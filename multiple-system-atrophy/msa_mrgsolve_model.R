## =====================================================================
##  Multiple System Atrophy (MSA) — QSP model for mrgsolve
##  ---------------------------------------------------------------------
##  ORGANISING IDEA
##  ---------------
##  Upright blood pressure is written as a PRODUCT of independently
##  breakable gains, not as a rate:
##
##      dMAP(upright) = [baroreceptor afferent]           (intact in MSA)
##                    x [central integration NTS/RVLM]    (G_CENT: dying)
##                    x [preganglionic IML output]        (G_CENT: dying)
##                    x [postganglionic releasable NE]    (POSTG: SPARED)
##                    x [vascular alpha1 density]         (A1R: UP-regulated)
##                    x [effective circulating volume]    (ECF: eroded at night)
##                    x [venous capacitance]              (CAPREL)
##
##  MSA breaks the CENTRAL terms and spares the POSTGANGLIONIC one.
##  Pure autonomic failure (PAF) and Parkinson disease (PD) do the reverse.
##  Every phenotype in this file is produced by changing only the VULN_*
##  vulnerability parameters — never by changing an equation.  Four clinical
##  facts then EMERGE rather than being coded:
##
##  (1) DRUG SELECTIVITY IS A MAP OF *WHERE* THE PRODUCT IS BROKEN.
##      Atomoxetine enters only as a multiplier on NET clearance of NEs, so
##      its pressor effect is proportional to NEREL = SNA x NEVES x POSTG.
##      In MSA (POSTG high, SNA low-but-present) it works; in the PAF
##      phenotype (POSTG -> 0) the same dose does essentially nothing.
##      Midodrine's active metabolite is added to the agonist term DIRECTLY,
##      downstream of every lesion, so it works in both.  Nothing in the
##      code says "atomoxetine is MSA-selective".
##
##  (2) CARBIDOPA ANTAGONISES DROXIDOPA.  Droxidopa becomes NE only through
##      DROXNE = KAADC*CDRX*POSTG*(1-CARBI).  Carbidopa, given for the
##      parkinsonism, occupies the SAME peripheral AADC term.  One shared
##      enzyme, written once, produces the interaction.
##
##  (3) THE PRESSOR CEILING IS SET AT NIGHT, NOT BY DAYTIME POTENCY.
##      Supine hypertension -> pressure natriuresis (UNAV, UV both rise
##      exponentially in MAP) -> overnight loss of NAB/ECF -> smaller Vc on
##      standing next morning -> WORSE orthostatic hypotension.  Because the
##      loop is closed, bedtime desmopressin and head-up-tilt sleeping raise
##      DAYTIME standing BP with no daytime pressor at all, and a long-acting
##      pressor at bedtime makes the next morning worse.
##
##  (4) LEVODOPA FAILURE IN MSA-P IS POSTSYNAPTIC.  Striatal output is
##      STRIAT = G_POST * DA/(EC50DA+DA) with G_POST = NMSN^GEXP_MSN.
##      Levodopa can only raise DA.  In the PD phenotype NMSN stays at 1 and
##      levodopa nearly normalises STRIAT; in MSA-P NMSN falls and the SAME
##      brain exposure buys a shrinking benefit that wanes over 1-2 years.
##
##  TWO TIMESCALES, ONE MODEL
##  -------------------------
##  Time unit is HOURS.  Neurodegeneration has time constants of years;
##  the baroreflex has time constants of seconds (TAU_MAP = 0.002 h ~ 7 s).
##  Rather than fight the stiffness in one long run, the driver below uses a
##  two-stage architecture:
##      msa_history()  — POSMODE=2 (circadian), 12 years, 12-h output
##      msa_state_at() — freeze the state vector at a chosen disease year
##      msa_tilt()     — restart from that state, POSMODE=3, 40-min tilt at
##                       30-s resolution  -> DeltaSBP(3 min), DeltaHR/DeltaSBP
##      msa_day()      — restart from that state, 48-h circadian run with
##                       real dosing  -> supine HTN, nocturnal natriuresis
##  The tilt test therefore becomes a READOUT OF DISEASE STAGE computed from
##  the same 61 equations that generate the 12-year trajectory.
##
##  UNITS
##  -----
##  time h · volumes L · MAP/SBP mmHg · HR bpm · sodium mmol · plasma drug
##  concentrations ng/mL (levodopa/brain levodopa ug/mL) · antibody mg/L ·
##  NES/NEVES/SNA/TPR/A1R dimensionless (1 = healthy supine reference) ·
##  neuronal states are SURVIVING FRACTIONS (1 = intact).
##
##  CALIBRATION TARGETS (see msa_references.md for the source of each)
##  ------------------------------------------------------------------
##   * nOH definition: standing SBP fall >= 30 mmHg (or DBP >= 15) at 3 min
##   * neurogenic signature: DeltaHR/DeltaSBP < 0.5 bpm/mmHg
##   * supine plasma NE normal (~250 pg/mL) with a blunted upright increment
##     (< ~60%, vs > 100% in health) — the MSA "central lesion" fingerprint
##   * UMSARS-II progresses ~5-8 points/year (natural-history cohorts)
##   * median survival 6-10 years from motor onset
##   * midodrine 10 mg: +15 to +25 mmHg standing SBP, peak ~1 h, ~4 h duration
##   * droxidopa 300 mg t.i.d.: ~ +8 to +12 mmHg standing SBP, OHSA item-1
##     improvement of ~0.9-1.0 units
##   * atomoxetine 18 mg: +25 to +50 mmHg seated SBP in MSA, ~nil in PAF
##   * fludrocortisone 0.1 mg: volume-mediated, days to plateau
##   * supine hypertension in roughly half of treated patients
##   * ~30% of MSA-P show a partial, waning levodopa response; no dyskinesia
##   * verdiperstat phase 3 (M-STAR) NEGATIVE — the model must reproduce a
##     null, not a benefit, when the intervention starts at diagnosis
##
##  RUN
##  ---
##      Rscript msa_mrgsolve_model.R          # builds, self-tests, 19 scenarios
##
##  Educational / research model.  NOT for clinical use.
## =====================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
})

msa_code <- '
$PROB
# Multiple System Atrophy — 61-ODE QSP model (time unit: hours)

$GLOBAL
#define SIGM(x, k) (1.0/(1.0 + exp(-(k)*(x))))
// smooth Gaussian meal bumps at 08:00, 13:00, 19:00
#define MEALB(t) (exp(-pow(((t)-8.0)/1.1,2)) + exp(-pow(((t)-13.0)/1.1,2)) + exp(-pow(((t)-19.0)/1.1,2)))

$PARAM @annotated
// ---------- protocol / posture ----------
POSMODE : 2    : 0 supine, 1 upright, 2 circadian, 3 tilt protocol
UPFRAC  : 0.55 : fraction of waking hours spent upright (circadian mode)
HUTF    : 0.0  : night-time orthostatic load from head-up tilt sleeping (0-0.3)
WAKE_H  : 7    : wake time (h of day)
SLEEP_H : 22   : bedtime (h of day)
TILT_T  : 0.25 : time of head-up tilt in tilt mode (h)
MEALAMP : 1.0  : post-prandial splanchnic shunt amplitude (0 = fasting)
BINDER  : 0    : abdominal binder / compression garment (0-1)
COUNTERM: 0    : physical counter-manoeuvres + slow rising (0-1)

// ---------- phenotype: regional vulnerability (this IS the phenotype) ----
KND     : 0.13 : neurodegeneration drive scale (per year at unit pathology)
VULN_IML : 1.00 : IML preganglionic sympathetic vulnerability
VULN_MED : 0.85 : NTS/CVLM/RVLM medullary cardiovascular vulnerability
VULN_VAG : 0.45 : cardiovagal (nucleus ambiguus) vulnerability
VULN_SN  : 0.70 : substantia nigra (presynaptic) vulnerability
VULN_MSN : 0.75 : striatal medium spiny neuron (POSTsynaptic) vulnerability
VULN_OPC : 0.35 : olivopontocerebellar vulnerability
VULN_ONUF: 1.05 : Onuf nucleus vulnerability
VULN_RESP: 0.40 : arcuate/raphe/pre-Botzinger respiratory vulnerability
VULN_PG  : 0.08 : POSTGANGLIONIC vulnerability — near-zero in MSA, high in PAF
COQ2F    : 1.0  : CoQ10 biosynthetic capacity (COQ2 variant carriers < 1)

// ---------- glial alpha-synuclein pathology (per-year rates) -------------
KSYN_A  : 1.00 : oligodendroglial alpha-syn monomer supply (per yr)
KUPT    : 2.20 : uptake of extracellular seed into oligodendrocyte (per yr)
KNUC    : 0.85 : p25alpha-templated nucleation (per yr)
KOXA    : 0.45 : oxidative promotion of oligomerisation (per yr)
KFIB    : 1.10 : oligomer -> fibrillar GCI (per yr)
KCLR    : 0.05 : GCI clearance (per yr)
KDEGM   : 1.60 : monomer degradation via autophagy/UPS (per yr)
KDEGO   : 0.55 : oligomer degradation via autophagy (per yr)
KAUT    : 1.30 : oligomer clogging of autophagic flux
KREL    : 0.30 : seed release from GCI (per yr)
KRELO   : 0.55 : seed release of oligomers (per yr)
KCLRS   : 3.20 : glymphatic/microglial seed clearance (per yr)
KP25    : 1.40 : p25alpha relocalisation rate (per yr)
KP25OFF : 0.35 : p25alpha return to myelin (per yr)
HP25    : 0.45 : baseline p25alpha availability floor for nucleation

// ---------- myelin, inflammation, mitochondria ---------------------------
KMYEREP : 0.30 : myelin repair (per yr)
KMYELOSS: 0.30 : myelin loss per unit GCI + oxidant (per yr)
KHOCL   : 0.55 : myeloperoxidase/HOCl contribution to myelin loss
KMGL    : 2.50 : microglial activation gain (per yr)
KDEB    : 0.60 : myelin-debris contribution to microglial activation
KMGLOFF : 2.20 : microglial deactivation (per yr)
KMPO    : 1.60 : myeloperoxidase induction by activated microglia (per yr)
KMPOOFF : 1.50 : myeloperoxidase turnover (per yr)
KOX     : 1.30 : oxidant generation from MPO (per yr)
KOXCQ   : 0.90 : oxidant generation from CoQ10 deficit
KOXOFF  : 2.00 : antioxidant handling (per yr)
TAU_CQ  : 0.50 : CoQ10 pool time constant (yr)
KCQOX   : 0.35 : oxidative consumption of CoQ10
WOXS    : 0.55 : weight of oxidative stress in neurodegeneration drive
WINF    : 0.40 : weight of neuroinflammation in the drive
WTROPH  : 0.85 : weight of LOST OLIGODENDROGLIAL TROPHIC SUPPORT in the drive
KNFL    : 55.0 : neurofilament-light release per unit neuronal loss rate
KNFLEL  : 0.35 : NfL elimination (per h)
NFL0    : 12.0 : baseline plasma NfL (pg/mL)

// ---------- central autonomic gain --------------------------------------
// The reflex arc is a SERIAL CHAIN (NTS -> CVLM -> RVLM -> IML), so a partial
// loss at each relay multiplies.  Exponents > 1 encode that compounding; with
// exponents of 1 a patient who has lost 40% of each relay still retains ~35%
// of reflex gain, which is far too much function to match the clinic.
GEXP_MED: 2.20 : exponent on medullary survival in central gain
GEXP_IML: 1.60 : exponent on IML survival in central gain
FTON    : 0.22 : lesion-independent (spinal/local) fraction of tonic drive
SNA0    : 1.00 : reference tonic sympathetic activity
SNAMAX  : 2.50 : ceiling on preganglionic sympathetic drive
GBARO   : 0.160: baroreflex gain (per mmHg error)
MAP_SET : 92   : baroreflex operating point (mmHg)
TAU_SNA : 0.003: sympathetic activity time constant (h ~ 11 s)

// ---------- postganglionic NE handling ----------------------------------
KREL_NE : 43.0 : NE release rate constant
KA2     : 6.00 : alpha2 autoreceptor affinity (NES units)
KUP     : 30.0 : NET reuptake clearance (per h)
KNEDEG  : 6.00 : MAO/COMT clearance of synaptic NE (per h)
KVES    : 2.40 : vesicular NE repletion (per h)
KUSE    : 0.022: vesicular depletion per unit release
KRECYC  : 0.035: fraction of reuptaken NE recycled into vesicles
NEPL0   : 250  : supine plasma NE reference (pg/mL) at reference release
KNEREF  : 38.5 : healthy reference NE release rate (spillover denominator)
KSPILL  : 1.00 : extra plasma spillover per unit NET blockade
KAADC   : 0.007: peripheral AADC conversion of droxidopa (per ng/mL)
IC50CBD : 400  : carbidopa concentration for 50% AADC block (ng/mL)
ENET    : 0.92 : maximal fractional NET blockade by atomoxetine
EC50ATX : 20   : atomoxetine EC50 for NET blockade (ng/mL)
EYOH    : 0.70 : maximal alpha2 blockade by yohimbine
EC50YOH : 55   : yohimbine EC50 (ng/mL)

// ---------- vascular effector ------------------------------------------
KSUPER  : 0.60 : denervation supersensitivity gain on alpha1 density
TAU_A1  : 336  : alpha1 receptor up-regulation time constant (h, ~2 wk)
EC50A1  : 1.00 : agonist occupancy EC50 (NES units)
KDMM_EQ : 45   : ng/mL desglymidodrine equivalent to 1 NES unit
ETPR    : 1.10 : maximal fractional rise in TPR from alpha1 occupancy
TPR0    : 1.00 : reference peripheral resistance
TAU_TPR : 0.05 : vascular tone time constant (h ~ 3 min)
EWB     : 0.28 : water-bolus osmopressor maximal TPR effect
EC50WB  : 210  : water-bolus effect-site EC50 (mL equivalent)
KVEN    : 2.60 : venous capacitance reduction per unit alpha1 occupancy
KBIND   : 0.55 : venous capacitance reduction from abdominal binder
KOCTV   : 0.60 : splanchnic constriction from octreotide (max fraction)
EC50OCT : 0.35 : octreotide effect-site EC50 (ug equivalent)
VPOOL   : 1.05 : maximal orthostatic splanchnic/dependent pooling (L)
VMEAL   : 0.30 : post-prandial additional splanchnic pooling (L)
TAU_POOL: 0.008: pooling time constant (h ~ 30 s)

// ---------- cardiac and haemodynamics ----------------------------------
BV0     : 5.00 : reference blood volume (L)
KBV     : 0.55 : blood volume sensitivity to ECF
ECF0    : 14.0 : reference extracellular fluid volume (L)
FCENT   : 0.28 : fraction of blood volume that is cardiopulmonary/central
VC0     : 1.40 : reference central blood volume (L) = FCENT*BV0
TAU_VC  : 0.006: central volume equilibration (h)
SV0     : 72   : reference stroke volume (mL)
HSV     : 0.85 : Frank-Starling exponent on central volume
KINO    : 0.12 : sympathetic inotropy gain
HR0     : 66   : intrinsic heart rate reference (bpm)
KCHR    : 0.30 : sympathetic chronotropy gain
KVAG    : 0.28 : vagal bradycardic gain
VAG0    : 1.00 : reference vagal tone
KVAGP   : 0.85 : vagal withdrawal on standing
TAU_VAG : 0.004: vagal time constant (h)
TAU_HR  : 0.005: heart-rate time constant (h)
KMAP    : 11.1 : MAP per unit cardiac output x resistance (mmHg / (L/min))
TAU_MAP : 0.002: arterial pressure time constant (h ~ 7 s)
KPP     : 0.62 : pulse pressure per mL stroke volume (mmHg/mL)
CBFFLOOR: 58   : MAP below which cerebral autoregulation fails (mmHg)

// ---------- renal / volume / endocrine ---------------------------------
NAIN    : 6.25 : dietary sodium intake (mmol/h; 150 mmol/day)
WIN     : 0.094: water intake (L/h; ~2.25 L/day)
UNA0    : 6.25 : reference sodium excretion (mmol/h)
MAPNAT  : 92   : MAP at which natriuresis equals intake (mmHg)
KPN     : 0.30 : PRESSURE-natriuresis steepness (per 10 mmHg) — deliberately
KVOLN   : 6.0  : VOLUME-error natriuresis steepness (tight homeostasis)
ECFTGT  : 14.0 : renal extracellular-volume target (L)
KMR     : 0.80 : exponent of mineralocorticoid sodium retention
UV0     : 0.130: reference urine output scale (L/h)
KPNW    : 0.22 : pressure-diuresis steepness (per 10 mmHg)
KVOLW   : 4.0  : volume-error diuresis steepness
KTHIRST : 1.50 : osmotic thirst gain (fractional water intake per 10% Na rise)
FAQ     : 0.62 : maximal fractional water reabsorption from AQP2
KAQ     : 0.55 : AQP2 EC50 for water reabsorption
EAQ     : 1.00 : maximal AQP2 trafficking
EC50AVP : 2.60 : AVP EC50 for V2R/AQP2 (pg/mL)
KDDA    : 3.10 : desmopressin AVP-equivalents per ng/mL
TAU_AQP : 1.20 : AQP2 trafficking time constant (h)
AVP0    : 2.10 : basal plasma AVP (pg/mL)
KOSM    : 60.0 : osmotic AVP gain (per fractional plasma-sodium change)
KBAVP   : 0.16 : BARORECEPTOR AVP gain (per mmHg) — gated by G_CENT
TAU_AVP : 0.35 : AVP time constant (h)
PRA0    : 1.00 : reference plasma renin activity
FRTON   : 0.30 : sympathetic-independent renin fraction
KBRENIN : 0.70 : beta1-sympathetic renin gain
KPRAP   : 1.30 : renal perfusion-pressure renin gain
TAU_PRA : 0.60 : renin time constant (h)
ALD0    : 1.00 : reference aldosterone
TAU_ALD : 1.50 : aldosterone time constant (h)
EFLU    : 2.30 : maximal fludrocortisone mineralocorticoid effect
EC50FLU : 1.40 : fludrocortisone effect-site EC50 (ng/mL)
EPYR    : 0.85 : maximal pyridostigmine amplification of BAROREFLEX traffic
EC50PYR : 26   : pyridostigmine EC50 (ng/mL)

// ---------- dopaminergic ------------------------------------------------
KDAEND  : 1.00 : endogenous striatal dopamine per unit nigral survival
KDALD   : 0.85 : dopamine produced per ug/mL brain levodopa
FGLIA   : 0.30 : fraction of levodopa decarboxylation independent of terminals
EC50DA  : 0.55 : striatal dopamine EC50 for postsynaptic signalling
GEXP_MSN: 1.15 : exponent on MSN survival in postsynaptic gain
TAU_DA  : 0.30 : striatal dopamine time constant (h)
PMAX    : 26   : maximal parkinsonism contribution to UMSARS-II
AMAX    : 22   : maximal ataxia contribution to UMSARS-II
PVRMAX  : 320  : maximal post-void residual (mL)
SMAX    : 1.00 : maximal stridor severity index

// ---------- clinical endpoints -----------------------------------------
SBPTH   : 96   : standing SBP at which OH symptoms are half-maximal (mmHg)
OHSLOPE : 9.0  : steepness of the symptom-vs-standing-SBP curve (mmHg)
OHMAX   : 10.0 : maximal OHSA item-1 style score
TAU_OHSA: 24   : symptom-score time constant (h)
TAU_U   : 720  : UMSARS time constant (h, ~1 month)
WPARK   : 1.00 : weight of parkinsonism in UMSARS-II
WATAX   : 1.00 : weight of ataxia in UMSARS-II
WBULB   : 6.0  : weight of bulbar/respiratory involvement in UMSARS-II
WOH     : 0.55 : weight of OH symptoms in UMSARS-I
WURO    : 3.5  : weight of urogenital failure in UMSARS-I
WS_IML  : 12   : SCOPA-AUT weight, sympathetic failure
WS_URO  : 10   : SCOPA-AUT weight, urogenital
WS_GI   : 8    : SCOPA-AUT weight, gastrointestinal (via DMV proxy NRESP)
WS_OH   : 1.4  : SCOPA-AUT weight, orthostatic symptoms
H0      : 0.0000017 : baseline hazard (per h)
B_U2    : 0.42 : log-hazard per 10 UMSARS-II points
B_STR   : 0.85 : log-hazard for stridor
B_OH    : 0.22 : log-hazard per 5 OHSA points
B_RESP  : 1.15 : log-hazard for respiratory nucleus loss

// ---------- drug pharmacokinetics --------------------------------------
KA_MID  : 3.00 : midodrine absorption (per h)
KM_MID  : 1.85 : midodrine -> desglymidodrine conversion (per h)
KE_MID  : 0.45 : midodrine non-converting elimination (per h)
VMID    : 60   : midodrine central volume (L)
KE_DMM  : 0.231: desglymidodrine elimination (per h; t1/2 3 h)
VDMM    : 95   : desglymidodrine volume (L)
KA_DRX  : 0.90 : droxidopa absorption (per h)
KE_DRX  : 0.277: droxidopa elimination (per h; t1/2 2.5 h)
VDRX    : 35   : droxidopa volume (L)
KE_FLU  : 0.029: fludrocortisone effective elimination (per h; t1/2 24 h)
VFLU    : 45   : fludrocortisone volume (L)
KA_FLU  : 1.60 : fludrocortisone absorption (per h)
KE_ATX  : 0.133: atomoxetine elimination (per h; t1/2 5.2 h extensive metaboliser)
CYP2D6  : 1.0  : CYP2D6 activity multiplier (0.24 = poor metaboliser)
VATX    : 62   : atomoxetine volume (L)
KA_ATX  : 2.00 : atomoxetine absorption (per h)
KE_PYR  : 0.198: pyridostigmine elimination (per h; t1/2 3.5 h)
VPYR    : 70   : pyridostigmine apparent volume (L)
KA_PYR  : 1.00 : pyridostigmine absorption (per h)
KE_DDA  : 0.231: desmopressin elimination (per h; t1/2 3 h)
VDDA    : 25   : desmopressin volume (L)
KE_LD   : 0.462: levodopa elimination (per h; t1/2 1.5 h with carbidopa)
VLD     : 40   : levodopa central volume (L)
KA_LD   : 2.00 : levodopa absorption (per h, gastroparesis-gated)
KGP     : 0.55 : fractional reduction of KA_LD by gastroparesis
KIN_LDB : 0.35 : levodopa brain influx (LAT1) (per h)
KOUT_LDB: 0.95 : brain levodopa efflux + decarboxylation (per h)
VLDB    : 1.40 : brain levodopa distribution volume (L)
KE_CBD  : 0.347: carbidopa elimination (per h; t1/2 2 h)
VCBD    : 50   : carbidopa volume (L)
KA_CBD  : 1.50 : carbidopa absorption (per h)
KE_UBQ  : 0.020: ubiquinol elimination (per h; t1/2 ~34 h)
VUBQ    : 120  : ubiquinol volume (L)
KA_UBQ  : 0.25 : ubiquinol absorption (per h, low bioavailability)
EUBQ    : 0.55 : maximal CoQ10 repletion effect
EC50UBQ : 2200 : ubiquinol EC50 (ng/mL)
KE_VRD  : 0.058: verdiperstat elimination (per h; t1/2 12 h)
VVRD    : 55   : verdiperstat volume (L)
KA_VRD  : 1.20 : verdiperstat absorption (per h)
EVRD    : 0.85 : maximal myeloperoxidase inhibition
EC50VRD : 900  : verdiperstat EC50 (ng/mL)
KE_MAB  : 0.00138 : anti-alpha-syn mAb elimination (per h; t1/2 ~21 d)
VMAB    : 3.60 : antibody central volume (L)
KMABSEED: 0.42 : seed neutralisation per mg/L antibody (per yr)
KE_YOH  : 0.462: yohimbine elimination (per h)
VYOH    : 80   : yohimbine volume (L)
KA_YOH  : 2.20 : yohimbine absorption (per h)
KE_OCT  : 0.385: octreotide elimination (per h; t1/2 1.8 h)
VOCT    : 20   : octreotide volume (L)
KE_WB   : 0.92 : water-bolus effect-site decay (per h)
IASO    : 0.0  : SNCA-lowering ASO fractional knockdown (0-0.8)
GLYMPH  : 0.0  : additional glymphatic seed clearance (0-1)

// ---------- multiscale housekeeping ------------------------------------
// TSCALE stretches ONLY the fast cardiovascular time constants (seconds to
// minutes).  The disease layer, the renal/endocrine layer and every
// quasi-steady VALUE are untouched, so a multi-year run with TSCALE = 40
// (TAU_MAP 7 s -> 5 min) reproduces the same trajectory as TSCALE = 1 at a
// small fraction of the cost: a 61-state stiff system with a numerical
// Jacobian otherwise spends essentially all of its time resolving a
// 7-second transient that no annual endpoint depends on.  Tilt runs and
// 48-h runs use TSCALE = 1, where those seconds ARE the measurement.
// msa_selftest() asserts the two agree to better than 1%.
TSCALE  : 1    : stretch factor on the fast haemodynamic time constants
// INITMAIN gates the explicit compartment initialisation in $MAIN.  It must
// be 1 for a fresh healthy start and 0 whenever a run is RESTARTED from a
// saved state via init(): mrgsolve applies $MAIN after init(), so leaving
// the X_0 assignments unguarded silently resets a diseased patient back to
// healthy physiology and every tilt/48-h experiment then reports the SAME
// numbers regardless of disease stage or phenotype.
INITMAIN: 1    : 1 = initialise from parameters, 0 = keep init() state
HPY     : 8766 : hours per year (converts per-year rates to per-hour)

$CMT @annotated
// --- glial proteinopathy (10) ---
P25A  : p25alpha/TPPP mislocalised to soma (fraction)
ASYNM : oligodendroglial alpha-syn monomer (a.u.)
ASYNO : alpha-syn oligomer (a.u.)
GCI   : glial cytoplasmic inclusion burden (a.u.)
SEED  : extracellular seeding-competent alpha-syn (a.u.)
MGL   : activated microglia (1 = resting baseline)
MPO   : myeloperoxidase activity (a.u.)
OXS   : oxidative stress (a.u.)
CQ    : coenzyme Q10 pool (1 = normal)
MYE   : myelin integrity (1 = normal)
// --- neuronal populations, surviving fraction (9) ---
NIML  : spinal IML preganglionic sympathetic neurons
NMED  : NTS/CVLM/RVLM medullary cardiovascular neurons
NVAG  : cardiovagal preganglionic neurons
NSN   : substantia nigra dopaminergic neurons (PREsynaptic)
NMSN  : striatal medium spiny neurons (POSTsynaptic)
NOPC  : olivopontocerebellar neurons
NONUF : Onuf nucleus motor neurons
NRESP : arcuate/raphe/pre-Botzinger respiratory neurons
POSTG : POSTGANGLIONIC sympathetic terminal integrity
// --- biomarker (1) ---
NFL   : plasma neurofilament light (pg/mL)
// --- cardiovascular fast loop (10) ---
SNA   : sympathetic nerve activity (1 = healthy supine)
NES   : synaptic norepinephrine at vascular alpha1 (a.u.)
NEVES : releasable vesicular NE pool (1 = full)
A1R   : vascular alpha1-adrenoceptor density (1 = normal)
TPR   : total peripheral resistance (1 = reference)
VSP   : splanchnic/dependent pooled volume (L)
VC    : central blood volume (L)
HR    : heart rate (bpm)
VAG   : cardiovagal tone (1 = reference)
MAP   : mean arterial pressure (mmHg)
// --- renal / volume / endocrine (6) ---
NAB   : total body exchangeable sodium (mmol)
ECF   : extracellular fluid volume (L)
AQP2  : collecting-duct AQP2 apical trafficking (1 = maximal)
AVP   : plasma arginine vasopressin (pg/mL)
PRA   : plasma renin activity (1 = reference)
ALD   : aldosterone (1 = reference)
// --- dopaminergic (1) ---
DA    : striatal synaptic dopamine (a.u.)
// --- drug PK (19) ---
MIDG  : midodrine gut (mg)
MIDC  : midodrine plasma (mg)
DMM   : desglymidodrine plasma (mg)
DRXG  : droxidopa gut (mg)
DRXC  : droxidopa plasma (mg)
FLUG  : fludrocortisone gut (mg)
FLU   : fludrocortisone plasma (mg)
ATXG  : atomoxetine gut (mg)
ATX   : atomoxetine plasma (mg)
PYRG  : pyridostigmine gut (mg)
PYR   : pyridostigmine plasma (mg)
DDA   : desmopressin plasma (ug)
LDG   : levodopa gut (mg)
LDC   : levodopa plasma (mg)
LDB   : brain levodopa (mg)
CBDG  : carbidopa gut (mg)
CBD   : carbidopa plasma (mg)
UBQG  : ubiquinol gut (mg)
UBQ   : ubiquinol plasma (mg)
// --- drug PK, continued (5) ---
VRDG  : verdiperstat gut (mg)
VRD   : verdiperstat plasma (mg)
MAB   : anti-alpha-syn antibody plasma (mg)
YOHG  : yohimbine gut (mg)
YOH   : yohimbine plasma (mg)
OCT   : octreotide plasma (ug)
WB    : water-bolus effect site (mL equivalent)
// --- clinical endpoints (5) ---
OHSA  : orthostatic-hypotension symptom score (0-10)
U1    : UMSARS Part I (activities of daily living)
U2    : UMSARS Part II (motor examination)
SCOPA : SCOPA-AUT autonomic symptom score
CUMH  : cumulative mortality hazard

$MAIN
// ---- initial conditions: healthy physiology, zero pathology ----
// Guarded by INITMAIN so that init() can restart from a frozen disease state.
if(INITMAIN > 0.5) {
P25A_0  = HP25;
ASYNM_0 = KSYN_A/KDEGM;
ASYNO_0 = 0.0;
GCI_0   = 0.0;
SEED_0  = 0.0;
MGL_0   = 1.0;
MPO_0   = 1.0;
OXS_0   = 1.0;
CQ_0    = COQ2F;
MYE_0   = 1.0;

NIML_0 = 1.0;  NMED_0  = 1.0;  NVAG_0 = 1.0;  NSN_0   = 1.0;
NMSN_0 = 1.0;  NOPC_0  = 1.0;  NONUF_0= 1.0;  NRESP_0 = 1.0;
POSTG_0= 1.0;
NFL_0  = NFL0;

SNA_0   = SNA0;
NES_0   = 1.0;
NEVES_0 = 1.0;
A1R_0   = 1.0;
TPR_0   = TPR0;
VSP_0   = 0.0;
VC_0    = FCENT*BV0;
HR_0    = HR0;
VAG_0   = VAG0;
MAP_0   = MAP_SET;

ECF_0   = ECF0;
NAB_0   = 140.0*ECF0;
AQP2_0  = 0.45;
AVP_0   = AVP0;
PRA_0   = PRA0;
ALD_0   = ALD0;

DA_0    = KDAEND;

OHSA_0  = 0.0;
U1_0    = 0.0;
U2_0    = 0.0;
SCOPA_0 = 0.0;
CUMH_0  = 0.0;
}

F_LDG = 1.0;
ALAG_LDG = 0.0;

$ODE
// =====================================================================
// 0. posture / meal drivers (smooth, so the solver never sees a step)
// =====================================================================
double tod  = fmod(SOLVERTIME, 24.0);
double wake = SIGM(tod - WAKE_H, 8.0) * SIGM(SLEEP_H - tod, 8.0);
double POS;
if(POSMODE < 0.5)      POS = 0.0;
else if(POSMODE < 1.5) POS = 1.0;
else if(POSMODE < 2.5) POS = UPFRAC*wake + HUTF*(1.0 - wake);
else                   POS = SIGM(SOLVERTIME - TILT_T, 120.0);
POS = POS*(1.0 - 0.25*COUNTERM);            // counter-manoeuvres offload the pool
double meal = MEALAMP*wake*MEALB(tod);

// =====================================================================
// 1. glial alpha-synuclein pathology  (rates given per year -> per hour)
// =====================================================================
double AUTOF = 1.0/(1.0 + KAUT*ASYNO);               // autophagic flux, clogged by oligomers
double p25   = HP25 + (1.0 - HP25)*P25A;             // effective nucleation template
double CMAB  = MAB/VMAB;                             // mg/L

dxdt_P25A  = (KP25*(1.0 - MYE)*(1.0 - P25A) - KP25OFF*P25A)/HPY;
dxdt_ASYNM = (KSYN_A*(1.0 - IASO) + KUPT*SEED
              - KNUC*ASYNM*p25 - KOXA*OXS*ASYNM*0.5
              - KDEGM*ASYNM*AUTOF)/HPY;
dxdt_ASYNO = (KNUC*ASYNM*p25 + KOXA*OXS*ASYNM*0.5
              - KFIB*ASYNO - KDEGO*ASYNO*AUTOF)/HPY;
dxdt_GCI   = (KFIB*ASYNO - KCLR*GCI)/HPY;
dxdt_SEED  = (KREL*GCI + KRELO*ASYNO
              - KCLRS*(1.0 + GLYMPH)*SEED
              - KMABSEED*CMAB*SEED)/HPY;

// =====================================================================
// 2. myelin, neuroinflammation, mitochondria
// =====================================================================
double MGLX = (MGL > 1.0) ? (MGL - 1.0) : 0.0;
double CUBQ = 1000.0*UBQ/VUBQ;                       // ng/mL
double CVRD = 1000.0*VRD/VVRD;                       // ng/mL
double IVRD = EVRD*CVRD/(EC50VRD + CVRD);
double CQss = COQ2F*(1.0 + EUBQ*CUBQ/(EC50UBQ + CUBQ))/(1.0 + KCQOX*(OXS - 1.0 > 0 ? OXS - 1.0 : 0.0));

dxdt_MGL = (KMGL*(SEED + KDEB*(1.0 - MYE)) - KMGLOFF*MGLX)/HPY;
dxdt_MPO = (KMPO*MGLX*(1.0 - IVRD) - KMPOOFF*(MPO - 1.0))/HPY;
dxdt_OXS = (KOX*(MPO - 1.0) + KOXCQ*(1.0 - CQ) - KOXOFF*(OXS - 1.0))/HPY;
dxdt_CQ  = (CQss - CQ)/(TAU_CQ*HPY);
dxdt_MYE = (KMYEREP*(1.0 - MYE) - KMYELOSS*(GCI + KHOCL*(MPO - 1.0))*MYE)/HPY;

// =====================================================================
// 3. neurodegeneration — ONE drive, many vulnerabilities
//    The trophic term is what makes MSA a GLIAL disease: neurons die
//    because their oligodendrocytes stopped supporting them.
// =====================================================================
double NDRIVE = KND*(GCI + WOXS*(OXS - 1.0) + WINF*MGLX + WTROPH*(1.0 - MYE))/HPY;
if(NDRIVE < 0.0) NDRIVE = 0.0;

dxdt_NIML  = -NDRIVE*VULN_IML *NIML;
dxdt_NMED  = -NDRIVE*VULN_MED *NMED;
dxdt_NVAG  = -NDRIVE*VULN_VAG *NVAG;
dxdt_NSN   = -NDRIVE*VULN_SN  *NSN;
dxdt_NMSN  = -NDRIVE*VULN_MSN *NMSN;
dxdt_NOPC  = -NDRIVE*VULN_OPC *NOPC;
dxdt_NONUF = -NDRIVE*VULN_ONUF*NONUF;
dxdt_NRESP = -NDRIVE*VULN_RESP*NRESP;
dxdt_POSTG = -NDRIVE*VULN_PG  *POSTG;

double LOSSRATE = NDRIVE*(VULN_IML*NIML + VULN_MED*NMED + VULN_SN*NSN
                        + VULN_MSN*NMSN + VULN_OPC*NOPC + VULN_RESP*NRESP)/6.0;
dxdt_NFL = KNFL*LOSSRATE*HPY - KNFLEL*(NFL - NFL0);

// =====================================================================
// 4. THE PRODUCT — central gain vs postganglionic integrity
// =====================================================================
double G_CENT = pow(NMED > 1e-6 ? NMED : 1e-6, GEXP_MED)
              * pow(NIML > 1e-6 ? NIML : 1e-6, GEXP_IML);

double CPYR  = 1000.0*PYR/VPYR;
double PYRAMP= 1.0 + EPYR*CPYR/(EC50PYR + CPYR);     // amplifies REFLEX traffic only
double berr  = MAP_SET - MAP;

double SNAss = SNA0*(FTON + (1.0 - FTON)*G_CENT) + GBARO*G_CENT*PYRAMP*berr;
if(SNAss < 0.02) SNAss = 0.02;
if(SNAss > SNAMAX) SNAss = SNAMAX;   // preganglionic drive cannot rise without limit
dxdt_SNA = (SNAss - SNA)/(TAU_SNA*TSCALE);

// =====================================================================
// 5. norepinephrine handling at the postganglionic terminal
// =====================================================================
double CYOH  = 1000.0*YOH/VYOH;
double IYOH  = EYOH*CYOH/(EC50YOH + CYOH);
double A2INH = 1.0/(1.0 + (NES/KA2)*(1.0 - IYOH));
double NEREL = KREL_NE*SNA*NEVES*POSTG*A2INH;

double CDRX  = 1000.0*DRXC/VDRX;
double CCBD  = 1000.0*CBD/VCBD;
double CARBI = CCBD/(IC50CBD + CCBD);                // <-- shared AADC term
double DROXNE= KAADC*CDRX*POSTG*(1.0 - CARBI);       // needs terminals AND no carbidopa

double CATX  = 1000.0*ATX/VATX;
double NETI  = ENET*CATX/(EC50ATX + CATX);
// NET is expressed ON the postganglionic terminal, so reuptake capacity is
// gated by POSTG.  This is what makes NET inhibition phenotype-selective:
// where the terminal is gone there is no transporter left to block and
// atomoxetine has nothing to do.  Without the gate a PAF patient gets almost
// the same atomoxetine response as an MSA patient, which is the opposite of
// what the crossover trials show.
double NEUP  = KUP*POSTG*(1.0 - NETI)*NES;

dxdt_NES   = NEREL + DROXNE - NEUP - KNEDEG*NES;
dxdt_NEVES = KVES*POSTG*(1.0 - NEVES) - KUSE*NEREL + KRECYC*NEUP;

// =====================================================================
// 6. vascular effector — denervation supersensitivity multiplies EVERY
//    alpha1 agonist, endogenous or pharmacological
// =====================================================================
double A1ss = 1.0 + KSUPER*(1.0 - G_CENT);
dxdt_A1R = (A1ss - A1R)/TAU_A1;

double CDMM = 1000.0*DMM/VDMM;
double AGON_N = NES;
double AGON = NES + CDMM/KDMM_EQ;
double OCC  = AGON/(EC50A1 + AGON);
double WBEF = EWB*WB/(EC50WB + WB);
double TPRss= TPR0*(1.0 + ETPR*A1R*OCC + WBEF);
dxdt_TPR = (TPRss - TPR)/(TAU_TPR*TSCALE);

double COCT   = OCT/VOCT;
double OCTEF  = KOCTV*COCT/(EC50OCT + COCT);
double CAPREL = 1.0/(1.0 + KVEN*A1R*OCC + KBIND*BINDER);
double POOLss = (POS*VPOOL + meal*VMEAL*(1.0 - OCTEF))*CAPREL*(1.0 - 0.55*OCTEF*0.0);
dxdt_VSP = (POOLss - VSP)/(TAU_POOL*TSCALE);

// =====================================================================
// 7. cardiac output and arterial pressure
// =====================================================================
double BV   = BV0*(1.0 + KBV*(ECF/ECF0 - 1.0));
double VCss = FCENT*BV - VSP;
if(VCss < 0.15) VCss = 0.15;
dxdt_VC = (VCss - VC)/(TAU_VC*TSCALE);

double SV_  = SV0*pow(VC/VC0 > 0.05 ? VC/VC0 : 0.05, HSV)*(1.0 + KINO*SNA*POSTG);
double GMED = pow(NMED > 1e-6 ? NMED : 1e-6, GEXP_MED);
double VAGss= VAG0*pow(NVAG > 1e-6 ? NVAG : 1e-6, 1.0)*(1.0 - KVAGP*POS*GMED);
dxdt_VAG = (VAGss - VAG)/(TAU_VAG*TSCALE);
double HRss = HR0*(1.0 + KCHR*SNA*POSTG)/(1.0 + KVAG*VAG);
dxdt_HR = (HRss - HR)/(TAU_HR*TSCALE);

double CO_  = SV_*HR/1000.0;                          // L/min
double MAPss= KMAP*CO_*TPR;
dxdt_MAP = (MAPss - MAP)/(TAU_MAP*TSCALE);

// =====================================================================
// 8. renal sodium / water — THE NOCTURNAL LOOP
//    Both natriuresis and diuresis rise exponentially with MAP, so supine
//    hypertension at night spends the very volume the patient needs at 07:00.
// =====================================================================
double CFLU = 1000.0*FLU/VFLU;
double MR   = ALD/ALD0 + EFLU*CFLU/(EC50FLU + CFLU);
if(MR < 0.05) MR = 0.05;
// Sodium excretion carries TWO terms with deliberately different gains:
// a STEEP volume-error term (long-run homeostasis, keeps total body sodium
// bounded) and a SHALLOW pressure term (~1.5-2x over a 25 mmHg supine rise).
// It is the shallow term that closes the nocturnal loop; making it steep
// would make the loop implausibly violent, making it zero would abolish the
// central clinical trade-off of the disease.
double VOLERR = ECF/ECFTGT;
if(VOLERR < 0.3) VOLERR = 0.3;
double UNAV = UNA0*pow(VOLERR, KVOLN)*exp(KPN*(MAP - MAPNAT)/10.0)/pow(MR, KMR);
dxdt_NAB = NAIN - UNAV;
if(NAB < 200.0 && dxdt_NAB < 0.0) dxdt_NAB = 0.0;

double CDDA = 1000.0*DDA/VDDA;                        // pg/mL-equivalent scale
double AQPss= EAQ*(AVP + KDDA*CDDA)/(EC50AVP + AVP + KDDA*CDDA);
dxdt_AQP2 = (AQPss - AQP2)/TAU_AQP;

double UVOL = UV0*pow(VOLERR, KVOLW)*exp(KPNW*(MAP - MAPNAT)/10.0)
              *(1.0 - FAQ*AQP2/(KAQ + AQP2));
if(UVOL < 0.005) UVOL = 0.005;
double NACONC = NAB/(ECF > 1.0 ? ECF : 1.0);
// osmotic thirst — the intact half of volume defence in MSA
double WINE = WIN*(1.0 + KTHIRST*(NACONC/140.0 - 1.0)*10.0);
if(WINE < 0.3*WIN) WINE = 0.3*WIN;
dxdt_ECF = WINE - UVOL;
double AVPss  = AVP0*(1.0 + KOSM*(NACONC/140.0 - 1.0))
              + KBAVP*G_CENT*(berr > 0.0 ? berr : 0.0);   // baroreflex arm is GATED
if(AVPss > 40.0) AVPss = 40.0;
if(AVPss < 0.1) AVPss = 0.1;
dxdt_AVP = (AVPss - AVP)/TAU_AVP;

double PRAss = PRA0*(FRTON + KBRENIN*SNA*POSTG)
               *pow(MAPNAT/(MAP > 20.0 ? MAP : 20.0), KPRAP);
dxdt_PRA = (PRAss - PRA)/TAU_PRA;
double ALDss = ALD0*pow(PRA/PRA0 > 0.05 ? PRA/PRA0 : 0.05, 0.8);
dxdt_ALD = (ALDss - ALD)/TAU_ALD;

// =====================================================================
// 9. dopaminergic: PREsynaptic supply vs POSTsynaptic target
// =====================================================================
double CLDB = LDB/VLDB;                               // mg/L = ug/mL
double DAss = KDAEND*NSN + KDALD*CLDB*(FGLIA + (1.0 - FGLIA)*NSN);
dxdt_DA = (DAss - DA)/TAU_DA;

// =====================================================================
// 10. drug pharmacokinetics
// =====================================================================
dxdt_MIDG = -KA_MID*MIDG;
dxdt_MIDC =  KA_MID*MIDG - (KM_MID + KE_MID)*MIDC;
dxdt_DMM  =  KM_MID*MIDC - KE_DMM*DMM;

dxdt_DRXG = -KA_DRX*DRXG;
dxdt_DRXC =  KA_DRX*DRXG - KE_DRX*DRXC;

dxdt_FLUG = -KA_FLU*FLUG;
dxdt_FLU  =  KA_FLU*FLUG - KE_FLU*FLU;

dxdt_ATXG = -KA_ATX*ATXG;
dxdt_ATX  =  KA_ATX*ATXG - (KE_ATX*CYP2D6)*ATX;

dxdt_PYRG = -KA_PYR*PYRG;
dxdt_PYR  =  KA_PYR*PYRG - KE_PYR*PYR;

dxdt_DDA  = -KE_DDA*DDA;

double KA_LDe = KA_LD*(1.0 - KGP*(1.0 - NRESP));      // gastroparesis gates absorption
dxdt_LDG = -KA_LDe*LDG;
dxdt_LDC =  KA_LDe*LDG - KE_LD*LDC - KIN_LDB*LDC;
dxdt_LDB =  KIN_LDB*LDC - KOUT_LDB*LDB;

dxdt_CBDG = -KA_CBD*CBDG;
dxdt_CBD  =  KA_CBD*CBDG - KE_CBD*CBD;

dxdt_UBQG = -KA_UBQ*UBQG;
dxdt_UBQ  =  KA_UBQ*UBQG - KE_UBQ*UBQ;

dxdt_VRDG = -KA_VRD*VRDG;
dxdt_VRD  =  KA_VRD*VRDG - KE_VRD*VRD;

dxdt_MAB = -KE_MAB*MAB;

dxdt_YOHG = -KA_YOH*YOHG;
dxdt_YOH  =  KA_YOH*YOHG - KE_YOH*YOH;

dxdt_OCT = -KE_OCT*OCT;
dxdt_WB  = -KE_WB*WB;

// =====================================================================
// 11. clinical endpoints
// =====================================================================
// Plasma norepinephrine is measured SPILLOVER — a function of RELEASE and of
// how much escapes reuptake — not of the synaptic concentration that drives
// the receptor.  Reading it off NES makes a phenotype with dead terminals but
// intact central drive look near-normal, inverting the single most useful
// bedside discriminator between MSA (normal supine NE, blunted increment) and
// PAF (low supine NE).
double NEPLA = NEPL0*(NEREL/KNEREF)*(1.0 + KSPILL*NETI);

double PPr  = KPP*SV_;
double SBPn = MAP + PPr/2.0;
double OHSt = OHMAX/(1.0 + exp((SBPn - SBPTH)/OHSLOPE));
dxdt_OHSA = (OHSt - OHSA)/TAU_OHSA;

double G_POST = pow(NMSN > 1e-6 ? NMSN : 1e-6, GEXP_MSN);
double STRIAT = G_POST*DA/(EC50DA + DA);
double STR0   = KDAEND/(EC50DA + KDAEND);
double PARK   = PMAX*(1.0 - STRIAT/STR0);
if(PARK < 0.0) PARK = 0.0;
double ATAX   = AMAX*(1.0 - NOPC);
double PVRES  = PVRMAX*(1.0 - NONUF);
double STRID  = SMAX*(1.0 - NRESP);

double U2t = WPARK*PARK + WATAX*ATAX + WBULB*STRID;
dxdt_U2 = (U2t - U2)/TAU_U;
double U1t = 0.80*U2 + WOH*OHSA + WURO*PVRES/PVRMAX*4.0;
dxdt_U1 = (U1t - U1)/TAU_U;
double SCt = WS_IML*(1.0 - NIML) + WS_URO*(1.0 - NONUF)
           + WS_GI*(1.0 - NRESP) + WS_OH*OHSA;
dxdt_SCOPA = (SCt - SCOPA)/TAU_U;

double HAZ = H0*exp(B_U2*U2/10.0 + B_STR*STRID + B_OH*OHSA/5.0 + B_RESP*(1.0 - NRESP));
dxdt_CUMH = HAZ;

$TABLE
// ---- RECOMPUTED at the exact output time --------------------------------
// mrgsolve hoists every `double` declaration from $ODE/$TABLE into one
// shared scope, so the algebra below deliberately RE-ASSIGNS (no `double`)
// the quantities that $ODE also computes.  That matters: the last derivative
// evaluation the solver makes is not generally at the output time, and with
// TAU_MAP = 0.002 h a stale MAP-derived value would corrupt the tilt test.
tod  = fmod(TIME, 24.0);
wake = SIGM(tod - WAKE_H, 8.0)*SIGM(SLEEP_H - tod, 8.0);
double POSo;
if(POSMODE < 0.5)      POSo = 0.0;
else if(POSMODE < 1.5) POSo = 1.0;
else if(POSMODE < 2.5) POSo = UPFRAC*wake + HUTF*(1.0 - wake);
else                   POSo = SIGM(TIME - TILT_T, 120.0);
POSo = POSo*(1.0 - 0.25*COUNTERM);

G_CENT = pow(NMED > 1e-6 ? NMED : 1e-6, GEXP_MED)
       * pow(NIML > 1e-6 ? NIML : 1e-6, GEXP_IML);
G_POST = pow(NMSN > 1e-6 ? NMSN : 1e-6, GEXP_MSN);

SV_  = SV0*pow(VC/VC0 > 0.05 ? VC/VC0 : 0.05, HSV)*(1.0 + KINO*SNA*POSTG);
CO_  = SV_*HR/1000.0;
PPr  = KPP*SV_;
double SBP = MAP + PPr/2.0;
double DBP = MAP - PPr/3.0;

CDMM = 1000.0*DMM/VDMM;
CDRX = 1000.0*DRXC/VDRX;
CATX = 1000.0*ATX/VATX;
CCBD = 1000.0*CBD/VCBD;
CFLU = 1000.0*FLU/VFLU;
CPYR = 1000.0*PYR/VPYR;
CLDB = LDB/VLDB;
double CLDC = 1000.0*LDC/VLD;
CUBQ = 1000.0*UBQ/VUBQ;
CVRD = 1000.0*VRD/VVRD;
CMAB = MAB/VMAB;

AGON = NES + CDMM/KDMM_EQ;
OCC  = AGON/(EC50A1 + AGON);
CARBI= CCBD/(IC50CBD + CCBD);
NETI = ENET*CATX/(EC50ATX + CATX);
DROXNE = KAADC*CDRX*POSTG*(1.0 - CARBI);

MR = ALD/ALD0 + EFLU*CFLU/(EC50FLU + CFLU);
if(MR < 0.05) MR = 0.05;
VOLERR = ECF/ECFTGT;
if(VOLERR < 0.3) VOLERR = 0.3;
UNAV = UNA0*pow(VOLERR, KVOLN)*exp(KPN*(MAP - MAPNAT)/10.0)/pow(MR, KMR);
UVOL = UV0*pow(VOLERR, KVOLW)*exp(KPNW*(MAP - MAPNAT)/10.0)
       *(1.0 - FAQ*AQP2/(KAQ + AQP2));
double NAPL = NAB/(ECF > 1.0 ? ECF : 1.0);

AGON_N = NES;
double IYOHt  = EYOH*(1000.0*YOH/VYOH)/(EC50YOH + 1000.0*YOH/VYOH);
double A2INHt = 1.0/(1.0 + (NES/KA2)*(1.0 - IYOHt));
NEREL = KREL_NE*SNA*NEVES*POSTG*A2INHt;
double NEPL = NEPL0*(NEREL/KNEREF)*(1.0 + KSPILL*NETI);
STRIAT = G_POST*DA/(EC50DA + DA);
STR0   = KDAEND/(EC50DA + KDAEND);
PARK   = PMAX*(1.0 - STRIAT/STR0);
if(PARK < 0.0) PARK = 0.0;
ATAX   = AMAX*(1.0 - NOPC);
PVRES  = PVRMAX*(1.0 - NONUF);
STRID  = SMAX*(1.0 - NRESP);
double U4      = 1.0 + 4.0*(U2/(U2 + 28.0));
double SURV    = exp(-CUMH);
double CBFMARG = MAP - CBFFLOOR;
double YEARS   = TIME/8766.0;

$CAPTURE @annotated
YEARS  : disease years since pathology onset
POSo   : orthostatic load (0 supine, 1 fully upright)
SBP    : systolic blood pressure (mmHg)
DBP    : diastolic blood pressure (mmHg)
CO_    : cardiac output (L/min)
SV_    : stroke volume (mL)
G_CENT : CENTRAL autonomic gain (product of medullary x IML survival)
G_POST : POSTsynaptic striatal gain
OCC    : vascular alpha1 occupancy
NETI   : fractional NET blockade by atomoxetine
CARBI  : fractional peripheral AADC blockade by carbidopa
DROXNE : NE generated from droxidopa (a.u./h)
NEPL   : plasma norepinephrine, spillover-based (pg/mL)
NEREL  : postganglionic NE release rate (a.u./h)
MR     : mineralocorticoid receptor activation (1 = reference)
UNAV   : urinary sodium excretion (mmol/h)
NAPL   : plasma sodium concentration (mmol/L)
UVOL   : urine output (L/h)
PARK   : parkinsonism severity (UMSARS-II points)
ATAX   : ataxia severity (UMSARS-II points)
STRIAT : striatal output (postsynaptic gain x dopamine occupancy)
PVRES  : post-void residual (mL)
STRID  : stridor severity index
U4     : UMSARS Part IV global disability (1-5)
SURV   : survival probability
CBFMARG: cerebral autoregulation margin (mmHg above floor)
CDMM   : desglymidodrine (ng/mL)
CDRX   : droxidopa (ng/mL)
CATX   : atomoxetine (ng/mL)
CLDC   : plasma levodopa (ng/mL)
CLDB   : brain levodopa (ug/mL)
CFLU   : fludrocortisone (ng/mL)
CPYR   : pyridostigmine (ng/mL)
CUBQ   : ubiquinol (ng/mL)
CVRD   : verdiperstat (ng/mL)
CMAB   : anti-alpha-syn antibody (mg/L)
'

mod <- mcode("msa_qsp", msa_code, soloc = tempdir())

## =====================================================================
##  PHENOTYPES — set ONLY by regional vulnerability, never by equations
## =====================================================================
## MSA-P  striatonigral predominant (about 60-70% in Western series)
## MSA-C  olivopontocerebellar predominant (predominant in Japanese series)
## PAF    pure autonomic failure — the postganglionic mirror image
## PD     Parkinson disease — presynaptic loss, POSTsynaptic target intact
pheno <- list(
  MSA_P = c(VULN_SN = 0.70, VULN_MSN = 0.75, VULN_OPC = 0.35,
            VULN_IML = 1.00, VULN_MED = 0.85, VULN_PG = 0.08),
  MSA_C = c(VULN_SN = 0.22, VULN_MSN = 0.20, VULN_OPC = 1.05,
            VULN_IML = 0.95, VULN_MED = 0.80, VULN_PG = 0.08),
  PAF   = c(VULN_SN = 0.05, VULN_MSN = 0.05, VULN_OPC = 0.05,
            VULN_IML = 0.10, VULN_MED = 0.08, VULN_PG = 1.90,
            VULN_ONUF = 0.10, VULN_RESP = 0.05, KND = 0.34),
  PD    = c(VULN_SN = 1.10, VULN_MSN = 0.02, VULN_OPC = 0.05,
            VULN_IML = 0.28, VULN_MED = 0.10, VULN_PG = 0.55,
            VULN_ONUF = 0.12, VULN_RESP = 0.10, KND = 0.30)
)

## An empty list is a legal "no override" — mrgsolve rejects it, so gate it.
pp <- function(m, p) if (length(p)) mrgsolve::param(m, p) else m

## =====================================================================
##  STAGE 1 — natural history (years).  Circadian posture, no drug.
## =====================================================================
msa_history <- function(years = 12, phenotype = "MSA_P", p = list(), delta = 24,
                        tscale = 40) {
  mod %>%
    param(as.list(pheno[[phenotype]])) %>%
    pp(p) %>%
    param(POSMODE = 2, TSCALE = tscale) %>%
    mrgsim(end = years*8766, delta = delta, hmax = 12,
           atol = 1e-7, rtol = 1e-6) %>%
    as.data.frame()
}

## Freeze the full state vector at a chosen disease year -> the initial
## condition for every fast-timescale experiment below.
msa_state_at <- function(year, phenotype = "MSA_P", p = list(), tscale = 40) {
  out <- mod %>%
    param(as.list(pheno[[phenotype]])) %>%
    pp(p) %>%
    param(POSMODE = 2, TSCALE = tscale) %>%
    mrgsim(end = year*8766, delta = 24, hmax = 12, atol = 1e-7, rtol = 1e-6) %>%
    as.data.frame()
  last <- out[nrow(out), ]
  cmts <- as.character(outvars(mod)$cmt)
  st <- as.list(last[, cmts])
  names(st) <- cmts
  st
}

## =====================================================================
##  STAGE 2a — head-up tilt / active stand from a frozen disease state
## =====================================================================
msa_tilt <- function(state, phenotype = "MSA_P", p = list(),
                     dose = NULL, settle = 1.0, dur = 0.75) {
  m <- mod %>%
    param(as.list(pheno[[phenotype]])) %>%
    pp(p) %>%
    param(POSMODE = 3, TILT_T = settle, MEALAMP = 0, TSCALE = 1, INITMAIN = 0) %>%
    init(state)
  if (!is.null(dose)) m <- m %>% ev(dose)
  m %>% mrgsim(end = settle + dur, delta = 1/240, hmax = 0.002,
               atol = 1e-9, rtol = 1e-7, recsort = 3) %>% as.data.frame()
}

## Derive the two diagnostic numbers that define neurogenic OH
tilt_metrics <- function(tt, settle = 1.0) {
  base <- tt[which.min(abs(tt$time - (settle - 0.01))), ]
  at3  <- tt[which.min(abs(tt$time - (settle + 3/60))), ]
  dSBP <- at3$SBP - base$SBP
  dHR  <- at3$HR  - base$HR
  data.frame(SBP_supine = base$SBP, SBP_3min = at3$SBP, dSBP = dSBP,
             HR_supine = base$HR, HR_3min = at3$HR, dHR = dHR,
             HR_SBP_ratio = ifelse(dSBP < 0, dHR/(-dSBP), NA_real_),
             nOH = dSBP <= -20,
             nOH_strict = dSBP <= -30,
             NE_supine = base$NEPL, NE_3min = at3$NEPL,
             NE_rise_pct = 100*(at3$NEPL - base$NEPL)/base$NEPL)
}

## =====================================================================
##  STAGE 2b — 48-h circadian run with real dosing, from a frozen state
## =====================================================================
msa_day <- function(state, phenotype = "MSA_P", p = list(), dose = NULL,
                    days = 2, delta = 0.1) {
  m <- mod %>%
    param(as.list(pheno[[phenotype]])) %>%
    pp(p) %>%
    param(POSMODE = 2, TSCALE = 1, INITMAIN = 0) %>%
    init(state)
  if (!is.null(dose)) m <- m %>% ev(dose)
  m %>% mrgsim(end = days*24, delta = delta, hmax = 0.05,
               atol = 1e-9, rtol = 1e-7, recsort = 3) %>% as.data.frame()
}

## Morning (07:30) standing SBP and the night-time supine peak — the two
## numbers the nocturnal loop trades off against each other.
day_metrics <- function(dd, day = 2) {
  off <- (day - 1)*24
  night <- dd[dd$time >= off + 0 & dd$time <= off + 6, ]
  morn  <- dd[dd$time >= off + 7.4 & dd$time <= off + 7.6, ]
  aft   <- dd[dd$time >= off + 14 & dd$time <= off + 18, ]
  data.frame(SBP_supine_peak_night = max(night$SBP),
             SBP_morning = mean(morn$SBP),
             SBP_afternoon = mean(aft$SBP),
             UNa_night_mmol = sum(night$UNAV)*mean(diff(night$time)),
             UV_night_L = sum(night$UVOL)*mean(diff(night$time)),
             ECF_dawn = morn$ECF[1],
             OHSA = mean(aft$OHSA))
}

## Standard dosing regimens (mg unless noted).  Day-0 = 07:00 reference.
rx <- list(
  midodrine_tid  = ev(amt = 10, cmt = "MIDG", time = 1, ii = 5, addl = 2),
  midodrine_hs   = ev(amt = 10, cmt = "MIDG", time = 22),   # 22:00 — the wrong time
  droxidopa_tid  = ev(amt = 300, cmt = "DRXG", time = 1, ii = 5, addl = 2),
  droxidopa_hi   = ev(amt = 600, cmt = "DRXG", time = 1, ii = 5, addl = 2),
  fludro         = ev(amt = 0.1, cmt = "FLUG", time = 1, ii = 24, addl = 30),
  atomoxetine    = ev(amt = 18, cmt = "ATXG", time = 1, ii = 12, addl = 5),
  pyridostigmine = ev(amt = 60, cmt = "PYRG", time = 1, ii = 5, addl = 2),
  desmopressin   = ev(amt = 2, cmt = "DDA", time = 22, ii = 24, addl = 30),
  water_bolus    = ev(amt = 500, cmt = "WB", time = 0.75),
  levodopa       = ev(amt = 200, cmt = "LDG", time = 1, ii = 5, addl = 2),
  carbidopa      = ev(amt = 50, cmt = "CBDG", time = 1, ii = 5, addl = 2),
  octreotide     = ev(amt = 50, cmt = "OCT", time = 7.5),
  yohimbine      = ev(amt = 5.4, cmt = "YOHG", time = 1, ii = 8, addl = 2),
  ubiquinol      = ev(amt = 1500, cmt = "UBQG", time = 0, ii = 24, addl = 3650),
  verdiperstat   = ev(amt = 600, cmt = "VRDG", time = 0, ii = 12, addl = 3650),
  mab_asyn       = ev(amt = 1200, cmt = "MAB", time = 0, ii = 672, addl = 60)
)

## Chronic (multi-year) regimens for disease-modification scenarios
rx_chronic <- function(what, start_yr = 0, years = 12) {
  t0 <- start_yr*8766
  switch(what,
    ubiquinol    = ev(amt = 1500, cmt = "UBQG", time = t0, ii = 24,
                      addl = floor((years - start_yr)*365)),
    verdiperstat = ev(amt = 600, cmt = "VRDG", time = t0, ii = 12,
                      addl = floor((years - start_yr)*730)),
    mab_asyn     = ev(amt = 1200, cmt = "MAB", time = t0, ii = 672,
                      addl = floor((years - start_yr)*13)),
    levodopa     = ev(amt = 200, cmt = "LDG", time = t0, ii = 5,
                      addl = floor((years - start_yr)*365*3)),
    stop("unknown chronic regimen")
  )
}

msa_history_rx <- function(years = 12, phenotype = "MSA_P", p = list(),
                           dose = NULL, delta = 24, tscale = 40) {
  m <- mod %>%
    param(as.list(pheno[[phenotype]])) %>%
    pp(p) %>%
    param(POSMODE = 2, TSCALE = tscale)
  if (!is.null(dose)) m <- m %>% ev(dose)
  m %>% mrgsim(end = years*8766, delta = delta, hmax = 12,
               atol = 1e-7, rtol = 1e-6, recsort = 3) %>% as.data.frame()
}

## =====================================================================
##  19 SCENARIOS
## =====================================================================
run_scenarios <- function(verbose = TRUE) {
  say <- function(...) if (verbose) cat(...)
  res <- list()

  say("\n================ MSA QSP MODEL — 19 SCENARIOS ================\n")

  ## ---- 1. natural history, MSA-P -----------------------------------
  h_p <- msa_history(12, "MSA_P")
  res$s01 <- h_p
  say(sprintf("\n[1] Natural history MSA-P: UMSARS-II yr4 %.1f -> yr8 %.1f (%.1f pts/yr)\n",
              h_p$U2[which.min(abs(h_p$YEARS - 4))],
              h_p$U2[which.min(abs(h_p$YEARS - 8))],
              (h_p$U2[which.min(abs(h_p$YEARS - 8))] -
               h_p$U2[which.min(abs(h_p$YEARS - 4))])/4))
  say(sprintf("    survival at 8 yr %.2f · G_CENT yr6 %.2f · POSTG yr6 %.2f · NfL yr6 %.0f pg/mL\n",
              h_p$SURV[which.min(abs(h_p$YEARS - 8))],
              h_p$G_CENT[which.min(abs(h_p$YEARS - 6))],
              h_p$POSTG[which.min(abs(h_p$YEARS - 6))],
              h_p$NFL[which.min(abs(h_p$YEARS - 6))]))

  ## ---- 2. natural history, MSA-C ------------------------------------
  h_c <- msa_history(12, "MSA_C")
  res$s02 <- h_c
  say(sprintf("\n[2] Natural history MSA-C: ataxia yr6 %.1f vs parkinsonism yr6 %.1f (MSA-P: %.1f / %.1f)\n",
              h_c$ATAX[which.min(abs(h_c$YEARS - 6))],
              h_c$PARK[which.min(abs(h_c$YEARS - 6))],
              h_p$ATAX[which.min(abs(h_p$YEARS - 6))],
              h_p$PARK[which.min(abs(h_p$YEARS - 6))]))

  ## ---- 3. tilt test as a stage readout ------------------------------
  say("\n[3] Head-up tilt at successive disease years (MSA-P, no drug):\n")
  tilt_tab <- do.call(rbind, lapply(c(1, 3, 5, 7, 9), function(y) {
    st <- msa_state_at(y, "MSA_P")
    m  <- tilt_metrics(msa_tilt(st, "MSA_P"))
    cbind(year = y, m)
  }))
  res$s03 <- tilt_tab
  print(tilt_tab[, c("year", "SBP_supine", "SBP_3min", "dSBP", "dHR",
                     "HR_SBP_ratio", "NE_rise_pct")], digits = 3)

  ## ---- 4. midodrine 10 mg t.i.d. -----------------------------------
  st5 <- msa_state_at(7, "MSA_P")
  t_base <- tilt_metrics(msa_tilt(st5, "MSA_P"))
  t_mido <- tilt_metrics(msa_tilt(st5, "MSA_P",
                                  dose = ev(amt = 10, cmt = "MIDG", time = 0)))
  res$s04 <- rbind(cbind(rx = "none", t_base), cbind(rx = "midodrine 10 mg", t_mido))
  say(sprintf("\n[4] Midodrine 10 mg (yr7): standing SBP %.0f -> %.0f mmHg (%+.0f)\n",
              t_base$SBP_3min, t_mido$SBP_3min, t_mido$SBP_3min - t_base$SBP_3min))

  ## ---- 5. droxidopa 300 mg -----------------------------------------
  t_drox <- tilt_metrics(msa_tilt(st5, "MSA_P",
                                  dose = ev(amt = 300, cmt = "DRXG", time = 0)))
  res$s05 <- t_drox
  say(sprintf("[5] Droxidopa 300 mg  (yr7): standing SBP %.0f -> %.0f mmHg (%+.0f)\n",
              t_base$SBP_3min, t_drox$SBP_3min, t_drox$SBP_3min - t_base$SBP_3min))

  ## ---- 6. droxidopa + carbidopa: shared-enzyme antagonism -----------
  t_dc <- tilt_metrics(msa_tilt(st5, "MSA_P",
                                dose = c(ev(amt = 300, cmt = "DRXG", time = 0),
                                         ev(amt = 50, cmt = "CBDG", time = 0))))
  res$s06 <- t_dc
  say(sprintf("[6] Droxidopa + carbidopa 50 mg: %+.0f mmHg — %.0f%% of the droxidopa effect lost\n",
              t_dc$SBP_3min - t_base$SBP_3min,
              100*(1 - (t_dc$SBP_3min - t_base$SBP_3min)/
                     max(t_drox$SBP_3min - t_base$SBP_3min, 1e-6))))

  ## ---- 7/8. atomoxetine: MSA vs the PAF mirror image ----------------
  st5_paf <- msa_state_at(5, "PAF")
  t_paf_b <- tilt_metrics(msa_tilt(st5_paf, "PAF"))
  t_msa_a <- tilt_metrics(msa_tilt(st5, "MSA_P",
                                   dose = ev(amt = 18, cmt = "ATXG", time = 0)))
  t_paf_a <- tilt_metrics(msa_tilt(st5_paf, "PAF",
                                   dose = ev(amt = 18, cmt = "ATXG", time = 0)))
  res$s07 <- data.frame(
    phenotype = c("MSA-P", "MSA-P", "PAF", "PAF"),
    atomoxetine = c("no", "yes", "no", "yes"),
    SBP_3min = c(t_base$SBP_3min, t_msa_a$SBP_3min,
                 t_paf_b$SBP_3min, t_paf_a$SBP_3min),
    POSTG = c(st5$POSTG, st5$POSTG, st5_paf$POSTG, st5_paf$POSTG))
  say(sprintf("\n[7/8] Atomoxetine 18 mg — the POSTG-dependence test:\n"))
  say(sprintf("      MSA-P (POSTG %.2f): %+.0f mmHg   |   PAF (POSTG %.2f): %+.0f mmHg\n",
              st5$POSTG, t_msa_a$SBP_3min - t_base$SBP_3min,
              st5_paf$POSTG, t_paf_a$SBP_3min - t_paf_b$SBP_3min))

  ## ---- 9. non-pharmacological bundle --------------------------------
  t_bundle <- tilt_metrics(msa_tilt(st5, "MSA_P",
                                    p = list(BINDER = 1, COUNTERM = 1),
                                    dose = ev(amt = 500, cmt = "WB", time = 0)))
  res$s09 <- t_bundle
  say(sprintf("\n[9] Binder + counter-manoeuvres + 500 mL water bolus: %+.0f mmHg\n",
              t_bundle$SBP_3min - t_base$SBP_3min))

  ## ---- 10-12. THE NOCTURNAL LOOP ------------------------------------
  say("\n[10-12] The nocturnal loop — same daytime pressor exposure, different nights:\n")
  w <- function(p = list(), dose = NULL)
    day_metrics(msa_day(st5, "MSA_P", p = p, dose = dose, days = 7, delta = 0.2), day = 7)
  night_tab <- rbind(
    cbind(arm = "untreated", w()),
    cbind(arm = "midodrine t.i.d. (correct timing)", w(dose = rx$midodrine_tid)),
    cbind(arm = "midodrine at bedtime (wrong timing)",
          w(dose = ev(amt = 10, cmt = "MIDG", time = 22, ii = 24, addl = 7))),
    cbind(arm = "head-up tilt sleeping", w(p = list(HUTF = 0.22))),
    cbind(arm = "bedtime desmopressin", w(dose = rx$desmopressin)),
    cbind(arm = "head-up tilt + desmopressin",
          w(p = list(HUTF = 0.22), dose = rx$desmopressin)))
  res$s10 <- night_tab
  print(night_tab, digits = 4, row.names = FALSE)

  ## ---- 13. pyridostigmine: reflex-only amplification ----------------
  t_pyr <- tilt_metrics(msa_tilt(st5, "MSA_P",
                                 dose = ev(amt = 60, cmt = "PYRG", time = 0)))
  d_pyr <- day_metrics(msa_day(st5, "MSA_P", dose = rx$pyridostigmine))
  d_mid <- day_metrics(msa_day(st5, "MSA_P", dose = rx$midodrine_tid))
  res$s13 <- list(tilt = t_pyr, day = d_pyr)
  say(sprintf("\n[13] Pyridostigmine 60 mg: standing %+.0f mmHg, night supine peak %.0f mmHg\n",
              t_pyr$SBP_3min - t_base$SBP_3min, d_pyr$SBP_supine_peak_night))
  say(sprintf("     (midodrine for comparison: night supine peak %.0f mmHg)\n",
              d_mid$SBP_supine_peak_night))

  ## ---- 14. fludrocortisone (volume arm) -----------------------------
  d_flu <- day_metrics(msa_day(st5, "MSA_P", dose = rx$fludro, days = 2))
  h_flu <- msa_day(st5, "MSA_P", dose = rx$fludro, days = 14, delta = 1)
  res$s14 <- h_flu
  say(sprintf("\n[14] Fludrocortisone 0.1 mg: ECF %.2f -> %.2f L over 14 d, afternoon SBP %+.0f mmHg\n",
              h_flu$ECF[1], h_flu$ECF[nrow(h_flu)],
              mean(tail(h_flu$SBP[h_flu$POSo > 0.4], 40)) -
                mean(head(h_flu$SBP[h_flu$POSo > 0.4], 40))))

  ## ---- 15/16. levodopa: MSA-P vs PD, and the waning response --------
  say("\n[15/16] Levodopa 200 mg t.i.d. — postsynaptic gain decides the answer:\n")
  ld_tab <- do.call(rbind, lapply(c(2, 4, 6, 8), function(y) {
    st_m <- msa_state_at(y, "MSA_P"); st_d <- msa_state_at(y, "PD")
    dm <- msa_day(st_m, "MSA_P", dose = rx$levodopa, days = 1, delta = 0.25)
    dd <- msa_day(st_d, "PD",    dose = rx$levodopa, days = 1, delta = 0.25)
    bm <- msa_day(st_m, "MSA_P", days = 1, delta = 0.25)
    bd <- msa_day(st_d, "PD",    days = 1, delta = 0.25)
    data.frame(year = y,
               MSA_NMSN = round(st_m$NMSN, 3), PD_NMSN = round(st_d$NMSN, 3),
               MSA_park_off = round(mean(bm$PARK), 1),
               MSA_park_on  = round(mean(dm$PARK), 1),
               MSA_gain_pct = round(100*(mean(bm$PARK) - mean(dm$PARK))/
                                      max(mean(bm$PARK), 1e-6), 1),
               PD_gain_pct  = round(100*(mean(bd$PARK) - mean(dd$PARK))/
                                      max(mean(bd$PARK), 1e-6), 1))
  }))
  res$s15 <- ld_tab
  print(ld_tab, row.names = FALSE)

  ## ---- 17. high-dose ubiquinol -------------------------------------
  h_ubq <- msa_history_rx(12, "MSA_P", dose = rx_chronic("ubiquinol", 3, 12))
  h_ubq_coq2 <- msa_history_rx(12, "MSA_P", p = list(COQ2F = 0.55),
                               dose = rx_chronic("ubiquinol", 3, 12))
  h_coq2 <- msa_history(12, "MSA_P", p = list(COQ2F = 0.55))
  res$s17 <- list(wt = h_ubq, coq2 = h_ubq_coq2, coq2_untreated = h_coq2)
  say(sprintf("\n[17] Ubiquinol 1500 mg/d from yr3 — UMSARS-II at yr9:\n"))
  say(sprintf("     COQ2 wild-type: %.1f (untreated %.1f) | COQ2-deficient: %.1f (untreated %.1f)\n",
              h_ubq$U2[which.min(abs(h_ubq$YEARS - 9))],
              h_p$U2[which.min(abs(h_p$YEARS - 9))],
              h_ubq_coq2$U2[which.min(abs(h_ubq_coq2$YEARS - 9))],
              h_coq2$U2[which.min(abs(h_coq2$YEARS - 9))]))

  ## ---- 18. verdiperstat — must reproduce a NULL ---------------------
  h_vrd <- msa_history_rx(12, "MSA_P", dose = rx_chronic("verdiperstat", 3, 12))
  res$s18 <- h_vrd
  say(sprintf("\n[18] Verdiperstat 600 mg b.i.d. from yr3 (M-STAR reproduction):\n"))
  say(sprintf("     UMSARS-II at yr4 %.1f vs untreated %.1f (48-week delta %.2f pts)\n",
              h_vrd$U2[which.min(abs(h_vrd$YEARS - 3.92))],
              h_p$U2[which.min(abs(h_p$YEARS - 3.92))],
              h_p$U2[which.min(abs(h_p$YEARS - 3.92))] -
                h_vrd$U2[which.min(abs(h_vrd$YEARS - 3.92))]))

  ## ---- 19. anti-alpha-syn antibody: early vs late -------------------
  h_mab_e <- msa_history_rx(12, "MSA_P", dose = rx_chronic("mab_asyn", 1, 12))
  h_mab_l <- msa_history_rx(12, "MSA_P", dose = rx_chronic("mab_asyn", 6, 12))
  res$s19 <- list(early = h_mab_e, late = h_mab_l)
  g <- function(d, y) d$U2[which.min(abs(d$YEARS - y))]
  say(sprintf("\n[19] Anti-alpha-syn mAb — the therapeutic window is EMERGENT:\n"))
  say(sprintf("     UMSARS-II at yr10: untreated %.1f | start yr1 %.1f | start yr6 %.1f\n",
              g(h_p, 10), g(h_mab_e, 10), g(h_mab_l, 10)))
  say(sprintf("     survival at yr10:  untreated %.2f | start yr1 %.2f | start yr6 %.2f\n",
              h_p$SURV[which.min(abs(h_p$YEARS - 10))],
              h_mab_e$SURV[which.min(abs(h_mab_e$YEARS - 10))],
              h_mab_l$SURV[which.min(abs(h_mab_l$YEARS - 10))]))

  say("\n==============================================================\n")
  invisible(res)
}

## =====================================================================
##  SELF-TEST — hard assertions on the behaviours the model must show.
##  Every threshold below is a clinical fact, not a fitted output: if an
##  edit breaks one of them the model has stopped describing MSA.
##  Reference stage for drug tests is disease YEAR 7 ("established MSA"),
##  which is where the model puts a patient meeting the nOH definition.
## =====================================================================
msa_selftest <- function() {
  ok <- TRUE
  .stc <- new.env(parent = emptyenv())
  state <- function(year, phenotype = "MSA_P") {
    k <- paste0(phenotype, "|", year)
    if (is.null(.stc[[k]])) .stc[[k]] <- msa_state_at(year, phenotype)
    .stc[[k]]
  }
  .tlc <- new.env(parent = emptyenv())
  tm <- function(year, ph = "MSA_P") {
    k <- paste0(ph, "|", year)
    if (is.null(.tlc[[k]])) .tlc[[k]] <- tilt_metrics(msa_tilt(state(year, ph), ph))
    .tlc[[k]]
  }
  gain <- function(year, ph, dose) {
    tilt_metrics(msa_tilt(state(year, ph), ph, dose = dose))$SBP_3min - tm(year, ph)$SBP_3min
  }
  chk <- function(name, cond, detail = "") {
    cond <- isTRUE(cond)
    cat(sprintf("  [%s] %s%s\n", if (cond) "PASS" else "FAIL", name,
                if (nzchar(detail)) paste0(" — ", detail) else ""))
    if (!cond) ok <<- FALSE
    invisible(cond)
  }
  cat("\n--- self-test ---\n")

  ## ---------------- 1. the healthy control -------------------------
  healthy <- mod %>% param(KND = 0, KSYN_A = 0, KNUC = 0) %>%
    param(POSMODE = 2, TSCALE = 1) %>%
    mrgsim(end = 24*15, delta = 1, hmax = 0.05) %>% as.data.frame()
  st_h <- as.list(healthy[nrow(healthy), as.character(outvars(mod)$cmt)])
  th <- tilt_metrics(msa_tilt(st_h))
  chk("healthy control has no nOH", th$dSBP > -20, sprintf("dSBP %+.1f mmHg", th$dSBP))
  chk("healthy control mounts a brisk tachycardia", th$dHR > 10, sprintf("dHR %+.1f bpm", th$dHR))
  chk("healthy supine plasma NE is 150-400 pg/mL",
      th$NE_supine > 150 && th$NE_supine < 400, sprintf("%.0f pg/mL", th$NE_supine))
  chk("healthy control mounts a plasma NE rise on standing", th$NE_rise_pct > 10,
      sprintf("%+.0f%%", th$NE_rise_pct))
  chk("healthy supine SBP/DBP are physiological",
      th$SBP_supine > 105 && th$SBP_supine < 135, sprintf("SBP %.0f mmHg", th$SBP_supine))

  ## ---------------- 2. the MSA fingerprint -------------------------
  s5 <- state(5); s7 <- state(7); s9 <- state(9)
  t5 <- tm(5); t7 <- tm(7); t9 <- tm(9)
  gc <- function(st) st$NMED^2.20 * st$NIML^1.60
  chk("MSA-P year 7 meets the nOH definition", t7$dSBP <= -20,
      sprintf("dSBP %+.1f mmHg", t7$dSBP))
  chk("MSA-P shows the neurogenic dHR/dSBP < 0.5 signature",
      !is.na(t7$HR_SBP_ratio) && t7$HR_SBP_ratio < 0.5, sprintf("ratio %.2f", t7$HR_SBP_ratio))
  chk("orthostatic fall worsens monotonically with disease stage",
      t5$dSBP > t7$dSBP && t7$dSBP > t9$dSBP,
      sprintf("yr5 %+.0f -> yr7 %+.0f -> yr9 %+.0f mmHg", t5$dSBP, t7$dSBP, t9$dSBP))
  chk("established MSA keeps a NORMAL supine plasma NE",
      t5$NE_supine > 100 && t5$NE_supine < 400, sprintf("yr5 %.0f pg/mL", t5$NE_supine))
  chk("late MSA has a BLUNTED upright NE increment",
      t9$NE_rise_pct < th$NE_rise_pct/3,
      sprintf("yr9 %+.0f%% vs healthy %+.0f%%", t9$NE_rise_pct, th$NE_rise_pct))
  chk("postganglionic terminals are comparatively SPARED in MSA", s7$POSTG > 0.85,
      sprintf("POSTG %.2f", s7$POSTG))
  chk("central gain is nearly abolished in MSA by year 7", gc(s7) < 0.10,
      sprintf("G_CENT %.3f", gc(s7)))
  chk("alpha1 receptors are UP-regulated as drive falls", s7$A1R > 1.3*state(2)$A1R,
      sprintf("A1R %.2f (yr7) vs %.2f (yr2)", s7$A1R, state(2)$A1R))

  ## ---------------- 3. the PAF mirror image ------------------------
  p5 <- state(5, "PAF"); tp5 <- tm(5, "PAF")
  chk("PAF phenotype loses the POSTGANGLIONIC terminal instead", p5$POSTG < 0.15,
      sprintf("POSTG %.3f (central gain %.2f, i.e. INTACT)", p5$POSTG, gc(p5)))
  chk("PAF has nOH of comparable severity (so the split is not just severity)",
      tp5$dSBP <= -20, sprintf("dSBP %+.1f vs MSA yr7 %+.1f mmHg", tp5$dSBP, t7$dSBP))
  chk("PAF has a LOW supine plasma NE (the bedside discriminator)",
      tp5$NE_supine < 0.5*t5$NE_supine,
      sprintf("PAF %.0f vs MSA %.0f pg/mL", tp5$NE_supine, t5$NE_supine))

  ## ---------------- 4. drug selectivity by lesion site ------------
  mido7 <- gain(7, "MSA_P", ev(amt = 10, cmt = "MIDG", time = 0))
  midoP <- gain(5, "PAF",   ev(amt = 10, cmt = "MIDG", time = 0))
  chk("midodrine 10 mg raises standing SBP by 15-60 mmHg in established MSA",
      mido7 > 15 && mido7 < 60, sprintf("%+.1f mmHg", mido7))
  chk("midodrine works in BOTH phenotypes (it is downstream of every lesion)",
      midoP > 10, sprintf("PAF %+.1f mmHg", midoP))

  drox7 <- gain(7, "MSA_P", ev(amt = 300, cmt = "DRXG", time = 0))
  droxP <- gain(5, "PAF",   ev(amt = 300, cmt = "DRXG", time = 0))
  droxc <- gain(7, "MSA_P", c(ev(amt = 300, cmt = "DRXG", time = 0),
                              ev(amt = 50, cmt = "CBDG", time = 0)))
  chk("droxidopa raises standing SBP in MSA", drox7 > 10, sprintf("%+.1f mmHg", drox7))
  chk("droxidopa REQUIRES the postganglionic terminal (fails in PAF)",
      droxP < drox7/3, sprintf("PAF %+.1f vs MSA %+.1f mmHg", droxP, drox7))
  chk("carbidopa antagonises droxidopa via the shared AADC term (>20% lost)",
      droxc < 0.8*drox7, sprintf("%+.1f vs %+.1f mmHg (%.0f%% lost)",
                                 droxc, drox7, 100*(1 - droxc/drox7)))

  atx7 <- gain(7, "MSA_P", ev(amt = 18, cmt = "ATXG", time = 0))
  atxP <- gain(5, "PAF",   ev(amt = 18, cmt = "ATXG", time = 0))
  chk("atomoxetine works in MSA (central lesion, intact terminals)", atx7 > 12,
      sprintf("%+.1f mmHg", atx7))
  chk("atomoxetine FAILS in PAF at matched OH severity (<1/4 of the MSA effect)",
      atxP < atx7/4, sprintf("PAF %+.1f vs MSA %+.1f mmHg", atxP, atx7))

  pyr5 <- gain(5, "MSA_P", ev(amt = 60, cmt = "PYRG", time = 0))
  pyr9 <- gain(9, "MSA_P", ev(amt = 60, cmt = "PYRG", time = 0))
  pyrP <- gain(5, "PAF",   ev(amt = 60, cmt = "PYRG", time = 0))
  chk("pyridostigmine works while a baroreflex remains", pyr5 > 2,
      sprintf("yr5 %+.1f mmHg", pyr5))
  chk("pyridostigmine fades as the reflex it amplifies dies", pyr9 < 1 && pyr9 < pyr5,
      sprintf("yr9 %+.1f mmHg", pyr9))
  chk("pyridostigmine does nothing when the effector is gone (PAF)", pyrP < 1.5,
      sprintf("PAF %+.1f mmHg", pyrP))

  g2 <- gain(2, "MSA_P", ev(amt = 5, cmt = "MIDG", time = 0))
  g7 <- gain(7, "MSA_P", ev(amt = 5, cmt = "MIDG", time = 0))
  chk("denervation supersensitivity: 5 mg midodrine buys more at yr7 than yr2",
      g7 > 2*g2, sprintf("%+.1f (yr7) vs %+.1f mmHg (yr2)", g7, g2))

  ## ---------------- 5. THE NOCTURNAL LOOP (7-day runs) ------------
  dm <- function(p = list(), dose = NULL)
    day_metrics(msa_day(s7, "MSA_P", p = p, dose = dose, days = 7, delta = 0.2), day = 7)
  d_un <- dm()
  ## t = 0 is midnight in these runs, so a BEDTIME dose is time = 22.
  MIDO_HS <- ev(amt = 10, cmt = "MIDG", time = 22, ii = 24, addl = 7)
  d_hs <- dm(dose = MIDO_HS)
  ## Night 1 is where the natriuretic TRANSIENT lives.  By day 7 the patient
  ## has settled at a lower volume, at which point sodium balance is restored
  ## and nightly excretion is back to intake -- so the loop does not bleed
  ## sodium forever, it moves the volume SETPOINT down.  Both facts are worth
  ## asserting, and asserting only the day-7 rate would look like a failure.
  ## the first night AFTER a 22:00 dose is the 24-30 h window, i.e. day = 2
  n1 <- function(p = list(), dose = NULL)
    day_metrics(msa_day(s7, "MSA_P", p = p, dose = dose, days = 3, delta = 0.2), day = 2)
  n1_un <- n1(); n1_hs <- n1(dose = MIDO_HS)
  d_hu <- dm(p = list(HUTF = 0.22))
  d_dd <- dm(dose = rx$desmopressin)
  d_bo <- dm(p = list(HUTF = 0.22), dose = rx$desmopressin)
  chk("bedtime pressor raises the nocturnal supine peak",
      d_hs$SBP_supine_peak_night > d_un$SBP_supine_peak_night + 3,
      sprintf("%.1f vs %.1f mmHg", d_hs$SBP_supine_peak_night, d_un$SBP_supine_peak_night))
  chk("bedtime pressor increases sodium loss on the FIRST night (the transient)",
      n1_hs$UNa_night_mmol > n1_un$UNa_night_mmol,
      sprintf("night 1: %.1f vs %.1f mmol", n1_hs$UNa_night_mmol, n1_un$UNa_night_mmol))
  ## A pressor that is still on board at 02:00 keeps the nocturnal pressure
  ## term elevated EVERY night, so natriuresis does not settle back to intake
  ## and the volume deficit keeps widening.  That is the loop, and it is why
  ## the >4 h-before-recumbency rule is a rule and not a nicety.
  chk("a bedtime pressor keeps overnight sodium loss elevated through day 7",
      d_hs$UNa_night_mmol > d_un$UNa_night_mmol + 5 &&
        d_hs$ECF_dawn < d_un$ECF_dawn - 0.3,
      sprintf("day 7 Na %.1f vs %.1f mmol, dawn ECF %.2f vs %.2f L",
              d_hs$UNa_night_mmol, d_un$UNa_night_mmol, d_hs$ECF_dawn, d_un$ECF_dawn))
  chk("bedtime pressor leaves the patient volume-depleted by dawn",
      d_hs$ECF_dawn < d_un$ECF_dawn,
      sprintf("%.3f vs %.3f L", d_hs$ECF_dawn, d_un$ECF_dawn))
  chk("head-up tilt sleeping LOWERS the nocturnal supine peak",
      d_hu$SBP_supine_peak_night < d_un$SBP_supine_peak_night,
      sprintf("%.1f vs %.1f mmHg", d_hu$SBP_supine_peak_night, d_un$SBP_supine_peak_night))
  chk("head-up tilt sleeping SPARES volume overnight (no daytime drug at all)",
      d_hu$ECF_dawn > d_un$ECF_dawn,
      sprintf("%.3f vs %.3f L", d_hu$ECF_dawn, d_un$ECF_dawn))
  chk("bedtime desmopressin cuts overnight urine output",
      d_dd$UV_night_L < d_un$UV_night_L,
      sprintf("%.3f vs %.3f L", d_dd$UV_night_L, d_un$UV_night_L))
  chk("bedtime desmopressin raises MORNING standing SBP",
      d_dd$SBP_morning > d_un$SBP_morning,
      sprintf("%.2f vs %.2f mmHg", d_dd$SBP_morning, d_un$SBP_morning))
  chk("head-up tilt + desmopressin is additive on dawn volume",
      d_bo$ECF_dawn > d_dd$ECF_dawn && d_bo$ECF_dawn > d_hu$ECF_dawn,
      sprintf("%.3f L vs %.3f / %.3f L", d_bo$ECF_dawn, d_dd$ECF_dawn, d_hu$ECF_dawn))
  d_py <- dm(dose = rx$pyridostigmine)
  d_mi <- dm(dose = rx$midodrine_tid)
  ## Correctly timed midodrine (last dose ~11:00) is cleared by bedtime, so it
  ## is NOT the arm that causes supine hypertension -- that is the whole point
  ## of the >4 h rule.  The meaningful contrast for pyridostigmine, whose
  ## amplification is gated on baroreflex traffic and therefore silent supine,
  ## is against a pressor that IS on board at night.
  chk("pyridostigmine spares the supine peak relative to a night-time pressor",
      d_py$SBP_supine_peak_night < d_hs$SBP_supine_peak_night - 3,
      sprintf("%.1f vs %.1f mmHg (correctly timed midodrine: %.1f)",
              d_py$SBP_supine_peak_night, d_hs$SBP_supine_peak_night,
              d_mi$SBP_supine_peak_night))
  chk("correctly timed midodrine does NOT raise the nocturnal supine peak",
      d_mi$SBP_supine_peak_night < d_hs$SBP_supine_peak_night - 3,
      sprintf("%.1f vs %.1f mmHg", d_mi$SBP_supine_peak_night, d_hs$SBP_supine_peak_night))

  ## ---------------- 6. levodopa: postsynaptic gain decides --------
  ldgain <- function(year, ph) {
    st <- state(year, ph)
    off <- msa_day(st, ph, days = 1, delta = 0.25)
    on  <- msa_day(st, ph, dose = rx$levodopa, days = 1, delta = 0.25)
    100*(mean(off$PARK) - mean(on$PARK))/max(mean(off$PARK), 1e-6)
  }
  gm7 <- ldgain(7, "MSA_P"); gd7 <- ldgain(7, "PD"); gm2 <- ldgain(2, "MSA_P")
  chk("levodopa buys much less in MSA-P than in PD at equal brain exposure",
      gm7 < gd7/2, sprintf("MSA %.0f%% vs PD %.0f%%", gm7, gd7))
  chk("the MSA-P levodopa response WANES with disease stage", gm2 > gm7,
      sprintf("%.0f%% (yr2) -> %.0f%% (yr7)", gm2, gm7))

  ## ---------------- 7. natural history and survival ---------------
  h <- msa_history(14, "MSA_P", delta = 48)
  at <- function(y, col) h[[col]][which.min(abs(h$YEARS - y))]
  r <- (at(8, "U2") - at(4, "U2"))/4
  chk("UMSARS-II progresses 3-11 points/year", r > 3 && r < 11, sprintf("%.1f pts/yr", r))
  sv <- approx(h$SURV, h$YEARS, xout = 0.5)$y
  chk("median survival is 6-11 years from pathology onset",
      !is.na(sv) && sv > 6 && sv < 11, sprintf("%.1f yr", sv))
  chk("post-void residual exceeds 100 mL by year 5", at(5, "PVRES") > 100,
      sprintf("%.0f mL", at(5, "PVRES")))
  chk("plasma NfL rises to at least twice the healthy baseline",
      max(h$NFL) > 2*h$NFL[1], sprintf("%.0f -> %.0f pg/mL", h$NFL[1], max(h$NFL)))
  chk("stridor emerges and contributes to the hazard",
      at(8, "STRID") > 0.2, sprintf("index %.2f at yr8", at(8, "STRID")))
  hc <- msa_history(8, "MSA_C", delta = 48)
  chk("MSA-C is ataxia-dominant while MSA-P is parkinsonism-dominant",
      tail(hc$ATAX, 1) > tail(hc$PARK, 1) && at(8, "PARK") > at(8, "ATAX"),
      sprintf("MSA-C %.0f/%.0f vs MSA-P %.0f/%.0f (ataxia/park)",
              tail(hc$ATAX, 1), tail(hc$PARK, 1), at(8, "ATAX"), at(8, "PARK")))

  ## ---------------- 8. disease-modification arms ------------------
  hv <- msa_history_rx(14, "MSA_P", dose = rx_chronic("verdiperstat", 3, 14), delta = 48)
  d48 <- at(3.92, "U2") - hv$U2[which.min(abs(hv$YEARS - 3.92))]
  chk("a myeloperoxidase inhibitor started at diagnosis is clinically NULL at 48 weeks",
      abs(d48) < 1.5, sprintf("UMSARS-II delta %.2f pts", d48))
  he <- msa_history_rx(14, "MSA_P", dose = rx_chronic("mab_asyn", 1, 14), delta = 48)
  hl <- msa_history_rx(14, "MSA_P", dose = rx_chronic("mab_asyn", 6, 14), delta = 48)
  ge <- at(10, "U2") - he$U2[which.min(abs(he$YEARS - 10))]
  gl <- at(10, "U2") - hl$U2[which.min(abs(hl$YEARS - 10))]
  chk("anti-alpha-syn antibody started EARLY beats the same drug started late",
      ge > gl + 1, sprintf("%.1f vs %.1f UMSARS-II points saved", ge, gl))

  ## ---------------- 9. numerics and conservation ------------------
  q40 <- msa_history(4, "MSA_P", delta = 168, tscale = 40)
  q01 <- msa_history(4, "MSA_P", delta = 168, tscale = 1)
  relv <- function(v) abs(tail(q40[[v]], 1) - tail(q01[[v]], 1))/max(abs(tail(q01[[v]], 1)), 1e-9)
  worst <- max(sapply(c("GCI", "MYE", "G_CENT", "POSTG", "NMSN", "NOPC", "NONUF",
                        "U1", "U2", "SCOPA", "ECF", "NAB", "NFL", "SURV"), relv))
  chk("TSCALE=40 reproduces the stiff TSCALE=1 solve to <1%", worst < 0.01,
      sprintf("worst relative difference %.1e", worst))
  chk("plasma sodium stays physiological across 14 years",
      all(h$NAPL > 132) && all(h$NAPL < 148),
      sprintf("range %.1f-%.1f mmol/L", min(h$NAPL), max(h$NAPL)))
  chk("extracellular volume stays physiological across 14 years",
      all(h$ECF > 10) && all(h$ECF < 20), sprintf("range %.1f-%.1f L", min(h$ECF), max(h$ECF)))
  chk("no state goes negative or non-finite anywhere in the 14-year run",
      all(is.finite(as.matrix(h))) && min(h$MYE, h$NIML, h$POSTG, h$NAB, h$ECF, h$CQ) > 0,
      "all finite and positive")

  cat(sprintf("--- self-test %s ---\n\n", if (ok) "PASSED" else "FAILED"))
  invisible(ok)
}

if (!interactive() && identical(environment(), globalenv())) {
  ok <- msa_selftest()
  run_scenarios()
  if (!isTRUE(ok)) quit(status = 1)
}
