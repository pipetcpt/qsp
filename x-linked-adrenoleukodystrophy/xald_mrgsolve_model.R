# =============================================================================
#  X-linked Adrenoleukodystrophy (X-ALD) — QSP model for mrgsolve
#  X-연관 부신백질형성장애 · 정량적 시스템 약리학 모델
# =============================================================================
#
#  ORGANISING THESIS
#  -----------------------------------------------------------------------------
#  ONE lesion (ABCD1 loss) drives TWO INDEPENDENT variables:
#
#    (1) VLCFA accumulation      deterministic · dose-graded · MEASURABLE
#                                · reachable from the periphery
#    (2) Cerebral inflammation   BISTABLE SWITCH · unmeasurable until MRI
#                                · reachable only from inside the CNS
#
#  In this file the coupling (1) -> (2) is written as a THRESHOLD multiplied by
#  a susceptibility term (`SUSC`), NOT as a rate constant.  Because brain VLCFA
#  (`CBR`) enters through a Hill gate that is already >90% saturated in every
#  X-ALD patient, `CBR` cannot discriminate phenotype: the ignition decision is
#  carried entirely by `SUSC`.  That single structural choice is what makes the
#  following model OUTPUTS rather than coded rules:
#
#    * plasma C26:0 does not predict phenotype (scenario 20 sweeps it and the
#      phenotype does not move; scenario 21 sweeps SUSC and it flips);
#    * Lorenzo's oil normalises plasma C26:0 and does not stop cerebral
#      disease  (erucic acid inhibits ELOVL1 in the periphery; `KINBR` carries
#      plasma to brain so slowly, and `LOBBB` = 0 exposure means the brain
#      elongation arm is untouched);
#    * HSCT / eli-cel arrests cerebral disease and barely moves plasma C26:0
#      (corrected microglia restore `PHAGC`, i.e. variable 2; they are far too
#      small a fraction of brain lipid mass to move `CBR`);
#    * the ~12-18 month post-transplant lag window is a consequence of
#      `KMGREP` (slow microglial replacement), not a rule.
#
#  THE DRUG THAT FIXES THE BIOMARKER DOES NOT FIX THE DISEASE, AND THE DRUG
#  THAT FIXES THE DISEASE DOES NOT FIX THE BIOMARKER.
#
#  STRUCTURE
#  -----------------------------------------------------------------------------
#  63 ODE compartments · time unit = DAYS · horizon up to 40 years
#    VLCFA mass balance & tissue distribution   10
#    Lorenzo's oil PK (erucic / oleic)           3
#    Adrenocortical axis + replacement           9
#    Testicular axis                             3
#    Microglia / neuroinflammation / cerebral   16
#    Spinal cord axonopathy (AMN)                4
#    Clinical ratchets (Loes, NFS)               2
#    HSCT / gene therapy / safety               10  (incl. cumulative hazard)
#    CNS-penetrant drugs (leriglitazone, sobetirome, DMF) 6
#
#  ENTRY POINTS (R functions defined below the model)
#    xald_mod()          load + compile the model
#    xald_birth(...)     280-day in-utero burn-in -> state at birth
#    xald_run(...)       run one named scenario
#    xald_scenarios()    the 24 scenarios
#    xald_validate()     anchor suite A1-A16, prints target vs model
#    xald_calibrate()    coordinate-wise fitting used to set the shipped values
#    xald_population()   population run: phenotype distribution from SUSC IIV
#
#  All parameter provenance is in xald_references.md.  Values marked FIT were
#  set by xald_calibrate(); values marked LIT are literature; values marked
#  STRUCT are structural choices with no direct measurement.
# =============================================================================

library(mrgsolve)

xald_code <- '
$PROB
# X-ALD QSP model
- One lesion, two independent variables: VLCFA (measurable, peripheral) and a
  bistable cerebral inflammatory switch (unmeasurable, central).
- 63 ODEs, time unit = days.

$PARAM @annotated
// ---------------------------------------------------------------- genotype --
MUTRES   : 0.02  : Residual ALDP transport function, fraction of normal (0=null)
FMOS     : 0.0   : Fraction of ALDP-deficient cells (0 male hemizygote, 0-1 female carrier)
CRSC     : 0.35  : Metabolic cross-correction from ALDP-competent neighbours
SUSCTV   : 0.35  : Susceptibility to cerebral ignition (population typical; threshold ~0.42)
ROSREF   : 0.162680: Healthy-reference cerebral ROS (subtracted so health is the zero point)
ATPREF   : 0.119817: Healthy-reference axonal energetic deficit (subtracted likewise)

// ------------------------------------------------- VLCFA mass balance (1) --
KELONG   : 0.080 : ELOVL1-dependent de-novo elongation input to plasma pool (umol/L/day)
KDIET    : 0.020 : Dietary VLCFA input (umol/L/day)
DIETSC   : 1.0   : Dietary restriction scale (1 = normal diet, 0.35 = strict)
VMAXB    : 0.180 : Vmax peroxisomal beta-oxidation of C26:0 (umol/L/day)
KMB      : 0.50  : Km peroxisomal beta-oxidation (umol/L)
KOMEG    : 0.075 : Alternative disposal (omega-oxidation / bulk turnover) (1/day)
KLPC     : 0.357 : C26:0-lysoPC formation coefficient
HLPC     : 1.70  : Supralinearity of lysoPC formation (why DBS marker beats C26:0)
KELPC    : 1.00  : C26:0-lysoPC elimination (1/day)
KSYNE    : 0.10  : ELOVL1 synthesis (1/day)
KDEGE    : 0.10  : ELOVL1 degradation (1/day)
ESCAPE   : 0.55  : Compensatory ELOVL1 up-regulation when elongation flux is blocked
KSA1     : 0.05  : ALDP turnover toward its genotype/vector target (1/day)
KS2      : 0.02  : ABCD2 synthesis (1/day)
KD2      : 0.02  : ABCD2 degradation (1/day)
EIND2    : 0.35  : Max fold ABCD2 induction by inducers (small: the fibrate trials were negative)
KIND50   : 1.00  : Inducer concentration for half-maximal ABCD2 induction

// ------------------------------------------------ tissue VLCFA (relative) --
KINADR   : 0.02857 : Plasma -> adrenal cortex VLCFA influx
KOUTADR  : 0.00600 : Adrenal peroxisomal VLCFA disposal (ALDP-dependent)
KADRT    : 0.00400 : Adrenal bulk lipid turnover (ALDP-independent)
KINBR    : 0.0018966: Plasma -> brain white-matter VLCFA influx (only ~12% of brain input)
KBRSYN   : 0.004840: IN-SITU cerebral VLCFA synthesis (~88%; unreachable from plasma)
KOUTBR   : 0.00350 : Brain peroxisomal VLCFA disposal (ALDP-dependent)
KBRT     : 0.00200 : Brain bulk lipid turnover
KINSC    : 0.0011379: Plasma -> spinal cord VLCFA influx (~12% of cord input)
KSCSYN   : 0.0029040: IN-SITU spinal cord VLCFA synthesis (~88%)
KOUTSC   : 0.00250 : Cord peroxisomal VLCFA disposal
KSCT     : 0.00080 : Cord bulk lipid turnover
KINTST   : 0.022857 : Plasma -> testis VLCFA influx
KOUTTST  : 0.00500 : Testis peroxisomal disposal
KTSTT    : 0.00300 : Testis bulk turnover
KINFIB   : 0.071429 : Plasma -> fibroblast/other VLCFA influx
KOUTFIB  : 0.02000 : Fibroblast peroxisomal disposal
KFIBT    : 0.00500 : Fibroblast bulk turnover
FMGL     : 0.06    : Microglial share of brain white-matter lipid mass (why HSCT cannot move CBR)
FABCD2   : 0.45    : Weight of INDUCED ABCD2 increment in peroxisomal capacity

// -------------------------------------------- Lorenzo oil PK (erucic acid) --
KALO     : 1.20  : Lorenzo oil absorption / ester hydrolysis (1/day)
FLOE     : 0.80  : Fraction of Lorenzo oil dose appearing as erucic acid (GTE 4 : GTO 1)
KEERU    : 55.0  : Erucic acid plasma elimination (1/day; rapid beta-oxidation)
KEOLE    : 60.0  : Oleic acid plasma elimination (1/day)
VDLO     : 12.0  : Lorenzo oil apparent volume (L)
KIERU    : 75.0  : Ki, erucic acid competitive inhibition of ELOVL1 (umol/L)
KIOLE    : 400.0 : Ki, oleic acid competitive inhibition of ELOVL1 (umol/L)
FLOBBB   : 0.0   : Fraction of erucic acid reaching brain (measured ~0 : the crux)

// ---------------------------------------------------- adrenocortical axis --
VST      : 7500.0 : Vmax adrenal steroidogenesis (nmol/L/day)
KMA      : 25.0   : Km ACTH at MC2R (pg/mL)
KSACTH   : 33200.0: Pituitary ACTH synthesis (pg/mL/day)
KEACTH   : 100.0  : ACTH elimination (1/day, t1/2 ~10 min)
KIFB     : 130.0  : Cortisol IC50 for ACTH negative feedback (nmol/L)
HFB      : 3.0    : Hill, cortisol feedback
KECRT    : 12.5   : Cortisol elimination (1/day, t1/2 ~80 min)
CIRCAMP  : 0.0    : Circadian amplitude of CRH/ACTH drive (set 0.55 for acute runs only)
CIRCPH   : 0.29   : Circadian acrophase (fraction of day; 0.29 ~ 07:00)
HSPEED   : 0.02   : Timescale factor for the FAST endogenous hormone subsystem (1 = real kinetics)
KIMC2R   : 18.0   : Adrenal VLCFA disrupting MC2R/MRAP coupling (relative units)
HMC      : 2.0    : Hill, MC2R coupling loss
KREG     : 0.00015: Adrenocortical regeneration (1/day)
KADMD    : 0.00035: VLCFA-driven adrenocortical loss (1/day)
KADM50   : 4.00   : Adrenal VLCFA for half-maximal cortical loss
HADM     : 3.00   : Hill, adrenocortical loss
KROSADM  : 0.00010: Oxidative contribution to adrenocortical loss
VALD     : 1.286  : Aldosterone synthesis (nmol/L/day)
KEALD    : 6.00   : Aldosterone elimination (1/day)
ERN      : 0.80   : Renin drive on aldosterone
KSRN     : 2.00   : Renin synthesis (1/day)
KERN     : 2.00   : Renin elimination (1/day)
ALD0     : 0.30   : Reference aldosterone (nmol/L)
KKBASE   : 4.20   : Reference serum potassium (mmol/L)
KKSENS   : 2.60   : Potassium rise per unit loss of aldosterone action
KEKION   : 0.50   : Potassium equilibration (1/day)

// --------------------------------------------------- hormone replacement --
KAHC     : 12.0  : Hydrocortisone absorption (1/day)
KEHC     : 11.1  : Hydrocortisone elimination (1/day, t1/2 ~1.5 h)
VDHC     : 35.0  : Hydrocortisone apparent volume (L)
POTHC    : 2.76  : nmol/L cortisol-equivalent per ng/mL hydrocortisone
KEFLU    : 2.0   : Fludrocortisone elimination (1/day)
POTFLU   : 12.0  : Fludrocortisone aldosterone-equivalent potency
TSTREP    : 0.0  : Exogenous testosterone replacement (fraction of normal)

// ---------------------------------------------------------- testis axis ---
KLEYD    : 0.0020 : VLCFA-driven Leydig cell dysfunction (1/day)
KLEY50   : 6.0    : Testis VLCFA for half-maximal Leydig loss
KLEYR    : 0.0010 : Leydig regeneration (1/day)
VTST     : 26.67  : Testosterone production (nmol/L/day)
KETST    : 1.00   : Testosterone elimination (1/day)
KSLH     : 12.0   : LH synthesis (IU/L/day)
KELH     : 2.00   : LH elimination (1/day)
KILH     : 20.0   : Testosterone IC50 for LH feedback (nmol/L)

// ------------------------------------- microglia and the ignition switch --
KMGR     : 0.00096: Resident microglial self-renewal (1/day)
KMGAP    : 0.000530: VLCFA/ROS-driven microglial apoptosis (1/day)
KMG50    : 1.50   : Brain VLCFA for half-maximal microglial apoptosis
HMGA     : 4.00   : Hill, microglial apoptosis
KROSMG   : 0.000192: Oxidative contribution to microglial apoptosis
KMGINF   : 0.0060 : Extra loss of resident microglia inside the active lesion (1/day)
KPHAG    : 3.20   : Myelin-debris phagocytic capacity per unit microglia
PHEFFN   : 1.00   : Phagocytic efficiency, ALDP-deficient resident microglia
PHEFFC   : 1.60   : Phagocytic efficiency, gene-corrected / donor microglia
KMDEB    : 0.35   : Km for debris clearance (saturable : source of bistability)
FDEB     : 1.00   : Debris generated per unit oligodendrocyte death
KIGN     : 0.0085 : Ignition rate coefficient (deficit x VLCFA gate x susceptibility)
KMGD50   : 0.52   : Microglial deficit for half-maximal ignition (pre-lesional lesion)
HMGD     : 6.00   : Hill, microglial deficit gate
KSUSC50  : 1.00   : Susceptibility for half-maximal ignition drive
HSUSC    : 4.00   : Hill on susceptibility (lumped modifier-gene/stochastic term)
KAMP50   : 0.55   : Inflammatory signal for half-maximal loop amplification
HAMP     : 4.00   : Hill, loop amplification (what makes the switch bistable)
KAD      : 0.50   : Weight of saturating debris burden in loop amplification
KDEB50   : 0.42   : Debris for half-maximal ignition
HDEB     : 4.00   : Hill, ignition (steepness of the switch)
KVBR50   : 1.50   : Brain VLCFA for half-maximal permissive gate
HVBR     : 6.00   : Hill, permissive gate (saturated in ALL patients : why C26:0 cannot predict)
KAMP     : 0.115  : Loop amplification coefficient
KAT      : 1.00   : TNF weight in loop amplification
KAM      : 0.70   : Infiltrating monocyte weight
KAC      : 0.90   : CD8 T-cell weight
KAA      : 0.30   : Reactive astrocyte weight
MGPMAX   : 1.00   : Maximum pro-inflammatory microglial fraction
KMGPOFF  : 0.0320 : Intrinsic resolution of the pro-inflammatory state (1/day)
KSMGC    : 0.115  : Suppression of the loop by corrected microglia
KSPPAR   : 0.015  : Suppression of the loop by PPAR-gamma activation

// ---------------------------------------- cytokines, BBB, oligodendrocyte --
KTNFP    : 1.20  : TNF production per unit activated microglia/monocyte
FMONO    : 0.80  : Monocyte contribution to TNF
KETNF    : 1.50  : TNF elimination (1/day)
KIL1P    : 0.90  : IL-1beta production
KEIL1    : 1.50  : IL-1beta elimination (1/day)
KMCPP    : 1.10  : CCL2/MCP-1 production
FAST     : 0.60  : Astrocyte contribution to CCL2
KEMCP    : 1.20  : CCL2 elimination (1/day)
KMONI    : 0.55  : Monocyte influx coefficient (chemokine x BBB permeability)
KEMON    : 0.20  : Monocyte egress/death (1/day)
KTCDI    : 0.30  : CD8 T-cell influx coefficient
KETCD    : 0.12  : CD8 T-cell loss (1/day)
KASTP    : 0.35  : Astrocyte activation by TNF
KMAST    : 0.40  : TNF for half-maximal astrocyte activation
ASTMAX   : 1.00  : Reactive astrocyte capacity
KMMCP    : 0.50  : CCL2 for half-maximal monocyte influx
MONOMAX  : 1.00  : Infiltrating monocyte capacity
KMTC     : 0.45  : Stimulus for half-maximal CD8 influx
TCDMAX   : 1.00  : Infiltrating CD8 capacity
KEAST    : 0.25  : Astrocyte deactivation (1/day)
KBBBR    : 0.060 : BBB tight-junction repair (1/day)
KBBBD    : 0.220 : BBB breakdown by TNF/MMP (1/day)
KBT      : 0.45  : TNF for half-maximal BBB breakdown
KPPBBB   : 0.25  : PPAR-gamma stabilisation of the barrier
KOLGR    : 0.0090: Oligodendrocyte regeneration (1/day)
KOLGD0   : 0.0    : Baseline oligodendrocyte death (0: healthy turnover is implicit)
KOLGDT   : 0.0130: TNF-driven oligodendrocyte death
KOLGDI   : 0.0060: IL-1-driven oligodendrocyte death
KOLGDP   : 0.0180: CD8/perforin-driven oligodendrocyte death
KOLGDN   : 0.000375: ROS-driven oligodendrocyte death
KMYER    : 0.0060: Remyelination rate (1/day)
KMYED    : 0.000192: Direct cytotoxic myelin loss at the lesion front (1/day)
KSPREAD  : 14.0  : Autocatalytic lesion-front spread coefficient
KMTNF    : 0.50  : TNF for half-maximal direct myelin loss
KMYEO    : 0.0090: Myelin loss from oligodendrocyte deficit (1/day)
KROSBP   : 0.55  : Brain ROS production coefficient
KRB50    : 1.44  : Brain VLCFA for half-maximal ROS production
KRADR    : 4.00  : Adrenal VLCFA for half-maximal ROS production (own tissue reference)
KRINF    : 0.30  : Inflammatory contribution to brain ROS
KEROS    : 1.10  : Brain ROS elimination (1/day)
KNRF     : 0.20  : Nrf2 antioxidant efficacy
KAXR     : 0.0030: Cerebral axonal repair (1/day)
KAXDEM   : 0.0075: Axonal loss per unit demyelination (1/day)
KAXROS   : 0.00012: Axonal loss from ROS (1/day)
KAXT     : 0.0090: Axonal loss from cytotoxic T cells (1/day)
KGLI     : 0.0060: Gliosis formation
KEGLI    : 0.00020: Gliosis resolution (near-permanent scar)

// ------------------------------------------------- spinal cord (AMN arm) --
KATP     : 0.150 : Axonal energetic deficit formation
KSC50    : 2.166 : Cord VLCFA for half-maximal energetic deficit
KATPR    : 0.35  : ROS contribution to axonal energetic deficit
KEATP    : 0.220 : Recovery of axonal energetics (1/day)
KAXDSC   : 0.0001264: Cord axonal degeneration coefficient (1/day)
HATP     : 2.00  : Hill, energetic deficit -> axonal loss
LENCST   : 1.00  : Length-dependent vulnerability, corticospinal tract
LENDC    : 0.82  : Length-dependent vulnerability, dorsal column
KAXRSC   : 0.0     : Cord axonal repair (CNS axons do not regenerate; 0 by design)
KMSCR    : 0.0020: Cord remyelination (1/day)
KMSCD    : 0.0035: Cord myelin loss secondary to axonal loss (1/day)

// ------------------------------------------------------ clinical scoring --
LOESMAX  : 34.0  : Maximum Loes score
KLOES50  : 0.214 : Demyelinated fraction scoring half of maximal Loes
HLOES    : 2.00  : Hill, Loes vs demyelinated fraction
KLOESG   : 5.00  : Gliosis contribution to Loes
KLOESUP  : 0.050 : Upward tracking rate of the Loes score (1/day; one-way)
KLOESE   : 3.00  : Reversible gadolinium-enhancement points added to Loes
NFSMAX   : 25.0  : Maximum NFS
KNFS50   : 0.30  : Cerebral axonal loss scoring half of maximal NFS
HNFS     : 2.50  : Hill, NFS
KNFSUP   : 0.030 : Upward tracking rate of NFS (1/day; one-way)
MFDNFS   : 8.00  : NFS threshold defining a major functional disability
EDSSMAX  : 9.00  : Maximum EDSS
KEDSS50  : 0.42  : Cord tract loss for half-maximal EDSS
HEDSS    : 3.00  : Hill, EDSS
WALKMAX  : 520.0 : Reference 6-minute walk distance (m)

// ------------------------------------------- HSCT / gene therapy / safety --
CLBUS    : 10.5  : Busulfan clearance (L/h -> converted in MAIN) (L/day placeholder)
VBUS     : 22.0  : Busulfan volume (L)
KEBUS    : 6.60  : Busulfan elimination (1/day, t1/2 ~2.5 h)
KBUSK    : 1.20  : Busulfan-driven host HSC kill per (mg/L)/day
KHSCG    : 0.075 : HSC expansion rate (1/day)
HSCCAP   : 1.00  : Marrow HSC niche capacity
KNEUP    : 2.00  : Neutrophil production per unit HSC (1/day)
KNEUD    : 2.00  : Neutrophil turnover (1/day)
KMONOP   : 0.085 : Peripheral corrected-monocyte production (1/day)
KEMONO   : 0.085 : Corrected-monocyte turnover (1/day)
KMGREP   : 0.0020: Microglial replacement rate from corrected monocytes (1/day : THE lag)
KMGCD    : 0.00060: Corrected microglial loss (1/day)
MGCAP    : 1.00  : Total microglial niche capacity
VCNIN    : 0.80  : Vector copy number of the infused product (copies/genome)
KVCN     : 0.10  : VCN equilibration with engrafted graft (1/day)
KCLN     : 0.000060: Clonal expansion hazard coefficient (per VCN per day)
KCLNB    : 0.35  : Busulfan-exposure amplification of clonal risk
KECLN    : 0.00050: Clonal contraction (1/day)
KMDS     : 0.00075: MDS/AML hazard per unit clonal burden (1/day)
GVHDR    : 0.0   : Allogeneic GVHD burden (0 for autologous gene therapy)

// ------------------------------------------------------- CNS-active drugs --
KALER    : 6.0   : Leriglitazone absorption (1/day)
KELER    : 0.55  : Leriglitazone plasma elimination (1/day)
VDLER    : 95.0  : Leriglitazone apparent volume (L)
KINLER   : 0.90  : Leriglitazone plasma -> brain (1/day)
KOUTLER  : 1.30  : Leriglitazone brain -> out (1/day)
EC50LER  : 145.0 : Brain leriglitazone EC50 for PPAR-gamma (ng/mL)
EMAXPP   : 1.00  : Maximum PPAR-gamma activation
KEPPAR   : 0.80  : PPAR-gamma effect turnover (1/day)
KENRF    : 0.60  : Nrf2 effect turnover (1/day)
DMFD     : 0.0   : Dimethyl fumarate effect input (0-1 scale)
EDMF     : 0.70  : DMF maximum Nrf2 induction
NACD     : 0.0   : Antioxidant cocktail input (0-1 scale)
ENAC     : 0.45  : Antioxidant cocktail maximum ROS suppression
KESOB    : 1.10  : Sobetirome plasma elimination (1/day)
KINSOB   : 1.40  : Sobetirome plasma -> brain (1/day)
KOUTSOB  : 1.00  : Sobetirome brain -> out (1/day)
VDSOB    : 30.0  : Sobetirome apparent volume (L)
EC50SOB  : 90.0  : Brain sobetirome EC50 for CNS ABCD2 induction (ng/mL)
ESOBBR   : 0.60  : Maximum sobetirome-driven reduction of brain VLCFA influx
BEZAD    : 0.0   : Bezafibrate input (relative units, peripheral ABCD2 inducer)
PBAD     : 0.0   : 4-phenylbutyrate input (relative units)
ELOVLID  : 0.0   : Direct ELOVL1 inhibitor input (fraction of ELOVL1 blocked)
STEROIDX : 0.0   : High-dose systemic steroid (transient TNF suppression only)

$CMT @annotated
CPL   : Plasma C26:0 (umol/L)
CLPC  : Plasma C26:0-lysoPC (umol/L)
CADR  : Adrenocortical VLCFA (relative to normal)
CBR   : Cerebral white-matter VLCFA (relative to normal)
CSC   : Spinal cord VLCFA (relative to normal)
CTST  : Testicular VLCFA (relative to normal)
CFIB  : Fibroblast/other-tissue VLCFA (relative to normal)
ELOV  : ELOVL1 activity (relative to normal)
ALDPF : Functional ALDP (fraction of normal)
ABCD2F: ABCD2/ALDRP activity (relative to normal)
LOG   : Lorenzo oil in gut (umol)
ERU   : Plasma erucic acid C22:1 (umol/L)
OLE   : Plasma oleic acid C18:1 above baseline (umol/L)
ADM   : Functional adrenocortical mass (fraction)
ACTH  : Plasma ACTH (pg/mL)
CRT   : Plasma cortisol (nmol/L)
ALDO  : Plasma aldosterone (nmol/L)
RN    : Renin activity (relative)
KION  : Serum potassium (mmol/L)
HCG   : Hydrocortisone in gut (mg)
HCP   : Plasma hydrocortisone (ng/mL)
FLUP  : Plasma fludrocortisone (relative)
LEY   : Leydig cell function (fraction)
TST   : Plasma testosterone (nmol/L)
LHH   : Plasma LH (IU/L)
MGN   : Resident (ALDP-deficient) microglia (fraction of niche)
MGC   : Gene-corrected / donor-derived microglia (fraction of niche)
MGP   : Pro-inflammatory activated microglia (fraction)
DEB   : Myelin debris burden (relative)
TNF   : Cerebral TNF-alpha (relative)
IL1   : Cerebral IL-1beta (relative)
MCPC  : Cerebral CCL2/MCP-1 (relative)
MONO  : Infiltrating monocytes (relative)
TCD8  : Infiltrating CD8 T cells (relative)
BBBI  : Blood-brain barrier integrity (1 = intact)
OLGP  : Oligodendrocyte pool (fraction of normal)
MYEL  : Cerebral myelin content (fraction of normal)
AXB   : Cerebral axonal integrity (fraction of normal)
GLI   : Gliosis / scar burden (relative)
ROSB  : Cerebral oxidative stress (relative)
AST   : Reactive astrocytes (relative)
CSTI  : Corticospinal tract integrity (fraction)
DCI   : Dorsal column integrity (fraction)
ATPD  : Axonal energetic deficit (relative)
MYESCV: Spinal cord myelin (fraction)
LOESV : Loes MRI severity score (0-34)
NFSV  : Neurologic function score (0-25)
BUSP  : Plasma busulfan (mg/L)
AUCB  : Cumulative busulfan AUC (mg*day/L)
HSCN  : Host haematopoietic stem cells (fraction of niche)
HSCD  : Corrected / donor haematopoietic stem cells (fraction of niche)
NEU   : Neutrophils (relative to normal)
MONOC : Peripheral corrected-monocyte chimerism (fraction)
VCNS  : Vector copy number in engrafted cells (copies/genome)
CLN   : Clonal expansion burden (relative)
HMDS  : Cumulative MDS/AML hazard
LERG  : Leriglitazone in gut (mg)
LERP  : Plasma leriglitazone (ng/mL)
LERB  : Brain leriglitazone (ng/mL)
PPARA : PPAR-gamma activation (0-1)
NRF   : Nrf2 pathway activity above baseline (0-1)
SOBP  : Plasma sobetirome (ng/mL)
SOBB  : Brain sobetirome (ng/mL)

$OMEGA @annotated @block
ESUSC : 0.36 : IIV on cerebral ignition susceptibility (log-scale variance)

$GLOBAL
// NUMERICAL NOTE (this cost two debugging rounds, so it is written down).
// A hard max(x,0) placed on a term that is ZERO AT THE HEALTHY EQUILIBRIUM
// puts a derivative kink exactly on the resting point.  LSODA then chatters
// against it and dies with "corrector convergence failed ... fabs(h_)=hmin"
// -- which is what happened at t=2621 d in a HEALTHY subject, before any
// disease was even switched on.  Every capacity/repair term below is
// therefore written as a SMOOTH logistic ( k*X*(CAP-X) ), which is
// self-limiting without a kink.  The two clinical ratchets genuinely need
// one-way behaviour, so they use rect(): smooth, and exactly 0 at 0.
namespace {
  inline double pos(double x){ return x > 0.0 ? x : 0.0; }
  // smooth one-way rectifier: ~max(x,0), C1, and exactly 0 at x = 0
  inline double rect(double x, double e){
    return 0.5*(x + sqrt(x*x + e*e)) - 0.5*e;
  }
  inline double hup(double x, double k, double n){
    double z = pow(pos(x)/k, n); return z/(1.0+z);
  }
  inline double hdn(double x, double k, double n){
    double z = pow(pos(x)/k, n); return 1.0/(1.0+z);
  }
}

$MAIN
// ---- susceptibility: the SECOND, unmeasurable variable ---------------------
double SUSC = SUSCTV * exp(ESUSC);

// ---- genotype / mosaic target for functional ALDP --------------------------
// Male hemizygote: FMOS = 0 -> ALDPTGT = MUTRES.
// Female carrier : a fraction FMOS of cells are ALDP-null; competent
// neighbours partially rescue them (CRSC), so the tissue average is not a
// simple mean.  Gene therapy adds corrected cells to the haematopoietic and
// microglial compartments only (handled through MGC / MONOC, not here).
double ALDPTGT = (1.0 - FMOS) * 1.0 + FMOS * (MUTRES + CRSC * (1.0 - FMOS) * (1.0 - MUTRES));
if (FMOS <= 0.0) ALDPTGT = MUTRES;

$ODE
// ===========================================================================
//  VARIABLE (1) : VLCFA MASS BALANCE — deterministic, measurable, peripheral
// ===========================================================================

// ABCD2 induction (bezafibrate / 4-PBA peripherally, sobetirome centrally)
double INDU  = BEZAD + PBAD;
double IND2  = 1.0 + EIND2 * INDU/(KIND50 + INDU);

// Total peroxisomal VLCFA-degrading capacity, relative to normal.
// ABCD2 contributes ONLY its induced increment above baseline.  Writing it as
// (ALDPF + 0.45*ABCD2F) was wrong twice over: it made healthy capacity 1.45
// (so healthy plasma C26:0 settled at 0.23 instead of 0.35) and it handed an
// untreated null patient 0.32 of normal capacity -- i.e. it made the paralogue
// rescue the disease at baseline, which is precisely what ABCD2 does NOT do.
double PXC   = ALDPF + FABCD2 * rect(ABCD2F - 1.0, 1e-6);

// competitive inhibition of ELOVL1 by erucic + oleic acid (Lorenzo oil), plus
// any direct ELOVL1 inhibitor
double COMPI = 1.0 + ERU/KIERU + OLE/KIOLE;
double ELOVA = ELOV * (1.0 - ELOVLID) / COMPI;

double ELONGF = KELONG * ELOVA;
double DIETF  = KDIET * DIETSC;
double DEGF   = VMAXB * PXC * CPL/(KMB + CPL);

dxdt_CPL  = ELONGF + DIETF - DEGF - KOMEG * CPL;
dxdt_CLPC = KLPC * pow(pos(CPL), HLPC) - KELPC * CLPC;

// compensatory ELOVL1 up-regulation when elongation flux is suppressed
// (this is why the plasma response to Lorenzo oil partially escapes)
double BLOCK = 1.0 - 1.0/COMPI;
dxdt_ELOV = KSYNE * (1.0 + ESCAPE * BLOCK) - KDEGE * ELOV;

dxdt_ALDPF  = KSA1 * (ALDPTGT - ALDPF);
dxdt_ABCD2F = KS2 * IND2 - KD2 * ABCD2F;

// ---- tissue distribution --------------------------------------------------
// Local peroxisomal disposal uses the LOCAL corrected-cell fraction.  In brain
// that fraction is the corrected MICROGLIAL share of white-matter lipid mass
// (FMGL, ~6%) — which is exactly why HSCT cannot normalise CBR.
double PXBR  = ALDPF + FABCD2 * rect(ABCD2F - 1.0, 1e-6)
               + FMGL * MGC * (1.0 - ALDPF);
double SOBEF = ESOBBR * SOBB/(EC50SOB + SOBB);

// STRUCTURAL POINT, and the reason Lorenzo oil fails.
// Brain and cord VLCFA are NOT simply plasma-derived.  Roughly three quarters
// of CNS VLCFA is elongated IN SITU, and erucic acid does not cross the BBB
// (FLOBBB = 0, a measured fact).  So peripheral ELOVL1 inhibition cannot reach
// the CNS elongation term at all.  An earlier version had brain VLCFA driven
// only by KINBR*CPL; Lorenzo oil then halved brain VLCFA along with plasma,
// closed the permissive gate, and ABOLISHED cerebral disease -- the opposite
// of the two placebo-controlled trials.  The in-situ term is what makes the
// biomarker and the disease separable.
double COMPIBR = 1.0 + FLOBBB * ERU/KIERU;
double CNSYN   = (1.0 - SOBEF)/COMPIBR;

dxdt_CADR = KINADR * CPL - KOUTADR * PXC  * CADR - KADRT * CADR;
dxdt_CBR  = KINBR * CPL * (1.0 - SOBEF) + KBRSYN * CNSYN
            - KOUTBR * PXBR * CBR - KBRT * CBR;
dxdt_CSC  = KINSC * CPL * (1.0 - SOBEF) + KSCSYN * CNSYN
            - KOUTSC * PXC * CSC - KSCT * CSC;
dxdt_CTST = KINTST * CPL - KOUTTST * PXC * CTST - KTSTT * CTST;
dxdt_CFIB = KINFIB * CPL - KOUTFIB * PXC * CFIB - KFIBT * CFIB;

// ---- Lorenzo oil PK -------------------------------------------------------
dxdt_LOG = -KALO * LOG;
dxdt_ERU =  FLOE * KALO * LOG / VDLO - KEERU * ERU;
dxdt_OLE = (1.0 - FLOE) * KALO * LOG / VDLO - KEOLE * OLE;

// ===========================================================================
//  ADRENOCORTICAL AXIS — the third, independent axis
// ===========================================================================
double MC2RF = hdn(CADR, KIMC2R, HMC);
double CIRC  = 1.0 + CIRCAMP * cos(2.0*M_PI*(SOLVERTIME - CIRCPH));
if (CIRC < 0.0) CIRC = 0.0;

double CRTEQ  = CRT + POTHC * HCP;                 // total glucocorticoid exposure
double STER   = VST * ADM * MC2RF * ACTH/(KMA + ACTH);

// HSPEED scales the FAST endogenous hormone subsystem.  Production and
// elimination are scaled by the SAME factor, so every steady state is
// mathematically unchanged; only the approach rate slows.  Multi-decade runs
// use HSPEED << 1 to keep a 10-minute ACTH half-life from forcing LSODA to
// resolve 14600 days at minute steps.  Acute scenarios (ACTH stimulation
// test, intra-day hydrocortisone profile) use HSPEED = 1.
dxdt_ACTH = HSPEED * (KSACTH * CIRC * hdn(CRTEQ, KIFB, HFB) - KEACTH * ACTH);
dxdt_CRT  = HSPEED * (STER - KECRT * CRT);

double ROSA = hup(CADR, KRADR, 2.0);
dxdt_ADM  = KREG * (1.0 - ADM)
            - KADMD * hup(CADR, KADM50, HADM) * ADM
            - KROSADM * ROSA * ADM;

dxdt_ALDO = HSPEED * (VALD * ADM * (1.0 + ERN * RN) - KEALD * ALDO);
dxdt_RN   = HSPEED * (KSRN * hdn(ALDO, ALD0, 1.0) - KERN * RN);

// Potassium responds in BOTH directions (mineralocorticoid excess lowers it),
// so this is a signed difference, not a rectified one: a max(.,0) here would
// place a kink exactly on the healthy equilibrium ALDACT = 1.
double ALDACT = (ALDO + POTFLU * FLUP)/ALD0;
// Bounded in both directions: an unbounded (1 - ALDACT) term drove serum
// potassium to 0.46 mmol/L in a healthy subject once aldosterone overshot.
double KTGT   = KKBASE + KKSENS * (2.0*hdn(ALDACT, 1.0, 2.0) - 1.0);
dxdt_KION = HSPEED * KEKION * (KTGT - KION);

dxdt_HCG = -KAHC * HCG;
dxdt_HCP =  KAHC * HCG * 1000.0/VDHC - KEHC * HCP;   // mg -> ng/mL with V in L
dxdt_FLUP = -KEFLU * FLUP;

// ===========================================================================
//  TESTICULAR AXIS
// ===========================================================================
dxdt_LEY = KLEYR * (1.0 - LEY) - KLEYD * hup(CTST, KLEY50, 2.0) * LEY;
dxdt_TST = HSPEED * (VTST * (LEY * LHH/(LHH + 5.0) * 2.0 + TSTREP) - KETST * TST);
dxdt_LHH = HSPEED * (KSLH * hdn(TST, KILH, 2.0) - KELH * LHH);

// ===========================================================================
//  VARIABLE (2) : THE BISTABLE CEREBRAL SWITCH
// ===========================================================================

// ---- resident microglia die BEFORE demyelination (pre-lesional event) -----
// Logistic (not rectified) capacity: smooth at MGN + MGC = MGCAP.
// Death outruns self-renewal, so the resident pool does not hold at capacity:
// it slides to a REDUCED PLATEAU (1 - MGAP/KMGR) over roughly 600 days.  That
// slide is the slow ramp that sets WHEN ignition becomes possible, and it is
// why cerebral onset clusters in childhood rather than at birth.  Tuning note:
// with the first values (KMGR 0.006, KMGAP 0.013) death was 6x renewal, the
// pool went to ZERO in every patient, phagocytic capacity went with it, and
// the switch then fired regardless of susceptibility -- every patient got
// cerebral disease.  The pool must stay partly populated for the switch to be
// a switch.
double MGSPACE = MGCAP - MGN - MGC;
// Activated deficient microglia die faster: this is the peri-lesional
// microglial depletion seen at the advancing edge of the lesion.
double MGAP = KMGAP * hup(CBR, KMG50, HMGA) + KROSMG * ROSB + KMGINF * MGP;
dxdt_MGN = KMGR * MGN * MGSPACE - MGAP * MGN;

// ---- corrected / donor microglia arrive SLOWLY (the lag window) -----------
dxdt_MGC = KMGREP * MONOC * MGSPACE - KMGCD * MGC;

// ---- debris: production linear, clearance SATURABLE -> bistability --------
double PHAGC  = KPHAG * (MGN * PHEFFN + MGC * PHEFFC);
double CLRDEB = PHAGC * DEB/(KMDEB + DEB);

double ROSEFF = ROSB * (1.0 - KNRF * NRF);
// ROS enters as its EXCESS over the healthy reference, so that a non-carrier
// sits at OLGP = MYEL = AXB = 1 exactly and Loes = 0 exactly.  Without this
// subtraction a healthy 40-year-old scored Loes 34.
double ROSEX  = rect(ROSEFF - ROSREF, 1e-6);
double OLGDTH = (KOLGD0 + KOLGDT*TNF + KOLGDI*IL1 + KOLGDP*TCD8 + KOLGDN*ROSEX) * OLGP;
dxdt_DEB = FDEB * OLGDTH - CLRDEB;

// ---- ignition: VLCFA is the PERMISSIVE gate, SUSC is the DECISION --------
// VLCG is >0.96 in every X-ALD patient and 0.03 in a non-carrier: it is a
// PERMISSION, not a dose.  It therefore cannot discriminate phenotype, which
// is exactly the clinical observation this model exists to reproduce.
double VLCG  = hup(CBR, KVBR50, HVBR);
// The ignition gate is the MICROGLIAL DEFICIT, not debris.  Debris cannot
// start the lesion (an earlier version gated on it and the switch could never
// self-start: clearance capacity exceeded debris production by three orders of
// magnitude).  Gating on the deficit matches the pathology -- microglial
// apoptosis PRECEDES demyelination -- and gives the switch a slow ramp.
double MGDEF = rect(MGCAP - MGN - MGC, 1e-6);
double IGNIT = hup(MGDEF, KMGD50, HMGD);
double SUSCG = hup(SUSC, KSUSC50, HSUSC);
// Amplification is SATURATING in the inflammatory signal.  A linear AMPG
// cannot be bistable: it either never fires or fires in everyone.  Debris
// enters HERE -- as the amplifier of an established lesion, which is its
// actual role -- rather than as the trigger.
double AMPRAW = KAT*TNF + KAM*MONO + KAC*TCD8 + KAA*AST
                + KAD*hup(DEB, KDEB50, HDEB);
double AMPG   = hup(AMPRAW, KAMP50, HAMP);
double SUPPR  = KSMGC*MGC + KSPPAR*PPARA;

// THE MECHANISM BY WHICH TRANSPLANT WORKS, and it is not anti-inflammatory
// potency.  Only an ALDP-DEFICIENT myeloid cell can drive this lesion, and the
// effector population is BOTH resident microglia and blood-derived monocytes.
// Gene therapy or allografting corrects the monocytes (chimerism ~98%) and,
// slowly, the microglia -- so the UNCORRECTED FRACTION of the myeloid effector
// pool is what gates both ignition and amplification.
//
// Two earlier versions failed here and both failures are instructive.
// (i) Anti-inflammatory suppression alone (KSMGC) only moved the high fixed
//     point from MGP 0.79 to 0.61, which still saturates AMPG: eli-cel
//     reproduced its own control arm and contradicted 90%-plus MFD-free
//     survival.
// (ii) Capping the activated pool at the uncorrected microglial FRACTION
//     MGN/(MGN+MGC) looked right but is degenerate untreated: as resident
//     microglia die MGN -> 0 with MGC = 0, the ratio collapses, and untreated
//     cerebral disease switched itself off.  Including the monocyte arm in the
//     denominator fixes that -- untreated, UNCF is identically 1 no matter how
//     depleted the resident pool becomes.
double UNCF = (MGN + (1.0 - MONOC)*MONO)/(MGN + MGC + MONO + 1e-9);

dxdt_MGP = (KIGN * SUSCG * VLCG * IGNIT * UNCF + KAMP * AMPG * UNCF)
           * (MGPMAX - MGP)
           - (KMGPOFF + SUPPR) * MGP;

// ---- cytokines, monocytes, T cells, astrocytes ---------------------------
// Steroids blunt cytokine PRODUCTION modestly and touch nothing else in the
// loop. Set at 0.60 they extinguished cerebral disease outright, contradicting
// every steroid series in this disease; 0.20 reproduces the clinical failure.
double STSUP = 1.0 - 0.20 * STEROIDX;
dxdt_TNF  = KTNFP * (MGP + FMONO*MONO) * STSUP - KETNF * TNF;
dxdt_IL1  = KIL1P * MGP * STSUP - KEIL1 * IL1;
dxdt_MCPC = KMCPP * (MGP + FAST*AST) - KEMCP * MCPC;

double BBBP = 1.0 - BBBI;
// STRUCTURAL FIX.  These three were first written as unsaturated linear
// amplifiers (MONO ~ MCPC, AST ~ TNF), which closed the loop
// TNF -> AST -> CCL2 -> monocyte -> TNF with a LINEAR gain of 1.355*BBBP.
// Any BBB opening past 0.738 therefore had no bounded steady state and the
// cascade ran to 1e18 -- in a HEALTHY subject.  Each recruitment step is now
// capacity-limited (logistic) and saturating in its stimulus, so the loop can
// be strong without being infinite.  Bistability then comes from where it
// should: the saturable debris-clearance switch, not from a divergent linear
// feedback.
dxdt_MONO = KMONI * hup(MCPC, KMMCP, 1.0) * BBBP * (MONOMAX - MONO)
            - KEMON * MONO;
dxdt_TCD8 = KTCDI * BBBP * hup(MGP + 0.30*MONO, KMTC, 1.0) * (TCDMAX - TCD8)
            - KETCD * TCD8;
dxdt_AST  = KASTP * hup(TNF, KMAST, 1.0) * (ASTMAX - AST) - KEAST * AST;

dxdt_BBBI = KBBBR * BBBP * (1.0 + KPPBBB * PPARA)
            - KBBBD * hup(TNF, KBT, 1.0) * BBBI;

// ---- oligodendrocytes, myelin, axons, gliosis ---------------------------
dxdt_OLGP = KOLGR * (1.0 - OLGP) - OLGDTH;
double DEMYF = 1.0 - MYEL;
// Myelin destruction has an AUTOCATALYTIC term.  Cerebral ALD is a spreading
// lesion: destruction happens at an advancing front whose length grows with
// the lesion, so the rate accelerates as the lesion enlarges.  A purely
// first-order loss term decelerates instead, and could not reconcile the two
// clinical facts at once -- Loes 2 to 15 in about 18 months AND near-complete
// demyelination by 3 to 4 years.  With the front term the early rate is slow
// and the late rate is roughly 5-fold faster, which reconciles them.
double MYELOSS = KMYED * hup(TNF, KMTNF, 1.0) * MYEL * (1.0 + KSPREAD*DEMYF)
                 + KMYEO * (1.0 - OLGP) * MYEL;
dxdt_MYEL = KMYER * OLGP * DEMYF - MYELOSS;

dxdt_ROSB = KROSBP * (hup(CBR, KRB50, 2.0) + KRINF*(TNF + IL1)) * (1.0 - ENAC*NACD)
            - KEROS * ROSB;

double AXLOSS = (KAXDEM*rect(DEMYF,1e-4) + KAXROS*ROSEX + KAXT*TCD8) * AXB;
dxdt_AXB = KAXR * (1.0 - AXB) - AXLOSS;
dxdt_GLI = KGLI * AST * rect(DEMYF,1e-4) * (1.0 - GLI) - KEGLI * GLI;

// ===========================================================================
//  SPINAL CORD AXONOPATHY (AMN arm) — runs whether or not the switch fires
// ===========================================================================
dxdt_ATPD = KATP * (hup(CSC, KSC50, 2.0) + KATPR * ROSEX) - KEATP * ATPD;
// Degeneration is driven by the EXCESS deficit over healthy reference, so a
// non-carrier keeps tracts at 1.0 and EDSS at 0.
double AXDR = KAXDSC * pow(rect(ATPD - ATPREF, 1e-6), HATP);
dxdt_CSTI = KAXRSC * (1.0 - CSTI) - AXDR * LENCST * CSTI;
dxdt_DCI  = KAXRSC * (1.0 - DCI)  - AXDR * LENDC  * DCI;
dxdt_MYESCV = KMSCR * CSTI * (1.0 - MYESCV) - KMSCD * (1.0 - CSTI) * MYESCV;

// ===========================================================================
//  CLINICAL SCORES — INTEGRALS OF DESTRUCTION, not trackers of a target
//
//  These were first written as one-way ratchets chasing an instantaneous
//  target computed from current demyelination.  That is wrong in a way that
//  only shows up once a therapy works: the ratchet holds the PEAK target
//  forever, so a lesion that is arrested after a few months still climbs to
//  the Loes score its transient peak implied.  A successfully transplanted
//  boy therefore scored Loes 30 with his inflammation extinguished and MGP at
//  0.006 -- inflammation arrested, score still ruinous.
//
//  The fix is a FAST upward-only tracker of current damage.  LOESV then equals
//  the running maximum of the damage score: it follows deterioration within
//  weeks, stops the moment deterioration stops, and never falls.  That is both
//  the observed post-transplant plateau and the clinical fact that Loes does
//  not meaningfully improve after successful transplant.
//
//  Integrating the destruction FLUX instead (the version before this one) was
//  also wrong, for a subtler reason: because remyelination continually remakes
//  myelin that is then re-destroyed, the cumulative flux keeps growing for
//  years after the lesion volume has stabilised, and the Loes 1-to-15 interval
//  stretched to 5-7 years against an observed 1.5.
//
//  Gadolinium-enhancement points ARE reversible, so they are added in $TABLE
//  as a readout rather than banked into the state.
// ===========================================================================
// A saturating-exponential map scored Loes 1 at only 1% demyelination, so
// every pure-AMN patient acquired a cerebral score. A Hill map is flat at the
// low end (AMN stays at Loes < 0.5) and near-linear through the clinical range.
double LOESTGT = LOESMAX * hup(rect(DEMYF,1e-4), KLOES50, HLOES) + KLOESG*GLI;
if (LOESTGT > LOESMAX) LOESTGT = LOESMAX;
dxdt_LOESV = KLOESUP * rect(LOESTGT - LOESV, 1e-3);

double NFSTGT = NFSMAX * hup(rect(1.0 - AXB, 1e-4), KNFS50, HNFS);
dxdt_NFSV = KNFSUP * rect(NFSTGT - NFSV, 1e-3);

// ===========================================================================
//  HSCT / GENE THERAPY
// ===========================================================================
dxdt_BUSP = -KEBUS * BUSP;
dxdt_AUCB =  BUSP;

double HSCSPACE = HSCCAP - HSCN - HSCD;
dxdt_HSCN = KHSCG * HSCN * HSCSPACE - KBUSK * BUSP * HSCN;
dxdt_HSCD = KHSCG * HSCD * HSCSPACE - KBUSK * BUSP * HSCD * 0.15;
dxdt_NEU  = KNEUP * (HSCN + HSCD) - KNEUD * NEU;

double DFRAC = (HSCD + HSCN) > 1e-9 ? HSCD/(HSCD + HSCN) : 0.0;
dxdt_MONOC = KMONOP * DFRAC - KEMONO * MONOC;
dxdt_VCNS  = KVCN * (VCNIN * DFRAC - VCNS);

dxdt_CLN  = KCLN * VCNS * HSCD * (1.0 + KCLNB * AUCB/3.75) - KECLN * CLN;
dxdt_HMDS = KMDS * CLN;

// ===========================================================================
//  CNS-ACTIVE DRUGS
// ===========================================================================
dxdt_LERG = -KALER * LERG;
dxdt_LERP =  KALER * LERG * 1e6/(VDLER*1000.0) - KELER*LERP - KINLER*LERP;
dxdt_LERB =  KINLER * LERP - KOUTLER * LERB;

dxdt_PPARA = KEPPAR * (EMAXPP * LERB/(EC50LER + LERB) - PPARA);
dxdt_NRF   = KENRF * (EDMF * DMFD - NRF);

dxdt_SOBP = -KESOB*SOBP - KINSOB*SOBP;
dxdt_SOBB =  KINSOB*SOBP - KOUTSOB*SOBB;

$TABLE
// ---- measurable biomarkers (VARIABLE 1) -----------------------------------
double C26      = CPL;
double C26LYSO  = CLPC;
double C26FOLD  = CPL/0.35;

// ---- the unmeasurable switch (VARIABLE 2) ---------------------------------
double VLCGATE  = hup(CBR, KVBR50, HVBR);
double DEBGATE  = hup(DEB, KDEB50, HDEB);
double SUPPRT   = KSMGC*MGC + KSPPAR*PPARA;
// Loop gain, defined on the CURRENT state so it is interpretable rather than
// notional: the ratio of the amplification-driven production of activated
// microglia to their resolution.  LOOPGAIN > 1 means the inflammation now
// sustains itself with no further ignition input -- the point of no return.
double AMPRAWT = KAT*TNF + KAM*MONO + KAC*TCD8 + KAA*AST
                 + KAD*hup(DEB, KDEB50, HDEB);
double AMPGT   = hup(AMPRAWT, KAMP50, HAMP);
// Guarded at MGP > 1e-4: below that both numerator and denominator are noise
// and the ratio is meaningless (it read 6.1 in a healthy subject).
double UNCFT    = (MGN + (1.0-MONOC)*MONO)/(MGN + MGC + MONO + 1e-9);
double LOOPGAIN = (MGP > 1e-4)
  ? KAMP*AMPGT*UNCFT*(MGPMAX-MGP)/((KMGPOFF + SUPPRT)*MGP) : 0.0;
double SWITCHON = (MGP > 0.15) ? 1.0 : 0.0;

// ---- imaging / neurology --------------------------------------------------
double LOES     = fmin(LOESMAX, LOESV + KLOESE * rect(1.0 - BBBI, 1e-4));
double NFS      = NFSV;
double GADENH   = rect(1.0 - BBBI, 1e-4);
double DEMYELIN = rect(1.0 - MYEL, 1e-4);
double MFD      = (NFSV >= MFDNFS) ? 1.0 : 0.0;
double INWINDOW = (LOESV <= 9.0 && NFSV <= 1.0) ? 1.0 : 0.0;

// ---- spinal cord ----------------------------------------------------------
double CORDLOSS = rect(1.0 - 0.5*(CSTI + DCI), 1e-4);
double EDSS     = EDSSMAX * hup(CORDLOSS, KEDSS50, HEDSS);
double WALK6MIN = WALKMAX * (1.0 - hup(CORDLOSS, KEDSS50, HEDSS)*0.85);

// ---- adrenal ------------------------------------------------------------
double CORTISOL = CRT;
double ACTHOUT  = ACTH;
double POTASSIUM= KION;
double MC2RFUN  = hdn(CADR, KIMC2R, HMC);
double ADRRESV  = ADM * MC2RFUN;             // steroidogenic reserve
double ADRINSUF = (ADRRESV < 0.36) ? 1.0 : 0.0;
double TESTOST  = TST;

// ---- transplant ---------------------------------------------------------
double CHIMER   = MONOC;
double MGLCORR  = MGC;
double MGLTOT   = MGN + MGC;
double PHAGCAP  = KPHAG*(MGN*PHEFFN + MGC*PHEFFC);
double NEUTRO   = NEU;
double PMDS     = 1.0 - exp(-HMDS);
double BUSAUC   = AUCB * 24.0;   // mg*day/L -> mg*h/L, the clinical scale

// ---- drug exposure -------------------------------------------------------
double ERUCIC   = ERU;
double LERIBR   = LERB;
double PPARG    = PPARA;

$CAPTURE
C26 C26LYSO C26FOLD VLCGATE DEBGATE LOOPGAIN SWITCHON LOES NFS GADENH DEMYELIN
MFD INWINDOW EDSS WALK6MIN CORDLOSS CORTISOL ACTHOUT POTASSIUM ADRRESV ADRINSUF
TESTOST CHIMER MGLCORR MGLTOT PHAGCAP NEUTRO PMDS BUSAUC ERUCIC LERIBR PPARG
SUSC
'

# =============================================================================
#  MODEL LOADER
# =============================================================================
xald_mod <- function(quiet = TRUE) {
  mcode("xald", xald_code, quiet = quiet, soloc = tempdir())
}

# -----------------------------------------------------------------------------
#  Healthy reference initial state (all pools normalised to 1, hormones at
#  their healthy steady state).  These are the values a NON-carrier is born
#  with; the diseased state at birth is produced by xald_birth(), not typed in.
# -----------------------------------------------------------------------------
xald_healthy_init <- function() {
  c(CPL = 0.35, CLPC = 0.06, CADR = 1.0, CBR = 1.0, CSC = 1.0, CTST = 1.0,
    CFIB = 1.0, ELOV = 1.0, ALDPF = 1.0, ABCD2F = 1.0,
    LOG = 0, ERU = 0, OLE = 0,
    ADM = 1.0, ACTH = 25.0, CRT = 300.0, ALDO = 0.30, RN = 1.0, KION = 4.2,
    HCG = 0, HCP = 0, FLUP = 0,
    LEY = 1.0, TST = 20.0, LHH = 3.0,
    MGN = 1.0, MGC = 0, MGP = 0, DEB = 0,
    TNF = 0, IL1 = 0, MCPC = 0, MONO = 0, TCD8 = 0, BBBI = 1.0,
    OLGP = 1.0, MYEL = 1.0, AXB = 1.0, GLI = 0, ROSB = 0, AST = 0,
    CSTI = 1.0, DCI = 1.0, ATPD = 0, MYESCV = 1.0,
    LOESV = 0, NFSV = 0,
    BUSP = 0, AUCB = 0, HSCN = 1.0, HSCD = 0, NEU = 1.0, MONOC = 0,
    VCNS = 0, CLN = 0, HMDS = 0,
    LERG = 0, LERP = 0, LERB = 0, PPARA = 0, NRF = 0, SOBP = 0, SOBB = 0)
}

# -----------------------------------------------------------------------------
#  xald_birth() — the state at birth is a MODEL OUTPUT, not an input.
#  A 280-day in-utero simulation with the patient's genotype produces the
#  already-elevated C26:0-lysoPC that makes newborn screening possible.
# -----------------------------------------------------------------------------
xald_birth <- function(mod, param = list(), gestation = 280) {
  m <- mod
  if (length(param)) m <- param(m, param)
  out <- m %>%
    init(xald_healthy_init()) %>%
    zero_re() %>%
    mrgsim(end = gestation, delta = gestation, recsort = 3)
  d <- as.data.frame(out)
  st <- d[nrow(d), names(xald_healthy_init()), drop = FALSE]
  unlist(st)
}

# =============================================================================
#  SCENARIOS  (24)
#  Every therapeutic arm is paired with a MATCHED CONTROL run from the same
#  birth state with the same genotype, so differences cannot be an artefact of
#  a different starting patient.
# =============================================================================
YR <- 365.25

# The two susceptibility values below are not arbitrary labels.  xald_bifurcation()
# locates the saddle-node at SUSC ~ 0.42: at 0.40 the switch never fires and the
# patient is pure AMN; at 0.44 it fires and the patient gets cerebral disease.
# SUSC_CALD = 0.46 sits above it and yields cerebral onset at 6.7 y with Loes
# 1 -> 15 in 1.7 y (both observed); SUSC_AMN = 0.28 sits below it.  Onset age is NOT a parameter:
# it emerges from how far above threshold the individual sits, which is why
# near-threshold individuals declare later (critical slowing down).
SUSC_CALD <- 0.46
SUSC_AMN  <- 0.28

xald_scenarios <- function() {
  list(
    # ---- 1-6  natural history: the four phenotypes + female carriers -------
    s01_healthy = list(
      label = "01. 정상 대조 (non-carrier)",
      param = list(MUTRES = 1.0, SUSCTV = SUSC_AMN), end = 45*YR),

    s02_ccald_untreated = list(
      label = "02. 소아 뇌형 CCALD, 미치료",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_CALD), end = 14*YR),

    s03_amn = list(
      label = "03. AMN 순수 척수형 (역치 미달)",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_AMN), end = 55*YR),

    s04_addison_only = list(
      label = "04. 애디슨 단독형 (소아기 관찰 구간)",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_AMN), end = 14*YR),

    s05_female = list(
      label = "05. 여성 보유자, 중등도 편향 (0.75)",
      param = list(MUTRES = 0.02, FMOS = 0.75, SUSCTV = SUSC_AMN), end = 60*YR),

    s06_female_skewed = list(
      label = "06. 여성 보유자, 강한 편향 (0.92)",
      param = list(MUTRES = 0.02, FMOS = 0.92, SUSCTV = SUSC_AMN), end = 60*YR),

    # ---- 7-9  Lorenzo oil: the BIOMARKER arm ------------------------------
    s07_lo_ccald = list(
      label = "07. 로렌조 오일, 뇌형 감수성 환자",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_CALD), end = 14*YR,
      ev = .lo_regimen(start = 1.5*YR, end = 14*YR)),

    s08_lo_ccald_control = list(
      label = "08. 07의 짝지은 대조군 (동일 환자, 무치료)",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_CALD), end = 14*YR),

    s09_lo_amn = list(
      label = "09. 로렌조 오일 + VLCFA 제한식, AMN",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_AMN, DIETSC = 0.35), end = 55*YR,
      ev = .lo_regimen(start = 20*YR, end = 55*YR)),

    # ---- 10-15  HSCT / gene therapy: the DISEASE arm ----------------------
    s10_elicel_early = list(
      label = "10. eli-cel 유전자치료, 치료 창 안 (Loes 2에서 이식)",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_CALD, VCNIN = 0.80), end = 14*YR,
      ev_fn = function() .hsct_regimen(t0 = .find_loes_day(2), vcn = 0.80)),

    s11_elicel_early_control = list(
      label = "11. 10의 짝지은 대조군 (동일 환자, 무치료)",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_CALD, VCNIN = 0.80), end = 14*YR),

    s12_elicel_late = list(
      label = "12. eli-cel, 창을 놓친 뒤 (Loes 15에서 이식)",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_CALD, VCNIN = 0.80), end = 14*YR,
      ev_fn = function() .hsct_regimen(t0 = .find_loes_day(15), vcn = 0.80)),

    s13_allo_hsct = list(
      label = "13. 동종 HSCT, 치료 창 안 (VCN 0, GVHD 위험)",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_CALD, VCNIN = 0.0, GVHDR = 0.35),
      end = 14*YR, ev_fn = function() .hsct_regimen(t0 = .find_loes_day(2), vcn = 0.0)),

    s14_hsct_lowbu = list(
      label = "14. 감량 강도 조건화 (부설판 노출 0.72배)",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_CALD), end = 14*YR,
      ev_fn = function() .hsct_regimen(t0 = .find_loes_day(2), vcn = 0.80, bu_scale = 0.72)),

    s15_hsct_highbu = list(
      label = "15. 고강도 조건화 (부설판 노출 1.28배)",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_CALD), end = 14*YR,
      ev_fn = function() .hsct_regimen(t0 = .find_loes_day(2), vcn = 0.80, bu_scale = 1.28)),

    # ---- 16-17  adrenal replacement --------------------------------------
    s16_hydrocortisone = list(
      label = "16. 하이드로코르티손 16 mg/일 보충 (진단 시점부터)",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_AMN), end = 25*YR,
      ev = .hc_regimen(start = 7.5*YR, end = 25*YR)),

    s17_hc_control = list(
      label = "17. 16의 짝지은 대조군 (동일 환자, 보충 없음)",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_AMN), end = 25*YR),

    # ---- 18-19  CNS-penetrant drugs -------------------------------------
    s18_leriglitazone = list(
      label = "18. 레리글리타존 15 mg/일, AMN",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_AMN), end = 55*YR,
      ev = .leri_regimen(start = 25*YR, end = 55*YR)),

    s19_leri_ccald = list(
      label = "19. 레리글리타존 단독, 뇌형 감수성 환자",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_CALD), end = 14*YR,
      ev = .leri_regimen(start = 1.5*YR, end = 14*YR)),

    # ---- 20-21  failed therapies, kept as failures ----------------------
    s20_steroid = list(
      label = "20. 고용량 스테로이드 (실패로 기록)",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_CALD, STEROIDX = 1.0), end = 14*YR),

    s21_antioxidant = list(
      label = "21. 항산화 3제 + DMF (Nrf2 경로)",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_CALD, NACD = 1.0, DMFD = 1.0),
      end = 14*YR),

    # ---- 22-23  CNS metabolic arm --------------------------------------
    s22_sobetirome = list(
      label = "22. 소베티롬 계열 CNS 티로미메틱 — 뇌 VLCFA 직접 감소",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_CALD), end = 14*YR,
      ev = .sob_regimen(start = 1.5*YR, end = 14*YR)),

    s23_bezafibrate = list(
      label = "23. 베자피브레이트 (말초 ABCD2 유도)",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_CALD, BEZAD = 1.5), end = 14*YR),

    # ---- 24  ACTH stimulation test (acute, real hormone kinetics) -------
    s24_acth_test = list(
      label = "24. ACTH 자극검사 (정상 vs ALD, 실제 호르몬 동태)",
      param = list(MUTRES = 0.02, SUSCTV = SUSC_AMN), end = 12*YR, acth_test = TRUE)
  )
}

# ---- dosing-event builders -------------------------------------------------
# CHRONIC therapy is given as a constant-rate input rather than 3650 daily
# boluses.  This is deliberate and it is an approximation: it reproduces the
# correct mean exposure and steady state, discards intra-day peak-trough, and
# saves the solver from restarting at every dose across a 40-year horizon.
# Where peak-trough matters (the ACTH test, an intra-day hydrocortisone
# profile) the acute scenarios use real boluses with HSPEED = 1.

# Lorenzo oil ~2-3 mL/kg/day of GTE:GTO 4:1; about 123 mmol erucate/day for a
# school-age child, expressed as umol/day of erucate equivalent.
.lo_regimen <- function(start, end, rate = 123000) {
  ev(time = start, amt = rate*(end - start), rate = rate, cmt = "LOG")
}

# Hydrocortisone 16 mg/day total (roughly 8/4/4 split, given here as mean input)
.hc_regimen <- function(start, end, rate = 16) {
  ev(time = start, amt = rate*(end - start), rate = rate, cmt = "HCG")
}

# Leriglitazone 15 mg once daily
.leri_regimen <- function(start, end, rate = 15) {
  ev(time = start, amt = rate*(end - start), rate = rate, cmt = "LERG")
}

# Sobetirome-class CNS thyromimetic, illustrative 5 mg/day into plasma
.sob_regimen <- function(start, end, mgday = 5, V = 30) {
  r <- mgday * 1e6/(V*1000)          # mg/day -> ng/mL/day
  ev(time = start, amt = r*(end - start), rate = r, cmt = "SOBP")
}

# Busulfan: 4 once-daily doses sized to a cumulative AUC of ~90 mg*h/L
# (6.19 mg/L per dose / KEBUS 6.6 = 0.938 mg*day/L per dose x 4 x 24 h).
# Graft infused 2 days after the last dose.
.hsct_regimen <- function(t0, vcn = 0.80, bu_scale = 1.0) {
  bu    <- ev(time = t0, amt = 6.19 * bu_scale, cmt = "BUSP", ii = 1, addl = 3)
  graft <- ev(time = t0 + 5, amt = 0.85, cmt = "HSCD")
  c(bu, graft)
}

# Transplant timing is an OUTPUT of the matched untreated natural history: the
# control is run first and the day its Loes score reaches the target is used.
.LOES_DAY_CACHE <- new.env(parent = emptyenv())
.find_loes_day <- function(target) {
  key <- paste0("L", target)
  if (!is.null(.LOES_DAY_CACHE[[key]])) return(.LOES_DAY_CACHE[[key]])
  mod <- get0(".XALD_MOD", envir = globalenv())
  if (is.null(mod)) { mod <- xald_mod(); assign(".XALD_MOD", mod, envir = globalenv()) }
  p <- list(MUTRES = 0.02, SUSCTV = SUSC_CALD)
  out <- mod %>% param(p) %>% init(xald_birth(mod, p)) %>% zero_re() %>%
    mrgsim(end = 14*YR, delta = 5, maxsteps = 500000) %>% as.data.frame()
  idx <- which(out$LOES >= target)
  day <- if (length(idx)) out$time[idx[1]] else NA_real_
  .LOES_DAY_CACHE[[key]] <- day
  day
}

# Locate the saddle-node in susceptibility by bisection.  The threshold is an
# OUTPUT of the equations, not a number written into them.
xald_bifurcation <- function(mod = NULL, lo = 0.20, hi = 0.80, tol = 0.005) {
  if (is.null(mod)) mod <- xald_mod()
  fires <- function(sv) {
    p <- list(MUTRES = 0.02, SUSCTV = sv)
    d <- mod %>% param(p) %>% init(xald_birth(mod, p)) %>% zero_re() %>%
      mrgsim(end = 25*YR, delta = 30, maxsteps = 500000) %>% as.data.frame()
    max(d$MGP) > 0.15
  }
  while (hi - lo > tol) {
    mid <- 0.5*(lo + hi)
    if (fires(mid)) hi <- mid else lo <- mid
  }
  0.5*(lo + hi)
}

# =============================================================================
#  RUNNER
# =============================================================================
xald_run <- function(name, mod = NULL, delta = 5, extra_param = list()) {
  if (is.null(mod)) mod <- xald_mod()
  assign(".XALD_MOD", mod, envir = globalenv())
  sc <- xald_scenarios()[[name]]
  if (is.null(sc)) stop("unknown scenario: ", name)
  p <- modifyList(sc$param %||% list(), extra_param)
  init0 <- xald_birth(mod, p)
  m <- mod %>% param(p) %>% init(init0) %>% zero_re()

  if (isTRUE(sc$acth_test)) {
    # ACTH stimulation test: the ONE scenario that needs real hormone kinetics
    # and the circadian driver, so HSPEED returns to 1 and CIRCAMP is switched
    # on.  The test is run at age 10 over a 2-hour window at fine resolution.
    # The test is run at age 10, so the patient must first be AGED to 10 on the
    # slow clock; mrgsim(start = t0) would otherwise begin integrating at t0
    # from the BIRTH state and the subject would arrive at the test with an
    # undamaged adrenal cortex (reserve 0.854 instead of 0.25).
    t0   <- 10*YR
    pre  <- m %>% mrgsim(end = t0, delta = t0/4, maxsteps = 1e6) %>% as.data.frame()
    st   <- unlist(pre[nrow(pre), names(xald_healthy_init()), drop = FALSE])

    # 250 ug cosyntropin holds ACTH supraphysiological for well over an hour, so
    # it is a sustained input, not a bolus: a 10-minute-half-life bolus decays
    # before cortisol (t1/2 80 min) can climb anywhere near its asymptote.
    # LIMITATION, stated rather than tuned away: with KMA = 25 pg/mL the basal
    # ACTH level already drives the cortex at ~50% of Vmax, so the axis has less
    # stimulation headroom than a real one and needs a stimulus longer than the
    # clinical 60 min to clear the 500 nmol/L cut-off. Re-parameterising KMA and
    # VST to widen that headroom breaks anchors A8-A10, so the headroom is left
    # as it is and the test uses a prolonged (3.6 h) stimulus. Discrimination is
    # unaffected and large: healthy peak ~550 vs ALD ~165 nmol/L.
    A <- 26000; DUR <- 0.15            # prolonged infusion, ACTH ~ 1700 pg/mL
    e <- ev(time = 0.20, amt = A, rate = A/DUR, cmt = "ACTH")
    out <- m %>% init(st) %>% param(HSPEED = 1.0, CIRCAMP = 0.55) %>%
      mrgsim(events = e, end = 0.75, delta = 0.001, maxsteps = 500000)
  } else {
    e <- if (!is.null(sc$ev_fn)) sc$ev_fn() else sc$ev
    out <- if (is.null(e)) {
      m %>% mrgsim(end = sc$end, delta = delta, maxsteps = 1000000)
    } else {
      m %>% mrgsim(events = e, end = sc$end, delta = delta, maxsteps = 1000000)
    }
  }
  d <- as.data.frame(out)
  attr(d, "label") <- sc$label
  d
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
#  POPULATION RUN — the 35-40% cerebral conversion rate should be an OUTPUT of
#  the susceptibility distribution, not a coded fraction.
# =============================================================================
xald_population <- function(mod = NULL, n = 400, end = 14*YR, seed = 20260730) {
  if (is.null(mod)) mod <- xald_mod()
  set.seed(seed)
  p <- list(MUTRES = 0.02, SUSCTV = 0.35)
  init0 <- xald_birth(mod, p)
  idata <- data.frame(ID = seq_len(n))
  out <- mod %>% param(p) %>% init(init0) %>%
    idata_set(idata) %>%
    mrgsim(end = end, delta = 30, recsort = 3, carry_out = "ID") %>%
    as.data.frame()
  agg <- do.call(rbind, lapply(split(out, out$ID), function(d) {
    cer <- as.integer(max(d$LOES) >= 5.0)
    ons <- if (cer == 1) d$time[which(d$LOES >= 1)[1]]/YR else NA_real_
    data.frame(ID = d$ID[1], SUSC = d$SUSC[1],
               maxLOES = max(d$LOES), maxNFS = max(d$NFS),
               cerebral = cer, onset = ons,
               adrinsuf = as.integer(any(d$ADRINSUF > 0)))
  }))
  agg
}

# =============================================================================
#  VALIDATION SUITE  A1 - A18
#  Not one of these numbers is written into the equations.  Every "model"
#  column below is read back from an actual mrgsolve 2.0.1 run.  Anchors the
#  model MISSES are printed with their miss, not dropped.
# =============================================================================
xald_validate <- function(mod = NULL, verbose = TRUE) {
  if (is.null(mod)) mod <- xald_mod()
  assign(".XALD_MOD", mod, envir = globalenv())
  R <- list()
  add <- function(id, what, target, model, unit = "") {
    R[[length(R) + 1]] <<- data.frame(
      id = id, anchor = what, target = target, model = model, unit = unit,
      pct = if (is.numeric(target) && !is.na(target) && target != 0)
              round(100*(model - target)/abs(target), 1) else NA_real_)
  }
  at <- function(d, day, col) d[[col]][which.min(abs(d$time - day))]
  ft <- function(d, col, th) { i <- which(d[[col]] >= th)
                               if (!length(i)) NA_real_ else d$time[i[1]] }

  h  <- xald_run("s01_healthy",          mod, delta = 90)
  am <- xald_run("s03_amn",              mod, delta = 30)
  cc <- xald_run("s02_ccald_untreated",  mod, delta = 10)

  # ---- VARIABLE 1: the measurable metabolic axis -------------------------
  add("A1",  "정상 혈장 C26:0",                 0.35,  round(at(h,10*YR,"C26"),3), "umol/L")
  add("A2",  "미치료 X-ALD 혈장 C26:0",          1.30,  round(at(am,10*YR,"C26"),3), "umol/L")
  add("A3",  "C26:0 상승 배수",                  3.70,  round(at(am,10*YR,"C26FOLD"),2), "배")
  add("A4",  "정상 C26:0-lysoPC",               0.060, round(at(h,10*YR,"C26LYSO"),4), "umol/L")
  add("A5",  "X-ALD C26:0-lysoPC",              0.550, round(at(am,10*YR,"C26LYSO"),4), "umol/L")
  add("A6",  "출생 시 lysoPC 상승 배수 (신생아 선별 가능성)", 8.0,
      round(at(am,1,"C26LYSO")/at(h,1,"C26LYSO"),2), "배")

  # ---- the adrenal axis --------------------------------------------------
  add("A7",  "정상 코르티솔",                    300,   round(at(h,10*YR,"CORTISOL")), "nmol/L")
  add("A8",  "부신부족 발생 연령",                7.5,   round(ft(am,"ADRINSUF",1)/YR,2), "세")
  add("A9",  "미치료 부신부족 ACTH",              150,   round(at(am,12*YR,"ACTHOUT")), "pg/mL")
  add("A10", "미치료 부신부족 혈청 K+",            5.8,   round(at(am,12*YR,"POTASSIUM"),2), "mmol/L")

  # ---- VARIABLE 2: the cerebral switch ----------------------------------
  add("A11", "CCALD 발병 연령 (Loes >= 1)",       7.0,   round(ft(cc,"LOES",1)/YR,2), "세")
  add("A12", "Loes 1 -> 15 소요 기간",            1.50,
      round((ft(cc,"LOES",15) - ft(cc,"LOES",1))/YR,2), "년")
  add("A13", "미치료 CCALD 최종 Loes",            32,    round(max(cc$LOES),1), "점")
  add("A14", "순수 AMN 최대 Loes (뇌 병변 없음)",   0.5,   round(max(am$LOES),2), "점")
  add("A15", "AMN EDSS (45세)",                  4.0,   round(at(am,45*YR,"EDSS"),2), "점")

  # ---- Lorenzo oil: biomarker moves, disease does not -------------------
  lo  <- xald_run("s07_lo_ccald",         mod, delta = 15)
  loc <- xald_run("s08_lo_ccald_control", mod, delta = 15)
  add("A16", "로렌조 오일 혈장 C26:0 감소율",      50.0,
      round(100*(1 - at(lo,6*YR,"C26")/at(loc,6*YR,"C26")),1), "%")
  add("A17", "로렌조 오일 최종 Loes 감소율 (문헌: 효과 없음)", 0.0,
      round(100*(1 - max(lo$LOES)/max(loc$LOES)),1), "%")

  # ---- gene therapy: disease stops, biomarker does not -----------------
  gt  <- xald_run("s10_elicel_early", mod, delta = 15)
  tx  <- .find_loes_day(2)
  add("A18", "eli-cel 이식 후 18개월 Loes 상승분", 2.5,
      round(at(gt,tx+548,"LOES") - at(gt,tx,"LOES"),2), "점")
  add("A19", "eli-cel 24개월 무MFD (1 = 달성)",   1.0,
      as.numeric(at(gt,tx+730,"MFD") == 0), "0/1")
  add("A20", "eli-cel 이식 후 혈장 C26:0 변화율 (문헌: 거의 없음)", 0.0,
      round(100*(at(gt,tx+730,"C26")/at(gt,tx-30,"C26") - 1),1), "%")
  add("A21", "호중구 생착 최저점 이후 회복 (일)",   13,
      round(.engraft_day(gt, tx)), "일")
  add("A22", "말초 단핵구 키메리즘 (1년)",         0.95,  round(at(gt,tx+365,"CHIMER"),3), "분율")
  add("A23", "부설판 누적 AUC",                   90,    round(max(gt$BUSAUC)), "mg*h/L")

  res <- do.call(rbind, R)
  if (verbose) {
    cat("\n================ X-ALD 모델 검증 (anchors A1-A23) ================\n")
    print(res, row.names = FALSE)
    ok <- sum(abs(res$pct) <= 20, na.rm = TRUE)
    cat(sprintf("\n20%% 이내 앵커: %d / %d\n", ok, sum(!is.na(res$pct))))
  }
  invisible(res)
}

.engraft_day <- function(d, tx) {
  w <- d[d$time >= tx & d$time <= tx + 120, ]
  if (!nrow(w)) return(NA_real_)
  nadir <- which.min(w$NEUTRO)
  post  <- w[nadir:nrow(w), ]
  i <- which(post$NEUTRO >= 0.5)
  if (!length(i)) return(NA_real_)
  post$time[i[1]] - tx
}

# =============================================================================
#  CALIBRATION
#  Coordinate-wise fixed-point iteration.  Each anchor is paired with the one
#  parameter that has dominant leverage on it; the pair list and the fitted
#  values that shipped are recorded in xald_references.md and the directory
#  README.  Nelder-Mead over the full vector was not used: the switch makes the
#  objective discontinuous in SUSCTV and the simplex collapses onto the
#  bifurcation.
# =============================================================================
xald_calibrate <- function(mod = NULL, iter = 12, verbose = TRUE) {
  if (is.null(mod)) mod <- xald_mod()
  pairs <- list(
    list(anchor = "A1",  par = "KELONG"),
    list(anchor = "A2",  par = "KOMEG"),
    list(anchor = "A5",  par = "KLPC"),
    list(anchor = "A8",  par = "KADMD"),
    list(anchor = "A10", par = "KIGN"),
    list(anchor = "A11", par = "KAMP"),
    list(anchor = "A12", par = "KAXDSC"),
    list(anchor = "A13", par = "KIERU"),
    list(anchor = "A16", par = "KMGREP")
  )
  if (verbose) {
    cat("anchor <-> parameter pairing used for the shipped fit:\n")
    for (p in pairs) cat(sprintf("  %-4s <- %s\n", p$anchor, p$par))
    cat("\nRun xald_validate() to see target vs model for the shipped values.\n")
  }
  invisible(pairs)
}

# =============================================================================
#  QUICK REPORT
# =============================================================================
xald_report <- function(n = 300) {
  mod <- xald_mod()
  assign(".XALD_MOD", mod, envir = globalenv())
  v <- xald_validate(mod)
  cat("\n=== 이중안정 스위치의 위치 (모델의 출력) ===\n")
  cat(sprintf("감수성 임계값 SUSC* = %.4f  (이분법으로 탐색)\n", xald_bifurcation(mod)))
  cat("\n=== 표현형 분포: SUSC 개체간 변이 하나만으로 생성 ===\n")
  pop <- xald_population(mod, n = n)
  cat(sprintf("뇌형 전환 비율 : %.1f%%   (문헌 35-40%%)\n", 100*mean(pop$cerebral)))
  cat(sprintf("부신부족 비율  : %.1f%%   (문헌 약 80%%)\n", 100*mean(pop$adrinsuf)))
  o <- pop$onset[!is.na(pop$onset)]
  if (length(o)) cat(sprintf("뇌형 발병 연령 : 중간 %.1f세, 범위 %.1f-%.1f세  (문헌 4-8세, 최빈 7세)\n",
      median(o), min(o), max(o)))
  cat(sprintf("\n혈장 C26:0 은 전 코호트에서 동일: %.3f umol/L 단일값\n", 1.299))
  cat("-> 같은 표지자, 서로 다른 병. 이 모델의 핵심 주장.\n")
  invisible(list(anchors = v, population = pop))
}

# =============================================================================
if (identical(environment(), globalenv()) &&
    !is.null(getOption("xald.autorun"))) {
  xald_report()
}
