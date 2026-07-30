## =============================================================================
##  zes_mrgsolve_model.R
##  Zollinger-Ellison Syndrome (ZES) / gastrinoma
##  Quantitative Systems Pharmacology model for mrgsolve
##
##  59 ODE compartments.  Time unit = HOURS, because this disease forces a
##  6-minute gastrin half-life, a 50-hour proton-pump half-life, a 3-week
##  parietal-cell trophic time constant and a 2-year progression-free survival
##  to live in the same system.
##
##  ---------------------------------------------------------------------------
##  THE ONE STRUCTURAL COMMITMENT
##  ---------------------------------------------------------------------------
##
##  Gastric acid output is not a state variable.  It is a PRODUCT:
##
##      BAO  =  BAOCAP  x  PCM  x  PUMPA  x  ACTP
##                         ^^^     ^^^^^     ^^^^
##                       FACTOR 1  FACTOR 2  FACTOR 3
##                       parietal  fraction  per-pump activation by the
##                       cell mass of pumps  secretagogue drive at
##                       (weeks)   at the    CCK2R / H2 / M3 / CaSR
##                                 membrane  (minutes-hours)
##                                 (hours-days)
##
##  Hypergastrinaemia raises FACTOR 1 (trophic, over weeks) and FACTOR 3
##  (acute, over minutes).  Every acid-suppressing drug ever licensed acts on
##  FACTOR 2 and only on FACTOR 2.  That asymmetry -- not a special resistance
##  of gastrinoma parietal cells -- is what this file is built to express, and
##  it is why the following are OUTPUTS here rather than coded rules:
##
##    * the ZES maintenance PPI dose is several times the reflux dose.
##      Removing the same FRACTION of the pump pool leaves a proportionally
##      larger ABSOLUTE residual when factors 1 and 3 are both raised;
##
##    * the DOSING INTERVAL fails before the daily dose does.  KSYNP (pump
##      synthesis) is up-regulated by the trophic gastrin signal, so the pool
##      refills FASTER in ZES, so the trough arrives earlier.  Nocturnal acid
##      breakthrough is a consequence of pump-pool kinetics plus the circadian
##      histamine peak, not a separate phenomenon;
##
##    * an H2-receptor antagonist is out-competed.  DRIVE is a PRODUCT of
##      three potentiating arms, so blocking one leaves the other two
##      multiplying;  a covalent pump blocker sits downstream of all three and
##      cannot be out-competed;
##
##    * a somatostatin analogue, cytoreductive surgery or PRRT is additive
##      with a PPI rather than redundant, because those agents lower SPHEN and
##      therefore factors 1 and 3, which no PPI touches;
##
##    * rebound hypersecretion on withdrawal is the enlarged PUMPI pool
##      discharging onto the enlarged PCM.
##
##  ---------------------------------------------------------------------------
##  TWO SIGN INVERSIONS, WRITTEN AS SIGN INVERSIONS
##  ---------------------------------------------------------------------------
##
##  1. SECRETIN.  SECN = 1/(1 + CSEC/KSECN) multiplies the ANTRAL G-cell
##     target (inhibition).  SECT = 1 + ESECT*CSEC/(KSECT+CSEC) multiplies the
##     TUMOUR secretion rate (stimulation).  Same compartment, opposite sign,
##     two targets.  The secretin test is therefore an output of the model, and
##     so is its behaviour in the two situations that matter clinically -- a
##     healthy subject and a healthy subject on a PPI.
##
##  2. LUMINAL ACID.  BRK = 1/(1+(HC/HBRK50)^HBRKN) gates the antral G cell
##     and NOTHING gates the tumour.  That single missing edge -- not a higher
##     secretion rate per cell -- is what makes the hypersecretion unremitting,
##     and the SAME missing edge is why a PPI raises gastrin two- to threefold
##     in a person WITHOUT a gastrinoma (BRK rises towards 1) and thereby
##     manufactures the commonest false-positive diagnosis in this disease.
##
##  ---------------------------------------------------------------------------
##  THE MEN1 ARM IS A LOOP, NOT A LABEL
##  ---------------------------------------------------------------------------
##
##  GLANDM (parathyroid functional mass) -> PTH -> CAION -> CAF, a CaSR-mediated
##  secretagogue multiplier on BOTH tumour and antral gastrin release, and ->
##  the calcium arm of DRIVE.  Nothing in this file says "parathyroidectomy
##  lowers gastrin".  Scenario 17 sets GLANDM back to 1.0 and the model says it.
##
##  ---------------------------------------------------------------------------
##  WHAT IS AN INPUT AND WHAT IS AN OUTPUT
##  ---------------------------------------------------------------------------
##
##  INPUTS (things a user sets):  tumour burden TUM0/MET0, grade via GRADEF,
##  SSTR2 density SSTR2D, MEN1 status (GLANDM, MENFLG), CYP2C19 phenotype CYPF,
##  renal function RENF, NSAID / H. pylori exposure, drug doses and schedules.
##
##  OUTPUTS (things the model must produce and is never told):  basal and
##  maximal acid output, the BAO/MAO ratio, fasting serum gastrin, intragastric
##  pH and its 24-h holding time, the secretin-test increment and its false
##  positive, the required maintenance dose, nocturnal acid breakthrough, ulcer
##  and oesophageal injury, diarrhoea, withdrawal rebound, ECL hyperplasia and
##  gastric-NET emergence, tumour progression, renal dosimetry, and the three
##  headroom fractions.
##
##  ---------------------------------------------------------------------------
##  ONE DELIBERATE PHENOMENOLOGY, DECLARED
##  ---------------------------------------------------------------------------
##
##  Intragastric pH is NOT computed as -log10[H+] of the luminal pool.  Doing
##  that makes the achlorhydric limit unreachable: with first-order emptying
##  the residual H+ pool can never fall by the four orders of magnitude that
##  separate pH 1.7 from pH 6.3, so a fully blocked stomach comes out at pH 3.
##  The luminal contents are a buffered system, so pH follows a TITRATION
##  CURVE, and the model uses a one-line Hill titration curve
##
##      pH = PHMIN + (PHMAX - PHMIN) / (1 + (HC/HC50)^PHN)
##
##  fitted to the two endpoints that are measured directly -- fasting normal
##  pH 1.7 and complete achlorhydria pH ~6.3.  This is the single frankly
##  empirical relation in the file and it is flagged here rather than buried.
##
##  ---------------------------------------------------------------------------
##  UNITS
##  ---------------------------------------------------------------------------
##    time                    hours
##    acid output (BAO/MAO)   mEq H+ per hour
##    luminal H+ (HLUM)       mEq          gastric volume: L
##    intragastric pH         pH units
##    gastrin G17, G34        pg/mL (G17-equivalent immunoreactivity)
##    PCM, ECL, HIS, GANT     fold of the healthy reference (1.0 = normal)
##    PUMPI/A/B/R             fraction of the healthy TOTAL pump pool
##    tumour TUM, MET, GNET   cm^3
##    drug amounts            mg    (PRRT: GBq)   absorbed dose: Gy
##    calcium, magnesium      mmol/L            PTH: pmol/L
##    injury indices          0-1
##
##  Author: QSP Disease Model Library (Claude Code Routine)
##  Licence: see repository LICENSE.  EDUCATIONAL / RESEARCH USE ONLY --
##  not validated for clinical or regulatory use.
## =============================================================================

library(mrgsolve)

zes_code <- r"---(
$PROB
# Zollinger-Ellison syndrome / gastrinoma QSP model
- BAO = BAOCAP * PCM * PUMPA * ACTP  (three independently breakable factors)
- 59 ODEs, time unit = hours

$PARAM @annotated
// ---------------------------------------------------------------- tumour ----
TUM0    :  3.0   : primary gastrinoma volume at t=0 (cm3)
MET0    :  0.0   : hepatic metastatic volume at t=0 (cm3)
TMAX    :  300   : Gompertz carrying capacity (cm3)
KG0     :  3.54e-5 : Gompertz growth constant at reference grade (1/h)
GRADEF  :  1.0   : grade multiplier on growth (Ki-67 about 5 percent = 1.0)
KGMET   :  1.35  : metastatic growth relative to primary (-)
KSECT0  :  562.45 : gastrin secretion per cm3 per unit SPHEN (pg/mL/h/cm3)
KSPH    :  0.01  : secretory-phenotype turnover rate (1/h)
KRES    :  6.0e-5: resistant-clone emergence rate under mTOR blockade (1/h)
SSTR2D  :  1.0   : SSTR2 surface density, fold of reference (-)
MGMTM   :  0.0   : MGMT promoter methylated (1 = methylated, alkylator-sensitive)
KDIN    :  0.02  : damage-signal formation rate (1/h)
KDOUT   :  0.004 : damage-signal resolution rate (1/h)
KKILL   :  1.6e-4: tumour kill per unit alkylator damage signal (1/h)
KRADK   :  0.016 : effective radiation alpha for protracted PRRT delivery (1/Gy)

// --------------------------------------------------------------- gastrin ----
KEL17   :  6.93  : G17 elimination (1/h) t1/2 = 6 min
KEL34   :  0.99  : G34 elimination (1/h) t1/2 = 42 min
KCONV   :  0.15  : peripheral G34 to G17 conversion (1/h)
FR17T   :  0.50  : fraction of tumour gastrin released as G17 (-)
FR17A   :  0.75  : fraction of antral gastrin released as G17 (-)
KSECA   :  901.01 : antral gastrin secretion scale (pg/mL/h per unit GANT)
KGANT   :  0.35  : antral G-cell output turnover rate (1/h)
HBRK50  :  1.0   : luminal H+ giving half the D-cell acid brake (mmol/L)
HBRKN   :  0.70  : Hill exponent of the acid brake (-)
EGRP    :  1.6   : meal / GRP stimulation of the antral G cell (-)
ATRF    :  1.0   : antral G-cell mass factor (atrophy < 1)
RENF    :  1.0   : renal + hepatic gastrin clearance factor (CKD < 1)
GBNORM  :  27.121 : reference potency-weighted gastrin (pg/mL-equiv)
KDG     :  60    : CCK2R KD on ECL and trophic targets (pg/mL-equiv)
KDGP    :  400   : CCK2R KD on the parietal cell, low affinity (pg/mL-equiv)
PENTA   :  0.0   : pentagastrin-equivalent occupancy added to CCK2R (-)
ICFLAG :  1.0   : apply the built-in healthy initial conditions (0 = use init())

// -------------------------------------------------- trophic cell masses ----
EMAXPCM :  2.0674 : maximal fold increase of parietal cell mass (-)
EC50PCM :  3.0   : trophic-drive excess giving half-maximal PCM growth (-)
KPCM    :  0.00138: parietal-mass turnover rate (1/h) t1/2 = 3 weeks
EMAXECL :  3.0   : maximal fold increase of ECL mass (-)
EC50ECL :  3.5   : trophic-drive excess giving half-maximal ECL growth (-)
KECL    :  0.00206: ECL-mass turnover rate (1/h) t1/2 = 2 weeks
KHIS    :  0.05  : histamine store turnover toward ECL mass (1/h)

// ------------------------------------------------------- pump lifecycle ----
KDEGP   :  0.01386: pump degradation (1/h) t1/2 = 50 h
GSYN    :  0.70  : maximal gastrin up-regulation of pump synthesis (-)
KGSYN   :  2.0   : trophic-drive excess for half-maximal synthesis induction (-)
KACT0   :  0.02  : constitutive trafficking of pumps to the canaliculus (1/h)
KACT1   :  1.20  : drive-dependent trafficking (1/h)
KINACT  :  0.35  : endocytic retrieval of active pumps (1/h)
KGSH    :  0.0015: glutathione-mediated reactivation of the PPI adduct (1/h)

// -------------------------------------------------------- acid and drive ---
// BAOCAP and KDRV are FITTED (see the CALIBRATION section for all nine)
BAOCAP  :  78.879 : acid capacity at PCM=1, PUMPA=1, ACTP=1 (mEq/h)
KDRV    :  45.447 : drive giving half-maximal per-pump activation (-)
AG      :  3.5   : weight of the direct parietal CCK2R arm (-)
AH      :  6.0   : weight of the histamine H2 arm (-)
ACHW    :  2.0   : weight of the cholinergic M3 arm (-)
ACAW    :  2.5   : weight of the ionised-calcium arm (per mmol/L above 1.25)
ACHBAS  :  0.30  : basal vagal cholinergic tone (-)
ACHMEAL :  1.20  : meal-evoked cholinergic tone (-)
CIRCA   :  0.35  : circadian amplitude of histamine drive (-)
CIRCT   :  1.0   : clock hour of peak nocturnal drive (h)
KHREL   :  1.0   : histamine release scale (-)
KH2R    :  1.0   : histamine signal for half-maximal H2 activation (-)
KRH2    :  0.02  : H2-receptor density turnover (1/h)
TACHM   :  0.90  : maximal H2-receptor up-regulation (tachyphylaxis) (-)
KTACH   :  0.15  : H2RA concentration for half-maximal up-regulation (mg/L)
DRVREF  :  7.12  : reference basal drive of a healthy stomach (-)

// ------------------------------------------------- lumen, meals, volume ----
MEALS   :  1.0   : meals on (1) or fasted study (0)
BUFMEAL : 107.10  : titratable buffer delivered per meal (mEq)
VG0     :  0.040 : fasting residual gastric volume (L)
KVACID  :  256   : mEq H+ per hour per litre of secretory volume retained (mEq/h/L)
KVBUF   :  0.008 : gastric volume added per mEq of meal buffer (L/mEq)
KEMPT   :  2.6218 : gastric emptying of buffer (1/h)
KEMPTH  :  1.6   : gastric emptying of H+ (1/h)
KNEUT   :  1.0   : buffer-limited neutralisation capacity (1/h per mEq buffer)
KHN     :  1.5   : luminal H+ at which neutralisation is half-saturated (mEq)
HCO3    :  1.2   : salivary + duodenogastric alkali flux (mEq/h)
KHSAT   :  0.15  : luminal H+ at which alkali neutralisation is half-saturated (mEq)
PHMIN   :  1.20  : titration-curve floor of intragastric pH (-)
PHMAX   :  7.00  : titration-curve ceiling (complete achlorhydria) (-)
HC50    :  7.00  : luminal H+ concentration at the curve mid-point (mmol/L)
PHN     :  1.45  : titration-curve steepness (-)

// --------------------------------------------- mucosa, injury, endpoints ---
KMUC    :  0.02  : mucosal-defence turnover (1/h)
ENSAID  :  0.55  : NSAID suppression of mucus-bicarbonate defence (-)
NSAID   :  0.0   : NSAID / aspirin exposure (0-1)
EHP     :  0.25  : H. pylori suppression of defence (-)
HPYL    :  0.0   : H. pylori co-infection (0-1)
ATHR    :  14.0  : acid output a normal mucosa tolerates (mEq/h)
KUL50   :  12.0  : excess acid load for half-maximal ulcerogenesis (mEq/h)
KULIN   :  0.012 : ulcer formation rate (1/h)
KULHEAL :  0.006 : ulcer healing rate (1/h)
KBLEED  :  1.2e-4: haemorrhage / perforation hazard per unit ulcer squared (1/h)
KESIN   :  0.006 : oesophageal injury rate (1/h)
KES50   :  18.0  : excess acid load for half-maximal oesophageal injury (mEq/h)
KESHEAL :  0.003 : oesophageal healing rate (1/h)
KDIA    :  0.06  : diarrhoea index turnover (1/h)
KDIA50  :  30.0  : excess acid load giving half-maximal diarrhoea (mEq/h)

// ------------------------------------ long-term consequences of high pH ----
KB12OUT :  1.5e-4: cobalamin store turnover (1/h) t1/2 about 190 days
PHB12   :  4.0   : pH at which B12 liberation is half-blocked (-)
NB12    :  4.0   : Hill exponent for acid-dependent B12 liberation (-)
KMGOUT  :  0.002 : magnesium turnover (1/h)
MGREF   :  0.85  : reference serum magnesium (mmol/L)
FMGFL   :  0.75  : acid-independent floor of magnesium absorption (-)
PHMG    :  4.5   : pH at which Mg absorption is half-blocked (-)
KBMD    :  3.0e-7: bone loss rate at complete achlorhydria (1/h)
KGNET   :  1.4e-5: gastric-NET formation per unit excess ECL mass (cm3/h)
MENFLG  :  0.0   : MEN1 syndrome present (1) - permissive for type-2 gastric NET
KGNETR  :  2.0e-6: gastric-NET regression rate (1/h)

// -------------------------------------------------------- MEN1 Ca-PTH loop -
GLANDM  :  1.0   : parathyroid functional mass (1 = normal, 1.8 = MEN1 pHPT)
PTHBAS  :  4.0   : reference PTH (pmol/L)
KPTH    :  0.35  : PTH turnover (1/h)
CASET   :  1.25  : calcium set-point of the parathyroid (mmol/L)
CAHILL  :  4.0   : steepness of calcium feedback on PTH (-)
KCAIN   :  0.030 : calcium influx per unit PTH (mmol/L/h per pmol/L)
KCAOUT  :  0.095 : calcium clearance (1/h)
VITDF   :  1.0   : vitamin-D-dependent intestinal calcium absorption factor (-)
ECA     :  1.60  : CaSR secretagogue effect on gastrin release (per mmol/L)

// -------------------------------------------------------------- PPI PK/PD --
KAPPI   :  1.10  : PPI absorption rate (1/h)
FPPI    :  0.55  : PPI oral bioavailability (-)
VPPI    :  30.0  : PPI central volume (L)
CLPPI0  :  21.0  : PPI clearance in a CYP2C19 normal metaboliser (L/h)
CYPF    :  1.0   : CYP2C19 clearance factor (UM 1.9 / NM 1.0 / IM 0.5 / PM 0.28)
KACTP   :  4.0   : sulfenamide formation rate from plasma prodrug (1/h per mg/L)
KMPA    :  0.035 : active-pump pool giving half-maximal drug activation (-)
KOUTCAN :  1.60  : canalicular activated-species disappearance (1/h)
KBIND   : 40.0    : covalent pump inactivation per unit activated species (1/h)

// ------------------------------------------------------------ P-CAB PK/PD --
KAVON   :  0.85  : vonoprazan absorption rate (1/h)
FVON    :  0.60  : vonoprazan oral bioavailability (-)
VVON    :  85.0  : vonoprazan central volume (L)
CLVON   :  7.6   : vonoprazan clearance (L/h)
KONV    :  32.0  : K+-competitive on-rate, ion-trapping-adjusted (1/h per mg/L)
KOFFV   :  0.55  : K+-competitive off-rate (1/h)
FRESTV  :  0.45  : accessibility of the RESTING pump pool to a P-CAB (-)

// ------------------------------------------------------------- H2RA PK/PD --
KAH2    :  1.30  : famotidine absorption rate (1/h)
FH2     :  0.42  : famotidine oral bioavailability (-)
VH2     :  90.0  : famotidine central volume (L)
CLH2    :  27.0  : famotidine clearance (L/h)
KIH2    :  0.030 : famotidine H2 inhibition constant (mg/L)

// ------------------------------------------------- somatostatin analogue ---
KRELD   :  0.0016: octreotide LAR microsphere release (1/h)
CLOCT   :  16.0  : octreotide clearance (L/h)
VOCT    :  20.0  : octreotide volume (L)
KDSST   :  0.35  : SSTR2 KD for octreotide (ng/mL)
ESSAS   :  0.62  : maximal suppression of tumour gastrin secretion by SSA (-)
ESSAE   :  0.55  : maximal suppression of ECL mass by SSA (-)
ESSAP   :  0.30  : maximal suppression of the parietal histamine arm by SSA (-)
ESSAT   :  0.50  : maximal antiproliferative effect of SSA (-)

// --------------------------------------------------------- antitumour PK ---
KAEVE   :  1.20  : everolimus absorption (1/h)
FEVE    :  0.16  : everolimus bioavailability (-)
VEVE    :  200   : everolimus volume (L)
CLEVE   :  4.6   : everolimus clearance (L/h)
IC50EVE :  0.005 : everolimus concentration for half-maximal mTORC1 block (mg/L)
EMAXEVE :  0.85  : maximal net-growth suppression by everolimus (-)
KASUN   :  0.55  : sunitinib absorption (1/h)
FSUN    :  0.50  : sunitinib bioavailability (-)
VSUN    :  2030  : sunitinib volume (L)
CLSUN   :  17.0  : sunitinib clearance, parent + SU12662 (L/h)
IC50SUN :  0.015 : sunitinib concentration for half-maximal VEGFR block (mg/L)
EMAXSUN :  0.86  : maximal net-growth suppression by sunitinib (-)
KACAP   :  1.00  : capecitabine absorption (1/h)
VCAP    :  100   : capecitabine (5-FU surrogate) volume (L)
CLCAP   :  180   : capecitabine clearance (L/h)
KATMZ   :  1.50  : temozolomide absorption (1/h)
VTMZ    :  30.0  : temozolomide volume (L)
CLTMZ   :  11.0  : temozolomide clearance (L/h)
KALK    :  2.6   : alkylator damage per unit (CAP x TMZ) exposure (-)
EMGMT   :  0.72  : loss of alkylator damage when MGMT is unmethylated (-)
KANET   :  1.40  : netazepide absorption (1/h)
VNET    :  60.0  : netazepide volume (L)
CLNET   :  30.0  : netazepide clearance (L/h)
KINET   :  0.020 : netazepide CCK2R inhibition constant (mg/L)
KELSEC  :  13.86 : secretin elimination (1/h) t1/2 = 3 min
VSEC    :  5.0   : secretin distribution volume (L)
KSECN   :  0.60  : secretin concentration inhibiting the antral G cell (U/L)
KSECT   :  0.35  : secretin concentration stimulating the tumour (U/L)
ESECT   :  3.0   : maximal secretin stimulation of tumour gastrin release (-)

// ----------------------------------------------------------------- PRRT ----
KELPRC  :  0.55  : 177Lu-DOTATATE plasma disappearance (1/h)
KUPT    :  0.055 : tumour uptake rate per unit SSTR2 (1/h)
KOUTT   :  0.0055: tumour washout plus physical decay (1/h)
KUPK    :  0.030 : renal cortical uptake (1/h)
KOUTK   :  0.010 : renal washout plus decay (1/h)
AALYS   :  1.0   : amino-acid renal protection on (1) or off (0)
FPROT   :  0.55  : fractional reduction of renal uptake by amino acids (-)
SDOSET  :  0.167 : absorbed dose per GBq-h of tumour activity (Gy/GBq/h)
SDOSEK  :  0.279 : absorbed dose per GBq-h of renal activity (Gy/GBq/h)

$CMT @annotated
TUM    : primary gastrinoma volume (cm3)
MET    : hepatic metastatic volume (cm3)
RESCL  : resistant clone fraction (-)
SPHEN  : gastrin secretion per unit tumour volume (fold)
DAMT   : alkylator damage signal (-)
G17    : little gastrin, plasma (pg/mL)
G34    : big gastrin, plasma (pg/mL)
GANT   : antral G-cell output (fold of normal)
ECL    : ECL cell mass (fold of normal)
HIS    : mucosal histamine store (fold of normal)
PCM    : parietal cell mass (fold of normal)
PUMPI  : inactive tubulovesicular pump pool (fraction of normal total)
PUMPA  : active canalicular pump pool (fraction of normal total)
PUMPB  : covalently blocked pump pool (fraction of normal total)
PUMPR  : reversibly K+-competitively blocked pump pool (fraction)
HLUM   : luminal H+ (mEq)
BUF    : intragastric titratable buffer (mEq)
MUCUS  : mucus-bicarbonate defence index (0-1)
ULCD   : duodenal ulcer burden index (0-1)
ESOPH  : oesophageal injury index (0-1)
DIARR  : diarrhoea index (0-1)
BLEED  : cumulative haemorrhage / perforation hazard (-)
B12    : cobalamin store (fold of normal)
MGS    : serum magnesium (mmol/L)
BMD    : bone mineral density (fold of normal)
GNET   : type-2 gastric NET burden (cm3)
PTH    : parathyroid hormone (pmol/L)
CAION  : ionised calcium (mmol/L)
RH2    : H2-receptor density (fold of normal)
PPIA   : PPI enteric-coated gut depot (mg)
PPIC   : PPI plasma prodrug (mg)
PPICAN : canalicular activated sulfenamide (arbitrary units)
VONA   : vonoprazan gut depot (mg)
VONC   : vonoprazan plasma (mg)
H2A    : famotidine gut depot (mg)
H2C    : famotidine plasma (mg)
OCTD   : octreotide LAR intramuscular depot (mg)
OCTC   : octreotide plasma (microgram)
EVEA   : everolimus gut depot (mg)
EVEC   : everolimus plasma (mg)
SUNA   : sunitinib gut depot (mg)
SUNC   : sunitinib plasma (mg)
CAPC   : capecitabine / 5-FU plasma (mg)
TMZC   : temozolomide plasma (mg)
NETA   : netazepide gut depot (mg)
NETC   : netazepide plasma (mg)
SECRC  : secretin plasma (U)
PRRTC  : 177Lu-DOTATATE plasma activity (GBq)
PRRTT  : 177Lu-DOTATATE tumour-bound activity (GBq)
PRRTK  : 177Lu-DOTATATE renal cortical activity (GBq)
DOSET  : cumulative absorbed tumour dose (Gy)
DOSEK  : cumulative absorbed renal dose (Gy)
AUCBAO : cumulative acid output (mEq)
AUCACD : cumulative acid load above the treatment target (mEq)
TPH4   : cumulative time with intragastric pH above 4 (h)
TBAO10 : cumulative time with BAO below 10 mEq/h (h)
HDPCM  : counterfactual cumulative acid output with a NORMAL parietal mass (mEq)
HDDRV  : counterfactual cumulative acid output with a NORMAL drive (mEq)
HDPMP  : counterfactual cumulative acid output with NO pump inhibition (mEq)

$GLOBAL
#define TWOPI 6.283185307179586
// non-negative floor helper (mrgsolve 2.x does not export posf)
inline double pz(double x) { return x > 0.0 ? x : 0.0; }
// values shared between $MAIN, $ODE and $TABLE
double gBAO, gPH, gFSG, gGB, gDRIVE, gACTP, gTD, gPOOL, gSSA, gHOCC,
       gMEAL, gCIRC, gBAOP, gBAOD, gBAOU, gSGE, gSGP, gBRK, gKSYN,
       gCPPI, gCVON, gCH2, gCOCT, gCEVE, gCSUN, gALOAD, gKACT, gACHT,
       gGANTS, gEPCM, gESEC, gEECL, gALK, gDRT, gDRK, gHC;

$MAIN
// between-patient variability enters the three gains that matter most:
// the parietal trophic gain, the tumour secretion rate, the ECL trophic gain
gEPCM = EMAXPCM * exp(ETA(1));
gESEC = KSECT0  * exp(ETA(2));
gEECL = EMAXECL * exp(ETA(3));

// The _0 assignments below are the HEALTHY reference initial condition.  They
// are guarded, because in mrgsolve a CMT_0 assignment in $MAIN silently
// overrides whatever init() supplied -- which would throw away the run-in
// state and make every trophic loop start from a healthy stomach again.
if (ICFLAG > 0.5) {
TUM_0    = TUM0;
MET_0    = MET0;
RESCL_0  = 0.0;
SPHEN_0  = 1.0;
DAMT_0   = 0.0;
// gastrin and the pump pool are initialised at the HEALTHY fasting steady
// state; the chronic ZES phenotype is not an initial condition, it is what
// the trophic loops settle to (see zes_init() in the R section)
G17_0    = 10.0;
G34_0    = 19.4;
GANT_0   = 0.140;
ECL_0    = 1.0;
HIS_0    = 1.0;
PCM_0    = 1.0;
PUMPI_0  = 0.6452;
PUMPA_0  = 0.3548;
PUMPB_0  = 0.0;
PUMPR_0  = 0.0;
HLUM_0   = 1.20;
BUF_0    = 0.0;
MUCUS_0  = 1.0;
ULCD_0   = 0.0;
ESOPH_0  = 0.0;
DIARR_0  = 0.0;
BLEED_0  = 0.0;
B12_0    = 1.0;
MGS_0    = MGREF;
BMD_0    = 1.0;
GNET_0   = 0.0;
PTH_0    = PTHBAS * GLANDM;
CAION_0  = 1.25;
RH2_0    = 1.0;
}

$ODE
// ===========================================================================
// 0 - clock: three smooth meal bumps and the circadian histamine peak
// ===========================================================================
double tod = fmod(SOLVERTIME, 24.0);
gMEAL = MEALS * ( exp(-pow(tod -  8.0, 2.0) / 0.605)
                + exp(-pow(tod - 13.0, 2.0) / 0.605)
                + exp(-pow(tod - 19.0, 2.0) / 0.605) );
gCIRC = 1.0 + CIRCA * cos(TWOPI * (tod - CIRCT) / 24.0);

// ===========================================================================
// 1 - gastrin: species, potency weighting, receptor occupancy
// ===========================================================================
gGB  = G17 + G34 / 6.0;                        // potency-weighted bioactivity
double NETFR = 1.0 + (NETC / VNET) / KINET;    // competitive CCK2R blockade
gSGE = gGB / (gGB + KDG  * NETFR) + PENTA;     // ECL / trophic occupancy
gSGP = gGB / (gGB + KDGP * NETFR) + PENTA;     // parietal occupancy
if (gSGE > 0.98) gSGE = 0.98;
if (gSGP > 0.98) gSGP = 0.98;
double SGE0 = GBNORM / (GBNORM + KDG);
gTD  = gSGE / SGE0;                            // normalised trophic drive
gFSG = G17 + G34;                              // assay read-out

// ===========================================================================
// 2 - DRIVE (FACTOR 3): three potentiating arms plus calcium
// ===========================================================================
gCOCT = OCTC / VOCT;                                          // ng/mL
gSSA  = SSTR2D * gCOCT / (gCOCT + KDSST);                     // SSTR2 occupancy
gCH2  = H2C / VH2;
double HREL = KHREL * HIS * gTD * gCIRC;                      // histamine signal
gHOCC = HREL / (HREL + KH2R * (1.0 + gCH2 / KIH2) / RH2);
gHOCC = gHOCC * (1.0 - ESSAP * gSSA);                         // SST2-Gi on parietal
gACHT = ACHBAS + ACHMEAL * gMEAL;
gDRIVE = (1.0 + AG * gSGP) * (1.0 + AH * gHOCC)
       * (1.0 + ACHW * gACHT) * (1.0 + ACAW * pz(CAION - 1.25));
gACTP = gDRIVE / (gDRIVE + KDRV);                             // per-pump activation

// ===========================================================================
// 3 - pump lifecycle (FACTOR 2) - the only factor a drug reaches
// ===========================================================================
// EXC is the trophic-drive EXCESS.  It is deliberately allowed to go negative
// (floored at -0.9) rather than rectified with pz().  Rectifying it biases the
// trophic loops upward, because gTD oscillates around 1 with the meal cycle
// and only the positive half of the oscillation would then be counted: a
// healthy stomach drifted to a parietal mass of 1.22 and an ECL mass of 1.27
// with no tumour anywhere in the system.
double EXC = gTD - 1.0;
if (EXC < -0.9) EXC = -0.9;
gKSYN = KDEGP * (1.0 + GSYN * EXC / (KGSYN + EXC));
gKACT = KACT0 + KACT1 * gACTP;
gCPPI = PPIC / VPPI;
gCVON = VONC / VVON;
double KB     = KBIND * PPICAN;
double VON_ON = KONV * gCVON;

dxdt_PUMPI = gKSYN - KDEGP * PUMPI - gKACT * PUMPI + KINACT * PUMPA
           + KGSH * PUMPB - VON_ON * FRESTV * PUMPI;
dxdt_PUMPA = gKACT * PUMPI - KINACT * PUMPA - KDEGP * PUMPA
           - KB * PUMPA - VON_ON * PUMPA + KOFFV * PUMPR;
dxdt_PUMPB = KB * PUMPA - KDEGP * PUMPB - KGSH * PUMPB;
dxdt_PUMPR = VON_ON * (PUMPA + FRESTV * PUMPI) - KOFFV * PUMPR - KDEGP * PUMPR;
gPOOL = PUMPI + PUMPA + PUMPB + PUMPR;

// ===========================================================================
// 4 - THE PRODUCT and its three counterfactuals
// ===========================================================================
gBAO  = BAOCAP * PCM * PUMPA * gACTP;
// A: the same stomach with a NORMAL parietal cell mass
gBAOP = BAOCAP * 1.0 * PUMPA * gACTP;
// B: the same pump pool driven at the NORMAL basal secretagogue tone
double ACTP0 = DRVREF / (DRVREF + KDRV);
gBAOD = BAOCAP * PCM * PUMPA * ACTP0;
// C: no pump inhibition - the blocked pools returned to service and split
//    active:inactive in the same ratio the UNBLOCKED pumps currently show.
//    Written this way it collapses exactly onto gBAO when no drug is on
//    board (an earlier version used the trafficking equilibrium instead and
//    reported a headroom of -0.67 for an untreated patient).
double FREE_A = PUMPA + PUMPI;
double FA_NOW = (FREE_A > 1.0e-9) ? PUMPA / FREE_A : 0.0;
gBAOU = BAOCAP * PCM * (FREE_A + PUMPB + PUMPR) * FA_NOW * gACTP;

// ===========================================================================
// 5 - gastric volume, luminal H+ concentration, the titration curve and the
//     D-cell acid brake.  These are placed AFTER the product because the
//     secretory volume arrives with the acid, so VG depends on BAO.
// ===========================================================================
// gastric volume: fasting residual + meal volume + the isotonic secretory
// volume that comes WITH the acid (litres a day in severe ZES)
double VG = VG0 + KVBUF * BUF + gBAO / KVACID;
gHC  = pz(HLUM) / VG;                                          // mmol/L
gPH  = PHMIN + (PHMAX - PHMIN) / (1.0 + pow(gHC / HC50, PHN));
gBRK = 1.0 / (1.0 + pow(gHC / HBRK50, HBRKN));                 // D-cell brake

// luminal balance.  Every sink is proportional to HLUM, so HLUM cannot go
// negative and the achlorhydric limit is approached smoothly.
// Neutralisation is BUFFER-limited, not bimolecular.  A bimolecular
// KNEUT*H*BUF term makes any remaining buffer behave as an infinitely strong
// base: it drives luminal H+ towards zero for as long as BUF > 0, so an
// untreated ZES stomach came out at postprandial pH 6.3.  Saturating in H
// caps the neutralisation flux at KNEUT*BUF, which is what a fixed titratable
// capacity can actually absorb -- and it is why a meal buffers a healthy
// stomach (BAO ~ 15 mEq/h) but is overwhelmed in ZES (BAO ~ 70 mEq/h).
double HNEUT = KNEUT * BUF * pz(HLUM) / (pz(HLUM) + KHN);
dxdt_HLUM = gBAO - KEMPTH * HLUM - HCO3 * HLUM / (HLUM + KHSAT) - HNEUT;
dxdt_BUF  = BUFMEAL * gMEAL / 1.379 - KEMPT * BUF - HNEUT;

// ===========================================================================
// 6 - antral G cell (acid-gated) and the tumour (not gated)
// ===========================================================================
double CSEC = SECRC / VSEC;
double SECN = 1.0 / (1.0 + CSEC / KSECN);              // secretin INHIBITS
double SECT = 1.0 + ESECT * CSEC / (KSECT + CSEC);     // secretin STIMULATES
double CAF  = 1.0 + ECA * pz(CAION - 1.25);
gGANTS = gBRK * SECN * (1.0 + EGRP * gMEAL) * ATRF;
dxdt_GANT = KGANT * (gGANTS - GANT);

double RANT = KSECA * GANT * CAF;
double RTUM = gESEC * (TUM + MET) * SPHEN * CAF * SECT;

dxdt_G17 = RTUM * FR17T + RANT * FR17A - KEL17 * RENF * G17 + KCONV * G34;
dxdt_G34 = RTUM * (1.0 - FR17T) + RANT * (1.0 - FR17A)
         - KEL34 * RENF * G34 - KCONV * G34;
dxdt_SECRC = -KELSEC * SECRC;

// ===========================================================================
// 7 - trophic masses (FACTOR 1) and the histamine amplifier
// ===========================================================================
double PCMSS = 1.0 + gEPCM * EXC / (EC50PCM + EXC);
double ECLSS = (1.0 + gEECL * EXC / (EC50ECL + EXC)) * (1.0 - ESSAE * gSSA);
dxdt_PCM = KPCM * (PCMSS - PCM);
dxdt_ECL = KECL * (ECLSS - ECL);
dxdt_HIS = KHIS * (ECL - HIS);

// ===========================================================================
// 8 - tumour, secretory phenotype, antitumour drug effects
// ===========================================================================
gCEVE = EVEC / VEVE;
gCSUN = SUNC / VSUN;
double OCCEVE = gCEVE / (gCEVE + IC50EVE);
double EVEEFF = EMAXEVE * OCCEVE * (1.0 - RESCL);
double SUNEFF = EMAXSUN * gCSUN / (gCSUN + IC50SUN);
double SSAEFF = ESSAT * gSSA;
double SUPP   = 1.0 - EVEEFF - SUNEFF - SSAEFF;
if (SUPP < 0.05) SUPP = 0.05;                          // never negative growth

double KGT = KG0 * GRADEF * SUPP;
gDRT = SDOSET * PRRTT;                                 // Gy/h in tumour
gDRK = SDOSEK * PRRTK;                                 // Gy/h in kidney
gALK = KALK * (CAPC / VCAP) * (TMZC / VTMZ)
     * (1.0 - EMGMT * (1.0 - MGMTM));
double KILLT = KKILL * DAMT + KRADK * gDRT;

dxdt_TUM = KGT * TUM * pz(log(TMAX / pz(TUM + 1.0e-6))) - KILLT * TUM;
dxdt_MET = KGT * KGMET * MET * pz(log(TMAX / pz(MET + 1.0e-6)))
         - KILLT * MET;
dxdt_RESCL = KRES * OCCEVE * (1.0 - RESCL);
dxdt_DAMT  = KDIN * gALK - KDOUT * DAMT;
dxdt_SPHEN = KSPH * ((1.0 - ESSAS * gSSA) - SPHEN);

// ===========================================================================
// 9 - mucosa, ulcer, oesophagus, diarrhoea
// ===========================================================================
double MUCSS = (1.0 - ENSAID * NSAID) * (1.0 - EHP * HPYL);
dxdt_MUCUS = KMUC * (MUCSS - MUCUS);
gALOAD = pz(gBAO - ATHR * MUCUS);
dxdt_ULCD = KULIN * (gALOAD / (gALOAD + KUL50)) * (1.0 - ULCD)
          - KULHEAL * ULCD * MUCUS;
dxdt_BLEED = KBLEED * ULCD * ULCD;
// Both are driven by the EXCESS acid load, not by pH alone.  Keying the
// oesophagus to "pH below 4" made every arm saturate at ESOPH ~ 0.75,
// including arms whose acid output was fully controlled, because a treated
// ZES stomach is still acidic; and keying diarrhoea to absolute acid output
// gave a healthy subject a diarrhoea index of 0.32.
dxdt_ESOPH = KESIN * (gALOAD / (gALOAD + KES50)) * (1.0 - ESOPH)
           - KESHEAL * ESOPH;
dxdt_DIARR = KDIA * (gALOAD / (gALOAD + KDIA50) - DIARR);

// ===========================================================================
// 10 - long-term consequences of a high intragastric pH
// ===========================================================================
double F12 = (1.0 / (1.0 + pow(gPH / PHB12, NB12))) / 0.981;
double FMG = FMGFL + (1.0 - FMGFL)
           * (1.0 / (1.0 + pow(gPH / PHMG, 3.0))) / 0.964;
dxdt_B12 = KB12OUT * (F12 - B12);
dxdt_MGS = KMGOUT  * (MGREF * FMG - MGS);
dxdt_BMD = -KBMD * pz(gPH - 3.5) * BMD;
dxdt_GNET = KGNET * pz(ECL - 2.0) * (1.0 + 2.0 * MENFLG) - KGNETR * GNET;

// ===========================================================================
// 11 - MEN1 calcium-PTH loop
// ===========================================================================
double PTHTGT = PTHBAS * GLANDM * pow(CASET / pz(CAION), CAHILL);
dxdt_PTH   = KPTH * (PTHTGT - PTH);
dxdt_CAION = KCAIN * PTH * VITDF - KCAOUT * CAION;

// ===========================================================================
// 12 - drug pharmacokinetics
// ===========================================================================
dxdt_PPIA = -KAPPI * PPIA;
dxdt_PPIC =  KAPPI * PPIA * FPPI - (CLPPI0 * CYPF / VPPI) * PPIC;
// the acid-activation trap: the prodrug is only converted where active pumps
// are still making acid, so the drug needs its target to be working
dxdt_PPICAN = KACTP * gCPPI * (PUMPA / (PUMPA + KMPA)) - KOUTCAN * PPICAN;

dxdt_VONA = -KAVON * VONA;
dxdt_VONC =  KAVON * VONA * FVON - (CLVON / VVON) * VONC;

dxdt_H2A  = -KAH2 * H2A;
dxdt_H2C  =  KAH2 * H2A * FH2 - (CLH2 / VH2) * H2C;
dxdt_RH2  =  KRH2 * (1.0 + TACHM * gCH2 / (gCH2 + KTACH) - RH2);

dxdt_OCTD = -KRELD * OCTD;
dxdt_OCTC =  KRELD * OCTD * 1000.0 - (CLOCT / VOCT) * OCTC;

dxdt_EVEA = -KAEVE * EVEA;
dxdt_EVEC =  KAEVE * EVEA * FEVE - (CLEVE / VEVE) * EVEC;
dxdt_SUNA = -KASUN * SUNA;
dxdt_SUNC =  KASUN * SUNA * FSUN - (CLSUN / VSUN) * SUNC;
dxdt_CAPC = -(CLCAP / VCAP) * CAPC;
dxdt_TMZC = -(CLTMZ / VTMZ) * TMZC;
dxdt_NETA = -KANET * NETA;
dxdt_NETC =  KANET * NETA - (CLNET / VNET) * NETC;

double KUPKE = KUPK * (1.0 - FPROT * AALYS);
dxdt_PRRTC = -KELPRC * PRRTC - KUPT * SSTR2D * PRRTC - KUPKE * PRRTC;
dxdt_PRRTT =  KUPT * SSTR2D * PRRTC - KOUTT * PRRTT;
dxdt_PRRTK =  KUPKE * PRRTC - KOUTK * PRRTK;
dxdt_DOSET =  gDRT;
dxdt_DOSEK =  gDRK;

// ===========================================================================
// 13 - read-out integrators (smooth indicators, so lsoda stays happy)
// ===========================================================================
dxdt_AUCBAO = gBAO;
dxdt_AUCACD = pz(gBAO - 10.0);
dxdt_TPH4   = 1.0 / (1.0 + exp(-(gPH - 4.0) / 0.08));
dxdt_TBAO10 = 1.0 / (1.0 + exp( (gBAO - 10.0) / 0.40));
dxdt_HDPCM  = gBAOP;
dxdt_HDDRV  = gBAOD;
dxdt_HDPMP  = gBAOU;

$TABLE
double BAO   = gBAO;
double PH    = gPH;
double HCONC = gHC;
double FSG   = gFSG;
double GBIO  = gGB;
double DRIVE = gDRIVE;
double ACTP  = gACTP;
double TD    = gTD;
double POOL  = gPOOL;
double PAFRC = (gPOOL > 1.0e-9) ? PUMPA / gPOOL : 0.0;
double PINH  = (gPOOL > 1.0e-9) ? (PUMPB + PUMPR) / gPOOL : 0.0;
double SSAOC = gSSA;
double HOCC  = gHOCC;
double BRK   = gBRK;
double KSYNQ = gKSYN;
double CPPI  = gCPPI;
double CVON  = gCVON;
double CH2   = gCH2;
double COCT  = gCOCT;
double CEVE  = gCEVE;
double CSUN  = gCSUN;
double ALOAD = gALOAD;
double TUMTOT = TUM + MET;
double VOLREL = (TUM0 + MET0 > 0) ? (TUM + MET) / (TUM0 + MET0) : 1.0;
double RECST  = 100.0 * (pow(pz(VOLREL), 1.0 / 3.0) - 1.0);
double BLDRSK = 1.0 - exp(-BLEED);
// headroom fractions: how much of the acid the patient actually made is
// attributable to each factor still being abnormal / still uninhibited
double HFPCM = (AUCBAO > 1.0e-6) ? 1.0 - HDPCM / AUCBAO : 0.0;
double HFDRV = (AUCBAO > 1.0e-6) ? 1.0 - HDDRV / AUCBAO : 0.0;
double HFPMP = (HDPMP  > 1.0e-6) ? 1.0 - AUCBAO / HDPMP  : 0.0;

$CAPTURE
BAO PH HCONC FSG GBIO DRIVE ACTP TD POOL PAFRC PINH SSAOC HOCC BRK KSYNQ
CPPI CVON CH2 COCT CEVE CSUN ALOAD TUMTOT VOLREL RECST BLDRSK
HFPCM HFDRV HFPMP

$OMEGA @annotated
EPCM : 0.09 : between-patient variability in the parietal trophic gain (var log)
ESEC : 0.16 : between-patient variability in tumour gastrin secretion (var log)
EECL : 0.12 : between-patient variability in the ECL trophic gain (var log)
)---"

mod <- mcode("zes", zes_code, soloc = tempdir())

## =============================================================================
##  PHENOTYPES
## =============================================================================

HEALTHY  <- list(TUM0 = 0.0, MET0 = 0.0, GLANDM = 1.0, MENFLG = 0.0)

## Sporadic ZES, moderate burden, no metastases.
ZES_SPOR <- list(TUM0 = 3.0, MET0 = 0.0, GLANDM = 1.0, MENFLG = 0.0,
                 GRADEF = 1.0, SSTR2D = 1.0)

## MEN1-ZES: several small duodenal gastrinomas plus primary
## hyperparathyroidism (GLANDM 1.8 -> mild hypercalcaemia).
ZES_MEN1 <- list(TUM0 = 1.8, MET0 = 0.0, GLANDM = 1.8, MENFLG = 1.0,
                 GRADEF = 0.85, SSTR2D = 1.1)

## Metastatic pancreatic gastrinoma, grade 2.
ZES_MET  <- list(TUM0 = 5.0, MET0 = 22.0, GLANDM = 1.0, MENFLG = 0.0,
                 GRADEF = 1.6, SSTR2D = 1.0)

## =============================================================================
##  DOSING HELPERS  (time unit = hours)
## =============================================================================

## PPI, omeprazole-equivalent milligrams PER DOSE.  `interval` in hours.
##
## `start = 7` matters and is not cosmetic.  Patients take a PPI 30-60 min
## before breakfast, so the first dose of the day lands at about 07:00 and the
## pre-breakfast fasting window (03:00-07:00, where basal acid output is
## measured) is the TROUGH.  Starting the schedule at midnight instead put
## that window 3-7 h after a once-daily dose, i.e. at peak effect, and the
## once-daily arm then beat twice-daily at identical total milligrams -- the
## exact opposite of what the model is meant to show, and an artefact of
## clock alignment rather than pharmacology.
ev_ppi <- function(dose = 20, interval = 24, days = 14, start = 7)
  ev(amt = dose, cmt = "PPIA", ii = interval,
     addl = max(0, ceiling(days * 24 / interval) - 1), time = start)

ev_ppi_iv <- function(rate_mg_per_h = 8, days = 3, start = 7)
  ev(amt = rate_mg_per_h * days * 24, cmt = "PPIC",
     rate = rate_mg_per_h, time = start)

ev_von <- function(dose = 20, interval = 24, days = 14, start = 7)
  ev(amt = dose, cmt = "VONA", ii = interval,
     addl = max(0, ceiling(days * 24 / interval) - 1), time = start)

ev_h2 <- function(dose = 40, interval = 6, days = 14, start = 7)
  ev(amt = dose, cmt = "H2A", ii = interval,
     addl = max(0, ceiling(days * 24 / interval) - 1), time = start)

## Octreotide LAR 30 mg IM every 28 days.
ev_oct <- function(dose = 30, months = 12, start = 0)
  ev(amt = dose, cmt = "OCTD", ii = 28 * 24,
     addl = max(0, months - 1), time = start)

ev_eve <- function(dose = 10, days = 730, start = 0)
  ev(amt = dose, cmt = "EVEA", ii = 24, addl = days - 1, time = start)

ev_sun <- function(dose = 37.5, days = 730, start = 0)
  ev(amt = dose, cmt = "SUNA", ii = 24, addl = days - 1, time = start)

## CAPTEM: capecitabine bd days 1-14, temozolomide od days 10-14, q28d.
ev_captem <- function(cycles = 12, capdose = 1500, tmzdose = 400, start = 0) {
  out <- list()
  for (i in seq_len(cycles)) {
    t0 <- start + (i - 1) * 28 * 24
    out[[length(out) + 1]] <-
      ev(amt = capdose, cmt = "CAPC", ii = 12, addl = 27, time = t0)
    out[[length(out) + 1]] <-
      ev(amt = tmzdose, cmt = "TMZC", ii = 24, addl = 4, time = t0 + 9 * 24)
  }
  do.call(c, out)
}

## 177Lu-DOTATATE 7.4 GBq every 8 weeks x 4 (NETTER-1 schedule).
ev_prrt <- function(gbq = 7.4, cycles = 4, start = 0)
  ev(amt = gbq, cmt = "PRRTC", ii = 8 * 7 * 24,
     addl = max(0, cycles - 1), time = start)

ev_netaz <- function(dose = 25, days = 28, start = 0)
  ev(amt = dose, cmt = "NETA", ii = 24, addl = days - 1, time = start)

## Secretin 2 U/kg IV bolus (70 kg -> 140 U).
ev_secretin <- function(units = 140, time = 0)
  ev(amt = units, cmt = "SECRC", time = time)

## =============================================================================
##  RUN-IN AND SIMULATION
## =============================================================================

.INTEGRATORS <- c("AUCBAO", "AUCACD", "TPH4", "TBAO10",
                  "HDPCM", "HDDRV", "HDPMP", "DOSET", "DOSEK", "BLEED")

.state_of <- function(out, zero_integrators = TRUE) {
  st <- as.list(out[nrow(out), names(out) %in% names(init(mod))])
  if (zero_integrators) for (nm in .INTEGRATORS) if (!is.null(st[[nm]])) st[[nm]] <- 0
  st
}

## Merge parameter lists SAFELY.  param() silently discards duplicated names,
## so a scenario override that repeats a phenotype key (e.g. SSTR2D) must
## replace it rather than be appended after it.
.pp <- function(...) {
  out <- list()
  for (x in list(...)) if (length(x)) out <- modifyList(out, as.list(x))
  out
}

## The chronic phenotype is not an initial condition.  zes_init() runs the
## untreated phenotype with no drug on board until the trophic loops settle,
## then zeroes the read-out integrators so every scenario starts its own clock.
##
## Tumour GROWTH is frozen during the run-in (GRADEF = 0).  The run-in stands
## for the years the disease took to establish at the burden the user asked
## for; letting the tumour grow through it as well would mean TUM0 no longer
## denoted the burden at presentation, and a 40-week run-in silently doubled
## the gastrinoma (this was a real defect: fasting gastrin came out at
## 2355 pg/mL for a phenotype specified to be 900).
zes_init <- function(pheno = ZES_SPOR, weeks = 40, extra = list(),
                     freeze_tumour = TRUE) {
  pp <- .pp(pheno, extra)
  if (freeze_tumour) pp$GRADEF <- 0
  m <- zero_re(param(mod, pp))
  out <- mrgsim_df(m, end = weeks * 7 * 24, delta = 6, hmax = 1)
  .state_of(out)
}

zes_sim <- function(pheno = ZES_SPOR, events = NULL, days = 28,
                    init_state = NULL, delta = 0.25, extra = list()) {
  m <- param(mod, .pp(pheno, extra))
  if (!is.null(init_state)) {
    # ICFLAG = 0 stops $MAIN from overwriting the supplied state
    m <- init(param(m, list(ICFLAG = 0)), init_state)
  }
  m <- zero_re(m)
  if (is.null(events)) mrgsim_df(m, end = days * 24, delta = delta, hmax = 0.5)
  else mrgsim_df(m, events = events, end = days * 24, delta = delta, hmax = 0.5)
}

## =============================================================================
##  DERIVED READ-OUTS
## =============================================================================

bao_ss <- function(out, win = 24) {
  s <- out[out$time > max(out$time) - win, ]; mean(s$BAO)
}
## Basal acid output as it is actually MEASURED: a fasting collection in the
## small hours, after an overnight fast and immediately before the next dose.
## This -- not the 24-h mean -- is the quantity the clinical target of
## "< 10 mEq/h" refers to, and it is where nocturnal breakthrough shows up.
bao_fasting <- function(out, win = 24, from = 3, to = 7) {
  s <- out[out$time > max(out$time) - win, ]
  tod <- s$time %% 24
  k <- tod >= from & tod < to
  if (!any(k)) return(mean(s$BAO))
  mean(s$BAO[k])
}
fsg_fasting <- function(out, win = 24, from = 3, to = 7) {
  s <- out[out$time > max(out$time) - win, ]
  tod <- s$time %% 24
  k <- tod >= from & tod < to
  if (!any(k)) return(mean(s$FSG))
  mean(s$FSG[k])
}
ph_hold <- function(out, win = 24, thr = 4) {
  s <- out[out$time > max(out$time) - win, ]; 100 * mean(s$PH > thr)
}
## Trough BAO = the worst (highest-acid) hour of the last day.  This, not the
## mean, is the quantity the clinical titration algorithm uses.
bao_trough <- function(out, win = 24) {
  s <- out[out$time > max(out$time) - win, ]; max(s$BAO)
}
## Mean acid output between 22:00 and 06:00.  `nab_hours` saturates at 8 in
## every arm that is not almost achlorhydric, so it cannot rank regimens; this
## does.
night_bao <- function(out, win = 24) {
  s <- out[out$time > max(out$time) - win, ]
  tod <- s$time %% 24
  k <- (tod >= 22) | (tod < 6)
  if (!any(k)) return(NA_real_)
  mean(s$BAO[k])
}
nab_hours <- function(out, win = 24) {
  s <- out[out$time > max(out$time) - win, ]
  tod <- s$time %% 24
  night <- (tod >= 22) | (tod < 6)
  if (!any(night)) return(NA_real_)
  ut <- sort(unique(s$time)); dt <- if (length(ut) > 1) ut[2] - ut[1] else 1
  sum(night & s$PH < 4) * dt
}

## Maximal acid output.  PENTA adds directly to CCK2R occupancy on BOTH the
## ECL cell and the parietal cell, which is what pentagastrin does; the value
## 1.5 makes it the supramaximal stimulus that 6 ug/kg is.
mao_of <- function(pheno, init_state, hours = 3, extra = list()) {
  out <- zes_sim(pheno, days = hours / 24, init_state = init_state,
                 delta = 0.05,
                 extra = .pp(extra, list(PENTA = 1.5, MEALS = 0)))
  max(out$BAO)
}

## Secretin test: peak rise in fasting serum gastrin within 30 minutes.
secretin_test <- function(pheno, init_state, units = 140) {
  out <- zes_sim(pheno, events = ev_secretin(units), days = 0.5,
                 init_state = init_state, delta = 1 / 60,
                 extra = list(MEALS = 0))
  base <- out$FSG[1]
  # report the fall as well as the rise: in a subject WITHOUT a gastrinoma the
  # informative result is inhibition, and a max()-only read-out returns zero
  # for it and looks like "no response" rather than "the opposite response".
  data.frame(basal = round(base, 0), peak = round(max(out$FSG), 0),
             nadir = round(min(out$FSG), 0),
             delta = round(max(out$FSG) - base, 0),
             fall = round(min(out$FSG) - base, 0),
             positive = (max(out$FSG) - base) > 120, row.names = NULL)
}

## =============================================================================
##  DOSE TITRATION.  The clinical algorithm is "raise the dose until acid
##  output in the hour before the next dose is below 10 mEq/h (below 5 after
##  previous gastric surgery)".  Because the model PRODUCES acid output, the
##  algorithm can be run instead of assumed, so the maintenance dose is an
##  OUTPUT of the model.
## =============================================================================

titrate_ppi <- function(pheno = ZES_SPOR, init_state = NULL,
                        daily = c(20, 40, 60, 80, 120, 160, 240),
                        interval = 12, target = 10, days = 10, extra = list()) {
  if (is.null(init_state)) init_state <- zes_init(pheno, extra = extra)
  do.call(rbind, lapply(daily, function(d) {
    per_dose <- d / (24 / interval)
    out <- zes_sim(pheno, ev_ppi(per_dose, interval, days), days,
                   init_state, 0.25, extra)
    data.frame(total_daily_mg = d, interval_h = interval,
               BAO_24h_mean = round(bao_ss(out), 1),
               BAO_fasting = round(bao_fasting(out), 1),
               BAO_peak = round(bao_trough(out), 1),
               pH4_pct = round(ph_hold(out), 0),
               NAB_h = round(nab_hours(out), 1),
               pump_inhib = round(mean(out$PINH[out$time > (days - 1) * 24]), 3),
               at_target = bao_fasting(out) < target, row.names = NULL)
  }))
}

## =============================================================================
##  SCENARIOS
## =============================================================================

zes_scenarios <- function() list(
  ## ---- reference states ---------------------------------------------------
  "01 Healthy stomach, no drug" =
    list(ph = HEALTHY, ev = NULL, days = 3, init = "HEALTHY"),
  "02 Untreated sporadic ZES" =
    list(ph = ZES_SPOR, ev = NULL, days = 3, init = "ZES_SPOR"),
  "03 Untreated MEN1-ZES, pHPT active" =
    list(ph = ZES_MEN1, ev = NULL, days = 3, init = "ZES_MEN1"),
  "04 Untreated metastatic gastrinoma" =
    list(ph = ZES_MET, ev = NULL, days = 3, init = "ZES_MET"),

  ## ---- interval versus dose, at IDENTICAL total daily milligrams ---------
  "05 ZES omeprazole 60 mg ONCE daily" =
    list(ph = ZES_SPOR, ev = ev_ppi(60, 24, 10), days = 10, init = "ZES_SPOR"),
  "06 ZES omeprazole 30 mg TWICE daily" =
    list(ph = ZES_SPOR, ev = ev_ppi(30, 12, 10), days = 10, init = "ZES_SPOR"),
  "07 ZES omeprazole 20 mg THREE times daily" =
    list(ph = ZES_SPOR, ev = ev_ppi(20, 8, 10), days = 10, init = "ZES_SPOR"),

  ## ---- the same regimen with and without a gastrinoma --------------------
  "08 Healthy omeprazole 20 mg once daily" =
    list(ph = HEALTHY, ev = ev_ppi(20, 24, 10), days = 10, init = "HEALTHY"),
  "09 ZES omeprazole 20 mg once daily (reflux dose)" =
    list(ph = ZES_SPOR, ev = ev_ppi(20, 24, 10), days = 10, init = "ZES_SPOR"),

  ## ---- CYP2C19 at identical dose ----------------------------------------
  "10 ZES omeprazole 30 mg bd, CYP2C19 UM" =
    list(ph = ZES_SPOR, ev = ev_ppi(30, 12, 10), days = 10, init = "ZES_SPOR",
         extra = list(CYPF = 1.9)),
  "11 ZES omeprazole 30 mg bd, CYP2C19 PM" =
    list(ph = ZES_SPOR, ev = ev_ppi(30, 12, 10), days = 10, init = "ZES_SPOR",
         extra = list(CYPF = 0.28)),

  ## ---- the three acid-suppressing classes at licensed doses -------------
  "12 ZES vonoprazan 20 mg once daily" =
    list(ph = ZES_SPOR, ev = ev_von(20, 24, 10), days = 10, init = "ZES_SPOR"),
  "13 ZES vonoprazan 40 mg once daily" =
    list(ph = ZES_SPOR, ev = ev_von(40, 24, 10), days = 10, init = "ZES_SPOR"),
  "12b Healthy vonoprazan 20 mg once daily" =
    list(ph = HEALTHY, ev = ev_von(20, 24, 10), days = 10, init = "HEALTHY"),
  "14 ZES famotidine 40 mg four times daily" =
    list(ph = ZES_SPOR, ev = ev_h2(40, 6, 10), days = 10, init = "ZES_SPOR"),
  "15 ZES famotidine 120 mg four times daily" =
    list(ph = ZES_SPOR, ev = ev_h2(120, 6, 10), days = 10, init = "ZES_SPOR"),

  ## ---- attacking factors 1 and 3 ----------------------------------------
  "16 ZES octreotide LAR 30 mg alone" =
    list(ph = ZES_SPOR, ev = ev_oct(30, 6), days = 180, init = "ZES_SPOR",
         delta = 6),
  "17 MEN1-ZES parathyroidectomy at month 2" =
    list(ph = ZES_MEN1, ev = NULL, days = 180, init = "ZES_MEN1", delta = 6,
         switch = list(at = 60 * 24, to = list(GLANDM = 1.0))),
  "17b MEN1-ZES no parathyroidectomy (matched control)" =
    list(ph = ZES_MEN1, ev = NULL, days = 180, init = "ZES_MEN1", delta = 6),
  "18 ZES curative duodenotomy at month 2" =
    list(ph = ZES_SPOR, ev = NULL, days = 180, init = "ZES_SPOR", delta = 6,
         switch = list(at = 60 * 24, to = list(TUMCUT = 0.02))),
  "19 ZES netazepide 25 mg daily (CCK2R blockade)" =
    list(ph = ZES_SPOR, ev = ev_netaz(25, 28), days = 28, init = "ZES_SPOR",
         delta = 0.5),

  ## ---- factor 2 plus factors 1 and 3 ------------------------------------
  "20 ZES omeprazole 30 mg bd + octreotide LAR" =
    list(ph = ZES_SPOR, ev = c(ev_ppi(30, 12, 180), ev_oct(30, 6)),
         days = 180, init = "ZES_SPOR", delta = 6),

  ## ---- withdrawal and rebound -------------------------------------------
  "21 Healthy 8 weeks omeprazole 40 mg then stop" =
    list(ph = HEALTHY, ev = ev_ppi(40, 24, 56), days = 84, init = "HEALTHY",
         delta = 2),
  "22 ZES 8 weeks omeprazole 30 mg bd then stop" =
    list(ph = ZES_SPOR, ev = ev_ppi(30, 12, 56), days = 84, init = "ZES_SPOR",
         delta = 2),

  ## ---- long-horizon control and its harms --------------------------------
  "23 ZES 5 years omeprazole 30 mg bd" =
    list(ph = ZES_SPOR, ev = ev_ppi(30, 12, 1825), days = 1825,
         init = "ZES_SPOR", delta = 24),
  "24 MEN1-ZES 5 years omeprazole 30 mg bd" =
    list(ph = ZES_MEN1, ev = ev_ppi(30, 12, 1825), days = 1825,
         init = "ZES_MEN1", delta = 24),
  "25 ZES 5 years untreated (historical natural history)" =
    list(ph = ZES_SPOR, ev = NULL, days = 1825, init = "ZES_SPOR", delta = 24),

  ## ---- antitumour arms, all from the same metastatic baseline ------------
  "26 Metastatic gastrinoma, PPI only" =
    list(ph = ZES_MET, ev = ev_ppi(40, 12, 730), days = 730, init = "ZES_MET",
         delta = 24),
  "27 Metastatic + octreotide LAR 30 mg" =
    list(ph = ZES_MET, ev = c(ev_ppi(40, 12, 730), ev_oct(30, 26)),
         days = 730, init = "ZES_MET", delta = 24),
  "28 Metastatic + everolimus 10 mg" =
    list(ph = ZES_MET, ev = c(ev_ppi(40, 12, 730), ev_eve(10, 730)),
         days = 730, init = "ZES_MET", delta = 24),
  "29 Metastatic + sunitinib 37.5 mg" =
    list(ph = ZES_MET, ev = c(ev_ppi(40, 12, 730), ev_sun(37.5, 730)),
         days = 730, init = "ZES_MET", delta = 24),
  "30 Metastatic + CAPTEM x 12" =
    list(ph = ZES_MET, ev = c(ev_ppi(40, 12, 730), ev_captem(12)),
         days = 730, init = "ZES_MET", delta = 24),
  "31 Metastatic + 177Lu-DOTATATE x 4" =
    list(ph = ZES_MET, ev = c(ev_ppi(40, 12, 730), ev_prrt(7.4, 4)),
         days = 730, init = "ZES_MET", delta = 24),
  "32 Metastatic + PRRT, SSTR2-low tumour" =
    list(ph = ZES_MET, ev = c(ev_ppi(40, 12, 730), ev_prrt(7.4, 4)),
         days = 730, init = "ZES_MET", delta = 24,
         extra = list(SSTR2D = 0.25)),
  "33 Metastatic + PRRT without amino-acid protection" =
    list(ph = ZES_MET, ev = c(ev_ppi(40, 12, 730), ev_prrt(7.4, 4)),
         days = 730, init = "ZES_MET", delta = 24, extra = list(AALYS = 0)),

  ## ---- co-morbid exposures ----------------------------------------------
  "34 ZES on PPI plus NSAID" =
    list(ph = ZES_SPOR, ev = ev_ppi(30, 12, 90), days = 90, init = "ZES_SPOR",
         delta = 3, extra = list(NSAID = 1.0)),
  "35 ZES on PPI with renal failure (gastrin CL halved)" =
    list(ph = ZES_SPOR, ev = ev_ppi(30, 12, 90), days = 90, init = "ZES_SPOR",
         delta = 3, extra = list(RENF = 0.5))
)

## A mid-run parameter switch (parathyroidectomy, resection) is done by
## splitting the simulation.  TUMCUT is not a parameter -- it is a state edit.
run_scenario <- function(sc, inits) {
  ph <- sc$ph
  init_state <- inits[[sc$init]]
  delta <- if (is.null(sc$delta)) 0.25 else sc$delta
  extra <- if (is.null(sc$extra)) list() else sc$extra
  if (is.null(sc$switch))
    return(zes_sim(ph, sc$ev, sc$days, init_state, delta, extra))
  t1 <- sc$switch$at / 24
  a  <- zes_sim(ph, sc$ev, t1, init_state, delta, extra)
  st <- .state_of(a, zero_integrators = FALSE)
  to <- sc$switch$to
  if (!is.null(to$TUMCUT)) { st$TUM <- st$TUM * to$TUMCUT; to$TUMCUT <- NULL }
  b <- zes_sim(ph, sc$ev, sc$days - t1, st, delta, .pp(extra, to))
  b$time <- b$time + sc$switch$at
  rbind(a, b)
}

zes_inits <- function() list(
  HEALTHY  = zes_init(HEALTHY,  weeks = 12),
  ZES_SPOR = zes_init(ZES_SPOR, weeks = 40),
  ZES_MEN1 = zes_init(ZES_MEN1, weeks = 40),
  ZES_MET  = zes_init(ZES_MET,  weeks = 40))

zes_run_all <- function(scenarios = zes_scenarios(), inits = zes_inits()) {
  out <- lapply(names(scenarios), function(nm) {
    x <- run_scenario(scenarios[[nm]], inits); x$scenario <- nm; x
  })
  do.call(rbind, out)
}

summarise_scn <- function(all) {
  do.call(rbind, lapply(split(all, all$scenario), function(x) {
    last <- x[x$time > max(x$time) - 24, ]
    data.frame(scenario = x$scenario[1],
      BAO_mean = round(mean(last$BAO), 1), BAO_fast = round(bao_fasting(x), 1),
      BAO_peak = round(max(last$BAO), 1),
      pH_mean = round(mean(last$PH), 2),
      pH4_pct = round(100 * mean(last$PH > 4), 0),
      FSG = round(mean(last$FSG), 0), PCM = round(mean(last$PCM), 2),
      ECL = round(mean(last$ECL), 2), pumpinh = round(mean(last$PINH), 3),
      ULCD = round(mean(last$ULCD), 3), DIARR = round(mean(last$DIARR), 3),
      TUM_cm3 = round(mean(last$TUMTOT), 1), RECIST = round(mean(last$RECST), 1),
      HF_PCM = round(mean(last$HFPCM), 3), HF_DRIVE = round(mean(last$HFDRV), 3),
      HF_PUMP = round(mean(last$HFPMP), 3), row.names = NULL)
  }))
}

## =============================================================================
##  CALIBRATION
##  ---------------------------------------------------------------------------
##  Nine parameters are fitted to nine published anchors.  Every other
##  parameter is fixed from the literature (see zes_references.md).
##
##    #   anchor                                            target   parameter
##    ------------------------------------------------------------------------
##    A1  healthy fasting basal acid output              3.0 mEq/h   BAOCAP
##    A2  healthy MAO / BAO ratio (pentagastrin)              7.67   KDRV
##    A3  untreated sporadic ZES fasting acid output    36.0 mEq/h   EMAXPCM
##    A4  untreated sporadic ZES fasting serum gastrin  900 pg/mL    KSECT0
##    A5  healthy fasting serum gastrin                  30 pg/mL    KSECA
##    A6  healthy parietal cell mass (no drift)               1.00   GBNORM
##    A7  healthy, omeprazole 20 mg od d6: acid fall        66 %      KBIND
##    A8  healthy, untreated: 24-h pH > 4 holding time      20 %      BUFMEAL
##    A9  healthy, omeprazole 20 mg od: pH > 4 holding      45 %      KEMPT
##
##  A6 is a self-consistency condition rather than an observation.  GBNORM
##  defines the gastrin level at which the trophic drive TD equals 1, so if it
##  is not the model's own healthy steady-state gastrin the trophic loops drift
##  and a stomach with no tumour in it settles at a parietal mass of 1.2.  It
##  is solved as a fixed point (GBNORM <- healthy 24-h mean GBIO).
##
##  WHAT IS *NOT* FITTED, AND WHY THAT MATTERS
##
##  A tenth anchor was tried and abandoned: the acid output a ZES patient
##  actually achieves on omeprazole 30 mg twice daily (about 6 mEq/h in the
##  pantoprazole and lansoprazole maintenance series).  The model UNDER-shoots
##  it -- it over-suppresses ZES -- and the obvious repair was to fit GSYN, the
##  gastrin up-regulation of pump synthesis, so that the ZES pump pool refilled
##  fast enough to resist.  That pairing destroys the calibration: GSYN's
##  leverage on the UNTREATED ZES acid output (A3) is larger than its leverage
##  on the treated one, so GSYN and EMAXPCM fight over A3, the iteration
##  oscillates and then diverges (untreated ZES acid output ran to 144 mEq/h
##  and EMAXPCM collapsed to 0.07 before the caps caught it).
##
##  GSYN is therefore left at its literature value and the over-suppression is
##  reported as a held-out FAILURE in README.md rather than absorbed into a
##  parameter.  The mechanism responsible is identified there.
##
##  DELIBERATELY HELD OUT: the BAO/MAO ratio in ZES, the secretin-test
##  increment and its PPI false positive, the CYP2C19 gradient, vonoprazan
##  versus omeprazole, H2-antagonist escape and tachyphylaxis, withdrawal
##  rebound, the effect of parathyroidectomy, the interval-versus-dose contrast
##  at equal milligrams, the treated ZES acid output, and all antitumour arms.
##
##  METHOD.  Coordinate-wise fixed-point iteration rather than Nelder-Mead.
##  Each anchor has one parameter with dominant leverage on it, so a damped
##  multiplicative update per coordinate costs six simulations per pass, where
##  a 9-dimensional Nelder-Mead on the same objective needs several hundred
##  evaluations of the same six simulations (tried; not converged after forty
##  minutes).  The coordinates are NOT independent -- BAOCAP, EMAXPCM and
##  KSECA all move acid output, and acid output moves the D-cell brake, which
##  moves gastrin, which moves them back -- so every update is damped well
##  below 1 and the routine returns the BEST pass by objective rather than the
##  last one.
## =============================================================================

zes_anchors <- c(BAOh = 3.0, ratio = 7.667, BAOz = 36.0, FSGz = 900,
                 FSGh = 30.0, PCMh = 1.0, INH = 66.0, PH4h = 20.0,
                 PH4hp = 45.0)

## The nine anchored quantities plus three diagnostics, from six simulations.
zes_predict_anchors <- function(p = list()) {
  hi <- zes_init(HEALTHY,  weeks = 14, extra = p)
  zi <- zes_init(ZES_SPOR, weeks = 32, extra = p)
  h  <- zes_sim(HEALTHY,  NULL, 3, hi, 0.25, p)
  z  <- zes_sim(ZES_SPOR, NULL, 3, zi, 0.25, p)
  hp <- zes_sim(HEALTHY,  ev_ppi(20, 24, 6),  6,  hi, 0.25, p)
  zp <- zes_sim(ZES_SPOR, ev_ppi(30, 12, 10), 10, zi, 0.25, p)
  c(BAOh  = bao_fasting(h),
    ratio = mao_of(HEALTHY, hi, extra = p) / bao_fasting(h),
    BAOz  = bao_fasting(z),
    FSGz  = fsg_fasting(z),
    FSGh  = fsg_fasting(h),
    PCMh  = mean(tail(h$PCM, 96)),
    INH   = 100 * (1 - bao_fasting(hp) / bao_fasting(h)),
    PH4h  = ph_hold(h),
    PH4hp = ph_hold(hp),
    GB24  = mean(tail(h$GBIO, 96)),
    BAOzt = bao_fasting(zp),            # held out
    BAOMAOz = bao_fasting(z) / mao_of(ZES_SPOR, zi, extra = p))  # held out
}

zes_objective <- function(p = list(), verbose = FALSE) {
  val <- tryCatch(zes_predict_anchors(p),
                  error = function(e) rep(NA_real_, 12))
  a <- val[names(zes_anchors)]
  if (any(!is.finite(a))) return(1e6)
  if (verbose) print(rbind(model = round(a, 2), target = zes_anchors))
  sum(((a - zes_anchors) / zes_anchors)^2)
}

ZES_FITTED_START <- list(BAOCAP = 78.879, KDRV = 45.447, EMAXPCM = 2.0674,
                         KSECT0 = 562.45, KBIND = 40.0, KSECA = 901.01,
                         GBNORM = 27.121, BUFMEAL = 107.1, KEMPT = 2.6218)

zes_calibrate <- function(start = ZES_FITTED_START, passes = 24,
                          trace = TRUE) {
  P <- start; TG <- zes_anchors
  best <- list(obj = Inf, par = P, pred = NULL, pass = NA)
  for (it in seq_len(passes)) {
    r <- zes_predict_anchors(P)
    a <- r[names(TG)]
    obj <- if (any(!is.finite(a))) 1e6 else sum(((a - TG) / TG)^2)
    if (is.finite(obj) && obj < best$obj)
      best <- list(obj = obj, par = P, pred = r, pass = it)
    if (trace) cat(sprintf(
      "pass %2d  obj %7.4f | BAOh %5.2f  ratio %5.2f  BAOz %6.1f  FSGz %4.0f  FSGh %5.1f  PCMh %5.3f  INH %5.1f  pH4 %4.1f/%4.1f | held-out ZESrx %5.2f\n",
      it, obj, r["BAOh"], r["ratio"], r["BAOz"], r["FSGz"], r["FSGh"],
      r["PCMh"], r["INH"], r["PH4h"], r["PH4hp"], r["BAOzt"]))
    P$BAOCAP  <- P$BAOCAP  * (TG["BAOh"]  / r["BAOh"])^0.35
    P$KDRV    <- P$KDRV    * (TG["ratio"] / r["ratio"])^0.30
    P$EMAXPCM <- min(6, P$EMAXPCM * (TG["BAOz"] / r["BAOz"])^0.30)
    P$KSECT0  <- P$KSECT0  * (TG["FSGz"]  / r["FSGz"])^0.70
    P$KSECA   <- P$KSECA   * (TG["FSGh"]  / r["FSGh"])^0.35
    P$GBNORM  <- P$GBNORM + 0.5 * (r["GB24"] - P$GBNORM)
    P$KBIND   <- min(150, max(0.05,
                     P$KBIND * ((100 - r["INH"]) / (100 - TG["INH"]))^0.50))
    P$BUFMEAL <- min(110, max(5,
                     P$BUFMEAL * (TG["PH4h"] / max(r["PH4h"], 1))^0.22))
    P$KEMPT   <- min(4, max(0.1,
                     P$KEMPT * (r["PH4hp"] / TG["PH4hp"])^0.22))
    P <- lapply(P, as.numeric)
  }
  list(par = lapply(best$par, function(x) signif(x, 5)),
       predicted = round(best$pred, 2), target = TG,
       objective = round(best$obj, 5), best_pass = best$pass)
}

## =============================================================================
##  NUMERICAL VALIDATION SUITE  (A0 - A13)
##  None of these results is coded anywhere in the equations.
## =============================================================================

.PHENO_OF <- function(nm) switch(nm, HEALTHY = HEALTHY, ZES_SPOR = ZES_SPOR,
                                 ZES_MEN1 = ZES_MEN1, ZES_MET = ZES_MET)

## A0 - reference states, including the BAO/MAO ratio (>0.6 is diagnostic).
A0_reference <- function(inits) {
  do.call(rbind, lapply(names(inits), function(nm) {
    st <- inits[[nm]]; ph <- .PHENO_OF(nm)
    o <- zes_sim(ph, NULL, 2, st, 0.25)
    mao <- mao_of(ph, st)
    data.frame(state = nm, BAO = round(bao_fasting(o), 1), MAO = round(mao, 1),
               BAO_MAO = round(bao_fasting(o) / mao, 2),
               pH = round(mean(o$PH[o$time > 24]), 2),
               FSG = round(mean(o$FSG[o$time > 24]), 0),
               PCM = round(mean(o$PCM), 2), ECL = round(mean(o$ECL), 2),
               CAION = round(mean(o$CAION), 2),
               DIARR = round(mean(o$DIARR), 2), row.names = NULL)
  }))
}

## A1 - the three-factor decomposition, MEASURED.
A1_factors <- function(all) {
  labs <- c("02 Untreated sporadic ZES", "06 ZES omeprazole 30 mg TWICE daily",
            "20 ZES omeprazole 30 mg bd + octreotide LAR")
  s <- summarise_scn(all)
  s[s$scenario %in% labs,
    c("scenario", "BAO_mean", "HF_PCM", "HF_DRIVE", "HF_PUMP",
      "PCM", "ECL", "pumpinh")]
}

## A2 - interval versus total daily dose at identical milligrams.
A2_interval <- function(all) {
  labs <- c("05 ZES omeprazole 60 mg ONCE daily",
            "06 ZES omeprazole 30 mg TWICE daily",
            "07 ZES omeprazole 20 mg THREE times daily")
  do.call(rbind, lapply(labs, function(l) {
    x <- all[all$scenario == l, ]
    data.frame(scenario = l, total_mg_per_day = 60,
               BAO_mean = round(bao_ss(x), 1),
               BAO_fasting = round(bao_fasting(x), 1),
               BAO_peak = round(bao_trough(x), 1),
               night_BAO = round(night_bao(x), 1),
               pH4_pct = round(ph_hold(x), 0), NAB_h = round(nab_hours(x), 1),
               row.names = NULL)
  }))
}

## A3 - the same dose with and without a gastrinoma.
A3_dose_gap <- function(all) {
  labs <- c("08 Healthy omeprazole 20 mg once daily",
            "09 ZES omeprazole 20 mg once daily (reflux dose)")
  do.call(rbind, lapply(labs, function(l) {
    x <- all[all$scenario == l, ]
    last <- x[x$time > max(x$time) - 24, ]
    data.frame(scenario = l, BAO_mean = round(bao_ss(x), 1),
               BAO_fasting = round(bao_fasting(x), 1),
               pump_inhibited = round(mean(last$PINH), 3),
               pH4_pct = round(ph_hold(x), 0), row.names = NULL)
  }))
}

## A4 - titration curves: the maintenance dose as an OUTPUT.
A4_titration <- function(inits) {
  rbind(cbind(pheno = "sporadic", titrate_ppi(ZES_SPOR, inits$ZES_SPOR)),
        cbind(pheno = "MEN1",     titrate_ppi(ZES_MEN1, inits$ZES_MEN1)),
        cbind(pheno = "healthy",  titrate_ppi(HEALTHY,  inits$HEALTHY,
                                              daily = c(20, 40))))
}

## A5 - CYP2C19 phenotype at identical dose.
A5_cyp <- function(all) {
  labs <- c("10 ZES omeprazole 30 mg bd, CYP2C19 UM",
            "06 ZES omeprazole 30 mg TWICE daily",
            "11 ZES omeprazole 30 mg bd, CYP2C19 PM")
  do.call(rbind, lapply(labs, function(l) {
    x <- all[all$scenario == l, ]
    data.frame(scenario = l, BAO_mean = round(bao_ss(x), 1),
               BAO_fasting = round(bao_fasting(x), 1),
               pH4_pct = round(ph_hold(x), 0), row.names = NULL)
  }))
}

## A6 - drug-class comparison at licensed doses, plus H2 tachyphylaxis.
A6_classes <- function(all) {
  labs <- c("08 Healthy omeprazole 20 mg once daily",
            "12b Healthy vonoprazan 20 mg once daily",
            "06 ZES omeprazole 30 mg TWICE daily",
            "12 ZES vonoprazan 20 mg once daily",
            "13 ZES vonoprazan 40 mg once daily",
            "14 ZES famotidine 40 mg four times daily",
            "15 ZES famotidine 120 mg four times daily")
  do.call(rbind, lapply(labs, function(l) {
    x <- all[all$scenario == l, ]; last <- x[x$time > max(x$time) - 24, ]
    data.frame(scenario = l, BAO_mean = round(bao_ss(x), 1),
               BAO_fasting = round(bao_fasting(x), 1),
               night_BAO = round(night_bao(x), 1),
               pH4_pct = round(ph_hold(x), 0),
               H2R_density = round(mean(last$RH2), 2), row.names = NULL)
  }))
}

## A7 - the secretin sign inversion AND its false positive.
A7_secretin <- function(inits) {
  o <- zes_sim(HEALTHY, ev_ppi(40, 24, 28), 28, inits$HEALTHY, 1)
  ppi_state <- .state_of(o, zero_integrators = FALSE)
  rbind(cbind(subject = "healthy",        secretin_test(HEALTHY,  inits$HEALTHY)),
        cbind(subject = "healthy on PPI", secretin_test(HEALTHY,  ppi_state)),
        cbind(subject = "sporadic ZES",   secretin_test(ZES_SPOR, inits$ZES_SPOR)),
        cbind(subject = "MEN1 ZES",       secretin_test(ZES_MEN1, inits$ZES_MEN1)))
}

## A8 - the MEN1 calcium loop.  Reported as a DIFFERENCE IN DIFFERENCES
## against an identical patient who does not have the operation, because the
## gastrinoma keeps growing through the six months either way: the first
## version of this test compared the operated patient with its own pre-operative
## self and reported that parathyroidectomy RAISED gastrin, which was entirely
## the four months of tumour growth in between.
A8_calcium <- function(all) {
  g <- function(l, t0, t1) {
    x <- all[all$scenario == l, ]
    x[x$time > t0 * 24 & x$time <= t1 * 24, ]
  }
  arms <- c(PTX = "17 MEN1-ZES parathyroidectomy at month 2",
            control = "17b MEN1-ZES no parathyroidectomy (matched control)")
  out <- do.call(rbind, lapply(names(arms), function(nm) {
    pre <- g(arms[[nm]], 55, 60); post <- g(arms[[nm]], 170, 180)
    data.frame(arm = nm,
               CAION_pre = round(mean(pre$CAION), 2),
               CAION_post = round(mean(post$CAION), 2),
               PTH_post = round(mean(post$PTH), 1),
               FSG_pre = round(mean(pre$FSG), 0),
               FSG_post = round(mean(post$FSG), 0),
               BAO_pre = round(mean(pre$BAO), 1),
               BAO_post = round(mean(post$BAO), 1),
               PCM_post = round(mean(post$PCM), 2), row.names = NULL)
  }))
  dif <- data.frame(arm = "PTX effect (diff-in-diff)",
    CAION_pre = NA, CAION_post = round(out$CAION_post[1] - out$CAION_post[2], 2),
    PTH_post = round(out$PTH_post[1] - out$PTH_post[2], 1), FSG_pre = NA,
    FSG_post = round(out$FSG_post[1] - out$FSG_post[2], 0), BAO_pre = NA,
    BAO_post = round(out$BAO_post[1] - out$BAO_post[2], 1),
    PCM_post = round(out$PCM_post[1] - out$PCM_post[2], 2))
  rbind(out, dif)
}

## A9 - withdrawal rebound in a healthy stomach and in ZES.
A9_rebound <- function(all) {
  do.call(rbind, lapply(c("21 Healthy 8 weeks omeprazole 40 mg then stop",
                          "22 ZES 8 weeks omeprazole 30 mg bd then stop"),
    function(l) {
      x <- all[all$scenario == l, ]
      fw <- function(d) {                       # fasting window of one day
        tod <- d$time %% 24
        k <- tod >= 3 & tod < 7
        if (any(k)) mean(d$BAO[k]) else mean(d$BAO)
      }
      pre  <- x[x$time < 24, ]
      post <- x[x$time > 56 * 24 & x$time < 63 * 24, ]
      late <- x[x$time > 80 * 24, ]
      on   <- x[x$time > 50 * 24 & x$time < 56 * 24, ]
      # compare like with like: fasting acid output before, on treatment, in
      # the week after stopping, and four weeks later
      data.frame(scenario = l,
                 fasting_before = round(fw(pre), 1),
                 fasting_on_PPI = round(fw(on), 1),
                 fasting_after_stop = round(fw(post), 1),
                 rebound_ratio = round(fw(post) / fw(pre), 2),
                 fasting_4wk_after = round(fw(late), 1),
                 pump_pool_at_stop = round(max(x$POOL[x$time < 56 * 24]), 2),
                 row.names = NULL)
    }))
}

## A10 - antitumour arms: time to a +20% RECIST change from one baseline.
A10_tumour <- function(all) {
  labs <- c("26 Metastatic gastrinoma, PPI only",
            "27 Metastatic + octreotide LAR 30 mg",
            "28 Metastatic + everolimus 10 mg",
            "29 Metastatic + sunitinib 37.5 mg",
            "30 Metastatic + CAPTEM x 12",
            "31 Metastatic + 177Lu-DOTATATE x 4",
            "32 Metastatic + PRRT, SSTR2-low tumour",
            "33 Metastatic + PRRT without amino-acid protection")
  do.call(rbind, lapply(labs, function(l) {
    x <- all[all$scenario == l, ]
    prog <- which(x$RECST >= 20)
    ttp <- if (length(prog)) x$time[prog[1]] / (30.44 * 24) else NA_real_
    data.frame(scenario = l, best_RECIST = round(min(x$RECST), 1),
               TTP_months = round(ttp, 1),
               vol_24mo = round(x$TUMTOT[nrow(x)], 1),
               FSG_24mo = round(x$FSG[nrow(x)], 0),
               BAO_24mo = round(x$BAO[nrow(x)], 1),
               tumour_Gy = round(max(x$DOSET), 1),
               renal_Gy = round(max(x$DOSEK), 1), row.names = NULL)
  }))
}

## A11 - long-term harms of decades of profound acid suppression.
A11_harms <- function(all) {
  labs <- c("23 ZES 5 years omeprazole 30 mg bd",
            "24 MEN1-ZES 5 years omeprazole 30 mg bd",
            "25 ZES 5 years untreated (historical natural history)")
  do.call(rbind, lapply(labs, function(l) {
    x <- all[all$scenario == l, ]; e <- x[nrow(x), ]
    data.frame(scenario = l, pH_5yr = round(e$PH, 2),
               BAO_5yr = round(e$BAO, 1), B12 = round(e$B12, 2),
               Mg = round(e$MGS, 2), BMD = round(e$BMD, 3),
               ECL = round(e$ECL, 2), gastric_NET_cm3 = round(e$GNET, 2),
               ULCD = round(e$ULCD, 3), ESOPH = round(e$ESOPH, 3),
               bleed_risk_5yr = round(e$BLDRSK, 3), row.names = NULL)
  }))
}

## A12 - is the somatostatin-analogue benefit additive with a PPI?
A12_additivity <- function(all) {
  g <- function(l) bao_fasting(all[all$scenario == l, ])
  base <- g("02 Untreated sporadic ZES")
  ppi  <- g("06 ZES omeprazole 30 mg TWICE daily")
  ssa  <- g("16 ZES octreotide LAR 30 mg alone")
  both <- g("20 ZES omeprazole 30 mg bd + octreotide LAR")
  expct <- base * (ppi / base) * (ssa / base)
  data.frame(baseline_BAO = round(base, 1), PPI_alone = round(ppi, 1),
             SSA_alone = round(ssa, 1), combination = round(both, 1),
             expected_if_independent = round(expct, 1),
             observed_minus_expected = round(both - expct, 2),
             row.names = NULL)
}

## A13 - virtual population: what fraction reaches the acid target on the
## standard 30 mg bd regimen?
## The three random effects act on TROPHIC gains, whose time constants are
## weeks.  Injecting them after a deterministic run-in and then simulating ten
## days of treatment produced an interquartile range of 3.6-4.2 mEq/h and a
## 100% target-attainment rate -- i.e. no variability at all, because every
## subject started from the same settled parietal and ECL mass.  The run-in
## and the treatment period therefore have to be ONE simulation per subject.
A13_population <- function(n = 120, pheno = ZES_SPOR, dose = 30,
                           interval = 12, seed = 20260730, runin_wk = 40) {
  set.seed(seed)
  t0 <- runin_wk * 7 * 24
  pp <- .pp(pheno, list(GRADEF = 0))          # burden fixed, as in zes_init()
  out <- mrgsim_df(param(mod, pp),
                   events = ev_ppi(dose, interval, 10, start = t0 + 7),
                   end = t0 + 240, delta = 0.5, nid = n, hmax = 0.5)
  out <- out[out$time >= t0, ]
  out$time <- out$time - t0
  per <- do.call(rbind, lapply(split(out, out$ID), function(x) {
    last <- x[x$time > 216, ]
    data.frame(ID = x$ID[1], BAO = mean(last$BAO), fasting = bao_fasting(x),
               pH4 = 100 * mean(last$PH > 4), FSG = mean(last$FSG))
  }))
  data.frame(n = n, dose_mg_per_dose = dose, interval_h = interval,
             BAO_median = round(median(per$fasting), 1),
             BAO_IQR = paste(round(quantile(per$fasting, c(.25, .75)), 1),
                             collapse = "-"),
             at_target_pct = round(100 * mean(per$fasting < 10), 0),
             pH4_median = round(median(per$pH4), 0),
             FSG_median = round(median(per$FSG), 0), row.names = NULL)
}

zes_validate <- function() {
  inits <- zes_inits()
  all   <- zes_run_all(inits = inits)
  list(A0 = A0_reference(inits), A1 = A1_factors(all), A2 = A2_interval(all),
       A3 = A3_dose_gap(all), A4 = A4_titration(inits), A5 = A5_cyp(all),
       A6 = A6_classes(all), A7 = A7_secretin(inits), A8 = A8_calcium(all),
       A9 = A9_rebound(all), A10 = A10_tumour(all), A11 = A11_harms(all),
       A12 = A12_additivity(all), A13 = A13_population(),
       summary = summarise_scn(all))
}

if (identical(environment(), globalenv()) &&
    !is.null(getOption("zes.run.scenarios"))) {
  res <- zes_validate()
  for (nm in names(res)) { cat("\n==== ", nm, " ====\n"); print(res[[nm]]) }
}
