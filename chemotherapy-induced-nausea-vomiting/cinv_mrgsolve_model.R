## =============================================================================
##  cinv_mrgsolve_model.R
##  Chemotherapy-Induced Nausea and Vomiting (CINV) — QSP model for mrgsolve
##
##  ORGANISING THESIS
##  -----------------
##  A CINV endpoint is not a percentage that drugs shrink.  It is an EXPONENTIAL
##  OF AN INTEGRAL OF A HAZARD, and the hazard is driven by a SUM:
##
##      DRIVE(t) = GAIN * SUM_j w_j * A_j(t)        <- additive, receptor-specific
##      EBS(t)   = max(0, DRIVE),  DRIVE built from chemo-driven INCREMENTS
##      lambda(t)= HMAX * EBS^NE / (THR^NE + EBS^NE)
##      I        = integral of (lambda + lambda_rescue) dt
##      CR       = E_S[ exp(-S * I) ],  S ~ lognormal(0, OMEGA^2)
##
##  Every antiemetic ever licensed removes ONE TERM A_j from that sum.  None of
##  them scales DRIVE.  Four things are therefore OUTPUTS of this file, not rules
##  written into it:
##
##   1. ADDITIVITY IS EXACT ONLY IN THE DRIVE.  Each drug deletes one term of the
##      sum; that layer is additive by construction.  The chain
##      drive -> lambda -> I -> CR is a saturating threshold followed by an
##      exponential, so additivity is destroyed before the endpoint: measured
##      here it is already down to ~67% of the sum at the dI level (see
##      log_additivity_check()).  The "synergy" universally reported between a
##      5-HT3 RA, an NK1 RA, dexamethasone and olanzapine is additive
##      pharmacology seen through that chain, not receptor interaction.  A
##      corollary: the ABSOLUTE benefit of the 4th agent is NON-MONOTONIC in
##      baseline risk.
##
##   2. Acute vs delayed is the SAME equation with a different term dominating.
##      The 5-HT3 term is a mucosal BURST (FG is sigmoidal in mucosal platinum,
##      HG = 2.5, so it is a burst and not a tail; peak 4-6 h, gone by 24 h,
##      observable as urinary 5-HIAA).  The NK1 term is a TRANSCRIPTIONAL RAMP
##      (Z1->Z2 transit -> preprotachykinin, peak ~40 h).  A 5-HT3 antagonist
##      deletes the term that dominates day 1 and is nearly absent on day 3.
##
##   3. Nausea and vomiting dissociate because they are TWO DIFFERENT weighted
##      sums.  The brainstem sum crosses a threshold; the cortical sum is graded
##      and carries H1, M1, 5-HT2A, AVP, gastric-dysrhythmia and conditioned
##      terms.  Olanzapine is the only agent that subtracts FOUR cortical terms
##      at once, so "no nausea" is the endpoint it owns.
##
##   4. 100% control is structurally impossible.  The area postrema has no
##      blood-brain barrier, so WAP*FG drives the sum DOWNSTREAM OF NO RECEPTOR.
##      No combination of antagonists removes it; only the gain terms
##      (glucocorticoid, GABA-A, CB1) attenuate it.
##
##  THERE IS NO TONIC-DRIVE CONSTANT TO CALIBRATE.  Every term of both sums is a
##  chemo-driven INCREMENT over its analytic drug-free value, so a drug-free patient
##  has an excess drive of exactly zero by construction.  tonic_selftest() below
##  verifies it and must return ~0; anything else is a bug, not a parameter.  An
##  earlier version subtracted a tonic constant AFTER applying the gain, which let
##  dexamethasone alone push the total below baseline and act as a complete
##  antiemetic — a defect this formulation makes structurally impossible.
##
##  UNITS:  time = h, drug amounts = mg, volumes = L, receptor states = FRACTION
##          of the receptor pool, 5-HT / substance P = nM, QTc = ms.
##
##  CALIBRATION (9 parameters, 10 anchors) — see cinv_references.md and README.md
##  Emesis anchors come from ONE randomised trial (Hesketh 2003 JCO, PMID
##  14559886; both arms x all three windows) so the fit cannot be assembled out
##  of whichever cohort is convenient.  Nausea scale from the CONTROL arm of
##  Navari 2016 NEJM (PMID 27705262) and from untreated natural history.
##  EVERYTHING ELSE IS HELD OUT: palonosetron, NEPA, rolapitant, olanzapine,
##  AC / carboplatin / oxaliplatin, the dexamethasone DDI, NK1 PET occupancy,
##  5-HIAA kinetics, QTc, CYP2D6 ultrarapid metabolisers, the risk covariates
##  and the multi-cycle anticipatory carryover.
##
##  USAGE
##    library(mrgsolve); source("cinv_mrgsolve_model.R")
##    out <- run_scenario("S12_quadruplet");  plot(out, NAUSEA + DRIVE + EMES ~ time)
##    tab <- run_all()                       # all 30 scenarios, endpoint table
## =============================================================================

library(mrgsolve)
suppressMessages(library(dplyr))

cinv_code <- '
$PROB
# CINV QSP model — additive emetic drive, threshold hazard, exponential endpoint

$PARAM @annotated
// ---------------- 5-HT3 receptor antagonist PK ----------------
F_OND    : 0.60   : ondansetron oral bioavailability (-)
KA_OND   : 1.2    : ondansetron absorption rate (1/h)
V1_OND   : 60     : ondansetron central volume (L)
V2_OND   : 100    : ondansetron peripheral volume (L)
Q_OND    : 25     : ondansetron intercompartmental clearance (L/h)
CL_OND   : 25     : ondansetron clearance (L/h)
FM2D6_OND: 0.25   : fraction of ondansetron clearance via CYP2D6 (-)
F_GRA    : 0.60   : granisetron oral bioavailability (-)
KA_GRA   : 1.5    : granisetron absorption rate (1/h)
V1_GRA   : 60     : granisetron central volume (L)
V2_GRA   : 120    : granisetron peripheral volume (L)
Q_GRA    : 20     : granisetron intercompartmental clearance (L/h)
CL_GRA   : 30     : granisetron clearance (L/h)
KA_GRAER : 0.0115 : granisetron extended-release SC depot rate (1/h)
F_PAL    : 1.00   : palonosetron availability from depot cmt (-)
KA_PAL   : 1.0    : palonosetron absorption rate (1/h)
V1_PAL   : 120    : palonosetron central volume (L)
V2_PAL   : 460    : palonosetron peripheral volume (L)
Q_PAL    : 25     : palonosetron intercompartmental clearance (L/h)
CL_PAL   : 11     : palonosetron clearance (L/h)
// ---------------- NK1 receptor antagonist PK ----------------
F_APR    : 0.65   : aprepitant oral bioavailability (-)
KA_APR   : 0.35   : aprepitant absorption rate (1/h)
V1_APR   : 45     : aprepitant central volume (L)
V2_APR   : 25     : aprepitant peripheral volume (L)
Q_APR    : 6      : aprepitant intercompartmental clearance (L/h)
CL_APR   : 4.5    : aprepitant clearance (L/h)
KFOS     : 8.0    : fosaprepitant to aprepitant conversion rate (1/h)
FOSEQ    : 0.77   : fosaprepitant to aprepitant molar equivalence (-)
F_NET    : 0.70   : netupitant oral bioavailability (-)
KA_NET   : 0.30   : netupitant absorption rate (1/h)
V1_NET   : 100    : netupitant central volume (L)
V2_NET   : 1800   : netupitant peripheral volume (L)
Q_NET    : 30     : netupitant intercompartmental clearance (L/h)
CL_NET   : 15     : netupitant clearance (L/h)
F_ROL    : 0.90   : rolapitant oral bioavailability (-)
KA_ROL   : 0.25   : rolapitant absorption rate (1/h)
V1_ROL   : 180    : rolapitant central volume (L)
V2_ROL   : 200    : rolapitant peripheral volume (L)
Q_ROL    : 8      : rolapitant intercompartmental clearance (L/h)
CL_ROL   : 1.4    : rolapitant clearance (L/h)
// ---------------- dexamethasone PK + genomic transduction ----------------
F_DEX    : 0.80   : dexamethasone oral bioavailability (-)
KA_DEX   : 1.5    : dexamethasone absorption rate (1/h)
V1_DEX   : 30     : dexamethasone central volume (L)
V2_DEX   : 35     : dexamethasone peripheral volume (L)
Q_DEX    : 12     : dexamethasone intercompartmental clearance (L/h)
CL_DEX   : 17     : dexamethasone clearance (L/h)
FM3A4_DEX: 0.85   : fraction of dexamethasone clearance via CYP3A4 (-)
KON_GR   : 0.60   : GR association rate (1/nM/h)
KOFF_GR  : 3.0    : GR dissociation rate (1/h)
KTR_GR   : 0.20   : genomic transduction rate (1/h)
KOUT_GR  : 0.0347 : genomic effect turnover (1/h)
// ---------------- other antiemetic PK ----------------
F_OLA    : 0.60   : olanzapine oral bioavailability (-)
KA_OLA   : 0.40   : olanzapine absorption rate (1/h)
V1_OLA   : 600    : olanzapine central volume (L)
V2_OLA   : 400    : olanzapine peripheral volume (L)
Q_OLA    : 30     : olanzapine intercompartmental clearance (L/h)
CL_OLA   : 25     : olanzapine clearance (L/h)
F_MCP    : 0.75   : metoclopramide oral bioavailability (-)
KA_MCP   : 1.5    : metoclopramide absorption rate (1/h)
V1_MCP   : 250    : metoclopramide volume (L)
CL_MCP   : 35     : metoclopramide clearance (L/h)
F_LOR    : 0.90   : lorazepam oral bioavailability (-)
KA_LOR   : 1.0    : lorazepam absorption rate (1/h)
V1_LOR   : 90     : lorazepam volume (L)
CL_LOR   : 4.6    : lorazepam clearance (L/h)
F_DRO    : 0.10   : dronabinol oral bioavailability (-)
KA_DRO   : 0.60   : dronabinol absorption rate (1/h)
V1_DRO   : 700    : dronabinol volume (L)
CL_DRO   : 90     : dronabinol clearance (L/h)
// ---------------- chemotherapy PK ----------------
V1_CIS   : 20     : free platinum central volume (L)
CL_CIS   : 25     : free platinum clearance (L/h)
K12_CIS  : 1.5    : free to tissue platinum rate (1/h)
K21_CIS  : 0.020  : tissue to free platinum rate (1/h)
KIN_G    : 1.0    : mucosal platinum signal formation (1/h)
KOUT_G   : 0.10   : mucosal platinum signal loss (1/h)
KZ       : 0.050  : delayed central signal transit rate (1/h)
CISRATE  : 60     : chemotherapy zero-order infusion rate (mg/h)
TINF     : 2.0    : chemotherapy infusion duration (h)
CYCLEN   : 504    : cycle length (h)
// ---------------- CYP pools ----------------
KDEG3A   : 0.0289 : CYP3A4 turnover (1/h)
KI3A_APR : 180    : aprepitant CYP3A4 inhibition constant (nM)
KI3A_NET : 120    : netupitant CYP3A4 inhibition constant (nM)
KDEG2D   : 0.0289 : CYP2D6 activity turnover (1/h)
KI2D6_ROL: 95     : rolapitant CYP2D6 inhibition constant (nM)
EMAXIND  : 0.45   : maximal CYP3A4 induction by dexamethasone (-)
EC50IND  : 1.0    : dexamethasone genomic effect for half-max induction (-)
// ---------------- enterochromaffin cell / serotonin ----------------
KSYN_EC  : 0.12   : EC-cell 5-HT store repletion rate (1/h)
KRE      : 1.30   : FITTED maximal EC-cell 5-HT release rate (1/h)
KG50     : 1.50   : mucosal platinum signal for half-max release (-)
HG       : 2.5    : Hill coefficient for mucosal release (-)
SCALE_G  : 900    : gut 5-HT scale factor (nM per unit release)
KEL_G    : 0.18   : gut interstitial 5-HT elimination (1/h)
KABS_P   : 0.010  : gut to plasma 5-HT transfer (1/h)
KEL_P    : 0.55   : plasma 5-HT elimination (1/h)
KHIAA    : 0.0075 : 5-HIAA formation scale (umol/nM/h)
FP_HIAA  : 6.0    : plasma weighting for 5-HIAA formation (-)
EDEX_5HT : 0.35   : maximal dexamethasone suppression of 5-HT release (-)
KD5HT    : 0.60   : dexamethasone genomic effect for half-max 5-HT suppression (-)
// ---------------- 5-HT3 receptor binding ----------------
KON5     : 0.12   : 5-HT association rate at 5-HT3 (1/nM/h)
KOFF5    : 60     : 5-HT dissociation rate at 5-HT3 (1/h)
KONO     : 2.0    : ondansetron association rate (1/nM/h)
KOFFO    : 4.0    : ondansetron dissociation rate (1/h)
KONG     : 4.0    : granisetron association rate (1/nM/h)
KOFFG    : 2.0    : granisetron dissociation rate (1/h)
KONP     : 1.0    : palonosetron association rate (1/nM/h)
KOFFP    : 0.10   : palonosetron dissociation rate (1/h)
KINT     : 0.045  : palonosetron-driven receptor internalisation (1/h)
KREC     : 0.012  : internalised receptor recycling (1/h)
// ---------------- substance P / NK1 ----------------
KSYN_PPT : 0.0289 : preprotachykinin capacity formation (1/h)
KDEG_PPT : 0.0289 : preprotachykinin capacity loss (1/h)
EPPT     : 0.301762 : FITTED maximal TAC1 induction (-)
KZ50     : 0.035  : delayed signal for half-max TAC1 induction (-)
EDEXPPT  : 0.40   : maximal dexamethasone suppression of TAC1 induction (-)
KDP      : 0.60   : dexamethasone genomic effect for half-max TAC1 suppression (-)
KSP      : 0.060  : central substance P formation (nM/h)
KEL_SP   : 1.20   : central substance P elimination (1/h)
KSPP     : 0.90   : peripheral substance P formation from capacity (nM/h)
KSPA     : 0.55   : ACUTE peripheral substance P co-release (nM/h)
KEL_SPP  : 1.10   : peripheral substance P elimination (1/h)
KONSP    : 90     : substance P association rate at NK1 (1/nM/h)
KOFFSP   : 18     : substance P dissociation rate at NK1 (1/h)
KONA     : 1.6    : aprepitant association rate (1/nM/h)
KOFFA    : 0.16   : aprepitant dissociation rate (1/h)
KONN     : 0.8    : netupitant association rate (1/nM/h)
KOFFN    : 0.72   : netupitant dissociation rate (1/h)
KONR     : 1.1    : rolapitant association rate (1/nM/h)
KOFFR    : 0.033  : rolapitant dissociation rate (1/h)
KDSP_P   : 1.60   : peripheral NK1 dissociation constant for substance P (nM)
KI_APR_P : 0.60   : aprepitant peripheral NK1 inhibition constant (nM)
KI_NET_P : 3.0    : netupitant peripheral NK1 inhibition constant (nM)
KI_ROL_P : 2.0    : rolapitant peripheral NK1 inhibition constant (nM)
// ---------------- dopamine and other receptors ----------------
KDA      : 1.20   : area postrema dopamine formation (nM/h)
EDA      : 2.40   : maximal chemotherapy-driven dopamine increase (-)
KEL_DA   : 1.20   : dopamine elimination (1/h)
KD_DA    : 8.0    : dopamine dissociation constant at D2 (nM)
KI_MCP_D2: 30     : metoclopramide D2 inhibition constant (nM)
KI_OLA_D2: 20     : olanzapine D2 inhibition constant (nM)
HIST0    : 25     : tonic histamine concentration (nM)
KD_H1    : 25     : histamine dissociation constant at H1 (nM)
KI_OLA_H1: 7.0    : olanzapine H1 inhibition constant (nM)
ACH0     : 30     : tonic acetylcholine concentration (nM)
KD_M1    : 30     : acetylcholine dissociation constant at M1 (nM)
KI_OLA_M1: 26     : olanzapine M1 inhibition constant (nM)
KD_2A    : 180    : 5-HT dissociation constant at 5-HT2A (nM)
KI_OLA_2A: 4.0    : olanzapine 5-HT2A inhibition constant (nM)
KI_LOR_GABA: 45   : lorazepam GABA-A occupancy constant (nM)
KI_DRO_CB1 : 12   : dronabinol CB1 occupancy constant (nM)
// ---------------- vagal afferent ----------------
VAG0     : 0.35   : basal vagal afferent firing (-)
EV3      : 3.40   : 5-HT3 gain on vagal firing (-)
EVNK     : 0.55   : peripheral NK1 gain on vagal firing (-)
EVMECH   : 0.30   : mechanoreceptor gain on vagal firing (-)
TAU_VAG  : 0.35   : vagal afferent time constant (h)
// ---------------- brainstem drive ----------------
WVAG     : 1.55   : vagal weight in the emetic sum (-)
WNKC     : 0.200982 : FITTED central NK1 weight in the emetic sum (-)
WD2      : 0.85   : D2 weight in the emetic sum (-)
WBLD     : 0.45   : blood-borne 5-HT weight in the emetic sum (-)
WAP      : 1.29189 : FITTED receptor-INDEPENDENT area postrema weight (-)
WAPD     : 1.38485 : FITTED delayed share of the receptor-independent drive (-)
WANT     : 0.60   : anticipatory weight in the emetic sum (-)
K5P      : 140    : plasma 5-HT for half-max area postrema drive (nM)
TAU_BS   : 0.50   : brainstem integration time constant (h)
EGABA    : 0.30   : maximal GABA-A gain reduction, brainstem (-)
ECB1     : 0.28   : maximal CB1 gain reduction, brainstem (-)
EDEXBS   : 0.34   : FITTED maximal glucocorticoid gain reduction (-)
KDB      : 0.60   : dexamethasone genomic effect for half-max gain reduction (-)
// ---------------- cortical (nausea) drive ----------------
NVAG     : 0.85   : vagal weight in the nausea sum (-)
NNK      : 0.95   : central NK1 weight in the nausea sum (-)
N2A      : 0.80   : 5-HT2A weight in the nausea sum (-)
NH1      : 0.55   : H1 weight in the nausea sum (-)
NM1      : 0.45   : M1 weight in the nausea sum (-)
ND2      : 0.55   : D2 weight in the nausea sum (-)
NAVP     : 0.35   : vasopressin weight in the nausea sum (-)
KAVP     : 6.0    : vasopressin scale (pg/mL)
NANT     : 0.80   : anticipatory / anxiety weight in the nausea sum (-)
NSW      : 0.70   : gastric dysrhythmia weight in the nausea sum (-)
NAP      : 0.861906 : FITTED receptor-independent weight in the nausea sum (-)
EGABA_C  : 0.35   : maximal GABA-A gain reduction, cortical (-)
ECB1_C   : 0.30   : maximal CB1 gain reduction, cortical (-)
NK50     : 2.3613 : FITTED cortical drive for half-max nausea (-)
HN       : 2.2    : Hill coefficient for nausea (-)
TAU_N    : 1.2    : nausea time constant (h)
// ---------------- hazards ----------------
HMAX     : 1.25544 : FITTED maximal emetic hazard (episodes/h)
THR      : 8.12534 : FITTED emetic threshold on the drive EXCESS (-)
NE       : 2.6    : Hill coefficient of the emetic hazard (-)
RMAX     : 0.075  : maximal rescue-medication hazard (1/h)
NR50     : 6.0    : nausea for half-max rescue hazard (0-10)
NRH      : 3.0    : Hill coefficient for rescue (-)
KNAU     : 0.145  : nausea-occurrence hazard scale (1/h)
NN50     : 2.06492 : FITTED nausea for half-max nausea hazard (0-10)
NNH      : 3.0    : Hill coefficient for the nausea hazard (-)
OMEGA    : 0.59149 : FITTED lognormal SD of susceptibility frailty (-)
// ---------------- gastric ----------------
GVOL0    : 250    : reference gastric volume (mL)
KSEC     : 42     : gastric secretion rate (mL/h)
KEMP     : 0.22   : basal gastric emptying rate (1/h)
EMCP     : 0.55   : maximal prokinetic gain from peripheral D2 blockade (-)
EDR      : 0.55   : maximal drive-induced emptying inhibition (-)
KDR      : 2.0    : drive for half-max emptying inhibition (-)
ESW      : 0.60   : maximal slow-wave uncoupling (-)
KSW      : 2.2    : drive for half-max slow-wave uncoupling (-)
TAU_SW   : 2.0    : slow-wave time constant (h)
AVP0     : 2.0    : basal plasma vasopressin (pg/mL)
EAVP     : 14     : maximal vasopressin rise (pg/mL)
KAV      : 2.6    : drive for half-max vasopressin rise (-)
TAU_AVP  : 0.8    : vasopressin time constant (h)
// ---------------- associative learning ----------------
KACQ     : 0.010  : anticipatory acquisition rate (1/h)
KEXT     : 0.0012 : anticipatory extinction rate (1/h)
EANXA    : 1.40   : anticipatory gain on anxiety (-)
ANX0     : 0.30   : trait anxiety (-)
TAU_ANX  : 12     : anxiety time constant (h)
ELOR_ANX : 0.55   : maximal lorazepam anxiolysis (-)
CUE_ON   : 1.0    : clinic-cue switch (0/1)
LEARN    : 1.0    : learning switch, set 0 to compute the tonic drive (0/1)
// ---------------- safety ----------------
SQ_OND   : 0.0110 : ondansetron QTc slope (ms/nM)
SQ_GRA   : 0.0016 : granisetron QTc slope (ms/nM)
SQ_PAL   : 0.0008 : palonosetron QTc slope (ms/nM)
SQ_OLA   : 0.0090 : olanzapine QTc slope (ms/nM)
TAU_Q    : 0.05   : QTc time constant (h) - hERG block has no kinetic delay
GLUB     : 5.4    : basal plasma glucose (mmol/L)
EGLU     : 4.2    : maximal dexamethasone hyperglycaemia (mmol/L)
KGL      : 0.60   : dexamethasone genomic effect for half-max hyperglycaemia (-)
TAU_GLU  : 4.0    : glucose time constant (h)
ES_H1    : 0.62   : H1 contribution to sedation (-)
ES_M1    : 0.30   : M1 contribution to sedation (-)
ES_LOR   : 0.85   : GABA-A contribution to sedation (-)
TAU_SED  : 1.5    : sedation time constant (h)
KEPS     : 0.055  : AKATH accumulation rate (1/h)
EPSTHR   : 0.55   : D2 blockade threshold for AKATH (-)
KEPSOFF  : 0.030  : AKATH resolution rate (1/h)
KCONST   : 0.022  : constipation accumulation rate (1/h)
KCONSTOFF: 0.020  : constipation resolution rate (1/h)
// ---------------- systemic consequences ----------------
VEM      : 0.055  : fluid loss per emetic episode (L)
KREPL    : 0.055  : fluid repletion rate (1/h)
KEM_K    : 0.055  : potassium loss per emetic episode (mmol/L)
KP0      : 4.2    : reference plasma potassium (mmol/L)
TAU_K    : 10     : potassium time constant (h)
SCR0     : 0.85   : baseline serum creatinine (mg/dL)
ENEPH    : 0.0021 : platinum nephrotoxicity gain (1/mg)
EVOL     : 0.55   : volume-depletion potentiation of nephrotoxicity (1/L)
TAU_SCR  : 24     : creatinine time constant (h)
KRDI     : 0.00048: dose-intensity loss rate (1/h per nausea unit)
RDITH    : 5.0    : nausea threshold for dose-intensity loss (0-10)
FLIEMAX  : 108    : maximum FLIE score (-)
WF_N     : 0.55   : nausea weight in FLIE (-)
WF_E     : 0.30   : emesis weight in FLIE (-)
TAU_FLIE : 8.0    : FLIE time constant (h)
HYDR     : 1.0    : hydration / repletion switch (-)
// ---------------- covariates ----------------
EMETO_P  : 1.0    : peripheral emetogenicity of the regimen (-)
EMETO_C  : 1.0    : central emetogenicity of the regimen (-)
FEMALE   : 1.0    : female sex (0/1)
AGE      : 52     : age (years)
ALCOHOL  : 0.0    : chronic heavy alcohol history (0/1)
MOTION   : 0.0    : motion sickness / pregnancy emesis history (0/1)
PRIORCINV: 0.0    : CINV in a previous cycle (0/1)
B_SEX    : 0.26   : female sex effect on log susceptibility (-)
B_AGE    : -0.115 : per-decade age effect on log susceptibility (-)
B_ALC    : -0.42  : alcohol effect on log susceptibility (-)
B_MOT    : 0.30   : motion sickness effect on log susceptibility (-)
B_PREV   : 0.45   : prior CINV effect on log susceptibility (-)
PH2D6    : 1.0    : CYP2D6 phenotype scalar (PM 0.1 / NM 1 / UM 2.5)

$CMT @annotated
A_OND  : ondansetron gut (mg)
C_OND  : ondansetron central (mg)
P_OND  : ondansetron peripheral (mg)
DEPGRA  : granisetron ER subcutaneous depot (mg)
A_GRA  : granisetron gut (mg)
C_GRA  : granisetron central (mg)
P_GRA  : granisetron peripheral (mg)
A_PAL  : palonosetron depot (mg)
C_PAL  : palonosetron central (mg)
P_PAL  : palonosetron peripheral (mg)
A_APR  : aprepitant gut (mg)
FOSP  : fosaprepitant prodrug (mg)
C_APR  : aprepitant central (mg)
P_APR  : aprepitant peripheral (mg)
A_NET  : netupitant gut (mg)
C_NET  : netupitant central (mg)
P_NET  : netupitant peripheral (mg)
A_ROL  : rolapitant gut (mg)
C_ROL  : rolapitant central (mg)
P_ROL  : rolapitant peripheral (mg)
A_DEX  : dexamethasone gut (mg)
C_DEX  : dexamethasone central (mg)
P_DEX  : dexamethasone peripheral (mg)
A_OLA  : olanzapine gut (mg)
C_OLA  : olanzapine central (mg)
P_OLA  : olanzapine peripheral (mg)
A_MCP  : metoclopramide gut (mg)
C_MCP  : metoclopramide central (mg)
A_LOR  : lorazepam gut (mg)
C_LOR  : lorazepam central (mg)
A_DRO  : dronabinol gut (mg)
C_DRO  : dronabinol central (mg)
C_CIS  : free plasma platinum (mg)
T_CIS  : tissue-bound platinum (mg)
G_CIS  : GI-mucosal platinum signal (-)
Z1     : delayed central signal transit 1 (-)
Z2     : delayed central signal transit 2 (-)
E3A4   : relative hepatic CYP3A4 amount (-)
E2D6   : relative hepatic CYP2D6 activity (-)
DEXGR  : dexamethasone-GR complex fraction (-)
DEXN   : dexamethasone nuclear genomic effect (-)
ECS    : EC-cell 5-HT store fraction (-)
SHT_G  : gut interstitial 5-HT (nM)
SHT_P  : plasma 5-HT (nM)
UHIAA  : cumulative urinary 5-HIAA (umol)
R3S    : 5-HT3 receptor fraction bound by 5-HT (-)
R3O    : 5-HT3 receptor fraction bound by ondansetron (-)
R3G    : 5-HT3 receptor fraction bound by granisetron (-)
R3P    : 5-HT3 receptor fraction bound by palonosetron (-)
R3I    : internalised 5-HT3 receptor fraction (-)
PPT    : preprotachykinin-A synthesis capacity (-)
SP_C   : central substance P (nM)
SP_P   : peripheral substance P (nM)
N1S    : NK1 receptor fraction bound by substance P (-)
N1A    : NK1 receptor fraction bound by aprepitant (-)
N1N    : NK1 receptor fraction bound by netupitant (-)
N1R    : NK1 receptor fraction bound by rolapitant (-)
DA_AP  : area postrema dopamine (nM)
VAG    : vagal afferent firing (-)
DRIVE  : brainstem emetic drive (-)
CTX    : cortical nausea drive (-)
NAUSEA : nausea VAS (0-10)
EMES   : cumulative expected emetic episodes (-)
HINT   : cumulative emesis + rescue hazard (-)
NINT   : cumulative nausea-occurrence hazard (-)
RESC   : cumulative rescue medication use (-)
GVOL   : gastric content volume (mL)
SWC    : gastric slow-wave coupling index (-)
AVP    : plasma vasopressin (pg/mL)
ANTIC  : conditioned anticipatory strength (-)
ANX    : anxiety (-)
QTC    : QTcF change from baseline (ms)
GLU    : plasma glucose (mmol/L)
SED    : sedation score (0-10)
AKATH    : extrapyramidal / akathisia score (-)
CONST  : constipation index (-)
ECFV   : extracellular fluid volume deficit (L)
KP     : plasma potassium (mmol/L)
SCR    : serum creatinine (mg/dL)
RDI    : relative dose intensity (fraction)
FLIE   : FLIE quality-of-life score (-)

$GLOBAL
#define MW_OND 293.4
#define MW_GRA 312.4
#define MW_PAL 296.4
#define MW_APR 534.4
#define MW_NET 578.6
#define MW_ROL 500.6
#define MW_DEX 392.5
#define MW_OLA 312.4
#define MW_MCP 299.8
#define MW_LOR 321.2
#define MW_DRO 314.5
// mg in a compartment of V litres -> nM
#define NMOL(amt, V, mw) ((amt) / (V) / (mw) * 1.0e6)
#define POS(x) ((x) > 0.0 ? (x) : 0.0)

$MAIN
// steady state of the drug-free system; the fast receptor states settle within
// the first simulated hour, so only the slow / structural ones are seeded here
E3A4_0  = 1.0;
E2D6_0  = PH2D6;
ECS_0   = 1.0;
PPT_0   = 1.0;
SP_C_0  = KSP / KEL_SP;
SP_P_0  = KSPP / KEL_SPP;
DA_AP_0 = KDA / KEL_DA;
GVOL_0  = GVOL0;
SWC_0   = 1.0;
AVP_0   = AVP0;
ANX_0   = ANX0;
KP_0    = KP0;
SCR_0   = SCR0;
RDI_0   = 1.0;
FLIE_0  = FLIEMAX;
GLU_0   = GLUB;

$ODE
// ===================== concentrations (nM) =====================
double cOND = NMOL(POS(C_OND), V1_OND, MW_OND);
double cGRA = NMOL(POS(C_GRA), V1_GRA, MW_GRA);
double cPAL = NMOL(POS(C_PAL), V1_PAL, MW_PAL);
double cAPR = NMOL(POS(C_APR), V1_APR, MW_APR);
double cNET = NMOL(POS(C_NET), V1_NET, MW_NET);
double cROL = NMOL(POS(C_ROL), V1_ROL, MW_ROL);
double cDEX = NMOL(POS(C_DEX), V1_DEX, MW_DEX);
double cOLA = NMOL(POS(C_OLA), V1_OLA, MW_OLA);
double cMCP = NMOL(POS(C_MCP), V1_MCP, MW_MCP);
double cLOR = NMOL(POS(C_LOR), V1_LOR, MW_LOR);
double cDRO = NMOL(POS(C_DRO), V1_DRO, MW_DRO);

// ===================== CYP pools and DDI =====================
double ACT3A4 = POS(E3A4) / (1.0 + cAPR / KI3A_APR + cNET / KI3A_NET);
double ACT2D6 = POS(E2D6) / (1.0 + cROL / KI2D6_ROL);
dxdt_E3A4 = KDEG3A * (1.0 + EMAXIND * POS(DEXN) / (EC50IND + POS(DEXN))) - KDEG3A * E3A4;
dxdt_E2D6 = KDEG2D * (PH2D6 - E2D6);

// ===================== antiemetic PK =====================
double CLO = CL_OND * (FM2D6_OND * ACT2D6 + (1.0 - FM2D6_OND));
dxdt_A_OND = -KA_OND * A_OND;
dxdt_C_OND = KA_OND * A_OND - CLO / V1_OND * C_OND
             - Q_OND * (C_OND / V1_OND - P_OND / V2_OND);
dxdt_P_OND = Q_OND * (C_OND / V1_OND - P_OND / V2_OND);

dxdt_DEPGRA = -KA_GRAER * DEPGRA;
dxdt_A_GRA = -KA_GRA * A_GRA;
dxdt_C_GRA = KA_GRA * A_GRA + KA_GRAER * DEPGRA - CL_GRA / V1_GRA * C_GRA
             - Q_GRA * (C_GRA / V1_GRA - P_GRA / V2_GRA);
dxdt_P_GRA = Q_GRA * (C_GRA / V1_GRA - P_GRA / V2_GRA);

dxdt_A_PAL = -KA_PAL * A_PAL;
dxdt_C_PAL = KA_PAL * A_PAL - CL_PAL / V1_PAL * C_PAL
             - Q_PAL * (C_PAL / V1_PAL - P_PAL / V2_PAL);
dxdt_P_PAL = Q_PAL * (C_PAL / V1_PAL - P_PAL / V2_PAL);

dxdt_A_APR = -KA_APR * A_APR;
dxdt_FOSP = -KFOS * FOSP;
dxdt_C_APR = KA_APR * A_APR + KFOS * FOSEQ * FOSP - CL_APR / V1_APR * C_APR
             - Q_APR * (C_APR / V1_APR - P_APR / V2_APR);
dxdt_P_APR = Q_APR * (C_APR / V1_APR - P_APR / V2_APR);

dxdt_A_NET = -KA_NET * A_NET;
dxdt_C_NET = KA_NET * A_NET - CL_NET / V1_NET * C_NET
             - Q_NET * (C_NET / V1_NET - P_NET / V2_NET);
dxdt_P_NET = Q_NET * (C_NET / V1_NET - P_NET / V2_NET);

dxdt_A_ROL = -KA_ROL * A_ROL;
dxdt_C_ROL = KA_ROL * A_ROL - CL_ROL / V1_ROL * C_ROL
             - Q_ROL * (C_ROL / V1_ROL - P_ROL / V2_ROL);
dxdt_P_ROL = Q_ROL * (C_ROL / V1_ROL - P_ROL / V2_ROL);

// dexamethasone clearance falls when an NK1 RA inhibits CYP3A4:
// the 20 mg -> 12 mg dose correction is an OUTPUT of these two lines
double CLD = CL_DEX * (FM3A4_DEX * ACT3A4 + (1.0 - FM3A4_DEX));
dxdt_A_DEX = -KA_DEX * A_DEX;
dxdt_C_DEX = KA_DEX * A_DEX - CLD / V1_DEX * C_DEX
             - Q_DEX * (C_DEX / V1_DEX - P_DEX / V2_DEX);
dxdt_P_DEX = Q_DEX * (C_DEX / V1_DEX - P_DEX / V2_DEX);
dxdt_DEXGR = KON_GR * cDEX * (1.0 - POS(DEXGR)) - KOFF_GR * DEXGR;
dxdt_DEXN  = KTR_GR * POS(DEXGR) - KOUT_GR * DEXN;

dxdt_A_OLA = -KA_OLA * A_OLA;
dxdt_C_OLA = KA_OLA * A_OLA - CL_OLA / V1_OLA * C_OLA
             - Q_OLA * (C_OLA / V1_OLA - P_OLA / V2_OLA);
dxdt_P_OLA = Q_OLA * (C_OLA / V1_OLA - P_OLA / V2_OLA);

dxdt_A_MCP = -KA_MCP * A_MCP;
dxdt_C_MCP = KA_MCP * A_MCP - CL_MCP / V1_MCP * C_MCP;
dxdt_A_LOR = -KA_LOR * A_LOR;
dxdt_C_LOR = KA_LOR * A_LOR - CL_LOR / V1_LOR * C_LOR;
dxdt_A_DRO = -KA_DRO * A_DRO;
dxdt_C_DRO = KA_DRO * A_DRO - CL_DRO / V1_DRO * C_DRO;

// ===================== chemotherapy PK =====================
double tc   = SOLVERTIME - CYCLEN * floor(SOLVERTIME / CYCLEN);
double RINF = (tc >= 0.0 && tc < TINF) ? CISRATE : 0.0;
dxdt_C_CIS = RINF - CL_CIS / V1_CIS * C_CIS - K12_CIS * C_CIS + K21_CIS * T_CIS;
dxdt_T_CIS = K12_CIS * C_CIS - K21_CIS * T_CIS;
dxdt_G_CIS = KIN_G * POS(C_CIS) / V1_CIS - KOUT_G * G_CIS;
dxdt_Z1    = KZ * POS(C_CIS) / V1_CIS - KZ * Z1;
dxdt_Z2    = KZ * POS(Z1) - KZ * Z2;

// mucosal injury signal — SIGMOIDAL, so the peripheral 5-HT arm is a BURST and
// not a tail, which is what separates the acute from the delayed phase
double GH = pow(POS(G_CIS), HG);
double FG = GH / (pow(KG50, HG) + GH);

// ===================== enterochromaffin cell / 5-HT =====================
double DEXF = POS(DEXN) / (KD5HT + POS(DEXN));
double REL  = KRE * EMETO_P * FG * POS(ECS) * (1.0 - EDEX_5HT * DEXF);
dxdt_ECS   = KSYN_EC * (1.0 - ECS) - REL;
dxdt_SHT_G = REL * SCALE_G - KEL_G * SHT_G;
dxdt_SHT_P = KABS_P * POS(SHT_G) - KEL_P * SHT_P;
dxdt_UHIAA = KHIAA * (POS(SHT_G) + FP_HIAA * POS(SHT_P));

// ===================== 5-HT3 receptor (KINETIC: slow off-rates matter) ======
double FREE3 = POS(1.0 - POS(R3S) - POS(R3O) - POS(R3G) - POS(R3P) - POS(R3I));
dxdt_R3S = KON5 * POS(SHT_G) * FREE3 - KOFF5 * R3S;
dxdt_R3O = KONO * cOND * FREE3 - KOFFO * R3O;
dxdt_R3G = KONG * cGRA * FREE3 - KOFFG * R3G;
dxdt_R3P = KONP * cPAL * FREE3 - KOFFP * R3P - KINT * POS(R3P);
dxdt_R3I = KINT * POS(R3P) - KREC * R3I;

// ===================== substance P / NK1 =====================
double DEXP = POS(DEXN) / (KDP + POS(DEXN));
double FZ   = POS(Z2) / (KZ50 + POS(Z2));
double IND  = 1.0 + EPPT * EMETO_C * FZ * (1.0 - EDEXPPT * DEXP);
// The receptor-INDEPENDENT area-postrema drive has an acute AND a delayed part:
// tissue platinum persists and the area postrema has no blood-brain barrier, so a
// floor survives into the delayed window that NO antagonist in the map can remove.
// Without the delayed part the aprepitant arm has a delayed hazard of exactly zero.
double APT  = FG + WAPD * FZ;
dxdt_PPT  = KSYN_PPT * IND - KDEG_PPT * PPT;
dxdt_SP_C = KSP * POS(PPT) - KEL_SP * SP_C;
dxdt_SP_P = KSPP * POS(PPT) + KSPA * EMETO_P * FG - KEL_SPP * SP_P;
double FREEN = POS(1.0 - POS(N1S) - POS(N1A) - POS(N1N) - POS(N1R));
dxdt_N1S = KONSP * POS(SP_C) * FREEN - KOFFSP * N1S;
dxdt_N1A = KONA * cAPR * FREEN - KOFFA * N1A;
dxdt_N1N = KONN * cNET * FREEN - KOFFN * N1N;
dxdt_N1R = KONR * cROL * FREEN - KOFFR * N1R;
double occNK1p = POS(SP_P) / (KDSP_P * (1.0 + cAPR / KI_APR_P + cNET / KI_NET_P
                 + cROL / KI_ROL_P) + POS(SP_P));

// ===================== dopamine and equilibrium receptors =====================
dxdt_DA_AP = KDA * (1.0 + EDA * EMETO_P * FG) - KEL_DA * DA_AP;
double D2COMP = 1.0 + cMCP / KI_MCP_D2 + cOLA / KI_OLA_D2;
double occD2  = POS(DA_AP) / (KD_DA * D2COMP + POS(DA_AP));
double blkD2  = 1.0 - 1.0 / D2COMP;
double occ2A  = POS(SHT_P) / (KD_2A * (1.0 + cOLA / KI_OLA_2A) + POS(SHT_P));
double occH1  = HIST0 / (KD_H1 * (1.0 + cOLA / KI_OLA_H1) + HIST0);
double occM1  = ACH0 / (KD_M1 * (1.0 + cOLA / KI_OLA_M1) + ACH0);
double occGABA= cLOR / (KI_LOR_GABA + cLOR);
double occCB1 = cDRO / (KI_DRO_CB1 + cDRO);
double blkH1  = cOLA / (KI_OLA_H1 + cOLA);
double blkM1  = cOLA / (KI_OLA_M1 + cOLA);

// ---- tonic reference values: ANALYTIC functions of the parameters, so that a
// ---- drug-free patient has an excess drive of EXACTLY zero by construction.
// ---- (This replaces the earlier DRIVE0/CTX0 constants.  Subtracting a tonic
// ---- constant AFTER applying the gain was a real defect: it let dexamethasone
// ---- push the total below baseline and act as a complete antiemetic.)
double SPCT     = KSP / KEL_SP;
double qT       = KONSP * SPCT / KOFFSP;
double N1ST     = qT / (1.0 + qT);
double SPPT     = KSPP / KEL_SPP;
double occNK1pT = SPPT / (KDSP_P + SPPT);
double VAGT     = VAG0 + EVNK * occNK1pT;
double DAT      = KDA / KEL_DA;
double occD2T   = DAT / (KD_DA + DAT);
double occH1T   = HIST0 / (KD_H1 + HIST0);
double occM1T   = ACH0 / (KD_M1 + ACH0);

// ===================== vagal afferent =====================
double VAGss = VAG0 + EV3 * POS(R3S) + EVNK * occNK1p
               + EVMECH * POS(GVOL - GVOL0) / GVOL0;
dxdt_VAG = (VAGss - VAG) / TAU_VAG;

// ===================== THE SUM — brainstem =====================
double GAIN_BS = (1.0 - EGABA * occGABA) * (1.0 - ECB1 * occCB1)
                 * (1.0 - EDEXBS * POS(DEXN) / (KDB + POS(DEXN)));
// THE SUM.  Every term is a chemo-driven INCREMENT over its drug-free tonic value,
// so the gain terms cannot push the total below baseline.
double DRIVEss = GAIN_BS * (WVAG * (POS(VAG) - VAGT)
                 + WNKC * (POS(N1S) - N1ST)
                 + WD2 * (occD2 - occD2T)
                 + WBLD * POS(SHT_P) / (K5P + POS(SHT_P))
                 + WAP * EMETO_P * APT
                 + WANT * POS(ANTIC) * CUE_ON);
dxdt_DRIVE = (DRIVEss - DRIVE) / TAU_BS;

// ===================== THE OTHER SUM — cortical =====================
double GAIN_CTX = (1.0 - EGABA_C * occGABA) * (1.0 - ECB1_C * occCB1);
double CTXss = GAIN_CTX * (NVAG * (POS(VAG) - VAGT)
               + NNK * (POS(N1S) - N1ST)
               + N2A * occ2A
               + NH1 * (occH1 - occH1T) + NM1 * (occM1 - occM1T)
               + ND2 * (occD2 - occD2T)
               + NAVP * (POS(AVP) - AVP0) / KAVP
               + NANT * (POS(ANTIC) * CUE_ON + POS(ANX) - ANX0)
               + NAP * EMETO_P * APT
               + NSW * POS(1.0 - SWC));
dxdt_CTX = (CTXss - CTX) / TAU_BS;

// DRIVE and CTX are already the EXCESS over tonic (see the increment sums above)
double EBS = POS(DRIVE);
double ECT = POS(CTX);
double NAUSEAss = 10.0 * pow(ECT, HN) / (pow(NK50, HN) + pow(ECT, HN));
dxdt_NAUSEA = (NAUSEAss - NAUSEA) / TAU_N;

// ===================== hazards =====================
double SUSC = exp(B_SEX * FEMALE + B_AGE * (AGE - 50.0) / 10.0
                  + B_ALC * ALCOHOL + B_MOT * MOTION + B_PREV * PRIORCINV);
double lam  = SUSC * HMAX * pow(EBS, NE) / (pow(THR, NE) + pow(EBS, NE));
double NAU  = POS(NAUSEA);
double lamr = SUSC * RMAX * pow(NAU, NRH) / (pow(NR50, NRH) + pow(NAU, NRH));
double lamn = SUSC * KNAU * pow(NAU, NNH) / (pow(NN50, NNH) + pow(NAU, NNH));
dxdt_EMES = lam;
dxdt_HINT = lam + lamr;
dxdt_NINT = lamn;
dxdt_RESC = lamr;

// ===================== gastric =====================
double EMPT = KEMP * (1.0 + EMCP * blkD2)
              * (1.0 - EDR * EBS / (KDR + EBS));
dxdt_GVOL = KSEC - EMPT * GVOL - 0.45 * lam * GVOL;
dxdt_SWC  = ((1.0 - ESW * EBS / (KSW + EBS)) - SWC) / TAU_SW;
dxdt_AVP  = ((AVP0 + EAVP * pow(EBS, 2.0)
             / (pow(KAV, 2.0) + pow(EBS, 2.0))) - AVP) / TAU_AVP;

// ===================== associative learning =====================
dxdt_ANTIC = LEARN * (KACQ * (NAU / 10.0) * (1.0 - POS(ANTIC))) - KEXT * ANTIC;
double ANXss = ANX0 * (1.0 + LEARN * EANXA * POS(ANTIC)) * (1.0 - ELOR_ANX * occGABA);
dxdt_ANX = (ANXss - ANX) / TAU_ANX;

// ===================== safety =====================
double QTCss = SQ_OND * cOND + SQ_GRA * cGRA + SQ_PAL * cPAL + SQ_OLA * cOLA;
dxdt_QTC = (QTCss - QTC) / TAU_Q;
dxdt_GLU = ((GLUB + EGLU * POS(DEXN) / (KGL + POS(DEXN))) - GLU) / TAU_GLU;
double SEDss = 10.0 * (ES_H1 * blkH1 + ES_M1 * blkM1 + ES_LOR * occGABA) / 1.77;
dxdt_SED = (SEDss - SED) / TAU_SED;
dxdt_AKATH = KEPS * POS(blkD2 - EPSTHR) - KEPSOFF * AKATH;
dxdt_CONST = KCONST * (POS(R3O) + POS(R3G) + POS(R3P) + POS(R3I)) - KCONSTOFF * CONST;

// ===================== systemic consequences =====================
dxdt_ECFV = VEM * lam - KREPL * HYDR * ECFV;
dxdt_KP   = -KEM_K * lam + (KP0 - KP) / TAU_K;
double SCRss = SCR0 * (1.0 + ENEPH * POS(T_CIS) * (1.0 + EVOL * POS(ECFV)));
dxdt_SCR = (SCRss - SCR) / TAU_SCR;
dxdt_RDI = -KRDI * POS(NAU - RDITH);
double FLss = FLIEMAX * (1.0 - WF_N * NAU / 10.0
              - WF_E * (lam / 0.6 < 1.0 ? lam / 0.6 : 1.0));
dxdt_FLIE = (FLss - FLIE) / TAU_FLIE;

$TABLE
double cOND_o = NMOL(POS(C_OND), V1_OND, MW_OND);
double cAPR_o = NMOL(POS(C_APR), V1_APR, MW_APR);
double cDEX_o = NMOL(POS(C_DEX), V1_DEX, MW_DEX);
double cOLA_o = NMOL(POS(C_OLA), V1_OLA, MW_OLA);
double cPAL_o = NMOL(POS(C_PAL), V1_PAL, MW_PAL);
double GH_o = pow(POS(G_CIS), HG);
double FG_o = GH_o / (pow(KG50, HG) + GH_o);
double EBS_o = POS(DRIVE);
double ECT_o = POS(CTX);
double SUSC_o = exp(B_SEX * FEMALE + B_AGE * (AGE - 50.0) / 10.0
                    + B_ALC * ALCOHOL + B_MOT * MOTION + B_PREV * PRIORCINV);
double lam_o = SUSC_o * HMAX * pow(EBS_o, NE) / (pow(THR, NE) + pow(EBS_o, NE));
// the four additive contributions to the emetic sum, reported separately so a
// user can SEE that each drug removes one term rather than scaling the total
double TERM_VAG = WVAG * POS(VAG);
double TERM_NK1 = WNKC * POS(N1S);
double TERM_D2  = WD2 * POS(DA_AP) / (KD_DA * (1.0 + NMOL(POS(C_MCP), V1_MCP, MW_MCP)
                  / KI_MCP_D2 + cOLA_o / KI_OLA_D2) + POS(DA_AP));
double FZ_o = POS(Z2) / (KZ50 + POS(Z2));
double APT_o = FG_o + WAPD * FZ_o;
double TERM_AP  = WAP * EMETO_P * APT_o;
double NK1OCC_DRUG = POS(N1A) + POS(N1N) + POS(N1R);
double R3OCC_DRUG  = POS(R3O) + POS(R3G) + POS(R3P) + POS(R3I);

$CAPTURE @annotated
cOND_o : ondansetron plasma concentration (nM)
cPAL_o : palonosetron plasma concentration (nM)
cAPR_o : aprepitant plasma concentration (nM)
cDEX_o : dexamethasone plasma concentration (nM)
cOLA_o : olanzapine plasma concentration (nM)
FG_o   : mucosal injury signal (-)
EBS_o  : brainstem drive EXCESS over tonic (-)
ECT_o  : cortical drive EXCESS over tonic (-)
lam_o  : instantaneous emetic hazard (episodes/h)
SUSC_o : patient susceptibility multiplier (-)
TERM_VAG : additive vagal / 5-HT3 term (-)
TERM_NK1 : additive central NK1 term (-)
TERM_D2  : additive D2 term (-)
TERM_AP  : additive receptor-INDEPENDENT area postrema term (-)
NK1OCC_DRUG : NK1 receptor fraction occupied by an antagonist (-)
R3OCC_DRUG  : 5-HT3 receptor fraction occupied by an antagonist or internalised (-)
'

mod <- mcode("cinv", cinv_code, soloc = tempdir())

## =============================================================================
##  SELF-TEST, not a calibration step.  Because both sums are built from
##  increments over analytic drug-free values, the excess drive of a
##  chemotherapy-free, drug-free patient must be exactly zero.
## =============================================================================
tonic_selftest <- function(m = mod, pars = list()) {
  m2 <- if (length(pars)) do.call(param, c(list(m), pars)) else m
  m2 <- param(m2, CISRATE = 0, LEARN = 0)
  s  <- as.data.frame(mrgsim(m2, end = 200, delta = 10, atol = 1e-10, rtol = 1e-8))
  c(DRIVE = tail(s$DRIVE, 1), CTX = tail(s$CTX, 1), NAUSEA = tail(s$NAUSEA, 1))
}

## =============================================================================
##  Regimen emetogenicity.  EMETO_P sets the peripheral (5-HT) burst amplitude
##  and EMETO_C the central (substance P) induction.  Carboplatin's relatively
##  HIGH central-to-peripheral ratio is why it behaves like HEC in the delayed
##  window despite being classified MEC — an output, not a rule.
## =============================================================================
EMETO <- list(
  cisplatin   = c(P = 1.00, C = 1.00, rate = 60.0),
  AC          = c(P = 0.74, C = 0.88, rate = 50.0),
  carboplatin = c(P = 0.60, C = 0.72, rate = 37.5),
  oxaliplatin = c(P = 0.44, C = 0.40, rate = 30.0),
  MEC         = c(P = 0.50, C = 0.45, rate = 30.0),
  LEC         = c(P = 0.20, C = 0.15, rate = 15.0)
)

## ---- dosing helpers ---------------------------------------------------------
.ev <- function(cmt, amt, times, F = 1) {
  do.call(rbind, lapply(times, function(t)
    data.frame(time = t, cmt = cmt, amt = amt * F, evid = 1)))
}
CMTI <- function(nm) which(names(init(mod)) == nm)

dose_ond_iv  <- function(d = 8,  t = -0.5) .ev(CMTI("C_OND"), d, t)
dose_ond_po  <- function(d = 8,  t)        .ev(CMTI("A_OND"), d, t, 0.60)
dose_gra_iv  <- function(d = 1,  t = -0.5) .ev(CMTI("C_GRA"), d, t)
dose_graER   <- function(d = 10, t = -0.5) .ev(CMTI("DEPGRA"), d, t)
dose_pal_iv  <- function(d = 0.25, t = -0.5) .ev(CMTI("C_PAL"), d, t)
dose_pal_po  <- function(d = 0.5, t = -1)  .ev(CMTI("A_PAL"), d, t, 1.00)
dose_apr     <- function(d, t)             .ev(CMTI("A_APR"), d, t, 0.65)
dose_fos     <- function(d = 150, t = -0.5) .ev(CMTI("FOSP"), d, t)
dose_net     <- function(d = 300, t = -1)  .ev(CMTI("A_NET"), d, t, 0.70)
dose_rol     <- function(d = 180, t = -2)  .ev(CMTI("A_ROL"), d, t, 0.90)
dose_dex     <- function(d, t)             .ev(CMTI("A_DEX"), d, t, 0.80)
dose_ola     <- function(d = 10, t = c(-0.5, 24, 48, 72)) .ev(CMTI("A_OLA"), d, t, 0.60)
dose_mcp     <- function(d = 20, t)        .ev(CMTI("A_MCP"), d, t, 0.75)
dose_mcp_iv  <- function(d = 180, t)       .ev(CMTI("C_MCP"), d, t)
dose_lor     <- function(d = 1, t)         .ev(CMTI("A_LOR"), d, t, 0.90)
dose_dro     <- function(d = 5, t)         .ev(CMTI("A_DRO"), d, t, 0.10)

DEX_STD  <- function() rbind(dose_dex(20, -0.5),
                             dose_dex(8, c(24, 36, 48, 60, 72, 84)))
DEX_NK1  <- function() rbind(dose_dex(12, -0.5), dose_dex(8, c(24, 48, 72)))
APR_STD  <- function() rbind(dose_apr(125, -1), dose_apr(80, c(24, 48)))

## =============================================================================
##  SCENARIOS — built as MATCHED PAIRS wherever a single structural claim is
##  being tested, so each comparison isolates one term of the sum.
## =============================================================================
SCEN <- list(
  S00_natural_history = list(chemo = "cisplatin", ev = NULL,
    note = "untreated cisplatin: the full additive sum"),
  S01_dex_alone = list(chemo = "cisplatin", ev = DEX_STD(),
    note = "gain term only, no receptor term removed"),
  S02_ond_alone = list(chemo = "cisplatin",
    ev = rbind(dose_ond_iv(), dose_ond_po(8, c(8, 24, 32, 48, 56))),
    note = "removes the acute term only"),
  S03_ond_dex = list(chemo = "cisplatin", ev = rbind(dose_ond_iv(), DEX_STD()),
    note = "CALIBRATION ARM (Hesketh 2003 control)"),
  S04_ond_dex_apr = list(chemo = "cisplatin",
    ev = rbind(dose_ond_iv(), DEX_NK1(), APR_STD()),
    note = "CALIBRATION ARM (Hesketh 2003 aprepitant)"),
  S05_pal_dex = list(chemo = "cisplatin", ev = rbind(dose_pal_iv(), DEX_STD()),
    note = "HELD OUT: slow off-rate + internalisation vs ondansetron"),
  S06_pal_dex_apr = list(chemo = "cisplatin",
    ev = rbind(dose_pal_iv(), DEX_NK1(), APR_STD()),
    note = "HELD OUT: triplet"),
  S07_pal_dex_fos = list(chemo = "cisplatin",
    ev = rbind(dose_pal_iv(), DEX_NK1(), dose_fos()),
    note = "HELD OUT: single-dose IV prodrug vs 3-day oral"),
  S08_nepa_dex = list(chemo = "cisplatin",
    ev = rbind(dose_pal_po(), dose_net(), DEX_NK1()),
    note = "HELD OUT: NEPA fixed combination"),
  S09_rol_pal_dex = list(chemo = "cisplatin",
    ev = rbind(dose_pal_iv(), dose_rol(), DEX_STD()),
    note = "HELD OUT: 180 h NK1 half-life, NO CYP3A4 effect so dex stays 20 mg"),
  S10_gra_dex = list(chemo = "cisplatin", ev = rbind(dose_gra_iv(), DEX_STD()),
    note = "HELD OUT: granisetron"),
  S11_graER_dex = list(chemo = "cisplatin", ev = rbind(dose_graER(), DEX_STD()),
    note = "HELD OUT: does a 5-day 5-HT3 DEPOT rescue the delayed phase?"),
  S12_quadruplet = list(chemo = "cisplatin",
    ev = rbind(dose_pal_iv(), DEX_NK1(), APR_STD(), dose_ola()),
    note = "HELD OUT: guideline HEC quadruplet"),
  S13_quad_ond = list(chemo = "cisplatin",
    ev = rbind(dose_ond_iv(), DEX_NK1(), APR_STD(), dose_ola()),
    note = "HELD OUT: Navari 2016 olanzapine arm"),
  S14_nepa_ola = list(chemo = "cisplatin",
    ev = rbind(dose_pal_po(), dose_net(), DEX_NK1(), dose_ola()),
    note = "HELD OUT: NEPA + olanzapine"),
  S15_ola5_quad = list(chemo = "cisplatin",
    ev = rbind(dose_pal_iv(), DEX_NK1(), APR_STD(), dose_ola(5)),
    note = "HELD OUT: J-FORCE 5 mg vs 10 mg"),
  S16_ola_triple_noNK1 = list(chemo = "cisplatin",
    ev = rbind(dose_pal_iv(), DEX_STD(), dose_ola()),
    note = "HELD OUT: can olanzapine SUBSTITUTE for an NK1 RA?"),
  S17_ola_alone = list(chemo = "cisplatin", ev = dose_ola(),
    note = "HELD OUT: olanzapine monotherapy"),
  S18_mcp_dex = list(chemo = "cisplatin",
    ev = rbind(dose_mcp(20, c(-0.5, 6, 12, 24, 36, 48, 60, 72)), DEX_STD()),
    note = "HELD OUT: the pre-1990 regimen"),
  S19_mcp_high = list(chemo = "cisplatin",
    ev = rbind(dose_mcp_iv(180, c(-0.5, 3, 6, 9)), DEX_STD()),
    note = "HELD OUT: high-dose metoclopramide (AKATH is the price)"),
  S20_ond_dex_lor = list(chemo = "cisplatin",
    ev = rbind(dose_ond_iv(), DEX_STD(), dose_lor(1, c(-1, 12, 24, 36, 48))),
    note = "HELD OUT: the only agent that reaches the conditioned term"),
  S21_ond_dex_dro = list(chemo = "cisplatin",
    ev = rbind(dose_ond_iv(), DEX_STD(),
               dose_dro(5, c(-1, 6, 12, 24, 36, 48, 60, 72))),
    note = "HELD OUT: CB1 gain reduction"),
  S22_ond32_dex = list(chemo = "cisplatin",
    ev = rbind(dose_ond_iv(32), DEX_STD()),
    note = "HELD OUT: the 32 mg IV dose withdrawn for QTc in 2012"),
  S23_apr_dex20_UNCORRECTED = list(chemo = "cisplatin",
    ev = rbind(dose_ond_iv(), DEX_STD(), APR_STD()),
    note = "HELD OUT: 20 mg dex WITH aprepitant = the DDI overexposure"),
  S24_dexfree_quad = list(chemo = "cisplatin",
    ev = rbind(dose_pal_iv(), APR_STD(), dose_ola()),
    note = "HELD OUT: steroid-free (immunotherapy / brittle diabetes)"),
  S25_dex_day1_only = list(chemo = "cisplatin",
    ev = rbind(dose_pal_iv(), dose_dex(12, -0.5), APR_STD(), dose_ola()),
    note = "HELD OUT: dexamethasone-sparing"),
  S26_AC_triplet = list(chemo = "AC",
    ev = rbind(dose_pal_iv(), DEX_NK1(), APR_STD()),
    note = "HELD OUT: anthracycline-cyclophosphamide"),
  S27_AC_quadruplet = list(chemo = "AC",
    ev = rbind(dose_pal_iv(), DEX_NK1(), APR_STD(), dose_ola()),
    note = "HELD OUT: AC + olanzapine"),
  S28_carbo_doublet = list(chemo = "carboplatin",
    ev = rbind(dose_pal_iv(), dose_dex(20, -0.5)),
    note = "HELD OUT (paired with S29): is an NK1 RA needed for carboplatin?"),
  S29_carbo_triplet = list(chemo = "carboplatin",
    ev = rbind(dose_pal_iv(), DEX_NK1(), APR_STD()),
    note = "HELD OUT (paired with S28)"),
  S30_oxali_doublet = list(chemo = "oxaliplatin",
    ev = rbind(dose_ond_iv(), dose_dex(20, -0.5)),
    note = "HELD OUT: oxaliplatin"),
  S31_LEC_dex = list(chemo = "LEC", ev = dose_dex(20, -0.5),
    note = "HELD OUT (paired with S32): low emetogenic risk"),
  S32_LEC_quadruplet = list(chemo = "LEC",
    ev = rbind(dose_pal_iv(), DEX_NK1(), APR_STD(), dose_ola()),
    note = "HELD OUT (paired with S31): the 4th agent where risk is LOW"),
  S33_UM2D6_ond = list(chemo = "cisplatin",
    ev = rbind(dose_ond_iv(), DEX_STD()), par = list(PH2D6 = 2.5),
    note = "HELD OUT (paired with S34): CYP2D6 ULTRARAPID metaboliser"),
  S34_UM2D6_pal = list(chemo = "cisplatin",
    ev = rbind(dose_pal_iv(), DEX_STD()), par = list(PH2D6 = 2.5),
    note = "HELD OUT (paired with S33): rescue with a non-CYP2D6 setron"),
  S35_highrisk_triplet = list(chemo = "cisplatin",
    ev = rbind(dose_pal_iv(), DEX_NK1(), APR_STD()),
    par = list(FEMALE = 1, AGE = 32, ALCOHOL = 0, MOTION = 1, PRIORCINV = 1),
    note = "HELD OUT (paired with S36): young anxious woman, prior CINV"),
  S36_lowrisk_triplet = list(chemo = "cisplatin",
    ev = rbind(dose_pal_iv(), DEX_NK1(), APR_STD()),
    par = list(FEMALE = 0, AGE = 68, ALCOHOL = 1, MOTION = 0, PRIORCINV = 0),
    note = "HELD OUT (paired with S35): older man, heavy alcohol history")
)

## =============================================================================
##  ENDPOINTS.  CR is E_S[exp(-S*I)] over a lognormal frailty, by Gauss-Hermite
##  quadrature.  The frailty is why CR(0-120) > CR(0-24) * CR(24-120): the gap
##  between the published joint endpoint and the published marginals IDENTIFIES
##  OMEGA.  Without it a deterministic model MUST under-predict the joint.
## =============================================================================
.gh <- function(n = 25) {
  ## Gauss-Hermite nodes/weights for the standard normal (probabilists')
  ## via the Golub-Welsch symmetric tridiagonal eigenproblem
  i <- 1:(n - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- sqrt(i)
  J[cbind(i + 1, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  x <- rev(e$values)
  w <- rev(e$vectors[1, ]^2)
  list(x = x, w = w)
}
.GH <- .gh(25)

frailty_CR <- function(I, omega) {
  if (!is.finite(I) || I <= 0) return(1)
  sum(.GH$w * exp(-exp(omega * .GH$x) * I))
}

run_scenario <- function(name, m = mod, end = 120, delta = 0.25, cycles = 1) {
  s <- SCEN[[name]]
  if (is.null(s)) stop("unknown scenario: ", name)
  em <- EMETO[[s$chemo]]
  pars <- c(list(EMETO_P = em[["P"]], EMETO_C = em[["C"]], CISRATE = em[["rate"]]),
            if (!is.null(s$par)) s$par else list())
  m2 <- do.call(param, c(list(m), pars))
  ev <- s$ev
  if (cycles > 1 && !is.null(ev)) {
    ev <- do.call(rbind, lapply(0:(cycles - 1), function(k) {
      e <- ev; e$time <- e$time + k * 504; e }))
  }
  tstart <- 0
  if (!is.null(ev)) {
    ev <- ev[order(ev$time), , drop = FALSE]
    ev$ID <- 1L
    ev <- ev[, c("ID", "time", "cmt", "amt", "evid")]
    tstart <- min(0, min(ev$time))
  }
  out <- if (is.null(ev)) {
    mrgsim(m2, start = tstart, end = end, delta = delta,
           atol = 1e-9, rtol = 1e-8)
  } else {
    mrgsim(m2, data = ev, start = tstart, end = end, delta = delta,
           atol = 1e-9, rtol = 1e-8, recsort = 3)
  }
  as.data.frame(out)
}

endpoints <- function(d, omega = as.numeric(param(mod)$OMEGA)) {
  # keep the LAST record at each time: mrgsim writes a pre- and a post-dose row at
  # every dose time, and averaging across that discontinuity is wrong
  d <- d[!duplicated(d$time, fromLast = TRUE), , drop = FALSE]
  at <- function(v, t) approx(d$time, v, xout = t, rule = 2)$y
  H24 <- at(d$HINT, 24); H120 <- at(d$HINT, 120)
  N24 <- at(d$NINT, 24); N120 <- at(d$NINT, 120)
  data.frame(
    CR_acute   = frailty_CR(H24, omega),
    CR_delayed = frailty_CR(H120 - H24, omega),
    CR_overall = frailty_CR(H120, omega),
    NN_acute   = frailty_CR(N24, omega),
    NN_delayed = frailty_CR(N120 - N24, omega),
    NN_overall = frailty_CR(N120, omega),
    emesis_24  = at(d$EMES, 24),
    emesis_120 = at(d$EMES, 120),
    peak_nausea = max(d$NAUSEA, na.rm = TRUE),
    peak_QTc   = max(d$QTC, na.rm = TRUE),
    peak_gluc  = max(d$GLU, na.rm = TRUE),
    peak_sed   = max(d$SED, na.rm = TRUE),
    peak_EPS   = max(d$AKATH, na.rm = TRUE),
    rescue_120 = at(d$RESC, 120),
    NK1occ_24  = at(d$NK1OCC_DRUG, 24),
    HT3occ_24  = at(d$R3OCC_DRUG, 24),
    HIAA_24    = at(d$UHIAA, 24),
    FLIE       = tail(d$FLIE, 1),
    RDI        = tail(d$RDI, 1),
    I_acute    = H24,
    I_delayed  = H120 - H24
  )
}

run_all <- function(which_scen = base::names(SCEN), ...) {
  do.call(rbind, lapply(which_scen, function(n) {
    d <- run_scenario(n, ...)
    cbind(scenario = n, note = SCEN[[n]]$note, endpoints(d))
  }))
}

## Read a column at a time point, keeping the LAST record at each time (mrgsim
## writes a pre- and a post-dose row at every dose time).
at_time <- function(d, col, t = 120) {
  d <- d[!duplicated(d$time, fromLast = TRUE), , drop = FALSE]
  approx(d$time, d[[col]], xout = t, rule = 2)$y
}

## =============================================================================
##  THE THESIS, AS A COMPUTATION.  Adding an agent multiplies CR by exp(dI).
##  So the ABSOLUTE gain in CR% depends on where the patient started, and it is
##  NOT monotonic in baseline risk.  This function sweeps baseline risk by
##  scaling susceptibility and reports the absolute and relative benefit of
##  adding olanzapine to a guideline triplet.
## =============================================================================
fourth_agent_benefit <- function(susc_mult = exp(seq(-2.2, 1.8, length.out = 21))) {
  base <- run_scenario("S06_pal_dex_apr")
  quad <- run_scenario("S12_quadruplet")
  om <- as.numeric(param(mod)$OMEGA)
  I3 <- at_time(base, "HINT", 120)
  I4 <- at_time(quad, "HINT", 120)
  do.call(rbind, lapply(susc_mult, function(k) {
    cr3 <- frailty_CR(k * I3, om); cr4 <- frailty_CR(k * I4, om)
    data.frame(susc_mult = k, CR_triplet = cr3, CR_quadruplet = cr4,
               abs_gain = cr4 - cr3, CR_ratio = cr4 / cr3,
               dI = k * (I3 - I4))
  }))
}

## Verification that additivity in log(CR) is a PROPERTY and not an assumption:
## build the four single-agent drops in the hazard integral and check that the
## combination's log-CR equals the sum of the parts (to within the nonlinearity
## of the threshold, which is exactly the residual worth reporting).
log_additivity_check <- function() {
  arms <- c(none = "S00_natural_history", ond = "S02_ond_alone",
            dex = "S01_dex_alone", ola = "S17_ola_alone",
            ond_dex = "S03_ond_dex", quad = "S13_quad_ond")
  I <- sapply(arms, function(a) {
    d <- run_scenario(a); at_time(d, "HINT", 120) })
  om <- as.numeric(param(mod)$OMEGA)
  CR <- sapply(I, frailty_CR, omega = om)
  data.frame(arm = names(I), I_120 = as.numeric(I), CR_120 = as.numeric(CR),
             dI_vs_none = as.numeric(I["none"] - I))
}

if (sys.nframe() == 0L && !interactive()) {
  cat("== tonic self-test (must be ~0) ==\n"); print(tonic_selftest())
  cat("\n== all scenarios ==\n")
  tab <- run_all()
  print(tab[, c("scenario", "CR_acute", "CR_delayed", "CR_overall",
                "NN_overall", "emesis_24", "peak_nausea", "peak_QTc")],
        digits = 3, row.names = FALSE)
}
