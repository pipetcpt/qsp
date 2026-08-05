## =============================================================================
##  sbe_mrgsolve_model.R
##  Snakebite envenoming + antivenom — 50-ODE QSP model
##  독사 교상(뱀독 중독) + 항독소 — 50-ODE 정량적 시스템 약리학 모델
## =============================================================================
##
##  THE ONE STRUCTURAL CLAIM
##  ------------------------
##      Antivenom BINDS.  It does not UNDO.
##
##  An antivenom antibody can only remove venom that is still free in plasma or
##  lymph.  It cannot pull a three-finger toxin off an acetylcholine receptor,
##  cannot un-cleave a fibrinogen molecule, and cannot un-destroy a motor nerve
##  terminal.  The model is therefore built with two structurally separate
##  clocks:
##
##      clock 1  (hours)  the VENOM clock    — antivenom acts here
##      clock 2  (days)   the SUBSTRATE clock — resynthesis / regeneration
##
##  and every clinical endpoint is read off clock 2, driven by the time-integral
##  of clock 1.  Four consequences follow arithmetically:
##
##   1  TIME-TO-ANTIVENOM SETS THE DEPTH OF THE NADIR, NOT THE SLOPE OF
##      RECOVERY.  Fibrinogen climbs back at the rate the liver can make it.
##      The hepatic synthesis term contains no antivenom variable at all.
##   2  RECURRENCE IS A RATIO OF TWO HALF-LIVES.  Venom keeps arriving from the
##      bite-site depot and the deep tissue compartment with a terminal half-life
##      of ~50 h.  Ovine Fab is cleared with a half-life of ~22 h.  When the
##      antivenom half-life is the shorter one, free venom MUST reappear, and no
##      larger front-loaded dose prevents it.
##   3  PRESYNAPTIC AND POSTSYNAPTIC PARALYSIS LOOK IDENTICAL AND ANSWER
##      OPPOSITE QUESTIONS.  Postsynaptic block is an occupancy that decays once
##      free toxin is gone and can be out-competed by raising acetylcholine.
##      Presynaptic PLA2 neurotoxicity is destruction: neither antivenom nor
##      neostigmine appears anywhere in the dTERM/dt equation.
##   4  ACUTE KIDNEY INJURY IS AN INTEGRAL, NOT A LEVEL.  Nephron loss is
##      driven by accumulated fibrin flux and pigment cast burden, so creatinine
##      peaks on day 5 — long after venom is undetectable.
##
##  UNITS
##  -----
##  Time      hours
##  Venom     mg of toxin protein, per toxin class
##  Antivenom mgNE = milligrams of venom-neutralising equivalents, which is the
##            unit antivenom is actually LABELLED in (Indian polyvalent ASV:
##            0.6 mg Daboia russelii venom neutralised per mL = 6 mgNE per
##            10-mL vial).  Working in the labelled unit makes antivenom potency
##            a regulatory fact rather than a fitted parameter.
##
##  PROVENANCE / VERIFICATION
##  -------------------------
##  The build environment for this repository has no R runtime, so every
##  equation below was independently re-implemented in Python
##  (`sbe_reference_model.py`) and actually integrated.  All numbers quoted in
##  README.md come from that run (`sbe_reference_output.txt`).  The two
##  implementations are term-for-term identical; the Python file carries inline
##  NOTE(defect n) comments at each of the places integration exposed a real
##  error in the first draft, and those same notes are reproduced here.
##
##  Requires: mrgsolve, dplyr, tidyr, ggplot2
## =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

# =============================================================================
#  MODEL CODE
# =============================================================================
sbe_code <- '
$PARAM @annotated
// ---------------------------------------------------------------- venom dose
VDOSE   :  63.0 : Injected venom protein mass (mg)
F_SVMP  :  0.28 : Mass fraction, snake venom metalloproteinase
F_SVSP  :  0.14 : Mass fraction, snake venom serine protease
F_PLA2  :  0.34 : Mass fraction, secreted phospholipase A2
F_TFTX  :  0.01 : Mass fraction, three-finger toxin
F_PRE   :  0.09 : Fraction of PLA2 mass that is presynaptically neurotoxic
F_MYO   :  0.50 : Relative systemic myotoxicity of this venoms PLA2
F_LOC   :  0.30 : Relative local necrotising potency
KOFFM   :  1.00 : Multiplier on alpha-neurotoxin dissociation rate
F_LEAK  :  1.00 : Relative capillary-leak potency

// ------------------------------------------------- venom class disposition
// Two parallel depots and ONE systemic compartment.  ka falls with molecular
// weight because absorption from a subcutaneous bite is lymphatic for the large
// enzymes and partly capillary for the small ones.
//
// NOTE(defect 6).  The first version gave venom a conventional two-compartment
// SYSTEMIC disposition with a computed terminal half-life of 49.5 h, and CLAIM 2
// WOULD NOT REPRODUCE: the bite-site depot was empty by 24 h, the deep
// compartment held almost nothing, and no antivenom fragment — Fab included —
// could run out while venom was still arriving.  The defect was structural, not
// numerical.  Prolonged venom antigenaemia after viperid bites is an ABSORPTION
// phenomenon (sequestration in a deep depot that releases over days), not a
// distribution phenomenon, so the peripheral compartment was rebuilt as a SLOW
// DEPOT.  The terminal half-life of the venom is now flip-flop-limited:
//   SVMP 44.7 h · SVSP 36.5 h · PLA2 19.3 h · 3FTx 6.0 h
// which is 5.5x longer than its systemic elimination half-life, and it is the
// number an antivenom fragment actually has to outlast.
KA_SVMP :  0.250 : Fast-depot absorption rate constant, SVMP (1/h)
KA_SVSP :  0.350 : Fast-depot absorption rate constant, SVSP (1/h)
KA_PLA2 :  0.500 : Fast-depot absorption rate constant, PLA2 (1/h)
KA_TFTX :  1.200 : Fast-depot absorption rate constant, 3FTx (1/h)
KAS_SVMP: 0.0120 : SLOW-depot release rate constant, SVMP (1/h)
KAS_SVSP: 0.0150 : SLOW-depot release rate constant, SVSP (1/h)
KAS_PLA2: 0.0300 : SLOW-depot release rate constant, PLA2 (1/h)
KAS_TFTX: 0.1000 : SLOW-depot release rate constant, 3FTx (1/h)
FSL_SVMP:  0.350 : Fraction of the SVMP mass sequestered in the slow depot
FSL_SVSP:  0.320 : Fraction of the SVSP mass sequestered in the slow depot
FSL_PLA2:  0.150 : Fraction of the PLA2 mass sequestered in the slow depot
FSL_TFTX:  0.050 : Fraction of the 3FTx mass sequestered in the slow depot
KDD_SVMP:  0.010 : Fast-depot local proteolysis, SVMP (1/h)
KDD_SVSP:  0.012 : Fast-depot local proteolysis, SVSP (1/h)
KDD_PLA2:  0.015 : Fast-depot local proteolysis, PLA2 (1/h)
KDD_TFTX:  0.020 : Fast-depot local proteolysis, 3FTx (1/h)
KDS_SVMP: 0.0035 : Slow-depot local proteolysis, SVMP (1/h)
KDS_SVSP: 0.0040 : Slow-depot local proteolysis, SVSP (1/h)
KDS_PLA2: 0.0060 : Slow-depot local proteolysis, PLA2 (1/h)
KDS_TFTX: 0.0150 : Slow-depot local proteolysis, 3FTx (1/h)
VC_SVMP :   3.50 : Central volume, SVMP (L)
VC_SVSP :   3.50 : Central volume, SVSP (L)
VC_PLA2 :   4.00 : Central volume, PLA2 (L)
VC_TFTX :   8.00 : Central volume, 3FTx (L)
CL_SVMP :   0.30 : Systemic clearance, SVMP (L/h)
CL_SVSP :   0.40 : Systemic clearance, SVSP (L/h)
CL_PLA2 :   0.90 : Systemic clearance, PLA2 (L/h)
CL_TFTX :   3.00 : Systemic clearance, 3FTx (L/h)

// eps = mg of toxin neutralised per mgNE of antivenom.  Antivenoms are raised
// against whole venom, whose immunodominant components are the LARGE enzymes,
// so per unit of labelled neutralising capacity they cover SVMP and SVSP well
// and the small, weakly immunogenic PLA2 and three-finger toxins badly.  This
// single vector is why a neurotoxic elapid bite needs 2-3x the vial count of a
// viperid bite carrying MORE venom protein.
// It is a MODEL PARAMETER, not a measured potency — see sbe_references.md.
EPS_SVMP:   1.20 : Relative neutralising potency, SVMP (mg toxin per mgNE)
EPS_SVSP:   1.00 : Relative neutralising potency, SVSP
EPS_PLA2:   0.60 : Relative neutralising potency, PLA2
EPS_TFTX:   0.35 : Relative neutralising potency, 3FTx

// -------------------------------------------------------- antivenom product
// Terminal half-lives implied by these values (computed, not asserted):
//   Fab      V1 4.0 V2 12.0 Q 2.00 CL 0.600  ->  t1/2 beta  21.8 h
//   F(ab)2   V1 3.5 V2  4.5 Q 0.50 CL 0.045  ->  t1/2 beta 126.8 h
//   IgG      V1 3.2 V2  2.6 Q 0.30 CL 0.030  ->  t1/2 beta 136.8 h
// matching the published 12-23 h (ovine Fab), ~130 h (equine F(ab)2) and
// 90-200 h (equine whole IgG).
V1A     :   3.2 : Antivenom central volume (L)
V2A     :   2.6 : Antivenom peripheral volume (L)
QA      :  0.30 : Antivenom intercompartmental clearance (L/h)
CLA     : 0.030 : Antivenom clearance (L/h)
KTIS    : 0.012 : Antivenom transfer into bite-site tissue (1/h)
RHO     : 1.300 : Product reactogenicity (anaphylactoid index per mgNE/min)
ICF     :  1.00 : Relative immune-complex burden per mgNE of protein

// ---------------------------------------------------- antivenom infusion(s)
AVT0    :   4.0 : Time of first antivenom infusion (h post-bite)
AVAMT   :  60.0 : First antivenom dose (mgNE)
AVDUR   :   1.0 : Infusion duration (h)
RPT_T   : 1e6   : Time of a single ad-hoc repeat dose (h)
RPT_AMT :   0.0 : Ad-hoc repeat dose (mgNE)
MSTART  : 1e6   : Start of maintenance schedule (h)
MINT    :   6.0 : Maintenance dosing interval (h)
MAMT    :   0.0 : Maintenance dose per administration (mgNE)
MN      :   0.0 : Number of maintenance doses
FLS     : 1e6   : Crystalloid infusion start (h)
FLE     : 1e6   : Crystalloid infusion stop (h)
FLR     :  0.00 : Crystalloid infusion rate (L/h)

// ---------------------------------------------------- binding / complexes
KB0     :  12.0 : Plasma antivenom-toxin association rate (L/mgNE/h)
KB_LOC  :  0.45 : Bite-site antivenom-toxin association rate (L/mgNE/h)
KEL_CPLX: 0.060 : Reticuloendothelial clearance of complex (1/h)
// Antibody-antigen complexes do dissociate.  The DEFAULT is ZERO, so that every
// recurrence result rests on absorption alone and nothing is smuggled in; the
// sensitivity analysis sweeps this to ask how fast dissociation would have to be
// to explain the clinical phenomenon that absorption cannot.
KOFF_AV : 0.0   : Dissociation of the antivenom-toxin complex (1/h)
CV_LOD  : 0.0020: Limit of detection of a venom antigen EIA (mg/L, ~2 ng/mL)
KOUT_T  : 0.030 : Antivenom loss from bite-site tissue (1/h)
VLOC    :  0.15 : Bite-site tissue distribution volume (L)

// ------------------------------------------------------------- haemostasis
FG0     :  2.80 : Baseline fibrinogen (g/L)
KDEG_FG :0.00722: Fibrinogen turnover rate, t1/2 96 h (1/h)
KM_FG   :  4.00 : Michaelis constant of venom fibrinogenolysis (g/L)
KCAT_SVSP: 2.50 : Thrombin-like fibrinogen consumption ((g/L/h) per (mg/L))
KCAT_SVMP: 0.30 : SVMP fibrinogen consumption ((g/L/h) per (mg/L))
F_PLASMIN: 0.05 : Plasmin-mediated fraction of fibrinogen loss
KUP_SYN : 0.060 : Acute-phase upregulation rate of hepatic synthesis (1/h)
KDOWN_SYN:0.020 : Decay of acute-phase upregulation (1/h)
SYNUP_MAX: 3.20 : Ceiling on hepatic fibrinogen synthesis upregulation
FX0     :  1.00 : Baseline factor X activity (fraction)
K_FX_REC:0.01733: Factor X resynthesis, t1/2 40 h (1/h)
KFX     :  0.35 : FX activation/consumption by SVMP (1/h per mg/L)
PLT0    : 250.0 : Baseline platelet count (1e9/L)
KDEG_PLT:0.005775: Platelet turnover, t1/2 120 h (1/h)
KPLT_BOOST: 1.50: Marrow response gain
KPLT_SVMP: 0.045: Platelet consumption by SVMP (1/h per mg/L)
KPLT_SVSP: 0.070: Platelet consumption by SVSP (1/h per mg/L)
XDP0    :  0.25 : Baseline D-dimer (ug/mL)
KEL_XDP : 0.0866: D-dimer elimination, t1/2 8 h (1/h)
G_XDP   :   9.0 : D-dimer generated per (g/L/h) of fibrin turnover

// --------------------------------------------------- neuromuscular junction
SF0     :  4.00 : Normal neuromuscular safety factor
KOFF_R  : 0.020 : alpha-neurotoxin dissociation from nAChR, t1/2 35 h (1/h)
KD_R    : 0.060 : Apparent perijunctional Kd of alpha-neurotoxin (mg/L)
KDES_TERM: 0.550: Presynaptic terminal destruction (1/h per mg/L PLA2)
KRGN_TERM:0.00578: Motor terminal regeneration, t1/2 5 d (1/h)
EMAX_NEO:  1.20 : Maximal neostigmine effect on cleft acetylcholine
EC50_NEO:0.0050 : Neostigmine EC50 (mg/L)
// NOTE(defect 5).  Neostigmine was first written as a free-standing multiplier
// on the safety factor: SF = SF0*TERM*(1-BR)*neo.  That made an
// anticholinesterase rescue a PRESYNAPTIC krait bite as efficiently as a
// postsynaptic cobra bite (ventilation 28 h -> 4 h), which is the opposite of
// the clinical fact and destroyed claim 3.  The error was MECHANISTIC, not
// numerical: raising cleft acetylcholine works by OUT-COMPETING a competitive
// antagonist at the receptor, so it belongs inside the occupancy term, where a
// small BR leaves it almost nothing to do.  W_PRE_NEO retains the small, real,
// non-decisive effect of prolonging the dwell time of whatever transmitter a
// damaged terminal still releases.
W_PRE_NEO: 0.050: Weight of the neostigmine effect on a damaged terminal

// ------------------------------------------------------------ local tissue
// NOTE(defect 1).  The bite site is 0.15 L, so local toxin concentrations are
// two orders of magnitude above plasma (a Bothrops bite starts near 167 mg/L
// of SVMP locally).  The first calibration set KNEC by analogy with the plasma
// rate constants and drove every venomous species to >50% limb necrosis inside
// an hour.  KNEC must be ~40x SMALLER than a plasma rate constant precisely
// BECAUSE the local compartment is small: the concentration does the work.
KNEC    :6.50e-4: Myonecrosis rate (1/h per mg/L local myotoxin)
A_SVMP_NEC: 1.00: Relative necrotising weight, SVMP
A_PLA2_NEC: 1.60: Relative necrotising weight, myotoxic PLA2
A_TFTX_NEC: 2.50: Relative necrotising weight, cytotoxic 3FTx
KHEAL_NEC:0.00150: Muscle healing, t1/2 19 d (1/h)
KEDEMA  :6.0e-4 : Local oedema formation (L/h per mg/L)
KDRAIN  : 0.0200: Oedema drainage (1/h)
CPB     :   8.0 : Baseline compartment pressure (mmHg)
CP_MAX  :  58.0 : Maximal compartment pressure increment (mmHg)
CP_K50  :  0.35 : Oedema volume at half-maximal pressure (L)
// NOTE(defect 2).  CK release proportional to the necrosis FLUX alone put the
// CK peak at 4-5 h and at 2e5 U/L.  Real CK peaks at 12-48 h: an enzyme leaves
// a necrotic fibre for as long as the fibre stays necrotic, so release must be
// driven by the necrosis STOCK as well as its flux.
G_CK_LOC: 1.60e4: CK release per unit local necrosis flux (U/L per 1/h)
G_CK_STOCK: 380.0: CK release per unit standing necrosis (U/L/h)
G_CK_SYS: 7.0e3 : CK release per unit systemic myotoxicity (U/L per 1/h)
KEL_CK  :0.01925: CK elimination, t1/2 36 h (1/h)
KMYO_SYS: 0.055 : Systemic myotoxicity (1/h per mg/L PLA2)
// NOTE(defect 9, found by the baseline-stability check).  Myoglobin had a
// release term and two clearance terms but NO basal production, so plasma
// myoglobin decayed from 0.030 to 0 over 14 days in a patient with no
// envenoming at all: the model baseline was not a fixed point, the classic
// signature of a missing zero-order production term.
G_MB    :  95.0 : Myoglobin release coefficient (ug/mL per 1/h)
MB0     : 0.030 : Baseline plasma myoglobin (ug/mL)
KEL_MB_REN: 0.260: Renal myoglobin clearance at full nephron mass (1/h)
KEL_MB_NONREN: 0.045: Non-renal myoglobin clearance (1/h)

// ------------------------------------------------------------------ kidney
KCAST   : 0.0180: Pigment cast formation (1/h per ug/mL myoglobin)
KCLR_CAST:0.0350: Cast clearance (1/h)
// NOTE(defect 3).  KN_FIB was first sized as if FGLOSS were dimensionless.  Its
// time integral is the total GRAMS PER LITRE of fibrinogen consumed (~5 g/L in
// an untreated Russell viper bite), so a rate constant of 0.013 contributed
// 0.05 of a nephron log-loss, i.e. nothing.  The intravascular-coagulation arm
// of venom AKI was silently ABSENT from the first version of the model.
KN_FIB  : 0.2000: Nephron loss per (g/L) of fibrinogen consumed
KN_CAST :0.00300: Nephron loss per unit cast burden (1/h)
KN_ISCH : 0.0900: Nephron loss from renal hypoperfusion (1/h)
KN_DIR  :0.00800: Direct venom nephrotoxicity (1/h per mg/L SVMP)
KRECOV_NEPH:0.00800: Nephron repair (1/h)
KFIBROSIS: 0.100: Fraction of nephron loss that becomes permanent scar
SCR0    :  0.90 : Baseline serum creatinine (mg/dL)
KEL_CR  : 0.060 : Creatinine elimination rate constant at full GFR (1/h)

// ----------------------------------------------------------- haemodynamics
// NOTE(defect 4).  The first leak calibration took plasma volume to 0.9 L and
// MAP to 20 mmHg in every viperid scenario, killing 100% of untreated patients
// from shock before coagulopathy or AKI could be read at all.  The
// transcapillary refill term was an order of magnitude too weak.
PV0     :  3.00 : Baseline plasma volume (L)
KLEAK   : 0.0300: Transcapillary plasma loss (1/h)
KREFILL : 0.0600: Transcapillary refill (1/h)
KGLX    : 0.0300: Glycocalyx shedding (1/h per mg/L)
KGLXR   : 0.0200: Glycocalyx repair (1/h)
MAP0    :  92.0 : Baseline mean arterial pressure (mmHg)
MAP_EXP :  1.60 : Exponent of the volume-pressure relation
MAP_ANAPHYL: 0.38: Maximal fractional MAP fall from anaphylaxis

// ------------------------------------------------------------ inflammation
G_IL6   : 170.0 : IL-6 generation coefficient (pg/mL/h per mg/L)
KEL_IL6 : 0.347 : IL-6 elimination, t1/2 2 h (1/h)
G_TNF   :  55.0 : TNF-alpha generation coefficient (pg/mL/h per mg/L)
KEL_TNF : 0.990 : TNF-alpha elimination (1/h)

// -------------------------------------------------- antivenom adverse events
KEL_MCA : 0.800 : Decay of anaphylactoid activation (1/h)
MCA_ANAPHYL: 1.00: Anaphylactoid index defining a clinical reaction
K_IC    :0.00090: Immune-complex formation (1/h per mgNE)
KEL_IC  : 0.0200: Immune-complex clearance (1/h)
K_SS    :0.01800: Serum sickness generation (1/h)
KEL_SS  : 0.0120: Serum sickness resolution (1/h)

// ----------------------------------------------------------------- hazards
H_BLEED :0.00830: Systemic bleeding hazard with unclottable blood (1/h)
FG_UNCLOT: 0.50 : Fibrinogen at which the 20WBCT turns positive (g/L)
PLT_LOW :  50.0 : Platelet count defining the high-risk term (1e9/L)
HD_BLEED:0.00180: Death hazard from haemorrhage (1/h)
HD_RESP_NOVENT:0.09000: Death hazard, respiratory failure, no ventilator (1/h)
HD_RESP_VENT:0.000400: Death hazard, respiratory failure, ventilated (1/h)
HD_AKI_NODIAL:0.01000: Death hazard from AKI without dialysis (1/h)
HD_AKI_DIAL:0.000600: Death hazard from AKI with dialysis (1/h)
NEPH_HAZ:  0.500: Nephron fraction below which the AKI hazard is graded up
HD_SHOCK:0.00600: Death hazard from shock (1/h)
HD_ANAPHYL:0.00800: Death hazard from anaphylaxis (1/h)
ICU     :   1.0 : Mechanical ventilation available (1 = yes)
DIALYSIS:   1.0 : Renal replacement therapy available (1 = yes)

// -------------------------------------------------------- adjunct drug PK
NEO_KA  : 1.400 : Neostigmine IM absorption (1/h)
NEO_KEL : 0.835 : Neostigmine elimination, t1/2 50 min (1/h)
NEO_V   :  50.0 : Neostigmine volume of distribution (L)
VAR_KA  : 1.500 : Varespladib-methyl absorption (1/h)
VAR_KEL : 0.231 : Varespladib elimination, t1/2 3 h (1/h)
VAR_V   :  12.0 : Varespladib volume of distribution (L)
VAR_F   : 0.500 : Varespladib oral bioavailability
VAR_IC50: 0.0100: Varespladib IC50 against secreted PLA2 (mg/L)
TXA_KEL : 0.347 : Tranexamic acid elimination, t1/2 2 h (1/h)
TXA_V   :  12.0 : Tranexamic acid volume of distribution (L)
TXA_IC50:  5.00 : Tranexamic acid IC50 against plasmin (mg/L)

$CMT @annotated
D_SVMP  : Bite-site depot, SVMP (mg)
D_SVSP  : Bite-site depot, SVSP (mg)
D_PLA2  : Bite-site depot, PLA2 (mg)
D_TFTX  : Bite-site depot, 3FTx (mg)
V_SVMP  : Central free SVMP (mg)
V_SVSP  : Central free SVSP (mg)
V_PLA2  : Central free PLA2 (mg)
V_TFTX  : Central free 3FTx (mg)
P_SVMP  : SLOW / sequestered bite-site depot, SVMP (mg)
P_SVSP  : SLOW / sequestered bite-site depot, SVSP (mg)
P_PLA2  : SLOW / sequestered bite-site depot, PLA2 (mg)
P_TFTX  : Deep tissue 3FTx (mg)
C_SVMP  : Antivenom-SVMP complex (mg toxin)
C_SVSP  : Antivenom-SVSP complex (mg toxin)
C_PLA2  : Antivenom-PLA2 complex (mg toxin)
C_TFTX  : Antivenom-3FTx complex (mg toxin)
A_C     : Antivenom, central (mgNE)
A_P     : Antivenom, peripheral (mgNE)
A_T     : Antivenom, bite-site tissue (mgNE)
NEO_A   : Neostigmine IM depot (mg)
NEO_C   : Neostigmine central (mg)
VAR_A   : Varespladib gut depot (mg)
VAR_C   : Varespladib central (mg)
TXA_C   : Tranexamic acid central (mg)
FG      : Fibrinogen (g/L)
FX      : Factor X activity (fraction)
PLT     : Platelet count (1e9/L)
XDP     : D-dimer / FDP (ug/mL)
SYNUP   : Hepatic fibrinogen synthesis upregulation (x normal)
BR      : nAChR fractional occupancy by 3FTx
TERM    : Presynaptic motor terminal integrity (0-1)
NEC     : Local myonecrosis fraction (0-1)
EDEMA   : Local oedema volume (L)
CK      : Creatine kinase (U/L)
MB      : Plasma myoglobin (ug/mL)
NEPH    : Functional nephron fraction (0-1)
SCR     : Serum creatinine (mg/dL)
CAST    : Tubular pigment cast burden (AU)
FIBR    : Permanent renal interstitial fibrosis (0-1)
PV      : Plasma volume (L)
GLX     : Endothelial glycocalyx integrity (0-1)
IL6     : Interleukin-6 (pg/mL)
TNF     : TNF-alpha (pg/mL)
MCA     : Mast-cell / anaphylactoid activation index (AU)
IC      : Circulating immune complex burden (AU)
SS      : Serum sickness score (AU)
HBLD    : Cumulative systemic-bleeding hazard
HDTH    : Cumulative death hazard
VAUC    : Free venom AUC (mg.h/L)
AVCUM   : Cumulative antivenom administered (mgNE)

$MAIN
// ---- the bite: partition the injected mass across the four toxin classes ----
// the injected mass splits at t = 0 between a fast, superficial depot and a
// slow, sequestered one.  The slow one is the terminal phase of the venom.
D_SVMP_0 = VDOSE * F_SVMP * (1.0 - FSL_SVMP);
D_SVSP_0 = VDOSE * F_SVSP * (1.0 - FSL_SVSP);
D_PLA2_0 = VDOSE * F_PLA2 * (1.0 - FSL_PLA2);
D_TFTX_0 = VDOSE * F_TFTX * (1.0 - FSL_TFTX);
P_SVMP_0 = VDOSE * F_SVMP * FSL_SVMP;
P_SVSP_0 = VDOSE * F_SVSP * FSL_SVSP;
P_PLA2_0 = VDOSE * F_PLA2 * FSL_PLA2;
P_TFTX_0 = VDOSE * F_TFTX * FSL_TFTX;

FG_0    = FG0;
FX_0    = FX0;
PLT_0   = PLT0;
XDP_0   = XDP0;
SYNUP_0 = 1.0;
BR_0    = 0.0;
TERM_0  = 1.0;
NEC_0   = 0.0;
EDEMA_0 = 0.0;
CK_0    = 90.0;
MB_0    = MB0;
NEPH_0  = 1.0;
SCR_0   = SCR0;
// the steady state of the pigment-cast compartment at baseline myoglobin and
// full nephron mass; starting it at zero left CAST drifting for 14 days
CAST_0  = KCAST*MB0*0.30/KCLR_CAST;
FIBR_0  = 0.0;
PV_0    = PV0;
GLX_0   = 1.0;
IL6_0   = 3.0;
TNF_0   = 4.0;

$ODE
// =========================================================================
//  infusion schedules (parameter driven, so that the anaphylactoid term can
//  see the infusion RATE — reactions are rate driven, not dose driven)
// =========================================================================
double AVRATE = 0.0;
if(SOLVERTIME >= AVT0 && SOLVERTIME < AVT0 + AVDUR)     AVRATE += AVAMT / AVDUR;
if(SOLVERTIME >= RPT_T && SOLVERTIME < RPT_T + AVDUR)   AVRATE += RPT_AMT / AVDUR;
if(MN > 0.5 && SOLVERTIME >= MSTART && SOLVERTIME < MSTART + MN * MINT) {
  double tm = SOLVERTIME - MSTART;
  double ph = tm - MINT * floor(tm / MINT);
  if(ph < AVDUR) AVRATE += MAMT / AVDUR;
}
double FLRATE = (SOLVERTIME >= FLS && SOLVERTIME < FLE) ? FLR : 0.0;

// =========================================================================
//  derived concentrations
// =========================================================================
double cSVMP = (V_SVMP > 0 ? V_SVMP : 0) / VC_SVMP;
double cSVSP = (V_SVSP > 0 ? V_SVSP : 0) / VC_SVSP;
double cPLA2 = (V_PLA2 > 0 ? V_PLA2 : 0) / VC_PLA2;
double cTFTX = (V_TFTX > 0 ? V_TFTX : 0) / VC_TFTX;
double cVfree = cSVMP + cSVSP + cPLA2 + cTFTX;

double ed   = (EDEMA > 0 ? EDEMA : 0);
double vloc = VLOC + ed;
double lSVMP = (D_SVMP > 0 ? D_SVMP : 0) / vloc;
double lSVSP = (D_SVSP > 0 ? D_SVSP : 0) / vloc;
double lPLA2 = (D_PLA2 > 0 ? D_PLA2 : 0) / vloc;
double lTFTX = (D_TFTX > 0 ? D_TFTX : 0) / vloc;
double cAt   = (A_T > 0 ? A_T : 0) / vloc;
double cA    = (A_C > 0 ? A_C : 0) / V1A;

double cNEO = (NEO_C > 0 ? NEO_C : 0) / NEO_V;
double cVAR = (VAR_C > 0 ? VAR_C : 0) / VAR_V;
double cTXA = (TXA_C > 0 ? TXA_C : 0) / TXA_V;

// varespladib is a catalytic-site PLA2 inhibitor with a 12 L volume of
// distribution: it reaches the bite site and the nerve terminal, where an IgG
// cannot.  INH_PLA2 multiplies EVERY PLA2-driven term in the model.
double INH_PLA2    = 1.0 / (1.0 + cVAR / VAR_IC50);
double INH_PLASMIN = 1.0 / (1.0 + cTXA / TXA_IC50);

// =========================================================================
//  VENOM PK + ANTIVENOM BINDING
//  Binding in plasma runs at ~150/h (a 15-second half-life) so the system is
//  effectively STOICHIOMETRIC: free venom is ~0 while a molar excess of
//  antivenom exists and returns the instant that excess is gone.  That is the
//  whole of claim 2 and it is emergent, not switched.
// =========================================================================
double fg_s = (FG > 1e-9 ? FG : 1e-9);

// local neutralisation.  KB_LOC/KB0 = 0.038: a 27-fold penalty, because an
// antibody penetrates oedematous, poorly-perfused, necrotic tissue badly.
// This is why intravenous antivenom does not prevent local necrosis.
// The fast depot is the superficial, well-drained tissue; the slow depot is the
// deep intramuscular and fascial sequestration that releases over DAYS.  Local
// neutralisation gets a 27-fold penalty in the fast depot (KB_LOC/KB0 = 0.038)
// and a further 4-fold penalty in the slow one — which is exactly why draining
// the depot is not an available way to stop recurrence.
double bl_SVMP = KB_LOC * cAt * (D_SVMP > 0 ? D_SVMP : 0);
double bl_SVSP = KB_LOC * cAt * (D_SVSP > 0 ? D_SVSP : 0);
double bl_PLA2 = KB_LOC * cAt * (D_PLA2 > 0 ? D_PLA2 : 0);
double bl_TFTX = KB_LOC * cAt * (D_TFTX > 0 ? D_TFTX : 0);
double bs2_SVMP = 0.25 * KB_LOC * cAt * (P_SVMP > 0 ? P_SVMP : 0);
double bs2_SVSP = 0.25 * KB_LOC * cAt * (P_SVSP > 0 ? P_SVSP : 0);
double bs2_PLA2 = 0.25 * KB_LOC * cAt * (P_PLA2 > 0 ? P_PLA2 : 0);
double bs2_TFTX = 0.25 * KB_LOC * cAt * (P_TFTX > 0 ? P_TFTX : 0);

double bs_SVMP = KB0 * EPS_SVMP * cA * (V_SVMP > 0 ? V_SVMP : 0);
double bs_SVSP = KB0 * EPS_SVSP * cA * (V_SVSP > 0 ? V_SVSP : 0);
double bs_PLA2 = KB0 * EPS_PLA2 * cA * (V_PLA2 > 0 ? V_PLA2 : 0);
double bs_TFTX = KB0 * EPS_TFTX * cA * (V_TFTX > 0 ? V_TFTX : 0);

// complex dissociation returns both toxin and neutralising capacity
double ds_SVMP = KOFF_AV * C_SVMP;
double ds_SVSP = KOFF_AV * C_SVSP;
double ds_PLA2 = KOFF_AV * C_PLA2;
double ds_TFTX = KOFF_AV * C_TFTX;

// mgNE consumed per mg of toxin bound = 1/eps
double bind_sys = (bs_SVMP - ds_SVMP)/EPS_SVMP + (bs_SVSP - ds_SVSP)/EPS_SVSP
                + (bs_PLA2 - ds_PLA2)/EPS_PLA2 + (bs_TFTX - ds_TFTX)/EPS_TFTX;
double bind_loc = (bl_SVMP + bs2_SVMP)/EPS_SVMP + (bl_SVSP + bs2_SVSP)/EPS_SVSP
                + (bl_PLA2 + bs2_PLA2)/EPS_PLA2 + (bl_TFTX + bs2_TFTX)/EPS_TFTX;

dxdt_D_SVMP = -KA_SVMP*D_SVMP - KDD_SVMP*D_SVMP - bl_SVMP;
dxdt_D_SVSP = -KA_SVSP*D_SVSP - KDD_SVSP*D_SVSP - bl_SVSP;
dxdt_D_PLA2 = -KA_PLA2*D_PLA2 - KDD_PLA2*D_PLA2 - bl_PLA2;
dxdt_D_TFTX = -KA_TFTX*D_TFTX - KDD_TFTX*D_TFTX - bl_TFTX;

dxdt_P_SVMP = -KAS_SVMP*P_SVMP - KDS_SVMP*P_SVMP - bs2_SVMP;
dxdt_P_SVSP = -KAS_SVSP*P_SVSP - KDS_SVSP*P_SVSP - bs2_SVSP;
dxdt_P_PLA2 = -KAS_PLA2*P_PLA2 - KDS_PLA2*P_PLA2 - bs2_PLA2;
dxdt_P_TFTX = -KAS_TFTX*P_TFTX - KDS_TFTX*P_TFTX - bs2_TFTX;

dxdt_V_SVMP = KA_SVMP*D_SVMP + KAS_SVMP*P_SVMP - CL_SVMP*V_SVMP/VC_SVMP
              - bs_SVMP + ds_SVMP;
dxdt_V_SVSP = KA_SVSP*D_SVSP + KAS_SVSP*P_SVSP - CL_SVSP*V_SVSP/VC_SVSP
              - bs_SVSP + ds_SVSP;
dxdt_V_PLA2 = KA_PLA2*D_PLA2 + KAS_PLA2*P_PLA2 - CL_PLA2*V_PLA2/VC_PLA2
              - bs_PLA2 + ds_PLA2;
dxdt_V_TFTX = KA_TFTX*D_TFTX + KAS_TFTX*P_TFTX - CL_TFTX*V_TFTX/VC_TFTX
              - bs_TFTX + ds_TFTX;

dxdt_C_SVMP = bs_SVMP + bl_SVMP + bs2_SVMP - KEL_CPLX*C_SVMP - ds_SVMP;
dxdt_C_SVSP = bs_SVSP + bl_SVSP + bs2_SVSP - KEL_CPLX*C_SVSP - ds_SVSP;
dxdt_C_PLA2 = bs_PLA2 + bl_PLA2 + bs2_PLA2 - KEL_CPLX*C_PLA2 - ds_PLA2;
dxdt_C_TFTX = bs_TFTX + bl_TFTX + bs2_TFTX - KEL_CPLX*C_TFTX - ds_TFTX;

// =========================================================================
//  ANTIVENOM PK
// =========================================================================
double tis = KTIS * (A_C > 0 ? A_C : 0) * (1.0 + 1.5*ed);
dxdt_A_C = AVRATE - CLA*A_C/V1A - QA*A_C/V1A + QA*A_P/V2A - tis - bind_sys;
dxdt_A_P = QA*A_C/V1A - QA*A_P/V2A;
dxdt_A_T = tis - KOUT_T*A_T - bind_loc;
dxdt_AVCUM = AVRATE;

// =========================================================================
//  ADJUNCT DRUG PK
// =========================================================================
dxdt_NEO_A = -NEO_KA*NEO_A;
dxdt_NEO_C =  NEO_KA*NEO_A - NEO_KEL*NEO_C;
dxdt_VAR_A = -VAR_KA*VAR_A;
dxdt_VAR_C =  VAR_F*VAR_KA*VAR_A - VAR_KEL*VAR_C;
dxdt_TXA_C = -TXA_KEL*TXA_C;

// =========================================================================
//  HAEMOSTASIS — clock 2, the substrate clock
//  Fibrinogen falls at the rate the venom enzyme works and rises at the rate
//  the LIVER works.  Those two rates have nothing to do with each other, and
//  antivenom can only touch the first.  Claim 1 is this pair of lines.
// =========================================================================
double mm  = fg_s / (KM_FG + fg_s);
double enz = KCAT_SVSP*cSVSP + KCAT_SVMP*cSVMP;
// 95% of the loss is direct enzymatic cleavage, untouchable by tranexamic acid;
// 5% is secondary plasmin-mediated fibrinogenolysis.
double FGLOSS = enz * mm * ((1.0 - F_PLASMIN) + F_PLASMIN*INH_PLASMIN);

dxdt_FG    = KDEG_FG*FG0*SYNUP - KDEG_FG*FG - FGLOSS;
double defc = 1.0 - FG/FG0;  if(defc < 0) defc = 0.0;
dxdt_SYNUP = KUP_SYN*defc*(SYNUP_MAX - SYNUP) - KDOWN_SYN*(SYNUP - 1.0);
dxdt_FX    = K_FX_REC*(FX0 - FX) - KFX*cSVMP*FX;
double pdef = 1.0 - PLT/PLT0;  if(pdef < 0) pdef = 0.0;
dxdt_PLT   = PLT0*KDEG_PLT*(1.0 + KPLT_BOOST*pdef) - KDEG_PLT*PLT
             - (KPLT_SVMP*cSVMP + KPLT_SVSP*cSVSP)*PLT;
dxdt_XDP   = G_XDP*FGLOSS - KEL_XDP*(XDP - XDP0);

// =========================================================================
//  NEUROMUSCULAR JUNCTION — the two mechanisms that look the same
// =========================================================================
// POSTSYNAPTIC: an occupancy.  Remove free toxin and it decays; raise
// acetylcholine and you out-compete what is left.  Reversible.
double konR = KOFF_R / KD_R;
double freeR = 1.0 - BR;  if(freeR < 0) freeR = 0.0;
dxdt_BR = konR*cTFTX*freeR - KOFF_R*KOFFM*BR;

// PRESYNAPTIC: destruction.  The only recovery term is regeneration of the
// terminal on a 5-day half-life.  NEITHER antivenom NOR neostigmine appears
// anywhere in this equation — that is the entirety of claim 3.
dxdt_TERM = -KDES_TERM*F_PRE*cPLA2*INH_PLA2*TERM + KRGN_TERM*(1.0 - TERM);

// =========================================================================
//  LOCAL TISSUE
//  The three-finger-toxin term is what makes an elapid bite necrotic despite
//  its carrying almost no metalloproteinase.
// =========================================================================
double freeN = 1.0 - NEC;  if(freeN < 0) freeN = 0.0;
double NECFLUX = KNEC * (A_SVMP_NEC*lSVMP + A_PLA2_NEC*lPLA2*INH_PLA2
                         + A_TFTX_NEC*lTFTX) * F_LOC * freeN;
double MYOSYS  = KMYO_SYS * cPLA2 * F_MYO * INH_PLA2;

dxdt_NEC   = NECFLUX - KHEAL_NEC*NEC;
dxdt_EDEMA = KEDEMA*(lSVMP + 0.5*lPLA2*INH_PLA2)*F_LOC - KDRAIN*EDEMA;
dxdt_CK    = G_CK_LOC*NECFLUX + G_CK_STOCK*NEC + G_CK_SYS*MYOSYS
             - KEL_CK*(CK - 90.0);
dxdt_MB    = MB0*(KEL_MB_REN + KEL_MB_NONREN)
             + G_MB*(NECFLUX*0.35 + NEC*0.010 + MYOSYS)
             - KEL_MB_REN*(NEPH > 0 ? NEPH : 0)*MB - KEL_MB_NONREN*MB;

// =========================================================================
//  HAEMODYNAMICS / CAPILLARY LEAK
// =========================================================================
dxdt_GLX = -KGLX*(cSVMP + 0.3*cPLA2*INH_PLA2)*F_LEAK*GLX + KGLXR*(1.0 - GLX);
dxdt_PV  = FLRATE - KLEAK*(1.0 - GLX)*PV + KREFILL*(PV0 - PV);

double pvr = (PV > 0.10 ? PV : 0.10) / PV0;
double MAP = MAP0 * pow(pvr, MAP_EXP) * (1.0 - MAP_ANAPHYL*MCA/(MCA + 0.50));
if(MAP < 20.0) MAP = 20.0;

// =========================================================================
//  KIDNEY — an integral, not a level (claim 4)
// =========================================================================
double uflow = (MAP - 45.0)/45.0;
if(uflow < 0) uflow = 0.0;  if(uflow > 1) uflow = 1.0;
uflow = uflow * (NEPH > 0 ? NEPH : 0);

dxdt_CAST = KCAST*MB*(1.0 - 0.7*uflow) - KCLR_CAST*CAST;

double isch = (65.0 - MAP)/65.0;  if(isch < 0) isch = 0.0;
double nloss = (KN_FIB*FGLOSS + KN_CAST*CAST + KN_ISCH*isch + KN_DIR*cSVMP)
               * (NEPH > 0 ? NEPH : 0);
double ceiling = 1.0 - FIBR;  if(ceiling < 0) ceiling = 0.0;
double headroom = ceiling - NEPH;  if(headroom < 0) headroom = 0.0;
dxdt_NEPH = -nloss + KRECOV_NEPH*headroom;
dxdt_FIBR = KFIBROSIS*nloss;
dxdt_SCR  = SCR0*KEL_CR - KEL_CR*(NEPH > 0.02 ? NEPH : 0.02)*SCR;

// =========================================================================
//  INFLAMMATION
// =========================================================================
dxdt_IL6 = G_IL6*(cSVMP + cPLA2 + 12.0*NECFLUX) - KEL_IL6*(IL6 - 3.0);
dxdt_TNF = G_TNF*(cSVMP + cPLA2) - KEL_TNF*(TNF - 4.0);

// =========================================================================
//  ANTIVENOM ADVERSE REACTIONS
//  The acute reaction is driven by the INFUSION RATE of foreign protein, not
//  by the cumulative dose — which is why slowing the infusion works.
// =========================================================================
dxdt_MCA = RHO*AVRATE/60.0 - KEL_MCA*MCA;
dxdt_IC  = K_IC*ICF*((A_C > 0 ? A_C : 0) + (A_P > 0 ? A_P : 0)) - KEL_IC*IC;
dxdt_SS  = K_SS*IC - KEL_SS*SS;

// =========================================================================
//  NEUROMUSCULAR SAFETY FACTOR AND CUMULATIVE HAZARDS
// =========================================================================
// competition at the receptor: elevated acetylcholine effectively dilutes a
// occupancy of a competitive antagonist.  A destroyed terminal releases little
// transmitter, so prolonging its dwell time buys only W_PRE_NEO of the effect.
double neo_amp = 1.0 + EMAX_NEO*cNEO/(EC50_NEO + cNEO);
double brc = (BR < 1.0 ? BR : 1.0);  if(brc < 0) brc = 0.0;
double R_eff = 1.0 - brc/neo_amp;  if(R_eff < 0) R_eff = 0.0;
double tm2 = (TERM > 0 ? TERM : 0) * (1.0 + W_PRE_NEO*(neo_amp - 1.0));
double SF  = SF0 * tm2 * R_eff;
double VCFRAC = (SF - 1.0)/2.0;
if(VCFRAC < 0) VCFRAC = 0.0;  if(VCFRAC > 1) VCFRAC = 1.0;

double unclot = 1.0/(1.0 + pow(fg_s/FG_UNCLOT, 4.0));
double pl_s   = (PLT > 1e-6 ? PLT : 1e-6);
double lowplt = 1.0/(1.0 + pow(pl_s/PLT_LOW, 3.0));

dxdt_HBLD = H_BLEED*unclot*(1.0 + 2.0*lowplt);

double hd = HD_BLEED*unclot*(1.0 + 2.0*lowplt);
if(VCFRAC < 0.15) hd += (ICU > 0.5 ? HD_RESP_VENT : HD_RESP_NOVENT);
double nh = (NEPH_HAZ - (NEPH > 0 ? NEPH : 0))/NEPH_HAZ;
if(nh < 0) nh = 0.0;
hd += (DIALYSIS > 0.5 ? HD_AKI_DIAL : HD_AKI_NODIAL) * nh;
double sh = (60.0 - MAP)/60.0;  if(sh < 0) sh = 0.0;
hd += HD_SHOCK * sh * 10.0;
if(MCA > MCA_ANAPHYL) hd += HD_ANAPHYL;
dxdt_HDTH = hd;

dxdt_VAUC = cVfree;

$TABLE
double cSVMPo = (V_SVMP > 0 ? V_SVMP : 0)/VC_SVMP;
double cSVSPo = (V_SVSP > 0 ? V_SVSP : 0)/VC_SVSP;
double cPLA2o = (V_PLA2 > 0 ? V_PLA2 : 0)/VC_PLA2;
double cTFTXo = (V_TFTX > 0 ? V_TFTX : 0)/VC_TFTX;
capture CVFREE = cSVMPo + cSVSPo + cPLA2o + cTFTXo;
// a venom antigen EIA measures free AND antivenom-bound toxin: antigen can be
// high while the patient is completely protected
capture CVTOTAL = CVFREE + C_SVMP/VC_SVMP + C_SVSP/VC_SVSP
                  + C_PLA2/VC_PLA2 + C_TFTX/VC_TFTX;
capture CA_OUT  = (A_C > 0 ? A_C : 0)/V1A;
capture DETECT  = (CVFREE > CV_LOD) ? 1.0 : 0.0;

double edo   = (EDEMA > 0 ? EDEMA : 0);
capture CPRESS = CPB + CP_MAX*edo/(edo + CP_K50);

double cNEOo = (NEO_C > 0 ? NEO_C : 0)/NEO_V;
double neoo  = 1.0 + EMAX_NEO*cNEOo/(EC50_NEO + cNEOo);
double brco  = (BR < 1.0 ? BR : 1.0);  if(brco < 0) brco = 0.0;
double R_effo = 1.0 - brco/neoo;  if(R_effo < 0) R_effo = 0.0;
capture SFOUT  = SF0*(TERM > 0 ? TERM : 0)*(1.0 + W_PRE_NEO*(neoo - 1.0))*R_effo;
capture VCFRACO = (SFOUT - 1.0)/2.0 < 0 ? 0.0 : ((SFOUT - 1.0)/2.0 > 1 ? 1.0 : (SFOUT - 1.0)/2.0);
capture SBC     = 42.0*VCFRACO;
capture PTOSIS  = (SFOUT < 2.20) ? 1.0 : 0.0;
capture VENT    = (VCFRACO < 0.15) ? 1.0 : 0.0;

double fgo = (FG > 1e-9 ? FG : 1e-9);
capture WBCT20 = (fgo < FG_UNCLOT) ? 1.0 : 0.0;

double pvro = (PV > 0.10 ? PV : 0.10)/PV0;
capture MAPOUT = MAP0*pow(pvro, MAP_EXP)*(1.0 - MAP_ANAPHYL*MCA/(MCA + 0.50));
capture HCT    = (42.0*PV0/(PV > 0.30 ? PV : 0.30) > 72.0) ? 72.0 : 42.0*PV0/(PV > 0.30 ? PV : 0.30);
capture EGFR   = 100.0*(NEPH > 0 ? NEPH : 0);
capture P_BLEED = 1.0 - exp(-HBLD);
capture P_DEATH = 1.0 - exp(-HDTH);
capture ANAPHYL = (MCA > MCA_ANAPHYL) ? 1.0 : 0.0;
'

mod <- mcode("sbe", sbe_code, end = 336, delta = 0.25, atol = 1e-8, rtol = 1e-6)

# =============================================================================
#  SNAKE ARCHETYPES
#  `dose` is the mass of venom protein delivered by a median significant bite,
#  not the maximum milkable gland yield.  Class fractions do not sum to 1
#  because every venom also contains L-amino-acid oxidase, CRISPs,
#  disintegrins, lectins and non-toxic protein this model does not resolve.
# =============================================================================
SNAKES <- list(
  "Daboia russelii (Sri Lanka)" = list(
    ko = "러셀살무사", VDOSE = 63, F_SVMP = 0.28, F_SVSP = 0.14, F_PLA2 = 0.34,
    F_TFTX = 0.01, F_PRE = 0.09, F_MYO = 0.50, F_LOC = 0.30, KOFFM = 1.00, F_LEAK = 1.00),
  "Echis ocellatus" = list(
    ko = "서아프리카카펫바이퍼", VDOSE = 22, F_SVMP = 0.68, F_SVSP = 0.05, F_PLA2 = 0.12,
    F_TFTX = 0.00, F_PRE = 0.00, F_MYO = 0.10, F_LOC = 0.55, KOFFM = 1.00, F_LEAK = 0.60),
  "Naja naja" = list(
    ko = "인도코브라", VDOSE = 40, F_SVMP = 0.10, F_SVSP = 0.02, F_PLA2 = 0.28,
    F_TFTX = 0.55, F_PRE = 0.05, F_MYO = 0.60, F_LOC = 1.30, KOFFM = 1.00, F_LEAK = 0.35),
  "Bungarus caeruleus" = list(
    ko = "인도크레이트", VDOSE = 10, F_SVMP = 0.03, F_SVSP = 0.01, F_PLA2 = 0.66,
    F_TFTX = 0.20, F_PRE = 1.00, F_MYO = 0.05, F_LOC = 0.03, KOFFM = 0.15, F_LEAK = 0.05),
  "Crotalus atrox" = list(
    ko = "서부다이아몬드방울뱀", VDOSE = 55, F_SVMP = 0.42, F_SVSP = 0.14, F_PLA2 = 0.08,
    F_TFTX = 0.00, F_PRE = 0.03, F_MYO = 0.30, F_LOC = 0.70, KOFFM = 1.00, F_LEAK = 0.70),
  "Bothrops asper" = list(
    ko = "중미창머리독사", VDOSE = 50, F_SVMP = 0.50, F_SVSP = 0.06, F_PLA2 = 0.30,
    F_TFTX = 0.00, F_PRE = 0.02, F_MYO = 0.80, F_LOC = 1.50, KOFFM = 1.00, F_LEAK = 0.85),
  "Bitis arietans" = list(
    ko = "퍼프애더", VDOSE = 90, F_SVMP = 0.62, F_SVSP = 0.04, F_PLA2 = 0.16,
    F_TFTX = 0.00, F_PRE = 0.00, F_MYO = 0.50, F_LOC = 1.30, KOFFM = 1.00, F_LEAK = 1.10),
  "Dry bite" = list(
    ko = "건성교상", VDOSE = 0, F_SVMP = 0, F_SVSP = 0, F_PLA2 = 0,
    F_TFTX = 0, F_PRE = 0, F_MYO = 0, F_LOC = 0, KOFFM = 1, F_LEAK = 0)
)

# =============================================================================
#  ANTIVENOM PRODUCTS
#  mgNE_vial for the Indian polyvalent IgG is the LABELLED potency
#  (0.6 mg Daboia venom per mL x 10 mL = 6 mgNE per vial); the Fab and F(ab')2
#  per-vial figures are calibrated so that the labelled initial dose achieves
#  control in the median bite of the species each product is licensed for.
# =============================================================================
AV <- list(
  "Fab"     = list(label = "ovine Fab (CroFab-like)",
                   V1A = 4.0, V2A = 12.0, QA = 2.00, CLA = 0.600,
                   KTIS = 0.050, mgNE_vial = 8, RHO = 0.105, ICF = 0.10),
  "F(ab')2" = list(label = "equine F(ab')2 (ANAVIP-like)",
                   V1A = 3.5, V2A = 4.5, QA = 0.50, CLA = 0.045,
                   KTIS = 0.020, mgNE_vial = 8, RHO = 0.320, ICF = 0.45),
  "IgG"     = list(label = "equine whole IgG (Indian polyvalent ASV-like)",
                   V1A = 3.2, V2A = 2.6, QA = 0.30, CLA = 0.030,
                   KTIS = 0.012, mgNE_vial = 6, RHO = 1.300, ICF = 1.00)
)

EPS <- c(SVMP = 1.20, SVSP = 1.00, PLA2 = 0.60, TFTX = 0.35)

#' mgNE of antivenom needed to neutralise one bite's toxin load exactly.
#' This is the arithmetic behind "how many vials is enough".
stoich_need <- function(snake, venom_mult = 1) {
  s <- SNAKES[[snake]]
  d <- s$VDOSE * venom_mult
  d * s$F_SVMP/EPS["SVMP"] + d * s$F_SVSP/EPS["SVSP"] +
    d * s$F_PLA2/EPS["PLA2"] + d * s$F_TFTX/EPS["TFTX"]
}

#' Apparent TERMINAL half-life of each venom class: absorption-limited by the
#' slow depot (flip-flop), which is 5.5x longer than systemic elimination and is
#' the number an antivenom fragment must outlast.
venom_terminal_thalf <- function() {
  tibble(class = c("SVMP", "SVSP", "PLA2", "3FTx"),
         systemic_thalf_h = log(2) * c(3.5, 3.5, 4.0, 8.0) / c(0.30, 0.40, 0.90, 3.00),
         terminal_thalf_h = log(2) / (c(0.0120, 0.0150, 0.0300, 0.1000) +
                                      c(0.0035, 0.0040, 0.0060, 0.0150)))
}

#' Terminal (beta) half-life of a two-compartment model — computed, not asserted.
tbeta <- function(V1, V2, Q, CL) {
  k10 <- CL/V1; k12 <- Q/V1; k21 <- Q/V2
  s <- k10 + k12 + k21; p <- k10*k21
  log(2) / ((s - sqrt(max(s*s - 4*p, 0)))/2)
}

# =============================================================================
#  SIMULATION DRIVER
# =============================================================================
run_sbe <- function(snake = "Daboia russelii (Sri Lanka)",
                    venom_mult = 1,
                    product = NULL, vials = 0, t_av = NA,
                    av_infuse = 1,
                    repeat_dose = NULL,          # c(time, vials)
                    maintenance = NULL,          # c(start, interval, vials, n)
                    neostigmine = NULL,          # c(t0, interval, n) - 0.5 mg IM
                    varespladib = NULL,          # c(t0, interval, n) - 500 mg PO
                    txa = NULL,                  # c(t0, interval, n) - 1 g IV
                    cryo = NULL,                 # data.frame(time, dFG)
                    fluids = NULL,               # c(t0, t1, L_per_h)
                    icu = TRUE, dialysis = TRUE,
                    tmax = 336, delta = 0.25) {

  s <- SNAKES[[snake]]
  p <- list(VDOSE = s$VDOSE * venom_mult,
            F_SVMP = s$F_SVMP, F_SVSP = s$F_SVSP, F_PLA2 = s$F_PLA2, F_TFTX = s$F_TFTX,
            F_PRE = s$F_PRE, F_MYO = s$F_MYO, F_LOC = s$F_LOC,
            KOFFM = s$KOFFM, F_LEAK = s$F_LEAK,
            ICU = as.numeric(icu), DIALYSIS = as.numeric(dialysis),
            AVDUR = av_infuse,
            AVT0 = 1e6, AVAMT = 0, RPT_T = 1e6, RPT_AMT = 0,
            MSTART = 1e6, MINT = 6, MAMT = 0, MN = 0,
            FLS = 1e6, FLE = 1e6, FLR = 0)

  if (!is.null(product)) {
    a <- AV[[product]]
    p$V1A <- a$V1A; p$V2A <- a$V2A; p$QA <- a$QA; p$CLA <- a$CLA
    p$KTIS <- a$KTIS; p$RHO <- a$RHO; p$ICF <- a$ICF
    if (!is.na(t_av) && vials > 0) {
      p$AVT0 <- t_av; p$AVAMT <- vials * a$mgNE_vial
    }
    if (!is.null(repeat_dose)) {
      p$RPT_T <- repeat_dose[1]; p$RPT_AMT <- repeat_dose[2] * a$mgNE_vial
    }
    if (!is.null(maintenance)) {
      p$MSTART <- maintenance[1]; p$MINT <- maintenance[2]
      p$MAMT <- maintenance[3] * a$mgNE_vial; p$MN <- maintenance[4]
    }
  }
  if (!is.null(fluids)) { p$FLS <- fluids[1]; p$FLE <- fluids[2]; p$FLR <- fluids[3] }

  ev_all <- NULL
  add_ev <- function(e) if (is.null(ev_all)) e else c(ev_all, e)
  if (!is.null(neostigmine))
    ev_all <- add_ev(ev(amt = 0.5, cmt = "NEO_A", time = neostigmine[1],
                        ii = neostigmine[2], addl = neostigmine[3] - 1))
  if (!is.null(varespladib))
    ev_all <- add_ev(ev(amt = 500, cmt = "VAR_A", time = varespladib[1],
                        ii = varespladib[2], addl = varespladib[3] - 1))
  if (!is.null(txa))
    ev_all <- add_ev(ev(amt = 1000, cmt = "TXA_C", time = txa[1],
                        ii = txa[2], addl = txa[3] - 1))
  if (!is.null(cryo))
    for (i in seq_len(nrow(cryo)))
      ev_all <- add_ev(ev(amt = cryo$dFG[i], cmt = "FG", time = cryo$time[i]))

  m <- mod %>% param(p) %>% update(end = tmax, delta = delta)
  if (is.null(ev_all)) {
    # no discrete doses at all: seed a zero-amount record so mrgsim has input
    out <- m %>% ev(amt = 0, cmt = "TXA_C", time = 0) %>% mrgsim()
  } else {
    out <- m %>% ev(ev_all) %>% mrgsim()
  }
  as_tibble(out)
}

# =============================================================================
#  SUMMARY METRICS
#  The recurrence detector is DERIVED from the trajectory (free venom rises
#  again after having been suppressed below 2% of its own peak for >= 4 h).
#  Nothing in the model switches it on.
# =============================================================================
summarise_sbe <- function(d, snake, venom_mult = 1, product = NULL, vials = 0,
                          t_av = NA, name = "") {
  need <- as.numeric(stoich_need(snake, venom_mult))
  mgNE <- if (is.null(product)) 0 else vials * AV[[product]]$mgNE_vial
  tt <- d$time

  unclot <- d$WBCT20 > 0.5
  hrs_unclot <- if (any(unclot)) sum(diff(tt) * head(unclot, -1)) else 0
  t_restore <- if (any(unclot)) {
    last <- max(which(unclot)); if (last < length(tt)) tt[last + 1] else NA_real_
  } else NA_real_

  # fibrinogen rise rate measured at a FIXED substrate level (1.0 g/L) on the
  # way up.  Comparing slopes at the same FG value rather than at each trace's
  # own nadir is what isolates the hepatic synthesis term — see claim 1.
  j0 <- which.min(d$FG); rise1 <- NA_real_
  if (min(d$FG) < 1 && j0 < length(tt)) {
    up <- d$FG[j0:length(tt)]
    k <- which(up >= 1)[1]
    if (!is.na(k) && k > 1) {
      t1 <- tt[j0 + k - 1]; t2 <- min(t1 + 12, max(tt))
      rise1 <- (approx(tt, d$FG, t2)$y - 1) / (t2 - t1) * 24
    }
  }

  # RECURRENT VENOM ANTIGENAEMIA, defined as a rebound: after the acute peak,
  # free venom reaches its FIRST local minimum and a later maximum exceeds that
  # trough >= 2-fold, with the rebound peak above the assay limit of detection.
  # Referencing the rebound to the trough rather than to the acute pre-antivenom
  # peak is the whole point: the acute peak belongs to a different regime.
  # (Defect 7: the first detector compared the rebound with 5% of the ACUTE peak
  # and reported "no recurrence" for a Fab arm whose free venom demonstrably rose
  # 3.6-fold from its trough.  Defect 8: without the LOD gate a trough of
  # 1e-9 mg/L produced a spurious 20000-fold "rebound".)
  cv <- d$CVFREE; tt2 <- d$time
  recur_t <- NA_real_; rebound <- 1; trough <- NA_real_
  LOD <- 0.0020
  if (max(cv) > 0) {
    jpk <- which.max(cv)
    tail_v <- cv[jpk:length(cv)]; tail_t <- tt2[jpk:length(tt2)]
    jmin <- NA_integer_
    if (length(tail_v) > 3) {
      for (m in 2:(length(tail_v) - 1)) {
        if (tail_v[m] <= tail_v[m + 1] && (tail_t[m] - tail_t[1]) >= 4) { jmin <- m; break }
      }
    }
    if (!is.na(jmin) && tail_v[jmin] > 0) {
      trough <- tail_v[jmin]
      after <- tail_v[jmin:length(tail_v)]
      rebound <- max(after) / trough
      if (rebound >= 2 && max(after) > LOD)
        recur_t <- tail_t[jmin + which.max(after) - 1]
    }
  }

  # LATE free-venom exposure (>48 h) is the clean, threshold-free separator
  # between fragment formats: it is the integral the short-lived fragment fails
  # to cover.
  late <- tt2 >= 48
  vauc_late <- if (sum(late) > 1) sum(diff(tt2[late]) * head(cv[late], -1)) else 0
  hrs_detect_late <- if (sum(late) > 1)
    sum(diff(tt2[late]) * head(as.numeric(cv[late] > LOD), -1)) else 0

  # RECURRENT COAGULOPATHY: fibrinogen recovers >= 0.5 g/L from its nadir and
  # then falls again by >= 0.4 g/L.  This is the endpoint the clinical literature
  # reports, and like the venom detector it is read off the trajectory.
  recoag_t <- NA_real_
  jn <- which.min(d$FG); seg <- d$FG[jn:length(d$FG)]
  if (length(seg) > 2) {
    runmax <- cummax(seg)
    hit <- which(runmax - seg[1] >= 0.5 & runmax - seg >= 0.4)
    if (length(hit)) recoag_t <- tt2[jn + hit[1] - 1]
  }

  tibble(
    name = name, snake = snake,
    product = if (is.null(product)) "no antivenom" else AV[[product]]$label,
    vials = vials, t_av = t_av, av_mgNE = mgNE, need_mgNE = need,
    molar_margin = if (need > 0) mgNE/need else Inf,
    FG_nadir = min(d$FG), t_FG_nadir = tt[which.min(d$FG)],
    FG_d7 = approx(tt, d$FG, 168)$y, FG_rise_at_1 = rise1,
    hrs_unclottable = hrs_unclot, t_clottable = t_restore,
    never_unclottable = !any(unclot),
    PLT_nadir = min(d$PLT), FX_nadir = min(d$FX), XDP_peak = max(d$XDP),
    SF_nadir = min(d$SFOUT), SBC_nadir = min(d$SBC),
    TERM_nadir = min(d$TERM), BR_peak = max(d$BR),
    hrs_ptosis = sum(diff(tt) * head(d$PTOSIS, -1)),
    hrs_ventilated = sum(diff(tt) * head(d$VENT, -1)),
    ventilated = any(d$VENT > 0.5),
    NEC_final = tail(d$NEC, 1), CP_peak = max(d$CPRESS),
    CK_peak = max(d$CK), t_CK_peak = tt[which.max(d$CK)], MB_peak = max(d$MB),
    SCR_peak = max(d$SCR), t_SCR_peak = tt[which.max(d$SCR)],
    eGFR_nadir = min(d$EGFR), NEPH_final = tail(d$NEPH, 1),
    FIBR_final = tail(d$FIBR, 1), AKI3 = any(d$SCR >= 2.7),
    MAP_nadir = min(d$MAPOUT), Hct_peak = max(d$HCT), GLX_nadir = min(d$GLX),
    IL6_peak = max(d$IL6), MCA_peak = max(d$MCA), anaphylaxis = any(d$ANAPHYL > 0.5),
    SS_peak = max(d$SS), t_SS_peak = tt[which.max(d$SS)],
    VAUC = tail(d$VAUC, 1),
    p_bleed = tail(d$P_BLEED, 1), p_death = tail(d$P_DEATH, 1),
    cV_peak = max(cv), cV_trough = trough, rebound_factor = rebound,
    VAUC_late = vauc_late, hrs_detectable_late = hrs_detect_late,
    recur_t = recur_t, recurrence = !is.na(recur_t),
    recoag_t = recoag_t, recoagulopathy = !is.na(recoag_t)
  )
}

# =============================================================================
#  27 TREATMENT SCENARIOS
#  Grouped by which of the four claims each one is there to test.
# =============================================================================
RV    <- "Daboia russelii (Sri Lanka)"
CA    <- "Crotalus atrox"
NN    <- "Naja naja"
BC    <- "Bungarus caeruleus"
FLUID <- c(1, 8, 0.30)

# The 27 scenarios, grouped by which claim each one exists to test.  These are
# the same 27 that `sbe_reference_model.py` integrates; every number quoted in
# README.md comes from that run.
SCENARIOS <- list(
  # ================= claim 1: depth versus slope ===========================
  list(name = "01 Russell's viper — no antivenom",
       args = list(snake = RV, fluids = FLUID),
       note = "untreated reference; VICC + AKI + capillary leak"),
  list(name = "02 Russell's viper — ASV 10 v @ 1 h",
       args = list(snake = RV, product = "IgG", vials = 10, t_av = 1, fluids = FLUID),
       note = "claim 1: earliest realistic antivenom"),
  list(name = "03 Russell's viper — ASV 10 v @ 4 h",
       args = list(snake = RV, product = "IgG", vials = 10, t_av = 4, fluids = FLUID),
       note = "claim 1: median real-world delay"),
  list(name = "04 Russell's viper — ASV 10 v @ 12 h",
       args = list(snake = RV, product = "IgG", vials = 10, t_av = 12, fluids = FLUID),
       note = "claim 1: late antivenom, substrate already at the floor"),
  list(name = "05 Russell's viper — ASV 20 v @ 4 h",
       args = list(snake = RV, product = "IgG", vials = 20, t_av = 4, fluids = FLUID),
       note = "claim 1: does doubling the dose buy back the lost time?"),
  list(name = "06 Russell's viper — ASV 10 v @ 4 h + cryoprecipitate",
       args = list(snake = RV, product = "IgG", vials = 10, t_av = 4, fluids = FLUID,
                   cryo = data.frame(time = c(4.5, 10), dFG = c(0.9, 0.9))),
       note = "claim 1: substitute for the liver instead of waiting for it"),
  list(name = "07 Russell's viper — ASV 10 v @ 4 h + tranexamic acid",
       args = list(snake = RV, product = "IgG", vials = 10, t_av = 4, fluids = FLUID,
                   txa = c(4, 8, 3)),
       note = "claim 1: only 5% of the fibrinogen loss is plasmin-mediated"),

  # ================= claim 2: recurrence ====================================
  # Crotalus atrox + ovine Fab is the setting in which recurrent venom
  # antigenaemia is actually documented, and the CroFab label's mandatory
  # 2-vials-q6h-x3 maintenance schedule is the fix the inequality predicts.
  # Arms 10-12 hold the MOLAR MARGIN constant while changing the bite size,
  # which is what separates "the margin matters" from "the vial count matters".
  list(name = "08 Crotalus atrox — Fab 6 v @ 2 h (margin 1.4x)",
       args = list(snake = CA, product = "Fab", vials = 6, t_av = 2, fluids = FLUID),
       note = "claim 2: short-lived fragment against a multi-day venom tail"),
  list(name = "09 Crotalus atrox — F(ab')2 6 v @ 2 h (same mgNE, margin 1.4x)",
       args = list(snake = CA, product = "F(ab')2", vials = 6, t_av = 2, fluids = FLUID),
       note = "claim 2: identical molar dose, 5.8x the half-life"),
  list(name = "10 Crotalus atrox — Fab 12 v @ 2 h (margin 2.8x)",
       args = list(snake = CA, product = "Fab", vials = 12, t_av = 2, fluids = FLUID),
       note = "claim 2: doubling works — at THIS bite size"),
  list(name = "11 Crotalus atrox 2x venom — Fab 12 v @ 2 h (margin 1.4x again)",
       args = list(snake = CA, venom_mult = 2, product = "Fab", vials = 12,
                   t_av = 2, fluids = FLUID),
       note = "claim 2: same margin, twice the tail — recurrence returns"),
  list(name = "12 Crotalus atrox 2x venom — F(ab')2 12 v @ 2 h (margin 1.4x)",
       args = list(snake = CA, venom_mult = 2, product = "F(ab')2", vials = 12,
                   t_av = 2, fluids = FLUID),
       note = "claim 2: the long fragment does not care about the tail"),
  list(name = "13 Crotalus atrox — Fab 6 v @ 2 h + 2 v q6h x 3 (label schedule)",
       args = list(snake = CA, product = "Fab", vials = 6, t_av = 2, fluids = FLUID,
                   maintenance = c(8, 6, 2, 3)),
       note = "claim 2: the fix is a schedule, not a bigger bolus"),

  # ================= claim 3: two paralyses =================================
  list(name = "14 Naja naja — no antivenom",
       args = list(snake = NN),
       note = "claim 3: postsynaptic paralysis, untreated"),
  list(name = "15 Naja naja — neostigmine q4h only, no antivenom",
       args = list(snake = NN, neostigmine = c(2, 4, 10)),
       note = "claim 3: the competition effect, isolated"),
  list(name = "16 Naja naja — ASV 10 v @ 4 h (label dose, margin 0.70x)",
       args = list(snake = NN, product = "IgG", vials = 10, t_av = 4),
       note = "claim 3: the label dose is under-stoichiometric for a cobra"),
  list(name = "17 Naja naja — ASV 10 v @ 4 h + neostigmine q4h",
       args = list(snake = NN, product = "IgG", vials = 10, t_av = 4,
                   neostigmine = c(4, 4, 10)),
       note = "claim 3: acetylcholine out-competes a reversible block"),
  list(name = "18 Naja naja — ASV 20 v @ 4 h (margin 1.40x)",
       args = list(snake = NN, product = "IgG", vials = 20, t_av = 4),
       note = "claim 3: what adequate neutralisation looks like"),
  list(name = "19 Bungarus caeruleus — ASV 10 v @ 6 h",
       args = list(snake = BC, product = "IgG", vials = 10, t_av = 6),
       note = "claim 3: presynaptic destruction already done"),
  list(name = "20 Bungarus caeruleus — ASV 10 v @ 6 h + neostigmine q4h",
       args = list(snake = BC, product = "IgG", vials = 10, t_av = 6,
                   neostigmine = c(6, 4, 10)),
       note = "claim 3: nothing to out-compete"),
  list(name = "21 Bungarus caeruleus — ASV 10 v @ 0.5 h (prevention)",
       args = list(snake = BC, product = "IgG", vials = 10, t_av = 0.5),
       note = "claim 3: the same drug, before the destruction"),
  list(name = "22 Bungarus caeruleus — ASV 10 v @ 6 h, NO ventilator",
       args = list(snake = BC, product = "IgG", vials = 10, t_av = 6, icu = FALSE),
       note = "claim 3: what actually saves a krait bite"),

  # ================= where antivenom cannot go ==============================
  list(name = "23 Bothrops asper — F(ab')2 10 v @ 4 h",
       args = list(snake = "Bothrops asper", product = "F(ab')2", vials = 10,
                   t_av = 4, fluids = FLUID),
       note = "systemic control, local failure"),
  list(name = "24 Bothrops asper — F(ab')2 10 v @ 4 h + varespladib @ 0.5 h",
       args = list(snake = "Bothrops asper", product = "F(ab')2", vials = 10,
                   t_av = 4, fluids = FLUID, varespladib = c(0.5, 12, 8)),
       note = "a small molecule reaches tissue an IgG cannot"),

  # ================= other settings =========================================
  list(name = "25 Echis ocellatus — polyvalent IgG 3 v @ 3 h",
       args = list(snake = "Echis ocellatus", product = "IgG", vials = 3,
                   t_av = 3, fluids = FLUID),
       note = "West African pure-SVMP coagulopathy"),
  list(name = "26 Dry bite — ASV 10 v @ 1 h",
       args = list(snake = "Dry bite", product = "IgG", vials = 10, t_av = 1),
       note = "all of the risk, none of the benefit"),
  list(name = "27 Russell's viper 2.5x venom — ASV 10 v @ 4 h",
       args = list(snake = RV, venom_mult = 2.5, product = "IgG", vials = 10,
                   t_av = 4, fluids = FLUID),
       note = "the label dose against a large bite: margin 0.39x")
)

run_all <- function() {
  res <- lapply(SCENARIOS, function(sc) {
    d <- do.call(run_sbe, sc$args)
    s <- summarise_sbe(d, sc$args$snake,
                       venom_mult = if (is.null(sc$args$venom_mult)) 1 else sc$args$venom_mult,
                       product = sc$args$product,
                       vials = if (is.null(sc$args$vials)) 0 else sc$args$vials,
                       t_av = if (is.null(sc$args$t_av)) NA else sc$args$t_av,
                       name = sc$name)
    s$note <- sc$note
    list(sim = d, sum = s)
  })
  names(res) <- vapply(SCENARIOS, `[[`, "", "name")
  res
}

# =============================================================================
#  ANALYSIS FUNCTIONS — one per claim
# =============================================================================

#' Claim 1.  Antivenom timing moves the DEPTH of the fibrinogen nadir.  The
#' rise rate measured at a fixed fibrinogen level does not move, because it is
#' the hepatic synthesis rate and no antivenom term appears in it.
analyse_depth_vs_slope <- function(res) {
  k <- c("01 Russell's viper — no antivenom",
         "02 Russell's viper — ASV 10 v @ 1 h",
         "03 Russell's viper — ASV 10 v @ 4 h",
         "04 Russell's viper — ASV 10 v @ 12 h",
         "05 Russell's viper — ASV 20 v @ 4 h",
         "06 Russell's viper — ASV 10 v @ 4 h + cryoprecipitate",
         "07 Russell's viper — ASV 10 v @ 4 h + tranexamic acid")
  bind_rows(lapply(res[k], `[[`, "sum")) %>%
    select(name, FG_nadir, hrs_unclottable, FG_rise_at_1, FG_d7, SCR_peak,
           p_bleed, AKI3)
}

#' Claim 2.  Recurrence is a property of the fragment.  Compare the terminal
#' half-life of the antivenom with that of the venom input.
analyse_recurrence <- function(res) {
  k <- c("08 Crotalus atrox — Fab 6 v @ 2 h (margin 1.4x)",
         "09 Crotalus atrox — F(ab')2 6 v @ 2 h (same mgNE, margin 1.4x)",
         "10 Crotalus atrox — Fab 12 v @ 2 h (margin 2.8x)",
         "11 Crotalus atrox 2x venom — Fab 12 v @ 2 h (margin 1.4x again)",
         "12 Crotalus atrox 2x venom — F(ab')2 12 v @ 2 h (margin 1.4x)",
         "13 Crotalus atrox — Fab 6 v @ 2 h + 2 v q6h x 3 (label schedule)")
  bind_rows(lapply(res[k], `[[`, "sum")) %>%
    select(name, av_mgNE, molar_margin, rebound_factor, recurrence, recur_t,
           VAUC_late, hrs_detectable_late, FG_nadir)
}

#' Claim 3.  Two paralyses, one bedside appearance, opposite pharmacology.
analyse_two_paralyses <- function(res) {
  k <- c("14 Naja naja — no antivenom",
         "15 Naja naja — neostigmine q4h only, no antivenom",
         "16 Naja naja — ASV 10 v @ 4 h (label dose, margin 0.70x)",
         "17 Naja naja — ASV 10 v @ 4 h + neostigmine q4h",
         "18 Naja naja — ASV 20 v @ 4 h (margin 1.40x)",
         "19 Bungarus caeruleus — ASV 10 v @ 6 h",
         "20 Bungarus caeruleus — ASV 10 v @ 6 h + neostigmine q4h",
         "21 Bungarus caeruleus — ASV 10 v @ 0.5 h (prevention)",
         "22 Bungarus caeruleus — ASV 10 v @ 6 h, NO ventilator")
  bind_rows(lapply(res[k], `[[`, "sum")) %>%
    select(name, BR_peak, TERM_nadir, SF_nadir, hrs_ptosis, hrs_ventilated, p_death)
}

#' Claim 4.  The kidney is an integral: creatinine peaks days after the venom
#' has gone, and the size of the peak tracks the free-venom AUC.
analyse_kidney_integral <- function(res) {
  k <- c("01 Russell's viper — no antivenom",
         "02 Russell's viper — ASV 10 v @ 1 h",
         "03 Russell's viper — ASV 10 v @ 4 h",
         "04 Russell's viper — ASV 10 v @ 12 h",
         "27 Russell's viper 2.5x venom — ASV 10 v @ 4 h")
  bind_rows(lapply(res[k], `[[`, "sum")) %>%
    select(name, VAUC, SCR_peak, t_SCR_peak, eGFR_nadir, FIBR_final)
}

#' Where antivenom cannot go, and what can get there.
analyse_local <- function(res) {
  k <- c("23 Bothrops asper — F(ab')2 10 v @ 4 h",
         "24 Bothrops asper — F(ab')2 10 v @ 4 h + varespladib @ 0.5 h")
  bind_rows(lapply(res[k], `[[`, "sum")) %>%
    select(name, NEC_final, CK_peak, CP_peak, TERM_nadir, VAUC)
}

#' How many vials is "enough" — the stoichiometric table.
analyse_stoichiometry <- function() {
  sn <- names(SNAKES)[vapply(SNAKES, function(s) s$VDOSE > 0, TRUE)]
  tibble(
    snake = sn,
    venom_mg = vapply(sn, function(s) SNAKES[[s]]$VDOSE, 0),
    mgNE_need = vapply(sn, function(s) as.numeric(stoich_need(s)), 0)
  ) %>%
    mutate(IgG_vials = mgNE_need/AV$IgG$mgNE_vial,
           Fab_vials = mgNE_need/AV$Fab$mgNE_vial,
           margin_at_10_IgG_vials = 10 * AV$IgG$mgNE_vial / mgNE_need)
}

#' The recurrence inequality, evaluated.
analyse_halflives <- function() {
  tibble(
    entity = c("ovine Fab", "equine F(ab')2", "equine whole IgG",
               "venom SVMP (input)", "venom SVSP (input)",
               "venom PLA2 (input)", "venom 3FTx (input)"),
    t_half_beta_h = c(
      tbeta(4.0, 12.0, 2.00, 0.600),
      tbeta(3.5,  4.5, 0.50, 0.045),
      tbeta(3.2,  2.6, 0.30, 0.030),
      log(2)/(0.0120 + 0.0035),
      log(2)/(0.0150 + 0.0040),
      log(2)/(0.0300 + 0.0060),
      log(2)/(0.1000 + 0.0150))
  )
}

# =============================================================================
#  VIRTUAL POPULATION
#  A trial arm is a MIXTURE.  A single median patient cannot reproduce a trial
#  event rate, because every endpoint here is a threshold and the mean of a
#  threshold is not the threshold of a mean.  Between-subject variability is
#  dominated by injected venom mass (a bite is not a calibrated injection) and
#  by presentation delay.
# =============================================================================
virtual_population <- function(n = 300, seed = 20260805) {
  set.seed(seed)
  arms <- list(
    list(label = "Daboia + ASV 10 v, delay ~ 4 h", snake = RV, product = "IgG", vials = 10, delay = 4),
    list(label = "Daboia + ASV 10 v, delay ~ 1 h", snake = RV, product = "IgG", vials = 10, delay = 1),
    list(label = "Daboia + ASV 20 v, delay ~ 4 h", snake = RV, product = "IgG", vials = 20, delay = 4),
    list(label = "Daboia, no antivenom",           snake = RV, product = NULL,  vials = 0,  delay = NA),
    list(label = "C. atrox + Fab 6 v @ 2 h",       snake = "Crotalus atrox", product = "Fab",     vials = 6, delay = 2),
    list(label = "C. atrox + F(ab')2 6 v @ 2 h",   snake = "Crotalus atrox", product = "F(ab')2", vials = 6, delay = 2)
  )
  bind_rows(lapply(arms, function(a) {
    vm <- pmin(pmax(rlnorm(n, log(1) - 0.5*0.55^2, 0.55), 0.10), 3.5)
    td <- if (is.na(a$delay)) rep(NA_real_, n) else
      pmax(0.3, rlnorm(n, log(a$delay) - 0.5*0.5^2, 0.5))
    bind_rows(lapply(seq_len(n), function(i) {
      d <- run_sbe(snake = a$snake, venom_mult = vm[i], product = a$product,
                   vials = a$vials, t_av = td[i], fluids = FLUID,
                   tmax = 240, delta = 0.5)
      s <- summarise_sbe(d, a$snake, vm[i], a$product, a$vials, td[i], a$label)
      s$arm <- a$label; s
    }))
  })) %>%
    group_by(arm) %>%
    summarise(
      n = n(),
      pct_unclottable   = 100*mean(hrs_unclottable > 0),
      median_hrs_unclot = median(hrs_unclottable),
      pct_recurrence    = 100*mean(recurrence),
      pct_detect_late   = 100*mean(hrs_detectable_late > 0),
      pct_AKI_any       = 100*mean(SCR_peak >= 1.35),
      pct_AKI3          = 100*mean(AKI3),
      pct_ventilated    = 100*mean(ventilated),
      mean_p_bleed      = 100*mean(p_bleed),
      mean_p_death      = 100*mean(p_death),
      # the labelled 10-vial dose is a 0.98x margin against the MEDIAN Daboia
      # bite, so roughly half of a real bite distribution is under-dosed BY
      # CONSTRUCTION.  Nothing about the antivenom is failing.
      pct_underdosed    = 100*mean(molar_margin < 1),
      .groups = "drop")
}

# =============================================================================
#  PLOTS
# =============================================================================
plot_claim1 <- function(res) {
  k <- c("01 Russell's viper — no antivenom",
         "02 Russell's viper — ASV 10 v @ 1 h",
         "03 Russell's viper — ASV 10 v @ 4 h",
         "04 Russell's viper — ASV 10 v @ 12 h")
  bind_rows(lapply(k, function(n) res[[n]]$sim %>% mutate(scenario = n))) %>%
    filter(time <= 120) %>%
    ggplot(aes(time, FG, colour = scenario)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    geom_line(linewidth = 0.8) +
    labs(x = "hours post-bite", y = "fibrinogen (g/L)",
         title = "Antivenom timing changes the DEPTH of the nadir",
         subtitle = "the rise rate out of the nadir is the liver's, not the antivenom's; dashed line = 20WBCT threshold") +
    theme_bw() + theme(legend.position = "bottom")
}

plot_claim2 <- function(res) {
  k <- c("08 Crotalus atrox — Fab 6 v @ 2 h (margin 1.4x)",
         "09 Crotalus atrox — F(ab')2 6 v @ 2 h (same mgNE, margin 1.4x)",
         "11 Crotalus atrox 2x venom — Fab 12 v @ 2 h (margin 1.4x again)",
         "13 Crotalus atrox — Fab 6 v @ 2 h + 2 v q6h x 3 (label schedule)")
  bind_rows(lapply(k, function(n) res[[n]]$sim %>% mutate(scenario = n))) %>%
    ggplot(aes(time, CVFREE, colour = scenario)) +
    geom_line(linewidth = 0.8) +
    scale_y_log10() +
    labs(x = "hours post-bite", y = "free venom (mg/L, log scale)",
         title = "Recurrence is a ratio of two half-lives",
         subtitle = "doubling the front-loaded Fab dose delays recurrence by one Fab half-life; it does not abolish it") +
    theme_bw() + theme(legend.position = "bottom")
}

plot_claim3 <- function(res) {
  k <- c("16 Naja naja — ASV 10 v @ 4 h (label dose, margin 0.70x)",
         "17 Naja naja — ASV 10 v @ 4 h + neostigmine q4h",
         "19 Bungarus caeruleus — ASV 10 v @ 6 h",
         "20 Bungarus caeruleus — ASV 10 v @ 6 h + neostigmine q4h")
  bind_rows(lapply(k, function(n) res[[n]]$sim %>% mutate(scenario = n))) %>%
    select(time, scenario, SFOUT, BR, TERM) %>%
    pivot_longer(c(SFOUT, BR, TERM)) %>%
    ggplot(aes(time, value, colour = scenario)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~name, scales = "free_y") +
    labs(x = "hours post-bite", y = NULL,
         title = "Neostigmine moves an occupancy (BR) and cannot move a destruction (TERM)") +
    theme_bw() + theme(legend.position = "bottom")
}

# =============================================================================
#  DEMONSTRATION
# =============================================================================
if (interactive()) {
  print(analyse_halflives())
  print(venom_terminal_thalf())
  print(analyse_stoichiometry())
  res <- run_all()
  print(bind_rows(lapply(res, `[[`, "sum")), n = 27, width = Inf)
  print(analyse_depth_vs_slope(res))
  print(analyse_recurrence(res))
  print(analyse_two_paralyses(res))
  print(analyse_kidney_integral(res))
  print(analyse_local(res))
  print(virtual_population())
  print(plot_claim1(res)); print(plot_claim2(res)); print(plot_claim3(res))
}
