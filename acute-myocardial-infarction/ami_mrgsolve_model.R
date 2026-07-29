# =============================================================================
#  ami_mrgsolve_model.R
#  Acute Myocardial Infarction (STEMI) — Quantitative Systems Pharmacology model
# =============================================================================
#
#  THE CLAIM THIS MODEL IS BUILT TO MAKE
#  -------------------------------------
#  An infarct is not "the artery closed, therefore the muscle died".
#  It is TWO RACES and ONE BIFURCATION:
#
#    Race 1 (minutes-hours).  A wavefront of necrosis crawls from endocardium to
#      epicardium.  Its speed is set by how far residual COLLATERAL flow falls
#      short of BASAL metabolic demand -- not of resting demand.  Tissue that
#      can still pay its basal bill stops contracting and survives indefinitely
#      (hibernation); tissue that cannot, dies.  The clock is time-to-reperfusion.
#
#    Race 2 (the first minutes of reflow).  Ischemia leaves two things behind:
#      ACID, which holds the mitochondrial permeability transition pore SHUT,
#      and SUCCINATE, which becomes an oxidant burst the instant oxygen returns.
#      Restoring flow washes out the acid and ignites the succinate at the same
#      moment.  Reperfusion is simultaneously the treatment and a second injury.
#
#    The bifurcation (weeks-months).  Laplace wall stress -> dilation -> more
#      wall stress is a POSITIVE loop.  Concentric hypertrophy is the NEGATIVE
#      loop opposing it.  Whether a given infarct settles into a stable
#      ventricle or runs away into heart failure is a bifurcation, not a dose.
#
#  There is no parameter called golden_hour, critical_infarct_size,
#  no_reflow_switch, reperfusion_injury_fraction, therapeutic_window or
#  severity.  All of those are computed results.  See README.md for the numbers
#  and ami_reference_check.py for a scipy re-implementation that regenerates
#  every one of them (it is the numerical source of truth; this file must
#  agree with it).
#
#  STRUCTURE — 82 ODEs
#    5 transmural layers x 8 states = 40   (E G H C P NI NR SUC)
#    29 global states                      (injury, inflammation, mechanics,
#                                           neurohormonal, biomarkers)
#    13 pharmacokinetic compartments
#
#  UNITS.  Time is HOURS from symptom onset throughout.  Flows are fractions of
#  normal resting myocardial blood flow (1 = normal).  Volumes are mL, mass g,
#  pressures mmHg, drug amounts mg.
#
#  USAGE
#    library(mrgsolve); library(dplyr)
#    mod <- mread("ami_mrgsolve_model.R")
#    out <- mod %>% param(T_PCI = 1.5, T_ASA = 1.0, T_TIC = 1.0, T_HEP = 1.0,
#                         T_RAM = 24, T_MET_PO = 24, T_EPL = 48, T_EMP = 24) %>%
#                   mrgsim(end = 24*180, delta = 0.05)
#    plot(out, IS + MVO + EF + EDV ~ time)
#
#  NOTE ON DOSING.  Every intervention is driven by a TIME PARAMETER
#  (T_PCI, T_LYT, T_RAM, ...) evaluated inside $ODE as a smooth switch, NOT by
#  mrgsolve event records.  This is deliberate: it makes the R model and the
#  Python reference bit-for-bit comparable, and it lets a scenario be specified
#  entirely through param().  Set a time parameter to a large value (1e6) to
#  withhold that intervention.  If you prefer real ev() dosing for the PK
#  compartments, the cmt names are exposed and additive.
#
#  Author: QSP disease-model library (see repository README).
# =============================================================================

$PROB
# Acute Myocardial Infarction (STEMI) — QSP
# 82 ODEs | 5 transmural layers | two races and one bifurcation

$PARAM @annotated
// ---- anatomy of the ischemic bed -------------------------------------------
AAR      : 0.35   : Area at risk, fraction of LV (proximal LAD)
COLL     : 0.10   : Mean collateral flow, fraction of normal resting flow
W1       : 0.20   : Transmural mass weight, layer 1 (subendocardial)
W2       : 0.20   : Transmural mass weight, layer 2
W3       : 0.20   : Transmural mass weight, layer 3 (midwall)
W4       : 0.20   : Transmural mass weight, layer 4
W5       : 0.20   : Transmural mass weight, layer 5 (subepicardial)
GC1      : 0.05   : Collateral flow multiplier, layer 1 - THE key asymmetry
GC2      : 0.35   : Collateral flow multiplier, layer 2
GC3      : 0.75   : Collateral flow multiplier, layer 3
GC4      : 1.35   : Collateral flow multiplier, layer 4
GC5      : 2.50   : Collateral flow multiplier, layer 5
GM1      : 1.60   : Susceptibility to microvascular plugging, layer 1
GM2      : 1.30   : Susceptibility to microvascular plugging, layer 2
GM3      : 1.00   : Susceptibility to microvascular plugging, layer 3
GM4      : 0.70   : Susceptibility to microvascular plugging, layer 4
GM5      : 0.40   : Susceptibility to microvascular plugging, layer 5

// ---- energetics: the survival criterion is BASAL demand --------------------
KE       : 30.0   : 1/h fall of energy charge per unit energy deficit
KREC     : 0.35   : Recovery of energy charge per unit surplus
EFLOOR   : 0.02   : Soft floor keeping energy charge non-negative
DEM      : 1.00   : Demand scale at baseline rate-pressure product
DEM_BAS  : 0.20   : Basal (non-contractile) fraction of demand
F_ANAERO : 0.18   : Maximum anaerobic contribution to supply
HI50     : 0.50   : Product inhibition of glycolysis by its own H+
KDIAST   : 1.00   : Collateral perfusion is diastolic (slower heart = more flow)
KGLY     : 5.00   : 1/h glycogen consumption
KGLYREP  : 0.25   : 1/h glycogen repletion (requires flow)

// ---- acidosis: made anaerobically, removed only by flow -------------------
KH       : 8.00   : 1/h H+ generation per unit anaerobic flux
KHW      : 4.00   : 1/h H+ washout - REQUIRES FLOW (this is the pH paradox)

// ---- ischemic succinate and the reperfusion oxidant burst -----------------
KSUC     : 2.50   : 1/h succinate accumulation per unit anaerobic flux
KSUCO    : 6.00   : 1/h succinate oxidation once oxygen returns
KROS     : 4.00   : 1/h ROS per unit succinate burned (reverse electron transport)
KROS_N   : 0.60   : 1/h ROS from neutrophil NADPH oxidase
KROSD    : 3.00   : 1/h oxidant scavenging

// ---- calcium ---------------------------------------------------------------
KCA      : 1.20   : 1/h Ca2+ gain from failing SERCA and NCX
NCA      : 1.50   : Exponent on energy deficit for Ca2+ gain
KCA_RP   : 2.50   : 1/h reflow Ca2+ surge, proportional to flow x accumulated acid
KCAO     : 3.00   : 1/h Ca2+ extrusion, requires energy

// ---- mitochondrial permeability transition pore ---------------------------
KP_ON    : 8.00   : 1/h pore opening rate
CA50     : 0.55   : Ca2+ load for half-maximal pore opening
NP_CA    : 3.00   : Hill coefficient on Ca2+
R50      : 0.50   : Oxidant level for half-maximal pore opening
NP_R     : 2.00   : Hill coefficient on oxidant
HP50     : 0.30   : Acid gate constant - low pH KEEPS THE PORE SHUT
KP_OFF   : 2.00   : 1/h pore resealing, requires energy

// ---- necrosis --------------------------------------------------------------
KNI      : 0.50   : 1/h ischemic necrosis rate constant
NNI      : 4.00   : Power law on energy deficit (deliberately no threshold)
KOSM     : 0.60   : Osmotic / contracture acceleration by acid
KNR      : 2.50   : 1/h mPTP-mediated necrosis rate constant

// ---- microvascular obstruction --------------------------------------------
KMVO_E   : 1.20   : Endothelial swelling from energy failure
KMVO_N   : 0.020  : Neutrophil capillary plugging
KMVO_P   : 0.40   : Distal atherothrombotic embolisation
KMVO_R   : 0.012  : 1/h resolution of microvascular obstruction (weeks)

// ---- epicardial thrombus ---------------------------------------------------
KREG     : 0.30   : 1/h platelet-driven thrombus (re)growth
KLYS0    : 0.02   : 1/h endogenous lysis
KLYS_PL  : 3.50   : 1/h plasmin-driven lysis
KPL50    : 0.35   : Plasmin for half-maximal lysis
NPL      : 1.50   : Hill coefficient on plasmin
KPCI     : 60.0   : 1/h mechanical thrombus removal at PCI
PCI_FAIL : 0.05   : Residual thrombus burden after a successful PCI

// ---- fibrinolytic pharmacology (tenecteplase) -----------------------------
TNK_V1   : 4.2    : Tenecteplase central volume (L)
TNK_CL   : 6.3    : Tenecteplase clearance (L/h)
TNK_Q    : 1.5    : Tenecteplase intercompartmental clearance (L/h)
TNK_V2   : 6.6    : Tenecteplase peripheral volume (L)
KPG      : 8.00   : 1/h plasminogen to plasmin per unit tenecteplase
KPLD     : 10.0   : 1/h plasmin inhibition by alpha2-antiplasmin
KI_PAI   : 1.20   : PAI-1 inhibition constant
KPGD     : 1.50   : 1/h systemic plasminogen consumption
KPGS     : 0.06   : 1/h plasminogen resynthesis
KFIBD    : 0.50   : 1/h fibrinogenolysis per unit plasmin
KFIBS    : 0.02   : 1/h fibrinogen resynthesis
KPAID    : 0.35   : 1/h PAI-1 turnover
KPAI_I   : 0.60   : IL-6 driven PAI-1 induction

// ---- antiplatelet and anticoagulant ---------------------------------------
TIC_KA   : 1.20   : Ticagrelor absorption rate (1/h)
TIC_F    : 0.36   : Ticagrelor bioavailability
TIC_V    : 88.0   : Ticagrelor volume (L)
TIC_CL   : 13.5   : Ticagrelor clearance (L/h)
TIC_EC50 : 0.15   : Ticagrelor EC50 for P2Y12 inhibition (mg/L)
TIC_EMAX : 0.95   : Maximum P2Y12 inhibition
ASA_EFF  : 0.45   : Aspirin contribution to platelet inhibition
HEP_EFF  : 0.50   : Heparin contribution to blocking thrombus regrowth
KPLT_PLN : 1.80   : Plasmin ACTIVATES platelets - why a lytic undoes itself

// ---- beta blockade ---------------------------------------------------------
MET_KA   : 1.10   : Metoprolol absorption rate (1/h)
MET_F    : 0.40   : Metoprolol bioavailability
MET_V    : 250.0  : Metoprolol volume (L)
MET_CL   : 60.0   : Metoprolol clearance (L/h)
MET_EC50 : 0.030  : Metoprolol EC50 (mg/L)
MET_EMAX : 1.00   : Maximum beta-blockade
BB_HR    : 0.28   : Maximum fractional heart-rate reduction
BB_NE    : 0.55   : Maximum blockade of the norepinephrine effect

// ---- RAAS, SGLT2, anti-inflammatory, mPTP inhibitor -----------------------
RAM_KA   : 1.50   : Ramipril absorption rate (1/h)
RAM_F    : 0.28   : Ramipril bioavailability
RAM_V    : 90.0   : Ramipril volume (L)
RAM_CL   : 8.00   : Ramipril clearance (L/h)
RAM_KD   : 0.010  : Ramipril concentration for half-maximal ACE inhibition (mg/L)
ACEI_EMAX: 0.72   : Maximum ACE inhibition
EPL_KEL  : 0.14   : Eplerenone elimination rate (1/h)
EPL_KD   : 30.0   : Eplerenone amount for half-maximal MR blockade (mg)
MRA_EMAX : 0.68   : Maximum mineralocorticoid receptor blockade
EMP_KEL  : 0.055  : Empagliflozin elimination rate (1/h)
EMP_KD   : 6.00   : Empagliflozin amount for half-maximal effect (mg)
SGLT_EMAX: 1.00   : Maximum SGLT2 effect
COLC_KEL : 0.035  : Colchicine elimination rate (1/h)
COLC_KD  : 0.20   : Colchicine amount for half-maximal NLRP3 inhibition (mg)
COLC_EMAX: 0.45   : Maximum colchicine inhibition of IL-1beta production
CAN_KEL  : 0.0012 : Canakinumab elimination rate (1/h)
CAN_KD   : 40.0   : Canakinumab amount for half-maximal IL-1beta neutralisation (mg)
CAN_EMAX : 0.85   : Maximum IL-1beta neutralisation
CSA_KEL  : 0.058  : Cyclosporine elimination rate (1/h)
CSA_KD   : 120.0  : Cyclosporine amount for half-maximal mPTP inhibition (mg)
CSA_EMAX : 0.70   : Maximum mPTP inhibition
ARNI     : 0.0    : Flag - sacubitril/valsartan in place of a plain ACE inhibitor

// ---- inflammation and healing ---------------------------------------------
KD1      : 1.00   : DAMP release per unit necrosis rate
KD2      : 0.35   : 1/h DAMP clearance
KN1      : 0.25   : Neutrophil influx per unit DAMP x perfusion
KN2      : 0.075  : 1/h neutrophil clearance
KM1      : 0.050  : M1 macrophage recruitment
KM1D     : 0.020  : 1/h M1 clearance
KSW      : 0.055  : M1 to M2 conversion rate
KM2D     : 0.014  : 1/h M2 clearance
KSW50    : 1.20   : Constant of the neutrophil-clearance switch
KI1      : 0.10   : IL-1beta production
KI1D     : 0.09   : 1/h IL-1beta clearance
KI6      : 0.16   : IL-6 production
KI6D     : 0.12   : 1/h IL-6 clearance
KC1      : 0.09   : CRP production per unit IL-6
KC2      : 0.035  : 1/h CRP clearance
KT1      : 0.045  : TGF-beta production per unit M2
KT2      : 0.020  : 1/h TGF-beta clearance
KMY      : 0.030  : Myofibroblast activation
KMYD     : 0.012  : 1/h myofibroblast loss
KCO      : 0.090  : Collagen deposition per unit myofibroblast
KCOD     : 0.012  : Collagen degradation per unit MMP
COLMAX   : 8.00   : Saturating collagen content
KMM      : 0.050  : MMP activation
KMMD     : 0.055  : 1/h MMP inactivation
KSC      : 1.10   : Collagen giving half-maximal scar tensile strength
SCARSET  : 0.85   : Scar strength at which the infarct stops stretching

// ---- LV mechanics ----------------------------------------------------------
EDV0     : 110.0  : Baseline end-diastolic volume (mL)
ESV0     : 42.0   : Baseline end-systolic volume (mL)
MASS0    : 150.0  : Baseline LV wall mass (g)
V0E      : 12.0   : ESPVR volume intercept (mL)
PED0     : 8.00   : Baseline filling pressure (mmHg)
NPED     : 1.30   : Diastolic pressure-volume exponent
KPES     : 0.50   : Incomplete emptying to filling pressure gain
PESSAT   : 1.50   : Saturation of that coupling (unbounded, it destabilises)
KSTIF_COL: 0.50   : Fibrosis stiffening of the diastolic PV relation
SBP0     : 120.0  : Baseline systolic pressure (mmHg)
HR0      : 72.0   : Baseline heart rate (1/min)
EDVMAX   : 450.0  : Numerical ceiling on chamber growth (mL) - NOT biology
BNPMAX   : 20.0   : Cap on the wall-stress ratio driving natriuretic peptide
KDIL     : 2.5e-4 : 1/h eccentric dilation per unit overstress (POSITIVE loop)
KDIL_A   : 0.60   : Angiotensin II amplification of dilation
KDIL_L   : 0.40   : Aldosterone amplification of dilation
KEXP     : 8.0e-4 : 1/h early infarct expansion while the scar is immature
KREVR    : 1.0e-4 : 1/h reverse remodelling below the wall-stress set point
KHYP     : 3.0e-4 : 1/h CONCENTRIC hypertrophy per unit systolic overstress
KHYP_E   : 8.0e-4 : 1/h ECCENTRIC hypertrophy per unit diastolic overstress
MASSMAX  : 300.0  : Ceiling on LV wall mass (g) - hypertrophy is not unlimited
KATR     : 1.2e-4 : 1/h atrophy below the set point
KTH      : 3.0e-3 : 1/h infarct wall thinning
THINMAX  : 0.55   : A healed infarct thins to ~45% of wall, not to nothing
KTHR     : 4.0e-4 : 1/h recovery of infarct thickness once scarred
KST_ON   : 0.55   : Stunning induction by oxidant and Ca2+ load
KST_OFF  : 0.012  : 1/h stunning recovery (tau about 3.5 days)
EMASS_N  : 1.00   : Exponent mapping wall mass to end-systolic elastance
EDIL_N   : 2.00   : Exponent mapping chamber size to elastance (from Laplace)
NE_DESENS: 0.18   : Chronic norepinephrine to beta-receptor downregulation

// ---- neurohormonal ---------------------------------------------------------
KNE      : 0.15   : 1/h sympathetic relaxation
GNE_IS   : 6.00   : Sympathetic drive from infarct size (independent of output)
GANG_IS  : 4.00   : Angiotensin II drive from infarct size
GNE      : 6.00   : Sympathetic gain on cardiac output deficit
KANG     : 0.30   : 1/h angiotensin II relaxation
GANG     : 4.00   : Angiotensin II gain on cardiac output deficit
GANG_NE  : 0.60   : Angiotensin II gain on sympathetic drive
KALD     : 0.20   : 1/h aldosterone relaxation
GALD     : 1.10   : Aldosterone gain on angiotensin II
KBNP     : 0.08   : 1/h natriuretic peptide relaxation
NBNP     : 2.20   : Exponent mapping wall stress to natriuretic peptide
KVOL     : 0.05   : 1/h volume relaxation
GVOL_A   : 0.25   : Aldosterone gain on plasma volume
GVOL_N   : 0.15   : Natriuretic peptide effect on plasma volume
GVOL_S   : 0.10   : SGLT2 inhibitor effect on plasma volume

// ---- biomarkers ------------------------------------------------------------
KTN      : 900.0  : cTnI release per unit necrosis rate
KTND     : 0.035  : 1/h cTnI clearance (t1/2 about 20 h)
KCK      : 700.0  : CK-MB release per unit necrosis rate
KCKD     : 0.058  : 1/h CK-MB clearance (t1/2 about 12 h)
WASH_FLOOR: 0.15  : Marker release fraction achievable without perfusion

// ---- host covariates -------------------------------------------------------
AGE      : 62.0   : Age (years) - reporting covariate only
PRECOND  : 0.0    : Flag - pre-infarction angina (ischemic preconditioning)
PRECOND_EFF: 0.35 : Reduction in mPTP opening rate if preconditioned

// ---- intervention times, hours from symptom onset (1e6 = not given) -------
T_OCC    : 0.0    : Time of coronary occlusion
T_PCI    : 1e6    : Time of mechanical reperfusion (primary or rescue PCI)
T_LYT    : 1e6    : Time of the tenecteplase bolus
TNK_DOSE : 40.0   : Tenecteplase dose (mg)
T_ASA    : 1e6    : Time aspirin is given
T_TIC    : 1e6    : Time the ticagrelor loading dose is given
T_HEP    : 1e6    : Time anticoagulation starts
TIC_LOAD : 180.0  : Ticagrelor loading dose (mg)
TIC_MAINT: 90.0   : Ticagrelor maintenance dose (mg)
TIC_TAU  : 12.0   : Ticagrelor dosing interval (h)
T_MET_IV : 1e6    : Time of IV metoprolol
MET_IV_DOSE: 15.0 : IV metoprolol dose (mg)
T_MET_PO : 1e6    : Time oral metoprolol starts
MET_PO   : 100.0  : Oral metoprolol daily dose (mg/day)
T_RAM    : 1e6    : Time ramipril starts
RAM_PO   : 10.0   : Ramipril daily dose (mg/day)
T_EPL    : 1e6    : Time eplerenone starts
EPL_PO   : 50.0   : Eplerenone daily dose (mg/day)
T_EMP    : 1e6    : Time empagliflozin starts
EMP_PO   : 10.0   : Empagliflozin daily dose (mg/day)
T_COLC   : 1e6    : Time colchicine starts
COLC_PO  : 0.50   : Colchicine daily dose (mg/day)
T_CAN    : 1e6    : Time of the anti-IL-1beta dose
CAN_DOSE : 150.0  : Anti-IL-1beta dose (mg)
T_CSA    : 1e6    : Time of the IV cyclosporine (mPTP inhibitor) dose
CSA_DOSE : 175.0  : Cyclosporine dose (mg)

$CMT @annotated
// ---- layer 1 = subendocardium ... layer 5 = subepicardium -----------------
E1   : Energy charge, layer 1
E2   : Energy charge, layer 2
E3   : Energy charge, layer 3
E4   : Energy charge, layer 4
E5   : Energy charge, layer 5
G1   : Glycogen, layer 1
G2   : Glycogen, layer 2
G3   : Glycogen, layer 3
G4   : Glycogen, layer 4
G5   : Glycogen, layer 5
H1   : Acid load, layer 1
H2   : Acid load, layer 2
H3   : Acid load, layer 3
H4   : Acid load, layer 4
H5   : Acid load, layer 5
C1   : Calcium load, layer 1
C2   : Calcium load, layer 2
C3   : Calcium load, layer 3
C4   : Calcium load, layer 4
C5   : Calcium load, layer 5
P1   : mPTP open fraction, layer 1
P2   : mPTP open fraction, layer 2
P3   : mPTP open fraction, layer 3
P4   : mPTP open fraction, layer 4
P5   : mPTP open fraction, layer 5
NI1  : Ischemic necrosis fraction, layer 1
NI2  : Ischemic necrosis fraction, layer 2
NI3  : Ischemic necrosis fraction, layer 3
NI4  : Ischemic necrosis fraction, layer 4
NI5  : Ischemic necrosis fraction, layer 5
NR1  : Reperfusion necrosis fraction, layer 1
NR2  : Reperfusion necrosis fraction, layer 2
NR3  : Reperfusion necrosis fraction, layer 3
NR4  : Reperfusion necrosis fraction, layer 4
NR5  : Reperfusion necrosis fraction, layer 5
SU1  : Ischemic succinate pool, layer 1
SU2  : Ischemic succinate pool, layer 2
SU3  : Ischemic succinate pool, layer 3
SU4  : Ischemic succinate pool, layer 4
SU5  : Ischemic succinate pool, layer 5
// ---- global ---------------------------------------------------------------
ROS  : Oxidant burst
MVO  : Microvascular obstruction
THR  : Occlusive thrombus burden
PLN  : Plasmin activity
PLG  : Plasminogen, fraction of normal
FIB  : Fibrinogen, fraction of normal
PAI  : PAI-1, fraction of normal
DAMP : Danger-associated molecular patterns
NEU  : Infarct neutrophils
M1   : Pro-inflammatory macrophages
M2   : Reparative macrophages
IL1  : IL-1beta
IL6  : IL-6
CRP  : C-reactive protein
TGF  : TGF-beta
MYOF : Myofibroblasts
COL  : Infarct collagen
MMP  : MMP-2/9 activity
EDVS : Structural (unstressed) chamber volume, mL
MASS : LV wall mass, g
THIN : Infarct wall thinning, 0-1
STUN : Stunning, 0-1
NE   : Sympathetic drive, 1 = normal
ANG  : Angiotensin II, 1 = normal
ALD  : Aldosterone, 1 = normal
BNP  : Natriuretic peptide, 1 = normal
VOL  : Plasma volume, 1 = normal
TNI  : Circulating cTnI
CKMB : Circulating CK-MB
// ---- pharmacokinetic ------------------------------------------------------
TNK1 : Tenecteplase central (mg)
TNK2 : Tenecteplase peripheral (mg)
TICd : Ticagrelor depot (mg)
TICc : Ticagrelor central (mg)
METd : Metoprolol depot (mg)
METc : Metoprolol central (mg)
RAMd : Ramipril depot (mg)
RAMc : Ramipril central (mg)
EPL  : Eplerenone effect compartment (mg)
EMP  : Empagliflozin effect compartment (mg)
COLC : Colchicine effect compartment (mg)
CAN  : Anti-IL-1beta effect compartment (mg)
CSA  : Cyclosporine effect compartment (mg)

$MAIN
// The healthy ventricle must be an EXACT fixed point, so the mechanical set
// points are DERIVED from the baseline geometry rather than typed in.  If they
// were typed in, an uninjured heart would drift and every remodelling result
// below would be contaminated by that drift.
double r_ed0 = pow(3.0 * EDV0 / (4.0 * M_PI), 1.0 / 3.0);
double h_ed0 = MASS0 / (1.05 * 4.0 * M_PI * r_ed0 * r_ed0);
double r_es0 = pow(3.0 * ESV0 / (4.0 * M_PI), 1.0 / 3.0);
double h_es0 = MASS0 / (1.05 * 4.0 * M_PI * r_es0 * r_es0);
double Pes0  = 0.90 * SBP0;

WS_SET    = PED0 * r_ed0 / (2.0 * h_ed0);
WS_SET_ES = Pes0 * r_es0 / (2.0 * h_es0);
EMAX0     = Pes0 / (ESV0 - V0E);
CO0       = (EDV0 - ESV0) * HR0 / 1000.0;

E1_0 = 1.0; E2_0 = 1.0; E3_0 = 1.0; E4_0 = 1.0; E5_0 = 1.0;
G1_0 = 1.0; G2_0 = 1.0; G3_0 = 1.0; G4_0 = 1.0; G5_0 = 1.0;
PLG_0 = 1.0; FIB_0 = 1.0; PAI_0 = 1.0;
NE_0 = 1.0; ANG_0 = 1.0; ALD_0 = 1.0; BNP_0 = 1.0; VOL_0 = 1.0;
EDVS_0 = EDV0; MASS_0 = MASS0;
// The culprit lesion is already thrombosed at t = 0.  Plaque rupture is the
// initial condition, not something the model has to grow.  Whether that
// thrombus obstructs flow is governed separately by T_OCC, so a run with
// T_OCC = 1e6 is a healthy heart with a mural thrombus and normal flow.
THR_0 = 1.0;

$GLOBAL
// Quantities computed in $ODE and reported in $TABLE must live at file scope:
// mrgsolve compiles $ODE and $TABLE into separate function bodies, so a
// variable declared locally in $ODE is NOT visible in $TABLE.
double is_lv, ni_lv, nr_lv;
double EDVv, ESVv, SVv, COv, HRb, SBPb, Ped, Pes;
double ws_ed, ws_es, ws_inf, h_ed, h_es, h_inf;
double contractile, col_n, pat, emax, rpp, qavg;
double NIv[5], NRv[5];
// Derived in $MAIN, used in $ODE and $TABLE -- so also file scope.
double WS_SET, WS_SET_ES, EMAX0, CO0;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define RELU(x)     ((x) > 0.0 ? (x) : 0.0)
#define CLIP01(x)   ((x) < 0.0 ? 0.0 : ((x) > 1.0 ? 1.0 : (x)))
#define HILLF(x,k,n) ((x) <= 0.0 ? 0.0 : pow(x,n) / (pow(k,n) + pow(x,n)))
// Smooth 0->1 step: a tanh rather than an if(), so the stiff solver never has
// to integrate across a true discontinuity.  Width is in hours.
#define SSTEP(t,t0,w) (0.5 * (1.0 + tanh(((t) - (t0)) / (w))))
#define WINDOW(t,t0,dur,w) (SSTEP(t,t0,w) * (1.0 - SSTEP(t,(t0)+(dur),w)))

$ODE
// ===========================================================================
//  A.  drug input rates (mg/h), from time parameters
// ===========================================================================
double r_tnk    = TNK_DOSE / (5.0/60.0) * WINDOW(SOLVERTIME, T_LYT, 5.0/60.0, 0.03);
double r_tic    = TIC_LOAD / (5.0/60.0) * WINDOW(SOLVERTIME, T_TIC, 5.0/60.0, 0.03)
                + TIC_MAINT / TIC_TAU * SSTEP(SOLVERTIME, T_TIC + TIC_TAU, 0.2);
double r_met_iv = MET_IV_DOSE / (10.0/60.0) * WINDOW(SOLVERTIME, T_MET_IV, 10.0/60.0, 0.03);
double r_met_po = MET_PO / 24.0 * SSTEP(SOLVERTIME, T_MET_PO, 0.2);
double r_ram    = RAM_PO  / 24.0 * SSTEP(SOLVERTIME, T_RAM,  0.2);
double r_epl    = EPL_PO  / 24.0 * SSTEP(SOLVERTIME, T_EPL,  0.2);
double r_emp    = EMP_PO  / 24.0 * SSTEP(SOLVERTIME, T_EMP,  0.2);
double r_colc   = COLC_PO / 24.0 * SSTEP(SOLVERTIME, T_COLC, 0.2);
double r_can    = CAN_DOSE / 6.0 * WINDOW(SOLVERTIME, T_CAN, 6.0, 0.2);
double r_csa    = CSA_DOSE / (10.0/60.0) * WINDOW(SOLVERTIME, T_CSA, 10.0/60.0, 0.03);

// ===========================================================================
//  B.  pharmacokinetics
// ===========================================================================
double c_tnk  = RELU(TNK1) / TNK_V1;
double c_tnk2 = RELU(TNK2) / TNK_V2;
dxdt_TNK1 = r_tnk - TNK_CL * c_tnk - TNK_Q * c_tnk + TNK_Q * c_tnk2;
dxdt_TNK2 = TNK_Q * c_tnk - TNK_Q * c_tnk2;

double c_tic = RELU(TICc) / TIC_V;
double c_met = RELU(METc) / MET_V;
double c_ram = RELU(RAMc) / RAM_V;
dxdt_TICd = r_tic * TIC_F - TIC_KA * RELU(TICd);
dxdt_TICc = TIC_KA * RELU(TICd) - TIC_CL * c_tic;
dxdt_METd = r_met_po * MET_F - MET_KA * RELU(METd);
dxdt_METc = r_met_iv + MET_KA * RELU(METd) - MET_CL * c_met;
dxdt_RAMd = r_ram * RAM_F - RAM_KA * RELU(RAMd);
dxdt_RAMc = RAM_KA * RELU(RAMd) - RAM_CL * c_ram;
dxdt_EPL  = r_epl  - EPL_KEL  * RELU(EPL);
dxdt_EMP  = r_emp  - EMP_KEL  * RELU(EMP);
dxdt_COLC = r_colc - COLC_KEL * RELU(COLC);
dxdt_CAN  = r_can  - CAN_KEL  * RELU(CAN);
dxdt_CSA  = r_csa  - CSA_KEL  * RELU(CSA);

// ---- pharmacodynamic occupancies ------------------------------------------
double ipi     = TIC_EMAX * c_tic / (TIC_EC50 + c_tic);
double asa     = ASA_EFF * SSTEP(SOLVERTIME, T_ASA, 0.2);
double hep     = HEP_EFF * SSTEP(SOLVERTIME, T_HEP, 0.2);
double plt_inh = 1.0 - (1.0 - asa) * (1.0 - 0.85 * ipi);
double bb      = MET_EMAX  * c_met / (MET_EC50 + c_met);
double acei    = ACEI_EMAX * c_ram / (RAM_KD + c_ram);
double mra     = MRA_EMAX  * RELU(EPL)  / (EPL_KD  + RELU(EPL));
double sglt    = SGLT_EMAX * RELU(EMP)  / (EMP_KD  + RELU(EMP));
double colc    = COLC_EMAX * RELU(COLC) / (COLC_KD + RELU(COLC));
double can     = CAN_EMAX  * RELU(CAN)  / (CAN_KD  + RELU(CAN));
double csa     = CSA_EMAX  * RELU(CSA)  / (CSA_KD  + RELU(CSA));

// ===========================================================================
//  C.  haemodynamics.  Computed BEFORE the ischemia block, because the
//      rate-pressure product is what sets the wavefront's speed -- that is
//      the entire mechanism by which an IV beta-blocker can shrink an infarct.
// ===========================================================================
double ne_eff  = 1.0 + (RELU(NE) - 1.0) * (1.0 - BB_NE * bb);
double ald_eff = RELU(ALD) * (1.0 - mra);          // MRA blocks the RECEPTOR
double np_eff  = RELU(BNP) * (1.0 + 0.80 * ARNI);

HRb = HR0 * (1.0 + 0.30 * (ne_eff - 1.0)) * (1.0 - BB_HR * bb);
if (HRb < 45.0)  HRb = 45.0;
if (HRb > 160.0) HRb = 160.0;
SBPb = SBP0 * (1.0 + 0.12 * (RELU(ANG) - 1.0) + 0.08 * (ne_eff - 1.0)
                      - 0.10 * acei - 0.06 * ARNI);
if (SBPb < 70.0)  SBPb = 70.0;
if (SBPb > 200.0) SBPb = 200.0;
Pes = 0.90 * SBPb;
rpp = (HRb / HR0) * (SBPb / SBP0);

// ---- epicardial patency ----------------------------------------------------
double occ = SSTEP(SOLVERTIME, T_OCC, 0.02);
pat = 1.0 - occ * CLIP01(THR);
// BASAL metabolism is rate-INDEPENDENT: by definition it is the
// non-contractile cost of ion homeostasis and does not scale with the
// rate-pressure product.  Only CONTRACTILE demand does, which is where rpp
// appears below.  Collateral flow, by contrast, is DIASTOLIC, so it rises when
// the heart slows -- and that, not the basal bill, is how rate control reaches
// the wavefront.
double basal = DEM_BAS * DEM;
double diast = pow(HR0 / (HRb > 1.0 ? HRb : 1.0), KDIAST);

// ---- per-layer flow, supply and AFFORDABLE contraction ---------------------
// The most consequential lines in the model.  A layer spends what it can
// afford on contraction and no more, so a layer whose supply still covers
// basal metabolism goes AKINETIC BUT STAYS ALIVE (hibernating myocardium)
// instead of dying.  That is why a subepicardial rim survives a permanently
// occluded artery, and it is not written anywhere as an assumption.
double GCv[5] = {GC1, GC2, GC3, GC4, GC5};
double GMv[5] = {GM1, GM2, GM3, GM4, GM5};
double Wv[5]  = {W1, W2, W3, W4, W5};
double Ev[5]  = {CLIP01(E1), CLIP01(E2), CLIP01(E3), CLIP01(E4), CLIP01(E5)};
double Gv[5]  = {RELU(G1), RELU(G2), RELU(G3), RELU(G4), RELU(G5)};
double Hv[5]  = {RELU(H1), RELU(H2), RELU(H3), RELU(H4), RELU(H5)};
double Cv[5]  = {RELU(C1), RELU(C2), RELU(C3), RELU(C4), RELU(C5)};
double Pv[5]  = {CLIP01(P1), CLIP01(P2), CLIP01(P3), CLIP01(P4), CLIP01(P5)};
NIv[0] = CLIP01(NI1); NIv[1] = CLIP01(NI2); NIv[2] = CLIP01(NI3);
NIv[3] = CLIP01(NI4); NIv[4] = CLIP01(NI5);
NRv[0] = CLIP01(NR1); NRv[1] = CLIP01(NR2); NRv[2] = CLIP01(NR3);
NRv[3] = CLIP01(NR4); NRv[4] = CLIP01(NR5);
double SUv[5] = {RELU(SU1), RELU(SU2), RELU(SU3), RELU(SU4), RELU(SU5)};

double q_ante[5], qc[5], anaer[5], supply[5], deficit[5], surplus[5], CF[5], VIA[5];
int i;
for (i = 0; i < 5; i++) {
  q_ante[i] = RELU(pat * (1.0 - GMv[i] * CLIP01(MVO)));
  double q  = q_ante[i] + (1.0 - pat) * COLL * GCv[i] * diast;
  qc[i]     = q < 1.0 ? q : 1.0;
  anaer[i]  = F_ANAERO * RELU(1.0 - qc[i]) * Gv[i] / (1.0 + Hv[i] / HI50);
  supply[i] = q + anaer[i];
  deficit[i] = RELU(basal - supply[i]);
  surplus[i] = RELU(supply[i] - basal);
  double afford = CLIP01((supply[i] / (DEM * rpp) - DEM_BAS) / (1.0 - DEM_BAS));
  CF[i]  = afford * Ev[i];
  VIA[i] = RELU(1.0 - NIv[i] - NRv[i]);
}

// ---- infarct size, stunning and contractile capacity ----------------------
is_lv = 0.0; ni_lv = 0.0; nr_lv = 0.0;
double contr_ar = 0.0, ca_avg = 0.0;
for (i = 0; i < 5; i++) {
  is_lv += Wv[i] * (NIv[i] + NRv[i]);
  ni_lv += Wv[i] * NIv[i];
  nr_lv += Wv[i] * NRv[i];
  contr_ar += Wv[i] * VIA[i] * CF[i] * (1.0 - CLIP01(STUN));
  ca_avg += Wv[i] * Cv[i];
}
is_lv *= AAR; ni_lv *= AAR; nr_lv *= AAR;
contractile = (1.0 - AAR) + AAR * contr_ar;

// ---- closed-loop mechanics -------------------------------------------------
// Order of causation matters:  contractile mass -> Emax -> ESV;  ESV
// (incomplete emptying) plus volume status -> filling pressure;  filling
// pressure displaces the STRUCTURAL chamber size along the diastolic PV curve
// to give the OBSERVED EDV.  One state (EDVS) therefore carries both the day-1
// dilation and the six-month dilation without conflating them.
// Scar strength is MATURITY, not amount.  A small infarct lays down less
// collagen in absolute terms but its scar is just as cross-linked, so both the
// half-maximal constant and the saturating cap are scaled by how much of the
// at-risk bed needs a scar.  Without this, a small infarct's normalised scar
// never matures, infarct expansion never switches off, and every infarct
// diverges however small it is.
double inj = is_lv / AAR; if (inj < 0.02) inj = 0.02;
col_n = RELU(COL) / (RELU(COL) + KSC * inj);
double EDVSb0 = RELU(EDVS) > 25.0 ? RELU(EDVS) : 25.0;
emax = EMAX0 * (contractile > 0.02 ? contractile : 0.02)
              * pow(RELU(MASS) / MASS0, EMASS_N)
              * (1.0 - NE_DESENS * (RELU(NE) - 1.0) / 2.0 > 0.0
                 ? (1.0 - NE_DESENS * (RELU(NE) - 1.0) / 2.0 < 1.0
                    ? 1.0 - NE_DESENS * (RELU(NE) - 1.0) / 2.0 : 1.0) : 0.0);
if (emax < 0.15) emax = 0.15;
// The ESPVR volume intercept is NOT fixed: it shifts rightward as the chamber
// remodels.  Without this, end-systolic volume stays small however far the
// ventricle dilates, ejection fraction rises instead of falling, end-systolic
// wall stress never rises, and the concentric-hypertrophy brake never engages
// -- so every infarct diverges.  One term, and it is what gives the
// remodelling bifurcation a stable branch at all.
double v0 = V0E * (EDVSb0 / EDV0);
// From Laplace: generated pressure ~ sigma*h/r ~ sigma*MASS/V, so elastance
// dP/dV ~ sigma*MASS/V^2.  Without this size term, elastance is
// size-independent, end-systolic volume stays small however far the ventricle
// dilates, and ejection fraction RISES during remodelling -- the opposite of
// what happens in patients.
emax = emax / pow(EDVSb0 / EDV0, EDIL_N);
if (emax < 0.15) emax = 0.15;
ESVv = v0 + Pes / emax;
// Filling pressure is set by regulated volume status plus a SATURATING
// contribution from incomplete emptying.  Diastolic stiffening does not raise
// the filling pressure -- it reduces the volume reached at that pressure,
// which is where it belongs (restrictive physiology).
double xes = RELU(ESVv / ESV0 - 1.0);
Ped = PED0 * (RELU(VOL) > 0.3 ? RELU(VOL) : 0.3)
             * (1.0 + KPES * xes / (1.0 + xes / PESSAT));
EDVv = EDVSb0 * pow((Ped / PED0) / (1.0 + KSTIF_COL * col_n), 1.0 / NPED);
if (ESVv > EDVv - 5.0) ESVv = EDVv - 5.0;
SVv = RELU(EDVv - ESVv); if (SVv < 5.0) SVv = 5.0;
COv = SVv * HRb / 1000.0;

double r_ed = pow(3.0 * EDVv / (4.0 * M_PI), 1.0 / 3.0);
h_ed = RELU(MASS) / (1.05 * 4.0 * M_PI * r_ed * r_ed);
double r_es = pow(3.0 * (ESVv > 15.0 ? ESVv : 15.0) / (4.0 * M_PI), 1.0 / 3.0);
h_es = RELU(MASS) / (1.05 * 4.0 * M_PI * r_es * r_es);
h_inf = h_ed * (1.0 - CLIP01(THIN)); if (h_inf < 0.15) h_inf = 0.15;
ws_ed  = Ped * r_ed / (2.0 * h_ed);
ws_es  = Pes * r_es / (2.0 * h_es);
ws_inf = Ped * r_ed / (2.0 * h_inf);

// ===========================================================================
//  D.  epicardial thrombus, lysis and the fibrinolytic system
// ===========================================================================
double lysis = KLYS0 + KLYS_PL * HILLF(RELU(PLN), KPL50, NPL);
// Plasmin activates platelets, so a lytic partly undoes itself.  The whole
// clinical benefit of adding a P2Y12 inhibitor to fibrinolysis comes out of
// this one term (Fitzgerald 1988; CLARITY-TIMI 28).
double plt_drive = (1.0 - plt_inh) * (1.0 + KPLT_PLN * RELU(PLN)) * (1.0 - hep);
double pci = KPCI * SSTEP(SOLVERTIME, T_PCI, 0.05);
dxdt_THR = KREG * plt_drive * (1.0 - CLIP01(THR)) - lysis * CLIP01(THR)
           - pci * RELU(CLIP01(THR) - PCI_FAIL);

double tnk_n = c_tnk / 5.0;
dxdt_PLN = KPG * tnk_n * CLIP01(PLG) / (1.0 + RELU(PAI) / KI_PAI) - KPLD * RELU(PLN);
dxdt_PLG = -KPGD * tnk_n * CLIP01(PLG) + KPGS * (1.0 - CLIP01(PLG));
dxdt_FIB = -KFIBD * RELU(PLN) + KFIBS * (1.0 - RELU(FIB));
dxdt_PAI = KPAID * (1.0 + KPAI_I * RELU(IL6) - RELU(PAI));

// ===========================================================================
//  E.  per-layer ischemia and reperfusion  ==  THE WAVEFRONT
// ===========================================================================
double dNI[5], dNR[5], suc_ox[5];
double dNI_tot = 0.0, dNR_tot = 0.0, reox_flux = 0.0, mvo_endo = 0.0;
qavg = 0.0;
double kon = KP_ON * (1.0 - csa) * (1.0 - PRECOND * PRECOND_EFF);
for (i = 0; i < 5; i++) {
  suc_ox[i] = KSUCO * qc[i] * SUv[i];
  double ph_gate = 1.0 / (1.0 + Hv[i] / HP50);
  double dP = kon * HILLF(Cv[i], CA50, NP_CA) * HILLF(RELU(ROS), R50, NP_R)
              * ph_gate * (1.0 - Pv[i]) - KP_OFF * Pv[i] * Ev[i];
  dNI[i] = KNI * pow(RELU(1.0 - Ev[i]), NNI) * (1.0 + KOSM * Hv[i]) * VIA[i];
  dNR[i] = KNR * Pv[i] * VIA[i];
  dNI_tot += Wv[i] * dNI[i];
  dNR_tot += Wv[i] * dNR[i];
  reox_flux += Wv[i] * VIA[i] * suc_ox[i];
  mvo_endo  += Wv[i] * GMv[i] * pow(RELU(1.0 - Ev[i]), 3.0)
               * (q_ante[i] < 1.0 ? q_ante[i] : 1.0);
  qavg += Wv[i] * qc[i];

  // stash the pore derivative; assigned to the right compartment below
  if (i == 0) dxdt_P1 = dP;
  if (i == 1) dxdt_P2 = dP;
  if (i == 2) dxdt_P3 = dP;
  if (i == 3) dxdt_P4 = dP;
  if (i == 4) dxdt_P5 = dP;
}

// Energy charge falls only when supply cannot cover BASAL demand.  The
// E/(E+EFLOOR) factor is not biology, it is a soft floor: it keeps the fall
// near-linear while there is ATP left and stops the state going negative once
// there is not.
#define DEDT(k) (KE * (-deficit[k] * Ev[k] / (Ev[k] + EFLOOR) \
                       + KREC * surplus[k] * (1.0 - Ev[k])))
dxdt_E1 = DEDT(0); dxdt_E2 = DEDT(1); dxdt_E3 = DEDT(2);
dxdt_E4 = DEDT(3); dxdt_E5 = DEDT(4);

#define DGDT(k) (-KGLY * anaer[k] + KGLYREP * qc[k] * (1.0 - Gv[k]))
dxdt_G1 = DGDT(0); dxdt_G2 = DGDT(1); dxdt_G3 = DGDT(2);
dxdt_G4 = DGDT(3); dxdt_G5 = DGDT(4);

// Acid is made anaerobically and removed ONLY by flow.  This is the whole of
// the pH paradox: the tissue's own protection against the pore is washed away
// by the treatment.
#define DHDT(k) (KH * anaer[k] - KHW * qc[k] * Hv[k])
dxdt_H1 = DHDT(0); dxdt_H2 = DHDT(1); dxdt_H3 = DHDT(2);
dxdt_H4 = DHDT(3); dxdt_H5 = DHDT(4);

// Succinate accumulates while anaerobic and is burned the instant oxygen
// returns.  Its oxidation IS the oxidant burst, which is why the burst is a
// reperfusion event and is exactly zero in a heart that was never ischemic.
#define DSDT(k) (KSUC * anaer[k] - suc_ox[k])
dxdt_SU1 = DSDT(0); dxdt_SU2 = DSDT(1); dxdt_SU3 = DSDT(2);
dxdt_SU4 = DSDT(3); dxdt_SU5 = DSDT(4);

// Calcium is gained as energy fails, plus a surge at reflow proportional to
// accumulated acid (NHE-1 -> Na+ -> reverse-mode NCX).
#define DCDT(k) (KCA * pow(RELU(1.0 - Ev[k]), NCA) + KCA_RP * qc[k] * Hv[k] \
                 - KCAO * Ev[k] * Cv[k])
dxdt_C1 = DCDT(0); dxdt_C2 = DCDT(1); dxdt_C3 = DCDT(2);
dxdt_C4 = DCDT(3); dxdt_C5 = DCDT(4);

// Necrosis by two mechanisms, integrated separately so that the split between
// them is a RESULT rather than an assumption.
dxdt_NI1 = dNI[0]; dxdt_NI2 = dNI[1]; dxdt_NI3 = dNI[2];
dxdt_NI4 = dNI[3]; dxdt_NI5 = dNI[4];
dxdt_NR1 = dNR[0]; dxdt_NR2 = dNR[1]; dxdt_NR3 = dNR[2];
dxdt_NR4 = dNR[3]; dxdt_NR5 = dNR[4];

dxdt_ROS = KROS * reox_flux + KROS_N * RELU(NEU) - KROSD * RELU(ROS);

// ===========================================================================
//  F.  microvascular obstruction, and its positive feedback onto flow
// ===========================================================================
// Lysing or crushing a clot sends some of it downstream, so distal
// embolisation is proportional to the RATE of proximal clot removal.
double embol = KMVO_P * (1.0 - plt_inh) * RELU(-dxdt_THR);
dxdt_MVO = (KMVO_E * mvo_endo + KMVO_N * RELU(NEU) + embol) * (1.0 - CLIP01(MVO))
           - KMVO_R * CLIP01(MVO);

// ===========================================================================
//  G.  inflammation and healing
// ===========================================================================
double il1_act = RELU(IL1) * (1.0 - can);
double inv = 1.0 / (RELU(NEU) + 0.05);
double swch = inv / (KSW50 + inv);   // M1 -> M2 follows neutrophil clearance

dxdt_DAMP = KD1 * (dNI_tot + dNR_tot) - KD2 * RELU(DAMP);
dxdt_NEU  = KN1 * RELU(DAMP) * qavg - KN2 * RELU(NEU);
dxdt_M1   = KM1 * (0.6 * RELU(DAMP) + 0.4 * RELU(NEU))
            - KSW * RELU(M1) * swch - KM1D * RELU(M1);
dxdt_M2   = KSW * RELU(M1) * swch - KM2D * RELU(M2);
dxdt_IL1  = KI1 * (RELU(M1) + 0.5 * RELU(NEU)) * (1.0 - colc) - KI1D * RELU(IL1);
dxdt_IL6  = KI6 * (0.7 * il1_act + 0.3 * RELU(ne_eff - 1.0)) - KI6D * RELU(IL6);
dxdt_CRP  = KC1 * RELU(IL6) - KC2 * RELU(CRP);
dxdt_TGF  = KT1 * RELU(M2) - KT2 * RELU(TGF);
dxdt_MYOF = KMY * RELU(TGF) * (1.0 + 0.4 * RELU(ald_eff - 1.0)) * (1.0 - CLIP01(MYOF))
            - KMYD * CLIP01(MYOF) * (1.0 + 2.0 * col_n);
dxdt_COL  = KCO * CLIP01(MYOF) * RELU(1.0 - RELU(COL) / (COLMAX * inj))
            - KCOD * RELU(MMP) * RELU(COL);
dxdt_MMP  = KMM * (0.5 * RELU(NEU) + 0.3 * RELU(M1) + 0.2 * il1_act)
            - KMMD * RELU(MMP);

// ===========================================================================
//  H.  mechanics  ==  THE BIFURCATION
// ===========================================================================
double over_ed  = RELU(ws_ed / WS_SET - 1.0);
double over_es  = RELU(ws_es / WS_SET_ES - 1.0);
double over_inf = RELU(ws_inf / WS_SET - 1.0);
double under_ed = RELU(1.0 - ws_ed / WS_SET);
double npe = RELU(np_eff - 1.0) / 3.0; if (npe > 1.0) npe = 1.0;
double amp = (1.0 + KDIL_A * RELU(RELU(ANG) - 1.0) + KDIL_L * RELU(ald_eff - 1.0))
             * (1.0 - 0.15 * sglt) * (1.0 - 0.10 * npe);

double EDVSb = RELU(EDVS) > 25.0 ? RELU(EDVS) : 25.0;
// `room` is a numerical ceiling, not a mechanism.  The model contains no
// death and no saturating mechanism, so a divergent trajectory would run to
// infinity and overflow.  A run that reaches EDVMAX should be read as
// "diverged", never as a volume prediction.
double room = RELU(1.0 - EDVSb / EDVMAX);
double unhealed = RELU(1.0 - col_n / SCARSET);
dxdt_EDVS = (KDIL * EDVSb * over_ed * amp
             + KEXP * EDVSb * over_inf * unhealed * is_lv / AAR) * room
            - KREVR * EDVSb * under_ed;
// SCAR CANNOT HYPERTROPHY.  Only surviving myocardium can add mass in
// response to wall stress, so the compensatory brake is scaled by the living
// fraction.  That single factor is what makes the remodelling bifurcation
// exist: below some infarct size the surviving wall can still add mass fast
// enough to hold wall stress at its set point, and above it, it cannot.  No
// critical infarct size is written anywhere -- it is where two rates cross.
double alive = RELU(1.0 - is_lv);
// Hypertrophy has a ceiling.  Once the brake saturates, wall stress keeps
// rising, end-systolic volume rises with it, and ejection fraction FALLS --
// which is what a decompensating ventricle actually does.
double headroom = RELU(1.0 - RELU(MASS) / MASSMAX);
dxdt_MASS = RELU(MASS) * alive * headroom * (KHYP * over_es + KHYP_E * over_ed)
            - KATR * RELU(MASS) * RELU(1.0 - ws_es / WS_SET_ES);
// Bounded thinning: without THINMAX the infarct-zone thickness reaches its
// floor, Laplace stress there becomes arbitrarily large, and infarct expansion
// runs away.  A healed infarct retains roughly 45% of the original wall.
dxdt_THIN = KTH * over_inf * unhealed * (1.0 + 0.5 * RELU(MMP))
            * RELU(THINMAX - CLIP01(THIN)) - KTHR * CLIP01(THIN) * col_n;
dxdt_STUN = KST_ON * (RELU(ROS) + 0.5 * ca_avg) * (1.0 - CLIP01(STUN))
            - KST_OFF * CLIP01(STUN);

// ===========================================================================
//  I.  neurohormonal
// ===========================================================================
double co_n = COv / CO0;
dxdt_NE = KNE * (1.0 + GNE * RELU(1.0 - co_n) + GNE_IS * is_lv - RELU(NE));
double ang_tgt = (1.0 + GANG * RELU(1.0 - co_n) + GANG_IS * is_lv
                  + GANG_NE * RELU(ne_eff - 1.0)) * (1.0 - acei);
dxdt_ANG = KANG * (ang_tgt - RELU(ANG));
dxdt_ALD = KALD * (1.0 + GALD * RELU(RELU(ANG) - 1.0) - RELU(ALD));
double ws_ratio = ws_ed / WS_SET; if (ws_ratio > BNPMAX) ws_ratio = BNPMAX;
dxdt_BNP = KBNP * (pow(ws_ratio, NBNP) - RELU(BNP));
// the natriuretic-peptide effect on plasma volume saturates
double npv = RELU(np_eff - 1.0); if (npv > 1.5) npv = 1.5;
double vol_tgt = 1.0 + GVOL_A * RELU(ald_eff - 1.0) - GVOL_N * npv
                 - GVOL_S * sglt;
if (vol_tgt < 0.6) vol_tgt = 0.6;
dxdt_VOL = KVOL * (vol_tgt - (RELU(VOL) > 0.3 ? RELU(VOL) : 0.3));

// ===========================================================================
//  J.  biomarkers.  Release into blood requires perfusion, which is why a
//      reperfused infarct produces a tall EARLY troponin peak and an occluded
//      one a late blunted peak at a LARGER true infarct size.
// ===========================================================================
double wash = WASH_FLOOR + (1.0 - WASH_FLOOR) * qavg;
double ndot = (dNI_tot + dNR_tot) * AAR;
dxdt_TNI  = KTN * ndot * wash - KTND * RELU(TNI);
dxdt_CKMB = KCK * ndot * wash - KCKD * RELU(CKMB);

$TABLE
// ---- recompute the reported quantities from the same expressions -----------
double IS      = 100.0 * is_lv;
double IS_ISCH = 100.0 * ni_lv;
double IS_REP  = 100.0 * nr_lv;
double RPFRAC  = is_lv > 1e-9 ? 100.0 * nr_lv / is_lv : 0.0;
double SALV    = 100.0 * (AAR - is_lv) / AAR;
double MVO_LV  = 100.0 * CLIP01(MVO) * AAR;
double EDV     = EDVv;
double ESV     = ESVv;
double EF      = 100.0 * SVv / EDVv;
double CO      = COv;
double HR      = HRb;
double SBP     = SBPb;
double WS      = ws_ed;
double WS_INF  = ws_inf;
double PED     = Ped;
double H_ED    = h_ed;
double H_INF   = h_inf;
double CONTR   = contractile;
double SCAR    = col_n;
double THINPCT = 100.0 * CLIP01(THIN);
double STUNPCT = 100.0 * CLIP01(STUN);
double PATENCY = 100.0 * pat;
double NECRO1  = 100.0 * (NIv[0] + NRv[0]);
double NECRO2  = 100.0 * (NIv[1] + NRv[1]);
double NECRO3  = 100.0 * (NIv[2] + NRv[2]);
double NECRO4  = 100.0 * (NIv[3] + NRv[3]);
double NECRO5  = 100.0 * (NIv[4] + NRv[4]);
double CTNI    = RELU(TNI);
double CKMBP   = RELU(CKMB);
double NTPROBNP = RELU(BNP);
double CRPX    = RELU(CRP);
double FIBX    = RELU(FIB);

$CAPTURE @annotated
IS      : Infarct size (% of LV)
IS_ISCH : Infarct size attributable to ischemic necrosis (% of LV)
IS_REP  : Infarct size attributable to reperfusion necrosis (% of LV)
RPFRAC  : Reperfusion-injury fraction of the final infarct (%)
SALV    : Myocardial salvage index (% of the area at risk saved)
MVO_LV  : Microvascular obstruction (% of LV)
NECRO1  : Necrosis in layer 1, subendocardium (% of layer)
NECRO2  : Necrosis in layer 2 (% of layer)
NECRO3  : Necrosis in layer 3, midwall (% of layer)
NECRO4  : Necrosis in layer 4 (% of layer)
NECRO5  : Necrosis in layer 5, subepicardium (% of layer)
EDV     : Observed LV end-diastolic volume (mL)
ESV     : LV end-systolic volume (mL)
EF      : LV ejection fraction (%)
CO      : Cardiac output (L/min)
HR      : Heart rate (1/min)
SBP     : Systolic blood pressure (mmHg)
PED     : LV filling pressure (mmHg)
WS      : Diastolic wall stress, remote zone (arbitrary Laplace units)
WS_INF  : Diastolic wall stress, infarct zone
H_ED    : Diastolic wall thickness (cm)
H_INF   : Infarct-zone wall thickness (cm)
CONTR   : Contractile capacity, fraction of normal
SCAR    : Scar tensile strength, 0-1
THINPCT : Infarct wall thinning (%)
STUNPCT : Stunning (%)
PATENCY : Epicardial patency (%)
CTNI    : Circulating cTnI (arbitrary units)
CKMBP   : Circulating CK-MB (arbitrary units)
NTPROBNP: Natriuretic peptide, multiples of normal
CRPX    : C-reactive protein, multiples of normal
FIBX    : Fibrinogen, fraction of normal

// =============================================================================
//  REFERENCE SCENARIOS
//  ------------------------------------------------------------------------
//  Every scenario is specified purely through param().  The numbers each one
//  should produce are tabulated in README.md and regenerated by
//  ami_reference_check.py (experiment 15).  If R and Python disagree, the
//  Python file is the source of truth and this file has a transcription bug.
//
//  library(mrgsolve); library(dplyr)
//  mod <- mread("ami_mrgsolve_model.R")
//
//  # ---- S1  late presenter, artery never opened, no guideline therapy -----
//  s1 <- mod %>% mrgsim(end = 24*180, delta = 0.05)
//
//  # ---- S2  primary PCI at 90 min + DAPT, no guideline therapy -----------
//  s2 <- mod %>% param(T_PCI = 1.5, T_ASA = 1.0, T_TIC = 1.0, T_HEP = 1.0) %>%
//                mrgsim(end = 24*180, delta = 0.05)
//
//  # ---- S3  primary PCI at 90 min + full GDMT from 24 h ------------------
//  s3 <- mod %>% param(T_PCI = 1.5, T_ASA = 1.0, T_TIC = 1.0, T_HEP = 1.0,
//                      T_RAM = 24, T_MET_PO = 24, T_EPL = 48, T_EMP = 24) %>%
//                mrgsim(end = 24*180, delta = 0.05)
//
//  # ---- S4  prehospital lysis at 45 min, rescue PCI at 4 h, GDMT --------
//  s4 <- mod %>% param(T_LYT = 0.75, T_ASA = 0.65, T_TIC = 0.65, T_HEP = 0.65,
//                      T_PCI = 4.0,
//                      T_RAM = 24, T_MET_PO = 24, T_EPL = 48, T_EMP = 24) %>%
//                mrgsim(end = 24*180, delta = 0.05)
//
//  # ---- S5  PCI 90 min + GDMT + IV metoprolol + CsA at reflow + colchicine
//  s5 <- mod %>% param(T_PCI = 1.5, T_ASA = 1.0, T_TIC = 1.0, T_HEP = 1.0,
//                      T_MET_IV = 0.6, T_CSA = 1.45, T_COLC = 24,
//                      T_RAM = 24, T_MET_PO = 24, T_EPL = 48, T_EMP = 24) %>%
//                mrgsim(end = 24*180, delta = 0.05)
//
//  # ---- S6  large area at risk (0.45), PCI at 4 h, full GDMT ------------
//  s6 <- mod %>% param(AAR = 0.45, T_PCI = 4.0, T_ASA = 3.5, T_TIC = 3.5,
//                      T_HEP = 3.5,
//                      T_RAM = 24, T_MET_PO = 24, T_EPL = 48, T_EMP = 24) %>%
//                mrgsim(end = 24*180, delta = 0.05)
//
//  # ---- S7  small area at risk (0.18), PCI at 3 h, full GDMT -----------
//  s7 <- mod %>% param(AAR = 0.18, T_PCI = 3.0, T_ASA = 2.5, T_TIC = 2.5,
//                      T_HEP = 2.5,
//                      T_RAM = 24, T_MET_PO = 24, T_EPL = 48, T_EMP = 24) %>%
//                mrgsim(end = 24*180, delta = 0.05)
//
//  # ---- the salvage curve: sweep reperfusion time only ------------------
//  idata <- data.frame(ID = 1:14,
//                      T_PCI = c(15,30,45,60,75,90,120,150,180,240,300,360,480,720)/60)
//  sweep <- mod %>% idata_set(idata) %>% param(T_ASA = 0.05, T_TIC = 0.05,
//                                              T_HEP = 0.05) %>%
//                   mrgsim(end = 48, delta = 0.05)
//
//  # ---- the remodelling bifurcation: sweep the area at risk -------------
//  ib <- data.frame(ID = 1:11,
//                   AAR = c(0.10,0.15,0.20,0.25,0.28,0.30,0.32,0.35,0.40,0.45,0.50))
//  bif <- mod %>% idata_set(ib) %>% param(T_PCI = 6, T_ASA = 5.5, T_TIC = 5.5,
//                                         T_HEP = 5.5) %>%
//                 mrgsim(end = 24*365, delta = 1)
//
//  # ---- structural sensitivity: which loop is load-bearing? -------------
//  #   KP_ON = 0     no mitochondrial permeability transition pore
//  #   KCA_RP = 0    no reflow calcium surge
//  #   HP50 = 1e6    acid gate removed (the pore stops being pH-gated)
//  #   KSUC = 0      no ischemic succinate, hence no reperfusion oxidant burst
//  #   GC1..GC5 = 1  flat collateral gradient (the wavefront stops being transmural)
//  #   COLL = 0      no collateral flow at all
//  #   GM1..GM5 = 0  no microvascular-obstruction feedback onto flow
//  #   KDIL = 0      no dilation loop (every infarct becomes benign at 6 months)
//  #   KHYP = 0      no hypertrophy brake
//  #   PRECOND = 1   preconditioned host
// =============================================================================
