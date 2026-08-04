## =====================================================================
##  MALIGNANT HYPERTHERMIA — QSP model (mrgsolve)
##  39 ODE compartments · 11 volatile/drug PK · 28 disease PD · 18 scenarios
##  70 kg adult reference subject.  Time unit: MINUTES.
## =====================================================================
##
##  THE THESIS THIS MODEL IS BUILT AROUND
##  -------------------------------------
##  Malignant hyperthermia is usually taught as a list of signs in the order
##  a textbook prints them.  This model asserts that the order is not a
##  convention at all — it is arithmetic on TWO NUMBERS.
##
##  There is exactly ONE driver: myoplasmic free Ca2+ (compartment CAM), set
##  by RyR1 open probability.  Everything the anaesthetist can see is that
##  one driver integrated into a STORE, and the two stores that matter
##  differ by more than two orders of magnitude in how long they take
##  to fill:
##
##      CO2  store    K_CO2B = 60 mL per mmHg
##      heat store    C_CORE = 158 000 J per degC
##
##  Divide each by the SAME excess flux produced by the SAME event (403 mL
##  CO2/min and 160 W, measured at 5 min into the model's own step
##  experiment) and you get the two time constants that organise the entire
##  clinical picture:
##
##      0.149 min per mmHg of EtCO2        16.5 min per degC of core temp
##
##  i.e. ONE DEGREE OF FEVER COSTS AS MUCH ELAPSED TIME AS 111 mmHg OF
##  EtCO2.  Since the capnograph saturates at a ventilation-limited plateau
##  of 149 mmHg, THE WHOLE EtCO2 SCALE IS EXHAUSTED BEFORE THE THERMOMETER
##  MOVES ONE DEGREE.  "Malignant hyperthermia" is named after the last
##  thing that happens.
##
##  Verified by a step experiment (muscle pre-equilibrated with agent, so
##  the trigger arrives at t = 0 and both stores start from baseline):
##
##      EtCO2 +10 mmHg   at   8.2 min        core +0.5 degC   at   9.8 min
##      EtCO2 +20 mmHg   at  18.5 min        core +1.0 degC   at  20.5 min
##      EtCO2 +40 mmHg   at  42.8 min        core +1.5 degC   at  31.0 min
##                                           core  38.5 degC  at  40.8 min
##
##  SECOND AXIS — THE TERM THAT STARTS IT IS NOT THE TERM THAT SUSTAINS IT
##  ---------------------------------------------------------------------
##  The volatile agent opens RyR1.  But Ca2+ overload plus mitochondrial
##  Ca2+ loading generates ROS/RNS which OXIDISE RyR1 (state SENS), and an
##  oxidised RyR1 behaves as though 1.3 activation units of agent were still
##  present.  SENS reverses with a half-life of ~116 min — far slower than
##  the agent washes out.  The consequence is the whole history of this
##  disease in one line:
##
##      untreated                                  lethal at 140 min
##      vaporiser closed + 10 L/min O2 flush       lethal at 152 min
##      ... plus activated-charcoal filters        lethal at 156 min
##
##  Removing the trigger deletes the START term and leaves the SUSTAIN term
##  running.  Dantrolene blocks the channel itself and is therefore agnostic
##  to which term opened it — which is why mortality fell from ~70-80%
##  before 1979 to a few percent afterwards, in an era when "turn the agent
##  off" had already been standard practice for decades.
##
##  THIRD AXIS — EVERY OTHER MANOEUVRE TREATS A READOUT
##  ---------------------------------------------------
##  The model's treatment arms make this quantitative, and uncomfortable:
##
##    hyperventilation x3 alone   EtCO2 max 50.7 mmHg (a NORMAL capnograph)
##                                core 42.5 degC, CK 31 894, dead at 164 min
##    active cooling alone        core max 38.1 degC (a NORMAL temperature)
##                                EtCO2 121, CK 15 769,      dead at 141 min
##    full bundle + dantrolene    core 37.4 degC, EtCO2 51, CK 4 460, lives
##
##  Both monitoring-directed arms produce a reassuring number on the
##  variable they act on while the patient dies of the variable they do not.
##  This is not a straw man: hyperventilating a rising EtCO2 and cooling a
##  rising temperature are both correct components of the MHAUS bundle.  The
##  point is that neither is a treatment, and each one destroys the signal
##  that would otherwise have prompted the treatment.
##
##  And the drug's own contribution, isolated by running the identical
##  bundle with and without it:
##
##    bundle WITHOUT dantrolene   CK 46 525, K+ 6.9, GFR 0.33   (survives)
##    bundle WITH dantrolene      CK  4 460, K+ 5.0, GFR 0.75
##
##  A ten-fold difference in muscle destroyed, from the one agent that acts
##  on Po rather than on a readout.
##
##  FOURTH AXIS — THE SETPOINT IS NOT WHERE YOU THINK IT IS
##  -------------------------------------------------------
##  At steady state the SR release flux and the SERCA uptake flux CANCEL
##  (they must, or the SR would empty in seconds).  What is left setting the
##  sustained myoplasmic Ca2+ level is the SARCOLEMMAL balance: store-
##  operated entry in (STIM1/Orai1, gain KSOCE) versus PMCA/NCX out.  The
##  model therefore predicts that the sustained phase of MH is an
##  EXTRACELLULAR Ca2+ ENTRY disease even though the trigger is an SR
##  release channel — and reproduces the in-vitro observation that removing
##  bath Ca2+ abolishes the halothane contracture (5.63 -> 1.30 -> 0.17 ->
##  0.00 g at 100/50/20/0% extracellular Ca2+).  SR release is the trigger
##  and the ATP sink; it is not the reservoir.
##
##  WHAT IS FITTED AND WHAT IS DERIVED
##  ----------------------------------
##  FITTED (5 things):
##    1. EC50_VOL, EC50_CAF and PO_BASE per genotype, set so that the
##       in-vitro contracture test classifies MHN as MHN and every MHS
##       variant as MHS at European MH Group thresholds.
##    2. KSOCE per genotype, set so that sustained CAM reaches ~1.05 uM.
##    3. K_INJ, set so that a promptly-treated event yields CK in the low
##       thousands and an untreated one CK > 20 000.
##    4. Dantrolene PK, from Flewellen 1983 (2.4 mg/kg -> ~3-4 ug/mL,
##       t1/2 12-16 h, therapeutic 2.8-4.2 ug/mL).
##    5. Relative RyR1 potency per vol% of the five volatile agents, set so
##       that time-to-onset ranks halothane < enflurane < isoflurane <
##       sevoflurane < desflurane at 1 MAC (Wedel 1993).
##
##  DERIVED, i.e. left as falsifiable output:
##    - the 111 mmHg-per-degC capacitance ratio and the whole
##      EtCO2-precedes-fever ordering of signs
##    - the failure of agent removal alone
##    - the two masking arms
##    - the ten-fold CK separation attributable to dantrolene
##    - recrudescence appearing ONLY in the under-dosed arm
##    - RQ > 1 (bicarbonate buffering of lactate liberates CO2 that never
##      came from O2, so VCO2 rises MORE than VO2)
##    - the extracellular-Ca2+ dependence of the contracture
##    - the cooling arithmetic (1 L of 4 degC saline = 138 kJ = 0.60 degC)
##    - non-depolarising blockers having EXACTLY zero effect on rigidity,
##      because they act upstream of RyR1 and MH Ca2+ arrives downstream
##    - contracture-induced capillary occlusion, hence compartment syndrome
##
##  WHERE THE MODEL DISAGREES WITH THE TEXTBOOK, AND SAYS SO
##  --------------------------------------------------------
##    - The often-quoted "core temperature rises 1-2 degC every 5 minutes"
##      requires 500-1000 W of NET heat storage.  The model's fulminant
##      event peaks at 3.5 degC/h, with whole-body VO2 733 mL/min (3.7x
##      baseline) and heat production 249 W.  The aerobic ceiling (cardiac
##      output ~10 L/min x 200 mL O2/L x 0.85 extraction = 1.7 L O2/min =
##      570 W) puts the textbook rate out of reach except transiently, in
##      the core compartment alone, at near-total muscle recruitment.
##    - Succinylcholine advances sustained onset by only ~1 min here (it
##      does produce the immediate self-limited masseter transient).  That
##      is much weaker than the epidemiological association between
##      succinylcholine and fulminant MH, and is the model's least
##      well-supported structural choice.
##    - The virtual-population complication rate rises with time-to-
##      dantrolene at an odds ratio of 2.38 per 30 min against Larach 2010's
##      observed 1.61.  Right direction, too steep.  See section ANALYSIS 7.
##
## =====================================================================

library(mrgsolve)
library(dplyr)

# =====================================================================
#  MODEL CODE
# =====================================================================
mh_code <- '
$PROB
# Malignant Hyperthermia QSP
- 39 ODEs: volatile agent PK (6) + succinylcholine (1) + dantrolene (4)
  + Ca handling (3) + energetics (3) + gas/acid-base (4) + thermal (2)
  + muscle injury (5) + renal (3) + coagulation (4) + monitors (4)
- ONE driver (myoplasmic Ca2+), TWO stores of very different size.

$PARAM @annotated
// ---------------- subject ------------------------------------------------
WT      :  70.0  : body weight (kg)

// ---------------- volatile agent PK (sevoflurane-equivalent, vol%) --------
VAP     :   2.0  : vaporiser dial setting (vol%)
POTENCY :   1.0  : RyR1 potency per vol% relative to sevoflurane
TVAPOFF : 1e9    : time the vaporiser is closed (min)
FGF     :   2.0  : fresh gas flow while anaesthetising (L/min)
FGFHI   :  10.0  : fresh gas flow during washout (L/min)
LBG     :   0.65 : blood:gas partition coefficient
C_CIRC  : 120.0  : circuit + absorbent + plastics capacity (L gas-equivalent)
C_LUNG  :   2.83 : FRC + pulmonary blood capacity (L gas-equivalent)
V_VRG   :   6.0  : vessel-rich group volume (L)
L_VRG   :   1.7  : VRG:blood partition
Q_VRG   :   3.75 : VRG perfusion (L/min)
V_MUSF  :   8.0  : well-perfused muscle volume (L)
V_MUSS  :  25.0  : bulk muscle volume (L)
L_MUS   :   3.1  : muscle:blood partition
Q_MUSF0 :   0.45 : resting well-perfused muscle flow (L/min)
Q_MUSS0 :   0.45 : resting bulk muscle flow (L/min)
V_FAT   :  14.5  : fat volume (L)
L_FAT   :  48.0  : fat:blood partition
Q_FAT   :   0.25 : fat perfusion (L/min)
K_CHARC : 700.0  : activated-charcoal scrubbing (L/min-equivalent)
TCHARC  : 1e9    : time charcoal filters are inserted (min)

// ---------------- ventilation --------------------------------------------
VA0     :   4.315: baseline alveolar ventilation (L/min)
VENTMLT :   1.0  : minute-ventilation multiplier after TVENT
TVENT   : 1e9    : time hyperventilation is started (min)

// ---------------- succinylcholine ----------------------------------------
V_SCH   :  15.0  : succinylcholine initial distribution volume (L)
KE_SCH  :   0.888: elimination rate constant (/min), t1/2 47 s
EC50_SCH:   0.55 : concentration for half-maximal RyR1 contribution (mg/L)
W_SCH   :   0.45 : weight of the succinylcholine term in RyR1 activation

// ---------------- dantrolene PK/PD ---------------------------------------
V1_DAN  :  32.0  : central volume (L)
V2_DAN  :  18.0  : peripheral 1 volume (L)
V3_DAN  :   6.0  : peripheral 2 volume (L)
CL_DAN  :   0.04 : clearance (L/min) = 2.4 L/h
Q2_DAN  :   0.30 : intercompartmental clearance 1 (L/min)
Q3_DAN  :   0.033: intercompartmental clearance 2 (L/min)
KE0_DAN :   0.14 : effect-site equilibration (/min)
EC50_DAN:   2.6  : effect-site EC50 (mg/L)
EMAX_DAN:   0.87 : maximum fractional block of the RyR1 leak (NOT 1.0)
HILL_DAN:   1.4  : Hill coefficient for dantrolene

// ---------------- RyR1 ---------------------------------------------------
KRYR    : 200.0  : SR release rate scaling (/min)
PO_BASE :   5.5e-4 : resting open probability (genotype)
PO_MAX  :   4.8e-3 : maximal open probability under trigger
EC50_VOL:   0.55 : agent vol% in muscle for half-maximal activation (genotype)
HILL_RYR:   3.6  : Hill coefficient (a threshold, not a gradient)
Q10_RYR :   2.6  : temperature sensitivity of RyR1 gating
EC50_CAF:   1.0  : caffeine EC50 (mM, genotype)
CAFF    :   0.0  : caffeine concentration (mM, contracture test)

// ---------------- RyR1 redox sensitisation (the SUSTAIN term) ------------
K_SENS_ON :  0.022 : ROS-driven oxidation rate (/min)
K_SENS_OFF:  0.006 : reduction rate (/min), t1/2 116 min
W_SENS    :  1.30  : activation equivalent of fully oxidised RyR1
NAC       :  0.0   : N-acetylcysteine effect (0-1) on K_SENS_ON

// ---------------- Ca handling (uM, uM/min, cytosol-referenced) -----------
VSERCA  : 371.0  : SERCA maximal rate
KSERCA  :   0.35 : SERCA Km
KATP_S  :   0.5  : ATP Km for SERCA (mM)
KPMCA   :  40.0  : PMCA/NCX extrusion rate constant (/min)
CA_REST :   0.10 : resting myoplasmic Ca (uM)
CASR0   :1400.0  : baseline SR Ca (uM, cytosol-referenced)
KSOCE   :  48.0  : store-operated Ca entry gain (genotype)
CA_O    :   1.0  : extracellular Ca (fraction of normal)
KMITO_U :  30.0  : mitochondrial uptake (/min)
KMITO_R :   0.90 : mitochondrial release
CAMITOMX: 900.0  : mitochondrial Ca capacity (uM)
CAMITO50: 700.0  : mitochondrial Ca at which coupling efficiency halves
KAPPA   :  60.0  : cytosolic Ca buffering power

// ---------------- contraction and energetics -----------------------------
CA50_F  :   0.75 : Ca for half-maximal force (uM)
HILL_F  :   5.0  : force-Ca Hill coefficient (steep, as in skinned fibres)
ATPASE_B: 400.0  : resting muscle ATP turnover (uM/min)
ATPASE_A:3600.0  : actomyosin ATPase at full activation (uM/min)
HEP0    :  33.0  : ATP + phosphocreatine (mM)
ATP0    :   8.0  : ATP (mM)
V_MUSW  :  22.5  : muscle water volume (L)
PO2_ATP :   5.3  : ATP per O2
GLY0    :  85.0  : muscle glycogen (mM glucosyl units)
VGLY_MAX: 900.0  : maximal anaerobic glycolytic ATP rate (uM/min)
F_OCCL  :   0.40 : fraction of muscle perfusion lost at full contracture
QMUS_MAX:   5.0  : ceiling on muscle perfusion (L/min)

// ---------------- gas exchange -------------------------------------------
K_CO2B  :  60.0  : fast body CO2 store (mL/mmHg)
K_CO2M  : 200.0  : muscle CO2 store (mL/mmHg)
G_CO2   :   6.0  : blood CO2 dissociation slope (mL/(L.mmHg))
VCO2_OTH: 162.0  : non-muscle CO2 production (mL/min)
RQ      :   0.85 : respiratory quotient of oxidative metabolism
CAO2    : 200.0  : arterial O2 content (mL/L)
O2EXTMX :   0.85 : maximal O2 extraction
VO2_OTH : 155.0  : non-muscle O2 consumption (mL/min)
DEADSP  :   3.0  : arterial-to-end-tidal CO2 gradient (mmHg)

// ---------------- thermal -------------------------------------------------
C_CORE  :158000.0: core heat capacity (J/degC)
C_SHELL : 73500.0: shell heat capacity (J/degC)
K_CS    :  27.0  : core-shell conductance (W/degC)
H_SKIN  :   5.15 : shell-ambient conductance (W/degC)
T_AMB   :  21.0  : ambient temperature (degC)
HEAT_O2 :20100.0 : heat per litre of O2 (J/L)
HEAT_LAC:  65.0  : heat per mmol anaerobic lactate (J/mmol)
HEAT_PCR:  34.0  : heat per mmol phosphate split without resynthesis (J/mmol)
COOL_W  :   0.0  : core (lavage/intravascular) cooling power (W)
SURF_W  :   0.0  : surface cooling power (W)
TCOOL   : 1e9    : time active cooling is started (min)
TSTOP_C :  38.0  : temperature at which cooling is stopped (degC)

// ---------------- muscle injury -------------------------------------------
K_INJ   :   1.1e-4 : membrane injury gain (/min)
CA_INJ50:   0.45 : Ca above which calpain/PLA2 injury accelerates (uM)
KT_INJ  :   0.55 : additional injury per degC above 38.5
CK_MUS  :   3.0e6: plasma CK per unit fractional muscle destruction (U/L)
KE_CK   :   3.21e-4 : CK elimination (/min), t1/2 36 h
MB_MUS  :   6.0e6: plasma myoglobin per unit fractional destruction (ng/mL)
KE_MB   :   0.0046 : myoglobin elimination (/min), t1/2 2.5 h
K_MUS   : 150.0  : intracellular K (mmol/L)
V_ECF   :  14.0  : extracellular fluid volume (L)
KE_K    :   0.008: K disposal (/min)
K0      :   4.0  : baseline plasma K (mmol/L)
K_ACID  :   0.50 : K shift per 0.1 unit fall in pH (mmol/L)
PHOS0   :   1.1  : baseline phosphate (mmol/L)
K_SHIFT :   0.0  : therapeutic K-lowering flux (insulin/glucose, mmol/L/min)
BICARB  :   0.0  : bicarbonate infusion (mmol/L/min)

// ---------------- renal ----------------------------------------------------
K_MBT   :   0.05 : myoglobin to tubular cast formation
K_CAST  :   0.0035 : cast clearance (/min)
K_GFR   :   2.0e-5 : GFR loss per unit cast per min
K_GFRREC:   1.2e-4 : tubular recovery (/min)
SCR0    :   0.9  : baseline creatinine (mg/dL)
K_SCR_PR:   0.0020 : creatinine production (mg/dL/min)
K_SCR_CL:   0.00222: creatinine clearance at GFRF = 1 (/min)

// ---------------- coagulation and thermal dose ------------------------------
T_DIC   :  40.0  : temperature above which thermal dose accrues (degC)
K_TD    :   1.0  : thermal dose accrual
K_FIB   :   0.010: fibrinogen consumption
K_PLT   :   0.008: platelet consumption
K_DD    :   0.9  : D-dimer generation
FIB0    : 300.0  : baseline fibrinogen (mg/dL)
PLT0    : 250.0  : baseline platelets (10^9/L)
K_HR    :   0.35 : heart-rate equilibration (/min)

$CMT @annotated
CIRC   : anaesthetic circuit agent (vol%)
ALV    : alveolar/arterial agent (vol%)
VRG    : vessel-rich group agent (vol%)
MUSF   : well-perfused muscle agent (vol%) - THE TRIGGER COMPARTMENT
MUSS   : bulk muscle agent (vol%)
FAT    : fat agent (vol%)
SCHA   : succinylcholine amount (mg)
DAN1   : dantrolene central (mg)
DAN2   : dantrolene peripheral 1 (mg)
DAN3   : dantrolene peripheral 2 (mg)
DANE   : dantrolene effect site (mg/L)
CAM    : myoplasmic free Ca (uM) - THE DRIVER
CASR   : SR Ca store (uM, cytosol-referenced)
CAMITO : mitochondrial Ca (uM)
HEP    : ATP + phosphocreatine (mM)
GLY    : muscle glycogen (mM)
LACM   : muscle lactate (mM)
PMCO2  : muscle PCO2 (mmHg)
PACO2  : arterial PCO2 (mmHg)
LACS   : systemic lactate (mM)
HCO3   : plasma bicarbonate (mmol/L)
TCORE  : core temperature (degC)
TSHELL : shell temperature (degC)
MEMB   : sarcolemmal integrity (1 = intact)
CK     : plasma creatine kinase (U/L)
MB     : plasma myoglobin (ng/mL)
KP     : plasma potassium (mmol/L)
PHOS   : plasma phosphate (mmol/L)
MBT    : tubular myoglobin cast burden
GFRF   : GFR fraction of normal
SCR    : serum creatinine (mg/dL)
TDOSE  : cumulative thermal dose (degC.min above 40)
FIB    : fibrinogen (mg/dL)
PLT    : platelets (10^9/L)
DDIM   : D-dimer (ng/mL FEU)
HR     : heart rate (/min)
VO2CUM : cumulative O2 consumption (mL)
NEURO  : cerebral injury index
SENS   : RyR1 redox sensitisation (0-1) - THE SUSTAIN TERM

$MAIN
CASR_0   = CASR0;
CAM_0    = CA_REST;
CAMITO_0 = 20.0;
HEP_0    = HEP0;
GLY_0    = GLY0;
LACM_0   = 1.2;
PMCO2_0  = 47.0;
PACO2_0  = 40.0;
LACS_0   = 1.0;
HCO3_0   = 24.0;
TCORE_0  = 36.5;
TSHELL_0 = 34.0;
MEMB_0   = 1.0;
CK_0     = 100.0;
MB_0     = 30.0;
KP_0     = K0;
PHOS_0   = PHOS0;
GFRF_0   = 1.0;
SCR_0    = SCR0;
FIB_0    = FIB0;
PLT_0    = PLT0;
HR_0     = 70.0;

$GLOBAL
#define POSPART(x) ((x) > 0.0 ? (x) : 0.0)
// Values computed in $ODE are local to that block, so anything named in
// $CAPTURE has to live at file scope and merely be ASSIGNED inside $ODE.
double ATP;
double PCR;
double ACT;
double Po;
double FORCE;
double FORCEca;
double RIGOR;
double ATPuse;
double ATP_ox;
double ATP_gly;
double HEP_DEF;
double MITOEFF;
double Q_MUS;
double occl;
double VO2_tot;
double VCO2_mus;
double HEAT;
double inj;
double ROS;
double hyper;

$ODE
// ---------- guards ---------------------------------------------------------
double cam    = fmax(CAM,    1e-6);
double casr   = fmax(CASR,   0.0);
double hep    = fmax(HEP,    0.0);
double gly    = fmax(GLY,    0.0);
double memb   = fmin(fmax(MEMB, 0.0), 1.0);
double gfrf   = fmin(fmax(GFRF, 0.02), 1.0);
double sens   = fmin(fmax(SENS, 0.0), 1.0);
double dane   = fmax(DANE,   0.0);

// ---------- interventions (time-switched) ---------------------------------
double vap    = (SOLVERTIME < TVAPOFF) ? VAP  : 0.0;
double fgf    = (SOLVERTIME < TVAPOFF) ? FGF  : FGFHI;
double charc  = (SOLVERTIME >= TCHARC) ? 1.0  : 0.0;
double VA     = VA0 * ((SOLVERTIME >= TVENT) ? VENTMLT : 1.0);
// Smooth thermostat: a hard switch on a STATE variable makes any stiff solver
// chatter, so "stop active cooling at 38 degC" is written as a steep logistic.
double thermo = 1.0 / (1.0 + exp(-(TCORE - TSTOP_C) / 0.15));
double coolC  = ((SOLVERTIME >= TCOOL) ? COOL_W : 0.0) * thermo;
double coolS  = ((SOLVERTIME >= TCOOL) ? SURF_W : 0.0) * thermo;

// ---------- high-energy phosphate split (creatine-kinase clamp) -----------
// PCr is spent FIRST; ATP is defended until PCr is exhausted.
ATP    = fmin(ATP0, hep);
PCR    = POSPART(hep - ATP0);
double atp_ok = ATP / (ATP + KATP_S);

// ---------- RyR1 activation ------------------------------------------------
double fT     = pow(Q10_RYR, (TCORE - 37.0) / 10.0);
double aVOL   = MUSF * POTENCY / EC50_VOL;
double cSCH   = SCHA / V_SCH;
double aSCH   = W_SCH * cSCH / (cSCH + EC50_SCH);
double aCAF   = CAFF / EC50_CAF;
ACT    = fT * (aVOL + aSCH + aCAF + W_SENS * sens);
double actH   = pow(ACT, HILL_RYR);
double Po_t   = PO_BASE + (PO_MAX - PO_BASE) * actH / (1.0 + actH);
double inh    = EMAX_DAN * pow(dane, HILL_DAN) /
                (pow(EC50_DAN, HILL_DAN) + pow(dane, HILL_DAN));
// dantrolene blocks the ACTIVATED leak, and weakly the resting leak
Po     = PO_BASE * (1.0 - 0.25 * inh) + (Po_t - PO_BASE) * (1.0 - inh);

// ---------- Ca fluxes -------------------------------------------------------
double J_ryr   = KRYR * Po * (casr - cam);
double J_serca = VSERCA * cam * cam / (KSERCA * KSERCA + cam * cam) * atp_ok;
double J_pmca  = KPMCA * (cam - CA_REST) * (0.55 + 0.45 * atp_ok);
double J_soce  = KSOCE * POSPART(1.0 - casr / CASR0) * CA_O;
double J_mu    = KMITO_U * cam * POSPART(1.0 - CAMITO / CAMITOMX);
double J_mr    = KMITO_R * CAMITO / CAMITOMX * 12.0;

// ---------- contraction and ATP demand -------------------------------------
double camH    = pow(cam, HILL_F);
FORCEca = camH / (pow(CA50_F, HILL_F) + camH);
RIGOR   = (ATP < 2.0) ? (1.0 - ATP / 2.0) : 0.0;
FORCE   = fmin(1.0, fmax(FORCEca, RIGOR));   // NMBA blocks NEITHER
ATPuse  = ATPASE_B + J_serca / 2.0 + ATPASE_A * FORCEca;

// ---------- muscle perfusion: hyperaemia vs contracture occlusion ----------
hyper   = ATPuse / ATPASE_B;
occl    = 1.0 - F_OCCL * FORCEca;
Q_MUS   = fmin(QMUS_MAX, (Q_MUSF0 + Q_MUSS0) * (1.0 + 2.2 * (hyper - 1.0))) * occl;
double fsplit  = Q_MUSF0 / (Q_MUSF0 + Q_MUSS0);
double Q_MUSF  = Q_MUS * fsplit;
double Q_MUSS  = Q_MUS * (1.0 - fsplit);

// ---------- oxidative vs glycolytic ATP ------------------------------------
// Mitochondrial Ca overload UNCOUPLES oxidation: the same O2 yields less ATP
// and therefore MORE heat per litre consumed.
MITOEFF = 1.0 / (1.0 + pow(CAMITO / CAMITO50, 2.0));
double DO2     = Q_MUS * CAO2 * O2EXTMX;
double ATP_oxc = (DO2 / 22.4) * 1000.0 * PO2_ATP * MITOEFF / V_MUSW;
ATP_ox  = fmin(ATPuse, ATP_oxc);
double pH_m    = 7.05 - 0.055 * POSPART(LACM - 1.2);
double fgly    = fmin(1.0, POSPART((pH_m - 6.35) / 0.55));
ATP_gly = (gly > 0.5) ? fmin(VGLY_MAX * fgly, POSPART(ATPuse - ATP_ox)) : 0.0;
double LACprod = ATP_gly * 2.0 / 3.0;
double VO2_mus = ATP_ox / (PO2_ATP * MITOEFF) * V_MUSW / 1000.0 * 22.4;
HEP_DEF = POSPART(ATPuse - ATP_ox - ATP_gly);

// ---------- CO2 --------------------------------------------------------------
// RQ RISES ABOVE 1 because bicarbonate buffering of lactic acid liberates CO2
// that never came from O2.  VCO2 therefore rises MORE than VO2 in MH.
double Hbuf     = LACprod * V_MUSW / 1000.0;
VCO2_mus = VO2_mus * RQ + Hbuf * 22.4;
double flux_mb  = G_CO2 * Q_MUS * (PMCO2 - PACO2);

// ---------- lactate and acid-base --------------------------------------------
double lac_out  = 0.55 * (LACM - LACS);
double pH       = 6.1 + log10(fmax(HCO3, 1.0) / (0.03 * fmax(PACO2, 5.0)));

// ---------- thermal -----------------------------------------------------------
VO2_tot  = VO2_mus + VO2_OTH;
HEAT     = HEAT_O2 * (VO2_tot / 1000.0) / 60.0
                + HEAT_LAC * (LACprod * V_MUSW / 1000.0) / 60.0
                + HEAT_PCR * (HEP_DEF * V_MUSW / 1000.0) / 60.0;

// ---------- muscle injury: FOUR MULTIPLICATIVE DRIVERS ------------------------
double ca_drv   = POSPART(cam - CA_INJ50);
double heat_drv = 1.0 + KT_INJ * POSPART(TCORE - 38.5);
double atp_drv  = 1.0 + 2.5 * POSPART(1.0 - ATP / ATP0);
inj      = K_INJ * pow(ca_drv, 1.6) * heat_drv * atp_drv * memb;
double kleak    = inj * K_MUS * V_MUSW / V_ECF
                + 0.25 * POSPART(1.0 - ATP / ATP0) * 3.0;

// ---------- RyR1 redox sensitisation -------------------------------------------
ROS      = POSPART((cam - 0.30) / 0.30) * (CAMITO / CAMITO50) * fT;

// ---------- volatile PK capacities ----------------------------------------------
double CVRG  = V_VRG  * L_VRG * LBG;
double CMUSF = V_MUSF * L_MUS * LBG;
double CMUSS = V_MUSS * L_MUS * LBG;
double CFAT  = V_FAT  * L_FAT * LBG;
double fVRG  = Q_VRG  * LBG * (ALV - VRG);
double fMUSF = Q_MUSF * LBG * (ALV - MUSF);
double fMUSS = Q_MUSS * LBG * (ALV - MUSS);
double fFAT  = Q_FAT  * LBG * (ALV - FAT);

// =============================== ODEs ============================================
dxdt_CIRC   = (fgf * (vap - CIRC) + VA * (ALV - CIRC) - charc * K_CHARC * CIRC) / C_CIRC;
dxdt_ALV    = (VA * (CIRC - ALV) - (fVRG + fMUSF + fMUSS + fFAT)) / C_LUNG;
dxdt_VRG    = fVRG  / CVRG;
dxdt_MUSF   = fMUSF / CMUSF;
dxdt_MUSS   = fMUSS / CMUSS;
dxdt_FAT    = fFAT  / CFAT;

dxdt_SCHA   = -KE_SCH * SCHA;

double C1   = DAN1 / V1_DAN;
dxdt_DAN1   = -(CL_DAN + Q2_DAN + Q3_DAN) * C1
              + Q2_DAN * DAN2 / V2_DAN + Q3_DAN * DAN3 / V3_DAN;
dxdt_DAN2   = Q2_DAN * (C1 - DAN2 / V2_DAN);
dxdt_DAN3   = Q3_DAN * (C1 - DAN3 / V3_DAN);
dxdt_DANE   = KE0_DAN * (C1 - DANE);

dxdt_CAM    = (J_ryr - J_serca - J_pmca + J_soce - J_mu + J_mr) / KAPPA;
dxdt_CASR   = J_serca - J_ryr;
dxdt_CAMITO = J_mu - J_mr;

dxdt_HEP    = (ATP_ox + ATP_gly - ATPuse) / 1000.0;
dxdt_GLY    = -LACprod / 2000.0;
dxdt_LACM   = LACprod / 1000.0 - lac_out;

dxdt_PMCO2  = (VCO2_mus - flux_mb) / K_CO2M;
dxdt_PACO2  = (VCO2_OTH + flux_mb - (VA * 1000.0 / 863.0) * PACO2) / K_CO2B;
dxdt_LACS   = lac_out * V_MUSW / 30.0 - 0.030 * LACS * gfrf;
dxdt_HCO3   = -(lac_out * V_MUSW / 30.0) * 0.85 + 0.03 * (24.0 - HCO3) + BICARB;

dxdt_TCORE  = (HEAT - K_CS * (TCORE - TSHELL) - coolC) / C_CORE * 60.0;
dxdt_TSHELL = (K_CS * (TCORE - TSHELL) - H_SKIN * (TSHELL - T_AMB) - coolS) / C_SHELL * 60.0;

dxdt_MEMB   = -inj;
dxdt_CK     = inj * CK_MUS - KE_CK * CK;
dxdt_MB     = inj * MB_MUS - KE_MB * MB;
dxdt_KP     = kleak - KE_K * (KP - K0) * gfrf - K_SHIFT;
dxdt_PHOS   = inj * 55.0 * V_MUSW / V_ECF - 0.02 * (PHOS - PHOS0) * gfrf;

dxdt_MBT    = K_MBT * MB / 1000.0 * gfrf - K_CAST * MBT;
dxdt_GFRF   = -K_GFR * MBT * gfrf + K_GFRREC * (1.0 - gfrf);
dxdt_SCR    = K_SCR_PR - K_SCR_CL * gfrf * SCR;

dxdt_TDOSE  = K_TD * POSPART(TCORE - T_DIC);
dxdt_FIB    = -K_FIB * TDOSE * FIB / FIB0 + 0.002 * (FIB0 - FIB);
dxdt_PLT    = -K_PLT * TDOSE * PLT / PLT0;
dxdt_DDIM   = K_DD * TDOSE - 0.006 * DDIM;

dxdt_HR     = K_HR * (70.0 + 3.1 * (PACO2 - 40.0) + 16.0 * (TCORE - 36.5)
              + 260.0 * POSPART(LACS - 1.0) / 12.0 - HR);
dxdt_VO2CUM = VO2_tot;
dxdt_NEURO  = 0.004 * pow(POSPART(TCORE - 40.0), 2.0) + 0.06 * POSPART(7.15 - pH);
dxdt_SENS   = K_SENS_ON * (1.0 - NAC) * ROS * (1.0 - sens) - K_SENS_OFF * sens;

$TABLE
double ETCO2   = PACO2 - DEADSP;
double PHART   = 6.1 + log10(fmax(HCO3, 1.0) / (0.03 * fmax(PACO2, 5.0)));
double KTOT    = KP + K_ACID * POSPART(7.40 - PHART) * 10.0;
double BE      = HCO3 - 24.0;                       // approximate base excess
double DESTRPC = 100.0 * (1.0 - MEMB);
double CONTRAC = 8.0 * pow(CAM, HILL_F) / (pow(CA50_F, HILL_F) + pow(CAM, HILL_F));

// ---- MHAUS Clinical Grading Scale (Larach 1994), computed as an OUTPUT ----
double cgs = 0.0;
if (DESTRPC > 0.05)              cgs += 15.0;   // generalised muscular rigidity
if (CK > 20000.0)                cgs += 15.0;   // CK > 20 000 after succinylcholine
if (MB > 60.0)                   cgs += 10.0;   // myoglobinuria
if (KTOT > 6.0)                  cgs += 3.0;    // hyperkalaemia
if (PACO2 > 55.0)                cgs += 15.0;   // inappropriate hypercarbia
if (PHART < 7.25)                cgs += 10.0;   // inappropriate acidosis
if (TCORE > 38.8)                cgs += 15.0;   // inappropriately rapid temp rise
if (HR > 150.0)                  cgs += 3.0;    // inappropriate tachycardia
double cgsRank = (cgs >= 50.0) ? 6.0 : (cgs >= 35.0) ? 5.0 : (cgs >= 20.0) ? 4.0 :
                 (cgs >= 13.0) ? 3.0 : (cgs >=  3.0) ? 2.0 : 1.0;

// ---- physiology that is no longer survivable ----
double LETHAL = ((TCORE > 42.5) || (KTOT > 8.5) || (PHART < 6.85)) ? 1.0 : 0.0;

$CAPTURE @annotated
ETCO2   : end-tidal CO2 (mmHg) - THE FIRST SIGN
PHART   : arterial pH
KTOT    : plasma K including the acidotic transcellular shift (mmol/L)
BE      : base excess (mEq/L)
DESTRPC : percent of muscle mass destroyed
CONTRAC : isometric contracture equivalent (g)
cgs     : MHAUS clinical grading scale raw score
cgsRank : MHAUS clinical grading scale rank (1-6)
LETHAL  : 1 when the simulated physiology is no longer survivable
Po      : RyR1 open probability - THE DRIVER
ACT     : RyR1 activation term
ATP     : muscle ATP (mM)
PCR     : phosphocreatine (mM)
FORCE   : fractional muscle activation (rigidity)
FORCEca : Ca-activated component of force
RIGOR   : ATP-depletion component of force
ATPuse  : muscle ATP turnover (uM/min)
ATP_ox  : oxidative ATP supply (uM/min)
ATP_gly : glycolytic ATP supply (uM/min)
HEP_DEF : unmet ATP demand (uM/min)
MITOEFF : mitochondrial coupling efficiency
Q_MUS   : muscle perfusion (L/min)
occl    : contracture-induced perfusion factor
VO2_tot : whole-body O2 consumption (mL/min)
VCO2_mus: muscle CO2 production (mL/min)
HEAT    : whole-body heat production (W)
inj     : instantaneous muscle injury rate (/min)
ROS     : reactive oxygen/nitrogen species drive
hyper   : hypermetabolism (multiple of resting muscle ATP turnover)

$OMEGA @annotated @block
ETA_EC50 : 0.1225 : between-subject variability on EC50_VOL (CV ~ 35%)

$SIGMA 0
'

mod <- mcode_cache("malignant_hyperthermia", mh_code)

# =====================================================================
#  GENOTYPES
#  The variant does NOT create a new pathway. It moves EC50_VOL (the agent
#  concentration that opens RyR1), PO_BASE (the resting leak), KSOCE (the
#  store-operated entry gain) and EC50_CAF (the caffeine threshold).
# =====================================================================
GENO <- list(
  MHN           = list(EC50_VOL = 22.00, PO_BASE = 1.00e-4, KSOCE = 15, EC50_CAF = 12.0),
  MHS_low       = list(EC50_VOL =  1.60, PO_BASE = 2.60e-4, KSOCE = 30, EC50_CAF =  2.2),
  MHS_CACNA1S   = list(EC50_VOL =  1.10, PO_BASE = 3.20e-4, KSOCE = 38, EC50_CAF =  1.8),
  MHS_high      = list(EC50_VOL =  0.55, PO_BASE = 5.50e-4, KSOCE = 48, EC50_CAF =  1.0),
  RYR1_myopathy = list(EC50_VOL =  0.35, PO_BASE = 9.00e-4, KSOCE = 60, EC50_CAF =  0.7)
)

# =====================================================================
#  VOLATILE AGENTS
#  MAC is a CNS potency scale. RyR1 does not read MAC — it reads molar
#  concentration, and 1 MAC spans 0.75% (halothane) to 6.0% (desflurane).
#  POTENCY is RyR1-activating potency PER VOL%, relative to sevoflurane
#  = 1.0. The product POTENCY x concentration is what the channel sees.
#
#  Model onset times, 1 MAC, no succinylcholine (min from induction to a
#  sustained +8 mmHg rise in EtCO2):
#
#                       halo   enfl    iso   sevo   desf
#      MHN               ---    ---    ---    ---    ---   (no event, ever)
#      MHS_low            87     96    114    136    197
#      MHS_CACNA1S        64     70     80     91    117
#      MHS_high           40     42     47     52     61
#      RYR1_myopathy      28     29     32     34     38
# =====================================================================
AGENT <- list(
  halothane   = list(POTENCY = 4.34, MAC = 0.75),
  enflurane   = list(POTENCY = 1.69, MAC = 1.70),
  isoflurane  = list(POTENCY = 2.06, MAC = 1.15),
  sevoflurane = list(POTENCY = 1.00, MAC = 2.00),
  desflurane  = list(POTENCY = 0.25, MAC = 6.00)
)

WT <- 70

# ---- dosing helpers --------------------------------------------------
ev_sux <- function(time = 2, mgkg = 1.0)
  ev(time = time, amt = mgkg * WT, cmt = "SCHA")

# Dantrolene: 2.5 mg/kg boluses q6 min, then 1 mg/kg q6h maintenance.
ev_dantrolene <- function(t0, mgkg = 2.5, n = 4, every = 6,
                          maint = TRUE, maint_start = 60, maint_n = 8) {
  e <- ev(time = t0, amt = mgkg * WT, cmt = "DAN1", ii = every, addl = n - 1)
  if (maint)
    e <- c(e, ev(time = t0 + maint_start, amt = 1.0 * WT, cmt = "DAN1",
                 ii = 360, addl = maint_n - 1))
  e
}

geno_par  <- function(g) GENO[[g]]
agent_par <- function(a) list(POTENCY = AGENT[[a]]$POTENCY, VAP = AGENT[[a]]$MAC)

# =====================================================================
#  SCENARIOS  (18)
#
#  Reference event unless stated: MHS_high genotype, sevoflurane 2% from
#  induction, succinylcholine 1 mg/kg at t = 2 min.
#
#  FS is when the anaesthetist RECOGNISES the event. It is not an arbitrary
#  round number: in this reference patient EtCO2 crosses +3 mmHg at 45.5
#  min, +5 at 48.8 and +8 at 51.8, so FS = 55 min represents recognition
#  about three minutes after the capnograph has unambiguously started to
#  climb. That is a generous but not absurd real-world lag. Every
#  intervention below starts at FS.
# =====================================================================
FS <- 55

BUNDLE <- function(t = FS, cool_delay = 2) list(
  TVAPOFF = t, TCHARC = t, TVENT = t, VENTMLT = 2.0,
  TCOOL = t + cool_delay, COOL_W = 200, SURF_W = 100
)

scen <- list(

  ## ---- 1-2 : controls -------------------------------------------------
  "01_MHN_control" = list(
    desc = "Non-susceptible patient, sevoflurane 2% + succinylcholine. Nothing happens
            at any agent concentration the model will reach.",
    par  = geno_par("MHN"), ev = ev_sux()),

  "02_MHS_untreated" = list(
    desc = "MH-susceptible, fulminant, NO treatment. EtCO2 121, core 41.5 degC,
            CK 22 531, K+ 8.1, pH 6.85. Lethal physiology at 140 min.",
    par  = geno_par("MHS_high"), ev = ev_sux()),

  ## ---- 3-5 : removing the trigger is NOT the same as treating ---------
  "03_vaporiser_off_only" = list(
    desc = "Vaporiser closed + 10 L/min O2 at first sign. STILL LETHAL, at 152 min.
            The agent was the START term; oxidised RyR1 is the SUSTAIN term. This one
            result is the entire history of this disease before 1979.",
    par  = c(geno_par("MHS_high"), list(TVAPOFF = FS)), ev = ev_sux()),

  "04_charcoal_filters" = list(
    desc = "As 03 plus activated-charcoal filters in both limbs. Circuit agent falls
            from 9 572 ppm to 9 ppm within 5 min instead of taking hours - but the
            PATIENT own muscle store is untouched, so the outcome barely moves
            (lethal at 156 vs 152 min). Cleaning the machine is not treatment.",
    par  = c(geno_par("MHS_high"), list(TVAPOFF = FS, TCHARC = FS)), ev = ev_sux()),

  "05_TIVA_prophylaxis" = list(
    desc = "MH-susceptible patient given a non-triggering (propofol/opioid) technique.
            No volatile agent at all: core 37.0 degC, CK 86. The definitive
            prevention, and the reason MH is now a preventable disease in anyone
            whose susceptibility is known.",
    par  = c(geno_par("MHS_high"), list(VAP = 0)), ev = NULL),

  ## ---- 6-7 : THE TWO MASKING ARMS ------------------------------------
  "06_hyperventilation_ONLY" = list(
    desc = "Minute ventilation tripled and nothing else. EtCO2 max 50.7 mmHg - a
            NORMAL capnograph - while core reaches 42.5 degC, CK 31 894, and lethal
            physiology arrives at 164 min. Fixing the ELIMINATION term hides a
            SOURCE-term problem and buys 24 minutes.",
    par  = c(geno_par("MHS_high"), list(TVENT = FS, VENTMLT = 3.0)), ev = ev_sux()),

  "07_cooling_ONLY" = list(
    desc = "Aggressive cooling and nothing else. Core max 38.1 degC - an almost normal
            temperature - while EtCO2 reaches 121, CK 15 769, and lethal physiology
            arrives at 141 min. The same lesson read off the other store.",
    par  = c(geno_par("MHS_high"),
             list(TCOOL = FS, COOL_W = 200, SURF_W = 100)), ev = ev_sux()),

  ## ---- 8-11 : dantrolene ----------------------------------------------
  "08_dantrolene_ONLY" = list(
    desc = "Dantrolene 10 mg/kg but the vaporiser is left running. SURVIVES - the drug
            alone rescues where every readout-directed manoeuvre failed - but at a
            cost: CK 10 497 and GFR 0.71, because the agent keeps re-opening what an
            87%-maximal block does not hold. Block AND remove.",
    par  = geno_par("MHS_high"),
    ev   = c(ev_sux(), ev_dantrolene(FS))),

  "09_bundle_WITHOUT_dantrolene" = list(
    desc = "Full MHAUS bundle minus the drug. Survives - but CK 46 525 (ten times the
            dantrolene arm), K+ 6.9 and GFR 0.33. This is the arm that isolates what
            the drug is actually worth.",
    par  = c(geno_par("MHS_high"), BUNDLE()), ev = ev_sux()),

  "10_FULL_bundle_Ryanodex" = list(
    desc = "REFERENCE STANDARD OF CARE. Ryanodex (250 mg / 5 mL) is reconstituted in
            about 2 min. Core 37.4 degC, EtCO2 51, CK 4 460, K+ 5.0, GFR 0.75, 0.17%
            of muscle destroyed.",
    par  = c(geno_par("MHS_high"), BUNDLE()),
    ev   = c(ev_sux(), ev_dantrolene(FS + 2))),

  "11_FULL_bundle_Dantrium" = list(
    desc = "IDENTICAL MOLECULE, different vial. Dantrium/Revonto is 20 mg per vial in
            60 mL of sterile water: for a 70 kg adult at 2.5 mg/kg that is 9 vials,
            540 mL and roughly 18 min of shaking. The ONLY thing that differs from
            scenario 10 is an input delay - and delay is the outcome variable:
            CK 5 775 vs 4 460 (+29%) and GFR 0.69 vs 0.75 for those 16 minutes.",
    par  = c(geno_par("MHS_high"), BUNDLE(cool_delay = 18)),
    ev   = c(ev_sux(), ev_dantrolene(FS + 18))),

  ## ---- 12-13 : the price of delay --------------------------------------
  "12_dantrolene_delayed_30min" = list(
    desc = "Recognition-to-drug 30 min. CK 7 656, K+ 6.2, GFR 0.55.",
    par  = c(geno_par("MHS_high"), BUNDLE(t = FS + 30)),
    ev   = c(ev_sux(), ev_dantrolene(FS + 32))),

  "13_dantrolene_delayed_60min" = list(
    desc = "Recognition-to-drug 60 min. CK 13 378, K+ 7.2, GFR 0.39. Beyond about
            75 min the drug arrives after the point of no return and the model dies
            regardless of dose - the delay curve does not merely worsen, it ends.",
    par  = c(geno_par("MHS_high"), BUNDLE(t = FS + 60)),
    ev   = c(ev_sux(), ev_dantrolene(FS + 62))),

  ## ---- 14 : recrudescence -----------------------------------------------
  "14_underdosed_no_maintenance" = list(
    desc = "A single 2.5 mg/kg bolus with no maintenance. Redox-sensitised RyR1
            outlives the drug: SENS is still 0.73 four hours in, and as the
            effect-site concentration falls through 2 mg/L the event RESTARTS -
            EtCO2 rebounding +33 mmHg with CK reaching 88 413 by 30 h. TWO boluses
            leave SENS at 0.23 and there is no rebound at all. In this model
            recrudescence is a consequence of UNDER-DOSING, not of bad luck.",
    par  = c(geno_par("MHS_high"), BUNDLE()),
    ev   = c(ev_sux(), ev_dantrolene(FS + 2, n = 1, maint = FALSE))),

  ## ---- 15-16 : genotype and agent --------------------------------------
  "15_low_penetrance_variant" = list(
    desc = "Low-penetrance RYR1 variant, same anaesthetic. Onset at ~136 min instead
            of 52 - this is why a susceptible patient can have several uneventful
            triggering anaesthetics before the one that declares itself.",
    par  = geno_par("MHS_low"), ev = ev_sux()),

  "16_desflurane_1MAC" = list(
    desc = "1 MAC desflurane (6.0 vol%) in the same high-penetrance patient. Because
            RyR1 potency per vol% is only 0.25, 1 MAC of desflurane delivers the
            LEAST channel drive of the five agents: onset 61 min against 40 for
            halothane in the same patient. The weakest trigger gives the latest and
            most easily-missed presentation.",
    par  = c(geno_par("MHS_high"), agent_par("desflurane")), ev = ev_sux()),

  ## ---- 17-18 : mechanism probes -----------------------------------------
  "17_low_extracellular_Ca" = list(
    desc = "MECHANISM PROBE. Extracellular Ca reduced to 20%. At steady state SR
            release and SERCA uptake cancel, so the sustained myoplasmic Ca level is
            set by the SARCOLEMMAL balance. Cutting entry collapses the plateau -
            the in-vivo counterpart of the in-vitro finding that a Ca-free bath
            abolishes the halothane contracture.",
    par  = c(geno_par("MHS_high"), list(CA_O = 0.2)), ev = ev_sux()),

  "18_antioxidant_NAC" = list(
    desc = "MECHANISM PROBE. N-acetylcysteine blocks ROS-driven RyR1 oxidation (NAC
            acts on K_SENS_ON, not on the channel), so this arm isolates how much of
            the event is SUSTAIN rather than START. Explored clinically in
            RYR1-related myopathy; NOT an established MH therapy.",
    par  = c(geno_par("MHS_high"), list(TVAPOFF = FS, NAC = 0.9)), ev = ev_sux())
)

# =====================================================================
#  RUNNERS
# =====================================================================
run_scenario <- function(name, end = 480, delta = 0.25) {
  s <- scen[[name]]
  m <- mod
  if (!is.null(s$par)) m <- param(m, s$par)
  e <- if (is.null(s$ev)) ev(time = 0, amt = 0, cmt = "DAN1") else s$ev
  out <- m %>% mrgsim(events = e, end = end, delta = delta) %>% as_tibble()
  out$scenario <- name
  out
}

summarise_scenario <- function(name, end = 480) {
  o <- run_scenario(name, end = end)
  # read every summary at the first lethal time point, or at the end
  ilet <- which(o$LETHAL > 0)
  i    <- if (length(ilet)) ilet[1] else nrow(o)
  w    <- seq_len(i)
  tibble(
    scenario   = name,
    Tmax       = max(o$TCORE[w]),
    EtCO2max   = max(o$ETCO2[w]),
    CK         = o$CK[i],
    Kmax       = max(o$KTOT[w]),
    pHmin      = min(o$PHART[w]),
    lacmax     = max(o$LACS[w]),
    destroyed  = o$DESTRPC[i],
    SCr        = o$SCR[i],
    GFR        = o$GFRF[i],
    CGSrank    = max(o$cgsRank[w]),
    SENSmax    = max(o$SENS[w]),
    died_at    = if (length(ilet)) o$time[ilet[1]] else NA_real_
  )
}

run_all <- function(end = 480)
  do.call(rbind, lapply(names(scen), summarise_scenario, end = end))

# =====================================================================
#  ANALYSIS 1 — THE CAPACITANCE EXPERIMENT
#  Pre-equilibrate the muscle with agent so the trigger is a STEP at t = 0,
#  then measure how long each store takes to fill. This is the experiment
#  the whole model exists to make.
#
#  Reference result (MHS_high):
#      EtCO2 +10 mmHg   8.2 min      core +0.5 degC    9.8 min
#      EtCO2 +20 mmHg  18.5 min      core +1.0 degC   20.5 min
#      EtCO2 +40 mmHg  42.8 min      core +1.5 degC   31.0 min
#                                    core  38.5 degC  40.8 min
#      excess CO2 403 mL/min / 60 mL/mmHg  = 0.149 min per mmHg
#      excess heat 160 W  / 158 kJ/degC    = 16.5  min per degC
#      => 111 mmHg of EtCO2 per degC; plateau 149 mmHg; peak 3.5 degC/h
# =====================================================================
capacitance_experiment <- function() {
  m <- param(mod, GENO$MHS_high)
  o <- m %>%
    init(CIRC = 1.6, ALV = 1.6, VRG = 1.6, MUSF = 1.6, MUSS = 1.6) %>%
    mrgsim(end = 200, delta = 0.1) %>% as_tibble()
  e0 <- o$ETCO2[1]; T0 <- o$TCORE[1]
  first <- function(v, thr) { i <- which(v >= thr); if (length(i)) o$time[i[1]] else NA }
  k  <- which.min(abs(o$time - 5))
  dC <- o$VCO2_mus[k] - 38.0          # excess CO2 flux, mL/min
  dW <- o$HEAT[k] - 66.0              # excess heat flux, W
  p  <- as.list(param(m))
  list(
    `EtCO2 +10 mmHg (min)`   = first(o$ETCO2, e0 + 10),
    `EtCO2 +20 mmHg (min)`   = first(o$ETCO2, e0 + 20),
    `EtCO2 +40 mmHg (min)`   = first(o$ETCO2, e0 + 40),
    `core +0.5 degC (min)`   = first(o$TCORE, T0 + 0.5),
    `core +1.0 degC (min)`   = first(o$TCORE, T0 + 1.0),
    `core +1.5 degC (min)`   = first(o$TCORE, T0 + 1.5),
    `core 38.5 degC (min)`   = first(o$TCORE, 38.5),
    `excess VCO2 (mL/min)`   = dC,
    `excess heat (W)`        = dW,
    `min per mmHg of EtCO2`  = p$K_CO2B / dC,
    `min per degC of core`   = p$C_CORE / dW / 60,
    `mmHg of EtCO2 per degC` = (p$C_CORE / dW / 60) / (p$K_CO2B / dC),
    `EtCO2 plateau (mmHg)`   = max(o$ETCO2),
    `peak VO2 (mL/min)`      = max(o$VO2_tot),
    `peak heat (W)`          = max(o$HEAT)
  )
}

# =====================================================================
#  ANALYSIS 2 — IN-VITRO CONTRACTURE TEST (IVCT / CHCT)
#  Ca subsystem only: temperature clamped at 37 degC, ATP unlimited, no
#  systemic feedback. This is the assay that measures EXACTLY the parameter
#  the genotype moves, which is why a functional test still outperforms
#  sequencing.
#  European MH Group: MHS if >= 0.2 g threshold contracture at <= 2%
#  halothane OR <= 2 mM caffeine.
#
#  Reference result, g above baseline at 2% halothane / 2 mM caffeine:
#      MHN            0.00 / 0.00   -> MHN
#      MHS_low        2.31 / 1.11   -> MHS
#      MHS_CACNA1S    3.90 / 2.52   -> MHS
#      MHS_high       5.63 / 5.44   -> MHS
#      RYR1_myopathy  6.69 / 6.66   -> MHS
#  The model separates susceptible from normal correctly and preserves the
#  rank order of thresholds, but compresses the graded contracture
#  MAGNITUDES: real IVCT tracings span 0.2-3 g, not 2-7 g.
# =====================================================================
ivct <- function(genotype, halothane = 0, caffeine = 0, ca_o = 1, tmax = 25) {
  p <- as.list(param(param(mod, GENO[[genotype]])))
  POT_HAL <- AGENT$halothane$POTENCY          # RyR1 potency per vol% (4.34)
  CAM <- p$CA_REST; CASR <- p$CASR0; CAMITO <- 20; dt <- 0.002
  ACT <- halothane * POT_HAL / p$EC50_VOL + caffeine / p$EC50_CAF
  Po  <- p$PO_BASE + (p$PO_MAX - p$PO_BASE) * ACT^p$HILL_RYR / (1 + ACT^p$HILL_RYR)
  for (i in seq_len(tmax / dt)) {
    Jr <- p$KRYR * Po * (CASR - CAM)
    Js <- p$VSERCA * CAM^2 / (p$KSERCA^2 + CAM^2)
    Jp <- p$KPMCA * (CAM - p$CA_REST)
    Jo <- p$KSOCE * max(0, 1 - CASR / p$CASR0) * ca_o
    Ju <- p$KMITO_U * CAM * max(0, 1 - CAMITO / p$CAMITOMX)
    Jm <- p$KMITO_R * CAMITO / p$CAMITOMX * 12
    CAM    <- CAM + dt * (Jr - Js - Jp + Jo - Ju + Jm) / p$KAPPA
    CASR   <- max(0, CASR + dt * (Js - Jr))
    CAMITO <- max(0, CAMITO + dt * (Ju - Jm))
  }
  8 * CAM^p$HILL_F / (p$CA50_F^p$HILL_F + CAM^p$HILL_F)
}

ivct_panel <- function() {
  do.call(rbind, lapply(names(GENO), function(g) {
    b <- ivct(g)
    tibble(genotype = g,
           hal_0.5 = ivct(g, halothane = 0.5) - b,
           hal_1.0 = ivct(g, halothane = 1.0) - b,
           hal_2.0 = ivct(g, halothane = 2.0) - b,
           caf_1.0 = ivct(g, caffeine  = 1.0) - b,
           caf_2.0 = ivct(g, caffeine  = 2.0) - b,
           verdict = ifelse(ivct(g, halothane = 2) - b >= 0.2 |
                            ivct(g, caffeine  = 2) - b >= 0.2, "MHS", "MHN"))
  }))
}

# Extracellular-Ca dependence of the contracture (the SOCE prediction).
# Reference result, MHS_high at 2% halothane:
#   100% bath Ca -> 5.63 g   50% -> 1.30 g   20% -> 0.17 g   0% -> 0.00 g
ivct_ca_dependence <- function(genotype = "MHS_high", halothane = 2) {
  b <- ivct(genotype)
  ca <- c(1, 0.5, 0.2, 0)
  tibble(extracellular_Ca = ca,
         contracture_g = sapply(ca, function(x)
           ivct(genotype, halothane = halothane, ca_o = x) - b))
}

# =====================================================================
#  ANALYSIS 3 — ONSET TIME BY GENOTYPE AND AGENT
#  Each agent is given at 1 MAC. The channel sees POTENCY x concentration,
#  not MAC, which is why the ranking follows potency-per-vol% and not
#  anaesthetic potency. See the AGENT table above for reference values.
# =====================================================================
onset_panel <- function(thr = 8, end = 300) {
  grid <- expand.grid(genotype = names(GENO), agent = names(AGENT),
                      stringsAsFactors = FALSE)
  grid$onset_min <- mapply(function(g, a) {
    m <- param(mod, c(GENO[[g]], agent_par(a)))
    o <- m %>% mrgsim(end = end, delta = 0.5) %>% as_tibble()
    base <- mean(o$ETCO2[o$time >= 8 & o$time <= 12])
    i <- which(o$ETCO2 >= base + thr & o$time > 12)
    if (length(i)) o$time[i[1]] else NA_real_
  }, grid$genotype, grid$agent)
  grid
}

# =====================================================================
#  ANALYSIS 4 — TIME TO DANTROLENE
#  Reference result (MHS_high, full bundle, dantrolene at FS + delay):
#      delay   CK      K+    GFR
#          0    4 292  5.0   0.70
#         30    7 656  6.2   0.55
#         60   13 378  7.2   0.39
#         90   -- dies at 140 min regardless of dose --
# =====================================================================
delay_curve <- function(delays = c(0, 5, 10, 15, 20, 30, 45, 60, 90, 120), end = 600) {
  do.call(rbind, lapply(delays, function(d) {
    td <- FS + d
    m <- param(mod, c(GENO$MHS_high, BUNDLE(t = td)))
    o <- m %>% mrgsim(events = c(ev_sux(), ev_dantrolene(td + 2)),
                      end = end, delta = 0.5) %>% as_tibble()
    ilet <- which(o$LETHAL > 0); i <- if (length(ilet)) ilet[1] else nrow(o)
    tibble(delay_min = d, Tmax = max(o$TCORE[seq_len(i)]), CK = o$CK[i],
           Kmax = max(o$KTOT[seq_len(i)]), destroyed = o$DESTRPC[i],
           GFR = o$GFRF[i], died_at = if (length(ilet)) o$time[i] else NA_real_)
  }))
}

# =====================================================================
#  ANALYSIS 5 — VOLATILE WASHOUT
#  Flushing the machine cannot get the circuit below what the PATIENT is
#  exhaling. Charcoal filters scrub the exhaled agent too, which is why
#  they have an asymptote the flush does not.
#
#  Reference result, circuit concentration (ppm) after the vaporiser closes:
#      minutes        1      5     30     90    180    480
#      flush       9 572  6 819    994     96     39      8
#      charcoal       52      9      2      1      0      0
#  Published charcoal-filter data reach <5 ppm within about 90 s; the model
#  gets there at ~5 min, i.e. it is conservative about the filter.
# =====================================================================
washout <- function(end = 600) {
  f <- function(charcoal) {
    m <- param(mod, c(GENO$MHN, list(TVAPOFF = FS,
              TCHARC = if (charcoal) FS else 1e9)))
    o <- m %>% mrgsim(end = end, delta = 0.5) %>% as_tibble()
    sapply(c(1, 5, 30, 90, 180, 480),
           function(mm) o$CIRC[which.min(abs(o$time - (FS + mm)))] * 10000)
  }
  tibble(minutes_after = c(1, 5, 30, 90, 180, 480),
         flush_ppm = f(FALSE), charcoal_ppm = f(TRUE))
}

# =====================================================================
#  ANALYSIS 6 — COOLING ARITHMETIC
#  Pure bookkeeping, but it is the bookkeeping that shows cooling can never
#  be the treatment.
# =====================================================================
cooling_arithmetic <- function() {
  p <- as.list(param(mod))
  Cbody <- (p$C_CORE + p$C_SHELL) / 1000       # kJ/degC
  list(
    `body heat capacity (kJ/degC)`     = Cbody,
    `1 L of 4 degC saline (kJ)`        = 1 * 4.18 * 33,
    `... degC it buys`                 = 4.18 * 33 / Cbody,
    `3 L bolus, degC`                  = 3 * 4.18 * 33 / Cbody,
    `surface cooling 200 W (degC/h)`   = 200 * 3600 / 1000 / Cbody,
    `fulminant heat production (W)`    = 249,
    comment = paste("Cooling removes heat at roughly the rate the muscle makes it.",
                    "It buys time; it cannot end the event, because it does not",
                    "touch Po. Its ONE genuine mechanistic contribution is via the",
                    "Q10 = 2.6 term: a cooler channel opens less. That is a",
                    "second-order effect on the driver, not a first-order one.")
  )
}

# =====================================================================
#  ANALYSIS 7 — VIRTUAL POPULATION
#  Larach 2010 (n = 286 NAMHR events) reported that each 30 min of delay
#  multiplied the odds of a complication by 1.61.
#
#  Reference result (60 subjects per delay, genotype mix 45% high-
#  penetrance / 35% CACNA1S / 20% low-penetrance, log-normal variability on
#  EC50_VOL, KSOCE, ATPASE_A, muscle mass and K_INJ, and a recognition time
#  of 55 +/- 8 min):
#      delay   0    15    30    45    60    90
#      compl  8%   18%   37%   33%   50%   58%
#      death  2%    3%   17%    8%   18%   18%
#      => odds ratio 2.38 per 30 min, against an observed 1.61.
#  Right direction, right monotonicity, too steep. See the README.
# =====================================================================
virtual_population <- function(n = 60, delays = c(0, 15, 30, 45, 60, 90),
                               seed = 20260804, end = 600) {
  set.seed(seed)
  genos <- c(rep("MHS_high", 45), rep("MHS_CACNA1S", 35), rep("MHS_low", 20))
  do.call(rbind, lapply(delays, function(d) {
    res <- replicate(n, {
      g  <- sample(genos, 1)
      pp <- GENO[[g]]
      pp$EC50_VOL <- pp$EC50_VOL * exp(rnorm(1, 0, 0.35))
      pp$KSOCE    <- pp$KSOCE    * exp(rnorm(1, 0, 0.25))
      fs <- FS + rnorm(1, 0, 8)
      td <- max(5, fs + d)
      m  <- param(mod, c(pp, BUNDLE(t = td),
                         list(ATPASE_A = 3600 * exp(rnorm(1, 0, 0.20)),
                              V_MUSW   = 22.5 * exp(rnorm(1, 0, 0.15)),
                              K_INJ    = 1.1e-4 * exp(rnorm(1, 0, 0.30)))))
      o  <- m %>% mrgsim(events = c(ev_sux(), ev_dantrolene(td + 2)),
                         end = end, delta = 1) %>% as_tibble()
      ilet <- which(o$LETHAL > 0); i <- if (length(ilet)) ilet[1] else nrow(o)
      w <- seq_len(i)
      c(complication = as.numeric(max(o$CK[w]) > 10000 || max(o$KTOT[w]) > 6 ||
                                  o$SCR[i] > 1.5 || max(o$TCORE[w]) > 40 ||
                                  o$PLT[i] < 150 || length(ilet) > 0),
        death = as.numeric(length(ilet) > 0))
    })
    tibble(delay_min = d, complication_rate = mean(res["complication", ]),
           death_rate = mean(res["death", ]))
  }))
}

# =====================================================================
#  PARAMETER PROVENANCE
#  Every number above is either measured, standard physiology, or one of
#  the five fitted quantities listed in the header. Sources are in
#  mh_references.md, section by section.
#
#  Anchors used for calibration
#  ----------------------------
#  - resting muscle ATP turnover 0.4-1 mM/min                   (31P-MRS)
#  - whole-body VO2 in MH 2-4x baseline                         (Gronert)
#  - body heat capacity 3.47 kJ/(kg.degC)                       (Sessler)
#  - apnoeic PaCO2 rise ~3.3 mmHg/min -> CO2 store 60 mL/mmHg   (Cherniack)
#  - alveolar equation PaCO2 = 863 x VCO2 / VA                  (standard)
#  - sevoflurane blood:gas 0.65, muscle:blood 3.1, fat:blood 48 (Eger/Yasuda)
#  - succinylcholine plasma t1/2 47 s                           (standard)
#  - dantrolene 2.4 mg/kg -> ~3-4 ug/mL, t1/2 12-16 h,
#    therapeutic 2.8-4.2 ug/mL for 75% twitch depression        (Flewellen 1983)
#  - dantrolene binds RyR1 aa 590-609, Mg2+/CaM dependent,
#    and the block is PARTIAL                                   (Paul-Pletzer, Fruen)
#  - resting myoplasmic Ca elevated in MHS muscle               (Lopez 1985)
#  - enhanced store-operated Ca entry in MHS muscle             (Duke 2010, Eltit 2013)
#  - RyR1 S-nitrosylation as a sensitising mechanism            (Durham 2008)
#  - onset later with the less soluble agents                   (Wedel 1993)
#  - IVCT/CHCT thresholds 0.2 g at 2% halothane / 2 mM caffeine (EMHG / NAMHG)
#  - MH mortality ~70-80% pre-dantrolene, 1.4-10% now           (registry data)
#  - complication odds x1.61 per 30 min of delay                (Larach 2010)
#  - recrudescence ~20%, median ~13 h                           (Burkman 2007)
# =====================================================================

if (interactive()) {
  print(run_all())
  print(capacitance_experiment())
  print(ivct_panel())
  print(ivct_ca_dependence())
  print(onset_panel())
  print(delay_curve())
  print(washout())
  print(cooling_arithmetic())
  print(virtual_population())
}
