$PROB
# Generalized Anxiety Disorder (GAD) -- Quantitative Systems Pharmacology model
#
# ---------------------------------------------------------------------------
# THE CENTRAL OBJECT
# ---------------------------------------------------------------------------
# Anxiety is written as ONE dimensionless corticolimbic gain
#
#        Phi = (E_amy * S_glu) / (C_pfc * (1 + e_ex*EXPECT) * I_gaba)
#
# a RATIO with two numerator factors and two denominator factors.  Each of the
# four factors is owned by a different drug class and each moves on its own
# time constant:
#
#   I_gaba   GABA-A inhibitory efficacy     benzodiazepines      hours
#   S_glu    glutamate release probability  pregabalin (a2d-1)   ~2 days
#   C_pfc    prefrontal control capacity    SSRI/SNRI, CBT       ~13 d gate
#                                                                + ~30 d growth
#   E_amy    amygdala / BNST drive          CBT, chronic 5-HT    weeks
#
# and the whole thing is read out through HAM-A, an instrument that carries a
# FIFTH clock of its own (expectancy + regression to the mean) which is
# subtracted from every arm equally.  Because HAM-A is a SATURATING function of
# Phi, several things that are usually asserted become arithmetic:
#
#   * "onset of action" is a property of WHICH FACTOR a drug moves;
#   * the flat SSRI dose-response is the SERT occupancy hyperbola, not a fit;
#   * drugs are MULTIPLICATIVE in Phi and SUB-ADDITIVE on HAM-A;
#   * a high-expectancy site has a smaller MEASURABLE delta with identical
#     pharmacology (assay sensitivity, derived);
#   * one benzodiazepine adaptation state (DEPEND) gives partial tolerance on
#     drug and rebound off it -- the same parameter with opposite sign.
#
# ---------------------------------------------------------------------------
# WHAT WAS FITTED
# ---------------------------------------------------------------------------
# FIVE parameters, to five numbers:
#   FLUCT0, dvisit, e_ex  <- Khan 2011 (PMID 21694613) placebo arm
#                            (week 1 -5.94, week 8 -11.10) and Rickels 2005
#                            (PMID 16143734) placebo arm week 4 (-8.4)
#   emax_pgb              <- Rickels 2005 pregabalin 300 mg week 4 (-12.2)
#   emax_a2               <- Rickels 2005 alprazolam 1.5 mg week 4 (-10.9)
# Everything else is a measured PK or occupancy constant, or structure.
# All other trial numbers reproduced by this model are out-of-sample.
#
# ---------------------------------------------------------------------------
# RELATIONSHIP TO THE PYTHON REFERENCE
# ---------------------------------------------------------------------------
# gad_python_reference.py is the executable reference; every number quoted in
# README.md comes from it.  This file mirrors it equation for equation with ONE
# deliberate difference: the Python code applies the study-visit expectancy
# update as an exact discrete map EXPECT <- EXPECT + dvisit*(Emax - EXPECT),
# which mrgsolve cannot express as a dose.  Here the same map is produced by a
# fast decaying cue compartment VISCUE (bolus of 1 unit at each visit, kvis =
# 48/d) driving
#       dEXPECT/dt = kexp_on*VISCUE*(Emax - EXPECT) - kexp_dec*EXPECT
# with kexp_on = -log(1 - dvisit)*kvis, which reproduces the discrete map in
# the limit of a fast cue.  NOTE ALSO, honestly: no R toolchain was available
# in the environment where this model was built, so this file has been checked
# by inspection against the Python reference but has NOT itself been executed.
#
# Compartments: 50 (49 ODE states of the reference + VISCUE).
#
# Dosing (all amounts in ug; 1 mg = 1000 ug):
#   cmt 1  ESC_GUT   escitalopram        QD
#   cmt 5  VEN_GUT   venlafaxine ER      QD
#   cmt 8  DLX_GUT   duloxetine          QD
#   cmt 10 PGB_GUT   pregabalin          BID
#   cmt 12 BZD_GUT   lorazepam/alprazolam TID
#   cmt 15 BUS_GUT   buspirone           TID
#   cmt 18 QTP_GUT   quetiapine XR       QD
#   cmt 50 VISCUE    study visit         amt = 1 at each visit time
#
# Alprazolam: set BZD_POT = 1.53, F_bzd = 0.88, V_bzd = 65, CL_bzd = 90.
# CBT: set CBT = 1 during the therapy window (via a time-varying data item).

$PLUGIN autodec

$PARAM @annotated
// ---- escitalopram --------------------------------------------------------
F_esc      : 0.80   : escitalopram bioavailability
ka_esc     : 9.6    : escitalopram absorption rate (1/d)
V1_esc     : 350    : central volume (L)
V2_esc     : 400    : peripheral volume (L)
Q_esc      : 250    : intercompartmental clearance (L/d)
CL_esc     : 400    : clearance (L/d) -- Css 20 ng/mL at 10 mg QD
keo_esc    : 8.0    : brain effect-site equilibration (1/d)
EC50_sert_esc : 5.0 : SERT EC50 (ng/mL) -- 80 pct occupancy at 10 mg
// ---- venlafaxine ER + O-desmethylvenlafaxine ------------------------------
F_ven      : 0.42   : venlafaxine bioavailability
ka_ven     : 2.4    : ER absorption rate (1/d)
V_ven      : 300    : venlafaxine volume (L)
CL_ven     : 999    : venlafaxine total clearance (L/d)
V_odv      : 250    : ODV volume (L)
CL_odv     : 378    : ODV clearance (L/d)
mw_ratio_odv : 0.95 : molar correction parent to metabolite
fm2d6      : 0.70   : CYP2D6 metabolic fraction (EM 0.70 / PM 0.15 / UM 0.85)
EC50_sert_ven : 22.5 : SERT EC50 for total active moiety (ng/mL)
EC50_net_ven  : 260  : NET EC50 for total active moiety (ng/mL)
// ---- duloxetine -----------------------------------------------------------
F_dlx      : 0.50   : duloxetine bioavailability
ka_dlx     : 4.0    : absorption rate (1/d)
V_dlx      : 700    : volume (L)
CL_dlx     : 970    : clearance (L/d)
EC50_sert_dlx : 7.7 : SERT EC50 (ng/mL)
EC50_net_dlx  : 46   : NET EC50 (ng/mL)
// ---- pregabalin -----------------------------------------------------------
F_pgb      : 0.90   : pregabalin bioavailability
ka_pgb     : 16.0   : absorption rate (1/d)
V_pgb      : 40     : volume (L)
CL_pgb     : 105    : clearance at CrCl 100 (L/d) -- t1/2 6.3 h
crcl       : 100    : creatinine clearance (mL/min)
EC50_a2d   : 2100   : alpha2delta EC50 (ng/mL)
ktraf      : 0.50   : alpha2delta trafficking rate (1/d) -- tau 2 d
emax_pgb   : @@emax_pgb@@ : maximum fractional suppression of S_glu   [FITTED]
// ---- benzodiazepine --------------------------------------------------------
F_bzd      : 0.90   : bioavailability (lorazepam)
ka_bzd     : 12.0   : absorption rate (1/d)
V_bzd      : 90     : volume (L)
CL_bzd     : 108    : clearance (L/d)
keo_bzd    : 14.0   : brain effect-site equilibration (1/d)
BZD_POT    : 1.0    : lorazepam-equivalent potency (alprazolam 1.53)
EC50_bz    : 96     : BZ-site EC50 (ng/mL lorazepam-equivalent, Atack 2007)
emax_a2    : @@emax_a2@@ : maximum alpha2/3 potentiation of I_gaba    [FITTED]
wadapt     : 0.20   : fraction of the effect the adaptation state subtracts
kdep_on    : 0.05   : adaptation rate (1/d) -- tau 20 d
ktol1      : 1.20   : alpha1 pool tolerance rate (1/d per unit occupancy)
krec1      : 0.15   : alpha1 pool recovery rate (1/d)
ktol2      : 0.012  : alpha2/3 pool tolerance rate (1/d per unit occupancy)
krec2      : 0.020  : alpha2/3 pool recovery rate (1/d)
// ---- buspirone + 1-PP -------------------------------------------------------
F_bus      : 0.05   : buspirone bioavailability (high first pass)
ka_bus     : 14.0   : absorption rate (1/d)
V_bus      : 350    : volume (L)
CL_bus     : 1800   : clearance (L/d)
f_pp       : 0.30   : fraction of dose appearing as 1-PP
V_pp       : 200    : 1-PP volume (L)
CL_pp      : 554    : 1-PP clearance (L/d)
EC50_1a    : 1.20   : 5-HT1A EC50 for buspirone (ng/mL)
ia_bus     : 0.35   : intrinsic activity at the autoreceptor
ia_bus_post : 0.55  : postsynaptic contribution coefficient
EC50_pp1   : 45     : alpha2-adrenergic EC50 for 1-PP (ng/mL)
kpp1       : 0.25   : 1-PP disinhibition of LC firing
// ---- quetiapine XR + norquetiapine -------------------------------------------
F_qtp      : 1.0    : apparent bioavailability (CL and V are apparent)
ka_qtp     : 3.5    : XR absorption rate (1/d)
V_qtp      : 700    : apparent volume (L)
CL_qtp     : 2000   : apparent clearance (L/d) -- Css 75 ng/mL at 150 mg
f_nqt      : 0.30   : fraction of clearance to norquetiapine
V_nqt      : 600    : norquetiapine volume (L)
CL_nqt     : 1500   : norquetiapine clearance (L/d)
EC50_h1    : 20     : H1 EC50 (ng/mL)
EC50_net_nqt : 60   : NET EC50 for norquetiapine (ng/mL)
// ---- serotonin ------------------------------------------------------------------
kin_sht    : 24.0   : 5-HT appearance rate (1/d)
kout_sht   : 24.0   : 5-HT removal rate (1/d)
sert_floor : 0.15   : non-SERT fraction of 5-HT clearance
gamma_auto : 6.0    : autoreceptor feedback gain
kdes_auto  : 0.24   : autoreceptor desensitisation rate
krec_auto  : 0.020  : autoreceptor recovery rate -- tau 13 d net
// ---- noradrenaline ----------------------------------------------------------------
kin_ne     : 18.0   : NE appearance rate (1/d)
kout_ne    : 18.0   : NE removal rate (1/d)
net_floor  : 0.25   : non-NET fraction of NE clearance
gamma_a2   : 5.0    : alpha2 autoreceptor feedback gain
kdes_a2    : 0.18   : alpha2 autoreceptor desensitisation rate
krec_a2    : 0.025  : alpha2 autoreceptor recovery rate
// ---- plasticity / PFC / amygdala ------------------------------------------------------
kb         : 0.20   : BDNF formation rate (1/d)
kdb        : 0.20   : BDNF loss rate (1/d)
w_5ht      : 0.55   : weight of 5-HT on BDNF
w_ne       : 0.25   : weight of NE on BDNF
w_tonic    : 0.20   : tonic BDNF drive
kcort_bdnf : 0.35   : cortisol suppression of BDNF
cbt_bdnf   : 0.30   : CBT contribution to BDNF
kp_pfc     : 0.030  : C_pfc growth rate (1/d)
kl_pfc     : 0.0175 : C_pfc loss rate (1/d)
cpfc_max   : 2.20   : C_pfc ceiling
kc_pfc     : 0.45   : cortisol acceleration of C_pfc loss
kdis_pfc   : 1.14   : disease acceleration of C_pfc loss
e_ex       : @@e_ex@@ : expectancy coefficient on C_pfc              [FITTED]
ke_in      : 0.050  : E_amy formation rate (1/d)
ke_out     : 0.050  : E_amy loss rate (1/d)
EAMAX      : 3.0    : E_amy logistic ceiling
kcrh_amy   : 0.15   : CRH drive on E_amy
kslp_amy   : 0.18   : sleep deficit drive on E_amy
kdis_amy   : 0.80   : disease drive on E_amy
a5ht_amy   : 0.45   : 5-HT acceleration of E_amy decay
acbt_amy   : 0.55   : CBT acceleration of E_amy decay
kglu_stress : 0.35  : E_amy drive on glutamate release probability
kglu_cort  : 0.20   : cortisol drive on glutamate release probability
igaba_base0 : 1.0   : healthy GABAergic tone
kdis_gaba  : 0.12   : disease reduction of baseline GABAergic tone
// ---- bounded excess gain -----------------------------------------------------------------
ZMAX       : 3.0    : ceiling on the saturating excess z
kcrh_phi   : 0.12   : weight of z on CRH drive
ks_phi     : 1.00   : weight of z on sympathetic tone
ksl_phi    : 1.00   : weight of z on sleep deficit
kw_phi     : 1.00   : weight of z on the worry engine
SLMAX      : 4.0    : sleep deficit ceiling
// ---- HPA axis --------------------------------------------------------------------------
kcrh       : 1.395  : CRH drive (normalises healthy cortisol to 1)
kcout      : 0.9    : CRH elimination (1/d)
kgr_fb     : 0.55   : strength of GR-mediated negative feedback
ka_acth    : 6.0    : ACTH formation (1/d)
ke_acth    : 6.0    : ACTH elimination (1/d)
kc_cort    : 4.0    : cortisol formation (1/d)
kec_cort   : 4.0    : cortisol elimination (1/d)
kgr_in     : 0.08   : GR recovery (1/d)
kgr_dn     : 0.10   : GR downregulation by excess cortisol
// ---- autonomic and symptom layer ----------------------------------------------------------
ks_in      : 0.60   : sympathetic formation (1/d)
ks_out     : 0.60   : sympathetic loss (1/d)
w_ne_sns   : 0.55   : NE drive on sympathetic tone
kau_in     : 0.45   : autonomic symptom formation (1/d)
kau_out    : 0.45   : autonomic symptom loss (1/d)
b_gaba     : 0.85   : GABA relief of autonomic symptoms
b_pgb      : 0.75   : alpha2delta relief of autonomic symptoms
ksl_in     : 0.35   : sleep deficit formation (1/d)
ksl_out    : 0.35   : sleep deficit resolution (1/d)
w_cort_sl  : 0.30   : cortisol drive on sleep deficit
h_h1       : 1.30   : H1 blockade relief of sleep deficit
h_gaba     : 1.60   : benzodiazepine relief of sleep deficit
h_pgb      : 0.90   : alpha2delta relief of sleep deficit
kw_in      : 0.16   : worry formation (1/d)
kw_out     : 0.16   : worry resolution (1/d)
rho_worry  : 0.55   : self-recruitment gain of the worry engine
kw_m       : 1.20   : half-saturation of self-recruitment
acbt_worry : 0.80   : CBT acceleration of worry resolution
a5ht_worry : 0.35   : 5-HT acceleration of worry resolution
// ---- trial machinery ------------------------------------------------------------------------
kexp_dec   : 0.010  : expectancy decay (1/d)
dvisit     : @@dvisit@@ : per-visit expectancy increment              [FITTED]
emax_exp   : 1.00   : expectancy ceiling
kvis       : 48.0   : visit-cue decay rate (1/d); see header note
kfl        : 0.075  : enrolment-peak decay (1/d) -- tau 13 d
FLUCT0     : @@fluct0@@ : enrolment peak, HAM-A points               [FITTED]
// ---- HAM-A link ---------------------------------------------------------------------------------
c1         : 0.62   : worry weight, psychic cluster
c2         : 0.30   : Phi_n weight, psychic cluster
c3         : 0.26   : sleep deficit weight, psychic cluster
Kpsy       : 3.02   : psychic half-saturation
psy0       : 1.098  : healthy psychic raw score
psy_floor  : 1.70   : healthy psychic HAM-A
v1         : 0.42   : autonomic weight, somatic cluster
v2         : 0.30   : sympathetic weight, somatic cluster
v3         : 0.30   : muscle tension (S_glu) weight, somatic cluster
v4         : 0.22   : sleep weight, somatic cluster
Ksom       : 2.73   : somatic half-saturation
som0       : 1.020  : healthy somatic raw score
som_floor  : 1.30   : healthy somatic HAM-A
// ---- adverse effects -------------------------------------------------------------------------------
ktol_nau   : 1.40   : nausea tolerance rate
krec_nau   : 0.16   : nausea pool recovery
ktol_dizz  : 0.55   : dizziness tolerance rate
krec_dizz  : 0.09   : dizziness pool recovery
ktol_h1    : 0.55   : H1 sedation tolerance rate
krec_h1    : 0.06   : H1 pool recovery
ks_sex     : 0.22   : sexual dysfunction build-up rate
kr_sex     : 0.10   : sexual dysfunction resolution rate
ktol_act   : 1.10   : activation tolerance rate
krec_act   : 0.14   : activation pool recovery
kwt        : 0.030  : weight gain rate (kg/d per unit drive)
kwt_off    : 0.010  : weight loss rate (1/d)
ae_sed_h1  : 2.20   : H1 contribution to sedation
ae_sed_bz  : 1.60   : benzodiazepine contribution to sedation
ae_sed_pgb : 1.10   : alpha2delta contribution to sedation
ae_dizz_pgb : 1.50  : alpha2delta contribution to dizziness
ae_dizz_bz : 0.80   : benzodiazepine contribution to dizziness
ae_nau     : 2.00   : SERT contribution to nausea
ae_act     : 1.40   : NE contribution to activation
ae_sex     : 1.00   : sexual dysfunction weight
w_ae_sed   : 0.45   : sedation weight in AE burden
w_ae_dizz  : 0.35   : dizziness weight in AE burden
w_ae_nau   : 0.50   : nausea weight in AE burden
w_ae_act   : 0.25   : activation weight in AE burden
w_ae_sex   : 0.15   : sexual dysfunction weight in AE burden
h_drop0    : 0.0016 : baseline dropout hazard (1/d)
b_ae       : 0.0070 : AE contribution to dropout hazard
// ---- comorbid depression --------------------------------------------------------------------------
km_in      : 0.05   : MADRS state formation (1/d)
km_out     : 0.05   : MADRS state resolution (1/d)
km_phi     : 0.80   : weight of z on the depressive state
m5ht       : 0.75   : 5-HT acceleration of depressive resolution
mne        : 0.35   : NE acceleration of depressive resolution
mcbt       : 0.45   : CBT acceleration of depressive resolution
madrs_scale : 8.0   : MADRS points per unit state above 1
madrs_floor : 3.0   : healthy MADRS
// ---- normalising constants and covariates ------------------------------------------------------------
phi_healthy  : 0.98484849 : Phi at the DIS = 0 attractor
sglu_healthy : 1.00000000 : S_glu at the DIS = 0 attractor
DIS        : 1.0    : disease severity (1 = typical trial patient)
CBT        : 0.0    : cognitive behavioural therapy indicator (0/1)
$END

$CMT @annotated
ESC_GUT  : escitalopram depot (ug)
ESC_C    : escitalopram central (ug)
ESC_P    : escitalopram peripheral (ug)
ESC_E    : escitalopram brain effect site (ng/mL)
VEN_GUT  : venlafaxine depot (ug)
VEN_C    : venlafaxine central (ug)
ODV_C    : O-desmethylvenlafaxine central (ug)
DLX_GUT  : duloxetine depot (ug)
DLX_C    : duloxetine central (ug)
PGB_GUT  : pregabalin depot (ug)
PGB_C    : pregabalin central (ug)
BZD_GUT  : benzodiazepine depot (ug)
BZD_C    : benzodiazepine central (ug)
BZD_E    : benzodiazepine brain effect site (ng/mL)
BUS_GUT  : buspirone depot (ug)
BUS_C    : buspirone central (ug)
PP1_C    : 1-pyrimidinylpiperazine central (ug)
QTP_GUT  : quetiapine depot (ug)
QTP_C    : quetiapine central (ug)
NQT_C    : norquetiapine central (ug)
SHT      : extracellular forebrain 5-HT (fraction of healthy)
AUTO     : 5-HT1A somatodendritic autoreceptor availability
NE       : extracellular noradrenaline (fraction of healthy)
A2AUTO   : alpha2-adrenergic autoreceptor availability
BDNF     : BDNF / plasticity signal
CPFC     : prefrontal control capacity (denominator of Phi)
EAMY     : amygdala / BNST drive (numerator of Phi)
TRAF     : alpha2delta trafficking effect state
RA2      : alpha2/alpha3 GABA-A pool (anxiolysis)
RA1      : alpha1 GABA-A pool (sedation)
DEPEND   : benzodiazepine adaptation state
CRH      : hypothalamic CRH
ACTH     : pituitary ACTH
CORT     : cortisol
GR       : glucocorticoid receptor function
SNS      : sympathetic tone
AUTON    : autonomic symptom load
SLEEPD   : sleep-continuity deficit
WORRY    : worry engine
EXPECT   : expectancy state
FLUCT    : enrolment-peak / regression-to-the-mean component
RNAU     : nausea tolerance pool
RDIZZ    : dizziness tolerance pool
RH1      : H1 sedation tolerance pool
SEXD     : sexual dysfunction
WT       : weight change (kg)
RACT     : activation tolerance pool
MADRSS   : depressive symptom state
CUMHAZ   : cumulative dropout hazard
VISCUE   : study-visit cue (dose amt = 1 at each visit)
$END

$MAIN
// ---------------------------------------------------------------------------
// Initial conditions: the UNTREATED attractor as a function of DIS.
//
// The attractor is the fixed point of the drug-free system.  It cannot be
// written in closed form, so it is supplied here as a cubic in DIS fitted to
// the Python reference's relaxation over DIS in [0.40, 1.75] (worst relative
// error 6.2e-3, on SLEEPD).  Outside that range, run a drug-free burn-in of
// >= 400 days instead of trusting the polynomial.
// ---------------------------------------------------------------------------
double d1 = DIS;
double d2 = DIS*DIS;
double d3 = d2*DIS;

@@INITBLOCK@@

FLUCT_0 = FLUCT0;

$ODE
// ===========================================================================
// 1. PHARMACOKINETICS
// ===========================================================================
double Cesc  = ESC_C / V1_esc;
double Cescp = ESC_P / V2_esc;
dxdt_ESC_GUT = -ka_esc * ESC_GUT;
dxdt_ESC_C   = ka_esc*ESC_GUT*F_esc - CL_esc*Cesc - Q_esc*(Cesc - Cescp);
dxdt_ESC_P   = Q_esc*(Cesc - Cescp);
dxdt_ESC_E   = keo_esc*(Cesc - ESC_E);

double Cven = VEN_C / V_ven;
double Codv = ODV_C / V_odv;
double CLm  = fm2d6 * CL_ven;
dxdt_VEN_GUT = -ka_ven * VEN_GUT;
dxdt_VEN_C   = ka_ven*VEN_GUT*F_ven - CL_ven*Cven;
dxdt_ODV_C   = mw_ratio_odv*CLm*Cven - CL_odv*Codv;

double Cdlx = DLX_C / V_dlx;
dxdt_DLX_GUT = -ka_dlx * DLX_GUT;
dxdt_DLX_C   = ka_dlx*DLX_GUT*F_dlx - CL_dlx*Cdlx;

double CLpgb = CL_pgb * (crcl / 100.0);
double Cpgb  = PGB_C / V_pgb;
dxdt_PGB_GUT = -ka_pgb * PGB_GUT;
dxdt_PGB_C   = ka_pgb*PGB_GUT*F_pgb - CLpgb*Cpgb;

double Cbzd = BZD_C / V_bzd;
dxdt_BZD_GUT = -ka_bzd * BZD_GUT;
dxdt_BZD_C   = ka_bzd*BZD_GUT*F_bzd - CL_bzd*Cbzd;
dxdt_BZD_E   = keo_bzd*(Cbzd - BZD_E);

double Cbus = BUS_C / V_bus;
double Cpp1 = PP1_C / V_pp;
dxdt_BUS_GUT = -ka_bus * BUS_GUT;
dxdt_BUS_C   = ka_bus*BUS_GUT*F_bus - CL_bus*Cbus;
dxdt_PP1_C   = f_pp*ka_bus*BUS_GUT - CL_pp*Cpp1;

double Cqtp = QTP_C / V_qtp;
double Cnqt = NQT_C / V_nqt;
dxdt_QTP_GUT = -ka_qtp * QTP_GUT;
dxdt_QTP_C   = ka_qtp*QTP_GUT*F_qtp - CL_qtp*Cqtp;
dxdt_NQT_C   = f_nqt*CL_qtp*Cqtp - CL_nqt*Cnqt;

// ===========================================================================
// 2. RECEPTOR OCCUPANCY  (competitive, multi-ligand: X/(1+X))
// ===========================================================================
double Xsert = ESC_E/EC50_sert_esc + (Cven + Codv)/EC50_sert_ven + Cdlx/EC50_sert_dlx;
double occ_sert = Xsert/(1.0 + Xsert);
double Xnet = (Cven + Codv)/EC50_net_ven + Cdlx/EC50_net_dlx + Cnqt/EC50_net_nqt;
double occ_net = Xnet/(1.0 + Xnet);
double occ_a2d = Cpgb/(Cpgb + EC50_a2d);
double Cbz_eq  = BZD_E * BZD_POT;
double occ_bz  = Cbz_eq/(Cbz_eq + EC50_bz);
double occ_h1  = Cqtp/(Cqtp + EC50_h1);
double occ_1a  = Cbus/(Cbus + EC50_1a);
double occ_a2adr = Cpp1/(Cpp1 + EC50_pp1);

// ===========================================================================
// 3. SEROTONIN -- occupancy is fast, the GATE is slow
// ===========================================================================
double ex_sht = (SHT > 1.0) ? (SHT - 1.0) : 0.0;
double fire = (1.0/(1.0 + gamma_auto*AUTO*ex_sht)) * (1.0 - ia_bus*occ_1a*AUTO);
dxdt_SHT  = kin_sht*fire - kout_sht*SHT*(sert_floor + (1.0 - sert_floor)*(1.0 - occ_sert));
dxdt_AUTO = -kdes_auto*AUTO*ex_sht + krec_auto*(1.0 - AUTO);
double sht_post = SHT + ia_bus_post*occ_1a;
double ex_shtp  = (sht_post > 1.0) ? (sht_post - 1.0) : 0.0;

// ===========================================================================
// 4. NORADRENALINE
// ===========================================================================
double ex_ne = (NE > 1.0) ? (NE - 1.0) : 0.0;
double fire_ne = (1.0/(1.0 + gamma_a2*A2AUTO*ex_ne)) * (1.0 + kpp1*occ_a2adr);
dxdt_NE     = kin_ne*fire_ne - kout_ne*NE*(net_floor + (1.0 - net_floor)*(1.0 - occ_net));
dxdt_A2AUTO = -kdes_a2*A2AUTO*ex_ne + krec_a2*(1.0 - A2AUTO);

// ===========================================================================
// 5. CORTISOL ENTERS EVERYWHERE THROUGH A SATURATING TRANSFORM
//    Without it, CRH -> cortisol -> GR down -> weaker feedback -> more CRH has
//    loop gain > 1 and the axis diverges.  This was an actual defect found
//    during construction, not a cosmetic choice.
// ===========================================================================
double exc = (CORT > 1.0) ? (CORT - 1.0) : 0.0;
double ex_cort = exc/(1.0 + 0.5*exc);

// ===========================================================================
// 6. PLASTICITY, PREFRONTAL CONTROL, AMYGDALA
// ===========================================================================
double bdnf_in = (w_5ht*sht_post + w_ne*NE + w_tonic + cbt_bdnf*CBT)*(1.0 - kcort_bdnf*ex_cort);
if (bdnf_in < 0.0) bdnf_in = 0.0;
dxdt_BDNF = kb*bdnf_in - kdb*BDNF;

double hill = BDNF*BDNF/(BDNF*BDNF + 1.0);
dxdt_CPFC = kp_pfc*hill*(cpfc_max - CPFC) - kl_pfc*CPFC*(1.0 + kc_pfc*ex_cort + kdis_pfc*DIS);

double sleepd = (SLEEPD > 0.0) ? SLEEPD : 0.0;
double ex_crh = (CRH > 1.0) ? (CRH - 1.0) : 0.0;
double amy_drive = 1.0 + kcrh_amy*ex_crh + kslp_amy*sleepd + kdis_amy*DIS;
double amy_room  = (EAMAX - EAMY > 0.0) ? (EAMAX - EAMY) : 0.0;
dxdt_EAMY = ke_in*amy_drive*amy_room/(EAMAX - 1.0)
            - ke_out*EAMY*(1.0 + a5ht_amy*ex_shtp + acbt_amy*CBT);

dxdt_TRAF = ktraf*(occ_a2d - TRAF);

// ===========================================================================
// 7. THE FOUR FACTORS AND Phi
// ===========================================================================
double ex_amy = (EAMY > 1.0) ? (EAMY - 1.0) : 0.0;
double sglu = (1.0 + kglu_stress*ex_amy + kglu_cort*ex_cort)*(1.0 - emax_pgb*TRAF);
if (sglu < 0.02) sglu = 0.02;

dxdt_RA1 = -ktol1*occ_bz*RA1 + krec1*(1.0 - RA1);
dxdt_RA2 = -ktol2*occ_bz*RA2 + krec2*(1.0 - RA2);
dxdt_DEPEND = kdep_on*(occ_bz - DEPEND);

double igaba_base = igaba_base0 - kdis_gaba*DIS;
double igaba = igaba_base*(1.0 + emax_a2*RA2*(occ_bz - wadapt*DEPEND));
if (igaba < 0.05) igaba = 0.05;

double cpfc_eff = CPFC*(1.0 + e_ex*EXPECT);
if (cpfc_eff < 0.05) cpfc_eff = 0.05;

double phi   = (EAMY*sglu)/(cpfc_eff*igaba);
double phi_n = phi/phi_healthy;

// bounded excess gain: EVERY downstream consequence reads z, not phi_n
double exphi = phi_n - 1.0;
double z = (exphi > 0.0) ? (exphi/(1.0 + exphi/ZMAX)) : exphi;
if (z < -0.90) z = -0.90;
if (z > ZMAX)  z = ZMAX;
double zp = (z > 0.0) ? z : 0.0;

// ===========================================================================
// 8. HPA AXIS
// ===========================================================================
dxdt_CRH  = kcrh*(1.0 + kcrh_phi*zp) - kcout*CRH*(1.0 + kgr_fb*GR*CORT);
dxdt_ACTH = ka_acth*CRH - ke_acth*ACTH;
dxdt_CORT = kc_cort*ACTH - kec_cort*CORT;
dxdt_GR   = kgr_in*(1.0 - GR) - kgr_dn*ex_cort*GR;

// ===========================================================================
// 9. AUTONOMIC, SLEEP AND WORRY
// ===========================================================================
dxdt_SNS = ks_in*(1.0 + ks_phi*zp + w_ne_sns*ex_ne) - ks_out*SNS;

double gaba_rel = (igaba/igaba_base - 1.0);
if (gaba_rel < 0.0) gaba_rel = 0.0;
dxdt_AUTON = kau_in*SNS - kau_out*AUTON*(1.0 + b_gaba*gaba_rel + b_pgb*TRAF);

double sl_room = (SLMAX - SLEEPD > 0.0) ? (SLMAX - SLEEPD) : 0.0;
dxdt_SLEEPD = ksl_in*(ksl_phi*zp + w_cort_sl*ex_cort)*sl_room
              - ksl_out*SLEEPD*(1.0 + h_h1*occ_h1*RH1 + h_gaba*occ_bz*RA1 + h_pgb*TRAF);

dxdt_WORRY = kw_in*(1.0 + kw_phi*zp)*(1.0 + rho_worry*WORRY/(kw_m + WORRY))
             - kw_out*WORRY*(1.0 + acbt_worry*CBT + a5ht_worry*ex_shtp);

// ===========================================================================
// 10. TRIAL MACHINERY -- the fifth clock
// ===========================================================================
double kexp_on = -log(1.0 - dvisit)*kvis;
dxdt_VISCUE = -kvis*VISCUE;
dxdt_EXPECT = kexp_on*VISCUE*(emax_exp - EXPECT) - kexp_dec*EXPECT;
dxdt_FLUCT  = -kfl*FLUCT;

// ===========================================================================
// 11. ADVERSE EFFECTS
// ===========================================================================
dxdt_RNAU  = -ktol_nau*occ_sert*RNAU + krec_nau*(1.0 - RNAU);
dxdt_RDIZZ = -ktol_dizz*TRAF*RDIZZ + krec_dizz*(1.0 - RDIZZ);
dxdt_RH1   = -ktol_h1*occ_h1*RH1 + krec_h1*(1.0 - RH1);
dxdt_SEXD  = ks_sex*occ_sert*(1.0 - SEXD) - kr_sex*SEXD;
dxdt_RACT  = -ktol_act*ex_ne*RACT + krec_act*(1.0 - RACT);
dxdt_WT    = kwt*(TRAF + 0.8*occ_h1) - kwt_off*WT;

double sedation  = ae_sed_h1*occ_h1*RH1 + ae_sed_bz*occ_bz*RA1 + ae_sed_pgb*TRAF*RDIZZ;
double dizziness = ae_dizz_pgb*TRAF*RDIZZ + ae_dizz_bz*occ_bz*RA1;
double nausea    = ae_nau*occ_sert*RNAU;
double activation = ae_act*ex_ne*RACT;
double ae_burden = w_ae_sed*sedation + w_ae_dizz*dizziness + w_ae_nau*nausea
                   + w_ae_act*activation + w_ae_sex*ae_sex*SEXD;

// ===========================================================================
// 12. COMORBID DEPRESSION AND DROPOUT HAZARD
// ===========================================================================
dxdt_MADRSS = km_in*(1.0 + km_phi*zp)
              - km_out*MADRSS*(1.0 + m5ht*ex_shtp + mne*ex_ne + mcbt*CBT);
dxdt_CUMHAZ = h_drop0 + b_ae*ae_burden;

$TABLE
double Cesc_o  = ESC_C/V1_esc;
double Cven_o  = VEN_C/V_ven;
double Codv_o  = ODV_C/V_odv;
double Cdlx_o  = DLX_C/V_dlx;
double Cpgb_o  = PGB_C/V_pgb;
double Cbzd_o  = BZD_C/V_bzd;
double Cbus_o  = BUS_C/V_bus;
double Cqtp_o  = QTP_C/V_qtp;
double Cnqt_o  = NQT_C/V_nqt;

double Xs = ESC_E/EC50_sert_esc + (Cven_o + Codv_o)/EC50_sert_ven + Cdlx_o/EC50_sert_dlx;
double OCC_SERT = Xs/(1.0 + Xs);
double Xn = (Cven_o + Codv_o)/EC50_net_ven + Cdlx_o/EC50_net_dlx + Cnqt_o/EC50_net_nqt;
double OCC_NET  = Xn/(1.0 + Xn);
double OCC_A2D  = Cpgb_o/(Cpgb_o + EC50_a2d);
double CbzeqT   = BZD_E*BZD_POT;
double OCC_BZ   = CbzeqT/(CbzeqT + EC50_bz);
double OCC_H1   = Cqtp_o/(Cqtp_o + EC50_h1);
double OCC_1A   = Cbus_o/(Cbus_o + EC50_1a);

double SHTP = SHT + ia_bus_post*OCC_1A;
double excT = (CORT > 1.0) ? (CORT - 1.0) : 0.0;
double EXCORT = excT/(1.0 + 0.5*excT);
double exAmy = (EAMY > 1.0) ? (EAMY - 1.0) : 0.0;
double SGLU = (1.0 + kglu_stress*exAmy + kglu_cort*EXCORT)*(1.0 - emax_pgb*TRAF);
if (SGLU < 0.02) SGLU = 0.02;
double IGB0 = igaba_base0 - kdis_gaba*DIS;
double IGABA = IGB0*(1.0 + emax_a2*RA2*(OCC_BZ - wadapt*DEPEND));
if (IGABA < 0.05) IGABA = 0.05;
double CPFCE = CPFC*(1.0 + e_ex*EXPECT);
if (CPFCE < 0.05) CPFCE = 0.05;
double PHI   = (EAMY*SGLU)/(CPFCE*IGABA);
double PHI_N = PHI/phi_healthy;

double SLD = (SLEEPD > 0.0) ? SLEEPD : 0.0;
double psy_raw = c1*WORRY + c2*PHI_N + c3*SLD;
double psy_eff = (psy_raw > psy0) ? (psy_raw - psy0) : 0.0;
double HAMA_PSY = psy_floor + (28.0 - psy_floor)*psy_eff/(psy_eff + Kpsy);
double som_raw = v1*AUTON + v2*SNS + v3*(SGLU/sglu_healthy) + v4*SLD;
double som_eff = (som_raw > som0) ? (som_raw - som0) : 0.0;
double HAMA_SOM = som_floor + (28.0 - som_floor)*som_eff/(som_eff + Ksom);
double HAMA = HAMA_PSY + HAMA_SOM + FLUCT;

double SEDATION  = ae_sed_h1*OCC_H1*RH1 + ae_sed_bz*OCC_BZ*RA1 + ae_sed_pgb*TRAF*RDIZZ;
double DIZZINESS = ae_dizz_pgb*TRAF*RDIZZ + ae_dizz_bz*OCC_BZ*RA1;
double NAUSEA    = ae_nau*OCC_SERT*RNAU;
double exNE      = (NE > 1.0) ? (NE - 1.0) : 0.0;
double ACTIVATION = ae_act*exNE*RACT;
double AE_BURDEN = w_ae_sed*SEDATION + w_ae_dizz*DIZZINESS + w_ae_nau*NAUSEA
                   + w_ae_act*ACTIVATION + w_ae_sex*ae_sex*SEXD;
double MADRS = madrs_floor + madrs_scale*((MADRSS > 1.0) ? (MADRSS - 1.0) : 0.0);
double SURVIVAL = exp(-CUMHAZ);

$CAPTURE @annotated
Cesc_o   : escitalopram plasma (ng/mL)
Cven_o   : venlafaxine plasma (ng/mL)
Codv_o   : ODV plasma (ng/mL)
Cdlx_o   : duloxetine plasma (ng/mL)
Cpgb_o   : pregabalin plasma (ng/mL)
Cbzd_o   : benzodiazepine plasma (ng/mL)
Cqtp_o   : quetiapine plasma (ng/mL)
OCC_SERT : serotonin transporter occupancy
OCC_NET  : noradrenaline transporter occupancy
OCC_A2D  : alpha2delta occupancy
OCC_BZ   : benzodiazepine-site occupancy
OCC_H1   : histamine H1 occupancy
SHTP     : postsynaptic 5-HT signal
SGLU     : S_glu, numerator factor
IGABA    : I_gaba, denominator factor
CPFCE    : C_pfc including expectancy
PHI      : corticolimbic gain
PHI_N    : corticolimbic gain, normalised to healthy
HAMA     : HAM-A total (0-56)
HAMA_PSY : HAM-A psychic cluster
HAMA_SOM : HAM-A somatic cluster
MADRS    : comorbid depressive symptoms
SEDATION : sedation burden
DIZZINESS : dizziness burden
NAUSEA   : nausea burden
ACTIVATION : activation / jitteriness burden
AE_BURDEN : weighted adverse-effect burden
SURVIVAL : probability of remaining in the trial
$END
