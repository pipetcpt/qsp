## =============================================================================
##  Severe Traumatic Brain Injury (sTBI) -- mrgsolve QSP model
##  중증 외상성 뇌손상 정량적 시스템 약리학 모델
## =============================================================================
##
##  47 ODEs.  Craniospinal mechanics (Monro-Kellie on an exponential
##  pressure-volume curve), cerebral haemodynamics with a saturating
##  autoregulatory sigmoid, CO2 reactivity with choroid-plexus bicarbonate
##  adaptation, two-region brain water with REGIONAL reflection coefficients,
##  regional pericontusional perfusion, the excitotoxic/mitochondrial cascade,
##  neuroinflammation with a biphasic barrier, haematoma expansion with
##  fibrinolysis, four biomarkers, and the PK/PD of propofol, thiopental,
##  noradrenaline, mannitol, hypertonic saline and tranexamic acid.
##
##  THE STRUCTURAL IDEA.  ICP is not a state with its own biology.  It is the
##  residual of a volume balance:
##
##      dICP/dt * [PVI/(ICP*ln10) + Ca - dVv/dICP]
##            = dCa/dt*(MAP-ICP) + Ca*dMAP/dt + (If - Io - Qevd)
##              + Jw_intact + Jw_injured + dVhaematoma/dt + mass
##
##  Everything a drug can do, it does by changing one term on the right.
##
##  PROVENANCE.  Every equation and parameter here mirrors, line for line, the
##  pure-Python reference implementation `tbi_reference_model.py`, whose
##  committed output `tbi_reference_output.txt` is the source of every number
##  quoted in README.md.  The Python version is the one that has been run; if
##  the two ever disagree, the Python file is the specification.
##
##  Units: time min | volume mL | pressure mmHg | flow mL/min
##         resistance mmHg*min/mL | osmolality mOsm/kg | CBF per 100 g
##
##  Usage:
##      library(mrgsolve); library(dplyr)
##      mod <- mread("tbi_mrgsolve_model", "tbi_mrgsolve_model.R")
##      out <- mod %>% mrgsim(end = 4320, delta = 5)
##      plot(out, ICP + CPP + PbtO2 + CBF_pen ~ time)
## =============================================================================

$PROB
# Severe traumatic brain injury: intracranial pressure, perfusion and metabolism
# 중증 외상성 뇌손상 QSP 모델

$PARAM @annotated
// ---- craniospinal mechanics ------------------------------------------------
PVI       : 26.0  : Pressure-volume index (mL)
PVI_crani : 95.0  : PVI after decompressive craniectomy (mL)
CRANI_T   : -1    : Time of craniectomy, negative = never (min)
Pvs       : 6.0   : Sagittal sinus pressure (mmHg)
V_csf_min : 62.0  : CSF volume at which the cisterns are effaced (mL)
MASS_RATE : 0.0   : Extra intracranial mass growth, a probe (mL/min)

// ---- cerebrovascular resistance and compliance -----------------------------
Ran       : 0.0800 : Nominal arteriolar resistance (mmHg*min/mL)
Rpv       : 0.0200 : Post-capillary resistance (mmHg*min/mL)
Rv        : 0.0120 : Cerebral venous resistance (mmHg*min/mL)
Rf        : 63.53  : CSF formation resistance (mmHg*min/mL)
Ro        : 11.43  : CSF outflow resistance (mmHg*min/mL)
R_evd     : 8.0    : External ventricular drain resistance (mmHg*min/mL)
Can       : 0.200  : Nominal arteriolar compliance (mL/mmHg)
Ca_lo_f   : 0.45   : Minimum compliance as a fraction of Can
Ca_hi_f   : 1.55   : Maximum compliance as a fraction of Can
k_sig     : 0.070  : Steepness of the compliance sigmoid
tau_aut   : 0.333  : Autoregulatory time constant (min)
G_aut     : 3.00   : Autoregulatory loop gain
AUTOREG   : 1.0    : Multiplier on G_aut, 1 intact and 0 abolished
Vv_max    : 26.0   : Collapsible venous blood volume (mL)
Pv_half   : 25.0   : ICP at which half the venous buffer is gone (mmHg)
kv        : 9.0    : Steepness of venous collapse (mmHg)
kR        : 486.72 : Poiseuille constant, back-solved at the operating point
q0        : 697.07 : Baseline CBF (mL/min)
Pc0       : 32.23  : Baseline capillary pressure (mmHg)
hyd0      : 22.23  : Baseline capillary-to-tissue gradient (mmHg)

// ---- CO2 -------------------------------------------------------------------
fco2_lo   : 0.35  : CO2 response floor
fco2_hi   : 2.05  : CO2 response ceiling
co2_center: 45.0  : CO2 sigmoid centre (mmHg)
co2_k     : 10.0  : CO2 sigmoid width (mmHg)
k_co2_cal : 0.45  : Direct CO2 action on achievable caliber
HCO3_n    : 24.0  : Normal CSF bicarbonate (mmol/L)
tau_hco3  : 360.0 : Bicarbonate adaptation time constant (min)
hco3_slope: 0.40  : Bicarbonate change per mmHg PaCO2
PACO2     : 38.0  : Set arterial PaCO2 (mmHg)
SAO2      : 0.98  : Arterial oxygen saturation
PAO2      : 95.0  : Arterial oxygen tension (mmHg)

// ---- metabolism and oxygen -------------------------------------------------
CMRO2n    : 3.30  : Normal CMRO2 (mL O2/100g/min)
Q10       : 2.30  : Temperature coefficient
met_alpha : 0.90  : Flow-metabolism coupling exponent
OEF_max   : 0.85  : Maximum oxygen extraction fraction
CBFn100   : 49.79 : Nominal CBF (mL/100g/min)
brain_g   : 1400  : Brain mass (g)
Hb        : 12.5  : Haemoglobin (g/dL)
k_pbto2   : 0.75  : PbtO2 per unit venous PO2
k_diff_ed : 1.4   : Oedema penalty on O2 diffusion
k_diff_mv : 1.6   : Microvascular collapse penalty
k_pen_res : 0.55  : Pericontusional resistance penalty

// ---- brain water -----------------------------------------------------------
V_int0    : 950.0   : Intact-region brain water (mL)
V_inj0    : 170.0   : Injured-region brain water (mL)
LpS_int   : 0.00190 : Osmotic conductance, intact (mL/min/mmHg)
LpS_inj   : 0.00260 : Osmotic conductance, injured (mL/min/mmHg)
Lh_int    : 0.0035  : Hydraulic conductance, intact (mL/min/mmHg)
Lh_inj    : 0.0300  : Hydraulic conductance, injured (mL/min/mmHg)
mmHg_mOsm : 19.3    : mmHg per mOsm/kg
sig_prot_open : 0.95 : Oncotic reflection lost at full barrier opening
k_glym_int: 0.0130  : ISF clearance, intact (1/min)
k_glym_inj: 0.0062  : ISF clearance, injured (1/min)
k_leak_int: 0.00030 : Solute equilibration across an intact barrier (1/min)
k_leak_bbb: 0.0100  : Solute leak per unit barrier opening (1/min)
k_mann_bbb: 0.0055  : Mannitol entry per unit barrier opening (1/min)
k_idio    : 0.115   : Idiogenic osmole generation
k_idio_clr: 0.0055  : Idiogenic osmole clearance (1/min)
Osm_n     : 293.0   : Reference osmolality (mOsm/kg)

// ---- systemic osmotic and renal --------------------------------------------
Na_n      : 140.0   : Reference plasma sodium (mmol/L)
V_ecf0    : 15000   : Reference ECF volume (mL)
Glu_p     : 8.0     : Plasma glucose (mmol/L)
BUN       : 5.0     : Plasma urea (mmol/L)
GFR       : 0.110   : Glomerular filtration rate (L/min)
k_na_corr : 0.0055  : Proportional renal sodium correction
k_vol_corr: 0.0012  : Proportional renal volume correction
Mann_V1   : 7000    : Mannitol central volume (mL)
Mann_V2   : 8000    : Mannitol peripheral volume (mL)
Mann_Q    : 0.90    : Mannitol intercompartmental clearance (L/min)
Osm_aki   : 320.0   : Plasma osmolality above which GFR degrades
FLUID     : 1.6     : Maintenance isotonic fluid (mL/min)

// ---- mass lesion and coagulopathy -------------------------------------------
k_hem_exp : 0.0072  : Haematoma expansion rate constant
tau_hem   : 210.0   : Haematoma expansion window (min)
k_hem_res : 0.00016 : Haematoma resorption (1/min)
k_fib     : 0.0140  : Fibrinolysis activation
k_fib_off : 0.0035  : Fibrinolysis decay (1/min)

// ---- excitotoxic cascade ------------------------------------------------------
Glu_n     : 2.0     : Resting ECF glutamate (uM)
k_glu_rel : 520.0   : Glutamate release on ionic failure
k_glu_upt : 0.34    : Glutamate uptake (1/min)
K_n       : 3.0     : Resting ECF potassium (mmol/L)
k_k_rel   : 52.0    : Potassium efflux on ionic failure
k_k_upt   : 0.30    : Potassium reuptake (1/min)
k_ca_in   : 0.055   : Calcium influx
Km_glu    : 28.0    : Half-maximal glutamate for calcium influx (uM)
k_ca_out  : 0.052   : Calcium extrusion (1/min)
Ca_n      : 0.10    : Resting intracellular calcium (arbitrary)
k_mito    : 0.0135  : Mitochondrial injury rate
k_mito_rep: 0.0021  : Mitochondrial repair rate (1/min)
k_ros     : 0.030   : ROS generation
k_ros_clr : 0.055   : ROS clearance (1/min)
k_lac     : 0.155   : Lactate generation
k_lac_clr : 0.030   : Lactate clearance (1/min)
Lac_n     : 1.6     : Resting brain lactate (mmol/L)
Glc_br_n  : 1.8     : Resting brain glucose (mmol/L)
k_glc_in  : 0.055   : Brain glucose transport (1/min)
k_glc_use : 0.075   : Brain glucose utilisation (1/min)

// ---- tissue fate and spreading depolarisation ----------------------------------
k_pen2core: 0.00105 : Penumbra to core conversion
k_pen_rec : 0.00034 : Penumbra recovery
k_sd_base : 0.55    : Spreading depolarisation rate scale
sd_cost   : 0.42    : CMRO2 surcharge per unit SD burden

// ---- neuroinflammation and barrier ----------------------------------------------
k_micro   : 0.0075   : Microglial activation
d_micro   : 0.00095  : Microglial decay (1/min)
k_cyto    : 0.0125   : Cytokine production
d_cyto    : 0.0042   : Cytokine decay (1/min)
k_neut    : 0.0035   : Neutrophil recruitment
d_neut    : 0.0016   : Neutrophil decay (1/min)
k_mmp     : 0.0060   : MMP-9 production
d_mmp     : 0.0021   : MMP-9 decay (1/min)
k_bbb_mech_off : 0.0038 : Mechanical barrier resealing (1/min)
k_bbb_infl: 0.00092  : Inflammatory barrier opening
k_bbb_seal: 0.00055  : Inflammatory barrier resealing (1/min)

// ---- biomarkers ---------------------------------------------------------------
k_gfap    : 145.0    : GFAP release
ke_gfap   : 0.00048  : GFAP elimination, t1/2 ~24 h (1/min)
k_uchl1   : 62.0     : UCH-L1 release
ke_uchl1  : 0.00165  : UCH-L1 elimination, t1/2 ~7 h (1/min)
k_nfl     : 9.0      : NfL release
ke_nfl    : 0.000023 : NfL elimination, t1/2 ~21 d (1/min)
k_s100b   : 17.0     : S100B release
ke_s100b  : 0.0116   : S100B elimination, t1/2 ~1 h (1/min)

// ---- temperature -----------------------------------------------------------------
T_n       : 37.0   : Normal core temperature (C)
k_fever   : 1.35   : Maximum cytokine-driven fever (C)
tau_T     : 95.0   : Fever time constant (min)
tau_cool  : 42.0   : Active cooling time constant (min)
T_TARGET  : -1     : Target temperature, negative = no active control (C)

// ---- drugs ------------------------------------------------------------------------
prop_V1   : 15.9  : Propofol V1 (L)
prop_V2   : 32.6  : Propofol V2 (L)
prop_V3   : 203.0 : Propofol V3 (L)
prop_Cl1  : 1.94  : Propofol elimination clearance (L/min)
prop_Cl2  : 1.75  : Propofol Q2 (L/min)
prop_Cl3  : 0.79  : Propofol Q3 (L/min)
prop_ke0  : 0.26  : Propofol effect-site rate constant (1/min)
prop_Emax : 0.55  : Maximum propofol CMRO2 suppression
prop_EC50 : 3.40  : Propofol CMRO2 EC50 (ug/mL)
prop_MAPmax : 26.0 : Maximum propofol MAP reduction (mmHg)
prop_MAP50: 3.10  : Propofol MAP EC50 (ug/mL)
PROP_RATE : 3.0   : Propofol infusion (mg/kg/h)
thio_V1   : 28.0  : Thiopental V1 (L)
thio_V2   : 110.0 : Thiopental V2 (L)
thio_Cl   : 0.21  : Thiopental clearance (L/min)
thio_Q    : 1.05  : Thiopental Q (L/min)
thio_Emax : 0.55  : Maximum thiopental CMRO2 suppression
thio_EC50 : 26.0  : Thiopental CMRO2 EC50 (ug/mL)
thio_MAPmax : 24.0 : Maximum thiopental MAP reduction (mmHg)
thio_MAP50: 32.0  : Thiopental MAP EC50 (ug/mL)
THIO_RATE : 0.0   : Thiopental infusion (mg/kg/h)
supp_floor: 0.60  : Combined anaesthetic CMRO2 suppression ceiling
ne_ke0    : 0.55  : Noradrenaline effect-site rate constant (1/min)
ne_Emax   : 55.0  : Maximum noradrenaline MAP effect (mmHg)
ne_EC50   : 0.300 : Noradrenaline EC50 (ug/kg/min)
NE_RATE   : 0.0   : Noradrenaline infusion (ug/kg/min)
txa_V1    : 12.0  : TXA V1 (L)
txa_V2    : 15.0  : TXA V2 (L)
txa_Cl    : 0.115 : TXA clearance (L/min)
txa_Q     : 0.16  : TXA Q (L/min)
txa_IC50  : 8.0   : TXA antifibrinolytic IC50 (mg/L)

// ---- clinical controls ----------------------------------------------------------
MAP_BASE  : 88.0  : Baseline mean arterial pressure (mmHg)
EVD_OPEN  : 0     : External ventricular drain open (0/1)
EVD_SET   : 12.0  : EVD drainage threshold (mmHg)
SEV       : 0.62  : Injury severity (0-1)
AGE       : 42    : Age (years)
WT        : 75    : Body weight (kg)

$CMT @annotated
Pic     : Intracranial pressure (mmHg)
x_aut   : Autoregulatory tone (dimensionless)
MAPs    : Mean arterial pressure (mmHg)
V_csf   : CSF volume (mL)
HCO3    : CSF bicarbonate (mmol/L)
V_int   : Intact-region brain water (mL)
V_inj   : Injured-region brain water (mL)
Osm_int : Intact-region osmolality (mOsm/kg)
Osm_inj : Injured-region osmolality, mannitol excluded (mOsm/kg)
Mann_inj: Mannitol inside the lesion (mmol)
Na_ecf  : Exchangeable sodium (mmol)
V_ecf   : Extracellular fluid volume (mL)
Mann_c  : Mannitol, central (mmol)
Mann_p  : Mannitol, peripheral (mmol)
V_hem   : Haematoma volume (mL)
Fib     : Fibrinolytic activity (dimensionless)
Glu     : ECF glutamate (uM)
K_ec    : ECF potassium (mmol/L)
Ca_i    : Intracellular calcium (arbitrary)
MitoD   : Mitochondrial dysfunction (0-1)
ROS     : Reactive oxygen species (arbitrary)
Lac     : Brain lactate (mmol/L)
Glc_br  : Brain glucose (mmol/L)
F_core  : Non-viable tissue fraction
F_pen   : Penumbral tissue fraction
Micro   : Microglial activation (0-1)
Cyto    : IL-1beta/TNF composite (arbitrary)
Neut    : Neutrophil infiltration (arbitrary)
MMP9    : MMP-9 activity (arbitrary)
BBB_mech: Mechanical barrier disruption (0-1)
BBB_infl: Inflammatory barrier opening (0-1)
GFAP    : Serum GFAP (arbitrary)
UCHL1   : Serum UCH-L1 (arbitrary)
NfL     : Serum neurofilament light (arbitrary)
S100B   : Serum S100B (arbitrary)
Temp    : Core temperature (C)
Prop1   : Propofol central (mg)
Prop2   : Propofol peripheral 1 (mg)
Prop3   : Propofol peripheral 2 (mg)
Prop_e  : Propofol effect site (ug/mL)
Thio1   : Thiopental central (mg)
Thio2   : Thiopental peripheral (mg)
NE_e    : Noradrenaline effect site (ug/kg/min)
TXA_c   : Tranexamic acid central (mg)
TXA_p   : Tranexamic acid peripheral (mg)
D_icp   : Cumulative ICP dose (mmHg*h)
D_cpp   : Cumulative CPP deficit dose (mmHg*h)

$GLOBAL
#define LN10 2.302585092994046

// Numerically safe logistic
double sigf(double z) {
  if (z < -60.0) return 0.0;
  if (z >  60.0) return 1.0;
  return 1.0/(1.0+exp(-z));
}
double clampf(double v, double lo, double hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}
// x^e guarded against the transient negative intermediates a solver can produce
double pwf(double b, double e) { return b > 0.0 ? pow(b, e) : 0.0; }

// Smooth min(a,b): accurate to <0.1% when the arguments differ two-fold, and
// continuously differentiable, which min() is not.
double softminf(double a, double b) {
  if (a <= 0.0 || b <= 0.0) return 0.0;
  return a / pow(1.0 + pow(a/b, 8.0), 0.125);
}

// Arteriolar compliance as a function of tone, with the CO2 caliber multiplier.
// CO2 enters the model twice on purpose: through the flow the servo is aiming
// at, and through the CEILING the smooth muscle can reach.  Without the second
// action a maximally dilated vessel could not respond to hyperventilation.
double Ca_of(double x, double Can, double lo_f, double hi_f, double ks,
             double m_co2) {
  double lo = lo_f*Can*m_co2, hi = hi_f*Can*m_co2;
  return lo + (hi-lo)*sigf(x/ks);
}
double dCa_dx_of(double x, double Can, double lo_f, double hi_f, double ks,
                 double m_co2) {
  double lo = lo_f*Can*m_co2, hi = hi_f*Can*m_co2;
  double s = sigf(x/ks);
  return (hi-lo)*s*(1.0-s)/ks;
}
// Inverse oxyhaemoglobin dissociation curve (Severinghaus)
double pvo2_of(double sv) {
  double ss = clampf(sv, 0.02, 0.995);
  return clampf(pow(23400.0/(1.0/ss - 1.0), 1.0/3.0), 1.0, 120.0);
}

$MAIN
// ---- initial conditions scale with injury severity -------------------------
if (NEWIND <= 1) {
  Pic_0      = 10.0;
  x_aut_0    = 0.0;
  MAPs_0     = MAP_BASE;
  // Compensatory reserve is not the same in every patient and is not a free
  // parameter: effacement of the basal cisterns is a Marshall/Rotterdam CT
  // criterion and the strongest single CT predictor of intracranial
  // hypertension, while cerebral atrophy is why an eighty-year-old tolerates a
  // subdural that would kill a twenty-year-old.  Both are the SAME quantity --
  // how much CSF there was to give up -- so both enter here.
  V_csf_0    = clampf(140.0 - 62.0*SEV + 0.85*(AGE - 40.0), 70.0, 190.0);
  HCO3_0     = HCO3_n;
  V_int_0    = V_int0;
  V_inj_0    = V_inj0;
  Osm_int_0  = Osm_n;
  Osm_inj_0  = Osm_n;
  Na_ecf_0   = Na_n * V_ecf0 / 1000.0;
  V_ecf_0    = V_ecf0;
  V_hem_0    = 14.0 * SEV;      // clot volume already present on arrival
  Fib_0      = 0.10;
  Glu_0      = Glu_n;
  K_ec_0     = K_n;
  Ca_i_0     = Ca_n;
  MitoD_0    = 0.02 + 0.10*SEV;
  ROS_0      = 0.02 + 0.10*SEV;
  Lac_0      = Lac_n;
  Glc_br_0   = Glc_br_n;
  // The impact-instant injury: the core is already dead on arrival, and the
  // penumbra is the entire target of everything that follows.
  F_core_0   = 0.055 * SEV;
  F_pen_0    = 0.300 * SEV;
  Micro_0    = 0.04 * SEV;
  Cyto_0     = 0.02 * SEV;
  Neut_0     = 0.01 * SEV;
  MMP9_0     = 0.05 * SEV;
  BBB_mech_0 = 0.55 * SEV;
  BBB_infl_0 = 0.02;
  Temp_0     = 37.0;
}

$ODE
// ===========================================================================
//  1. TISSUE FATE AND METABOLIC DEMAND
// ===========================================================================
double Fc  = clampf(F_core, 0.0, 1.0);
double Fp  = clampf(F_pen,  0.0, 1.0-Fc);
double viable = 1.0 - Fc;

double Ce_prop = Prop_e;
double C_thio  = Thio1/thio_V1;
double supp_raw = prop_Emax*Ce_prop/(prop_EC50+Ce_prop)
                + thio_Emax*C_thio/(thio_EC50+C_thio);
double supp = supp_floor*(1.0 - exp(-supp_raw/supp_floor));

double q10f = pow(Q10, (Temp-37.0)/10.0);
double sd_drive = sigf((K_ec-5.4)/0.85) * (Fp/0.30);
double SD_rate  = k_sd_base*12.0*sd_drive;          // events per hour

double CMRO2_dem = CMRO2n*q10f*(1.0-supp)*(viable*0.86+0.14)
                 * (1.0 + sd_cost*sd_drive);

// ===========================================================================
//  2. CO2 -> PERIVASCULAR pH -> EFFECTIVE PaCO2
// ===========================================================================
double PaCO2e = PACO2 * (HCO3_n / (HCO3 > 6.0 ? HCO3 : 6.0));
double f_co2  = fco2_lo + (fco2_hi-fco2_lo)*sigf((PaCO2e-co2_center)/co2_k);
double m_co2  = clampf(1.0 + k_co2_cal*(f_co2-1.0), 0.65, 1.35);

double f_met   = pwf(CMRO2_dem/CMRO2n, met_alpha);
double q_target = CBFn100 * brain_g/100.0 * f_co2 * f_met;

// ===========================================================================
//  3. VASCULAR MECHANICS AND THE STARLING RESISTOR
// ===========================================================================
int crani = (CRANI_T >= 0.0 && SOLVERTIME >= CRANI_T) ? 1 : 0;
double PVIe = crani ? PVI_crani : PVI;

double Ca    = Ca_of(x_aut, Can, Ca_lo_f, Ca_hi_f, k_sig, m_co2);
double dCadx = dCa_dx_of(x_aut, Can, Ca_lo_f, Ca_hi_f, k_sig, m_co2);
double Va    = Ca*(MAPs - Pic); if (Va < 1.0) Va = 1.0;
double Ra    = kR*Can*Can/(Va*Va);

// Cerebral venous outflow is ICP-referenced, not sinus-referenced: the bridging
// veins are collapsible and ICP normally exceeds sinus pressure.
double Pdown = Pvs + (Pic - Pvs)*sigf((Pic - Pvs)/1.0);
double Rtot  = Ra + Rpv + Rv;
double q     = (MAPs - Pdown)/Rtot; if (q < 1.0) q = 1.0;
double Pc    = MAPs - q*Ra;
double CBF100 = q/(brain_g/100.0);
double CPP    = MAPs - Pic;

double sv_    = sigf((Pic - Pv_half)/kv);
double Vv     = Vv_max*(1.0 - sv_);
double dVv_dP = -Vv_max*sv_*(1.0-sv_)/kv;

// ===========================================================================
//  4. OXYGEN: GLOBAL AND REGIONAL
// ===========================================================================
double CaO2 = 1.34*Hb*SAO2 + 0.003*PAO2;
double DO2  = CBF100*CaO2/100.0;
double CMRO2_act = softminf(CMRO2_dem, DO2*OEF_max);
double OEF  = CMRO2_act/(DO2 > 1e-6 ? DO2 : 1e-6);
double SvO2 = SAO2*(1.0-OEF);

// The pericontusional compartment sits behind a much higher local resistance.
// This is why a jugular bulb catheter can read normally over dying tissue.
double ed_inj_frac = (V_inj - V_inj0)/V_inj0; if (ed_inj_frac < 0.0) ed_inj_frac = 0.0;
double mv_col  = sigf((Pic-28.0)/6.0)*0.7 + 0.3*Fp;
double BBB     = clampf(BBB_mech + BBB_infl, 0.0, 1.0);
double pen_res = 1.0 + k_pen_res*(mv_col + 0.6*BBB + 5.0*ed_inj_frac
                                  + 0.6*sd_drive);
double CBF_pen = CBF100/pen_res;
double DO2_pen = CBF_pen*CaO2/100.0;
double CMRO2_pen = softminf(CMRO2_dem, DO2_pen*OEF_max);
double r_pen  = CMRO2_pen/(CMRO2_dem > 1e-6 ? CMRO2_dem : 1e-6);
double OEF_pen = CMRO2_pen/(DO2_pen > 1e-6 ? DO2_pen : 1e-6);
double SvO2_pen = SAO2*(1.0-OEF_pen);

double PvO2_pen = pvo2_of(SvO2_pen);
double Dfac  = 1.0/(1.0 + k_diff_ed*ed_inj_frac + k_diff_mv*mv_col);
double PbtO2 = k_pbto2*PvO2_pen*Dfac;

// Reduced pump CAPACITY and frank IONIC FAILURE are different thresholds.
double E_pump  = sigf((r_pen-0.62)/0.045);
double ionfail = sigf((0.55-r_pen)/0.035);
double LPR = 15.0*(1.0 + 6.0*pwf(1.0-r_pen,1.3) + 2.5*(MitoD > 0 ? MitoD : 0));

// ===========================================================================
//  5. OSMOTIC STATE
// ===========================================================================
double Na_p    = Na_ecf*1000.0/(V_ecf > 1000.0 ? V_ecf : 1000.0);
double Mann_pl = Mann_c/Mann_V1*1000.0;
double Osm_p   = 2.0*Na_p + Glu_p + BUN + Mann_pl;
double Mann_inj_mM = Mann_inj/V_inj*1000.0;
double Osm_int_t = Osm_int;
double Osm_inj_t = Osm_inj + Mann_inj_mM;

// ===========================================================================
//  6. AUTOREGULATORY TONE
// ===========================================================================
// Loss of autoregulation is loss of the ability to CHANGE caliber, not a
// command to return to some neutral caliber.  A vasoparalytic vessel keeps
// whatever tone it had and thereafter behaves as a passive distensible tube,
// so the servo target is blended toward the CURRENT tone.  At AUTOREG = 0 this
// gives dx/dt = 0 exactly, which is what makes the abolished-autoregulation
// arm a fair comparator rather than a different patient.
double G = G_aut;
double x_reg = clampf(-G*(q - q_target)/q_target, -8.0, 8.0);
double x_ss = AUTOREG*x_reg + (1.0 - AUTOREG)*x_aut;
double dx = (x_ss - x_aut)/tau_aut;
dxdt_x_aut = dx;

// ===========================================================================
//  7. ARTERIAL PRESSURE, WITH THE CUSHING REFLEX AS AN EMERGENT LOOP
// ===========================================================================
double map_ne   = ne_Emax*NE_e/(ne_EC50+NE_e);
double map_prop = prop_MAPmax*Ce_prop/(prop_MAP50+Ce_prop);
double map_thio = thio_MAPmax*C_thio/(thio_MAP50+C_thio);
double cush     = 35.0*sigf((30.0-CPP)/5.0);
double MAP_tgt  = MAP_BASE + map_ne + cush - map_prop - map_thio;
dxdt_MAPs = (MAP_tgt - MAPs)/0.50;
double HR = 88.0 - 0.62*cush + 6.0*(Temp-37.0);

// ===========================================================================
//  8. CSF: A FINITE BUFFER
// ===========================================================================
double If = (Pc > Pic) ? (Pc - Pic)/Rf : 0.0;
double csf_avail = sigf((V_csf - V_csf_min)/6.0);
double Io = ((Pic > Pvs) ? (Pic - Pvs)/Ro : 0.0) * csf_avail;
double evd = 0.0;
if (EVD_OPEN > 0.5 && Pic > EVD_SET)
  evd = (Pic - EVD_SET)/R_evd * sigf((V_csf - V_csf_min + 8.0)/5.0);
dxdt_V_csf = If - Io - evd;

dxdt_HCO3 = (HCO3_n + hco3_slope*(PACO2-40.0) - HCO3)/tau_hco3;

// ===========================================================================
//  9. BRAIN WATER: TWO REGIONS, TWO REFLECTION COEFFICIENTS
// ---------------------------------------------------------------------------
//  Water follows the osmole, so a HIGHER plasma osmolality drives water OUT.
//  Vasogenic oedema is not a rise in capillary pressure: an intact barrier
//  holds net filtration at zero DESPITE a 22 mmHg gradient, because its
//  reflection coefficient is ~1.  Opening the barrier does not push harder --
//  it stops pulling back, and the pre-existing gradient does the rest.
// ===========================================================================
double sig_int = 0.97;
double sig_inj = 0.97*(1.0 - 0.92*BBB);
double sig_prot_inj = 1.0 - sig_prot_open*BBB;
double hyd_int = Lh_int*((Pc - Pic) - 1.0*hyd0);
double hyd_inj = Lh_inj*((Pc - Pic) - sig_prot_inj*hyd0);

double Jw_int = LpS_int*mmHg_mOsm*sig_int*(Osm_int_t - Osm_p) + hyd_int
              - k_glym_int*(V_int - V_int0);
double Jw_inj = LpS_inj*mmHg_mOsm*sig_inj*(Osm_inj_t - Osm_p) + hyd_inj
              - k_glym_inj*(V_inj - V_inj0);
dxdt_V_int = Jw_int;
dxdt_V_inj = Jw_inj;

double idio_int = k_idio*ionfail*0.25;
double idio_inj = k_idio*ionfail*1.00 + 0.020*Fp;
double leak_int = k_leak_int*(Osm_p - Osm_int_t);
double leak_inj = k_leak_bbb*BBB*(Osm_p - Mann_pl - Osm_inj);
dxdt_Osm_int = idio_int + leak_int - k_idio_clr*(Osm_int_t - Osm_n)
             - Osm_int_t*Jw_int/V_int;
dxdt_Osm_inj = idio_inj + leak_inj - k_idio_clr*0.55*(Osm_inj - Osm_n)
             - Osm_inj*Jw_inj/V_inj;

// The rebound mechanism in one line: mannitol that crossed the broken barrier
// stays there, and when plasma mannitol clears the gradient goes NEGATIVE.
double mann_flux = k_mann_bbb*BBB*(Mann_pl - Mann_inj_mM)*V_inj/1000.0;
dxdt_Mann_inj = mann_flux - 0.0022*Mann_inj;

// ===========================================================================
//  10. SYSTEMIC SODIUM, VOLUME AND MANNITOL
// ===========================================================================
double na_in  = FLUID*0.154;   // isotonic maintenance carries its salt
double vol_in = FLUID;
// Hypertonic saline and mannitol arrive as dosing events into Na_ecf / Mann_c.
double gfr = GFR*(1.0 - 0.75*sigf((Osm_p - Osm_aki)/7.0));
double na_out = na_in + k_na_corr*(Na_p - Na_n)*V_ecf/1000.0;
double mann_filt = Mann_pl*gfr;
double osm_diuresis = mann_filt/0.35;
double urine = FLUID + osm_diuresis + k_vol_corr*(V_ecf - V_ecf0);
dxdt_Na_ecf = na_in - na_out;
dxdt_V_ecf  = vol_in - urine;

double Cp_mann = Mann_p/Mann_V2*1000.0;
dxdt_Mann_c = -mann_filt - Mann_Q*(Mann_pl - Cp_mann) - mann_flux;
dxdt_Mann_p =  Mann_Q*(Mann_pl - Cp_mann);

// ===========================================================================
//  11. HAEMATOMA AND COAGULOPATHY
// ===========================================================================
double txa_c   = TXA_c/txa_V1;
double txa_inh = txa_c/(txa_IC50 + txa_c);
dxdt_V_hem = k_hem_exp*Fib*exp(-SOLVERTIME/tau_hem)*(MAPs/88.0)*SEV*60.0
           - k_hem_res*V_hem;
dxdt_Fib   = k_fib*SEV*exp(-SOLVERTIME/240.0)*(1.0-txa_inh) - k_fib_off*Fib;

// ===========================================================================
//  12. EXCITOTOXIC AND METABOLIC CASCADE
// ===========================================================================
double sd_burst = SD_rate/12.0;
dxdt_Glu = k_glu_rel*ionfail*viable + 34.0*sd_burst
         - k_glu_upt*(0.25+0.75*E_pump)*(Glu - Glu_n);
dxdt_K_ec = k_k_rel*ionfail*viable*0.045 + 0.55*sd_burst
          - k_k_upt*(0.20+0.80*E_pump)*(K_ec - K_n);
double dglu = (Glu - Glu_n) > 0 ? (Glu - Glu_n) : 0.0;
dxdt_Ca_i = k_ca_in*dglu/(Km_glu + dglu) - k_ca_out*E_pump*(Ca_i - Ca_n);
double dca = (Ca_i - Ca_n) > 0 ? (Ca_i - Ca_n) : 0.0;
dxdt_ROS  = k_ros*(dca*0.8 + MitoD + 0.4*Neut) - k_ros_clr*ROS;
dxdt_MitoD = k_mito*dca*(1.0+1.6*ROS)*(1.0-MitoD) - k_mito_rep*MitoD*E_pump;
dxdt_Lac  = k_lac*(1.0-r_pen > 0 ? 1.0-r_pen : 0.0)*12.0 - k_lac_clr*(Lac - Lac_n);
dxdt_Glc_br = k_glc_in*(Glu_p*0.32 - Glc_br)
            - k_glc_use*(0.35+0.65*(1.0-E_pump))*Glc_br + 0.128;

// ===========================================================================
//  13. TISSUE FATE
//  Two thresholds, because one cannot do both jobs: a graded arm from the
//  electrical-failure threshold, and a steep arm at ionic-pump failure.
// ===========================================================================
double kill = k_pen2core*(pwf(MitoD,1.8)*9.0
                          + 12.0*pwf(0.80-r_pen,1.2)
                          + 55.0*pwf(0.35-r_pen,2.0)
                          + 0.35*sd_burst);
double rescue = k_pen_rec*pwf(E_pump,3.0)*sigf((PbtO2-20.0)/4.0);
dxdt_F_pen  = -(kill + rescue)*F_pen;
dxdt_F_core =  kill*F_pen;

// ===========================================================================
//  14. NEUROINFLAMMATION AND THE BIPHASIC BARRIER
// ===========================================================================
double damp = kill*F_pen*400.0 + 0.30*F_core + 0.5*ROS;
dxdt_Micro = k_micro*damp*(1.0-Micro) - d_micro*Micro;
dxdt_Cyto  = k_cyto*Micro - d_cyto*Cyto;
dxdt_Neut  = k_neut*Cyto - d_neut*Neut;
dxdt_MMP9  = k_mmp*(Cyto + 0.8*Neut) - d_mmp*MMP9;
dxdt_BBB_mech = -k_bbb_mech_off*BBB_mech;
dxdt_BBB_infl =  k_bbb_infl*MMP9*(1.0-BBB_infl) - k_bbb_seal*BBB_infl;

// ===========================================================================
//  15. BIOMARKERS -- one release process, four elimination half-lives
// ===========================================================================
double rel = kill*F_pen;
dxdt_GFAP  = k_gfap*rel*100.0  - ke_gfap*GFAP;
dxdt_UCHL1 = k_uchl1*rel*100.0 - ke_uchl1*UCHL1;
dxdt_NfL   = k_nfl*rel*100.0   - ke_nfl*NfL;
dxdt_S100B = k_s100b*(rel*100.0 + 0.02*BBB) - ke_s100b*S100B;

// ===========================================================================
//  16. TEMPERATURE
// ===========================================================================
if (T_TARGET >= 0.0) {
  dxdt_Temp = (T_TARGET - Temp)/tau_cool;
} else {
  dxdt_Temp = (T_n + k_fever*Cyto/(0.6+Cyto) - Temp)/tau_T;
}

// ===========================================================================
//  17. DRUG PHARMACOKINETICS
// ===========================================================================
double rate_prop = PROP_RATE*WT/60.0;      // mg/min
double C1 = Prop1/prop_V1, C2 = Prop2/prop_V2, C3 = Prop3/prop_V3;
dxdt_Prop1 = rate_prop - prop_Cl1*C1 - prop_Cl2*(C1-C2) - prop_Cl3*(C1-C3);
dxdt_Prop2 = prop_Cl2*(C1-C2);
dxdt_Prop3 = prop_Cl3*(C1-C3);
dxdt_Prop_e = prop_ke0*(C1 - Prop_e);

double rate_thio = THIO_RATE*WT/60.0;
double T1 = Thio1/thio_V1, T2 = Thio2/thio_V2;
dxdt_Thio1 = rate_thio - thio_Cl*T1 - thio_Q*(T1-T2);
dxdt_Thio2 = thio_Q*(T1-T2);

dxdt_NE_e = ne_ke0*(NE_RATE - NE_e);

double Tc = TXA_c/txa_V1, Tp = TXA_p/txa_V2;
dxdt_TXA_c = -txa_Cl*Tc - txa_Q*(Tc-Tp);
dxdt_TXA_p =  txa_Q*(Tc-Tp);

// ===========================================================================
//  18. THE MONRO-KELLIE BALANCE
//  Solved for dICP/dt rather than book-kept, so the constraint cannot drift.
// ===========================================================================
double denom = PVIe/((Pic > 1.0 ? Pic : 1.0)*LN10) + Ca - dVv_dP;
double dCa_dt = dCadx*dx;
double net = If - Io - evd + Jw_int + Jw_inj + MASS_RATE
           + dxdt_V_hem + dCa_dt*(MAPs - Pic) + Ca*dxdt_MAPs;
dxdt_Pic = net/denom;

// cumulative pressure-time doses, in mmHg*hour
dxdt_D_icp = (Pic > 20.0 ? (Pic-20.0) : 0.0)/60.0;
dxdt_D_cpp = (CPP <  60.0 ? (60.0-CPP) : 0.0)/60.0;

$TABLE
double ICP   = Pic;
double CPP_o = MAPs - Pic;
double MAP_o = MAPs;
double SjvO2 = SvO2*100.0;
double CBV   = Va + Vv;
double V_ed  = (V_int - V_int0) + (V_inj - V_inj0);
double osm_gap = Mann_pl;
double Cp_prop = Prop1/prop_V1;
double Cp_thio = Thio1/thio_V1;
// A cheap PRx surrogate.  The reference implementation solves the quasi-steady
// tone by bisection; here the blend of the passive response (always positive)
// and the regulated response is approximated from the local sigmoid slope.
double dVa_pass = Ca;
double dVa_reg  = Ca - dCadx*(G_aut*AUTOREG)*3.2*(MAPs-Pic)/q_target;
double dVa_eff  = dVa_pass*0.20 + dVa_reg*0.80;
double PRx = tanh(4.0*dVa_eff/(denom > 0.05 ? denom : 0.05));

$CAPTURE @annotated
ICP     : Intracranial pressure (mmHg)
CPP_o   : Cerebral perfusion pressure (mmHg)
MAP_o   : Mean arterial pressure (mmHg)
CBF100  : Global cerebral blood flow (mL/100g/min)
CBF_pen : Pericontusional cerebral blood flow (mL/100g/min)
PbtO2   : Brain tissue oxygen tension (mmHg)
SjvO2   : Jugular bulb oxygen saturation (%)
LPR     : Microdialysis lactate/pyruvate ratio
CMRO2_act : Cerebral metabolic rate for oxygen (mL O2/100g/min)
PRx     : Pressure reactivity index
CBV     : Cerebral blood volume (mL)
Va      : Arterial-arteriolar volume (mL)
Vv      : Venous blood volume (mL)
Ca      : Arteriolar compliance (mL/mmHg)
Pc      : Capillary pressure (mmHg)
V_ed    : Excess brain water (mL)
Osm_p   : Plasma osmolality (mOsm/kg)
Na_p    : Plasma sodium (mmol/L)
osm_gap : Osmolar gap, i.e. plasma mannitol (mOsm/kg)
Mann_inj_mM : Mannitol trapped inside the lesion (mmol/L)
csf_avail : CSF buffer remaining (0-1)
E_pump  : Na/K-ATPase capacity (0-1)
r_pen   : Pericontusional supply adequacy
ionfail : Ionic pump failure (0-1)
BBB     : Blood-brain barrier opening (0-1)
SD_rate : Spreading depolarisation rate (per hour)
PaCO2e  : Effective PaCO2 seen by the vessel (mmHg)
f_co2   : CO2 flow multiplier
cush    : Cushing reflex pressure contribution (mmHg)
HR      : Heart rate (per min)
If      : CSF formation (mL/min)
Io      : CSF absorption (mL/min)
evd     : External ventricular drain flow (mL/min)
urine   : Urine output (mL/min)
gfr     : Glomerular filtration rate (L/min)
Cp_prop : Propofol plasma concentration (ug/mL)
Cp_thio : Thiopental plasma concentration (ug/mL)
kill    : Penumbra-to-core conversion rate (per min)
denom   : Effective craniospinal compliance (mL/mmHg)

## =============================================================================
##  SCENARIOS
## -----------------------------------------------------------------------------
##  Everything below is ordinary R and runs after mread().  The numbers these
##  scenarios produce are reproduced, with commentary, in tbi_reference_output.txt
##  by the Python reference implementation.
## =============================================================================
##
## library(mrgsolve); library(dplyr); library(ggplot2)
## mod <- mread("tbi", "tbi_mrgsolve_model.R")
##
## ---------------------------------------------------------------------------
## Scenario 1 -- natural history: no ICP-directed therapy
## ---------------------------------------------------------------------------
## s1 <- mod %>% param(SEV = 0.55, PROP_RATE = 3) %>% mrgsim(end = 4320, delta = 5)
##
## ---------------------------------------------------------------------------
## Scenario 2 -- tier 1: external ventricular drain plus hypertonic saline
##   3% NaCl 250 mL contains 128 mmol of sodium; give it into Na_ecf over 30 min,
##   and the accompanying volume into V_ecf.
## ---------------------------------------------------------------------------
## hts <- bind_rows(
##   ev(time = 360, amt = 128.2, cmt = "Na_ecf", tinf = 30),
##   ev(time = 360, amt = 250,   cmt = "V_ecf",  tinf = 30))
## s2 <- mod %>% param(SEV = 0.55, EVD_OPEN = 1, EVD_SET = 12) %>%
##   ev(hts) %>% mrgsim(end = 1440, delta = 2)
##
## ---------------------------------------------------------------------------
## Scenario 3 -- mannitol 0.5 g/kg, and the rebound that follows it
##   0.5 g/kg of a 75 kg patient = 37.5 g = 205.9 mmol, in ~188 mL of 20% solution
## ---------------------------------------------------------------------------
## mann <- bind_rows(
##   ev(time = 360, amt = 205.9, cmt = "Mann_c", tinf = 20),
##   ev(time = 360, amt = 188,   cmt = "V_ecf",  tinf = 20))
## s3 <- mod %>% param(SEV = 0.55) %>% ev(mann) %>% mrgsim(end = 1440, delta = 2)
## # watch Mann_inj_mM keep rising while osm_gap falls: that is the rebound
##
## ---------------------------------------------------------------------------
## Scenario 4 -- hyperventilation is a loan
##   PaCO2 30 for 24 h, then back to 38.  Watch HCO3 adapt, the ICP benefit
##   evaporate, and the restoration produce a rebound because the bicarbonate
##   has not come back yet.
## ---------------------------------------------------------------------------
## s4 <- mod %>% param(SEV = 0.55) %>%
##   ev(bind_rows(ev(time =   60, amt = 0, cmt = 0, evid = 2, PACO2 = 30),
##                ev(time = 1500, amt = 0, cmt = 0, evid = 2, PACO2 = 38))) %>%
##   mrgsim(end = 1800, delta = 2, recsort = 3)
##
## ---------------------------------------------------------------------------
## Scenario 5 -- CPP augmentation, and the sign flip
##   Identical noradrenaline in a patient with intact and with abolished
##   autoregulation.  The ICP goes in opposite directions.
## ---------------------------------------------------------------------------
## s5a <- mod %>% param(SEV = 0.55, NE_RATE = 0.18, AUTOREG = 1.0) %>% mrgsim(end = 720)
## s5b <- mod %>% param(SEV = 0.55, NE_RATE = 0.18, AUTOREG = 0.0) %>% mrgsim(end = 720)
##
## ---------------------------------------------------------------------------
## Scenario 6 -- barbiturate coma
## ---------------------------------------------------------------------------
## s6 <- mod %>% param(SEV = 0.62, THIO_RATE = 5, PROP_RATE = 6) %>% mrgsim(end = 2880)
##
## ---------------------------------------------------------------------------
## Scenario 7 -- therapeutic hypothermia with a controlled rewarm
## ---------------------------------------------------------------------------
## s7 <- mod %>% param(SEV = 0.62) %>%
##   ev(bind_rows(ev(time =  120, amt = 0, cmt = 0, evid = 2, T_TARGET = 33.5),
##                ev(time = 3000, amt = 0, cmt = 0, evid = 2, T_TARGET = 37.0))) %>%
##   mrgsim(end = 4320, delta = 5, recsort = 3)
##
## ---------------------------------------------------------------------------
## Scenario 8 -- decompressive craniectomy, early versus late
## ---------------------------------------------------------------------------
## s8a <- mod %>% param(SEV = 0.62, CRANI_T =  480) %>% mrgsim(end = 4320, delta = 5)
## s8b <- mod %>% param(SEV = 0.62, CRANI_T = 1800) %>% mrgsim(end = 4320, delta = 5)
##
## ---------------------------------------------------------------------------
## Scenario 9 -- tranexamic acid, and the CRASH-3 time window
##   1 g over 10 min then 1 g over 8 h, given at 1 h and at 6 h post-injury
## ---------------------------------------------------------------------------
## txa <- function(t0) bind_rows(ev(time = t0, amt = 1000, cmt = "TXA_c", tinf = 10),
##                               ev(time = t0 + 10, amt = 1000, cmt = "TXA_c", tinf = 480))
## s9a <- mod %>% param(SEV = 0.62) %>% ev(txa(60))  %>% mrgsim(end = 2880, delta = 5)
## s9b <- mod %>% param(SEV = 0.62) %>% ev(txa(360)) %>% mrgsim(end = 2880, delta = 5)
##
## ---------------------------------------------------------------------------
## Scenario 10 -- a secondary insult, and why the hour matters
##   Transient hypotension and hypoxaemia together, delivered early and late.
## ---------------------------------------------------------------------------
## insult <- function(t0) bind_rows(
##   ev(time = t0,      amt = 0, cmt = 0, evid = 2, MAP_BASE = 55, SAO2 = 0.75),
##   ev(time = t0 + 30, amt = 0, cmt = 0, evid = 2, MAP_BASE = 88, SAO2 = 0.98))
## s10a <- mod %>% param(SEV = 0.55) %>% ev(insult(60))   %>% mrgsim(end = 2880, recsort = 3)
## s10b <- mod %>% param(SEV = 0.55) %>% ev(insult(2160)) %>% mrgsim(end = 2880, recsort = 3)
##
## ---------------------------------------------------------------------------
## Scenario 11 -- compensatory reserve: the point of the whole model
##   Grow a mass at 0.5 mL/min from different starting points and record how
##   many millilitres go in before ICP reaches 30.  ICP barely moves for hours
##   while the answer to that question collapses.
## ---------------------------------------------------------------------------
## s11 <- mod %>% param(SEV = 0.50, MASS_RATE = 0.5) %>% mrgsim(end = 900, delta = 1)
##
## ---------------------------------------------------------------------------
## Scenario 12 -- the emergent Lassen curve
##   Hold MAP and let the system settle; CBF is flat between CPP 60 and 130 and
##   linear outside it, and nothing in the model was told to do that.
## ---------------------------------------------------------------------------
## lassen <- lapply(seq(40, 170, 10), function(m)
##   mod %>% param(SEV = 0, MAP_BASE = m, PROP_RATE = 0) %>%
##     mrgsim(end = 25, delta = 25) %>% as_tibble() %>% slice_tail(n = 1) %>%
##     mutate(MAP_set = m)) %>% bind_rows()
## # plot(lassen$CPP_o, lassen$CBF100, type = "b")
