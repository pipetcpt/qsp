## =============================================================================
##  Carbon Monoxide Poisoning -- Quantitative Systems Pharmacology model
##  mrgsolve implementation, 45 ODEs
## =============================================================================
##
##  ORGANISING IDEA
##  ---------------
##  Carbon monoxide creates two occupancies with two different time constants:
##
##    (i)  carboxyhaemoglobin -- fast (t1/2 = 313 min on air, 72 min on a mask),
##         responds to oxygen within minutes, and is the ONLY one measured; and
##    (ii) CO on myoglobin and on the reduced a3 haem of cytochrome c oxidase --
##         slow (off-rate tau ~455 min), protected by tissue PO2 in a way COHb is
##         not, and the pool that actually writes the injury.
##
##  Everything downstream -- a reassuring pulse oximeter, a normal COHb in a
##  patient who will deteriorate, delayed neurological sequelae weeks later, and
##  two randomised trials of hyperbaric oxygen that disagree -- follows from
##  monitoring and treating (i) while the disease is being written by (ii).
##
##  STRUCTURE
##  ---------
##   1. Environment          C_env (ppm), source rate, air changes
##   2. Pulmonary exchange   Coburn-Forster-Kane mass balance on COHb
##   3. Oxygen transport     content, Haldane left shift, Fick, autoregulation
##   4. Tissue CO            4 perfusion-limited compartments + myoglobin
##   5. Terminal oxidase     CO on reduced a3, competitively protected by O2
##   6. Energetics           brain, watershed (globus pallidus), myocardium
##   7. Reoxygenation injury xanthine oxidase -> ROS -> lipid peroxidation
##   8. Innate               beta-2 integrin -> neutrophil -> myeloperoxidase
##   9. Adaptive             MBP adduct -> autoreactive clone -> demyelination
##                           (a bistable switch with a computable threshold)
##  10. Organ injury         necrosis, oedema, ICP, troponin, EF, CK, creatinine
##  11. Cyanide              fire smoke; converges on the same terminal oxidase
##  12. Drugs                O2, HBO, hydroxocobalamin, NAC, allopurinol
##  13. Fetus                placental transfer with a slow fetal compartment
##
##  CALIBRATION (see co_references.md for the sources of each anchor)
##  -----------------------------------------------------------------
##   * CFK constants are resting physiology, NOT fitted: VA 4.2 L/min (VE 6
##     minus 30% dead space), DLco 25 mL/min/mmHg, Haldane M = 245, Vb 5.5 L.
##     They reproduce Weaver's measured half-lives to within 2-3%:
##          model 313 min on room air   vs observed 320
##          model  72 min on a 100% O2 mask vs observed  74
##     The single fitted quantity in the whole gas-exchange block is `shunt`
##     (0.72), which folds venous admixture and mask leak into the mean
##     pulmonary capillary PO2, and it was set ONCE on the 74 min figure.
##   * Baseline COHb 0.618% (endogenous haem catabolism 0.007 mL/min plus 1.5
##     ppm ambient) against a measured normal of 0.4-0.8%.
##   * The healthy unexposed subject is an EXACT fixed point of all 45 states.
##   * DNS incidence is anchored on Weaver 2002 (NEJM 347:1057): cognitive
##     sequelae at 6 weeks 46% with normobaric oxygen, 25% with HBO x3.
##
##  KNOWN MISSES -- stated here rather than buried
##  ---------------------------------------------
##   1. CFK CANNOT REPRODUCE THE HYPERBARIC HALF-LIFE, and the reason is
##      structural, not a bad parameter.  Raising chamber pressure raises the
##      driving force P_c,O2 and the transfer resistance P_L/VA in exact
##      proportion, so P_c,O2/B -> VA and the rate constant tends to a CEILING
##      k_max = VA/(M [O2Hb] Vb).  That puts a floor of 31.3 min on the
##      half-life at any pressure whatsoever.  The observed half-life at
##      2.5-3 ATA is ~20 min, i.e. 36% BELOW a floor the equation says cannot
##      be crossed.  Either an extrapulmonary elimination route exists or the
##      Haldane ratio M is itself pressure-dependent.  The model reports 38.6
##      min at 3 ATA and is simply wrong there; no attempt has been made to
##      hide this by re-fitting `shunt` per pressure.
##   2. Cyanide is represented with a lumped tissue compartment and a single
##      rhodanese term.  It reproduces the direction and rough magnitude of the
##      CO/CN interaction but has not been calibrated against fire-victim data.
##   3. N-acetylcysteine and allopurinol are included because they are the
##      mechanistically obvious adjuncts, not because there is human outcome
##      evidence.  Their predicted effects should be read as hypotheses.
##   4. The globus pallidus watershed is one lumped compartment.  Real CO
##      lesions are patchy and the model cannot speak to their distribution.
##
##  USAGE
##  -----
##    library(mrgsolve); library(dplyr)
##    mod <- mread("co_mrgsolve_model.R")
##    # a severe accidental exposure treated with a non-rebreather mask
##    mod %>% param(ppm_fix = 2298, texp = 60, to2_start = 90, to2_stop = 450) %>%
##            mrgsim(end = 45*1440, delta = 5) %>% plot(COHb_pct + CcOb + MBPad + Cog ~ time)
##
##  Every equation in this file has been independently re-implemented in
##  Python/scipy (copoison.py) and the reference output is committed as
##  co_reference_output.txt.  That cross-check found and fixed six defects; they
##  are listed in README.md.
## =============================================================================

$PROB
# Carbon monoxide poisoning: 45-ODE QSP model
- Coburn-Forster-Kane pulmonary exchange coupled to a bistable autoimmune
  demyelination switch.
- Author: QSP disease model library
- Units: time in MINUTES throughout.

$PARAM @annotated
// ---------------- subject -------------------------------------------------
Hb        : 15.0   : Haemoglobin (g/dL)
Vb        : 5500   : Blood volume (mL)
WT        : 70     : Body weight (kg)

// ---------------- Coburn-Forster-Kane -------------------------------------
DL        : 25.0   : Pulmonary diffusing capacity for CO (mL/min/mmHg)
VA        : 4200   : Alveolar ventilation at rest (mL/min)
MHald     : 245    : Haldane ratio, CO/O2 haemoglobin affinity (-)
Vco_endo  : 0.007  : Endogenous CO production from haem oxygenase (mL/min)
shunt     : 0.72   : Pulmonary capillary efficiency, venous admixture + mask leak (-)
ppm_amb   : 1.5    : Ambient background CO (ppm)

// ---------------- exposure ------------------------------------------------
ppm_fix   : -1     : Fixed room CO; negative uses the dynamic Cenv compartment (ppm)
Rsrc      : 0      : CO source rate into the room (mL/min)
Vroom     : 30000  : Enclosed space volume (L equivalent, mL basis)
ACH       : 0.5    : Air changes per hour (1/h)
texp      : 0      : Exposure duration (min)

// ---------------- therapy -------------------------------------------------
FiO2      : 0.21   : Baseline inspired O2 fraction (-)
to2_start : 1e9    : Start of normobaric oxygen (min)
to2_stop  : 1e9    : Stop of normobaric oxygen (min)
FiO2_trt  : 0.85   : Inspired O2 fraction on a tight non-rebreather (-)
VAtrt     : 4200   : Alveolar ventilation during treatment (mL/min)
FiCO2     : 0      : Inspired CO2 fraction, carbogen (-)
thbo_start: 1e9    : First hyperbaric session start (min)
thbo_dur  : 90     : Duration of each hyperbaric session (min)
nhbo      : 1      : Number of hyperbaric sessions (-)
hbo_gap   : 480    : Interval between hyperbaric sessions (min)
ATA       : 1.0    : Chamber pressure during a session (atmospheres absolute)
hypotherm : 0      : Targeted temperature management switch, 0/1 (-)

// ---------------- tissue CO distribution ----------------------------------
kappa     : 1750   : Tissue CO tension per unit COHb fraction (milli-mmHg)
tau_br    : 9      : Brain CO equilibration time constant (min)
tau_ht    : 45     : Myocardial CO equilibration time constant (min)
tau_mu    : 320    : Skeletal muscle CO equilibration time constant (min)
tau_rest  : 120    : Splanchnic and other CO time constant (min)
Mmb       : 36     : Myoglobin Haldane ratio, CO/O2 affinity (-)

// ---------------- terminal oxidase ----------------------------------------
Kc        : 0.55   : CO dissociation constant at reduced haem a3 (mmHg)
Ko        : 0.55   : O2 tension that half-protects haem a3 by oxidising it (mmHg)
kon_cco   : 0.10   : CO on-rate onto reduced haem a3 (1/min)
koff_cco  : 0.0022 : CO off-rate from haem a3, tau ~455 min (1/min)
CcO_base  : 0.0095 : Baseline haem a3 occupancy from endogenous CO (-)
hbo_koff  : 2.6    : Fractional increase in CcO off-rate inside the chamber (-)
Kc_cn     : 6.0    : Cyanide concentration half-inhibiting cytochrome c oxidase (uM)

// ---------------- oxygen transport ----------------------------------------
P50       : 26.8   : Oxyhaemoglobin P50 (mmHg)
nHill     : 2.7    : Hill coefficient of the oxyhaemoglobin curve (-)
alphaCO2  : 0.75   : Strength of the CO-induced allosteric left shift (-)
CMRO2     : 3.30   : Cerebral metabolic rate for O2 (mL/100g/min)
CBF0      : 50     : Baseline cerebral blood flow (mL/100g/min)
tau_cbf   : 1.5    : Cerebrovascular response time constant (min)
PaCO2b   : 40     : Baseline arterial CO2 tension (mmHg)
CBFmax    : 2.4    : Maximum autoregulatory flow reserve, multiple of baseline (-)

// ---------------- energetics ----------------------------------------------
tau_atp   : 6      : Energy charge equilibration time constant (min)
Km_o2     : 12.0   : TISSUE-level O2 constant at the diffusion field edge (mmHg)
PtO2ref   : 33.6   : Reference healthy cerebral tissue PO2 (mmHg)
ATPthr    : 0.92   : Energy reserve before the oxidative cascade ignites (-)
f_wshed   : 0.55   : Watershed tissue PO2 as a fraction of mean cerebral (-)
kLac      : 0.55   : Lactate production gain (mM/min)
kLacEl    : 0.030  : Lactate elimination rate (1/min)
Lac0      : 1.0    : Baseline arterial lactate (mM)

// ---------------- reoxygenation injury ------------------------------------
kXOon     : 0.028  : Xanthine dehydrogenase to oxidase conversion rate (1/min)
kXOoff    : 0.0016 : Reversal of the conversion (1/min)
kROS      : 0.185  : ROS generation gain (1/min)
kROSel    : 0.10   : ROS clearance rate (1/min)
kNO       : 0.225  : NO / peroxynitrite generation gain (1/min)
kNOel     : 0.09   : NO / peroxynitrite clearance rate (1/min)
kGSH      : 0.0032 : Glutathione resynthesis rate (1/min)
kGSHox    : 0.0075 : Glutathione consumption by ROS (1/min)
EC50_nac  : 42     : NAC concentration half-maximally boosting glutathione (mg/L)
IC50_oxy  : 6.0    : Oxypurinol concentration half-inhibiting xanthine oxidase (mg/L)

// ---------------- innate immune amplification -----------------------------
kadh      : 0.0075 : Neutrophil adhesion gain (1/min)
kdeadh    : 0.0075 : Neutrophil de-adhesion rate (1/min)
hbo_adh   : 0.82   : Maximal HBO blockade of beta-2 integrin adhesion (-)
tau_hbo_adh: 420   : Decay of the HBO anti-adhesion effect after a session (min)
kMPO      : 0.0045 : Myeloperoxidase release gain (1/min)
kMPOel    : 0.0045 : Myeloperoxidase clearance rate (1/min)
kLPO      : 0.00114: Lipid peroxidation gain (1/min)
kLPOel    : 0.0021 : Lipid peroxide clearance rate (1/min)
wMPO      : 0.85   : Weight of myeloperoxidase relative to ROS in peroxidation (-)

// ---------------- adaptive autoimmunity, the DNS switch -------------------
kAd       : 0.00050: MBP adduct formation gain (1/min)
kSpread   : 0.00040: Epitope spreading, antigen released by demyelination (1/min)
kAdClr    : 4.0e-5 : Adduct clearance on the myelin turnover timescale (1/min)
theta     : 0.450  : Adduct burden at half-maximal clonal proliferation (-)
nT        : 6.0    : Hill coefficient of the tolerance threshold (-)
Tprol     : 0.00040: Autoreactive clone proliferation rate (1/min)
Tdeath    : 0.00019: Autoreactive clone contraction rate (1/min)
Tmax      : 1.0    : Clone carrying capacity (-)
T0        : 0.004  : Naive autoreactive clone frequency (-)
kMicro    : 0.0010 : Microglial activation gain (1/min)
kMicroEl  : 0.0012 : Microglial deactivation rate (1/min)
kDemy     : 2.2e-5 : Demyelination rate per unit clone (1/min)
kDemyDir  : 8.0e-6 : Direct microglial demyelination rate (1/min)
kRepair   : 1.5e-5 : Remyelination rate (1/min)

// ---------------- adaptive protective -------------------------------------
kHO1      : 0.0016 : Haem oxygenase-1 induction gain (1/min)
kHO1el    : 0.00035: Haem oxygenase-1 decay rate (1/min)
kHO1co    : 1.9    : Extra endogenous CO per unit HO-1 induction (-)
kHIF      : 0.006  : HIF-1 alpha stabilisation gain (1/min)
kHIFel    : 0.0016 : HIF-1 alpha degradation rate (1/min)

// ---------------- organ injury --------------------------------------------
ATPcrit   : 0.42   : Watershed energy charge below which necrosis proceeds (-)
kNec      : 0.0075 : Necrosis rate (1/min)
kEdema    : 0.020  : Cerebral oedema formation gain (1/min)
kEdemaEl  : 0.0018 : Cerebral oedema resolution rate (1/min)
ICP0      : 10     : Baseline intracranial pressure (mmHg)
kICP      : 42     : Intracranial pressure per unit oedema (mmHg)
ATPcrit_h : 0.72   : Myocardial energy charge below which troponin leaks (-)
kTrop     : 1.5    : Troponin release gain (ng/mL/min)
kTropEl   : 0.0060 : Troponin elimination rate (1/min)
kEF       : 0.85   : Ejection fraction sensitivity to energy deficit (-)
kEFrec    : 0.0016 : Ejection fraction recovery rate (1/min)
EF0       : 0.62   : Baseline left ventricular ejection fraction (-)
CK0       : 100    : Baseline creatine kinase (U/L)
MbCK      : 0.15   : Myoglobin-CO occupancy above which muscle leaks CK (-)
kCK       : 300    : Creatine kinase release gain (U/L/min)
kCKel     : 0.0009 : Creatine kinase elimination rate (1/min)
Cr0       : 0.9    : Baseline serum creatinine (mg/dL)
kCr       : 0.00025: Creatinine rise per unit myoglobin load (mg/dL/min)
kCrEl     : 0.0013 : Creatinine elimination rate (1/min)

// ---------------- outcome -------------------------------------------------
wDemy     : 0.62   : Cognitive weight of demyelination (-)
wNec      : 0.55   : Cognitive weight of necrosis (-)
tau_cog   : 4320   : Cognitive score time constant, 3 days (min)
kHaz      : 0.02   : DNS hazard accrual per unit demyelination (1/min)

// ---------------- cyanide (fire smoke) ------------------------------------
CN_rate   : 0      : Cyanide absorption rate during smoke exposure (uM/min)
CN_tau    : 25     : Duration of cyanide absorption (min)
CL_cn     : 0.55   : Cyanide clearance (L/min)
V_cn      : 25     : Cyanide central volume (L)
k_rhod    : 0.0075 : Rhodanese-mediated detoxification rate (1/min)
k_cnt     : 0.045  : Cyanide blood-tissue transfer rate (1/min)
lam_cn    : 2.2    : Cyanide tissue-to-blood partition coefficient (-)
k_ohcbl   : 0.085  : Second-order hydroxocobalamin-cyanide binding (1/uM/min)
CL_ohc    : 0.30   : Hydroxocobalamin clearance (L/min)
V_ohc     : 18     : Hydroxocobalamin volume (L)

// ---------------- adjunct drug PK -----------------------------------------
CL_nac    : 13.0   : N-acetylcysteine clearance (L/min)
V1_nac    : 12.0   : NAC central volume (L)
Q_nac     : 8.0    : NAC intercompartmental clearance (L/min)
V2_nac    : 20.0   : NAC peripheral volume (L)
CL_oxy    : 1.2    : Oxypurinol clearance (L/min)
V_oxy     : 25     : Oxypurinol volume (L)

// ---------------- fetus ---------------------------------------------------
fet_ratio : 1.8    : Fetal-to-maternal COHb ratio at equilibrium (-)
fet_tau   : 250    : Fetal COHb equilibration time constant (min)

$CMT @annotated
Cenv   : Room CO concentration (ppm)
ACO    : CO bound to haemoglobin, whole body (mL CO)
Pbr    : Brain tissue CO tension (milli-mmHg)
Pht    : Myocardial tissue CO tension (milli-mmHg)
Pmu    : Skeletal muscle CO tension (milli-mmHg)
Prest  : Splanchnic and other tissue CO tension (milli-mmHg)
MbHt   : Cardiac myoglobin CO-occupied fraction (-)
MbMu   : Skeletal myoglobin CO-occupied fraction (-)
CcOb   : Brain cytochrome c oxidase CO-inhibited fraction (-)
CcOh   : Myocardial cytochrome c oxidase CO-inhibited fraction (-)
ATPb   : Brain energy charge (-)
ATPh   : Myocardial energy charge (-)
Lac    : Arterial lactate (mM)
CBF    : Cerebral blood flow (mL/100g/min)
PaCO2  : Arterial CO2 tension (mmHg)
XO     : Xanthine oxidase converted fraction (-)
ROS    : Brain reactive oxygen species index (-)
NOx    : NO / peroxynitrite index (-)
Neut   : Adherent neutrophil index (-)
MPO    : Myeloperoxidase activity index (-)
LPO    : Lipid peroxidation product index (-)
MBPad  : Myelin basic protein adduct burden (-)
Tcell  : Autoreactive lymphocyte clone (-)
Micro  : Activated microglia index (-)
Demy   : Demyelination burden (-)
HO1    : Haem oxygenase-1 induction (-)
HIF    : HIF-1 alpha index (-)
Necb   : Watershed brain necrosis burden (-)
Edema  : Cerebral oedema index (-)
ICP    : Intracranial pressure (mmHg)
TnI    : Cardiac troponin I (ng/mL)
EF     : Left ventricular ejection fraction (-)
CK     : Creatine kinase (U/L)
Cr     : Serum creatinine (mg/dL)
Cog    : Neurocognitive composite, fraction of baseline (-)
HazD   : Cumulative DNS hazard (-)
CNbl   : Blood cyanide (uM)
CNtis  : Tissue cyanide (uM)
OHCbl  : Hydroxocobalamin (uM)
NACc   : N-acetylcysteine, central (mg/L)
NACp   : N-acetylcysteine, peripheral (mg/L)
GSH    : Brain glutathione, fraction of baseline (-)
Oxy    : Oxypurinol (mg/L)
FetCO  : Fetal COHb fraction (-)
ATPgp  : Watershed (globus pallidus) energy charge (-)

$GLOBAL
#define PH2O 47.0

// Hill function, guarded at zero
double hillf(double x, double k, double n) {
  double xx = (x > 0.0) ? x : 0.0;
  double xn = pow(xx, n);
  return xn / (xn + pow(k, n));
}
// oxyhaemoglobin saturation and its inverse
double sat_o2(double po2, double p50, double n) {
  double p = (po2 > 1e-9) ? po2 : 1e-9;
  return pow(p, n) / (pow(p, n) + pow(p50, n));
}
double po2_of_sat(double s, double p50, double n) {
  double ss = (s < 1e-9) ? 1e-9 : ((s > 1.0-1e-9) ? 1.0-1e-9 : s);
  return p50 * pow(ss/(1.0-ss), 1.0/n);
}
double dmax(double a, double b) { return (a > b) ? a : b; }
double dmin(double a, double b) { return (a < b) ? a : b; }

$MAIN
// The healthy subject must be an EXACT fixed point.  COHb starts at the value
// where endogenous production plus ambient uptake balances pulmonary loss:
//     Vco_endo = (100 F / (M (1-F)) - PIamb) / B      solved for F
double capH  = 1.34*Hb/100.0;
double Bres0 = 1.0/DL + (760.0-PH2O)/VA;
double PIamb = ppm_amb*1e-6*(760.0-PH2O);
double rhs0  = (Vco_endo*Bres0 + PIamb)*MHald/100.0;   // = F/(1-F)
double F0    = rhs0/(1.0 + rhs0);

ACO_0   = F0*capH*Vb;
Pbr_0   = kappa*F0;
Pht_0   = kappa*F0;
Pmu_0   = kappa*F0;
Prest_0 = kappa*F0;
ATPb_0  = 1.0;
ATPh_0  = 1.0;
ATPgp_0 = 1.0;
Lac_0   = Lac0;
CBF_0   = CBF0;
PaCO2_0 = PaCO2b;
Tcell_0 = T0;
ICP_0   = ICP0;
EF_0    = EF0;
CK_0    = CK0;
Cr_0    = Cr0;
Cog_0   = 1.0;
GSH_0   = 1.0;
FetCO_0 = fet_ratio*F0;

$ODE
// =========================================================================
// therapy state: which gas, at what pressure, with what ventilation
// =========================================================================
int    inO2    = (SOLVERTIME >= to2_start && SOLVERTIME < to2_stop) ? 1 : 0;
double FiO2_c  = inO2 ? FiO2_trt : FiO2;
double VA_c    = inO2 ? VAtrt    : VA;
double ata     = 1.0;
double chamber = 0.0;
double since   = 1e9;
for (int i = 0; i < (int) nhbo; ++i) {
  double t0 = thbo_start + i*hbo_gap;
  if (SOLVERTIME >= t0) { double s = SOLVERTIME - t0; if (s < since) since = s; }
  if (SOLVERTIME >= t0 && SOLVERTIME < t0 + thbo_dur) {
    ata = ATA; chamber = 1.0; FiO2_c = 1.0;
  }
}
double PL = 760.0*ata - PH2O;

// =========================================================================
// 1. environment and inspired CO
// =========================================================================
double src = (SOLVERTIME < texp) ? Rsrc : 0.0;
dxdt_Cenv  = src/Vroom*1e6 - ACH/60.0*Cenv;
double Cenv_eff = (ppm_fix >= 0.0 && SOLVERTIME < texp) ? ppm_fix : Cenv;
double PIco = (chamber > 0.5) ? 0.0 : (Cenv_eff + ppm_amb)*1e-6*PL;

// =========================================================================
// ventilation sets PaCO2; inspired CO2 (carbogen) clamps it back up
// =========================================================================
double PaCO2ss = PaCO2b*VA/dmax(VA_c, 500.0) + FiCO2*(760.0*ata - PH2O)*0.55;
dxdt_PaCO2 = (PaCO2ss - PaCO2)/2.0;

// =========================================================================
// 2. Coburn-Forster-Kane mass balance on carboxyhaemoglobin
// =========================================================================
double cap   = 1.34*Hb/100.0;                 // mL CO (or O2) per mL blood
double COHb  = ACO/Vb;                        // mL CO per mL blood
double FCOHb = dmin(dmax(COHb/cap, 0.0), 0.98);
double O2Hb  = dmax(cap - COHb, 1e-9);

// mean pulmonary capillary O2 tension: room air by construction, otherwise the
// alveolar value degraded by one efficiency factor for shunt and mask leak
double PAO2c = FiO2_c*(760.0*ata - PH2O) - PaCO2/0.8;
double PcO2  = (FiO2_c <= 0.25) ? 100.0 : shunt*PAO2c;

double Bres = 1.0/DL + PL/dmax(VA_c, 500.0);  // CO transfer resistance
double Vco  = Vco_endo*(1.0 + kHO1co*HO1);    // HO-1 makes CO: positive feedback
dxdt_ACO    = Vco - (COHb*PcO2/(MHald*O2Hb) - PIco)/Bres;

// =========================================================================
// 3. tissue CO distribution
// -------------------------------------------------------------------------
// NOTE (a defect found by the independent Python cross-check): tissue loading
// is driven by blood CO CONTENT, not by the Haldane free tension
// COHb*PcO2/(M*O2Hb).  That tension is a lung-exchange construct and rises
// ~20-fold inside a hyperbaric chamber at fixed COHb, so using it here pumped
// CO INTO muscle during treatment and the model produced HBO-induced
// rhabdomyolysis (CK 4803 with HBO against 100 without).
// =========================================================================
double Ptgt = kappa*FCOHb;
dxdt_Pbr    = (Ptgt - Pbr  )/tau_br;
dxdt_Pht    = (Ptgt - Pht  )/tau_ht;
dxdt_Pmu    = (Ptgt - Pmu  )/tau_mu;
dxdt_Prest  = (Ptgt - Prest)/tau_rest;

// =========================================================================
// 4. oxygen transport: the two hits of one ligand
// =========================================================================
double PaO2 = FiO2_c*(760.0*ata - PH2O) - PaCO2/0.8 - 10.0;
double SaO2 = sat_o2(PaO2, P50, nHill);
double CaO2 = 1.34*Hb*SaO2*(1.0 - FCOHb) + 0.003*dmax(PaO2, 0.0);   // hit 1
double P50e = P50*(1.0 - alphaCO2*FCOHb);                            // hit 2

// cerebral flow: CO2 reactivity times a delivery-defending autoregulatory term
double fc  = 0.28 + 0.72/(1.0 + exp(-(PaCO2   - 32.0)/6.0));
double fc0 = 0.28 + 0.72/(1.0 + exp(-(PaCO2b - 32.0)/6.0));
double fco2 = fc/fc0;
double DO2n = CBF*CaO2/(CBF0*1.34*Hb*0.97);
double CBFss = dmin(CBF0*fco2*(1.0 + 1.25*dmax(0.0, 1.0 - DO2n)), CBFmax*CBF0);
dxdt_CBF = (CBFss - CBF)/tau_cbf;

// tissue (end-capillary) PO2 by Fick on the LEFT-SHIFTED curve
double CMRO2d = CMRO2*(1.0 - 0.35*hypotherm);
double cn_inh = hillf(CNtis, Kc_cn, 1.0);
double CcOtb  = 1.0 - (1.0 - CcOb)*(1.0 - cn_inh);   // one axis, two ligands
double CvO2   = CaO2 - CMRO2d*100.0/dmax(CBF, 1.0);
double Sv     = (CvO2 - 0.003*40.0)/dmax(1.34*Hb*(1.0 - FCOHb), 1e-6);
Sv = dmin(dmax(Sv, 1e-6), 1.0-1e-6);
double PtO2 = dmax(po2_of_sat(Sv, P50e, nHill), 0.05);

// =========================================================================
// 5. myoglobin
// =========================================================================
double MbHtss = 1.0/(1.0 + PtO2/dmax(Mmb*Pht*1e-3, 1e-9));
double MbMuss = 1.0/(1.0 + 40.0 /dmax(Mmb*Pmu*1e-3, 1e-9));
dxdt_MbHt = (MbHtss - MbHt)/6.0;
dxdt_MbMu = (MbMuss - MbMu)/25.0;

// =========================================================================
// 6. terminal oxidase.  CO binds ONLY the reduced a3 haem, so oxygen is a
//    competitive protector: occupancy drive falls as 1/(1 + PtO2/Ko).
//    This is the whole reason the tissue pool behaves differently from COHb.
// =========================================================================
double drive_b = (Pbr*1e-3)/(Kc*(1.0 + PtO2/Ko));
double drive_h = (Pht*1e-3)/(Kc*(1.0 + 40.0/Ko));
dxdt_CcOb = kon_cco*drive_b*(1.0 - CcOb) - koff_cco*(1.0 + hbo_koff*chamber)*CcOb;
dxdt_CcOh = kon_cco*drive_h*(1.0 - CcOh) - koff_cco*(1.0 + hbo_koff*chamber)*CcOh;

// =========================================================================
// 7. energetics.  Km_o2 is a TISSUE-level constant (the far edge of the Krogh
//    diffusion field, ~12 mmHg), not the mitochondrial Km (~1 mmHg): cells die
//    at the edge of the field long before the enzyme runs out of substrate.
//    Supplies are normalised so the healthy state is EXACTLY 1.
// =========================================================================
double fo2_ref = PtO2ref/(PtO2ref + Km_o2);
double fo2     = PtO2/(PtO2 + Km_o2);
double supply_b = dmin((1.0 - CcOtb)/(1.0 - CcO_base)*dmin(1.0, fo2/fo2_ref), 1.0);
dxdt_ATPb = (supply_b - ATPb)/tau_atp;

// the watershed is not handicapped at baseline; it is STEEPER, because it sits
// lower on its own supply curve, so the same fractional fall in PtO2 costs more
double fo2_gp     = (f_wshed*PtO2)   /(f_wshed*PtO2    + Km_o2);
double fo2_gp_ref = (f_wshed*PtO2ref)/(f_wshed*PtO2ref + Km_o2);
double supply_gp  = dmin((1.0 - CcOtb)/(1.0 - CcO_base)*dmin(1.0, fo2_gp/fo2_gp_ref), 1.0);
dxdt_ATPgp = (supply_gp - ATPgp)/tau_atp;

double CcOth = 1.0 - (1.0 - CcOh)*(1.0 - cn_inh);
double supply_h = dmin((1.0 - CcOth)/(1.0 - CcO_base)*(1.0 - 0.55*MbHt), 1.0);
dxdt_ATPh = (supply_h - ATPh)/tau_atp;

dxdt_Lac = kLac*(dmax(0.0, ATPthr - supply_b) + dmax(0.0, ATPthr - supply_h))
           - kLacEl*(Lac - Lac0);

// =========================================================================
// 8. reoxygenation injury.  The ROS burst requires BOTH a prior energy failure
//    (to convert the dehydrogenase) AND restored oxygen (as co-substrate), so
//    the burst is generated by the treatment, not by the poison.
// =========================================================================
double Edef = dmax(0.0, ATPthr - ATPb);        // exactly 0 in health
dxdt_XO = kXOon*Edef*(1.0 - XO) - kXOoff*XO;

double o2avail = PtO2/(PtO2 + 8.0) + 0.55*chamber;
double I_oxy   = 1.0 - Oxy/(Oxy + IC50_oxy);
dxdt_ROS = kROS*XO*o2avail*I_oxy*(2.0 - GSH) - kROSel*ROS;
dxdt_NOx = kNO*(dmax(0.0, FCOHb - 0.01) + 0.5*dmax(0.0, CcOb - CcO_base)) - kNOel*NOx;

// =========================================================================
// 9. innate immune amplification.  HBO blocks beta-2 integrin adhesion
//    DIRECTLY; that action is neither CO clearance nor oxygen delivery, and it
//    persists for hours after the patient leaves the chamber.
// =========================================================================
double hboadh = (since < 1e8) ? hbo_adh*exp(-since/tau_hbo_adh) : 0.0;
dxdt_Neut = kadh*(ROS + 0.6*NOx)*(1.0 - dmin(hboadh, 0.95)) - kdeadh*Neut;
dxdt_MPO  = kMPO*Neut - kMPOel*MPO;
dxdt_LPO  = kLPO*(ROS + wMPO*MPO)*(2.0 - GSH) - kLPOel*LPO;

// =========================================================================
// 10. the DNS switch.  Peroxidation charge-modifies myelin basic protein into
//     an antigen; the clone expands only above a tolerance threshold; and
//     demyelination liberates further antigen, so above the threshold the loop
//     LATCHES.  The weeks-long latency and the all-or-none character are
//     consequences of this structure, not assumptions:
//        clone grows iff Tprol*H(MBPad) > Tdeath
//        =>  MBPad_crit = theta*((Tdeath/Tprol)/(1 - Tdeath/Tprol))^(1/nT)
//     Adduct clearance is on the MYELIN TURNOVER timescale (weeks).  With a
//     days-long clearance the primed clone cannot outrun antigen loss and DNS
//     becomes structurally unreachable at any dose -- a defect the Python
//     cross-check exposed.
// =========================================================================
dxdt_MBPad = (kAd*LPO + kSpread*Demy)*(1.0 - MBPad) - kAdClr*MBPad;
dxdt_Tcell = Tprol*Tcell*hillf(MBPad, theta, nT)*(1.0 - Tcell/Tmax)
             - Tdeath*(Tcell - T0);
dxdt_Micro = kMicro*(dmax(0.0, Tcell - T0) + 0.4*LPO) - kMicroEl*Micro;
dxdt_Demy  = (kDemy*dmax(0.0, Tcell - T0) + kDemyDir*Micro)*(1.0 - Demy)
             - kRepair*Demy;

// =========================================================================
// 11. adaptive protective responses
// =========================================================================
dxdt_HO1 = kHO1*(ROS + 4.0*dmax(0.0, FCOHb - 0.01)) - kHO1el*HO1;
dxdt_HIF = kHIF*Edef - kHIFel*HIF;

// =========================================================================
// 12. organ injury
// =========================================================================
dxdt_Necb  = kNec*dmax(0.0, ATPcrit - ATPgp)*(1.0 - Necb);
dxdt_Edema = kEdema*(Necb + 0.35*dmax(0.0, ATPcrit - ATPgp)) - kEdemaEl*Edema;
dxdt_ICP   = (ICP0 + kICP*Edema - ICP)/8.0;

dxdt_TnI = kTrop*dmax(0.0, ATPcrit_h - ATPh) - kTropEl*TnI;
double EFss = EF0*(1.0 - kEF*dmax(0.0, ATPcrit_h - ATPh));
dxdt_EF  = (EFss - EF)/25.0 + kEFrec*(EF0 - EF);

dxdt_CK = kCK*dmax(0.0, MbMu - MbCK) - kCKel*(CK - CK0);
dxdt_Cr = kCr*dmax(0.0, CK - 5000.0)/1000.0 - kCrEl*(Cr - Cr0);

// =========================================================================
// 13. outcome
// =========================================================================
double Cogss = 1.0 - wDemy*Demy - wNec*Necb;
dxdt_Cog  = (Cogss - Cog)/tau_cog;
dxdt_HazD = kHaz*Demy + 0.5*kHaz*Necb;

// =========================================================================
// 14. cyanide (fire smoke) and its antidote
// =========================================================================
double cn_in = (SOLVERTIME < CN_tau) ? CN_rate : 0.0;
dxdt_CNbl = cn_in - (CL_cn/V_cn)*CNbl - k_rhod*CNbl - k_ohcbl*OHCbl*CNbl
            - k_cnt*(CNbl - CNtis/lam_cn);
dxdt_CNtis = k_cnt*(CNbl - CNtis/lam_cn)*(V_cn/40.0) - 0.010*CNtis;
dxdt_OHCbl = -(CL_ohc/V_ohc)*OHCbl - k_ohcbl*OHCbl*CNbl;

// =========================================================================
// 15. adjunct drug PK and glutathione
// =========================================================================
dxdt_NACc = -(CL_nac/V1_nac)*NACc - (Q_nac/V1_nac)*NACc + (Q_nac/V2_nac)*NACp;
dxdt_NACp =  (Q_nac/V2_nac)*NACc - (Q_nac/V2_nac)*NACp;
dxdt_GSH  = kGSH*(1.0 + 1.6*NACc/(NACc + EC50_nac))*(1.0 - GSH) - kGSHox*ROS*GSH;
dxdt_Oxy  = -(CL_oxy/V_oxy)*Oxy;

// =========================================================================
// 16. fetus.  Fetal haemoglobin binds CO more avidly and the fetal compartment
//     equilibrates slowly, so the fetus lags, peaks higher and clears far more
//     slowly than the mother.  The bedside rule "continue oxygen for about five
//     times the time needed to normalise the mother" is a CONSEQUENCE of this
//     time constant rather than an independent instruction.
// =========================================================================
dxdt_FetCO = (fet_ratio*FCOHb - FetCO)/fet_tau;

$TABLE
double capT   = 1.34*Hb/100.0;
double FCOHbT = ACO/(capT*Vb);
capture COHb_pct  = 100.0*FCOHbT;
capture FetCO_pct = 100.0*FetCO;

// oxygen transport readouts
int    inO2T   = (TIME >= to2_start && TIME < to2_stop) ? 1 : 0;
double FiO2_t  = inO2T ? FiO2_trt : FiO2;
double ataT    = 1.0;
for (int i = 0; i < (int) nhbo; ++i) {
  double t0 = thbo_start + i*hbo_gap;
  if (TIME >= t0 && TIME < t0 + thbo_dur) { ataT = ATA; FiO2_t = 1.0; }
}
double PaO2T = FiO2_t*(760.0*ataT - PH2O) - PaCO2/0.8 - 10.0;
double SaO2T = sat_o2(PaO2T, P50, nHill);
capture PaO2_mmHg  = PaO2T;
capture SaO2_pct   = 100.0*SaO2T*(1.0 - FCOHbT);
capture CaO2       = 1.34*Hb*SaO2T*(1.0 - FCOHbT) + 0.003*dmax(PaO2T, 0.0);
capture O2_diss    = 0.003*dmax(PaO2T, 0.0);
capture P50_eff    = P50*(1.0 - alphaCO2*FCOHbT);

// tissue PO2, recomputed for output
double CvO2T = CaO2 - CMRO2*(1.0 - 0.35*hypotherm)*100.0/dmax(CBF, 1.0);
double SvT   = dmin(dmax((CvO2T - 0.12)/dmax(1.34*Hb*(1.0-FCOHbT), 1e-6), 1e-6), 1.0-1e-6);
capture PtO2_mmHg = dmax(po2_of_sat(SvT, P50_eff, nHill), 0.05);

// ---- MEASUREMENT MODEL: what a two-wavelength pulse oximeter displays ----
// COHb absorbs at 660 nm almost exactly as O2Hb does and is invisible at 940,
// so the ratio-of-ratios inversion attributes the CO-occupied fraction to
// oxygenated haemoglobin.  SpO2 therefore reports (O2Hb + COHb).
double e660_O2 = 0.081, e940_O2 = 0.290;
double e660_HH = 0.845, e940_HH = 0.170;
double e660_CO = 0.083, e940_CO = 0.010;
double fO2 = (1.0 - FCOHbT)*SaO2T;
double fHH = (1.0 - FCOHbT)*(1.0 - SaO2T);
double A660 = fO2*e660_O2 + fHH*e660_HH + FCOHbT*e660_CO;
double A940 = fO2*e940_O2 + fHH*e940_HH + FCOHbT*e940_CO;
double Rr = A660/A940;
double sp = (e660_HH - Rr*e940_HH)/((e660_HH - e660_O2) - Rr*(e940_HH - e940_O2));
sp = dmin(dmax(sp, 0.0), 1.0);
capture SpO2_displayed = 100.0*sp;
capture sat_gap        = 100.0*(sp - fO2);

// composite CcO inhibition (CO and cyanide on one axis) and clinical indices
capture CcO_total_pct = 100.0*(1.0 - (1.0-CcOb)*(1.0 - hillf(CNtis, Kc_cn, 1.0)));
capture DNS_flag      = (Demy > 0.05) ? 1.0 : 0.0;
capture Cog_pct       = 100.0*Cog;

## =============================================================================
##  SCENARIO LIBRARY  (14 scenarios; see co_reference_output.txt for the numbers
##  produced by the independent Python implementation of this same system)
## =============================================================================
##
##  The exposure scale was calibrated by root-finding on peak COHb:
##      mild      538 ppm x 60 min -> COHb 10%
##      moderate 1409 ppm x 60 min -> COHb 25%
##      severe   2298 ppm x 60 min -> COHb 40%
##      critical 4198 ppm x 45 min -> COHb 55%
##  For reference: OSHA 8-h PEL 50 ppm, NIOSH IDLH 1200 ppm, a petrol engine in
##  a closed garage exceeds 30 000 ppm.
##
## ---- 01  severe exposure, no treatment --------------------------------------
##   param(ppm_fix = 2298, texp = 60)
##
## ---- 02  severe, nasal cannula 6 L/min --------------------------------------
##   param(ppm_fix = 2298, texp = 60, to2_start = 90, to2_stop = 450, FiO2_trt = 0.44)
##
## ---- 03  severe, tight non-rebreather mask for 6 h --------------------------
##   param(ppm_fix = 2298, texp = 60, to2_start = 90, to2_stop = 450)
##
## ---- 04  severe, mask for 24 h ----------------------------------------------
##   param(ppm_fix = 2298, texp = 60, to2_start = 90, to2_stop = 1530)
##   (tests whether prolonging normobaric oxygen past COHb clearance adds
##    anything -- in the model it does not, except in pregnancy)
##
## ---- 05  severe, intubated with deliberate hyperventilation -----------------
##   param(ppm_fix = 2298, texp = 60, to2_start = 90, to2_stop = 450, VAtrt = 8000)
##
## ---- 06  severe, carbogen: hyperventilation with PaCO2 clamped --------------
##   param(ppm_fix = 2298, texp = 60, to2_start = 90, to2_stop = 450,
##         VAtrt = 8000, FiCO2 = 0.05)
##
## ---- 07  severe, ONE hyperbaric session at 2 h, 3.0 ATA ---------------------
##   param(ppm_fix = 2298, texp = 60, to2_start = 90, to2_stop = 450,
##         thbo_start = 120, ATA = 3.0, nhbo = 1)
##
## ---- 08  severe, THREE sessions in 24 h starting at 2 h (Weaver protocol) ---
##   param(ppm_fix = 2298, texp = 60, to2_start = 90, to2_stop = 2970,
##         thbo_start = 120, ATA = 3.0, nhbo = 3, hbo_gap = 480)
##
## ---- 09  severe, the same three sessions delayed to 20 h --------------------
##   param(ppm_fix = 2298, texp = 60, to2_start = 90, to2_stop = 2970,
##         thbo_start = 1200, ATA = 3.0, nhbo = 3, hbo_gap = 480)
##   (arms 08 and 09 are the Weaver / Scheinkestel contrast)
##
## ---- 10  severe, three sessions at 2.0 ATA ----------------------------------
##   param(ppm_fix = 2298, texp = 60, to2_start = 90, to2_stop = 2970,
##         thbo_start = 120, ATA = 2.0, nhbo = 3)
##
## ---- 11  critical exposure, mask only --------------------------------------
##   param(ppm_fix = 4198, texp = 45, to2_start = 60, to2_stop = 1500)
##
## ---- 12  critical exposure, early hyperbaric oxygen ------------------------
##   param(ppm_fix = 4198, texp = 45, to2_start = 60, to2_stop = 2940,
##         thbo_start = 90, ATA = 3.0, nhbo = 3)
##
## ---- 13  severe exposure on a background of anaemia (Hb 9 g/dL) ------------
##   param(ppm_fix = 2298, texp = 60, to2_start = 90, to2_stop = 450, Hb = 9)
##
## ---- 14  severe exposure with N-acetylcysteine 150 mg/kg -------------------
##   mod %>% param(ppm_fix = 2298, texp = 60, to2_start = 90, to2_stop = 450) %>%
##       ev(time = 90, amt = 150*70, cmt = "NACc") %>% mrgsim(end = 45*1440)
##
## ---- 15  fire smoke: CO together with cyanide ------------------------------
##   param(ppm_fix = 1600, texp = 60, CN_rate = 4.0, CN_tau = 25,
##         to2_start = 90, to2_stop = 1530)
##   plus hydroxocobalamin:  ev(time = 95, amt = 205, cmt = "OHCbl")
##   (5 g in ~18 L is about 205 uM; the point of the scenario is that oxygen
##    cannot displace cyanide and cobalamin cannot displace CO, so a lactate
##    that will not fall on 100% oxygen is the only bedside discriminator)
##
## ---- 16  pregnancy ---------------------------------------------------------
##   Any of the above; read FetCO_pct rather than COHb_pct.  Maternal COHb
##   reaches 5% at 204 min on a mask while the fetal compartment does not until
##   755 min, a ratio of 3.7 -- which is where the "five times" rule comes from.
##
## ---- 17  hyperbaric half-life floor (a structural demonstration) -----------
##   Sweep ATA from 1.5 to 20 with thbo_start = 0, nhbo = 1, thbo_dur = 1e6 and
##   watch the half-life stop improving: as pressure rises, P_c,O2 and P_L/VA
##   rise together, P_c,O2/B -> VA, and the half-life approaches
##   ln2*M*[O2Hb]*Vb/VA = 31.3 min.  The observed ~20 min at 2.5-3 ATA is below
##   that floor, so CFK is falsified at hyperbaric pressure.  This is reported
##   as a miss, not fitted away.
## =============================================================================
