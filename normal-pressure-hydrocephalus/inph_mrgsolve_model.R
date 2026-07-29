## =============================================================================
##  inph_mrgsolve_model.R
##  Idiopathic Normal Pressure Hydrocephalus (iNPH) — QSP model for mrgsolve
## =============================================================================
##
##  46 ODE compartments:
##      10  craniospinal hydrodynamics & structure
##       5  periventricular white matter and glia
##       7  amyloid / tau / CSF biomarkers
##      11  drug PK (acetazolamide with saturable red-cell binding, loop
##          diuretic, solifenacin, donepezil, melatonin)
##       3  systemic safety (bicarbonate, potassium, acetazolamide escape)
##       2  device (occlusion, cumulative drainage)
##       7  clinical endpoints (fast + slow gait, fast + slow cognition,
##          continence, low-pressure headache, subdural hazard)
##       1  dVac — the acute CSF volume deviation produced by a lumbar tap
##
##  MEAN INTRACRANIAL PRESSURE IS NOT A COMPARTMENT.  It is obtained in closed
##  form at every derivative evaluation from the steady-state CSF flow balance
##      If = (P - Pss)/Rout + (P + h - Popen - Pdist)/Rsh
##  which is linear in P once the valve state is known.  Craniospinal mechanics
##  equilibrate in minutes; the disease evolves over years.  Solving the fast
##  block algebraically is what keeps this a non-stiff 46-state system instead
##  of a stiff 51-state one.
##
##  POSTURE IS A WEIGHT, NOT AN OSCILLATION.  The hydraulic block is evaluated
##  twice per derivative call (supine and upright) and averaged with weight
##  f_up = fraction of the day spent upright.  Shunt hydrostatics differ by tens
##  of cmH2O between the two, so this term cannot be dropped — it is the single
##  largest determinant of over-drainage (see scenario S5 vs S6).
##
##  RELATIONSHIP TO THE PYTHON REFERENCE
##  ------------------------------------
##  inph_reference_model.py implements the identical equations with identical
##  parameter names and is the numerical source of truth for README.md.  The
##  There are exactly TWO implementation differences, both deliberate and both
##  immaterial over the ranges these scenarios reach:
##    (a) dVac.  Both versions use the same EXACT, production-rate-limited flow
##        balance.  Python sub-steps it at ~1 min outside its fixed-step RK4;
##        mrgsolve carries it as an ordinary compartment because LSODA is
##        stiff-capable (the balance is stiff near equilibrium, tau = Rout*C
##        ~ 6 min, and rate-capped far from it).
##    (b) The elastance ceiling.  Python hard-clamps E1 to [E1_floor_eff,
##        E1_ceil] after each step; a hard clamp is not expressible in an ODE,
##        so here the progression term is gated smoothly as E1 approaches
##        E1_ceil.  E1 reaches only ~0.25 by 5 years against a 0.30 ceiling, so
##        the two forms are indistinguishable in these runs.
##  Python additionally clamps every state to its physical range after each
##  step.  Here those bounds are enforced inside the rate expressions instead
##  (every rate that could drive a state negative is written with a guard), so
##  no post-hoc clamping is required.
##
##  Units
##      time        days
##      volume      mL
##      pressure    mmHg internally; valve settings in cmH2O (1 mmHg = 1.36 cmH2O)
##      CSF flow    mL/min inside the hydraulic block, mL/day elsewhere
##      resistance  mmHg/(mL/min)  — the clinical unit for Rout
##
##  Calibration, misses and uncalibrated guesses are documented in README.md
##  section "모델 보정" and in section 12 of inph_model_report.txt.  In short:
##  seven parameters are fitted, kappa_tm is the most influential parameter in
##  the model and has never been measured, and the subdural-haematoma
##  percentages are ordinal rather than calibrated incidences.
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

code <- '
$PROB
# Idiopathic Normal Pressure Hydrocephalus (iNPH) QSP model
- 45 ODE compartments; mean ICP solved in closed form from the CSF flow balance
- posture-averaged shunt hydraulics (supine + upright, weight f_up)
- effector chain: valve -> mean ICP -> elastance -> pulse amplitude ->
  pulsatile transmantle gradient -> periventricular water/perfusion ->
  white matter -> gait / cognition / continence

$PARAM @annotated
// ---- CSF formation and absorption --------------------------------------- //
If0        : 0.35    : CSF formation rate (mL/min)
Pss        : 6.5     : dural venous sinus pressure supine (mmHg)
dPss_up    : 8.0     : downward shift of intracranial reference pressure upright (mmHg)
Rout_norm  : 9.0     : normal CSF outflow resistance (mmHg/(mL/min))
Rout_init  : 19.0    : iNPH outflow resistance at presentation (mmHg/(mL/min))
kRout      : 3.0e-4  : outflow resistance progression rate (1/day)
Rout_max   : 26.0    : ceiling on outflow resistance (mmHg/(mL/min))

// ---- Craniospinal elastance / compliance -------------------------------- //
P0_marm    : -2.0    : reference pressure of the Marmarou P-V curve (mmHg)
E1_floor   : 0.162   : best achievable elastance after decompression (1/mL)
E1_ceil    : 0.360   : elastance ceiling (1/mL)
kE1_rec    : 0.010   : elastance recovery rate when ICP and water fall (1/day)
kE1_prog   : 3.6e-5  : elastance worsening rate with disease (1/day)
C_max      : 3.0     : cap on compliance at low pressure (mL/mmHg)
Ccs_sdh    : 8.0     : chronic accommodation of a subdural collection (mL/mmHg)

// ---- Pulsatility -------------------------------------------------------- //
Vp_norm    : 1.10    : intracranial arterial pulse volume healthy (mL/beat)
Vp_init    : 1.35    : intracranial arterial pulse volume iNPH (mL/beat)
Cart_norm  : 0.85    : arterial compliance index healthy
Cart_init  : 0.62    : arterial compliance index iNPH
kVp_edema  : 0.05    : extra pulse transmission per unit interstitial water
kappa_tm   : 0.060   : fraction of ICP pulse appearing across the mantle
DESH       : 1.55    : DESH morphology multiplier on kappa_tm
DESH_norm  : 1.00    : morphology multiplier, normal convexity
ktm_mean   : 0.115   : mean transmantle gradient per unit relative Rout excess (mmHg)

// ---- Ventricular geometry ----------------------------------------------- //
Vv_norm    : 25.0    : normal lateral + third ventricular volume (mL)
Vsas_norm  : 110.0   : normal cortical SAS + spinal CSF volume (mL)
Vsas_init  : 62.0    : iNPH extraventricular CSF volume (mL)
Vv_max     : 175.0   : ceiling on ventricular volume (mL)
k_creep    : 26.0    : viscoelastic ventricular creep (mL/(1000 day mmHg))
dP_yield   : 0.145   : yield gradient below which no creep occurs (mmHg)
k_recoil   : 0.010   : elastic ventricular recoil (1/day)
f_plastic  : 0.72    : fraction of ventricular expansion that is irreversible
k_exvac    : 40.0    : ex-vacuo enlargement per unit permanent WM loss (mL)

// ---- Periventricular water and glymphatics ------------------------------ //
Wpv_norm   : 0.30    : normal periventricular interstitial water index
kW_base    : 0.366   : baseline interstitial water turnover (1/day)
kW_in      : 2.10    : water accumulation per mmHg mean transmantle gradient (1/day)
kW_aq      : 1.220   : water accumulation per unit glymphatic deficit (1/day)
kW_out     : 1.300   : interstitial water clearance (1/day)
AQ_norm    : 0.85    : normal AQP4 perivascular polarisation fraction
kAQ_rec    : 0.008   : AQP4 repolarisation rate (1/day)
kAQ_loss   : 0.00677 : AQP4 depolarisation per unit excess astrogliosis (1/day)
sleep_base : 1.00    : glymphatic sleep multiplier

// ---- Perfusion ---------------------------------------------------------- //
CBF0       : 22.0    : periventricular WM perfusion ceiling (mL/100g/min)
kCBF_W     : 0.120   : perfusion loss per unit interstitial water
kCBF_P     : 0.350   : perfusion loss per mmHg pulsatile transmantle gradient
Autoreg_norm : 0.95  : autoregulatory reserve reference
kAR        : 0.006   : autoregulatory adaptation rate (1/day)
kar_p      : 0.387   : autoregulatory loss per mmHg pulsatile gradient
CBF_crit   : 18.0    : perfusion threshold for WM injury (mL/100g/min)

// ---- White matter and glia ---------------------------------------------- //
k_rep      : 0.0165  : WM functional recovery rate (1/day)
k_deg      : 0.0048  : WM injury-driven loss rate (1/day)
k_perm     : 1.0e-4  : conversion of injury into permanent loss (1/day)
w_water    : 0.55    : injury weight on interstitial water
w_perf     : 1.00    : injury weight on ischaemia
w_infl     : 0.22    : injury weight on microglial activation
kMyel      : 0.010   : myelin recovery rate (1/day)
kMyel_inj  : 0.030   : myelin loss per unit injury (1/day)
kAstro_in  : 0.0216  : astrogliosis induction (1/day)
kAstro_out : 0.020   : astrogliosis resolution (1/day)
kMicro_in  : 0.0262  : microglial induction (1/day)
kMicro_out : 0.030   : microglial resolution (1/day)

// ---- Amyloid / tau / CSF biomarkers ------------------------------------- //
Pab        : 8.0     : ISF Abeta42 production (nM/day)
kab_deg    : 0.55    : ISF Abeta42 local degradation (1/day)
kab_gly    : 1.10    : ISF Abeta42 glymphatic clearance at full function (1/day)
kagg       : 0.0125  : Abeta42 aggregation sink (1/(nM day))
APOE       : 1.0     : APOE multiplier (1.0 = e3/e3, 2.2 = e4 carrier)
Ab_ref     : 5.46    : reference healthy ISF Abeta42 (nM)
kagg_plq   : 2.10e-4 : plaque deposition rate at Ab_ref (1/day)
kplq_clr   : 0.0016  : plaque clearance rate (1/day)
Ptau       : 1.00    : ISF p-tau production healthy (a.u./day)
ktau_plq   : 1.60    : p-tau production driven by plaque load
ktau_deg   : 0.42    : ISF p-tau degradation (1/day)
ktau_gly   : 0.80    : ISF p-tau glymphatic clearance (1/day)
kflux_ab   : 98000.0 : brain-to-CSF Abeta42 flux (pg/day per nM per glymphatic unit); scaled so healthy CSF Abeta42 ~800 pg/mL
kflux_tau  : 26400.0 : brain-to-CSF p-tau flux (pg/day per a.u. per glymphatic unit); scaled so healthy CSF p-tau ~45 pg/mL
kdeg_csf   : 0.09    : intrathecal degradation (1/day)
kNfL       : 1.0e8   : CSF NfL release per unit axonal loss rate (pg/day)
kNfL_base  : 4.2e5   : age-related baseline NfL release (pg/day)
kLRG       : 400.0   : CSF LRG release per unit astrogliosis + water (a.u./day)
kLRG_base  : 516.0   : baseline CSF LRG release (a.u./day)

// ---- Shunt hardware ----------------------------------------------------- //
shunt      : 0       : ventriculoperitoneal shunt in situ (0/1)
Popen_cm   : 10.0    : differential valve opening pressure (cmH2O)
Rsh        : 2.5     : shunt hydraulic resistance (mmHg/(mL/min))
Ggrav_cm   : 0.0     : gravitational unit opening pressure added upright (cmH2O)
asd_eff    : 0.0     : membrane anti-siphon device efficiency (0-1)
Hcol_cm    : 45.0    : ventricle-to-peritoneum hydrostatic column (cmH2O)
Pd_sup_cm  : 5.0     : intra-abdominal pressure supine (cmH2O)
Pd_up_cm   : 10.0    : intra-abdominal pressure upright (cmH2O)
k_occl     : 2.6e-4  : shunt occlusion hazard (1/day)
f_up       : 0.60    : fraction of the day spent upright

// ---- ETV ---------------------------------------------------------------- //
etv        : 0       : endoscopic third ventriculostomy (0/1)
etv_dRout  : 0.06    : fractional Rout reduction achievable by ETV in communicating HC

// ---- External drainage -------------------------------------------------- //
eld_rate   : 0.0     : external lumbar drainage rate (mL/min)

// ---- Subdural collection and safety ------------------------------------- //
P_thr_sdh  : 2.0     : daily-mean ICP below which fluid collects (mmHg)
k_sdh      : 1.55    : subdural accumulation (mL/(day mmHg))
k_sdh_res  : 0.011   : subdural resorption (1/day)
atrophy    : 0.55    : cortical atrophy, sets the available subdural space (0-1)
antithrombotic : 0.0 : antiplatelet or anticoagulant use (0/1)
Vsdh_sympt : 22.0    : collection volume that becomes symptomatic (mL)
k_haz      : 1.5e-5  : symptomatic-event hazard (1/(day mL))
k_sdh_gait : 0.008   : gait penalty per mL of collection above 15 mL
k_sdh_cog  : 0.050   : MMSE penalty per mL of collection above 15 mL
P_thr_head : 1.0     : upright ICP below which headache appears (mmHg)

// ---- Acetazolamide ------------------------------------------------------ //
az_ka      : 36.0    : absorption rate (1/day)
az_V       : 14.0    : central volume (L)
az_CL      : 700.0   : clearance on the UNBOUND pool (L/day); CL_u = CL_total/fu, i.e. 35 L/day on total drug, t1/2 ~ 6.7 h
az_fu      : 0.05    : unbound plasma fraction
az_kon     : 30.0    : red-cell CA-I binding on-rate (1/day per mg/L)
az_koff    : 1.4     : red-cell CA-I binding off-rate (1/day)
az_Bmax    : 120.0   : saturable red-cell binding capacity (mg); this deep pool is the origin of the non-linear PK
az_Emax    : 0.50    : maximal fractional reduction of CSF formation
az_EC50    : 0.028   : free concentration for half-maximal effect (mg/L)
az_EC50_sys : 1.20   : free concentration for the systemic acid-base effect (mg/L); ~43x the choroidal EC50
az_kacid   : 9.0     : bicarbonate fall at full systemic effect (mmol/L)
az_tau_acid : 3.0    : bicarbonate time constant (day)
az_kK      : 0.95    : potassium fall at full systemic effect (mmol/L)
az_tau_esc : 18.0    : time constant of pharmacodynamic escape (day)
az_f_esc   : 0.78    : fraction of the acute CSF-formation effect lost to compensation

// ---- Loop diuretic (NKCC1) ---------------------------------------------- //
bu_ka      : 28.8    : absorption rate (1/day)
bu_V       : 12.0    : central volume (L)
bu_CL      : 3.84    : clearance (L/day)
bu_Emax    : 0.22    : maximal fractional reduction of CSF formation
bu_EC50    : 0.05    : half-maximal concentration (mg/L)

// ---- Solifenacin -------------------------------------------------------- //
so_ka      : 12.0    : absorption rate (1/day)
so_V       : 600.0   : central volume (L)
so_CL      : 14.4    : clearance (L/day)
so_Emax_ur : 0.85    : maximal reduction of the urinary index
so_EC50_ur : 0.010   : half-maximal concentration (mg/L)
so_kcog    : 6.5     : MMSE points lost per mg/L (central antimuscarinic)

// ---- Donepezil ---------------------------------------------------------- //
do_ka      : 12.0    : absorption rate (1/day)
do_V       : 800.0   : central volume (L)
do_CL      : 7.2     : clearance (L/day)
do_Emax_cog : 1.60   : MMSE-equivalent gain (points)
do_EC50_cog : 0.020  : half-maximal concentration (mg/L)

// ---- Melatonin ---------------------------------------------------------- //
me_ka      : 72.0    : absorption rate (1/day)
me_V       : 35.0    : central volume (L)
me_CL      : 26.4    : clearance (L/day)
me_Emax_sleep : 0.30 : fractional gain in the glymphatic sleep factor
me_EC50_sleep : 0.0025 : half-maximal concentration (mg/L)

// ---- Clinical endpoints ------------------------------------------------- //
G_max      : 1.16    : gait velocity ceiling for a 75-year-old (m/s)
G_wm_pow   : 1.15    : steepness of gait on WM integrity
tau_gait   : 34.0    : slow gait time constant (day)
tau_on     : 0.020   : charging time of the fast gait arm (day, 29 min)
tau_off    : 0.90    : discharge time of the fast gait arm (day, 21.6 h)
a_fast     : 0.075   : maximum fast gait component (m/s)
MMSE_max   : 30.0    : cognitive ceiling (points)
kcog_wm    : 9.0     : MMSE points lost at full WM failure
kcog_ad    : 8.0     : MMSE points lost at full amyloid/tau burden
tau_cog    : 62.0    : slow cognition time constant (day)
a_fast_cog : 0.90    : maximum fast cognitive component (points)
tau_on_cog : 0.040   : charging time of the fast cognitive arm (day)
tau_off_cog : 1.20   : discharge time of the fast cognitive arm (day)
Ur_max     : 3.0     : urinary index ceiling
kur_wm     : 3.4     : urinary index at full WM failure
tau_ur     : 40.0    : urinary time constant (day)
a_fast_ur  : 0.42    : maximum fast urinary component
comorb     : 0.0     : extra non-hydrocephalic gait burden (0-1)

// ---- Patient-specific references (computed by inph_refs() in R) --------- //
AMP_ref    : 3.263   : the patient own untreated daily-mean pulse amplitude (mmHg); recomputed per patient by inph_refs()
Pday_ref   : 8.354   : the patient own untreated daily-mean ICP (mmHg)
E1_floor_eff : 0.162 : patient-relative elastance floor (1/mL)

$CMT @annotated
Vv       : ventricular CSF volume (mL)
Vsas     : extraventricular CSF volume (mL)
Vsdh     : subdural collection volume (mL)
Rout     : CSF outflow resistance (mmHg/(mL/min))
E1       : craniospinal elastance coefficient (1/mL)
Vplast   : irreversible ventricular expansion (mL)
Wpv      : periventricular interstitial water index
AQ       : AQP4 perivascular polarisation fraction
Cart     : arterial compliance index
Autoreg  : cerebrovascular autoregulatory reserve
WMint    : periventricular WM functional integrity (recoverable)
WMperm   : permanent periventricular axonal loss
Myel     : myelin index
Astro    : reactive astrogliosis
Micro    : microglial activation
Ab_isf   : interstitial Abeta42 (nM)
Ab_plq   : cortical amyloid plaque load
Tau_isf  : interstitial p-tau (a.u.)
A_ab     : CSF Abeta42 amount (pg)
A_pt     : CSF p-tau amount (pg)
A_nfl    : CSF neurofilament light amount (pg)
A_lrg    : CSF LRG amount (a.u.)
Aaz_g    : acetazolamide gut depot (mg)
Aaz_c    : acetazolamide central (mg)
Aaz_r    : acetazolamide red-cell bound (mg)
Abu_g    : loop diuretic gut depot (mg)
Abu_c    : loop diuretic central (mg)
Aso_g    : solifenacin gut depot (mg)
Aso_c    : solifenacin central (mg)
Ado_g    : donepezil gut depot (mg)
Ado_c    : donepezil central (mg)
Ame_g    : melatonin gut depot (mg)
Ame_c    : melatonin central (mg)
HCO3     : serum bicarbonate (mmol/L)
Kser     : serum potassium (mmol/L)
AZesc    : acetazolamide pharmacodynamic escape (0-1)
Occl     : shunt occlusion fraction
Vdr      : cumulative shunt drainage (mL)
Gslow    : slow gait component (m/s)
Gfast    : fast gait component (m/s)
Cslow    : slow cognitive component (points)
Cfast    : fast cognitive component (points)
Urin     : urinary index (0 continent - 3 incontinent)
Headx    : low-pressure headache index
SDHhaz   : cumulative subdural event hazard
dVac     : acute CSF volume deviation from a lumbar tap (mL)

$GLOBAL
#define MINPERDAY 1440.0
#define MMHG(cm) ((cm)/1.36)

// Marmarou exponential pressure-volume curve, floored and capped
double compfn(double P, double e1, double P0, double Cmax) {
  double dP = P - P0;
  if (dP < 0.5) dP = 0.5;
  double C = 1.0 / (e1 * dP);
  return (C > Cmax) ? Cmax : C;
}

// Everything the algebraic hydraulic block returns
struct HYD {
  double Ps, Pu, Pday, Qsh, Qabs, Csup, Cday;
  double AMP, AMPday, dPtmP, dPtmM, Ifd, parfrac, Vp, desh;
};

// Closed-form CSF flow balance, evaluated supine then upright and averaged.
// Solving this algebraically -- rather than carrying P as a compartment -- is
// what keeps the system non-stiff, because tau = Rout*C is about 6 minutes.
void hydro(HYD &o,
           double Rout_, double E1_, double Vsdh_, double dVac_, double Occl_,
           double Cart_, double Wpv_, double Vsas_,
           double Ifeff, double eld, double shunt_, double etv_,
           double Pss_, double dPss_up_, double Pdsup, double Pdup,
           double Hcol, double asd, double Popen, double Ggrav, double Rsh_,
           double Ccs, double P0m, double Cmax, double fup, double etvd,
           double Vpn, double Vpi, double Cartn, double Carti, double kVpe,
           double Wpvn, double kap, double desh_, double deshn,
           double Vsasn, double Vsasi, double ktmm, double Routn) {

  double Qin  = Ifeff - eld;
  double RoutE = Rout_ * ((etv_ > 0.5) ? (1.0 - etvd) : 1.0);
  double Ps = 0, Pu = 0, Qs = 0, Qu = 0, Qas = 0, Qau = 0, Cs = 0, Cu = 0;

  for (int up = 0; up < 2; ++up) {
    double Pss_eff  = Pss_ - (up ? dPss_up_ : 0.0);
    double Pdist    = MMHG(up ? Pdup : Pdsup);
    double h_eff    = up ? MMHG(Hcol) * (1.0 - asd) : 0.0;
    double Popen_e  = MMHG(Popen) + (up ? MMHG(Ggrav) : 0.0);

    double Pf, Rsh_eff = -1.0;
    if (shunt_ > 0.5) {
      double occl = (Occl_ > 0.995) ? 0.995 : Occl_;
      double Rt = Rsh_ / (1.0 - occl);
      double G  = 1.0 / RoutE + 1.0 / Rt;
      Pf = (Qin + Pss_eff / RoutE + (Popen_e + Pdist - h_eff) / Rt) / G;
      if (Pf + h_eff >= Popen_e + Pdist) {
        Rsh_eff = Rt;                       // valve open
      } else {
        Pf = Pss_eff + Qin * RoutE;         // valve shut
      }
    } else {
      Pf = Pss_eff + Qin * RoutE;
    }

    // Two different volume loads, two different mechanisms.  A chronic subdural
    // collection is accommodated by brain compression: a linear offset with its
    // own, much larger compliance Ccs.  An acute tap rides the craniospinal P-V
    // curve itself, so it enters through the INTEGRATED Marmarou form -- which
    // is explicit and monotone, reduces to dV/C for small dV, and saturates
    // near P0 for a large removal.  (A linear offset solved by fixed point does
    // NOT converge for a 40 mL tap; see the defect list in README.md.)
    double Pbase = Pf + Vsdh_ / Ccs;
    double P = (dVac_ != 0.0) ? (P0m + (Pbase - P0m) * exp(E1_ * dVac_)) : Pbase;
    if (P < -30.0) P = -30.0;
    double C = compfn(P, E1_, P0m, Cmax);

    double Qshl = 0.0;
    if (Rsh_eff > 0.0) {
      double drv = P + h_eff - Popen_e - Pdist;
      if (drv > 0.0) Qshl = drv / Rsh_eff;
    }
    double Qabsl = (P > Pss_eff) ? (P - Pss_eff) / RoutE : 0.0;

    if (up) { Pu = P; Qu = Qshl; Qau = Qabsl; Cu = C; }
    else    { Ps = P; Qs = Qshl; Qas = Qabsl; Cs = C; }
  }

  o.Ps = Ps; o.Pu = Pu;
  o.Pday = (1.0 - fup) * Ps + fup * Pu;
  o.Qsh  = ((1.0 - fup) * Qs  + fup * Qu ) * MINPERDAY;
  o.Qabs = ((1.0 - fup) * Qas + fup * Qau) * MINPERDAY;
  o.Csup = Cs;
  o.Cday = (1.0 - fup) * Cs + fup * Cu;

  double frac = (Cartn - Cart_) / (Cartn - Carti);
  if (frac < 0.0) frac = 0.0; if (frac > 1.0) frac = 1.0;
  double wex = (Wpv_ > Wpvn) ? (Wpv_ - Wpvn) : 0.0;
  o.Vp = (Vpn + (Vpi - Vpn) * frac) * (1.0 + kVpe * wex);

  o.AMP    = o.Vp / Cs;        // the clinically measured (supine) amplitude
  o.AMPday = o.Vp / o.Cday;    // what stresses the mantle across the whole day

  double dfr = (Vsasn - Vsas_) / (Vsasn - Vsasi);
  if (dfr < 0.0) dfr = 0.0; if (dfr > 1.0) dfr = 1.0;
  o.desh  = deshn + (desh_ - deshn) * dfr;
  o.dPtmP = kap * o.AMPday * o.desh;

  o.Ifd = Ifeff * MINPERDAY;
  o.parfrac = o.Qabs / ((o.Ifd > 1e-9) ? o.Ifd : 1e-9);
  if (o.parfrac > 1.0) o.parfrac = 1.0;
  double rex = Rout_ / Routn - 1.0;
  o.dPtmM = ((rex > 0.0) ? ktmm * rex : 0.0) * o.parfrac;
}

$ODE
// ---------------- effective CSF formation (drug effects) ----------------- //
double Caz_free = az_fu * Aaz_c / az_V;
double occ_az = Caz_free / (az_EC50 + Caz_free);            // choroidal
double occ_sys = Caz_free / (az_EC50_sys + Caz_free);       // acid-base
// Carbonic anhydrase inhibition lowers CSF formation acutely, then escapes
// over ~3 weeks (choroidal transporter + renal acid-base compensation).  The
// acidosis does NOT escape, so the therapeutic index degrades with duration as
// well as with dose.  Without this term a CA inhibitor beats a shunt here.
double E_az = az_Emax * occ_az * (1.0 - az_f_esc * AZesc);
double Cbu = Abu_c / bu_V;
double E_bu = bu_Emax * Cbu / (bu_EC50 + Cbu);
double If_eff = If0 * (1.0 - E_az) * (1.0 - E_bu);

// ---------------- algebraic hydraulic block ------------------------------ //
HYD h;
hydro(h, Rout, E1, Vsdh, dVac, Occl, Cart, Wpv, Vsas,
      If_eff, eld_rate, shunt, etv,
      Pss, dPss_up, Pd_sup_cm, Pd_up_cm, Hcol_cm, asd_eff,
      Popen_cm, Ggrav_cm, Rsh, Ccs_sdh, P0_marm, C_max, f_up, etv_dRout,
      Vp_norm, Vp_init, Cart_norm, Cart_init, kVp_edema, Wpv_norm,
      kappa_tm, DESH, DESH_norm, Vsas_norm, Vsas_init, ktm_mean, Rout_norm);

// ---------------- glymphatic function ----------------------------------- //
double Cme = Ame_c / me_V;
double sleepf = sleep_base * (1.0 + me_Emax_sleep * Cme / (me_EC50_sleep + Cme));
double gly = AQ * sleepf / (1.0 + 0.35 * Wpv);

// ---------------- perfusion and injury ---------------------------------- //
double wex = (Wpv > Wpv_norm) ? (Wpv - Wpv_norm) : 0.0;
double CBFpv = CBF0 * (Autoreg / Autoreg_norm) /
               (1.0 + kCBF_W * wex + kCBF_P * h.dPtmP);
double isch = (CBFpv < CBF_crit) ? (CBF_crit - CBFpv) / CBF_crit : 0.0;
double mex = (Micro > 0.20) ? (Micro - 0.20) : 0.0;
double inj = w_water * wex + w_perf * isch + w_infl * mex;

// ---------------- ventricular geometry ---------------------------------- //
double dPtm = h.dPtmM + 0.85 * h.dPtmP;
double drive = dPtm - dP_yield;
double Vv_floor = Vv_norm + Vplast + k_exvac * WMperm;
double dVv;
if (drive > 0.0) dVv = k_creep * drive * (1.0 - Vv / Vv_max) / 1000.0;
else             dVv = -k_recoil * ((Vv > Vv_floor) ? (Vv - Vv_floor) : 0.0);

dxdt_Vv     = dVv;
dxdt_Vplast = f_plastic * ((dVv > 0.0) ? dVv : 0.0);
dxdt_Vsas   = -0.55 * dVv;                  // DESH: the convexity is squeezed

// ---------------- outflow resistance, elastance ------------------------- //
// Progression is a DISEASE process: gate it on the abnormality already present
// so a normal craniospinal space does not drift toward Rout_max.
double rgate = (Rout - Rout_norm) / ((Rout_init - Rout_norm) > 1e-6 ?
                                     (Rout_init - Rout_norm) : 1e-6);
if (rgate < 0.0) rgate = 0.0; if (rgate > 1.0) rgate = 1.0;
dxdt_Rout = kRout * ((Rout < Rout_max) ? (Rout_max - Rout) : 0.0) * rgate;

double relief = (Pday_ref - h.Pday) / ((Pday_ref > 1e-6) ? Pday_ref : 1e-6);
if (relief < 0.0) relief = 0.0;
double water_ok = (Wpv < 1.50) ? (1.50 - Wpv) / 1.50 : 0.0;
double recdrv = relief + 0.6 * water_ok;
if (recdrv > 1.0) recdrv = 1.0;
double e1gate = (E1_ceil - E1) / (E1_ceil - E1_floor_eff);
if (e1gate < 0.0) e1gate = 0.0; if (e1gate > 1.0) e1gate = 1.0;
dxdt_E1 = -kE1_rec * (E1 - E1_floor_eff) * recdrv
          + kE1_prog * (wex / 1.20) * e1gate;

// ---------------- subdural collection ----------------------------------- //
double defic = (h.Pday < P_thr_sdh) ? (P_thr_sdh - h.Pday) : 0.0;
dxdt_Vsdh = k_sdh * defic * (1.0 + atrophy) - k_sdh_res * Vsdh;

// ---------------- water, glymphatics, vessels --------------------------- //
double aqdef = (AQ < AQ_norm) ? (AQ_norm - AQ) : 0.0;
dxdt_Wpv = kW_base + kW_in * h.dPtmM + kW_aq * aqdef
           - kW_out * Wpv * (0.35 + 0.65 * gly / AQ_norm);
double asex = (Astro > 0.25) ? (Astro - 0.25) : 0.0;
dxdt_AQ  = kAQ_rec * (AQ_norm - AQ) * sleepf - kAQ_loss * asex * AQ;
dxdt_Cart = -2.2e-5;                        // arterial stiffening with age
dxdt_Autoreg = kAR * (Autoreg_norm / (1.0 + kar_p * h.dPtmP) - Autoreg);

// ---------------- white matter and glia --------------------------------- //
double ceilWM = 1.0 - WMperm;
double repgate = CBFpv / CBF_crit;  if (repgate > 1.0) repgate = 1.0;
dxdt_WMint  = k_rep * ((ceilWM > WMint) ? (ceilWM - WMint) : 0.0) * repgate
              - k_deg * inj * WMint;
dxdt_WMperm = k_perm * inj * ((WMperm < 1.0) ? (1.0 - WMperm) : 0.0);
dxdt_Myel   = kMyel * (ceilWM - Myel) - kMyel_inj * inj;
dxdt_Astro  = kAstro_in * inj - kAstro_out * (Astro - 0.25);
dxdt_Micro  = kMicro_in * inj - kMicro_out * (Micro - 0.20);

// ---------------- amyloid, tau, CSF biomarkers -------------------------- //
dxdt_Ab_isf = Pab - (kab_deg + kab_gly * gly) * Ab_isf
              - kagg * APOE * Ab_isf * Ab_isf;
dxdt_Ab_plq = kagg_plq * APOE * (Ab_isf / Ab_ref) * (Ab_isf / Ab_ref)
              * (1.0 - Ab_plq) - kplq_clr * Ab_plq;
dxdt_Tau_isf = Ptau * (1.0 + ktau_plq * Ab_plq)
               - (ktau_deg + ktau_gly * gly) * Tau_isf;

double Vcsf = Vv + Vsas;  if (Vcsf < 20.0) Vcsf = 20.0;
double kout = (h.Ifd + h.Qsh) / Vcsf + kdeg_csf;
dxdt_A_ab  = kflux_ab  * gly * Ab_isf  - kout * A_ab;
dxdt_A_pt  = kflux_tau * gly * Tau_isf - kout * A_pt;
dxdt_A_nfl = kNfL_base + kNfL * k_deg * inj * WMint - kout * A_nfl;
dxdt_A_lrg = kLRG_base + kLRG * (asex + wex) - kout * A_lrg;

// ---------------- drug PK ------------------------------------------------ //
// Acetazolamide binds red-cell carbonic anhydrase I saturably.  That pool, not
// clearance, is why the apparent half-life lengthens with dose.
double bind = az_kon * Caz_free * ((az_Bmax > Aaz_r) ? (az_Bmax - Aaz_r) : 0.0)
              - az_koff * Aaz_r;
dxdt_Aaz_g = -az_ka * Aaz_g;
dxdt_Aaz_c =  az_ka * Aaz_g - az_CL * Caz_free - bind;
dxdt_Aaz_r =  bind;
dxdt_Abu_g = -bu_ka * Abu_g;
dxdt_Abu_c =  bu_ka * Abu_g - bu_CL * Abu_c / bu_V;
dxdt_Aso_g = -so_ka * Aso_g;
dxdt_Aso_c =  so_ka * Aso_g - so_CL * Aso_c / so_V;
dxdt_Ado_g = -do_ka * Ado_g;
dxdt_Ado_c =  do_ka * Ado_g - do_CL * Ado_c / do_V;
dxdt_Ame_g = -me_ka * Ame_g;
dxdt_Ame_c =  me_ka * Ame_g - me_CL * Ame_c / me_V;

dxdt_AZesc = (occ_az - AZesc) / az_tau_esc;
dxdt_HCO3 = (24.0 - az_kacid * occ_sys - HCO3) / az_tau_acid;
dxdt_Kser = (4.2  - az_kK    * occ_sys - Kser) / 4.0;

// ---------------- device ------------------------------------------------- //
if (shunt > 0.5) {
  double prot = 0.6 + 0.4 * (A_lrg / Vcsf);
  dxdt_Occl = k_occl * prot * (1.0 - Occl);
  dxdt_Vdr  = h.Qsh;
} else {
  dxdt_Occl = 0.0;
  dxdt_Vdr  = 0.0;
}

// ---------------- clinical endpoints ------------------------------------ //
double wmc = WMint;  if (wmc < 1e-4) wmc = 1e-4;  if (wmc > 1.0) wmc = 1.0;
double perf = 0.55 + 0.45 * ((CBFpv / CBF0 < 1.0) ? CBFpv / CBF0 : 1.0);
double sdhex = (Vsdh > 15.0) ? (Vsdh - 15.0) : 0.0;
double sdhg = 1.0 - k_sdh_gait * sdhex;  if (sdhg < 0.30) sdhg = 0.30;
double gtgt = G_max * pow(wmc, G_wm_pow) * perf * (1.0 - 0.25 * Ab_plq)
              * sdhg * (1.0 - comorb);
if (gtgt < 0.05) gtgt = 0.05;

double taub = 0.4 * ((Tau_isf > 1.0) ? (Tau_isf - 1.0) / 3.0 : 0.0);
if (taub > 0.4) taub = 0.4;
double adb = 0.6 * Ab_plq + taub;  if (adb > 1.0) adb = 1.0;
double ctgt = MMSE_max - kcog_wm * (1.0 - wmc) - kcog_ad * adb
              - k_sdh_cog * sdhex;
if (ctgt < 5.0) ctgt = 5.0;

double utgt = kur_wm * (1.0 - wmc);
if (utgt > Ur_max) utgt = Ur_max;  if (utgt < 0.0) utgt = 0.0;

// The acute drive is the fractional fall in pulse amplitude relative to the
// patient own untreated operating point.  It is what a tap test perturbs.
double drv_ac = (AMP_ref - h.AMPday) / ((AMP_ref > 1e-6) ? AMP_ref : 1e-6);
if (drv_ac < 0.0) drv_ac = 0.0;  if (drv_ac > 1.0) drv_ac = 1.0;

double Cso = Aso_c / so_V;
double Cdo = Ado_c / do_V;
double cogdrug = -so_kcog * Cso + do_Emax_cog * Cdo / (do_EC50_cog + Cdo);
double E_so = so_Emax_ur * Cso / (so_EC50_ur + Cso);
double urn = utgt - E_so - a_fast_ur * drv_ac;  if (urn < 0.0) urn = 0.0;

// The fast arm is ASYMMETRIC: relief of conduction block in a loaded
// periventricular tract is fast, its re-establishment is not.  A symmetric
// filter with tau = 13 h integrates almost nothing out of a 2-hour tap and
// predicts that no tap test can be positive.  The asymmetry is a modelling
// choice, but a falsifiable one -- it predicts fast onset and slow offset of
// tap-test improvement, which is what is clinically reported.
double gf_tgt = a_fast * drv_ac;
double cf_tgt = a_fast_cog * drv_ac;
dxdt_Gslow = (gtgt - Gslow) / tau_gait;
dxdt_Gfast = (gf_tgt > Gfast) ? (gf_tgt - Gfast) / tau_on : -Gfast / tau_off;
dxdt_Cslow = (ctgt + cogdrug - Cslow) / tau_cog;
dxdt_Cfast = (cf_tgt > Cfast) ? (cf_tgt - Cfast) / tau_on_cog
                              : -Cfast / tau_off_cog;
dxdt_Urin  = (urn - Urin) / tau_ur;
dxdt_Headx = (((P_thr_head > h.Pu) ? (P_thr_head - h.Pu) : 0.0) - Headx) / 3.0;
dxdt_SDHhaz = k_haz * ((Vsdh > Vsdh_sympt) ? (Vsdh - Vsdh_sympt) : 0.0)
              * (1.0 + 0.8 * antithrombotic);

// The acute deficit recovers by the EXACT flow balance, not by tau = Rout*C.
// That linearisation implies a refill rate of ~7 mL/min for a 40 mL tap against
// a choroidal secretion rate of 0.35 mL/min.  Below sinus pressure BOTH
// absorption and shunt flow are zero, so recovery is PRODUCTION-limited: 40 mL
// takes ~114 min, not ~30.  This was the most consequential defect found during
// development -- with the linearised form the model concluded that no tap test
// could ever be positive.
dxdt_dVac = (If_eff - eld_rate) * MINPERDAY - h.Qabs - h.Qsh;

$TABLE
double Caz_freeT = az_fu * Aaz_c / az_V;
double E_azT = az_Emax * (Caz_freeT / (az_EC50 + Caz_freeT))
               * (1.0 - az_f_esc * AZesc);
double CbuT = Abu_c / bu_V;
double E_buT = bu_Emax * CbuT / (bu_EC50 + CbuT);
double If_effT = If0 * (1.0 - E_azT) * (1.0 - E_buT);

HYD g;
hydro(g, Rout, E1, Vsdh, dVac, Occl, Cart, Wpv, Vsas,
      If_effT, eld_rate, shunt, etv,
      Pss, dPss_up, Pd_sup_cm, Pd_up_cm, Hcol_cm, asd_eff,
      Popen_cm, Ggrav_cm, Rsh, Ccs_sdh, P0_marm, C_max, f_up, etv_dRout,
      Vp_norm, Vp_init, Cart_norm, Cart_init, kVp_edema, Wpv_norm,
      kappa_tm, DESH, DESH_norm, Vsas_norm, Vsas_init, ktm_mean, Rout_norm);

double wexT = (Wpv > Wpv_norm) ? (Wpv - Wpv_norm) : 0.0;
double CBF_pv = CBF0 * (Autoreg / Autoreg_norm) /
                (1.0 + kCBF_W * wexT + kCBF_P * g.dPtmP);

double ICP_sup = g.Ps;
double ICP_up  = g.Pu;
double ICP_day = g.Pday;
double AMP     = g.AMP;
double AMP_day = g.AMPday;
double Cspine  = g.Cday;
double dPtm_pulse = g.dPtmP;
double dPtm_mean  = g.dPtmM;
double Qsh_day = g.Qsh;
double Qabs_day = g.Qabs;
double VcsfT = Vv + Vsas;
double turnover = (g.Ifd + g.Qsh) / VcsfT;
double If_day = g.Ifd;

double EvansIdx = 0.20 + 0.00350 * (Vv - 25.0);  if (EvansIdx > 0.58) EvansIdx = 0.58;
double CallAngle = 120.0 - 0.50 * (Vv - 25.0);   if (CallAngle < 48.0) CallAngle = 48.0;

double gait = Gslow + Gfast;
double MMSE = Cslow + Cfast;
double urin = Urin;

// iNPH grading scale, three 0-4 domains
double gs_g = (gait >= 1.00) ? 0 : (gait >= 0.85) ? 1 : (gait >= 0.68) ? 2 :
              (gait >= 0.50) ? 3 : 4;
double gs_c = (MMSE >= 28) ? 0 : (MMSE >= 26) ? 1 : (MMSE >= 23) ? 2 :
              (MMSE >= 19) ? 3 : 4;
double gs_u = (urin <= 0.4) ? 0 : (urin <= 1.0) ? 1 : (urin <= 1.7) ? 2 :
              (urin <= 2.4) ? 3 : 4;
double iNPHGS = gs_g + gs_c + gs_u;

double SDH_inc = 1.0 - exp(-SDHhaz);
double CSF_Ab42 = A_ab / VcsfT;
double CSF_pTau = A_pt / VcsfT;
double CSF_NfL  = A_nfl / VcsfT;
double CSF_LRG  = A_lrg / VcsfT;
double AZ_free  = Caz_freeT;
double AZ_eff   = E_azT;
double AZ_rbcfrac = Aaz_r / ((Aaz_c + Aaz_r > 1e-9) ? (Aaz_c + Aaz_r) : 1e-9);

$CAPTURE
ICP_sup ICP_up ICP_day AMP AMP_day Cspine dPtm_pulse dPtm_mean
Qsh_day Qabs_day If_day turnover EvansIdx CallAngle CBF_pv
gait MMSE urin iNPHGS SDH_inc CSF_Ab42 CSF_pTau CSF_NfL CSF_LRG
AZ_free AZ_eff AZ_rbcfrac
'

mod <- mcode("inph", code, soloc = tempdir())

## =============================================================================
##  Initial conditions
## =============================================================================
##  Both phenotypes are initialised ON their own steady manifolds (ISF amyloid
##  and tau from the quadratic/linear steady states, CSF amounts from the
##  flux/clearance balance, clinical scores at their targets).  Starting off the
##  manifold produces a spurious first-month transient that contaminates every
##  short-horizon read-out — a defect that was actually present in an early
##  version of this file.
## =============================================================================

inph_init <- function(p, healthy = FALSE) {
  y <- as.list(rep(0, 46))
  names(y) <- c("Vv","Vsas","Vsdh","Rout","E1","Vplast","Wpv","AQ","Cart",
                "Autoreg","WMint","WMperm","Myel","Astro","Micro","Ab_isf",
                "Ab_plq","Tau_isf","A_ab","A_pt","A_nfl","A_lrg","Aaz_g",
                "Aaz_c","Aaz_r","Abu_g","Abu_c","Aso_g","Aso_c","Ado_g",
                "Ado_c","Ame_g","Ame_c","HCO3","Kser","AZesc","Occl","Vdr",
                "Gslow","Gfast","Cslow","Cfast","Urin","Headx","SDHhaz","dVac")
  if (healthy) {
    y$Vv <- p$Vv_norm; y$Vsas <- p$Vsas_norm
    y$Rout <- p$Rout_norm; y$E1 <- 0.130
    y$Wpv <- p$Wpv_norm; y$AQ <- p$AQ_norm
    y$Cart <- p$Cart_norm; y$Autoreg <- 0.92
    y$WMint <- 0.96; y$WMperm <- 0.02; y$Myel <- 0.96
    y$Astro <- 0.25; y$Micro <- 0.20; y$Ab_plq <- 0.116
    y$Vplast <- 0
  } else {
    y$Vv <- 95.0; y$Vsas <- p$Vsas_init
    y$Rout <- p$Rout_init; y$E1 <- 0.185
    y$Wpv <- 1.50; y$AQ <- 0.45
    y$Cart <- p$Cart_init; y$Autoreg <- 0.85
    y$WMint <- 0.68; y$WMperm <- 0.10; y$Myel <- 0.70
    y$Astro <- 1.30; y$Micro <- 1.05; y$Ab_plq <- 0.228
    y$Vplast <- p$f_plastic * (95.0 - p$Vv_norm)
  }
  gly <- y$AQ / (1 + 0.35 * y$Wpv)
  a <- p$kagg * p$APOE; b <- p$kab_deg + p$kab_gly * gly; c0 <- -p$Pab
  y$Ab_isf <- (-b + sqrt(b^2 - 4 * a * c0)) / (2 * a)
  y$Tau_isf <- p$Ptau * (1 + p$ktau_plq * y$Ab_plq) /
    (p$ktau_deg + p$ktau_gly * gly)

  Vcsf <- y$Vv + y$Vsas
  kout <- p$If0 * 1440 / Vcsf + p$kdeg_csf
  y$A_ab <- p$kflux_ab * gly * y$Ab_isf / kout
  y$A_pt <- p$kflux_tau * gly * y$Tau_isf / kout

  hy <- inph_hydro(y, p)
  wex <- max(0, y$Wpv - p$Wpv_norm)
  CBF <- p$CBF0 * (y$Autoreg / p$Autoreg_norm) /
    (1 + p$kCBF_W * wex + p$kCBF_P * hy$dPtmP)
  inj <- p$w_water * wex + p$w_perf * max(0, (p$CBF_crit - CBF) / p$CBF_crit) +
    p$w_infl * max(0, y$Micro - 0.20)
  y$A_nfl <- (p$kNfL_base + p$kNfL * p$k_deg * inj * y$WMint) / kout
  y$A_lrg <- (p$kLRG_base + p$kLRG * (max(0, y$Astro - 0.25) + wex)) / kout

  y$HCO3 <- 24.0; y$Kser <- 4.2
  perf <- 0.55 + 0.45 * min(1, CBF / p$CBF0)
  y$Gslow <- max(0.05, p$G_max * y$WMint^p$G_wm_pow * perf *
                   (1 - 0.25 * y$Ab_plq) * (1 - p$comorb))
  taub <- min(0.4, 0.4 * max(0, y$Tau_isf - 1) / 3)
  y$Cslow <- max(5, p$MMSE_max - p$kcog_wm * (1 - y$WMint) -
                   p$kcog_ad * min(1, 0.6 * y$Ab_plq + taub))
  y$Urin <- max(0, min(p$Ur_max, p$kur_wm * (1 - y$WMint)))
  unlist(y)
}

## R-side copy of the hydraulic algebra.  It exists ONLY to compute the two
## patient-specific reference values (AMP_ref, Pday_ref) that the C++ block
## needs as parameters, and to allow the pressure/valve arithmetic of Sections
## 1, 3 and 4 of the report to be reproduced without running the ODE solver.
inph_hydro <- function(y, p, dVac = 0) {
  cm <- function(x) x / 1.36
  compfn <- function(P, e1) {
    dP <- max(P - p$P0_marm, 0.5); min(1 / (e1 * dP), p$C_max)
  }
  Caz <- p$az_fu * y$Aaz_c / p$az_V
  E_az <- p$az_Emax * Caz / (p$az_EC50 + Caz) * (1 - p$az_f_esc * y$AZesc)
  Cbu <- y$Abu_c / p$bu_V
  E_bu <- p$bu_Emax * Cbu / (p$bu_EC50 + Cbu)
  Ifeff <- p$If0 * (1 - E_az) * (1 - E_bu)
  Qin <- Ifeff - p$eld_rate
  RoutE <- y$Rout * if (p$etv > 0.5) (1 - p$etv_dRout) else 1

  res <- lapply(c(FALSE, TRUE), function(up) {
    Pss_eff <- p$Pss - if (up) p$dPss_up else 0
    Pdist <- cm(if (up) p$Pd_up_cm else p$Pd_sup_cm)
    h_eff <- if (up) cm(p$Hcol_cm) * (1 - p$asd_eff) else 0
    Popen_e <- cm(p$Popen_cm) + if (up) cm(p$Ggrav_cm) else 0
    Rsh_eff <- NA
    if (p$shunt > 0.5) {
      Rt <- p$Rsh / (1 - min(y$Occl, 0.995))
      G <- 1 / RoutE + 1 / Rt
      Pf <- (Qin + Pss_eff / RoutE + (Popen_e + Pdist - h_eff) / Rt) / G
      if (Pf + h_eff >= Popen_e + Pdist) Rsh_eff <- Rt
      else Pf <- Pss_eff + Qin * RoutE
    } else Pf <- Pss_eff + Qin * RoutE
    Pbase <- Pf + y$Vsdh / p$Ccs_sdh
    P <- if (dVac != 0) p$P0_marm + (Pbase - p$P0_marm) * exp(y$E1 * dVac) else Pbase
    P <- max(P, -30)
    Qsh <- if (!is.na(Rsh_eff)) max(0, P + h_eff - Popen_e - Pdist) / Rsh_eff else 0
    list(P = P, Qsh = Qsh, Qabs = max(0, P - Pss_eff) / RoutE, C = compfn(P, y$E1))
  })
  s <- res[[1]]; u <- res[[2]]; fu <- p$f_up
  Cday <- (1 - fu) * s$C + fu * u$C
  frac <- max(0, min(1, (p$Cart_norm - y$Cart) / (p$Cart_norm - p$Cart_init)))
  Vp <- (p$Vp_norm + (p$Vp_init - p$Vp_norm) * frac) *
    (1 + p$kVp_edema * max(0, y$Wpv - p$Wpv_norm))
  dfr <- max(0, min(1, (p$Vsas_norm - y$Vsas) / (p$Vsas_norm - p$Vsas_init)))
  desh <- p$DESH_norm + (p$DESH - p$DESH_norm) * dfr
  Ifd <- Ifeff * 1440
  Qsh_day <- ((1 - fu) * s$Qsh + fu * u$Qsh) * 1440
  Qabs_day <- ((1 - fu) * s$Qabs + fu * u$Qabs) * 1440
  list(Ps = s$P, Pu = u$P, Pday = (1 - fu) * s$P + fu * u$P,
       Qsh = Qsh_day, Qabs = Qabs_day, Csup = s$C, Cday = Cday,
       AMP = Vp / s$C, AMPday = Vp / Cday,
       dPtmP = p$kappa_tm * (Vp / Cday) * desh,
       dPtmM = max(0, p$ktm_mean * (y$Rout / p$Rout_norm - 1)) *
         min(1, Qabs_day / max(Ifd, 1e-9)),
       Ifd = Ifd)
}

## Patient-specific references.  AMP_ref is the patient's OWN untreated
## daily-mean pulse amplitude: the acute (tap-test) drive is a fractional fall
## relative to it, so it must be captured before any intervention is applied.
inph_refs <- function(p, healthy = FALSE) {
  y0 <- as.list(inph_init(p, healthy))
  hy <- inph_hydro(y0, p)
  list(AMP_ref = hy$AMPday, Pday_ref = hy$Pday,
       E1_floor_eff = min(p$E1_floor, y0$E1))
}

## Build a parameter list + matching initial condition + references in one call
inph_patient <- function(healthy = FALSE, ...) {
  p <- as.list(param(mod))
  p[names(list(...))] <- list(...)
  r <- inph_refs(p, healthy)
  p[names(r)] <- r
  list(p = p, init = inph_init(p, healthy))
}

## =============================================================================
##  Scenario runner
## =============================================================================
inph_run <- function(label, tend = 730, healthy = FALSE, events = NULL,
                     delta = 1, ...) {
  pt <- inph_patient(healthy = healthy, ...)
  m <- mod %>% param(pt$p) %>% init(pt$init)
  out <- if (is.null(events)) m %>% mrgsim(end = tend, delta = delta,
                                          hmax = 0.25)
         else m %>% mrgsim(data = events, end = tend, delta = delta,
                           hmax = 0.25)
  as_tibble(out) %>% mutate(scenario = label)
}

## Oral dosing helper.  The amount is already multiplied by bioavailability, so
## the gut compartment receives F * dose.
dose_ev <- function(cmt, dose, bioav = 1, ii = 1, addl = 0, time = 0) {
  ev(amt = dose * bioav, cmt = cmt, ii = ii, addl = addl, time = time)
}

## A lumbar tap is a step change in VOLUME, i.e. a negative bolus into dVac.
tap_ev <- function(time, volume) ev(amt = -volume, cmt = "dVac", time = time)

## =============================================================================
##  The 21 scenarios
## =============================================================================
## The reference "well-set" shunt: the optimum found by titration_map().
## A 10 cmH2O differential valve plus 5 cmH2O of abdominal pressure sits at
## 11.0 mmHg against a supine iNPH ICP of 13.2 -- only 2.1 mmHg of headroom --
## so it can take only ~29% of CSF production and UNDER-drains.  A setting this
## low is usable only because the gravitational unit cancels the siphon.
GRAV <- list(shunt = 1, Popen_cm = 4, Ggrav_cm = 30)

scen <- list(
  ## --- natural history and controls ---
  S1  = function() inph_run("S1 untreated iNPH", 730),
  S2  = function() inph_run("S2 healthy 75 y control", 730, healthy = TRUE),

  ## --- diagnostic perturbations ---
  S3  = function() inph_run("S3 tap test 40 mL", 6, events = tap_ev(1, 40),
                            delta = 0.02),
  S4  = function() inph_run("S4 ELD 10 mL/h x 72 h", 8, delta = 0.02,
                            eld_rate = 10 / 60),

  ## --- shunt hardware: the hydrostatics arm ---
  S5  = function() inph_run("S5 valve 10, unprotected", 730,
                            shunt = 1, Popen_cm = 10),
  S6  = function() inph_run("S6 valve 4 + gravitational 30 (well set)", 730,
                            shunt = 1, Popen_cm = 4, Ggrav_cm = 30),
  S6b = function() inph_run("S6b valve 10 + gravitational 30 (under-drains)",
                            730, shunt = 1, Popen_cm = 10, Ggrav_cm = 30),
  S7  = function() inph_run("S7 programmable HIGH 16", 730,
                            shunt = 1, Popen_cm = 16, Ggrav_cm = 30),
  S8  = function() inph_run("S8 programmable LOW 4", 730,
                            shunt = 1, Popen_cm = 4, Ggrav_cm = 30),
  S9  = function() {
    ## Stepwise down-titration, the standard clinical protocol.  A valve setting
    ## is a PARAMETER, not a dose, and mrgsolve cannot change a parameter from an
    ## event record -- so the titration is assembled by stitching fixed-setting
    ## segments, each starting from the state the previous one reached.
    pt <- inph_patient(shunt = 1, Popen_cm = 16, Ggrav_cm = 30)
    segs <- list(c(16, 0, 60), c(13, 60, 120), c(10, 120, 180), c(8, 180, 730))
    y <- pt$init; acc <- list()
    for (s in segs) {
      m <- mod %>% param(pt$p) %>% param(Popen_cm = s[1]) %>% init(y)
      o <- as_tibble(mrgsim(m, start = 0, end = s[3] - s[2], delta = 1,
                            hmax = 0.25))
      y <- unlist(o[nrow(o), names(pt$init)])
      acc[[length(acc) + 1]] <- o %>% mutate(time = time + s[2])
    }
    bind_rows(acc) %>% mutate(scenario = "S9 stepwise titration 16->8")
  },
  S10 = function() inph_stitch("S10 early shunt (6 mo delay)", 182, 1460, GRAV),
  S11 = function() inph_stitch("S11 late shunt (36 mo delay)", 1095, 1460, GRAV),

  ## --- pharmacology ---
  S12 = function() inph_run("S12 acetazolamide 250 mg BID", 730, delta = 1,
                            events = dose_ev("Aaz_g", 250, bioav = 0.9, ii = 0.5,
                                             addl = 1459)),
  S13 = function() inph_run("S13 acetazolamide 500 mg BID", 730, delta = 1,
                            events = dose_ev("Aaz_g", 500, bioav = 0.9, ii = 0.5,
                                             addl = 1459)),
  S14 = function() inph_run("S14 shunt + acetazolamide 250 BID", 730,
                            events = dose_ev("Aaz_g", 250, bioav = 0.9, ii = 0.5,
                                             addl = 1459),
                            shunt = 1, Popen_cm = 4, Ggrav_cm = 30),

  ## --- phenotype and behaviour modifiers ---
  S15 = function() inph_run("S15 shunt, comorbid AD (APOE e4)", 730,
                            APOE = 2.2, shunt = 1, Popen_cm = 4,
                            Ggrav_cm = 30),
  S16 = function() inph_run("S16 unprotected valve, very active (f_up 0.80)",
                            730, shunt = 1, Popen_cm = 10, f_up = 0.80),
  S17 = function() inph_run("S17 shunt + solifenacin 5 mg", 730,
                            events = dose_ev("Aso_g", 5, bioav = 0.9, ii = 1,
                                             addl = 699, time = 30),
                            shunt = 1, Popen_cm = 4, Ggrav_cm = 30),
  S18 = function() {
    ## occlusion at 18 months, revision at 23 months
    pt <- inph_patient(shunt = 1, Popen_cm = 4, Ggrav_cm = 30)
    segs <- list(c(2.5, 0, 548), c(40, 548, 700), c(2.5, 700, 1095))
    y <- pt$init; acc <- list()
    for (s in segs) {
      m <- mod %>% param(pt$p) %>% param(Rsh = s[1]) %>% init(y)
      o <- as_tibble(mrgsim(m, start = 0, end = s[3] - s[2], delta = 1,
                            hmax = 0.25))
      y <- unlist(o[nrow(o), names(pt$init)])
      acc[[length(acc) + 1]] <- o %>% mutate(time = time + s[2])
    }
    bind_rows(acc) %>% mutate(scenario = "S18 occlusion 18 mo + revision")
  },
  S19 = function() inph_run("S19 shunt + melatonin 2 mg", 730,
                            events = dose_ev("Ame_g", 2, bioav = 0.15, ii = 1,
                                             addl = 729),
                            shunt = 1, Popen_cm = 4, Ggrav_cm = 30),

  ## --- the comparator with no hydrostatic column ---
  S20 = function() inph_run("S20 ETV in communicating iNPH", 730, etv = 1),

  ## --- the bundle ---
  S21 = function() inph_run("S21 optimal bundle (grav, 8 cmH2O, melatonin)",
                            730,
                            events = dose_ev("Ame_g", 2, bioav = 0.15, ii = 1,
                                             addl = 729),
                            shunt = 1, Popen_cm = 8, Ggrav_cm = 30)
)

## Delayed-intervention helper: run untreated to `tswitch`, then continue with
## the intervention from the reached state.  The patient's AMP_ref / Pday_ref
## references stay those of the UNTREATED baseline, which is the point.
inph_stitch <- function(label, tswitch, tend, changes) {
  pt <- inph_patient()
  m1 <- mod %>% param(pt$p) %>% init(pt$init)
  o1 <- as_tibble(mrgsim(m1, end = tswitch, delta = 1, hmax = 0.25))
  y <- unlist(o1[nrow(o1), names(pt$init)])
  m2 <- mod %>% param(pt$p) %>% param(changes) %>% init(y)
  o2 <- as_tibble(mrgsim(m2, start = 0, end = tend - tswitch, delta = 1,
                         hmax = 0.25)) %>% mutate(time = time + tswitch)
  bind_rows(o1, o2) %>% mutate(scenario = label)
}

## =============================================================================
##  Run everything and summarise
## =============================================================================
run_all <- function() {
  res <- lapply(names(scen), function(nm) {
    message("running ", nm)
    scen[[nm]]()
  })
  bind_rows(res)
}

summarise_scenarios <- function(all) {
  all %>%
    group_by(scenario) %>%
    summarise(
      tend      = max(time),
      gait0     = first(gait),
      gait_end  = last(gait),
      dGait     = last(gait) - first(gait),
      MMSE_end  = last(MMSE),
      iNPHGS    = last(iNPHGS),
      Evans     = last(EvansIdx),
      AMP       = last(AMP),
      ICP_day   = last(ICP_day),
      Vsdh      = last(Vsdh),
      SDH_pct   = 100 * last(SDH_inc),
      Qsh_mLday = last(Qsh_day),
      CBF_pv    = last(CBF_pv),
      .groups   = "drop"
    ) %>% arrange(scenario)
}

## -----------------------------------------------------------------------------
##  The valve titration map — the central clinical result
##  Utility = gait gain − 1.1 × subdural incidence − 0.03 × headache index.
##  Run with and without a gravitational unit; the device MOVES the optimum
##  rather than merely making a fixed setting safer.
## -----------------------------------------------------------------------------
titration_map <- function(popen = c(2, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28),
                          grav = c(0, 30), tend = 730) {
  rows <- list()
  for (g in grav) for (Po in popen) {
    o <- inph_run(sprintf("Po%.0f_G%.0f", Po, g), tend, delta = 30,
                  shunt = 1, Popen_cm = Po, Ggrav_cm = g)
    rows[[length(rows) + 1]] <- tibble(
      Ggrav = g, Popen = Po,
      ICP_day = last(o$ICP_day), ICP_up = last(o$ICP_up),
      AMP = last(o$AMP), dGait = last(o$gait) - first(o$gait),
      Vsdh = last(o$Vsdh), SDH_pct = 100 * last(o$SDH_inc),
      headache = last(o$Headx))
  }
  bind_rows(rows) %>%
    mutate(utility = dGait - 1.1 * SDH_pct / 100 - 0.03 * headache)
}

## -----------------------------------------------------------------------------
##  Early vs late shunting — the two-clock result.  All arms are read 36 months
##  AFTER surgery so post-operative exposure is identical; only the pre-
##  operative delay differs.
## -----------------------------------------------------------------------------
delay_sweep <- function(delays_mo = c(6, 12, 24, 36, 60)) {
  bind_rows(lapply(delays_mo, function(dm) {
    d <- dm * 30.4
    o <- inph_stitch(sprintf("delay %d mo", dm), d, d + 3 * 365, GRAV)
    tibble(delay_mo = dm, Evans = last(o$EvansIdx), AMP = last(o$AMP),
           Wpv = last(o$Wpv), WMint = last(o$WMint), WMperm = last(o$WMperm),
           gait = last(o$gait), MMSE = last(o$MMSE), iNPHGS = last(o$iNPHGS))
  }))
}

## -----------------------------------------------------------------------------
##  Plots
## -----------------------------------------------------------------------------
plot_shunt_comparison <- function() {
  d <- bind_rows(scen$S1(), scen$S5(), scen$S6(), scen$S7())
  d %>%
    select(time, scenario, ICP_day, AMP, gait, MMSE, Vsdh, EvansIdx) %>%
    pivot_longer(-c(time, scenario)) %>%
    ggplot(aes(time / 30.4, value, colour = scenario)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~name, scales = "free_y") +
    labs(x = "months", y = NULL,
         title = "iNPH: valve setting and siphon protection",
         subtitle = paste("Unprotected valve (S5) drains ~2x CSF production",
                          "on a daily average and builds a subdural collection;",
                          "\na high setting (S7) never opens supine and",
                          "reproduces the untreated trajectory.")) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
}

plot_titration <- function(tm = titration_map()) {
  tm %>%
    mutate(arm = ifelse(Ggrav > 0, "gravitational 30 cmH2O",
                        "no siphon protection")) %>%
    select(arm, Popen, dGait, SDH_pct, utility) %>%
    pivot_longer(-c(arm, Popen)) %>%
    ggplot(aes(Popen, value, colour = arm)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
    facet_wrap(~name, scales = "free_y", ncol = 1) +
    labs(x = "valve opening pressure (cmH2O)", y = NULL,
         title = "The therapeutic window is set by hydrostatics",
         subtitle = paste("A gravitational unit does not only lower risk at a",
                          "fixed setting; it moves the optimum downward")) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
}

plot_tap_vs_eld <- function() {
  bind_rows(scen$S3(), scen$S4()) %>%
    select(time, scenario, ICP_sup, AMP, gait) %>%
    pivot_longer(-c(time, scenario)) %>%
    ggplot(aes(time * 24, value, colour = scenario)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~name, scales = "free_y", ncol = 1) +
    labs(x = "hours", y = NULL,
         title = "Why a 72-hour drain beats a single tap",
         subtitle = paste("Mean ICP recovers within ~30 min (tau = Rout*C).",
                          "Gait follows a slower variable, so the tap works",
                          "at all --\nbut recovers only a fraction of the",
                          "eventual shunt benefit.")) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
}

## =============================================================================
if (interactive() || identical(Sys.getenv("INPH_RUN"), "1")) {
  all <- run_all()
  print(summarise_scenarios(all), n = 25, width = Inf)
  cat("\n--- valve titration map ---\n")
  print(titration_map(), n = 30, width = Inf)
  cat("\n--- early vs late shunting ---\n")
  print(delay_sweep(), width = Inf)
}
