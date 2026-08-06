## =============================================================================
## Hemolytic Disease of the Fetus and Newborn (HDFN)
## mrgsolve implementation -- 45 differential equations
## =============================================================================
##
## THE ORGANISING IDENTITY
## ----------------------
## Fetal red-cell destruction is a PRODUCT of three factors, and every therapy
## in this disease owns exactly one of them:
##
##     destruction = kops * A * f_ag * (Ro + Rr),  cleared with saturation
##
##     A     = Fa1 + 1.7*Fa3           antibody DELIVERED to the fetus.
##                                     Owned by nipocalimab, plasmapheresis,
##                                     IVIG (FcRn competition) and -- upstream
##                                     of all of them -- anti-D prophylaxis.
##     f_ag  = (Ro+Rr)/(Ro+Rr+Rd)      the antigen-POSITIVE fraction.  Owned by
##                                     intrauterine transfusion: donor cells are
##                                     antigen-negative, so a transfusion is a
##                                     substrate dilution and not merely a
##                                     haemoglobin top-up.
##     M     = Res*(1-FcgammaR block)  clearance capacity.  Owned by high-dose
##                                     IVIG and by splenic mass.
##
## WHAT IS FITTED (eight numbers, all published summary statistics)
## ---------------------------------------------------------------
##   v0, g_pl        placental conveyor   <- Malek 1996 (PMID 8955500): fetal:
##                                          maternal IgG 0.075 at 19.5 wk and
##                                          1.25 at 39 wk
##   kops            destruction potency  <- anti-D 15 IU/mL needs its first IUT
##                                          at 26 wk (Nishie 2012 cohort mean
##                                          26.1 wk, PMID 22949399)
##   ksens, fmh_ante sensitisation        <- 16% unprophylaxed, 1.6% with
##                                          postpartum 300 ug only
##   ugt_birth, ugt_t50  conjugation      <- HEALTHY term newborn peaks at
##                                          8 mg/dL TSB at ~4 d
##   emh_thresh      extramedullary switch<- overt ascites at ~5-6 g/dL
##
## WHAT IS PREDICTED (never shown to the fitter)
## --------------------------------------------
##   * the flattening post-transfusion decline, 0.29 -> 0.18 -> 0.17 g/dL/day
##     (reported 0.40, SD 0.25, between the first two procedures -- the model
##     under-predicts the LEVEL by a third and reproduces the ORDER; see
##     hdfn_reference_output.txt section 3, which also shows that matching the
##     level requires donor red cells to survive 30-45 d rather than 70)
##   * MCA-PSV 1.5 MoM is the arithmetic image of Hb 0.65 MoM, and the 12%
##     false-positive rate of Mari 2000 follows from a 10% measurement CV
##   * hydrops appears at Hb MoM = 1/co_max = 0.45
##   * 0.1-0.4% residual sensitisation under antenatal + postpartum prophylaxis
##   * the UNITY nipocalimab result (7/13 = 54% IUT-free live birth >= 32 wk)
##   * anti-K anaemia with LOW bilirubin and low dOD450; anti-D with high
##   * why a fetus with 8 g/dL of haemolysis is not jaundiced and a newborn is
##
## NOTE ON PROVENANCE.  hdfn_python_reference.py is the executed reference
## implementation; this file mirrors it equation for equation.  No R toolchain
## was available in the environment where the model was built, so THIS FILE HAS
## NOT ITSELF BEEN RUN.  Where the two disagree, the Python file is correct.
##
## Units: time in DAYS (gestational age = TIME/7 weeks, TIME starts at ga0*7).
##        antibody amounts in IU, IgG in g, nipocalimab in mg, haemoglobin
##        masses in g, volumes in mL, bilirubin in mg, concentrations derived.
## =============================================================================

library(mrgsolve)
library(dplyr)

code <- '
$PROB
# HDFN: destruction = A x f_ag x M

$PARAM @annotated
// ---- maternal IgG homeostasis ---------------------------------------------
igg0      : 10.0   : maternal total IgG at 12 wk (g/L)
k_int     : 0.120  : IgG internalisation rate (1/d)
phi_rescue: 0.967  : fraction of internalised IgG rescued by free FcRn
K_igg     : 30.0   : IgG concentration at half-saturation of FcRn (g/L)
K_nip     : 8.0    : nipocalimab at half-saturation of FcRn (mg/L)
// ---- anti-D humoral response ----------------------------------------------
ad_iu_per_ug : 5.0 : IU of anti-D per ug (300 ug = 1500 IU)
k_pc      : 0.05   : plasma-cell return to set point (1/d)
sec_pc    : 1.0    : anti-D secretion per unit plasma cell (IU/d)
mpc_set   : 0.0    : plasma-cell set point sustaining the baseline titre
k_bm      : 0.002  : memory B-cell decay (1/d)
k_prime   : 0.60   : memory priming per mL.d of uncoated fetal RBC exposure
k_boost   : 0.35   : memory B to plasma cell differentiation (1/d)
// ---- fetomaternal haemorrhage and prophylaxis -----------------------------
k_sen_free  : 0.0099 : senescence of UNCOATED fetal D+ RBC in mother (1/d)
k_clear_coat: 6.0    : splenic clearance of coated fetal RBC (1/d)
iu_per_ml_rbc: 100.0 : IU anti-D needed to coat 1 mL fetal RBC (= 20 ug/mL)
ksens     : 0.00778122 : FITTED priming per mL.d of exposure
dev       : 0.98   : STRUCTURAL immune-deviation efficiency of coating
k_a_rhig  : 0.35   : absorption of the intramuscular RhIG depot (1/d)
// ---- nipocalimab ----------------------------------------------------------
nip_V1    : 3.0    : central volume (L)
nip_V2    : 4.0    : peripheral volume (L)
nip_Q     : 0.80   : intercompartmental clearance (L/d)
nip_CLlin : 0.85   : non-saturable clearance (L/d)
nip_CLfcrn: 0.25   : FcRn-mediated clearance at low concentration (L/d)
nip_pl_pen: 0.15   : placental interstitial penetration (fraction of plasma)
// ---- IVIG -----------------------------------------------------------------
ivig_compete: 1.0  : does IVIG-derived IgG compete for FcRn? (0 isolates the
//                     FcgammaR mechanism; zeroing the IVIG POOL instead removes
//                     both mechanisms at once, which is how the first version of
//                     that decomposition reached the wrong conclusion)
ivig_fcgr : 0.55   : maximal fetal FcgammaR blockade
ivig_K    : 6.0    : IVIG-derived IgG for half-maximal blockade (g/L)
// ---- placental conveyor (FITTED) ------------------------------------------
v0        : 0.000206568 : conveyor capacity at 20 wk
g_pl      : 0.2677 : exponential growth of conveyor capacity (1/wk)
sub3      : 0.60   : IgG3 transfer relative to IgG1
// ---- fetal IgG ------------------------------------------------------------
vd_igg_f  : 2.0    : fetal IgG distribution volume (x plasma volume)
k_int_f   : 0.10   : fetal IgG internalisation (1/d)
// ---- red cells ------------------------------------------------------------
kbv       : 0.105  : fetoplacental blood volume per g EFW (mL/g)
bv_neo    : 92.0   : blood volume per kg after birth (mL/kg)
k_vol     : 0.60   : plasma volume regulation (1/d)
vpl_expand: 0.50   : plasma volume expansion at zero haematocrit
mchc      : 32.8   : mean corpuscular haemoglobin concentration (g/dL)
kops      : 0.106972 : FITTED opsonisation rate per (IU/mL x site density)
Kres      : 0.30   : STRUCTURAL opsonised Hb at half-maximal clearance (g)
vmax_res  : 1.0    : maximal clearance per unit RES (g Hb/d)
k_res_grow: 0.10   : RES expansion (1/d)
res_max   : 4.0    : maximal RES capacity (x normal)
t_own     : 75.0   : fetal red-cell lifespan (d)
t_don     : 70.0   : donor red-cell lifespan in the fetus (d)
k_mat     : 0.5    : reticulocyte maturation (1/d)
site_D    : 1.0    : relative antigen density (R1r 1.0, R2R2 3.0, K 0.2)
pot3      : 1.7    : ADCC potency of IgG3 relative to IgG1 per IU
// ---- erythropoiesis -------------------------------------------------------
epo0      : 15.0   : fetal EPO at reference Hb (mU/mL)
epo_n     : 5.1    : EPO exponent on 1/HbMoM (BIDIRECTIONAL)
epo_cap   : 40.0   : maximal fold rise of EPO
epo_pow   : 1.2    : EPO to progenitor recruitment exponent
k_epo     : 6.0    : EPO turnover (1/d)
k_prog    : 0.55   : progenitor recruitment at reference EPO (1/d)
prog_max  : 5.0    : maximal marrow progenitor expansion
k_prog_d  : 0.481  : progenitor loss (1/d) -- set so Prog=1 is the steady state
k_prod    : 1.0    : production as a multiple of the requirement
k_emh     : 0.09   : recruitment of hepatic extramedullary erythropoiesis (1/d)
k_emh_off : 0.06   : resolution of EMH (1/d)
emh_thresh: 4.0    : FITTED progenitor expansion at which EMH is recruited
emh_prod  : 0.60   : red cells made per unit EMH
kell_kill : 0.0    : anti-K killing of PROGENITORS per (IU/mL) (1/d)
// ---- liver, albumin, hydrops ---------------------------------------------
alb0      : 1.90   : fetal albumin at 20 wk (g/dL)
alb_slope : 0.070  : rise in albumin per week (g/dL/wk)
k_alb_cat : 0.115  : albumin turnover (1/d)
emh_alb   : 0.40   : maximal fractional suppression of albumin synthesis
k_portal  : 0.80   : umbilical venous pressure per unit EMH (mmHg)
cvp0      : 4.5    : baseline umbilical venous pressure (mmHg)
lp_s      : 6.0    : filtration coefficient (mL/d/mmHg/kg)
pi_frac   : 0.55   : interstitial oncotic pressure as fraction of plasma
asc_compl : 150.0  : ascites per mmHg of interstitial pressure (mL/mmHg)
p_i       : 1.0    : baseline interstitial hydrostatic pressure (mmHg)
sigma     : 0.85   : oncotic reflection coefficient
alb_pi    : 4.5    : plasma oncotic pressure per g/dL albumin (mmHg)
k_lymph_p : 0.50   : lymphatic safety factor per mmHg
lymph_cvp : 0.16   : suppression of lymph return per mmHg of venous pressure
asc_hydrops: 25.0  : ascites that is sonographically overt at 900 g (mL)
k_perm    : 0.35   : permeability rise per unit hypoxia
// ---- oxygen and cardiac ---------------------------------------------------
co_ref    : 450.0  : combined fetal cardiac output (mL/min/kg)
sao2      : 0.55   : fetal arterial saturation
do2_alpha : 1.0    : DERIVED exponent of cerebral flow compensation
visc_k    : 0.0    : viscosity term -- ZERO ON PURPOSE (see header)
co_max    : 2.20   : maximal cardiac output reserve (x baseline)
k_fail    : 0.22   : accumulation of cardiac decompensation (1/d)
k_recov   : 0.12   : recovery of cardiac reserve (1/d)
do2_crit  : 0.55   : relative oxygen delivery below which lactate rises
k_lac     : 2.2    : lactate accumulation at zero reserve (mmol/L/d)
k_lac_out : 1.1    : lactate clearance (1/d)
// ---- bilirubin ------------------------------------------------------------
bil_per_g : 34.0   : mg bilirubin per g haemoglobin catabolised
bil_shunt : 2.0    : shunt bilirubin (mg/kg/d)
k_bil_pl  : 14.0   : PLACENTAL bilirubin clearance (1/d)
vd_bp     : 0.09   : rapidly mixing bilirubin volume (L/kg)
vd_bt     : 0.11   : slowly mixing tissue bilirubin volume (L/kg)
k_bx      : 10.0   : plasma to tissue bilirubin exchange (1/d)
k_amn     : 0.55   : appearance in amniotic fluid (1/d)
k_amn_out : 0.42   : amniotic fluid bilirubin turnover (1/d)
od450_per_mg: 0.052: dOD450 per mg/dL amniotic bilirubin
ugt_birth : 0.060375 : FITTED UGT1A1 activity at birth (fraction of adult)
ugt_t50   : 20.4664  : FITTED half-maturation of UGT1A1 (d)
ugt_n     : 1.5    : Hill coefficient of UGT1A1 maturation
vmax_ugt  : 95.0   : conjugation at full adult activity (mg/d/kg)
km_bil    : 8.0    : Michaelis constant of conjugation (mg/dL)
k_ehc     : 0.22   : enterohepatic recirculation (1/d)
k_gut     : 0.55   : faecal loss of gut bilirubin (1/d)
k_photo   : 0.022  : photoisomerisation per (uW/cm2/nm) (1/d)
k_iso_out : 1.6    : biliary loss of photoisomers (1/d)
k_iso_back: 0.35   : reversion of photoisomer (1/d)
alb_bind  : 0.72   : mg bilirubin bound per g/dL albumin per dL
bind_thresh: 0.65  : B/A ratio at which neurotoxic exposure accumulates
// ---- hazards --------------------------------------------------------------
h0        : 0.00012 : background fetal loss (1/d)
h_hydrops : 0.010   : additional hazard with overt hydrops (1/d)
h_severe  : 0.020   : additional hazard below Hb 0.35 MoM (1/d)
h_acid    : 0.030   : additional hazard below pH 7.15 (1/d)
h_iut     : 0.012   : procedure-related loss per IUT (Zwiers 2017)
ph0       : 7.35    : reference fetal pH
ph_lac    : 0.035   : pH units per mmol/L lactate
// ---- protocol switches (set by the R driver, not by the ODE) --------------
BORN      : 0.0    : 1 after delivery
GA_BIRTH  : 40.0   : gestational age at delivery (wk)
PHOTO     : 0.0    : phototherapy irradiance (uW/cm2/nm)
EHC_ON    : 0.0    : enterohepatic recirculation switch (feeding)

$CMT @annotated
Mig   : maternal total IgG (g)
Ma1   : maternal anti-D IgG1 (IU)
Ma3   : maternal anti-D IgG3 (IU)
Mbm   : maternal memory B cells (rel)
Mpc   : maternal anti-D plasma cells (rel)
Mfree : uncoated fetal D+ RBC in maternal blood (mL)
Mcoat : coated fetal D+ RBC in maternal blood (mL)
Mrhig : passive prophylactic anti-D (IU)
Mdep  : intramuscular RhIG depot (IU)
Sens  : cumulative sensitisation signal (-log survival)
Nc    : nipocalimab central (mg)
Np    : nipocalimab peripheral (mg)
Miv   : IVIG-derived maternal IgG (g)
Fig   : fetal total IgG (g)
Fa1   : fetal anti-D IgG1 (IU)
Fa3   : fetal anti-D IgG3 (IU)
Fnip  : fetal nipocalimab (mg)
EFW   : estimated fetal weight (g)
Vpl   : fetal/neonatal plasma volume (mL)
Ro    : own antigen-POSITIVE red-cell Hb mass (g)
Rd    : donor antigen-NEGATIVE red-cell Hb mass (g)
Rr    : own reticulocyte Hb mass (g)
Rop   : opsonised red-cell Hb mass awaiting clearance (g)
Prog  : erythroid progenitor pool (rel)
EMH   : hepatic extramedullary erythropoiesis (rel)
Epo   : fetal erythropoietin (mU/mL)
Res   : reticuloendothelial capacity (rel)
Alb   : albumin concentration (g/dL)
Asc   : ascites and interstitial oedema (mL)
Card  : cardiac decompensation index (0-1)
Lac   : lactate (mmol/L)
Bil   : plasma bilirubin (mg)
Bex   : extravascular tissue bilirubin (mg)
Bgut  : intestinal bilirubin (mg)
Biso  : photoisomer pool (mg)
Amn   : amniotic-fluid bilirubin (mg)
Ugt   : UGT1A1 activity (fraction of adult)
Bind  : cumulative free-bilirubin neurotoxic exposure
Hzd   : cumulative perinatal-loss hazard
Fer   : transfusional iron load (mg Fe)
Nrbc  : circulating nucleated red cells (rel)
Wgt   : postnatal weight (kg)
Adose : cumulative anti-D delivered to the fetus (IU)
Dest  : cumulative antigen-positive Hb destroyed (g)
SenC  : cumulative senescent Hb loss (g)

$GLOBAL
#define EPS 1e-12

// Hadlock 50th-centile estimated fetal weight, log-cubic fit
double efw_ref(double ga) {
  double g = ga < 14.0 ? 14.0 : (ga > 42.0 ? 42.0 : ga);
  double lp = -0.00013746*g*g*g + 0.00650447*g*g + 0.15310042*g + 1.35293281;
  return exp(lp);
}
// Nicolaides reference fetal haemoglobin
double hb_ref(double ga) {
  double g = ga < 17.0 ? 17.0 : (ga > 42.0 ? 42.0 : ga);
  return 10.9 + 0.1565*(g - 17.0);
}
// postnatal haemoglobin reference: cord Hb is 14% above intrauterine venous
double hb_ref_pn(double ga_birth, double pna) {
  double gb = ga_birth > 40.0 ? 40.0 : ga_birth;
  double hb0 = hb_ref(ga_birth)*1.14;
  double nadir = 9.0 + 0.25*(gb - 32.0);
  double t = pna < 0.0 ? 0.0 : pna;
  return nadir + (hb0 - nadir)*exp(-t/22.0);
}
double alb_ref_f(double ga, double a0, double sl) {
  double g = ga > 40.0 ? 40.0 : ga;
  return a0 + sl*(g - 20.0);
}
double ugt_ref_f(double pna, double u0, double t50, double n) {
  double t = pna < 0.0 ? 0.0 : pna;
  double x = pow(t/t50, n);
  return u0 + (1.0 - u0)*x/(1.0 + x);
}
// Mari 2000 median MCA peak systolic velocity (cm/s)
double psv_med(double ga) {
  double g = ga < 16.0 ? 16.0 : (ga > 42.0 ? 42.0 : ga);
  return exp(2.31 + 0.046*g);
}
// maternal plasma volume (L)
double vm_plasma(double ga) {
  return 2.5*(1.0 + 0.45/(1.0 + exp(-(ga - 24.0)/4.5)));
}

$MAIN
// nothing to initialise here: the R driver sets the initial condition, because
// the initial state depends on the gestational age at which a case is entered
// and on the maternal titre.

$ODE
double ga   = SOLVERTIME/7.0;
double born = BORN > 0.5 ? 1.0 : 0.0;
double pna  = born > 0.5 ? (ga - GA_BIRTH)*7.0 : 0.0;

// ---------------- maternal concentrations ----------------------------------
double Vm   = vm_plasma(born > 0.5 ? GA_BIRTH : ga);
double Cig  = Mig/Vm;                       // g/L
double Civ  = Miv/Vm;                       // g/L
double Cnip = Nc/nip_V1;                    // mg/L
double Ca1  = Ma1/(Vm*1000.0);              // IU/mL
double Ca3  = Ma3/(Vm*1000.0);
double Crh  = Mrhig/(Vm*1000.0);
double fu_m = 1.0/(1.0 + (Cig + ivig_compete*Civ)/K_igg + Cnip/K_nip);

// ---------------- the fetal circulation: plasma PLUS cells ----------------
double efw   = EFW > 80.0 ? EFW : 80.0;
double wt    = born > 0.5 ? Wgt : efw/1000.0;
double hbmass= Ro + Rd + Rr;
double Vrbc  = hbmass/(mchc/100.0);
double Vplx  = Vpl > 1.0 ? Vpl : 1.0;
double Vfp   = Vplx + Vrbc;
double hb    = 100.0*hbmass/Vfp;
double hct   = Vrbc/Vfp;
double Vd_ig = vd_igg_f*Vplx/1000.0;        // L
double Cf_ig = Fig/(Vd_ig + EPS);
double Cf_a1 = Fa1/(Vd_ig*1000.0 + EPS);
double Cf_a3 = Fa3/(Vd_ig*1000.0 + EPS);
double Cf_np = Fnip/(Vd_ig + EPS);

// ---------------- THE THREE FACTORS ---------------------------------------
double A_eff = (Cf_a1 + pot3*Cf_a3)*site_D;
double f_ag  = (Ro + Rr)/(hbmass + EPS);
double fcgr  = ivig_fcgr*Civ/(ivig_K + Civ);
double M_eff = Res*(1.0 - fcgr);

// ---------------- oxygen delivery and the Doppler -------------------------
double hbr    = born > 0.5 ? hb_ref_pn(GA_BIRTH, pna) : hb_ref(ga);
double hbmom  = hb/hbr;
double hbm_fl = hbmom < 0.12 ? 0.12 : hbmom;
double dem_br = pow(1.0/hbm_fl, do2_alpha);          // cerebral, NOT capped
double demand = dem_br < co_max ? dem_br : co_max;    // cardiac, capped
double do2    = demand*hbmom;
double hct_rf = hbr/mchc;
double visc   = exp(-visc_k*(hct - hct_rf));
double psv_mom= dem_br*visc;

// ---------------- Starling balance ----------------------------------------
double cvp    = cvp0 + k_portal*EMH + 6.0*Card;
double pi_p   = alb_pi*Alb;
double hypox  = do2 < do2_crit ? (1.0 - do2/do2_crit) : 0.0;
double kg     = wt > 0.05 ? wt : 0.05;
double lp     = lp_s*kg*(1.0 + k_perm*hypox);
double pc     = 0.55*cvp + 8.0;
double p_i_ef = p_i + Asc/asc_compl;
double jv     = lp*((pc - p_i_ef) - sigma*pi_p*(1.0 - pi_frac));
double pc0    = 0.55*cvp0 + 8.0;
double pi_rf  = alb_pi*alb_ref_f(born > 0.5 ? GA_BIRTH : ga, alb0, alb_slope);
double jv_ref = lp_s*kg*((pc0 - p_i) - sigma*pi_rf*(1.0 - pi_frac));
double lymph  = jv_ref*(1.0 + k_lymph_p*(p_i_ef - p_i))*
                exp(-lymph_cvp*(cvp > cvp0 ? cvp - cvp0 : 0.0));

// ---------------- bilirubin readouts --------------------------------------
double vd_b = vd_bp*wt;
double tsb  = Bil/(vd_b*10.0 + EPS);
double ctis = Bex/(vd_bt*wt*10.0 + EPS);
double ba   = tsb/(Alb > 0.3 ? Alb : 0.3)/8.5;

// ---------------- maternal IgG and anti-D ---------------------------------
double kel_igg = k_int*(1.0 - phi_rescue*fu_m);
double syn_igg = igg0*vm_plasma(12.0)*k_int*(1.0 - phi_rescue/(1.0 + igg0/K_igg));

// ---------------- the placental conveyor ----------------------------------
double cap   = born > 0.5 ? 0.0 : v0*exp(g_pl*(ga - 20.0));
double fu_pl = 1.0/(1.0 + (Cig + ivig_compete*Civ)/K_igg + nip_pl_pen*Cnip/K_nip);
double J_ig  = cap*fu_pl*Cig;
double J_iv  = cap*fu_pl*Civ;
double J_a1  = cap*fu_pl*Ca1*1000.0;
double J_a3  = sub3*cap*fu_pl*Ca3*1000.0;
double J_rh  = cap*fu_pl*Crh*1000.0;
double J_nip = 0.55*cap*fu_pl*Cnip;

// ---------------- fetomaternal haemorrhage and priming --------------------
double need   = iu_per_ml_rbc*(Mfree + Mcoat);
double avail  = Mrhig + Ma1 + Ma3;
double coat_f = (avail + need) > 0.0 ? avail/(avail + need) : 0.0;
double to_coat= 4.0*coat_f*Mfree;
// PRIMING IS THE EXPOSURE INTEGRAL, not a clearance flux.  This is what makes
// the 72-hour window 72 hours: an uncoated fetal cell survives for months.
double ag_pr  = Mfree + (1.0 - dev)*Mcoat;

dxdt_Mig   = syn_igg - kel_igg*Mig - J_ig;
dxdt_Miv   = -kel_igg*Miv - J_iv;
dxdt_Ma1   = sec_pc*Mpc - kel_igg*Ma1 - J_a1;
dxdt_Ma3   = 0.18*sec_pc*Mpc - kel_igg*Ma3 - J_a3;
dxdt_Mrhig = k_a_rhig*Mdep - kel_igg*Mrhig - k_clear_coat*Mcoat*iu_per_ml_rbc
             - J_rh;
dxdt_Mdep  = -k_a_rhig*Mdep;
dxdt_Mfree = -to_coat - k_sen_free*Mfree;
dxdt_Mcoat = to_coat - k_clear_coat*Mcoat;
dxdt_Sens  = ksens*ag_pr;
dxdt_Mbm   = k_prime*ag_pr*(1.0 - Mbm/5.0) - k_bm*Mbm;
// long-lived plasma cells decay back to a SET POINT, not to zero: otherwise
// every sensitised mother's titre collapses during the pregnancy
dxdt_Mpc   = k_boost*Mbm*ag_pr - k_pc*(Mpc - mpc_set);

// ---------------- nipocalimab ---------------------------------------------
double cl_nip = nip_CLlin + nip_CLfcrn*K_nip/(K_nip + Cnip);
dxdt_Nc = -cl_nip*Cnip - nip_Q*(Cnip - Np/nip_V2) - J_nip;
dxdt_Np =  nip_Q*(Cnip - Np/nip_V2);

// ---------------- fetal antibody ------------------------------------------
double kel_f = k_int_f*(1.0 - phi_rescue/(1.0 + Cf_ig/K_igg + Cf_np/K_nip));
dxdt_Fig  = J_ig + J_iv - kel_f*Fig;
dxdt_Fa1  = J_a1 + J_rh - kel_f*Fa1;
dxdt_Fa3  = J_a3 - kel_f*Fa3;
dxdt_Fnip = J_nip - (nip_CLlin/nip_V1*1.5)*Fnip;
dxdt_Adose= J_a1 + J_a3;

// ---------------- growth and plasma volume --------------------------------
double target = efw_ref(ga);
if (born > 0.5) {
  dxdt_EFW = 0.0;
  dxdt_Wgt = 0.018*Wgt*(1.0 - 0.3*Card);
} else {
  dxdt_EFW = (target - EFW)*0.35 +
             (efw_ref(ga + 0.01) - target)/0.07*(1.0 - 0.45*Card);
  dxdt_Wgt = 0.0;
}
double vb_ref = born > 0.5 ? bv_neo*wt : kbv*efw;
double vpl_t  = vb_ref*(1.0 - hct_rf)*
                (1.0 + vpl_expand*(hct < hct_rf ? 1.0 - hct/hct_rf : 0.0));
dxdt_Vpl = k_vol*(vpl_t - Vpl);

// ---------------- THE PRODUCT: destruction --------------------------------
double ops  = kops*A_eff*(Ro + Rr);
double clr  = vmax_res*M_eff*Rop/(Kres + Rop + EPS);
double sen_o= log(2.0)/t_own*Ro;
double sen_d= log(2.0)/t_don*Rd;

// production is a multiple of the REQUIREMENT (senescence + blood-volume
// growth), which is why k_prod = 1 is not a fitted number
double mref, dmref;
if (born > 0.5) {
  mref  = hb_ref_pn(GA_BIRTH, pna)*bv_neo*wt/100.0;
  dmref = (hb_ref_pn(GA_BIRTH, pna + 0.05)*bv_neo*wt/100.0 - mref)/0.05;
} else {
  mref  = hb_ref(ga)*kbv*efw_ref(ga)/100.0;
  dmref = (hb_ref(ga + 0.01)*kbv*efw_ref(ga + 0.01)/100.0 - mref)/0.07;
}
double req  = log(2.0)/t_own*mref + (dmref > 0.0 ? dmref : 0.0);
double prod = k_prod*(Prog + emh_prod*EMH)*req;
double shr  = (Ro + Rr) > EPS ? Ro/(Ro + Rr) : 0.0;

dxdt_Rr  = prod - k_mat*Rr - ops*(1.0 - shr);
dxdt_Ro  = k_mat*Rr - sen_o - ops*shr;
dxdt_Rd  = -sen_d;
dxdt_Rop = ops - clr;
dxdt_Dest= clr;
dxdt_SenC= sen_o + sen_d;
dxdt_Res = k_res_grow*(Rop/(Kres + Rop + EPS))*(res_max - Res)
           - 0.05*(Res - 1.0);

// ---------------- erythropoiesis ------------------------------------------
double epo_r = pow(1.0/(hbmom < 0.20 ? 0.20 : hbmom), epo_n);
double epo_t = epo0*(epo_r < epo_cap ? epo_r : epo_cap);
dxdt_Epo = k_epo*(epo_t - Epo);
double drive = pow(Epo/epo0, epo_pow);
dxdt_Prog = k_prog*drive*(1.0 - Prog/prog_max) - k_prog_d*Prog
            - kell_kill*A_eff*Prog;
dxdt_EMH  = k_emh*(Prog > emh_thresh ? Prog - emh_thresh : 0.0) - k_emh_off*EMH;
dxdt_Nrbc = 0.6*(Prog - Nrbc);

// ---------------- albumin, hydrops, cardiac -------------------------------
double alb_t = alb_ref_f(born > 0.5 ? GA_BIRTH : ga, alb0, alb_slope)*
               (1.0 - emh_alb*EMH/(1.0 + EMH));
dxdt_Alb = k_alb_cat*(alb_t - Alb);
if (born > 0.5) {
  dxdt_Asc = -0.25*Asc;
} else {
  double net = jv - lymph;
  dxdt_Asc = net > -0.25*Asc ? net : -0.25*Asc;
}
dxdt_Card = k_fail*(dem_br/co_max > 1.0 ? dem_br/co_max - 1.0 : 0.0)*
            (1.0 - Card) - k_recov*Card;
dxdt_Lac  = k_lac*hypox - k_lac_out*Lac;

// ---------------- bilirubin -----------------------------------------------
double bil_prod = bil_per_g*(clr + sen_o + sen_d) + bil_shunt*wt;
double conj = 0.0, bil_out = 0.0, amn_in = 0.0, photo = 0.0, ehc = 0.0;
if (born > 0.5) {
  conj    = vmax_ugt*wt*Ugt*tsb/(tsb + km_bil);
  bil_out = conj;
  photo   = k_photo*PHOTO*Bil;
  ehc     = k_ehc*Bgut*EHC_ON;
} else {
  bil_out = k_bil_pl*Bil;
  amn_in  = k_amn*Bil*0.06;
}
double xchg = k_bx*(tsb - ctis)*vd_bp*wt*10.0;
dxdt_Bil  = bil_prod - bil_out - photo + k_iso_back*Biso + ehc - xchg;
dxdt_Bex  = xchg;
dxdt_Biso = photo - k_iso_out*Biso - k_iso_back*Biso;
dxdt_Bgut = conj - k_gut*Bgut - ehc;
dxdt_Amn  = amn_in - k_amn_out*Amn;
dxdt_Ugt  = born > 0.5 ? 3.0*(ugt_ref_f(pna, ugt_birth, ugt_t50, ugt_n) - Ugt)
                       : 0.0;
dxdt_Bind = ba > bind_thresh ? ba - bind_thresh : 0.0;

// ---------------- hazard ---------------------------------------------------
double ph  = ph0 - ph_lac*Lac;
double thr = asc_hydrops*efw/900.0;
dxdt_Hzd = h0 + (Asc > thr ? h_hydrops : 0.0)
              + (hbmom < 0.35 ? h_severe : 0.0)
              + (ph < 7.15 ? h_acid : 0.0);
dxdt_Fer = 0.0;

$TABLE
double ga_o   = TIME/7.0;
double born_o = BORN > 0.5 ? 1.0 : 0.0;
double pna_o  = born_o > 0.5 ? (ga_o - GA_BIRTH)*7.0 : 0.0;
double Vm_o   = vm_plasma(born_o > 0.5 ? GA_BIRTH : ga_o);
double efw_o  = EFW > 80.0 ? EFW : 80.0;
double wt_o   = born_o > 0.5 ? Wgt : efw_o/1000.0;
double mass_o = Ro + Rd + Rr;
double Vrbc_o = mass_o/(mchc/100.0);
double Vfp_o  = (Vpl > 1.0 ? Vpl : 1.0) + Vrbc_o;
double HB     = 100.0*mass_o/Vfp_o;
double HCT    = Vrbc_o/Vfp_o;
double HBREF  = born_o > 0.5 ? hb_ref_pn(GA_BIRTH, pna_o) : hb_ref(ga_o);
double HBMOM  = HB/HBREF;
double PSVMOM = pow(1.0/(HBMOM < 0.12 ? 0.12 : HBMOM), do2_alpha)*
                exp(-visc_k*(HCT - HBREF/mchc));
double PSV    = PSVMOM*psv_med(ga_o);
double Vd_o   = vd_igg_f*(Vpl > 1.0 ? Vpl : 1.0)/1000.0;
double CFA1   = Fa1/(Vd_o*1000.0 + EPS);
double CFA3   = Fa3/(Vd_o*1000.0 + EPS);
double AEFF   = (CFA1 + pot3*CFA3)*site_D;
double FAG    = (Ro + Rr)/(mass_o + EPS);
double TSB    = Bil/(vd_bp*wt_o*10.0 + EPS);
double OD450  = od450_per_mg*Amn/(0.4*efw_o/1000.0 > 0.05 ?
                                  0.4*efw_o/1000.0 : 0.05);
double CA1    = Ma1/(Vm_o*1000.0);
double CIG    = Mig/Vm_o;
double CNIP   = Nc/nip_V1;
double HYDROPS= Asc > asc_hydrops*efw_o/900.0 ? 1.0 : 0.0;
double SURV   = exp(-Hzd);
double BA     = TSB/(Alb > 0.3 ? Alb : 0.3)/8.5;
double PSENS  = 1.0 - exp(-Sens);

$CAPTURE ga_o HB HCT HBMOM PSV PSVMOM AEFF FAG TSB OD450
$CAPTURE CA1 CIG CNIP CFA1 HYDROPS SURV BA PSENS Vfp_o efw_o wt_o
'

hdfn <- mcode("hdfn", code)

## =============================================================================
## INITIAL CONDITIONS
## The state at entry depends on the gestational age and the maternal titre, so
## it is built in R rather than in $MAIN.
## =============================================================================
hdfn_init <- function(mod, ga0 = 13, anti_d_iu = 15, hb_frac = 1.0) {
  p    <- as.list(param(mod))
  Vm   <- 2.5*(1 + 0.45/(1 + exp(-(ga0 - 24)/4.5)))
  efw  <- exp(-0.00013746*ga0^3 + 0.00650447*ga0^2 + 0.15310042*ga0 + 1.35293281)
  hbr  <- 10.9 + 0.1565*(max(ga0, 17) - 17)
  Vb   <- p$kbv*efw
  mass <- hbr*hb_frac*Vb/100
  kel  <- p$k_int*(1 - p$phi_rescue/(1 + p$igg0/p$K_igg))
  Ma1  <- anti_d_iu*Vm*1000*0.85
  list(
    init = c(Mig = p$igg0*Vm, Ma1 = Ma1, Ma3 = anti_d_iu*Vm*1000*0.15,
             Mbm = ifelse(anti_d_iu > 0, 2, 0),
             Mpc = ifelse(anti_d_iu > 0, kel*Ma1/p$sec_pc, 0),
             EFW = efw, Vpl = Vb - mass/(p$mchc/100),
             Ro = mass*0.92, Rr = mass*0.08, Prog = 1, Epo = p$epo0, Res = 1,
             Alb = p$alb0 + p$alb_slope*(min(ga0, 40) - 20),
             Ugt = p$ugt_birth, Nrbc = 1, Wgt = efw/1000,
             Bil = 0.4*p$vd_bp*(efw/1000)*10,
             Bex = 0.4*p$vd_bt*(efw/1000)*10),
    mpc_set = ifelse(anti_d_iu > 0, kel*Ma1/p$sec_pc, 0))
}

## =============================================================================
## INTRAUTERINE TRANSFUSION -- mass balance, not a nomogram
##   (V_rbc + Hct_don*V) / (V_b + V) = Hct_target
## Only 20% of a packed unit is plasma, which is why a transfusion does NOT
## dilute the antibody that is already inside the fetus.
## =============================================================================
iut_volume <- function(state, p, target_hct = 0.45, donor_hct = 0.80) {
  mass <- state[["Ro"]] + state[["Rd"]] + state[["Rr"]]
  Vrbc <- mass/(p$mchc/100)
  Vb   <- max(state[["Vpl"]], 1) + Vrbc
  if (Vrbc/Vb >= target_hct) return(0)
  V <- (target_hct*Vb - Vrbc)/(donor_hct - target_hct)
  min(V, 0.60*Vb)
}

apply_iut <- function(state, p, target_hct = 0.45, donor_hct = 0.80) {
  V <- iut_volume(state, p, target_hct, donor_hct)
  if (V <= 0) return(list(state = state, volume = 0))
  state[["Rd"]]  <- state[["Rd"]] + donor_hct*V*p$mchc/100
  state[["Vpl"]] <- state[["Vpl"]] + (1 - donor_hct)*V
  state[["Fer"]] <- state[["Fer"]] + V*1.08
  state[["Mfree"]] <- state[["Mfree"]] + 0.35   # every IUT is itself an FMH
  state[["Hzd"]] <- state[["Hzd"]] + p$h_iut
  list(state = state, volume = V)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## =============================================================================
## THE DRIVER
## Day-by-day, because the obstetric protocol is a discrete-time controller:
## weekly MCA-PSV, cordocentesis to confirm a first transfusion, then repeat
## transfusions on the calendar (2 weeks, then 3), delivery, then neonatal
## phototherapy / exchange / top-up thresholds.
## =============================================================================
simulate_hdfn <- function(mod, ga0 = 13, ga_deliver = 36, postnatal_days = 70,
                          anti_d_iu = 15, protocol = c("mca", "fixed", "none"),
                          nip_dose = 0, nip_start = 99, nip_stop = 35,
                          ivig_dose = 0, ivig_start = 99, mat_wt = 70,
                          mca_cut = 1.50, iut_target = 0.45, photo_irr = 30,
                          seed = 7) {
  protocol <- match.arg(protocol)
  set.seed(seed)
  p  <- as.list(param(mod))
  I  <- hdfn_init(mod, ga0, anti_d_iu)
  mod <- param(mod, mpc_set = I$mpc_set, BORN = 0, GA_BIRTH = ga_deliver,
               PHOTO = 0, EHC_ON = 0)
  st <- init(mod, as.list(I$init)) %>% as.numeric()
  names(st) <- names(init(mod))
  t   <- ga0*7
  end <- ga_deliver*7 + postnatal_days
  out <- NULL
  iut_log <- list()
  next_iut <- NA_real_
  born <- FALSE
  n_exch <- 0L; n_topup <- 0L; last_exch <- -99
  while (t < end - 1e-9) {
    ga <- t/7
    ## ---- weekly drug administration ----------------------------------------
    if (!born && nip_dose > 0 && ga >= nip_start && ga <= nip_stop &&
        abs(((ga - nip_start)*7) %% 7) < 0.5)
      st[["Nc"]] <- st[["Nc"]] + nip_dose*mat_wt
    if (!born && ivig_dose > 0 && ga >= ivig_start &&
        abs(((ga - ivig_start)*7) %% 7) < 0.5)
      st[["Miv"]] <- st[["Miv"]] + ivig_dose*mat_wt
    ## ---- surveillance and transfusion --------------------------------------
    if (!born && protocol != "none" && ga >= 17.5 && ga <= 35) {
      mass <- st[["Ro"]] + st[["Rd"]] + st[["Rr"]]
      Vrbc <- mass/(p$mchc/100); Vb <- max(st[["Vpl"]], 1) + Vrbc
      hb   <- 100*mass/Vb
      hbr  <- 10.9 + 0.1565*(max(ga, 17) - 17)
      psvm <- 1/max(hb/hbr, 0.12)
      trig <- FALSE
      weekly <- (((ga - ga0)*7) %% 7) < 0.5
      if (protocol == "mca" && weekly) {
        trig <- psvm*(1 + 0.10*rnorm(1)) >= mca_cut
      } else if (protocol == "fixed" && weekly && is.na(next_iut)) {
        ## "fixed" transfuses on the calendar, but something still has to start
        ## the course: the first cordocentesis below 0.85 MoM.
        trig <- hb < 0.85*hbr
      }
      sched <- !is.na(next_iut) && ga >= next_iut
      if (sched) trig <- TRUE
      if (trig && (sched || protocol == "fixed" ||
                   hb*(1 + 0.03*rnorm(1)) < 0.75*hbr)) {
        res <- apply_iut(st, p, iut_target)
        if (res$volume > 0) {
          st <- res$state
          iut_log[[length(iut_log) + 1]] <-
            data.frame(ga = ga, hb_pre = hb, volume = res$volume)
          next_iut <- ga + ifelse(protocol == "fixed" ||
                                  length(iut_log) == 1, 2, 3)
          if (next_iut > 34.5) next_iut <- NA_real_
        }
      }
    }
    ## ---- delivery -----------------------------------------------------------
    if (!born && t >= ga_deliver*7 - 1e-9) {
      born <- TRUE
      st[["Wgt"]] <- max(st[["EFW"]], 80)/1000
      mod <- param(mod, BORN = 1, EHC_ON = 1)
    }
    ## ---- neonatal management ------------------------------------------------
    if (born) {
      hours <- (t - ga_deliver*7)*24
      tsb   <- st[["Bil"]]/(p$vd_bp*st[["Wgt"]]*10)
      mass  <- st[["Ro"]] + st[["Rd"]] + st[["Rr"]]
      Vb    <- max(st[["Vpl"]], 1) + mass/(p$mchc/100)
      mod   <- param(mod, PHOTO = ifelse(tsb > 5 + min(hours, 96)*0.10,
                                         photo_irr, 0))
      ## A double-volume exchange is VOLUME-NEUTRAL: it removes ~85% of
      ## everything in the circulation -- bilirubin, antibody, the infant's own
      ## cells AND any donor cells already there -- and then restores the red
      ## cell volume.  Adding donor cells without removing the old ones makes
      ## the donor pool, its senescence and therefore bilirubin production grow
      ## with every procedure: a positive feedback loop.
      if (tsb > 12 + min(hours, 96)*0.09 && (t - last_exch) > 0.5) {
        fx <- 0.85; keep <- 1 - fx
        st[["Bil"]] <- st[["Bil"]]*keep;  st[["Bex"]] <- st[["Bex"]]*0.65
        st[["Fa1"]] <- st[["Fa1"]]*keep;  st[["Fa3"]] <- st[["Fa3"]]*keep
        st[["Fig"]] <- st[["Fig"]]*keep
        st[["Ro"]]  <- st[["Ro"]]*keep;   st[["Rr"]]  <- st[["Rr"]]*keep
        st[["Rd"]]  <- st[["Rd"]]*keep
        mass_now <- st[["Ro"]] + st[["Rd"]] + st[["Rr"]]
        mass_t   <- 0.50/(1 - 0.50)*max(st[["Vpl"]], 1)*p$mchc/100
        if (mass_t > mass_now) {
          st[["Rd"]]  <- st[["Rd"]] + (mass_t - mass_now)
          st[["Fer"]] <- st[["Fer"]] + (mass_t - mass_now)/(p$mchc/100)*1.08
        }
        n_exch <- n_exch + 1; last_exch <- t
      }
      if (100*mass/Vb < 7.5 && hours > 24) {
        st[["Rd"]]  <- st[["Rd"]] + 0.60*0.15*Vb*p$mchc/100
        st[["Vpl"]] <- st[["Vpl"]] + 0.40*0.15*Vb
        n_topup <- n_topup + 1L
      }
    }
    ## ---- integrate one day --------------------------------------------------
    df <- mod %>% init(as.list(st)) %>%
      mrgsim(start = t, end = t + 1, delta = 1) %>% as.data.frame()
    out <- rbind(out, df[nrow(df), , drop = FALSE])
    st  <- unlist(df[nrow(df), names(st)])
    st[st < 0] <- 0
    t <- t + 1
  }
  list(out = out, iut = do.call(rbind, iut_log), n_iut = length(iut_log),
       n_exch = n_exch, n_topup = n_topup)
}

## =============================================================================
## SIXTEEN SCENARIOS
## The numbers each of these produces are printed in hdfn_reference_output.txt,
## section 13, from the Python reference implementation.
## =============================================================================
hdfn_scenarios <- function(mod) {
  S <- list(
    list(id = 1,  lab = "no antibody (reference fetus)",
         args = list(anti_d_iu = 0, protocol = "none")),
    list(id = 2,  lab = "anti-D 4 IU/mL, no intervention",
         args = list(anti_d_iu = 4, protocol = "none")),
    list(id = 3,  lab = "anti-D 15 IU/mL, no intervention",
         args = list(anti_d_iu = 15, protocol = "none")),
    list(id = 4,  lab = "anti-D 60 IU/mL, no intervention",
         args = list(anti_d_iu = 60, protocol = "none")),
    list(id = 5,  lab = "anti-D 15, MCA-PSV surveillance + IUT",
         args = list(anti_d_iu = 15)),
    list(id = 6,  lab = "anti-D 60, MCA-PSV surveillance + IUT",
         args = list(anti_d_iu = 60)),
    list(id = 7,  lab = "anti-D 60, fixed 2-weekly IUT",
         args = list(anti_d_iu = 60, protocol = "fixed")),
    list(id = 8,  lab = "anti-D 60 + IVIG 1 g/kg/wk from 12 wk",
         args = list(anti_d_iu = 60, ivig_dose = 1, ivig_start = 12)),
    list(id = 9,  lab = "anti-D 60 + IVIG 2 g/kg/wk from 12 wk",
         args = list(anti_d_iu = 60, ivig_dose = 2, ivig_start = 12)),
    list(id = 10, lab = "anti-D 60 + nipocalimab 30 mg/kg from 14 wk",
         args = list(anti_d_iu = 60, nip_dose = 30, nip_start = 14)),
    list(id = 11, lab = "anti-D 60 + nipocalimab 45 mg/kg from 14 wk",
         args = list(anti_d_iu = 60, nip_dose = 45, nip_start = 14)),
    list(id = 12, lab = "anti-D 120 + nipocalimab 30 mg/kg from 14 wk",
         args = list(anti_d_iu = 120, nip_dose = 30, nip_start = 14)),
    list(id = 13, lab = "anti-D 120 + nipocalimab 30 + IVIG 1 g/kg",
         args = list(anti_d_iu = 120, nip_dose = 30, nip_start = 14,
                     ivig_dose = 1, ivig_start = 14)),
    list(id = 14, lab = "anti-D 120 + nipocalimab from 24 wk (too late)",
         args = list(anti_d_iu = 120, nip_dose = 30, nip_start = 24)),
    list(id = 15, lab = "anti-D 15, deliver at 34 wk",
         args = list(anti_d_iu = 15, ga_deliver = 34)),
    list(id = 16, lab = "anti-D 15, deliver at 37 wk",
         args = list(anti_d_iu = 15, ga_deliver = 37))
  )
  do.call(rbind, lapply(S, function(s) {
    r  <- do.call(simulate_hdfn, c(list(mod = mod), s$args))
    o  <- r$out
    nb <- o[o$time >= (s$args$ga_deliver %||% 36)*7, ]
    data.frame(scenario = s$id, label = s$lab, n_iut = r$n_iut,
               first_iut = ifelse(is.null(r$iut), NA, min(r$iut$ga)),
               hydrops = as.integer(any(o$HYDROPS > 0.5)),
               survival = min(o$SURV),
               hb_birth = ifelse(nrow(nb) > 0, nb$HB[1], NA),
               tsb_peak = ifelse(nrow(nb) > 0, max(nb$TSB), NA),
               n_exch = r$n_exch, n_topup = r$n_topup)
  }))
}

## Sensitisation / prophylaxis regimens are a maternal-side-only problem and are
## driven with the same model by setting anti_d_iu = 0 and dosing Mdep:
##   300 ug  = 1500 IU into the Mdep compartment
##   an FMH of V mL of fetal whole blood = 0.42*V mL into Mfree
## The population risk is the expectation of 1 - exp(-Sens) over the log-normal
## FMH distribution (median 0.20 mL, sigma 1.82), which is what
## hdfn_reference_output.txt section 9 reports:
##      none 15.6% | postpartum 300 ug 1.35% | antenatal + postpartum ~0.2%
## and the 300 ug column turns over at exactly 30 mL of fetal whole blood,
## which is the stated coverage of that dose.

## Run the scenario table only when explicitly asked, so that sourcing this
## file from hdfn_shiny_app.R does not launch sixteen simulations:
##     HDFN_RUN_SCENARIOS=1 Rscript hdfn_mrgsolve_model.R
if (identical(Sys.getenv("HDFN_RUN_SCENARIOS"), "1")) {
  print(hdfn_scenarios(hdfn))
}
