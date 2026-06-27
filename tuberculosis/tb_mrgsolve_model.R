##############################################################################
# Tuberculosis QSP Model — mrgsolve Implementation
# Disease: Pulmonary Tuberculosis (Mtb infection)
# Drug:    RIPE (Rifampicin, Isoniazid, Pyrazinamide, Ethambutol) ±
#          Bedaquiline (BDQ) for DR-TB
#
# Model structure:
#   • PK:  1-compartment (RIF, INH, PZA, EMB, BDQ) with oral absorption
#          RIF auto-induction of CYP3A4; INH NAT2 acetylation polymorphism
#   • Host: Macrophage dynamics (resting→infected→activated)
#           T-cell & cytokine (Th1, IFN-γ, IL-10, TNF-α)
#   • PD:  Three bacterial subpopulations (AR, SR, NR)
#          Emax killing for each drug vs each population
#
# Calibration references:
#   • Gumbo et al. (2007) Antimicrob Agents Chemother 51:3633–3641 (RIF EKmax)
#   • Jayaram et al. (2003) Antimicrob Agents Chemother 47:2118–2124 (INH PKPD)
#   • Srivastava et al. (2011) Antimicrob Agents Chemother 55:1743–1748 (PZA)
#   • Magombedze et al. (2010) Theor Biol Med Model 7:41 (immune ODE)
#   • Wigginton & Kirschner (2001) J Immunol 166:1951–1967 (granuloma model)
#
# Scenarios:
#   1. Natural history (no treatment) — Mtb growth to clearance/death
#   2. Standard RIPE 6-month (2HRZE/4HR) — drug-susceptible TB
#   3. Poor adherence (30 % missed doses) — subtherapeutic exposure
#   4. MDR-TB regimen (BDQ-Pa-LZD substitute) — DR-TB
#   5. HIV co-infected (CD4<200) — impaired Th1 / macrophage function
#   6. Diabetic host (HbA1c 9%) — impaired macrophage activation
##############################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

tb_model_code <- '
$PLUGIN autodiff

$PARAM @annotated
//------------------------------------------------------------------
// Bacterial parameters
//------------------------------------------------------------------
kgAR   : 0.6931 : Active replicating growth rate (/day; doubling ~1d)
kgSR   : 0.0693 : Slow replicating growth rate (/day; doubling ~10d)
kAR2SR : 0.050  : AR→SR transition rate (/day; hypoxia/stress)
kSR2NR : 0.020  : SR→NR transition rate (/day; full dormancy)
kNR2AR : 0.0005 : NR reactivation rate (/day; immunosuppression)
kdAR   : 0.040  : Natural AR death rate (/day)
kdSR   : 0.010  : Natural SR death rate (/day)
kdNR   : 0.002  : Natural NR death rate (/day)
Bmax   : 1.0e9  : Maximum bacterial burden (CFU; ~10^9 at end-stage)

//------------------------------------------------------------------
// Innate / Adaptive immune parameters
//------------------------------------------------------------------
UM_ss  : 500.0  : Steady-state uninfected macrophage count (AU)
kUM_p  : 50.0   : Macrophage production rate (AU/day)
kUM_d  : 0.10   : Macrophage natural death rate (/day)
kInf   : 1.0e-4 : Macrophage infection rate (AR_intra uptake; /CFU/day)
kActv  : 0.50   : Macrophage activation rate by IFNg (/AU/day)
kIM_d  : 0.15   : Infected macrophage death rate (/day; necrosis+apoptosis)
kAM_d  : 0.08   : Activated macrophage death rate (/day)
kAMkill: 2.0    : Activated Mφ bacterial kill rate (/day)

// T-cell and cytokine parameters
kTh1_p : 5.0    : Th1 production rate driven by IL-12/IFN-γ (AU/day)
kTh1_d : 0.10   : Th1 death rate (/day)
kTh1_half: 50.0 : Bacterial load (CFU, half-maximum) for Th1 stim.
kIFNg_p: 2.0    : IFN-γ production per Th1 cell (AU/AU/day)
kIFNg_d: 1.50   : IFN-γ clearance rate (/day)
kTNFa_p: 1.0    : TNF-α production rate by AM_inf (AU/day/AU)
kTNFa_d: 2.0    : TNF-α clearance (/day)
kIL10_p: 0.30   : IL-10 production by Treg/AM (/day)
kIL10_d: 1.0    : IL-10 clearance (/day)
kIL10_sup: 0.5  : IL-10 suppression of IFN-γ (half-max, AU)

// Host modifier (normal = 1.0; reduced in HIV/DM scenarios)
phi_Th1  : 1.0  : Th1 functional scaling factor (HIV: 0.2, normal: 1.0)
phi_Mact : 1.0  : Macrophage activation scaling (DM: 0.4, normal: 1.0)

//------------------------------------------------------------------
// RIF PK: 1-compartment + gut depot
//   Reference: Wilkins et al. (2008) Antimicrob Agents Chemother
//   Typical pop: ka=1.5/h, V=49L/70kg, CL=9L/h/70kg, F=0.68
//------------------------------------------------------------------
RIF_dose_mg : 600.0   : Rifampicin daily dose (mg)
RIF_ka      : 1.50    : RIF absorption rate (/h)
RIF_V       : 49.0    : RIF volume of distribution (L)
RIF_CL_base : 9.0     : RIF clearance at baseline (L/h)
RIF_F       : 0.68    : RIF bioavailability
RIF_ind_max : 0.60    : Maximum auto-induction of RIF CL (fractional)
RIF_ind_ec50: 4.0     : EC50 for auto-induction (mg/L)
RIF_lung_f  : 0.25    : RIF lung/lesion penetration ratio

//------------------------------------------------------------------
// INH PK: 1-compartment, NAT2 acetylation
//   Reference: Kinzig-Schippers et al. (2005) Clin Pharmacokinet
//   Slow: CL~3.5 L/h, Fast: CL~14 L/h
//------------------------------------------------------------------
INH_dose_mg : 300.0   : Isoniazid daily dose (mg)
INH_ka      : 1.20    : INH absorption rate (/h)
INH_V       : 38.0    : INH volume (L)
INH_CL      : 5.0     : INH acetylation CL (L/h; weighted mean slow/fast)
INH_F       : 0.90    : INH bioavailability
NAT2_factor : 1.0     : NAT2 acetylation factor (slow: 0.4, fast: 1.8)

//------------------------------------------------------------------
// PZA PK: 1-compartment
//   Reference: Zhu et al. (2002) Clin Pharmacol Ther
//   Cmax~35 mg/L, t1/2~9h, V=36L, CL=2.8 L/h
//------------------------------------------------------------------
PZA_dose_mg : 1500.0  : Pyrazinamide daily dose (mg)
PZA_ka      : 1.0     : PZA absorption rate (/h)
PZA_V       : 36.0    : PZA volume (L)
PZA_CL      : 2.8     : PZA clearance (L/h)
PZA_F       : 0.98    : PZA bioavailability
PZA_pH_act  : 1.0     : PZA activity scaling (acid pH 5.5 enhancement)

//------------------------------------------------------------------
// EMB PK: 1-compartment
//   Reference: Peloquin et al. (1999) Pharmacotherapy
//   Cmax~2.5 mg/L, t1/2~7h, V=60L, CL=9.4 L/h (renal)
//------------------------------------------------------------------
EMB_dose_mg : 1200.0  : Ethambutol daily dose (mg)
EMB_ka      : 1.1     : EMB absorption rate (/h)
EMB_V       : 60.0    : EMB volume (L)
EMB_CL      : 9.4     : EMB clearance (L/h)
EMB_F       : 0.78    : EMB bioavailability

//------------------------------------------------------------------
// BDQ PK: 1-compartment (for DR-TB scenario)
//   Reference: Svensson et al. (2013) Antimicrob Agents Chemother
//   Loading 400 mg QD x14d → 200 mg TIW; t1/2 ~5.5 months
//------------------------------------------------------------------
BDQ_dose_mg : 400.0   : Bedaquiline loading dose (mg)
BDQ_ka      : 0.50    : BDQ absorption rate (/h)
BDQ_V       : 164000. : BDQ volume (L; very high distribution)
BDQ_CL      : 91.0    : BDQ clearance (L/h)
BDQ_F       : 0.45    : BDQ bioavailability
BDQ_on      : 0.0     : BDQ flag (1=DR-TB regimen, 0=off)

//------------------------------------------------------------------
// PD: Emax killing parameters
//   AR = actively replicating bacteria
//   SR = slowly replicating bacteria
//   NR = non-replicating persister
//------------------------------------------------------------------
// RIF (bactericidal on AR; sterilizing on SR)
RIF_Emax_AR : 3.50    : RIF max kill rate on AR (/day)
RIF_EC50_AR : 0.60    : RIF EC50 on AR (mg/L; ≈MIC for susceptible)
RIF_Emax_SR : 1.50    : RIF max kill rate on SR (/day)
RIF_EC50_SR : 1.50    : RIF EC50 on SR (mg/L)
RIF_hill    : 1.0     : RIF Hill coefficient

// INH (bactericidal on AR; less active on SR)
INH_Emax_AR : 4.00    : INH max kill rate on AR (/day)
INH_EC50_AR : 0.10    : INH EC50 on AR (mg/L; MIC ≤0.2 mg/L)
INH_Emax_SR : 0.50    : INH max kill rate on SR (/day)
INH_EC50_SR : 0.50    : INH EC50 on SR (mg/L)
INH_hill    : 1.0     : INH Hill coefficient

// PZA (sterilizing; active on SR and NR at low pH)
PZA_Emax_SR : 2.00    : PZA max kill rate on SR (/day)
PZA_EC50_SR : 20.0    : PZA EC50 on SR (mg/L; higher MIC)
PZA_Emax_NR : 1.00    : PZA max kill rate on NR (/day)
PZA_EC50_NR : 30.0    : PZA EC50 on NR (mg/L)
PZA_hill    : 1.5     : PZA Hill coefficient

// EMB (bacteriostatic on AR)
EMB_Emax_AR : 0.80    : EMB max inhibition of AR growth (/day)
EMB_EC50_AR : 1.0     : EMB EC50 on AR (mg/L)
EMB_hill    : 1.0     : EMB Hill coefficient

// BDQ (bactericidal/sterilizing; active on SR and NR)
BDQ_Emax_SR : 2.50    : BDQ max kill rate on SR (/day)
BDQ_EC50_SR : 0.06    : BDQ EC50 on SR (mg/L)
BDQ_Emax_NR : 2.00    : BDQ max kill rate on NR (/day)
BDQ_EC50_NR : 0.10    : BDQ EC50 on NR (mg/L)
BDQ_hill    : 1.0     : BDQ Hill coefficient

// Resistance flags (1=resistant, 0=susceptible)
rpoB_mut   : 0.0      : rpoB mutation (1=RIF resistant)
katG_mut   : 0.0      : katG/inhA mutation (1=INH resistant)

//------------------------------------------------------------------
// Treatment duration flags (1=drug present, 0=absent)
//------------------------------------------------------------------
on_RIF : 1.0  : Rifampicin on/off
on_INH : 1.0  : Isoniazid on/off
on_PZA : 1.0  : Pyrazinamide on/off (first 2 months)
on_EMB : 1.0  : Ethambutol on/off (first 2 months)
adhere : 1.0  : Adherence factor (1=full, 0.7=30% missed)

$CMT @annotated
// Drug PK compartments
RIF_gut  : RIF gut depot (mg)
RIF_c    : RIF central plasma (mg)
INH_c    : INH central plasma (mg)
PZA_c    : PZA central plasma (mg)
EMB_c    : EMB central plasma (mg)
BDQ_c    : BDQ central plasma (mg)

// Bacterial subpopulations (log10 CFU/mL)
AR       : Actively replicating bacteria (CFU)
SR       : Slowly replicating bacteria (CFU)
NR       : Non-replicating persisters (CFU)

// Immune compartments
UM       : Uninfected macrophages (AU)
IM       : Infected macrophages (AU)
AM       : Activated macrophages (AU)
Th1      : Th1 CD4+ T cells (AU)
IFNg     : IFN-γ concentration (AU)
TNFa     : TNF-α concentration (AU)
IL10     : IL-10 concentration (AU)

$GLOBAL
double RIF_conc, INH_conc, PZA_conc, EMB_conc, BDQ_conc;
double kRIF_AR, kRIF_SR, kINH_AR, kINH_SR;
double kPZA_SR, kPZA_NR, kEMB_AR_inh, kBDQ_SR, kBDQ_NR;
double RIF_CL_eff;

$MAIN
// Initial conditions — simulate initiation of infection
AR_0  = 10.0;    // initial inoculum (10 CFU)
SR_0  = 0.0;
NR_0  = 0.0;
UM_0  = UM_ss;
IM_0  = 0.0;
AM_0  = 0.0;
Th1_0 = 0.0;
IFNg_0 = 0.001;
TNFa_0 = 0.01;
IL10_0 = 0.01;
RIF_gut_0 = 0.0;
RIF_c_0   = 0.0;
INH_c_0   = 0.0;
PZA_c_0   = 0.0;
EMB_c_0   = 0.0;
BDQ_c_0   = 0.0;

// Effective RIF clearance with auto-induction
RIF_CL_eff = RIF_CL_base * (1.0 + RIF_ind_max * RIF_c/(RIF_ind_ec50*RIF_V + RIF_c));

$ODE
//==========================================================
// RIF PK (absorbed → plasma → clearance)
//==========================================================
double RIF_rate_in = on_RIF * adhere * RIF_dose_mg * RIF_F * RIF_ka;
double RIF_ka_out  = RIF_ka * RIF_gut;
dxdt_RIF_gut = RIF_rate_in - RIF_ka_out;

RIF_CL_eff = RIF_CL_base * (1.0 + RIF_ind_max * (RIF_c/RIF_V) /
             (RIF_ind_ec50 + RIF_c/RIF_V));
dxdt_RIF_c = RIF_ka_out - (RIF_CL_eff / RIF_V) * RIF_c;

//==========================================================
// INH PK (accounting for NAT2 phenotype)
//==========================================================
double INH_rate_in = on_INH * adhere * INH_dose_mg * INH_F * INH_ka;
double INH_CL_eff  = INH_CL * NAT2_factor;
dxdt_INH_c = INH_rate_in - (INH_CL_eff / INH_V) * INH_c;

//==========================================================
// PZA PK
//==========================================================
double PZA_rate_in = on_PZA * adhere * PZA_dose_mg * PZA_F * PZA_ka;
dxdt_PZA_c = PZA_rate_in - (PZA_CL / PZA_V) * PZA_c;

//==========================================================
// EMB PK
//==========================================================
double EMB_rate_in = on_EMB * adhere * EMB_dose_mg * EMB_F * EMB_ka;
dxdt_EMB_c = EMB_rate_in - (EMB_CL / EMB_V) * EMB_c;

//==========================================================
// BDQ PK
//==========================================================
double BDQ_rate_in = BDQ_on * adhere * BDQ_dose_mg * BDQ_F * BDQ_ka;
dxdt_BDQ_c = BDQ_rate_in - (BDQ_CL / BDQ_V) * BDQ_c;

//==========================================================
// Derived drug concentrations (mg/L = mg / L volume)
//==========================================================
RIF_conc = (RIF_c / RIF_V) * RIF_lung_f;  // lung penetration
INH_conc = INH_c / INH_V;
PZA_conc = (PZA_c / PZA_V) * PZA_pH_act;
EMB_conc = EMB_c / EMB_V;
BDQ_conc = BDQ_c / BDQ_V;

//==========================================================
// Drug PD: Emax killing rates
// (set to 0 if resistant or drug absent)
//==========================================================
double HILLrif = pow(RIF_conc, RIF_hill);
double HILLinh = pow(INH_conc, INH_hill);
double HILLpza = pow(PZA_conc, PZA_hill);
double HILLemb = pow(EMB_conc, EMB_hill);
double HILLbdq = pow(BDQ_conc, BDQ_hill);

// RIF effect (abolished if rpoB mutant)
kRIF_AR = (1.0 - rpoB_mut) * RIF_Emax_AR * HILLrif /
          (pow(RIF_EC50_AR, RIF_hill) + HILLrif);
kRIF_SR = (1.0 - rpoB_mut) * RIF_Emax_SR * HILLrif /
          (pow(RIF_EC50_SR, RIF_hill) + HILLrif);

// INH effect (abolished if katG mutant)
kINH_AR = (1.0 - katG_mut) * INH_Emax_AR * HILLinh /
          (pow(INH_EC50_AR, INH_hill) + HILLinh);
kINH_SR = (1.0 - katG_mut) * INH_Emax_SR * HILLinh /
          (pow(INH_EC50_SR, INH_hill) + HILLinh);

// PZA effect
kPZA_SR  = PZA_Emax_SR * HILLpza /
           (pow(PZA_EC50_SR, PZA_hill) + HILLpza);
kPZA_NR  = PZA_Emax_NR * HILLpza /
           (pow(PZA_EC50_NR, PZA_hill) + HILLpza);

// EMB growth inhibition on AR
kEMB_AR_inh = EMB_Emax_AR * HILLemb /
              (pow(EMB_EC50_AR, EMB_hill) + HILLemb);

// BDQ effect (DR-TB sterilizing)
kBDQ_SR = BDQ_on * BDQ_Emax_SR * HILLbdq /
          (pow(BDQ_EC50_SR, BDQ_hill) + HILLbdq);
kBDQ_NR = BDQ_on * BDQ_Emax_NR * HILLbdq /
          (pow(BDQ_EC50_NR, BDQ_hill) + HILLbdq);

//==========================================================
// Immune system ODEs
//==========================================================
double total_B = AR + SR + NR;
double total_B_safe = total_B < 1.0 ? 1.0 : total_B;

// Uninfected macrophage (UM): produced at constant rate, die,
// become infected (IM) when AR bacteria present
double UM_pos = UM < 0.0 ? 0.0 : UM;
double AR_pos = AR < 0.0 ? 0.0 : AR;
double SR_pos = SR < 0.0 ? 0.0 : SR;
double NR_pos = NR < 0.0 ? 0.0 : NR;

dxdt_UM = kUM_p - kUM_d * UM_pos - kInf * UM_pos * AR_pos;

// Infected macrophage (IM): gain from UM infection; killed by CTL/activated Mφ
double IFNg_pos = IFNg < 0.0 ? 0.0 : IFNg;
double AM_pos = AM < 0.0 ? 0.0 : AM;

dxdt_IM = kInf * UM_pos * AR_pos - kIM_d * IM < 0.0 ? 0.0 :
          kInf * UM_pos * AR_pos - kIM_d * (IM < 0.0 ? 0.0 : IM);

// Activated macrophage (AM): IFN-γ drives activation; kills bacteria
double IFNg_act = IFNg_pos / (kIL10_sup + IL10 < 0.0 ? 0.0 : IL10);
dxdt_AM = kActv * phi_Mact * IFNg_act * (IM < 0.0 ? 0.0 : IM) -
          kAM_d * (AM < 0.0 ? 0.0 : AM);

// Th1 cells: stimulated by bacterial load; inhibited by IL-10
double Th1_pos = Th1 < 0.0 ? 0.0 : Th1;
double Th1_stim = total_B_safe / (kTh1_half + total_B_safe);
dxdt_Th1 = phi_Th1 * kTh1_p * Th1_stim - kTh1_d * Th1_pos;

// IFN-γ: produced by Th1; cleared; suppressed by IL-10
dxdt_IFNg = kIFNg_p * Th1_pos - kIFNg_d * IFNg_pos;

// TNF-α: produced by IM; cleared
double IM_pos = IM < 0.0 ? 0.0 : IM;
dxdt_TNFa = kTNFa_p * IM_pos - kTNFa_d * (TNFa < 0.0 ? 0.0 : TNFa);

// IL-10: produced; cleared
dxdt_IL10 = kIL10_p * Th1_pos - kIL10_d * (IL10 < 0.0 ? 0.0 : IL10);

//==========================================================
// Bacterial ODE: AR (Actively Replicating)
//==========================================================
// Growth: logistic; Death: natural + drug killing + AM killing
// AR → SR: stress/hypoxia
double AR_growth = kgAR * AR_pos * (1.0 - total_B_safe / Bmax) -
                   kEMB_AR_inh * AR_pos;
double AR_kill   = (kRIF_AR + kINH_AR + kdAR) * AR_pos +
                   kAMkill * AM_pos * AR_pos / (100.0 + AR_pos);
dxdt_AR = AR_growth - AR_kill - kAR2SR * AR_pos;

//==========================================================
// Bacterial ODE: SR (Slowly Replicating)
//==========================================================
double SR_growth = kgSR * SR_pos * (1.0 - total_B_safe / Bmax);
double SR_kill   = (kRIF_SR + kINH_SR + kPZA_SR + kBDQ_SR + kdSR) * SR_pos;
dxdt_SR = SR_growth - SR_kill + kAR2SR * AR_pos - kSR2NR * SR_pos;

//==========================================================
// Bacterial ODE: NR (Non-Replicating Persisters)
//==========================================================
double NR_kill = (kPZA_NR + kBDQ_NR + kdNR) * NR_pos;
dxdt_NR = kSR2NR * SR_pos - NR_kill - kNR2AR * NR_pos;

$TABLE
// Derived quantities for output
double total_bact   = AR + SR + NR;
double log10_AR     = log10(AR  + 0.001);
double log10_SR     = log10(SR  + 0.001);
double log10_NR     = log10(NR  + 0.001);
double log10_total  = log10(total_bact + 0.001);
double RIF_Cplasma  = RIF_c / RIF_V;
double INH_Cplasma  = INH_c / INH_V;
double PZA_Cplasma  = PZA_c / PZA_V;
double EMB_Cplasma  = EMB_c / EMB_V;
double BDQ_Cplasma  = BDQ_c / BDQ_V;
// Culture conversion: total bacteria < 10 CFU/mL
double cult_conv    = (total_bact < 10.0) ? 1.0 : 0.0;
// Sputum smear: AR > 1e4 CFU/mL (roughly)
double smear_pos    = (AR > 1.0e4) ? 1.0 : 0.0;
double immune_index = AM / (UM + 0.1);   // ratio activated/resting

$CAPTURE
total_bact log10_total log10_AR log10_SR log10_NR
RIF_Cplasma INH_Cplasma PZA_Cplasma EMB_Cplasma BDQ_Cplasma
UM IM AM Th1 IFNg TNFa IL10
cult_conv smear_pos immune_index
'

#------------------------------------------------------------------------------
# Compile the model
#------------------------------------------------------------------------------
mod <- mcode("tb_qsp", tb_model_code)

#==============================================================================
# SCENARIO DEFINITIONS
#==============================================================================
# Helper: build once-daily dosing event table (days 1-180 for PZA/EMB, 1-180 for RIF/INH)
make_events <- function(start_day = 0, end_rife = 60, end_ri = 180,
                        rif_mg = 600, inh_mg = 300, pza_mg = 1500, emb_mg = 1200,
                        dose_hrs = 0,   # time of day (hours)
                        adhere_prob = 1.0, seed_val = 42) {
  set.seed(seed_val)
  days_ri   <- seq(start_day, end_ri   - 1)
  days_rife <- seq(start_day, end_rife - 1)
  # missed doses flag
  dose_RIF <- sapply(days_ri,   function(d) rbinom(1, 1, adhere_prob) * rif_mg)
  dose_INH <- sapply(days_ri,   function(d) rbinom(1, 1, adhere_prob) * inh_mg)
  dose_PZA <- sapply(days_rife, function(d) rbinom(1, 1, adhere_prob) * pza_mg)
  dose_EMB <- sapply(days_rife, function(d) rbinom(1, 1, adhere_prob) * emb_mg)

  ev_RIF <- ev(time = days_ri   * 24 + dose_hrs, amt = dose_RIF, cmt = "RIF_gut", addl = 0)
  ev_INH <- ev(time = days_ri   * 24 + dose_hrs, amt = dose_INH, cmt = "INH_c",   addl = 0)
  ev_PZA <- ev(time = days_rife * 24 + dose_hrs, amt = dose_PZA, cmt = "PZA_c",   addl = 0)
  ev_EMB <- ev(time = days_rife * 24 + dose_hrs, amt = dose_EMB, cmt = "EMB_c",   addl = 0)
  as_data_frame(ev_RIF) %>%
    bind_rows(as_data_frame(ev_INH)) %>%
    bind_rows(as_data_frame(ev_PZA)) %>%
    bind_rows(as_data_frame(ev_EMB)) %>%
    arrange(time)
}

#------------------------------------------------------------------------------
# Scenario 1: Natural History (no treatment)
#------------------------------------------------------------------------------
run_nat <- function(end_hrs = 360*24) {
  mod %>%
    param(on_RIF=0, on_INH=0, on_PZA=0, on_EMB=0, adhere=1) %>%
    mrgsim(end = end_hrs, delta = 24, hmax = 6) %>%
    as_tibble() %>%
    mutate(scenario = "1. Natural History")
}

#------------------------------------------------------------------------------
# Scenario 2: Standard RIPE (2HRZE / 4HR) — full adherence
#------------------------------------------------------------------------------
run_ripe <- function() {
  ev_data <- make_events(end_rife = 60, end_ri = 180, adhere_prob = 1.0)
  mod %>%
    param(on_RIF=1, on_INH=1, on_PZA=1, on_EMB=1,
          adhere=1, phi_Th1=1, phi_Mact=1) %>%
    data_set(ev_data) %>%
    mrgsim(end = 210*24, delta = 24, hmax = 6) %>%
    as_tibble() %>%
    mutate(scenario = "2. Standard RIPE (Full Adherence)")
}

#------------------------------------------------------------------------------
# Scenario 3: Poor Adherence (30% missed doses)
#------------------------------------------------------------------------------
run_poor_adhere <- function() {
  ev_data <- make_events(end_rife = 60, end_ri = 180, adhere_prob = 0.70)
  mod %>%
    param(on_RIF=1, on_INH=1, on_PZA=1, on_EMB=1,
          adhere=1, phi_Th1=1, phi_Mact=1) %>%
    data_set(ev_data) %>%
    mrgsim(end = 210*24, delta = 24, hmax = 6) %>%
    as_tibble() %>%
    mutate(scenario = "3. Poor Adherence (30% Missed)")
}

#------------------------------------------------------------------------------
# Scenario 4: MDR-TB → BDQ-based regimen (RIF + INH resistant)
#------------------------------------------------------------------------------
run_mdr <- function() {
  # BDQ loading 400 mg x14d then 200 mg x3/wk for 24 weeks
  ev_bdq_load <- ev(time  = seq(0,  13) * 24, amt = 400, cmt = "BDQ_c")
  ev_bdq_maint <- ev(time = seq(14, 168, by = 7/3) * 24, amt = 200, cmt = "BDQ_c")
  ev_bdq <- bind_rows(as_data_frame(ev_bdq_load), as_data_frame(ev_bdq_maint))

  mod %>%
    param(on_RIF=0, on_INH=0, on_PZA=1, on_EMB=0,
          BDQ_on=1, adhere=1,
          rpoB_mut=1, katG_mut=1,   # MDR-TB: RIF+INH resistant
          phi_Th1=1, phi_Mact=1) %>%
    data_set(ev_bdq) %>%
    mrgsim(end = 180*24, delta = 24, hmax = 6) %>%
    as_tibble() %>%
    mutate(scenario = "4. MDR-TB → BDQ-Based Regimen")
}

#------------------------------------------------------------------------------
# Scenario 5: HIV Co-Infected (CD4 < 200; impaired Th1 phi_Th1=0.2)
#------------------------------------------------------------------------------
run_hiv <- function() {
  ev_data <- make_events(end_rife = 60, end_ri = 180, adhere_prob = 0.90)
  mod %>%
    param(on_RIF=1, on_INH=1, on_PZA=1, on_EMB=1,
          adhere=1,
          phi_Th1 = 0.20,   # severe Th1 impairment
          phi_Mact = 0.60) %>%
    data_set(ev_data) %>%
    mrgsim(end = 210*24, delta = 24, hmax = 6) %>%
    as_tibble() %>%
    mutate(scenario = "5. HIV Co-Infected (CD4<200)")
}

#------------------------------------------------------------------------------
# Scenario 6: Diabetic Host (HbA1c ~9%; impaired macrophage activation)
#------------------------------------------------------------------------------
run_dm <- function() {
  ev_data <- make_events(end_rife = 60, end_ri = 180, adhere_prob = 1.0)
  mod %>%
    param(on_RIF=1, on_INH=1, on_PZA=1, on_EMB=1,
          adhere=1,
          phi_Th1 = 0.80,   # mild Th1 impairment
          phi_Mact = 0.40) %>%  # significant macrophage activation impairment
    data_set(ev_data) %>%
    mrgsim(end = 210*24, delta = 24, hmax = 6) %>%
    as_tibble() %>%
    mutate(scenario = "6. Diabetic Host (Impaired Mφ)")
}

#==============================================================================
# RUN ALL SCENARIOS
#==============================================================================
cat("Running Scenario 1: Natural History...\n")
out1 <- run_nat()
cat("Running Scenario 2: Standard RIPE...\n")
out2 <- run_ripe()
cat("Running Scenario 3: Poor Adherence...\n")
out3 <- run_poor_adhere()
cat("Running Scenario 4: MDR-TB...\n")
out4 <- run_mdr()
cat("Running Scenario 5: HIV Co-infected...\n")
out5 <- run_hiv()
cat("Running Scenario 6: Diabetic...\n")
out6 <- run_dm()

all_out <- bind_rows(out1, out2, out3, out4, out5, out6) %>%
  mutate(day = time / 24)

#==============================================================================
# PLOTTING
#==============================================================================
theme_tb <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#E8F5E9"),
        panel.grid.minor = element_blank())

sc_cols <- c(
  "1. Natural History"             = "#B71C1C",
  "2. Standard RIPE (Full Adherence)" = "#1B5E20",
  "3. Poor Adherence (30% Missed)" = "#F57F17",
  "4. MDR-TB → BDQ-Based Regimen"  = "#4A148C",
  "5. HIV Co-Infected (CD4<200)"   = "#880E4F",
  "6. Diabetic Host (Impaired Mφ)" = "#1565C0"
)

# --- Plot 1: Total bacterial burden over time ---------------------------------
p1 <- all_out %>%
  filter(day <= 210) %>%
  ggplot(aes(day, log10_total, color = scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  annotate("text", x=5, y=1.3, label="Culture conversion threshold", size=3) +
  scale_color_manual(values = sc_cols) +
  labs(
    title    = "Tuberculosis QSP Model — Total Bacterial Burden",
    subtitle = "RIPE standard therapy vs. impaired host / MDR-TB scenarios",
    x = "Day", y = "log₁₀ Total Bacterial Burden (CFU/mL)",
    color = "Scenario"
  ) +
  theme_tb +
  guides(color = guide_legend(ncol = 2))

# --- Plot 2: Bacterial subpopulations (Scenario 2 only) ----------------------
p2 <- out2 %>%
  filter(day <= 210) %>%
  select(day, log10_AR, log10_SR, log10_NR) %>%
  pivot_longer(-day, names_to = "pop", values_to = "log10_cfu") %>%
  mutate(pop = recode(pop,
    log10_AR = "Active Replicating (AR)",
    log10_SR = "Slow Replicating (SR)",
    log10_NR = "Non-Replicating (NR)")) %>%
  ggplot(aes(day, log10_cfu, color = pop, linetype = pop)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = c(0, 60), linetype = "dotted", color = "grey40") +
  annotate("text", x=30, y=9.5, label="Intensive Phase\n(2HRZE)", size=3.5) +
  annotate("text", x=120, y=9.5, label="Continuation\n(4HR)", size=3.5) +
  scale_color_manual(values = c("#E65100","#1565C0","#6A1B9A")) +
  labs(
    title    = "Bacterial Subpopulation Dynamics — Standard RIPE",
    x = "Day", y = "log₁₀ Bacterial Subpopulation (CFU/mL)",
    color = "Population", linetype = "Population"
  ) +
  theme_tb

# --- Plot 3: Drug PK profiles (Scenario 2, first 48h) ------------------------
p3 <- out2 %>%
  filter(day <= 14) %>%
  select(day, RIF_Cplasma, INH_Cplasma, PZA_Cplasma, EMB_Cplasma) %>%
  pivot_longer(-day, names_to = "drug", values_to = "conc_mgL") %>%
  mutate(drug = recode(drug,
    RIF_Cplasma = "Rifampicin (target >2 mg/L)",
    INH_Cplasma = "Isoniazid (target >3 mg/L)",
    PZA_Cplasma = "Pyrazinamide (target >20 mg/L)",
    EMB_Cplasma = "Ethambutol (target >2 mg/L)")) %>%
  ggplot(aes(day, conc_mgL, color = drug)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = c("#D32F2F","#1565C0","#2E7D32","#6D4C41")) +
  facet_wrap(~drug, scales = "free_y") +
  labs(
    title = "RIPE Drug PK Profiles (First 14 Days)",
    x = "Day", y = "Plasma Concentration (mg/L)",
    color = "Drug"
  ) +
  theme_tb +
  theme(legend.position = "none")

# --- Plot 4: Immune dynamics (Scenario 2) ------------------------------------
p4 <- out2 %>%
  filter(day <= 210) %>%
  select(day, UM, IM, AM, Th1, IFNg, TNFa) %>%
  pivot_longer(-day, names_to = "variable", values_to = "value") %>%
  ggplot(aes(day, value, color = variable)) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~variable, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c(
    UM="#4CAF50", IM="#F44336", AM="#FF9800",
    Th1="#9C27B0", IFNg="#2196F3", TNFa="#FF5722")) +
  labs(
    title = "Immune Compartment Dynamics — Standard RIPE",
    x = "Day", y = "Compartment Value (AU)",
    color = "Variable"
  ) +
  theme_tb +
  theme(legend.position = "none")

# --- Plot 5: Culture conversion comparison ------------------------------------
p5 <- all_out %>%
  filter(day <= 210) %>%
  ggplot(aes(day, cult_conv, color = scenario)) +
  geom_step(linewidth = 1.0) +
  scale_color_manual(values = sc_cols) +
  scale_y_continuous(labels = scales::percent, limits = c(-0.05, 1.1)) +
  labs(
    title    = "Sputum Culture Conversion Over Time",
    subtitle = "1 = culture negative; 0 = culture positive",
    x = "Day", y = "Culture Conversion (binary)",
    color = "Scenario"
  ) +
  theme_tb +
  guides(color = guide_legend(ncol = 2))

# Print all plots
print(p1); print(p2); print(p3); print(p4); print(p5)

#==============================================================================
# SUMMARY TABLE
#==============================================================================
summary_tbl <- all_out %>%
  group_by(scenario) %>%
  summarise(
    peak_log10_bact = max(log10_total, na.rm = TRUE),
    day_cult_conv   = suppressWarnings(
      min(day[cult_conv == 1 & !is.na(cult_conv)], na.rm = TRUE)),
    min_log10_bact  = min(log10_total, na.rm = TRUE),
    final_log10_bact = last(log10_total),
    .groups = "drop"
  ) %>%
  mutate(
    day_cult_conv = ifelse(is.infinite(day_cult_conv), NA_real_, day_cult_conv)
  )

cat("\n========== TB QSP Model — Scenario Summary ==========\n")
print(summary_tbl, n = Inf)
cat("=====================================================\n")
