## ---------------------------------------------------------------------------
## icic_mrgsolve_model.R
## ===========================================================================
## IMMUNE CHECKPOINT INHIBITOR-INDUCED COLITIS (ICI colitis / immune-mediated
## diarrhoea and colitis, IMDC) — a 45-state QSP model for mrgsolve.
##
## This file is a line-by-line translation of icic_reference_model.py, which is
## the numerical reference implementation.  Every parameter value here appears
## there, the baseline steady state is exact in both, and the two are compared
## in icic_cross_validation.txt.
##
## ---------------------------------------------------------------------------
## THE ORGANISING IDEA
## ---------------------------------------------------------------------------
## Write the colonic inflammatory drive as a RATIO, not a level:
##
##      Phi_col = (a_eff * Teff + a_trm * Trm) / (a_reg * Treg + kappa)
##
## The two checkpoint drugs enter that fraction in structurally different
## places, and that single fact explains the two things the clinical data
## insist on and a "more immune activation" story cannot:
##
##  NUMERATOR — anti-PD-1 multiplies per-cell activity by an occupancy term
##    O = C_t/(C_t + K_D).  K_D ~ 1 nM, tissue concentrations 2-200 nM, so O
##    runs 0.854 -> 0.995 across a THIRTYFOLD dose range, and the numerator
##    MULTIPLIER only 1.68 -> 1.80 (6.7%).  Saturated at the
##    lowest clinical dose => NO dose-response, which is what nivolumab
##    0.1-10 mg/kg and pembrolizumab 2 vs 10 mg/kg actually show.  PD-1
##    blockade is also EFFERENT: it unleashes clones that already exist, so it
##    enters expansion only weakly (pd1_prol_pow = 0.3) and does not raise the
##    licensed clone count at all.
##
##  DENOMINATOR — anti-CTLA-4 depletes CTLA-4-high colonic Tregs by
##    FcgammaRIIIa-dependent ADCC.  Fc-mediated killing needs immune-complex
##    DENSITY, not mere occupancy, so EC50_ADCC is 10-100x the binding K_D and
##    tissue free drug is ~13% of plasma: EC50_ADCC lands in the MIDDLE of the
##    clinical range.  A_adcc = 0.23 / 0.47 / 0.75 at 1 / 3 / 10 mg/kg.
##    The release factor G = Reg0/Reg is HYPERBOLIC, so it keeps climbing after
##    every occupancy term has flattened.  What caps it is kappa, the non-Treg
##    regulatory floor — the most influential and least measurable parameter
##    in the model (see the sensitivity analysis).
##
## ---------------------------------------------------------------------------
## THE CENSORING RESERVE (why grade is a bad instrument)
## ---------------------------------------------------------------------------
##      Stool = max(Stool_min, L + Sec - A_max * S_eff)
##      S*    = L / A_max = 1750 / 4500 = 0.389
## 61% of colonic absorptive function can be destroyed with a normal stool
## chart.  Calprotectin and endoscopy read the lesion and are NOT censored.
## Symptom-triggered treatment is late by construction.
##
## A caveat stated up front.  Threshold censoring was ALSO meant to explain why
## ipi+nivo colitis is about the same as ipi monotherapy (CheckMate 067: 7.7%
## vs 8.7%).  It does not: in this parameterisation Ppd1 and G multiply, and
## the virtual population gives a clearly SUPER-additive combination.  That is
## a documented failure of the model, not a result of it — see analysis K in
## icic_reference_output.txt.  The structural claims below do not depend on it.
##
## ---------------------------------------------------------------------------
## CASCADE DEPTH IS ONSET TIME (the model is not told these numbers)
## ---------------------------------------------------------------------------
##   infliximab   neutralises TNF-alpha, the terminal executor      -> 2-3 d
##   corticosteroid  suppresses the transcriptional programme       -> 4-6 d
##   vedolizumab  blocks alpha4beta7 TRAFFICKING and cannot touch a
##                cell already resident in the mucosa               -> 14-21 d
##
## ---------------------------------------------------------------------------
## SELECTIVITY IS THE CURRENCY
## ---------------------------------------------------------------------------
## The clones that damage the colon are the clones that kill the tumour, so
## rescue costs tumour control in proportion to non-gut-selectivity chi:
## prednisolone 1.00, infliximab 0.60, vedolizumab 0.10.
##
## ---------------------------------------------------------------------------
## THE DISEASE EATS ITS OWN ANTIDOTE
## ---------------------------------------------------------------------------
## CL_ifx = 0.40 * (ALB/4.0)^-0.9 * (1 + 0.006*CRP), and severe colitis causes
## a protein-losing enteropathy.  A fixed 5 mg/kg dose therefore delivers the
## LOWEST exposure to the SICKEST patient.
##
## ---------------------------------------------------------------------------
## A NOTE ON THE MUCOSAL TOLERANCE BAND (dead_*)
## ---------------------------------------------------------------------------
## Every amplifying loop (IFN-g -> CXCL10 -> macrophage -> TNF -> MAdCAM ->
## homing -> more effectors; and injury -> IL-15 -> Trm self-renewal -> more
## injury) is engaged only once its input leaves a tolerance band.  Without
## these deadbands the resting mucosa has open-loop gain > 1, the healthy state
## is not a steady state, and colitis becomes the inevitable fate of any
## fluctuation.  With them the healthy state is stable and ICI colitis is a
## threshold event — which is what it clinically is.
##
## ---------------------------------------------------------------------------
## USAGE
## ---------------------------------------------------------------------------
##   library(mrgsolve); library(dplyr); library(ggplot2)
##   mod <- mread("icic", "icic_mrgsolve_model.R")
##   # 1. verify the drug-free steady state
##   mod %>% mrgsim(end = 365, delta = 5) %>% as_tibble() %>% tail(1)
##   # 2. the representative susceptible patient on ipilimumab 3 mg/kg
##   mod %>% param(CASE_PARAMS) %>% init(CASE_INIT) %>%
##       ev(amt = 210, cmt = 1, ii = 21, addl = 3) %>%
##       mrgsim(end = 210, delta = 0.5) %>% plot(grade + Calpro + S_eff ~ time)
##
## Dosing compartments (amt in mg):
##   1 A_ipi1  ipilimumab IV        3 A_pd11  anti-PD-1 IV
##   5 A_ifx1  infliximab IV        7 A_vdz1  vedolizumab IV
##   9 A_toc1  tocilizumab IV      11 A_predg prednisolone ORAL
##  14 A_jak   JAK inhibitor ORAL
## Prednisolone is normally given as a continuous mg/day rate through
## PRED_RATE (set by a taper schedule) rather than as discrete oral doses.
## ---------------------------------------------------------------------------

[ PROB ]
# ICI colitis: a ratio with a saturated numerator and a depletable denominator

[ SET ] end = 210, delta = 0.5, atol = 1e-10, rtol = 1e-8

[ CMT ] @annotated
A_ipi1  : ipilimumab central (mg)
A_ipi2  : ipilimumab peripheral (mg)
A_pd11  : anti-PD-1 central (mg)
A_pd12  : anti-PD-1 peripheral (mg)
A_ifx1  : infliximab central (mg)
A_ifx2  : infliximab peripheral (mg)
A_vdz1  : vedolizumab central (mg)
A_vdz2  : vedolizumab peripheral (mg)
A_toc1  : tocilizumab central (mg)
A_toc2  : tocilizumab peripheral (mg)
A_predg : prednisolone gut depot (mg)
A_pred1 : prednisolone central (mg)
Ce_pred : prednisolone effect site (mg/L)
A_jak   : JAK inhibitor central (mg)
Tn      : naive/precursor clone pool (rel)
Nclone  : licensed clone count (priming integral)
Teff    : colonic pathogenic effectors
Trm     : colonic tissue-resident memory (hysteretic state)
Treg    : colonic regulatory capacity (depletable denominator)
Th17    : colonic Th17
Mac     : M1 macrophage (TNF source and ADCC effector)
Neut    : colonic neutrophils
IFNg    : interferon-gamma (rel)
TNFa    : TNF-alpha (rel)
IL6     : interleukin-6 (rel)
IL23    : interleukin-23 (rel)
IL15    : epithelial IL-15 (rel) — the Trm fuel
CXCL10  : CXCL10 (rel)
MADCAM  : MAdCAM-1 (rel)
IL10    : interleukin-10 (rel)
ISC     : LGR5+ stem/crypt units (rel)
Ent     : colonocyte absorptive mass S (rel)
TJ      : tight-junction integrity (rel)
Muc     : mucus / goblet mass (rel)
Dv      : microbial alpha-diversity (rel)
Bu      : utilisable butyrate (rel)
LPS     : translocated microbial product
Alb     : serum albumin (g/dL)
CRP     : C-reactive protein (mg/L)
Calpro  : faecal calprotectin (ug/g)
Ulcer   : deep ulceration index
Tumor   : tumour burden
TumTeff : intratumoural effectors
SterCum : cumulative prednisone-equivalent (mg)
InfHaz  : cumulative opportunistic infection hazard

[ PARAM ] @annotated
// ---- monoclonal antibody PK (population estimates, 70 kg) --------------
CL_ipi : 0.367 : ipilimumab CL (L/d) — 15.3 mL/h, Feng 2013
V1_ipi : 4.4 : ipilimumab central volume (L)
Q_ipi : 0.6 : ipilimumab intercompartmental clearance (L/d)
V2_ipi : 2.8 : ipilimumab peripheral volume (L)
MW_ipi : 148 : ipilimumab molecular weight (kDa)
CL_pd1 : 0.228 : anti-PD-1 CL (L/d) — nivolumab 9.5 mL/h
V1_pd1 : 8 : anti-PD-1 central volume (L)
Q_pd1 : 0.6 : anti-PD-1 intercompartmental clearance (L/d)
V2_pd1 : 4 : anti-PD-1 peripheral volume (L)
MW_pd1 : 146 : anti-PD-1 molecular weight (kDa)
MW_ifx : 149 : infliximab molecular weight (kDa)
MW_vdz : 147 : vedolizumab molecular weight (kDa)
MW_toc : 148 : tocilizumab molecular weight (kDa)
CL_ifx : 0.4 : infliximab reference CL (L/d)
V1_ifx : 5.5 : infliximab central volume (L)
Q_ifx : 0.3 : infliximab intercompartmental clearance (L/d)
V2_ifx : 2 : infliximab peripheral volume (L)
alb_exp_ifx : -0.9 : albumin exponent on infliximab CL
crp_slope_ifx : 0.006 : CRP slope on infliximab CL (per mg/L)
ALB_REF : 4 : reference albumin for CL covariate (g/dL)
CL_vdz : 0.157 : vedolizumab CL (L/d), Rosario 2015
V1_vdz : 4.2 : vedolizumab central volume (L)
Q_vdz : 0.12 : vedolizumab intercompartmental clearance (L/d)
V2_vdz : 1.7 : vedolizumab peripheral volume (L)
alb_exp_vdz : -0.9 : albumin exponent on vedolizumab CL
CL_toc : 0.23 : tocilizumab CL (L/d), linear approximation
V1_toc : 3.5 : tocilizumab central volume (L)
Q_toc : 0.25 : tocilizumab intercompartmental clearance (L/d)
V2_toc : 2 : tocilizumab peripheral volume (L)
// ---- prednisolone and JAK inhibitor -----------------------------------
F_pred : 0.8 : prednisolone oral bioavailability
ka_pred : 18 : prednisolone absorption rate (1/d)
CL_pred : 200 : prednisolone CL (L/d) — t1/2 2.9 h
V_pred : 35 : prednisolone volume (L)
keo_pred : 1.4 : prednisolone effect-site rate (1/d) — genomic effect
EC50_pred : 0.1 : prednisolone effect-site EC50 (mg/L)
Emax_pred : 0.92 : prednisolone maximal effect
ke_jak : 5.2 : JAK inhibitor elimination (1/d) — t1/2 3.2 h
V_jak : 87 : JAK inhibitor apparent volume (L)
EC50_jak : 0.025 : JAK inhibitor EC50 (mg/L)
Emax_jak : 0.85 : JAK inhibitor maximal effect
// ---- tissue penetration and target engagement -------------------------
BDC : 0.13 : plasma->colon biodistribution coefficient
BDC_infl : 0.6 : inflammation-driven increase in penetration
mac_pow : 0.35 : macrophage exponent on ADCC (self-amplification)
KD_PD1_nM : 1 : PD-1 binding KD (nM)
KD_CTLA4_nM : 5.25 : CTLA-4 binding KD (nM), ipilimumab
EC50_ADCC : 4 : ADCC EC50 (tissue ug/mL) — 10-100x KD by design
EC50_IFX : 0.1 : infliximab TNF-neutralisation EC50 (tissue ug/mL)
Emax_IFX : 0.93 : infliximab maximal TNF neutralisation
KD_VDZ : 0.1 : vedolizumab alpha4beta7 KD (tissue ug/mL)
Emax_VDZ : 0.92 : vedolizumab maximal homing blockade
EC50_TOC : 0.3 : tocilizumab EC50 (tissue ug/mL)
Emax_TOC : 0.9 : tocilizumab maximal IL-6 blockade
// ---- THE RATIO --------------------------------------------------------
a_reg : 1 : per-Treg regulatory weight
kappa : 0.15 : NON-Treg regulatory floor — caps the hyperbola G
a_trm : 1.6 : per-cell potency of Trm relative to Teff
Emax_pd1 : 0.8 : efferent-phase per-cell release by PD-1 blockade
Emax_ctla : 2.5 : afferent-phase clone licensing by CTLA-4 blockade
pd1_prol_pow : 0.3 : weak PD-1 effect on expansion (efferent, not licensing)
// ---- priming / clone accounting ---------------------------------------
kp_prime : 0.0015 : baseline priming rate
kout_ncl : 0.03 : licensed-clone loss (1/d) — 33 d memory
lam_tn : 0.05 : naive pool replenishment (1/d)
c_tn : 0.5 : naive pool consumption per unit excess priming
ag_scale : 1 : antigen / susceptibility scale (patient covariate)
a_ff : 0.6 : cap on the barrier->antigen feed-forward
K_ff : 1 : half-constant of the feed-forward
// ---- colonic effector / regulatory dynamics ---------------------------
k_hom : 0.3 : effector homing rate into colon
k_prol : 0.2 : effector proliferation rate (1/d)
Teff_max : 1.2 : lamina propria carrying capacity
k_out_eff : 0.473167 : effector egress (1/d) — SOLVED for steady state
k_conv : 0.006 : Teff -> Trm residency conversion (1/d)
k_dtrm : 0.015 : Trm loss (1/d) — t1/2 46 d
k_sr : 0.05 : Trm self-renewal rate (1/d)
K15 : 1.5 : IL-15 half-constant for Trm self-renewal
Trm_max : 0.35 : Trm carrying capacity
dead_il15 : 1.9 : IL-15 THRESHOLD for the residency latch
k_in_treg : 0.06 : Treg supply (1/d)
k_out_treg : 0.06 : Treg loss (1/d)
b_but : 0.2 : butyrate dependence of pTreg supply
k_adcc : 0.085 : ADCC depletion rate constant (1/d)
phi_FCGR : 1 : FCGR3A genotype factor (V/V 1.60, V/F 1.00, F/F 0.55)
kd_th17 : 0.15 : Th17 turnover (1/d)
E_th17 : 2.5 : maximal Th17 induction
K_th17 : 1.5 : Th17 induction half-constant
kd_mac : 0.2 : macrophage turnover (1/d)
m_cx : 0.9 : CXCL10 -> macrophage gain
m_lps : 0.7 : LPS -> macrophage gain
K_cx_mac : 1 : macrophage induction half-constant
K_lps : 2 : LPS half-constant
kd_neut : 1.2 : neutrophil turnover (1/d)
E_n_th17 : 0.8 : Th17 -> neutrophil gain
E_n_tnf : 0.8 : TNF -> neutrophil gain
K_neut : 1.5 : neutrophil induction half-constant
// ---- mucosal tolerance band -------------------------------------------
dead_cx : 0.2 : CXCL10 tolerance band
dead_tnf : 0.25 : TNF tolerance band
dead_inj : 0.15 : injury tolerance band
dead_th17 : 0.2 : Th17 tolerance band
dead_ifn : 0.1 : IFN-gamma tolerance band
dead_tj : 0.15 : tight-junction tolerance band
dead_neut : 0.1 : neutrophil tolerance band
hinge_s : 0.02 : C1 smoothing width of the hinge function
g_amp : 1 : global recruitment-amplifier gain (sets the separatrix)
// ---- cytokines (saturating targets) -----------------------------------
kd_ifng : 6 : IFN-gamma turnover (1/d)
kd_tnf : 8 : TNF-alpha turnover (1/d)
kd_il6 : 6 : IL-6 turnover (1/d)
kd_il23 : 4 : IL-23 turnover (1/d)
kd_il15 : 0.8 : IL-15 turnover (1/d)
kd_cxcl : 3 : CXCL10 turnover (1/d)
kd_mad : 1.5 : MAdCAM-1 turnover (1/d)
kd_il10 : 4 : IL-10 turnover (1/d)
K_cyt : 90 : drive->cytokine transduction half-constant (excess drive)
E_ifn : 2.2 : maximal IFN-gamma induction
E_tnf_d : 1.5 : drive -> TNF gain
E_tnf_m : 1.2 : macrophage -> TNF gain
E_tnf_l : 1.5 : LPS -> TNF gain
E_il6_d : 1.5 : drive -> IL-6 gain
E_il6_m : 1.2 : macrophage -> IL-6 gain
E_il23 : 2 : macrophage -> IL-23 gain
b_il15 : 3 : injury -> IL-15 gain
K_inj15 : 1 : injury half-constant for IL-15
b_cxcl : 3 : IFN-gamma -> CXCL10 gain
K_ifn_cx : 4 : CXCL10 induction half-constant
b_mad : 1.2 : TNF -> MAdCAM-1 gain
K_tnf_mad : 3 : MAdCAM-1 induction half-constant
mad_max : 2.5 : MAdCAM-1 ceiling
c_hom : 0.5 : CXCL10 -> homing gain
// ---- epithelium -------------------------------------------------------
k_isc : 0.035 : stem compartment recovery (1/d)
k_isck : 0.3 : crypt-base killing rate (1/d)
inj_deep : 1.2 : injury threshold for crypt-base killing
k_prod_e : 0.25 : colonocyte production (1/d)
k_turn : 0.25 : colonocyte turnover (1/d) — 3-5 d lifespan
k_inj : 0.11 : injury-driven colonocyte loss (1/d per unit Ix)
Ix_max : 4 : ceiling of the saturating injury variable
K_ix : 2 : half-constant of the injury saturation
k_tj : 0.5 : tight-junction repair (1/d)
k_tjd : 0.35 : tight-junction damage (1/d)
k_muc : 0.2 : mucus turnover (1/d)
k_mucd : 0.15 : mucus damage (1/d)
w_cyt : 0.55 : cytokine weight in the injury index
w_gzmb : 0.3 : cytotoxic weight in the injury index
w_neut : 0.15 : neutrophil weight in the injury index
E_gzmb : 1.5 : maximal cytotoxic contribution
f_nsaid : 1 : NSAID factor on epithelial restitution (0.85 if exposed)
// ---- lumen ------------------------------------------------------------
k_dv : 0.03 : diversity recovery (1/d)
k_abx : 0.5 : antibiotic effect on diversity (1/d)
k_dys : 0.01 : inflammation-driven dysbiosis (1/d)
k_bu : 0.35 : butyrate turnover (1/d)
k_buc : 0.1 : loss of butyrate utilisation with injury (1/d)
LPS0 : 0.02 : baseline translocated product
k_lps : 2 : LPS equilibration (1/d)
b_lps : 3 : barrier-loss gain on translocation
// ---- water and stool --------------------------------------------------
L_pres : 1750 : ileal effluent presented to colon (mL/d)
A_max : 4500 : maximal colonic absorptive capacity (mL/d)
Stool_min : 100 : normal stool water (mL/d)
mL_per_stool : 155 : incremental stool volume per stool (mL)
freq0 : 1.2 : baseline stool frequency (per day)
f_tj_w : 0.15 : tight-junction weight in absorptive function
f_tr_w : 0.15 : transporter weight in absorptive function
nhe_tnf : 0.9 : TNF effect on NHE3
nhe_ifn : 0.6 : IFN-gamma effect on NHE3
sec_lps : 180 : LPS-driven secretion (mL/d)
sec_bile : 250 : bile-acid-driven secretion (mL/d)
// ---- biomarkers and complications -------------------------------------
Calpro0 : 30 : baseline faecal calprotectin (ug/g)
kd_cal : 1 : calprotectin equilibration (1/d)
cal_neut : 20 : neutrophil gain on calprotectin
cal_ulc : 15 : ulceration gain on calprotectin
CRP0 : 3 : baseline CRP (mg/L)
kd_crp : 0.6 : CRP turnover (1/d) — t1/2 19 h
k_ulc : 0.25 : ulceration formation (1/d)
k_ulch : 0.12 : ulceration healing (1/d)
ulc_thr : 0.8 : injury threshold for ulceration
ulc_max : 3 : ulceration ceiling
Alb0 : 4.2 : baseline albumin (g/dL)
kcat_alb : 0.04 : albumin catabolism (1/d)
k_ple : 0.03 : protein-losing enteropathy coefficient
alb_crp : 0.012 : CRP suppression of albumin synthesis
h_inf : 0.004 : infection hazard rate
// ---- steroid pharmacodynamic actions ----------------------------------
s_cyt : 0.85 : steroid suppression of cytokine production
s_prol : 0.7 : steroid suppression of effector proliferation
s_kill : 0.35 : steroid-induced effector apoptosis (1/d)
s_hom : 0.5 : steroid suppression of homing
s_kill_trm : 0.06 : steroid effect on resident memory (small — the point)
// ---- tumour -----------------------------------------------------------
k_hom_t : 0.3 : effector homing into tumour (NOT alpha4beta7-dependent)
k_prol_t : 0.2 : intratumoural effector proliferation (1/d)
k_out_t : 0.45 : intratumoural effector egress (1/d) — SOLVED for SS
TumTeff_max : 0.5 : intratumoural effector carrying capacity
kg_tum : 0.02 : tumour growth rate (1/d)
Tmax_tum : 1000 : tumour carrying capacity
kkill_tum : 45 : immune kill rate constant
Km_tum : 400 : kill saturation constant
chi_pred : 1 : prednisolone non-gut-selectivity
chi_ifx : 0.6 : infliximab non-gut-selectivity
chi_vdz : 0.1 : vedolizumab non-gut-selectivity — the gut ADDRESS
chi_toc : 0.5 : tocilizumab non-gut-selectivity
// ---- external inputs (set by the driver / Shiny app) ------------------
PRED_RATE : 0.0  : prednisolone input (mg/day, continuous)
JAK_RATE  : 0.0  : JAK inhibitor input (mg/day, continuous)
ABX_ON    : 0.0  : antibiotic exposure indicator (0/1)

[ GLOBAL ]
#define TEFF0 0.05
#define TRM0  0.02
#define TREG0 1.00
#define MAC0  0.20
#define NEUT0 0.05
#define TH170 0.10
#define DRIVE0 (TEFF0 + a_trm * TRM0)

// C1 max(0,x): EXACTLY zero for x <= 0 (so the drug-free mucosa is an exact
// steady state, not a slow drift), quadratically smoothed over the first s
// (so the Jacobian stays continuous when the tolerance band is crossed).
#define HINGE(x) ((x) <= 0.0 ? 0.0 : ((x) < hinge_s ? (x)*(x)/(2.0*hinge_s) \
                                                     : (x) - 0.5*hinge_s))
#define MM(x, k) ((x) / ((x) + (k)))

[ MAIN ]
// Two egress rates are not free parameters: they are whatever makes the
// drug-free colon and the untreated intratumoural pool stationary.  Solving
// them here is what makes the steady-state self-test pass exactly.
double ncl0 = kp_prime / kout_ncl;
double sp0  = 1.0 - TEFF0 / Teff_max;
double K_OUT_EFF = ((k_hom * ncl0 + k_prol * TEFF0) * sp0
                    - k_conv * TEFF0) / TEFF0;
double K_OUT_T   = ((k_hom_t * ncl0 + k_prol_t * TEFF0)
                    * (1.0 - TEFF0 / TumTeff_max)) / TEFF0;

[ ODE ]
// =========================================================================
// 0.  Plasma concentrations (mg/L == ug/mL)
// =========================================================================
double ALBs = Alb  > 1.2 ? Alb  : 1.2;
double CRPs = CRP  > 0.1 ? CRP  : 0.1;
double C_ipi = A_ipi1 / V1_ipi;
double C_pd1 = A_pd11 / V1_pd1;
double C_ifx = A_ifx1 / V1_ifx;
double C_vdz = A_vdz1 / V1_vdz;
double C_toc = A_toc1 / V1_toc;
double C_pred = A_pred1 / V_pred;
double C_jak = A_jak / V_jak;

// ---- THE FIFTH IDEA: disease severity raises the clearance of its antidote
double CLifx = CL_ifx * pow(ALBs / ALB_REF, alb_exp_ifx)
                      * (1.0 + crp_slope_ifx * CRPs);
double CLvdz = CL_vdz * pow(ALBs / ALB_REF, alb_exp_vdz);

// =========================================================================
// 1.  Tissue penetration and target engagement
// =========================================================================
// Inflammation raises vascular permeability.  Ulcer and calprotectin are used
// as the smooth monotone surrogate for mucosal inflammation so that the
// penetration term does not depend on quantities computed further down.
double infl = Ulcer / 1.5 + Calpro / 600.0;
if (infl > 1.0) infl = 1.0;
double bdc = BDC * (1.0 + BDC_infl * infl);

double Ct_ipi = C_ipi * bdc;
double Ct_pd1 = C_pd1 * bdc;
double Ct_ifx = C_ifx * bdc * 0.77;
double Ct_vdz = C_vdz * bdc * 0.77;
double Ct_toc = C_toc * bdc;

double nM_pd1 = Ct_pd1 * 1000.0 / MW_pd1;
double nM_ipi = Ct_ipi * 1000.0 / MW_ipi;
double O_PD1   = MM(nM_pd1, KD_PD1_nM);      // saturated across all doses
double O_CTLA4 = MM(nM_ipi, KD_CTLA4_nM);    // also saturated
double A_adcc  = MM(Ct_ipi, EC50_ADCC);      // NOT saturated — the dose dial

double E_ifx  = Emax_IFX * MM(Ct_ifx, EC50_IFX);
double O_vdz  = MM(Ct_vdz, KD_VDZ);
double E_vdz  = Emax_VDZ * O_vdz;
double E_toc  = Emax_TOC * MM(Ct_toc, EC50_TOC);
double E_pred = Emax_pred * MM(Ce_pred, EC50_pred);
double E_jak  = Emax_jak * MM(C_jak, EC50_jak);

// =========================================================================
// 2.  THE RATIO
// =========================================================================
double Tregp = Treg > 0.0 ? Treg : 0.0;
double Reg  = a_reg * Tregp + kappa;
double Reg0 = a_reg * TREG0 + kappa;
double G    = Reg0 / (Reg > 1e-6 ? Reg : 1e-6);   // hyperbolic, capped by kappa
double Ppd1 = 1.0 + Emax_pd1 * O_PD1;             // efferent-phase release
double drive = (Teff + a_trm * Trm) * Ppd1 * G;
double drive_rel = drive / DRIVE0;

// Saturating transduction: everything downstream of the T cells enters
// through this one bounded quantity.  Without it the model is a cascade of
// linear amplifiers and diverges the moment the tolerance band is crossed.
double drx = HINGE(drive_rel - 1.0);
double h_drive = MM(drx, K_cyt);

// =========================================================================
// 3.  Antigen availability (bounded feed-forward)
// =========================================================================
double LPSex = HINGE(LPS / LPS0 - 1.0);
double breach = LPSex + 2.0 * HINGE(1.0 - Muc);
double Agx = ag_scale * (1.0 + a_ff * MM(breach, K_ff));

// =========================================================================
// 4.  Signalling modifiers
// =========================================================================
double fs      = 1.0 - s_cyt * E_pred;            // steroid on cytokines
double TNF_eff = TNFa * (1.0 - E_ifx);            // infliximab neutralisation
double IL6_eff = IL6 * (1.0 - E_toc) * (1.0 - 0.6 * E_jak);
double IFN_sig = IFNg * (1.0 - 0.80 * E_jak);

// =========================================================================
// 5.  Injury: gated by the tolerance band, then saturated
// =========================================================================
double IFNp = IFN_sig > 1e-6 ? IFN_sig : 1e-6;
double cyt_kill = TNF_eff * pow(IFNp, 0.70);
double Injury = w_cyt * cyt_kill
              + w_gzmb * (1.0 + E_gzmb * h_drive)
              + w_neut * Neut / NEUT0;
double raw = HINGE(Injury - 1.0 - dead_inj);
double Ix  = Ix_max * MM(raw, K_ix);

// =========================================================================
// 6.  Epithelium, transport and THE CENSORING THRESHOLD
// =========================================================================
double Bup = Bu > 0.0 ? Bu : 0.0;
double Rep = (0.30 + 0.70 * Bup) * f_nsaid;
double f_tj = (1.0 - f_tj_w) + f_tj_w * (TJ > 0.0 ? TJ : 0.0);
double nhe = 1.0 / (1.0 + nhe_tnf * HINGE(TNF_eff - 1.0 - dead_tnf)
                        + nhe_ifn * HINGE(IFN_sig - 1.0 - dead_ifn));
double f_tr = (1.0 - f_tr_w) + f_tr_w * nhe;
double S_eff = (Ent > 0.0 ? Ent : 0.0) * f_tj * f_tr;
double S_star = L_pres / A_max;                   // = 0.389
double Sec = sec_lps * (LPSex / (1.0 + LPSex)) + sec_bile * HINGE(1.0 - Dv);
double Stool = L_pres + Sec - A_max * S_eff;
if (Stool < Stool_min) Stool = Stool_min;
double dfreq = (Stool - Stool_min) / mL_per_stool;
double freq = freq0 + dfreq;
double grade = dfreq < 1.0 ? 0.0 : (dfreq < 4.0 ? 1.0 : (dfreq < 7.0 ? 2.0 : 3.0));

// =========================================================================
// 7.  Systemic immunosuppression weighted by gut-selectivity
// =========================================================================
double Imm_sys = chi_pred * E_pred + chi_ifx * E_ifx
               + chi_vdz * E_vdz + chi_toc * E_toc;
if (Imm_sys > 0.95) Imm_sys = 0.95;

// =========================================================================
// 8.  PK ODEs
// =========================================================================
dxdt_A_ipi1 = -(CL_ipi/V1_ipi + Q_ipi/V1_ipi)*A_ipi1 + (Q_ipi/V2_ipi)*A_ipi2;
dxdt_A_ipi2 =  (Q_ipi/V1_ipi)*A_ipi1 - (Q_ipi/V2_ipi)*A_ipi2;
dxdt_A_pd11 = -(CL_pd1/V1_pd1 + Q_pd1/V1_pd1)*A_pd11 + (Q_pd1/V2_pd1)*A_pd12;
dxdt_A_pd12 =  (Q_pd1/V1_pd1)*A_pd11 - (Q_pd1/V2_pd1)*A_pd12;
dxdt_A_ifx1 = -(CLifx/V1_ifx + Q_ifx/V1_ifx)*A_ifx1 + (Q_ifx/V2_ifx)*A_ifx2;
dxdt_A_ifx2 =  (Q_ifx/V1_ifx)*A_ifx1 - (Q_ifx/V2_ifx)*A_ifx2;
dxdt_A_vdz1 = -(CLvdz/V1_vdz + Q_vdz/V1_vdz)*A_vdz1 + (Q_vdz/V2_vdz)*A_vdz2;
dxdt_A_vdz2 =  (Q_vdz/V1_vdz)*A_vdz1 - (Q_vdz/V2_vdz)*A_vdz2;
dxdt_A_toc1 = -(CL_toc/V1_toc + Q_toc/V1_toc)*A_toc1 + (Q_toc/V2_toc)*A_toc2;
dxdt_A_toc2 =  (Q_toc/V1_toc)*A_toc1 - (Q_toc/V2_toc)*A_toc2;

double abs_p = ka_pred * A_predg;
dxdt_A_predg = PRED_RATE - abs_p;
dxdt_A_pred1 = F_pred * abs_p - (CL_pred/V_pred) * A_pred1;
dxdt_Ce_pred = keo_pred * (C_pred - Ce_pred);
dxdt_A_jak   = JAK_RATE - ke_jak * A_jak;

// =========================================================================
// 9.  Priming — the AFFERENT phase, where anti-CTLA-4 acts
// =========================================================================
double Prime  = kp_prime * Agx * Tn * (1.0 + Emax_ctla * O_CTLA4)
                * sqrt(G > 1e-9 ? G : 1e-9);
double Prime0 = kp_prime;
dxdt_Tn     = lam_tn * (1.0 - Tn) - c_tn * (Prime - Prime0);
dxdt_Nclone = Prime - kout_ncl * Nclone;

// =========================================================================
// 10.  Colonic effector and regulatory compartments
// =========================================================================
double homing = MADCAM * (1.0 + g_amp * c_hom * HINGE(CXCL10 - 1.0 - dead_cx))
                * (1.0 - E_vdz) * (1.0 - s_hom * E_pred);
if (homing < 0.0) homing = 0.0;
double space = 1.0 - Teff / Teff_max;
if (space < 0.0) space = 0.0;
double influx = k_hom * Nclone * homing * space;

// NOTE.  Ppd1 enters expansion only with exponent pd1_prol_pow.  PD-1
// blockade is efferent: it raises the activity of clones that already exist
// (that is in `drive`), it does not license new ones.  Giving it full weight
// here makes anti-PD-1 monotherapy as colitogenic as ipilimumab, which is
// exactly what the clinical data say it is not.
double prolif = k_prol * Teff * Agx * G * pow(Ppd1, pd1_prol_pow) * space
                * (1.0 - s_prol * E_pred) * (1.0 - 0.4 * E_jak);
double conv = k_conv * Teff;
dxdt_Teff = influx + prolif - conv - K_OUT_EFF * Teff
            - s_kill * E_pred * Teff;

// Trm self-renewal needs a THRESHOLD level of epithelial IL-15, not merely
// any elevation.  Without this deadband the residency latch closes in every
// regimen and there is no such thing as a patient who was exposed and did
// not get colitis.
double il15ex = HINGE(IL15 - 1.0 - dead_il15);
double selfren = k_sr * Trm * MM(il15ex, K15) * (1.0 - Trm / Trm_max);
dxdt_Trm = conv + selfren - k_dtrm * Trm - s_kill_trm * E_pred * Trm;

// ADCC is SELF-AMPLIFYING: the inflammation it causes raises antibody
// penetration (bdc), recruits the FcgammaRIIIa+ effector cells (Mac) and
// depletes the butyrate that supplies pTreg.
double adcc = k_adcc * phi_FCGR * A_adcc
              * pow(Mac / MAC0 > 1e-9 ? Mac / MAC0 : 1e-9, mac_pow) * Treg;
dxdt_Treg = k_in_treg * (1.0 + b_but * (Bu - 1.0))
            - k_out_treg * Treg - adcc;

// ---- myeloid and Th17 (saturating; steroid suppresses only the excess) ---
double hcx = HINGE(CXCL10 - 1.0 - dead_cx);
double mac_t = MAC0 * (1.0 + fs * (g_amp * m_cx * MM(hcx, K_cx_mac)
                                   + m_lps * MM(LPSex, K_lps)));
dxdt_Mac = kd_mac * (mac_t - Mac);

double h23 = HINGE((IL23 - 1.0) * G);
double th17_t = TH170 * (1.0 + fs * E_th17 * MM(h23, K_th17));
dxdt_Th17 = kd_th17 * (th17_t - Th17);

double hn1 = HINGE(Th17 / TH170 - 1.0 - dead_th17);
double hn2 = HINGE(TNF_eff - 1.0 - dead_tnf);
double neut_t = NEUT0 * (1.0 + fs * (E_n_th17 * MM(hn1, K_neut)
                                     + E_n_tnf * MM(hn2, K_neut)));
dxdt_Neut = kd_neut * (neut_t - Neut);

// =========================================================================
// 11.  Cytokines — saturating targets, exact baseline of 1
// =========================================================================
double hmac = HINGE(Mac / MAC0 - 1.0);
double hmac_s = MM(hmac, K_cx_mac);
double hlps_s = MM(LPSex, K_lps);

dxdt_IFNg = kd_ifng * ((1.0 + fs * E_ifn * h_drive) - IFNg);
dxdt_TNFa = kd_tnf  * ((1.0 + fs * (E_tnf_m * hmac_s + E_tnf_d * h_drive
                                    + E_tnf_l * hlps_s)) - TNFa);
dxdt_IL6  = kd_il6  * ((1.0 + fs * (E_il6_m * hmac_s + E_il6_d * h_drive)) - IL6);
dxdt_IL23 = kd_il23 * ((1.0 + fs * E_il23 * hmac_s) - IL23);
dxdt_IL15 = kd_il15 * ((1.0 + b_il15 * MM(Ix, K_inj15)) - IL15);
double hcx_in = HINGE(IFN_sig - 1.0 - dead_ifn);
dxdt_CXCL10 = kd_cxcl * ((1.0 + b_cxcl * MM(hcx_in, K_ifn_cx)) - CXCL10);
double mad_t = 1.0 + g_amp * b_mad * MM(hn2, K_tnf_mad);
if (mad_t > mad_max) mad_t = mad_max;
dxdt_MADCAM = kd_mad * (mad_t - MADCAM);
dxdt_IL10 = kd_il10 * (Treg - IL10);

// =========================================================================
// 12.  Epithelium
// =========================================================================
dxdt_ISC = k_isc * (1.0 - ISC) - k_isck * HINGE(Ix - inj_deep) * ISC;
dxdt_Ent = k_prod_e * ISC * Rep - (k_turn + k_inj * Ix) * Ent;
double tj_dam = HINGE((TNF_eff + 0.5 * IFN_sig) / 1.5 - 1.0 - dead_tj);
dxdt_TJ = k_tj * (1.0 - TJ) - k_tjd * tj_dam * TJ;
dxdt_Muc = k_muc * (ISC * Rep - Muc) - k_mucd * Ix * Muc;

// =========================================================================
// 13.  Lumen
// =========================================================================
dxdt_Dv = k_dv * (1.0 - Dv) - k_abx * ABX_ON * Dv - k_dys * Ix * Dv;
dxdt_Bu = k_bu * (pow(Dv > 0.0 ? Dv : 0.0, 1.3) - Bu) - k_buc * Ix * Bu;
double Mucs = Muc > 0.25 ? Muc : 0.25;
double lps_t = LPS0 * (1.0 + b_lps * (1.0 - TJ) / Mucs);
dxdt_LPS = k_lps * (lps_t - LPS);

// =========================================================================
// 14.  Biomarkers and complications
// =========================================================================
double cal_t = Calpro0 * (1.0 + cal_neut * HINGE(Neut/NEUT0 - 1.0 - dead_neut)
                          + cal_ulc * Ulcer);
dxdt_Calpro = kd_cal * (cal_t - Calpro);
dxdt_CRP = kd_crp * CRP0 * IL6_eff - kd_crp * CRP;
dxdt_Ulcer = k_ulc * HINGE(Ix - ulc_thr) * (1.0 - Ulcer / ulc_max)
             - k_ulch * Ulcer * Rep;
double Entc = Ent < 1.0 ? Ent : 1.0;
double ksyn = kcat_alb * Alb0 * (1.0 + alb_crp * CRP0);
dxdt_Alb = ksyn / (1.0 + alb_crp * CRP)
           - (kcat_alb + k_ple * Ulcer * (1.0 - Entc)) * Alb;

// =========================================================================
// 15.  Tumour — THE FOURTH IDEA: selectivity is the currency
// =========================================================================
// Vedolizumab blocks alpha4beta7:MAdCAM-1, a GUT address.  The tumour does
// not use it, so vedolizumab enters this term with chi = 0.10 while
// prednisolone enters with 1.00.
double traffic_block = chi_pred * E_pred + chi_ifx * E_ifx + chi_vdz * E_vdz;
if (traffic_block > 0.95) traffic_block = 0.95;
double sp_t = 1.0 - TumTeff / TumTeff_max;
if (sp_t < 0.0) sp_t = 0.0;
dxdt_TumTeff = (k_hom_t * Nclone * (1.0 - traffic_block)
                + k_prol_t * TumTeff * Ppd1 * G * (1.0 - Imm_sys)) * sp_t
               - K_OUT_T * TumTeff;
double Tp = Tumor > 0.0 ? Tumor : 0.0;
dxdt_Tumor = kg_tum * Tp * (1.0 - Tp / Tmax_tum)
             - kkill_tum * TumTeff * Tp / (Tp + Km_tum);

dxdt_SterCum = PRED_RATE;
dxdt_InfHaz = h_inf * (E_pred + 0.6 * E_ifx + 0.15 * E_vdz + 0.3 * E_toc);

[ TABLE ]
capture O_PD1_o    = O_PD1;
capture O_CTLA4_o  = O_CTLA4;
capture A_adcc_o   = A_adcc;
capture G_o        = G;
capture Ppd1_o     = Ppd1;
capture drive_rel_o = drive_rel;
capture Injury_o   = Injury;
capture Ix_o       = Ix;
capture S_eff_o    = S_eff;
capture S_star_o   = S_star;
capture Stool_o    = Stool;
capture dfreq_o    = dfreq;
capture freq_o     = freq;
capture grade_o    = grade;
capture E_pred_o   = E_pred;
capture E_ifx_o    = E_ifx;
capture E_vdz_o    = E_vdz;
capture Imm_sys_o  = Imm_sys;
capture CL_ifx_o   = CLifx;
capture C_ipi_o    = C_ipi;
capture C_pd1_o    = C_pd1;
capture C_ifx_o    = C_ifx;
capture C_vdz_o    = C_vdz;

[ CAPTURE ] @annotated
O_PD1_o     : PD-1 receptor occupancy (saturated at all clinical doses)
A_adcc_o    : ADCC drive (the dose-dependent step)
G_o         : regulatory-release factor Reg0/Reg (hyperbolic)
S_eff_o     : effective colonic absorptive function
S_star_o    : censoring threshold L/A_max
dfreq_o     : increase in stools per day over baseline
grade_o     : CTCAE diarrhoea grade (0-3)
CL_ifx_o    : infliximab clearance (rises as albumin falls)

## ---------------------------------------------------------------------------
## DRIVER CODE — the treatment scenarios
## ---------------------------------------------------------------------------
## Everything below is ordinary R and is ignored by mread().  Source this file
## after mread() (or copy the block) to reproduce the scenarios reported in
## README.md and icic_reference_output.txt.
## ---------------------------------------------------------------------------
if (FALSE) {

library(mrgsolve); library(dplyr); library(tidyr); library(ggplot2)
mod <- mread("icic", "icic_mrgsolve_model.R")
BW <- 70

## --- the representative SUSCEPTIBLE patient -----------------------------
## The MEDIAN patient does not develop colitis on ipilimumab 3 mg/kg -- only
## ~7-12% do, and a model in which the median patient gets colitis would be
## wrong.  Deterministic scenarios therefore run on a named susceptible
## phenotype, every element of which is a published risk factor.
CASE_PARAMS <- list(ag_scale = 1.45,   # higher commensal antigen drive
                    kappa    = 0.11,   # thin non-Treg regulatory floor
                    phi_FCGR = 1.60)   # FCGR3A 158 V/V, efficient ADCC
CASE_INIT   <- list(Dv = 0.70, Bu = 0.70^1.3,   # antibiotic-thinned microbiome
                    Trm = 0.06)                # latent primed / resident pool
## severe (steroid-refractory) phenotype used for the rescue comparisons
SEVERE_PARAMS <- list(ag_scale = 1.75, kappa = 0.095, phi_FCGR = 1.60)
SEVERE_INIT   <- list(Dv = 0.60, Bu = 0.60^1.3, Trm = 0.12)

case_mod <- mod %>% param(CASE_PARAMS) %>% init(CASE_INIT)

## --- 0. steady-state check (median patient, no drug) --------------------
mod %>% mrgsim(end = 365, delta = 5) %>% as_tibble() %>%
  summarise(across(c(Ent, TJ, Treg, Calpro, S_eff_o), list(first = first,
                                                           last = last)))

## --- SCENARIO 1-3. the ipilimumab dose-response -------------------------
## Occupancy is saturated at all three doses; ADCC is not.  The dose-response
## lives entirely in the denominator of the ratio.
ipi_dose <- function(mgkg, n = 4)
  ev(amt = mgkg * BW, cmt = 1, ii = 21, addl = n - 1)

s_ipi <- bind_rows(lapply(c(1, 3, 10), function(d)
  case_mod %>% ev(ipi_dose(d)) %>% mrgsim(end = 210, delta = 0.5) %>%
    as_tibble() %>% mutate(arm = paste0("ipilimumab ", d, " mg/kg"))))

## --- SCENARIO 4-5. anti-PD-1 is FLAT ------------------------------------
## 3 vs 10 mg/kg changes occupancy from 0.9985 to 0.9995.  There is no
## dose-response to find because the drug is on a plateau everywhere.
s_pd1 <- bind_rows(lapply(c(3, 10), function(d)
  case_mod %>% ev(amt = d * BW, cmt = 3, ii = 14, addl = 12) %>%
    mrgsim(end = 210, delta = 0.5) %>% as_tibble() %>%
    mutate(arm = paste0("anti-PD-1 ", d, " mg/kg"))))

## --- SCENARIO 6. combination (ipi 3 + nivo 1 -> nivo 480 q4w) -----------
e_combo <- ev(amt = 3 * BW, cmt = 1, ii = 21, addl = 3) +
           ev(amt = 1 * BW, cmt = 3, ii = 21, addl = 3) +
           ev(time = 84, amt = 480, cmt = 3, ii = 28, addl = 3)
s_combo <- case_mod %>% ev(e_combo) %>% mrgsim(end = 210, delta = 0.5)

## --- SCENARIO 7. steroid, started at grade 2 ----------------------------
## PRED_RATE is a continuous mg/day input: prednisolone has a 3 h half-life
## but a genomic effect with a ~12 h effect half-life, so a daily dose and a
## continuous rate are pharmacodynamically equivalent here.
## 1 mg/kg for 7 days, then a linear taper over `taper_wk` weeks.
pred_sched <- function(start_d, mgkg = 1, plateau = 7, taper_wk = 4,
                       end = 210) {
  top <- mgkg * BW
  tt <- seq(0, end, by = 1)
  rate <- sapply(tt - start_d, function(tau) {
    if (tau < 0) return(0)
    if (tau < plateau) return(top)
    wk <- (tau - plateau) / 7
    if (wk >= taper_wk) return(0)
    top * (1 - wk / taper_wk)
  })
  data.frame(time = tt, PRED_RATE = rate)
}
## grade 2 is reached on day 53 in this patient (see icic_reference_output.txt)
## Pass the taper as a time-varying covariate data set.
d_ster <- pred_sched(51) %>% mutate(ID = 1, evid = 0, cmt = 0, amt = 0) %>%
  bind_rows(as.data.frame(ipi_dose(3)) %>% mutate(ID = 1, PRED_RATE = NA)) %>%
  arrange(time) %>% fill(PRED_RATE, .direction = "downup")
s_ster <- case_mod %>% data_set(d_ster) %>% mrgsim(end = 210, delta = 0.5)

## --- SCENARIO 8. steroid-refractory -> infliximab 5 mg/kg at day 5 ------
## TNF is the TERMINAL executor of epithelial apoptosis, so this is the
## fastest-acting option (2-3 days).  Note CL_ifx_o in the output: the
## hypoalbuminaemia of a severe colitis has already raised the clearance of
## the drug you are about to give.
d_ifx <- d_ster %>% bind_rows(
  data.frame(ID = 1, time = 56, amt = 5 * BW, cmt = 5, evid = 1),
  data.frame(ID = 1, time = 70, amt = 5 * BW, cmt = 5, evid = 1),
  data.frame(ID = 1, time = 98, amt = 5 * BW, cmt = 5, evid = 1)) %>%
  arrange(time) %>% fill(PRED_RATE, .direction = "downup")
s_ifx <- case_mod %>% data_set(d_ifx) %>% mrgsim(end = 210, delta = 0.5)

## --- SCENARIO 9. steroid-refractory -> vedolizumab 300 mg wk 0/2/6 ------
## alpha4beta7 blockade acts on TRAFFICKING and cannot touch a cell already
## resident in the mucosa: slowest onset (14-21 d), most durable, and by far
## the cheapest in tumour control (chi = 0.10 vs 1.00 for steroid).
d_vdz <- d_ster %>% bind_rows(
  data.frame(ID = 1, time = c(56, 70, 98), amt = 300, cmt = 7, evid = 1)) %>%
  arrange(time) %>% fill(PRED_RATE, .direction = "downup")
s_vdz <- case_mod %>% data_set(d_vdz) %>% mrgsim(end = 210, delta = 0.5)

## --- SCENARIO 10. calprotectin-guided start (the uncensored instrument) --
## Calprotectin crosses 200 ug/g on day 45; grade 2 arrives on day 51.
## Starting on the biomarker treats a smaller lesion, leaves less deep
## ulceration and a smaller residual Trm pool -- it buys durability, not
## just speed.  It costs steroid given to patients who would never have
## become symptomatic, and the model quantifies both sides.
d_early <- pred_sched(45) %>% mutate(ID = 1, evid = 0, cmt = 0, amt = 0) %>%
  bind_rows(as.data.frame(ipi_dose(3)) %>% mutate(ID = 1, PRED_RATE = NA)) %>%
  arrange(time) %>% fill(PRED_RATE, .direction = "downup")
s_early <- case_mod %>% data_set(d_early) %>% mrgsim(end = 210, delta = 0.5)

## --- SCENARIO 11. FCGR3A genotype at a fixed ipilimumab dose ------------
s_fcgr <- bind_rows(lapply(c(V_V = 1.60, V_F = 1.00, F_F = 0.55), function(phi)
  case_mod %>% param(phi_FCGR = phi) %>% ev(ipi_dose(3)) %>%
    mrgsim(end = 210, delta = 1) %>% as_tibble()), .id = "genotype")

## --- SCENARIO 12. antibiotic pre-exposure (diversity/butyrate axis) -----
s_abx <- bind_rows(lapply(c(none = 1.00, recent = 0.55, severe = 0.35),
  function(dv) case_mod %>% init(Dv = dv, Bu = dv^1.3) %>% ev(ipi_dose(3)) %>%
    mrgsim(end = 210, delta = 1) %>% as_tibble()), .id = "microbiome")

## --- SCENARIO 13. the oncological cost of each rescue -------------------
## Same colitis control, different currency.  Compare Tumor at day 180.
bind_rows(
  s_ster  %>% as_tibble() %>% mutate(arm = "steroid only"),
  s_ifx   %>% as_tibble() %>% mutate(arm = "steroid + infliximab"),
  s_vdz   %>% as_tibble() %>% mutate(arm = "steroid + vedolizumab")) %>%
  group_by(arm) %>% filter(time == max(time)) %>%
  select(arm, Tumor, SterCum, InfHaz, Trm)

## --- plots --------------------------------------------------------------
## The censoring reserve, drawn: S_eff falls for weeks with a flat stool chart
s_ipi %>% ggplot(aes(time, S_eff_o, colour = arm)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = unique(s_ipi$S_star_o), linetype = 2) +
  annotate("text", x = 5, y = 0.42, hjust = 0, size = 3,
           label = "S* = L/A_max = 0.389 — below this, and only below this,\nthe stool chart moves") +
  labs(x = "day", y = "effective absorptive function S_eff",
       title = "61% of colonic absorptive function is spent before symptom 1")

## The two exposure-response shapes side by side
s_ipi %>% ggplot(aes(time)) +
  geom_line(aes(y = O_CTLA4_o, colour = "CTLA-4 occupancy (saturated)")) +
  geom_line(aes(y = A_adcc_o, colour = "ADCC drive (unsaturated)")) +
  facet_wrap(~arm) +
  labs(x = "day", y = "fraction of maximum", colour = NULL,
       title = "Why ipilimumab has a dose-response: it is not the occupancy")
}
