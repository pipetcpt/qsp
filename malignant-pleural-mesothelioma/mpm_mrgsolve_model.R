## =============================================================================
##  mpm_mrgsolve_model.R
##  Malignant pleural mesothelioma (MPM)
##  Quantitative Systems Pharmacology model for mrgsolve
##
##  53 ODE compartments.  TIME UNIT = DAYS, because this disease forces a
##  40-minute free-platinum half-life, a 3-week cycle, a 21-day necrotic-debris
##  clearance, a 460-day collagen turnover and a 3-year survival horizon to live
##  in one system.
##
##  ---------------------------------------------------------------------------
##  THE ONE STRUCTURAL COMMITMENT
##  ---------------------------------------------------------------------------
##
##  Mesothelioma is a RIND, not a mass.  Tumour lives as a sheet of thickness h
##  spread over a pleural surface of area A, so
##
##          V = A h        and        dV/dh = A , a CONSTANT.
##
##  Tumour burden is therefore NOT a state variable in this model.  The state
##  variables are MASSES on three surfaces (parietal, visceral, fissural), and
##  thickness is computed from them:
##
##          h_i = (T_i + M_i + N_i) / (rho A_i)
##                 ^^^^  ^^^^  ^^^^
##                 viable  matrix  necrotic
##
##  Every clinical quantity is then a different FUNCTIONAL of that same sheet:
##
##      mRECIST      = 4 h_par + 2 h_vis          -> LINEAR in h
##      drug at the
##      cell         = (lam/h)(1 - exp(-h/lam))   -> SUBLINEAR, saturating
##      lung trapping= stiffness x h_vis           -> VISCERAL LEAF ONLY
##      effusion     = f(parietal rind sitting on the lymphatic stomata)
##
##  Because dV/dh is a constant, a 30% fall in mRECIST is a 30% fall in tumour
##  VOLUME.  In the spherical tumour that RECIST was designed for, a 30% fall in
##  diameter is a 66% fall in volume.  The phrase "partial response" therefore
##  denotes a 2.2-fold different degree of cytoreduction in a mesothelioma trial
##  than in the trial next door.  This is arithmetic, and the model reports it.
##
##  ---------------------------------------------------------------------------
##  ONE COLLAGEN STATE, THREE CONSEQUENCES
##  ---------------------------------------------------------------------------
##
##  phi = M / (T + M + N), the collagenous fraction of the rind, enters the
##  model in exactly three places:
##
##    (1) in h, because matrix mass does not disappear when tumour cells die.
##        Running the calibrated model through six cycles of pemetrexed +
##        cisplatin kills 60.3% of the viable cell mass and moves the mRECIST
##        sum by only -30.1%.  The missing 30 points are collagen and
##        not-yet-cleared necrotic debris, and the gap is remarkably stable
##        across regimens: 27.7 points for pemetrexed/cisplatin, 28.2 with
##        bevacizumab, 28.4 for nivolumab + ipilimumab.  Every thickness
##        endpoint in this disease is biased downward by about the same amount.
##
##    (2) in lambda(phi) = lambda_0 (1-phi)^1.5, because collagen collapses
##        microvessels and raises interstitial pressure.  So the tumour becomes
##        harder to reach EXACTLY AS IT RESPONDS: phi rises from 0.265 to 0.499
##        across those six cycles, and the depth-averaged small-molecule
##        exposure falls from 0.257 to 0.215 even though the rind got thinner.
##
##    (3) in the added elastance of the visceral leaf, E_add = k phi h_vis.
##        This is the reason the model predicts something initially surprising
##        and clinically familiar: across a 30% radiographic RESPONSE, E_add
##        RISES from 8.15 to 10.01 cmH2O/L and FVC does not improve (2.27 ->
##        2.11 L).  A responding mesothelioma does not breathe better.  That is
##        why large radiographic responses in this disease buy small survival
##        gains, and the model derives it rather than assuming it.
##
##  ---------------------------------------------------------------------------
##  TWO DELIVERY ROUTES WITH OPPOSITE GRADIENTS
##  ---------------------------------------------------------------------------
##
##  Systemic drug enters at the VASCULARISED BASE of the rind and decays inward.
##  Intrapleural drug enters at the FREE SURFACE and decays toward the base.
##  Both use the same depth-average, with different lengths:
##
##      lambda_smallmolecule ~ 2.5 mm      lambda_IgG ~ 0.6 mm
##      lambda_Tcell         ~ 0.5 mm      lambda_intrapleural ~ 3.5 mm
##
##  At phi = 0.25 a 6 mm rind sees a depth-averaged interstitial antibody
##  concentration of 6.5% of the plasma-equilibrium value -- 18% of what a 1 mm
##  rind sees at the same plasma exposure.  No dose escalation fixes that: f is
##  bounded above by lambda/h, so doubling the dose doubles a number that is
##  already 0.065.  It is also why every
##  positive intracavitary result in this disease is reported AFTER macroscopic
##  complete resection and never instead of it.
##
##  ---------------------------------------------------------------------------
##  HISTOLOGY AS ONE CONTINUOUS AXIS
##  ---------------------------------------------------------------------------
##
##  EMT in [0,1]: 0 epithelioid, 0.5 biphasic, 1 sarcomatoid.  One parameter
##  raises proliferation (x1.9 at x=1), raises chemoresistance (kill x 0.45),
##  raises collagen deposition (x4.2) and raises PD-L1 (0.20 -> 0.60) and
##  T-cell infiltration (x2).  The first three hurt chemotherapy and the fourth
##  helps checkpoint blockade, so the SIGN of OS(IO) - OS(chemo) along x is an
##  OUTPUT.  Only ONE point on that axis was fitted (the CheckMate 743 trial
##  population); the epithelioid / non-epithelioid split was left as a
##  prediction -- AND THE PREDICTION FAILED.  The model has the IO advantage
##  SHRINKING along the axis (median ratio 0.70 epithelioid -> 0.87
##  non-epithelioid) where CheckMate 743 has it GROWING (HR 0.86 -> 0.46).
##  This is reported rather than tuned away; mpm_emt_sensitivity.py shows that
##  neither the chemoresistance slope nor PD-L1 can repair it, and locates the
##  error in the fact that this model applies the histology penalty to BOTH
##  arms while the trial says it falls on chemotherapy alone.
##
##  ---------------------------------------------------------------------------
##  PROVENANCE
##  ---------------------------------------------------------------------------
##  Every equation below was integrated first in mpm_reference_model.py (pure
##  Python, RK4, no dependencies) because this build environment has no R
##  runtime.  mpm_calibration.py fits five constants by one-dimensional
##  bisection, one endpoint each, and prints mpm_calibration_output.txt.  The
##  numbers quoted in the comments above come from that run.
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

code <- '
$PROB
# Malignant pleural mesothelioma QSP model
# Tumour as a rind: thickness is derived, not integrated.

$PARAM @annotated
// ---------------- patient ------------------------------------------------
BSA      : 1.80   : Body surface area (m2)
WT       : 72     : Body weight (kg)
FVC_PRED : 3.30   : Predicted FVC (L)
EMT      : 0.25   : Histology axis 0=epithelioid 0.5=biphasic 1=sarcomatoid

// ---------------- rind geometry ------------------------------------------
A_PAR_ANAT : 1000 : Parietal pleural surface one hemithorax (cm2)
A_VIS_ANAT : 900  : Visceral (lung) surface (cm2)
A_FIS_ANAT : 400  : Fissural + mediastinal + diaphragmatic surface (cm2)
COV_PAR  : 0.45   : Fraction of parietal surface carrying tumour
COV_VIS  : 0.35   : Fraction of visceral surface carrying tumour
COV_FIS  : 0.50   : Fraction of fissural surface carrying tumour
RHO      : 1.05   : Tissue density (g/cm3)
HPAR0    : 0.60   : Parietal rind thickness at diagnosis (cm)
HVIS0    : 0.40   : Visceral rind thickness at diagnosis (cm)
HFIS0    : 0.50   : Fissural rind thickness at diagnosis (cm)

// ---------------- growth --------------------------------------------------
KG0      : 0.0032 : Intrinsic proliferation rate of a thin rind (1/d)
EMT_KG   : 0.90   : Proliferation multiplier per unit EMT
LAM_G    : 0.42   : Depth of the vascularised proliferative zone (cm)
KSEED    : 0.0006 : Cross-seeding between pleural surfaces (1/d)

// ---------------- matrix / collagen ---------------------------------------
KCOL0    : 0.00125: Matrix deposited per gram of tumour (1/d)
EMT_KCOL : 3.20   : Collagen multiplier per unit EMT
KMDEG    : 0.0015 : Matrix turnover (1/d, t1/2 ~ 460 d)
KNCLR    : 0.0476 : Clearance of killed mass (1/d, t1/2 ~ 15 d)

// ---------------- penetration ---------------------------------------------
LAM_SM0  : 0.25   : Small-molecule penetration length at phi=0 (cm)
LAM_AB0  : 0.060  : IgG penetration length at phi=0 (cm)
LAM_T0   : 0.050  : T-cell infiltration length at phi=0 (cm)
LAM_IP0  : 0.35   : Intrapleural penetration length at phi=0 (cm)
PHI_EXP  : 1.5    : Exponent in lambda = lambda0 (1-phi)^PHI_EXP

// ---------------- cisplatin PK --------------------------------------------
CIS_V1   : 25     : Central volume free platinum (L)
CIS_V2   : 40     : Peripheral volume (L)
CIS_CL   : 600    : Clearance free platinum (L/d)
CIS_Q    : 200    : Intercompartmental clearance (L/d)
CIS_KIN  : 1.00   : Tumour uptake / adduct formation (1/d)
CIS_KOUT : 0.35   : Pt-DNA adduct repair by NER (1/d)

// ---------------- pemetrexed PK -------------------------------------------
PEM_V1   : 16     : Central volume (L)
PEM_V2   : 8      : Peripheral volume (L)
PEM_CL0  : 130    : Clearance at eGFR 90 (L/d)
PEM_Q    : 30     : Intercompartmental clearance (L/d)
PEM_KIN  : 1.00   : RFC/SLC19A1-mediated uptake (1/d)
PEM_KOUT : 0.277  : Polyglutamate loss (1/d, t1/2 2.5 d)

// ---------------- gemcitabine PK ------------------------------------------
GEM_V    : 50     : Volume (L)
GEM_CL   : 2500   : Clearance (L/d)
GEM_KIN  : 1.00   : Tumour uptake (1/d)
GEM_KOUT : 1.00   : dFdCTP loss (1/d)

// ---------------- monoclonal antibody PK ----------------------------------
MAB_V1   : 3.6    : Central volume (L)
MAB_V2   : 2.8    : Peripheral volume (L)
MAB_CL   : 0.21   : Clearance (L/d)
MAB_Q    : 0.45   : Intercompartmental clearance (L/d)
MAB_KIN  : 0.55   : Into rind interstitium (1/d)
MAB_KOUT : 0.55   : Out of rind interstitium (1/d)

// ---------------- chemotherapy PD -----------------------------------------
EMAX_CIS : 0.148438 : Maximal cisplatin kill of a thin rind (1/d) FITTED
EC50_CIS : 0.10   : Adduct-equivalent EC50 (mg/L)
EMAX_PEM : 0.028359 : Maximal pemetrexed kill of a thin rind (1/d) FITTED
EC50_PEM : 0.35   : Polyglutamate-equivalent EC50 (mg/L)
SYN_PEMCIS : 0.30 : Fractional synergy of the doublet
EMT_CHEMO: 0.55   : Chemoresistance per unit EMT
EMAX_GEM : 0.030  : Maximal gemcitabine kill (1/d)
EC50_GEM : 0.35   : EC50 (mg/L)

// ---------------- immune ---------------------------------------------------
TEFF0    : 0.050  : Baseline intratumoural CD8 (arbitrary units)
KPRIME   : 0.055  : Priming rate from tumour antigen (1/d)
TMB_FACT : 0.42   : Antigenicity scaling for a low-TMB tumour
KTDEATH  : 0.070  : Effector loss (1/d)
KTREG    : 0.020  : Treg turnover (1/d)
KIFN     : 0.9    : IFN-gamma production
KIFN_OUT : 1.2    : IFN-gamma loss (1/d)
PDL1_BASE: 0.20   : Constitutive PD-L1 at EMT = 0
EMT_PDL1 : 2.0    : PD-L1 multiplier per unit EMT (0.20 -> 0.60)
EMT_TINF : 1.00   : Infiltration multiplier per unit EMT
PDL1_IND : 0.55   : IFN-gamma inducible PD-L1 component
EMAX_IO  : 0.624424 : Maximal immune kill of a thin rind (1/d) FITTED
KT50     : 0.080  : T_eff for half-maximal immune kill (its operating range)
KD_PD1   : 0.30   : Nivolumab conc for 50% PD-1 occupancy (mg/L)
KD_CTLA4 : 1.20   : Ipilimumab conc for 50% CTLA-4 occupancy (mg/L)
CTLA4_REL: 0.25   : Kill released by CTLA-4 blockade alone
IPI_PRIME: 2.30   : Priming gain at full CTLA-4 blockade
TGFB_SUPP: 0.85   : TGF-beta suppression of effector function

// ---------------- signalling ----------------------------------------------
KVEGF_OUT: 0.7    : VEGF turnover (1/d)
KTGFB_OUT: 0.8    : TGF-beta turnover (1/d)
KIL6_OUT : 0.9    : IL-6 turnover (1/d)
BEV_KD   : 0.9    : Bevacizumab conc for 50% VEGF neutralisation (mg/L)
BEV_LAM  : 0.295625 : Maximal fractional rise in lambda (normalisation) FITTED

// ---------------- pleural space -------------------------------------------
QFORM0   : 250    : Normal pleural fluid formation (mL/d)
QFORM_VEGF: 5.5   : Formation multiplier at maximal VEGF
QDRAIN_MAX: 5000  : Stomatal lymphatic reserve (mL/d)
H50_STOMA: 0.22   : Parietal rind for 50% stomatal obstruction (cm)
PLV_MAX  : 3500   : Hemithorax fluid capacity (mL)
PLV0     : 628    : Pleural effusion at diagnosis (mL)
KSYMPH   : 0.10   : Pleurodesis symphysis formation (1/d)
KSYMPH_LOSS: 0.010: Symphysis loss (1/d)

// ---------------- mechanics -----------------------------------------------
KELAST   : 76     : Added elastance per (phi x cm) of visceral rind
E50      : 14.5   : Trapped-lung elastance threshold (cmH2O/L)
KREC     : 0.045  : Lung re-expansion rate (1/d)
PLV_FVC  : 0.55   : Fraction of FVC lost at maximal effusion
KDYSP    : 0.20   : Dyspnoea equilibration (1/d)

// ---------------- toxicity ------------------------------------------------
CIRC0    : 4.0    : Baseline ANC (10^9/L)
MTT      : 4.6    : Marrow mean transit time (d)
GAM_FB   : 0.16   : Rebound feedback exponent
SLOPE_CIS_ANC : 0.90 : Cisplatin marrow slope on the retained adduct (L/mg)
SLOPE_PEM_ANC : 0.115: Pemetrexed marrow slope on the polyglutamate (L/mg)
FOLATE   : 1      : 1 = folic acid + B12 supplemented
FOLATE_PROT : 0.60: Reduction in pemetrexed marrow effect with supplementation
KPT_UP   : 35.0   : Renal cortex platinum uptake
KPT_EL   : 0.020  : Renal platinum elimination (1/d)
GFR0     : 90     : Baseline eGFR (mL/min/1.73m2)
KGFR_DAM : 0.400  : Platinum nephrotoxicity slope
KGFR_REP : 0.0045 : Tubular repair (1/d)
KNEURO   : 0.750  : Cumulative neuropathy slope
KIRAE_PD1: 0.020  : irAE hazard from PD-1 occupancy
KIRAE_CTLA4: 0.075: irAE hazard from CTLA-4 occupancy
KIRAE_OUT: 0.035  : irAE resolution (1/d)

// ---------------- cachexia ------------------------------------------------
KCACHEX  : 0.055  : Weight loss at maximal IL-6 (%/d)
KCACHEX_REC : 0.010 : Weight regain (1/d)

// ---------------- biomarker -----------------------------------------------
KSMRP    : 0.00090: SMRP shedding per gram of viable tumour (nM/d/g)
KSMRP_CL : 0.60   : SMRP renal clearance at eGFR 90 (1/d)

// ---------------- survival ------------------------------------------------
HZ0      : 0.00045625 : Baseline hazard at the reference state (1/d) FITTED
HZ_VOL   : 0.30   : Coefficient on ln(V / V_ref)
HZ_VREF  : 500    : Reference tumour volume (cm3)
HZ_ECOG  : 0.60   : Coefficient per ECOG point
HZ_EMT   : 0.55   : Coefficient per unit of the histology axis
HZ_FVC   : 1.20   : Coefficient per unit fractional FVC loss

// ---------------- procedures / switches -----------------------------------
ON_ADI   : 0      : Arginine deprivation (ADI-PEG20) on
ASS1_NEG : 0      : 1 = ASS1-methylated tumour (arginine auxotroph)
ADI_EFF  : 0.026  : Arginine-deprivation kill rate (1/d)
ON_IP    : 0      : Intrapleural chemotherapy on
SURG_DAY : -1     : Day of cytoreduction (-1 = none)
SURG_RESID : 0.10 : Residual thickness after macroscopic complete resection (cm)
SURG_FIS_SPARE : 0.65 : Fraction of fissural disease NOT resectable
SURG_FVC_HIT : 0  : 0 = P/D, 0.35 = EPP (a whole lung removed)
TALC_DAY : -1     : Day of talc pleurodesis (-1 = none)
IPC_DAY  : -1     : Day of indwelling pleural catheter (-1 = none)
IPC_RATE : 700    : Ambulatory drainage (mL/d)
TTF      : 0      : Tumour-treating fields on
TTF_EFF  : 0.0045 : TTFields kill rate (1/d)

$CMT @annotated
// ---- PK ----
CIS1  : Free platinum central (mg)
CIS2  : Free platinum peripheral (mg)
CIST  : Pt-DNA adduct equivalent in tumour (mg/L)
PEM1  : Pemetrexed central (mg)
PEM2  : Pemetrexed peripheral (mg)
PEMT  : Pemetrexed polyglutamate in tumour (mg/L)
GEM1  : Gemcitabine central (mg)
GEMT  : dFdCTP in tumour (mg/L)
BEV1  : Bevacizumab central (mg)
BEV2  : Bevacizumab peripheral (mg)
NIV1  : Nivolumab (or pembrolizumab) central (mg)
NIV2  : Nivolumab peripheral (mg)
NIVT  : Nivolumab in rind interstitium (mg/L)
IPI1  : Ipilimumab central (mg)
IPI2  : Ipilimumab peripheral (mg)
IPPL  : Intrapleural drug in the effusion (mg)
IPT   : Intrapleural drug in the rind (mg/L)
// ---- rind ----
TPAR  : Viable tumour mass parietal leaf (g)
TVIS  : Viable tumour mass visceral leaf (g)
TFIS  : Viable tumour mass fissural sanctuary (g)
MPAR  : Matrix mass parietal (g)
MVIS  : Matrix mass visceral (g)
MFIS  : Matrix mass fissural (g)
NPAR  : Necrotic mass parietal (g)
NVIS  : Necrotic mass visceral (g)
NFIS  : Necrotic mass fissural (g)
// ---- immune ----
TEFF  : Intratumoural effector CD8 (a.u.)
TREG  : Regulatory T cells (a.u.)
PRIME : Lymph-node primed T-cell pool (a.u.)
IFNG  : Intratumoural IFN-gamma (a.u.)
PDL1  : Tumour PD-L1 fraction
TAM   : M2 tumour-associated macrophages (a.u.)
// ---- signalling ----
VEGF  : VEGF-A (a.u.)
TGFB  : TGF-beta (a.u.)
IL6   : IL-6 (a.u.)
// ---- pleural space ----
PLV   : Pleural effusion volume (mL)
SYMPH : Pleural symphysis fraction
// ---- mechanics ----
VEXP  : Expandable lung volume fraction
DYSP  : Dyspnoea score
// ---- toxicity ----
PROL  : Proliferating marrow pool (10^9/L)
TR1   : Marrow transit 1
TR2   : Marrow transit 2
TR3   : Marrow transit 3
ANC   : Absolute neutrophil count (10^9/L)
PTK   : Renal cortex platinum (a.u.)
GFR   : eGFR (mL/min/1.73m2)
IRAE  : Immune-related adverse event driver
NEURO : Cumulative neuropathy index
// ---- host ----
CACHEX: Weight loss (%)
SMRP  : Serum soluble mesothelin-related peptide (nM)
ARG   : Plasma arginine (fraction of normal)
CH    : Cumulative hazard

$GLOBAL
#define A_PAR (A_PAR_ANAT * COV_PAR)
#define A_VIS (A_VIS_ANAT * COV_VIS)
#define A_FIS (A_FIS_ANAT * COV_FIS)

// Depth-averaged relative concentration through a slab of thickness h fed from
// ONE face with decay length lam:  <C>/C0 = (lam/h)(1 - exp(-h/lam)).
// This one function makes drug delivery, oxygen supply and T-cell infiltration
// all SUBLINEAR in rind thickness.
double fpen(double h, double lam) {
  if (h <= 1e-9) return 1.0;
  double r = h / lam;
  if (r > 60.0) return lam / h;
  return (lam / h) * (1.0 - exp(-r));
}
double lam_of(double lam0, double phi, double ex, double norm) {
  double f = pow(1.0 - (phi > 0.95 ? 0.95 : phi), ex) * (1.0 + norm);
  return (f * lam0 < 1e-4) ? 1e-4 : f * lam0;
}
double pos(double v) { return v > 0.0 ? v : 0.0; }

$MAIN
// ---- initial rind masses from the prescribed THICKNESSES -----------------
// phi at diagnosis is the steady state implied by deposition vs turnover at
// the intrinsic growth rate -- it is not a free parameter.
double kg_i   = KG0   * (1.0 + EMT_KG   * EMT);
double kcol_i = KCOL0 * (1.0 + EMT_KCOL * EMT);
double phi0   = kcol_i / (kcol_i + KMDEG + kg_i);

TPAR_0 = HPAR0 * RHO * A_PAR * (1.0 - phi0);
MPAR_0 = HPAR0 * RHO * A_PAR * phi0;
TVIS_0 = HVIS0 * RHO * A_VIS * (1.0 - phi0);
MVIS_0 = HVIS0 * RHO * A_VIS * phi0;
TFIS_0 = HFIS0 * RHO * A_FIS * (1.0 - phi0);
MFIS_0 = HFIS0 * RHO * A_FIS * phi0;

TEFF_0  = TEFF0 * (1.0 + EMT_TINF * EMT);
TREG_0  = 1.0;
PRIME_0 = 1.0;
IFNG_0  = 0.35;
PDL1_0  = PDL1_BASE * (1.0 + EMT_PDL1 * EMT);
TAM_0   = 1.0;
VEGF_0  = 1.0;
TGFB_0  = 1.0;
IL6_0   = 1.0;
PLV_0   = PLV0;

// the fast mechanical states start at the value implied by the initial rind
double Ead0  = KELAST * phi0 * HVIS0;
double vex0  = 1.0 / (1.0 + pow(Ead0 / E50, 2.0));
VEXP_0 = vex0;
double fvc0  = FVC_PRED * vex0 * (1.0 - PLV_FVC * PLV0 / PLV_MAX);
DYSP_0 = 2.4 * (1.0 - fvc0 / FVC_PRED) + 1.4 * PLV0 / PLV_MAX;

PROL_0 = CIRC0; TR1_0 = CIRC0; TR2_0 = CIRC0; TR3_0 = CIRC0; ANC_0 = CIRC0;
GFR_0  = GFR0;
ARG_0  = 1.0;

// Macroscopic complete resection is NOT applied here.  It is a state RESET,
// delivered by evid = 8 replacement records built in surgery_records() below,
// because surgery sets THICKNESS -- it does not reset phi and it cannot reach
// the fissural sanctuary.

$ODE
// ======================= rind geometry ====================================
double tp = pos(TPAR), tv = pos(TVIS), tf = pos(TFIS);
double mp = pos(MPAR), mv = pos(MVIS), mf = pos(MFIS);
double np_ = pos(NPAR), nv = pos(NVIS), nf = pos(NFIS);

double totp = tp + mp + np_;
double totv = tv + mv + nv;
double totf = tf + mf + nf;

double h_par = totp / (RHO * A_PAR);
double h_vis = totv / (RHO * A_VIS);
double h_fis = totf / (RHO * A_FIS);

double phi_par = (totp > 1e-9) ? mp / totp : 0.0;
double phi_vis = (totv > 1e-9) ? mv / totv : 0.0;
double phi_fis = (totf > 1e-9) ? mf / totf : 0.0;

double V_tot = A_PAR * h_par + A_VIS * h_vis + A_FIS * h_fis;
double T_tot = tp + tv + tf;
double N_tot = np_ + nv + nf;

// ======================= penetration ======================================
double c_bev = pos(BEV1) / MAB_V1;
double norm  = BEV_LAM * c_bev / (BEV_KD + c_bev);   // vascular normalisation

double f_sm_p = fpen(h_par, lam_of(LAM_SM0, phi_par, PHI_EXP, norm));
double f_sm_v = fpen(h_vis, lam_of(LAM_SM0, phi_vis, PHI_EXP, norm));
double f_sm_f = fpen(h_fis, lam_of(LAM_SM0, phi_fis, PHI_EXP, norm));
double f_ab_p = fpen(h_par, lam_of(LAM_AB0, phi_par, PHI_EXP, norm));
double f_T_p  = fpen(h_par, lam_of(LAM_T0,  phi_par, PHI_EXP, norm));
double f_T_v  = fpen(h_vis, lam_of(LAM_T0,  phi_vis, PHI_EXP, norm));
double f_T_f  = fpen(h_fis, lam_of(LAM_T0,  phi_fis, PHI_EXP, norm));
double f_ip_p = fpen(h_par, lam_of(LAM_IP0, phi_par, PHI_EXP, 0.0));
double f_ip_v = fpen(h_vis, lam_of(LAM_IP0, phi_vis, PHI_EXP, 0.0));
double f_ip_f = fpen(h_fis, lam_of(LAM_IP0, phi_fis, PHI_EXP, 0.0));
double f_g_p  = fpen(h_par, LAM_G);
double f_g_v  = fpen(h_vis, LAM_G);
double f_g_f  = fpen(h_fis, LAM_G);

// ======================= drug concentrations ==============================
double c_cis = pos(CIS1) / CIS_V1;
double c_pem = pos(PEM1) / PEM_V1;
double c_gem = pos(GEM1) / GEM_V;
double c_niv = pos(NIV1) / MAB_V1;
double c_ipi = pos(IPI1) / MAB_V1;
double plv   = pos(PLV);  if (plv > PLV_MAX) plv = PLV_MAX;
double c_ippl = pos(IPPL) / ((plv / 1000.0 < 0.05) ? 0.05 : plv / 1000.0);

double ro_pd1   = pos(NIVT) / (KD_PD1 + pos(NIVT));
double ro_ctla4 = c_ipi / (KD_CTLA4 + c_ipi);   // the lymph node has no barrier

// ======================= kill rates =======================================
double e_cis = EMAX_CIS * pos(CIST) / (EC50_CIS + pos(CIST));
double e_pem = EMAX_PEM * pos(PEMT) / (EC50_PEM + pos(PEMT));
double e_gem = EMAX_GEM * pos(GEMT) / (EC50_GEM + pos(GEMT));
double syn   = SYN_PEMCIS * sqrt(e_cis * e_pem);
double kill_chem = (e_cis + e_pem + e_gem + syn) * (1.0 - EMT_CHEMO * EMT);

double kill_arg = 0.0;
if (ON_ADI > 0.5) {
  double defc = 1.0 - pos(ARG) / 0.35;
  kill_arg = ADI_EFF * ASS1_NEG * (defc > 0.0 ? defc : 0.0);
}

double pdl1v = pos(PDL1); if (pdl1v > 1.0) pdl1v = 1.0;
// RELEASED kill, not scaled kill: exactly zero with no checkpoint inhibitor on
// board.  Baseline immune surveillance lives inside KG0 (fitted on untreated
// disease); EMAX_IO can therefore only move the checkpoint arms.  What blockade
// releases is bounded by how much PD-L1 was braking -- and PD-L1 is one of the
// things the EMT axis raises, which is where claim 3 comes from.
double release = pdl1v * ro_pd1 + CTLA4_REL * ro_ctla4;
double kill_io = EMAX_IO * pos(TEFF) / (KT50 + pos(TEFF)) * release;

double e_ip = (ON_IP > 0.5) ? 0.055 * c_ippl / (2.0 + c_ippl) : 0.0;
double ttf  = (TTF > 0.5) ? TTF_EFF : 0.0;

double k_par = kill_chem * f_sm_p + kill_io * f_T_p + kill_arg * f_sm_p + e_ip * f_ip_p + ttf;
double k_vis = kill_chem * f_sm_v + kill_io * f_T_v + kill_arg * f_sm_v + e_ip * f_ip_v + ttf;
double k_fis = kill_chem * f_sm_f + kill_io * f_T_f + kill_arg * f_sm_f + e_ip * f_ip_f + 0.4 * ttf;

// ======================= PK ================================================
dxdt_CIS1 = -CIS_CL / CIS_V1 * CIS1 - CIS_Q / CIS_V1 * CIS1 + CIS_Q / CIS_V2 * CIS2;
dxdt_CIS2 =  CIS_Q / CIS_V1 * CIS1 - CIS_Q / CIS_V2 * CIS2;
dxdt_CIST =  CIS_KIN * c_cis - CIS_KOUT * CIST;

double pem_cl = PEM_CL0 * ((GFR > 15.0 ? GFR : 15.0) / GFR0);
dxdt_PEM1 = -pem_cl / PEM_V1 * PEM1 - PEM_Q / PEM_V1 * PEM1 + PEM_Q / PEM_V2 * PEM2;
dxdt_PEM2 =  PEM_Q / PEM_V1 * PEM1 - PEM_Q / PEM_V2 * PEM2;
dxdt_PEMT =  PEM_KIN * c_pem - PEM_KOUT * PEMT;

dxdt_GEM1 = -GEM_CL / GEM_V * GEM1;
dxdt_GEMT =  GEM_KIN * c_gem - GEM_KOUT * GEMT;

dxdt_BEV1 = -MAB_CL / MAB_V1 * BEV1 - MAB_Q / MAB_V1 * BEV1 + MAB_Q / MAB_V2 * BEV2;
dxdt_BEV2 =  MAB_Q / MAB_V1 * BEV1 - MAB_Q / MAB_V2 * BEV2;
dxdt_NIV1 = -MAB_CL / MAB_V1 * NIV1 - MAB_Q / MAB_V1 * NIV1 + MAB_Q / MAB_V2 * NIV2;
dxdt_NIV2 =  MAB_Q / MAB_V1 * NIV1 - MAB_Q / MAB_V2 * NIV2;
dxdt_IPI1 = -MAB_CL / MAB_V1 * IPI1 - MAB_Q / MAB_V1 * IPI1 + MAB_Q / MAB_V2 * IPI2;
dxdt_IPI2 =  MAB_Q / MAB_V1 * IPI1 - MAB_Q / MAB_V2 * IPI2;

// the binding-site barrier lives HERE: antibody entering the rind is gated by
// the IgG penetration fraction of the thickest leaf
dxdt_NIVT = MAB_KIN * c_niv * f_ab_p - MAB_KOUT * NIVT;

// ======================= pleural space =====================================
double q_form  = QFORM0 * (1.0 + (QFORM_VEGF - 1.0) * pos(VEGF) / (1.0 + pos(VEGF)));
double blockf  = h_par / (H50_STOMA + h_par);
double symph   = pos(SYMPH); if (symph > 1.0) symph = 1.0;
double q_drain = QDRAIN_MAX * (1.0 - blockf);
// Pleurodesis removes the SPACE fluid can filter into; it does not disable the
// parietal lymphatics.  The (1 - symphysis) factor therefore belongs on
// FORMATION, not on drainage.
double q_form_eff = q_form * (1.0 - 0.95 * symph);

double ipc = (IPC_DAY >= 0 && SOLVERTIME >= IPC_DAY) ? IPC_RATE : 0.0;
dxdt_PLV = q_form_eff * (1.0 - plv / PLV_MAX)
           - q_drain * plv / (400.0 + plv)
           - ipc * plv / (200.0 + plv);

// intrapleural drug: diluted by the effusion, cleared by residual lymphatics
dxdt_IPPL = -(0.9 + q_drain / 1000.0) * IPPL;
dxdt_IPT  = 0.0;

// pleurodesis needs apposition; it fails in proportion to lung trapping
double E_add = KELAST * phi_vis * h_vis;
double vexp_t = 1.0 / (1.0 + pow(E_add / E50, 2.0));
double talc_on = (TALC_DAY >= 0 && SOLVERTIME >= TALC_DAY) ? 1.0 : 0.0;
dxdt_SYMPH = KSYMPH * talc_on * vexp_t * (1.0 - symph) - KSYMPH_LOSS * symph;

// ======================= rind dynamics =====================================
double kg   = KG0   * (1.0 + EMT_KG   * EMT);
double kcol = KCOL0 * (1.0 + EMT_KCOL * EMT);

double seed_p = KSEED * (tv + tf) * (1.0 - (h_par / 1.5 > 0.95 ? 0.95 : h_par / 1.5));
double seed_v = KSEED * (tp + tf) * (1.0 - (h_vis / 1.5 > 0.95 ? 0.95 : h_vis / 1.5));
double seed_f = KSEED * (tp + tv) * (1.0 - (h_fis / 1.5 > 0.95 ? 0.95 : h_fis / 1.5));

dxdt_TPAR = kg * f_g_p * tp - k_par * tp + seed_p;
dxdt_TVIS = kg * f_g_v * tv - k_vis * tv + seed_v;
dxdt_TFIS = kg * f_g_f * tf - k_fis * tf + seed_f;

dxdt_MPAR = kcol * tp - KMDEG * mp;
dxdt_MVIS = kcol * tv - KMDEG * mv;
dxdt_MFIS = kcol * tf - KMDEG * mf;

dxdt_NPAR = k_par * tp - KNCLR * np_;
dxdt_NVIS = k_vis * tv - KNCLR * nv;
dxdt_NFIS = k_fis * tf - KNCLR * nf;

// ======================= immune ============================================
double antigen = TMB_FACT * (0.15 * T_tot / 300.0 + N_tot / 60.0);
double prime   = KPRIME * antigen * (1.0 + IPI_PRIME * ro_ctla4);
dxdt_PRIME = prime - 0.25 * PRIME;

double supp = 1.0 + 0.9 * pos(TREG) * (1.0 - 0.65 * ro_ctla4)
              + TGFB_SUPP * pos(TGFB) + 0.5 * pos(TAM);
dxdt_TEFF = 0.25 * pos(PRIME) * (1.0 + EMT_TINF * EMT) / supp - KTDEATH * TEFF;
dxdt_TREG = KTREG * (1.0 + 0.8 * pos(TGFB)) - KTREG * TREG - 0.045 * ro_ctla4 * TREG;
dxdt_IFNG = KIFN * 0.30 * pos(TEFF) - KIFN_OUT * IFNG;

double pdl1_ss = PDL1_BASE * (1.0 + EMT_PDL1 * EMT) * (1.0 + PDL1_IND * pos(IFNG));
if (pdl1_ss > 0.95) pdl1_ss = 0.95;
dxdt_PDL1 = 0.5 * (pdl1_ss - PDL1);
dxdt_TAM  = 0.10 * (1.0 + 0.5 * pos(IL6) - TAM);

// ======================= signalling ========================================
double hyp    = 1.0 - f_g_p;
double veg_ss = 0.7 + 2.2 * hyp + 0.5 * T_tot / 400.0;
double bevn   = c_bev / (BEV_KD + c_bev);
dxdt_VEGF = KVEGF_OUT * (veg_ss * (1.0 - 0.90 * bevn) - VEGF);

double tgf_ss = 0.6 + 1.4 * (1.0 + EMT_KCOL * EMT) * T_tot / 400.0;
dxdt_TGFB = KTGFB_OUT * (tgf_ss - TGFB);

double il6_ss = 0.5 + 1.8 * V_tot / 600.0 + 0.6 * N_tot / 60.0;
dxdt_IL6 = KIL6_OUT * (il6_ss - IL6);

// ======================= mechanics =========================================
dxdt_VEXP = KREC * (vexp_t - VEXP);
double vexp = pos(VEXP); if (vexp > 1.0) vexp = 1.0;
double fvc  = FVC_PRED * vexp * (1.0 - PLV_FVC * plv / PLV_MAX);
if (SURG_DAY >= 0 && SOLVERTIME >= SURG_DAY) fvc = fvc * (1.0 - SURG_FVC_HIT);
double dysp_ss = 2.4 * (1.0 - fvc / FVC_PRED) + 1.4 * plv / PLV_MAX;
dxdt_DYSP = KDYSP * (dysp_ss - DYSP);

// ======================= myelosuppression (Friberg) ========================
// driven by the RETAINED species (Pt-DNA adducts, polyglutamates), not by the
// plasma concentrations, whose half-lives are under two hours
double edrug = SLOPE_CIS_ANC * pos(CIST)
             + SLOPE_PEM_ANC * pos(PEMT) * (1.0 - FOLATE_PROT * FOLATE);
if (edrug > 0.98) edrug = 0.98;
double ktr = 4.0 / MTT;
double ancv = (ANC > 0.05) ? ANC : 0.05;
double fb   = pow(CIRC0 / ancv, GAM_FB);
dxdt_PROL = ktr * PROL * (1.0 - edrug) * fb - ktr * PROL;
dxdt_TR1  = ktr * (PROL - TR1);
dxdt_TR2  = ktr * (TR1 - TR2);
dxdt_TR3  = ktr * (TR2 - TR3);
dxdt_ANC  = ktr * (TR3 - ANC);

// ======================= renal, neuro, irAE ================================
dxdt_PTK   = KPT_UP * pos(CIST) - KPT_EL * PTK;
dxdt_GFR   = -KGFR_DAM * pos(PTK) / 100.0 + KGFR_REP * (GFR0 - GFR);
dxdt_NEURO = KNEURO * pos(CIST) - 0.0025 * NEURO;
dxdt_IRAE  = KIRAE_PD1 * ro_pd1 + KIRAE_CTLA4 * ro_ctla4 - KIRAE_OUT * IRAE;

// ======================= cachexia, biomarker, arginine =====================
dxdt_CACHEX = KCACHEX * pos(IL6) / (1.0 + pos(IL6)) - KCACHEX_REC * CACHEX;

// SMRP is SHED from viable cells and cleared renally, so it tracks the quantity
// the drug actually changes -- and it moves EARLIER than thickness does
dxdt_SMRP = KSMRP * T_tot - KSMRP_CL * ((GFR > 10.0 ? GFR : 10.0) / GFR0) * SMRP;
dxdt_ARG  = 0.35 * (1.0 - ARG) - (2.2 * ON_ADI) * ARG;

// ======================= hazard ============================================
double ecog = 0.20 + 1.9 * (1.0 - fvc / FVC_PRED) + 0.050 * pos(CACHEX)
              + 0.30 * (pos(DYSP) > 3.0 ? 3.0 : pos(DYSP));
if (ecog > 4.0) ecog = 4.0;

dxdt_CH = HZ0 * exp(HZ_VOL * log((V_tot > 1.0 ? V_tot : 1.0) / HZ_VREF))
              * exp(HZ_ECOG * ecog)
              * exp(HZ_EMT * EMT)
              * exp(HZ_FVC * (1.0 - fvc / FVC_PRED));

$TABLE
double totp_o = pos(TPAR) + pos(MPAR) + pos(NPAR);
double totv_o = pos(TVIS) + pos(MVIS) + pos(NVIS);
double totf_o = pos(TFIS) + pos(MFIS) + pos(NFIS);
double hpar_o = totp_o / (RHO * A_PAR);
double hvis_o = totv_o / (RHO * A_VIS);
double hfis_o = totf_o / (RHO * A_FIS);

capture h_par_mm  = 10.0 * hpar_o;
capture h_vis_mm  = 10.0 * hvis_o;
capture h_fis_mm  = 10.0 * hfis_o;
capture mRECIST   = 10.0 * (4.0 * hpar_o + 2.0 * hvis_o);
capture Vtumour   = A_PAR * hpar_o + A_VIS * hvis_o + A_FIS * hfis_o;
capture Tviable   = pos(TPAR) + pos(TVIS) + pos(TFIS);
capture phi       = (pos(MPAR)+pos(MVIS)+pos(MFIS)) / ((totp_o+totv_o+totf_o) > 1e-9 ? (totp_o+totv_o+totf_o) : 1e-9);
capture phi_vis   = (totv_o > 1e-9) ? pos(MVIS) / totv_o : 0.0;
capture Eadd      = KELAST * phi_vis * hvis_o;
capture vexp_out  = pos(VEXP) > 1.0 ? 1.0 : pos(VEXP);
capture plv_out   = pos(PLV) > PLV_MAX ? PLV_MAX : pos(PLV);
capture FVC       = FVC_PRED * vexp_out * (1.0 - PLV_FVC * plv_out / PLV_MAX)
                    * ((SURG_DAY >= 0 && TIME >= SURG_DAY) ? (1.0 - SURG_FVC_HIT) : 1.0);
capture FVCpct    = 100.0 * FVC / FVC_PRED;
capture ECOG      = (0.20 + 1.9 * (1.0 - FVC / FVC_PRED) + 0.050 * pos(CACHEX)
                     + 0.30 * (pos(DYSP) > 3.0 ? 3.0 : pos(DYSP)));
capture Surv      = exp(-CH);
capture f_sm      = fpen(hpar_o, lam_of(LAM_SM0, phi, PHI_EXP, 0.0));
capture f_ab      = fpen(hpar_o, lam_of(LAM_AB0, phi, PHI_EXP, 0.0));
capture f_Tcell   = fpen(hpar_o, lam_of(LAM_T0,  phi, PHI_EXP, 0.0));
capture lam_sm_mm = 10.0 * lam_of(LAM_SM0, phi, PHI_EXP, 0.0);
capture lam_ab_mm = 10.0 * lam_of(LAM_AB0, phi, PHI_EXP, 0.0);
capture Ccis      = pos(CIS1) / CIS_V1;
capture Cpem      = pos(PEM1) / PEM_V1;
capture Cniv      = pos(NIV1) / MAB_V1;
capture Cipi      = pos(IPI1) / MAB_V1;
capture Cbev      = pos(BEV1) / MAB_V1;
capture RO_PD1    = pos(NIVT) / (KD_PD1 + pos(NIVT));
capture RO_CTLA4  = Cipi / (KD_CTLA4 + Cipi);
capture stomablk  = hpar_o / (H50_STOMA + hpar_o);
'

mod <- mrgsolve::mcode("mpm", code, end = 1100, delta = 1)


## =============================================================================
##  DOSING REGIMENS
##  Doses are real labelled doses; nothing here is a fitted scaling factor.
## =============================================================================
BSA <- 1.80
WT  <- 72

cis_ev  <- function(n = 6) ev(amt = 75 * BSA,  cmt = "CIS1", ii = 21, addl = n - 1)
carb_ev <- function(n = 6) ev(amt = 0.62 * 5 * (90 + 25), cmt = "CIS1", ii = 21, addl = n - 1)
pem_ev  <- function(n = 6) ev(amt = 500 * BSA, cmt = "PEM1", ii = 21, addl = n - 1)
gem_ev  <- function(n = 6) ev(amt = 1000 * BSA, cmt = "GEM1", ii = 7, addl = 2 * n - 1)
bev_ev  <- function(n = 6, maint = TRUE) ev(amt = 15 * WT, cmt = "BEV1", ii = 21,
                                            addl = if (maint) 34 else n - 1)
niv_ev  <- function(months = 24) ev(amt = 360, cmt = "NIV1", ii = 21,
                                    addl = floor(months * 30.44 / 21) - 1)
ipi_ev  <- function(months = 24) ev(amt = 1 * WT, cmt = "IPI1", ii = 42,
                                    addl = floor(months * 30.44 / 42) - 1)
pembro_ev <- function(n = 6) ev(amt = 200, cmt = "NIV1", ii = 21, addl = n - 1)
ip_ev   <- function(times = c(63.5), amt = 400)
  Reduce(`+`, lapply(times, function(tt) ev(time = tt, amt = amt, cmt = "IPPL")))

## Surgery is a state RESET, not a dose.  It is applied by re-initialising the
## rind compartments at SURG_DAY with mrgsolve's `evid = 8` (replace) records:
## the fissural sanctuary keeps SURG_FIS_SPARE of its mass, and the collagen
## FRACTION of what remains equals the fraction of what was removed.
surgery_records <- function(mod, out_pre, surg_day, resid_cm = 0.10,
                            fis_spare = 0.65) {
  r <- out_pre %>% filter(time == surg_day) %>% slice(1)
  A_PAR <- 1000 * 0.45; A_VIS <- 900 * 0.35; A_FIS <- 400 * 0.50; RHO <- 1.05
  mk <- function(cmtT, cmtM, cmtN, Tm, Mm, Nm, A, spare) {
    tot <- Tm + Mm + Nm
    if (tot <= 0) return(NULL)
    keep <- max(resid_cm * RHO * A, spare * tot)
    f <- min(1, keep / tot)
    ev(time = surg_day, amt = f * Tm, cmt = cmtT, evid = 8) +
      ev(time = surg_day, amt = f * Mm, cmt = cmtM, evid = 8) +
      ev(time = surg_day, amt = f * Nm, cmt = cmtN, evid = 8)
  }
  Reduce(`+`, Filter(Negate(is.null), list(
    mk("TPAR", "MPAR", "NPAR", r$TPAR, r$MPAR, r$NPAR, A_PAR, 0),
    mk("TVIS", "MVIS", "NVIS", r$TVIS, r$MVIS, r$NVIS, A_VIS, 0),
    mk("TFIS", "MFIS", "NFIS", r$TFIS, r$MFIS, r$NFIS, A_FIS, fis_spare))))
}


## =============================================================================
##  SCENARIOS  (20 of them)
## =============================================================================
SCENARIOS <- list(

  ## ---- 1. natural history --------------------------------------------------
  bsc = list(
    label = "Best supportive care",
    param = list(),
    ev    = NULL),

  ## ---- 2-6. first-line systemic therapy -----------------------------------
  cis = list(
    label = "Cisplatin 75 mg/m2 q3w x6 (EMPHACIS control arm)",
    param = list(),
    ev    = cis_ev()),

  pemcis = list(
    label = "Pemetrexed 500 + cisplatin 75 mg/m2 q3w x6 (EMPHACIS)",
    param = list(),
    ev    = cis_ev() + pem_ev()),

  pemcarbo = list(
    label = "Pemetrexed + carboplatin AUC5 q3w x6 (unfit / elderly)",
    param = list(),
    ev    = carb_ev() + pem_ev()),

  gemcis = list(
    label = "Gemcitabine 1000 mg/m2 d1,8 + cisplatin q3w x6",
    param = list(),
    ev    = cis_ev() + gem_ev()),

  pemcisbev = list(
    label = "Pemetrexed + cisplatin + bevacizumab 15 mg/kg, bev maintenance (MAPS)",
    param = list(),
    ev    = cis_ev() + pem_ev() + bev_ev()),

  ## ---- 7-10. immunotherapy -------------------------------------------------
  nivoipi = list(
    label = "Nivolumab 360 mg q3w + ipilimumab 1 mg/kg q6w x2 yr (CheckMate 743)",
    param = list(),
    ev    = niv_ev() + ipi_ev()),

  nivo = list(
    label = "Nivolumab monotherapy (CONFIRM, second line)",
    param = list(),
    ev    = niv_ev()),

  pembrochemo = list(
    label = "Pembrolizumab 200 mg + pemetrexed/cisplatin (IND227 / DREAM3R)",
    param = list(),
    ev    = cis_ev() + pem_ev() + pembro_ev()),

  chemo_then_io = list(
    label = "4 cycles pem/cis then nivolumab + ipilimumab",
    param = list(),
    ev    = cis_ev(4) + pem_ev(4) + ev(time = 84, amt = 360, cmt = "NIV1",
                                       ii = 21, addl = 30) +
            ev(time = 84, amt = 1 * WT, cmt = "IPI1", ii = 42, addl = 14)),

  ## ---- 11-12. biology-directed --------------------------------------------
  adi_pemcis = list(
    label = "ADI-PEG20 + pemetrexed/cisplatin in ASS1-deficient tumour (ATOMIC-Meso)",
    param = list(ON_ADI = 1, ASS1_NEG = 1),
    ev    = cis_ev() + pem_ev()),

  ttf_chemo = list(
    label = "Tumour-treating fields + pemetrexed/cisplatin (STELLAR)",
    param = list(TTF = 1),
    ev    = cis_ev() + pem_ev()),

  ## ---- 13-16. local therapy -------------------------------------------------
  pd_chemo = list(
    label = "Extended pleurectomy/decortication at day 63 + chemotherapy (MARS2 arm)",
    param = list(SURG_DAY = 63, SURG_FVC_HIT = 0, SURG_FIS_SPARE = 0.65),
    ev    = cis_ev() + pem_ev(),
    surgery = TRUE),

  epp_chemo = list(
    label = "Extrapleural pneumonectomy at day 63 + chemotherapy (MARS arm)",
    param = list(SURG_DAY = 63, SURG_FVC_HIT = 0.35, SURG_FIS_SPARE = 0.15),
    ev    = cis_ev() + pem_ev(),
    surgery = TRUE),

  surg_hithoc = list(
    label = "Cytoreduction + hyperthermic intrathoracic chemoperfusion + chemo",
    param = list(SURG_DAY = 63, ON_IP = 1),
    ev    = cis_ev() + pem_ev() + ip_ev(63.5, 400),
    surgery = TRUE),

  ip_only = list(
    label = "Intrapleural chemotherapy WITHOUT cytoreduction (the control case)",
    param = list(ON_IP = 1),
    ev    = ip_ev(c(7, 28, 49), 400)),

  ## ---- 17-19. pleural management -------------------------------------------
  talc = list(
    label = "Talc pleurodesis day 7, no systemic therapy",
    param = list(TALC_DAY = 7),
    ev    = NULL),

  ipc = list(
    label = "Indwelling pleural catheter day 7, no systemic therapy",
    param = list(IPC_DAY = 7),
    ev    = NULL),

  talc_chemo = list(
    label = "Talc pleurodesis day 7 + pemetrexed/cisplatin",
    param = list(TALC_DAY = 7),
    ev    = cis_ev() + pem_ev()),

  ## ---- 20. toxicity stress test --------------------------------------------
  pemcis_nofolate = list(
    label = "Pemetrexed/cisplatin WITHOUT folate + B12 (pre-amendment EMPHACIS)",
    param = list(FOLATE = 0),
    ev    = cis_ev() + pem_ev())
)


run_scenario <- function(key, emt = 0.25, end = 1100, delta = 1, extra = list()) {
  s <- SCENARIOS[[key]]
  if (is.null(s)) stop("unknown scenario: ", key)
  p <- modifyList(modifyList(s$param, list(EMT = emt)), extra)
  m <- mod %>% param(p) %>% update(end = end, delta = delta)
  if (isTRUE(s$surgery)) {
    ## two-pass: run to the day of surgery, read the rind, then re-run with the
    ## resection applied as replacement records
    pre <- m %>% mrgsim_df(events = if (is.null(s$ev)) ev() else s$ev,
                           end = p$SURG_DAY, delta = 1)
    sev <- surgery_records(m, pre, p$SURG_DAY, p$SURG_RESID %||% 0.10,
                           p$SURG_FIS_SPARE %||% 0.65)
    e <- if (is.null(s$ev)) sev else s$ev + sev
  } else {
    e <- s$ev
  }
  out <- if (is.null(e)) m %>% mrgsim_df() else m %>% mrgsim_df(events = e)
  out$scenario <- key
  out$label <- s$label
  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a


## =============================================================================
##  ENDPOINT EXTRACTION
## =============================================================================
median_os <- function(out) {
  i <- which(out$CH >= log(2))[1]
  if (is.na(i) || i < 2) return(NA_real_)
  t0 <- out$time[i - 1]; t1 <- out$time[i]
  c0 <- out$CH[i - 1];   c1 <- out$CH[i]
  (t0 + (log(2) - c0) / (c1 - c0) * (t1 - t0)) / 30.44
}

best_response <- function(out, months = 12) {
  b <- out$mRECIST[1]
  100 * (min(out$mRECIST[out$time <= months * 30.44]) - b) / b
}

pfs_months <- function(out) {
  nad <- cummin(out$mRECIST)
  i <- which(out$mRECIST >= 1.20 * nad & (out$mRECIST - nad) >= 5)[1]
  if (is.na(i)) return(max(out$time) / 30.44)
  out$time[i] / 30.44
}

endpoint_table <- function(keys = names(SCENARIOS), emt = 0.25) {
  do.call(rbind, lapply(keys, function(k) {
    o <- run_scenario(k, emt = emt)
    data.frame(
      scenario   = k,
      label      = SCENARIOS[[k]]$label,
      medOS_mo   = round(median_os(o), 1),
      PFS_mo     = round(pfs_months(o), 1),
      bestmRECIST= round(best_response(o), 1),
      viable_d180= round(100 * (o$Tviable[o$time == 180] - o$Tviable[1]) / o$Tviable[1], 1),
      phi_d180   = round(o$phi[o$time == 180], 3),
      FVC_d180   = round(o$FVC[o$time == 180], 2),
      ANCnadir   = round(min(o$ANC), 2),
      GFR_end    = round(min(o$GFR), 0),
      stringsAsFactors = FALSE)
  }))
}


## =============================================================================
##  THE HISTOLOGY-AXIS CROSSOVER
##  The only pre-registered prediction in this model: the SIGN of
##  OS(IO) - OS(chemo) flips somewhere on the epithelioid-to-sarcomatoid axis.
## =============================================================================
crossover_curve <- function(xs = seq(0, 1, by = 0.05)) {
  do.call(rbind, lapply(xs, function(x) {
    oc <- run_scenario("pemcis",  emt = x)
    oi <- run_scenario("nivoipi", emt = x)
    data.frame(EMT = x,
               OS_chemo = median_os(oc),
               OS_io    = median_os(oi),
               delta    = median_os(oi) - median_os(oc),
               median_ratio = median_os(oc) / median_os(oi))
  }))
}


## =============================================================================
##  GEOMETRY: what a 30% response actually means
## =============================================================================
geometry_note <- function() {
  A_par <- 1000 * 0.45
  cat(sprintf("dV/dh on the parietal leaf   : %.0f mL per mm of rind\n", A_par / 10))
  cat(sprintf("RIND   30%% thickness drop     : %.1f%% volume kill\n", 30))
  cat(sprintf("SPHERE 30%% diameter drop      : %.1f%% volume kill\n",
              100 * (1 - 0.7^3)))
  cat(sprintf("ratio of implied cytoreduction: %.2f x\n", (1 - 0.7^3) / 0.30))
}

penetration_table <- function(phis = c(0.25, 0.55), hs_mm = c(1, 2, 3, 5, 8, 12, 20)) {
  fp <- function(h, lam) if (h <= 0) 1 else (lam / h) * (1 - exp(-h / lam))
  do.call(rbind, lapply(phis, function(phi) {
    do.call(rbind, lapply(hs_mm, function(h_mm) {
      h <- h_mm / 10
      ls <- 0.25 * (1 - phi)^1.5; la <- 0.060 * (1 - phi)^1.5; lt <- 0.050 * (1 - phi)^1.5
      data.frame(h_mm = h_mm, phi = phi,
                 lam_sm_mm = round(10 * ls, 2), f_sm = round(fp(h, ls), 3),
                 lam_ab_mm = round(10 * la, 3), f_ab = round(fp(h, la), 3),
                 f_Tcell   = round(fp(h, lt), 3))
    }))
  }))
}


## =============================================================================
##  PLOTS
## =============================================================================
plot_rind <- function(keys = c("bsc", "pemcis", "pemcisbev", "nivoipi"), emt = 0.25) {
  d <- do.call(rbind, lapply(keys, run_scenario, emt = emt))
  d %>%
    select(time, scenario, mRECIST, Tviable, phi, FVC) %>%
    pivot_longer(-c(time, scenario)) %>%
    ggplot(aes(time, value, colour = scenario)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~name, scales = "free_y") +
    labs(x = "days", y = NULL,
         title = "Viable cell mass falls much further than measured thickness",
         subtitle = "the gap is collagen (phi) that does not resolve") +
    theme_bw()
}

plot_crossover <- function() {
  cc <- crossover_curve()
  ggplot(cc, aes(EMT)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(aes(y = delta), linewidth = 1) +
    labs(x = "EMT axis  (0 epithelioid  -  1 sarcomatoid)",
         y = "median OS (IO) - median OS (chemo), months",
         title = "The sign of the treatment comparison is an OUTPUT of one parameter",
         subtitle = "CheckMate 743: epithelioid HR 0.86, non-epithelioid HR 0.46") +
    theme_bw()
}


## =============================================================================
##  DEMO
## =============================================================================
if (interactive() && !isTRUE(getOption("mpm.no_demo"))) {
  geometry_note()
  print(penetration_table())
  print(endpoint_table())
  print(crossover_curve(seq(0, 1, 0.1)))
  print(plot_rind())
  print(plot_crossover())
}
