# =============================================================================
#  DENGUE / SEVERE DENGUE (DSS) -- QSP MODEL FOR mrgsolve
# =============================================================================
#  45 ODEs.  Antibody-dependent enhancement, target-cell-limited replication,
#  a two-bed Starling system (systemic + serosal), baroreflex haemodynamics,
#  immune thrombocytopenia, shock liver, fever, and eighteen therapeutic
#  scenarios including the WHO fluid ladder.
#
#  THE ORGANISING IDEA
#  -------------------
#  Dengue is written here as ONE antibody response read through TWO channels
#  of opposite sign:
#
#    channel 1 (protective)   neutralising IgG clears virions, viraemia falls,
#                             interferon falls, the fever breaks;
#    channel 2 (destructive)  the SAME IgG binds NS1 and virion, forming immune
#                             complexes that activate complement and mast
#                             cells and arm cross-reactive T cells; the
#                             glycocalyx is stripped; the reflection
#                             coefficient sigma falls; plasma leaves the
#                             circulation.
#
#  Nothing in the code links the temperature equation to the leak equations
#  after the antibody rise -- temperature is computed from hypothalamic PGE2
#  driven by interferon, the leak from a glycocalyx driven by NS1, TNF, VEGF
#  and chymase -- yet the clinical nadir and the fever break land within a day
#  of each other.  That is the model's central claim and it is a computed
#  result, not a rule.
#
#  THE SECOND IDEA: the leak is self-limiting only through shock.
#      J = Kf * [ (Pc - Pi) - sigma*(PIp - PIi) ]
#  and the only term that falls fast enough to switch filtration off is Pc,
#  which falls because the plasma volume that feeds it has gone.  Resuscitation
#  therefore restores perfusion by restarting the leak, which is why the fluid
#  dose-response has an interior optimum instead of being monotone.
#
#  VERIFICATION
#  ------------
#  Every equation below was independently re-implemented in Python/scipy
#  (denv_reference_model.py) and solved with LSODA.  All quoted numbers are
#  that solver's output.  The re-implementation found and fixed five defects
#  in earlier drafts of this model; they are listed in README.md.
#
#  CALIBRATION ANCHORS (each is a published measurement, not a fitted value)
#  ------------------------------------------------------------------------
#   * Salje 2018 Nature 557:719  -- hospitalised-dengue risk is raised at
#     pre-infection titres 1:21-1:80 and suppressed above 1:1280.  The model's
#     entry factor peaks at 1:57 (E = 12.8) and crosses 1.0 at 1:696.
#   * Kliks 1988 / Chau 2009     -- infant DHF peaks at 6-8 months of age.
#     The model predicts 6.3 months from a cord titre and an IgG half-life.
#   * Guyton, Textbook of Medical Physiology -- Pc 17.3, Pi -3.0, Pi_p 27.4,
#     Pi_i 8.0 mmHg and 140 mL/h baseline lymph flow.  Kf and Kfs are then
#     NOT free: they are whatever makes the uninfected host stationary.
#   * Miserocchi 1997 -- parietal pleural lymphatic absorption ceiling
#     ~0.65 mL/kg/h.  This is the reason effusions accumulate at all.
#   * WHO 2009 haemoconcentration criterion -- a 20 % haematocrit rise.
#     The model computes what that costs: 800 mL, 28.6 % of the compartment.
#   * Wills 2005 NEJM 353:877, Lye 2017 Lancet 389:1611 (AAPT),
#     Tam 2012 Clin Infect Dis 55:1216, Nguyen 2013 J Infect Dis 207:1442
#     (balapiravir) -- four negative trials the model reproduces.
#
#  Usage:
#    library(mrgsolve); library(dplyr); library(ggplot2)
#    source("denv_mrgsolve_model.R")      # defines `mod` via mcode_cache
#    out <- run_scenario("S02_secondary_untreated")
#    compare_table(run_all())
#
#  UNITS: time in hours from FEVER ONSET (illness day 0).  The 4-7 day
#  incubation is deliberately outside the model; every scenario starts from
#  the same viraemia so that the antibody state, and not a different starting
#  point, is the only thing that differs between arms.
# =============================================================================

library(mrgsolve)

code <- '
$PROB
# DENGUE / SEVERE DENGUE QSP MODEL
# 45 ODEs | 18 scenarios | ADE bell curve | two-bed Starling | WHO fluid ladder

$PARAM @annotated
// ---------- host -----------------------------------------------------------
WT       :  70.0   : body weight (kg)
VP0      : 2800.0  : baseline plasma volume (mL)
VI0      : 10500.0 : baseline interstitial volume (mL)
RBC0     : 2000.0  : red-cell volume (mL)
TP_P0    :  7.30   : baseline plasma total protein (g/dL)
TP_I0    :  2.00   : baseline interstitial total protein (g/dL)

// ---------- virology -------------------------------------------------------
TGT0     : 2.0e4   : productively infectable FcgR+ target cells (/mL)
dT       : 4.1667e-4 : target-cell turnover (1/h)
kECL     : 0.25    : eclipse exit rate (1/h)
dINF0    : 0.041667 : infected-cell death rate (1/h)
kCTLkill : 0.10    : CTL killing per unit CTL (1/h)
BETA0    : 1.19e-9 : FcgR-independent entry rate constant (mL/copy/h)
PPROD0   : 1050.0  : IFN-free burst rate (copies/cell/h)
PINT     : 3.00    : intrinsic-ADE gain on burst size (-)
IC50IFN  : 400.0   : IFN halving viral production (pg/mL)
KVSAT    : 1.5e7   : FcgR uptake saturation (copies/mL)
cV       : 0.145833 : virion clearance (1/h)
kNeutCl  : 2.0e-4  : opsonic virion clearance per titre unit (1/h)
pNS1     : 0.0155  : NS1 secretion ((ng/mL)/(cell/mL)/h)
cNS1     : 0.098   : NS1 elimination (1/h)
kNS1ab   : 2.2e-4  : immune-complex NS1 clearance per titre unit (1/h)

// ---------- antibody-dependent enhancement --------------------------------
FCROSS   : 0.12    : neutralising weight of heterotypic IgG (-)
AVIDN    : 3.00    : avidity advantage of homotypic de-novo IgG (-)
NT50     : 30.0    : reciprocal titre for 50 percent neutralisation (-)
HNEUT    : 2.50    : Hill coefficient for multi-hit neutralisation (-)
KOPS     : 4.00    : reciprocal titre for 50 percent FcgR ligation (-)
PHI      : 14.0    : FcgR amplification of monocyte entry (-)
KINT     : 0.85    : intrinsic-ADE suppression of interferon (-)

// ---------- adaptive response ---------------------------------------------
kPBL     : 0.055   : naive B-cell recruitment (1/h)
kPBLM    : 0.185   : memory B-cell recruitment (1/h)
KVPBL    : 3.0e6   : half-max antigen drive (copies/mL)
dPBL     : 0.030   : plasmablast loss (1/h)
kMAT     : 0.030   : B cell to antibody-secreting cell maturation (1/h)
dASC     : 0.021   : ASC loss (1/h)
kAB      : 3.30    : IgG production per ASC (titre/h)
dAB      : 1.375e-3 : IgG catabolism, t-half 21 d (1/h)
dABH     : 6.716e-4 : heterotypic IgG waning, t-half 43 d (1/h)
kCTL     : 0.058   : CD8 activation (1/h)
KICTL    : 40.0    : half-max infected-cell drive (cells/mL)
dCTL     : 0.025   : CD8 contraction (1/h)
MEMT     : 2.40    : memory T-cell amplification (-)
XAVID    : 0.45    : cross-reactive T-cell killing efficiency (-)
XTNF     : 3.10    : cross-reactive T-cell TNF multiple (-)

// ---------- mediators ------------------------------------------------------
kIFN     : 0.0085  : interferon production ((pg/mL)/(cell/mL)/h)
dIFN     : 0.10    : interferon elimination (1/h)
kTNFm    : 34.0    : TNF from immune complexes ((pg/mL)/h)
kTNFt    : 3.20    : TNF from activated T cells ((pg/mL)/h)
dTNF     : 0.29    : TNF elimination (1/h)
kIL10    : 620.0   : IL-10 production ((pg/mL)/h)
dIL10    : 0.21    : IL-10 elimination (1/h)
kVEGF    : 78.0    : VEGF-A release ((pg/mL)/h)
dVEGF    : 0.14    : VEGF-A elimination (1/h)
kCHYM    : 13.5    : mast-cell chymase release ((ng/mL)/h)
dCHYM    : 0.11    : chymase elimination (1/h)

// ---------- endothelium ----------------------------------------------------
kGSYN    : 0.0065  : glycocalyx regeneration, t-half 107 h (1/h)
kGDEG    : 0.082   : maximal glycocalyx stripping (1/h)
wNS1     : 0.20    : weight of direct NS1 damage (-)
KNS1     : 1100.0  : NS1 half-effect (ng/mL)
wTNF     : 0.36    : weight of TNF damage (-)
KTNF     : 110.0   : TNF half-effect (pg/mL)
wVEG     : 0.12    : weight of VEGF damage (-)
KVEG     : 520.0   : VEGF half-effect (pg/mL)
wCHY     : 0.32    : weight of chymase damage (-)
KCHY     : 11.0    : chymase half-effect (ng/mL)
SIG0     : 0.97    : baseline reflection coefficient (-)
SIGMIN   : 0.55    : reflection coefficient floor (-)
NSIG     : 1.30    : exponent linking glycocalyx to sigma (-)
ALPHLP   : 1.20    : maximal rise in hydraulic conductance (-)
tauLP    : 3.00    : conductance time constant (h)

// ---------- Starling: systemic bed -----------------------------------------
KF0      : 91.77324 : systemic filtration coefficient (mL/h/mmHg)
PI0      : -3.00   : baseline interstitial hydrostatic pressure (mmHg)
PIMAX    :  2.60   : interstitial pressure ceiling (mmHg)
VSC      : 1600.0  : interstitial compliance scale (mL)
JL0      : 140.0   : baseline lymph flow (mL/h)
kLYMPH   : 154.0   : lymph recruitment (mL/h/mmHg)
JLMAX    : 1100.0  : lymph pump ceiling (mL/h)
PSPROT   : 47.0    : diffusive protein permeability-surface (mL/h)
PRSYN0   : 1.04    : plasma protein synthesis (g/h)
PRSUP    : 0.55    : acute-phase suppression of protein synthesis (-)
KPRSUP   : 60.0    : TNF half-effect on synthesis (pg/mL)
kPRDEG   : 0.00509 : plasma protein catabolism (1/h)

// ---------- Starling: serosal bed ------------------------------------------
KFS      : 13.88890 : serosal filtration coefficient (mL/h/mmHg)
PSER0    : -5.00   : resting pleural/peritoneal pressure (mmHg)
CSER     : 150.0   : serosal compliance (mL/mmHg)
TP_S0    :  1.50   : resting serosal protein (g/dL)
TP_SMAX  :  3.20   : protein of an exudative dengue effusion (g/dL)
DSER0    : 20.0    : baseline serosal turnover (mL/h)
VSERB    : 20.0    : resting serosal fluid volume (mL)
DSERMAX  : 46.0    : parietal lymphatic absorption ceiling (mL/h)
SIGCOL_F : 1.02    : colloid sieving advantage over native protein (-)
COLLCOP  : 5.05    : colloid oncotic pressure per g/dL (mmHg)
COLLEL   : 0.010   : colloid metabolic loss (1/h)

// ---------- haemodynamics --------------------------------------------------
SVMAX    : 140.0   : maximal stroke volume (mL)
VPUN     : 1800.0  : unstressed plasma volume (mL)
KSV      : 1000.0  : Frank-Starling half-constant (mL)
kINO     : 0.35    : sympathetic inotropy (-)
kMYOC    : 0.28    : myocardial depression at saturating TNF (-)
HR0      : 72.0    : baseline heart rate (bpm)
kHRB     : 46.0    : baroreflex chronotropy (bpm)
kHRT     : 7.50    : fever chronotropy (bpm/degC)
tauHR    : 0.10    : heart-rate time constant (h)
SVR0     : 17.46   : baseline systemic vascular resistance (mmHg*min/L)
BARGAIN  : 5.20    : baroreflex gain (-)
SVRMAX   : 2.55    : maximal SVR multiplier (-)
tauSVR   : 0.16    : SVR time constant (h)
MAPSET   : 88.0    : baroreflex set point (mmHg)
kVASOPL  : 0.55    : cytokine erosion of the SVR ceiling (-)
KVASOPL  : 140.0   : TNF half-effect on vasoplegia (pg/mL)
CART     : 1.82    : arterial compliance (mL/mmHg)
PV0      : 8.00    : baseline central venous pressure (mmHg)
CVEN     : 200.0   : systemic venous compliance (mL/mmHg)
FPLAS    : 0.583   : plasma fraction of blood volume (-)
PVMIN    : 1.00    : floor on central venous pressure (mmHg)
RRATIO0  : 7.60    : baseline pre/post-capillary resistance ratio (-)
kRRAT    : 2.50    : sympathetic rise in that ratio (-)

// ---------- renal and intake -----------------------------------------------
UO0      : 62.0    : baseline urine output (mL/h)
KMAPU    : 62.0    : renal autoregulation breakpoint (mmHg)
NMAPU    : 3.00    : renal pressure exponent (-)
kUOVOL   : 4.00    : volume diuresis exponent (-)
INSENS   : 43.0    : insensible loss (mL/h)
ORAL0    : 105.0   : baseline oral intake (mL/h)
ANOREX   : 0.12    : fall in oral intake at peak fever (-)

// ---------- haematology ----------------------------------------------------
kPLTP    : 1.05    : marrow platelet output (1e9/L/h)
dPLT     : 0.00417 : platelet senescence (1/h)
kPLTI    : 0.265   : immune platelet destruction (1/h)
kPLTC    : 0.030   : platelet consumption on damaged endothelium (1/h)
kMKR     : 0.011   : marrow recovery (1/h)
kMKS     : 0.055   : marrow suppression by virus (1/h)
KVMK     : 2.0e6   : half-max marrow suppression (copies/mL)
kAPLT    : 0.0125  : anti-platelet antibody generation (1/h)
dAPLT    : 0.011   : anti-platelet antibody loss (1/h)
kFIBC    : 0.020   : fibrinogen consumption (1/h)
kFIBS    : 0.0090  : fibrinogen synthesis (1/h)
FIBAPR   : 1.45    : acute-phase fibrinogen multiple (-)
kWBCS    : 0.075   : leukopenia rate (1/h)
kWBCR    : 0.030   : leukocyte recovery (1/h)

// ---------- liver ----------------------------------------------------------
kHEPV    : 0.0042  : maximal virus-driven hepatocyte injury (1/h)
KVHEP    : 1.0e8   : half-max viral hepatotoxicity (copies/mL)
kHEPS    : 0.020   : shock-liver injury (1/h)
kHEPAP   : 1.1e-4  : paracetamol hepatotoxicity (1/h per mg/L)
APAPTH   : 22.0    : paracetamol hepatotoxic threshold (mg/L)
kHEPR    : 0.0055  : hepatocyte regeneration (1/h)
kAST     : 7400.0  : AST release per unit injury rate (U/L)
kALT     : 3100.0  : ALT release per unit injury rate (U/L)
dAST     : 0.041   : AST elimination (1/h)
dALT     : 0.0154  : ALT elimination (1/h)

// ---------- fever ----------------------------------------------------------
kPGE     : 0.155   : hypothalamic PGE2 production (1/h)
KPYR     : 62.0    : half-max pyrogen (pg/mL)
dPGE     : 0.115   : PGE2 clearance (1/h)
TSETMAX  : 3.10    : maximal set-point rise (degC)
KPGE     : 0.52    : PGE2 half-effect on set point (-)
tauT     : 1.60    : thermal time constant (h)
APAPIC50 : 9.00    : paracetamol COX IC50 (mg/L)

// ---------- drugs ----------------------------------------------------------
kaAV     : 1.10    : antiviral absorption (1/h)
kelAV    : 0.077   : antiviral elimination (1/h)
VdAV     : 92.0    : antiviral volume of distribution (L)
EMAXAV   : 0.965   : maximal block of viral production (-)
EC50AV   : 0.42    : antiviral EC50 (mg/L)
kelAPAP  : 0.315   : paracetamol elimination (1/h)
VdAPAP   : 49.0    : paracetamol volume of distribution (L)
kelSTER  : 0.231   : methylprednisolone elimination (1/h)
VdSTER   : 84.0    : methylprednisolone volume of distribution (L)
STEREM   : 0.72    : steroid suppression of TNF (-)
STEREC   : 0.55    : steroid EC50 (mg/L)
STERVIR  : 0.70    : steroid suppression of CD8 response (-)
STERAB   : 0.45    : steroid suppression of ASC maturation (-)

// ---------- scenario switches (set per scenario, not per patient) ----------
ABH0     : 0.0     : pre-existing heterotypic IgG titre (-)
V0       : 4.0e5   : viraemia at fever onset (copies/mL)
RCRYS    : 0.0     : crystalloid infusion rate (mL/h)
RCOLL    : 0.0     : colloid infusion rate (mL/h)
RALB     : 0.0     : albumin infusion rate (g/h)
RORAL    : 0.0     : supplemental oral rehydration (mL/h)
NS1MAB   : 0.0     : anti-NS1 monoclonal clearance (1/h)

$CMT @annotated
TGT  : susceptible FcgR+ target cells (cells/mL)
ECL  : eclipse-phase infected cells (cells/mL)
INF  : productively infected cells (cells/mL)
V    : plasma viraemia (copies/mL)
NS1  : circulating NS1 antigen (ng/mL)
ABH  : pre-existing heterotypic IgG (reciprocal titre)
ABN  : de-novo neutralising IgG (reciprocal titre)
PBL  : plasmablasts (relative)
ASC  : antibody-secreting cells (relative)
CTL  : activated CD8 T cells (relative)
IFN  : type-I interferon (pg/mL)
TNF  : TNF-alpha (pg/mL)
IL10 : IL-10 (pg/mL)
VEGF : free VEGF-A (pg/mL)
CHYM : mast-cell chymase (ng/mL)
GLX  : endothelial glycocalyx integrity (0-1)
LPX  : hydraulic conductance multiplier (-)
VP   : plasma volume (mL)
VI   : interstitial volume (mL)
PRP  : plasma protein mass (g)
PRI  : interstitial protein mass (g)
VSER : serosal fluid volume (mL)
RBCV : red-cell volume (mL)
PLT  : platelets (1e9/L)
MKC  : marrow platelet output capacity (relative)
APLT : anti-platelet antibody (relative)
FIB  : fibrinogen (mg/dL)
WBC  : leukocytes (1e9/L)
HEP  : viable hepatocyte fraction (0-1)
AST  : aspartate aminotransferase (U/L)
ALT  : alanine aminotransferase (U/L)
PGE  : hypothalamic PGE2 (relative)
TEMP : core temperature (degC)
SVRR : systemic vascular resistance multiplier (-)
HR   : heart rate (bpm)
LAC  : arterial lactate (mmol/L)
COLL : synthetic colloid mass in plasma (g)
AVD  : antiviral gut depot (mg)
AVC  : antiviral plasma concentration (mg/L)
APAP : paracetamol concentration (mg/L)
STER : methylprednisolone concentration (mg/L)
CIN  : cumulative fluid in (mL)
COUT : cumulative urine out (mL)
NAUC : cumulative NS1 exposure (ng/mL*h)
ICAUC: cumulative immune-complex exposure (h)

$GLOBAL
#define hill2(x, k) (((x)*(x)) / ((x)*(x) + (k)*(k)))
// Landis-Pappenheimer on TOTAL plasma protein (g/dL).  Applying this to
// albumin alone -- a mistake that survived one draft of this model -- puts
// baseline net filtration at 13 mmHg instead of 1.5 and makes every
// downstream number meaningless.
#define copP(c) (2.1*(c) + 0.16*(c)*(c) + 0.009*(c)*(c)*(c))
// Interstitial colloid osmotic pressure, calibrated to the Guyton value of 8 mmHg at
// 2.0 g/dL; interstitial fluid is not plasma and does not share its curve.
#define copI(c) (3.6*(c) + 0.20*(c)*(c))

$MAIN
if (NEWIND < 2) {
  // ---- t = 0 is FEVER ONSET.  ECL and INF are seeded at the values
  // ---- consistent with pre-symptomatic exponential growth at lam, so the
  // ---- trajectory is continuous rather than starting with a kink.
  double lam = 0.080;
  TGT_0  = TGT0;
  V_0    = V0;
  INF_0  = (lam + cV) * V0 / PPROD0;
  ECL_0  = (lam + dINF0) * INF_0 / kECL;
  TGT_0  = TGT0 - INF_0 - ECL_0;
  ABH_0  = ABH0;
  IFN_0  = (V0 > 0.0) ? 28.0  : 0.0;
  PGE_0  = (V0 > 0.0) ? 0.49  : 0.0;
  TEMP_0 = (V0 > 0.0) ? 38.4  : 37.0;
  GLX_0  = 1.0;
  LPX_0  = 1.0;
  VP_0   = VP0;
  VI_0   = VI0;
  PRP_0  = TP_P0 / 100.0 * VP0;
  PRI_0  = TP_I0 / 100.0 * VI0;
  VSER_0 = VSERB;
  RBCV_0 = RBC0;
  PLT_0  = 250.0;
  MKC_0  = 1.0;
  FIB_0  = 300.0;
  WBC_0  = 6.5;
  HEP_0  = 1.0;
  AST_0  = 24.0;
  ALT_0  = 22.0;
  SVRR_0 = 1.0;
  HR_0   = HR0;
  LAC_0  = 1.0;
}

$ODE
// ===========================================================================
// 1.  ANTIBODY-DEPENDENT ENHANCEMENT
// ===========================================================================
// Two occupancies are read off the SAME antibody pool.  Neutralisation needs
// a high, multi-hit stoichiometric occupancy and is driven by the
// neutralising-equivalent titre; opsonisation needs one or two IgG per virion
// to ligate FcgR and is driven by the total titre.  Writing O for the
// opsonised-but-not-neutralised fraction,
//     E = (1 - N - O) + PHI*O = (1 - N)*[1 + (PHI-1)*Aopz/(Aopz+KOPS)]
// which is exactly 1 in a naive host, peaks at an intermediate titre and
// falls back through 1 once neutralisation takes over.  The bell shape is the
// product of a rising saturating term and a falling Hill term -- it is not
// imposed on the model, it comes out of it.
double aeff = AVIDN * ABN + FCROSS * ABH;
double aopz = ABN + ABH;
double Nneut = (aeff > 0.0)
  ? pow(aeff, HNEUT) / (pow(aeff, HNEUT) + pow(NT50, HNEUT)) : 0.0;
double xops  = aopz / (aopz + KOPS);
double Oops  = (1.0 - Nneut) * xops;
double Eade  = (1.0 - Nneut) * (1.0 + (PHI - 1.0) * xops);
double memory = (ABH > 1.0) ? 1.0 : 0.0;

// ===========================================================================
// 2.  VIROLOGY -- target-cell limited, with FcgR uptake saturation
// ===========================================================================
double av_block = EMAXAV * AVC / (AVC + EC50AV);
double p_eff = PPROD0 * (1.0 + PINT * Oops) / (1.0 + IFN / IC50IFN)
             * (1.0 - av_block);
double infect = BETA0 * Eade * V * TGT / (1.0 + V / KVSAT);
double dINFrate = dINF0 + kCTLkill * CTL * (memory > 0.5 ? XAVID : 1.0);

dxdt_TGT = dT * (TGT0 - TGT) - infect;
dxdt_ECL = infect - kECL * ECL;
dxdt_INF = kECL * ECL - dINFrate * INF;
dxdt_V   = p_eff * INF - cV * V - kNeutCl * (ABN + 0.35 * ABH) * V;
dxdt_NS1 = pNS1 * (1.0 - av_block) * INF - (cNS1 + NS1MAB) * NS1
         - kNS1ab * ABN * NS1;

// ===========================================================================
// 3.  ADAPTIVE RESPONSE
// ===========================================================================
double antigen  = V / (V + KVPBL);
double ster_imm = 1.0 - STERVIR * STER / (STER + STEREC);
double ster_ab  = 1.0 - STERAB  * STER / (STER + STEREC);

dxdt_PBL = (kPBL + kPBLM * memory) * antigen - (dPBL + kMAT) * PBL;
dxdt_ASC = kMAT * PBL * ster_ab - dASC * ASC;
dxdt_ABN = kAB * ASC - dAB * ABN;
dxdt_ABH = -dABH * ABH;
dxdt_CTL = kCTL * INF / (INF + KICTL) * (1.0 + MEMT * memory) * ster_imm
         - dCTL * CTL;

// ===========================================================================
// 4.  IMMUNE COMPLEXES AND MEDIATORS
// ===========================================================================
// IC is a product of a RISING function of antibody and a FALLING function of
// antigen, so it peaks strictly between the antibody rise and the antigen
// fall.  That is why the mediator storm sits at the fever break rather than
// at the viraemia peak.
double fAb = ABN / (ABN + 400.0);
double IC  = fAb * (NS1 / (NS1 + 240.0) + 0.55 * V / (V + 4.0e6));
double ster_tnf = 1.0 - STEREM * STER / (STER + STEREC);

dxdt_IFN  = kIFN * INF / (1.0 + KINT * Oops) - dIFN * IFN;
dxdt_TNF  = (kTNFm * IC + kTNFt * CTL * (memory > 0.5 ? XTNF : 1.0)) * ster_tnf
          - dTNF * TNF;
dxdt_IL10 = kIL10 * (IC + 0.02 * CTL) - dIL10 * IL10;
dxdt_VEGF = kVEGF * (IC + 0.006 * TNF) - dVEGF * VEGF;
dxdt_CHYM = kCHYM * IC - dCHYM * CHYM;

// ===========================================================================
// 5.  ENDOTHELIUM -- the glycocalyx is a STRUCTURE, so it integrates
// ===========================================================================
// Each driver enters as a Hill-2 term: the glycocalyx tolerates a mediator
// load and then gives way.  That threshold is what separates dengue fever
// from dengue haemorrhagic fever along a continuous driver axis, and the slow
// time constant (20-30 h) is what puts the leak at illness day 4-6 instead of
// at the cytokine peak on day 2.
double damage = wNS1 * hill2(NS1, KNS1) + wTNF * hill2(TNF, KTNF)
              + wVEG * hill2(VEGF, KVEG) + wCHY * hill2(CHYM, KCHY);
dxdt_GLX = kGSYN * (1.0 - GLX) - kGDEG * GLX * damage;
dxdt_LPX = (1.0 + ALPHLP * (1.0 - GLX) - LPX) / tauLP;

double sigma     = SIGMIN + (SIG0 - SIGMIN) * pow(GLX, NSIG);
double sigma_col = fmin(0.995, sigma * SIGCOL_F);

// ===========================================================================
// 6.  HAEMODYNAMICS
// ===========================================================================
double stressed = fmax(VP - VPUN, 0.0);
double myoc = 1.0 - kMYOC * TNF / (TNF + 150.0);
double SV   = SVMAX * stressed / (KSV + stressed)
            * (1.0 + kINO * (SVRR - 1.0)) * myoc;
double CO   = SV * HR / 1000.0;
double MAP  = fmax(CO * SVR0 * SVRR, 5.0);
double PP   = SV / CART;

double ceiling = fmax(1.0, SVRMAX * (1.0 - kVASOPL * TNF / (TNF + KVASOPL)));
double svr_t = 1.0 + BARGAIN * (MAPSET - MAP) / MAPSET;
svr_t = fmin(fmax(svr_t, 0.85), ceiling);
dxdt_SVRR = (svr_t - SVRR) / tauSVR;
double hr_t = HR0 + kHRB * fmax(SVRR - 1.0, 0.0) + kHRT * fmax(TEMP - 37.0, 0.0);
dxdt_HR = (hr_t - HR) / tauHR;

// Capillary pressure.  Sympathetic arteriolar constriction raises the
// pre/post ratio and so DROPS Pc -- the autotransfusion that is the leak
// only physiological brake, and the reason the untreated patient stops
// leaking exactly when the circulation fails.
double PVcvp  = fmax(PVMIN, PV0 + (VP - VP0) / FPLAS / CVEN);
double rratio = RRATIO0 * (1.0 + kRRAT * fmax(SVRR - 1.0, 0.0) / SVRMAX);
double Pc = (MAP + rratio * PVcvp) / (1.0 + rratio);

// ===========================================================================
// 7.  STARLING EXCHANGE -- two beds with very different drains
// ===========================================================================
double Cp   = PRP / VP * 100.0;
double Ci   = PRI / VI * 100.0;
double Ccol = COLL / VP * 100.0;
double PIp  = copP(Cp) + COLLCOP * Ccol;
double PIi  = copI(Ci);
double Pi   = PIMAX - (PIMAX - PI0) * exp(-fmax(VI - VI0, 0.0) / VSC);

double Jv = KF0 * LPX * ((Pc - Pi) - sigma * (PIp - PIi));
double Jl = fmax(0.0, fmin(JL0 + kLYMPH * (Pi - PI0), JLMAX));

// The serosal bed.  Same sigma, but the resting pressure is -5 mmHg and the
// drain is parietal pleural and diaphragmatic lymphatics whose ceiling is two
// orders of magnitude below systemic lymph flow.  Fluid therefore accumulates
// here and is stopped only by the pressure it generates -- which is why
// dengue produces effusions rather than generalised oedema.
double Pser = PSER0 + fmax(VSER - VSERB, 0.0) / CSER;
double Cser = TP_S0 + (TP_SMAX - TP_S0) * (SIG0 - sigma)
            / fmax(SIG0 - SIGMIN, 1e-9);
double Jser = KFS * LPX * ((Pc - Pser) - sigma * (PIp - copI(Cser)));
double Dser = fmin(DSERMAX, DSER0 * sqrt(fmax(VSER, 0.0) / VSERB));

double anorexia = 1.0 - ANOREX * fmin(fmax(TEMP - 37.0, 0.0) / 2.5, 1.0);
double oral   = (ORAL0 * anorexia + RORAL) * (MAP < 60.0 ? 0.6 : 1.0);
double insens = INSENS * (1.0 + 0.13 * fmax(TEMP - 37.0, 0.0));
double f_map  = (pow(MAP, NMAPU) / (pow(MAP, NMAPU) + pow(KMAPU, NMAPU)))
              / (pow(MAPSET, NMAPU) / (pow(MAPSET, NMAPU) + pow(KMAPU, NMAPU)));
double f_vol  = fmin(3.0, fmax(0.15, pow(VP / VP0, kUOVOL)));
double uo     = UO0 * f_map * f_vol;

dxdt_VP   = -Jv + Jl - Jser + Dser + RCRYS + RCOLL + oral - uo - insens;
dxdt_VI   = Jv - Jl;
dxdt_VSER = Jser - Dser;

double Js     = (1.0 - sigma) * fmax(Jv, 0.0) * (PRP / VP)
              + PSPROT * (Cp - Ci) / 100.0;
double back   = Jl * (PRI / VI);
double Jsprot = (Jser - Dser) * Cser / 100.0;
double prsyn  = PRSYN0 * (1.0 - PRSUP * TNF / (TNF + KPRSUP));

dxdt_PRP  = prsyn - Js + back + RALB - Jsprot - kPRDEG * PRP;
dxdt_PRI  = Js - back;
dxdt_COLL = RCOLL * 0.06 - (1.0 - sigma_col) * fmax(Jv, 0.0) * (COLL / VP)
          - COLLEL * COLL;
dxdt_RBCV = 0.0;
dxdt_CIN  = RCRYS + RCOLL + oral;
dxdt_COUT = uo;

// ===========================================================================
// 8.  HAEMATOLOGY
// ===========================================================================
dxdt_MKC  = kMKR * (1.0 - MKC) - kMKS * V / (V + KVMK) * MKC;
dxdt_APLT = kAPLT * IC * 100.0 - dAPLT * APLT;
dxdt_PLT  = kPLTP * MKC - dPLT * PLT - kPLTI * APLT * PLT / 250.0
          - kPLTC * (1.0 - GLX) * PLT;
double fib_target = 300.0 * FIBAPR * (1.0 - 0.55 * (1.0 - HEP));
dxdt_FIB  = kFIBS * (fib_target - FIB) - kFIBC * FIB * (1.0 - GLX);
dxdt_WBC  = kWBCR * (6.5 - WBC) - kWBCS * WBC * V / (V + KVMK);

// ===========================================================================
// 9.  LIVER
// ===========================================================================
double perf_def = fmax(0.0, 1.0 - CO / 4.6);
double apap_tox = kHEPAP * fmax(APAP - APAPTH, 0.0);
double inj = (kHEPV * V / (V + KVHEP) + kHEPS * perf_def * perf_def + apap_tox) * HEP;
dxdt_HEP = kHEPR * (1.0 - HEP) - inj;
dxdt_AST = kAST * inj - dAST * (AST - 24.0);
dxdt_ALT = kALT * inj - dALT * (ALT - 22.0);

// ===========================================================================
// 10. FEVER -- the other channel, and it shares no equation with the leak
// ===========================================================================
double pyrogen = IFN + 0.06 * TNF;
double cox = 1.0 / (1.0 + APAP / APAPIC50);
dxdt_PGE = kPGE * pyrogen / (pyrogen + KPYR) * cox - dPGE * PGE;
double tset = 37.0 + TSETMAX * PGE / (PGE + KPGE);
dxdt_TEMP = (tset - TEMP) / tauT;

// ===========================================================================
// 11. OXYGEN DELIVERY, DRUG PK, EXPOSURE INTEGRALS
// ===========================================================================
double Hct = RBCV / (RBCV + VP) * 100.0;
double DO2 = CO * (Hct / 3.0) * 1.34 * 0.97 * 10.0 / 1000.0;
dxdt_LAC = 0.24 + 2.6 * fmax(0.0, 1.0 - DO2 / 0.72) - 0.24 * LAC * fmax(HEP, 0.2);

dxdt_AVD  = -kaAV * AVD;
dxdt_AVC  = kaAV * AVD / VdAV - kelAV * AVC;
dxdt_APAP = -kelAPAP * APAP;
dxdt_STER = -kelSTER * STER;

dxdt_NAUC  = NS1;
dxdt_ICAUC = IC;

$TABLE
double aeff_o = AVIDN * ABN + FCROSS * ABH;
double aopz_o = ABN + ABH;
double Nn_o = (aeff_o > 0.0)
  ? pow(aeff_o, HNEUT) / (pow(aeff_o, HNEUT) + pow(NT50, HNEUT)) : 0.0;
double E_o  = (1.0 - Nn_o) * (1.0 + (PHI - 1.0) * aopz_o / (aopz_o + KOPS));

double sig_o = SIGMIN + (SIG0 - SIGMIN) * pow(GLX, NSIG);
double str_o = fmax(VP - VPUN, 0.0);
double myo_o = 1.0 - kMYOC * TNF / (TNF + 150.0);
double SV_o  = SVMAX * str_o / (KSV + str_o) * (1.0 + kINO * (SVRR - 1.0)) * myo_o;
double CO_o  = SV_o * HR / 1000.0;
double MAP_o = fmax(CO_o * SVR0 * SVRR, 5.0);
double PP_o  = SV_o / CART;
double SBP_o = MAP_o + 2.0 * PP_o / 3.0;
double DBP_o = MAP_o - PP_o / 3.0;
double PV_o  = fmax(PVMIN, PV0 + (VP - VP0) / FPLAS / CVEN);
double rr_o  = RRATIO0 * (1.0 + kRRAT * fmax(SVRR - 1.0, 0.0) / SVRMAX);
double Pc_o  = (MAP_o + rr_o * PV_o) / (1.0 + rr_o);
double Cp_o  = PRP / VP * 100.0;
double Ci_o  = PRI / VI * 100.0;
double PIp_o = copP(Cp_o) + COLLCOP * (COLL / VP * 100.0);
double PIi_o = copI(Ci_o);
double Pi_o  = PIMAX - (PIMAX - PI0) * exp(-fmax(VI - VI0, 0.0) / VSC);
double Jv_o  = KF0 * LPX * ((Pc_o - Pi_o) - sig_o * (PIp_o - PIi_o));
double Jl_o  = fmax(0.0, fmin(JL0 + kLYMPH * (Pi_o - PI0), JLMAX));
double Pser_o = PSER0 + fmax(VSER - VSERB, 0.0) / CSER;
double Cser_o = TP_S0 + (TP_SMAX - TP_S0) * (SIG0 - sig_o) / fmax(SIG0 - SIGMIN, 1e-9);
double Jser_o = KFS * LPX * ((Pc_o - Pser_o) - sig_o * (PIp_o - copI(Cser_o)));
double Dser_o = fmin(DSERMAX, DSER0 * sqrt(fmax(VSER, 0.0) / VSERB));
double Hct_o  = RBCV / (RBCV + VP) * 100.0;
double Hct0_o = RBC0 / (RBC0 + VP0) * 100.0;

// Haemostasis is a PRODUCT of four independent requirements, so the worst
// term sets the answer.  This is the whole reason prophylactic platelet
// transfusion does nothing: at the nadir the worst term is the vessel wall.
double f_plt = 1.0 / (1.0 + pow(fmax(PLT, 1e-6) / 18.0, 2.2));
double f_coa = 1.0 / (1.0 + pow(fmax(FIB, 1e-6) / 90.0, 2.5));
double f_ves = (1.0 - GLX) * (1.0 - GLX);
double f_shk = 1.0 / (1.0 + pow(MAP_o / 52.0, 5.0));
double BLEED = 1.0 - (1.0 - 0.95 * f_plt) * (1.0 - 0.75 * f_coa)
             * (1.0 - 0.60 * f_ves) * (1.0 - 0.80 * f_shk);

double WARN = (Hct_o >= Hct0_o * 1.20) + (PLT < 100.0)
            + ((VSER - VSERB) > 250.0) + (AST > 200.0)
            + (PP_o <= 20.0) + (LAC > 2.5);

capture ILLDAY = self.time / 24.0;
capture EADE   = E_o;
capture NEUTF  = Nn_o;
capture SIGMA  = sig_o;
capture Jv     = Jv_o;
capture Jlymph = Jl_o;
capture Jser   = Jser_o;
capture Dser   = Dser_o;
capture NETSER = Jser_o - Dser_o;
capture Pcap   = Pc_o;
capture Pint   = Pi_o;
capture Pser   = Pser_o;
capture COPp   = PIp_o;
capture COPi   = PIi_o;
capture Hct    = Hct_o;
capture HCTRISE = (Hct_o - Hct0_o) / Hct0_o * 100.0;
capture SV     = SV_o;
capture CO     = CO_o;
capture MAP    = MAP_o;
capture SBP    = SBP_o;
capture DBP    = DBP_o;
capture PP     = PP_o;
capture TPROT  = Cp_o;
capture ALBUM  = 0.575 * Cp_o;
capture EFFUS  = fmax(VSER - VSERB, 0.0);
capture LOG10V = log10(fmax(V, 1e-9));
capture BLEEDIX = BLEED;
capture WARNSIGN = WARN;
capture SHOCK  = ((PP_o <= 20.0) || (MAP_o < 60.0)) ? 1.0 : 0.0;
capture FLUIDBAL = CIN - COUT;
capture DO2out = CO_o * (Hct_o / 3.0) * 1.34 * 0.97 * 10.0 / 1000.0;

$CAPTURE ABH0 V0 RCRYS RCOLL
'

mod <- mcode_cache("denv_qsp", code)


## ============================================================================
##  R DRIVER -- scenarios, protocols and the analyses reported in README.md
## ============================================================================
##  The WHO 2009 compensated-shock fluid ladder, as a data_set of infusion-rate
##  changes.  Rates are mL/kg/h, anchored at the moment the model ITSELF first
##  satisfies the WHO admission rule of three warning signs -- t = 48 h in the
##  untreated secondary case, which is a model output and not a chosen hour.
## ----------------------------------------------------------------------------

T_PRESENT <- 48

who_ladder <- function(t0, wt = 70, scale = 1) {
  a <- c(0, 1, 3, 5, 12, 36, 48)
  r <- c(10, 6, 4, 2.5, 1.5, 1, 0)
  data.frame(time = t0 + a, RCRYS = r * wt * scale)
}

colloid_bolus <- function(t0, wt = 70, mlkg = 15, over_h = 1, n = 2, gap = 6) {
  do.call(rbind, lapply(seq_len(n) - 1, function(k)
    data.frame(time  = t0 + gap * k + c(0, over_h),
               RCOLL = c(mlkg * wt / over_h, 0))))
}

## Compartment indices used for bolus events (order of $CMT):
##   24 = PLT, 38 = AVD (antiviral depot), 40 = APAP, 41 = STER
CMT_PLT <- 24; CMT_AVD <- 38; CMT_APAP <- 40; CMT_STER <- 41

scenarios <- list(
  S01_primary_naive           = list(ABH0 = 0),
  S02_secondary_untreated     = list(ABH0 = 55),
  S03_secondary_WHO_fluids    = list(ABH0 = 55, crys = 1.00),
  S04_under_resuscitation_50  = list(ABH0 = 55, crys = 0.50),
  S05_over_resuscitation_200  = list(ABH0 = 55, crys = 2.00),
  S06_colloid_rescue          = list(ABH0 = 55, crys = 0.55, colloid = TRUE),
  S07_albumin_rescue          = list(ABH0 = 55, crys = 0.55, alb = 15),
  S08_prophylactic_platelets  = list(ABH0 = 55, crys = 1.00, plt = 38),
  S09_steroid_at_presentation = list(ABH0 = 55, crys = 1.00, ster = T_PRESENT),
  S09b_steroid_at_fever_onset = list(ABH0 = 55, crys = 1.00, ster = 0),
  S10_antiviral_day2_illness  = list(ABH0 = 55, crys = 1.00, av = 48),
  S11_antiviral_at_present    = list(ABH0 = 55, crys = 1.00, av = T_PRESENT),
  S12_vaccine_seronegative    = list(ABH0 = 48),
  S13_vaccine_seropositive    = list(ABH0 = 3100),
  S14_tertiary_high_titre     = list(ABH0 = 2600),
  S15_paracetamol_high_dose   = list(ABH0 = 55, crys = 1.00, apap = TRUE),
  S16_anti_NS1_mab            = list(ABH0 = 55, crys = 1.00, NS1MAB = 0.35),
  S17_early_oral_rehydration  = list(ABH0 = 55, RORAL = 95),
  S18_delayed_presentation_8h = list(ABH0 = 55, crys = 1.00, delay = 8))

## Build one mrgsolve data_set.  Infusions are expressed as PARAMETER changes
## (RCRYS, RCOLL, RALB are rates read directly by $ODE), boluses as ordinary
## evid = 1 records.  Every row carries every column, because a parameter
## column left NA on a bolus row would silently blank the infusion.
build_events <- function(s) {
  tp <- T_PRESENT + if (is.null(s$delay)) 0 else s$delay

  rate_rows <- data.frame(time = 0, RCRYS = 0, RCOLL = 0, RALB = 0)
  if (!is.null(s$crys)) {
    w <- who_ladder(tp, scale = s$crys)
    rate_rows <- rbind(rate_rows,
                       data.frame(time = w$time, RCRYS = w$RCRYS, RCOLL = 0, RALB = 0))
  }
  if (isTRUE(s$colloid)) {
    cb <- colloid_bolus(tp)
    rate_rows <- rbind(rate_rows,
                       data.frame(time = cb$time, RCRYS = NA, RCOLL = cb$RCOLL, RALB = 0))
  }
  if (!is.null(s$alb)) {
    rate_rows <- rbind(rate_rows,
                       data.frame(time = tp + c(0, 2), RCRYS = NA, RCOLL = NA,
                                  RALB = c(s$alb, 0)))
  }
  rate_rows <- rate_rows[order(rate_rows$time), ]
  # carry each rate forward from the last row that actually set it
  for (nm in c("RCRYS", "RCOLL", "RALB")) {
    v <- rate_rows[[nm]]
    for (i in seq_along(v)) if (is.na(v[i])) v[i] <- if (i > 1) v[i - 1] else 0
    rate_rows[[nm]] <- v
  }
  rate_rows$cmt <- 0; rate_rows$evid <- 2; rate_rows$amt <- 0

  bol <- NULL
  addb <- function(times, cmt, amt) data.frame(time = times, cmt = cmt,
                                               evid = 1, amt = amt)
  if (!is.null(s$plt))  bol <- rbind(bol, addb(tp + c(0, 24, 48), CMT_PLT, s$plt))
  if (!is.null(s$ster)) bol <- rbind(bol, addb(s$ster + c(0, 24, 48), CMT_STER, 500 / 84))
  if (!is.null(s$av))   bol <- rbind(bol, addb(s$av + seq(0, 13) * 12, CMT_AVD, 450))
  if (isTRUE(s$apap))   bol <- rbind(bol, addb(120 + seq(0, 27) * 6, CMT_APAP, 1000 / 49))

  if (!is.null(bol)) {
    # a bolus must inherit whatever infusion rates are running at that instant
    idx <- findInterval(bol$time, rate_rows$time)
    idx[idx < 1] <- 1
    bol$RCRYS <- rate_rows$RCRYS[idx]
    bol$RCOLL <- rate_rows$RCOLL[idx]
    bol$RALB  <- rate_rows$RALB[idx]
  }

  ev <- rbind(rate_rows[, c("time", "cmt", "evid", "amt", "RCRYS", "RCOLL", "RALB")],
              if (is.null(bol)) NULL
              else bol[, c("time", "cmt", "evid", "amt", "RCRYS", "RCOLL", "RALB")])
  ev <- ev[order(ev$time, ev$evid), ]
  ev$ID <- 1
  ev
}

run_scenario <- function(name, tend = 336, delta = 0.5) {
  s  <- scenarios[[name]]
  ev <- build_events(s)
  pars <- s[intersect(names(s), c("ABH0", "RORAL", "NS1MAB"))]
  out <- mod
  if (length(pars)) out <- mrgsolve::param(out, pars)
  out <- mrgsolve::mrgsim(mrgsolve::data_set(out, ev), end = tend, delta = delta)
  d <- as.data.frame(out)
  d$scenario <- name
  d
}

run_all <- function() do.call(rbind, lapply(names(scenarios), run_scenario))

compare_table <- function(all) {
  do.call(rbind, lapply(split(all, all$scenario), function(d) data.frame(
    scenario      = d$scenario[1],
    peakV_log10   = round(max(d$LOG10V), 2),
    peak_NS1      = round(max(d$NS1)),
    defervesce_d  = d$ILLDAY[which(d$TEMP < 37.8 & d$ILLDAY > 1)[1]],
    min_sigma     = round(min(d$SIGMA), 3),
    hct_rise_pct  = round(max(d$HCTRISE), 1),
    min_PP        = round(min(d$PP), 1),
    min_MAP       = round(min(d$MAP), 1),
    effusion_mL   = round(max(d$EFFUS)),
    min_PLT       = round(min(d$PLT)),
    min_albumin   = round(min(d$ALBUM), 2),
    max_AST       = round(max(d$AST)),
    max_lactate   = round(max(d$LAC), 2),
    shock_hours   = sum(d$SHOCK) * 0.5,
    max_warnsigns = max(d$WARNSIGN),
    row.names = NULL)))
}

## ----------------------------------------------------------------------------
##  Analyses that need no ODE solver at all
## ----------------------------------------------------------------------------

## The ADE bell curve.  Peak E = 12.83 at 1:57; E returns through 1.0 at 1:696;
## E(1:1280) = 0.231.  Salje 2018 measured raised risk at 1:21-1:80 and
## protection above 1:1280.
ade_curve <- function(A, FCROSS = 0.12, NT50 = 30, HNEUT = 2.5,
                      KOPS = 4, PHI = 14) {
  aeff <- FCROSS * A
  N <- ifelse(aeff > 0, aeff^HNEUT / (aeff^HNEUT + NT50^HNEUT), 0)
  (1 - N) * (1 + (PHI - 1) * A / (A + KOPS))
}

ade_peak_titre <- function(...) {
  A <- 10^seq(-1, 4.4, length.out = 4000)
  A[which.max(ade_curve(A, ...))]
}

## The infant prediction: a cord titre, an IgG half-life and one logarithm.
## Returns 6.3 months; observed peak of infant DHF is 6-8 months, and no
## infant datum entered the calculation.
infant_age_months <- function(cord = 1280, t_half_d = 43) {
  log2(cord / ade_peak_titre()) * t_half_d / 30.44
}

## What the WHO 20 % haemoconcentration criterion actually costs.
## 41.7 % -> 50.0 %, i.e. 800 mL and 28.6 % of the plasma compartment gone.
hct_cost <- function(RBC = 2000, VP = 2800) {
  h0 <- RBC / (RBC + VP) * 100
  h1 <- 1.20 * h0
  c(hct0 = h0, threshold = h1, plasma_lost_mL = VP - RBC * (100 - h1) / h1,
    plasma_lost_pct = (VP - RBC * (100 - h1) / h1) / VP * 100)
}

## Net filtration pressure and J_v as sigma falls, everything else at baseline.
starling_table <- function(sigmas = c(0.97, 0.90, 0.80, 0.70, 0.60, 0.55)) {
  copP <- function(c) 2.1 * c + 0.16 * c^2 + 0.009 * c^3
  copI <- function(c) 3.6 * c + 0.20 * c^2
  Pc <- (88 + 7.6 * 8) / 8.6
  dPI <- copP(7.3) - copI(2.0)
  nfp <- (Pc + 3) - sigmas * dPI
  data.frame(sigma = sigmas, NFP_mmHg = round(nfp, 3),
             Jv_mL_h = round(91.77324 * nfp),
             minus_baseline_lymph = round(91.77324 * nfp - 140))
}

## ============================================================================
##  EXAMPLE
## ============================================================================
if (interactive() && !exists("DENV_NO_DEMO")) {
  all <- run_all()
  print(compare_table(all), row.names = FALSE)

  cat(sprintf("\nADE peak titre 1:%.0f, infant peak age %.1f months\n",
              ade_peak_titre(), infant_age_months()))
  print(hct_cost())
  print(starling_table())

  ## The headline figure: both channels of the same antibody rise, one axis.
  d <- run_scenario("S02_secondary_untreated")
  op <- par(mfrow = c(2, 2))
  plot(d$ILLDAY, d$TEMP,  type = "l", xlab = "illness day", ylab = "core temp (C)")
  plot(d$ILLDAY, d$SIGMA, type = "l", xlab = "illness day", ylab = "sigma")
  plot(d$ILLDAY, d$PP,    type = "l", xlab = "illness day", ylab = "pulse pressure")
  abline(h = 20, lty = 2, col = "red")
  plot(d$ILLDAY, d$EFFUS, type = "l", xlab = "illness day", ylab = "effusion (mL)")
  par(op)
}


## ============================================================================
##  HEADLINE RESULTS (all from denv_reference_model.py, LSODA, rtol 1e-6)
# ---------------------------------------------------------------------------
# 1. THE ENHANCEMENT CURVE IS A BELL BECAUSE IT IS A PRODUCT.
#    E(A) = (1-N(A)) * [1 + (PHI-1)*A/(A+KOPS)] peaks at E = 12.83 at a
#    reciprocal titre of 1:57 and crosses back through 1.0 at 1:696.
#    Salje 2018 measured raised risk at 1:21-1:80 and protection above 1:1280.
#
# 2. THE AGE OF PEAK INFANT DENGUE IS A DIVISION PROBLEM.
#    log2(1280 / 57) = 4.49 maternal-IgG half-lives; x 43 d = 193 d =
#    6.3 months.  Observed peak of infant DHF: 6-8 months.  No infant datum
#    entered the model.
#
# 3. A 20 % HAEMATOCRIT RISE MEANS 800 mL HAS ALREADY GONE (28.6 % of the
#    plasma compartment), which is why the criterion fires late and why any
#    fluid already given erases it.
#
# 4. DEFERVESCENCE AND THE CLINICAL NADIR ARE ONE EVENT.  Temperature is
#    computed from hypothalamic PGE2 driven by interferon; the leak from a
#    glycocalyx driven by NS1, TNF, VEGF and chymase.  The two chains share no
#    equation after the antibody rise, yet the narrowest pulse pressure falls
#    21.5 h before defervescence and the haematocrit peak 20.5 h before it.
#
# 5. THE FLUID DOSE-RESPONSE HAS AN INTERIOR OPTIMUM at 0.75 x the WHO ladder
#    (~5.0 L).  Zero fluid: 66 h of shock, lactate 4.35, severity 0.372.
#    3x: 4.5 h of shock but 2737 mL of effusion, severity 0.297.  Optimum
#    0.249.  The curve is U-shaped because crystalloid raises Pc and dilutes
#    Pi_p, and both terms increase J_v.
#
# 6. RESUSCITATION EFFICIENCY COLLAPSES THREEFOLD.  With sigma held at 0.97 it
#    takes 5.5 mL of crystalloid to add 1 mL to the plasma volume; at the leak
#    nadir it takes 16.2 mL, and 29.7 mL at 3x the ladder.
#
# 7. THE ANTIVIRAL WINDOW CLOSES BEFORE ANYONE PRESENTS.  Benefit is 25 % at
#    fever onset, 6.2 % at 24 h, 0.2 % at 48 h and nil thereafter, because the
#    NS1 exposure integral is 84 % spent by 24 h.  Balapiravir enrolled at
#    <=72 h and celgosivir at <=48 h; both were negative.
#
# 8. CORTICOSTEROIDS HAVE THE SAME CLIFF.  Methylprednisolone from fever onset
#    cuts severity 0.249 -> 0.195 but raises the NS1 AUC from 15.6k to 21.1k
#    by delaying viral clearance; from 72 h it is worth 5 %.  Tam 2012 enrolled
#    within 72 h and found nothing.
#
# 9. PROPHYLACTIC PLATELETS RAISE THE COUNT AND NOT THE HAEMOSTASIS.  Three
#    pools move the nadir 32 -> 44 x10^9/L and the bleeding index 0.490 ->
#    0.444, because haemostasis is a product and the worst term at the nadir
#    is the vessel wall.  AAPT (Lye 2017) found the same.
#
# 10. COLLOID'S ADVANTAGE IS PROPORTIONAL TO sigma.  Margin over crystalloid
#    grows -0.016 -> -0.027 -> -0.030 as the sigma defect deepens, because the
#    oncotic term colloid adds enters the Starling equation multiplied by
#    sigma.  Wills 2005 found no overall difference and a benefit confined to
#    the narrowest-pulse-pressure stratum.
#
# 11. LESION-BY-LESION, sigma IS THE DISEASE.  Holding sigma at 0.97 cuts
#    plasma loss 30.6 % -> 6.6 % and effusion 1825 -> 77 mL.  Holding Kf at
#    baseline changes plasma loss only 30.6 % -> 27.6 %.  Permeability is a
#    detail; sieving is the disease.
# ---------------------------------------------------------------------------
