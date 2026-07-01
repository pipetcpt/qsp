## =============================================================================
## Retinitis Pigmentosa (RP) — mrgsolve QSP model
## RHO/RPGR/USH2A/PDE6/RPE65 genotype-dependent primary rod photoreceptor
##   apoptosis (rhodopsin misfolding/ER-stress, PDE6-cGMP-Ca2+ excitotoxicity,
##   RPGR-ciliopathy transport failure) -> secondary cone death (RdCVF loss,
##   oxidative stress) -> microglial/gliotic feedback -> clinical vision loss,
##   coupled to voretigene neparvovec (AAV2-RPE65 subretinal gene therapy),
##   investigational RPGR gene augmentation (AAV8/AAV5-RPGR), MCO-010
##   optogenetic gene therapy (AAV2 multi-characteristic opsin, intravitreal),
##   CNTF encapsulated-cell neuroprotection, N-acetylcysteine (antioxidant),
##   and vitamin A palmitate/DHA supplementation PK/PD.
##
## Time unit: hours. Disease horizon simulated over years (end = 24*365*N).
##
## Calibration anchors (see rp_references.md for full PMID list):
##   - Natural history ERG/VF decline (~18%/yr exponential): Berson 1985 Arch
##     Ophthalmol; Sandberg 2005 IOVS; Grover 1997 Ophthalmology genotype-
##     specific progression rates (RPE65/LCA fastest, USH2A slower)
##   - Voretigene neparvovec phase 3 (Luxturna): Russell 2017 Lancet (MLMT,
##     FST, multi-luminance improvement sustained to 1yr); Maguire 2019/2021
##     Ophthalmology (4-yr durability, stable transgene expression)
##   - RPGR gene therapy (XLRP): Cehajic-Kapetanovic 2020 Nat Med (AAV8-RPGR
##     phase 1/2, dose-dependent retinal thinning at high dose); Pennesi 2022
##     XIRIUS trial (AAV5-RPGR, cotoretigene toliparvovec)
##   - MCO-010 optogenetic therapy: Nanoscope RESTORE trial topline (2023),
##     Sahel 2021 Nat Med (ChrimsonR partial vision recovery case report),
##     Busskamp 2010 Science (ChR2 proof-of-concept)
##   - CNTF encapsulated cell implant: Birch 2013 IOVS (paradoxical ERG
##     suppression despite structural rescue signal), Sieving 2006 PNAS
##   - Vitamin A palmitate + DHA: Berson 1993 Arch Ophthalmol (DBA trial,
##     ~20% slower ERG decline), Berson 2004 Arch Ophthalmol (DHA substudy)
##   - N-acetylcysteine: Campochiaro 2020 IOVS phase 1 dose escalation
##   - Cystoid macular edema / CAI: Fishman 1989/Chen 2016 (dorzolamide,
##     acetazolamide reduce CST in RP-CME, ~50% partial responders)
##   - Retinal prosthesis: Humayun 2012 Ophthalmology (Argus II), da Cruz
##     2016 Ophthalmology (5-yr outcomes)
## =============================================================================

$PROB
# Retinitis Pigmentosa QSP model (23-compartment PK/PD + disease + safety)

$PARAM
// ---- Genotype / severity ----
SEVERITY     = 1.0     // 0.75=USH2A, 1.0=RHO-adRP, 1.6=RPGR-XLRP, 2.5=RPE65-LCA/EOSRD
IS_RPE65     = 0       // 1 = biallelic RPE65 genotype (eligible for voretigene neparvovec)
IS_XLRP      = 0       // 1 = RPGR-XLRP genotype (eligible for RPGR gene therapy)

// ---- Voretigene neparvovec (AAV2-RPE65) subretinal gene therapy ----
GT65_KTRANS   = 0.010   // 1/h, subretinal vector-to-RPE transduction rate
GT65_YIELD    = 0.85    // fraction of transduced vector producing functional protein
GT65_KLOSS    = 6.0e-7  // 1/h, transgene expression loss (~19-yr half-life, durable per Maguire)
VC_BASE_DEFICIT = 0.05  // residual visual-cycle flux in untreated biallelic RPE65-null (near-abolished)
VC_EC50       = 0.30    // RPE65 transgene expression level for half-max visual-cycle restoration
APO_OPSIN_GAIN = 1.8    // extra rod-death drive from constitutive opsin signaling when visual cycle flux is low

// ---- Investigational RPGR gene therapy (AAV8/AAV5-RPGR) ----
GTRPGR_KTRANS = 0.008
GTRPGR_YIELD  = 0.80
GTRPGR_KLOSS  = 8.0e-7  // slightly less durable than RPE65 GT (earlier-generation vector data)
RPGR_EC50     = 0.30
RPGR_EMAX     = 0.65    // max fractional correction of RPGR-driven excess rod death

// ---- MCO-010 optogenetic gene therapy (AAV2, intravitreal) ----
OP_KTRANS     = 0.006
OP_YIELD      = 0.55    // lower transduction efficiency, intravitreal route to inner retina
OP_KLOSS      = 5.0e-7
OP_EC50       = 0.35
OP_BYPASS_EMAX = 0.55   // max fractional restoration of light-sensitivity signal (bypass pathway)

// ---- CNTF encapsulated cell technology implant (NT-501, intravitreal) ----
CNTF_PROD    = 0.050    // relative units/h zero-order release while device implanted
CNTF_KOUT    = 0.010    // 1/h clearance
CNTF_EC50    = 1.2
CNTF_EMAX_ROD = 0.25     // max fractional slowing of rod apoptosis
CNTF_EMAX_CONE = 0.20    // max fractional slowing of cone apoptosis
CNTF_ERG_SUPPRESS_EMAX = 0.30  // paradoxical reversible ERG amplitude suppression (Birch 2013)
CNTF_ERG_SUPPRESS_EC50 = 1.0

// ---- N-acetylcysteine (NAC) oral, antioxidant ----
NAC_KA   = 0.35     // 1/h
NAC_F    = 0.10     // low oral bioavailability
NAC_V1   = 40        // L
NAC_CL   = 12        // L/h
NAC_EC50 = 3.0        // mg/L, half-max ROS-scavenging effect
NAC_EMAX = 0.45       // max fractional ROS clearance enhancement
NAC_CONE_PROTECT_EMAX = 0.25  // max fractional slowing of secondary cone death

// ---- Vitamin A palmitate + DHA oral supplementation ----
VITA_KIN   = 0.0020   // 1/h, buildup of chromophore-support storage pool while dosed
VITA_KOUT  = 0.0015   // 1/h, turnover
VITA_EC50  = 1.0
VITA_EMAX_ERG = 0.20   // max ~20% slowing of ERG amplitude decline (Berson 1993)

// ---- Recombinant RdCVF protein (investigational, simplified static modifier) ----
RDCVF_THERAPY_ON   = 0
RDCVF_PROTECT_FRAC = 0.35  // fractional reduction of rod-loss-driven cone death signal

// ---- Rod photoreceptor apoptosis ----
K_ROD0        = 2.1e-5  // 1/h baseline rate (~ -17%/yr amplitude decline, Berson/Sandberg natural history)
MG_FEEDBACK_GAIN = 0.9   // amplification of rod death by microglial activation

// ---- Cone secondary death ----
K_CONE0        = 1.0e-5
RDCVF_GAIN     = 1.3     // cone-death drive per unit rod loss (RdCVF depletion)
ROS_CONE_GAIN  = 0.6     // cone-death drive per unit oxidative stress

// ---- Oxidative stress (ROS) ----
ROS_KIN   = 0.012
ROS_KOUT  = 0.010

// ---- Microglial activation ----
MG_KIN    = 6.0
MG_KOUT   = 0.020

// ---- Retinal ganglion cell (relatively preserved, late decline) ----
K_RGC0    = 4.0e-6

// ---- Clinical endpoint relaxation (assessment/tissue-response lag) ----
K_CLIN    = 0.0035    // 1/h (~12-day time constant)
ERG_ROD_MAX  = 250     // uV, normal rod b-wave
ERG_ROD_HILL = 1.4
ERG_CONE_MAX = 80      // uV, normal cone amplitude
ERG_CONE_HILL = 1.2
VF_MAX       = 1800    // deg^2, normal Goldmann III4e field area
BCVA_BEST    = 0.0     // logMAR (20/20)
BCVA_WORST   = 1.6     // logMAR (~20/800)
FST_MIN      = -5      // dB, normal (very sensitive)
FST_MAX      = 55      // dB, dark/unresponsive
MLMT_MIN     = -1      // MLMT score, worst navigable light level (per Russell 2017 scale)
MLMT_MAX     = 4       // best (navigates at lowest luminance)
CME_BASELINE = 250     // um, normal central subfield thickness
CME_GAIN     = 25       // um increase per unit microglial/inflammatory index
CAI_ON       = 0        // topical dorzolamide / oral acetazolamide flag
CAI_EFFECT   = 0.40     // fractional reduction of edema-driven thickening

$CMT @annotated
GT65_VG      : Voretigene neparvovec subretinal vector genome (relative units)
GT65_EXPR    : RPE65 transgene protein expression (relative units)
GTRPGR_VG    : RPGR gene-therapy vector genome (relative units)
GTRPGR_EXPR  : RPGR transgene expression (relative units)
OP_VG        : MCO-010 optogenetic vector genome (relative units)
OP_EXPR      : Multi-characteristic opsin expression (relative units)
CNTF_DEV     : CNTF encapsulated-implant tissue level (relative units)
NAC_GUT      : N-acetylcysteine gut depot (mg)
NAC_CENT     : N-acetylcysteine plasma concentration compartment (mg)
VITA_STORE   : Vitamin A/DHA chromophore-support storage pool (relative units)
ROD_FRAC     : Viable rod photoreceptor fraction (0-1)
CONE_FRAC    : Viable cone photoreceptor fraction (0-1)
ROS          : Oxidative stress index (relative units)
MICROGLIA    : Microglial activation index (relative units)
RGC_FRAC     : Viable retinal ganglion cell fraction (0-1)
ERG_ROD      : Full-field ERG rod b-wave amplitude (uV)
ERG_CONE     : Full-field ERG cone amplitude (uV)
VF_AREA      : Goldmann kinetic visual field area (deg^2)
BCVA         : Best-corrected visual acuity (logMAR)
FST          : Full-field stimulus threshold (dB)
MLMT         : Multi-luminance mobility test score
CME_CST      : Central subfield thickness, cystoid macular edema (um)
FX_YEARS     : Elapsed simulation time tracker (years, for plotting convenience)

$MAIN
double GT65_Cp   = GT65_EXPR;
double GTRPGR_Cp = GTRPGR_EXPR;
double OP_Cp     = OP_EXPR;
double NAC_Conc  = NAC_CENT / NAC_V1;

// ---- Visual cycle restoration (RPE65 genotype only) ----
double visual_cycle = 1.0;
if (IS_RPE65 == 1) {
  visual_cycle = VC_BASE_DEFICIT + (1.0-VC_BASE_DEFICIT) * (GT65_EXPR/(VC_EC50+GT65_EXPR+1e-9));
}
double apo_opsin_toxicity = (IS_RPE65==1) ? APO_OPSIN_GAIN*(1.0-visual_cycle) : 0.0;

// ---- RPGR ciliary correction (XLRP genotype only) ----
double cilia_correction = (IS_XLRP==1) ? RPGR_EMAX*(GTRPGR_EXPR/(RPGR_EC50+GTRPGR_EXPR+1e-9)) : 0.0;

// ---- CNTF neuroprotection ----
double cntf_frac = CNTF_DEV/(CNTF_EC50+CNTF_DEV+1e-9);
double cntf_protect_rod  = CNTF_EMAX_ROD*cntf_frac;
double cntf_protect_cone = CNTF_EMAX_CONE*cntf_frac;
double cntf_erg_suppress = CNTF_ERG_SUPPRESS_EMAX*(CNTF_DEV/(CNTF_ERG_SUPPRESS_EC50+CNTF_DEV+1e-9));

// ---- NAC antioxidant effect ----
double nac_ros_effect  = NAC_EMAX*(NAC_Conc/(NAC_EC50+NAC_Conc+1e-9));
double nac_cone_protect = NAC_CONE_PROTECT_EMAX*(NAC_Conc/(NAC_EC50+NAC_Conc+1e-9));

// ---- Vitamin A / DHA chromophore support ----
double vita_erg_support = VITA_EMAX_ERG*(VITA_STORE/(VITA_EC50+VITA_STORE+1e-9));

// ---- Optogenetic bypass signal (requires surviving RGC) ----
double op_bypass = OP_BYPASS_EMAX*(OP_EXPR/(OP_EC50+OP_EXPR+1e-9))*RGC_FRAC;

$ODE
// ---------------- Voretigene neparvovec (RPE65) subretinal gene therapy ----------------
dxdt_GT65_VG   = -GT65_KTRANS*GT65_VG;
dxdt_GT65_EXPR =  GT65_KTRANS*GT65_VG*GT65_YIELD - GT65_KLOSS*GT65_EXPR;

// ---------------- Investigational RPGR gene therapy ----------------
dxdt_GTRPGR_VG   = -GTRPGR_KTRANS*GTRPGR_VG;
dxdt_GTRPGR_EXPR =  GTRPGR_KTRANS*GTRPGR_VG*GTRPGR_YIELD - GTRPGR_KLOSS*GTRPGR_EXPR;

// ---------------- MCO-010 optogenetic gene therapy ----------------
dxdt_OP_VG   = -OP_KTRANS*OP_VG;
dxdt_OP_EXPR =  OP_KTRANS*OP_VG*OP_YIELD - OP_KLOSS*OP_EXPR;

// ---------------- CNTF encapsulated implant ----------------
dxdt_CNTF_DEV = CNTF_PROD - CNTF_KOUT*CNTF_DEV;

// ---------------- N-acetylcysteine PK ----------------
dxdt_NAC_GUT  = -NAC_KA*NAC_GUT;
dxdt_NAC_CENT =  NAC_KA*NAC_GUT*NAC_F - (NAC_CL/NAC_V1)*NAC_CENT;

// ---------------- Vitamin A / DHA storage pool ----------------
dxdt_VITA_STORE = VITA_KIN*1.0 - VITA_KOUT*VITA_STORE;   // dosing controlled via VITA_KIN reset in scenario param set

// ---------------- Oxidative stress ----------------
dxdt_ROS = ROS_KIN*(1.0-ROD_FRAC) - ROS_KOUT*(1.0+nac_ros_effect)*ROS;

// ---------------- Microglial activation (driven by instantaneous photoreceptor loss flux) ----------------
double k_rod_raw = K_ROD0*SEVERITY*(1.0+apo_opsin_toxicity)*(1.0+MG_FEEDBACK_GAIN*MICROGLIA)
                    *(1.0-cilia_correction)*(1.0-cntf_protect_rod)*(1.0-vita_erg_support*0.5);
if (k_rod_raw < 0) k_rod_raw = 0;
double k_cone_raw = K_CONE0*(RDCVF_GAIN*(1.0-ROD_FRAC)*(1.0-RDCVF_THERAPY_ON*RDCVF_PROTECT_FRAC)
                    + ROS_CONE_GAIN*ROS)*(1.0-cntf_protect_cone)*(1.0-nac_cone_protect);
if (k_cone_raw < 0) k_cone_raw = 0;
double death_flux = k_rod_raw*ROD_FRAC + k_cone_raw*CONE_FRAC;
dxdt_MICROGLIA = MG_KIN*death_flux - MG_KOUT*MICROGLIA;

// ---------------- Rod / cone photoreceptor survival ----------------
dxdt_ROD_FRAC  = -k_rod_raw*ROD_FRAC;
dxdt_CONE_FRAC = -k_cone_raw*CONE_FRAC;

// ---------------- Retinal ganglion cell (late, slow decline) ----------------
dxdt_RGC_FRAC = -K_RGC0*(1.0-ROD_FRAC)*(1.0-CONE_FRAC)*RGC_FRAC;

// ---------------- Clinical endpoints (relaxation toward mechanistic targets) ----------------
double erg_rod_target  = ERG_ROD_MAX*pow(ROD_FRAC, ERG_ROD_HILL)*(1.0-cntf_erg_suppress)*(1.0+vita_erg_support);
double erg_cone_target = ERG_CONE_MAX*pow(CONE_FRAC, ERG_CONE_HILL)*(1.0-cntf_erg_suppress*0.5);
double photoreceptor_signal = ROD_FRAC + 0.3*CONE_FRAC;
if (photoreceptor_signal > 1.0) photoreceptor_signal = 1.0;
double vf_target   = VF_MAX*pow(ROD_FRAC,0.7)*pow(CONE_FRAC,0.3);
double bcva_target = BCVA_BEST + (BCVA_WORST-BCVA_BEST)*(1.0-CONE_FRAC);
double fst_target  = FST_MAX - (FST_MAX-FST_MIN)*photoreceptor_signal - op_bypass*(FST_MAX-FST_MIN)*0.5;
if (fst_target < FST_MIN) fst_target = FST_MIN;
double mlmt_target = MLMT_MIN + (MLMT_MAX-MLMT_MIN)*(photoreceptor_signal*0.7 + op_bypass*0.3);
if (mlmt_target > MLMT_MAX) mlmt_target = MLMT_MAX;
double cme_target  = CME_BASELINE + CME_GAIN*MICROGLIA*(1.0-CAI_ON*CAI_EFFECT);

dxdt_ERG_ROD  = K_CLIN*(erg_rod_target-ERG_ROD);
dxdt_ERG_CONE = K_CLIN*(erg_cone_target-ERG_CONE);
dxdt_VF_AREA  = K_CLIN*(vf_target-VF_AREA);
dxdt_BCVA     = K_CLIN*(bcva_target-BCVA);
dxdt_FST      = K_CLIN*(fst_target-FST);
dxdt_MLMT     = K_CLIN*(mlmt_target-MLMT);
dxdt_CME_CST  = K_CLIN*(cme_target-CME_CST);
dxdt_FX_YEARS = 1.0/(24.0*365.0);

$INIT
GT65_VG = 0, GT65_EXPR = 0, GTRPGR_VG = 0, GTRPGR_EXPR = 0,
OP_VG = 0, OP_EXPR = 0, CNTF_DEV = 0, NAC_GUT = 0, NAC_CENT = 0, VITA_STORE = 0,
ROD_FRAC = 1.0, CONE_FRAC = 1.0, ROS = 0, MICROGLIA = 0, RGC_FRAC = 1.0,
ERG_ROD = 250, ERG_CONE = 80, VF_AREA = 1800, BCVA = 0, FST = -5, MLMT = 4,
CME_CST = 250, FX_YEARS = 0

$CAPTURE visual_cycle apo_opsin_toxicity cilia_correction cntf_protect_rod cntf_protect_cone
$CAPTURE nac_ros_effect vita_erg_support op_bypass k_rod_raw k_cone_raw photoreceptor_signal

## =============================================================================
## Treatment scenarios (see rp_shiny_app.R for interactive dosing UI)
##
## 1. Natural history, RHO-adRP (SEVERITY=1.0, untreated)
## 2. Natural history, RPGR-XLRP (SEVERITY=1.6, faster progression)
## 3. Natural history, RPE65-LCA/EOSRD (SEVERITY=2.5, IS_RPE65=1, untreated)
## 4. Voretigene neparvovec (Luxturna) subretinal gene therapy, RPE65-biallelic
##    (IS_RPE65=1, one-time bolus into GT65_VG; Russell 2017 Lancet regimen)
## 5. Investigational RPGR gene therapy (AAV8/AAV5-RPGR), XLRP (IS_XLRP=1,
##    one-time bolus into GTRPGR_VG; Cehajic-Kapetanovic 2020 / XIRIUS dosing)
## 6. MCO-010 optogenetic gene therapy, genotype-agnostic late-stage disease
##    (intravitreal bolus into OP_VG; requires RGC_FRAC > ~0.5 for meaningful
##    bypass signal; Nanoscope RESTORE trial regimen)
## 7. CNTF encapsulated-cell implant (continuous, CNTF_PROD "on" from t=0)
## 8. N-acetylcysteine oral (chronic BID/TID dosing into NAC_GUT)
## 9. Vitamin A palmitate 15,000 IU/day + DHA (VITA_KIN "on" from t=0;
##    Berson 1993/2004 regimen)
## 10. Combination: voretigene neparvovec + adjunct NAC neuroprotection
##
## Example mrgsolve event code (voretigene neparvovec, RPE65 genotype):
##   mod <- mread("rp_mrgsolve_model") %>% param(SEVERITY = 2.5, IS_RPE65 = 1)
##   e_gt <- ev(amt = 100, cmt = "GT65_VG", time = 0)  # one-time subretinal dose
##   out <- mod %>% ev(e_gt) %>% mrgsim(end = 24*365*10, delta = 24)  # 10-yr horizon
## =============================================================================
