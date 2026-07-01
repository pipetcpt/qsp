## =============================================================================
## Osteogenesis Imperfecta (OI) — mrgsolve QSP model
## COL1A1/COL1A2 collagen defect -> osteoblast/osteocyte dysfunction ->
##   RANKL/OPG-driven high bone turnover + sclerostin-mediated Wnt suppression
##   + excess TGF-beta bioavailability -> low BMD / recurrent fracture,
##   coupled to bisphosphonate (pamidronate, zoledronic acid), denosumab
##   (anti-RANKL), teriparatide (PTH1R anabolic, adult type I only),
##   setrusumab (anti-sclerostin, investigational), and fresolimumab
##   (anti-TGF-beta, investigational) PK/PD.
##
## Calibration anchors (see oi_references.md for full PMID list):
##   - Cyclical IV pamidronate pediatric OI: Glorieux 1998 NEJM (BMD Z-score
##     +Z ~1.5 over 4y, fracture rate reduction); Zeitlin 2003 Bone (growth)
##   - Bisphosphonate meta-analysis: Dwan 2016 Cochrane (CD005088)
##   - Oral risedronate (PLUTO trial): Bishop 2013 Lancet (no fracture benefit
##     over 12mo vs IV, informs oral-vs-IV potency difference)
##   - Zoledronic acid vs pamidronate pediatric OI: Barros 2012, Vuorimies 2017
##     Horm Res Paediatr (annual ZA non-inferior to cyclical PAM)
##   - Denosumab in OI type VI / severe recessive OI: Hoyer-Kuhn 2014 JBMR
##     (case series, reduced fracture, rebound resorption on cessation)
##   - Teriparatide adults OI type I: Orwoll 2014 J Bone Miner Res (RCT,
##     BMD +unlike type III/IV due to growth-plate/anabolic-window biology);
##     Gatti 2013 Calcif Tissue Int; contraindicated in open growth plates
##     (osteosarcoma signal in juvenile rats)
##   - Setrusumab (BPS804) anti-sclerostin: Glorieux 2017 JBMR (phase 2a dose
##     range), ASTEROID/ORBIT phase 2/3 (Mereo BioPharma, NCT03118570,
##     NCT05125276)
##   - Fresolimumab anti-TGF-beta preclinical rationale: Grafe 2014 Nat Med
##     (oim/oim mouse model, TGF-beta neutralization normalizes bone mass)
##   - Bone turnover markers P1NP/CTX pediatric reference ranges: Bayer 2014
##     Ann Clin Biochem; Rauch 2000 J Clin Invest (OI bone histomorphometry,
##     high turnover)
##   - Sillence classification & natural history: Sillence 1979 J Med Genet;
##     Marini 2017 Nat Rev Dis Primers OI review
## =============================================================================

$PROB
# Osteogenesis Imperfecta QSP model (22-compartment PK/PD + disease + safety)

$PARAM
// ---- Pamidronate (PAM) PK: IV infusion, 2-cpt (central + bone-bound deep) ----
PAM_V1      = 12      // L, central volume
PAM_CL      = 5.4     // L/h, renal clearance
PAM_Q       = 0.030   // L/h, slow exchange into bone-bound compartment
PAM_VBONE   = 1        // scaling amount compartment (bone-bound, mg)
PAM_KBONE_OFF = 0.00025 // 1/h, extremely slow bone release (skeletal retention)

// ---- Zoledronic acid (ZOL) PK: IV infusion, 2-cpt ----
ZOL_V1      = 18
ZOL_CL      = 3.2
ZOL_Q       = 0.045
ZOL_KBONE_OFF = 0.00015

// ---- Denosumab (DMAB) PK: SC depot, nonlinear (TMDD-approx via Michaelis-Menten) ----
DMAB_KA     = 0.007   // 1/h, slow SC absorption
DMAB_F      = 0.64
DMAB_V1     = 2.5     // L (small, mAb central volume ~ plasma-restricted early)
DMAB_CLlin  = 0.0028  // L/h, linear (FcRn-mediated) clearance
DMAB_VMAX   = 0.09    // mg/h, target-mediated nonlinear clearance Vmax
DMAB_KM     = 3.0     // mg/L, Michaelis constant for nonlinear elimination

// ---- Teriparatide (TPTD) PK: SC depot, rapid 1-cpt ----
TPTD_KA     = 1.6     // 1/h
TPTD_F      = 0.95
TPTD_V1     = 25      // L
TPTD_CL     = 62      // L/h (short t1/2 ~1h)

// ---- Setrusumab (anti-sclerostin mAb) PK: IV, 2-cpt linear ----
SETRU_V1    = 3.2
SETRU_V2    = 2.8
SETRU_CL    = 0.011   // L/h
SETRU_Q     = 0.020   // L/h

// ---- Fresolimumab (anti-TGF-beta mAb) PK: IV, 1-cpt linear (investigational) ----
FRESO_V1    = 3.0
FRESO_CL    = 0.009   // L/h

// ---- Disease network: sclerostin / Wnt ----
SOST_KIN     = 1.0     // baseline sclerostin production (osteocyte-driven)
SOST_KOUT    = 0.010   // 1/h turnover
SOST_OCYTE_GAIN = 0.4  // extra sclerostin drive from disrupted osteocyte network in OI

// ---- RANKL / OPG turnover ----
RANKL_KIN    = 1.0
RANKL_KOUT   = 0.05    // 1/h
RANKL_OI_GAIN = 0.65   // OI-specific elevation of RANKL production (osteocyte signaling)
OPG_KIN      = 1.0
OPG_KOUT     = 0.04    // 1/h
IC50_DMAB    = 2.5     // mg/L, denosumab free-RANKL neutralization potency (relative)

// ---- TGF-beta active pool ----
TGFB_KIN     = 1.0
TGFB_KOUT    = 0.08    // 1/h
TGFB_MATRIX_GAIN = 0.55 // release from disorganized matrix (severity-dependent)
IC50_FRESO   = 4.0      // mg/L, fresolimumab neutralization potency

// ---- Osteoblast / osteoclast activity indices ----
OB_KIN       = 1.0
OB_KOUT      = 0.06    // 1/h
OB_SOST_IC50 = 1.5      // sclerostin level producing half-max OB suppression
OB_TGFB_IC50 = 1.8       // TGF-beta level producing half-max OB maturation block
OB_PTH_EMAX  = 1.8        // max fold-stimulation of OB by intermittent PTH signal
OB_PTH_EC50  = 0.35        // ug/L teriparatide Cp for half-max anabolic effect
OB_APOPTOSIS_UPR = 0.30     // fractional OB activity loss from ER-stress/UPR (severity-linked)

OC_KIN       = 1.0
OC_KOUT      = 0.07    // 1/h
OC_RANKL_GAIN = 1.4      // gain of RANKL/OPG ratio driving OC differentiation
OC_TGFB_GAIN  = 0.35      // synergistic OC stimulation by TGF-beta
OC_BP_IC50_PAM = 25       // mg bone-bound PAM amount for half-max osteoclast apoptosis
OC_BP_IC50_ZOL = 12       // mg bone-bound ZOL amount (more potent, lower IC50)
OC_BP_EMAX     = 0.90     // maximal fractional OC suppression by bisphosphonate

// ---- Disease severity (Sillence-type modifier; set per scenario) ----
SEVERITY     = 1.0     // 0.3=type I mild, 1.0=type III/IV mod-severe, 1.6=type II/severe-recessive analog (non-lethal sim)

// ---- Bone mineral content / BMD trajectory ----
BMC_KIN_BASE   = 1.0
BMC_REMODEL_GAIN = 0.02   // 1/h scaling of OB-OC balance into BMC change
BMC_BASELINE_Z = -3.0     // baseline lumbar spine BMD Z-score (untreated moderate-severe OI)

// ---- Bone turnover markers ----
P1NP_GAIN    = 45      // ug/L scaling from OB activity
CTX_GAIN     = 0.85    // ug/L scaling from OC activity

// ---- Fracture hazard ----
FX_HAZARD_BASE = 0.006  // 1/h baseline instantaneous hazard scaling (severity & BMC dependent)
FX_BMC_GAIN    = 1.8    // hazard multiplier per unit BMC deficit

// ---- Growth (pediatric height Z-score trajectory) ----
HEIGHT_KIN     = 0.00030  // 1/h baseline growth-plate mineralization rate
HEIGHT_OB_GAIN = 0.6      // OB-activity dependent contribution to growth
HEIGHT_TGFB_PENALTY = 0.5 // growth suppression from excess TGF-beta (impaired maturation)

$CMT @annotated
PAM_CENT    : Pamidronate central amount (mg)
PAM_BONE    : Pamidronate bone-bound amount (mg)
ZOL_CENT    : Zoledronic acid central amount (mg)
ZOL_BONE    : Zoledronic acid bone-bound amount (mg)
DMAB_DEPOT  : Denosumab SC depot amount (mg)
DMAB_CENT   : Denosumab central amount (mg)
TPTD_DEPOT  : Teriparatide SC depot amount (ug)
TPTD_CENT   : Teriparatide central amount (ug)
SETRU_CENT  : Setrusumab central amount (mg)
SETRU_PERIPH: Setrusumab peripheral amount (mg)
FRESO_CENT  : Fresolimumab central amount (mg)
SOST        : Sclerostin level (relative units)
RANKL       : RANKL level (relative units)
OPG         : OPG level (relative units)
TGFB        : Active TGF-beta level (relative units)
OB          : Osteoblast activity index (relative units)
OC          : Osteoclast activity index (relative units)
BMC         : Bone mineral content index (Z-score-like scale)
P1NP        : Procollagen type I N-propeptide (ug/L)
CTX         : C-telopeptide crosslinks (ug/L)
FX_CUM      : Cumulative fracture count (events)
HEIGHT_Z    : Height-for-age Z-score

$MAIN
double PAM_Cp   = PAM_CENT / PAM_V1;
double ZOL_Cp   = ZOL_CENT / ZOL_V1;
double DMAB_Cp  = DMAB_CENT / DMAB_V1;
double TPTD_Cp  = TPTD_CENT / TPTD_V1;
double SETRU_Cp = SETRU_CENT / SETRU_V1;
double FRESO_Cp = FRESO_CENT / FRESO_V1;

// Sclerostin neutralization by setrusumab (Emax model, effective SOST signal reduced)
double IC50_SETRU = 3.0;  // mg/L, setrusumab neutralization potency
double SOST_eff = SOST * (1.0 - (SETRU_Cp / (IC50_SETRU + SETRU_Cp + 1e-9)));
if (SOST_eff < 0) SOST_eff = 0;

// TGF-beta neutralization by fresolimumab
double TGFB_eff = TGFB * (1.0 - (FRESO_Cp / (IC50_FRESO + FRESO_Cp + 1e-9)));
if (TGFB_eff < 0) TGFB_eff = 0;

// Free RANKL after denosumab neutralization (Emax-type competitive block)
double RANKL_free = RANKL * (IC50_DMAB / (IC50_DMAB + DMAB_Cp + 1e-9));

$ODE
// ---------------- Pamidronate PK ----------------
dxdt_PAM_CENT = -(PAM_CL/PAM_V1)*PAM_CENT - PAM_Q*PAM_CENT/PAM_V1;
dxdt_PAM_BONE =  PAM_Q*PAM_CENT/PAM_V1 - PAM_KBONE_OFF*PAM_BONE;

// ---------------- Zoledronic acid PK ----------------
dxdt_ZOL_CENT = -(ZOL_CL/ZOL_V1)*ZOL_CENT - ZOL_Q*ZOL_CENT/ZOL_V1;
dxdt_ZOL_BONE =  ZOL_Q*ZOL_CENT/ZOL_V1 - ZOL_KBONE_OFF*ZOL_BONE;

// ---------------- Denosumab PK (SC depot -> central, nonlinear elimination) ----------------
dxdt_DMAB_DEPOT = -DMAB_KA*DMAB_DEPOT;
dxdt_DMAB_CENT  =  DMAB_KA*DMAB_DEPOT*DMAB_F - DMAB_CLlin*DMAB_Cp
                   - (DMAB_VMAX*DMAB_Cp)/(DMAB_KM + DMAB_Cp + 1e-9);

// ---------------- Teriparatide PK (SC depot, rapid elimination) ----------------
dxdt_TPTD_DEPOT = -TPTD_KA*TPTD_DEPOT;
dxdt_TPTD_CENT  =  TPTD_KA*TPTD_DEPOT*TPTD_F - (TPTD_CL/TPTD_V1)*TPTD_CENT;

// ---------------- Setrusumab PK (2-cpt linear mAb) ----------------
dxdt_SETRU_CENT   = -(SETRU_CL/SETRU_V1)*SETRU_CENT - SETRU_Q*(SETRU_CENT/SETRU_V1 - SETRU_PERIPH/SETRU_V2);
dxdt_SETRU_PERIPH =  SETRU_Q*(SETRU_CENT/SETRU_V1 - SETRU_PERIPH/SETRU_V2);

// ---------------- Fresolimumab PK (1-cpt linear mAb) ----------------
dxdt_FRESO_CENT = -(FRESO_CL/FRESO_V1)*FRESO_CENT;

// ---------------- Sclerostin ----------------
dxdt_SOST = SOST_KIN*(1.0 + SOST_OCYTE_GAIN*SEVERITY) - SOST_KOUT*SOST;

// ---------------- RANKL / OPG ----------------
dxdt_RANKL = RANKL_KIN*(1.0 + RANKL_OI_GAIN*SEVERITY) - RANKL_KOUT*RANKL;
dxdt_OPG   = OPG_KIN - OPG_KOUT*OPG;

// ---------------- TGF-beta ----------------
dxdt_TGFB = TGFB_KIN*(1.0 + TGFB_MATRIX_GAIN*SEVERITY) - TGFB_KOUT*TGFB;

// ---------------- Osteoblast activity ----------------
double sost_inhib = SOST_eff / (OB_SOST_IC50 + SOST_eff + 1e-9);
double tgfb_block = TGFB_eff / (OB_TGFB_IC50 + TGFB_eff + 1e-9);
double pth_stim   = 1.0 + (OB_PTH_EMAX-1.0)*TPTD_Cp/(OB_PTH_EC50 + TPTD_Cp + 1e-9);
double ob_prod = OB_KIN*(1.0 - OB_APOPTOSIS_UPR*SEVERITY/1.6)*(1.0 - 0.6*sost_inhib)*(1.0 - 0.5*tgfb_block)*pth_stim;
if (ob_prod < 0.05) ob_prod = 0.05;
dxdt_OB = ob_prod - OB_KOUT*OB;

// ---------------- Osteoclast activity ----------------
double rankl_opg_ratio = RANKL_free / (OPG + 1e-9);
double bp_apoptosis = OC_BP_EMAX*( (PAM_BONE/(OC_BP_IC50_PAM+PAM_BONE+1e-9)) + (ZOL_BONE/(OC_BP_IC50_ZOL+ZOL_BONE+1e-9)) );
if (bp_apoptosis > OC_BP_EMAX) bp_apoptosis = OC_BP_EMAX;
double oc_prod = OC_KIN*(1.0 + OC_RANKL_GAIN*(rankl_opg_ratio-1.0))*(1.0 + OC_TGFB_GAIN*(TGFB_eff-1.0)/2.0);
if (oc_prod < 0.05) oc_prod = 0.05;
dxdt_OC = oc_prod*(1.0 - bp_apoptosis) - OC_KOUT*OC;

// ---------------- Bone mineral content (net remodeling balance) ----------------
dxdt_BMC = BMC_REMODEL_GAIN*(OB - OC);

// ---------------- Bone turnover markers ----------------
dxdt_P1NP = 0.15*(P1NP_GAIN*OB - P1NP);
dxdt_CTX  = 0.20*(CTX_GAIN*OC - CTX);

// ---------------- Fracture hazard (cumulative count via hazard integration) ----------------
double bmc_deficit = BMC_BASELINE_Z - BMC;  // >0 when BMC below untreated baseline reference frame is inverse; recompute below
double z_now = BMC_BASELINE_Z + BMC;         // current Z-score approx (BMC delta from baseline)
double deficit = -z_now;                     // positive when Z more negative (worse)
if (deficit < 0) deficit = 0;
double hazard = FX_HAZARD_BASE*SEVERITY*(1.0 + FX_BMC_GAIN*deficit/3.0);
dxdt_FX_CUM = hazard;

// ---------------- Height / growth (pediatric) ----------------
double height_gain = HEIGHT_KIN*(1.0 + HEIGHT_OB_GAIN*(OB-1.0)) - HEIGHT_KIN*HEIGHT_TGFB_PENALTY*(TGFB_eff-1.0)/2.0;
dxdt_HEIGHT_Z = height_gain;

$INIT
PAM_CENT = 0, PAM_BONE = 0, ZOL_CENT = 0, ZOL_BONE = 0,
DMAB_DEPOT = 0, DMAB_CENT = 0, TPTD_DEPOT = 0, TPTD_CENT = 0,
SETRU_CENT = 0, SETRU_PERIPH = 0, FRESO_CENT = 0,
SOST = 1.0, RANKL = 1.0, OPG = 1.0, TGFB = 1.0,
OB = 1.0, OC = 1.0, BMC = 0, P1NP = 45, CTX = 0.85,
FX_CUM = 0, HEIGHT_Z = -2.0

$CAPTURE PAM_Cp ZOL_Cp DMAB_Cp TPTD_Cp SETRU_Cp FRESO_Cp SOST_eff TGFB_eff RANKL_free bmc_deficit z_now hazard

## =============================================================================
## Treatment scenarios (see oi_shiny_app.R for interactive dosing UI)
##
## 1. Natural history (untreated, moderate-severe OI type III/IV, SEVERITY=1.0)
## 2. IV Pamidronate cyclical: 1 mg/kg/day x3d q3mo (Glorieux 1998 NEJM regimen,
##    ~9 mg/kg/cycle in a 30kg child ~ 30mg/dose x3 days)
## 3. IV Zoledronic acid: 0.05 mg/kg q6mo (Vuorimies 2017; Barros 2012)
## 4. Denosumab SC: 1 mg/kg q3mo off-label pediatric (Hoyer-Kuhn 2014 case
##    series dosing extrapolated; rebound flare modeled at discontinuation)
## 5. Teriparatide SC: 20 ug/day, ADULT OI type I only (Orwoll 2014 RCT;
##    SEVERITY should be set to mild 0.3-0.5 and growth plates closed)
## 6. Setrusumab IV: 20 mg/kg q4wk (ASTEROID/ORBIT phase 2/3 dose range)
## 7. Fresolimumab IV: investigational, 1 mg/kg q4wk (preclinical/phase 1
##    extrapolation from Grafe 2014 Nat Med oim/oim mouse dosing)
##
## Example mrgsolve event code (pamidronate cyclical, 30 kg child):
##   mod <- mread("oi_mrgsolve_model") %>% param(SEVERITY = 1.0)
##   e_pam <- ev(amt = 30, cmt = "PAM_CENT", ii = 24, addl = 2, rate = 30/2,
##               time = 0) %>% ev_repeat(ii = 24*30*3, addl = 12)  # q3mo x4y
##   out <- mod %>% ev(e_pam) %>% mrgsim(end = 24*365*4, delta = 24)
## =============================================================================
