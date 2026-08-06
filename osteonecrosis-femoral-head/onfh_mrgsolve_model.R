## =============================================================================
##  Osteonecrosis of the Femoral Head (ONFH)
##  Quantitative Systems Pharmacology model -- mrgsolve implementation, 49 ODEs
## =============================================================================
##
##  ORGANISING IDEA
##  ---------------
##  The femoral head is a closed compartment in which THREE CLOCKS run at once.
##
##    CLOCK 1  PERFUSION   hours to days.  Marrow adipogenesis and micro-
##       thrombosis raise intraosseous pressure until the Starling resistor
##       closes; the head infarcts.  This clock has stopped before the patient
##       has a symptom, let alone a diagnosis.  Everything aimed at it -- statin,
##       anticoagulant, prostacyclin, and the venting half of a core
##       decompression -- can therefore only work as PROPHYLAXIS.
##
##    CLOCK 2  REPAIR      months.  Creeping substitution.  Osteoclasts excavate
##       the dead trabeculae before osteoblasts refill them, so the interface
##       between dead and living bone is MECHANICALLY WEAKER while it heals than
##       it was while it was merely dead.  And because the front keeps moving,
##       the active zone is continuously reset to mid-cycle: the weakness lasts
##       as long as the front is still crossing the lesion.
##
##    CLOCK 3  FATIGUE     months to years.  Microdamage.  Living bone erases it
##       by targeted remodelling, which is signalled by osteocytes.  Necrotic
##       bone has none, so its damage is permanent.  This clock alone sets the
##       endpoint.
##
##  Collapse is a clock-3 event whose RATE is set by clock 2 through exactly one
##  number: the load-bearing bone volume fraction of the reparative interface,
##
##       beta = 1 - CAV - 0.65 * NB * (1 - MINZ)
##
##  The necrotic bone itself never fails.  It is fully mineralised and as stiff
##  as it ever was.  The hip collapses BECAUSE IT STARTED TO HEAL.
##
##  THE MECHANICAL CORE
##  -------------------
##  The lesion is a cone with its apex at the centre of the head, half-angle
##  alpha, axis tilted phi from the hip joint resultant force.  Two areas follow
##  from that and nothing else:
##
##      loaded cap        A_cap = 2 pi R^2 (1 - cos alpha)     ~ alpha^2
##      conical interface A_int =   pi R^2   sin alpha         ~ alpha
##
##  Every newton that lands on the necrotic cap has to LEAVE through the conical
##  interface, because the apex is a point.  Hence
##
##      sigma_int = L_eff / A_int          (an equilibrium traction, NOT a
##                                          contact stress -- it does not get
##                                          the load-spreading factor that the
##                                          contact probes get)
##      S_int     = S_tr0 * beta^2
##
##  and the whole of ONFH staging is that quotient.  Kerboull's combined
##  necrotic angle measures alpha (CNA = 4 alpha for a circular cone).  The JIC
##  types measure WHERE the cap sits relative to the acetabular contact patch,
##  which is what sets L_eff.  Steinberg's percent-volume measures a mixture of
##  the two, which is why it predicts collapse worse than either.
##
##  WHAT IS FITTED
##  --------------
##  Three numbers, each flagged FITTED in $PARAM:
##     k_dmg    microdamage rate constant  -- ONE anchor, the 5-year collapse
##                                            rate of untreated JIC C1 hips
##     h0_tha   baseline arthroplasty hazard
##     k_pain   pain scale
##  Everything else is geometry, material properties, or published physiology.
##  The JIC A / B / C2 collapse rates, the Kerboull threshold, the timing of the
##  hazard peak, the ceiling on bisphosphonate benefit and the sign flip of core
##  decompression are all PREDICTIONS.
##
##  VERIFICATION
##  ------------
##  Every equation below is also implemented, independently, in
##  onfh_python_reference.py (numpy/scipy).  The cross-check found and fixed
##  TEN substantive defects -- among them a crescent term with a non-zero floor
##  that collapsed every hip including type A, an oedema-pressure loop with gain
##  above one that let a NORMAL femoral head infarct itself (and so made every
##  steroid schedule produce the identical maximal lesion), and incident
##  scenarios that never fed the necrotic volume back into the lesion geometry.
##  All ten are listed in README.md.  NOTE: no R
##  toolchain was available in the environment where this model was written, so
##  this file mirrors the executed Python reference equation for equation but
##  has not itself been run.
##
##  USAGE
##  -----
##    library(mrgsolve); library(dplyr)
##    mod <- mread("onfh_mrgsolve_model.R")
##
##    ## an established JIC C1 lesion, five years, untreated
##    mod %>% param(alpha_deg = 45, phi_deg = -12, lesion_preset = 1) %>%
##            mrgsim(end = 1826, delta = 1) %>% plot(DEPR + sr_int + beta_eff ~ time)
##
##    ## alendronate 70 mg weekly for two years (see SCENARIOS at the bottom)
##    ev_aln <- ev(amt = 16.0, cmt = "ALN_c", ii = 7, addl = 103)
##    mod %>% param(alpha_deg = 45, phi_deg = -12) %>% mrgsim(ev_aln, end = 1826)
##
##  Time unit: DAYS.  Length: mm.  Force: N.  Stress: MPa.  Pressure: mmHg.
## =============================================================================

$PROB
# ONFH: osteonecrosis of the femoral head, 49-ODE QSP model

$PARAM @annotated
// ---- anatomy and materials (literature; not fitted) -----------------------
R_head      : 22.5   : femoral head radius (mm)
t_sp        : 1.10   : subchondral bone plate thickness (mm)
E_sp0       : 5000   : apparent modulus of the subchondral plate (MPa)
nu          : 0.30   : Poisson ratio (-)
E_tr0       : 620    : apparent modulus of femoral head trabecular bone (MPa)
S_tr0       : 9.5    : compressive strength of the same at BV/TV = 1 (MPa)
S_pl0       : 130    : bending strength of the subchondral plate (MPa)
h_found     : 6.0    : Winkler foundation depth (mm)
theta_c     : 50.0   : acetabular contact half-angle (deg)
BW          : 70.0   : body weight (kg)
f_hip       : 2.80   : peak hip contact force (x body weight)
Ncyc        : 5000   : gait cycles per hip per day (1/d)
m_fat       : 4.0    : structural fatigue exponent; specimens give 12-16 and a structure that sheds load from failing elements is flatter (-)
f_spread    : 0.30   : contact-to-trabecular stress spreading factor (-)

// ---- lesion geometry (per patient) ---------------------------------------
alpha_deg   : 45.0   : lesion cone half-angle (deg); Kerboull CNA = 4*alpha
phi_deg     : 0.0    : lesion axis tilt, positive = lateral (deg)
lesion_preset : 1    : 1 = lesion established at t=0; 0 = incident model
XF0         : 0.5    : front penetration already achieved at presentation (mm)

// ---- glucocorticoid PK (prednisolone equivalents) -------------------------
GC_F        : 0.80   : oral bioavailability (-)
GC_ka       : 48.0   : absorption rate (1/d)
GC_CL       : 210    : clearance (L/d)
GC_V        : 35.0   : volume (L)
GC_ke0      : 0.15   : marrow transcriptional integration rate (1/d)
GC_EC50     : 0.220  : EC50 on adipogenic drive; placed at the exposure a 40 mg/d schedule produces (mg/L)
GC_gam      : 2.5    : Hill coefficient; >1 is why PEAK dose matters (-)
GC_Emax     : 3.2    : maximum fold rise of adipogenic drive (-)

// ---- alendronate ----------------------------------------------------------
ALN_CL      : 110    : renal clearance (L/d)
ALN_V       : 28.0   : central volume (L)
ALN_kupt    : 3.4    : plasma to bone surface (1/d)
ALN_kbur    : 0.0035 : bone surface to buried, t1/2 ~ 200 d (1/d)
ALN_koff    : 0.010  : desorption from bone surface (1/d)
ALN_IC50    : 0.55   : IC50 on osteoclast activity (umol/kg)
ALN_Imax    : 0.85   : maximum osteoclast inhibition (-)

// ---- statin, enoxaparin, iloprost, teriparatide, denosumab ---------------
STA_kel     : 0.9    : statin exposure-index turnover (1/d)
STA_IC50    : 0.45   : IC50 on PPARgamma drive (-)
STA_Imax    : 0.55   : maximum inhibition of adipogenic drive (-)
STA_Erunx   : 0.45   : BMP-2 mediated RUNX2 induction (-)
ENX_kel     : 1.8    : enoxaparin elimination (1/d)
ENX_EC50    : 0.35   : EC50 on thrombus dissolution (IU/mL)
ENX_Emax    : 2.2    : maximum effect (-)
ILO_kel     : 0.55   : iloprost effect-site decay (1/d)
ILO_Eendo   : 0.55   : endothelial effect (-)
ILO_Eedema  : 1.30   : oedema resolution effect (-)
TPT_ka      : 26.0   : teriparatide absorption (1/d)
TPT_CL      : 1900   : teriparatide clearance (L/d)
TPT_V       : 100    : teriparatide volume (L)
TPT_EC50    : 12.0   : EC50 on the anabolic index (pg/mL)
TPT_Eob     : 1.35   : osteoblast activity effect (-)
TPT_Eocl    : 0.35   : coupled osteoclast rise (-)
TPT_Emin    : 0.55   : mineralisation rate effect (-)
DMB_kel     : 0.028  : denosumab elimination, t1/2 ~ 25 d (1/d)
DMB_IC50    : 0.55   : IC50 on RANKL (ug/mL)
DMB_Imax    : 0.92   : maximum RANKL inhibition (-)

// ---- marrow lineage -------------------------------------------------------
kin_pparg   : 0.10   : PPARgamma production (1/d)
kout_pparg  : 0.10   : PPARgamma loss (1/d)
kin_runx2   : 0.12   : RUNX2 production (1/d)
kout_runx2  : 0.12   : RUNX2 loss (1/d)
pparg_runx2_I : 0.60 : reciprocal lineage inhibition (-)
k_adipo     : 0.020  : adipocyte expansion (1/d)
ADIPO0      : 0.55   : baseline marrow adipocyte volume fraction (-)
ADIPO_max   : 0.90   : maximum adipocyte volume fraction (-)
k_lipo      : 0.020  : lipolysis back to baseline (1/d)
ETOH        : 0.0    : alcohol adipogenic drive, 0-2 (-)
k_msc       : 0.05   : stromal pool turnover (1/d)
k_opro      : 0.06   : osteoprogenitor turnover (1/d)

// ---- vascular and coagulation --------------------------------------------
kin_pai1    : 0.25   : PAI-1 production (1/d)
kout_pai1   : 0.25   : PAI-1 loss (1/d)
PAI1_gc     : 1.40   : glucocorticoid induction of PAI-1 (-)
THROMBOPHILIA : 0.0  : inherited hypofibrinolysis multiplier, 0-2 (-)
k_thr_on    : 0.045  : thrombus formation (1/d)
k_thr_off   : 0.030  : thrombus lysis (1/d)
kin_endo    : 0.12   : endothelial competence turnover (1/d)
ENDO_gc     : 0.45   : glucocorticoid impairment of eNOS (-)
k_edema_on  : 0.10   : oedema formation (1/d)
k_edema_off : 0.035  : oedema resolution (1/d)
k_angio     : 0.012  : neovascularisation of the front (1/d)

// ---- compartment and perfusion -------------------------------------------
PIO0        : 16.0   : normal femoral head marrow pressure (mmHg)
k_pio_fat   : 55.0   : pressure per unit normalised fat excess (mmHg)
k_pio_ed    : 25.0   : pressure per unit oedema fraction; above ~31 the oedema-pressure loop has gain > 1 (mmHg)
k_pio_thr   : 30.0   : pressure per unit venous occlusion (mmHg)
tau_pio     : 3.0    : pressure equilibration time (d)
k_cd        : 0.55   : venting through a patent channel (1/d)
PIO_cd      : 9.0    : pressure a patent channel vents to (mmHg)
k_close     : 0.0115 : channel closure, t1/2 ~ 60 d (1/d)
P_art       : 65.0   : subchondral arteriolar pressure (mmHg)
P_ven       : 12.0   : venous pressure (mmHg)
R_vasc      : 1.0    : normalised vascular resistance (-)
thr_R       : 4.0    : resistance rise per unit occlusion (-)
pO2_0       : 35.0   : normal marrow tissue oxygen tension (mmHg)
pO2_exp     : 0.75   : flow-to-pO2 exponent (-)

// ---- viability ------------------------------------------------------------
pO2_crit_hem : 15.0  : haematopoietic marrow critical pO2 (mmHg)
n_hem       : 4.0    : Hill coefficient (-)
k_hem       : 3.0    : haematopoietic death rate (1/d)
pO2_crit_ocy : 9.0   : osteocyte critical pO2 (mmHg)
n_ocy       : 4.0    : Hill coefficient (-)
k_ocy       : 0.35   : osteocyte death rate (1/d)
k_ocy_rec   : 0.05   : recovery of stunned osteocytes (1/d)
k_nec       : 0.28   : irreversible commitment to necrosis (1/d)
NEC_MAX     : 0.42   : largest infarct the head geometry allows; alpha 80.8 deg, Kerboull CNA 323 (-)

// ---- reparative front -----------------------------------------------------
v_front     : 0.0247 : creeping substitution ~0.75 mm/month (mm/d)
L_eng       : 3.0    : front depth at which the interface is fully reparative (mm)
W_front     : 2.0    : thickness of the active reparative zone (mm)
k_resp      : 0.055  : osteoclastic excavation at the interface (1/d)
RESP_max    : 0.55   : maximum total turnover of interface bone volume (-)
k_fill      : 0.030  : osteoblastic refill of cavities (1/d)
k_minz      : 0.020  : primary mineralisation, t1/2 ~ 35 d (1/d)
k_fibr_on   : 0.045  : fibrous tissue formation (1/d)
k_fibr_off  : 0.018  : fibrous tissue replacement by bone (1/d)
FIBR_max    : 0.30   : maximum fibrous fraction (-)
k_rim       : 0.010  : reactive sclerosis outside the lesion (1/d)

// ---- remodelling regulators ----------------------------------------------
kin_rankl   : 0.20   : RANKL production (1/d)
kout_rankl  : 0.20   : RANKL loss (1/d)
RANKL_nec   : 2.4    : necrotic debris drive on RANKL (-)
kin_opg     : 0.15   : OPG production (1/d)
kout_opg    : 0.15   : OPG loss (1/d)
OPG_runx2   : 0.60   : RUNX2 drive on OPG (-)
k_ocl       : 0.30   : osteoclast activity turnover (1/d)
k_obl       : 0.22   : osteoblast activity turnover (1/d)

// ---- structure, damage, collapse -----------------------------------------
k_dmg       : 0.90   : FITTED microdamage rate constant; the ONE outcome-fitted number (1/d)
k_dmg_rep   : 0.0015 : targeted removal of microdamage, t1/2 ~ 460 d (1/d)
c_Dself     : 0.8    : crack coalescence; kept weak so that time-to-collapse carries information (-)
zeta_max    : 0.25   : maximum lateral buttress load shedding (-)
dtheta_sh   : 20.0   : angle over which the buttress engages (deg)
chi_rim0    : 1.35   : acetabular edge stress concentration (-)
c_rim_soft  : 0.55   : extra edge loading when the rim sits on the plug (-)
c_depr      : 0.55   : incongruity pressure rise per mm of depression (1/mm)
k_cresc     : 0.030  : crescent propagation (1/d)
k_depr      : 0.055  : head depression rate (1/d)
DEPR_max    : 8.0    : maximum depression (mm)
COLLAPSE_MM : 2.0    : radiographic collapse threshold, ARCO IIIA/IIIB (mm)
k_subpl     : 0.010  : plate degradation once damaged (1/d)

// ---- joint and clinical ---------------------------------------------------
k_cart      : 0.0022 : cartilage loss rate (1/d)
k_synv_on   : 0.030  : synovitis onset (1/d)
k_synv_off  : 0.020  : synovitis resolution (1/d)
k_pain      : 0.16   : FITTED pain scale (-)
tau_pain    : 12.0   : pain equilibration time (d)
h0_tha      : 4.20e-4 : FITTED arthroplasty hazard scale; the hazard is PROPORTIONAL to disease, so a spherical painless hip accrues none (1/d)
a_pain      : 1.5    : pain term in the arthroplasty hazard (-)
a_depr      : 2.6    : depression term in the arthroplasty hazard (-)
a_cart      : 2.2    : cartilage term in the arthroplasty hazard (-)

// ---- interventions --------------------------------------------------------
CD_time     : -1     : core decompression time; <0 = never (d)
CD_ntrack   : 1      : 1 = single 8-10 mm track, 3 = multiple small drillings
BMAC_dose   : 0.0    : bone marrow concentrate, x10^4 CFU-F
SUPPORT_GRAFT : 0.0  : structural support, 0-1 (tantalum rod / fibular graft)
PWB         : 1.0    : protected weight bearing multiplier on gait cycles (-)

$CMT @annotated
GC_gut  : prednisolone-equivalent, gut depot (mg)
GC_cen  : prednisolone-equivalent, central (mg)
GC_eff  : marrow transcriptional exposure index (mg/L-eq)
ALN_c   : alendronate, plasma (ug/L)
ALN_s   : alendronate, bone SURFACE, the active pool (umol/kg)
ALN_d   : alendronate, buried in bone, inactive (umol/kg)
STA_c   : statin exposure index (-)
ENX     : enoxaparin anti-Xa (IU/mL)
ILO     : iloprost effect-site (-)
TPT_sc  : teriparatide subcutaneous depot (ug)
TPT_c   : teriparatide plasma (pg/mL)
DMB_c   : denosumab plasma (ug/mL)
PPARG   : adipogenic transcriptional drive (-)
RUNX2   : osteogenic transcriptional drive (-)
ADIPO   : marrow adipocyte volume fraction (-)
MSC     : marrow stromal progenitor pool (-)
OPRO    : osteoprogenitors available at the front (-)
PAI1    : plasminogen activator inhibitor-1 (-)
THROMB  : microvascular occlusion fraction (-)
ENDO    : endothelial and NO competence (-)
EDEMA   : marrow interstitial oedema fraction (-)
ANGIO   : neovascular density at the reparative front (-)
PIO     : intraosseous pressure of the head (mmHg)
CDCH    : core-decompression channel patency (-)
HEMV    : haematopoietic marrow viability (-)
OCYV    : osteocyte viability in the at-risk zone (-)
NECF    : necrotic volume fraction of the head (-)
XF      : creeping-substitution penetration depth (mm)
CAV     : OPEN resorption cavity fraction at the interface (-)
NB      : new bone that has refilled cavities (-)
MINZ    : mineralisation degree of that new bone (-)
RIMD    : reactive sclerotic rim density index (-)
FIBR    : fibrous / granulation fraction at the interface (-)
RANKL   : RANKL (-)
OPG     : osteoprotegerin (-)
OCL     : osteoclast activity at the interface (-)
OBL     : osteoblast activity at the interface (-)
SUBPL   : subchondral plate integrity (-)
D_int   : microdamage, reparative interface (-)
D_nec   : microdamage, untouched necrotic bone (-)
D_plate : microdamage, subchondral plate at the lesion rim (-)
D_rim   : microdamage under the acetabular rim (-)
D_liv   : microdamage, adjacent living bone (negative control) (-)
CRESC   : subchondral fracture (crescent) extent (-)
DEPR    : femoral head depression (mm)
CART    : articular cartilage integrity (-)
SYNV    : synovitis and effusion (-)
PAINS   : slow pain state (0-10)
H_THA   : cumulative hazard of arthroplasty (-)

$GLOBAL
#include <cmath>

#define CLAMP(x, lo, hi) ((x) < (lo) ? (lo) : ((x) > (hi) ? (hi) : (x)))
#define DEGR 0.017453292519943295

// fatigue microdamage production: a power-law S-N curve in stress/strength,
// with a weak crack-coalescence term
#define DMGPROD(K, NR, SR, MM, CS, DD) \
  ((K) * (NR) * pow(fmin((SR), 1.6), (MM)) * (1.0 + (CS) * fmin((DD), 1.0)))

// ---------------------------------------------------------------------------
// Static lesion geometry.  Computed once per individual in $MAIN, because it
// depends only on alpha, phi, body weight and the contact half-angle.
// ---------------------------------------------------------------------------
double g_ACAP, g_AINT, g_LNEC, g_LTOT, g_THETAL, g_P0, g_RIMON, g_RIMSOFT,
       g_FLOAD, g_XFMAX;

void onfh_geometry(double R, double alpha_deg, double phi_deg,
                   double thetac_deg, double BW, double fhip)
{
  const double a  = alpha_deg  * DEGR;
  const double ph = phi_deg    * DEGR;
  const double tc = thetac_deg * DEGR;

  // peak contact pressure from equilibrium of the cosine pressure field:
  //   F = 2 pi R^2 p0 (1 - cos^3 theta_c) / 3          [R in metres]
  const double Rm = R * 1e-3;
  const double F  = fhip * BW * 9.80665;
  const double geom = 2.0 * M_PI * Rm * Rm * (1.0 - pow(cos(tc), 3.0)) / 3.0;
  g_P0 = (F / geom) * 1e-6;                                        // MPa

  g_ACAP = 2.0 * M_PI * R * R * (1.0 - cos(a));                    // mm^2
  g_AINT =       M_PI * R * R * sin(a);                            // mm^2

  // integrate p(theta) cos(theta) dA over the contact patch, and over the
  // intersection of the contact patch with the lesion cone
  const int NT = 121, NP = 242;
  const double dth = tc / (NT - 1);
  const double dps = 2.0 * M_PI / NP;
  double Ltot = 0.0, Lnec = 0.0;
  for (int i = 0; i < NT; ++i) {
    const double th = i * dth;
    const double pr = g_P0 * cos(th);
    const double dA = R * R * sin(th);
    const double contrib = pr * cos(th) * dA * dth * dps;
    for (int j = 0; j < NP; ++j) {
      const double ps = j * dps;
      Ltot += contrib;
      const double cosd = cos(th) * cos(ph) + sin(th) * sin(ph) * cos(ps);
      if (cosd >= cos(a)) Lnec += contrib;
    }
  }
  g_LTOT  = Ltot;
  g_LNEC  = Lnec;
  g_FLOAD = (Ltot > 0.0) ? Lnec / Ltot : 0.0;

  g_THETAL   = phi_deg + alpha_deg;                    // lateral boundary
  g_RIMON    = (g_THETAL >= thetac_deg) ? 1.0 : 0.0;
  g_RIMSOFT  = CLAMP((g_THETAL - thetac_deg) / 15.0, 0.0, 1.0);
  g_XFMAX    = fmax(2.0, R * sin(a));
}

$MAIN
if (NEWIND <= 1) {
  onfh_geometry(R_head, alpha_deg, phi_deg, theta_c, BW, f_hip);

  PPARG_0  = 1.0;  RUNX2_0 = 1.0;  ADIPO_0 = ADIPO0;  MSC_0 = 1.0;
  OPRO_0   = 1.0;  PAI1_0  = 1.0;  ENDO_0  = 1.0;     PIO_0 = PIO0;
  HEMV_0   = 1.0;  OCYV_0  = 1.0;  RANKL_0 = 1.0;     OPG_0 = 1.0;
  OCL_0    = 1.0;  OBL_0   = 1.0;  SUBPL_0 = 1.0;     CART_0 = 1.0;

  if (lesion_preset > 0.5) {
    // established lesion: necrosis complete, creeping substitution beginning
    NECF_0  = (1.0 - cos(alpha_deg * DEGR)) / 2.0;
    OCYV_0  = 0.0;
    HEMV_0  = 0.0;
    ANGIO_0 = 0.25;
    XF_0    = XF0;
  }
}

$ODE
// ===========================================================================
// ALGEBRAIC LAYER -- re-evaluated at every derivative call
// ===========================================================================

// ---- drug effects ---------------------------------------------------------
double aln   = fmax(ALN_s, 0.0);
double I_aln = ALN_Imax * aln / (ALN_IC50 + aln);
double sta   = fmax(STA_c, 0.0);
double I_sta = STA_Imax * sta / (STA_IC50 + sta);
double E_sta_runx = STA_Erunx * sta / (STA_IC50 + sta);
double enx   = fmax(ENX, 0.0);
double E_enx = ENX_Emax * enx / (ENX_EC50 + enx);
double ilo   = fmax(ILO, 0.0);
double E_ilo_endo = ILO_Eendo  * ilo / (0.5 + ilo);
double E_ilo_ed   = ILO_Eedema * ilo / (0.5 + ilo);
double tpt   = fmax(TPT_c, 0.0);
double E_tpt = tpt / (TPT_EC50 + tpt);
double dmb   = fmax(DMB_c, 0.0);
double I_dmb = DMB_Imax * dmb / (DMB_IC50 + dmb);
double gce   = fmax(GC_eff, 0.0);
double E_gc  = GC_Emax * pow(gce, GC_gam)
             / (pow(GC_EC50, GC_gam) + pow(gce, GC_gam));

// ---- perfusion: Starling resistor / vascular waterfall --------------------
double pio   = fmax(PIO, 0.0);
double p_out = fmax(P_ven, pio);
double p_in  = P_art * (0.55 + 0.45 * CLAMP(ENDO, 0.0, 1.5));
double Rv    = R_vasc * (1.0 + thr_R * CLAMP(THROMB, 0.0, 1.0));
double Qf    = fmax(0.0, (p_in - p_out) / Rv);
// Normalise to the BASELINE operating point, not to venous pressure.  With
// the venous normalisation a healthy hip sat at Qrel = 0.925, which fed the
// oedema term, which raised P_io, which lowered Qrel further: the loop gain
// exceeded one and a normal femoral head slowly infarcted itself.
double Q0    = (P_art - fmax(P_ven, PIO0)) / R_vasc;
double Qrel  = Qf / Q0;
double pO2   = pO2_0 * pow(Qrel, pO2_exp);

// ---- reparative interface material ---------------------------------------
// Volume bookkeeping.  Three phases sum to one:
//   ORIG = 1 - CAV - NB  original dead, fully mineralised trabeculae
//   CAV                  open resorption cavity, no load-bearing value
//   NB                   new bone, relative stiffness 0.35 + 0.65*MINZ
double cav  = CLAMP(CAV,  0.0, 0.95);
double nb   = CLAMP(NB,   0.0, 0.95);
double fibr = CLAMP(FIBR, 0.0, 0.95);
double minz = CLAMP(MINZ, 0.0, 1.0);
double beta_rep = 1.0 - cav - 0.65 * nb * (1.0 - minz) - 0.25 * fibr * cav;
beta_rep = CLAMP(beta_rep, 0.05, 1.25);

// how much of the conical interface has been converted from "dead bone still
// welded to living bone" into reparative tissue
double f_eng    = CLAMP(XF / L_eng, 0.0, 1.0);
double beta_eff = (1.0 - f_eng) + beta_rep * f_eng;
double E_int    = E_tr0 * pow(beta_eff, 2.5);
double S_int    = S_tr0 * pow(beta_eff, 2.0);

// ---- load transfer --------------------------------------------------------
double depr      = fmax(DEPR, 0.0);
double chi_depr  = 1.0 + c_depr * depr;
double dtheta    = theta_c - g_THETAL;
double zeta      = zeta_max * CLAMP(dtheta / dtheta_sh, 0.0, 1.0);
double supp      = CLAMP(SUPPORT_GRAFT, 0.0, 1.0);
double L_eff     = g_LNEC * (1.0 - zeta) * chi_depr * (1.0 - 0.35 * supp);

double cd_defect = 0.0;
if (CD_time >= 0.0 && SOLVERTIME >= CD_time)
  cd_defect = 0.075 / fmax(CD_ntrack, 1.0);
double A_int_eff = g_AINT * (1.0 - cd_defect);

double sigma_int = L_eff / fmax(A_int_eff, 1.0);
double sr_int    = sigma_int / fmax(S_int, 1e-6);

// ---- untouched necrotic bone: strong, but nothing clears its damage -------
double sigma_nec = g_LNEC * chi_depr * f_spread / fmax(g_ACAP, 1.0);
double sr_nec    = sigma_nec / S_tr0;

// ---- subchondral plate bending over the sinking plug ----------------------
double h_eff  = CLAMP(XF, 1.0, 6.0);
double k_seg  = E_int / h_eff;
double w_sink = sigma_int / fmax(k_seg, 1e-6);
double D_pl   = E_sp0 * pow(t_sp, 3.0) / (12.0 * (1.0 - nu * nu));
double k_h    = E_tr0 / h_found;
double lam_h  = pow(4.0 * D_pl / k_h, 0.25);
double lam_s  = pow(4.0 * D_pl / fmax(k_seg, 1e-3), 0.25);
double lam    = sqrt(lam_h * lam_s);
double sigma_pl = E_sp0 * t_sp * w_sink / (lam * lam);
double sr_pl    = sigma_pl / fmax(S_pl0 * CLAMP(SUBPL, 0.02, 1.0), 1e-6);

// ---- acetabular rim -------------------------------------------------------
double chi_rim   = chi_rim0 + c_rim_soft * g_RIMSOFT;
double sigma_rim = g_P0 * cos(theta_c * DEGR) * chi_rim * chi_depr * f_spread;
double beta_rim  = (g_RIMON > 0.5) ? beta_eff : 1.0;
double sr_rim    = sigma_rim / fmax(S_tr0 * beta_rim * beta_rim, 1e-6);
double rim_living = (g_RIMON > 0.5) ? 0.0 : 1.0;

// ---- living-bone negative control ----------------------------------------
// negative control: bone outside the lesion and outside the collapsing
// segment, so it does not see the incongruity concentration
double sr_liv = g_P0 * f_spread / S_tr0;

// ---- remodelling availability --------------------------------------------
double REM_int = fmin(CLAMP(OCL, 0.0, 5.0) / 2.5, 1.2) * f_eng;

// ===========================================================================
// DIFFERENTIAL EQUATIONS
// ===========================================================================

// ---- PK -------------------------------------------------------------------
dxdt_GC_gut = -GC_ka * GC_gut;
dxdt_GC_cen =  GC_F * GC_ka * GC_gut - (GC_CL / GC_V) * GC_cen;
dxdt_GC_eff =  GC_ke0 * (GC_cen / GC_V - GC_eff);

dxdt_ALN_c  = -(ALN_CL / ALN_V) * ALN_c - ALN_kupt * ALN_c
              + ALN_koff * ALN_s * 10.0;
dxdt_ALN_s  =  ALN_kupt * ALN_c / 10.0 - ALN_kbur * ALN_s - ALN_koff * ALN_s;
dxdt_ALN_d  =  ALN_kbur * ALN_s;

dxdt_STA_c  = -STA_kel * STA_c;
dxdt_ENX    = -ENX_kel * ENX;
dxdt_ILO    = -ILO_kel * ILO;
dxdt_TPT_sc = -TPT_ka * TPT_sc;
dxdt_TPT_c  =  TPT_ka * TPT_sc * 1e6 / TPT_V / 1000.0 - (TPT_CL / TPT_V) * TPT_c;
dxdt_DMB_c  = -DMB_kel * DMB_c;

// ---- marrow lineage -------------------------------------------------------
double drive_ad = 1.0 + E_gc + ETOH;
dxdt_PPARG = kin_pparg * drive_ad * (1.0 - I_sta) - kout_pparg * PPARG;
dxdt_RUNX2 = kin_runx2 * (1.0 + E_sta_runx + TPT_Eob * E_tpt)
             / (1.0 + pparg_runx2_I * fmax(PPARG - 1.0, 0.0))
             - kout_runx2 * RUNX2;
dxdt_ADIPO = k_adipo * fmax(PPARG - 1.0, 0.0) * (ADIPO_max - ADIPO)
             - k_lipo * (ADIPO - ADIPO0);
dxdt_MSC   = k_msc * (1.0 - MSC) - 0.02 * fmax(PPARG - 1.0, 0.0);
dxdt_OPRO  = k_opro * (MSC * RUNX2 - OPRO);

// ---- vascular and coagulation --------------------------------------------
dxdt_PAI1 = kin_pai1 * (1.0 + PAI1_gc * E_gc + THROMBOPHILIA) - kout_pai1 * PAI1;
double thr = CLAMP(THROMB, 0.0, 1.0);
dxdt_THROMB = k_thr_on * fmax(PAI1 - 1.0, 0.0) * (1.0 - thr)
              * (2.0 - CLAMP(ENDO, 0.0, 1.0))
              - k_thr_off * (1.0 + E_enx) * thr;
dxdt_ENDO = kin_endo * (1.0 + E_ilo_endo + 0.35 * E_sta_runx
                        - ENDO_gc * E_gc - ENDO);
double ed = CLAMP(EDEMA, 0.0, 1.0);
// oedema needs a real perfusion deficit, not any deficit at all
double isch_ed = fmax(1.0 - CLAMP(Qrel, 0.0, 1.0) - 0.10, 0.0) / 0.90;
dxdt_EDEMA = k_edema_on * isch_ed * (0.55 - ed)
             - k_edema_off * (1.0 + E_ilo_ed) * ed;
double isch = 1.0 - CLAMP(Qrel, 0.0, 1.0);
dxdt_ANGIO = k_angio * CLAMP(NECF / 0.01, 0.0, 1.0) * (1.0 + 0.8 * OPRO)
             * (1.0 - CLAMP(ANGIO, 0.0, 1.0))
             - 0.004 * CLAMP(ANGIO, 0.0, 1.0) * isch;

// ---- compartment ----------------------------------------------------------
double fat_excess = (ADIPO - ADIPO0) / fmax(1.0 - ADIPO0, 1e-3);
double pio_t = PIO0 + k_pio_fat * fmax(fat_excess, 0.0) + k_pio_ed * ed
             + k_pio_thr * thr;
dxdt_PIO  = (pio_t - PIO) / tau_pio
            - k_cd * CLAMP(CDCH, 0.0, 1.0) * (PIO - PIO_cd);
dxdt_CDCH = -k_close * CDCH;

// ---- viability ------------------------------------------------------------
double h_hem = 1.0 / (1.0 + pow(pO2 / pO2_crit_hem, n_hem));
double h_ocy = 1.0 / (1.0 + pow(pO2 / pO2_crit_ocy, n_ocy));
dxdt_HEMV = -k_hem * h_hem * HEMV + 0.08 * (1.0 - h_hem) * (1.0 - HEMV);
dxdt_OCYV = -k_ocy * h_ocy * OCYV
            + k_ocy_rec * (1.0 - h_ocy) * (1.0 - NECF / 0.5) * (1.0 - OCYV);
if (lesion_preset < 0.5) {
  // h_ocy is the fraction of osteocyte territory below the critical oxygen
  // tension, so NEC_MAX * h_ocy is the volume fraction the current depth of
  // ischaemia can commit.  The ratchet makes commitment irreversible:
  // reperfusion does not resurrect an infarct.  Making the CEILING depend on
  // severity, rather than the RATE, is what lets a mild insult produce a
  // small lesion instead of the same maximal one.
  dxdt_NECF = k_nec * fmax(NEC_MAX * h_ocy - CLAMP(NECF, 0.0, NEC_MAX), 0.0);
} else {
  dxdt_NECF = 0.0;
}

// ---- reparative front -----------------------------------------------------
double ocl = CLAMP(OCL,   0.0, 6.0);
double obl = CLAMP(OBL,   0.0, 6.0);
double ang = CLAMP(ANGIO, 0.0, 1.0);
double adv = CLAMP(1.0 - XF / g_XFMAX, 0.0, 1.0);
dxdt_XF = v_front * ang * (0.25 + 0.75 * ocl) * adv;

double exc  = k_resp * ocl * ang * adv * fmax(RESP_max - cav - nb, 0.0);
double fill = k_fill * obl * CLAMP(OPRO, 0.0, 4.0)
              * (1.0 - 0.8 * CLAMP(FIBR, 0.0, 1.0)) * cav;
// RENEWAL.  Creeping substitution is a MOVING front: it leaves finished bone
// behind and invades fresh necrotic bone ahead, so the active zone is
// continuously reset to "not yet excavated".  Without this term the interface
// is a single synchronised remodelling cohort that heals in weeks and the
// disease disappears.
double nu_ren = dxdt_XF / W_front;
dxdt_CAV  = exc - fill - nu_ren * cav;
dxdt_NB   = fill - nu_ren * nb;
dxdt_MINZ = k_minz * (1.0 + TPT_Emin * E_tpt) * (1.0 - CLAMP(MINZ, 0.0, 1.0))
            - nu_ren * CLAMP(MINZ, 0.0, 1.0);
dxdt_RIMD = k_rim * ang * (1.0 - CLAMP(RIMD, 0.0, 1.0));
dxdt_FIBR = k_fibr_on * ang * (1.0 - 0.5 * CLAMP(OPRO, 0.0, 2.0))
            * (FIBR_max - CLAMP(FIBR, 0.0, 1.0))
            - k_fibr_off * (0.2 + nb) * CLAMP(FIBR, 0.0, 1.0);

// ---- remodelling regulators ----------------------------------------------
double nec_debris = CLAMP(NECF / 0.20, 0.0, 1.5) * f_eng;
dxdt_RANKL = kin_rankl * (1.0 + RANKL_nec * nec_debris) * (1.0 - I_dmb)
             - kout_rankl * RANKL;
dxdt_OPG   = kin_opg * (1.0 + OPG_runx2 * fmax(RUNX2 - 1.0, 0.0))
             - kout_opg * OPG;
double ratio = CLAMP(RANKL, 0.01, 20.0) / CLAMP(OPG, 0.05, 20.0);
double ocl_t = ratio * (1.0 - I_aln) * (1.0 + TPT_Eocl * E_tpt);
dxdt_OCL = k_ocl * (ocl_t - OCL);
// coupling: osteoblast recruitment follows resorbed surface, which is why an
// antiresorptive cannot be a pure benefit -- it removes the stimulus for the
// very formation it is trying to protect
double coup  = 0.30 + 0.70 * CLAMP(CAV, 0.0, 1.0) / fmax(RESP_max, 1e-3);
double obl_t = CLAMP(RUNX2, 0.0, 5.0) * CLAMP(OPRO, 0.0, 5.0) * coup
               * (1.0 + TPT_Eob * E_tpt);
dxdt_OBL = k_obl * (obl_t - OBL);

// ---- microdamage ----------------------------------------------------------
// Production follows a power-law S-N curve in stress/strength.  Removal is
// TARGETED remodelling, which requires living osteocytes to signal it and is
// therefore identically zero inside necrotic bone.
double Nrel = (Ncyc * PWB) / 5000.0;
#define DP(SR, DD) DMGPROD(k_dmg, Nrel, (SR), m_fat, c_Dself, (DD))

dxdt_D_int   = DP(sr_int, fmax(D_int, 0.0))
               - k_dmg_rep * REM_int * fmax(D_int, 0.0);
dxdt_D_nec   = DP(sr_nec, fmax(D_nec, 0.0));            // NO removal term:
                                                        // no osteocytes, no
                                                        // targeted remodelling
dxdt_D_plate = DP(sr_pl,  fmax(D_plate, 0.0))
               - k_dmg_rep * 0.25 * REM_int * fmax(D_plate, 0.0);
dxdt_D_rim   = DP(sr_rim, fmax(D_rim, 0.0))
               - k_dmg_rep * rim_living * fmax(D_rim, 0.0);
dxdt_D_liv   = DP(sr_liv, fmax(D_liv, 0.0))
               - k_dmg_rep * 1.0 * fmax(D_liv, 0.0);

// ---- collapse -------------------------------------------------------------
double Dmax = fmax(fmax(D_int, D_nec), fmax(D_plate, D_rim));
Dmax = fmax(Dmax, 0.0);
double phi_fail = (Dmax > 0.0) ? pow(Dmax, 8.0) / (1.0 + pow(Dmax, 8.0)) : 0.0;
double cr = CLAMP(CRESC, 0.0, 1.0);
dxdt_CRESC = k_cresc * fmax(phi_fail - cr, 0.0) * 4.0;
dxdt_DEPR  = k_depr * cr * (0.35 + sr_int) * (DEPR_max - fmax(DEPR, 0.0));
dxdt_SUBPL = -k_subpl * cr * CLAMP(SUBPL, 0.0, 1.0);

// ---- joint and clinical ---------------------------------------------------
// Cartilage is lost because the head has become incongruent, not because time
// passes.  A baseline term here silently destroyed healthy hips.
dxdt_CART = -k_cart * cr * (1.0 + depr / 2.0) * CLAMP(CART, 0.0, 1.0);
dxdt_SYNV = k_synv_on * (0.5 * cr + 0.5 * (1.0 - CART))
            * (1.0 - CLAMP(SYNV, 0.0, 1.0))
            - k_synv_off * CLAMP(SYNV, 0.0, 1.0);
double pain_t = 10.0 * (1.0 - exp(-k_pain * (
                  0.09 * fmax(PIO - PIO0, 0.0)
                + 3.0  * ed
                + 12.0 * cr
                + 2.2  * depr
                + 6.0  * (1.0 - CLAMP(CART, 0.0, 1.0))
                + 4.0  * CLAMP(SYNV, 0.0, 1.0))));
dxdt_PAINS = (pain_t - PAINS) / tau_pain;
// Arthroplasty hazard is PROPORTIONAL to disease.  With an additive constant
// an asymptomatic, spherical, painless hip accrued a 34% five-year probability
// of being replaced.
dxdt_H_THA = h0_tha * (a_pain * PAINS / 10.0 + a_depr * depr / 4.0
                       + a_cart * (1.0 - CLAMP(CART, 0.0, 1.0)));

$TABLE
// Every probe is recomputed here EXACTLY as in $ODE, so that captured output
// corresponds to the requested output time rather than to whatever internal
// step the solver last took.
double beta_out = CLAMP(1.0 - CAV - 0.65 * NB * (1.0 - MINZ)
                        - 0.25 * FIBR * CAV, 0.05, 1.25);
double feng_out = CLAMP(XF / L_eng, 0.0, 1.0);
double BETA_EFF = (1.0 - feng_out) + beta_out * feng_out;
double E_INT    = E_tr0 * pow(BETA_EFF, 2.5);
double S_INT    = S_tr0 * pow(BETA_EFF, 2.0);
double CHI      = 1.0 + c_depr * fmax(DEPR, 0.0);
double ZETA     = zeta_max * CLAMP((theta_c - g_THETAL) / dtheta_sh, 0.0, 1.0);
double CDDEF    = (CD_time >= 0.0 && TIME >= CD_time)
                  ? 0.075 / fmax(CD_ntrack, 1.0) : 0.0;
double SIGMA_INT = g_LNEC * (1.0 - ZETA) * CHI
                   * (1.0 - 0.35 * CLAMP(SUPPORT_GRAFT, 0.0, 1.0))
                   / fmax(g_AINT * (1.0 - CDDEF), 1.0);
double SR_INT   = SIGMA_INT / fmax(S_INT, 1e-6);

double SR_NEC   = g_LNEC * CHI * f_spread / fmax(g_ACAP, 1.0) / S_tr0;

double KSEG     = E_INT / CLAMP(XF, 1.0, 6.0);
double WSINK    = SIGMA_INT / fmax(KSEG, 1e-6);
double DPLT     = E_sp0 * pow(t_sp, 3.0) / (12.0 * (1.0 - nu * nu));
double LAMH     = pow(4.0 * DPLT / (E_tr0 / h_found), 0.25);
double LAMS     = pow(4.0 * DPLT / fmax(KSEG, 1e-3), 0.25);
double SR_PL    = E_sp0 * t_sp * WSINK / (LAMH * LAMS)
                  / fmax(S_pl0 * CLAMP(SUBPL, 0.02, 1.0), 1e-6);

double BRIM     = (g_RIMON > 0.5) ? BETA_EFF : 1.0;
double SR_RIM   = g_P0 * cos(theta_c * DEGR)
                  * (chi_rim0 + c_rim_soft * g_RIMSOFT) * CHI * f_spread
                  / fmax(S_tr0 * BRIM * BRIM, 1e-6);
double SR_LIV   = g_P0 * f_spread / S_tr0;

double COLLAPSED = (DEPR >= COLLAPSE_MM) ? 1.0 : 0.0;
double CNA      = 4.0 * alpha_deg;
double THETA_L  = g_THETAL;
double F_LOAD   = g_FLOAD;
double A_CAP    = g_ACAP;
double A_INT    = g_AINT;
double AMPLIF   = g_ACAP / fmax(g_AINT, 1e-6);        // = 2 tan(alpha/2)
double P0       = g_P0;
double HHS      = 44.0 * (1.0 - PAINS / 10.0)
                + 47.0 * CLAMP(1.0 - DEPR / 6.0, 0.0, 1.0) * CLAMP(CART, 0.0, 1.0)
                +  9.0 * CLAMP(1.0 - DEPR / 8.0, 0.0, 1.0);
double P_THA    = 1.0 - exp(-H_THA);
double VAS      = PAINS;

$CAPTURE @annotated
BETA_EFF  : load-bearing bone volume fraction of the interface (-)
E_INT     : interface modulus (MPa)
S_INT     : interface strength (MPa)
SIGMA_INT : interface traction (MPa)
SR_INT    : interface stress / strength (-)
SR_NEC    : necrotic core stress / strength (-)
SR_PL     : subchondral plate stress / strength (-)
SR_RIM    : acetabular rim stress / strength (-)
SR_LIV    : adjacent living bone stress / strength, negative control (-)
COLLAPSED : femoral head depression has reached the radiographic threshold
CNA       : Kerboull combined necrotic angle (deg)
THETA_L   : lateral lesion boundary relative to the load axis (deg)
F_LOAD    : fraction of the joint load landing on the lesion (-)
A_CAP     : loaded cap area (mm2)
A_INT     : conical interface area (mm2)
AMPLIF    : A_cap / A_int = 2 tan(alpha/2) (-)
P0        : peak contact pressure (MPa)
HHS       : Harris Hip Score (0-100)
VAS       : pain, visual analogue scale (0-10)
P_THA     : cumulative probability of arthroplasty (-)

## =============================================================================
##  SCENARIOS
##  ---------
##  All of these are executed in onfh_python_reference.py; the numbers they
##  produce are tabulated in onfh_reference_output.txt and discussed in
##  README.md.  JIC presets are chosen so that theta_L = phi + alpha lands
##  mid-class:  A (28, -55) · B (34, -34) · C1 (45, -12) · C2 (58, +4).
##
##    library(mrgsolve); library(dplyr)
##    mod <- mread("onfh_mrgsolve_model.R")
##    jic <- function(m, a, ph, ...) param(m, alpha_deg = a, phi_deg = ph, ...)
##
##  S01  healthy hip                 jic(mod, 0.001, 0) |> param(lesion_preset=0)
##  S02  high-dose steroid, incident lesion
##         ev(amt = 60, cmt = "GC_gut", ii = 1, addl = 29) then a 90-day taper;
##         param(lesion_preset = 0)
##  S03  same cumulative dose at a low peak: amt = 32.1, ii = 1, addl = 144
##  S04  JIC A untreated             jic(mod, 28, -55)
##  S05  JIC B untreated             jic(mod, 34, -34)
##  S06  JIC C1 untreated            jic(mod, 45, -12)      <- the ONE anchor
##  S07  JIC C2 untreated            jic(mod, 58,   4)
##  S08  C1 + alendronate            ev(amt = 16.0, cmt = "ALN_c", ii = 7, addl = 103)
##  S09  C1 + core decompression     param(CD_time = 30, CD_ntrack = 1, PWB = 0.55)
##  S10  C1 + multiple drilling + BMAC
##                                   param(CD_time = 30, CD_ntrack = 3,
##                                         BMAC_dose = 2, PWB = 0.55)
##  S11  C1 + teriparatide           ev(amt = 20, cmt = "TPT_sc", ii = 1, addl = 547)
##  S12  C1 + alendronate 12 mo then teriparatide 18 mo
##  S13  steroid + rosuvastatin prophylaxis
##                                   ev(amt = 0.9, cmt = "STA_c", ii = 1, addl = 364)
##  S14  steroid + enoxaparin x 12 weeks, THROMBOPHILIA = 0.8
##  S15  JIC C2 + core decompression        <- where the sign flips
##  S16  JIC B  + core decompression
##  S17  C1 + iloprost 5-day infusion
##  S18  C1 + denosumab 60 mg q6m
##  S19  C1 + structural graft       param(SUPPORT_GRAFT = 1, CD_ntrack = 3,
##                                         BMAC_dose = 2, PWB = 0.50)
##  S20  C1 + protected weight bearing alone   param(PWB = 0.45)
##
##  Collapse is DEPR >= COLLAPSE_MM (2 mm), i.e. the ARCO IIIA/IIIB boundary.
## =============================================================================
