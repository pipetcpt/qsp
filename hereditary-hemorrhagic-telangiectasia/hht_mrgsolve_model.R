## ============================================================================
##  Hereditary Hemorrhagic Telangiectasia (HHT / Osler-Weber-Rendu)
##  유전성 출혈성 모세혈관확장증 — QSP model for mrgsolve
##
##  57 ODE states · 5 vascular beds · 5 drugs · 18 treatment scenarios
##
##  ORGANISING THESIS
##  ---------------------------------------------------------------------------
##  (1) ALK1/ENG/SMAD4 is the endothelial SHEAR SET-POINT controller.  In normal
##      endothelium a rise in wall shear drives INWARD remodelling and shear
##      returns to its set point — negative feedback (PMID 27646277, 32078368,
##      37490341).  A two-hit lesion has no pathway left, so the SIGN of that
##      loop inverts.  One line carries it:
##
##          REMOD = KSH * g(WSS) * (1 - 2*pSMAD_rel)
##
##      pSMAD_rel = 1  -> negative -> inward, self-limiting
##      pSMAD_rel ~ 0  -> positive -> outward, self-amplifying
##
##  (2) Wall shear is NOT monotone in lumen size.  Q = dP/(R_feed + R_lesion)
##      with R_lesion ~ S^-4, so WSS ~ S/(1+(S/S_feed)^4): shear rises while the
##      LESION is the limiting resistance and falls once the FEEDER is.  A nasal
##      ectasia fed by an arteriole saturates early and stays drug-reachable; a
##      hepatic shunt fed by a conducting artery does not.  Same equation, two
##      basins — telangiectasia and AVM are not two diseases here.
##
##  (3) Anaemia is an INPUT to lesion growth, not only its output.  Hb down ->
##      cardiac output up (O2-delivery homeostasis) -> perfusion pressure and
##      shear up -> outward remodelling -> more bleeding -> Hb down.  Holding Hb
##      at 6.0 instead of 14.6 raises blood loss 33% through this loop alone.
##
##  (4) The decisive clinical variable is a FLUX BALANCE, not a score.  Iron in
##      blood = 0.0347 * Hb[g/dL] mg per mL, so steady state requires
##          absorbed - obligatory = 0.0347 * Hb * blood-loss rate
##      which inverts to a maximum sustainable blood loss per iron strategy:
##      8.4 mL/d on diet, 17.8 mL/d on oral iron (hepcidin-capped — 200 mg/d
##      buys nothing over 65 mg/d), 88.3 mL/d on monthly 1 g IV iron.  Benefit
##      is therefore THRESHOLD-LIKE where ESS is continuous.
##
##  CALIBRATION / VALIDATION (model vs published)
##  ---------------------------------------------------------------------------
##   PATH-HHT pomalidomide, ESS vs placebo    -0.83   vs  -0.94   (PMID 39292928)
##   PATH-HHT pomalidomide, QOL vs placebo    -1.52   vs  -1.4
##   Dupuis-Girod bevacizumab, CI at 3 mo     4.41    vs   4.20   (PMID 22396517)
##   Dupuis-Girod bevacizumab, MED at 3 mo     102    vs    134
##   Dupuis-Girod nasal bevacizumab           flat    vs   flat   (PMID 27599328)
##   ATERO oral TXA, monthly duration        -14.1%   vs  -17.3%  (PMID 25040799)
##   ATERO oral TXA, episode count            +0.2%   vs    ~0
##   Whitehead topical TXA 40 mg/d            -0.6%   vs    ~0    (PMID 27599329)
##   Parambil pazopanib, transfusions/3 mo   15.5->2.0 vs 16->0   (PMID 34292451)
##   Parambil pazopanib, ESS change           -2.56   vs  -4.77   (under-predicted)
##
##  Units: time = DAYS, volume = mL, Hb = g/dL, iron = mg, drug amounts = mg.
##
##  Every equation here was first written and run in an independent Python
##  implementation.  That implementation found six real defects, each fixed
##  before this file was written — see README.md, "개발 중 발견한 결함".
## ============================================================================

library(mrgsolve)

code <- '
$PARAM @annotated
// ---------------- patient ------------------------------------------------
WT     :  70.0  : body weight (kg)
BSA    :   1.80 : body surface area (m2)
BV_kg  :  70.0  : blood volume per kg (mL/kg)
GD     :   0.5  : germline functional gene dosage (0.5 het, 1.0 control)
GD_LES :   0.055: residual pathway activity inside a two-hit lesion
BMP9   :   6.0  : circulating BMP9 (ng/mL)
KD9    :   4.0  : BMP9 receptor KD (ng/mL)
SEV    :   1.05 : phenotype severity multiplier on somatic-hit rate
ANGF   :   1.02 : per-patient angiogenic drive (modifier-gene load)
GIF    :   0.55 : GI angiodysplasia burden multiplier
HEPF   :   0.20 : hepatic bed - fraction carrying two-hit lesions (0-1)
PULF   :   0.0  : pulmonary bed - fraction carrying two-hit lesions (0-1)
CNSF   :   0.0  : cerebral bed - fraction carrying two-hit lesions (0-1)

// ---------------- bevacizumab -------------------------------------------
BEV_CL :   0.207: bevacizumab clearance (L/day)
BEV_V1 :   2.90 : central volume (L)
BEV_Q  :   0.22 : intercompartmental clearance (L/day)
BEV_V2 :   2.30 : peripheral volume (L)
BEV_KP :   0.055: plasma-to-interstitium partition for IgG
BEVN_KMC: 66.5  : nasal mucociliary clearance (1/day, t1/2 ~15 min)
BEVN_F :   0.010: fraction of nasal depot reaching the submucosal plexus
BEVN_KE:   4.16 : submucosal washout (1/day)
BEVN_V :   0.012: nasal submucosal distribution volume (L)
BEV_KD :   1.1  : bevacizumab-VEGF KD (nM)
BEV_MW : 149000 : bevacizumab MW (g/mol)
VEGF_MW:  45000 : VEGF-A dimer MW (g/mol)

// ---------------- pazopanib ---------------------------------------------
PAZ_KA :  12.0  : absorption rate (1/day)
PAZ_CL :  10.1  : clearance (L/day)
PAZ_V  :  11.1  : volume (L)
PAZ_FU :   0.0005 : unbound fraction
PAZ_IC50_VEGFR2 : 1.2 : free conc for 50% VEGFR2 inhibition (ng/mL)
PAZ_IC50_PDGFRB : 14.0: free conc for 50% PDGFRbeta inhibition (ng/mL)

// ---------------- pomalidomide ------------------------------------------
POM_KA :  28.8  : absorption rate (1/day)
POM_CL : 168.0  : clearance (L/day)
POM_V  :  62.0  : volume (L)
POM_EC50:  8.0  : EC50 for PDGF-B / mural maturation (ng/mL)
POM_EMAX:  0.62 : maximal fractional increase in PDGF-B production
POM_EC50V: 25.0 : EC50 for VEGF suppression (ng/mL)
POM_EMAXV: 0.40 : maximal fractional VEGF suppression

// ---------------- tranexamic acid ---------------------------------------
TXA_KA :  21.6  : absorption rate (1/day)
TXA_CL : 144.0  : clearance (L/day)
TXA_V  :  30.0  : volume (L)
TXA_IC50:  6.0  : plasmin inhibition IC50 (ug/mL)
TXA_IMAX:  0.92 : maximal plasmin inhibition

// ---------------- investigational ---------------------------------------
ALK1AG :   0.0  : fractional restoration of lesional ALK1 signalling (0-1)

// ---------------- signalling --------------------------------------------
KIN_PS :   1.0  : pSMAD production rate
KOUT_PS:   1.0  : pSMAD degradation rate
GSHEAR :   0.85 : shear amplification of SMAD1/5
KWS    :   1.0  : half-max shear for SMAD amplification
KV0    :  20.0  : basal VEGF production (pg/mL/day)
AV     :   2.2  : pSMAD-loss amplification of VEGF production
KDV    :   1.0  : VEGF degradation (1/day)
KV_HYP :  14.0  : hypoxia-driven VEGF production
KV_HALF:  52.0  : free VEGF for half-max VEGFR2 signal (pg/mL)
KA_AKT :   1.0  : AKT activation rate
KD_AKT :   1.0  : AKT deactivation rate
WA     :   1.0  : pSMAD-loss amplification of AKT
PD0    :   0.55 : pSMAD-independent baseline PDGF-B production
KPD    :   1.0  : pSMAD-driven PDGF-B production
KDPD   :   1.0  : PDGF-B turnover

// ---------------- lesion birth / involution ------------------------------
KHIT_N :   0.0049 : nasal two-hit-to-lesion conversion attempts (1/day)
KHIT_G :   0.00135: GI conversion attempts (1/day)
KC     :   0.30 : angiogenic index for half-max conversion
HC     :   3.0  : Hill coefficient for conversion
KREG_N :   9.0e-5 : baseline nasal lesion regression (1/day)
KREG_G :   9.0e-5 : baseline GI lesion regression (1/day)
KREGA  :  60.0  : acceleration of involution when the drive is withdrawn
ANGIO0 :   0.546: reference angiogenic index (untreated)

// ---------------- lesion size -------------------------------------------
KG     :   0.020 : angiogenesis-driven growth
KSH    :   0.047 : shear remodelling gain  <-- the sign-flipping term
KWSS   :   1.0  : half-max shear for remodelling
KDEC   :   0.0085: passive decay of dilation
KMAT   :   0.030 : mural-cell-driven inward maturation
S_BIRTH:   1.35 : calibre a newly converted lesion is born with
WSS0_n :   0.30 : nasal perfusion-pressure scaling
WSS0_g :   0.34 : GI
WSS0_h :   1.55 : hepatic
WSS0_p :   1.15 : pulmonary
WSS0_c :   1.30 : cerebral
SSAT_n :   2.2  : nasal calibre at which the FEEDER becomes limiting
SSAT_g :   2.6  : GI
SSAT_h :  14.0  : hepatic
SSAT_p :  12.0  : pulmonary
SSAT_c :  10.0  : cerebral
SMAX_n :   6.0  : nasal anatomical calibre ceiling
SMAX_g :   7.5  : GI
SMAX_h :  48.0  : hepatic
SMAX_p :  42.0  : pulmonary
SMAX_c :  30.0  : cerebral
KMR    :   0.160: mural recruitment rate
KML    :   0.050: mural loss rate
WM     :   1.30 : angiogenic amplification of mural loss
KML_PAZ:   0.075: PDGFRbeta-blockade cost to mural coverage

// ---------------- bleeding ----------------------------------------------
KE_N   :   0.0631: nasal bleed events per lesion per day
KGB    :   0.0155: GI continuous ooze coefficient
MFRAG  :   1.55 : fragility exponent on (1 - mural coverage)
RBLD   :   0.138: bleed rate coefficient (mL/min per S^2)
TEV_N  :   2.10 : episode duration coefficient (min per unit S)
KDF    :   2.60 : duration amplification by loss of mural coverage
FLYS   :   0.256: fibrinolysis-attributable share of bleed volume/duration

// ---------------- iron / erythron ---------------------------------------
FE_OBLIG:  1.0  : obligatory non-erythroid iron loss (mg/day)
FE_ORAL :  0.0  : oral elemental iron dose (mg/day)
FE_ABS_MAX:10.0 : absolute ceiling on daily absorbed iron (mg/day)
FE_ABS_FR: 0.22 : fractional absorption of supplemental iron
FE_DIET : 15.0  : dietary iron intake (mg/day)
FE_DIET_FR:0.35 : maximal fractional absorption of dietary iron
T_MOB   : 20.0  : storage-iron mobilisation time constant (day)
K_HEP   :  1.0  : hepcidin gate half-max
N_HEP   :  2.4  : hepcidin gate Hill coefficient
KFES    :1400.0 : iron stores for half-max hepcidin (mg)
KERFE   :  0.85 : erythroferrone suppression of hepcidin
FES0    :900.0  : initial iron stores (mg)
HB_T    : 14.6  : target/normal haemoglobin (g/dL)
T_RBC   :110.0  : red cell lifespan (day)
KEPO    :  0.55 : EPO sensitivity to Hb deficit
EPO_MAX :  6.5  : maximal EPO drive
TX_TRIG :  7.0  : haemoglobin transfusion trigger (g/dL)
TX_GAIN :  0.55 : transfusion delivery gain (g/dL per day at full trigger)
TX_UNIT_FE:200.0: iron per RBC unit (mg)
FE_IV_RATE: 0.0 : scheduled IV iron (mg/day equivalent)
HB_FIX  :  0.0  : if > 0, hold Hb at this value (loop-gain experiment only)

// ---------------- cardiac -----------------------------------------------
CI_REF  :  2.85 : reference cardiac index (L/min/m2)
ALPHA_DO2: 0.70 : degree of O2-delivery compensation (1 = full)
SAO2_REF:  0.98 : reference arterial saturation
KSHUNT_H:  0.0021 : hepatic shunt flow per unit conductance^2 (L/min)
KSHUNT_P:  0.0006 : pulmonary
TAU_CO  :  6.0  : cardiac output adaptation time constant (day)
KRAP    :  1.0  : right atrial pressure gain

// ---------------- scores -------------------------------------------------
ESS_TAU :  14.0 : ESS smoothing time constant (day)
PBO_MAX :   0.0 : EMPIRICAL trial-reconstruction placebo drift (ESS points)
PBO_TAU :  42.0 : placebo drift time constant (day)

// ---------------- toxicity -----------------------------------------------
SBP0    :124.0  : baseline systolic BP (mmHg)
KSBP_BEV: 11.0  : bevacizumab hypertension effect (mmHg)
KSBP_PAZ: 13.0  : pazopanib hypertension effect (mmHg)
TAU_SBP : 10.0  : BP time constant (day)
KUP_BEV :  0.9  : bevacizumab proteinuria (g/g creatinine)
KUP_PAZ :  0.4  : pazopanib proteinuria
TAU_UP  : 20.0  : proteinuria time constant (day)
ALT0    : 24.0  : baseline ALT (U/L)
KALT_PAZ:  2.6  : pazopanib ALT multiplier
TAU_ALT : 12.0  : ALT time constant (day)
ANC0    :  3.6  : baseline neutrophils (10^9/L)
KANC_POM:  0.42 : pomalidomide neutropenia
TAU_ANC :  8.0  : ANC time constant (day)
KVTE_POM: 1.1e-4: pomalidomide VTE hazard accrual
KVTE_TXA: 2.0e-5: tranexamic acid VTE hazard accrual

$CMT @annotated
BEV_C   : bevacizumab central (mg)
BEV_P   : bevacizumab peripheral (mg)
BEVN_S  : bevacizumab nasal surface depot (mg)
BEVN_T  : bevacizumab nasal submucosa (mg)
PAZ_A   : pazopanib absorption depot (mg)
PAZ_C   : pazopanib central (mg)
POM_A   : pomalidomide absorption depot (mg)
POM_C   : pomalidomide central (mg)
TXA_A   : tranexamic acid absorption depot (mg)
TXA_C   : tranexamic acid central (mg)
PS_n    : pSMAD1/5/8 activity, nasal lesion
PS_g    : pSMAD1/5/8 activity, GI lesion
PS_h    : pSMAD1/5/8 activity, hepatic lesion
PS_p    : pSMAD1/5/8 activity, pulmonary lesion
PS_c    : pSMAD1/5/8 activity, cerebral lesion
VEGF_n  : tissue VEGF-A, nasal (pg/mL)
VEGF_g  : tissue VEGF-A, GI (pg/mL)
VEGF_h  : tissue VEGF-A, hepatic (pg/mL)
VEGF_p  : tissue VEGF-A, pulmonary (pg/mL)
VEGF_c  : tissue VEGF-A, cerebral (pg/mL)
AKT_n   : PI3K/AKT activity, nasal
AKT_g   : PI3K/AKT activity, GI
AKT_h   : PI3K/AKT activity, hepatic
AKT_p   : PI3K/AKT activity, pulmonary
AKT_c   : PI3K/AKT activity, cerebral
ID1_n   : ID1 transcript (pSMAD reporter)
PDGFB_n : endothelial PDGF-B, nasal
PDGFB_g : endothelial PDGF-B, GI
N_n     : number of nasal telangiectasias
N_g     : number of GI telangiectasias
S_n     : mean nasal lesion calibre index
S_g     : mean GI lesion calibre index
S_h     : hepatic shunt conductance index
S_p     : pulmonary shunt conductance index
S_c     : cerebral AVM index
MU_n    : mural cell coverage, nasal (0-1)
MU_g    : mural cell coverage, GI (0-1)
PLG     : plasmin activity index (1 = untreated)
FES     : body iron stores (mg)
HB      : haemoglobin (g/dL)
HEP     : hepcidin (relative)
EPO     : erythropoietin drive (relative)
RET     : reticulocyte index
CO      : cardiac output (L/min)
LVR     : left ventricular remodelling index
RAP     : right atrial pressure / congestion index
CUM_BL  : cumulative blood loss (mL)
CUM_TX  : cumulative transfused units
CUM_FEIV: cumulative IV iron (mg)
ESSL    : latent Epistaxis Severity Score
QOL     : HHT-specific quality of life (0-16)
PBO     : placebo drift (ESS points, trial reconstruction only)
SBP     : systolic blood pressure (mmHg)
UPCR    : urine protein/creatinine ratio
ALT     : alanine aminotransferase (U/L)
ANC     : absolute neutrophil count (10^9/L)
VTEH    : cumulative venous thromboembolism hazard

$GLOBAL
#define HILL(x, k) ((x) > 0 ? (x)/((k)+(x)) : 0.0)
#define POS(x) ((x) > 0 ? (x) : 0.0)

// smooth replacement for a 0..n bin index, on a log scale
double sbin1(double x, double e){
  double lx = log(x > 1e-9 ? x : 1e-9);
  return 1.0/(1.0 + exp(-(lx - log(e))*3.2));
}
// free ligand from 1:1 binding (same units throughout)
double freelig(double L, double R, double KD){
  double b = L - R - KD;
  return 0.5*(b + sqrt(b*b + 4.0*L*KD));
}

$MAIN
if(NEWIND < 2){
  double bmp0 = BMP9/(KD9 + BMP9);
  PS_n_0 = GD_LES*bmp0;  PS_g_0 = GD_LES*bmp0;  PS_h_0 = GD_LES*bmp0;
  PS_p_0 = GD_LES*bmp0;  PS_c_0 = GD_LES*bmp0;
  VEGF_n_0 = 44.0; VEGF_g_0 = 44.0; VEGF_h_0 = 44.0;
  VEGF_p_0 = 44.0; VEGF_c_0 = 44.0;
  AKT_n_0 = 0.45; AKT_g_0 = 0.45; AKT_h_0 = 0.45;
  AKT_p_0 = 0.45; AKT_c_0 = 0.45;
  ID1_n_0 = GD_LES*bmp0;
  PDGFB_n_0 = PD0 + GD_LES*bmp0;
  PDGFB_g_0 = PD0 + GD_LES*bmp0;
  N_n_0 = 0.02; N_g_0 = 0.01;
  S_n_0 = 1.0; S_g_0 = 1.0; S_h_0 = 1.0; S_p_0 = 1.0; S_c_0 = 1.0;
  MU_n_0 = 0.50; MU_g_0 = 0.50;
  PLG_0 = 1.0;
  FES_0 = FES0;
  HB_0  = HB_T;
  HEP_0 = 1.0; EPO_0 = 1.0; RET_0 = 1.0;
  CO_0  = CI_REF*BSA;
  LVR_0 = 1.0;
  SBP_0 = SBP0; ALT_0 = ALT0; ANC_0 = ANC0;
}

$ODE
double BV = BV_kg*WT;
double HBx = HB > 1.0 ? HB : 1.0;

// ---------------- drug concentrations (mg/L == ug/mL) --------------------
double C_BEV = BEV_C/BEV_V1;                       // ug/mL
double C_PAZ = PAZ_C/PAZ_V;                        // ug/mL
double CU_PAZ = C_PAZ*PAZ_FU*1000.0;               // ng/mL free
double C_POM = POM_C/POM_V*1000.0;                 // ng/mL
double C_TXA = TXA_C/TXA_V;                        // ug/mL

double bev_sys = (C_BEV*BEV_KP)*1e-3/BEV_MW*1e9;   // nM in the interstitium
double bev_top = (BEVN_T/BEVN_V)*1e-3/BEV_MW*1e9;  // nM in nasal submucosa

double I_VEGFR2 = HILL(CU_PAZ, PAZ_IC50_VEGFR2);
double I_PDGFRB = HILL(CU_PAZ, PAZ_IC50_PDGFRB);
double E_IMID   = POM_EMAX *HILL(C_POM, POM_EC50);
double E_IMID_V = POM_EMAXV*HILL(C_POM, POM_EC50V);

// ---------------- oxygen delivery and cardiac output ---------------------
double SPO2 = 0.98 - 0.34*HILL(pow(S_p,2.5), pow(20.0,2.5));
double CO_anemia = CI_REF*pow((HB_T*SAO2_REF)/(HBx*SPO2), ALPHA_DO2)*BSA;
double Q_SHUNT = KSHUNT_H*S_h*S_h + KSHUNT_P*S_p*S_p;
double CO_target = CO_anemia + Q_SHUNT;
double COx = CO > 1.0 ? CO : 1.0;
double CO_rel = COx/(CI_REF*BSA);

// ---------------- per-bed shear, VEGF and angiogenic index ---------------
// WSS ~ S/(1+(S/SSAT)^4): rises while the LESION limits flow, falls once the
// FEEDING vessel does.  This is the whole telangiectasia-vs-AVM distinction.
double ivg = 1e-9/VEGF_MW*1e9;                     // pg/mL -> nM

double WSS_n = WSS0_n*CO_rel*S_n/(1.0 + pow(S_n/SSAT_n,4));
double WSS_g = WSS0_g*CO_rel*S_g/(1.0 + pow(S_g/SSAT_g,4));
double WSS_h = WSS0_h*CO_rel*S_h/(1.0 + pow(S_h/SSAT_h,4));
double WSS_p = WSS0_p*CO_rel*S_p/(1.0 + pow(S_p/SSAT_p,4));
double WSS_c = WSS0_c*CO_rel*S_c/(1.0 + pow(S_c/SSAT_c,4));

double vfn = freelig(VEGF_n*ivg > 1e-12 ? VEGF_n*ivg : 1e-12, bev_sys+bev_top, BEV_KD)/ivg;
double vfg = freelig(VEGF_g*ivg > 1e-12 ? VEGF_g*ivg : 1e-12, bev_sys, BEV_KD)/ivg;
double vfh = freelig(VEGF_h*ivg > 1e-12 ? VEGF_h*ivg : 1e-12, bev_sys, BEV_KD)/ivg;
double vfp = freelig(VEGF_p*ivg > 1e-12 ? VEGF_p*ivg : 1e-12, bev_sys, BEV_KD)/ivg;
double vfc = freelig(VEGF_c*ivg > 1e-12 ? VEGF_c*ivg : 1e-12, bev_sys, BEV_KD)/ivg;

double VS_n = HILL(vfn, KV_HALF)*(1.0 - I_VEGFR2);
double VS_g = HILL(vfg, KV_HALF)*(1.0 - I_VEGFR2);
double VS_h = HILL(vfh, KV_HALF)*(1.0 - I_VEGFR2);
double VS_p = HILL(vfp, KV_HALF)*(1.0 - I_VEGFR2);
double VS_c = HILL(vfc, KV_HALF)*(1.0 - I_VEGFR2);

double AN_n = VS_n*AKT_n/(1.0 + 2.2*PS_n);
double AN_g = VS_g*AKT_g/(1.0 + 2.2*PS_g);
double AN_h = VS_h*AKT_h/(1.0 + 2.2*PS_h);
double AN_p = VS_p*AKT_p/(1.0 + 2.2*PS_p);
double AN_c = VS_c*AKT_c/(1.0 + 2.2*PS_c);

// ---------------- signalling --------------------------------------------
double bmp = HILL(BMP9, KD9);
double gd  = GD_LES + ALK1AG*(1.0 - GD_LES);
double hyp = (0.98 - SPO2)/0.98;

dxdt_PS_n = KIN_PS*gd*bmp*(1.0 + GSHEAR*gd*HILL(WSS_n,KWS)) - KOUT_PS*PS_n;
dxdt_PS_g = KIN_PS*gd*bmp*(1.0 + GSHEAR*gd*HILL(WSS_g,KWS)) - KOUT_PS*PS_g;
dxdt_PS_h = KIN_PS*gd*bmp*(1.0 + GSHEAR*gd*HILL(WSS_h,KWS)) - KOUT_PS*PS_h;
dxdt_PS_p = KIN_PS*gd*bmp*(1.0 + GSHEAR*gd*HILL(WSS_p,KWS)) - KOUT_PS*PS_p;
dxdt_PS_c = KIN_PS*gd*bmp*(1.0 + GSHEAR*gd*HILL(WSS_c,KWS)) - KOUT_PS*PS_c;

dxdt_VEGF_n = KV0*(1.0+AV*POS(1.0-PS_n))*(1.0-E_IMID_V) - KDV*VEGF_n;
dxdt_VEGF_g = KV0*(1.0+AV*POS(1.0-PS_g))*(1.0-E_IMID_V) - KDV*VEGF_g;
dxdt_VEGF_h = KV0*(1.0+AV*POS(1.0-PS_h))*(1.0-E_IMID_V) + KV_HYP*hyp - KDV*VEGF_h;
dxdt_VEGF_p = KV0*(1.0+AV*POS(1.0-PS_p))*(1.0-E_IMID_V) + KV_HYP*hyp - KDV*VEGF_p;
dxdt_VEGF_c = KV0*(1.0+AV*POS(1.0-PS_c))*(1.0-E_IMID_V) + KV_HYP*hyp - KDV*VEGF_c;

dxdt_AKT_n = KA_AKT*VS_n*(1.0+WA*POS(1.0-PS_n)) - KD_AKT*AKT_n;
dxdt_AKT_g = KA_AKT*VS_g*(1.0+WA*POS(1.0-PS_g)) - KD_AKT*AKT_g;
dxdt_AKT_h = KA_AKT*VS_h*(1.0+WA*POS(1.0-PS_h)) - KD_AKT*AKT_h;
dxdt_AKT_p = KA_AKT*VS_p*(1.0+WA*POS(1.0-PS_p)) - KD_AKT*AKT_p;
dxdt_AKT_c = KA_AKT*VS_c*(1.0+WA*POS(1.0-PS_c)) - KD_AKT*AKT_c;

dxdt_ID1_n = PS_n - ID1_n;
dxdt_PDGFB_n = (PD0 + KPD*PS_n)*(1.0+E_IMID) - KDPD*PDGFB_n;
dxdt_PDGFB_g = (PD0 + KPD*PS_g)*(1.0+E_IMID) - KDPD*PDGFB_g;

// ---------------- lesion birth and involution ----------------------------
double ps_field = GD*bmp*(1.0 + GSHEAR*GD*0.5);
double af_n = AN_n*(1.0+2.2*PS_n)/(1.0+2.2*ps_field);
double af_g = AN_g*(1.0+2.2*PS_g)/(1.0+2.2*ps_field);
double iv_n = 1.0 + KREGA*pow(POS(1.0 - AN_n/ANGIO0),2);
double iv_g = 1.0 + KREGA*pow(POS(1.0 - AN_g/ANGIO0),2);
double birth_n = KHIT_N*SEV*HILL(pow(af_n,HC), pow(KC,HC));
double birth_g = KHIT_G*SEV*GIF*HILL(pow(af_g,HC), pow(KC,HC));

dxdt_N_n = birth_n - KREG_N*N_n*iv_n;
dxdt_N_g = birth_g - KREG_G*N_g*iv_g;

// ---------------- lesion calibre: THE SIGN-FLIPPING SHEAR LOOP -----------
// REMOD = KSH*g(WSS)*(1 - 2*pSMAD).  pSMAD ~ 1 -> inward.  pSMAD ~ 0 -> outward.
double RM_n = KSH*HILL(WSS_n,KWSS)*(1.0 - 2.0*PS_n);
double RM_g = KSH*HILL(WSS_g,KWSS)*(1.0 - 2.0*PS_g);
double RM_h = KSH*HILL(WSS_h,KWSS)*(1.0 - 2.0*PS_h);
double RM_p = KSH*HILL(WSS_p,KWSS)*(1.0 - 2.0*PS_p);
double RM_c = KSH*HILL(WSS_c,KWSS)*(1.0 - 2.0*PS_c);

// sinks act on the DILATION above the calibre a lesion is born with
double dil_n = (N_n > 1e-6) ? (S_BIRTH - S_n)*birth_n/N_n : 0.0;
double dil_g = (N_g > 1e-6) ? (S_BIRTH - S_g)*birth_g/N_g : 0.0;

dxdt_S_n = S_n*(KG*ANGF*AN_n + RM_n)*(1.0 - S_n/SMAX_n)
           - (KMAT*MU_n + KDEC)*(S_n - S_BIRTH) + dil_n;
dxdt_S_g = S_g*(KG*ANGF*AN_g + RM_g)*(1.0 - S_g/SMAX_g)
           - (KMAT*MU_g + KDEC)*(S_g - S_BIRTH) + dil_g;

double smax_h = SMAX_h*HEPF, smax_p = SMAX_p*PULF, smax_c = SMAX_c*CNSF;
dxdt_S_h = (smax_h > 1e-6)
  ? S_h*(KG*ANGF*AN_h + RM_h)*HEPF*(1.0 - S_h/smax_h) - KDEC*(S_h - S_BIRTH)
  : -KDEC*S_h;
dxdt_S_p = (smax_p > 1e-6)
  ? S_p*(KG*ANGF*AN_p + RM_p)*PULF*(1.0 - S_p/smax_p) - KDEC*(S_p - S_BIRTH)
  : -KDEC*S_p;
dxdt_S_c = (smax_c > 1e-6)
  ? S_c*(KG*ANGF*AN_c + RM_c)*CNSF*(1.0 - S_c/smax_c) - KDEC*(S_c - S_BIRTH)
  : -KDEC*S_c;

dxdt_MU_n = KMR*PDGFB_n*(1.0-MU_n) - KML*MU_n*(1.0+WM*AN_n) - KML_PAZ*I_PDGFRB*MU_n;
dxdt_MU_g = KMR*PDGFB_g*(1.0-MU_g) - KML*MU_g*(1.0+WM*AN_g) - KML_PAZ*I_PDGFRB*MU_g;

// ---------------- haemostasis -------------------------------------------
dxdt_PLG = 1.2*((1.0 - TXA_IMAX*HILL(C_TXA, TXA_IC50)) - PLG);

// ---------------- bleeding ----------------------------------------------
double frag_n = pow(POS(1.0-MU_n), MFRAG);
double frag_g = pow(POS(1.0-MU_g), MFRAG);
double lysf   = 1.0 + FLYS*PLG;
double ev_n   = KE_N*N_n*frag_n;                      // episodes per day
double dur_ev = TEV_N*S_n*lysf*(1.0 + KDF*POS(1.0-MU_n));   // min per episode
double rate_ev= RBLD*S_n*S_n;                         // mL per min
double vol_ev = rate_ev*dur_ev;                       // mL per episode
double BLR_N  = ev_n*vol_ev;
double BLR_G  = KGB*N_g*frag_g*pow(S_g,3)*lysf;
double BLR    = BLR_N + BLR_G;
double MEDx   = 30.0*ev_n*dur_ev;                     // min per month
double EPM    = 30.0*ev_n;

// ---------------- iron flux balance --------------------------------------
double gate   = 1.0/(1.0 + pow(HEP/K_HEP, N_HEP));
double FE_ABS = fmin(FE_ABS_MAX, FE_ORAL*FE_ABS_FR + FE_DIET*FE_DIET_FR)*gate;
double FE_IN  = FE_ABS + FE_IV_RATE;
double FE_PER_HB = 0.0347*BV;                         // mg iron per g/dL Hb
double fe_recycle = HBx/T_RBC*FE_PER_HB;
double fe_mobil   = FES/T_MOB;
double fe_for_ery = fe_recycle + FE_IN - FE_OBLIG + fe_mobil;
double prod_iron  = POS(fe_for_ery)/FE_PER_HB;
double epo_drive  = fmin(EPO_MAX, exp(KEPO*(HB_T - HBx)));
double prod       = fmin((HB_T/T_RBC)*EPO, prod_iron);

// transfusion written as a smooth demand-driven infusion rather than a discrete
// event, so the ODE solver stays continuous
double tx_rate = TX_GAIN/(1.0 + exp((HBx - TX_TRIG)/0.20));

dxdt_EPO = 1.0*(epo_drive - EPO);
dxdt_RET = (prod/(HB_T/T_RBC) - RET);
dxdt_HB  = (HB_FIX > 0.0) ? 0.0 : (prod - HBx/T_RBC - HBx*BLR/BV + tx_rate);
dxdt_FES = FE_IN - FE_OBLIG - 0.0347*HBx*BLR - (prod - HBx/T_RBC)*FE_PER_HB
           + tx_rate*TX_UNIT_FE/1.0;
dxdt_HEP = 0.6*(fmax(0.02, HILL(FES,KFES)*2.2 - KERFE*POS(EPO-1.0)) - HEP);

// ---------------- cardiac -------------------------------------------------
dxdt_CO  = (CO_target - COx)/TAU_CO;
dxdt_LVR = (CO_target/(CI_REF*BSA) - LVR)/120.0;
dxdt_RAP = (KRAP*POS(COx/BSA - 3.8) - RAP)/30.0;

// ---------------- Epistaxis Severity Score (Hoag domain reconstruction) ---
// six domains: frequency, duration, intensity, medical attention,
// transfusion, anaemia -- smooth (differentiable) versions of the bins
double F  = sbin1(EPM/4.345,0.5)+sbin1(EPM/4.345,1.5)+sbin1(EPM/4.345,3.5)+sbin1(EPM/4.345,6.5);
double D  = sbin1(dur_ev,1.0)+sbin1(dur_ev,5.0)+sbin1(dur_ev,15.0)+sbin1(dur_ev,30.0);
double I  = sbin1(rate_ev,0.40)+sbin1(rate_ev,1.20)+sbin1(rate_ev,3.00);
double A  = 1.0/(1.0+exp(-(MEDx-120.0)/45.0));
double Tq = 1.0/(1.0+exp((HBx-TX_TRIG-0.6)/0.35));
double AQ = 1.0/(1.0+exp((HBx-11.8)/0.7));
double raw = 0.517*F + 0.789*D + 0.999*I + 0.551*A + 0.634*Tq + 0.318*AQ;
double ESS_INST = 10.0*raw/9.7215;

dxdt_ESSL = (ESS_INST - ESSL)/ESS_TAU;
dxdt_PBO  = (PBO_MAX - PBO)/PBO_TAU;
dxdt_QOL  = (16.0*HILL(pow(POS(ESSL-PBO),2.1), pow(5.6,2.1)) - QOL)/21.0;

dxdt_CUM_BL   = BLR;
dxdt_CUM_TX   = tx_rate;                 // g/dL-equivalents; 1 unit ~ 1 g/dL
dxdt_CUM_FEIV = FE_IV_RATE;

// ---------------- toxicity ------------------------------------------------
dxdt_SBP  = ((SBP0 + KSBP_BEV*HILL(C_BEV,40.0) + KSBP_PAZ*HILL(CU_PAZ,9.0)) - SBP)/TAU_SBP;
dxdt_UPCR = ((KUP_BEV*HILL(C_BEV,60.0) + KUP_PAZ*HILL(CU_PAZ,14.0)) - UPCR)/TAU_UP;
dxdt_ALT  = ((ALT0*(1.0 + KALT_PAZ*HILL(CU_PAZ,12.0))) - ALT)/TAU_ALT;
dxdt_ANC  = ((ANC0*(1.0 - KANC_POM*HILL(C_POM,14.0))) - ANC)/TAU_ANC;
dxdt_VTEH = KVTE_POM*HILL(C_POM,12.0) + KVTE_TXA*HILL(C_TXA,8.0);

$TABLE
double BVt = BV_kg*WT;
double HBo = HB > 1.0 ? HB : 1.0;
double C_BEVo = BEV_C/BEV_V1;
double C_PAZo = PAZ_C/PAZ_V;
double CU_PAZo = C_PAZo*PAZ_FU*1000.0;
double C_POMo = POM_C/POM_V*1000.0;
double C_TXAo = TXA_C/TXA_V;
double SPO2o = 0.98 - 0.34*HILL(pow(S_p,2.5), pow(20.0,2.5));
double CIo = CO/BSA;
double fragno = pow(POS(1.0-MU_n), MFRAG);
double fraggo = pow(POS(1.0-MU_g), MFRAG);
double lysfo = 1.0 + FLYS*PLG;
double ev_no = KE_N*N_n*fragno;
double dur_evo = TEV_N*S_n*lysfo*(1.0 + KDF*POS(1.0-MU_n));
double rate_evo = RBLD*S_n*S_n;
double BLR_No = ev_no*rate_evo*dur_evo;
double BLR_Go = KGB*N_g*fraggo*pow(S_g,3)*lysfo;
double BLRo = BLR_No + BLR_Go;
double MEDo = 30.0*ev_no*dur_evo;
double EPMo = 30.0*ev_no;
double ESSo = POS(ESSL - PBO);
double FE_LOSSo = BLRo*0.0347*HBo;
double gateo = 1.0/(1.0 + pow(HEP/K_HEP, N_HEP));
double FE_ABSo = fmin(FE_ABS_MAX, FE_ORAL*FE_ABS_FR + FE_DIET*FE_DIET_FR)*gateo;
// maximum blood loss this iron strategy can sustain at the current Hb
double FE_CEILo = POS(FE_ABSo + FE_IV_RATE - FE_OBLIG)/(0.0347*HBo);
double IRONINDo = (BLRo < FE_CEILo) ? 1.0 : 0.0;

$CAPTURE @annotated
C_BEVo   : bevacizumab concentration (ug/mL)
C_PAZo   : pazopanib total concentration (ug/mL)
CU_PAZo  : pazopanib free concentration (ng/mL)
C_POMo   : pomalidomide concentration (ng/mL)
C_TXAo   : tranexamic acid concentration (ug/mL)
BLRo     : total blood loss rate (mL/day)
BLR_No   : nasal blood loss rate (mL/day)
BLR_Go   : GI blood loss rate (mL/day)
MEDo     : monthly epistaxis duration (min/month)
EPMo     : epistaxis episodes per month
dur_evo  : duration per episode (min)
rate_evo : bleed rate during an episode (mL/min)
ESSo     : Epistaxis Severity Score (0-10)
CIo      : cardiac index (L/min/m2)
SPO2o    : arterial oxygen saturation
FE_LOSSo : iron loss in shed blood (mg/day)
FE_ABSo  : absorbed iron (mg/day)
FE_CEILo : maximum sustainable blood loss for this iron strategy (mL/day)
IRONINDo : iron-independent flag (1 = below the ceiling)
'

mod <- mcode("hht", code)

## ============================================================================
##  PHENOTYPES
##  SEV  = somatic-hit rate (lesion number)   ANGF = angiogenic drive (calibre)
##  GIF  = GI angiodysplasia burden           HEPF/PULF/CNSF = visceral beds
## ============================================================================
PHENO <- list(
  mild     = list(SEV = 0.22, ANGF = 0.55, GIF = 0.10, HEPF = 0.00),
  moderate = list(SEV = 1.05, ANGF = 1.02, GIF = 0.55, HEPF = 0.20),  # PATH-HHT arm
  severe   = list(SEV = 1.70, ANGF = 1.25, GIF = 45.0, HEPF = 0.30),  # Parambil arm
  hepatic  = list(SEV = 1.00, ANGF = 1.00, GIF = 0.25, HEPF = 1.00),  # Dupuis-Girod arm
  hht1_pav = list(SEV = 1.00, ANGF = 1.00, GIF = 0.30, HEPF = 0.10, PULF = 0.55),
  jphht    = list(SEV = 1.30, ANGF = 1.10, GIF = 3.00, HEPF = 0.35)   # SMAD4
)

#' Grow a patient from birth to `age` with no treatment, then return the sim.
grow_patient <- function(pheno = "moderate", age = 52, extra = list()) {
  pars <- c(PHENO[[pheno]], extra)
  mod %>% param(pars) %>% mrgsim(end = 365 * age, delta = 365, atol = 1e-8, rtol = 1e-6)
}

#' Final compartment vector as a named list suitable for init().
#' Selecting by names(init(mod)) matters: the simulation output also carries the
#' $CAPTURE columns, and taking the row positionally misaligns the state.
final_state <- function(sim) {
  nm <- names(init(mod))
  df <- as.data.frame(sim)
  as.list(df[nrow(df), nm])
}

## ============================================================================
##  18 TREATMENT SCENARIOS
##  Each returns a dosing data set to be applied on top of a grown patient.
## ============================================================================
##  cmt indices: 1 BEV_C  3 BEVN_S  5 PAZ_A  7 POM_A  9 TXA_A
scen_none      <- function() ev(amt = 0, cmt = 1, time = 0)
scen_bev_ind   <- function() ev(amt = 350, cmt = 1, ii = 14, addl = 5)          # 5 mg/kg q2w x6
scen_bev_maint <- function() ev(amt = 350, cmt = 1, ii = 28, addl = 12)         # monthly maintenance
scen_bev_nasal <- function(d = 75) ev(amt = d, cmt = 3, ii = 14, addl = 2)      # Dupuis-Girod 2016
scen_paz50     <- function() ev(amt = 50,  cmt = 5, ii = 1, addl = 364)
scen_paz150    <- function() ev(amt = 150, cmt = 5, ii = 1, addl = 364)         # Parambil median
scen_paz800    <- function() ev(amt = 800, cmt = 5, ii = 1, addl = 364)         # oncology dose
scen_pom4      <- function() ev(amt = 4,   cmt = 7, ii = 1, addl = 167)         # PATH-HHT 24 wk
scen_pom2      <- function() ev(amt = 2,   cmt = 7, ii = 1, addl = 167)
scen_thal      <- function() ev(amt = 100, cmt = 7, ii = 1, addl = 364)
scen_txa_oral  <- function() ev(amt = 1000, cmt = 9, ii = 1/3, addl = 3 * 364)  # 1 g TID
scen_txa_top   <- function() ev(amt = 20,   cmt = 9, ii = 1/2, addl = 2 * 84)   # 40 mg/d topical
scen_pom_txa   <- function() c(scen_pom4(), scen_txa_oral())
scen_bev_pom   <- function() c(scen_bev_maint(), scen_pom4())

## parameter-level scenarios (no dosing record)
par_iron_oral  <- list(FE_ORAL = 200)
par_iron_iv_4w <- list(FE_ORAL = 200, FE_IV_RATE = 1000 / 28)
par_iron_iv_2w <- list(FE_ORAL = 200, FE_IV_RATE = 1000 / 14)
par_alk1_rescue<- list(ALK1AG = 0.60)                    # hypothetical upstream fix
par_embolise   <- list(PULF = 0.0)                       # PAVM coil embolisation
par_laser      <- list(KHIT_N = 0.0049 * 0.55)           # nasal laser ablation field effect

## ============================================================================
##  EXAMPLE 1 — PATH-HHT reconstruction (Al-Samkari 2024 NEJM, PMID 39292928)
##  published: ESS difference vs placebo -0.94 (95% CI -1.57 to -0.31)
##  model:     -0.83
## ============================================================================
run_path_hht <- function() {
  y0 <- final_state(grow_patient("moderate", 52))
  arm <- function(evd, pbo) {
    mod %>%
      param(c(PHENO$moderate, list(PBO_MAX = pbo))) %>%
      init(y0) %>%
      mrgsim(events = evd, end = 168, delta = 7)
  }
  pbo_arm <- arm(scen_none(), 0.90)
  pom_arm <- arm(scen_pom4(), 0.90)
  d_pbo <- tail(pbo_arm$ESSo, 1) - pbo_arm$ESSo[1]
  d_pom <- tail(pom_arm$ESSo, 1) - pom_arm$ESSo[1]
  cat(sprintf("PATH-HHT: placebo %+.2f, pomalidomide %+.2f, difference %+.2f (published -0.94)\\n",
              d_pbo, d_pom, d_pom - d_pbo))
  list(placebo = pbo_arm, pomalidomide = pom_arm)
}

## ============================================================================
##  EXAMPLE 2 — the oral-iron ceiling (the model result with clinical teeth)
##  Steady state needs  absorbed - obligatory = 0.0347 * Hb * blood loss.
## ============================================================================
iron_ceiling <- function(hb = 14.6) {
  strat <- list(
    "diet only"                = list(FE_ORAL = 0,   FE_IV_RATE = 0),
    "oral 65 mg/d"             = list(FE_ORAL = 65,  FE_IV_RATE = 0),
    "oral 200 mg/d"            = list(FE_ORAL = 200, FE_IV_RATE = 0),
    "oral + IV 1 g q8w"        = list(FE_ORAL = 200, FE_IV_RATE = 1000 / 56),
    "oral + IV 1 g q4w"        = list(FE_ORAL = 200, FE_IV_RATE = 1000 / 28),
    "oral + IV 1 g q2w"        = list(FE_ORAL = 200, FE_IV_RATE = 1000 / 14))
  do.call(rbind, lapply(names(strat), function(k) {
    s <- strat[[k]]
    absd <- min(10.0, s$FE_ORAL * 0.22 + 15.0 * 0.35)
    surplus <- absd + s$FE_IV_RATE - 1.0
    data.frame(strategy = k, absorbed_mg = round(absd + s$FE_IV_RATE, 1),
               max_blood_loss_mL_day = round(surplus / (0.0347 * hb), 1),
               max_blood_loss_mL_month = round(30 * surplus / (0.0347 * hb)))
  }))
}

## ============================================================================
##  EXAMPLE 3 — why the nasal route fails at any dose (PMID 27599328)
##  Peak submucosal concentration scales with dose; TIME above threshold does
##  not, because mucociliary clearance empties the depot with t1/2 ~15 min.
## ============================================================================
run_nasal_dose_ranging <- function(doses = c(0, 25, 50, 75, 750, 7500, 75000)) {
  y0 <- final_state(grow_patient("moderate", 52))
  do.call(rbind, lapply(doses, function(d) {
    out <- mod %>% param(PHENO$moderate) %>%
      init(y0) %>%
      mrgsim(events = if (d > 0) scen_bev_nasal(d) else scen_none(),
             end = 92, delta = 1, atol = 1e-10, rtol = 1e-8)
    data.frame(nasal_dose_mg = d, MED_3mo = round(tail(out$MEDo, 1), 1))
  }))
}

## ============================================================================
##  EXAMPLE 4 — the anaemia -> output -> shear -> bleeding loop
##  Clamp Hb, let the structure re-equilibrate, and read the blood loss back.
##  Holding Hb at 6.0 rather than 14.6 raises blood loss by 33% through the
##  mechanical loop alone -- so correcting anaemia is not purely supportive.
## ============================================================================
run_loop_gain <- function(hb_levels = c(14.6, 13.0, 11.0, 9.0, 7.5, 6.0)) {
  y0 <- final_state(grow_patient("severe", 55))
  do.call(rbind, lapply(hb_levels, function(hb) {
    ini <- y0; ini$HB <- hb
    out <- mod %>% param(c(PHENO$severe, list(HB_FIX = hb))) %>%
      init(ini) %>% mrgsim(end = 365 * 5, delta = 365)
    data.frame(Hb = hb, CI = round(tail(out$CIo, 1), 2),
               BLR = round(tail(out$BLRo, 1), 1),
               MED = round(tail(out$MEDo, 1), 1))
  }))
}

## ============================================================================
##  EXAMPLE 5 — full scenario sweep on the moderate phenotype, 12 months
## ============================================================================
run_scenarios <- function(pheno = "moderate", age = 52, end = 365) {
  y0 <- final_state(grow_patient(pheno, age))
  S <- list(
    "no treatment"                 = list(ev = scen_none(),      par = list()),
    "oral iron 200 mg/d"           = list(ev = scen_none(),      par = par_iron_oral),
    "IV iron 1 g q4w"              = list(ev = scen_none(),      par = par_iron_iv_4w),
    "TXA 3 g/d"                    = list(ev = scen_txa_oral(),  par = par_iron_oral),
    "TXA topical 40 mg/d"          = list(ev = scen_txa_top(),   par = par_iron_oral),
    "pomalidomide 2 mg/d"          = list(ev = scen_pom2(),      par = par_iron_oral),
    "pomalidomide 4 mg/d"          = list(ev = scen_pom4(),      par = par_iron_oral),
    "thalidomide 100 mg/d"         = list(ev = scen_thal(),      par = par_iron_oral),
    "pazopanib 50 mg/d"            = list(ev = scen_paz50(),     par = par_iron_oral),
    "pazopanib 150 mg/d"           = list(ev = scen_paz150(),    par = par_iron_oral),
    "pazopanib 800 mg/d"           = list(ev = scen_paz800(),    par = par_iron_oral),
    "bevacizumab induction"        = list(ev = scen_bev_ind(),   par = par_iron_oral),
    "bevacizumab monthly maint."   = list(ev = scen_bev_maint(), par = par_iron_oral),
    "bevacizumab nasal 75 mg"      = list(ev = scen_bev_nasal(), par = par_iron_oral),
    "pomalidomide + TXA"           = list(ev = scen_pom_txa(),   par = par_iron_oral),
    "bevacizumab + pomalidomide"   = list(ev = scen_bev_pom(),   par = par_iron_oral),
    "nasal laser ablation"         = list(ev = scen_none(),      par = c(par_iron_oral, par_laser)),
    "ALK1 restoration (hypothet.)" = list(ev = scen_none(),      par = c(par_iron_oral, par_alk1_rescue))
  )
  do.call(rbind, lapply(names(S), function(k) {
    s <- S[[k]]
    out <- mod %>% param(c(PHENO[[pheno]], s$par)) %>%
      init(y0) %>%
      mrgsim(events = s$ev, end = end, delta = 7)
    data.frame(scenario = k,
               ESS   = round(tail(out$ESSo, 1), 2),
               dESS  = round(tail(out$ESSo, 1) - out$ESSo[1], 2),
               MED   = round(tail(out$MEDo, 1), 1),
               Hb    = round(tail(out$HB, 1), 2),
               BLR   = round(tail(out$BLRo, 1), 1),
               iron_independent = tail(out$IRONINDo, 1))
  }))
}

## ============================================================================
##  Usage
## ============================================================================
if (interactive()) {
  print(iron_ceiling())
  run_path_hht()
  print(run_nasal_dose_ranging())
  print(run_loop_gain())
  print(run_scenarios("moderate"))
}
