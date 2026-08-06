## =====================================================================
##  Hereditary Spherocytosis (HS) — mrgsolve QSP model
##  유전성 구상적혈구증 — 72-ODE 정량적 시스템 약리학 모델
##
##  This file mirrors hsph_python_reference.py equation for equation.  The
##  Python file is the model of record: it is the one that was executed,
##  calibrated and used to produce hsph_reference_output.txt.  No R
##  toolchain was available in the build environment, so THIS FILE HAS NOT
##  BEEN RUN.  Treat any discrepancy with the Python reference as a bug in
##  this file.
##
##  ------------------------------------------------------------------
##  THE IDEA
##  A red cell is two numbers: membrane area S (um^2) and volume V (um^3).
##  Their combination gives the minimum cylindrical diameter
##
##      V_sph = S^1.5 / (6 sqrt(pi))          (critical haemolytic volume)
##      s     = V / V_sph                     (sphericity, 0..1)
##      D_c   = 2 sqrt(S/pi) cos(acos(-s)/3 - 2pi/3)
##
##  and D_c is the ONLY property of a red cell that the splenic sinus wall
##  can measure.  The disease is the loop
##
##      area loss  ->  D_c up  ->  longer cordal residence  ->  area loss
##
##  closed inside one organ, which is why splenectomy works and why it
##  works differently depending on which membrane protein is missing.
##
##  Usage
##    library(mrgsolve); library(dplyr)
##    mod <- mread("hsph_mrgsolve_model", ".")
##    mod %>% mrgsim(end = 1400, delta = 1) %>% plot(Hb + RETpct + TBIL ~ time)
##
##    ## moderate HS, then splenectomy at day 400
##    mod %>% param(fdef = 0.30) %>%
##      ev(ev(time = 400, amt = 0, evid = 8, cmt = 1)) %>%   # see $PREAMBLE
##      mrgsim(end = 1400)
## =====================================================================

$PROB
# Hereditary spherocytosis: a geometry-driven haemolysis model
# 72 ODEs — 9 red cell age cohorts x 5 states, plus erythropoiesis,
# spleen, haem catabolism, iron and drug PK/PD.

$PARAM @annotated
// ---------------- whole body
BV       :    5   : blood volume (L)
// ---------------- genotype
fdef     :     0  : vertical-linkage deficit, 1 - spectrin fraction (-)
a_ent_def:     1  : fraction of the deficit already present at release (-)
f_b3ves  :     1  : band 3 content of the shed vesicle vs parent membrane (-)
f_hbves  :  0.05  : haemoglobin content of the shed vesicle vs parent (-)
// ---------------- reticulocyte entry geometry
A0       :   140  : membrane area at release (um2)
V0       :    94  : cell volume at release (um3)
H0       :    30  : haemoglobin per cell at release (pg)
B0       :    1.2 : band 3 copies at release (1e6/cell)
V_RET_MULT: 1.22  : reticulocyte volume excess (-)
// ---------------- membrane loss  [kv_base CALIBRATED: Waugh area loss]
kv_base  : 0.000923326 : basal fractional area loss (1/day)
kv_def   : 0.00307772 : extra fractional area loss per unit fdef (1/day)
cordamp  :    1   : area loss inside the cord, multiple of circulating (-)
p_atp    :  1.5   : ATP exponent on vesiculation (-)
r_ves3   :  0.025 : vesicle volume/area ratio = radius/3 (um)
// ---------------- dehydration  [kd_base CALIBRATED: Waugh volume loss]
kd_base  : 0.00366211 : basal dehydration rate (1/day)
kd_def   : 0.0244142 : constitutive cation leak per unit fdef; NOT fitted, = 6.667 x kd_base (1/day)
kd_cord  : 0.150146 : extra dehydration at full cordal residence; NOT fitted, = 41 x kd_base (1/day)
MCHC_MAX :    41   : maximum MCHC, sets the volume floor (g/dL)
// ---------------- splenic geometric filter
spl_flow_exp: 0.35 : splenic flow scales as mass^0.35 not mass^1; with linear scaling the hazard-workload-mass-flow loop gain exceeds 1 and the steady state stops being unique (-)
spl_flow_cap: 1.6  : splenic flow never exceeds 1.6x normal (-)
f_pass0  :    72   : splenic passes per cell per day (1/day)
D50      :    3.27 : D_c at which half the passes go cordal (um)
wD       :    0.09 : steepness of the cordal filter (um)
tau0     :   0.00104 : cordal dwell of a normal cell (day)
Dc_ref   :   2.856 : reference D_c for the dwell law (um)
w_esc    :     0.3 : e-folding of dwell time with D_c (um)
visc_k   :   0.231 : dehydration retention coefficient  [CALIBRATED]; the published bulk cytoplasmic-viscosity slope is ln2/3 = 0.231 per g/dL of MCHC (1/(g/dL))
MCHC_ref :    33   : reference MCHC for the viscosity law (g/dL)
R_MAX    :    0.6  : cap on the cordal time fraction (-)
k_ph     :  0.0302021 : phagocytic hazard while held in the cord (1/day)
// ---------------- opsonic arm (band 3 clustering -> natural IgG)
Bden_ref : 0.00857143 : band 3 copies per um2, normal (1e6/um2)
IgG_bg   :    45   : background red-cell-bound IgG (molecules/cell)
IgG_span :   180   : maximum excess IgG (molecules/cell)
K_cl     : 0.17589  : band 3 density excess at half-maximal IgG (-)
IgG_thr  :    70   : IgG below which nothing is eaten (molecules/cell)
IgG50    :    38   : excess IgG for half-maximal phagocytosis (molecules/cell)
m_ops    :     3   : Hill coefficient, opsonic phagocytosis (-)
k_ops    : 0.00906356  : maximal opsonic clearance rate (1/day)
spl_ops_cap: 1.5   : splenic opsonic capacity scales with mass only to 1.5x (-)
w_spl_ops:   0.72  : splenic share of opsonic clearance (-)
w_liv_ops:   0.28  : hepatic share of opsonic clearance (-)
// ---------------- other clearance
k_liv_geo:   0.035 : hepatic geometric clearance rate (1/day)
D_liv    :   3.45  : D_c at which hepatic retention is half-maximal (um)
w_liv    :   0.16  : steepness of hepatic retention (um)
k_lys    :   0.06  : intravascular lysis rate at s = 1 (1/day)
s_lys    :    0.97 : sphericity above which lysis begins (-)
k_sen    :    0.8  : senescence hazard scale (1/day)
tau50    : 211.238  : senescence half-time  [CALIBRATED: lifespan 120 d] (day)
m_sen    :     8   : senescence Hill exponent (-)
// ---------------- erythropoiesis
k_epo    :    24   : EPO turnover (1/day)
EPO_a    :   4.17  : intercept of log10(EPO) vs Hb (-)
EPO_b    :   0.211 : slope of log10(EPO) vs Hb (dL/g)
EPO_norm :    10   : normal EPO (mIU/mL)
Emax_mar :     8   : maximum erythropoietic expansion (-)
K_epo    :    48   : EPO excess at half-maximal expansion (mIU/mL)
k_prog   : 0.333333: progenitor transit (1/day)
k_erb    :   0.25  : erythroblast transit (1/day)
k_retm   :    0.4  : marrow reticulocyte release (1/day)
tau_ret0 :     1   : blood reticulocyte maturation time (day)
prod0    :   0.0417: baseline red cell production (1e12/L/day)
shift_max:   1.6   : maximal marrow release acceleration (-)
// ---------------- spleen
SPL_base :   150   : normal spleen volume (mL)
SPL_max  :   1200  : maximum spleen volume (mL)
k_spl    : 0.00454545 : spleen remodelling rate (1/day)
g_spl    :  0.146  : spleen volume gain per unit workload (-)
spl_frac :     1   : surviving splenic fraction, 0 = splenectomised (-)
k_regrow :     0   : remnant regrowth rate after partial splenectomy (1/day)
k_mac    : 0.0714286: macrophage activation turnover (1/day)
a_mac    :   0.8   : macrophage activation gain (-)
b_mac    :     2   : macrophage activation saturation (-)
f_sen_spl:    0.3  : splenic share of senescent-cell clearance (-)
f_lys_spl:     0   : splenic share of intravascular lysis (-)
W0       :   0.0125: reference splenic erythrophagocytic workload (1e12/L/d)
// ---------------- haem catabolism
BR_PER_G :    34   : bilirubin produced per g of haemoglobin (mg/g)
Vd_bil   :     5   : bilirubin distribution volume (L)
Vmax_ugt :  12750  : UGT1A1 maximal conjugation (mg/day)
Km_ugt   :    30   : UGT1A1 Km (mg/dL)
ugt_f    :     1   : UGT1A1 activity factor; 0.287 = UGT1A1*28 homozygote (-)
k_bilc   :    20   : conjugated bilirubin elimination (1/day)
k_bile   :     1   : biliary bilirubin turnover (1/day)
f_ucb_bile:  0.012 : unconjugated fraction of biliary bilirubin (-)
k_stone  :  2.3e-06: stone nucleation rate constant (1/mg/day)
k_hp_syn :    14   : haptoglobin synthesis (mg/dL/day)
k_hp_deg :   0.14  : haptoglobin turnover (1/day)
k_hp_bind:     9   : haptoglobin consumption by free Hb (1/day)
k_fhb    :     8   : free haemoglobin clearance (1/day)
LDH_per  :     1   : LDH scaling placeholder (-)
k_ldh    :   2.2   : LDH turnover (1/day)
LDH0     :   170   : baseline LDH (U/L)
// ---------------- iron
k_hepc   :   0.5   : hepcidin turnover (1/day)
HEPC0    :     1   : normal hepcidin (relative)
erfe_k   :   0.45  : erythroferrone suppression of hepcidin (-)
k_ferr   : 0.0333333: ferritin turnover (1/day)
fe_per_cell: 1.05  : iron per 1e12 red cells (mg)
k_fe_abs :   0.0012: duodenal iron absorption scale (mg/day)
k_fe_liv :   0.0025: liver iron turnover (1/day)
tx_fe    :   200   : iron per transfused unit (mg)
k_chel   :     0   : iron chelation rate; deferasirox ~0.0016 (1/day)
// ---------------- mitapivat PK/PD
ka_m     :   1.6   : absorption rate (1/h)
CL_m     :   4.2   : clearance (L/h)
Vc_m     :    38   : central volume (L)
Q_m      :     3   : intercompartmental clearance (L/h)
Vp_m     :    25   : peripheral volume (L)
EC50_m   :  1150   : concentration for half-maximal PKR activation (ng/mL)
Emax_atp :   0.62  : maximal fractional ATP increase (-)
Emax_dpg :   0.45  : maximal fractional 2,3-DPG decrease (-)
k_atp    :     2   : ATP turnover (1/day)
k_dpg    :   1.2   : 2,3-DPG turnover (1/day)
dose_m   :     0   : mitapivat dose per administration, BID (mg)
// ---------------- insults and cofactors
parvo_t  :    -1   : parvovirus B19 onset; <0 = never (day)
parvo_dur:     8   : duration of erythroid arrest (day)
parvo_supp:  0.97  : fractional suppression of progenitor input (-)
fol_ok   :     1   : folate supply, 1 = replete (-)
k_fol    :   0.05  : folate turnover (1/day)
// ---------------- transfusion
tx_start :    -1   : first transfusion; <0 = never (day)
tx_interval:   28  : transfusion interval (day)
tx_units :     2   : units per transfusion (-)
tx_cells :     2   : red cells per unit (1e12)
tau_dose :   0.5   : infusion duration (day)
k_don    : 0.0181818 : donor red cell loss rate (1/day)
 @annotated

$CMT @annotated
// --- 9 red cell age cohorts, edges 0/4/9/16/26/42/64/92/130/200 days
N1 : cohort 1 cell count (1e12/L)
N2 : cohort 2 cell count (1e12/L)
N3 : cohort 3 cell count (1e12/L)
N4 : cohort 4 cell count (1e12/L)
N5 : cohort 5 cell count (1e12/L)
N6 : cohort 6 cell count (1e12/L)
N7 : cohort 7 cell count (1e12/L)
N8 : cohort 8 cell count (1e12/L)
N9 : cohort 9 cell count (1e12/L)
NA1 : cohort 1 count x membrane area (1e12 um2/L)
NA2 : cohort 2 count x membrane area (1e12 um2/L)
NA3 : cohort 3 count x membrane area (1e12 um2/L)
NA4 : cohort 4 count x membrane area (1e12 um2/L)
NA5 : cohort 5 count x membrane area (1e12 um2/L)
NA6 : cohort 6 count x membrane area (1e12 um2/L)
NA7 : cohort 7 count x membrane area (1e12 um2/L)
NA8 : cohort 8 count x membrane area (1e12 um2/L)
NA9 : cohort 9 count x membrane area (1e12 um2/L)
NV1 : cohort 1 count x volume (1e12 um3/L)
NV2 : cohort 2 count x volume (1e12 um3/L)
NV3 : cohort 3 count x volume (1e12 um3/L)
NV4 : cohort 4 count x volume (1e12 um3/L)
NV5 : cohort 5 count x volume (1e12 um3/L)
NV6 : cohort 6 count x volume (1e12 um3/L)
NV7 : cohort 7 count x volume (1e12 um3/L)
NV8 : cohort 8 count x volume (1e12 um3/L)
NV9 : cohort 9 count x volume (1e12 um3/L)
NH1 : cohort 1 count x haemoglobin (1e12 pg/L)
NH2 : cohort 2 count x haemoglobin (1e12 pg/L)
NH3 : cohort 3 count x haemoglobin (1e12 pg/L)
NH4 : cohort 4 count x haemoglobin (1e12 pg/L)
NH5 : cohort 5 count x haemoglobin (1e12 pg/L)
NH6 : cohort 6 count x haemoglobin (1e12 pg/L)
NH7 : cohort 7 count x haemoglobin (1e12 pg/L)
NH8 : cohort 8 count x haemoglobin (1e12 pg/L)
NH9 : cohort 9 count x haemoglobin (1e12 pg/L)
NB1 : cohort 1 count x band 3 copies (1e18/L)
NB2 : cohort 2 count x band 3 copies (1e18/L)
NB3 : cohort 3 count x band 3 copies (1e18/L)
NB4 : cohort 4 count x band 3 copies (1e18/L)
NB5 : cohort 5 count x band 3 copies (1e18/L)
NB6 : cohort 6 count x band 3 copies (1e18/L)
NB7 : cohort 7 count x band 3 copies (1e18/L)
NB8 : cohort 8 count x band 3 copies (1e18/L)
NB9 : cohort 9 count x band 3 copies (1e18/L)
// --- erythropoiesis
EPO  : plasma erythropoietin (mIU/mL)
PROG : erythroid progenitor pool (1e12/L)
ERB  : erythroblast pool (1e12/L)
RETM : marrow reticulocytes (1e12/L)
RETB : blood reticulocytes (1e12/L)
// --- transfused cells
NDON : transfused donor red cells (1e12/L)
// --- spleen
SPLV : spleen volume (mL)
CORD : cordal red cell pool (1e12/L)
MAC  : red pulp macrophage activation (relative)
// --- haem catabolism
HPT  : haptoglobin (mg/dL)
FHB  : free plasma haemoglobin (mg/dL)
BILU : unconjugated bilirubin (mg/dL)
BILC : conjugated bilirubin (mg/dL)
BILE : biliary bilirubin pool (mg)
STONE: gallstone index (0-1)
LDH  : lactate dehydrogenase (U/L)
// --- iron
HEPC : hepcidin (relative)
FERR : serum ferritin (ng/mL)
FELIV: liver iron concentration (mg/g dw)
FESPL: splenic macrophage iron (mg)
// --- mitapivat
MGUT : mitapivat gut depot (mg)
MCEN : mitapivat central (mg)
MPER : mitapivat peripheral (mg)
ATP  : red cell ATP (relative)
DPG  : red cell 2,3-DPG (relative)
// --- insults
PARVO: parvovirus B19 erythroid arrest (0-1)
FOL  : folate status (relative)

$GLOBAL
#define NCOH 9
static const double AGEMID[NCOH] = {2.0, 6.5, 12.5, 21.0, 34.0, 53.0,
                                    78.0, 111.0, 165.0};
static const double ADV[NCOH]    = {1.0/4.0, 1.0/5.0, 1.0/7.0, 1.0/10.0,
                                    1.0/16.0, 1.0/22.0, 1.0/28.0, 1.0/38.0,
                                    1.0/70.0};
#define SQRTPI 1.7724538509055159
#define TWOPI3 2.0943951023931953

// closed-form minimum cylindrical diameter (Canham & Burton geometry)
inline double dcrit(double A, double V) {
  if (A < 1e-6) return 0.0;
  double Vs = pow(A, 1.5) / (6.0 * SQRTPI);
  double s  = V / Vs;
  if (s > 1.0) s = 1.0;
  if (s < 1e-9) s = 1e-9;
  double th = acos(-s);
  return 2.0 * sqrt(A / M_PI) * cos(th / 3.0 - TWOPI3);
}
inline double vsph(double A) { return pow(A, 1.5) / (6.0 * SQRTPI); }

$MAIN
if (NEWIND < 2) {
  // steady-ish initial condition; run 900-1400 days to the fixed point
  double Aent = A0 * (1.0 - a_ent_def * fdef * 0.42);
  double wsum = 0.0, w[NCOH];
  double width[NCOH] = {4,5,7,10,16,22,28,38,70};
  for (int i = 0; i < NCOH; ++i) { w[i] = width[i]*exp(-AGEMID[i]/110.0);
                                   wsum += w[i]; }
  N_1  = 5.0*w[0]/wsum;  N_2 = 5.0*w[1]/wsum;  N_3 = 5.0*w[2]/wsum;
  N_4  = 5.0*w[3]/wsum;  N_5 = 5.0*w[4]/wsum;  N_6 = 5.0*w[5]/wsum;
  N_7  = 5.0*w[6]/wsum;  N_8 = 5.0*w[7]/wsum;  N_9 = 5.0*w[8]/wsum;
  NA1_0 = N_1*Aent; NA2_0 = N_2*Aent; NA3_0 = N_3*Aent;
  NA4_0 = N_4*Aent; NA5_0 = N_5*Aent; NA6_0 = N_6*Aent;
  NA7_0 = N_7*Aent; NA8_0 = N_8*Aent; NA9_0 = N_9*Aent;
  NV1_0 = N_1*V0; NV2_0 = N_2*V0; NV3_0 = N_3*V0;
  NV4_0 = N_4*V0; NV5_0 = N_5*V0; NV6_0 = N_6*V0;
  NV7_0 = N_7*V0; NV8_0 = N_8*V0; NV9_0 = N_9*V0;
  NH1_0 = N_1*H0; NH2_0 = N_2*H0; NH3_0 = N_3*H0;
  NH4_0 = N_4*H0; NH5_0 = N_5*H0; NH6_0 = N_6*H0;
  NH7_0 = N_7*H0; NH8_0 = N_8*H0; NH9_0 = N_9*H0;
  double Bent = B0*(1.0 - f_b3ves*(1.0 - Aent/A0));
  NB1_0 = N_1*Bent; NB2_0 = N_2*Bent; NB3_0 = N_3*Bent;
  NB4_0 = N_4*Bent; NB5_0 = N_5*Bent; NB6_0 = N_6*Bent;
  NB7_0 = N_7*Bent; NB8_0 = N_8*Bent; NB9_0 = N_9*Bent;
  N1_0 = N_1; N2_0 = N_2; N3_0 = N_3; N4_0 = N_4; N5_0 = N_5;
  N6_0 = N_6; N7_0 = N_7; N8_0 = N_8; N9_0 = N_9;
  EPO_0 = 10.0;  PROG_0 = prod0/k_prog;  ERB_0 = prod0/k_erb;
  RETM_0 = prod0/k_retm;  RETB_0 = prod0*tau_ret0;
  SPLV_0 = SPL_base;  MAC_0 = 1.0;  HPT_0 = 100.0;
  BILU_0 = 0.6;  BILC_0 = 0.15;  LDH_0 = LDH0;
  HEPC_0 = HEPC0;  FERR_0 = 90.0;  FELIV_0 = 0.8;  FESPL_0 = 300.0;
  ATP_0 = 1.0;  DPG_0 = 1.0;  FOL_0 = 1.0;
}

$ODE
double Nc[NCOH]  = {N1,N2,N3,N4,N5,N6,N7,N8,N9};
double NAc[NCOH] = {NA1,NA2,NA3,NA4,NA5,NA6,NA7,NA8,NA9};
double NVc[NCOH] = {NV1,NV2,NV3,NV4,NV5,NV6,NV7,NV8,NV9};
double NHc[NCOH] = {NH1,NH2,NH3,NH4,NH5,NH6,NH7,NH8,NH9};
double NBc[NCOH] = {NB1,NB2,NB3,NB4,NB5,NB6,NB7,NB8,NB9};

double A[NCOH], V[NCOH], H[NCOH], Bq[NCOH], Dc[NCOH], sph[NCOH];
double pslow[NCOH], tauc[NCOH], Rcord[NCOH], kv[NCOH];
double vA[NCOH], vV[NCOH], vH[NCOH], vB[NCOH];
double hgeom[NCOH], hops[NCOH], hliv[NCOH], hlys[NCOH], hsen[NCOH], hz[NCOH];
double IgG[NCOH];

// ---- spleen: functional mass and traffic
double spl_rel = spl_frac * SPLV / SPL_base;
double spl_c   = spl_rel > spl_ops_cap ? spl_ops_cap : spl_rel;
double spl_fl  = pow(spl_rel, spl_flow_exp);
if (spl_fl > spl_flow_cap) spl_fl = spl_flow_cap;
double f_pass  = f_pass0 * spl_fl;
double atp     = ATP < 0.25 ? 0.25 : ATP;
double atp_f   = pow(atp, -p_atp);

double Ntot = 0.0, HBt = 0.0, VOLt = 0.0;
double destN = 0.0, destHb = 0.0, lysHb = 0.0, Wspl = 0.0, Rpool = 0.0;

for (int i = 0; i < NCOH; ++i) {
  double n = Nc[i] > 1e-14 ? Nc[i] : 1e-14;
  A[i] = NAc[i]/n;  if (A[i] < 20.0) A[i] = 20.0;
  V[i] = NVc[i]/n;  if (V[i] < 15.0) V[i] = 15.0;
  H[i] = NHc[i]/n;
  Bq[i] = NBc[i]/n;

  Dc[i]  = dcrit(A[i], V[i]);
  sph[i] = V[i]/vsph(A[i]);

  pslow[i] = 1.0/(1.0 + exp(-(Dc[i]-D50)/wD));
  double ex = (Dc[i]-Dc_ref)/w_esc;
  if (ex >  8.0) ex =  8.0;
  if (ex < -8.0) ex = -8.0;
  double vk = visc_k*(100.0*H[i]/V[i] - MCHC_ref);
  if (vk >  4.0) vk =  4.0;
  if (vk < -3.0) vk = -3.0;
  tauc[i]  = tau0*exp(ex)*exp(vk);
  Rcord[i] = f_pass*pslow[i]*tauc[i];
  if (Rcord[i] > R_MAX) Rcord[i] = R_MAX;

  // membrane loss: intrinsic + cordal amplification, ATP-modulated
  kv[i] = (kv_base + kv_def*fdef) * (1.0 + (cordamp-1.0)*Rcord[i]) * atp_f;
  vA[i] = A[i]*kv[i];
  double Vfloor = 100.0*H[i]/MCHC_MAX;
  double dv = V[i]-Vfloor; if (dv < 0.0) dv = 0.0;
  double kd = (kd_base + kd_def*fdef + kd_cord*Rcord[i])
              / (atp < 0.3 ? 0.3 : atp);
  vV[i] = vA[i]*r_ves3 + dv*kd;
  vH[i] = vA[i]*r_ves3*(H[i]/V[i])*f_hbves;
  vB[i] = (vA[i]/A[i])*Bq[i]*f_b3ves;

  // clearance 1: geometric, splenic only
  double pidest = 1.0 - exp(-k_ph*MAC*tauc[i]*exp(vk));
  hgeom[i] = f_pass*pslow[i]*pidest;
  // clearance 2: opsonic, band 3 density -> natural anti-band-3 IgG
  double rho = (Bq[i]/A[i])/Bden_ref;
  double ce  = rho > 1.0 ? rho-1.0 : 0.0;
  IgG[i] = IgG_bg + IgG_span*ce*ce/(K_cl*K_cl + ce*ce);
  double ig = IgG[i]-IgG_thr; if (ig < 0.0) ig = 0.0;
  double fo = pow(ig/IgG50, m_ops);
  hops[i] = k_ops*(w_spl_ops*spl_c*MAC + w_liv_ops)*fo/(1.0+fo);
  // clearance 3: hepatic geometric
  hliv[i] = k_liv_geo/(1.0 + exp(-(Dc[i]-D_liv)/w_liv));
  // clearance 4: intravascular lysis of the nearly spherical
  double sl = sph[i]-s_lys; if (sl < 0.0) sl = 0.0;
  hlys[i] = k_lys*sl/(1.0-s_lys);
  // clearance 5: senescence -- molecular tagging, geometry-blind
  hsen[i] = k_sen*pow(AGEMID[i]/tau50, m_sen);

  hz[i] = hgeom[i]+hops[i]+hliv[i]+hlys[i]+hsen[i];

  Ntot  += Nc[i];
  HBt   += NHc[i];
  VOLt  += NVc[i];
  destN += hz[i]*Nc[i];
  destHb+= hz[i]*NHc[i];
  lysHb += hlys[i]*NHc[i];
  Wspl  += (hgeom[i] + f_sen_spl*hsen[i] + f_lys_spl*hlys[i])*Nc[i]
           + (w_spl_ops*spl_c/(w_spl_ops+w_liv_ops))*hops[i]*Nc[i];
  Rpool += Rcord[i]*Nc[i];
}
destN  += ADV[NCOH-1]*Nc[NCOH-1];
destHb += ADV[NCOH-1]*NHc[NCOH-1];

double RBCt = Ntot + RETB + NDON;
double Hbg  = (HBt + RETB*H0 + NDON*30.0)/10.0;

// ---- reticulocyte input
double tr = EPO/EPO_norm/12.0; if (tr > 1.7) tr = 1.7;
double tau_ret = tau_ret0*(1.0 + 0.9*tr);
double prod_in = RETB/tau_ret;
double Aent    = A0*(1.0 - a_ent_def*fdef*0.42);

// ---- cohort ODEs (upwind advection of the extensive quantities)
double inN[NCOH], inA[NCOH], inV[NCOH], inH[NCOH], inB[NCOH];
inN[0]=prod_in;         inA[0]=prod_in*Aent;  inV[0]=prod_in*V0;
inH[0]=prod_in*H0;
// band 3 at release tracks the membrane already missing at release
double B_ent = B0*(1.0 - f_b3ves*(1.0 - Aent/A0));
inB[0]=prod_in*B_ent;
for (int i = 1; i < NCOH; ++i) {
  inN[i]=ADV[i-1]*Nc[i-1];  inA[i]=ADV[i-1]*NAc[i-1];
  inV[i]=ADV[i-1]*NVc[i-1]; inH[i]=ADV[i-1]*NHc[i-1];
  inB[i]=ADV[i-1]*NBc[i-1];
}
double dN[NCOH], dNA[NCOH], dNV[NCOH], dNH[NCOH], dNB[NCOH];
for (int i = 0; i < NCOH; ++i) {
  dN[i]  = inN[i] - ADV[i]*Nc[i]  - hz[i]*Nc[i];
  dNA[i] = inA[i] - ADV[i]*NAc[i] - hz[i]*NAc[i] - Nc[i]*vA[i];
  dNV[i] = inV[i] - ADV[i]*NVc[i] - hz[i]*NVc[i] - Nc[i]*vV[i];
  dNH[i] = inH[i] - ADV[i]*NHc[i] - hz[i]*NHc[i] - Nc[i]*vH[i];
  dNB[i] = inB[i] - ADV[i]*NBc[i] - hz[i]*NBc[i] - Nc[i]*vB[i];
}
dxdt_N1=dN[0]; dxdt_N2=dN[1]; dxdt_N3=dN[2]; dxdt_N4=dN[3]; dxdt_N5=dN[4];
dxdt_N6=dN[5]; dxdt_N7=dN[6]; dxdt_N8=dN[7]; dxdt_N9=dN[8];
dxdt_NA1=dNA[0]; dxdt_NA2=dNA[1]; dxdt_NA3=dNA[2]; dxdt_NA4=dNA[3];
dxdt_NA5=dNA[4]; dxdt_NA6=dNA[5]; dxdt_NA7=dNA[6]; dxdt_NA8=dNA[7];
dxdt_NA9=dNA[8];
dxdt_NV1=dNV[0]; dxdt_NV2=dNV[1]; dxdt_NV3=dNV[2]; dxdt_NV4=dNV[3];
dxdt_NV5=dNV[4]; dxdt_NV6=dNV[5]; dxdt_NV7=dNV[6]; dxdt_NV8=dNV[7];
dxdt_NV9=dNV[8];
dxdt_NH1=dNH[0]; dxdt_NH2=dNH[1]; dxdt_NH3=dNH[2]; dxdt_NH4=dNH[3];
dxdt_NH5=dNH[4]; dxdt_NH6=dNH[5]; dxdt_NH7=dNH[6]; dxdt_NH8=dNH[7];
dxdt_NH9=dNH[8];
dxdt_NB1=dNB[0]; dxdt_NB2=dNB[1]; dxdt_NB3=dNB[2]; dxdt_NB4=dNB[3];
dxdt_NB5=dNB[4]; dxdt_NB6=dNB[5]; dxdt_NB7=dNB[6]; dxdt_NB8=dNB[7];
dxdt_NB9=dNB[8];

// ---- erythropoiesis
double hbe = Hbg < 2.0 ? 2.0 : Hbg;
double EPO_t = pow(10.0, EPO_a - EPO_b*hbe);
dxdt_EPO = k_epo*(EPO_t - EPO);
double e    = EPO - EPO_norm; if (e < 0.0) e = 0.0;
double stim = 1.0 + (Emax_mar-1.0)*e/(K_epo+e);
double parv = 1.0 - parvo_supp*PARVO; if (parv < 0.0) parv = 0.0;
double folf = (FOL/(0.35+FOL))*1.35; if (folf > 1.0) folf = 1.0;
dxdt_PROG = prod0*stim*parv*folf - k_prog*PROG;
dxdt_ERB  = k_prog*PROG - k_erb*ERB;
double shift = 1.0 + (shift_max-1.0)*e/(K_epo+e);
dxdt_RETM = k_erb*ERB - k_retm*shift*RETM;
dxdt_RETB = k_retm*shift*RETM - prod_in;
double sx = stim-1.0; if (sx < 0.0) sx = 0.0;
dxdt_FOL  = k_fol*(fol_ok - FOL) - 0.010*sx*FOL;

// ---- transfusion
double txin = 0.0;
if (tx_start >= 0.0 && SOLVERTIME >= tx_start) {
  double ph = fmod(SOLVERTIME - tx_start, tx_interval);
  if (ph < tau_dose) txin = tx_units*tx_cells/tau_dose/BV;
}
dxdt_NDON = txin - k_don*NDON;
destHb += k_don*NDON*30.0;
double destHb_g = destHb*BV;
double lysHb_g  = lysHb*BV;

// ---- spleen
double load  = Wspl/W0;
double lx    = load-1.0; if (lx < 0.0) lx = 0.0;
double SPL_t = SPL_base*(1.0+g_spl*lx); if (SPL_t > SPL_max) SPL_t = SPL_max;
double gr    = 1.0 - SPLV/(SPL_t > 1.0 ? SPL_t : 1.0); if (gr < 0.0) gr = 0.0;
dxdt_SPLV = k_spl*(SPL_t - SPLV)*spl_frac + k_regrow*SPLV*gr;
dxdt_CORD = Rpool - CORD;
dxdt_MAC  = k_mac*(1.0 + a_mac*lx/(1.0+lx/b_mac) - MAC);

// ---- haem catabolism and bilirubin
double BR_prod = BR_PER_G*destHb_g;
dxdt_FHB  = lysHb_g*1000.0/BV/10.0 - k_fhb*FHB;
dxdt_HPT  = k_hp_syn - k_hp_deg*HPT - k_hp_bind*FHB*HPT/(10.0+HPT);
double conj = Vmax_ugt*ugt_f*BILU/(Km_ugt+BILU);
dxdt_BILU = (BR_prod - conj)/(10.0*Vd_bil);
dxdt_BILC = conj/(10.0*Vd_bil) - k_bilc*BILC;
dxdt_BILE = k_bilc*BILC*10.0*Vd_bil - k_bile*BILE;
double ucb = f_ucb_bile*BILE*(1.0 + 0.5*(BILU/0.6 - 1.0));
if (ucb < 0.0) ucb = 0.0;
double stfree = 1.0 - (STONE > 1.0 ? 1.0 : STONE);
dxdt_STONE = k_stone*ucb*stfree;
dxdt_LDH = k_ldh*(LDH0 + 210.0*(destN/(5.0/120.0)-1.0) + 900.0*lysHb_g - LDH);

// ---- iron
double erfe = sx;
double tx_fe_d = 0.0;
if (tx_start >= 0.0 && SOLVERTIME >= tx_start)
  tx_fe_d = tx_units*tx_fe/tx_interval;
dxdt_HEPC = k_hepc*(HEPC0*(1.0+0.9*(FELIV/0.8-1.0))/(1.0+erfe_k*erfe) - HEPC);
double hp = HEPC < 0.05 ? 0.05 : HEPC;
double abs_fe = k_fe_abs*1000.0/hp;
dxdt_FELIV = (abs_fe + 0.55*tx_fe_d - k_chel*FELIV*12.0 - 1.0)*k_fe_liv;
dxdt_FESPL = 0.45*tx_fe_d + 0.02*destN*fe_per_cell*BV - 0.02*FESPL
             - k_chel*FESPL*0.6;
dxdt_FERR  = k_ferr*(30.0 + 120.0*FELIV/0.8 - FERR);

// ---- mitapivat PK/PD (BID)
double dose = 0.0;
if (dose_m > 0.0) { double ph = fmod(SOLVERTIME, 0.5);
                    if (ph < 0.02) dose = dose_m/0.02; }
double conc = MCEN/Vc_m*1000.0;
dxdt_MGUT = dose - ka_m*24.0*MGUT;
dxdt_MCEN = ka_m*24.0*MGUT - CL_m*24.0*MCEN/Vc_m
            - Q_m*24.0*(MCEN/Vc_m - MPER/Vp_m);
dxdt_MPER = Q_m*24.0*(MCEN/Vc_m - MPER/Vp_m);
double drv = conc/(EC50_m+conc);
double rfrac = Rpool/(RBCt > 1e-6 ? RBCt : 1e-6); if (rfrac > 1.0) rfrac = 1.0;
dxdt_ATP = k_atp*(1.0 + Emax_atp*drv - 0.10*rfrac - ATP);
dxdt_DPG = k_dpg*(1.0 - Emax_dpg*drv - DPG);

// ---- parvovirus B19
double onp = 0.0;
if (parvo_t >= 0.0 && SOLVERTIME >= parvo_t
    && SOLVERTIME < parvo_t+parvo_dur) onp = 1.0;
dxdt_PARVO = (parvo_t >= 0.0) ? 3.0*(onp-PARVO) : -3.0*PARVO;

$TABLE
double Nt = 0.0, HBs = 0.0, VLs = 0.0, dN_ = 0.0, dHb_ = 0.0;
double wA = 0.0, wV = 0.0, wD_ = 0.0, wS = 0.0, wI = 0.0, wR = 0.0;
double gG = 0.0, gO = 0.0, gS = 0.0, gL = 0.0, gY = 0.0;
double Nq[NCOH]  = {N1,N2,N3,N4,N5,N6,N7,N8,N9};
double NAq[NCOH] = {NA1,NA2,NA3,NA4,NA5,NA6,NA7,NA8,NA9};
double NVq[NCOH] = {NV1,NV2,NV3,NV4,NV5,NV6,NV7,NV8,NV9};
double NHq[NCOH] = {NH1,NH2,NH3,NH4,NH5,NH6,NH7,NH8,NH9};
double NBq[NCOH] = {NB1,NB2,NB3,NB4,NB5,NB6,NB7,NB8,NB9};
double sr = spl_frac*SPLV/SPL_base;
double sc = sr > spl_ops_cap ? spl_ops_cap : sr;
double sf = pow(sr, spl_flow_exp); if (sf > spl_flow_cap) sf = spl_flow_cap;
double fp = f_pass0*sf;
for (int i = 0; i < NCOH; ++i) {
  double n = Nq[i] > 1e-14 ? Nq[i] : 1e-14;
  double Ai = NAq[i]/n, Vi = NVq[i]/n, Hi = NHq[i]/n, Bi = NBq[i]/n;
  if (Ai < 20.0) Ai = 20.0;
  if (Vi < 15.0) Vi = 15.0;
  double Di = dcrit(Ai, Vi), si = Vi/vsph(Ai);
  double ps = 1.0/(1.0+exp(-(Di-D50)/wD));
  double ex = (Di-Dc_ref)/w_esc; if (ex>8.0) ex=8.0; if (ex<-8.0) ex=-8.0;
  double vk = visc_k*(100.0*Hi/Vi - MCHC_ref);
  if (vk >  4.0) vk =  4.0;
  if (vk < -3.0) vk = -3.0;
  double tc = tau0*exp(ex)*exp(vk);
  double Rc = fp*ps*tc; if (Rc > R_MAX) Rc = R_MAX;
  double pd = 1.0-exp(-k_ph*MAC*tc*exp(vk));
  double hg = fp*ps*pd;
  double rho = (Bi/Ai)/Bden_ref; double ce = rho>1.0 ? rho-1.0 : 0.0;
  double Ig = IgG_bg + IgG_span*ce*ce/(K_cl*K_cl+ce*ce);
  double ig = Ig-IgG_thr; if (ig<0.0) ig=0.0;
  double fo = pow(ig/IgG50, m_ops);
  double ho = k_ops*(w_spl_ops*sc*MAC+w_liv_ops)*fo/(1.0+fo);
  double hl = k_liv_geo/(1.0+exp(-(Di-D_liv)/w_liv));
  double sl = si-s_lys; if (sl<0.0) sl=0.0;
  double hy = k_lys*sl/(1.0-s_lys);
  double hs = k_sen*pow(AGEMID[i]/tau50, m_sen);
  double ht = hg+ho+hl+hy+hs;
  Nt += Nq[i]; HBs += NHq[i]; VLs += NVq[i];
  dN_ += ht*Nq[i]; dHb_ += ht*NHq[i];
  wA += Nq[i]*Ai; wV += Nq[i]*Vi; wD_ += Nq[i]*Di; wS += Nq[i]*si;
  wI += Nq[i]*Ig; wR += Nq[i]*Rc;
  gG += hg*Nq[i]; gO += ho*Nq[i]; gS += hs*Nq[i];
  gL += hl*Nq[i]; gY += hy*Nq[i];
}
dN_  += ADV[NCOH-1]*N9;
dHb_ += ADV[NCOH-1]*NH9;
double RBC  = Nt + RETB + NDON;
double Hb   = (HBs + RETB*H0 + NDON*30.0)/10.0;
double VOL  = VLs + RETB*V0*V_RET_MULT + NDON*88.0;
double Hct  = VOL/10.0;                       // %
double MCV  = VOL/(RBC > 1e-9 ? RBC : 1e-9);
double MCH  = (HBs + RETB*H0 + NDON*30.0)/(RBC > 1e-9 ? RBC : 1e-9);
double MCHC = 100.0*MCH/(MCV > 1e-9 ? MCV : 1e-9);
double RETpct = 100.0*RETB/(RBC > 1e-9 ? RBC : 1e-9);
double RETabs = 1000.0*RETB;
double LIFE = Nt/(dN_ > 1e-12 ? dN_ : 1e-12);
double AREA = wA/(Nt > 1e-12 ? Nt : 1e-12);
double VOLc = wV/(Nt > 1e-12 ? Nt : 1e-12);
double DCRIT= wD_/(Nt > 1e-12 ? Nt : 1e-12);
double SPHER= wS/(Nt > 1e-12 ? Nt : 1e-12);
double IGGC = wI/(Nt > 1e-12 ? Nt : 1e-12);
double RCORD= wR/(Nt > 1e-12 ? Nt : 1e-12);
double EMA  = AREA/A0;
double TBIL = BILU + BILC;
double BRPROD = BR_PER_G*dHb_*BV;
double FGEOM = gG/(dN_ > 1e-12 ? dN_ : 1e-12);
double FOPS  = gO/(dN_ > 1e-12 ? dN_ : 1e-12);
double FSEN  = gS/(dN_ > 1e-12 ? dN_ : 1e-12);
double FLIV  = gL/(dN_ > 1e-12 ? dN_ : 1e-12);
double FLYS  = gY/(dN_ > 1e-12 ? dN_ : 1e-12);
// osmotic fragility of the mean cell: %NaCl at which V reaches V_sph
double vs_  = 0.75*MCH + 0.04*VOLc;           // osmotically inactive volume
double vw_  = VOLc - vs_;
double vmax_= vsph(AREA);
double MCF  = (vmax_ > vs_) ? 290.0*vw_/(vmax_-vs_)*(0.9/308.0) : 0.9;
double SPLEEN = spl_frac*SPLV;
double STONEPCT = 100.0*STONE;
double MITA = MCEN/Vc_m*1000.0;

$CAPTURE @annotated
Hb      : haemoglobin (g/dL)
Hct     : haematocrit (%)
RBC     : red cell count (1e12/L)
MCV     : mean cell volume (fL)
MCH     : mean cell haemoglobin (pg)
MCHC    : mean cell haemoglobin concentration (g/dL)
RETpct  : reticulocytes (%)
RETabs  : absolute reticulocytes (1e9/L)
LIFE    : mean red cell lifespan (day)
AREA    : mean membrane area (um2)
VOLc    : mean cell volume of mature cells (um3)
DCRIT   : mean minimum cylindrical diameter (um)
SPHER   : mean sphericity V/V_sph (-)
EMA     : EMA binding relative to normal (-)
MCF     : 50% osmotic lysis point (%NaCl)
IGGC    : red-cell-bound IgG (molecules/cell)
RCORD   : fraction of time spent in the splenic cords (-)
FGEOM   : share of destruction that is geometric (-)
FOPS    : share of destruction that is opsonic (-)
FSEN    : share of destruction that is senescence (-)
FLIV    : share of destruction that is hepatic (-)
FLYS    : share of destruction that is intravascular lysis (-)
TBIL    : total bilirubin (mg/dL)
BRPROD  : bilirubin production (mg/day)
SPLEEN  : spleen volume (mL)
STONEPCT: gallstone index (%)
MITA    : mitapivat plasma concentration (ng/mL)
