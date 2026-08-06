# =====================================================================
#  Bile Acid Diarrhoea / Bile Acid Malabsorption (BAM)
#  Quantitative Systems Pharmacology model  —  mrgsolve
# ---------------------------------------------------------------------
#  ORGANISING IDEA
#  ---------------
#  The enterohepatic circulation is a negative-feedback loop whose
#  SENSOR (ileal FXR -> FGF19) measures the ABSORBED flux, while the
#  SYMPTOM is produced by the SPILLED flux.  At steady state, spillage
#  equals hepatic synthesis exactly, so
#
#         colonic bile-acid load  =  S  =  hepatic synthesis rate
#
#  and malabsorption BY ITSELF cannot raise the colonic load.  It raises
#  it only through the feedback it disinhibits.  Two independent lesions
#  are therefore needed to write the disease:
#
#         phi    surviving fraction of ileal ASBT absorptive capacity
#         kappa  gain of the ileal FGF19 sensor
#
#  and the two routine tests read them SEPARATELY:
#         SeHCAT 7-day retention  ->  phi   alone
#         serum C4 / faecal BA    ->  S     (phi AND kappa together)
#
#  Set FEEDBACK = 0 to FREEZE the loop at its normal operating point.
#  That single switch is the model's central experiment.  Verified result:
#  with the loop frozen, dropping phi from 1.00 to 0.05 moves the colonic
#  bile-acid load by +0.2 % and the stool water by -0.0 %, while SeHCAT
#  still collapses from 50.2 % to 0.00 %.  With the loop live the same
#  lesion multiplies the colonic load 5.56x and the stool water 6.81x.
#  The diarrhoea is made by the answer, not by the injury.
#
#  THERAPY CLASSES (see $PARAM dose_* below)
#    A  sequestrants  (colesevelam, colestyramine) — bind the luminal
#       load, but also remove the sensor's ligand: C4 RISES, and the
#       dose-response is therefore shallow.  Structural escape.
#    B  FXR agonists  (obeticholic acid, tropifexor) — ADDITIVE at the
#       receptor, independent of phi; the only class that lowers S.
#    C  ASBT inhibitors (elobixibat, odevixibat) — a pharmacological
#       phi lesion; the model must reproduce iatrogenic BAM.
#    D  transit agents (loperamide, ondansetron) — slow colonic transit.
#       NOTE: the model REFUTED the hypothesis that these break a vicious
#       cycle.  Measuring the opened loop gives a NEGATIVE loop gain
#       (L = -0.18 to -1.6): faster transit dilutes and evacuates the
#       secretagogue faster than the extra ileal spillage replaces it, so
#       bile acid diarrhoea cannot latch.  Loperamide works by buying
#       colonic water-absorption time.  What DID survive is the ileal
#       half: longer ileal contact raises f, so faecal bile acid falls
#       too (2048 -> 1924 umol/day), not just stool water.
#    E  microbiome    (rifaximin) — lowers the POTENCY of a given load
#       by suppressing 7alpha-dehydroxylation (CA w=0.15 -> DCA w=1.15).
#
#  UNITS  time = hours, bile acid = umol, water = mL, drugs as noted.
#
#  Companion files
#    README.md              verified results, and an honest account of the
#                           two places the model misses (the 100 cm rule
#                           lands at 34 cm; drug effect sizes overshoot)
#    bam_qsp_model.dot      mechanistic map (137 nodes, 12 clusters)
#    bam_verify_python.py   independent re-implementation; produces
#                           bam_verification_output.txt, the source of
#                           every number quoted in README.md
#    bam_shiny_app.R        interactive dashboard (10 tabs)
# =====================================================================

library(mrgsolve)
library(dplyr)

code <- '
$PROB
# Bile Acid Diarrhoea / Malabsorption QSP model
- 51 ODE compartments (enterohepatic bile acid, FXR-FGF19-CYP7A1 axis,
  colonic microbial transformation, colonic secretion/motility, five
  drug classes, organ-system endpoints, SeHCAT tracer)
- Two lesion axes: phi (ileal ASBT capacity) and kappa (FGF19 gain)

$PLUGIN autodec

$GLOBAL
#define MEALAUC 1.8            // h/day, integral of the meal signal (3 x 1.2/2)
#define BAMPI   3.14159265358979323846

// Raised-cosine meal / CCK signal, 24 h period, peak 1.0
double meal_signal(double t) {
  double tod = fmod(t, 24.0);
  if (tod < 0) tod += 24.0;
  double m[3] = {8.0, 13.0, 19.0};
  double dur = 1.2, s = 0.0;
  for (int i = 0; i < 3; i++) {
    double d = tod - m[i];
    if (d < -12.0) d += 24.0;
    else if (d > 12.0) d -= 24.0;
    if (d >= -dur/2.0 && d <= dur/2.0)
      s += 0.5 * (1.0 + cos(2.0 * BAMPI * d / dur));
  }
  return s > 1.0 ? 1.0 : s;
}

$PARAM @annotated
// ================= THE TWO LESION AXES =============================
phi       : 1.0   : surviving ileal ASBT absorptive capacity (0-1)
kappa     : 1.0   : ileal FGF19 sensor gain (type-2 BAD < 1)
GBfun     : 1.0   : 1 = gallbladder present, 0 = cholecystectomy
ENTMASS   : 1.0   : functional ileal enterocyte mass (Crohn/radiation)
FEEDBACK  : 1.0   : 1 = FGF19->CYP7A1 loop live, 0 = FROZEN at normal

// ================= hepatic synthesis / CYP7A1 ======================
Vs        : 29.17 : umol/h synthesis at CYP7A1 = 1 (700 umol/day)
kout7     : 0.15  : /h CYP7A1 protein turnover (t1/2 4.6 h)
kin7      : 0.90  : /h max CYP7A1 transcription (ceiling 6x basal)
IC50f     : 118.049 : pg/mL FGF19 potency on CYP7A1 [CALIBRATED]
hf        : 1.5   : Hill coefficient, FGF19 arm
IC50s     : 1.0   : relative SHP potency on CYP7A1
hs        : 1.5   : Hill coefficient, SHP arm
circA     : 0.35  : circadian amplitude of CYP7A1 transcription
circPk    : 13.0  : h acrophase of CYP7A1 transcription
Rfnorm    : 0.2334 : frozen-loop clamp, day-averaged FGF19 repression
Rsnorm    : 0.7140 : frozen-loop clamp, day-averaged SHP repression

// ================= bile acid trafficking ===========================
kbsep     : 3.00  : /h hepatocyte -> canaliculus (BSEP)
fGB       : 0.72  : fasting fraction of bile diverted to gallbladder
kejF      : 0.130 : /h interdigestive gallbladder emptying
kejM      : 1.90  : /h CCK-driven gallbladder emptying
kduo      : 1.05  : /h duodenum -> jejunum
kjej      : 0.95  : /h jejunum -> ileum
kpass     : 0.020 : /h passive jejunal absorption
kile      : 1.00  : /h transit between ileal tanks and out
kost      : 6.00  : /h enterocyte -> portal (OSTa/b)
kpor      : 12.0  : /h portal -> hepatocyte (NTCP/OATP)
fspill    : 0.11  : portal fraction escaping hepatic first pass
ksys      : 0.90  : /h systemic -> liver re-uptake

// ================= ASBT (the phi axis) =============================
Vmax_asbt : 1164.94 : umol/h maximal uptake per ileal tank [CALIBRATED]
Km_asbt   : 6.0    : mM luminal half-saturating concentration
V_ile     : 75.0   : mL luminal volume per ileal tank
V_ent     : 60.0   : mL ileal enterocyte water
V_duo     : 250.0  : mL duodenal / meal volume
V_col     : 300.0  : mL colonic luminal water reference

// ================= FXR / FGF19 =====================================
fu_ent    : 0.010  : free fraction of enterocyte BA (IBABP buffers the rest)
EC50_ba   : 26.0   : uM FREE enterocyte bile acid at FXR (CDCA-equivalent)
EC50_oca  : 0.16   : uM obeticholic acid at FXR (~100x CDCA)
EC50_tro  : 0.0003 : uM tropifexor at FXR (sub-nM, non-steroidal)
kf19      : 805.08 : pg/mL/h FGF19 production at FXRi=1 [CALIBRATED]
kef19     : 1.386  : /h FGF19 elimination (t1/2 30 min)
kshp      : 0.50   : /h SHP formation
kdshp     : 0.50   : /h SHP decay
EC50_hep  : 900.0  : umol hepatic bile acid for half-max SHP
kc4       : 0.0295424 : (ng/mL)/(umol/h)/h C4 appearance [CALIBRATED]
kec4      : 0.0417 : /h C4 elimination (t1/2 16.6 h)

// ================= CYP8B1 / pool composition =======================
kin8      : 0.20   : /h CYP8B1 transcription
kout8     : 0.10   : /h CYP8B1 turnover
fCA0      : 0.55   : cholate fraction of new synthesis at set-point
E8ref     : 0.4385 : CYP8B1 level at the normal set-point

// ================= adaptive ASBT ===================================
kASBT     : 0.0060 : /h adaptation rate (t1/2 ~5 days)
ASBTmax   : 1.35   : ceiling of compensatory ASBT up-regulation

// ================= colonic microbial transformation ================
kbai      : 0.16   : /h 7alpha-dehydroxylation at BAI = 1
kdcaabs   : 0.0471678 : /h passive colonic absorption of DCA
kcdcaabs  : 0.032  : /h passive colonic absorption of CDCA
klcaout   : 0.030  : /h LCA sulfation / precipitation
kfec      : 0.0420822 : /h colonic emptying at TRANS = 1 [CALIBRATED]
kbsh      : 0.02   : /h BSH recovery rate
kbai_t    : 0.02   : /h BAI recovery rate

// ================= colonic secretion / motility ====================
wCA       : 0.15   : secretory potency weight, cholate
wCDCA     : 1.00   : secretory potency weight, chenodeoxycholate
wDCA      : 1.15   : secretory potency weight, deoxycholate
wLCA      : 0.30   : secretory potency weight, lithocholate
EC50_G    : 3.00   : mM-eq colonic secretory threshold (Mekjian 1971)
hG        : 2.5    : Hill coefficient of the secretory switch
Jsec      : 190.0  : mL/h maximal TGR5/CFTR-driven secretion
Qin       : 62.5   : mL/h ileocaecal water inflow (1500 mL/day)
kwabs     : 0.528619 : /h colonic water absorption [CALIBRATED]
kstool    : 0.030  : /h stool water evacuation
tauTR     : 3.0    : h transit-state time constant
aTR       : 1.55   : LOOP GAIN secretory drive -> transit acceleration
TRANSclamp: 0.0    : >0 pins TRANS, i.e. OPENS the transit loop
tauHT     : 1.0    : h 5-HT signal time constant
aHT       : 1.0    : gain of the EC-cell arm
tauPERM   : 12.0   : h permeability time constant
aPERM     : 0.55   : gain of the permeability arm

// ================= fat absorption / CMC ============================
CMC       : 1.5    : mM effective critical micellar concentration
hCMC      : 4.0    : Hill coefficient of micellar solubilisation
fatmax    : 0.97   : maximal fractional fat absorption
fatIn     : 100.0  : g/day dietary fat

// ================= systemic endpoints ==============================
LDL0      : 120.0  : mg/dL LDL-C at normal synthesis
expLDL    : 0.35   : exponent, LDL-C ~ (S0/S)^expLDL
tauLDL    : 336.0  : h LDL-C time constant
UOX0      : 30.0   : mg/day baseline urinary oxalate
aUOXfat   : 3.6    : mg/day per g faecal fat above 3 g
aUOXperm  : 9.0    : mg/day per unit permeability
tauUOX    : 48.0   : h urinary oxalate time constant
VITD0     : 30.0   : ng/mL 25-OH vitamin D at normal fat absorption
tauVITD   : 720.0  : h vitamin D time constant

// ================= drugs ===========================================
seq_cap   : 0.22   : umol bile acid bound per mg resin, IN VIVO effective
kbind     : 0.020  : /(mM.h) second-order binding onto free resin capacity
kseq12    : 1.30   : /h resin stomach -> small intestine
kseq23    : 0.55   : /h resin small intestine -> colon
kseq3out  : 0.05   : /h resin colon -> stool (residence ~20 h)
ka_oca    : 0.55   : /h obeticholic acid gut absorption
foca_pass : 0.40   : fraction of OCA uptake independent of ASBT
kOCAep    : 1.10   : /h OCA enterocyte -> plasma
kOCAout   : 0.055  : /h OCA systemic elimination
kOCAgo    : 0.15   : /h OCA luminal transit loss
IC50_elo  : 0.35   : uM luminal elobixibat (ASBT inhibition)
kelo_out  : 0.55   : /h elobixibat luminal clearance
ka_tro    : 0.75   : /h tropifexor absorption
ke_tro    : 0.06   : /h tropifexor elimination
V_trop    : 55.0   : L tropifexor volume of distribution
ftro_ile  : 0.55   : plasma -> ileal enterocyte partition
ka_lop    : 0.90   : /h loperamide absorption to effect site
ke_lop    : 0.14   : /h loperamide effect-site loss
IC50_lop  : 2.50   : mg loperamide at effect site
ke_ond    : 0.19   : /h ondansetron elimination
IC50_ond  : 3.2    : mg ondansetron
fOND      : 0.45   : max fraction of the 5-HT3 arm blocked
ke_rif    : 0.35   : /h rifaximin colonic clearance
IC50_rif  : 900.0  : mg rifaximin in colon
fRIF      : 0.60   : max suppression of 7alpha-dehydroxylation

$CMT @annotated
LIV     : hepatocyte conjugated bile acid (umol)
GB      : gallbladder bile acid (umol)
DUO     : duodenum + proximal jejunum lumen (umol)
JEJ     : mid small intestine lumen (umol)
ILE1    : proximal terminal ileum lumen (umol)
ILE2    : mid terminal ileum lumen (umol)
ILE3    : distal terminal ileum lumen (umol)
ENT     : ileal enterocyte cytosolic bile acid (umol)
POR     : portal blood bile acid (umol)
SYS     : systemic plasma bile acid (umol)
fCAp    : circulating pool fraction that is cholate
fDCAp   : circulating pool fraction that is deoxycholate
CCA     : colonic cholate (umol)
CCDCA   : colonic chenodeoxycholate (umol)
CDCA    : colonic deoxycholate (umol)
CLCA    : colonic lithocholate (umol)
BND     : bile acid bound to sequestrant, in transit (umol)
FEC     : cumulative faecal bile acid (umol)
FGF19   : plasma FGF19 (pg/mL)
SHP     : hepatic SHP (relative)
CYP7A1  : hepatic CYP7A1 protein (relative, 1 = normal)
C4      : serum 7a-hydroxy-4-cholesten-3-one (ng/mL)
CYP8B1  : hepatic CYP8B1 (relative)
ASBTx   : ileal ASBT protein (relative, adaptive)
WCOL    : colonic luminal water (mL)
TRANS   : colonic transit rate state (relative)
HT5     : mucosal 5-HT / EC-cell signal (relative)
PERM    : colonic permeability (relative)
STOOL   : cumulative stool water (mL)
BSHx    : bacterial bile salt hydrolase activity (relative)
BAI     : 7alpha-dehydroxylating capacity (relative)
SEQ1    : sequestrant, stomach (mg)
SEQ2    : sequestrant, small intestine (mg)
SEQ3    : sequestrant, colon (mg)
OCAg    : obeticholic acid, gut lumen (umol)
OCAe    : obeticholic acid, ileal enterocyte (umol)
OCAp    : obeticholic acid, plasma (umol)
ELOg    : elobixibat, ileal lumen (umol)
LOPg    : loperamide, gut depot (mg)
LOPe    : loperamide, myenteric effect site (mg)
RIFg    : rifaximin, colon (mg)
ONDc    : ondansetron, plasma (mg)
TROg    : tropifexor, gut (umol)
TROp    : tropifexor, plasma (umol)
HCHOL   : hepatic cholesterol pool (relative)
LDLC    : plasma LDL-C (mg/dL)
FECFAT  : cumulative faecal fat (g)
UOX     : urinary oxalate (mg/day)
VITD    : 25-OH vitamin D (ng/mL)
SEH     : SeHCAT tracer remaining in body (fraction of dose)
FATF    : faecal fat appearance rate, slow filter (g/day)

$MAIN
// ---- initial conditions: a healthy adult at the normal set point ----
LIV_0    = 500.0;  GB_0   = 2200.0; DUO_0  = 700.0;  JEJ_0  = 900.0;
ILE1_0   = 500.0;  ILE2_0 = 200.0;  ILE3_0 = 80.0;   ENT_0  = 200.0;
POR_0    = 100.0;  SYS_0  = 30.0;
fCAp_0   = 0.42;   fDCAp_0 = 0.21;
CCA_0    = 150.0;  CCDCA_0 = 120.0; CDCA_0 = 260.0;  CLCA_0 = 150.0;
FGF19_0  = 230.0;  SHP_0  = 1.0;    CYP7A1_0 = 1.0;  C4_0   = 18.0;
CYP8B1_0 = 0.4385;  ASBTx_0 = 1.0;
WCOL_0   = 180.0;  TRANS_0 = 1.0;   HT5_0  = 0.05;   PERM_0 = 1.0;
BSHx_0   = 1.0;    BAI_0  = 1.0;
HCHOL_0  = 1.0;    LDLC_0 = 120.0;  UOX_0  = 30.0;   VITD_0 = 30.0;
SEH_0    = 1.0;    FATF_0 = 3.5;

$ODE
double M = meal_signal(SOLVERTIME);

// ===================== hepatic synthesis ===========================
double circ = 1.0 + circA * cos(2.0*BAMPI*(fmod(SOLVERTIME,24.0) - circPk)/24.0);
double Rf, Rs;
if (FEEDBACK >= 0.5) {
  Rf = 1.0/(1.0 + pow(FGF19/IC50f, hf));
  Rs = 1.0/(1.0 + pow(SHP/IC50s, hs));
} else {
  // Loop FROZEN.  CYP7A1 integrates the DAY-AVERAGE repression, not the
  // fasting trough, so the clamp is the time-averaged Rf and Rs of the
  // healthy run -- otherwise the phi = 1 frozen run would not reproduce
  // the healthy state and the central comparison would be confounded.
  Rf = Rfnorm;
  Rs = Rsnorm;
}
dxdt_CYP7A1 = kin7*Rf*Rs*circ - kout7*CYP7A1;
double S = Vs * CYP7A1 * HCHOL;                       // umol/h

dxdt_CYP8B1 = kin8*Rf - kout8*CYP8B1;
double fCAnew = fCA0 * CYP8B1 / E8ref;
if (fCAnew < 0.25) fCAnew = 0.25;
if (fCAnew > 0.80) fCAnew = 0.80;

// ===================== bile secretion / gallbladder =================
double Jbsep = kbsep * LIV;
double to_gb = 0.0, Jeject = 0.0;
if (GBfun >= 0.5) {
  to_gb   = (1.0 - M) * fGB * Jbsep;
  Jeject  = (kejF + kejM*M) * GB;
}
double Jduo_in = Jbsep - to_gb + Jeject;

// ===================== small-intestinal transit =====================
double Jd2j  = kduo * (1.0 + 0.8*M) * DUO;
double Jj2i  = kjej * (1.0 + 0.8*M) * JEJ;
double Jpass = kpass * JEJ;

// ---- sequestrant: a CONSUMABLE capacity, not an equilibrium ---------
// Binding removes bile acid irreversibly and uses up resin capacity.  This
// is what makes the dose-response shallow: capacity is finite and most of
// it is spent in the small intestine against the 25 mmol/day cycling flux
// long before the resin reaches the colon.
double ILEtot = ILE1 + ILE2 + ILE3;
double SEQtot = SEQ2 + SEQ3;
double cap_free = seq_cap*SEQtot - BND;
if (cap_free < 0.0) cap_free = 0.0;
double Jbind_ile = 0.0;
if (cap_free > 0.0 && SEQtot > 1e-9 && ILEtot > 1e-9) {
  Jbind_ile = kbind * cap_free * (SEQ2/SEQtot) * (ILEtot/(3.0*V_ile));
  if (Jbind_ile > 0.60*ILEtot) Jbind_ile = 0.60*ILEtot;
}

// ---- ASBT uptake (the phi axis), three tanks in series -------------
double C_elo = ELOg/(3.0*V_ile) * 1000.0;             // uM
double Ielo  = 1.0/(1.0 + C_elo/IC50_elo);
double Vmax_eff = Vmax_asbt * phi * ASBTx * ENTMASS * Ielo;
double C1 = ILE1/V_ile;
double C2 = ILE2/V_ile;
double C3 = ILE3/V_ile;
double Ja1 = Vmax_eff * C1/(Km_asbt + C1);
double Ja2 = Vmax_eff * C2/(Km_asbt + C2);
double Ja3 = Vmax_eff * C3/(Km_asbt + C3);
double Jasbt = Ja1 + Ja2 + Ja3;

double kt  = kile * (1.0 + 0.6*M);
double J12 = kt*ILE1, J23 = kt*ILE2, J3c = kile*TRANS*ILE3;
double fb1 = ILE1/ILEtot, fb2 = ILE2/ILEtot, fb3 = ILE3/ILEtot;
if (ILEtot < 1e-9) { fb1 = fb2 = fb3 = 0.0; }
dxdt_ILE1 = Jj2i - Ja1 - J12 - Jbind_ile*fb1;
dxdt_ILE2 = J12  - Ja2 - J23 - Jbind_ile*fb2;
dxdt_ILE3 = J23  - Ja3 - J3c - Jbind_ile*fb3;

double Ji2c_free = J3c;

// ===================== enterocyte / FXR ==============================
double C_ent   = ENT/V_ent * 1000.0;                  // uM
double C_oca_e = OCAe/V_ent * 1000.0;                 // uM
double C_tro_e = ftro_ile * TROp / V_trop;            // uM
// Additive receptor occupancy.  Endogenous bile acid and OCA are both
// IBABP substrates and share the same free fraction; tropifexor is not a
// bile acid and reaches FXR from the blood, so it is NOT buffered -- which
// is exactly why its effect is independent of phi.
double occ  = fu_ent*C_ent/EC50_ba + fu_ent*C_oca_e/EC50_oca
              + C_tro_e/EC50_tro;
double FXRi = occ/(1.0 + occ);
double Jost = kost * ENT;
dxdt_ENT = Jasbt - Jost;

dxdt_FGF19 = kf19*kappa*ENTMASS*FXRi - kef19*FGF19;
double hepFXR = LIV/(EC50_hep + LIV);
dxdt_SHP = kshp*(hepFXR/0.5) - kdshp*SHP;
dxdt_C4  = kc4*S - kec4*C4;

// ===================== colon =========================================
double fCDCAp = 1.0 - fCAp - fDCAp;
if (fCDCAp < 0.0) fCDCAp = 0.0;
double in_CA   = Ji2c_free * fCAp;
double in_DCA  = Ji2c_free * fDCAp;
double in_CDCA = Ji2c_free * fCDCAp;

double RIFeff = 1.0 - fRIF * RIFg/(IC50_rif + RIFg);
double kdh    = kbai * BAI * RIFeff;

double tot_col = CCA + CCDCA + CDCA + CLCA;
double Jbind_col = 0.0;
if (cap_free > 0.0 && SEQtot > 1e-9 && tot_col > 1e-9) {
  Jbind_col = kbind * cap_free * (SEQ3/SEQtot) * (tot_col/V_col);
  if (Jbind_col > 0.60*tot_col) Jbind_col = 0.60*tot_col;
}
double fc1 = CCA/tot_col, fc2 = CCDCA/tot_col;
double fc3 = CDCA/tot_col, fc4 = CLCA/tot_col;
if (tot_col < 1e-9) { fc1 = fc2 = fc3 = fc4 = 0.0; }

double kfe = kfec * TRANS;
double Jabs_dca  = kdcaabs  * CDCA;
double Jabs_cdca = kcdcaabs * CCDCA;
dxdt_CCA   = in_CA   - kdh*CCA   - kfe*CCA              - Jbind_col*fc1;
dxdt_CCDCA = in_CDCA - kdh*CCDCA - kfe*CCDCA - Jabs_cdca - Jbind_col*fc2;
dxdt_CDCA  = in_DCA  + kdh*CCA   - kfe*CDCA  - Jabs_dca  - Jbind_col*fc3;
dxdt_CLCA  = kdh*CCDCA - kfe*CLCA - klcaout*CLCA        - Jbind_col*fc4;
dxdt_BND   = Jbind_ile + Jbind_col - kseq3out*BND;

double Jfec = kfe*(CCA + CCDCA + CDCA + CLCA) + klcaout*CLCA
              + kseq3out*BND;
dxdt_FEC = Jfec;
double Jcolabs = Jabs_dca + Jabs_cdca;

// ---- secretory drive (mM CDCA-equivalent, FREE bile acid only) ------
double Ds = (wCA*CCA + wCDCA*CCDCA + wDCA*CDCA + wLCA*CLCA)/V_col;
double G  = pow(Ds,hG)/(pow(EC50_G,hG) + pow(Ds,hG));

double OND = 1.0 - fOND*ONDc/(IC50_ond + ONDc);
dxdt_HT5 = (aHT*G*OND - HT5)/tauHT;
double LOP = 1.0/(1.0 + LOPe/IC50_lop);
double TRANSss = (1.0 + aTR*(0.55*G + 0.45*HT5)) * LOP;
if (TRANSclamp > 0.0) dxdt_TRANS = (TRANSclamp - TRANS)/0.25;  // loop OPENED
else                  dxdt_TRANS = (TRANSss - TRANS)/tauTR;
dxdt_PERM  = (1.0 + aPERM*G - PERM)/tauPERM;

double Jsecw  = Jsec*(0.6*G + 0.4*HT5);
double Jabsw  = kwabs*WCOL/TRANS;
double Jstool = kstool*TRANS*WCOL;
dxdt_WCOL  = Qin + Jsecw - Jabsw - Jstool;
dxdt_STOOL = Jstool;

dxdt_BSHx = kbsh*(1.0 - BSHx);
dxdt_BAI  = kbai_t*(1.0 - BAI);

// ===================== bile acid mass balance =========================
double Jpor_liv = kpor*POR*(1.0 - fspill);
double Jpor_sys = kpor*POR*fspill;
double Jsys_liv = ksys*SYS;
dxdt_LIV = S + Jpor_liv + Jsys_liv - Jbsep;
dxdt_GB  = to_gb - Jeject;
dxdt_DUO = Jduo_in - Jd2j;
dxdt_JEJ = Jd2j - Jj2i - Jpass;
dxdt_POR = Jost + Jpass + Jcolabs - kpor*POR;
dxdt_SYS = Jpor_sys - Jsys_liv;

// ---- circulating-pool composition (fractions, no mass leak) ---------
double Pool = LIV + GB + DUO + JEJ + ILE1 + ILE2 + ILE3 + ENT + POR + SYS;
double Jin  = S + Jcolabs;
dxdt_fCAp  = 0.0;
dxdt_fDCAp = 0.0;
if (Pool > 1e-6 && Jin > 1e-12) {
  dxdt_fCAp  = Jin*(S*fCAnew/Jin - fCAp)/Pool;
  dxdt_fDCAp = Jin*(Jabs_dca/Jin - fDCAp)/Pool;
}

// ===================== fat absorption / CMC ============================
double C_duo   = DUO/V_duo;                            // mM
double fat_eff = fatmax*pow(C_duo,hCMC)/(pow(CMC,hCMC) + pow(C_duo,hCMC));
double fat_app = fatIn*M/MEALAUC*(1.0 - fat_eff);      // g/h into stool
dxdt_FECFAT = fat_app;
// FATF low-pass filters the meal-pulsed appearance rate into a smooth
// g/day figure.  Downstream endpoints (oxalate, vitamin D) must see the
// DAILY loss, not the between-meal instantaneous value, which is
// meaningless because nothing is being eaten then.
dxdt_FATF = (24.0*fat_app - FATF)/24.0;
double fatrate = FATF;                                 // g/day

// ===================== systemic endpoints ===============================
dxdt_HCHOL = (1.0 - HCHOL)/48.0;
double Sx = (S > 1e-6) ? S : 1e-6;
double LDLss = LDL0*pow(Vs/Sx, expLDL);
dxdt_LDLC = (LDLss - LDLC)/tauLDL;
double exfat = fatrate - 3.0; if (exfat < 0.0) exfat = 0.0;
double UOXss = UOX0 + aUOXfat*exfat + aUOXperm*(PERM - 1.0);
dxdt_UOX  = (UOXss - UOX)/tauUOX;
double fabs_day = (fatIn - fatrate)/fatIn/fatmax;
if (fabs_day > 1.0) fabs_day = 1.0;
if (fabs_day < 0.0) fabs_day = 0.0;
dxdt_VITD = (VITD0*fabs_day - VITD)/tauVITD;

// ===================== SeHCAT tracer =====================================
dxdt_SEH = -(Jfec/(Pool > 1e-9 ? Pool : 1e-9))*SEH;

// ===================== drug PK ============================================
dxdt_SEQ1 = -kseq12*SEQ1;
dxdt_SEQ2 =  kseq12*SEQ1 - kseq23*SEQ2;
dxdt_SEQ3 =  kseq23*SEQ2 - kseq3out*SEQ3;
double oca_asbt = (1.0 - foca_pass)*phi*ASBTx*ENTMASS;
double kaO = ka_oca*(foca_pass + oca_asbt);
dxdt_OCAg = -(kaO + kOCAgo)*OCAg;
dxdt_OCAe =  kaO*OCAg - kOCAep*OCAe;
dxdt_OCAp =  kOCAep*OCAe - kOCAout*OCAp;
dxdt_ELOg = -kelo_out*ELOg;
dxdt_LOPg = -ka_lop*LOPg;
dxdt_LOPe =  ka_lop*LOPg - ke_lop*LOPe;
dxdt_RIFg = -ke_rif*RIFg;
dxdt_ONDc = -ke_ond*ONDc;
dxdt_TROg = -ka_tro*TROg;
dxdt_TROp =  ka_tro*TROg - ke_tro*TROp;

double tgt = 1.0 + (ASBTmax - 1.0)*(1.0 - FXRi);
dxdt_ASBTx = kASBT*(tgt - ASBTx);

$TABLE
double POOLtot = LIV + GB + DUO + JEJ + ILE1 + ILE2 + ILE3 + ENT + POR + SYS;
double SYNTH   = Vs * CYP7A1 * HCHOL * 24.0;                   // umol/day
double CDUO    = DUO/V_duo;                                    // mM
double DS = (wCA*CCA + wCDCA*CCDCA + wDCA*CDCA + wLCA*CLCA)/V_col;
double GDRIVE = pow(DS,hG)/(pow(EC50_G,hG) + pow(DS,hG));
double SERUMBA = SYS/12.0;                                     // uM (12 L)
double FATEFF = fatmax*pow(CDUO,hCMC)/(pow(CMC,hCMC) + pow(CDUO,hCMC));
double POOL_G = POOLtot*500.0/1e6;                             // g

$CAPTURE @annotated
POOLtot : total circulating bile acid pool (umol)
POOL_G  : total circulating bile acid pool (g)
SYNTH   : instantaneous hepatic synthesis rate (umol/day)
DS      : colonic secretory drive (mM CDCA-equivalent)
GDRIVE  : fractional activation of the colonic secretory switch
CDUO    : duodenal bile acid concentration (mM)
FATEFF  : fractional fat absorption
SERUMBA : systemic serum bile acid (uM)
'

mod <- mcode_cache("bam_qsp", code)

# =====================================================================
#  DOSING HELPERS
# =====================================================================
#  All regimens are given in the units of the corresponding compartment.
#  Molecular weights: obeticholic acid 420.6, elobixibat 514.8,
#  tropifexor 604.7 g/mol.
# =====================================================================

ev_seq <- function(mg = 1875, n_per_day = 2, days = 56) {
  h <- switch(as.character(n_per_day),
              "1" = 8, "2" = c(8, 19), "3" = c(8, 13, 19), c(8, 19))
  e <- NULL
  for (tt in h) {
    ei <- ev(amt = mg, cmt = "SEQ1", time = tt, ii = 24, addl = days - 1)
    e <- if (is.null(e)) ei else e + ei
  }
  e
}
ev_oca <- function(mg = 25, days = 56)
  ev(amt = mg / 420.6 * 1000, cmt = "OCAg", time = 8, ii = 24, addl = days - 1)
ev_tro <- function(mg = 0.090, days = 56)
  ev(amt = mg / 604.7 * 1000, cmt = "TROg", time = 8, ii = 24, addl = days - 1)
ev_elo <- function(mg = 10, days = 56)
  ev(amt = mg / 514.8 * 1000, cmt = "ELOg", time = 7.5, ii = 24, addl = days - 1)
ev_lop <- function(mg = 2, days = 56)
  ev(amt = mg, cmt = "LOPg", time = 7.5, ii = 24, addl = days - 1) +
  ev(amt = mg, cmt = "LOPg", time = 18.5, ii = 24, addl = days - 1)
ev_ond <- function(mg = 4, days = 56)
  ev(amt = mg, cmt = "ONDc", time = 7.5, ii = 24, addl = days - 1) +
  ev(amt = mg, cmt = "ONDc", time = 18.5, ii = 24, addl = days - 1)
ev_rif <- function(mg = 550, days = 56)
  ev(amt = mg, cmt = "RIFg", time = 8, ii = 24, addl = days - 1) +
  ev(amt = mg, cmt = "RIFg", time = 13, ii = 24, addl = days - 1) +
  ev(amt = mg, cmt = "RIFg", time = 19, ii = 24, addl = days - 1)

# ---------------------------------------------------------------------
#  phi from resection length: ASBT density falls exponentially in the
#  proximal direction with a decay length of 60 cm (anatomical input,
#  NOT fitted to resection outcome data).
# ---------------------------------------------------------------------
phi_from_resection <- function(L_cm, lambda = 60) exp(-L_cm / lambda)

# ---------------------------------------------------------------------
#  Day-averaged readouts from the last simulated 24 h.
# ---------------------------------------------------------------------
bam_readout <- function(out) {
  d <- as.data.frame(out)
  t_end <- max(d$time)
  d <- dplyr::filter(d, time >= t_end - 24)
  span <- max(d$time) - min(d$time)
  tint <- function(x) sum(diff(d$time) * (head(x, -1) + tail(x, -1)) / 2) / span
  fec_day  <- (tail(d$FEC, 1)    - head(d$FEC, 1))    / span * 24
  stool    <- (tail(d$STOOL, 1)  - head(d$STOOL, 1))  / span * 24
  fat      <- (tail(d$FECFAT, 1) - head(d$FECFAT, 1)) / span * 24
  pool     <- tint(d$POOLtot)
  kturn    <- fec_day / pool
  fasting  <- dplyr::filter(d, (time %% 24) >= 5, (time %% 24) <= 7.5)
  data.frame(
    pool_g      = pool * 500 / 1e6,
    synth_umold = tint(d$SYNTH),
    fecBA_umold = fec_day,
    C4          = mean(fasting$C4),
    FGF19       = mean(fasting$FGF19),
    SeHCAT_pct  = 100 * exp(-7 * kturn),
    kturn_day   = kturn,
    Ds          = tint(d$DS),
    stool_mL    = stool,
    BM_per_day  = 1.10 * (stool / 130)^1.05,
    bristol     = pmin(7, pmax(1, 3.6 + 3 * (1 - exp(-1.05 * (stool / 130 - 1))))),
    fecfat_g    = fat,
    Cduo_peak   = max(d$CDUO),
    LDLC        = tint(d$LDLC),
    UOX         = tint(d$UOX),
    TRANS       = tint(d$TRANS)
  )
}

sim_bam <- function(days = 56, events = NULL, ...) {
  p <- list(...)
  m <- if (length(p)) do.call(mrgsolve::param, c(list(mod), p)) else mod
  e <- if (is.null(events)) ev(amt = 0, cmt = "SEQ1", time = 0) else events
  mrgsolve::mrgsim(m, events = e, end = days * 24, delta = 0.25,
                   hmax = 0.35, rtol = 1e-7, atol = 1e-9)
}

# =====================================================================
#  SCENARIO LIBRARY  —  12 scenarios
#  (numbers reproduced by bam_verify_python.py; see
#   bam_verification_output.txt)
# =====================================================================
bam_scenarios <- function() {
  S <- list(
    # ---- 1. healthy control --------------------------------------
    "01 healthy control" =
      list(par = list(), ev = NULL),

    # ---- 2. THE CENTRAL EXPERIMENT: loop live vs loop frozen -----
    #  Identical phi lesion; only the feedback switch differs.
    "02a severe malabsorption, feedback LIVE" =
      list(par = list(phi = 0.05, FEEDBACK = 1), ev = NULL),
    "02b severe malabsorption, feedback FROZEN" =
      list(par = list(phi = 0.05, FEEDBACK = 0), ev = NULL),

    # ---- 3. the four aetiological types ---------------------------
    "03 type 1 - 60 cm ileal resection" =
      list(par = list(phi = phi_from_resection(60)), ev = NULL),
    "04 type 1 - ileal Crohn disease" =
      list(par = list(phi = 0.25, ENTMASS = 0.70), ev = NULL),
    "05 type 2 - idiopathic (primary) BAD" =
      list(par = list(phi = 0.95, kappa = 0.40), ev = NULL),
    "06 type 3 - post-cholecystectomy" =
      list(par = list(GBfun = 0, phi = 0.75), ev = NULL),
    "07 type 4 - metformin-associated" =
      list(par = list(phi = 0.55, kbai = 0.16 * 1.4), ev = NULL),

    # ---- 4. therapy, all in the same reference patient ------------
    "08 moderate BAM - untreated (phi 0.20)" =
      list(par = list(phi = 0.20), ev = NULL),
    "09 moderate BAM + colesevelam 1.875 g BID" =
      list(par = list(phi = 0.20), ev = ev_seq(1875, 2)),
    "10 moderate BAM + obeticholic acid 25 mg" =
      list(par = list(phi = 0.20), ev = ev_oca(25)),
    "11 moderate BAM + colesevelam + OCA" =
      list(par = list(phi = 0.20), ev = ev_seq(1875, 2) + ev_oca(25)),
    "12 moderate BAM + loperamide 2 mg BID" =
      list(par = list(phi = 0.20), ev = ev_lop(2)),

    # ---- 5. iatrogenic BAM: the mirror image ----------------------
    "13 normal gut + elobixibat 10 mg (constipation Rx)" =
      list(par = list(), ev = ev_elo(10)),

    # ---- 6. beyond the compensation ceiling -----------------------
    "14 120 cm resection - steatorrhoea territory" =
      list(par = list(phi = phi_from_resection(120)), ev = NULL)
  )
  res <- lapply(names(S), function(nm) {
    s <- S[[nm]]
    out <- do.call(sim_bam, c(list(days = 56, events = s$ev), s$par))
    cbind(scenario = nm, bam_readout(out))
  })
  dplyr::bind_rows(res)
}

# =====================================================================
#  PARAMETER PROVENANCE / CALIBRATION NOTES
# ---------------------------------------------------------------------
#  Everything below is either (i) a directly measured physiological
#  quantity, (ii) fixed by a single published experiment, or
#  (iii) calibrated to the healthy steady state.  Nothing was fitted to
#  bile-acid-diarrhoea outcome data — the disease behaviour is a
#  prediction, not a fit.
#
#  MEASURED / FIXED BY THE LITERATURE
#    total pool 3 g, 4-6 enterohepatic cycles/day, duodenal delivery
#      20-40 mmol/day, synthesis 0.35 g/day .......... PMID 6682120, 7017051
#    per-pass ileal conservation f ~ 0.97 ............ PMID 5029077, 19498215
#    ASBT Km ~ 6 mM luminal, distally weighted ....... PMID 19498215, 9109432
#    colonic secretory threshold 3 mM di-OH BA ....... PMID 4938344  <= EC50_G
#    potency order CA << CDCA < DCA .................. PMID 4938344, 23041323
#    fasting FGF19 ~ 230 pg/mL, t1/2 ~ 30 min ........ PMID 17116003, 23238290
#    fasting C4 ~ 18 ng/mL (BAD cut-off ~ 48) ........ PMID 23238290, 28691284
#    SeHCAT grading 15/10/5 % ........................ PMID 24351663
#    OCA is ~100x more potent than CDCA at FXR ....... PMID 12166927
#    CYP7A1 protein t1/2 ~ 4 h, strong circadian ..... PMID 26579436, 17116003
#    ASBT decay length 60 cm (anatomy, NOT fitted) ... PMID 5029077, 30116880
#
#  CALIBRATED TO THE HEALTHY STEADY STATE (5 numbers only)
#    Vmax_asbt  -> total pool 6000 umol
#    kf19       -> fasting FGF19 230 pg/mL
#    IC50f      -> CYP7A1 = 1.0 (i.e. synthesis 700 umol/day)
#    kc4        -> fasting C4 18 ng/mL
#    kfec       -> colonic secretory drive 1.2 mM-eq
#    kwabs      -> stool water 130 mL/day
#
#  NOT CALIBRATED AT ALL (pure predictions)
#    every disease scenario, every drug effect, the 100 cm rule,
#    the SeHCAT/C4 dissociation, the sequestrant C4 escape, and the
#    superadditivity of sequestrant + FXR agonist.
# =====================================================================

if (sys.nframe() == 0L) {
  print(as.data.frame(bam_scenarios()), digits = 3)
}
