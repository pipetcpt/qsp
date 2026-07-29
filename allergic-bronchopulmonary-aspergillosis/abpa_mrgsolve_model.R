## ===========================================================================
##  Allergic Bronchopulmonary Aspergillosis (ABPA) — mrgsolve QSP model
##  ---------------------------------------------------------------------------
##  43 ODE compartments.  The $PARAM block is the SAME parameter block as
##  abpa_reference_implementation.py — all 163 shared parameters agree to
##  machine precision (the R file adds only the two counterfactual switches
##  DDI_OFF and AF_OFF), so any number in README.md can be reproduced from
##  either file.  The Python reference carries 42 ODE states rather than 43
##  because it integrates the omalizumab SC depot analytically between events
##  instead of as a compartment; the two are equivalent, the depot having no
##  feedback.
##
##  THE STRUCTURAL CLAIM THE MODEL EXISTS TO TEST
##  ---------------------------------------------
##  The Aspergillus population is partitioned into two compartments:
##       FLUM   luminal / mucosa-adherent  -> sees unbound plasma drug
##       FPLG   embedded in a mucus plug   -> sees f_pen x unbound plasma drug
##  coupled by entrapment k_in*PLUG*FLUM and release k_out(PLUG)*FPLG.
##  Linearising at low burden gives a 2x2 Jacobian whose stability condition
##  contains the dose only through f_pen in the sanctuary row, so there is a
##  plug-clearance rate below which NO antifungal exposure clears the organism.
##  k_out is not an antifungal parameter — it belongs to the type-2 mucin axis.
##  Therefore the azole and the biologic are not two ways of doing one thing:
##  one lowers a threshold the other has to clear.
##
##  Two further structural traps are implemented explicitly:
##   * TOTAL IgE RISES on omalizumab while free IgE falls, because the
##     omalizumab:IgE complex is cleared more slowly than free IgE.  The ABPA
##     response criterion ("total IgE down 35-50%") therefore moves OPPOSITE to
##     the pharmacology.  Outputs IGE_TOTAL and IGE_FREE are both reported and
##     must never be plotted on one axis without saying which is which.
##   * ITRACONAZOLE INHIBITS CYP3A4 at two sequential sites (gut wall and
##     liver), so part of its apparent corticosteroid-sparing effect is a drug
##     interaction whose size depends on WHICH steroid the protocol specified.
##
##  CALIBRATION ANCHORS  (see abpa_references.md for the sources)
##  -------------------------------------------------------------
##   untreated endotype   total IgE 2000 IU/mL, blood eos 800/uL, FEV1 ~68% pred,
##                        mucus-plug score ~6/18, bronchiectasis 3/18
##   itraconazole PK      200 mg BID -> parent ~1.28 mg/L, OH-itra ~2.20 mg/L
##                        (OH/parent 1.7), nonlinear via CYP3A4 autoinhibition
##   voriconazole PK      200 mg BID -> ~2.2 mg/L (CYP2C19 NM), ~4 mg/L (PM)
##   CYP3A4 inhibition    fractional activity I ~ 0.29 at itraconazole steady
##                        state; steroid AUC fold-changes 1.12 / 2.57 / 3.92
##                        for prednisolone / methylprednisolone / budesonide
##                        against observed 1.24 / 2.6 / 4.2
##   omalizumab           free IgE -96%, total IgE x3.1, complex t1/2 8 d
##   corticosteroid       ISHAM medium-dose taper; cortisol 14 -> ~9 ug/dL on
##                        prednisolone 10 mg/d maintenance
##
##  USAGE
##  -----
##    library(mrgsolve); library(dplyr)
##    mod <- mread("abpa_mrgsolve_model", "path/to/dir")
##    out <- mod %>% ev(abpa_scenario("itraconazole")) %>% mrgsim(end = 364)
##    plot(out, IGE_TOTAL + IGE_FREE + PLUG_SCORE + FPLG ~ time)
##
##  Scenario builders are at the bottom of this file (abpa_scenario()).
## ===========================================================================

$PROB
# ABPA QSP — sanctuary partition, IgE TMDD, CYP3A4 interaction

$PARAM @annotated
// ---- antigen / immune drive ------------------------------------------------
kag     : 1.00    : antigen release per unit fungal burden (1/d)
kagd    : 0.49    : antigen clearance (1/d)
ag_plug : 0.80    : relative antigen release from plug-embedded hyphae (-)
kin2    : 0.2937  : Th2 recruitment rate (1/d)
Kag2    : 0.35    : antigen EC50 for Th2 drive (-)
kout2   : 0.12    : Th2 loss rate (1/d)
th2_base: 0.45    : antigen-independent memory Th2 drive (-)
s13     : 0.90    : IL-13 secretion per Th2 (1/d)
kd13    : 1.20    : IL-13 clearance (1/d)
s5      : 0.90    : IL-5 secretion per Th2 (1/d)
kd5     : 1.20    : IL-5 clearance (1/d)

// ---- plasma cell / IgE -----------------------------------------------------
spc     : 0.03718 : antigen-driven plasma-cell recruitment (1/d)
pc0     : 0.0538  : long-lived bone-marrow plasma-cell floor (1/d)
Kagpc   : 0.40    : antigen EC50 for plasma-cell drive (-)
a4      : 0.85    : IL-4/IL-13 amplification of class switching (-)
kpc     : 0.020   : plasma-cell loss rate (1/d)
kige    : 1.30    : IgE secretion per plasma-cell unit (nM/d)
kdegE   : 0.277   : free IgE elimination, t1/2 2.5 d (1/d)
kelCX   : 0.0866  : omalizumab:IgE complex elimination, t1/2 8 d (1/d)
Kd_oma  : 1.0     : omalizumab:IgE dissociation constant (nM)
eta_oma : 0.50    : IgE neutralised per omalizumab, 2:1 average (-)
ka_oma  : 0.55    : omalizumab SC absorption (1/d)
F_oma   : 0.62    : omalizumab SC bioavailability (-)
kel_oma : 0.0267  : omalizumab elimination, t1/2 26 d (1/d)
Kd_fcer : 0.10    : FcepsilonRI affinity for IgE (nM)
ksyn_fcer : 0.050 : FcepsilonRI synthesis (1/d)
kdeg_fcer : 0.050 : FcepsilonRI turnover (1/d)
b_fcer  : 6.0     : IgE-driven FcepsilonRI up-regulation (-)
Kup_fcer: 0.35    : free IgE EC50 for receptor up-regulation (nM)

// ---- eosinophils -----------------------------------------------------------
kin_e   : 120.0   : marrow eosinophil output (cells/uL/d)
a5      : 1.20    : IL-5 amplification of eosinophilopoiesis (-)
kout_e  : 0.42    : blood eosinophil loss (1/d)
cs_eos  : 1.60    : corticosteroid enhancement of eosinophil loss (-)
kADCC   : 0.90    : benralizumab ADCC-mediated depletion (1/d)
EC50_benr : 0.30  : benralizumab EC50 for ADCC (mg/L)
ktr_e   : 0.170   : blood-to-airway eosinophil egress (1/d per 1000/uL)
aeot    : 0.80    : IL-13/eotaxin amplification of egress (-)
kout_ea : 0.30    : airway eosinophil loss (1/d)
cs_apop : 1.30    : corticosteroid enhancement of eosinophil apoptosis (-)
dup_egress : 0.55 : fraction of egress blocked by IL-4Ra blockade (-)
sepx    : 1.00    : EPX release per airway eosinophil (1/d)
kepx    : 0.625   : EPX clearance (1/d)

// ---- mucus plug ------------------------------------------------------------
smuc    : 0.352   : IL-13-driven mucin secretion (1/d)
Kmuc    : 1.50    : IL-13 EC50 for mucin secretion (-)
muc0    : 0.030   : IL-13-independent mucin secretion (1/d)
smuc_ige: 0.060   : mast-cell/IgE-driven mucin secretion (1/d)
cs_muc  : 0.55    : corticosteroid suppression of mucin secretion (-)
kout0   : 0.550   : baseline plug clearance (1/d)
cs_plug : 1.90    : corticosteroid enhancement of plug clearance (-)
g_epx   : 0.55    : EPX cross-linking penalty on plug clearance (-)
g_br    : 0.60    : bronchiectasis penalty on plug clearance (-)

// ---- THE SANCTUARY PARTITION ----------------------------------------------
gl      : 0.85    : luminal fungal growth rate (1/d)
gp      : 0.350   : INTRA-PLUG fungal growth rate (1/d) -- never measured in a
                  : human airway; A2(e) shows the model turns on sign(gp-kout)
kin_f   : 0.063   : entrapment rate per unit plug (1/d)
seed    : 0.0020  : inhaled conidial deposition (1/d) -- makes true
                  : sterilisation impossible by construction
k_host  : 1.265   : host (phagocyte) killing of luminal fungus (1/d)
cs_imm  : 0.55    : corticosteroid suppression of host killing (-)
f_pen   : 0.10    : azole penetration into a plug, unbound-fraction-equivalent
f_pen_amb : 0.03  : liposomal amphotericin penetration into a plug (-)
Emax_af : 1.55    : maximum azole kill rate = the WHOLE CLASS ceiling (1/d)
EC50_af : 0.0038  : unbound itraconazole-equivalent EC50 (mg/L)
Emax_amb: 1.20    : maximum amphotericin kill rate (1/d)
EC50_amb: 0.90    : amphotericin EC50 in epithelial lining fluid (mg/L)

// ---- itraconazole PK (nonlinear, autoinhibiting) --------------------------
ka_itra : 6.0     : absorption (1/d)
F_itra  : 0.55    : bioavailability, capsule with food (-)
V_itra  : 700.0   : volume of distribution (L)
kel_itra: 0.79    : elimination at low concentration (1/d)
Ki_auto : 0.55    : autoinhibition constant (mg/L)
fm_oh   : 1.85    : fraction of cleared parent appearing as OH-itraconazole (-)
kel_oh  : 0.26    : OH-itraconazole elimination (1/d)
fu_itra : 0.0020  : itraconazole unbound fraction (99.8% bound)
fu_oh   : 0.0016  : OH-itraconazole unbound fraction (-)
pot_oh  : 1.00    : OH-itraconazole antifungal potency vs parent (-)
Ki3A4_itra : 0.0015 : unbound CYP3A4 Ki, itraconazole (mg/L)
Ki3A4_oh   : 0.0050 : unbound CYP3A4 Ki, OH-itraconazole (mg/L)

// ---- voriconazole PK (mixed saturable + linear) ---------------------------
ka_vori : 8.0     : absorption (1/d)
F_vori  : 0.90    : bioavailability (-)
V_vori  : 280.0   : volume of distribution (L)
Vmax_vori : 400.0 : saturable clearance capacity (mg/d)
Km_vori : 1.60    : Michaelis constant (mg/L)
CLlin_vori : 60.0 : non-saturable clearance (L/d)
cyp2c19 : 1.0     : CYP2C19 activity scalar; 1.0 NM, 0.40 PM, 1.6 UM (-)
fu_vori : 0.42    : unbound fraction (-)
pot_vori: 0.0066  : potency vs unbound itraconazole (-)
Ki3A4_vori : 0.35 : unbound CYP3A4 Ki, voriconazole (mg/L)
kel_amb : 0.14    : nebulised liposomal amphotericin airway decay (1/d)

// ---- corticosteroid PK, SEQUENTIAL gut + hepatic CYP3A4 -------------------
// A single-site model cannot reproduce the observed itraconazole-budesonide
// 4.2x rise: the ceiling of any single-site model is 1/I = 3.4 at this I.
ka_pred : 8.0     : prednisolone absorption (1/d)
F_pred  : 0.80    : prednisolone bioavailability (-)
V_pred  : 45.0    : prednisolone volume (L)
kel_pred: 5.55    : prednisolone elimination, t1/2 3 h (1/d)
s3A4g_pred : 0.00 : prednisolone gut-wall CYP3A4 share (-)
s3A4h_pred : 0.15 : prednisolone hepatic CYP3A4 share (-)
ka_mpred: 8.0     : methylprednisolone absorption (1/d)
F_mpred : 0.85    : methylprednisolone bioavailability (-)
V_mpred : 90.0    : methylprednisolone volume (L)
kel_mpred : 4.60  : methylprednisolone elimination (1/d)
s3A4g_mpred : 0.15 : methylprednisolone gut-wall CYP3A4 share (-)
s3A4h_mpred : 0.80 : methylprednisolone hepatic CYP3A4 share (fitted to 2.6x)
ka_bud  : 8.0     : swallowed budesonide absorption (1/d)
F_bud   : 0.11    : budesonide systemic bioavailability (-)
V_bud   : 210.0   : budesonide volume (L)
kel_bud : 9.20    : budesonide elimination (1/d)
s3A4g_bud : 0.70  : budesonide gut-wall CYP3A4 share (fitted with the next)
s3A4h_bud : 0.70  : budesonide hepatic CYP3A4 share (fitted to 4.2x)
fu_pred : 0.25    : prednisolone unbound fraction (-)
fu_mpred: 0.23    : methylprednisolone unbound fraction (-)
fu_bud  : 0.12    : budesonide unbound fraction (-)
pot_mpred : 1.25  : methylprednisolone GR potency vs prednisolone (-)
pot_bud : 40.0    : budesonide GR potency vs prednisolone (-)
Emax_cs : 1.00    : maximum GR effect (-)
EC50_cs : 0.0090  : unbound prednisolone-equivalent EC50 (mg/L)
w_t2    : 0.85    : GR weight on Th2 recruitment (-)
w_cyt   : 0.80    : GR weight on cytokine secretion (-)
w_csr   : 0.60    : GR weight on class switching (-)

// ---- HPA axis and toxicity ------------------------------------------------
kin_c   : 2.10    : cortisol production (ug/dL/d)
kout_c  : 0.150   : cortisol elimination (1/d)
Imax_c  : 0.96    : maximum HPA suppression (-)
kbmd    : 0.00042 : GR-driven bone loss (1/d)
krep_bmd: 0.00020 : bone recovery toward baseline (1/d)
tau_a1c : 45.0    : HbA1c turnover time (d)
a1c0    : 5.40    : baseline HbA1c (%)
a1c_gain: 2.00    : maximum steroid-driven HbA1c rise (%)

// ---- biologics PK/PD ------------------------------------------------------
ka_mep  : 0.28    : mepolizumab SC absorption (1/d)
F_mep   : 0.80    : mepolizumab bioavailability (-)
V_mep   : 3.60    : mepolizumab volume (L)
kel_mep : 0.0347  : mepolizumab elimination, t1/2 20 d (1/d)
Emax_mep: 0.90    : maximum IL-5 neutralisation (-)
EC50_mep: 0.55    : mepolizumab EC50 (mg/L)
ka_ben  : 0.25    : benralizumab SC absorption (1/d)
F_ben   : 0.58    : benralizumab bioavailability (-)
V_ben   : 3.10    : benralizumab volume (L)
kel_ben : 0.0462  : benralizumab elimination, t1/2 15 d (1/d)
Emax_ben_block : 0.85 : maximum IL-5Ra blockade (-)
EC50_ben: 0.30    : benralizumab EC50 (mg/L)
ka_dup  : 0.26    : dupilumab SC absorption (1/d)
F_dup   : 0.64    : dupilumab bioavailability (-)
V_dup   : 4.80    : dupilumab volume (L)
kel_dup : 0.0630  : dupilumab elimination (1/d)
Emax_dup: 0.82    : maximum IL-4Ra blockade (-)
EC50_dup: 1.90    : dupilumab EC50 (mg/L)
ka_tez  : 0.20    : tezepelumab SC absorption (1/d)
F_tez   : 0.77    : tezepelumab bioavailability (-)
V_tez   : 3.90    : tezepelumab volume (L)
kel_tez : 0.0315  : tezepelumab elimination (1/d)
Emax_tez: 0.50    : maximum TSLP-driven alarmin blockade (-)
EC50_tez: 2.00    : tezepelumab EC50 (mg/L)

// ---- structure and physiology --------------------------------------------
kb1     : 0.00060 : plug-driven bronchiectasis accrual (score/d)
Pthr    : 0.35    : plug threshold below which no structural damage accrues (-)
kb2     : 0.340   : exacerbation-driven bronchiectasis accrual (score/event)
BRONmax : 18.0    : maximum bronchiectasis score (-)
FEV1_0  : 100.0   : FEV1 in the absence of all disease terms (% predicted)
cp_fev  : 0.34    : plug contribution to airflow obstruction (-)
Kp_fev  : 0.85    : plug EC50 for obstruction (-)
cb_fev  : 0.42    : bronchiectasis contribution to fixed obstruction (-)
ca_fev  : 0.14    : hyper-responsiveness contribution (-)
Ka_fev  : 1.00    : IL-13 EC50 for hyper-responsiveness (-)
w_ahr_il13 : 0.60 : IL-13 weight in the AHR term (-)
w_ahr_ige  : 0.40 : mast-cell/IgE weight in the AHR term (-)
h0      : 0.00021 : baseline exacerbation hazard (1/d)
b1_eos  : 0.85    : airway eosinophil weight in the hazard (-)
b2_plug : 1.15    : mucus plug weight in the hazard (-)
b3_fun  : 0.55    : fungal burden weight in the hazard (-)
b4_ige  : 0.55    : IgE effector weight in the hazard (-)

// ---- switches used for counterfactual arms -------------------------------
DDI_OFF : 0       : set to 1 to abolish the CYP3A4 interaction (counterfactual)
AF_OFF  : 0       : set to 1 to abolish antifungal killing (counterfactual)

$CMT @annotated
AG    : airway antigen / allergen load (-)
TH2   : Af-specific Th2 pool (-)
IL13  : IL-13 (-)
IL5   : IL-5 (-)
PC    : IgE-secreting plasma-cell pool (-)
TI    : TOTAL IgE, free + complexed (nM)
TO    : TOTAL omalizumab, free + complexed (nM)
OMAD  : omalizumab SC depot (mg)
FCER  : FcepsilonRI density on mast cells / basophils (relative)
EOSB  : blood eosinophils (cells/uL)
EOSA  : airway eosinophils (-)
EPX   : eosinophil peroxidase cross-linking capacity (-)
PLUG  : mucus plug mass (-)
FLUM  : luminal fungal burden -- DRUG-ACCESSIBLE (-)
FPLG  : plug-embedded fungal burden -- DRUG-SANCTUARY (-)
AITR  : itraconazole gut depot (mg)
ITRA  : itraconazole plasma (mg/L)
OHIT  : OH-itraconazole plasma (mg/L)
AVOR  : voriconazole gut depot (mg)
VORI  : voriconazole plasma (mg/L)
AMB   : nebulised amphotericin in epithelial lining fluid (mg/L)
APRD  : prednisolone gut depot (mg)
PRED  : prednisolone plasma (mg/L)
AMPD  : methylprednisolone gut depot (mg)
MPRD  : methylprednisolone plasma (mg/L)
ABUD  : swallowed budesonide depot (mg)
BUD   : budesonide plasma (mg/L)
CORT  : endogenous cortisol (ug/dL)
BMD   : lumbar bone mineral density (fraction of baseline)
HBA1C : HbA1c (%)
CUMO  : cumulative prednisolone-equivalent dose (mg)
MEPD  : mepolizumab SC depot (mg)
MEPO  : mepolizumab plasma (mg/L)
BEND  : benralizumab SC depot (mg)
BENR  : benralizumab plasma (mg/L)
DUPD  : dupilumab SC depot (mg)
DUPI  : dupilumab plasma (mg/L)
TEZD  : tezepelumab SC depot (mg)
TEZE  : tezepelumab plasma (mg/L)
BRON  : bronchiectasis extent, 0-18 -- AN INTEGRATOR, never falls
CHAZ  : cumulative exacerbation hazard (expected events)
AUCP  : integral of unbound prednisolone-equivalent (mg/L.d)
AUCCS : integral of the GR effect actually delivered (d)

$GLOBAL
#define IGE_IU_TO_NM  (0.01263)     // 1 IU/mL = 2.4 ng/mL, MW 190 kDa
#define IGE_NM_TO_IU  (79.1765)
#define IGE_NM_TO_NGML (190.024)    // nM -> ng/mL
#define OMA_MG_TO_NM  (1.22930)     // 1 mg absorbed -> nM in 5.46 L

// rapid-equilibrium 1:1 binding; returns the complex concentration
double bind11(double B, double L, double Kd) {
  if (B <= 0.0 || L <= 0.0) return 0.0;
  double s = B + L + Kd;
  double disc = s*s - 4.0*B*L;
  if (disc < 0.0) disc = 0.0;
  return 0.5*(s - sqrt(disc));
}
double pos(double x) { return x > 0.0 ? x : 0.0; }

$MAIN
// steady-state-consistent initial conditions (the untreated ABPA endotype)
AG_0    = 0.90;   TH2_0  = 2.00;  IL13_0 = 1.50;  IL5_0  = 1.50;
PC_0    = 5.38;   TI_0   = 2000.0*IGE_IU_TO_NM;   TO_0   = 0.0;
FCER_0  = 1.00;   EOSB_0 = 800.0; EOSA_0 = 1.00;  EPX_0  = 1.60;
PLUG_0  = 1.00;   FLUM_0 = 0.25;  FPLG_0 = 0.35;
CORT_0  = 14.0;   BMD_0  = 1.00;  HBA1C_0 = 5.40; BRON_0 = 3.00;

$ODE
// ===========================================================================
//  algebra shared by the whole system
// ===========================================================================
double ige_c   = bind11(eta_oma*pos(TO), pos(TI), Kd_oma);   // complexed IgE
double ige_f   = pos(TI) - ige_c;                            // FREE IgE (nM)
double oma_f   = pos(TO) - ige_c/eta_oma;                    // free omalizumab

// unbound prednisolone-equivalent, and the GR effect it delivers
double peq = fu_pred*pos(PRED)
           + fu_mpred*pos(MPRD)*pot_mpred
           + fu_bud*pos(BUD)*pot_bud;
double CS  = (peq > 0.0) ? Emax_cs*peq/(peq + EC50_cs) : 0.0;

// fractional CYP3A4 activity remaining (the drug interaction, in one line)
double I3A4 = 1.0/(1.0 + fu_itra*pos(ITRA)/Ki3A4_itra
                       + fu_oh  *pos(OHIT)/Ki3A4_oh
                       + fu_vori*pos(VORI)/Ki3A4_vori);
if (DDI_OFF > 0.5) I3A4 = 1.0;

// biologic occupancies
double E_dup = Emax_dup       *pos(DUPI)/(pos(DUPI) + EC50_dup);
double E_mep = Emax_mep       *pos(MEPO)/(pos(MEPO) + EC50_mep);
double E_ben = Emax_ben_block *pos(BENR)/(pos(BENR) + EC50_ben);
double E_tez = Emax_tez       *pos(TEZE)/(pos(TEZE) + EC50_tez);
double IL13e = pos(IL13)*(1.0 - E_dup);
double IL5e  = pos(IL5) *(1.0 - (E_mep > E_ben ? E_mep : E_ben));

// mast-cell / basophil effector activation: occupancy x receptor density
double EFF = (ige_f/(ige_f + Kd_fcer))*pos(FCER);

// plug clearance -- the parameter the biologic moves and the azole cannot
double kout_plug = kout0*(1.0 + cs_plug*CS)
                   / ((1.0 + g_epx*pos(EPX))*(1.0 + g_br*pos(BRON)/BRONmax));

// ===========================================================================
//  antigen, Th2, cytokines
// ===========================================================================
dxdt_AG = kag*(pos(FLUM) + ag_plug*pos(FPLG)) - kagd*AG;

double drive = th2_base + (1.0 - th2_base)*pos(AG)/(pos(AG) + Kag2)*(1.0 - E_tez);
dxdt_TH2  = kin2*drive*(1.0 - w_t2*CS) - kout2*TH2;
dxdt_IL13 = s13*pos(TH2)*(1.0 - w_cyt*CS) - kd13*IL13;
dxdt_IL5  = s5 *pos(TH2)*(1.0 - w_cyt*CS) - kd5 *IL5;

// ===========================================================================
//  plasma cells and the IgE / omalizumab TMDD (rapid-equilibrium form)
// ===========================================================================
dxdt_PC = spc*pos(AG)/(pos(AG) + Kagpc)*(1.0 + a4*IL13e)*(1.0 - w_csr*CS)
          + pc0 - kpc*PC;

// TI and TO are TOTALS; the two elimination routes differ, which is the whole
// reason total IgE rises while free IgE falls.
dxdt_TI   = kige*pos(PC) - kdegE*ige_f - kelCX*ige_c;
dxdt_OMAD = -ka_oma*OMAD;
dxdt_TO   = ka_oma*OMAD*F_oma*OMA_MG_TO_NM - kel_oma*pos(oma_f)
            - kelCX*(ige_c/eta_oma);

dxdt_FCER = ksyn_fcer*(1.0 + b_fcer*ige_f/(ige_f + Kup_fcer))
            - kdeg_fcer*(1.0 + b_fcer)*pos(FCER);

// ===========================================================================
//  eosinophils and granule chemistry
// ===========================================================================
dxdt_EOSB = kin_e*(1.0 + a5*IL5e)
            - kout_e*pos(EOSB)*(1.0 + cs_eos*CS)
            - kADCC*pos(BENR)/(pos(BENR) + EC50_benr)*pos(EOSB);
double egress = ktr_e*(pos(EOSB)/1000.0)*(1.0 + aeot*IL13e)*(1.0 - dup_egress*E_dup);
dxdt_EOSA = egress - kout_ea*pos(EOSA)*(1.0 + cs_apop*CS);
dxdt_EPX  = sepx*pos(EOSA) - kepx*EPX;

// ===========================================================================
//  the mucus plug
// ===========================================================================
double mdrive = muc0 + smuc*IL13e/(IL13e + Kmuc) + smuc_ige*EFF;
dxdt_PLUG = mdrive*(1.0 - cs_muc*CS) - kout_plug*PLUG;

// ===========================================================================
//  THE SANCTUARY PARTITION -- the structural core
// ===========================================================================
double cu_az = fu_itra*pos(ITRA) + fu_oh*pos(OHIT)*pot_oh + fu_vori*pos(VORI)*pot_vori;
double E_az  = (AF_OFF > 0.5) ? 0.0 : Emax_af *cu_az   /(cu_az   + EC50_af);
double E_amb = (AF_OFF > 0.5) ? 0.0 : Emax_amb*pos(AMB)/(pos(AMB) + EC50_amb);
double E_lum = E_az + E_amb;
double E_plg = f_pen*E_az + f_pen_amb*E_amb;        // the sanctuary discount
double host  = k_host*(1.0 - cs_imm*CS);           // steroids REDUCE this
double entrap = kin_f*pos(PLUG)*pos(FLUM)*pos(1.0 - pos(FPLG));

dxdt_FLUM = gl*pos(FLUM)*pos(1.0 - pos(FLUM) - pos(FPLG))
            + kout_plug*pos(FPLG) - entrap - host*pos(FLUM) - E_lum*pos(FLUM) + seed;
dxdt_FPLG = gp*pos(FPLG)*pos(1.0 - pos(FPLG)) + entrap
            - kout_plug*pos(FPLG) - E_plg*pos(FPLG);

// ===========================================================================
//  antifungal PK
// ===========================================================================
dxdt_AITR = -ka_itra*AITR;
double cl_itra = kel_itra/(1.0 + pos(ITRA)/Ki_auto);      // autoinhibition
dxdt_ITRA = ka_itra*AITR*F_itra/V_itra - cl_itra*pos(ITRA);
dxdt_OHIT = fm_oh*cl_itra*pos(ITRA) - kel_oh*pos(OHIT);

dxdt_AVOR = -ka_vori*AVOR;
dxdt_VORI = ka_vori*AVOR*F_vori/V_vori
            - (Vmax_vori*cyp2c19*pos(VORI)/(Km_vori + pos(VORI))
               + CLlin_vori*pos(VORI))/V_vori;

dxdt_AMB  = -kel_amb*pos(AMB);

// ===========================================================================
//  corticosteroid PK with the SEQUENTIAL gut + hepatic interaction
// ===========================================================================
double gutP = (1.0 - s3A4g_pred)  + s3A4g_pred *I3A4;
double hepP = (1.0 - s3A4h_pred)  + s3A4h_pred *I3A4;
double gutM = (1.0 - s3A4g_mpred) + s3A4g_mpred*I3A4;
double hepM = (1.0 - s3A4h_mpred) + s3A4h_mpred*I3A4;
double gutB = (1.0 - s3A4g_bud)   + s3A4g_bud  *I3A4;
double hepB = (1.0 - s3A4h_bud)   + s3A4h_bud  *I3A4;
double FP = F_pred /gutP;  if (FP > 1.0) FP = 1.0;
double FM = F_mpred/gutM;  if (FM > 1.0) FM = 1.0;
double FB = F_bud  /gutB;  if (FB > 1.0) FB = 1.0;

dxdt_APRD = -ka_pred*APRD;
dxdt_PRED = ka_pred*APRD*FP/V_pred - kel_pred*hepP*pos(PRED);
dxdt_AMPD = -ka_mpred*AMPD;
dxdt_MPRD = ka_mpred*AMPD*FM/V_mpred - kel_mpred*hepM*pos(MPRD);
dxdt_ABUD = -ka_bud*ABUD;
dxdt_BUD  = ka_bud*ABUD*FB/V_bud - kel_bud*hepB*pos(BUD);

// ===========================================================================
//  HPA axis, bone, glucose, cumulative steroid
// ===========================================================================
dxdt_CORT  = kin_c*(1.0 - Imax_c*CS) - kout_c*pos(CORT);
dxdt_BMD   = -kbmd*CS + krep_bmd*(1.0 - BMD);
dxdt_HBA1C = ((a1c0 + a1c_gain*CS) - HBA1C)/tau_a1c;
dxdt_CUMO  = 0.0;                       // incremented by dosing records

// ===========================================================================
//  biologic PK
// ===========================================================================
dxdt_MEPD = -ka_mep*MEPD;  dxdt_MEPO = ka_mep*MEPD*F_mep/V_mep - kel_mep*pos(MEPO);
dxdt_BEND = -ka_ben*BEND;  dxdt_BENR = ka_ben*BEND*F_ben/V_ben - kel_ben*pos(BENR);
dxdt_DUPD = -ka_dup*DUPD;  dxdt_DUPI = ka_dup*DUPD*F_dup/V_dup - kel_dup*pos(DUPI);
dxdt_TEZD = -ka_tez*TEZD;  dxdt_TEZE = ka_tez*TEZD*F_tez/V_tez - kel_tez*pos(TEZE);

// ===========================================================================
//  structure, hazard, exposure integrals
// ===========================================================================
double haz = h0*exp(b1_eos*pos(EOSA) + b2_plug*pos(PLUG)
                    + b3_fun*(pos(FLUM) + pos(FPLG)) + b4_ige*EFF);
dxdt_CHAZ  = haz;
dxdt_AUCP  = peq;
dxdt_AUCCS = CS;
dxdt_BRON  = (BRON < BRONmax) ? (kb1*pos(pos(PLUG) - Pthr) + kb2*haz) : 0.0;

$TABLE
double ige_c_o = bind11(eta_oma*pos(TO), pos(TI), Kd_oma);
double ige_f_o = pos(TI) - ige_c_o;
double peq_o = fu_pred*pos(PRED) + fu_mpred*pos(MPRD)*pot_mpred
             + fu_bud*pos(BUD)*pot_bud;
double CS_o  = (peq_o > 0.0) ? Emax_cs*peq_o/(peq_o + EC50_cs) : 0.0;
double EFF_o = (ige_f_o/(ige_f_o + Kd_fcer))*pos(FCER);
double IL13e_o = pos(IL13)*(1.0 - Emax_dup*pos(DUPI)/(pos(DUPI) + EC50_dup));
double ahr = w_ahr_il13*IL13e_o/(IL13e_o + Ka_fev) + w_ahr_ige*EFF_o;

capture IGE_TOTAL  = pos(TI)*IGE_NM_TO_IU;                  // IU/mL -- RISES on omalizumab
capture IGE_FREE   = ige_f_o*IGE_NM_TO_IU;                  // IU/mL -- FALLS on omalizumab
capture IGE_FREE_NG = ige_f_o*IGE_NM_TO_NGML;               // ng/mL, the 25 ng/mL target
capture IGE_CPLX   = ige_c_o*IGE_NM_TO_IU;
capture FCER_OCC   = ige_f_o/(ige_f_o + Kd_fcer);
capture EFFECTOR   = EFF_o;
capture PLUG_SCORE = (6.0*pos(PLUG) < 18.0) ? 6.0*pos(PLUG) : 18.0;
capture FTOT       = pos(FLUM) + pos(FPLG);
capture SANCT_FRAC = (FTOT > 1e-12) ? pos(FPLG)/FTOT : 0.0;  // sanctuary share
capture KOUT_PLUG  = kout0*(1.0 + cs_plug*CS_o)
                     /((1.0 + g_epx*pos(EPX))*(1.0 + g_br*pos(BRON)/BRONmax));
capture E_LUM_CAP  = Emax_af*(fu_itra*pos(ITRA) + fu_oh*pos(OHIT)*pot_oh
                              + fu_vori*pos(VORI)*pot_vori)
                     /((fu_itra*pos(ITRA) + fu_oh*pos(OHIT)*pot_oh
                        + fu_vori*pos(VORI)*pot_vori) + EC50_af);
// the analytic sanctuary floor: the luminal kill rate below which the plug
// population sustains itself no matter what the drug does
capture E_FLOOR    = (gp - KOUT_PLUG > 0.0) ? (gp - KOUT_PLUG)/f_pen : 0.0;
capture SANCT_OK   = (E_LUM_CAP > E_FLOOR) ? 1.0 : 0.0;
capture I3A4_CAP   = 1.0/(1.0 + fu_itra*pos(ITRA)/Ki3A4_itra
                              + fu_oh*pos(OHIT)/Ki3A4_oh
                              + fu_vori*pos(VORI)/Ki3A4_vori);
capture CS_CAP     = CS_o;
capture FEV1       = FEV1_0*(1.0 - cp_fev*pos(PLUG)/(pos(PLUG) + Kp_fev))
                          *(1.0 - cb_fev*pos(BRON)/BRONmax)
                          *(1.0 - ca_fev*ahr);
capture HAZ_YR     = 365.0*h0*exp(b1_eos*pos(EOSA) + b2_plug*pos(PLUG)
                                  + b3_fun*FTOT + b4_ige*EFF_o);
capture IGE_RESP   = 0.0;   // filled by post-processing: total-IgE response
                            // criterion; see the warning in the header

$CAPTURE AG TH2 IL13 IL5 PC EOSB EOSA EPX PLUG FLUM FPLG ITRA OHIT VORI AMB
$CAPTURE PRED MPRD BUD CORT BMD HBA1C CUMO MEPO BENR DUPI TEZE BRON CHAZ
$CAPTURE AUCP AUCCS


## ===========================================================================
##  SCENARIO BUILDERS  (R code -- source this file's tail after mread())
## ===========================================================================
## Not parsed by mrgsolve; keep below the model spec.
##
## abpa_scenario(name, days = 364, wt = 70)
##   "untreated"            no therapy
##   "pred_isham"           ISHAM medium-dose prednisolone taper
##   "pred_maint"           prednisolone 10 mg/d maintenance
##   "itraconazole"         itraconazole 200 mg BID
##   "pred_itra"            ISHAM taper + itraconazole
##   "voriconazole"         voriconazole 200 mg BID (set cyp2c19 for PM)
##   "neb_amb"              nebulised liposomal amphotericin B 10 mg q48h
##   "omalizumab"           omalizumab 375 mg q2w  (WATCH IGE_TOTAL vs IGE_FREE)
##   "omalizumab_flux"      omalizumab dosed to the model's flux requirement
##   "mepolizumab"          mepolizumab 100 mg q4w
##   "benralizumab"         benralizumab 30 mg q4w x3 then q8w
##   "dupilumab"            dupilumab 300 mg q2w
##   "dupi_itra"            dupilumab + itraconazole (the synergy arm)
##   "bud_itra"             inhaled budesonide + itraconazole (the trap)
##   "lavage_itra"          bronchoscopic plug removal at day 0 + itraconazole
##   "triple"               ISHAM taper + itraconazole + dupilumab
##
## Example:
##   library(mrgsolve); library(dplyr)
##   mod <- mread("abpa_mrgsolve_model", ".")
##   out <- mod %>% ev(abpa_scenario("dupi_itra")) %>% mrgsim(end = 364, delta = 1)
##   plot(out, FPLG + SANCT_FRAC + PLUG_SCORE + FEV1 ~ time)
##
## The counterfactual arms used in the report are parameter overrides, not
## events -- e.g. the DDI-off arm is
##   mod %>% param(DDI_OFF = 1) %>% ev(abpa_scenario("pred_itra")) %>% mrgsim()
## and the antifungal-off arm is  param(AF_OFF = 1).
##
## ---------------------------------------------------------------------------
## if (FALSE) {   ## <- unwrap to use
##
## abpa_scenario <- function(name, days = 364, wt = 70) {
##   ev0 <- function(...) mrgsolve::ev(...)
##   bind <- function(...) do.call(c, Filter(Negate(is.null), list(...)))
##
##   itra <- ev0(time = 0, amt = 400, cmt = "AITR", ii = 1, addl = days - 1)
##   vori <- ev0(time = 0, amt = 400, cmt = "AVOR", ii = 1, addl = days - 1)
##   amb  <- ev0(time = 0, amt = 10,  cmt = "AMB",  ii = 2, addl = days/2 - 1)
##   bud  <- ev0(time = 0, amt = 1.6, cmt = "ABUD", ii = 1, addl = days - 1)
##   oma  <- function(mg) ev0(time = 0, amt = mg, cmt = "OMAD", ii = 14,
##                            addl = floor(days/14) - 1)
##   mepo <- ev0(time = 0, amt = 100, cmt = "MEPD", ii = 28, addl = floor(days/28) - 1)
##   dupi <- ev0(time = 0, amt = 300, cmt = "DUPD", ii = 14, addl = floor(days/14) - 1)
##   benr <- c(ev0(time = 0,  amt = 30, cmt = "BEND", ii = 28, addl = 2),
##             ev0(time = 84, amt = 30, cmt = "BEND", ii = 56,
##                 addl = floor((days - 84)/56) - 1))
##   ## ISHAM medium-dose taper: 0.5 mg/kg/d x 2 wk, alternate-day x 8 wk,
##   ## then -5 mg every 2 wk to zero by ~5 months
##   isham <- function() {
##     d0 <- 0.5*wt; out <- NULL
##     for (day in 0:(days - 1)) {
##       dose <- if (day < 14) d0
##               else if (day < 70) { if ((day - 14) %% 2 == 0) d0 else 0 }
##               else if (day < 154) max(d0 - 5*((day - 70) %/% 14 + 1), 0)
##               else 0
##       if (dose > 0) out <- c(out, ev0(time = day, amt = dose, cmt = "APRD"))
##     }
##     out
##   }
##   maint <- ev0(time = 0, amt = 10, cmt = "APRD", ii = 1, addl = days - 1)
##
##   switch(name,
##     untreated       = ev0(time = 0, amt = 0, cmt = "AITR"),
##     pred_isham      = isham(),
##     pred_maint      = maint,
##     itraconazole    = itra,
##     pred_itra       = c(isham(), itra),
##     voriconazole    = vori,
##     neb_amb         = amb,
##     omalizumab      = oma(375),
##     omalizumab_flux = oma(1200),
##     mepolizumab     = mepo,
##     benralizumab    = benr,
##     dupilumab       = dupi,
##     dupi_itra       = c(dupi, itra),
##     bud_itra        = c(bud, itra),
##     lavage_itra     = itra,   ## plug removal is a state reset, not an event:
##                               ## use  init(FPLG = 0.05*FPLG, PLUG = 0.05*PLUG)
##     triple          = c(isham(), itra, dupi),
##     stop("unknown scenario: ", name))
## }
##
## }  ## end if (FALSE)
## ---------------------------------------------------------------------------
##
## READ THIS BEFORE PLOTTING ANYTHING FROM AN OMALIZUMAB ARM
## --------------------------------------------------------
## IGE_TOTAL and IGE_FREE move in OPPOSITE directions on omalizumab.  Total IgE
## rises about three-fold while free IgE falls about 96%.  Any figure, table or
## response criterion built on total IgE will score a working drug as a failure.
## This is not a bug and it is not a parameter choice: it follows from
## kelCX < kdegE, which is guaranteed as long as an IgG-like complex is cleared
## more slowly than free IgE.  Sweep kelCX over any plausible range and the
## direction never changes; only the size does.
##
## READ THIS BEFORE INTERPRETING ANY ITRACONAZOLE ARM
## -------------------------------------------------
## If the arm also contains a corticosteroid, part of what looks like antifungal
## efficacy is the CYP3A4 interaction raising steroid exposure.  Run the same
## arm with param(DDI_OFF = 1) to see how much.  On a prednisolone backbone the
## AUC rise is small (1.12x); on inhaled budesonide it is 4x, and there the
## interaction accounts for essentially the whole apparent benefit -- and can
## also produce adrenal suppression on its own, with nothing in the disease
## model having changed.
## ===========================================================================
