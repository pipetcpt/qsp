## =====================================================================
##  icp_mrgsolve_model.R
##  Intrahepatic Cholestasis of Pregnancy (ICP) — QSP model for mrgsolve
##  임신성 간내 담즙정체 — mrgsolve QSP 모델
##
##  78 ODEs = 5 bile acid species x 7 compartments + 43 global states
##
##  THE ONE STRUCTURAL IDEA
##  -----------------------
##  ICP is measured with a single number — total serum bile acid from a
##  3alpha-hydroxysteroid dehydrogenase assay — and that number is not
##  what harms either patient.  So this model refuses to use it as a
##  state variable.  It carries FIVE bile acid species (cholic,
##  chenodeoxycholic, deoxycholic, lithocholic, ursodeoxycholic) through
##  hepatocyte, bile, small intestine, colon, maternal plasma, fetal
##  plasma and fetal liver, and it reconstructs the assay at the end as
##  an observation — including the ursodeoxycholate and hyocholate that
##  the assay also counts but nobody attributes to it.
##
##  TWO CAUSAL CHAINS OUT OF ONE LESION
##  -----------------------------------
##  BILE ACID AXIS   sex steroids -> BSEP inhibition -> maternal pool ->
##                   placental carrier -> fetal hydrophobic load ->
##                   Cx43 uncoupling + Ca2+ overload -> arrhythmia ->
##                   stillbirth.  Moved by UDCA, cholestyramine, an IBAT
##                   inhibitor.
##  ITCH AXIS        sex steroids -> autotaxin -> LPA -> MRGPRX4/TRPV1 ->
##                   spinal GRPR -> pruritus.  Coupled to bile acid only
##                   weakly and saturably.  Moved by rifampicin.
##  The two axes barely share a parameter, which is why a trial powered
##  on one cannot rank drugs for the other.  See section F of
##  icp_calibration.py.
##
##  WHERE THE 100 umol/L THRESHOLD COMES FROM
##  -----------------------------------------
##  Not from the dose-response and not from the placenta.  Re-fitting the
##  Ovadia 2019 stillbirth strata (0.13% / 0.28% / 3.44% for <40 / 40-99
##  / >=100 umol/L) against three different driver variables:
##      hazard on maternal total bile acid   needs exponent 2.95
##      hazard on fetal cytotoxic load       needs exponent 2.55
##      hazard on myocardial arrhythmia idx  needs exponent 1.55
##  The arrhythmia index spans 37-fold across the same three pregnancies
##  that span 8.7-fold in maternal bile acid.  The steepness lives in
##  connexin-43 uncoupling, which is a threshold process.  Ablating it
##  (or linearising it) destroys the threshold; ablating placental
##  carrier saturation does not.  The map was drawn expecting the
##  opposite; the fit corrected it.
##
##  VALIDATION
##  ----------
##  No R runtime was available where this was written, so every equation
##  below was first implemented and integrated in dependency-free Python
##  (icp_reference_model.py, fixed-step RK4), and this file is a
##  line-by-line port of what that produced.  The Python runs exposed
##  NINE real defects; each is flagged [DEFECT n] here and documented at
##  the point of the fix in the Python file.  The initial conditions in
##  $INIT below are the converged pre-third-trimester state of that
##  reference implementation, not guesses.
##
##  Two deliberate differences from the Python file, both harmless for
##  every endpoint reported (all of which integrate over weeks):
##    * dosing here uses proper discrete mrgsolve events; the reference
##      implementation used equivalent daily input rates.
##    * mrgsolve integrates with LSODA, which handles the fast maternal
##      plasma compartment natively; the RK4 reference needed the
##      apparent-volume formulation described at VSP below.
##
##  USAGE
##      library(mrgsolve); library(dplyr)
##      mod <- mread("icp_mrgsolve_model", "path/to/dir")
##      out <- mod %>% mrgsim(end = 154, delta = 1) %>% as_tibble()
##  TIME is days from GA 20+0. Delivery happens at TDEL (days from GA 20).
##  Scenarios (35 of them) are at the bottom of this file.
## =====================================================================

$PROB
# Intrahepatic Cholestasis of Pregnancy (ICP)
# 5 species x 7 compartments + 43 global = 78 ODEs

$PARAM @annotated
// ---- gestational frame ----------------------------------------------
GA0     :  20.0  : Gestational age at TIME = 0 (weeks)
TDEL    : 133.0  : Days from TIME=0 to delivery (133 d = GA 39+0)

// ---- constitutional / genetic susceptibility ------------------------
// These five multipliers are the model's patient descriptors.  Severity
// strata are DEFINED by the bile acid band they produce, exactly as the
// literature defines them, not by a genotype fixed in advance.
GBSEP   :  1.00  : ABCB11/BSEP transport capacity (1 = wild type)
GMDR3   :  1.00  : ABCB4/MDR3 phosphatidylcholine floppase capacity
GFIC1   :  1.00  : ATP8B1/FIC1 capacity (multiplies BSEP)
GSULT   :  1.00  : SULT2A1 steroid sulfotransferase activity
GFXR    :  1.00  : NR1H4/FXR responsiveness
GP      :  1.00  : Placental fetal-to-maternal transport capacity
TWIN    :  1.00  : Sex-steroid load multiplier (1.55 for twins)

// ---- volumes ---------------------------------------------------------
VHL     : 15.0   : Effective hepatocyte distribution volume (L)
VIL     :  1.0   : Small-intestinal luminal volume (L)
VCO     :  0.5   : Colonic luminal volume (L)
// VSP is the APPARENT volume of distribution of total bile acid, not the
// plasma volume: >95% of circulating bile acid is albumin-bound and the
// bound pool equilibrates with interstitial fluid within minutes, so the
// measured serum concentration behaves as amount/Vd with Vd ~ 4x plasma
// volume.  [DEFECT 1] Using plasma volume put a 7-minute half-life on a
// 4 L compartment; the RK4 reference went negative at dt = 0.02 d.
VSP     : 18.0   : Apparent Vd of total bile acid, maternal (L)

// ---- synthesis and its feedback -------------------------------------
VSYN    : 520.0  : De-novo bile acid synthesis at CYP7A1 = 1 (umol/d)
FRCA    :  0.55  : CYP8B1 branch fraction to cholic acid
KS_C7   :  1.0   : CYP7A1 synthesis rate constant
KD_C7   :  1.0   : CYP7A1 degradation rate constant (1/d)
KSHP    :  1.0   : SHP potency on CYP7A1
KF19    :  1.2   : FGF19 potency on CYP7A1
KS_SHP  :  2.2   : SHP synthesis rate constant
KD_SHP  :  2.2   : SHP degradation rate constant (1/d)
EC50FXR : 120.0  : Hepatocyte FXR EC50, FXR-weighted bile acid (umol/L)
KS_F19  :  3.0   : FGF19 synthesis rate constant
KD_F19  :  3.0   : FGF19 degradation rate constant (1/d)
EC50FXI : 900.0  : Ileal FXR EC50, luminal FXR-weighted bile acid (umol/L)

// ---- MRP2: the canalicular route for UDCA ---------------------------
// [DEFECT 2] UDCA was first modelled as a BSEP substrate inside the
// shared competitive denominator.  That forced a choice between UDCA not
// reaching bile (plasma UDCA 60+ umol/L, never observed) and UDCA
// displacing the endogenous species so that the drug RAISED their serum
// level.  MRP2/ABCC2 is the correct route and it does not compete.
VMRP2   : 24000. : MRP2 maximal UDCA canalicular export (umol/d)
KMMRP2  : 55.0   : MRP2 Km for UDCA (umol/L, effective hepatocyte conc)
EMRP2P  :  0.60  : PXR induction of MRP2 per unit CYP3A above baseline

// ---- BSEP (ABCB11) --------------------------------------------------
// [DEFECT 3] BSEP was initially placed near Vmax (denB ~ 4.4).  A 1.3x
// fall in transport capacity then produced a 6x rise in hepatocyte load
// and the genotype-to-phenotype map overshot the clinical spread of ICP
// by an order of magnitude (term TBA 200-860 umol/L).
VBSEP   : 74000. : BSEP maximal canalicular export (umol/d)
KMB1    : 52.0   : BSEP Km, cholate (umol/L)
KMB2    : 42.0   : BSEP Km, chenodeoxycholate (umol/L)
KMB3    : 42.0   : BSEP Km, deoxycholate (umol/L)
KMB4    : 28.0   : BSEP Km, lithocholate (umol/L)
ETRAF   :  0.14  : UDCA-driven canalicular carrier insertion, Emax
KTRAF   :  6.0   : UDCA concentration for half-maximal insertion (umol/L)
ESAM    :  0.18  : SAMe effect on canalicular export, Emax
KSAM    :  9.0   : SAMe concentration for half-maximal effect (mg/L eq)

// ---- sex-steroid inhibition of BSEP: the lesion ---------------------
KI_E2G  : 11.0   : Estradiol-17b-glucuronide Ki, cis-inhibition (nmol/L)
KI_P4S  : 140.0  : Progesterone sulfate Ki, trans-inhibition (nmol/L)

// ---- basolateral escape valve (MRP3/MRP4/OSTa-b) --------------------
// [DEFECT 4] Written as a proportional FXR term the valve was already
// 76% engaged in normal pregnancy, so it had no headroom to express the
// 5-20x induction human cholestatic liver actually shows.  It is now a
// Hill in FXR occupancy.
KBL     : 39.0   : Basolateral efflux clearance at MRP4 = 1 (L/d)
SBL1    :  1.00  : Basolateral selectivity, cholate
SBL2    :  0.90  : Basolateral selectivity, chenodeoxycholate
SBL3    :  0.85  : Basolateral selectivity, deoxycholate
SBL4    :  0.60  : Basolateral selectivity, lithocholate
SBL5    :  0.18  : Basolateral selectivity, ursodeoxycholate
KS_M4   :  1.0   : MRP4 synthesis rate constant
KD_M4   :  1.0   : MRP4 degradation rate constant (1/d)
EM4     :  7.0   : Maximal cholestatic induction of the escape valve
KM4     :  0.62  : FXR occupancy for half-maximal valve induction
HM4     :  2.5   : Hill coefficient, valve induction
EM4P    :  0.25  : PXR induction of the valve per unit CYP3A

// ---- hepatic uptake / first pass ------------------------------------
EUP1    :  0.990 : Intrinsic first-pass extraction, cholate
EUP2    :  0.985 : Intrinsic first-pass extraction, chenodeoxycholate
EUP3    :  0.980 : Intrinsic first-pass extraction, deoxycholate
EUP4    :  0.970 : Intrinsic first-pass extraction, lithocholate
EUP5    :  0.930 : Intrinsic first-pass extraction, ursodeoxycholate
JHALF   : 800000.: Absorbed load halving first-pass efficiency (umol/d)
KHLSAT  : 300.0  : Hepatocyte load at which NTCP is down-regulated (umol/L)
HHLSAT  :  2.0   : Hill: NTCP down-regulation is a threshold, not linear
KIRIFO  : 12.0   : Rifampicin Ki at OATP1B1/NTCP, total conc (mg/L)
CLHSYS  : 706.0  : Hepatic clearance of recirculating systemic BA (L/d)

// ---- bile flow and MDR3 --------------------------------------------
KBILE   :  7.0   : Biliary emptying rate constant (1/d)
VPC     : 6000.  : Biliary phosphatidylcholine flux at MDR3 = 1 (umol/d)
BSPC0   :  2.8   : Normal biliary bile salt : phosphatidylcholine ratio

// ---- intestine -----------------------------------------------------
VASBT   : 20500. : ASBT maximal ileal reabsorption (umol/d)
KMA1    : 300.0  : ASBT Km, cholate (umol/L)
KMA2    : 260.0  : ASBT Km, chenodeoxycholate (umol/L)
KMA3    : 260.0  : ASBT Km, deoxycholate (umol/L)
KMA4    : 200.0  : ASBT Km, lithocholate (umol/L)
KMA5    : 5000.  : ASBT Km, ursodeoxycholate (umol/L)
KPAS1   :  0.004 : Jejunal passive absorption, cholate (1/d)
KPAS2   :  0.010 : Jejunal passive absorption, chenodeoxycholate (1/d)
KPAS3   :  0.026 : Jejunal passive absorption, deoxycholate (1/d)
KPAS4   :  0.030 : Jejunal passive absorption, lithocholate (1/d)
KPAS5   :  0.007 : Jejunal passive absorption, ursodeoxycholate (1/d)
KTR     :  0.42  : Ileum-to-colon transit rate constant (1/d)
KSEQ    :  0.45  : Maximal resin sequestration rate constant (1/d)
KCHOL   :  3.0   : Cholestyramine luminal amount, half-maximal (g)
FBI1    :  1.00  : Resin binding selectivity, cholate
FBI2    :  1.05  : Resin binding selectivity, chenodeoxycholate
FBI3    :  1.05  : Resin binding selectivity, deoxycholate
FBI4    :  1.10  : Resin binding selectivity, lithocholate
FBI5    :  0.95  : Resin binding selectivity, ursodeoxycholate
KIODEV  :  0.9   : IBAT inhibitor luminal amount for 50% ASBT block (umol)
KDH_CA  :  1.10  : Microbial 7a-dehydroxylation, CA to DCA (1/d)
KDH_CD  :  1.25  : Microbial 7a-dehydroxylation, CDCA to LCA (1/d)
KDH_UD  :  0.20  : Microbial 7b-dehydroxylation, UDCA to LCA (1/d)
KCOL1   :  0.06  : Colonic passive absorption, cholate (1/d)
KCOL2   :  0.10  : Colonic passive absorption, chenodeoxycholate (1/d)
KCOL3   :  0.55  : Colonic passive absorption, deoxycholate (1/d)
KCOL4   :  0.60  : Colonic passive absorption, lithocholate (1/d)
KCOL5   :  0.08  : Colonic passive absorption, ursodeoxycholate (1/d)
KFEC    :  1.05  : Faecal elimination rate constant (1/d)

// ---- hepatic lithocholate detoxification ---------------------------
// [DEFECT 5] Without this route the model had no way to dispose of the
// lithocholate that gut bacteria make from UDCA, so UDCA therapy shifted
// the fetal pool 43% more hydrophobic and the drug came out harmful.
VLCAS   : 900.0  : SULT2A1 maximal LCA sulfation (umol/d)
KMLCAS  :  6.0   : SULT2A1 Km for LCA (umol/L)

// ---- renal ---------------------------------------------------------
CLR1    :  0.170 : Renal clearance, cholate (L/d)
CLR2    :  0.160 : Renal clearance, chenodeoxycholate (L/d)
CLR3    :  0.145 : Renal clearance, deoxycholate (L/d)
CLR4    :  0.095 : Renal clearance, lithocholate (L/d)
// [DEFECT 6] UDCA glucuronides are renally cleared far better than
// endogenous conjugates.  Given the endogenous value, UDCA accumulated
// behind the cholestatic block to plasma concentrations of 100+ umol/L.
CLR5    :  3.000 : Renal clearance, ursodeoxycholate (L/d)
ESULT   :  4.0   : SULT2A1 induction effect on renal BA clearance
KS_S    :  1.0   : SULT2A1 synthesis rate constant
KD_S    :  1.0   : SULT2A1 degradation rate constant (1/d)
EPXR_S  :  0.85  : PXR induction of SULT2A1 per unit CYP3A

// ---- sex steroids --------------------------------------------------
E2BASE  : 20.0   : Estradiol at GA 20 (nmol/L)
E2K     :  0.0627: Estradiol exponential rise constant (1/wk)
KOUTE2  :  2.77  : Estradiol elimination rate constant (1/d)
P4BASE  : 200.0  : Progesterone at GA 20 (nmol/L)
P4K     :  0.0550: Progesterone exponential rise constant (1/wk)
KOUTP4  :  1.39  : Progesterone elimination rate constant (1/d)
KFE2G   :  0.185 : Estradiol-17b-glucuronide formation rate constant (1/d)
KOUTE2G :  3.10  : Estradiol-17b-glucuronide elimination (1/d)
EUGT    :  0.15  : Rifampicin co-induction of UGT1A1
KFP4S   :  0.088 : Progesterone sulfate formation rate constant (1/d)
KOUTP4S :  0.55  : Progesterone sulfate elimination rate constant (1/d)

// ---- placental transfer --------------------------------------------
VP      : 88.0   : Maximal fetal-to-maternal carrier capacity (umol/d)
KMP1    : 26.0   : Placental carrier Km, cholate (umol/L)
KMP2    : 20.0   : Placental carrier Km, chenodeoxycholate (umol/L)
KMP3    : 20.0   : Placental carrier Km, deoxycholate (umol/L)
KMP4    : 15.0   : Placental carrier Km, lithocholate (umol/L)
KMP5    : 55.0   : Placental carrier Km, ursodeoxycholate (umol/L)
TRANSIN :  0.30  : Maternal-side occupancy of the placental carrier
// [DEFECT 7] The diffusive component was first written one-way
// (maternal -> fetal only), as several published ICP models write it.
// The fetal compartment then had no concentration-independent escape
// once the carrier saturated and the fetal pool grew without bound
// (cord TBA 700 umol/L).  Diffusion is bidirectional.
PS1     :  0.550 : Diffusive placental permeability, cholate (L/d)
PS2     :  0.750 : Diffusive placental permeability, chenodeoxycholate (L/d)
PS3     :  0.850 : Diffusive placental permeability, deoxycholate (L/d)
PS4     :  0.950 : Diffusive placental permeability, lithocholate (L/d)
PS5     :  0.150 : Diffusive placental permeability, ursodeoxycholate (L/d)
SF      : 13.0   : Fetal de-novo bile acid synthesis (umol/d)
FSF1    :  0.60  : Fetal synthesis fraction to cholate
FSF2    :  0.40  : Fetal synthesis fraction to chenodeoxycholate
CLFL    :  0.55  : Fetal plasma-to-liver clearance (L/d)
KFLO    :  2.20  : Fetal liver-to-plasma rate constant (1/d)
KFB     :  0.55  : Fetal biliary secretion rate constant (1/d)
KUR     :  0.200 : Fetal urinary bile acid clearance (L/d per L Vfp)
KSW     :  1.30  : Fetal swallowing clearance of amniotic fluid (1/d)
KIMA    :  0.55  : Intramembranous absorption of amniotic fluid (1/d)

// ---- fetal myocardium: where the threshold lives -------------------
IMAXGJ  :  0.85  : Maximal connexin-43 uncoupling
IC50GJ  : 30.0   : Fetal cytotoxic load for half-maximal uncoupling
HGJ     :  1.60  : Hill coefficient, connexin-43 uncoupling
KS_GJ   :  1.0   : Gap junction synthesis rate constant
KD_GJ   :  1.0   : Gap junction turnover rate constant (1/d)
ECA     :  2.50  : Maximal calcium overload above baseline
KCA     : 45.0   : Fetal cytotoxic load for half-maximal Ca overload
KS_CA   :  1.0   : Calcium index synthesis rate constant
KD_CA   :  1.0   : Calcium index turnover rate constant (1/d)

// ---- placental vasculature / fetal hypoxia -------------------------
EV      :  1.0   : Maximal chorionic vasoconstriction
KV      : 40.0   : Fetal cytotoxic load for half-maximal vasoconstriction
KS_V    :  1.0   : Vasoconstriction synthesis rate constant
KD_V    :  1.0   : Vasoconstriction resolution rate constant (1/d)
KT      : 55.0   : Fetal cytotoxic load for half-maximal trophoblast damage
KS_T    :  1.0   : Trophoblast damage rate constant
KD_T    :  1.0   : Trophoblast repair rate constant (1/d)
KS_H    :  1.0   : Hypoxia index formation rate constant
KD_H    :  1.0   : Hypoxia index resolution rate constant (1/d)
WT      :  0.65  : Weight of trophoblast damage in the hypoxia index

// ---- outcome hazards (the fitted constants; see icp_calibration.py) -
HSB0    :  1.0887e-5 : Background stillbirth hazard (1/d) [FITTED]
HSBSC   :  3.1896e-3 : Scale on the arrhythmia index (1/d) [FITTED]
HN      :  1.55  : Exponent on the arrhythmia index [FITTED]
KHY     :  0.55  : Amplification of the hazard by fetal hypoxia
HPT0    :  5.4e-4: Background spontaneous preterm birth hazard (1/d)
KOT     :  2.20  : Oxytocin receptor effect on the preterm hazard
KUT     :  1.00  : Uterine activity effect on the preterm hazard
HM0     :  3.7e-3: Background meconium passage hazard (1/d)
HMSC    :  2.52e-2 : Bile-acid-driven meconium hazard scale (1/d)
KMG     : 25.0   : Fetal cytotoxic load for half-maximal meconium hazard

// ---- myometrium ----------------------------------------------------
KS_O    :  1.0   : Oxytocin receptor synthesis rate constant
KD_O    :  1.0   : Oxytocin receptor turnover rate constant (1/d)
EO      :  1.60  : Maximal bile-acid induction of oxytocin receptor
KO      : 22.0   : Maternal cytotoxic load for half-maximal OTR induction
KS_U    :  1.0   : Uterine activity formation rate constant
KD_U    :  1.0   : Uterine activity resolution rate constant (1/d)

// ---- pruritus: the second axis -------------------------------------
// Autotaxin is driven mainly by the sex-steroid load and only weakly,
// and SATURABLY, by circulating bile acid.  [DEFECT 8] Driven from
// maternal bile acid with KA_CH = 17 umol/L, UDCA reduced itch by 3.3 cm
// on a 10 cm scale.  The trial figure is 0.7 cm, serum autotaxin does
// not fall with UDCA at all, and in a large minority of women itch
// precedes the biochemical abnormality by weeks.
KS_A    :  1.0   : Autotaxin synthesis rate constant
KD_A    :  1.0   : Autotaxin turnover rate constant (1/d)
EA_E2   :  1.05  : Estradiol drive on autotaxin
KA_E2   : 48.0   : Estradiol for half-maximal autotaxin drive (nmol/L)
EA_P4S  :  3.00  : Progesterone sulfate drive on autotaxin
KA_P4S  : 200.0  : Progesterone sulfate, half-maximal drive (nmol/L)
HP4S    :  2.0   : Hill coefficient, progesterone sulfate drive
EA_CH   :  0.35  : Bile acid drive on autotaxin (weak)
KA_CH   :  3.0   : Maternal cytotoxic load, half-maximal (umol/L)
ERIFA   :  0.62  : Rifampicin (PXR) suppression of autotaxin
KS_L    :  1.0   : LPA formation rate constant
KD_L    :  1.0   : LPA turnover rate constant (1/d)
IMAXI   :  1.0   : Maximal LPA-driven itch signal
ILP50   :  3.20  : LPA for half-maximal itch signal
HLP     :  4.0   : Hill coefficient, LPA to itch (LPAR-TRPV1-GRPR cascade)
WBA     :  0.30  : Direct bile acid (MRGPRX4) contribution to itch
KBA     : 25.0   : Maternal cytotoxic load, half-maximal MRGPRX4 (umol/L)
KSENS   :  0.90  : Scratch-itch central sensitisation gain
KSC     :  0.45  : Central itch state for half-maximal sensitisation
KS_I    :  1.0   : Central itch state formation rate constant
KD_I    :  1.0   : Central itch state resolution rate constant (1/d)
ENTX    :  0.30  : Naltrexone maximal effect on the central itch state
EAH     :  0.22  : Antihistamine maximal effect (sedation)
VSCALE  :  4.81  : cm on a 10 cm VAS per unit central itch state [FITTED]
KS_SL   :  1.0   : Sleep disturbance formation rate constant
KD_SL   :  1.0   : Sleep disturbance resolution rate constant (1/d)

// ---- hepatocellular injury -----------------------------------------
KS_R    :  1.0   : Oxidative stress formation rate constant
KD_R    :  1.0   : Oxidative stress resolution rate constant (1/d)
KR      : 65.0   : Hepatocyte hydrophobic load, half-maximal stress
HR      :  3.0   : Hill coefficient, hydrophobic load to stress
WBSPC   :  0.62  : Weight of biliary bile salt:PC excess in injury
EUDROS  :  0.42  : UDCA cytoprotection, Emax
KUDROS  :  7.0   : UDCA concentration for half-maximal cytoprotection
ESAM_R  :  0.55  : SAMe effect on oxidative stress resolution
WRIFROS :  0.28  : Rifampicin idiosyncratic hepatotoxic contribution
KRIFROS :  3.0   : Rifampicin concentration, half-maximal (mg/L)
KS_ALT  : 20.0   : Baseline ALT production (U/L per day)
KD_ALT  :  1.0   : ALT elimination rate constant (1/d)
EALT    :  9.0   : Oxidative stress amplification of ALT

// ---- fat-soluble vitamins and coagulation --------------------------
MIC50   : 1136.  : Luminal bile acid for half-maximal micelle capacity
KIN_K   :  1.0   : Vitamin K absorption rate constant
KD_K    :  1.0   : Vitamin K turnover rate constant (1/d)
DIETK   :  1.0   : Dietary vitamin K intake (relative)
VKDOSE  :  0.0   : Supplemental vitamin K (relative units, 1.6 = 10 mg/d)
VKSTART : 98.0   : Vitamin K start time (days from TIME 0)
KS_PI   :  1.0   : PIVKA-II formation rate constant
KD_PI   :  1.0   : PIVKA-II elimination rate constant (1/d)
KVK     :  0.62  : Vitamin K status for half-maximal PIVKA-II formation
WINR    :  0.42  : INR increment per unit PIVKA-II above normal
PIVKA0  :  0.278 : PIVKA-II in an uncomplicated pregnancy

// ---- fetal lung maturity -------------------------------------------
KS_SURF :  1.0   : Surfactant index formation rate constant
KD_SURF :  1.0   : Surfactant index turnover rate constant (1/d)
EBET    :  1.30  : Betamethasone maximal effect on lung maturity
KBET    :  1.8   : Betamethasone concentration, half-maximal (mg/L)

// ---- drug disposition ----------------------------------------------
KAUD    :  8.0   : UDCA absorption rate constant from the gut depot (1/d)
KTRUD   :  3.2   : UDCA gut depot transit to colon (1/d)
FUD     :  0.42  : UDCA fraction absorbed and presented to the liver
KA_RIF  : 15.0   : Rifampicin absorption rate constant (1/d)
F_RIF   :  0.90  : Rifampicin oral bioavailability
CL_RIF  : 170.0  : Rifampicin clearance (L/d)
VC_RIF  : 45.0   : Rifampicin central volume (L)
Q_RIF   : 30.0   : Rifampicin intercompartmental clearance (L/d)
VP_RIF  : 20.0   : Rifampicin peripheral volume (L)
EMAX3A  :  4.2   : Maximal CYP3A4 induction by rifampicin
EC503A  :  1.6   : Rifampicin concentration for half-maximal induction (mg/L)
KS_3A   :  0.75  : CYP3A4 synthesis rate constant
KD_3A   :  0.75  : CYP3A4 degradation rate constant (1/d)
V6A     : 260.0  : CYP3A4 maximal 6a-hydroxylation capacity (umol/d)
F6A1    :  0.15  : 6a-hydroxylation selectivity, cholate
F6A2    :  1.00  : 6a-hydroxylation selectivity, chenodeoxycholate
F6A3    :  0.30  : 6a-hydroxylation selectivity, deoxycholate
F6A4    :  0.90  : 6a-hydroxylation selectivity, lithocholate
F6A5    :  0.05  : 6a-hydroxylation selectivity, ursodeoxycholate
CLHCA   : 60.0   : Hyocholate clearance (L/d)
KCLRCH  :  4.5   : Cholestyramine luminal transit rate constant (1/d)
CL_SAM  : 48.0   : SAMe clearance (L/d)
V_SAM   : 22.0   : SAMe volume of distribution (L)
KODEV   :  6.0   : IBAT inhibitor luminal transit rate constant (1/d)
CL_NTX  : 190.0  : Naltrexone clearance (L/d)
V_NTX   : 1350.  : Naltrexone volume of distribution (L)
EC50NTX :  1.3   : Naltrexone concentration for 50% MOR occupancy (ug/L)
CL_AH   : 95.0   : Antihistamine clearance (L/d)
V_AH    : 340.0  : Antihistamine volume of distribution (L)
EC50AH  :  9.0   : Antihistamine concentration for 50% effect (ug/L)
CL_BET  : 13.0   : Betamethasone clearance (L/d)
V_BET   : 42.0   : Betamethasone volume of distribution (L)

$CMT @annotated
// ---- 5 bile acid species x 7 compartments (umol) --------------------
HLCA  : Hepatocyte cholate
HLCD  : Hepatocyte chenodeoxycholate
HLDC  : Hepatocyte deoxycholate
HLLC  : Hepatocyte lithocholate
HLUD  : Hepatocyte ursodeoxycholate
BICA  : Biliary cholate
BICD  : Biliary chenodeoxycholate
BIDC  : Biliary deoxycholate
BILC  : Biliary lithocholate
BIUD  : Biliary ursodeoxycholate
ILCA  : Small-intestinal cholate
ILCD  : Small-intestinal chenodeoxycholate
ILDC  : Small-intestinal deoxycholate
ILLC  : Small-intestinal lithocholate
ILUD  : Small-intestinal ursodeoxycholate
COCA  : Colonic cholate
COCD  : Colonic chenodeoxycholate
CODC  : Colonic deoxycholate
COLC  : Colonic lithocholate
COUD  : Colonic ursodeoxycholate
SPCA  : Maternal plasma cholate
SPCD  : Maternal plasma chenodeoxycholate
SPDC  : Maternal plasma deoxycholate
SPLC  : Maternal plasma lithocholate
SPUD  : Maternal plasma ursodeoxycholate
FPCA  : Fetal plasma cholate
FPCD  : Fetal plasma chenodeoxycholate
FPDC  : Fetal plasma deoxycholate
FPLC  : Fetal plasma lithocholate
FPUD  : Fetal plasma ursodeoxycholate
FLCA  : Fetal liver cholate
FLCD  : Fetal liver chenodeoxycholate
FLDC  : Fetal liver deoxycholate
FLLC  : Fetal liver lithocholate
FLUD  : Fetal liver ursodeoxycholate
// ---- global states --------------------------------------------------
AF    : Amniotic fluid bile acid (umol)
MEC   : Meconium bile acid, cumulative (umol)
CYP7A1: Cholesterol 7a-hydroxylase (relative)
SHP   : Small heterodimer partner (relative)
FGF19 : Ileal FGF19 (relative)
E2    : Estradiol (nmol/L)
P4    : Progesterone (nmol/L)
P4S   : Progesterone sulfates (nmol/L)
E2G   : Estradiol-17b-glucuronide (nmol/L)
ATX   : Autotaxin (relative, 1 = non-pregnant)
LPA   : Lysophosphatidic acid (relative)
ITCHC : Central itch state (relative)
SLEEPD: Sleep disturbance index
ROS   : Hepatocyte oxidative stress (relative)
ALT   : Alanine aminotransferase (U/L)
GJ    : Fetal cardiac connexin-43 coupling (relative)
FCA   : Fetal cardiomyocyte calcium index (relative)
VASO  : Chorionic vasoconstriction index
TROPH : Trophoblast damage index
HYP   : Fetal hypoxia index
OTR   : Myometrial oxytocin receptor density (relative)
UTA   : Uterine activity index
VITK  : Maternal vitamin K status (relative)
PIVKA : Under-carboxylated prothrombin (relative)
HSB   : Cumulative stillbirth hazard
HPT   : Cumulative spontaneous preterm birth hazard
HMEC  : Cumulative meconium passage hazard
UDDEP : UDCA gut depot (umol)
RIFDEP: Rifampicin gut depot (mg)
RIFC  : Rifampicin central (mg)
RIFP  : Rifampicin peripheral (mg)
CYP3A : CYP3A4 induction (relative)
HCA   : Hyocholate / hyodeoxycholate pool (umol)
CHOLL : Cholestyramine luminal (g)
SAMC  : SAMe central (mg)
ODEV  : IBAT inhibitor luminal (umol)
NTXC  : Naltrexone central (mg)
AHC   : Antihistamine central (mg)
BETC  : Betamethasone central (mg)
SURF  : Fetal surfactant / lung maturity (relative)
SULT  : SULT2A1 induction (relative)
MRP4  : Basolateral MRP3/MRP4/OSTa-b (relative)
TDOSE : Cumulative UDCA dose (umol, bookkeeping)

$INIT @annotated
// Converged pre-third-trimester state of icp_reference_model.py after a
// 140-day burn-in at frozen GA 20 on the wild-type susceptibility
// vector.  Starting an ICP genotype from here is deliberate: at GA 20
// the disease has not manifested, and the pool re-equilibrates within
// about two weeks of simulated time.
HLCA  : 154.139   : umol
HLCD  : 107.812   : umol
HLDC  :  35.4695  : umol
HLLC  :   2.8481  : umol
HLUD  :   0.0     : umol
BICA  : 1091.18   : umol
BICD  :  944.938  : umol
BIDC  :  310.880  : umol
BILC  :   37.4442 : umol
BIUD  :    0.0    : umol
ILCA  :  522.028  : umol
ILCD  :  393.169  : umol
ILDC  :  129.228  : umol
ILLC  :   12.0445 : umol
ILUD  :    0.0    : umol
COCA  :   99.209  : umol
COCD  :   68.8046 : umol
CODC  :  102.129  : umol
COLC  :   55.1905 : umol
COUD  :    0.0    : umol
SPCA  :   20.9479 : umol
SPCD  :   16.0279 : umol
SPDC  :    5.47396: umol
SPLC  :    0.62179: umol
SPUD  :    0.0    : umol
FPCA  :    0.158479 : umol
FPCD  :    0.0855994: umol
FPDC  :    0.00448707 : umol
FPLC  :    0.000442971 : umol
FPUD  :    0.0    : umol
FLCA  :    3.21932: umol
FLCD  :    2.09775: umol
FLDC  :    0.0108424 : umol
FLLC  :    0.00107033 : umol
FLUD  :    0.0    : umol
AF    :    0.0269183 : umol
MEC   :    0.0    : umol
CYP7A1:    0.739536 : relative
SHP   :    0.0894186: relative
FGF19 :    0.413045 : relative
E2    :   20.0018 : nmol/L
P4    :  200.016  : nmol/L
P4S   :   32.0025 : nmol/L
E2G   :    1.19366: nmol/L
ATX   :    1.47258: relative
LPA   :    1.47258: relative
ITCHC :    0.0605206 : relative
SLEEPD:    0.0605206 : relative
ROS   :    0.00204795: relative
ALT   :   20.3686 : U/L
GJ    :    0.996051 : relative
FCA   :    1.05689: relative
VASO  :    0.0255278 : relative
TROPH :    0.0186958 : relative
HYP   :    0.0376801 : relative
OTR   :    1.07093: relative
UTA   :    0.0709295 : relative
VITK  :    0.963726 : relative
PIVKA :    0.292727 : relative
HSB   :    0.0     : cumulative hazard
HPT   :    0.0     : cumulative hazard
HMEC  :    0.0     : cumulative hazard
UDDEP :    0.0     : umol
RIFDEP:    0.0     : mg
RIFC  :    0.0     : mg
RIFP  :    0.0     : mg
CYP3A :    1.0     : relative
HCA   :   15.7773  : umol
CHOLL :    0.0     : g
SAMC  :    0.0     : mg
ODEV  :    0.0     : umol
NTXC  :    0.0     : mg
AHC   :    0.0     : mg
BETC  :    0.0     : mg
SURF  :    1.0     : relative
SULT  :    1.0     : relative
MRP4  :    1.05486 : relative
TDOSE :    0.0     : umol

$GLOBAL
#define NSP 5
// Cytotoxicity / detergency weights.  Ordering and spacing follow the
// hydrophobicity (critical micellar concentration) series and the fetal
// cardiomyocyte data: LCA > DCA > CDCA >> CA >>> UDCA.  UDCA is not zero
// (taurine-conjugated UDCA is weakly detergent) but it is two orders
// below LCA, and that gap is why UDCA can lower the assay by a third
// without lowering the toxic load by a third.
#define W_CA   0.20
#define W_CDCA 0.60
#define W_DCA  0.72
#define W_LCA  1.00
#define W_UDCA 0.02
// FXR agonist potency; CDCA is the physiological ligand, UDCA is inert
#define FX_CA   0.30
#define FX_CDCA 1.00
#define FX_DCA  0.60
#define FX_LCA  0.50
#define FX_UDCA 0.02

$MAIN
double GA = GA0 + SOLVERTIME / 7.0;
double PL = (SOLVERTIME < TDEL) ? 1.0 : 0.0;   // placental function switch

$ODE
// ====================================================================
// gestational covariates
// ====================================================================
double ga = GA0 + SOLVERTIME / 7.0;
double pl = (SOLVERTIME < TDEL) ? 1.0 : 0.0;

// Hadlock-style log-quadratic fetal weight (kg)
double fw = exp(0.578 + 0.332 * ga - 0.00354 * ga * ga) / 1000.0;
double VFP = (0.25 * fw > 0.05) ? 0.25 * fw : 0.05;
double VFL = (0.035 * fw > 0.01) ? 0.035 * fw : 0.01;

// ====================================================================
// concentrations
// ====================================================================
double CHL[NSP], CIL[NSP], CSP[NSP], CFP[NSP];
double HLv[NSP] = {HLCA, HLCD, HLDC, HLLC, HLUD};
double BIv[NSP] = {BICA, BICD, BIDC, BILC, BIUD};
double ILv[NSP] = {ILCA, ILCD, ILDC, ILLC, ILUD};
double COv[NSP] = {COCA, COCD, CODC, COLC, COUD};
double SPv[NSP] = {SPCA, SPCD, SPDC, SPLC, SPUD};
double FPv[NSP] = {FPCA, FPCD, FPDC, FPLC, FPUD};
double FLv[NSP] = {FLCA, FLCD, FLDC, FLLC, FLUD};
int i;
for (i = 0; i < NSP; i++) {
  if (HLv[i] < 0.0) HLv[i] = 0.0;
  if (BIv[i] < 0.0) BIv[i] = 0.0;
  if (ILv[i] < 0.0) ILv[i] = 0.0;
  if (COv[i] < 0.0) COv[i] = 0.0;
  if (SPv[i] < 0.0) SPv[i] = 0.0;
  if (FPv[i] < 0.0) FPv[i] = 0.0;
  if (FLv[i] < 0.0) FLv[i] = 0.0;
  CHL[i] = HLv[i] / VHL;
  CIL[i] = ILv[i] / VIL;
  CSP[i] = SPv[i] / VSP;
  CFP[i] = FPv[i] / VFP;
}
double KMB[NSP] = {KMB1, KMB2, KMB3, KMB4, 1e9};
double SBL[NSP] = {SBL1, SBL2, SBL3, SBL4, SBL5};
double EUP[NSP] = {EUP1, EUP2, EUP3, EUP4, EUP5};
double KMA[NSP] = {KMA1, KMA2, KMA3, KMA4, KMA5};
double KPAS[NSP] = {KPAS1, KPAS2, KPAS3, KPAS4, KPAS5};
double FBI[NSP] = {FBI1, FBI2, FBI3, FBI4, FBI5};
double KCOL[NSP] = {KCOL1, KCOL2, KCOL3, KCOL4, KCOL5};
double CLR[NSP] = {CLR1, CLR2, CLR3, CLR4, CLR5};
double KMP[NSP] = {KMP1, KMP2, KMP3, KMP4, KMP5};
double PSv[NSP] = {PS1, PS2, PS3, PS4, PS5};
double F6A[NSP] = {F6A1, F6A2, F6A3, F6A4, F6A5};
double WTOX[NSP] = {W_CA, W_CDCA, W_DCA, W_LCA, W_UDCA};
double WFXR[NSP] = {FX_CA, FX_CDCA, FX_DCA, FX_LCA, FX_UDCA};
double FSFv[NSP] = {FSF1, FSF2, 0.0, 0.0, 0.0};

// ====================================================================
// sex steroids: the trigger and the clock
// ====================================================================
double E2tgt = E2BASE * exp(E2K * (ga - 20.0)) * TWIN;
double P4tgt = P4BASE * exp(P4K * (ga - 20.0)) * TWIN;
double cyp3a = (CYP3A > 1e-6) ? CYP3A : 1e-6;
double sult  = (SULT  > 1e-6) ? SULT  : 1e-6;

dxdt_E2 = KOUTE2 * (E2tgt * pl - E2);
dxdt_P4 = KOUTP4 * (P4tgt * pl - P4);
// Rifampicin co-induces UGT1A1, so it makes MORE of the BSEP cis-
// inhibitor while suppressing autotaxin.  The tension is kept rather
// than hidden; at 300 mg bid the autotaxin effect dominates.
dxdt_E2G = KFE2G * E2 * (1.0 + EUGT * (cyp3a - 1.0)) - KOUTE2G * E2G;
dxdt_P4S = KFP4S * P4 * GSULT * sult - KOUTP4S * P4S;

double SIP = E2G / KI_E2G + P4S / KI_P4S;   // steroid inhibitory pressure
double IBSEP = 1.0 / (1.0 + SIP);

// ====================================================================
// hepatic synthesis and its negative feedback
// ====================================================================
double FXRW = 0.0, CHLtot = 0.0, HSTRESS = 0.0, ILBA = 0.0, ILtot = 0.0;
for (i = 0; i < NSP; i++) {
  FXRW    += WFXR[i] * CHL[i];
  CHLtot  += CHL[i];
  HSTRESS += WTOX[i] * CHL[i];
  ILBA    += WFXR[i] * CIL[i];
  ILtot   += CIL[i];
}
double FXRH = GFXR * FXRW / (EC50FXR + FXRW);
double FXRI = GFXR * ILBA / (EC50FXI + ILBA);
double shp = (SHP > 0.0) ? SHP : 0.0;
double fgf = (FGF19 > 0.0) ? FGF19 : 0.0;
double c7  = (CYP7A1 > 0.0) ? CYP7A1 : 0.0;

dxdt_SHP   = KS_SHP * FXRH - KD_SHP * shp;
dxdt_FGF19 = KS_F19 * FXRI - KD_F19 * fgf;
dxdt_CYP7A1 = KS_C7 / (1.0 + pow(shp / KSHP, 2.0) + fgf / KF19) - KD_C7 * c7;

double SYN = VSYN * c7;
double Ssyn[NSP] = {SYN * FRCA, SYN * (1.0 - FRCA), 0.0, 0.0, 0.0};

// ====================================================================
// canalicular export: BSEP for the four endogenous species, MRP2 for UDCA
// ====================================================================
double samc = (SAMC > 0.0) ? SAMC : 0.0;
double TRAF  = 1.0 + ETRAF * CHL[4] / (KTRAF + CHL[4]);
double SAMEF = 1.0 + ESAM * samc / (KSAM + samc);
double denB = 1.0;
for (i = 0; i < 4; i++) denB += CHL[i] / KMB[i];
double gB = GBSEP * GFIC1;
double JB[NSP];
for (i = 0; i < 4; i++)
  JB[i] = VBSEP * gB * IBSEP * TRAF * SAMEF * (CHL[i] / KMB[i]) / denB;
JB[4] = VMRP2 * (1.0 + EMRP2P * (cyp3a - 1.0)) * CHL[4] / (KMMRP2 + CHL[4]);

// ---- lithocholate detoxification (SULT2A1) --------------------------
double JLCAS = VLCAS * sult * CHL[3] / (KMLCAS + CHL[3]);

// ---- 6a-hydroxylation (CYP3A4) -> hyocholate ------------------------
double J6A[NSP], J6Atot = 0.0;
for (i = 0; i < NSP; i++) {
  J6A[i] = V6A * cyp3a * F6A[i] * CHL[i] / 100.0;
  J6Atot += J6A[i];
}

// ---- basolateral escape valve ---------------------------------------
double mrp4 = (MRP4 > 0.0) ? MRP4 : 0.0;
double rm4 = pow(FXRH / KM4, HM4);
dxdt_MRP4 = KS_M4 * (1.0 + EM4 * rm4 / (1.0 + rm4) + EM4P * (cyp3a - 1.0))
            - KD_M4 * mrp4;
double JBL[NSP];
for (i = 0; i < NSP; i++) JBL[i] = KBL * mrp4 * SBL[i] * CHL[i];

// ====================================================================
// intestinal handling
// ====================================================================
double choll = (CHOLL > 0.0) ? CHOLL : 0.0;
double odev  = (ODEV  > 0.0) ? ODEV  : 0.0;
double denA = 1.0 + odev / KIODEV;
for (i = 0; i < NSP; i++) denA += CIL[i] / KMA[i];
double JASBT[NSP], JSEQ[NSP];
double fseq = KSEQ * choll / (KCHOL + choll);
for (i = 0; i < NSP; i++) {
  JASBT[i] = VASBT * (CIL[i] / KMA[i]) / denA + KPAS[i] * ILv[i];
  JSEQ[i]  = fseq * FBI[i] * ILv[i];
}

// ---- colonic microbial transformation -------------------------------
double rCA = KDH_CA * COv[0];
double rCD = KDH_CD * COv[1];
double rUD = KDH_UD * COv[4];
double JCOL[NSP];
for (i = 0; i < NSP; i++) JCOL[i] = KCOL[i] * COv[i];

// ====================================================================
// hepatic first pass and systemic re-uptake
// ====================================================================
double ABS_[NSP], ABStot = 0.0;
for (i = 0; i < NSP; i++) { ABS_[i] = JASBT[i] + JCOL[i]; ABStot += ABS_[i]; }
double CRIF = ((RIFC > 0.0) ? RIFC : 0.0) / VC_RIF;
// Sinusoidal uptake efficiency.  NTCP/OATP1B1 are transcriptionally
// DOWN-regulated once the hepatocyte is loaded, so the term is a
// threshold in CHLtot, not a proportionality: first pass is essentially
// complete in health and degrades only once cholestasis is established.
// This is what turns a linear fall in BSEP capacity into a supralinear
// rise in maternal plasma bile acid.
double updep = 1.0 + ABStot / JHALF + CRIF / KIRIFO
               + pow(CHLtot / KHLSAT, HHLSAT);
double JU[NSP], JSPILL[NSP], JUS[NSP];
for (i = 0; i < NSP; i++) {
  double E = EUP[i] / updep;
  JU[i]     = E * ABS_[i];
  JSPILL[i] = (1.0 - E) * ABS_[i];
  JUS[i]    = CLHSYS / updep * CSP[i];
}

// ---- renal ----------------------------------------------------------
double JREN[NSP];
for (i = 0; i < NSP; i++)
  JREN[i] = CLR[i] * (1.0 + ESULT * (sult - 1.0)) * CSP[i];
dxdt_SULT = KS_S * (1.0 + EPXR_S * (cyp3a - 1.0)) - KD_S * sult;

// ====================================================================
// placental transfer
// ====================================================================
double denP = 1.0;
for (i = 0; i < NSP; i++) denP += CFP[i] / KMP[i] + TRANSIN * CSP[i] / KMP[i];
double JF2M[NSP], JPASS[NSP];
for (i = 0; i < NSP; i++) {
  JF2M[i]  = VP * GP * pl * (CFP[i] / KMP[i]) / denP;
  // The non-carrier component of placental transfer is DIFFUSIVE and
  // therefore bidirectional; the gradient term is what bounds the fetal
  // pool once the carrier saturates.
  JPASS[i] = PSv[i] * pl * (CSP[i] - CFP[i]);
}

// ---- fetal liver, urine, amniotic fluid, meconium -------------------
double JFLIN[NSP], JFLOUT[NSP], JFBILE[NSP];
double FPtot = 0.0, JFBtot = 0.0, FGUT = 0.0;
for (i = 0; i < NSP; i++) {
  JFLIN[i]  = CLFL * CFP[i];
  JFLOUT[i] = KFLO * FLv[i];
  JFBILE[i] = KFB * FLv[i];
  FPtot  += CFP[i];
  JFBtot += JFBILE[i];
  FGUT   += WTOX[i] * FLv[i];
}
FGUT = FGUT / VFL;
double af = (AF > 0.0) ? AF : 0.0;
double JUR = KUR * FPtot * VFP;
dxdt_AF  = JUR - (KSW + KIMA) * af;
dxdt_MEC = JFBtot * pl;

// ====================================================================
// species mass balances
// ====================================================================
double dHL[NSP], dBI[NSP], dIL[NSP], dCO[NSP], dSP[NSP], dFP[NSP], dFL[NSP];
for (i = 0; i < NSP; i++) {
  dHL[i] = Ssyn[i] + JU[i] + JUS[i] - JB[i] - JBL[i] - J6A[i];
  dBI[i] = JB[i] - KBILE * BIv[i];
  dIL[i] = KBILE * BIv[i] - JASBT[i] - JSEQ[i] - KTR * ILv[i];
  dCO[i] = KTR * ILv[i] - JCOL[i] - KFEC * COv[i];
  dSP[i] = JSPILL[i] + JBL[i] - JUS[i] - JREN[i] + JF2M[i] - JPASS[i];
  dFP[i] = JPASS[i] - JF2M[i] + JFLOUT[i] - JFLIN[i] - KUR * CFP[i] * VFP;
  dFL[i] = SF * FSFv[i] * pl + JFLIN[i] - JFLOUT[i] - JFBILE[i];
}
dHL[3] -= JLCAS;                       // hepatic LCA sulfation
dCO[0] -= rCA;  dCO[2] += rCA;         // CA  -> DCA
dCO[1] -= rCD;  dCO[3] += rCD;         // CDCA -> LCA
dCO[4] -= rUD;  dCO[3] += rUD;         // UDCA -> LCA (the UDCA liability)

// ---- oral UDCA: absorbed to liver, unabsorbed to colon --------------
double uddep = (UDDEP > 0.0) ? UDDEP : 0.0;
dxdt_UDDEP = -(KAUD + KTRUD) * uddep;
dHL[4] += KAUD * uddep * FUD;
dSP[4] += KAUD * uddep * FUD * (1.0 - FUD) * 0.25;
dCO[4] += KTRUD * uddep;
dxdt_TDOSE = 0.0;

dxdt_HLCA = dHL[0]; dxdt_HLCD = dHL[1]; dxdt_HLDC = dHL[2];
dxdt_HLLC = dHL[3]; dxdt_HLUD = dHL[4];
dxdt_BICA = dBI[0]; dxdt_BICD = dBI[1]; dxdt_BIDC = dBI[2];
dxdt_BILC = dBI[3]; dxdt_BIUD = dBI[4];
dxdt_ILCA = dIL[0]; dxdt_ILCD = dIL[1]; dxdt_ILDC = dIL[2];
dxdt_ILLC = dIL[3]; dxdt_ILUD = dIL[4];
dxdt_COCA = dCO[0]; dxdt_COCD = dCO[1]; dxdt_CODC = dCO[2];
dxdt_COLC = dCO[3]; dxdt_COUD = dCO[4];
dxdt_SPCA = dSP[0]; dxdt_SPCD = dSP[1]; dxdt_SPDC = dSP[2];
dxdt_SPLC = dSP[3]; dxdt_SPUD = dSP[4];
dxdt_FPCA = dFP[0]; dxdt_FPCD = dFP[1]; dxdt_FPDC = dFP[2];
dxdt_FPLC = dFP[3]; dxdt_FPUD = dFP[4];
dxdt_FLCA = dFL[0]; dxdt_FLCD = dFL[1]; dxdt_FLDC = dFL[2];
dxdt_FLLC = dFL[3]; dxdt_FLUD = dFL[4];

dxdt_HCA = J6Atot + JLCAS - CLHCA * ((HCA > 0.0) ? HCA : 0.0) / VSP;

// ====================================================================
// derived exposure indices
// ====================================================================
double FCL = 0.0, MCL = 0.0;   // fetal / maternal cytotoxic load
for (i = 0; i < NSP; i++) { FCL += WTOX[i] * CFP[i]; MCL += WTOX[i] * CSP[i]; }

// ====================================================================
// fetal myocardium: where the 100 umol/L threshold lives
// ====================================================================
double gj  = (GJ  > 0.0) ? GJ  : 0.0;
double fca = (FCA > 0.0) ? FCA : 0.0;
double rgj = (FCL > 0.0) ? pow(FCL / IC50GJ, HGJ) : 0.0;
double GJinh = IMAXGJ * rgj / (1.0 + rgj);
dxdt_GJ  = KS_GJ * (1.0 - GJinh) - KD_GJ * gj;
dxdt_FCA = KS_CA * (1.0 + ECA * FCL / (KCA + FCL)) - KD_CA * fca;
double ARRI = (1.0 - gj > 0.0 ? 1.0 - gj : 0.0) * fca;

// ====================================================================
// placental vasculature and fetal hypoxia
// ====================================================================
double vaso  = (VASO  > 0.0) ? VASO  : 0.0;
double troph = (TROPH > 0.0) ? TROPH : 0.0;
double hyp   = (HYP   > 0.0) ? HYP   : 0.0;
dxdt_VASO  = KS_V * EV * FCL / (KV + FCL) - KD_V * vaso;
dxdt_TROPH = KS_T * FCL / (KT + FCL) - KD_T * troph;
dxdt_HYP   = KS_H * (vaso + WT * troph) - KD_H * hyp;

// ====================================================================
// outcome hazards
// ====================================================================
// Gating: "stillbirth" is only defined from 24 weeks, so the hazard is
// not allowed to accrue before then, or the fitted background hazard
// silently absorbs four weeks of pre-viable gestation and stops meaning
// what the epidemiology means by it.
double gate_sb  = (ga >= 24.0) ? 1.0 : 0.0;
double gate_pt  = (ga >= 24.0 && ga < 37.0) ? 1.0 : 0.0;
// Meconium passage in utero is a near-term phenomenon (gut motility and
// TGR5 expression both mature late).
double gate_mec = 1.0 / (1.0 + exp(-1.2 * (ga - 35.0)));

double otr = (OTR > 0.0) ? OTR : 0.0;
double uta = (UTA > 0.0) ? UTA : 0.0;
dxdt_OTR = KS_O * (1.0 + EO * MCL / (KO + MCL)) - KD_O * otr;
dxdt_UTA = KS_U * ((otr - 1.0 > 0.0) ? otr - 1.0 : 0.0) - KD_U * uta;

double hSB = HSB0 + HSBSC * pow(ARRI, HN) * (1.0 + KHY * hyp);
double hPT = HPT0 * (1.0 + KOT * ((otr - 1.0 > 0.0) ? otr - 1.0 : 0.0))
                  * (1.0 + KUT * uta);
double hMEC = (HM0 + HMSC * FCL / (KMG + FCL) * (1.0 + hyp)) * gate_mec;
dxdt_HSB  = hSB  * pl * gate_sb;
dxdt_HPT  = hPT  * pl * gate_pt;
dxdt_HMEC = hMEC * pl;

// ====================================================================
// pruritus: the second axis
// ====================================================================
double atx   = (ATX   > 0.0) ? ATX   : 0.0;
double lpa   = (LPA   > 0.0) ? LPA   : 0.0;
double itchc = (ITCHC > 0.0) ? ITCHC : 0.0;
double rifatx = 1.0 / (1.0 + ERIFA * (cyp3a - 1.0));
double rp4s = (P4S > 0.0) ? pow(P4S / KA_P4S, HP4S) : 0.0;
dxdt_ATX = KS_A * (1.0 + EA_E2 * E2 / (KA_E2 + E2)
                   + EA_P4S * rp4s / (1.0 + rp4s)
                   + EA_CH * MCL / (KA_CH + MCL)) * rifatx - KD_A * atx;
dxdt_LPA = KS_L * atx - KD_L * lpa;

double ntxc = ((NTXC > 0.0) ? NTXC : 0.0) / V_NTX * 1000.0;
double NTXocc = ntxc / (EC50NTX + ntxc);
double ahc = ((AHC > 0.0) ? AHC : 0.0) / V_AH * 1000.0;
double AHocc = ahc / (EC50AH + ahc);
double rlp = (lpa > 0.0) ? pow(lpa / ILP50, HLP) : 0.0;
double itch_raw = IMAXI * rlp / (1.0 + rlp) + WBA * MCL / (KBA + MCL);
double sens = 1.0 + KSENS * itchc / (KSC + itchc);
dxdt_ITCHC = KS_I * itch_raw * sens * (1.0 - ENTX * NTXocc)
             * (1.0 - EAH * AHocc) - KD_I * itchc;
dxdt_SLEEPD = KS_SL * itchc - KD_SL * ((SLEEPD > 0.0) ? SLEEPD : 0.0);

// ====================================================================
// maternal hepatocellular injury
// ====================================================================
double PCFLUX = VPC * GMDR3;
// UDCA is excluded from the biliary bile-salt : phosphatidylcholine
// ratio: it is not detergent at biliary concentrations, and this is
// precisely why it protects the ABCB4-deficient canaliculus.  Counting
// it made the drug look as if it worsened the lesion it treats.
double JBendo = JB[0] + JB[1] + JB[2] + JB[3];
double BSPC = JBendo / ((PCFLUX > 1.0) ? PCFLUX : 1.0);
double ros = (ROS > 0.0) ? ROS : 0.0;
double prot_ud = 1.0 - EUDROS * CHL[4] / (KUDROS + CHL[4]);
double rr = (HSTRESS > 0.0) ? pow(HSTRESS / KR, HR) : 0.0;
double bsx = (BSPC / BSPC0 - 1.0 > 0.0) ? BSPC / BSPC0 - 1.0 : 0.0;
// Rifampicin carries its own recognised hepatotoxic signal, so it is not
// allowed to drive ALT to normal purely by depleting CDCA.
dxdt_ROS = KS_R * (rr / (1.0 + rr) + WBSPC * bsx
                   + WRIFROS * CRIF / (KRIFROS + CRIF)) * prot_ud
           - KD_R * ros * (1.0 + ESAM_R * samc / (KSAM + samc));
dxdt_ALT = KS_ALT * (1.0 + EALT * ros) - KD_ALT * ((ALT > 0.0) ? ALT : 0.0);

// ====================================================================
// fat-soluble vitamins and coagulation
// ====================================================================
double MICELLE = ILtot / (MIC50 + ILtot);
double vitk = (VITK > 1e-6) ? VITK : 1e-6;
double vkin = (SOLVERTIME >= VKSTART && SOLVERTIME < TDEL) ? VKDOSE : 0.0;
dxdt_VITK = KIN_K * (DIETK + vkin) * MICELLE / 0.5 - KD_K * vitk;
dxdt_PIVKA = KS_PI / (1.0 + pow(vitk / KVK, 2.0))
             - KD_PI * ((PIVKA > 0.0) ? PIVKA : 0.0);

// ====================================================================
// fetal lung maturity
// ====================================================================
double betc = (BETC > 0.0) ? BETC : 0.0;
double cbet = betc / V_BET;
dxdt_SURF = KS_SURF * (1.0 + EBET * cbet / (KBET + cbet))
            - KD_SURF * ((SURF > 0.0) ? SURF : 0.0);

// ====================================================================
// drug disposition
// ====================================================================
double rifdep = (RIFDEP > 0.0) ? RIFDEP : 0.0;
dxdt_RIFDEP = -KA_RIF * rifdep;
dxdt_RIFC = KA_RIF * rifdep * F_RIF - CL_RIF / VC_RIF * ((RIFC > 0.0) ? RIFC : 0.0)
            - Q_RIF / VC_RIF * ((RIFC > 0.0) ? RIFC : 0.0)
            + Q_RIF / VP_RIF * ((RIFP > 0.0) ? RIFP : 0.0);
dxdt_RIFP = Q_RIF / VC_RIF * ((RIFC > 0.0) ? RIFC : 0.0)
            - Q_RIF / VP_RIF * ((RIFP > 0.0) ? RIFP : 0.0);
dxdt_CYP3A = KS_3A * (1.0 + EMAX3A * CRIF / (EC503A + CRIF)) - KD_3A * cyp3a;

dxdt_CHOLL = -KCLRCH * choll;
dxdt_SAMC  = -CL_SAM / V_SAM * samc;
dxdt_ODEV  = -KODEV * odev;
dxdt_NTXC  = -CL_NTX / V_NTX * ((NTXC > 0.0) ? NTXC : 0.0);
dxdt_AHC   = -CL_AH / V_AH * ((AHC > 0.0) ? AHC : 0.0);
dxdt_BETC  = -CL_BET / V_BET * betc;

$TABLE
double VFPt = (0.25 * exp(0.578 + 0.332 * GA - 0.00354 * GA * GA) / 1000.0);
if (VFPt < 0.05) VFPt = 0.05;
double VFLt = (0.035 * exp(0.578 + 0.332 * GA - 0.00354 * GA * GA) / 1000.0);
if (VFLt < 0.01) VFLt = 0.01;

double cCA = SPCA / VSP, cCD = SPCD / VSP, cDC = SPDC / VSP;
double cLC = SPLC / VSP, cUD = SPUD / VSP, cHCA = HCA / VSP;
double fCA = FPCA / VFPt, fCD = FPCD / VFPt, fDC = FPDC / VFPt;
double fLC = FPLC / VFPt, fUD = FPUD / VFPt;

// The routine clinical assay is 3a-hydroxysteroid dehydrogenase based:
// it counts ursodeoxycholate and hyocholate along with the endogenous
// species.  Reporting both the assay value and the endogenous fraction
// is the point of carrying species at all.
capture TBA      = cCA + cCD + cDC + cLC + cUD + cHCA;
capture TBA_ENDO = cCA + cCD + cDC + cLC;
capture UDCA_P   = cUD;
capture HCA_P    = cHCA;
capture CORD     = fCA + fCD + fDC + fLC + fUD;
capture FCLo = W_CA*fCA + W_CDCA*fCD + W_DCA*fDC + W_LCA*fLC + W_UDCA*fUD;
capture MCLo = W_CA*cCA + W_CDCA*cCD + W_DCA*cDC + W_LCA*cLC + W_UDCA*cUD;
capture WBAR     = (CORD > 1e-9) ? FCLo / CORD : 0.0;
capture ARRIo    = ((1.0 - GJ) > 0.0 ? (1.0 - GJ) : 0.0) * FCA;
capture VAS      = (VSCALE * ITCHC < 10.0) ? VSCALE * ITCHC : 10.0;
capture INR      = 1.0 + WINR * ((PIVKA - PIVKA0 > 0.0) ? PIVKA - PIVKA0 : 0.0);
capture SBRISK   = 1.0 - exp(-((HSB > 0.0) ? HSB : 0.0));
capture PTBRISK  = 1.0 - exp(-((HPT > 0.0) ? HPT : 0.0));
capture MECRISK  = 1.0 - exp(-((HMEC > 0.0) ? HMEC : 0.0));
// Respiratory distress and any-neonatal-unit-admission are functions of
// gestational age at delivery, NOT of bile acid.  That is exactly why
// they cannot be moved by a bile-acid-lowering drug, and it is the whole
// arithmetic of the PITCHES null result.
capture RDSRISK  = (0.0028 + 1.0 / (1.0 + exp(1.15 * (GA - 33.6))))
                   / (1.0 + 1.55 * ((SURF - 1.0 > 0.0) ? SURF - 1.0 : 0.0));
capture NICURISK = (0.020 + 1.0 / (1.0 + exp(0.82 * (GA - 35.3))))
                   / (1.0 + 0.85 * ((SURF - 1.0 > 0.0) ? SURF - 1.0 : 0.0))
                   + 0.16 * MECRISK + 0.22 * ((HYP / 1.2 < 1.0) ? HYP / 1.2 : 1.0);
capture GAw      = GA;
capture SEVERE_FLAG = (TBA >= 40.0) ? 1.0 : 0.0;
capture VSEV_FLAG   = (TBA >= 100.0) ? 1.0 : 0.0;

## =====================================================================
##  SCENARIOS
##  ---------
##  library(mrgsolve); library(dplyr); library(tidyr)
##  mod <- mread("icp_mrgsolve_model", ".")
##
##  Susceptibility vectors.  Severity strata are DEFINED by the bile acid
##  band they produce, exactly as the literature defines them (>=10
##  diagnostic, >=40 severe, >=100 very severe), not by a genotype fixed
##  in advance.  These four were found by scanning GBSEP/GMDR3/GSULT.
##
##  WT     <- list(GBSEP = 1.000, GMDR3 = 1.00, GSULT = 1.00)
##  MILD1  <- list(GBSEP = 0.795, GMDR3 = 0.695, GSULT = 1.72)  # TBA ~16
##  MILD   <- list(GBSEP = 0.735, GMDR3 = 0.620, GSULT = 1.95)  # TBA ~26
##  SEV    <- list(GBSEP = 0.650, GMDR3 = 0.500, GSULT = 2.40)  # TBA ~58
##  VSEV   <- list(GBSEP = 0.520, GMDR3 = 0.350, GSULT = 3.10)  # TBA ~138
##
##  Dosing helpers.  TIME is days from GA 20+0, so GA 30+0 is day 70.
##  ga2d <- function(ga) (ga - 20) * 7
##
##  udca <- function(mg_per_dose = 500, from_ga = 30, to_ga = 39)
##    ev(amt = mg_per_dose / 392.57 * 1000, cmt = "UDDEP", ii = 0.5,
##       addl = (to_ga - from_ga) * 14 - 1, time = ga2d(from_ga))
##  rif  <- function(mg = 300, from_ga = 32, to_ga = 39)
##    ev(amt = mg, cmt = "RIFDEP", ii = 0.5,
##       addl = (to_ga - from_ga) * 14 - 1, time = ga2d(from_ga))
##  chol <- function(g = 4, from_ga = 32, to_ga = 39)
##    ev(amt = g, cmt = "CHOLL", ii = 0.25,
##       addl = (to_ga - from_ga) * 28 - 1, time = ga2d(from_ga))
##  same <- function(mg = 1000, from_ga = 32, to_ga = 39)
##    ev(amt = mg, cmt = "SAMC", ii = 1, addl = (to_ga - from_ga) * 7 - 1,
##       time = ga2d(from_ga))
##  ibat <- function(umol = 4, from_ga = 32, to_ga = 39)
##    ev(amt = umol, cmt = "ODEV", ii = 1, addl = (to_ga - from_ga) * 7 - 1,
##       time = ga2d(from_ga))
##  ntx  <- function(mg = 50, from_ga = 32, to_ga = 39)
##    ev(amt = mg, cmt = "NTXC", ii = 1, addl = (to_ga - from_ga) * 7 - 1,
##       time = ga2d(from_ga))
##  anti <- function(mg = 4, from_ga = 32, to_ga = 39)
##    ev(amt = mg, cmt = "AHC", ii = 1/3, addl = (to_ga - from_ga) * 21 - 1,
##       time = ga2d(from_ga))
##  beta <- function(from_ga = 35)
##    ev(amt = 12, cmt = "BETC", ii = 1, addl = 1, time = ga2d(from_ga))
##
##  run <- function(gen = list(), dose = NULL, del_ga = 39, extra = list())
##    mod %>% param(c(gen, extra, list(TDEL = ga2d(del_ga)))) %>%
##      (\(m) if (is.null(dose)) m else mrgsolve::ev(m, dose)) %>%
##      mrgsim(end = ga2d(del_ga) + 21, delta = 0.5) %>% as_tibble()
##
##  --- 1-4  natural history by stratum -------------------------------
##  s01 <- run(WT)                                  # normal pregnancy
##  s02 <- run(MILD1)                               # TBA <40 band
##  s03 <- run(MILD)                                # TBA ~26
##  s04 <- run(SEV)                                 # TBA 40-99 band
##  s05 <- run(VSEV)                                # TBA >=100 band
##
##  --- 6-11  UDCA dose and timing ------------------------------------
##  s06 <- run(SEV,  udca(500))                     # 500 mg bid from 30 wk
##  s07 <- run(SEV,  udca(750))                     # 750 mg bid
##  s08 <- run(SEV,  udca(1000))                    # 1000 mg bid
##  s09 <- run(SEV,  udca(500, from_ga = 34))       # late start
##  s10 <- run(VSEV, udca(500))
##  s11 <- run(VSEV, udca(750))
##
##  --- 12-17  the itch axis ------------------------------------------
##  s12 <- run(SEV,  rif(300))                      # rifampicin 300 bid
##  s13 <- run(SEV,  rif(150))                      # low dose
##  s14 <- run(SEV,  c(udca(500), rif(300)))        # combination
##  s15 <- run(VSEV, c(udca(500), rif(300)))
##  s16 <- run(SEV,  ntx(50))                       # naltrexone
##  s17 <- run(SEV,  anti(4))                       # antihistamine
##
##  --- 18-22  other and hypothetical agents --------------------------
##  s18 <- run(SEV,  chol(4))                       # cholestyramine 16 g/d
##  s19 <- run(SEV,  chol(4), extra = list(VKDOSE = 1.6))
##  s20 <- run(SEV,  same(1000))                    # SAMe
##  s21 <- run(SEV,  c(udca(500), same(1000)))
##  s22 <- run(VSEV, ibat(4))                       # IBAT inhibitor
##  s23 <- run(SEV,  ibat(4))
##
##  --- 24-28  genotype subgroups -------------------------------------
##  s24 <- run(list(GBSEP = 0.55, GMDR3 = 0.95, GSULT = 2.2))  # ABCB11 hom
##  s25 <- run(list(GBSEP = 0.92, GMDR3 = 0.30, GSULT = 2.2))  # ABCB4 severe
##  s26 <- run(list(GFIC1 = 0.68, GMDR3 = 0.80, GSULT = 2.4))  # ATP8B1
##  s27 <- run(c(SEV, list(GFXR = 0.55)))                      # NR1H4 hypo
##  s28 <- run(c(SEV, list(GP = 0.60)))                        # low placenta
##
##  --- 29-30  twins --------------------------------------------------
##  s29 <- run(c(MILD1, list(TWIN = 1.55)))
##  s30 <- run(c(SEV,   list(TWIN = 1.55)), udca(500))
##
##  --- 31-35  delivery timing ----------------------------------------
##  s31 <- run(VSEV, del_ga = 36)
##  s32 <- run(VSEV, beta(35), del_ga = 36)
##  s33 <- run(VSEV, del_ga = 38)
##  s34 <- run(SEV,  del_ga = 37)
##  s35 <- run(VSEV, c(udca(750), rif(300), ibat(4), beta(35)), del_ga = 36)
##
##  --- reading the output --------------------------------------------
##  The four columns that matter, and why they disagree with each other:
##    TBA        what the clinic measures (counts UDCA and hyocholate)
##    TBA_ENDO   the endogenous fraction (what the clinic thinks it measures)
##    FCLo       fetal cytotoxic load  -> drives ARRIo -> drives SBRISK
##    VAS        itch, which tracks ATX and not TBA
##  Compare s08 with s12: the first halves TBA and leaves VAS alone, the
##  second abolishes VAS and leaves TBA mid-table.  Compare SBRISK with
##  NICURISK across s31/s33: the two move in opposite directions with
##  gestational age, and their sum is what the delivery decision is
##  actually minimising.
## =====================================================================
