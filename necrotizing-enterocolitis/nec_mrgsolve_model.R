## =============================================================================
##  nec_mrgsolve_model.R
##  Necrotising Enterocolitis (NEC) of the preterm infant — QSP / PK-PD model
##  신생아 괴사성 장염 — 정량적 시스템 약리학 모델 (mrgsolve)
##
##  39 ODEs.  Mucosal barrier · luminal ecology · innate signalling · splanchnic
##  perfusion · the lesion · systemic response · distant organs · drug PK.
##
##  -----------------------------------------------------------------------------
##  THE ONE STRUCTURAL CLAIM THIS FILE IMPLEMENTS
##  -----------------------------------------------------------------------------
##  NEC is written as a CLOSED POSITIVE LOOP on a single mucosal state variable
##  E (enterocyte / villus integrity), whose loop gain is set by the luminal
##  pathobiont load B:
##
##      E -> BI = E·(TJ)·(MUC)                       barrier integrity
##        -> Pb = Pmin + (Pmax-Pmin)·(1-BI)^hP       permeability
##        -> Jtr = Pb · B                            translocation flux
##        -> TLR4s = TLR4expr · Jtr/(Jtr + Ktlr)     innate signal
##        -> (a) epithelial loss  kA·INJ·E/(E+KEap)  UP
##           (b) crypt proliferation  1/(1+(INJ/Ki)^2)  DOWN
##        -> E DOWN  ... and round again
##
##  Both injury arms act on E, and E is what sets Pb.  That is the loop.  The
##  1-D reduced field g(E) has an odd number of roots; a saddle-node pair
##  appears at B = B_lo (NEC becomes POSSIBLE) and the healthy root is
##  annihilated at B = B_hi (NEC becomes INEVITABLE).  Between the two, fate is
##  decided by a PRECIPITATING INSULT, not by the load.
##
##  Two further structural facts are built in on purpose:
##
##   * FEEDING IS TWO-SIGNED.  Enteral substrate appears with a + sign in the
##     trophic drive Ftroph, a + sign in dB/dt (it is bacterial substrate), and
##     a - sign in the splanchnic O2 balance (it raises demand).  Gut rest is
##     therefore not a safe default: Ftroph falls to its atrophy floor.
##
##   * THRESHOLD-MOVERS AND STATE-MOVERS ARE DIFFERENT DRUGS.  HMO and sIgA
##     enter ONLY through Ktlr; probiotics enter ONLY through B.  Human milk is
##     the only intervention here that is both, which is why nothing on the
##     formulary replaces it.
##
##  -----------------------------------------------------------------------------
##  IMPORTANT — HOW THESE EQUATIONS WERE VERIFIED
##  -----------------------------------------------------------------------------
##  The build environment for this repository has no R runtime, so this file
##  could not be executed here.  Every equation below is re-implemented term for
##  term in `nec_reference_model.py` (pure-Python RK4) and actually integrated;
##  the captured log is `nec_reference_output.txt` and every number quoted in
##  `README.md` comes from that run.  If you change an equation here, change it
##  there too, or the two will drift apart.
##
##  -----------------------------------------------------------------------------
##  CALIBRATION TARGETS (aggregate, not patient-level)
##  -----------------------------------------------------------------------------
##   * NEC incidence by GA band in a real-world feeding mix (~60 % human milk,
##     ~55 % exposed to early empirical antibiotics):
##         24-25 wk ~12-15 % · 26-27 wk ~8-11 % · 28-29 wk ~4-6 % · 30-32 wk ~1-3 %
##   * Onset: postnatal day ~10-25, later at lower gestational age
##   * Pneumatosis requires a breached mucosa (Pb > PbPNE), never intact bowel
##   * Bell II -> III in a fulminant case: hours, not days
##   * Post-NEC stricture ~10-25 % of survivors
##   * Neonatal PK: ampicillin t1/2 ~3 h, gentamicin ~8 h, metronidazole ~33 h,
##     indomethacin ~19 h, ibuprofen ~59 h; CL scaled by (PMA/40)^1.3
##
##  DISCLAIMER: teaching / hypothesis-generating model.  Not validated for any
##  clinical, prescribing, or regulatory use.
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

nec_code <- '
$PARAM @annotated
// ---- gestational / postnatal maturation ------------------------------------
GA       :  28.0 : Gestational age at birth (wk)
BW       :  1.05 : Birth weight (kg)
PMA50    :  30.5 : PMA at half-maximal mucosal maturation (wk)
nMAT     :   6.0 : Hill coefficient of maturation
phiTLR   :  0.55 : Fraction of TLR4 over-expression removed by maturation
TLR4max  :  1.00 : Immature-gut TLR4 signalling capacity

// ---- barrier: enterocyte / villus integrity --------------------------------
kE       :  0.45 : Crypt renewal rate constant (1/d)
kA       :  0.90 : Injury-driven enterocyte loss (1/d)
KEap     :  0.08 : E at which injury-driven loss is half-maximal
Ki       :  0.25 : Injury at which crypt proliferation is half-blocked
nKi      :  2.00 : Steepness of the proliferation block
trophfl  :  0.25 : Ftroph floor when completely NPO (mucosal atrophy)
Kfeed    :  40.0 : Feed volume at half-maximal trophic drive (mL/kg/d)
eEGF     :  0.10 : Extra trophic drive from milk EGF / HB-EGF
eSCFA    :  0.35 : Extra trophic drive from butyrate
KsSCFA   :   8.0 : SCFA at half-maximal trophic effect (mM)

// ---- barrier accessories ---------------------------------------------------
kTJ      :  1.20 : Tight-junction assembly (1/d)
dTJ      :  1.10 : TJ disassembly by TNF / IL-1b / NO (1/d)
kMUC     :  1.00 : Mucus production (1/d)
dMUC     :  0.90 : Mucus degradation by neutrophil elastase (1/d)
kIgAm    :  1.60 : sIgA delivered per unit milk bioactivity (1/d)
kIgAe    :  0.10 : Endogenous sIgA, maturation-limited (1/d)
dIgA     :  1.50 : sIgA loss (1/d)
kDEF     :  0.80 : Paneth-cell defensin output (1/d)
dDEF     :  1.00 : Defensin turnover (1/d)

// ---- permeability / translocation -----------------------------------------
Pmin     : 0.050 : Paracellular permeability of intact barrier
Pmax     : 1.200 : Permeability of denuded mucosa
hP       :   3.0 : Steepness of permeability vs barrier loss
fP       :  0.15 : Probiotic relative translocation weight

// ---- innate signalling -----------------------------------------------------
Ktlr0    :  60.0 : Translocation flux for half-maximal TLR4 signal
thrboost :  1.00 : THRESHOLD-MOVER lever (oral 2-FL / anti-TLR4 / rPAF-AH)
hHMO     :  0.15 : HMO decoy-receptor raising of Ktlr
KHMO     :  0.40 : HMO pool for half of that effect (g/kg)
hIgA     :  0.10 : sIgA immune-exclusion raising of Ktlr
ktlr_on  :  24.0 : NF-kB activation (1/d)
ktlr_off :  24.0 : NF-kB deactivation (1/d)
iIL10    :  0.90 : IL-10 damping of NF-kB
iDEX     :  0.70 : Glucocorticoid damping of NF-kB (Emax)
EC50DEX  :  12.0 : Dexamethasone EC50 (ug/L)

// ---- cytokines -------------------------------------------------------------
kIL1     :   5.0 : IL-1b production
dIL1     :   6.0 : IL-1b turnover (1/d)
kTNF     :   6.0 : TNF-a production
dTNF     :   8.0 : TNF-a turnover (1/d)
kIL8     :   5.0 : IL-8 production
dIL8     :   4.0 : IL-8 turnover (1/d)
aTNF8    :   0.6 : TNF drive on IL-8
aIL18    :   0.6 : IL-1b drive on IL-8
kIL10    :   3.0 : IL-10 production (maturation-limited)
dIL10    :   3.0 : IL-10 turnover (1/d)
kPAF     :   4.0 : PAF production
dPAF0    :   1.2 : PAF-acetylhydrolase floor (preterm deficiency)
dPAFm    :   6.0 : Maturation-dependent PAF-AH capacity
kNEU     :   3.0 : Neutrophil recruitment
dNEU     :   2.0 : Neutrophil clearance (1/d)
NEUmax   :   4.0 : Saturation of the mucosal neutrophil pool
kNOx     :   2.5 : iNOS-derived NO production
dNOx     :   5.0 : NO / peroxynitrite turnover (1/d)
wNO      :  0.35 : Weight of nitrosative injury in INJ
wISCH    :  0.90 : Weight of ischaemia in INJ
wNEC     :  1.80 : Weight of the lesion itself in INJ (DAMPs)
SDRcrit  :  0.80 : O2 supply/demand ratio below which injury starts

// ---- luminal ecology -------------------------------------------------------
muB      :   4.5 : Pathobiont max growth (1/d)
KsB      :  0.20 : Substrate for half-maximal B growth (g/kg)
muC      :   3.0 : Commensal anaerobe max growth (1/d)
KsC      :  0.30 : Substrate for half-maximal C growth (g/kg)
muCH     :  2.60 : Commensal growth on HMO (a niche B cannot enter)
muP      :   3.4 : Probiotic max growth (1/d)
KsP      :  0.30 : Substrate for half-maximal P growth (g/kg)
Ktot     :  26.0 : Total colonisable niche (1e9 CFU/g)
alphaX   :  0.55 : Cross-taxon competition (<1 permits coexistence)
binfant  :  1.00 : HMO-utilising strain carried (1) or not (0.05)
kwash    :  1.50 : Stool washout (1/d)
kIgAB    :  0.12 : sIgA-mediated clearance of B (1/d)
kDEFB    :  0.25 : Defensin-mediated clearance of B (1/d)
KpH      :  10.0 : SCFA (mM) halving Enterobacteriaceae growth
seedChm  : 0.030 : Bifidobacterium seeding from milk (1e9/g/d)
seedCenv : 0.004 : Environmental anaerobe seeding
seedBenv : 0.010 : Environmental Enterobacteriaceae seeding
dysbio   :  1.00 : Multiplier on environmental B seeding
ySCFA    :  32.0 : SCFA yield per unit anaerobe fermentation
dSCFA    :   2.2 : SCFA absorption / washout (1/d)
yGAS     :  0.55 : Gas yield per unit fermentation
dGAS     :  1.60 : Gas clearance (1/d)

// ---- substrate handling ----------------------------------------------------
gPerMl   : 0.070 : Fermentable macronutrient per mL feed (g/mL)
Smuc     :  2.20 : Host-derived substrate: mucin O-glycans (g/kg/d)
Vabs     :  14.0 : Absorptive Vmax at E = 1 (g/kg/d)
Kabs     :  0.35 : Michaelis constant of absorption (g/kg)
ktrans   :  2.00 : Substrate loss by transit (1/d)
qferm    : 0.075 : Fermentation per unit biomass (g/kg/d per 1e9)
Ksf      :  0.30 : Michaelis constant of fermentation (g/kg)
hmoPerMl : 0.0130 : HMO in mother-s own milk (g/mL)
kHMOferm : 0.055 : HMO consumption per unit anaerobe (1/d per 1e9)
dHMO     :  1.50 : HMO washout (1/d)

// ---- perfusion -------------------------------------------------------------
kPERF    :   6.0 : Perfusion relaxation (1/d)
aPAF     :  0.45 : PAF-driven splanchnic vasoconstriction
KaPAF    :  0.60 : PAF for half-maximal vasoconstriction
aCOX     :  0.35 : Perfusion lost when vasodilator prostanoid tone is removed
aNEC     :  0.40 : Perfusion lost inside the lesion
aDEM     :  0.40 : Extra postprandial O2 demand at 200 mL/kg/d
kPGE     :   4.0 : Prostanoid synthesis (1/d)
EC50IND  :  0.35 : Indomethacin EC50 for COX inhibition (mg/L)
EC50IBU  :  12.0 : Ibuprofen EC50 for COX inhibition (mg/L)

// ---- lesion ----------------------------------------------------------------
INJth    :  0.35 : Injury threshold for transmural necrosis
knec     :  4.00 : Necrosis propagation (1/d)
nNEC     :   1.5 : Steepness of necrosis propagation
krep     :  0.30 : Granulation / re-epithelialisation (1/d)
NECrev   :  0.10 : Only injury below this can re-epithelialise
kPNE     :  3.20 : Pneumatosis formation (1/d)
PbPNE    :  0.20 : Permeability floor below which no gas dissects
dPNE     :  0.55 : Pneumatosis resolution (1/d)

// ---- systemic --------------------------------------------------------------
kLPS     :  0.35 : Endotoxaemia from translocation
dLPS     :   4.0 : Endotoxin clearance (1/d)
wNECLPS  :   8.0 : Amplification of endotoxaemia by the lesion
PLT0     : 250.0 : Baseline platelet count (1e9/L)
kPLT     :  0.35 : Platelet production (1/d)
kPLTc    :  0.55 : Platelet consumption
kCRP     :  22.0 : CRP synthesis
dCRP     :  0.60 : CRP elimination (1/d)
kLAC     :   6.0 : Lactate from hypoperfusion
kLAC2    :   5.0 : Lactate from the lesion
dLAC     :   3.0 : Lactate clearance (1/d)
LAC0     :   1.0 : Baseline lactate (mmol/L)
kNIN     :  0.55 : Neuro-inflammation accrual
dNIN     :  0.30 : Neuro-inflammation resolution (1/d)
kBIL     :  0.70 : TPN cholestasis accrual
dBIL     :  0.22 : Bilirubin clearance (1/d)
kcal_ml  :  0.68 : kcal per mL of feed
kcal_tpn :  62.0 : kcal/kg/d from parenteral nutrition
kcal_mnt :  60.0 : Maintenance kcal/kg/d
kcal_g   :   5.0 : kcal per g of weight gain
kKIN     :  0.06 : Aminoglycoside tubular injury accrual
dKIN     :  0.20 : Tubular recovery (1/d)

// ---- drug PK ---------------------------------------------------------------
Vamp     :  0.50 : Ampicillin V (L/kg)
CLamp    :  2.60 : Ampicillin CL (L/kg/d)
Vgen     :  0.50 : Gentamicin central V (L/kg)
CLgen    :  1.05 : Gentamicin CL (L/kg/d)
Qgen     :  0.35 : Gentamicin intercompartmental Q (L/kg/d)
Vgep     :  0.30 : Gentamicin peripheral V (L/kg)
Vmtz     :  0.72 : Metronidazole V (L/kg)
CLmtz    :  0.36 : Metronidazole CL (L/kg/d)
Vind     :  0.36 : Indomethacin V (L/kg)
CLind    :  0.31 : Indomethacin CL (L/kg/d)
Vibu     :  0.32 : Ibuprofen V (L/kg)
CLibu    :  0.09 : Ibuprofen CL (L/kg/d)
Vdex     :  1.10 : Dexamethasone V (L/kg)
CLdex    :  3.60 : Dexamethasone CL (L/kg/d)
// LUMINAL AVAILABILITY.  The microbiome is reshaped by the concentration that
// reaches the colonic lumen, not by the plasma concentration.  Ampicillin is
// biliary/renally excreted and arrives in quantity; metronidazole distributes
// almost freely; an intravenous aminoglycoside barely gets there at all.  This
// asymmetry — not any assumption about resistance — is what makes empirical
// ampicillin + gentamicin destroy the anaerobes while comparatively sparing
// the Enterobacteriaceae, and it is what makes a long course change sign.
fl_amp   :  0.45 : Ampicillin luminal fraction of plasma concentration
fl_gen   :  0.05 : Gentamicin luminal fraction (intravenous aminoglycoside)
fl_mtz   :  0.80 : Metronidazole luminal fraction
EmxBamp  :  1.10 : Ampicillin Emax on B (partial: resistance)
EC5Bamp  :  40.0 : Ampicillin EC50 on B (mg/L)
EmxCamp  :  3.40 : Ampicillin Emax on C (anaerobes)
EC5Camp  :   6.0 : Ampicillin EC50 on C (mg/L)
EmxBgen  :  4.20 : Gentamicin Emax on B
EC5Bgen  :   2.0 : Gentamicin EC50 on B (mg/L)
EmxCgen  :  0.35 : Gentamicin Emax on C
EC5Cgen  :   4.0 : Gentamicin EC50 on C (mg/L)
EmxBmtz  :  0.10 : Metronidazole Emax on B
EC5Bmtz  :   8.0 : Metronidazole EC50 on B (mg/L)
EmxCmtz  :  4.60 : Metronidazole Emax on C (anaerobes)
EC5Cmtz  :   4.0 : Metronidazole EC50 on C (mg/L)
EmxPamp  :  2.20 : Ampicillin Emax on the probiotic strain
EmxPgen  :  0.30 : Gentamicin Emax on the probiotic strain
EmxPmtz  :  1.30 : Metronidazole Emax on the probiotic strain

// ---- nutrition / regimen ---------------------------------------------------
fMOM     :  1.00 : Fraction of enteral volume that is mother-s own milk
fDM      :  0.00 : Fraction that is pasteurised donor milk
fFORM    :  0.00 : Fraction that is preterm formula
fortify  :  0.00 : Bovine fortifier fraction of volume
tfeed0   :  2.00 : Postnatal day enteral feeding starts
feedini  :  20.0 : Starting feed volume (mL/kg/d)
feedrate :  20.0 : Advance rate (mL/kg/d per day)
feedmax  : 160.0 : Target feed volume (mL/kg/d)
NPO_ON   :  0.00 : Set to 1 to hold gut rest (overrides the feed schedule)

// ---- clinical response to a diagnosis of NEC --------------------------------
// The response is NOT instantaneous: suspicion -> abdominal film -> nil by
// mouth actually in force.  The lesion keeps growing during that lag, and how
// deep the collapse gets during the lag is exactly what separates medical from
// surgical NEC in this model.
dx_lag     : 1.20 : Suspicion -> gut rest actually started (d)
npo_on_nec : 1.00 : Set to 1 for gut rest once Bell >= II
nec_abx    : 1.00 : Set to 1 for therapeutic triple antibiotics once Bell >= II
rx_days    : 10.0 : Duration of the therapeutic course (d)

// ---- precipitating insults (up to three boxcars) ---------------------------
// Inside the bistable window an infant sits on the healthy branch and stays
// there unless something transiently pushes E below E*.  Apnoea-bradycardia,
// hypotension, PDA diastolic steal, a packed-cell transfusion, late-onset
// sepsis.  DEP = splanchnic perfusion depth, NFK = extra NF-kB drive.
INS1A    : 1e6   : Insult 1 start (d)
INS1B    : 1e6   : Insult 1 end (d)
INS1DEP  : 0.00  : Insult 1 perfusion depth
INS1NFK  : 0.00  : Insult 1 extra NF-kB drive
INS2A    : 1e6   : Insult 2 start (d)
INS2B    : 1e6   : Insult 2 end (d)
INS2DEP  : 0.00  : Insult 2 perfusion depth
INS2NFK  : 0.00  : Insult 2 extra NF-kB drive
INS3A    : 1e6   : Insult 3 start (d)
INS3B    : 1e6   : Insult 3 end (d)
INS3DEP  : 0.00  : Insult 3 perfusion depth
INS3NFK  : 0.00  : Insult 3 extra NF-kB drive

$CMT @annotated
// ---- mucosal barrier -------------------------------------------------------
E     : Enterocyte / villus integrity (0-1 index)
TJ    : Tight-junction (ZO-1 / claudin-3) index
MUC   : Mucus layer adequacy (0-1)
IgA   : Luminal secretory IgA (0-1 index)
DEF   : Paneth-cell defensin index (0-1)
// ---- luminal ecology -------------------------------------------------------
B     : Pathobiont (Enterobacteriaceae) load (1e9 CFU/g)
C     : Commensal obligate anaerobe load (1e9 CFU/g)
P     : Administered probiotic load (1e9 CFU/g)
SUB   : Unabsorbed fermentable substrate (g/kg)
SCFA  : Luminal butyrate-equivalent SCFA (mM)
GAS   : Luminal / intramural gas index
HMO   : Human-milk oligosaccharide pool (g/kg)
// ---- innate signalling -----------------------------------------------------
TLRS  : TLR4-MyD88-NF-kB activity (0-1 index)
IL1B  : Mucosal IL-1beta (index)
TNFA  : Mucosal TNF-alpha (index)
IL8   : Mucosal IL-8 / CXCL8 (index)
IL10  : Mucosal IL-10 (index)
PAF   : Platelet-activating factor (index)
NEU   : Mucosal neutrophil infiltrate (index)
NOX   : iNOS-derived NO / peroxynitrite (index)
// ---- perfusion and lesion --------------------------------------------------
PERF  : Mesenteric perfusion (1 = normal)
NECA  : Transmural necrotic bowel area fraction (0-1)
PNEU  : Pneumatosis intestinalis radiographic index (0-1)
// ---- systemic --------------------------------------------------------------
LPSP  : Plasma endotoxin / bacteraemia index
PLT   : Platelet count (1e9/L)
CRP   : C-reactive protein (mg/L)
LAC   : Arterial lactate (mmol/L)
NIN   : Neuro-inflammation (NDI driver) index
BILI  : Conjugated bilirubin (mg/dL)
WT    : Body weight (kg)
KIN   : Aminoglycoside tubular injury index
// ---- drug PK ---------------------------------------------------------------
AMPC  : Ampicillin plasma (mg/L)
GENC  : Gentamicin plasma central (mg/L)
GENP  : Gentamicin peripheral / renal cortex (mg/L-eq)
MTZC  : Metronidazole plasma (mg/L)
INDC  : Indomethacin plasma (mg/L)
IBUC  : Ibuprofen plasma (mg/L)
DEXC  : Dexamethasone plasma (ug/L)
PGE   : Mucosal vasodilator prostanoid index (0-1)
// ---- latch: time elapsed since the FIRST radiographic diagnosis ------------
DXT   : Days since Bell II was first reached (latching clock)

$MAIN
// initial conditions are gestational-age dependent
double MAT0 = pow(GA/PMA50, nMAT) / (1.0 + pow(GA/PMA50, nMAT));
E_0    = fmin(0.985, 0.70 + 0.026*(GA - 24.0));
TJ_0   = 0.30 + 0.65*MAT0;
MUC_0  = 0.30 + 0.60*MAT0;
IgA_0  = 0.02;
DEF_0  = 0.15 + 0.55*MAT0;
B_0    = 0.010;
C_0    = 0.005;
SUB_0  = 0.02;
SCFA_0 = 0.50;
GAS_0  = 0.05;
PERF_0 = 1.00;
PLT_0  = PLT0;
CRP_0  = 1.00;
LAC_0  = LAC0;
BILI_0 = 0.40;
WT_0   = BW;
PGE_0  = 1.00;

$ODE
// ===========================================================================
//  0.  maturation, milk bioactivity, regimen
// ===========================================================================
double PMA = GA + SOLVERTIME/7.0;
double MAT = pow(PMA/PMA50, nMAT) / (1.0 + pow(PMA/PMA50, nMAT));

// pasteurised donor milk keeps HMO but loses sIgA / lactoferrin / EGF
double bio = fMOM + 0.45*fDM;

// ---- clinical response: a latching clock started by the first diagnosis ----
// DXT counts days since Bell II was first reached and never resets, so the
// response fires ONCE rather than being re-triggered at every step.
double bell2now = (PNEU > 0.30 || NECA > 0.10) ? 1.0 : 0.0;
dxdt_DXT = (bell2now > 0.5 || DXT > 0.0) ? 1.0 : 0.0;
double onRx  = (DXT > dx_lag && DXT < dx_lag + rx_days) ? 1.0 : 0.0;
double restD = onRx*npo_on_nec;                    // gut rest in force

// feed schedule.  After a course of gut rest, feeds restart from feedini.
double onNPO = (NPO_ON > 0.5 || restD > 0.5) ? 1.0 : 0.0;
double trestart = (DXT > dx_lag + rx_days) ? (SOLVERTIME - (DXT - dx_lag - rx_days))
                                           : tfeed0;
double tref = fmax(tfeed0, trestart);
double FEED  = 0.0;
if (SOLVERTIME >= tfeed0 && onNPO < 0.5) {
  FEED = feedini + feedrate*(SOLVERTIME - tref);
  FEED = fmin(feedmax, fmax(0.0, FEED));
}
double enteral = FEED / (feedmax > 1e-9 ? feedmax : 1e-9);
if (enteral > 1.0) enteral = 1.0;

// precipitating insults
double insP = 0.0, insN = 0.0;
if (SOLVERTIME >= INS1A && SOLVERTIME < INS1B) { insP += INS1DEP; insN += INS1NFK; }
if (SOLVERTIME >= INS2A && SOLVERTIME < INS2B) { insP += INS2DEP; insN += INS2NFK; }
if (SOLVERTIME >= INS3A && SOLVERTIME < INS3B) { insP += INS3DEP; insN += INS3NFK; }

// ===========================================================================
//  1.  barrier integrity -> permeability -> translocation   (THE LOOP)
// ===========================================================================
double BI  = E * (0.35 + 0.65*TJ) * (0.55 + 0.45*MUC);
BI = fmin(1.0, fmax(0.0, BI));
double Pb  = Pmin + (Pmax - Pmin) * pow(1.0 - BI, hP);
double Jtr = Pb * (B + fP*P);

// ===========================================================================
//  2.  innate signalling threshold and NF-kB
// ===========================================================================
double Ktlr = Ktlr0 * thrboost
              * (1.0 + hHMO * HMO/(HMO + KHMO))
              * (1.0 + hIgA * IgA);
double TLR4expr = TLR4max * (1.0 - phiTLR*MAT);
double dexE = iDEX * DEXC/(DEXC + EC50DEX);
dxdt_TLRS = ktlr_on * TLR4expr * (1.0 - dexE) * (Jtr/(Jtr + Ktlr) + insN)
            - ktlr_off * TLRS * (1.0 + iIL10*IL10);

// ===========================================================================
//  3.  splanchnic O2 balance.  Feeding raises DEMAND, so the same perfusion
//      buys less reserve on full feeds than on gut rest.
// ===========================================================================
double demand = 1.0 + aDEM * FEED/200.0;
double ISCH   = fmax(0.0, 1.0 - (PERF/demand)/SDRcrit);

// total injury.  The lesion is itself an inflammatory stimulus (DAMPs), which
// is what makes it irreversible above NECa_crit = INJth / wNEC.
double INJ = TLRS + wNO*NOX + wISCH*ISCH + wNEC*NECA;

// ===========================================================================
//  4.  enterocyte balance
// ===========================================================================
double Ftroph = (trophfl + (1.0 - trophfl)*FEED/(FEED + Kfeed))
                * (1.0 + eEGF*bio)
                * (1.0 + eSCFA*SCFA/(SCFA + KsSCFA));
double block  = 1.0 / (1.0 + pow(INJ/Ki, nKi));
dxdt_E = kE*(1.0 - E)*Ftroph*block - kA*INJ*E/(E + KEap);

dxdt_TJ  = kTJ*MAT*(1.0 - TJ)*(1.0 + 0.5*SCFA/(SCFA + KsSCFA))
           - dTJ*TJ*(0.45*TNFA + 0.45*IL1B + 0.35*NOX);
dxdt_MUC = kMUC*MAT*(1.0 - MUC)*(1.0 + 0.4*bio)
           - dMUC*MUC*(0.6*NEU + 0.2*TLRS);
dxdt_IgA = kIgAm*bio*enteral + kIgAe*MAT - dIgA*IgA - 0.05*IgA*B;
dxdt_DEF = kDEF*MAT*E*(1.0 - DEF) - dDEF*DEF;

// ===========================================================================
//  5.  antibiotic kill rates (Emax per taxon).  Note the asymmetry that drives
//      the sign change: ampicillin and metronidazole hit C much harder than B.
// ===========================================================================
double aL = fl_amp*AMPC;        // luminal, not plasma, concentrations
double gL = fl_gen*GENC;
double mL = fl_mtz*MTZC;
double killB = EmxBamp*aL/(aL + EC5Bamp)
             + EmxBgen*gL/(gL + EC5Bgen)
             + EmxBmtz*mL/(mL + EC5Bmtz);
double killC = EmxCamp*aL/(aL + EC5Camp)
             + EmxCgen*gL/(gL + EC5Cgen)
             + EmxCmtz*mL/(mL + EC5Cmtz);
double killP = EmxPamp*aL/(aL + EC5Camp)
             + EmxPgen*gL/(gL + EC5Bgen)
             + EmxPmtz*mL/(mL + EC5Cmtz);

// ===========================================================================
//  6.  luminal ecology.  alphaX < 1 means the two guilds are not competing for
//      one interchangeable niche, so coexistence is possible.
// ===========================================================================
double occB = fmin(1.0, (B + alphaX*(C + P))/Ktot);
double occC = fmin(1.0, (C + P + alphaX*B)/Ktot);
double fpH  = 1.0 / (1.0 + pow(SCFA/KpH, 2.0));
double gB   = muB*SUB/(SUB + KsB)*(1.0 - occB)*fpH;
double gC   = (muC*SUB/(SUB + KsC) + muCH*binfant*HMO/(HMO + KHMO))*(1.0 - occC);
double gP   = muP*SUB/(SUB + KsP)*(1.0 - occC);

dxdt_B = B*(gB - kwash - killB - kIgAB*IgA - kDEFB*DEF) + seedBenv*dysbio;
dxdt_C = C*(gC - kwash - killC) + seedChm*bio*enteral + seedCenv;
dxdt_P = P*(gP - kwash - killP);

// ===========================================================================
//  7.  substrate, SCFA, gas, HMO.  Note the SECOND positive feedback:
//      E down -> absorption down -> SUB up -> B up -> E down.
// ===========================================================================
double Sin    = FEED*gPerMl*(1.0 + 0.25*fortify) + Smuc*(0.5 + 0.5*E);
double absorb = Vabs*E*SUB/(SUB + Kabs);
double ferm   = qferm*(B + C + P)*SUB/(SUB + Ksf);
dxdt_SUB  = Sin - absorb - ferm - ktrans*SUB;
dxdt_SCFA = ySCFA*qferm*(C + P)*SUB/(SUB + Ksf) - dSCFA*SCFA;
dxdt_GAS  = yGAS*ferm*(1.0 + 1.5*(B/(B + 5.0))) - dGAS*GAS;
dxdt_HMO  = FEED*hmoPerMl*fMOM + FEED*hmoPerMl*0.95*fDM
            - kHMOferm*(C + P)*HMO/(HMO + KHMO) - dHMO*HMO;

// ===========================================================================
//  8.  cytokines and lipid mediators
// ===========================================================================
double dexf = 1.0 - dexE;
double NEUs = NEU/(1.0 + NEU);      // saturating feed-forward
double PAFs = PAF/(1.0 + PAF);
dxdt_IL1B = kIL1*dexf*(TLRS + 0.35*NEUs) - dIL1*IL1B;
dxdt_TNFA = kTNF*dexf*(TLRS + 0.20*NEUs) - dTNF*TNFA;
dxdt_IL8  = kIL8*dexf*(TLRS + aTNF8*TNFA + aIL18*IL1B) - dIL8*IL8;
dxdt_IL10 = kIL10*MAT*(TNFA + IL1B) - dIL10*IL10;
double PAFAH = dPAF0 + dPAFm*MAT;   // preterm PAF-AH deficiency
dxdt_PAF  = kPAF*(TLRS + 0.5*TNFA) - PAFAH*PAF;
dxdt_NEU  = kNEU*IL8*(1.0 + 0.8*PAFs)*(1.0 - NEU/NEUmax) - dNEU*NEU;
dxdt_NOX  = kNOx*(TNFA + IL1B) - dNOx*NOX;

// ===========================================================================
//  9.  perfusion
// ===========================================================================
double coxI = fmin(1.0, INDC/(INDC + EC50IND) + IBUC/(IBUC + EC50IBU));
dxdt_PGE = kPGE*((1.0 - coxI) - PGE);
double target = 1.0 - aPAF*PAF/(PAF + KaPAF) - aCOX*(1.0 - PGE)
                - aNEC*NECA - insP;
target = fmin(1.0, fmax(0.15, target));
dxdt_PERF = kPERF*(target - PERF);

// ===========================================================================
// 10.  the lesion.  Repair applies ONLY below NECrev: transmural necrosis
//      cannot re-epithelialise.
// ===========================================================================
double over = fmax(0.0, INJ - INJth);
dxdt_NECA = knec*(1.0 - NECA)*pow(over, nNEC)
            - krep*fmin(NECA, NECrev)*E*(INJ < INJth ? 1.0 : 0.0);
// gas dissects into the wall only through an already-breached mucosa
dxdt_PNEU = kPNE*GAS*fmax(0.0, Pb - PbPNE)*(1.0 - PNEU) - dPNE*PNEU;

// ===========================================================================
// 11.  systemic, liver, brain, kidney, growth
// ===========================================================================
dxdt_LPSP = kLPS*Jtr*(1.0 + wNECLPS*NECA) - dLPS*LPSP;
dxdt_PLT  = kPLT*(PLT0 - PLT) - kPLTc*PLT*(0.5*LPSP + 3.0*NECA);
dxdt_CRP  = kCRP*(0.6*IL1B + 0.4*TNFA) - dCRP*CRP;
dxdt_LAC  = kLAC*pow(fmax(0.0, 1.0 - PERF), 2.0) + kLAC2*NECA - dLAC*(LAC - LAC0);
dxdt_NIN  = kNIN*(0.4*TNFA + 0.4*IL1B + 0.5*LPSP) - dNIN*NIN;

double tpn = fmax(0.0, 1.0 - enteral);
dxdt_BILI = kBIL*tpn*(1.0 + 0.5*LPSP) - dBIL*BILI*(1.0 + enteral);
double kcal_in = (absorb/gPerMl)*kcal_ml + kcal_tpn*tpn;
dxdt_WT   = fmax(-0.004, (kcal_in - kcal_mnt*(1.0 + 0.30*NECA))/kcal_g/1000.0);
dxdt_KIN  = kKIN*GENP - dKIN*KIN;

// ===========================================================================
// 12.  drug PK, clearance scaled by post-menstrual age
// ===========================================================================
double sc = pow(PMA/40.0, 1.3);
// therapeutic amp + gent + metronidazole while Bell >= II, entered as the
// continuous-rate equivalent of q12h / q36h / q24h dosing
double rx = onRx*nec_abx;
dxdt_AMPC = -CLamp*sc/Vamp*AMPC + rx*(100.0/Vamp)/0.5;
dxdt_GENC = -CLgen*sc/Vgen*GENC - Qgen/Vgen*(GENC - GENP) + rx*(4.0/Vgen)/1.5;
dxdt_GENP =  Qgen/Vgep*(GENC - GENP);
dxdt_MTZC = -CLmtz*sc/Vmtz*MTZC + rx*(7.5/Vmtz)/1.0;
dxdt_INDC = -CLind*sc/Vind*INDC;
dxdt_IBUC = -CLibu*sc/Vibu*IBUC;
dxdt_DEXC = -CLdex*sc/Vdex*DEXC;

$TABLE
double MATo   = pow((GA + TIME/7.0)/PMA50, nMAT) /
                (1.0 + pow((GA + TIME/7.0)/PMA50, nMAT));
double BIo    = E*(0.35 + 0.65*TJ)*(0.55 + 0.45*MUC);
double Pbo    = Pmin + (Pmax - Pmin)*pow(1.0 - fmin(1.0, fmax(0.0, BIo)), hP);
double Jtro   = Pbo*(B + fP*P);
double Ktlro  = Ktlr0*thrboost*(1.0 + hHMO*HMO/(HMO + KHMO))*(1.0 + hIgA*IgA);
double onRxo  = (DXT > dx_lag && DXT < dx_lag + rx_days) ? 1.0 : 0.0;
double FEEDo  = 0.0;
if (TIME >= tfeed0 && NPO_ON < 0.5 && onRxo*npo_on_nec < 0.5) {
  double tro = (DXT > dx_lag + rx_days) ? (TIME - (DXT - dx_lag - rx_days)) : tfeed0;
  FEEDo = fmin(feedmax, fmax(0.0, feedini + feedrate*(TIME - fmax(tfeed0, tro))));
}
double demo   = 1.0 + aDEM*FEEDo/200.0;
double ISCHo  = fmax(0.0, 1.0 - (PERF/demo)/SDRcrit);
double INJo   = TLRS + wNO*NOX + wISCH*ISCHo + wNEC*NECA;
double NECcrt = INJth/wNEC;                       // the point of no return
// modified Bell stage read out of the state vector
double BELL = 0;
if (GAS > 0.55 || CRP > 12.0 || NECA > 0.01)          BELL = 1;
if (PNEU > 0.30 || NECA > 0.10)                        BELL = 2;
if (NECA > 0.28 || (PLT < 100.0 && LAC > 4.0))         BELL = 3;

capture Barrier      = BIo;
capture Permeability = Pbo;
capture Translocat   = Jtro;
capture Threshold    = Ktlro;
capture Injury       = INJo;
capture Ischaemia    = ISCHo;
capture FeedVol      = FEEDo;
capture Maturation   = MATo;
capture NECa_crit    = NECcrt;
capture BellStage    = BELL;
capture OnTreatment  = onRxo;
capture WtGain_g     = (WT - BW)*1000.0;
'

mod <- mcode("nec", nec_code, end = 56, delta = 0.05)

## =============================================================================
##  DOSING AND REGIMEN HELPERS
## =============================================================================

## Empirical / therapeutic antibiotics, PDA treatment, probiotics.
## Doses are entered as bolus additions to the plasma compartment expressed as
## amount/V, i.e. directly in mg/L, because the compartments are concentrations.
nec_events <- function(abx_start = 0, abx_days = 0, mtz_days = 0,
                       ind_start = 3, ind_days = 0,
                       ibu_start = 3, ibu_days = 0,
                       dex_days = 0,
                       prob_dose = 0, prob_start = 2, prob_stop = 56,
                       Vamp = 0.50, Vgen = 0.50, Vmtz = 0.72,
                       Vind = 0.36, Vibu = 0.32, Vdex = 1.10) {
  ev <- NULL
  add <- function(ev, times, cmt, amt) {
    if (!length(times)) return(ev)
    dplyr::bind_rows(ev, data.frame(ID = 1, time = times, cmt = cmt,
                                    amt = amt, evid = 1))
  }
  if (abx_days > 0) {
    ev <- add(ev, seq(abx_start, abx_start + abx_days - 1e-9, by = 0.5),
              "AMPC", 100 / Vamp)                       # 100 mg/kg q12h
    ev <- add(ev, seq(abx_start, abx_start + abx_days - 1e-9, by = 1.5),
              "GENC", 4 / Vgen)                         # 4 mg/kg q36h
  }
  if (mtz_days > 0) {
    ev <- add(ev, seq(abx_start, abx_start + mtz_days - 1e-9, by = 1.0),
              "MTZC", 7.5 / Vmtz)                       # 7.5 mg/kg q24h
  }
  if (ind_days > 0) {
    ev <- add(ev, seq(ind_start, ind_start + ind_days - 1e-9, by = 1.0),
              "INDC", 0.2 / Vind)                       # 0.2 mg/kg q24h x3
  }
  if (ibu_days > 0) {
    tt <- seq(ibu_start, ibu_start + ibu_days - 1e-9, by = 1.0)
    ev <- add(ev, tt, "IBUC", c(10, rep(5, max(0, length(tt) - 1))) / Vibu)
  }
  if (dex_days > 0) {
    ev <- add(ev, seq(0, dex_days - 1e-9, by = 0.5), "DEXC", 150 / Vdex)
  }
  if (prob_dose > 0) {
    ev <- add(ev, seq(prob_start, prob_stop - 1e-9, by = 1.0), "P", prob_dose)
  }
  if (is.null(ev)) ev <- data.frame(ID = 1, time = 0, cmt = "AMPC",
                                   amt = 0, evid = 1)
  dplyr::arrange(ev, time)
}

## The high-risk index-patient insult schedule used by the scenarios below:
## two transfusions and two hypoperfusion episodes, the second of them deep and
## landing on day 14, when the luminal load has finished building.
std_insults <- list(INS1A = 6.0,  INS1B = 6.3,  INS1DEP = 0.22, INS1NFK = 0.000,
                    INS2A = 14.0, INS2B = 14.6, INS2DEP = 0.40, INS2NFK = 0.010,
                    INS3A = 20.0, INS3B = 20.5, INS3DEP = 0.16, INS3NFK = 0.030)

milk_mix <- function(kind) {
  switch(kind,
    MOM     = list(fMOM = 1.0, fDM = 0.0, fFORM = 0.0),
    DONOR   = list(fMOM = 0.0, fDM = 1.0, fFORM = 0.0),
    FORMULA = list(fMOM = 0.0, fDM = 0.0, fFORM = 1.0),
    MIXED   = list(fMOM = 0.5, fDM = 0.2, fFORM = 0.3),
    stop("unknown milk kind: ", kind))
}

run_nec <- function(GA = 28, BW = 1.05, milk = "MOM", feedrate = 20,
                    abx_days = 0, mtz_days = 0, ind_days = 0, ibu_days = 0,
                    prob_dose = 0, dysbio = 1.5, binfant = 1.0,
                    thrboost = 1.0, insults = std_insults,
                    NPO_ON = 0, end = 56, ...) {
  pars <- c(list(GA = GA, BW = BW, feedrate = feedrate, dysbio = dysbio,
                 binfant = binfant, thrboost = thrboost, NPO_ON = NPO_ON),
            milk_mix(milk), insults, list(...))
  ev <- nec_events(abx_days = abx_days, mtz_days = mtz_days,
                   ind_days = ind_days, ibu_days = ibu_days,
                   prob_dose = prob_dose)
  mod %>%
    param(pars) %>%
    mrgsim(data = ev, end = end, delta = 0.05) %>%
    as_tibble()
}

## =============================================================================
##  SCENARIOS  (11)
##  These describe a HIGH-RISK INDEX PATIENT, not a population average: a
##  dysbiotic infant who receives one deep hypotensive episode on day 14.  A
##  median infant on the same regimens does not tip at all — comparing arms in
##  an infant who never approaches the separatrix would show nothing.
## =============================================================================
scenarios <- list(
  S1  = list(label = "28주 · 모유 · 항생제 없음",
             args  = list(GA = 28, BW = 1.05, milk = "MOM")),
  S2  = list(label = "28주 · 분유 · 항생제 없음",
             args  = list(GA = 28, BW = 1.05, milk = "FORMULA")),
  S3  = list(label = "28주 · 분유 · 경험적 항생제 7일",
             args  = list(GA = 28, BW = 1.05, milk = "FORMULA", abx_days = 7)),
  S4  = list(label = "28주 · 분유 · 프로바이오틱스",
             args  = list(GA = 28, BW = 1.05, milk = "FORMULA", prob_dose = 1.10)),
  S5  = list(label = "25주 · 분유 · 항생제 7일 + 인도메타신",
             args  = list(GA = 25, BW = 0.68, milk = "FORMULA", abx_days = 7,
                          ind_days = 3, dysbio = 1.8)),
  S6  = list(label = "25주 · 모유 + 프로바이오틱스 · 항생제 2일",
             args  = list(GA = 25, BW = 0.68, milk = "MOM", prob_dose = 1.10,
                          abx_days = 2, dysbio = 1.8)),
  S7  = list(label = "32주 · 분유",
             args  = list(GA = 32, BW = 1.75, milk = "FORMULA", feedrate = 30)),
  S8  = list(label = "28주 · 기증모유(저온살균)",
             args  = list(GA = 28, BW = 1.05, milk = "DONOR")),
  S9  = list(label = "28주 · 완전 금식 + TPN",
             args  = list(GA = 28, BW = 1.05, milk = "MOM", NPO_ON = 1)),
  S10 = list(label = "28주 · 분유 · 항생제 10일 + 메트로니다졸",
             args  = list(GA = 28, BW = 1.05, milk = "FORMULA", abx_days = 10,
                          mtz_days = 10)),
  S11 = list(label = "27주 · 분유 · 중증 (구조 개입 없음)",
             args  = list(GA = 27, BW = 0.92, milk = "FORMULA", feedrate = 30,
                          abx_days = 7, mtz_days = 5, ind_days = 3,
                          dysbio = 1.8))
)

run_all_scenarios <- function() {
  dplyr::bind_rows(lapply(names(scenarios), function(k) {
    s <- scenarios[[k]]
    do.call(run_nec, s$args) %>% mutate(scenario = k, label = s$label)
  }))
}

## Endpoint summary: peak Bell stage, onset day, whether the lesion passed the
## point of no return, and the extra-intestinal readouts.
summarise_scenarios <- function(sim) {
  sim %>%
    group_by(scenario, label) %>%
    summarise(
      peak_bell   = max(BellStage),
      onset_day   = ifelse(any(BellStage >= 2), min(time[BellStage >= 2]), NA),
      bell3_day   = ifelse(any(BellStage >= 3), min(time[BellStage >= 3]), NA),
      passed_crit = any(NECA > first(NECa_crit)),
      NECa_max    = max(NECA),
      E_min       = min(E),
      B_max       = max(B),
      SCFA_end    = last(SCFA),
      bili_end    = last(BILI),
      NDI_index   = last(NIN),
      tubular     = last(KIN),
      wt_gain_g   = last(WtGain_g),
      .groups     = "drop") %>%
    mutate(bell2_to_3_h = (bell3_day - onset_day)*24)
}

## =============================================================================
##  THE TWO CRITICAL LOADS — computed, not asserted
##  Reproduces section 1 of nec_reference_output.txt.  Holds B, HMO, IgA and the
##  fast block quasi-static and counts the roots of g(E) as B is raised.
## =============================================================================
reduced_field <- function(E, B, pp, ctx) {
  Pb  <- pp$Pmin + (pp$Pmax - pp$Pmin) * (1 - E)^pp$hP
  INJ <- ctx$injury_of(Pb * B)
  blk <- 1/(1 + (INJ/pp$Ki)^pp$nKi)
  pp$kE*(1 - E)*ctx$Ftroph*blk - pp$kA*INJ*E/(E + pp$KEap)
}

## Steady state of the fast block at a held translocation flux — the same nine
## equations as the ODE above, solved by damped fixed-point iteration.
quasistatic_injury <- function(Jtr, pp, ctx, niter = 250) {
  Ktlr  <- pp$Ktlr0*pp$thrboost*(1 + pp$hHMO*ctx$HMO/(ctx$HMO + pp$KHMO))*
           (1 + pp$hIgA*ctx$IgA)
  drive <- ctx$TLR4expr * Jtr/(Jtr + Ktlr)
  PAFAH <- pp$dPAF0 + pp$dPAFm*ctx$MAT
  x <- drive; TNF <- IL1 <- IL8 <- NEU <- PAF <- NOx <- IL10 <- 0; PERF <- 1
  w <- 0.25
  for (i in seq_len(niter)) {
    xn   <- drive/(1 + pp$iIL10*IL10)
    NEUs <- NEU/(1 + NEU); PAFs <- PAF/(1 + PAF)
    TNFn <- pp$kTNF*(xn + 0.20*NEUs)/pp$dTNF
    IL1n <- pp$kIL1*(xn + 0.35*NEUs)/pp$dIL1
    IL8n <- pp$kIL8*(xn + pp$aTNF8*TNF + pp$aIL18*IL1)/pp$dIL8
    I10n <- pp$kIL10*ctx$MAT*(TNF + IL1)/pp$dIL10
    PAFn <- pp$kPAF*(xn + 0.5*TNF)/PAFAH
    a    <- pp$kNEU*IL8*(1 + 0.8*PAFs)
    NEUn <- a/(pp$dNEU + a/pp$NEUmax)
    NOxn <- pp$kNOx*(TNF + IL1)/pp$dNOx
    tgt  <- 1 - pp$aPAF*PAF/(PAF + pp$KaPAF)
    PERFn <- min(1, max(0.15, tgt))
    x <- x + w*(xn - x);       TNF <- TNF + w*(TNFn - TNF)
    IL1 <- IL1 + w*(IL1n - IL1); IL8 <- IL8 + w*(IL8n - IL8)
    IL10 <- IL10 + w*(I10n - IL10); PAF <- PAF + w*(PAFn - PAF)
    NEU <- NEU + w*(NEUn - NEU); NOx <- NOx + w*(NOxn - NOx)
    PERF <- PERF + w*(PERFn - PERF)
  }
  ISCH <- max(0, 1 - (PERF/ctx$demand)/pp$SDRcrit)
  x + pp$wNO*NOx + pp$wISCH*ISCH
}

nec_context <- function(pp, milk = "FORMULA", GA = 28, feed = 150,
                        probiotic = FALSE) {
  mm  <- milk_mix(milk)
  bio <- mm$fMOM + 0.45*mm$fDM
  MAT <- ((GA + 2)/pp$PMA50)^pp$nMAT / (1 + ((GA + 2)/pp$PMA50)^pp$nMAT)
  HMO <- feed*pp$hmoPerMl*(mm$fMOM + 0.95*mm$fDM)/(pp$dHMO + 0.25)
  IgA <- min(1, (pp$kIgAm*bio*feed/pp$feedmax + pp$kIgAe*MAT)/pp$dIgA)
  SCFA <- if (bio > 0.5 || probiotic) 12.0 else 4.5
  Ftroph <- (pp$trophfl + (1 - pp$trophfl)*feed/(feed + pp$Kfeed)) *
            (1 + pp$eEGF*bio) * (1 + pp$eSCFA*SCFA/(SCFA + pp$KsSCFA))
  ctx <- list(HMO = HMO, IgA = IgA, SCFA = SCFA, Ftroph = Ftroph, MAT = MAT,
              TLR4expr = pp$TLR4max*(1 - pp$phiTLR*MAT),
              demand = 1 + pp$aDEM*feed/200)
  ## tabulate INJ(Jtr) once and interpolate — INJ depends on E and B only
  ## through Jtr, so this is exact, not an approximation
  js <- c(0, 400*((1:500)/500)^2.5)
  vs <- vapply(js, quasistatic_injury, numeric(1), pp = pp, ctx = ctx)
  ctx$injury_of <- stats::approxfun(js, vs, rule = 2)
  ctx
}

count_roots <- function(B, pp, ctx, ngrid = 1200) {
  xs <- seq(1/ngrid, 1 - 1/ngrid, length.out = ngrid)
  vs <- vapply(xs, reduced_field, numeric(1), B = B, pp = pp, ctx = ctx)
  idx <- which(vs[-length(vs)]*vs[-1] < 0)
  vapply(idx, function(i)
    stats::uniroot(reduced_field, c(xs[i], xs[i + 1]), B = B, pp = pp,
                   ctx = ctx, tol = 1e-10)$root, numeric(1))
}

## B_lo = the load at which bistability appears  (NEC becomes POSSIBLE)
## B_hi = the load at which the healthy root dies (NEC becomes INEVITABLE)
bifurcation_loads <- function(pp, ctx, Bmax = 240, nB = 240) {
  grid <- Bmax*(1:nB)/nB
  ks   <- vapply(grid, function(B) length(count_roots(B, pp, ctx)), numeric(1))
  hi_r <- vapply(grid, function(B) { r <- count_roots(B, pp, ctx)
                                    if (length(r)) max(r) else NA }, numeric(1))
  i_lo <- which(ks >= 3)[1]
  i_hi <- which(ks == 1 & !is.na(hi_r) & hi_r < 0.5)[1]
  refine <- function(a, b, pred) {
    for (i in 1:34) { m <- (a + b)/2; if (pred(m)) b <- m else a <- m }
    (a + b)/2
  }
  B_lo <- if (!is.na(i_lo) && i_lo > 1)
    refine(grid[i_lo - 1], grid[i_lo],
           function(B) length(count_roots(B, pp, ctx)) >= 3) else NA
  B_hi <- if (!is.na(i_hi) && i_hi > 1)
    refine(grid[i_hi - 1], grid[i_hi],
           function(B) { r <- count_roots(B, pp, ctx)
                         length(r) > 0 && max(r) < 0.5 }) else NA
  c(B_lo = B_lo, B_hi = B_hi)
}

par_list <- function() as.list(param(mod))

critical_load_table <- function() {
  pp <- par_list()
  out <- expand.grid(GA = c(25, 28, 32),
                     milk = c("FORMULA", "DONOR", "MOM"),
                     stringsAsFactors = FALSE)
  res <- t(mapply(function(g, m) {
    bifurcation_loads(pp, nec_context(pp, milk = m, GA = g))
  }, out$GA, out$milk))
  cbind(out, as.data.frame(res)) %>%
    mutate(window_ratio = B_hi/B_lo,
           reachable    = B_lo < pp$Ktot)
}

## =============================================================================
##  POPULATION SIMULATION
## =============================================================================
draw_insults <- function(GA, tmax = 42) {
  rate <- 0.09*(1 + 0.45*max(0, 30 - GA))
  tt <- numeric(0); t <- 0.5
  repeat {
    t <- t + rexp(1, rate); if (t >= tmax) break
    if (runif(1) <= exp(-t/20)) tt <- c(tt, t)      # front-loaded instability
  }
  ntx <- max(0, round(3.2 - 0.28*(GA - 24) + 0.8*rnorm(1)))
  ins <- list()
  if (length(tt)) ins$hyp <- data.frame(a = tt, b = tt + runif(length(tt), .1, .3),
                                        dep = runif(length(tt), .10, .30), nfk = 0)
  if (ntx > 0) { s <- 6 + runif(ntx)*(tmax - 12)
                 ins$tx <- data.frame(a = s, b = s + 0.5, dep = 0.16, nfk = 0.030) }
  d <- dplyr::bind_rows(ins)
  if (!nrow(d)) return(std_insults)
  d <- d[order(-d$dep), ][seq_len(min(3, nrow(d))), , drop = FALSE]
  z <- std_insults
  for (i in seq_len(nrow(d))) {
    z[[paste0("INS", i, "A")]]   <- d$a[i]
    z[[paste0("INS", i, "B")]]   <- d$b[i]
    z[[paste0("INS", i, "DEP")]] <- d$dep[i]
    z[[paste0("INS", i, "NFK")]] <- d$nfk[i]
  }
  z
}

nec_population <- function(n = 200, arm = list(), seed = 20260805, end = 42) {
  set.seed(seed)
  dplyr::bind_rows(lapply(seq_len(n), function(i) {
    GA <- 24 + 8*runif(1)^0.75
    BW <- max(0.42, 0.10*(GA - 22) + 0.35 + 0.09*rnorm(1))
    iiv <- list(Ktlr0 = exp(0.22*rnorm(1)), kA = exp(0.20*rnorm(1)),
                kE = exp(0.18*rnorm(1)), muB = exp(0.20*rnorm(1)),
                Pmax = exp(0.15*rnorm(1)), seedBenv = exp(0.45*rnorm(1)),
                kPAF = exp(0.25*rnorm(1)), dPAF0 = exp(0.25*rnorm(1)),
                seedChm = exp(1.00*rnorm(1)), seedCenv = exp(1.00*rnorm(1)))
    pp <- par_list()
    mult <- lapply(names(iiv), function(k) pp[[k]]*iiv[[k]])
    names(mult) <- names(iiv)
    a <- c(list(GA = GA, BW = BW,
                binfant = if (runif(1) < 0.35) 1.0 else 0.05,
                dysbio  = exp(0.35*rnorm(1)),
                insults = draw_insults(GA, end), end = end),
           mult, arm)
    do.call(run_nec, a) %>%
      summarise(ID = i, GA = GA,
                peak_bell = max(BellStage),
                nec = as.integer(max(BellStage) >= 2),
                surgical = as.integer(max(BellStage) >= 3),
                onset_day = ifelse(any(BellStage >= 2),
                                   min(time[BellStage >= 2]), NA_real_))
  }))
}

population_summary <- function(pop) {
  pop %>% summarise(n = n(), nec_pct = 100*mean(nec),
                    surgical_pct = 100*mean(surgical),
                    median_onset = median(onset_day, na.rm = TRUE))
}

## =============================================================================
##  THE FIVE EXPERIMENTS THE MODEL EXISTS TO RUN
## =============================================================================

## (1) two critical loads by gestational age and milk type
##       critical_load_table()

## (2) feeding advance is two-signed: the optimum is interior, and the slope is
##     flat for a 30-weeker and steep for a 25-weeker
feeding_sweep <- function(GA = 25, n = 150, rates = c(0, 10, 20, 30, 40)) {
  dplyr::bind_rows(lapply(rates, function(r) {
    pop <- nec_population(n, arm = list(GA = GA, milk = "FORMULA",
                                        feedrate = r), seed = 700 + r)
    population_summary(pop) %>%
      mutate(feed_rate = r,
             days_to_full = ifelse(r > 0, (160 - 20)/r + 2, Inf))
  }))
}

## (3) empirical antibiotic duration changes sign.  Run this in HUMAN-MILK-fed
##     infants on purpose: the harm is the loss of C, so it can only appear
##     where C exists.  In a formula-fed infant there is no commensal community
##     to destroy and antibiotics look protective — a falsifiable prediction.
abx_duration_sweep <- function(n = 200, days = c(0, 2, 3, 5, 7, 10, 14)) {
  dplyr::bind_rows(lapply(days, function(d) {
    pop <- nec_population(n, arm = list(GA = 26, milk = "MOM", abx_days = d),
                          seed = 4400 + d)
    population_summary(pop) %>% mutate(abx_days = d)
  }))
}

## (4) threshold-movers x state-movers multiply.  thrboost is a pure
##     threshold-mover; prob_dose is a pure state-mover; human milk is both.
factorial_arms <- function(n = 250) {
  arms <- list(
    "formula"                  = list(milk = "FORMULA"),
    "formula + probiotic"      = list(milk = "FORMULA", prob_dose = 1.10),
    "formula + threshold agent"= list(milk = "FORMULA", thrboost = 1.5),
    "formula + both"           = list(milk = "FORMULA", prob_dose = 1.10,
                                      thrboost = 1.5),
    "mother's own milk"        = list(milk = "MOM"))
  res <- dplyr::bind_rows(lapply(names(arms), function(k)
    population_summary(nec_population(n, arm = c(arms[[k]],
                                                 list(abx_days = 3)),
                                      seed = 9091)) %>% mutate(arm = k)))
  ref <- res$nec_pct[res$arm == "formula"]
  res %>% mutate(RR = nec_pct/ref)
}

## (5) the point of no return, verified by simulation rather than asserted
point_of_no_return <- function() {
  pp <- par_list()
  crit <- pp$INJth/pp$wNEC
  sim  <- do.call(run_nec, scenarios$S11$args)
  list(NECa_crit = crit,
       crossed_at = ifelse(any(sim$NECA > crit), min(sim$time[sim$NECA > crit]), NA),
       bell2 = ifelse(any(sim$BellStage >= 2), min(sim$time[sim$BellStage >= 2]), NA),
       bell3 = ifelse(any(sim$BellStage >= 3), min(sim$time[sim$BellStage >= 3]), NA))
}

## =============================================================================
##  PLOTS
## =============================================================================
plot_core_loop <- function(sim) {
  sim %>%
    select(time, scenario, E, Permeability, Translocat, Injury, NECA) %>%
    pivot_longer(-c(time, scenario)) %>%
    mutate(name = factor(name, levels = c("E", "Permeability", "Translocat",
                                          "Injury", "NECA"))) %>%
    ggplot(aes(time, value, colour = scenario)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~name, scales = "free_y", ncol = 1) +
    labs(x = "postnatal day", y = NULL,
         title = "The closed loop: E -> permeability -> translocation -> injury -> E") +
    theme_bw(base_size = 10)
}

plot_ecology <- function(sim) {
  sim %>%
    select(time, scenario, B, C, P, SCFA, HMO) %>%
    pivot_longer(-c(time, scenario)) %>%
    ggplot(aes(time, value, colour = scenario)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~name, scales = "free_y") +
    labs(x = "postnatal day", y = NULL, title = "Luminal ecology") +
    theme_bw(base_size = 10)
}

plot_field <- function(GA = 28, milk = "FORMULA", Bs = c(6, 12, 18, 26)) {
  pp <- par_list(); ctx <- nec_context(pp, milk = milk, GA = GA)
  grid <- expand.grid(E = seq(0.02, 0.99, by = 0.005), B = Bs)
  grid$g <- mapply(reduced_field, grid$E, grid$B,
                   MoreArgs = list(pp = pp, ctx = ctx))
  ggplot(grid, aes(E, g, colour = factor(B))) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(linewidth = 0.7) +
    labs(x = "villus integrity E", y = "g(E)  (dE/dt)",
         colour = "B (1e9/g)",
         title = sprintf("Reduced field at GA %g wk, %s: roots are the fates",
                         GA, milk)) +
    theme_bw(base_size = 10)
}

## =============================================================================
##  EXAMPLE SESSION
## =============================================================================
if (interactive()) {
  sim <- run_all_scenarios()
  print(summarise_scenarios(sim), n = 20)
  print(critical_load_table())
  print(point_of_no_return())
  print(feeding_sweep(GA = 25))
  print(feeding_sweep(GA = 30))
  print(abx_duration_sweep())
  print(factorial_arms())
  print(plot_field())
  print(plot_core_loop(dplyr::filter(sim, scenario %in% c("S1", "S2", "S5", "S11"))))
  print(plot_ecology(dplyr::filter(sim, scenario %in% c("S1", "S2", "S4", "S8"))))
}
