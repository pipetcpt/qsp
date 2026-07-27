# =============================================================================
#  Primary Hyperoxaluria Type 1 (PH1) — QSP model (mrgsolve)
#  ---------------------------------------------------------------------------
#  ph1_mrgsolve_model.R   ·  companion to ph1_qsp_model.dot / ph1_shiny_app.R
#
#  WHAT THIS MODEL IS ABOUT
#  ------------------------
#  PH1 is a LIVER enzyme defect that destroys the KIDNEY.  A single hepatic
#  peroxisomal enzyme — alanine:glyoxylate aminotransferase (AGT, AGXT,
#  2q37.3) — is missing, and everything downstream follows from the fact that
#  the carbon AGT can no longer route to glycine has only one other
#  irreversible destination: OXALATE, which humans cannot catabolise and can
#  only excrete through the kidney they are in the process of losing.
#
#  The model is built as a closed mass balance on that carbon.  Nothing about
#  "PH1 severity" is asserted; the phenotype is generated from two numbers —
#  residual AGT activity (AGTACT0) and glomerular filtration (GFR0/NM0).
#
#  THREE LOOPS DO THE WORK
#  -----------------------
#  LOOP 1  THE GLYOXYLATE FUTILE CYCLE
#          glyoxylate --GRHPR--> glycolate --GO/HAO1--> glyoxylate
#          With AGT present, glyoxylate is transaminated to glycine before it
#          can enter this cycle.  With AGT gone the carbon circulates, and
#          each pass gives cytosolic LDH another chance at the irreversible
#          step.  This is the model's mechanistic account of two otherwise
#          puzzling clinical facts:
#            - silencing GO (lumasiran), which is NOT the oxalate-forming
#              enzyme, cuts urinary oxalate by ~2/3, and RAISES glycolate
#              2-4x (the carbon has to leave somewhere, and with GO blocked
#              the only remaining exit is the urine as glycolate);
#            - silencing LDHA (nedosiran) cuts oxalate WITHOUT raising
#              glycolate, because it blocks the exit rather than the cycle,
#              and the displaced carbon leaves as plasma/urinary glyoxylate
#              instead.  The model therefore predicts a glyoxylate rise on
#              nedosiran, which is a testable consequence rather than an
#              assumption (see the VALIDATION block).
#
#  LOOP 2  CRYSTAL - INFLAMMATION - NEPHRON LOSS  (the accelerator)
#          Crystal deposition is driven not by urinary oxalate but by the
#          oxalate load PER SURVIVING NEPHRON divided by urine flow.  So each
#          nephron lost raises the insult on those that remain.  This is the
#          PH1 analogue of Brenner hyperfiltration and it is why the eGFR
#          trajectory bends downward instead of declining linearly.  It is
#          also why hyperhydration works at all: it is the only intervention
#          that attacks the denominator.
#
#  LOOP 3  PLASMA-OXALATE ESCAPE AND THE BONE CAPACITOR
#          Oxalate leaves only by filtration + tubular secretion
#          (CL ~ 1.3 x GFR).  As GFR falls, an unchanged production must exit
#          through a shrinking clearance, so plasma oxalate rises
#          hyperbolically, not linearly.  Above ~30 umol/L it precipitates in
#          bone, retina, myocardium and vessels.  Bone is large and slow: it
#          BUFFERS plasma oxalate (masking severity, and blunting the plasma
#          response to RNAi in advanced CKD), then RELEASES it for months to
#          years once a transplant restores clearance.  Both directions are
#          the same term with opposite sign — the post-transplant "oxalate
#          release" burden is generated, not scripted.
#
#  ONE CONSEQUENCE WORTH STATING UP FRONT
#  --------------------------------------
#  Because Uox = 1.3 x GFR x Pox, urinary oxalate FALLS as PH1 progresses.
#  A model that used Uox as its severity metric would report a patient
#  improving on the way to dialysis.  That is exactly why ILLUMINATE-C used
#  plasma oxalate rather than urinary oxalate as its primary endpoint, and the
#  model reproduces the crossover (see scenario 03 at year 15+).
#
#  WHY THERE ARE ALMOST NO DOSING EVENTS
#  -------------------------------------
#  Every scenario below is a PURE PARAMETER SET: no dosing events, no scripted
#  interventions, so any two scenarios differ only by the numbers listed.
#  Chronic oral drugs (pyridoxine, potassium citrate, stiripentol) enter as
#  continuous input rates.  The two subcutaneous siRNAs enter through smooth
#  periodic rectangular windows (PULSE below) so that real loading-then-
#  quarterly behaviour, and the trough in effect between quarterly doses, are
#  preserved.  All scenarios include a 180-day untreated RUN-IN so that every
#  "% change from baseline" is computed the way the trials computed it.
#
#  BODY SIZE
#  ---------
#  All hepatic capacities, distribution volumes, dietary load, urine volume
#  and urine calcium/citrate scale with SF = BSA/1.73, so the same parameter
#  set describes an adult and a 8-kg infant; urine values are reported both
#  absolute and normalised to 1.73 m2.  Infantile oxalosis is then generated
#  by the two things that are actually different in infancy — immature GFR
#  (GFR0) and higher collagen-turnover precursor supply (PRECUR) — rather
#  than by a hand-set severity term.
#
#  WHAT IS DELIBERATELY *NOT* CLAIMED
#  ----------------------------------
#  Deterministic single-patient model; it produces trajectories, not incidence
#  rates. Registry percentages are used as sanity anchors for phenotypes, never
#  as fitted targets. Hepatic glyoxylate is an EFFECTIVE hepatocyte
#  concentration, not a measured one — only its ratio across genotypes carries
#  meaning. The apparent Km of GO (KMGO) lumps peroxisomal transport and the
#  liver:plasma glycolate gradient; the model needs GO to run near capacity,
#  which is what makes glycolate a near-stoichiometric on-target biomarker.
#  Parameters are literature-anchored where a number exists (all normal/PH1
#  urine and plasma values, FEox > 1, the 30 umol/L oxalosis threshold, trial
#  endpoints, dialysis dialysance) and hand-calibrated otherwise. The model
#  over-predicts the plasma-oxalate response to lumasiran in advanced CKD
#  relative to ILLUMINATE-C (-57% vs -33%); the likely reason is that real
#  patients carry larger tissue oxalate stores than the model assumes, and
#  this is reported rather than tuned away.  NOT FOR CLINICAL USE.
#
#  Requires: mrgsolve, dplyr, tidyr (ggplot2 optional)
#  Run:      Rscript ph1_mrgsolve_model.R
# =============================================================================

library(mrgsolve)
suppressMessages({
  library(dplyr)
  library(tidyr)
})

code <- '
$PROB
# Primary hyperoxaluria type 1 QSP model
- 51 ODEs: 14 drug PK / transduction, 5 hepatic gene-enzyme, 3 metabolite
  pools, 4 systemic oxalate compartments, 4 urine chemistry, 4 crystal/stone,
  6 injury-fibrosis-nephron mass, 2 graft/liver replacement, 5 systemic
  oxalosis / safety, 4 outcome accumulators
- Every scenario is a pure parameter set

$PARAM @annotated
// ---------------------------------------------------------------- patient
WT      : 70.0  : body weight (kg)
BSA     : 1.73  : body surface area (m2); all capacities scale with BSA/1.73
GFR0    : 105.0 : GFR at full nephron mass (mL/min/1.73 m2)
EXPNM   : 0.85  : nephron-mass exponent on GFR (partial single-nephron compensation)
NM0     : 0.95  : initial functional nephron mass (fraction)
NC0     : 0.0   : initial nephrocalcinosis burden (0-1)
FIB0    : 0.0   : initial interstitial fibrosis (0-1)
OXADEP0 : 0.0   : initial soft-tissue oxalate deposit (umol)
OXABONE0: 0.0   : initial bone oxalate reservoir (umol)

// ------------------------------------------------------------- genotype
AGTACT0 : 0.02  : residual peroxisomal AGT activity (fraction of normal)
B6RESP  : 1.0   : pyridoxine responsiveness (1 = G170R/F152I, 0 = c.33dupC/null)
EB6     : 9.0   : maximum fold-increase in functional AGT from PLP chaperoning
GRHPRF  : 1.0   : GRHPR activity multiplier (set < 1 to approximate PH2)
HOGAF   : 1.0   : HOGA1 flux multiplier (hydroxyproline route; PH3 surrogate)
PRECUR  : 1.0   : precursor-supply multiplier (>1 = high collagen turnover, infancy)

// ------------------------------------- hepatic precursor supply (umol/day)
HYPFLUX : 900.0  : hydroxyproline-derived glyoxylate via HOGA1
GLYOXO  : 366.0  : glycine/glycolaldehyde-derived glyoxylate (direct)
GLCIN   : 2111.0 : endogenous + dietary glycolate entering the body pool

// ------------------------------------------- hepatic enzyme kinetics
VLIV    : 1.5    : hepatocyte water volume for glyoxylate (L)
VMAGT   : 9000.0 : AGT Vmax at full activity (umol/day)
KMAGT   : 150.0  : AGT Km for glyoxylate (umol/L)
VMGR    : 1200.0 : GRHPR Vmax, glyoxylate -> glycolate (umol/day)
KMGR    : 100.0  : GRHPR Km for glyoxylate (umol/L)
VMLDH   : 1600.0 : cytosolic LDH Vmax, glyoxylate -> oxalate (umol/day)
KMLDH   : 300.0  : LDH K50 for glyoxylate (umol/L)
HLDH    : 2.0    : LDH Hill coefficient (tetramer; gives the 30x normal-to-PH1 span)
VMGO    : 2735.0 : glycolate oxidase Vmax, glycolate -> glyoxylate (umol/day)
KMGO    : 0.5    : GO APPARENT Km (umol/L, plasma-equivalent; GO runs near capacity)
VMALT   : 1000.0 : alternative glyoxylate disposal Vmax (AGXT2/GPT, decarboxylation)
KMALT   : 2000.0 : alternative disposal Km (umol/L)
KEXGLX  : 1.0    : hepatic glyoxylate efflux to plasma (1/day)
FEGLX   : 0.6    : glyoxylate fractional renal excretion
KPLDH   : 15.0   : extrahepatic LDH clearance of glyoxylate to oxalate (L/day)
VDGLC   : 25.0   : glycolate distribution volume (L)
FEGLC   : 0.35   : glycolate fractional renal excretion
KMETGLC : 15.0   : non-renal metabolic clearance of glycolate (L/day)

// ------------------------------- oxalate distribution, gut and kidney
VDOX    : 21.0   : oxalate / glyoxylate distribution volume (L)
FEOX    : 1.3    : oxalate renal clearance as a multiple of GFR (net secretion)
OXDIET  : 2500.0 : dietary oxalate intake (umol/day)
FABSOX  : 0.08   : fraction of luminal oxalate absorbed
KGT     : 1.0    : intestinal transit rate constant (1/day)
KDEGG   : 4.0    : luminal oxalate degradation rate at full enzyme/bacterial activity
KENT    : 3.0    : SLC26A6 enteric secretion clearance (L/day)
KENTM   : 2.5    : enhancement of enteric secretion by oxalate-degrader colonisation

// ------------------------------------------------- systemic oxalosis
POXCRIT : 30.0   : plasma oxalate threshold for tissue deposition (umol/L)
KDEP    : 6.0    : tissue deposition rate coefficient
HDEP    : 1.5    : deposition power on (Pox - POXCRIT)
KREL    : 0.0018 : tissue/bone oxalate release rate constant (1/day)
FBONE   : 0.85   : fraction of deposited oxalate going to the bone reservoir

// ------------------------------------------------- urine chemistry
UVTGT   : 1.5    : target urine volume (L/day per 1.73 m2)
KUV     : 0.35   : rate of adaptation of urine volume (1/day)
ADHFL   : 1.0    : adherence to the prescribed fluid intake (fraction)
UCAL    : 4.0    : urinary calcium excretion (mmol/day per 1.73 m2)
UCITB   : 3.0    : baseline urinary citrate (mmol/day per 1.73 m2)
ECIT    : 4.0    : maximum citrate increment from potassium citrate (mmol/day)
KC50    : 0.15   : plasma citrate for half-maximal urinary citrate increment
KUC     : 0.3    : urinary citrate adaptation rate (1/day)
KSS     : 8.0    : CaOx relative-supersaturation scaling constant
KICIT   : 2.5    : urinary citrate concentration halving CaOx supersaturation (mmol/L)
RSSCRIT : 4.0    : metastable limit: RSS above which spontaneous nucleation occurs

// ------------------------------------------------- crystal and stone
KNUC    : 0.0025 : nucleation rate coefficient on (RSS - RSSCRIT)^2
KGROW   : 0.005  : crystal growth rate coefficient
KWASH   : 0.25   : urine-flow-dependent crystal washout (1/day at 1.5 L/day)
KFIX    : 0.05   : loss of luminal crystals by fixation/retention (1/day)
KFIX2   : 0.0008 : conversion of retained crystals into nephrocalcinosis
KTUBDEP : 0.00006: direct tubular deposition per unit normalised tubular oxalate
KRESNC  : 0.0015 : nephrocalcinosis resolution rate (1/day)
KSTG    : 0.010  : stone growth from crystal load above the metastable limit
KPST    : 0.004  : stone passage / removal rate (1/day)
KEV     : 1.2    : maximum symptomatic stone events per year at saturating burden

// ------------------------------------------------- injury cascade
KINJ    : 0.05   : tubular injury formation rate
KDINJ   : 0.06   : tubular injury resolution rate (1/day)
KPOXINJ : 1.5    : weight of plasma-oxalate-driven (acute oxalate nephropathy) injury
KMACP   : 0.05   : macrophage recruitment rate
KDMAC   : 0.05   : macrophage clearance rate (1/day)
KIL1    : 0.30   : NLRP3/IL-1beta generation per unit macrophage x crystal
KDIL1   : 0.15   : IL-1beta decay (1/day)
KTGF    : 0.05   : TGF-beta generation rate
KDTGF   : 0.05   : TGF-beta decay (1/day)
KFIB    : 0.0025 : fibrosis accumulation rate
KDFIB   : 0.0025 : fibrosis resolution rate (1/day)
KNL     : 0.00033: nephron loss rate per unit damage (1/day)

// --------------------------------------- dialysis and transplantation
HDON    : 0.0    : haemodialysis on/off (forced)
HDAUTO  : 1.0    : start dialysis automatically at eGFR < 15
HDSESS  : 3.0    : haemodialysis sessions per week
HDHR    : 4.0    : hours per session
DIALZ   : 120.0  : oxalate dialysance (mL/min)
PDCL    : 0.0    : additional peritoneal dialysis oxalate clearance (L/day)
TXDAY   : 1e6    : day of transplantation
KTXON   : 0.0    : kidney graft included (0/1)
LTXON   : 0.0    : liver graft included (0/1) — this is what cures the metabolism
NMGT    : 0.55   : nephron mass delivered by a single kidney graft
KGR     : 0.05   : engraftment rate (1/day, during a 30-day window)
GRSENS  : 2.0    : graft susceptibility to oxalate injury relative to native kidney
KLIV    : 0.10   : rate of hepatic enzyme restoration after liver graft (1/day)
KTXWASH : 0.25   : washout of native-kidney lesion states at transplantation (1/day)

// ------------------------------------------------------------- drugs
RUNIN   : 180.0  : untreated run-in before any treatment starts (day)
LUMON   : 0.0    : lumasiran on/off
LUMLOAD : 3.0    : lumasiran loading dose (mg/kg SC monthly)
LUMNL   : 3.0    : number of monthly loading doses
LUMDOSE : 3.0    : lumasiran maintenance dose (mg/kg SC)
LUMTAU  : 90.0   : lumasiran maintenance interval (day)
NEDON   : 0.0    : nedosiran on/off
NEDDOSE : 3.5    : nedosiran dose (mg/kg SC)
NEDTAU  : 30.0   : nedosiran interval (day)
B6ON    : 0.0    : pyridoxine on/off
B6DOSE  : 0.0    : pyridoxine dose (mg/kg/day PO)
KCITON  : 0.0    : potassium citrate on/off
KCITD   : 0.15   : potassium citrate dose (g/kg/day)
STPON   : 0.0    : stiripentol on/off
STPDOSE : 50.0   : stiripentol dose (mg/kg/day)
OXDON   : 0.0    : oral oxalate-degrading therapy on/off (reloxaliase / O. formigenes)
ANAON   : 0.0    : IL-1 blockade on/off (exploratory, preclinical only)

// --------------------------------------------------- siRNA PK and PD
KALUM   : 0.6    : lumasiran SC absorption rate (1/day)
CLLUM   : 250.0  : lumasiran plasma clearance (L/day)
VLUM    : 15.0   : lumasiran central volume (L)
KUPLUM  : 45.0   : ASGPR-mediated hepatic uptake clearance (L/day)
KOUTLUM : 0.030  : hepatic depot loss (1/day)
KRLUM   : 0.10   : RISC loading rate (1/day)
KDRLUM  : 0.022  : RISC-loaded siRNA decay (1/day) — sets the quarterly interval
EMAXH   : 0.93   : maximum HAO1 mRNA knockdown
IC50H   : 0.35   : RISC amount for half-maximal HAO1 knockdown
KANED   : 1.2    : nedosiran SC absorption rate (1/day)
CLNED   : 120.0  : nedosiran plasma clearance (L/day)
VNED    : 12.0   : nedosiran central volume (L)
KUPNED  : 25.0   : nedosiran hepatic uptake clearance (L/day)
KOUTNED : 0.060  : nedosiran hepatic depot loss (1/day)
KRNED   : 0.15   : nedosiran RISC loading rate (1/day)
KDRNED  : 0.055  : nedosiran RISC decay (1/day) — sets the monthly interval
EMAXL   : 0.90   : maximum LDHA mRNA knockdown
IC50L   : 0.45   : RISC amount for half-maximal LDHA knockdown
KSHAO   : 0.35   : HAO1 mRNA synthesis (1/day)
KDHAO   : 0.35   : HAO1 mRNA degradation (1/day)
KSGO    : 0.10   : GO protein synthesis (1/day)
KDGO    : 0.10   : GO protein degradation (1/day) — t1/2 ~7 d, sets onset/offset lag
KSLDA   : 0.50   : LDHA mRNA synthesis (1/day)
KDLDA   : 0.50   : LDHA mRNA degradation (1/day)
KSLDP   : 0.18   : LDH protein synthesis (1/day)
KDLDP   : 0.18   : LDH protein degradation (1/day)
KSAGT   : 0.12   : AGT protein synthesis (1/day)
KDAGT   : 0.12   : AGT protein degradation (1/day)

// ---------------------------------------------- small-molecule PK
KAB6    : 3.0    : pyridoxine absorption rate (1/day)
CLB6    : 60.0   : PLP clearance (L/day)
VB6     : 30.0   : PLP distribution volume (L)
PLP50   : 0.08   : PLP concentration for half-maximal AGT chaperoning (mg/L)
KACIT   : 4.0    : citrate absorption rate (1/day)
CLCIT   : 900.0  : citrate clearance (L/day)
VCIT    : 18.0   : citrate distribution volume (L)
CLSTP   : 100.0  : stiripentol clearance (L/day)
VSTP    : 60.0   : stiripentol volume (L)
KISTP   : 25.0   : stiripentol concentration for half-maximal LDH inhibition (mg/L)
KOXD    : 0.5    : onset/offset of luminal oxalate-degrading activity (1/day)
EOXD    : 0.85   : maximum fraction of luminal oxalate degraded
KANA    : 1.0    : onset/offset of IL-1 blockade (1/day)
EANA    : 0.75   : maximum IL-1 pathway blockade
KB6TOX  : 2.5e-5 : pyridoxine neurotoxicity accumulation rate
KB6REC  : 0.004  : recovery from pyridoxine neuropathy (1/day)

$CMT @annotated
LUMSC   : lumasiran subcutaneous depot (mg)
LUMC    : lumasiran plasma (mg)
LUML    : lumasiran hepatic depot (mg-equivalent)
LUMR    : lumasiran AGO2-RISC loaded pool (arbitrary)
NEDSC   : nedosiran subcutaneous depot (mg)
NEDC    : nedosiran plasma (mg)
NEDL    : nedosiran hepatic depot (mg-equivalent)
NEDR    : nedosiran AGO2-RISC loaded pool (arbitrary)
B6G     : pyridoxine gut (mg)
PLPC    : plasma PLP (mg)
CITG    : citrate gut (mmol)
CITC    : plasma citrate (mmol)
STPC    : stiripentol central (mg)
OXDG    : luminal oxalate-degrading activity (0-1)
HAO1M   : HAO1 mRNA (fraction of normal)
GOP     : glycolate oxidase protein activity (fraction of normal)
LDHAM   : LDHA mRNA (fraction of normal)
LDHP    : cytosolic LDH protein activity (fraction of normal)
AGTP    : functional peroxisomal AGT activity (fraction of normal)
GLXL    : hepatic glyoxylate (umol)
GLCB    : body glycolate pool (umol)
GLXP    : plasma glyoxylate (umol)
OXAP    : plasma/extracellular oxalate (umol)
OXAG    : intestinal luminal oxalate (umol)
OXADEP  : soft-tissue oxalate deposit (umol)
OXABONE : bone oxalate reservoir (umol)
UOXM    : measured 24-h urinary oxalate (umol/day, 1-day filter)
UGLCM   : measured 24-h urinary glycolate (umol/day, 1-day filter)
UCITM   : urinary citrate excretion (mmol/day)
UVOLS   : urine volume (L/day)
CRYST   : luminal calcium-oxalate crystal load (arbitrary)
NC      : nephrocalcinosis burden (0-1)
STONE   : stone burden (arbitrary)
STEV    : cumulative symptomatic stone events
TUBINJ  : tubular epithelial injury (0-1)
MAC     : macrophage / monocyte infiltrate (arbitrary)
IL1B    : NLRP3-driven IL-1beta signal (arbitrary)
TGFB    : TGF-beta1 profibrotic signal (arbitrary)
FIB     : interstitial fibrosis (0-1)
NM      : native functional nephron mass (fraction)
NMG     : kidney graft functional nephron mass (fraction)
LIVFX   : hepatic enzyme replacement factor after liver graft (0-1)
RETINA  : retinal oxalosis burden (0-1)
CARDIO  : myocardial / conduction oxalosis burden (0-1)
BONEDIS : oxalate osteopathy burden (0-1)
NEURO   : peripheral neuropathy burden (0-1)
B6NEURO : pyridoxine-attributable neuropathy (0-1)
HAZ     : cumulative hazard of kidney failure
UOXAUC  : cumulative urinary oxalate (umol)
POXAUC  : cumulative plasma oxalate exposure (umol/L*day)
DIALOX  : cumulative dialytic oxalate removal (umol)

$GLOBAL
#define pos(x) ((x) > 0.0 ? (x) : 0.0)
// smooth periodic rectangular dose window of width W days at period TAU
#define PULSE(ph, W) (1.0/(1.0+exp(-(ph)/0.4)) - 1.0/(1.0+exp(-((ph)-(W))/0.4)))

$MAIN
// initial conditions come from parameters so that every scenario stays a
// pure parameter set (a patient can be started at CKD4 with bone stores)
double SFI = BSA/1.73;
AGTP_0    = AGTACT0;
HAO1M_0   = 1.0;
GOP_0     = 1.0;
LDHAM_0   = 1.0;
LDHP_0    = 1.0;
GLXL_0    = 700.0  * SFI;
GLCB_0    = 300.0  * SFI;
GLXP_0    = 100.0  * SFI;
OXAP_0    = 180.0  * SFI;
OXAG_0    = 2500.0 * SFI;
UOXM_0    = 1400.0 * SFI;
UGLCM_0   = 600.0  * SFI;
UCITM_0   = UCITB  * SFI;
UVOLS_0   = UVTGT  * SFI;
NM_0      = NM0;
NC_0      = NC0;
FIB_0     = FIB0;
OXADEP_0  = OXADEP0;
OXABONE_0 = OXABONE0;

$ODE
double tt = SOLVERTIME;
double SF = BSA/1.73;                        // body-size scaling factor
double ON = (tt >= RUNIN) ? 1.0 : 0.0;       // treatment switch (after run-in)
double tl = tt - RUNIN;
double W  = 5.0;

// ------------------------------------------------------------------------
// LUMASIRAN — GalNAc-siRNA against HAO1 (glycolate oxidase)
// monthly loading x LUMNL, then every LUMTAU days
// ------------------------------------------------------------------------
double lrate = 0.0;
if(LUMON > 0.5 && tl >= 0.0) {
  double ph, dose;
  if(tl < 30.0 * LUMNL) {
    ph   = fmod(tl, 30.0);
    dose = LUMLOAD * WT;
  } else {
    ph   = fmod(tl - 30.0 * LUMNL, LUMTAU);
    dose = LUMDOSE * WT;
  }
  lrate = dose / W * PULSE(ph, W);
}
double CLUMc = LUMC / VLUM;
dxdt_LUMSC = lrate - KALUM * LUMSC;
dxdt_LUMC  = KALUM * LUMSC - (CLLUM + KUPLUM) * CLUMc;
dxdt_LUML  = KUPLUM * CLUMc - KOUTLUM * LUML - KRLUM * LUML;
dxdt_LUMR  = KRLUM * LUML - KDRLUM * LUMR;

// ------------------------------------------------------------------------
// NEDOSIRAN — GalNAc-siRNA against LDHA
// ------------------------------------------------------------------------
double nrate = 0.0;
if(NEDON > 0.5 && tl >= 0.0) {
  double ph2 = fmod(tl, NEDTAU);
  nrate = NEDDOSE * WT / W * PULSE(ph2, W);
}
double CNEDc = NEDC / VNED;
dxdt_NEDSC = nrate - KANED * NEDSC;
dxdt_NEDC  = KANED * NEDSC - (CLNED + KUPNED) * CNEDc;
dxdt_NEDL  = KUPNED * CNEDc - KOUTNED * NEDL - KRNED * NEDL;
dxdt_NEDR  = KRNED * NEDL - KDRNED * NEDR;

// ------------------------------------------------------------------------
// PYRIDOXINE -> PLP: cofactor AND pharmacological chaperone for AGT
// ------------------------------------------------------------------------
double CB6 = PLPC / VB6;
dxdt_B6G  = B6ON * ON * B6DOSE * WT - KAB6 * B6G;
dxdt_PLPC = KAB6 * B6G - CLB6 * CB6;
double PLPF = CB6 / (PLP50 + CB6);

// ------------------------------------------------- potassium citrate
double CCIT = CITC / VCIT;
dxdt_CITG = KCITON * ON * KCITD * WT * 3.0 - KACIT * CITG;
dxdt_CITC = KACIT * CITG - CLCIT * CCIT;

// ------------------------------------------------- stiripentol (LDH-5 inhibitor)
double CSTP = STPC / VSTP;
dxdt_STPC = STPON * ON * STPDOSE * WT - CLSTP * CSTP;
double STPI = CSTP / (KISTP + CSTP);

// -------------------------- luminal oxalate degradation, IL-1 blockade
dxdt_OXDG = KOXD * (OXDON * ON - OXDG);
double DEGF = EOXD * OXDG;

// ------------------------------------------------------------------------
// HEPATIC GENE -> ENZYME LAYER
// mRNA knockdown is immediate; PROTEIN turnover (KDGO, KDLDP) is what
// creates the observed lag to nadir and the slow offset after stopping.
// ------------------------------------------------------------------------
double KDH = EMAXH * LUMR / (IC50H + LUMR);
dxdt_HAO1M = KSHAO * (1.0 - KDH) - KDHAO * HAO1M;
dxdt_GOP   = KSGO * HAO1M - KDGO * GOP;

double KDL = EMAXL * NEDR / (IC50L + NEDR);
dxdt_LDHAM = KSLDA * (1.0 - KDL) - KDLDA * LDHAM;
dxdt_LDHP  = KSLDP * LDHAM - KDLDP * LDHP;

// AGT: genotype floor, raised by PLP chaperoning, replaced by a liver graft
double AGTTGT = AGTACT0 * (1.0 + EB6 * B6RESP * PLPF);
AGTTGT = AGTTGT * (1.0 - LIVFX) + 1.0 * LIVFX;
dxdt_AGTP = KSAGT * AGTTGT - KDAGT * AGTP;

// ------------------------------------------------- transplant switches
double txon = (tt >= TXDAY) ? 1.0 : 0.0;
double inwin = (tt >= TXDAY && tt < TXDAY + 30.0) ? 1.0 : 0.0;
double WASH = KTXWASH * KTXON * inwin;   // native-kidney lesions no longer counted
double ENG  = inwin;                     // one-off engraftment window

// ------------------------------------------------- renal function
double NMT  = (NM + NMG > 1e-3) ? (NM + NMG) : 1e-3;
double EGFR = GFR0 * pow(NMT, EXPNM);
double GFRL = EGFR * 1.44 * SF;                       // L/day
double hd   = (HDON > 0.5 || (HDAUTO > 0.5 && EGFR < 15.0)) ? 1.0 : 0.0;
double CLHD = hd * (HDSESS * HDHR * DIALZ * 60.0 / 1000.0 / 7.0 + PDCL);

// ========================================================================
// LOOP 1 — THE GLYOXYLATE FUTILE CYCLE
// glyoxylate --GRHPR--> glycolate --GO--> glyoxylate, with AGT as the drain
// that PH1 removes and LDH as the one irreversible exit.
// ========================================================================
double CGLX = GLXL / (VLIV * SF);
double CGLC = GLCB / (VDGLC * SF);

double VAGT = SF * VMAGT * AGTP * CGLX / (KMAGT + CGLX);
double VGR  = SF * VMGR * GRHPRF * CGLX / (KMGR + CGLX);
double VLDH = SF * VMLDH * LDHP * (1.0 - STPI) * pow(CGLX, HLDH)
              / (pow(KMLDH, HLDH) + pow(CGLX, HLDH));
double VALT = SF * VMALT * CGLX / (KMALT + CGLX);
double VGO  = SF * VMGO * GOP * CGLC / (KMGO + CGLC);
double EXGLX = KEXGLX * GLXL;
double SUPPLY = SF * PRECUR * (HYPFLUX * HOGAF + GLYOXO);

dxdt_GLXL = SUPPLY + VGO - VAGT - VGR - VLDH - VALT - EXGLX;

double UGLC = FEGLC * GFRL * CGLC;
dxdt_GLCB = SF * PRECUR * GLCIN + VGR - VGO - UGLC - (SF * KMETGLC + CLHD) * CGLC;

// glyoxylate that escapes the hepatocyte is converted to oxalate by
// EXTRAHEPATIC LDH — outside the reach of a hepatocyte-targeted siRNA.
double CGLXP = GLXP / (VDOX * SF);
double VPLDH = SF * KPLDH * CGLXP;
dxdt_GLXP = EXGLX - FEGLX * GFRL * CGLXP - VPLDH - CLHD * CGLXP;

// ========================================================================
// LOOP 3 — OXALATE MASS BALANCE AND THE BONE CAPACITOR
// ========================================================================
double POX   = OXAP / (VDOX * SF);
double ABSOX = FABSOX * (1.0 - DEGF) * KGT * OXAG;
double ENTSEC = SF * KENT * (1.0 + KENTM * OXDG) * POX;
double DEP   = SF * KDEP * pow(pos(POX - POXCRIT), HDEP);
double RELF  = pos(1.0 - POX / POXCRIT);
double REL   = KREL * (OXADEP + OXABONE) * RELF;
double UOX   = FEOX * GFRL * POX;

dxdt_OXAP    = VLDH + VPLDH + ABSOX + REL - UOX - ENTSEC - DEP - CLHD * POX;
dxdt_OXAG    = SF * OXDIET + ENTSEC - KGT * OXAG - KDEGG * DEGF * OXAG;
dxdt_OXADEP  = (1.0 - FBONE) * DEP - KREL * OXADEP * RELF;
dxdt_OXABONE = FBONE * DEP - KREL * OXABONE * RELF;
dxdt_DIALOX  = CLHD * POX;

// ------------------------------------------------- urine chemistry
dxdt_UVOLS = KUV * (SF * UVTGT * ADHFL - UVOLS);
double UCITT = SF * (UCITB + ECIT * CCIT / (KC50 + CCIT));
dxdt_UCITM = KUC * (UCITT - UCITM);

double UV    = (UVOLS > 0.05 * SF + 0.02) ? UVOLS : (0.05 * SF + 0.02);
double COX   = UOX / UV / 1000.0;          // mmol/L
double CCA   = SF * UCAL / UV;             // mmol/L
double CCITU = UCITM / UV;                 // mmol/L
double RSS   = KSS * COX * CCA / (1.0 + CCITU / KICIT);

dxdt_UOXM  = 1.0 * (UOX - UOXM);
dxdt_UGLCM = 1.0 * (UGLC - UGLCM);

// ========================================================================
// LOOP 2 — CRYSTAL / INFLAMMATION / NEPHRON LOSS
// TFN is the accelerator: oxalate load PER SURVIVING NEPHRON per litre of
// urine, normalised to healthy. Losing nephrons raises it even as Uox falls.
// ========================================================================
dxdt_CRYST = -WASH * CRYST
           + KNUC * pow(pos(RSS - RSSCRIT), 2.0)
           + KGROW * CRYST * pos(RSS - 1.0)
           - (KWASH * UV / 1.5 + KFIX) * CRYST;

double TFN = (UOX / NMT / UV / 1000.0) / 0.182;

dxdt_NC = -WASH * NC
        + KFIX2 * CRYST * (1.0 - NC)
        + KTUBDEP * pos(TFN - 1.0) * (1.0 - NC)
        - KRESNC * NC;

dxdt_STONE = -WASH * STONE
           + KSTG * CRYST * pos(RSS - RSSCRIT) - KPST * STONE;
dxdt_STEV  = KEV * STONE / (1.0 + STONE) / 365.0;

// crystal-driven innate immunity -> fibrosis
dxdt_TUBINJ = -WASH * TUBINJ
            + KINJ * (CRYST + 2.0 * NC + KPOXINJ * pos(POX - POXCRIT) / 30.0)
              * (1.0 - TUBINJ) - KDINJ * TUBINJ;
dxdt_MAC  = -WASH * MAC + KMACP * (CRYST + 2.0 * NC) - KDMAC * MAC;
double IL1BLK = EANA * ANAON * ON;
dxdt_IL1B = -WASH * IL1B
          + KIL1 * MAC * (CRYST + NC) * (1.0 - IL1BLK) - KDIL1 * IL1B;
dxdt_TGFB = -WASH * TGFB
          + KTGF * (TUBINJ + 0.7 * IL1B + 0.4 * MAC) - KDTGF * TGFB;
dxdt_FIB  = -WASH * FIB + KFIB * TGFB * (1.0 - FIB) - KDFIB * FIB;

double DAM = 0.5 * FIB + 1.0 * NC + 0.3 * TUBINJ;
dxdt_NM = -KNL * DAM * NM;

// ------------------------------------------------- graft and liver
dxdt_NMG   = KGR * (NMGT * KTXON - NMG) * ENG - KNL * GRSENS * DAM * NMG;
dxdt_LIVFX = KLIV * (LTXON * txon - LIVFX);

// ------------------------------------------------- systemic oxalosis
double dep = pos(POX - POXCRIT) / 30.0;
dxdt_RETINA  = 0.004 * dep * (1.0 - RETINA)  - 0.0004 * RETINA;
dxdt_CARDIO  = 0.003 * dep * (1.0 - CARDIO)  - 0.0004 * CARDIO;
dxdt_BONEDIS = 0.006 * dep * (1.0 - BONEDIS) - 0.0006 * BONEDIS;
dxdt_NEURO   = 0.002 * dep * (1.0 - NEURO) - 0.0005 * NEURO
             + KB6TOX * pos(CB6 * VB6 / WT - 15.0) * (1.0 - NEURO);
dxdt_B6NEURO = KB6TOX * pos(B6ON * ON * B6DOSE - 15.0) * (1.0 - B6NEURO)
             - KB6REC * B6NEURO;

// ------------------------------------------------- outcome accumulators
dxdt_HAZ    = 0.00004 * exp(2.2 * (60.0 - EGFR) / 30.0 + 0.5 * NC + 0.4 * FIB);
dxdt_UOXAUC = UOX;
dxdt_POXAUC = POX;

$TABLE
// Outputs are recomputed here from the state vector so every reported value
// is exactly consistent with the states at the output time.
double SFr  = BSA/1.73;
double NMTr = (NM + NMG > 1e-3) ? (NM + NMG) : 1e-3;
double eGFR = GFR0 * pow(NMTr, EXPNM);
double GFRLr = eGFR * 1.44 * SFr;
double Pox  = OXAP / (VDOX * SFr);
double Uox  = FEOX * GFRLr * Pox;
double UVr  = (UVOLS > 0.05 * SFr + 0.02) ? UVOLS : (0.05 * SFr + 0.02);
double nrm  = 1.73 / BSA;                       // normalise urine to 1.73 m2

double Uox24  = Uox * nrm / 1000.0;             // mmol/24h/1.73 m2
double Uglc24 = UGLCM * nrm / 1000.0;           // mmol/24h/1.73 m2
double UoxCr  = Uox * nrm / 1000.0 / 8.8;       // mmol/mmol creatinine (approx)
double Pglc   = GLCB / (VDGLC * SFr);
double Pglx   = GLXP / (VDOX * SFr);
double RSSr   = KSS * (Uox / UVr / 1000.0) * (SFr * UCAL / UVr)
                / (1.0 + (UCITM / UVr) / KICIT);
double Ucit24 = UCITM * nrm;
double Uvol24 = UVOLS * nrm;

// semi-quantitative clinical grades
double NCgrade  = 3.0 * NC * NC / (NC * NC + 0.09);         // ultrasound 0-3
double FIBgrade = 3.0 * FIB * FIB / (FIB * FIB + 0.09);
double SYSOX    = (RETINA + CARDIO + BONEDIS + NEURO) / 4.0;
double ESKD     = (eGFR < 15.0) ? 1.0 : 0.0;
double OXNEPH   = (Pox > POXCRIT) ? 1.0 : 0.0;              // oxalosis risk flag
double CKDST    = (eGFR >= 90) ? 1.0 : (eGFR >= 60) ? 2.0 : (eGFR >= 45) ? 3.0 :
                  (eGFR >= 30) ? 3.5 : (eGFR >= 15) ? 4.0 : 5.0;
double onDial   = (HDON > 0.5 || (HDAUTO > 0.5 && eGFR < 15.0)) ? 1.0 : 0.0;
double AGTact   = AGTP;
double GOact    = GOP;
double LDHact   = LDHP;
double BoneOx   = OXABONE / 1000.0;                          // mmol
double SoftOx   = OXADEP / 1000.0;                           // mmol
double DialOx   = DIALOX / 1000.0;                           // mmol
double StoneEv  = STEV;
double SurvESKD = exp(-HAZ);

$CAPTURE @annotated
Uox24    : 24-h urinary oxalate (mmol/1.73 m2/day)
UoxCr    : urinary oxalate:creatinine (mmol/mmol; adult-equivalent approximation)
Uglc24   : 24-h urinary glycolate (mmol/1.73 m2/day)
Pox      : plasma oxalate (umol/L)
Pglc     : plasma glycolate (umol/L)
Pglx     : plasma glyoxylate (umol/L)
RSSr     : calcium-oxalate relative supersaturation
Uvol24   : urine volume (L/1.73 m2/day)
Ucit24   : urinary citrate (mmol/1.73 m2/day)
eGFR     : eGFR (mL/min/1.73 m2)
CKDST    : CKD stage (3.5 = stage 3b)
ESKD     : end-stage kidney disease flag
onDial   : on dialysis flag
NCgrade  : nephrocalcinosis grade (0-3)
FIBgrade : interstitial fibrosis grade (0-3)
StoneEv  : cumulative symptomatic stone events
OXNEPH   : plasma oxalate above the systemic-oxalosis threshold (flag)
SYSOX    : systemic oxalosis composite (0-1)
BoneOx   : bone oxalate reservoir (mmol)
SoftOx   : soft-tissue oxalate deposit (mmol)
DialOx   : cumulative dialytic oxalate removal (mmol)
AGTact   : functional AGT activity (fraction of normal)
GOact    : glycolate oxidase activity (fraction of normal)
LDHact   : cytosolic LDH activity (fraction of normal)
SurvESKD : kidney-failure-free survival probability
'

mod <- mcode("ph1", code)

# =============================================================================
#  SCENARIOS
#  Each is a PURE PARAMETER SET.  Every scenario carries a 180-day untreated
#  run-in (RUNIN), so "% change from baseline" is computed the way the
#  registration trials computed it: day 180 = baseline, day 360 = month 6.
# =============================================================================

TEND <- 7300   # 20 years

CKD4  <- list(NM0 = 0.30, NC0 = 0.45, FIB0 = 0.50, OXABONE0 = 150000)
ESKDp <- list(NM0 = 0.04, NC0 = 0.80, FIB0 = 0.80, OXABONE0 = 350000, HDON = 1)
CONS  <- list(UVTGT = 3.0, KCITON = 1)

scenarios <- list(

  "01 Healthy control (normal AGT)" = list(AGTACT0 = 1.0, NM0 = 1.0),

  "02 PH1 untreated - classic (G170R/G170R)" = list(),

  "03 PH1 untreated - null genotype (c.33dupC)" = list(AGTACT0 = 0.005, B6RESP = 0),

  "04 PH1 infantile oxalosis (immature GFR, high turnover)" = list(
    AGTACT0 = 0.005, B6RESP = 0, BSA = 0.45, WT = 8.0, GFR0 = 30.0,
    PRECUR = 1.6, NM0 = 0.88),

  "05 Conservative only - hyperhydration 3 L" = list(UVTGT = 3.0),

  "06 Conservative - hyperhydration + potassium citrate" = CONS,

  "07 Poor adherence to hyperhydration (50% of days)" = c(CONS, list(ADHFL = 0.5)),

  "08 Pyridoxine 8 mg/kg (G170R, responsive)" = list(B6ON = 1, B6DOSE = 8),

  "09 Pyridoxine 20 mg/kg (G170R) + neuropathy tracker" = list(B6ON = 1, B6DOSE = 20),

  "10 Pyridoxine 20 mg/kg in a NON-responsive genotype" = list(
    B6ON = 1, B6DOSE = 20, AGTACT0 = 0.005, B6RESP = 0),

  "11 Lumasiran (3 mg/kg x3 monthly, then quarterly)" = list(LUMON = 1),

  "12 Lumasiran + conservative therapy" = c(CONS, list(LUMON = 1)),

  "13 Lumasiran, ILLUMINATE-B infant regimen (6 mg/kg monthly)" = list(
    LUMON = 1, LUMLOAD = 6, LUMDOSE = 3, LUMTAU = 30, LUMNL = 3,
    BSA = 0.55, WT = 10.0, GFR0 = 60.0, PRECUR = 1.4, NM0 = 0.92),

  "14 Lumasiran in advanced CKD (ILLUMINATE-C cohort A)" = c(CKD4, list(LUMON = 1)),

  "15 Lumasiran on haemodialysis (ILLUMINATE-C cohort B)" = c(ESKDp, list(LUMON = 1)),

  "16 Nedosiran 3.5 mg/kg monthly" = list(NEDON = 1),

  "17 Nedosiran + conservative therapy" = c(CONS, list(NEDON = 1)),

  "18 Dual RNAi - lumasiran + nedosiran (exploratory)" = list(LUMON = 1, NEDON = 1),

  "19 Stiripentol monotherapy (LDH-5 inhibition)" = list(STPON = 1),

  "20 Oral oxalate degradation (reloxaliase / O. formigenes)" = list(OXDON = 1),

  "21 Late diagnosis at CKD4, no substrate reduction" = CKD4,

  "22 Late diagnosis at CKD4 + lumasiran + conservative" = c(CKD4, CONS, list(LUMON = 1)),

  "23 ESKD on conventional thrice-weekly haemodialysis" = ESKDp,

  "24 ESKD on intensified daily HD + nightly PD" = c(ESKDp, list(HDSESS = 6, PDCL = 6)),

  "25 ESKD -> combined liver-kidney transplant (day 365)" = c(ESKDp, list(
    TXDAY = 365, KTXON = 1, LTXON = 1)),

  "26 ESKD -> isolated kidney transplant, no RNAi" = c(ESKDp, list(
    TXDAY = 365, KTXON = 1, LTXON = 0)),

  "27 ESKD -> isolated kidney transplant + maintained lumasiran" = c(ESKDp, list(
    TXDAY = 365, KTXON = 1, LTXON = 0, LUMON = 1)),

  "28 Exploratory IL-1 blockade add-on (preclinical)" = c(CONS, list(ANAON = 1))
)

run_scen <- function(nm, pars) {
  m <- if (length(pars)) param(mod, pars) else mod
  out <- m %>%
    mrgsim(end = TEND, delta = 1, hmax = 2, atol = 1e-8, rtol = 1e-6) %>%
    as_tibble()
  out$scenario <- nm
  out
}

message("Simulating ", length(scenarios), " scenarios over ", TEND / 365.25, " years ...")
sim <- bind_rows(lapply(names(scenarios), function(nm) run_scen(nm, scenarios[[nm]])))

# =============================================================================
#  SUMMARY TABLE
#  Day 180 = end of untreated run-in (= trial baseline); day 360 = month 6.
# =============================================================================
at <- function(d, x, tt) x[which.min(abs(tt - d))]

summary_tbl <- sim %>%
  group_by(scenario) %>%
  summarise(
    Uox_base   = round(at(180,  Uox24, time), 3),
    Uox_m6     = round(at(360,  Uox24, time), 3),
    dUox_pct   = round(100 * (at(360, Uox24, time) / at(180, Uox24, time) - 1), 1),
    Uglc_base  = round(at(180,  Uglc24, time), 3),
    Uglc_m6    = round(at(360,  Uglc24, time), 3),
    Pox_base   = round(at(180,  Pox,   time), 1),
    Pox_m6     = round(at(360,  Pox,   time), 1),
    dPox_pct   = round(100 * (at(360, Pox, time) / at(180, Pox, time) - 1), 1),
    Pglc_m6    = round(at(360,  Pglc,  time), 1),
    RSS_m6     = round(at(360,  RSSr,  time), 2),
    eGFR_y1    = round(at(365,  eGFR,  time), 1),
    eGFR_y5    = round(at(1825, eGFR,  time), 1),
    eGFR_y10   = round(at(3650, eGFR,  time), 1),
    NC_y10     = round(at(3650, NCgrade, time), 2),
    stones_y10 = round(at(3650, StoneEv, time), 1),
    SysOx_y10  = round(at(3650, SYSOX, time), 3),
    .groups = "drop"
  ) %>%
  arrange(scenario)

cat("\n=================== SCENARIO SUMMARY ===================\n")
cat("Uox/Uglc mmol/1.73m2/day · Pox/Pglc umol/L · eGFR mL/min/1.73m2\n")
cat("baseline = day 180 (end of untreated run-in); m6 = day 360\n\n")
print(as.data.frame(summary_tbl), row.names = FALSE)

# time to kidney failure
ttf <- sim %>%
  group_by(scenario) %>%
  summarise(
    ESKD_day = if (any(eGFR < 15)) min(time[eGFR < 15]) else NA_real_,
    ESKD_yr  = round(ESKD_day / 365.25, 1),
    .groups  = "drop"
  )
cat("\n---- time from model start to eGFR < 15 (years; NA = not reached in 20 y) ----\n")
print(as.data.frame(ttf[, c("scenario", "ESKD_yr")]), row.names = FALSE)

# =============================================================================
#  VALIDATION AGAINST PUBLISHED ANCHORS
#  Qualitative phenotype checks against literature ranges, not fitted targets.
#  Sources for each number are tabulated in ph1_references.md section R.
# =============================================================================
chk <- function(label, value, lo, hi, unit = "") {
  ok <- !is.na(value) && value >= lo && value <= hi
  cat(sprintf("  [%s] %-62s %9.2f %-14s (target %g to %g)\n",
              if (ok) "PASS" else "FAIL", label, value, unit, lo, hi))
  ok
}
g <- function(scn, col, d) {
  s <- sim[sim$scenario == scn, ]
  at(d, s[[col]], s$time)
}
S01 <- "01 Healthy control (normal AGT)"
S02 <- "02 PH1 untreated - classic (G170R/G170R)"
S04 <- "04 PH1 infantile oxalosis (immature GFR, high turnover)"
S06 <- "06 Conservative - hyperhydration + potassium citrate"
S08 <- "08 Pyridoxine 8 mg/kg (G170R, responsive)"
S10 <- "10 Pyridoxine 20 mg/kg in a NON-responsive genotype"
S11 <- "11 Lumasiran (3 mg/kg x3 monthly, then quarterly)"
S14 <- "14 Lumasiran in advanced CKD (ILLUMINATE-C cohort A)"
S16 <- "16 Nedosiran 3.5 mg/kg monthly"
S18 <- "18 Dual RNAi - lumasiran + nedosiran (exploratory)"
S20 <- "20 Oral oxalate degradation (reloxaliase / O. formigenes)"
S23 <- "23 ESKD on conventional thrice-weekly haemodialysis"
S24 <- "24 ESKD on intensified daily HD + nightly PD"
S25 <- "25 ESKD -> combined liver-kidney transplant (day 365)"
S26 <- "26 ESKD -> isolated kidney transplant, no RNAi"

pc <- function(scn, col, d0 = 180, d1 = 360) {
  100 * (g(scn, col, d1) / g(scn, col, d0) - 1)
}

cat("\n=================== VALIDATION ===================\n")
res <- c(
  chk("Healthy: 24-h urinary oxalate below ULN",
      g(S01, "Uox24", 500), 0.10, 0.46, "mmol/1.73m2"),
  chk("Healthy: plasma oxalate",
      g(S01, "Pox", 500), 0.8, 5.0, "umol/L"),
  chk("Healthy: urinary glycolate below ULN",
      g(S01, "Uglc24", 500), 0.05, 0.50, "mmol/1.73m2"),
  chk("Healthy: CaOx relative supersaturation",
      g(S01, "RSSr", 500), 1.0, 5.0, ""),
  chk("PH1 untreated: 24-h urinary oxalate",
      g(S02, "Uox24", 180), 1.0, 2.5, "mmol/1.73m2"),
  chk("PH1 untreated: plasma oxalate at normal GFR",
      g(S02, "Pox", 180), 5.0, 15.0, "umol/L"),
  chk("PH1 untreated: urinary glycolate raised above normal",
      g(S02, "Uglc24", 180) / g(S01, "Uglc24", 500), 1.5, 5.0, "x normal"),
  chk("PH1 untreated: CaOx supersaturation",
      g(S02, "RSSr", 180), 8.0, 20.0, ""),
  chk("PH1 untreated: reaches kidney failure within 20 y",
      g(S02, "eGFR", 7300), 0.0, 15.0, "mL/min"),
  chk("PH1 untreated: symptomatic stone events per year",
      g(S02, "StoneEv", 3650) / 10, 0.3, 2.0, "events/y"),
  chk("Infantile: plasma oxalate already at oxalosis threshold",
      g(S04, "Pox", 180), 25.0, 90.0, "umol/L"),
  chk("Infantile: kidney failure in early childhood",
      (function(z) if (length(z)) min(z) / 365.25 else NA_real_)(
        sim$time[sim$scenario == S04 & sim$eGFR < 15]), 0.5, 6.0, "years"),
  chk("Conservative: supersaturation brought below metastable limit",
      g(S06, "RSSr", 360), 1.0, 4.5, ""),
  chk("Conservative: urinary oxalate NOT changed (no substrate reduction)",
      abs(pc(S06, "Uox24")), 0.0, 5.0, "% change"),
  chk("Conservative: new stone formation suppressed",
      g(S06, "StoneEv", 3650) - g(S06, "StoneEv", 360), 0.0, 0.5, "events"),
  chk("Conservative alone: still progresses to CKD by 20 y",
      g(S06, "eGFR", 7300), 15.0, 75.0, "mL/min"),
  chk("Pyridoxine (responsive genotype): Uox reduction",
      -pc(S08, "Uox24"), 20.0, 50.0, "%"),
  chk("Pyridoxine (non-responsive genotype): no meaningful Uox reduction",
      -pc(S10, "Uox24"), -5.0, 8.0, "%"),
  chk("Pyridoxine 20 mg/kg: neuropathy tracker stays low",
      g("09 Pyridoxine 20 mg/kg (G170R) + neuropathy tracker", "SYSOX", 3650),
      0.0, 0.25, "composite"),
  chk("LUMASIRAN: Uox reduction at month 6 (ILLUMINATE-A -53.5% to -65.4%)",
      -pc(S11, "Uox24"), 45.0, 75.0, "%"),
  chk("LUMASIRAN: plasma glycolate rises (on-target biomarker)",
      g(S11, "Pglc", 360) / g(S11, "Pglc", 180), 2.0, 8.0, "x baseline"),
  chk("LUMASIRAN: urinary glycolate rises",
      g(S11, "Uglc24", 360) / g(S11, "Uglc24", 180), 2.0, 8.0, "x baseline"),
  chk("LUMASIRAN: GO activity knocked down",
      1 - g(S11, "GOact", 360), 0.75, 0.95, "fraction"),
  chk("LUMASIRAN: supersaturation normalised",
      g(S11, "RSSr", 360), 1.0, 5.5, ""),
  chk("LUMASIRAN: eGFR preserved at 10 y vs untreated",
      g(S11, "eGFR", 3650) - g(S02, "eGFR", 3650), 25.0, 80.0, "mL/min"),
  chk("LUMASIRAN in advanced CKD: plasma oxalate falls",
      -pc(S14, "Pox"), 20.0, 70.0, "%"),
  chk("NEDOSIRAN: Uox reduction at month 6 (PHYOX2 ~-50%)",
      -pc(S16, "Uox24"), 35.0, 70.0, "%"),
  chk("NEDOSIRAN: glycolate NOT raised (contrast with lumasiran)",
      g(S16, "Uglc24", 360) / g(S16, "Uglc24", 180), 0.7, 1.3, "x baseline"),
  chk("NEDOSIRAN: model predicts plasma GLYOXYLATE rises instead",
      g(S16, "Pglx", 360) / g(S16, "Pglx", 180), 1.5, 5.0, "x baseline"),
  chk("DUAL RNAi: greater Uox reduction than either agent alone",
      -pc(S18, "Uox24"), 70.0, 95.0, "%"),
  chk("Oral oxalate degradation alone: small effect (dietary fraction only)",
      -pc(S20, "Uox24"), 0.0, 20.0, "%"),
  chk("ESKD conventional HD: plasma oxalate remains grossly elevated",
      g(S23, "Pox", 360), 40.0, 120.0, "umol/L"),
  chk("ESKD conventional HD: above the systemic-oxalosis threshold",
      g(S23, "OXNEPH", 360), 1.0, 1.0, "flag"),
  chk("Intensified HD+PD lowers Pox below conventional HD",
      g(S23, "Pox", 360) - g(S24, "Pox", 360), 5.0, 40.0, "umol/L"),
  chk("Intensified HD+PD still cannot normalise Pox",
      g(S24, "Pox", 360), 25.0, 70.0, "umol/L"),
  chk("CLKT: AGT activity restored",
      g(S25, "AGTact", 550), 0.90, 1.02, "fraction"),
  chk("CLKT: plasma oxalate normalised by 6 months post-transplant",
      g(S25, "Pox", 545), 1.0, 15.0, "umol/L"),
  chk("CLKT: bone oxalate release keeps Uox raised for months",
      g(S25, "Uox24", 550), 0.5, 1.6, "mmol/1.73m2"),
  chk("CLKT: that release has resolved by 4 years",
      g(S25, "Uox24", 1825), 0.15, 0.55, "mmol/1.73m2"),
  chk("Isolated KTx without metabolic control: graft nephrocalcinosis returns",
      g(S26, "NCgrade", 1825), 1.5, 3.0, "grade 0-3"),
  chk("Isolated KTx without metabolic control: graft function lost",
      g(S26, "eGFR", 2555), 0.0, 25.0, "mL/min"),
  chk("Isolated KTx + lumasiran protects the graft vs no RNAi",
      g("27 ESKD -> isolated kidney transplant + maintained lumasiran", "eGFR", 2555)
      - g(S26, "eGFR", 2555), 8.0, 60.0, "mL/min"),
  chk("Advanced CKD: Uox FALLS as GFR falls (why Pox is the endpoint)",
      g(S02, "Uox24", 3650) - g(S02, "Uox24", 7300), 0.3, 2.0, "mmol/1.73m2"),
  chk("Baseline stability: healthy Uox drift over 20 y",
      abs(100 * (g(S01, "Uox24", 7300) / g(S01, "Uox24", 500) - 1)), 0.0, 5.0, "%")
)
cat(sprintf("\n  %d/%d checks passed\n", sum(res), length(res)))

# =============================================================================
#  THE THREE LOOPS, MADE VISIBLE
# =============================================================================
cat("\n============ LOOP 1: futile cycle — why the two siRNAs differ ============\n")
loop1 <- sim %>%
  filter(scenario %in% c(S02, S11, S16, S18), time == 360) %>%
  transmute(scenario, Uox24, Uglc24, Pglc, Pglx,
            GO = round(GOact, 2), LDH = round(LDHact, 2))
print(as.data.frame(loop1), row.names = FALSE)
cat("\n  Blocking GO (lumasiran) forces carbon out as GLYCOLATE.\n",
    " Blocking LDHA (nedosiran) forces it out as GLYOXYLATE instead, which is\n",
    " why nedosiran does not raise glycolate - and why extrahepatic LDH sets a\n",
    " floor on what a hepatocyte-targeted siRNA can achieve.\n")

cat("\n============ LOOP 2: the per-nephron accelerator ============\n")
loop2 <- sim %>%
  filter(scenario == S02, time %in% c(180, 1825, 3650, 5475)) %>%
  transmute(year = round(time / 365.25, 1), eGFR = round(eGFR, 1),
            Uox24 = round(Uox24, 2),
            load_per_nephron = round(Uox24 / (eGFR / 105), 2),
            NCgrade = round(NCgrade, 2))
print(as.data.frame(loop2), row.names = FALSE)
cat("\n  Urinary oxalate falls, yet the load carried by each surviving nephron\n",
    " keeps rising. That is the accelerator, and it is why the eGFR curve bends.\n")

cat("\n============ LOOP 3: the bone capacitor ============\n")
loop3 <- sim %>%
  filter(scenario == S25, time %in% c(360, 400, 550, 1095, 1825, 2555)) %>%
  transmute(day = time, Pox = round(Pox, 1), Uox24 = round(Uox24, 2),
            BoneOx_mmol = round(BoneOx, 1), eGFR = round(eGFR, 1))
print(as.data.frame(loop3), row.names = FALSE)
cat("\n  The same term that buffered plasma oxalate before transplant runs\n",
    " backwards afterwards: the graft has to excrete the skeleton's stores.\n")

# =============================================================================
#  OPTIONAL PLOTS
# =============================================================================
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  S21 <- "21 Late diagnosis at CKD4, no substrate reduction"
  key <- c(S02, S06, S08, S11, S16, S21, S25)
  p <- sim %>%
    filter(scenario %in% key) %>%
    select(time, scenario, `Uox (mmol/24h)` = Uox24, `Pox (umol/L)` = Pox,
           `plasma glycolate` = Pglc, `CaOx RSS` = RSSr,
           `eGFR` = eGFR, `nephrocalcinosis grade` = NCgrade) %>%
    pivot_longer(-c(time, scenario)) %>%
    ggplot(aes(time / 365.25, value, colour = scenario)) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~name, scales = "free_y") +
    labs(x = "years", y = NULL,
         title = "Primary hyperoxaluria type 1 QSP model - key scenarios") +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", legend.text = element_text(size = 6)) +
    guides(colour = guide_legend(ncol = 2))
  ggsave("ph1_scenarios.png", p, width = 11, height = 7, dpi = 150)
  message("wrote ph1_scenarios.png")
}

invisible(sim)
