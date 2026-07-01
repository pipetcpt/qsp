## =============================================================================
## Alport Syndrome (AS) — mrgsolve QSP model
## COL4A3/A4/A5 genotype-dependent type-IV collagen network failure ->
##   glomerular basement membrane (GBM) mechanical instability -> podocyte
##   foot-process effacement/apoptosis -> proteinuria/hematuria -> compensatory
##   hyperfiltration & RAAS/endothelin-driven glomerular hypertension ->
##   TGF-beta1/CTGF/miR-21 fibrotic cascade -> nephron loss -> eGFR decline ->
##   ESRD, coupled to cochlear (progressive SNHL) and ocular (lenticonus)
##   basement-membrane phenotypes, and to RAAS blockade (ACEi/ARB), sparsentan
##   (dual ETA/AT1 antagonist), bardoxolone methyl (Nrf2 activator), lademirsen
##   (anti-miR-21 ASO), dapagliflozin (SGLT2i), and finerenone (nonsteroidal
##   MRA) PK/PD.
##
## Time unit: hours. Disease horizon simulated over years (end = 24*365*N).
##
## Calibration anchors (see alp_references.md for full PMID list):
##   - Natural history eGFR decline & genotype-specific ESRD age: Jais 2000/2003
##     JASN (XLAS male ESRD median ~25yr, female ~60s); Gross 2012 KI (ARAS ESRD
##     childhood-adolescence); Kashtan 2018 Pediatr Nephrol clinical guidance
##   - ACEi/ARB delay-of-ESRD natural history comparator: Gross 2012 Kidney Int
##     (ramipril cohort, ~13-yr ESRD delay vs untreated historical controls);
##     EARLY PRO-TECT Alport (Gross 2020 Kidney Int, presymptomatic ramipril
##     initiation, albuminuria/GBM benefit)
##   - Sparsentan dual ETA/AT1 antagonism: DUPLEX FSGS trial PK/PD bridging
##     (Komers 2022 Kidney Int Rep) + dedicated pediatric/adult Alport program
##     rationale (Gross 2022 Nephrol Dial Transplant perspective)
##   - Bardoxolone methyl Nrf2 activation, Alport-specific: CARDINAL trial
##     (Chertow 2021 CJASN; eGFR acute rise, chronic-phase slope benefit signal,
##     cardiovascular fluid-retention caveat noted in label)
##   - Lademirsen (RG-012) anti-miR-21 ASO: HERA trial (Gomez 2022 topline;
##     phase 2 discontinued for lack of efficacy signal at interim -> modeled
##     with attenuated EMAX reflecting negative/inconclusive trial)
##   - SGLT2i tubuloglomerular feedback / hyperfiltration reduction: extrapolated
##     from DAPA-CKD (Heerspink 2020 NEJM) proteinuric CKD subgroup physiology
##   - Finerenone MRA anti-fibrotic: FIDELIO-DKD/FIGARO-DKD (Bakris 2020/Pitt
##     2021 NEJM) mechanistic extrapolation to non-diabetic proteinuric CKD
##   - Cochlear SNHL natural history: Merchant 2004 Otol Neurotol (temporal
##     bone histopathology); Gubbels 2009 (audiometric progression cohort)
##   - Ocular lenticonus/retinopathy: Colville 1997 Ophthalmology (prevalence
##     by genotype/sex); Kato 2008 Am J Ophthalmol
## =============================================================================

$PROB
# Alport Syndrome QSP model (24-compartment PK/PD + disease + safety)

$PARAM
// ---- Genotype / severity ----
SEVERITY     = 1.0     // 0.4=ADAS(het thin-BM), 1.0=XLAS-male/ARAS(reference), 0.55=XLAS-female(mosaic), 1.15=ARAS-truncating
IS_MALE_XL   = 1       // 1 = XLAS hemizygous male (fastest progression), 0 = female/autosomal

// ---- Ramipril (ACEi) oral PK ----
RAM_KA   = 0.90    // 1/h
RAM_F    = 0.60
RAM_V1   = 60       // L (ramiprilat)
RAM_CL   = 8.0      // L/h
ACEI_EC50 = 4.0      // relative conc for half-max AngII suppression
ACEI_EMAX = 0.55     // max fractional AngII generation suppression

// ---- Losartan (ARB) oral PK (E-3174 active metabolite) ----
LOS_KA   = 0.70
LOS_F    = 0.33      // fraction converted to active carboxylic-acid metabolite
LOS_V1   = 50
LOS_CL   = 6.5
ARB_EC50 = 3.5
ARB_EMAX = 0.55      // max fractional AT1 blockade

// ---- Sparsentan (dual ETA/AT1 antagonist) oral PK ----
SPAR_KA  = 0.55
SPAR_F   = 0.70
SPAR_V1  = 40
SPAR_CL  = 2.5
SPAR_EC50 = 2.0
SPAR_EMAX_AT1 = 0.60   // AT1 component
SPAR_EMAX_ETA = 0.50   // ETA component
SPAR_EMAX_FIB = 0.35   // direct additive anti-fibrotic (beyond hemodynamic)

// ---- Bardoxolone methyl (Nrf2 activator) oral PK ----
BARD_KA  = 0.30
BARD_F   = 0.45
BARD_V1  = 70
BARD_CL  = 3.2
BARD_EC50 = 1.5
BARD_EMAX_ROS  = 0.45   // antioxidant/Nrf2 gene induction
BARD_EMAX_FIB  = 0.20   // modest chronic anti-fibrotic
BARD_ACUTE_EGFR_BUMP = 4.0  // mL/min/1.73m2, creatinine-independent tubular-secretion pseudo-rise
BARD_CV_RISK_SIGNAL = 1     // flag: fluid-retention/CV safety caveat (CARDINAL, informational only)

// ---- Lademirsen (RG-012, anti-miR-21 ASO) SC PK ----
LAD_KA   = 0.05     // slow SC absorption
LAD_F    = 0.50
LAD_V1   = 8
LAD_CL   = 0.35
LAD_KTISSUE = 0.02   // 1/h, plasma -> kidney-tissue distribution
LAD_KOUT_TIS = 0.010 // 1/h, tissue clearance
LAD_EC50 = 1.0
LAD_EMAX_MIR21 = 0.40   // attenuated per HERA negative/inconclusive phase-2 signal

// ---- Dapagliflozin (SGLT2i) oral PK ----
DAPA_KA  = 1.1
DAPA_F   = 0.78
DAPA_V1  = 118
DAPA_CL  = 8.5
SGLT2_EC50 = 0.06
SGLT2_EMAX_SNGFR = 0.30   // fractional single-nephron hyperfiltration reduction via TGF restoration

// ---- Finerenone (nonsteroidal MRA) oral PK ----
FIN_KA   = 0.60
FIN_F    = 0.44
FIN_V1   = 45
FIN_CL   = 6.0
FIN_EC50 = 0.30
FIN_EMAX_ALDO_FIB = 0.35   // fractional aldosterone-driven fibrosis suppression

// ---- Disease progression rate constants ----
K_GBM0      = 0.00030   // 1/h baseline GBM structural failure rate
K_NEPHRON0  = 0.00012   // 1/h baseline nephron-loss rate (sclerosis accrual)
K_PODO0     = 0.00020   // 1/h baseline podocyte depletion rate
HYPERFILT_GAIN = 1.8    // amplification of SNGFR from nephron loss
RAAS_GAIN   = 1.2       // AngII generation gain per unit glomerular pressure
ET1_GAIN    = 0.9       // ET-1 upregulation gain per unit mesangial activation
FIB_GAIN    = 1.4        // TGF-b1/CTGF drive per unit AT1+ETA signaling
MIR21_GAIN  = 0.8        // miR-21 amplification of CTGF
VICIOUS_GAIN = 0.55       // feedback gain, glomerular pressure -> further GBM/podocyte stress
EGFR_MAX    = 100         // mL/min/1.73m2 healthy reference
EGFR_HALFLIFE_K = 0.60     // relaxation rate of eGFR toward mechanistic target
UACR_MAX    = 3000         // mg/g ceiling
COCHLEA_K   = 0.00006      // 1/h baseline cochlear BM/hair-cell decline
OCULAR_K    = 0.00004      // 1/h baseline lens/retinal collagen decline
ESRD_EGFR_THRESHOLD = 15    // mL/min/1.73m2

$CMT @annotated
RAM_GUT     : Ramipril gut depot (mg)
RAM_CENT    : Ramiprilat plasma concentration compartment (mg)
LOS_GUT     : Losartan gut depot (mg)
LOS_CENT    : Losartan carboxylic-acid (E-3174) plasma compartment (mg)
SPAR_GUT    : Sparsentan gut depot (mg)
SPAR_CENT   : Sparsentan plasma compartment (mg)
BARD_GUT    : Bardoxolone methyl gut depot (mg)
BARD_CENT   : Bardoxolone methyl plasma compartment (mg)
LAD_SC      : Lademirsen SC depot (mg)
LAD_PLASMA  : Lademirsen plasma compartment (mg)
LAD_TISSUE  : Lademirsen kidney-tissue compartment (relative units)
DAPA_GUT    : Dapagliflozin gut depot (mg)
DAPA_CENT   : Dapagliflozin plasma compartment (mg)
FIN_GUT     : Finerenone gut depot (mg)
FIN_CENT    : Finerenone plasma compartment (mg)
GBM_INTEG   : GBM structural integrity index (1=normal, 0=failed)
PODO_FRAC   : Viable podocyte fraction (0-1)
NEPHRON_FRAC: Surviving functional nephron fraction (0-1)
FIBROSIS    : Interstitial fibrosis / IFTA index (0-1)
MIR21       : Relative miR-21 activity index
UACR        : Urine albumin-to-creatinine ratio (mg/g)
EGFR        : Estimated GFR, CKD-EPI (mL/min/1.73m2)
HEARING_LOSS: Hearing threshold shift composite (dB HL)
OCULAR_SCORE: Ocular (lenticonus/retinopathy) severity index (0-1)
FX_YEARS    : Elapsed simulation time tracker (years)

$MAIN
double RAM_Cp  = RAM_CENT/RAM_V1;
double LOS_Cp  = LOS_CENT/LOS_V1;
double SPAR_Cp = SPAR_CENT/SPAR_V1;
double BARD_Cp = BARD_CENT/BARD_V1;
double LAD_Cp  = LAD_TISSUE;
double DAPA_Cp = DAPA_CENT/DAPA_V1;
double FIN_Cp  = FIN_CENT/FIN_V1;

// ---- RAAS blockade (ACEi / ARB) ----
double acei_effect = ACEI_EMAX*(RAM_Cp/(ACEI_EC50+RAM_Cp+1e-9));
double arb_effect  = ARB_EMAX*(LOS_Cp/(ARB_EC50+LOS_Cp+1e-9));
double spar_at1    = SPAR_EMAX_AT1*(SPAR_Cp/(SPAR_EC50+SPAR_Cp+1e-9));
double spar_eta    = SPAR_EMAX_ETA*(SPAR_Cp/(SPAR_EC50+SPAR_Cp+1e-9));
double spar_fib    = SPAR_EMAX_FIB*(SPAR_Cp/(SPAR_EC50+SPAR_Cp+1e-9));
double at1_block_total = 1.0 - (1.0-acei_effect)*(1.0-arb_effect)*(1.0-spar_at1);
if (at1_block_total > 0.95) at1_block_total = 0.95;
double eta_block_total = spar_eta;

// ---- Bardoxolone Nrf2 effect ----
double bard_ros_effect = BARD_EMAX_ROS*(BARD_Cp/(BARD_EC50+BARD_Cp+1e-9));
double bard_fib_effect = BARD_EMAX_FIB*(BARD_Cp/(BARD_EC50+BARD_Cp+1e-9));
double bard_acute_bump = BARD_ACUTE_EGFR_BUMP*(BARD_Cp/(BARD_EC50+BARD_Cp+1e-9));

// ---- Lademirsen anti-miR-21 ----
double lad_mir21_effect = LAD_EMAX_MIR21*(LAD_Cp/(LAD_EC50+LAD_Cp+1e-9));

// ---- SGLT2i tubuloglomerular feedback ----
double sglt2_sngfr_effect = SGLT2_EMAX_SNGFR*(DAPA_Cp/(SGLT2_EC50+DAPA_Cp+1e-9));

// ---- Finerenone MRA anti-fibrotic ----
double fin_fib_effect = FIN_EMAX_ALDO_FIB*(FIN_Cp/(FIN_EC50+FIN_Cp+1e-9));

// ---- Composite hemodynamic & fibrotic drive ----
double sngfr = HYPERFILT_GAIN*(1.0-NEPHRON_FRAC+1e-6)*(1.0-sglt2_sngfr_effect);
double glomerular_pressure = RAAS_GAIN*sngfr*(1.0-at1_block_total);
double et1_activity = ET1_GAIN*(1.0-NEPHRON_FRAC)*(1.0-eta_block_total);
double fibrotic_drive = FIB_GAIN*(glomerular_pressure+et1_activity)*(1.0+MIR21_GAIN*MIR21)
                          *(1.0-spar_fib)*(1.0-bard_fib_effect)*(1.0-fin_fib_effect);
if (fibrotic_drive < 0) fibrotic_drive = 0;
double vicious_feedback = VICIOUS_GAIN*glomerular_pressure;

$ODE
// ---------------- Ramipril / ramiprilat PK ----------------
dxdt_RAM_GUT  = -RAM_KA*RAM_GUT;
dxdt_RAM_CENT =  RAM_KA*RAM_GUT*RAM_F - (RAM_CL/RAM_V1)*RAM_CENT;

// ---------------- Losartan / E-3174 PK ----------------
dxdt_LOS_GUT  = -LOS_KA*LOS_GUT;
dxdt_LOS_CENT =  LOS_KA*LOS_GUT*LOS_F - (LOS_CL/LOS_V1)*LOS_CENT;

// ---------------- Sparsentan PK ----------------
dxdt_SPAR_GUT  = -SPAR_KA*SPAR_GUT;
dxdt_SPAR_CENT =  SPAR_KA*SPAR_GUT*SPAR_F - (SPAR_CL/SPAR_V1)*SPAR_CENT;

// ---------------- Bardoxolone methyl PK ----------------
dxdt_BARD_GUT  = -BARD_KA*BARD_GUT;
dxdt_BARD_CENT =  BARD_KA*BARD_GUT*BARD_F - (BARD_CL/BARD_V1)*BARD_CENT;

// ---------------- Lademirsen PK (SC -> plasma -> kidney tissue) ----------------
dxdt_LAD_SC     = -LAD_KA*LAD_SC;
dxdt_LAD_PLASMA =  LAD_KA*LAD_SC*LAD_F - (LAD_CL/LAD_V1)*LAD_PLASMA - LAD_KTISSUE*LAD_PLASMA;
dxdt_LAD_TISSUE =  LAD_KTISSUE*LAD_PLASMA - LAD_KOUT_TIS*LAD_TISSUE;

// ---------------- Dapagliflozin PK ----------------
dxdt_DAPA_GUT  = -DAPA_KA*DAPA_GUT;
dxdt_DAPA_CENT =  DAPA_KA*DAPA_GUT*DAPA_F - (DAPA_CL/DAPA_V1)*DAPA_CENT;

// ---------------- Finerenone PK ----------------
dxdt_FIN_GUT  = -FIN_KA*FIN_GUT;
dxdt_FIN_CENT =  FIN_KA*FIN_GUT*FIN_F - (FIN_CL/FIN_V1)*FIN_CENT;

// ---------------- GBM structural integrity ----------------
double k_gbm = K_GBM0*SEVERITY*(1.0+vicious_feedback)*(1.0-bard_ros_effect*0.3);
if (k_gbm < 0) k_gbm = 0;
dxdt_GBM_INTEG = -k_gbm*GBM_INTEG;

// ---------------- Podocyte depletion ----------------
double k_podo = K_PODO0*SEVERITY*(1.0-GBM_INTEG+1.0)*(1.0+vicious_feedback)*(1.0-bard_ros_effect*0.2);
if (k_podo < 0) k_podo = 0;
dxdt_PODO_FRAC = -k_podo*PODO_FRAC;

// ---------------- miR-21 pro-fibrotic amplifier ----------------
dxdt_MIR21 = 0.02*(fibrotic_drive) - 0.02*MIR21*(1.0+lad_mir21_effect) ;

// ---------------- Fibrosis / IFTA accrual ----------------
dxdt_FIBROSIS = 0.00015*fibrotic_drive*(1.0-FIBROSIS) - 0.00002*FIBROSIS;
if (FIBROSIS > 1.0) FIBROSIS = 1.0;

// ---------------- Nephron loss ----------------
double k_nephron = K_NEPHRON0*SEVERITY*(1.0+2.0*FIBROSIS)*(1.0+0.5*(1.0-PODO_FRAC))*(1.0-at1_block_total*0.4);
if (k_nephron < 0) k_nephron = 0;
dxdt_NEPHRON_FRAC = -k_nephron*NEPHRON_FRAC;

// ---------------- Proteinuria (UACR) ----------------
double uacr_target = UACR_MAX*(1.0-PODO_FRAC)*(1.0-GBM_INTEG*0.3)*(1.0+0.5*glomerular_pressure);
if (uacr_target > UACR_MAX) uacr_target = UACR_MAX;
if (uacr_target < 8) uacr_target = 8;
dxdt_UACR = 0.01*(uacr_target-UACR);

// ---------------- eGFR ----------------
double egfr_target = EGFR_MAX*NEPHRON_FRAC*PODO_FRAC*(1.0-FIBROSIS*0.5) + bard_acute_bump*(BARD_CV_RISK_SIGNAL);
if (egfr_target < 2) egfr_target = 2;
if (egfr_target > EGFR_MAX) egfr_target = EGFR_MAX;
dxdt_EGFR = EGFR_HALFLIFE_K*0.002*(egfr_target-EGFR);

// ---------------- Cochlear hearing loss ----------------
double hearing_target = 90.0*(1.0-pow(GBM_INTEG,0.5));
dxdt_HEARING_LOSS = COCHLEA_K*SEVERITY*(hearing_target-HEARING_LOSS);

// ---------------- Ocular severity ----------------
double ocular_target = 1.0*(1.0-pow(GBM_INTEG,0.6));
dxdt_OCULAR_SCORE = OCULAR_K*SEVERITY*(ocular_target-OCULAR_SCORE);

dxdt_FX_YEARS = 1.0/(24.0*365.0);

$INIT
RAM_GUT = 0, RAM_CENT = 0, LOS_GUT = 0, LOS_CENT = 0, SPAR_GUT = 0, SPAR_CENT = 0,
BARD_GUT = 0, BARD_CENT = 0, LAD_SC = 0, LAD_PLASMA = 0, LAD_TISSUE = 0,
DAPA_GUT = 0, DAPA_CENT = 0, FIN_GUT = 0, FIN_CENT = 0,
GBM_INTEG = 1.0, PODO_FRAC = 1.0, NEPHRON_FRAC = 1.0, FIBROSIS = 0, MIR21 = 0,
UACR = 15, EGFR = 100, HEARING_LOSS = 0, OCULAR_SCORE = 0, FX_YEARS = 0

$CAPTURE acei_effect arb_effect spar_at1 spar_eta spar_fib bard_ros_effect bard_fib_effect
$CAPTURE lad_mir21_effect sglt2_sngfr_effect fin_fib_effect sngfr glomerular_pressure
$CAPTURE et1_activity fibrotic_drive at1_block_total eta_block_total

## =============================================================================
## Treatment scenarios (see alp_shiny_app.R for interactive dosing UI)
##
## 1. Natural history, XLAS male (SEVERITY=1.0, IS_MALE_XL=1, untreated)
## 2. Natural history, ARAS (SEVERITY=1.15, untreated)
## 3. Natural history, ADAS/heterozygous thin-BM (SEVERITY=0.4, untreated)
## 4. Ramipril (ACEi) chronic daily oral, started at proteinuria onset
##    (Gross 2012 KI cohort dosing ~5-10 mg/day titrated)
## 5. Presymptomatic early ramipril initiation (EARLY PRO-TECT Alport regimen,
##    started before overt proteinuria)
## 6. Losartan (ARB) chronic daily oral, alternative/add-on RAAS blockade
## 7. Sparsentan chronic daily oral (dual ETA/AT1 antagonism)
## 8. Bardoxolone methyl chronic daily oral (Nrf2 activation, CARDINAL regimen;
##    note acute creatinine-independent eGFR bump vs chronic slope effect)
## 9. Lademirsen (RG-012) chronic subcutaneous dosing (anti-miR-21 ASO, HERA
##    regimen; attenuated effect per inconclusive phase 2 data)
## 10. Combination: ACEi/ARB maximal RAAS blockade + dapagliflozin (SGLT2i) +
##     sparsentan (multi-target regimen)
##
## Example mrgsolve event code (ramipril chronic dosing, XLAS male):
##   mod <- mread("alp_mrgsolve_model") %>% param(SEVERITY = 1.0, IS_MALE_XL = 1)
##   e_ram <- ev(amt = 5, cmt = "RAM_GUT", time = 0, ii = 24, addl = 365*15-1)
##   out <- mod %>% ev(e_ram) %>% mrgsim(end = 24*365*15, delta = 24)  # 15-yr horizon
## =============================================================================
