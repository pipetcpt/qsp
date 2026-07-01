## =============================================================================
## C3 Glomerulopathy (C3G: Dense Deposit Disease [DDD] & C3 Glomerulonephritis
## [C3GN]) -- mrgsolve QSP model
##
## Genetic/autoantibody alternative-pathway (AP) dysregulation (CFH/CFI/CFB/C3
## variants, C3 Nephritic Factor [C3NeF]) -> uncontrolled C3 convertase
## (C3bBb) amplification -> glomerular C3b deposition (mesangial/subendothelial
## in C3GN, intramembranous dense-ribbon in DDD) -> mesangial proliferation,
## podocyte injury, proteinuria -> TGF-beta-driven interstitial fibrosis ->
## nephron loss -> eGFR decline -> ESKD/transplant recurrence, coupled to
## upstream alternative-pathway inhibitor (iptacopan [Factor B], pegcetacoplan
## [C3/C3b], danicopan [Factor D]) and terminal-pathway inhibitor (eculizumab,
## ravulizumab [C5]) PK/PD.
##
## Time unit: hours. Disease horizon simulated over years (end = 24*365*N).
##
## Calibration anchors (see c3g_references.md for full PMID list):
##   - AP amplification loop & C3NeF stabilization kinetics: Pickering 2013
##     Kidney Int; Zhang 2012 JASN
##   - KDIGO/consensus C3G histologic classification (DDD vs C3GN,
##     C3-dominant IF criterion, EM ribbon deposits): Pickering 2013 Kidney
##     Int; Sethi & Fervenza 2012 Semin Nephrol; Goodship 2017 Kidney Int
##     (C3G Consensus Conference)
##   - Natural-history progression to ESKD (~50% by 10yr): Servais 2012
##     Kidney Int; Iatropoulos 2018 Mol Immunol registry data
##   - Iptacopan Factor B inhibition, oral PK/PD & AP-blockade: APPEAR-C3G
##     phase 3 (Bomback 2025 NEJM/ASN abstract topline; proteinuria reduction
##     ~35% vs placebo); iptacopan PNH-program PK bridging (Jang 2022 Clin
##     Pharmacokinet)
##   - Pegcetacoplan C3/C3b inhibition: Wong 2023 Kidney Int Rep phase 2 (C3G
##     & IC-MPGN, UPCR reduction, C3 normalization); PEGASUS PNH-program PK
##     bridging (Hillmen 2021 NEJM, elimination half-life ~8 days)
##   - Eculizumab/ravulizumab C5 blockade in C3G, modest/heterogeneous benefit
##     with persistent upstream C3 deposition: Bomback 2012 CJASN (eculizumab
##     case series, sC5b-9 responders); Le Quintrec 2018 Am J Transplant;
##     ravulizumab PK from PNH program (Lee 2019 Blood, near-complete C5
##     saturation with FcRn-recycling extended half-life)
##   - Danicopan Factor D add-on for breakthrough AP activity atop C5
##     inhibitor background: PNH extravascular-hemolysis program (Risitano
##     2021 Blood; ALPHA trial) mechanistic extrapolation to C3G breakthrough
##     deposition
##   - Kidney transplant recurrence (up to 50-90% histologic recurrence,
##     allograft loss): Zand 2014 CJASN; Regunathan-Shenk 2019 AJKD
##   - Serum C3/C4/sC5b-9/Factor B biomarker natural history: Zhang 2012 CJASN;
##     Corvillo 2019 Front Immunol (AP functional assay)
## =============================================================================

$PROB
# C3 Glomerulopathy (C3G) QSP model (21-compartment PK/PD + disease + biomarker)

$PARAM
// ---- Genotype / phenotype severity ----
SEVERITY     = 1.0     // 0.7=milder C3GN, 1.0=reference C3GN, 1.35=DDD (intramembranous, faster progression)
CFH_LOF      = 0       // 1 = CFH loss-of-function variant present (raises baseline AP drive)
C3NEF_TITER  = 1.0      // relative C3NeF autoantibody titer (0=absent/seronegative, 1=reference positive, 2=high-titer persistent)
GENETIC_GAIN = 1.0      // composite gain from CFB/C3 gain-of-function variants (1=reference, >1=stronger convertase-stabilizing variant)

// ---- Iptacopan (oral Factor B inhibitor) PK ----
IPTA_KA   = 0.80    // 1/h
IPTA_F    = 0.85
IPTA_V1   = 25       // L
IPTA_CL   = 3.0      // L/h
IPTA_EC50 = 3.0      // relative plasma conc for half-max AP blockade
IPTA_EMAX = 0.92     // near-complete Factor B blockade at steady trough (APPEAR-C3G regimen)

// ---- Pegcetacoplan (SC C3/C3b inhibitor) PK ----
PEG_KA    = 0.02     // slow SC depot absorption
PEG_F     = 0.75
PEG_V1    = 6        // L (small Vd, large peptide)
PEG_CL    = 0.045    // L/h (long half-life ~8 days)
PEG_EC50  = 40
PEG_EMAX  = 0.95     // near-complete C3/AP+CP convertase blockade

// ---- Eculizumab (IV C5 inhibitor) PK (2-compartment mAb) ----
ECU_V1    = 3.5      // L central
ECU_VP    = 3.0      // L peripheral
ECU_Q     = 0.015    // L/h inter-compartmental
ECU_CL    = 0.0055   // L/h
ECU_EC50  = 15
ECU_EMAX  = 0.98     // near-complete C5/terminal blockade at trough

// ---- Ravulizumab (IV C5 inhibitor, engineered FcRn recycling) PK ----
RAVU_V1   = 3.2
RAVU_VP   = 2.8
RAVU_Q    = 0.010
RAVU_CL   = 0.0011   // ~4x lower CL than eculizumab -> q8w dosing
RAVU_EC50 = 12
RAVU_EMAX = 0.98

// ---- Danicopan (oral Factor D inhibitor, add-on) PK ----
DANI_KA   = 1.1
DANI_F    = 0.70
DANI_V1   = 18
DANI_CL   = 2.4
DANI_EC50 = 2.5
DANI_EMAX = 0.55     // partial add-on AP suppression on top of C5i background

// ---- Alternative-pathway amplification & convertase kinetics ----
AP_BASELINE   = 0.05     // baseline healthy AP tick-over
AP_GAIN       = 1.0      // scaling for disease amplification loop
AP_K_IN       = 0.08     // 1/h rate constant driving AP_ACTIVITY toward its target
C3_SYNTH      = 1.2      // mg/dL/h hepatic C3 synthesis rate
C3_NORMAL     = 110      // mg/dL reference normal serum C3
C3_CONSUME_K  = 0.014    // consumption rate constant scaling with AP_ACTIVITY

// ---- Glomerular deposition / injury kinetics ----
K_DEPOSIT0    = 0.00060  // deposit accrual rate constant
K_DEPOSIT_CLR = 0.00030  // deposit clearance/turnover rate
K_MESANGIAL0  = 0.00045
K_PODO0       = 0.00035
K_FIB0        = 0.00020
K_NEPHRON0    = 0.00018

// ---- Clinical/biomarker scaling ----
UPCR_MAX      = 6000     // mg/g ceiling
EGFR_MAX      = 100      // mL/min/1.73m2
EGFR_K        = 0.0025
SC5B9_MAX     = 900      // ng/mL ceiling (~10x ULN)
SC5B9_NORMAL  = 250

// ---- Transplant recurrence flag ----
TRANSPLANT    = 0        // 1 = post-transplant state (allograft, residual systemic AP dysregulation persists)

$CMT @annotated
IPTA_GUT    : Iptacopan gut depot (mg)
IPTA_CENT   : Iptacopan plasma compartment (mg)
PEG_SC      : Pegcetacoplan SC depot (mg)
PEG_CENT    : Pegcetacoplan plasma compartment (mg)
ECU_CENT    : Eculizumab central compartment (mg)
ECU_PERIPH  : Eculizumab peripheral compartment (mg)
RAVU_CENT   : Ravulizumab central compartment (mg)
RAVU_PERIPH : Ravulizumab peripheral compartment (mg)
DANI_GUT    : Danicopan gut depot (mg)
DANI_CENT   : Danicopan plasma compartment (mg)
AP_ACTIVITY : Alternative-pathway convertase activity index (0-1+, relative)
C3_LEVEL    : Serum C3 concentration (mg/dL)
SC5B9       : Soluble C5b-9 terminal complex (ng/mL)
GLOM_DEPOSIT: Glomerular C3b/deposit burden (relative units, 0-1)
MESANGIAL   : Mesangial proliferation / histologic activity index (0-1)
PODO_FRAC   : Viable podocyte fraction (0-1)
FIBROSIS    : Interstitial fibrosis / IF-TA index (0-1)
NEPHRON_FRAC: Surviving functional nephron fraction (0-1)
UPCR        : Urine protein-to-creatinine ratio (mg/g)
EGFR        : Estimated GFR, CKD-EPI (mL/min/1.73m2)
FX_YEARS    : Elapsed simulation time tracker (years)

$MAIN
double IPTA_Cp = IPTA_CENT/IPTA_V1;
double PEG_Cp  = PEG_CENT/PEG_V1;
double ECU_Cp  = ECU_CENT/ECU_V1;
double RAVU_Cp = RAVU_CENT/RAVU_V1;
double DANI_Cp = DANI_CENT/DANI_V1;

// ---- Upstream AP-pathway inhibitor effects (act on convertase formation) ----
double ipta_effect = IPTA_EMAX*(IPTA_Cp/(IPTA_EC50+IPTA_Cp+1e-9));
double peg_effect  = PEG_EMAX*(PEG_Cp/(PEG_EC50+PEG_Cp+1e-9));
double dani_effect = DANI_EMAX*(DANI_Cp/(DANI_EC50+DANI_Cp+1e-9));
double ap_upstream_block = 1.0 - (1.0-ipta_effect)*(1.0-peg_effect)*(1.0-dani_effect);
if (ap_upstream_block > 0.99) ap_upstream_block = 0.99;

// ---- Terminal-pathway (C5) inhibitor effects (act on MAC/sC5b-9 only) ----
double ecu_effect  = ECU_EMAX*(ECU_Cp/(ECU_EC50+ECU_Cp+1e-9));
double ravu_effect = RAVU_EMAX*(RAVU_Cp/(RAVU_EC50+RAVU_Cp+1e-9));
double c5_block = 1.0 - (1.0-ecu_effect)*(1.0-ravu_effect);
if (c5_block > 0.99) c5_block = 0.99;

// ---- Composite AP disease drive ----
double genetic_drive = GENETIC_GAIN*(1.0 + 0.6*CFH_LOF)*(1.0 + 0.5*C3NEF_TITER)*SEVERITY;
double ap_target = AP_BASELINE + AP_GAIN*genetic_drive*(1.0-ap_upstream_block);
if (ap_target < 0) ap_target = 0;

$ODE
// ---------------- Iptacopan PK ----------------
dxdt_IPTA_GUT  = -IPTA_KA*IPTA_GUT;
dxdt_IPTA_CENT =  IPTA_KA*IPTA_GUT*IPTA_F - (IPTA_CL/IPTA_V1)*IPTA_CENT;

// ---------------- Pegcetacoplan PK ----------------
dxdt_PEG_SC   = -PEG_KA*PEG_SC;
dxdt_PEG_CENT =  PEG_KA*PEG_SC*PEG_F - (PEG_CL/PEG_V1)*PEG_CENT;

// ---------------- Eculizumab PK (2-cpt mAb) ----------------
dxdt_ECU_CENT   = -(ECU_CL/ECU_V1)*ECU_CENT - (ECU_Q/ECU_V1)*ECU_CENT + (ECU_Q/ECU_VP)*ECU_PERIPH;
dxdt_ECU_PERIPH =  (ECU_Q/ECU_V1)*ECU_CENT - (ECU_Q/ECU_VP)*ECU_PERIPH;

// ---------------- Ravulizumab PK (2-cpt mAb, FcRn-recycled) ----------------
dxdt_RAVU_CENT   = -(RAVU_CL/RAVU_V1)*RAVU_CENT - (RAVU_Q/RAVU_V1)*RAVU_CENT + (RAVU_Q/RAVU_VP)*RAVU_PERIPH;
dxdt_RAVU_PERIPH =  (RAVU_Q/RAVU_V1)*RAVU_CENT - (RAVU_Q/RAVU_VP)*RAVU_PERIPH;

// ---------------- Danicopan PK ----------------
dxdt_DANI_GUT  = -DANI_KA*DANI_GUT;
dxdt_DANI_CENT =  DANI_KA*DANI_GUT*DANI_F - (DANI_CL/DANI_V1)*DANI_CENT;

// ---------------- Alternative-pathway convertase activity ----------------
dxdt_AP_ACTIVITY = AP_K_IN*(ap_target-AP_ACTIVITY);

// ---------------- Serum C3 (consumed by AP activity, synthesized by liver) ----------------
double c3_consumption = C3_CONSUME_K*AP_ACTIVITY*C3_LEVEL;
dxdt_C3_LEVEL = C3_SYNTH*(1.0-C3_LEVEL/(C3_NORMAL*1.3)) - c3_consumption;

// ---------------- Soluble C5b-9 (terminal pathway, blocked by C5 inhibitors only) ----------------
double sc5b9_target = SC5B9_NORMAL + (SC5B9_MAX-SC5B9_NORMAL)*AP_ACTIVITY*(1.0-c5_block);
dxdt_SC5B9 = 0.05*(sc5b9_target-SC5B9);

// ---------------- Glomerular C3b/deposit burden (upstream: NOT reduced by C5i alone) ----------------
double k_deposit = K_DEPOSIT0*AP_ACTIVITY;
dxdt_GLOM_DEPOSIT = k_deposit*(1.0-GLOM_DEPOSIT) - K_DEPOSIT_CLR*GLOM_DEPOSIT;

// ---------------- Mesangial proliferation / histologic activity ----------------
double k_mesangial = K_MESANGIAL0*GLOM_DEPOSIT*(1.0+0.3*SC5B9/SC5B9_NORMAL);
dxdt_MESANGIAL = k_mesangial*(1.0-MESANGIAL) - 0.00020*MESANGIAL;

// ---------------- Podocyte injury ----------------
double k_podo = K_PODO0*(GLOM_DEPOSIT+0.5*MESANGIAL)*(1.0-0.15*c5_block);
if (k_podo < 0) k_podo = 0;
dxdt_PODO_FRAC = -k_podo*PODO_FRAC;

// ---------------- Interstitial fibrosis ----------------
double fibrotic_drive = K_FIB0*(MESANGIAL+(1.0-PODO_FRAC))*(1.0+0.4*TRANSPLANT);
dxdt_FIBROSIS = fibrotic_drive*(1.0-FIBROSIS) - 0.00003*FIBROSIS;
if (FIBROSIS > 1.0) FIBROSIS = 1.0;

// ---------------- Nephron loss ----------------
double k_nephron = K_NEPHRON0*(1.0+2.0*FIBROSIS)*(1.0+0.5*(1.0-PODO_FRAC));
dxdt_NEPHRON_FRAC = -k_nephron*NEPHRON_FRAC;

// ---------------- Proteinuria (UPCR) ----------------
double upcr_target = UPCR_MAX*(1.0-PODO_FRAC)*(1.0+0.5*MESANGIAL);
if (upcr_target > UPCR_MAX) upcr_target = UPCR_MAX;
if (upcr_target < 80) upcr_target = 80;
dxdt_UPCR = 0.012*(upcr_target-UPCR);

// ---------------- eGFR ----------------
double egfr_target = EGFR_MAX*NEPHRON_FRAC*PODO_FRAC*(1.0-0.5*FIBROSIS);
if (egfr_target < 2) egfr_target = 2;
if (egfr_target > EGFR_MAX) egfr_target = EGFR_MAX;
dxdt_EGFR = EGFR_K*(egfr_target-EGFR);

dxdt_FX_YEARS = 1.0/(24.0*365.0);

$INIT
IPTA_GUT = 0, IPTA_CENT = 0, PEG_SC = 0, PEG_CENT = 0,
ECU_CENT = 0, ECU_PERIPH = 0, RAVU_CENT = 0, RAVU_PERIPH = 0,
DANI_GUT = 0, DANI_CENT = 0,
AP_ACTIVITY = 0.05, C3_LEVEL = 100, SC5B9 = 260, GLOM_DEPOSIT = 0,
MESANGIAL = 0, PODO_FRAC = 1.0, FIBROSIS = 0, NEPHRON_FRAC = 1.0,
UPCR = 100, EGFR = 95, FX_YEARS = 0

$CAPTURE ipta_effect peg_effect dani_effect ap_upstream_block ecu_effect ravu_effect
$CAPTURE c5_block genetic_drive ap_target c3_consumption sc5b9_target k_deposit
$CAPTURE k_mesangial k_podo fibrotic_drive k_nephron upcr_target egfr_target

## =============================================================================
## Treatment scenarios (see c3g_shiny_app.R for interactive dosing UI)
##
## 1. Natural history, C3GN phenotype (SEVERITY=1.0, C3NEF_TITER=1.0, untreated)
## 2. Natural history, DDD phenotype (SEVERITY=1.35, C3NEF_TITER=2.0, high-
##    titer persistent C3NeF, untreated -> fastest progression)
## 3. Iptacopan 200 mg oral BID (Factor B inhibitor; APPEAR-C3G regimen)
## 4. Pegcetacoplan 1080 mg SC twice-weekly (C3/C3b inhibitor)
## 5. Eculizumab 900 mg IV q2w maintenance (off-label; high-sC5b-9 responder
##    subgroup, upstream C3 deposition persists)
## 6. Ravulizumab weight-based IV q8w maintenance (engineered FcRn recycling,
##    extended terminal-pathway blockade)
## 7. Danicopan 150 mg oral TID add-on to eculizumab background (breakthrough
##    residual AP activity)
## 8. CFH loss-of-function genotype (CFH_LOF=1) + iptacopan
## 9. Post-transplant recurrence (TRANSPLANT=1) with pegcetacoplan prophylaxis
## 10. Iptacopan dose reduction/discontinuation after 2 years -> relapse
##
## Example mrgsolve event code (iptacopan chronic BID dosing, C3GN):
##   mod <- mread("c3g_mrgsolve_model") %>% param(SEVERITY = 1.0, C3NEF_TITER = 1.0)
##   e_ipta <- ev(amt = 200, cmt = "IPTA_GUT", time = 0, ii = 12, addl = 2*365*2-1)
##   out <- mod %>% ev(e_ipta) %>% mrgsim(end = 24*365*10, delta = 24)  # 10-yr horizon
## =============================================================================
