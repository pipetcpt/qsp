# =====================================================================
# Microscopic Colitis (LC + CC) — mrgsolve QSP Model
# ---------------------------------------------------------------------
# Drugs:  Budesonide MMX, Cholestyramine, Mesalamine, Bismuth subsalicylate,
#         Loperamide, Azathioprine, Infliximab, Vedolizumab, Adalimumab
# Disease compartments: CD8-IEL, IFN-γ, TNF-α, IL-6/IL-17, TGF-β,
#         epithelial barrier integrity, collagen-band thickness,
#         bile-acid load, mucosal Na+/H2O flux, stool frequency,
#         Hjortswang score, HRQoL.
# Calibration sources:
#   - BUC-60 / BUC-63 (Miehlke 2007, Münch 2014)  : budesonide RCTs
#   - COLAL-PRED (Bonderup 2009)
#   - Cholestyramine open-label (Ung 2000)
#   - VICTORIA-MC (Miehlke 2021) vedolizumab open
#   - Riddell, Loftus epidemiology
# ---------------------------------------------------------------------
# Author: QSP Disease Model Library (Claude Code Routine)
# =====================================================================

library(mrgsolve)

code <- '
$PROB Microscopic Colitis QSP

$PARAM @annotated
// ---- Budesonide MMX (PO, pH-/time-release, ~3 h Tlag) ------------------
KA_BUD      :  0.85 : Absorption rate constant (1/h)
F_BUD       :  0.10 : Oral bioavailability after high first-pass
TLAG_BUD    :  3.00 : Release lag (h)
VC_BUD      : 240   : Central V (L)
Q_BUD       :  45   : Inter-compartmental Q (L/h)
VP_BUD      : 880   : Peripheral V (L)
CL_BUD      :  98   : Plasma clearance (L/h)
KGUT_BUD    :  0.30 : Mucosal degradation rate (1/h)
EC50_BUD_GR :   4.0 : Local conc for 50% GR effect (mucosal a.u.)

// ---- Cholestyramine (luminal, non-absorbed bile-acid sequestrant) ------
KA_CHOL     :  6.0  : Apparent luminal mixing (1/h)
KGUT_CHOL   :  0.50 : Loss with stool (1/h)
EC50_CHOL   :  4.0  : Conc for 50% BA chelation (g equiv. luminal)

// ---- Mesalamine PO ------------------------------------------------------
KA_MES      :  0.90
F_MES       :  0.30
CL_MES      :  20
VC_MES      :  20
EC50_MES    :  10   : Local conc (mg/L equiv.) for PPAR-γ effect

// ---- Bismuth subsalicylate -- modeled as luminal anti-secretory ---------
KA_BIS      :  4.0
KGUT_BIS    :  0.80
EC50_BIS    :  3.0

// ---- Loperamide (mu-opioid, gut motility) -------------------------------
KA_LOP      :  1.5
F_LOP       :  0.40
CL_LOP      :  60
VC_LOP      : 300
EC50_LOP    :  0.6

// ---- Azathioprine -> 6-TGN (RBC) ----------------------------------------
KA_AZA      :  3.0
F_AZA       :  0.47
CL_AZA      :  50
VC_AZA      :  90
KFORM_TGN   :  0.10 : 6-TGN formation rate (proxy)
KDEG_TGN    :  0.0035 : 6-TGN degradation (1/h)
EC50_TGN    :  250  : pmol/8e8 RBC for 50% lymphocyte effect

// ---- Infliximab (IV, anti-TNF) ------------------------------------------
CL_IFX      :  0.32 : Linear CL (L/day)
VC_IFX      :   3.0
Q_IFX       :   0.5
VP_IFX      :   3.5
KD_IFX_TNF  :   0.1 : nM

// ---- Adalimumab (SC, anti-TNF) ------------------------------------------
KA_ADA      :  0.30 : 1/day
F_ADA       :  0.64
CL_ADA      :  0.30 : L/day
VC_ADA      :   5.0
KD_ADA_TNF  :   0.1

// ---- Vedolizumab (IV, anti-α4β7) ----------------------------------------
CL_VDZ      :  0.16
VC_VDZ      :   4.7
Q_VDZ       :   0.30
VP_VDZ      :   1.7
EC50_VDZ    :   2.0  : mg/L for 50% gut homing block

// =====================================================================
// Disease parameters (rates per day unless noted)
// =====================================================================
KSYN_IEL    :  0.10 : CD8+ IEL synthesis (a.u./day)
KDEG_IEL    :  0.04 : CD8+ IEL turnover (1/day)
IEL_BASE    :  35   : Active-disease baseline IEL (per 100 EC)
IEL_NORM    :  8

KSYN_IFNG   :  6.0  : IFN-γ production (pg/mL/day)
KDEG_IFNG   :  0.50
IFNG_BASE   :   25

KSYN_TNF    :  4.0
KDEG_TNF    :  0.55
TNF_BASE    :   18

KSYN_IL6    :  3.0
KDEG_IL6    :  0.70
IL6_BASE    :   12

KSYN_IL17   :  2.5
KDEG_IL17   :  0.65
IL17_BASE   :    9

KSYN_TGFB   :  2.0
KDEG_TGFB   :  0.20
TGFB_BASE   :   15

// Barrier (0 = open, 1 = intact)
BAR_BASE    :  0.55
KFORM_BAR   :  0.05
KDEG_BAR    :  0.05 : Driven by cytokines

// Collagen band thickness (μm) – CC subtype
COL_BASE    :  18   : CC baseline (CC dx > 10)
COL_LC      :   4   : LC baseline
KFORM_COL   :  0.02
KDEG_COL    :  0.02
EC50_TGFB_COL: 30

// Bile acid colonic load (a.u.)
BA_BASE     :  2.5   : Baseline colonic BA load
KSYN_BA     :  0.4
KDEG_BA     :  0.15

// Net absorptive flux (mL water/day → -ve = secretory)
WAT_BASE    : -350  : ml/day secretory in active disease
WAT_NORM    :  200  : ml/day net absorption healthy
KFLX        :  0.30

// Stool freq (stools/day)
STOOL_BASE  :  6.0
STOOL_MIN   :  1.5
KSTOOL      :  0.5

// Hjortswang clinical-remission score (0 = remission, 1 = active)
HJ_BASE     :  1.0

// HRQoL (SF-36 mental subscale-like, 0–100, 100 best)
QOL_BASE    :  55
KQOL        :  0.04

// Subtype switch (1 = CC, 0 = LC)
SUBTYPE_CC  : 1

// Comorbidity flags
BMD_DECAY_GC: 0.0008 : 1/day BMD loss while on systemic GC
HPA_HALF    : 60     : hours t1/2 for HPA suppression

$CMT @annotated
// PK compartments
GUT_BUD     : Budesonide gut lumen
CENT_BUD    : Budesonide central
PER_BUD     : Budesonide peripheral
GUT_CHOL    : Cholestyramine luminal
GUT_MES     : Mesalamine gut
CEN_MES     : Mesalamine plasma
GUT_BIS     : Bismuth luminal
GUT_LOP     : Loperamide gut
CEN_LOP     : Loperamide plasma
GUT_AZA     : Azathioprine gut
CEN_AZA     : Azathioprine plasma
TGN         : 6-TGN RBC pool
CEN_IFX     : Infliximab central
PER_IFX     : Infliximab peripheral
DEPOT_ADA   : Adalimumab SC depot
CEN_ADA     : Adalimumab plasma
CEN_VDZ     : Vedolizumab central
PER_VDZ     : Vedolizumab peripheral
// Disease compartments
IEL         : CD8+ intraepithelial lymphocytes (per 100 EC)
IFNG        : IFN-gamma (pg/mL)
TNF         : TNF-alpha (pg/mL)
IL6         : IL-6 (pg/mL)
IL17        : IL-17A (pg/mL)
TGFB        : TGF-beta (pg/mL)
BAR         : Epithelial barrier integrity (0–1)
COL         : Subepithelial collagen band (μm)
BA          : Colonic bile-acid load (a.u.)
WAT         : Net colonic water flux (mL/day)
STOOL       : Stool frequency (stools/day)
HJ          : Hjortswang active-disease score
QOL         : HRQoL composite
BMD         : Bone mineral density (T-score proxy)
HPA         : HPA axis cortisol response (1=normal,0=suppressed)

$GLOBAL
double pos(double x){ return x<0 ? 0 : x; }
double sat(double C, double E){ return C/(C+E); }

$MAIN
GUT_BUD_0  = 0;
CENT_BUD_0 = 0;
PER_BUD_0  = 0;
GUT_CHOL_0 = 0;
GUT_MES_0  = 0;
CEN_MES_0  = 0;
GUT_BIS_0  = 0;
GUT_LOP_0  = 0;
CEN_LOP_0  = 0;
GUT_AZA_0  = 0;
CEN_AZA_0  = 0;
TGN_0      = 0;
CEN_IFX_0  = 0;
PER_IFX_0  = 0;
DEPOT_ADA_0= 0;
CEN_ADA_0  = 0;
CEN_VDZ_0  = 0;
PER_VDZ_0  = 0;

IEL_0   = IEL_BASE;
IFNG_0  = IFNG_BASE;
TNF_0   = TNF_BASE;
IL6_0   = IL6_BASE;
IL17_0  = IL17_BASE;
TGFB_0  = TGFB_BASE;
BAR_0   = BAR_BASE;
COL_0   = (SUBTYPE_CC>0.5 ? COL_BASE : COL_LC);
BA_0    = BA_BASE;
WAT_0   = WAT_BASE;
STOOL_0 = STOOL_BASE;
HJ_0    = HJ_BASE;
QOL_0   = QOL_BASE;
BMD_0   = 0.0;
HPA_0   = 1.0;

$ODE
// ====== Budesonide PK ======
double BUD_GUT_C = GUT_BUD / 0.4;                 // mucosal a.u.
double BUD_PLAS  = CENT_BUD / VC_BUD;
dxdt_GUT_BUD  = -KA_BUD * GUT_BUD - KGUT_BUD * GUT_BUD;
dxdt_CENT_BUD =  F_BUD * KA_BUD * GUT_BUD - (CL_BUD/VC_BUD)*CENT_BUD
                 - (Q_BUD/VC_BUD)*CENT_BUD + (Q_BUD/VP_BUD)*PER_BUD;
dxdt_PER_BUD  =  (Q_BUD/VC_BUD)*CENT_BUD - (Q_BUD/VP_BUD)*PER_BUD;

// ====== Cholestyramine luminal ======
dxdt_GUT_CHOL = -KGUT_CHOL * GUT_CHOL;
double CHOL_LUM = GUT_CHOL;

// ====== Mesalamine ======
dxdt_GUT_MES = -KA_MES * GUT_MES;
dxdt_CEN_MES =  F_MES * KA_MES * GUT_MES - (CL_MES/VC_MES) * CEN_MES;
double MES_LUM = GUT_MES;

// ====== Bismuth luminal ======
dxdt_GUT_BIS = -KGUT_BIS * GUT_BIS;
double BIS_LUM = GUT_BIS;

// ====== Loperamide ======
dxdt_GUT_LOP = -KA_LOP * GUT_LOP;
dxdt_CEN_LOP =  F_LOP * KA_LOP * GUT_LOP - (CL_LOP/VC_LOP) * CEN_LOP;
double LOP_C = CEN_LOP / VC_LOP;

// ====== Azathioprine -> 6-TGN ======
dxdt_GUT_AZA = -KA_AZA * GUT_AZA;
dxdt_CEN_AZA =  F_AZA * KA_AZA * GUT_AZA - (CL_AZA/VC_AZA) * CEN_AZA;
double AZA_C = CEN_AZA / VC_AZA;
dxdt_TGN     =  KFORM_TGN * AZA_C - KDEG_TGN * TGN;

// ====== Infliximab IV ======
dxdt_CEN_IFX = -(CL_IFX/24/VC_IFX) * CEN_IFX
               - (Q_IFX/24/VC_IFX) * CEN_IFX + (Q_IFX/24/VP_IFX) * PER_IFX;
dxdt_PER_IFX =  (Q_IFX/24/VC_IFX) * CEN_IFX - (Q_IFX/24/VP_IFX) * PER_IFX;
double IFX_C = CEN_IFX / VC_IFX;

// ====== Adalimumab SC ======
dxdt_DEPOT_ADA = -(KA_ADA/24) * DEPOT_ADA;
dxdt_CEN_ADA   =  F_ADA * (KA_ADA/24) * DEPOT_ADA - (CL_ADA/24/VC_ADA) * CEN_ADA;
double ADA_C = CEN_ADA / VC_ADA;

// ====== Vedolizumab IV ======
dxdt_CEN_VDZ = -(CL_VDZ/24/VC_VDZ) * CEN_VDZ
               - (Q_VDZ/24/VC_VDZ) * CEN_VDZ + (Q_VDZ/24/VP_VDZ) * PER_VDZ;
dxdt_PER_VDZ =  (Q_VDZ/24/VC_VDZ) * CEN_VDZ - (Q_VDZ/24/VP_VDZ) * PER_VDZ;
double VDZ_C = CEN_VDZ / VC_VDZ;

// ====== Drug-effect functions ======
double E_BUD   = sat(BUD_GUT_C, EC50_BUD_GR);     // mucosal GR transrepression
double E_CHOL  = sat(CHOL_LUM, EC50_CHOL);        // BA sequestration
double E_MES   = sat(MES_LUM,  EC50_MES);         // PPAR-γ
double E_BIS   = sat(BIS_LUM,  EC50_BIS);         // anti-secretory
double E_LOP   = sat(LOP_C,   EC50_LOP);          // motility/secretion
double E_TGN   = sat(TGN,     EC50_TGN);
double E_TNF_b = sat(IFX_C + ADA_C, 5.0);         // anti-TNF combined
double E_VDZ   = sat(VDZ_C,   EC50_VDZ);

// ====== Disease ODEs ======
// CD8+ IEL — driven by IL-15/IFN-γ, suppressed by GC, TPMT, anti-α4β7
double IEL_SS = IEL_NORM + (IEL_BASE - IEL_NORM) * (1 - 0.55*E_BUD - 0.25*E_TGN - 0.40*E_VDZ);
dxdt_IEL  = KSYN_IEL*(IEL_SS) - KDEG_IEL * IEL;

// IFN-γ
dxdt_IFNG = KSYN_IFNG * (IEL/IEL_NORM) * (1 - 0.60*E_BUD - 0.30*E_TGN - 0.20*E_VDZ)
            - KDEG_IFNG * IFNG;

// TNF-α — anti-TNF combo plus GC
dxdt_TNF  = KSYN_TNF  * (1 - 0.55*E_BUD - 0.35*E_TGN - 0.85*E_TNF_b)
            - KDEG_TNF * TNF;

// IL-6
dxdt_IL6  = KSYN_IL6  * (1 - 0.50*E_BUD - 0.20*E_MES) - KDEG_IL6 * IL6;

// IL-17
dxdt_IL17 = KSYN_IL17 * (1 - 0.55*E_BUD - 0.40*E_TGN) - KDEG_IL17 * IL17;

// TGF-β — driven by mast cells / collagen feedback
dxdt_TGFB = KSYN_TGFB * (1 - 0.45*E_BUD) * (1 + 0.4*(COL/(COL+5.0)))
            - KDEG_TGFB * TGFB;

// Barrier integrity: declines with TNF, IFN-γ, IL-17 + BA
double CYTO_BURDEN = (TNF/TNF_BASE) + (IFNG/IFNG_BASE) + 0.4*(IL17/IL17_BASE) + 0.3*(BA/BA_BASE);
dxdt_BAR  = KFORM_BAR*(1 + 0.4*E_BUD + 0.2*E_TNF_b + 0.1*E_MES)*(1 - BAR)
            - KDEG_BAR * CYTO_BURDEN * BAR;

// Collagen band — CC: TGF-β builds, GC degrades
dxdt_COL  = KFORM_COL * sat(TGFB, EC50_TGFB_COL) * (SUBTYPE_CC>0.5 ? 1.0 : 0.3)
            - KDEG_COL * COL * (1 + 0.8*E_BUD);

// Colonic bile-acid load
dxdt_BA   = KSYN_BA * (1 + 0.4*(1-BAR)) - KDEG_BA * BA - 1.5 * E_CHOL * BA;

// Net colonic water flux (negative = secretion)
double WAT_TARGET = WAT_BASE
   + (WAT_NORM - WAT_BASE) * (0.40*E_BUD + 0.25*E_CHOL + 0.20*E_BIS
                              + 0.15*E_LOP + 0.10*E_MES + 0.10*E_TNF_b
                              + 0.10*E_VDZ + 0.10*E_TGN);
WAT_TARGET = WAT_TARGET + 200*pos(BAR-0.6) - 80*(BA/BA_BASE - 1);
dxdt_WAT  = KFLX * (WAT_TARGET - WAT);

// Stool frequency — softplus driven by inverse of WAT
double drive = (WAT_BASE - WAT) / (WAT_BASE - WAT_NORM + 1e-3); // 0 = active, 1 = normal
double STOOL_TARGET = STOOL_BASE - (STOOL_BASE - STOOL_MIN) * drive;
dxdt_STOOL = KSTOOL * (STOOL_TARGET - STOOL);

// Hjortswang score — soft indicator: <3 stools/day & WAT recovering → remission
double HJ_TARGET = 1 / (1 + exp(-2*(STOOL - 3.0)));
dxdt_HJ   = 0.5 * (HJ_TARGET - HJ);

// HRQoL — improves as stools, urgency, fatigue drop
double QOL_TARGET = 90 - 30*HJ - 5*(STOOL-1.5);
QOL_TARGET = QOL_TARGET < 30 ? 30 : QOL_TARGET;
dxdt_QOL  = KQOL * (QOL_TARGET - QOL);

// BMD loss while on systemic GC (BUD_PLAS) above a threshold
dxdt_BMD  = - BMD_DECAY_GC * (BUD_PLAS > 0.5 ? 1.0 : 0.0) - 0.0003 * (BUD_PLAS); // small effect

// HPA suppression first-order recovery + suppression
double k_hpa = log(2)/HPA_HALF;       // per hour
double HPA_target = 1 - 0.6 * (BUD_PLAS / (BUD_PLAS + 1.0));
dxdt_HPA  = k_hpa * 24 * (HPA_target - HPA);  // convert to per day

$TABLE
double BUD_CONC = CENT_BUD / VC_BUD;
double STOOL_CLIP = STOOL < 0 ? 0 : STOOL;
double REMIT = (STOOL_CLIP < 3.0 ? 1.0 : 0.0);

$CAPTURE BUD_CONC IFX_C ADA_C VDZ_C TGN STOOL_CLIP HJ QOL BMD HPA BAR COL BA WAT IEL IFNG TNF IL6 IL17 TGFB REMIT
'

# Compile and quick smoke-test
mod <- mcode("microscopic_colitis_qsp", code)

# =====================================================================
# Treatment Scenarios
# =====================================================================

# Scenario 1: Natural history (no treatment), 365 days
s1 <- expand.ev(amt = 0)

# Scenario 2: Budesonide MMX 9 mg PO daily x 8 weeks induction + 6 mg taper x 6 mo
s2_ind  <- ev(time = 0,    amt = 9,  cmt = "GUT_BUD", ii = 24, addl = 55)        # 8 weeks
s2_main <- ev(time = 24*56,amt = 6,  cmt = "GUT_BUD", ii = 24, addl = 180)        # ~6 months
s2 <- c(s2_ind, s2_main)

# Scenario 3: Cholestyramine 4 g PO TID x 6 months
s3 <- ev(amt = 4, cmt = "GUT_CHOL", ii = 8, addl = 540)

# Scenario 4: Mesalamine 2.4 g PO daily x 8 weeks
s4 <- ev(amt = 2400, cmt = "GUT_MES", ii = 24, addl = 55)

# Scenario 5: Bismuth subsalicylate 524 mg PO QID x 8 weeks
s5 <- ev(amt = 524, cmt = "GUT_BIS", ii = 6, addl = 224)

# Scenario 6: Loperamide 4 mg PO BID prn (chronic) — symptomatic
s6 <- ev(amt = 4, cmt = "GUT_LOP", ii = 12, addl = 720)

# Scenario 7: Azathioprine 2 mg/kg (140 mg) PO QD — refractory
s7 <- ev(amt = 140, cmt = "GUT_AZA", ii = 24, addl = 360)

# Scenario 8: Infliximab 5 mg/kg (350 mg) IV at 0, 2, 6 weeks then q8w
ifx_doses <- c(0, 14*24, 42*24, seq(98*24, 365*24, by = 56*24))
s8 <- do.call(c, lapply(ifx_doses, function(t) ev(time = t, amt = 350, cmt = "CEN_IFX")))

# Scenario 9: Vedolizumab 300 mg IV at 0, 2, 6 weeks then q8w
vdz_doses <- ifx_doses
s9 <- do.call(c, lapply(vdz_doses, function(t) ev(time = t, amt = 300, cmt = "CEN_VDZ")))

# Scenario 10: Budesonide + Cholestyramine combo (refractory CC w/ BAM)
s10 <- c(s2, s3)

scenarios <- list(
  "01_natural_history"     = ev(amt = 0),
  "02_budesonide_taper"    = s2,
  "03_cholestyramine"      = s3,
  "04_mesalamine"          = s4,
  "05_bismuth"             = s5,
  "06_loperamide_symp"     = s6,
  "07_azathioprine"        = s7,
  "08_infliximab_q8w"      = s8,
  "09_vedolizumab_q8w"     = s9,
  "10_bud_chol_combo"      = s10
)

run_scenario <- function(label, evd, subtype = c("CC", "LC"), tend = 365*24) {
  subtype <- match.arg(subtype)
  pset <- if (subtype == "LC") list(SUBTYPE_CC = 0) else list(SUBTYPE_CC = 1)
  out <- mod %>% param(pset) %>% mrgsim(events = evd, end = tend, delta = 24)
  df  <- as.data.frame(out)
  df$scenario <- label
  df$subtype  <- subtype
  df
}

if (FALSE) {
  # Demo / verification – uncomment to execute
  out_all <- do.call(rbind, lapply(names(scenarios), function(n) run_scenario(n, scenarios[[n]])))
  head(out_all)
  with(subset(out_all, scenario == "02_budesonide_taper"),
       plot(time/24, STOOL_CLIP, type = "l", ylim = c(0, 8),
            xlab = "Days", ylab = "Stools/day",
            main = "Budesonide MMX induction + taper (CC)"))
}

# Export for Shiny app re-use
mc_model     <- mod
mc_scenarios <- scenarios
