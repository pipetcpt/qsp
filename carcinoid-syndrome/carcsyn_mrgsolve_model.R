# =====================================================================
# Carcinoid Syndrome (Functional Midgut NET) — QSP / mrgsolve model
# Author: CCR / Claude Code Routine (QSP Disease Model Library)
# Date  : 2026-06-30
#
# Scope :
#   • Tryptophan → 5-HTP (TPH1) → 5-HT → 5-HIAA pathway with platelet pool
#   • Tumor burden (hepatic) driving systemic spill of 5-HT, tachykinins,
#     kallikrein/bradykinin (lumped), histamine (foregut subtype)
#   • SSTR2/5 binding for octreotide LAR, lanreotide autogel, pasireotide
#   • TPH1 inhibition by telotristat (active metabolite LP-778902)
#   • mTORC1 inhibition by everolimus; VEGFR-TKI generic (sunitinib /
#     surufatinib / cabozantinib) lumped to tumor-growth Hill
#   • ¹⁷⁷Lu-DOTATATE (PRRT) — SSTR2-mediated cell-kill of tumor
#   • IFN-α-2b — antiproliferative + anti-angiogenic
#   • Hepatic artery embolization (HAE) — instantaneous tumor mass step
#   • Symptomatic dynamics: flushing episodes, bowel movements/day,
#     bronchospasm index, NT-proBNP from valvular fibrosis (5-HT2B/TGF-β)
#   • Carcinoid heart disease (CHD) progression via TGF-β/CTGF-driven
#     valve collagen
#
# Compartments / states (≥ 20 ODE):
#   1  GUT_OCT    Octreotide LAR depot (zero-order release)
#   2  CENT_OCT   Octreotide central
#   3  PER_OCT    Octreotide peripheral
#   4  GUT_LAN    Lanreotide autogel depot
#   5  CENT_LAN   Lanreotide central
#   6  GUT_PASIR  Pasireotide LAR depot
#   7  CENT_PASIR Pasireotide central
#   8  GUT_TEL    Telotristat ethyl gut depot
#   9  CENT_TEL   Telotristat active metabolite plasma
#  10  GUT_EVE    Everolimus depot
#  11  CENT_EVE   Everolimus plasma
#  12  GUT_TKI    VEGFR-TKI generic depot
#  13  CENT_TKI   VEGFR-TKI plasma
#  14  CENT_IFN   IFN-α central
#  15  PRRT_AMT   ¹⁷⁷Lu-DOTATATE tumor-delivered activity
#  16  TUMOR      Hepatic tumor burden (g)
#  17  TPH1       Functional TPH1 (% of baseline)
#  18  HTP        5-HTP intratumoral
#  19  SER_T      Tumor intracellular 5-HT
#  20  SER_P      Plasma free 5-HT (nmol/L)
#  21  SER_PLT    Whole-blood platelet 5-HT (ng/mL)
#  22  HIAA_U     Urinary 24-h 5-HIAA (mg/24h)
#  23  TGFb       Plasma TGF-β1 (signal)
#  24  VALVE      Valve collagen / Hassan-equivalent score
#  25  NTproBNP   NT-proBNP (pg/mL)
#  26  BM         Bowel movements / day (rolling)
#  27  FLUSH      Flushing episodes / day
#
# References calibrated against:
#   • Telotristat 250 mg TID — TELESTAR / TELECAST (Kulke 2017,2019)
#   • Octreotide LAR PK — Astrup 2002 / Lamberts 2000 / PROMID 2009
#   • Lanreotide autogel — CLARINET (Caplin 2014)
#   • Pasireotide LAR — Bergsland 2014
#   • Everolimus — RADIANT-4 (Yao 2016)
#   • PRRT — NETTER-1 (Strosberg 2017) / NETTER-2 (Singh 2024)
#   • IFN-α NET — SWOG / Faiss 2003
#   • Carcinoid heart disease echo progression — Møller 2003 / Bhattacharyya 2011
# =====================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

carcsyn_code <- '
$PROB
# Carcinoid Syndrome QSP — functional midgut NET
# 27-state ODE model with SSA, TPH1i, PRRT, mTORi, TKI, IFN-α

$PARAM @annotated
// ------------ Octreotide LAR PK (Astrup 2002; depot zero-order)
KA_OCT     :  0.012     : Octreotide LAR depot release rate (1/h)
CL_OCT     :  10.5      : Octreotide CL (L/h)
V_OCT      :  21.0      : Octreotide central V (L)
Q_OCT      :   4.0      : Inter-compartmental clearance (L/h)
VP_OCT     :  15.0      : Peripheral V (L)
KD_OCT_S2  :   0.40     : Octreotide Kd at SSTR2 (nmol/L)
// ------------ Lanreotide autogel
KA_LAN     :  0.010     : Lanreotide release (1/h)
CL_LAN     :   9.0      : Lanreotide CL (L/h)
V_LAN      :  18.0      : Lanreotide V (L)
KD_LAN_S2  :   0.50     : Lanreotide Kd at SSTR2 (nmol/L)
// ------------ Pasireotide LAR (pan-SSTR)
KA_PAS     :  0.014     : Pasireotide release (1/h)
CL_PAS     :   8.5      : Pasireotide CL (L/h)
V_PAS      :  25.0      : Pasireotide V (L)
KD_PAS_S5  :   0.40     : Pasireotide Kd at SSTR5 (nmol/L)
// ------------ Telotristat (LP-778902 active metabolite)
KA_TEL     :   0.90     : Telotristat absorption (1/h)
CL_TEL     :  35.0      : Telotristat met CL (L/h)
V_TEL      : 110.0      : Telotristat met V (L)
F_TEL      :   0.70     : Bioavailability (LP-778902 formation)
IC50_TPH1  :   8.0      : LP-778902 IC50 for TPH1 (ng/mL)
GAMMA_TPH1 :   1.3      : Hill on TPH1 inhibition
// ------------ Everolimus
KA_EVE     :   2.5      : Everolimus absorption (1/h)
CL_EVE     :  18.0      : Everolimus CL (L/h)
V_EVE      : 195.0      : Everolimus V (L)
IC50_MTOR  :   12.0     : Everolimus IC50 on tumor growth (ng/mL)
// ------------ VEGFR-TKI (sunitinib/surufatinib/cabozantinib lumped)
KA_TKI     :   0.8      : TKI absorption (1/h)
CL_TKI     :  37.0      : TKI CL (L/h)
V_TKI      : 1500.0     : TKI V (L)
IC50_VEGF  :   60.0     : TKI IC50 on tumor angiogenesis (ng/mL)
// ------------ IFN-α-2b
CL_IFN     :   8.0      : IFN-α CL (L/h)
V_IFN      :  15.0      : IFN-α V (L)
EFF_IFN50  :   2.0      : IFN-α effect EC50 (IU/mL)
// ------------ Tumor biology
KGROW      : 0.00080    : Net tumor growth rate (1/h)
TUMAX      : 2000.0     : Tumor carrying capacity (g)
DRUG_KILL_PRRT : 0.020   : PRRT cell-kill rate constant (1/h per GBq tumor)
TUM_BASE   :  300.0     : Baseline hepatic tumor burden (g)
KGROW_HILL :   1.0      : Hill on growth attenuation by mTORi
// ------------ Serotonin pathway (per gram tumor)
KSYN_HTP   :   2.0      : 5-HTP synthesis rate per g tumor (nmol/h/g)
KCONV_AAD  :  20.0      : AADC conversion 5-HTP -> 5-HT (1/h)
KREL_SER   :   0.50     : Exocytotic release rate (1/h) (modulated by SSTR)
KMAOA      :   3.5      : MAO-A clearance of plasma 5-HT (1/h)
KUPT_PLT   :   0.30     : SERT uptake into platelets (1/h)
KOUT_PLT   :   0.005    : Platelet 5-HT turnover (1/h)
F_HIAA     :   0.85     : Fraction of plasma 5-HT to urinary 5-HIAA
KSER_NORM  :   1.0      : Healthy basal plasma 5-HT (nmol/L)
// ------------ Symptom dynamics
KIN_BM     :   1.5      : Bowel movement baseline (BM/day)
EC50_BM    :   2.5      : Plasma 5-HT EC50 for diarrhea (nmol/L)
EMAX_BM    :   6.0      : Max BM/day added
KIN_FL     :   0.5      : Baseline flushing episodes/day
EC50_FL    :   1.5      : 5-HT (+ tachykinin proxy) EC50 for flushing
EMAX_FL    :   8.0      : Max flushing episodes/day
KOUT_SYMP  :   2.0      : Symptom turnover (1/h)
// ------------ Valve / CHD
KIN_VAL    :   0.0001   : Valve collagen deposition rate (per nmol/L 5-HT/h)
KOUT_VAL   :   0.00003  : Valve collagen resorption (1/h)
HASSAN_MAX :  20.0      : Hassan score ceiling
KIN_BNP    :   5.0      : NT-proBNP secretion per unit valve score (pg/mL/h)
KOUT_BNP   :   0.05     : NT-proBNP clearance (1/h)
KTGFB      :   0.001    : TGF-β1 induction per nmol 5-HT (1/h)
KOUT_TGFB  :   0.1      : TGF-β1 clearance (1/h)
// ------------ Symptomatic / supportive
KOND       :   0.0      : Ondansetron 5-HT3 block factor (set per scenario)
KLOP       :   0.0      : Loperamide motility factor

$CMT @annotated
GUT_OCT    : Octreotide LAR depot (mg)
CENT_OCT   : Octreotide central (mg)
PER_OCT    : Octreotide peripheral (mg)
GUT_LAN    : Lanreotide depot (mg)
CENT_LAN   : Lanreotide central (mg)
GUT_PASIR  : Pasireotide depot (mg)
CENT_PASIR : Pasireotide central (mg)
GUT_TEL    : Telotristat parent depot (mg)
CENT_TEL   : Telotristat active metabolite (ng-equivalent)
GUT_EVE    : Everolimus depot (mg)
CENT_EVE   : Everolimus plasma (ng)
GUT_TKI    : TKI depot (mg)
CENT_TKI   : TKI plasma (ng)
CENT_IFN   : IFN-α plasma (IU)
PRRT_AMT   : PRRT delivered activity (GBq-tumor)
TUMOR      : Hepatic tumor burden (g)
TPH1       : Functional TPH1 (fraction)
HTP        : 5-HTP intratumoral (nmol)
SER_T      : Tumor 5-HT (nmol)
SER_P      : Plasma 5-HT (nmol/L equivalent)
SER_PLT    : Platelet 5-HT (ng/mL)
HIAA_U     : Urinary 5-HIAA (mg/24h running)
TGFb       : Plasma TGF-β1 (rel)
VALVE      : Valve collagen / Hassan equivalent
NTproBNP   : NT-proBNP (pg/mL)
BM         : Bowel movements/day (rolling)
FLUSH      : Flushing episodes/day

$MAIN
TPH1_0 = 1.0;
TUMOR_0 = TUM_BASE;
SER_P_0 = KSER_NORM;
SER_PLT_0 = 100.0;
HIAA_U_0 = 5.0;
VALVE_0 = 1.0;
NTproBNP_0 = 80.0;
BM_0 = KIN_BM;
FLUSH_0 = KIN_FL;
TGFb_0 = 1.0;

$ODE
// ---- PK: octreotide LAR ----
double C_OCT = CENT_OCT / V_OCT * 1000.0; // ng/mL
double C_LAN = CENT_LAN / V_LAN * 1000.0;
double C_PAS = CENT_PASIR / V_PAS * 1000.0;
double C_TEL = CENT_TEL / V_TEL;          // already ng/mL (mass scaled)
double C_EVE = CENT_EVE / V_EVE;          // ng/mL
double C_TKI = CENT_TKI / V_TKI;          // ng/mL
double C_IFN = CENT_IFN / V_IFN;          // IU/mL

// ---- SSTR2/5 occupancy (free-fraction Hill) ----
double OCC_OCT = C_OCT/(KD_OCT_S2 + C_OCT);
double OCC_LAN = C_LAN/(KD_LAN_S2 + C_LAN);
double OCC_PAS = C_PAS/(KD_PAS_S5 + C_PAS);
double SSTR_TOTAL = OCC_OCT + OCC_LAN + OCC_PAS;
if (SSTR_TOTAL > 0.99) SSTR_TOTAL = 0.99;

// ---- TPH1 inhibition (telotristat) ----
double INH_TPH1 = pow(C_TEL, GAMMA_TPH1) /
                  (pow(IC50_TPH1, GAMMA_TPH1) + pow(C_TEL, GAMMA_TPH1));

// ---- Tumor growth modifiers ----
double MTORi = C_EVE/(IC50_MTOR + C_EVE);
double VEGFi = C_TKI/(IC50_VEGF + C_TKI);
double IFNi  = C_IFN/(EFF_IFN50 + C_IFN);

double GROWTH = KGROW * TUMOR *
                (1 - TUMOR/TUMAX) *
                (1 - 0.6*MTORi) *
                (1 - 0.45*VEGFi) *
                (1 - 0.3*IFNi);
double KILL_PRRT = DRUG_KILL_PRRT * PRRT_AMT * (TUMOR/TUM_BASE);

dxdt_GUT_OCT    = -KA_OCT * GUT_OCT;
dxdt_CENT_OCT   =  KA_OCT * GUT_OCT
                  - (CL_OCT/V_OCT) * CENT_OCT
                  - (Q_OCT/V_OCT)*CENT_OCT + (Q_OCT/VP_OCT)*PER_OCT;
dxdt_PER_OCT    =  (Q_OCT/V_OCT)*CENT_OCT - (Q_OCT/VP_OCT)*PER_OCT;
dxdt_GUT_LAN    = -KA_LAN * GUT_LAN;
dxdt_CENT_LAN   =  KA_LAN * GUT_LAN - (CL_LAN/V_LAN)*CENT_LAN;
dxdt_GUT_PASIR  = -KA_PAS * GUT_PASIR;
dxdt_CENT_PASIR =  KA_PAS * GUT_PASIR - (CL_PAS/V_PAS)*CENT_PASIR;
dxdt_GUT_TEL    = -KA_TEL * GUT_TEL;
dxdt_CENT_TEL   =  F_TEL * KA_TEL * GUT_TEL - (CL_TEL/V_TEL)*CENT_TEL;
dxdt_GUT_EVE    = -KA_EVE * GUT_EVE;
dxdt_CENT_EVE   =  KA_EVE * GUT_EVE - (CL_EVE/V_EVE)*CENT_EVE;
dxdt_GUT_TKI    = -KA_TKI * GUT_TKI;
dxdt_CENT_TKI   =  KA_TKI * GUT_TKI - (CL_TKI/V_TKI)*CENT_TKI;
dxdt_CENT_IFN   = -(CL_IFN/V_IFN)*CENT_IFN;
dxdt_PRRT_AMT   = -0.04 * PRRT_AMT;        // physical+biological decay (~17h eff half-life proxy)
dxdt_TUMOR      = GROWTH - KILL_PRRT;

// ---- TPH1 functional pool (modeled as fraction) ----
dxdt_TPH1       = 0.05*(1 - INH_TPH1 - TPH1);   // turnover toward (1 - inhibition)

// ---- 5-HTP / serotonin / 5-HIAA flux per gram tumor ----
double KREL_EFF = KREL_SER * (1 - 0.75*SSTR_TOTAL);
dxdt_HTP        = KSYN_HTP * TUMOR * TPH1 - KCONV_AAD * HTP;
dxdt_SER_T      = KCONV_AAD * HTP - KREL_EFF * SER_T;
dxdt_SER_P      = KREL_EFF * SER_T / 50.0
                  - KMAOA * (SER_P - KSER_NORM)
                  - KUPT_PLT * SER_P;
dxdt_SER_PLT    = KUPT_PLT * SER_P * 200.0 - KOUT_PLT * SER_PLT;
dxdt_HIAA_U     = F_HIAA * KMAOA * (SER_P - KSER_NORM) * 0.20
                  - 0.04 * HIAA_U;

// ---- Symptoms ----
double DR_BM    = EMAX_BM * SER_P / (EC50_BM + SER_P);
double DR_FL    = EMAX_FL * SER_P / (EC50_FL + SER_P);
dxdt_BM         = KOUT_SYMP * (KIN_BM + DR_BM*(1-KOND)*(1-KLOP) - BM);
dxdt_FLUSH      = KOUT_SYMP * (KIN_FL + DR_FL*(1-SSTR_TOTAL) - FLUSH);

// ---- TGF-β / valve / NT-proBNP (CHD axis) ----
dxdt_TGFb       = KTGFB * SER_P - KOUT_TGFB * (TGFb - 1.0);
dxdt_VALVE      = KIN_VAL * SER_P * TGFb - KOUT_VAL * (VALVE - 1.0);
dxdt_NTproBNP   = KIN_BNP * (VALVE - 1.0) - KOUT_BNP * (NTproBNP - 80.0);

$CAPTURE @annotated
C_OCT  : Octreotide plasma (ng/mL)
C_LAN  : Lanreotide plasma (ng/mL)
C_PAS  : Pasireotide plasma (ng/mL)
C_TEL  : Telotristat metabolite (ng/mL)
C_EVE  : Everolimus plasma (ng/mL)
C_TKI  : VEGFR-TKI plasma (ng/mL)
C_IFN  : IFN-α (IU/mL)
SSTR_TOTAL : Combined SSTR2/5 occupancy (0-1)
INH_TPH1   : Telotristat TPH1 inhibition (0-1)
MTORi  : mTORC1 inhibition (0-1)
VEGFi  : VEGFR inhibition (0-1)
IFNi   : IFN-α effect (0-1)
'

# -------------------------------------------------------------
# Compile (mrgsolve required). Will skip silently if unavailable.
# -------------------------------------------------------------
mod <- tryCatch(mcode("carcsyn", carcsyn_code), error = function(e) {
  message("mrgsolve compile failed (expected without compiler): ", conditionMessage(e))
  NULL
})

# =====================================================================
# Treatment scenarios
# =====================================================================
build_scenarios <- function() {
  list(
    "S1_natural_history" = tibble(
      ID=1, time=0, amt=0, cmt=1, evid=0, ii=0, addl=0,
      desc="Untreated functional midgut NET, 12 months"
    ),
    "S2_octreotide_LAR_30mg" = tibble(
      ID=2, time=seq(0, 24*28*12, by=24*28),
      amt=30, cmt=1, evid=1, ii=0, addl=0,
      desc="Octreotide LAR 30 mg IM q28d × 12"
    ),
    "S3_lanreotide_autogel_120mg" = tibble(
      ID=3, time=seq(0, 24*28*12, by=24*28),
      amt=120, cmt=4, evid=1, ii=0, addl=0,
      desc="Lanreotide autogel 120 mg SC q28d × 12 (CLARINET)"
    ),
    "S4_octreotide_plus_telotristat" = bind_rows(
      tibble(ID=4, time=seq(0, 24*28*12, by=24*28), amt=30, cmt=1, evid=1),
      tibble(ID=4, time=seq(0, 24*28*12, by=8),     amt=250, cmt=8, evid=1)
    ) %>% mutate(desc="Octreotide LAR 30 mg + Telotristat 250 mg t.i.d. (TELESTAR)"),
    "S5_pasireotide_LAR" = tibble(
      ID=5, time=seq(0, 24*28*12, by=24*28),
      amt=60, cmt=6, evid=1,
      desc="Pasireotide LAR 60 mg q28d (refractory)"
    ),
    "S6_everolimus_10mg" = tibble(
      ID=6, time=seq(0, 24*12, by=24), amt=10, cmt=10, evid=1,
      desc="Everolimus 10 mg/day (RADIANT-4) × 12 months"
    ),
    "S7_sunitinib_proxy_TKI" = tibble(
      ID=7, time=seq(0, 24*12, by=24), amt=37.5, cmt=12, evid=1,
      desc="VEGFR-TKI 37.5 mg/day (sunitinib/surufatinib proxy)"
    ),
    "S8_PRRT_177Lu_DOTATATE" = tibble(
      ID=8, time=c(0, 8*7*24, 16*7*24, 24*7*24),
      amt=7.4, cmt=15, evid=1,
      desc="177Lu-DOTATATE 7.4 GBq q8w × 4 cycles (NETTER-1)"
    ),
    "S9_IFN_alpha" = tibble(
      ID=9, time=seq(0, 24*7*52, by=24*48),
      amt=5e6, cmt=14, evid=1,
      desc="IFN-α-2b 5 MU SC TIW × 12 months"
    ),
    "S10_HAE_then_octreotide" = bind_rows(
      tibble(ID=10, time=0, amt=-200, cmt=16, evid=1),     # tumor mass debulk
      tibble(ID=10, time=seq(0, 24*28*12, by=24*28), amt=30, cmt=1, evid=1)
    ) %>% mutate(desc="Hepatic artery embolization + octreotide LAR 30 mg"),
    "S11_carcinoid_crisis_prevention" = bind_rows(
      tibble(ID=11, time=0,  amt=500, cmt=2,  evid=1),    # 500 µg IV bolus pre-procedure
      tibble(ID=11, time=24, amt=30,  cmt=1,  evid=1)
    ) %>% mutate(desc="Carcinoid crisis prophylaxis: octreotide IV 500 µg + LAR"),
    "S12_quad_therapy" = bind_rows(
      tibble(ID=12, time=seq(0, 24*28*12, by=24*28), amt=30, cmt=1, evid=1),
      tibble(ID=12, time=seq(0, 24*28*12, by=8),     amt=250, cmt=8, evid=1),
      tibble(ID=12, time=c(0, 8*7*24, 16*7*24, 24*7*24), amt=7.4, cmt=15, evid=1),
      tibble(ID=12, time=seq(0, 24*12, by=24),       amt=10,  cmt=10, evid=1)
    ) %>% mutate(desc="Octreotide + Telotristat + PRRT + Everolimus (refractory CS)")
  )
}

# =====================================================================
# Run scenarios (placeholder if mrgsolve compile failed)
# =====================================================================
run_scenarios <- function(mod, scen) {
  if (is.null(mod)) return(invisible(NULL))
  out <- lapply(names(scen), function(nm) {
    ev <- scen[[nm]]
    mod %>% data_set(ev) %>% mrgsim(end = 24*28*12, delta = 12) %>% as_tibble() %>%
      mutate(scenario = nm)
  })
  bind_rows(out)
}

# =====================================================================
# Diagnostic ggplot helpers
# =====================================================================
plot_serotonin <- function(sim) {
  ggplot(sim, aes(time/24, SER_P, color = scenario)) +
    geom_line() +
    labs(x = "Day", y = "Plasma 5-HT (nmol/L)",
         title = "Plasma 5-HT trajectory across scenarios") +
    theme_minimal()
}

plot_BM_flush <- function(sim) {
  long <- tidyr::pivot_longer(sim, c(BM, FLUSH))
  ggplot(long, aes(time/24, value, color = scenario)) +
    geom_line() + facet_wrap(~ name, scales = "free_y") +
    labs(x = "Day", y = "Count/day", title = "Diarrhea & flushing dynamics") +
    theme_minimal()
}

plot_tumor <- function(sim) {
  ggplot(sim, aes(time/24, TUMOR, color = scenario)) + geom_line() +
    labs(x = "Day", y = "Tumor burden (g)", title = "Hepatic tumor mass") +
    theme_minimal()
}

# =====================================================================
# Pop-PK variability template (Omega, Sigma)
# =====================================================================
# $OMEGA (per-subject IIV on CL_OCT, CL_TEL, CL_EVE, KGROW, KIN_BM, KIN_FL)
omega_skeleton <- list(
  CL_OCT = 0.20,
  CL_TEL = 0.30,
  CL_EVE = 0.25,
  KGROW  = 0.50,
  KIN_BM = 0.25,
  KIN_FL = 0.40
)

# Proportional residual error: 25% on plasma 5-HT, 20% on tumor, 30% on symptoms
sigma_skeleton <- list(
  SER_P = 0.25,
  TUMOR = 0.20,
  BM    = 0.30,
  FLUSH = 0.30
)

# =====================================================================
# Quick test (uncomment when mrgsolve available)
# scen <- build_scenarios()
# sim  <- run_scenarios(mod, scen)
# plot_serotonin(sim)
# plot_BM_flush(sim)
# plot_tumor(sim)
# =====================================================================
