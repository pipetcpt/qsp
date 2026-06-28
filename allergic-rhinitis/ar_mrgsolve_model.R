################################################################################
# Allergic Rhinitis (AR) – Quantitative Systems Pharmacology (QSP) Model
# mrgsolve implementation  |  R ≥ 4.2, mrgsolve ≥ 1.0
#
# Model scope:
#   • IgE-mast cell axis (sensitization + allergen challenge)
#   • Th2 cytokine cascade (IL-4, IL-5, IL-13)
#   • Eosinophil dynamics (blood & nasal tissue)
#   • Histamine & cysteinyl-leukotriene release
#   • Drug PK/PD:
#       – Cetirizine       10 mg q.d. (oral, 1-comp)
#       – Fluticasone FP   200 μg/d   (intranasal, 1-comp)
#       – Montelukast      10 mg q.d. (oral, 1-comp)
#       – Omalizumab       300 mg q4w (SC, 2-comp)
#   • Clinical endpoints: TNSS, eosinophilia, specific IgE
#
# Calibration targets (literature):
#   • Cetirizine H1-RO ≥80% at steady state  [Yanai 1995 JACI]
#   • FP reduces TNSS ~35-40% vs placebo     [Meltzer 2005 JACI]
#   • Montelukast reduces TNSS ~20-25%       [Philip 2002 JACI]
#   • Omalizumab reduces free IgE >95%       [Fahy 1997 AJRCCM]
#   • Nasal eosinophilia reduced ~50% by FP  [Holgate 2003 Allergy]
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ============================================================
# Model code block
# ============================================================

AR_model_code <- '
$PROB
  Allergic Rhinitis QSP Model
  IgE/Mast Cell · Th2 Cytokines · Eosinophils · PK/PD
  Cetirizine | Fluticasone | Montelukast | Omalizumab

$PARAM @annotated
  // --- Allergen ---
  ALLERGEN_SS  : 1.0   : Steady-state nasal allergen load (AU/mL, normalized)
  K_ALLERGEN   : 0.1   : Allergen clearance rate (1/h)

  // --- IgE dynamics ---
  KSY_IGE      : 0.002 : IgE synthesis rate (AU/h, baseline)
  KDEG_IGE     : 0.005 : Free IgE degradation rate (1/h, t1/2~138h)
  K_BIND_MAST  : 0.1   : IgE binding to mast cell FcεRI (1/h)
  K_OFF_MAST   : 0.001 : IgE dissociation from FcεRI (1/h)
  MAST_TOTAL   : 1.0   : Total mast cell FcεRI capacity (AU, normalized)

  // --- Mast cell activation & degranulation ---
  EC50_CROSS   : 0.5   : Allergen EC50 for crosslinking IgE-Ag (AU)
  HILL_CROSS   : 2.0   : Hill coefficient for crosslinking
  KDEG_HIST    : 2.0   : Histamine degradation (1/h, DAO/HNMT)
  KDEG_LT      : 1.0   : CysLT degradation (1/h)
  KHIST_PROD   : 5.0   : Max histamine production per unit mast activation
  KLT_PROD     : 2.0   : Max CysLT production per unit mast activation
  KMAST_REC    : 0.05  : Mast cell recharging rate after degranulation (1/h)

  // --- Th2 cytokines (IL-4, IL-5, IL-13) ---
  KSY_TH2      : 0.05  : Th2 cell activation rate driven by IgE/allergen
  KDEG_TH2     : 0.1   : Th2 cell decay (1/h)
  KSY_IL4      : 1.0   : IL-4 synthesis per Th2 unit (pg/mL/h)
  KSY_IL5      : 0.8   : IL-5 synthesis per Th2 unit
  KSY_IL13     : 1.2   : IL-13 synthesis per Th2 unit
  KDEG_IL4     : 0.5   : IL-4 clearance (1/h)
  KDEG_IL5     : 0.3   : IL-5 clearance (1/h)
  KDEG_IL13    : 0.4   : IL-13 clearance (1/h)
  TH2_BASE     : 0.5   : Baseline Th2 activity (atopic, normalized)

  // --- Eosinophil dynamics ---
  EOS_BLOOD_0  : 300   : Baseline blood eosinophil (cells/μL)
  KEO_PROD     : 0.05  : Eosinophil production stimulated by IL-5 (cells/μL/h)
  KEO_SURV     : 0.02  : IL-5-mediated eosinophil survival boost (1/h)
  KEO_DEATH    : 0.08  : Eosinophil apoptosis rate (1/h, t1/2~8.7h)
  KEO_MIGRATE  : 0.01  : Blood→nasal tissue migration (1/h) per chemokine unit
  KCHEMOKINE   : 2.0   : Eotaxin/CCL11 driven chemotaxis coefficient
  KEO_TISSUE0  : 10.0  : Baseline nasal tissue eosinophil (cells/μL tissue)
  KEO_TIS_DEATH: 0.05  : Nasal tissue eosinophil apoptosis (1/h)

  // --- Symptom model ---
  HIST_EC50    : 1.0   : Histamine EC50 for symptom generation
  LT_EC50      : 0.5   : CysLT EC50 for congestion
  EOS_EC50     : 50.0  : Nasal eosinophil EC50 for late symptoms
  SNEEZE_MAX   : 3.0   : Maximum sneezing score
  RHINO_MAX    : 3.0   : Maximum rhinorrhea score
  CONG_MAX     : 3.0   : Maximum congestion score
  PRUR_MAX     : 3.0   : Maximum pruritus score

  // =====================================================
  // PK Parameters – Cetirizine (oral, 1-comp)
  // Tmax ~1h, t1/2 ~10h, Vd 0.56 L/kg
  // =====================================================
  KA_CETI      : 0.9   : Cetirizine absorption rate (1/h)
  CL_CETI      : 7.0   : Cetirizine clearance (L/h)
  VD_CETI      : 70.0  : Volume of distribution (L)
  F_CETI       : 0.70  : Oral bioavailability cetirizine
  H1_IC50_CETI : 15.0  : Cetirizine IC50 for H1R occupancy (ng/mL)
  H1_HILL      : 1.0   : Hill coeff H1 block

  // =====================================================
  // PK Parameters – Fluticasone Propionate (intranasal)
  // Local mucosal concentration, systemic exposure <2%
  // =====================================================
  KA_FP        : 1.5   : FP nasal mucosal absorption (1/h)
  CL_FP_LOCAL  : 2.0   : FP local clearance from nasal mucosa (1/h)
  VD_FP_LOCAL  : 1.0   : FP mucosal volume (L, local compartment)
  GR_IC50_FP   : 0.5   : FP IC50 for GR-mediated cytokine suppression (nM)
  GR_HILL      : 1.2   : Hill coeff for GR inhibition
  FP_DOSE_BIOCONV: 200.0: Intranasal dose converted to nM (200 μg ~ 450 nM local)

  // =====================================================
  // PK Parameters – Montelukast (oral, 1-comp)
  // Tmax 3-4h, t1/2 5.5h, BA ~64%
  // =====================================================
  KA_MLKT      : 0.5   : Montelukast absorption (1/h)
  CL_MLKT      : 45.0  : Montelukast clearance (L/h)
  VD_MLKT      : 10.0  : Volume of distribution (L)
  F_MLKT       : 0.64  : Oral bioavailability montelukast
  CYSLTR1_IC50 : 2.0   : Montelukast IC50 CysLT1 (ng/mL)

  // =====================================================
  // PK Parameters – Omalizumab (SC, 2-comp)
  // t1/2 ~26d, ka ~0.004/h, Vss ~78 mL/kg
  // =====================================================
  KA_OMA       : 0.004 : Omalizumab SC absorption (1/h)
  CL_OMA       : 0.14  : Omalizumab central clearance (mL/h)
  VC_OMA       : 3140  : Omalizumab central volume (mL, ~45 mL/kg × 70kg)
  VP_OMA       : 2360  : Omalizumab peripheral volume (mL)
  Q_OMA        : 0.40  : Intercompartmental clearance (mL/h)
  KON_IGE      : 1.0   : Omalizumab-IgE association rate
  KOFF_IGE     : 0.0001: Omalizumab-IgE dissociation (KD ~0.1 nM)

$INIT @annotated
  // Allergen
  AG       = 0.0   : Nasal allergen (AU)

  // IgE compartments
  IGE_FREE = 50.0  : Free IgE (IU/mL equivalent)
  IGE_MAST = 0.5   : Mast cell-bound IgE (normalized)

  // Mast cell state
  MAST_ACT = 0.0   : Mast cell activation level (0-1)
  MAST_CHG = 1.0   : Mast cell granule charge (0-1)

  // Mediators
  HISTAMINE = 0.0  : Nasal histamine (normalized AU)
  CYS_LT   = 0.0   : Cysteinyl leukotrienes (AU)

  // Th2 / cytokines
  TH2      = 0.5   : Th2 cell activity (normalized, atopic baseline)
  IL4      = 5.0   : IL-4 (pg/mL)
  IL5      = 3.0   : IL-5 (pg/mL)
  IL13     = 8.0   : IL-13 (pg/mL)

  // Eosinophils
  EOS_B    = 300.0 : Blood eosinophil (cells/μL)
  EOS_N    = 10.0  : Nasal tissue eosinophil (cells/μL tissue)

  // Cetirizine PK (depot + central)
  CETI_D   = 0.0   : Cetirizine depot (mg)
  CETI_C   = 0.0   : Cetirizine central (mg)

  // Fluticasone PK (nasal local)
  FP_LOC   = 0.0   : Fluticasone local nasal conc (nM-normalized)

  // Montelukast PK
  MLKT_D   = 0.0   : Montelukast depot (mg)
  MLKT_C   = 0.0   : Montelukast central (mg)

  // Omalizumab PK (depot + 2-comp + IgE complex)
  OMA_D    = 0.0   : Omalizumab SC depot (mg)
  OMA_C    = 0.0   : Omalizumab central (mg/mL × volume → mg)
  OMA_P    = 0.0   : Omalizumab peripheral (mg)
  OMA_IGE  = 0.0   : Omalizumab-IgE complex (IU-equivalents)

$MAIN
  // Allergen challenge: set AG from input or steady-state
  // Dose event (compartment 1) delivers allergen pulse

$ODE
  // ==============================================================
  // 1. Allergen kinetics
  // ==============================================================
  double Ag = AG;
  dxdt_AG = -K_ALLERGEN * AG;  // cleared by mucociliary; replenished via dose event

  // ==============================================================
  // 2. IgE dynamics
  // ==============================================================
  double degen_free = KDEG_IGE * IGE_FREE;
  double bind_rate  = K_BIND_MAST * IGE_FREE * (MAST_TOTAL - IGE_MAST);
  double unbind_rate= K_OFF_MAST  * IGE_MAST;

  // Omalizumab complex: removes free IgE
  double oma_cp = OMA_C / VC_OMA;  // mg/mL → concentration
  double ige_capture = KON_IGE * oma_cp * IGE_FREE - KOFF_IGE * OMA_IGE;

  dxdt_IGE_FREE = KSY_IGE * 200.0 - degen_free - bind_rate + unbind_rate - ige_capture;
  dxdt_IGE_MAST = bind_rate - unbind_rate;
  dxdt_OMA_IGE  = ige_capture;

  // ==============================================================
  // 3. Mast cell crosslinking & activation
  // ==============================================================
  double crosslink_frac = pow(Ag, HILL_CROSS) / (pow(EC50_CROSS, HILL_CROSS) + pow(Ag, HILL_CROSS));
  double mast_trigger   = crosslink_frac * IGE_MAST * MAST_CHG;  // depends on bound IgE and granule charge

  dxdt_MAST_ACT = mast_trigger - 0.5 * MAST_ACT;
  dxdt_MAST_CHG = KMAST_REC * (1.0 - MAST_CHG) - mast_trigger * MAST_CHG;

  // ==============================================================
  // 4. Histamine & CysLT release/degradation
  // ==============================================================
  dxdt_HISTAMINE = KHIST_PROD * MAST_ACT - KDEG_HIST * HISTAMINE;
  dxdt_CYS_LT   = KLT_PROD   * MAST_ACT - KDEG_LT   * CYS_LT;

  // ==============================================================
  // 5. Th2 / cytokine cascade
  // ==============================================================
  double th2_drive = TH2_BASE + KSY_TH2 * Ag * IGE_FREE / (1.0 + IGE_FREE);
  dxdt_TH2  = th2_drive - KDEG_TH2 * TH2;
  dxdt_IL4  = KSY_IL4  * TH2 - KDEG_IL4  * IL4;
  dxdt_IL5  = KSY_IL5  * TH2 - KDEG_IL5  * IL5;
  dxdt_IL13 = KSY_IL13 * TH2 - KDEG_IL13 * IL13;

  // ==============================================================
  // 6. Eosinophil dynamics
  // ==============================================================
  double eos_prod    = KEO_PROD * IL5 + KEO_SURV * EOS_B * IL5 / (1.0 + IL5);
  double eos_death_b = KEO_DEATH * EOS_B;
  double chemokine   = KCHEMOKINE * IL5 * IL13;
  double eos_migrate = KEO_MIGRATE * chemokine * EOS_B;
  double eos_death_n = KEO_TIS_DEATH * EOS_N;

  dxdt_EOS_B = eos_prod - eos_death_b - eos_migrate;
  dxdt_EOS_N = eos_migrate - eos_death_n;

  // ==============================================================
  // 7. Drug PK
  // ==============================================================

  // 7a. Cetirizine
  double ke_ceti = CL_CETI / VD_CETI;
  dxdt_CETI_D = -KA_CETI * CETI_D;
  dxdt_CETI_C =  KA_CETI * F_CETI * CETI_D - ke_ceti * CETI_C;
  double ceti_cp_ngmL = CETI_C / VD_CETI * 1000.0;  // mg/L → ng/mL

  // 7b. Fluticasone Propionate (intranasal local)
  dxdt_FP_LOC = KA_FP * FP_DOSE_BIOCONV / VD_FP_LOCAL - CL_FP_LOCAL * FP_LOC;
  // FP_LOC replenished by dose events (compartment 3 delivers dose in μg→nM)

  // 7c. Montelukast
  double ke_mlkt = CL_MLKT / VD_MLKT;
  dxdt_MLKT_D = -KA_MLKT * MLKT_D;
  dxdt_MLKT_C =  KA_MLKT * F_MLKT * MLKT_D - ke_mlkt * MLKT_C;
  double mlkt_cp_ngmL = MLKT_C / VD_MLKT * 1000.0;

  // 7d. Omalizumab (2-comp)
  dxdt_OMA_D = -KA_OMA * OMA_D;
  dxdt_OMA_C =  KA_OMA * OMA_D - (CL_OMA/VC_OMA + Q_OMA/VC_OMA) * OMA_C + Q_OMA/VP_OMA * OMA_P;
  dxdt_OMA_P =  Q_OMA/VC_OMA * OMA_C - Q_OMA/VP_OMA * OMA_P;

  // ==============================================================
  // 8. Drug pharmacodynamic effects (inhibitory)
  // ==============================================================

  // Cetirizine: H1R occupancy → attenuate histamine effect
  double h1ro = pow(ceti_cp_ngmL, H1_HILL) / (pow(H1_IC50_CETI, H1_HILL) + pow(ceti_cp_ngmL, H1_HILL));

  // Fluticasone: GR occupancy → suppress cytokines (IL-4, IL-5, IL-13) & eosinophil
  double fp_nM = FP_LOC;
  double gr_occ = pow(fp_nM, GR_HILL) / (pow(GR_IC50_FP, GR_HILL) + pow(fp_nM, GR_HILL));

  // Apply FP effect on cytokines (transrepression)
  dxdt_IL4  -= gr_occ * 0.6 * IL4;
  dxdt_IL5  -= gr_occ * 0.7 * IL5;
  dxdt_IL13 -= gr_occ * 0.6 * IL13;
  dxdt_EOS_N -= gr_occ * 0.5 * EOS_N;

  // Montelukast: CysLT1 block → reduce LT-driven congestion
  double cysltr1_inh = pow(mlkt_cp_ngmL, 1.0) / (CYSLTR1_IC50 + mlkt_cp_ngmL);

  // Modify mediator-driven outputs via drug PD variables (used in $TABLE)
  // (PD captured as secondary vars below)

$TABLE
  // ==============================================================
  // Secondary PK variables
  // ==============================================================
  double CETI_CP     = CETI_C / VD_CETI * 1000.0;       // ng/mL
  double FP_LOCAL_NM = FP_LOC;                           // nM
  double MLKT_CP     = MLKT_C / VD_MLKT * 1000.0;       // ng/mL
  double OMA_CP      = OMA_C  / (VC_OMA/1000.0);        // μg/mL

  // H1 receptor occupancy
  double H1_RO = pow(CETI_CP, H1_HILL) / (pow(H1_IC50_CETI, H1_HILL) + pow(CETI_CP, H1_HILL)) * 100.0;  // %

  // GR occupancy (FP)
  double GR_OCC_FP = pow(FP_LOCAL_NM, GR_HILL) / (pow(GR_IC50_FP, GR_HILL) + pow(FP_LOCAL_NM, GR_HILL)) * 100.0; // %

  // CysLT1 receptor inhibition (montelukast)
  double CYSLTR1_INH = MLKT_CP / (CYSLTR1_IC50 + MLKT_CP) * 100.0;  // %

  // Free IgE reduction (omalizumab; vs baseline 50)
  double IGE_REDUCTION_PCT = (1.0 - IGE_FREE/50.0) * 100.0;

  // ==============================================================
  // Drug-modified mediator outputs
  // ==============================================================
  double H1_EFF_HIST = HISTAMINE * (1.0 - H1_RO/100.0);     // Histamine × (1-H1RO)
  double LT_EFF      = CYS_LT   * (1.0 - CYSLTR1_INH/100.0); // LT × (1-CysLT1 block)

  // ==============================================================
  // Symptom scores (0-3 scale, Emax model)
  // ==============================================================
  // Sneezing: mainly histamine + sensory nerve
  double SNEEZE = SNEEZE_MAX * H1_EFF_HIST / (HIST_EC50 + H1_EFF_HIST);

  // Rhinorrhea: histamine + LTs + glandular
  double RHINORRHEA = RHINO_MAX * (0.6 * H1_EFF_HIST + 0.4 * LT_EFF) /
                      (HIST_EC50 + 0.6 * H1_EFF_HIST + 0.4 * LT_EFF);

  // Congestion: mainly LTs + PGD2 (approximated by CYS_LT) + eosinophil
  double EOS_CONG = EOS_N / (EOS_EC50 + EOS_N);
  double CONGESTION = CONG_MAX * (0.5 * LT_EFF / (LT_EC50 + LT_EFF) + 0.3 * EOS_CONG + 0.2);

  // Pruritus: histamine + Th2 cytokines (IL-31 proxy → IL-13)
  double IL13_NORM = IL13 / (IL13 + 8.0);
  double PRURITUS = PRUR_MAX * (0.7 * H1_EFF_HIST / (HIST_EC50 + H1_EFF_HIST) + 0.3 * IL13_NORM);

  // TNSS (0-12)
  double TNSS = SNEEZE + RHINORRHEA + CONGESTION + PRURITUS;

  // ==============================================================
  // Biomarker endpoints
  // ==============================================================
  double FREE_IGE_IU  = IGE_FREE;          // IU/mL
  double BLOOD_EOS_UL = EOS_B;             // cells/μL
  double NASAL_EOS_UL = EOS_N;             // cells/μL tissue
  double TRYPTASE_MCG = MAST_ACT * 15.0;  // serum tryptase proxy (ng/mL)

  capture CETI_CP MLKT_CP OMA_CP FP_LOCAL_NM H1_RO GR_OCC_FP CYSLTR1_INH
  capture SNEEZE RHINORRHEA CONGESTION PRURITUS TNSS
  capture FREE_IGE_IU BLOOD_EOS_UL NASAL_EOS_UL TRYPTASE_MCG
  capture IL4 IL5 IL13 HISTAMINE CYS_LT MAST_ACT EOS_N
  capture H1_EFF_HIST LT_EFF IGE_REDUCTION_PCT

'

# ============================================================
# Compile model
# ============================================================
AR_model <- mcode("AllergicRhinitisQSP", AR_model_code)

# ============================================================
# Helper: build dosing regimens
# ============================================================
build_doses <- function(
    cetirizine  = FALSE,  # 10 mg QD oral
    fluticasone = FALSE,  # 200 μg/day intranasal (split as 100 μg BID)
    montelukast = FALSE,  # 10 mg QD oral
    omalizumab  = FALSE,  # 300 mg q4w SC
    allergen_pulse = FALSE, # allergen challenge on day 28
    sim_duration_d = 84
) {
  ev_list <- list()

  if (cetirizine)  ev_list[["ceti"]]  <- ev(cmt = "CETI_D",  amt = 10,   time = 0, ii = 24,   addl = sim_duration_d - 1)
  if (fluticasone) ev_list[["fp"]]    <- ev(cmt = "FP_LOC",  amt = 450,  time = 0, ii = 24,   addl = sim_duration_d - 1)
  if (montelukast) ev_list[["mlkt"]]  <- ev(cmt = "MLKT_D",  amt = 10,   time = 0, ii = 24,   addl = sim_duration_d - 1)
  if (omalizumab)  ev_list[["oma"]]   <- ev(cmt = "OMA_D",   amt = 300,  time = 0, ii = 28*24, addl = 2)

  if (allergen_pulse) {
    ev_list[["ag"]] <- ev(cmt = "AG", amt = 5.0, time = 28*24)
  }

  if (length(ev_list) == 0) return(ev(cmt = "AG", amt = 5.0, time = 28*24))
  Reduce(c, ev_list)
}

# ============================================================
# Scenario definitions
# ============================================================
scenarios <- list(
  list(
    name        = "1. Natural History\n(Allergen only, no Tx)",
    cetirizine  = FALSE, fluticasone = FALSE,
    montelukast = FALSE, omalizumab  = FALSE,
    allergen_pulse = TRUE, color = "#E53935"
  ),
  list(
    name        = "2. Cetirizine 10 mg QD",
    cetirizine  = TRUE,  fluticasone = FALSE,
    montelukast = FALSE, omalizumab  = FALSE,
    allergen_pulse = TRUE, color = "#1E88E5"
  ),
  list(
    name        = "3. Fluticasone FP 200 μg/d",
    cetirizine  = FALSE, fluticasone = TRUE,
    montelukast = FALSE, omalizumab  = FALSE,
    allergen_pulse = TRUE, color = "#43A047"
  ),
  list(
    name        = "4. Montelukast 10 mg QD",
    cetirizine  = FALSE, fluticasone = FALSE,
    montelukast = TRUE,  omalizumab  = FALSE,
    allergen_pulse = TRUE, color = "#FB8C00"
  ),
  list(
    name        = "5. Cetirizine + Fluticasone\n(Combination)",
    cetirizine  = TRUE,  fluticasone = TRUE,
    montelukast = FALSE, omalizumab  = FALSE,
    allergen_pulse = TRUE, color = "#8E24AA"
  ),
  list(
    name        = "6. Omalizumab 300 mg q4w\n(Anti-IgE)",
    cetirizine  = FALSE, fluticasone = FALSE,
    montelukast = FALSE, omalizumab  = TRUE,
    allergen_pulse = TRUE, color = "#00897B"
  ),
  list(
    name        = "7. Triple Therapy\n(Ceti + FP + MLKT)",
    cetirizine  = TRUE,  fluticasone = TRUE,
    montelukast = TRUE,  omalizumab  = FALSE,
    allergen_pulse = TRUE, color = "#6D4C41"
  )
)

# ============================================================
# Simulation function
# ============================================================
run_scenario <- function(sc) {
  evs <- build_doses(
    cetirizine  = sc$cetirizine,
    fluticasone = sc$fluticasone,
    montelukast = sc$montelukast,
    omalizumab  = sc$omalizumab,
    allergen_pulse = sc$allergen_pulse,
    sim_duration_d = 84
  )
  AR_model %>%
    ev(evs) %>%
    mrgsim(end = 84 * 24, delta = 1) %>%
    as.data.frame() %>%
    mutate(scenario = sc$name, color = sc$color, time_d = time / 24)
}

# Run all scenarios
message("Running ", length(scenarios), " scenarios...")
results <- bind_rows(lapply(scenarios, run_scenario))

# ============================================================
# Visualization functions
# ============================================================

# 1. TNSS over time
plot_tnss <- function(data) {
  ggplot(data, aes(x = time_d, y = TNSS, color = scenario)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = setNames(
      unique(data$color), unique(data$scenario)
    )) +
    geom_vline(xintercept = 28, linetype = "dashed", color = "grey40", alpha = 0.7) +
    annotate("text", x = 28.5, y = 11, label = "Allergen\nChallenge", size = 3, hjust = 0) +
    labs(
      title = "Total Nasal Symptom Score (TNSS) – All Scenarios",
      subtitle = "0-12 scale; allergen challenge on Day 28",
      x = "Time (days)", y = "TNSS (0-12)",
      color = "Scenario"
    ) +
    ylim(0, 12) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right", legend.text = element_text(size = 8))
}

# 2. Individual symptom scores
plot_symptoms <- function(data) {
  sym_df <- data %>%
    select(time_d, scenario, color, SNEEZE, RHINORRHEA, CONGESTION, PRURITUS) %>%
    pivot_longer(c(SNEEZE, RHINORRHEA, CONGESTION, PRURITUS), names_to = "symptom", values_to = "score")

  ggplot(sym_df, aes(x = time_d, y = score, color = scenario)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~symptom, scales = "free_y", ncol = 2) +
    scale_color_manual(values = setNames(
      unique(data$color), unique(data$scenario)
    )) +
    geom_vline(xintercept = 28, linetype = "dashed", color = "grey40", alpha = 0.5) +
    labs(
      title = "Individual Symptom Scores",
      x = "Time (days)", y = "Score (0-3)",
      color = "Scenario"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.text = element_text(size = 7))
}

# 3. PK – H1R occupancy (cetirizine)
plot_pk_ceti <- function(data) {
  df <- data %>% filter(grepl("Cetirizine|Combination|Triple", scenario))
  ggplot(df, aes(x = time_d, y = H1_RO, color = scenario)) +
    geom_line(linewidth = 0.9) +
    labs(
      title = "H1 Receptor Occupancy – Cetirizine",
      subtitle = "Steady-state target ≥80% for clinical efficacy",
      x = "Time (days)", y = "H1 Receptor Occupancy (%)",
      color = "Scenario"
    ) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "blue") +
    annotate("text", x = 5, y = 82, label = "80% target", color = "blue", size = 3) +
    ylim(0, 100) +
    theme_bw(base_size = 12)
}

# 4. Biomarkers – Free IgE & Omalizumab
plot_ige <- function(data) {
  df_oma <- data %>% filter(grepl("Omalizumab", scenario) | grepl("Natural", scenario))
  p1 <- ggplot(df_oma, aes(x = time_d, y = FREE_IGE_IU, color = scenario)) +
    geom_line(linewidth = 0.9) +
    labs(title = "Free IgE (IU/mL)", x = "Time (days)", y = "Free IgE (IU/mL)") +
    theme_bw(base_size = 11)
  print(p1)
}

# 5. Eosinophil dynamics
plot_eos <- function(data) {
  df <- data %>%
    select(time_d, scenario, color, BLOOD_EOS_UL, NASAL_EOS_UL) %>%
    pivot_longer(c(BLOOD_EOS_UL, NASAL_EOS_UL), names_to = "compartment", values_to = "eos")
  ggplot(df, aes(x = time_d, y = eos, color = scenario)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~compartment, scales = "free_y") +
    scale_color_manual(values = setNames(unique(data$color), unique(data$scenario))) +
    geom_vline(xintercept = 28, linetype = "dashed", color = "grey40", alpha = 0.5) +
    labs(
      title = "Eosinophil Dynamics",
      x = "Time (days)", y = "Eosinophils (cells/μL)",
      color = "Scenario"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.text = element_text(size = 7))
}

# 6. Summary table: mean TNSS at peak (day 28-35) & end (day 77-84)
summary_table <- function(data) {
  data %>%
    group_by(scenario) %>%
    summarise(
      TNSS_baseline    = mean(TNSS[time_d < 28], na.rm = TRUE),
      TNSS_peak        = max(TNSS[time_d >= 28 & time_d <= 35], na.rm = TRUE),
      TNSS_wk12        = mean(TNSS[time_d >= 77 & time_d <= 84], na.rm = TRUE),
      Free_IgE_wk12    = mean(FREE_IGE_IU[time_d >= 77 & time_d <= 84], na.rm = TRUE),
      Blood_Eos_wk12   = mean(BLOOD_EOS_UL[time_d >= 77 & time_d <= 84], na.rm = TRUE),
      Nasal_Eos_wk12   = mean(NASAL_EOS_UL[time_d >= 77 & time_d <= 84], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      TNSS_pct_chg = round((TNSS_wk12 - TNSS_peak[1]) / TNSS_peak[1] * 100, 1)
    )
}

# ============================================================
# Generate plots
# ============================================================
p_tnss <- plot_tnss(results)
p_sym  <- plot_symptoms(results)
p_eos  <- plot_eos(results)

print(p_tnss)
print(p_sym)
print(p_eos)

tbl <- summary_table(results)
print(tbl)

message("\nScenario summary:")
message("  TNSS at Week 12 (days 77-84):")
for (i in seq_len(nrow(tbl))) {
  message(sprintf("    %-45s  %.2f", tbl$scenario[i], tbl$TNSS_wk12[i]))
}
