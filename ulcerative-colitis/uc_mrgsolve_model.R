# ============================================================
# Ulcerative Colitis (UC) - Quantitative Systems Pharmacology
# mrgsolve ODE Model
# ============================================================
# References:
#   Fasanmade AA et al. Clin Pharmacokinet 2010 (Infliximab PK)
#   Rosario M et al. Clin Pharmacokinet 2015 (Vedolizumab PK)
#   Dowty ME et al. Drug Metab Dispos 2014 (Tofacitinib PK)
#   Xu Z et al. J Crohns Colitis 2017 (Ustekinumab UNIFI PK)
#   D'Ambrosio D et al. Br J Pharmacol 2021 (Ozanimod PK)
# ============================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# ============================================================
# MODEL DEFINITION
# ============================================================
uc_model_code <- '
$PROB Ulcerative Colitis QSP Model v1.0
  Compartments: IFX (2-cmt TMDD), VDZ (1-cmt TMDD),
                TOF (1-cmt oral), OZA+M1 (2-cmt oral),
                UST (depot+central), Cytokines (TNFa/IL17/IL13/IL10),
                Immune cells (Th2/Th17/Treg/Neutrophil),
                Disease state (MayoScore/MucosalHealing/CRP/FC)
  Total ODE compartments: 25

$PARAM
// --------------------------------------------------------
// Infliximab PK (Fasanmade et al. Clin Pharmacokinet 2010)
// --------------------------------------------------------
CL_IFX    = 0.407   // L/day; total clearance
V1_IFX    = 3.28    // L; central volume of distribution
V2_IFX    = 3.57    // L; peripheral volume
Q_IFX     = 0.484   // L/day; inter-compartmental clearance
kon_IFX   = 0.097   // 1/(nM·day); association rate (TMDD)
koff_IFX  = 0.001   // 1/day; dissociation rate
kdeg_IFX  = 0.0693  // 1/day; drug-target complex degradation
Rbase_TNF = 0.5     // nM; baseline free TNF-alpha

// --------------------------------------------------------
// Vedolizumab PK (Rosario et al. Clin Pharmacokinet 2015)
// --------------------------------------------------------
CL_VDZ    = 0.271   // L/day
V1_VDZ    = 5.24    // L
kon_VDZ   = 0.217   // 1/(nM·day)
koff_VDZ  = 0.033   // 1/day
kdeg_VDZ  = 0.0528  // 1/day
Rbase_a4b7 = 2.0    // nM; baseline alpha4beta7+ lymphocytes (proxy)

// --------------------------------------------------------
// Tofacitinib PK (Dowty et al. Drug Metab Dispos 2014)
// --------------------------------------------------------
ka_TOF    = 83.28   // 1/day (3.47/h converted)
CL_TOF    = 1027.2  // L/day (42.8 L/h converted)
V_TOF     = 2088.0  // L (87 L is wrong unit - actually 87 L correct)
F_TOF     = 0.74    // oral bioavailability

// --------------------------------------------------------
// Ozanimod PK + Active Metabolite CC112273
// (D'Ambrosio et al. Br J Pharmacol 2021)
// --------------------------------------------------------
ka_OZA    = 28.8    // 1/day (1.20/h)
CL_OZA    = 1056.0  // L/day (44.0 L/h)
V_OZA     = 3760.0  // L
k_OZA_M1  = 0.432   // 1/day (0.018/h); conversion to CC112273
CL_M1     = 252.0   // L/day (10.5 L/h)
V_M1      = 4850.0  // L

// --------------------------------------------------------
// Ustekinumab PK (Xu et al. J Crohns Colitis 2017; UNIFI)
// --------------------------------------------------------
ka_UST    = 0.44    // 1/day; SC absorption rate
CL_UST    = 0.192   // L/day
V1_UST    = 4.62    // L
F_UST_SC  = 0.615   // SC bioavailability

// --------------------------------------------------------
// Disease Biology Parameters
// --------------------------------------------------------
// Cytokine turnover
kin_TNF   = 0.020   // nM/day; TNF-alpha production baseline
kout_TNF  = 0.040   // 1/day; TNF-alpha clearance
kin_IL17  = 0.015   // nM/day; IL-17A production
kout_IL17 = 0.030   // 1/day
kin_IL13  = 0.010   // nM/day; IL-13 (Th2 dominant in UC)
kout_IL13 = 0.025   // 1/day
kin_IL10  = 0.012   // nM/day; IL-10 (regulatory)
kout_IL10 = 0.030   // 1/day

// Immune cell turnover (cells/uL/day)
kin_Th2   = 10.0    // Th2 production
kout_Th2  = 0.030   // 1/day; Th2 elimination
kin_Th17  = 5.0     // Th17 production
kout_Th17 = 0.030   // 1/day
kin_Treg  = 8.0     // Treg production
kout_Treg = 0.025   // 1/day
kin_Neut  = 50.0    // Neutrophil production (% baseline units)
kout_Neut = 0.15    // 1/day (fast turnover)

// Disease activity dynamics
kdam      = 0.002   // damage rate coefficient
krep      = 0.008   // mucosal repair rate
kCRP      = 0.50    // CRP production scaling
kFC       = 2.00    // Fecal calprotectin scaling
Hill_n    = 2.0     // Hill coefficient
EC50_TNF  = 0.30    // nM; EC50 TNF-driven damage
EC50_IL13 = 0.20    // nM; EC50 IL-13 barrier disruption

// Drug PD parameters (Emax model)
Emax_IFX  = 0.92    // max anti-TNF effect (infliximab)
EC50_IFX  = 1.5     // ug/mL; EC50 for IFX PD effect
Emax_VDZ  = 0.85    // max anti-homing effect (vedolizumab)
EC50_VDZ  = 10.0    // ug/mL
Emax_TOF  = 0.88    // max JAK inhibition (tofacitinib)
EC50_TOF  = 50.0    // ng/mL
Emax_OZA  = 0.75    // max lymphocyte retention (ozanimod)
EC50_OZA  = 0.5     // ng/mL (CC112273)
Emax_UST  = 0.80    // max IL-12/23 inhibition (ustekinumab)
EC50_UST  = 1.0     // ug/mL

$CMT
// Infliximab PK (2-compartment + TMDD)
IFX_C1   // central compartment (mg)
IFX_C2   // peripheral compartment (mg)
IFX_RC   // drug-TNF complex (nM)

// Vedolizumab PK (1-compartment + TMDD)
VDZ_C1   // central compartment (mg)
VDZ_RC   // drug-alpha4beta7 complex (nM)

// Tofacitinib PK (oral 1-compartment)
TOF_GI   // GI absorption compartment (mg)
TOF_C1   // plasma compartment (mg)

// Ozanimod PK + active metabolite
OZA_GI   // GI compartment (mg)
OZA_C1   // ozanimod plasma (mg)
OZA_M1   // CC112273 metabolite (mg)

// Ustekinumab PK (SC depot + central)
UST_DEPOT // SC depot (mg)
UST_C1    // plasma compartment (mg)

// Cytokine dynamics (nM)
TNFa      // TNF-alpha
IL17      // IL-17A
IL13      // IL-13 (Th2 dominant in UC)
IL10      // IL-10 (regulatory)

// Immune cell populations (cells/uL or %)
Th2       // Th2 lymphocytes
Th17      // Th17 lymphocytes
Treg      // Regulatory T cells
Neutrophil // Neutrophil activity index

// Disease state variables
MayoScore     // Mayo total score (0-12, continuous)
MucosalHealing // Mucosal healing index (0-1)
CRP           // C-reactive protein (ratio to baseline)
FC            // Fecal calprotectin (ratio to baseline)

$INIT
// PK compartments start at 0
IFX_C1 = 0,  IFX_C2 = 0,  IFX_RC = 0
VDZ_C1 = 0,  VDZ_RC = 0
TOF_GI = 0,  TOF_C1 = 0
OZA_GI = 0,  OZA_C1 = 0,  OZA_M1 = 0
UST_DEPOT = 0, UST_C1 = 0

// Disease at baseline (moderate-severe UC)
TNFa = 0.5     // elevated TNF (nM)
IL17 = 0.5     // elevated IL-17
IL13 = 0.4     // elevated IL-13
IL10 = 0.4     // reduced IL-10

// Immune cells at disease baseline
Th2  = 400.0   // cells/uL (elevated gut-homing Th2)
Th17 = 200.0   // cells/uL
Treg = 300.0   // cells/uL (reduced suppression)
Neutrophil = 50.0  // elevated

// Disease activity
MayoScore = 9.0      // moderate-severe UC
MucosalHealing = 0.1 // minimal healing at baseline
CRP = 4.0            // elevated 4x baseline
FC  = 10.0           // elevated 10x baseline

$ODE
// ============================================================
// PHARMACOKINETIC DIFFERENTIAL EQUATIONS
// ============================================================

// --- INFLIXIMAB (2-cmt IV + TMDD with TNFa) ---
// Units: IFX_C1/C2 in mg; conc = mg/V(L) = mg/L = ug/mL (x1000 for mg/3.28L)
double IFX_conc_ugmL = (IFX_C1 / V1_IFX) * 1000.0; // ug/mL

// Emax PD: inhibition of TNF-driven inflammation
double E_IFX = Emax_IFX * pow(IFX_conc_ugmL, Hill_n) /
               (pow(EC50_IFX, Hill_n) + pow(IFX_conc_ugmL, Hill_n));

dxdt_IFX_C1 = -(CL_IFX + Q_IFX) / V1_IFX * IFX_C1
              + Q_IFX / V2_IFX * IFX_C2
              - kon_IFX * (IFX_C1 / V1_IFX) * TNFa
              + koff_IFX * IFX_RC;
dxdt_IFX_C2 = Q_IFX / V1_IFX * IFX_C1 - Q_IFX / V2_IFX * IFX_C2;
dxdt_IFX_RC = kon_IFX * (IFX_C1 / V1_IFX) * TNFa
              - koff_IFX * IFX_RC - kdeg_IFX * IFX_RC;

// --- VEDOLIZUMAB (1-cmt IV + TMDD with alpha4beta7) ---
double VDZ_conc_ugmL = (VDZ_C1 / V1_VDZ) * 1000.0;
double E_VDZ = Emax_VDZ * pow(VDZ_conc_ugmL, Hill_n) /
               (pow(EC50_VDZ, Hill_n) + pow(VDZ_conc_ugmL, Hill_n));

dxdt_VDZ_C1 = -CL_VDZ / V1_VDZ * VDZ_C1
              - kon_VDZ * (VDZ_C1 / V1_VDZ) * Rbase_a4b7
              + koff_VDZ * VDZ_RC;
dxdt_VDZ_RC = kon_VDZ * (VDZ_C1 / V1_VDZ) * Rbase_a4b7
              - koff_VDZ * VDZ_RC - kdeg_VDZ * VDZ_RC;

// --- TOFACITINIB (oral 1-cmt) ---
// Dose in mg; V_TOF = 87 L; convert to ng/mL
double TOF_conc_ngmL = (TOF_C1 / 87.0) * 1e6;  // mg/L -> ng/mL

double E_TOF = Emax_TOF * pow(TOF_conc_ngmL, Hill_n) /
               (pow(EC50_TOF, Hill_n) + pow(TOF_conc_ngmL, Hill_n));

dxdt_TOF_GI = -ka_TOF * TOF_GI;
dxdt_TOF_C1 = ka_TOF * F_TOF * TOF_GI - (CL_TOF / 87.0) * TOF_C1;

// --- OZANIMOD + CC112273 METABOLITE (oral) ---
double OZA_M1_conc_ngmL = (OZA_M1 / V_M1) * 1e6;  // mg/L -> ng/mL (V_M1 in L)
double E_OZA = Emax_OZA * pow(OZA_M1_conc_ngmL, Hill_n) /
               (pow(EC50_OZA, Hill_n) + pow(OZA_M1_conc_ngmL, Hill_n));

dxdt_OZA_GI = -ka_OZA * OZA_GI;
dxdt_OZA_C1 = ka_OZA * OZA_GI - (CL_OZA / V_OZA + k_OZA_M1) * OZA_C1;
dxdt_OZA_M1 = k_OZA_M1 * (V_OZA / V_M1) * OZA_C1 - (CL_M1 / V_M1) * OZA_M1;

// --- USTEKINUMAB (SC depot -> plasma) ---
double UST_conc_ugmL = (UST_C1 / V1_UST) * 1000.0;
double E_UST = Emax_UST * pow(UST_conc_ugmL, Hill_n) /
               (pow(EC50_UST, Hill_n) + pow(UST_conc_ugmL, Hill_n));

dxdt_UST_DEPOT = -ka_UST * UST_DEPOT;
dxdt_UST_C1   = ka_UST * F_UST_SC * UST_DEPOT - CL_UST / V1_UST * UST_C1;

// ============================================================
// COMBINED DRUG EFFECTS FOR DISEASE MODEL
// ============================================================
double inh_TNF     = E_IFX;           // IFX neutralizes TNFa
double inh_homing  = E_VDZ;           // VDZ blocks gut lymphocyte homing
double inh_JAK     = E_TOF;           // TOF inhibits JAK1/3 -> STAT signaling
double inh_S1P     = E_OZA;           // OZA retains lymphocytes in lymph nodes
double inh_IL12_23 = E_UST;           // UST neutralizes IL-12/23

// ============================================================
// CYTOKINE DYNAMICS
// ============================================================

// TNF-alpha: elevated in active UC, suppressed by anti-TNF
dxdt_TNFa = kin_TNF * (1.0 + 2.0 * MayoScore / 12.0) * (1.0 - inh_TNF)
            - kout_TNF * TNFa;

// IL-17A: Th17-driven; reduced by JAK inhibition and IL-12/23 blockade
dxdt_IL17 = kin_IL17 * (Th17 / 200.0) * (1.0 - inh_JAK * 0.6) * (1.0 - inh_IL12_23 * 0.4)
            - kout_IL17 * IL17;

// IL-13: Th2-driven (dominant in UC epithelial barrier disruption)
// Reduced by JAK inhibition (STAT6) and gut lymphocyte homing blockade
dxdt_IL13 = kin_IL13 * (Th2 / 400.0) * (1.0 - inh_JAK * 0.7) * (1.0 - inh_homing * 0.5)
            - kout_IL13 * IL13;

// IL-10: regulatory cytokine from Treg cells
// Treg relative expansion under VDZ (VDZ preferentially spares Tregs)
dxdt_IL10 = kin_IL10 * (Treg / 300.0) * (1.0 + inh_homing * 0.2)
            - kout_IL10 * IL10;

// ============================================================
// IMMUNE CELL POPULATION DYNAMICS
// ============================================================

// Th2 cells (gut-homing alpha4beta7+; dominant in UC)
// Reduced by VDZ (gut homing block) and OZA (lymphocyte retention)
dxdt_Th2 = kin_Th2 * (1.0 - inh_homing) * (1.0 - inh_S1P * 0.6)
           - kout_Th2 * Th2;

// Th17 cells
// Reduced by JAK inhibition (STAT3/IL-6 axis) and S1P modulation
dxdt_Th17 = kin_Th17 * (1.0 - inh_JAK * 0.5) * (1.0 - inh_S1P * 0.5) * (1.0 - inh_IL12_23 * 0.3)
            - kout_Th17 * Th17;

// Regulatory T cells (Treg)
// VDZ can relatively spare Tregs; steroid effect not modeled here
dxdt_Treg = kin_Treg * (1.0 + inh_homing * 0.3)
            - kout_Treg * Treg;

// Neutrophil activity index (driven by IL-17 and CXCL8)
dxdt_Neutrophil = kin_Neut * (1.0 + IL17 / 0.5) * (1.0 - inh_JAK * 0.4)
                  - kout_Neut * Neutrophil;

// ============================================================
// DISEASE STATE DYNAMICS
// ============================================================

// Combined damage driver: TNF, IL-13 (barrier disruption), IL-17 (neutrophil influx)
// Modulated by IL-10 (anti-inflammatory)
double IL10_ratio  = IL10 / 0.4;  // ratio to target level
double TNF_ratio   = TNFa / Rbase_TNF;
double IL17_ratio  = IL17 / 0.5;
double IL13_ratio  = IL13 / 0.4;

double damage_driver = TNF_ratio * (1.0 - IL10_ratio * 0.3)
                       + 0.5 * IL13_ratio
                       + 0.3 * IL17_ratio;

// Mayo Score (0-12 continuous; disease worsens with inflammation, heals with IL-10)
dxdt_MayoScore = kdam * damage_driver * (12.0 - MayoScore)
                 - krep * IL10_ratio * MayoScore;

// Mucosal Healing Index (0=unhealed, 1=fully healed)
dxdt_MucosalHealing = krep * IL10_ratio * (1.0 - MucosalHealing)
                      - kdam * damage_driver * MucosalHealing;

// CRP (ratio to baseline; log-linear production driven by TNF/IL-17/IL-6)
dxdt_CRP = kCRP * (TNF_ratio + 0.5 * IL17_ratio - 0.5 * IL10_ratio - 1.0)
           - 0.15 * CRP;

// Fecal calprotectin (ratio to baseline; driven by neutrophil influx and disease activity)
dxdt_FC = kFC * (Neutrophil / 50.0 - 1.0 + 0.5 * MayoScore / 12.0)
          - 0.10 * FC;

$TABLE
// ============================================================
// DERIVED PK CONCENTRATIONS FOR OUTPUT
// ============================================================
double IFX_conc   = (IFX_C1 / V1_IFX) * 1000.0;   // ug/mL
double VDZ_conc   = (VDZ_C1 / V1_VDZ) * 1000.0;   // ug/mL
double TOF_conc   = (TOF_C1 / 87.0) * 1e6;         // ng/mL
double UST_conc   = (UST_C1 / V1_UST) * 1000.0;   // ug/mL
double OZA_M1_conc = (OZA_M1 / V_M1) * 1e6;       // ng/mL

// ============================================================
// DERIVED CLINICAL ENDPOINTS
// ============================================================
double Mayo_total     = MayoScore;
double pMayo          = MayoScore * 0.75;  // partial Mayo (3-component)
double MH_index       = MucosalHealing;
double CRP_val        = CRP + 1.0;         // add baseline offset
double FC_val         = FC + 50.0;         // ug/g (baseline ~50)

// Response thresholds (binary indicators)
double clin_remission = (MayoScore <= 2.0) ? 1.0 : 0.0;    // Mayo <= 2
double clin_response  = (MayoScore <= 6.0) ? 1.0 : 0.0;    // Mayo reduction >= 3 (proxy)
double MH_resp        = (MucosalHealing >= 0.7) ? 1.0 : 0.0;  // mucosal healing threshold
double deep_remission = ((MayoScore <= 2.0) && (MucosalHealing >= 0.7)) ? 1.0 : 0.0;

$CAPTURE
IFX_conc VDZ_conc TOF_conc UST_conc OZA_M1_conc
Mayo_total pMayo MH_index CRP_val FC_val
TNFa IL17 IL13 IL10 Th2 Th17 Treg Neutrophil
clin_remission clin_response MH_resp deep_remission
'

# Compile model
cat("Compiling UC QSP mrgsolve model...\n")
mod <- mcode("UC_QSP", uc_model_code, quiet = TRUE)
cat("Model compiled successfully.\n")
cat(sprintf("  Compartments: %d\n", length(init(mod))))
cat(sprintf("  Parameters:   %d\n", length(param(mod))))

# ============================================================
# DOSING REGIMENS - 6 TREATMENT SCENARIOS
# ============================================================
# Simulation duration: 365 days (~52 weeks)
TEND   <- 365  # days
DELTA  <- 1    # output every 1 day

# Helper function to build event objects
# Weights assumed: 70 kg body weight for weight-based dosing

BW <- 70  # kg

# ---- SCENARIO 1: PLACEBO ----
ev_placebo <- ev(time = 0, amt = 0, cmt = "IFX_C1", rate = 0)

# ---- SCENARIO 2: INFLIXIMAB (ACT1 regimen: 5 mg/kg IV at wk0,2,6 then Q8W) ----
# 5 mg/kg x 70 kg = 350 mg IV
IFX_dose <- 5 * BW  # mg
ev_IFX <- ev(time = c(0, 14, 42, 98, 154, 210, 266, 322),
              amt  = IFX_dose,
              cmt  = "IFX_C1",
              rate = -2)  # infuse over 2h -> rate=-2 in mrgsolve means rate computed from amt/dur

# ---- SCENARIO 3: VEDOLIZUMAB (GEMINI1: 300 mg IV wk0,2,6 then Q8W) ----
ev_VDZ <- ev(time = c(0, 14, 42, 98, 154, 210, 266, 322),
              amt  = 300,
              cmt  = "VDZ_C1",
              rate = -2)

# ---- SCENARIO 4: TOFACITINIB (OCTAVE: 10 mg BID x56d induction, then 5 mg BID) ----
# BID dosing approximated as once-daily with double amount for simplicity
# 10 mg BID = 20 mg/day for 0-56 days; then 5 mg BID = 10 mg/day
ev_TOF_induction <- ev(time = seq(0, 55, by = 1),
                        amt  = 10,  # 10 mg per dose
                        cmt  = "TOF_GI",
                        ii   = 12/24,  # every 12h expressed as fraction of day
                        addl = 1)      # 2 doses per day
# For mrgsolve simple approach: add doses at 0, 0.5, 1, 1.5... days
tof_ind_times <- sort(c(seq(0, 55.5, by = 0.5)))
tof_main_times <- sort(c(seq(56, 364.5, by = 0.5)))
ev_TOF <- ev(time = c(tof_ind_times, tof_main_times),
              amt  = c(rep(10, length(tof_ind_times)),
                       rep(5, length(tof_main_times))),
              cmt  = "TOF_GI")

# ---- SCENARIO 5: USTEKINUMAB (UNIFI: ~520 mg IV induction, then 90 mg SC Q8W) ----
# IV induction: weight-tiered; use 520 mg for 70 kg patient
ev_UST_IV <- ev(time = 0, amt = 520, cmt = "UST_C1")  # IV goes directly to central
ev_UST_SC <- ev(time = c(56, 112, 168, 224, 280, 336),
                 amt  = 90,
                 cmt  = "UST_DEPOT")  # SC goes to depot
ev_UST <- ev_UST_IV + ev_UST_SC

# ---- SCENARIO 6: OZANIMOD (TRUE NORTH: 0.92 mg QD oral) ----
# Dose in mg; daily oral dosing
ev_OZA <- ev(time = seq(0, 364, by = 1),
              amt  = 0.92,
              cmt  = "OZA_GI")

# ============================================================
# RUN SIMULATIONS
# ============================================================
cat("\nRunning treatment scenario simulations...\n")

sim_time <- seq(0, TEND, by = DELTA)

# Common simulation settings
run_sim <- function(model, events, scenario_name) {
  tryCatch({
    out <- model %>%
      ev(events) %>%
      mrgsim(end = TEND, delta = DELTA) %>%
      as.data.frame()
    out$scenario <- scenario_name
    out
  }, error = function(e) {
    cat(sprintf("  Warning in %s: %s\n", scenario_name, e$message))
    NULL
  })
}

# Simulate each scenario
sim_placebo <- run_sim(mod, ev_placebo,  "1_Placebo")
sim_IFX     <- run_sim(mod, ev_IFX,     "2_Infliximab")
sim_VDZ     <- run_sim(mod, ev_VDZ,     "3_Vedolizumab")
sim_TOF     <- run_sim(mod, ev_TOF,     "4_Tofacitinib")
sim_UST     <- run_sim(mod, ev_UST,     "5_Ustekinumab")
sim_OZA     <- run_sim(mod, ev_OZA,     "6_Ozanimod")

# Combine results
sim_all <- bind_rows(
  sim_placebo, sim_IFX, sim_VDZ, sim_TOF, sim_UST, sim_OZA
)
sim_all$scenario <- factor(sim_all$scenario,
  levels = c("1_Placebo","2_Infliximab","3_Vedolizumab",
             "4_Tofacitinib","5_Ustekinumab","6_Ozanimod"))

cat("Simulations complete.\n")

# ============================================================
# CLINICAL RESPONSE CALCULATIONS
# ============================================================
# Extract responses at key timepoints: Week 8 (day 56), Week 16 (112), Week 52 (365)
timepoints <- c(wk0 = 0, wk8 = 56, wk16 = 112, wk52 = 365)

response_table <- sim_all %>%
  filter(time %in% timepoints) %>%
  mutate(week = case_when(
    time == 0   ~ "Baseline",
    time == 56  ~ "Week 8",
    time == 112 ~ "Week 16",
    time == 365 ~ "Week 52"
  )) %>%
  group_by(scenario, week) %>%
  summarise(
    Mayo_mean     = round(mean(Mayo_total, na.rm=TRUE), 2),
    ClinRemission = round(mean(clin_remission, na.rm=TRUE) * 100, 1),
    ClinResponse  = round(mean(clin_response, na.rm=TRUE) * 100, 1),
    MH_rate       = round(mean(MH_resp, na.rm=TRUE) * 100, 1),
    DeepRemission = round(mean(deep_remission, na.rm=TRUE) * 100, 1),
    MH_index_mean = round(mean(MH_index, na.rm=TRUE), 3),
    CRP_mean      = round(mean(CRP_val, na.rm=TRUE), 2),
    FC_mean       = round(mean(FC_val, na.rm=TRUE), 1),
    .groups = "drop"
  )

cat("\n============================================================\n")
cat("CLINICAL RESPONSE TABLE\n")
cat("============================================================\n")
print(as.data.frame(response_table), row.names = FALSE)

# ============================================================
# PLOTS
# ============================================================

# Color palette for scenarios
scen_colors <- c(
  "1_Placebo"      = "#999999",
  "2_Infliximab"   = "#E41A1C",
  "3_Vedolizumab"  = "#377EB8",
  "4_Tofacitinib"  = "#4DAF4A",
  "5_Ustekinumab"  = "#984EA3",
  "6_Ozanimod"     = "#FF7F00"
)

scen_labels <- c(
  "1_Placebo"      = "Placebo",
  "2_Infliximab"   = "Infliximab (IFX)",
  "3_Vedolizumab"  = "Vedolizumab (VDZ)",
  "4_Tofacitinib"  = "Tofacitinib (TOF)",
  "5_Ustekinumab"  = "Ustekinumab (UST)",
  "6_Ozanimod"     = "Ozanimod (OZA)"
)

theme_uc <- theme_bw() +
  theme(
    legend.position  = "bottom",
    legend.title     = element_blank(),
    strip.background = element_rect(fill = "lightblue"),
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", hjust = 0.5)
  )

# ---- Plot 1: Mayo Score over time ----
p1 <- ggplot(sim_all, aes(x = time / 7, y = Mayo_total, color = scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black", linewidth = 0.7) +
  annotate("text", x = 52, y = 2.3, label = "Remission threshold (Mayo <= 2)",
           hjust = 1, size = 3, fontface = "italic") +
  scale_x_continuous(breaks = c(0, 8, 16, 24, 32, 40, 52), name = "Time (weeks)") +
  scale_y_continuous(name = "Mayo Total Score (0-12)", limits = c(0, 12)) +
  scale_color_manual(values = scen_colors, labels = scen_labels) +
  labs(title = "Mayo Score Over Time by Treatment",
       subtitle = "UC QSP Model Simulation - Moderate-Severe UC at Baseline") +
  theme_uc

# ---- Plot 2: Mucosal Healing Index ----
p2 <- ggplot(sim_all, aes(x = time / 7, y = MH_index, color = scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "darkgreen", linewidth = 0.7) +
  annotate("text", x = 52, y = 0.73, label = "Mucosal healing threshold (>=0.7)",
           hjust = 1, size = 3, fontface = "italic", color = "darkgreen") +
  scale_x_continuous(breaks = c(0, 8, 16, 24, 32, 40, 52), name = "Time (weeks)") +
  scale_y_continuous(name = "Mucosal Healing Index (0-1)", limits = c(0, 1)) +
  scale_color_manual(values = scen_colors, labels = scen_labels) +
  labs(title = "Mucosal Healing Index Over Time",
       subtitle = "UC QSP Model Simulation") +
  theme_uc

# ---- Plot 3: Cytokines (TNFa and IL-13) ----
cyto_long <- sim_all %>%
  select(time, scenario, TNFa, IL17, IL13, IL10) %>%
  pivot_longer(cols = c(TNFa, IL17, IL13, IL10),
               names_to = "cytokine", values_to = "conc")

p3 <- ggplot(cyto_long, aes(x = time / 7, y = conc, color = scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ cytokine, scales = "free_y", nrow = 2,
             labeller = as_labeller(c(TNFa = "TNF-alpha (nM)", IL17 = "IL-17A (nM)",
                                      IL13 = "IL-13 (nM)", IL10 = "IL-10 (nM)"))) +
  scale_x_continuous(breaks = c(0, 8, 16, 24, 52), name = "Time (weeks)") +
  scale_color_manual(values = scen_colors, labels = scen_labels) +
  labs(title = "Cytokine Dynamics by Treatment",
       y = "Concentration (nM)") +
  theme_uc

# ---- Plot 4: Biomarkers (CRP and Fecal Calprotectin) ----
biom_long <- sim_all %>%
  select(time, scenario, CRP_val, FC_val) %>%
  pivot_longer(cols = c(CRP_val, FC_val),
               names_to = "biomarker", values_to = "value")

p4 <- ggplot(biom_long, aes(x = time / 7, y = value, color = scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ biomarker, scales = "free_y",
             labeller = as_labeller(c(CRP_val = "CRP (relative to baseline)",
                                      FC_val = "Fecal Calprotectin (ug/g)"))) +
  scale_x_continuous(breaks = c(0, 8, 16, 24, 52), name = "Time (weeks)") +
  scale_color_manual(values = scen_colors, labels = scen_labels) +
  labs(title = "Biomarker Kinetics Under Treatment",
       y = "Biomarker Level") +
  theme_uc

# ---- Plot 5: PK Concentrations for Biologics ----
pk_long <- sim_all %>%
  select(time, scenario, IFX_conc, VDZ_conc, UST_conc) %>%
  pivot_longer(cols = c(IFX_conc, VDZ_conc, UST_conc),
               names_to = "drug", values_to = "conc_ugmL")

p5 <- ggplot(pk_long %>% filter(conc_ugmL > 0),
             aes(x = time / 7, y = conc_ugmL, color = scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ drug, scales = "free_y",
             labeller = as_labeller(c(IFX_conc = "Infliximab (ug/mL)",
                                      VDZ_conc  = "Vedolizumab (ug/mL)",
                                      UST_conc  = "Ustekinumab (ug/mL)"))) +
  scale_x_continuous(breaks = c(0, 8, 16, 24, 52), name = "Time (weeks)") +
  scale_color_manual(values = scen_colors, labels = scen_labels) +
  labs(title = "Biologic Drug Concentrations (PK Profiles)",
       y = "Concentration (ug/mL)") +
  theme_uc

# ---- Plot 6: Immune Cell Dynamics ----
immune_long <- sim_all %>%
  select(time, scenario, Th2, Th17, Treg, Neutrophil) %>%
  pivot_longer(cols = c(Th2, Th17, Treg, Neutrophil),
               names_to = "cell_type", values_to = "count")

p6 <- ggplot(immune_long, aes(x = time / 7, y = count, color = scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ cell_type, scales = "free_y", nrow = 2,
             labeller = as_labeller(c(Th2 = "Th2 Cells (cells/uL)",
                                      Th17 = "Th17 Cells (cells/uL)",
                                      Treg = "Treg Cells (cells/uL)",
                                      Neutrophil = "Neutrophil Index (%)"))) +
  scale_x_continuous(breaks = c(0, 8, 16, 24, 52), name = "Time (weeks)") +
  scale_color_manual(values = scen_colors, labels = scen_labels) +
  labs(title = "Immune Cell Population Dynamics",
       y = "Cell Count / Activity") +
  theme_uc

# ---- Bar plot: Response rates at Week 8 and Week 52 ----
resp_rates <- response_table %>%
  filter(week %in% c("Week 8", "Week 52")) %>%
  select(scenario, week, ClinRemission, MH_rate, DeepRemission) %>%
  pivot_longer(cols = c(ClinRemission, MH_rate, DeepRemission),
               names_to = "endpoint", values_to = "rate_pct")

p7 <- ggplot(resp_rates, aes(x = scenario, y = rate_pct, fill = scenario)) +
  geom_col(position = "dodge", color = "white") +
  facet_grid(endpoint ~ week,
             labeller = as_labeller(c(
               ClinRemission = "Clinical Remission (%)",
               MH_rate       = "Mucosal Healing (%)",
               DeepRemission = "Deep Remission (%)",
               "Week 8"  = "Week 8",
               "Week 52" = "Week 52"
             ))) +
  geom_text(aes(label = paste0(rate_pct, "%")), vjust = -0.4, size = 2.5) +
  scale_fill_manual(values = scen_colors, labels = scen_labels) +
  scale_x_discrete(labels = scen_labels) +
  scale_y_continuous(limits = c(0, 105), name = "Response Rate (%)") +
  labs(title = "Treatment Response Rates at Week 8 and Week 52",
       x = NULL) +
  theme_uc +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Print plots
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)

# ============================================================
# SUMMARY OUTPUT
# ============================================================
cat("\n============================================================\n")
cat("UC QSP MODEL - SUMMARY\n")
cat("============================================================\n")
cat("Model: Ulcerative Colitis QSP (mrgsolve)\n")
cat(sprintf("Compartments: %d ODE compartments\n", length(init(mod))))
cat("Treatments simulated: 6\n")
cat("  1. Placebo\n")
cat("  2. Infliximab 5mg/kg IV wk0,2,6 then Q8W (ACT1 regimen)\n")
cat("  3. Vedolizumab 300mg IV wk0,2,6 then Q8W (GEMINI1)\n")
cat("  4. Tofacitinib 10mg BID x8wk, then 5mg BID (OCTAVE)\n")
cat("  5. Ustekinumab ~520mg IV induction, 90mg SC Q8W (UNIFI)\n")
cat("  6. Ozanimod 0.92mg QD oral (TRUE NORTH)\n")
cat("\nKey clinical endpoints tracked:\n")
cat("  - Mayo total score (0-12)\n")
cat("  - Partial Mayo score\n")
cat("  - Mucosal healing index (0-1)\n")
cat("  - CRP (serum, relative to baseline)\n")
cat("  - Fecal calprotectin (ug/g)\n")
cat("  - Clinical remission (Mayo <= 2)\n")
cat("  - Mucosal healing response (MH index >= 0.7)\n")
cat("  - Deep remission (Mayo <= 2 + MH >= 0.7)\n")
cat("  - Cytokines: TNF-alpha, IL-17A, IL-13, IL-10\n")
cat("  - Immune cells: Th2, Th17, Treg, Neutrophil\n")

# Return key objects for use in Shiny app
invisible(list(
  model = mod,
  sim_all = sim_all,
  response_table = response_table,
  plots = list(p1=p1, p2=p2, p3=p3, p4=p4, p5=p5, p6=p6, p7=p7)
))
