## ============================================================
## AMD (Age-related Macular Degeneration) QSP Model
## mrgsolve-compatible R file
## ============================================================
## Compartments (20 ODEs):
##   [1] Drug (anti-VEGF) — vitreous
##   [2] Drug — retina/RPE (effect site)
##   [3] Drug — systemic plasma
##   [4] VEGF_free (vitreous/retinal)
##   [5] VEGF_bound (Drug:VEGF complex)
##   [6] VEGFR2_active (signaling)
##   [7] ANG2_free
##   [8] ANG2_bound (Faricimab Ang-2 arm)
##   [9] Complement_C3 (local retinal)
##  [10] Complement_C5 (local retinal)
##  [11] MAC (membrane attack complex)
##  [12] RPE_normal (cells, fraction)
##  [13] RPE_damaged
##  [14] Lipofuscin
##  [15] Drusen (area, mm²)
##  [16] CNV_area (mm²)
##  [17] Fluid_total (CST proxy, μm above baseline)
##  [18] GA_area (geographic atrophy, mm²)
##  [19] BCVA (ETDRS letters)
##  [20] Photoreceptor (fraction surviving)
##
## Key references:
##   - Lanzetta et al. Graefes Arch 2013 (anti-VEGF PK vitreous)
##   - Xu et al. Invest Ophthalmol 2013 (ranibizumab ocular PK)
##   - Holz et al. Ophthalmology 2014 (VIEW AFL)
##   - Heier et al. Ophthalmology 2012 (VIEW1)
##   - Schmidt-Erfurth et al. Br J Ophthalmol 2014 (HARBOR)
##   - Dugel et al. Ophthalmology 2020 (HAWK/HARRIER brolucizumab)
##   - Khanani et al. NEJM 2022 (TENAYA/LUCERNE faricimab)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

code <- '
$PROB AMD QSP Model - Anti-VEGF PK/PD + Disease Progression

$PARAM
// === Drug PK parameters ===
// Ranibizumab (default drug)
Dose_vit = 0.5          // mg IVT dose
Fabs_ret  = 0.3         // fraction drug reaching retina from vitreous
kel_vit   = 0.0965      // /day vitreous elimination (t1/2 ~7.2d Fab)
kel_sys   = 0.693       // /day systemic elimination (t1/2~1d Fab)
ktr_vit2ret = 0.02      // /day vitreous->retina transfer
ktr_vit2sys = 0.005     // /day vitreous->systemic drainage (AH)
ktr_ret2sys = 0.01      // /day retina->systemic clearance
V_vit     = 0.0046      // L vitreous volume (~4.6 mL)

// === VEGF kinetics ===
k_VEGF_syn  = 0.5       // nM/day VEGF production by RPE/hypoxia
k_VEGF_deg  = 0.2       // /day VEGF degradation (t1/2 ~3.5d)
kon_VEGF    = 5.0       // nM^-1 day^-1 drug-VEGF binding
koff_VEGF   = 0.002     // /day drug-VEGF dissociation
Kd_VEGF_RNB = 0.04      // nM VEGF:ranibizumab Kd
Kd_VEGF_AFL = 0.0005    // nM VEGF:aflibercept Kd
Kd_VEGF_BEV = 0.2       // nM VEGF:bevacizumab Kd
Kd_VEGF_FAR = 0.0003    // nM VEGF:faricimab (VEGF arm) Kd
Kd_VEGF_BRO = 0.06      // nM VEGF:brolucizumab Kd
VEGF_baseline = 2.5     // nM baseline vitreous VEGF (wet AMD ~10x normal)

// === VEGFR2 activation ===
k_R2_act  = 0.5         // /day VEGF->VEGFR2 activation rate
k_R2_inact = 0.3        // /day VEGFR2 inactivation
EC50_R2   = 0.5         // nM VEGF for half-max VEGFR2 activation
Hill_R2   = 1.2

// === Ang-2 / Tie2 (faricimab Ang-2 arm) ===
ANG2_baseline = 1.0     // nM vitreous Ang-2 (elevated in AMD)
k_ANG2_syn  = 0.15      // nM/day Ang-2 production
k_ANG2_deg  = 0.1       // /day Ang-2 degradation
kon_ANG2    = 2.0       // nM^-1 day^-1 faricimab Ang-2 arm binding
koff_ANG2   = 0.00087   // /day (Kd~0.9pM for faricimab Ang-2)

// === Complement activation ===
C3_baseline = 5.0       // AU C3 local level (drusen-associated)
C5_baseline = 2.0       // AU C5 local level
k_C3_syn    = 0.3       // AU/day C3 synthesis
k_C3_deg    = 0.1       // /day C3 degradation
k_C5_syn    = 0.15      // AU/day C5 synthesis
k_C5_deg    = 0.08      // /day C5 degradation
k_MAC_form  = 0.05      // /day MAC formation from C5
k_MAC_deg   = 0.1       // /day MAC clearance
E_RPE_MAC   = 0.15      // MAC-to-RPE damage rate coefficient

// === RPE cell dynamics ===
RPE0 = 1.0              // initial RPE fraction (normalized)
k_RPE_death = 0.0002    // /day basal RPE loss rate
E_RPE_ox    = 0.0005    // oxidative stress amplification factor
E_RPE_MAC2  = 0.0008    // MAC-driven RPE loss
k_RPE_repair = 0.00005  // /day marginal RPE self-repair
GA_threshold = 0.6      // RPE fraction below which GA progresses

// === Lipofuscin / drusen ===
k_LF_accum  = 0.001     // AU/day lipofuscin accumulation
k_LF_clear  = 0.0001    // /day lipofuscin clearance (minimal)
k_Drusen_grow = 0.003   // mm^2/yr->day drusen growth (~0.001)
k_Drusen_base = 0.0003  // baseline drusen growth /day

// === CNV dynamics (wet AMD) ===
k_CNV_init  = 0.0005    // /day baseline CNV initiation rate
E_VEGF_CNV  = 0.002     // VEGF amplification of CNV growth
EC50_VEGF_CNV = 1.0     // nM VEGF for half-max CNV growth
k_CNV_reg   = 0.01      // /day CNV regression (VEGF-dependent)
CNV_max     = 25.0      // mm^2 maximum CNV area

// === Fluid dynamics (CST above baseline) ===
k_Fluid_in  = 10.0      // μm/day per unit VEGFR2 activity
k_Fluid_out = 0.15      // /day fluid resorption rate
CST_baseline = 280.0    // μm baseline CST (normal retina)

// === GA progression ===
k_GA_grow   = 0.00055   // mm^2/day ~1.9 mm^2/yr (natural history)
E_MAC_GA    = 0.0003    // MAC-driven GA acceleration
E_RPE_GA    = 0.001     // RPE loss -> GA

// === Photoreceptor dynamics ===
PR0 = 1.0               // normalized photoreceptor count
k_PR_death  = 0.0001    // /day basal loss
E_RPE_PR    = 0.002     // RPE-death -> PR loss coupling
E_fluid_PR  = 0.00005   // excess fluid -> PR damage

// === Visual acuity (ETDRS letters) ===
BCVA0 = 55.0            // letters baseline (wet AMD ~20/80)
k_BCVA_CNV  = 0.8       // letters/mm^2 CNV effect
k_BCVA_GA   = 1.2       // letters/mm^2 GA effect (central)
k_BCVA_fluid = 0.02     // letters/μm excess CST effect
k_BCVA_PR   = 15.0      // letters from PR loss
BCVA_max    = 85.0      // maximum attainable BCVA

// === Disease type flags ===
WET_AMD     = 1         // 1=wet AMD (CNV-driven), 0=dry AMD (GA-driven)
DRUG_TYPE   = 1         // 1=RNB, 2=AFL, 3=BEV, 4=FAR, 5=BRO

// === Dosing schedule ===
n_loading   = 3         // number of loading doses (monthly)
inter_maint = 56.0      // days between maintenance injections (q8w=56)

$CMT
// Compartments
DRUG_VIT    // Drug in vitreous [nM]
DRUG_RET    // Drug in retina/RPE [nM]
DRUG_SYS    // Drug in systemic circulation [mg]
VEGF_FREE   // Free VEGF in retina [nM]
VEGF_BOUND  // Drug:VEGF complex [nM]
VEGFR2_ACT  // Active VEGFR-2 [normalized 0-1]
ANG2_FREE   // Free Ang-2 [nM]
ANG2_BOUND  // Drug:Ang-2 complex [nM]
C3_LOCAL    // Local complement C3 [AU]
C5_LOCAL    // Local complement C5 [AU]
MAC_LOCAL   // Local MAC [AU]
RPE_NORM    // Normal RPE fraction [0-1]
RPE_DAM     // Damaged RPE fraction [0-1]
LIPOFUSCIN  // Lipofuscin accumulation [AU]
DRUSEN      // Drusen area [mm^2]
CNV_AREA    // CNV lesion area [mm^2]
FLUID_EX    // Excess fluid over CST baseline [μm]
GA_AREA     // Geographic atrophy area [mm^2]
BCVA_SCORE  // Best-corrected VA [ETDRS letters]
PR_FRAC     // Photoreceptor fraction [0-1]

$INIT
DRUG_VIT   = 0
DRUG_RET   = 0
DRUG_SYS   = 0
VEGF_FREE  = 2.5      // elevated in wet AMD
VEGF_BOUND = 0
VEGFR2_ACT = 0.3      // partially active at baseline (wet AMD)
ANG2_FREE  = 1.0
ANG2_BOUND = 0
C3_LOCAL   = 5.0
C5_LOCAL   = 2.0
MAC_LOCAL  = 0.5
RPE_NORM   = 1.0
RPE_DAM    = 0.0
LIPOFUSCIN = 0.2
DRUSEN     = 0.8      // pre-existing drusen (intermediate AMD)
CNV_AREA   = 2.0      // small existing CNV at start (wet AMD)
FLUID_EX   = 120.0    // excess fluid (CST~400μm)
GA_AREA    = 0.0      // no GA at start (wet AMD model)
BCVA_SCORE = 55.0     // 20/80 equivalent
PR_FRAC    = 0.95

$ODE

// ============================================================
// Drug PK
// ============================================================
double Kd_eff = (DRUG_TYPE==1) ? Kd_VEGF_RNB :
               (DRUG_TYPE==2) ? Kd_VEGF_AFL  :
               (DRUG_TYPE==3) ? Kd_VEGF_BEV  :
               (DRUG_TYPE==4) ? Kd_VEGF_FAR  :
               Kd_VEGF_BRO;

double kon_eff = 0.693 / (Kd_eff * 1.0);   // derived kon for each drug
// kon = koff/Kd; koff ~ 0.001-0.01/day
double koff_eff = (DRUG_TYPE==1) ? 0.0039 :
                 (DRUG_TYPE==2) ? 0.0005 :
                 (DRUG_TYPE==3) ? 0.002 :
                 (DRUG_TYPE==4) ? 0.0003 :
                 0.006;
double kon_used = koff_eff / Kd_eff;

// VEGF binding
double R_bind_VEGF = kon_used * DRUG_RET * VEGF_FREE - koff_eff * VEGF_BOUND;

// Ang-2 binding (only faricimab DRUG_TYPE==4)
double R_bind_ANG2 = (DRUG_TYPE==4) ?
    (kon_ANG2 * DRUG_RET * ANG2_FREE - koff_ANG2 * ANG2_BOUND) : 0.0;

dxdt_DRUG_VIT = -kel_vit * DRUG_VIT
                - ktr_vit2ret * DRUG_VIT
                - ktr_vit2sys * DRUG_VIT;

dxdt_DRUG_RET = ktr_vit2ret * DRUG_VIT
                - ktr_ret2sys * DRUG_RET
                - kon_used * DRUG_RET * VEGF_FREE
                + koff_eff * VEGF_BOUND
                - R_bind_ANG2;

dxdt_DRUG_SYS = ktr_vit2sys * DRUG_VIT * V_vit
                + ktr_ret2sys * DRUG_RET
                - kel_sys * DRUG_SYS;

// ============================================================
// VEGF dynamics
// ============================================================
double VEGF_upregulation = 1.0 + 0.5 * VEGFR2_ACT + 0.3 * (RPE_DAM);
double VEGF_syn_eff = k_VEGF_syn * VEGF_upregulation;

dxdt_VEGF_FREE = VEGF_syn_eff
                 - k_VEGF_deg * VEGF_FREE
                 - R_bind_VEGF;

dxdt_VEGF_BOUND = R_bind_VEGF
                  - (k_VEGF_deg + kel_vit) * VEGF_BOUND;  // complex cleared

// ============================================================
// VEGFR-2 activation
// ============================================================
double Hill_num = pow(VEGF_FREE, Hill_R2);
double Hill_den = pow(EC50_R2, Hill_R2) + Hill_num;
double VEGFR2_target = (Hill_num / Hill_den);

dxdt_VEGFR2_ACT = k_R2_act * (VEGFR2_target - VEGFR2_ACT);

// ============================================================
// Ang-2 dynamics
// ============================================================
double ANG2_feedback = 1.0 + 0.3 * VEGFR2_ACT;   // VEGF stimulates Ang-2
dxdt_ANG2_FREE = k_ANG2_syn * ANG2_feedback
                 - k_ANG2_deg * ANG2_FREE
                 - R_bind_ANG2;

dxdt_ANG2_BOUND = R_bind_ANG2
                  - k_ANG2_deg * ANG2_BOUND;

// ============================================================
// Complement system (local retinal)
// ============================================================
double C3_stim = 1.0 + 0.2 * LIPOFUSCIN + 0.1 * MAC_LOCAL;
dxdt_C3_LOCAL = k_C3_syn * C3_stim - k_C3_deg * C3_LOCAL;

double C5_stim = 1.0 + 0.15 * C3_LOCAL;
dxdt_C5_LOCAL = k_C5_syn * C5_stim - k_C5_deg * C5_LOCAL;

dxdt_MAC_LOCAL = k_MAC_form * C5_LOCAL - k_MAC_deg * MAC_LOCAL;

// ============================================================
// RPE cell dynamics
// ============================================================
double RPE_ox_damage = E_RPE_ox * (MAC_LOCAL * E_RPE_MAC2 + LIPOFUSCIN * 0.5);
double RPE_death_rate = k_RPE_death + RPE_ox_damage + E_RPE_MAC2 * MAC_LOCAL;

// RPE total = RPE_NORM + RPE_DAM; loss from RPE_NORM to RPE_DAM then death
dxdt_RPE_NORM = -RPE_death_rate * RPE_NORM + k_RPE_repair * (1.0 - RPE_NORM - RPE_DAM);
dxdt_RPE_DAM  = RPE_death_rate * RPE_NORM * 0.5 - 0.05 * RPE_DAM;  // 50% to damaged, 50% lost

// ============================================================
// Lipofuscin accumulation
// ============================================================
double LF_accumulation = k_LF_accum * (2.0 - RPE_NORM);  // more accum as RPE fails
dxdt_LIPOFUSCIN = LF_accumulation - k_LF_clear * LIPOFUSCIN;

// ============================================================
// Drusen dynamics
// ============================================================
double Drusen_growth = k_Drusen_base + k_Drusen_grow * (1.0 - RPE_NORM) * LIPOFUSCIN;
dxdt_DRUSEN = Drusen_growth * (1.0 - DRUSEN / 20.0);  // logistic ceiling 20mm^2

// ============================================================
// CNV dynamics (VEGF-driven)
// ============================================================
double VEGF_eff_CNV = VEGF_FREE / (EC50_VEGF_CNV + VEGF_FREE);
double CNV_growth   = (k_CNV_init + E_VEGF_CNV * VEGF_eff_CNV)
                      * (1.0 - CNV_AREA / CNV_max) * WET_AMD;
// CNV regression with anti-VEGF: proportional to VEGFR2 suppression
double VEGFR2_suppress = 1.0 - VEGFR2_ACT;
double CNV_regress  = k_CNV_reg * VEGFR2_suppress * CNV_AREA;

dxdt_CNV_AREA = (CNV_growth - CNV_regress) * CNV_AREA;

// ============================================================
// Fluid dynamics
// ============================================================
double Fluid_inflow = k_Fluid_in * VEGFR2_ACT * (1.0 + 0.5 * ANG2_FREE / (1.0 + ANG2_FREE));
double Fluid_outflow = k_Fluid_out * FLUID_EX;

dxdt_FLUID_EX = Fluid_inflow - Fluid_outflow;
if(FLUID_EX < 0) FLUID_EX = 0;

// ============================================================
// Geographic Atrophy (dry AMD / late dry component)
// ============================================================
double RPE_loss_rate = fmax(0.0, (1.0 - RPE_NORM));
double GA_growth = (k_GA_grow + E_MAC_GA * MAC_LOCAL + E_RPE_GA * RPE_loss_rate);
dxdt_GA_AREA = GA_growth * (1.0 + GA_AREA / 5.0);  // expanding rim growth

// ============================================================
// Photoreceptors
// ============================================================
double PR_loss_RPE = E_RPE_PR * fmax(0.0, (RPE_death_rate - k_RPE_death));
double PR_loss_fluid = E_fluid_PR * fmax(0.0, FLUID_EX - 50.0);
dxdt_PR_FRAC = -k_PR_death * PR_FRAC - PR_loss_RPE - PR_loss_fluid;
if(PR_FRAC < 0) PR_FRAC = 0;

// ============================================================
// Visual Acuity (ETDRS letters) — clinical endpoint
// ============================================================
double CNV_BCVA_loss  = k_BCVA_CNV  * CNV_AREA;
double GA_BCVA_loss   = k_BCVA_GA   * GA_AREA * 0.5;   // central GA only (x0.5)
double Fluid_BCVA_loss = k_BCVA_fluid * fmax(0.0, FLUID_EX - 30.0);
double PR_BCVA_loss   = k_BCVA_PR   * (1.0 - PR_FRAC);

double BCVA_target = BCVA_max
                     - CNV_BCVA_loss
                     - GA_BCVA_loss
                     - Fluid_BCVA_loss
                     - PR_BCVA_loss;
BCVA_target = fmax(0.0, fmin(BCVA_max, BCVA_target));

dxdt_BCVA_SCORE = 0.05 * (BCVA_target - BCVA_SCORE);  // sluggish VA dynamics

$TABLE
capture VEGF_FREE_out  = VEGF_FREE;
capture VEGF_BOUND_out = VEGF_BOUND;
capture DRUG_VIT_out   = DRUG_VIT;
capture DRUG_RET_out   = DRUG_RET;
capture VEGFR2_out     = VEGFR2_ACT;
capture CNV_out        = CNV_AREA;
capture CST_out        = CST_baseline + FLUID_EX;
capture GA_out         = GA_AREA;
capture BCVA_out       = BCVA_SCORE;
capture RPE_out        = RPE_NORM;
capture PR_out         = PR_FRAC;
capture MAC_out        = MAC_LOCAL;
capture Drusen_out     = DRUSEN;
capture ANG2_out       = ANG2_FREE;
capture BCVA_change    = BCVA_SCORE - 55.0;

$CAPTURE VEGF_FREE_out VEGF_BOUND_out DRUG_VIT_out DRUG_RET_out
        VEGFR2_out CNV_out CST_out GA_out BCVA_out
        RPE_out PR_out MAC_out Drusen_out ANG2_out BCVA_change
'

mod <- mcode("AMD_QSP", code)

## ============================================================
## DOSING FUNCTIONS
## ============================================================

make_events <- function(drug_type = 1,
                        n_load = 3,
                        load_interval = 28,
                        maint_interval = 56,
                        n_maint = 10,
                        cmt = 1) {
  # Loading doses
  load_times <- seq(0, (n_load - 1) * load_interval, by = load_interval)
  # Maintenance doses
  maint_start <- max(load_times) + maint_interval
  maint_times <- seq(maint_start, maint_start + (n_maint - 1) * maint_interval,
                     by = maint_interval)
  all_times <- c(load_times, maint_times)

  # Dose amounts (converted to approximate nM in vitreous ~4.6 mL)
  # MW: RNB=48kDa, AFL=115kDa, BEV=149kDa, FAR=150kDa, BRO=26kDa
  dose_mg <- c(0.5, 2.0, 1.25, 6.0, 6.0)[drug_type]
  MW_kDa  <- c(48,  115,  149,  150,  26)[drug_type]
  V_vit_L <- 0.0046
  dose_nM <- (dose_mg / MW_kDa) / V_vit_L * 1000

  ev <- ev(time = all_times, amt = dose_nM, cmt = cmt, addl = 0)
  return(ev)
}

## ============================================================
## SCENARIO 1: Ranibizumab q4w × 3 loading → q8w maintenance
## ============================================================
cat("\n--- SCENARIO 1: Ranibizumab (0.5mg) q4w×3 → q8w ---\n")

ev_rnb <- make_events(drug_type = 1, n_load = 3, load_interval = 28,
                      maint_interval = 56, n_maint = 10)

out_rnb <- mod %>%
  param(DRUG_TYPE = 1, WET_AMD = 1, BCVA0 = 55) %>%
  ev(ev_rnb) %>%
  mrgsim(end = 730, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Ranibizumab q4w→q8w")

## ============================================================
## SCENARIO 2: Aflibercept q4w × 3 loading → q8w maintenance
## ============================================================
cat("--- SCENARIO 2: Aflibercept (2mg) q4w×3 → q8w ---\n")

ev_afl <- make_events(drug_type = 2, n_load = 3, load_interval = 28,
                      maint_interval = 56, n_maint = 10)

out_afl <- mod %>%
  param(DRUG_TYPE = 2, WET_AMD = 1,
        kel_vit = 0.099,  # AFL t1/2 ~7d in vitreous
        BCVA0 = 55) %>%
  ev(ev_afl) %>%
  mrgsim(end = 730, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Aflibercept q4w→q8w")

## ============================================================
## SCENARIO 3: Faricimab (dual anti-VEGF/Ang-2) q4w × 4 → q16w
## ============================================================
cat("--- SCENARIO 3: Faricimab (6mg) q4w×4 → q16w ---\n")

ev_far <- make_events(drug_type = 4, n_load = 4, load_interval = 28,
                      maint_interval = 112, n_maint = 7)

out_far <- mod %>%
  param(DRUG_TYPE = 4, WET_AMD = 1,
        kel_vit = 0.099,
        BCVA0 = 55) %>%
  ev(ev_far) %>%
  mrgsim(end = 730, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Faricimab q4w→q16w (T&E)")

## ============================================================
## SCENARIO 4: Brolucizumab q6w × 3 → q8w or q12w
## ============================================================
cat("--- SCENARIO 4: Brolucizumab (6mg) q6w×3 → q12w ---\n")

ev_bro <- make_events(drug_type = 5, n_load = 3, load_interval = 42,
                      maint_interval = 84, n_maint = 8)

out_bro <- mod %>%
  param(DRUG_TYPE = 5, WET_AMD = 1,
        kel_vit = 0.173,  # BRO t1/2 ~4d (small scFv, faster)
        BCVA0 = 55) %>%
  ev(ev_bro) %>%
  mrgsim(end = 730, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Brolucizumab q6w→q12w")

## ============================================================
## SCENARIO 5: Untreated natural history (wet AMD)
## ============================================================
cat("--- SCENARIO 5: Natural history (no treatment) ---\n")

out_noTx <- mod %>%
  param(DRUG_TYPE = 1, WET_AMD = 1, BCVA0 = 55) %>%
  mrgsim(end = 730, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Natural History (No Treatment)")

## ============================================================
## SCENARIO 6: Dry AMD natural history + AREDS2 supplements
## ============================================================
cat("--- SCENARIO 6: Dry AMD with AREDS2 supplements ---\n")

out_dry_areds <- mod %>%
  param(DRUG_TYPE = 1, WET_AMD = 0,
        BCVA0 = 70,        # better baseline in dry AMD
        CNV_AREA = 0,      # no CNV
        FLUID_EX = 0,      # no fluid
        GA_AREA = 1.5,     # existing GA
        k_RPE_death  = 0.00015,  # AREDS reduces RPE loss by ~25%
        k_GA_grow    = 0.00045,  # AREDS reduces GA growth ~18%
        E_RPE_ox     = 0.0003)   %>%
  mrgsim(end = 1460, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Dry AMD + AREDS2")

## ============================================================
## COMBINE RESULTS AND PLOT
## ============================================================
results_wet <- bind_rows(out_rnb, out_afl, out_far, out_bro, out_noTx)

## Plot 1: BCVA over time
p1 <- ggplot(results_wet, aes(x = time, y = BCVA_out, color = scenario)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = c(55, 70), linetype = "dashed", color = "gray50") +
  annotate("text", x = 700, y = 71, label = "+15 letters\n(clinically meaningful)", size = 3) +
  labs(title = "AMD QSP Model — BCVA Outcomes (Wet AMD)",
       subtitle = "Comparison of Anti-VEGF Regimens vs Natural History",
       x = "Time (days)", y = "BCVA (ETDRS letters)",
       color = "Treatment") +
  theme_bw(base_size = 12) +
  scale_x_continuous(breaks = seq(0, 730, 90)) +
  coord_cartesian(ylim = c(20, 80))

print(p1)

## Plot 2: CNV Area over time
p2 <- ggplot(results_wet, aes(x = time, y = CNV_out, color = scenario)) +
  geom_line(size = 1.2) +
  labs(title = "CNV Lesion Area over Time",
       x = "Time (days)", y = "CNV Area (mm²)", color = "Treatment") +
  theme_bw(base_size = 12)
print(p2)

## Plot 3: CST (fluid) over time
p3 <- ggplot(results_wet, aes(x = time, y = CST_out, color = scenario)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 280, linetype = "dashed", color = "black") +
  annotate("text", x = 650, y = 270, label = "Normal CST ~280μm", size = 3) +
  labs(title = "Central Subfield Thickness (CST) — OCT Proxy",
       x = "Time (days)", y = "CST (μm)", color = "Treatment") +
  theme_bw(base_size = 12)
print(p3)

## Plot 4: VEGF free vs drug
p4 <- ggplot(results_wet %>% filter(scenario != "Natural History (No Treatment)"),
             aes(x = time, y = VEGF_FREE_out, color = scenario)) +
  geom_line(size = 1.0) +
  labs(title = "Free VEGF-A in Vitreous/Retina",
       x = "Time (days)", y = "Free VEGF (nM)", color = "Treatment") +
  theme_bw(base_size = 12)
print(p4)

## Plot 5: PK in vitreous
p5 <- ggplot(results_wet %>% filter(scenario != "Natural History (No Treatment)"),
             aes(x = time, y = DRUG_VIT_out, color = scenario)) +
  geom_line(size = 1.0) +
  scale_y_log10() +
  labs(title = "Drug Concentration in Vitreous (log scale)",
       x = "Time (days)", y = "Drug [nM] (log)", color = "Treatment") +
  theme_bw(base_size = 12)
print(p5)

## Plot 6: Dry AMD progression
p6 <- ggplot(out_dry_areds, aes(x = time, y = GA_out)) +
  geom_line(size = 1.2, color = "#E74C3C") +
  geom_area(fill = "#E74C3C", alpha = 0.2) +
  labs(title = "Dry AMD — Geographic Atrophy Progression (with AREDS2)",
       x = "Time (days)", y = "GA Area (mm²)") +
  theme_bw(base_size = 12)
print(p6)

## ============================================================
## SUMMARY TABLE: 1-year outcomes
## ============================================================
summary_1yr <- results_wet %>%
  filter(time == 365) %>%
  select(scenario, BCVA_out, BCVA_change, CNV_out, CST_out,
         VEGF_FREE_out, RPE_out, GA_out) %>%
  rename(
    `BCVA (letters)` = BCVA_out,
    `ΔBCVA from BL`  = BCVA_change,
    `CNV Area (mm²)` = CNV_out,
    `CST (μm)`       = CST_out,
    `Free VEGF (nM)` = VEGF_FREE_out,
    `RPE fraction`   = RPE_out,
    `GA Area (mm²)`  = GA_out
  )

cat("\n=== 1-YEAR OUTCOMES SUMMARY ===\n")
print(summary_1yr, digits = 3)

## ============================================================
## SENSITIVITY ANALYSIS: Kd effect on BCVA gain
## ============================================================
cat("\n--- Sensitivity: Kd impact on 1yr BCVA ---\n")

Kd_values <- c(0.001, 0.01, 0.04, 0.1, 0.5, 1.0)  # nM
sens_res <- lapply(Kd_values, function(kd) {
  out <- mod %>%
    param(DRUG_TYPE = 1, WET_AMD = 1,
          Kd_VEGF_RNB = kd) %>%
    ev(ev_rnb) %>%
    mrgsim(end = 365, delta = 7) %>%
    as.data.frame()
  data.frame(Kd_nM = kd,
             BCVA_yr1 = tail(out$BCVA_out, 1),
             BCVA_change = tail(out$BCVA_out, 1) - 55)
})
sens_df <- bind_rows(sens_res)
cat("Kd sensitivity (Ranibizumab, 1-year):\n")
print(sens_df)

cat("\n=== AMD QSP Model — COMPLETE ===\n")
cat("Deliverables: 20-ODE model, 6 clinical scenarios, PK/PD visualization\n")
cat("Key endpoints: BCVA, CNV area, CST, GA area, free VEGF, RPE fraction\n")
