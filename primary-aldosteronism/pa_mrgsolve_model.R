## =============================================================================
## Primary Aldosteronism (Conn's Syndrome) — mrgsolve QSP Model
## =============================================================================
## Model: RAAS cascade + adrenal steroidogenesis + renal tubule signaling +
##        cardiovascular hemodynamics + target organ damage + drug PK/PD
##
## Compartments (23 ODEs):
##   PK:  [1]  C_spiro     — Spironolactone plasma (ng/mL)
##        [2]  C_canrenone — Canrenone (active metabolite) (ng/mL)
##        [3]  C_eple      — Eplerenone plasma (ng/mL)
##        [4]  C_fine      — Finerenone plasma (ng/mL)
##        [5]  C_acei      — ACEi (ramiprilat) plasma (ng/mL)
##        [6]  C_ccb       — CCB (amlodipine) plasma (ng/mL)
##   RAAS/Steroids:
##        [7]  Renin       — Plasma renin activity (ng/mL/h)
##        [8]  AngII       — Angiotensin II (pmol/L)
##        [9]  Aldo        — Aldosterone plasma (ng/dL)
##   Renal/Ions:
##        [10] ENaC_act    — ENaC activity (relative units, 1=basal)
##        [11] Na_plasma   — Plasma sodium (mEq/L)
##        [12] K_plasma    — Plasma potassium (mEq/L)
##        [13] HCO3        — Plasma bicarbonate (mEq/L)
##        [14] Vol_plasma  — Plasma volume (relative units, 1=normal)
##   Cardiovascular:
##        [15] MAP         — Mean arterial pressure (mmHg)
##        [16] TPR         — Total peripheral resistance (relative units)
##   Renal function:
##        [17] GFR_c       — GFR (mL/min/1.73m²)
##   Target organ damage:
##        [18] CardFibrosis— Cardiac fibrosis index (0–1)
##        [19] LVMi        — LV mass index (g/m²)
##   Adrenal:
##        [20] APA_activity— APA autonomous aldosterone secretion (relative)
##        [21] CYP11B2_act — Aldosterone synthase activity
##   Biomarkers:
##        [22] ARR_model   — Modeled ARR (ng/dL per ng/mL/h)
##        [23] HOMA_proxy  — Metabolic insulin resistance proxy
##
## References:
##   Funder JW et al. J Clin Endocrinol Metab 2016 (Endocrine Society Guidelines)
##   Rossi GP et al. J Hypertens 2016 (PAPY Study); Williams B et al. Lancet 2018
##   Bakris GL et al. NEJM 2020 (FIDELIO-DKD — finerenone); Conn JW 1955
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ---------------------------------------------------------------------------
## Model code
## ---------------------------------------------------------------------------
pa_code <- '
$PROB Primary Aldosteronism (PA) QSP mrgsolve Model v1.0

$PARAM
// ---- Drug doses ----
dose_spiro = 0       // Spironolactone daily dose (mg)
dose_eple  = 0       // Eplerenone daily dose (mg)
dose_fine  = 0       // Finerenone daily dose (mg)
dose_acei  = 0       // ACEi (ramipril equivalent) daily dose (mg)
dose_ccb   = 0       // CCB (amlodipine) daily dose (mg)
surgery    = 0       // Adrenalectomy (0=no, 1=yes, effect applied at t=0)

// ---- PK parameters ----
ka_spiro   = 1.2    // Spironolactone absorption rate (1/h)
CL_spiro   = 22.0   // Spironolactone clearance (L/h)
Vd_spiro   = 80.0   // Spironolactone volume (L)
km_spiro   = 0.4    // Spironolactone→canrenone conversion (1/h)
CL_canr    = 0.25   // Canrenone clearance (1/h) — T½~20h
Vd_canr    = 160.0

ka_eple    = 1.0
CL_eple    = 18.0
Vd_eple    = 50.0

ka_fine    = 0.8
CL_fine    = 25.0
Vd_fine    = 70.0

ka_acei    = 1.5
CL_acei    = 15.0
Vd_acei    = 100.0

ka_ccb     = 0.5
CL_ccb     = 7.0
Vd_ccb     = 2100.0

// ---- Disease parameters ----
// Baseline PCOS values (PA phenotype)
APA_severity = 1.0  // 0=BAH, 1.0=typical APA, 2.0=severe APA
BAH_severity = 0.5  // BAH additive aldosterone excess (0–1)

// RAAS parameters
kout_Renin   = 2.0    // Renin clearance (1/h)
kin_Renin    = 0.6    // Renin basal production
IC50_Aldo_Renin = 15  // Aldosterone IC50 for renin suppression (ng/dL)
IC50_MAP_Renin  = 100 // MAP IC50 for renin (pressure suppression)
EC50_AngII_Renin = 20 // Ang II IC50 for renin negative feedback

kout_AngII   = 6.0    // Ang II clearance (1/h) — T½~2-4 min
kin_AngII    = 5.0    // Ang II basal production from renin
kRenin_AngII = 3.0    // Renin-driven Ang II production coefficient
IC50_ACEI_ACE = 5.0   // ACEi IC50 for ACE inhibition (ng/mL)

kout_Aldo    = 0.7    // Aldosterone clearance (1/h) — T½~20 min
kin_Aldo_bas = 2.0    // Basal ZG aldosterone production
kAngII_Aldo  = 1.5    // AngII-stimulated aldosterone production
kAPA_Aldo    = 8.0    // APA autonomous aldosterone production
kBAH_Aldo    = 3.0    // BAH excess aldosterone

// Renal tubule parameters
kout_ENaC    = 0.05   // ENaC activity clearance (1/h) — slow turnover
kin_ENaC     = 0.05   // ENaC basal activity
kAldo_ENaC   = 0.08   // Aldo-driven ENaC upregulation
EC50_Aldo_ENaC = 12   // Aldo EC50 for ENaC (ng/dL)

kout_Na      = 0.1    // Na+ equilibration rate
Na_set       = 140    // Normal serum Na+ (mEq/L)
kENaC_Na     = 0.5    // ENaC-driven Na+ retention

kout_K       = 0.08   // K+ equilibration rate
K_set        = 4.2    // Normal serum K+ (mEq/L)
kENaC_K_loss = 0.6    // ENaC-driven K+ excretion (ROMK)

kout_HCO3    = 0.04   // HCO3- equilibration
HCO3_set     = 24     // Normal HCO3- (mEq/L)
kK_HCO3      = 0.2    // Hypokalemia-driven alkalosis

kout_Vol     = 0.02   // Plasma volume equilibration
Vol_set      = 1.0    // Normal plasma volume
kNa_Vol      = 0.3    // Na+ retention → volume expansion

// Cardiovascular parameters
kout_MAP     = 0.1    // MAP equilibration (1/h)
MAP_basal    = 90     // Normal MAP (mmHg)
kVol_MAP     = 15     // Volume expansion → MAP increase
kTPR_MAP     = 20     // TPR → MAP contribution
kAngII_TPR   = 0.3    // AngII → vasoconstriction

kout_TPR     = 0.05   // TPR equilibration
TPR_basal    = 1.0    // Normal TPR (relative)
kAngII_TPR_stim = 0.4 // AngII stimulates TPR

kout_GFR     = 0.005  // GFR slow turnover
GFR_set      = 90     // Normal GFR (mL/min/1.73m²)
kMap_GFR     = 0.5    // MAP → pressure-induced hyperfiltration
kKdrop_GFR   = 0.02   // Hypokalemia → renal tubular damage

// Target organ damage parameters
kout_CardFib = 0.001  // Cardiac fibrosis resolution (very slow)
kAldo_CardFib = 0.0005 // Aldo-driven fibrosis rate
kMAP_CardFib  = 0.0003 // Pressure-driven fibrosis

kout_LVMi    = 0.001  // LVM index regression
LVMi_normal  = 80     // Normal LVMi (g/m²)
kMAP_LVMi    = 0.5    // MAP-driven LVH
kAldo_LVMi   = 0.2    // Aldo-driven LVH (excess over pressure)

// Adrenal
kout_APA     = 0.0001 // APA growth (very slow)
kout_CYP11B2 = 0.1    // CYP11B2 regulation
surgery_effect = 1.0  // 1.0 = complete cure (0 = no effect)

// Drug PD parameters
IC50_spiro_MR  = 30   // Spironolactone+canrenone IC50 for MR blockade (ng/mL)
IC50_eple_MR   = 25   // Eplerenone IC50 for MR blockade
IC50_fine_MR   = 8    // Finerenone IC50 for MR blockade (high affinity)
Emax_MR_ENaC   = 0.9  // Max fraction ENaC reduction by MR blockade
IC50_ACEI_AngII = 8   // ACEi IC50 for AngII reduction
EC50_ccb_TPR   = 20   // CCB EC50 for TPR reduction (ng/mL)
Emax_ccb_TPR   = 0.35 // Max TPR reduction by CCB

$CMT
C_spiro C_canrenone C_eple C_fine C_acei C_ccb
Renin_c AngII_c Aldo_c
ENaC_act Na_c K_c HCO3_c Vol_c
MAP_c TPR_c GFR_c
CardFib LVMi_c
APA_act CYP11B2_c
ARR_c HOMA_proxy

$MAIN
// Compute drug-receptor occupancy at each time step
double MR_block_spiro = (C_spiro + C_canrenone) / (IC50_spiro_MR + C_spiro + C_canrenone);
double MR_block_eple  = C_eple / (IC50_eple_MR + C_eple);
double MR_block_fine  = C_fine / (IC50_fine_MR + C_fine);
double MR_block_total = 1.0 - (1.0 - MR_block_spiro) * (1.0 - MR_block_eple) * (1.0 - MR_block_fine);
if (MR_block_total > 0.98) MR_block_total = 0.98;

double ACE_block = C_acei / (IC50_ACEI_ACE + C_acei);
double CCB_effect = Emax_ccb_TPR * C_ccb / (EC50_ccb_TPR + C_ccb);

// Surgery effect on APA
double APA_remaining = (surgery > 0.5) ? (1.0 - surgery_effect) : 1.0;

// Safe compartment values
double Aldo_safe = (Aldo_c < 0.1) ? 0.1 : Aldo_c;
double K_safe    = (K_c < 1.0) ? 1.0 : K_c;
double AngII_safe = (AngII_c < 0.01) ? 0.01 : AngII_c;
double Renin_safe = (Renin_c < 0.01) ? 0.01 : Renin_c;

$ODE
// =====================================================================
// PK: Drug concentrations
// =====================================================================
dxdt_C_spiro     = ka_spiro   * dose_spiro/24.0  - (CL_spiro/Vd_spiro)*C_spiro - km_spiro*C_spiro;
dxdt_C_canrenone = km_spiro*C_spiro*Vd_spiro/Vd_canr - (CL_canr/Vd_canr)*C_canrenone;
dxdt_C_eple      = ka_eple    * dose_eple/24.0   - (CL_eple/Vd_eple)*C_eple;
dxdt_C_fine      = ka_fine    * dose_fine/24.0   - (CL_fine/Vd_fine)*C_fine;
dxdt_C_acei      = ka_acei    * dose_acei/24.0   - (CL_acei/Vd_acei)*C_acei;
dxdt_C_ccb       = ka_ccb     * dose_ccb/24.0    - (CL_ccb/Vd_ccb)*C_ccb;

// =====================================================================
// Renin (ng/mL/h)
// Suppressed by aldosterone (volume feedback) and AngII (short feedback)
// Recovered after adrenalectomy
// =====================================================================
double Aldo_supp_Renin = IC50_Aldo_Renin / (IC50_Aldo_Renin + Aldo_safe);
double MAP_supp_Renin  = IC50_MAP_Renin / (IC50_MAP_Renin + MAP_c);
double AngII_supp_Renin = EC50_AngII_Renin / (EC50_AngII_Renin + AngII_safe);
double Renin_surgery_release = (surgery > 0.5) ? 3.0 : 1.0; // renin rebounds after APA removal
dxdt_Renin_c = kin_Renin * Aldo_supp_Renin * MAP_supp_Renin * AngII_supp_Renin * Renin_surgery_release
              - kout_Renin * Renin_c;

// =====================================================================
// Angiotensin II (pmol/L)
// Produced from Ang I by ACE (inhibited by ACEi)
// ACE2 counterbalances
// =====================================================================
double ACE_activity = 1.0 - ACE_block;
dxdt_AngII_c = kRenin_AngII * Renin_safe * ACE_activity * kin_AngII
             - kout_AngII * AngII_c;

// =====================================================================
// Aldosterone (ng/dL)
// Sources: physiological (AngII-driven from ZG) + APA autonomous + BAH excess
// MR antagonists do NOT reduce aldosterone; ACEi partially reduces via AngII
// Surgery eliminates APA source
// =====================================================================
double AngII_stim_Aldo = kAngII_Aldo * AngII_safe / (5.0 + AngII_safe);
double APA_source = kAPA_Aldo * APA_remaining * APA_severity;
double BAH_source = kBAH_Aldo * BAH_severity;
double ACEi_Aldo_supp = 1.0 - ACE_block * 0.4; // partial suppression via Ang II
dxdt_Aldo_c = (kin_Aldo_bas + AngII_stim_Aldo + APA_source + BAH_source) * ACEi_Aldo_supp
            - kout_Aldo * Aldo_c;

// =====================================================================
// ENaC activity (relative, 1=basal)
// Driven by aldosterone via MR-SGK1-Nedd4-2 pathway
// Blocked by MR antagonists
// =====================================================================
double Aldo_ENaC_stim = Aldo_safe / (EC50_Aldo_ENaC + Aldo_safe);
double MR_block_ENaC  = Emax_MR_ENaC * MR_block_total;
dxdt_ENaC_act = kin_ENaC * (1.0 + kAldo_ENaC * Aldo_ENaC_stim) * (1.0 - MR_block_ENaC)
             - kout_ENaC * ENaC_act;

// =====================================================================
// Plasma Sodium (mEq/L)
// ENaC excess → Na+ retention → hypernatremia (often high-normal)
// =====================================================================
double ENaC_safe = (ENaC_act < 0.01) ? 0.01 : ENaC_act;
dxdt_Na_c = kENaC_Na * (ENaC_safe - 1.0) * 0.5 - kout_Na * (Na_c - Na_set);

// =====================================================================
// Plasma Potassium (mEq/L) — KEY biomarker of PA
// ENaC↑ → lumenal electronegativity → ROMK drives K+ secretion
// =====================================================================
dxdt_K_c = -(kENaC_K_loss * ENaC_safe * 0.4) - kout_K * (K_c - K_set);

// =====================================================================
// Plasma HCO3- (mEq/L) — metabolic alkalosis from H+ excretion
// Also driven by K+ depletion (H+ shifts intracellularly as K+ exits cells)
// =====================================================================
double K_deficit = (K_safe < K_set) ? (K_set - K_safe) : 0.0;
dxdt_HCO3_c = kK_HCO3 * K_deficit - kout_HCO3 * (HCO3_c - HCO3_set);

// =====================================================================
// Plasma Volume (relative, 1=normal)
// Na+ retention → volume expansion → hypertension
// =====================================================================
dxdt_Vol_c = kNa_Vol * (Na_c - Na_set) * 0.01 - kout_Vol * (Vol_c - Vol_set);

// =====================================================================
// Mean Arterial Pressure (mmHg)
// Driven by volume expansion and increased TPR
// Reduced by antihypertensives (MR antagonists, CCB, ACEi)
// =====================================================================
double Vol_effect = kVol_MAP * (Vol_c - Vol_set);
double TPR_effect = kTPR_MAP * (TPR_c - TPR_basal);
double MAP_target = MAP_basal + Vol_effect + TPR_effect;
dxdt_MAP_c = kout_MAP * (MAP_target - MAP_c);

// =====================================================================
// Total Peripheral Resistance (relative)
// AngII → AT1R → VSMC contraction; CCB relaxes VSMC; MR antagonists reduce over time
// =====================================================================
double AngII_TPR_stim = kAngII_TPR_stim * AngII_safe / (10.0 + AngII_safe);
double TPR_target = TPR_basal + AngII_TPR_stim - CCB_effect - MR_block_total * 0.15;
dxdt_TPR_c = kout_TPR * (TPR_target - TPR_c);

// =====================================================================
// GFR (mL/min/1.73m²)
// Hyperfiltration from elevated MAP; hypokalemia → tubular injury → GFR loss
// =====================================================================
double MAP_GFR_effect = (MAP_c > MAP_basal) ? kMap_GFR * (MAP_c - MAP_basal) * 0.02 : 0.0;
double K_GFR_injury   = kKdrop_GFR * K_deficit;
dxdt_GFR_c = kout_GFR * (GFR_set + MAP_GFR_effect - K_GFR_injury - GFR_c);

// =====================================================================
// Cardiac Fibrosis Index (0–1)
// Excess aldosterone drives cardiac MR → TGF-β → collagen beyond pressure effects
// Finerenone preferentially reduces cardiac fibrosis
// =====================================================================
double fine_cardio_protect = 1.0 - 0.7 * C_fine / (IC50_fine_MR + C_fine); // extra cardiac effect
double Aldo_fibrosis_drive = kAldo_CardFib * (Aldo_safe - 8.0) * (Aldo_safe > 8.0 ? 1.0 : 0.0);
double MAP_fibrosis_drive  = kMAP_CardFib * (MAP_c - MAP_basal) * (MAP_c > MAP_basal ? 1.0 : 0.0);
dxdt_CardFib = (Aldo_fibrosis_drive + MAP_fibrosis_drive) * fine_cardio_protect
            - kout_CardFib * CardFib;

// =====================================================================
// LV Mass Index (g/m²)
// Driven by MAP (pressure hypertrophy) + aldosterone (direct MR-mediated LVH)
// Reduces slowly with effective treatment
// =====================================================================
double LVMi_target = LVMi_normal + kMAP_LVMi * (MAP_c - MAP_basal) + kAldo_LVMi * (Aldo_safe - 8.0) * MR_block_total * (-1.0) + kAldo_LVMi * (Aldo_safe - 8.0);
dxdt_LVMi_c = kout_LVMi * (LVMi_target - LVMi_c);

// =====================================================================
// APA activity (relative, 1=baseline, 0=post-surgery)
// Surgery sets to 0; slowly grows without treatment (autonomous)
// =====================================================================
double APA_reduction = (surgery > 0.5) ? surgery_effect * APA_act : 0.0;
dxdt_APA_act = -kout_APA * APA_reduction;

// CYP11B2 activity (tracks aldosterone synthase)
dxdt_CYP11B2_c = kout_CYP11B2 * (APA_act * APA_severity - CYP11B2_c);

// =====================================================================
// Modeled ARR (ng/dL per ng/mL/h)
// Key diagnostic: > 30 = screening positive for PA
// =====================================================================
dxdt_ARR_c = 0.1 * (Aldo_safe / (Renin_safe < 0.1 ? 0.1 : Renin_safe) - ARR_c);

// =====================================================================
// Metabolic proxy (HOMA-like) — hyperaldosteronism affects glucose metabolism
// =====================================================================
dxdt_HOMA_proxy = 0.01 * (K_deficit * 0.5 + (Aldo_safe - 8.0) * 0.1 - HOMA_proxy);

$TABLE
double SBP_model = MAP_c * 1.33; // Approximate SBP from MAP
double DBP_model = MAP_c * 0.85;
double K_status  = K_c;          // For easy access
double Aldo_PAC  = Aldo_c;       // PAC equivalent
double Renin_PRA = Renin_c;      // PRA equivalent

$CAPTURE SBP_model DBP_model K_status Aldo_PAC Renin_PRA
'

## ---------------------------------------------------------------------------
## Compile model
## ---------------------------------------------------------------------------
pa_mod <- mcode("PA_QSP", pa_code)

## ---------------------------------------------------------------------------
## Initial conditions — PA baseline (APA phenotype)
## ---------------------------------------------------------------------------
pa_init_APA <- list(
  C_spiro=0, C_canrenone=0, C_eple=0, C_fine=0, C_acei=0, C_ccb=0,
  Renin_c  = 0.3,    # suppressed renin in PA (PRA < 1.0)
  AngII_c  = 12.0,   # relatively normal Ang II
  Aldo_c   = 28.0,   # elevated PAC (normal: 3-16 ng/dL; PA: >15)
  ENaC_act = 2.2,    # upregulated ENaC
  Na_c     = 142.5,  # high-normal Na+
  K_c      = 3.2,    # hypokalemia (30% of PA patients)
  HCO3_c   = 27.0,   # mild metabolic alkalosis
  Vol_c    = 1.12,   # modestly expanded plasma volume
  MAP_c    = 110.0,  # elevated MAP (SBP ~146 mmHg)
  TPR_c    = 1.3,    # elevated TPR
  GFR_c    = 95.0,   # mildly elevated (hyperfiltration)
  CardFib  = 0.25,   # early cardiac fibrosis (PA excess vs essential HT)
  LVMi_c   = 115.0,  # LVH (normal <100 g/m²)
  APA_act  = 1.0,    # APA present and active
  CYP11B2_c = 1.0,
  ARR_c    = 100.0,  # ARR >> 30 → screening positive
  HOMA_proxy = 1.5
)

pa_init_BAH <- modifyList(pa_init_APA, list(
  Aldo_c  = 20.0,    # BAH: less severe aldosterone excess
  K_c     = 3.6,     # BAH often normokalemic
  MAP_c   = 105.0,
  CardFib = 0.18,
  LVMi_c  = 108.0,
  ARR_c   = 60.0,
  APA_act = 0.0      # No APA in BAH
))

## Normal control
pa_init_normal <- list(
  C_spiro=0, C_canrenone=0, C_eple=0, C_fine=0, C_acei=0, C_ccb=0,
  Renin_c=1.2, AngII_c=8.0, Aldo_c=8.0, ENaC_act=1.0,
  Na_c=140, K_c=4.2, HCO3_c=24, Vol_c=1.0, MAP_c=90, TPR_c=1.0,
  GFR_c=90, CardFib=0.05, LVMi_c=82, APA_act=0, CYP11B2_c=0,
  ARR_c=6.7, HOMA_proxy=0.5
)

## ---------------------------------------------------------------------------
## Helper: run simulation
## ---------------------------------------------------------------------------
run_pa <- function(init, params, dur_days = 365, label = "unnamed") {
  pa_mod %>%
    init(init) %>%
    param(params) %>%
    mrgsim(end = dur_days * 24, delta = 12) %>%
    as.data.frame() %>%
    mutate(Treatment = label, time_days = time / 24)
}

## ---------------------------------------------------------------------------
## Scenario 1: Untreated APA — 24-month natural history
## ---------------------------------------------------------------------------
sim1 <- run_pa(pa_init_APA,
               list(APA_severity=1.0, BAH_severity=0.0, surgery=0),
               dur_days=730, label="Untreated APA")

## ---------------------------------------------------------------------------
## Scenario 2: Adrenalectomy (curative surgery for APA)
## ---------------------------------------------------------------------------
sim2 <- run_pa(pa_init_APA,
               list(APA_severity=1.0, BAH_severity=0.0, surgery=1,
                    surgery_effect=0.95),
               dur_days=730, label="Adrenalectomy")

## ---------------------------------------------------------------------------
## Scenario 3: Spironolactone 100mg/day (BAH first-line)
## ---------------------------------------------------------------------------
sim3 <- run_pa(pa_init_BAH,
               list(APA_severity=0.0, BAH_severity=0.5, surgery=0,
                    dose_spiro=100),
               dur_days=730, label="Spironolactone 100mg")

## ---------------------------------------------------------------------------
## Scenario 4: Eplerenone 50mg BID (selective MR antagonist)
## ---------------------------------------------------------------------------
sim4 <- run_pa(pa_init_BAH,
               list(APA_severity=0.0, BAH_severity=0.5, surgery=0,
                    dose_eple=100),
               dur_days=730, label="Eplerenone 100mg/d")

## ---------------------------------------------------------------------------
## Scenario 5: Finerenone 20mg/day (non-steroidal, cardio-renal protective)
## ---------------------------------------------------------------------------
sim5 <- run_pa(pa_init_BAH,
               list(APA_severity=0.0, BAH_severity=0.5, surgery=0,
                    dose_fine=20),
               dur_days=730, label="Finerenone 20mg")

## ---------------------------------------------------------------------------
## Scenario 6: Spironolactone + CCB (combination for resistant hypertension)
## ---------------------------------------------------------------------------
sim6 <- run_pa(pa_init_APA,
               list(APA_severity=1.0, BAH_severity=0.0, surgery=0,
                    dose_spiro=100, dose_ccb=10),
               dur_days=730, label="Spiro+CCB (No Surgery)")

## ---------------------------------------------------------------------------
## Scenario 7: Normal control (reference)
## ---------------------------------------------------------------------------
sim7 <- run_pa(pa_init_normal,
               list(APA_severity=0.0, BAH_severity=0.0, surgery=0),
               dur_days=730, label="Normal Control")

## ---------------------------------------------------------------------------
## Scenario 8: ACEi (ramipril 10mg/day) — limited benefit in PA
## ---------------------------------------------------------------------------
sim8 <- run_pa(pa_init_BAH,
               list(APA_severity=0.0, BAH_severity=0.5, surgery=0,
                    dose_acei=10),
               dur_days=730, label="ACEi (Ramipril 10mg)")

## ---------------------------------------------------------------------------
## Combine all
## ---------------------------------------------------------------------------
all_sims <- bind_rows(sim1, sim2, sim3, sim4, sim5, sim6, sim7, sim8) %>%
  mutate(Treatment = factor(Treatment, levels = c(
    "Normal Control", "Untreated APA", "Adrenalectomy",
    "Spironolactone 100mg", "Eplerenone 100mg/d", "Finerenone 20mg",
    "Spiro+CCB (No Surgery)", "ACEi (Ramipril 10mg)"
  )))

pa_colors <- c(
  "Normal Control"          = "#2ECC71",
  "Untreated APA"           = "#E74C3C",
  "Adrenalectomy"           = "#1A5276",
  "Spironolactone 100mg"    = "#3498DB",
  "Eplerenone 100mg/d"      = "#9B59B6",
  "Finerenone 20mg"         = "#E91E63",
  "Spiro+CCB (No Surgery)"  = "#F39C12",
  "ACEi (Ramipril 10mg)"    = "#7F8C8D"
)

## ---------------------------------------------------------------------------
## Plots
## ---------------------------------------------------------------------------
plot_pa <- function(data, yvar, ylabel, title_str) {
  ggplot(data, aes(time_days, .data[[yvar]], color = Treatment)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = pa_colors) +
    labs(x = "Time (days)", y = ylabel, title = title_str) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right", plot.title = element_text(face = "bold"))
}

p1 <- plot_pa(all_sims, "Aldo_c",   "Aldosterone (ng/dL)", "Plasma Aldosterone (PAC)")
p2 <- plot_pa(all_sims, "Renin_c",  "PRA (ng/mL/h)",       "Renin Activity (PRA)")
p3 <- plot_pa(all_sims, "ARR_c",    "ARR",                  "Aldosterone-to-Renin Ratio")
p4 <- plot_pa(all_sims, "SBP_model","SBP (mmHg)",           "Systolic Blood Pressure")
p5 <- plot_pa(all_sims, "K_status", "Serum K+ (mEq/L)",     "Plasma Potassium")
p6 <- plot_pa(all_sims, "HCO3_c",  "HCO3- (mEq/L)",        "Bicarbonate (Alkalosis)")
p7 <- plot_pa(all_sims, "LVMi_c",  "LVMi (g/m²)",          "LV Mass Index")
p8 <- plot_pa(all_sims, "CardFib", "Fibrosis Index (0–1)",  "Cardiac Fibrosis")
p9 <- plot_pa(all_sims, "GFR_c",   "GFR (mL/min/1.73m²)",  "GFR")

## ---------------------------------------------------------------------------
## 12-month endpoint summary
## ---------------------------------------------------------------------------
summary_12mo <- all_sims %>%
  filter(abs(time_days - 365) < 0.6) %>%
  group_by(Treatment) %>%
  slice(1) %>%
  ungroup() %>%
  select(Treatment, Aldo_PAC, Renin_PRA, ARR_c, SBP_model, K_status,
         HCO3_c, LVMi_c, CardFib, GFR_c) %>%
  rename(
    `PAC (ng/dL)` = Aldo_PAC, `PRA (ng/mL/h)` = Renin_PRA, ARR = ARR_c,
    `SBP (mmHg)` = SBP_model, `K+ (mEq/L)` = K_status,
    `HCO3- (mEq/L)` = HCO3_c, `LVMi (g/m²)` = LVMi_c,
    `Cardiac Fibrosis` = CardFib, `GFR` = GFR_c
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

cat("\n===== PA QSP Model — 12-Month Endpoint Summary =====\n")
print(as.data.frame(summary_12mo))

## ---------------------------------------------------------------------------
## Dose-response: Spironolactone
## ---------------------------------------------------------------------------
spiro_doses <- c(12.5, 25, 50, 100, 200, 400)
spiro_dr <- do.call(bind_rows, lapply(spiro_doses, function(d) {
  run_pa(pa_init_BAH, list(APA_severity=0, BAH_severity=0.5, dose_spiro=d),
         dur_days=365, label=paste("Spiro", d, "mg")) %>%
    filter(abs(time_days - 365) < 0.6) %>%
    slice(1) %>%
    mutate(Dose_mg = d)
}))

dr_plot <- ggplot(spiro_dr, aes(Dose_mg, SBP_model)) +
  geom_line(color = "#3498DB", linewidth = 1.2) +
  geom_point(size = 3, color = "#3498DB") +
  labs(x = "Spironolactone Dose (mg/day)", y = "SBP at 12 months (mmHg)",
       title = "Spironolactone Dose–Response in PA (BAH)") +
  theme_bw(base_size = 12)

cat("\nSpironolactone dose-response:\n")
print(spiro_dr %>% select(Dose_mg, SBP_model, K_status, ARR_c) %>%
      mutate(across(where(is.numeric), ~round(.x, 2))))

message("PA QSP model loaded. All scenario simulations available in all_sims.")
