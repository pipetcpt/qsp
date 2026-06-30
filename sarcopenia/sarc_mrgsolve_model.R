################################################################################
##  Sarcopenia QSP Model — mrgsolve implementation
##
##  Scope:
##    - Age-related skeletal muscle wasting (primary sarcopenia)
##    - MPS (mTORC1 / Akt) vs MPB (UPS + autophagy + myostatin/SMAD)
##    - Drug PK/PD: Bimagrumab (ActRII mAb), Apitegromab (latent MSTN mAb),
##      Testosterone IM, Vitamin D3 + 25(OH)D pool, Leucine/HMB, Creatine,
##      Resistance exercise (covariate driver), Anamorelin (ghrelin/GHS-R).
##    - Aging covariates: AGE, SEX, BMI, baseline activity, frailty.
##    - Outputs: Appendicular Lean Mass (ALM, kg), ALMI (kg/m²), Grip strength
##      (kg), Gait speed (m/s), SPPB (0-12), and a frailty hazard.
##
##  Calibration anchors (informational):
##    - Bimagrumab 30 mg/kg IV q4w: BMR-101 (Rooks 2017 J Frailty Aging) ALM
##      +7.1% @ 24 wk vs +0.5% placebo, gait +0.14 m/s.
##    - Apitegromab 20 mg/kg IV q4w: TOPAZ (Day 2021) Type II SMA HFMSE +5.3
##      points; applied here as proof-of-mechanism for latent-MSTN inhibition.
##    - Testosterone 100 mg/wk IM: Bhasin 1996 NEJM ΔLBM +6.1% @10 wk.
##    - Vitamin D 800 IU/d: Bischoff-Ferrari 2009 BMJ falls RR 0.81.
##    - Resistance training: Peterson 2011 MSSE +1.1 kg lean / 12 wk in elderly.
##    - Natural ALM decline: Janssen 2000 J Appl Physiol -0.46 %/yr after age 50.
##
##  Usage:
##    library(mrgsolve); library(tidyverse)
##    mod <- mread("sarc_mrgsolve_model.R", project = "sarcopenia")
##    out <- mod %>% ev(amt = 30, ii = 28, addl = 5, cmt = "BIMA_DEPOT") %>%
##                   mrgsim(end = 168, delta = 0.5)
################################################################################

[ PROB ]
Sarcopenia QSP: aging + anabolic/catabolic balance + 5 drug arms.

[ PLUGIN ] Rcpp

[ PARAM ] @annotated
// === Patient / aging covariates ===
AGE          : 75   : Chronological age (yr)
SEX          : 0    : Sex (0=Male, 1=Female)
WT           : 70   : Body weight (kg)
HT_M         : 1.70 : Height (m)
ACTIVITY     : 0.5  : Activity level (0=bedrest, 1=active elderly)
BASE_ALM     : 18.0 : Baseline appendicular lean mass (kg)
BASE_GRIP    : 25.0 : Baseline grip strength (kg)
BASE_GAIT    : 0.95 : Baseline gait speed (m/s)
BASE_SPPB    : 9    : Baseline SPPB score (0-12)
BASE_IGF1    : 130  : Baseline serum IGF-1 (ng/mL)
BASE_VITD    : 22   : Baseline 25(OH)D (ng/mL)
BASE_MSTN    : 6.0  : Baseline serum myostatin (ng/mL)
BASE_GDF15   : 1200 : Baseline serum GDF-15 (pg/mL)
BASE_IL6     : 3.5  : Baseline IL-6 (pg/mL)

// === Aging-driven decline (per-year) ===
KAGE_ALM     : 0.0046 : Natural ALM loss fraction per year (0.46%/yr)
KAGE_GRIP    : 0.014  : Grip strength loss /yr after 60
KAGE_IGF1    : 0.012  : IGF-1 decline /yr
KAGE_MSTN    : 0.010  : Myostatin rise /yr
KAGE_INFL    : 0.020  : Inflammaging IL-6 rise /yr

// === Anabolic pathway (mTORC1 / IGF-1 / Akt) ===
KMPS_BASAL   : 0.011  : Basal MPS rate constant (fraction/day, ALM units)
EMAX_LEU     : 0.40   : Max MPS stimulation by leucine
EC50_LEU     : 2.5    : Leucine concentration EC50 (mM)
EMAX_IGF     : 0.55   : Max anabolic effect of IGF-1
EC50_IGF     : 140    : IGF-1 EC50 (ng/mL)
EMAX_TEST    : 0.60   : Max testosterone anabolic effect
EC50_TEST    : 6.0    : Testosterone EC50 (ng/mL serum)
EMAX_EX      : 0.30   : Resistance exercise MPS effect
KEX_HALF     : 0.6    : Exercise saturation half-point
ANABRES_AGE  : 0.0035 : Anabolic resistance gain /yr after 60

// === Catabolic pathway (UPS / autophagy / SMAD) ===
KMPB_BASAL   : 0.0105 : Basal MPB rate (fraction/day)
EMAX_MSTN    : 0.50   : Max catabolic stimulation by myostatin
EC50_MSTN    : 5.0    : Myostatin EC50 (ng/mL)
EMAX_INFL    : 0.45   : Max IL-6/TNF catabolic effect
EC50_INFL    : 5.0    : IL-6 EC50 (pg/mL)
EMAX_INACT   : 0.55   : Disuse-driven MPB
SOCS3_GAIN   : 0.30   : IL-6 → SOCS3-mediated anabolic resistance

// === Bimagrumab (BYM338) — ActRIIA/B antagonist mAb ===
//   IV q4w 30 mg/kg; reported t½ ~21 d, V1 ~3 L
BIMA_F       : 1.0    : IV bioavailability
KA_BIMA      : 0      : (IV; not used)
BIMA_CL      : 0.005  : Bimagrumab clearance (L/h)
BIMA_V1      : 3.0    : Central V (L)
BIMA_V2      : 4.5    : Peripheral V (L)
BIMA_Q       : 0.02   : Inter-compartmental Q (L/h)
EMAX_BIMA    : 0.55   : Max muscle mass effect (relative MPS increase)
EC50_BIMA    : 5.0    : Plasma concentration EC50 (mg/L)
KBIMA_OFF    : 0.001  : Receptor occupancy decay /h

// === Apitegromab (SRK-015) — pro/latent-myostatin mAb ===
//   IV q4w 20 mg/kg; t½ ~30 d, V ~3.5 L
APIT_CL      : 0.0042 : L/h
APIT_V1      : 3.5
APIT_V2      : 5.0
APIT_Q       : 0.018
EMAX_APIT    : 0.55   : Max muscle gain via latent MSTN sequestration
EC50_APIT    : 6.0    : EC50 (mg/L)

// === Testosterone IM (cypionate/enanthate 100 mg/wk) ===
//   1-cpt with first-order release from oil depot
KA_TEST      : 0.012  : 1/h release (t½ ~3 d)
TEST_CL      : 50     : L/h apparent
TEST_V       : 30     : L
TEST_F       : 1.0
TEST_TO_AR   : 1.0    : Direct mapping to anabolic effect

// === Vitamin D3 (calcifediol/cholecalciferol) ===
KA_VITD      : 0.05   : /h absorption
VITD_CL      : 0.02   : Conversion D3 → 25(OH)D (1/d)
VITD_V       : 8      : L
KOFF_25OHD   : 0.012  : 25(OH)D elimination /d (t½ ~3 wk)
EMAX_VITD    : 0.20   : Strength / fall risk effect
EC50_VITD    : 30     : Target 25(OH)D (ng/mL)

// === Leucine / HMB / EAA ===
KA_LEU       : 1.2    : /h absorption
LEU_CL       : 18     : L/h
LEU_V        : 15     : L
LEU_DOSE_MMOL : 30    : nominal per intake

// === Anamorelin (ghrelin / GHS-R agonist) ===
KA_ANA       : 0.45   : /h
ANA_CL       : 17     : L/h
ANA_V        : 80     : L
EMAX_ANA     : 0.35   : GH/IGF-1 stimulation
EC50_ANA     : 25     : ng/mL

// === Resistance exercise (covariate / driver, sessions/wk) ===
RT_SESSIONS  : 0      : Sessions per week (0 = none, 3 = ACSM)
RT_INTENSITY : 0.7    : Fraction of 1RM
RT_DAYS      : 0      : Total RT-adapted days (state)

// === Hazard / clinical translation ===
KFALL_BASE   : 0.0008 : Baseline fall hazard /day
GAIT_HALF    : 0.8    : m/s threshold for sarcopenia
KFRAIL       : 0.0003 : Frailty progression /day

[ CMT ] @annotated
BIMA_CENT    : Bimagrumab central (mg)
BIMA_PERI    : Bimagrumab peripheral (mg)
APIT_CENT    : Apitegromab central (mg)
APIT_PERI    : Apitegromab peripheral (mg)
TEST_DEPOT   : Testosterone oil depot (mg)
TEST_CENT    : Testosterone central (ng·L convention via amt mg)
VITD_GUT     : Vitamin D3 gut (IU)
VITD_25OH    : 25(OH)D pool (IU equivalent)
LEU_GUT      : Leucine gut (mmol)
LEU_CENT     : Leucine central (mmol)
ANA_GUT      : Anamorelin gut (mg)
ANA_CENT     : Anamorelin central (mg)
IGF1         : Serum IGF-1 (ng/mL)  // dynamic
MSTN_TOT     : Total serum myostatin (ng/mL)  // pool
IL6          : Serum IL-6 (pg/mL)
GDF15        : Serum GDF-15 (pg/mL)
ALM          : Appendicular lean mass (kg)
GRIP         : Grip strength (kg)
GAIT         : Gait speed (m/s)
SPPB_S       : SPPB total (0-12)
CUM_FALL     : Cumulative fall events
CUM_FRAIL    : Cumulative frailty index (0-1)

[ MAIN ]
// Initial conditions
IGF1_0     = BASE_IGF1;
MSTN_TOT_0 = BASE_MSTN;
IL6_0      = BASE_IL6;
GDF15_0    = BASE_GDF15;
ALM_0      = BASE_ALM;
GRIP_0     = BASE_GRIP;
GAIT_0     = BASE_GAIT;
SPPB_S_0   = BASE_SPPB;
VITD_25OH_0 = BASE_VITD * VITD_V;   // store as "amount"

// Sex/age adjustments
double sex_factor = (SEX==1 ? 0.85 : 1.0);
double age_excess = (AGE - 50.0 > 0 ? AGE - 50.0 : 0);
double anabres    = ANABRES_AGE * age_excess; // additive resistance

[ ODE ]
// --- Bimagrumab 2-cpt PK ---
double CBIMA   = BIMA_CENT / BIMA_V1;          // mg/L
dxdt_BIMA_CENT = - (BIMA_CL/BIMA_V1) * BIMA_CENT
                 - (BIMA_Q /BIMA_V1) * BIMA_CENT
                 + (BIMA_Q /BIMA_V2) * BIMA_PERI;
dxdt_BIMA_PERI =   (BIMA_Q /BIMA_V1) * BIMA_CENT
                 - (BIMA_Q /BIMA_V2) * BIMA_PERI;

// --- Apitegromab 2-cpt PK ---
double CAPIT   = APIT_CENT / APIT_V1;
dxdt_APIT_CENT = - (APIT_CL/APIT_V1) * APIT_CENT
                 - (APIT_Q /APIT_V1) * APIT_CENT
                 + (APIT_Q /APIT_V2) * APIT_PERI;
dxdt_APIT_PERI =   (APIT_Q /APIT_V1) * APIT_CENT
                 - (APIT_Q /APIT_V2) * APIT_PERI;

// --- Testosterone IM depot ---
double CTEST   = TEST_CENT / TEST_V * 1000.0;  // ng/mL from mg
dxdt_TEST_DEPOT = - KA_TEST * TEST_DEPOT;
dxdt_TEST_CENT  =   KA_TEST * TEST_DEPOT - (TEST_CL/TEST_V) * TEST_CENT;

// --- Vitamin D3 absorption → 25(OH)D ---
dxdt_VITD_GUT  = - KA_VITD * VITD_GUT;
dxdt_VITD_25OH =   KA_VITD * VITD_GUT
                 - KOFF_25OHD/24.0 * VITD_25OH;  // /h
double CVITD = VITD_25OH / VITD_V;                // ng/mL approx

// --- Leucine absorption ---
double CLEU = LEU_CENT / LEU_V;                   // mmol/L
dxdt_LEU_GUT  = - KA_LEU * LEU_GUT;
dxdt_LEU_CENT =   KA_LEU * LEU_GUT - (LEU_CL/LEU_V) * LEU_CENT;

// --- Anamorelin ---
double CANA  = ANA_CENT / ANA_V * 1000.0;         // ng/mL
dxdt_ANA_GUT  = - KA_ANA * ANA_GUT;
dxdt_ANA_CENT =   KA_ANA * ANA_GUT - (ANA_CL/ANA_V) * ANA_CENT;

// --- Endocrine / cytokine slow turnover (per hour) ---
// IGF-1 dynamics: baseline + drug stimulation - age decline
double IGF1_stim = EMAX_ANA * CANA / (EC50_ANA + CANA);
double IGF1_ss   = BASE_IGF1 * (1.0 - KAGE_IGF1*age_excess/100.0)
                              * (1.0 + IGF1_stim);
dxdt_IGF1   = 0.06 * (IGF1_ss - IGF1) / 24.0;     // h^-1 (slow)

// Myostatin total pool: latent included; affected by Apit (sequestration) and Bima (no direct)
double MSTN_ss   = BASE_MSTN * (1.0 + KAGE_MSTN*age_excess/100.0);
double APIT_SEQ  = EMAX_APIT * CAPIT / (EC50_APIT + CAPIT);
double MSTN_FREE = MSTN_TOT * (1.0 - APIT_SEQ);
dxdt_MSTN_TOT    = 0.05 * (MSTN_ss - MSTN_TOT) / 24.0;

// IL-6 (inflammaging)
double IL6_ss   = BASE_IL6 * (1.0 + KAGE_INFL*age_excess/100.0);
dxdt_IL6        = 0.08 * (IL6_ss - IL6) / 24.0;

// GDF-15
double GDF_ss   = BASE_GDF15 * (1.0 + 0.015*age_excess);
dxdt_GDF15      = 0.04 * (GDF_ss - GDF15) / 24.0;

// --- ANABOLIC EFFECTS (fractional, multiplicative) ---
double leu_eff   = EMAX_LEU * CLEU / (EC50_LEU + CLEU);
double igf_eff   = EMAX_IGF * IGF1 / (EC50_IGF + IGF1);
double test_eff  = EMAX_TEST * CTEST / (EC50_TEST + CTEST);
double bima_eff  = EMAX_BIMA * CBIMA / (EC50_BIMA + CBIMA);
double apit_eff  = EMAX_APIT * MSTN_FREE / (MSTN_FREE + EC50_MSTN);  // surrogate
double ex_eff    = EMAX_EX  * RT_SESSIONS / (RT_SESSIONS + KEX_HALF);
double vitd_eff  = EMAX_VITD * CVITD / (EC50_VITD + CVITD);

// Anabolic resistance attenuates additive anabolic drive
double ANABDRIVE = (1.0 + leu_eff + igf_eff + test_eff + bima_eff + apit_eff
                       + ex_eff + 0.5*vitd_eff) / (1.0 + anabres + SOCS3_GAIN*(IL6/EC50_INFL));

// --- CATABOLIC EFFECTS ---
double mstn_eff  = EMAX_MSTN * MSTN_FREE / (EC50_MSTN + MSTN_FREE);
double infl_eff  = EMAX_INFL * IL6 / (EC50_INFL + IL6);
double disuse    = EMAX_INACT * (1.0 - ACTIVITY);
double bima_anti = 0.40 * CBIMA / (EC50_BIMA + CBIMA);  // mAb reduces MPB
double CATDRIVE  = (1.0 + mstn_eff + infl_eff + disuse) * (1.0 - bima_anti);

// --- ALM dynamics (kg) ---
double MPS = KMPS_BASAL/24.0 * ALM * ANABDRIVE;
double MPB = KMPB_BASAL/24.0 * ALM * CATDRIVE;
double AGE_LOSS = KAGE_ALM/24.0/365.0 * ALM * age_excess/10.0;
dxdt_ALM   = MPS - MPB - AGE_LOSS;

// --- Grip strength (kg) ---
double GRIP_target = BASE_GRIP * (ALM/BASE_ALM) * (1.0 + vitd_eff*0.5 + test_eff*0.4);
dxdt_GRIP  = 0.04 * (GRIP_target - GRIP) / 24.0;

// --- Gait speed (m/s) ---
double GAIT_target = BASE_GAIT * pow(ALM/BASE_ALM, 0.6) * (1.0 + ex_eff*0.3 + bima_eff*0.2);
dxdt_GAIT  = 0.03 * (GAIT_target - GAIT) / 24.0;

// --- SPPB (composite 0-12) ---
double SPPB_target = 12 - 12 * (1.0 - GAIT/1.2) - (BASE_SPPB<10 ? 1 : 0);
if (SPPB_target < 0) SPPB_target = 0;
if (SPPB_target > 12) SPPB_target = 12;
dxdt_SPPB_S = 0.02 * (SPPB_target - SPPB_S) / 24.0;

// --- Falls hazard ---
double fall_rate = KFALL_BASE *
   pow((GAIT_HALF / (GAIT + 1e-3)), 2.0) *
   (1.0 + 0.5*(1 - vitd_eff));
dxdt_CUM_FALL = fall_rate;

// --- Frailty index ---
double frail_rate = KFRAIL * (1.0 - ANABDRIVE/2.5) * (1.0 + infl_eff);
if (frail_rate < 0) frail_rate = 0;
dxdt_CUM_FRAIL = frail_rate;

[ TABLE ]
double ALMI = ALM / (HT_M*HT_M);
double SARC_DX = (ALMI < (SEX==1 ? 6.0 : 7.0) && GRIP < (SEX==1 ? 16 : 27)) ? 1 : 0;
double SEVERE  = (SARC_DX==1 && GAIT < 0.8) ? 1 : 0;

[ CAPTURE ] @annotated
CBIMA    : Bimagrumab concentration (mg/L)
CAPIT    : Apitegromab concentration (mg/L)
CTEST    : Testosterone (ng/mL)
CVITD    : 25(OH)D (ng/mL)
CLEU     : Leucine plasma (mmol/L)
ANABDRIVE: Net anabolic drive (dimensionless)
CATDRIVE : Net catabolic drive
ALMI     : ALM Index (kg/m²)
SARC_DX  : Sarcopenia diagnosis (EWGSOP2)
SEVERE   : Severe sarcopenia flag

################################################################################
##  Five+ treatment scenarios (driver script — pseudocode)
##
##    library(mrgsolve); library(tibble); library(dplyr); library(ggplot2)
##    mod <- mread("sarc_mrgsolve_model.R", project = "sarcopenia")
##    sim <- function(label, evt, end = 24*7*52) {
##      mod %>% ev(evt) %>% mrgsim(end = end, delta = 1) %>%
##        as_tibble() %>% mutate(arm = label)
##    }
##    arm1 <- sim("No treatment", ev(amt = 0, cmt = "LEU_GUT", time = 0))
##    arm2 <- sim("Bimagrumab 30 mg/kg IV q4w",
##                ev(amt = 30*70, ii = 24*28, addl = 12, cmt = "BIMA_CENT"))
##    arm3 <- sim("Apitegromab 20 mg/kg IV q4w",
##                ev(amt = 20*70, ii = 24*28, addl = 12, cmt = "APIT_CENT"))
##    arm4 <- sim("Testosterone 100 mg IM q1w",
##                ev(amt = 100, ii = 24*7,  addl = 51, cmt = "TEST_DEPOT"))
##    arm5 <- sim("Vit D + EAA + RT 3x/wk",
##                ev(amt = 2000, ii = 24,   addl = 364, cmt = "VITD_GUT") %>%
##                  mutate(RT_SESSIONS = 3))
##    arm6 <- sim("Combo: Bima + RT + EAA", ...)
##
##  Calibration / Sensitivity:
##    - Vary EC50_MSTN over 1-15 ng/mL to bracket assay heterogeneity.
##    - Test SOCS3_GAIN 0-1 for anabolic resistance under inflammaging.
##    - Run age = 65/75/85 cohorts; sex contrast for grip/ALMI cutoffs.
##
################################################################################
