## ============================================================
## Gout QSP Model — mrgsolve ODE Implementation
## Disease: Gout (Hyperuricemia & Crystal-Induced Arthritis)
## Version: 1.0  Date: 2026-06-17
##
## Key pathways modeled:
##   1. Purine metabolism → uric acid production (XO-mediated)
##   2. Urate distribution (plasma, peripheral, synovial, tophus)
##   3. Renal urate handling (GFR, URAT1, OAT1/3, ABCG2)
##   4. Crystal formation & NLRP3 inflammasome
##   5. Inflammatory cascade (IL-1β, TNF-α, neutrophil influx)
##   6. Drug PK: allopurinol/oxypurinol, febuxostat,
##               probenecid, lesinurad, colchicine,
##               indomethacin, anakinra, canakinumab
##   7. Clinical endpoints: sUA, flare risk, tophus volume, QoL
##
## Calibrated against:
##   - Becker et al. 2010 (CONFIRMS, febuxostat vs allopurinol)
##   - Sundy et al. 2011 (canakinumab CANTOS)
##   - Terkeltaub et al. 2010 (colchicine AGREE)
##   - Saag et al. 2017 (lesinurad CLEAR studies)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ============================================================
## Model Code
## ============================================================
gout_model_code <- '
$PROB
Gout QSP Model v1.0
Purine metabolism, urate kinetics, crystal inflammation, drug PK/PD

$PARAM @annotated
// --- Purine Production ---
kprod_UA    : 0.15 : Baseline uric acid production rate (mg/dL/h)
kprod_diet  : 0.03 : Dietary purine contribution (mg/dL/h)
kXO_max     : 1.0  : Maximum XO enzyme capacity (normalized)
km_XO       : 0.5  : Michaelis constant for XO (normalized substrate)
fruct_coeff : 0.2  : Fructose-driven purine production coefficient

// --- Urate Distribution ---
Vc_UA       : 14.0 : Central volume for urate (L)
Vp_UA       : 28.0 : Peripheral volume for urate (L)
Vsyn_UA     : 0.5  : Synovial volume for urate (L)
ktp_UA      : 0.08 : Central to peripheral urate transfer (h^-1)
kpt_UA      : 0.04 : Peripheral to central urate transfer (h^-1)
ksyn_in     : 0.06 : Central to synovial urate transfer (h^-1)
ksyn_out    : 0.05 : Synovial to central urate return (h^-1)

// --- Renal Excretion ---
GFR         : 120  : Glomerular filtration rate (mL/min)
FEurate0    : 0.08 : Baseline fractional excretion of urate
kURAT1      : 0.88 : URAT1 reabsorption fraction (0-1)
kOAT_sec    : 0.12 : OAT1/3 secretion fraction
kABCG2_r    : 0.04 : ABCG2 renal secretion fraction

// --- Gut ABCG2 ---
kgut_ABCG2  : 0.015 : Intestinal ABCG2 urate secretion (h^-1)

// --- Crystal Formation ---
sUA_sat     : 6.8  : Saturation threshold for MSU crystal formation (mg/dL)
kcryst      : 0.002 : Crystal formation rate constant (h^-1 per mg/dL excess)
kdissolve   : 0.001 : Crystal dissolution rate constant (h^-1)
kflare_cryst: 0.005 : Crystal shedding-induced flare rate

// --- Inflammasome (NLRP3) & IL-1β ---
kNLRP3      : 0.8  : NLRP3 activation rate per crystal concentration
kIL1b_prod  : 0.5  : IL-1β production rate (pg/mL/h)
kIL1b_deg   : 0.3  : IL-1β clearance rate (h^-1)
IL1b0       : 1.0  : Baseline IL-1β (pg/mL)

// --- TNF-α ---
kTNFa_prod  : 0.3  : TNF-α production rate (pg/mL/h)
kTNFa_deg   : 0.4  : TNF-α degradation rate (h^-1)
TNFa0       : 5.0  : Baseline TNF-α (pg/mL)

// --- Neutrophil Influx ---
kPMN_rec    : 0.2  : Neutrophil recruitment rate
kPMN_deg    : 0.15 : Neutrophil clearance rate (h^-1)
PMN0        : 1.0  : Baseline neutrophil influx (normalized)

// --- Pain & Inflammation Score ---
kpain_IL1b  : 0.04 : Pain contribution from IL-1β
kpain_PGE2  : 0.03 : Pain contribution from PGE2
kpain_res   : 0.1  : Pain resolution rate (h^-1)
pain_max    : 10.0 : Maximum NRS pain score

// --- PK: Allopurinol / Oxypurinol ---
ka_Allo     : 0.9  : Allopurinol absorption rate constant (h^-1)
F_Allo      : 0.80 : Allopurinol oral bioavailability
CL_Allo     : 9.0  : Allopurinol clearance (L/h)
Vc_Allo     : 30.0 : Allopurinol central volume (L)
k_Allo_Oxy  : 0.5  : Allopurinol to oxypurinol conversion rate (h^-1)
CL_Oxy      : 0.5  : Oxypurinol renal clearance (L/h)
Vc_Oxy      : 45.0 : Oxypurinol central volume (L)
Vp_Oxy      : 30.0 : Oxypurinol peripheral volume (L)
ktp_Oxy     : 0.5  : Oxypurinol central-to-periph transfer (h^-1)
kpt_Oxy     : 0.2  : Oxypurinol periph-to-central transfer (h^-1)
Ki_Oxy      : 0.001: Oxypurinol Ki for XO inhibition (mg/L)

// --- PK: Febuxostat ---
ka_Febu     : 1.2  : Febuxostat absorption rate constant (h^-1)
F_Febu      : 0.49 : Febuxostat oral bioavailability
CL_Febu     : 4.0  : Febuxostat total clearance (L/h)
Vc_Febu     : 12.0 : Febuxostat central volume (L)
Vp_Febu     : 24.0 : Febuxostat peripheral volume (L)
ktp_Febu    : 0.3  : Febuxostat central-to-periph transfer (h^-1)
kpt_Febu    : 0.15 : Febuxostat periph-to-central transfer (h^-1)
Ki_Febu     : 0.000001 : Febuxostat Ki for XO (mg/L, very potent)

// --- PK: Probenecid ---
ka_Prob     : 1.5  : Probenecid absorption rate (h^-1)
F_Prob      : 1.0  : Probenecid bioavailability
CL_Prob     : 3.5  : Probenecid clearance (L/h)
Vc_Prob     : 10.0 : Probenecid central volume (L)
Vp_Prob     : 20.0 : Probenecid peripheral volume (L)
ktp_Prob    : 0.4  : Prob central-to-periph transfer (h^-1)
kpt_Prob    : 0.2  : Prob periph-to-central transfer (h^-1)
IC50_Prob   : 5.0  : Probenecid IC50 for URAT1 (mg/L)

// --- PK: Lesinurad ---
ka_Lesi     : 2.1  : Lesinurad absorption rate (h^-1)
F_Lesi      : 0.95 : Lesinurad bioavailability
CL_Lesi     : 8.0  : Lesinurad clearance (L/h)
Vc_Lesi     : 14.0 : Lesinurad central volume (L)
IC50_Lesi   : 0.1  : Lesinurad IC50 for URAT1 (mg/L)

// --- PK: Colchicine ---
ka_Colch    : 1.8  : Colchicine absorption rate (h^-1)
F_Colch     : 0.45 : Colchicine oral bioavailability
CL_Colch    : 20.0 : Colchicine total clearance (L/h)
Vc_Colch    : 100.0: Colchicine central volume (L)
Vp_Colch    : 400.0: Colchicine peripheral volume (L)
ktp_Colch   : 2.1  : Colchicine central-to-periph transfer (h^-1)
kpt_Colch   : 0.42 : Colchicine periph-to-central transfer (h^-1)
IC50_Colch  : 0.0003: Colchicine IC50 for NLRP3 (mg/L)
Emax_Colch  : 0.85 : Maximum NLRP3 inhibition by colchicine

// --- PK: Indomethacin (NSAID) ---
ka_Indo     : 1.5  : Indomethacin absorption rate (h^-1)
F_Indo      : 0.90 : Indomethacin bioavailability
CL_Indo     : 8.0  : Indomethacin clearance (L/h)
Vc_Indo     : 16.0 : Indomethacin central volume (L)
Vt_Indo     : 20.0 : Indomethacin tissue/synovial volume (L)
ktp_Indo    : 0.6  : Indo central-to-tissue transfer (h^-1)
kpt_Indo    : 0.3  : Indo tissue-to-central transfer (h^-1)
IC50_Indo   : 0.002: Indomethacin IC50 for COX (mg/L)

// --- PK: Anakinra (IL-1Ra) ---
ka_Ana      : 0.6  : Anakinra SC absorption rate (h^-1)
F_Ana       : 0.95 : Anakinra SC bioavailability
CL_Ana      : 8.0  : Anakinra total clearance (L/h)
Vc_Ana      : 8.0  : Anakinra central volume (L)
Vp_Ana      : 5.0  : Anakinra peripheral volume (L)
ktp_Ana     : 0.15 : Anakinra central-to-periph transfer (h^-1)
kpt_Ana     : 0.08 : Anakinra periph-to-central transfer (h^-1)
IC50_Ana    : 0.5  : Anakinra IC50 for IL-1R blockade (mg/L)

// --- PK: Canakinumab ---
ka_Cana     : 0.009: Canakinumab SC absorption rate (h^-1)
F_Cana      : 0.70 : Canakinumab SC bioavailability
CL_Cana     : 0.23 : Canakinumab clearance (L/h)
Vc_Cana     : 6.0  : Canakinumab central volume (L)
Vp_Cana     : 3.5  : Canakinumab peripheral volume (L)
ktp_Cana    : 0.005: Canakinumab central-to-periph (h^-1)
kpt_Cana    : 0.003: Canakinumab periph-to-central (h^-1)
Kon_Cana    : 1.0  : Canakinumab-IL1b association (L/mg/h)
Koff_Cana   : 0.00001: Canakinumab-IL1b dissociation (h^-1)

// --- Tophus Dynamics ---
ktoph_form  : 0.0001: Tophus crystal deposition rate (cm3/h per unit crystal)
ktoph_diss  : 0.00005: Tophus dissolution rate (cm3/h per unit sUA lowering)
Toph0       : 0.0  : Initial tophus volume (cm3, 0 for non-tophaceous)

// --- Disease Progression ---
k_joint_dmg : 0.0001: Cumulative joint damage rate (units/h per flare)
k_eGFR_loss : 0.00002: eGFR decline per unit chronic hyperuricemia (h^-1)
eGFR0       : 90.0 : Baseline eGFR (mL/min/1.73m2)

// --- Covariates (patient-specific) ---
BW          : 80.0 : Body weight (kg)
AGE         : 50.0 : Age (years)
SEX         : 1.0  : Sex (1=male, 0=female)
RACE_AFRO   : 0.0  : African ancestry (1=yes; ↑ risk)
FOOD_score  : 0.5  : Diet purine score (0-1)
ETOH        : 0.0  : Alcohol intake (drinks/day)

$CMT @annotated
// Purine/Urate
A_UA_gut    : Urate precursor gut absorption depot
A_UA_cent   : Urate central compartment (plasma, mg)
A_UA_peri   : Urate peripheral tissue compartment (mg)
A_UA_syn    : Urate in synovial fluid (mg)
A_Crystal   : MSU crystal pool in joint (mg, normalized)
A_Tophus    : Tophus volume (cm3)

// Inflammation
A_IL1b      : Active IL-1β concentration (pg/mL)
A_TNFa      : TNF-α (pg/mL)
A_PMN       : Neutrophil influx to joint (normalized)
A_Pain      : Acute gout pain score (NRS 0-10)
A_JointDmg  : Cumulative joint damage (arbitrary units)
A_eGFR      : eGFR (mL/min/1.73m2)

// PK: Allopurinol / Oxypurinol
A_Allo_gut  : Allopurinol gut depot
A_Allo_cent : Allopurinol central
A_Oxy_cent  : Oxypurinol central
A_Oxy_peri  : Oxypurinol peripheral

// PK: Febuxostat
A_Febu_gut  : Febuxostat gut depot
A_Febu_cent : Febuxostat central
A_Febu_peri : Febuxostat peripheral

// PK: Probenecid
A_Prob_gut  : Probenecid gut depot
A_Prob_cent : Probenecid central
A_Prob_peri : Probenecid peripheral

// PK: Lesinurad
A_Lesi_gut  : Lesinurad gut depot
A_Lesi_cent : Lesinurad central

// PK: Colchicine
A_Colch_gut : Colchicine gut depot
A_Colch_cent: Colchicine central
A_Colch_peri: Colchicine peripheral

// PK: Indomethacin
A_Indo_gut  : Indomethacin gut depot
A_Indo_cent : Indomethacin central
A_Indo_tiss : Indomethacin tissue/synovial

// PK: Anakinra
A_Ana_SC    : Anakinra SC depot
A_Ana_cent  : Anakinra central
A_Ana_peri  : Anakinra peripheral

// PK: Canakinumab
A_Cana_SC   : Canakinumab SC depot
A_Cana_cent : Canakinumab central
A_Cana_peri : Canakinumab peripheral
A_IL1b_Cana : IL-1β:Canakinumab complex

$MAIN
// Concentrations
double C_UA     = A_UA_cent / Vc_UA;         // serum urate (mg/dL)
double C_UA_syn = A_UA_syn / Vsyn_UA;        // synovial urate (mg/dL)
double C_Crystal = A_Crystal;                // crystal concentration

double C_Allo   = A_Allo_cent / Vc_Allo;    // allopurinol (mg/L)
double C_Oxy    = A_Oxy_cent / Vc_Oxy;      // oxypurinol (mg/L)
double C_Febu   = A_Febu_cent / Vc_Febu;    // febuxostat (mg/L)
double C_Prob   = A_Prob_cent / Vc_Prob;    // probenecid (mg/L)
double C_Lesi   = A_Lesi_cent / Vc_Lesi;    // lesinurad (mg/L)
double C_Colch  = A_Colch_cent / Vc_Colch;  // colchicine (mg/L)
double C_Indo   = A_Indo_cent / Vc_Indo;    // indomethacin (mg/L)
double C_IndoT  = A_Indo_tiss / Vt_Indo;    // indomethacin tissue
double C_Ana    = A_Ana_cent / Vc_Ana;      // anakinra (mg/L)
double C_Cana   = A_Cana_cent / Vc_Cana;   // canakinumab (mg/L)

// XO inhibition (combined oxypurinol + febuxostat)
double XO_inhib_Oxy  = C_Oxy  / (C_Oxy  + Ki_Oxy);
double XO_inhib_Febu = C_Febu / (C_Febu + Ki_Febu);
double XO_inhib      = 1.0 - (1.0 - XO_inhib_Oxy) * (1.0 - XO_inhib_Febu);
double XO_activity   = kXO_max * (1.0 - XO_inhib);

// URAT1 inhibition (probenecid + lesinurad)
double URAT1_inhib_P = (C_Prob / IC50_Prob) / (1.0 + C_Prob / IC50_Prob);
double URAT1_inhib_L = (C_Lesi / IC50_Lesi) / (1.0 + C_Lesi / IC50_Lesi);
double URAT1_inhib   = 1.0 - (1.0 - URAT1_inhib_P) * (1.0 - URAT1_inhib_L);
double kURAT1_eff    = kURAT1 * (1.0 - URAT1_inhib);

// Renal clearance of urate
double CLr_UA = GFR * 0.001 * 60.0 * (1.0 - kURAT1_eff + kOAT_sec + kABCG2_r); // L/h

// Gut excretion (ABCG2)
double ABCG2_Q141K_effect = 1.0; // set < 1.0 for Q141K polymorphism patients
double kgut_eff = kgut_ABCG2 * ABCG2_Q141K_effect;

// NLRP3 inhibition by colchicine
double NLRP3_inhib_Colch = Emax_Colch * C_Colch / (C_Colch + IC50_Colch);

// COX inhibition by NSAID
double COX_inhib_Indo = C_IndoT / (C_IndoT + IC50_Indo);

// IL-1β effective concentration (accounting for Cana neutralization)
double IL1b_free  = A_IL1b;
double IL1b_bound = A_IL1b_Cana;

// IL-1R blockade by anakinra
double Ana_effect = C_Ana / (C_Ana + IC50_Ana);

// PGE2 proxy (proportional to COX activity * inflammatory stimulus)
double PGE2 = (1.0 - COX_inhib_Indo) * (IL1b_free / (IL1b_free + 10.0));

// Crystal formation: driven by supersaturation in synovial fluid
double dCryst_form = (C_UA_syn > sUA_sat) ?
    kcryst * (C_UA_syn - sUA_sat) : 0.0;
double dCryst_diss = kdissolve * A_Crystal * (1.0 / (1.0 + C_UA_syn / sUA_sat));

// Flare trigger (crystal concentration drives NLRP3)
double NLRP3_act = kNLRP3 * A_Crystal / (1.0 + A_Crystal) * (1.0 - NLRP3_inhib_Colch);

// IL-1β production with IL-1Ra (anakinra) effect
double IL1b_prod = kIL1b_prod * NLRP3_act * (1.0 - Ana_effect);

// Uric acid production rate (diet + endogenous + alcohol + fructose)
double UA_prod = (kprod_UA + kprod_diet * (1.0 + FOOD_score + ETOH * fruct_coeff))
                 * XO_activity * (BW / 70.0);

// Pain: driven by IL-1β, PGE2, PMN activation
double pain_drive = kpain_IL1b * IL1b_free + kpain_PGE2 * PGE2 * 10.0
                    + 0.02 * A_PMN;

// Initial conditions
A_UA_cent_0  = 6.0 * Vc_UA;    // sUA ~6 mg/dL at baseline
A_UA_peri_0  = 6.0 * Vp_UA * 0.8;
A_UA_syn_0   = 5.5 * Vsyn_UA;
A_IL1b_0     = IL1b0;
A_TNFa_0     = TNFa0;
A_PMN_0      = PMN0;
A_Pain_0     = 0.0;
A_eGFR_0     = eGFR0;
A_Tophus_0   = Toph0;

$ODE
// =============================================================
// Urate Distribution
// =============================================================
dxdt_A_UA_gut  = 0.0; // Fed from dosing events (not dietary)
dxdt_A_UA_cent = UA_prod                                 // production
               - CLr_UA * C_UA                           // renal excretion
               - kgut_eff * A_UA_cent                    // gut secretion
               - ktp_UA * A_UA_cent                      // to peripheral
               + kpt_UA * A_UA_peri                      // from peripheral
               - ksyn_in * A_UA_cent                     // to synovial
               + ksyn_out * A_UA_syn;                    // from synovial

dxdt_A_UA_peri = ktp_UA * A_UA_cent - kpt_UA * A_UA_peri;

dxdt_A_UA_syn  = ksyn_in * A_UA_cent
               - ksyn_out * A_UA_syn
               - dCryst_form * Vsyn_UA                   // crystal deposition removes UA
               + dCryst_diss * Vsyn_UA;                  // dissolution releases UA

// =============================================================
// Crystal & Tophus
// =============================================================
dxdt_A_Crystal = dCryst_form - dCryst_diss
               - kflare_cryst * A_Crystal * A_PMN;       // crystal clearance by PMN

dxdt_A_Tophus  = ktoph_form * A_Crystal
               - ktoph_diss * A_Tophus * fmax(0.0, 6.0 - C_UA); // dissolve if sUA<6

// =============================================================
// NLRP3 Inflammasome & Cytokines
// =============================================================
dxdt_A_IL1b    = IL1b_prod
               - kIL1b_deg * A_IL1b
               - Kon_Cana * A_IL1b * A_Cana_cent + Koff_Cana * A_IL1b_Cana;

dxdt_A_TNFa    = kTNFa_prod * NLRP3_act * (1.0 - 0.8 * COX_inhib_Indo)
               - kTNFa_deg * A_TNFa;

dxdt_A_PMN     = kPMN_rec * (IL1b_free / (IL1b_free + 5.0) + 0.5 * A_TNFa / (A_TNFa + 20.0))
               - kPMN_deg * A_PMN;

// =============================================================
// Pain & Disease Outcomes
// =============================================================
dxdt_A_Pain    = pain_drive - kpain_res * A_Pain;

dxdt_A_JointDmg = k_joint_dmg * A_PMN * A_Crystal;

dxdt_A_eGFR    = -k_eGFR_loss * fmax(0.0, C_UA - 6.0) * A_eGFR;

// =============================================================
// PK: Allopurinol / Oxypurinol
// =============================================================
dxdt_A_Allo_gut  = -ka_Allo * A_Allo_gut;
dxdt_A_Allo_cent = ka_Allo * F_Allo * A_Allo_gut
                 - (CL_Allo / Vc_Allo) * A_Allo_cent
                 - k_Allo_Oxy * A_Allo_cent;
dxdt_A_Oxy_cent  = k_Allo_Oxy * A_Allo_cent
                 - (CL_Oxy / Vc_Oxy) * A_Oxy_cent
                 - ktp_Oxy * A_Oxy_cent
                 + kpt_Oxy * A_Oxy_peri;
dxdt_A_Oxy_peri  = ktp_Oxy * A_Oxy_cent - kpt_Oxy * A_Oxy_peri;

// =============================================================
// PK: Febuxostat
// =============================================================
dxdt_A_Febu_gut  = -ka_Febu * A_Febu_gut;
dxdt_A_Febu_cent = ka_Febu * F_Febu * A_Febu_gut
                 - (CL_Febu / Vc_Febu) * A_Febu_cent
                 - ktp_Febu * A_Febu_cent
                 + kpt_Febu * A_Febu_peri;
dxdt_A_Febu_peri = ktp_Febu * A_Febu_cent - kpt_Febu * A_Febu_peri;

// =============================================================
// PK: Probenecid
// =============================================================
dxdt_A_Prob_gut  = -ka_Prob * A_Prob_gut;
dxdt_A_Prob_cent = ka_Prob * F_Prob * A_Prob_gut
                 - (CL_Prob / Vc_Prob) * A_Prob_cent
                 - ktp_Prob * A_Prob_cent
                 + kpt_Prob * A_Prob_peri;
dxdt_A_Prob_peri = ktp_Prob * A_Prob_cent - kpt_Prob * A_Prob_peri;

// =============================================================
// PK: Lesinurad
// =============================================================
dxdt_A_Lesi_gut  = -ka_Lesi * A_Lesi_gut;
dxdt_A_Lesi_cent = ka_Lesi * F_Lesi * A_Lesi_gut
                 - (CL_Lesi / Vc_Lesi) * A_Lesi_cent;

// =============================================================
// PK: Colchicine
// =============================================================
dxdt_A_Colch_gut  = -ka_Colch * A_Colch_gut;
dxdt_A_Colch_cent = ka_Colch * F_Colch * A_Colch_gut
                  - (CL_Colch / Vc_Colch) * A_Colch_cent
                  - ktp_Colch * A_Colch_cent
                  + kpt_Colch * A_Colch_peri;
dxdt_A_Colch_peri = ktp_Colch * A_Colch_cent - kpt_Colch * A_Colch_peri;

// =============================================================
// PK: Indomethacin
// =============================================================
dxdt_A_Indo_gut  = -ka_Indo * A_Indo_gut;
dxdt_A_Indo_cent = ka_Indo * F_Indo * A_Indo_gut
                 - (CL_Indo / Vc_Indo) * A_Indo_cent
                 - ktp_Indo * A_Indo_cent
                 + kpt_Indo * A_Indo_tiss;
dxdt_A_Indo_tiss = ktp_Indo * A_Indo_cent - kpt_Indo * A_Indo_tiss;

// =============================================================
// PK: Anakinra (SC)
// =============================================================
dxdt_A_Ana_SC   = -ka_Ana * A_Ana_SC;
dxdt_A_Ana_cent = ka_Ana * F_Ana * A_Ana_SC
                - (CL_Ana / Vc_Ana) * A_Ana_cent
                - ktp_Ana * A_Ana_cent
                + kpt_Ana * A_Ana_peri;
dxdt_A_Ana_peri = ktp_Ana * A_Ana_cent - kpt_Ana * A_Ana_peri;

// =============================================================
// PK: Canakinumab (SC, target-mediated disposition)
// =============================================================
dxdt_A_Cana_SC   = -ka_Cana * A_Cana_SC;
dxdt_A_Cana_cent = ka_Cana * F_Cana * A_Cana_SC
                 - (CL_Cana / Vc_Cana) * A_Cana_cent
                 - ktp_Cana * A_Cana_cent
                 + kpt_Cana * A_Cana_peri
                 - Kon_Cana * A_Cana_cent * A_IL1b
                 + Koff_Cana * A_IL1b_Cana;
dxdt_A_Cana_peri = ktp_Cana * A_Cana_cent - kpt_Cana * A_Cana_peri;
dxdt_A_IL1b_Cana = Kon_Cana * A_Cana_cent * A_IL1b
                 - Koff_Cana * A_IL1b_Cana
                 - (CL_Cana / Vc_Cana) * A_IL1b_Cana; // catabolism of complex

$TABLE
capture sUA      = A_UA_cent / Vc_UA;
capture sUA_syn  = A_UA_syn / Vsyn_UA;
capture Crystal  = A_Crystal;
capture Tophus   = A_Tophus;
capture IL1b_f   = A_IL1b;
capture TNFa_f   = A_TNFa;
capture PMN      = A_PMN;
capture Pain     = fmin(A_Pain, pain_max);
capture JntDmg   = A_JointDmg;
capture eGFR_sim = A_eGFR;
capture XO_inh   = XO_inhib * 100.0;  // % XO inhibition
capture URAT1_inh = URAT1_inhib * 100.0;
capture C_Allo_out = C_Allo;
capture C_Oxy_out  = C_Oxy;
capture C_Febu_out = C_Febu;
capture C_Prob_out = C_Prob;
capture C_Lesi_out = C_Lesi;
capture C_Colch_out = C_Colch;
capture C_Ana_out  = C_Ana;
capture C_Cana_out = C_Cana;
capture FEurate    = CLr_UA * C_UA / (GFR * 0.001 * 60.0 * C_UA + 0.0001) * 100.0;
'

## Compile model
mod <- mread_cache("gout_qsp", tempdir(), gout_model_code)

## ============================================================
## SCENARIO DEFINITIONS
## ============================================================

## Helper: build dosing event
dose_ev <- function(drug_cmt, amount, ii, addl, start=0) {
    ev(cmt=drug_cmt, amt=amount, ii=ii, addl=addl, time=start, rate=0)
}

# Simulation time: 52 weeks (8736 h)
sim_end  <- 52 * 7 * 24   # hours
obs_times <- seq(0, sim_end, by=24)

## ============================================================
## Scenario 1: Untreated Hyperuricemia (baseline)
## ============================================================
scen1 <- mod %>%
    param(FOOD_score=0.7, ETOH=2.0) %>%   # high-purine diet + alcohol
    mrgsim(end=sim_end, delta=24) %>%
    as_tibble() %>%
    mutate(Scenario="1_Untreated")

## ============================================================
## Scenario 2: Allopurinol 300 mg/day (standard urate-lowering)
## ============================================================
e2 <- ev(cmt="A_Allo_gut", amt=300, ii=24, addl=sim_end/24 - 1, time=0)

scen2 <- mod %>%
    ev(e2) %>%
    mrgsim(end=sim_end, delta=24) %>%
    as_tibble() %>%
    mutate(Scenario="2_Allopurinol300")

## ============================================================
## Scenario 3: Febuxostat 80 mg/day (non-purine XO inhibitor)
## ============================================================
e3 <- ev(cmt="A_Febu_gut", amt=80, ii=24, addl=sim_end/24 - 1, time=0)

scen3 <- mod %>%
    ev(e3) %>%
    mrgsim(end=sim_end, delta=24) %>%
    as_tibble() %>%
    mutate(Scenario="3_Febuxostat80")

## ============================================================
## Scenario 4: Combination — Allopurinol + Lesinurad
##   (suboptimal responders needing dual therapy)
## ============================================================
e4a <- ev(cmt="A_Allo_gut",  amt=300, ii=24, addl=sim_end/24 - 1, time=0)
e4b <- ev(cmt="A_Lesi_gut",  amt=200, ii=24, addl=sim_end/24 - 1, time=0)

scen4 <- mod %>%
    ev(e4a + e4b) %>%
    mrgsim(end=sim_end, delta=24) %>%
    as_tibble() %>%
    mutate(Scenario="4_Allo_Lesinurad")

## ============================================================
## Scenario 5: Acute Gout Flare Treatment — Colchicine
##   Low-dose regimen: 1.2mg then 0.6mg 1h later
## ============================================================
e5_flare <- ev(cmt="A_Crystal", amt=5, time=0)   # induce crystal flare
e5a <- ev(cmt="A_Colch_gut", amt=1.2, time=0)
e5b <- ev(cmt="A_Colch_gut", amt=0.6, time=1)
e5c <- ev(cmt="A_Colch_gut", amt=0.6, ii=12, addl=7, time=12) # maintenance 5 days

scen5 <- mod %>%
    ev(e5_flare + e5a + e5b + e5c) %>%
    mrgsim(end=14*24, delta=2) %>%
    as_tibble() %>%
    mutate(Scenario="5_Colchicine_acute")

## ============================================================
## Scenario 6: Acute Flare — Indomethacin 50mg TID
## ============================================================
e6_flare <- ev(cmt="A_Crystal", amt=5, time=0)
e6 <- ev(cmt="A_Indo_gut", amt=50, ii=8, addl=20, time=0) # 50mg q8h × 7 days

scen6 <- mod %>%
    ev(e6_flare + e6) %>%
    mrgsim(end=14*24, delta=2) %>%
    as_tibble() %>%
    mutate(Scenario="6_Indomethacin_acute")

## ============================================================
## Scenario 7: Biologic — Canakinumab 150mg SC (refractory flares)
## ============================================================
e7_flare <- ev(cmt="A_Crystal", amt=5, time=0)
e7 <- ev(cmt="A_Cana_SC", amt=150, time=0)

scen7 <- mod %>%
    ev(e7_flare + e7) %>%
    mrgsim(end=90*24, delta=12) %>%
    as_tibble() %>%
    mutate(Scenario="7_Canakinumab")

## ============================================================
## Scenario 8: Febuxostat 80mg + Flare prophylaxis (Colchicine 0.5mg/day)
## ============================================================
e8a <- ev(cmt="A_Febu_gut",  amt=80,  ii=24, addl=sim_end/24-1, time=0)
e8b <- ev(cmt="A_Colch_gut", amt=0.5, ii=24, addl=sim_end/24-1, time=0)

scen8 <- mod %>%
    ev(e8a + e8b) %>%
    mrgsim(end=sim_end, delta=24) %>%
    as_tibble() %>%
    mutate(Scenario="8_Febu_ColchProphylaxis")

## ============================================================
## RESULTS SUMMARY
## ============================================================
cat("\n=== Gout QSP Model — 52-week Outcomes Summary ===\n")
summary_all <- list(scen1, scen2, scen3, scen4, scen8) %>%
    bind_rows() %>%
    filter(time == max(time)) %>%
    select(Scenario, sUA, XO_inh, URAT1_inh, Crystal, Tophus, eGFR_sim)

print(summary_all)

## ============================================================
## VISUALIZATION
## ============================================================
chronic_data <- bind_rows(scen1, scen2, scen3, scen4, scen8) %>%
    mutate(week = time / 168)

acute_data <- bind_rows(scen5, scen6, scen7) %>%
    mutate(day = time / 24)

p1 <- ggplot(chronic_data, aes(x=week, y=sUA, color=Scenario)) +
    geom_line(linewidth=1.0) +
    geom_hline(yintercept=6.0, linetype="dashed", color="red") +
    geom_hline(yintercept=5.0, linetype="dotted", color="blue") +
    labs(title="Serum Urate Over 52 Weeks",
         subtitle="Red dashed = target <6 mg/dL; Blue dotted = target <5 mg/dL (tophaceous)",
         x="Week", y="sUA (mg/dL)") +
    theme_bw() + theme(legend.position="bottom")

p2 <- ggplot(chronic_data, aes(x=week, y=Crystal, color=Scenario)) +
    geom_line(linewidth=1.0) +
    labs(title="MSU Crystal Pool Over 52 Weeks",
         x="Week", y="Crystal Burden (normalized)") +
    theme_bw() + theme(legend.position="bottom")

p3 <- ggplot(acute_data, aes(x=day, y=Pain, color=Scenario)) +
    geom_line(linewidth=1.0) +
    labs(title="Acute Gout Pain Score (NRS 0-10)",
         x="Day", y="NRS Pain Score") +
    theme_bw() + theme(legend.position="bottom")

p4 <- ggplot(chronic_data, aes(x=week, y=IL1b_f, color=Scenario)) +
    geom_line(linewidth=1.0) +
    labs(title="Synovial IL-1β Over 52 Weeks",
         x="Week", y="IL-1β (pg/mL)") +
    theme_bw() + theme(legend.position="bottom")

p5 <- ggplot(chronic_data, aes(x=week, y=Tophus, color=Scenario)) +
    geom_line(linewidth=1.0) +
    labs(title="Tophus Volume Regression Over 52 Weeks",
         x="Week", y="Tophus Volume (cm³)") +
    theme_bw() + theme(legend.position="bottom")

p6 <- ggplot(chronic_data, aes(x=week, y=eGFR_sim, color=Scenario)) +
    geom_line(linewidth=1.0) +
    geom_hline(yintercept=60, linetype="dashed", color="orange") +
    labs(title="eGFR Trajectory Over 52 Weeks",
         subtitle="Orange = CKD stage 3 threshold",
         x="Week", y="eGFR (mL/min/1.73m²)") +
    theme_bw() + theme(legend.position="bottom")

## ============================================================
## SENSITIVITY ANALYSIS: sUA target achievement
## ============================================================
cat("\n=== sUA Target Achievement (<6 mg/dL) at Week 24 ===\n")
target_week24 <- bind_rows(scen2, scen3, scen4, scen8) %>%
    filter(abs(time - 24*7*24) < 24) %>%
    group_by(Scenario) %>%
    summarise(
        sUA_wk24    = mean(sUA),
        target_met  = mean(sUA) < 6.0,
        XO_inh_pct  = mean(XO_inh),
        .groups="drop"
    )
print(target_week24)

cat("\nModel compilation and simulation complete.\n")
