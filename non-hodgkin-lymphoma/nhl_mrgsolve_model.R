##############################################################################
# Non-Hodgkin Lymphoma (DLBCL) — QSP mrgsolve Model
# Disease: Diffuse Large B-Cell Lymphoma (DLBCL)
# Subtypes: GCB, ABC (Activated B-Cell), Double-Hit Lymphoma (DHL)
#
# Model Structure (22 compartments):
#   PK (10 cmt): Rituximab 2-cmt + TMDD, Cyclophosphamide, Doxorubicin,
#                Venetoclax, Ibrutinib
#   PD (12 cmt): Tumor burden, CD20 surface density, BCR/NF-kB signal,
#                BCL-2 occupancy, Effector cells (NK, CD8+),
#                CRS/toxicity markers, ANC
#
# Treatment Scenarios (6):
#   1. Untreated DLBCL (natural progression)
#   2. R-CHOP x6 cycles (standard first-line)
#   3. R-CHOP + Polatuzumab Vedotin (Pola-R-CHP) — POLARIX regimen
#   4. R-CHOP + Ibrutinib (ABC-DLBCL, PHOENIX trial)
#   5. R2 (Rituximab + Lenalidomide, relapsed/refractory)
#   6. Venetoclax + R-CHOP (BCL-2 high, CAVALLI study)
#
# Key clinical calibration:
#   - R-CHOP 6cy: ORR ~75-85%, CR ~65-70%, 2yr PFS ~60% (Coiffier 2002 NEJM)
#   - Pola-R-CHP vs R-CHOP: PFS 76.7% vs 70.2% at 2yr (Tilly 2022 NEJM)
#   - PHOENIX (ibrutinib+R-CHOP ABC): EFS HR 0.934 (not significant overall,
#     benefit in younger non-GCB) (Younes 2019 NEJM)
#   - CAVALLI venetoclax+R-CHOP: BCL2+ ORR 88% (Morschhauser 2021 JCO)
#
# Units: time in days, concentration in mg/L (µg/mL), doses in mg
# References: see nhl_references.md
##############################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

##############################################################################
# MODEL CODE
##############################################################################
code <- '
$PROB
Non-Hodgkin Lymphoma (DLBCL) QSP Model
Rituximab TMDD + CHOP PK + Tumor Dynamics + Immune Effectors

$PARAM @annotated
// --- Rituximab PK (2-cmt + TMDD) ---
CLR   : 0.008 : Rituximab linear clearance (L/h) [population mean]
VcR   : 3.1   : Rituximab central volume (L)
QR    : 0.12  : Rituximab intercompartmental clearance (L/h)
VpR   : 3.7   : Rituximab peripheral volume (L)
konR  : 0.27  : RTX-CD20 association rate (1/nM/h)
koffR : 0.0003: RTX-CD20 dissociation rate (1/h)
kintR : 0.05  : RTX-CD20 internalization rate (1/h)
ksynCD20: 0.012: CD20 synthesis rate (nM/h per cell)
kdegCD20: 0.004: CD20 baseline degradation rate (1/h)
Bmax_CD20: 0.35: Total CD20 capacity (nM, referenced to central volume)
RTX_dose: 0.0 : Rituximab dose (mg); set per scenario
RTX_MW  : 144000.0: Rituximab molecular weight (g/mol)

// --- Cyclophosphamide PK (1-cmt, active metabolite) ---
F_CPP  : 0.9  : Bioavailability (IV=1 assumed)
kaCPP  : 2.5  : CYP2B6 activation rate constant (1/h) → 4-OH-CPP
kelCPP : 0.18 : Active metabolite elimination (1/h) [T1/2~3.8h]
VdCPP  : 50.0 : Apparent volume 4-OH-CPP (L)
CPP_dose: 0.0 : Cyclophosphamide dose (mg)

// --- Doxorubicin PK (2-cmt) ---
CLDox  : 45.0 : Doxorubicin clearance (L/h) [population mean]
VcDox  : 25.0 : Central volume (L)
QDox   : 15.0 : Intercompartmental CL (L/h)
VpDox  : 400.0: Peripheral volume (L) [extensive tissue binding]
Dox_dose: 0.0 : Doxorubicin dose (mg)

// --- Venetoclax PK (1-cmt, oral) ---
kaVEN  : 0.2  : Absorption rate (1/h)
FVEN   : 0.5  : Oral bioavailability (with food)
kelVEN : 0.035: Elimination rate (1/h) [T1/2~19h]
VdVEN  : 256.0: Volume of distribution (L)
VEN_dose: 0.0 : Venetoclax dose (mg/day, oral)

// --- Ibrutinib PK (1-cmt, oral) ---
kaIBR  : 0.5  : Absorption rate (1/h)
FIBR   : 0.03 : Oral bioavailability (extensive first-pass ~3%)
kelIBR : 0.22 : Elimination rate (1/h) [T1/2~3h]
VdIBR  : 10000.0: Volume (L) [widely distributed]
IBR_dose: 0.0 : Ibrutinib dose (mg/day, oral)

// --- Tumor Dynamics ---
kg_tumor: 0.012: Tumor growth rate (1/day) [doubling ~58 days]
kd_base : 0.003: Baseline tumor death rate (1/day)
K_carry : 1000.0: Tumor carrying capacity (arbitrary units)
T0      : 100.0 : Initial tumor burden (au)
EC50_RTX: 0.15  : RTX concentration for 50% max tumor kill (mg/L)
Emax_RTX: 0.85  : Max fractional tumor kill by RTX
EC50_CPP: 0.8   : 4-OH-CPP EC50 (mg/L)
Emax_CPP: 0.80  : Max fractional kill by cyclophosphamide
EC50_Dox: 0.05  : Doxorubicin EC50 (mg/L)
Emax_Dox: 0.75  : Max fractional kill by doxorubicin
EC50_VEN: 0.3   : Venetoclax EC50 (mg/L)
Emax_VEN: 0.80  : Max fractional kill by venetoclax
EC50_IBR: 0.12  : Ibrutinib EC50 (nM, converted) effective ABC subtype
Emax_IBR: 0.50  : Max fractional kill by ibrutinib (partial in ABC)
Hill_n  : 1.5   : Hill coefficient for drug effects

// --- BCR/NF-kB Signaling (simplified PD transduction) ---
kBCR_base: 0.02 : Baseline BCR signal activation rate (1/h)
kBCR_deg : 0.15 : BCR signal decay rate (1/h)
BCR_max  : 1.0  : Maximum signal (normalized)
IBR_BCR_IC50: 0.08: Ibrutinib IC50 for BTK inhibition (mg/L equivalent)

// --- BCL-2 Occupancy (Venetoclax PD) ---
BCL2_Kd  : 0.2  : Venetoclax BCL-2 binding affinity (mg/L)
BCL2_0   : 1.0  : Baseline BCL-2 expression (normalized)

// --- NK Cell Dynamics (ADCC effectors) ---
kNK_in  : 0.5  : NK cell influx rate (normalized, 1/day)
kNK_out : 0.1  : NK cell turnover (1/day)
NK0     : 5.0  : Baseline NK cells (normalized units)
RTX_NK_stim: 0.3: RTX ADCC amplification of NK activity

// --- CD8+ T Cell Dynamics ---
kCD8_in : 0.3  : CD8+ T cell influx (1/day)
kCD8_out: 0.05 : CD8+ T cell turnover (1/day)
CD8_0   : 10.0 : Baseline CD8+ cells (normalized)
kexhaust: 0.01 : Exhaustion rate (1/day, PD-L1 mediated)

// --- ANC (Absolute Neutrophil Count) Dynamics ---
ANC_0   : 5.0  : Baseline ANC (×10⁹/L)
kANC_rec: 0.08 : ANC recovery rate (1/day)
kANC_CPP_kill: 0.25: Cyclophosphamide neutropenia rate constant
kANC_Dox_kill: 0.15: Doxorubicin neutropenia rate constant
ANC_min : 0.1  : Minimum ANC floor (×10⁹/L)

// --- Resistance ---
kresist : 0.0003: Rate of resistance development (1/day × tumor)
resist_0: 0.01  : Initial resistance fraction

// --- Disease subtype flags (0=GCB, 1=ABC) ---
is_ABC  : 0.0   : Subtype: 0=GCB-DLBCL, 1=ABC-DLBCL
is_DHL  : 0.0   : Double-Hit flag (MYC+BCL2; worse prognosis)
DHL_growth_mult: 1.5: DHL growth rate multiplier

$CMT @annotated
// PK compartments
RuxCent   : Rituximab central (mg)
RuxPeriph : Rituximab peripheral (mg)
CD20_free : Free CD20 antigen on tumor surface (nM·cell equiv)
CD20_RTX  : CD20-RTX complex (nM·cell equiv)
CPP_active: Cyclophosphamide active metabolite (mg)
DoxCent   : Doxorubicin central (mg)
DoxPeriph : Doxorubicin peripheral (mg)
VEN_gut   : Venetoclax gut absorption compartment (mg)
VEN_cent  : Venetoclax central (mg)
IBR_cent  : Ibrutinib central (mg)

// PD compartments
Tumor     : Tumor burden (arbitrary units, normalized to T0=100)
BCR_signal: BCR/NF-kB activation (normalized 0-1)
BCL2_occ  : BCL-2 occupancy by Venetoclax (0-1)
NK_cells  : NK cell pool (normalized)
CD8_cells : CD8+ T-cell pool (normalized)
ANC       : Absolute neutrophil count (×10⁹/L)
Resistance: Resistance fraction (0-1)
Cum_RTX_dose: Cumulative rituximab exposure (mg·day/L)
Tumor_resp  : Tumor response indicator (1=initial, decreasing=response)

// Additional toxicity marker
CRS_risk  : CRS risk index (0-1, for novel agents)

$MAIN
// --- Rituximab dose as bolus (mg converted for TMDD in nM) ---
// Tumor_growth_rate accounts for DHL aggressiveness
double kg_eff = kg_tumor * (1 + is_DHL * (DHL_growth_mult - 1));

// Rituximab concentration (mg/L)
double CRux = RuxCent / VcR;
// 4-OH-cyclophosphamide conc (mg/L)
double CCPP = CPP_active / VdCPP;
// Doxorubicin central conc (mg/L)
double CDox = DoxCent / VcDox;
// Venetoclax conc (mg/L)
double CVEN = VEN_cent / VdVEN;
// Ibrutinib conc (mg/L) → convert to nM for IC50
double CIBR = IBR_cent / VdIBR;

$ODE
// ====================================================================
// RITUXIMAB PK (2-cmt + TMDD)
// ====================================================================
double koffR_eff = koffR + kintR;  // effective off rate
dxdt_RuxCent   = -CLR/VcR * RuxCent
                 - QR/VcR * RuxCent + QR/VpR * RuxPeriph
                 - konR * CRux * CD20_free / VcR
                 + koffR * CD20_RTX / VcR;
dxdt_RuxPeriph = QR/VcR * RuxCent - QR/VpR * RuxPeriph;

// CD20 dynamics: synthesis on tumor, degradation, RTX binding
double CD20_prod = ksynCD20 * (Tumor / 100.0);  // proportional to tumor cells
dxdt_CD20_free = CD20_prod - kdegCD20 * CD20_free
                 - konR * CRux * CD20_free + koffR * CD20_RTX;
dxdt_CD20_RTX  = konR * CRux * CD20_free - koffR_eff * CD20_RTX;

// ====================================================================
// CYCLOPHOSPHAMIDE PK (1-cmt active metabolite)
// ====================================================================
dxdt_CPP_active = -kelCPP * CPP_active;
// Dose added via $TABLE via event system (IV bolus at t=0 on Day 1 of each cycle)

// ====================================================================
// DOXORUBICIN PK (2-cmt)
// ====================================================================
dxdt_DoxCent   = -CLDox/VcDox * DoxCent - QDox/VcDox * DoxCent + QDox/VpDox * DoxPeriph;
dxdt_DoxPeriph = QDox/VcDox * DoxCent - QDox/VpDox * DoxPeriph;

// ====================================================================
// VENETOCLAX PK (1-cmt oral)
// ====================================================================
dxdt_VEN_gut  = -kaVEN * VEN_gut;
dxdt_VEN_cent = FVEN * kaVEN * VEN_gut - kelVEN * VEN_cent;

// ====================================================================
// IBRUTINIB PK (1-cmt oral)
// ====================================================================
dxdt_IBR_cent = FIBR * (-kelIBR * IBR_cent);  // simplified 1-cmt with bioav

// ====================================================================
// BCR/NF-kB SIGNAL (simplified transduction chain)
// ====================================================================
double BCR_drive = kBCR_base;   // tonic BCR signal
// Ibrutinib inhibits BTK -> reduces BCR signal in ABC
double IBR_BCR_inh = (is_ABC > 0.5) ?
    (CIBR / (CIBR + IBR_BCR_IC50)) * 0.8 : 0.0;
dxdt_BCR_signal = BCR_drive * (BCR_max - BCR_signal) * (1 - IBR_BCR_inh)
                  - kBCR_deg * BCR_signal;

// ====================================================================
// BCL-2 OCCUPANCY (Venetoclax PD)
// ====================================================================
double VEN_occ = CVEN / (CVEN + BCL2_Kd);
dxdt_BCL2_occ = kaVEN * (VEN_occ - BCL2_occ);  // first-order approach to equil

// ====================================================================
// NK CELL DYNAMICS
// ====================================================================
// RTX enhances NK via ADCC signaling
double NK_RTX_boost = 1.0 + RTX_NK_stim * (CD20_RTX / (CD20_RTX + 0.05));
dxdt_NK_cells = kNK_in * NK_RTX_boost - kNK_out * NK_cells;

// ====================================================================
// CD8+ T CELL DYNAMICS
// ====================================================================
// Exhaustion increases with tumor burden and PD-L1 (simplified)
double PD1_exhaust = kexhaust * (Tumor / 100.0);
dxdt_CD8_cells = kCD8_in - (kCD8_out + PD1_exhaust) * CD8_cells;

// ====================================================================
// TUMOR DYNAMICS
// ====================================================================
// Drug kill effects (additive, capped at 1)
double E_RTX  = Emax_RTX * pow(CRux, Hill_n) /
                (pow(CRux, Hill_n) + pow(EC50_RTX, Hill_n));
// ADCC component from NK cells
double E_ADCC = E_RTX * (NK_cells / NK0) * 1.2;
E_ADCC = (E_ADCC > Emax_RTX) ? Emax_RTX : E_ADCC;

double E_CPP  = Emax_CPP * pow(CCPP, Hill_n) /
                (pow(CCPP, Hill_n) + pow(EC50_CPP, Hill_n));
double E_Dox  = Emax_Dox * pow(CDox, Hill_n) /
                (pow(CDox, Hill_n) + pow(EC50_Dox, Hill_n));
double E_VEN  = Emax_VEN * BCL2_occ;  // proportional to BCL-2 occupancy
// Ibrutinib effect in ABC-DLBCL via BCR signal suppression
double E_IBR  = (is_ABC > 0.5) ?
    Emax_IBR * CIBR / (CIBR + EC50_IBR) * BCR_signal : 0.0;

// NK and CD8 killing contribution
double E_immune = 0.01 * (NK_cells / NK0 + CD8_cells / CD8_0) * (Tumor / 100.0);

// Combined kill (Bliss independence approximation)
double E_total = 1 - (1 - E_ADCC) * (1 - E_CPP) * (1 - E_Dox) *
                     (1 - E_VEN) * (1 - E_IBR);
if(E_total > 0.98) E_total = 0.98;

// Resistance effect: reduces drug efficacy
double resist_factor = 1 - Resistance * 0.8;
E_total = E_total * resist_factor;

// Logistic tumor growth with drug kill
double growth_term = kg_eff * Tumor * (1 - Tumor / K_carry);
double kill_term = (kd_base + E_total) * Tumor + E_immune;
dxdt_Tumor = growth_term - kill_term;
if(Tumor < 0.01) dxdt_Tumor = 0;

// ====================================================================
// RESISTANCE EVOLUTION
// ====================================================================
dxdt_Resistance = kresist * Tumor * (1 - Resistance);

// ====================================================================
// ANC DYNAMICS (myelosuppression)
// ====================================================================
double ANC_drive = kANC_rec * (ANC_0 - ANC);
double ANC_CPP_kill = kANC_CPP_kill * CCPP * ANC / ANC_0;
double ANC_Dox_kill = kANC_Dox_kill * CDox * ANC / ANC_0;
dxdt_ANC = ANC_drive - ANC_CPP_kill - ANC_Dox_kill;
if(ANC < ANC_min) dxdt_ANC = 0;

// ====================================================================
// CUMULATIVE RTX EXPOSURE
// ====================================================================
dxdt_Cum_RTX_dose = CRux;

// ====================================================================
// TUMOR RESPONSE TRACKER
// ====================================================================
dxdt_Tumor_resp = 0;  // updated in $TABLE

// ====================================================================
// CRS RISK INDEX (novel agents / CAR-T surrogate)
// ====================================================================
// Rapid tumor kill rate drives CRS
double tumor_kill_rate = kill_term;
dxdt_CRS_risk = 0.1 * tumor_kill_rate * (1 - CRS_risk) - 0.05 * CRS_risk;

$TABLE
// Derived concentrations for output
capture CRux_mgL  = RuxCent / VcR;          // Rituximab central conc (mg/L)
capture CCPP_mgL  = CPP_active / VdCPP;     // Active CPP metabolite (mg/L)
capture CDox_mgL  = DoxCent / VcDox;        // Doxorubicin central (mg/L)
capture CVEN_mgL  = VEN_cent / VdVEN;       // Venetoclax (mg/L)
capture CIBR_mgL  = IBR_cent / VdIBR;       // Ibrutinib (mg/L)

// Tumor response classification (Lugano-like)
// CR: <5% of T0; PR: 5-50% of T0; SD: 50-150%; PD: >150%
double resp_frac = Tumor / T0;
capture Resp_frac = resp_frac;              // Fraction of initial tumor
capture CR_flag   = (resp_frac < 0.05) ? 1 : 0;
capture PR_flag   = (resp_frac >= 0.05 && resp_frac < 0.50) ? 1 : 0;
capture PD_flag   = (resp_frac > 1.50) ? 1 : 0;

// ANC categories
capture ANC_grade  = (ANC < 0.5) ? 4 : (ANC < 1.0) ? 3 :
                     (ANC < 1.5) ? 2 : (ANC < 2.0) ? 1 : 0;
capture febrile_risk = (ANC < 0.5) ? 1 : 0;

// BCR/NF-kB signal output
capture BCR_sig_out = BCR_signal;

// Relative SPD (Sum of Product Diameters, proportional to tumor burden)
capture SPD_relative = Tumor / T0 * 100;    // % of baseline SPD

$INIT
RuxCent    = 0
RuxPeriph  = 0
CD20_free  = 0.30   // near Bmax at baseline
CD20_RTX   = 0
CPP_active = 0
DoxCent    = 0
DoxPeriph  = 0
VEN_gut    = 0
VEN_cent   = 0
IBR_cent   = 0
Tumor      = 100    // T0 = 100 normalized units
BCR_signal = 0.133  // baseline tonic signal (= kBCR_base/kBCR_deg)
BCL2_occ   = 0
NK_cells   = 5      // NK0
CD8_cells  = 10     // CD8_0
ANC        = 5.0    // normal ANC (×10⁹/L)
Resistance = 0.01   // 1% initial resistance
Cum_RTX_dose = 0
Tumor_resp = 100
CRS_risk   = 0

$OMEGA @labels ECL EVCD20 ETumor EANCR
0.04 0.09 0.16 0.04

$SIGMA
0.02

$CAPTURE CRux_mgL CCPP_mgL CDox_mgL CVEN_mgL CIBR_mgL
         Resp_frac CR_flag PR_flag PD_flag
         ANC_grade febrile_risk BCR_sig_out SPD_relative
         Resistance NK_cells CD8_cells CRS_risk
'

# Compile the model
dlbcl_mod <- mcode("dlbcl_qsp", code)

##############################################################################
# DOSING EVENTS — Helper function for R-CHOP cycles
##############################################################################

make_rchop_events <- function(n_cycles = 6,
                               cycle_len = 21,       # days
                               start_day = 0,
                               rtx_dose_mg  = 1350,  # ~375 mg/m² × 3.6 m² BSA
                               cpp_dose_mg  = 2700,  # ~750 mg/m² × 3.6 m²
                               dox_dose_mg  = 180,   # ~50 mg/m² × 3.6 m²
                               vin_note     = TRUE,  # Vincristine noted only
                               pred_note    = TRUE) {
  ev_list <- vector("list", n_cycles * 3)
  idx <- 1
  for (i in seq_len(n_cycles)) {
    day <- start_day + (i - 1) * cycle_len
    # RTX (Day 1 of each cycle) — 2-hour IV bolus → model as rapid bolus
    ev_list[[idx]] <- ev(cmt = 1, time = day * 24,
                         amt  = rtx_dose_mg, rate = -2)  # -2 = infusion 2h
    idx <- idx + 1
    # Cyclophosphamide IV bolus (Day 1)
    ev_list[[idx]] <- ev(cmt = 5, time = day * 24,
                         amt  = cpp_dose_mg, rate = -2)
    idx <- idx + 1
    # Doxorubicin IV bolus (Day 1)
    ev_list[[idx]] <- ev(cmt = 6, time = day * 24,
                         amt  = dox_dose_mg, rate = -1)
    idx <- idx + 1
  }
  do.call(c, ev_list)
}

# Venetoclax daily oral (continuous dosing)
make_venetoclax_events <- function(dose_mg = 800, days = 180, interval = 24) {
  ev(cmt = 8, time = 0, amt = dose_mg, ii = interval, addl = days - 1)
}

# Ibrutinib daily oral
make_ibrutinib_events <- function(dose_mg = 560, days = 365, interval = 24) {
  ev(cmt = 10, time = 0, amt = dose_mg, ii = interval, addl = days - 1)
}

##############################################################################
# SCENARIO DEFINITIONS
##############################################################################

scenarios <- list(
  list(
    name    = "1. Untreated DLBCL (GCB)",
    params  = list(is_ABC = 0, is_DHL = 0),
    events  = ev(cmt = 1, time = 0, amt = 0),
    end_day = 365, label = "Untreated GCB"
  ),
  list(
    name    = "2. R-CHOP x6 (Standard, GCB)",
    params  = list(is_ABC = 0, is_DHL = 0),
    events  = make_rchop_events(n_cycles = 6),
    end_day = 365, label = "R-CHOP×6 GCB"
  ),
  list(
    name    = "3. Pola-R-CHP x6 (POLARIX, GCB)",
    params  = list(is_ABC = 0, is_DHL = 0,
                   Emax_RTX = 0.90, EC50_CPP = 0.7),  # Pola adds MMAE effect (proxy)
    events  = make_rchop_events(n_cycles = 6, rtx_dose_mg = 1350,
                                 cpp_dose_mg = 2700, dox_dose_mg = 180),
    end_day = 365, label = "Pola-R-CHP×6"
  ),
  list(
    name    = "4. R-CHOP + Ibrutinib (ABC-DLBCL, PHOENIX)",
    params  = list(is_ABC = 1, is_DHL = 0),
    events  = c(make_rchop_events(n_cycles = 6),
                make_ibrutinib_events(dose_mg = 560, days = 126)),
    end_day = 365, label = "R-CHOP+Ibrutinib ABC"
  ),
  list(
    name    = "5. Venetoclax + R-CHOP (BCL-2 high, CAVALLI)",
    params  = list(is_ABC = 0, is_DHL = 0,
                   EC50_VEN = 0.25, Emax_VEN = 0.85),
    events  = c(make_rchop_events(n_cycles = 6),
                make_venetoclax_events(dose_mg = 800, days = 126)),
    end_day = 365, label = "Venetoclax+R-CHOP"
  ),
  list(
    name    = "6. R-CHOP x6 (Double-Hit Lymphoma, DHL)",
    params  = list(is_ABC = 0, is_DHL = 1),
    events  = make_rchop_events(n_cycles = 6),
    end_day = 365, label = "R-CHOP DHL (high-risk)"
  )
)

##############################################################################
# SIMULATION
##############################################################################

run_scenario <- function(scenario, mod, delta = 6) {
  params_list <- scenario$params
  param_vec   <- as.numeric(unlist(params_list))
  param_names <- names(params_list)
  p <- do.call(param, c(list(mod), params_list))
  out <- mrgsim(p, events = scenario$events,
                end = scenario$end_day * 24, delta = delta,
                carry.out = "evid")
  df <- as.data.frame(out)
  df$Scenario <- scenario$label
  df$time_day <- df$time / 24
  df
}

# Run all 6 scenarios
cat("Running 6 DLBCL treatment scenarios...\n")
set.seed(42)
results <- lapply(scenarios, run_scenario, mod = dlbcl_mod)
all_results <- bind_rows(results)

##############################################################################
# POPULATION SIMULATION (N=100 virtual patients, Scenario 2: R-CHOP)
##############################################################################

cat("Population simulation (N=100 virtual patients)...\n")
idata_pop <- data.frame(
  ID   = 1:100,
  kg_tumor = rlnorm(100, log(0.012), 0.4),   # between-patient variability
  T0       = rlnorm(100, log(100), 0.5),
  EC50_RTX = rlnorm(100, log(0.15), 0.3)
)

pop_sim <- mrgsim(dlbcl_mod,
                  idata = idata_pop,
                  events = make_rchop_events(n_cycles = 6),
                  end    = 365 * 24,
                  delta  = 24,
                  carry.out = "evid")
pop_df <- as.data.frame(pop_sim)
pop_df$time_day <- pop_df$time / 24

##############################################################################
# SUMMARY STATISTICS
##############################################################################

summary_tab <- all_results %>%
  filter(time_day %in% c(0, 42, 84, 126, 168, 252, 365)) %>%
  group_by(Scenario, time_day) %>%
  summarise(
    Tumor_median = round(median(Tumor), 1),
    SPD_pct      = round(median(SPD_relative), 1),
    CR_rate      = round(mean(CR_flag), 3),
    ANC_mean     = round(mean(ANC), 2),
    CRS_risk_max = round(max(CRS_risk), 3),
    .groups = "drop"
  )

cat("\n=== Tumor Response Summary at Key Timepoints ===\n")
print(summary_tab, n = Inf)

# CR and PD rates at end of treatment (~Day 126 = 18 weeks)
eot_resp <- all_results %>%
  filter(abs(time_day - 126) < 3) %>%
  group_by(Scenario) %>%
  summarise(
    CR_rate = round(mean(CR_flag) * 100, 1),
    PR_rate = round(mean(PR_flag) * 100, 1),
    PD_rate = round(mean(PD_flag) * 100, 1),
    Median_SPD_pct = round(median(SPD_relative), 1)
  )
cat("\n=== End-of-Treatment Response Rates (Day 126) ===\n")
print(eot_resp)

# Population-level CR rate at Day 126
pop_eot <- pop_df %>%
  filter(abs(time_day - 126) < 1) %>%
  summarise(
    CR_rate_pop = round(mean(CR_flag) * 100, 1),
    median_SPD  = round(median(SPD_relative), 1),
    p25_SPD     = round(quantile(SPD_relative, 0.25), 1),
    p75_SPD     = round(quantile(SPD_relative, 0.75), 1)
  )
cat("\n=== Population R-CHOP: Day-126 CR Rate (N=100) ===\n")
print(pop_eot)

##############################################################################
# VISUALIZATION
##############################################################################

theme_dlbcl <- theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#E3F2FD"),
        panel.grid.minor = element_blank())

# Figure 1: Tumor Burden Over Time — All Scenarios
p1 <- ggplot(all_results, aes(x = time_day, y = SPD_relative,
                               color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = c(5, 50), linetype = "dashed",
             color = c("#4CAF50", "#FF9800"), linewidth = 0.8) +
  annotate("text", x = 360, y = 7,  label = "CR threshold (5%)",  size = 3,
           color = "#4CAF50") +
  annotate("text", x = 360, y = 52, label = "PR threshold (50%)", size = 3,
           color = "#FF9800") +
  scale_y_continuous(limits = c(0, 250)) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "DLBCL Tumor Burden Dynamics — 6 Treatment Scenarios",
       subtitle = "SPD = Sum of Product Diameters (% of baseline)",
       x = "Time (days)", y = "SPD (% of baseline)",
       color = "Treatment") +
  theme_dlbcl
print(p1)

# Figure 2: Rituximab PK — Concentration-Time Profile
rtx_pk <- all_results %>% filter(grepl("R-CHOP", Scenario) | Scenario == "Venetoclax+R-CHOP")
p2 <- ggplot(rtx_pk, aes(x = time_day, y = CRux_mgL, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_y_log10() +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Rituximab PK — Central Compartment Concentration",
       subtitle = "2-compartment model with TMDD (CD20 binding)",
       x = "Time (days)", y = "Rituximab [mg/L] (log scale)") +
  theme_dlbcl
print(p2)

# Figure 3: ANC Over Time (Myelosuppression)
p3 <- ggplot(all_results %>% filter(!grepl("Untreated", Scenario)),
             aes(x = time_day, y = ANC, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = c(0.5, 1.0), linetype = "dashed",
             color = c("red", "orange"), linewidth = 0.7) +
  annotate("text", x = 360, y = 0.6, label = "Grade 4 threshold (0.5)",
           size = 3, color = "red") +
  annotate("text", x = 360, y = 1.1, label = "Grade 3 threshold (1.0)",
           size = 3, color = "orange") +
  scale_color_brewer(palette = "Set2") +
  labs(title = "ANC Dynamics — Myelosuppression Profile",
       subtitle = "Cyclophosphamide + Doxorubicin-induced neutropenia",
       x = "Time (days)", y = "ANC (×10⁹/L)") +
  theme_dlbcl
print(p3)

# Figure 4: NK and CD8+ T cell dynamics
imm_df <- all_results %>%
  select(time_day, Scenario, NK_cells, CD8_cells) %>%
  pivot_longer(cols = c(NK_cells, CD8_cells),
               names_to = "cell_type", values_to = "count")
p4 <- ggplot(imm_df, aes(x = time_day, y = count, color = Scenario,
                          linetype = cell_type)) +
  geom_line(linewidth = 0.9) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Immune Effector Cell Dynamics",
       subtitle = "NK cell ADCC enhancement by Rituximab; CD8+ T cell exhaustion",
       x = "Time (days)", y = "Cell count (normalized)",
       linetype = "Cell type") +
  theme_dlbcl
print(p4)

# Figure 5: Population Variability (R-CHOP N=100)
pop_summary <- pop_df %>%
  group_by(time_day) %>%
  summarise(
    p50 = median(SPD_relative),
    p25 = quantile(SPD_relative, 0.25),
    p75 = quantile(SPD_relative, 0.75),
    p05 = quantile(SPD_relative, 0.05),
    p95 = quantile(SPD_relative, 0.95)
  )
p5 <- ggplot(pop_summary, aes(x = time_day)) +
  geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.15, fill = "#1565C0") +
  geom_ribbon(aes(ymin = p25, ymax = p75), alpha = 0.30, fill = "#1565C0") +
  geom_line(aes(y = p50), color = "#0D47A1", linewidth = 1.5) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "#4CAF50") +
  labs(title = "Population Variability — R-CHOP (N = 100 Virtual Patients)",
       subtitle = "Median (solid) ± IQR (dark band) ± 90th percentile (light band)",
       x = "Time (days)", y = "SPD (% of baseline)") +
  theme_dlbcl
print(p5)

# Figure 6: Response Waterfall at Day 126
waterfall_df <- pop_df %>%
  filter(abs(time_day - 126) < 1) %>%
  arrange(SPD_relative) %>%
  mutate(rank = row_number(),
         response = case_when(
           SPD_relative < 5   ~ "CR",
           SPD_relative < 50  ~ "PR",
           SPD_relative > 150 ~ "PD",
           TRUE               ~ "SD"
         ))
p6 <- ggplot(waterfall_df, aes(x = rank, y = (SPD_relative - 100),
                                 fill = response)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = c(-95, -50, 50), linetype = "dashed", linewidth = 0.6) +
  scale_fill_manual(values = c("CR" = "#4CAF50", "PR" = "#8BC34A",
                                 "SD" = "#FFC107", "PD" = "#F44336")) +
  labs(title = "Waterfall Plot — Tumor Response at Day 126 (N=100)",
       subtitle = "R-CHOP regimen; % change from baseline SPD",
       x = "Patient rank", y = "SPD change from baseline (%)",
       fill = "Response") +
  theme_dlbcl
print(p6)

cat("\n=== DLBCL QSP Model Simulation Complete ===\n")
cat("Calibration targets:\n")
cat("  R-CHOP CR rate ~65-70% → Model Day-126 CR:", pop_eot$CR_rate_pop, "%\n")
cat("  Clinical ORR (CR+PR) ~75-85% for R-CHOP × 6 cycles\n")
cat("  CAVALLI venetoclax+R-CHOP: BCL2+ ORR 88% (Morschhauser 2021)\n")
cat("  POLARIX Pola-R-CHP: 2yr PFS 76.7% vs R-CHOP 70.2% (Tilly 2022 NEJM)\n")
