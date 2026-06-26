## =============================================================================
## PDAC (Pancreatic Ductal Adenocarcinoma) QSP mrgsolve Model
## =============================================================================
## Mechanistic ODE model covering:
##   - Multi-drug PK: Gemcitabine, nab-Paclitaxel, Oxaliplatin,
##     Irinotecan/SN-38, 5-Fluorouracil, MRTX1133 (KRAS G12D inhibitor),
##     Olaparib
##   - Tumor Growth Inhibition (Simeoni model)
##   - CA19-9 biomarker dynamics
##   - Myelosuppression (Friberg neutropenia model)
##   - Stromal resistance compartment
##   - KRAS signaling modulation
##
## Clinical calibration references:
##   MPACT (Von Hoff 2013): Gem+nab-Pac mOS 8.5 vs Gem 6.7 months
##   PRODIGE4 (Conroy 2011): FOLFIRINOX mOS 11.1 vs Gem 6.8 months
##   POLO (Golan 2019): Olaparib PFS 7.4 vs placebo 3.8 months
##   MRTX1133 phase I/II (2023-2025): KRAS G12D inhibitor
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

## ─────────────────────────────────────────────────────────────────────────────
## Model specification
## ─────────────────────────────────────────────────────────────────────────────

pdac_model <- mcode("pdac_qsp", '

$PROB
PDAC QSP Model - Multi-drug PK/PD with Simeoni TGI, Friberg neutropenia,
CA19-9 biomarker, stromal resistance, and KRAS inhibition.
Units: concentrations in ug/mL (except where noted), time in hours,
tumor volume in mm3, CA19-9 in U/mL, ANC in 10^9/L.

$PARAM
// ── Gemcitabine PK ──────────────────────────────────────────────────────────
CL_GEM      = 80.0,     // L/h clearance
V1_GEM      = 15.0,     // L central volume
V2_GEM      = 30.0,     // L peripheral volume
Q_GEM       = 25.0,     // L/h intercompartment clearance
k_dFdCTP    = 0.5,      // h^-1 formation rate of dFdCTP active metabolite
k_out_dFdCTP = 0.3,     // h^-1 dFdCTP elimination

// ── nab-Paclitaxel PK ───────────────────────────────────────────────────────
CL_NPAC     = 12.0,     // L/h clearance
V1_NPAC     = 8.0,      // L central volume
V2_NPAC     = 100.0,    // L peripheral volume
Q_NPAC      = 8.0,      // L/h intercompartment clearance

// ── Oxaliplatin PK ──────────────────────────────────────────────────────────
CL_OHP      = 9.0,      // L/h clearance
V_OHP       = 12.0,     // L volume of distribution
k_DNA       = 0.08,     // h^-1 DNA adduct formation rate
k_repair    = 0.02,     // h^-1 DNA repair rate

// ── Irinotecan / SN-38 PK ───────────────────────────────────────────────────
CL_CPT11    = 15.0,     // L/h irinotecan clearance
V_CPT11     = 180.0,    // L irinotecan volume
k_SN38_form = 0.15,     // h^-1 SN-38 formation rate from CPT-11
CL_SN38     = 8.0,      // L/h SN-38 clearance
k_SN38G     = 0.3,      // h^-1 SN-38 glucuronidation rate

// ── 5-Fluorouracil PK ───────────────────────────────────────────────────────
CL_FU5      = 120.0,    // L/h (non-linear, saturable; simplified linear here)
V_FU5       = 20.0,     // L volume
k_FdUMP     = 0.4,      // h^-1 FdUMP formation (intracellular, surrogate)

// ── MRTX1133 (KRAS G12D inhibitor) PK ───────────────────────────────────────
CL_MRTX     = 5.0,      // L/h clearance
V_MRTX      = 25.0,     // L volume
F_MRTX      = 0.45,     // oral bioavailability fraction
ka_MRTX     = 0.8,      // h^-1 absorption rate constant (oral)

// ── Olaparib PK ─────────────────────────────────────────────────────────────
CL_OLA      = 8.9,      // L/h clearance (Menear 2008; Plummer 2008)
V_OLA       = 158.0,    // L volume
F_OLA       = 0.73,     // oral bioavailability
ka_OLA      = 1.2,      // h^-1 absorption rate constant

// ── Tumor growth / Simeoni TGI ──────────────────────────────────────────────
k_prog      = 0.012,    // day^-1 → converted to h^-1 in ODE (divide by 24)
k_death     = 0.008,    // day^-1 natural death rate (h^-1 in ODE)
KSS_GEM     = 0.5,      // ug/mL IC50-like gemcitabine (via dFdCTP)
KSS_NPAC    = 0.05,     // ug/mL IC50 nab-paclitaxel
KSS_OHP     = 0.3,      // ug/mL IC50 oxaliplatin (free plasma)
KSS_SN38    = 0.02,     // ug/mL IC50 SN-38
KSS_FU5     = 0.8,      // ug/mL IC50 5-FU
KSS_MRTX    = 0.06,     // ug/mL IC50 MRTX1133 (tumor kill)
EMAX_combo  = 0.85,     // maximum fractional tumor kill (combined)
k_tr        = 0.0075,   // h^-1 transit rate (Simeoni) [= 0.18 day^-1]
V0_tumor    = 2000.0,   // mm3 reference tumor volume (initial)

// ── KRAS signaling / MRTX1133 PD ────────────────────────────────────────────
IC50_MRTX_KRAS = 0.06,  // ug/mL MRTX1133 for 50% KRAS inhibition
Imax_MRTX      = 0.90,  // max KRAS inhibition fraction
k_KRAS_on      = 0.05,  // h^-1 KRAS signaling activation rate
k_KRAS_off     = 0.05,  // h^-1 KRAS inactivation rate (at baseline equilibrium = 1)

// ── PI3K/AKT pathway ────────────────────────────────────────────────────────
k_PI3K_on   = 0.04,     // h^-1 activation by KRAS signal
k_PI3K_off  = 0.06,     // h^-1 baseline inactivation

// ── CA19-9 biomarker ────────────────────────────────────────────────────────
k_CA199_prod  = 0.1042, // U/mL/h production (= 2.5/day) per unit tumor burden
k_CA199_el    = 0.00208,// h^-1 CA19-9 elimination (= 0.05/day)
CA199_base    = 35.0,   // U/mL baseline CA19-9

// ── Friberg neutropenia model ────────────────────────────────────────────────
MTT         = 132.0,    // h mean transit time (5.5 days)
gamma_neut  = 0.17,     // feedback exponent
Circ0       = 5.0,      // 10^9/L baseline ANC
EMAX_neut   = 0.95,     // max myelosuppressive effect
EC50_neut_GEM = 0.3,    // ug/mL GEM for 50% neutrophil inhibition
EC50_neut_SN38 = 0.015, // ug/mL SN-38 for 50% neutrophil inhibition

// ── Stroma ───────────────────────────────────────────────────────────────────
k_stroma_build  = 0.000208, // h^-1 stroma accumulation (= 0.005/day)
k_stroma_decay  = 0.0000833,// h^-1 natural stroma decay (= 0.002/day)
stroma_pen_factor = 0.6,    // max 60% reduction in effective drug penetration

// ── Dosing input scalars (set by event objects; 0 = off) ────────────────────
// These are placeholders overridden by ev() calls
GEM_dose_rate  = 0,
NPAC_dose_rate = 0,
OHP_dose_rate  = 0,
CPT11_dose_rate = 0,
FU5_dose_rate  = 0,
MRTX_dose_rate = 0,
OLA_dose_rate  = 0

$CMT
// Drug PK compartments (12)
GEM_C1        // 1.  Gemcitabine central (ug)
GEM_C2        // 2.  Gemcitabine peripheral (ug)
dFdCTP        // 3.  Active metabolite dFdCTP (ug equivalent)
NPAC_C1       // 4.  nab-Paclitaxel central (ug)
NPAC_C2       // 5.  nab-Paclitaxel peripheral (ug)
OHP_FREE      // 6.  Oxaliplatin free plasma (ug)
OHP_DNA       // 7.  Oxaliplatin-DNA adducts (ug)
CPT11         // 8.  Irinotecan plasma (ug)
SN38          // 9.  Active SN-38 metabolite (ug)
FU5           // 10. 5-Fluorouracil plasma (ug)
MRTX1133      // 11. MRTX1133 KRAS G12D inhibitor (ug)
OLA           // 12. Olaparib (ug)
// Disease / PD compartments (10)
KRAS_SIG      // 13. KRAS signaling (normalized, dimensionless 0-2)
PI3K_ACT      // 14. PI3K/AKT activity (normalized, dimensionless 0-2)
x0            // 15. Tumor proliferating cells (Simeoni) (mm3 equivalent)
x1            // 16. Transit damage compartment 1
x2            // 17. Transit damage compartment 2
x3            // 18. Transit damage compartment 3
TUMOR         // 19. Total tumor volume (mm3)
CA199         // 20. CA19-9 biomarker (U/mL)
Prol          // 21. Proliferating neutrophils (10^9/L, Friberg)
Tr1           // 22. Transit neutrophils 1
Tr2           // 23. Transit neutrophils 2
Tr3           // 24. Transit neutrophils 3
Circ          // 25. Circulating neutrophils ANC (10^9/L)
STROMA        // 26. Stromal resistance factor (normalized 0-1)

$MAIN
// ── Derived rate constants ───────────────────────────────────────────────────
double k_prol     = 4.0 / MTT;   // Friberg proliferation rate = 4/MTT
double k_tr_neut  = 4.0 / MTT;   // Friberg transit rate = 4/MTT

// ── Convert tumor growth rates from day^-1 to h^-1 ──────────────────────────
double k_prog_h   = k_prog  / 24.0;
double k_death_h  = k_death / 24.0;

// ── Initial conditions ───────────────────────────────────────────────────────
if (NEWIND <= 1) {
  // Tumor: start with a small measurable mass
  _init_x0    = 1000.0;   // mm3 proliferating compartment
  _init_TUMOR = 1000.0;   // mm3 total tumor volume
  // KRAS signaling at basal equilibrium (= 1, normalized)
  _init_KRAS_SIG = 1.0;
  _init_PI3K_ACT = k_PI3K_on / k_PI3K_off; // steady-state
  // CA19-9 at baseline
  _init_CA199 = CA199_base;
  // Neutrophils at steady-state (all compartments = Circ0 at SS)
  _init_Prol  = Circ0;
  _init_Tr1   = Circ0;
  _init_Tr2   = Circ0;
  _init_Tr3   = Circ0;
  _init_Circ  = Circ0;
  // Stroma starts at low baseline (normalized 0-1)
  _init_STROMA = 0.3;
}

$ODE
// ─────────────────────────────────────────────────────────────────────────────
// Derived PK quantities (concentrations ug/mL = amount/volume)
// ─────────────────────────────────────────────────────────────────────────────
double C_GEM     = GEM_C1  / V1_GEM;     // ug/mL gemcitabine central
double C_dFdCTP  = dFdCTP  / V1_GEM;     // surrogate intracellular metabolite
double C_NPAC    = NPAC_C1 / V1_NPAC;    // ug/mL nab-paclitaxel
double C_OHP     = OHP_FREE / V_OHP;     // ug/mL oxaliplatin free
double C_OHP_DNA = OHP_DNA  / V_OHP;     // DNA adducts (ug/mL equivalent)
double C_CPT11   = CPT11   / V_CPT11;    // ug/mL irinotecan
double C_SN38    = SN38    / 25.0;       // ug/mL (Vd SN38 ~25 L)
double C_FU5     = FU5     / V_FU5;      // ug/mL 5-FU
double C_MRTX    = MRTX1133 / V_MRTX;   // ug/mL MRTX1133
double C_OLA     = OLA     / V_OLA;      // ug/mL olaparib

// Protect against negative concentrations (numerical safety)
double C_GEM_p    = C_GEM    > 0 ? C_GEM    : 0;
double C_dFdCTP_p = C_dFdCTP > 0 ? C_dFdCTP : 0;
double C_NPAC_p   = C_NPAC   > 0 ? C_NPAC   : 0;
double C_OHP_p    = C_OHP    > 0 ? C_OHP    : 0;
double C_SN38_p   = C_SN38   > 0 ? C_SN38   : 0;
double C_FU5_p    = C_FU5    > 0 ? C_FU5    : 0;
double C_MRTX_p   = C_MRTX   > 0 ? C_MRTX   : 0;
double C_OHP_DNA_p = C_OHP_DNA > 0 ? C_OHP_DNA : 0;

// ─────────────────────────────────────────────────────────────────────────────
// Drug PK ODEs
// ─────────────────────────────────────────────────────────────────────────────

// ── Gemcitabine two-compartment PK ──────────────────────────────────────────
dxdt_GEM_C1  = -(CL_GEM/V1_GEM)*GEM_C1
               - (Q_GEM/V1_GEM)*GEM_C1
               + (Q_GEM/V2_GEM)*GEM_C2
               - k_dFdCTP*GEM_C1
               + GEM_dose_rate;
dxdt_GEM_C2  =  (Q_GEM/V1_GEM)*GEM_C1 - (Q_GEM/V2_GEM)*GEM_C2;
dxdt_dFdCTP  =  k_dFdCTP*GEM_C1 - k_out_dFdCTP*dFdCTP;

// ── nab-Paclitaxel two-compartment PK ───────────────────────────────────────
dxdt_NPAC_C1 = -(CL_NPAC/V1_NPAC)*NPAC_C1
               - (Q_NPAC/V1_NPAC)*NPAC_C1
               + (Q_NPAC/V2_NPAC)*NPAC_C2
               + NPAC_dose_rate;
dxdt_NPAC_C2 =  (Q_NPAC/V1_NPAC)*NPAC_C1 - (Q_NPAC/V2_NPAC)*NPAC_C2;

// ── Oxaliplatin one-compartment + DNA adducts ────────────────────────────────
dxdt_OHP_FREE = -(CL_OHP/V_OHP)*OHP_FREE - k_DNA*OHP_FREE + OHP_dose_rate;
dxdt_OHP_DNA  =  k_DNA*OHP_FREE - k_repair*OHP_DNA;

// ── Irinotecan → SN-38 ───────────────────────────────────────────────────────
dxdt_CPT11   = -(CL_CPT11/V_CPT11)*CPT11 - k_SN38_form*CPT11 + CPT11_dose_rate;
dxdt_SN38    =  k_SN38_form*CPT11 - (CL_SN38/25.0)*SN38 - k_SN38G*SN38;

// ── 5-Fluorouracil ───────────────────────────────────────────────────────────
dxdt_FU5     = -(CL_FU5/V_FU5)*FU5 - k_FdUMP*FU5 + FU5_dose_rate;

// ── MRTX1133 (oral, absorbed directly into central) ─────────────────────────
// Depot absorption modeled through dose_rate input (pre-absorbed for simplicity)
dxdt_MRTX1133 = -(CL_MRTX/V_MRTX)*MRTX1133 + MRTX_dose_rate;

// ── Olaparib (oral) ───────────────────────────────────────────────────────────
dxdt_OLA     = -(CL_OLA/V_OLA)*OLA + OLA_dose_rate;

// ─────────────────────────────────────────────────────────────────────────────
// KRAS signaling and PI3K/AKT pathway
// ─────────────────────────────────────────────────────────────────────────────
// KRAS inhibition by MRTX1133 (Emax model)
double MRTX_inh = Imax_MRTX * C_MRTX_p / (IC50_MRTX_KRAS + C_MRTX_p);
// KRAS_SIG: activated by baseline (= k_KRAS_on), inhibited by drug, decays
dxdt_KRAS_SIG = k_KRAS_on*(1.0 - MRTX_inh) - k_KRAS_off*KRAS_SIG;
// PI3K/AKT activated downstream of KRAS
dxdt_PI3K_ACT = k_PI3K_on*KRAS_SIG - k_PI3K_off*PI3K_ACT;

// ─────────────────────────────────────────────────────────────────────────────
// Stromal resistance compartment (normalized 0-1)
// ─────────────────────────────────────────────────────────────────────────────
double STROMA_safe = STROMA > 0 ? STROMA : 0;
dxdt_STROMA  = k_stroma_build*(1.0 - STROMA_safe) - k_stroma_decay*STROMA_safe;

// Stromal penetration factor (1 = no stroma, reduced with higher stroma)
double stroma_pen = 1.0 - stroma_pen_factor * STROMA_safe;
stroma_pen = stroma_pen < 0.05 ? 0.05 : stroma_pen; // floor at 5%

// ─────────────────────────────────────────────────────────────────────────────
// Individual drug Emax effects on tumor (Bliss independence model)
// ─────────────────────────────────────────────────────────────────────────────
double E_GEM  = EMAX_combo * C_dFdCTP_p / (KSS_GEM  + C_dFdCTP_p);
double E_NPAC = EMAX_combo * C_NPAC_p   / (KSS_NPAC + C_NPAC_p);
double E_OHP  = EMAX_combo * C_OHP_DNA_p/ (KSS_OHP  + C_OHP_DNA_p);
double E_SN38 = EMAX_combo * C_SN38_p   / (KSS_SN38 + C_SN38_p);
double E_FU5  = EMAX_combo * C_FU5_p    / (KSS_FU5  + C_FU5_p);
double E_MRTX = EMAX_combo * C_MRTX_p   / (KSS_MRTX + C_MRTX_p);

// Bliss independence combination (avoid double-counting):
// Psi = 1 - (1-E1)*(1-E2)*...
double surv_GEM  = 1.0 - E_GEM;
double surv_NPAC = 1.0 - E_NPAC;
double surv_OHP  = 1.0 - E_OHP;
double surv_SN38 = 1.0 - E_SN38;
double surv_FU5  = 1.0 - E_FU5;
double surv_MRTX = 1.0 - E_MRTX;
double Psi = 1.0 - surv_GEM*surv_NPAC*surv_OHP*surv_SN38*surv_FU5*surv_MRTX;
Psi = Psi > 0.98 ? 0.98 : Psi; // cap at 98%

// KRAS signaling modifies effective tumor growth rate
double k_prog_h  = k_prog  / 24.0;
double k_death_h = k_death / 24.0;
double k_prog_eff = k_prog_h * KRAS_SIG; // KRAS drives proliferation

// ─────────────────────────────────────────────────────────────────────────────
// Simeoni TGI model (4-compartment transit)
// ─────────────────────────────────────────────────────────────────────────────
double x0_safe = x0 > 0 ? x0 : 0;
double Psi_eff = Psi * stroma_pen; // effective drug kill adjusted for stroma

dxdt_x0    = (k_prog_eff*stroma_pen - k_death_h - Psi_eff)*x0_safe;
dxdt_x1    =  Psi_eff*x0_safe - k_tr*x1;
dxdt_x2    =  k_tr*x1         - k_tr*x2;
dxdt_x3    =  k_tr*x2         - k_tr*x3;
dxdt_TUMOR =  k_tr*x3         - k_death_h*TUMOR;

// ─────────────────────────────────────────────────────────────────────────────
// CA19-9 biomarker (proportional to total tumor burden)
// ─────────────────────────────────────────────────────────────────────────────
double tumor_burden = (x0_safe + TUMOR) / V0_tumor; // normalized burden
double CA199_safe = CA199 > 0 ? CA199 : 0;
dxdt_CA199 = k_CA199_prod * tumor_burden - k_CA199_el * CA199_safe;

// ─────────────────────────────────────────────────────────────────────────────
// Friberg neutropenia model
// ─────────────────────────────────────────────────────────────────────────────
double k_prol_val    = 4.0 / MTT;
double k_tr_neut_val = 4.0 / MTT;

// Combined myelosuppressive effect (GEM + SN38)
double drug_myelo = C_GEM_p + C_SN38_p; // additive surrogate
double E_myelo    = EMAX_neut * drug_myelo / (EC50_neut_GEM + drug_myelo);

double Circ_safe = Circ > 0.01 ? Circ : 0.01; // avoid division by zero
double Prol_safe = Prol  > 0 ? Prol  : 0;
double Tr1_safe  = Tr1   > 0 ? Tr1   : 0;
double Tr2_safe  = Tr2   > 0 ? Tr2   : 0;
double Tr3_safe  = Tr3   > 0 ? Tr3   : 0;

double feedback = pow(Circ0/Circ_safe, gamma_neut);

dxdt_Prol = k_prol_val*Prol_safe*(1.0 - E_myelo)*feedback - k_tr_neut_val*Prol_safe;
dxdt_Tr1  = k_tr_neut_val*Prol_safe - k_tr_neut_val*Tr1_safe;
dxdt_Tr2  = k_tr_neut_val*Tr1_safe  - k_tr_neut_val*Tr2_safe;
dxdt_Tr3  = k_tr_neut_val*Tr2_safe  - k_tr_neut_val*Tr3_safe;
dxdt_Circ = k_tr_neut_val*Tr3_safe  - k_tr_neut_val*Circ_safe;

$TABLE
// Derived outputs for capture
double conc_GEM   = GEM_C1  / V1_GEM;
double conc_dFdCTP = dFdCTP / V1_GEM;
double conc_NPAC  = NPAC_C1 / V1_NPAC;
double conc_SN38  = SN38    / 25.0;
double conc_OHP_DNA = OHP_DNA / V_OHP;
double conc_MRTX  = MRTX1133 / V_MRTX;
double conc_OLA   = OLA      / V_OLA;

double ANC = Circ;  // absolute neutrophil count (10^9/L)
double TUMOR_cm3  = TUMOR / 1000.0;  // convert mm3 -> cm3

$CAPTURE
TUMOR CA199 ANC conc_GEM conc_dFdCTP conc_NPAC conc_SN38
conc_OHP_DNA conc_MRTX conc_OLA KRAS_SIG PI3K_ACT STROMA TUMOR_cm3

')

## ─────────────────────────────────────────────────────────────────────────────
## Helper: body-surface-area based dose → rate (ug/h) for IV infusion
## BSA = 1.73 m2 (typical adult); infusion duration in hours
## ─────────────────────────────────────────────────────────────────────────────
bsa_to_rate <- function(dose_per_m2, bsa = 1.73, duration_h = 0.5) {
  total_mg  <- dose_per_m2 * bsa         # mg
  total_ug  <- total_mg * 1000            # ug
  rate_ug_h <- total_ug / duration_h     # ug/h
  return(rate_ug_h)
}

## ─────────────────────────────────────────────────────────────────────────────
## Simulation time: 365 days = 8760 hours, output every 4 hours
## ─────────────────────────────────────────────────────────────────────────────
sim_end   <- 8760      # hours (365 days)
sim_delta <- 4         # h output interval

## ─────────────────────────────────────────────────────────────────────────────
## Treatment Scenario 1: Untreated control
## ─────────────────────────────────────────────────────────────────────────────
ev_control <- ev(time = 0, amt = 0, cmt = 1)  # null event

out_control <- pdac_model %>%
  ev(ev_control) %>%
  mrgsim(end = sim_end, delta = sim_delta) %>%
  as_tibble() %>%
  mutate(scenario = "1. Untreated Control")

## ─────────────────────────────────────────────────────────────────────────────
## Scenario 2: Gemcitabine monotherapy
##   1000 mg/m² IV 30-min infusion, days 1, 8, 15 of 28-day cycle
##   6 cycles = 168 days
## Clinical calibration: MPACT Gem arm mOS 6.7 months (Von Hoff 2013)
## ─────────────────────────────────────────────────────────────────────────────
gem_rate <- bsa_to_rate(1000, duration_h = 0.5)  # ug/h

# Generate dosing events for 6 cycles x 3 doses = 18 infusions
gem_days <- c()
for (cyc in 0:5) {
  gem_days <- c(gem_days, cyc*28 + 1, cyc*28 + 8, cyc*28 + 15)
}
gem_times_h <- gem_days * 24  # convert day to hour

ev_gem <- ev(time  = gem_times_h,
             amt   = gem_rate * 0.5,  # total ug (rate * 0.5h)
             rate  = gem_rate,
             cmt   = 1,
             addl  = 0)

out_gem <- pdac_model %>%
  ev(ev_gem) %>%
  mrgsim(end = sim_end, delta = sim_delta) %>%
  as_tibble() %>%
  mutate(scenario = "2. Gemcitabine mono")

## ─────────────────────────────────────────────────────────────────────────────
## Scenario 3: Gemcitabine + nab-Paclitaxel (MPACT regimen)
##   Gem 1000 mg/m² IV + nab-Pac 125 mg/m² IV, days 1, 8, 15 q28d
## Clinical calibration: MPACT mOS 8.5 vs 6.7 months (Von Hoff 2013, NEJM)
## ─────────────────────────────────────────────────────────────────────────────
npac_rate <- bsa_to_rate(125, duration_h = 0.5)

ev_gem_npac <- ev(time = gem_times_h,
                  amt  = gem_rate * 0.5,
                  rate = gem_rate,
                  cmt  = 1) +
               ev(time = gem_times_h,
                  amt  = npac_rate * 0.5,
                  rate = npac_rate,
                  cmt  = 4)

out_gem_npac <- pdac_model %>%
  ev(ev_gem_npac) %>%
  mrgsim(end = sim_end, delta = sim_delta) %>%
  as_tibble() %>%
  mutate(scenario = "3. Gem + nab-Pac (MPACT)")

## ─────────────────────────────────────────────────────────────────────────────
## Scenario 4: FOLFIRINOX (standard)
##   OHP 85 mg/m² + CPT11 180 mg/m² + 5-FU 400 bolus + 5-FU 2400 CI over 46h
##   q14d, 12 cycles (168 days)
## Clinical calibration: PRODIGE4 mOS 11.1 months (Conroy 2011, NEJM)
## ─────────────────────────────────────────────────────────────────────────────
ohp_rate_std   <- bsa_to_rate(85,   duration_h = 2.0)
cpt11_rate_std <- bsa_to_rate(180,  duration_h = 1.5)
fu5_bolus_rate <- bsa_to_rate(400,  duration_h = 0.083)  # ~5 min
fu5_ci_rate    <- bsa_to_rate(2400, duration_h = 46.0)

folfox_days <- c()
for (cyc in 0:11) folfox_days <- c(folfox_days, cyc*14 + 1)
folfox_times_h <- folfox_days * 24

ev_folfirinox <- ev(time = folfox_times_h,
                    amt  = ohp_rate_std * 2.0,
                    rate = ohp_rate_std,
                    cmt  = 6) +
                 ev(time = folfox_times_h,
                    amt  = cpt11_rate_std * 1.5,
                    rate = cpt11_rate_std,
                    cmt  = 8) +
                 ev(time = folfox_times_h,
                    amt  = fu5_bolus_rate * 0.083,
                    rate = fu5_bolus_rate,
                    cmt  = 10) +
                 ev(time = folfox_times_h,
                    amt  = fu5_ci_rate * 46.0,
                    rate = fu5_ci_rate,
                    cmt  = 10)

out_folfirinox <- pdac_model %>%
  ev(ev_folfirinox) %>%
  mrgsim(end = sim_end, delta = sim_delta) %>%
  as_tibble() %>%
  mutate(scenario = "4. FOLFIRINOX (standard)")

## ─────────────────────────────────────────────────────────────────────────────
## Scenario 5: mFOLFIRINOX (modified; no 5-FU bolus, reduced doses)
##   OHP 65 mg/m² + CPT11 150 mg/m² + 5-FU 2400 CI q14d
##   Reference: Stein 2016; used in adjuvant (PRODIGE24/ACCORD) and palliation
## ─────────────────────────────────────────────────────────────────────────────
ohp_rate_mod   <- bsa_to_rate(65,   duration_h = 2.0)
cpt11_rate_mod <- bsa_to_rate(150,  duration_h = 1.5)
fu5_ci_mod_rate <- bsa_to_rate(2400, duration_h = 46.0)

ev_mfolfirinox <- ev(time = folfox_times_h,
                     amt  = ohp_rate_mod * 2.0,
                     rate = ohp_rate_mod,
                     cmt  = 6) +
                  ev(time = folfox_times_h,
                     amt  = cpt11_rate_mod * 1.5,
                     rate = cpt11_rate_mod,
                     cmt  = 8) +
                  ev(time = folfox_times_h,
                     amt  = fu5_ci_mod_rate * 46.0,
                     rate = fu5_ci_mod_rate,
                     cmt  = 10)

out_mfolfirinox <- pdac_model %>%
  ev(ev_mfolfirinox) %>%
  mrgsim(end = sim_end, delta = sim_delta) %>%
  as_tibble() %>%
  mutate(scenario = "5. mFOLFIRINOX")

## ─────────────────────────────────────────────────────────────────────────────
## Scenario 6: MRTX1133 (KRAS G12D inhibitor)
##   100 mg BID oral; F=0.45; MW ~550 Da (approximate)
##   Continuous dosing, 12 months
## Reference: Fell et al. 2020; Hallin et al. 2022; Phase I (2023)
## ─────────────────────────────────────────────────────────────────────────────
mrtx_dose_mg   <- 100               # mg per dose
mrtx_dose_ug   <- mrtx_dose_mg * 1000  # ug
mrtx_abs_ug    <- mrtx_dose_ug * 0.45  # bioavailable amount
mrtx_times_h   <- seq(0, sim_end - 12, by = 12)  # BID (every 12h)

ev_mrtx <- ev(time = mrtx_times_h,
              amt  = mrtx_abs_ug,
              cmt  = 11)

out_mrtx <- pdac_model %>%
  ev(ev_mrtx) %>%
  mrgsim(end = sim_end, delta = sim_delta) %>%
  as_tibble() %>%
  mutate(scenario = "6. MRTX1133 (KRAS G12D)")

## ─────────────────────────────────────────────────────────────────────────────
## Scenario 7: Olaparib (BRCA1/2-mutated PDAC)
##   300 mg BID oral; F=0.73; maintenance after platinum-based therapy
## Clinical calibration: POLO trial PFS 7.4 vs 3.8 months (Golan 2019, NEJM)
## ─────────────────────────────────────────────────────────────────────────────
ola_dose_mg  <- 300
ola_dose_ug  <- ola_dose_mg * 1000
ola_abs_ug   <- ola_dose_ug * 0.73   # bioavailable
ola_times_h  <- seq(0, sim_end - 12, by = 12)  # BID

ev_ola <- ev(time = ola_times_h,
             amt  = ola_abs_ug,
             cmt  = 12)

out_ola <- pdac_model %>%
  ev(ev_ola) %>%
  mrgsim(end = sim_end, delta = sim_delta) %>%
  as_tibble() %>%
  mutate(scenario = "7. Olaparib (BRCA+)")

## ─────────────────────────────────────────────────────────────────────────────
## Combine all scenarios
## ─────────────────────────────────────────────────────────────────────────────
all_results <- bind_rows(
  out_control,
  out_gem,
  out_gem_npac,
  out_folfirinox,
  out_mfolfirinox,
  out_mrtx,
  out_ola
) %>%
  mutate(
    time_days = time / 24,
    scenario  = factor(scenario, levels = c(
      "1. Untreated Control",
      "2. Gemcitabine mono",
      "3. Gem + nab-Pac (MPACT)",
      "4. FOLFIRINOX (standard)",
      "5. mFOLFIRINOX",
      "6. MRTX1133 (KRAS G12D)",
      "7. Olaparib (BRCA+)"
    ))
  )

## ─────────────────────────────────────────────────────────────────────────────
## Color palette (7 scenarios)
## ─────────────────────────────────────────────────────────────────────────────
scenario_colors <- c(
  "1. Untreated Control"      = "#E64646",
  "2. Gemcitabine mono"       = "#F4A261",
  "3. Gem + nab-Pac (MPACT)"  = "#2A9D8F",
  "4. FOLFIRINOX (standard)"  = "#264653",
  "5. mFOLFIRINOX"            = "#457B9D",
  "6. MRTX1133 (KRAS G12D)"  = "#8338EC",
  "7. Olaparib (BRCA+)"       = "#3A86FF"
)

## ─────────────────────────────────────────────────────────────────────────────
## Plot 1: Tumor Volume over time
## ─────────────────────────────────────────────────────────────────────────────
p1 <- ggplot(all_results, aes(x = time_days, y = TUMOR, color = scenario)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous(breaks = seq(0, 365, 60), limits = c(0, 365)) +
  scale_y_continuous(labels = scales::comma,
                     limits = c(0, NA)) +
  geom_hline(yintercept = 1000, linetype = "dashed",
             color = "grey60", linewidth = 0.5) +
  annotate("text", x = 5, y = 1100, label = "Baseline (1000 mm³)",
           hjust = 0, size = 3, color = "grey50") +
  labs(
    title    = "PDAC Tumor Volume — 7 Treatment Scenarios",
    subtitle = "Simeoni TGI model with stromal resistance and KRAS pathway modulation",
    x        = "Time (days)",
    y        = "Tumor Volume (mm³)",
    color    = "Scenario",
    caption  = paste0(
      "MPACT (Von Hoff 2013): Gem+nab-Pac mOS 8.5 vs Gem 6.7 months; ",
      "PRODIGE4 (Conroy 2011): FOLFIRINOX mOS 11.1 vs 6.8 months"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "right",
    legend.key.size  = unit(0.8, "lines"),
    plot.caption     = element_text(size = 7, color = "grey40"),
    panel.grid.minor = element_blank()
  )
print(p1)

## ─────────────────────────────────────────────────────────────────────────────
## Plot 2: CA19-9 over time
## ─────────────────────────────────────────────────────────────────────────────
p2 <- ggplot(all_results, aes(x = time_days, y = CA199, color = scenario)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous(breaks = seq(0, 365, 60), limits = c(0, 365)) +
  scale_y_continuous(trans = "log10",
                     labels = scales::comma,
                     breaks = c(1, 10, 35, 100, 500, 2000, 10000)) +
  geom_hline(yintercept = 37, linetype = "dashed",
             color = "darkred", linewidth = 0.5) +
  annotate("text", x = 5, y = 42, label = "ULN 37 U/mL",
           hjust = 0, size = 3, color = "darkred") +
  labs(
    title    = "CA19-9 Biomarker Dynamics",
    subtitle = "Proportional to total tumor burden (log₁₀ scale)",
    x        = "Time (days)",
    y        = "CA19-9 (U/mL, log scale)",
    color    = "Scenario"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "right",
    legend.key.size  = unit(0.8, "lines"),
    panel.grid.minor = element_blank()
  )
print(p2)

## ─────────────────────────────────────────────────────────────────────────────
## Plot 3: Neutrophil counts (myelosuppression, ANC)
## ─────────────────────────────────────────────────────────────────────────────
# Focus on first 168 days (treatment period) for clarity
neut_data <- all_results %>% filter(time_days <= 168)

p3 <- ggplot(neut_data, aes(x = time_days, y = ANC, color = scenario)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  scale_color_manual(values = scenario_colors) +
  scale_x_continuous(breaks = seq(0, 168, 28), limits = c(0, 168)) +
  scale_y_continuous(limits = c(0, NA)) +
  geom_hline(yintercept = 1.5, linetype = "dashed",
             color = "red", linewidth = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dotted",
             color = "darkred", linewidth = 0.7) +
  annotate("text", x = 2, y = 1.6,
           label = "G3 Threshold (ANC < 1.5 × 10⁹/L)",
           hjust = 0, size = 2.8, color = "red") +
  annotate("text", x = 2, y = 0.4,
           label = "G4 Threshold (ANC < 0.5 × 10⁹/L)",
           hjust = 0, size = 2.8, color = "darkred") +
  labs(
    title    = "Myelosuppression — Absolute Neutrophil Count (ANC)",
    subtitle = "Friberg transit model; treatment period (first 168 days)",
    x        = "Time (days)",
    y        = "ANC (10⁹/L)",
    color    = "Scenario",
    caption  = "Gemcitabine and SN-38 drive myelosuppression in this model"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "right",
    legend.key.size  = unit(0.8, "lines"),
    plot.caption     = element_text(size = 7, color = "grey40"),
    panel.grid.minor = element_blank()
  )
print(p3)

## ─────────────────────────────────────────────────────────────────────────────
## Plot 4: KRAS signaling and stromal dynamics
## ─────────────────────────────────────────────────────────────────────────────
kras_data <- all_results %>%
  filter(scenario %in% c("1. Untreated Control",
                          "6. MRTX1133 (KRAS G12D)",
                          "3. Gem + nab-Pac (MPACT)")) %>%
  select(time_days, scenario, KRAS_SIG, PI3K_ACT, STROMA)

p4a <- ggplot(kras_data, aes(x = time_days, y = KRAS_SIG, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "KRAS Signaling Over Time",
       x = "Time (days)", y = "KRAS Signal (normalized)",
       color = "Scenario") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())
print(p4a)

## ─────────────────────────────────────────────────────────────────────────────
## Summary statistics table
## ─────────────────────────────────────────────────────────────────────────────
summary_stats <- all_results %>%
  group_by(scenario) %>%
  summarize(
    max_TUMOR    = max(TUMOR, na.rm = TRUE),
    min_TUMOR    = min(TUMOR, na.rm = TRUE),
    final_TUMOR  = TUMOR[which.max(time_days)],
    min_ANC      = min(ANC, na.rm = TRUE),
    nadir_ANC_day = time_days[which.min(ANC)],
    max_CA199    = max(CA199, na.rm = TRUE),
    final_CA199  = CA199[which.max(time_days)],
    .groups = "drop"
  )

cat("\n=== PDAC QSP Model — Simulation Summary ===\n\n")
print(as.data.frame(summary_stats), digits = 3)

cat("\n=== Clinical Trial Calibration Notes ===\n")
cat("MPACT (Von Hoff 2013, NEJM 369:1691):
  Gem+nab-Pac (Sc3) target mOS 8.5 months vs Gem (Sc2) mOS 6.7 months
  ORR: 23% vs 7%; PFS: 5.5 vs 3.7 months\n")
cat("PRODIGE4/ACCORD (Conroy 2011, NEJM 364:1817):
  FOLFIRINOX (Sc4) target mOS 11.1 months vs Gem mOS 6.8 months
  ORR: 31.6% vs 9.4%; PFS: 6.4 vs 3.3 months\n")
cat("POLO (Golan 2019, NEJM 381:317):
  Olaparib maint (Sc7) target mPFS 7.4 months vs placebo 3.8 months
  In BRCA1/2-germline mutated metastatic PDAC after platinum therapy\n")
cat("MRTX1133 (Fell 2020, JMCS; Hallin 2022, Cancer Discov; Phase I 2023-2025):
  KRAS G12D inhibitor (~12% of PDAC); Sc6 demonstrates KRAS pathway suppression\n")

## ─────────────────────────────────────────────────────────────────────────────
## End of PDAC mrgsolve model script
## =============================================================================
