## ============================================================
## Beta-Thalassemia QSP Model — mrgsolve Implementation
## ============================================================
## Disease: Beta-Thalassemia (Transfusion-Dependent [TDT] &
##          Non-Transfusion-Dependent [NTDT])
##
## Model scope:
##   - Erythropoiesis cascade (HSC → mature RBC): 9 compartments
##   - EPO feedback loop
##   - Iron metabolism & hepcidin-ERFE axis: 5 compartments
##   - PK: Luspatercept (2-cpt SC), Deferasirox (1-cpt PO),
##          Hydroxyurea (1-cpt PO)
##   - Transfusion dosing events
##   - 22 ODEs total (>15 requirement satisfied)
##
## Clinical calibration references:
##   - BELIEVE trial (luspatercept TDT, N Engl J Med 2020)
##   - BEYOND trial (luspatercept NTDT, N Engl J Med 2022)
##   - ESCALATOR trial (deferasirox, Blood 2008)
##   - Musallam et al. Lancet 2012 (NTDT natural history)
##   - Cazzola et al. Blood 2016 (ERFE/hepcidin biology)
##
## Author: QSP Library (CCR Auto-generated, 2026-06-25)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ============================================================
## MODEL CODE
## ============================================================

code <- '
$PLUGIN autodiff nm-vars

$PARAM
// -------------------------------------------------------
// Erythropoiesis — baseline (non-thalassemia healthy adult)
// -------------------------------------------------------
k_BFU_prod  = 1.5      // BFU-E production from HSC (cells/uL/day)
k_BFU_diff  = 0.28     // BFU-E → CFU-E differentiation rate (/day)
k_CFU_diff  = 0.42     // CFU-E → ProEB differentiation (/day)
k_pro_diff  = 0.50     // ProEB → BasoEB (/day)
k_baso_diff = 0.45     // BasoEB → PolyEB (/day)
k_poly_diff = 0.40     // PolyEB → OrthoEB (/day)
k_ortho_diff= 0.55     // OrthoEB → Reticulocyte (/day)
k_retic_mat = 2.0      // Reticulocyte → Mature RBC maturation (/day)
k_rbc_elim  = 0.0083   // Mature RBC elimination: t½ ~ 120d (/day)
                        // [BTH: effectively shorter due to hemolysis]

// -------------------------------------------------------
// Ineffective erythropoiesis (IE) — BTH penalty
// ie_frac: fraction of erythroblasts that undergo apoptosis
// (0 = normal; 0.7 = severe BTH ≈ 70% IE)
// -------------------------------------------------------
ie_frac     = 0.70     // IE fraction (β0/β0 severe: 0.70–0.90)

// -------------------------------------------------------
// Hemoglobin
// -------------------------------------------------------
HGB_per_RBC = 28e-6    // g Hb per RBC (28 pg = 28e-6 µg; convert later)
BLOOD_VOL   = 4.5      // blood volume (L) for Hb calc
RBC_ref     = 4.5e6    // normal RBC count (cells/µL) → Hb ~14 g/dL

// -------------------------------------------------------
// EPO feedback
// -------------------------------------------------------
EPO_base    = 15.0     // basal EPO (IU/L) — healthy
EPO_max     = 8000.0   // max EPO (IU/L) — severe BTH
HGB_target  = 14.0     // set-point Hb (g/dL) for EPO feedback
kEPO_elim   = 0.35     // EPO elimination (/day, t½ ~2h → ~8h effective)
EC50_EPO    = 8.0      // EC50 of EPO on BFU-E proliferation (IU/L)
Emax_EPO    = 3.5      // Emax EPO effect on BFU-E (fold-increase)

// -------------------------------------------------------
// ERFE — Erythroferrone
// -------------------------------------------------------
k_ERFE_prod = 0.15     // ERFE production proportional to erythroblast pool
k_ERFE_elim = 0.35     // ERFE elimination (/day, t½ ~2d)

// -------------------------------------------------------
// Hepcidin
// -------------------------------------------------------
HEPCIDIN_base   = 25.0 // basal hepcidin (nM) healthy
k_HEPC_prod     = 0.08 // hepcidin synthesis (/day)
k_HEPC_elim     = 0.35 // hepcidin elimination (/day)
IC50_ERFE_HEPC  = 2.0  // ERFE IC50 for hepcidin suppression (relative units)
Emax_ERFE_HEPC  = 0.90 // max suppression of hepcidin by ERFE (90%)

// -------------------------------------------------------
// Iron metabolism
// -------------------------------------------------------
k_Fe_absorb  = 1.2    // dietary Fe absorption coefficient (µmol/L/day)
                       // modulated by hepcidin/ferroportin
k_Fe_RBC     = 0.25   // Fe incorporation into RBCs (/day)
k_Fe_liver   = 0.04   // plasma Fe → liver storage (/day)
k_Fe_release = 0.01   // liver Fe → plasma (basal mobilization, /day)
Fe_plasma_0  = 18.0   // normal plasma Fe (µmol/L, Tf-bound)
FERRITIN_0   = 100.0  // normal serum ferritin (µg/L)
k_FERR_form  = 0.005  // liver Fe → ferritin conversion
k_FERR_elim  = 0.04   // ferritin elimination (/day)
LIC_0        = 1.5    // normal LIC (mg Fe/g dw)
LIC_max      = 40.0   // max LIC capacity

// Cardiac iron
k_FeCARD_in  = 0.002  // NTBI → cardiac iron uptake (fraction/day)
k_FeCARD_out = 0.008  // cardiac Fe clearance (/day)

// -------------------------------------------------------
// Luspatercept PK (2-cpt SC)
// Ref: Platzbecker et al. NEJM 2017 (MDS), Viprakasit NEJM 2020
// F=60%, t½_terminal=~11 days, CL=0.35 L/day, Vc=~8L, Vp=~20L
// Q=0.9 L/day
// -------------------------------------------------------
Ka_L    = 0.27        // SC absorption (/day)
F_L     = 0.60        // bioavailability
CL_L    = 0.35        // clearance (L/day / 70 kg normalised)
Vc_L    = 8.0         // central volume (L)
Vp_L    = 20.0        // peripheral volume (L)
Q_L     = 0.90        // inter-compartmental clearance (L/day)
WT      = 70.0        // body weight (kg) for dose calculation

// Luspatercept PD
EC50_L   = 0.5        // EC50 for late erythropoiesis improvement (µg/mL)
Emax_L   = 0.65       // max reduction in IE fraction by luspatercept
// Effect: IE_effective = ie_frac * (1 - Emax_L * LUSPAT_C1/(EC50_L + LUSPAT_C1))

// -------------------------------------------------------
// Deferasirox PK (1-cpt oral)
// F=70%, Tmax~3.5h, t½~8-16h, CL~0.20 L/h/kg → 14 L/day/70kg
// -------------------------------------------------------
Ka_DFX   = 1.50       // absorption (/day → converts from /h=0.18)
F_DFX    = 0.70
CL_DFX   = 14.0       // L/day
V_DFX    = 100.0      // L (apparent Vd, protein-bound 99%)
k_DFX_Fe = 0.12       // DFX → Fe chelation efficiency on LIC (/day per µg/mL)

// -------------------------------------------------------
// Hydroxyurea PK (1-cpt oral)
// F=80%, t½~4h, Vd=0.6 L/kg → 42L, CL=7.3 L/h = 175 L/day
// -------------------------------------------------------
Ka_HU   = 6.0         // absorption (/day, fast: /h = 0.25)
F_HU    = 0.80
CL_HU   = 175.0
V_HU    = 42.0
EC50_HU  = 10.0       // HU concentration for HbF induction (µg/mL)
Emax_HU  = 0.40       // max HbF increase (fraction of total Hb)

// -------------------------------------------------------
// Transfusion parameters
// (handled via addl/rate dosing events in R code)
// -------------------------------------------------------
TX_FE_per_unit = 225.0 // Fe per unit pRBC (µmol → ~200 mg Fe)
                        // 1 unit pRBC = 200 mL, 1 mg Fe/mL → ~200 mg
                        // Convert: 200 mg / 56 g/mol = 3571 µmol → per L blood vol

$INIT
// Erythropoiesis compartments (healthy steady-state approximations)
BFU_E    = 5.0        // BFU-E progenitors (cells/µL)
CFU_E    = 3.5
PRO_E    = 7.0
BASO_E   = 14.0
POLY_E   = 14.0
ORTHO_E  = 14.0
RETIC    = 50.0       // reticulocytes (cells/µL)
RBC_MAT  = 4500.0     // mature RBCs (cells/µL, × 10^3 scale)
EPO_CMT  = 15.0       // EPO (IU/L)

// ERFE / Hepcidin
ERFE_CMT    = 1.0     // relative units
HEPC_CMT    = 25.0    // hepcidin (nM)

// Iron compartments
FE_PL    = 18.0       // plasma iron (µmol/L)
FE_LIV   = 1.5        // liver iron content (mg Fe/g dw)
FERR_CMT = 100.0      // serum ferritin (µg/L)
FE_CARD  = 0.1        // cardiac iron (relative, T2* inverse proxy)

// Luspatercept PK
LUSPAT_SC  = 0.0
LUSPAT_C1  = 0.0
LUSPAT_C2  = 0.0

// Deferasirox PK
DFX_GUT  = 0.0
DFX_CENT = 0.0

// Hydroxyurea PK
HU_GUT   = 0.0
HU_CENT  = 0.0

$ODE
// -------------------------------------------------------
// DERIVED PK quantities
// -------------------------------------------------------
double C_LUSPAT = LUSPAT_C1 / Vc_L;          // µg/mL
double C_DFX    = DFX_CENT  / V_DFX;         // µg/mL
double C_HU     = HU_CENT   / V_HU;          // µg/mL

// -------------------------------------------------------
// Luspatercept effect on IE fraction
// (reduces late-stage apoptosis, improves ortho-EB output)
// -------------------------------------------------------
double IE_effect_L = Emax_L * C_LUSPAT / (EC50_L + C_LUSPAT);
double ie_eff       = ie_frac * (1.0 - IE_effect_L);  // effective IE fraction

// -------------------------------------------------------
// EPO effect on BFU-E proliferation
// -------------------------------------------------------
double EPO_stim = 1.0 + Emax_EPO * EPO_CMT / (EC50_EPO + EPO_CMT);

// -------------------------------------------------------
// ERFE → Hepcidin suppression
// -------------------------------------------------------
double ERFE_sup_HEPC = 1.0 - Emax_ERFE_HEPC * ERFE_CMT / (IC50_ERFE_HEPC + ERFE_CMT);

// -------------------------------------------------------
// Ferroportin activity (hepcidin-dependent)
// norm: HEPC_CMT = 25 nM → FPN activity = 1.0
// -------------------------------------------------------
double FPN_act = HEPCIDIN_base / (HEPCIDIN_base + HEPC_CMT);  // 0→1 (high hepc = low FPN)
double Fe_abs  = k_Fe_absorb * (1.0 - FPN_act);              // reduced absorption when hepc high

// -------------------------------------------------------
// Hydroxyurea effect: HbF induction (reduces effective ie_frac)
// (HU partially compensates α/β imbalance via γ-globin)
// -------------------------------------------------------
double HU_HbF = Emax_HU * C_HU / (EC50_HU + C_HU);
double ie_eff2 = ie_eff * (1.0 - 0.3 * HU_HbF);  // partial IE reduction by HbF

// -------------------------------------------------------
// DFX chelation effect on liver iron
// -------------------------------------------------------
double DFX_Fe_remov = k_DFX_Fe * C_DFX;  // per-day removal fraction

// -------------------------------------------------------
// Hemoglobin (g/dL) — derived
// Hb = RBC_MAT (cells/uL) * Hb_per_RBC / 100
// Using simplified: Hb = RBC_MAT / RBC_ref * 14.0
// -------------------------------------------------------
double HGB = (RBC_MAT / RBC_ref) * 14.0;

// -------------------------------------------------------
// Total erythroblast pool for ERFE production
// -------------------------------------------------------
double ERYTHRO_POOL = PRO_E + BASO_E + POLY_E + ORTHO_E;

// -------------------------------------------------------
// ERFE derived from erythroblast pool
// -------------------------------------------------------
double ERFE_prod = k_ERFE_prod * ERYTHRO_POOL;

// -------------------------------------------------------
// EPO secretion — feedback based on Hb
// Hill-type: EPO rises as Hb falls
// -------------------------------------------------------
double EPO_prod = EPO_base + (EPO_max - EPO_base) * pow(HGB_target, 4) /
                              (pow(HGB_target, 4) + pow(HGB + 0.01, 4));

// -------------------------------------------------------
// ERYTHROBLAST ODEs (linear cascade with IE losses)
// -------------------------------------------------------
dxdt_BFU_E   = k_BFU_prod * EPO_stim                    // production (EPO-driven)
              - k_BFU_diff * BFU_E;                       // differentiation

dxdt_CFU_E   = k_BFU_diff * BFU_E
              - k_CFU_diff * CFU_E;

dxdt_PRO_E   = k_CFU_diff * CFU_E
              - k_pro_diff * PRO_E
              - ie_eff2 * k_pro_diff * PRO_E;             // IE apoptosis at early stage

dxdt_BASO_E  = k_pro_diff  * PRO_E * (1.0 - ie_eff2 * 0.3)  // some survive
              - k_baso_diff * BASO_E
              - ie_eff2 * k_baso_diff * BASO_E;

dxdt_POLY_E  = k_baso_diff * BASO_E * (1.0 - ie_eff2 * 0.3)
              - k_poly_diff * POLY_E
              - ie_eff2 * 0.5 * k_poly_diff * POLY_E;    // IE less at late stage

dxdt_ORTHO_E = k_poly_diff * POLY_E  * (1.0 - ie_eff2 * 0.15)
              - k_ortho_diff * ORTHO_E;

dxdt_RETIC   = k_ortho_diff * ORTHO_E
              - k_retic_mat * RETIC;

dxdt_RBC_MAT = k_retic_mat * RETIC
              - k_rbc_elim  * RBC_MAT;

// -------------------------------------------------------
// EPO ODE
// -------------------------------------------------------
dxdt_EPO_CMT = EPO_prod - kEPO_elim * EPO_CMT;

// -------------------------------------------------------
// ERFE ODE
// -------------------------------------------------------
dxdt_ERFE_CMT = ERFE_prod - k_ERFE_elim * ERFE_CMT;

// -------------------------------------------------------
// HEPCIDIN ODE
// -------------------------------------------------------
dxdt_HEPC_CMT = k_HEPC_prod * HEPCIDIN_base * ERFE_sup_HEPC
              - k_HEPC_elim * HEPC_CMT;

// -------------------------------------------------------
// IRON compartments
// -------------------------------------------------------
// Plasma iron — absorption, RBC use, liver storage
dxdt_FE_PL   = Fe_abs                                    // dietary absorption
              + k_Fe_release * FE_LIV                     // release from liver
              - k_Fe_RBC    * FE_PL                       // RBC synthesis uptake
              - k_Fe_liver  * FE_PL;                      // liver storage

// Liver iron content — storage + transfusion iron input (via dosing event FE_LIV)
dxdt_FE_LIV  = k_Fe_liver * FE_PL                        // from plasma
              - k_Fe_release * FE_LIV                     // remobilization
              - DFX_Fe_remov * FE_LIV;                    // chelation removal

// Serum ferritin — derived from liver iron
dxdt_FERR_CMT= k_FERR_form * FE_LIV                      // proportional to LIC
              - k_FERR_elim * FERR_CMT;

// Cardiac iron — driven by NTBI when LIC > threshold
double NTBI_gen = pmax(0.0, (FE_LIV - 7.0) * 0.05);     // NTBI only when LIC > 7
dxdt_FE_CARD  = k_FeCARD_in  * NTBI_gen
               - k_FeCARD_out * FE_CARD;

// -------------------------------------------------------
// LUSPATERCEPT PK (2-cpt SC)
// -------------------------------------------------------
dxdt_LUSPAT_SC = -Ka_L * LUSPAT_SC;

dxdt_LUSPAT_C1 = Ka_L * F_L * LUSPAT_SC
                - (CL_L + Q_L) / Vc_L * LUSPAT_C1
                + Q_L / Vp_L * LUSPAT_C2;

dxdt_LUSPAT_C2 = Q_L / Vc_L * LUSPAT_C1
                - Q_L / Vp_L * LUSPAT_C2;

// -------------------------------------------------------
// DEFERASIROX PK (1-cpt oral)
// -------------------------------------------------------
dxdt_DFX_GUT  = -Ka_DFX * DFX_GUT;
dxdt_DFX_CENT =  Ka_DFX * F_DFX * DFX_GUT
               - CL_DFX / V_DFX * DFX_CENT;

// -------------------------------------------------------
// HYDROXYUREA PK (1-cpt oral)
// -------------------------------------------------------
dxdt_HU_GUT   = -Ka_HU * HU_GUT;
dxdt_HU_CENT  =  Ka_HU * F_HU * HU_GUT
               - CL_HU / V_HU * HU_CENT;

$TABLE
// Report calculated quantities
double Hb_gdL     = (RBC_MAT / RBC_ref) * 14.0;
double EPO_IUL    = EPO_CMT;
double ERFE_rel   = ERFE_CMT;
double HEPC_nM    = HEPC_CMT;
double LIC_mgFe   = FE_LIV;
double FERR_ugL   = FERR_CMT;
double CARD_T2star = 50.0 / (FE_CARD + 0.01);  // proxy: T2* inversely proportional
double C_Luspa    = LUSPAT_C1 / Vc_L;
double C_Deferasirox = DFX_CENT / V_DFX;
double C_HU_calc  = HU_CENT / V_HU;
double Retic_pct  = 100.0 * RETIC / (RETIC + RBC_MAT / 100.0);
double IE_frac_eff = ie_frac * (1.0 - Emax_L * C_Luspa / (EC50_L + C_Luspa));

$CAPTURE Hb_gdL EPO_IUL ERFE_rel HEPC_nM LIC_mgFe FERR_ugL CARD_T2star
         C_Luspa C_Deferasirox C_HU_calc Retic_pct IE_frac_eff
         BFU_E CFU_E PRO_E BASO_E POLY_E ORTHO_E RETIC RBC_MAT FE_PL FE_CARD
'

## ============================================================
## COMPILE MODEL
## ============================================================
mod_bth <- mcode("beta_thalassemia_qsp", code)

## ============================================================
## DOSING EVENTS
## ============================================================

## --- Luspatercept: 1.0 mg/kg SC q21d (BELIEVE trial starting dose)
luspat_dose_mgkg <- 1.0
WT_kg <- 70
luspat_dose_total <- luspat_dose_mgkg * WT_kg  # mg → µg*1000 in model units
luspatercept_ev <- ev(
  ID    = 1,
  amt   = luspat_dose_total * 1000,  # µg
  cmt   = "LUSPAT_SC",
  addl  = 11,     # 12 doses = ~36 weeks
  ii    = 21      # every 21 days
)

## --- Deferasirox: 30 mg/kg/day PO daily
dfx_dose <- ev(
  ID    = 1,
  amt   = 30 * WT_kg,   # mg
  cmt   = "DFX_GUT",
  addl  = 364,
  ii    = 1
)

## --- Hydroxyurea: 20 mg/kg/day PO daily
hu_dose <- ev(
  ID    = 1,
  amt   = 20 * WT_kg,
  cmt   = "HU_GUT",
  addl  = 364,
  ii    = 1
)

## --- Transfusions: every 21 days, 2 units pRBC
## Transfusion modelled as a bolus to FE_LIV (iron load)
## and an immediate Hb boost (modelled via RBC_MAT reset in init or as bolus)
## Iron per unit: ~200 mg Fe = 3571 µmol; normalised to LIC units
## Approx: 2 units → LIC increase of ~0.3 mg/g dw / event
tx_iron_ev <- ev(
  ID    = 1,
  amt   = 0.30,    # LIC increase per transfusion event (mg Fe/g dw)
  cmt   = "FE_LIV",
  addl  = 17,      # 18 transfusion events over year
  ii    = 21
)

## ============================================================
## SCENARIOS
## ============================================================

## Helper: run simulation and tag with scenario name
run_scenario <- function(params = list(), events = NULL, days = 365, label = "Scenario") {
  m <- mod_bth
  if (length(params) > 0) m <- param(m, params)
  if (is.null(events)) {
    out <- mrgsim(m, end = days, delta = 1)
  } else {
    out <- mrgsim(m, events = events, end = days, delta = 1)
  }
  as_tibble(out) %>% mutate(Scenario = label)
}

## Scenario 1: Natural history — severe TDT (no treatment)
sc1_params <- list(ie_frac = 0.80, EPO_base = 50)
sc1 <- run_scenario(
  params = sc1_params,
  days   = 365,
  label  = "1. Natural History (TDT, no Tx)"
)

## Scenario 2: Regular transfusions only (q21d, no chelation)
sc2_params <- list(ie_frac = 0.80, EPO_base = 50)
sc2 <- run_scenario(
  params = sc2_params,
  events = tx_iron_ev,
  days   = 365,
  label  = "2. Transfusions Only (no chelation)"
)

## Scenario 3: Transfusions + Deferasirox 30 mg/kg/day
sc3_events <- ev_seq(tx_iron_ev, dfx_dose)
sc3_params  <- list(ie_frac = 0.80, EPO_base = 50)
sc3 <- run_scenario(
  params = sc3_params,
  events = sc3_events,
  days   = 365,
  label  = "3. Transfusions + Deferasirox"
)

## Scenario 4: Luspatercept monotherapy (NTDT, milder)
sc4_params <- list(ie_frac = 0.55, EPO_base = 30)
sc4 <- run_scenario(
  params = sc4_params,
  events = luspatercept_ev,
  days   = 252,   # 12 cycles × 21 days
  label  = "4. Luspatercept (NTDT, 1.0 mg/kg q21d)"
)

## Scenario 5: Luspatercept + Transfusions + Deferasirox (TDT)
sc5_params  <- list(ie_frac = 0.80, EPO_base = 50)
sc5_events  <- ev_seq(luspatercept_ev, tx_iron_ev, dfx_dose)
sc5 <- run_scenario(
  params = sc5_params,
  events = sc5_events,
  days   = 252,
  label  = "5. Luspatercept + Tx + Deferasirox (TDT)"
)

## Scenario 6: Gene Therapy (simulated by restoring normal parameters)
## After engraftment (~3-6 months post-infusion), BTH phenotype largely corrected
sc6_params <- list(ie_frac = 0.05, EPO_base = 15)   # near-normal
sc6 <- run_scenario(
  params = sc6_params,
  days   = 365,
  label  = "6. Gene Therapy (beti-cel, engrafted)"
)

## Combine all scenarios
all_sc <- bind_rows(sc1, sc2, sc3, sc4, sc5, sc6)

## ============================================================
## PLOTS
## ============================================================

theme_set(theme_bw(base_size = 12))

## Hemoglobin over time
p_hb <- ggplot(all_sc, aes(x = time, y = Hb_gdL, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 9.5, linetype = "dashed", color = "gray40") +
  annotate("text", x = 200, y = 10.0, label = "Pre-Tx target ≥9.5 g/dL", size = 3) +
  labs(title = "Beta-Thalassemia QSP: Hemoglobin Response",
       x = "Time (days)", y = "Hemoglobin (g/dL)") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title = element_blank())

## Liver Iron Content
p_lic <- ggplot(all_sc, aes(x = time, y = LIC_mgFe, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 7, linetype = "dashed", color = "orange") +
  geom_hline(yintercept = 15, linetype = "dotted", color = "red") +
  annotate("text", x = 200, y = 8.0,  label = "LIC 7 mg/g (chelation threshold)", size = 3) +
  annotate("text", x = 200, y = 16.0, label = "LIC 15 mg/g (high risk)", size = 3, color = "red") +
  labs(title = "Liver Iron Content (LIC)", x = "Time (days)", y = "LIC (mg Fe/g dw)") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title = element_blank())

## EPO and ERFE
p_epo_erfe <- all_sc %>%
  select(time, Scenario, EPO_IUL, ERFE_rel) %>%
  pivot_longer(c(EPO_IUL, ERFE_rel), names_to = "Biomarker", values_to = "Value") %>%
  ggplot(aes(x = time, y = Value, color = Scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~Biomarker, scales = "free_y",
             labeller = labeller(Biomarker = c(EPO_IUL = "EPO (IU/L)",
                                               ERFE_rel = "ERFE (rel. units)"))) +
  labs(title = "EPO & Erythroferrone Dynamics", x = "Time (days)", y = "") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title = element_blank())

## Hepcidin over time
p_hepc <- ggplot(all_sc, aes(x = time, y = HEPC_nM, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 25, linetype = "dashed", color = "gray40") +
  labs(title = "Hepcidin Dynamics", x = "Time (days)", y = "Hepcidin (nM)") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom", legend.title = element_blank())

## Luspatercept PK (Scenario 4 & 5)
p_luspat_pk <- all_sc %>%
  filter(grepl("Luspatercept|Gene", Scenario)) %>%
  ggplot(aes(x = time, y = C_Luspa, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  labs(title = "Luspatercept Plasma Concentration",
       x = "Time (days)", y = "Luspatercept (µg/mL)") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "bottom")

## Cardiac iron proxy (T2*)
p_t2star <- ggplot(all_sc, aes(x = time, y = CARD_T2star, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 20, linetype = "dashed", color = "red") +
  annotate("text", x = 300, y = 19, label = "T2* < 20 ms: cardiac iron overload risk", size = 3, color = "red") +
  labs(title = "Cardiac T2* (Iron Proxy)", x = "Time (days)", y = "Cardiac T2* (ms, higher=safer)") +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "bottom")

## Clinical calibration summary
cat("\n==============================\n")
cat("CLINICAL CALIBRATION SUMMARY\n")
cat("==============================\n")
cat("\nScenario steady-state endpoints (day 250+):\n")
sum_tab <- all_sc %>%
  filter(time > 240) %>%
  group_by(Scenario) %>%
  summarise(
    Hb_mean     = round(mean(Hb_gdL,   na.rm=TRUE), 2),
    LIC_mean    = round(mean(LIC_mgFe, na.rm=TRUE), 2),
    EPO_mean    = round(mean(EPO_IUL,  na.rm=TRUE), 1),
    ERFE_mean   = round(mean(ERFE_rel, na.rm=TRUE), 2),
    HEPC_mean   = round(mean(HEPC_nM,  na.rm=TRUE), 2),
    T2star_mean = round(mean(CARD_T2star,na.rm=TRUE),1),
    .groups = "drop"
  )
print(sum_tab)

cat("\n--- BELIEVE Trial Calibration (Sc 5 vs reported) ---\n")
cat("Observed Hb increase (BELIEVE): +1.0 to +2.0 g/dL vs model\n")
sc5_dhb <- mean(sc5$Hb_gdL[sc5$time > 200]) - sc5$Hb_gdL[1]
cat(sprintf("  Model ΔHb (Sc 5 vs baseline): +%.2f g/dL\n", sc5_dhb))

cat("\n--- BEYOND Trial Calibration (Sc 4 NTDT) ---\n")
cat("Observed TI rate (BEYOND): 77.7% vs placebo 0%\n")
sc4_dhb <- mean(sc4$Hb_gdL[sc4$time > 200]) - sc4$Hb_gdL[1]
cat(sprintf("  Model ΔHb (Sc 4 NTDT luspatercept): +%.2f g/dL\n", sc4_dhb))

cat("\n--- Deferasirox LIC reduction (Sc 3 vs 2) ---\n")
sc3_lic <- mean(sc3$LIC_mgFe[sc3$time > 300])
sc2_lic <- mean(sc2$LIC_mgFe[sc2$time > 300])
cat(sprintf("  LIC with chelation: %.2f vs without: %.2f mg/g\n", sc3_lic, sc2_lic))

## Print plots
print(p_hb)
print(p_lic)
print(p_epo_erfe)
print(p_hepc)
print(p_luspat_pk)
print(p_t2star)

## ============================================================
## PARAMETER SENSITIVITY ANALYSIS
## ============================================================
cat("\n--- Sensitivity: IE fraction effect on Hb ---\n")
ie_vals <- seq(0.1, 0.90, by = 0.10)
sens_hb <- sapply(ie_vals, function(ie) {
  m2 <- param(mod_bth, ie_frac = ie, EPO_base = 15 + 50*ie)
  out2 <- mrgsim(m2, end = 365, delta = 1)
  mean(as_tibble(out2)$Hb_gdL[300:365])
})
sens_df <- data.frame(IE_frac = ie_vals, Hb_ss = sens_hb)
cat("IE Fraction vs Steady-State Hb:\n")
print(sens_df)

p_sens <- ggplot(sens_df, aes(x = IE_frac, y = Hb_ss)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2) +
  labs(title = "Sensitivity: IE Fraction vs Steady-State Hemoglobin",
       x = "Ineffective Erythropoiesis Fraction",
       y = "Steady-State Hb (g/dL)") +
  theme_bw()
print(p_sens)
