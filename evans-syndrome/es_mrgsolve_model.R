## ============================================================
## Evans Syndrome (ES) – mrgsolve QSP Model
## Combined AIHA + Immune Thrombocytopenia
## ============================================================
## Author : Claude Code (CCR)
## Date   : 2026-06-20
## References:
##   - Audia et al. (2020) Blood Rev. Pathophysiology of AIHA.
##   - Michel et al. (2021) Lancet Haematol. Evans syndrome review.
##   - Weksler (2022) Am J Hematol. Rituximab in Evans syndrome.
##   - Fattizzo et al. (2023) Haematologica. ES treatment.
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ============================================================
## 1.  mrgsolve Model Code Block
## ============================================================
ES_code <- '
$PROB Evans Syndrome QSP Model (AIHA + ITP)

$PARAM
// ---- Disease baseline parameters ----
B_RBC0    = 500      // Anti-RBC B cell baseline (cells/uL)
B_Plt0    = 500      // Anti-Plt B cell baseline
Ab_RBC0   = 2.0      // Anti-RBC IgG baseline (mg/L; normal ~0.1)
Ab_Plt0   = 2.0      // Anti-Plt IgG baseline (mg/L)
Treg0     = 80       // Treg cell baseline (cells/uL)
Treg_min  = 20       // Minimum Treg (disease nadir)

// ---- Erythrocyte / Hemoglobin ----
RBC0      = 5.0      // Baseline RBC (x10^12/L)
Hgb0      = 14.0     // Baseline Hgb (g/dL)
Hgb_per_RBC = 2.8   // Hgb g/dL per x10^12 RBC/L
Retic0    = 1.5      // Reticulocyte % baseline
BM_Ery0   = 100      // BM erythroid progenitor pool (AU)

// ---- Platelet ----
Plt0      = 250      // Baseline platelet count (x10^9/L)
MK0       = 100      // Megakaryocyte pool baseline (AU)

// ---- Kinetic parameters – B cells ----
kprol_B   = 0.025    // B cell proliferation rate (1/day)
kdeg_B    = 0.020    // B cell apoptosis rate (1/day)
kB_Ab_RBC = 0.004    // B cell -> Ab_RBC production (mg/L/day per cell/uL)
kB_Ab_Plt = 0.004    // B cell -> Ab_Plt production
kclear_Ab = 0.10     // Antibody catabolism rate (1/day)

// ---- Kinetic parameters – Treg ----
kprol_Treg = 0.015   // Treg proliferation (1/day)
kdeg_Treg  = 0.012   // Treg apoptosis (1/day)
ISuppression = 0.4   // Max immunosuppression by Treg (fraction)
EC50_Treg  = 60      // Treg EC50 for B cell suppression (cells/uL)

// ---- Complement ----
kform_C3b  = 0.05    // C3b formation rate constant (1/mg/L/day)
kclear_C3b = 0.30    // C3b clearance from RBC surface (1/day)
C3b0       = 0.1     // Baseline C3b (AU)

// ---- RBC / hemolysis ----
kops_RBC   = 0.015   // Opsonization rate constant (1/mg/L/day)
kphago_RBC = 0.30    // Phagocytosis rate of opsonized RBC (1/day)
kdeg_RBC   = 0.012   // Normal RBC senescence (1/day; t1/2 ~90 days)
kinput_RBC = 0.060   // Baseline RBC influx from BM (x10^12/L/day)
EPO_EC50   = 7.0     // Hgb (g/dL) at which EPO = 2x baseline
kEPO_stim  = 3.0     // Maximum EPO stimulation fold

// ---- Reticulocyte / BM erythroid ----
kprol_BM_Ery = 0.050 // BM erythroid progenitor proliferation (1/day)
kdeg_BM_Ery  = 0.040 // BM progenitor death/export (1/day)
kmat_Retic   = 0.50  // Reticulocyte maturation rate (1/day)

// ---- Platelet / MK ----
kops_Plt   = 0.018   // Platelet opsonization (1/mg/L/day)
kphago_Plt = 0.35    // Phagocytosis of opsonized Plt (1/day)
kdeg_Plt   = 0.043   // Normal Plt lifespan (1/day; t1/2~10d)
kinput_Plt = 10.8    // Baseline Plt production (x10^9/L/day)
kprol_MK   = 0.040   // MK proliferation (1/day)
kdeg_MK    = 0.030   // MK turnover (1/day)
krel_Plt   = 0.108   // Plt release rate from MK pool (AU/day)

// ---- Drug indicators (0/1 or dose) ----
PRED_dose  = 0       // Prednisone dose (mg/kg)
IVIG_dose  = 0       // IVIg dose flag (0/1)
RTX_dose   = 0       // Rituximab dose flag (0/1)
MMF_dose   = 0       // MMF dose (g/day)
ELT_dose   = 0       // Eltrombopag dose (mg/day)
SIRO_dose  = 0       // Sirolimus dose (mg/day)
SPLEN      = 0       // Splenectomy (0/1)

// ---- Drug PK parameters ----
// Prednisone -> Prednisolone
CL_Pred    = 15.0    // CL (L/day)
Vc_Pred    = 45.0    // Central Vd (L)
F_Pred     = 0.82    // Bioavailability
ka_Pred    = 8.0     // Absorption rate (1/day)

// IVIg
CL_IVIg    = 0.30    // CL (L/day) - FcRn mediated
Vc_IVIg    = 3.5     // Vd (L)

// Rituximab
CL_Rtx     = 0.33    // CL (L/day; non-linear Michaelis)
Vc_Rtx     = 4.0     // Vd (L)
Km_Rtx     = 5.0     // Michaelis constant for target-mediated

// MPA (active form of MMF)
CL_MPA     = 25.0    // CL (L/day)
Vc_MPA     = 15.0    // Vd (L)
F_MPA      = 0.94    // Bioavailability
ka_MPA     = 6.0     // Absorption (1/day)

// Eltrombopag
CL_Elt     = 12.0    // CL (L/day)
Vc_Elt     = 20.0    // Vd (L)
F_Elt      = 0.52    // Bioavailability
ka_Elt     = 3.0     // Absorption (1/day)

// Sirolimus
CL_Siro    = 15.0    // CL (L/day) (whole blood)
Vc_Siro    = 480.0   // Vd – large due to RBC partitioning
F_Siro     = 0.14    // Bioavailability
ka_Siro    = 2.0     // Absorption (1/day)

// ---- Drug PD parameters ----
// Corticosteroid
Emax_CS_B  = 0.70    // Max B cell apoptosis by CS (fraction)
EC50_CS_B  = 0.08    // Prednisolone EC50 for B cell (mg/L)
Emax_CS_Mac = 0.60   // Max FcgR downreg
EC50_CS_Mac = 0.06

// IVIg
Emax_IVIg_FcR = 0.75 // Max FcgR blockade
EC50_IVIg     = 10.0 // IVIg EC50 (mg/L)
kcat_IgG      = 0.05 // Extra IgG catabolism by IVIg (1/day)

// Rituximab
Emax_Rtx_B  = 0.95   // Max B cell depletion
EC50_Rtx    = 1.0    // Rtx EC50 (mg/L)

// MMF/MPA
Emax_MPA_Lym = 0.65  // Max lymphocyte antiproliferative
EC50_MPA     = 1.5   // MPA EC50 (mg/L)

// Eltrombopag
Emax_Elt_MK = 2.5    // Max-fold MK stimulation
EC50_Elt    = 0.5    // Elt EC50 (mg/L)

// Sirolimus
Emax_Siro_Treg = 2.0 // Max Treg fold expansion
EC50_Siro_Treg = 1.5 // Sirolimus EC50 (ng/mL)

// Splenectomy
Splen_Phago_red = 0.80 // Phagocytosis reduction after splenectomy

$CMT
// Drug PK compartments
DEPOT_PRED  // Prednisone depot (oral)
CENT_PRED   // Prednisolone central (mg)
DEPOT_MPA   // MMF depot (oral)
CENT_MPA    // MPA central (mg)
DEPOT_ELT   // Eltrombopag depot (oral)
CENT_ELT    // Eltrombopag central (mg)
DEPOT_SIRO  // Sirolimus depot (oral)
CENT_SIRO   // Sirolimus central (mg)
CENT_IVIG   // IVIg central (mg)
CENT_RTX    // Rituximab central (mg)

// Immunology
B_RBC       // Anti-RBC B cells (cells/uL)
AB_RBC      // Anti-RBC IgG (mg/L)
B_PLT       // Anti-Plt B cells (cells/uL)
AB_PLT      // Anti-Plt IgG (mg/L)
TREG        // Regulatory T cells (cells/uL)
C3B         // Complement C3b on RBC (AU)

// Hematology – RBC compartment
ORBC        // Opsonized RBC (x10^12/L)
RBC_CIRC    // Circulating (non-opsonized) RBC (x10^12/L)
RETIC       // Reticulocyte (% of RBC)
BM_ERY      // BM erythroid precursors (AU)

// Hematology – Platelet compartment
OPLT        // Opsonized platelets (x10^9/L)
PLT_CIRC    // Circulating platelets (x10^9/L)
MK          // Megakaryocyte pool (AU)

$MAIN
// Initial conditions
DEPOT_PRED_0  = 0;
CENT_PRED_0   = 0;
DEPOT_MPA_0   = 0;
CENT_MPA_0    = 0;
DEPOT_ELT_0   = 0;
CENT_ELT_0    = 0;
DEPOT_SIRO_0  = 0;
CENT_SIRO_0   = 0;
CENT_IVIG_0   = 0;
CENT_RTX_0    = 0;

B_RBC_0    = B_RBC0;
AB_RBC_0   = Ab_RBC0;
B_PLT_0    = B_Plt0;
AB_PLT_0   = Ab_Plt0;
TREG_0     = Treg0;
C3B_0      = C3b0;
ORBC_0     = 0.1;
RBC_CIRC_0 = RBC0;
RETIC_0    = Retic0;
BM_ERY_0   = BM_Ery0;
OPLT_0     = 5.0;
PLT_CIRC_0 = Plt0;
MK_0       = MK0;

$ODE
// ===== Drug PK =====
// --- Prednisone/Prednisolone ---
double dose_PRED_daily = PRED_dose * 70.0; // mg/day (70kg patient)
double kin_PRED  = F_Pred * dose_PRED_daily * ka_Pred;  // mg/day -> depot
double C_Pred    = CENT_PRED / Vc_Pred;    // mg/L (prednisolone)
dxdt_DEPOT_PRED  = kin_PRED - ka_Pred * DEPOT_PRED;
dxdt_CENT_PRED   = ka_Pred * DEPOT_PRED - CL_Pred * C_Pred;

// --- MPA (MMF) ---
double dose_MPA_daily = MMF_dose * 1000.0; // g -> mg
double kin_MPA = F_MPA * dose_MPA_daily * ka_MPA;
double C_MPA   = CENT_MPA / Vc_MPA;
dxdt_DEPOT_MPA = kin_MPA - ka_MPA * DEPOT_MPA;
dxdt_CENT_MPA  = ka_MPA * DEPOT_MPA - CL_MPA * C_MPA;

// --- Eltrombopag ---
double kin_ELT = F_Elt * ELT_dose * ka_Elt;
double C_Elt   = CENT_ELT / Vc_Elt;
dxdt_DEPOT_ELT = kin_ELT - ka_Elt * DEPOT_ELT;
dxdt_CENT_ELT  = ka_Elt * DEPOT_ELT - CL_Elt * C_Elt;

// --- Sirolimus ---
double kin_SIRO = F_Siro * SIRO_dose * ka_Siro;
double C_Siro   = CENT_SIRO / Vc_Siro;  // ng/mL whole blood (approx)
dxdt_DEPOT_SIRO = kin_SIRO - ka_Siro * DEPOT_SIRO;
dxdt_CENT_SIRO  = ka_Siro * DEPOT_SIRO - CL_Siro * C_Siro;

// --- IVIg (IV bolus handled by event) ---
double C_IVIg = CENT_IVIG / Vc_IVIg;
// FcRn saturation: at high IVIg, catabolism accelerates (non-linear)
double IVIg_cat_extra = kcat_IgG * C_IVIg / (EC50_IVIg + C_IVIg);
dxdt_CENT_IVIG = -CL_IVIg * C_IVIg;

// --- Rituximab (IV doses via event) ---
double C_Rtx = CENT_RTX / Vc_Rtx;
// Target-mediated disposition (simplified Michaelis-Menten)
dxdt_CENT_RTX = -CL_Rtx * C_Rtx - (CL_Rtx * Km_Rtx * C_Rtx) / (Km_Rtx + C_Rtx);

// ===== Drug PD effects =====
// Corticosteroid effects
double E_CS_B   = Emax_CS_B   * C_Pred / (EC50_CS_B   + C_Pred);  // B apoptosis
double E_CS_Mac = Emax_CS_Mac * C_Pred / (EC50_CS_Mac + C_Pred);  // FcgR downreg

// IVIg effects
double E_IVIg_FcR = Emax_IVIg_FcR * C_IVIg / (EC50_IVIg + C_IVIg); // FcgR blockade
double E_IVIg_cat = IVIg_cat_extra;                                    // Ab catabolism

// Rituximab effects
double E_Rtx_B  = Emax_Rtx_B * C_Rtx / (EC50_Rtx + C_Rtx);  // B depletion

// MPA effects
double E_MPA_Lym = Emax_MPA_Lym * C_MPA / (EC50_MPA + C_MPA);  // Lymphoproliferation↓

// Eltrombopag MK stimulation
double E_Elt_MK = 1.0 + (Emax_Elt_MK - 1.0) * C_Elt / (EC50_Elt + C_Elt);

// Sirolimus – Treg expansion
double E_Siro_Treg = 1.0 + (Emax_Siro_Treg - 1.0) * C_Siro / (EC50_Siro_Treg + C_Siro);

// Splenectomy – reduce phagocytosis
double Splen_factor = 1.0 - SPLEN * Splen_Phago_red;

// Combined FcgR suppression (CS + IVIg + splenectomy)
double FcgR_RBC = Splen_factor * (1.0 - E_CS_Mac) * (1.0 - E_IVIg_FcR);
double FcgR_Plt = Splen_factor * (1.0 - E_CS_Mac) * (1.0 - E_IVIg_FcR);

// ===== Treg dynamics =====
// Treg suppression of B cells: Hill function
double Treg_suppression = ISuppression * TREG / (EC50_Treg + TREG);
// Sirolimus expands Treg
dxdt_TREG = kprol_Treg * E_Siro_Treg * Treg0 - kdeg_Treg * TREG;

// ===== B cell dynamics =====
// Net B cell proliferation reduced by Treg suppression, Rtx, and MPA
double B_net_suppression = (1.0 - Treg_suppression) * (1.0 - E_Rtx_B) * (1.0 - E_MPA_Lym);
dxdt_B_RBC = kprol_B * B_net_suppression * B_RBC0 - kdeg_B * B_RBC - E_CS_B * kdeg_B * B_RBC;
dxdt_B_PLT = kprol_B * B_net_suppression * B_Plt0 - kdeg_B * B_PLT - E_CS_B * kdeg_B * B_PLT;

// ===== Autoantibody dynamics =====
// Ab catabolism increased by IVIg (FcRn saturation)
double kclear_Ab_RBC = kclear_Ab + E_IVIg_cat;
double kclear_Ab_Plt = kclear_Ab + E_IVIg_cat;
dxdt_AB_RBC = kB_Ab_RBC * B_RBC * (1.0 - E_CS_B) - kclear_Ab_RBC * AB_RBC;
dxdt_AB_PLT = kB_Ab_Plt * B_PLT * (1.0 - E_CS_B) - kclear_Ab_Plt * AB_PLT;

// ===== Complement C3b on RBC =====
// Classical pathway activated proportional to Ab_RBC
double C3b_form = kform_C3b * AB_RBC * RBC_CIRC;
dxdt_C3B = C3b_form - kclear_C3b * C3B;

// ===== RBC hemolysis =====
// Opsonization rate depends on Ab_RBC and C3b
double kops_eff  = kops_RBC * (AB_RBC + 0.5 * C3B);
double phago_RBC = kphago_RBC * FcgR_RBC * ORBC;
dxdt_ORBC     = kops_eff * RBC_CIRC - phago_RBC - kdeg_RBC * ORBC;
dxdt_RBC_CIRC = kinput_RBC - kops_eff * RBC_CIRC - kdeg_RBC * RBC_CIRC;

// ===== EPO / BM erythropoiesis =====
double Hgb_calc    = RBC_CIRC * Hgb_per_RBC;  // Hgb g/dL
double EPO_fold    = 1.0 + (kEPO_stim - 1.0) * EPO_EC50 / (EPO_EC50 + Hgb_calc);
dxdt_BM_ERY = kprol_BM_Ery * EPO_fold * BM_Ery0 - kdeg_BM_Ery * BM_ERY;
// Reticulocytes produced from BM
double Retic_prod = kdeg_BM_Ery * BM_ERY * 0.10;  // fraction exported as retics
dxdt_RETIC = Retic_prod - kmat_Retic * RETIC;
// Mature RBC influx from BM (reticulocytes mature → RBC)
double influx_RBC = kmat_Retic * RETIC * RBC0 / Retic0;
dxdt_RBC_CIRC = influx_RBC - kops_eff * RBC_CIRC - kdeg_RBC * RBC_CIRC;  // overrides above

// ===== Platelet destruction =====
double kops_Plt_eff = kops_Plt * AB_PLT;
double phago_Plt    = kphago_Plt * FcgR_Plt * OPLT;
dxdt_OPLT     = kops_Plt_eff * PLT_CIRC - phago_Plt - kdeg_Plt * OPLT;
dxdt_PLT_CIRC = kinput_Plt - kops_Plt_eff * PLT_CIRC - kdeg_Plt * PLT_CIRC;

// ===== Megakaryocyte pool =====
// Plt count feeds back on TPO -> MK (negative feedback)
double TPO_fold = 1.0 + (3.0 - 1.0) * Plt0 / (Plt0 + PLT_CIRC);
dxdt_MK = kprol_MK * TPO_fold * E_Elt_MK * MK0 - kdeg_MK * MK;
// Platelet production from MK
double Plt_prod = krel_Plt * MK;
dxdt_PLT_CIRC = Plt_prod - kops_Plt_eff * PLT_CIRC - kdeg_Plt * PLT_CIRC;  // overrides above

$TABLE
capture Hgb    = RBC_CIRC * Hgb_per_RBC;         // g/dL
capture Plt    = PLT_CIRC;                         // x10^9/L
capture Retic_pct = RETIC;                         // %
capture C_Prednisolone = CENT_PRED / Vc_Pred;     // mg/L
capture C_MPA_plasma   = CENT_MPA / Vc_MPA;       // mg/L
capture C_Rtx_plasma   = CENT_RTX / Vc_Rtx;       // mg/L
capture C_Elt_plasma   = CENT_ELT / Vc_Elt;       // mg/L
capture C_Siro_blood   = CENT_SIRO / Vc_Siro;     // ng/mL
capture C_IVIg_plasma  = CENT_IVIG / Vc_IVIg;    // mg/L
capture Treg_cells     = TREG;                     // cells/uL
capture Ab_RBC_level   = AB_RBC;                   // mg/L
capture Ab_Plt_level   = AB_PLT;                   // mg/L
capture B_RBC_count    = B_RBC;                    // cells/uL
capture B_Plt_count    = B_PLT;                    // cells/uL
capture C3b_level      = C3B;                      // AU
'

## ============================================================
## 2.  Compile Model
## ============================================================
es_mod <- mcode("Evans_Syndrome_QSP", ES_code)

## ============================================================
## 3.  Helper – Make Event Table
## ============================================================
make_events <- function(scenario) {
  ev <- ev(time=0, amt=0, cmt=1)  # null event
  if (scenario == "no_treatment") return(ev)

  if (scenario == "prednisolone_mono") {
    return(ev(time=0, amt=0, cmt=1))  # handled via PARAM
  }
  if (scenario == "rituximab") {
    # 4 weekly doses of rituximab (375 mg/m² × 1.7m² ≈ 640mg)
    return(ev(time=c(1, 8, 15, 22), amt=640, cmt="CENT_RTX", ii=0, addl=0))
  }
  if (scenario == "ivig") {
    # IVIg 1g/kg = 70g = 70,000mg, Vc=3.5L => ~20,000 mg in central?
    return(ev(time=c(1, 2), amt=35000, cmt="CENT_IVIG", ii=0, addl=0))
  }
  if (scenario == "combo_rtx_ivig") {
    ev_rtx  <- ev(time=c(1, 8, 15, 22), amt=640, cmt="CENT_RTX")
    ev_ivig <- ev(time=c(1, 2), amt=35000, cmt="CENT_IVIG")
    return(ev_rtx + ev_ivig)
  }
  return(ev(time=0, amt=0, cmt=1))
}

## ============================================================
## 4.  Simulation Parameters for 5 Scenarios
## ============================================================
params_list <- list(

  # Scenario 1: No treatment (natural disease course)
  no_treatment = list(
    PRED_dose=0, IVIG_dose=0, RTX_dose=0, MMF_dose=0,
    ELT_dose=0, SIRO_dose=0, SPLEN=0,
    label="No Treatment (Natural Course)"
  ),

  # Scenario 2: Corticosteroid monotherapy (first-line)
  prednisolone = list(
    PRED_dose=1.5, IVIG_dose=0, RTX_dose=0, MMF_dose=0,
    ELT_dose=0, SIRO_dose=0, SPLEN=0,
    label="Prednisone 1.5 mg/kg/day (1st-line)"
  ),

  # Scenario 3: Rituximab + prednisolone (2nd-line)
  pred_rtx = list(
    PRED_dose=0.5, IVIG_dose=0, RTX_dose=1, MMF_dose=0,
    ELT_dose=0, SIRO_dose=0, SPLEN=0,
    label="Prednisone 0.5 mg/kg + Rituximab"
  ),

  # Scenario 4: MMF + low-dose pred (3rd-line)
  pred_mmf = list(
    PRED_dose=0.25, IVIG_dose=0, RTX_dose=0, MMF_dose=2.0,
    ELT_dose=0, SIRO_dose=0, SPLEN=0,
    label="Prednisone 0.25mg/kg + MMF 2g/day"
  ),

  # Scenario 5: Sirolimus + eltrombopag (refractory)
  siro_elt = list(
    PRED_dose=0.1, IVIG_dose=0, RTX_dose=0, MMF_dose=0,
    ELT_dose=50, SIRO_dose=4, SPLEN=0,
    label="Sirolimus 4mg + Eltrombopag 50mg (refractory)"
  ),

  # Scenario 6: Splenectomy (surgical option)
  splenectomy = list(
    PRED_dose=0.1, IVIG_dose=0, RTX_dose=0, MMF_dose=0,
    ELT_dose=0, SIRO_dose=0, SPLEN=1,
    label="Splenectomy + Low-dose Prednisone"
  )
)

## ============================================================
## 5.  Run Simulations
## ============================================================
sim_time <- seq(0, 365, by=1)  # 1 year simulation (days)

run_scenario <- function(pname, pvals, sim_t = sim_time) {
  p_update <- pvals[names(pvals) != "label"]
  ev_obj <- make_events(pname)

  out <- es_mod %>%
    param(p_update) %>%
    ev(ev_obj) %>%
    mrgsim(end=max(sim_t), delta=1, obsonly=TRUE) %>%
    as.data.frame() %>%
    mutate(scenario = pvals$label)
  out
}

results <- bind_rows(lapply(names(params_list), function(nm) {
  run_scenario(nm, params_list[[nm]])
}))

## ============================================================
## 6.  Quick Diagnostic Plots
## ============================================================
scenario_colors <- c(
  "No Treatment (Natural Course)"          = "#E53935",
  "Prednisone 1.5 mg/kg/day (1st-line)"   = "#FB8C00",
  "Prednisone 0.5 mg/kg + Rituximab"      = "#1E88E5",
  "Prednisone 0.25mg/kg + MMF 2g/day"    = "#43A047",
  "Sirolimus 4mg + Eltrombopag 50mg (refractory)" = "#8E24AA",
  "Splenectomy + Low-dose Prednisone"     = "#00897B"
)

p_hgb <- ggplot(results, aes(time, Hgb, color=scenario)) +
  geom_line(size=1.1) +
  geom_hline(yintercept=12, linetype="dashed", color="gray40") +
  labs(title="Hemoglobin over Time", x="Day", y="Hgb (g/dL)",
       color="Scenario") +
  scale_color_manual(values=scenario_colors) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", legend.text=element_text(size=8))

p_plt <- ggplot(results, aes(time, Plt, color=scenario)) +
  geom_line(size=1.1) +
  geom_hline(yintercept=100, linetype="dashed", color="gray40") +
  labs(title="Platelet Count over Time", x="Day", y="Plt (×10⁹/L)",
       color="Scenario") +
  scale_color_manual(values=scenario_colors) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", legend.text=element_text(size=8))

p_ab_rbc <- ggplot(results, aes(time, Ab_RBC_level, color=scenario)) +
  geom_line(size=1.1) +
  labs(title="Anti-RBC IgG Autoantibody", x="Day", y="Ab-RBC (mg/L)",
       color="Scenario") +
  scale_color_manual(values=scenario_colors) +
  theme_bw(base_size=11) +
  theme(legend.position="none")

p_treg <- ggplot(results, aes(time, Treg_cells, color=scenario)) +
  geom_line(size=1.1) +
  labs(title="Regulatory T cells (Treg)", x="Day", y="Treg (cells/µL)",
       color="Scenario") +
  scale_color_manual(values=scenario_colors) +
  theme_bw(base_size=11) +
  theme(legend.position="none")

combined_plot <- (p_hgb + p_plt) / (p_ab_rbc + p_treg)
print(combined_plot)

## ============================================================
## 7.  Response Summary at Day 180
## ============================================================
summary_d180 <- results %>%
  filter(time == 180) %>%
  select(scenario, Hgb, Plt, Retic_pct, Ab_RBC_level, Ab_Plt_level,
         Treg_cells, B_RBC_count, B_Plt_count) %>%
  mutate(
    AIHA_Response = case_when(
      Hgb >= 12 ~ "Complete",
      Hgb >= 10 ~ "Partial",
      TRUE       ~ "No Response"
    ),
    ITP_Response = case_when(
      Plt >= 100 ~ "Complete",
      Plt >= 30  ~ "Partial",
      TRUE        ~ "No Response"
    )
  ) %>%
  arrange(desc(Hgb))

cat("\n===== Evans Syndrome QSP: Treatment Response at Day 180 =====\n")
print(summary_d180, n=20)

## ============================================================
## 8.  PK Profile Plots (Rituximab + Prednisolone scenarios)
## ============================================================
pk_rtx <- results %>%
  filter(scenario == "Prednisone 0.5 mg/kg + Rituximab", time <= 60)

p_pk_rtx <- ggplot(pk_rtx, aes(time, C_Rtx_plasma)) +
  geom_line(color="#1E88E5", size=1.2) +
  labs(title="Rituximab PK (Scenario 3)", x="Day", y="Rtx (mg/L)") +
  theme_bw(base_size=11)

pk_pred <- results %>%
  filter(scenario == "Prednisone 1.5 mg/kg/day (1st-line)", time <= 30)

p_pk_pred <- ggplot(pk_pred, aes(time, C_Prednisolone)) +
  geom_line(color="#FB8C00", size=1.2) +
  labs(title="Prednisolone PK (Scenario 2)", x="Day", y="Prednisolone (mg/L)") +
  theme_bw(base_size=11)

print(p_pk_rtx + p_pk_pred)
