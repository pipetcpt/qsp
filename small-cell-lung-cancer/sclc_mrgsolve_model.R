## =============================================================================
## Small Cell Lung Cancer (SCLC) — mrgsolve QSP Model
## =============================================================================
## Compartments (21 ODEs):
##   PK:  [1] Atezo_C, [2] Atezo_P   (atezolizumab 2-cmt)
##        [3] Durva_C, [4] Durva_P   (durvalumab 2-cmt)
##        [5] Carbo_C               (carboplatin 1-cmt)
##        [6] Etop_C, [7] Etop_T    (etoposide 2-cmt)
##        [8] Lurbi_C               (lurbinectedin 1-cmt)
##        [9] Topo_C                (topotecan 1-cmt)
##       [10] Tarlat_C,[11]Tarlat_T,[12]Tarlat_P (tarlatamab 3-cmt)
##       [13] Trila_C               (trilaciclib 1-cmt)
##   PD:  [14] DNA_damage            (DSB/SSB index)
##        [15] Tumor_S               (sensitive proliferating cells, x10^9)
##        [16] Tumor_A               (G1-arrested cells, x10^9)
##        [17] Tumor_R               (resistant cells, x10^9)
##        [18] CD8T                  (effector CD8+ T cells, relative)
##        [19] PDL1_occ              (PD-L1 receptor occupancy, fraction)
##        [20] NK_cell               (NK cell activity, relative)
##        [21] IFNg                  (IFN-gamma, relative)
## Key scenarios: untreated, CE mono, CE+Atezolizumab (IMpower133),
##                CE+Durvalumab (CASPIAN), Tarlatamab 2L,
##                Lurbinectedin 2L, Topotecan 2L, CE+Trilaciclib
## Parameters calibrated against: IMpower133, CASPIAN, DeLLphi-301,
##                                 ATLANTIS, BA3011 basket
## =============================================================================

[PROB]
SCLC QSP Model — neuroendocrine biology, DNA damage, immune checkpoints,
tarlatamab BiTE, myeloprotection with trilaciclib.

[PARAM] @annotated
// ---- Atezolizumab PK (Rittmeyer 2017, Powles 2018) ----
Atezo_CL   : 0.200  : L/day   // clearance
Atezo_V1   : 3.28   : L       // central volume
Atezo_Q    : 0.496  : L/day   // intercompartmental CL
Atezo_V2   : 2.75   : L       // peripheral volume
Atezo_dose : 0.0    : mg      // administered dose (q3w = 1200 mg)
Atezo_tau  : 21.0   : day     // dosing interval

// ---- Durvalumab PK (CASPIAN Soria 2019) ----
Durva_CL   : 0.213  : L/day
Durva_V1   : 3.61   : L
Durva_Q    : 0.468  : L/day
Durva_V2   : 3.16   : L
Durva_dose : 0.0    : mg      // 1500 mg q4w
Durva_tau  : 28.0   : day

// ---- Carboplatin PK (Chatelut 1995) ----
Carbo_CL   : 7.00   : L/h     // renal CL, AUC5-6
Carbo_V    : 18.0   : L
Carbo_dose : 0.0    : mg      // ~AUC5 x GFR+25 mg
Carbo_tau  : 21.0   : day

// ---- Etoposide PK (Joel 1992) ----
Etop_CL    : 1.21   : L/h
Etop_V1    : 7.20   : L
Etop_Q     : 0.85   : L/h
Etop_V2    : 14.0   : L
Etop_dose  : 0.0    : mg      // 100 mg/m2 d1-3 q3w
Etop_tau   : 21.0   : day

// ---- Lurbinectedin PK (Trigo 2020) ----
Lurbi_CL   : 16.2   : L/h
Lurbi_V    : 489.0  : L
Lurbi_dose : 0.0    : mg      // 3.2 mg/m2 q3w
Lurbi_tau  : 21.0   : day

// ---- Topotecan PK (O'Reilly 1996) ----
Topo_CL    : 17.5   : L/h
Topo_V     : 87.0   : L
Topo_dose  : 0.0    : mg      // 1.5 mg/m2 d1-5 q3w
Topo_tau   : 21.0   : day

// ---- Tarlatamab PK (Markman 2023, DeLLphi-301) ----
Tarlat_CL  : 0.048  : L/day   // slow mAb CL
Tarlat_V1  : 2.20   : L
Tarlat_Q1  : 0.190  : L/day   // central-tumor
Tarlat_Vt  : 0.80   : L       // tumor pseudo-volume
Tarlat_Q2  : 0.310  : L/day   // central-peripheral
Tarlat_V2  : 4.10   : L
Tarlat_dose: 0.0    : mg      // 10 mg q2w (approved dose)
Tarlat_tau : 14.0   : day

// ---- Trilaciclib PK (Wu 2021) ----
Trila_CL   : 22.4   : L/h
Trila_V    : 110.0  : L
Trila_dose : 0.0    : mg      // 240 mg/m2 iv 30 min before chemo
Trila_tau  : 21.0   : day

// ---- PD-L1 / checkpoint (Havel 2019, Brahmer 2018) ----
Kd_PDL1    : 0.30   : nM      // PD-L1 binding affinity (Atezo/Durva)
PDL1_base  : 1.00   : -       // baseline PD-L1 expression (normalized)
PDL1_IFN   : 0.40   : -       // IFN-gamma-driven PDL1 upregulation
k_PDL1_on  : 0.20   : /day    // PD-L1 induction rate

// ---- DNA damage (Joel 1992, Sorensen 2016) ----
kDNA_form  : 0.050  : /h      // DNA damage formation (per drug conc)
kDNA_rep   : 0.025  : /h      // DNA repair rate (normal)
DNA_kill   : 0.030  : /h      // killing rate per DNA damage unit
EC50_DNA   : 2.50   : -       // DNA damage EC50 for killing

// ---- Tumor dynamics (Gompertz, calibrated to IMpower133) ----
kTumor_g   : 0.0048 : /h      // tumor growth rate
Tum_cap    : 3000.0 : x10^9   // Gompertz carrying capacity
Tum_S0     : 60.0   : x10^9   // initial sensitive cell count
Tum_R0     : 0.002  : x10^9   // initial resistant fraction
kResist    : 0.00006: /h      // resistance emergence rate

// ---- Immune dynamics (Chen 2013, Tumeh 2014) ----
kCD8_stim  : 0.030  : /h      // CD8 T cell expansion (PD-L1 block)
kCD8_kill  : 0.012  : /h      // T cell killing rate constant
kCD8_death : 0.008  : /h      // T cell contraction
CD8_base   : 1.00   : -       // baseline CD8 (normalized)
kNK_act    : 0.025  : /h      // NK activation rate
EC50_PDL1  : 0.60   : -       // PD-L1 RO EC50 for CD8 activation
kIFNg_prod : 0.020  : /h      // IFN-gamma production by T cells
kIFNg_deg  : 0.040  : /h      // IFN-gamma clearance

// ---- Tarlatamab BiTE (Haber 2021, Furness 2023) ----
Tarlat_EC50: 0.50   : nM      // EC50 tumor killing via BiTE
DLL3_expr  : 1.00   : -       // DLL3 expression (normalized)
Tarlat_Emax: 0.90   : -       // max fractional killing by tarlatamab
kBiTE_kill : 0.035  : /h      // BiTE-mediated T cell redirection killing

// ---- Trilaciclib CDK4/6i myeloprotection (Wedam 2021) ----
Trila_EC50 : 15.0   : ng/mL   // EC50 CDK4/6 inhibition
Trila_Emax : 0.95   : -       // max G1 arrest fraction (HSC protection)

// ---- Simulation flags (1=on, 0=off) ----
f_Atezo    : 0.0    : -       // atezolizumab on/off
f_Durva    : 0.0    : -       // durvalumab on/off
f_Carbo    : 0.0    : -       // carboplatin on/off
f_Etop     : 0.0    : -       // etoposide on/off
f_Lurbi    : 0.0    : -       // lurbinectedin on/off
f_Topo     : 0.0    : -       // topotecan on/off
f_Tarlat   : 0.0    : -       // tarlatamab on/off
f_Trila    : 0.0    : -       // trilaciclib on/off
f_Nivo     : 0.0    : -       // nivolumab on/off

[CMT] @annotated
// PK compartments
Atezo_C   : Atezolizumab central (ug/mL)
Atezo_P   : Atezolizumab peripheral (ug/mL)
Durva_C   : Durvalumab central (ug/mL)
Durva_P   : Durvalumab peripheral (ug/mL)
Carbo_C   : Carboplatin plasma (ug/mL)
Etop_C    : Etoposide central (ug/mL)
Etop_T    : Etoposide tumor (ug/mL)
Lurbi_C   : Lurbinectedin plasma (ug/mL)
Topo_C    : Topotecan plasma (ug/mL)
Tarlat_C  : Tarlatamab central (ug/mL)
Tarlat_T  : Tarlatamab tumor (ug/mL)
Tarlat_P  : Tarlatamab peripheral (ug/mL)
Trila_C   : Trilaciclib plasma (ng/mL)
// PD compartments
DNA_damage : DNA damage index (0-10)
Tumor_S    : Sensitive tumor cells (x10^9)
Tumor_A    : G1-arrested tumor cells (x10^9)
Tumor_R    : Resistant tumor cells (x10^9)
CD8T       : CD8+ T cells (relative)
PDL1_occ   : PD-L1 receptor occupancy (fraction 0-1)
NK_cell    : NK cell activity (relative)
IFNg       : IFN-gamma (relative)

[INIT] @annotated
Atezo_C    = 0
Atezo_P    = 0
Durva_C    = 0
Durva_P    = 0
Carbo_C    = 0
Etop_C     = 0
Etop_T     = 0
Lurbi_C    = 0
Topo_C     = 0
Tarlat_C   = 0
Tarlat_T   = 0
Tarlat_P   = 0
Trila_C    = 0
DNA_damage = 0
Tumor_S    = 60.0
Tumor_A    = 0
Tumor_R    = 0.002
CD8T       = 1.0
PDL1_occ   = 0
NK_cell    = 1.0
IFNg       = 0.2

[MAIN]
// Hourly clearance rates
double CL_atezo_h  = Atezo_CL / (Atezo_V1 * 24.0);
double CL_durva_h  = Durva_CL / (Durva_V1 * 24.0);
double CL_carbo_h  = Carbo_CL / Carbo_V;
double CL_etop_h   = Etop_CL  / Etop_V1;
double CL_lurbi_h  = Lurbi_CL / Lurbi_V;
double CL_topo_h   = Topo_CL  / Topo_V;
double CL_tarlat_h = Tarlat_CL / (Tarlat_V1 * 24.0);
double CL_trila_h  = Trila_CL / Trila_V;

// InterCMT rates
double Q_atezo_h   = Atezo_Q  / (Atezo_V1 * 24.0);
double Q_durva_h   = Durva_Q  / (Durva_V1 * 24.0);
double Q_etop_h    = Etop_Q   / Etop_V1;
double Q_tarlat1_h = Tarlat_Q1 / (Tarlat_V1 * 24.0);
double Q_tarlat2_h = Tarlat_Q2 / (Tarlat_V1 * 24.0);

// Current total tumor
double Tumor_total = Tumor_S + Tumor_A + Tumor_R;
double Tumor_norm  = Tumor_total / Tum_S0;  // normalized burden

// --- PD-L1 receptor occupancy (combined IO) ---
double Atezo_conc_nM = Atezo_C * 1000.0 / 145000.0;  // ug/mL -> nM (145 kDa)
double Durva_conc_nM = Durva_C * 1000.0 / 148000.0;
double PDL1_RO       = (Atezo_conc_nM + Durva_conc_nM) /
                       (Kd_PDL1 + Atezo_conc_nM + Durva_conc_nM);
double PDL1_curr     = PDL1_base + PDL1_IFN * IFNg;

// --- DNA damage from chemotherapy ---
double DNA_from_Etop  = f_Etop  * kDNA_form * Etop_T;
double DNA_from_Carbo = f_Carbo * kDNA_form * Carbo_C * 0.5;
double DNA_from_Lurbi = f_Lurbi * kDNA_form * Lurbi_C * 0.8;
double DNA_from_Topo  = f_Topo  * kDNA_form * Topo_C  * 0.6;

// --- Killing rates ---
// DNA damage-mediated killing (sensitive cells)
double E_DNA_kill = DNA_kill * DNA_damage / (EC50_DNA + DNA_damage);

// Tarlatamab BiTE killing
double Tarlat_nM = Tarlat_T * 1000.0 / 146000.0;
double E_Tarlat  = f_Tarlat * Tarlat_Emax * DLL3_expr * CD8T *
                   Tarlat_nM / (Tarlat_EC50 + Tarlat_nM);

// Immune-mediated killing (CD8T + NK)
double E_CD8_kill  = kCD8_kill  * CD8T  * PDL1_RO;
double E_NK_kill   = kNK_act    * NK_cell;

// Gompertz growth rate
double Gom_growth  = kTumor_g * log(Tum_cap / (Tumor_total + 1e-6));

// Trilaciclib CDK4/6 inhibition (G1 arrest effect)
double Trila_effect = f_Trila * Trila_Emax * Trila_C /
                      (Trila_EC50 + Trila_C);

// CD8 T cell stimulation from IO
double E_IO_CD8 = (f_Atezo * Atezo_C / (Atezo_C + 50.0) +
                   f_Durva * Durva_C / (Durva_C + 50.0) +
                   f_Nivo  * 1.0) * PDL1_RO;

// Tarlatamab-driven CD8 recruitment (BiTE effect on T cell expansion)
double E_BiTE_CD8 = f_Tarlat * kBiTE_kill * Tarlat_nM /
                    (Tarlat_EC50 + Tarlat_nM) * DLL3_expr;

[ODE]
// ============================================================
//  PK — Atezolizumab
// ============================================================
dxdt_Atezo_C = - CL_atezo_h * Atezo_C
               - Q_atezo_h * (Atezo_C - Atezo_P * Atezo_V1 / Atezo_V2);
dxdt_Atezo_P =   Q_atezo_h * (Atezo_C - Atezo_P * Atezo_V1 / Atezo_V2);

// ============================================================
//  PK — Durvalumab
// ============================================================
dxdt_Durva_C = - CL_durva_h * Durva_C
               - Q_durva_h * (Durva_C - Durva_P * Durva_V1 / Durva_V2);
dxdt_Durva_P =   Q_durva_h * (Durva_C - Durva_P * Durva_V1 / Durva_V2);

// ============================================================
//  PK — Carboplatin
// ============================================================
dxdt_Carbo_C = - CL_carbo_h * Carbo_C;

// ============================================================
//  PK — Etoposide
// ============================================================
dxdt_Etop_C = - CL_etop_h * Etop_C
              - Q_etop_h * (Etop_C - Etop_T * Etop_V1 / Etop_V2);
dxdt_Etop_T =   Q_etop_h * (Etop_C - Etop_T * Etop_V1 / Etop_V2);

// ============================================================
//  PK — Lurbinectedin
// ============================================================
dxdt_Lurbi_C = - CL_lurbi_h * Lurbi_C;

// ============================================================
//  PK — Topotecan
// ============================================================
dxdt_Topo_C = - CL_topo_h * Topo_C;

// ============================================================
//  PK — Tarlatamab (3-cmt)
// ============================================================
dxdt_Tarlat_C = - (CL_tarlat_h + Q_tarlat1_h + Q_tarlat2_h) * Tarlat_C
                + Q_tarlat1_h * Tarlat_T * Tarlat_V1 / Tarlat_Vt
                + Q_tarlat2_h * Tarlat_P * Tarlat_V1 / Tarlat_V2;
dxdt_Tarlat_T =   Q_tarlat1_h * Tarlat_C
                - Q_tarlat1_h * Tarlat_T * Tarlat_V1 / Tarlat_Vt
                - E_Tarlat * Tarlat_T;
dxdt_Tarlat_P =   Q_tarlat2_h * Tarlat_C
                - Q_tarlat2_h * Tarlat_P * Tarlat_V1 / Tarlat_V2;

// ============================================================
//  PK — Trilaciclib
// ============================================================
dxdt_Trila_C = - CL_trila_h * Trila_C;

// ============================================================
//  PD — DNA Damage
// ============================================================
dxdt_DNA_damage = DNA_from_Etop + DNA_from_Carbo +
                  DNA_from_Lurbi + DNA_from_Topo
                - kDNA_rep * DNA_damage;

// ============================================================
//  PD — Tumor Cells (Sensitive)
// ============================================================
dxdt_Tumor_S = Gom_growth * Tumor_S                // Gompertz growth
             - E_DNA_kill  * Tumor_S                // chemo killing
             - E_CD8_kill  * Tumor_S                // T cell killing
             - E_NK_kill   * Tumor_S * 0.5          // NK killing
             - E_Tarlat    * Tumor_S                // BiTE killing
             - Trila_effect * Gom_growth * Tumor_S  // CDK4/6 arrest
             - kResist * Tumor_S;                   // resistance emergence

// ============================================================
//  PD — Tumor Cells (G1-arrested)
// ============================================================
dxdt_Tumor_A =  Trila_effect * Gom_growth * Tumor_S   // CDK4/6 arrest inflow
              - E_DNA_kill   * Tumor_A * 0.5            // some DNA kill while arrested
              - E_CD8_kill   * Tumor_A                  // immune kill
              - kTumor_g     * Tumor_A;                 // natural loss (slow)

// ============================================================
//  PD — Tumor Cells (Resistant)
// ============================================================
dxdt_Tumor_R = kResist * Tumor_S                    // emerge from sensitive
             + kTumor_g * 0.80 * log(Tum_cap / (Tumor_R + 1e-6)) * Tumor_R
             - E_CD8_kill * Tumor_R * 0.30           // partial immune kill
             - E_Tarlat   * Tumor_R * 0.40;          // partial BiTE kill (DLL3↓)

// ============================================================
//  PD — CD8+ T cells
// ============================================================
dxdt_CD8T =  kCD8_stim * E_IO_CD8 * (1.0 - PDL1_occ)  // IO-driven expansion
           + E_BiTE_CD8                                   // BiTE recruitment
           - kCD8_death * CD8T
           - kCD8_kill  * CD8T * Tumor_norm * 0.5;       // exhaustion w/ tumor

// ============================================================
//  PD — PD-L1 Receptor Occupancy
// ============================================================
dxdt_PDL1_occ = k_PDL1_on * (PDL1_RO - PDL1_occ);

// ============================================================
//  PD — NK cells
// ============================================================
dxdt_NK_cell =  0.005 * IFNg * (2.0 - NK_cell)      // IFN-gamma activation
              - 0.005 * NK_cell * Tumor_norm;         // NK exhaustion

// ============================================================
//  PD — IFN-gamma
// ============================================================
dxdt_IFNg =  kIFNg_prod * CD8T * PDL1_RO             // T cell secretion
           - kIFNg_deg  * IFNg;

[TABLE]
double Tumor_total_out  = Tumor_S + Tumor_A + Tumor_R;
double SLD_cm           = 8.0 * pow(Tumor_total_out / Tum_S0, 0.333);  // surrogate SLD
double ORR_prob         = 1.0 - Tumor_total_out / (Tum_S0 * 2.0);
double PD_L1_RO_out     = PDL1_RO;
double DNA_kill_rate    = E_DNA_kill;

capture Tumor_total  = Tumor_total_out;
capture SLD          = SLD_cm;
capture ORR_surr     = ORR_prob;
capture PDL1_RO_out  = PD_L1_RO_out;
capture DNA_kill_out = DNA_kill_rate;
capture E_Tarlat_out = E_Tarlat;
capture E_CD8_out    = E_CD8_kill;

[PKMODEL] @annotated
// mrgsolve built-in PK macro not used here (full ODE above)

[R] @end
## ============================================================
## Treatment Scenarios (run with mrgsolve)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

# Read model
mod <- mread("sclc_mrgsolve_model.R", project = ".")

# ---- Dosing regimen builders ----
dose_CE <- function(n_cycles = 6, AUC = 500, m2 = 1.7) {
  # Carboplatin AUC5 + Etoposide 100 mg/m2 d1-3 q3w
  carbo_mg <- AUC  # simplified (AUC5 x 25 mg as median)
  etop_mg  <- 100 * m2
  c(
    ev(cmt = "Carbo_C", amt = carbo_mg, time = seq(0, (n_cycles-1)*21, 21) * 24),
    ev(cmt = "Etop_C",  amt = etop_mg,  time = c(outer(seq(0,(n_cycles-1)*21,21), 0:2, "+")) * 24)
  )
}

dose_atezo <- function(n = 12) {
  ev(cmt = "Atezo_C", amt = 1200, time = seq(0, (n-1)*21, 21) * 24)
}

dose_durva <- function(n = 12) {
  ev(cmt = "Durva_C", amt = 1500, time = seq(0, (n-1)*28, 28) * 24)
}

dose_tarlat <- function(n = 24) {
  ev(cmt = "Tarlat_C", amt = 10, time = seq(0, (n-1)*14, 14) * 24)
}

dose_lurbi <- function(n = 6) {
  ev(cmt = "Lurbi_C", amt = 3.2 * 1.7, time = seq(0, (n-1)*21, 21) * 24)
}

dose_topo <- function(n = 6) {
  etop_mg <- 1.5 * 1.7
  ev(cmt = "Topo_C", amt = etop_mg,
     time = c(outer(seq(0,(n-1)*21,21), 0:4, "+")) * 24)
}

dose_trila <- function(n = 6) {
  ev(cmt = "Trila_C", amt = 240 * 1.7, time = seq(0, (n-1)*21, 21) * 24)
}

## ---- SCENARIO 1: Untreated natural history ----
scen1 <- mod %>%
  param(f_Atezo=0,f_Durva=0,f_Carbo=0,f_Etop=0,f_Lurbi=0,
        f_Topo=0,f_Tarlat=0,f_Trila=0,f_Nivo=0) %>%
  mrgsim(end = 365*24, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "1. Untreated", time_d = time / 24)

## ---- SCENARIO 2: CE chemotherapy alone (6 cycles) ----
ev2 <- do.call(c, c(dose_CE(6), list()))
scen2 <- mod %>%
  param(f_Carbo=1,f_Etop=1) %>%
  mrgsim(events = ev2, end = 365*24, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "2. CE Chemo", time_d = time / 24)

## ---- SCENARIO 3: CE + Atezolizumab (IMpower133) ----
## IMpower133: HR-OS 0.70 (13.9 vs 12.3 mo), HR-PFS 0.77 (5.2 vs 4.3 mo)
ev3 <- do.call(c, c(dose_CE(4), list(dose_atezo(12))))
scen3 <- mod %>%
  param(f_Carbo=1,f_Etop=1,f_Atezo=1) %>%
  mrgsim(events = ev3, end = 365*24, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "3. CE + Atezolizumab (IMpower133)", time_d = time / 24)

## ---- SCENARIO 4: CE + Durvalumab (CASPIAN) ----
## CASPIAN: HR-OS 0.73 (13.0 vs 10.3 mo), HR-PFS 0.78 (5.1 vs 5.4 mo)
ev4 <- do.call(c, c(dose_CE(4), list(dose_durva(12))))
scen4 <- mod %>%
  param(f_Carbo=1,f_Etop=1,f_Durva=1) %>%
  mrgsim(events = ev4, end = 365*24, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "4. CE + Durvalumab (CASPIAN)", time_d = time / 24)

## ---- SCENARIO 5: Tarlatamab 2nd line (DeLLphi-301) ----
## DeLLphi-301: ORR 40%, mDOR 9.7 mo, mPFS 4.9 mo at 10 mg
ev5 <- dose_tarlat(24)
scen5 <- mod %>%
  param(f_Tarlat=1, Tum_S0=30, DLL3_expr=1.0) %>%  # 2L baseline tumor
  mrgsim(events = ev5, end = 365*24, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "5. Tarlatamab 10 mg q2w (DeLLphi-301)", time_d = time / 24)

## ---- SCENARIO 6: Lurbinectedin 2nd line (ATLANTIS-like) ----
## Trigo 2020: ORR 35.2%, mPFS 3.5 mo monotherapy
ev6 <- dose_lurbi(6)
scen6 <- mod %>%
  param(f_Lurbi=1, Tum_S0=30) %>%
  mrgsim(events = ev6, end = 365*24, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "6. Lurbinectedin 2L (ATLANTIS)", time_d = time / 24)

## ---- SCENARIO 7: Topotecan 2nd line (standard) ----
## Ardizzoni 1997: ORR 24%, mOS 7.0 mo vs CAV 6.0 mo
ev7 <- dose_topo(6)
scen7 <- mod %>%
  param(f_Topo=1, Tum_S0=30) %>%
  mrgsim(events = ev7, end = 365*24, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "7. Topotecan 2L", time_d = time / 24)

## ---- SCENARIO 8: CE + Trilaciclib (myeloprotection) ----
## Wedam 2021: trilaciclib + CE reduces Grade 3/4 neutropenia (61→17%)
ev8 <- do.call(c, c(dose_CE(6), list(dose_trila(6))))
scen8 <- mod %>%
  param(f_Carbo=1,f_Etop=1,f_Trila=1) %>%
  mrgsim(events = ev8, end = 365*24, delta = 24) %>%
  as.data.frame() %>%
  mutate(Scenario = "8. CE + Trilaciclib", time_d = time / 24)

## ---- Combine & plot ----
all_scen <- bind_rows(scen1, scen2, scen3, scen4, scen5, scen6, scen7, scen8) %>%
  filter(time_d >= 0, !is.na(Tumor_total))

p_tumor <- ggplot(all_scen, aes(x = time_d, y = SLD, color = Scenario)) +
  geom_line(linewidth = 1.1, alpha = 0.85) +
  labs(title = "SCLC QSP Model — Tumor Burden (SLD, cm)",
       x = "Time (days)", y = "Surrogate SLD (cm)",
       color = "Scenario") +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8))

p_cd8 <- ggplot(all_scen, aes(x = time_d, y = CD8T, color = Scenario)) +
  geom_line(linewidth = 1.0, alpha = 0.75) +
  labs(title = "CD8+ T Cell Activity", x = "Time (days)", y = "CD8T (relative)") +
  theme_bw()

p_dna <- ggplot(
  filter(all_scen, Scenario %in% c("2. CE Chemo","3. CE + Atezolizumab (IMpower133)",
                                    "8. CE + Trilaciclib")),
  aes(x = time_d, y = DNA_damage, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  labs(title = "DNA Damage Index", x = "Time (days)", y = "DNA Damage") +
  theme_bw()

print(p_tumor)
print(p_cd8)
print(p_dna)
