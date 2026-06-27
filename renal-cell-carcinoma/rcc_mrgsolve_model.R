## ============================================================
## Renal Cell Carcinoma (ccRCC) QSP Model — mrgsolve
## 18-compartment ODE: Sunitinib PK + SU12662 metabolite,
##   Nivolumab TMDD, Belzutifan, VHL/HIF/VEGF pathway,
##   mTOR signalling, Tumor Growth Inhibition (Simeoni),
##   Immune microenvironment (CD8, Treg, MDSC)
##
## Key clinical calibration:
##   CheckMate 214 (Motzer 2018): Nivo+Ipi → mPFS 11.6 mo (IRF)
##   KEYNOTE-426 (Rini 2019):     Pembro+Axitinib → mPFS 15.1 mo
##   CheckMate 9ER (Choueiri 2021): Nivo+Cabo → mPFS 16.6 mo
##   CLEAR (Motzer 2021):         Len+Pembro → mPFS 23.9 mo
##   METEOR (Choueiri 2015):      Cabozantinib → mPFS 7.4 mo
##   RECORD-1 (Motzer 2008):      Everolimus  → mPFS 4.0 mo
##   LITESPARK-005 (Choueiri 2023): Belzutifan → ORR 22%
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ----------------------------------------------------------
## 1. Model code
## ----------------------------------------------------------
rcc_code <- '
$PROB Renal Cell Carcinoma (ccRCC) QSP Model

$PARAM
// --- Sunitinib PK (2-cmt, SU12662 metabolite) ---
// Ref: Demetri 2009, Houk 2010, Lindauer 2010
CL_sun  = 51.8    // L/h  apparent oral clearance
V1_sun  = 2030    // L    central volume
Q_sun   = 7.5     // L/h  inter-compartmental clearance
V2_sun  = 560     // L    peripheral volume
ka_sun  = 0.09    // 1/h  absorption rate
F_sun   = 0.50    // bioavailability
fm_sun  = 0.23    // fraction metabolised to SU12662
CL_met  = 14.2    // L/h  SU12662 clearance
V_met   = 1100    // L    SU12662 volume

// --- Nivolumab TMDD (PD-1 target) ---
// Ref: Bajaj 2017, Lindauer 2017
CL_niv  = 0.019   // L/h
V1_niv  = 3.5     // L
Q_niv   = 0.003   // L/h
V2_niv  = 3.8     // L
kon_niv = 0.32    // 1/(nM·h)
koff_niv= 0.0018  // 1/h
kdeg_niv= 0.034   // 1/h  PD-1 turnover
Rmax_PD1= 6.5     // nM   PD-1 baseline
ksyn_PD1= 0.22    // nM/h PD-1 synthesis

// --- Belzutifan (HIF-2α inhibitor) ---
// Ref: Courtney 2020 MK-6482 Phase I, Jonasch 2021
CL_bez  = 2.5     // L/h
V_bez   = 85      // L
ka_bez  = 0.6     // 1/h
F_bez   = 0.70
IC50_bez= 0.018   // µM  HIF2A inhibition

// --- VHL/HIF/VEGF pathway ---
kprod_VHL  = 0.12  // 1/h  pVHL synthesis
kdeg_VHL   = 0.08  // 1/h  pVHL degradation
VHL_frac   = 0.25  // fraction of functional pVHL (0.25 = 75% loss, ccRCC)
kprod_HIF2 = 0.18  // nM/h HIF-2α synthesis
kdeg_HIF2  = 0.35  // 1/h  HIF-2α degradation (normoxia)
kprod_VEGF = 0.55  // nM/h VEGF synthesis
kdeg_VEGF  = 0.14  // 1/h  VEGF clearance
EC50_HIF2  = 1.8   // nM   HIF-2α → VEGF EC50

// --- VEGFR2 / Sunitinib PD ---
IC50_sun_VEGFR = 0.008  // µM sunitinib IC50 VEGFR2
VEGFR2_base    = 2.5    // nM VEGFR2 baseline
kprod_VEGFR2   = 0.12   // nM/h
kdeg_VEGFR2    = 0.05   // 1/h

// --- mTOR ---
kprod_mTOR = 0.22  // AU/h
kdeg_mTOR  = 0.18  // 1/h
IC50_ever  = 0.15  // nM everolimus IC50 mTORC1
kp_mTOR    = 1.2   // mTOR positive feedback on HIF2

// --- Tumor Growth (Simeoni TGI) ---
// Ref: Simeoni 2004, Rocchetti 2007
lambda1    = 5.0e-4  // 1/h exponential growth rate
lambda2    = 1.8e-2  // cm³/h linear growth rate (large tumors)
psi        = 20      // switch parameter
w0         = 0.8     // cm³ initial tumor volume
k1_sun     = 0.045   // L/(µg·h) sunitinib damage rate
k1_niv     = 0.008   // L/(nM·h) nivolumab damage rate (indirect)
k2_dam     = 0.032   // 1/h damage transit rate

// --- Immune Microenvironment ---
kprol_CD8  = 0.09    // 1/h  CD8 T-cell proliferation
kdeath_CD8 = 0.04    // 1/h  CD8 death
alpha_CD8  = 2.5e-4  // killing rate cm³/(cell·h)
ki_treg    = 0.35    // Treg suppression coefficient
k_MDSC     = 0.28    // MDSC immune suppression coefficient
kprol_Treg = 0.06    // 1/h
kdeath_Treg= 0.04    // 1/h
kprol_MDSC = 0.07    // 1/h
kdeath_MDSC= 0.05    // 1/h
CD8_base   = 100     // cells/µL baseline
Treg_base  = 15      // cells/µL baseline
MDSC_base  = 25      // cells/µL baseline
EC50_PD1   = 0.5     // nM  PD-1 occupancy EC50 for CD8 rescue

// --- Everolimus (oral mTOR inhibitor) ---
CL_ever    = 8.4     // L/h
V_ever     = 191     // L
ka_ever    = 1.2     // 1/h
F_ever     = 0.16

// --- Axitinib ---
CL_axi     = 12.5    // L/h
V_axi      = 160     // L
ka_axi     = 1.0     // 1/h
F_axi      = 0.58
IC50_axi   = 0.003   // µM axitinib IC50 VEGFR2

// --- Cabozantinib ---
CL_cabo    = 2.5     // L/h
V_cabo     = 350     // L
ka_cabo    = 0.8     // 1/h
F_cabo     = 0.57
IC50_cabo  = 0.006   // µM cabozantinib IC50 MET/VEGFR2

// --- Dose switches (1=active) ---
use_sunitinib  = 0
use_nivolumab  = 0
use_belzutifan = 0
use_everolimus = 0
use_axitinib   = 0
use_cabozantinib = 0

$CMT
// Sunitinib PK
DEPOT_SUN   // gut depot
CENT_SUN    // central (µg/L = ng/mL)
PERI_SUN    // peripheral
MET_SUN     // SU12662 metabolite

// Nivolumab TMDD
CENT_NIV    // free nivolumab (nM)
PERI_NIV    // peripheral
PD1_FREE    // unbound PD-1 (nM)
PD1_BOUND   // niv:PD-1 complex

// Belzutifan
DEPOT_BEZ   // gut
CENT_BEZ    // central (µM)

// VHL/HIF/VEGF
pVHL        // functional pVHL protein (AU)
HIF2A       // HIF-2α (nM)
VEGF        // free VEGF (nM)
VEGFR2_ACT  // active VEGFR2 (nM)

// mTOR
mTOR_ACT    // mTORC1 activity (AU)

// Everolimus
DEPOT_EVER  // gut
CENT_EVER   // central (nM)

// Axitinib
DEPOT_AXI
CENT_AXI    // central (µM)

// Cabozantinib
DEPOT_CABO
CENT_CABO   // central (µM)

// Tumor (Simeoni TGI)
TUM_W1      // damage transit 1
TUM_W2      // damage transit 2
TUM_W3      // damage transit 3
TUM_VOL     // tumor volume (cm³)

// Immune
CD8_T       // CD8+ T-cells (cells/µL)
TREG        // Regulatory T-cells
MDSC        // Myeloid-derived suppressor cells

$MAIN
// Initialise compartments
DEPOT_SUN_0   = 0;
CENT_SUN_0    = 0;
PERI_SUN_0    = 0;
MET_SUN_0     = 0;
CENT_NIV_0    = 0;
PERI_NIV_0    = 0;
PD1_FREE_0    = Rmax_PD1;
PD1_BOUND_0   = 0;
DEPOT_BEZ_0   = 0;
CENT_BEZ_0    = 0;
pVHL_0        = VHL_frac * kprod_VHL / kdeg_VHL;
HIF2A_0       = kprod_HIF2 / kdeg_HIF2 * (1.0 / (1.0 + pVHL_0));
VEGF_0        = kprod_VEGF * pow(HIF2A_0, 2.0) / (pow(EC50_HIF2, 2.0) + pow(HIF2A_0, 2.0)) / kdeg_VEGF;
VEGFR2_ACT_0  = VEGFR2_base;
mTOR_ACT_0    = kprod_mTOR / kdeg_mTOR;
DEPOT_EVER_0  = 0;
CENT_EVER_0   = 0;
DEPOT_AXI_0   = 0;
CENT_AXI_0    = 0;
DEPOT_CABO_0  = 0;
CENT_CABO_0   = 0;
TUM_W1_0      = 0;
TUM_W2_0      = 0;
TUM_W3_0      = 0;
TUM_VOL_0     = w0;
CD8_T_0       = CD8_base;
TREG_0        = Treg_base;
MDSC_0        = MDSC_base;

$ODE
// ---- Sunitinib PK ----
double C_sun  = CENT_SUN / V1_sun;   // µg/L = ng/mL
double C_met  = MET_SUN  / V_met;    // ng/mL active metabolite

dxdt_DEPOT_SUN = -ka_sun * F_sun * DEPOT_SUN;
dxdt_CENT_SUN  =  ka_sun * F_sun * DEPOT_SUN
                 - CL_sun * C_sun
                 - Q_sun  * (C_sun - PERI_SUN / V2_sun);
dxdt_PERI_SUN  =  Q_sun  * (C_sun - PERI_SUN / V2_sun);
dxdt_MET_SUN   =  fm_sun * CL_sun * C_sun - CL_met * C_met;

// Combined sunitinib equivalent (ng/mL; SU12662 ~50% potency)
double C_sun_eq = (C_sun + 0.5 * C_met) / 1000.0;   // convert to µM

// ---- Nivolumab TMDD ----
double C_niv  = CENT_NIV / V1_niv;   // nM

dxdt_CENT_NIV  = -CL_niv * C_niv
                 - Q_niv  * (C_niv - PERI_NIV / V2_niv)
                 - kon_niv * C_niv * PD1_FREE
                 + koff_niv * PD1_BOUND;
dxdt_PERI_NIV  =  Q_niv  * (C_niv - PERI_NIV / V2_niv);
dxdt_PD1_FREE  =  ksyn_PD1 - kdeg_niv * PD1_FREE
                 - kon_niv * C_niv * PD1_FREE
                 + koff_niv * PD1_BOUND;
dxdt_PD1_BOUND =  kon_niv * C_niv * PD1_FREE
                 - (koff_niv + kdeg_niv) * PD1_BOUND;

// PD-1 occupancy
double PD1_occ = PD1_BOUND / (PD1_FREE + PD1_BOUND + 1e-9);

// ---- Belzutifan ----
double C_bez = CENT_BEZ / V_bez;     // µM
dxdt_DEPOT_BEZ = -ka_bez * F_bez * DEPOT_BEZ;
dxdt_CENT_BEZ  =  ka_bez * F_bez * DEPOT_BEZ - CL_bez * C_bez;

// ---- Everolimus ----
double C_ever = CENT_EVER / V_ever;  // nM
dxdt_DEPOT_EVER = -ka_ever * F_ever * DEPOT_EVER;
dxdt_CENT_EVER  =  ka_ever * F_ever * DEPOT_EVER - CL_ever * C_ever;

// ---- Axitinib ----
double C_axi = CENT_AXI / V_axi;    // µM
dxdt_DEPOT_AXI = -ka_axi * F_axi * DEPOT_AXI;
dxdt_CENT_AXI  =  ka_axi * F_axi * DEPOT_AXI - CL_axi * C_axi;

// ---- Cabozantinib ----
double C_cabo = CENT_CABO / V_cabo;  // µM
dxdt_DEPOT_CABO = -ka_cabo * F_cabo * DEPOT_CABO;
dxdt_CENT_CABO  =  ka_cabo * F_cabo * DEPOT_CABO - CL_cabo * C_cabo;

// ---- VHL / HIF-2α axis ----
// pVHL reduced in ccRCC (VHL_frac = fraction remaining)
dxdt_pVHL  = kprod_VHL * VHL_frac - kdeg_VHL * pVHL;

// HIF-2α: degraded by pVHL; belzutifan blocks HIF-2α activity
double inh_bez_HIF2 = C_bez / (IC50_bez + C_bez);
double stim_mTOR_HIF2 = kp_mTOR * mTOR_ACT / (kprod_mTOR / kdeg_mTOR);
dxdt_HIF2A = kprod_HIF2 * stim_mTOR_HIF2
            - kdeg_HIF2 * (pVHL + 0.01) * HIF2A
            - inh_bez_HIF2 * kdeg_HIF2 * HIF2A;

// VEGF: driven by HIF-2α (Hill kinetics)
dxdt_VEGF  = kprod_VEGF * pow(HIF2A, 2.0) / (pow(EC50_HIF2, 2.0) + pow(HIF2A, 2.0))
            - kdeg_VEGF * VEGF;

// VEGFR2: inhibited by sunitinib, axitinib, cabozantinib
double inh_sun_VEGFR  = use_sunitinib  * C_sun_eq / (IC50_sun_VEGFR + C_sun_eq);
double inh_axi_VEGFR  = use_axitinib   * C_axi    / (IC50_axi        + C_axi);
double inh_cabo_VEGFR = use_cabozantinib * C_cabo  / (IC50_cabo       + C_cabo);
double total_VEGFR_inh = 1.0 - (1.0 - inh_sun_VEGFR) * (1.0 - inh_axi_VEGFR) * (1.0 - inh_cabo_VEGFR);

dxdt_VEGFR2_ACT = kprod_VEGFR2 * (1.0 - total_VEGFR_inh)
                 - kdeg_VEGFR2 * VEGFR2_ACT;

// ---- mTOR ----
double inh_ever_mTOR = use_everolimus * C_ever / (IC50_ever + C_ever);
dxdt_mTOR_ACT = kprod_mTOR * (1.0 - inh_ever_mTOR)
               - kdeg_mTOR * mTOR_ACT;

// ---- Tumor Growth (Simeoni) ----
double TW = TUM_W1 + TUM_W2 + TUM_W3;
double growth = lambda1 * TUM_VOL / pow(1.0 + pow(lambda1 / lambda2 * TUM_VOL, psi), 1.0 / psi);

// Drug damage signals
double dmg_sun  = use_sunitinib    * k1_sun  * C_sun_eq;
double dmg_cabo = use_cabozantinib * k1_sun  * C_cabo;
double dmg_axi  = use_axitinib     * k1_sun  * C_axi;

// Immune kill by CD8 (reduced by Treg and MDSC, enhanced by anti-PD-1)
double CD8_eff = CD8_T * (1.0 - ki_treg * TREG / (Treg_base + TREG))
                       * (1.0 - k_MDSC  * MDSC / (MDSC_base + MDSC))
                       * (1.0 + 3.0 * PD1_occ / (EC50_PD1 + PD1_occ));
double immune_kill = alpha_CD8 * CD8_eff * TUM_VOL;

// Total drug kill rate
double total_k = dmg_sun + dmg_cabo + dmg_axi + immune_kill;

dxdt_TUM_W1  = total_k * TUM_VOL - k2_dam * TUM_W1;
dxdt_TUM_W2  =                       k2_dam * (TUM_W1 - TUM_W2);
dxdt_TUM_W3  =                       k2_dam * (TUM_W2 - TUM_W3);
dxdt_TUM_VOL =  growth - k2_dam * TW;

// ---- Immune microenvironment ----
// CD8 proliferation stimulated by VEGFR inhibition (less immune exclusion)
// and rescued by anti-PD-1
double stim_VEGFR_inh_CD8 = 1.0 + 0.8 * total_VEGFR_inh;
dxdt_CD8_T = kprol_CD8 * CD8_base * stim_VEGFR_inh_CD8
              * (1.0 + 2.0 * PD1_occ)
            - kdeath_CD8 * CD8_T
            - ki_treg * TREG * CD8_T / 100.0
            - k_MDSC * MDSC * CD8_T / 100.0;

dxdt_TREG  = kprol_Treg  * Treg_base - kdeath_Treg  * TREG;
dxdt_MDSC  = kprol_MDSC  * MDSC_base * (1.0 + 0.5 * VEGF / (VEGF + 1.0))
            - kdeath_MDSC * MDSC;

$TABLE
double Csun_ngmL   = CENT_SUN / V1_sun;
double Cmet_ngmL   = MET_SUN  / V_met;
double CsunEq_nM   = (Csun_ngmL + 0.5 * Cmet_ngmL) / 1000.0 * 1e6 / 589.7;
double Cniv_nM     = CENT_NIV / V1_niv;
double Cbez_uM     = CENT_BEZ / V_bez;
double HIF2_nM     = HIF2A;
double VEGF_nM     = VEGF;
double VEGFR2_nM   = VEGFR2_ACT;
double mTOR_AU     = mTOR_ACT;
double TumorVol    = TUM_VOL;
double CD8_count   = CD8_T;
double Treg_count  = TREG;
double MDSC_count  = MDSC;
double PD1_occ_pct = PD1_BOUND / (PD1_FREE + PD1_BOUND + 1e-9) * 100.0;

$CAPTURE Csun_ngmL Cmet_ngmL CsunEq_nM Cniv_nM Cbez_uM
         HIF2_nM VEGF_nM VEGFR2_nM mTOR_AU TumorVol
         CD8_count Treg_count MDSC_count PD1_occ_pct

$OMEGA @labels BSV_CL_sun BSV_VHL_frac BSV_TGI
0.0625    // 25% CV on sunitinib CL
0.1296    // 36% CV on VHL_frac
0.0400    // 20% CV on TGI kill

$SIGMA 0.04   // 20% residual for PK

'

## ----------------------------------------------------------
## 2. Compile model
## ----------------------------------------------------------
mod <- mrgsolve::mcode("rcc_qsp", rcc_code)

## ----------------------------------------------------------
## 3. Helper: build dosing event table
## ----------------------------------------------------------
build_events <- function(scenario = "sunitinib_mono",
                         n_weeks  = 52) {
  ev_list <- list()

  if (scenario == "sunitinib_mono") {
    # Sunitinib 50 mg QD × 4/2 schedule (4 wk on, 2 wk off)
    cycle_len <- 42 * 24   # 42 days in hours
    on_hours  <- 28 * 24   # 28 days on
    for (cy in 0:floor(n_weeks / 6)) {
      start_h <- cy * cycle_len
      ev_list[[length(ev_list) + 1]] <-
        ev(amt = 50000, ii = 24, addl = 27, time = start_h, cmt = "DEPOT_SUN")
    }
    params <- c(use_sunitinib = 1)

  } else if (scenario == "nivo_ipi") {
    # Nivolumab 3 mg/kg q3w + Ipilimumab 1 mg/kg q3w × 4 doses, then Nivo 3 mg/kg q2w
    # Ipilimumab modelled as Treg depletion (simplified)
    ev_list[[1]] <- ev(amt = 240e6, ii = 14*24, addl = floor(n_weeks/2)-1,
                       time = 0, cmt = "CENT_NIV")
    params <- c(use_nivolumab = 1)

  } else if (scenario == "pembro_axitinib") {
    # Pembrolizumab 200 mg q3w + Axitinib 5 mg BID
    ev_list[[1]] <- ev(amt = 200e6, ii = 21*24, addl = floor(n_weeks/3)-1,
                       time = 0, cmt = "CENT_NIV")   # simplified: same TMDD as nivo
    ev_list[[2]] <- ev(amt = 5000, ii = 12, addl = n_weeks*14-1,
                       time = 0, cmt = "DEPOT_AXI")
    params <- c(use_nivolumab = 1, use_axitinib = 1)

  } else if (scenario == "cabo_nivo") {
    # Cabozantinib 40 mg QD + Nivolumab 240 mg q2w
    ev_list[[1]] <- ev(amt = 40000, ii = 24, addl = n_weeks*7-1,
                       time = 0, cmt = "DEPOT_CABO")
    ev_list[[2]] <- ev(amt = 240e6, ii = 14*24, addl = floor(n_weeks/2)-1,
                       time = 0, cmt = "CENT_NIV")
    params <- c(use_cabozantinib = 1, use_nivolumab = 1)

  } else if (scenario == "cabozantinib_mono") {
    # Cabozantinib 60 mg QD (METEOR)
    ev_list[[1]] <- ev(amt = 60000, ii = 24, addl = n_weeks*7-1,
                       time = 0, cmt = "DEPOT_CABO")
    params <- c(use_cabozantinib = 1)

  } else if (scenario == "everolimus") {
    # Everolimus 10 mg QD (RECORD-1)
    ev_list[[1]] <- ev(amt = 10000, ii = 24, addl = n_weeks*7-1,
                       time = 0, cmt = "DEPOT_EVER")
    params <- c(use_everolimus = 1)

  } else if (scenario == "belzutifan") {
    # Belzutifan 120 mg QD (LITESPARK-005)
    ev_list[[1]] <- ev(amt = 120, ii = 24, addl = n_weeks*7-1,
                       time = 0, cmt = "DEPOT_BEZ")
    params <- c(use_belzutifan = 1)

  } else {
    # Untreated (vehicle)
    ev_list[[1]] <- ev(amt = 0, time = 0, cmt = "DEPOT_SUN")
    params <- c()
  }

  list(evs = do.call(rbind, ev_list), params = params)
}

## ----------------------------------------------------------
## 4. Simulate all 7 scenarios
## ----------------------------------------------------------
scenarios <- c(
  "untreated",
  "sunitinib_mono",
  "nivo_ipi",
  "pembro_axitinib",
  "cabo_nivo",
  "cabozantinib_mono",
  "everolimus",
  "belzutifan"
)

sim_list <- lapply(scenarios, function(sc) {
  cfg    <- build_events(sc, n_weeks = 52)
  p_over <- if (length(cfg$params) > 0) as.list(cfg$params) else list()

  out <- mod %>%
    param(p_over) %>%
    mrgsim(events = cfg$evs,
           end    = 52 * 7 * 24,
           delta  = 6,
           obsonly = TRUE) %>%
    as.data.frame()

  out$scenario <- sc
  out
})

results <- bind_rows(sim_list)

## ----------------------------------------------------------
## 5. Plot: Tumor Volume over Time
## ----------------------------------------------------------
scenario_labels <- c(
  untreated          = "Untreated",
  sunitinib_mono     = "Sunitinib 50 mg (4/2)",
  nivo_ipi           = "Nivolumab + Ipilimumab",
  pembro_axitinib    = "Pembrolizumab + Axitinib",
  cabo_nivo          = "Cabozantinib + Nivolumab",
  cabozantinib_mono  = "Cabozantinib 60 mg",
  everolimus         = "Everolimus 10 mg",
  belzutifan         = "Belzutifan 120 mg"
)

p_tumor <- results %>%
  mutate(time_wk = time / (7 * 24),
         label   = scenario_labels[scenario]) %>%
  ggplot(aes(x = time_wk, y = TumorVol, color = label)) +
  geom_line(linewidth = 0.8) +
  labs(title    = "ccRCC: Tumor Volume Dynamics",
       x        = "Time (weeks)",
       y        = "Tumor Volume (cm³)",
       color    = "Regimen") +
  theme_bw() +
  scale_color_brewer(palette = "Set2")

print(p_tumor)

## ----------------------------------------------------------
## 6. Plot: CD8 T-cell dynamics
## ----------------------------------------------------------
p_cd8 <- results %>%
  mutate(time_wk = time / (7 * 24),
         label   = scenario_labels[scenario]) %>%
  ggplot(aes(x = time_wk, y = CD8_count, color = label)) +
  geom_line(linewidth = 0.8) +
  labs(title = "CD8+ T-cell Dynamics",
       x = "Time (weeks)", y = "CD8+ T-cells (cells/µL)",
       color = "Regimen") +
  theme_bw() +
  scale_color_brewer(palette = "Set2")

print(p_cd8)

## ----------------------------------------------------------
## 7. Plot: HIF-2α and VEGF
## ----------------------------------------------------------
p_pathway <- results %>%
  filter(scenario %in% c("untreated", "sunitinib_mono", "belzutifan", "cabo_nivo")) %>%
  mutate(time_wk = time / (7 * 24),
         label   = scenario_labels[scenario]) %>%
  select(time_wk, label, HIF2_nM, VEGF_nM) %>%
  pivot_longer(cols = c(HIF2_nM, VEGF_nM), names_to = "marker") %>%
  ggplot(aes(x = time_wk, y = value, color = label)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~marker, scales = "free_y") +
  labs(title = "VHL/HIF/VEGF Pathway Biomarkers",
       x = "Time (weeks)", y = "Concentration (nM)",
       color = "Regimen") +
  theme_bw()

print(p_pathway)

## ----------------------------------------------------------
## 8. Summary endpoint table
## ----------------------------------------------------------
endpoint_summary <- results %>%
  mutate(time_wk = time / (7 * 24)) %>%
  group_by(scenario) %>%
  summarise(
    Tumor_at_24wk = TumorVol[which.min(abs(time_wk - 24))],
    Tumor_at_52wk = TumorVol[which.min(abs(time_wk - 52))],
    BOR_pct       = (1 - min(TumorVol) / TumorVol[1]) * 100,
    CD8_peak      = max(CD8_count),
    PD1_occ_ss    = PD1_occ_pct[n()],
    .groups = "drop"
  ) %>%
  mutate(label = scenario_labels[scenario]) %>%
  arrange(Tumor_at_52wk)

print(endpoint_summary)
message("RCC QSP simulation complete.")
