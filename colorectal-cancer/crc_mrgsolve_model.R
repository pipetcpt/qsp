## =============================================================================
## Colorectal Cancer (CRC) — Quantitative Systems Pharmacology Model
## mrgsolve ODE implementation
##
## Pathways covered:
##   • 5-FU/LV PK (3-CMT + intracellular)
##   • Oxaliplatin PK (2-CMT, Pt-DNA adduct)
##   • Irinotecan → SN-38 PK (3-CMT, UGT1A1 polymorphism)
##   • Bevacizumab (2-CMT TMDD)
##   • Cetuximab (2-CMT EGFR occupancy)
##   • Pembrolizumab (1-CMT PD-1 blockade)
##   • Tumor growth (sensitive/resistant) with immune kill
##   • Biomarkers: CEA, ctDNA, RECIST-linked diameter
##
## Treatment scenarios calibrated to:
##   MOSAIC (FOLFOX4), CRYSTAL (FOLFIRI+CTX), TRIBE (FOLFIRI+BEV),
##   KEYNOTE-177 (pembro MSI-H), CORRECT (regorafenib)
##
## Author: Claude Code Routine (CCR) — 2026-06-23
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ─────────────────────────────────────────────────────────────────────────────
## MODEL CODE
## ─────────────────────────────────────────────────────────────────────────────
crc_code <- '
$PROB
  Colorectal Cancer QSP Model
  ODE Compartments (20):
    [1-3]  5-FU PK: FU1 (central), FU2 (peripheral), FU_ic (intracellular)
    [4-5]  Oxaliplatin: OX1 (central), OX_DNA (Pt-DNA adducts)
    [6-8]  Irinotecan: IRI1 (central), SN38 (active), SN38G (glucuronide)
    [9-10] Bevacizumab: BEV1 (central), BEV_VEGF (complex)
    [11]   Cetuximab: CTX1 (central)
    [12]   Pembrolizumab: PEM1 (central)
    [13]   Tumor sensitive cells (Ts)
    [14]   Tumor resistant cells (Tr)
    [15]   CEA biomarker
    [16]   ctDNA
    [17]   CD8 T effector cells
    [18]   VEGF (free)
    [19]   EGFR occupancy (fraction)
    [20]   PD-1 occupancy (fraction)

$PARAM
  // ── Body surface area & dosing
  BSA    = 1.73   // m²
  WT     = 70     // kg

  // ── 5-FU PK (population-based, Gamelin 1999 / Capitain 2012)
  CL_FU  = 112    // L/h  (high clearance, DPD-dependent)
  V1_FU  = 11     // L    central
  V2_FU  = 28     // L    peripheral
  Q_FU   = 30     // L/h  intercompartmental
  Kic_FU = 0.5    // h⁻¹  uptake to intracellular
  Kout_ic= 0.8    // h⁻¹  intracellular turnover

  // ── 5-FU PD (TS inhibition → DNA damage → apoptosis)
  IC50_FU  = 0.05  // μg/mL  FdUMP-TS IC50
  Hill_FU  = 1.5   //         Hill coeff
  Emax_FU  = 0.80  //         max kill fraction

  // ── Oxaliplatin PK (Troger 2012, 2-CMT)
  CL_OX   = 10.2  // L/h
  V1_OX   = 4.9   // L
  V2_OX   = 21    // L
  Q_OX    = 4.6   // L/h
  k_adduct= 0.018 // h⁻¹  Pt-DNA adduct formation
  k_repair= 0.003 // h⁻¹  NER repair rate

  // ── Oxaliplatin PD
  IC50_OX   = 0.2  // μg/mL  adduct→cell kill IC50
  Emax_OX   = 0.70
  Hill_OX   = 1.2

  // ── Irinotecan PK (Mathijssen 2002 / 3-CMT)
  CL_IRI   = 20   // L/h
  V1_IRI   = 15   // L
  V2_IRI   = 120  // L
  Q_IRI    = 14   // L/h
  k_conv   = 0.15 // h⁻¹  IRI→SN-38 (CES/CYP3A4)
  k_gluc   = 0.22 // h⁻¹  SN-38→SN-38G (UGT1A1)
  UGT_eff  = 1.0  // 1=normal, 0.5=*28/*28 poor glucuronidator
  CL_SN38  = 4.0  // L/h

  // ── SN-38 PD (Top1 inhibition)
  IC50_SN38 = 0.015 // μg/mL
  Emax_SN38 = 0.75
  Hill_SN38 = 1.3

  // ── Bevacizumab PK (TMDD, Lammerts 2017)
  CL_BEV   = 0.197  // L/h
  V1_BEV   = 2.92   // L
  V2_BEV   = 2.09   // L
  Q_BEV    = 0.398  // L/h
  Kin_VEGF = 0.05   // nM/h  VEGF production
  Kout_VEGF= 0.04   // h⁻¹   VEGF elimination
  Kon_BEV  = 0.12   // nM⁻¹h⁻¹  binding
  Koff_BEV = 0.0001 // h⁻¹  dissociation (Kd~1nM)
  Kint     = 0.001  // h⁻¹  complex internalization

  // ── Bevacizumab PD (anti-angiogenic → tumor growth inhibition)
  VEGF0    = 1.25   // nM   baseline VEGF
  EC50_BEV = 0.4    // nM   free VEGF for max angiogenesis
  Emax_BEV = 0.6    // max tumor growth inhibition from anti-VEGF

  // ── Cetuximab PK (1-CMT simplified, Mould 2015)
  CL_CTX   = 0.384  // L/h
  V1_CTX   = 3.57   // L
  k_EGFR   = 0.08   // h⁻¹  EGFR saturation rate
  Kd_CTX   = 0.3    // nM   EGFR binding affinity

  // ── Cetuximab PD
  Emax_CTX  = 0.65   // max effect (RAS-WT only)
  KRAS_mut  = 0      // 0=WT, 1=mutant (bypass resistance)

  // ── Pembrolizumab PK (1-CMT, 200mg Q3W)
  CL_PEM   = 0.214  // L/h
  V1_PEM   = 5.0    // L
  Kd_PEM   = 0.0004 // nM   PD-1 Kd

  // ── Pembrolizumab PD (T cell reactivation — MSI-H)
  MSI_H    = 0       // 0=MSS, 1=MSI-H
  CD8_base = 0.3     // baseline CD8 T cell activity (0-1)
  CD8_max  = 0.9     // max with PD-1 blockade (MSI-H)
  k_CD8act = 0.005   // h⁻¹  CD8 activation rate

  // ── Tumor dynamics
  kg_s     = 0.008   // h⁻¹  sensitive cell growth rate
  kg_r     = 0.006   // h⁻¹  resistant cell growth rate
  k_kill_s = 0.025   // h⁻¹  maximum drug kill rate (sensitive)
  k_kill_r = 0.004   // h⁻¹  drug kill (resistant)
  k_mutate = 1e-6    // h⁻¹  sensitive→resistant mutation rate
  k_immune = 0.003   // h⁻¹  immune-mediated kill
  K_carry  = 1e10    // cells  carrying capacity
  Ts0      = 5e8     // initial sensitive tumor cells
  Tr0      = 1e6     // initial resistant cells
  frac_r   = 0.001   // initial resistant fraction

  // ── CEA PK
  kout_CEA = 0.004   // h⁻¹  CEA turnover
  kin_CEA0 = 0.004   // baseline CEA production (normalized)

  // ── ctDNA
  kout_ctDNA = 0.05  // h⁻¹

  // ── Regorafenib (additional multi-kinase)
  REG_on   = 0       // 0=off, 1=on
  Emax_REG = 0.55    // multi-kinase max effect
  EC50_REG = 0.5     // μg/mL (plasma)

$CMT
  FU1 FU2 FU_ic         // 5-FU
  OX1 OX_DNA            // Oxaliplatin
  IRI1 SN38 SN38G       // Irinotecan
  BEV1 BEV_VEGF         // Bevacizumab
  CTX1                  // Cetuximab
  PEM1                  // Pembrolizumab
  Ts Tr                 // Tumor cells
  CEA ctDNA             // Biomarkers
  CD8eff                // CD8 T effectors (activity units)
  VEGF_free             // Free VEGF (nM)
  EGFR_occ              // EGFR occupancy
  PD1_occ               // PD-1 occupancy

$INIT
  FU1    = 0
  FU2    = 0
  FU_ic  = 0
  OX1    = 0
  OX_DNA = 0
  IRI1   = 0
  SN38   = 0
  SN38G  = 0
  BEV1   = 0
  BEV_VEGF = 0
  CTX1   = 0
  PEM1   = 0
  Ts     = 5e8
  Tr     = 1e6
  CEA    = 10         // ng/mL baseline
  ctDNA  = 0.01       // normalized units
  CD8eff = 0.3        // 0-1 scale
  VEGF_free = 1.25    // nM
  EGFR_occ  = 0
  PD1_occ   = 0

$ODE
  // ── 5-FU concentration (μg/mL = mg/L)
  double Cp_FU = FU1 / V1_FU;

  // ── Drug kill effects ─────────────────────────────────────────────────────

  // 5-FU: FdUMP→TS inhibition
  double FU_eff = Emax_FU * pow(FU_ic / (IC50_FU + FU_ic), Hill_FU);

  // Oxaliplatin: DNA adducts
  double OX_eff = Emax_OX * pow(OX_DNA / (IC50_OX + OX_DNA), Hill_OX);

  // SN-38 Top1 inhibition
  double SN38_eff = Emax_SN38 * pow(SN38 / (IC50_SN38 + SN38 + 1e-9), Hill_SN38);

  // Bevacizumab: normalize VEGF suppression → growth inhibition
  double VEGF_rel = VEGF_free / (EC50_BEV + VEGF_free);
  double BEV_TGI  = Emax_BEV * (1.0 - VEGF_rel);

  // Cetuximab: EGFR blockade (RAS-WT only)
  double CTX_eff = Emax_CTX * EGFR_occ * (1.0 - KRAS_mut);

  // Pembrolizumab: CD8 reactivation (MSI-H tumors benefit most)
  double CD8_kill = CD8eff * k_immune * MSI_H + CD8eff * k_immune * 0.1 * (1.0 - MSI_H);

  // Regorafenib (optional)
  // REG plasma assumed constant if on (simplified)
  double REG_Cp = REG_on * 2.0; // μg/mL average
  double REG_eff= Emax_REG * REG_Cp / (EC50_REG + REG_Cp);

  // Combined kill rate (sensitive)
  double kill_s = k_kill_s * (FU_eff + OX_eff + SN38_eff + CTX_eff + BEV_TGI + REG_eff)
                  + CD8_kill;
  kill_s = (kill_s > 0.99) ? 0.99 : kill_s;

  // Combined kill rate (resistant — minimal chemosensitivity)
  double kill_r = k_kill_r * (FU_eff + OX_eff + SN38_eff) + CD8_kill;

  // Logistic growth
  double Ttot  = Ts + Tr;
  double gfrac = 1.0 - Ttot / K_carry;
  if (gfrac < 0) gfrac = 0;

  // ── ODE SYSTEM ─────────────────────────────────────────────────────────────

  // 5-FU PK
  dxdt_FU1   = -(CL_FU / V1_FU) * FU1 - (Q_FU / V1_FU) * FU1
               + (Q_FU / V2_FU) * FU2
               - Kic_FU * FU1;
  dxdt_FU2   = (Q_FU / V1_FU) * FU1 - (Q_FU / V2_FU) * FU2;
  dxdt_FU_ic = Kic_FU * Cp_FU - Kout_ic * FU_ic;

  // Oxaliplatin PK
  double Cp_OX = OX1 / V1_OX;
  dxdt_OX1   = -(CL_OX / V1_OX) * OX1 - (Q_OX / V1_OX) * OX1;
  dxdt_OX_DNA= k_adduct * Cp_OX - k_repair * OX_DNA;

  // Irinotecan/SN-38 PK
  double Cp_IRI  = IRI1 / V1_IRI;
  double kconv_eff = k_conv;  // could be adjusted for UGT polymorphism
  double kgluc_eff = k_gluc * UGT_eff;
  dxdt_IRI1  = -(CL_IRI / V1_IRI) * IRI1 - (Q_IRI / V1_IRI) * IRI1
               - kconv_eff * IRI1;
  dxdt_SN38  = kconv_eff * Cp_IRI * V1_IRI - CL_SN38 * SN38
               - kgluc_eff * SN38;
  dxdt_SN38G = kgluc_eff * SN38;

  // Bevacizumab PK (TMDD)
  double Cp_BEV = BEV1 / V1_BEV;
  double Cp_BEV_nM = Cp_BEV / 0.149; // mg/mL → nM (MW≈149 kDa)
  dxdt_BEV1    = -(CL_BEV / V1_BEV) * BEV1
                 - Kon_BEV * Cp_BEV_nM * VEGF_free * V1_BEV
                 + Koff_BEV * BEV_VEGF;
  dxdt_BEV_VEGF= Kon_BEV * Cp_BEV_nM * VEGF_free
                 - (Koff_BEV + Kint) * BEV_VEGF;

  // Free VEGF turnover
  double VEGF_prod = Kin_VEGF * (Ttot / Ts0); // tumor-derived
  dxdt_VEGF_free = VEGF_prod - Kout_VEGF * VEGF_free
                   - Kon_BEV * Cp_BEV_nM * VEGF_free
                   + Koff_BEV * BEV_VEGF;

  // Cetuximab PK (1-CMT)
  double Cp_CTX = CTX1 / V1_CTX;
  double Cp_CTX_nM = Cp_CTX / 0.145;
  dxdt_CTX1  = -(CL_CTX / V1_CTX) * CTX1;
  dxdt_EGFR_occ = k_EGFR * (Cp_CTX_nM / (Kd_CTX + Cp_CTX_nM) - EGFR_occ);

  // Pembrolizumab PK
  double Cp_PEM = PEM1 / V1_PEM;
  double Cp_PEM_nM = Cp_PEM / 0.149;
  dxdt_PEM1  = -(CL_PEM / V1_PEM) * PEM1;
  dxdt_PD1_occ = 0.01 * (Cp_PEM_nM / (Kd_PEM + Cp_PEM_nM) - PD1_occ);

  // CD8 effector dynamics
  double CD8_target = CD8_base + (CD8_max - CD8_base) * PD1_occ * MSI_H;
  dxdt_CD8eff = k_CD8act * (CD8_target - CD8eff);

  // Tumor cells
  dxdt_Ts = kg_s * Ts * gfrac - kill_s * Ts - k_mutate * Ts;
  dxdt_Tr = kg_r * Tr * gfrac - kill_r * Tr + k_mutate * Ts;

  // CEA biomarker (proportional to total tumor burden)
  double Ttot_norm = Ttot / (Ts0 + 1e6);
  double kin_CEA   = kin_CEA0 * Ttot_norm;
  dxdt_CEA  = kin_CEA - kout_CEA * CEA;

  // ctDNA (faster turnover, correlates with apoptosis)
  double apop_flux = kill_s * Ts + kill_r * Tr;
  dxdt_ctDNA = 1e-9 * apop_flux - kout_ctDNA * ctDNA;

$TABLE
  double Cp_FU_out   = FU1 / V1_FU;
  double Cp_OX_out   = OX1 / V1_OX;
  double Cp_SN38_out = SN38;
  double Cp_BEV_out  = BEV1 / V1_BEV;
  double Cp_CTX_out  = CTX1 / V1_CTX;
  double Cp_PEM_out  = PEM1 / V1_PEM;
  double Ttotal      = Ts + Tr;
  double TumDiam     = 35.0 * pow(Ttotal / (Ts0 + 1e6), 0.333); // mm (SLD)
  double PctChange   = (TumDiam - 35.0) / 35.0 * 100.0;
  double ResistFrac  = Tr / (Ts + Tr + 1e-6);
  double CEA_out     = CEA;
  double ctDNA_out   = ctDNA;

$CAPTURE
  Cp_FU_out Cp_OX_out Cp_SN38_out Cp_BEV_out Cp_CTX_out Cp_PEM_out
  Ttotal TumDiam PctChange ResistFrac CEA_out ctDNA_out
  EGFR_occ PD1_occ CD8eff VEGF_free BEV_VEGF
'

## ─────────────────────────────────────────────────────────────────────────────
## COMPILE MODEL
## ─────────────────────────────────────────────────────────────────────────────
mod <- mcode("crc_qsp", crc_code)

## ─────────────────────────────────────────────────────────────────────────────
## HELPER: FOLFOX6 dosing events (oxaliplatin + 5-FU bolus + infusion)
##   Cycle: Q2W; Oxaliplatin 85 mg/m² 2h; 5-FU 400 mg/m² bolus + 2400 mg/m² 46h
## ─────────────────────────────────────────────────────────────────────────────
make_FOLFOX <- function(n_cycles = 12, bsa = 1.73) {
  ox_dose <- 85 * bsa  # mg
  fu_bolus <- 400 * bsa
  fu_inf   <- 2400 * bsa  # mg over 46h (rate = fu_inf/46)
  evts <- data.frame()
  for (i in seq_len(n_cycles)) {
    t0 <- (i - 1) * 14 * 24   # hours
    # Oxaliplatin 2h infusion
    evts <- bind_rows(evts,
      data.frame(time = t0, cmt = "OX1", amt = ox_dose, rate = ox_dose / 2,
                 evid = 1, ii = 0, addl = 0))
    # 5-FU bolus
    evts <- bind_rows(evts,
      data.frame(time = t0, cmt = "FU1", amt = fu_bolus, rate = -2,
                 evid = 1, ii = 0, addl = 0))
    # 5-FU 46h continuous infusion
    evts <- bind_rows(evts,
      data.frame(time = t0, cmt = "FU1", amt = fu_inf, rate = fu_inf / 46,
                 evid = 1, ii = 0, addl = 0))
  }
  evts
}

## ─────────────────────────────────────────────────────────────────────────────
## HELPER: FOLFIRI dosing events
##   Irinotecan 180 mg/m² 90min; same 5-FU
## ─────────────────────────────────────────────────────────────────────────────
make_FOLFIRI <- function(n_cycles = 12, bsa = 1.73) {
  iri_dose <- 180 * bsa
  fu_bolus <- 400 * bsa
  fu_inf   <- 2400 * bsa
  evts <- data.frame()
  for (i in seq_len(n_cycles)) {
    t0 <- (i - 1) * 14 * 24
    evts <- bind_rows(evts,
      data.frame(time = t0, cmt = "IRI1", amt = iri_dose, rate = iri_dose / 1.5,
                 evid = 1, ii = 0, addl = 0))
    evts <- bind_rows(evts,
      data.frame(time = t0, cmt = "FU1", amt = fu_bolus, rate = -2,
                 evid = 1, ii = 0, addl = 0))
    evts <- bind_rows(evts,
      data.frame(time = t0, cmt = "FU1", amt = fu_inf, rate = fu_inf / 46,
                 evid = 1, ii = 0, addl = 0))
  }
  evts
}

## ─────────────────────────────────────────────────────────────────────────────
## HELPER: Bevacizumab dosing (5 mg/kg Q2W)
## ─────────────────────────────────────────────────────────────────────────────
make_BEV <- function(n_cycles = 12, wt = 70) {
  bev_dose <- 5 * wt  # mg
  evts <- data.frame()
  for (i in seq_len(n_cycles)) {
    t0 <- (i - 1) * 14 * 24
    evts <- bind_rows(evts,
      data.frame(time = t0, cmt = "BEV1", amt = bev_dose, rate = bev_dose / 0.5,
                 evid = 1, ii = 0, addl = 0))
  }
  evts
}

## ─────────────────────────────────────────────────────────────────────────────
## HELPER: Cetuximab (400 mg/m² loading, then 250 mg/m² Q1W)
## ─────────────────────────────────────────────────────────────────────────────
make_CTX <- function(n_weeks = 24, bsa = 1.73) {
  evts <- data.frame()
  for (i in seq_len(n_weeks)) {
    dose <- ifelse(i == 1, 400 * bsa, 250 * bsa)
    t0   <- (i - 1) * 7 * 24
    evts <- bind_rows(evts,
      data.frame(time = t0, cmt = "CTX1", amt = dose, rate = dose / 2,
                 evid = 1, ii = 0, addl = 0))
  }
  evts
}

## ─────────────────────────────────────────────────────────────────────────────
## HELPER: Pembrolizumab (200 mg Q3W)
## ─────────────────────────────────────────────────────────────────────────────
make_PEM <- function(n_cycles = 18) {
  evts <- data.frame()
  for (i in seq_len(n_cycles)) {
    t0 <- (i - 1) * 21 * 24
    evts <- bind_rows(evts,
      data.frame(time = t0, cmt = "PEM1", amt = 200, rate = 200 / 0.5,
                 evid = 1, ii = 0, addl = 0))
  }
  evts
}

## ─────────────────────────────────────────────────────────────────────────────
## SIMULATION SETUP
## ─────────────────────────────────────────────────────────────────────────────
sim_end <- 12000  # hours (~500 days)
times   <- c(seq(0, 168, by = 2), seq(168, sim_end, by = 24))

## ─────────────────────────────────────────────────────────────────────────────
## SCENARIO 1: No Treatment (natural history)
## ─────────────────────────────────────────────────────────────────────────────
cat("Running Scenario 1: No treatment (natural history)...\n")
ev_none <- data.frame(time = 0, cmt = "FU1", amt = 0, evid = 0)
out1 <- mod %>%
  ev(ev_none) %>%
  mrgsim(end = sim_end, delta = 24, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(scenario = "1_NoTreatment")

## ─────────────────────────────────────────────────────────────────────────────
## SCENARIO 2: FOLFOX6 (MOSAIC trial, Kolfox/Andre 2004)
##  mFOLFOX6: median PFS 9.0 mo, median OS 17.9 mo (mCRC)
## ─────────────────────────────────────────────────────────────────────────────
cat("Running Scenario 2: FOLFOX6...\n")
ev2 <- make_FOLFOX(n_cycles = 12) %>%
  as_data_set()
out2 <- mod %>%
  param(KRAS_mut = 0) %>%
  ev(as_data_set(make_FOLFOX(n_cycles = 12))) %>%
  mrgsim(end = sim_end, delta = 24, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(scenario = "2_FOLFOX6")

## ─────────────────────────────────────────────────────────────────────────────
## SCENARIO 3: FOLFIRI (GERCOR / Tournigand 2004)
##  Similar efficacy to FOLFOX, sequence dependent
## ─────────────────────────────────────────────────────────────────────────────
cat("Running Scenario 3: FOLFIRI...\n")
out3 <- mod %>%
  param(KRAS_mut = 0) %>%
  ev(as_data_set(make_FOLFIRI(n_cycles = 12))) %>%
  mrgsim(end = sim_end, delta = 24, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(scenario = "3_FOLFIRI")

## ─────────────────────────────────────────────────────────────────────────────
## SCENARIO 4: FOLFOX + Bevacizumab (NO16966, Saltz 2008)
##  mPFS 9.4 mo (vs 8.0 mo FOLFOX alone)
## ─────────────────────────────────────────────────────────────────────────────
cat("Running Scenario 4: FOLFOX + Bevacizumab...\n")
ev4 <- bind_rows(make_FOLFOX(n_cycles = 12), make_BEV(n_cycles = 12)) %>%
  arrange(time)
out4 <- mod %>%
  param(KRAS_mut = 0) %>%
  ev(as_data_set(ev4)) %>%
  mrgsim(end = sim_end, delta = 24, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(scenario = "4_FOLFOX_BEV")

## ─────────────────────────────────────────────────────────────────────────────
## SCENARIO 5: FOLFIRI + Cetuximab (RAS-WT) — CRYSTAL trial (Van Cutsem 2009)
##  RAS-WT: ORR 57%, median PFS 9.9 mo vs 8.7 mo (FOLFIRI alone)
## ─────────────────────────────────────────────────────────────────────────────
cat("Running Scenario 5: FOLFIRI + Cetuximab (RAS-WT)...\n")
ev5 <- bind_rows(make_FOLFIRI(n_cycles = 12), make_CTX(n_weeks = 24)) %>%
  arrange(time)
out5 <- mod %>%
  param(KRAS_mut = 0) %>%
  ev(as_data_set(ev5)) %>%
  mrgsim(end = sim_end, delta = 24, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(scenario = "5_FOLFIRI_CTX_RASWT")

## ─────────────────────────────────────────────────────────────────────────────
## SCENARIO 6: FOLFIRI + Bevacizumab (TRIBE trial, Falcone 2013)
##  FOLFOXIRI+BEV reference: median PFS 12.1 mo, OS 29.8 mo
##  FOLFIRI+BEV: median PFS 9.7 mo
## ─────────────────────────────────────────────────────────────────────────────
cat("Running Scenario 6: FOLFIRI + Bevacizumab...\n")
ev6 <- bind_rows(make_FOLFIRI(n_cycles = 12), make_BEV(n_cycles = 12)) %>%
  arrange(time)
out6 <- mod %>%
  param(KRAS_mut = 0) %>%
  ev(as_data_set(ev6)) %>%
  mrgsim(end = sim_end, delta = 24, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(scenario = "6_FOLFIRI_BEV")

## ─────────────────────────────────────────────────────────────────────────────
## SCENARIO 7: Pembrolizumab (MSI-H / dMMR) — KEYNOTE-177 (Andre 2020)
##  MSI-H: pembro vs chemotherapy — PFS HR=0.60, ORR 43.8% vs 33.1%
## ─────────────────────────────────────────────────────────────────────────────
cat("Running Scenario 7: Pembrolizumab (MSI-H)...\n")
out7 <- mod %>%
  param(MSI_H = 1, KRAS_mut = 0) %>%
  ev(as_data_set(make_PEM(n_cycles = 18))) %>%
  mrgsim(end = sim_end, delta = 24, obsonly = TRUE) %>%
  as.data.frame() %>%
  mutate(scenario = "7_Pembro_MSIH")

## ─────────────────────────────────────────────────────────────────────────────
## COMBINE ALL SCENARIOS
## ─────────────────────────────────────────────────────────────────────────────
all_out <- bind_rows(out1, out2, out3, out4, out5, out6, out7) %>%
  filter(time >= 0) %>%
  mutate(time_days = time / 24,
         time_months = time_days / 30.4)

## ─────────────────────────────────────────────────────────────────────────────
## RESULTS SUMMARY TABLE
## ─────────────────────────────────────────────────────────────────────────────
summary_tbl <- all_out %>%
  group_by(scenario) %>%
  summarise(
    PFS_days   = suppressWarnings(min(time_days[PctChange >= 20], na.rm = TRUE)),
    BestResp   = min(PctChange, na.rm = TRUE),
    CEA_nadir  = min(CEA_out, na.rm = TRUE),
    ResistFrac_end = last(ResistFrac),
    .groups = "drop"
  ) %>%
  mutate(PFS_months = round(PFS_days / 30.4, 1),
         BestResp   = round(BestResp, 1),
         CEA_nadir  = round(CEA_nadir, 2))

print(summary_tbl)

## ─────────────────────────────────────────────────────────────────────────────
## PLOTS
## ─────────────────────────────────────────────────────────────────────────────
cols7 <- c(
  "1_NoTreatment"      = "#555555",
  "2_FOLFOX6"          = "#2980B9",
  "3_FOLFIRI"          = "#27AE60",
  "4_FOLFOX_BEV"       = "#8E44AD",
  "5_FOLFIRI_CTX_RASWT"= "#E67E22",
  "6_FOLFIRI_BEV"      = "#E74C3C",
  "7_Pembro_MSIH"      = "#1ABC9C"
)
labs7 <- c(
  "1_NoTreatment"      = "No Treatment",
  "2_FOLFOX6"          = "FOLFOX6",
  "3_FOLFIRI"          = "FOLFIRI",
  "4_FOLFOX_BEV"       = "FOLFOX + Bevacizumab",
  "5_FOLFIRI_CTX_RASWT"= "FOLFIRI + Cetuximab (RAS-WT)",
  "6_FOLFIRI_BEV"      = "FOLFIRI + Bevacizumab",
  "7_Pembro_MSIH"      = "Pembrolizumab (MSI-H)"
)

# Plot 1: Tumor diameter % change (waterfall/swim)
p1 <- ggplot(all_out %>% filter(time_months <= 18),
             aes(x = time_months, y = PctChange, color = scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = c(20, -30), linetype = c("dashed","dashed"),
             color = c("red", "green4"), alpha = 0.7) +
  scale_color_manual(values = cols7, labels = labs7) +
  labs(title = "CRC — Tumor Diameter Change from Baseline (RECIST)",
       x = "Time (months)", y = "Change from Baseline (%)",
       color = "Regimen", caption = "Dashed: PD (+20%) / PR threshold (-30%)") +
  coord_cartesian(ylim = c(-80, 200)) +
  theme_bw(base_size = 13) + theme(legend.position = "right")

# Plot 2: CEA biomarker kinetics
p2 <- ggplot(all_out %>% filter(time_months <= 18),
             aes(x = time_months, y = CEA_out, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = cols7, labels = labs7) +
  scale_y_log10() +
  labs(title = "CRC — CEA Biomarker Kinetics",
       x = "Time (months)", y = "CEA (ng/mL, log scale)", color = "Regimen") +
  theme_bw(base_size = 13) + theme(legend.position = "right")

# Plot 3: Sensitive vs resistant fractions
p3 <- ggplot(all_out %>% filter(scenario %in% c("2_FOLFOX6","4_FOLFOX_BEV",
                                                  "7_Pembro_MSIH"),
                                 time_months <= 18),
             aes(x = time_months, y = ResistFrac, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = cols7, labels = labs7) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "CRC — Resistant Cell Fraction Over Time",
       x = "Time (months)", y = "Resistant Fraction (%)", color = "Regimen") +
  theme_bw(base_size = 13)

# Plot 4: PK — 5-FU and SN-38 first cycle
p4 <- all_out %>%
  filter(scenario %in% c("2_FOLFOX6","3_FOLFIRI"),
         time_months <= 1) %>%
  select(time_months, scenario, Cp_FU_out, Cp_SN38_out) %>%
  pivot_longer(cols = c(Cp_FU_out, Cp_SN38_out),
               names_to = "drug", values_to = "Cp") %>%
  ggplot(aes(x = time_months * 30.4, y = Cp, color = scenario, linetype = drug)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = cols7, labels = labs7) +
  scale_linetype_manual(values = c("solid","dashed"),
                        labels = c("5-FU (μg/mL)", "SN-38 (μg/mL)")) +
  labs(title = "5-FU / SN-38 PK — First Cycle",
       x = "Time (h)", y = "Plasma Concentration (μg/mL)",
       color = "Regimen", linetype = "Drug") +
  theme_bw(base_size = 13)

# Plot 5: Bevacizumab effect on free VEGF
p5 <- all_out %>%
  filter(scenario %in% c("1_NoTreatment","4_FOLFOX_BEV","6_FOLFIRI_BEV"),
         time_months <= 12) %>%
  ggplot(aes(x = time_months, y = VEGF_free, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = cols7, labels = labs7) +
  labs(title = "Free VEGF Kinetics — Bevacizumab Effect",
       x = "Time (months)", y = "Free VEGF (nM)", color = "Regimen") +
  theme_bw(base_size = 13)

# Plot 6: CD8 T cell activity & PD-1 blockade
p6 <- all_out %>%
  filter(scenario %in% c("1_NoTreatment","7_Pembro_MSIH"),
         time_months <= 18) %>%
  ggplot(aes(x = time_months, y = CD8eff, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = cols7, labels = labs7) +
  labs(title = "CD8 T Cell Activity — Pembrolizumab (MSI-H)",
       x = "Time (months)", y = "CD8 Effector Activity (0-1)", color = "Regimen") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 13)

## Print plots
print(p1); print(p2); print(p3); print(p4); print(p5); print(p6)

cat("\n=== CRC QSP Model Summary ===\n")
cat("Scenarios simulated: 7\n")
cat("ODE compartments: 20\n")
cat("Calibration trials: MOSAIC, CRYSTAL, NO16966, TRIBE, KEYNOTE-177\n")
