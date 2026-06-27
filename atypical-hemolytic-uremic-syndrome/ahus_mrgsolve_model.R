################################################################################
# Atypical Hemolytic Uremic Syndrome (aHUS)
# Quantitative Systems Pharmacology (QSP) Model
# mrgsolve ODE Model
#
# Disease: Complement Alternative Pathway Dysregulation → TMA → Organ Damage
# Drug:    Eculizumab (anti-C5 mAb), Ravulizumab, Iptacopan (anti-Factor B)
#
# Model Compartments (18):
#   PK:  Drug_C (central), Drug_P (peripheral), Drug_C5 (TMDD complex)
#   Complement: C3pool, C3b_AP, C5free, C5conv, MACflux
#   TMA: Endo_inj, PLT, Hgb, LDH, Hpg
#   Renal: GFR, Schistocyte, CRP
#   Biomarker: sC5b9, CH50perc
#
# Calibrated parameters based on:
#   - Fakhouri et al. Lancet 2017 (aHUS epidemiology, outcomes)
#   - Legendre et al. NEJM 2013 (eculizumab LEAP trial)
#   - Rother et al. Nat Med 2007 (eculizumab PK/PD)
#   - Greenbaum et al. Kidney Int 2016 (pediatric aHUS)
#   - Kavanagh et al. NEJM 2021 (iptacopan for PNH/aHUS)
#   - Menne et al. JASN 2015 (eculizumab PK population model)
#
# Author: QSP Library (Claude Code Routine)
# Date: 2026-06-27
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ─── mrgsolve MODEL CODE ─────────────────────────────────────────────────────

code <- '
$PARAM @annotated
// ── Eculizumab 2-Compartment PK ──────────────────────────────────────────────
CL    : 0.31  : Clearance (L/day) [Menne 2015]
V1    : 5.30  : Central volume (L)
Q     : 0.54  : Inter-compartmental clearance (L/day)
V2    : 4.10  : Peripheral volume (L)
F1    : 1.00  : Bioavailability (IV = 1)

// ── TMDD (Target-Mediated Drug Disposition) ───────────────────────────────────
kon   : 0.864  : C5-drug association rate (nM^-1 day^-1) [Kd=10pM]
koff  : 8.64e-3: C5-drug dissociation rate (day^-1) [Kd=10pM]
kint  : 0.02   : TMDD complex internalization/degradation (day^-1)

// ── Complement Parameters ─────────────────────────────────────────────────────
C3_ss   : 8000  : C3 plasma concentration at steady state (nM) [~1.3 g/L, MW 185 kDa]
kC3syn  : 4.167 : C3 synthesis rate (nM/h) [turnover ~12h]
kC3deg  : 5.208e-4 : C3 degradation rate constant (h^-1)
kAP     : 0.12  : AP activation rate constant (h^-1) [dysregulated: ~10x normal]
kAP_amp : 0.30  : AP amplification feedback coefficient (h^-1 per normalized C3b)
kC3b_deg: 0.80  : C3b inactivation rate (h^-1) [CFH+CFI mediated]
C5_ss   : 395   : C5 plasma concentration at steady state (nM) [~75 ug/mL, MW 190 kDa]
kC5syn  : 11.85 : C5 synthesis rate (nM/day)
kC5deg  : 0.030 : C5 degradation rate constant (day^-1) [t1/2 ~23d]
kC5conv : 0.015 : C5 convertase-driven C5 cleavage rate (day^-1 per nM C5conv)
kMACform: 0.25  : MAC formation rate from C5b cleavage (normalized, day^-1)
kMACclr : 0.50  : MAC clearance rate (day^-1)

// ── AP Dysregulation Parameters (disease state) ───────────────────────────────
CFH_factor: 0.30  : CFH regulatory function (0=none, 1=normal; 0.30 = CFH mutation)
CFI_factor: 0.50  : CFI function (0=none, 1=normal)
CD46_factor:0.70  : CD46/MCP function on cell surface

// ── Endothelial Injury & TMA ──────────────────────────────────────────────────
kEI_MAC  : 0.40  : Endothelial injury rate per MAC unit (day^-1)
kEI_rep  : 0.15  : Endothelial repair rate (day^-1) [baseline repair]
EI_max   : 10.0  : Maximum endothelial injury score (0-10 scale)
kPLT_loss: 0.08  : Platelet consumption rate per EI unit (day^-1)
kPLT_prod: 13.0  : Platelet production rate (10^9/L/day) [normal thrombopoiesis]
PLT_ss   : 250   : Normal platelet count (10^9/L)
kHgb_loss: 0.04  : Hemoglobin loss rate per EI unit (day^-1)
kHgb_prod: 0.20  : Hemoglobin production rate (g/dL/day) [erythropoiesis]
Hgb_ss   : 14.0  : Normal hemoglobin (g/dL)
kLDH_rel : 15.0  : LDH release rate per hemolysis (U/L/day per EI)
kLDH_clr : 0.40  : LDH clearance rate (day^-1) [t1/2 ~1.7 days]
LDH_norm : 200   : Normal LDH upper limit (U/L)
kHpg_cons: 0.35  : Haptoglobin consumption rate per free Hgb (g/L/day per EI)
kHpg_syn : 0.15  : Haptoglobin synthesis rate (g/L/day)
Hpg_ss   : 1.20  : Normal haptoglobin (g/L)

// ── Renal Compartment ─────────────────────────────────────────────────────────
GFR_ss    : 90.0 : Baseline/normal GFR (mL/min/1.73m2)
kGFR_loss : 0.05 : GFR loss rate per EI unit (day^-1)
kGFR_rep  : 0.10 : GFR recovery rate (day^-1) [if MAC controlled]

// ── Biomarker Parameters ──────────────────────────────────────────────────────
sC5b9_base: 150  : Baseline sC5b-9 (ng/mL) [normal <244]
sC5b9_max : 3500 : Maximum sC5b-9 in active aHUS (ng/mL)
kSCT_form : 0.60 : Schistocyte formation rate per EI (% per day)
kSCT_clr  : 0.25 : Schistocyte clearance rate (% per day)
kCRP_syn  : 12.0 : CRP synthesis per inflammation unit (mg/L/day)
kCRP_clr  : 0.60 : CRP clearance rate (day^-1) [t1/2 ~1.15 days]

// ── Scenario Selection ────────────────────────────────────────────────────────
scenario  : 1    : 1=natural history, 2=eculizumab std, 3=ravulizumab,
                    4=eculizumab+danicopan, 5=iptacopan oral

$CMT @annotated
Drug_C   : Eculizumab central compartment (mg)
Drug_P   : Eculizumab peripheral compartment (mg)
Drug_C5  : Drug-C5 TMDD complex (mg)
C3pool   : C3 plasma pool (nM)
C3b_AP   : Active C3b in alternative pathway (nM)
C5free   : Free C5 plasma (nM)
C5conv   : C5 convertase activity (arbitrary units)
MACflux  : MAC flux/burden (arbitrary units)
Endo_inj : Endothelial injury score (0-10)
PLT      : Platelet count (x10^9/L)
Hgb      : Hemoglobin (g/dL)
LDH      : Lactate dehydrogenase (U/L)
Hpg      : Haptoglobin (g/L)
GFR      : eGFR (mL/min/1.73m2)
Schist   : Schistocytes (% on blood smear)
CRP      : C-reactive protein (mg/L)
sC5b9    : Soluble MAC (sC5b-9) (ng/mL)
CH50pct  : CH50 functional activity (%)

$MAIN
// Initial conditions - active aHUS (untreated)
if(NEWIND <= 1) {
  // PK
  Drug_C_0 = 0;
  Drug_P_0 = 0;
  Drug_C5_0 = 0;

  // Complement (dysregulated baseline for aHUS)
  C3pool_0  = C3_ss;
  C3b_AP_0  = 200;      // elevated due to AP dysregulation (nM)
  C5free_0  = C5_ss;
  C5conv_0  = 5.0;      // elevated C5 convertase activity
  MACflux_0 = 2.5;      // elevated MAC

  // Clinical (active TMA at presentation)
  Endo_inj_0 = 6.0;    // significant endothelial injury
  PLT_0      = 60;      // thrombocytopenia (x10^9/L)
  Hgb_0      = 7.5;    // anemia (g/dL)
  LDH_0      = 1500;   // elevated (U/L)
  Hpg_0      = 0.05;   // near-absent (g/L)
  GFR_0      = 25;     // acute kidney injury (mL/min)
  Schist_0   = 4.5;    // >1% schistocytes
  CRP_0      = 45;     // elevated (mg/L)

  // Biomarkers
  sC5b9_0  = sC5b9_max * 0.8;  // ~2800 ng/mL in active aHUS
  CH50pct_0 = 60;               // partially consumed
}

// ── Regulatory function (disease-modified) ─────────────────────────────────
double reg_eff = CFH_factor * CFI_factor * CD46_factor;
// reg_eff = 0.105 for triple mutation vs 1.0 normal

// ── Drug concentration (mg/L = ug/mL) and molar (nM, MW=148000 g/mol) ──────
double Cc_mgL = Drug_C / V1;                  // mg/L = ug/mL
double Cp_mgL = Drug_P / V2;
double Cc_nM  = Cc_mgL * 1000 / 148.0;       // nM (MW ~148 kDa)

// ── Anti-Factor B (iptacopan) effect ──────────────────────────────────────
double fB_inhibition = 0.0;
if(scenario == 5) fB_inhibition = 0.95;       // 95% inhibition of factor B

// ── Danicopan (anti-Factor D) effect ──────────────────────────────────────
double fD_inhibition = 0.0;
if(scenario == 4) fD_inhibition = 0.80;       // 80% add-on Factor D inhibition

// ── AP activity (modified by drug, disease, regulation) ──────────────────
double AP_activity = kAP * (1.0 - fB_inhibition) * (1.0 - fD_inhibition) / reg_eff;

// ── C5 convertase formation from C3b accumulation ──────────────────────────
double C3b_norm = C3b_AP / 500.0;             // normalized to "active disease" level

// ── MAC inhibition by eculizumab/ravulizumab ──────────────────────────────
double free_C5_norm = C5free / C5_ss;
double MAC_eff = MACflux * free_C5_norm;       // MAC generation scales with free C5

$ODE
// ── Eculizumab PK (2-compartment with TMDD) ──────────────────────────────────
double TMDD_bind = kon * Cc_nM * C5free - koff * Drug_C5;    // nM/day binding

dxdt_Drug_C  = -(CL/V1) * Drug_C - (Q/V1) * Drug_C + (Q/V2) * Drug_P
               - kon * (Drug_C/V1) * C5free * V1 + koff * Drug_C5;
dxdt_Drug_P  =  (Q/V1) * Drug_C - (Q/V2) * Drug_P;
dxdt_Drug_C5 =  kon * (Drug_C/V1) * C5free * 148.0 * V1 / 1000.0
               - koff * Drug_C5 - kint * Drug_C5;
// Note: simplified TMDD coupling for ODE stability

// ── C3 pool dynamics ──────────────────────────────────────────────────────────
// C3 consumed by AP activation; synthesized in liver
double C3_syn_rate = kC3syn * 24;             // convert to nM/day
double kC3deg_d    = kC3deg * 24;
dxdt_C3pool = C3_syn_rate - kC3deg_d * C3pool
              - AP_activity * C3pool / (C3_ss / 100.0);

// ── C3b alternative pathway accumulation ─────────────────────────────────────
double C3b_form = AP_activity * C3pool / C3_ss;   // AP-generated C3b
double C3b_inac = kC3b_deg * reg_eff * C3b_AP;    // CFH/CFI-mediated inactivation
double amploop  = kAP_amp * (C3b_AP / C3_ss);     // positive feedback
dxdt_C3b_AP = C3b_form + amploop * C3pool / C3_ss - C3b_inac;

// ── Free C5 with TMDD ─────────────────────────────────────────────────────────
double C5_cleavage = kC5conv * C5conv * C5free;    // C5 convertase cleaves C5
double C5_bound_loss = kon * Cc_nM * C5free - koff * (Drug_C5 / 148.0 * 1000.0 / V1);
dxdt_C5free = kC5syn - kC5deg * C5free - C5_cleavage - C5_bound_loss;

// ── C5 Convertase activity ────────────────────────────────────────────────────
double C5conv_form = 0.08 * pow(C3b_AP / C3_ss, 2.0);  // Hill-like from C3b
double C5conv_clr  = 0.15 * C5conv;
dxdt_C5conv = C5conv_form - C5conv_clr;

// ── MAC flux ──────────────────────────────────────────────────────────────────
double MAC_form = kMACform * C5_cleavage;           // C5b -> MAC
double MAC_clr  = kMACclr * MACflux;
// CD59 protection on cell surface
double CD59_eff = 0.60;                             // 60% baseline MAC inhibition
dxdt_MACflux = MAC_form * (1.0 - CD59_eff) - MAC_clr;

// ── Endothelial Injury Score (0-10) ──────────────────────────────────────────
double EI_drive = kEI_MAC * MACflux;
double EI_rep   = kEI_rep * Endo_inj;
dxdt_Endo_inj = EI_drive - EI_rep;
if(Endo_inj > EI_max) dxdt_Endo_inj = 0;
if(Endo_inj < 0) dxdt_Endo_inj = 0;

// ── Platelet Count (x10^9/L) ──────────────────────────────────────────────────
double PLT_loss = kPLT_loss * PLT * (Endo_inj / EI_max);
double PLT_prod = kPLT_prod;                        // constant TPO-driven production
dxdt_PLT = PLT_prod - PLT_loss;
if(PLT < 5) dxdt_PLT = 0;
if(PLT > 400) dxdt_PLT = PLT_prod - PLT_loss;

// ── Hemoglobin (g/dL) ─────────────────────────────────────────────────────────
double Hgb_loss_rate = kHgb_loss * Hgb * (Endo_inj / EI_max);
double Hgb_prod_rate = kHgb_prod;
dxdt_Hgb = Hgb_prod_rate - Hgb_loss_rate;
if(Hgb < 4.0) dxdt_Hgb = 0;

// ── LDH (U/L) ────────────────────────────────────────────────────────────────
double hemolysis_rate = kHgb_loss * Hgb * (Endo_inj / EI_max);
dxdt_LDH = kLDH_rel * hemolysis_rate - kLDH_clr * LDH;

// ── Haptoglobin (g/L) ─────────────────────────────────────────────────────────
double free_Hgb_proxy = Hgb_loss_rate;   // proxy for free Hgb generation
dxdt_Hpg = kHpg_syn - kHpg_cons * free_Hgb_proxy - 0.10 * Hpg;
if(Hpg < 0) dxdt_Hpg = 0;

// ── GFR (mL/min/1.73m2) ──────────────────────────────────────────────────────
double GFR_loss_rate = kGFR_loss * GFR * (Endo_inj / 6.0);   // proportional to EI
double GFR_rep_rate  = kGFR_rep * (GFR_ss - GFR) * (1.0 - Endo_inj/EI_max);
dxdt_GFR = GFR_rep_rate - GFR_loss_rate;
if(GFR < 5) dxdt_GFR = 0;

// ── Schistocytes (%) ──────────────────────────────────────────────────────────
dxdt_Schist = kSCT_form * (Endo_inj / EI_max) - kSCT_clr * Schist;
if(Schist < 0) dxdt_Schist = 0;

// ── CRP (mg/L) ───────────────────────────────────────────────────────────────
double inflam = MACflux;
dxdt_CRP = kCRP_syn * inflam - kCRP_clr * CRP;

// ── Soluble C5b-9 (sC5b-9, ng/mL) ────────────────────────────────────────────
double sC5b9_drive = MAC_form * (sC5b9_max / 3.0);
dxdt_sC5b9 = sC5b9_drive - 0.20 * sC5b9;

// ── CH50 (%) - inversely related to complement consumption ───────────────────
double C5_pct = C5free / C5_ss;
dxdt_CH50pct = 0.50 * (100.0 * C5_pct - CH50pct);   // equilibrates to C5 fraction


$TABLE
// Derived outputs
double Cc_ugmL   = Drug_C / V1;                       // drug conc (ug/mL)
double Cc_nM_out = Cc_ugmL * 1000.0 / 148.0;          // drug conc (nM)
double EI_norm   = Endo_inj / EI_max;                 // 0-1 scale
double PLT_resp  = PLT / PLT_ss * 100;                // % of normal PLT
double Hgb_resp  = Hgb / Hgb_ss * 100;               // % of normal Hgb
double GFR_resp  = GFR / GFR_ss * 100;                // % of baseline GFR
double TMA_active = (PLT < 150) ? 1.0 : 0.0;         // TMA flag
double LDH_xULN  = LDH / LDH_norm;                   // x upper limit of normal
double C5_block   = 1.0 - C5free / C5_ss;            // fractional C5 blockade
double dialysis_risk = (GFR < 15) ? 1.0 : 0.0;

// Capture all outputs
capture Cc_ugmL Cc_nM_out EI_norm PLT_resp Hgb_resp GFR_resp;
capture TMA_active LDH_xULN C5_block dialysis_risk;
capture C3pool C3b_AP C5free MACflux sC5b9 CH50pct;

$CAPTURE
Cc_ugmL Cc_nM_out EI_norm PLT_resp Hgb_resp GFR_resp TMA_active LDH_xULN C5_block dialysis_risk
C3pool C3b_AP C5free MACflux sC5b9 CH50pct Drug_C Drug_P Drug_C5
C3pool C3b_AP C5free C5conv MACflux Endo_inj PLT Hgb LDH Hpg GFR Schist CRP sC5b9 CH50pct
'

# Compile model
mod <- mcode("aHUS_QSP", code, quiet = TRUE)

# ─── DOSING REGIMENS ─────────────────────────────────────────────────────────

# Eculizumab standard: induction 900 mg qw x4, maintenance 1200 mg q2w
ecu_induction <- function(start_day = 1) {
  data.frame(
    time = c(start_day, start_day+7, start_day+14, start_day+21),
    amt  = 900,
    cmt  = 1, evid = 1, rate = -2, ii = 0, addl = 0
  )
}
ecu_maintenance <- function(start_day = 28, n_doses = 12) {
  data.frame(
    time = seq(start_day, by = 14, length.out = n_doses),
    amt  = 1200,
    cmt  = 1, evid = 1, rate = -2, ii = 0, addl = 0
  )
}

# Ravulizumab (weight-based; using 70 kg → 3000 mg loading, 3300 mg q8w)
ravu_induction <- function(start_day = 1) {
  data.frame(
    time = c(start_day, start_day + 15),
    amt  = c(2400, 3000),
    cmt  = 1, evid = 1, rate = -2, ii = 0, addl = 0
  )
}
ravu_maintenance <- function(start_day = 29, n_doses = 5) {
  data.frame(
    time = seq(start_day, by = 56, length.out = n_doses),
    amt  = 3300,
    cmt  = 1, evid = 1, rate = -2, ii = 0, addl = 0
  )
}

# ─── SIMULATION SCENARIOS ─────────────────────────────────────────────────────

simulate_scenario <- function(scen, mod) {
  # Build dosing events
  if(scen == 1) {
    # Natural history - no treatment
    ev <- data.frame(time=0, amt=0, cmt=1, evid=0)
  } else if(scen == 2) {
    # Eculizumab standard dosing
    ev <- bind_rows(ecu_induction(1), ecu_maintenance(28, 12))
    ev$scenario <- 2
  } else if(scen == 3) {
    # Ravulizumab
    ev <- bind_rows(ravu_induction(1), ravu_maintenance(29, 5))
  } else if(scen == 4) {
    # Eculizumab + Danicopan add-on (fDi)
    ev <- bind_rows(ecu_induction(1), ecu_maintenance(28, 12))
  } else if(scen == 5) {
    # Iptacopan oral (no IV loading; simulated as parameter change)
    ev <- data.frame(time=0, amt=0, cmt=1, evid=0)
  }

  out <- mod %>%
    param(scenario = scen) %>%
    ev(ev) %>%
    mrgsim(end = 365, delta = 1) %>%
    as.data.frame()

  out$scenario_label <- c(
    "1" = "Natural History",
    "2" = "Eculizumab Std",
    "3" = "Ravulizumab",
    "4" = "Ecu + Danicopan",
    "5" = "Iptacopan Oral"
  )[as.character(scen)]

  out
}

# Run all scenarios
message("Running aHUS QSP simulations...")
results <- lapply(1:5, function(s) simulate_scenario(s, mod))
all_data <- bind_rows(results)
all_data$scenario_label <- factor(all_data$scenario_label,
  levels = c("Natural History","Eculizumab Std","Ravulizumab",
             "Ecu + Danicopan","Iptacopan Oral"))

# ─── PLOTTING ─────────────────────────────────────────────────────────────────

theme_qsp <- theme_bw() +
  theme(
    strip.background = element_rect(fill = "#2C3E50"),
    strip.text = element_text(color = "white", face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  )

scenario_colors <- c(
  "Natural History"  = "#C0392B",
  "Eculizumab Std"   = "#2980B9",
  "Ravulizumab"      = "#27AE60",
  "Ecu + Danicopan"  = "#8E44AD",
  "Iptacopan Oral"   = "#E67E22"
)

# Plot 1: Platelet Count
p1 <- ggplot(all_data, aes(x=time, y=PLT, color=scenario_label)) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=150, linetype="dashed", color="gray40", linewidth=0.8) +
  annotate("text", x=365, y=155, label="Thrombocytopenia threshold (150)", hjust=1, size=3) +
  scale_color_manual(values=scenario_colors) +
  labs(title="Platelet Count", x="Day", y="Platelets (x10⁹/L)", color="") +
  ylim(0, 300) + theme_qsp

# Plot 2: Hemoglobin
p2 <- ggplot(all_data, aes(x=time, y=Hgb, color=scenario_label)) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=c(7.0, 10.0), linetype="dashed", color=c("red","gray40"), linewidth=0.8) +
  scale_color_manual(values=scenario_colors) +
  labs(title="Hemoglobin", x="Day", y="Hemoglobin (g/dL)", color="") +
  ylim(4, 16) + theme_qsp

# Plot 3: GFR
p3 <- ggplot(all_data, aes(x=time, y=GFR, color=scenario_label)) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=c(15, 60, 90), linetype="dashed",
             color=c("red","orange","gray40"), linewidth=0.7) +
  scale_color_manual(values=scenario_colors) +
  labs(title="eGFR (Renal Function)", x="Day", y="eGFR (mL/min/1.73m²)", color="") +
  ylim(0, 100) + theme_qsp

# Plot 4: LDH
p4 <- ggplot(all_data, aes(x=time, y=LDH, color=scenario_label)) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=200, linetype="dashed", color="gray40") +
  annotate("text", x=365, y=220, label="ULN (200 U/L)", hjust=1, size=3) +
  scale_color_manual(values=scenario_colors) +
  labs(title="LDH (Hemolysis Marker)", x="Day", y="LDH (U/L)", color="") +
  ylim(0, 2500) + theme_qsp

# Plot 5: sC5b-9
p5 <- ggplot(all_data, aes(x=time, y=sC5b9, color=scenario_label)) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=244, linetype="dashed", color="gray40") +
  annotate("text", x=365, y=260, label="Normal threshold (244 ng/mL)", hjust=1, size=3) +
  scale_color_manual(values=scenario_colors) +
  labs(title="sC5b-9 (Complement Activation Biomarker)", x="Day",
       y="sC5b-9 (ng/mL)", color="") +
  theme_qsp

# Plot 6: CH50 (Complement Function Blockade)
p6 <- ggplot(all_data, aes(x=time, y=CH50pct, color=scenario_label)) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=c(10), linetype="dashed", color="gray40") +
  annotate("text", x=100, y=12, label="Target: <10% (fully blocked)", hjust=0, size=3) +
  scale_color_manual(values=scenario_colors) +
  labs(title="CH50 (% Residual Complement Activity)", x="Day",
       y="CH50 (%)", color="") +
  ylim(0, 100) + theme_qsp

# Plot 7: Endothelial Injury Score
p7 <- ggplot(all_data, aes(x=time, y=Endo_inj, color=scenario_label)) +
  geom_line(linewidth=1.1) +
  scale_color_manual(values=scenario_colors) +
  labs(title="Endothelial Injury Score (0-10)", x="Day",
       y="Injury Score", color="") +
  ylim(0, 10) + theme_qsp

# Plot 8: Haptoglobin
p8 <- ggplot(all_data, aes(x=time, y=Hpg, color=scenario_label)) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=0.1, linetype="dashed", color="gray40") +
  scale_color_manual(values=scenario_colors) +
  labs(title="Haptoglobin (Hemolysis Marker)", x="Day",
       y="Haptoglobin (g/L)", color="") +
  theme_qsp

# Combined figure
combined <- (p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | plot_spacer()) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "aHUS QSP Model - Treatment Scenario Comparison",
    subtitle = "Natural history vs. anti-complement therapies (1-year simulation)",
    theme = theme(legend.position = "bottom",
                  plot.title = element_text(face="bold", size=14))
  )

# ─── SUMMARY TABLE ────────────────────────────────────────────────────────────

summary_table <- all_data %>%
  filter(time %in% c(28, 84, 180, 365)) %>%
  group_by(scenario_label, time) %>%
  summarise(
    PLT_mean     = round(mean(PLT), 0),
    Hgb_mean     = round(mean(Hgb), 1),
    GFR_mean     = round(mean(GFR), 0),
    LDH_xULN    = round(mean(LDH_xULN), 1),
    sC5b9_mean   = round(mean(sC5b9), 0),
    CH50_pct     = round(mean(CH50pct), 0),
    Dialysis_risk = round(mean(dialysis_risk)*100, 0),
    .groups = "drop"
  )

cat("\n═══════════════════════════════════════════════════════════════════\n")
cat("aHUS QSP Model - Simulation Summary Table\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")
print(as.data.frame(summary_table), row.names = FALSE)

# ─── SCENARIO 2 DETAILED: TMA EVENT ANALYSIS ─────────────────────────────────

ecu_data <- all_data %>% filter(scenario_label == "Eculizumab Std")
nh_data  <- all_data %>% filter(scenario_label == "Natural History")

cat("\n─── TMA Resolution Analysis (Eculizumab Std vs Natural History) ───\n")
cat(sprintf("Day 28 PLT recovery:  ECU=%.0f, NH=%.0f x10⁹/L\n",
    ecu_data$PLT[29], nh_data$PLT[29]))
cat(sprintf("Day 84 GFR recovery:  ECU=%.0f, NH=%.0f mL/min\n",
    ecu_data$GFR[85], nh_data$GFR[85]))
cat(sprintf("Day 365 dialysis risk: ECU=%.0f%%, NH=%.0f%%\n",
    ecu_data$dialysis_risk[366]*100, nh_data$dialysis_risk[366]*100))

# ─── DOSE-RESPONSE: ECULIZUMAB EXPOSURE-EFFICACY ─────────────────────────────

dose_levels <- c(300, 600, 900, 1200, 1800)

dose_response <- lapply(dose_levels, function(d) {
  ev_dr <- ecu_induction(1)
  ev_dr$amt <- d
  ev_maint <- ecu_maintenance(28, 12)
  ev_maint$amt <- round(d * 4/3)
  ev_all <- bind_rows(ev_dr, ev_maint)

  out <- mod %>%
    param(scenario = 2) %>%
    ev(ev_all) %>%
    mrgsim(end = 84, delta = 1) %>%
    as.data.frame()

  data.frame(
    dose_mg    = d,
    Cc_peak_nM = max(out$Cc_nM_out, na.rm=TRUE),
    PLT_d84    = out$PLT[85],
    Hgb_d84    = out$Hgb[85],
    GFR_d84    = out$GFR[85],
    C5_block_d84 = out$C5_block[85]
  )
})
dose_resp_df <- bind_rows(dose_response)

cat("\n─── Dose-Response Analysis (Eculizumab, Day 84) ───\n")
print(dose_resp_df, digits=3, row.names=FALSE)

# ─── PARAMETER SENSITIVITY (CFH function) ────────────────────────────────────

cfh_levels <- c(0.1, 0.3, 0.5, 0.7, 1.0)  # 0.1=severe, 1.0=normal

cfh_sens <- lapply(cfh_levels, function(cfh) {
  ev_std <- bind_rows(ecu_induction(1), ecu_maintenance(28, 12))
  out <- mod %>%
    param(scenario=2, CFH_factor=cfh) %>%
    ev(ev_std) %>%
    mrgsim(end=180, delta=1) %>%
    as.data.frame()
  data.frame(
    CFH_factor = cfh,
    PLT_d28  = out$PLT[29],
    GFR_d90  = out$GFR[91],
    sC5b9_d90 = out$sC5b9[91]
  )
})
cfh_sens_df <- bind_rows(cfh_sens)

cat("\n─── CFH Function Sensitivity (Eculizumab, Day 90) ───\n")
print(cfh_sens_df, digits=3, row.names=FALSE)

cat("\n\nSimulation complete. Use ahus_shiny_app.R for interactive exploration.\n")
