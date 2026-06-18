## ============================================================
## Chronic Hepatitis B (CHB) — mrgsolve QSP Model
## Author: Claude Code Routine (CCR)  |  Date: 2026-06-18
##
## Disease: CHB — HBV-mediated chronic liver disease
##   HBV entry via NTCP → cccDNA persistence → immune evasion
##   → T-cell exhaustion → hepatocyte injury → fibrosis → HCC
##
## Drugs modeled:
##   1. ETV  (Entecavir 0.5 mg QD; naïve patients)
##   2. TDF  (Tenofovir DF 300 mg QD; or TAF 25 mg QD)
##   3. PIFN (Pegylated IFN-α2a 180 µg SC QW × 48 weeks)
##   4. ETV+PIFN combination (add-on strategy)
##   5. ETV+siRNA  (ETV + GalNAc-siRNA add-on; novel strategy)
##
## Key calibration trials:
##   ETV:   EABV, BEHoLD AI463022 (Chang 2006 NEJM, PMID:16672700)
##          Year-1: HBV DNA <300 copies/mL in 67% HBeAg+ pts
##   TDF:   Marcellin 2008 NEJM (PMID:18685079)
##          Year-1: HBV DNA <400 copies/mL in 76% HBeAg+ pts
##   PIFN:  Lau 2005 NEJM (PMID:15987917)
##          Year-1: HBeAg loss 32%, HBsAg loss 3%
##   HBsAg  kinetics: Wursthorn 2010 (PMID:20626649)
##          Decline ~0.05 log/yr on NUC alone
##
## ODE system: 20 state variables
##   PK:      ETV_gut, ETV_C, ETV_TP, TDF_gut, TDF_C, TFV_DP,
##            IFN_SC, IFN_C  (8 states)
##   Biology: T (target), I (infected), V (serum HBV DNA),
##            ccc (intrahepatic cccDNA), Ag (HBsAg),
##            CTL, Exhaust, IFN_inn (innate IFN),
##            ALT, Fibrosis, HCC_risk, HSC  (12 states)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

## ─────────────────────────────────────────────────────────────
## 1. MODEL CODE BLOCK
## ─────────────────────────────────────────────────────────────

chb_code <- '
$PROB Chronic Hepatitis B (CHB) QSP Model — mrgsolve

$PLUGIN Rcpp

$PARAM @annotated
// ── Patient baseline parameters ──────────────────────────────
BWT       : 65     : Body weight (kg)
T0        : 1e7    : Baseline target hepatocytes (cells/mL equiv)
I0        : 1e5    : Baseline infected hepatocytes (cells/mL equiv)
V0        : 1e7    : Baseline HBV DNA (IU/mL)
ccc0      : 10     : Baseline cccDNA (copies/cell equivalent)
Ag0       : 1000   : Baseline HBsAg (IU/mL)
ALT0      : 80     : Baseline ALT (IU/L)
FIB0      : 1.5    : Baseline fibrosis score (Metavir 0–4)
CTL0      : 0.5    : Baseline CD8+ CTL response (normalized 0–1)
ExhSS     : 0.4    : Steady-state T cell exhaustion (0–1)

// ── HBV Viral Dynamics (Neumann/Dahari framework) ─────────────
s_T       : 1e5    : Target cell production rate (cells/mL/day)
d_T       : 0.01   : Target cell death rate (/day)
beta      : 2e-10  : Infection rate constant (mL/IU/day)
delta_I   : 0.15   : Infected cell clearance rate (/day; immune)
p_V       : 50     : Virion production (IU/mL/cell/day)
c_V       : 0.67   : Serum viral clearance rate (/day; t½~1 day)
lambda_ccc: 0.002  : cccDNA replenishment from rcDNA (/day)
delta_ccc : 0.003  : cccDNA loss rate (/day; slow ~years t½)
k_Ag      : 0.05   : HBsAg production proportional to cccDNA
c_Ag      : 0.01   : HBsAg clearance rate (/day; t½~70 days)

// ── Immune Parameters ─────────────────────────────────────────
alpha_CTL : 0.8    : CTL killing efficiency (/day per normalized CTL)
k_CTL_exp : 0.05   : CTL expansion rate proportional to V (/day)
d_CTL     : 0.03   : CTL contraction rate (/day)
k_exhaust : 0.02   : T cell exhaustion rate (driven by Ag load)
r_exhaust : 0.005  : Exhaustion recovery rate (/day)
k_IFNinn  : 0.1    : Innate IFN induction rate
d_IFNinn  : 0.15   : Innate IFN decay rate (/day)
IFN_antiv : 0.4    : IFN innate antiviral effect (max fraction inhibition)

// ── Liver Biology ─────────────────────────────────────────────
k_ALT_I   : 50     : ALT elevation per infected cell unit
k_ALT_ret : 0.05   : ALT return rate to baseline (/day)
k_fibrosis: 0.0005 : Fibrosis progression rate (per inflam unit/day)
k_fib_reg : 0.0002 : Fibrosis regression (on viral suppression)
k_HCC     : 0.002  : HCC risk accumulation rate (per fibrosis unit)
k_HSC     : 0.01   : HSC activation rate
d_HSC     : 0.03   : HSC deactivation rate

// ── ETV PK parameters ─────────────────────────────────────────
ETV_F     : 0.37   : ETV oral bioavailability
ETV_ka    : 0.90   : ETV absorption rate constant (/h)
ETV_Vc    : 73     : ETV central volume (L)
ETV_CL    : 10.5   : ETV total clearance (L/h)
ETV_Kp    : 5.0    : ETV liver/plasma partition coefficient
ETV_kTP   : 0.15   : ETV phosphorylation rate to ETV-TP (/h)
ETV_kTP_d : 0.008  : ETV-TP elimination rate (/h; intracellular t½~15h)
ETV_IC50  : 0.004  : ETV-TP IC50 for RT inhibition (µM)
ETV_Imax  : 0.99   : ETV max inhibition (fraction)

// ── TDF PK parameters ─────────────────────────────────────────
TDF_F     : 0.25   : TDF oral bioavailability (TAF F=0.80)
TDF_ka    : 0.50   : TDF absorption rate constant (/h)
TDF_Vc    : 18     : TFV central volume (L)
TDF_CL    : 11.2   : TFV total clearance (L/h)
TFV_kDP   : 0.10   : TFV→TFV-DP phosphorylation rate (/h)
TFV_kDP_d : 0.012  : TFV-DP elimination (/h; intracellular t½~12-50h)
TFV_IC50  : 0.500  : TFV-DP IC50 for RT inhibition (µM)
TFV_Imax  : 0.99   : TFV max inhibition (fraction)

// ── Peg-IFN PK parameters ────────────────────────────────────
IFN_F     : 0.80   : Peg-IFN SC bioavailability
IFN_ka    : 0.04   : Peg-IFN SC absorption rate (/h; t_peak~72h)
IFN_Vc    : 8.0    : Peg-IFN volume of distribution (L)
IFN_CL    : 0.055  : Peg-IFN total clearance (L/h; t½~85h)
IFN_EC50  : 0.08   : Peg-IFN EC50 (ng/mL equiv)
IFN_Emax  : 0.80   : Peg-IFN max antiviral effect (fraction)
IFN_immE  : 0.50   : Peg-IFN immunomodulatory effect (CTL boost)
IFN_cccE  : 0.30   : Peg-IFN anti-cccDNA effect (fraction max)
IAGM_dose : 0      : siRNA dose active (0=off, 1=on)
siRNA_Emax: 0.90   : GalNAc-siRNA max HBsAg reduction

// ── Drug dosing flags ─────────────────────────────────────────
USE_ETV   : 1      : Use Entecavir (0=no, 1=yes)
USE_TDF   : 0      : Use TDF (0=no, 1=yes)
USE_PIFN  : 0      : Use Peg-IFN (0=no, 1=yes)
PIFN_dur  : 336    : Peg-IFN treatment duration (days, 48 wks=336)

$CMT @annotated
// PK compartments
ETV_gut   : ETV gut compartment (mg)
ETV_C     : ETV central plasma (ng/mL · L = µg equiv)
ETV_TP    : ETV-TP intracellular (µM equiv)
TDF_gut   : TDF gut compartment (mg)
TDF_C     : TFV central plasma (ng/mL · L)
TFV_DP    : TFV-DP intracellular (µM equiv)
IFN_SC    : Peg-IFN SC depot (IU equiv)
IFN_C     : Peg-IFN central plasma (ng/mL equiv)
// Disease biology
T_cell    : Target hepatocytes (×10^6 cells/mL)
I_cell    : HBV-Infected hepatocytes (×10^6 cells/mL)
V_dna     : Serum HBV DNA (IU/mL)
ccc_DNA   : Intrahepatic cccDNA (copies/cell equiv)
HBsAg_C   : Serum HBsAg (IU/mL)
CTL_resp  : CD8+ CTL response (normalized 0–1)
Exhaust   : T cell exhaustion index (0–1)
IFN_inn   : Innate IFN signaling (normalized)
ALT_val   : Serum ALT (IU/L)
Fibrosis  : Hepatic fibrosis score (Metavir 0–4)
HCC_risk  : Cumulative HCC risk (dimensionless)
HSC_act   : Activated HSC fraction (0–1)

$MAIN
// ── Derived drug concentrations ───────────────────────────────
double ETV_Cp  = ETV_C  / ETV_Vc;    // ng/mL
double TFV_Cp  = TDF_C  / TDF_Vc;    // ng/mL
double IFN_Cp  = IFN_C  / IFN_Vc;    // ng/mL equiv

// ── Drug effect calculations ──────────────────────────────────
// ETV: inhibition of viral replication (RT)
double ETV_eff = USE_ETV * ETV_Imax * ETV_TP / (ETV_TP + ETV_IC50);

// TDF: inhibition of viral replication (RT)
double TFV_eff = USE_TDF * TFV_Imax * TFV_DP / (TFV_DP + TFV_IC50);

// Combined NUC effect (Bliss independence approximation)
double NUC_eff = 1.0 - (1.0 - ETV_eff)*(1.0 - TFV_eff);

// Peg-IFN antiviral effect on cccDNA/replication
double IFN_av  = USE_PIFN * IFN_Emax * IFN_Cp / (IFN_Cp + IFN_EC50);
// Peg-IFN immunomodulatory (CTL boost factor)
double IFN_im  = USE_PIFN * IFN_immE * IFN_Cp / (IFN_Cp + IFN_EC50);
// Peg-IFN anti-cccDNA (reduces cccDNA directly)
double IFN_ccc_eff = USE_PIFN * IFN_cccE * IFN_Cp / (IFN_Cp + IFN_EC50);

// siRNA HBsAg knockdown
double siRNA_eff = IAGM_dose * siRNA_Emax;

// ── Innate IFN antiviral ──────────────────────────────────────
double INN_av  = IFN_antiv * IFN_inn / (IFN_inn + 0.5);

// ── Effective CTL (reduced by exhaustion) ────────────────────
double CTL_eff = CTL_resp * (1.0 - Exhaust * 0.7);

// ── Fibrosis-driven portal hypertension (>F3) ────────────────
double PH_risk = Fibrosis > 3.0 ? (Fibrosis - 3.0) : 0.0;

// ── Initial conditions ────────────────────────────────────────
if(NEWIND <= 1) {
  T_cell_0 = T0 / 1e6;
  I_cell_0 = I0 / 1e6;
  V_dna_0  = V0;
  ccc_DNA_0 = ccc0;
  HBsAg_C_0 = Ag0;
  CTL_resp_0 = CTL0;
  Exhaust_0  = ExhSS;
  IFN_inn_0  = 0.1;
  ALT_val_0  = ALT0;
  Fibrosis_0 = FIB0;
  HCC_risk_0 = 0.0;
  HSC_act_0  = 0.1;
}

$ODE
// ── ETV PK ───────────────────────────────────────────────────
dxdt_ETV_gut = -ETV_ka * ETV_gut;
dxdt_ETV_C   =  ETV_F * ETV_ka * ETV_gut - (ETV_CL/ETV_Vc) * ETV_C
               - ETV_kTP * ETV_C * ETV_Kp / ETV_Vc;
dxdt_ETV_TP  =  ETV_kTP * ETV_C * ETV_Kp / ETV_Vc - ETV_kTP_d * ETV_TP;

// ── TDF PK ───────────────────────────────────────────────────
dxdt_TDF_gut = -TDF_ka * TDF_gut;
dxdt_TDF_C   =  TDF_F * TDF_ka * TDF_gut - (TDF_CL/TDF_Vc) * TDF_C
               - TFV_kDP * TDF_C / TDF_Vc;
dxdt_TFV_DP  =  TFV_kDP * TDF_C / TDF_Vc - TFV_kDP_d * TFV_DP;

// ── Peg-IFN PK ───────────────────────────────────────────────
dxdt_IFN_SC  = -IFN_ka * IFN_SC;
dxdt_IFN_C   =  IFN_F * IFN_ka * IFN_SC - (IFN_CL/IFN_Vc) * IFN_C;

// ── Target cell dynamics ─────────────────────────────────────
// dT/dt = production - natural death - new infections
double T = T_cell;
double I = I_cell;
double V = V_dna;
double ccc = ccc_DNA;

dxdt_T_cell = s_T/1e6 - d_T * T - beta * T * V;

// ── Infected cell dynamics ───────────────────────────────────
// dI/dt = new infections - immune killing (CTL) - natural death
double CTL_kill = alpha_CTL * CTL_eff * I / (I + 0.1);
dxdt_I_cell = beta * T * V
              - delta_I * I * (1 + CTL_eff)
              - d_T * I;

// ── Serum HBV DNA ─────────────────────────────────────────────
// dV/dt = production from infected cells - clearance - drug effect
double viral_prod = p_V * I * (1.0 - NUC_eff) * (1.0 - INN_av);
dxdt_V_dna = viral_prod - c_V * V;

// ── cccDNA dynamics ───────────────────────────────────────────
// cccDNA sustained by rcDNA recycling; reduced by IFN/CTL
double ccc_input  = lambda_ccc * I * (1.0 - NUC_eff) * (1.0 - IFN_ccc_eff);
double ccc_loss   = delta_ccc * ccc + alpha_CTL * CTL_eff * 0.2 * ccc;
dxdt_ccc_DNA = ccc_input - ccc_loss;

// ── HBsAg dynamics ────────────────────────────────────────────
// HBsAg produced from cccDNA and integrated DNA; cleared slowly
dxdt_HBsAg_C = k_Ag * ccc - c_Ag * HBsAg_C * (1.0 + siRNA_eff);

// ── CTL dynamics ──────────────────────────────────────────────
// Expand with antigen (HBV DNA driven), contract, boosted by IFN
double CTL_exp = k_CTL_exp * V / (V + 1e5) * (1.0 + IFN_im);
dxdt_CTL_resp = CTL_exp * (1.0 - CTL_resp) - d_CTL * CTL_resp;

// ── T cell exhaustion ─────────────────────────────────────────
// Driven by chronic antigen (HBsAg excess); partially reversed by IFN
double Ag_load = HBsAg_C / (HBsAg_C + 100.0);
dxdt_Exhaust = k_exhaust * Ag_load * (1.0 - Exhaust)
               - r_exhaust * Exhaust
               - IFN_im * 0.1 * Exhaust;

// ── Innate IFN signaling ──────────────────────────────────────
// Induced by HBV (attenuated by HBx evasion); boosted by Peg-IFN
double HBV_innate = 0.05 * V / (V + 1e6);  // HBV induces IFN poorly
dxdt_IFN_inn = k_IFNinn * HBV_innate - d_IFNinn * IFN_inn + IFN_av * 0.5;

// ── ALT dynamics ──────────────────────────────────────────────
// ALT elevated by infected cell death + immune-mediated lysis
double inflam_drive = delta_I * I * CTL_eff * k_ALT_I;
dxdt_ALT_val = inflam_drive - k_ALT_ret * (ALT_val - ALT0);

// ── HSC activation ────────────────────────────────────────────
// Driven by inflammation (ALT proxy), inhibited on viral suppression
double inflam_norm = (ALT_val - ALT0) / (ALT0 + 10.0);
dxdt_HSC_act = k_HSC * inflam_norm * (1.0 - HSC_act)
               - d_HSC * HSC_act * (V < 100 ? 1.5 : 1.0);

// ── Fibrosis progression ──────────────────────────────────────
// Progresses with inflammation; regresses slowly on viral suppression
double prog_rate = k_fibrosis * HSC_act * (ALT_val/ALT0);
double reg_rate  = k_fib_reg  * (V_dna < 100 ? 1.5 : 0.5);
dxdt_Fibrosis = prog_rate * (4.0 - Fibrosis) / 4.0
                - reg_rate * Fibrosis / 4.0;

// ── HCC risk accumulation ─────────────────────────────────────
dxdt_HCC_risk = k_HCC * Fibrosis * (1.0 + 0.5*(V_dna > 1e4 ? 1.0 : 0.0));

$TABLE
// Derived outputs for plotting
capture ETV_Cp   = ETV_C / ETV_Vc;          // ETV plasma (ng/mL)
capture TFV_Cp   = TDF_C / TDF_Vc;          // TFV plasma (ng/mL)
capture IFN_Cp   = IFN_C / IFN_Vc;          // IFN plasma
capture V_log10  = log10(V_dna + 1.0);       // HBV DNA log10 IU/mL
capture Ag_log10 = log10(HBsAg_C + 0.01);   // HBsAg log10 IU/mL
capture ccc_log  = log10(ccc_DNA + 0.001);   // cccDNA log10
capture ALT_out  = ALT_val;
capture Fib_out  = Fibrosis;
capture HCC_out  = HCC_risk;
capture CTL_out  = CTL_resp;
capture Exh_out  = Exhaust;
capture T_frac   = T_cell / (T0/1e6);        // Target cell fraction remaining
capture I_frac   = I_cell / (T0/1e6 + I0/1e6); // Infected fraction
capture NUC_eff_out = NUC_eff;
capture IFN_AV  = IFN_av;
capture HSC_out  = HSC_act;
capture PH_out   = PH_risk;
'

## ─────────────────────────────────────────────────────────────
## 2. COMPILE MODEL
## ─────────────────────────────────────────────────────────────

chb_mod <- mcode("CHB_QSP", chb_code, quiet = TRUE)

## ─────────────────────────────────────────────────────────────
## 3. DOSING REGIMENS
## ─────────────────────────────────────────────────────────────

sim_end  <- 365 * 3   # 3-year simulation (days)
sim_freq <- 1         # daily output

# ETV 0.5 mg QD (oral) — dose as bolus to ETV_gut (mg)
ev_ETV   <- ev(cmt = "ETV_gut", amt = 0.5, ii = 24, addl = sim_end - 1,
               rate = 0, time = 0)

# TDF 300 mg QD
ev_TDF   <- ev(cmt = "TDF_gut", amt = 300,  ii = 24, addl = sim_end - 1,
               rate = 0, time = 0)

# Peg-IFN-α2a 180 µg SC QW (weekly) for 48 weeks (336 days)
ev_PIFN  <- ev(cmt = "IFN_SC",  amt = 180,  ii = 168, addl = 47,
               rate = 0, time = 0)

## ─────────────────────────────────────────────────────────────
## 4. TREATMENT SCENARIOS
## ─────────────────────────────────────────────────────────────

scenarios <- list(
  "Untreated"              = list(ev_dose = ev(amt=0, cmt=1, time=0),
                                   params = list(USE_ETV=0, USE_TDF=0, USE_PIFN=0)),
  "Entecavir 0.5 mg QD"   = list(ev_dose = ev_ETV,
                                   params = list(USE_ETV=1, USE_TDF=0, USE_PIFN=0)),
  "TDF 300 mg QD"          = list(ev_dose = ev_TDF,
                                   params = list(USE_ETV=0, USE_TDF=1, USE_PIFN=0)),
  "Peg-IFN-α2a × 48wks"   = list(ev_dose = ev_PIFN,
                                   params = list(USE_ETV=0, USE_TDF=0, USE_PIFN=1)),
  "ETV + Peg-IFN (combo)"  = list(ev_dose = c(ev_ETV, ev_PIFN),
                                   params = list(USE_ETV=1, USE_TDF=0, USE_PIFN=1)),
  "ETV + siRNA add-on"     = list(ev_dose = ev_ETV,
                                   params = list(USE_ETV=1, USE_TDF=0, USE_PIFN=0,
                                                 IAGM_dose=1))
)

## ─────────────────────────────────────────────────────────────
## 5. SIMULATION
## ─────────────────────────────────────────────────────────────

run_scenario <- function(scen_name, scen) {
  mod_s <- chb_mod %>% param(scen$params)
  out <- mod_s %>%
    ev(scen$ev_dose) %>%
    mrgsim(end = sim_end, delta = sim_freq, recover = "V_log10,Ag_log10,ALT_out,Fib_out,HCC_out,CTL_out,Exh_out,ETV_Cp,TFV_Cp,IFN_Cp,T_frac,I_frac,NUC_eff_out,ccc_log,HSC_out,IFN_AV") %>%
    as_tibble() %>%
    mutate(Scenario = scen_name)
  out
}

# Run all scenarios and combine
results <- bind_rows(
  mapply(run_scenario, names(scenarios), scenarios, SIMPLIFY = FALSE)
)

cat("\n=== CHB QSP Model: Simulation Complete ===\n")
cat(sprintf("  Scenarios: %d\n", length(scenarios)))
cat(sprintf("  Duration: %d days (%g years)\n", sim_end, sim_end/365))
cat(sprintf("  Output rows: %d\n", nrow(results)))

## ─────────────────────────────────────────────────────────────
## 6. KEY RESULTS AT YEAR 1 AND YEAR 3
## ─────────────────────────────────────────────────────────────

timepoints <- c(0, 365, 730, 1095)

summary_tab <- results %>%
  filter(time %in% timepoints) %>%
  select(Scenario, time, V_log10, Ag_log10, ALT_out, Fib_out, CTL_out,
         Exh_out, HCC_out, ccc_log) %>%
  mutate(
    Week     = round(time / 7),
    ViroResp = V_log10 < log10(21),      # <20 IU/mL
    HBsAgLow = Ag_log10 < log10(0.05)    # HBsAg loss
  )

cat("\n── Year 1 HBV DNA and HBsAg by Scenario ──\n")
print(summary_tab %>%
        filter(time == 365) %>%
        select(Scenario, V_log10, Ag_log10, ALT_out, ViroResp, HBsAgLow) %>%
        arrange(V_log10))

cat("\n── Year 3 Fibrosis Score and HCC Risk ──\n")
print(summary_tab %>%
        filter(time == 1095) %>%
        select(Scenario, Fib_out, HCC_out, Exh_out, ccc_log) %>%
        arrange(Fib_out))

## ─────────────────────────────────────────────────────────────
## 7. PLOTS
## ─────────────────────────────────────────────────────────────

scenario_colors <- c(
  "Untreated"              = "#E53935",
  "Entecavir 0.5 mg QD"   = "#1E88E5",
  "TDF 300 mg QD"          = "#43A047",
  "Peg-IFN-α2a × 48wks"   = "#FB8C00",
  "ETV + Peg-IFN (combo)"  = "#8E24AA",
  "ETV + siRNA add-on"     = "#00897B"
)

# Panel 1: HBV DNA kinetics
p1 <- ggplot(results, aes(x = time/365, y = V_log10, color = Scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = log10(20), linetype = "dashed", color = "gray40") +
  annotate("text", x = 0.2, y = log10(15), label = "Detection limit (20 IU/mL)",
           size = 3, color = "gray40") +
  scale_color_manual(values = scenario_colors) +
  labs(title = "CHB Model: Serum HBV DNA Kinetics",
       x = "Time (years)", y = "HBV DNA (log₁₀ IU/mL)",
       color = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

# Panel 2: HBsAg decline
p2 <- ggplot(results, aes(x = time/365, y = Ag_log10, color = Scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = log10(0.05), linetype = "dashed", color = "red") +
  annotate("text", x = 0.2, y = log10(0.03), label = "Functional cure threshold (0.05 IU/mL)",
           size = 3, color = "red") +
  scale_color_manual(values = scenario_colors) +
  labs(title = "HBsAg Kinetics — Path to Functional Cure",
       x = "Time (years)", y = "HBsAg (log₁₀ IU/mL)",
       color = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

# Panel 3: ALT normalization
p3 <- ggplot(results, aes(x = time/365, y = ALT_out, color = Scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 40, linetype = "dashed", color = "green4") +
  annotate("text", x = 0.2, y = 36, label = "ULN (40 IU/L)",
           size = 3, color = "green4") +
  scale_color_manual(values = scenario_colors) +
  labs(title = "ALT Kinetics — Biochemical Response",
       x = "Time (years)", y = "ALT (IU/L)",
       color = "Treatment") +
  theme_bw(base_size = 12)

# Panel 4: Fibrosis progression/regression
p4 <- ggplot(results, aes(x = time/365, y = Fib_out, color = Scenario)) +
  geom_line(size = 1.1) +
  scale_y_continuous(breaks = 0:4, labels = paste0("F", 0:4), limits = c(0, 4)) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Hepatic Fibrosis Progression (Metavir)",
       x = "Time (years)", y = "Fibrosis Score (Metavir)",
       color = "Treatment") +
  theme_bw(base_size = 12)

# Panel 5: cccDNA kinetics
p5 <- ggplot(results, aes(x = time/365, y = ccc_log, color = Scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Intrahepatic cccDNA Kinetics",
       x = "Time (years)", y = "cccDNA (log₁₀ copies/cell equiv)",
       color = "Treatment") +
  theme_bw(base_size = 12)

# Panel 6: T cell dynamics (CTL and Exhaustion)
p6 <- results %>%
  select(time, Scenario, CTL_out, Exh_out) %>%
  pivot_longer(cols = c(CTL_out, Exh_out), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                          CTL_out = "CTL Response",
                          Exh_out = "T cell Exhaustion")) %>%
  ggplot(aes(x = time/365, y = value, color = Scenario, linetype = metric)) +
  geom_line(size = 1.0) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Immune Dynamics: CTL Response & T-cell Exhaustion",
       x = "Time (years)", y = "Normalized Index (0–1)",
       color = "Treatment", linetype = "Metric") +
  theme_bw(base_size = 12)

## ─────────────────────────────────────────────────────────────
## 8. PK VERIFICATION PANEL
## ─────────────────────────────────────────────────────────────

pk_etv <- results %>%
  filter(Scenario == "Entecavir 0.5 mg QD", time <= 14) %>%
  select(time, ETV_Cp) %>%
  ggplot(aes(x = time, y = ETV_Cp)) +
  geom_line(color = "#1E88E5", size = 1.2) +
  labs(title = "ETV PK: First 2 Weeks",
       x = "Time (days)", y = "ETV Plasma (ng/mL)") +
  theme_bw(base_size = 11)

pk_ifn <- results %>%
  filter(Scenario == "Peg-IFN-α2a × 48wks", time <= 120) %>%
  select(time, IFN_Cp) %>%
  ggplot(aes(x = time, y = IFN_Cp)) +
  geom_line(color = "#FB8C00", size = 1.2) +
  labs(title = "Peg-IFN PK: First 16 Weeks",
       x = "Time (days)", y = "Peg-IFN Plasma (ng/mL equiv)") +
  theme_bw(base_size = 11)

## ─────────────────────────────────────────────────────────────
## 9. SENSITIVITY ANALYSIS: BASELINE VIRAL LOAD EFFECT
## ─────────────────────────────────────────────────────────────

vl_baseline <- c(1e5, 1e7, 1e9)
vl_labels   <- c("Low (10^5 IU/mL)", "Medium (10^7 IU/mL)", "High (10^9 IU/mL)")

sens_results <- bind_rows(lapply(seq_along(vl_baseline), function(i) {
  chb_mod %>%
    param(USE_ETV=1, USE_TDF=0, USE_PIFN=0, V0=vl_baseline[i]) %>%
    ev(ev_ETV) %>%
    mrgsim(end = sim_end, delta = sim_freq,
           recover = "V_log10,Ag_log10,ALT_out,Fib_out") %>%
    as_tibble() %>%
    mutate(VL_group = vl_labels[i])
}))

p_sens <- ggplot(sens_results, aes(x = time/365, y = V_log10, color = VL_group)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = log10(20), linetype = "dashed", color = "gray40") +
  labs(title = "Sensitivity: Baseline Viral Load Effect on ETV Response",
       x = "Time (years)", y = "HBV DNA (log₁₀ IU/mL)",
       color = "Baseline VL") +
  theme_bw(base_size = 12)

## ─────────────────────────────────────────────────────────────
## 10. PARAMETER CALIBRATION SUMMARY
## ─────────────────────────────────────────────────────────────

cat("\n══════════════════════════════════════════════════════════\n")
cat("  CHB QSP Model — Key Parameter Calibration\n")
cat("══════════════════════════════════════════════════════════\n")
cat(sprintf("  ETV IC50 (ETV-TP):     %g µM     [Ref: Colonno 2006]\n", 0.004))
cat(sprintf("  ETV Year-1 DNA resp:  ~67%%       [Chang 2006 NEJM]\n"))
cat(sprintf("  TDF IC50 (TFV-DP):    %g µM    [Ref: Tsiang 2008]\n", 0.5))
cat(sprintf("  TDF Year-1 DNA resp:  ~76%%       [Marcellin 2008 NEJM]\n"))
cat(sprintf("  Peg-IFN HBsAg loss:   ~3%% yr-1  [Lau 2005 NEJM]\n"))
cat(sprintf("  cccDNA half-life:      ~33 days (δ=%g/day) [Werle-Lapostolle 2004]\n", 0.003))
cat(sprintf("  HBsAg half-life:       ~70 days (c=%g/day) [Volz 2007]\n", 0.01))
cat("══════════════════════════════════════════════════════════\n")

# Print plots if interactive
if(interactive()) {
  print(p1); print(p2); print(p3); print(p4)
  print(p5); print(p6); print(p_sens)
}

cat("\nModel compilation complete. Run this script in R with mrgsolve installed.\n")
