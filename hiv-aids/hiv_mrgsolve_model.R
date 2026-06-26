## ============================================================
## HIV/AIDS QSP Model — mrgsolve ODE Simulation
## ============================================================
## Viral dynamics: Perelson et al. (1996, 1997) Science;
##   Nowak & May (2000) "Virus Dynamics"
## PK: FDA CPR reviews for TDF, FTC, DTG, BIC, EFV, DRV/r
## PD: IC50 from published in vitro / protein-adjusted data
## 18 compartments · 8 treatment scenarios
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ─────────────────────────────────────────────────────────────
# mrgsolve MODEL CODE
# ─────────────────────────────────────────────────────────────
code <- '
$PROB HIV/AIDS QSP Model (Perelson viral dynamics + ART PK/PD)

$PARAM
  // Drug use switches (0=off, 1=on)
  use_TDF = 0, use_FTC = 0, use_DTG = 0,
  use_BIC = 0, use_EFV = 0, use_DRV = 0,

  // TDF/TFV PK (FDA NDA 21-356 review)
  ka_TDF   = 1.0,    // absorption rate (h-1)
  CL_TDF   = 16.0,   // plasma CL (L/h)
  Vc_TDF   = 45.0,   // Vd central (L)
  k_IC_TDF = 0.010,  // TFV→TFV-DP intracell formation (h-1)
  k_eg_TDF = 0.00462,// TFV-DP elimination (h-1; t1/2=150 h)

  // FTC PK (FDA NDA 21-500 review)
  ka_FTC   = 1.5,    // h-1
  CL_FTC   = 13.7,   // L/h
  Vc_FTC   = 213.0,  // L
  k_IC_FTC = 0.020,  // FTC→FTC-TP (h-1)
  k_eg_FTC = 0.01782,// FTC-TP elimination (h-1; t1/2=39 h)

  // DTG PK (FDA NDA 204790 review; Edelman 2014 JAIDS)
  ka_DTG   = 0.9,    // h-1
  CL_DTG   = 1.0,    // L/h (UGT1A1/CYP3A4)
  Vc_DTG   = 17.4,   // L

  // BIC PK (FDA NDA 210251; Gallant 2017 Lancet)
  ka_BIC   = 0.5,    // h-1
  CL_BIC   = 0.41,   // L/h
  Vc_BIC   = 25.0,   // L

  // EFV PK (FDA NDA 20-972; extensive CYP2B6 metabolizer)
  ka_EFV   = 0.5,    // h-1
  CL_EFV   = 4.5,    // L/h
  Vc_EFV   = 270.0,  // L

  // DRV/r PK (FDA NDA 22-341; RTV-boosted)
  ka_DRV   = 1.2,    // h-1
  CL_DRV   = 1.0,    // L/h (CYP3A4 inhibited by RTV)
  Vc_DRV   = 88.4,   // L

  // PD IC50 (protein-adjusted; nM)
  IC50_TFV = 100.0,  // TFV-DP IC50 (nM; Robbins 2003 AAC)
  IC50_FTC = 1000.0, // FTC-TP IC50 (nM)
  IC50_DTG = 2.0,    // DTG IC50 (nM; Veltri 2010 AAC)
  IC50_BIC = 1.7,    // BIC IC50 (nM; Tsiang 2016 AAC)
  IC50_EFV = 3000.0, // EFV IC50 (nM; Ren 2000 Structure)
  IC50_DRV = 0.5,    // DRV IC50 (nM; Koh 2003 JBC)
  n_H      = 1.5,    // Hill coefficient

  // HIV viral dynamics (Perelson 1996 Science; t in days)
  s_T      = 10.0,   // CD4 naive source (cells/uL/d)
  d_T      = 0.01,   // CD4 natural death (d-1)
  beta     = 5e-8,   // infection rate const. (mL/copy/d)
  delta_I  = 1.0,    // infected CD4 death (d-1)
  p_V      = 3000.0, // viral production per infected cell (copies/mL/d per cell/uL)
  c_V      = 23.0,   // viral clearance (d-1)
  k_lat    = 0.02,   // fraction → latent reservoir
  r_lat    = 1e-5,   // latent reactivation (d-1)

  // CTL dynamics (Nowak & May 2000)
  s_E      = 5.0,    // CTL source (cells/uL/d)
  p_E      = 0.5,    // CTL expansion (d-1)
  d_E      = 0.05,   // CTL death (d-1)
  k_kill   = 0.004,  // CTL killing rate (uL/cell/d)

  // Inflammation
  k_IL6    = 0.01,   // IL-6 prod from infected cells (pg/mL per cell/uL per d)
  d_IL6    = 0.10,   // IL-6 decay (d-1)
  IL6_0    = 2.0,    // baseline IL-6 (pg/mL)

  // Resistance
  k_mut    = 1e-6,   // mutation rate (d-1 × copies)
  k_rd     = 0.001   // resistance decay with fitness cost (d-1)

$CMT
  // Drug PK
  TDF_GUT      // TDF gut depot (mg)
  TDF_PLASMA   // TFV plasma (ug/mL)
  TFV_DP       // TFV-DP intracellular (nM)
  FTC_PLASMA   // FTC plasma (ug/mL)
  FTC_TP       // FTC-TP intracellular (nM)
  DTG_PLASMA   // DTG plasma (ug/mL)
  BIC_PLASMA   // BIC plasma (ug/mL)
  EFV_PLASMA   // EFV plasma (ug/mL)
  DRV_PLASMA   // DRV plasma (ug/mL)
  // HIV Viral dynamics
  T_CELL       // Uninfected CD4+ T cells (cells/uL)
  I_CELL       // Productively infected CD4+ (cells/uL)
  V_FREE       // Free virus (copies/mL ×1e-3)
  L_CELL       // Latently infected cells (cells/uL)
  E_CELL       // CD8+ CTL effectors (cells/uL)
  // Biomarkers / state
  INFLAM       // IL-6 systemic inflammation (pg/mL)
  RESIST       // Resistance score (0–1)
  VL_LOG       // log10(VL copies/mL) smoothed
  CD4_SMOOTH   // Smoothed CD4 count (cells/uL)

$INIT
  TDF_GUT=0, TDF_PLASMA=0, TFV_DP=0,
  FTC_PLASMA=0, FTC_TP=0,
  DTG_PLASMA=0, BIC_PLASMA=0,
  EFV_PLASMA=0, DRV_PLASMA=0,
  T_CELL=1000, I_CELL=0.01, V_FREE=1.0,
  L_CELL=0.001, E_CELL=100,
  INFLAM=2.0, RESIST=0,
  VL_LOG=3.0, CD4_SMOOTH=1000

$ODE
  // ── Drug PK ────────────────────────────────────────────────
  // TDF: gut→plasma (1-comp first-order abs) + intracell TFV-DP
  double ke_TDF = CL_TDF / Vc_TDF;       // h-1
  dxdt_TDF_GUT    = -ka_TDF * TDF_GUT;
  dxdt_TDF_PLASMA =  ka_TDF * TDF_GUT / Vc_TDF - ke_TDF * TDF_PLASMA;
  // TFV plasma→TFV-DP intracell (TFV MW=287.2; µg/mL → nM × 3484)
  double TFV_nM = TDF_PLASMA * 1000.0 / 287.2;
  dxdt_TFV_DP   = k_IC_TDF * TFV_nM - k_eg_TDF * TFV_DP;

  // FTC: 1-comp + intracell FTC-TP (FTC MW=247.2; µg/mL → nM × 4044)
  double ke_FTC = CL_FTC / Vc_FTC;
  dxdt_FTC_PLASMA = -ke_FTC * FTC_PLASMA;
  double FTC_nM = FTC_PLASMA * 1000.0 / 247.2;
  dxdt_FTC_TP   = k_IC_FTC * FTC_nM - k_eg_FTC * FTC_TP;

  // DTG (MW=419.4)
  double ke_DTG = CL_DTG / Vc_DTG;
  dxdt_DTG_PLASMA = -ke_DTG * DTG_PLASMA;

  // BIC (MW=449.4)
  double ke_BIC = CL_BIC / Vc_BIC;
  dxdt_BIC_PLASMA = -ke_BIC * BIC_PLASMA;

  // EFV (MW=315.7)
  double ke_EFV = CL_EFV / Vc_EFV;
  dxdt_EFV_PLASMA = -ke_EFV * EFV_PLASMA;

  // DRV/r (MW=547.7)
  double ke_DRV = CL_DRV / Vc_DRV;
  dxdt_DRV_PLASMA = -ke_DRV * DRV_PLASMA;

  // ── Drug PD: antiviral efficacy ─────────────────────────────
  // All concentrations in nM
  double DTG_nM = DTG_PLASMA * 1000.0 / 419.4;
  double BIC_nM = BIC_PLASMA * 1000.0 / 449.4;
  double EFV_nM = EFV_PLASMA * 1000.0 / 315.7;
  double DRV_nM = DRV_PLASMA * 1000.0 / 547.7;

  // Emax Hill model: η = C^n / (IC50^n + C^n)
  double eTFV = pow(TFV_DP, n_H) / (pow(IC50_TFV, n_H) + pow(TFV_DP, n_H));
  double eFTC = pow(FTC_TP,  n_H) / (pow(IC50_FTC, n_H) + pow(FTC_TP,  n_H));
  double eDTG = pow(DTG_nM,  n_H) / (pow(IC50_DTG, n_H) + pow(DTG_nM,  n_H));
  double eBIC = pow(BIC_nM,  n_H) / (pow(IC50_BIC, n_H) + pow(BIC_nM,  n_H));
  double eEFV = pow(EFV_nM,  n_H) / (pow(IC50_EFV, n_H) + pow(EFV_nM,  n_H));
  double eDRV = pow(DRV_nM,  n_H) / (pow(IC50_DRV, n_H) + pow(DRV_nM,  n_H));

  // Combined NRTI: TFV-DP & FTC-TP (independent inhibition)
  double eNRTI = 1.0 - (1.0 - eTFV * use_TDF) * (1.0 - eFTC * use_FTC);
  // RT inhibition: NRTI + NNRTI
  double eRT   = 1.0 - (1.0 - eNRTI) * (1.0 - eEFV * use_EFV);
  // INSTI: DTG or BIC
  double eINSTI = (use_DTG > 0.5) ? eDTG : ((use_BIC > 0.5) ? eBIC : 0.0);
  // Overall efficacy (three-class: RT, INSTI, PI)
  double eta = 1.0 - (1.0 - eRT) * (1.0 - eINSTI) * (1.0 - eDRV * use_DRV);
  if (eta > 0.9999) eta = 0.9999;
  if (eta < 0.0)    eta = 0.0;

  // ── HIV Viral Dynamics (ODEs in days) ──────────────────────
  // Guard against negative state
  double T = (T_CELL > 0) ? T_CELL : 0.0;
  double I = (I_CELL > 0) ? I_CELL : 0.0;
  double V = (V_FREE > 0) ? V_FREE : 0.0;
  double L = (L_CELL > 0) ? L_CELL : 0.0;
  double E = (E_CELL > 0) ? E_CELL : 0.0;
  double R = (RESIST > 0) ? RESIST : 0.0;

  // Effective infection rate elevated by resistance
  double beta_eff = beta * (1.0 - eta) * (1.0 + 2.0 * R);

  // Uninfected CD4+ (cells/uL/day)
  dxdt_T_CELL = s_T - d_T * T - beta_eff * V * T * (1.0 - k_lat) + r_lat * L;

  // Productively infected CD4+ (cells/uL/day)
  dxdt_I_CELL = beta_eff * V * T * (1.0 - k_lat) - delta_I * I - k_kill * E * I;

  // Free virus (copies/mL ×1e-3 / day)
  // Virion production partially suppressed by PI (half-efficacy approximation)
  double pi_supp = 1.0 - eDRV * use_DRV * 0.5;
  dxdt_V_FREE = p_V * I * pi_supp - c_V * V + p_V * r_lat * L * 0.1;

  // Latent reservoir (cells/uL/day)
  dxdt_L_CELL = beta_eff * V * T * k_lat - r_lat * L - k_kill * E * L * 0.01;

  // CTL effectors (cells/uL/day; Nowak & May)
  dxdt_E_CELL = s_E + p_E * I / (0.5 + I) - d_E * E;

  // Inflammation – IL-6 (pg/mL/day)
  dxdt_INFLAM = k_IL6 * I - d_IL6 * (INFLAM - IL6_0);

  // Resistance score (0–1)
  dxdt_RESIST = k_mut * V * (1.0 - eta) - k_rd * R;

  // Smoothed observables (first-order filter, half-life ~1.4 days)
  double VL_now = (V * 1000.0 > 1.0) ? log10(V * 1000.0) : 0.0;
  dxdt_VL_LOG    = (VL_now - VL_LOG) * 0.5;
  dxdt_CD4_SMOOTH = ((T - I) - CD4_SMOOTH) * 0.5;

$TABLE
  // Plasma concentrations (µg/mL)
  double C_TFV  = TDF_PLASMA;
  double C_FTC  = FTC_PLASMA;
  double C_DTG  = DTG_PLASMA;
  double C_BIC  = BIC_PLASMA;
  double C_EFV  = EFV_PLASMA;
  double C_DRV  = DRV_PLASMA;

  // Intracellular active metabolites (nM)
  double TFV_DP_nM = TFV_DP;
  double FTC_TP_nM = FTC_TP;

  // Clinical outputs
  double CD4_count  = (T_CELL - I_CELL > 0) ? (T_CELL - I_CELL) : 0.0;
  double VL_copies  = V_FREE * 1000.0;
  double VL_log10   = (VL_copies > 1.0) ? log10(VL_copies) : 0.0;
  double VL_supp    = (VL_copies < 50.0) ? 1.0 : 0.0;   // Virologic suppression
  double AIDS_risk  = (CD4_count < 200.0) ? 1.0 : 0.0;  // CD4 <200
  double Lat_IUPM   = L_CELL * 1e3;                      // IUPM proxy
  double IL6_conc   = INFLAM;
  double Res_score  = RESIST;
  double CTL_count  = E_CELL;

  // Recompute efficacy for TABLE output
  double TFV_n = TDF_PLASMA * 1000.0 / 287.2;
  double FTC_n = FTC_PLASMA * 1000.0 / 247.2;
  double DTG_n = DTG_PLASMA * 1000.0 / 419.4;
  double BIC_n = BIC_PLASMA * 1000.0 / 449.4;
  double EFV_n = EFV_PLASMA * 1000.0 / 315.7;
  double DRV_n = DRV_PLASMA * 1000.0 / 547.7;

  double eTFV2  = pow(TFV_DP,n_H)/(pow(IC50_TFV,n_H)+pow(TFV_DP,n_H));
  double eFTC2  = pow(FTC_TP,n_H)/(pow(IC50_FTC,n_H)+pow(FTC_TP,n_H));
  double eDTG2  = pow(DTG_n,n_H)/(pow(IC50_DTG,n_H)+pow(DTG_n,n_H));
  double eBIC2  = pow(BIC_n,n_H)/(pow(IC50_BIC,n_H)+pow(BIC_n,n_H));
  double eEFV2  = pow(EFV_n,n_H)/(pow(IC50_EFV,n_H)+pow(EFV_n,n_H));
  double eDRV2  = pow(DRV_n,n_H)/(pow(IC50_DRV,n_H)+pow(DRV_n,n_H));
  double eNRTI2 = 1.0-(1.0-eTFV2*use_TDF)*(1.0-eFTC2*use_FTC);
  double eRT2   = 1.0-(1.0-eNRTI2)*(1.0-eEFV2*use_EFV);
  double eI2    = (use_DTG>0.5)?eDTG2:((use_BIC>0.5)?eBIC2:0.0);
  double Eta    = 1.0-(1.0-eRT2)*(1.0-eI2)*(1.0-eDRV2*use_DRV);
  if(Eta>0.9999) Eta=0.9999;
  if(Eta<0)      Eta=0.0;

$CAPTURE
  C_TFV C_FTC C_DTG C_BIC C_EFV C_DRV
  TFV_DP_nM FTC_TP_nM
  CD4_count VL_copies VL_log10 VL_supp AIDS_risk
  Lat_IUPM IL6_conc Res_score Eta CTL_count
'

# ─────────────────────────────────────────────────────────────
# Compile
# ─────────────────────────────────────────────────────────────
mod <- mrgsolve::mcode("hiv_qsp", code)

# ─────────────────────────────────────────────────────────────
# Dosing event builder
# ─────────────────────────────────────────────────────────────
make_ev <- function(drugs, t_start, t_end) {
  # Doses converted to plasma compartment initial bolus via ka absorption
  # TDF: amt in gut depot (mg);  others: direct plasma (µg/mL approximation)
  dose_params <- list(
    TDF = list(cmt="TDF_GUT",    amt=300,             ii=24),
    FTC = list(cmt="FTC_PLASMA", amt=200/213.0*0.93,  ii=24),
    DTG = list(cmt="DTG_PLASMA", amt=50/17.4*0.53,    ii=24),
    BIC = list(cmt="BIC_PLASMA", amt=50/25.0*0.95,    ii=24),
    EFV = list(cmt="EFV_PLASMA", amt=600/270.0*0.42,  ii=24),
    DRV = list(cmt="DRV_PLASMA", amt=800/88.4*0.82,   ii=24)
  )
  ev_list <- lapply(drugs, function(d) {
    dp   <- dose_params[[d]]
    addl <- max(0L, as.integer(floor((t_end - t_start) / dp$ii)))
    ev(cmt=dp$cmt, amt=dp$amt, ii=dp$ii, addl=addl, time=t_start)
  })
  do.call(c, ev_list)
}

# ─────────────────────────────────────────────────────────────
# Scenario definitions
# ─────────────────────────────────────────────────────────────
SIM_DAYS <- 5 * 365  # 5 years

scenarios <- list(
  # 1. No ART
  list(id=1,
       name="① No ART (Natural Progression)",
       params=list(use_TDF=0,use_FTC=0,use_DTG=0,use_BIC=0,use_EFV=0,use_DRV=0),
       events=NULL),

  # 2. TDF/FTC/DTG (Triumeq equivalent — first-line)
  list(id=2,
       name="② TDF/FTC/DTG (Triumeq-type 1st-line)",
       params=list(use_TDF=1,use_FTC=1,use_DTG=1,use_BIC=0,use_EFV=0,use_DRV=0),
       events=make_ev(c("TDF","FTC","DTG"), 0, SIM_DAYS)),

  # 3. TAF/FTC/BIC (Biktarvy — contemporary 1st-line)
  # Simulate TAF by increasing k_IC_TDF (10× TFV-DP vs TDF)
  list(id=3,
       name="③ TAF/FTC/BIC (Biktarvy 1st-line)",
       params=list(use_TDF=1,use_FTC=1,use_DTG=0,use_BIC=1,use_EFV=0,use_DRV=0,
                   k_IC_TDF=0.05, ka_TDF=2.5),  # TAF kinetics
       events=make_ev(c("TDF","FTC","BIC"), 0, SIM_DAYS)),

  # 4. TDF/FTC/EFV (WHO LMIC first-line)
  list(id=4,
       name="④ TDF/FTC/EFV (WHO LMIC 1st-line)",
       params=list(use_TDF=1,use_FTC=1,use_DTG=0,use_BIC=0,use_EFV=1,use_DRV=0),
       events=make_ev(c("TDF","FTC","EFV"), 0, SIM_DAYS)),

  # 5. DRV/r + DTG (Salvage / 2nd-line)
  list(id=5,
       name="⑤ DRV/r + DTG (Salvage 2nd-line)",
       params=list(use_TDF=0,use_FTC=0,use_DTG=1,use_BIC=0,use_EFV=0,use_DRV=1),
       events=make_ev(c("DRV","DTG"), 0, SIM_DAYS)),

  # 6. Delayed ART (treatment start at Day 180)
  list(id=6,
       name="⑥ Delayed ART (TDF/FTC/DTG start Day 180)",
       params=list(use_TDF=1,use_FTC=1,use_DTG=1,use_BIC=0,use_EFV=0,use_DRV=0),
       events=make_ev(c("TDF","FTC","DTG"), 180, SIM_DAYS)),

  # 7. Treatment Interruption (STI: off Day 365–730, restart Day 730)
  list(id=7,
       name="⑦ Structured Treatment Interruption (off 1y–2y)",
       params=list(use_TDF=1,use_FTC=1,use_DTG=1,use_BIC=0,use_EFV=0,use_DRV=0),
       events=c(make_ev(c("TDF","FTC","DTG"), 0, 365),
                make_ev(c("TDF","FTC","DTG"), 730, SIM_DAYS))),

  # 8. PrEP simulation (TDF/FTC, low viral challenge V0=10)
  list(id=8,
       name="⑧ PrEP (TDF/FTC — HIV Prevention)",
       params=list(use_TDF=1,use_FTC=1,use_DTG=0,use_BIC=0,use_EFV=0,use_DRV=0),
       events=make_ev(c("TDF","FTC"), 0, SIM_DAYS),
       init_override=list(V_FREE=0.01, I_CELL=1e-6, L_CELL=1e-9))
)

# ─────────────────────────────────────────────────────────────
# Simulation runner
# ─────────────────────────────────────────────────────────────
run_scenario <- function(mod, scen) {
  m <- mod %>% param(scen$params)

  if (!is.null(scen$init_override)) {
    m <- m %>% init(scen$init_override)
  }

  if (!is.null(scen$events)) {
    out <- mrgsim(m, events=scen$events, end=SIM_DAYS, delta=1, obsonly=TRUE)
  } else {
    out <- mrgsim(m, end=SIM_DAYS, delta=1, obsonly=TRUE)
  }

  as_tibble(out) %>% mutate(scenario_id=scen$id, scenario=scen$name)
}

# ─────────────────────────────────────────────────────────────
# Run all 8 scenarios
# ─────────────────────────────────────────────────────────────
message("Running 8 HIV/ART scenarios ...")
results <- lapply(scenarios, function(s) {
  tryCatch(run_scenario(mod, s),
           error = function(e) { message("Error in '", s$name, "': ", e$message); NULL })
})
results <- bind_rows(Filter(Negate(is.null), results))
message("Done. Rows: ", nrow(results))

# ─────────────────────────────────────────────────────────────
# Summary statistics
# ─────────────────────────────────────────────────────────────
summary_tbl <- results %>%
  group_by(scenario_id, scenario) %>%
  summarise(
    VL_w48   = round(VL_log10[which.min(abs(time - 336))], 2),
    CD4_w48  = round(CD4_count[which.min(abs(time - 336))]),
    CD4_nadir= round(min(CD4_count, na.rm=TRUE)),
    VL_supp_pct = round(mean(VL_supp, na.rm=TRUE)*100, 1),
    AIDS_events = max(AIDS_risk, na.rm=TRUE),
    Resist_max  = round(max(Res_score, na.rm=TRUE), 4),
    .groups = "drop"
  )

message("\n=== HIV/AIDS QSP Summary ===")
print(as.data.frame(summary_tbl))

# ─────────────────────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────────────────────

pal8 <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#999999")

# Plot 1: Viral Load (log10 copies/mL)
p1 <- ggplot(results, aes(x=time/365, y=VL_log10, color=scenario)) +
  geom_line(linewidth=0.8) +
  geom_hline(yintercept=log10(50), linetype="dashed", color="red", alpha=0.7) +
  annotate("text", x=4.8, y=log10(50)+0.2, label="VL < 50 c/mL",
           size=3, color="red", hjust=1) +
  scale_x_continuous("Time (Years)", breaks=0:5) +
  scale_y_continuous("Viral Load (log₁₀ copies/mL)", limits=c(0, 8)) +
  scale_color_manual(values=pal8) +
  labs(title="HIV Viral Load — 8 ART Scenarios", color=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", legend.text=element_text(size=7),
        legend.key.width=unit(0.5,"cm"))

# Plot 2: CD4 Count
p2 <- ggplot(results, aes(x=time/365, y=CD4_count, color=scenario)) +
  geom_line(linewidth=0.8) +
  geom_hline(yintercept=200, linetype="dashed", color="red",   alpha=0.6) +
  geom_hline(yintercept=500, linetype="dashed", color="green4",alpha=0.6) +
  annotate("text", x=4.8, y=215, label="CD4=200 (AIDS)", size=3, color="red",  hjust=1) +
  annotate("text", x=4.8, y=515, label="CD4=500 (normal)", size=3, color="green4", hjust=1) +
  scale_x_continuous("Time (Years)", breaks=0:5) +
  scale_y_continuous("CD4⁺ Count (cells/µL)", limits=c(0, 1300)) +
  scale_color_manual(values=pal8) +
  labs(title="CD4⁺ T Cell Dynamics", color=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", legend.text=element_text(size=7),
        legend.key.width=unit(0.5,"cm"))

# Plot 3: Latent Reservoir
p3 <- ggplot(results %>% filter(scenario_id > 1),
             aes(x=time/365, y=log10(pmax(Lat_IUPM, 1e-6)), color=scenario)) +
  geom_line(linewidth=0.8) +
  scale_x_continuous("Time (Years)", breaks=0:5) +
  scale_y_continuous("Latent Reservoir (log₁₀ IUPM)") +
  scale_color_manual(values=pal8[-1]) +
  labs(title="HIV Latent Reservoir Dynamics on ART",
       subtitle="IUPM = Infectious Units Per Million CD4⁺ T cells", color=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", legend.text=element_text(size=7),
        legend.key.width=unit(0.5,"cm"))

# Plot 4: Drug PK (scenario 2 — TDF/FTC/DTG, first 30 days)
pk_long <- results %>%
  filter(scenario_id == 2, time <= 30) %>%
  select(time, C_TFV, C_FTC, C_DTG) %>%
  pivot_longer(-time, names_to="Drug", values_to="Conc_ug_mL") %>%
  mutate(Drug = recode(Drug, C_TFV="TFV (plasma)", C_FTC="FTC (plasma)",
                       C_DTG="DTG (plasma)"))

p4 <- ggplot(pk_long, aes(x=time, y=Conc_ug_mL, color=Drug)) +
  geom_line(linewidth=0.9) +
  scale_x_continuous("Time (Days)", breaks=seq(0,30,5)) +
  scale_y_continuous("Plasma Conc. (µg/mL)") +
  scale_color_brewer(palette="Set1") +
  labs(title="TDF/FTC/DTG Plasma PK — Days 1–30", color=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom")

# Plot 5: Intracellular active metabolites
ic_long <- results %>%
  filter(scenario_id == 2, time <= 30) %>%
  select(time, TFV_DP_nM, FTC_TP_nM) %>%
  pivot_longer(-time, names_to="Metabolite", values_to="Conc_nM") %>%
  mutate(Metabolite = recode(Metabolite,
                             TFV_DP_nM = "TFV-DP (intracell.)",
                             FTC_TP_nM = "FTC-TP (intracell.)"))

ic50_lines <- data.frame(
  Metabolite = c("TFV-DP (intracell.)","FTC-TP (intracell.)"),
  IC50 = c(100, 1000)
)

p5 <- ggplot(ic_long, aes(x=time, y=Conc_nM, color=Metabolite)) +
  geom_line(linewidth=0.9) +
  geom_hline(data=ic50_lines, aes(yintercept=IC50, color=Metabolite),
             linetype="dashed", alpha=0.7) +
  scale_x_continuous("Time (Days)", breaks=seq(0,30,5)) +
  scale_y_continuous("Intracellular Conc. (nM)") +
  scale_color_manual(values=c("TFV-DP (intracell.)"="#E66101",
                               "FTC-TP (intracell.)"="#5E3C99")) +
  labs(title="Intracellular Active Metabolites: TFV-DP & FTC-TP",
       subtitle="Dashed = IC50", color=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom")

# Plot 6: Antiviral Efficacy (η)
p6 <- ggplot(results %>% filter(scenario_id > 1),
             aes(x=time/365, y=Eta*100, color=scenario)) +
  geom_line(linewidth=0.8) +
  scale_x_continuous("Time (Years)", breaks=0:5) +
  scale_y_continuous("Overall ART Efficacy η (%)", limits=c(0,101)) +
  scale_color_manual(values=pal8[-1]) +
  labs(title="Combined ART Antiviral Efficacy Over Time", color=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", legend.text=element_text(size=7),
        legend.key.width=unit(0.5,"cm"))

# Plot 7: CTL dynamics
p7 <- ggplot(results, aes(x=time/365, y=CTL_count, color=scenario)) +
  geom_line(linewidth=0.8) +
  scale_x_continuous("Time (Years)", breaks=0:5) +
  scale_y_continuous("CD8⁺ CTL Count (cells/µL)") +
  scale_color_manual(values=pal8) +
  labs(title="CD8⁺ CTL Immune Dynamics", color=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", legend.text=element_text(size=7),
        legend.key.width=unit(0.5,"cm"))

# Plot 8: Resistance Score
p8 <- ggplot(results, aes(x=time/365, y=Res_score*100, color=scenario)) +
  geom_line(linewidth=0.8) +
  scale_x_continuous("Time (Years)", breaks=0:5) +
  scale_y_continuous("Drug Resistance Score (%)", limits=c(0,100)) +
  scale_color_manual(values=pal8) +
  labs(title="Drug Resistance Score (Mutation Burden)", color=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", legend.text=element_text(size=7),
        legend.key.width=unit(0.5,"cm"))

# ─────────────────────────────────────────────────────────────
# Return
# ─────────────────────────────────────────────────────────────
invisible(list(
  model   = mod,
  results = results,
  summary = summary_tbl,
  plots   = list(p1, p2, p3, p4, p5, p6, p7, p8)
))
