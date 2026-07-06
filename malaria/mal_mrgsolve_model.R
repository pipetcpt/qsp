## ============================================================
## Malaria (Plasmodium falciparum / P. vivax) QSP Model
## mrgsolve ODE Simulation
## ============================================================
## Parasite life-cycle: age-structured ring→trophozoite→schizont
##   compartments (White NJ 1997 Parasitology; Simpson JA 2000
##   Br J Clin Pharmacol — parasite multiplication/clearance kinetics)
## Artemisinin PK/PD: Morris CA 2011 Malar J (AS/DHA population PK);
##   Saralamba S 2011 PNAS (ring-stage killing kinetics)
## Partner-drug PK: Tarning J 2012 (lumefantrine); Tarning J 2008 AAC
##   (piperaquine); Hietala SF 2007 AAC (amodiaquine/desethyl-AQ)
## Resistance: Ashley EA 2014 NEJM (K13/ring-stage survival assay);
##   WWARN Parasite Clearance Estimator methodology
## 22 compartments (13 PK + 9 disease/biomarker) · 8 treatment scenarios
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ─────────────────────────────────────────────────────────────
# mrgsolve MODEL CODE  (time unit = HOURS)
# ─────────────────────────────────────────────────────────────
code <- '
$PROB Malaria QSP Model (age-structured erythrocytic cycle + ACT PK/PD)

$PARAM
  // Drug-active switches (0/1)
  use_AS = 0, use_LUM = 0, use_PPQ = 0, use_AQ = 0, use_PQ = 0,
  K13_RES = 0,  // 0=wild-type, 1=K13 propeller mutant (artemisinin resistance)

  // ── Artesunate (AS) -> Dihydroartemisinin (DHA) PK (h-1, L) ──
  // Morris 2011 Malar J population PK (oral/IV artesunate)
  ka_AS     = 4.0,     // fast oral absorption (h-1)
  ke_AS     = 1.386,   // AS elimination / hydrolysis to DHA (t1/2=0.5h)
  Vc_AS     = 30.0,    // L
  ke_DHA    = 0.90,    // DHA elimination (h-1; t1/2~0.77h)
  Vc_DHA    = 70.0,    // L
  frac_AS_DHA = 0.90,  // fraction of AS converted to DHA

  // ── Lumefantrine (LUM) PK — Tarning 2012 CPT ──
  ka_LUM    = 0.30,    // slow, food-dependent absorption (h-1)
  CL_LUM    = 0.14,    // L/h (very low CL, long t1/2 3-6d)
  Vc_LUM    = 21.0,    // L
  Q_LUM     = 0.05,    // L/h inter-compartmental
  Vp_LUM    = 500.0,   // L (large peripheral -> long terminal phase)

  // ── Piperaquine (PPQ) PK — Tarning 2008 AAC ──
  ka_PPQ    = 0.10,    // h-1
  CL_PPQ    = 0.90,    // L/h
  Vc_PPQ    = 150.0,   // L
  Q_PPQ     = 0.60,    // L/h
  Vp_PPQ    = 5000.0,  // L (huge Vd -> t1/2 20-28 days)

  // ── Amodiaquine -> Desethylamodiaquine (DEAQ) PK — Hietala 2007 AAC ──
  ka_AQ     = 0.50,    // h-1
  ke_DEAQ   = 0.0032,  // h-1 (t1/2 ~9 days)
  Vc_DEAQ   = 1500.0,  // L

  // ── Primaquine (PQ) PK ──
  ka_PQ     = 0.60,    // h-1
  ke_PQ     = 0.11,    // h-1 (t1/2 ~6.3h; drives daily dosing)
  Vc_PQ     = 220.0,   // L

  // ── Drug PD: Emax/Hill kill parameters ──
  // Artemisinin: dominant, broad (ring+troph) killer; PRR~1e4/cycle
  Emax_AS    = 0.42,   // max additional log-kill rate (h-1) at high DHA
  EC50_AS    = 5.0,    // ng/mL (approx in vitro IC50-equivalent)
  hill_AS    = 1.8,
  Kres_shift = 6.0,    // EC50 multiplier under K13 resistance (ring-stage survival)

  // Partner drugs: slower onset, broad troph/schizont kill + prophylaxis
  Emax_LUM = 0.10, EC50_LUM = 280.0, hill_LUM = 1.0,   // ng/mL
  Emax_PPQ = 0.10, EC50_PPQ = 30.0,  hill_PPQ = 1.0,   // ng/mL
  Emax_AQ  = 0.10, EC50_AQ  = 40.0,  hill_AQ  = 1.0,   // ng/mL (as DEAQ)

  // Primaquine: hypnozoiticidal + gametocytocidal (CYP2D6-dependent, simplified)
  Emax_PQ_liver = 0.06, EC50_PQ = 15.0, hill_PQ = 1.0,
  Emax_PQ_gam   = 0.30,
  G6PD_deficient = 0,    // 0=normal, 1=deficient (oxidative hemolysis risk)
  k_PQ_hemolysis = 0.0006,

  // ── Parasite life-cycle kinetics (P. falciparum, 48h cycle) ──
  k_RT      = 0.0556,  // ring->troph (1/18h)
  k_TS      = 0.0833,  // troph->schizont (1/12h)
  k_SR      = 0.0556,  // schizont->rupture (1/18h)
  burst     = 16.0,    // merozoites released per rupture
  inv_eff0  = 0.45,    // baseline successful invasion fraction
  immune_block = 0.85, // max fractional reduction of invasion by acquired immunity
  gam_commit   = 0.015,// fraction of rings committing to gametocytogenesis
  k_gam_mat    = 0.0052,// gametocyte maturation+loss (h-1; ~8d lifespan)

  // Liver stage / hypnozoite (P. vivax relapse) kinetics
  k_liver_egress = 0.0060, // liver schizont -> blood merozoites (h-1; ~7d prepatent)
  relapse_on     = 0,      // 1 = P. vivax relapse-prone strain
  k_relapse      = 0.00028,// hypnozoite reactivation source (h-1; weeks-months)

  // RBC / hemoglobin dynamics
  RBC_base   = 5.0e6,  // cells/uL
  k_RBC_sen  = 3.47e-4,// senescence (h-1; 120-day lifespan)
  k_ery_comp = 0.0020, // compensatory erythropoiesis (h-1)
  HB_base    = 13.5,   // g/dL
  hb_per_rupture = 2.5e-6, // g/dL lost per (schizont rupture)/uL/h
  k_dysery   = 3.0e-6, // dyserythropoiesis suppression per parasite/uL

  // Fever / cytokine proxy
  Tbase      = 37.0,
  Tmax_rise  = 3.2,
  K_fever    = 2.0e4,  // parasites/uL for half-max fever response
  k_temp_decay = 0.25, // h-1

  // Acquired immunity
  k_immune_gain  = 3.0e-7,
  k_immune_decay = 0.0010,
  K_immune       = 5.0e4  // parasites/uL for half-max immune stimulation

$INIT
  // Drug PK (13 compartments)
  AS_GUT=0, AS_PLASMA=0, DHA_PLASMA=0,
  LUM_GUT=0, LUM_CENTRAL=0, LUM_PERIPH=0,
  PPQ_GUT=0, PPQ_CENTRAL=0, PPQ_PERIPH=0,
  AQ_GUT=0, DEAQ_PLASMA=0,
  PQ_GUT=0, PQ_PLASMA=0,
  // Disease / biomarkers (9 compartments)
  RBC_U=5.0e6, PRBC_RING=10.0, PRBC_TROPH=6.0, PRBC_SCHIZONT=4.0,
  LIVER_PARASITE=1.0e4, GAMETOCYTE=0.1,
  HB=13.5, TEMP=37.0, IMMUNITY=0.05

$ODE
  // ── PK ──────────────────────────────────────────────────────
  dxdt_AS_GUT     = -ka_AS * AS_GUT;
  dxdt_AS_PLASMA  =  ka_AS * AS_GUT / Vc_AS - ke_AS * AS_PLASMA;
  dxdt_DHA_PLASMA =  frac_AS_DHA * ke_AS * AS_PLASMA * Vc_AS / Vc_DHA - ke_DHA * DHA_PLASMA;

  // 2-compartment PK expressed directly in concentrations (mg/L):
  // dC1/dt = in - k10*C1 - (Q/V1)*(C1-C2);  dC2/dt = (Q/V2)*(C1-C2)
  double keLUM = CL_LUM / Vc_LUM;
  double kcpLUM = Q_LUM / Vc_LUM;
  double kpcLUM = Q_LUM / Vp_LUM;
  dxdt_LUM_GUT     = -ka_LUM * LUM_GUT;
  dxdt_LUM_CENTRAL =  ka_LUM * LUM_GUT / Vc_LUM - keLUM * LUM_CENTRAL - kcpLUM * (LUM_CENTRAL - LUM_PERIPH);
  dxdt_LUM_PERIPH  =  kpcLUM * (LUM_CENTRAL - LUM_PERIPH);

  double kePPQ = CL_PPQ / Vc_PPQ;
  double kcpPPQ = Q_PPQ / Vc_PPQ;
  double kpcPPQ = Q_PPQ / Vp_PPQ;
  dxdt_PPQ_GUT     = -ka_PPQ * PPQ_GUT;
  dxdt_PPQ_CENTRAL =  ka_PPQ * PPQ_GUT / Vc_PPQ - kePPQ * PPQ_CENTRAL - kcpPPQ * (PPQ_CENTRAL - PPQ_PERIPH);
  dxdt_PPQ_PERIPH  =  kpcPPQ * (PPQ_CENTRAL - PPQ_PERIPH);

  dxdt_AQ_GUT      = -ka_AQ * AQ_GUT;
  dxdt_DEAQ_PLASMA =  ka_AQ * AQ_GUT / Vc_DEAQ - ke_DEAQ * DEAQ_PLASMA;

  dxdt_PQ_GUT     = -ka_PQ * PQ_GUT;
  dxdt_PQ_PLASMA  =  ka_PQ * PQ_GUT / Vc_PQ - ke_PQ * PQ_PLASMA;

  // ── PD: Hill/Emax drug effects (concentrations in ng/mL) ──────
  // All *_PLASMA / *_CENTRAL states are concentrations already (mg/L); ng/mL = mg/L * 1000
  double C_DHA   = (DHA_PLASMA > 0) ? DHA_PLASMA * 1000.0 : 0.0;
  double Clum_ng = (LUM_CENTRAL > 0) ? LUM_CENTRAL * 1000.0 : 0.0;
  double Cppq_ng = (PPQ_CENTRAL > 0) ? PPQ_CENTRAL * 1000.0 : 0.0;
  double Cdeaq_ng= (DEAQ_PLASMA > 0) ? DEAQ_PLASMA * 1000.0 : 0.0;
  double Cpq_ng  = (PQ_PLASMA   > 0) ? PQ_PLASMA   * 1000.0 : 0.0;

  double EC50_AS_eff = EC50_AS * (1.0 + (Kres_shift - 1.0) * K13_RES);
  double E_AS  = use_AS  * Emax_AS  * pow(C_DHA,  hill_AS)  / (pow(EC50_AS_eff, hill_AS)  + pow(C_DHA,  hill_AS)+1e-9);
  double E_LUM = use_LUM * Emax_LUM * pow(Clum_ng, hill_LUM) / (pow(EC50_LUM, hill_LUM) + pow(Clum_ng, hill_LUM)+1e-9);
  double E_PPQ = use_PPQ * Emax_PPQ * pow(Cppq_ng, hill_PPQ) / (pow(EC50_PPQ, hill_PPQ) + pow(Cppq_ng, hill_PPQ)+1e-9);
  double E_AQ  = use_AQ  * Emax_AQ  * pow(Cdeaq_ng,hill_AQ)  / (pow(EC50_AQ,  hill_AQ)  + pow(Cdeaq_ng,hill_AQ)+1e-9);
  double E_partner = 1.0 - (1.0 - E_LUM) * (1.0 - E_PPQ) * (1.0 - E_AQ); // independent-action combo
  double E_PQ_liver = use_PQ * Emax_PQ_liver * pow(Cpq_ng, hill_PQ) / (pow(EC50_PQ, hill_PQ) + pow(Cpq_ng, hill_PQ)+1e-9);
  double E_PQ_gam   = use_PQ * Emax_PQ_gam   * pow(Cpq_ng, hill_PQ) / (pow(EC50_PQ, hill_PQ) + pow(Cpq_ng, hill_PQ)+1e-9);

  // ── Parasite life cycle (age-structured, guarded non-negative) ──
  double RU = (RBC_U > 0) ? RBC_U : 0.0;
  double R  = (PRBC_RING > 0) ? PRBC_RING : 0.0;
  double T  = (PRBC_TROPH > 0) ? PRBC_TROPH : 0.0;
  double S  = (PRBC_SCHIZONT > 0) ? PRBC_SCHIZONT : 0.0;
  double LP = (LIVER_PARASITE > 0) ? LIVER_PARASITE : 0.0;
  double IMM = (IMMUNITY > 0) ? IMMUNITY : 0.0;

  double total_para = R + T + S;
  double inv_eff = inv_eff0 * (1.0 - immune_block * IMM);
  double liver_flux = k_liver_egress * LP;                 // primary liver schizogony -> blood
  double relapse_source = relapse_on * k_relapse * 5.0e3; // continuous slow hypnozoite reactivation source
  double rupture_flux = k_SR * S;                           // parasites rupturing per hour
  double new_infections = (liver_flux + rupture_flux * burst) * inv_eff * (RU / RBC_base);

  // Drug kill terms (additional first-order loss, h-1)
  double kill_ring  = E_AS;                    // artemisinin: dominant ring-stage action
  double kill_troph = E_AS * 0.6 + E_partner;  // partner drugs: broad troph/schizont action
  double kill_schiz = E_partner * 0.8;

  dxdt_LIVER_PARASITE = relapse_source - liver_flux - E_PQ_liver * LP;

  dxdt_RBC_U = k_RBC_sen * (RBC_base - RU) + k_ery_comp * (HB_base - HB > 0 ? (HB_base - HB) : 0.0) * 1.0e5
               / (1.0 + k_dysery * total_para)
               - new_infections;

  dxdt_PRBC_RING     = new_infections - k_RT * R - kill_ring  * R;
  dxdt_PRBC_TROPH    = k_RT * R - k_TS * T - kill_troph * T - gam_commit * k_RT * R;
  dxdt_PRBC_SCHIZONT = k_TS * T - k_SR * S - kill_schiz * S;

  dxdt_GAMETOCYTE = gam_commit * k_RT * R - k_gam_mat * GAMETOCYTE - E_PQ_gam * GAMETOCYTE;

  // ── Hemoglobin / fever / immunity ─────────────────────────────
  // Hemolysis proportional to schizont rupture events (1 host RBC destroyed per rupture)
  double hemolysis = hb_per_rupture * rupture_flux
                    + (G6PD_deficient > 0.5 ? k_PQ_hemolysis * Cpq_ng * (HB/HB_base) : 0.0);
  dxdt_HB = -hemolysis + 0.0015 * (HB_base - HB) / (1.0 + k_dysery * total_para);

  double fever_drive = Tmax_rise * rupture_flux * burst / (rupture_flux * burst + K_fever + 1e-6);
  dxdt_TEMP = (Tbase + fever_drive - TEMP) * k_temp_decay;

  dxdt_IMMUNITY = k_immune_gain * total_para / (total_para + K_immune) * (1.0 - IMM) - k_immune_decay * IMM;

$TABLE
  double C_DHA_out  = (DHA_PLASMA  > 0) ? DHA_PLASMA  * 1000.0 : 0.0;  // mg/L -> ng/mL
  double C_LUM_out  = (LUM_CENTRAL > 0) ? LUM_CENTRAL * 1000.0 : 0.0;
  double C_PPQ_out  = (PPQ_CENTRAL > 0) ? PPQ_CENTRAL * 1000.0 : 0.0;
  double C_DEAQ_out = (DEAQ_PLASMA > 0) ? DEAQ_PLASMA * 1000.0 : 0.0;
  double C_PQ_out   = (PQ_PLASMA   > 0) ? PQ_PLASMA   * 1000.0 : 0.0;

  double PERIPH_PARASITEMIA = PRBC_RING + 0.10 * PRBC_TROPH;  // visible on smear (sequestration)
  double TOTAL_PARASITEMIA  = PRBC_RING + PRBC_TROPH + PRBC_SCHIZONT;
  double log10_para = (PERIPH_PARASITEMIA > 0.05) ? log10(PERIPH_PARASITEMIA) : -1.301;
  double Hb_gdl     = HB;
  double Temp_C     = TEMP;
  double Gam_uL      = GAMETOCYTE;
  double Immunity_pct= IMMUNITY * 100.0;
  double Liver_burden= LIVER_PARASITE;

$CAPTURE
  C_DHA_out C_LUM_out C_PPQ_out C_DEAQ_out C_PQ_out
  PERIPH_PARASITEMIA TOTAL_PARASITEMIA log10_para
  Hb_gdl Temp_C Gam_uL Immunity_pct Liver_burden
'

# ─────────────────────────────────────────────────────────────
# Compile
# ─────────────────────────────────────────────────────────────
mod <- mrgsolve::mcode("malaria_qsp", code)

# ─────────────────────────────────────────────────────────────
# Dosing event builder (70 kg adult; doses in mg -> compartment "amount")
# ─────────────────────────────────────────────────────────────
make_ev <- function(regimen, t_start) {
  switch(regimen,
    # Artemether-lumefantrine: here approximated with artesunate-equivalent
    # DHA exposure + lumefantrine twice daily x3d (WHO 6-dose regimen)
    AL = c(
      ev(cmt="AS_GUT",  amt=100, ii=12, addl=5, time=t_start),
      ev(cmt="LUM_GUT", amt=480, ii=12, addl=5, time=t_start)
    ),
    ASAQ = c(
      ev(cmt="AS_GUT", amt=200, ii=24, addl=2, time=t_start),
      ev(cmt="AQ_GUT", amt=600, ii=24, addl=2, time=t_start)
    ),
    DP = c(
      ev(cmt="AS_GUT",  amt=240, ii=24, addl=2, time=t_start),
      ev(cmt="PPQ_GUT", amt=960, ii=24, addl=2, time=t_start)
    ),
    IV_AS = c(
      # 2.4 mg/kg IV bolus, 70kg adult; AS_PLASMA is a concentration compartment
      # (Vc_AS = 30 L), so bolus amount is dose_mg / Vc_AS -> mg/L bump
      ev(cmt="AS_PLASMA", amt=2.4*70/30.0, ii=12, addl=2, time=t_start)
    ),
    PQ14 = ev(cmt="PQ_GUT", amt=30, ii=24, addl=13, time=t_start),
    NONE = NULL
  )
}

SIM_H <- 40 * 24   # 40 days for standard treatment scenarios
SIM_H_VIVAX <- 200 * 24  # 200 days to show relapse biology

# ─────────────────────────────────────────────────────────────
# Scenario definitions (>=5 required; 8 provided)
# ─────────────────────────────────────────────────────────────
scenarios <- list(
  list(id=1, name="① Untreated (Non-Immune, Natural History)",
       params=list(use_AS=0,use_LUM=0,use_PPQ=0,use_AQ=0,use_PQ=0,K13_RES=0),
       events=NULL, end=14*24, init=list(IMMUNITY=0.02)),

  list(id=2, name="② Untreated (Semi-Immune, High Transmission Setting)",
       params=list(use_AS=0,use_LUM=0,use_PPQ=0,use_AQ=0,use_PQ=0,K13_RES=0),
       events=NULL, end=14*24, init=list(IMMUNITY=0.55)),

  list(id=3, name="③ Artemether-Lumefantrine (AL, 3-day ACT)",
       params=list(use_AS=1,use_LUM=1,use_PPQ=0,use_AQ=0,use_PQ=0,K13_RES=0),
       events=make_ev("AL", 0), end=SIM_H, init=list()),

  list(id=4, name="④ Artesunate-Amodiaquine (ASAQ, 3-day ACT)",
       params=list(use_AS=1,use_LUM=0,use_PPQ=0,use_AQ=1,use_PQ=0,K13_RES=0),
       events=make_ev("ASAQ", 0), end=SIM_H, init=list()),

  list(id=5, name="⑤ Dihydroartemisinin-Piperaquine (DP, 3-day ACT)",
       params=list(use_AS=1,use_LUM=0,use_PPQ=1,use_AQ=0,use_PQ=0,K13_RES=0),
       events=make_ev("DP", 0), end=SIM_H, init=list()),

  list(id=6, name="⑥ Severe Malaria: IV Artesunate → Oral AL Follow-On",
       params=list(use_AS=1,use_LUM=1,use_PPQ=0,use_AQ=0,use_PQ=0,K13_RES=0),
       events=c(make_ev("IV_AS", 0), make_ev("AL", 72)),
       end=SIM_H, init=list(PRBC_RING=400, PRBC_TROPH=250, PRBC_SCHIZONT=150, HB=9.0)),

  list(id=7, name="⑦ Artemisinin-Resistant (K13 C580Y) + Standard AL",
       params=list(use_AS=1,use_LUM=1,use_PPQ=0,use_AQ=0,use_PQ=0,K13_RES=1),
       events=make_ev("AL", 0), end=SIM_H, init=list()),

  list(id=8, name="⑧ P. vivax: ACT + 14-day Primaquine Radical Cure",
       params=list(use_AS=1,use_LUM=1,use_PPQ=0,use_AQ=0,use_PQ=1,K13_RES=0,relapse_on=1),
       events=c(make_ev("AL", 0), make_ev("PQ14", 0)),
       end=SIM_H_VIVAX, init=list()),

  list(id=9, name="⑨ P. vivax: ACT WITHOUT Radical Cure (Relapse)",
       params=list(use_AS=1,use_LUM=1,use_PPQ=0,use_AQ=0,use_PQ=0,K13_RES=0,relapse_on=1),
       events=make_ev("AL", 0), end=SIM_H_VIVAX, init=list())
)

# ─────────────────────────────────────────────────────────────
# Simulation runner
# ─────────────────────────────────────────────────────────────
run_scenario <- function(mod, scen) {
  m <- mod %>% param(scen$params)
  if (length(scen$init) > 0) m <- m %>% init(scen$init)
  if (!is.null(scen$events)) {
    out <- mrgsim(m, events=scen$events, end=scen$end, delta=1, obsonly=TRUE, maxsteps=100000)
  } else {
    out <- mrgsim(m, end=scen$end, delta=1, obsonly=TRUE, maxsteps=100000)
  }
  as_tibble(out) %>% mutate(scenario_id=scen$id, scenario=scen$name)
}

message("Running 9 malaria treatment/natural-history scenarios ...")
results <- bind_rows(lapply(scenarios, function(s) {
  tryCatch(run_scenario(mod, s), error=function(e) { message("Error in '", s$name, "': ", e$message); NULL })
}))
message("Done. Rows: ", nrow(results))

# ─────────────────────────────────────────────────────────────
# Clinical/WWARN-style summary endpoints
# ─────────────────────────────────────────────────────────────
summary_tbl <- results %>%
  group_by(scenario_id, scenario) %>%
  summarise(
    Peak_parasitemia   = round(max(PERIPH_PARASITEMIA, na.rm=TRUE)),
    PCT_hours          = { d <- PERIPH_PARASITEMIA; tt <- time
                            below <- tt[which(d < 10)]
                            if (length(below) > 0) round(min(below[below > 12])) else NA },
    Day3_positive      = ifelse(any(abs(time-72) < 1 & PERIPH_PARASITEMIA >= 10), 1, 0),
    Hb_nadir           = round(min(Hb_gdl, na.rm=TRUE), 1),
    Tmax_C             = round(max(Temp_C, na.rm=TRUE), 1),
    Day28_recrud       = ifelse(any(time > 28*24 & PERIPH_PARASITEMIA >= 10), 1, 0),
    .groups="drop"
  )

message("\n=== Malaria QSP Summary (WWARN-style endpoints) ===")
print(as.data.frame(summary_tbl))

# ─────────────────────────────────────────────────────────────
# Plots
# ─────────────────────────────────────────────────────────────
pal9 <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF","#999999","#66C2A5")

p1 <- ggplot(results, aes(x=time/24, y=log10_para, color=scenario)) +
  geom_line(linewidth=0.8) +
  geom_hline(yintercept=log10(10), linetype="dashed", color="red", alpha=0.6) +
  scale_x_continuous("Time (Days)") +
  scale_y_continuous("Peripheral Parasitemia (log10 parasites/µL)") +
  scale_color_manual(values=pal9) +
  labs(title="Peripheral Parasitemia — 9 Scenarios", color=NULL) +
  theme_bw(base_size=11) + theme(legend.position="bottom", legend.text=element_text(size=7))

p2 <- ggplot(results %>% filter(scenario_id %in% c(3,4,5), time<=96),
             aes(x=time, y=C_DHA_out, color=scenario)) +
  geom_line(linewidth=0.9) +
  scale_x_continuous("Time (Hours)") + scale_y_continuous("DHA Plasma (ng/mL)") +
  labs(title="Dihydroartemisinin (Active Metabolite) PK — First 4 Days", color=NULL) +
  theme_bw(base_size=11) + theme(legend.position="bottom")

p3 <- ggplot(results %>% filter(scenario_id==3),
             aes(x=time/24, y=C_LUM_out)) +
  geom_line(linewidth=0.9, color="#984EA3") +
  scale_x_continuous("Time (Days)") + scale_y_continuous("Lumefantrine Plasma (ng/mL)") +
  labs(title="Lumefantrine PK — Long Terminal Phase (AL Regimen)") +
  theme_bw(base_size=11)

p4 <- ggplot(results, aes(x=time/24, y=Hb_gdl, color=scenario)) +
  geom_line(linewidth=0.8) +
  scale_x_continuous("Time (Days)") + scale_y_continuous("Hemoglobin (g/dL)") +
  scale_color_manual(values=pal9) +
  labs(title="Hemoglobin Trajectory (Malarial Anemia & Recovery)", color=NULL) +
  theme_bw(base_size=11) + theme(legend.position="bottom", legend.text=element_text(size=7))

p5 <- ggplot(results, aes(x=time/24, y=Temp_C, color=scenario)) +
  geom_line(linewidth=0.8) +
  geom_hline(yintercept=37.5, linetype="dashed", color="grey40") +
  scale_x_continuous("Time (Days)") + scale_y_continuous("Body Temperature (°C)") +
  scale_color_manual(values=pal9) +
  labs(title="Fever Course (Paroxysmal, Rupture-Driven)", color=NULL) +
  theme_bw(base_size=11) + theme(legend.position="bottom", legend.text=element_text(size=7))

p6 <- ggplot(results %>% filter(scenario_id %in% c(8,9)),
             aes(x=time/24, y=Liver_burden, color=scenario)) +
  geom_line(linewidth=0.9) +
  scale_x_continuous("Time (Days)") + scale_y_continuous("Liver-Stage/Hypnozoite Burden") +
  labs(title="P. vivax Hypnozoite Reservoir: Radical Cure vs No Radical Cure", color=NULL) +
  theme_bw(base_size=11) + theme(legend.position="bottom")

p7 <- ggplot(results %>% filter(scenario_id %in% c(8,9)),
             aes(x=time/24, y=log10_para, color=scenario)) +
  geom_line(linewidth=0.9) +
  scale_x_continuous("Time (Days)") + scale_y_continuous("Parasitemia (log10/µL)") +
  labs(title="P. vivax Relapse Pattern (Radical Cure vs Not)", color=NULL) +
  theme_bw(base_size=11) + theme(legend.position="bottom")

p8 <- ggplot(results, aes(x=time/24, y=Gam_uL, color=scenario)) +
  geom_line(linewidth=0.8) +
  scale_x_continuous("Time (Days)") + scale_y_continuous("Mature Gametocyte Density (/µL)") +
  scale_color_manual(values=pal9) +
  labs(title="Gametocyte Carriage (Transmission Potential)", color=NULL) +
  theme_bw(base_size=11) + theme(legend.position="bottom", legend.text=element_text(size=7))

p9 <- ggplot(results %>% filter(scenario_id %in% c(3,7)),
             aes(x=time/24, y=log10_para, color=scenario)) +
  geom_line(linewidth=1.0) +
  scale_x_continuous("Time (Days)", limits=c(0,10)) + scale_y_continuous("Parasitemia (log10/µL)") +
  scale_color_manual(values=c("#4DAF4A","#F781BF")) +
  labs(title="Artemisinin Resistance: Delayed Parasite Clearance (K13 Mutant)", color=NULL) +
  theme_bw(base_size=11) + theme(legend.position="bottom")

# ─────────────────────────────────────────────────────────────
# Return
# ─────────────────────────────────────────────────────────────
invisible(list(
  model   = mod,
  results = results,
  summary = summary_tbl,
  plots   = list(p1,p2,p3,p4,p5,p6,p7,p8,p9)
))
