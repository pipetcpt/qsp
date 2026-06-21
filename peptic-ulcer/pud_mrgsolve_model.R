## ============================================================
## Peptic Ulcer Disease (PUD) – QSP Model in mrgsolve
## ============================================================
## Disease: Peptic Ulcer Disease (소화성 궤양)
## Pathophysiology: Imbalance between aggressive factors (acid, H. pylori,
##   NSAIDs) and mucosal defense (mucus, bicarbonate, prostaglandins)
## Treatments modeled: PPI (omeprazole), H2RA (famotidine),
##   H. pylori eradication triple therapy (PPI + amoxicillin + clarithromycin),
##   NSAIDs (ibuprofen), misoprostol (PGE1 analog)
##
## Parameter calibration notes:
##   - PPI PK: Shin JM et al. Pharmacol Ther 2009
##   - H. pylori eradication rates: Graham DY, Yamaoka Y Curr Opin Gastroenterol 2019
##   - Acid secretion model: Feldman M, et al. Gastroenterology 1984
##   - Ulcer healing rates: Laine L, et al. NEJM 2008
##   - NSAID gastropathy model: Scheiman JM Ann Int Med 1995
##   - Mucosal defense: Wallace JL Physiol Rev 2008
## ============================================================

library(mrgsolve)

code <- '
$PROB
Peptic Ulcer Disease (PUD) QSP Model
Version 1.0 - 2026-06-19
Compartments: PPI PK (3), H2RA PK (2), Amoxicillin PK (2),
              Clarithromycin PK (2), NSAID PK (2),
              H. pylori dynamics (1), Acid secretion (2),
              Mucosal defense (2), Inflammation (2),
              Ulcer dynamics (2) = 20 ODEs total

$PARAM
// ── PPI (Omeprazole) PK Parameters ──────────────────────────
F_PPI       = 0.65   // oral bioavailability (pH-dependent)
ka_PPI      = 1.2    // absorption rate (/hr)
Vc_PPI      = 14.0   // central volume (L)
Vp_PPI      = 8.0    // peripheral volume (L)
CL_PPI      = 28.0   // total clearance (L/hr) - CYP2C19 dependent
Q_PPI       = 6.0    // intercompartmental clearance (L/hr)
EC50_PPI    = 0.3    // EC50 for pump inhibition (mg/L)
Emax_PPI    = 0.98   // maximum pump inhibition fraction

// ── H2RA (Famotidine) PK Parameters ─────────────────────────
F_H2RA      = 0.45   // oral bioavailability
ka_H2RA     = 0.8    // absorption rate (/hr)
Vc_H2RA     = 25.0   // central volume (L)
CL_H2RA     = 10.0   // renal clearance (L/hr)
EC50_H2RA   = 0.05   // EC50 for H2R inhibition (mg/L)
Emax_H2RA   = 0.90   // maximum H2R inhibition

// ── Amoxicillin (AMX) PK Parameters ─────────────────────────
F_AMX       = 0.85   // oral bioavailability
ka_AMX      = 1.5    // absorption rate (/hr)
Vc_AMX      = 22.0   // central volume (L)
CL_AMX      = 18.0   // renal clearance (L/hr)
MIC_AMX     = 0.125  // MIC against H. pylori (mg/L)
Emax_AMX    = 0.85   // max kill fraction of H. pylori

// ── Clarithromycin (CLR) PK Parameters ──────────────────────
F_CLR       = 0.55   // oral bioavailability
ka_CLR      = 1.0    // absorption rate (/hr)
Vc_CLR      = 72.0   // central volume (L)
CL_CLR      = 30.0   // CYP3A4-mediated clearance (L/hr)
MIC_CLR     = 0.25   // MIC against H. pylori (mg/L)
Emax_CLR    = 0.90   // max kill fraction
ResistanceCLR = 0.15 // baseline resistance rate (15%)

// ── NSAID (Ibuprofen) PK Parameters ─────────────────────────
F_NSAID     = 0.80   // oral bioavailability
ka_NSAID    = 2.0    // absorption rate (/hr)
Vc_NSAID    = 8.0    // central volume (L)
CL_NSAID    = 4.0    // hepatic clearance (L/hr)
IC50_COX1   = 0.5    // IC50 for COX-1 inhibition (mg/L)
IC50_COX2   = 2.0    // IC50 for COX-2 inhibition (mg/L) - ibuprofen
Imax_COX    = 0.95   // max COX inhibition

// ── H. pylori Dynamics ───────────────────────────────────────
HP0         = 1e7    // initial H. pylori load (CFU/mL)
k_grow_HP   = 0.02   // H. pylori growth rate (/hr)
K_HP        = 1e9    // carrying capacity (CFU/mL)
k_clear_HP  = 0.001  // natural immune clearance (/hr)
pH_kill_HP  = 3.0    // pH below which HP killed (pH unit threshold)
k_pH_kill   = 0.05   // rate constant for acid killing of HP

// ── Acid Secretion Parameters ────────────────────────────────
MaxAcidOutput = 40.0  // maximal acid output (mEq/hr)
BasalAcid   = 4.0    // basal acid output (mEq/hr)
k_pumpTurn  = 0.12   // proton pump turnover/synthesis rate (/hr)
pHbasal     = 1.5    // basal intragastric pH
pHmax       = 7.0    // maximum intragastric pH
HP_acid_stim = 0.3   // H. pylori stimulation of acid (fraction)

// ── Mucosal Defense Parameters ───────────────────────────────
MucusMax    = 1.0    // maximum mucus layer (normalized)
k_mucusProd = 0.10   // mucus production rate (/hr)
k_mucusDeg  = 0.05   // mucus degradation rate (/hr)
PG_protect  = 0.8    // prostaglandin protection (fraction)
k_PG_base   = 0.10   // basal prostaglandin synthesis rate (/hr)
k_PG_deg    = 0.08   // prostaglandin degradation rate (/hr)
MucusIC50   = 4.0    // pH IC50 for mucus disruption

// ── Inflammation Parameters ──────────────────────────────────
k_inflam    = 0.05   // HP-driven inflammation onset rate (/hr)
k_inflam_res = 0.02  // inflammation resolution rate (/hr)
InflamMax   = 10.0   // maximum inflammation score (0-10)
IL8_HP      = 0.8    // HP contribution to IL-8/neutrophil infiltration

// ── Ulcer Dynamics ───────────────────────────────────────────
k_damage    = 0.15   // mucosal damage rate (per unit acid*inflammation)
k_heal      = 0.08   // baseline healing rate (/hr)
k_recur     = 0.005  // ulcer recurrence risk (/hr, HP-dependent)
UlcerMax    = 100.0  // max ulcer area (mm²)
EGF_heal    = 0.3    // EGF contribution to healing (fraction boost)

$CMT
// PPI compartments
PPI_GUT PPI_CENT PPI_PERI

// H2RA compartments
H2RA_GUT H2RA_CENT

// Amoxicillin compartments
AMX_GUT AMX_CENT

// Clarithromycin compartments
CLR_GUT CLR_CENT

// NSAID compartments
NSAID_GUT NSAID_CENT

// H. pylori load (log10 scale state)
HP_LOAD

// Acid secretion: active pump fraction and intragastric pH
PUMP_ACTIVE pH_INTRAGASTRIC

// Mucosal defense: mucus layer and prostaglandin level
MUCUS PG_LEVEL

// Inflammation score and ulcer area
INFLAM ULCER_AREA

$MAIN
// ── Derived PK quantities ────────────────────────────────────
double C_PPI    = PPI_CENT  / Vc_PPI;
double C_H2RA   = H2RA_CENT / Vc_H2RA;
double C_AMX    = AMX_CENT  / Vc_AMX;
double C_CLR    = CLR_CENT  / Vc_CLR;
double C_NSAID  = NSAID_CENT / Vc_NSAID;

// ── PPI Pharmacodynamics ─────────────────────────────────────
// Acid-activated PPI binds proton pumps covalently → irreversible inhibition
// New pumps synthesized at rate k_pumpTurn; inhibited fraction:
double inh_PPI  = Emax_PPI * C_PPI / (EC50_PPI + C_PPI);

// ── H2RA Pharmacodynamics ────────────────────────────────────
double inh_H2RA = Emax_H2RA * C_H2RA / (EC50_H2RA + C_H2RA);

// Combined acid suppression (PPI dominant due to irreversible binding)
double total_acid_inh = 1.0 - (1.0 - inh_PPI) * (1.0 - inh_H2RA);

// ── NSAID COX inhibition ─────────────────────────────────────
double inh_COX1 = Imax_COX * C_NSAID / (IC50_COX1 + C_NSAID);
double inh_COX2 = Imax_COX * C_NSAID / (IC50_COX2 + C_NSAID);

// Prostaglandin reduction due to COX-1 inhibition
double PG_NSAID_factor = 1.0 - 0.7 * inh_COX1; // COX-1 is main PG source

// ── H. pylori killing by antibiotics ─────────────────────────
double kill_AMX = (C_AMX > MIC_AMX) ? Emax_AMX * (C_AMX - MIC_AMX) / (C_AMX - MIC_AMX + MIC_AMX) : 0.0;
double kill_CLR = (C_CLR > MIC_CLR && ResistanceCLR < 0.5) ?
                  Emax_CLR * (C_CLR - MIC_CLR) / (C_CLR - MIC_CLR + MIC_CLR) : 0.0;
double kill_combo = 1.0 - (1.0 - kill_AMX) * (1.0 - kill_CLR);

// ── Current intragastric pH-based HP suppression ─────────────
double pH = pH_INTRAGASTRIC;
double acid_HP_kill = (pH < pH_kill_HP) ? k_pH_kill * (pH_kill_HP - pH) : 0.0;

$ODE
// ── PPI PK (3-compartment oral) ──────────────────────────────
dxdt_PPI_GUT  = -ka_PPI * PPI_GUT;
dxdt_PPI_CENT =  ka_PPI * F_PPI * PPI_GUT
                 - (CL_PPI + Q_PPI) * C_PPI
                 + Q_PPI * (PPI_PERI / Vp_PPI);
dxdt_PPI_PERI =  Q_PPI * C_PPI - Q_PPI * (PPI_PERI / Vp_PPI);

// ── H2RA PK (2-compartment oral) ─────────────────────────────
dxdt_H2RA_GUT  = -ka_H2RA * H2RA_GUT;
dxdt_H2RA_CENT =  ka_H2RA * F_H2RA * H2RA_GUT - CL_H2RA * C_H2RA;

// ── Amoxicillin PK (2-compartment oral) ──────────────────────
dxdt_AMX_GUT  = -ka_AMX * AMX_GUT;
dxdt_AMX_CENT =  ka_AMX * F_AMX * AMX_GUT - CL_AMX * C_AMX;

// ── Clarithromycin PK (2-compartment oral) ───────────────────
dxdt_CLR_GUT  = -ka_CLR * CLR_GUT;
dxdt_CLR_CENT =  ka_CLR * F_CLR * CLR_GUT - CL_CLR * C_CLR;

// ── NSAID PK (2-compartment oral) ────────────────────────────
dxdt_NSAID_GUT  = -ka_NSAID * NSAID_GUT;
dxdt_NSAID_CENT =  ka_NSAID * F_NSAID * NSAID_GUT - CL_NSAID * C_NSAID;

// ── H. pylori Dynamics ───────────────────────────────────────
// Logistic growth - antibiotic killing - acid killing - immune clearance
double HP = HP_LOAD;
double HP_growth = k_grow_HP * HP * (1.0 - HP / K_HP);
double HP_kill_total = (kill_combo + acid_HP_kill + k_clear_HP) * HP;
dxdt_HP_LOAD = HP_growth - HP_kill_total;

// ── Proton Pump Active Fraction ──────────────────────────────
// New pumps synthesized at k_pumpTurn; PPI binding reduces active fraction
// HP stimulates gastrin → more pumps activated
double HP_stim = HP_acid_stim * (HP / (HP + 1e6));  // normalize
double pump_activation = k_pumpTurn * (1.0 + HP_stim) * (1.0 - PUMP_ACTIVE);
double pump_inhibition_rate = inh_PPI * k_pumpTurn;  // PPI irreversible
dxdt_PUMP_ACTIVE = pump_activation - pump_inhibition_rate * PUMP_ACTIVE;

// ── Intragastric pH ──────────────────────────────────────────
// pH is driven by acid output (pump activity) and buffering
// Low pH = more acid; Total inhibition raises pH toward 7
double acid_rate = (MaxAcidOutput - BasalAcid) * PUMP_ACTIVE * (1.0 - total_acid_inh) + BasalAcid;
double target_pH = pHbasal + (pHmax - pHbasal) * (1.0 - acid_rate / MaxAcidOutput);
// pH relaxes toward target with half-life ~0.5 hr
dxdt_pH_INTRAGASTRIC = 1.5 * (target_pH - pH_INTRAGASTRIC);

// ── Mucosal Defense – Mucus Layer ───────────────────────────
// PGE2 stimulates mucus; acid and H. pylori (VacA) degrade mucus
double mucus_prod = k_mucusProd * PG_LEVEL * PG_NSAID_factor;
double mucus_acid_deg = k_mucusDeg * (1.0 + (1.0 - PG_NSAID_factor)) *
                        exp(-0.5 * (pH - 2.0));  // more degradation at low pH
double mucus_HP_deg = 0.03 * (HP / K_HP);  // HP protease/VacA degrades mucus
dxdt_MUCUS = mucus_prod * (MucusMax - MUCUS) - (mucus_acid_deg + mucus_HP_deg) * MUCUS;
if(MUCUS < 0) MUCUS = 0;

// ── Prostaglandin Level ──────────────────────────────────────
// COX-1 → constitutive PG; NSAID reduces COX-1; inflammation boosts COX-2 PG
double PG_prod = k_PG_base * PG_NSAID_factor + 0.02 * INFLAM;  // COX-2 adds some
double PG_deg  = k_PG_deg * PG_LEVEL;
dxdt_PG_LEVEL = PG_prod - PG_deg;
if(PG_LEVEL < 0) PG_LEVEL = 0;

// ── Inflammation Score ───────────────────────────────────────
// HP-driven (via CagA/IL-8) and NSAID-driven; resolves with HP clearance
double inflam_drive = k_inflam * (HP / K_HP) * InflamMax;
double NSAID_inflam = 0.1 * inh_COX1 * InflamMax;  // COX-1 inhib reduces anti-inflam PGs
double inflam_res   = k_inflam_res * INFLAM * (MUCUS + 0.1);  // mucus helps resolution
dxdt_INFLAM = inflam_drive + NSAID_inflam - inflam_res;
if(INFLAM > InflamMax) INFLAM = InflamMax;
if(INFLAM < 0) INFLAM = 0;

// ── Ulcer Area Dynamics ──────────────────────────────────────
// Damage driven by acid (low pH) × inflammation × mucosal defense failure
// Healing driven by PG, mucus, EGF (simplified)
double mucosal_integrity = MUCUS * PG_LEVEL;  // combined defense
double acid_damage_factor = (pH < 3.0) ? (3.0 - pH) : 0.0;
double ulcer_damage = k_damage * acid_damage_factor * INFLAM / (mucosal_integrity + 0.1);
double EGF_boost = 1.0 + EGF_heal * (1.0 - ULCER_AREA / UlcerMax);
double ulcer_healing = k_heal * EGF_boost * PG_LEVEL * MUCUS * ULCER_AREA;
double ulcer_recur   = k_recur * (HP / K_HP) * (UlcerMax - ULCER_AREA);
dxdt_ULCER_AREA = ulcer_damage * (UlcerMax - ULCER_AREA) + ulcer_recur - ulcer_healing;
if(ULCER_AREA < 0) ULCER_AREA = 0;
if(ULCER_AREA > UlcerMax) ULCER_AREA = UlcerMax;

$TABLE
// ── Derived clinical outputs ─────────────────────────────────
double C_PPI_out     = PPI_CENT / Vc_PPI;
double C_H2RA_out    = H2RA_CENT / Vc_H2RA;
double C_AMX_out     = AMX_CENT / Vc_AMX;
double C_CLR_out     = CLR_CENT / Vc_CLR;
double C_NSAID_out   = NSAID_CENT / Vc_NSAID;

// Intragastric pH
double pH_out        = pH_INTRAGASTRIC;

// Pump inhibition
double pump_inh_pct  = (1.0 - PUMP_ACTIVE) * 100.0;

// HP load as log10
double log10_HP      = (HP_LOAD > 1) ? log10(HP_LOAD) : 0.0;

// HP eradication (defined as < 100 CFU)
double HP_erad_flag  = (HP_LOAD < 100) ? 1.0 : 0.0;

// Mucosal protection score (0-1)
double mucosal_prot  = MUCUS * PG_LEVEL;

// Pain score (0-10 VAS, correlated with ulcer area and inflammation)
double pain_VAS      = 10.0 * (ULCER_AREA / UlcerMax) * (0.5 + 0.5 * INFLAM / 10.0);

// Bleeding risk (increases with ulcer depth proxy)
double bleeding_risk = (ULCER_AREA > 20.0) ?
                       0.01 * (ULCER_AREA - 20.0) / (UlcerMax - 20.0) : 0.0;

// Endoscopy score (Lanza 0-4 scale, simplified)
double endoscopy_score = (ULCER_AREA < 1)  ? 0.0 :
                         (ULCER_AREA < 10) ? 1.0 :
                         (ULCER_AREA < 25) ? 2.0 :
                         (ULCER_AREA < 50) ? 3.0 : 4.0;

// COX inhibition percentages
double COX1_inh_pct  = (1.0 - (1.0 - Imax_COX * C_NSAID_out / (IC50_COX1 + C_NSAID_out))) * 100.0;
double COX2_inh_pct  = (1.0 - (1.0 - Imax_COX * C_NSAID_out / (IC50_COX2 + C_NSAID_out))) * 100.0;

capture C_PPI_out C_H2RA_out C_AMX_out C_CLR_out C_NSAID_out
capture pH_out pump_inh_pct log10_HP HP_erad_flag
capture mucosal_prot pain_VAS bleeding_risk endoscopy_score
capture COX1_inh_pct COX2_inh_pct MUCUS PG_LEVEL INFLAM ULCER_AREA

$INIT
// PK compartments start empty (doses via event objects)
PPI_GUT = 0, PPI_CENT = 0, PPI_PERI = 0
H2RA_GUT = 0, H2RA_CENT = 0
AMX_GUT = 0, AMX_CENT = 0
CLR_GUT = 0, CLR_CENT = 0
NSAID_GUT = 0, NSAID_CENT = 0

// Baseline disease state (HP positive, established ulcer)
HP_LOAD = 1e7          // ~10^7 CFU/mL H. pylori
PUMP_ACTIVE = 0.95     // ~95% pump active at baseline
pH_INTRAGASTRIC = 1.5  // fasting intragastric pH
MUCUS = 0.6            // partially degraded mucus
PG_LEVEL = 0.7         // reduced prostaglandins (HP effect)
INFLAM = 4.0           // moderate inflammation (HP driven)
ULCER_AREA = 25.0      // active peptic ulcer (25 mm²)
'

# Compile the model
mod <- mcode("PUD_QSP", code)

# ============================================================
# Helper functions
# ============================================================

make_events <- function(scenario, dur = 168) {
  # Returns an event data frame
  evts <- list()

  if (scenario == "ppi_bid") {
    # PPI 20mg BID (omeprazole-equivalent) for 4 weeks
    evts[["ppi"]] <- ev(amt = 20, cmt = "PPI_GUT", ii = 12, addl = 27, time = 0)

  } else if (scenario == "h2ra_bid") {
    # Famotidine 20mg BID
    evts[["h2ra"]] <- ev(amt = 20, cmt = "H2RA_GUT", ii = 12, addl = 27, time = 0)

  } else if (scenario == "triple_therapy") {
    # Standard triple therapy (PPI + AMX + CLR, 14 days)
    evts[["ppi"]]  <- ev(amt = 20, cmt = "PPI_GUT",  ii = 12, addl = 27, time = 0)
    evts[["amx"]]  <- ev(amt = 1000, cmt = "AMX_GUT", ii = 12, addl = 27, time = 0)
    evts[["clr"]]  <- ev(amt = 500, cmt = "CLR_GUT",  ii = 12, addl = 27, time = 0)

  } else if (scenario == "nsaid_only") {
    # Ibuprofen 400mg TID (NSAID gastropathy model)
    evts[["nsaid"]] <- ev(amt = 400, cmt = "NSAID_GUT", ii = 8, addl = 41, time = 0)

  } else if (scenario == "nsaid_ppi") {
    # Ibuprofen + PPI co-therapy (gastroprotection)
    evts[["nsaid"]] <- ev(amt = 400, cmt = "NSAID_GUT", ii = 8, addl = 41, time = 0)
    evts[["ppi"]]   <- ev(amt = 20,  cmt = "PPI_GUT",   ii = 12, addl = 27, time = 0)
  }

  if (length(evts) == 0) return(ev(amt = 0, cmt = 1, time = 999))
  do.call(c, evts)
}

# ============================================================
# Scenario 1: No treatment (natural history)
# ============================================================
sim_no_tx <- mod %>%
  param(HP0 = 1e7) %>%
  mrgsim(end = 672, delta = 1) %>%  # 28 days
  as.data.frame()

# ============================================================
# Scenario 2: PPI monotherapy (omeprazole 20mg BID)
# ============================================================
sim_ppi <- mod %>%
  ev(make_events("ppi_bid")) %>%
  mrgsim(end = 672, delta = 1) %>%
  as.data.frame()

# ============================================================
# Scenario 3: H2RA monotherapy (famotidine 20mg BID)
# ============================================================
sim_h2ra <- mod %>%
  ev(make_events("h2ra_bid")) %>%
  mrgsim(end = 672, delta = 1) %>%
  as.data.frame()

# ============================================================
# Scenario 4: H. pylori Triple Therapy (standard 14-day)
# ============================================================
sim_triple <- mod %>%
  ev(make_events("triple_therapy")) %>%
  mrgsim(end = 672, delta = 1) %>%
  as.data.frame()

# ============================================================
# Scenario 5: NSAID-induced gastropathy (no protection)
# ============================================================
sim_nsaid <- mod %>%
  init(HP_LOAD = 100, ULCER_AREA = 0, INFLAM = 1.0) %>%  # HP-negative, healthy start
  ev(make_events("nsaid_only")) %>%
  mrgsim(end = 336, delta = 1) %>%  # 14 days
  as.data.frame()

# ============================================================
# Scenario 5b: NSAID + PPI co-therapy (gastroprotection)
# ============================================================
sim_nsaid_ppi <- mod %>%
  init(HP_LOAD = 100, ULCER_AREA = 0, INFLAM = 1.0) %>%
  ev(make_events("nsaid_ppi")) %>%
  mrgsim(end = 336, delta = 1) %>%
  as.data.frame()

# ============================================================
# Population variability simulation (n=100, PPI BID)
# ============================================================
set.seed(42)
n_pop <- 100

# CYP2C19 genotype variability affects PPI PK
pop_param <- data.frame(
  ID = 1:n_pop,
  # Poor metabolizers (PM) ~3%, Intermediate (IM) ~30%, Normal (NM) ~67%
  CYP_type = sample(c("PM", "IM", "NM"), n_pop,
                    prob = c(0.03, 0.30, 0.67), replace = TRUE)
)
pop_param$CL_PPI <- ifelse(pop_param$CYP_type == "PM", 8.0,
                    ifelse(pop_param$CYP_type == "IM", 18.0, 28.0))
pop_param$F_PPI  <- ifelse(pop_param$CYP_type == "PM", 0.90,
                    ifelse(pop_param$CYP_type == "IM", 0.75, 0.65))

sim_pop <- mod %>%
  idata_set(pop_param) %>%
  ev(make_events("ppi_bid")) %>%
  mrgsim(end = 336, delta = 2) %>%
  as.data.frame()

# ============================================================
# Summary Statistics
# ============================================================
library(dplyr)

summary_at_end <- function(df, scenario_name) {
  df %>%
    filter(time == max(time)) %>%
    summarise(
      scenario     = scenario_name,
      pH_mean      = mean(pH_out),
      pump_inh_pct = mean(pump_inh_pct),
      log10_HP     = mean(log10_HP),
      HP_erad      = mean(HP_erad_flag),
      ulcer_area   = mean(ULCER_AREA),
      pain_VAS     = mean(pain_VAS),
      endoscopy    = mean(endoscopy_score)
    )
}

results_summary <- bind_rows(
  summary_at_end(sim_no_tx,    "No Treatment"),
  summary_at_end(sim_ppi,      "PPI BID"),
  summary_at_end(sim_h2ra,     "H2RA BID"),
  summary_at_end(sim_triple,   "Triple Therapy"),
  summary_at_end(sim_nsaid,    "NSAID Only"),
  summary_at_end(sim_nsaid_ppi,"NSAID + PPI")
)

print(results_summary)

# ============================================================
# Key Visualization Plots
# ============================================================
library(ggplot2)
library(patchwork)

# Combined scenario data
all_sims <- bind_rows(
  sim_no_tx    %>% mutate(scenario = "No Treatment", time_day = time/24),
  sim_ppi      %>% mutate(scenario = "PPI 20mg BID", time_day = time/24),
  sim_h2ra     %>% mutate(scenario = "H2RA 20mg BID", time_day = time/24),
  sim_triple   %>% mutate(scenario = "Triple Therapy", time_day = time/24)
)

# 1. Intragastric pH over time
p1 <- ggplot(all_sims, aes(time_day, pH_out, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 4.0, linetype = "dashed", color = "gray50") +
  labs(title = "Intragastric pH Over Time",
       x = "Time (days)", y = "Intragastric pH",
       subtitle = "Target pH > 4 for ulcer healing") +
  theme_bw() + theme(legend.position = "bottom")

# 2. H. pylori Load over time
p2 <- ggplot(all_sims, aes(time_day, log10_HP, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = log10(100), linetype = "dashed", color = "red") +
  labs(title = "H. pylori Load",
       x = "Time (days)", y = "log10(CFU/mL)",
       subtitle = "Red dashed = eradication threshold") +
  theme_bw() + theme(legend.position = "bottom")

# 3. Ulcer Area over time
p3 <- ggplot(all_sims, aes(time_day, ULCER_AREA, color = scenario)) +
  geom_line(linewidth = 1) +
  labs(title = "Ulcer Area Dynamics",
       x = "Time (days)", y = "Ulcer Area (mm²)") +
  theme_bw() + theme(legend.position = "bottom")

# 4. Pain VAS over time
p4 <- ggplot(all_sims, aes(time_day, pain_VAS, color = scenario)) +
  geom_line(linewidth = 1) +
  labs(title = "Epigastric Pain (VAS 0-10)",
       x = "Time (days)", y = "Pain Score") +
  scale_y_continuous(limits = c(0, 10)) +
  theme_bw() + theme(legend.position = "bottom")

# NSAID comparison
nsaid_sims <- bind_rows(
  sim_nsaid     %>% mutate(scenario = "NSAID Only",   time_day = time/24),
  sim_nsaid_ppi %>% mutate(scenario = "NSAID + PPI", time_day = time/24)
)

# 5. NSAID gastropathy: ulcer area
p5 <- ggplot(nsaid_sims, aes(time_day, ULCER_AREA, color = scenario)) +
  geom_line(linewidth = 1.2) +
  labs(title = "NSAID Gastropathy Prevention",
       x = "Time (days)", y = "Ulcer Area (mm²)") +
  theme_bw() + theme(legend.position = "bottom")

# 6. Population variability (PPI PK by CYP2C19)
p6 <- sim_pop %>%
  mutate(time_day = time/24) %>%
  left_join(pop_param, by = "ID") %>%
  filter(time_day %% 1 == 0) %>%
  ggplot(aes(time_day, pH_out, group = ID, color = CYP_type)) +
  geom_line(alpha = 0.4) +
  stat_summary(aes(group = CYP_type), fun = mean, geom = "line", linewidth = 1.5) +
  labs(title = "PPI pH Response by CYP2C19 Phenotype",
       x = "Time (days)", y = "Intragastric pH",
       color = "CYP2C19") +
  theme_bw() + theme(legend.position = "right")

combined_plot <- (p1 + p2) / (p3 + p4) / (p5 + p6)
print(combined_plot)

cat("\n=== Peptic Ulcer Disease QSP Model Summary ===\n")
cat("Compartments:", 20, "ODEs\n")
cat("Scenarios simulated:", 6, "\n")
cat("H. pylori eradication rate (triple therapy):",
    round(mean(sim_triple$HP_erad_flag[sim_triple$time == max(sim_triple$time)]) * 100, 1), "%\n")
cat("Ulcer healing at 4 weeks (PPI BID):",
    round((1 - mean(sim_ppi$ULCER_AREA[sim_ppi$time == max(sim_ppi$time)]) / 25) * 100, 1), "% reduction\n")
