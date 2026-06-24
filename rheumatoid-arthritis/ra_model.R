## ============================================================
## Rheumatoid Arthritis — QSP mrgsolve Model
## ============================================================
## Model: Tocilizumab (IL-6Rα mAb) +/- Methotrexate in RA
##
## Architecture:
##   PK  : Tocilizumab 2-compartment + TMDD (sIL6R/mIL6R)
##         Methotrexate 1-compartment oral
##   PD  : TNF-α dynamics (MTX-sensitive)
##         IL-6 dynamics (tocilizumab-sensitive via TMDD)
##         CRP production (STAT3/IL-6 driven)
##         RANKL dynamics (bone erosion driver)
##         DAS28-CRP, ACR response rates
##
## Key references:
##   Gibiansky et al. J Clin Pharmacol 2012 (TCZ PK)
##   Levi et al. JPP 2018 (TCZ CRP PK/PD)
##   Simeoni et al. JPP 2019 (cytokine QSP)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ── mrgsolve model code ────────────────────────────────────────
code <- '
$PROB
Rheumatoid Arthritis QSP Model
Tocilizumab + Methotrexate

$PARAM @annotated
// ── Tocilizumab PK (2-compartment, SC/IV) ──
CL_TCZ   : 0.224  : L/day  | TCZ central clearance
V1_TCZ   : 3.72   : L      | TCZ central volume
Q_TCZ    : 0.697  : L/day  | TCZ intercompartmental clearance
V2_TCZ   : 2.91   : L      | TCZ peripheral volume
F_SC_TCZ : 0.80   : .      | SC bioavailability
ka_SC_TCZ: 0.29   : 1/day  | SC absorption rate constant

// ── TMDD — sIL6R / mIL6R ──
kon_TCZ  : 0.0272  : 1/(nM*day) | TCZ-IL6R association rate
koff_TCZ : 0.0096  : 1/day      | TCZ-IL6R dissociation rate
ksyn_R   : 0.132   : nM/day     | sIL6R synthesis rate
kdeg_R   : 0.066   : 1/day      | free sIL6R degradation rate
kdeg_RC  : 0.016   : 1/day      | TCZ-IL6R complex degradation rate

// ── Methotrexate PK (1-compartment, oral) ──
F_MTX    : 0.70   : .      | MTX oral bioavailability
ka_MTX   : 1.50   : 1/hr   | MTX absorption rate (1/hr)
CL_MTX   : 0.09   : L/hr   | MTX total clearance
V_MTX    : 0.80   : L/kg   | MTX volume of distribution (scaled by WT)
WT       : 70     : kg     | body weight

// ── TNF-α PD ──
ksyn_TNF   : 1.20   : ng/mL/day | TNF-α zero-order synthesis (RA baseline)
kout_TNF   : 0.50   : 1/day     | TNF-α first-order elimination
TNF0       : 8.0    : ng/mL     | healthy baseline TNF-α
RA_TNF_fold: 3.0    : .         | fold increase in RA vs healthy
Emax_MTX   : 0.50   : .         | MTX maximum TNF reduction (Emax)
EC50_MTX   : 0.05   : mg/L      | MTX plasma EC50 for TNF inhibition

// ── IL-6 PD ──
ksyn_IL6   : 6.0    : pg/mL/day | IL-6 zero-order synthesis in RA
kout_IL6   : 0.30   : 1/day     | IL-6 first-order elimination
IL60       : 5.0    : pg/mL     | healthy baseline IL-6
kfb_TNF_IL6: 0.08   : .         | TNF→IL-6 positive feedback coefficient

// ── CRP PD (STAT3/IL-6 driven) ──
ksyn_CRP   : 8.0    : mg/L/day  | CRP maximum synthesis rate
kout_CRP   : 0.693  : 1/day     | CRP elimination (t½~1 day)
EC50_IL6_CRP: 8.0   : pg/mL    | IL-6 EC50 for CRP synthesis
CRP0       : 2.0    : mg/L      | healthy baseline CRP
CRP_RA0    : 20.0   : mg/L      | RA baseline CRP (pre-treatment)

// ── RANKL/OPG dynamics (bone erosion risk) ──
ksyn_RANKL  : 0.50  : pmol/L/day | RANKL synthesis (baseline RA)
kout_RANKL  : 0.20  : 1/day      | RANKL clearance
RANKL0      : 2.0   : pmol/L     | healthy baseline RANKL
TNF_RANKL_slope: 0.30: .         | TNF→RANKL production slope

// ── DAS28 parameters ──
DAS28_base  : 5.8   : .          | pre-treatment DAS28-CRP
SJC28_base  : 14    : .          | pre-treatment swollen joint count
TJC28_base  : 18    : .          | pre-treatment tender joint count
PatGlob_base: 65    : mm         | pre-treatment patient global (VAS 0-100)

// ── Joint adaptation kinetics ──
k_SJC  : 0.15   : 1/day | SJC28 adaptation rate to inflammation
k_TJC  : 0.10   : 1/day | TJC28 adaptation rate (pain-driven)

$CMT @annotated
// TCZ PK
DEPOT_TCZ : SC depot for tocilizumab
C1_TCZ    : Central plasma compartment (mg)
C2_TCZ    : Peripheral tissue compartment (mg)
R_FREE    : Free soluble IL-6 receptor (nM)
RC        : TCZ-IL6R complex (nM)

// MTX PK
GI_MTX    : GI absorption compartment (mg)
C1_MTX    : Plasma MTX (mg)

// Biomarker/PD states
TNFa      : TNF-alpha (ng/mL)
IL6       : Interleukin-6 (pg/mL)
CRP       : C-reactive protein (mg/L)
RANKL_pd  : RANKL (pmol/L)

// Clinical scores (ODE-driven for smoothness)
SJC28_ode  : Swollen joint count (ODE)
TJC28_ode  : Tender joint count (ODE)

$MAIN
// Steady-state initial conditions
R_FREE_0 = ksyn_R / kdeg_R;         // free sIL6R at SS (~2 nM)

// Disease-state baseline for PD states
TNFa_0 = ksyn_TNF / kout_TNF;       // RA TNF-α baseline
IL6_0  = ksyn_IL6 / kout_IL6;       // RA IL-6 baseline
CRP_0  = CRP_RA0;                   // RA CRP baseline

// Effective MTX volume
double Vd_MTX = V_MTX * WT;         // L (body weight scaled)

$ODE
// ── TCZ PK ──
double C_TCZ = C1_TCZ / V1_TCZ;     // TCZ plasma conc [mg/L] = [nM * MW/1000]
double C_TCZ_nM = C_TCZ / 0.148;    // convert mg/L to nM (MW 148 kDa)

dxdt_DEPOT_TCZ = -ka_SC_TCZ * DEPOT_TCZ;

dxdt_C1_TCZ = ka_SC_TCZ * F_SC_TCZ * DEPOT_TCZ
              - (CL_TCZ / V1_TCZ) * C1_TCZ
              - (Q_TCZ  / V1_TCZ) * C1_TCZ
              + (Q_TCZ  / V2_TCZ) * C2_TCZ
              - kon_TCZ * C_TCZ_nM * R_FREE * V1_TCZ   // drug-target binding
              + koff_TCZ * RC * V1_TCZ;

dxdt_C2_TCZ = (Q_TCZ / V1_TCZ) * C1_TCZ - (Q_TCZ / V2_TCZ) * C2_TCZ;

// ── TMDD: IL-6 receptor ──
dxdt_R_FREE = ksyn_R
              - kdeg_R * R_FREE
              - kon_TCZ * C_TCZ_nM * R_FREE
              + koff_TCZ * RC;

dxdt_RC = kon_TCZ * C_TCZ_nM * R_FREE
          - koff_TCZ * RC
          - kdeg_RC * RC;

// ── MTX PK ──
double CL_MTX_total = CL_MTX;       // L/hr
dxdt_GI_MTX = -ka_MTX * GI_MTX;
dxdt_C1_MTX = ka_MTX * F_MTX * GI_MTX
              - (CL_MTX_total / Vd_MTX) * C1_MTX;

double C_MTX_mgL = C1_MTX / Vd_MTX; // mg/L

// MTX effect on TNF (sigmoid Emax)
double MTX_Eff = Emax_MTX * C_MTX_mgL / (EC50_MTX + C_MTX_mgL);

// ── TNF-α dynamics ──
dxdt_TNFa = ksyn_TNF * (1.0 - MTX_Eff)
            - kout_TNF * TNFa;

// ── IL-6 dynamics ──
// TNF-driven feedback; receptor blockade reduces effective IL-6 action
// (IL6 production unaffected — IL6R blockade changes downstream, not IL6 itself)
double TNF_norm = TNFa / (RA_TNF_fold * TNF0);
dxdt_IL6 = ksyn_IL6 * (1.0 + kfb_TNF_IL6 * (TNF_norm - 1.0))
           - kout_IL6 * IL6;

// ── CRP dynamics ──
// CRP synthesis driven by free IL-6 signaling
// IL-6R blockade (high RC) reduces effective signaling
double R_total = R_FREE + RC + 1e-6;
double fR_blocked = RC / R_total;           // fraction IL-6R occupied by TCZ
double IL6_effective = IL6 * (1.0 - fR_blocked); // reduced signal to liver

double CRP_stim = ksyn_CRP * IL6_effective / (EC50_IL6_CRP + IL6_effective);
dxdt_CRP = CRP_stim - kout_CRP * CRP;

// ── RANKL dynamics ──
// TNF and IL-6 both upregulate RANKL on FLS and osteoblasts
double RANKL_drive = ksyn_RANKL * (1.0 + TNF_RANKL_slope * (TNF_norm - 1.0));
dxdt_RANKL_pd = RANKL_drive - kout_RANKL * RANKL_pd;

// ── Clinical score adaptation ──
// Inflammation signal → synovitis → joint counts
double inflam_idx = 0.5 * (TNFa / (RA_TNF_fold * TNF0))
                  + 0.5 * (IL6  / (ksyn_IL6 / kout_IL6));
inflam_idx = inflam_idx < 0.01 ? 0.01 : inflam_idx;

double SJC28_ss = SJC28_base * inflam_idx;
double TJC28_ss = TJC28_base * (0.6 * inflam_idx + 0.4 * (CRP / CRP_RA0));

dxdt_SJC28_ode = k_SJC * (SJC28_ss - SJC28_ode);
dxdt_TJC28_ode = k_TJC * (TJC28_ss - TJC28_ode);

$TABLE
// Concentrations
capture C_TCZ_mgL   = C1_TCZ / V1_TCZ;
capture C_TCZ_nM    = (C1_TCZ / V1_TCZ) / 0.148;
capture C_MTX_pl    = C1_MTX / (V_MTX * WT);
capture fR_blocked  = RC / (R_FREE + RC + 1e-6);

// Biomarkers
capture TNFa_out    = TNFa;
capture IL6_out     = IL6;
capture CRP_out     = CRP;
capture RANKL_out   = RANKL_pd;
capture R_FREE_out  = R_FREE;
capture RC_out      = RC;

// DAS28-CRP formula:
// DAS28-CRP = 0.56*sqrt(TJC28) + 0.28*sqrt(SJC28) + 0.36*ln(CRP+1) + 0.014*PatGlob + 0.96
double PatGlob_dyn  = PatGlob_base * (CRP / CRP_RA0) * 0.6
                    + PatGlob_base * (TJC28_ode / TJC28_base) * 0.4;
PatGlob_dyn = PatGlob_dyn < 5.0 ? 5.0 : PatGlob_dyn;

capture SJC28_out   = SJC28_ode;
capture TJC28_out   = TJC28_ode;
capture PatGlob_out = PatGlob_dyn;

double DAS28_CRP_val = 0.56 * sqrt(TJC28_ode)
                     + 0.28 * sqrt(SJC28_ode)
                     + 0.36 * log(CRP + 1.0)
                     + 0.014 * PatGlob_dyn
                     + 0.96;
capture DAS28_CRP_out = DAS28_CRP_val;

// ACR response (% improvement vs baseline DAS28)
double pct_impr = (DAS28_base - DAS28_CRP_val) / DAS28_base * 100.0;
pct_impr = pct_impr < 0.0 ? 0.0 : pct_impr;
capture pct_improvement = pct_impr;
capture ACR20 = (pct_impr >= 20.0) ? 1.0 : 0.0;
capture ACR50 = (pct_impr >= 50.0) ? 1.0 : 0.0;
capture ACR70 = (pct_impr >= 70.0) ? 1.0 : 0.0;

// EULAR response categories (based on absolute DAS28 and change)
double dDAS28 = DAS28_base - DAS28_CRP_val;
capture EULAR_good     = (DAS28_CRP_val < 3.2 && dDAS28 > 1.2) ? 1.0 : 0.0;
capture EULAR_moderate = ((DAS28_CRP_val < 5.1 && dDAS28 > 1.2) ||
                          (DAS28_CRP_val < 3.2 && dDAS28 >= 0.6)) ? 1.0 : 0.0;
capture EULAR_remission = (DAS28_CRP_val < 2.6) ? 1.0 : 0.0;

// HAQ-DI approximation: correlates with disease activity
capture HAQ_DI = 0.2 * (DAS28_CRP_val / DAS28_base) * 2.5;  // scaled 0-2.5

$CAPTURE C_TCZ_mgL C_TCZ_nM C_MTX_pl fR_blocked
         TNFa_out IL6_out CRP_out RANKL_out R_FREE_out RC_out
         SJC28_out TJC28_out PatGlob_out DAS28_CRP_out pct_improvement
         ACR20 ACR50 ACR70 EULAR_good EULAR_moderate EULAR_remission HAQ_DI
'

## ── Compile model ─────────────────────────────────────────────
mod <- mcode("ra_qsp", code)
mod <- update(mod, end = 168, delta = 1, outvars = "all")  # 24 weeks (168 days)


## ── Helper: dosing event builder ─────────────────────────────
make_events <- function(
    tcz_dose_mg   = 162,   # mg (SC; e.g. 162 mg q2w)
    tcz_interval  = 14,    # days between TCZ doses
    tcz_nDose     = 12,    # number of TCZ doses
    tcz_route     = "SC",  # "SC" or "IV"
    mtx_dose_mg   = 15,    # mg (oral weekly)
    mtx_weekly    = TRUE,  # dose MTX weekly
    mtx_nDose     = 24,    # doses over 24 weeks
    tcz_only      = FALSE, # TRUE = no MTX
    sim_dur       = 168    # days
) {

  ev_list <- list()

  # TCZ dosing
  tcz_times <- seq(0, by = tcz_interval, length.out = tcz_nDose)
  tcz_times <- tcz_times[tcz_times <= sim_dur]
  if (tcz_route == "SC") {
    ev_list[[1]] <- ev(cmt = "DEPOT_TCZ", amt = tcz_dose_mg,
                       time = tcz_times, rate = 0)
  } else {
    ev_list[[1]] <- ev(cmt = "C1_TCZ", amt = tcz_dose_mg,
                       time = tcz_times, rate = tcz_dose_mg / (1/24))  # 1-hr infusion
  }

  # MTX dosing (weekly, oral)
  if (!tcz_only) {
    mtx_times <- seq(0, by = 7, length.out = mtx_nDose)
    mtx_times <- mtx_times[mtx_times <= sim_dur]
    ev_list[[2]] <- ev(cmt = "GI_MTX", amt = mtx_dose_mg,
                       time = mtx_times, rate = 0)
  }

  # Combine
  do.call(c, ev_list[!sapply(ev_list, is.null)])
}


## ── Initial conditions (RA disease state pre-treatment) ──────
init_ra <- function(m) {
  params <- param(m)
  ksyn_R  <- params$ksyn_R
  kdeg_R  <- params$kdeg_R
  ksyn_TNF <- params$ksyn_TNF
  kout_TNF <- params$kout_TNF
  ksyn_IL6 <- params$ksyn_IL6
  kout_IL6 <- params$kout_IL6
  CRP_RA0  <- params$CRP_RA0
  ksyn_RANKL <- params$ksyn_RANKL
  kout_RANKL <- params$kout_RANKL

  m <- init(m,
    R_FREE    = ksyn_R / kdeg_R,          # ~2 nM sIL6R
    TNFa      = ksyn_TNF / kout_TNF,      # RA baseline TNF
    IL6       = ksyn_IL6 / kout_IL6,      # RA baseline IL-6
    CRP       = CRP_RA0,                  # RA baseline CRP
    RANKL_pd  = ksyn_RANKL / kout_RANKL,  # RA baseline RANKL
    SJC28_ode = params$SJC28_base,
    TJC28_ode = params$TJC28_base,
    RC        = 0,
    C1_TCZ    = 0, C2_TCZ = 0, DEPOT_TCZ = 0,
    GI_MTX    = 0, C1_MTX = 0
  )
  m
}


## ── Scenario definitions ────────────────────────────────────
scenarios <- list(
  "TCZ_162mg_q2w_SC_plus_MTX" = list(
    tcz_dose_mg = 162, tcz_interval = 14, tcz_route = "SC",
    tcz_only = FALSE, mtx_dose_mg = 15
  ),
  "TCZ_162mg_qw_SC_monotherapy" = list(
    tcz_dose_mg = 162, tcz_interval = 7, tcz_route = "SC",
    tcz_only = TRUE
  ),
  "TCZ_8mgkg_q4w_IV_plus_MTX" = list(
    tcz_dose_mg = 560, tcz_interval = 28, tcz_route = "IV",   # ~8 mg/kg x 70 kg
    tcz_only = FALSE, mtx_dose_mg = 15
  ),
  "MTX_15mg_weekly_monotherapy" = list(
    tcz_dose_mg = 0, tcz_only = TRUE,   # no TCZ
    mtx_dose_mg = 15
  )
)


## ── Run simulations ─────────────────────────────────────────
mod_ra <- init_ra(mod)

run_scenario <- function(scen_params, n_subj = 1) {
  ev_tcz <- make_events(
    tcz_dose_mg  = scen_params$tcz_dose_mg,
    tcz_interval = ifelse(!is.null(scen_params$tcz_interval), scen_params$tcz_interval, 14),
    tcz_route    = scen_params$tcz_route,
    tcz_only     = scen_params$tcz_only,
    mtx_dose_mg  = ifelse(!is.null(scen_params$mtx_dose_mg), scen_params$mtx_dose_mg, 0)
  )
  mrgsim(mod_ra, events = ev_tcz, end = 168, delta = 0.5, carry_out = "evid")
}

## Run all scenarios
results <- lapply(scenarios, run_scenario)
names(results) <- names(scenarios)


## ── Plots ───────────────────────────────────────────────────
theme_set(theme_bw(base_size = 12))

# Helper to extract tidy data from mrgsim output
tidy_sim <- function(sim_out, scenario_name) {
  as.data.frame(sim_out) %>%
    mutate(Scenario = scenario_name) %>%
    filter(evid == 0)
}

# Combine all scenarios
all_res <- bind_rows(
  mapply(tidy_sim, results, names(results), SIMPLIFY = FALSE)
)

# ── Plot 1: Tocilizumab PK ──
p1 <- all_res %>%
  filter(grepl("TCZ", Scenario)) %>%
  ggplot(aes(time, C_TCZ_mgL, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous("Time (days)", breaks = seq(0, 168, 14)) +
  scale_y_continuous("Tocilizumab Plasma Conc. (mg/L)") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Tocilizumab Plasma Concentration",
       subtitle = "2-Compartment + TMDD Model (SC/IV)") +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

# ── Plot 2: IL-6 Receptor Blockade ──
p2 <- all_res %>%
  filter(grepl("TCZ", Scenario)) %>%
  ggplot(aes(time, fR_blocked * 100, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 80, linetype = 2, color = "gray40") +
  scale_x_continuous("Time (days)", breaks = seq(0, 168, 14)) +
  scale_y_continuous("IL-6 Receptor Blockade (%)", limits = c(0, 100)) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "IL-6 Receptor Occupancy by Tocilizumab",
       subtitle = "TMDD: sIL6Rα fraction occupied")

# ── Plot 3: CRP Kinetics ──
p3 <- all_res %>%
  ggplot(aes(time, CRP_out, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 5, linetype = 2, color = "gray40") +
  annotate("text", x = 168, y = 5.5, label = "Normal CRP (<5 mg/L)", size = 3, hjust = 1) +
  scale_x_continuous("Time (days)", breaks = seq(0, 168, 14)) +
  scale_y_continuous("CRP (mg/L)") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "C-Reactive Protein Kinetics",
       subtitle = "IL-6/STAT3-driven CRP synthesis")

# ── Plot 4: TNF-α and IL-6 ──
p4a <- all_res %>%
  ggplot(aes(time, TNFa_out, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("TNF-α (ng/mL)") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "TNF-α", subtitle = "MTX-sensitive pathway")

p4b <- all_res %>%
  ggplot(aes(time, IL6_out, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("IL-6 (pg/mL)") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "IL-6", subtitle = "TNF-driven feedback")

# ── Plot 5: DAS28-CRP over time ──
p5 <- all_res %>%
  ggplot(aes(time, DAS28_CRP_out, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 2.6, linetype = 2, color = "green4",  linewidth = 0.8) +
  geom_hline(yintercept = 3.2, linetype = 2, color = "orange",  linewidth = 0.8) +
  geom_hline(yintercept = 5.1, linetype = 2, color = "red3",    linewidth = 0.8) +
  annotate("text", x=5, y=2.4, label="Remission (<2.6)", size=3, color="green4", hjust=0) +
  annotate("text", x=5, y=3.0, label="Low DA (<3.2)",   size=3, color="orange",  hjust=0) +
  annotate("text", x=5, y=4.9, label="Moderate DA (<5.1)", size=3, color="red3", hjust=0) +
  scale_x_continuous("Time (days)", breaks = seq(0, 168, 28)) +
  scale_y_continuous("DAS28-CRP Score", limits = c(1, 7)) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "DAS28-CRP Disease Activity Score",
       subtitle = "Clinical response categories shown")

# ── Plot 6: ACR response probability over time ──
acr_df <- all_res %>%
  select(time, Scenario, ACR20, ACR50, ACR70) %>%
  pivot_longer(ACR20:ACR70, names_to = "Threshold", values_to = "Response")

p6 <- acr_df %>%
  ggplot(aes(time, Response, color = Scenario, linetype = Threshold)) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous("Time (days)", breaks = seq(0, 168, 28)) +
  scale_y_continuous("ACR Response (1 = achieved)", limits = c(-0.05, 1.1)) +
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(values = c("solid","dashed","dotted")) +
  labs(title = "ACR20/50/70 Response",
       subtitle = "1 = Response achieved, 0 = Not achieved")

# ── Plot 7: RANKL (bone erosion driver) ──
p7 <- all_res %>%
  ggplot(aes(time, RANKL_out, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous("Time (days)", breaks = seq(0, 168, 28)) +
  scale_y_continuous("RANKL (pmol/L)") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "RANKL Dynamics",
       subtitle = "TNF-driven bone erosion signal")

# ── Summary table at Week 24 ────────────────────────────────
summary_w24 <- all_res %>%
  filter(abs(time - 168) < 0.6) %>%
  group_by(Scenario) %>%
  slice(1) %>%
  select(Scenario, CRP_out, DAS28_CRP_out, ACR20, ACR50, ACR70,
         EULAR_remission, HAQ_DI, RANKL_out) %>%
  rename(
    `CRP (mg/L)` = CRP_out,
    `DAS28-CRP`  = DAS28_CRP_out,
    `EULAR Rem.` = EULAR_remission
  )

cat("\n=== Week 24 Summary ===\n")
print(as.data.frame(summary_w24), digits = 3)

# Return plots and summary for use in Shiny
list(
  plots = list(pk = p1, receptor = p2, crp = p3,
               tnf = p4a, il6 = p4b, das28 = p5,
               acr = p6, rankl = p7),
  summary = summary_w24,
  data    = all_res,
  model   = mod_ra
)
