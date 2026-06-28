## ============================================================
## Acute Kidney Injury (AKI) — QSP mrgsolve ODE Model
## ============================================================
## Disease: Acute Kidney Injury (급성 신손상)
## Subtypes: Ischemia-Reperfusion (IRI), Nephrotoxic (NTX),
##           Sepsis-Associated (SA-AKI)
## Framework: mrgsolve + R
## Author: Claude Code Routine (CCR) | Date: 2026-06-28
##
## Key References:
##  - Rabb et al. JASN 2016 (IRI mechanisms)
##  - Mehta et al. Lancet 2015 (KDIGO staging)
##  - Bellomo et al. Crit Care 2012 (CRRT outcomes)
##  - Meersch et al. Lancet 2017 (TIMP-2·IGFBP7)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ============================================================
## MODEL CODE BLOCK
## ============================================================

code <- '
$PROB
AKI QSP Model: Ischemia-Reperfusion / Nephrotoxic / Sepsis
20-compartment ODE system
Drug PK: Furosemide (2-cmpt), Norepinephrine (1-cmpt), NAC (1-cmpt)
Disease PD: GFR, Tubular cell viability, ROS, Inflammation (IL6/TNFa),
            Oxidative stress, Biomarkers (NGAL, KIM1, Cystatin C, sCr),
            Repair pathway, Fibrosis (TGF-b / AKI-to-CKD)

$PARAM
// ─── AKI Trigger Parameters ───
AKI_type   = 1      // 1=IRI, 2=NTX, 3=Sepsis
isch_sev   = 0.6    // Ischemia severity [0-1]
tox_dose   = 0      // Nephrotoxin dose (mg/kg cisplatin equiv)
lps_amt    = 0      // LPS burden (ng/mL equiv)

// ─── Baseline Physiology ───
GFR0       = 100    // Baseline GFR (mL/min/1.73m2)
TCV0       = 1.0    // Tubular Cell Viability (normalized)
RBF0       = 1.0    // Renal Blood Flow (normalized, 1=normal)
SCR0       = 0.9    // Baseline serum creatinine (mg/dL)

// ─── Tubular Injury Kinetics ───
k_inj      = 0.15   // Rate of tubular injury (h-1)
k_apop     = 0.08   // Apoptosis rate constant (h-1)
k_nec      = 0.04   // Necrosis rate constant (h-1)
k_ferr     = 0.03   // Ferroptosis rate constant (h-1)
EC50_ros   = 0.5    // ROS EC50 for injury (normalized)
EC50_atp   = 0.4    // ATP depletion EC50 for injury
n_hill     = 2      // Hill coefficient

// ─── GFR Dynamics ───
k_gfr_rec  = 0.005  // GFR recovery rate (h-1) after injury resolves
k_gfr_inj  = 0.10   // GFR decline rate from tubular injury (h-1)
k_backleak = 0.02   // Backleak contribution to effective GFR loss
k_obstruct = 0.015  // Tubular obstruction contribution

// ─── Inflammatory Mediators ───
k_il6_prod = 0.3    // IL-6 production rate (pg/mL/h)
k_il6_deg  = 0.1    // IL-6 degradation (h-1)
k_tnfa_prod= 0.2    // TNF-α production rate (pg/mL/h)
k_tnfa_deg = 0.15   // TNF-α degradation (h-1)
IL6_base   = 5      // Baseline IL-6 (pg/mL)
TNFa_base  = 3      // Baseline TNF-α (pg/mL)
k_nfkb     = 0.5    // NF-κB activation by injury/ROS

// ─── Oxidative Stress ───
k_ros_prod = 0.2    // ROS production (h-1)
k_ros_deg  = 0.15   // ROS degradation/antioxidant (h-1)
k_gsh_prod = 0.05   // GSH replenishment rate (h-1)
k_gsh_dep  = 0.1    // GSH depletion by ROS (h-1)
GSH0       = 1.0    // Baseline GSH (normalized)
ROS0       = 0.1    // Baseline ROS (normalized)

// ─── ATP Dynamics ───
k_atp_prod = 0.2    // ATP production via OXPHOS (h-1)
k_atp_dep  = 0.15   // ATP consumption (h-1)
ATP0       = 1.0    // Baseline ATP (normalized, 1=normal)

// ─── Biomarker Kinetics ───
k_ngal_prod= 0.5    // NGAL production from injured cells (ng/mL/h)
k_ngal_deg = 0.12   // NGAL elimination (h-1)
k_kim1_prod= 0.3    // KIM-1 shedding rate (ng/mL/h)
k_kim1_deg = 0.08   // KIM-1 elimination (h-1)
k_cysC_prod= 0.05   // Cystatin C production rate (mg/L/h)
k_cysC_deg = 0.03   // Cystatin C elimination (h-1)
CysC0      = 0.8    // Baseline cystatin C (mg/L)
k_scr_prod = 0.005  // Creatinine endogenous production (mg/dL/h)

// ─── Repair Pathway ───
k_repair   = 0.003  // Tubular repair rate (h-1)
k_egf      = 0.1    // EGF/HGF-driven proliferation (h-1)
EGF0       = 1.0    // Baseline EGF signaling
k_wnt      = 0.05   // Wnt/β-catenin repair signal

// ─── Fibrosis (AKI-to-CKD) ───
k_tgfb_prod= 0.02   // TGF-β1 production rate (h-1)
k_tgfb_deg = 0.05   // TGF-β1 degradation (h-1)
k_fibrosis = 0.001  // Fibrosis rate from TGF-β (h-1)
TGFb0      = 1.0    // Baseline TGF-β (pg/mL)
k_myo_act  = 0.015  // Myofibroblast activation (h-1)
k_myo_deg  = 0.02   // Myofibroblast turnover (h-1)

// ─── Furosemide PK ───
Vd_fur1    = 8      // Central volume (L)
Vd_fur2    = 14     // Peripheral volume (L)
CL_fur     = 8      // Clearance (L/h)
Q_fur      = 4      // Inter-compartment clearance (L/h)
F_fur      = 0.6    // Oral bioavailability
OAT_eff    = 1.0    // OAT1/3 secretion efficiency (reduced in AKI)
NKCC2_IC50 = 0.5    // Furosemide IC50 for NKCC2 (mg/L)

// ─── Norepinephrine PK ───
CL_ne      = 250    // Clearance NE (L/h, very fast t1/2 ~2.5 min)
Vd_ne      = 10     // Volume NE (L)
Emax_map   = 0.4    // Maximum MAP increase fraction
EC50_ne    = 0.01   // NE EC50 for MAP (mcg/mL)

// ─── N-Acetylcysteine PK ───
CL_nac     = 15     // NAC clearance (L/h)
Vd_nac     = 30     // NAC volume (L)
Emax_gsh   = 0.6    // Max GSH replenishment by NAC
EC50_nac   = 10     // NAC EC50 for GSH (mg/L)
Emax_no    = 0.3    // Max NO enhancement by NAC

// ─── CRRT Parameters ───
CRRT_on    = 0      // 1=CRRT active
CL_crrt    = 3      // CRRT creatinine clearance (L/h, ~50 mL/min)
CL_crrt_il6= 0.5    // CRRT IL-6 clearance (L/h, adsorption)

$CMT
// Drug PK compartments
FUR_GUT     // Furosemide gut (oral)
FUR_CENT    // Furosemide central
FUR_PERI    // Furosemide peripheral
NE_CENT     // Norepinephrine central
NAC_CENT    // NAC central

// AKI Pathophysiology
ATP         // ATP level (normalized, 0=depleted, 1=normal)
ROS         // Reactive oxygen species (normalized)
GSH         // Glutathione (normalized)
TCV         // Tubular cell viability (0-1)
GFR         // Glomerular filtration rate (mL/min/1.73m2)
IL6         // IL-6 (pg/mL)
TNFa        // TNF-alpha (pg/mL)

// Biomarkers
NGAL        // Urine/plasma NGAL (ng/mL)
KIM1        // Urine KIM-1 (ng/mL)
SCR         // Serum creatinine (mg/dL)
CysC        // Serum cystatin C (mg/L)

// Repair & Fibrosis
REPAIR_CAP  // Repair capacity (EGF/HGF signaling, normalized)
TGFb        // TGF-β1 (pg/mL)
MYO         // Myofibroblast burden (normalized)
FIBROSIS    // Cumulative fibrosis index (0-1)

$MAIN
// Trigger injury stimulus
double inj_stim = 0;
if(AKI_type == 1) inj_stim = isch_sev;       // IRI
if(AKI_type == 2) inj_stim = tox_dose / 20;  // NTX (20 mg/kg ~ full injury)
if(AKI_type == 3) inj_stim = lps_amt / 10;   // Sepsis (10 ng/mL equiv)

// Furosemide effect on urine output (OAT-mediated tubular secretion)
double FUR_TUB = FUR_CENT * OAT_eff * TCV; // effective tubular conc
double E_NKCC2 = FUR_TUB / (FUR_TUB + NKCC2_IC50);
double diuresis = 1.0 + 3.0 * E_NKCC2;    // diuresis fold over baseline

// Norepinephrine: MAP and RBF effect
double E_NE = Emax_map * NE_CENT / (NE_CENT + EC50_ne);
double RBF  = RBF0 * (1 + E_NE);           // MAP↑ → RBF autoregulation
double E_RBF_gfr = RBF / RBF0;             // relative RBF effect on GFR

// NAC effects
double E_NAC_gsh = Emax_gsh * NAC_CENT / (NAC_CENT + EC50_nac);
double E_NAC_no  = Emax_no  * NAC_CENT / (NAC_CENT + EC50_nac);

// Injury composite (ATP depletion + ROS-driven Hill)
double ROS_frac = ROS / (ROS + EC50_ros);
double ATP_frac = (1.0 - ATP) / (1.0 - ATP + EC50_atp);
double inj_composite = (0.5 * ROS_frac + 0.5 * ATP_frac);

// Repair stimulus from surviving cells
double repair_stim = TCV * REPAIR_CAP * k_egf;

// NF-κB driver
double nfkb = k_nfkb * (inj_stim + ROS + (1.0 - TCV));

// GPx4 activity (GSH-dependent ferroptosis gate)
double GPx4_act = GSH;
double ferr_rate = k_ferr * (1 - GPx4_act) * inj_composite;

// Creatinine production (endogenous, independent of GFR)
double CR_prod  = k_scr_prod;
// Creatinine elimination driven by GFR + CRRT
double GFR_norm = GFR / GFR0;
double CR_elim  = (GFR_norm * 0.08 + CRRT_on * CL_crrt / 100.0) * SCR;

// AKI staging (KDIGO) based on sCr fold change
double scr_fold = SCR / SCR0;
int AKI_STAGE;
if(scr_fold < 1.5)      AKI_STAGE = 0;
else if(scr_fold < 2.0) AKI_STAGE = 1;
else if(scr_fold < 3.0) AKI_STAGE = 2;
else                     AKI_STAGE = 3;

$ODE
// ─── Furosemide PK ───
dxdt_FUR_GUT  = -KA_FUR * FUR_GUT;  // KA defined by user or fixed
dxdt_FUR_CENT = F_fur * KA_FUR * FUR_GUT
                - (CL_fur / Vd_fur1) * FUR_CENT
                - (Q_fur / Vd_fur1) * FUR_CENT
                + (Q_fur / Vd_fur2) * FUR_PERI;
dxdt_FUR_PERI = (Q_fur / Vd_fur1) * FUR_CENT
                - (Q_fur / Vd_fur2) * FUR_PERI;

// ─── Norepinephrine PK ───
dxdt_NE_CENT  = - (CL_ne / Vd_ne) * NE_CENT;  // continuous infusion via RATE

// ─── NAC PK ───
dxdt_NAC_CENT = - (CL_nac / Vd_nac) * NAC_CENT;

// ─── ATP ───
dxdt_ATP = k_atp_prod * (1 - inj_stim) * (1 - FIBROSIS)
           - k_atp_dep * (1 + inj_stim + 0.5 * (1 - TCV));
// clamp implicitly via init

// ─── Reactive Oxygen Species ───
dxdt_ROS = k_ros_prod * (1 + 2 * inj_stim + (1 - TCV) + (1 - ATP))
           - k_ros_deg * GSH * ROS
           - 0.2 * E_NAC_no * ROS;  // NAC scavenging

// ─── Glutathione ───
dxdt_GSH = k_gsh_prod * (1 + E_NAC_gsh)
           - k_gsh_dep * ROS * GSH;

// ─── Tubular Cell Viability ───
double inj_rate_apop = k_apop * inj_composite * TCV * pow(ROS, n_hill);
double inj_rate_nec  = k_nec  * (1 - ATP) * TCV;
double inj_rate_ferr = ferr_rate * TCV;
dxdt_TCV = -inj_rate_apop - inj_rate_nec - inj_rate_ferr
           + repair_stim * (1 - TCV) * (1 - FIBROSIS);

// ─── GFR ───
// GFR declines with tubular injury (obstruction + backleak), recovers as TCV recovers
double gfr_inj  = k_gfr_inj * (1 - TCV) * (1 + k_backleak + k_obstruct);
double gfr_rec  = k_gfr_rec * TCV * E_RBF_gfr * (GFR0 - GFR);
dxdt_GFR = gfr_rec - gfr_inj * GFR;

// ─── IL-6 ───
dxdt_IL6 = k_il6_prod * (1 + 5 * nfkb + lps_amt / 5)
           - k_il6_deg * IL6
           - CRRT_on * CL_crrt_il6 * IL6 / 5.0;

// ─── TNF-alpha ───
dxdt_TNFa = k_tnfa_prod * (1 + 3 * nfkb + lps_amt / 8)
            - k_tnfa_deg * TNFa;

// ─── Biomarkers ───
// NGAL (injury marker, rises within 2h)
dxdt_NGAL = k_ngal_prod * (1 - TCV) * (1 + 3 * inj_stim)
            - k_ngal_deg * NGAL;

// KIM-1 (brush border shedding, slower)
dxdt_KIM1 = k_kim1_prod * (1 - TCV) * (1 + 2 * inj_stim)
            - k_kim1_deg * KIM1;

// Serum Creatinine
dxdt_SCR = CR_prod - CR_elim;

// Cystatin C
double GFR_n2 = GFR / GFR0;
dxdt_CysC = k_cysC_prod - k_cysC_deg * GFR_n2 * CysC;

// ─── Repair Capacity ───
// EGF/HGF signaling declines with injury, recovers when inflammation resolves
dxdt_REPAIR_CAP = k_repair * (1 - FIBROSIS) * (1 + k_wnt)
                  - 0.05 * (IL6 / (IL6 + 20)) * REPAIR_CAP
                  - 0.03 * (1 - TCV) * REPAIR_CAP;

// ─── TGF-β1 ───
// Produced by macrophages, hypoxic tubular cells; driven by prolonged injury
dxdt_TGFb = k_tgfb_prod * (1 + 3 * (1 - TCV) + nfkb)
             - k_tgfb_deg * TGFb;

// ─── Myofibroblast Burden ───
dxdt_MYO = k_myo_act * TGFb * (1 - MYO) - k_myo_deg * MYO;

// ─── Cumulative Fibrosis ───
dxdt_FIBROSIS = k_fibrosis * MYO * (1 - FIBROSIS);

$TABLE
double AKI_stage_out = 0;
double scr_ratio = SCR / SCR0;
if(scr_ratio >= 1.5) AKI_stage_out = 1;
if(scr_ratio >= 2.0) AKI_stage_out = 2;
if(scr_ratio >= 3.0) AKI_stage_out = 3;

double eGFR = GFR;
double FE_Na = 100.0 * (1.0 - TCV * 0.9);  // FENa proxy (tubular NaK-ATPase loss)
double UO_mLkgh = 1.0 * diuresis * (TCV * 0.8 + 0.2) * (GFR / GFR0);  // mL/kg/h approx

$CAPTURE
AKI_stage_out eGFR UO_mLkgh FE_Na
FUR_CENT NE_CENT NAC_CENT
ATP ROS GSH TCV
IL6 TNFa NGAL KIM1 SCR CysC
TGFb MYO FIBROSIS REPAIR_CAP
diuresis E_NKCC2 E_NE E_NAC_gsh
'

## ============================================================
## Auxiliary: KA for furosemide (defined in $PARAM via macro)
## We inject KA_FUR as a parameter in the model
## ============================================================

code2 <- gsub("dxdt_FUR_GUT  = -KA_FUR", "dxdt_FUR_GUT  = -1.5", code)
code2 <- gsub("dxdt_FUR_CENT = F_fur \\* KA_FUR", "dxdt_FUR_CENT = F_fur * 1.5", code2)

## ============================================================
## Compile Model
## ============================================================

mod <- mcode("AKI_QSP", code2)

cat("=== AKI QSP Model compiled successfully ===\n")
cat("Compartments:", length(Init(mod)), "\n")

## ============================================================
## INITIAL CONDITIONS
## ============================================================

mod <- mod %>%
  init(
    FUR_GUT  = 0,
    FUR_CENT = 0,
    FUR_PERI = 0,
    NE_CENT  = 0,
    NAC_CENT = 0,
    ATP      = 1.0,
    ROS      = 0.1,
    GSH      = 1.0,
    TCV      = 1.0,
    GFR      = 100,
    IL6      = 5,
    TNFa     = 3,
    NGAL     = 5,    # baseline ~5 ng/mL
    KIM1     = 0.5,  # baseline
    SCR      = 0.9,
    CysC     = 0.8,
    REPAIR_CAP = 1.0,
    TGFb     = 1.0,
    MYO      = 0.0,
    FIBROSIS = 0.0
  )

## ============================================================
## SCENARIO DEFINITIONS
## ============================================================

## Scenario 1: Natural course — IRI (moderate, no treatment)
scen1 <- mod %>%
  param(AKI_type=1, isch_sev=0.6) %>%
  mrgsim(end=168, delta=0.5) %>%  # 7 days
  as_tibble() %>%
  mutate(Scenario="IRI - No Treatment")

## Scenario 2: IRI + Furosemide (40mg IV q12h)
ev2 <- ev(amt=40, ii=12, addl=13, cmt="FUR_GUT", time=2)  # start at 2h
scen2 <- mod %>%
  param(AKI_type=1, isch_sev=0.6) %>%
  mrgsim(events=ev2, end=168, delta=0.5) %>%
  as_tibble() %>%
  mutate(Scenario="IRI + Furosemide 40mg IV q12h")

## Scenario 3: IRI + NE support (MAP target) + Fluid
ev3 <- ev(amt=0.1, rate=0.1, cmt="NE_CENT", time=0, evid=1)  # NE continuous
scen3 <- mod %>%
  param(AKI_type=1, isch_sev=0.6) %>%
  mrgsim(events=ev3, end=168, delta=0.5) %>%
  as_tibble() %>%
  mutate(Scenario="IRI + Vasopressor (NE) Support")

## Scenario 4: IRI + NAC prophylaxis (IV 150 mg/kg loading then maintenance)
ev4 <- ev(amt=10000, rate=1000, cmt="NAC_CENT", time=-2)  # pre-injury NAC
scen4 <- mod %>%
  param(AKI_type=1, isch_sev=0.6) %>%
  mrgsim(events=ev4, end=168, delta=0.5) %>%
  as_tibble() %>%
  mutate(Scenario="IRI + NAC Prophylaxis")

## Scenario 5: Sepsis-associated AKI — natural course
scen5 <- mod %>%
  param(AKI_type=3, lps_amt=8, isch_sev=0.3) %>%
  mrgsim(end=168, delta=0.5) %>%
  as_tibble() %>%
  mutate(Scenario="SA-AKI - No Treatment")

## Scenario 6: Sepsis-AKI + NE + Furosemide + Early CRRT
ev6a <- ev(amt=0.1, rate=0.1, cmt="NE_CENT", time=0)
ev6b <- ev(amt=40, ii=12, addl=13, cmt="FUR_GUT", time=6)
ev6  <- ev6a + ev6b
scen6 <- mod %>%
  param(AKI_type=3, lps_amt=8, isch_sev=0.3, CRRT_on=1) %>%
  mrgsim(events=ev6, end=168, delta=0.5) %>%
  as_tibble() %>%
  mutate(Scenario="SA-AKI + NE + Furosemide + CRRT")

## Scenario 7: Nephrotoxic AKI (cisplatin 75 mg/m² equiv)
scen7 <- mod %>%
  param(AKI_type=2, tox_dose=75) %>%
  mrgsim(end=336, delta=1) %>%  # 14 days
  as_tibble() %>%
  mutate(Scenario="Nephrotoxic AKI (Cisplatin)")

## ============================================================
## COMBINE & PLOT
## ============================================================

all_scen <- bind_rows(scen1, scen2, scen3, scen4, scen5, scen6)

## Plot 1: GFR over time
p1 <- ggplot(all_scen, aes(time, eGFR, color=Scenario)) +
  geom_line(size=1.2) +
  geom_hline(yintercept=60, linetype="dashed", color="gray40") +
  annotate("text", x=5, y=62, label="CKD threshold (eGFR 60)", size=3, color="gray40") +
  scale_x_continuous(breaks=seq(0,168,24), labels=paste0(seq(0,168,24)/24, "d")) +
  labs(title="AKI: GFR Dynamics by Treatment Scenario",
       x="Time (days)", y="eGFR (mL/min/1.73m²)", color="Scenario") +
  theme_bw(base_size=12) +
  theme(legend.position="bottom", legend.text=element_text(size=8)) +
  guides(color=guide_legend(ncol=2))

## Plot 2: Serum Creatinine + AKI staging
p2 <- ggplot(all_scen, aes(time, SCR, color=Scenario)) +
  geom_line(size=1.2) +
  geom_hline(yintercept=c(0.9*1.5, 0.9*2.0, 0.9*3.0),
             linetype=c("dotted","dashed","solid"), color="red", alpha=0.6) +
  annotate("text", x=160, y=c(1.35+0.06, 1.8+0.06, 2.7+0.06),
           label=c("AKI Stage 1","Stage 2","Stage 3"), color="red", size=3, hjust=1) +
  scale_x_continuous(breaks=seq(0,168,24), labels=paste0(seq(0,168,24)/24, "d")) +
  labs(title="Serum Creatinine Trajectory (KDIGO AKI Staging)",
       x="Time", y="sCr (mg/dL)", color="Scenario") +
  theme_bw(base_size=12) +
  theme(legend.position="bottom") +
  guides(color=guide_legend(ncol=2))

## Plot 3: NGAL & KIM-1 (early biomarkers)
bm_long <- all_scen %>%
  select(time, Scenario, NGAL, KIM1) %>%
  pivot_longer(cols=c(NGAL, KIM1), names_to="Biomarker", values_to="Level")

p3 <- ggplot(filter(bm_long, time <= 48), aes(time, Level, color=Scenario, linetype=Biomarker)) +
  geom_line(size=1.1) +
  labs(title="Early AKI Biomarkers (First 48h)",
       x="Time (hours)", y="Biomarker Level (ng/mL)", color="Scenario") +
  theme_bw(base_size=12) +
  theme(legend.position="bottom") +
  guides(color=guide_legend(ncol=2))

## Plot 4: Inflammatory cascade (IL-6, TNF-α)
inf_long <- all_scen %>%
  select(time, Scenario, IL6, TNFa) %>%
  pivot_longer(cols=c(IL6, TNFa), names_to="Cytokine", values_to="Level")

p4 <- ggplot(inf_long, aes(time, Level, color=Scenario, linetype=Cytokine)) +
  geom_line(size=1.1) +
  scale_x_continuous(breaks=seq(0,168,24), labels=paste0(seq(0,168,24)/24, "d")) +
  labs(title="Inflammatory Cytokines (IL-6, TNF-α)",
       x="Time", y="Cytokine (pg/mL)", color="Scenario", linetype="Cytokine") +
  theme_bw(base_size=12) +
  theme(legend.position="bottom") +
  guides(color=guide_legend(ncol=2))

## Plot 5: Tubular cell viability + ROS + GSH
cell_long <- all_scen %>%
  select(time, Scenario, TCV, ROS, GSH) %>%
  pivot_longer(cols=c(TCV, ROS, GSH), names_to="Variable", values_to="Value")

p5 <- ggplot(cell_long, aes(time, Value, color=Scenario, linetype=Variable)) +
  geom_line(size=1) +
  scale_x_continuous(breaks=seq(0,168,24), labels=paste0(seq(0,168,24)/24, "d")) +
  labs(title="Tubular Cell Biology: Viability, ROS, GSH",
       x="Time", y="Normalized Level", color="Scenario", linetype="Variable") +
  theme_bw(base_size=12) +
  theme(legend.position="bottom") +
  guides(color=guide_legend(ncol=2))

## Plot 6: Cisplatin NTX fibrosis trajectory
p6 <- ggplot(scen7, aes(time/24)) +
  geom_line(aes(y=FIBROSIS*100, color="Fibrosis Index (%)"), size=1.2) +
  geom_line(aes(y=TGFb*5, color="TGF-β1 (×5)"), size=1) +
  geom_line(aes(y=MYO*100, color="Myofibroblast (%)"), size=1, linetype="dashed") +
  labs(title="Nephrotoxic AKI: AKI-to-CKD Transition (Cisplatin)",
       x="Time (days)", y="Level", color="Variable") +
  theme_bw(base_size=12) +
  theme(legend.position="right")

## Print all plots
print(p1); print(p2); print(p3); print(p4); print(p5); print(p6)

## ============================================================
## CLINICAL SUMMARY TABLE
## ============================================================

summary_tbl <- all_scen %>%
  group_by(Scenario) %>%
  summarise(
    Peak_sCr    = round(max(SCR), 2),
    Max_AKI_Stage = max(AKI_stage_out),
    Nadir_GFR   = round(min(eGFR), 1),
    Peak_NGAL   = round(max(NGAL), 1),
    Peak_IL6    = round(max(IL6), 1),
    Min_TCV     = round(min(TCV), 3),
    Final_GFR   = round(last(eGFR), 1),
    Final_Fibr  = round(last(FIBROSIS), 4),
    .groups     = "drop"
  )

cat("\n=== Clinical Summary Table ===\n")
print(as.data.frame(summary_tbl))

cat("\n=== AKI QSP Model: Run Complete ===\n")
cat("Scenarios: 7 (IRI alone, IRI+Furosemide, IRI+NE, IRI+NAC,",
    "SA-AKI, SA-AKI+combo, NTX)\n")
cat("Endpoint coverage: GFR, sCr, NGAL, KIM-1, IL-6, TNF-α, TGF-β, fibrosis\n")
