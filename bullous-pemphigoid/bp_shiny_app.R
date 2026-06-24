## =============================================================================
## Bullous Pemphigoid QSP — Interactive Shiny Dashboard
## Author: Claude Code Routine (CCR)
## Date: 2026-06-19
##
## Tabs:
##   1. Patient Profile & Disease Overview
##   2. Drug PK (prednisolone, dupilumab, omalizumab, rituximab, doxycycline)
##   3. PD Key Biomarkers (anti-BP180 IgG/IgE, eosinophils, mast cells, C5a)
##   4. Clinical Endpoints (BPDAI, Itch NRS, BSA, DEJ integrity)
##   5. Treatment Scenario Comparison
##   6. Mechanistic Biomarker Explorer (heatmap + DEJ pathway)
## =============================================================================

library(shiny)
library(bslib)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(DT)

## ── Inline mrgsolve model (abbreviated for Shiny; same ODEs) ─────────────────
bp_code <- '
$PARAM
ka_pred=2.4 F_pred=0.82 Vd_pred=45 CL_pred=6.5 Vp_pred=15 Q_pred=1.2
ka_dup=0.022 F_dup=0.64 Vd_dup=4.8 CL_dup=0.0064
ka_oma=0.010 F_oma=0.62 Vd_oma=7.7 CL_oma=0.0112
Vd_ritu=3.6 CL_ritu=0.0162
ka_doxy=0.58 F_doxy=0.93 Vd_doxy=150 CL_doxy=1.8
kBn_prod=0.006 kBn_die=0.0012 kBact=0.0003 kGC=0.0008
kMem=0.0004 kSLPC=0.0005 kLLPC=0.00004
kBact_die=0.004 kMem_die=0.00005 kSLPC_die=0.012 kLLPC_die=0.00014
kTh2_prod=0.0008 kTh2_die=0.002 Th2_0=80
kIgG_prod=0.00018 kIgG_deg=0.0013 kIgE_prod=0.00006 kIgE_deg=0.0046
kEos_prod=0.0004 kEos_die=0.0017 kEos_skin=0.0008 kEos_sk_die=0.004
kMast_act=0.0030 kMast_base=0.0005 kMast_die=0.0050
kC5a_prod=0.0020 kC5a_deg=0.030
kDEJ_dam=0.0015 kDEJ_repair=0.0006
EC50_pred_IL4=0.015 EC50_pred_Eos=0.010 EC50_dup_IL4=0.008
EC50_oma_IgE=0.030 EC50_ritu_B=0.015 EC50_doxy_MMP=2.00 Emax=0.90
BP180_Ag_level=1.0 IL4_stim=1.5 Th2_bias=2.0

$CMT
PRED_GUT PRED_C PRED_P DUP_DEPOT DUP_C OMA_DEPOT OMA_C RITU_C DOXY_GUT DOXY_C
B_NAIVE B_ACT B_MEM SLPC LLPC TH2 IGG_BP180 IGE_BP180
EOS_BLOOD EOS_SKIN MAST_ACT C5A DEJ

$INIT
B_NAIVE=200 B_ACT=10 B_MEM=20 SLPC=5 LLPC=2 TH2=80
IGG_BP180=1.0 IGE_BP180=0.5 EOS_BLOOD=0.35 EOS_SKIN=0.20
MAST_ACT=0.10 C5A=0.05 DEJ=0.50

$ODE
dxdt_PRED_GUT = -ka_pred*PRED_GUT;
dxdt_PRED_C = ka_pred*F_pred*PRED_GUT - (CL_pred/Vd_pred+Q_pred/Vd_pred)*PRED_C + (Q_pred/Vp_pred)*PRED_P;
dxdt_PRED_P = (Q_pred/Vd_pred)*PRED_C - (Q_pred/Vp_pred)*PRED_P;
dxdt_DUP_DEPOT = -ka_dup*DUP_DEPOT;
dxdt_DUP_C = ka_dup*F_dup*DUP_DEPOT - (CL_dup/Vd_dup)*DUP_C;
dxdt_OMA_DEPOT = -ka_oma*OMA_DEPOT;
dxdt_OMA_C = ka_oma*F_oma*OMA_DEPOT - (CL_oma/Vd_oma)*OMA_C;
dxdt_RITU_C = -(CL_ritu/Vd_ritu)*RITU_C;
dxdt_DOXY_GUT = -ka_doxy*DOXY_GUT;
dxdt_DOXY_C = ka_doxy*F_doxy*DOXY_GUT - (CL_doxy/Vd_doxy)*DOXY_C;
double Cpred = PRED_C/Vd_pred;
double Cdup  = DUP_C/Vd_dup;
double Coma  = OMA_C/Vd_oma;
double Critu = RITU_C/Vd_ritu;
double Cdoxy = DOXY_C/Vd_doxy;
double E_pred_IL4 = Emax*Cpred/(EC50_pred_IL4+Cpred);
double E_pred_Eos = Emax*Cpred/(EC50_pred_Eos+Cpred);
double E_dup_IL4  = Emax*Cdup/(EC50_dup_IL4+Cdup);
double E_oma_IgE  = Emax*Coma/(EC50_oma_IgE+Coma);
double E_ritu_B   = Emax*Critu/(EC50_ritu_B+Critu);
double E_doxy_MMP = Emax*Cdoxy/(EC50_doxy_MMP+Cdoxy);
double E_Th2_total = 1.0-(1.0-E_pred_IL4)*(1.0-E_dup_IL4);
double Th2_drive = kTh2_prod*Th2_bias*IL4_stim*BP180_Ag_level*(1.0-E_Th2_total);
dxdt_TH2 = Th2_drive*B_NAIVE/(B_NAIVE+50.0) - kTh2_die*TH2;
double BCR_stim = IGG_BP180*BP180_Ag_level;
double Bcell_dep = E_ritu_B;
dxdt_B_NAIVE = kBn_prod - kBn_die*B_NAIVE - kBact*BCR_stim*B_NAIVE*(1.0-Bcell_dep);
dxdt_B_ACT = kBact*BCR_stim*B_NAIVE*(1.0-Bcell_dep) - (kGC+kBact_die)*B_ACT*(1.0-Bcell_dep);
dxdt_B_MEM = kMem*B_ACT - kMem_die*B_MEM*(1.0-Bcell_dep);
dxdt_SLPC = kSLPC*B_ACT - (kSLPC_die+kLLPC)*SLPC;
dxdt_LLPC = kLLPC*SLPC - kLLPC_die*LLPC;
double IgG_prod_rate = kIgG_prod*(LLPC+0.3*SLPC)*(1.0-0.6*E_ritu_B);
dxdt_IGG_BP180 = IgG_prod_rate - kIgG_deg*IGG_BP180;
double IgE_prod_rate = kIgE_prod*SLPC*(1.0-E_Th2_total);
double IgE_neutralise = E_oma_IgE*0.8*IGE_BP180;
dxdt_IGE_BP180 = IgE_prod_rate - kIgE_deg*IGE_BP180 - IgE_neutralise;
double Eos_stim = (TH2/Th2_0)*IL4_stim*(1.0-E_pred_Eos);
dxdt_EOS_BLOOD = kEos_prod*Eos_stim - kEos_die*EOS_BLOOD - kEos_skin*IGG_BP180*EOS_BLOOD;
double Skin_recruit = kEos_skin*IGG_BP180*EOS_BLOOD*(1.0-E_dup_IL4*0.5);
dxdt_EOS_SKIN = Skin_recruit - kEos_sk_die*EOS_SKIN;
double Mast_drive = kMast_act*IGE_BP180*(1.0-E_oma_IgE) + kMast_base;
dxdt_MAST_ACT = Mast_drive - kMast_die*MAST_ACT;
double C5a_drive = kC5a_prod*IGG_BP180*(1.0-E_pred_IL4*0.3);
dxdt_C5A = C5a_drive - kC5a_deg*C5A;
double DEJ_damage_rate = kDEJ_dam*(EOS_SKIN+0.5*C5A+0.3*MAST_ACT)*(1.0-E_doxy_MMP);
double DEJ_repair_rate = kDEJ_repair*(1.0-DEJ);
dxdt_DEJ = DEJ_repair_rate - DEJ_damage_rate*DEJ;

$TABLE
double BPDAI_act = 30.0*(1.0-DEJ)*(IGG_BP180/(IGG_BP180+0.5));
if(BPDAI_act < 0) BPDAI_act = 0;
if(BPDAI_act > 90) BPDAI_act = 90;
double itch_raw = 6.0*MAST_ACT/0.10*IGE_BP180/0.5*(1.0/(1.0+(DUP_C/Vd_dup)/EC50_dup_IL4));
double Itch = (itch_raw < 0) ? 0 : (itch_raw > 10 ? 10 : itch_raw);
double BSA = 20.0*(1.0-DEJ)*(IGG_BP180/(IGG_BP180+1.0));
double new_blisters = BPDAI_act/10.0;
double in_remission = (new_blisters < 0.3) ? 1.0 : 0.0;
double Cpred_ngml = PRED_C/Vd_pred*1000.0;
double Cdup_mgL = DUP_C/Vd_dup;
double Coma_mgL = OMA_C/Vd_oma;
double Critu_mgL = RITU_C/Vd_ritu;
double Cdoxy_mcg = DOXY_C/Vd_doxy*1000.0;

$CAPTURE BPDAI_act Itch BSA new_blisters in_remission
Cpred_ngml Cdup_mgL Coma_mgL Critu_mgL Cdoxy_mcg
IGG_BP180 IGE_BP180 EOS_BLOOD EOS_SKIN MAST_ACT C5A DEJ
B_NAIVE B_ACT LLPC TH2
'

mod_shiny <- mread_cache("bp_shiny", code = bp_code)

## ── Helper: Run Single Scenario ───────────────────────────────────────────────
run_bp_sim <- function(
    pred_dose = 0.5,   # mg/kg/d
    bw = 65,
    use_dup = FALSE,
    use_oma = FALSE,
    use_ritu = FALSE,
    use_doxy = FALSE,
    weeks = 52) {

  hours <- weeks * 7 * 24

  ev_list <- list()

  # Prednisolone taper
  if (pred_dose > 0) {
    dose_mg <- pred_dose * bw
    ev_list[[length(ev_list)+1]] <- ev(cmt="PRED_GUT", amt=dose_mg, ii=24,
                                        addl=4*7-1, time=0)
    dose2 <- dose_mg * 0.5
    ev_list[[length(ev_list)+1]] <- ev(cmt="PRED_GUT", amt=dose2, ii=24,
                                        addl=4*7-1, time=4*7*24)
    ev_list[[length(ev_list)+1]] <- ev(cmt="PRED_GUT", amt=max(dose2*0.5,5), ii=24,
                                        addl=weeks*7-1, time=8*7*24)
  }

  if (use_dup) {
    ev_list[[length(ev_list)+1]] <- ev(cmt="DUP_DEPOT", amt=600, time=0)
    ev_list[[length(ev_list)+1]] <- ev(cmt="DUP_DEPOT", amt=300, ii=14*24,
                                        addl=weeks/2-1, time=14*24)
  }

  if (use_oma) {
    ev_list[[length(ev_list)+1]] <- ev(cmt="OMA_DEPOT", amt=300, ii=28*24,
                                        addl=floor(weeks/4), time=0)
  }

  if (use_ritu) {
    ev_list[[length(ev_list)+1]] <- ev(cmt="RITU_C", amt=1000, time=0)
    ev_list[[length(ev_list)+1]] <- ev(cmt="RITU_C", amt=1000, time=14*24)
    if (weeks > 26)
      ev_list[[length(ev_list)+1]] <- ev(cmt="RITU_C", amt=1000, time=26*7*24)
  }

  if (use_doxy) {
    ev_list[[length(ev_list)+1]] <- ev(cmt="DOXY_GUT", amt=200, ii=24,
                                        addl=weeks*7-1, time=0)
  }

  # Combine events
  total_ev <- Reduce(`+`, ev_list)

  out <- mod_shiny %>%
    mrgsim_df(events = total_ev, end = hours, delta = 24, hmax = 1) %>%
    mutate(time_weeks = time / (7 * 24))

  return(out)
}

## ── UI ────────────────────────────────────────────────────────────────────────
ui <- page_navbar(
  title = "Bullous Pemphigoid QSP Model",
  theme = bs_theme(bootswatch = "flatly", base_font = font_google("Inter")),

  ## ── TAB 1: Patient Profile ──────────────────────────────────────────────────
  nav_panel(
    "Patient Profile",
    icon = icon("user-md"),
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        h5("Patient Parameters"),
        sliderInput("bw", "Body Weight (kg)", 40, 120, 65, step = 5),
        sliderInput("weeks", "Simulation Duration (weeks)", 12, 104, 52, step = 4),
        hr(),
        h5("Disease Status"),
        selectInput("severity", "Initial Disease Severity",
                    c("Mild (BPDAI<20)"="mild",
                      "Moderate (BPDAI 20-56)"="moderate",
                      "Severe (BPDAI>56)"="severe"),
                    selected = "moderate"),
        checkboxInput("neuro_comorbid", "Neurological Comorbidity", FALSE),
        checkboxInput("drug_induced", "Drug-induced BP (Checkpoint inhibitor/Gliptin)", FALSE),
        hr(),
        h5("Disease Pathology Info"),
        tags$small(tags$b("Bullous Pemphigoid (BP)"), " is the most common autoimmune blistering disease, affecting primarily elderly patients (>65 years). Key features:"),
        tags$ul(
          tags$li("Th2-dominant immune response"),
          tags$li("Autoantibodies against BP180 (Col XVII) and BP230"),
          tags$li("Subepidermal blister formation"),
          tags$li("Intense pruritus (itch)"),
          tags$li("Annual incidence: 21-66 per million"),
          tags$li("1-year mortality: 10-30% in elderly")
        )
      ),
      layout_columns(
        col_widths = c(4, 4, 4),
        value_box(
          title = "Disease Mechanism",
          value = "Th2/IgE-mediated",
          showcase = icon("dna"),
          theme = "danger"
        ),
        value_box(
          title = "Primary Target",
          value = "BP180 NC16A",
          showcase = icon("bullseye"),
          theme = "warning"
        ),
        value_box(
          title = "Peak Incidence",
          value = ">70 years",
          showcase = icon("person-cane"),
          theme = "info"
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Pathomechanism Summary"),
          card_body(
            tags$table(
              class = "table table-sm table-bordered",
              tags$thead(tags$tr(
                tags$th("Phase"), tags$th("Key Events"), tags$th("Biomarkers")
              )),
              tags$tbody(
                tags$tr(tags$td("Sensitisation"), tags$td("BP180 antigen presentation, Th2 polarisation"), tags$td("IL-4↑, TSLP↑")),
                tags$tr(tags$td("B Cell Activation"), tags$td("GC reaction, IgG/IgE class switch"), tags$td("Anti-BP180 IgG/IgE↑")),
                tags$tr(tags$td("Effector Phase"), tags$td("Complement activation, eosinophil/mast cell recruitment"), tags$td("C5a↑, Eosinophils↑")),
                tags$tr(tags$td("Tissue Damage"), tags$td("DEJ cleavage, blister formation"), tags$td("BPDAI↑, BSA↑")),
                tags$tr(tags$td("Pruritus"), tags$td("IL-31R signaling, histamine release"), tags$td("Itch NRS↑")),
                tags$tr(tags$td("Remission"), tags$td("Antibody clearance, DEJ repair"), tags$td("BPDAI=0"))
              )
            )
          )
        ),
        card(
          card_header("Clinical Endpoints"),
          card_body(
            tags$dl(
              tags$dt("BPDAI"), tags$dd("Bullous Pemphigoid Disease Area Index (0–90): counts blisters, erosions, urticaria, erythema weighted by BSA"),
              tags$dt("Itch NRS"), tags$dd("Numerical Rating Scale for pruritus (0–10)"),
              tags$dt("IGA-BP"), tags$dd("Investigator Global Assessment for BP (0–5)"),
              tags$dt("DLQI"), tags$dd("Dermatology Life Quality Index (0–30)"),
              tags$dt("Anti-BP180 AU"), tags$dd("ELISA titer against NC16A domain; >9 AU/mL diagnostic")
            )
          )
        )
      )
    )
  ),

  ## ── TAB 2: Drug PK ──────────────────────────────────────────────────────────
  nav_panel(
    "Drug PK",
    icon = icon("pills"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("Treatment Selection"),
        numericInput("pred_dose_pk", "Prednisolone (mg/kg/d)", 0.5, 0, 1.5, 0.1),
        checkboxInput("use_dup_pk", "Dupilumab 300mg SC Q2W", FALSE),
        checkboxInput("use_oma_pk", "Omalizumab 300mg SC Q4W", FALSE),
        checkboxInput("use_ritu_pk", "Rituximab 1g IV (protocol)", FALSE),
        checkboxInput("use_doxy_pk", "Doxycycline 200mg/d", FALSE),
        sliderInput("bw_pk", "Body Weight (kg)", 40, 120, 65, 5),
        sliderInput("weeks_pk", "Duration (weeks)", 12, 104, 52, 4),
        actionButton("run_pk", "Run PK Simulation", class = "btn-primary btn-block"),
        hr(),
        tags$small(
          tags$b("PK Parameters:"),
          tags$br(), "Prednisolone: F=82%, t½=3.5h",
          tags$br(), "Dupilumab: F=64%, t½=21d",
          tags$br(), "Omalizumab: F=62%, t½=26d",
          tags$br(), "Rituximab: IV, t½=22d",
          tags$br(), "Doxycycline: F=93%, t½=16h"
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Prednisolone Plasma Concentration"),
             card_body(plotlyOutput("pk_pred", height = "300px"))),
        card(card_header("Dupilumab Concentration"),
             card_body(plotlyOutput("pk_dup", height = "300px")))
      ),
      layout_columns(
        col_widths = c(4, 4, 4),
        card(card_header("Omalizumab"),
             card_body(plotlyOutput("pk_oma", height = "260px"))),
        card(card_header("Rituximab"),
             card_body(plotlyOutput("pk_ritu", height = "260px"))),
        card(card_header("Doxycycline"),
             card_body(plotlyOutput("pk_doxy", height = "260px")))
      )
    )
  ),

  ## ── TAB 3: PD Biomarkers ────────────────────────────────────────────────────
  nav_panel(
    "PD Biomarkers",
    icon = icon("flask"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("Treatment Settings"),
        numericInput("pred_pd", "Prednisolone (mg/kg/d)", 0.5, 0, 1.5, 0.1),
        checkboxInput("dup_pd", "Add Dupilumab", FALSE),
        checkboxInput("oma_pd", "Add Omalizumab", FALSE),
        checkboxInput("ritu_pd", "Add Rituximab", FALSE),
        checkboxInput("doxy_pd", "Add Doxycycline", FALSE),
        sliderInput("bw_pd", "Body Weight (kg)", 40, 120, 65, 5),
        sliderInput("wk_pd", "Duration (weeks)", 12, 104, 52, 4),
        actionButton("run_pd", "Run PD Simulation", class = "btn-success btn-block"),
        hr(),
        checkboxGroupInput("biomarkers_show",
                           "Biomarkers to display:",
                           c("Anti-BP180 IgG" = "IGG_BP180",
                             "Anti-BP180 IgE" = "IGE_BP180",
                             "Blood Eosinophils" = "EOS_BLOOD",
                             "Skin Eosinophils" = "EOS_SKIN",
                             "Mast Cells (Skin)" = "MAST_ACT",
                             "C5a Complement" = "C5A",
                             "DEJ Integrity" = "DEJ",
                             "Th2 Cells" = "TH2",
                             "LLPC (Plasma Cells)" = "LLPC"),
                           selected = c("IGG_BP180","IGE_BP180","EOS_BLOOD","DEJ"))
      ),
      layout_columns(
        col_widths = c(8, 4),
        card(card_header("Biomarker Trajectories"),
             card_body(plotlyOutput("pd_biomarkers", height = "400px"))),
        card(card_header("B Cell Dynamics"),
             card_body(plotlyOutput("pd_bcell", height = "400px")))
      ),
      card(
        card_header("Biomarker Reference Ranges"),
        card_body(
          DT::dataTableOutput("pd_ref_table")
        )
      )
    )
  ),

  ## ── TAB 4: Clinical Endpoints ────────────────────────────────────────────────
  nav_panel(
    "Clinical Endpoints",
    icon = icon("chart-line"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("Treatment Inputs"),
        numericInput("pred_ep", "Prednisolone (mg/kg/d)", 0.5, 0, 1.5, 0.1),
        checkboxInput("dup_ep", "Add Dupilumab", FALSE),
        checkboxInput("oma_ep", "Add Omalizumab", FALSE),
        checkboxInput("ritu_ep", "Add Rituximab", FALSE),
        checkboxInput("doxy_ep", "Add Doxycycline", FALSE),
        sliderInput("bw_ep", "Body Weight (kg)", 40, 120, 65, 5),
        sliderInput("wk_ep", "Duration (weeks)", 12, 104, 52, 4),
        actionButton("run_ep", "Run Endpoint Simulation", class = "btn-warning btn-block"),
        hr(),
        h6("Remission Criteria"),
        tags$small(
          tags$b("Complete Remission:"), " BPDAI = 0 for ≥8 weeks off therapy",
          tags$br(), tags$b("Partial Remission:"), " ≤3 new blisters/day",
          tags$br(), tags$b("Relapse:"), " ≥3 new blisters/day after CR"
        )
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("BPDAI Activity Score"),
             card_body(plotlyOutput("ep_bpdai", height = "280px"))),
        card(card_header("Itch NRS Score"),
             card_body(plotlyOutput("ep_itch", height = "280px")))
      ),
      layout_columns(
        col_widths = c(4, 4, 4),
        card(card_header("BSA Affected (%)"),
             card_body(plotlyOutput("ep_bsa", height = "250px"))),
        card(card_header("New Blisters/Day"),
             card_body(plotlyOutput("ep_blisters", height = "250px"))),
        card(card_header("DEJ Integrity"),
             card_body(plotlyOutput("ep_dej", height = "250px")))
      ),
      card(
        card_header("Weekly Summary Table"),
        card_body(DT::dataTableOutput("ep_table"))
      )
    )
  ),

  ## ── TAB 5: Scenario Comparison ──────────────────────────────────────────────
  nav_panel(
    "Scenario Comparison",
    icon = icon("balance-scale"),
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        h5("Select Scenarios to Compare"),
        checkboxInput("sc1", "Untreated (no therapy)", TRUE),
        checkboxInput("sc2", "High-dose Prednisolone (0.75 mg/kg taper)", TRUE),
        checkboxInput("sc3", "Low-dose Pred + Doxycycline", TRUE),
        checkboxInput("sc4", "Dupilumab + Low-dose Pred", TRUE),
        checkboxInput("sc5", "Omalizumab + Low-dose Pred", TRUE),
        checkboxInput("sc6", "Rituximab + Short Pred course", TRUE),
        hr(),
        sliderInput("bw_sc", "Body Weight (kg)", 40, 120, 65, 5),
        sliderInput("wk_sc", "Duration (weeks)", 12, 104, 52, 4),
        actionButton("run_sc", "Run All Scenarios", class = "btn-primary btn-block")
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("BPDAI Comparison"),
             card_body(plotlyOutput("sc_bpdai", height = "320px"))),
        card(card_header("Itch NRS Comparison"),
             card_body(plotlyOutput("sc_itch", height = "320px")))
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Anti-BP180 IgG Titer Comparison"),
             card_body(plotlyOutput("sc_igg", height = "300px"))),
        card(card_header("Blood Eosinophils Comparison"),
             card_body(plotlyOutput("sc_eos", height = "300px")))
      ),
      card(
        card_header("Efficacy Summary Table — All Scenarios"),
        card_body(DT::dataTableOutput("sc_table"))
      )
    )
  ),

  ## ── TAB 6: Biomarker Explorer ───────────────────────────────────────────────
  nav_panel(
    "Biomarker Explorer",
    icon = icon("microscope"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280,
        h5("Biomarker Analysis"),
        numericInput("pred_bm", "Prednisolone (mg/kg/d)", 0.5, 0, 1.5, 0.1),
        checkboxInput("dup_bm", "Add Dupilumab", FALSE),
        checkboxInput("oma_bm", "Add Omalizumab", FALSE),
        checkboxInput("ritu_bm", "Add Rituximab", FALSE),
        checkboxInput("doxy_bm", "Add Doxycycline", FALSE),
        sliderInput("bw_bm", "Body Weight (kg)", 40, 120, 65, 5),
        sliderInput("wk_bm", "Duration (weeks)", 12, 104, 52, 4),
        actionButton("run_bm", "Run Simulation", class = "btn-info btn-block"),
        hr(),
        h6("Biomarker Correlation"),
        selectInput("bm_x", "X-axis Biomarker",
                    c("Anti-BP180 IgG"="IGG_BP180",
                      "Blood Eosinophils"="EOS_BLOOD",
                      "Skin Eosinophils"="EOS_SKIN",
                      "C5a"="C5A",
                      "Mast Cells"="MAST_ACT",
                      "Th2"="TH2"),
                    selected = "IGG_BP180"),
        selectInput("bm_y", "Y-axis Biomarker",
                    c("BPDAI Score"="BPDAI_act",
                      "Itch NRS"="Itch",
                      "BSA Affected"="BSA",
                      "DEJ Integrity"="DEJ",
                      "Anti-BP180 IgE"="IGE_BP180"),
                    selected = "BPDAI_act")
      ),
      layout_columns(
        col_widths = c(7, 5),
        card(card_header("All Biomarkers — Heatmap (% change from baseline)"),
             card_body(plotlyOutput("bm_heatmap", height = "420px"))),
        card(card_header("Biomarker Correlation Plot"),
             card_body(plotlyOutput("bm_corr", height = "420px")))
      ),
      card(
        card_header("Mechanistic Pathway Activation Status"),
        card_body(
          layout_columns(
            col_widths = c(4, 4, 4),
            card(
              card_header("Th2 / B Cell Axis"),
              card_body(
                div(style="font-size:13px",
                    p("• IL-4 / IL-13 → IgE class switch"),
                    p("• IL-5 → Eosinophil proliferation"),
                    p("• Tfh cells → GC B cell help"),
                    p("• BAFF/APRIL → Plasma cell survival"),
                    p("➡ Targets: dupilumab (IL-4Rα), mepolizumab (IL-5)")
                )
              )
            ),
            card(
              card_header("IgE / Mast Cell Axis"),
              card_body(
                div(style="font-size:13px",
                    p("• Anti-BP180 IgE → FcεRI cross-linking"),
                    p("• Mast cell degranulation → histamine, tryptase"),
                    p("• Histamine → pruritus (H1R)"),
                    p("• IL-31 → JAK1/STAT3 itch"),
                    p("➡ Targets: omalizumab (anti-IgE), nemolizumab (IL-31R)")
                )
              )
            ),
            card(
              card_header("Complement / Eosinophil Axis"),
              card_body(
                div(style="font-size:13px",
                    p("• Anti-BP180 IgG → C1q → C5a"),
                    p("• C5a → eosinophil/neutrophil recruitment"),
                    p("• MBP, ECP → DEJ proteolysis"),
                    p("• MMP-9 → collagen XVII cleavage"),
                    p("➡ Targets: avacopan (C5aR1), doxycycline (MMP inh.)")
                )
              )
            )
          )
        )
      )
    )
  )
)

## ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## ─ Reactive Simulation Data ─────────────────────────────────────────────────

  # PK Tab
  pk_data <- eventReactive(input$run_pk, {
    run_bp_sim(
      pred_dose = input$pred_dose_pk,
      bw = input$bw_pk,
      use_dup = input$use_dup_pk,
      use_oma = input$use_oma_pk,
      use_ritu = input$use_ritu_pk,
      use_doxy = input$use_doxy_pk,
      weeks = input$weeks_pk
    )
  }, ignoreNULL = FALSE)

  # PD Tab
  pd_data <- eventReactive(input$run_pd, {
    run_bp_sim(input$pred_pd, input$bw_pd, input$dup_pd, input$oma_pd,
               input$ritu_pd, input$doxy_pd, input$wk_pd)
  }, ignoreNULL = FALSE)

  # Endpoints Tab
  ep_data <- eventReactive(input$run_ep, {
    run_bp_sim(input$pred_ep, input$bw_ep, input$dup_ep, input$oma_ep,
               input$ritu_ep, input$doxy_ep, input$wk_ep)
  }, ignoreNULL = FALSE)

  # Biomarker Explorer
  bm_data <- eventReactive(input$run_bm, {
    run_bp_sim(input$pred_bm, input$bw_bm, input$dup_bm, input$oma_bm,
               input$ritu_bm, input$doxy_bm, input$wk_bm)
  }, ignoreNULL = FALSE)

  # Scenario Comparison
  sc_data <- eventReactive(input$run_sc, {
    scen_list <- list()
    bw <- input$bw_sc; wks <- input$wk_sc
    if (input$sc1) scen_list[["Untreated"]] <- run_bp_sim(0, bw, F, F, F, F, wks)
    if (input$sc2) scen_list[["High-dose Pred\n(0.75 mg/kg)"]] <- run_bp_sim(0.75, bw, F, F, F, F, wks)
    if (input$sc3) scen_list[["Low-dose Pred\n+ Doxy"]] <- run_bp_sim(0.30, bw, F, F, F, T, wks)
    if (input$sc4) scen_list[["Dupilumab\n+ Low Pred"]] <- run_bp_sim(0.25, bw, T, F, F, F, wks)
    if (input$sc5) scen_list[["Omalizumab\n+ Low Pred"]] <- run_bp_sim(0.30, bw, F, T, F, F, wks)
    if (input$sc6) scen_list[["Rituximab\n+ Short Pred"]] <- run_bp_sim(0.50, bw, F, F, T, F, wks)
    lapply(names(scen_list), function(nm) {
      scen_list[[nm]] %>% mutate(Scenario = nm)
    }) %>% bind_rows()
  }, ignoreNULL = FALSE)

  ## ─ PK Outputs ───────────────────────────────────────────────────────────────
  make_pk_plot <- function(df, y_col, y_lab, color) {
    plot_ly(df, x = ~time_weeks, y = as.formula(paste0("~", y_col)),
            type = "scatter", mode = "lines",
            line = list(color = color, width = 2)) %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = y_lab),
             hovermode = "x unified")
  }

  output$pk_pred <- renderPlotly({
    make_pk_plot(pk_data(), "Cpred_ngml", "Prednisolone (ng/mL)", "#FF7043")
  })
  output$pk_dup <- renderPlotly({
    make_pk_plot(pk_data(), "Cdup_mgL", "Dupilumab (mg/L)", "#42A5F5")
  })
  output$pk_oma <- renderPlotly({
    make_pk_plot(pk_data(), "Coma_mgL", "Omalizumab (mg/L)", "#26C6DA")
  })
  output$pk_ritu <- renderPlotly({
    make_pk_plot(pk_data(), "Critu_mgL", "Rituximab (mg/L)", "#AB47BC")
  })
  output$pk_doxy <- renderPlotly({
    make_pk_plot(pk_data(), "Cdoxy_mcg", "Doxycycline (µg/mL)", "#8D6E63")
  })

  ## ─ PD Outputs ───────────────────────────────────────────────────────────────
  output$pd_biomarkers <- renderPlotly({
    df <- pd_data()
    bm <- input$biomarkers_show
    if (length(bm) == 0) return(plotly_empty())
    df_long <- df %>% select(time_weeks, all_of(bm)) %>%
      pivot_longer(-time_weeks, names_to = "Biomarker", values_to = "Value")
    plot_ly(df_long, x = ~time_weeks, y = ~Value, color = ~Biomarker,
            type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = "Value (normalised)"),
             hovermode = "x unified",
             legend = list(orientation = "h"))
  })

  output$pd_bcell <- renderPlotly({
    df <- pd_data() %>%
      select(time_weeks, B_NAIVE, B_ACT, LLPC) %>%
      pivot_longer(-time_weeks, names_to = "Cell", values_to = "Count")
    plot_ly(df, x = ~time_weeks, y = ~Count, color = ~Cell,
            type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Weeks"),
             yaxis = list(title = "Count (AU)"),
             hovermode = "x unified")
  })

  output$pd_ref_table <- DT::renderDataTable({
    data.frame(
      Biomarker = c("Anti-BP180 IgG", "Anti-BP180 IgE", "Blood Eosinophils",
                    "Skin Eosinophils", "C5a", "DEJ Integrity", "BPDAI"),
      `Normal Range` = c("<9 AU/mL", "Undetectable", "<0.4×10⁹/L",
                          "Absent", "Undetectable", "1.0 (intact)",
                          "0 (normal)"),
      `Active BP` = c(">20 AU/mL", "Elevated", ">0.5×10⁹/L",
                       "Moderate-high", "Elevated", "<0.5 (disrupted)",
                       ">20")
    )
  }, options = list(dom = "t", pageLength = 10))

  ## ─ Endpoint Outputs ──────────────────────────────────────────────────────────
  output$ep_bpdai <- renderPlotly({
    plot_ly(ep_data(), x = ~time_weeks, y = ~BPDAI_act,
            type = "scatter", mode = "lines",
            line = list(color = "#E53935", width = 2.5)) %>%
      layout(xaxis = list(title = "Weeks"),
             yaxis = list(title = "BPDAI Activity (0–90)", range = c(0, 90)))
  })
  output$ep_itch <- renderPlotly({
    plot_ly(ep_data(), x = ~time_weeks, y = ~Itch,
            type = "scatter", mode = "lines",
            line = list(color = "#FF7043", width = 2.5)) %>%
      layout(xaxis = list(title = "Weeks"),
             yaxis = list(title = "Itch NRS (0–10)", range = c(0, 10)))
  })
  output$ep_bsa <- renderPlotly({
    plot_ly(ep_data(), x = ~time_weeks, y = ~BSA,
            type = "scatter", mode = "lines",
            line = list(color = "#FFA726", width = 2)) %>%
      layout(xaxis = list(title = "Weeks"), yaxis = list(title = "BSA (%)"))
  })
  output$ep_blisters <- renderPlotly({
    plot_ly(ep_data(), x = ~time_weeks, y = ~new_blisters,
            type = "scatter", mode = "lines",
            line = list(color = "#AB47BC", width = 2)) %>%
      layout(xaxis = list(title = "Weeks"), yaxis = list(title = "Blisters/day"))
  })
  output$ep_dej <- renderPlotly({
    plot_ly(ep_data(), x = ~time_weeks, y = ~DEJ,
            type = "scatter", mode = "lines",
            line = list(color = "#66BB6A", width = 2)) %>%
      layout(xaxis = list(title = "Weeks"),
             yaxis = list(title = "DEJ Integrity (0–1)", range = c(0, 1)))
  })
  output$ep_table <- DT::renderDataTable({
    ep_data() %>%
      filter(time_weeks %% 4 < 0.1) %>%
      select(time_weeks, BPDAI_act, Itch, BSA, new_blisters, IGG_BP180, DEJ) %>%
      mutate(across(where(is.numeric), ~round(.x, 2))) %>%
      rename("Week"="time_weeks", "BPDAI"="BPDAI_act",
             "Itch NRS"="Itch", "BSA%"="BSA",
             "Blisters/d"="new_blisters",
             "Anti-BP180 IgG"="IGG_BP180",
             "DEJ Integrity"="DEJ")
  }, options = list(pageLength = 13, dom = "tip"))

  ## ─ Scenario Comparison Outputs ───────────────────────────────────────────────
  sc_plot <- function(y_col, y_lab) {
    df <- sc_data()
    if (nrow(df) == 0) return(plotly_empty())
    plot_ly(df, x = ~time_weeks, y = as.formula(paste0("~", y_col)),
            color = ~Scenario, type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Weeks"),
             yaxis = list(title = y_lab),
             hovermode = "x unified",
             legend = list(orientation = "h", y = -0.3))
  }

  output$sc_bpdai <- renderPlotly({ sc_plot("BPDAI_act", "BPDAI Activity (0–90)") })
  output$sc_itch  <- renderPlotly({ sc_plot("Itch",       "Itch NRS (0–10)") })
  output$sc_igg   <- renderPlotly({ sc_plot("IGG_BP180",  "Anti-BP180 IgG (AU)") })
  output$sc_eos   <- renderPlotly({ sc_plot("EOS_BLOOD",  "Blood Eosinophils (×10⁹/L)") })

  output$sc_table <- DT::renderDataTable({
    req(nrow(sc_data()) > 0)
    sc_data() %>%
      group_by(Scenario) %>%
      summarise(
        `BPDAI Wk4`   = round(BPDAI_act[which.min(abs(time_weeks-4))],  1),
        `BPDAI Wk12`  = round(BPDAI_act[which.min(abs(time_weeks-12))], 1),
        `BPDAI Wk52`  = round(BPDAI_act[which.min(abs(time_weeks-52))], 1),
        `Itch Wk12`   = round(Itch[which.min(abs(time_weeks-12))],      1),
        `Anti-BP180 Wk52` = round(IGG_BP180[which.min(abs(time_weeks-52))]*100, 0),
        `Remission%`  = round(mean(in_remission[time_weeks > 12])*100,  0),
        .groups = "drop"
      )
  }, options = list(dom = "t"))

  ## ─ Biomarker Explorer Outputs ────────────────────────────────────────────────
  output$bm_heatmap <- renderPlotly({
    df <- bm_data()
    cols <- c("IGG_BP180","IGE_BP180","EOS_BLOOD","EOS_SKIN",
              "MAST_ACT","C5A","DEJ","TH2","LLPC","BPDAI_act","Itch","BSA")
    # Baseline values at time 0
    baseline <- df %>% filter(time_weeks == 0) %>% select(all_of(cols)) %>%
      summarise(across(everything(), mean)) %>% unlist()
    # Weekly snapshots at 4, 12, 26, 52 wks
    wk_snaps <- c(4, 12, 26, 52)
    pct_change <- do.call(rbind, lapply(wk_snaps, function(w) {
      row <- df %>% filter(abs(time_weeks - w) == min(abs(time_weeks - w))) %>%
        head(1) %>% select(all_of(cols)) %>% unlist()
      (row - baseline) / (abs(baseline) + 1e-6) * 100
    }))
    rownames(pct_change) <- paste0("Wk", wk_snaps)

    plot_ly(
      x = cols, y = rownames(pct_change),
      z = pct_change,
      type = "heatmap",
      colorscale = list(c(0,"#2196F3"), c(0.5,"#FFFFFF"), c(1,"#F44336")),
      zmid = 0,
      colorbar = list(title = "% Change\nfrom BL")
    ) %>%
      layout(
        xaxis = list(title = "", tickangle = 45),
        yaxis = list(title = "Time Point")
      )
  })

  output$bm_corr <- renderPlotly({
    df <- bm_data()
    plot_ly(df, x = as.formula(paste0("~", input$bm_x)),
            y = as.formula(paste0("~", input$bm_y)),
            type = "scatter", mode = "markers",
            marker = list(color = ~time_weeks, colorscale = "Viridis",
                          size = 5, showscale = TRUE,
                          colorbar = list(title = "Weeks"))) %>%
      layout(xaxis = list(title = input$bm_x),
             yaxis = list(title = input$bm_y))
  })
}

## ── Launch ────────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
