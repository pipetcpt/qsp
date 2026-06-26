################################################################################
# Giant Cell Arteritis (GCA) — Interactive QSP Shiny Dashboard
# 거대세포 동맥염 QSP 인터랙티브 대시보드
#
# Author : Claude Code Routine (CCR) | 2026-06-17
# Package: shiny, mrgsolve, ggplot2, plotly, DT, dplyr
#
# Tabs (6+):
#  1. Patient Profile       — demographics, GCA phenotype, baseline parameters
#  2. Pharmacokinetics (PK) — drug concentration-time profiles
#  3. PD Biomarkers         — CRP, ESR, IL-6, sIL-6R, VEGF dynamics
#  4. Disease Activity      — GCA-DAS, remission, relapse probability
#  5. Treatment Comparison  — head-to-head scenario comparison table & plots
#  6. Safety Monitoring     — cumulative GC dose, BMD loss, risk assessment
#  7. Mechanistic Map       — embedded QSP diagram viewer
################################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(patchwork)
library(tidyr)

# ── Embed the mrgsolve model code ─────────────────────────────────────────────
gca_model_code <- '
$PROB GCA QSP Shiny Model — Giant Cell Arteritis

$PARAM
ka_pred=1.35 F_pred=0.82 kconv_pred=0.925 CL_predsl=1.67 Vd_predsl=32.0
k12_predsl=0.12 k21_predsl=0.06
TCZ_CL=0.0108 TCZ_V1=3.72 TCZ_V2=3.44 TCZ_Q=0.315
TCZ_SC_F=0.80 TCZ_SC_ka=0.0031
TCZ_TMDD_Kin=0.035 TCZ_TMDD_Kout=0.012 TCZ_TMDD_kon=0.094 TCZ_TMDD_koff=0.00014
IL6_kin_base=0.42 IL6_kout=0.693 IL6_DA_amp=8.0
CRP_kin_max=2.4 CRP_EC50_IL6=12.5 CRP_kout=0.0365
ESR_kin=3.0 ESR_kout=0.02 ESR_IL6_slope=0.8
DA_baseline=1.0 DA_kout=0.004
DA_GC_Emax=0.90 DA_GC_EC50=8.0
DA_TCZ_Emax=0.80 DA_TCZ_EC50=0.15
Mac_kin=0.01 Mac_kout=0.008 Mac_GC_Emax=0.75 Mac_GC_EC50=5.0
Th17_kin=0.005 Th17_kout=0.004 Th17_TCZ_Emax=0.60 Th17_TCZ_EC50=0.25
VEGF_kin=0.5 VEGF_kout=0.046 VEGF_Th17_slope=40.0
BMD_kout=0.00009 BMD_GC_slope=0.000005
TCZ_route=0 use_pred=1
bw=70

$INIT
Pred_GI=0 Pred_plasma=0 Prednisolone=0 Pred_tissue=0
TCZ_depot=0 TCZ_central=0 TCZ_periph=0
sIL6R_free=3.0 TCZ_sIL6R=0
IL6_free=95.0 IL6_bound=0.5
CRP=80.0 ESR_state=85.0
DA=1.0 MacActiv=1.0 Th17_rel=1.0 VEGF_state=320.0
BMD_rel=1.0 CumGC=0

$ODE
double GC_effect = Prednisolone / (Prednisolone + DA_GC_EC50);
double GC_effect_mac = (Mac_GC_Emax * Prednisolone) / (Prednisolone + Mac_GC_EC50);
double total_sIL6R = sIL6R_free + TCZ_sIL6R;
double IL6R_block_frac = (total_sIL6R > 0.01) ? TCZ_sIL6R / total_sIL6R : 0.0;
double TCZ_DA_effect = (DA_TCZ_Emax * IL6R_block_frac) / (IL6R_block_frac + DA_TCZ_EC50);
double DA_suppress = 1.0 - (1.0 - DA_GC_Emax * GC_effect) * (1.0 - TCZ_DA_effect);
double DA_suppress_capped = (DA_suppress > 0.98) ? 0.98 : DA_suppress;

dxdt_Pred_GI = -ka_pred * Pred_GI;
double kconv_rate = 0.5 * ka_pred;
dxdt_Pred_plasma = ka_pred * Pred_GI - kconv_rate * Pred_plasma - (CL_predsl/Vd_predsl)*0.1*Pred_plasma;
dxdt_Prednisolone = kconv_pred * kconv_rate * Pred_plasma - (CL_predsl/Vd_predsl)*Prednisolone
                    - k12_predsl*Prednisolone + k21_predsl*Pred_tissue;
dxdt_Pred_tissue = k12_predsl * Prednisolone - k21_predsl * Pred_tissue;
dxdt_TCZ_depot = -TCZ_SC_ka * TCZ_depot;
double TCZ_SC_influx = (TCZ_route == 2) ? TCZ_SC_F * TCZ_SC_ka * TCZ_depot : 0.0;
double TCZ_bind_rate = TCZ_TMDD_kon * (TCZ_central/TCZ_V1) * sIL6R_free;
double TCZ_unbind_rate = TCZ_TMDD_koff * TCZ_sIL6R;
dxdt_TCZ_central = TCZ_SC_influx - TCZ_CL*(TCZ_central/TCZ_V1)
                   - TCZ_Q*(TCZ_central/TCZ_V1 - TCZ_periph/TCZ_V2)
                   - TCZ_bind_rate*TCZ_V1 + TCZ_unbind_rate*TCZ_V1;
dxdt_TCZ_periph = TCZ_Q*(TCZ_central/TCZ_V1 - TCZ_periph/TCZ_V2);
dxdt_sIL6R_free = TCZ_TMDD_Kin - TCZ_TMDD_Kout*sIL6R_free - TCZ_bind_rate;
dxdt_TCZ_sIL6R = TCZ_bind_rate - TCZ_unbind_rate - TCZ_TMDD_Kout*TCZ_sIL6R;
double IL6_prod = IL6_kin_base * (1.0 + IL6_DA_amp * DA * MacActiv);
double IL6_bind = 0.05 * IL6_free * (sIL6R_free/(sIL6R_free+1.0));
dxdt_IL6_free = IL6_prod - IL6_kout*IL6_free - IL6_bind;
dxdt_IL6_bound = IL6_bind - 0.3*IL6_bound;
double CRP_prod = CRP_kin_max * (IL6_free/(IL6_free + CRP_EC50_IL6));
dxdt_CRP = CRP_prod - CRP_kout*CRP;
double ESR_il6_effect = 1.0 + ESR_IL6_slope*(IL6_free/10.0);
dxdt_ESR_state = ESR_kin * ESR_il6_effect - ESR_kout*ESR_state;
double DA_inflow = DA_kout * DA_baseline * (1.0-DA_suppress_capped) * (MacActiv+Th17_rel)/2.0;
dxdt_DA = DA_inflow - DA_kout*DA;
if(DA < 0.0) dxdt_DA = -DA*10.0;
double Mac_inflow = Mac_kin * DA;
double Mac_suppress_total = GC_effect_mac + (1.0-GC_effect_mac)*IL6R_block_frac*0.5;
dxdt_MacActiv = Mac_inflow*(1.0-Mac_suppress_total) - Mac_kout*MacActiv;
if(MacActiv < 0.0) dxdt_MacActiv = -MacActiv*10.0;
double Th17_TCZ_suppress = (Th17_TCZ_Emax*IL6R_block_frac)/(IL6R_block_frac+Th17_TCZ_EC50);
double Th17_GC_suppress = 0.4*GC_effect;
dxdt_Th17_rel = Th17_kin*DA - Th17_kout*Th17_rel*(1.0+Th17_TCZ_suppress+Th17_GC_suppress);
if(Th17_rel < 0.0) dxdt_Th17_rel = -Th17_rel*10.0;
dxdt_VEGF_state = VEGF_kin + VEGF_Th17_slope*Th17_rel*0.01 - VEGF_kout*VEGF_state;
dxdt_BMD_rel = -BMD_GC_slope * Prednisolone;
if(BMD_rel < 0.60) dxdt_BMD_rel = 0.0;
dxdt_CumGC = Prednisolone;

$TABLE
double CRP_obs  = CRP;
double ESR_obs  = ESR_state;
double IL6_obs  = IL6_free;
double VEGF_obs = VEGF_state;
double DA_obs   = DA * 100.0;
double BMD_obs  = BMD_rel * 100.0;
double CumGC_g  = CumGC / 24000.0;
double TCZ_Cp   = TCZ_central / TCZ_V1;
double Pred_Cp  = Prednisolone;
double IL6R_bl  = (sIL6R_free + TCZ_sIL6R > 0.01) ? TCZ_sIL6R/(sIL6R_free+TCZ_sIL6R)*100.0 : 0.0;
double sIL6R_o  = sIL6R_free;
double remiss   = (CRP_obs < 5.0 && ESR_obs < 20.0 && DA_obs < 15.0) ? 1.0 : 0.0;

$CAPTURE CRP_obs ESR_obs IL6_obs VEGF_obs DA_obs BMD_obs CumGC_g TCZ_Cp Pred_Cp IL6R_bl sIL6R_o remiss
'

# Compile model (cached globally)
gca_mod <- mcode("gca_shiny_v1", gca_model_code, quiet = TRUE)

# ── Taper helpers ─────────────────────────────────────────────────────────────
make_pred_events <- function(start_dose, taper_speed = "slow") {
  if (taper_speed == "slow") {
    sched <- data.frame(wk = c(0,2,4,8,12,16,20,26), mg = c(start_dose,40,30,25,20,15,10,0))
  } else if (taper_speed == "fast") {
    sched <- data.frame(wk = c(0,2,4,6,8,12,18,26), mg = c(start_dose,40,30,25,20,10,5,0))
  } else {
    sched <- data.frame(wk = c(0,4,8,12,26,52), mg = c(start_dose,30,20,15,10,5))
  }
  ev_list <- list()
  for (i in seq_len(nrow(sched) - 1)) {
    hrs <- seq(sched$wk[i]*168, (sched$wk[i+1]*168-24), by=24)
    ev_list[[i]] <- data.frame(time=hrs, amt=sched$mg[i]*0.82, cmt=1, evid=1, rate=0)
  }
  do.call(rbind, ev_list)
}

make_tcz_iv_events <- function(n_doses=13, bw=70) {
  data.frame(
    time = seq(0,(n_doses-1)*4*168, by=4*168),
    amt  = 8*bw, cmt=6, evid=1, rate=-2
  )
}

make_tcz_sc_events <- function(freq="qw", n_doses=52) {
  interval <- if(freq=="qw") 168 else 336
  data.frame(
    time = seq(0,(n_doses-1)*interval, by=interval),
    amt=162, cmt=5, evid=1, rate=0
  )
}

run_scenario <- function(mod, pred_ev=NULL, tcz_ev=NULL, tcz_route=0,
                          sim_wks=52, label="Scenario") {
  tgrid <- seq(0, sim_wks*168, by=12)
  ev    <- bind_rows(pred_ev, tcz_ev)
  if (nrow(ev) == 0) return(NULL)
  mod %>%
    param(TCZ_route = tcz_route) %>%
    data_set(ev) %>%
    mrgsim(tgrid = tgrid, obsonly = TRUE) %>%
    as.data.frame() %>%
    mutate(time_wk = time/168, scenario = label)
}

# ── Color palette ─────────────────────────────────────────────────────────────
pal_scenarios <- c(
  "GC Mono (Slow Taper)"     = "#E67E22",
  "GC Mono (Fast Taper)"     = "#E74C3C",
  "TCZ IV q4w + GC"          = "#2980B9",
  "TCZ SC qw + GC"           = "#27AE60",
  "TCZ SC q2w + GC"          = "#8E44AD",
  "Custom Regimen"            = "#1ABC9C"
)

theme_gca <- theme_bw(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face="bold", size=14))

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = HTML("GCA QSP Dashboard"),
    titleWidth = 260
  ),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      id = "tabs",
      menuItem("Patient Profile",    tabName = "profile",    icon = icon("user-md")),
      menuItem("Pharmacokinetics",   tabName = "pk",         icon = icon("flask")),
      menuItem("PD Biomarkers",      tabName = "pd",         icon = icon("chart-line")),
      menuItem("Disease Activity",   tabName = "disease",    icon = icon("heartbeat")),
      menuItem("Treatment Compare",  tabName = "compare",    icon = icon("balance-scale")),
      menuItem("Safety Monitoring",  tabName = "safety",     icon = icon("shield-alt")),
      menuItem("Mechanistic Map",    tabName = "map",        icon = icon("project-diagram"))
    ),

    hr(),
    h5("  Global Settings", style="color:#BDC3C7; padding-left:12px;"),

    sliderInput("sim_weeks", "Simulation Duration (weeks):",
                min=12, max=104, value=52, step=4),
    sliderInput("bw", "Body Weight (kg):", min=40, max=120, value=70, step=5),

    hr(),
    h5("  GC Regimen", style="color:#BDC3C7; padding-left:12px;"),
    sliderInput("pred_start", "Initial Prednisone (mg/day):",
                min=20, max=80, value=60, step=5),
    selectInput("pred_taper", "Taper Protocol:",
                choices=c("Slow (GiACTA-like)"="slow",
                          "Fast (Placebo-like)"="fast",
                          "Prolonged"="prolonged")),

    hr(),
    h5("  Biologic Therapy", style="color:#BDC3C7; padding-left:12px;"),
    selectInput("tcz_regimen", "Tocilizumab Regimen:",
                choices=c("None"="none",
                          "IV 8 mg/kg q4w"="iv",
                          "SC 162 mg qw"="sc_qw",
                          "SC 162 mg q2w"="sc_q2w")),
    numericInput("tcz_start_wk", "TCZ Start (week):", value=0, min=0, max=52)
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background-color: #F8F9FA; }
        .box { border-radius: 8px; }
        .info-box { border-radius: 8px; }
        .value-box { border-radius: 8px; }
      "))
    ),

    tabItems(
      # ── TAB 1: Patient Profile ──────────────────────────────────────────────
      tabItem(
        tabName = "profile",
        fluidRow(
          box(
            title = "Patient Demographics & GCA Phenotype",
            width = 5, status = "primary", solidHeader = TRUE,
            selectInput("sex", "Sex:", choices=c("Female"=1,"Male"=0), selected=1),
            sliderInput("age", "Age (years):", min=50, max=90, value=72, step=1),
            sliderInput("duration_sxs", "Symptom Duration (weeks):", min=1, max=52, value=6),
            selectInput("phenotype", "GCA Phenotype:",
                        choices=c("Cranial GCA (classic)"="cranial",
                                  "Large Vessel GCA (LV-GCA)"="lv_gca",
                                  "PMR + GCA overlap"="pmr_gca",
                                  "Occult GCA (biopsy-proven)"="occult")),
            checkboxInput("visual_sx", "Visual Symptoms at Presentation", value=FALSE),
            checkboxInput("jaw_claud", "Jaw Claudication", value=TRUE),
            checkboxInput("pmr_sx", "Polymyalgia Rheumatica Symptoms", value=TRUE),
            sliderInput("baseline_crp", "Baseline CRP (mg/L):", min=20, max=200, value=80),
            sliderInput("baseline_esr", "Baseline ESR (mm/hr):", min=30, max=120, value=85)
          ),
          box(
            title = "ACR/EULAR 2022 Classification Criteria",
            width = 4, status = "warning", solidHeader = TRUE,
            p("Scoring ≥6 = GCA classification", style="font-weight:bold;"),
            tableOutput("acr_score_tbl"),
            hr(),
            valueBoxOutput("acr_total_box", width=12)
          ),
          box(
            title = "Baseline Inflammatory Markers",
            width = 3, status = "danger", solidHeader = TRUE,
            valueBoxOutput("vb_crp", width=12),
            valueBoxOutput("vb_esr", width=12),
            valueBoxOutput("vb_il6", width=12),
            p(style="font-size:11px; color:#7F8C8D; margin-top:10px;",
              "IL-6 estimated from CRP/ESR correlation (Weyand 2003).
               Normal IL-6: <7 pg/mL. Active GCA: 50-300 pg/mL.")
          )
        ),
        fluidRow(
          box(
            title = "Key Clinical Features of Giant Cell Arteritis",
            width = 12, status = "info", solidHeader = TRUE, collapsible=TRUE,
            column(4,
              h4("Cranial Manifestations"),
              tags$ul(
                tags$li("New temporal headache (87%)"),
                tags$li("Jaw claudication (50%)"),
                tags$li("Scalp tenderness (47%)"),
                tags$li("Temporal artery abnormality on exam (53%)"),
                tags$li("Visual loss (15-20% if untreated)")
              )
            ),
            column(4,
              h4("Systemic & Rheumatic"),
              tags$ul(
                tags$li("Fever, weight loss, night sweats"),
                tags$li("PMR in 40-50% of GCA patients"),
                tags$li("Fatigue and malaise"),
                tags$li("Normochromic, normocytic anemia"),
                tags$li("Elevated CRP/ESR (near universal)")
              )
            ),
            column(4,
              h4("Large Vessel (LV-GCA)"),
              tags$ul(
                tags$li("Subclavian/axillary stenosis (>15%)"),
                tags$li("Aortic involvement (18-58%)"),
                tags$li("Arm claudication"),
                tags$li("Blood pressure asymmetry"),
                tags$li("Aortic aneurysm (long-term risk)")
              )
            )
          )
        )
      ),

      # ── TAB 2: Pharmacokinetics ─────────────────────────────────────────────
      tabItem(
        tabName = "pk",
        fluidRow(
          box(
            title = "Drug Concentration — Time Profiles",
            width = 12, status = "primary", solidHeader = TRUE,
            plotlyOutput("pk_plot", height="480px")
          )
        ),
        fluidRow(
          box(
            title = "PK Parameter Summary",
            width = 6, status = "info", solidHeader = TRUE,
            tableOutput("pk_param_tbl")
          ),
          box(
            title = "Dosing Event Log",
            width = 6, status = "info", solidHeader = TRUE,
            p("Prednisone oral taper schedule"),
            tableOutput("dose_schedule_tbl")
          )
        )
      ),

      # ── TAB 3: PD Biomarkers ────────────────────────────────────────────────
      tabItem(
        tabName = "pd",
        fluidRow(
          box(
            title = "CRP — C-Reactive Protein",
            width = 6, status = "danger", solidHeader = TRUE,
            plotlyOutput("crp_plot", height="300px")
          ),
          box(
            title = "ESR — Erythrocyte Sedimentation Rate",
            width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("esr_plot", height="300px")
          )
        ),
        fluidRow(
          box(
            title = "Serum IL-6",
            width = 4, status = "danger", solidHeader = TRUE,
            plotlyOutput("il6_plot", height="280px")
          ),
          box(
            title = "sIL-6R (IL-6R Saturation by TCZ)",
            width = 4, status = "primary", solidHeader = TRUE,
            plotlyOutput("sil6r_plot", height="280px")
          ),
          box(
            title = "VEGF (Angiogenesis Marker)",
            width = 4, status = "warning", solidHeader = TRUE,
            plotlyOutput("vegf_plot", height="280px")
          )
        ),
        fluidRow(
          box(
            title = "Biomarker Reference Ranges",
            width = 12, status = "info", solidHeader = TRUE,
            tableOutput("bm_ref_tbl")
          )
        )
      ),

      # ── TAB 4: Disease Activity ─────────────────────────────────────────────
      tabItem(
        tabName = "disease",
        fluidRow(
          valueBoxOutput("da_wk12_box", width=3),
          valueBoxOutput("da_wk26_box", width=3),
          valueBoxOutput("da_wk52_box", width=3),
          valueBoxOutput("remiss_box",  width=3)
        ),
        fluidRow(
          box(
            title = "GCA Disease Activity Index Over Time",
            width = 7, status = "primary", solidHeader = TRUE,
            plotlyOutput("da_plot", height="340px")
          ),
          box(
            title = "Macrophage & Th17 Dynamics",
            width = 5, status = "warning", solidHeader = TRUE,
            plotlyOutput("immune_plot", height="340px")
          )
        ),
        fluidRow(
          box(
            title = "Relapse Risk Assessment",
            width = 12, status = "danger", solidHeader = TRUE,
            column(4,
              h4("Risk Factors for GCA Relapse"),
              tags$ul(
                tags$li("Rapid GC taper (major risk)"),
                tags$li("High baseline IL-6 (>100 pg/mL)"),
                tags$li("Persistent elevated CRP on GC"),
                tags$li("LV-GCA phenotype (higher relapse rate)"),
                tags$li("HLA-DRB1*04 genotype")
              )
            ),
            column(4,
              h4("Remission Definition (GiACTA)"),
              p("Absence of GCA symptoms AND signs PLUS:"),
              tags$ul(
                tags$li("CRP < 5 mg/L"),
                tags$li("ESR < 20 mm/hr"),
                tags$li("Disease Activity Index < 15%")
              ),
              p("Primary endpoint: Sustained remission at week 52")
            ),
            column(4,
              plotlyOutput("relapse_gauge", height="200px")
            )
          )
        )
      ),

      # ── TAB 5: Treatment Comparison ─────────────────────────────────────────
      tabItem(
        tabName = "compare",
        fluidRow(
          box(
            title = "Select Scenarios to Compare",
            width = 12, status = "primary", solidHeader = TRUE,
            column(2, checkboxInput("s1_on", "GC Mono Slow", TRUE)),
            column(2, checkboxInput("s2_on", "GC Mono Fast", TRUE)),
            column(2, checkboxInput("s3_on", "TCZ IV q4w", TRUE)),
            column(2, checkboxInput("s4_on", "TCZ SC qw", TRUE)),
            column(2, checkboxInput("s5_on", "TCZ SC q2w", TRUE)),
            column(2, actionButton("run_compare", "Run Comparison",
                                   class="btn-primary btn-lg", icon=icon("play")))
          )
        ),
        fluidRow(
          box(
            title = "CRP & ESR Comparison",
            width = 6, status = "danger", solidHeader = TRUE,
            plotlyOutput("compare_crp_esr", height="320px")
          ),
          box(
            title = "Disease Activity Comparison",
            width = 6, status = "primary", solidHeader = TRUE,
            plotlyOutput("compare_da", height="320px")
          )
        ),
        fluidRow(
          box(
            title = "52-Week Endpoint Summary Table",
            width = 12, status = "success", solidHeader = TRUE,
            DTOutput("compare_table")
          )
        )
      ),

      # ── TAB 6: Safety Monitoring ────────────────────────────────────────────
      tabItem(
        tabName = "safety",
        fluidRow(
          valueBoxOutput("cum_gc_box",  width=4),
          valueBoxOutput("bmd_box",     width=4),
          valueBoxOutput("infect_box",  width=4)
        ),
        fluidRow(
          box(
            title = "Cumulative GC Dose Over Time",
            width = 6, status = "danger", solidHeader = TRUE,
            plotlyOutput("cumgc_plot", height="300px")
          ),
          box(
            title = "Bone Mineral Density (% baseline)",
            width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("bmd_plot", height="300px")
          )
        ),
        fluidRow(
          box(
            title = "GC Side-Effect Risk Assessment",
            width = 6, status = "danger", solidHeader = TRUE,
            tableOutput("gc_risk_tbl")
          ),
          box(
            title = "Recommended Monitoring & Prophylaxis",
            width = 6, status = "success", solidHeader = TRUE,
            h4("Mandatory at GC Initiation"),
            tags$ul(
              tags$li("DXA scan (DEXA) — bone mineral density baseline"),
              tags$li("Fasting glucose + HbA1c"),
              tags$li("Blood pressure monitoring"),
              tags$li("Ophthalmology referral (IOP)"),
              tags$li("Pneumococcal + influenza vaccination")
            ),
            h4("Prophylactic Therapy"),
            tags$ul(
              tags$li("Calcium 1200 mg/day + Vitamin D 800 IU/day"),
              tags$li("Bisphosphonate (if T-score < -1.5 or high fracture risk)"),
              tags$li("TMP-SMX PCP prophylaxis (if prednisone >20 mg/day >4 wks)"),
              tags$li("PPI if NSAID co-use or high GI risk"),
              tags$li("Low-dose aspirin 81 mg/day (reduces ischemic events by ~50%)")
            )
          )
        )
      ),

      # ── TAB 7: Mechanistic Map ──────────────────────────────────────────────
      tabItem(
        tabName = "map",
        fluidRow(
          box(
            title = "GCA QSP Mechanistic Map — Pathways & Drug Actions",
            width = 12, status = "primary", solidHeader = TRUE,
            p("This map depicts the full immunopathogenesis of giant cell arteritis
               including vascular dendritic cell activation, Th1/Th17 immune responses,
               cytokine networks (IL-6, IFN-γ, IL-17A), vascular pathology,
               drug mechanisms (glucocorticoids, tocilizumab, abatacept, JAK inhibitors),
               and clinical endpoints."),
            tags$div(
              style = "overflow:auto; max-height:800px;",
              tags$img(
                src   = "gca_qsp_model.png",
                style = "max-width:100%; border-radius:8px; border:1px solid #DEE2E6;"
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "Cluster Legend",
            width = 12, status = "info", solidHeader = TRUE,
            column(3,
              tags$div(style="background:#D5E8D4;padding:8px;border-radius:6px;margin:4px;",
                "Genetic Risk & Trigger")
            ),
            column(3,
              tags$div(style="background:#FDEBD0;padding:8px;border-radius:6px;margin:4px;",
                "Innate Immunity (Macrophages, DCs)")
            ),
            column(3,
              tags$div(style="background:#D5F5E3;padding:8px;border-radius:6px;margin:4px;",
                "Adaptive Immunity (Th1, Th17, Treg)")
            ),
            column(3,
              tags$div(style="background:#FEF9E7;padding:8px;border-radius:6px;margin:4px;",
                "Cytokine Network (IL-6, IFN-γ, IL-17A)")
            ),
            column(3,
              tags$div(style="background:#FADBD8;padding:8px;border-radius:6px;margin:4px;",
                "Vascular Pathology (Temporal Artery)")
            ),
            column(3,
              tags$div(style="background:#F4ECF7;padding:8px;border-radius:6px;margin:4px;",
                "Clinical Manifestations & Ischemia")
            ),
            column(3,
              tags$div(style="background:#D5F5E3;padding:8px;border-radius:6px;margin:4px;",
                "Biomarkers (CRP, ESR, IL-6, PET/CT)")
            ),
            column(3,
              tags$div(style="background:#D6EAF8;padding:8px;border-radius:6px;margin:4px;",
                "Drug PK (Prednisone, Tocilizumab)")
            ),
            column(3,
              tags$div(style="background:#F0E6FF;padding:8px;border-radius:6px;margin:4px;",
                "Drug PD (GR, NFκB, JAK-STAT3)")
            ),
            column(3,
              tags$div(style="background:#F2F3F4;padding:8px;border-radius:6px;margin:4px;",
                "GC Side Effects (BMD, DM, Infect.)")
            )
          )
        )
      )
    )
  )
)

# ── SERVER ─────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Run main simulation ────────────────────────────────────────────────────
  sim_result <- reactive({
    req(input$pred_start, input$pred_taper, input$tcz_regimen, input$sim_weeks)

    pred_ev <- make_pred_events(input$pred_start, input$pred_taper)

    tcz_route <- 0
    tcz_ev    <- NULL
    if (input$tcz_regimen == "iv") {
      tcz_route <- 1
      n_doses   <- ceiling(input$sim_weeks / 4)
      tcz_ev    <- make_tcz_iv_events(n_doses = n_doses, bw = input$bw)
      tcz_ev$time <- tcz_ev$time + input$tcz_start_wk * 168
    } else if (input$tcz_regimen == "sc_qw") {
      tcz_route <- 2
      tcz_ev    <- make_tcz_sc_events("qw", n_doses = input$sim_weeks)
      tcz_ev$time <- tcz_ev$time + input$tcz_start_wk * 168
    } else if (input$tcz_regimen == "sc_q2w") {
      tcz_route <- 2
      tcz_ev    <- make_tcz_sc_events("q2w", n_doses = ceiling(input$sim_weeks/2))
      tcz_ev$time <- tcz_ev$time + input$tcz_start_wk * 168
    }

    ev <- bind_rows(pred_ev, tcz_ev)

    tgrid <- seq(0, input$sim_weeks * 168, by = 12)

    gca_mod %>%
      param(TCZ_route = tcz_route, bw = input$bw) %>%
      data_set(ev) %>%
      mrgsim(tgrid = tgrid, obsonly = TRUE) %>%
      as.data.frame() %>%
      mutate(time_wk = time / 168)
  })

  # ── TAB 1: Patient profile outputs ────────────────────────────────────────
  output$vb_crp <- renderValueBox({
    valueBox(
      paste0(input$baseline_crp, " mg/L"),
      "Baseline CRP",
      icon  = icon("vial"),
      color = if (input$baseline_crp > 50) "red" else "orange"
    )
  })

  output$vb_esr <- renderValueBox({
    valueBox(
      paste0(input$baseline_esr, " mm/hr"),
      "Baseline ESR",
      icon  = icon("tint"),
      color = if (input$baseline_esr > 50) "red" else "orange"
    )
  })

  output$vb_il6 <- renderValueBox({
    est_il6 <- round(input$baseline_crp * 1.1 + 5, 0)
    valueBox(
      paste0("~", est_il6, " pg/mL"),
      "Estimated IL-6 (active GCA)",
      icon  = icon("dna"),
      color = "red"
    )
  })

  output$acr_score_tbl <- renderTable({
    criteria <- data.frame(
      Criterion = c(
        "Age ≥50 years",
        "New headache",
        "Temporal artery abnormality",
        "ESR ≥50 mm/hr or CRP ≥10 mg/L",
        "Jaw or tongue claudication",
        "Visual symptoms",
        "Temporal artery biopsy (positive)",
        "Ultrasound halo sign",
        "FDG PET/CT (large vessel)"
      ),
      Points = c("+2", "+2", "+2", "+3", "+2", "+1", "+5", "+5", "+4"),
      Present = c(
        ifelse(input$age >= 50, "Yes (+2)", "No"),
        "Yes (+2)",
        ifelse(input$jaw_claud, "Yes (+2)", "No"),
        ifelse(input$baseline_esr >= 50 | input$baseline_crp >= 10, "Yes (+3)", "No"),
        ifelse(input$jaw_claud, "Yes (+2)", "No"),
        ifelse(input$visual_sx, "Yes (+1)", "No"),
        "Pending",
        "Pending",
        "Pending"
      )
    )
    criteria
  }, striped = TRUE, bordered = TRUE)

  output$acr_total_box <- renderValueBox({
    score <- 2 + 2 + ifelse(input$jaw_claud, 2, 0) +
      ifelse(input$baseline_esr >= 50 | input$baseline_crp >= 10, 3, 0) +
      ifelse(input$visual_sx, 1, 0)
    valueBox(
      paste0(score, " / 21"),
      if (score >= 6) "GCA Classification: POSITIVE (≥6)" else "GCA Classification: Uncertain",
      icon  = icon("check-circle"),
      color = if (score >= 6) "green" else "yellow"
    )
  })

  # ── TAB 2: PK plots ─────────────────────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    df <- sim_result()
    p1 <- ggplot(df, aes(time_wk, Pred_Cp)) +
      geom_line(color = "#E67E22", size = 1) +
      labs(x = "Weeks", y = "Prednisolone (mg/L)", title = "Prednisolone Plasma [Cp]") +
      theme_gca

    p2 <- ggplot(df, aes(time_wk, TCZ_Cp)) +
      geom_line(color = "#2980B9", size = 1) +
      labs(x = "Weeks", y = "TCZ (mg/L)", title = "Tocilizumab Plasma [Cp]") +
      theme_gca

    subplot(ggplotly(p1), ggplotly(p2), nrows = 2, shareX = TRUE)
  })

  output$pk_param_tbl <- renderTable({
    data.frame(
      Parameter = c("Prednisone ka", "Prednisolone CL", "Prednisolone Vd",
                    "TCZ V1 (central)", "TCZ V2 (periph)", "TCZ CL",
                    "TCZ T1/2 (α)", "TCZ SC F", "TCZ SC ka"),
      Value     = c("1.35 h⁻¹", "1.67 L/h", "32 L",
                    "3.72 L", "3.44 L", "0.011 L/h",
                    "~11 days", "80%", "0.0031 h⁻¹"),
      Reference = c("Bergrem 2005", "Rose 1981", "Rose 1981",
                    "Rau 2014", "Rau 2014", "Rau 2014",
                    "GiACTA PK", "GiACTA PK", "GiACTA PK")
    )
  }, striped = TRUE, bordered = TRUE)

  output$dose_schedule_tbl <- renderTable({
    data.frame(
      Week  = c(0, 2, 4, 8, 12, 16, 20, 26),
      Dose  = c(input$pred_start, 40, 30, 25, 20, 15, 10, 0),
      Route = rep("Oral (mg/day)", 8)
    )
  }, striped = TRUE, bordered = TRUE)

  # ── TAB 3: PD Biomarkers ──────────────────────────────────────────────────
  output$crp_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time_wk, CRP_obs)) +
      geom_line(color="#E74C3C", size=1.1) +
      geom_hline(yintercept=5, linetype="dashed", color="gray50") +
      annotate("text", x=max(df$time_wk)*0.9, y=8, label="ULN 5 mg/L", size=3) +
      labs(x="Weeks", y="CRP (mg/L)") + theme_gca
    ggplotly(p)
  })

  output$esr_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time_wk, ESR_obs)) +
      geom_line(color="#F39C12", size=1.1) +
      geom_hline(yintercept=20, linetype="dashed", color="gray50") +
      labs(x="Weeks", y="ESR (mm/hr)") + theme_gca
    ggplotly(p)
  })

  output$il6_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time_wk, IL6_obs)) +
      geom_line(color="#E74C3C", size=1.1) +
      geom_hline(yintercept=7, linetype="dashed", color="gray50") +
      labs(x="Weeks", y="IL-6 (pg/mL)") + theme_gca
    ggplotly(p)
  })

  output$sil6r_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time_wk, IL6R_bl)) +
      geom_line(color="#2980B9", size=1.1) +
      geom_hline(yintercept=80, linetype="dashed") +
      annotate("text", x=max(df$time_wk)*0.7, y=83, label="80% target", size=3) +
      labs(x="Weeks", y="IL-6R Saturation (%)") + theme_gca
    ggplotly(p)
  })

  output$vegf_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time_wk, VEGF_obs)) +
      geom_line(color="#9B59B6", size=1.1) +
      geom_hline(yintercept=150, linetype="dashed", color="gray50") +
      annotate("text", x=max(df$time_wk)*0.7, y=160, label="Active threshold", size=3) +
      labs(x="Weeks", y="VEGF (pg/mL)") + theme_gca
    ggplotly(p)
  })

  output$bm_ref_tbl <- renderTable({
    data.frame(
      Biomarker = c("CRP", "ESR", "IL-6", "sIL-6R", "VEGF", "Platelet"),
      Normal    = c("<5 mg/L", "<20 mm/hr (F), <15 (M)", "<7 pg/mL",
                    "20-60 ng/mL", "<150 pg/mL", "150-400 ×10⁹/L"),
      Active_GCA= c("50-200 mg/L", "70-120 mm/hr", "50-300 pg/mL",
                    "Elevated (60-150 ng/mL)", "150-500 pg/mL", "400-600 ×10⁹/L"),
      Remission = c("<5 mg/L", "<20 mm/hr", "<12.5 pg/mL",
                    "↑↑↑ on TCZ", "<150 pg/mL", "Normal")
    )
  }, striped = TRUE, bordered = TRUE)

  # ── TAB 4: Disease Activity ───────────────────────────────────────────────
  get_da_at_wk <- function(df, wk) {
    df %>% filter(abs(time_wk - wk) < 0.1) %>% slice(1) %>% pull(DA_obs)
  }

  output$da_wk12_box <- renderValueBox({
    df <- sim_result(); da <- round(get_da_at_wk(df, 12), 1)
    valueBox(paste0(da, "%"), "DA at Week 12", icon=icon("thermometer-half"),
             color = if(da < 30) "green" else if(da < 60) "yellow" else "red")
  })
  output$da_wk26_box <- renderValueBox({
    df <- sim_result(); da <- round(get_da_at_wk(df, 26), 1)
    valueBox(paste0(da, "%"), "DA at Week 26", icon=icon("thermometer-half"),
             color = if(da < 20) "green" else if(da < 40) "yellow" else "red")
  })
  output$da_wk52_box <- renderValueBox({
    df <- sim_result(); da <- round(get_da_at_wk(df, min(52, input$sim_weeks)), 1)
    valueBox(paste0(da, "%"), "DA at Week 52", icon=icon("thermometer-half"),
             color = if(da < 15) "green" else "red")
  })
  output$remiss_box <- renderValueBox({
    df <- sim_result()
    in_remiss <- df %>%
      filter(time_wk >= (input$sim_weeks - 4)) %>%
      summarise(r = mean(remiss, na.rm=TRUE)) %>% pull(r)
    valueBox(
      if(in_remiss > 0.5) "YES" else "NO",
      "Sustained Remission (last 4wk)",
      icon  = icon("check-circle"),
      color = if(in_remiss > 0.5) "green" else "red"
    )
  })

  output$da_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time_wk, DA_obs)) +
      geom_line(color="#2980B9", size=1.3) +
      geom_ribbon(aes(ymin=0, ymax=DA_obs), fill="#AED6F1", alpha=0.3) +
      geom_hline(yintercept=15, linetype="dashed", color="green4") +
      annotate("text", x=max(df$time_wk)*0.75, y=17,
               label="Remission (<15%)", color="green4", size=3.5) +
      labs(x="Weeks", y="Disease Activity (%)") + theme_gca
    ggplotly(p)
  })

  output$immune_plot <- renderPlotly({
    df <- sim_result() %>%
      select(time_wk, MacAct_obs, Th17_obs) %>%
      pivot_longer(c(MacAct_obs, Th17_obs), names_to="cell", values_to="value") %>%
      mutate(cell = ifelse(cell=="MacAct_obs","Macrophage Activation","Th17 Population"))
    p <- ggplot(df, aes(time_wk, value, color=cell)) +
      geom_line(size=1.1) +
      scale_color_manual(values=c("#E74C3C","#27AE60")) +
      labs(x="Weeks", y="Relative Level (normalized)") + theme_gca
    ggplotly(p)
  })

  output$relapse_gauge <- renderPlotly({
    df <- sim_result()
    da_end <- get_da_at_wk(df, min(52, input$sim_weeks))
    crp_end <- df %>% filter(abs(time_wk-min(52,input$sim_weeks))<0.1) %>%
      slice(1) %>% pull(CRP_obs)
    risk_pct <- min(100, round(da_end + crp_end * 0.3, 0))

    plot_ly(
      type = "indicator", mode = "gauge+number",
      value = risk_pct,
      title = list(text = "Relapse Risk Score"),
      gauge = list(
        axis = list(range = c(0, 100)),
        bar = list(color = "#E74C3C"),
        steps = list(
          list(range=c(0,30), color="#2ECC71"),
          list(range=c(30,60), color="#F39C12"),
          list(range=c(60,100), color="#E74C3C")
        ),
        threshold = list(line=list(color="black",width=4), thickness=0.75, value=30)
      )
    ) %>% layout(height=200, margin=list(t=30,b=10,l=10,r=10))
  })

  # ── TAB 5: Treatment Comparison ───────────────────────────────────────────
  compare_data <- eventReactive(input$run_compare, {
    tgrid <- seq(0, 52 * 168, by = 12)
    results <- list()

    if (input$s1_on) {
      ev <- make_pred_events(60, "slow")
      results[["GC Mono (Slow Taper)"]] <- gca_mod %>%
        param(TCZ_route=0) %>% data_set(ev) %>%
        mrgsim(tgrid=tgrid, obsonly=TRUE) %>% as.data.frame() %>%
        mutate(time_wk=time/168, scenario="GC Mono (Slow Taper)")
    }
    if (input$s2_on) {
      ev <- make_pred_events(60, "fast")
      results[["GC Mono (Fast Taper)"]] <- gca_mod %>%
        param(TCZ_route=0) %>% data_set(ev) %>%
        mrgsim(tgrid=tgrid, obsonly=TRUE) %>% as.data.frame() %>%
        mutate(time_wk=time/168, scenario="GC Mono (Fast Taper)")
    }
    if (input$s3_on) {
      ev <- bind_rows(make_pred_events(60,"slow"), make_tcz_iv_events(13, 70))
      results[["TCZ IV q4w + GC"]] <- gca_mod %>%
        param(TCZ_route=1) %>% data_set(ev) %>%
        mrgsim(tgrid=tgrid, obsonly=TRUE) %>% as.data.frame() %>%
        mutate(time_wk=time/168, scenario="TCZ IV q4w + GC")
    }
    if (input$s4_on) {
      ev <- bind_rows(make_pred_events(60,"slow"), make_tcz_sc_events("qw", 52))
      results[["TCZ SC qw + GC"]] <- gca_mod %>%
        param(TCZ_route=2) %>% data_set(ev) %>%
        mrgsim(tgrid=tgrid, obsonly=TRUE) %>% as.data.frame() %>%
        mutate(time_wk=time/168, scenario="TCZ SC qw + GC")
    }
    if (input$s5_on) {
      ev <- bind_rows(make_pred_events(60,"slow"), make_tcz_sc_events("q2w", 26))
      results[["TCZ SC q2w + GC"]] <- gca_mod %>%
        param(TCZ_route=2) %>% data_set(ev) %>%
        mrgsim(tgrid=tgrid, obsonly=TRUE) %>% as.data.frame() %>%
        mutate(time_wk=time/168, scenario="TCZ SC q2w + GC")
    }
    bind_rows(results)
  }, ignoreNULL = FALSE)

  output$compare_crp_esr <- renderPlotly({
    df <- compare_data()
    if (is.null(df) || nrow(df)==0) return(plotly_empty())
    p <- ggplot(df, aes(time_wk, CRP_obs, color=scenario)) +
      geom_line(size=0.9) +
      geom_hline(yintercept=5, linetype="dashed") +
      scale_color_manual(values=pal_scenarios, drop=FALSE) +
      labs(x="Weeks", y="CRP (mg/L)", title="CRP Over Time") + theme_gca
    ggplotly(p) %>% layout(legend=list(orientation="h", y=-0.3))
  })

  output$compare_da <- renderPlotly({
    df <- compare_data()
    if (is.null(df) || nrow(df)==0) return(plotly_empty())
    p <- ggplot(df, aes(time_wk, DA_obs, color=scenario)) +
      geom_line(size=0.9) +
      geom_hline(yintercept=15, linetype="dashed", color="green4") +
      scale_color_manual(values=pal_scenarios, drop=FALSE) +
      labs(x="Weeks", y="Disease Activity (%)", title="Disease Activity") + theme_gca
    ggplotly(p) %>% layout(legend=list(orientation="h", y=-0.3))
  })

  output$compare_table <- renderDT({
    df <- compare_data()
    if (is.null(df) || nrow(df)==0) return(datatable(data.frame()))
    df %>%
      filter(abs(time_wk - 52) < 0.15) %>%
      group_by(scenario) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(
        `CRP (mg/L)`    = round(CRP_obs, 1),
        `ESR (mm/hr)`   = round(ESR_obs, 1),
        `DA (%)`        = round(DA_obs, 1),
        `IL-6 (pg/mL)`  = round(IL6_obs, 1),
        `BMD (%)`       = round(BMD_obs, 1),
        `Cum GC (g)`    = round(CumGC_g, 2),
        `Remission`     = ifelse(remiss > 0.5, "YES", "NO")
      ) %>%
      select(scenario, `CRP (mg/L)`, `ESR (mm/hr)`, `DA (%)`,
             `IL-6 (pg/mL)`, `BMD (%)`, `Cum GC (g)`, `Remission`) %>%
      datatable(
        options  = list(pageLength=10, dom='t'),
        rownames = FALSE
      ) %>%
      formatStyle("Remission",
                  backgroundColor = styleEqual(c("YES","NO"), c("#D5F5E3","#FADBD8"))) %>%
      formatStyle("CRP (mg/L)",
                  backgroundColor = styleInterval(c(5,50), c("#D5F5E3","#FEF9E7","#FADBD8")))
  })

  # ── TAB 6: Safety outputs ─────────────────────────────────────────────────
  output$cum_gc_box <- renderValueBox({
    df <- sim_result()
    cg <- df %>% filter(time_wk==max(time_wk)) %>% pull(CumGC_g)
    valueBox(paste0(round(cg, 1), " g"), "Cumulative GC (prednisolone-eq.)",
             icon=icon("capsules"),
             color = if(cg < 3) "green" else if(cg < 8) "yellow" else "red")
  })
  output$bmd_box <- renderValueBox({
    df <- sim_result()
    bmd <- df %>% filter(time_wk==max(time_wk)) %>% pull(BMD_obs)
    valueBox(paste0(round(bmd, 1), "%"), "BMD at End (% baseline)",
             icon=icon("bone"),
             color = if(bmd > 97) "green" else if(bmd > 93) "yellow" else "red")
  })
  output$infect_box <- renderValueBox({
    df <- sim_result()
    cg <- df %>% filter(time_wk==max(time_wk)) %>% pull(CumGC_g)
    risk <- if(cg > 8) "High" else if(cg > 4) "Moderate" else "Low"
    valueBox(risk, "Infection Risk (GC-related)",
             icon=icon("biohazard"),
             color = if(risk=="Low") "green" else if(risk=="Moderate") "yellow" else "red")
  })

  output$cumgc_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time_wk, CumGC_g)) +
      geom_line(color="#E74C3C", size=1.2) +
      geom_hline(yintercept=3, linetype="dashed", color="green4") +
      geom_hline(yintercept=8, linetype="dashed", color="orange") +
      annotate("text", x=max(df$time_wk)*0.6, y=3.2, label="Low risk", color="green4", size=3) +
      annotate("text", x=max(df$time_wk)*0.6, y=8.2, label="High risk", color="orange", size=3) +
      labs(x="Weeks", y="Cumulative GC (g)") + theme_gca
    ggplotly(p)
  })

  output$bmd_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time_wk, BMD_obs)) +
      geom_line(color="#8E44AD", size=1.2) +
      geom_hline(yintercept=95, linetype="dashed", color="orange") +
      labs(x="Weeks", y="BMD (% baseline)") + theme_gca
    ggplotly(p)
  })

  output$gc_risk_tbl <- renderTable({
    df <- sim_result()
    cg <- df %>% filter(time_wk==max(time_wk)) %>% pull(CumGC_g)
    data.frame(
      `Side Effect`   = c("Osteoporosis/Fracture", "Steroid DM", "Hypertension",
                          "Infection (PCP, etc.)", "Cataract", "Adrenal Insufficiency",
                          "Steroid Myopathy", "Weight Gain"),
      `Risk Level`    = c(
        if(cg>6) "High (Bisphosphonate needed)" else "Moderate",
        if(cg>5) "Moderate (monitor HbA1c)" else "Low",
        "Moderate (monitor BP weekly)",
        if(cg>8) "High (TMP-SMX prophylaxis)" else "Moderate",
        if(cg>10) "Moderate (annual eye exam)" else "Low",
        if(cg>6) "Moderate (slow taper)" else "Low",
        if(cg>8) "Moderate" else "Low",
        "Moderate (diet counseling)"
      )
    )
  }, striped=TRUE, bordered=TRUE)
}

# ── Run app ───────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
