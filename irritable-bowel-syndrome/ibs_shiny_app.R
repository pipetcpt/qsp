## =============================================================================
## IBS QSP Shiny Dashboard — Irritable Bowel Syndrome
## 6 Tabs: Patient Profile · PK Profiles · PD Biomarkers ·
##         Clinical Endpoints · Scenario Comparison · Sensitivity Analysis
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(scales)

# ─────────────────────────────────────────────────────────────────────────────
# Inline model (same ODE code as ibs_mrgsolve_model.R — abbreviated for Shiny)
# ─────────────────────────────────────────────────────────────────────────────

ibs_model_code <- '
$PARAM
  KA_ALO=0.70, CL_ALO=50.0, VC_ALO=77.0, VP_ALO=52.0, Q_ALO=18.0, F_ALO=0.60, IC50_ALO=0.002,
  KA_LIN=0.10, CL_LIN=200.0, VC_LIN=5.0, F_LIN=0.02, EC50_LIN=0.05,
  KA_RIF=0.05, CL_RIF=300.0, VC_RIF=10.0, F_RIF=0.004, KLUM_RIF=0.08, EC50_RIF=0.10,
  KA_LOP=0.30, CL_LOP=14.0, VC_LOP=100.0, VP_LOP=80.0, Q_LOP=12.0, F_LOP=0.40, IC50_LOP=0.001,
  KA_TCA=0.35, CL_TCA=22.0, VC_TCA=1400.0, F_TCA=0.50, EC50_TCA=0.08,
  KA_TEGA=0.55, CL_TEGA=77.0, VC_TEGA=368.0, F_TEGA=0.10, EC50_TEGA=0.015,
  KSY_5HT=0.50, KRU_5HT=1.20, KDG_5HT=0.30, HT3R_BASE=0.20, HT4R_BASE=0.30,
  CTT_BASE=36.0, K_CTT=0.05, SEC_BASE=1.0, K_SEC=0.10,
  P_BASE=5.5, K_PAIN=0.15, VH_BASE=0.60, K_VH=0.008, K_CSS=0.004,
  MC_BASE=0.40, K_MC=0.12, LGI_BASE=0.45, K_LGI=0.05, PERM_BASE=0.50, K_PERM=0.03,
  DYS_BASE=0.55, K_DYS=0.02, STRESS_BASE=0.60, K_STRESS=0.04,
  IBSD_FLAG=1, IBS_SEVERITY=2,
  USE_ALO=0, USE_LIN=0, USE_RIF=0, USE_LOP=0, USE_TCA=0, USE_TEGA=0

$CMT
  GUT_ALO, CENT_ALO, PERIPH_ALO,
  LUM_LIN, CENT_LIN,
  LUM_RIF, CENT_RIF,
  GUT_LOP, CENT_LOP, PERIPH_LOP,
  CENT_TCA, CENT_TEGA,
  HT5_EC, HT3R_OCC, HT4R_TON,
  CTT, SEC_IDX, MMC_IDX,
  VH_IDX, SP_SENS, CENT_SENS, PAIN_NRS,
  MC_ACT, LGI_IDX, PERM_IDX,
  DYS_IDX, STRESS_IDX

$MAIN
  double CTT_init = IBSD_FLAG==1 ? CTT_BASE*0.72 : CTT_BASE*1.35;
  double PAIN_init = P_BASE + (IBS_SEVERITY-1)*1.5;
  if(NEWIND<=1){
    CTT_0=CTT_init; PAIN_NRS_0=PAIN_init;
    STRESS_IDX_0=STRESS_BASE+(IBS_SEVERITY-1)*0.10;
    SEC_IDX_0=IBSD_FLAG==1?SEC_BASE*1.4:SEC_BASE*0.7;
    MC_ACT_0=MC_BASE+(IBS_SEVERITY-1)*0.08;
    LGI_IDX_0=LGI_BASE+(IBS_SEVERITY-1)*0.08;
    VH_IDX_0=VH_BASE+(IBS_SEVERITY-1)*0.10;
    DYS_IDX_0=DYS_BASE; PERM_IDX_0=PERM_BASE;
    SP_SENS_0=0.45; CENT_SENS_0=0.35;
    HT5_EC_0=IBSD_FLAG==1?1.2:0.8;
    HT3R_OCC_0=HT3R_BASE; HT4R_TON_0=HT4R_BASE; MMC_IDX_0=0.60;
  }

$ODE
  double Cp_ALO=CENT_ALO/VC_ALO; double Cp_LOP=CENT_LOP/VC_LOP;
  double Cp_TCA=CENT_TCA/VC_TCA; double Cp_TEGA=CENT_TEGA/VC_TEGA;
  double E_ALO=USE_ALO*(Cp_ALO/(Cp_ALO+IC50_ALO));
  double E_LIN=USE_LIN*(LUM_LIN/(LUM_LIN+EC50_LIN));
  double E_RIF=USE_RIF*(LUM_RIF/(LUM_RIF+EC50_RIF));
  double E_LOP=USE_LOP*(Cp_LOP/(Cp_LOP+IC50_LOP));
  double E_TCA=USE_TCA*(Cp_TCA/(Cp_TCA+EC50_TCA));
  double E_TEGA=USE_TEGA*(Cp_TEGA/(Cp_TEGA+EC50_TEGA))*0.6;

  dxdt_GUT_ALO=-KA_ALO*GUT_ALO;
  dxdt_CENT_ALO=KA_ALO*F_ALO*GUT_ALO-(CL_ALO+Q_ALO)/VC_ALO*CENT_ALO+Q_ALO/VP_ALO*PERIPH_ALO;
  dxdt_PERIPH_ALO=Q_ALO/VC_ALO*CENT_ALO-Q_ALO/VP_ALO*PERIPH_ALO;
  dxdt_LUM_LIN=-KA_LIN*LUM_LIN-0.30*LUM_LIN;
  dxdt_CENT_LIN=KA_LIN*F_LIN*LUM_LIN-CL_LIN/VC_LIN*CENT_LIN;
  dxdt_LUM_RIF=-KLUM_RIF*LUM_RIF-KA_RIF*LUM_RIF;
  dxdt_CENT_RIF=KA_RIF*F_RIF*LUM_RIF-CL_RIF/VC_RIF*CENT_RIF;
  dxdt_GUT_LOP=-KA_LOP*GUT_LOP;
  dxdt_CENT_LOP=KA_LOP*F_LOP*GUT_LOP-(CL_LOP+Q_LOP)/VC_LOP*CENT_LOP+Q_LOP/VP_LOP*PERIPH_LOP;
  dxdt_PERIPH_LOP=Q_LOP/VC_LOP*CENT_LOP-Q_LOP/VP_LOP*PERIPH_LOP;
  dxdt_CENT_TCA=0-CL_TCA/VC_TCA*CENT_TCA;
  dxdt_CENT_TEGA=0-CL_TEGA/VC_TEGA*CENT_TEGA;

  double TPH1_up=1.0+0.6*LGI_IDX;
  double SERT_fun=1.0-0.85*E_TCA;
  dxdt_HT5_EC=KSY_5HT*TPH1_up-KRU_5HT*SERT_fun*HT5_EC-KDG_5HT*HT5_EC;
  double HT3R_drive=HT5_EC/(HT5_EC+1.0);
  dxdt_HT3R_OCC=0.30*(HT3R_drive-HT3R_OCC)-E_ALO*HT3R_OCC*1.5;
  double HT4R_drive=HT5_EC/(HT5_EC+1.5);
  dxdt_HT4R_TON=0.20*(HT4R_drive-HT4R_TON)+E_TEGA*(1.0-HT4R_TON)*0.5;

  double CTT_init2=IBSD_FLAG==1?CTT_BASE*0.72:CTT_BASE*1.35;
  double CTT_target=CTT_init2*(1.0-0.40*E_ALO)*(1.0-0.25*E_LOP)*(1.0-0.10*E_TCA);
  double CTT_lina_eff=(USE_LIN==1&&IBSD_FLAG==0)?-0.15*E_LIN*CTT:0;
  dxdt_CTT=K_CTT*(CTT_target-CTT)+CTT_lina_eff;

  double SEC_base2=IBSD_FLAG==1?SEC_BASE*1.4:SEC_BASE*0.7;
  double SEC_target=SEC_base2*(1.0+0.50*E_LIN*(1-IBSD_FLAG))*(1.0-0.35*E_LOP);
  dxdt_SEC_IDX=K_SEC*(SEC_target-SEC_IDX);

  dxdt_MMC_IDX=0.04*(0.60+0.25*E_RIF-MMC_IDX);

  double VH_target=VH_BASE+0.20*HT3R_OCC-0.30*E_ALO-0.15*E_TCA-0.12*E_LIN;
  VH_target=VH_target<0.10?0.10:VH_target>1.0?1.0:VH_target;
  dxdt_VH_IDX=K_VH*(VH_target-VH_IDX);

  double SP_target=0.40*VH_IDX+0.20*LGI_IDX+0.10*STRESS_IDX-0.20*E_TCA;
  SP_target=SP_target<0.05?0.05:SP_target>0.95?0.95:SP_target;
  dxdt_SP_SENS=K_CSS*(SP_target-SP_SENS);

  double CS_target=0.50*SP_SENS+0.25*STRESS_IDX-0.25*E_TCA;
  CS_target=CS_target<0.05?0.05:CS_target>0.95?0.95:CS_target;
  dxdt_CENT_SENS=K_CSS*0.8*(CS_target-CENT_SENS);

  double P_target=P_BASE+3.0*VH_IDX+2.5*SP_SENS+2.0*CENT_SENS+1.5*MC_ACT
                  +(IBS_SEVERITY-1)*1.5-3.5*E_ALO-2.0*E_TCA-1.5*E_LIN;
  P_target=P_target<0?0:P_target>10?10:P_target;
  dxdt_PAIN_NRS=K_PAIN*(P_target-PAIN_NRS);

  double MC_target=MC_BASE+0.35*STRESS_IDX+0.20*DYS_IDX-0.10*E_RIF;
  MC_target=MC_target<0.10?0.10:MC_target>1.0?1.0:MC_target;
  dxdt_MC_ACT=K_MC*(MC_target-MC_ACT);

  double LGI_target=LGI_BASE+0.30*MC_ACT+0.25*DYS_IDX+0.15*PERM_IDX-0.15*E_RIF;
  LGI_target=LGI_target<0.10?0.10:LGI_target>1.0?1.0:LGI_target;
  dxdt_LGI_IDX=K_LGI*(LGI_target-LGI_IDX);

  double PERM_target=PERM_BASE+0.25*LGI_IDX+0.20*STRESS_IDX-0.15*E_RIF;
  PERM_target=PERM_target<0.10?0.10:PERM_target>1.0?1.0:PERM_target;
  dxdt_PERM_IDX=K_PERM*(PERM_target-PERM_IDX);

  double DYS_target=DYS_BASE+0.20*STRESS_IDX-0.45*E_RIF;
  DYS_target=DYS_target<0.10?0.10:DYS_target>1.0?1.0:DYS_target;
  dxdt_DYS_IDX=K_DYS*(DYS_target-DYS_IDX);

  dxdt_STRESS_IDX=K_STRESS*(STRESS_BASE-0.10*E_TCA-STRESS_IDX);

$TABLE
  capture Cp_ALO=CENT_ALO/VC_ALO;
  capture Cp_LOP=CENT_LOP/VC_LOP;
  capture LumRIF=LUM_RIF;
  capture LumLIN=LUM_LIN;
  capture pain=PAIN_NRS;
  capture ctransit=CTT;
  capture secretion=SEC_IDX;
  capture mmc=MMC_IDX;
  capture viscHyp=VH_IDX;
  capture spinal=SP_SENS;
  capture central=CENT_SENS;
  capture mast=MC_ACT;
  capture inflam=LGI_IDX;
  capture permeab=PERM_IDX;
  capture dysbiosis=DYS_IDX;
  capture stress=STRESS_IDX;
  capture ht3r=HT3R_OCC;
  capture ht4r=HT4R_TON;
  capture serotonin=HT5_EC;
  capture ibs_sss=fmin(500,(pain/10.0)*100+(viscHyp*150)+(fabs(ctransit-36.0)/36.0)*100+(inflam*75)+(dysbiosis*75));
'

mod <- mcode("ibs_shiny", ibs_model_code, quiet = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# Helper: run one scenario
# ─────────────────────────────────────────────────────────────────────────────
run_sim <- function(input_list) {
  p <- do.call(param, c(list(mod), input_list$params))
  out <- mrgsim(p, events = input_list$events,
                tgrid = seq(0, 91*24, 6), obsonly = TRUE)
  df <- as.data.frame(out)
  df$time_wk <- df$time / 168
  df
}

# ─────────────────────────────────────────────────────────────────────────────
# Shiny UI
# ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "black",
  dashboardHeader(title = span("IBS QSP Dashboard", style = "font-size:16px")),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "tab_patient",  icon = icon("user")),
      menuItem("PK Profiles",         tabName = "tab_pk",       icon = icon("chart-line")),
      menuItem("PD Biomarkers",       tabName = "tab_pd",       icon = icon("flask")),
      menuItem("Clinical Endpoints",  tabName = "tab_clinical", icon = icon("stethoscope")),
      menuItem("Scenario Comparison", tabName = "tab_scenario", icon = icon("table")),
      menuItem("Sensitivity Analysis",tabName = "tab_sens",     icon = icon("sliders-h"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      body, .content-wrapper, .right-side { background-color: #0a0a14; }
      .box { background-color: #121230; border: 1px solid #2a2a50; }
      .box-header { background-color: #1a1a40; color: white; }
      label { color: #aaaacc !important; }
      .small-box { background-color: #1a1a40 !important; }
      h3, h4 { color: white; }
      .nav-tabs-custom > .nav-tabs > li.active { border-top-color: #2a9aff; }
    "))),

    tabItems(

      # ── Tab 1: Patient Profile ────────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Profile", width = 4, solidHeader = TRUE,
              status = "primary",
              selectInput("ibs_type",   "IBS Subtype",
                          choices = c("IBS-D (Diarrhea)" = 1,
                                      "IBS-C (Constipation)" = 0),
                          selected = 1),
              selectInput("severity",   "Symptom Severity",
                          choices = c("Mild (SSS 75–175)" = 1,
                                      "Moderate (SSS 175–300)" = 2,
                                      "Severe (SSS >300)" = 3),
                          selected = 2),
              sliderInput("stress_base","Baseline Stress Level (HPA)",
                          min=0.2, max=0.95, value=0.60, step=0.05),
              sliderInput("dys_base",  "Baseline Dysbiosis Score",
                          min=0.2, max=0.90, value=0.55, step=0.05),
              sliderInput("perm_base", "Baseline Gut Permeability",
                          min=0.2, max=0.90, value=0.50, step=0.05),
              actionButton("run_btn",  "Run Simulation",
                           class = "btn-primary", style="width:100%")
          ),
          box(title = "Drug Treatment", width = 4, solidHeader = TRUE,
              status = "info",
              checkboxInput("use_alo",  "Alosetron (IBS-D) 1 mg BID", FALSE),
              conditionalPanel("input.use_alo",
                sliderInput("alo_dose","Alosetron dose (mg)", 0.5, 1.0, 1.0, 0.25)),
              checkboxInput("use_lin",  "Linaclotide (IBS-C) 290 µg QD", FALSE),
              checkboxInput("use_rif",  "Rifaximin 550 mg TID × 14 d", FALSE),
              checkboxInput("use_lop",  "Loperamide 2 mg BID (IBS-D)", FALSE),
              checkboxInput("use_tca",  "Amitriptyline 25 mg QD (TCA)", FALSE),
              checkboxInput("use_tega", "Tegaserod 6 mg BID (IBS-C)", FALSE)
          ),
          box(title = "ROME IV Diagnostic Criteria", width = 4,
              solidHeader = TRUE, status = "warning",
              tableOutput("rome_table")
          )
        ),
        fluidRow(
          box(title = "IBS Symptom Radar (Baseline)", width = 6,
              plotlyOutput("radar_plot", height = "350px")),
          box(title = "Baseline Biomarker Summary", width = 6,
              tableOutput("baseline_table"))
        )
      ),

      # ── Tab 2: PK Profiles ───────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Drug Plasma / Luminal Concentrations (0–7 days)", width = 12,
              solidHeader = TRUE, status = "primary",
              plotlyOutput("pk_plot", height = "420px"))
        ),
        fluidRow(
          box(title = "PK Parameter Summary", width = 12,
              DT::dataTableOutput("pk_table"))
        )
      ),

      # ── Tab 3: PD Biomarkers ─────────────────────────────────────────────
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(title = "5-HT Signaling & Receptor Occupancy (13 weeks)",
              width = 6, solidHeader = TRUE, status = "success",
              plotlyOutput("pd_5ht", height = "350px")),
          box(title = "Visceral Hypersensitivity & Pain Sensitization",
              width = 6, solidHeader = TRUE, status = "danger",
              plotlyOutput("pd_pain_sens", height = "350px"))
        ),
        fluidRow(
          box(title = "Mucosal Immunity: Mast Cells & Low-Grade Inflammation",
              width = 6, solidHeader = TRUE, status = "warning",
              plotlyOutput("pd_immune", height = "350px")),
          box(title = "Gut Barrier (Permeability) & Dysbiosis",
              width = 6, solidHeader = TRUE, status = "info",
              plotlyOutput("pd_barrier", height = "350px"))
        )
      ),

      # ── Tab 4: Clinical Endpoints ─────────────────────────────────────────
      tabItem(tabName = "tab_clinical",
        fluidRow(
          box(title = "Abdominal Pain NRS (0–10)", width = 6,
              solidHeader = TRUE, status = "danger",
              plotlyOutput("clin_pain", height = "320px")),
          box(title = "Colonic Transit Time (h) & Secretion Index",
              width = 6, solidHeader = TRUE, status = "info",
              plotlyOutput("clin_motility", height = "320px"))
        ),
        fluidRow(
          box(title = "IBS-SSS Composite Score (0–500)", width = 6,
              solidHeader = TRUE, status = "primary",
              plotlyOutput("clin_sss", height = "320px")),
          box(title = "Clinical Endpoint Summary Table", width = 6,
              DT::dataTableOutput("clin_table"))
        )
      ),

      # ── Tab 5: Scenario Comparison ───────────────────────────────────────
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "Scenario Configuration", width = 3,
              solidHeader = TRUE, status = "primary",
              helpText("Select scenarios to compare (multi-select):"),
              checkboxGroupInput("sc_select", NULL,
                choices = c("No Treatment" = "NT",
                            "Alosetron 1mg BID" = "ALO",
                            "Linaclotide 290µg QD" = "LIN",
                            "Rifaximin 550mg TID×14d" = "RIF",
                            "Loperamide 2mg BID" = "LOP",
                            "Alosetron + TCA" = "ALO_TCA",
                            "Linaclotide + Tegaserod" = "LIN_TEGA"),
                selected = c("NT","ALO","LIN","RIF")),
              selectInput("sc_endpoint", "Endpoint",
                          choices = c("Pain NRS" = "pain",
                                      "IBS-SSS" = "ibs_sss",
                                      "CTT (h)" = "ctransit",
                                      "Inflammation" = "inflam",
                                      "Dysbiosis" = "dysbiosis"),
                          selected = "pain"),
              actionButton("run_sc", "Compare Scenarios",
                           class = "btn-success", style="width:100%")
          ),
          box(title = "Scenario Comparison Plot", width = 9,
              solidHeader = TRUE, status = "primary",
              plotlyOutput("sc_plot", height = "450px"))
        ),
        fluidRow(
          box(title = "13-Week Endpoint Summary Table", width = 12,
              DT::dataTableOutput("sc_table"))
        )
      ),

      # ── Tab 6: Sensitivity Analysis ──────────────────────────────────────
      tabItem(tabName = "tab_sens",
        fluidRow(
          box(title = "Stress Level Sensitivity", width = 6,
              solidHeader = TRUE, status = "warning",
              sliderInput("sens_stress_range", "Stress Range",
                          min=0.1, max=1.0, value=c(0.2, 0.9), step=0.05),
              actionButton("run_sens", "Run Sensitivity",
                           class = "btn-warning", style="width:100%"),
              br(),
              plotlyOutput("sens_plot_stress", height = "380px")),
          box(title = "Alosetron Dose–Response", width = 6,
              solidHeader = TRUE, status = "info",
              sliderInput("alo_dr_max", "Max Alosetron Dose (mg)",
                          min=0.5, max=4.0, value=2.0, step=0.5),
              actionButton("run_dr", "Run Dose–Response",
                           class = "btn-info", style="width:100%"),
              br(),
              plotlyOutput("sens_dr_plot", height = "380px"))
        ),
        fluidRow(
          box(title = "Sensitivity Heatmap: Stress × Dysbiosis → Pain NRS",
              width = 12, solidHeader = TRUE, status = "primary",
              plotlyOutput("sens_heatmap", height = "400px"))
        )
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# Shiny Server
# ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  dark_theme <- function() {
    list(
      paper_bgcolor = "#0a0a14",
      plot_bgcolor  = "#121230",
      font          = list(color = "white"),
      xaxis = list(gridcolor="#2a2a50", color="white"),
      yaxis = list(gridcolor="#2a2a50", color="white"),
      legend = list(bgcolor="#0a0a20", bordercolor="#2a2a50", font=list(color="white"))
    )
  }

  # ── Build param list from UI ────────────────────────────────────────────
  build_params <- function() {
    list(
      IBSD_FLAG    = as.integer(input$ibs_type),
      IBS_SEVERITY = as.integer(input$severity),
      STRESS_BASE  = input$stress_base,
      DYS_BASE     = input$dys_base,
      PERM_BASE    = input$perm_base,
      USE_ALO      = as.integer(input$use_alo),
      USE_LIN      = as.integer(input$use_lin),
      USE_RIF      = as.integer(input$use_rif),
      USE_LOP      = as.integer(input$use_lop),
      USE_TCA      = as.integer(input$use_tca),
      USE_TEGA     = as.integer(input$use_tega)
    )
  }

  build_events <- reactive({
    evs <- list()
    if (input$use_alo)  evs[[length(evs)+1]] <- ev(amt=1,    cmt="GUT_ALO",  ii=12, addl=13*7*2-1)
    if (input$use_lin)  evs[[length(evs)+1]] <- ev(amt=0.29, cmt="LUM_LIN",  ii=24, addl=13*7-1)
    if (input$use_rif)  evs[[length(evs)+1]] <- ev(amt=550,  cmt="LUM_RIF",  ii=8,  addl=3*14-1)
    if (input$use_lop)  evs[[length(evs)+1]] <- ev(amt=2,    cmt="GUT_LOP",  ii=12, addl=13*7*2-1)
    if (input$use_tca)  evs[[length(evs)+1]] <- ev(amt=25,   cmt="CENT_TCA", ii=24, addl=13*7-1)
    if (input$use_tega) evs[[length(evs)+1]] <- ev(amt=6,    cmt="CENT_TEGA",ii=12, addl=13*7*2-1)
    if (length(evs)==0) return(ev(amt=0, cmt="GUT_ALO", time=0))
    do.call(c, evs)
  })

  # ── Reactive simulation result ──────────────────────────────────────────
  sim_result <- eventReactive(input$run_btn, {
    withProgress(message = "Running ODE simulation...", {
      p_mod  <- do.call(param, c(list(mod), build_params()))
      out    <- mrgsim(p_mod, events = build_events(),
                       tgrid = seq(0, 91*24, 6), obsonly = TRUE)
      df     <- as.data.frame(out)
      df$time_wk <- df$time / 168
      df
    })
  }, ignoreNULL = FALSE)

  # ── Tab 1: ROME IV table ────────────────────────────────────────────────
  output$rome_table <- renderTable({
    tibble(
      `ROME IV Criterion` = c(
        "Recurrent abdominal pain",
        "≥1 day/week × 3 months",
        "Related to defecation",
        "Change in stool frequency",
        "Change in stool consistency",
        "Onset ≥6 months prior"
      ),
      Status = c("✓","✓","✓","✓","✓","✓")
    )
  }, striped=TRUE, bordered=TRUE, hover=TRUE)

  output$baseline_table <- renderTable({
    sev <- as.integer(input$severity)
    tibble(
      Biomarker = c("Pain NRS","CTT (h)","VH Index","IBS-SSS (est.)","HPA Stress",
                    "Dysbiosis","Gut Permeability","Mast Cell Activity",
                    "Low-Grade Inflammation","EC-cell 5-HT"),
      `Baseline Value` = c(
        round(5.5 + (sev-1)*1.5, 1),
        ifelse(input$ibs_type==1, round(36*0.72,1), round(36*1.35,1)),
        round(0.60 + (sev-1)*0.10, 2),
        round(200 + (sev-1)*75),
        input$stress_base, input$dys_base, input$perm_base,
        round(0.40 + (sev-1)*0.08, 2),
        round(0.45 + (sev-1)*0.08, 2),
        ifelse(input$ibs_type==1, 1.2, 0.8)
      )
    )
  }, striped=TRUE, bordered=TRUE)

  output$radar_plot <- renderPlotly({
    sev <- as.integer(input$severity)
    categories <- c("Pain","Motility\nAbnormality","Visceral\nHypersens.",
                    "Inflammation","Dysbiosis","Gut\nPermeability",
                    "Stress","EC-5HT")
    values <- c(
      (5.5 + (sev-1)*1.5) / 10,
      ifelse(input$ibs_type==1, 0.7, 0.65),
      0.60 + (sev-1)*0.10,
      0.45 + (sev-1)*0.08,
      input$dys_base,
      input$perm_base,
      input$stress_base,
      0.75
    )
    plot_ly(type='scatterpolar', fill='toself',
            r=c(values, values[1]), theta=c(categories, categories[1]),
            line=list(color='#2a9aff'), fillcolor='rgba(42,154,255,0.3)') %>%
      layout(title=list(text="IBS Disease Burden Radar", font=list(color="white")),
             paper_bgcolor="#0a0a14", plot_bgcolor="#121230",
             polar=list(bgcolor="#121230",
                        radialaxis=list(color="white", range=c(0,1)),
                        angularaxis=list(color="white")),
             font=list(color="white"))
  })

  # ── Tab 2: PK Profiles ──────────────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    df <- sim_result() %>% filter(time <= 7*24) %>%
      mutate(time_h = time)
    p <- plot_ly() %>%
      layout(title=list(text="Drug PK Profiles (first 7 days)", font=list(color="white")))

    if (input$use_alo) {
      p <- add_trace(p, data=df, x=~time_h, y=~Cp_ALO*1000, name="Alosetron (ng/mL)",
                     type='scatter', mode='lines', line=list(color='#2a9aff'))
    }
    if (input$use_lop) {
      p <- add_trace(p, data=df, x=~time_h, y=~Cp_LOP*1000, name="Loperamide (ng/mL)",
                     type='scatter', mode='lines', line=list(color='#cc44cc'))
    }
    if (input$use_lin) {
      p <- add_trace(p, data=df, x=~time_h, y=~LumLIN*1000, name="Linaclotide Lumen (ng equiv.)",
                     type='scatter', mode='lines', line=list(color='#ff8844'))
    }
    if (input$use_rif) {
      p <- add_trace(p, data=df, x=~time_h, y=~LumRIF, name="Rifaximin Lumen (µg/mL)",
                     type='scatter', mode='lines', line=list(color='#44cc44'))
    }
    p %>% layout(
      xaxis=list(title="Time (h)", color="white", gridcolor="#2a2a50"),
      yaxis=list(title="Concentration", color="white", gridcolor="#2a2a50"),
      paper_bgcolor="#0a0a14", plot_bgcolor="#121230",
      font=list(color="white")
    )
  })

  output$pk_table <- DT::renderDataTable({
    tibble(
      Drug = c("Alosetron","Linaclotide","Rifaximin","Loperamide",
               "Amitriptyline (TCA)","Tegaserod"),
      `Dose` = c("1 mg BID","290 µg QD","550 mg TID×14d","2 mg PRN",
                  "25 mg QD","6 mg BID"),
      `F (%)` = c("60","<2","<0.4","40","50","10"),
      `t½ (h)` = c("1.5","~0.5 (lumen)","~0.5 (lumen)","10","20","11"),
      `Vc (L)` = c("77","5","10","100","1400","368"),
      `Target` = c("5-HT3R","GC-C receptor","RNA polymerase β","µ-OR",
                   "SERT/NET","5-HT4R"),
      `IC50/EC50` = c("1.6 nM","0.05 nM (GC-C)","~0.1 µg/mL (MIC)","1 nM",
                      "0.08 µg/mL","45 nM"),
      `Key Effect` = c("↓motility, ↓pain","↑secretion, ↓pain","↓SIBO/dysbiosis",
                       "↓transit, ↓secretion","↓pain centrally","↑motility (IBS-C)")
    ) %>%
      DT::datatable(options=list(pageLength=8, dom='t'),
                    class='table-dark table-striped')
  })

  # ── Tab 3: PD Biomarkers ────────────────────────────────────────────────
  output$pd_5ht <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_wk) %>%
      add_trace(y=~serotonin, name="EC-5HT (a.u.)", type='scatter', mode='lines',
                line=list(color='#ff8844')) %>%
      add_trace(y=~ht3r, name="5-HT3R Occupancy", type='scatter', mode='lines',
                line=list(color='#ffcc44')) %>%
      add_trace(y=~ht4r, name="5-HT4R Tone", type='scatter', mode='lines',
                line=list(color='#44ff44')) %>%
      layout(title=list(text="5-HT Signaling", font=list(color="white")),
             xaxis=list(title="Week", color="white", gridcolor="#2a2a50"),
             yaxis=list(title="Value (a.u.)", color="white", gridcolor="#2a2a50"),
             paper_bgcolor="#0a0a14", plot_bgcolor="#121230",
             font=list(color="white"), legend=list(bgcolor="#0a0a20"))
  })

  output$pd_pain_sens <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_wk) %>%
      add_trace(y=~viscHyp, name="Visceral Hypersens.", type='scatter', mode='lines',
                line=list(color='#ff4444')) %>%
      add_trace(y=~spinal, name="Spinal Sensitization", type='scatter', mode='lines',
                line=list(color='#ff8888')) %>%
      add_trace(y=~central, name="Central Sensitization", type='scatter', mode='lines',
                line=list(color='#ffaaaa')) %>%
      layout(title=list(text="Pain Sensitization Cascade", font=list(color="white")),
             xaxis=list(title="Week", color="white", gridcolor="#2a2a50"),
             yaxis=list(title="Index (0–1)", color="white", gridcolor="#2a2a50"),
             paper_bgcolor="#0a0a14", plot_bgcolor="#121230",
             font=list(color="white"))
  })

  output$pd_immune <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_wk) %>%
      add_trace(y=~mast, name="Mast Cell Activity", type='scatter', mode='lines',
                line=list(color='#cc2222')) %>%
      add_trace(y=~inflam, name="Low-Grade Inflammation", type='scatter', mode='lines',
                line=list(color='#ff6644')) %>%
      layout(title=list(text="Mucosal Immune Activation", font=list(color="white")),
             xaxis=list(title="Week", color="white", gridcolor="#2a2a50"),
             yaxis=list(title="Index (0–1)", color="white", gridcolor="#2a2a50"),
             paper_bgcolor="#0a0a14", plot_bgcolor="#121230",
             font=list(color="white"))
  })

  output$pd_barrier <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_wk) %>%
      add_trace(y=~permeab, name="Gut Permeability", type='scatter', mode='lines',
                line=list(color='#9a6a2a')) %>%
      add_trace(y=~dysbiosis, name="Dysbiosis Score", type='scatter', mode='lines',
                line=list(color='#6a8a2a')) %>%
      add_trace(y=~stress, name="HPA Stress", type='scatter', mode='lines',
                line=list(color='#9a44cc')) %>%
      layout(title=list(text="Gut Barrier & Microbiome", font=list(color="white")),
             xaxis=list(title="Week", color="white", gridcolor="#2a2a50"),
             yaxis=list(title="Index (0–1)", color="white", gridcolor="#2a2a50"),
             paper_bgcolor="#0a0a14", plot_bgcolor="#121230",
             font=list(color="white"))
  })

  # ── Tab 4: Clinical Endpoints ────────────────────────────────────────────
  output$clin_pain <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_wk, y=~pain, type='scatter', mode='lines',
            line=list(color='#ff4444', width=2)) %>%
      add_hline(y=3.5, line=list(dash='dash', color='#ffcc00')) %>%
      layout(title=list(text="Abdominal Pain NRS", font=list(color="white")),
             xaxis=list(title="Week", color="white", gridcolor="#2a2a50"),
             yaxis=list(title="Pain NRS (0–10)", range=c(0,10), color="white",
                        gridcolor="#2a2a50"),
             paper_bgcolor="#0a0a14", plot_bgcolor="#121230",
             font=list(color="white"))
  })

  output$clin_motility <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_wk) %>%
      add_trace(y=~ctransit, name="Colonic Transit (h)", type='scatter', mode='lines',
                line=list(color='#2a9aff', width=2)) %>%
      add_trace(y=~secretion*20, name="Secretion ×20", type='scatter', mode='lines',
                line=list(color='#44ccff', width=2)) %>%
      add_hline(y=30, line=list(dash='dot', color='#44ff44')) %>%
      add_hline(y=42, line=list(dash='dot', color='#44ff44')) %>%
      layout(title=list(text="Colonic Motility & Secretion", font=list(color="white")),
             xaxis=list(title="Week", color="white", gridcolor="#2a2a50"),
             yaxis=list(title="CTT (h) / Secretion Index ×20",
                        color="white", gridcolor="#2a2a50"),
             paper_bgcolor="#0a0a14", plot_bgcolor="#121230",
             font=list(color="white"))
  })

  output$clin_sss <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_wk, y=~ibs_sss, type='scatter', mode='lines',
            fill='tozeroy', fillcolor='rgba(42,154,255,0.2)',
            line=list(color='#2a9aff', width=2)) %>%
      add_hline(y=175, line=list(dash='dash', color='#44ff44')) %>%
      add_hline(y=300, line=list(dash='dash', color='#ffcc00')) %>%
      layout(title=list(text="IBS-SSS Composite Score", font=list(color="white")),
             xaxis=list(title="Week", color="white", gridcolor="#2a2a50"),
             yaxis=list(title="IBS-SSS (0–500)", color="white", gridcolor="#2a2a50"),
             paper_bgcolor="#0a0a14", plot_bgcolor="#121230",
             font=list(color="white"))
  })

  output$clin_table <- DT::renderDataTable({
    df <- sim_result() %>%
      mutate(week = round(time_wk)) %>%
      filter(week %in% c(0,4,8,13)) %>%
      group_by(week) %>%
      summarise(
        `Pain NRS`   = round(mean(pain, na.rm=TRUE), 1),
        `CTT (h)`    = round(mean(ctransit, na.rm=TRUE), 1),
        `IBS-SSS`    = round(mean(ibs_sss, na.rm=TRUE)),
        `Inflammation`= round(mean(inflam, na.rm=TRUE), 3),
        `Dysbiosis`  = round(mean(dysbiosis, na.rm=TRUE), 3),
        `VH Index`   = round(mean(viscHyp, na.rm=TRUE), 3),
        .groups = "drop"
      ) %>%
      rename(Week = week)
    DT::datatable(df, options=list(dom='t', pageLength=6),
                  class='table-dark table-striped')
  })

  # ── Tab 5: Scenario Comparison ───────────────────────────────────────────
  sc_data <- eventReactive(input$run_sc, {
    scenario_map <- list(
      NT      = list(params=list(IBSD_FLAG=1,IBS_SEVERITY=2,USE_ALO=0,USE_LIN=0,USE_RIF=0,USE_LOP=0,USE_TCA=0,USE_TEGA=0), ev=ev(amt=0,cmt="GUT_ALO"), col="#aaaaaa", nm="No Treatment"),
      ALO     = list(params=list(IBSD_FLAG=1,IBS_SEVERITY=2,USE_ALO=1,USE_LIN=0,USE_RIF=0,USE_LOP=0,USE_TCA=0,USE_TEGA=0), ev=ev(amt=1,cmt="GUT_ALO",ii=12,addl=13*7*2-1), col="#2a9aff", nm="Alosetron 1mg BID"),
      LIN     = list(params=list(IBSD_FLAG=0,IBS_SEVERITY=2,USE_ALO=0,USE_LIN=1,USE_RIF=0,USE_LOP=0,USE_TCA=0,USE_TEGA=0), ev=ev(amt=0.29,cmt="LUM_LIN",ii=24,addl=13*7-1), col="#ff8844", nm="Linaclotide 290µg QD"),
      RIF     = list(params=list(IBSD_FLAG=1,IBS_SEVERITY=2,USE_ALO=0,USE_LIN=0,USE_RIF=1,USE_LOP=0,USE_TCA=0,USE_TEGA=0), ev=ev(amt=550,cmt="LUM_RIF",ii=8,addl=3*14-1), col="#44cc44", nm="Rifaximin 550mg TID"),
      LOP     = list(params=list(IBSD_FLAG=1,IBS_SEVERITY=2,USE_ALO=0,USE_LIN=0,USE_RIF=0,USE_LOP=1,USE_TCA=0,USE_TEGA=0), ev=ev(amt=2,cmt="GUT_LOP",ii=12,addl=13*7*2-1), col="#cc44cc", nm="Loperamide 2mg BID"),
      ALO_TCA = list(params=list(IBSD_FLAG=1,IBS_SEVERITY=3,USE_ALO=1,USE_LIN=0,USE_RIF=0,USE_LOP=0,USE_TCA=1,USE_TEGA=0), ev=c(ev(amt=1,cmt="GUT_ALO",ii=12,addl=13*7*2-1), ev(amt=25,cmt="CENT_TCA",ii=24,addl=13*7-1)), col="#ff4444", nm="Alosetron+TCA"),
      LIN_TEGA= list(params=list(IBSD_FLAG=0,IBS_SEVERITY=3,USE_ALO=0,USE_LIN=1,USE_RIF=0,USE_LOP=0,USE_TCA=0,USE_TEGA=1), ev=c(ev(amt=0.29,cmt="LUM_LIN",ii=24,addl=13*7-1), ev(amt=6,cmt="CENT_TEGA",ii=12,addl=13*7*2-1)), col="#ffcc00", nm="Linaclotide+Tegaserod")
    )
    selected_keys <- input$sc_select
    bind_rows(lapply(selected_keys, function(k) {
      sc <- scenario_map[[k]]
      p_mod <- do.call(param, c(list(mod), sc$params,
                                list(STRESS_BASE=input$stress_base,
                                     DYS_BASE=input$dys_base)))
      out <- mrgsim(p_mod, events=sc$ev,
                    tgrid=seq(0,91*24,12), obsonly=TRUE)
      df <- as.data.frame(out)
      df$time_wk <- df$time/168
      df$scenario <- sc$nm
      df$color    <- sc$col
      df
    }))
  }, ignoreNULL = FALSE)

  output$sc_plot <- renderPlotly({
    df <- sc_data()
    ep <- input$sc_endpoint
    df_wk <- df %>% mutate(week=round(time_wk)) %>%
      group_by(scenario, color, week) %>%
      summarise(val = mean(.data[[ep]], na.rm=TRUE), .groups="drop")
    colors <- setNames(unique(df_wk$color), unique(df_wk$scenario))
    p <- plot_ly()
    for (sc_nm in unique(df_wk$scenario)) {
      sc_d <- df_wk %>% filter(scenario == sc_nm)
      p <- add_trace(p, data=sc_d, x=~week, y=~val, name=sc_nm,
                     type='scatter', mode='lines',
                     line=list(color=colors[[sc_nm]], width=2))
    }
    y_lab <- switch(ep, pain="Pain NRS (0–10)", ibs_sss="IBS-SSS (0–500)",
                    ctransit="CTT (h)", inflam="Inflammation (0–1)",
                    dysbiosis="Dysbiosis (0–1)")
    p %>% layout(
      title=list(text=paste("Scenario Comparison:", y_lab), font=list(color="white")),
      xaxis=list(title="Week", color="white", gridcolor="#2a2a50"),
      yaxis=list(title=y_lab, color="white", gridcolor="#2a2a50"),
      paper_bgcolor="#0a0a14", plot_bgcolor="#121230",
      font=list(color="white"), legend=list(bgcolor="#0a0a20")
    )
  })

  output$sc_table <- DT::renderDataTable({
    df <- sc_data() %>%
      mutate(week=round(time_wk)) %>%
      filter(week == 13) %>%
      group_by(scenario) %>%
      summarise(
        `Pain NRS (wk13)` = round(mean(pain, na.rm=TRUE), 1),
        `CTT (h)`         = round(mean(ctransit, na.rm=TRUE), 1),
        `IBS-SSS`         = round(mean(ibs_sss, na.rm=TRUE)),
        `Inflammation`    = round(mean(inflam, na.rm=TRUE), 3),
        `Dysbiosis`       = round(mean(dysbiosis, na.rm=TRUE), 3),
        .groups = "drop"
      ) %>%
      rename(Scenario = scenario)
    DT::datatable(df, options=list(dom='t', pageLength=10),
                  class='table-dark table-striped')
  })

  # ── Tab 6: Sensitivity Analysis ──────────────────────────────────────────
  sens_stress_data <- eventReactive(input$run_sens, {
    stress_range <- seq(input$sens_stress_range[1],
                        input$sens_stress_range[2], by=0.05)
    bind_rows(lapply(stress_range, function(sl) {
      p_s <- param(mod, IBSD_FLAG=1, IBS_SEVERITY=2,
                   STRESS_BASE=sl, DYS_BASE=input$dys_base,
                   USE_ALO=0, USE_LIN=0, USE_RIF=0,
                   USE_LOP=0, USE_TCA=0, USE_TEGA=0)
      out <- mrgsim(p_s, events=ev(amt=0,cmt="GUT_ALO"),
                    tgrid=seq(0,13*168,168), obsonly=TRUE)
      df <- tail(as.data.frame(out), 1)
      data.frame(stress=sl, pain=df$pain, ibs_sss=df$ibs_sss, inflam=df$inflam)
    }))
  })

  output$sens_plot_stress <- renderPlotly({
    df <- sens_stress_data()
    plot_ly(df, x=~stress) %>%
      add_trace(y=~pain, name="Pain NRS", type='scatter', mode='lines+markers',
                line=list(color='#ff4444')) %>%
      add_trace(y=~ibs_sss/50, name="IBS-SSS/50", type='scatter', mode='lines+markers',
                line=list(color='#ffcc00')) %>%
      add_trace(y=~inflam*10, name="Inflammation×10", type='scatter', mode='lines+markers',
                line=list(color='#44ccff', dash='dash')) %>%
      layout(title=list(text="Stress Sensitivity (Week 13)", font=list(color="white")),
             xaxis=list(title="Stress Index", color="white", gridcolor="#2a2a50"),
             yaxis=list(title="Scaled Endpoint", color="white", gridcolor="#2a2a50"),
             paper_bgcolor="#0a0a14", plot_bgcolor="#121230",
             font=list(color="white"))
  })

  dr_data <- eventReactive(input$run_dr, {
    doses <- seq(0.125, input$alo_dr_max, length.out=12)
    bind_rows(lapply(doses, function(d) {
      ev_d <- ev(amt=d, cmt="GUT_ALO", ii=12, addl=13*7*2-1)
      p_d  <- param(mod, IBSD_FLAG=1, IBS_SEVERITY=2,
                    STRESS_BASE=input$stress_base,
                    USE_ALO=1, USE_LIN=0, USE_RIF=0,
                    USE_LOP=0, USE_TCA=0, USE_TEGA=0)
      out  <- mrgsim(p_d, events=ev_d, tgrid=seq(0,13*168,168), obsonly=TRUE)
      df   <- tail(as.data.frame(out), 1)
      data.frame(dose=d, pain=df$pain, ht3r=df$ht3r, ibs_sss=df$ibs_sss)
    }))
  })

  output$sens_dr_plot <- renderPlotly({
    df <- dr_data()
    plot_ly(df, x=~dose) %>%
      add_trace(y=~pain, name="Pain NRS", type='scatter', mode='lines+markers',
                line=list(color='#ff4444')) %>%
      add_trace(y=~ht3r*10, name="5-HT3R Occ.×10", type='scatter', mode='lines+markers',
                line=list(color='#2a9aff', dash='dash')) %>%
      layout(title=list(text="Alosetron Dose–Response (Week 13)", font=list(color="white")),
             xaxis=list(title="Alosetron BID Dose (mg)", color="white", gridcolor="#2a2a50"),
             yaxis=list(title="Endpoint", color="white", gridcolor="#2a2a50"),
             paper_bgcolor="#0a0a14", plot_bgcolor="#121230",
             font=list(color="white"))
  })

  output$sens_heatmap <- renderPlotly({
    stress_v <- seq(0.2, 0.9, by=0.1)
    dys_v    <- seq(0.2, 0.9, by=0.1)
    mat <- matrix(NA, length(stress_v), length(dys_v))
    for (i in seq_along(stress_v))
      for (j in seq_along(dys_v)) {
        p_hm <- param(mod, IBSD_FLAG=1, IBS_SEVERITY=2,
                      STRESS_BASE=stress_v[i], DYS_BASE=dys_v[j],
                      USE_ALO=0, USE_LIN=0, USE_RIF=0,
                      USE_LOP=0, USE_TCA=0, USE_TEGA=0)
        out  <- mrgsim(p_hm, events=ev(amt=0,cmt="GUT_ALO"),
                       tgrid=seq(0,13*168,168), obsonly=TRUE)
        mat[i,j] <- tail(as.data.frame(out)$pain, 1)
      }
    plot_ly(x=dys_v, y=stress_v, z=mat, type='heatmap',
            colorscale='RdBu', reversescale=TRUE,
            colorbar=list(title="Pain NRS")) %>%
      layout(title=list(text="Pain NRS Heatmap: Stress × Dysbiosis", font=list(color="white")),
             xaxis=list(title="Dysbiosis Score", color="white"),
             yaxis=list(title="Stress Index", color="white"),
             paper_bgcolor="#0a0a14", plot_bgcolor="#121230",
             font=list(color="white"))
  })
}

shinyApp(ui, server)
