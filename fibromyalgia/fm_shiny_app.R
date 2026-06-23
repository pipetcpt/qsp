## ============================================================
##  Fibromyalgia QSP — Interactive Shiny Dashboard
##  6 Tabs:
##    1. Patient Profile & Disease Baseline
##    2. Drug PK (plasma concentrations)
##    3. PD Key Markers (central sensitization, transmitters)
##    4. Clinical Endpoints (pain, FIQ, fatigue, depression)
##    5. Scenario Comparison
##    6. Biomarkers (CSF SP, microglia, cortisol, sleep)
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

## ---- Inline model code (same as fm_mrgsolve_model.R) -------
fm_code <- '
$PARAM
ka_DUL=0.80 CL_DUL=54.0 V1_DUL=1640 Q_DUL=18.0 V2_DUL=820 F_DUL=0.50
ka_PRE=1.30 CL_PRE=6.8  V_PRE=42.0  F_PRE=0.90
ka_MIL=0.90 CL_MIL=50.0 V_MIL=300.0 F_MIL=0.85
ka_TCA=0.70 CL_TCA=40.0 V_TCA=1500  F_TCA=0.48
IC50_SERT_DUL=0.003 IC50_NET_DUL=0.012
IC50_SERT_MIL=0.015 IC50_NET_MIL=0.018
IC50_SERT_TCA=0.010 IC50_NET_TCA=0.025 Emax_SNRI=1.0
IC50_PRE_alpha2d=0.05 Emax_PRE=0.70
kprod_NGF=0.05 kdeg_NGF=0.10 kact_TRPV1=0.20 kdeg_PGE2=0.30 kprod_PGE2=0.08
kstim_DRG=0.40 kdeg_DRG=0.50
kprod_SP=0.15 kdeg_SP=0.20 kWU=0.05 kWU_decay=0.02
kLTP=0.08 kLTP_decay=0.005 kNMDA_act=0.12 kNMDA_decay=0.08 Emax_inhib=0.70
ksyn_NE=0.30 ksyn_5HT=0.25 kdeg_NE=0.40 kdeg_5HT=0.35 kdesc_NE=0.20 kdesc_5HT=0.15
kprod_CRH=0.50 kdeg_CRH=0.80 kprod_ACTH=0.40 kdeg_ACTH=0.60
kprod_CORT=0.30 kdeg_CORT=0.20 kfb_CORT=0.60
kSNS_base=0.60 kSNS_stress=0.20 kSNS_decay=0.40 kHRV_base=0.80
kaden_prod=0.15 kaden_clear=0.12 kSWS_drive=0.20 kSWS_decay=0.25
kSWS_pain_inh=0.30 kSWS_TCA=0.15
kprod_MG=0.10 kdeg_MG=0.12 kprod_IL1b=0.18 kdeg_IL1b=0.25 kMG_cortisol=0.15
k_pain_LTP=0.40 k_pain_SP=0.20 k_FIQ_pain=0.35 k_FIQ_sleep=0.25 k_FIQ_dep=0.20
k_fatigue=0.30 k_dep_LTP=0.15
pain_base=5.0 FIQ_base=55.0 fatigue_base=60.0
use_DUL=0 use_PRE=0 use_MIL=0 use_TCA=0

$CMT
DUL_gut DUL_cent DUL_peri PRE_gut PRE_cent MIL_gut MIL_cent TCA_gut TCA_cent
NGF PGE2 DRG_act SP_csf NMDA_state WindUp LTP_cs NE_syn SHT_syn DPMS
CRH ACTH CORT SNS_tone Adenosine SWS_depth MG_act IL1b_sp
Pain_score FIQ_score Fatigue_VAS Depression_score

$ODE
double ka_DUL_eff=use_DUL*ka_DUL;
dxdt_DUL_gut=-ka_DUL_eff*DUL_gut;
dxdt_DUL_cent=ka_DUL_eff*DUL_gut*F_DUL-(CL_DUL/V1_DUL)*DUL_cent-(Q_DUL/V1_DUL)*DUL_cent+(Q_DUL/V2_DUL)*DUL_peri;
dxdt_DUL_peri=(Q_DUL/V1_DUL)*DUL_cent-(Q_DUL/V2_DUL)*DUL_peri;
double Cp_DUL=DUL_cent/V1_DUL;
double ka_PRE_eff=use_PRE*ka_PRE;
dxdt_PRE_gut=-ka_PRE_eff*PRE_gut;
dxdt_PRE_cent=ka_PRE_eff*PRE_gut*F_PRE-(CL_PRE/V_PRE)*PRE_cent;
double Cp_PRE=PRE_cent/V_PRE;
double ka_MIL_eff=use_MIL*ka_MIL;
dxdt_MIL_gut=-ka_MIL_eff*MIL_gut;
dxdt_MIL_cent=ka_MIL_eff*MIL_gut*F_MIL-(CL_MIL/V_MIL)*MIL_cent;
double Cp_MIL=MIL_cent/V_MIL;
double ka_TCA_eff=use_TCA*ka_TCA;
dxdt_TCA_gut=-ka_TCA_eff*TCA_gut;
dxdt_TCA_cent=ka_TCA_eff*TCA_gut*F_TCA-(CL_TCA/V_TCA)*TCA_cent;
double Cp_TCA=TCA_cent/V_TCA;
double inh_SERT=Emax_SNRI*(Cp_DUL/(IC50_SERT_DUL+Cp_DUL)+Cp_MIL/(IC50_SERT_MIL+Cp_MIL)+Cp_TCA/(IC50_SERT_TCA+Cp_TCA));
if(inh_SERT>1.0) inh_SERT=1.0;
double inh_NET=Emax_SNRI*(Cp_DUL/(IC50_NET_DUL+Cp_DUL)+Cp_MIL/(IC50_NET_MIL+Cp_MIL)+Cp_TCA/(IC50_NET_TCA+Cp_TCA));
if(inh_NET>1.0) inh_NET=1.0;
double inh_Ca=Emax_PRE*Cp_PRE/(IC50_PRE_alpha2d+Cp_PRE);
double eff_TCA_sleep=Cp_TCA/(0.05+Cp_TCA);
dxdt_NGF=kprod_NGF*1.5-kdeg_NGF*NGF;
dxdt_PGE2=kprod_PGE2*(1.0+IL1b_sp*0.5)-kdeg_PGE2*PGE2;
double TRPV1_sens=kact_TRPV1*(PGE2+NGF*0.5);
dxdt_DRG_act=kstim_DRG*TRPV1_sens-kdeg_DRG*DRG_act;
dxdt_SP_csf=kprod_SP*DRG_act*(1.0-inh_Ca*0.6)-kdeg_SP*SP_csf;
dxdt_NMDA_state=kNMDA_act*(SP_csf+DRG_act*0.5)*(1.0-NMDA_state)-kNMDA_decay*NMDA_state;
dxdt_WindUp=kWU*DRG_act*NMDA_state-kWU_decay*WindUp;
dxdt_LTP_cs=kLTP*WindUp*(1.0+IL1b_sp*0.4)-kLTP_decay*LTP_cs-LTP_cs*DPMS*Emax_inhib;
dxdt_NE_syn=ksyn_NE*(1.0+inh_NET*2.0)-kdeg_NE*(1.0-inh_NET)*NE_syn;
dxdt_SHT_syn=ksyn_5HT*(1.0+inh_SERT*2.0)-kdeg_5HT*(1.0-inh_SERT)*SHT_syn;
dxdt_DPMS=kdesc_NE*NE_syn+kdesc_5HT*SHT_syn-0.30*DPMS-DPMS*SNS_tone*0.10;
dxdt_CRH=kprod_CRH*(1.0+LTP_cs*0.5+SNS_tone*0.3)-kdeg_CRH*CRH-kfb_CORT*CORT*CRH;
dxdt_ACTH=kprod_ACTH*CRH-kdeg_ACTH*ACTH;
dxdt_CORT=kprod_CORT*ACTH-kdeg_CORT*CORT;
dxdt_SNS_tone=kSNS_base+kSNS_stress*(LTP_cs+0.5*(1.0-CORT*0.5))-kSNS_decay*SNS_tone;
dxdt_Adenosine=kaden_prod-kaden_clear*Adenosine+SNS_tone*0.05;
dxdt_SWS_depth=kSWS_drive*Adenosine-kSWS_decay*SWS_depth-kSWS_pain_inh*LTP_cs*SWS_depth+kSWS_TCA*eff_TCA_sleep;
dxdt_MG_act=kprod_MG*DRG_act*(1.0+WindUp*0.5)-kdeg_MG*MG_act-kMG_cortisol*CORT*MG_act;
dxdt_IL1b_sp=kprod_IL1b*MG_act-kdeg_IL1b*IL1b_sp;
double pain_target=pain_base+k_pain_LTP*LTP_cs*6.0+k_pain_SP*SP_csf*2.0-DPMS*4.0-(inh_NET+inh_SERT)*2.0;
if(pain_target<0) pain_target=0; if(pain_target>10) pain_target=10;
dxdt_Pain_score=0.15*(pain_target-Pain_score);
double FIQ_target=FIQ_base+k_FIQ_pain*(Pain_score-5.0)*8.0+k_FIQ_sleep*(1.0-SWS_depth)*20.0+k_FIQ_dep*Depression_score*0.5;
if(FIQ_target<0) FIQ_target=0; if(FIQ_target>100) FIQ_target=100;
dxdt_FIQ_score=0.10*(FIQ_target-FIQ_score);
double fatigue_tgt=fatigue_base+k_fatigue*(1.0-SWS_depth)*30.0-inh_NET*20.0;
if(fatigue_tgt<0) fatigue_tgt=0; if(fatigue_tgt>100) fatigue_tgt=100;
dxdt_Fatigue_VAS=0.12*(fatigue_tgt-Fatigue_VAS);
double dep_target=8.0+k_dep_LTP*LTP_cs*12.0+(1.0-SWS_depth)*4.0-(inh_SERT+inh_NET)*5.0;
if(dep_target<0) dep_target=0; if(dep_target>27) dep_target=27;
dxdt_Depression_score=0.08*(dep_target-Depression_score);

$TABLE
double Cp_DUL_out=DUL_cent/V1_DUL;
double Cp_PRE_out=PRE_cent/V_PRE;
double Cp_MIL_out=MIL_cent/V_MIL;
double Cp_TCA_out=TCA_cent/V_TCA;
double inh_SERT_pct=100*(Cp_DUL_out/(IC50_SERT_DUL+Cp_DUL_out)+Cp_MIL_out/(IC50_SERT_MIL+Cp_MIL_out));
double inh_NET_pct=100*(Cp_DUL_out/(IC50_NET_DUL+Cp_DUL_out)+Cp_MIL_out/(IC50_NET_MIL+Cp_MIL_out));
double Ca_block_pct=100*Emax_PRE*Cp_PRE_out/(IC50_PRE_alpha2d+Cp_PRE_out);

$CAPTURE
Cp_DUL_out Cp_PRE_out Cp_MIL_out Cp_TCA_out
inh_SERT_pct inh_NET_pct Ca_block_pct
SP_csf NMDA_state WindUp LTP_cs DPMS NE_syn SHT_syn
CRH ACTH CORT SNS_tone SWS_depth Adenosine MG_act IL1b_sp
Pain_score FIQ_score Fatigue_VAS Depression_score
'

## Precompile once
fm_mod <- mcode("fm_shiny", fm_code, quiet = TRUE)

FM_init_vals <- c(
  NGF=1.5, PGE2=1.2, DRG_act=0.8, SP_csf=1.8, NMDA_state=0.4,
  WindUp=0.3, LTP_cs=0.55, NE_syn=0.7, SHT_syn=0.6, DPMS=0.35,
  CRH=0.9, ACTH=0.85, CORT=0.8, SNS_tone=0.75,
  Adenosine=1.2, SWS_depth=0.35, MG_act=0.65, IL1b_sp=0.55,
  Pain_score=6.5, FIQ_score=68, Fatigue_VAS=72, Depression_score=12
)

run_sim <- function(mod, use_DUL=0, use_PRE=0, use_MIL=0, use_TCA=0,
                    dose_DUL=60, dose_PRE=150, dose_MIL=50, dose_TCA=25,
                    weeks=12, severity=1.0) {
  init_adj <- FM_init_vals * severity
  evs <- NULL
  if (use_DUL) evs <- c(evs, list(ev(amt=dose_DUL, ii=24, addl=weeks*7-1, cmt="DUL_gut")))
  if (use_PRE) evs <- c(evs, list(ev(amt=dose_PRE, ii=12, addl=weeks*14-1, cmt="PRE_gut")))
  if (use_MIL) evs <- c(evs, list(ev(amt=dose_MIL, ii=12, addl=weeks*14-1, cmt="MIL_gut")))
  if (use_TCA) evs <- c(evs, list(ev(amt=dose_TCA, ii=24, addl=weeks*7-1, cmt="TCA_gut", time=22)))

  m <- mod %>%
    init(as.list(init_adj)) %>%
    param(use_DUL=use_DUL, use_PRE=use_PRE, use_MIL=use_MIL, use_TCA=use_TCA)

  if (!is.null(evs)) {
    combined_ev <- do.call(rbind, evs)
    m <- m %>% ev(combined_ev)
  }

  m %>%
    mrgsim(end=weeks*7*24, delta=6) %>%
    as.data.frame() %>%
    mutate(time_days = time / 24)
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "Fibromyalgia QSP Model"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "tab_patient",   icon = icon("user")),
      menuItem("Drug PK",            tabName = "tab_pk",        icon = icon("flask")),
      menuItem("PD Key Markers",     tabName = "tab_pd",        icon = icon("brain")),
      menuItem("Clinical Endpoints", tabName = "tab_clinical",  icon = icon("heartbeat")),
      menuItem("Scenario Comparison",tabName = "tab_scenario",  icon = icon("chart-bar")),
      menuItem("Biomarkers",         tabName = "tab_biomarker", icon = icon("vial"))
    ),
    hr(),
    h5("Drug Selection", style="padding-left:15px; color:#ddd;"),
    checkboxInput("use_DUL", "Duloxetine",   FALSE),
    conditionalPanel("input.use_DUL",
      sliderInput("dose_DUL", "Dose (mg QD)", 20, 120, 60, step=20)),
    checkboxInput("use_PRE", "Pregabalin",   FALSE),
    conditionalPanel("input.use_PRE",
      sliderInput("dose_PRE", "Dose (mg BID)", 75, 300, 150, step=75)),
    checkboxInput("use_MIL", "Milnacipran",  FALSE),
    conditionalPanel("input.use_MIL",
      sliderInput("dose_MIL", "Dose (mg BID)", 25, 100, 50, step=25)),
    checkboxInput("use_TCA", "Amitriptyline",FALSE),
    conditionalPanel("input.use_TCA",
      sliderInput("dose_TCA", "Dose (mg QHS)", 10, 75, 25, step=5)),
    hr(),
    sliderInput("weeks",    "Duration (weeks)", 4, 24, 12),
    sliderInput("severity", "Disease Severity", 0.5, 1.5, 1.0, step=0.1),
    actionButton("run_btn", "Run Simulation", class="btn-primary btn-block",
                 icon=icon("play"))
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .main-header .logo { font-weight: bold; font-size: 16px; }
      .box-title { font-weight: bold; }
      .kpi-box { text-align:center; padding:12px; border-radius:8px;
                 background:#fff; border:2px solid #9B59B6; margin:6px; }
      .kpi-val  { font-size:2em; font-weight:bold; color:#8E44AD; }
      .kpi-lbl  { font-size:0.9em; color:#555; }
    "))),

    tabItems(

      ## TAB 1 — Patient Profile
      tabItem("tab_patient",
        fluidRow(
          box(title="Fibromyalgia — Disease Overview", width=12, status="purple",
            fluidRow(
              column(4,
                h4("Diagnostic Criteria (2016 ACR)"),
                tags$ul(
                  tags$li("Widespread pain index (WPI) ≥ 7 + SS ≥ 5, or WPI 4-6 + SS ≥ 9"),
                  tags$li("Symptoms ≥ 3 months, generalized pain (≥4/5 regions)"),
                  tags$li("Rule out other diagnoses")
                ),
                h4("Key Pathophysiology"),
                tags$ul(
                  tags$li("Central sensitization (spinal LTP, wind-up)"),
                  tags$li("Descending pain inhibition deficit (DPMS ↓)"),
                  tags$li("Neuroinflammation (microglia, IL-1β, TNF-α)"),
                  tags$li("HPA axis dysregulation (blunted cortisol)"),
                  tags$li("Non-restorative sleep (alpha-delta intrusion)"),
                  tags$li("Autonomic imbalance (SNS↑, HRV↓)")
                )
              ),
              column(4,
                h4("QSP Model Structure"),
                tags$ul(
                  tags$li("30 ODE compartments"),
                  tags$li("Drug PK: 4 drugs, 2-compartment DUL, 1-cpt PRE/MIL/TCA"),
                  tags$li("Peripheral: NGF, PGE2, DRG afferent activity"),
                  tags$li("Spinal: SP (CSF), NMDA state, wind-up, LTP"),
                  tags$li("Brain: NE/5-HT synaptic pools, DPMS"),
                  tags$li("HPA: CRH-ACTH-cortisol loop"),
                  tags$li("ANS: SNS tone, HRV"),
                  tags$li("Sleep: adenosine, SWS depth"),
                  tags$li("Immune: microglia, IL-1β")
                )
              ),
              column(4,
                h4("Approved Pharmacotherapy"),
                tableOutput("drug_table")
              )
            )
          )
        ),
        fluidRow(
          box(title="Baseline Disease Metrics (Simulated)", width=12, status="warning",
            fluidRow(
              column(2, div(class="kpi-box",
                div(class="kpi-val", "6.5"), div(class="kpi-lbl", "Pain NRS"))),
              column(2, div(class="kpi-box",
                div(class="kpi-val", "68"), div(class="kpi-lbl", "FIQR"))),
              column(2, div(class="kpi-box",
                div(class="kpi-val", "72"), div(class="kpi-lbl", "Fatigue VAS"))),
              column(2, div(class="kpi-box",
                div(class="kpi-val", "12"), div(class="kpi-lbl", "PHQ-9"))),
              column(2, div(class="kpi-box",
                div(class="kpi-val", "0.35"), div(class="kpi-lbl", "SWS Depth"))),
              column(2, div(class="kpi-box",
                div(class="kpi-val", "0.55"), div(class="kpi-lbl", "LTP Index")))
            )
          )
        )
      ),

      ## TAB 2 — Drug PK
      tabItem("tab_pk",
        fluidRow(
          box(title="Plasma Concentration — Time Profiles", width=12, status="blue",
              plotOutput("pk_plot", height="450px"))
        ),
        fluidRow(
          box(title="PK Parameters Summary", width=6, status="info",
              tableOutput("pk_table")),
          box(title="Target Occupancy (PD)", width=6, status="info",
              plotOutput("occ_plot", height="250px"))
        )
      ),

      ## TAB 3 — PD Key Markers
      tabItem("tab_pd",
        fluidRow(
          box(title="Central Sensitization Index (Spinal LTP)", width=6, status="purple",
              plotOutput("ltp_plot", height="280px")),
          box(title="Descending Pain Modulation (DPMS)", width=6, status="purple",
              plotOutput("dpms_plot", height="280px"))
        ),
        fluidRow(
          box(title="Synaptic NE & 5-HT (CNS pools)", width=6, status="info",
              plotOutput("mono_plot", height="280px")),
          box(title="Wind-up & NMDA State", width=6, status="warning",
              plotOutput("wu_plot", height="280px"))
        )
      ),

      ## TAB 4 — Clinical Endpoints
      tabItem("tab_clinical",
        fluidRow(
          box(title="NRS Pain Score (0–10)", width=6, status="danger",
              plotOutput("pain_plot", height="280px")),
          box(title="FIQR Score (0–100)", width=6, status="warning",
              plotOutput("fiq_plot", height="280px"))
        ),
        fluidRow(
          box(title="Fatigue VAS (0–100)", width=6, status="olive",
              plotOutput("fatigue_plot", height="280px")),
          box(title="Depression PHQ-9 (0–27)", width=6, status="primary",
              plotOutput("dep_plot", height="280px"))
        ),
        fluidRow(
          box(title="Week-12 Outcome Summary", width=12, status="success",
              DTOutput("outcome_table"))
        )
      ),

      ## TAB 5 — Scenario Comparison
      tabItem("tab_scenario",
        fluidRow(
          box(title="Multi-Drug Scenario Comparison", width=12, status="purple",
            p("Compare all treatment scenarios at Week 12 (fixed doses):
               DUL 60mg QD | PRE 150mg BID | MIL 50mg BID | TCA 25mg QHS"),
            actionButton("run_all", "Run All Scenarios", class="btn-warning",
                         icon=icon("sync")),
            hr(),
            plotOutput("scenario_plot", height="500px")
          )
        ),
        fluidRow(
          box(title="Responder Analysis (≥30% / ≥50% Pain Reduction)", width=12,
              status="success", DTOutput("resp_table"))
        )
      ),

      ## TAB 6 — Biomarkers
      tabItem("tab_biomarker",
        fluidRow(
          box(title="CSF Substance P (central sensitization biomarker)", width=6,
              status="danger", plotOutput("sp_plot", height="260px")),
          box(title="Microglia Activation (neuroinflammation)", width=6,
              status="warning", plotOutput("mg_plot", height="260px"))
        ),
        fluidRow(
          box(title="HPA Axis: Cortisol (normalized)", width=6, status="info",
              plotOutput("cort_plot", height="260px")),
          box(title="Slow-Wave Sleep Depth", width=6, status="primary",
              plotOutput("sws_plot", height="260px"))
        ),
        fluidRow(
          box(title="Sympathetic Tone (ANS)", width=6, status="olive",
              plotOutput("sns_plot", height="260px")),
          box(title="IL-1β Spinal (neuroinflammation)", width=6, status="danger",
              plotOutput("il1b_plot", height="260px"))
        )
      )
    )
  )
)

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

  ## Drug table
  output$drug_table <- renderTable({
    data.frame(
      Drug       = c("Duloxetine","Pregabalin","Milnacipran","Amitriptyline"),
      Class      = c("SNRI","α2δ ligand","SNRI","TCA"),
      Mechanism  = c("SERT/NET","α2δ-1 Ca-ch","SERT/NET","SERT/NET/H1"),
      FDA_Status = c("Approved FM","Approved FM","Approved FM","Off-label"),
      stringsAsFactors = FALSE
    )
  })

  ## Reactive simulation
  sim_data <- eventReactive(input$run_btn, {
    run_sim(fm_mod,
            use_DUL=as.integer(input$use_DUL),
            use_PRE=as.integer(input$use_PRE),
            use_MIL=as.integer(input$use_MIL),
            use_TCA=as.integer(input$use_TCA),
            dose_DUL=input$dose_DUL, dose_PRE=input$dose_PRE,
            dose_MIL=input$dose_MIL, dose_TCA=input$dose_TCA,
            weeks=input$weeks, severity=input$severity)
  })

  ## PK plots
  output$pk_plot <- renderPlot({
    df <- sim_data()
    pk_long <- df %>%
      select(time_days, Cp_DUL_out, Cp_PRE_out, Cp_MIL_out, Cp_TCA_out) %>%
      pivot_longer(-time_days, names_to="Drug", values_to="Cp") %>%
      mutate(Drug = recode(Drug,
        Cp_DUL_out="Duloxetine", Cp_PRE_out="Pregabalin",
        Cp_MIL_out="Milnacipran", Cp_TCA_out="Amitriptyline"))
    ggplot(pk_long, aes(time_days, Cp, color=Drug)) +
      geom_line(linewidth=1) +
      labs(x="Time (days)", y="Plasma Conc (mg/L)", title="PK — Plasma Concentration Profiles") +
      theme_bw(13) + scale_color_brewer(palette="Set1")
  })

  output$pk_table <- renderTable({
    data.frame(
      Drug=c("Duloxetine","Pregabalin","Milnacipran","Amitriptyline"),
      Vd_L=c(1640,42,300,1500), CL_Lh=c(54,6.8,50,40),
      t_half_h=c(12,6,8,25), Bioavail=c("50%","90%","85%","48%"),
      stringsAsFactors=FALSE)
  })

  output$occ_plot <- renderPlot({
    df <- sim_data()
    occ <- df %>%
      select(time_days, inh_SERT_pct, inh_NET_pct, Ca_block_pct) %>%
      pivot_longer(-time_days, names_to="Target", values_to="Occupancy") %>%
      mutate(Target=recode(Target,
        inh_SERT_pct="SERT (5-HT)", inh_NET_pct="NET (NE)",
        Ca_block_pct="α2δ-1 (Ca-ch)"))
    ggplot(occ, aes(time_days, Occupancy, color=Target)) +
      geom_line(linewidth=1) +
      labs(x="Days", y="Target Occupancy (%)", title="Receptor Occupancy") +
      theme_bw(12) + ylim(0,100) + scale_color_brewer(palette="Set2")
  })

  ## PD plots
  output$ltp_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(time_days, LTP_cs)) +
      geom_line(color="#8E44AD", linewidth=1.2) +
      labs(x="Days", y="LTP index (0–1)", title="Spinal LTP / Central Sensitization") +
      theme_bw(12) + ylim(0,1)
  })

  output$dpms_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(time_days, DPMS)) +
      geom_line(color="#27AE60", linewidth=1.2) +
      labs(x="Days", y="DPMS (a.u.)", title="Descending Pain Modulation (↑=better)") +
      theme_bw(12)
  })

  output$mono_plot <- renderPlot({
    df <- sim_data()
    mono <- df %>% select(time_days, NE_syn, SHT_syn) %>%
      pivot_longer(-time_days, names_to="NT", values_to="Conc") %>%
      mutate(NT = recode(NT, NE_syn="Synaptic NE", SHT_syn="Synaptic 5-HT"))
    ggplot(mono, aes(time_days, Conc, color=NT)) +
      geom_line(linewidth=1.1) +
      scale_color_manual(values=c("Synaptic NE"="#E67E22","Synaptic 5-HT"="#9B59B6")) +
      labs(x="Days", y="Pool (a.u.)", title="CNS Monoamine Pools") +
      theme_bw(12)
  })

  output$wu_plot <- renderPlot({
    df <- sim_data()
    wu <- df %>% select(time_days, WindUp, NMDA_state) %>%
      pivot_longer(-time_days, names_to="Var", values_to="Val") %>%
      mutate(Var = recode(Var, WindUp="Wind-up", NMDA_state="NMDA Receptor State"))
    ggplot(wu, aes(time_days, Val, color=Var)) +
      geom_line(linewidth=1.1) +
      scale_color_manual(values=c("Wind-up"="#E74C3C","NMDA Receptor State"="#F39C12")) +
      labs(x="Days", y="Index (0–1)", title="Wind-up & NMDA Activation") +
      theme_bw(12)
  })

  ## Clinical endpoints
  output$pain_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(time_days, Pain_score)) +
      geom_line(color="#E74C3C", linewidth=1.2) +
      geom_hline(yintercept=c(3,5), linetype="dashed", color="grey60") +
      annotate("text", x=max(df$time_days)*0.8, y=3.3, label="Mild (3)", size=3) +
      annotate("text", x=max(df$time_days)*0.8, y=5.3, label="Moderate (5)", size=3) +
      labs(x="Days", y="NRS (0–10)", title="Pain Score") +
      theme_bw(13) + ylim(0,10)
  })

  output$fiq_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(time_days, FIQ_score)) +
      geom_line(color="#E67E22", linewidth=1.2) +
      labs(x="Days", y="FIQR (0–100)", title="Fibromyalgia Impact Questionnaire") +
      theme_bw(13)
  })

  output$fatigue_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(time_days, Fatigue_VAS)) +
      geom_line(color="#2ECC71", linewidth=1.2) +
      labs(x="Days", y="Fatigue VAS (0–100)", title="Fatigue (MFI proxy)") +
      theme_bw(13)
  })

  output$dep_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(time_days, Depression_score)) +
      geom_line(color="#3498DB", linewidth=1.2) +
      labs(x="Days", y="PHQ-9 (0–27)", title="Depression Score") +
      theme_bw(13)
  })

  output$outcome_table <- renderDT({
    df <- sim_data()
    wk <- df %>% filter(time_days >= (input$weeks - 0.5)) %>% slice_tail(n=4) %>%
      summarise(across(c(Pain_score,FIQ_score,Fatigue_VAS,Depression_score,SWS_depth,LTP_cs), mean)) %>%
      round(2)
    names(wk) <- c("Pain NRS","FIQR","Fatigue VAS","PHQ-9","SWS Depth","LTP Index")
    datatable(wk, rownames=FALSE, options=list(dom='t'))
  })

  ## Scenario comparison
  scenario_data <- eventReactive(input$run_all, {
    scenarios <- list(
      list(label="Untreated FM",          u=0,u_p=0,u_m=0,u_t=0),
      list(label="Duloxetine 60mg QD",    u=1,u_p=0,u_m=0,u_t=0),
      list(label="Pregabalin 150mg BID",  u=0,u_p=1,u_m=0,u_t=0),
      list(label="Milnacipran 50mg BID",  u=0,u_p=0,u_m=1,u_t=0),
      list(label="DUL + PRE Combo",       u=1,u_p=1,u_m=0,u_t=0),
      list(label="Amitriptyline 25mg QHS",u=0,u_p=0,u_m=0,u_t=1)
    )
    bind_rows(lapply(scenarios, function(s) {
      run_sim(fm_mod, use_DUL=s$u, use_PRE=s$u_p, use_MIL=s$u_m, use_TCA=s$u_t,
              weeks=12, severity=1.0) %>%
        mutate(scenario=s$label)
    }))
  })

  output$scenario_plot <- renderPlot({
    df <- scenario_data()
    cols <- c("Untreated FM"="#E74C3C","Duloxetine 60mg QD"="#3498DB",
              "Pregabalin 150mg BID"="#2ECC71","Milnacipran 50mg BID"="#9B59B6",
              "DUL + PRE Combo"="#E67E22","Amitriptyline 25mg QHS"="#1ABC9C")
    p1 <- ggplot(df, aes(time_days, Pain_score, color=scenario)) +
      geom_line(linewidth=0.9) + scale_color_manual(values=cols) +
      labs(x="Days",y="Pain NRS",title="Pain Score") + theme_bw(11) + ylim(0,10)
    p2 <- ggplot(df, aes(time_days, FIQ_score, color=scenario)) +
      geom_line(linewidth=0.9) + scale_color_manual(values=cols) +
      labs(x="Days",y="FIQR",title="FIQR") + theme_bw(11)
    p3 <- ggplot(df, aes(time_days, SWS_depth, color=scenario)) +
      geom_line(linewidth=0.9) + scale_color_manual(values=cols) +
      labs(x="Days",y="SWS Depth",title="Sleep Quality") + theme_bw(11)
    p4 <- ggplot(df, aes(time_days, LTP_cs, color=scenario)) +
      geom_line(linewidth=0.9) + scale_color_manual(values=cols) +
      labs(x="Days",y="LTP Index",title="Central Sensitization") + theme_bw(11)
    gridExtra::grid.arrange(p1, p2, p3, p4, nrow=2)
  })

  output$resp_table <- renderDT({
    df <- scenario_data()
    base_pain <- mean(filter(df, scenario=="Untreated FM", time_days<=1)$Pain_score)
    resp <- df %>%
      filter(time_days >= 83) %>%
      group_by(scenario) %>%
      summarise(
        `Pain NRS (wk12)` = round(mean(Pain_score),2),
        `% Change` = round(100*(mean(Pain_score)-base_pain)/base_pain,1),
        `≥30% Responder` = ifelse(mean(Pain_score) <= base_pain*0.70, "Yes","No"),
        `≥50% Responder` = ifelse(mean(Pain_score) <= base_pain*0.50, "Yes","No"),
        `FIQR (wk12)` = round(mean(FIQ_score),1),
        .groups="drop"
      ) %>% arrange(`Pain NRS (wk12)`)
    datatable(resp, rownames=FALSE)
  })

  ## Biomarker plots
  output$sp_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(time_days, SP_csf)) +
      geom_line(color="#E74C3C", linewidth=1.1) +
      labs(x="Days",y="SP (a.u.)",title="CSF Substance P") + theme_bw(12)
  })

  output$mg_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(time_days, MG_act)) +
      geom_line(color="#F39C12", linewidth=1.1) +
      labs(x="Days",y="Activation (a.u.)",title="Microglia Activation") + theme_bw(12)
  })

  output$cort_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(time_days, CORT)) +
      geom_line(color="#3498DB", linewidth=1.1) +
      labs(x="Days",y="Cortisol (norm.)",title="HPA Axis: Cortisol") + theme_bw(12)
  })

  output$sws_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(time_days, SWS_depth)) +
      geom_line(color="#9B59B6", linewidth=1.1) +
      labs(x="Days",y="SWS Depth (0–1)",title="Slow-Wave Sleep") + theme_bw(12) + ylim(0,1)
  })

  output$sns_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(time_days, SNS_tone)) +
      geom_line(color="#E67E22", linewidth=1.1) +
      labs(x="Days",y="SNS Tone (a.u.)",title="Sympathetic Nervous System Tone") + theme_bw(12)
  })

  output$il1b_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(time_days, IL1b_sp)) +
      geom_line(color="#C0392B", linewidth=1.1) +
      labs(x="Days",y="IL-1β (a.u.)",title="Spinal IL-1β (Neuroinflammation)") + theme_bw(12)
  })
}

shinyApp(ui, server)
