## ============================================================
## HCC QSP Shiny Dashboard
## Hepatocellular Carcinoma — Multi-drug Interactive Simulator
## Tabs: Patient Profile | Drug PK | Tumor Dynamics | PD/Biomarkers |
##       Scenario Comparison | Immune & Safety
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)

## ---- Inline mrgsolve model (compact version) ----------------------
hcc_code <- '
$PARAM
KA_Sora=0.25 F_Sora=0.38 CL_Sora=5.3 V2_Sora=7.5 Q_Sora=12 V3_Sora=99
KA_Lenva=1.20 F_Lenva=0.85 CL_Lenva=4.2 V2_Lenva=18 Q_Lenva=7.5 V3_Lenva=45
CL_Atez=0.200 V2_Atez=3.5 Q_Atez=0.41 V3_Atez=2.9
CL_Beva=0.170 V2_Beva=2.8 Q_Beva=0.19 V3_Beva=1.8
KA_Rego=0.20 F_Rego=0.70 CL_Rego=4.8 V2_Rego=6.0
lambda1=0.012 K_tumor=300 T0=1.0
Emax_Angio=0.75 IC50_Angio_S=2.5 IC50_Angio_L=0.08
Emax_MAPK=0.55 IC50_MAPK_S=3.5
Emax_VEGF=0.60 IC50_VEGF_B=15.0
Emax_ICB=0.70 IC50_ICB_A=10.0
k_kill=0.008 k_exhaust=0.002 k_Treg_sup=0.010
k_CD8_in=0.03 k_CD8_death=0.03
k_Treg_in=0.015 k_Treg_death=0.015
k_tumor_sup=0.005
k_Angio_in=0.02 k_Angio_deg=0.02 k_Angio_tum=0.01
VEGF_ss=1.0 k_VEGF=0.05 kd_VEGF=0.05
AFP_ss=100 k_AFP=0.05 k_AFP_cl=0.05
k_LF_decay=0.0005 k_LF_recover=0.001
use_sora=0 use_lenva=0 use_atez=0 use_beva=0 use_rego=0

$CMT GUT_S CENTRAL_S PERIPH_S MET_S
     GUT_L CENTRAL_L PERIPH_L
     CENTRAL_A PERIPH_A CENTRAL_B PERIPH_B
     GUT_R CENTRAL_R
     TUMOR CD8T TREG ANGIO VEGF_FREE AFP LF

$MAIN
double Cs=(CENTRAL_S>0)?CENTRAL_S/V2_Sora:0;
double Cl=(CENTRAL_L>0)?CENTRAL_L/V2_Lenva:0;
double Ca=(CENTRAL_A>0)?CENTRAL_A/V2_Atez:0;
double Cb=(CENTRAL_B>0)?CENTRAL_B/V2_Beva:0;
double Cr=(CENTRAL_R>0)?CENTRAL_R/V2_Rego:0;
double E_AS=Emax_Angio*Cs/(IC50_Angio_S+Cs)*use_sora;
double E_AL=Emax_Angio*Cl/(IC50_Angio_L+Cl)*use_lenva;
double E_VB=Emax_VEGF*Cb/(IC50_VEGF_B+Cb)*use_beva;
double E_Angio=1.0-(1.0-E_AS)*(1.0-E_AL)*(1.0-E_VB);
double E_MS=Emax_MAPK*Cs/(IC50_MAPK_S+Cs)*use_sora;
double E_MR=Emax_MAPK*Cr/(IC50_MAPK_S*0.8+Cr)*use_rego;
double E_MAPK=1.0-(1.0-E_MS)*(1.0-E_MR);
double E_ICB=Emax_ICB*Ca/(IC50_ICB_A+Ca)*use_atez;
double TGI=1.0-(1.0-E_Angio*0.6)*(1.0-E_MAPK*0.4);
F_GUT_S=F_Sora; F_GUT_L=F_Lenva; F_GUT_R=F_Rego;

$ODE
dxdt_GUT_S=   -KA_Sora*GUT_S;
dxdt_CENTRAL_S= KA_Sora*GUT_S-(CL_Sora+0.18+Q_Sora)/V2_Sora*CENTRAL_S+Q_Sora/V3_Sora*PERIPH_S;
dxdt_PERIPH_S=  Q_Sora/V2_Sora*CENTRAL_S-Q_Sora/V3_Sora*PERIPH_S;
dxdt_MET_S=    0.18/V2_Sora*CENTRAL_S-0.15*MET_S;
dxdt_GUT_L=   -KA_Lenva*GUT_L;
dxdt_CENTRAL_L= KA_Lenva*GUT_L-(CL_Lenva+Q_Lenva)/V2_Lenva*CENTRAL_L+Q_Lenva/V3_Lenva*PERIPH_L;
dxdt_PERIPH_L=  Q_Lenva/V2_Lenva*CENTRAL_L-Q_Lenva/V3_Lenva*PERIPH_L;
dxdt_CENTRAL_A=-(CL_Atez+Q_Atez)/V2_Atez*CENTRAL_A+Q_Atez/V3_Atez*PERIPH_A;
dxdt_PERIPH_A=  Q_Atez/V2_Atez*CENTRAL_A-Q_Atez/V3_Atez*PERIPH_A;
dxdt_CENTRAL_B=-(CL_Beva+Q_Beva)/V2_Beva*CENTRAL_B+Q_Beva/V3_Beva*PERIPH_B;
dxdt_PERIPH_B=  Q_Beva/V2_Beva*CENTRAL_B-Q_Beva/V3_Beva*PERIPH_B;
dxdt_GUT_R=   -KA_Rego*GUT_R;
dxdt_CENTRAL_R= KA_Rego*GUT_R-CL_Rego/V2_Rego*CENTRAL_R;
double T=(TUMOR>0)?TUMOR:0;
double lam_eff=lambda1*(1.0-T/K_tumor);
double kill_rate=k_kill*(1.0+E_ICB);
dxdt_TUMOR=lam_eff*T - TGI*lambda1*T - kill_rate*CD8T*T;
if(TUMOR<0.001) dxdt_TUMOR=0;
double CD8_infl=k_CD8_in*(1.0+0.5*T/(1.0+T));
double exhaust=k_exhaust*T*CD8T+k_Treg_sup*TREG*CD8T;
double reinvig=E_ICB*k_exhaust*T*CD8T;
dxdt_CD8T=CD8_infl-k_CD8_death*CD8T-exhaust+reinvig;
if(CD8T<0.001) dxdt_CD8T=0;
dxdt_TREG=k_Treg_in*(1.0+k_tumor_sup*T)-k_Treg_death*TREG;
if(TREG<0.001) dxdt_TREG=0;
dxdt_ANGIO=(k_Angio_in+k_Angio_tum*T)*(1.0-E_Angio)-k_Angio_deg*ANGIO;
if(ANGIO<0) dxdt_ANGIO=0;
dxdt_VEGF_FREE=k_VEGF*(1.0+0.5*T)*VEGF_ss-kd_VEGF*VEGF_FREE-E_VB*kd_VEGF*VEGF_FREE;
if(VEGF_FREE<0) dxdt_VEGF_FREE=0;
dxdt_AFP=k_AFP*T*AFP_ss-k_AFP_cl*AFP;
if(AFP<1) dxdt_AFP=0;
dxdt_LF=k_LF_recover*(1.0-LF)-k_LF_decay*T;

$TABLE
double CONC_Sora  =CENTRAL_S/V2_Sora;
double CONC_Lenva =CENTRAL_L/V2_Lenva;
double CONC_Atez  =CENTRAL_A/V2_Atez;
double CONC_Beva  =CENTRAL_B/V2_Beva;
double TGI_pct    =TGI*100;
double TumorRel   =TUMOR/T0;
double LF_pct     =LF*100;
double E_ICB_val  =E_ICB;
double pct_VEGFinh=E_VB*100;

$CAPTURE CONC_Sora CONC_Lenva CONC_Atez CONC_Beva
         TGI_pct TumorRel AFP LF_pct CD8T TREG ANGIO VEGF_FREE
         E_ICB_val pct_VEGFinh
'

mod_hcc <- mcode("HCC_Shiny", hcc_code, quiet = TRUE)

## ---- Helper: build events -------------------------------------------
build_events <- function(treatment, wt_kg, days) {
  evs <- list()
  if ("sorafenib" %in% treatment) {
    evs[[length(evs)+1]] <- ev(amt=400, ii=12, addl=days*2-1, cmt="GUT_S", time=0)
  }
  if ("lenvatinib" %in% treatment) {
    dose <- if (wt_kg >= 60) 12 else 8
    evs[[length(evs)+1]] <- ev(amt=dose, ii=24, addl=days-1, cmt="GUT_L", time=0)
  }
  if ("atez_beva" %in% treatment) {
    ncycles <- ceiling(days / 21)
    evs[[length(evs)+1]] <- ev(amt=1200,        ii=21*24, addl=ncycles-1, cmt="CENTRAL_A", time=0)
    evs[[length(evs)+1]] <- ev(amt=wt_kg*15,    ii=21*24, addl=ncycles-1, cmt="CENTRAL_B", time=0)
  }
  if ("regorafenib" %in% treatment) {
    ncycles <- ceiling(days / 28)
    for (i in seq_len(ncycles)) {
      evs[[length(evs)+1]] <- ev(amt=160, ii=24, addl=20, cmt="GUT_R",
                                  time=(i-1)*28*24)
    }
  }
  if (length(evs) == 0) return(ev(time=0, amt=0, cmt="GUT_S"))
  do.call(c, evs)
}

## ---- UI --------------------------------------------------------------
ui <- dashboardPage(
  skin = "black",
  dashboardHeader(title = "HCC QSP Simulator"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "profile",   icon = icon("user")),
      menuItem("Drug PK",             tabName = "pk",        icon = icon("pills")),
      menuItem("Tumor Dynamics",      tabName = "tumor",     icon = icon("chart-line")),
      menuItem("PD & Biomarkers",     tabName = "pd",        icon = icon("flask")),
      menuItem("Scenario Comparison", tabName = "compare",   icon = icon("balance-scale")),
      menuItem("Immune & Safety",     tabName = "immune",    icon = icon("shield-alt"))
    ),

    hr(),
    h5("  Patient Settings", style = "color:#aaa; margin-left:10px"),

    sliderInput("weight",   "Body Weight (kg)",  min=40, max=120, value=70, step=5),
    sliderInput("age",      "Age (years)",       min=20, max=85,  value=60, step=1),

    selectInput("bclc", "BCLC Stage",
      choices = c("B (Intermediate)" = "B", "C (Advanced)" = "C"),
      selected = "C"),

    selectInput("childpugh", "Child-Pugh Class",
      choices = c("A (5-6 pts)" = "A", "B (7-9 pts)" = "B"),
      selected = "A"),

    sliderInput("afp_base", "Baseline AFP (ng/mL)", min=10, max=5000, value=300, step=10),
    sliderInput("lf_base",  "Initial Liver Fn (%)", min=50, max=100, value=85, step=5),

    hr(),
    h5("  Treatment", style = "color:#aaa; margin-left:10px"),

    checkboxGroupInput("treatment", "Active Treatment(s):",
      choices = c("Sorafenib 400mg BID"        = "sorafenib",
                  "Lenvatinib 8/12mg QD"       = "lenvatinib",
                  "Atezo+Beva (q3w)"           = "atez_beva",
                  "Regorafenib 160mg QD (3/1)" = "regorafenib"),
      selected = "sorafenib"),

    sliderInput("sim_days", "Simulation Duration (days)", min=30, max=720, value=360, step=30),

    actionButton("run", "Run Simulation", icon=icon("play"),
                 class="btn-primary", style="margin:10px; width:90%")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #1a1a2e; }
      .box { background: #0d1117; border-top-color: #2ecc71; }
      .box-header { color: #2ecc71; font-weight: bold; }
      .nav-tabs-custom { background: #0d1117; }
    "))),

    tabItems(
      ## ---- TAB 1: Patient Profile -----------------------------------
      tabItem(tabName = "profile",
        fluidRow(
          valueBoxOutput("vbox_bclc",    width=3),
          valueBoxOutput("vbox_cp",      width=3),
          valueBoxOutput("vbox_afp",     width=3),
          valueBoxOutput("vbox_lf",      width=3)
        ),
        fluidRow(
          box(title="HCC Disease Staging Summary", width=6, solidHeader=TRUE,
              status="danger",
              plotlyOutput("plot_staging", height="350px")),
          box(title="Patient Baseline Parameters", width=6, solidHeader=TRUE,
              status="warning",
              DTOutput("tbl_params"))
        ),
        fluidRow(
          box(title="BCLC Treatment Algorithm", width=12, solidHeader=TRUE,
              status="primary",
              htmlOutput("bclc_algo"))
        )
      ),

      ## ---- TAB 2: Drug PK ------------------------------------------
      tabItem(tabName = "pk",
        fluidRow(
          box(title="PK Time-Concentration Profiles", width=8, solidHeader=TRUE,
              status="primary",
              plotlyOutput("plot_pk_all", height="450px")),
          box(title="PK Parameters", width=4, solidHeader=TRUE, status="info",
              DTOutput("tbl_pk"))
        ),
        fluidRow(
          box(title="Sorafenib Steady-State PK (Week 2-3)", width=6, solidHeader=TRUE,
              status="warning",
              plotlyOutput("plot_pk_sora_ss", height="300px")),
          box(title="Lenvatinib PK Profile", width=6, solidHeader=TRUE,
              status="warning",
              plotlyOutput("plot_pk_lenva", height="300px"))
        )
      ),

      ## ---- TAB 3: Tumor Dynamics -----------------------------------
      tabItem(tabName = "tumor",
        fluidRow(
          valueBoxOutput("vbox_tumor6mo",  width=3),
          valueBoxOutput("vbox_response",  width=3),
          valueBoxOutput("vbox_pfs_est",   width=3),
          valueBoxOutput("vbox_os_est",    width=3)
        ),
        fluidRow(
          box(title="Tumor Burden Over Time (Relative)", width=8, solidHeader=TRUE,
              status="danger",
              plotlyOutput("plot_tumor", height="400px")),
          box(title="Tumor Response Classification", width=4, solidHeader=TRUE,
              status="info",
              plotlyOutput("plot_waterfall", height="400px"))
        ),
        fluidRow(
          box(title="Tumor Growth Rate Analysis", width=12, solidHeader=TRUE,
              status="warning",
              plotlyOutput("plot_tgi", height="300px"))
        )
      ),

      ## ---- TAB 4: PD & Biomarkers ----------------------------------
      tabItem(tabName = "pd",
        fluidRow(
          box(title="AFP Biomarker Trajectory", width=6, solidHeader=TRUE,
              status="primary",
              plotlyOutput("plot_afp", height="350px")),
          box(title="Liver Function Reserve (%)", width=6, solidHeader=TRUE,
              status="warning",
              plotlyOutput("plot_lf", height="350px"))
        ),
        fluidRow(
          box(title="Angiogenesis State & VEGF-A", width=6, solidHeader=TRUE,
              status="danger",
              plotlyOutput("plot_vegf", height="300px")),
          box(title="TGI (%) per Drug Mechanism", width=6, solidHeader=TRUE,
              status="info",
              plotlyOutput("plot_mech", height="300px"))
        )
      ),

      ## ---- TAB 5: Scenario Comparison ------------------------------
      tabItem(tabName = "compare",
        fluidRow(
          box(title="All Scenarios: Relative Tumor Burden", width=12, solidHeader=TRUE,
              status="primary",
              plotlyOutput("plot_compare_tumor", height="400px"))
        ),
        fluidRow(
          box(title="All Scenarios: AFP", width=6, solidHeader=TRUE, status="warning",
              plotlyOutput("plot_compare_afp", height="300px")),
          box(title="All Scenarios: Liver Function", width=6, solidHeader=TRUE, status="danger",
              plotlyOutput("plot_compare_lf", height="300px"))
        ),
        fluidRow(
          box(title="Efficacy Summary Table (Month 6)", width=12, solidHeader=TRUE,
              status="success",
              DTOutput("tbl_compare"))
        )
      ),

      ## ---- TAB 6: Immune & Safety ----------------------------------
      tabItem(tabName = "immune",
        fluidRow(
          box(title="CD8+ T Cell & Treg Dynamics", width=6, solidHeader=TRUE,
              status="success",
              plotlyOutput("plot_immune", height="350px")),
          box(title="PD-L1 Blockade Effect (Atezolizumab)", width=6, solidHeader=TRUE,
              status="primary",
              plotlyOutput("plot_icb", height="350px"))
        ),
        fluidRow(
          box(title="Key Safety Signals", width=6, solidHeader=TRUE, status="danger",
              plotlyOutput("plot_safety", height="300px")),
          box(title="Adverse Event Profile by Drug", width=6, solidHeader=TRUE,
              status="warning",
              DTOutput("tbl_ae"))
        )
      )
    )
  )
)

## ---- Server ----------------------------------------------------------
server <- function(input, output, session) {

  ## ---- Reactive simulation -----------------------------------------
  sim_result <- eventReactive(input$run, {
    req(input$treatment)
    ev_all <- build_events(input$treatment, input$weight, input$sim_days)

    prm <- list(
      use_sora  = as.integer("sorafenib"    %in% input$treatment),
      use_lenva = as.integer("lenvatinib"   %in% input$treatment),
      use_atez  = as.integer("atez_beva"    %in% input$treatment),
      use_beva  = as.integer("atez_beva"    %in% input$treatment),
      use_rego  = as.integer("regorafenib"  %in% input$treatment),
      AFP_ss    = input$afp_base
    )

    ini <- list(TUMOR=1, CD8T=1, TREG=1, ANGIO=1, VEGF_FREE=1,
                AFP=input$afp_base, LF=input$lf_base/100)

    out <- mod_hcc %>%
      param(prm) %>%
      init(ini) %>%
      mrgsim(ev_all, end=input$sim_days, delta=1, hmax=0.05)

    df <- as.data.frame(out) %>%
      mutate(time_d = time, time_mo = time / 30,
             TumorVol_rel = TumorRel,
             AFP_ng = AFP,
             LF_pct_v = LF_pct)
    df
  }, ignoreNULL = FALSE)

  ## ---- All-scenario comparison ------------------------------------
  sim_all_scenarios <- reactive({
    scenarios <- list(
      list(name="Sorafenib 400mg BID",        trt=c("sorafenib")),
      list(name="Lenvatinib 8/12mg QD",        trt=c("lenvatinib")),
      list(name="Atezo+Beva q3w (IMbrave150)", trt=c("atez_beva")),
      list(name="Regorafenib 160mg QD (3/1)",  trt=c("regorafenib")),
      list(name="BSC (No Treatment)",           trt=character(0))
    )
    lapply(scenarios, function(s) {
      ev_all <- build_events(s$trt, input$weight, 360)
      prm <- list(
        use_sora  = as.integer("sorafenib"   %in% s$trt),
        use_lenva = as.integer("lenvatinib"  %in% s$trt),
        use_atez  = as.integer("atez_beva"   %in% s$trt),
        use_beva  = as.integer("atez_beva"   %in% s$trt),
        use_rego  = as.integer("regorafenib" %in% s$trt),
        AFP_ss    = input$afp_base
      )
      ini <- list(TUMOR=1,CD8T=1,TREG=1,ANGIO=1,VEGF_FREE=1,
                  AFP=input$afp_base, LF=input$lf_base/100)
      out <- mod_hcc %>% param(prm) %>% init(ini) %>%
        mrgsim(ev_all, end=360, delta=1, hmax=0.05)
      as.data.frame(out) %>% mutate(Scenario=s$name, time_mo=time/30)
    }) %>% bind_rows()
  })

  ## ---- Value boxes --------------------------------------------------
  output$vbox_bclc <- renderValueBox({
    valueBox(input$bclc, "BCLC Stage", icon=icon("sitemap"), color="red")
  })
  output$vbox_cp <- renderValueBox({
    valueBox(input$childpugh, "Child-Pugh", icon=icon("liver"), color="yellow")
  })
  output$vbox_afp <- renderValueBox({
    valueBox(paste0(input$afp_base, " ng/mL"), "Baseline AFP",
             icon=icon("vial"), color="blue")
  })
  output$vbox_lf <- renderValueBox({
    valueBox(paste0(input$lf_base, "%"), "Liver Function",
             icon=icon("heartbeat"), color="green")
  })

  output$vbox_tumor6mo <- renderValueBox({
    df <- sim_result()
    val <- df %>% filter(abs(time_d - 180) < 1) %>% pull(TumorRel) %>% first()
    pct <- if (!is.na(val)) round((val - 1) * 100, 1) else NA
    col <- if (!is.na(pct) && pct < -30) "green" else if (!is.na(pct) && pct > 25) "red" else "yellow"
    valueBox(paste0(pct, "%"), "Tumor Change (6mo)", icon=icon("chart-bar"), color=col)
  })

  output$vbox_response <- renderValueBox({
    df <- sim_result()
    val <- df %>% filter(abs(time_d - 180) < 1) %>% pull(TumorRel) %>% first()
    resp <- if (is.na(val)) "N/A" else if (val < 0.70) "PR" else if (val > 1.25) "PD" else "SD"
    col  <- if (resp == "PR") "green" else if (resp == "PD") "red" else "yellow"
    valueBox(resp, "Best Response (6mo)", icon=icon("stethoscope"), color=col)
  })

  output$vbox_pfs_est <- renderValueBox({
    df <- sim_result()
    pd_time <- df %>% filter(TumorRel > 1.25) %>% pull(time_d) %>% first()
    pfs_mo  <- if (is.na(pd_time)) paste0(">", round(input$sim_days/30)) else round(pd_time/30, 1)
    valueBox(paste0(pfs_mo, " mo"), "Est. PFS", icon=icon("clock"), color="blue")
  })

  output$vbox_os_est <- renderValueBox({
    df <- sim_result()
    lf_fail <- df %>% filter(LF_pct < 30) %>% pull(time_d) %>% first()
    os_mo <- if (is.na(lf_fail)) paste0(">", round(input$sim_days/30)) else round(lf_fail/30, 1)
    valueBox(paste0(os_mo, " mo"), "Est. OS (liver-fn limit)",
             icon=icon("heartbeat"), color="purple")
  })

  ## ---- Tab 1: Staging radar -----------------------------------------
  output$plot_staging <- renderPlotly({
    df_radar <- data.frame(
      parameter = c("Tumor Burden", "Liver Function", "VEGF Level",
                    "Immune Activity", "AFP Level"),
      value = c(0.8, input$lf_base/100, 0.7, 0.4,
                pmin(input$afp_base / 1000, 1))
    )
    plot_ly(df_radar, type="scatterpolar", r=~value, theta=~parameter,
            fill="toself", fillcolor="rgba(231,76,60,0.3)",
            line=list(color="#e74c3c")) %>%
      layout(polar=list(radialaxis=list(range=c(0,1))),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"))
  })

  output$tbl_params <- renderDT({
    df <- data.frame(
      Parameter = c("Weight (kg)", "Age (years)", "BCLC Stage",
                    "Child-Pugh", "Baseline AFP (ng/mL)", "Liver Fn (%)",
                    "Treatment(s)"),
      Value = c(input$weight, input$age, input$bclc,
                input$childpugh, input$afp_base, input$lf_base,
                paste(input$treatment, collapse=" + "))
    )
    datatable(df, options=list(dom="t", pageLength=10), rownames=FALSE)
  })

  output$bclc_algo <- renderUI({
    HTML("<div style='color:#aaa; padding:10px;'>
    <b style='color:#2ecc71'>BCLC Staging & Treatment Algorithm:</b><br/>
    <b>Stage 0 (Very Early)</b>: Single nodule &lt;2cm; Child-Pugh A → Resection / Ablation (OS ~5y)<br/>
    <b>Stage A (Early)</b>: Single or ≤3 nodules ≤3cm; Child-Pugh A-B → Resection / Transplant / Ablation (OS ~3y)<br/>
    <b>Stage B (Intermediate)</b>: Multinodular, no macrovascular invasion; Child-Pugh A-B → <b>TACE</b> (OS ~2.5y)<br/>
    <b>Stage C (Advanced)</b>: Macrovascular invasion or extrahepatic spread →
    <b>Atezolizumab+Bevacizumab</b> (1st-line, OS ~19mo) or
    <b>Sorafenib/Lenvatinib</b> (1st-line, OS ~13mo)<br/>
    <b>Stage D (Terminal)</b>: Child-Pugh C → Best Supportive Care (OS &lt;3mo)<br/>
    <i style='color:#f39c12'>2nd-line: Regorafenib, Cabozantinib, Ramucirumab (AFP≥400); Pembrolizumab, Nivolumab+Ipilimumab</i>
    </div>")
  })

  ## ---- Tab 2: Drug PK -----------------------------------------------
  output$plot_pk_all <- renderPlotly({
    df <- sim_result()
    pk_long <- df %>%
      select(time_mo, CONC_Sora, CONC_Lenva, CONC_Atez, CONC_Beva) %>%
      pivot_longer(-time_mo, names_to="Drug", values_to="Conc") %>%
      filter(Conc > 0.001) %>%
      mutate(Drug = recode(Drug,
        CONC_Sora="Sorafenib (μg/mL)", CONC_Lenva="Lenvatinib (μg/mL)",
        CONC_Atez="Atezolizumab (mg/L)", CONC_Beva="Bevacizumab (mg/L)"))
    plot_ly(pk_long, x=~time_mo, y=~Conc, color=~Drug, type="scatter",
            mode="lines", line=list(width=2)) %>%
      layout(xaxis=list(title="Time (months)", color="white"),
             yaxis=list(title="Concentration", color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$tbl_pk <- renderDT({
    df <- data.frame(
      Drug = c("Sorafenib","Lenvatinib","Atezolizumab","Bevacizumab","Regorafenib"),
      Route = c("Oral","Oral","IV","IV","Oral"),
      Dose = c("400mg BID","8/12mg QD","1200mg q3w","15mg/kg q3w","160mg QD"),
      `T½ (h)` = c("~27","~28","~27d","~20d","~28"),
      `F (%)` = c(38,85,100,100,70),
      check.names=FALSE
    )
    datatable(df, options=list(dom="t"), rownames=FALSE)
  })

  output$plot_pk_sora_ss <- renderPlotly({
    df <- sim_result() %>%
      filter(time_d >= 120, time_d <= 192)  # day 5-8 steady-state window
    plot_ly(df, x=~(time_d-120), y=~CONC_Sora, type="scatter", mode="lines",
            line=list(color="#e74c3c", width=2), name="Sorafenib") %>%
      layout(xaxis=list(title="Day (relative to day 120)", color="white"),
             yaxis=list(title="Conc (μg/mL)", color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"))
  })

  output$plot_pk_lenva <- renderPlotly({
    df <- sim_result() %>% filter(time_d <= 60)
    plot_ly(df, x=~time_d, y=~CONC_Lenva, type="scatter", mode="lines",
            line=list(color="#f39c12", width=2), name="Lenvatinib") %>%
      layout(xaxis=list(title="Time (days)", color="white"),
             yaxis=list(title="Conc (μg/mL)", color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"))
  })

  ## ---- Tab 3: Tumor -------------------------------------------------
  output$plot_tumor <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_mo, y=~TumorRel, type="scatter", mode="lines",
            line=list(color="#e74c3c", width=2.5), name="Tumor") %>%
      add_segments(x=0, xend=max(df$time_mo), y=0.7, yend=0.7,
                   line=list(dash="dash", color="#2ecc71"), name="PR threshold") %>%
      add_segments(x=0, xend=max(df$time_mo), y=1.25, yend=1.25,
                   line=list(dash="dot", color="#f39c12"), name="PD threshold") %>%
      layout(xaxis=list(title="Time (months)", color="white"),
             yaxis=list(title="Relative Tumor Burden", color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$plot_waterfall <- renderPlotly({
    df <- sim_result()
    resp_pts <- df %>%
      filter(time_d %in% c(28,56,84,112,140,168,196,224,252,280)) %>%
      mutate(pct_change = (TumorRel - 1) * 100,
             timepoint  = paste0("Day ", round(time_d)))
    plot_ly(resp_pts, x=~timepoint, y=~pct_change,
            type="bar",
            marker=list(color=ifelse(resp_pts$pct_change < -30, "#2ecc71",
                                     ifelse(resp_pts$pct_change > 25, "#e74c3c", "#f39c12")))) %>%
      add_segments(x=-0.5, xend=nrow(resp_pts)-0.5, y=-30, yend=-30,
                   line=list(dash="dash", color="#2ecc71"), name="PR") %>%
      add_segments(x=-0.5, xend=nrow(resp_pts)-0.5, y=25, yend=25,
                   line=list(dash="dash", color="#e74c3c"), name="PD") %>%
      layout(xaxis=list(title="Timepoint", color="white"),
             yaxis=list(title="Tumor Change (%)", color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"))
  })

  output$plot_tgi <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_mo, y=~TGI_pct, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(46,204,113,0.2)",
            line=list(color="#2ecc71", width=2)) %>%
      layout(xaxis=list(title="Time (months)", color="white"),
             yaxis=list(title="Tumor Growth Inhibition (%)", range=c(0,100), color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"))
  })

  ## ---- Tab 4: PD & Biomarkers ---------------------------------------
  output$plot_afp <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_mo, y=~AFP_ng, type="scatter", mode="lines",
            line=list(color="#9b59b6", width=2)) %>%
      add_segments(x=0, xend=max(df$time_mo), y=400, yend=400,
                   line=list(dash="dash", color="#e74c3c"), name="AFP 400 (ramucirumab threshold)") %>%
      layout(xaxis=list(title="Time (months)", color="white"),
             yaxis=list(title="AFP (ng/mL)", type="log", color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"))
  })

  output$plot_lf <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_mo, y=~LF_pct_v, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(52,152,219,0.2)",
            line=list(color="#3498db", width=2)) %>%
      add_segments(x=0, xend=max(df$time_mo), y=70, yend=70,
                   line=list(dash="dash", color="#f39c12"), name="Child-Pugh B") %>%
      add_segments(x=0, xend=max(df$time_mo), y=50, yend=50,
                   line=list(dash="dot", color="#e74c3c"), name="Child-Pugh C") %>%
      layout(xaxis=list(title="Time (months)", color="white"),
             yaxis=list(title="Liver Function (%)", range=c(0,105), color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"))
  })

  output$plot_vegf <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_mo, y=~VEGF_FREE, type="scatter", mode="lines",
            name="Free VEGF-A", line=list(color="#e74c3c", width=2)) %>%
      add_trace(y=~ANGIO, name="Angiogenesis Drive",
                line=list(color="#f39c12", width=2, dash="dash")) %>%
      layout(xaxis=list(title="Time (months)", color="white"),
             yaxis=list(title="Relative Level", color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$plot_mech <- renderPlotly({
    df <- sim_result() %>%
      mutate(AngioPct = pct_VEGFinh,
             MAPKPct  = pmin(TGI_pct, 55),
             ICBPct   = E_ICB_val * 70)
    mech_long <- df %>%
      select(time_mo, AngioPct, MAPKPct, ICBPct) %>%
      pivot_longer(-time_mo, names_to="Mechanism", values_to="Effect_pct")
    plot_ly(mech_long, x=~time_mo, y=~Effect_pct, color=~Mechanism,
            type="scatter", mode="lines") %>%
      layout(xaxis=list(title="Time (months)", color="white"),
             yaxis=list(title="Drug Effect (%)", range=c(0,100), color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  ## ---- Tab 5: Scenario Comparison -----------------------------------
  cols5 <- c("#e74c3c","#f39c12","#2ecc71","#9b59b6","#95a5a6")
  names(cols5) <- c("Sorafenib 400mg BID","Lenvatinib 8/12mg QD",
                    "Atezo+Beva q3w (IMbrave150)","Regorafenib 160mg QD (3/1)",
                    "BSC (No Treatment)")

  output$plot_compare_tumor <- renderPlotly({
    df <- sim_all_scenarios()
    plot_ly(df, x=~time_mo, y=~TumorRel, color=~Scenario,
            type="scatter", mode="lines", colors=cols5) %>%
      layout(xaxis=list(title="Time (months)", color="white"),
             yaxis=list(title="Relative Tumor Burden", color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$plot_compare_afp <- renderPlotly({
    df <- sim_all_scenarios()
    plot_ly(df, x=~time_mo, y=~AFP, color=~Scenario,
            type="scatter", mode="lines", colors=cols5) %>%
      layout(xaxis=list(title="Time (months)", color="white"),
             yaxis=list(title="AFP (ng/mL)", type="log", color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$plot_compare_lf <- renderPlotly({
    df <- sim_all_scenarios()
    plot_ly(df, x=~time_mo, y=~LF_pct, color=~Scenario,
            type="scatter", mode="lines", colors=cols5) %>%
      layout(xaxis=list(title="Time (months)", color="white"),
             yaxis=list(title="Liver Function (%)", color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$tbl_compare <- renderDT({
    df <- sim_all_scenarios()
    sum_tbl <- df %>%
      filter(abs(time_d - 180) < 1 | time == 180) %>%
      group_by(Scenario) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(
        `Tumor (6mo rel)` = round(TumorRel, 3),
        `AFP (6mo ng/mL)` = round(AFP, 0),
        `Liver Fn (6mo%)` = round(LF_pct, 1),
        `CD8 Level`       = round(CD8T, 3),
        Response          = case_when(
          TumorRel < 0.70  ~ "PR",
          TumorRel > 1.25  ~ "PD",
          TRUE             ~ "SD"
        )
      ) %>%
      select(Scenario, `Tumor (6mo rel)`, Response, `AFP (6mo ng/mL)`,
             `Liver Fn (6mo%)`, `CD8 Level`)
    datatable(sum_tbl, options=list(dom="t", pageLength=6), rownames=FALSE) %>%
      formatStyle("Response",
        backgroundColor = styleEqual(c("PR","SD","PD"), c("#27ae60","#f39c12","#c0392b")))
  })

  ## ---- Tab 6: Immune & Safety ---------------------------------------
  output$plot_immune <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_mo, y=~CD8T, name="CD8+ T cells",
            type="scatter", mode="lines", line=list(color="#2ecc71", width=2)) %>%
      add_trace(y=~TREG, name="Tregs",
                line=list(color="#e74c3c", width=2, dash="dash")) %>%
      layout(xaxis=list(title="Time (months)", color="white"),
             yaxis=list(title="Relative Cell Level", color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$plot_icb <- renderPlotly({
    df <- sim_result()
    plot_ly(df, x=~time_mo, y=~E_ICB_val, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(155,89,182,0.2)",
            line=list(color="#9b59b6", width=2)) %>%
      layout(xaxis=list(title="Time (months)", color="white"),
             yaxis=list(title="PD-L1 Blockade Effect (0–1)", range=c(0,1), color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"))
  })

  output$plot_safety <- renderPlotly({
    ae_df <- data.frame(
      AE = c("HFSR", "Hypertension", "Diarrhea", "Fatigue",
             "Hypothyroidism", "Proteinuria", "irAE (any)", "Hepatotoxicity"),
      Sorafenib   = c(30, 10, 43, 37, 7, 2, 0, 4),
      Lenvatinib  = c(27, 42, 39, 44, 57, 25, 0, 5),
      Atez_Beva   = c(3, 30, 25, 35, 12, 30, 35, 6),
      Regorafenib = c(53, 31, 40, 47, 8, 8, 0, 10)
    )
    plot_ly(ae_df, x=~AE, y=~Sorafenib, name="Sorafenib", type="bar",
            marker=list(color="#e74c3c")) %>%
      add_trace(y=~Lenvatinib,  name="Lenvatinib",  marker=list(color="#f39c12")) %>%
      add_trace(y=~Atez_Beva,   name="Atezo+Beva",  marker=list(color="#2ecc71")) %>%
      add_trace(y=~Regorafenib, name="Regorafenib", marker=list(color="#9b59b6")) %>%
      layout(barmode="group",
             xaxis=list(title="Adverse Event", color="white"),
             yaxis=list(title="Incidence (any grade, %)", color="white"),
             paper_bgcolor="#0d1117", plot_bgcolor="#0d1117",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$tbl_ae <- renderDT({
    ae_df <- data.frame(
      `Adverse Event` = c("HFSR","Hypertension","Diarrhea","Fatigue",
                          "Hypothyroidism","Proteinuria","irAE (any gr.)","Hepatotoxicity"),
      `Sorafenib (%)` = c(30,10,43,37,7,2,0,4),
      `Lenvatinib (%)` = c(27,42,39,44,57,25,0,5),
      `Atezo+Beva (%)` = c(3,30,25,35,12,30,35,6),
      `Regorafenib (%)` = c(53,31,40,47,8,8,0,10),
      check.names=FALSE
    )
    datatable(ae_df, options=list(dom="t", pageLength=10), rownames=FALSE)
  })
}

shinyApp(ui, server)
