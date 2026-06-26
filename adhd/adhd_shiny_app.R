## ============================================================
## ADHD QSP Shiny Dashboard
## Tabs: Patient Profile · PK · DA/NE · Executive Function
##       Clinical Endpoints · Scenario Comparison · Biomarkers
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(plotly)

## ---- Embed mrgsolve model ---------------------------------------------------
adhd_code <- '
$PARAM
ka1=1.2, CL1=31.5, V2c=448.0, Q1=75.6, V2p=560.0, F1=0.22
ka2=1.0, CL2=39.2, V3c=245.0, Q2=42.0, V3p=350.0, F2=0.75
ka3=2.0, CL3=24.5, V4c=59.5,  Q3=28.0, V4p=168.0, F3=0.63
ka4=0.6, CL4=3.5,  V5c=196.0, Q4=5.6,  V5p=1120.0,F4=0.80
ka5=0.5, CL5=19.6, V6c=105.0, F5=0.88
DA_base=1.0, krel_DA=0.5, kel_DA=2.0
NE_base=0.8, krel_NE=0.4, kel_NE=1.8
Ki_MPH_DAT=0.05, Ki_MPH_NET=0.34
Ki_AMP_DAT=0.10, Ki_AMP_NET=0.04
Ki_ATX_NET=0.002, Ki_GFN_A2A=0.001, Ki_VLX_NET=0.042
Emax_DA_WM=0.8, EC50_DA_WM=1.5, Emax_NE_WM=0.9, EC50_NE_WM=1.2, hill=1.5
DA_opt=2.0, DA_bw=1.5, NE_opt=1.5, NE_bw=1.2
ADHD_RS_base=32.0, k_symp=0.1, k_WM_effect=0.6
WT=70, AGE=10, CYP2D6=1

$CMT GUT1 CENT1 PER1 GUT2 CENT2 PER2 GUT3 CENT3 PER3
     GUT4 CENT4 PER4 GUT5 CENT5
     DA_syn NE_syn DAT_occ NET_occ PFC_DA PFC_NE WM_idx ExecFun
     ADHD_RS CGI_S QoL_idx

$INIT GUT1=0,CENT1=0,PER1=0, GUT2=0,CENT2=0,PER2=0,
      GUT3=0,CENT3=0,PER3=0, GUT4=0,CENT4=0,PER4=0,
      GUT5=0,CENT5=0,
      DA_syn=1.0, NE_syn=0.8, DAT_occ=0, NET_occ=0,
      PFC_DA=0.5, PFC_NE=0.5, WM_idx=0.35, ExecFun=0.30,
      ADHD_RS=32.0, CGI_S=4.0, QoL_idx=0.40

$ODE
double C1=CENT1/(V2c*WT/70.0), C2=CENT2/(V3c*WT/70.0);
double C3=CENT3/(V4c*WT/70.0), C4=CENT4/(V5c*WT/70.0);
double C5=CENT5/(V6c*WT/70.0);
double CL1a=CL1*pow(WT/70.0,0.75), CL2a=CL2*pow(WT/70.0,0.75);
double CL3a=CL3*CYP2D6*pow(WT/70.0,0.75), CL4a=CL4*pow(WT/70.0,0.75);
double CL5a=CL5*pow(WT/70.0,0.75);
dxdt_GUT1=-ka1*GUT1;
dxdt_CENT1=ka1*GUT1-(CL1a/V2c)*CENT1-(Q1/V2c)*CENT1+(Q1/V2p)*PER1;
dxdt_PER1=(Q1/V2c)*CENT1-(Q1/V2p)*PER1;
dxdt_GUT2=-ka2*GUT2;
dxdt_CENT2=ka2*GUT2-(CL2a/V3c)*CENT2-(Q2/V3c)*CENT2+(Q2/V3p)*PER2;
dxdt_PER2=(Q2/V3c)*CENT2-(Q2/V3p)*PER2;
dxdt_GUT3=-ka3*GUT3;
dxdt_CENT3=ka3*GUT3-(CL3a/V4c)*CENT3-(Q3/V4c)*CENT3+(Q3/V4p)*PER3;
dxdt_PER3=(Q3/V4c)*CENT3-(Q3/V4p)*PER3;
dxdt_GUT4=-ka4*GUT4;
dxdt_CENT4=ka4*GUT4-(CL4a/V5c)*CENT4-(Q4/V5c)*CENT4+(Q4/V5p)*PER4;
dxdt_PER4=(Q4/V5c)*CENT4-(Q4/V5p)*PER4;
dxdt_GUT5=-ka5*GUT5;
dxdt_CENT5=ka5*GUT5-(CL5a/V6c)*CENT5;
double MPH_uM=C1*1000.0/233.7, AMP_uM=C2*1000.0/135.2;
double ATX_uM=C3*1000.0/291.8, GFN_uM=C4*1000.0/246.7, VLX_uM=C5*1000.0/237.3;
double DAT_ss=(MPH_uM/Ki_MPH_DAT+AMP_uM/Ki_AMP_DAT)/
              (1.0+MPH_uM/Ki_MPH_DAT+AMP_uM/Ki_AMP_DAT);
double NET_ss=(MPH_uM/Ki_MPH_NET+AMP_uM/Ki_AMP_NET+ATX_uM/Ki_ATX_NET+VLX_uM/Ki_VLX_NET)/
              (1.0+MPH_uM/Ki_MPH_NET+AMP_uM/Ki_AMP_NET+ATX_uM/Ki_ATX_NET+VLX_uM/Ki_VLX_NET);
dxdt_DAT_occ=2.0*(DAT_ss-DAT_occ);
dxdt_NET_occ=2.0*(NET_ss-NET_occ);
double AEF=1.0+2.0*AMP_uM/(AMP_uM+Ki_AMP_DAT);
double GFN_occ=GFN_uM/(GFN_uM+Ki_GFN_A2A);
dxdt_DA_syn=krel_DA*AEF*DA_base-kel_DA*(1.0-0.85*DAT_occ)*DA_syn;
dxdt_NE_syn=krel_NE*(1.0-0.6*GFN_occ)*NE_base-kel_NE*(1.0-0.85*NET_occ)*NE_syn;
double PDA=exp(-pow(DA_syn-DA_opt,2.0)/(2.0*DA_bw*DA_bw));
double PNE=fmin(exp(-pow(NE_syn-NE_opt,2.0)/(2.0*NE_bw*NE_bw))*(1.0+0.5*GFN_occ),1.0);
dxdt_PFC_DA=0.5*(PDA-PFC_DA);
dxdt_PFC_NE=0.5*(PNE-PFC_NE);
double WMss=0.35+0.35*PFC_DA+0.30*PFC_NE;
double EFss=0.30+0.30*PFC_DA+0.40*PFC_NE;
dxdt_WM_idx=0.3*(WMss-WM_idx);
dxdt_ExecFun=0.3*(EFss-ExecFun);
double sr=(WM_idx-0.35)/0.65*k_WM_effect+(ExecFun-0.30)/0.70*(1.0-k_WM_effect);
dxdt_ADHD_RS=k_symp*(ADHD_RS_base*(1.0-sr)-ADHD_RS);
double cgi_ss=4.0-2.5*sr; if(cgi_ss<1.0) cgi_ss=1.0;
dxdt_CGI_S=0.05*(cgi_ss-CGI_S);
double qol_ss=0.40+0.55*sr; if(qol_ss>1.0) qol_ss=1.0;
dxdt_QoL_idx=0.05*(qol_ss-QoL_idx);

$TABLE
double Cp_MPH=CENT1/(V2c*WT/70.0), Cp_AMP=CENT2/(V3c*WT/70.0);
double Cp_ATX=CENT3/(V4c*WT/70.0), Cp_GFN=CENT4/(V5c*WT/70.0);
double Cp_VLX=CENT5/(V6c*WT/70.0);
double Cp_MPH_nM=Cp_MPH*1000.0/233.7, Cp_AMP_nM=Cp_AMP*1000.0/135.2;
double Cp_ATX_nM=Cp_ATX*1000.0/291.8, Cp_GFN_nM=Cp_GFN*1000.0/246.7;
double Cp_VLX_nM=Cp_VLX*1000.0/237.3;
double response_pct=100.0*(1.0-ADHD_RS/ADHD_RS_base);

$CAPTURE Cp_MPH Cp_AMP Cp_ATX Cp_GFN Cp_VLX
         Cp_MPH_nM Cp_AMP_nM Cp_ATX_nM Cp_GFN_nM Cp_VLX_nM
         DA_syn NE_syn DAT_occ NET_occ PFC_DA PFC_NE WM_idx ExecFun
         ADHD_RS CGI_S QoL_idx response_pct
'

mod <- mcode("adhd_shiny", adhd_code, quiet = TRUE)

## ---- UI ---------------------------------------------------------------------
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "ADHD QSP Dashboard", titleWidth = 280),
  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("① Patient Profile",    tabName = "patient",   icon = icon("user")),
      menuItem("② Drug PK",            tabName = "pk",        icon = icon("pills")),
      menuItem("③ DA/NE Dynamics",     tabName = "dane",      icon = icon("brain")),
      menuItem("④ PFC & Cognition",    tabName = "pfc",       icon = icon("sitemap")),
      menuItem("⑤ Clinical Endpoints", tabName = "endpoints", icon = icon("chart-line")),
      menuItem("⑥ Scenario Compare",   tabName = "scenarios", icon = icon("layer-group")),
      menuItem("⑦ Biomarker Panel",    tabName = "biomarker", icon = icon("microscope"))
    ),
    hr(),
    h4("  Patient Parameters", style = "color:#ccc; margin-left:10px"),
    sliderInput("wt",   "Body Weight (kg)", 20, 100, 70, step = 5),
    sliderInput("age",  "Age (years)",       6, 65,  10, step = 1),
    selectInput("cyp",  "CYP2D6 Phenotype",
                choices = c("EM (normal)"=1, "IM (intermediate)"=0.5, "PM (poor)"=0.1),
                selected = 1),
    selectInput("subtype", "ADHD Subtype",
                choices = c("Combined"="combined","Inattentive"="inattentive",
                            "Hyperactive/Impulsive"="hyperactive")),
    sliderInput("baseline_rs", "Baseline ADHD-RS (0–54)", 20, 54, 32, step = 1),
    hr(),
    h4("  Simulation Duration", style = "color:#ccc; margin-left:10px"),
    sliderInput("weeks", "Weeks to simulate", 4, 26, 12, step = 2)
  ),
  dashboardBody(
    tabItems(

      ## ── Tab 1: Patient Profile ────────────────────────────────────────────
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "Patient Demographics & Diagnosis", width = 6, solidHeader = TRUE, status = "primary",
            fluidRow(
              column(6,
                h4("Patient Summary"),
                tableOutput("pt_summary")
              ),
              column(6,
                h4("ADHD Diagnosis Criteria"),
                tableOutput("dx_criteria")
              )
            )
          ),
          box(title = "Comorbidities & Risk Factors", width = 6, solidHeader = TRUE, status = "warning",
            checkboxGroupInput("comorbid", "Select Comorbidities:",
              choices = c("Anxiety Disorder"="anxiety",
                          "Major Depression"="mdd",
                          "Oppositional Defiant Disorder (ODD)"="odd",
                          "Conduct Disorder"="cd",
                          "Learning Disability"="ld",
                          "Sleep Disturbance"="sleep",
                          "Tic Disorder / Tourette"="tic",
                          "Autism Spectrum"="asd"),
              selected = c("anxiety")),
            hr(),
            h5("Genetic Risk Profile:"),
            checkboxGroupInput("genetics", NULL,
              choices = c("DAT1 10-repeat"="dat1",
                          "DRD4 7-repeat"="drd4",
                          "COMT Val/Val"="comt",
                          "SNAP25 variant"="snap25"),
              selected = c("drd4"))
          )
        ),
        fluidRow(
          box(title = "ADHD-RS Symptom Profile", width = 12, solidHeader = TRUE, status = "info",
            plotOutput("pt_radar", height = "350px")
          )
        )
      ),

      ## ── Tab 2: Drug PK ───────────────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Drug Selection & Dose", width = 4, solidHeader = TRUE, status = "primary",
            radioButtons("drug_sel", "Select Drug:",
              choices = c(
                "Methylphenidate IR" = "mph_ir",
                "Methylphenidate ER" = "mph_er",
                "Amphetamine XR"     = "amp_xr",
                "Atomoxetine"        = "atx",
                "Guanfacine ER"      = "gfn",
                "Viloxazine ER"      = "vlx"
              ), selected = "mph_ir"),
            uiOutput("dose_slider"),
            hr(),
            h5("PK Parameters (estimated):"),
            tableOutput("pk_params_tbl")
          ),
          box(title = "Plasma Concentration-Time Profile (24 h)", width = 8, solidHeader = TRUE, status = "info",
            plotlyOutput("pk_plot_24h", height = "350px")
          )
        ),
        fluidRow(
          box(title = "DAT / NET Occupancy vs Time", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("occ_plot", height = "300px")
          ),
          box(title = "Steady-State PK Summary", width = 6, solidHeader = TRUE, status = "success",
            tableOutput("pk_ss_tbl")
          )
        )
      ),

      ## ── Tab 3: DA/NE Dynamics ────────────────────────────────────────────
      tabItem(tabName = "dane",
        fluidRow(
          box(title = "Drug & Dose", width = 3, solidHeader = TRUE, status = "primary",
            selectInput("dane_drug", "Drug:",
              choices = c("Untreated","MPH IR","MPH ER","AMP XR","ATX","GFN ER","VLX ER"),
              selected = "MPH IR"),
            uiOutput("dane_dose_ui"),
            sliderInput("dane_days", "Days to simulate:", 1, 28, 7)
          ),
          box(title = "Synaptic DA & NE Levels", width = 9, solidHeader = TRUE, status = "info",
            plotlyOutput("da_ne_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Inverted-U: DA Tone vs PFC Function", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("invU_plot", height = "300px")
          ),
          box(title = "DA / NE Neurotransmitter Panel", width = 6, solidHeader = TRUE, status = "success",
            plotlyOutput("nt_panel", height = "300px")
          )
        )
      ),

      ## ── Tab 4: PFC & Cognition ────────────────────────────────────────────
      tabItem(tabName = "pfc",
        fluidRow(
          box(title = "Drug Selection", width = 3, solidHeader = TRUE, status = "primary",
            selectInput("pfc_drug", "Drug:",
              choices = c("Untreated","MPH IR","MPH ER","AMP XR","ATX","GFN ER","VLX ER"),
              selected = "MPH IR"),
            uiOutput("pfc_dose_ui")
          ),
          box(title = "PFC Tone & Cognitive Index Over Time", width = 9, solidHeader = TRUE, status = "info",
            plotlyOutput("pfc_plot", height = "380px")
          )
        ),
        fluidRow(
          box(title = "Working Memory Index", width = 4, solidHeader = TRUE, status = "success",
            plotlyOutput("wm_gauge", height = "220px")
          ),
          box(title = "Executive Function Index", width = 4, solidHeader = TRUE, status = "success",
            plotlyOutput("ef_gauge", height = "220px")
          ),
          box(title = "Quality of Life Index", width = 4, solidHeader = TRUE, status = "success",
            plotlyOutput("qol_gauge", height = "220px")
          )
        )
      ),

      ## ── Tab 5: Clinical Endpoints ─────────────────────────────────────────
      tabItem(tabName = "endpoints",
        fluidRow(
          box(title = "Drug & Dose", width = 3, solidHeader = TRUE, status = "primary",
            selectInput("ep_drug", "Drug:",
              choices = c("Untreated","MPH IR","MPH ER","AMP XR","ATX","GFN ER","VLX ER"),
              selected = "MPH IR"),
            uiOutput("ep_dose_ui"),
            hr(),
            h5("Clinical Trial Benchmarks:"),
            tableOutput("ep_bench_tbl")
          ),
          box(title = "ADHD-RS-5 Score Over Time", width = 9, solidHeader = TRUE, status = "info",
            plotlyOutput("adhd_rs_plot", height = "380px")
          )
        ),
        fluidRow(
          box(title = "CGI-Severity", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("cgi_plot", height = "280px")
          ),
          box(title = "Response & Remission Rates", width = 6, solidHeader = TRUE, status = "success",
            h4("At Week 12:"),
            tableOutput("resp_tbl")
          )
        )
      ),

      ## ── Tab 6: Scenario Comparison ────────────────────────────────────────
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "Scenario Setup", width = 3, solidHeader = TRUE, status = "primary",
            checkboxGroupInput("sc_drugs", "Include in comparison:",
              choices = c("Untreated","MPH IR 10mg TID","MPH ER 36mg QD",
                          "AMP XR 20mg QD","ATX 80mg QD",
                          "GFN ER 4mg QD","VLX ER 400mg QD"),
              selected = c("Untreated","MPH IR 10mg TID","AMP XR 20mg QD","ATX 80mg QD"))
          ),
          box(title = "ADHD-RS-5 Trajectory — All Scenarios", width = 9, solidHeader = TRUE, status = "info",
            plotlyOutput("sc_adhd_rs", height = "380px")
          )
        ),
        fluidRow(
          box(title = "12-Week Summary Table", width = 12, solidHeader = TRUE, status = "success",
            DTOutput("sc_summary_dt")
          )
        )
      ),

      ## ── Tab 7: Biomarker Panel ────────────────────────────────────────────
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "Drug Selection", width = 3, solidHeader = TRUE, status = "primary",
            selectInput("bm_drug", "Drug:",
              choices = c("Untreated","MPH IR","MPH ER","AMP XR","ATX","GFN ER","VLX ER"),
              selected = "AMP XR"),
            uiOutput("bm_dose_ui"),
            hr(),
            h5("Biomarker Reference Ranges:"),
            tableOutput("bm_ref_tbl")
          ),
          box(title = "Multi-Biomarker Panel", width = 9, solidHeader = TRUE, status = "info",
            plotlyOutput("bm_panel_plot", height = "450px")
          )
        ),
        fluidRow(
          box(title = "DAT Occupancy & DA Levels", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("bm_dat_plot", height = "280px")
          ),
          box(title = "NET Occupancy & NE Levels", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("bm_net_plot", height = "280px")
          )
        )
      )
    )
  )
)

## ---- Server -----------------------------------------------------------------
server <- function(input, output, session) {

  ## ─── Helper: build event table ────────────────────────────────────────────
  build_ev <- function(drug, dose_mg, weeks) {
    n_days <- weeks * 7
    switch(drug,
      "Untreated"     = ev(amt=0, cmt="GUT1", time=0),
      "MPH IR"        = ev(amt=dose_mg, cmt="GUT1", ii=8,  addl=n_days*3-1, time=0),
      "MPH ER"        = ev(amt=dose_mg, cmt="GUT1", ii=24, addl=n_days-1,   time=0),
      "MPH IR 10mg TID" = ev(amt=10, cmt="GUT1", ii=8,  addl=n_days*3-1, time=0),
      "MPH ER 36mg QD"  = ev(amt=36, cmt="GUT1", ii=24, addl=n_days-1,   time=0),
      "AMP XR"        = ev(amt=dose_mg, cmt="GUT2", ii=24, addl=n_days-1, time=0),
      "AMP XR 20mg QD"= ev(amt=20, cmt="GUT2", ii=24, addl=n_days-1, time=0),
      "ATX"           = ev(amt=dose_mg, cmt="GUT3", ii=24, addl=n_days-1, time=0),
      "ATX 80mg QD"   = ev(amt=80, cmt="GUT3", ii=24, addl=n_days-1, time=0),
      "GFN ER"        = ev(amt=dose_mg, cmt="GUT4", ii=24, addl=n_days-1, time=0),
      "GFN ER 4mg QD" = ev(amt=4,  cmt="GUT4", ii=24, addl=n_days-1, time=0),
      "VLX ER"        = ev(amt=dose_mg, cmt="GUT5", ii=24, addl=n_days-1, time=0),
      "VLX ER 400mg QD"=ev(amt=400,cmt="GUT5", ii=24, addl=n_days-1, time=0)
    )
  }

  run_sim <- function(drug, dose, weeks, ka1_val = 1.2) {
    ev_obj <- build_ev(drug, dose, weeks)
    mrgsim(
      mod %>% param(WT = input$wt, AGE = input$age,
                    CYP2D6 = as.numeric(input$cyp),
                    ADHD_RS_base = input$baseline_rs,
                    ka1 = ka1_val),
      ev = ev_obj, delta = 0.5, end = weeks * 168
    ) %>% as_tibble()
  }

  ## ─── Tab 1: Patient profile ───────────────────────────────────────────────
  output$pt_summary <- renderTable({
    data.frame(
      Parameter = c("Age","Weight","CYP2D6","ADHD Subtype","Baseline ADHD-RS"),
      Value = c(paste(input$age, "yr"), paste(input$wt, "kg"),
                c("1"="EM","0.5"="IM","0.1"="PM")[input$cyp],
                input$subtype, input$baseline_rs)
    )
  })

  output$dx_criteria <- renderTable({
    data.frame(
      Criterion = c("Symptoms ≥6 inattentive", "Symptoms ≥6 hyperactive",
                    "Duration ≥6 months", "Present before age 12",
                    "≥2 settings", "DSM-5 diagnosis"),
      Met = c("Yes","Yes","Yes","Yes","Yes","Yes")
    )
  })

  output$pt_radar <- renderPlot({
    df <- data.frame(
      Domain = c("Inattention","Hyperactivity","Impulsivity",
                 "Working Memory","Executive Function","Emotional Regulation"),
      Severity = c(8, 6, 7, 4, 4, 6)
    )
    ggplot(df, aes(x = Domain, y = Severity)) +
      geom_col(fill = "#1565C0", alpha = 0.7) +
      coord_polar() +
      ylim(0, 9) +
      labs(title = "ADHD Symptom Severity Profile") +
      theme_minimal()
  })

  ## ─── Tab 2: Drug PK ──────────────────────────────────────────────────────
  output$dose_slider <- renderUI({
    limits <- list(
      mph_ir = c(5, 60, 10), mph_er = c(18, 72, 36),
      amp_xr = c(5, 40, 20), atx = c(10, 100, 80),
      gfn = c(1, 7, 4), vlx = c(100, 600, 400)
    )
    lims <- limits[[input$drug_sel]]
    sliderInput("pk_dose", "Dose (mg):", lims[1], lims[2], lims[3])
  })

  pk_data <- reactive({
    req(input$pk_dose)
    drug_map <- c(mph_ir="MPH IR", mph_er="MPH ER", amp_xr="AMP XR",
                  atx="ATX", gfn="GFN ER", vlx="VLX ER")
    ka1v <- ifelse(input$drug_sel == "mph_er", 0.4, 1.2)
    run_sim(drug_map[[input$drug_sel]], input$pk_dose, 1, ka1v) %>%
      filter(time <= 30)
  })

  cp_col <- reactive({
    c(mph_ir="Cp_MPH_nM", mph_er="Cp_MPH_nM", amp_xr="Cp_AMP_nM",
      atx="Cp_ATX_nM", gfn="Cp_GFN_nM", vlx="Cp_VLX_nM")[input$drug_sel]
  })

  output$pk_plot_24h <- renderPlotly({
    df <- pk_data(); col <- cp_col()
    p <- ggplot(df, aes(x = time, y = .data[[col]])) +
      geom_line(color = "#1565C0", linewidth = 1.2) +
      labs(x = "Time (h)", y = "Plasma Concentration (nM)",
           title = paste(input$drug_sel, "PK")) +
      theme_bw()
    ggplotly(p)
  })

  output$occ_plot <- renderPlotly({
    df <- pk_data()
    p <- df %>%
      select(time, DAT_occ, NET_occ) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x = time, y = value*100, color = name)) +
      geom_line(linewidth = 1) +
      labs(x = "Time (h)", y = "Occupancy (%)", color = "Target") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_params_tbl <- renderTable({
    drug_params <- list(
      mph_ir = data.frame(Param=c("t½","Tmax","F%","ki_DAT"),
                          Value=c("2.5h","1.5h","22%","34nM")),
      mph_er = data.frame(Param=c("t½","Tmax","F%","ki_DAT"),
                          Value=c("8h","6–8h","22%","34nM")),
      amp_xr = data.frame(Param=c("t½","Tmax","F%","ki_DAT"),
                          Value=c("9–14h","7h","75%","100nM")),
      atx    = data.frame(Param=c("t½","Tmax","F%","ki_NET"),
                          Value=c("5h(EM)","1–2h","63%","2nM")),
      gfn    = data.frame(Param=c("t½","Tmax","F%","ki_α2A"),
                          Value=c("17h","5h","80%","1nM")),
      vlx    = data.frame(Param=c("t½","Tmax","F%","ki_NET"),
                          Value=c("7h","5h","88%","42nM"))
    )
    drug_params[[input$drug_sel]]
  })

  output$pk_ss_tbl <- renderTable({
    data.frame(
      Metric = c("Cmax (nM)", "AUC24 (nM·h)", "DAT_occ_max (%)", "NET_occ_max (%)"),
      Value  = c("—", "—", "—", "—")
    )
  })

  ## ─── Tab 3: DA/NE Dynamics ───────────────────────────────────────────────
  output$dane_dose_ui <- renderUI({
    req(input$dane_drug)
    drug_defaults <- c("Untreated"=0,"MPH IR"=10,"MPH ER"=36,
                       "AMP XR"=20,"ATX"=80,"GFN ER"=4,"VLX ER"=400)
    sliderInput("dane_dose", "Dose (mg):", 0, 600,
                value = drug_defaults[[input$dane_drug]])
  })

  dane_sim <- reactive({
    req(input$dane_dose)
    run_sim(input$dane_drug, input$dane_dose, max(1, input$dane_days/7)) %>%
      filter(time <= input$dane_days * 24)
  })

  output$da_ne_plot <- renderPlotly({
    df <- dane_sim()
    p <- df %>%
      select(time, DA_syn, NE_syn) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x = time/24, y = value, color = name)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = c(1.0, 0.8), linetype = "dashed", alpha = 0.5) +
      labs(x = "Day", y = "Concentration (nM)", color = "Transmitter") +
      theme_bw()
    ggplotly(p)
  })

  output$invU_plot <- renderPlotly({
    da_seq <- seq(0, 5, by = 0.05)
    inv_u  <- data.frame(DA = da_seq,
                          PFC = exp(-(da_seq - 2.0)^2 / (2 * 1.5^2)))
    current_da <- if (nrow(dane_sim()) > 0) tail(dane_sim()$DA_syn, 1) else 1.0
    p <- ggplot(inv_u, aes(x = DA, y = PFC)) +
      geom_line(color = "steelblue", linewidth = 1.5) +
      geom_vline(xintercept = 1.0, color = "red", linetype = "dashed") +
      geom_vline(xintercept = current_da, color = "green", linewidth = 1) +
      annotate("text", x = 1.0, y = 0.9, label = "ADHD\nbaseline", color="red", size=3) +
      annotate("text", x = current_da + 0.1, y = 0.85,
               label = paste0("Drug\n[", round(current_da,2)," nM]"),
               color="green", size=3) +
      labs(x = "[DA] PFC (nM)", y = "PFC Function Index",
           title = "Inverted-U (Arnsten)") +
      theme_bw()
    ggplotly(p)
  })

  output$nt_panel <- renderPlotly({
    df <- dane_sim()
    p <- df %>%
      select(time, DAT_occ, NET_occ) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x = time/24, y = value*100, color = name)) +
      geom_line(linewidth = 1) +
      labs(x = "Day", y = "Occupancy (%)", color = "Transporter") +
      theme_bw()
    ggplotly(p)
  })

  ## ─── Tab 4: PFC & Cognition ──────────────────────────────────────────────
  output$pfc_dose_ui <- renderUI({
    defaults <- c("Untreated"=0,"MPH IR"=10,"MPH ER"=36,
                  "AMP XR"=20,"ATX"=80,"GFN ER"=4,"VLX ER"=400)
    sliderInput("pfc_dose", "Dose (mg):", 0, 600, defaults[[input$pfc_drug]])
  })

  pfc_sim <- reactive({
    req(input$pfc_dose)
    run_sim(input$pfc_drug, input$pfc_dose, input$weeks)
  })

  output$pfc_plot <- renderPlotly({
    df <- pfc_sim()
    p <- df %>%
      select(time, PFC_DA, PFC_NE, WM_idx, ExecFun) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x = time/168, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Week", y = "Index (0–1)", color = "Variable") +
      theme_bw()
    ggplotly(p)
  })

  make_gauge <- function(val, title) {
    plot_ly(type = "indicator", mode = "gauge+number",
      value = round(val * 100),
      title = list(text = title),
      gauge = list(
        axis = list(range = list(0, 100)),
        bar  = list(color = "#1565C0"),
        steps = list(
          list(range = c(0,40),  color = "#EF9A9A"),
          list(range = c(40,65), color = "#FFF9C4"),
          list(range = c(65,100),color = "#C8E6C9")
        )
      )
    ) %>% layout(margin = list(t=50))
  }

  output$wm_gauge  <- renderPlotly({ make_gauge(tail(pfc_sim()$WM_idx,1), "Working Memory (%)") })
  output$ef_gauge  <- renderPlotly({ make_gauge(tail(pfc_sim()$ExecFun,1), "Exec. Function (%)") })
  output$qol_gauge <- renderPlotly({ make_gauge(tail(pfc_sim()$QoL_idx,1), "Quality of Life (%)") })

  ## ─── Tab 5: Clinical Endpoints ───────────────────────────────────────────
  output$ep_dose_ui <- renderUI({
    defaults <- c("Untreated"=0,"MPH IR"=10,"MPH ER"=36,
                  "AMP XR"=20,"ATX"=80,"GFN ER"=4,"VLX ER"=400)
    sliderInput("ep_dose","Dose (mg):",0,600,defaults[[input$ep_drug]])
  })

  ep_sim <- reactive({
    req(input$ep_dose)
    run_sim(input$ep_drug, input$ep_dose, input$weeks)
  })

  output$adhd_rs_plot <- renderPlotly({
    df <- ep_sim()
    p <- ggplot(df, aes(x = time/168, y = ADHD_RS)) +
      geom_line(color="#1565C0", linewidth=1.3) +
      geom_hline(yintercept = input$baseline_rs * 0.5,
                 linetype = "dashed", color = "green") +
      annotate("text", x = 2, y = input$baseline_rs * 0.5 + 0.5,
               label = "50% reduction threshold", color = "green", size = 3) +
      labs(x = "Week", y = "ADHD-RS-5 Total Score",
           title = paste("ADHD-RS-5:", input$ep_drug)) +
      scale_y_reverse(limits = c(input$baseline_rs + 2, 5)) +
      theme_bw()
    ggplotly(p)
  })

  output$cgi_plot <- renderPlotly({
    df <- ep_sim()
    p <- ggplot(df, aes(x = time/168, y = CGI_S)) +
      geom_line(color="#E65100", linewidth=1.3) +
      labs(x = "Week", y = "CGI-Severity (1–7)", title = "CGI-S") +
      scale_y_reverse(limits = c(5.5, 1)) +
      theme_bw()
    ggplotly(p)
  })

  output$ep_bench_tbl <- renderTable({
    data.frame(
      Drug = c("MPH IR","AMP XR","ATX","GFN ER","VLX ER"),
      RS_reduction = c("-10","-12","-8","-7","-8"),
      Response_pct = c("60%","65%","50%","45%","47%"),
      Source = c("MTA","Biederman02","Michelson01","Sallee09","Nasser21")
    )
  })

  output$resp_tbl <- renderTable({
    df <- ep_sim() %>% filter(time == max(time))
    data.frame(
      Metric = c("ADHD-RS reduction","% Response (≥25%)","CGI-S final","QoL final"),
      Value  = c(round(input$baseline_rs - df$ADHD_RS, 1),
                 ifelse((input$baseline_rs - df$ADHD_RS) / input$baseline_rs >= 0.25, "Yes","No"),
                 round(df$CGI_S, 1),
                 round(df$QoL_idx, 2))
    )
  })

  ## ─── Tab 6: Scenario Comparison ──────────────────────────────────────────
  sc_sim_all <- reactive({
    weeks <- input$weeks
    scenarios <- list(
      "Untreated"       = list(ev=ev(amt=0,cmt="GUT1",time=0)),
      "MPH IR 10mg TID" = list(ev=ev(amt=10,cmt="GUT1",ii=8, addl=weeks*7*3-1,time=0)),
      "MPH ER 36mg QD"  = list(ev=ev(amt=36,cmt="GUT1",ii=24,addl=weeks*7-1,time=0)),
      "AMP XR 20mg QD"  = list(ev=ev(amt=20,cmt="GUT2",ii=24,addl=weeks*7-1,time=0)),
      "ATX 80mg QD"     = list(ev=ev(amt=80,cmt="GUT3",ii=24,addl=weeks*7-1,time=0)),
      "GFN ER 4mg QD"   = list(ev=ev(amt=4, cmt="GUT4",ii=24,addl=weeks*7-1,time=0)),
      "VLX ER 400mg QD" = list(ev=ev(amt=400,cmt="GUT5",ii=24,addl=weeks*7-1,time=0))
    )
    sel <- input$sc_drugs
    bind_rows(lapply(sel, function(s) {
      if (!s %in% names(scenarios)) return(NULL)
      mrgsim(mod %>% param(WT=input$wt, AGE=input$age, CYP2D6=as.numeric(input$cyp),
                            ADHD_RS_base=input$baseline_rs),
             ev = scenarios[[s]]$ev, delta=24, end=weeks*168) %>%
        as_tibble() %>% mutate(scenario=s, week=time/168)
    }))
  })

  output$sc_adhd_rs <- renderPlotly({
    df <- sc_sim_all()
    req(nrow(df) > 0)
    p <- ggplot(df, aes(x=week, y=ADHD_RS, color=scenario)) +
      geom_line(linewidth=1.2) +
      scale_y_reverse() +
      labs(x="Week", y="ADHD-RS-5", color="Treatment") +
      theme_bw()
    ggplotly(p)
  })

  output$sc_summary_dt <- renderDT({
    df <- sc_sim_all()
    req(nrow(df) > 0)
    df %>%
      filter(near(week, input$weeks, tol=0.5)) %>%
      group_by(scenario) %>%
      summarise(
        `ADHD-RS Final`    = round(mean(ADHD_RS), 1),
        `Reduction`        = round(mean(input$baseline_rs - ADHD_RS), 1),
        `Response (%)`     = round(mean(response_pct), 1),
        `CGI-S Final`      = round(mean(CGI_S), 1),
        `WM Index`         = round(mean(WM_idx), 2),
        `Exec Fun Index`   = round(mean(ExecFun), 2),
        `QoL Index`        = round(mean(QoL_idx), 2)
      ) %>%
      arrange(`ADHD-RS Final`) %>%
      datatable(options = list(pageLength = 10), rownames = FALSE)
  })

  ## ─── Tab 7: Biomarker Panel ──────────────────────────────────────────────
  output$bm_dose_ui <- renderUI({
    defaults <- c("Untreated"=0,"MPH IR"=10,"MPH ER"=36,
                  "AMP XR"=20,"ATX"=80,"GFN ER"=4,"VLX ER"=400)
    sliderInput("bm_dose","Dose (mg):",0,600,defaults[[input$bm_drug]])
  })

  bm_sim <- reactive({
    req(input$bm_dose)
    run_sim(input$bm_drug, input$bm_dose, input$weeks)
  })

  output$bm_panel_plot <- renderPlotly({
    df <- bm_sim()
    p <- df %>%
      select(time, DAT_occ, NET_occ, DA_syn, NE_syn, PFC_DA, PFC_NE, WM_idx, ExecFun) %>%
      pivot_longer(-time, names_to = "biomarker", values_to = "value") %>%
      ggplot(aes(x = time/168, y = value, color = biomarker)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~biomarker, scales = "free_y", ncol = 3) +
      labs(x = "Week", y = "Value") +
      theme_bw() + theme(legend.position = "none")
    ggplotly(p)
  })

  output$bm_dat_plot <- renderPlotly({
    df <- bm_sim()
    p <- df %>% select(time, DAT_occ, DA_syn) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x=time/168, y=value, color=name)) +
      geom_line(linewidth=1) +
      labs(x="Week", y="Value", color="Biomarker") + theme_bw()
    ggplotly(p)
  })

  output$bm_net_plot <- renderPlotly({
    df <- bm_sim()
    p <- df %>% select(time, NET_occ, NE_syn) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x=time/168, y=value, color=name)) +
      geom_line(linewidth=1) +
      labs(x="Week", y="Value", color="Biomarker") + theme_bw()
    ggplotly(p)
  })

  output$bm_ref_tbl <- renderTable({
    data.frame(
      Biomarker = c("DAT occupancy","NET occupancy","Synaptic DA","Synaptic NE"),
      ADHD = c("~0%","~0%","1.0 nM","0.8 nM"),
      Therapeutic = c("50–80%","70–90%","1.5–2.5 nM","1.2–2.0 nM")
    )
  })
}

## ---- Launch -----------------------------------------------------------------
shinyApp(ui, server)
