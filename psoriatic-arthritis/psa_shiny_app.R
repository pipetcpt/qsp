## ============================================================
## Psoriatic Arthritis (PsA) – Interactive Shiny QSP Dashboard
## 7 Tabs: Patient Profile · PK · Cytokines · Skin (PASI) ·
##         Joint (DAPSA) · Scenario Comparison · Biomarkers
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

## ---- Embed mrgsolve model code --------------------------
psa_code <- '
$PARAM
  BWT=85 BSA=2.0
  ka_ADA=0.0065 CL_ADA=0.0117 V1_ADA=2.75 V2_ADA=1.76 Q_ADA=0.0096 F_ADA=0.64
  kon_ADA=0.096 koff_ADA=0.00048 kdeg_ADA=0.0042
  ka_IXE=0.0055 CL_IXE=0.0119 V1_IXE=3.8 F_IXE=0.60
  kon_IXE=0.21 koff_IXE=0.00021 kdeg_IXE=0.0052
  ka_GUS=0.0072 CL_GUS=0.0082 V1_GUS=3.2 F_GUS=0.49
  kon_GUS=0.30 koff_GUS=0.00018 kdeg_GUS=0.0038
  ka_UPA=0.72 CL_UPA=38.4 V1_UPA=220 F_UPA=0.76
  ka_APR=0.65 CL_APR=10.1 V1_APR=84 F_APR=0.73
  kin_IL17=0.018 kout_IL17=0.014 IL17_base=1.286
  kin_TNF=0.022 kout_TNF=0.017 TNF_base=1.294
  kin_IL23=0.015 kout_IL23=0.012 IL23_base=1.25
  kin_Th17=0.008 kout_Th17=0.006 Th17_base=1.0
  kin_RANKL=0.010 kout_RANKL=0.008 RANKL_base=1.25
  kin_PASI=0.004 kout_PASI=0.003 PASI_base=18.0
  kin_DAPSA=0.005 kout_DAPSA=0.0038 DAPSA_base=28.0
  kin_CRP=0.030 kout_CRP=0.025 CRP_base=1.2
  kin_S100=0.012 kout_S100=0.009 S100_base=1.3
  kin_MTSS=0.00015 kout_MTSS=0.00005
  Emax_ADA=0.92 EC50_ADA=0.18
  Emax_IXE=0.96 EC50_IXE=0.12
  Emax_GUS=0.94 EC50_GUS=0.08
  Emax_UPA=0.88 EC50_UPA=0.082
  Emax_APR=0.78 EC50_APR=0.54
  n_hill=1.5
  fb_IL17_TNF=0.35 fb_TNF_IL17=0.28 fb_IL23_Th17=0.55 fb_Th17_IL17=0.70

$CMT
  ADA_DEPOT ADA_C1 ADA_C2 ADA_RC
  IXE_DEPOT IXE_C1 IXE_RC
  GUS_DEPOT GUS_C1 GUS_RC
  UPA_GI UPA_C1
  APR_GI APR_C1
  IL17 TNFa IL23 TH17 RANKL CRP_ratio S100 PASI DAPSA MTSS

$INIT
  ADA_DEPOT=0 ADA_C1=0 ADA_C2=0 ADA_RC=0
  IXE_DEPOT=0 IXE_C1=0 IXE_RC=0
  GUS_DEPOT=0 GUS_C1=0 GUS_RC=0
  UPA_GI=0 UPA_C1=0
  APR_GI=0 APR_C1=0
  IL17=1.286 TNFa=1.294 IL23=1.25 TH17=1.0 RANKL=1.25
  CRP_ratio=1.2 S100=1.3 PASI=18.0 DAPSA=28.0 MTSS=0.0

$MAIN
  double ADA_conc = ADA_C1/V1_ADA;
  double IXE_conc = IXE_C1/V1_IXE;
  double GUS_conc = GUS_C1/V1_GUS;
  double UPA_conc = (UPA_C1/V1_UPA)*1000.0;
  double APR_conc = (APR_C1/V1_APR)*1000.0;

  double Imax_ADA = Emax_ADA*pow(ADA_conc,n_hill)/(pow(EC50_ADA,n_hill)+pow(ADA_conc,n_hill));
  double Imax_IXE = Emax_IXE*pow(IXE_conc,n_hill)/(pow(EC50_IXE,n_hill)+pow(IXE_conc,n_hill));
  double Imax_GUS = Emax_GUS*pow(GUS_conc,n_hill)/(pow(EC50_GUS,n_hill)+pow(GUS_conc,n_hill));
  double Imax_UPA = Emax_UPA*pow(UPA_conc,n_hill)/(pow(EC50_UPA,n_hill)+pow(UPA_conc,n_hill));
  double Imax_APR = Emax_APR*pow(APR_conc,n_hill)/(pow(EC50_APR,n_hill)+pow(APR_conc,n_hill));

  double Itot_IL17  = 1.0-(1.0-Imax_IXE)*(1.0-Imax_UPA*0.75)*(1.0-Imax_APR*0.45);
  double Itot_TNF   = 1.0-(1.0-Imax_ADA)*(1.0-Imax_UPA*0.60)*(1.0-Imax_APR*0.55);
  double Itot_IL23  = 1.0-(1.0-Imax_GUS)*(1.0-Imax_UPA*0.80);
  double Itot_Th17  = Itot_IL23*0.6 + Itot_IL17*0.3;

  double TNF_drive  = 1.0 + 0.35*(TNFa/TNF_base-1.0);
  double IL23_drive = 1.0 + 0.55*(IL23/IL23_base-1.0);
  double Th17_drive = 1.0 + 0.70*(TH17/Th17_base-1.0);

$ODE
  dxdt_ADA_DEPOT = -ka_ADA*F_ADA*ADA_DEPOT;
  dxdt_ADA_C1    =  ka_ADA*F_ADA*ADA_DEPOT-(CL_ADA/V1_ADA)*ADA_C1-(Q_ADA/V1_ADA)*ADA_C1+(Q_ADA/V2_ADA)*ADA_C2-kon_ADA*ADA_conc*TNFa+koff_ADA*ADA_RC;
  dxdt_ADA_C2    =  (Q_ADA/V1_ADA)*ADA_C1-(Q_ADA/V2_ADA)*ADA_C2;
  dxdt_ADA_RC    =  kon_ADA*ADA_conc*TNFa-koff_ADA*ADA_RC-kdeg_ADA*ADA_RC;
  dxdt_IXE_DEPOT = -ka_IXE*F_IXE*IXE_DEPOT;
  dxdt_IXE_C1    =  ka_IXE*F_IXE*IXE_DEPOT-(CL_IXE/V1_IXE)*IXE_C1-kon_IXE*IXE_conc*IL17+koff_IXE*IXE_RC;
  dxdt_IXE_RC    =  kon_IXE*IXE_conc*IL17-koff_IXE*IXE_RC-kdeg_IXE*IXE_RC;
  dxdt_GUS_DEPOT = -ka_GUS*F_GUS*GUS_DEPOT;
  dxdt_GUS_C1    =  ka_GUS*F_GUS*GUS_DEPOT-(CL_GUS/V1_GUS)*GUS_C1-kon_GUS*GUS_conc*IL23+koff_GUS*GUS_RC;
  dxdt_GUS_RC    =  kon_GUS*GUS_conc*IL23-koff_GUS*GUS_RC-kdeg_GUS*GUS_RC;
  dxdt_UPA_GI    = -ka_UPA*F_UPA*UPA_GI;
  dxdt_UPA_C1    =  ka_UPA*F_UPA*UPA_GI-(CL_UPA/V1_UPA)*UPA_C1;
  dxdt_APR_GI    = -ka_APR*F_APR*APR_GI;
  dxdt_APR_C1    =  ka_APR*F_APR*APR_GI-(CL_APR/V1_APR)*APR_C1;
  dxdt_IL17  = kin_IL17*TNF_drive*IL23_drive*Th17_drive*(1.0-Itot_IL17)-kout_IL17*IL17;
  dxdt_TNFa  = kin_TNF*(1.0+0.35*(IL17/IL17_base-1.0))*(1.0-Itot_TNF)-kout_TNF*TNFa;
  dxdt_IL23  = kin_IL23*(1.0-Itot_IL23)-kout_IL23*IL23;
  dxdt_TH17  = kin_Th17*IL23_drive*(1.0-Itot_Th17)-kout_Th17*TH17;
  dxdt_RANKL = kin_RANKL*(IL17/IL17_base)*(TNFa/TNF_base)*(1.0-0.5*Itot_IL17-0.4*Itot_TNF)-kout_RANKL*RANKL;
  dxdt_CRP_ratio = kin_CRP*((IL17/IL17_base+TNFa/TNF_base+IL23/IL23_base)/3.0)*(1.0-0.4*Itot_IL17-0.4*Itot_TNF-0.2*Itot_IL23)-kout_CRP*CRP_ratio;
  dxdt_S100  = kin_S100*(IL17/IL17_base)*(TNFa/TNF_base)*(1.0-0.45*Itot_IL17-0.35*Itot_TNF-0.20*Itot_IL23)-kout_S100*S100;
  dxdt_PASI  = kin_PASI*(IL17/IL17_base)*(TNFa/TNF_base)*(1.0-0.65*Itot_IL17-0.20*Itot_IL23-0.10*Itot_TNF)-kout_PASI*PASI;
  dxdt_DAPSA = kin_DAPSA*(TNFa/TNF_base)*(IL17/IL17_base)*(1.0-0.45*Itot_TNF-0.35*Itot_IL17-0.15*Itot_IL23)-kout_DAPSA*DAPSA;
  dxdt_MTSS  = kin_MTSS*(RANKL/RANKL_base)*(1.0-0.7*(0.5*Itot_IL17+0.4*Itot_TNF));

$TABLE
  double CONC_ADA = ADA_C1/V1_ADA;
  double CONC_IXE = IXE_C1/V1_IXE;
  double CONC_GUS = GUS_C1/V1_GUS;
  double CONC_UPA = (UPA_C1/V1_UPA)*1000.0;
  double CONC_APR = (APR_C1/V1_APR)*1000.0;
  double CRP_mgL  = CRP_ratio*20.0;
  double Calprotectin = S100*3.5;
  double PASI_chg = (18.0>0) ? (18.0-PASI)/18.0*100.0 : 0.0;
  double DAPSA_chg= (28.0>0) ? (28.0-DAPSA)/28.0*100.0 : 0.0;
  int PASI75 = (PASI_chg>=75.0) ? 1:0;
  int PASI90 = (PASI_chg>=90.0) ? 1:0;
  int ACR20  = (DAPSA_chg>=20.0) ? 1:0;
  int ACR50  = (DAPSA_chg>=50.0) ? 1:0;
  int ACR70  = (DAPSA_chg>=70.0) ? 1:0;
  int DAPSA_REM = (DAPSA<=4.0) ? 1:0;
  capture CONC_ADA CONC_IXE CONC_GUS CONC_UPA CONC_APR
  capture IL17 TNFa IL23 TH17 RANKL CRP_mgL Calprotectin
  capture PASI DAPSA MTSS PASI_chg DAPSA_chg
  capture PASI75 PASI90 ACR20 ACR50 ACR70 DAPSA_REM
'

## ---- Dosing helper --------------------------------------
make_ev <- function(drug, weeks) {
  h <- weeks * 7 * 24
  switch(drug,
    "Adalimumab"   = ev(cmt=1, amt=40, ii=14*24, addl=floor(h/(14*24))),
    "Ixekizumab"   = {
      e1 <- ev(cmt=5, amt=c(160, rep(80,8)), time=seq(0,16*7*24,by=2*7*24))
      e2 <- ev(cmt=5, amt=80, time=seq(20*7*24,h,by=4*7*24))
      c(e1, e2)
    },
    "Guselkumab"   = ev(cmt=8, amt=100, time=c(0,4*7*24,seq(12*7*24,h,by=8*7*24))),
    "Upadacitinib" = ev(cmt=11, amt=15, ii=24, addl=floor(h/24)),
    "Apremilast"   = ev(cmt=13, amt=30, ii=12, addl=floor(h/12)),
    ev(cmt=1, amt=0, time=0)
  )
}

## ---- Compile model once ---------------------------------
mod <- mread("psa_shiny", tempdir(), psa_code)

## ========================================================
## UI
## ========================================================
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "PsA QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",  tabName = "profile",   icon = icon("user")),
      menuItem("Drug PK",          tabName = "pk",        icon = icon("pills")),
      menuItem("Cytokine Dynamics",tabName = "cytokines", icon = icon("dna")),
      menuItem("Skin (PASI)",      tabName = "skin",      icon = icon("band-aid")),
      menuItem("Joint (DAPSA)",    tabName = "joint",     icon = icon("bone")),
      menuItem("Scenario Compare", tabName = "scenario",  icon = icon("chart-bar")),
      menuItem("Biomarkers",       tabName = "biomarkers",icon = icon("vial"))
    )
  ),

  dashboardBody(
    tabItems(

      ## ---- Tab 1: Patient Profile -------------------------
      tabItem(tabName = "profile",
        fluidRow(
          box(title = "Patient Parameters", status = "primary", solidHeader = TRUE, width = 4,
            sliderInput("bwt",  "Body Weight (kg)", 40, 150, 85, step = 5),
            sliderInput("pasi_init", "Baseline PASI",  5, 50, 18, step = 1),
            sliderInput("dapsa_init","Baseline DAPSA", 5, 80, 28, step = 1),
            sliderInput("crp_init",  "Baseline CRP (mg/L)", 5, 100, 20, step = 1),
            selectInput("disease_duration", "Disease Duration",
                        choices = c("<1 yr","1-3 yrs","3-10 yrs",">10 yrs"),
                        selected = "1-3 yrs"),
            checkboxGroupInput("comorbid", "Comorbidities",
                               choices = c("Obesity","MetSyn","IBD","Uveitis","CV Risk"),
                               selected = c())
          ),
          box(title = "Disease Classification (CASPAR)", status = "info", solidHeader = TRUE, width = 4,
            checkboxInput("caspar_ps",     "Current psoriasis (2 pts)", TRUE),
            checkboxInput("caspar_nail",   "Nail dystrophy (1 pt)", FALSE),
            checkboxInput("caspar_rf_neg", "RF negative (1 pt)", TRUE),
            checkboxInput("caspar_dactyl", "Dactylitis (1 pt)", FALSE),
            checkboxInput("caspar_xray",   "Juxta-articular bone formation (1 pt)", FALSE),
            hr(),
            verbatimTextOutput("caspar_score")
          ),
          box(title = "Simulation Settings", status = "success", solidHeader = TRUE, width = 4,
            selectInput("treatment", "Treatment",
                        choices = c("Placebo","Adalimumab","Ixekizumab",
                                    "Guselkumab","Upadacitinib","Apremilast"),
                        selected = "Ixekizumab"),
            sliderInput("sim_weeks", "Simulation Duration (weeks)", 12, 104, 52, step = 4),
            actionButton("run_sim", "Run Simulation", class = "btn-success btn-lg", width = "100%")
          )
        ),
        fluidRow(
          valueBoxOutput("vbox_pasi"),
          valueBoxOutput("vbox_dapsa"),
          valueBoxOutput("vbox_crp")
        )
      ),

      ## ---- Tab 2: Drug PK ---------------------------------
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "PK Parameters", status = "primary", solidHeader = TRUE, width = 3,
            numericInput("dose_val", "Dose (mg)", value = 80, min = 1),
            selectInput("pk_drug", "Drug to Display",
                        choices = c("Adalimumab","Ixekizumab","Guselkumab",
                                    "Upadacitinib","Apremilast")),
            checkboxInput("show_trough", "Show trough lines", TRUE),
            checkboxInput("log_pk", "Log Y-axis", FALSE)
          ),
          box(title = "Concentration–Time Profile", status = "primary", solidHeader = TRUE, width = 9,
            plotOutput("pk_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "PK Parameters Table", status = "info", solidHeader = TRUE, width = 12,
            DTOutput("pk_params_table")
          )
        )
      ),

      ## ---- Tab 3: Cytokine Dynamics -----------------------
      tabItem(tabName = "cytokines",
        fluidRow(
          box(title = "IL-17A Dynamics", status = "primary", solidHeader = TRUE, width = 6,
            plotOutput("il17_plot", height = "300px")
          ),
          box(title = "TNF-α Dynamics", status = "warning", solidHeader = TRUE, width = 6,
            plotOutput("tnf_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "IL-23 Dynamics", status = "danger", solidHeader = TRUE, width = 6,
            plotOutput("il23_plot", height = "300px")
          ),
          box(title = "Th17 Cell Dynamics", status = "info", solidHeader = TRUE, width = 6,
            plotOutput("th17_plot", height = "300px")
          )
        )
      ),

      ## ---- Tab 4: Skin (PASI) -----------------------------
      tabItem(tabName = "skin",
        fluidRow(
          box(title = "PASI Over Time", status = "warning", solidHeader = TRUE, width = 8,
            plotOutput("pasi_plot", height = "380px")
          ),
          box(title = "Response Summary", status = "success", solidHeader = TRUE, width = 4,
            h4("At Week 12"),
            verbatimTextOutput("pasi_w12"),
            hr(),
            h4("At Week 24"),
            verbatimTextOutput("pasi_w24"),
            hr(),
            h4("At Week 52"),
            verbatimTextOutput("pasi_w52")
          )
        ),
        fluidRow(
          box(title = "PASI75/90 Achievement Over Time", status = "info", solidHeader = TRUE, width = 12,
            plotOutput("pasi_resp_plot", height = "280px")
          )
        )
      ),

      ## ---- Tab 5: Joint (DAPSA) ---------------------------
      tabItem(tabName = "joint",
        fluidRow(
          box(title = "DAPSA Over Time", status = "danger", solidHeader = TRUE, width = 8,
            plotOutput("dapsa_plot", height = "380px")
          ),
          box(title = "Joint Biomarkers", status = "info", solidHeader = TRUE, width = 4,
            plotOutput("rankl_plot", height = "180px"),
            plotOutput("mtss_plot",  height = "180px")
          )
        ),
        fluidRow(
          box(title = "ACR20/50/70 Over Time", status = "primary", solidHeader = TRUE, width = 12,
            plotOutput("acr_plot", height = "280px")
          )
        )
      ),

      ## ---- Tab 6: Scenario Comparison ---------------------
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "Select Scenarios to Compare", status = "primary", solidHeader = TRUE, width = 3,
            checkboxGroupInput("compare_scenarios", "Treatment Arms",
                               choices = c("Placebo","Adalimumab","Ixekizumab",
                                           "Guselkumab","Upadacitinib","Apremilast"),
                               selected = c("Placebo","Adalimumab","Ixekizumab","Guselkumab")),
            selectInput("compare_endpoint", "Primary Endpoint",
                        choices = c("PASI","DAPSA","IL17","CRP_mgL","MTSS"),
                        selected = "PASI"),
            actionButton("run_compare", "Compare All", class = "btn-primary", width = "100%")
          ),
          box(title = "Comparative Efficacy", status = "primary", solidHeader = TRUE, width = 9,
            plotOutput("compare_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Efficacy Summary Table (Week 24)", status = "info", solidHeader = TRUE, width = 12,
            DTOutput("summary_table")
          )
        )
      ),

      ## ---- Tab 7: Biomarkers ------------------------------
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "CRP (mg/L)", status = "warning", solidHeader = TRUE, width = 4,
            plotOutput("crp_plot", height = "260px")
          ),
          box(title = "Calprotectin (µg/mL)", status = "info", solidHeader = TRUE, width = 4,
            plotOutput("calp_plot", height = "260px")
          ),
          box(title = "RANKL Ratio", status = "danger", solidHeader = TRUE, width = 4,
            plotOutput("rankl_bm_plot", height = "260px")
          )
        ),
        fluidRow(
          box(title = "Biomarker Table (Week 4/12/24/52)", status = "success", solidHeader = TRUE, width = 12,
            DTOutput("bm_table")
          )
        )
      )
    )
  )
)

## ========================================================
## SERVER
## ========================================================
server <- function(input, output, session) {

  ## ---- CASPAR score ------------------------------------
  output$caspar_score <- renderText({
    pts <- 0
    if (input$caspar_ps)     pts <- pts + 2
    if (input$caspar_nail)   pts <- pts + 1
    if (input$caspar_rf_neg) pts <- pts + 1
    if (input$caspar_dactyl) pts <- pts + 1
    if (input$caspar_xray)   pts <- pts + 1
    cat_txt <- if (pts >= 3) "✓ PsA (CASPAR ≥3)" else "✗ Not PsA (CASPAR <3)"
    paste("CASPAR Score:", pts, "\n", cat_txt)
  })

  ## ---- Reactive simulation (single drug) ---------------
  sim_result <- eventReactive(input$run_sim, {
    req(input$treatment)
    ev_dose <- make_ev(input$treatment, input$sim_weeks)
    upd_mod <- param(mod,
                     PASI_base = input$pasi_init,
                     DAPSA_base= input$dapsa_init)
    upd_mod <- init(upd_mod,
                    PASI  = input$pasi_init,
                    DAPSA = input$dapsa_init,
                    CRP_ratio = input$crp_init / 20.0)

    out <- upd_mod %>%
      ev(ev_dose) %>%
      mrgsim(end = input$sim_weeks * 7 * 24, delta = 24, obsonly = TRUE) %>%
      as.data.frame()
    out$Week <- out$time / (7 * 24)
    out
  }, ignoreNULL = FALSE)

  ## Initialize with default run
  observe({
    if (!is.null(sim_result())) return()
    ev_dose <- make_ev("Ixekizumab", 52)
    out <- mod %>% ev(ev_dose) %>%
      mrgsim(end = 52*7*24, delta=24, obsonly=TRUE) %>% as.data.frame()
    out$Week <- out$time / (7*24)
  })

  ## ---- Value boxes ----------------------------------------
  output$vbox_pasi <- renderValueBox({
    df <- sim_result()
    if (is.null(df)) return(valueBox(18, "Baseline PASI", icon = icon("band-aid"), color = "orange"))
    last <- tail(df, 1)
    valueBox(round(last$PASI, 1),
             paste("PASI at wk", round(last$Week)),
             icon = icon("band-aid"),
             color = if (last$PASI < 3) "green" else if (last$PASI < 9) "yellow" else "red")
  })
  output$vbox_dapsa <- renderValueBox({
    df <- sim_result()
    if (is.null(df)) return(valueBox(28, "Baseline DAPSA", icon = icon("bone"), color = "red"))
    last <- tail(df, 1)
    valueBox(round(last$DAPSA, 1),
             paste("DAPSA at wk", round(last$Week)),
             icon = icon("bone"),
             color = if (last$DAPSA <= 4) "green" else if (last$DAPSA <= 14) "yellow" else "red")
  })
  output$vbox_crp <- renderValueBox({
    df <- sim_result()
    if (is.null(df)) return(valueBox(20, "Baseline CRP", icon = icon("vial"), color = "orange"))
    last <- tail(df, 1)
    valueBox(round(last$CRP_mgL, 1),
             paste("CRP (mg/L) at wk", round(last$Week)),
             icon = icon("vial"),
             color = if (last$CRP_mgL < 5) "green" else if (last$CRP_mgL < 10) "yellow" else "red")
  })

  ## ---- PK plot ------------------------------------------
  output$pk_plot <- renderPlot({
    df <- sim_result()
    if (is.null(df)) return()
    conc_col <- switch(input$pk_drug,
      "Adalimumab"   = "CONC_ADA",
      "Ixekizumab"   = "CONC_IXE",
      "Guselkumab"   = "CONC_GUS",
      "Upadacitinib" = "CONC_UPA",
      "Apremilast"   = "CONC_APR"
    )
    if (!conc_col %in% names(df)) {
      plot(0, type="n", main="No PK data for this drug in current scenario")
      return()
    }
    p <- ggplot(df, aes_string(x = "Week", y = conc_col)) +
      geom_line(color = "steelblue", size = 1.2) +
      labs(title = paste(input$pk_drug, "PK – Concentration–Time"),
           x = "Week", y = "Concentration (µg/mL)") +
      theme_bw(base_size = 13)
    if (input$log_pk) p <- p + scale_y_log10()
    print(p)
  })

  ## ---- PK params table -----------------------------------
  output$pk_params_table <- renderDT({
    pk_tbl <- data.frame(
      Drug = c("Adalimumab","Ixekizumab","Guselkumab","Upadacitinib","Apremilast"),
      Route = c("SC","SC","SC","Oral","Oral"),
      Dose  = c("40mg Q2W","80mg Q4W","100mg Q8W","15mg QD","30mg BID"),
      t_half_days = c(14, 13, 17, 0.5, 0.4),
      Vd_L  = c(2.75+1.76, 3.8, 3.2, 220, 84),
      F_pct = c(64, 60, 49, 76, 73),
      Target= c("TNF-α","IL-17A/F","IL-23p19","JAK1/2/3/TYK2","PDE4"),
      Trial = c("ADEPT","SPIRIT-P1","DISCOVER-1","SELECT-PsA 1","PALACE 1-3")
    )
    datatable(pk_tbl, options = list(pageLength = 5, scrollX = TRUE), rownames = FALSE)
  })

  ## ---- Cytokine plots ------------------------------------
  mk_cyto_plot <- function(df, col, title, color) {
    if (is.null(df) || !col %in% names(df)) return(ggplot() + theme_void())
    ggplot(df, aes_string(x="Week", y=col)) +
      geom_line(color = color, size = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dashed", alpha = 0.5) +
      labs(title = title, x = "Week", y = "Relative Concentration") +
      theme_bw(base_size = 11)
  }
  output$il17_plot <- renderPlot({ mk_cyto_plot(sim_result(), "IL17",  "IL-17A", "#E91E63") })
  output$tnf_plot  <- renderPlot({ mk_cyto_plot(sim_result(), "TNFa",  "TNF-α",  "#FF5722") })
  output$il23_plot <- renderPlot({ mk_cyto_plot(sim_result(), "IL23",  "IL-23",  "#9C27B0") })
  output$th17_plot <- renderPlot({ mk_cyto_plot(sim_result(), "TH17",  "Th17 Cells", "#3F51B5") })

  ## ---- PASI plots ----------------------------------------
  output$pasi_plot <- renderPlot({
    df <- sim_result()
    if (is.null(df)) return()
    base_pasi <- input$pasi_init
    ggplot(df, aes(x = Week, y = PASI)) +
      geom_line(color = "#FF9800", size = 1.3) +
      geom_hline(yintercept = base_pasi * 0.25, linetype="dashed", color="gray50") +
      geom_hline(yintercept = base_pasi * 0.10, linetype="dashed", color="gray70") +
      annotate("text", x = max(df$Week)*0.9, y = base_pasi*0.25+0.5,
               label = "PASI75 threshold", size = 3.5) +
      annotate("text", x = max(df$Week)*0.9, y = base_pasi*0.10+0.5,
               label = "PASI90 threshold", size = 3.5) +
      labs(title = paste("PASI –", input$treatment),
           x = "Week", y = "PASI Score (0–72)") +
      theme_bw(base_size = 13)
  })

  pasi_summary <- function(df, wk) {
    if (is.null(df)) return("No data")
    row <- df[which.min(abs(df$Week - wk)), ]
    paste0("PASI: ", round(row$PASI, 1), "\n",
           "PASI75: ", row$PASI75, " | PASI90: ", row$PASI90, "\n",
           "% Change: ", round(row$PASI_chg, 1), "%")
  }
  output$pasi_w12 <- renderText({ pasi_summary(sim_result(), 12) })
  output$pasi_w24 <- renderText({ pasi_summary(sim_result(), 24) })
  output$pasi_w52 <- renderText({ pasi_summary(sim_result(), 52) })

  output$pasi_resp_plot <- renderPlot({
    df <- sim_result()
    if (is.null(df)) return()
    df_long <- df %>%
      select(Week, PASI75, PASI90) %>%
      pivot_longer(c(PASI75, PASI90), names_to = "Criterion", values_to = "Achieved")
    ggplot(df_long, aes(x = Week, y = Achieved, color = Criterion)) +
      geom_line(size = 1.2) +
      scale_y_continuous(limits = c(-0.05, 1.05), breaks = c(0,1), labels = c("No","Yes")) +
      scale_color_manual(values = c(PASI75="#FF9800", PASI90="#E91E63")) +
      labs(title = "PASI75/90 Achievement Over Time",
           x = "Week", y = "Response Achieved") +
      theme_bw(base_size = 12)
  })

  ## ---- DAPSA / joint plots -------------------------------
  output$dapsa_plot <- renderPlot({
    df <- sim_result()
    if (is.null(df)) return()
    ggplot(df, aes(x = Week, y = DAPSA)) +
      geom_line(color = "#F44336", size = 1.3) +
      geom_hline(yintercept = 4,  linetype="dashed", color="green3") +
      geom_hline(yintercept = 14, linetype="dashed", color="orange") +
      geom_hline(yintercept = 28, linetype="dashed", color="gray60") +
      annotate("text", x=max(df$Week)*0.85, y=4.8,  label="Remission (≤4)",  size=3.5, color="green4") +
      annotate("text", x=max(df$Week)*0.85, y=14.8, label="LDA (≤14)",        size=3.5, color="darkorange") +
      labs(title = paste("DAPSA –", input$treatment),
           x = "Week", y = "DAPSA Score") +
      theme_bw(base_size = 13)
  })

  output$rankl_plot <- renderPlot({
    df <- sim_result()
    if (is.null(df)) return()
    ggplot(df, aes(x=Week, y=RANKL)) +
      geom_line(color="#795548", size=1) +
      labs(title="RANKL Ratio", x="Week", y="RANKL (relative)") +
      theme_bw(base_size=10)
  })

  output$mtss_plot <- renderPlot({
    df <- sim_result()
    if (is.null(df)) return()
    ggplot(df, aes(x=Week, y=MTSS)) +
      geom_line(color="#607D8B", size=1) +
      labs(title="Structural Damage (mTSS proxy)", x="Week", y="mTSS") +
      theme_bw(base_size=10)
  })

  output$acr_plot <- renderPlot({
    df <- sim_result()
    if (is.null(df)) return()
    df_long <- df %>%
      select(Week, ACR20, ACR50, ACR70) %>%
      pivot_longer(c(ACR20,ACR50,ACR70), names_to="Criterion", values_to="Achieved")
    ggplot(df_long, aes(x=Week, y=Achieved, color=Criterion)) +
      geom_line(size=1.2) +
      scale_y_continuous(limits=c(-0.05,1.05), breaks=c(0,1), labels=c("No","Yes")) +
      scale_color_manual(values=c(ACR20="#4CAF50",ACR50="#2196F3",ACR70="#FF9800")) +
      labs(title="ACR Response Achievement Over Time", x="Week", y="Response") +
      theme_bw(base_size=12)
  })

  ## ---- Scenario comparison --------------------------------
  compare_data <- eventReactive(input$run_compare, {
    scs <- input$compare_scenarios
    if (length(scs) == 0) return(NULL)
    results <- lapply(scs, function(sc) {
      ev_dose <- make_ev(sc, input$sim_weeks)
      out <- mod %>% ev(ev_dose) %>%
        mrgsim(end = input$sim_weeks*7*24, delta=24, obsonly=TRUE) %>%
        as.data.frame()
      out$Scenario <- sc
      out$Week <- out$time/(7*24)
      out
    })
    bind_rows(results)
  }, ignoreNULL = FALSE)

  output$compare_plot <- renderPlot({
    df <- compare_data()
    if (is.null(df)) return()
    ep <- input$compare_endpoint
    if (!ep %in% names(df)) return()
    ggplot(df, aes_string(x="Week", y=ep, color="Scenario")) +
      geom_line(size=1.2) +
      scale_color_brewer(palette="Set1") +
      labs(title=paste("Comparative:", ep, "over time"),
           x="Week", y=ep, color="Treatment") +
      theme_bw(base_size=13)
  })

  output$summary_table <- renderDT({
    df <- compare_data()
    if (is.null(df)) return(datatable(data.frame()))
    df %>%
      filter(abs(Week - 24) < 0.5) %>%
      group_by(Scenario) %>%
      summarise(
        PASI_24     = round(mean(PASI), 1),
        PASI75_pct  = round(mean(PASI75)*100, 0),
        PASI90_pct  = round(mean(PASI90)*100, 0),
        DAPSA_24    = round(mean(DAPSA), 1),
        ACR20_pct   = round(mean(ACR20)*100, 0),
        ACR50_pct   = round(mean(ACR50)*100, 0),
        ACR70_pct   = round(mean(ACR70)*100, 0),
        DAPSA_REM_pct = round(mean(DAPSA_REM)*100, 0),
        CRP_mgL     = round(mean(CRP_mgL), 1),
        MTSS_24     = round(mean(MTSS), 4),
        .groups="drop"
      ) %>%
      datatable(options=list(pageLength=8, scrollX=TRUE), rownames=FALSE)
  })

  ## ---- Biomarker plots ------------------------------------
  output$crp_plot <- renderPlot({
    df <- sim_result()
    if (is.null(df)) return()
    ggplot(df, aes(x=Week, y=CRP_mgL)) +
      geom_line(color="#FF5722", size=1.2) +
      geom_hline(yintercept=5, linetype="dashed", color="gray50") +
      annotate("text", x=max(df$Week)*0.8, y=6, label="Upper Normal (5 mg/L)", size=3.5) +
      labs(title="CRP (mg/L)", x="Week", y="CRP (mg/L)") +
      theme_bw(base_size=12)
  })

  output$calp_plot <- renderPlot({
    df <- sim_result()
    if (is.null(df)) return()
    ggplot(df, aes(x=Week, y=Calprotectin)) +
      geom_line(color="#2196F3", size=1.2) +
      geom_hline(yintercept=3.5, linetype="dashed", color="gray50") +
      labs(title="Calprotectin (µg/mL)", x="Week", y="S100A8/9 (µg/mL)") +
      theme_bw(base_size=12)
  })

  output$rankl_bm_plot <- renderPlot({
    df <- sim_result()
    if (is.null(df)) return()
    ggplot(df, aes(x=Week, y=RANKL)) +
      geom_line(color="#795548", size=1.2) +
      geom_hline(yintercept=1.0, linetype="dashed", color="gray50") +
      labs(title="RANKL Ratio", x="Week", y="RANKL/Baseline") +
      theme_bw(base_size=12)
  })

  output$bm_table <- renderDT({
    df <- sim_result()
    if (is.null(df)) return(datatable(data.frame()))
    df %>%
      filter(Week %in% c(0,4,12,24,52)) %>%
      mutate(Week = round(Week)) %>%
      select(Week, CRP_mgL, Calprotectin, IL17, TNFa, IL23, RANKL, MTSS) %>%
      mutate(across(where(is.numeric), ~round(.x, 3))) %>%
      datatable(options=list(pageLength=6), rownames=FALSE)
  })
}

## ========================================================
## Run app
## ========================================================
shinyApp(ui, server)
