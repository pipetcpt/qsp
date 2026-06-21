## =============================================================================
## Psoriasis QSP Interactive Shiny Dashboard
## =============================================================================
## 6 Tabs: Patient Profile · PK · IL-17/Cytokine PD · PASI Endpoints ·
##         Scenario Comparison · Biomarker Dashboard
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(DT)

## ─────────────────────────────────────────────────────────────────────────────
## EMBED MODEL CODE (same ODE as mrgsolve R file — simplified for app)
## ─────────────────────────────────────────────────────────────────────────────
pso_model_code <- '
$PARAM
ka_ADA=0.013 CL_ADA=0.247 V1_ADA=7.0 Q_ADA=0.30 V2_ADA=3.2 Kd_ADA=0.1
ka_SEC=0.015 CL_SEC=0.191 V1_SEC=7.1 Q_SEC=0.20 V2_SEC=3.8 Kd_SEC=0.08
ka_RSK=0.012 CL_RSK=0.078 V1_RSK=11.2 Kd_RSK=0.06
ka_APR=0.80  CL_APR=9.5   V1_APR=86.6 Q_APR=3.9  V2_APR=43.3 IC50_APR=74.0 Emax_APR=0.70
ka_TOF=1.20  CL_TOF=22.8  V1_TOF=87.0 IC50_TOF=1.0 Emax_TOF=0.80
ka_MTX=0.50  CL_MTX=4.8   V1_MTX=24.0 k_PG=0.02 k_PGelim=0.005 IC50_MTX=5.0
k_DC=0.05 k_DCd=0.02 k_IL23p=0.10 k_IL23d=0.50
k_Th17d=0.003 k_Th17x=0.005 EC50_IL23=50 Emax_IL23=5
k_IL17p=0.20 k_IL17d=0.30 EC50_IL17=80
k_TNFp=0.15 k_TNFd=0.40 k_IFNp=0.12 k_IFNd=0.35
k_KCb=0.008 k_KCx=0.004 k_KC17=0.006 k_KCtnf=0.004
k_PASIf=0.0010 k_PASIr=0.0050 PASI_ss=20
MW_ADA=148000 MW_SEC=147000 MW_RSK=153000

$CMT
DC IL23 Th17 IL17A TNFa IFNg KC PASI
ADA_SC ADA_C ADA_P
SEC_SC SEC_C SEC_P
RSK_SC RSK_C
APR_GI APR_C APR_P
TOF_GI TOF_C
MTX_GI MTX_C MTX_PG

$MAIN
double ADA_nM  = (ADA_C/V1_ADA)*1e6/MW_ADA;
double SEC_nM  = (SEC_C/V1_SEC)*1e6/MW_SEC;
double RSK_nM  = (RSK_C/V1_RSK)*1e6/MW_RSK;
double APR_nM  = (APR_C/V1_APR)*1e6/460.5;
double TOF_nM  = (TOF_C/V1_TOF)*1e6/312.4;
double MTX_nM  = MTX_PG/0.454;
double fADA    = ADA_nM/(ADA_nM+Kd_ADA);
double fSEC    = SEC_nM/(SEC_nM+Kd_SEC);
double fRSK    = RSK_nM/(RSK_nM+Kd_RSK);
double fAPR    = Emax_APR*APR_nM/(APR_nM+IC50_APR);
double fTOF    = Emax_TOF*TOF_nM/(TOF_nM+IC50_TOF);
double fMTX    = MTX_nM/(MTX_nM+IC50_MTX);
double IL23e   = IL23*(1-fRSK)*(1-fTOF*0.4);
double IL17e   = IL17A*(1-fSEC);
double TNFe    = TNFa*(1-fADA);
double Th17_dr = k_Th17d*(1+Emax_IL23*IL23e/(IL23e+EC50_IL23));
double KC_st   = k_KC17*IL17e/(IL17e+EC50_IL17)+k_KCtnf*TNFe/(TNFe+100.0);
double Th0ss   = 500.0;

$ODE
dxdt_DC    = k_DC - k_DCd*DC;
dxdt_IL23  = k_IL23p*DC*(1-fRSK) - k_IL23d*IL23;
dxdt_Th17  = Th17_dr*Th0ss*(1-fMTX*0.5)*(1-fTOF*0.6) - k_Th17x*Th17;
dxdt_IL17A = k_IL17p*Th17*(1-fAPR)*(1-fTOF*0.3) - k_IL17d*IL17A;
dxdt_TNFa  = (k_TNFp*DC+0.05*Th17)*(1-fAPR)*(1-fTOF*0.5) - k_TNFd*TNFa;
dxdt_IFNg  = k_IFNp*Th17 - k_IFNd*IFNg;
dxdt_KC    = KC_st*(100.0-KC) - k_KCx*(KC-100.0);
dxdt_PASI  = k_PASIf*KC - k_PASIr*PASI;
dxdt_ADA_SC = -ka_ADA*ADA_SC;
dxdt_ADA_C  =  ka_ADA*ADA_SC-(CL_ADA+Q_ADA)*(ADA_C/V1_ADA)+Q_ADA*(ADA_P/V2_ADA);
dxdt_ADA_P  =  Q_ADA*(ADA_C/V1_ADA)-Q_ADA*(ADA_P/V2_ADA);
dxdt_SEC_SC = -ka_SEC*SEC_SC;
dxdt_SEC_C  =  ka_SEC*SEC_SC-(CL_SEC+Q_SEC)*(SEC_C/V1_SEC)+Q_SEC*(SEC_P/V2_SEC);
dxdt_SEC_P  =  Q_SEC*(SEC_C/V1_SEC)-Q_SEC*(SEC_P/V2_SEC);
dxdt_RSK_SC = -ka_RSK*RSK_SC;
dxdt_RSK_C  =  ka_RSK*RSK_SC-CL_RSK*(RSK_C/V1_RSK);
dxdt_APR_GI = -ka_APR*APR_GI;
dxdt_APR_C  =  ka_APR*APR_GI-(CL_APR/V1_APR+Q_APR/V1_APR)*APR_C+Q_APR/V2_APR*APR_P;
dxdt_APR_P  =  Q_APR/V1_APR*APR_C-Q_APR/V2_APR*APR_P;
dxdt_TOF_GI = -ka_TOF*TOF_GI;
dxdt_TOF_C  =  ka_TOF*TOF_GI-CL_TOF/V1_TOF*TOF_C;
dxdt_MTX_GI = -ka_MTX*MTX_GI;
dxdt_MTX_C  =  ka_MTX*MTX_GI-CL_MTX/V1_MTX*MTX_C;
dxdt_MTX_PG =  k_PG*MTX_C-k_PGelim*MTX_PG;

$TABLE
double Th0ss = 500.0;
double ADA_nM = (ADA_C/V1_ADA)*1e6/MW_ADA;
double SEC_nM = (SEC_C/V1_SEC)*1e6/MW_SEC;
double RSK_nM = (RSK_C/V1_RSK)*1e6/MW_RSK;
double APR_nM = (APR_C/V1_APR)*1e6/460.5;
double TOF_nM = (TOF_C/V1_TOF)*1e6/312.4;
double PASI75  = (PASI <= 0.25*20.0) ? 1.0 : 0.0;
double PASI90  = (PASI <= 0.10*20.0) ? 1.0 : 0.0;
double PASI100 = (PASI <= 0.10) ? 1.0 : 0.0;
double IGA01   = (PASI <= 3.0) ? 1.0 : 0.0;

$CAPTURE
PASI PASI75 PASI90 PASI100 IGA01 IL17A TNFa Th17 IL23 KC IFNg
ADA_nM SEC_nM RSK_nM APR_nM TOF_nM
'

pso_mod <- mrgsolve::mcode("pso_shiny", pso_model_code)

## ─────────────────────────────────────────────────────────────────────────────
## HELPER: build dosing events
## ─────────────────────────────────────────────────────────────────────────────
build_events <- function(drug, dose_mg, weeks_total) {
  sim_end <- weeks_total * 7 * 24
  if (drug == "adalimumab") {
    mrgsolve::ev(
      amt  = c(80000, 40000, rep(40000, ceiling(weeks_total/2))),
      time = c(0, 168, seq(336, sim_end, by=336)),
      cmt  = "ADA_SC"
    )
  } else if (drug == "secukinumab") {
    mrgsolve::ev(
      amt  = c(rep(300000,5), rep(300000, ceiling(weeks_total/4))),
      time = c(0,168,336,504,672, seq(672+4*168, sim_end, by=4*168)),
      cmt  = "SEC_SC"
    )
  } else if (drug == "risankizumab") {
    mrgsolve::ev(
      amt  = c(150000, 150000, rep(150000, ceiling(weeks_total/12))),
      time = c(0, 4*168, seq(16*168, sim_end, by=12*168)),
      cmt  = "RSK_C"
    )
  } else if (drug == "apremilast") {
    mrgsolve::ev(
      amt  = rep(dose_mg * 1000, 2 * weeks_total * 7),
      time = seq(0, sim_end - 12, by=12),
      cmt  = "APR_GI"
    )
  } else if (drug == "tofacitinib") {
    mrgsolve::ev(
      amt  = rep(dose_mg * 1000, 2 * weeks_total * 7),
      time = seq(0, sim_end - 12, by=12),
      cmt  = "TOF_GI"
    )
  } else if (drug == "methotrexate") {
    mrgsolve::ev(
      amt  = rep(dose_mg * 1000, weeks_total),
      time = seq(0, sim_end - 168, by=168),
      cmt  = "MTX_GI"
    )
  } else {
    mrgsolve::ev()  # no treatment
  }
}

## ─────────────────────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "Psoriasis QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName="tab_patient",   icon=icon("user")),
      menuItem("PK Profile",         tabName="tab_pk",        icon=icon("pills")),
      menuItem("Cytokine PD",        tabName="tab_pd",        icon=icon("bacterium")),
      menuItem("PASI Endpoints",     tabName="tab_pasi",      icon=icon("chart-line")),
      menuItem("Scenario Comparison",tabName="tab_scenario",  icon=icon("exchange-alt")),
      menuItem("Biomarker Dashboard",tabName="tab_bm",        icon=icon("vial"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box-header { font-weight:bold; }
      .main-sidebar { background:#4A148C; }
    "))),

    tabItems(
      ## ─── TAB 1: PATIENT PROFILE ───────────────────────────────────────────
      tabItem("tab_patient",
        fluidRow(
          box(title="Patient & Disease Parameters", width=4, status="purple",
            sliderInput("base_pasi",  "Baseline PASI", 5, 72, 20, 1),
            sliderInput("base_il17",  "Baseline IL-17A (pg/mL)", 10, 500, 140, 10),
            sliderInput("base_tnf",   "Baseline TNF-α (pg/mL)",  10, 300, 80, 5),
            sliderInput("base_th17",  "Baseline Th17 (×1000/mL)", 50, 1000, 350, 25),
            selectInput("psoriasis_type","Psoriasis Phenotype",
                        choices=c("Plaque (moderate-severe)"="plaque",
                                  "Generalized Pustular (GPP)"="gpp",
                                  "Erythrodermic"="erythro",
                                  "Guttate (post-strep)"="guttate"))
          ),
          box(title="Treatment Selection", width=4, status="primary",
            selectInput("drug","Drug / Mechanism",
                        choices=c("No Treatment"="none",
                                  "Adalimumab (anti-TNF)"="adalimumab",
                                  "Secukinumab (anti-IL17A)"="secukinumab",
                                  "Risankizumab (anti-IL23p19)"="risankizumab",
                                  "Apremilast (PDE4i)"="apremilast",
                                  "Tofacitinib (JAKi)"="tofacitinib",
                                  "Methotrexate (DHFR)"="methotrexate")),
            numericInput("dose_mg","Dose (mg) [small molecules only]", 30, min=5, max=100),
            sliderInput("sim_weeks","Simulation Duration (weeks)", 12, 104, 52, 4),
            actionButton("run_sim","Run Simulation",
                         class="btn-primary btn-lg", icon=icon("play"))
          ),
          box(title="Disease Summary", width=4, status="warning",
            valueBoxOutput("vb_pasi",   width=12),
            valueBoxOutput("vb_sev",    width=12),
            valueBoxOutput("vb_drug",   width=12)
          )
        ),
        fluidRow(
          box(title="Psoriasis Pathophysiology — IL-23/IL-17 Axis", width=12,
              status="info",
              p("Psoriasis is driven by the IL-23/Th17/IL-17A axis. Myeloid DCs",
                "produce IL-23, which drives Th17 differentiation. Th17 cells secrete",
                "IL-17A (and IL-22), which activates keratinocytes via NF-κB to",
                "hyperproliferate, producing acanthosis and parakeratosis.",
                "TNF-α amplifies vascular changes and keratinocyte activation."),
              tableOutput("tbl_patient_summary")
          )
        )
      ),

      ## ─── TAB 2: PK PROFILE ───────────────────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(title="Biologic PK — Central Concentration", width=6, status="primary",
              plotlyOutput("plt_bio_pk", height="350px")),
          box(title="Small Molecule PK — Trough Concentration", width=6, status="success",
              plotlyOutput("plt_sm_pk", height="350px"))
        ),
        fluidRow(
          box(title="PK Parameters Reference", width=12, status="info",
              DTOutput("tbl_pk_params"))
        )
      ),

      ## ─── TAB 3: CYTOKINE PD ─────────────────────────────────────────────
      tabItem("tab_pd",
        fluidRow(
          box(title="IL-17A Dynamics (Free Serum)", width=6, status="warning",
              plotlyOutput("plt_il17", height="350px")),
          box(title="TNF-α Dynamics", width=6, status="danger",
              plotlyOutput("plt_tnf", height="350px"))
        ),
        fluidRow(
          box(title="Th17 Cell Count", width=6, status="info",
              plotlyOutput("plt_th17", height="350px")),
          box(title="IL-23 Dynamics", width=6, status="purple",
              plotlyOutput("plt_il23", height="350px"))
        )
      ),

      ## ─── TAB 4: PASI ENDPOINTS ──────────────────────────────────────────
      tabItem("tab_pasi",
        fluidRow(
          box(title="PASI Score over Time", width=8, status="success",
              plotlyOutput("plt_pasi", height="400px")),
          box(title="PASI Response Rate", width=4, status="primary",
              plotlyOutput("plt_pasi_resp", height="400px"))
        ),
        fluidRow(
          box(title="PASI75/90/100 over Time", width=8, status="warning",
              plotlyOutput("plt_pasi75", height="300px")),
          box(title="Key Wk12 & Wk16 Endpoints", width=4, status="info",
              tableOutput("tbl_endpoints"))
        )
      ),

      ## ─── TAB 5: SCENARIO COMPARISON ────────────────────────────────────
      tabItem("tab_scenario",
        fluidRow(
          box(title="Run All 7 Scenarios", width=12, status="primary",
              actionButton("run_all","Run All Scenarios", class="btn-success btn-lg",
                           icon=icon("sync")),
              helpText("Compares: No Tx · Adalimumab · Secukinumab · Risankizumab · Apremilast · Tofacitinib · MTX")
          )
        ),
        fluidRow(
          box(title="PASI — All Scenarios", width=6, status="success",
              plotlyOutput("plt_all_pasi", height="350px")),
          box(title="IL-17A — All Scenarios", width=6, status="warning",
              plotlyOutput("plt_all_il17", height="350px"))
        ),
        fluidRow(
          box(title="PASI75/90/100 at Week 16", width=12, status="info",
              plotlyOutput("plt_bar_resp", height="350px"))
        )
      ),

      ## ─── TAB 6: BIOMARKER DASHBOARD ────────────────────────────────────
      tabItem("tab_bm",
        fluidRow(
          valueBoxOutput("bm_pasi75",  width=3),
          valueBoxOutput("bm_pasi90",  width=3),
          valueBoxOutput("bm_pasi100", width=3),
          valueBoxOutput("bm_iga01",   width=3)
        ),
        fluidRow(
          box(title="Serum Biomarker Heatmap (relative change from baseline)", width=12,
              status="purple", plotlyOutput("plt_bm_heatmap", height="350px"))
        ),
        fluidRow(
          box(title="Biomarker Time Series", width=6, status="info",
              plotlyOutput("plt_bm_ts", height="300px")),
          box(title="Clinical Endpoint Summary Table", width=6, status="success",
              DTOutput("tbl_bm_summary"))
        )
      )
    )
  )
)

## ─────────────────────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## Reactive: run single scenario
  sim_result <- eventReactive(input$run_sim, {
    sim_end <- input$sim_weeks * 7 * 24
    init_vals <- list(
      DC=2.5, IL23=120, Th17=input$base_th17,
      IL17A=input$base_il17, TNFa=input$base_tnf, IFNg=50,
      KC=input$base_pasi * 10, PASI=input$base_pasi,
      ADA_SC=0, ADA_C=0, ADA_P=0,
      SEC_SC=0, SEC_C=0, SEC_P=0,
      RSK_SC=0, RSK_C=0,
      APR_GI=0, APR_C=0, APR_P=0,
      TOF_GI=0, TOF_C=0,
      MTX_GI=0, MTX_C=0, MTX_PG=0
    )
    ev_obj <- build_events(input$drug, input$dose_mg, input$sim_weeks)
    pso_mod %>%
      init(!!!init_vals) %>%
      ev(ev_obj) %>%
      mrgsim(end=sim_end, delta=12) %>%
      as.data.frame() %>%
      mutate(week=time/168)
  })

  ## Reactive: run all scenarios
  all_sim <- eventReactive(input$run_all, {
    init_vals <- list(
      DC=2.5, IL23=120, Th17=350,
      IL17A=140, TNFa=80, IFNg=50,
      KC=200, PASI=20,
      ADA_SC=0, ADA_C=0, ADA_P=0,
      SEC_SC=0, SEC_C=0, SEC_P=0,
      RSK_SC=0, RSK_C=0,
      APR_GI=0, APR_C=0, APR_P=0,
      TOF_GI=0, TOF_C=0,
      MTX_GI=0, MTX_C=0, MTX_PG=0
    )
    weeks_t <- 52
    drugs   <- c("none","adalimumab","secukinumab","risankizumab",
                 "apremilast","tofacitinib","methotrexate")
    labs    <- c("No Treatment","Adalimumab","Secukinumab","Risankizumab",
                 "Apremilast","Tofacitinib","Methotrexate")
    bind_rows(lapply(seq_along(drugs), function(i) {
      ev_obj <- build_events(drugs[i], 30, weeks_t)
      pso_mod %>%
        init(!!!init_vals) %>%
        ev(ev_obj) %>%
        mrgsim(end=weeks_t*7*24, delta=12) %>%
        as.data.frame() %>%
        mutate(week=time/168, scenario=labs[i])
    }))
  })

  ## Value boxes
  output$vb_pasi <- renderValueBox(
    valueBox(input$base_pasi, "Baseline PASI", icon=icon("ruler"), color="purple")
  )
  output$vb_sev <- renderValueBox({
    sev <- if(input$base_pasi < 10) "Mild" else if(input$base_pasi < 20) "Moderate" else "Severe"
    valueBox(sev, "Disease Severity", icon=icon("heartbeat"), color="orange")
  })
  output$vb_drug <- renderValueBox(
    valueBox(input$drug, "Selected Drug", icon=icon("pills"), color="blue")
  )

  ## Patient summary table
  output$tbl_patient_summary <- renderTable({
    data.frame(
      Parameter = c("Baseline PASI","Baseline IL-17A","Baseline TNF-α",
                    "Baseline Th17","Disease Type","Treatment"),
      Value     = c(input$base_pasi, paste0(input$base_il17," pg/mL"),
                    paste0(input$base_tnf," pg/mL"),
                    paste0(input$base_th17," ×1000/mL"),
                    input$psoriasis_type, input$drug)
    )
  })

  ## PK plots
  output$plt_bio_pk <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>%
      filter(time %% 168 == 0) %>%
      pivot_longer(c(ADA_nM, SEC_nM, RSK_nM), names_to="drug", values_to="Cp_nM") %>%
      mutate(drug=recode(drug, ADA_nM="Adalimumab",SEC_nM="Secukinumab",RSK_nM="Risankizumab"))
    p <- ggplot(d, aes(week, Cp_nM, color=drug)) + geom_line(size=1.1) +
      labs(x="Week",y="Cp (nM)",color="Biologic") + theme_bw()
    ggplotly(p)
  })

  output$plt_sm_pk <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>%
      filter(time %% 168 == 0) %>%
      pivot_longer(c(APR_nM, TOF_nM), names_to="drug", values_to="Cp_nM") %>%
      mutate(drug=recode(drug, APR_nM="Apremilast",TOF_nM="Tofacitinib"))
    p <- ggplot(d, aes(week, Cp_nM, color=drug)) + geom_line(size=1.1) +
      labs(x="Week",y="Cp (nM)",color="Drug") + theme_bw()
    ggplotly(p)
  })

  ## PK parameters table
  output$tbl_pk_params <- renderDT({
    data.frame(
      Drug=c("Adalimumab","Secukinumab","Risankizumab","Ustekinumab","Apremilast","Tofacitinib","Methotrexate"),
      Mechanism=c("anti-TNF-α","anti-IL-17A","anti-IL-23p19","anti-IL-12/23p40","PDE4 inhibitor","JAK1/3 inhibitor","DHFR inhibitor"),
      Route=c("SC","SC","SC","SC","Oral","Oral","Oral"),
      Regimen=c("40mg q2w","300mg wk0-4, q4w","150mg wk0,4, q12w","45mg wk0,4, q12w","30mg BID","10mg BID","20mg qweek"),
      F_pct=c("64","73","89","57","73","74","70"),
      t_half_d=c("14","27","28","21-23","9h","3h","3-10h"),
      CL_L_h=c("0.247","0.191","0.078","0.252","9.5","22.8","4.8")
    )
  }, options=list(pageLength=7, scrollX=TRUE))

  ## PD plots
  output$plt_il17 <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>% filter(time %% 168 == 0)
    p <- ggplot(d, aes(week, IL17A)) + geom_line(color="#2196F3",size=1.2) +
      geom_hline(yintercept=30, linetype="dashed", color="gray60") +
      labs(title="Serum IL-17A",x="Week",y="IL-17A (pg/mL)") + theme_bw() +
      annotate("text",x=2,y=32,label="Normal ~30 pg/mL",size=3,color="gray40")
    ggplotly(p)
  })

  output$plt_tnf <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>% filter(time %% 168 == 0)
    p <- ggplot(d, aes(week, TNFa)) + geom_line(color="#F44336",size=1.2) +
      labs(title="Serum TNF-α",x="Week",y="TNF-α (pg/mL)") + theme_bw()
    ggplotly(p)
  })

  output$plt_th17 <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>% filter(time %% 168 == 0)
    p <- ggplot(d, aes(week, Th17)) + geom_line(color="#9C27B0",size=1.2) +
      labs(title="Th17 Cell Count",x="Week",y="Th17 (×1000/mL)") + theme_bw()
    ggplotly(p)
  })

  output$plt_il23 <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>% filter(time %% 168 == 0)
    p <- ggplot(d, aes(week, IL23)) + geom_line(color="#FF9800",size=1.2) +
      labs(title="Serum IL-23",x="Week",y="IL-23 (pg/mL)") + theme_bw()
    ggplotly(p)
  })

  ## PASI plots
  output$plt_pasi <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>% filter(time %% 168 == 0)
    p <- ggplot(d, aes(week, PASI)) + geom_line(color="#4CAF50",size=1.3) +
      geom_hline(yintercept=c(0.25*input$base_pasi, 0.10*input$base_pasi),
                 linetype="dashed", color=c("#FF9800","#F44336")) +
      labs(title="PASI Score", x="Week", y="PASI") + theme_bw()
    ggplotly(p)
  })

  output$plt_pasi_resp <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>%
      filter(week %in% c(12, 16, 24, 52)) %>%
      group_by(week) %>%
      summarize(PASI75=mean(PASI75)*100, PASI90=mean(PASI90)*100, PASI100=mean(PASI100)*100) %>%
      pivot_longer(-week, names_to="endpoint", values_to="pct")
    p <- ggplot(d, aes(factor(week), pct, fill=endpoint)) +
      geom_col(position="dodge") +
      scale_fill_manual(values=c(PASI75="#4CAF50",PASI90="#2196F3",PASI100="#9C27B0")) +
      labs(title="Response Rate",x="Week",y="% Responders",fill="") + theme_bw()
    ggplotly(p)
  })

  output$plt_pasi75 <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>%
      filter(time %% 168 == 0) %>%
      pivot_longer(c(PASI75,PASI90,PASI100), names_to="resp", values_to="val") %>%
      group_by(week, resp) %>%
      summarize(pct=mean(val)*100, .groups="drop")
    p <- ggplot(d, aes(week, pct, color=resp)) + geom_line(size=1.1) +
      scale_color_manual(values=c(PASI75="#4CAF50",PASI90="#2196F3",PASI100="#9C27B0")) +
      labs(x="Week",y="% Achieving Response",color="") + theme_bw()
    ggplotly(p)
  })

  output$tbl_endpoints <- renderTable({
    req(sim_result())
    sim_result() %>%
      filter(week %in% c(12, 16, 52)) %>%
      group_by(week) %>%
      summarize(
        `PASI (mean)` = round(mean(PASI),1),
        `PASI75 (%)` = round(mean(PASI75)*100,1),
        `PASI90 (%)` = round(mean(PASI90)*100,1),
        `IGA 0/1 (%)` = round(mean(IGA01)*100,1)
      )
  })

  ## Scenario comparison plots
  scen_colors <- c("No Treatment"="#616161","Adalimumab"="#F44336",
                   "Secukinumab"="#2196F3","Risankizumab"="#4CAF50",
                   "Apremilast"="#FF9800","Tofacitinib"="#9C27B0",
                   "Methotrexate"="#795548")

  output$plt_all_pasi <- renderPlotly({
    req(all_sim())
    d <- all_sim() %>% filter(time %% 168 == 0)
    p <- ggplot(d, aes(week, PASI, color=scenario)) + geom_line(size=1.1) +
      scale_color_manual(values=scen_colors) +
      labs(title="PASI — All Scenarios",x="Week",y="PASI",color="") + theme_bw()
    ggplotly(p)
  })

  output$plt_all_il17 <- renderPlotly({
    req(all_sim())
    d <- all_sim() %>% filter(time %% 168 == 0)
    p <- ggplot(d, aes(week, IL17A, color=scenario)) + geom_line(size=1.1) +
      scale_color_manual(values=scen_colors) +
      labs(title="IL-17A — All Scenarios",x="Week",y="IL-17A (pg/mL)",color="") + theme_bw()
    ggplotly(p)
  })

  output$plt_bar_resp <- renderPlotly({
    req(all_sim())
    d <- all_sim() %>%
      filter(abs(week - 16) < 0.1) %>%
      group_by(scenario) %>%
      summarize(PASI75=mean(PASI75)*100, PASI90=mean(PASI90)*100, PASI100=mean(PASI100)*100) %>%
      pivot_longer(-scenario, names_to="endpoint", values_to="pct")
    p <- ggplot(d, aes(scenario, pct, fill=endpoint)) +
      geom_col(position="dodge") +
      scale_fill_manual(values=c(PASI75="#4CAF50",PASI90="#2196F3",PASI100="#9C27B0")) +
      labs(title="PASI Response at Week 16",x="",y="Responders (%)",fill="") +
      theme_bw() + theme(axis.text.x=element_text(angle=35, hjust=1))
    ggplotly(p)
  })

  ## Biomarker value boxes
  output$bm_pasi75 <- renderValueBox({
    req(sim_result())
    val <- sim_result() %>% filter(abs(week-16)<0.1) %>%
      summarize(v=mean(PASI75)*100) %>% pull(v)
    valueBox(paste0(round(val,1),"%"), "PASI75 at Wk16", icon=icon("check"), color="green")
  })
  output$bm_pasi90 <- renderValueBox({
    req(sim_result())
    val <- sim_result() %>% filter(abs(week-16)<0.1) %>%
      summarize(v=mean(PASI90)*100) %>% pull(v)
    valueBox(paste0(round(val,1),"%"), "PASI90 at Wk16", icon=icon("check-double"), color="blue")
  })
  output$bm_pasi100 <- renderValueBox({
    req(sim_result())
    val <- sim_result() %>% filter(abs(week-16)<0.1) %>%
      summarize(v=mean(PASI100)*100) %>% pull(v)
    valueBox(paste0(round(val,1),"%"), "PASI100 (clear) Wk16", icon=icon("star"), color="purple")
  })
  output$bm_iga01 <- renderValueBox({
    req(sim_result())
    val <- sim_result() %>% filter(abs(week-16)<0.1) %>%
      summarize(v=mean(IGA01)*100) %>% pull(v)
    valueBox(paste0(round(val,1),"%"), "IGA 0/1 at Wk16", icon=icon("thumbs-up"), color="yellow")
  })

  ## Biomarker heatmap
  output$plt_bm_heatmap <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>%
      filter(week %in% c(0, 4, 8, 12, 16, 24, 36, 52)) %>%
      group_by(week) %>%
      summarize(IL17A=mean(IL17A), TNFa=mean(TNFa), Th17=mean(Th17),
                IL23=mean(IL23), KC=mean(KC)) %>%
      mutate(across(c(IL17A,TNFa,Th17,IL23,KC), ~round(./.[1]*100-100, 1))) %>%
      pivot_longer(-week, names_to="biomarker", values_to="pct_chg")
    p <- ggplot(d, aes(factor(week), biomarker, fill=pct_chg)) +
      geom_tile(color="white") +
      geom_text(aes(label=paste0(round(pct_chg),"%")), size=3) +
      scale_fill_gradient2(low="#2196F3", mid="white", high="#F44336", midpoint=0) +
      labs(title="Biomarker % Change from Baseline",x="Week",y="",fill="% Change") +
      theme_bw()
    ggplotly(p)
  })

  output$plt_bm_ts <- renderPlotly({
    req(sim_result())
    d <- sim_result() %>%
      filter(time %% 168 == 0) %>%
      pivot_longer(c(IL17A, TNFa, Th17), names_to="bm", values_to="val")
    p <- ggplot(d, aes(week, val, color=bm)) + geom_line(size=1.1) +
      scale_color_manual(values=c(IL17A="#2196F3",TNFa="#F44336",Th17="#9C27B0")) +
      labs(x="Week",y="Value",color="Biomarker") + theme_bw()
    ggplotly(p)
  })

  output$tbl_bm_summary <- renderDT({
    req(sim_result())
    sim_result() %>%
      filter(week %in% c(0, 12, 16, 52)) %>%
      group_by(week) %>%
      summarize(
        `PASI`=round(mean(PASI),1), `IL-17A`=round(mean(IL17A),1),
        `TNF-α`=round(mean(TNFa),1), `Th17`=round(mean(Th17),0),
        `IL-23`=round(mean(IL23),1), `KC idx`=round(mean(KC),0)
      ) %>% ungroup()
  }, options=list(pageLength=5))
}

shinyApp(ui, server)
