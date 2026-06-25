## ============================================================
##  Aplastic Anemia QSP — Interactive Shiny Dashboard
##  6 Tabs: Patient Profile · PK · BM & Immunity ·
##          Clinical Endpoints · Scenario Comparison · Biomarkers
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(DT)

## ---- Embed model code ----------------------------------------
model_code <- '
$PARAM
  kCD8p=0.03, kCD8d=0.03, kTRp=3.0, kTRd=0.03,
  kCYp=0.5, kCYd=2.0, kHSCsr=0.08, kHSCd=0.02, kHSCdf=0.05,
  kCYkill=0.25, kPROGdf=0.15, kCFUEdf=0.20,
  kRETmat=0.14, kRBCd=0.0083, kMKPdf=0.15, kPLTd=0.10,
  kNEUPdf=0.33, kNEUd=2.0, kTELOs=0.0005,
  kATGel=0.116, kATG12=0.5, kATG21=0.2, kATGkill=0.008,
  kCSAabs=1.5, kCSAel=0.35, kCSA12=0.8, kCSA21=0.4,
  IC50CSA=150, VCSA=5.0,
  kEPGabs=0.55, kEPGel=0.77,
  EmaxHSC=0.8, EC50HSC=2.0, EmaxPLT=0.9, EC50PLT=1.5, EmaxTR=0.3,
  hATG_on=0, rATG_on=0, CSA_on=0, EPAG_on=0, GCSF_on=0,
  CD8_init=2000, TREG_init=20, HSC_sev=5,
  ATG_amt=2800, CSA_amt=350, EPAG_amt=150

$CMT CD8 TREG CYTO HSC PROG CFU_E RETIC RBC MKP PLT NEUP NEU TELO
     CATG PATG CGUT CCSA EGUT CEPG

$MAIN
  CD8_0  = CD8_init;
  TREG_0 = TREG_init;

$INIT
  CD8=2000, TREG=20, CYTO=3.0, HSC=5, PROG=5, CFU_E=5,
  RETIC=5, RBC=7.0, MKP=5, PLT=10, NEUP=3, NEU=0.35,
  TELO=0.90, CATG=0, PATG=0, CGUT=0, CCSA=0, EGUT=0, CEPG=0

$ODE
  double kATGel_eff = hATG_on*kATGel + rATG_on*(kATGel*0.23);
  double ATG_Tkill  = kATGkill * CATG;
  double CsA_inh    = CSA_on * CCSA / (CCSA + IC50CSA);
  double EPAG_HSCeff= EPAG_on * EmaxHSC * CEPG / (CEPG + EC50HSC);
  double EPAG_PLTeff= EPAG_on * EmaxPLT * CEPG / (CEPG + EC50PLT);
  double EPAG_TReff = EPAG_on * EmaxTR  * CEPG / (CEPG + EC50PLT);
  double GCSF_boost = GCSF_on * 3.0;
  double TReg_sup   = 1.0 / (1.0 + (TREG / 50.0));
  double CYTO_prod  = kCYp * (CD8/2000) * TReg_sup * (1 - CsA_inh*0.7);

  dxdt_CD8  = kCD8p*CD8*(1.0-CD8/10000) - kCD8d*CD8 - ATG_Tkill*CD8 - CsA_inh*0.4*kCD8p*CD8;
  dxdt_TREG = kTRp - kTRd*TREG - ATG_Tkill*0.7*TREG + EPAG_TReff*kTRp;
  dxdt_CYTO = CYTO_prod - kCYd*CYTO;
  dxdt_HSC  = kHSCsr*HSC*(1.0-HSC/100)*TELO*(1+EPAG_HSCeff) - kHSCd*HSC - kHSCdf*HSC - kCYkill*CYTO*HSC;
  dxdt_PROG = kHSCdf*HSC*10 - kPROGdf*PROG;
  dxdt_CFU_E= kPROGdf*PROG*0.35 - kCFUEdf*CFU_E;
  dxdt_RETIC= kCFUEdf*CFU_E*2  - kRETmat*RETIC;
  dxdt_RBC  = kRETmat*RETIC*0.07 - kRBCd*RBC;
  dxdt_MKP  = kPROGdf*PROG*0.20*(1+EPAG_PLTeff) - kMKPdf*MKP;
  dxdt_PLT  = kMKPdf*MKP*1.2 - kPLTd*PLT;
  dxdt_NEUP = kPROGdf*PROG*0.45*(1+GCSF_boost) - kNEUPdf*NEUP;
  dxdt_NEU  = kNEUPdf*NEUP*0.04 - kNEUd*NEU;
  dxdt_TELO = -kTELOs*(kHSCsr*HSC/100);
  dxdt_CATG = -(kATGel_eff+kATG12)*CATG + kATG21*PATG;
  dxdt_PATG = kATG12*CATG - kATG21*PATG;
  dxdt_CGUT = -kCSAabs*CGUT;
  dxdt_CCSA = kCSAabs*CGUT/VCSA - (kCSAel+kCSA12)*CCSA + kCSA21*(CCSA*kCSA12/kCSA21);
  dxdt_EGUT = -kEPGabs*EGUT;
  dxdt_CEPG = kEPGabs*EGUT - kEPGel*CEPG;

$TABLE
  double ANC  = NEU;
  double HGB  = RBC;
  double PLTc = PLT;
  double CR   = (ANC>1.0 && PLTc>100 && HGB>10) ? 1.0 : 0.0;
  double PR   = (ANC>0.5 && PLTc>20  && HGB>8 && CR==0) ? 1.0 : 0.0;
  capture ANC HGB PLTc CR PR HSC CYTO CD8 TREG TELO CATG CCSA CEPG
'

mod <- mcode("AA_Shiny", model_code, quiet=TRUE)

## ---- Helper: run simulation ----------------------------------
run_sim <- function(severity, hATG, rATG, csa, epag, gcsf,
                    atg_amt, csa_amt, epag_amt, duration) {
  # Initial conditions by severity
  cd8_init <- switch(severity,
    "sAA"  = 2000, "vsAA" = 2500, "nsAA" = 1200, "Normal" = 800)
  hsc_init <- switch(severity,
    "sAA"  = 5,    "vsAA" = 2,    "nsAA" = 20,   "Normal" = 100)

  m2 <- mod %>%
    init(CD8=cd8_init, HSC=hsc_init,
         PLT = ifelse(severity=="Normal",200,ifelse(severity=="nsAA",50,10)),
         NEU = ifelse(severity=="Normal",3.0, ifelse(severity=="nsAA",0.8,0.35)),
         RBC = ifelse(severity=="Normal",14,  ifelse(severity=="nsAA",9.0,7.0)),
         TREG= ifelse(severity=="Normal",100, 20)) %>%
    param(hATG_on=hATG, rATG_on=rATG, CSA_on=csa, EPAG_on=epag, GCSF_on=gcsf)

  ev_list <- list()
  if (hATG == 1) ev_list[[length(ev_list)+1]] <- ev(time=0:3,  cmt="CATG", amt=atg_amt,  rate=-2)
  if (rATG == 1) ev_list[[length(ev_list)+1]] <- ev(time=0:4,  cmt="CATG", amt=atg_amt*0.088, rate=-2)
  if (csa  == 1) ev_list[[length(ev_list)+1]] <- ev(time=0, cmt="CGUT", amt=csa_amt, addl=pmin(duration-1,179), ii=1)
  if (epag == 1) ev_list[[length(ev_list)+1]] <- ev(time=0, cmt="EGUT", amt=epag_amt, addl=duration-1, ii=1)

  if (length(ev_list) == 0) {
    mrgsim(m2, end=duration, delta=1) %>% as.data.frame()
  } else {
    ev_combined <- Reduce("+", ev_list)
    mrgsim(m2, events=ev_combined, end=duration, delta=1) %>% as.data.frame()
  }
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title="Aplastic Anemia QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName="tab_patient",   icon=icon("user-md")),
      menuItem("Pharmacokinetics",   tabName="tab_pk",        icon=icon("flask")),
      menuItem("BM & Immunity",      tabName="tab_bm",        icon=icon("dna")),
      menuItem("Clinical Endpoints", tabName="tab_clinical",  icon=icon("heartbeat")),
      menuItem("Scenario Comparison",tabName="tab_scenario",  icon=icon("chart-bar")),
      menuItem("Biomarkers",         tabName="tab_biomarker", icon=icon("vials"))
    ),
    hr(),
    h4("Treatment Settings", style="color:white;margin-left:15px;"),
    selectInput("severity", "Disease Severity:",
                choices=c("Severe AA (sAA)"="sAA","Very Severe AA (vsAA)"="vsAA",
                          "Non-Severe AA (nsAA)"="nsAA","Normal (reference)"="Normal")),
    checkboxInput("hATG",  "Horse ATG (hATG, 40 mg/kg×4d)", value=FALSE),
    checkboxInput("rATG",  "Rabbit ATG (rATG, 3.5 mg/kg×5d)", value=FALSE),
    checkboxInput("csa",   "Cyclosporine A (CsA)", value=FALSE),
    checkboxInput("epag",  "Eltrombopag (EPAG)", value=FALSE),
    checkboxInput("gcsf",  "G-CSF (Filgrastim)", value=FALSE),
    hr(),
    sliderInput("atg_amt",  "ATG dose (mg/admin):", min=500, max=5000, value=2800, step=100),
    sliderInput("csa_amt",  "CsA dose (mg/day):",   min=100, max=600,  value=350,  step=25),
    sliderInput("epag_amt", "EPAG dose (mg/day):",  min=25,  max=200,  value=150,  step=25),
    sliderInput("duration", "Simulation (days):",   min=90,  max=730,  value=365,  step=30),
    actionButton("run", "Run Simulation", icon=icon("play"),
                 style="color:#fff;background-color:#2196F3;border-color:#1565C0;width:90%;margin:5px 10px;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .skin-blue .main-header .logo { background-color:#1a237e; }
      .skin-blue .main-header .navbar { background-color:#283593; }
      .skin-blue .main-sidebar { background-color:#1a237e; }
      .value-box .inner h3 { font-size: 22px; }
    "))),

    tabItems(

      ## ---- TAB 1: Patient Profile --------------------------
      tabItem("tab_patient",
        h2("Patient Profile & Disease Classification"),
        fluidRow(
          valueBoxOutput("vb_anc",  width=3),
          valueBoxOutput("vb_hgb",  width=3),
          valueBoxOutput("vb_plt",  width=3),
          valueBoxOutput("vb_bm",   width=3)
        ),
        fluidRow(
          box(title="AA Diagnosis Criteria", width=6, solidHeader=TRUE, status="primary",
            tableOutput("diag_table")),
          box(title="Treatment Decision Framework", width=6, solidHeader=TRUE, status="warning",
            tableOutput("tx_table"))
        ),
        fluidRow(
          box(title="Response Classification at 6 Months", width=12, solidHeader=TRUE, status="success",
            plotlyOutput("response_gauge", height="250px"))
        )
      ),

      ## ---- TAB 2: Pharmacokinetics --------------------------
      tabItem("tab_pk",
        h2("Drug PK Profiles"),
        fluidRow(
          box(title="ATG Concentration–Time (Central Compartment)", width=12,
              solidHeader=TRUE, status="primary",
              plotlyOutput("pk_atg", height="280px"))
        ),
        fluidRow(
          box(title="Cyclosporine A Trough (Cmin)", width=6,
              solidHeader=TRUE, status="warning",
              plotlyOutput("pk_csa", height="280px")),
          box(title="Eltrombopag Plasma Concentration", width=6,
              solidHeader=TRUE, status="success",
              plotlyOutput("pk_epag", height="280px"))
        ),
        fluidRow(
          box(title="PK Summary at Key Timepoints", width=12, solidHeader=TRUE, status="info",
              DTOutput("pk_table"))
        )
      ),

      ## ---- TAB 3: BM & Immunity -----------------------------
      tabItem("tab_bm",
        h2("Bone Marrow & Immune Dynamics"),
        fluidRow(
          box(title="HSC Pool Recovery (% of Normal)", width=6,
              solidHeader=TRUE, status="primary",
              plotlyOutput("plot_hsc", height="300px")),
          box(title="Cytokine Index (IFN-γ/TNF-α)", width=6,
              solidHeader=TRUE, status="danger",
              plotlyOutput("plot_cyto", height="300px"))
        ),
        fluidRow(
          box(title="CD8+ T-Cell Dynamics", width=6,
              solidHeader=TRUE, status="warning",
              plotlyOutput("plot_cd8", height="300px")),
          box(title="Regulatory T-Cell (Treg)", width=6,
              solidHeader=TRUE, status="success",
              plotlyOutput("plot_treg", height="300px"))
        )
      ),

      ## ---- TAB 4: Clinical Endpoints -------------------------
      tabItem("tab_clinical",
        h2("Clinical Endpoints"),
        fluidRow(
          box(title="Absolute Neutrophil Count (ANC)", width=4,
              solidHeader=TRUE, status="danger",
              plotlyOutput("plot_anc", height="280px"),
              p("Threshold: CR >1.0 × 10⁹/L", style="color:red;font-size:11px;")),
          box(title="Hemoglobin (g/dL)", width=4,
              solidHeader=TRUE, status="warning",
              plotlyOutput("plot_hgb", height="280px"),
              p("Threshold: CR >10 g/dL", style="color:red;font-size:11px;")),
          box(title="Platelet Count (× 10⁹/L)", width=4,
              solidHeader=TRUE, status="primary",
              plotlyOutput("plot_plt", height="280px"),
              p("Threshold: CR >100 × 10⁹/L", style="color:red;font-size:11px;"))
        ),
        fluidRow(
          box(title="Response Status Over Time", width=12, solidHeader=TRUE, status="success",
              plotlyOutput("plot_response", height="250px"))
        ),
        fluidRow(
          box(title="Clinical Data Table", width=12, solidHeader=TRUE, status="info",
              DTOutput("clin_table"))
        )
      ),

      ## ---- TAB 5: Scenario Comparison -----------------------
      tabItem("tab_scenario",
        h2("Treatment Scenario Comparison"),
        fluidRow(
          box(title="Scenario Legend", width=12, status="info",
            p("All 6 predefined scenarios are simulated simultaneously for comparison."),
            tags$ul(
              tags$li("Scenario 1: Untreated severe AA (natural history)"),
              tags$li("Scenario 2: hATG + CsA (standard IST; Young 2012 NEJM target CR ~68%)"),
              tags$li("Scenario 3: rATG + CsA (inferior IST; Young 2012 NEJM target CR ~37%)"),
              tags$li("Scenario 4: hATG + CsA + EPAG (Townsley 2017 NEJM target CR ~58%)"),
              tags$li("Scenario 5: EPAG Monotherapy (Desmond 2013 JCI, ~44% response)"),
              tags$li("Scenario 6: CsA Monotherapy (non-severe or elderly)")
            )
          )
        ),
        fluidRow(
          box(title="ANC Comparison", width=6, solidHeader=TRUE, status="danger",
              plotlyOutput("sc_anc", height="300px")),
          box(title="Platelet Comparison", width=6, solidHeader=TRUE, status="primary",
              plotlyOutput("sc_plt", height="300px"))
        ),
        fluidRow(
          box(title="Hemoglobin Comparison", width=6, solidHeader=TRUE, status="warning",
              plotlyOutput("sc_hgb", height="300px")),
          box(title="HSC Recovery Comparison", width=6, solidHeader=TRUE, status="success",
              plotlyOutput("sc_hsc", height="300px"))
        ),
        fluidRow(
          box(title="6-Month Response Summary Table", width=12, solidHeader=TRUE, status="info",
              DTOutput("sc_table"))
        )
      ),

      ## ---- TAB 6: Biomarkers --------------------------------
      tabItem("tab_biomarker",
        h2("Biomarkers & Disease Progression"),
        fluidRow(
          box(title="Telomere Length Dynamics", width=6,
              solidHeader=TRUE, status="info",
              plotlyOutput("bm_telo", height="280px"),
              p("Short telomeres: poor IST response; consider Danazol if TL <10th percentile",
                style="color:blue;font-size:11px;"))  ,
          box(title="CD8:Treg Ratio (Immune Imbalance)", width=6,
              solidHeader=TRUE, status="warning",
              plotlyOutput("bm_ratio", height="280px"))
        ),
        fluidRow(
          box(title="Cytokine vs HSC (Phase Plot)", width=6,
              solidHeader=TRUE, status="danger",
              plotlyOutput("bm_phase", height="300px")),
          box(title="Biomarker Summary", width=6, solidHeader=TRUE, status="success",
              DTOutput("bm_table"))
        ),
        fluidRow(
          box(title="Infection & Bleeding Risk Indices", width=12,
              solidHeader=TRUE, status="primary",
              plotlyOutput("bm_risk", height="250px"))
        )
      )
    ) # end tabItems
  ) # end dashboardBody
) # end dashboardPage

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

  sim_data <- eventReactive(input$run, {
    run_sim(input$severity,
            as.integer(input$hATG), as.integer(input$rATG),
            as.integer(input$csa),  as.integer(input$epag),
            as.integer(input$gcsf),
            input$atg_amt, input$csa_amt, input$epag_amt,
            input$duration)
  }, ignoreNone=FALSE)

  # Pre-compute all 6 scenarios (for comparison tab)
  sc_data <- reactive({
    sc_params <- list(
      list(h=0,r=0,c=0,e=0,g=0, label="1: Untreated"),
      list(h=1,r=0,c=1,e=0,g=0, label="2: hATG+CsA"),
      list(h=0,r=1,c=1,e=0,g=0, label="3: rATG+CsA"),
      list(h=1,r=0,c=1,e=1,g=0, label="4: hATG+CsA+EPAG"),
      list(h=0,r=0,c=0,e=1,g=0, label="5: EPAG mono"),
      list(h=0,r=0,c=1,e=0,g=0, label="6: CsA mono")
    )
    bind_rows(lapply(sc_params, function(sc) {
      run_sim("sAA", sc$h, sc$r, sc$c, sc$e, sc$g,
              2800, 350, 150, 365) %>%
        mutate(Scenario=sc$label)
    }))
  })

  ## --- Value boxes ---
  last_vals <- reactive({
    d <- sim_data()
    d[nrow(d)-1,]
  })

  output$vb_anc <- renderValueBox({
    v <- round(sim_data() %>% filter(time==180) %>% pull(ANC), 2)
    valueBox(value=paste(v,"×10⁹/L"), subtitle="ANC at 6 mo",
             icon=icon("shield-alt"),
             color=ifelse(v>1.0,"green",ifelse(v>0.5,"yellow","red")))
  })
  output$vb_hgb <- renderValueBox({
    v <- round(sim_data() %>% filter(time==180) %>% pull(HGB), 1)
    valueBox(value=paste(v,"g/dL"), subtitle="Hemoglobin at 6 mo",
             icon=icon("tint"), color=ifelse(v>10,"green",ifelse(v>8,"yellow","red")))
  })
  output$vb_plt <- renderValueBox({
    v <- round(sim_data() %>% filter(time==180) %>% pull(PLTc), 0)
    valueBox(value=paste(v,"×10⁹/L"), subtitle="Platelets at 6 mo",
             icon=icon("circle"), color=ifelse(v>100,"green",ifelse(v>20,"yellow","red")))
  })
  output$vb_bm <- renderValueBox({
    v <- round(sim_data() %>% filter(time==180) %>% pull(HSC), 1)
    valueBox(value=paste(v,"%"), subtitle="BM Cellularity at 6 mo",
             icon=icon("bone"), color=ifelse(v>50,"green",ifelse(v>25,"yellow","red")))
  })

  ## --- Diagnostic table ---
  output$diag_table <- renderTable({
    data.frame(
      Criterion      = c("ANC threshold","PLT threshold","BM Cellularity","Reticulocytes"),
      `Non-Severe`   = c(">0.5×10⁹/L",  ">20×10⁹/L", "<35%",  "Variable"),
      `Severe (sAA)` = c("<0.5×10⁹/L",  "<20×10⁹/L", "<25%",  "<20×10⁹/L"),
      `Very Severe`  = c("<0.2×10⁹/L",  "<20×10⁹/L", "<25%",  "<20×10⁹/L")
    )
  }, striped=TRUE, bordered=TRUE)

  output$tx_table <- renderTable({
    data.frame(
      Category        = c("1st line sAA <40y","1st line sAA ≥40y","Refractory AA","Telomere disease"),
      Treatment       = c("Allo-HSCT (MSD)","hATG+CsA+EPAG","EPAG or 2nd IST","Danazol + IST"),
      Reference       = c("EBMT Guidelines 2022","Townsley 2017 NEJM","Desmond 2013 JCI","Townsley 2016 NEJM")
    )
  }, striped=TRUE, bordered=TRUE)

  output$response_gauge <- renderPlotly({
    d <- sim_data() %>% filter(time==180)
    cr <- d$CR; pr <- d$PR; nr <- d$NR
    plot_ly(type="bar", x=c("CR","PR","NR"), y=c(cr,pr,nr)*100,
            marker=list(color=c("#4CAF50","#FFC107","#F44336"))) %>%
      layout(title="Response at 6 Months (%)", yaxis=list(title="%"), xaxis=list(title="Response"))
  })

  ## --- PK plots ---
  output$pk_atg <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~CATG, type="scatter", mode="lines",
            line=list(color="#3F51B5",width=2)) %>%
      layout(title="ATG Concentration (μg/mL)", xaxis=list(title="Day"), yaxis=list(title="μg/mL"))
  })
  output$pk_csa <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~CCSA, type="scatter", mode="lines",
            line=list(color="#FF9800",width=2)) %>%
      add_hlines(y=list(100,200), line=list(dash="dot",color="red")) %>%
      layout(title="CsA Trough (ng/mL)", xaxis=list(title="Day"), yaxis=list(title="ng/mL"))
  })
  output$pk_epag <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~CEPG, type="scatter", mode="lines",
            line=list(color="#4CAF50",width=2)) %>%
      layout(title="Eltrombopag (μg/mL)", xaxis=list(title="Day"), yaxis=list(title="μg/mL"))
  })
  output$pk_table <- renderDT({
    d <- sim_data() %>% filter(time %in% c(1,4,7,14,28,90,180,270,365)) %>%
      select(time, CATG, CCSA, CEPG) %>%
      rename(Day=time, `ATG(μg/mL)`=CATG, `CsA Cmin(ng/mL)`=CCSA, `EPAG(μg/mL)`=CEPG) %>%
      mutate(across(where(is.numeric), ~round(.,3)))
    datatable(d, options=list(pageLength=10))
  })

  ## --- BM & Immunity ---
  output$plot_hsc <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~HSC, type="scatter", mode="lines",
            line=list(color="#1565C0",width=2)) %>%
      add_lines(y=rep(25,nrow(d)), x=d$time, line=list(dash="dot",color="red"), name="sAA threshold") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="% Normal"))
  })
  output$plot_cyto <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~CYTO, type="scatter", mode="lines",
            line=list(color="#D32F2F",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Cytokine Index"))
  })
  output$plot_cd8 <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~CD8, type="scatter", mode="lines",
            line=list(color="#F57F17",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="cells/μL"))
  })
  output$plot_treg <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~TREG, type="scatter", mode="lines",
            line=list(color="#388E3C",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="cells/μL"))
  })

  ## --- Clinical Endpoints ---
  output$plot_anc <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~ANC, type="scatter", mode="lines",
            line=list(color="#D32F2F",width=2)) %>%
      add_lines(y=rep(1.0,nrow(d)), x=d$time, line=list(dash="dot",color="darkred"), name="CR ≥1.0") %>%
      add_lines(y=rep(0.5,nrow(d)), x=d$time, line=list(dash="dash",color="red"), name="sAA <0.5") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="×10⁹/L"))
  })
  output$plot_hgb <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~HGB, type="scatter", mode="lines",
            line=list(color="#FF8F00",width=2)) %>%
      add_lines(y=rep(10,nrow(d)), x=d$time, line=list(dash="dot",color="darkorange"), name="CR ≥10") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="g/dL"))
  })
  output$plot_plt <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~PLTc, type="scatter", mode="lines",
            line=list(color="#1565C0",width=2)) %>%
      add_lines(y=rep(100,nrow(d)), x=d$time, line=list(dash="dot",color="darkblue"), name="CR ≥100") %>%
      add_lines(y=rep(20,nrow(d)),  x=d$time, line=list(dash="dash",color="blue"),    name="sAA <20") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="×10⁹/L"))
  })
  output$plot_response <- renderPlotly({
    d <- sim_data() %>%
      select(time, CR, PR) %>%
      mutate(NR=1-CR-PR) %>%
      pivot_longer(-time, names_to="Response", values_to="Flag") %>%
      filter(Flag==1)
    plot_ly(d, x=~time, y=~Response, type="scatter", mode="markers",
            color=~Response,
            colors=c("CR"="#4CAF50","PR"="#FFC107","NR"="#F44336")) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Response Status"))
  })
  output$clin_table <- renderDT({
    d <- sim_data() %>%
      filter(time %in% c(0,30,60,90,120,180,270,365)) %>%
      select(time, ANC, HGB, PLTc, HSC, CR, PR) %>%
      mutate(across(where(is.numeric), ~round(.,2)),
             Response = case_when(CR==1~"CR", PR==1~"PR", TRUE~"NR")) %>%
      select(-CR,-PR) %>%
      rename(Day=time, `ANC(×10⁹/L)`=ANC, `Hgb(g/dL)`=HGB,
             `PLT(×10⁹/L)`=PLTc, `BM(%)`=HSC)
    datatable(d, options=list(pageLength=8), rownames=FALSE)
  })

  ## --- Scenario Comparison ---
  output$sc_anc <- renderPlotly({
    d <- sc_data()
    plot_ly(d, x=~time, y=~ANC, color=~Scenario, type="scatter", mode="lines") %>%
      add_lines(y=rep(1.0,365), x=0:364, line=list(dash="dot",color="black"), name="CR threshold") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="ANC (×10⁹/L)"))
  })
  output$sc_plt <- renderPlotly({
    d <- sc_data()
    plot_ly(d, x=~time, y=~PLTc, color=~Scenario, type="scatter", mode="lines") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="PLT (×10⁹/L)"))
  })
  output$sc_hgb <- renderPlotly({
    d <- sc_data()
    plot_ly(d, x=~time, y=~HGB, color=~Scenario, type="scatter", mode="lines") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Hgb (g/dL)"))
  })
  output$sc_hsc <- renderPlotly({
    d <- sc_data()
    plot_ly(d, x=~time, y=~HSC, color=~Scenario, type="scatter", mode="lines") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="BM Cellularity (% normal)"))
  })
  output$sc_table <- renderDT({
    d <- sc_data() %>%
      filter(time==180) %>%
      group_by(Scenario) %>%
      summarise(
        ANC=round(mean(ANC),2), Hgb=round(mean(HGB),1),
        PLT=round(mean(PLTc),0), BM=round(mean(HSC),1),
        CR_pct=round(mean(CR)*100,0), PR_pct=round(mean(PR)*100,0)
      )
    datatable(d, rownames=FALSE, options=list(dom="t")) %>%
      formatStyle("CR_pct", background=styleColorBar(c(0,100),"#4CAF50"))
  })

  ## --- Biomarkers ---
  output$bm_telo <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time, y=~TELO, type="scatter", mode="lines",
            line=list(color="#00BCD4",width=2)) %>%
      add_lines(y=rep(0.75,nrow(d)), x=d$time, line=list(dash="dot",color="red"), name="Short TL threshold") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Telomere Length (relative)"))
  })
  output$bm_ratio <- renderPlotly({
    d <- sim_data() %>% mutate(ratio = CD8/pmax(TREG,1))
    plot_ly(d, x=~time, y=~ratio, type="scatter", mode="lines",
            line=list(color="#FF5722",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="CD8:Treg Ratio"))
  })
  output$bm_phase <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~CYTO, y=~HSC, type="scatter", mode="lines+markers",
            marker=list(size=4, color=~time, colorscale="Viridis", showscale=TRUE)) %>%
      layout(xaxis=list(title="Cytokine Index"), yaxis=list(title="HSC Pool (% normal)"))
  })
  output$bm_table <- renderDT({
    d <- sim_data() %>%
      filter(time %in% c(0,30,90,180,365)) %>%
      mutate(CD8_Treg_ratio = round(CD8/pmax(TREG,1),1)) %>%
      select(time, TELO, CD8, TREG, CD8_Treg_ratio, CYTO) %>%
      mutate(across(where(is.numeric), ~round(.,3))) %>%
      rename(Day=time, `TL (rel)`=TELO, `CD8(cells/μL)`=CD8,
             `Treg(cells/μL)`=TREG, `CD8:Treg`=CD8_Treg_ratio, Cytokine=CYTO)
    datatable(d, rownames=FALSE)
  })
  output$bm_risk <- renderPlotly({
    d <- sim_data() %>%
      mutate(
        Infect_risk = exp(-ANC/0.5) * 100,  # infection risk index (%)
        Bleed_risk  = exp(-PLTc/20) * 100   # bleeding risk index (%)
      )
    plot_ly(d, x=~time) %>%
      add_lines(y=~Infect_risk, name="Infection Risk Index", line=list(color="#F44336",width=2)) %>%
      add_lines(y=~Bleed_risk,  name="Bleeding Risk Index",  line=list(color="#9C27B0",width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Risk Index (%)"),
             title="Infection & Bleeding Risk Over Time")
  })

} # end server

shinyApp(ui, server)
