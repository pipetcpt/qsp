##########################################################################
# IgA Nephropathy (IgAN) QSP – Interactive Shiny Dashboard
#
# Tabs:
#   1. Patient Profile & Disease Staging
#   2. Disease Pathophysiology (biomarker cascade)
#   3. Drug PK Profiles
#   4. Proteinuria (UPCR) Endpoints
#   5. eGFR & Renal Function
#   6. Mechanistic Biomarkers
#   7. Scenario Comparison Table
#   8. Oxford MEST-C Visualisation
##########################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

# =========================================================================
# Inline mrgsolve model (same as igan_mrgsolve_model.R, trimmed for Shiny)
# =========================================================================
igan_code <- '
$PARAM
k_syn_IgA1=0.0578 k_deg_IgA1=0.0578
k_syn_AB=0.0198   k_deg_AB=0.0198
k_form_IC=0.10    k_clear_IC=0.12
k_syn_CP=0.09     k_deg_CP=0.18
k_act_MES=0.07    k_res_MES=0.04
k_inj_Pod=0.0045  k_rep_Pod=0.0018
k_syn_TIF=0.0028  k_deg_TIF=0.0004
k_syn_UPCR=2.0    k_deg_UPCR=0.28
k_loss_GFR=0.00018 k_RAAS_GFR=0.00014 eGFR_0=75
BP_0=135 k_BP_on=0.045 k_BP_off=0.045
ka_BUD=0.5  CL_BUD=90  V_BUD=200  F_BUD=0.15
Emax_BUD=0.55 EC50_BUD=8  Hill_BUD=1.5
ka_SPA=0.4  CL_SPA=14  V_SPA=120  F_SPA=0.85
Emax_SPA=0.62 EC50_SPA=200 Hill_SPA=1.2
ka_IPT=0.8  CL_IPT=8   V_IPT=60   F_IPT=0.90
Emax_IPT=0.88 EC50_IPT=50  Hill_IPT=2.0
ka_SIB=0.012 CL_SIB=0.008 V1_SIB=3.5 V2_SIB=2.5 Q_SIB=0.02 F_SIB=0.80
Emax_SIB=0.68 EC50_SIB=5   Hill_SIB=1.0
E_RAAS=0 Emax_RAAS_P=0.35 Emax_RAAS_G=0.22

$CMT
BUD_gut BUD_central SPA_gut SPA_central IPT_gut IPT_central
SIB_depot SIB_central SIB_periph
GdIgA1 AutoAb IC_mes CompAP Mesangial Podocyte TIF UPCR eGFR BP_sys

$MAIN
double C_BUD = BUD_central/V_BUD;
double C_SPA = SPA_central/V_SPA;
double C_IPT = IPT_central/V_IPT;
double C_SIB = SIB_central/V1_SIB;
double E_bud = Emax_BUD*pow(C_BUD,Hill_BUD)/(pow(EC50_BUD,Hill_BUD)+pow(C_BUD,Hill_BUD)+1e-12);
double E_spa = Emax_SPA*pow(C_SPA,Hill_SPA)/(pow(EC50_SPA,Hill_SPA)+pow(C_SPA,Hill_SPA)+1e-12);
double E_ipt = Emax_IPT*pow(C_IPT,Hill_IPT)/(pow(EC50_IPT,Hill_IPT)+pow(C_IPT,Hill_IPT)+1e-12);
double E_sib = Emax_SIB*pow(C_SIB,Hill_SIB)/(pow(EC50_SIB,Hill_SIB)+pow(C_SIB,Hill_SIB)+1e-12);
double E_mucosal = 1.0-(1.0-E_bud)*(1.0-E_sib);
double E_raas_p  = Emax_RAAS_P*E_RAAS;
double E_raas_g  = Emax_RAAS_G*E_RAAS;
double pod_ok    = (Podocyte>0)?Podocyte:0.0;
double tif_cap   = (TIF<1)?TIF:1.0;
double egfr_ok   = (eGFR>5)?eGFR:5.0;
double inj_drv   = CompAP*0.55+Mesangial*0.45;

$ODE
dxdt_BUD_gut     = -ka_BUD*BUD_gut;
dxdt_BUD_central =  ka_BUD*F_BUD*BUD_gut-(CL_BUD/V_BUD)*BUD_central;
dxdt_SPA_gut     = -ka_SPA*SPA_gut;
dxdt_SPA_central =  ka_SPA*F_SPA*SPA_gut-(CL_SPA/V_SPA)*SPA_central;
dxdt_IPT_gut     = -ka_IPT*IPT_gut;
dxdt_IPT_central =  ka_IPT*F_IPT*IPT_gut-(CL_IPT/V_IPT)*IPT_central;
dxdt_SIB_depot   = -ka_SIB*SIB_depot;
dxdt_SIB_central =  ka_SIB*F_SIB*SIB_depot-(CL_SIB/V1_SIB)*SIB_central
                    -(Q_SIB/V1_SIB)*SIB_central+(Q_SIB/V2_SIB)*SIB_periph;
dxdt_SIB_periph  =  (Q_SIB/V1_SIB)*SIB_central-(Q_SIB/V2_SIB)*SIB_periph;
dxdt_GdIgA1   = k_syn_IgA1*(1.0-E_mucosal)-k_deg_IgA1*GdIgA1;
dxdt_AutoAb   = k_syn_AB*GdIgA1-k_deg_AB*AutoAb;
dxdt_IC_mes   = k_form_IC*GdIgA1*AutoAb-k_clear_IC*IC_mes;
dxdt_CompAP   = k_syn_CP*IC_mes*(1.0-E_ipt)-k_deg_CP*CompAP;
dxdt_Mesangial= k_act_MES*(CompAP+IC_mes*0.3)-k_res_MES*Mesangial;
dxdt_Podocyte = -k_inj_Pod*inj_drv*pod_ok+k_rep_Pod*(1.0-pod_ok);
dxdt_TIF      = k_syn_TIF*(1.0-pod_ok)*(1.0+UPCR/3.0)*(1.0-tif_cap)-k_deg_TIF*tif_cap;
double UPCR_s = k_syn_UPCR*(1.0-pod_ok)*(1.0+Mesangial*0.4)*(1.0-E_spa)*(1.0-E_raas_p);
dxdt_UPCR     = UPCR_s-k_deg_UPCR*UPCR;
double GFR_loss=(k_loss_GFR*tif_cap+k_RAAS_GFR*(1.0-E_raas_g)*(1.0-pod_ok));
dxdt_eGFR     = -GFR_loss*egfr_ok;
dxdt_BP_sys   = k_BP_on*(( 1.0+Mesangial*0.12-E_raas_p*0.45)*BP_0-BP_sys)-k_BP_off*(BP_sys-BP_0);

$TABLE
double C_BUD_t = BUD_central/V_BUD;
double C_SPA_t = SPA_central/V_SPA;
double C_IPT_t = IPT_central/V_IPT;
double C_SIB_t = SIB_central/V1_SIB;
double UPCR_pct = 100.0*(UPCR-2.5)/2.5;
double CR50     = (UPCR<=1.25)?1.0:0.0;

$CAPTURE C_BUD_t C_SPA_t C_IPT_t C_SIB_t
         GdIgA1 AutoAb IC_mes CompAP Mesangial Podocyte TIF
         UPCR eGFR BP_sys UPCR_pct CR50
'

mod_base <- mcode("IgAN_Shiny", igan_code, quiet = TRUE)

init_default <- list(
  BUD_gut=0,BUD_central=0,SPA_gut=0,SPA_central=0,
  IPT_gut=0,IPT_central=0,SIB_depot=0,SIB_central=0,SIB_periph=0,
  GdIgA1=1.55, AutoAb=1.30, IC_mes=0.85, CompAP=0.65,
  Mesangial=0.75, Podocyte=0.82, TIF=0.12,
  UPCR=2.50, eGFR=75.0, BP_sys=136.0
)

# =========================================================================
# UI
# =========================================================================
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "IgAN QSP Simulator"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "tab_patient",  icon = icon("user")),
      menuItem("Pathophysiology",    tabName = "tab_patho",    icon = icon("project-diagram")),
      menuItem("Drug PK",            tabName = "tab_pk",       icon = icon("pills")),
      menuItem("Proteinuria (UPCR)", tabName = "tab_upcr",     icon = icon("tint")),
      menuItem("eGFR & Renal Fn",    tabName = "tab_egfr",     icon = icon("kidneys")),
      menuItem("Biomarkers",         tabName = "tab_bio",      icon = icon("flask")),
      menuItem("Scenario Table",     tabName = "tab_table",    icon = icon("table")),
      menuItem("Oxford MEST-C",      tabName = "tab_mestc",    icon = icon("microscope"))
    ),

    hr(),
    h5("  Patient Parameters", style = "color:#ecf0f1; margin-left:12px"),
    sliderInput("init_eGFR",    "Baseline eGFR (mL/min)", 20, 90, 75, step = 5),
    sliderInput("init_UPCR",    "Baseline UPCR (g/g)",    0.5, 5.0, 2.5, step = 0.1),
    sliderInput("init_Podocyte","Podocyte integrity",     0.5, 1.0, 0.82, step = 0.02),
    sliderInput("init_GdIgA1",  "Baseline Gd-IgA1 (norm)",1.0, 2.5, 1.55, step = 0.05),
    sliderInput("sim_dur",      "Simulation (months)",     6, 36, 24, step = 6),

    hr(),
    h5("  Drug Selection", style = "color:#ecf0f1; margin-left:12px"),
    checkboxInput("use_RAAS",  "RAAS Inhibitor (ACEi/ARB)", TRUE),
    checkboxInput("use_BUD",   "Budesonide TRF (16 mg/d)", FALSE),
    checkboxInput("use_SPA",   "Sparsentan (400 mg/d)",    FALSE),
    checkboxInput("use_IPT",   "Iptacopan (200 mg BID)",   FALSE),
    checkboxInput("use_SIB",   "Sibeprenlimab (500 mg Q4W)", FALSE),

    hr(),
    actionButton("run_sim", "  Run Simulation",
                 icon = icon("play"), class = "btn-success btn-lg",
                 style = "width:90%; margin:5px 5%")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .skin-blue .main-header .logo { background-color:#1565C0; }
      .skin-blue .main-header .navbar { background-color:#1976D2; }
      .skin-blue .main-sidebar { background-color:#1a237e; }
      .content-wrapper { background-color:#F5F5F5; }
      .box-header { background-color:#E3F2FD; }
    "))),

    tabItems(

      # ---- Tab 1: Patient Profile ----
      tabItem("tab_patient",
        fluidRow(
          valueBoxOutput("vb_eGFR",  width = 3),
          valueBoxOutput("vb_UPCR",  width = 3),
          valueBoxOutput("vb_GdIgA1",width = 3),
          valueBoxOutput("vb_CKD",   width = 3)
        ),
        fluidRow(
          box(title = "KDIGO Risk Category", width = 6, status = "primary",
            plotOutput("plot_kdigo_heat", height = "280px")
          ),
          box(title = "Disease Staging Summary", width = 6, status = "info",
            tableOutput("tbl_staging")
          )
        ),
        fluidRow(
          box(title = "IgAN Four-Hit Cascade — Baseline Severity", width = 12, status = "warning",
            plotOutput("plot_baseline_cascade", height = "200px")
          )
        )
      ),

      # ---- Tab 2: Pathophysiology ----
      tabItem("tab_patho",
        fluidRow(
          box(title = "Four-Hit Cascade Over Time", width = 12, status = "primary",
            plotlyOutput("plot_cascade", height = "420px")
          )
        ),
        fluidRow(
          box(title = "Mesangial Activation & Complement", width = 6, status = "warning",
            plotlyOutput("plot_mes_comp", height = "320px")
          ),
          box(title = "Podocyte Integrity Over Time", width = 6, status = "danger",
            plotlyOutput("plot_pod", height = "320px")
          )
        )
      ),

      # ---- Tab 3: Drug PK ----
      tabItem("tab_pk",
        fluidRow(
          box(title = "Budesonide TRF — Plasma Concentration (first 14 days)",
              width = 6, status = "success",
              plotlyOutput("pk_bud", height = "280px")),
          box(title = "Sparsentan — Plasma Concentration (first 14 days)",
              width = 6, status = "primary",
              plotlyOutput("pk_spa", height = "280px"))
        ),
        fluidRow(
          box(title = "Iptacopan — Plasma Concentration (first 14 days)",
              width = 6, status = "purple",
              plotlyOutput("pk_ipt", height = "280px")),
          box(title = "Sibeprenlimab (mAb) — Serum Concentration (6 months)",
              width = 6, status = "info",
              plotlyOutput("pk_sib", height = "280px"))
        )
      ),

      # ---- Tab 4: UPCR ----
      tabItem("tab_upcr",
        fluidRow(
          box(title = "UPCR Trajectory — All Treatments", width = 8, status = "danger",
            plotlyOutput("plot_upcr", height = "380px")
          ),
          box(title = "Week-36 Waterfall (% UPCR Change)", width = 4, status = "warning",
            plotlyOutput("plot_waterfall", height = "380px")
          )
        ),
        fluidRow(
          box(title = "Dose–Response: Sparsentan UPCR at Wk 36", width = 6, status = "primary",
            plotlyOutput("plot_dr_spa", height = "280px")
          ),
          box(title = "CR50 Achievement Over Time", width = 6, status = "success",
            plotlyOutput("plot_cr50", height = "280px")
          )
        )
      ),

      # ---- Tab 5: eGFR ----
      tabItem("tab_egfr",
        fluidRow(
          box(title = "eGFR Trajectory (CKD stage lines annotated)", width = 8, status = "info",
            plotlyOutput("plot_egfr", height = "380px")
          ),
          box(title = "eGFR Slope (mL/min/yr) at 2 Years", width = 4, status = "primary",
            plotlyOutput("plot_egfr_slope", height = "380px")
          )
        ),
        fluidRow(
          box(title = "Systolic BP Over Time", width = 6, status = "warning",
            plotlyOutput("plot_bp", height = "260px")
          ),
          box(title = "TIF Index Over Time", width = 6, status = "danger",
            plotlyOutput("plot_tif", height = "260px")
          )
        )
      ),

      # ---- Tab 6: Biomarkers ----
      tabItem("tab_bio",
        fluidRow(
          box(title = "Serum Gd-IgA1 (Hit 1 Biomarker)", width = 6, status = "warning",
            plotlyOutput("plot_gdiga1", height = "280px")
          ),
          box(title = "Anti-Gd-IgA1 IgG (Hit 2 Autoantibody)", width = 6, status = "danger",
            plotlyOutput("plot_autoab", height = "280px")
          )
        ),
        fluidRow(
          box(title = "Mesangial IC Deposits (Hit 3+4)", width = 6, status = "info",
            plotlyOutput("plot_ic", height = "280px")
          ),
          box(title = "Complement AP Activity", width = 6, status = "primary",
            plotlyOutput("plot_comp", height = "280px")
          )
        )
      ),

      # ---- Tab 7: Scenario Table ----
      tabItem("tab_table",
        fluidRow(
          box(title = "Scenario Comparison — Key Endpoints at Wk 36 & Wk 104",
              width = 12, status = "primary",
              DTOutput("dt_summary"))
        )
      ),

      # ---- Tab 8: Oxford MEST-C ----
      tabItem("tab_mestc",
        fluidRow(
          box(title = "Oxford MEST-C Score — Proxy from Model State",
              width = 7, status = "danger",
              plotlyOutput("plot_mestc", height = "380px")
          ),
          box(title = "MEST-C Score Legend & ESKD Risk", width = 5, status = "info",
            HTML("
              <table class='table table-condensed'>
                <thead><tr><th>Score</th><th>Lesion</th><th>Model proxy</th></tr></thead>
                <tbody>
                  <tr><td><b>M0/M1</b></td><td>Mesangial hypercellularity</td><td>Mesangial index</td></tr>
                  <tr><td><b>E0/E1</b></td><td>Endocapillary proliferation</td><td>CompAP × IC</td></tr>
                  <tr><td><b>S0/S1</b></td><td>Segmental glomerulosclerosis</td><td>1–Podocyte</td></tr>
                  <tr><td><b>T0/T1/T2</b></td><td>Tubular atrophy / TIF</td><td>TIF index</td></tr>
                  <tr><td><b>C0/C1/C2</b></td><td>Crescents</td><td>Mesangial × (1–Podocyte)</td></tr>
                </tbody>
              </table>
              <hr>
              <p><b>ESKD risk at 10 yr (KDIGO 2021):</b></p>
              <ul>
                <li>UPCR &lt; 1 g/g: ~5%</li>
                <li>UPCR 1–3 g/g: ~15–25%</li>
                <li>UPCR &gt; 3 g/g: ~40–60%</li>
                <li>T2 lesion: +25% additive risk</li>
              </ul>
            ")
          )
        )
      )
    )
  )
)

# =========================================================================
# SERVER
# =========================================================================
server <- function(input, output, session) {

  # ---- Reactive simulation ----
  sim_data <- eventReactive(input$run_sim, ignoreNULL = FALSE, {
    init <- modifyList(init_default, list(
      eGFR    = input$init_eGFR,
      UPCR    = input$init_UPCR,
      Podocyte= input$init_Podocyte,
      GdIgA1  = input$init_GdIgA1
    ))

    SIM_DAYS <- input$sim_dur * 30.44

    make_events <- function() {
      evs <- c()
      if (input$use_BUD) evs <- c(evs, list(ev(amt=16,  ii=24, addl=SIM_DAYS-1, cmt="BUD_gut")))
      if (input$use_SPA) evs <- c(evs, list(ev(amt=400, ii=24, addl=SIM_DAYS-1, cmt="SPA_gut")))
      if (input$use_IPT) evs <- c(evs, list(ev(amt=200, ii=12, addl=2*SIM_DAYS-1, cmt="IPT_gut")))
      if (input$use_SIB) evs <- c(evs, list(ev(amt=500, ii=28*24, addl=25, cmt="SIB_depot")))
      if (length(evs) == 0) return(NULL)
      do.call(c, evs)
    }

    raas_val <- ifelse(input$use_RAAS, 1.0, 0.0)
    m <- mod_base %>% init(init) %>% param(E_RAAS = raas_val)
    evs <- make_events()

    if (is.null(evs)) {
      out <- m %>% mrgsim(end = SIM_DAYS, delta = 1)
    } else {
      out <- m %>% mrgsim(events = evs, end = SIM_DAYS, delta = 1)
    }

    as_tibble(out) %>%
      mutate(time_mo = time / 30.44,
             time_wk = time / 7)
  })

  # All 7 scenarios for comparison
  all_scenarios <- eventReactive(input$run_sim, ignoreNULL = FALSE, {
    init <- modifyList(init_default, list(
      eGFR    = input$init_eGFR,
      UPCR    = input$init_UPCR,
      Podocyte= input$init_Podocyte,
      GdIgA1  = input$init_GdIgA1
    ))
    SIM_DAYS <- input$sim_dur * 30.44

    e_BUD <- ev(amt=16,  ii=24,    addl=SIM_DAYS-1,   cmt="BUD_gut")
    e_SPA <- ev(amt=400, ii=24,    addl=SIM_DAYS-1,   cmt="SPA_gut")
    e_IPT <- ev(amt=200, ii=12,    addl=2*SIM_DAYS-1, cmt="IPT_gut")
    e_SIB <- ev(amt=500, ii=28*24, addl=25,            cmt="SIB_depot")

    run_s <- function(evs, raas, lbl) {
      m <- mod_base %>% init(init) %>% param(E_RAAS = raas)
      out <- if (is.null(evs)) mrgsim(m, end=SIM_DAYS, delta=1)
             else               mrgsim(m, events=evs, end=SIM_DAYS, delta=1)
      as_tibble(out) %>% mutate(Scenario=lbl, time_mo=time/30.44, time_wk=time/7)
    }

    bind_rows(
      run_s(NULL,               0, "Untreated"),
      run_s(NULL,               1, "RAAS inh."),
      run_s(e_BUD,              1, "RAAS+BUD"),
      run_s(e_SPA,              0, "Sparsentan"),
      run_s(e_IPT,              1, "RAAS+IPT"),
      run_s(e_SIB,              1, "RAAS+SIB"),
      run_s(c(e_BUD,e_SPA,e_IPT),0, "Triple Combo")
    )
  })

  pal7 <- c(
    "Untreated"    = "#B71C1C", "RAAS inh."  = "#E65100",
    "RAAS+BUD"     = "#1B5E20", "Sparsentan" = "#0D47A1",
    "RAAS+IPT"     = "#4A148C", "RAAS+SIB"  = "#006064",
    "Triple Combo" = "#212121"
  )

  theme_shiny <- theme_bw(base_size=12) +
    theme(legend.position="bottom", legend.text=element_text(size=8))

  # ---- Value Boxes ----
  output$vb_eGFR <- renderValueBox({
    v <- input$init_eGFR
    color <- if (v >= 60) "green" else if (v >= 30) "yellow" else "red"
    valueBox(paste0(v, " mL/min"), "Baseline eGFR", icon=icon("tachometer-alt"), color=color)
  })
  output$vb_UPCR <- renderValueBox({
    v <- input$init_UPCR
    color <- if (v < 1) "green" else if (v < 3) "yellow" else "red"
    valueBox(paste0(v, " g/g"), "Baseline UPCR", icon=icon("tint"), color=color)
  })
  output$vb_GdIgA1 <- renderValueBox({
    valueBox(input$init_GdIgA1, "Gd-IgA1 (norm.)", icon=icon("dna"), color="orange")
  })
  output$vb_CKD <- renderValueBox({
    v <- input$init_eGFR
    stage <- if (v >= 90) "G1" else if (v >= 60) "G2" else if (v >= 45) "G3a" else if (v >= 30) "G3b" else if (v >= 15) "G4" else "G5"
    valueBox(stage, "CKD Stage", icon=icon("kidneys"), color="blue")
  })

  output$tbl_staging <- renderTable({
    data.frame(
      Parameter     = c("eGFR","UPCR","Gd-IgA1","Podocyte integrity"),
      `At Baseline` = c(
        paste0(input$init_eGFR, " mL/min/1.73m²"),
        paste0(input$init_UPCR, " g/g"),
        paste0(input$init_GdIgA1, " (norm.)"),
        paste0(round(input$init_Podocyte*100,0), "%")
      ),
      `Risk flag` = c(
        ifelse(input$init_eGFR < 60, "⚠ CKD ≥G3", "OK"),
        ifelse(input$init_UPCR >= 1, "⚠ High-risk", "Low-risk"),
        ifelse(input$init_GdIgA1 > 1.4, "⚠ Elevated", "Normal"),
        ifelse(input$init_Podocyte < 0.9, "⚠ Impaired", "Normal")
      )
    )
  }, rownames=FALSE)

  output$plot_baseline_cascade <- renderPlot({
    df <- data.frame(
      Hit = c("Gd-IgA1\n(Hit 1)","AutoAb\n(Hit 2)","IC Deposit\n(Hit 3+4)","CompAP","Mesangial","Podocyte\n(injured)"),
      Value = c(
        (input$init_GdIgA1-1)/1.5,
        0.30,
        0.85 * (input$init_GdIgA1-1)/1.5,
        0.65 * (input$init_GdIgA1-1)/1.5,
        0.75 * (input$init_GdIgA1-1)/1.5,
        1 - input$init_Podocyte
      )
    )
    ggplot(df, aes(x=reorder(Hit,Value), y=Value, fill=Value)) +
      geom_col(show.legend=FALSE) +
      coord_flip() +
      scale_fill_gradient(low="#FFF9C4", high="#B71C1C") +
      labs(x=NULL, y="Severity index (0–1)") +
      theme_minimal(base_size=12)
  })

  output$plot_kdigo_heat <- renderPlot({
    df <- expand.grid(
      eGFR_cat = c("G1\n(≥90)","G2\n(60-89)","G3a\n(45-59)","G3b\n(30-44)","G4\n(15-29)","G5\n(<15)"),
      UPCR_cat = c("A1\n(<0.3)","A2\n(0.3-1)","A3\n(≥1)")
    )
    df$risk <- c(1,1,1,2,2,3, 1,1,2,2,3,3, 1,2,3,3,3,3)
    ggplot(df, aes(x=UPCR_cat, y=eGFR_cat, fill=factor(risk))) +
      geom_tile(color="white", linewidth=1.5) +
      scale_fill_manual(values=c("1"="#4CAF50","2"="#FFC107","3"="#F44336"),
                        labels=c("Low","High","Very High"),name="ESKD Risk") +
      labs(x="Proteinuria (A category)", y="eGFR (G category)",
           title="KDIGO 2021 Risk Matrix") +
      theme_minimal(base_size=12)
  })

  # ---- Pathophysiology tab ----
  output$plot_cascade <- renderPlotly({
    d <- sim_data()
    df <- pivot_longer(d, cols=c("GdIgA1","AutoAb","IC_mes","CompAP"),
                       names_to="variable", values_to="value")
    p <- ggplot(df, aes(x=time_mo, y=value, color=variable)) +
      geom_line(linewidth=1.1) +
      labs(title="Four-Hit Cascade Biomarkers",
           x="Time (months)", y="Level (normalized)", color="") +
      theme_shiny
    ggplotly(p)
  })

  output$plot_mes_comp <- renderPlotly({
    d <- sim_data()
    df <- pivot_longer(d, cols=c("Mesangial","CompAP"),
                       names_to="variable", values_to="value")
    p <- ggplot(df, aes(x=time_mo, y=value, color=variable)) +
      geom_line(linewidth=1.1) +
      labs(x="Time (months)", y="Index (norm.)", color="") + theme_shiny
    ggplotly(p)
  })

  output$plot_pod <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=time_mo, y=Podocyte)) +
      geom_line(color="#B71C1C", linewidth=1.3) +
      geom_hline(yintercept=0.7, linetype="dashed", color="grey40") +
      labs(x="Time (months)", y="Podocyte integrity (0–1)") + theme_shiny
    ggplotly(p)
  })

  # ---- PK tabs ----
  pk_plot <- function(data_col, color, title, ylab, max_t=14) {
    renderPlotly({
      d <- sim_data() %>% filter(time <= max_t)
      p <- ggplot(d, aes_string(x="time", y=data_col)) +
        geom_line(color=color, linewidth=1.2) +
        labs(title=title, x="Time (days)", y=ylab) + theme_shiny
      ggplotly(p)
    })
  }
  output$pk_bud <- pk_plot("C_BUD_t","#1B5E20","Budesonide TRF","Conc. (ng/mL)")
  output$pk_spa <- pk_plot("C_SPA_t","#0D47A1","Sparsentan","Conc. (ng/mL)")
  output$pk_ipt <- pk_plot("C_IPT_t","#4A148C","Iptacopan","Conc. (ng/mL)")
  output$pk_sib <- pk_plot("C_SIB_t","#006064","Sibeprenlimab (SC)","Serum (μg/mL)", max_t=180)

  # ---- UPCR tab ----
  output$plot_upcr <- renderPlotly({
    d <- all_scenarios()
    p <- ggplot(d, aes(x=time_mo, y=UPCR, color=Scenario)) +
      geom_line(linewidth=1.0) +
      scale_color_manual(values=pal7) +
      geom_hline(yintercept=1.25, linetype="dashed", color="grey40") +
      annotate("text", x=1, y=1.15, label="CR50 (1.25 g/g)", size=3, hjust=0) +
      labs(x="Time (months)", y="UPCR (g/g)") + theme_shiny
    ggplotly(p)
  })

  output$plot_waterfall <- renderPlotly({
    d <- all_scenarios() %>%
      group_by(Scenario) %>%
      filter(time == max(time[time_wk <= 36])) %>%
      summarise(pct_chg = UPCR_pct[1], .groups="drop") %>%
      arrange(pct_chg)
    p <- ggplot(d, aes(x=reorder(Scenario,pct_chg), y=pct_chg, fill=pct_chg<0)) +
      geom_col(show.legend=FALSE) +
      coord_flip() +
      scale_fill_manual(values=c("TRUE"="#1B5E20","FALSE"="#B71C1C")) +
      geom_hline(yintercept=-50, linetype="dashed", color="black") +
      labs(x=NULL, y="% UPCR change vs baseline") + theme_shiny
    ggplotly(p)
  })

  output$plot_dr_spa <- renderPlotly({
    init <- modifyList(init_default, list(
      eGFR=input$init_eGFR, UPCR=input$init_UPCR, GdIgA1=input$init_GdIgA1))
    doses <- c(25,50,100,200,400,800,1200)
    res <- lapply(doses, function(d) {
      e <- ev(amt=d, ii=24, addl=251, cmt="SPA_gut")
      mod_base %>% init(init) %>% mrgsim(events=e, end=252, delta=1) %>%
        as_tibble() %>% filter(time==252) %>% mutate(Dose=d)
    }) %>% bind_rows()
    p <- ggplot(res, aes(x=Dose, y=UPCR)) +
      geom_line(color="#0D47A1",linewidth=1.2) + geom_point(size=3,color="#0D47A1") +
      scale_x_log10() +
      labs(x="Sparsentan dose (mg)", y="UPCR at Wk 36 (g/g)") + theme_shiny
    ggplotly(p)
  })

  output$plot_cr50 <- renderPlotly({
    d <- all_scenarios()
    p <- ggplot(d, aes(x=time_mo, y=CR50*100, color=Scenario)) +
      geom_line(linewidth=1.0) +
      scale_color_manual(values=pal7) +
      labs(x="Time (months)", y="CR50 Achieved (0/100%)") + theme_shiny
    ggplotly(p)
  })

  # ---- eGFR tab ----
  output$plot_egfr <- renderPlotly({
    d <- all_scenarios()
    p <- ggplot(d, aes(x=time_mo, y=eGFR, color=Scenario)) +
      geom_line(linewidth=1.0) +
      scale_color_manual(values=pal7) +
      geom_hline(yintercept=c(60,30,15), linetype="dashed", color="grey50") +
      labs(x="Time (months)", y="eGFR (mL/min/1.73 m²)") + theme_shiny
    ggplotly(p)
  })

  output$plot_egfr_slope <- renderPlotly({
    d <- all_scenarios() %>%
      filter(time_mo >= input$sim_dur * 0.8) %>%
      group_by(Scenario) %>%
      summarise(
        eGFR_final = mean(eGFR), eGFR_init = input$init_eGFR,
        months = input$sim_dur,
        slope = (eGFR_final - eGFR_init) / (months / 12),
        .groups="drop"
      )
    p <- ggplot(d, aes(x=reorder(Scenario,slope), y=slope, fill=slope>-3)) +
      geom_col(show.legend=FALSE) + coord_flip() +
      scale_fill_manual(values=c("TRUE"="#1B5E20","FALSE"="#B71C1C")) +
      labs(x=NULL, y="eGFR slope (mL/min/yr)") + theme_shiny
    ggplotly(p)
  })

  output$plot_bp <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=time_mo, y=BP_sys)) +
      geom_line(color="#00695C",linewidth=1.2) +
      geom_hline(yintercept=130, linetype="dashed", color="grey40") +
      annotate("text", x=0.5, y=131, label="Target BP 130 mmHg", size=3, hjust=0) +
      labs(x="Time (months)", y="Systolic BP (mmHg)") + theme_shiny
    ggplotly(p)
  })

  output$plot_tif <- renderPlotly({
    d <- all_scenarios()
    p <- ggplot(d, aes(x=time_mo, y=TIF, color=Scenario)) +
      geom_line(linewidth=1.0) + scale_color_manual(values=pal7) +
      labs(x="Time (months)", y="TIF Index (0–1)") + theme_shiny
    ggplotly(p)
  })

  # ---- Biomarker tab ----
  mk_bio_plot <- function(col, color, title) {
    renderPlotly({
      d <- all_scenarios()
      p <- ggplot(d, aes_string(x="time_mo", y=col, color="Scenario")) +
        geom_line(linewidth=0.9) + scale_color_manual(values=pal7) +
        labs(title=title, x="Time (months)", y=col) + theme_shiny
      ggplotly(p)
    })
  }
  output$plot_gdiga1 <- mk_bio_plot("GdIgA1","#E65100","Serum Gd-IgA1 (Hit 1)")
  output$plot_autoab <- mk_bio_plot("AutoAb","#880E4F","Anti-Gd-IgA1 IgG (Hit 2)")
  output$plot_ic     <- mk_bio_plot("IC_mes","#6A1B9A","Mesangial IC Deposit (Hit 3+4)")
  output$plot_comp   <- mk_bio_plot("CompAP","#1565C0","Complement AP Activity")

  # ---- Summary Table ----
  output$dt_summary <- renderDT({
    d <- all_scenarios()
    tps <- c(36*7, 52*7, min(input$sim_dur*30.44, 730))
    res <- lapply(tps, function(tp) {
      d %>%
        group_by(Scenario) %>%
        filter(time == max(time[time <= tp])) %>%
        summarise(
          Timepoint = paste0("Mo ", round(tp/30.44,0)),
          UPCR      = round(mean(UPCR),2),
          pct_UPCR  = round(mean(UPCR_pct),1),
          eGFR      = round(mean(eGFR),1),
          GdIgA1    = round(mean(GdIgA1),2),
          CompAP    = round(mean(CompAP),2),
          Podocyte  = round(mean(Podocyte),3),
          TIF       = round(mean(TIF),3),
          CR50      = round(mean(CR50),0),
          .groups="drop"
        )
    }) %>% bind_rows()
    datatable(res, options=list(pageLength=20, scrollX=TRUE),
              rownames=FALSE) %>%
      formatStyle("pct_UPCR",
        background=styleInterval(c(-50, 0), c("#C8E6C9","#FFF9C4","#FFCDD2"))) %>%
      formatStyle("CR50",
        backgroundColor=styleEqual(c(0,1), c("white","#C8E6C9")))
  })

  # ---- MEST-C tab ----
  output$plot_mestc <- renderPlotly({
    d <- sim_data() %>%
      mutate(
        M_score = pmin(Mesangial / 0.5, 1),
        E_score = pmin(CompAP * IC_mes / 0.4, 1),
        S_score = 1 - Podocyte,
        T_score = TIF,
        C_score = pmin(Mesangial * (1-Podocyte) / 0.3, 1)
      ) %>%
      pivot_longer(cols=c(M_score,E_score,S_score,T_score,C_score),
                   names_to="Score", values_to="Value")
    p <- ggplot(d, aes(x=time_mo, y=Value, color=Score)) +
      geom_line(linewidth=1.1) +
      scale_color_brewer(palette="Set1") +
      labs(x="Time (months)", y="Score index (0–1)") + theme_shiny
    ggplotly(p)
  })
}

# =========================================================================
# Launch
# =========================================================================
shinyApp(ui, server)
