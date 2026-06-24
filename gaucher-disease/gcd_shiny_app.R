################################################################################
## Gaucher Disease QSP — Interactive Shiny Dashboard
## Tabs: Overview · Patient Profile · PK · Enzyme & Substrate ·
##       Organ/Haematology · Bone · Scenario Comparison · Biomarkers ·
##       Virtual Population
################################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ──────────────────────────────────────────────────────────
## Inline model (abbreviated — same ODEs as gcd_mrgsolve_model.R)
## ──────────────────────────────────────────────────────────
gcd_model_code <- '
$PARAM
BW=70,GBA0=5,GC0=100,GL1_0=8.5,LYSO_0=18,CHITR0=4800
SV0=2200,LV0=1.65,HGB0=10.5,PLT0=65,BMD0=-1.8
k_gba_deg=0.0035,k_gba_syn=0.0175
k_gc_syn=0.028,k_gc_deg=0.014
k_gl1_out=0.055,k_lyso_out=0.045,k_chitr_out=0.008
k_ferr_out=0.012,FERRIT0=220
k_sv_in=0.0005,k_sv_out=0.0003
k_lv_in=0.0004,k_lv_out=0.0002
k_hgb_prod=0.0018,k_hgb_deg=0.0017
k_plt_prod=0.018,k_plt_deg=0.013
k_bmd_form=0.0005,k_bmd_res=0.0008
k_oc_stim=0.0012,k_oc_deg=0.0010
k_ob_stim=0.0004,k_ob_deg=0.0008
k_il6_prod=0.025,k_il6_deg=0.030
k_nfkb_on=0.018,k_nfkb_off=0.015
CL_ERT=1.4,V1_ERT=0.18,Q_ERT=0.35,V2_ERT=0.55
EMAX_ERT=0.85,EC50_ERT=0.6,VELA_MOD=1.0
KA_ELIS=0.80,F_ELIS=0.20,CL_ELIS=38.0,V_ELIS=106
IC50_ELIS=0.010,IMAX_ELIS=0.95
KA_MIGS=0.60,F_MIGS=0.97,CL_MIGS=4.5,V_MIGS=28
IC50_MIGS=50.0,IMAX_MIGS=0.80
USE_ERT=0,USE_ELIS=0,USE_MIGS=0,DOSE_ERT=60

$CMT ERT_C ERT_T ELIS_GUT ELIS_C MIGS_GUT MIGS_C
     GBA GC_MAC GC_SP GC_LV GC_BM
     GL1 LYSOGL1 CHITR FERRIT SV LV HGB PLT BMD OC OB IL6 NFKB

$INIT
ERT_C=0,ERT_T=0,ELIS_GUT=0,ELIS_C=0,MIGS_GUT=0,MIGS_C=0
GBA=5,GC_MAC=100,GC_SP=100,GC_LV=100,GC_BM=100
GL1=8.5,LYSOGL1=18,CHITR=4800,FERRIT=220
SV=2200,LV=1.65,HGB=10.5,PLT=65,BMD=-1.8,OC=2.5,OB=0.9,IL6=12,NFKB=1.0

$ODE
double k10_ert=CL_ERT/V1_ERT;
double k12_ert=Q_ERT/V1_ERT;
double k21_ert=Q_ERT/V2_ERT;
dxdt_ERT_C=-(k10_ert+k12_ert)*ERT_C+k21_ert*ERT_T;
dxdt_ERT_T=k12_ert*ERT_C-k21_ert*ERT_T;
dxdt_ELIS_GUT=-KA_ELIS*ELIS_GUT;
dxdt_ELIS_C=KA_ELIS*ELIS_GUT*F_ELIS/V_ELIS-(CL_ELIS/V_ELIS)*ELIS_C;
dxdt_MIGS_GUT=-KA_MIGS*MIGS_GUT;
dxdt_MIGS_C=KA_MIGS*MIGS_GUT*F_MIGS/V_MIGS-(CL_MIGS/V_MIGS)*MIGS_C;
double ERT_eff=USE_ERT*EMAX_ERT*VELA_MOD*ERT_T/(EC50_ERT+ERT_T);
double ELIS_inh=USE_ELIS*IMAX_ELIS*ELIS_C/(IC50_ELIS+ELIS_C);
double MIGS_inh=USE_MIGS*IMAX_MIGS*MIGS_C/(IC50_MIGS+MIGS_C);
double SRT_inh=1.0-(1.0-ELIS_inh)*(1.0-MIGS_inh);
double GBA_eff=GBA+ERT_eff*(30.0-GBA);
dxdt_GBA=k_gba_syn-k_gba_deg*GBA+ERT_eff*0.005;
double GCS_act=k_gc_syn*(1.0-SRT_inh);
double GC_cl=(k_gc_deg*GBA_eff/GBA0+ERT_eff*0.025);
dxdt_GC_MAC=GCS_act-GC_cl*GC_MAC;
dxdt_GC_SP=0.15*GC_MAC-0.018*(GBA_eff/GBA0)*GC_SP;
dxdt_GC_LV=0.12*GC_MAC-0.016*(GBA_eff/GBA0)*GC_LV;
dxdt_GC_BM=0.18*GC_MAC-0.014*(GBA_eff/GBA0)*GC_BM;
double GC_exc=(GC_MAC>GC0)?(GC_MAC-GC0)/GC0:0.0;
dxdt_NFKB=k_nfkb_on*GC_exc-k_nfkb_off*NFKB;
dxdt_IL6=k_il6_prod*NFKB-k_il6_deg*IL6;
dxdt_GL1=0.045*GC_MAC-k_gl1_out*GL1;
dxdt_LYSOGL1=0.025*GC_MAC-k_lyso_out*LYSOGL1;
dxdt_CHITR=3.5*IL6-k_chitr_out*CHITR;
dxdt_FERRIT=0.55*IL6-k_ferr_out*FERRIT;
double sf=(SV>300)?300/SV:1.0;
double bmf=(GC_BM>50)?exp(-0.008*(GC_BM-50)):1.0;
dxdt_SV=k_sv_in*GC_SP-k_sv_out*(SV-300.0)*(GBA_eff/GBA0);
dxdt_LV=k_lv_in*GC_LV-k_lv_out*(LV-1.0)*(GBA_eff/GBA0);
dxdt_HGB=k_hgb_prod*bmf*sf-k_hgb_deg*HGB;
dxdt_PLT=k_plt_prod*bmf*sf-k_plt_deg*PLT;
dxdt_OC=k_oc_stim*IL6-k_oc_deg*OC;
dxdt_OB=k_ob_stim/(1.0+0.3*IL6)-k_ob_deg*OB;
dxdt_BMD=k_bmd_form*OB-k_bmd_res*OC;

$TABLE
double GBA_pct=100.0*GBA/30.0;
double SRT_GCS_inh=(USE_ELIS*IMAX_ELIS*ELIS_C/(IC50_ELIS+ELIS_C)+
                    USE_MIGS*IMAX_MIGS*MIGS_C/(IC50_MIGS+MIGS_C))*100.0;

$CAPTURE GBA GBA_pct GC_MAC GC_SP GC_LV GC_BM
         GL1 LYSOGL1 CHITR FERRIT SV LV HGB PLT BMD OC OB IL6 NFKB
         ERT_C ERT_T ELIS_C MIGS_C SRT_GCS_inh
'

mod_global <- mcode("gcd_shiny", gcd_model_code)

## ──────────────────────────────────────────────────────────
## UI
## ──────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "Gaucher Disease QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview",           tabName = "tab_overview",  icon = icon("info-circle")),
      menuItem("Patient Profile",    tabName = "tab_patient",   icon = icon("user")),
      menuItem("Drug PK",            tabName = "tab_pk",        icon = icon("pills")),
      menuItem("Enzyme & Substrate", tabName = "tab_enzyme",    icon = icon("flask")),
      menuItem("Organ / Haematology",tabName = "tab_organ",     icon = icon("heart")),
      menuItem("Bone",               tabName = "tab_bone",      icon = icon("bone")),
      menuItem("Scenario Comparison",tabName = "tab_scenario",  icon = icon("chart-line")),
      menuItem("Biomarkers",         tabName = "tab_bio",       icon = icon("vial")),
      menuItem("Virtual Population", tabName = "tab_vp",        icon = icon("users"))
    )
  ),

  dashboardBody(
    tabItems(

      ## ── TAB 1: OVERVIEW ──────────────────────────────────
      tabItem("tab_overview",
        fluidRow(
          box(title = "Gaucher Disease — QSP Model Overview", width = 12, status = "primary",
            p("Gaucher disease is the most common lysosomal storage disorder, caused by biallelic",
              "loss-of-function variants in", strong("GBA1"), "(glucocerebrosidase gene). Deficient",
              "lysosomal GBA leads to accumulation of glucocerebroside (GC) and its toxic deacylated",
              "form lyso-glucocerebroside (lyso-GL1/GL-1) primarily in tissue macrophages (Gaucher cells)."),
            p("This QSP model integrates:"),
            tags$ul(
              tags$li("Two-compartment ERT PK (imiglucerase / velaglucerase α)"),
              tags$li("One-compartment oral SRT PK (eliglustat, miglustat)"),
              tags$li("GBA enzyme kinetics and M6P-receptor mediated ERT delivery"),
              tags$li("GC substrate dynamics across macrophage, spleen, liver, bone marrow compartments"),
              tags$li("Inflammation cascade (NF-κB → IL-6 → cytokines)"),
              tags$li("Organ volume, hematology, and bone PD endpoints"),
              tags$li("Virtual patient population (Monte Carlo)")
            ),
            br(),
            tags$img(src = "gcd_qsp_model.png",
                     style = "max-width:100%; border:1px solid #ccc;")
          )
        )
      ),

      ## ── TAB 2: PATIENT PROFILE ───────────────────────────
      tabItem("tab_patient",
        fluidRow(
          box(title = "Patient Parameters", width = 4, status = "warning",
            sliderInput("bw",    "Body Weight (kg)",   30, 120, 70, 1),
            selectInput("type",  "Gaucher Type",
                        choices = c("Type 1 (non-neuronopathic)"=1,
                                    "Type 3 (chronic neuronopathic)"=3),
                        selected = 1),
            sliderInput("gba0",  "Baseline GBA Activity (nmol/h/mg; normal~30)", 0.5, 20, 5, 0.5),
            sliderInput("sv0",   "Baseline Spleen Volume (mL)",  300, 5000, 2200, 50),
            sliderInput("hgb0",  "Baseline Hemoglobin (g/dL)",  6, 15, 10.5, 0.1),
            sliderInput("plt0",  "Baseline Platelets (×10⁹/L)", 10, 200, 65, 5),
            sliderInput("bmd0",  "Baseline Lumbar T-score",     -4, 0, -1.8, 0.1),
            sliderInput("sim_yrs","Simulation Duration (years)",1, 5, 2, 0.5)
          ),
          box(title = "Treatment", width = 4, status = "primary",
            h4("ERT"),
            checkboxInput("use_ert",  "Enable ERT", FALSE),
            selectInput("ert_type", "ERT Agent",
                        choices = c("Imiglucerase"="imig","Velaglucerase α"="vela"),
                        selected = "imig"),
            sliderInput("dose_ert", "ERT Dose (U/kg Q2W)", 15, 120, 60, 15),
            hr(),
            h4("SRT"),
            checkboxInput("use_elis","Enable Eliglustat (SRT)", FALSE),
            selectInput("cyp2d6","CYP2D6 Status",
                        choices=c("Extensive Metabolizer"="em","Poor Metabolizer"="pm"),
                        selected="em"),
            checkboxInput("use_migs","Enable Miglustat (SRT)", FALSE),
            sliderInput("dose_migs","Miglustat Dose (mg TID)", 50, 200, 100, 25)
          ),
          box(title = "Disease Severity", width = 4, status = "danger",
            p("Gaucher Severity Score Index (GSSI) — estimated from parameters:"),
            tableOutput("gssi_table"),
            br(),
            p("Zimran Severity Score Index (SSI):"),
            tableOutput("ssi_table")
          )
        )
      ),

      ## ── TAB 3: DRUG PK ───────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(title = "ERT Plasma Concentration", width = 6, status = "primary",
            plotlyOutput("pk_ert_plot", height = "350px")
          ),
          box(title = "SRT Plasma Concentration", width = 6, status = "primary",
            plotlyOutput("pk_srt_plot", height = "350px")
          )
        ),
        fluidRow(
          box(title = "PK Parameters Summary", width = 12, status = "info",
            DTOutput("pk_table")
          )
        )
      ),

      ## ── TAB 4: ENZYME & SUBSTRATE ────────────────────────
      tabItem("tab_enzyme",
        fluidRow(
          box(title = "GBA Enzyme Activity (% of Normal)", width = 6, status = "warning",
            plotlyOutput("gba_plot", height = "320px")
          ),
          box(title = "Glucocerebroside Burden", width = 6, status = "danger",
            plotlyOutput("gc_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "GCS Inhibition (SRT)", width = 6, status = "primary",
            plotlyOutput("gcs_plot", height = "280px")
          ),
          box(title = "NF-κB & Inflammation", width = 6, status = "warning",
            plotlyOutput("nfkb_plot", height = "280px")
          )
        )
      ),

      ## ── TAB 5: ORGAN / HAEMATOLOGY ───────────────────────
      tabItem("tab_organ",
        fluidRow(
          box(title = "Spleen Volume (mL)", width = 6, status = "primary",
            plotlyOutput("sv_plot", height = "300px"),
            p("Therapeutic goal: <5× normal (normal ≈ 300 mL)")
          ),
          box(title = "Liver Volume (×normal)", width = 6, status = "primary",
            plotlyOutput("lv_plot", height = "300px"),
            p("Therapeutic goal: <1.25×normal")
          )
        ),
        fluidRow(
          box(title = "Hemoglobin (g/dL)", width = 6, status = "danger",
            plotlyOutput("hgb_plot", height = "300px"),
            p("Therapeutic goal: ≥11 g/dL (women) / ≥12 g/dL (men)")
          ),
          box(title = "Platelets (×10⁹/L)", width = 6, status = "danger",
            plotlyOutput("plt_plot", height = "300px"),
            p("Therapeutic goal: >100 ×10⁹/L (avoid splenectomy if possible)")
          )
        )
      ),

      ## ── TAB 6: BONE ──────────────────────────────────────
      tabItem("tab_bone",
        fluidRow(
          box(title = "Lumbar Spine T-score (BMD)", width = 6, status = "warning",
            plotlyOutput("bmd_plot", height = "320px")
          ),
          box(title = "OB / OC Balance", width = 6, status = "info",
            plotlyOutput("oboc_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Bone Outcomes Summary", width = 12, status = "primary",
            p("Osteonecrosis risk is driven by GC bone marrow burden and osteoclast hyperactivity."),
            p("ERT reduces GC_BM → lowers RANKL-driven OC → improves BMD over 2–5 years."),
            p("Clinical evidence: Wenstrup 2007 (Clin Orthop Relat Res) showed +5.8% lumbar BMD at 2 years with imiglucerase.")
          )
        )
      ),

      ## ── TAB 7: SCENARIO COMPARISON ───────────────────────
      tabItem("tab_scenario",
        fluidRow(
          box(title = "Select Variables to Compare", width = 4, status = "primary",
            checkboxGroupInput("sc_vars", "Endpoints to display:",
              choices = c("GL-1 (GlucoCeramide)"="GL1",
                          "Lyso-GL1"="LYSOGL1",
                          "Chitotriosidase"="CHITR",
                          "Spleen Volume"="SV",
                          "Hemoglobin"="HGB",
                          "Platelets"="PLT",
                          "BMD (T-score)"="BMD"),
              selected = c("GL1","SV","HGB"))
          ),
          box(title = "Scenario Comparison Plots", width = 8, status = "primary",
            plotOutput("scenario_plot", height = "500px")
          )
        ),
        fluidRow(
          box(title = "Response Table (12-Month, all scenarios)", width = 12, status = "info",
            DTOutput("resp_table")
          )
        )
      ),

      ## ── TAB 8: BIOMARKERS ────────────────────────────────
      tabItem("tab_bio",
        fluidRow(
          box(title = "Plasma GL-1 (Glucosylceramide, μg/L)", width = 6, status = "primary",
            plotlyOutput("gl1_plot", height = "300px"),
            p("Normal: <1 μg/L | Untreated Type 1: 5–15 μg/L")
          ),
          box(title = "Lyso-GL1 (Lyso-GlucoCeramide, ng/mL)", width = 6, status = "primary",
            plotlyOutput("lyso_plot", height = "300px"),
            p("Most sensitive marker; normal: <1.5 ng/mL")
          )
        ),
        fluidRow(
          box(title = "Chitotriosidase (nmol/h/mL)", width = 6, status = "warning",
            plotlyOutput("chitr_plot", height = "300px"),
            p("Macrophage activation marker; normal: <100 nmol/h/mL (absent in ~6% of population)")
          ),
          box(title = "Serum Ferritin (μg/L) & IL-6", width = 6, status = "warning",
            plotlyOutput("ferr_plot", height = "300px"),
            p("Inflammation markers; elevated in active Gaucher disease")
          )
        )
      ),

      ## ── TAB 9: VIRTUAL POPULATION ────────────────────────
      tabItem("tab_vp",
        fluidRow(
          box(title = "Virtual Population Settings", width = 4, status = "primary",
            sliderInput("n_vp",  "Number of Patients", 50, 500, 200, 50),
            selectInput("vp_scenario", "Treatment Scenario",
                        choices = c("Natural History"="nh",
                                    "Imiglucerase ERT"="ert",
                                    "Eliglustat SRT"="srt"),
                        selected = "ert"),
            sliderInput("vp_time", "Timepoint (months)", 3, 24, 12, 3),
            actionButton("run_vp", "Run Simulation", class="btn-primary btn-lg")
          ),
          box(title = "Hb Response Distribution", width = 8, status = "info",
            plotlyOutput("vp_hgb_plot", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Spleen Volume % Reduction", width = 6, status = "primary",
            plotlyOutput("vp_sv_plot", height = "280px")
          ),
          box(title = "Response Rate Summary", width = 6, status = "success",
            tableOutput("vp_resp_summary")
          )
        )
      )
    )
  )
)

## ──────────────────────────────────────────────────────────
## SERVER
## ──────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## ── Reactive: run single simulation ─────────────────────
  sim_data <- reactive({
    sim_h <- input$sim_yrs * 365 * 24
    bw    <- input$bw
    p_ert <- if (input$ert_type == "vela") list(VELA_MOD = 1.05) else list(VELA_MOD = 1.0)
    f_elis <- if (input$cyp2d6 == "pm") 0.35 else 0.20
    cl_elis <- if (input$cyp2d6 == "pm") 8.0 else 38.0

    extra_p <- c(list(
      BW      = bw,
      GBA0    = input$gba0,
      SV0     = input$sv0,
      HGB0    = input$hgb0,
      PLT0    = input$plt0,
      BMD0    = input$bmd0,
      USE_ERT  = as.integer(input$use_ert),
      USE_ELIS = as.integer(input$use_elis),
      USE_MIGS = as.integer(input$use_migs),
      DOSE_ERT = input$dose_ert,
      F_ELIS   = f_elis,
      CL_ELIS  = cl_elis
    ), p_ert)

    m <- param(mod_global, extra_p)
    m <- init(m, list(
      GBA = input$gba0, SV = input$sv0,
      HGB = input$hgb0, PLT = input$plt0, BMD = input$bmd0
    ))

    events <- NULL
    if (input$use_ert) {
      events <- ev(cmt = "ERT_C", amt = input$dose_ert * bw / (0.18 * bw),
                   ii = 336, addl = floor(sim_h / 336))
    }
    if (input$use_elis) {
      e2 <- ev(cmt = "ELIS_GUT", amt = 84, ii = 12, addl = floor(sim_h / 12))
      events <- if (is.null(events)) e2 else c(events, e2)
    }
    if (input$use_migs) {
      e3 <- ev(cmt = "MIGS_GUT", amt = input$dose_migs, ii = 8, addl = floor(sim_h / 8))
      events <- if (is.null(events)) e3 else c(events, e3)
    }

    out <- if (is.null(events)) mrgsim(m, end=sim_h, delta=24) else
      mrgsim(m, events=events, end=sim_h, delta=24)
    as.data.frame(out) %>% mutate(time_days = time / 24)
  })

  ## PK plots
  output$pk_ert_plot <- renderPlotly({
    d <- sim_data() %>% filter(time_days <= 28)
    ggplotly(ggplot(d, aes(time_days, ERT_C)) + geom_line(color="steelblue",linewidth=1) +
      labs(x="Days",y="ERT Central Conc (U/mL)") + theme_bw())
  })
  output$pk_srt_plot <- renderPlotly({
    d <- sim_data() %>% filter(time_days <= 14) %>%
      pivot_longer(cols=c(ELIS_C,MIGS_C), names_to="Drug", values_to="Conc")
    ggplotly(ggplot(d, aes(time_days,Conc,color=Drug)) + geom_line(linewidth=1) +
      labs(x="Days",y="Plasma Conc (μg/mL)") + theme_bw())
  })
  output$pk_table <- renderDT({
    data.frame(
      Parameter = c("ERT CL","ERT V1","ERT V2","Eliglustat CL","Eliglustat F","Eliglustat V"),
      Value     = c("1.4 L/h/kg","0.18 L/kg","0.55 L/kg","38 L/h (EM)","20% (EM)","106 L"),
      Source    = c("Aerts 2003","Aerts 2003","Aerts 2003","Lukina 2014","Mistry 2015","Mistry 2015")
    )
  }, options=list(pageLength=6))

  ## Enzyme plots
  output$gba_plot <- renderPlotly({
    d <- sim_data()
    ggplotly(ggplot(d, aes(time_days, GBA_pct)) + geom_line(color="darkgreen",linewidth=1) +
      geom_hline(yintercept=100, linetype="dashed") +
      labs(x="Days", y="GBA Activity (% Normal)") + theme_bw())
  })
  output$gc_plot <- renderPlotly({
    d <- sim_data() %>%
      pivot_longer(cols=c(GC_MAC,GC_SP,GC_LV,GC_BM), names_to="Compartment", values_to="GC")
    ggplotly(ggplot(d, aes(time_days, GC, color=Compartment)) + geom_line(linewidth=1) +
      geom_hline(yintercept=100, linetype="dashed", color="gray") +
      labs(x="Days", y="GC Burden (AU)") + theme_bw())
  })
  output$gcs_plot <- renderPlotly({
    d <- sim_data()
    ggplotly(ggplot(d, aes(time_days, SRT_GCS_inh)) + geom_line(color="purple",linewidth=1) +
      labs(x="Days", y="GCS Inhibition (%)") + theme_bw())
  })
  output$nfkb_plot <- renderPlotly({
    d <- sim_data() %>% pivot_longer(cols=c(NFKB,IL6), names_to="Marker", values_to="Value")
    ggplotly(ggplot(d, aes(time_days, Value, color=Marker)) + geom_line(linewidth=1) +
      labs(x="Days", y="Activity / Concentration") + theme_bw())
  })

  ## Organ plots
  output$sv_plot <- renderPlotly({
    d <- sim_data()
    ggplotly(ggplot(d, aes(time_days, SV)) + geom_line(color="steelblue",linewidth=1) +
      geom_hline(yintercept=1500, linetype="dashed", color="orange") +
      labs(x="Days", y="Spleen Volume (mL)") + theme_bw())
  })
  output$lv_plot <- renderPlotly({
    d <- sim_data()
    ggplotly(ggplot(d, aes(time_days, LV)) + geom_line(color="darkblue",linewidth=1) +
      geom_hline(yintercept=1.25, linetype="dashed", color="orange") +
      labs(x="Days", y="Liver Volume (×normal)") + theme_bw())
  })
  output$hgb_plot <- renderPlotly({
    d <- sim_data()
    ggplotly(ggplot(d, aes(time_days, HGB)) + geom_line(color="firebrick",linewidth=1) +
      geom_hline(yintercept=12, linetype="dashed") +
      labs(x="Days", y="Hemoglobin (g/dL)") + theme_bw())
  })
  output$plt_plot <- renderPlotly({
    d <- sim_data()
    ggplotly(ggplot(d, aes(time_days, PLT)) + geom_line(color="darkorange",linewidth=1) +
      geom_hline(yintercept=100, linetype="dashed") +
      labs(x="Days", y="Platelets (×10⁹/L)") + theme_bw())
  })

  ## Bone
  output$bmd_plot <- renderPlotly({
    d <- sim_data()
    ggplotly(ggplot(d, aes(time_days, BMD)) + geom_line(color="brown",linewidth=1) +
      geom_hline(yintercept=-2.5, linetype="dashed", color="red") +
      labs(x="Days", y="Lumbar T-score") + theme_bw())
  })
  output$oboc_plot <- renderPlotly({
    d <- sim_data() %>% pivot_longer(cols=c(OB,OC), names_to="CellType", values_to="Activity")
    ggplotly(ggplot(d, aes(time_days, Activity, color=CellType)) + geom_line(linewidth=1) +
      labs(x="Days", y="Activity (AU)") + theme_bw())
  })

  ## Scenario comparison (use fixed scenarios)
  all_scenarios_data <- reactive({
    run_sc <- function(use_ert, use_elis, vela_mod=1.0, cl_elis=38, f_elis=0.20, label="") {
      m <- param(mod_global, list(USE_ERT=use_ert, USE_ELIS=use_elis,
                                  VELA_MOD=vela_mod, CL_ELIS=cl_elis, F_ELIS=f_elis))
      evts <- NULL
      if (use_ert == 1) evts <- ev(cmt="ERT_C", amt=60*70/(0.18*70), ii=336, addl=103)
      if (use_elis == 1) {
        e2 <- ev(cmt="ELIS_GUT", amt=84, ii=12, addl=1460)
        evts <- if (is.null(evts)) e2 else c(evts, e2)
      }
      out <- if (is.null(evts)) mrgsim(m, end=17520, delta=24) else
        mrgsim(m, events=evts, end=17520, delta=24)
      as.data.frame(out) %>% mutate(time_days = time / 24, scenario = label)
    }
    bind_rows(
      run_sc(0, 0, label = "S1: Natural History"),
      run_sc(1, 0, label = "S2: Imiglucerase 60U/kg"),
      run_sc(1, 0, vela_mod=1.05, label = "S3: Velaglucerase α"),
      run_sc(0, 1, f_elis=0.20, cl_elis=38, label = "S4: Eliglustat (EM)"),
      run_sc(0, 1, f_elis=0.35, cl_elis=8, label = "S5: Eliglustat (PM)"),
      run_sc(1, 1, label = "S6: Low ERT + Eliglustat")
    )
  })

  output$scenario_plot <- renderPlot({
    d <- all_scenarios_data()
    vars <- input$sc_vars
    if (length(vars) == 0) return(NULL)
    dl <- d %>% pivot_longer(cols = all_of(vars), names_to = "Variable", values_to = "Value")
    ggplot(dl, aes(time_days, Value, color=scenario)) +
      geom_line(linewidth=0.9) +
      facet_wrap(~Variable, scales="free_y", ncol=2) +
      labs(x="Days", y="Value", color="Scenario") +
      theme_bw(base_size=12) +
      theme(legend.position="bottom")
  })

  output$resp_table <- renderDT({
    d <- all_scenarios_data()
    bl <- d %>% filter(time_days==0) %>%
      select(scenario, GL1, SV, HGB, PLT) %>%
      rename_with(~paste0(.x,"_BL"), -scenario)
    d %>% filter(time_days==365) %>%
      left_join(bl, by="scenario") %>%
      transmute(
        Scenario = scenario,
        "GL-1 Δ%" = round((GL1-GL1_BL)/GL1_BL*100,1),
        "SV Δ%" = round((SV-SV_BL)/SV_BL*100,1),
        "Hb Δ g/dL" = round(HGB-HGB_BL,2),
        "PLT Δ%" = round((PLT-PLT_BL)/PLT_BL*100,1)
      )
  }, options=list(pageLength=8))

  ## Biomarkers
  output$gl1_plot <- renderPlotly({
    d <- sim_data()
    ggplotly(ggplot(d, aes(time_days,GL1)) + geom_line(color="dodgerblue",linewidth=1) +
      geom_hline(yintercept=1, linetype="dashed") + labs(x="Days",y="GL-1 (μg/L)") + theme_bw())
  })
  output$lyso_plot <- renderPlotly({
    d <- sim_data()
    ggplotly(ggplot(d, aes(time_days,LYSOGL1)) + geom_line(color="tomato",linewidth=1) +
      geom_hline(yintercept=1.5, linetype="dashed") + labs(x="Days",y="Lyso-GL1 (ng/mL)") + theme_bw())
  })
  output$chitr_plot <- renderPlotly({
    d <- sim_data()
    ggplotly(ggplot(d, aes(time_days,CHITR)) + geom_line(color="seagreen",linewidth=1) +
      geom_hline(yintercept=100, linetype="dashed") + labs(x="Days",y="Chitotriosidase") + theme_bw())
  })
  output$ferr_plot <- renderPlotly({
    d <- sim_data() %>% pivot_longer(cols=c(FERRIT,IL6), names_to="Marker", values_to="Val")
    ggplotly(ggplot(d, aes(time_days,Val,color=Marker)) + geom_line(linewidth=1) +
      labs(x="Days",y="Concentration") + theme_bw())
  })

  ## Disease severity tables
  output$gssi_table <- renderTable({
    sv_x <- input$sv0 / 300
    data.frame(
      Component=c("Spleen (×normal)","Hemoglobin","Platelets (×10⁹/L)","BMD T-score"),
      Value=c(round(sv_x,1), input$hgb0, input$plt0, input$bmd0),
      Severity=c(
        ifelse(sv_x>15,"Severe",ifelse(sv_x>8,"Moderate","Mild")),
        ifelse(input$hgb0<8,"Severe",ifelse(input$hgb0<10,"Moderate","Mild")),
        ifelse(input$plt0<50,"Severe",ifelse(input$plt0<100,"Moderate","Mild")),
        ifelse(input$bmd0< -2.5,"Osteoporosis",ifelse(input$bmd0< -1,"Osteopenia","Normal"))
      )
    )
  })
  output$ssi_table <- renderTable({
    data.frame(Component=c("Bone","Liver","Spleen","Lung","Neurological"),
               Score=c(2,3,3,0,as.integer(input$type)-1))
  })

  ## Virtual population
  vp_results <- eventReactive(input$run_vp, {
    n <- input$n_vp
    set.seed(42)
    bws   <- rnorm(n, 70, 12)
    gba0s <- pmax(rnorm(n, 5, 1.5), 0.5)
    sv0s  <- pmax(rnorm(n, 2200, 400), 400)
    hgb0s <- pmax(rnorm(n, 10.5, 1.2), 7)
    plt0s <- pmax(rnorm(n, 65, 18), 20)
    felis <- pmin(pmax(rnorm(n, 0.20, 0.06), 0.05), 0.50)

    use_ert_vp  <- if (input$vp_scenario == "ert")  1 else 0
    use_elis_vp <- if (input$vp_scenario == "srt")  1 else 0
    sim_h <- input$vp_time * 30 * 24

    withProgress(message="Simulating VPs…", value=0, {
      res <- lapply(seq_len(n), function(i) {
        incProgress(1/n)
        m <- param(mod_global, list(BW=bws[i], GBA0=gba0s[i],
          USE_ERT=use_ert_vp, USE_ELIS=use_elis_vp, F_ELIS=felis[i]))
        m <- init(m, list(GBA=gba0s[i], SV=sv0s[i], HGB=hgb0s[i], PLT=plt0s[i]))
        evts <- NULL
        if (use_ert_vp==1) evts <- ev(cmt="ERT_C", amt=60*bws[i]/(0.18*bws[i]),
                                       ii=336, addl=floor(sim_h/336))
        if (use_elis_vp==1) evts <- ev(cmt="ELIS_GUT",amt=84,ii=12,addl=floor(sim_h/12))
        out <- if (is.null(evts)) mrgsim(m, end=sim_h, delta=24) else
          mrgsim(m, events=evts, end=sim_h, delta=24)
        df <- as.data.frame(out)
        list(
          HGB_bl  = hgb0s[i],
          HGB_end = df$HGB[nrow(df)],
          SV_bl   = sv0s[i],
          SV_end  = df$SV[nrow(df)]
        )
      })
    })
    data.frame(
      HGB_delta = sapply(res, function(r) r$HGB_end - r$HGB_bl),
      SV_pct    = sapply(res, function(r) (r$SV_end - r$SV_bl) / r$SV_bl * 100)
    )
  })

  output$vp_hgb_plot <- renderPlotly({
    d <- vp_results()
    ggplotly(ggplot(d, aes(HGB_delta)) + geom_histogram(bins=30, fill="steelblue", color="white") +
      geom_vline(xintercept=1, linetype="dashed", color="red") +
      labs(x="ΔHb (g/dL)", y="Count", title="Hb Response Distribution") + theme_bw())
  })
  output$vp_sv_plot <- renderPlotly({
    d <- vp_results()
    ggplotly(ggplot(d, aes(SV_pct)) + geom_histogram(bins=30, fill="tomato", color="white") +
      geom_vline(xintercept=-30, linetype="dashed", color="darkblue") +
      labs(x="Spleen Volume % Change", y="Count") + theme_bw())
  })
  output$vp_resp_summary <- renderTable({
    d <- vp_results()
    data.frame(
      Endpoint=c("Hb ≥+1 g/dL","Hb ≥+2 g/dL","SV ≤-30%","SV ≤-50%"),
      Rate=c(
        sprintf("%.1f%%", mean(d$HGB_delta >= 1)*100),
        sprintf("%.1f%%", mean(d$HGB_delta >= 2)*100),
        sprintf("%.1f%%", mean(d$SV_pct <= -30)*100),
        sprintf("%.1f%%", mean(d$SV_pct <= -50)*100)
      )
    )
  })
}

shinyApp(ui, server)
