################################################################################
##  Lupus Nephritis QSP — Interactive Shiny Dashboard
##  Tabs: ① Patient Profile  ② Drug PK  ③ Immune Biomarkers
##        ④ Renal Function   ⑤ Clinical Endpoints  ⑥ Scenario Comparison
################################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

##──────────────────────────────────────────────────────────────────────────────
##  MODEL CODE (embedded)
##──────────────────────────────────────────────────────────────────────────────

model_code <- '
$PARAM
F_MMF=0.94, ka_MMF=1.2, CL_MPA=22.0, V_MPA=55.0,
ka_EHC=0.08, fr_MPAG=0.40, dose_MMF=0.0,
F_HCQ=0.74, ka_HCQ=0.38, CL_HCQ=0.095,
V_blood=7.0, V_tissue=800.0, k12_HCQ=0.03, k21_HCQ=0.0003,
dose_HCQ=0.0, WT=65.0,
F_VCS=0.27, ka_VCS=0.7, CL_VCS=31.6, V_VCS=167.0, dose_VCS=0.0,
CL_BEL=0.182, V1_BEL=4.1, V2_BEL=1.8, Q_BEL=0.12, dose_BEL=0.0,
CL_ANI=0.192, V1_ANI=5.6, V2_ANI=2.4, Q_ANI=0.15, dose_ANI=0.0,
EC50_MPA_B=2.5, EC50_MPA_T=4.0, Emax_MPA=0.85,
EC50_HCQ=0.50, Emax_HCQ=0.70,
EC50_VCS=0.015, Emax_VCS=0.80,
Kd_BEL=0.001, Emax_BEL=0.75,
EC50_ANI=0.005, Emax_ANI=0.90,
k_Bnaive_in=0.05, k_Bnaive_deg=0.014,
k_GC_stim=0.008, k_GC_deg=0.05,
k_PC_diff=0.10, k_PC_deg=0.003,
k_Tfh_stim=0.006, k_Tfh_deg=0.08,
k_Treg_base=0.003, k_Treg_deg=0.03,
BAFF_base=1.0,
k_Ab_prod=0.15, k_Ab_deg=0.004, Ab_baseline=100.0,
C3_baseline=0.80, k_C3_synth=0.01, k_C3_consume=0.04, k_C3_deg=0.005,
C4_baseline=0.25, k_C4_synth=0.003, k_C4_consume=0.06, k_C4_deg=0.006,
k_pod_inj=0.03, k_pod_repair=0.005, Pod_max=1.0,
k_prot_rise=0.15, k_prot_fall=0.04, UPCR_base=2.5,
eGFR_base=65.0, k_eGFR_loss=0.0008, k_eGFR_rec=0.001, eGFR_min=15.0,
IFN_base=1.5, k_IFN_BAFF=0.4

$CMT
MPA_gut MPA_plasma MPAG_gut HCQ_blood HCQ_tissue VCS_plasma
BEL_central BEL_periph ANI_central ANI_periph
B_naive B_GC Plasma_cell Tfh Treg
Anti_dsDNA C3 C4
Podocyte_inj Proteinuria eGFR_cmt

$INIT
MPA_gut=0, MPA_plasma=0, MPAG_gut=0,
HCQ_blood=0, HCQ_tissue=0, VCS_plasma=0,
BEL_central=0, BEL_periph=0, ANI_central=0, ANI_periph=0,
B_naive=250, B_GC=20, Plasma_cell=120, Tfh=15, Treg=30,
Anti_dsDNA=280, C3=0.55, C4=0.12,
Podocyte_inj=0.35, Proteinuria=2.5, eGFR_cmt=65

$ODE
double C_MPA   = MPA_plasma/V_MPA;
double C_HCQ   = HCQ_blood/(V_blood*WT);
double C_VCS   = VCS_plasma/V_VCS;
double C_BEL   = BEL_central/V1_BEL;
double C_ANI   = ANI_central/V1_ANI;
double E_MPA_B = Emax_MPA*pow(C_MPA,1.5)/(pow(EC50_MPA_B,1.5)+pow(C_MPA,1.5));
double E_MPA_T = Emax_MPA*C_MPA/(EC50_MPA_T+C_MPA);
double E_HCQ   = Emax_HCQ*C_HCQ/(EC50_HCQ+C_HCQ);
double E_VCS   = Emax_VCS*C_VCS/(EC50_VCS+C_VCS);
double E_BEL   = Emax_BEL*C_BEL/(Kd_BEL+C_BEL);
double E_ANI   = Emax_ANI*C_ANI/(EC50_ANI+C_ANI);
double IFN_cur = IFN_base*(1.0-E_ANI)*(1.0-0.5*E_HCQ);
double BAFF    = BAFF_base*(1.0+k_IFN_BAFF*(IFN_cur-1.0))*(1.0-E_BEL);
if(BAFF<0) BAFF=0.001;
double IC_load = (Anti_dsDNA/Ab_baseline)*(C3_baseline/(C3+0.01));
dxdt_MPA_gut    = dose_MMF*F_MMF - ka_MMF*MPA_gut;
dxdt_MPA_plasma = ka_MMF*MPA_gut + ka_EHC*MPAG_gut
                  - (CL_MPA/V_MPA)*MPA_plasma
                  - (fr_MPAG*CL_MPA/V_MPA)*MPA_plasma;
dxdt_MPAG_gut   = (fr_MPAG*CL_MPA/V_MPA)*MPA_plasma - ka_EHC*MPAG_gut;
dxdt_HCQ_blood  = dose_HCQ*F_HCQ - (ka_HCQ+CL_HCQ*WT/(V_blood*WT))*HCQ_blood
                  + k21_HCQ*HCQ_tissue;
dxdt_HCQ_tissue = k12_HCQ*HCQ_blood - k21_HCQ*HCQ_tissue;
dxdt_VCS_plasma = dose_VCS*F_VCS*ka_VCS - (CL_VCS/V_VCS)*VCS_plasma;
dxdt_BEL_central= dose_BEL - (CL_BEL+Q_BEL)/V1_BEL*BEL_central
                  + Q_BEL/V2_BEL*BEL_periph;
dxdt_BEL_periph = Q_BEL/V1_BEL*BEL_central - Q_BEL/V2_BEL*BEL_periph;
dxdt_ANI_central= dose_ANI - (CL_ANI+Q_ANI)/V1_ANI*ANI_central
                  + Q_ANI/V2_ANI*ANI_periph;
dxdt_ANI_periph = Q_ANI/V1_ANI*ANI_central - Q_ANI/V2_ANI*ANI_periph;
dxdt_B_naive    = k_Bnaive_in*BAFF - k_Bnaive_deg*B_naive;
dxdt_B_GC       = k_GC_stim*Tfh*B_naive*(1-E_MPA_B) - k_GC_deg*B_GC;
dxdt_Plasma_cell= k_PC_diff*B_GC - k_PC_deg*Plasma_cell;
dxdt_Tfh        = k_Tfh_stim*IFN_cur*B_naive*(1-E_MPA_T)*(1-E_VCS) - k_Tfh_deg*Tfh;
dxdt_Treg       = k_Treg_base*(1+0.3*E_VCS) - k_Treg_deg*Treg;
dxdt_Anti_dsDNA = k_Ab_prod*Plasma_cell - k_Ab_deg*Anti_dsDNA;
dxdt_C3         = k_C3_synth - k_C3_consume*IC_load*C3 - k_C3_deg*C3;
dxdt_C4         = k_C4_synth - k_C4_consume*IC_load*C4 - k_C4_deg*C4;
double Pod_stim = k_pod_inj*IC_load*(C3_baseline/(C3+0.001));
double Pod_rep  = k_pod_repair*(1+E_VCS);
dxdt_Podocyte_inj = Pod_stim*(1-Podocyte_inj) - Pod_rep*Podocyte_inj;
dxdt_Proteinuria  = k_prot_rise*Podocyte_inj - k_prot_fall*Proteinuria*(1-0.5*Podocyte_inj);
dxdt_eGFR_cmt     = k_eGFR_rec*(eGFR_base-eGFR_cmt) - k_eGFR_loss*Podocyte_inj*eGFR_cmt;

$TABLE
capture C_MPA=MPA_plasma/V_MPA;
capture C_HCQ_b=HCQ_blood/(V_blood*WT);
capture C_VCS=VCS_plasma/V_VCS;
capture C_BEL=BEL_central/V1_BEL;
capture C_ANI=ANI_central/V1_ANI;
capture BAFF_level=BAFF_base*(1.0+k_IFN_BAFF*(IFN_base*(1.0-E_ANI)*(1.0-0.5*E_HCQ)-1.0))*(1.0-E_BEL);
capture IC_load=(Anti_dsDNA/Ab_baseline)*(C3_baseline/(C3+0.01));
capture UPCR=Proteinuria;
capture eGFR=eGFR_cmt;
capture CRRP=(Proteinuria<0.5 && eGFR_cmt>=60)?1.0:0.0;
capture PRRP=(Proteinuria<1.0 && eGFR_cmt>=60)?1.0:0.0;
double E_MPA_B_out=Emax_MPA*pow(C_MPA,1.5)/(pow(EC50_MPA_B,1.5)+pow(C_MPA,1.5));
double E_VCS_out=Emax_VCS*C_VCS/(EC50_VCS+C_VCS);
double E_BEL_out=Emax_BEL*C_BEL/(Kd_BEL+C_BEL);
double E_ANI_out=Emax_ANI*C_ANI/(EC50_ANI+C_ANI);
'

##──────────────────────────────────────────────────────────────────────────────
##  HELPER: run simulation
##──────────────────────────────────────────────────────────────────────────────
run_sim <- function(mod, mmf_g, hcq_mg, vcs_bid_mg,
                    bel_mg_4w, ani_mg_4w,
                    wt_kg, ifn_level, sim_days, init_upcr, init_egfr) {
  # Convert doses
  dose_mmf <- mmf_g  * 1000 / 24
  dose_hcq <- hcq_mg / 24
  dose_vcs <- vcs_bid_mg * 2 / 24
  dose_bel <- bel_mg_4w / (4 * 7 * 24)
  dose_ani <- ani_mg_4w / (4 * 7 * 24)

  p <- param(mod,
    dose_MMF = dose_mmf, dose_HCQ = dose_hcq,
    dose_VCS = dose_vcs, dose_BEL = dose_bel,
    dose_ANI = dose_ani, WT = wt_kg, IFN_base = ifn_level
  )
  p <- init(p, Proteinuria = init_upcr, eGFR_cmt = init_egfr,
              Anti_dsDNA = ifelse(init_upcr > 2, 280, 150),
              C3 = ifelse(init_upcr > 2, 0.55, 0.75),
              C4 = ifelse(init_upcr > 2, 0.12, 0.18))

  out <- mrgsim(p, end = sim_days, delta = 0.5) %>% as.data.frame()
  out$week <- out$time / 7
  out
}

##──────────────────────────────────────────────────────────────────────────────
##  UI
##──────────────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(
    title = "Lupus Nephritis QSP",
    titleWidth = 260
  ),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "patient",   icon = icon("user")),
      menuItem("Drug PK",            tabName = "pk",        icon = icon("pills")),
      menuItem("Immune Biomarkers",  tabName = "immune",    icon = icon("dna")),
      menuItem("Renal Function",     tabName = "renal",     icon = icon("kidneys")),
      menuItem("Clinical Endpoints", tabName = "endpoint",  icon = icon("chart-line")),
      menuItem("Scenario Comparison",tabName = "scenario",  icon = icon("sliders"))
    ),

    tags$hr(),
    tags$div(style = "padding:10px 15px; font-size:12px; color:#aaa;",
      "Lupus Nephritis QSP", tags$br(),
      "mrgsolve v0.12 · 20 compartments", tags$br(),
      "AURORA 2021 · BLISS-LN 2020"
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box { border-top: 3px solid #3c8dbc; }
      .small-box .icon-large { font-size: 65px !important; }
    "))),

    tabItems(

      ##── Tab 1: Patient Profile ───────────────────────────────────────────────
      tabItem("patient",
        fluidRow(
          box(title = "Patient Parameters", status = "primary", solidHeader = TRUE,
              width = 4,
              numericInput("wt",     "Body Weight (kg)", 65, 40, 120, 5),
              selectInput("ln_class","LN Class (ISN/RPS)",
                          c("Class III" = "III", "Class IV" = "IV",
                            "Class V (pure)" = "V", "Class III/V or IV/V" = "mixed")),
              numericInput("init_upcr", "Baseline UPCR (g/g)", 2.5, 0.1, 15, 0.5),
              numericInput("init_egfr", "Baseline eGFR (mL/min)", 65, 10, 130, 5),
              sliderInput("ifn_level", "IFN Signature Level",
                          min = 0.5, max = 5.0, value = 1.5, step = 0.5),
              numericInput("sim_wks",  "Simulation Duration (weeks)", 52, 12, 104, 4)
          ),
          box(title = "Baseline Summary", status = "info", solidHeader = TRUE, width = 4,
              valueBoxOutput("vb_upcr",   width = 12),
              valueBoxOutput("vb_egfr",   width = 12),
              valueBoxOutput("vb_ifn",    width = 12)
          ),
          box(title = "LN Class Reference", status = "warning", solidHeader = TRUE,
              width = 4,
              tags$table(class = "table table-condensed",
                tags$thead(tags$tr(tags$th("Class"), tags$th("Description"), tags$th("Treatment"))),
                tags$tbody(
                  tags$tr(tags$td("I/II"), tags$td("Minimal/Mesangial"),
                          tags$td("HCQ ± low-dose GC")),
                  tags$tr(tags$td("III"), tags$td("Focal Proliferative"),
                          tags$td("MMF or CYC + GC")),
                  tags$tr(tags$td("IV"), tags$td("Diffuse Proliferative"),
                          tags$td("MMF or CYC + GC ± add-on")),
                  tags$tr(tags$td("V"), tags$td("Membranous"),
                          tags$td("MMF ± Voclosporin")),
                  tags$tr(tags$td("VI"), tags$td("Advanced Sclerotic"),
                          tags$td("RRT consideration"))
                )
              )
          )
        )
      ),

      ##── Tab 2: Drug PK ──────────────────────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(title = "Drug Regimen", status = "primary", solidHeader = TRUE, width = 3,
              sliderInput("dose_mmf",  "MMF Daily Dose (g)", 0, 3, 2,   0.5),
              sliderInput("dose_hcq",  "HCQ Daily Dose (mg)",0, 400, 400, 100),
              sliderInput("dose_vcs",  "Voclosporin BID (mg)", 0, 30, 0, 2.7),
              sliderInput("dose_bel",  "Belimumab q4w (mg)", 0, 1000, 0, 100),
              sliderInput("dose_ani",  "Anifrolumab q4w (mg)", 0, 300, 0, 300),
              actionButton("run_pk", "Run Simulation", class = "btn-primary btn-block",
                           icon = icon("play"))
          ),
          box(title = "MPA Concentration-Time", status = "info",
              solidHeader = TRUE, width = 9,
              plotlyOutput("pk_mpa", height = "220px")
          )
        ),
        fluidRow(
          box(title = "HCQ Blood Concentration", status = "info", solidHeader = TRUE, width = 4,
              plotlyOutput("pk_hcq", height = "220px")),
          box(title = "Voclosporin Plasma", status = "info", solidHeader = TRUE, width = 4,
              plotlyOutput("pk_vcs", height = "220px")),
          box(title = "Biologics (BEL / ANI)", status = "info", solidHeader = TRUE, width = 4,
              plotlyOutput("pk_bio", height = "220px"))
        )
      ),

      ##── Tab 3: Immune Biomarkers ─────────────────────────────────────────────
      tabItem("immune",
        fluidRow(
          box(title = "B Cell Compartments", status = "success",
              solidHeader = TRUE, width = 6,
              plotlyOutput("imm_bcell", height = "300px")),
          box(title = "T Cell Compartments", status = "success",
              solidHeader = TRUE, width = 6,
              plotlyOutput("imm_tcell", height = "300px"))
        ),
        fluidRow(
          box(title = "Anti-dsDNA Antibody", status = "warning",
              solidHeader = TRUE, width = 4,
              plotlyOutput("imm_dsdna", height = "260px")),
          box(title = "Complement C3 / C4", status = "warning",
              solidHeader = TRUE, width = 4,
              plotlyOutput("imm_comp", height = "260px")),
          box(title = "BAFF Level & IC Load", status = "warning",
              solidHeader = TRUE, width = 4,
              plotlyOutput("imm_baff", height = "260px"))
        )
      ),

      ##── Tab 4: Renal Function ────────────────────────────────────────────────
      tabItem("renal",
        fluidRow(
          box(title = "eGFR Trajectory", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("ren_egfr", height = "300px")),
          box(title = "Proteinuria (UPCR)", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("ren_upcr", height = "300px"))
        ),
        fluidRow(
          box(title = "Podocyte Injury Index", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("ren_pod", height = "260px")),
          box(title = "Renal Biomarkers Summary (Week 52)", status = "info",
              solidHeader = TRUE, width = 6,
              DTOutput("ren_table"))
        )
      ),

      ##── Tab 5: Clinical Endpoints ────────────────────────────────────────────
      tabItem("endpoint",
        fluidRow(
          valueBoxOutput("ep_crr",  width = 4),
          valueBoxOutput("ep_prr",  width = 4),
          valueBoxOutput("ep_egfr_w52", width = 4)
        ),
        fluidRow(
          box(title = "Renal Response Over Time", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("ep_response", height = "320px")),
          box(title = "SLEDAI Renal Score (proxy)", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("ep_sledai", height = "320px"))
        ),
        fluidRow(
          box(title = "Renal Response Waterfall (Weekly UPCR)", status = "warning",
              solidHeader = TRUE, width = 12,
              plotlyOutput("ep_waterfall", height = "250px"))
        )
      ),

      ##── Tab 6: Scenario Comparison ───────────────────────────────────────────
      tabItem("scenario",
        fluidRow(
          box(title = "Scenario Settings", status = "primary", solidHeader = TRUE, width = 3,
              h4("Scenario A (Reference)"),
              sliderInput("sa_mmf", "MMF (g/day)", 0, 3, 3, 0.5),
              sliderInput("sa_hcq", "HCQ (mg/day)", 0, 400, 400, 100),
              checkboxInput("sa_vcs", "Add Voclosporin (23.7 mg BID)", FALSE),
              checkboxInput("sa_bel", "Add Belimumab (650 mg q4w)", FALSE),
              checkboxInput("sa_ani", "Add Anifrolumab (300 mg q4w)", FALSE),
              tags$hr(),
              h4("Scenario B (Comparator)"),
              sliderInput("sb_mmf", "MMF (g/day)", 0, 3, 2, 0.5),
              sliderInput("sb_hcq", "HCQ (mg/day)", 0, 400, 400, 100),
              checkboxInput("sb_vcs", "Add Voclosporin", TRUE),
              checkboxInput("sb_bel", "Add Belimumab", FALSE),
              checkboxInput("sb_ani", "Add Anifrolumab", FALSE),
              actionButton("run_scen", "Compare Scenarios",
                           class = "btn-warning btn-block", icon = icon("balance-scale"))
          ),
          box(title = "UPCR Comparison", status = "info",
              solidHeader = TRUE, width = 9,
              plotlyOutput("sc_upcr", height = "300px"))
        ),
        fluidRow(
          box(title = "eGFR Comparison", status = "success", solidHeader = TRUE, width = 6,
              plotlyOutput("sc_egfr", height = "270px")),
          box(title = "Anti-dsDNA Comparison", status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("sc_ab", height = "270px"))
        ),
        fluidRow(
          box(title = "Response Rates at Week 52", status = "danger",
              solidHeader = TRUE, width = 12,
              DTOutput("sc_table"))
        )
      )
    )
  )
)

##──────────────────────────────────────────────────────────────────────────────
##  SERVER
##──────────────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # Compile model once
  mod <- mcode("LN_shiny", model_code, quiet = TRUE)

  # Reactive simulation (triggered by pk run button or on input changes)
  sim_data <- eventReactive(
    list(input$run_pk, input$wt, input$init_upcr, input$init_egfr,
         input$ifn_level, input$sim_wks, input$dose_mmf, input$dose_hcq,
         input$dose_vcs, input$dose_bel, input$dose_ani),
    ignoreNULL = FALSE,
    {
      run_sim(mod,
              mmf_g       = input$dose_mmf,
              hcq_mg      = input$dose_hcq,
              vcs_bid_mg  = input$dose_vcs,
              bel_mg_4w   = input$dose_bel,
              ani_mg_4w   = input$dose_ani,
              wt_kg       = input$wt,
              ifn_level   = input$ifn_level,
              sim_days    = input$sim_wks * 7,
              init_upcr   = input$init_upcr,
              init_egfr   = input$init_egfr)
    }
  )

  ## Value boxes
  output$vb_upcr <- renderValueBox({
    valueBox(input$init_upcr, "Baseline UPCR (g/g)", icon = icon("tint"),
             color = if (input$init_upcr > 2) "red" else "orange")
  })
  output$vb_egfr <- renderValueBox({
    valueBox(input$init_egfr, "Baseline eGFR", icon = icon("heartbeat"),
             color = if (input$init_egfr < 60) "yellow" else "green")
  })
  output$vb_ifn <- renderValueBox({
    lbl <- if (input$ifn_level >= 3) "High IFN Sig" else if (input$ifn_level >= 1.5) "Med IFN Sig" else "Low"
    valueBox(input$ifn_level, paste("Type I IFN Level —", lbl), icon = icon("shield"),
             color = if (input$ifn_level >= 3) "purple" else "blue")
  })

  ## ── PK Tab ──────────────────────────────────────────────────────────────────
  output$pk_mpa <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~week, y = ~C_MPA, type = "scatter", mode = "lines",
            line = list(color = "#1565c0", width = 2)) %>%
      layout(title = "MPA Free Concentration", xaxis = list(title = "Week"),
             yaxis = list(title = "C_MPA (mg/L)"))
  })
  output$pk_hcq <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~week, y = ~C_HCQ_b, type = "scatter", mode = "lines",
            line = list(color = "#00695c", width = 2)) %>%
      layout(title = "HCQ Blood", xaxis = list(title = "Week"),
             yaxis = list(title = "C_HCQ (mg/L)"))
  })
  output$pk_vcs <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~week, y = ~C_VCS, type = "scatter", mode = "lines",
            line = list(color = "#4a148c", width = 2)) %>%
      layout(title = "Voclosporin Plasma", xaxis = list(title = "Week"),
             yaxis = list(title = "C_VCS (mg/L)"))
  })
  output$pk_bio <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_trace(x = ~week, y = ~C_BEL, name = "Belimumab",
                type = "scatter", mode = "lines",
                line = list(color = "#e65100", width = 2)) %>%
      add_trace(x = ~week, y = ~C_ANI, name = "Anifrolumab",
                type = "scatter", mode = "lines",
                line = list(color = "#880e4f", width = 2)) %>%
      layout(title = "Biologic Concentrations", xaxis = list(title = "Week"),
             yaxis = list(title = "mg/L"))
  })

  ## ── Immune Tab ──────────────────────────────────────────────────────────────
  output$imm_bcell <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_trace(x=~week, y=~B_naive,    name="B Naive",    mode="lines", type="scatter",
                line=list(color="#42a5f5")) %>%
      add_trace(x=~week, y=~B_GC,       name="GC B",       mode="lines", type="scatter",
                line=list(color="#1e88e5")) %>%
      add_trace(x=~week, y=~Plasma_cell,name="Plasma Cell", mode="lines", type="scatter",
                line=list(color="#0d47a1")) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="Cells/µL"))
  })
  output$imm_tcell <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_trace(x=~week, y=~Tfh,  name="Tfh",  mode="lines", type="scatter",
                line=list(color="#43a047")) %>%
      add_trace(x=~week, y=~Treg, name="Treg", mode="lines", type="scatter",
                line=list(color="#ef6c00")) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="Cells/µL"))
  })
  output$imm_dsdna <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~week, y=~Anti_dsDNA, type="scatter", mode="lines",
            line=list(color="#c62828", width=2)) %>%
      add_trace(x=c(0, max(d$week)), y=c(200,200), mode="lines",
                line=list(dash="dash", color="grey"), showlegend=FALSE) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="IU/mL"),
             annotations=list(x=2, y=205, text="Upper normal", showarrow=FALSE, font=list(size=10)))
  })
  output$imm_comp <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_trace(x=~week, y=~C3, name="C3 (g/L)", yaxis="y1",
                mode="lines", type="scatter", line=list(color="#1565c0")) %>%
      add_trace(x=~week, y=~C4, name="C4 (g/L)", yaxis="y2",
                mode="lines", type="scatter", line=list(color="#e65100")) %>%
      layout(yaxis=list(title="C3 (g/L)"),
             yaxis2=list(title="C4 (g/L)", overlaying="y", side="right"),
             xaxis=list(title="Week"))
  })
  output$imm_baff <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_trace(x=~week, y=~BAFF_level, name="BAFF",    yaxis="y1",
                mode="lines", type="scatter", line=list(color="#7b1fa2")) %>%
      add_trace(x=~week, y=~IC_load,    name="IC Load", yaxis="y2",
                mode="lines", type="scatter", line=list(color="#c62828")) %>%
      layout(yaxis=list(title="BAFF (rel.)"),
             yaxis2=list(title="IC Load", overlaying="y", side="right"),
             xaxis=list(title="Week"))
  })

  ## ── Renal Tab ───────────────────────────────────────────────────────────────
  output$ren_egfr <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~week, y=~eGFR, type="scatter", mode="lines",
            line=list(color="#00695c", width=2)) %>%
      add_trace(x=c(0,max(d$week)), y=c(60,60), mode="lines",
                line=list(dash="dash", color="orange"), showlegend=FALSE) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="eGFR (mL/min/1.73m²)"))
  })
  output$ren_upcr <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~week, y=~UPCR, type="scatter", mode="lines",
            line=list(color="#c62828", width=2)) %>%
      add_trace(x=c(0,max(d$week)), y=c(0.5,0.5), mode="lines",
                line=list(dash="dash", color="blue"), name="CRR threshold",
                showlegend=TRUE) %>%
      add_trace(x=c(0,max(d$week)), y=c(1.0,1.0), mode="lines",
                line=list(dash="dot", color="green"), name="PRR threshold",
                showlegend=TRUE) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="UPCR (g/g)"))
  })
  output$ren_pod <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~week, y=~Podocyte_inj, type="scatter", mode="lines",
            line=list(color="#880e4f", width=2), fill="tozeroy",
            fillcolor="rgba(136,14,79,0.1)") %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="Podocyte Injury Index (0–1)"))
  })
  output$ren_table <- renderDT({
    d <- sim_data()
    w52 <- d %>% filter(abs(week - input$sim_wks) <= 0.6) %>% tail(1)
    tbl <- data.frame(
      Endpoint = c("UPCR (g/g)", "eGFR (mL/min)", "Podocyte Injury",
                   "Anti-dsDNA (IU/mL)", "C3 (g/L)", "C4 (g/L)",
                   "CRR achieved", "PRR achieved"),
      Value = c(round(w52$UPCR, 2), round(w52$eGFR, 1), round(w52$Podocyte_inj, 3),
                round(w52$Anti_dsDNA, 0), round(w52$C3, 3), round(w52$C4, 3),
                ifelse(w52$CRRP == 1, "Yes", "No"),
                ifelse(w52$PRRP == 1, "Yes", "No"))
    )
    datatable(tbl, options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  ## ── Endpoints Tab ───────────────────────────────────────────────────────────
  output$ep_crr <- renderValueBox({
    d <- sim_data()
    wk52 <- d %>% filter(abs(week - input$sim_wks) <= 0.6) %>% tail(1)
    col <- if (wk52$CRRP == 1) "green" else "red"
    valueBox(ifelse(wk52$CRRP == 1, "ACHIEVED", "NOT MET"),
             "Complete Renal Response (UPCR<0.5 + eGFR≥60)",
             icon = icon("check"), color = col)
  })
  output$ep_prr <- renderValueBox({
    d <- sim_data()
    wk52 <- d %>% filter(abs(week - input$sim_wks) <= 0.6) %>% tail(1)
    col <- if (wk52$PRRP == 1) "green" else "orange"
    valueBox(ifelse(wk52$PRRP == 1, "ACHIEVED", "NOT MET"),
             "Partial Renal Response (UPCR<1.0 + eGFR≥60)",
             icon = icon("check-circle"), color = col)
  })
  output$ep_egfr_w52 <- renderValueBox({
    d <- sim_data()
    wk52 <- d %>% filter(abs(week - input$sim_wks) <= 0.6) %>% tail(1)
    valueBox(round(wk52$eGFR, 1), "eGFR at Week 52",
             icon = icon("heartbeat"), color = "blue")
  })
  output$ep_response <- renderPlotly({
    d <- sim_data()
    wk_filtered <- d %>% filter(week >= 0)
    plot_ly(wk_filtered) %>%
      add_trace(x=~week, y=~UPCR, name="UPCR", mode="lines", type="scatter",
                line=list(color="#c62828", width=2)) %>%
      add_trace(x=c(0,max(d$week)), y=c(0.5,0.5), mode="lines",
                line=list(dash="dash", color="#1565c0"), name="CRR", showlegend=TRUE) %>%
      add_trace(x=c(0,max(d$week)), y=c(1.0,1.0), mode="lines",
                line=list(dash="dot", color="#43a047"), name="PRR", showlegend=TRUE) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="UPCR (g/g)"))
  })
  output$ep_sledai <- renderPlotly({
    d <- sim_data()
    # Simplified SLEDAI renal proxy: hematuria + proteinuria + eGFR decline
    d$SLEDAI_proxy <- pmax(0,
      2 * (d$UPCR > 0.5) +
      4 * (d$Podocyte_inj > 0.5) +
      2 * (d$eGFR < 60) +
      2 * (d$Anti_dsDNA > 200)
    )
    plot_ly(d, x=~week, y=~SLEDAI_proxy, type="scatter", mode="lines",
            line=list(color="#37474f", width=2), fill="tozeroy",
            fillcolor="rgba(55,71,79,0.15)") %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="SLEDAI Renal Proxy Score"))
  })
  output$ep_waterfall <- renderPlotly({
    d <- sim_data() %>%
      filter(week %% 4 < 0.4) %>%
      mutate(UPCR_change = (UPCR - input$init_upcr) / input$init_upcr * 100,
             wk_label = paste0("Wk", round(week)))
    colors <- ifelse(d$UPCR_change < 0, "#1565c0", "#c62828")
    plot_ly(d, x=~wk_label, y=~UPCR_change, type="bar",
            marker=list(color=colors)) %>%
      add_trace(x=c(d$wk_label[1], tail(d$wk_label,1)), y=c(-50,-50),
                mode="lines", type="scatter",
                line=list(dash="dash", color="orange"), name="≥50% reduction") %>%
      layout(xaxis=list(title="Visit", tickangle=-45),
             yaxis=list(title="UPCR % Change from Baseline"))
  })

  ## ── Scenario Comparison Tab ──────────────────────────────────────────────────
  sc_data <- eventReactive(input$run_scen, {
    simA <- run_sim(mod,
                    mmf_g=input$sa_mmf, hcq_mg=input$sa_hcq,
                    vcs_bid_mg=ifelse(input$sa_vcs,23.7,0),
                    bel_mg_4w=ifelse(input$sa_bel,650,0),
                    ani_mg_4w=ifelse(input$sa_ani,300,0),
                    wt_kg=input$wt, ifn_level=input$ifn_level,
                    sim_days=input$sim_wks*7,
                    init_upcr=input$init_upcr, init_egfr=input$init_egfr)
    simA$Scenario <- "Scenario A"

    simB <- run_sim(mod,
                    mmf_g=input$sb_mmf, hcq_mg=input$sb_hcq,
                    vcs_bid_mg=ifelse(input$sb_vcs,23.7,0),
                    bel_mg_4w=ifelse(input$sb_bel,650,0),
                    ani_mg_4w=ifelse(input$sb_ani,300,0),
                    wt_kg=input$wt, ifn_level=input$ifn_level,
                    sim_days=input$sim_wks*7,
                    init_upcr=input$init_upcr, init_egfr=input$init_egfr)
    simB$Scenario <- "Scenario B"
    bind_rows(simA, simB)
  })

  output$sc_upcr <- renderPlotly({
    d <- sc_data()
    plot_ly(d, x=~week, y=~UPCR, color=~Scenario,
            type="scatter", mode="lines",
            colors=c("#1565c0","#c62828")) %>%
      add_trace(x=c(0,max(d$week)), y=c(0.5,0.5), mode="lines", inherit=FALSE,
                line=list(dash="dash", color="grey"), showlegend=FALSE) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="UPCR (g/g)"))
  })
  output$sc_egfr <- renderPlotly({
    d <- sc_data()
    plot_ly(d, x=~week, y=~eGFR, color=~Scenario,
            type="scatter", mode="lines",
            colors=c("#1565c0","#c62828")) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="eGFR (mL/min)"))
  })
  output$sc_ab <- renderPlotly({
    d <- sc_data()
    plot_ly(d, x=~week, y=~Anti_dsDNA, color=~Scenario,
            type="scatter", mode="lines",
            colors=c("#1565c0","#c62828")) %>%
      layout(xaxis=list(title="Week"), yaxis=list(title="Anti-dsDNA (IU/mL)"))
  })
  output$sc_table <- renderDT({
    d <- sc_data()
    tbl <- d %>%
      group_by(Scenario) %>%
      filter(abs(week - input$sim_wks) <= 0.6) %>%
      summarise(
        `UPCR (g/g)` = round(last(UPCR), 2),
        `eGFR (mL/min)` = round(last(eGFR), 1),
        `Anti-dsDNA (IU/mL)` = round(last(Anti_dsDNA), 0),
        `C3 (g/L)` = round(last(C3), 3),
        `CRR` = ifelse(last(CRRP) == 1, "Yes", "No"),
        `PRR` = ifelse(last(PRRP) == 1, "Yes", "No"),
        .groups = "drop"
      )
    datatable(tbl, options = list(dom = "t"), rownames = FALSE)
  })
}

shinyApp(ui, server)
