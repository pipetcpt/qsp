## =============================================================================
## IPF QSP Shiny Dashboard
## Disease: Idiopathic Pulmonary Fibrosis (IPF)
## Drugs  : Pirfenidone · Nintedanib · Combination
## Author : Claude Code Routine (CCR) · 2026-06-17
## Tabs   : 1.Patient Profile · 2.Drug PK · 3.Lung Function (FVC/DLCO)
##          4.Biomarkers · 5.Scenario Comparison · 6.Mechanistic Map
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(tidyverse)
library(ggplot2)
library(plotly)
library(DT)

## ── Inline mrgsolve model (compact version for Shiny) ─────────────────────

ipf_code <- '
$PARAM
ka_P=1.74 F_P=0.81 CL_P=8.4 V1_P=20.4 V2_P=14.0 Q_P=3.6
ka_N=0.80 F_N=0.047 CL_N=22.0 V1_N=730.0 V2_N=900.0 Q_N=15.0
AEC2_ss=1.0 k_AEC=0.005 k_rep=0.003
kprod_TGFb=0.08 kdeg_TGFb=0.04 kact_TGFb=0.12
EC50_P_TGFb=30.0 Emax_P_TGFb=0.65
kact_F=0.015 kdiff_M=0.025 kapop_F=0.008 F_ss=1.0
EC50_N_F=15.0 Emax_N_F=0.70
kprod_Col=0.010 kdeg_Col=0.002 kprod_MMP=0.020 kdeg_MMP=0.030
kprod_TIMP=0.018 kdeg_TIMP=0.020 Col_ss=1.0
kM2_act=0.012 kdeg_M2=0.008 EC50_P_M2=25.0
kprod_ROS=0.04 kdeg_ROS=0.05 kfb_ROS=0.02 Emax_P_ROS=0.45
FVC_base=80.0 k_FVC_loss=0.0026 DLCO_base=65.0 k_DLCO_loss=0.0018
MW_P=185.22 MW_N=539.63

$CMT DEPOT_P CENT_P PERI_P DEPOT_N CENT_N PERI_N
     AEC2 TGFb M2 ROS FIBRO MYOFIB COLLAGEN MMP TIMP FVC_st DLCO_st

$INIT
AEC2=1.0 TGFb=1.0 M2=1.0 ROS=1.0 FIBRO=1.0 MYOFIB=1.0
COLLAGEN=1.0 MMP=1.0 TIMP=1.0 FVC_st=80.0 DLCO_st=65.0

$MAIN
double Cp_P = CENT_P / V1_P;
double Cn_nM = (CENT_N / V1_N) / MW_N * 1e3;
double inh_P_TGFb = Emax_P_TGFb * Cp_P / (EC50_P_TGFb + Cp_P);
double inh_P_M2   = Emax_P_TGFb * Cp_P / (EC50_P_M2   + Cp_P);
double inh_P_ROS  = Emax_P_ROS  * Cp_P / (EC50_P_TGFb + Cp_P);
double inh_N_F    = Emax_N_F    * Cn_nM / (EC50_N_F    + Cn_nM);

$ODE
dxdt_DEPOT_P = -ka_P * DEPOT_P;
dxdt_CENT_P  = ka_P*F_P*DEPOT_P - (CL_P+Q_P)/V1_P*CENT_P + Q_P/V2_P*PERI_P;
dxdt_PERI_P  = Q_P/V1_P*CENT_P - Q_P/V2_P*PERI_P;
dxdt_DEPOT_N = -ka_N * DEPOT_N;
dxdt_CENT_N  = ka_N*F_N*DEPOT_N - (CL_N+Q_N)/V1_N*CENT_N + Q_N/V2_N*PERI_N;
dxdt_PERI_N  = Q_N/V1_N*CENT_N - Q_N/V2_N*PERI_N;

double AEC2_damage_rate = k_AEC * ROS * AEC2;
double AEC2_repair_rate = k_rep * AEC2_ss * (1.0 - AEC2);
dxdt_AEC2 = AEC2_repair_rate - AEC2_damage_rate;

double TGFb_prod = kprod_TGFb + kact_TGFb*(AEC2_ss-AEC2) + 0.06*M2 + 0.04*MYOFIB;
dxdt_TGFb = TGFb_prod*(1.0-inh_P_TGFb) - kdeg_TGFb*TGFb;

dxdt_M2 = kM2_act*TGFb*(1.0-inh_P_M2) - kdeg_M2*M2;

double ROS_prod = kprod_ROS + kfb_ROS*(AEC2_ss-AEC2)*TGFb;
dxdt_ROS = ROS_prod*(1.0-inh_P_ROS) - kdeg_ROS*ROS;

dxdt_FIBRO  = kact_F*TGFb*F_ss*(1.0-inh_N_F) - kapop_F*FIBRO;
dxdt_MYOFIB = kdiff_M*FIBRO*TGFb*(1.0-inh_N_F*0.5) - 0.006*MYOFIB;

double Col_prod = kprod_Col * MYOFIB;
double Col_deg  = kdeg_Col * MMP / (TIMP + 0.1) * COLLAGEN;
dxdt_COLLAGEN = Col_prod - Col_deg;

dxdt_MMP  = (kprod_MMP*FIBRO*0.5 + kprod_MMP*M2*0.5) - kdeg_MMP*MMP;
dxdt_TIMP = kprod_TIMP*TGFb*MYOFIB - kdeg_TIMP*TIMP;

double col_excess = (COLLAGEN > 1.0) ? (COLLAGEN - 1.0) : 0.0;
dxdt_FVC_st  = -k_FVC_loss * col_excess * FVC_st;
dxdt_DLCO_st = -k_DLCO_loss * ((AEC2_ss-AEC2) + 0.5*col_excess) * DLCO_st;

$TABLE
double Cp_pirf = CENT_P / V1_P;
double Cn_nint = (CENT_N / V1_N) / MW_N * 1e3;
double FVC_pct  = FVC_st;
double DLCO_pct = DLCO_st;
double Col_norm = COLLAGEN;
double AEC2_norm = AEC2;
double TGFb_norm = TGFb;
double MMP_norm  = MMP;
double TIMP_norm = TIMP;
double Myofib_norm = MYOFIB;
double ROS_norm  = ROS;
double MMP_TIMP_ratio = (TIMP > 0) ? MMP / TIMP : 1.0;
double Periostin_proxy = COLLAGEN * 1.2 * MYOFIB;
double MMP7_proxy = MMP * TGFb * 0.8;
double KL6_proxy = (AEC2_ss - AEC2 + 0.5) * 800.0;

$CAPTURE Cp_pirf Cn_nint FVC_pct DLCO_pct Col_norm AEC2_norm TGFb_norm
         MMP_norm TIMP_norm Myofib_norm ROS_norm MMP_TIMP_ratio
         Periostin_proxy MMP7_proxy KL6_proxy
'

mod <- mcode("IPF_QSP_shiny", ipf_code, quiet=TRUE)

## ── Helper: run simulation ─────────────────────────────────────────────────

run_sim <- function(dur_weeks, fvc_base, dlco_base,
                    use_pirf, dose_pirf_mg,
                    use_nint, dose_nint_mg,
                    disease_severity) {

  dur_h   <- dur_weeks * 7 * 24
  obs_t   <- seq(0, dur_h, by=24)

  # Severity multipliers on key params
  sev_mult <- c(mild=0.7, moderate=1.0, severe=1.5)[disease_severity]
  params <- list(FVC_st = fvc_base, DLCO_st = dlco_base,
                 k_AEC = 0.005 * sev_mult,
                 kprod_TGFb = 0.08 * sev_mult,
                 kprod_ROS  = 0.04 * sev_mult)

  ev_list <- list()
  if (use_pirf && dose_pirf_mg > 0) {
    ev_list[["pirf"]] <- ev(cmt=1, amt=dose_pirf_mg*1e3, ii=8, addl=ceiling(dur_h/8))
  }
  if (use_nint && dose_nint_mg > 0) {
    ev_list[["nint"]] <- ev(cmt=4, amt=dose_nint_mg*1e6, ii=12, addl=ceiling(dur_h/12))
  }

  ev_obj <- if (length(ev_list)==0) ev(time=0, amt=0) else
    if (length(ev_list)==1) ev_list[[1]] else c(ev_list[[1]], ev_list[[2]])

  out <- mod %>%
    param(params) %>%
    init(FVC_st=fvc_base, DLCO_st=dlco_base) %>%
    mrgsim(events=ev_obj, end=dur_h, delta=24, digits=4) %>%
    as.data.frame()
  out$Week <- out$time / (7*24)
  out
}

## ─── UI ───────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = span(icon("lungs"), "IPF QSP Dashboard"),
    titleWidth = 320
  ),

  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("Patient Profile",     tabName="patient",   icon=icon("user-md")),
      menuItem("Drug PK",             tabName="pk",        icon=icon("pills")),
      menuItem("Lung Function",       tabName="lungfx",    icon=icon("lungs")),
      menuItem("Biomarkers",          tabName="biomarker", icon=icon("vials")),
      menuItem("Scenario Comparison", tabName="scenario",  icon=icon("chart-bar")),
      menuItem("Mechanistic Map",     tabName="mech",      icon=icon("project-diagram"))
    ),
    hr(),
    h4("  Patient Parameters", style="color:#ECF0F1; padding-left:10px"),
    sliderInput("fvc_base",   "Baseline FVC (% predicted):",
                min=40, max=95, value=75, step=1),
    sliderInput("dlco_base",  "Baseline DLCO (% predicted):",
                min=25, max=90, value=60, step=1),
    selectInput("severity", "Disease Severity:",
                choices=c("Mild"="mild","Moderate"="moderate","Severe"="severe"),
                selected="moderate"),
    sliderInput("dur_weeks", "Simulation Duration (weeks):",
                min=13, max=104, value=52, step=13),
    hr(),
    h4("  Treatment", style="color:#ECF0F1; padding-left:10px"),
    checkboxInput("use_pirf", "Pirfenidone", value=TRUE),
    conditionalPanel("input.use_pirf",
      selectInput("dose_pirf", "Pirfenidone Dose (mg TID):",
                  choices=c("267 mg"=267,"534 mg"=534,"801 mg"=801),
                  selected=801)
    ),
    checkboxInput("use_nint", "Nintedanib", value=FALSE),
    conditionalPanel("input.use_nint",
      selectInput("dose_nint", "Nintedanib Dose (mg BID):",
                  choices=c("100 mg"=100,"150 mg"=150),
                  selected=150)
    ),
    actionButton("run_sim", "Run Simulation",
                 class="btn-primary btn-block",
                 icon=icon("play"), style="margin:10px")
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background-color: #f5f7fa; }
        .box { border-radius: 8px; }
        .info-box { border-radius: 8px; }
        .nav-tabs-custom>.tab-content { padding:10px; }
      "))
    ),

    tabItems(

      ## ── TAB 1: Patient Profile ──────────────────────────────────────────
      tabItem(tabName="patient",
        fluidRow(
          valueBoxOutput("box_fvc",    width=3),
          valueBoxOutput("box_dlco",   width=3),
          valueBoxOutput("box_gap",    width=3),
          valueBoxOutput("box_stage",  width=3)
        ),
        fluidRow(
          box(title="Patient Summary", width=6, solidHeader=TRUE, status="primary",
            tableOutput("pt_summary_table")
          ),
          box(title="GAP Index Calculator", width=6, solidHeader=TRUE, status="warning",
            numericInput("pt_age",  "Age (years):", value=68, min=40, max=90),
            selectInput("pt_sex",  "Sex:", choices=c("Female"=0,"Male"=1), selected=1),
            sliderInput("pt_fvc",  "FVC (% predicted):", min=30, max=100, value=75),
            sliderInput("pt_dlco", "DLCO (% predicted):", min=20, max=100, value=60),
            actionButton("calc_gap", "Calculate GAP Index", class="btn-warning"),
            br(), br(),
            verbatimTextOutput("gap_result")
          )
        ),
        fluidRow(
          box(title="Disease Natural History — FVC Trajectory",
              width=12, solidHeader=TRUE, status="info",
            plotlyOutput("nat_history_plot", height="350px"),
            p("Gray band = typical IPF progression range (−100 to −300 mL/year ≈ −1.5 to −4% predicted/year)")
          )
        )
      ),

      ## ── TAB 2: Drug PK ─────────────────────────────────────────────────
      tabItem(tabName="pk",
        fluidRow(
          box(title="Pirfenidone PK Profile (72h, 801 mg TID)",
              width=6, solidHeader=TRUE, status="primary",
            plotlyOutput("pk_pirf_plot", height="320px"),
            p(em("F=81%, t½=2.4h, ka=1.74/h, CL=8.4 L/h; CYP1A2 major metabolism"))
          ),
          box(title="Nintedanib PK Profile (72h, 150 mg BID)",
              width=6, solidHeader=TRUE, status="success",
            plotlyOutput("pk_nint_plot", height="320px"),
            p(em("F=4.7%, t½=10h, ka=0.8/h, CL=22 L/h; P-gp efflux limits absorption"))
          )
        ),
        fluidRow(
          box(title="PK Parameters Reference", width=6, solidHeader=TRUE, status="info",
            tableOutput("pk_param_table")
          ),
          box(title="PK Variability (BSV) — Monte Carlo Preview",
              width=6, solidHeader=TRUE, status="warning",
            plotlyOutput("pk_variability_plot", height="280px"),
            p("50 simulated individuals with log-normal IIV (CV~30% on CL/V)")
          )
        )
      ),

      ## ── TAB 3: Lung Function ───────────────────────────────────────────
      tabItem(tabName="lungfx",
        fluidRow(
          valueBoxOutput("box_fvc52",  width=3),
          valueBoxOutput("box_dlco52", width=3),
          valueBoxOutput("box_fvc_decl", width=3),
          valueBoxOutput("box_resp",   width=3)
        ),
        fluidRow(
          box(title="FVC % Predicted — Time Course",
              width=6, solidHeader=TRUE, status="danger",
            plotlyOutput("fvc_plot", height="360px")
          ),
          box(title="DLCO % Predicted — Time Course",
              width=6, solidHeader=TRUE, status="warning",
            plotlyOutput("dlco_plot", height="360px")
          )
        ),
        fluidRow(
          box(title="Lung Function Decline Rate",
              width=6, solidHeader=TRUE, status="info",
            plotlyOutput("fvc_rate_plot", height="280px")
          ),
          box(title="Spirometry Reference Zones",
              width=6, solidHeader=TRUE, status="primary",
            tableOutput("spirometry_ref_table")
          )
        )
      ),

      ## ── TAB 4: Biomarkers ──────────────────────────────────────────────
      tabItem(tabName="biomarker",
        fluidRow(
          box(title="TGF-β1 & Fibroblast Activation",
              width=6, solidHeader=TRUE, status="danger",
            plotlyOutput("tgfb_fibro_plot", height="320px")
          ),
          box(title="ECM — Collagen & MMP:TIMP Balance",
              width=6, solidHeader=TRUE, status="warning",
            plotlyOutput("ecm_plot", height="320px")
          )
        ),
        fluidRow(
          box(title="Serum Biomarkers (Proxy)",
              width=6, solidHeader=TRUE, status="success",
            plotlyOutput("bm_serum_plot", height="320px"),
            p(em("MMP-7, KL-6, Periostin — scaled to approximate clinical ranges"))
          ),
          box(title="Oxidative Stress & AEC Damage",
              width=6, solidHeader=TRUE, status="primary",
            plotlyOutput("ros_aec_plot", height="320px")
          )
        ),
        fluidRow(
          box(title="Biomarker Correlation Heatmap (at Week 52)",
              width=12, solidHeader=TRUE, status="info",
            plotlyOutput("bm_heatmap", height="280px")
          )
        )
      ),

      ## ── TAB 5: Scenario Comparison ─────────────────────────────────────
      tabItem(tabName="scenario",
        fluidRow(
          box(title="Multi-Scenario FVC Comparison",
              width=12, solidHeader=TRUE, status="primary",
            plotlyOutput("scenario_fvc_plot", height="400px")
          )
        ),
        fluidRow(
          box(title="Scenario Results Table (Week 52)",
              width=7, solidHeader=TRUE, status="info",
            DTOutput("scenario_table")
          ),
          box(title="Treatment Effect (%FVC decline reduction vs. placebo)",
              width=5, solidHeader=TRUE, status="success",
            plotlyOutput("te_plot", height="300px")
          )
        ),
        fluidRow(
          box(title="Dose-Response: Pirfenidone",
              width=6, solidHeader=TRUE, status="warning",
            plotlyOutput("dr_pirf_plot", height="280px")
          ),
          box(title="Dose-Response: Nintedanib",
              width=6, solidHeader=TRUE, status="success",
            plotlyOutput("dr_nint_plot", height="280px")
          )
        )
      ),

      ## ── TAB 6: Mechanistic Map ─────────────────────────────────────────
      tabItem(tabName="mech",
        fluidRow(
          box(title="IPF QSP Mechanistic Map", width=12,
            solidHeader=TRUE, status="primary",
            imageOutput("mech_map", height="700px"),
            p(strong("Figure:"), "Full IPF QSP mechanistic map showing 130+ components across 12 biological subsystems.",
              "Rendered from Graphviz DOT source using fdp layout engine.",
              "Key: diamond nodes = biomarkers, hexagon = clinical endpoints,",
              "cylinder = dosing inputs. Blue/teal = drug PK compartments.",
              style="font-size:12px; color:#555")
          )
        ),
        fluidRow(
          box(title="Pathway Interaction Summary", width=6, solidHeader=TRUE, status="info",
            tableOutput("pathway_table")
          ),
          box(title="Drug Mechanism Summary", width=6, solidHeader=TRUE, status="warning",
            tableOutput("drug_mech_table")
          )
        )
      )

    ) # end tabItems
  )
)

## ─── SERVER ───────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  ## ── Reactive simulation ─────────────────────────────────────────────────
  sim_data <- eventReactive(input$run_sim, {
    run_sim(
      dur_weeks       = input$dur_weeks,
      fvc_base        = input$fvc_base,
      dlco_base       = input$dlco_base,
      use_pirf        = input$use_pirf,
      dose_pirf_mg    = as.numeric(input$dose_pirf),
      use_nint        = input$use_nint,
      dose_nint_mg    = as.numeric(input$dose_nint),
      disease_severity = input$severity
    )
  }, ignoreNULL=FALSE)

  ## All-scenario comparison (fixed scenarios)
  all_scenarios <- reactive({
    dur_h <- input$dur_weeks * 7 * 24
    sev_mult <- c(mild=0.7, moderate=1.0, severe=1.5)[input$severity]
    params <- list(k_AEC=0.005*sev_mult, kprod_TGFb=0.08*sev_mult,
                   FVC_st=input$fvc_base, DLCO_st=input$dlco_base)

    sc_list <- list(
      "Placebo"      = ev(time=0, amt=0),
      "Pirf 801 TID" = ev(cmt=1, amt=801e3, ii=8, addl=ceiling(dur_h/8)),
      "Nint 150 BID" = ev(cmt=4, amt=150e6, ii=12, addl=ceiling(dur_h/12)),
      "Combination"  = c(ev(cmt=1, amt=801e3, ii=8, addl=ceiling(dur_h/8)),
                         ev(cmt=4, amt=150e6, ii=12, addl=ceiling(dur_h/12))),
      "Pirf 534 TID" = ev(cmt=1, amt=534e3, ii=8, addl=ceiling(dur_h/8))
    )

    lapply(names(sc_list), function(nm) {
      out <- mod %>%
        param(params) %>%
        init(FVC_st=input$fvc_base, DLCO_st=input$dlco_base) %>%
        mrgsim(events=sc_list[[nm]], end=dur_h, delta=24, digits=4) %>%
        as.data.frame()
      out$Scenario <- nm
      out$Week     <- out$time / (7*24)
      out
    }) %>% bind_rows()
  })

  ## PK simulation (72h single dose)
  pk_sim <- reactive({
    mod %>%
      mrgsim(
        events = c(ev(cmt=1, amt=801e3, time=0),
                   ev(cmt=4, amt=150e6, time=0)),
        end=72, delta=0.25
      ) %>% as.data.frame()
  })

  ## ── TAB 1: Patient Profile ─────────────────────────────────────────────

  output$box_fvc  <- renderValueBox({
    valueBox(paste0(input$fvc_base, "%"), "Baseline FVC", icon=icon("lungs"), color="blue")
  })
  output$box_dlco <- renderValueBox({
    valueBox(paste0(input$dlco_base, "%"), "Baseline DLCO", icon=icon("wind"), color="yellow")
  })
  output$box_gap  <- renderValueBox({
    gap <- calc_gap_score(as.numeric(input$pt_age), as.numeric(input$pt_sex),
                          input$fvc_base, input$dlco_base)
    valueBox(gap$score, "GAP Score", icon=icon("chart-line"),
             color=ifelse(gap$score<=3,"green",ifelse(gap$score<=5,"yellow","red")))
  })
  output$box_stage <- renderValueBox({
    gap <- calc_gap_score(as.numeric(input$pt_age), as.numeric(input$pt_sex),
                          input$fvc_base, input$dlco_base)
    valueBox(gap$stage, "GAP Stage", icon=icon("flag"), color=gap$color)
  })

  output$pt_summary_table <- renderTable({
    data.frame(
      Parameter = c("Disease Severity","Baseline FVC","Baseline DLCO",
                    "Expected FVC Decline/yr","Expected 1-yr Mortality","Therapy"),
      Value = c(
        tools::toTitleCase(input$severity),
        paste0(input$fvc_base, "% predicted"),
        paste0(input$dlco_base, "% predicted"),
        ifelse(input$severity=="mild","~100 mL/yr",
               ifelse(input$severity=="moderate","~200 mL/yr","~300 mL/yr")),
        ifelse(input$fvc_base>70,"5-10%",ifelse(input$fvc_base>50,"15-30%","30-50%")),
        paste(c(if(input$use_pirf) paste0("Pirfenidone ",input$dose_pirf," mg TID"),
                if(input$use_nint) paste0("Nintedanib ",input$dose_nint," mg BID"),
                if(!input$use_pirf && !input$use_nint) "None (Placebo)"),
              collapse=" + ")
      )
    )
  })

  output$gap_result <- renderText({
    input$calc_gap
    isolate({
      gap <- calc_gap_score(input$pt_age, as.numeric(input$pt_sex),
                            input$pt_fvc, input$pt_dlco)
      paste0(
        "GAP Score: ", gap$score, "\n",
        "Stage: ", gap$stage, "\n",
        "Mortality (1-yr): ", gap$mort1yr, "\n",
        "Mortality (2-yr): ", gap$mort2yr, "\n",
        "Mortality (3-yr): ", gap$mort3yr
      )
    })
  })

  output$nat_history_plot <- renderPlotly({
    weeks <- seq(0, input$dur_weeks, by=1)
    fvc_mid  <- input$fvc_base - (2.5 / 52) * weeks
    fvc_low  <- input$fvc_base - (4.0 / 52) * weeks
    fvc_high <- input$fvc_base - (1.0 / 52) * weeks

    df <- data.frame(Week=weeks, Mid=fvc_mid, Low=fvc_low, High=fvc_high)

    plot_ly(df) %>%
      add_ribbons(x=~Week, ymin=~Low, ymax=~High, fillcolor="rgba(200,200,200,0.3)",
                  line=list(color="transparent"), name="Typical Range") %>%
      add_lines(x=~Week, y=~Mid, line=list(color="black", dash="dash", width=2),
                name="Median Natural History") %>%
      layout(title="Natural History — FVC Decline Trajectory",
             xaxis=list(title="Week"), yaxis=list(title="FVC % predicted"),
             shapes=list(list(type="line", x0=0, x1=max(weeks),
                              y0=50, y1=50, line=list(color="red", dash="dot"))))
  })

  ## ── TAB 2: Drug PK ─────────────────────────────────────────────────────

  output$pk_pirf_plot <- renderPlotly({
    pk <- pk_sim()
    plot_ly(pk, x=~time, y=~Cp_pirf, type="scatter", mode="lines",
            line=list(color="#2471A3", width=2.5), name="Pirfenidone") %>%
      add_segments(x=0, xend=72, y=30, yend=30,
                   line=list(color="red", dash="dash"), name="EC50") %>%
      layout(title="Pirfenidone: Plasma Cp (µg/mL)",
             xaxis=list(title="Time (h)"),
             yaxis=list(title="Concentration (µg/mL)"))
  })

  output$pk_nint_plot <- renderPlotly({
    pk <- pk_sim()
    plot_ly(pk, x=~time, y=~Cn_nint, type="scatter", mode="lines",
            line=list(color="#27AE60", width=2.5), name="Nintedanib") %>%
      add_segments(x=0, xend=72, y=20, yend=20,
                   line=list(color="red", dash="dash"), name="IC50") %>%
      layout(title="Nintedanib: Plasma Cn (nM)",
             xaxis=list(title="Time (h)"),
             yaxis=list(title="Concentration (nM)"))
  })

  output$pk_param_table <- renderTable({
    data.frame(
      Parameter = c("Pirfenidone F","Pirfenidone t½","Pirfenidone ka",
                    "Pirfenidone CL","Pirfenidone V1",
                    "Nintedanib F","Nintedanib t½","Nintedanib ka",
                    "Nintedanib CL","Nintedanib V1"),
      Value = c("81%","2.4 h","1.74 h⁻¹","8.4 L/h","20.4 L",
                "4.7%","10 h","0.80 h⁻¹","22.0 L/h","730 L"),
      Source = c(rep("Rubino 2009 ClinPK",5), rep("Stopfer 2011 ClinPK",5))
    )
  })

  output$pk_variability_plot <- renderPlotly({
    set.seed(42)
    n <- 50; t_seq <- seq(0, 24, by=0.5)
    pk_sims <- lapply(1:n, function(i) {
      cl_i <- 8.4 * exp(rnorm(1, 0, 0.3))
      v1_i <- 20.4 * exp(rnorm(1, 0, 0.25))
      cp <- (801e3 * 0.81 * 1.74) / (v1_i * (1.74 - cl_i/v1_i)) *
        (exp(-cl_i/v1_i * t_seq) - exp(-1.74 * t_seq))
      data.frame(time=t_seq, Cp=pmax(cp, 0), ID=i)
    }) %>% bind_rows()

    pk_sims %>%
      group_by(time) %>%
      summarise(p5=quantile(Cp,0.05), p50=median(Cp), p95=quantile(Cp,0.95)) ->
      pk_ci

    plot_ly(pk_ci) %>%
      add_ribbons(x=~time, ymin=~p5, ymax=~p95,
                  fillcolor="rgba(36,113,163,0.2)",
                  line=list(color="transparent"), name="5th-95th %ile") %>%
      add_lines(x=~time, y=~p50, line=list(color="#2471A3", width=2), name="Median") %>%
      layout(title="Pirfenidone Population PK Variability (n=50)",
             xaxis=list(title="Time (h)"), yaxis=list(title="Cp (µg/mL)"))
  })

  ## ── TAB 3: Lung Function ───────────────────────────────────────────────

  output$box_fvc52  <- renderValueBox({
    d <- sim_data() %>% filter(Week == max(Week))
    fvc_final <- mean(d$FVC_pct)
    valueBox(paste0(round(fvc_final,1),"%"),
             paste0("FVC at ", input$dur_weeks," wks"),
             icon=icon("lungs"),
             color=ifelse(fvc_final>70,"green",ifelse(fvc_final>50,"yellow","red")))
  })
  output$box_dlco52 <- renderValueBox({
    d <- sim_data() %>% filter(Week == max(Week))
    valueBox(paste0(round(mean(d$DLCO_pct),1),"%"),
             paste0("DLCO at ", input$dur_weeks," wks"),
             icon=icon("wind"), color="yellow")
  })
  output$box_fvc_decl <- renderValueBox({
    d0 <- sim_data() %>% filter(Week==0)
    dn <- sim_data() %>% filter(Week==max(Week))
    decl <- mean(d0$FVC_pct) - mean(dn$FVC_pct)
    valueBox(paste0("-",round(decl,2),"%"), "FVC Decline (total)",
             icon=icon("arrow-down"), color=ifelse(decl<3,"green",ifelse(decl<6,"yellow","red")))
  })
  output$box_resp <- renderValueBox({
    d0 <- sim_data() %>% filter(Week==0)
    dn <- sim_data() %>% filter(Week==max(Week))
    decl_treat <- mean(d0$FVC_pct) - mean(dn$FVC_pct)
    nat_decl    <- 2.5 * input$dur_weeks / 52
    resp_pct    <- max(0, (1 - decl_treat / nat_decl) * 100)
    valueBox(paste0(round(resp_pct,0),"%"), "Treatment Response",
             icon=icon("check-circle"),
             color=ifelse(resp_pct>40,"green",ifelse(resp_pct>15,"yellow","red")))
  })

  output$fvc_plot <- renderPlotly({
    d <- sim_data()
    p <- plot_ly(d, x=~Week, y=~FVC_pct, type="scatter", mode="lines",
                 line=list(color="#2471A3", width=2.5), name="Treatment") %>%
      add_segments(x=0, xend=max(d$Week), y=50, yend=50,
                   line=list(color="red", dash="dash"), name="Critical: 50%") %>%
      add_segments(x=0, xend=max(d$Week), y=70, yend=70,
                   line=list(color="orange", dash="dot"), name="Moderate: 70%") %>%
      layout(title="FVC % Predicted",
             xaxis=list(title="Week"), yaxis=list(title="FVC (%)"))
    p
  })

  output$dlco_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~Week, y=~DLCO_pct, type="scatter", mode="lines",
            line=list(color="#E67E22", width=2.5), name="DLCO") %>%
      add_segments(x=0, xend=max(d$Week), y=35, yend=35,
                   line=list(color="red", dash="dash"), name="Severe: 35%") %>%
      layout(title="DLCO % Predicted",
             xaxis=list(title="Week"), yaxis=list(title="DLCO (%)"))
  })

  output$fvc_rate_plot <- renderPlotly({
    d <- sim_data() %>% arrange(Week)
    rate <- diff(d$FVC_pct) / diff(d$Week) * 52  # annualized
    df_rate <- data.frame(Week=d$Week[-1], Rate=rate)
    plot_ly(df_rate, x=~Week, y=~Rate, type="scatter", mode="lines",
            line=list(color="#8E44AD", width=2), name="FVC change/yr") %>%
      layout(title="Annualized FVC Decline Rate",
             xaxis=list(title="Week"), yaxis=list(title="FVC %/year"))
  })

  output$spirometry_ref_table <- renderTable({
    data.frame(
      Category      = c("Normal","Mild IPF","Moderate IPF","Severe IPF","Very Severe"),
      `FVC %`       = c(">80%","66-80%","51-65%","<50%","<40%"),
      `DLCO %`      = c(">70%","56-70%","41-55%","<40%","<30%"),
      `1yr Mortality`= c("<5%","5-10%","15-30%","30-50%",">50%")
    )
  })

  ## ── TAB 4: Biomarkers ──────────────────────────────────────────────────

  output$tgfb_fibro_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_lines(x=~Week, y=~TGFb_norm, name="TGF-β1",
                line=list(color="#C0392B", width=2)) %>%
      add_lines(x=~Week, y=~Myofib_norm, name="Myofibroblasts",
                line=list(color="#6C3483", width=2)) %>%
      layout(title="TGF-β1 & Myofibroblast Activation",
             xaxis=list(title="Week"), yaxis=list(title="Normalized level"))
  })

  output$ecm_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_lines(x=~Week, y=~Col_norm, name="Collagen",
                line=list(color="#784212", width=2)) %>%
      add_lines(x=~Week, y=~MMP_norm, name="MMP Activity",
                line=list(color="#1A5276", width=2)) %>%
      add_lines(x=~Week, y=~TIMP_norm, name="TIMP Activity",
                line=list(color="#117A65", width=2)) %>%
      add_lines(x=~Week, y=~MMP_TIMP_ratio, name="MMP:TIMP Ratio",
                line=list(color="#E67E22", width=2, dash="dash")) %>%
      layout(title="ECM Remodeling: Collagen & MMP/TIMP",
             xaxis=list(title="Week"), yaxis=list(title="Normalized level"))
  })

  output$bm_serum_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_lines(x=~Week, y=~MMP7_proxy, name="Serum MMP-7 (proxy)",
                line=list(color="#27AE60", width=2)) %>%
      add_lines(x=~Week, y=~Periostin_proxy, name="Periostin (proxy)",
                line=list(color="#F39C12", width=2)) %>%
      add_lines(x=~Week, y=~KL6_proxy, name="KL-6 (proxy, U/mL)",
                line=list(color="#2980B9", width=2)) %>%
      layout(title="Serum Biomarker Trajectories",
             xaxis=list(title="Week"), yaxis=list(title="Level (proxy)"))
  })

  output$ros_aec_plot <- renderPlotly({
    d <- sim_data()
    plot_ly(d) %>%
      add_lines(x=~Week, y=~ROS_norm, name="ROS",
                line=list(color="#E74C3C", width=2)) %>%
      add_lines(x=~Week, y=~AEC2_norm, name="AEC-II Population",
                line=list(color="#2ECC71", width=2)) %>%
      layout(title="Oxidative Stress & Epithelial Integrity",
             xaxis=list(title="Week"), yaxis=list(title="Normalized level"))
  })

  output$bm_heatmap <- renderPlotly({
    d <- sim_data() %>% filter(Week==max(Week)) %>%
      select(TGFb_norm, Myofib_norm, Col_norm, MMP_norm, TIMP_norm,
             ROS_norm, AEC2_norm, MMP7_proxy, KL6_proxy)
    cormat <- cor(d)
    plot_ly(x=colnames(cormat), y=rownames(cormat), z=cormat,
            type="heatmap", colorscale="RdBu", zmin=-1, zmax=1) %>%
      layout(title="Biomarker Correlation (Week 52)")
  })

  ## ── TAB 5: Scenarios ───────────────────────────────────────────────────

  output$scenario_fvc_plot <- renderPlotly({
    d <- all_scenarios()
    sc_colors <- c("Placebo"="#E74C3C", "Pirf 801 TID"="#2471A3",
                   "Nint 150 BID"="#27AE60", "Combination"="#8E44AD",
                   "Pirf 534 TID"="#E67E22")
    p <- plot_ly()
    for (sc in unique(d$Scenario)) {
      dd <- d %>% filter(Scenario==sc)
      p  <- add_lines(p, data=dd, x=~Week, y=~FVC_pct, name=sc,
                      line=list(color=sc_colors[sc], width=2.2))
    }
    p %>% layout(title="Scenario Comparison: FVC % Predicted Over Time",
                 xaxis=list(title="Week"), yaxis=list(title="FVC (%)"),
                 legend=list(orientation="h"))
  })

  output$scenario_table <- renderDT({
    d <- all_scenarios()
    d %>% filter(Week %in% c(0, 13, 26, 39, 52)) %>%
      group_by(Scenario, Week) %>%
      summarise(FVC=round(mean(FVC_pct),2), DLCO=round(mean(DLCO_pct),2),
                TGFb=round(mean(TGFb_norm),3), Col=round(mean(Col_norm),3),
                .groups="drop") %>%
      datatable(options=list(pageLength=10, scrollX=TRUE),
                rownames=FALSE, class="compact")
  })

  output$te_plot <- renderPlotly({
    d    <- all_scenarios()
    base <- d %>% filter(Scenario=="Placebo", Week==max(Week))
    rest <- d %>% filter(Scenario!="Placebo", Week==max(Week))
    fvc0 <- input$fvc_base
    fvc_pl <- mean(base$FVC_pct)
    nat_decl <- fvc0 - fvc_pl

    te <- rest %>%
      group_by(Scenario) %>%
      summarise(FVC_final=mean(FVC_pct), .groups="drop") %>%
      mutate(Decline=fvc0-FVC_final,
             Reduction=pmax(0, (1 - Decline/nat_decl)*100))

    plot_ly(te, x=~Scenario, y=~Reduction, type="bar",
            marker=list(color=c("#2471A3","#27AE60","#8E44AD","#E67E22")),
            name="% Reduction") %>%
      layout(title="FVC Decline Reduction vs Placebo",
             xaxis=list(title=""), yaxis=list(title="Reduction (%)"),
             shapes=list(list(type="line", x0=-0.5, x1=3.5, y0=50, y1=50,
                              line=list(color="red", dash="dash"))))
  })

  output$dr_pirf_plot <- renderPlotly({
    dur_h <- input$dur_weeks * 7 * 24
    doses <- c(267, 534, 801, 1068)
    fvc_f <- sapply(doses, function(d) {
      out <- mod %>%
        init(FVC_st=input$fvc_base, DLCO_st=input$dlco_base) %>%
        mrgsim(events=ev(cmt=1, amt=d*1e3, ii=8, addl=ceiling(dur_h/8)),
               end=dur_h, delta=24, digits=4) %>%
        as.data.frame() %>% filter(time==max(time))
      mean(out$FVC_pct)
    })
    plot_ly(x=doses, y=fvc_f, type="scatter", mode="lines+markers",
            line=list(color="#2471A3", width=2), marker=list(size=8, color="#2471A3")) %>%
      layout(title="Pirfenidone Dose-Response (FVC at 52 wks)",
             xaxis=list(title="Dose (mg TID)"), yaxis=list(title="FVC % predicted"))
  })

  output$dr_nint_plot <- renderPlotly({
    dur_h <- input$dur_weeks * 7 * 24
    doses <- c(50, 100, 150, 200)
    fvc_f <- sapply(doses, function(d) {
      out <- mod %>%
        init(FVC_st=input$fvc_base, DLCO_st=input$dlco_base) %>%
        mrgsim(events=ev(cmt=4, amt=d*1e6, ii=12, addl=ceiling(dur_h/12)),
               end=dur_h, delta=24, digits=4) %>%
        as.data.frame() %>% filter(time==max(time))
      mean(out$FVC_pct)
    })
    plot_ly(x=doses, y=fvc_f, type="scatter", mode="lines+markers",
            line=list(color="#27AE60", width=2), marker=list(size=8, color="#27AE60")) %>%
      layout(title="Nintedanib Dose-Response (FVC at 52 wks)",
             xaxis=list(title="Dose (mg BID)"), yaxis=list(title="FVC % predicted"))
  })

  ## ── TAB 6: Mechanistic Map ─────────────────────────────────────────────

  output$mech_map <- renderImage({
    png_path <- file.path(".", "ipf_qsp_model.png")
    if (!file.exists(png_path)) {
      png_path <- system.file("ipf_qsp_model.png", package="base")
    }
    list(src=png_path, alt="IPF QSP Mechanistic Map",
         style="max-width:100%; height:auto")
  }, deleteFile=FALSE)

  output$pathway_table <- renderTable({
    data.frame(
      Cluster = c("Alveolar Epithelium","Macrophages","TGF-β Signaling",
                  "Fibroblast/Myofibroblast","ECM Remodeling",
                  "Growth Factors","Oxidative Stress","PK Pirfenidone",
                  "PK Nintedanib","Drug Targets","Clinical Endpoints","Comorbidities"),
      Nodes = c(15, 16, 17, 15, 17, 15, 15, 10, 10, 12, 17, 14),
      `Key Mediators` = c(
        "AEC-II, KL-6, SP-D, ER Stress",
        "M1/M2, NLRP3, IL-1β, TNF-α, MMP-12",
        "TGF-β1, SMAD2/3, CTGF, PAI-1",
        "α-SMA, PDGFR, FGFR, YAP/TAZ",
        "Collagen I/III, MMP-7, TIMP-1, Periostin",
        "PDGF, FGF-2, VEGF, ET-1, CXCL12",
        "NOX4, ROS, Nrf2, GSH, mTOR",
        "Gut→Plasma→Periph; CYP1A2 metabolism",
        "Gut→Plasma→Periph; P-gp efflux; biliary",
        "TGF-β/PDGF/FGF/VEGFR/ROS inhibition",
        "FVC, DLCO, MMP-7, KL-6, GAP index",
        "GERD, PH, Smoking, TERT/MUC5B"
      )
    )
  })

  output$drug_mech_table <- renderTable({
    data.frame(
      Drug  = c("Pirfenidone","Pirfenidone","Pirfenidone","Pirfenidone","Pirfenidone",
                "Nintedanib","Nintedanib","Nintedanib","Nintedanib"),
      Target = c("TGF-β1 production","TNF-α","PDGF","FGF","ROS/Antioxidant",
                 "FGFR-1/2/3","VEGFR-1/2/3","PDGFRα/β","Src/Lck/RET"),
      Mechanism = c("↓Transcription/Secretion","↓Cytokine","↓Growth Factor","↓Growth Factor","↑GSH",
                    "Competitive TKI","Competitive TKI","Competitive TKI","Kinase Inhibition"),
      EC50_IC50 = c("~30 µg/mL","~25 µg/mL","~30 µg/mL","~30 µg/mL","~25 µg/mL",
                    "~20 nM","~34 nM","~59 nM","~16 nM")
    )
  })

}

## ─── HELPER: GAP Score Calculator ─────────────────────────────────────────────

calc_gap_score <- function(age, sex, fvc, dlco) {
  g <- ifelse(sex==0, 0, 1)                        # female=0 (1 pt), male=1 (0 pt)
  g_pts <- ifelse(sex==0, 1, 0)

  a_pts <- ifelse(age > 70, 2, ifelse(age > 65, 1, 0))

  p_fvc <- ifelse(fvc > 75, 0, ifelse(fvc >= 50, 1, 2))

  p_dlco <- if (is.na(dlco) || dlco > 55) 0 else
    ifelse(dlco >= 36, 1, ifelse(dlco >= 21, 2, 3))

  score <- g_pts + a_pts + p_fvc + p_dlco
  stage <- ifelse(score <= 3, "I", ifelse(score <= 5, "II", "III"))
  color <- ifelse(score <= 3, "green", ifelse(score <= 5, "yellow", "red"))
  mort_ref <- list(
    I   = c("5.6%", "10.9%", "16.3%"),
    II  = c("16.2%", "29.9%", "42.1%"),
    III = c("39.2%", "62.1%", "76.8%")
  )
  m <- mort_ref[[stage]]
  list(score=score, stage=stage, color=color, mort1yr=m[1], mort2yr=m[2], mort3yr=m[3])
}

## ─── LAUNCH ───────────────────────────────────────────────────────────────────

shinyApp(ui, server)
