## ============================================================
## Osteoporosis QSP — Shiny Interactive Dashboard
## 6 tabs: Patient Profile · PK · PD Markers · Clinical Endpoints
##         Scenario Comparison · Biomarker Panel
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)

# ----------------------------------------------------------------
# Inline mrgsolve model (lightweight version for Shiny reactivity)
# ----------------------------------------------------------------
library(mrgsolve)

model_code <- '
$PARAM
USE_ALN=0, USE_ZOL=0, USE_DMAB=0, USE_TPTD=0, USE_ROMO=0
MENO=1, GIO=0
ALN_FPP_IC50=50, ZOL_FPP_IC50=10
DMAB_ka=0.005, DMAB_CL=0.0019, DMAB_Vc=2.46
DMAB_RANKL_IC50=0.003, DMAB_Emax=0.98
TPTD_ka=5.0, TPTD_ke=1.4, TPTD_Vc=1500, TPTD_PTH1R_EC50=100
TPTD_Emax_OB=0.8, TPTD_Emax_OC=0.3
ROMO_ka=0.01, ROMO_CL=0.015, ROMO_Vc=3.5, ROMO_Vp=2.5, ROMO_Q=0.04
ROMO_SOST_IC50=0.15, ROMO_Emax_OB=1.2, ROMO_Emax_AC=0.6
kOB_in=0.02, kOB_out=0.004, kOC_in=0.015, kOC_out=0.01
kPREOB_in=0.015, kPREOC_in=0.012
RANKL_ss=10, OPG_ss=40, SCLER_ss=35, PTH_ss=35
E2_ss_pre=80, E2_ss_post=12, Ca_ss=9.4
MENO_RANKL_fold=1.6, MENO_OPG_fold=0.7, MENO_OB_fold=0.85
kRANKL_in=0.10, kRANKL_out=0.10, kOPG_in=0.08, kOPG_out=0.08
kSCLER_in=0.05, kSCLER_out=0.05
BMD_form=0.002, BMD_resor=0.002
kCTX_in=0.07, kCTX_out=0.20, kP1NP_in=1.5, kP1NP_out=0.03
GIO_OB_supp=0.5, GIO_OC_stim=1.3, GIO_RANKL_up=1.4
FxRisk_k=0.001, FxRisk_BMD_slope=2.0

$CMT ALN_BONE ZOL_BONE DMAB_SC DMAB_C TPTD_C
     ROMO_SC ROMO_C ROMO_P
     OB OC PREOB PREOC
     BMD CTX P1NP RANKL OPG SCLER PTH_sys E2 Ca_sys FRACT_RISK

$INIT
ALN_BONE=0, ZOL_BONE=0, DMAB_SC=0, DMAB_C=0, TPTD_C=0,
ROMO_SC=0,  ROMO_C=0,   ROMO_P=0,
OB=1, OC=1, PREOB=2, PREOC=2,
BMD=1, CTX=0.35, P1NP=50,
RANKL=10, OPG=40, SCLER=35,
PTH_sys=35, E2=80, Ca_sys=9.4, FRACT_RISK=0

$ODE
double ALN_EFF  = (USE_ALN==1)  ? ALN_BONE/(ALN_BONE  + ALN_FPP_IC50)  : 0;
double ZOL_EFF  = (USE_ZOL==1)  ? ZOL_BONE/(ZOL_BONE  + ZOL_FPP_IC50)  : 0;
dxdt_ALN_BONE = -0.0000004*ALN_BONE;
dxdt_ZOL_BONE = -0.0001*ZOL_BONE;
double DMAB_abs = (USE_DMAB==1) ? DMAB_ka*DMAB_SC : 0;
dxdt_DMAB_SC = -DMAB_abs;
dxdt_DMAB_C  = DMAB_abs/DMAB_Vc - (DMAB_CL/DMAB_Vc)*DMAB_C;
double DMAB_inh = (USE_DMAB==1) ? DMAB_Emax*DMAB_C/(DMAB_C+DMAB_RANKL_IC50) : 0;
dxdt_TPTD_C  = -(TPTD_ke)*TPTD_C;
double TPTD_stim = (USE_TPTD==1) ? TPTD_C/(TPTD_C+TPTD_PTH1R_EC50) : 0;
dxdt_ROMO_SC = -ROMO_ka*ROMO_SC;
dxdt_ROMO_C  = ROMO_ka*ROMO_SC/ROMO_Vc-(ROMO_CL/ROMO_Vc)*ROMO_C
               -(ROMO_Q/ROMO_Vc)*ROMO_C+(ROMO_Q/ROMO_Vp)*ROMO_P;
dxdt_ROMO_P  = (ROMO_Q/ROMO_Vc)*ROMO_C-(ROMO_Q/ROMO_Vp)*ROMO_P;
double ROMO_inh_SOST=(USE_ROMO==1)?ROMO_Emax_OB*ROMO_C/(ROMO_C+ROMO_SOST_IC50):0;
double ROMO_inh_OC  =(USE_ROMO==1)?ROMO_Emax_AC*ROMO_C/(ROMO_C+ROMO_SOST_IC50):0;
double E2_target=(MENO==1)?E2_ss_post:E2_ss_pre;
dxdt_E2 = 0.1*(E2_target-E2);
dxdt_PTH_sys = 0.5*(PTH_ss*(Ca_ss/Ca_sys)-PTH_sys);
dxdt_Ca_sys  = -0.2*(Ca_sys-Ca_ss)+0.05*(OC-1)+0.01*(PTH_sys-PTH_ss);
double MENO_R=(MENO==1)?MENO_RANKL_fold:1; double GIO_R=(GIO==1)?GIO_RANKL_up:1;
double PTH_Rup=1+0.3*(PTH_sys/PTH_ss-1);
dxdt_RANKL=kRANKL_in*RANKL_ss*MENO_R*GIO_R*PTH_Rup/(E2/E2_ss_pre+0.001)
           -kRANKL_out*RANKL*(1-DMAB_inh);
double OPG_Mf=(MENO==1)?MENO_OPG_fold:1; double OPG_Gf=(GIO==1)?0.7:1;
dxdt_OPG = kOPG_in*OPG_ss*OPG_Mf*OPG_Gf*(E2/E2_ss_pre)*(1+0.5*ROMO_inh_OC)
           -kOPG_out*OPG;
double SOST_pi=(USE_TPTD==1)?(1-0.3*TPTD_stim):1;
dxdt_SCLER = kSCLER_in*SCLER_ss*SOST_pi - kSCLER_out*SCLER*(1-ROMO_inh_SOST);
double Wnt=1/(1+SCLER/SCLER_ss);
double TPTD_OBs=(USE_TPTD==1)?(1+TPTD_Emax_OB*TPTD_stim):1;
double GIO_Obf=(GIO==1)?GIO_OB_supp:1; double MENO_Obf=(MENO==1)?MENO_OB_fold:1;
dxdt_PREOB=kPREOB_in*2.0*Wnt*TPTD_OBs*(1+ROMO_inh_SOST)*GIO_Obf*MENO_Obf-kOB_in*PREOB;
dxdt_OB   =kOB_in*PREOB-kOB_out*OB;
double RROP=RANKL/(OPG+0.001)/(RANKL_ss/OPG_ss);
double GIO_Ocf=(GIO==1)?GIO_OC_stim:1;
double TPTD_OCs=(USE_TPTD==1)?(1+TPTD_Emax_OC*TPTD_stim):1;
dxdt_PREOC=kPREOC_in*2.0*RROP*GIO_Ocf*TPTD_OCs-kOC_in*PREOC;
dxdt_OC   =kOC_in*PREOC-(kOC_out*(1+(ALN_EFF+ZOL_EFF)*0.8+DMAB_inh+ROMO_inh_OC))*OC;
dxdt_BMD  =BMD_form*OB-BMD_resor*OC;
dxdt_CTX  =kCTX_in*0.35*OC-kCTX_out*CTX;
dxdt_P1NP =kP1NP_in*OB-kP1NP_out*P1NP;
dxdt_FRACT_RISK=FxRisk_k*exp(FxRisk_BMD_slope*(1-BMD)/0.1)*(1-FRACT_RISK);

$TABLE
capture BMD_gcm2   = BMD*0.85;
capture Tscore     = (BMD_gcm2-0.955)/0.120;
capture CTX_out    = CTX;
capture P1NP_out   = P1NP;
capture RANKL_out  = RANKL;
capture OPG_out    = OPG;
capture SCLER_out  = SCLER;
capture OB_out     = OB;
capture OC_out     = OC;
capture PTH_out    = PTH_sys;
capture E2_out     = E2;
capture Ca_out     = Ca_sys;
capture FxR        = 1-exp(-FxRisk_k*exp(FxRisk_BMD_slope*(1-BMD)/0.1)*10*365*24);
capture DMAB_Cout  = DMAB_C;
capture TPTD_Cout  = TPTD_C;
capture ROMO_Cout  = ROMO_C;
'

mod_shiny <- mrgsolve::mcode("op_shiny", model_code, quiet = TRUE)

# ----------------------------------------------------------------
# Helper: build dosing events
# ----------------------------------------------------------------
build_ev <- function(use_aln, use_zol, use_dmab, use_tptd, use_romo, dur_yr) {
  h_total <- dur_yr * 365 * 24
  evs <- list()
  if (use_aln)  evs$aln  <- ev(cmt="ALN_BONE", amt=1500,  time=seq(0, h_total-168, 168))
  if (use_zol)  evs$zol  <- ev(cmt="ZOL_BONE", amt=5000,  time=seq(0, h_total-8760, 8760))
  if (use_dmab) evs$dmab <- ev(cmt="DMAB_SC",  amt=60,    time=seq(0, h_total-4380, 4380))
  if (use_tptd) {
    daily_rate <- 20000 / 24
    evs$tptd <- ev(cmt="TPTD_C", amt=daily_rate*24*dur_yr*365,
                   rate=daily_rate, time=0, tinf=h_total)
  }
  if (use_romo) {
    mo_times <- seq(0, min(h_total-720, 11*720), 720)
    evs$romo <- ev(cmt="ROMO_SC", amt=210, time=mo_times)
  }
  if (length(evs) == 0) return(ev(time=0, amt=0, cmt=1))
  Reduce(c, evs)
}

# ----------------------------------------------------------------
# Run one scenario
# ----------------------------------------------------------------
run_sim <- function(params_list, ev_obj, init_list = list(), dur_yr = 3) {
  out <- mod_shiny %>%
    param(params_list) %>%
    init(init_list) %>%
    ev(ev_obj) %>%
    mrgsim(end = dur_yr*365*24, delta = 12) %>%
    as_tibble() %>%
    mutate(time_yr = time/(365*24))
  out
}

# ----------------------------------------------------------------
# UI
# ----------------------------------------------------------------
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Osteoporosis QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("환자 프로파일",    tabName = "profile",   icon = icon("user")),
      menuItem("약동학 (PK)",      tabName = "pk",        icon = icon("flask")),
      menuItem("PD 주요 지표",     tabName = "pd",        icon = icon("chart-line")),
      menuItem("임상 엔드포인트",  tabName = "endpoints", icon = icon("heartbeat")),
      menuItem("시나리오 비교",    tabName = "scenarios", icon = icon("sliders-h")),
      menuItem("바이오마커 패널",  tabName = "biomarkers",icon = icon("vial"))
    ),

    hr(),
    h5("  환자 설정", style="padding-left:10px; font-weight:bold;"),
    selectInput("meno_status", "Menopause Status:",
                choices = c("Postmenopausal (Default)" = 1, "Premenopausal" = 0), selected = 1),
    checkboxInput("gio", "Glucocorticoid-Induced OP", FALSE),
    sliderInput("dur_yr", "Simulation Duration (years):", min=1, max=5, value=3, step=0.5),

    hr(),
    h5("  치료 선택", style="padding-left:10px; font-weight:bold;"),
    checkboxInput("use_aln",  "Alendronate 70mg/wk",     FALSE),
    checkboxInput("use_zol",  "Zoledronate 5mg/yr (IV)",  FALSE),
    checkboxInput("use_dmab", "Denosumab 60mg/6mo",       FALSE),
    checkboxInput("use_tptd", "Teriparatide 20µg/day",    FALSE),
    checkboxInput("use_romo", "Romosozumab 210mg/mo",     FALSE),
    actionButton("run_btn", "Run Simulation", class="btn-primary btn-block",
                 icon=icon("play"))
  ),

  dashboardBody(
    tabItems(

      # ---- Tab 1: Patient Profile ----
      tabItem("profile",
        fluidRow(
          valueBoxOutput("vb_bmd",  width=3),
          valueBoxOutput("vb_tscore", width=3),
          valueBoxOutput("vb_ctx",  width=3),
          valueBoxOutput("vb_p1np", width=3)
        ),
        fluidRow(
          box(title="Disease Overview", width=6, solidHeader=TRUE, status="primary",
            HTML("
              <h4>골다공증 (Osteoporosis)</h4>
              <p><b>정의:</b> 골밀도 감소와 뼈 미세구조 악화로 골절 위험이 증가하는 대사성 골질환.
              T-score ≤ −2.5 (WHO 기준).</p>
              <p><b>유병률:</b> 전세계 2억명 이상, 50세 이상 여성의 30%, 남성의 8%.</p>
              <p><b>핵심 기전:</b></p>
              <ul>
                <li>RANK/RANKL/OPG 축: 파골세포 분화 조절</li>
                <li>Wnt/β-Catenin + Sclerostin: 조골세포 조절</li>
                <li>에스트로겐 결핍: 폐경 후 골소실 가속화</li>
                <li>PTH·비타민D·칼슘 항상성</li>
              </ul>
              <p><b>치료 목표:</b> 골절 위험 감소 (척추 50–70%, 고관절 30–40%)</p>
            ")
          ),
          box(title="WHO T-score Classification", width=3, solidHeader=TRUE,
            DT::dataTableOutput("tbl_who")
          ),
          box(title="Drug Mechanism Overview", width=3, solidHeader=TRUE,
            DT::dataTableOutput("tbl_drugs")
          )
        ),
        fluidRow(
          box(title="BMD Trajectory (Current Settings)", width=12, status="info",
              plotlyOutput("plot_bmd_profile", height="300px"))
        )
      ),

      # ---- Tab 2: PK ----
      tabItem("pk",
        fluidRow(
          box(title="Drug Plasma Concentrations", width=8, status="primary",
              plotlyOutput("plot_pk_conc", height="400px")),
          box(title="PK Parameters Summary", width=4, solidHeader=TRUE,
              DT::dataTableOutput("tbl_pk_params"))
        ),
        fluidRow(
          box(title="Bone Depot Kinetics (Bisphosphonates)", width=6, status="warning",
              plotlyOutput("plot_bone_depot", height="300px")),
          box(title="SC Depot Kinetics (Denosumab/Romosozumab)", width=6, status="success",
              plotlyOutput("plot_sc_depot", height="300px"))
        )
      ),

      # ---- Tab 3: PD ----
      tabItem("pd",
        fluidRow(
          box(title="Osteoblast / Osteoclast Activity", width=6, status="primary",
              plotlyOutput("plot_cells", height="350px")),
          box(title="RANKL / OPG / Sclerostin", width=6, status="info",
              plotlyOutput("plot_mediators", height="350px"))
        ),
        fluidRow(
          box(title="Wnt Signaling & Bone Formation Drive", width=6, status="success",
              plotlyOutput("plot_wnt", height="300px")),
          box(title="Calcium & PTH Homeostasis", width=6, status="warning",
              plotlyOutput("plot_cahormone", height="300px"))
        )
      ),

      # ---- Tab 4: Clinical Endpoints ----
      tabItem("endpoints",
        fluidRow(
          box(title="BMD — Lumbar Spine & T-score", width=6, status="primary",
              plotlyOutput("plot_bmd_ts", height="350px")),
          box(title="Bone Turnover Markers (CTX & P1NP)", width=6, status="info",
              plotlyOutput("plot_btm", height="350px"))
        ),
        fluidRow(
          box(title="Fracture Risk Accumulation", width=6, status="danger",
              plotlyOutput("plot_fx_risk", height="300px")),
          box(title="Key Timepoints Summary", width=6, solidHeader=TRUE,
              DT::dataTableOutput("tbl_timepoints"))
        )
      ),

      # ---- Tab 5: Scenario Comparison ----
      tabItem("scenarios",
        fluidRow(
          box(width=12,
            selectInput("scen_endpoint", "Select Endpoint:",
              choices = c("BMD (g/cm²)"="BMD_gcm2", "T-score"="Tscore",
                          "CTX (ng/mL)"="CTX_out",  "P1NP (µg/L)"="P1NP_out",
                          "Osteoblast (rel.)"="OB_out", "Osteoclast (rel.)"="OC_out",
                          "Sclerostin (pmol/L)"="SCLER_out",
                          "10-yr Fx Risk"="FxR"),
              selected = "BMD_gcm2")
          )
        ),
        fluidRow(
          box(title="All Scenarios — Selected Endpoint", width=8, status="primary",
              plotlyOutput("plot_scenarios", height="400px")),
          box(title="3-Year Comparison Table", width=4, solidHeader=TRUE,
              DT::dataTableOutput("tbl_scenarios"))
        )
      ),

      # ---- Tab 6: Biomarker Panel ----
      tabItem("biomarkers",
        fluidRow(
          box(title="Biomarker Status Panel", width=4, solidHeader=TRUE, status="info",
              DT::dataTableOutput("tbl_biomarkers")),
          box(title="CTX & P1NP Over Time", width=8, status="primary",
              plotlyOutput("plot_btm_detail", height="350px"))
        ),
        fluidRow(
          box(title="Estrogen & PTH Dynamics", width=6, status="warning",
              plotlyOutput("plot_hormones", height="300px")),
          box(title="Biomarker Interpretation Guide", width=6, solidHeader=TRUE,
              DT::dataTableOutput("tbl_biomarker_guide"))
        )
      )
    )
  )
)

# ----------------------------------------------------------------
# SERVER
# ----------------------------------------------------------------
server <- function(input, output, session) {

  # ---- Static reference tables ----
  who_tbl <- data.frame(
    Classification = c("Normal", "Osteopenia", "Osteoporosis", "Severe OP"),
    `T-score` = c("≥ −1.0", "−1.0 to −2.5", "≤ −2.5", "≤ −2.5 + fracture"),
    Action = c("Monitor", "Lifestyle + Ca/VitD", "Pharmacotherapy", "Urgent Rx + Fall prevention"),
    check.names = FALSE
  )
  drug_tbl <- data.frame(
    Drug = c("Alendronate","Zoledronate","Denosumab","Teriparatide","Romosozumab"),
    Class = c("N-BP","N-BP","Anti-RANKL mAb","PTH analog","Anti-Sclerostin mAb"),
    Action = c("Anti-resorptive","Anti-resorptive","Anti-resorptive","Anabolic","Dual"),
    `Vertebral Fx RR` = c("47%","70%","68%","65%","73%"),
    check.names = FALSE
  )
  pk_params <- data.frame(
    Drug = c("Alendronate","Zoledronate","Denosumab","Teriparatide","Romosozumab"),
    Route = c("PO","IV","SC","SC","SC"),
    `F (%)` = c("0.7","~100","~62","~95","~81"),
    `t½ (bone/plasma)` = c(">10yr",">1yr","~26d","~1h","~12d"),
    `Dose & Freq` = c("70mg/wk","5mg/yr","60mg/6mo","20µg/d","210mg/mo"),
    check.names = FALSE
  )

  output$tbl_who  <- DT::renderDataTable(who_tbl,  options=list(dom='t', pageLength=10), rownames=FALSE)
  output$tbl_drugs<- DT::renderDataTable(drug_tbl, options=list(dom='t', pageLength=10), rownames=FALSE)
  output$tbl_pk_params <- DT::renderDataTable(pk_params, options=list(dom='t'), rownames=FALSE)

  # ---- Reactive simulation ----
  sim_data <- eventReactive(input$run_btn, {
    withProgress(message="Simulating...", {
      meno_val <- as.integer(input$meno_status)
      gio_val  <- as.integer(input$gio)
      dur      <- input$dur_yr

      p_list <- list(
        USE_ALN=as.integer(input$use_aln), USE_ZOL=as.integer(input$use_zol),
        USE_DMAB=as.integer(input$use_dmab), USE_TPTD=as.integer(input$use_tptd),
        USE_ROMO=as.integer(input$use_romo), MENO=meno_val, GIO=gio_val
      )
      init_list <- if (meno_val == 1) {
        list(E2=12, RANKL=16, OPG=28, OC=1.4, OB=0.9, CTX=0.55, P1NP=38, BMD=0.96)
      } else {
        list()
      }
      ev_obj <- build_ev(input$use_aln, input$use_zol, input$use_dmab,
                         input$use_tptd, input$use_romo, dur)
      run_sim(p_list, ev_obj, init_list, dur)
    })
  }, ignoreNULL=FALSE)

  # ---- Scenario comparison data ----
  scenario_data <- reactive({
    dur <- input$dur_yr
    meno_val <- as.integer(input$meno_status)
    gio_val  <- as.integer(input$gio)
    init_meno <- if (meno_val == 1) list(E2=12, RANKL=16, OPG=28, OC=1.4, OB=0.9, CTX=0.55, P1NP=38, BMD=0.96) else list()

    scenarios_def <- list(
      list(name="S1: Untreated",       aln=0,zol=0,dmab=0,tptd=0,romo=0),
      list(name="S2: Alendronate",     aln=1,zol=0,dmab=0,tptd=0,romo=0),
      list(name="S3: Zoledronate",     aln=0,zol=1,dmab=0,tptd=0,romo=0),
      list(name="S4: Denosumab",       aln=0,zol=0,dmab=1,tptd=0,romo=0),
      list(name="S5: Teriparatide",    aln=0,zol=0,dmab=0,tptd=1,romo=0),
      list(name="S6: Romo→Denosumab",  aln=0,zol=0,dmab=1,tptd=0,romo=1)
    )
    all_sims <- lapply(scenarios_def, function(s) {
      p <- list(USE_ALN=s$aln, USE_ZOL=s$zol, USE_DMAB=s$dmab,
                USE_TPTD=s$tptd, USE_ROMO=s$romo, MENO=meno_val, GIO=gio_val)
      ev_obj <- build_ev(s$aln==1, s$zol==1, s$dmab==1, s$tptd==1, s$romo==1, dur)
      out <- run_sim(p, ev_obj, init_meno, dur)
      out$scenario <- s$name
      out
    })
    bind_rows(all_sims)
  })

  # ---- Value boxes ----
  output$vb_bmd <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    valueBox(round(last$BMD_gcm2, 3), "BMD (g/cm²) at end", icon=icon("bone"), color="blue")
  })
  output$vb_tscore <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    col <- if (last$Tscore > -1) "green" else if (last$Tscore > -2.5) "yellow" else "red"
    valueBox(round(last$Tscore, 2), "T-score at end", icon=icon("chart-bar"), color=col)
  })
  output$vb_ctx <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    valueBox(round(last$CTX_out, 3), "CTX (ng/mL) at end", icon=icon("vial"), color="orange")
  })
  output$vb_p1np <- renderValueBox({
    d <- sim_data(); last <- tail(d, 1)
    valueBox(round(last$P1NP_out, 1), "P1NP (µg/L) at end", icon=icon("flask"), color="green")
  })

  # ---- Tab 1: BMD profile ----
  output$plot_bmd_profile <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr, y=~BMD_gcm2, type='scatter', mode='lines',
            line=list(color='steelblue', width=2)) %>%
      add_lines(y=~Tscore, name="T-score", yaxis="y2",
                line=list(color='firebrick', dash='dash')) %>%
      layout(title="BMD & T-score Over Time",
             xaxis=list(title="Years"),
             yaxis=list(title="BMD (g/cm²)"),
             yaxis2=list(title="T-score", overlaying="y", side="right"),
             legend=list(orientation="h"))
  })

  # ---- Tab 2: PK ----
  output$plot_pk_conc <- renderPlotly({
    d <- sim_data()
    p <- plot_ly(d, x=~time_yr, type='scatter', mode='lines')
    if (any(d$DMAB_Cout > 0.001))
      p <- p %>% add_lines(y=~DMAB_Cout, name="Denosumab (µg/mL)")
    if (any(d$TPTD_Cout > 1))
      p <- p %>% add_lines(y=~TPTD_Cout/1000, name="Teriparatide (ng/mL×1000 scaled)")
    if (any(d$ROMO_Cout > 0.001))
      p <- p %>% add_lines(y=~ROMO_Cout, name="Romosozumab (µg/mL)")
    p %>% layout(xaxis=list(title="Years"), yaxis=list(title="Concentration"),
                 legend=list(orientation="h"))
  })
  output$plot_bone_depot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr, type='scatter', mode='lines') %>%
      add_lines(y=~ALN_BONE, name="Alendronate Bone Depot") %>%
      add_lines(y=~ZOL_BONE, name="Zoledronate Bone Depot") %>%
      layout(xaxis=list(title="Years"), yaxis=list(title="Bone Depot (nmol/g)"),
             legend=list(orientation="h"))
  })
  output$plot_sc_depot <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr, type='scatter', mode='lines') %>%
      add_lines(y=~DMAB_SC, name="Denosumab SC Depot (µg)") %>%
      add_lines(y=~ROMO_SC, name="Romosozumab SC Depot (µg)") %>%
      layout(xaxis=list(title="Years"), yaxis=list(title="SC Depot (µg)"),
             legend=list(orientation="h"))
  })

  # ---- Tab 3: PD ----
  output$plot_cells <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr, type='scatter', mode='lines') %>%
      add_lines(y=~OB_out, name="Osteoblast (rel.)", line=list(color='green')) %>%
      add_lines(y=~OC_out, name="Osteoclast (rel.)", line=list(color='red')) %>%
      layout(xaxis=list(title="Years"), yaxis=list(title="Relative Activity"),
             shapes=list(list(type='line',x0=0,x1=max(d$time_yr),y0=1,y1=1,
                              line=list(dash='dot',color='gray'))),
             legend=list(orientation="h"))
  })
  output$plot_mediators <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr, type='scatter', mode='lines') %>%
      add_lines(y=~RANKL_out, name="RANKL (pmol/L)", line=list(color='red')) %>%
      add_lines(y=~OPG_out,   name="OPG (pmol/L)",   line=list(color='green')) %>%
      add_lines(y=~SCLER_out, name="Sclerostin (pmol/L)", line=list(color='purple')) %>%
      layout(xaxis=list(title="Years"), yaxis=list(title="Concentration (pmol/L)"),
             legend=list(orientation="h"))
  })
  output$plot_wnt <- renderPlotly({
    d <- sim_data() %>% mutate(Wnt_index = 1/(1+SCLER_out/35))
    plot_ly(d, x=~time_yr, y=~Wnt_index, type='scatter', mode='lines',
            line=list(color='darkgreen')) %>%
      layout(title="Wnt Signaling Activity (1/(1+SOST/baseline))",
             xaxis=list(title="Years"), yaxis=list(title="Wnt Activity Index (0–1)"))
  })
  output$plot_cahormone <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr, type='scatter', mode='lines') %>%
      add_lines(y=~Ca_out,  name="Ca²⁺ (mg/dL)", line=list(color='orange')) %>%
      add_lines(y=~PTH_out/10, name="PTH/10 (pg/mL)", line=list(color='blue')) %>%
      add_lines(y=~E2_out/10,  name="E2/10 (pg/mL)", line=list(color='pink')) %>%
      layout(xaxis=list(title="Years"), yaxis=list(title="Value (scaled)"),
             legend=list(orientation="h"))
  })

  # ---- Tab 4: Endpoints ----
  output$plot_bmd_ts <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr, type='scatter', mode='lines') %>%
      add_lines(y=~BMD_gcm2, name="BMD (g/cm²)", line=list(color='steelblue')) %>%
      add_lines(y=~Tscore,   name="T-score", yaxis="y2", line=list(color='darkred',dash='dash')) %>%
      layout(xaxis=list(title="Years"),
             yaxis=list(title="BMD (g/cm²)"),
             yaxis2=list(title="T-score", overlaying="y", side="right"),
             shapes=list(list(type='line',x0=0,x1=max(d$time_yr),y0=-2.5,y1=-2.5,
                              line=list(dash='dot',color='red'),yref='y2')),
             legend=list(orientation="h"))
  })
  output$plot_btm <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr, type='scatter', mode='lines') %>%
      add_lines(y=~CTX_out,  name="CTX (ng/mL)",  line=list(color='orange')) %>%
      add_lines(y=~P1NP_out, name="P1NP (µg/L)", yaxis="y2", line=list(color='green')) %>%
      layout(xaxis=list(title="Years"),
             yaxis=list(title="CTX (ng/mL)"),
             yaxis2=list(title="P1NP (µg/L)", overlaying="y", side="right"),
             legend=list(orientation="h"))
  })
  output$plot_fx_risk <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr, y=~FxR*100, type='scatter', mode='lines',
            line=list(color='darkred', width=2), fill='tozeroy', fillcolor='rgba(200,0,0,0.1)') %>%
      layout(title="10-Year Fracture Probability",
             xaxis=list(title="Years"), yaxis=list(title="Probability (%)"))
  })
  output$tbl_timepoints <- DT::renderDataTable({
    d <- sim_data()
    tps <- c(0, 0.5, 1, 2, 3)
    tps <- tps[tps <= input$dur_yr]
    rows <- lapply(tps, function(t) {
      idx <- which.min(abs(d$time_yr - t))
      r <- d[idx, ]
      data.frame(
        Year=t, `BMD (g/cm²)`=round(r$BMD_gcm2,3), `T-score`=round(r$Tscore,2),
        `CTX (ng/mL)`=round(r$CTX_out,3), `P1NP (µg/L)`=round(r$P1NP_out,1),
        `OB (rel.)`=round(r$OB_out,3), `OC (rel.)`=round(r$OC_out,3),
        `Fx Risk 10yr (%)`=round(r$FxR*100,1), check.names=FALSE
      )
    })
    bind_rows(rows)
  }, options=list(dom='t', pageLength=10), rownames=FALSE)

  # ---- Tab 5: Scenario Comparison ----
  output$plot_scenarios <- renderPlotly({
    d <- scenario_data()
    ep <- input$scen_endpoint
    if (!ep %in% names(d)) return(plotly_empty())
    d$ep_val <- d[[ep]]
    p <- plot_ly(d, x=~time_yr, y=~ep_val, color=~scenario,
                 type='scatter', mode='lines') %>%
      layout(xaxis=list(title="Years"), yaxis=list(title=ep),
             legend=list(orientation="h"))
    p
  })
  output$tbl_scenarios <- DT::renderDataTable({
    d <- scenario_data() %>%
      group_by(scenario) %>%
      filter(abs(time_yr - input$dur_yr) == min(abs(time_yr - input$dur_yr))) %>%
      slice(1) %>%
      select(scenario, BMD_gcm2, Tscore, CTX_out, P1NP_out, FxR) %>%
      mutate(across(where(is.numeric), ~round(.x, 3))) %>%
      rename(`Scenario`=scenario, `BMD`=BMD_gcm2, `T-score`=Tscore,
             `CTX`=CTX_out, `P1NP`=P1NP_out, `Fx Risk`=FxR)
    d
  }, options=list(dom='t', pageLength=10), rownames=FALSE)

  # ---- Tab 6: Biomarkers ----
  output$tbl_biomarkers <- DT::renderDataTable({
    d <- sim_data(); last <- tail(d, 1)
    df <- data.frame(
      Biomarker = c("BMD (g/cm²)","T-score","CTX (ng/mL)","P1NP (µg/L)",
                    "RANKL (pmol/L)","OPG (pmol/L)","Sclerostin (pmol/L)",
                    "Serum Ca²⁺ (mg/dL)","PTH (pg/mL)","Estradiol (pg/mL)"),
      Value = c(round(last$BMD_gcm2,3), round(last$Tscore,2), round(last$CTX_out,3),
                round(last$P1NP_out,1), round(last$RANKL_out,1), round(last$OPG_out,1),
                round(last$SCLER_out,1), round(last$Ca_out,2), round(last$PTH_out,1),
                round(last$E2_out,1)),
      Reference = c("0.85–1.05","≥ −1.0","0.1–0.5","25–80","3–20","2–80","18–50",
                    "8.5–10.5","15–65","20–400 (pre) / 5–30 (post)"),
      Status = c(
        ifelse(last$BMD_gcm2 >= 0.80, "OK", "↓LOW"),
        ifelse(last$Tscore >= -2.5, ifelse(last$Tscore >= -1, "Normal","Osteopenia"), "OP"),
        ifelse(last$CTX_out <= 0.55, "OK", "↑HIGH"),
        ifelse(last$P1NP_out >= 25, "OK", "↓LOW"),
        "Info", "Info", "Info",
        ifelse(abs(last$Ca_out - 9.4) < 0.5, "OK", "Abnormal"),
        ifelse(last$PTH_out < 65, "OK", "↑HIGH"),
        ifelse(last$E2_out >= 20 | as.integer(input$meno_status)==1, "Info", "↓LOW")
      )
    )
    df
  }, options=list(dom='t', pageLength=15), rownames=FALSE)

  output$plot_btm_detail <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr, type='scatter', mode='lines') %>%
      add_lines(y=~CTX_out, name="CTX (ng/mL)", line=list(color='orange',width=2)) %>%
      add_lines(y=~P1NP_out/100, name="P1NP/100 (µg/L)", line=list(color='green',width=2,dash='dash')) %>%
      layout(title="Bone Turnover Markers",
             xaxis=list(title="Years"),
             yaxis=list(title="CTX (ng/mL) | P1NP/100"),
             legend=list(orientation="h"))
  })
  output$plot_hormones <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x=~time_yr, type='scatter', mode='lines') %>%
      add_lines(y=~E2_out,  name="Estradiol (pg/mL)", line=list(color='deeppink')) %>%
      add_lines(y=~PTH_out, name="PTH (pg/mL)",        line=list(color='royalblue')) %>%
      layout(xaxis=list(title="Years"), yaxis=list(title="pg/mL"),
             legend=list(orientation="h"))
  })
  output$tbl_biomarker_guide <- DT::renderDataTable({
    data.frame(
      Marker = c("CTX","P1NP","RANKL","OPG","Sclerostin","T-score","PTH","Ca²⁺","E2","FxRisk"),
      Meaning = c("Bone resorption","Bone formation","OC differentiation signal",
                  "RANKL decoy receptor","Wnt inhibitor (osteocyte)",
                  "Fracture risk threshold","Calcium regulator","Mineral balance",
                  "OB/OC balance regulator","10-year fragility fracture probability"),
      `Treatment Response` = c("BP/DMAB↓","TPTD/ROMO↑","DMAB blocks","E2/Treg↑","ROMO neutralizes",
                                "All agents↑","Monitor (secondary OP)","Ca+VitD supplement",
                                "HRT (selected)","Primary endpoint"),
      check.names = FALSE
    )
  }, options=list(dom='t', pageLength=12), rownames=FALSE)
}

shinyApp(ui, server)
