###############################################################################
# Chronic Pyelonephritis (CPN) — Interactive Shiny QSP Dashboard
# 만성 신우신염 인터랙티브 시뮬레이션 대시보드
# Tabs: 환자 프로파일 · PK · 세균 동태 · 신기능 · 시나리오 비교 · 바이오마커
###############################################################################

library(shiny)
library(mrgsolve)
library(tidyverse)
library(ggplot2)
library(plotly)
library(DT)
library(shinydashboard)

# --------------------------------------------------------------------------- #
# mrgsolve model (inline) — same as cpn_mrgsolve_model.R but streamlined
# --------------------------------------------------------------------------- #
cpn_code <- '
$PARAM
ka_cipro=1.5,F_cipro=0.70,Vc_cipro=140,Vp_cipro=100,k12_cipro=0.12,
k21_cipro=0.08,CL_cipro=25,fu_cipro=0.65,
ka_tmp=0.80,F_tmp=0.95,Vc_tmp=85,CL_tmp=3.5,
ka_nit=2.0,F_nit=0.75,ke_nit=0.4,fuUr_nit=0.4,
MIC_cipro=0.25,Emax_cipro=3.5,EC50_cipro=2.0,Hill_cipro=1.5,
MIC_tmp=2.0,Emax_tmp=2.5,EC50_tmp=1.5,Hill_tmp=2.0,
MIC_nit=32,Emax_nit=2.0,EC50_nit=64,Hill_nit=2.0,
kgrow=0.15,Bmax=9.0,kdeath=0.02,kbiofilm=0.04,kdisbio=0.01,
phi_biofilm=0.3,krecurr=0.0005,Resist_F=1.0,
kNin=2.0,kNout=0.3,kMin=0.5,kMout=0.05,
kIL6=0.8,kIL6_deg=0.2,kTGFb1=0.3,kTGFb1_deg=0.1,
neutKill=0.5,macKill=0.3,
kCol=0.001,kCol_deg=0.0002,kScar=0.0005,ScarMax=0.9,ColMax=5.0,
GFR0=100,kGFR_scar=5.0,kGFR_inf=0.1,kGFR_return=0.01,GFR_floor=5.0,
use_cipro=0,use_tmp=0,use_nit=0

$INIT
Cipro_gut=0,Cipro_C=0,Cipro_P=0,TMP_gut=0,TMP_C=0,NIT_gut=0,NIT_urine=0,
Bacteria=6,Biofilm=0.1,Neutrophil=1,Macrophage=1,IL6=1,TGFb1=1,
Collagen=1,RenalScar=0,GFR=100

$ODE
dxdt_Cipro_gut = -ka_cipro * Cipro_gut;
double Cp_cipro = Cipro_C;
dxdt_Cipro_C = ka_cipro*Cipro_gut/Vc_cipro
               -(CL_cipro/Vc_cipro+k12_cipro)*Cipro_C+k21_cipro*(Cipro_P/Vp_cipro);
dxdt_Cipro_P = k12_cipro*Vc_cipro*Cipro_C - k21_cipro*Cipro_P;
dxdt_TMP_gut = -ka_tmp*TMP_gut;
dxdt_TMP_C   = ka_tmp*TMP_gut/Vc_tmp - (CL_tmp/Vc_tmp)*TMP_C;
dxdt_NIT_gut   = -ka_nit*NIT_gut;
dxdt_NIT_urine = ka_nit*NIT_gut*fuUr_nit - ke_nit*NIT_urine;
double CiproFree = fu_cipro*Cp_cipro;
double MIC_c_eff = MIC_cipro*Resist_F;
double fAUC_MIC  = (CiproFree>0) ? CiproFree/MIC_c_eff : 0.0;
double Kill_cip  = (use_cipro>0.5) ? Emax_cipro*pow(fAUC_MIC,Hill_cipro)/
                   (pow(EC50_cipro,Hill_cipro)+pow(fAUC_MIC,Hill_cipro)) : 0.0;
double Kill_tmp_v = (use_tmp>0.5) ? Emax_tmp*pow(TMP_C,Hill_tmp)/
                    (pow(EC50_tmp*MIC_tmp*Resist_F,Hill_tmp)+pow(TMP_C,Hill_tmp)) : 0.0;
double Kill_nit_v = (use_nit>0.5) ? Emax_nit*pow(NIT_urine,Hill_nit)/
                    (pow(EC50_nit,Hill_nit)+pow(NIT_urine,Hill_nit)) : 0.0;
double Kill_drug = (Kill_cip+Kill_tmp_v+Kill_nit_v)*(1.0-Biofilm*(1.0-phi_biofilm));
double BactNorm = (Bacteria>0) ? Bacteria/Bmax : 0.0;
double GrowthRate = kgrow*(1.0-BactNorm);
double ImmuneKill = neutKill*(Neutrophil-1.0)+macKill*(Macrophage-1.0);
if(ImmuneKill<0.0) ImmuneKill=0.0;
double dBact = GrowthRate-Kill_drug-kdeath-ImmuneKill+krecurr;
dxdt_Bacteria = (Bacteria<=0.0 && dBact<0.0) ? 0.0 : dBact;
double dBiofilm = kbiofilm*Biofilm*(1.0-Biofilm)*BactNorm - kdisbio*Biofilm*Kill_drug;
if(Biofilm<=0.0 && dBiofilm<0.0) dBiofilm=0.0;
if(Biofilm>=1.0 && dBiofilm>0.0) dBiofilm=0.0;
dxdt_Biofilm = dBiofilm;
double BactSignal = (Bacteria>0.0) ? Bacteria/6.0 : 0.0;
dxdt_Neutrophil = kNin*BactSignal - kNout*Neutrophil;
dxdt_Macrophage = kMin*BactSignal+0.1*(Neutrophil-1.0) - kMout*Macrophage;
dxdt_IL6   = kIL6*BactSignal*(1.0+0.5*(Neutrophil-1.0)) - kIL6_deg*IL6;
double M2signal = (Macrophage>2.0) ? (Macrophage-2.0) : 0.0;
dxdt_TGFb1 = kTGFb1*(M2signal+0.5)+0.05*(IL6-1.0) - kTGFb1_deg*TGFb1;
dxdt_Collagen = kCol*TGFb1*(1.0-Collagen/ColMax) - kCol_deg*Collagen;
double ScarDrive = (Collagen>1.0) ? kScar*(Collagen-1.0) : 0.0;
dxdt_RenalScar  = ScarDrive*(1.0-RenalScar/ScarMax);
if(dxdt_RenalScar<0.0) dxdt_RenalScar=0.0;
double GFR_target = GFR0*(1.0-kGFR_scar/GFR0*RenalScar) - kGFR_inf*(IL6-1.0);
if(GFR_target<GFR_floor) GFR_target=GFR_floor;
dxdt_GFR = -kGFR_return*(GFR-GFR_target);
if(GFR<=GFR_floor && dxdt_GFR<0.0) dxdt_GFR=0.0;

$TABLE
double Creatinine  = 100.0/(GFR>5.0?GFR:5.0);
double CKD_Stage   = (GFR>=90)?1:(GFR>=60)?2:(GFR>=30)?3:(GFR>=15)?4:5;
double ScarPct     = RenalScar*100.0;
double fAUC_MIC_OB = CiproFree/(MIC_cipro*Resist_F);

$CAPTURE
Cipro_C,TMP_C,NIT_urine,Bacteria,Biofilm,Neutrophil,Macrophage,
IL6,TGFb1,Collagen,RenalScar,GFR,Creatinine,CKD_Stage,ScarPct,
fAUC_MIC_OB,Kill_drug
'

mod <- mcode("cpn_shiny", cpn_code, quiet = TRUE)

# --------------------------------------------------------------------------- #
# Helper: build events
# --------------------------------------------------------------------------- #
make_ev <- function(drug, amt_scaled, ii, ndays, start_h) {
  cmt_map <- c(cipro="Cipro_gut", tmp="TMP_gut", nit="NIT_gut")
  ev(amt  = amt_scaled,
     ii   = ii,
     addl = max(0, ndays * (24 / ii) - 1),
     cmt  = cmt_map[drug],
     time = start_h)
}

run_sim <- function(params_list, events_list, end_h, delta = 2) {
  m <- param(mod, params_list)
  evs <- if (length(events_list) > 0) do.call(c, events_list) else NULL
  if (!is.null(evs)) {
    out <- mrgsim(m, events = evs, end = end_h, delta = delta)
  } else {
    out <- mrgsim(m, end = end_h, delta = delta)
  }
  as.data.frame(out)
}

# --------------------------------------------------------------------------- #
# UI
# --------------------------------------------------------------------------- #
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "CPN QSP Dashboard — 만성 신우신염"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("① 환자 프로파일",   tabName = "patient",   icon = icon("user-circle")),
      menuItem("② PK 모니터링",     tabName = "pk",        icon = icon("chart-line")),
      menuItem("③ 세균 동태",       tabName = "bacteria",  icon = icon("bacteria")),
      menuItem("④ 신기능 (GFR)",    tabName = "renal",     icon = icon("heartbeat")),
      menuItem("⑤ 시나리오 비교",   tabName = "scenario",  icon = icon("sliders-h")),
      menuItem("⑥ 바이오마커 패널", tabName = "biomarker", icon = icon("vials")),
      menuItem("⑦ 참고문헌",        tabName = "refs",      icon = icon("book-open"))
    ),
    hr(),
    tags$div(style = "padding:10px; font-size:11px; color:#ccc;",
      "CPN QSP v1.0",
      br(), "mrgsolve · 2026")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color:#f4f6f9; }
      .box-header .box-title { font-weight: bold; }
      .val-box-label { font-size:12px; }
    "))),

    tabItems(

      # ====== TAB 1: Patient Profile =========================================
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "환자 파라미터 설정", width = 4, status = "primary",
              solidHeader = TRUE,
              sliderInput("GFR0",    "기저 GFR (mL/min/1.73m²)",  15, 120, 90, 5),
              selectInput("VUR_risk","VUR 중증도 (구조적 위험)",
                          c("없음 (정상)"="low", "경도 I–II"="mild",
                            "중등도 III"="mod", "고도 IV–V"="high")),
              selectInput("Bact_init","초기 세균 부하 (log₁₀ CFU/g)",
                          c("경증 10^4"="4", "중등도 10^6"="6",
                            "중증 10^8"="8"), selected="6"),
              checkboxInput("hasDM",   "당뇨병 (DM)",       FALSE),
              checkboxInput("hasImm",  "면역억제 상태",      FALSE),
              checkboxInput("hasPreg", "임신 (위험 3×↑)",   FALSE),
              numericInput("Resist_F", "내성 계수 (1=감수성, 4=중등내성)", 1, 1, 8, 0.5)
          ),
          box(title = "항생제 처방 설정", width = 4, status = "warning",
              solidHeader = TRUE,
              checkboxInput("use_cipro", "Ciprofloxacin 사용",    FALSE),
              conditionalPanel("input.use_cipro",
                sliderInput("cipro_dose", "용량 (mg/dose)", 250, 750, 500, 50),
                sliderInput("cipro_days", "처방 기간 (days)", 3, 28, 14, 1)
              ),
              hr(),
              checkboxInput("use_tmp",   "TMP-SMX 사용",         FALSE),
              conditionalPanel("input.use_tmp",
                sliderInput("tmp_days", "TMP-SMX 기간 (days)", 3, 28, 14, 1)
              ),
              hr(),
              checkboxInput("use_nit",   "Nitrofurantoin 예방 사용", FALSE),
              conditionalPanel("input.use_nit",
                sliderInput("nit_months", "예방 투여 기간 (months)", 1, 12, 6, 1)
              )
          ),
          box(title = "시뮬레이션 설정", width = 4, status = "success",
              solidHeader = TRUE,
              sliderInput("sim_months", "시뮬레이션 기간 (개월)", 1, 24, 12, 1),
              actionButton("run_btn", "▶  시뮬레이션 실행",
                           class = "btn-primary btn-lg", width = "100%"),
              hr(),
              htmlOutput("patient_summary_box")
          )
        ),
        fluidRow(
          valueBoxOutput("vb_gfr",    width = 3),
          valueBoxOutput("vb_stage",  width = 3),
          valueBoxOutput("vb_scar",   width = 3),
          valueBoxOutput("vb_bact",   width = 3)
        )
      ),

      # ====== TAB 2: PK Monitoring ===========================================
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Ciprofloxacin 혈중 농도 (Cp vs time)", width = 8,
              status = "info", solidHeader = TRUE,
              plotlyOutput("plot_pk_cipro", height = "340px")),
          box(title = "PK/PD 지표 (fAUC/MIC)", width = 4,
              status = "info", solidHeader = TRUE,
              plotlyOutput("plot_fAUC", height = "340px"))
        ),
        fluidRow(
          box(title = "TMP 혈중 농도", width = 6,
              status = "primary", solidHeader = TRUE,
              plotlyOutput("plot_pk_tmp", height = "280px")),
          box(title = "Nitrofurantoin 요중 농도", width = 6,
              status = "primary", solidHeader = TRUE,
              plotlyOutput("plot_pk_nit", height = "280px"))
        ),
        fluidRow(
          box(title = "PK 파라미터 요약", width = 12, status = "warning",
              DT::dataTableOutput("pk_param_table"))
        )
      ),

      # ====== TAB 3: Bacterial Dynamics ======================================
      tabItem(tabName = "bacteria",
        fluidRow(
          box(title = "신장 내 세균 부하 (log₁₀ CFU/g)", width = 8,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_bacteria", height = "360px")),
          box(title = "생물막 (Biofilm) 분율", width = 4,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_biofilm", height = "360px"))
        ),
        fluidRow(
          box(title = "항균 효과 (Kill Rate)", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_kill", height = "280px")),
          box(title = "면역세포 동원 (Neutrophil / Macrophage)", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_immune", height = "280px"))
        )
      ),

      # ====== TAB 4: Renal Function ==========================================
      tabItem(tabName = "renal",
        fluidRow(
          box(title = "GFR 변화 추이", width = 7,
              status = "primary", solidHeader = TRUE,
              plotlyOutput("plot_gfr", height = "360px")),
          box(title = "혈청 크레아티닌", width = 5,
              status = "primary", solidHeader = TRUE,
              plotlyOutput("plot_creat", height = "360px"))
        ),
        fluidRow(
          box(title = "신장 반흔 진행 (Scar %)", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_scar", height = "280px")),
          box(title = "CKD 단계 전환", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_ckd", height = "280px"))
        )
      ),

      # ====== TAB 5: Scenario Comparison =====================================
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "5가지 항생제 전략 비교 (GFR)", width = 12,
              status = "success", solidHeader = TRUE,
              plotlyOutput("plot_scenario_gfr", height = "380px"))
        ),
        fluidRow(
          box(title = "세균 부하 비교", width = 6,
              status = "success", solidHeader = TRUE,
              plotlyOutput("plot_scenario_bact", height = "300px")),
          box(title = "신장 반흔 비교", width = 6,
              status = "success", solidHeader = TRUE,
              plotlyOutput("plot_scenario_scar", height = "300px"))
        ),
        fluidRow(
          box(title = "시나리오 비교 요약표", width = 12,
              DT::dataTableOutput("scenario_table"))
        )
      ),

      # ====== TAB 6: Biomarker Panel =========================================
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "염증 및 섬유화 마커 (IL-6, TGF-β1, Collagen)", width = 8,
              status = "info", solidHeader = TRUE,
              plotlyOutput("plot_biomarker", height = "380px")),
          box(title = "바이오마커 임상 의미", width = 4,
              status = "info", solidHeader = TRUE,
              tableOutput("biomarker_legend"))
        ),
        fluidRow(
          box(title = "요로패혈증 위험 지수 (Bacteria > 10^7.5)", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_sepsis_risk", height = "280px")),
          box(title = "Collagen 침착 vs Scar 형성", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_collagen_scar", height = "280px"))
        )
      ),

      # ====== TAB 7: References ==============================================
      tabItem(tabName = "refs",
        box(title = "주요 참고문헌 (만성 신우신염 QSP 모델)", width = 12,
            htmlOutput("ref_html"))
      )

    )  # end tabItems
  )
)

# --------------------------------------------------------------------------- #
# SERVER
# --------------------------------------------------------------------------- #
server <- function(input, output, session) {

  # ---- Reactive: run main simulation
  sim_data <- eventReactive(input$run_btn, {
    end_h <- input$sim_months * 30 * 24

    # Comorbidity risk factor
    risk_mult <- 1.0
    if (input$hasDM)   risk_mult <- risk_mult * 1.8
    if (input$hasImm)  risk_mult <- risk_mult * 2.5
    if (input$hasPreg) risk_mult <- risk_mult * 1.5

    vur_krecurr <- switch(input$VUR_risk,
      low  = 0.0005,
      mild = 0.001,
      mod  = 0.002,
      high = 0.004
    )

    params <- list(
      GFR0     = input$GFR0,
      Bacteria = as.numeric(input$Bact_init),
      Resist_F = input$Resist_F,
      krecurr  = vur_krecurr * risk_mult,
      use_cipro = as.numeric(input$use_cipro),
      use_tmp   = as.numeric(input$use_tmp),
      use_nit   = as.numeric(input$use_nit)
    )

    evs <- list()
    if (input$use_cipro) {
      evs[["cipro"]] <- make_ev("cipro",
        input$cipro_dose * 0.70, 12, input$cipro_days, 0)
    }
    if (input$use_tmp) {
      evs[["tmp"]] <- make_ev("tmp",
        160 * 0.95, 12, input$tmp_days, 0)
    }
    if (input$use_nit) {
      start_nit <- if (input$use_cipro) input$cipro_days * 24 else 0
      evs[["nit"]] <- make_ev("nit",
        100 * 0.75, 24, input$nit_months * 30, start_nit)
    }

    m <- param(mod, params)
    init(m, Bacteria = as.numeric(input$Bact_init))

    run_sim(params, evs, end_h, delta = 6)
  }, ignoreNULL = FALSE)

  # ---- Scenario comparison reactive (runs all 5 preset scenarios)
  scenario_data <- eventReactive(input$run_btn, {
    end_h <- input$sim_months * 30 * 24
    base_p <- list(GFR0 = input$GFR0, Resist_F = input$Resist_F)
    sc_list <- list(
      list(label="未治療", use_cipro=0, use_tmp=0, use_nit=0, evs=list(), color="#cc0000"),
      list(label="Cipro 14d", use_cipro=1, use_tmp=0, use_nit=0,
           evs=list(cipro=make_ev("cipro",350,12,14,0)), color="#0066cc"),
      list(label="TMP-SMX 14d", use_cipro=0, use_tmp=1, use_nit=0,
           evs=list(tmp=make_ev("tmp",152,12,14,0)), color="#009900"),
      list(label="Cipro+NIT ppx", use_cipro=1, use_tmp=0, use_nit=1,
           evs=list(
             cipro=make_ev("cipro",350,12,14,0),
             nit=make_ev("nit",75,24,180,14*24)
           ), color="#990099"),
      list(label="NIT 예방 단독", use_cipro=0, use_tmp=0, use_nit=1,
           evs=list(nit=make_ev("nit",75,24,180,0)), color="#ff9900")
    )
    lapply(sc_list, function(sc) {
      p <- c(base_p, list(use_cipro=sc$use_cipro, use_tmp=sc$use_tmp,
                          use_nit=sc$use_nit))
      df <- run_sim(p, sc$evs, end_h, delta = 6)
      df$Scenario <- sc$label
      df$color    <- sc$color
      df
    }) |> bind_rows()
  }, ignoreNULL = FALSE)

  # ======== VALUE BOXES =====================================================
  output$vb_gfr <- renderValueBox({
    d <- sim_data()
    gfr_last <- tail(d$GFR, 1)
    valueBox(round(gfr_last, 1), "최종 GFR (mL/min)", icon = icon("tint"),
             color = if (gfr_last >= 60) "green" else if (gfr_last >= 30) "yellow" else "red")
  })
  output$vb_stage <- renderValueBox({
    d <- sim_data()
    st <- tail(d$CKD_Stage, 1)
    valueBox(paste0("G", round(st)), "CKD 단계", icon = icon("kidney"),
             color = if (st <= 2) "green" else if (st <= 3) "yellow" else "red")
  })
  output$vb_scar <- renderValueBox({
    d <- sim_data()
    scar <- round(tail(d$ScarPct, 1), 1)
    valueBox(paste0(scar, "%"), "신장 반흔", icon = icon("ban"),
             color = if (scar < 10) "green" else if (scar < 30) "yellow" else "red")
  })
  output$vb_bact <- renderValueBox({
    d <- sim_data()
    b <- round(tail(d$Bacteria, 1), 1)
    valueBox(paste0("10^", b), "최종 세균 부하 (log₁₀)", icon = icon("bacteria"),
             color = if (b < 3) "green" else if (b < 5) "yellow" else "red")
  })

  # ======== PATIENT SUMMARY =================================================
  output$patient_summary_box <- renderUI({
    tags$div(
      tags$b("현재 환자 설정:"),
      tags$ul(
        tags$li(paste0("기저 GFR: ", input$GFR0, " mL/min")),
        tags$li(paste0("VUR 위험: ", input$VUR_risk)),
        tags$li(paste0("초기 세균: 10^", input$Bact_init, " CFU/g")),
        tags$li(paste0("내성 계수: ", input$Resist_F))
      )
    )
  })

  # ======== TAB 2: PK ======================================================
  output$plot_pk_cipro <- renderPlotly({
    d <- sim_data()
    peak_cipro <- max(d$Cipro_C, na.rm = TRUE)
    p <- plot_ly(d, x = ~time, y = ~Cipro_C, type = "scatter", mode = "lines",
                 line = list(color = "#0066cc", width = 2),
                 name = "Cipro Cp (μg/mL)") |>
      add_segments(x = 0, xend = max(d$time), y = 0.25, yend = 0.25,
                   line = list(color = "red", dash = "dot"),
                   name = "MIC = 0.25 μg/mL") |>
      layout(title = "Ciprofloxacin 혈중 농도 프로파일",
             xaxis = list(title = "시간 (h)"),
             yaxis = list(title = "Cp (μg/mL)"))
    p
  })

  output$plot_fAUC <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time, y = ~fAUC_MIC_OB, type = "scatter", mode = "lines",
            line = list(color = "#009900", width = 2),
            name = "fAUC/MIC") |>
      add_segments(x = 0, xend = max(d$time), y = 125, yend = 125,
                   line = list(color = "orange", dash = "dot"),
                   name = "Target fAUC/MIC = 125") |>
      layout(title = "fAUC/MIC 달성도",
             xaxis = list(title = "시간 (h)"),
             yaxis = list(title = "fAUC/MIC"))
  })

  output$plot_pk_tmp <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time, y = ~TMP_C, type = "scatter", mode = "lines",
            line = list(color = "#009900", width = 2)) |>
      layout(xaxis = list(title = "시간 (h)"),
             yaxis = list(title = "TMP Cp (μg/mL)"))
  })

  output$plot_pk_nit <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time, y = ~NIT_urine, type = "scatter", mode = "lines",
            line = list(color = "#ff9900", width = 2)) |>
      add_segments(x = 0, xend = max(d$time), y = 32, yend = 32,
                   line = list(color = "red", dash = "dot"),
                   name = "MIC = 32 μg/mL") |>
      layout(xaxis = list(title = "시간 (h)"),
             yaxis = list(title = "Nitrofurantoin 요중 농도 (μg/mL)"))
  })

  output$pk_param_table <- DT::renderDataTable({
    tbl <- data.frame(
      항생제      = c("Ciprofloxacin", "TMP-SMX", "Nitrofurantoin"),
      생체이용률   = c("70%", "95% / 85%", "75%"),
      Vc_L       = c(140, 85, "-"),
      CL_Lh      = c(25, "3.5 / 1.5", "-"),
      t_half     = c("4h", "10h", "0.3-1h"),
      MIC_target = c("≤0.25 μg/mL", "≤2 mg/L TMP", "≤32 μg/mL (urine)"),
      PK_PD_index = c("fAUC/MIC > 125", "T>MIC > 40%", "Cmax/MIC > 4"),
      stringsAsFactors = FALSE
    )
    DT::datatable(tbl, options = list(dom = "t", pageLength = 5),
                  rownames = FALSE)
  })

  # ======== TAB 3: Bacteria =================================================
  output$plot_bacteria <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time / 24, y = ~Bacteria, type = "scatter", mode = "lines",
            line = list(color = "#cc0000", width = 2)) |>
      add_segments(x = 0, xend = max(d$time / 24), y = 3, yend = 3,
                   line = list(color = "green", dash = "dot")) |>
      layout(xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "log₁₀ CFU/g", range = c(0, 9.5)))
  })

  output$plot_biofilm <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time / 24, y = ~Biofilm, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(200,100,100,0.2)",
            line = list(color = "#cc3300", width = 2)) |>
      layout(xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "생물막 분율 (0–1)", range = c(0, 1)))
  })

  output$plot_kill <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time / 24, y = ~Kill_drug, type = "scatter", mode = "lines",
            line = list(color = "#007bff", width = 2)) |>
      layout(xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "항균 Kill Rate (log₁₀ CFU/h)"))
  })

  output$plot_immune <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time / 24, y = ~Neutrophil, type = "scatter", mode = "lines",
            name = "Neutrophil", line = list(color = "#ffd700", width = 2)) |>
      add_trace(y = ~Macrophage, name = "Macrophage",
                line = list(color = "#ffa500", width = 2)) |>
      layout(xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "정규화 세포 수 (1=기저)"))
  })

  # ======== TAB 4: Renal ====================================================
  output$plot_gfr <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time / 24, y = ~GFR, type = "scatter", mode = "lines",
            line = list(color = "#2255aa", width = 2.5)) |>
      add_segments(x = 0, xend = max(d$time / 24), y = 60, yend = 60,
                   line = list(color = "#ffa500", dash = "dot"), name = "CKD G3a") |>
      add_segments(x = 0, xend = max(d$time / 24), y = 30, yend = 30,
                   line = list(color = "#ff4400", dash = "dot"), name = "CKD G4") |>
      layout(xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "GFR (mL/min/1.73m²)", range = c(0, 115)))
  })

  output$plot_creat <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time / 24, y = ~Creatinine, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(204,0,0,0.1)",
            line = list(color = "#cc0000", width = 2)) |>
      layout(xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "혈청 크레아티닌 (mg/dL)"))
  })

  output$plot_scar <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time / 24, y = ~ScarPct, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(128,0,64,0.15)",
            line = list(color = "#880040", width = 2)) |>
      layout(xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "신장 반흔 (%)"))
  })

  output$plot_ckd <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time / 24, y = ~CKD_Stage, type = "scatter", mode = "lines",
            line = list(color = "#113388", width = 2, shape = "hv")) |>
      layout(xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "CKD 단계", range = c(0.5, 5.5),
                          tickvals = 1:5,
                          ticktext = c("G1 (≥90)", "G2 (60–89)",
                                       "G3 (30–59)", "G4 (15–29)", "G5 (<15)")))
  })

  # ======== TAB 5: Scenario =================================================
  output$plot_scenario_gfr <- renderPlotly({
    d <- scenario_data()
    plot_ly(d, x = ~time / 24, y = ~GFR, color = ~Scenario,
            type = "scatter", mode = "lines",
            colors = c("#cc0000","#0066cc","#009900","#990099","#ff9900")) |>
      layout(xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "GFR (mL/min/1.73m²)"),
             legend = list(orientation = "h"))
  })

  output$plot_scenario_bact <- renderPlotly({
    d <- scenario_data()
    plot_ly(d, x = ~time / 24, y = ~Bacteria, color = ~Scenario,
            type = "scatter", mode = "lines",
            colors = c("#cc0000","#0066cc","#009900","#990099","#ff9900")) |>
      layout(xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "log₁₀ CFU/g"))
  })

  output$plot_scenario_scar <- renderPlotly({
    d <- scenario_data()
    plot_ly(d, x = ~time / 24, y = ~ScarPct, color = ~Scenario,
            type = "scatter", mode = "lines",
            colors = c("#cc0000","#0066cc","#009900","#990099","#ff9900")) |>
      layout(xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "신장 반흔 (%)"))
  })

  output$scenario_table <- DT::renderDataTable({
    d <- scenario_data()
    d |>
      group_by(Scenario) |>
      summarise(
        `최종 GFR`       = round(tail(GFR, 1), 1),
        `GFR 변화`       = round(tail(GFR, 1) - head(GFR, 1), 1),
        `최종 CKD 단계`  = paste0("G", round(tail(CKD_Stage, 1))),
        `최종 반흔 (%)`  = round(tail(ScarPct, 1), 1),
        `세균 최저치`    = round(min(Bacteria), 2),
        `최종 세균 부하` = round(tail(Bacteria, 1), 2)
      ) |>
      DT::datatable(options = list(pageLength = 10), rownames = FALSE)
  })

  # ======== TAB 6: Biomarker ================================================
  output$plot_biomarker <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time / 24, y = ~IL6, type = "scatter", mode = "lines",
            name = "IL-6 (normalised)", line = list(color = "#44bb44", width = 2)) |>
      add_trace(y = ~TGFb1, name = "TGF-β1",
                line = list(color = "#228822", width = 2, dash = "dash")) |>
      add_trace(y = ~Collagen, name = "Collagen",
                line = list(color = "#880040", width = 2, dash = "dot")) |>
      layout(xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "정규화 수준 (1 = 기저)"),
             legend = list(orientation = "h"))
  })

  output$biomarker_legend <- renderTable({
    data.frame(
      마커      = c("IL-6", "TGF-β1", "Collagen", "Neutrophil", "Macrophage"),
      임상적의미 = c(
        "급성기 단백 유도, 발열 신호",
        "섬유화 핵심 매개체",
        "간질 섬유화 진행도",
        "급성 세균 감염 지표",
        "만성 염증 / M2→섬유화"
      ),
      정상범위  = c("< 5 pg/mL", "< 10 pg/mL", "기저=1", "기저=1", "기저=1")
    )
  })

  output$plot_sepsis_risk <- renderPlotly({
    d <- sim_data()
    d$SepsisRisk <- as.integer(d$Bacteria > 7.5)
    plot_ly(d, x = ~time / 24, y = ~SepsisRisk, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(255,0,0,0.2)",
            line = list(color = "#cc0000", width = 2)) |>
      layout(xaxis = list(title = "시간 (일)"),
             yaxis = list(title = "요로패혈증 위험 (1=High, 0=Low)",
                          range = c(-0.1, 1.2)))
  })

  output$plot_collagen_scar <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~Collagen, y = ~ScarPct, type = "scatter", mode = "markers",
            marker = list(color = ~time / max(d$time),
                          colorscale = "Viridis", size = 3,
                          colorbar = list(title = "시간 비율"))) |>
      layout(xaxis = list(title = "콜라겐 침착 (정규화)"),
             yaxis = list(title = "신장 반흔 (%)"))
  })

  # ======== TAB 7: References ===============================================
  output$ref_html <- renderUI({
    HTML("
    <div style='font-size:13px; line-height:1.8; padding:10px;'>
    <h4>핵심 병태생리 참고문헌</h4>
    <ol>
    <li>Nicolle LE (2008). Uncomplicated urinary tract infection in adults including uncomplicated pyelonephritis.
        <a href='https://pubmed.ncbi.nlm.nih.gov/19090895/' target='_blank'>Urol Clin North Am 35:1–12</a></li>
    <li>Hooton TM (2012). Uncomplicated urinary tract infection.
        <a href='https://pubmed.ncbi.nlm.nih.gov/22591447/' target='_blank'>N Engl J Med 366:1028–1037</a></li>
    <li>Foxman B (2010). The epidemiology of urinary tract infection.
        <a href='https://pubmed.ncbi.nlm.nih.gov/20141568/' target='_blank'>Nat Rev Urol 7:653–660</a></li>
    <li>Scholes D et al. (2000). Risk factors for recurrent urinary tract infection in young women.
        <a href='https://pubmed.ncbi.nlm.nih.gov/11011900/' target='_blank'>J Infect Dis 182:1177–1182</a></li>
    <li>Schappert SM, Rechtsteiner EA (2011). Ambulatory medical care utilization estimates.
        <a href='https://pubmed.ncbi.nlm.nih.gov/21976217/' target='_blank'>Vital Health Stat 13:1–32</a></li>
    </ol>
    <h4>VUR 및 신장 반흔</h4>
    <ol start='6'>
    <li>Weiss R et al. (1992). Clinical significance of primary vesicoureteral reflux.
        <a href='https://pubmed.ncbi.nlm.nih.gov/1580441/' target='_blank'>J Urol 148:1682–1685</a></li>
    <li>Craig JC et al. (2000). Vesico-ureteric reflux and timing of micturating cystourethrography.
        <a href='https://pubmed.ncbi.nlm.nih.gov/10768834/' target='_blank'>Lancet 356:1160–1161</a></li>
    <li>Salo J et al. (2011). Childhood urinary tract infections as a cause of CKD.
        <a href='https://pubmed.ncbi.nlm.nih.gov/21606568/' target='_blank'>Pediatrics 128:840–847</a></li>
    </ol>
    <h4>항생제 PK/PD</h4>
    <ol start='9'>
    <li>Craig WA (1998). PK/PD parameters: rationale for antibacterial dosing.
        <a href='https://pubmed.ncbi.nlm.nih.gov/9647345/' target='_blank'>Clin Infect Dis 26:1–12</a></li>
    <li>Rybak MJ (2006). PD: relation to antimicrobial resistance.
        <a href='https://pubmed.ncbi.nlm.nih.gov/16819721/' target='_blank'>Am J Infect Control 34:S38–45</a></li>
    <li>Turnidge J, Paterson DL (2007). Setting and revising antibacterial susceptibility breakpoints.
        <a href='https://pubmed.ncbi.nlm.nih.gov/17638701/' target='_blank'>Clin Microbiol Rev 20:391–408</a></li>
    <li>Gupta K et al. (2011). International clinical practice guidelines for the treatment of acute uncomplicated cystitis and pyelonephritis.
        <a href='https://pubmed.ncbi.nlm.nih.gov/21292654/' target='_blank'>Clin Infect Dis 52:e103–e120</a></li>
    </ol>
    <h4>신장 섬유화 / CKD 진행</h4>
    <ol start='13'>
    <li>Eddy AA (2014). Overview of the cellular and molecular basis of kidney fibrosis.
        <a href='https://pubmed.ncbi.nlm.nih.gov/24523596/' target='_blank'>Kidney Int Suppl 4:2–8</a></li>
    <li>Liu Y (2011). Cellular and molecular mechanisms of renal fibrosis.
        <a href='https://pubmed.ncbi.nlm.nih.gov/21655640/' target='_blank'>Nat Rev Nephrol 7:684–696</a></li>
    <li>Zeisberg M, Neilson EG (2010). Mechanisms of tubulointerstitial fibrosis.
        <a href='https://pubmed.ncbi.nlm.nih.gov/20360843/' target='_blank'>J Am Soc Nephrol 21:1819–1834</a></li>
    </ol>
    </div>
    ")
  })
}

shinyApp(ui, server)
