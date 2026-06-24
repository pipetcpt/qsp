## =============================================================================
## Ovarian Cancer (HGSOC) QSP Shiny Dashboard
## Interactive Simulation: Carboplatin/Paclitaxel, PARP Inhibitors, Bevacizumab
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ── mrgsolve model code (embedded) ─────────────────────────────────────────
oc_code <- '
$PARAM
CL_CAR=4.2, V1_CAR=15.0, Q_CAR=6.5, V2_CAR=35.0,
CL_PAC=13.2, V1_PAC=6.5, Q2_PAC=7.0, V2_PAC=113.0, Q3_PAC=2.2, V3_PAC=1088.0,
ka_OLA=1.74, F_OLA=0.77, CL_OLA=8.9, V1_OLA=67.0, Q_OLA=4.2, V2_OLA=100.0,
ka_NIRA=0.36, F_NIRA=0.73, CL_NIRA=16.2, V1_NIRA=537.0, Q_NIRA=5.0, V2_NIRA=537.0,
CL_BEV=0.207, V1_BEV=2.91, Q_BEV=0.469, V2_BEV=1.91,
VEGF0=0.15, ksyn_VEGF=1.5, kdeg_VEGF=10.0, kbind_BEV=50.0,
k_adduct=0.015, k_repair=0.08, HRD_sens=1.0,
k_HRD_in=0.1, k_HRD_out=0.02,
TV0=50.0, kg=0.008, TV_max=3000.0,
k_kill_Pt=0.004, k_kill_T=0.001, k_kill_Pi=0.003,
CA125_0=300.0, ksyn_CA125=3.0, kdeg_CA125=0.03,
CD8T_0=1.0, k_CD8_in=0.1, k_CD8_out=0.1, k_exhaust=0.3, k_ICI=2.0,
ICI_flag=0, BRCAmut=1, HRD_pos=1

$CMT CAR_C1 CAR_C2 PAC_C1 PAC_C2 PAC_C3 OLA_gut OLA_C1 OLA_C2
     NIRA_C1 NIRA_C2 BEV_C1 BEV_C2 VEGF TV CA125 Pt_DNA CD8T HRD

$MAIN
double eff_HRD = HRD_pos * HRD_sens + (1 - HRD_pos) * 0.3;

$ODE
dxdt_CAR_C1 = -(CL_CAR/V1_CAR)*CAR_C1 - (Q_CAR/V1_CAR)*CAR_C1 + (Q_CAR/V2_CAR)*CAR_C2;
dxdt_CAR_C2 = (Q_CAR/V1_CAR)*CAR_C1 - (Q_CAR/V2_CAR)*CAR_C2;
dxdt_PAC_C1 = -(CL_PAC/V1_PAC)*PAC_C1 - (Q2_PAC/V1_PAC)*PAC_C1 + (Q2_PAC/V2_PAC)*PAC_C2
              - (Q3_PAC/V1_PAC)*PAC_C1 + (Q3_PAC/V3_PAC)*PAC_C3;
dxdt_PAC_C2 = (Q2_PAC/V1_PAC)*PAC_C1 - (Q2_PAC/V2_PAC)*PAC_C2;
dxdt_PAC_C3 = (Q3_PAC/V1_PAC)*PAC_C1 - (Q3_PAC/V3_PAC)*PAC_C3;
dxdt_OLA_gut = -ka_OLA * OLA_gut;
dxdt_OLA_C1  = (F_OLA * ka_OLA * OLA_gut) / V1_OLA - (CL_OLA/V1_OLA)*OLA_C1
               - (Q_OLA/V1_OLA)*OLA_C1 + (Q_OLA/V2_OLA)*OLA_C2;
dxdt_OLA_C2  = (Q_OLA/V1_OLA)*OLA_C1 - (Q_OLA/V2_OLA)*OLA_C2;
dxdt_NIRA_C1 = (F_NIRA * ka_NIRA * NIRA_C1) / V1_NIRA - (CL_NIRA/V1_NIRA)*NIRA_C1
               - (Q_NIRA/V1_NIRA)*NIRA_C1 + (Q_NIRA/V2_NIRA)*NIRA_C2;
dxdt_NIRA_C2 = (Q_NIRA/V1_NIRA)*NIRA_C1 - (Q_NIRA/V2_NIRA)*NIRA_C2;
double BEV_effect_v = kbind_BEV * BEV_C1 * VEGF;
dxdt_BEV_C1 = -(CL_BEV/V1_BEV)*BEV_C1 - (Q_BEV/V1_BEV)*BEV_C1 + (Q_BEV/V2_BEV)*BEV_C2;
dxdt_BEV_C2 = (Q_BEV/V1_BEV)*BEV_C1 - (Q_BEV/V2_BEV)*BEV_C2;
dxdt_VEGF = ksyn_VEGF*(TV/TV0) - kdeg_VEGF*VEGF - BEV_effect_v;
dxdt_Pt_DNA = k_adduct*CAR_C1 - k_repair*(1 - 0.6*eff_HRD)*Pt_DNA;
double parp_trap = 0.0;
if(OLA_C1 > 0.1)  parp_trap += 0.8 * OLA_C1 / (OLA_C1 + 500);
if(NIRA_C1 > 0.1) parp_trap += 0.7 * NIRA_C1 / (NIRA_C1 + 2000);
parp_trap = parp_trap * eff_HRD;
dxdt_HRD = k_HRD_in*parp_trap - k_HRD_out*HRD;
double ICI_eff = 1.0 + ICI_flag*(k_ICI - 1.0);
double exhaust = k_exhaust * TV / TV_max;
dxdt_CD8T = k_CD8_in*ICI_eff - k_CD8_out*CD8T - exhaust*CD8T;
double grow_t   = kg * TV * log(TV_max / TV);
double kill_Pt  = k_kill_Pt * Pt_DNA * TV;
double pac_eff  = PAC_C1 / (PAC_C1 + 100.0);
double kill_pac = 0.6 * k_kill_Pt * pac_eff * TV;
double parp_kill= k_kill_Pi * HRD * TV;
double kill_CD8 = k_kill_T * CD8T * TV;
dxdt_TV = grow_t - kill_Pt - kill_pac - parp_kill - kill_CD8;
dxdt_CA125 = ksyn_CA125 * TV - kdeg_CA125 * CA125;

$TABLE
capture CAR_Conc = CAR_C1;
capture PAC_Conc = PAC_C1;
capture OLA_Conc = OLA_C1;
capture NIRA_Conc = NIRA_C1;
capture BEV_Conc = BEV_C1;
capture VEGF_f = VEGF;
capture TumorVol = TV;
capture CA125_lvl = CA125;
capture PtDNA_r = Pt_DNA;
capture HRD_d = HRD;
capture CD8T_r = CD8T;

$INIT
CAR_C1=0,CAR_C2=0,PAC_C1=0,PAC_C2=0,PAC_C3=0,
OLA_gut=0,OLA_C1=0,OLA_C2=0,NIRA_C1=0,NIRA_C2=0,
BEV_C1=0,BEV_C2=0,VEGF=0.15,TV=50,CA125=300,
Pt_DNA=0,CD8T=1.0,HRD=0
'

## Compile once at startup
mod_base <- mcode("oc_shiny", oc_code, quiet=TRUE)

## ── Shiny UI ───────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "난소암 QSP 대시보드 | Ovarian Cancer QSP"),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("① 환자 프로파일", tabName="patient",  icon=icon("user")),
      menuItem("② 약물 PK",       tabName="pk",       icon=icon("pills")),
      menuItem("③ PD 바이오마커", tabName="pd",       icon=icon("flask")),
      menuItem("④ 종양 반응",     tabName="tumor",    icon=icon("chart-line")),
      menuItem("⑤ 시나리오 비교", tabName="scenario", icon=icon("layer-group")),
      menuItem("⑥ 바이오마커 패널", tabName="biomarker", icon=icon("microscope"))
    ),
    hr(),
    h5("── 환자 특성 ──", style="color:#ccc;padding-left:15px"),
    sliderInput("TV0",    "초기 종양 크기 (cm³)",  10, 500, 50, step=10),
    sliderInput("CA125_0","초기 CA-125 (U/mL)",    35, 2000, 300, step=10),
    selectInput("BRCAmut","BRCA 변이 상태",
      choices=c("BRCA1/2 변이 (mut)"=1, "BRCA 야생형 (wt)"=0), selected=1),
    selectInput("HRD_pos","HRD 상태",
      choices=c("HRD 양성 (≥42)"=1, "HRD 음성 (<42)"=0), selected=1),
    sliderInput("GFR",    "GFR (mL/min)",    30, 150, 90, step=5),
    hr(),
    h5("── 치료 설정 ──", style="color:#ccc;padding-left:15px"),
    sliderInput("n_cycles","화학요법 사이클 수", 3, 8, 6, step=1),
    selectInput("maint_drug","유지요법 약물",
      choices=c("없음"="none","오라파립(BRCA+)"="ola","니라파립(HRD+)"="nira",
                "베바시주맙"="bev","오라파립+베바시주맙(PAOLA-1)"="ola_bev"), selected="ola"),
    checkboxInput("add_bev_front","1차 화학요법+베바시주맙", value=FALSE),
    checkboxInput("ICI_flag","면역관문억제제 (항PD-1) 추가", value=FALSE),
    sliderInput("sim_days","시뮬레이션 기간 (일)", 180, 1095, 730, step=30),
    actionButton("run_sim","▶ 시뮬레이션 실행", class="btn-success btn-block")
  ),

  dashboardBody(
    tabItems(
      ## Tab 1: Patient Profile
      tabItem("patient",
        fluidRow(
          valueBoxOutput("box_stage"), valueBoxOutput("box_brca"),
          valueBoxOutput("box_hrd")
        ),
        fluidRow(
          box(title="환자 특성 요약", width=6, status="danger",
            tableOutput("patient_summary")
          ),
          box(title="BRCA/HRD 치료 적합성", width=6, status="warning",
            tableOutput("eligibility_table")
          )
        ),
        fluidRow(
          box(title="HGSOC 병태생리 개요", width=12, status="info",
            HTML('
<table class="table table-bordered table-sm">
<thead><tr style="background:#C0392B;color:white">
  <th>특성</th><th>HGSOC (고등급 장액성)</th><th>LGSOC (저등급)</th>
</tr></thead>
<tbody>
<tr><td>빈도</td><td>~70% of all OC</td><td>~5-10%</td></tr>
<tr><td>주요 변이</td><td>TP53 (>96%), HRD (50%)</td><td>KRAS/BRAF/RAS pathway</td></tr>
<tr><td>BRCA1/2 변이</td><td>~15% germline, ~7% somatic</td><td>드묾</td></tr>
<tr><td>HRD 양성률</td><td>~50% (BRCA+포함)</td><td>~10%</td></tr>
<tr><td>백금 민감도</td><td>초기 민감; 내성 흔함</td><td>백금 불량 반응</td></tr>
<tr><td>PARP억제제 효과</td><td>HRD+에서 탁월 (SOLO-1, PRIMA)</td><td>낮음</td></tr>
<tr><td>5년 생존율</td><td>FIGO III: ~35%, IV: ~20%</td><td>III: ~50%</td></tr>
</tbody>
</table>'
            )
          )
        )
      ),

      ## Tab 2: Drug PK
      tabItem("pk",
        fluidRow(
          box(title="카보플라틴 PK (중심 구획)", width=6, status="warning",
            plotlyOutput("plot_carbo_pk", height=280)
          ),
          box(title="파클리탁셀 PK (중심 구획)", width=6, status="warning",
            plotlyOutput("plot_pacli_pk", height=280)
          )
        ),
        fluidRow(
          box(title="오라파립 PK (중심 구획)", width=6, status="danger",
            plotlyOutput("plot_ola_pk", height=280)
          ),
          box(title="니라파립 PK", width=6, status="danger",
            plotlyOutput("plot_nira_pk", height=280)
          )
        ),
        fluidRow(
          box(title="베바시주맙 PK + 유리 VEGF-A", width=12, status="primary",
            plotlyOutput("plot_bev_vegf", height=300)
          )
        )
      ),

      ## Tab 3: PD Biomarkers
      tabItem("pd",
        fluidRow(
          box(title="CA-125 혈청 농도 (U/mL)", width=6, status="danger",
            plotlyOutput("plot_ca125", height=300)
          ),
          box(title="백금-DNA 부가물 (상대값)", width=6, status="warning",
            plotlyOutput("plot_ptdna", height=300)
          )
        ),
        fluidRow(
          box(title="HRD 손상 축적 (PARP 억제제)", width=6, status="danger",
            plotlyOutput("plot_hrd", height=300)
          ),
          box(title="CD8+ T세포 (종양내 침윤 CTL)", width=6, status="success",
            plotlyOutput("plot_cd8t", height=300)
          )
        ),
        fluidRow(
          box(title="바이오마커 요약 (최종 시점)", width=12, status="info",
            tableOutput("bm_summary")
          )
        )
      ),

      ## Tab 4: Tumor Response
      tabItem("tumor",
        fluidRow(
          valueBoxOutput("box_best_resp"), valueBoxOutput("box_ca125_nadir"),
          valueBoxOutput("box_pfs")
        ),
        fluidRow(
          box(title="종양 부피 시간 경과 (cm³)", width=8, status="danger",
            plotlyOutput("plot_tv", height=380)
          ),
          box(title="RECIST 반응 분류", width=4, status="warning",
            tableOutput("recist_table")
          )
        ),
        fluidRow(
          box(title="종양 부피 변화율 (베이스라인 대비 %)", width=12, status="info",
            plotlyOutput("plot_tv_pct", height=250)
          )
        )
      ),

      ## Tab 5: Scenario Comparison
      tabItem("scenario",
        fluidRow(
          box(title="6가지 치료 시나리오 — 종양 부피 비교", width=12, status="danger",
            plotlyOutput("plot_scenario_tv", height=380)
          )
        ),
        fluidRow(
          box(title="6가지 치료 시나리오 — CA-125 비교", width=12, status="warning",
            plotlyOutput("plot_scenario_ca125", height=300)
          )
        ),
        fluidRow(
          box(title="시나리오 요약 테이블", width=12, status="info",
            DTOutput("scenario_table")
          )
        )
      ),

      ## Tab 6: Biomarker Panel
      tabItem("biomarker",
        fluidRow(
          box(title="종합 바이오마커 패널 (시간-농도)", width=12, status="primary",
            plotlyOutput("plot_bm_panel", height=500)
          )
        ),
        fluidRow(
          box(title="BRCA/HRD 진단 트리 — 치료 선택", width=6, status="danger",
            HTML('
<div style="font-family:monospace;line-height:1.8;font-size:13px">
<b>난소암 1차 치료 알고리즘 (BRCA/HRD)</b><br>
└─ <b>수술 후 잔류병변?</b><br>
   ├─ <b>R0 (no residue)</b><br>
   │  ├─ Carbo+Pacli ×6<br>
   │  └─ BRCA mut? → Olaparib 2년 (SOLO-1: mPFS NR)<br>
   └─ <b>R1/R2 또는 NACT</b><br>
      ├─ Carbo+Pacli+Bev ×6 → Bev 유지<br>
      └─ HRD+? (BRCA포함)<br>
         ├─ Ola+Bev 유지 (PAOLA-1: mPFS 22.1mo)<br>
         └─ Niraparib 유지 (PRIMA: mPFS 13.8mo)
</div>'
            )
          ),
          box(title="주요 임상시험 참조 수치", width=6, status="warning",
            tableOutput("trial_ref_table")
          )
        )
      )
    )
  )
)

## ── Server ─────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## Reactive: build dosing events and run simulation
  sim_result <- eventReactive(input$run_sim, {
    TV0_v   <- input$TV0
    CA125_0v <- input$CA125_0
    BRCA_v  <- as.numeric(input$BRCAmut)
    HRD_v   <- as.numeric(input$HRD_pos)
    GFR_v   <- input$GFR
    nc      <- input$n_cycles
    dur     <- input$sim_days
    ici_f   <- if(input$ICI_flag) 1 else 0

    ## Carboplatin dose (Calvert AUC=6)
    carbo_mg <- 6 * (GFR_v + 25)
    carbo_times <- seq(0, (nc-1)*21, by=21)

    ## Paclitaxel 300mg
    pacli_times <- carbo_times

    ## Maintenance start day
    maint_start <- nc * 21

    ## Build event table
    ev_list <- list(
      ## Carboplatin IV
      data.frame(time=carbo_times, cmt="CAR_C1", amt=carbo_mg/15, evid=1),
      ## Paclitaxel IV
      data.frame(time=pacli_times, cmt="PAC_C1", amt=300000/6.5, evid=1)
    )

    ## Bevacizumab in front-line
    if(input$add_bev_front) {
      bev_times <- carbo_times
      ev_list <- c(ev_list, list(
        data.frame(time=bev_times, cmt="BEV_C1", amt=1050/2.91, evid=1)
      ))
    }

    ## Maintenance drugs
    maint <- input$maint_drug
    if(maint == "ola" || maint == "ola_bev") {
      ola_times <- seq(maint_start*24, (dur)*24-12, by=12)
      ev_list <- c(ev_list, list(
        data.frame(time=ola_times, cmt="OLA_gut", amt=300, evid=1)
      ))
    }
    if(maint == "nira") {
      nira_times <- seq(maint_start*24, (dur)*24-24, by=24)
      ev_list <- c(ev_list, list(
        data.frame(time=nira_times, cmt="NIRA_C1", amt=300, evid=1)
      ))
    }
    if(maint == "bev" || maint == "ola_bev") {
      bev_m_times <- seq(maint_start, dur, by=21)
      ev_list <- c(ev_list, list(
        data.frame(time=bev_m_times, cmt="BEV_C1", amt=1050/2.91, evid=1)
      ))
    }

    ev_df <- bind_rows(ev_list) %>% arrange(time)

    mod_run <- mod_base %>%
      param(TV0=TV0_v, CA125_0=CA125_0v,
            BRCAmut=BRCA_v, HRD_pos=HRD_v,
            ICI_flag=ici_f)

    out <- mrgsim_df(mod_run, data=as_data_set(as.ev(ev_df)),
                     end=dur, delta=1)
    out
  }, ignoreNULL=FALSE)

  ## Six-scenario simulation (fixed for comparison tab)
  all_scenarios <- reactive({
    scenarios <- list(
      list(name="S1: 무치료", brca=1, hrd=1, maint="none", bev_front=FALSE, ici=0),
      list(name="S2: Carbo+Pacli ×6", brca=1, hrd=1, maint="none", bev_front=FALSE, ici=0),
      list(name="S3: Carbo+Pacli+Bev→Bev 유지", brca=0, hrd=1, maint="bev", bev_front=TRUE, ici=0),
      list(name="S4: →오라파립 유지 (BRCA+)", brca=1, hrd=1, maint="ola", bev_front=FALSE, ici=0),
      list(name="S5: →니라파립 유지 (HRD+)", brca=0, hrd=1, maint="nira", bev_front=FALSE, ici=0),
      list(name="S6: Carbo+Pacli+Bev→Ola+Bev (PAOLA-1)", brca=0, hrd=1, maint="ola_bev", bev_front=TRUE, ici=0)
    )

    lapply(scenarios, function(s) {
      ev_list <- list(
        data.frame(time=seq(0,5*21,by=21), cmt="CAR_C1", amt=690/15, evid=1),
        data.frame(time=seq(0,5*21,by=21), cmt="PAC_C1", amt=300000/6.5, evid=1)
      )
      maint_start <- 126
      if(s$bev_front)
        ev_list <- c(ev_list, list(data.frame(time=seq(0,5*21,by=21), cmt="BEV_C1", amt=1050/2.91, evid=1)))
      if(s$maint %in% c("ola","ola_bev")) {
        ev_list <- c(ev_list, list(data.frame(time=seq(maint_start*24, 730*24-12, by=12), cmt="OLA_gut", amt=300, evid=1)))
      }
      if(s$maint == "nira") {
        ev_list <- c(ev_list, list(data.frame(time=seq(maint_start*24, 730*24-24, by=24), cmt="NIRA_C1", amt=300, evid=1)))
      }
      if(s$maint %in% c("bev","ola_bev")) {
        ev_list <- c(ev_list, list(data.frame(time=seq(maint_start, 730, by=21), cmt="BEV_C1", amt=1050/2.91, evid=1)))
      }
      ev_df <- bind_rows(ev_list) %>% arrange(time)
      mod_run <- mod_base %>% param(BRCAmut=s$brca, HRD_pos=s$hrd, ICI_flag=s$ici)
      out <- mrgsim_df(mod_run, data=as_data_set(as.ev(ev_df)), end=730, delta=1)
      out$Scenario <- s$name
      out
    }) %>% bind_rows()
  })

  ## ── ValueBoxes ───────────────────────────────────────────────
  output$box_stage <- renderValueBox({
    valueBox("FIGO III-IV", "병기 (대부분 진단 시)", icon=icon("hospital"), color="red")
  })
  output$box_brca <- renderValueBox({
    label <- if(input$BRCAmut=="1") "BRCA 변이" else "BRCA 야생형"
    valueBox(label, "BRCA 상태", icon=icon("dna"), color="purple")
  })
  output$box_hrd <- renderValueBox({
    label <- if(input$HRD_pos=="1") "HRD 양성 (≥42)" else "HRD 음성 (<42)"
    valueBox(label, "HRD 점수", icon=icon("vial"), color="orange")
  })
  output$box_best_resp <- renderValueBox({
    df <- sim_result()
    best <- min(df$TumorVol, na.rm=TRUE)
    tv0  <- input$TV0
    pct  <- round((best - tv0)/tv0*100, 1)
    valueBox(paste0(pct, "%"), "최선 종양 반응 (RECIST)", icon=icon("chart-line"), color="red")
  })
  output$box_ca125_nadir <- renderValueBox({
    df <- sim_result()
    nadir <- round(min(df$CA125_lvl, na.rm=TRUE))
    valueBox(paste0(nadir, " U/mL"), "CA-125 최저치", icon=icon("flask"), color="purple")
  })
  output$box_pfs <- renderValueBox({
    df <- sim_result()
    tv_thresh <- input$TV0 * 2
    pfs_d <- df %>% filter(TumorVol > tv_thresh) %>% pull(time) %>% min()
    label <- if(is.infinite(pfs_d) || is.na(pfs_d))
      paste0(">", input$sim_days, " 일")
    else
      paste0(round(pfs_d), " 일 (", round(pfs_d/30.4, 1), " mo)")
    valueBox(label, "PFS 추정 (종양 2배 시점)", icon=icon("calendar"), color="green")
  })

  ## ── Tab 1: Patient Summary ────────────────────────────────────
  output$patient_summary <- renderTable({
    data.frame(
      항목 = c("초기 종양 크기","초기 CA-125","BRCA 변이","HRD 상태","GFR","화학요법 사이클"),
      값   = c(paste(input$TV0,"cm³"), paste(input$CA125_0,"U/mL"),
               ifelse(input$BRCAmut=="1","양성 (돌연변이)","음성 (야생형)"),
               ifelse(input$HRD_pos=="1","양성 (≥42)","음성 (<42)"),
               paste(input$GFR,"mL/min"),
               paste(input$n_cycles,"사이클"))
    )
  })

  output$eligibility_table <- renderTable({
    brca <- input$BRCAmut == "1"
    hrd  <- input$HRD_pos == "1"
    data.frame(
      치료 = c("오라파립 유지 (SOLO-1)", "니라파립 유지 (PRIMA)", "오라파립+베바시주맙 (PAOLA-1)", "베바시주맙 (ICON7/AVA)", "임상시험 (ICI)"),
      적합성 = c(
        ifelse(brca, "✅ 적합 (BRCA+)", "❌ 제한적"),
        ifelse(hrd,  "✅ 적합 (HRD+)", "⚠️ 제한적 (HRD-도 일부 이득)"),
        ifelse(hrd,  "✅ 최적 (HRD+)", "⚠️ 이득 제한"),
        "✅ 모든 환자 (고위험군)", "⚠️ 연구 중"
      )
    )
  })

  ## ── Tab 2: PK Plots ────────────────────────────────────────────
  pk_plot <- function(df, col, title, ylab, xlim=NULL) {
    p <- ggplot(df, aes(x=time, y=.data[[col]])) +
      geom_line(color="#C0392B") +
      labs(title=title, x="일(Day)", y=ylab) +
      theme_bw(base_size=10)
    if(!is.null(xlim)) p <- p + coord_cartesian(xlim=xlim)
    ggplotly(p)
  }

  output$plot_carbo_pk <- renderPlotly({
    df <- sim_result()
    pk_plot(df, "CAR_Conc", "Carboplatin PK", "농도 (µg/mL)", xlim=c(0,250))
  })
  output$plot_pacli_pk <- renderPlotly({
    df <- sim_result()
    pk_plot(df, "PAC_Conc", "Paclitaxel PK", "농도 (ng/mL)", xlim=c(0,250))
  })
  output$plot_ola_pk <- renderPlotly({
    df <- sim_result()
    pk_plot(df, "OLA_Conc", "오라파립 PK (유지요법 기간)", "농도 (ng/mL)")
  })
  output$plot_nira_pk <- renderPlotly({
    df <- sim_result()
    pk_plot(df, "NIRA_Conc", "니라파립 PK", "농도 (ng/mL)")
  })
  output$plot_bev_vegf <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df) +
      geom_line(aes(x=time, y=BEV_Conc*10, color="베바시주맙 (×10 mg/L)")) +
      geom_line(aes(x=time, y=VEGF_f*1000, color="유리 VEGF-A (×1000 ng/mL)")) +
      labs(title="베바시주맙 PK & VEGF 억제", x="일(Day)", y="농도") +
      scale_color_manual(values=c("베바시주맙 (×10 mg/L)"="#1A237E","유리 VEGF-A (×1000 ng/mL)"="#FF8A65")) +
      theme_bw(base_size=10) + theme(legend.position="bottom")
    ggplotly(p)
  })

  ## ── Tab 3: PD Biomarkers ───────────────────────────────────────
  output$plot_ca125 <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x=time, y=CA125_lvl)) +
      geom_line(color="#E91E63") +
      geom_hline(yintercept=35, linetype="dashed", color="darkgreen") +
      annotate("text", x=max(df$time)*0.8, y=50, label="정상 상한 35 U/mL", size=3, color="darkgreen") +
      scale_y_log10() +
      labs(title="CA-125 혈청 농도", x="일(Day)", y="CA-125 (U/mL, log scale)") +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$plot_ptdna <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x=time, y=PtDNA_r)) +
      geom_line(color="#FF6F00") +
      labs(title="백금-DNA 부가물 (상대값)", x="일(Day)", y="부가물 (0–1)") +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$plot_hrd <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x=time, y=HRD_d)) +
      geom_line(color="#880E4F") +
      labs(title="HRD 손상 축적 (PARPi 트래핑)", x="일(Day)", y="HRD 손상 (0–1)") +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$plot_cd8t <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(x=time, y=CD8T_r)) +
      geom_line(color="#2E7D32") +
      geom_hline(yintercept=1.0, linetype="dashed", color="gray50") +
      labs(title="CD8+ T세포 (종양내 CTL)", x="일(Day)", y="상대 수준 (기저치=1)") +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$bm_summary <- renderTable({
    df <- sim_result()
    last <- tail(df, 1)
    data.frame(
      바이오마커 = c("CA-125 (U/mL)","백금-DNA 부가물 (상대)","HRD 손상 (0-1)","CD8+ T세포 (상대)","유리 VEGF-A (ng/mL)","종양 부피 (cm³)"),
      최종값 = c(round(last$CA125_lvl,1), round(last$PtDNA_r,3),
                  round(last$HRD_d,3), round(last$CD8T_r,3),
                  round(last$VEGF_f,3), round(last$TumorVol,1))
    )
  })

  ## ── Tab 4: Tumor Response ─────────────────────────────────────
  output$plot_tv <- renderPlotly({
    df <- sim_result()
    thresh <- input$TV0 * 2
    p <- ggplot(df, aes(x=time, y=TumorVol)) +
      geom_line(color="#C0392B", size=1) +
      geom_hline(yintercept=thresh, linetype="dashed", color="gray40") +
      annotate("text", x=max(df$time)*0.7, y=thresh*1.05,
               label=paste("PD 기준 (2×기저, 2×", input$TV0, "cm³)"), size=3, color="gray40") +
      labs(title="종양 부피 (cm³)", x="일(Day)", y="종양 부피 (cm³)") +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$plot_tv_pct <- renderPlotly({
    df <- sim_result()
    tv0 <- input$TV0
    df$pct_chg <- (df$TumorVol - tv0) / tv0 * 100
    p <- ggplot(df, aes(x=time, y=pct_chg)) +
      geom_line(color="#E65100") +
      geom_hline(yintercept=c(-30, 0, 20), linetype="dashed",
                 color=c("#27AE60","black","#C0392B")) +
      labs(title="종양 부피 변화율 (베이스라인 대비 %)",
           x="일(Day)", y="변화율 (%)") +
      theme_bw(base_size=10)
    ggplotly(p)
  })
  output$recist_table <- renderTable({
    df <- sim_result()
    tv0 <- input$TV0
    min_tv <- min(df$TumorVol)
    pct <- (min_tv - tv0) / tv0 * 100
    cat_r <- if(pct < -100) "CR (완전 관해)" else if(pct < -30) "PR (부분 관해)" else if(pct < 20) "SD (안정)" else "PD (진행)"
    data.frame(
      구분   = c("최선 반응 %", "RECIST 범주", "CA-125 최저값"),
      값     = c(paste0(round(pct,1), "%"), cat_r, paste0(round(min(df$CA125_lvl)), " U/mL"))
    )
  })

  ## ── Tab 5: Scenario Comparison ────────────────────────────────
  output$plot_scenario_tv <- renderPlotly({
    df <- all_scenarios()
    p <- ggplot(df, aes(x=time, y=TumorVol, color=Scenario)) +
      geom_line(size=0.8) +
      geom_hline(yintercept=100, linetype="dashed", color="gray50") +
      labs(title="6가지 시나리오 — 종양 부피 비교 (2년)",
           x="일(Day)", y="종양 부피 (cm³)") +
      scale_color_brewer(palette="Set1") +
      theme_bw(base_size=10) + theme(legend.position="bottom")
    ggplotly(p) %>% layout(legend=list(orientation="h", y=-0.3))
  })
  output$plot_scenario_ca125 <- renderPlotly({
    df <- all_scenarios()
    p <- ggplot(df, aes(x=time, y=CA125_lvl, color=Scenario)) +
      geom_line(size=0.8) +
      geom_hline(yintercept=35, linetype="dashed", color="darkgreen") +
      scale_y_log10() +
      labs(title="6가지 시나리오 — CA-125 (log scale)",
           x="일(Day)", y="CA-125 (U/mL)") +
      scale_color_brewer(palette="Set1") +
      theme_bw(base_size=10) + theme(legend.position="bottom")
    ggplotly(p) %>% layout(legend=list(orientation="h", y=-0.3))
  })
  output$scenario_table <- renderDT({
    df <- all_scenarios()
    df %>% group_by(Scenario) %>% summarise(
      최소종양(cm3) = round(min(TumorVol, na.rm=TRUE), 1),
      최소CA125(U_mL) = round(min(CA125_lvl, na.rm=TRUE), 1),
      최대HRD손상 = round(max(HRD_d, na.rm=TRUE), 3),
      최종CD8T = round(last(CD8T_r), 3)
    ) %>% datatable(options=list(pageLength=6, dom="t"), rownames=FALSE)
  })

  ## ── Tab 6: Biomarker Panel ────────────────────────────────────
  output$plot_bm_panel <- renderPlotly({
    df <- sim_result()
    df_long <- df %>%
      select(time, CA125_lvl, TumorVol, PtDNA_r, HRD_d, CD8T_r, VEGF_f) %>%
      pivot_longer(-time, names_to="Biomarker", values_to="Value") %>%
      mutate(Biomarker = recode(Biomarker,
        CA125_lvl="CA-125 (U/mL)", TumorVol="종양 부피 (cm³)",
        PtDNA_r="Pt-DNA 부가물", HRD_d="HRD 손상",
        CD8T_r="CD8+ T세포", VEGF_f="유리 VEGF-A (ng/mL)"))
    p <- ggplot(df_long, aes(x=time, y=Value)) +
      geom_line(aes(color=Biomarker)) +
      facet_wrap(~Biomarker, scales="free_y", ncol=3) +
      labs(title="종합 바이오마커 패널", x="일(Day)", y="값") +
      theme_bw(base_size=9) + theme(legend.position="none")
    ggplotly(p)
  })

  output$trial_ref_table <- renderTable({
    data.frame(
      임상시험 = c("SOLO-1 (2018 NEJM)","PRIMA (2019 NEJM)","PAOLA-1 (2019 NEJM)","ICON7 (2011 NEJM)","AGO-OVAR16 (2014)"),
      치료 = c("Olaparib 유지 vs 위약 (BRCA+ 1차)","Niraparib 유지 vs 위약 (1차 CR/PR)","Ola+Bev vs 위약+Bev (HRD+ 1차)","Carbo+Pacli+Bev vs 표준","Pazopanib 유지 vs 위약"),
      mPFS = c("NR vs 13.8mo (HR 0.30)","13.8mo vs 8.2mo (HR 0.43, HRD+)","22.1mo vs 16.6mo (HR 0.33, HRD+)","19.0mo vs 17.4mo","17.9mo vs 12.3mo"),
      주_대상 = c("BRCA mut 1차","All-comer HRD+/−","HRD+ HRD−","전체 1차","1차 표준 후")
    )
  })
}

shinyApp(ui, server)
