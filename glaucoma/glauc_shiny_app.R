## ============================================================
## Glaucoma QSP — Shiny Dashboard
## 녹내장 QSP 인터랙티브 시뮬레이션 대시보드
##
## Tabs:
##   1. 환자 프로파일 (Patient Profile & Risk Assessment)
##   2. 약물 PK (Drug PK — Anterior Chamber Concentrations)
##   3. 방수 역학 & 안압 (Aqueous Dynamics & IOP)
##   4. 시신경 & RGC (Optic Nerve & RGC Survival)
##   5. 임상 엔드포인트 (Clinical Endpoints: RNFL, MD, VFI)
##   6. 치료 시나리오 비교 (Scenario Comparison)
##   7. 바이오마커 패널 (Biomarker Panel & Monitoring)
##   8. 민감도 분석 (Sensitivity Analysis — Tornado)
##
## Requirements:
##   install.packages(c("shiny","shinydashboard","mrgsolve",
##                      "ggplot2","dplyr","tidyr","plotly",
##                      "DT","patchwork","shinyWidgets"))
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)
library(shinyWidgets)

## ─── mrgsolve model (inline) ─────────────────────────
glauc_code <- '
$PARAM
ka_LAT=0.46, ke_LAT=0.35, ka_TIM=0.50, ke_TIM=0.25,
ka_DZL=0.42, ke_DZL=0.20, ka_BRI=0.48, ke_BRI=0.30,
ka_NET=0.45, ke_NET=0.28,
EC50_FP=0.5, Emax_FP=1.0, EC50_BET=50, Emax_BET=1.0,
EC50_ALP=1.5, Emax_ALP=1.0, EC50_CA=80, Emax_CA=0.9,
EC50_ROC=0.3, Emax_ROC=1.0,
kout_FP=0.15, kout_BET=0.20, kout_ALP=0.18,
kout_CA=0.10, kout_ROC=0.25,
F_base=2.50, C_base=0.22, Fu_base=0.40, EVP=9.0,
kECM_on=0.08, kECM_off=0.04, dFu_ECM=0.22,
dF_cAMP=0.60, dF_CA=0.18, dC_ROCK=0.07,
MAP=90.0, IOP_thr=14.0,
kRGC=1.46e-5, kONH=0.055,
kMD=0.045, tau_RNFL=2160,
RNFL_0=105.0, VF_MD0=-2.5

$CMT LAT_D TIM_D DZL_D BRI_D NET_D
    LAT_A TIM_A DZL_A BRI_A NET_A
    FP_R BETA_R ALP_R CA_I ROCK_I
    ECM cAMP_C AQ_P C_T FU_V
    IOP_C RGC_P RNFL_U VF_M OPP_C

$INIT LAT_D=0, TIM_D=0, DZL_D=0, BRI_D=0, NET_D=0,
      LAT_A=0, TIM_A=0, DZL_A=0, BRI_A=0, NET_A=0,
      FP_R=0, BETA_R=0, ALP_R=0, CA_I=0, ROCK_I=0,
      ECM=0, cAMP_C=1.0, AQ_P=2.50, C_T=0.22, FU_V=0.40,
      IOP_C=24.0, RGC_P=100.0, RNFL_U=105.0, VF_M=-2.5, OPP_C=36.0

$ODE
dxdt_LAT_D = -ka_LAT*LAT_D;
dxdt_LAT_A = ka_LAT*LAT_D*0.04 - ke_LAT*LAT_A;
dxdt_TIM_D = -ka_TIM*TIM_D;
dxdt_TIM_A = ka_TIM*TIM_D - ke_TIM*TIM_A;
dxdt_DZL_D = -ka_DZL*DZL_D;
dxdt_DZL_A = ka_DZL*DZL_D - ke_DZL*DZL_A;
dxdt_BRI_D = -ka_BRI*BRI_D;
dxdt_BRI_A = ka_BRI*BRI_D - ke_BRI*BRI_A;
dxdt_NET_D = -ka_NET*NET_D;
dxdt_NET_A = ka_NET*NET_D - ke_NET*NET_A;

double FP_SS=Emax_FP*LAT_A/(EC50_FP+LAT_A);
double BE_SS=Emax_BET*TIM_A/(EC50_BET+TIM_A);
double AL_SS=Emax_ALP*BRI_A/(EC50_ALP+BRI_A);
double CA_SS=Emax_CA*DZL_A/(EC50_CA+DZL_A);
double RO_SS=Emax_ROC*NET_A/(EC50_ROC+NET_A);
dxdt_FP_R  = kout_FP*(FP_SS-FP_R);
dxdt_BETA_R= kout_BET*(BE_SS-BETA_R);
dxdt_ALP_R = kout_ALP*(AL_SS-ALP_R);
dxdt_CA_I  = kout_CA*(CA_SS-CA_I);
dxdt_ROCK_I= kout_ROC*(RO_SS-ROCK_I);

dxdt_ECM   = kECM_on*FP_R*(1-ECM) - kECM_off*ECM;
double cT  = 1.0-dF_cAMP*(BETA_R+0.5*ALP_R)/(1+0.5*ALP_R);
if(cT<0.1) cT=0.1;
dxdt_cAMP_C= 0.2*(cT-cAMP_C);
double Ff  = cAMP_C*(1-dF_CA*CA_I);
if(Ff<0.30) Ff=0.30;
dxdt_AQ_P  = 0.10*(F_base*Ff - AQ_P);
dxdt_C_T   = 0.10*(C_base + dC_ROCK*ROCK_I - C_T);
dxdt_FU_V  = 0.10*(Fu_base + dFu_ECM*ECM + 0.06*ALP_R - FU_V);

double IOP_c=(AQ_P-FU_V)/C_T+EVP;
if(IOP_c<4) IOP_c=4;
dxdt_IOP_C = 6.0*(IOP_c-IOP_C);

double OPP_c=(2.0/3.0)*MAP-IOP_C;
if(OPP_c<0) OPP_c=0;
dxdt_OPP_C = 3.0*(OPP_c-OPP_C);

double Ie=IOP_C-IOP_thr; if(Ie<0) Ie=0;
double OPPd=0; if(OPP_C<40) OPPd=(40-OPP_C)/40.0;
dxdt_RGC_P = -kRGC*(Ie*Ie+0.3*OPPd)*RGC_P;
dxdt_RNFL_U= (RNFL_0*(RGC_P/100.0)-RNFL_U)/tau_RNFL;
double VFt = VF_MD0 - kMD*(RNFL_0-RNFL_U);
if(VFt<-30) VFt=-30;
dxdt_VF_M  = 0.003*(VFt-VF_M);
dxdt_RGC_P = dxdt_RGC_P * (RGC_P>0 ? 1 : 0);

$TABLE
double IOP_out=IOP_C, OPP_out=OPP_C, RGC_out=RGC_P;
double RNFL_out=RNFL_U, MD_out=VF_M;
double VFI_out=(VF_M>-30)?100.0*(1.0+VF_M/30.0):0.0;
double AQ_out=AQ_P, CT_out=C_T, FU_out=FU_V;
double ECM_out=ECM, cAMP_out=cAMP_C;
double FP_out=FP_R, BET_out=BETA_R, ALP_out=ALP_R;
double CA_out=CA_I, ROC_out=ROCK_I;

$CAPTURE IOP_out OPP_out RGC_out RNFL_out MD_out VFI_out
         AQ_out CT_out FU_out ECM_out cAMP_out
         FP_out BET_out ALP_out CA_out ROC_out
         LAT_A TIM_A DZL_A BRI_A NET_A
'

mod <- mcode("glauc_shiny", glauc_code, quiet = TRUE)

## ─── Simulation function ─────────────────────────────
run_sim <- function(input, n_years = 5) {
  n_hrs <- n_years * 365 * 24
  evs   <- list()

  if (input$use_lat) {
    evs[["lat"]] <- ev(amt = input$dose_lat, cmt = 1,
                       time = 20, ii = 24,
                       addl = n_years * 365 - 1)
  }
  if (input$use_tim) {
    evs[["tim"]] <- ev(amt = input$dose_tim, cmt = 2,
                       time = 0, ii = 12,
                       addl = n_years * 365 * 2 - 1)
  }
  if (input$use_dzl) {
    evs[["dzl"]] <- ev(amt = input$dose_dzl, cmt = 3,
                       time = 0, ii = 8,
                       addl = n_years * 365 * 3 - 1)
  }
  if (input$use_bri) {
    evs[["bri"]] <- ev(amt = input$dose_bri, cmt = 4,
                       time = 0, ii = 8,
                       addl = n_years * 365 * 3 - 1)
  }
  if (input$use_net) {
    evs[["net"]] <- ev(amt = input$dose_net, cmt = 5,
                       time = 20, ii = 24,
                       addl = n_years * 365 - 1)
  }

  if (length(evs) == 0) {
    e <- ev(amt = 0, cmt = 1, time = 0)
  } else {
    e <- Reduce(c, evs)
  }

  mod %>%
    param(MAP = input$pt_MAP,
          IOP_thr = input$pt_IOP_thr,
          VF_MD0  = input$pt_vf_init,
          RNFL_0  = input$pt_rnfl_init) %>%
    mrgsim(events = e, end = n_hrs, delta = 24) %>%
    as.data.frame() %>%
    mutate(time_yr = time / (365 * 24))
}

## ─── Scenario comparison ─────────────────────────────
scenarios_fixed <- list(
  "무치료"                = list(use_lat=F,use_tim=F,use_dzl=F,use_bri=F,use_net=F,
                                 pt_MAP=90, pt_IOP_thr=14, pt_vf_init=-2.5, pt_rnfl_init=105,
                                 dose_lat=2500,dose_tim=25000,dose_dzl=50000,
                                 dose_bri=5000,dose_net=500),
  "라타노프로스트 QD"     = list(use_lat=T,use_tim=F,use_dzl=F,use_bri=F,use_net=F,
                                 pt_MAP=90,pt_IOP_thr=14,pt_vf_init=-2.5,pt_rnfl_init=105,
                                 dose_lat=2500,dose_tim=25000,dose_dzl=50000,
                                 dose_bri=5000,dose_net=500),
  "티몰롤 BID"            = list(use_lat=F,use_tim=T,use_dzl=F,use_bri=F,use_net=F,
                                 pt_MAP=90,pt_IOP_thr=14,pt_vf_init=-2.5,pt_rnfl_init=105,
                                 dose_lat=2500,dose_tim=25000,dose_dzl=50000,
                                 dose_bri=5000,dose_net=500),
  "라타노+티몰롤 병용"    = list(use_lat=T,use_tim=T,use_dzl=F,use_bri=F,use_net=F,
                                 pt_MAP=90,pt_IOP_thr=14,pt_vf_init=-2.5,pt_rnfl_init=105,
                                 dose_lat=2500,dose_tim=25000,dose_dzl=50000,
                                 dose_bri=5000,dose_net=500),
  "3제 병용"              = list(use_lat=T,use_tim=T,use_dzl=T,use_bri=F,use_net=F,
                                 pt_MAP=90,pt_IOP_thr=14,pt_vf_init=-2.5,pt_rnfl_init=105,
                                 dose_lat=2500,dose_tim=25000,dose_dzl=50000,
                                 dose_bri=5000,dose_net=500),
  "네타수딜 QD"           = list(use_lat=F,use_tim=F,use_dzl=F,use_bri=F,use_net=T,
                                 pt_MAP=90,pt_IOP_thr=14,pt_vf_init=-2.5,pt_rnfl_init=105,
                                 dose_lat=2500,dose_tim=25000,dose_dzl=50000,
                                 dose_bri=5000,dose_net=500),
  "네타수딜+라타노(Rocklatan)" = list(use_lat=T,use_tim=F,use_dzl=F,use_bri=F,use_net=T,
                                 pt_MAP=90,pt_IOP_thr=14,pt_vf_init=-2.5,pt_rnfl_init=105,
                                 dose_lat=2500,dose_tim=25000,dose_dzl=50000,
                                 dose_bri=5000,dose_net=500)
)

## ─── UI ──────────────────────────────────────────────
ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "녹내장 QSP 대시보드 (GLAUC)"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("환자 프로파일", tabName = "tab_patient",   icon = icon("user-md")),
      menuItem("약물 PK",       tabName = "tab_pk",        icon = icon("capsules")),
      menuItem("방수 & 안압",   tabName = "tab_aq",        icon = icon("eye")),
      menuItem("시신경 & RGC",  tabName = "tab_rgc",       icon = icon("brain")),
      menuItem("임상 엔드포인트", tabName = "tab_ep",      icon = icon("chart-line")),
      menuItem("시나리오 비교", tabName = "tab_scen",      icon = icon("layer-group")),
      menuItem("바이오마커",    tabName = "tab_bio",        icon = icon("flask")),
      menuItem("민감도 분석",   tabName = "tab_sa",         icon = icon("sliders-h"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box { border-top: 3px solid #28a745; }
      .info-box { min-height: 80px; }
      .irs-grid-pol.small { height: 0px; }
    "))),

    tabItems(

      ## ── TAB 1: Patient Profile ──────────────────────
      tabItem("tab_patient",
        fluidRow(
          box(title = "환자 기본 정보", width = 4, status = "success",
            numericInput("pt_MAP",  "평균 동맥압 MAP (mmHg)", 90, min=50, max=130),
            numericInput("pt_IOP_thr", "IOP 손상 역치 (mmHg)", 14, min=10, max=20),
            numericInput("pt_vf_init","초기 VF MD (dB)", -2.5, min=-20, max=0),
            numericInput("pt_rnfl_init","초기 RNFL (µm)", 105, min=60, max=130),
            sliderInput("sim_years","시뮬레이션 기간 (년)", 1, 10, 5)
          ),
          box(title = "약물 선택 및 용량", width = 4, status = "success",
            checkboxInput("use_lat","라타노프로스트 QD", TRUE),
            conditionalPanel("input.use_lat",
              sliderInput("dose_lat","용량 (ng)", 500, 5000, 2500, step=100)),
            checkboxInput("use_tim","티몰롤 BID", FALSE),
            conditionalPanel("input.use_tim",
              sliderInput("dose_tim","용량 (ng)", 5000, 50000, 25000, step=1000)),
            checkboxInput("use_dzl","도르졸아마이드 TID", FALSE),
            conditionalPanel("input.use_dzl",
              sliderInput("dose_dzl","용량 (ng)", 10000, 100000, 50000, step=5000)),
            checkboxInput("use_bri","브리모니딘 TID", FALSE),
            conditionalPanel("input.use_bri",
              sliderInput("dose_bri","용량 (ng)", 1000, 10000, 5000, step=500)),
            checkboxInput("use_net","네타수딜 QD", FALSE),
            conditionalPanel("input.use_net",
              sliderInput("dose_net","용량 (ng)", 100, 1000, 500, step=50)),
            actionButton("run_sim", "시뮬레이션 실행", class = "btn-success btn-block")
          ),
          box(title = "OHTS 위험 계산기", width = 4, status = "warning",
            numericInput("ohts_age","나이 (세)", 55, 40, 80),
            numericInput("ohts_iop","IOP (mmHg)", 24, 21, 40),
            numericInput("ohts_cct","CCT (µm)", 545, 450, 650),
            numericInput("ohts_cd","C/D 비율", 0.5, 0.2, 0.9, step=0.05),
            numericInput("ohts_psd","HVF PSD (dB)", 2.0, 0.5, 6.0, step=0.1),
            verbatimTextOutput("ohts_score")
          )
        ),
        fluidRow(
          infoBoxOutput("ib_iop",  width=3),
          infoBoxOutput("ib_rnfl", width=3),
          infoBoxOutput("ib_md",   width=3),
          infoBoxOutput("ib_rgc",  width=3)
        )
      ),

      ## ── TAB 2: Drug PK ─────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(title = "전방 약물 농도 (AC Concentrations)", width = 12,
            plotlyOutput("plot_pk_ac", height = "400px")
          )
        ),
        fluidRow(
          box(title = "수용체/표적 점유율 (Receptor Occupancy)", width = 6,
            plotlyOutput("plot_ro", height = "350px")
          ),
          box(title = "PK 파라미터 요약", width = 6,
            DTOutput("tbl_pk_params")
          )
        )
      ),

      ## ── TAB 3: Aqueous Dynamics & IOP ──────────────
      tabItem("tab_aq",
        fluidRow(
          box(title = "안압 (IOP) 시간 경과", width = 6,
            plotlyOutput("plot_iop", height = "350px")
          ),
          box(title = "방수 역학 (Aqueous Humor Dynamics)", width = 6,
            plotlyOutput("plot_aq_dyn", height = "350px")
          )
        ),
        fluidRow(
          box(title = "소주 유출 용이도 & 포도막 유출", width = 6,
            plotlyOutput("plot_outflow", height = "300px")
          ),
          box(title = "Goldman 방정식 분석", width = 6,
            plotlyOutput("plot_goldman", height = "300px")
          )
        )
      ),

      ## ── TAB 4: Optic Nerve & RGC ───────────────────
      tabItem("tab_rgc",
        fluidRow(
          box(title = "RGC 생존율 (%)", width = 6,
            plotlyOutput("plot_rgc", height = "350px")
          ),
          box(title = "ONH 기계적 스트레스", width = 6,
            plotlyOutput("plot_onh", height = "350px")
          )
        ),
        fluidRow(
          box(title = "안구 관류압 (OPP)", width = 6,
            plotlyOutput("plot_opp", height = "300px")
          ),
          box(title = "치료 효과 — RGC 보존 (%)", width = 6,
            plotlyOutput("plot_rgc_5yr", height = "300px")
          )
        )
      ),

      ## ── TAB 5: Clinical Endpoints ──────────────────
      tabItem("tab_ep",
        fluidRow(
          box(title = "RNFL 두께 (OCT)", width = 6,
            plotlyOutput("plot_rnfl", height = "350px")
          ),
          box(title = "시야 MD (dB)", width = 6,
            plotlyOutput("plot_md", height = "350px")
          )
        ),
        fluidRow(
          box(title = "VFI (%) — Visual Field Index", width = 6,
            plotlyOutput("plot_vfi", height = "300px")
          ),
          box(title = "임상 엔드포인트 요약표", width = 6,
            DTOutput("tbl_ep")
          )
        )
      ),

      ## ── TAB 6: Scenario Comparison ─────────────────
      tabItem("tab_scen",
        fluidRow(
          box(title = "7가지 치료 시나리오 비교 — IOP", width = 6,
            plotlyOutput("plot_scen_iop", height = "380px")
          ),
          box(title = "7가지 치료 시나리오 비교 — RNFL", width = 6,
            plotlyOutput("plot_scen_rnfl", height = "380px")
          )
        ),
        fluidRow(
          box(title = "시야 MD 비교", width = 6,
            plotlyOutput("plot_scen_md", height = "350px")
          ),
          box(title = "5년 결과 비교표", width = 6,
            DTOutput("tbl_scen_5yr")
          )
        )
      ),

      ## ── TAB 7: Biomarker Panel ─────────────────────
      tabItem("tab_bio",
        fluidRow(
          box(title = "방수 생성률 & cAMP", width = 6,
            plotlyOutput("plot_bio_aq", height = "320px")
          ),
          box(title = "ECM 재형성 상태 (라타노프로스트)", width = 6,
            plotlyOutput("plot_bio_ecm", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Biomarker Radar Chart (1년 시점)", width = 6,
            plotlyOutput("plot_radar", height = "350px")
          ),
          box(title = "IOP Peak & Trough 분석", width = 6,
            plotlyOutput("plot_pk_pt", height = "350px")
          )
        )
      ),

      ## ── TAB 8: Sensitivity Analysis ────────────────
      tabItem("tab_sa",
        fluidRow(
          box(title = "파라미터 선택", width = 4, status = "info",
            sliderInput("sa_kRGC", "kRGC (×배수)", 0.5, 3.0, 1.0, step=0.1),
            sliderInput("sa_dFu",  "dFu_ECM 배수", 0.5, 2.0, 1.0, step=0.1),
            sliderInput("sa_dC",   "dC_ROCK 배수",  0.5, 2.0, 1.0, step=0.1),
            sliderInput("sa_EVP",  "EVP (mmHg)", 6, 14, 9),
            sliderInput("sa_MAP",  "MAP (mmHg)", 70, 110, 90),
            actionButton("run_sa", "민감도 분석 실행", class = "btn-info btn-block")
          ),
          box(title = "토네이도 차트 — 5년 IOP", width = 8,
            plotlyOutput("plot_tornado", height = "450px")
          )
        )
      )
    )
  )
)

## ─── Server ──────────────────────────────────────────
server <- function(input, output, session) {

  ## Reactive simulation
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "시뮬레이션 실행 중...", {
      run_sim(input, input$sim_years)
    })
  }, ignoreNULL = FALSE)

  ## Scenario comparison data (precomputed on load)
  scen_data <- reactive({
    withProgress(message = "시나리오 비교 계산 중...", {
      purrr::map2_dfr(names(scenarios_fixed), scenarios_fixed, function(nm, params) {
        s <- structure(params, class = "list")
        run_sim(s, 5) %>% mutate(scenario = nm)
      })
    })
  })

  ## ── OHTS Risk Score ──────────────────────────────
  output$ohts_score <- renderText({
    age  <- input$ohts_age
    iop  <- input$ohts_iop
    cct  <- input$ohts_cct
    cd   <- input$ohts_cd
    psd  <- input$ohts_psd
    # Simplified OHTS logistic model
    lp <- -7.10 + 0.10*age + 0.09*iop - 0.009*cct + 3.12*cd + 0.31*psd
    risk5yr <- 100 / (1 + exp(-lp))
    cat(sprintf("5년 POAG 발생 위험: %.1f%%\n", risk5yr),
        sprintf("IOP: %d mmHg | CCT: %d µm | C/D: %.2f\n", iop, cct, cd),
        if(risk5yr > 15) "치료 권고: 강력 권고\n"
        else if(risk5yr > 8) "치료 권고: 권고\n"
        else "치료 권고: 경과 관찰 가능\n")
  })

  ## ── Info Boxes ───────────────────────────────────
  output$ib_iop <- renderInfoBox({
    d <- sim_data(); last <- tail(d, 1)
    iop_val <- round(last$IOP_out, 1)
    color <- if(iop_val > 21) "red" else if(iop_val > 15) "yellow" else "green"
    infoBox("현재 IOP", paste0(iop_val, " mmHg"), icon=icon("eye"), color=color)
  })
  output$ib_rnfl <- renderInfoBox({
    d <- sim_data(); last <- tail(d, 1)
    v <- round(last$RNFL_out, 1)
    color <- if(v < 75) "red" else if(v < 90) "yellow" else "green"
    infoBox("RNFL 두께", paste0(v, " µm"), icon=icon("layer-group"), color=color)
  })
  output$ib_md <- renderInfoBox({
    d <- sim_data(); last <- tail(d, 1)
    v <- round(last$MD_out, 2)
    color <- if(v < -12) "red" else if(v < -6) "yellow" else "green"
    infoBox("시야 MD", paste0(v, " dB"), icon=icon("low-vision"), color=color)
  })
  output$ib_rgc <- renderInfoBox({
    d <- sim_data(); last <- tail(d, 1)
    v <- round(last$RGC_out, 1)
    color <- if(v < 70) "red" else if(v < 85) "yellow" else "green"
    infoBox("RGC 생존율", paste0(v, "%"), icon=icon("microscope"), color=color)
  })

  ## ── PK Plots ─────────────────────────────────────
  output$plot_pk_ac <- renderPlotly({
    d <- sim_data() %>% filter(time_yr < 0.1)
    p <- d %>%
      select(time_yr, LAT_A, TIM_A, DZL_A, BRI_A, NET_A) %>%
      pivot_longer(-time_yr) %>%
      ggplot(aes(time_yr*8760, value, color=name)) +
      geom_line() +
      labs(x="시간 (h)", y="전방 농도 (ng/mL)") +
      theme_minimal()
    ggplotly(p)
  })

  output$plot_ro <- renderPlotly({
    d <- sim_data() %>% filter(time_yr <= min(input$sim_years, 1))
    p <- d %>%
      select(time_yr, FP_out, BET_out, ALP_out, CA_out, ROC_out) %>%
      pivot_longer(-time_yr) %>%
      ggplot(aes(time_yr, value, color=name)) +
      geom_line() +
      ylim(0, 1) +
      labs(x="시간 (년)", y="점유율 (0–1)") +
      theme_minimal()
    ggplotly(p)
  })

  output$tbl_pk_params <- renderDT({
    data.frame(
      약물 = c("라타노프로스트","티몰롤","도르졸아마이드","브리모니딘","네타수딜"),
      표적 = c("FP 수용체","β2-AR","CA-II","α2-AR","ROCK1/2"),
      EC50_ngmL = c(0.5, 50, 80, 1.5, 0.3),
      IOP감소 = c("25–32%","20–25%","15–20%","18–20%","15–20%"),
      주요작용 = c("포도막유출↑","방수생성↓","방수생성↓","방수생성↓","소주유출↑")
    )
  }, options = list(pageLength = 5))

  ## ── Aqueous / IOP plots ──────────────────────────
  output$plot_iop <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_yr, IOP_out)) +
      geom_line(color="#2ecc71", size=1.2) +
      geom_hline(yintercept=21, linetype="dashed") +
      geom_hline(yintercept=15, linetype="dotted", color="#e74c3c") +
      labs(x="시간 (년)", y="IOP (mmHg)", title="안압 변화") +
      theme_minimal()
    ggplotly(p)
  })

  output$plot_aq_dyn <- renderPlotly({
    d <- sim_data()
    p <- d %>% select(time_yr, AQ_out, cAMP_out) %>%
      pivot_longer(-time_yr) %>%
      ggplot(aes(time_yr, value, color=name)) +
      geom_line(size=1) +
      labs(x="시간 (년)", y="정규화 값") +
      theme_minimal()
    ggplotly(p)
  })

  output$plot_outflow <- renderPlotly({
    d <- sim_data()
    p <- d %>% select(time_yr, CT_out, FU_out) %>%
      pivot_longer(-time_yr) %>%
      ggplot(aes(time_yr, value, color=name)) +
      geom_line(size=1) +
      labs(x="시간 (년)", y="유출 파라미터") +
      theme_minimal()
    ggplotly(p)
  })

  output$plot_goldman <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(AQ_out, IOP_out, color=time_yr)) +
      geom_path() + geom_point(size=0.8) +
      scale_color_viridis_c() +
      labs(x="방수 생성률 (µL/min)", y="IOP (mmHg)",
           title="Goldman 관계도") +
      theme_minimal()
    ggplotly(p)
  })

  ## ── RGC / ONH plots ──────────────────────────────
  output$plot_rgc <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_yr, RGC_out)) +
      geom_line(color="#e74c3c", size=1.2) +
      geom_hline(yintercept=50, linetype="dashed") +
      labs(x="시간 (년)", y="RGC 생존율 (%)", title="RGC 생존율") +
      theme_minimal()
    ggplotly(p)
  })

  output$plot_onh <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_yr, IOP_out - 14)) +
      geom_line(color="#9b59b6", size=1.2) +
      geom_hline(yintercept=0, linetype="dashed") +
      labs(x="시간 (년)", y="IOP 초과분 (mmHg)", title="사상판 압박 대리 지표") +
      theme_minimal()
    ggplotly(p)
  })

  output$plot_opp <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_yr, OPP_out)) +
      geom_line(color="#3498db", size=1.2) +
      geom_hline(yintercept=40, linetype="dashed", color="red") +
      labs(x="시간 (년)", y="OPP (mmHg)", title="안구 관류압") +
      theme_minimal()
    ggplotly(p)
  })

  output$plot_rgc_5yr <- renderPlotly({
    d <- scen_data() %>% filter(abs(time_yr - 5) < 0.1) %>%
      group_by(scenario) %>%
      summarise(RGC = mean(RGC_out), .groups="drop")
    p <- ggplot(d, aes(reorder(scenario, RGC), RGC, fill=scenario)) +
      geom_col() + coord_flip() +
      labs(x="", y="5년 RGC 생존율 (%)") +
      theme_minimal() + theme(legend.position="none")
    ggplotly(p)
  })

  ## ── Clinical Endpoint plots ──────────────────────
  output$plot_rnfl <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_yr, RNFL_out)) +
      geom_line(color="#f39c12", size=1.2) +
      geom_hline(yintercept=80, linetype="dashed") +
      labs(x="시간 (년)", y="RNFL 두께 (µm)", title="RNFL 두께") +
      theme_minimal()
    ggplotly(p)
  })

  output$plot_md <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_yr, MD_out)) +
      geom_line(color="#c0392b", size=1.2) +
      geom_hline(yintercept=-6, linetype="dashed") +
      geom_hline(yintercept=-12, linetype="dotted", color="red") +
      labs(x="시간 (년)", y="MD (dB)", title="시야 평균 편차 (MD)") +
      theme_minimal()
    ggplotly(p)
  })

  output$plot_vfi <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(time_yr, VFI_out)) +
      geom_line(color="#27ae60", size=1.2) +
      geom_hline(yintercept=70, linetype="dashed") +
      ylim(0, 100) +
      labs(x="시간 (년)", y="VFI (%)", title="시야 기능 지수 (VFI)") +
      theme_minimal()
    ggplotly(p)
  })

  output$tbl_ep <- renderDT({
    sim_data() %>%
      filter(time_yr %in% c(0, 1, 2, 3, 5)) %>%
      mutate(across(where(is.numeric), ~round(., 1))) %>%
      select(year = time_yr, IOP = IOP_out, RNFL = RNFL_out,
             MD = MD_out, VFI = VFI_out, RGC = RGC_out)
  }, options = list(pageLength=6))

  ## ── Scenario plots ───────────────────────────────
  output$plot_scen_iop <- renderPlotly({
    d <- scen_data()
    p <- ggplot(d, aes(time_yr, IOP_out, color=scenario)) +
      geom_line() +
      geom_hline(yintercept=21, linetype="dashed") +
      labs(x="년", y="IOP (mmHg)") + theme_minimal() +
      theme(legend.position="bottom")
    ggplotly(p)
  })

  output$plot_scen_rnfl <- renderPlotly({
    d <- scen_data()
    p <- ggplot(d, aes(time_yr, RNFL_out, color=scenario)) +
      geom_line() +
      geom_hline(yintercept=80, linetype="dashed") +
      labs(x="년", y="RNFL (µm)") + theme_minimal() +
      theme(legend.position="bottom")
    ggplotly(p)
  })

  output$plot_scen_md <- renderPlotly({
    d <- scen_data()
    p <- ggplot(d, aes(time_yr, MD_out, color=scenario)) +
      geom_line() +
      labs(x="년", y="MD (dB)") + theme_minimal() +
      theme(legend.position="bottom")
    ggplotly(p)
  })

  output$tbl_scen_5yr <- renderDT({
    scen_data() %>%
      filter(abs(time_yr - 5) < 0.1) %>%
      group_by(scenario) %>%
      summarise(
        IOP  = round(mean(IOP_out),1),
        RNFL = round(mean(RNFL_out),1),
        MD   = round(mean(MD_out),2),
        VFI  = round(mean(VFI_out),1),
        RGC  = round(mean(RGC_out),1),
        .groups="drop"
      ) %>%
      arrange(desc(IOP))
  })

  ## ── Biomarker plots ──────────────────────────────
  output$plot_bio_aq <- renderPlotly({
    d <- sim_data()
    p <- d %>% select(time_yr, AQ_out, cAMP_out, CA_out) %>%
      pivot_longer(-time_yr) %>%
      ggplot(aes(time_yr, value, color=name)) +
      geom_line(size=0.9) +
      labs(x="년", y="정규화 값") + theme_minimal()
    ggplotly(p)
  })

  output$plot_bio_ecm <- renderPlotly({
    d <- sim_data()
    p <- d %>% select(time_yr, ECM_out, FP_out, FU_out) %>%
      pivot_longer(-time_yr) %>%
      ggplot(aes(time_yr, value, color=name)) +
      geom_line(size=0.9) +
      labs(x="년", y="값") + theme_minimal()
    ggplotly(p)
  })

  output$plot_radar <- renderPlotly({
    d <- sim_data() %>%
      filter(abs(time_yr - 1) < 0.05) %>%
      slice(1)
    categories <- c("IOP 감소","RNFL 보존","VFI","RGC","OPP","방수생성감소")
    iop_base <- 24
    vals <- c(
      100*(iop_base - d$IOP_out)/iop_base,
      d$RNFL_out/105*100,
      d$VFI_out,
      d$RGC_out,
      min(d$OPP_out/50*100, 100),
      (1 - d$AQ_out/2.5)*100
    )
    vals <- pmax(vals, 0)
    plot_ly(type = "scatterpolar", mode = "lines+markers",
            r = c(vals, vals[1]), theta = c(categories, categories[1]),
            fill = "toself", fillcolor = "rgba(40,167,69,0.3)",
            line = list(color = "#28a745")) %>%
      layout(polar = list(radialaxis = list(range = c(0, 100))),
             title = "1년 치료 효과 레이더")
  })

  output$plot_pk_pt <- renderPlotly({
    d <- sim_data() %>% filter(time_yr < 0.5)
    p <- ggplot(d, aes(time_yr*365*24, IOP_out)) +
      geom_line(color="#2ecc71") +
      labs(x="시간 (h)", y="IOP (mmHg)", title="IOP 피크/트로프") +
      theme_minimal()
    ggplotly(p)
  })

  ## ── Sensitivity Analysis ─────────────────────────
  sa_results <- eventReactive(input$run_sa, {
    params_base <- list(kRGC=1.46e-5, dFu_ECM=0.22, dC_ROCK=0.07, EVP=9, MAP=90)
    params_vary <- list(
      kRGC   = 1.46e-5 * input$sa_kRGC,
      dFu_ECM= 0.22    * input$sa_dFu,
      dC_ROCK= 0.07    * input$sa_dC,
      EVP    = input$sa_EVP,
      MAP    = input$sa_MAP
    )
    withProgress(message = "민감도 분석...", {
      map_dfr(names(params_base), function(par) {
        plo <- params_vary; phi <- params_vary
        plo[[par]] <- params_base[[par]] * 0.5
        phi[[par]] <- params_base[[par]] * 2.0
        run_low <- mod %>% param(plo) %>%
          mrgsim(events=ev(amt=2500,cmt=1,time=20,ii=24,addl=5*365-1),
                 end=5*365*24, delta=24) %>%
          as.data.frame() %>% tail(1) %>% pull(IOP_out)
        run_high <- mod %>% param(phi) %>%
          mrgsim(events=ev(amt=2500,cmt=1,time=20,ii=24,addl=5*365-1),
                 end=5*365*24, delta=24) %>%
          as.data.frame() %>% tail(1) %>% pull(IOP_out)
        data.frame(param=par, IOP_low=run_low, IOP_high=run_high)
      })
    })
  })

  output$plot_tornado <- renderPlotly({
    d <- sa_results()
    base_iop <- 15.5
    d <- d %>%
      mutate(
        delta_low  = IOP_low  - base_iop,
        delta_high = IOP_high - base_iop
      ) %>%
      mutate(range = abs(delta_high - delta_low)) %>%
      arrange(range)
    p <- ggplot(d) +
      geom_segment(aes(x=delta_low, xend=delta_high,
                       y=reorder(param, range), yend=reorder(param, range)),
                   size=6, color="#3498db", alpha=0.7) +
      geom_vline(xintercept=0, color="red", linetype="dashed") +
      labs(x="IOP 변화 (기저치 대비 mmHg)", y="파라미터",
           title="토네이도 차트 — 5년 IOP 민감도") +
      theme_minimal()
    ggplotly(p)
  })
}

## Run app
shinyApp(ui = ui, server = server)
