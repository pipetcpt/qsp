## ============================================================
## Prurigo Nodularis QSP — Shiny Dashboard
## 6개 탭: 환자 프로파일 · PK · 면역 PD · 피부/신경 · 임상 엔드포인트 · 시나리오 비교 · 바이오마커
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)

## ─── mrgsolve model code ──────────────────────────────────────
model_code <- '
$PROB PN QSP Shiny Model

$PARAM
DUP_KA=0.0048 DUP_F1=0.64 DUP_CL=0.012 DUP_V1=3.0 DUP_Q=0.006 DUP_V2=5.5
DUP_KSYN=0.01 DUP_KDEG=0.004 DUP_KON=2.0 DUP_KOFF=0.003 DUP_KINT=0.002
NEM_KA=0.005 NEM_F1=0.72 NEM_CL=0.009 NEM_V1=2.8 NEM_Q=0.004 NEM_V2=4.5
NEM_KSYN=0.008 NEM_KDEG=0.003 NEM_KON=3.0 NEM_KOFF=0.002 NEM_KINT=0.0015
TRAL_KA=0.0042 TRAL_F1=0.77 TRAL_CL=0.010 TRAL_V1=3.2 TRAL_Q=0.005 TRAL_V2=5.0
CSA_KA=0.08 CSA_F1=0.40 CSA_CL=25.0 CSA_V=350.0
NALB_KA=0.03 NALB_F1=0.20 NALB_CL=140.0 NALB_V=800.0
TCS_KA_SK=0.02 TCS_KEL=0.08
TH2_0=100 TH2_KPROL=0.003 TH2_KDTH=0.003 TH2_EC50=2.0
IL4_0=5 IL4_KPROD=1.5 IL4_KDEG=0.30
IL13_0=8 IL13_KPROD=2.0 IL13_KDEG=0.25
IL31_0=12 IL31_KPROD=3.0 IL31_KDEG=0.20 IL31_KPROD_TH2=0.03
IGE_0=800 IGE_KPROD=0.002 IGE_KDEG=0.0001 IGE_IL4_EC50=5.0
MAST_0=100 MAST_KPROL=0.0005 MAST_KDTH=0.0005 MAST_IGE_EC50=500.0
EOS_0=400 EOS_KPROD=120 EOS_KDEG=0.30 EOS_IL13_EC50=5.0
TEWL_0=20 TEWL_MAX=80 TEWL_IL4_EC50=15 TEWL_KDEG=0.01
NERVE_0=150 NERVE_KGROWTH=0.004 NERVE_KDTH=0.003 NERVE_NGF_EC50=20
NODULE_0=20 NODULE_KFORM=0.002 NODULE_KHEAL=0.0015
DNRS_0=18
ITCH_0=7.5 ITCH_IL31_EC50=10 ITCH_KSEN=0.005 ITCH_KDES=0.003
CSA_IMAX=0.85 CSA_IC50=150
NALB_IMAX=0.70 NALB_IC50=8
TCS_IMAX=0.90 TCS_IC50=0.5
BW=75

$CMT
DUP_SC DUP_C1 DUP_C2 DUP_IL4Ra DUP_CMPLX
NEM_SC NEM_C1 NEM_C2 NEM_IL31Ra NEM_CMPLX
TRAL_SC TRAL_C1 TRAL_C2
CSA_GUT CSA_C1
NALB_GUT NALB_C1
TCS_DEPOT TCS_SKIN
TH2_CELLS IL4 IL13 IL31 IGE MAST_CELLS EOS_COUNT
TEWL NERVE_DEN NODULE_CNT ITCH_CS DNRS

$ODE
double DUP_RO = DUP_CMPLX > 0 ? DUP_CMPLX/(DUP_IL4Ra+DUP_CMPLX) : 0.0;
double INH_IL4_DUP=1.0-DUP_RO*0.95; double INH_IL13_DUP=1.0-DUP_RO*0.93;
double NEM_RO = NEM_CMPLX>0 ? NEM_CMPLX/(NEM_IL31Ra+NEM_CMPLX) : 0.0;
double INH_IL31_NEM=1.0-NEM_RO*0.90;
double TRAL_C=TRAL_C1; double INH_IL13_TRAL=1.0-(TRAL_C/(TRAL_C+2.5));
double CSA_C=CSA_C1; double INH_TH2_CSA=1.0-(CSA_IMAX*pow(CSA_C,2)/(pow(CSA_IC50,2)+pow(CSA_C,2)));
double NALB_C=NALB_C1; double INH_ITCH_NALB=1.0-(NALB_IMAX*NALB_C/(NALB_IC50+NALB_C));
double TCS_S=TCS_SKIN; double INH_INFLAM_TCS=1.0-(TCS_IMAX*TCS_S/(TCS_IC50+TCS_S));
double INH_IL13_TOTAL=INH_IL13_DUP*INH_IL13_TRAL;

dxdt_DUP_SC=-DUP_KA*DUP_SC;
dxdt_DUP_C1=DUP_KA*DUP_F1*DUP_SC/DUP_V1-DUP_CL/DUP_V1*DUP_C1-DUP_Q/DUP_V1*DUP_C1+DUP_Q/DUP_V2*DUP_C2-DUP_KON*DUP_C1*DUP_IL4Ra+DUP_KOFF*DUP_CMPLX;
dxdt_DUP_C2=DUP_Q/DUP_V1*DUP_C1-DUP_Q/DUP_V2*DUP_C2;
dxdt_DUP_IL4Ra=DUP_KSYN-DUP_KDEG*DUP_IL4Ra-DUP_KON*DUP_C1*DUP_IL4Ra+DUP_KOFF*DUP_CMPLX;
dxdt_DUP_CMPLX=DUP_KON*DUP_C1*DUP_IL4Ra-DUP_KOFF*DUP_CMPLX-DUP_KINT*DUP_CMPLX;

dxdt_NEM_SC=-NEM_KA*NEM_SC;
dxdt_NEM_C1=NEM_KA*NEM_F1*NEM_SC/NEM_V1-NEM_CL/NEM_V1*NEM_C1-NEM_Q/NEM_V1*NEM_C1+NEM_Q/NEM_V2*NEM_C2-NEM_KON*NEM_C1*NEM_IL31Ra+NEM_KOFF*NEM_CMPLX;
dxdt_NEM_C2=NEM_Q/NEM_V1*NEM_C1-NEM_Q/NEM_V2*NEM_C2;
dxdt_NEM_IL31Ra=NEM_KSYN-NEM_KDEG*NEM_IL31Ra-NEM_KON*NEM_C1*NEM_IL31Ra+NEM_KOFF*NEM_CMPLX;
dxdt_NEM_CMPLX=NEM_KON*NEM_C1*NEM_IL31Ra-NEM_KOFF*NEM_CMPLX-NEM_KINT*NEM_CMPLX;

dxdt_TRAL_SC=-TRAL_KA*TRAL_SC;
dxdt_TRAL_C1=TRAL_KA*TRAL_F1*TRAL_SC/TRAL_V1-TRAL_CL/TRAL_V1*TRAL_C1-TRAL_Q/TRAL_V1*TRAL_C1+TRAL_Q/TRAL_V2*TRAL_C2;
dxdt_TRAL_C2=TRAL_Q/TRAL_V1*TRAL_C1-TRAL_Q/TRAL_V2*TRAL_C2;

dxdt_CSA_GUT=-CSA_KA*CSA_GUT;
dxdt_CSA_C1=CSA_KA*CSA_F1*CSA_GUT/CSA_V-CSA_CL/CSA_V*CSA_C1;

dxdt_NALB_GUT=-NALB_KA*NALB_GUT;
dxdt_NALB_C1=NALB_KA*NALB_F1*NALB_GUT/NALB_V-NALB_CL/NALB_V*NALB_C1;

dxdt_TCS_DEPOT=-TCS_KA_SK*TCS_DEPOT;
dxdt_TCS_SKIN=TCS_KA_SK*TCS_DEPOT-TCS_KEL*TCS_SKIN;

double IL4_VAL=IL4>0?IL4:1e-6; double TH2_V=TH2_CELLS;
double TH2_stim=TH2_KPROL*IL4_VAL/(TH2_EC50+IL4_VAL);
dxdt_TH2_CELLS=(TH2_stim-TH2_KDTH)*TH2_CELLS*INH_TH2_CSA*INH_INFLAM_TCS;
dxdt_IL4=IL4_KPROD*(TH2_V/TH2_0)*INH_INFLAM_TCS-IL4_KDEG*IL4;
dxdt_IL13=IL13_KPROD*(TH2_V/TH2_0)*INH_INFLAM_TCS-IL13_KDEG*IL13;
dxdt_IL31=IL31_KPROD+IL31_KPROD_TH2*TH2_V-IL31_KDEG*IL31;
double IL4_NORM=IL4_VAL/(IGE_IL4_EC50+IL4_VAL);
dxdt_IGE=IGE_KPROD*(1.0+4.0*IL4_NORM)*INH_IL4_DUP*INH_IL13_TOTAL-IGE_KDEG*IGE;
double IGE_NORM=IGE/(MAST_IGE_EC50+IGE);
dxdt_MAST_CELLS=MAST_KPROL*(1.0+2.0*IGE_NORM)*MAST_CELLS-MAST_KDTH*MAST_CELLS;
double IL13_V=IL13>0?IL13:1e-6;
dxdt_EOS_COUNT=EOS_KPROD*(1.0+IL13_V/(EOS_IL13_EC50+IL13_V))*INH_IL13_TOTAL*INH_INFLAM_TCS-EOS_KDEG*EOS_COUNT;
double TEWL_stim=TEWL_MAX*IL4_NORM*INH_IL4_DUP*INH_IL13_TOTAL;
dxdt_TEWL=TEWL_KDEG*(TEWL_0+TEWL_stim-TEWL);
double NGF_proxy=20.0*(MAST_CELLS/MAST_0);
dxdt_NERVE_DEN=NERVE_KGROWTH*NGF_proxy/(NERVE_NGF_EC50+NGF_proxy)*NERVE_0-NERVE_KDTH*NERVE_DEN;
double SCRATCH_DRIVE=ITCH_CS>0?ITCH_CS:1e-6;
dxdt_NODULE_CNT=NODULE_KFORM*SCRATCH_DRIVE*(IL13_V/(5.0+IL13_V))-NODULE_KHEAL*INH_IL13_TOTAL*INH_INFLAM_TCS*NODULE_CNT;
double IL31_V=IL31>0?IL31:1e-6;
dxdt_ITCH_CS=ITCH_KSEN*(IL31_V/(ITCH_IL31_EC50+IL31_V))*(NERVE_DEN/NERVE_0)*INH_IL31_NEM*INH_ITCH_NALB-ITCH_KDES*ITCH_CS;
double DNRS_target=28.0*ITCH_CS*(NODULE_CNT/NODULE_0)*0.5;
dxdt_DNRS=0.02*(DNRS_target-DNRS);

$TABLE
double DUP_RO_PCT=DUP_CMPLX>0?100.0*DUP_CMPLX/(DUP_IL4Ra+DUP_CMPLX):0.0;
double NEM_RO_PCT=NEM_CMPLX>0?100.0*NEM_CMPLX/(NEM_IL31Ra+NEM_CMPLX):0.0;
double CSA_CONC_NG=CSA_C1*1202.0/1000.0;
double ITCH_NRS_OUT=10.0*ITCH_CS*(IL31/(ITCH_IL31_EC50+IL31));
if(ITCH_NRS_OUT>10.0)ITCH_NRS_OUT=10.0;
if(ITCH_NRS_OUT<0.0)ITCH_NRS_OUT=0.0;
double IGA_OUT=4.0*(NODULE_CNT/NODULE_0)*(TEWL/TEWL_MAX)*2.0;
if(IGA_OUT>4.0)IGA_OUT=4.0;
capture DUP_C1 NEM_C1 TRAL_C1 CSA_CONC_NG NALB_C1
capture DUP_RO_PCT NEM_RO_PCT
capture IL4 IL13 IL31 IGE EOS_COUNT MAST_CELLS TH2_CELLS
capture TEWL NERVE_DEN NODULE_CNT ITCH_CS DNRS
capture ITCH_NRS_OUT IGA_OUT

$INIT
DUP_SC=0 DUP_C1=0 DUP_C2=0 DUP_IL4Ra=2.5 DUP_CMPLX=0
NEM_SC=0 NEM_C1=0 NEM_C2=0 NEM_IL31Ra=2.667 NEM_CMPLX=0
TRAL_SC=0 TRAL_C1=0 TRAL_C2=0
CSA_GUT=0 CSA_C1=0
NALB_GUT=0 NALB_C1=0
TCS_DEPOT=0 TCS_SKIN=0
TH2_CELLS=100 IL4=5 IL13=8 IL31=12 IGE=800 MAST_CELLS=100 EOS_COUNT=400
TEWL=20 NERVE_DEN=150 NODULE_CNT=20 ITCH_CS=1.0 DNRS=18
'

mod_base <- mcode("pn_shiny", model_code, quiet=TRUE)

## ─── Simulation helper ─────────────────────────────────────────
run_pn <- function(scen_ids, duration_wk = 52, BW = 75,
                   base_il31 = 12, base_ige = 800,
                   dupilumab_dose = 300,
                   nemolizumab_dose = 60) {
  mod_run <- mod_base |>
    param(BW = BW, IL31_0 = base_il31, IGE_0 = base_ige)

  results <- lapply(scen_ids, function(sid) {
    ev_list <- list()
    n_d <- duration_wk * 7

    if (sid == 2) {
      amt <- (dupilumab_dose * 1e6 / 144000) * 3.0  # nmol in V1
      ev_list[["dup"]] <- ev(cmt="DUP_SC", amt=amt, ii=14*24,
                             addl=floor(duration_wk/2)-1, time=0)
    }
    if (sid == 3) {
      amt <- (nemolizumab_dose * 1e6 / 152000) * 2.8
      ev_list[["nem"]] <- ev(cmt="NEM_SC", amt=amt, ii=28*24,
                             addl=floor(duration_wk/4)-1, time=0)
    }
    if (sid == 4) {
      ev_list[["tral"]] <- ev(cmt="TRAL_SC", amt=625*3.2, ii=14*24,
                              addl=floor(duration_wk/2)-1, time=0)
    }
    if (sid == 5) {
      dose_mg <- BW * 5.0 / 2
      ev_list[["csa"]] <- ev(cmt="CSA_GUT", amt=dose_mg/2, ii=12,
                             addl=n_d*2-1, time=0)
    }
    if (sid == 6) {
      ev_list[["nalb"]] <- ev(cmt="NALB_GUT", amt=54, ii=12,
                              addl=n_d*2-1, time=0)
    }
    if (sid == 7) {
      amt <- (dupilumab_dose * 1e6 / 144000) * 3.0
      ev_list[["dup"]] <- ev(cmt="DUP_SC", amt=amt, ii=14*24,
                             addl=floor(duration_wk/2)-1, time=0)
      ev_list[["tcs"]] <- ev(cmt="TCS_DEPOT", amt=1.0, ii=12,
                             addl=n_d*2-1, time=0)
    }

    ev_obj <- if (length(ev_list) > 0) do.call(c, ev_list) else NULL
    out <- if (is.null(ev_obj)) {
      mrgsim(mod_run, end=duration_wk*7*24, delta=12)
    } else {
      mrgsim(mod_run, events=ev_obj, end=duration_wk*7*24, delta=12)
    }

    scen_names <- c("1"="Placebo","2"="Dupilumab","3"="Nemolizumab",
                    "4"="Tralokinumab","5"="Cyclosporine",
                    "6"="Nalbuphine ER","7"="Dupilumab+TCS")
    as.data.frame(out) |>
      mutate(scenario = scen_names[as.character(sid)],
             time_wk = time / (7*24))
  }) |> bind_rows()
  results
}

SCEN_COLORS <- c(
  "Placebo"         = "#999999",
  "Dupilumab"       = "#E41A1C",
  "Nemolizumab"     = "#377EB8",
  "Tralokinumab"    = "#4DAF4A",
  "Cyclosporine"    = "#FF7F00",
  "Nalbuphine ER"   = "#984EA3",
  "Dupilumab+TCS"   = "#A65628"
)

## ─── UI ─────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "PN QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("환자 프로파일", tabName="patient",  icon=icon("user")),
      menuItem("약동학 (PK)",   tabName="pk",       icon=icon("pills")),
      menuItem("면역 PD",       tabName="immune",   icon=icon("dna")),
      menuItem("피부/신경",     tabName="skin",     icon=icon("skin")),
      menuItem("임상 엔드포인트",tabName="clinical", icon=icon("chart-line")),
      menuItem("시나리오 비교", tabName="compare",  icon=icon("balance-scale")),
      menuItem("바이오마커",    tabName="bio",      icon=icon("vial"))
    )
  ),
  dashboardBody(
    tabItems(

      ## ── Tab 1: 환자 프로파일 ─────────────────────────────────
      tabItem("patient",
        fluidRow(
          box(title="환자 기본 정보", width=4, status="primary",
            sliderInput("bw",   "체중 (kg)", 40, 130, 75),
            sliderInput("age",  "나이 (세)",  18, 80, 50),
            radioButtons("sex", "성별", c("여성"=1,"남성"=0), inline=TRUE),
            sliderInput("base_il31", "기저 IL-31 (pg/mL)", 5, 50, 12),
            sliderInput("base_ige",  "기저 IgE (IU/mL)",  100, 3000, 800)
          ),
          box(title="치료 선택", width=4, status="warning",
            checkboxGroupInput("selected_scen",
              "치료 시나리오 (복수 선택 가능)",
              choices = c(
                "Placebo (무치료)" = 1,
                "Dupilumab 300mg SC Q2W" = 2,
                "Nemolizumab 60mg SC Q4W" = 3,
                "Tralokinumab 300mg SC Q2W" = 4,
                "Cyclosporine 5mg/kg/d" = 5,
                "Nalbuphine ER 54mg BID" = 6,
                "Dupilumab + TCS" = 7
              ),
              selected = c(1, 2, 3)
            ),
            sliderInput("duration", "시뮬레이션 기간 (주)", 16, 104, 52, step=4)
          ),
          box(title="용량 조정", width=4, status="info",
            sliderInput("dup_dose", "Dupilumab 용량 (mg)", 100, 600, 300, step=50),
            sliderInput("nem_dose", "Nemolizumab 용량 (mg)", 10, 120, 60, step=10),
            actionButton("run_sim", "시뮬레이션 실행", class="btn-success btn-lg",
                         icon=icon("play"))
          )
        ),
        fluidRow(
          box(title="환자 질환 중증도 프로파일", width=12,
            valueBoxOutput("vbox_itch",  width=3),
            valueBoxOutput("vbox_iga",   width=3),
            valueBoxOutput("vbox_il31",  width=3),
            valueBoxOutput("vbox_eos",   width=3)
          )
        )
      ),

      ## ── Tab 2: PK ─────────────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(title="약물 혈장 농도-시간 곡선", width=12,
            plotlyOutput("plot_pk", height="450px"))
        ),
        fluidRow(
          box(title="수용체 점유율 (Receptor Occupancy)", width=6,
            plotlyOutput("plot_ro", height="350px")),
          box(title="PK 파라미터 요약", width=6,
            tableOutput("tbl_pk"))
        )
      ),

      ## ── Tab 3: 면역 PD ────────────────────────────────────────
      tabItem("immune",
        fluidRow(
          box(title="사이토카인: IL-4, IL-13, IL-31", width=6,
            plotlyOutput("plot_cytokines", height="350px")),
          box(title="Th2 세포 및 IgE", width=6,
            plotlyOutput("plot_th2_ige", height="350px"))
        ),
        fluidRow(
          box(title="호산구 및 비만세포", width=6,
            plotlyOutput("plot_mast_eos", height="350px")),
          box(title="면역 PD 수치 요약 (Week 16)", width=6,
            tableOutput("tbl_immune"))
        )
      ),

      ## ── Tab 4: 피부/신경 ──────────────────────────────────────
      tabItem("skin",
        fluidRow(
          box(title="피부 장벽 기능 (TEWL)", width=6,
            plotlyOutput("plot_tewl", height="350px")),
          box(title="진피 신경 밀도", width=6,
            plotlyOutput("plot_nerve", height="350px"))
        ),
        fluidRow(
          box(title="결절 수 변화", width=6,
            plotlyOutput("plot_nodule", height="350px")),
          box(title="중추 감작 상태", width=6,
            plotlyOutput("plot_cs", height="350px"))
        )
      ),

      ## ── Tab 5: 임상 엔드포인트 ────────────────────────────────
      tabItem("clinical",
        fluidRow(
          box(title="가려움 NRS (0–10) 시간 경과", width=6,
            plotlyOutput("plot_itch_nrs", height="380px")),
          box(title="IGA 점수 변화", width=6,
            plotlyOutput("plot_iga", height="380px"))
        ),
        fluidRow(
          box(title="DNRS (Dynamic Neuropathic Response Score)", width=6,
            plotlyOutput("plot_dnrs", height="350px")),
          box(title="Week 16 임상 결과 요약", width=6,
            tableOutput("tbl_clinical"))
        )
      ),

      ## ── Tab 6: 시나리오 비교 ──────────────────────────────────
      tabItem("compare",
        fluidRow(
          box(title="Itch NRS — 전체 시나리오 비교", width=12,
            plotlyOutput("plot_compare_itch", height="420px"))
        ),
        fluidRow(
          box(title="Week 16 가려움 NRS % 감소", width=6,
            plotlyOutput("plot_waterfall", height="350px")),
          box(title="IGA 개선 (Week 16 vs Baseline)", width=6,
            plotlyOutput("plot_iga_bar", height="350px"))
        )
      ),

      ## ── Tab 7: 바이오마커 ────────────────────────────────────
      tabItem("bio",
        fluidRow(
          box(title="IL-31 혈청 수준", width=6,
            plotlyOutput("plot_bio_il31", height="350px")),
          box(title="혈청 IgE", width=6,
            plotlyOutput("plot_bio_ige", height="350px"))
        ),
        fluidRow(
          box(title="말초 호산구 수", width=6,
            plotlyOutput("plot_bio_eos", height="350px")),
          box(title="TEWL vs 치료 반응 산점도 (Week 16)", width=6,
            plotlyOutput("plot_bio_scatter", height="350px"))
        )
      )
    )
  )
)

## ─── Server ─────────────────────────────────────────────────────
server <- function(input, output, session) {

  sim_data <- eventReactive(input$run_sim, {
    req(input$selected_scen)
    withProgress(message="시뮬레이션 실행 중...", value=0.5, {
      run_pn(
        scen_ids     = as.integer(input$selected_scen),
        duration_wk  = input$duration,
        BW           = input$bw,
        base_il31    = input$base_il31,
        base_ige     = input$base_ige,
        dupilumab_dose = input$dup_dose,
        nemolizumab_dose = input$nem_dose
      )
    })
  }, ignoreNULL = FALSE)

  # Initialize with default scenarios on load
  observe({
    if (is.null(sim_data()) || nrow(sim_data()) == 0) {
      isolate(input$run_sim)
    }
  })

  ## Value boxes
  output$vbox_itch <- renderValueBox({
    valueBox(paste0(input$base_il31, " pg/mL"), "기저 IL-31", icon=icon("thermometer-half"), color="red")
  })
  output$vbox_iga  <- renderValueBox({
    valueBox("3-4", "기저 IGA", icon=icon("star"), color="orange")
  })
  output$vbox_il31 <- renderValueBox({
    valueBox("7-8/10", "기저 Itch NRS", icon=icon("arrows-v"), color="purple")
  })
  output$vbox_eos  <- renderValueBox({
    valueBox(paste0(input$base_ige, " IU/mL"), "기저 IgE", icon=icon("vial"), color="blue")
  })

  ## ── PK plots ─────────────────────────────────────────────────
  output$plot_pk <- renderPlotly({
    df <- sim_data(); req(df)
    df_long <- df |>
      select(time_wk, scenario, DUP_C1, NEM_C1, TRAL_C1) |>
      pivot_longer(c(DUP_C1, NEM_C1, TRAL_C1),
                   names_to="Drug", values_to="Conc_nM") |>
      mutate(Drug = recode(Drug,
        DUP_C1="Dupilumab (nM)",
        NEM_C1="Nemolizumab (nM)",
        TRAL_C1="Tralokinumab (nM)"))
    p <- ggplot(df_long, aes(x=time_wk, y=Conc_nM, color=scenario,
                              linetype=Drug)) +
      geom_line(linewidth=0.8, alpha=0.85) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="혈장 농도 (nM)", color="시나리오", linetype="약물") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_ro <- renderPlotly({
    df <- sim_data(); req(df)
    df_ro <- df |>
      select(time_wk, scenario, DUP_RO_PCT, NEM_RO_PCT) |>
      pivot_longer(c(DUP_RO_PCT, NEM_RO_PCT),
                   names_to="Drug", values_to="RO") |>
      mutate(Drug = recode(Drug,
        DUP_RO_PCT="Dupilumab (IL-4Rα RO%)",
        NEM_RO_PCT="Nemolizumab (IL-31Rα RO%)"))
    p <- ggplot(df_ro, aes(x=time_wk, y=RO, color=scenario, linetype=Drug)) +
      geom_line(linewidth=0.8) +
      scale_y_continuous(limits=c(0,100)) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="수용체 점유율 (%)", color="시나리오") +
      theme_bw()
    ggplotly(p)
  })

  output$tbl_pk <- renderTable({
    data.frame(
      약물 = c("Dupilumab","Nemolizumab","Tralokinumab","Cyclosporine","Nalbuphine ER"),
      용량 = c("300mg SC Q2W","60mg SC Q4W","300mg SC Q2W","5mg/kg/d PO","54mg PO BID"),
      t½_일 = c("21","25","14","0.8","3.5"),
      F_percent = c("64","72","77","40","20"),
      표적 = c("IL-4Rα","IL-31Rα","IL-13","칼시뉴린","κ-오피오이드")
    )
  })

  ## ── Immune PD plots ──────────────────────────────────────────
  output$plot_cytokines <- renderPlotly({
    df <- sim_data(); req(df)
    df_long <- df |>
      select(time_wk, scenario, IL4, IL13, IL31) |>
      pivot_longer(c(IL4,IL13,IL31), names_to="Cytokine", values_to="pg_mL")
    p <- ggplot(df_long, aes(x=time_wk, y=pg_mL, color=scenario,
                              linetype=Cytokine)) +
      geom_line(linewidth=0.8, alpha=0.85) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="사이토카인 (pg/mL)", color="시나리오") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_th2_ige <- renderPlotly({
    df <- sim_data(); req(df)
    p1 <- ggplot(df, aes(x=time_wk, y=TH2_CELLS, color=scenario)) +
      geom_line(linewidth=0.8) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="Th2 세포 (arb)", color="시나리오") +
      theme_bw() + theme(legend.position="none")

    p2 <- ggplot(df, aes(x=time_wk, y=IGE, color=scenario)) +
      geom_line(linewidth=0.8) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="IgE (IU/mL)", color="시나리오") +
      theme_bw()

    subplot(ggplotly(p1), ggplotly(p2), nrows=1, shareX=TRUE)
  })

  output$plot_mast_eos <- renderPlotly({
    df <- sim_data(); req(df)
    df_long <- df |>
      select(time_wk, scenario, MAST_CELLS, EOS_COUNT) |>
      pivot_longer(c(MAST_CELLS, EOS_COUNT), names_to="Cell", values_to="Count") |>
      mutate(Cell = recode(Cell, MAST_CELLS="비만세포", EOS_COUNT="호산구 (cells/μL)"))
    p <- ggplot(df_long, aes(x=time_wk, y=Count, color=scenario, linetype=Cell)) +
      geom_line(linewidth=0.8) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="세포 수", color="시나리오") +
      theme_bw()
    ggplotly(p)
  })

  output$tbl_immune <- renderTable({
    df <- sim_data(); req(df)
    df |>
      filter(abs(time_wk - 16) < 0.1) |>
      group_by(scenario) |>
      slice_head(n=1) |>
      transmute(
        시나리오 = scenario,
        IL31_pgmL = round(IL31, 1),
        IL4_pgmL = round(IL4, 1),
        IL13_pgmL = round(IL13, 1),
        IgE_IUmL = round(IGE, 0),
        Eos_cellsuL = round(EOS_COUNT, 0)
      )
  })

  ## ── Skin/Neuro plots ─────────────────────────────────────────
  output$plot_tewl <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time_wk, y=TEWL, color=scenario)) +
      geom_line(linewidth=0.8) +
      geom_hline(yintercept=10, linetype="dashed", color="green4", alpha=0.6) +
      annotate("text", x=2, y=11, label="정상 TEWL (10)", size=3, color="green4") +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="TEWL (g/m²/h)", color="시나리오") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_nerve <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time_wk, y=NERVE_DEN, color=scenario)) +
      geom_line(linewidth=0.8) +
      geom_hline(yintercept=100, linetype="dashed", color="green4", alpha=0.6) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="진피 신경 밀도 (arb)", color="시나리오") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_nodule <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time_wk, y=NODULE_CNT, color=scenario)) +
      geom_line(linewidth=0.8) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="결절 수 (arb units)", color="시나리오") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_cs <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time_wk, y=ITCH_CS, color=scenario)) +
      geom_line(linewidth=0.8) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="중추 감작 지수 (arb)", color="시나리오") +
      theme_bw()
    ggplotly(p)
  })

  ## ── Clinical endpoint plots ───────────────────────────────────
  output$plot_itch_nrs <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time_wk, y=ITCH_NRS_OUT, color=scenario)) +
      geom_line(linewidth=1.0) +
      scale_y_continuous(limits=c(0,10), breaks=0:10) +
      scale_x_continuous(breaks=seq(0,104,8)) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="가려움 NRS (0–10)", color="시나리오",
           title="Peak Pruritus NRS 시간 경과") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_iga <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time_wk, y=IGA_OUT, color=scenario)) +
      geom_line(linewidth=0.9) +
      scale_y_continuous(limits=c(0,4.5), breaks=0:4) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="IGA 점수 (0–4)", color="시나리오") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_dnrs <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time_wk, y=DNRS, color=scenario)) +
      geom_line(linewidth=0.9) +
      scale_y_continuous(limits=c(0,28)) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="DNRS (0–28)", color="시나리오") +
      theme_bw()
    ggplotly(p)
  })

  output$tbl_clinical <- renderTable({
    df <- sim_data(); req(df)
    df_bl <- df |>
      filter(time_wk < 0.1) |>
      group_by(scenario) |>
      slice_head(n=1) |>
      select(scenario, BL_NRS = ITCH_NRS_OUT)

    df |>
      filter(abs(time_wk - 16) < 0.1) |>
      group_by(scenario) |>
      slice_head(n=1) |>
      left_join(df_bl, by="scenario") |>
      transmute(
        시나리오 = scenario,
        `Itch NRS W16` = round(ITCH_NRS_OUT, 1),
        `IGA W16` = round(IGA_OUT, 1),
        `DNRS W16` = round(DNRS, 1),
        `결절수 변화(%)` = round((NODULE_CNT - 20)/20*100, 1)
      )
  })

  ## ── Scenario comparison plots ─────────────────────────────────
  output$plot_compare_itch <- renderPlotly({
    scens_all <- run_pn(1:7, duration_wk=input$duration,
                        BW=input$bw,
                        base_il31=input$base_il31,
                        base_ige=input$base_ige)
    p <- ggplot(scens_all, aes(x=time_wk, y=ITCH_NRS_OUT, color=scenario)) +
      geom_line(linewidth=1.0) +
      scale_y_continuous(limits=c(0,10)) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="Itch NRS (0–10)", color="시나리오",
           title="전체 7개 시나리오 Itch NRS 비교") +
      theme_bw(base_size=13)
    ggplotly(p)
  })

  output$plot_waterfall <- renderPlotly({
    df <- sim_data(); req(df)
    df_w16 <- df |>
      filter(abs(time_wk - 16) < 0.1) |>
      group_by(scenario) |>
      slice_head(n=1) |>
      mutate(NRS_pchg = (ITCH_NRS_OUT - 7.5) / 7.5 * 100) |>
      arrange(NRS_pchg)
    p <- ggplot(df_w16, aes(x=reorder(scenario, NRS_pchg), y=NRS_pchg,
                             fill=NRS_pchg < -20)) +
      geom_col(width=0.7) +
      geom_hline(yintercept=-20, linetype="dashed") +
      geom_hline(yintercept=-50, linetype="dashed", color="red") +
      coord_flip() +
      scale_fill_manual(values=c("TRUE"="#2196F3","FALSE"="#9E9E9E"),
                        guide="none") +
      labs(x="", y="Itch NRS % 변화 (Week 16)",
           title="Waterfall — Itch NRS % 감소") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_iga_bar <- renderPlotly({
    df <- sim_data(); req(df)
    df_w16 <- df |>
      filter(abs(time_wk - 16) < 0.1) |>
      group_by(scenario) |>
      slice_head(n=1)
    p <- ggplot(df_w16, aes(x=reorder(scenario, IGA_OUT), y=IGA_OUT,
                             fill=scenario)) +
      geom_col(width=0.7) +
      geom_hline(yintercept=4.0, linetype="dashed", alpha=0.5) +
      scale_fill_manual(values=SCEN_COLORS, guide="none") +
      coord_flip() +
      labs(x="", y="IGA (Week 16)", title="IGA 점수 — Week 16 비교") +
      theme_bw()
    ggplotly(p)
  })

  ## ── Biomarker plots ───────────────────────────────────────────
  output$plot_bio_il31 <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time_wk, y=IL31, color=scenario)) +
      geom_line(linewidth=0.8) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="IL-31 (pg/mL)", color="시나리오") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_bio_ige <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time_wk, y=IGE, color=scenario)) +
      geom_line(linewidth=0.8) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="IgE (IU/mL)", color="시나리오") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_bio_eos <- renderPlotly({
    df <- sim_data(); req(df)
    p <- ggplot(df, aes(x=time_wk, y=EOS_COUNT, color=scenario)) +
      geom_line(linewidth=0.8) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="시간 (주)", y="호산구 (cells/μL)", color="시나리오") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_bio_scatter <- renderPlotly({
    df <- sim_data(); req(df)
    df_w16 <- df |>
      filter(abs(time_wk - 16) < 0.1) |>
      group_by(scenario) |>
      slice_head(n=1)
    p <- ggplot(df_w16, aes(x=TEWL, y=ITCH_NRS_OUT,
                             color=scenario, size=IGE)) +
      geom_point(alpha=0.8) +
      scale_color_manual(values=SCEN_COLORS) +
      labs(x="TEWL (g/m²/h)", y="Itch NRS (0–10)",
           size="IgE (IU/mL)",
           title="TEWL vs Itch NRS (Week 16)") +
      theme_bw()
    ggplotly(p)
  })
}

## ─── Launch ─────────────────────────────────────────────────────
shinyApp(ui, server)
