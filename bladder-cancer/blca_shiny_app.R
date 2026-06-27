# =============================================================================
# Bladder Cancer (BLCA) QSP — Shiny Dashboard
# =============================================================================
# 7 tabs: Patient Profile | PK | Immune Biomarkers | Tumor Dynamics |
#         Clinical Endpoints | Scenario Comparison | Biomarker Analytics
# =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)

# ---- mrgsolve model ----
blca_model_code <- '
$PARAM
FGFR3_flag=0 Nectin4_flag=1 PDL1_CPS=0
kabs_BCG=0.5 kelim_BCG=0.12
kelim_Cis=0.35 kelim_Gem=1.2
CL_Pembro=0.22 V1_Pembro=3.4 Q_Pembro=0.55 V2_Pembro=4.8
IC50_Pembro=0.40 hill_Pembro=1.5
CL_Atezo=0.20 V1_Atezo=3.2 Q_Atezo=0.48 V2_Atezo=5.1
IC50_Atezo=0.60 hill_Atezo=1.4
ka_Erda=0.8 F_Erda=0.70 kelim_Erda=0.058
EC50_Erda=0.12 hill_Erda=1.8
kelim_EV=0.045 EC50_EV=0.08 hill_EV=1.6
Emax_BCG_CD8=3.0 EC50_BCG_CD8=0.5
Emax_BCG_kill=0.65 EC50_BCG_kill=0.6 hill_BCG=1.5
Emax_Cis=0.55 EC50_Cis=8.0 hill_Cis=1.2
Emax_Gem=0.45 EC50_Gem=2.5 hill_Gem=1.3
Emax_IO_kill=0.30 CD8_EC50_kill=0.8
Emax_Erda_kill=0.40 Emax_EV_kill=0.42
kg_tumor=0.0025 Kmax_tumor=1e10 T0_cells=1e8 SLD0=55.0
kin_CD8=0.02 kout_CD8=0.008 kin_Treg=0.005 kout_Treg=0.003
kin_MDSC=0.004 kout_MDSC=0.005
suppTreg=0.30 suppMDSC=0.25 Treg_EC50=1.0 MDSC_EC50=1.0
kprod_IFNg=0.015 kdeg_IFNg=0.10
kprod_PDL1=0.008 kdeg_PDL1=0.04 IFNg_stim=2.5
kprod_FGFR3=0.05 kdeg_FGFR3=0.02
kprod_NMP22=0.003 kdeg_NMP22=0.015 kg_SLD=0.00015

$CMT BCG_depot BCG_eff Cis_plasm Gem_plasm
Pembro_c Pembro_p Atezo_c Atezo_p
Erda_dep Erda_plasm EnFV_c
CD8_eff Treg MDSC_cmt TumBurd
FGFR3act PDL1_lvl IFNg_cmt NMP22_cmt SLD_cmt

$MAIN
if(NEWIND<=1){
  CD8_eff_0=kin_CD8/kout_CD8; Treg_0=kin_Treg/kout_Treg;
  MDSC_cmt_0=kin_MDSC/kout_MDSC; TumBurd_0=T0_cells;
  FGFR3act_0=(FGFR3_flag>0.5)?kprod_FGFR3/kdeg_FGFR3:0.1;
  PDL1_lvl_0=kprod_PDL1/kdeg_PDL1;
  IFNg_cmt_0=kprod_IFNg*(kin_CD8/kout_CD8)/kdeg_IFNg;
  NMP22_cmt_0=kprod_NMP22*T0_cells/kdeg_NMP22; SLD_cmt_0=SLD0;
}

$ODE
dxdt_BCG_depot=-kabs_BCG*BCG_depot;
dxdt_BCG_eff=kabs_BCG*BCG_depot-kelim_BCG*BCG_eff;
dxdt_Cis_plasm=-kelim_Cis*Cis_plasm;
dxdt_Gem_plasm=-kelim_Gem*Gem_plasm;
double k10P=CL_Pembro/V1_Pembro,k12P=Q_Pembro/V1_Pembro,k21P=Q_Pembro/V2_Pembro;
dxdt_Pembro_c=-(k10P+k12P)*Pembro_c+k21P*Pembro_p;
dxdt_Pembro_p=k12P*Pembro_c-k21P*Pembro_p;
double k10A=CL_Atezo/V1_Atezo,k12A=Q_Atezo/V1_Atezo,k21A=Q_Atezo/V2_Atezo;
dxdt_Atezo_c=-(k10A+k12A)*Atezo_c+k21A*Atezo_p;
dxdt_Atezo_p=k12A*Atezo_c-k21A*Atezo_p;
dxdt_Erda_dep=-ka_Erda*Erda_dep;
dxdt_Erda_plasm=F_Erda*ka_Erda*Erda_dep-kelim_Erda*Erda_plasm;
dxdt_EnFV_c=-kelim_EV*EnFV_c;
double E_BCG_CD8=Emax_BCG_CD8*pow(BCG_eff,hill_BCG)/(pow(EC50_BCG_CD8,hill_BCG)+pow(BCG_eff,hill_BCG));
double E_BCG_kill=Emax_BCG_kill*pow(BCG_eff,hill_BCG)/(pow(EC50_BCG_kill,hill_BCG)+pow(BCG_eff,hill_BCG));
double E_Cis=Emax_Cis*pow(Cis_plasm,hill_Cis)/(pow(EC50_Cis,hill_Cis)+pow(Cis_plasm,hill_Cis));
double E_Gem=Emax_Gem*pow(Gem_plasm,hill_Gem)/(pow(EC50_Gem,hill_Gem)+pow(Gem_plasm,hill_Gem));
double pc_ugmL=Pembro_c/V1_Pembro, ac_ugmL=Atezo_c/V1_Atezo;
double PRO=pow(pc_ugmL,hill_Pembro)/(pow(IC50_Pembro,hill_Pembro)+pow(pc_ugmL,hill_Pembro));
double ARO=pow(ac_ugmL,hill_Atezo)/(pow(IC50_Atezo,hill_Atezo)+pow(ac_ugmL,hill_Atezo));
double IO_RO=(PRO>ARO)?PRO:ARO;
double E_Erda=FGFR3_flag*Emax_Erda_kill*pow(Erda_plasm,hill_Erda)/(pow(EC50_Erda,hill_Erda)+pow(Erda_plasm,hill_Erda));
double E_EV=Nectin4_flag*Emax_EV_kill*pow(EnFV_c,hill_EV)/(pow(EC50_EV,hill_EV)+pow(EnFV_c,hill_EV));
double supp=1.0+suppTreg*Treg/(Treg_EC50+Treg)+suppMDSC*MDSC_cmt/(MDSC_EC50+MDSC_cmt);
dxdt_CD8_eff=kin_CD8*(1.0+E_BCG_CD8+3.0*IO_RO)-kout_CD8*supp*CD8_eff;
double tgfb=TumBurd/(1e9+TumBurd);
dxdt_Treg=kin_Treg*(1.0+2.0*tgfb)-kout_Treg*Treg;
double tdrive=TumBurd/(5e9+TumBurd);
dxdt_MDSC_cmt=kin_MDSC*(1.0+1.5*tdrive)-kout_MDSC*MDSC_cmt;
double CTL_kill=Emax_IO_kill*CD8_eff*IO_RO/(CD8_EC50_kill+CD8_eff*IO_RO);
double BCG_CTL=E_BCG_kill*CD8_eff/(0.5+CD8_eff);
double kd=BCG_CTL+(E_Cis+E_Gem)*0.012+CTL_kill*0.01+E_Erda*0.015+E_EV*0.014;
double Tg=kg_tumor*TumBurd*(1.0-TumBurd/Kmax_tumor);
dxdt_TumBurd=Tg-kd*TumBurd;
if(TumBurd<1.0) dxdt_TumBurd=0.0;
double kinh_F=E_Erda/(0.01+E_Erda);
dxdt_FGFR3act=kprod_FGFR3*FGFR3_flag-(kdeg_FGFR3+0.05*kinh_F)*FGFR3act;
double ifng_n=IFNg_cmt/(0.1+IFNg_cmt);
dxdt_PDL1_lvl=kprod_PDL1*(1.0+IFNg_stim*ifng_n)-kdeg_PDL1*PDL1_lvl;
dxdt_IFNg_cmt=kprod_IFNg*CD8_eff*(1.0+1.5*E_BCG_CD8)-kdeg_IFNg*IFNg_cmt;
dxdt_NMP22_cmt=kprod_NMP22*TumBurd-kdeg_NMP22*NMP22_cmt;
double Tn=TumBurd/T0_cells;
dxdt_SLD_cmt=kg_SLD*SLD_cmt*(Tn-1.0);

$TABLE
double pc_ug=Pembro_c/V1_Pembro, ac_ug=Atezo_c/V1_Atezo;
double PD1_RO=pow(pc_ug,hill_Pembro)/(pow(IC50_Pembro,hill_Pembro)+pow(pc_ug,hill_Pembro));
double PDL1_RO=pow(ac_ug,hill_Atezo)/(pow(IC50_Atezo,hill_Atezo)+pow(ac_ug,hill_Atezo));
double TumorRed=100.0*(T0_cells-TumBurd)/T0_cells;
if(TumorRed<-100.0) TumorRed=-100.0;
double SLD_chg=100.0*(SLD_cmt-SLD0)/SLD0;

$CAPTURE BCG_eff Cis_plasm Gem_plasm
pc_ug ac_ug PD1_RO PDL1_RO Erda_plasm EnFV_c
CD8_eff Treg MDSC_cmt TumBurd TumorRed SLD_cmt SLD_chg
FGFR3act PDL1_lvl IFNg_cmt NMP22_cmt
'

mod <- mcode("blca_shiny", blca_model_code, quiet = TRUE)

# ---- Dosing helpers ----
make_ev_blca <- function(drug, dose_mg, wt_kg = 70) {
  if (drug == "BCG") {
    induction <- ev(cmt = "BCG_depot", amt = 81,
                    time = seq(0, by = 168, length.out = 6))
    maint <- ev(cmt = "BCG_depot", amt = 81,
                time = seq(6*168, by = 3*168, length.out = 9))
    return(c(induction, maint))
  }
  if (drug == "GC") {
    t_q3w <- seq(0, by = 21*24, length.out = 6)
    return(c(ev(cmt="Cis_plasm", amt=dose_mg*1.73, time=t_q3w),
             ev(cmt="Gem_plasm", amt=1000*1.73, time=t_q3w)))
  }
  if (drug == "Pembrolizumab") {
    return(ev(cmt="Pembro_c", amt=200, time=seq(0, by=21*24, length.out=18)))
  }
  if (drug == "Atezolizumab") {
    return(ev(cmt="Atezo_c", amt=1200, time=seq(0, by=21*24, length.out=18)))
  }
  if (drug == "Erdafitinib") {
    return(ev(cmt="Erda_dep", amt=8, time=seq(0, by=24, length.out=365)))
  }
  if (drug == "Enfortumab Vedotin") {
    times_h <- unlist(lapply(0:7, function(c) c*28*24+c(0,7*24,14*24)))
    return(ev(cmt="EnFV_c", amt=1.25*wt_kg, time=times_h))
  }
  ev(time=0, amt=0, cmt=1)
}

run_sim <- function(drug, FGFR3=0, Nectin4=1, dose_mg=70, wt_kg=70, dur_d=365) {
  dose <- make_ev_blca(drug, dose_mg, wt_kg)
  pars <- list(FGFR3_flag = ifelse(drug=="Erdafitinib",1,FGFR3),
               Nectin4_flag = Nectin4)
  mrgsim(mod, ev=dose, param=pars,
         end=dur_d*24, delta=6, obsonly=TRUE) |>
    as.data.frame() |>
    mutate(Drug=drug, time_d=time/24)
}

drug_choices <- c("Untreated", "BCG", "GC (Cis+Gem)",
                  "Pembrolizumab", "Atezolizumab",
                  "Erdafitinib", "Enfortumab Vedotin")

drug_colors <- c("Untreated"="#999999","BCG"="#E41A1C",
                 "GC (Cis+Gem)"="#377EB8","Pembrolizumab"="#4DAF4A",
                 "Atezolizumab"="#984EA3","Erdafitinib"="#FF7F00",
                 "Enfortumab Vedotin"="#A65628")

# =====================================================================
# UI
# =====================================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "방광암 QSP 대시보드 (BLCA)"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("환자 프로파일", tabName = "patient", icon = icon("user")),
      menuItem("약동학 (PK)",   tabName = "pk",      icon = icon("flask")),
      menuItem("면역 바이오마커", tabName = "immune", icon = icon("shield-alt")),
      menuItem("종양 역학",     tabName = "tumor",   icon = icon("viruses")),
      menuItem("임상 엔드포인트", tabName = "endpoints", icon = icon("chart-line")),
      menuItem("시나리오 비교", tabName = "compare", icon = icon("balance-scale")),
      menuItem("바이오마커 분석", tabName = "biomarkers", icon = icon("dna"))
    )
  ),

  dashboardBody(
    tabItems(

      # ---- Tab 1: Patient Profile ----
      tabItem("patient",
        fluidRow(
          box(title="환자 설정", status="primary", solidHeader=TRUE, width=4,
            sliderInput("age", "나이 (세)", 40, 85, 65),
            radioButtons("sex", "성별", c("남성"="M","여성"="F"), inline=TRUE),
            sliderInput("wt", "체중 (kg)", 40, 120, 70),
            radioButtons("stage", "병기", c("NMIBC"="NMIBC","MIBC/전이"="MIBC"), inline=TRUE),
            checkboxInput("fgfr3", "FGFR3 변이 양성 (erdafitinib 대상)", FALSE),
            checkboxInput("pdl1", "PD-L1 CPS ≥ 10 (pembrolizumab 우선)", FALSE),
            selectInput("drug_main", "치료 선택", drug_choices),
            sliderInput("dur", "시뮬레이션 기간 (일)", 90, 730, 365),
            actionButton("run", "시뮬레이션 실행", class="btn-success btn-block")
          ),
          box(title="환자 요약", status="info", solidHeader=TRUE, width=8,
            DTOutput("patient_summary")
          )
        )
      ),

      # ---- Tab 2: Pharmacokinetics ----
      tabItem("pk",
        fluidRow(
          box(title="약물 혈중농도", status="primary", solidHeader=TRUE, width=12,
            plotlyOutput("pk_conc", height="350px"))
        ),
        fluidRow(
          box(title="수용체 점유율 (IO)", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("pk_ro", height="280px")),
          box(title="BCG 방광벽 효과 농도", status="danger", solidHeader=TRUE, width=6,
            plotlyOutput("pk_bcg", height="280px"))
        )
      ),

      # ---- Tab 3: Immune Biomarkers ----
      tabItem("immune",
        fluidRow(
          box(title="CD8+ 세포독성 T세포 (CTL)", status="success", solidHeader=TRUE, width=6,
            plotlyOutput("imm_cd8", height="280px")),
          box(title="조절 T세포 (Treg) & MDSC", status="danger", solidHeader=TRUE, width=6,
            plotlyOutput("imm_treg_mdsc", height="280px"))
        ),
        fluidRow(
          box(title="IFN-γ 동태", status="info", solidHeader=TRUE, width=6,
            plotlyOutput("imm_ifng", height="280px")),
          box(title="PD-L1 발현 수준", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("imm_pdl1", height="280px"))
        )
      ),

      # ---- Tab 4: Tumor Dynamics ----
      tabItem("tumor",
        fluidRow(
          box(title="종양 부담 (log scale)", status="primary", solidHeader=TRUE, width=6,
            plotlyOutput("tum_burden", height="300px")),
          box(title="종양 크기 변화 — SLD (mm)", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("tum_sld", height="300px"))
        ),
        fluidRow(
          box(title="종양 감소율 (%)", status="success", solidHeader=TRUE, width=6,
            plotlyOutput("tum_red", height="280px")),
          box(title="FGFR3 신호 활성도", status="info", solidHeader=TRUE, width=6,
            plotlyOutput("tum_fgfr3", height="280px"))
        )
      ),

      # ---- Tab 5: Clinical Endpoints ----
      tabItem("endpoints",
        fluidRow(
          box(title="SLD 변화 — 워터폴 플롯 (최종)", status="primary", solidHeader=TRUE, width=6,
            plotlyOutput("ep_waterfall", height="320px")),
          box(title="임상 결과 요약표", status="info", solidHeader=TRUE, width=6,
            DTOutput("ep_table"))
        ),
        fluidRow(
          box(title="RECIST 기반 종양 반응", status="success", solidHeader=TRUE, width=12,
            plotlyOutput("ep_recist", height="280px"))
        )
      ),

      # ---- Tab 6: Scenario Comparison ----
      tabItem("compare",
        fluidRow(
          box(title="비교 시나리오 선택", status="primary", solidHeader=TRUE, width=3,
            checkboxGroupInput("drugs_comp", "치료 선택 (복수 가능)",
              choices=drug_choices, selected=c("Untreated","BCG","Pembrolizumab")),
            checkboxInput("fgfr3_comp", "FGFR3 변이 양성", FALSE),
            actionButton("run_comp", "시나리오 비교 실행", class="btn-warning btn-block")
          ),
          box(title="종양 감소율 비교", status="info", solidHeader=TRUE, width=9,
            plotlyOutput("comp_tumor", height="380px"))
        ),
        fluidRow(
          box(title="CD8+ CTL 비교", status="success", solidHeader=TRUE, width=6,
            plotlyOutput("comp_cd8", height="280px")),
          box(title="NMP22 소변 마커 비교", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("comp_nmp22", height="280px"))
        )
      ),

      # ---- Tab 7: Biomarker Analytics ----
      tabItem("biomarkers",
        fluidRow(
          box(title="NMP22 소변 바이오마커", status="primary", solidHeader=TRUE, width=6,
            plotlyOutput("bio_nmp22", height="280px")),
          box(title="PD-L1 발현 vs IFN-γ", status="warning", solidHeader=TRUE, width=6,
            plotlyOutput("bio_pdl1_ifng", height="280px"))
        ),
        fluidRow(
          box(title="바이오마커 요약 테이블", status="info", solidHeader=TRUE, width=12,
            DTOutput("bio_table"))
        )
      )

    ) # end tabItems
  ) # end dashboardBody
)

# =====================================================================
# Server
# =====================================================================
server <- function(input, output, session) {

  # ---- Reactive simulation (single drug) ----
  sim_data <- eventReactive(input$run, {
    drug_name <- switch(input$drug_main,
      "Untreated" = "Untreated",
      "BCG" = "BCG",
      "GC (Cis+Gem)" = "GC",
      "Pembrolizumab" = "Pembrolizumab",
      "Atezolizumab" = "Atezolizumab",
      "Erdafitinib" = "Erdafitinib",
      "Enfortumab Vedotin" = "Enfortumab Vedotin",
      "Untreated"
    )
    withProgress(message="시뮬레이션 실행 중...", value=0.5, {
      run_sim(drug=drug_name,
              FGFR3=as.numeric(input$fgfr3),
              Nectin4=1,
              wt_kg=input$wt,
              dur_d=input$dur)
    })
  })

  # ---- Reactive comparison simulation ----
  comp_data <- eventReactive(input$run_comp, {
    withProgress(message="시나리오 비교 중...", value=0.3, {
      bind_rows(lapply(input$drugs_comp, function(d) {
        dn <- switch(d,
          "GC (Cis+Gem)"="GC",
          "Enfortumab Vedotin"="Enfortumab Vedotin", d)
        run_sim(drug=dn, FGFR3=as.numeric(input$fgfr3_comp), wt_kg=input$wt, dur_d=input$dur) |>
          mutate(Drug=d)
      }))
    })
  })

  # -- Patient summary --
  output$patient_summary <- renderDT({
    data.frame(
      항목=c("나이","성별","체중","병기","FGFR3 변이","PD-L1 CPS","선택 치료","시뮬레이션 기간"),
      값=c(paste0(input$age,"세"), input$sex, paste0(input$wt,"kg"),
           input$stage, ifelse(input$fgfr3,"양성(+)","음성(-)"),
           ifelse(input$pdl1,"CPS≥10","CPS<10"),
           input$drug_main, paste0(input$dur,"일"))
    ) |> datatable(options=list(dom='t', pageLength=10), rownames=FALSE)
  })

  # ---- Tab 2: PK ----
  output$pk_conc <- renderPlotly({
    df <- sim_data()
    d <- input$drug_main
    if (d=="Pembrolizumab") {
      p <- ggplot(df, aes(time_d, pc_ug)) + geom_line(color="#4DAF4A",size=1) +
        labs(x="시간 (일)", y="농도 (μg/mL)", title="펨브롤리주맙 혈중농도")
    } else if (d=="Atezolizumab") {
      p <- ggplot(df, aes(time_d, ac_ug)) + geom_line(color="#984EA3",size=1) +
        labs(x="시간 (일)", y="농도 (μg/mL)", title="아테졸리주맙 혈중농도")
    } else if (d=="Erdafitinib") {
      p <- ggplot(df, aes(time_d, Erda_plasm)) + geom_line(color="#FF7F00",size=1) +
        labs(x="시간 (일)", y="농도 (μg/mL)", title="에르다피티닙 혈중농도")
    } else if (d %in% c("GC (Cis+Gem)")) {
      p <- ggplot(df) +
        geom_line(aes(time_d, Cis_plasm, color="시스플라틴"), size=1) +
        geom_line(aes(time_d, Gem_plasm, color="젬시타빈"), size=1) +
        scale_color_manual(values=c("시스플라틴"="#377EB8","젬시타빈"="#E41A1C")) +
        labs(x="시간 (일)", y="농도 (μg/mL)", title="GC 화학요법 혈중농도", color="")
    } else if (d=="BCG") {
      p <- ggplot(df, aes(time_d, BCG_eff)) + geom_line(color="#E41A1C",size=1) +
        labs(x="시간 (일)", y="BCG 효과 농도 (a.u.)", title="BCG 방광벽 효과 농도")
    } else if (d=="Enfortumab Vedotin") {
      p <- ggplot(df, aes(time_d, EnFV_c)) + geom_line(color="#A65628",size=1) +
        labs(x="시간 (일)", y="농도 (μg/mL)", title="엔포르투맙 베도틴 혈중농도")
    } else {
      p <- ggplot(df, aes(time_d, TumBurd)) + geom_line(color="#999999",size=1) +
        labs(x="시간 (일)", y="종양 세포 수", title="무치료 — 종양 세포 수")
    }
    ggplotly(p + theme_bw())
  })

  output$pk_ro <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df) +
      geom_line(aes(time_d, PD1_RO, color="PD-1 RO (Pembro)"), size=1) +
      geom_line(aes(time_d, PDL1_RO, color="PD-L1 RO (Atezo)"), size=1) +
      scale_color_manual(values=c("PD-1 RO (Pembro)"="#4DAF4A","PD-L1 RO (Atezo)"="#984EA3")) +
      labs(x="시간 (일)", y="수용체 점유율 (0–1)", color="") + theme_bw()
    ggplotly(p)
  })

  output$pk_bcg <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time_d, BCG_eff)) + geom_line(color="#E41A1C",size=1) +
      labs(x="시간 (일)", y="BCG 방광벽 농도 (a.u.)") + theme_bw()
    ggplotly(p)
  })

  # ---- Tab 3: Immune Biomarkers ----
  output$imm_cd8 <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time_d, CD8_eff)) + geom_line(color="#4DAF4A",size=1) +
      labs(x="시간 (일)", y="CD8+ CTL (상대 단위)") + theme_bw()
    ggplotly(p)
  })

  output$imm_treg_mdsc <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df) +
      geom_line(aes(time_d, Treg, color="Treg"), size=1) +
      geom_line(aes(time_d, MDSC_cmt, color="MDSC"), size=1) +
      scale_color_manual(values=c("Treg"="#E41A1C","MDSC"="#FF7F00")) +
      labs(x="시간 (일)", y="세포 수 (상대 단위)", color="") + theme_bw()
    ggplotly(p)
  })

  output$imm_ifng <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time_d, IFNg_cmt)) + geom_line(color="#9C27B0",size=1) +
      labs(x="시간 (일)", y="IFN-γ 농도 (상대 단위)") + theme_bw()
    ggplotly(p)
  })

  output$imm_pdl1 <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time_d, PDL1_lvl)) + geom_line(color="#FF7F00",size=1) +
      labs(x="시간 (일)", y="PD-L1 발현 수준 (상대 단위)") + theme_bw()
    ggplotly(p)
  })

  # ---- Tab 4: Tumor Dynamics ----
  output$tum_burden <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time_d, log10(TumBurd+1))) + geom_line(color="#E41A1C",size=1) +
      labs(x="시간 (일)", y="log₁₀ 종양 세포 수") + theme_bw()
    ggplotly(p)
  })

  output$tum_sld <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time_d, SLD_cmt)) + geom_line(color="#377EB8",size=1) +
      labs(x="시간 (일)", y="SLD (mm)") + theme_bw()
    ggplotly(p)
  })

  output$tum_red <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time_d, TumorRed)) + geom_line(color="#4DAF4A",size=1) +
      geom_hline(yintercept=30, linetype="dashed", color="blue") +
      geom_hline(yintercept=-20, linetype="dashed", color="red") +
      labs(x="시간 (일)", y="종양 감소율 (%)") + theme_bw()
    ggplotly(p)
  })

  output$tum_fgfr3 <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time_d, FGFR3act)) + geom_line(color="#FF7F00",size=1) +
      labs(x="시간 (일)", y="FGFR3 신호 활성도 (a.u.)") + theme_bw()
    ggplotly(p)
  })

  # ---- Tab 5: Endpoints ----
  output$ep_waterfall <- renderPlotly({
    df <- sim_data()
    final_sld <- df |> filter(time_d == max(time_d)) |> slice(1)
    bar_df <- data.frame(
      Drug = input$drug_main,
      SLD_chg = final_sld$SLD_chg
    )
    p <- ggplot(bar_df, aes(x=Drug, y=SLD_chg, fill=SLD_chg>0)) +
      geom_col() +
      scale_fill_manual(values=c("TRUE"="#E41A1C","FALSE"="#4DAF4A"), guide="none") +
      geom_hline(yintercept=c(20,-30), linetype="dashed") +
      labs(x="", y="SLD 변화율 (%)", title="워터폴 플롯 (RECIST 기준선)") +
      theme_bw()
    ggplotly(p)
  })

  output$ep_table <- renderDT({
    df <- sim_data()
    final <- df |> filter(time_d == max(time_d)) |> slice(1)
    data.frame(
      지표=c("SLD 변화율 (%)","종양 감소율 (%)","CD8+ CTL (최종)","NMP22 (최종)","PD-1 RO (최종)","PD-L1 RO (최종)"),
      값=round(c(final$SLD_chg, final$TumorRed, final$CD8_eff,
                 final$NMP22_cmt, final$PD1_RO, final$PDL1_RO), 3)
    ) |> datatable(options=list(dom='t'), rownames=FALSE)
  })

  output$ep_recist <- renderPlotly({
    df <- sim_data()
    recist_df <- df |>
      mutate(Response = case_when(
        SLD_chg <= -30 ~ "부분 반응 (PR)",
        SLD_chg >= 20  ~ "진행 (PD)",
        TRUE           ~ "안정 (SD)"
      ))
    p <- ggplot(recist_df, aes(time_d, SLD_chg, color=Response)) +
      geom_line(size=1) +
      scale_color_manual(values=c("부분 반응 (PR)"="#4DAF4A","진행 (PD)"="#E41A1C","안정 (SD)"="#FF7F00")) +
      geom_hline(yintercept=c(20,-30), linetype="dashed", alpha=0.5) +
      labs(x="시간 (일)", y="SLD 변화율 (%)", color="RECIST") + theme_bw()
    ggplotly(p)
  })

  # ---- Tab 6: Scenario Comparison ----
  output$comp_tumor <- renderPlotly({
    df <- comp_data()
    p <- ggplot(df, aes(time_d, TumorRed, color=Drug)) + geom_line(size=1) +
      scale_color_manual(values=drug_colors) +
      geom_hline(yintercept=30, linetype="dashed", alpha=0.5) +
      labs(x="시간 (일)", y="종양 감소율 (%)", color="치료") + theme_bw()
    ggplotly(p)
  })

  output$comp_cd8 <- renderPlotly({
    df <- comp_data()
    p <- ggplot(df, aes(time_d, CD8_eff, color=Drug)) + geom_line(size=1) +
      scale_color_manual(values=drug_colors) +
      labs(x="시간 (일)", y="CD8+ CTL", color="치료") + theme_bw()
    ggplotly(p)
  })

  output$comp_nmp22 <- renderPlotly({
    df <- comp_data()
    p <- ggplot(df, aes(time_d, NMP22_cmt, color=Drug)) + geom_line(size=1) +
      scale_color_manual(values=drug_colors) +
      labs(x="시간 (일)", y="NMP22 (상대 단위)", color="치료") + theme_bw()
    ggplotly(p)
  })

  # ---- Tab 7: Biomarker Analytics ----
  output$bio_nmp22 <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time_d, NMP22_cmt)) + geom_line(color="#E41A1C",size=1) +
      labs(x="시간 (일)", y="NMP22 (상대 단위)", title="소변 NMP22 동태") + theme_bw()
    ggplotly(p)
  })

  output$bio_pdl1_ifng <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(IFNg_cmt, PDL1_lvl, color=time_d)) +
      geom_path(size=1) + geom_point(size=0.5) +
      scale_color_viridis_c(name="시간 (일)") +
      labs(x="IFN-γ (a.u.)", y="PD-L1 발현 (a.u.)", title="PD-L1 vs IFN-γ 위상 공간") +
      theme_bw()
    ggplotly(p)
  })

  output$bio_table <- renderDT({
    df <- sim_data()
    t_pts <- c(0, 30, 90, 180, 365) |> Filter(function(x) x <= input$dur, x=_)
    bio_sum <- df |>
      filter(round(time_d) %in% t_pts) |>
      group_by(time_d) |> slice(1) |> ungroup() |>
      select(time_d, NMP22_cmt, PDL1_lvl, IFNg_cmt, FGFR3act) |>
      mutate(across(where(is.numeric), ~round(.x, 4)))
    colnames(bio_sum) <- c("시간 (일)","NMP22","PD-L1","IFN-γ","FGFR3act")
    datatable(bio_sum, options=list(dom='t',pageLength=10), rownames=FALSE)
  })

}

shinyApp(ui, server)
