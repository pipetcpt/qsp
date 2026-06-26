## ============================================================
## Asymptomatic Hyperuricemia (AHY) — Shiny Interactive Dashboard
## 무증상 고요산혈증 인터랙티브 QSP 대시보드
## 6 Tabs: 환자프로파일 | PK | 요산역학 | 심혈관/신장 | 시나리오비교 | 바이오마커
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)

## ---- Compile mrgsolve model (inline) ----
ahy_code <- '
$PARAM
k_prod_base=700, k_dietary=200, fructose_load=0, alcohol_use=0,
k_fructose=0.15, k_alcohol=0.10,
XO_Vmax=1.0, XO_Km=50, k_XO_deg=0.02, k_XO_syn=0.02,
k_renal_base=400, k_URAT1=0.70, k_OAT_secr=0.50, k_ABCG2_ren=0.15,
GFR_0=90, k_gut_ABCG2=0.33,
F_allo=0.90, ka_allo=3.0, CL_allo=18, Vd_allo=1.6,
ke_oxy=0.033, k_form_oxy=0.60, IC50_oxy=8.0, n_oxy=1.5,
F_feb=0.84, ka_feb=2.5, CL_feb=5.2, Vd_feb=1.2,
ke_feb=0.12, IC50_feb=0.001, n_feb=1.0,
F_uric=0.95, ka_uric=2.0, CL_uric=15, Vd_uric=0.8,
ke_uric=0.14, IC50_URAT1=0.05, k_uricosuric_eff=0.35,
k_nucl=0.001, k_growth=0.005, k_dissolve=0.003, SUA_sat=6.8,
k_endo_dmg=0.005, k_endo_rep=0.02, NO_0=1.0, k_NO_UA=0.008,
k_BP_NO=5.0, BP_0=85,
k_GFR_UA=0.0008, k_GFR_BP=0.0005, k_GFR_rep=0.001,
k_IL1_MSU=0.02, k_IL1_deg=0.30, k_CRP_IL1=0.50, k_CRP_deg=0.15,
k_IR_UA=0.015, k_IR_base=1.0,
k_CV_UA=0.002, k_CV_BP=0.003, k_CV_CRP=0.001, k_CV_decay=0.0001,
BW=70, ABCG2_Q141K=0, k_ABCG2_inh=0.50,
DOSE_allo=0, DOSE_feb=0, DOSE_uric=0

$INIT
UA_plasma=7.5, UA_tissue=500, XO_free=1.0,
Oxypurinol=0, Febuxostat_C=0, Uricosuric_C=0, URAT1_free=1.0,
UrinaryUA=400, MSU_depot=0, Endothelial_fn=1.0, NO_level=1.0,
BP=85, GFR=90, IL1beta=5, CRP=2.0, InsulinResist=1.5,
CV_risk_score=0, Tophus_vol=0, ABCG2_frac=1.0

$ODE
double ABCG2_fn = ABCG2_frac*(1.0-ABCG2_Q141K*k_ABCG2_inh);
double UA_prod=(k_prod_base+k_dietary+k_fructose*fructose_load*k_prod_base+k_alcohol*alcohol_use*k_prod_base)*XO_free/24.0;
double inh_oxy=pow(Oxypurinol,n_oxy)/(pow(IC50_oxy,n_oxy)+pow(Oxypurinol,n_oxy));
double inh_feb=pow(Febuxostat_C,n_feb)/(pow(IC50_feb,n_feb)+pow(Febuxostat_C,n_feb));
double XO_active=XO_free*(1.0-std::max(inh_oxy,inh_feb));
dxdt_XO_free=k_XO_syn-k_XO_deg*XO_free;
double UA_prod_eff=UA_prod*XO_active;
double GFR_ratio=GFR/GFR_0;
double inh_URAT1=(Uricosuric_C>0)?Uricosuric_C/(IC50_URAT1+Uricosuric_C)*k_uricosuric_eff:0;
double URAT1_eff=k_URAT1*(1.0-inh_URAT1);
double ABCG2_ren_eff=k_ABCG2_ren*ABCG2_fn;
double FE_UA=std::min(1.0-URAT1_eff+k_OAT_secr+ABCG2_ren_eff,0.25);
double k_renal_eff=k_renal_base*GFR_ratio*FE_UA/12.0;
double k_gut_eff=k_gut_ABCG2*ABCG2_fn*k_renal_base/24.0;
dxdt_UA_plasma=UA_prod_eff/(0.6*BW*10.0)-k_renal_eff*UA_plasma-k_gut_eff*UA_plasma-0.5*UA_plasma+0.2*(UA_tissue/(0.4*BW*10.0));
dxdt_UA_tissue=0.5*UA_plasma*(0.6*BW*10.0)-0.2*UA_tissue;
dxdt_UrinaryUA=k_renal_eff*UA_plasma*24.0-0.01*UrinaryUA;
double k_oxy_form=(DOSE_allo>0)?DOSE_allo*F_allo*k_form_oxy/(Vd_allo*BW*24.0):0;
dxdt_Oxypurinol=k_oxy_form-ke_oxy*Oxypurinol;
double k_feb_in=(DOSE_feb>0)?DOSE_feb*F_feb/(Vd_feb*BW*24.0):0;
dxdt_Febuxostat_C=k_feb_in-ke_feb*Febuxostat_C;
double k_uric_in=(DOSE_uric>0)?DOSE_uric*F_uric/(Vd_uric*BW*24.0):0;
dxdt_Uricosuric_C=k_uric_in-ke_uric*Uricosuric_C;
dxdt_URAT1_free=0.05*(1.0-URAT1_free)-0.05*inh_URAT1;
double SUA_excess=std::max(UA_plasma-SUA_sat,0.0);
dxdt_MSU_depot=k_nucl*SUA_excess*(MSU_depot<1.0?1.0:0.1)+k_growth*SUA_excess*MSU_depot-k_dissolve*std::max(SUA_sat-UA_plasma,0.0)*MSU_depot;
dxdt_Tophus_vol=std::max(MSU_depot-10.0,0.0)*0.01-0.001*Tophus_vol;
dxdt_Endothelial_fn=k_endo_rep*(1.0-Endothelial_fn)-k_endo_dmg*UA_plasma*Endothelial_fn;
dxdt_NO_level=0.05*(1.0-NO_level)-k_NO_UA*std::max(UA_plasma-5.0,0.0)*std::max(UA_plasma-5.0,0.0)/10.0;
double BP_delta_NO=k_BP_NO*(1.0-NO_level);
dxdt_BP=0.01*(BP_0+BP_delta_NO-BP);
double GFR_dmg=k_GFR_UA*std::max(UA_plasma-6.0,0.0)+k_GFR_BP*std::max(BP-90.0,0.0);
double GFR_rep=(UA_plasma<6.0)?k_GFR_rep*(GFR_0-GFR):0;
dxdt_GFR=GFR_rep-GFR_dmg*GFR;
dxdt_IL1beta=k_IL1_MSU*MSU_depot-k_IL1_deg*IL1beta;
dxdt_CRP=k_CRP_IL1*IL1beta-k_CRP_deg*CRP;
dxdt_InsulinResist=0.001*(k_IR_base+k_IR_UA*UA_plasma-InsulinResist);
dxdt_CV_risk_score=k_CV_UA*std::max(UA_plasma-6.0,0.0)+k_CV_BP*std::max(BP-90.0,0.0)+k_CV_CRP*std::max(CRP-3.0,0.0)+k_CV_decay;
dxdt_ABCG2_frac=0;

$TABLE
double SUA=UA_plasma; double eGFR=GFR; double MAP=BP;
double hsCRP=CRP; double HOMA_IR=InsulinResist;
double GoutRisk=(UA_plasma>9)?3:(UA_plasma>8)?2:(UA_plasma>7)?1.5:1.0;
double XO_inh=(1.0-XO_free*(1.0-std::max(pow(Oxypurinol,n_oxy)/(pow(IC50_oxy,n_oxy)+pow(Oxypurinol,n_oxy)),pow(Febuxostat_C,n_feb)/(pow(IC50_feb,n_feb)+pow(Febuxostat_C,n_feb)))))*100.0;

$CAPTURE SUA eGFR MAP hsCRP HOMA_IR GoutRisk XO_inh
         Oxypurinol Febuxostat_C Uricosuric_C
         MSU_depot Tophus_vol IL1beta CV_risk_score UrinaryUA
'

mod <- mcode("ahy_shiny", ahy_code, quiet = TRUE)

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "AHY QSP Dashboard — 무증상 고요산혈증"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("환자 프로파일", tabName = "patient",
               icon = icon("user-md")),
      menuItem("약물 PK", tabName = "pk",
               icon = icon("pills")),
      menuItem("요산 역학", tabName = "urate",
               icon = icon("chart-line")),
      menuItem("심혈관·신장 영향", tabName = "organ",
               icon = icon("heartbeat")),
      menuItem("시나리오 비교", tabName = "scenario",
               icon = icon("balance-scale")),
      menuItem("바이오마커 대시보드", tabName = "biomarker",
               icon = icon("vial"))
    ),

    hr(),
    h5("  시뮬레이션 파라미터", style = "color:white; padding-left:10px"),

    ## Patient characteristics
    sliderInput("SUA_init",  "기저 SUA (mg/dL):", 5.0, 12.0, 7.5, 0.5),
    sliderInput("GFR_init",  "기저 eGFR (mL/min):", 30, 120, 90, 5),
    sliderInput("BW",        "체중 (kg):", 40, 150, 70, 5),
    sliderInput("sim_days",  "시뮬레이션 기간 (일):", 90, 1095, 730, 30),

    hr(),
    h5("  생활습관 요인", style = "color:white; padding-left:10px"),
    sliderInput("fructose",  "과당 부하 (0-1):", 0, 1, 0, 0.1),
    sliderInput("alcohol",   "알코올 사용 (0-1):", 0, 1, 0, 0.1),

    hr(),
    h5("  유전적 요인", style = "color:white; padding-left:10px"),
    checkboxInput("ABCG2_var", "ABCG2 Q141K 변이 보유", FALSE),

    hr(),
    h5("  약물 치료", style = "color:white; padding-left:10px"),
    numericInput("dose_allo",  "알로퓨리놀 (mg/day):", 0, 0, 900, 50),
    numericInput("dose_feb",   "페북소스타트 (mg/day):", 0, 0, 120, 20),
    numericInput("dose_uric",  "요산배설촉진제 (mg/day):", 0, 0, 400, 50),

    actionButton("run_sim", "시뮬레이션 실행", class = "btn-primary",
                 style = "margin:10px; width:90%")
  ),

  dashboardBody(
    tabItems(

      ## ---- TAB 1: 환자 프로파일 ----
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "무증상 고요산혈증 — 질환 개요", width = 12,
              status = "primary", solidHeader = TRUE,
              HTML("
              <div style='font-size:14px'>
              <h4>병태생리 요약</h4>
              <ul>
                <li><b>고요산혈증 정의</b>: 혈청 요산 ≥7.0mg/dL (남성), ≥6.0mg/dL (여성)</li>
                <li><b>무증상 고요산혈증</b>: 통풍 발작·토파이·요산 신석증 없는 고요산혈증</li>
                <li><b>유병률</b>: 전체 인구의 15-25%; 아시아인에서 더 높음</li>
                <li><b>주요 원인</b>: XO 과활성 (과생성형) 또는 URAT1/ABCG2 기능 이상 (배설 감소형)</li>
                <li><b>치료 결정 기준</b>: SUA ≥9mg/dL + 신장질환/고혈압/심혈관 위험 → 치료 고려</li>
              </ul>
              <h4>주요 합병증 위험</h4>
              <table border='1' style='width:100%; border-collapse:collapse'>
                <tr style='background:#1E88E5;color:white'><th>SUA 수준</th><th>통풍 발작</th><th>CKD 위험</th><th>CV 위험</th></tr>
                <tr><td>7.0-8.0 mg/dL</td><td>1.5배 증가</td><td>+18%</td><td>+12%</td></tr>
                <tr style='background:#FFF9C4'><td>8.0-9.0 mg/dL</td><td>3배 증가</td><td>+27%</td><td>+21%</td></tr>
                <tr style='background:#FFCDD2'><td>&gt;9.0 mg/dL</td><td>10배 증가</td><td>+38%</td><td>+31%</td></tr>
              </table>
              </div>")
          )
        ),
        fluidRow(
          valueBoxOutput("box_SUA"),
          valueBoxOutput("box_GFR"),
          valueBoxOutput("box_flare_risk")
        ),
        fluidRow(
          box(title = "SUA 시계열 (2년)", width = 8,
              plotlyOutput("plot_SUA_profile", height = "350px")),
          box(title = "고요산혈증 관련 위험 요약", width = 4,
              tableOutput("tbl_risk_summary"))
        )
      ),

      ## ---- TAB 2: 약물 PK ----
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "약물 농도-시간 곡선 (PK)", width = 12,
              status = "info", solidHeader = TRUE,
              plotlyOutput("plot_PK", height = "400px"))
        ),
        fluidRow(
          box(title = "XO 억제율 (%)", width = 6,
              plotlyOutput("plot_XO_inh", height = "300px")),
          box(title = "약물 PK 파라미터 요약", width = 6,
              tableOutput("tbl_PK_params"))
        )
      ),

      ## ---- TAB 3: 요산 역학 ----
      tabItem(tabName = "urate",
        fluidRow(
          box(title = "혈청 요산 (SUA) 시계열", width = 6,
              plotlyOutput("plot_SUA", height = "300px")),
          box(title = "요중 요산 (Urinary UA)", width = 6,
              plotlyOutput("plot_UUA", height = "300px"))
        ),
        fluidRow(
          box(title = "MSU 결정 침착량", width = 6,
              plotlyOutput("plot_MSU", height = "300px")),
          box(title = "통풍 발작 위험 지수", width = 6,
              plotlyOutput("plot_gout_risk", height = "300px"))
        )
      ),

      ## ---- TAB 4: 심혈관·신장 영향 ----
      tabItem(tabName = "organ",
        fluidRow(
          box(title = "사구체 여과율 (eGFR) 변화", width = 6,
              plotlyOutput("plot_GFR", height = "300px")),
          box(title = "혈압 (MAP) 변화", width = 6,
              plotlyOutput("plot_BP", height = "300px"))
        ),
        fluidRow(
          box(title = "누적 심혈관 위험 점수", width = 6,
              plotlyOutput("plot_CV", height = "300px")),
          box(title = "내피 기능 및 NO 수준", width = 6,
              plotlyOutput("plot_endo", height = "300px"))
        )
      ),

      ## ---- TAB 5: 시나리오 비교 ----
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "7개 치료 시나리오 비교", width = 12,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_scenario", height = "500px"))
        ),
        fluidRow(
          box(title = "2년 후 임상 결과 비교표", width = 12,
              DTOutput("tbl_scenario_summary"))
        )
      ),

      ## ---- TAB 6: 바이오마커 대시보드 ----
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "염증 바이오마커 (IL-1β, hs-CRP)", width = 6,
              plotlyOutput("plot_inflam", height = "300px")),
          box(title = "인슐린 저항성 (HOMA-IR)", width = 6,
              plotlyOutput("plot_HOMA", height = "300px"))
        ),
        fluidRow(
          box(title = "바이오마커 상관관계 (SUA vs. CRP/GFR)", width = 6,
              plotlyOutput("plot_scatter", height = "300px")),
          box(title = "바이오마커 요약표 (2년 최종)", width = 6,
              tableOutput("tbl_biomarker"))
        )
      )
    )  # end tabItems
  )    # end dashboardBody
)

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

  ## Reactive simulation
  sim_data <- eventReactive(input$run_sim, {
    p_update <- list(
      UA_plasma    = input$SUA_init,
      GFR          = input$GFR_init,
      GFR_0        = input$GFR_init,
      BW           = input$BW,
      fructose_load = input$fructose,
      alcohol_use   = input$alcohol,
      ABCG2_Q141K  = as.integer(input$ABCG2_var),
      DOSE_allo    = input$dose_allo,
      DOSE_feb     = input$dose_feb,
      DOSE_uric    = input$dose_uric
    )
    mod_upd <- param(mod, p_update)
    init_upd <- init(mod_upd, UA_plasma = input$SUA_init, GFR = input$GFR_init)
    out <- mrgsim(init_upd, end = input$sim_days, delta = 1, digits = 4)
    as.data.frame(out)
  }, ignoreNULL = FALSE)

  ## Scenario comparison (preset 7 scenarios)
  scenario_data <- reactive({
    sc_params <- list(
      list(label = "1. 미치료 AHY", p = list()),
      list(label = "2. 알로퓨리놀 300mg", p = list(DOSE_allo = 300)),
      list(label = "3. 페북소스타트 80mg", p = list(DOSE_feb = 80)),
      list(label = "4. 알로+요산배설촉진제", p = list(DOSE_allo = 300, DOSE_uric = 200)),
      list(label = "5. 고과당+알코올", p = list(fructose_load = 1.0, alcohol_use = 0.5)),
      list(label = "6. ABCG2Q141K+페북소120", p = list(ABCG2_Q141K = 1, DOSE_feb = 120)),
      list(label = "7. 적극치료 SUA<6", p = list(DOSE_allo = 600))
    )
    purrr::map_df(sc_params, function(sc) {
      mod_up <- param(mod, sc$p)
      out <- mrgsim(mod_up, end = 730, delta = 1, digits = 4)
      df <- as.data.frame(out)
      df$scenario <- sc$label
      df
    })
  })

  ## ---- Value boxes ----
  output$box_SUA <- renderValueBox({
    df <- sim_data()
    val <- round(tail(df$SUA, 1), 1)
    col <- if (val < 6) "green" else if (val < 7) "yellow" else "red"
    valueBox(paste0(val, " mg/dL"), "최종 SUA", icon = icon("tint"), color = col)
  })
  output$box_GFR <- renderValueBox({
    df <- sim_data()
    val <- round(tail(df$eGFR, 1), 0)
    col <- if (val >= 60) "green" else if (val >= 30) "yellow" else "red"
    valueBox(paste0(val, " mL/min"), "최종 eGFR", icon = icon("filter"), color = col)
  })
  output$box_flare_risk <- renderValueBox({
    df <- sim_data()
    sval <- tail(df$SUA, 1)
    risk <- if (sval > 9) "매우 높음 (×10)" else if (sval > 8) "높음 (×3)" else
            if (sval > 7) "중등도 (×1.5)" else "낮음"
    col  <- if (sval > 9) "red" else if (sval > 8) "orange" else
            if (sval > 7) "yellow" else "green"
    valueBox(risk, "통풍 발작 위험", icon = icon("exclamation-triangle"), color = col)
  })

  ## ---- Tab 1: SUA profile ----
  output$plot_SUA_profile <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, SUA)) +
      geom_line(color = "#E53935", size = 1) +
      geom_hline(yintercept = 6.8, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 6.0, linetype = "dotted", color = "blue") +
      labs(x = "시간 (일)", y = "SUA (mg/dL)", title = "혈청 요산 시계열") +
      theme_bw()
    ggplotly(p)
  })
  output$tbl_risk_summary <- renderTable({
    df <- sim_data()
    final <- tail(df, 1)
    data.frame(
      지표 = c("혈청 요산", "eGFR", "MAP", "hs-CRP", "HOMA-IR", "MSU 결정", "CV 위험 점수"),
      값   = round(c(final$SUA, final$eGFR, final$MAP,
                     final$hsCRP, final$HOMA_IR,
                     final$MSU_depot, final$CV_risk_score), 2),
      단위 = c("mg/dL", "mL/min", "mmHg", "mg/L", "", "mg", "rel.")
    )
  })

  ## ---- Tab 2: PK ----
  output$plot_PK <- renderPlotly({
    df <- sim_data()
    pk_long <- df %>%
      select(time, Oxypurinol, Febuxostat_C, Uricosuric_C) %>%
      pivot_longer(-time, names_to = "Drug", values_to = "Conc")
    p <- ggplot(pk_long, aes(time, Conc, color = Drug)) +
      geom_line(size = 1) +
      scale_color_manual(values = c(Oxypurinol = "#1E88E5",
                                    Febuxostat_C = "#43A047",
                                    Uricosuric_C = "#8E24AA")) +
      labs(x = "시간 (일)", y = "혈중 농도 (mg/L)", title = "약물 PK 시계열") +
      theme_bw()
    ggplotly(p)
  })
  output$plot_XO_inh <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, XO_inh)) +
      geom_line(color = "#F4511E", size = 1) +
      geom_hline(yintercept = 50, linetype = "dashed") +
      labs(x = "시간 (일)", y = "XO 억제율 (%)", title = "잔틴 산화효소 억제율") +
      theme_bw()
    ggplotly(p)
  })
  output$tbl_PK_params <- renderTable({
    data.frame(
      약물 = c("알로퓨리놀", "옥시퓨리놀", "페북소스타트", "레시누라드"),
      F    = c("90%", "-", "84%", "≈100%"),
      t_half = c("1-2h", "18-30h", "5-8h", "5h"),
      IC50_XO = c("(전구체)", "8.0mg/L", "0.001mg/L", "-"),
      IC50_URAT1 = c("-", "-", "-", "0.05mg/L")
    )
  })

  ## ---- Tab 3: 요산 역학 ----
  output$plot_SUA <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, SUA)) +
      geom_line(color = "#E53935", size = 1) +
      geom_hline(yintercept = 6.8, linetype = "dashed", color = "red") +
      labs(x = "일", y = "mg/dL", title = "혈청 요산") + theme_bw()
    ggplotly(p)
  })
  output$plot_UUA <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, UrinaryUA)) +
      geom_line(color = "#1E88E5", size = 1) +
      geom_hline(yintercept = 800, linetype = "dashed", color = "orange",
                 alpha = 0.7) +
      labs(x = "일", y = "mg/day", title = "요중 요산") + theme_bw()
    ggplotly(p)
  })
  output$plot_MSU <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, MSU_depot)) +
      geom_area(fill = "#EF9A9A", alpha = 0.5) +
      geom_line(color = "#B71C1C", size = 1) +
      labs(x = "일", y = "MSU (mg)", title = "MSU 결정 침착량") + theme_bw()
    ggplotly(p)
  })
  output$plot_gout_risk <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, GoutRisk)) +
      geom_line(color = "#F57F17", size = 1) +
      geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
      labs(x = "일", y = "위험 지수", title = "통풍 발작 위험 (상대적)") + theme_bw()
    ggplotly(p)
  })

  ## ---- Tab 4: 심혈관·신장 ----
  output$plot_GFR <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, eGFR)) +
      geom_line(color = "#0288D1", size = 1) +
      geom_hline(yintercept = 60, linetype = "dashed", color = "orange") +
      labs(x = "일", y = "mL/min/1.73m²", title = "eGFR 변화") + theme_bw()
    ggplotly(p)
  })
  output$plot_BP <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, MAP)) +
      geom_line(color = "#D32F2F", size = 1) +
      geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
      labs(x = "일", y = "MAP (mmHg)", title = "평균 동맥압 변화") + theme_bw()
    ggplotly(p)
  })
  output$plot_CV <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, CV_risk_score)) +
      geom_line(color = "#C62828", size = 1) +
      geom_area(fill = "#FFCDD2", alpha = 0.4) +
      labs(x = "일", y = "누적 점수", title = "누적 심혈관 위험") + theme_bw()
    ggplotly(p)
  })
  output$plot_endo <- renderPlotly({
    df <- sim_data() %>% mutate(NO = NO_level)
    p <- ggplot(df) +
      geom_line(aes(time, Endothelial_fn, color = "내피 기능"), size = 1) +
      geom_line(aes(time, NO, color = "NO 수준"), size = 1) +
      scale_color_manual(values = c("내피 기능" = "#00897B", "NO 수준" = "#1565C0")) +
      labs(x = "일", y = "정규화 값 (0-1)", title = "내피 기능 / NO", color = "") +
      theme_bw()
    ggplotly(p)
  })

  ## ---- Tab 5: 시나리오 비교 ----
  output$plot_scenario <- renderPlotly({
    df <- scenario_data()
    p <- ggplot(df, aes(time, SUA, color = scenario)) +
      geom_line(size = 0.9) +
      geom_hline(yintercept = 6.8, linetype = "dashed", color = "red", alpha = 0.7) +
      geom_hline(yintercept = 6.0, linetype = "dotted", color = "blue", alpha = 0.7) +
      labs(x = "시간 (일)", y = "SUA (mg/dL)",
           title = "7개 시나리오 — SUA 비교", color = "시나리오") +
      theme_bw(base_size = 11)
    ggplotly(p)
  })
  output$tbl_scenario_summary <- renderDT({
    df <- scenario_data() %>%
      filter(time == 730) %>%
      select(scenario, SUA, eGFR, MAP, hsCRP, HOMA_IR,
             MSU_depot, CV_risk_score, GoutRisk) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
    colnames(df) <- c("시나리오", "SUA(mg/dL)", "eGFR", "MAP(mmHg)",
                      "hs-CRP", "HOMA-IR", "MSU(mg)", "CV위험", "통풍위험")
    datatable(df, options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE) %>%
      DT::formatStyle("SUA(mg/dL)",
                      backgroundColor = DT::styleInterval(c(6.0, 7.0, 9.0),
                        c("#C8E6C9", "#FFFFFF", "#FFF9C4", "#FFCDD2")))
  })

  ## ---- Tab 6: 바이오마커 ----
  output$plot_inflam <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df) +
      geom_line(aes(time, IL1beta, color = "IL-1β (pg/mL)"), size = 1) +
      geom_line(aes(time, hsCRP, color = "hs-CRP (mg/L)"), size = 1) +
      scale_color_manual(values = c("IL-1β (pg/mL)" = "#8E24AA",
                                    "hs-CRP (mg/L)" = "#E53935")) +
      labs(x = "일", y = "농도", title = "염증 바이오마커", color = "") +
      theme_bw()
    ggplotly(p)
  })
  output$plot_HOMA <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(time, HOMA_IR)) +
      geom_line(color = "#F57F17", size = 1) +
      geom_hline(yintercept = 2.5, linetype = "dashed", color = "red") +
      labs(x = "일", y = "HOMA-IR", title = "인슐린 저항성 (HOMA-IR)") + theme_bw()
    ggplotly(p)
  })
  output$plot_scatter <- renderPlotly({
    df <- sim_data() %>% filter(time %% 30 == 0)
    p <- ggplot(df, aes(SUA, hsCRP, color = time)) +
      geom_point(size = 2) +
      geom_smooth(method = "loess", se = FALSE, color = "red") +
      scale_color_viridis_c() +
      labs(x = "SUA (mg/dL)", y = "hs-CRP (mg/L)",
           title = "SUA vs. hs-CRP 상관관계", color = "시간(일)") + theme_bw()
    ggplotly(p)
  })
  output$tbl_biomarker <- renderTable({
    df <- sim_data()
    final <- tail(df, 1)
    data.frame(
      바이오마커 = c("SUA", "hs-CRP", "IL-1β", "HOMA-IR",
                     "MSU 결정", "토파이 부피", "CV 위험", "eGFR"),
      값         = round(c(final$SUA, final$hsCRP, final$IL1beta,
                           final$HOMA_IR, final$MSU_depot,
                           final$Tophus_vol, final$CV_risk_score, final$eGFR), 2),
      단위       = c("mg/dL", "mg/L", "pg/mL", "",
                     "mg", "mm³", "rel.", "mL/min"),
      정상기준   = c("<7.0", "<3.0", "<10", "<2.5",
                     "0", "0", "-", ">60")
    )
  })
}

## Run
shinyApp(ui, server)
