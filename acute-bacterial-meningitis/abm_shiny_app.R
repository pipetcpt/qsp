## =============================================================================
##  abm_shiny_app.R
##  Acute Bacterial Meningitis (pneumococcal) — QSP 대시보드
##  급성 세균성 수막염 QSP 인터랙티브 대시보드 (12 탭)
##
##  Run with:
##      shiny::runApp("abm_shiny_app.R")
##  (abm_mrgsolve_model.R 가 같은 디렉토리에 있어야 한다)
##
##  이 앱이 보여주려는 것은 하나다: **손상 flux 는 곱이고, 그 곱은 첫 항생제
##  투여에서 최대가 된다.**  그래서 "덱사메타손 시점" 슬라이더가 이 앱의
##  중심 조작기이고, 다른 모든 탭은 그 조작이 어디로 전파되는지를 보여준다.
##
##      injury flux = k_kill(t) · N(t) · Y_lysis · (1 − E_dex(t))
##
##  두 번째로 보여주려는 것: 혈액-CSF 장벽은 양방향 문이다.  "장벽과 침투" 탭
##  에서 덱사메타손을 켜면 염증(그리고 부종·ICP)이 내려가는 동시에 반코마이신
##  CSF 농도도 내려간다 — 감수성 균주에서는 무해하고 내성 균주에서는 해롭다.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

source("abm_mrgsolve_model.R")

theme_abm <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

## 문헌 참조 범위 (진단 패널 탭의 음영 띠)
ref_bands <- list(
  PMN      = c(1000, 5000),
  PROT_CSF = c(100, 500),
  GLC_CSF  = c(0, 40),
  LAC_CSF  = c(3.5, 12),
  QALB     = c(30, 100),
  TNF      = c(100, 1000),
  MMP9     = c(100, 1000),
  ICP      = c(20, 40)
)

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("급성 세균성 수막염 (폐렴균) — QSP 대시보드 · 63-ODE 모델"),
  tags$p(style = "color:#444;",
         HTML("손상 flux = <b>k_kill(t) · N(t) · Y_lysis · (1 − E_dex(t))</b> — ",
              "이 곱은 <b>첫 항생제 투여 직후</b>에 최대가 된다. ",
              "그래서 스테로이드는 그 봉우리보다 <b>먼저</b> 켜져 있어야 한다.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("발현 시점의 환자"),
      sliderInput("logN0", "CSF 균 밀도 log₁₀(CFU/mL)", 4, 9, 7, step = 0.25),
      sliderInput("logNb0", "균혈증 log₁₀(CFU/mL)", 0, 6, 3, step = 0.5),
      sliderInput("host", "숙주방어 지수", 0.3, 1.0, 1.0, step = 0.05),
      sliderInput("wt", "체중 (kg)", 40, 120, 70, step = 5),

      hr(), h4("항생제"),
      checkboxInput("use_cef", "세프트리악손 2 g q12h", TRUE),
      checkboxInput("cef_ci", "  └ 지속주입으로 (4 g/일)", FALSE),
      sliderInput("cef_delay", "항생제 투여 지연 (h)", 0, 24, 0, step = 0.5),
      selectInput("mic_cef", "폐렴균 세프트리악손 MIC (mg/L)",
                  c("0.03 (감수성)" = 0.03, "0.5 (중간)" = 0.5,
                    "2 (내성)" = 2, "4 (고내성)" = 4), selected = 0.03),
      checkboxInput("use_van", "반코마이신 15 mg/kg q6h", FALSE),
      checkboxInput("use_rif", "리팜핀 600 mg q12h", FALSE),
      sliderInput("rif_lead", "  └ 리팜핀 선행 시간 (h)", 0, 6, 0, step = 0.5),

      hr(), h4("덱사메타손 — 이 앱의 중심 조작기"),
      checkboxInput("use_dex", "덱사메타손 0.15 mg/kg q6h", TRUE),
      sliderInput("dex_time", "투여 시점 (항생제 기준, h)", -2, 24, -0.33, step = 0.33),
      sliderInput("dex_days", "투여 기간 (일)", 1, 4, 4, step = 1),

      hr(), h4("보조요법"),
      checkboxInput("use_mann", "만니톨 0.5 g/kg q6h", FALSE),
      checkboxInput("use_gly", "글리세롤 1.5 g/kg q6h", FALSE),
      sliderInput("drain", "EVD 배액 (mL/h)", 0, 24, 0, step = 2),
      sliderInput("antipyr", "해열제 효과", 0, 1, 0, step = 0.1),
      sliderInput("anticonv", "항경련제 효과", 0, 1, 0, step = 0.1),

      hr(),
      sliderInput("tend", "시뮬레이션 기간 (h)", 48, 336, 336, step = 24),
      actionButton("go", "다시 계산", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ------------------------------------------------------------------ 1
        tabPanel("1 · 환자 프로파일",
                 h4("발현 시점의 상태와 즉시 예측되는 결과"),
                 fluidRow(column(6, tableOutput("tbl_profile")),
                          column(6, tableOutput("tbl_endpoints"))),
                 hr(),
                 h4("치료 타임라인"),
                 plotOutput("p_timeline", height = "220px"),
                 helpText("t = 0 은 감염 시작이 아니라 병원 도착(발현)이다. ",
                          "증상은 이미 12–36 시간 진행되어 있으므로 초기조건은 ",
                          "폐렴균 수막염의 보고된 CSF 패널에 맞춰져 있다.")),

        ## ------------------------------------------------------------------ 2
        tabPanel("2 · 항생제 PK",
                 h4("혈장과 CSF — 전달 단계가 치료지수를 결정한다"),
                 plotOutput("p_pk", height = "440px"),
                 fluidRow(column(6, h5("CSF 노출 요약"), tableOutput("tbl_pk")),
                          column(6, h5("침투율 (CSF / 총 혈장농도)"),
                                 plotOutput("p_pen", height = "220px"))),
                 helpText("세프트리악손의 단백결합은 포화성이다 — 유리분율이 ",
                          "30 mg/L 에서 0.076, 250 mg/L 에서 0.167 로 오른다. ",
                          "CSF 로 넘어가는 것은 유리약물뿐이다.")),

        ## ------------------------------------------------------------------ 3
        tabPanel("3 · 세균 동역학",
                 h4("CSF 유리균 · 부착 아집단 · 균혈증"),
                 plotOutput("p_bact", height = "340px"),
                 plotOutput("p_kill", height = "240px"),
                 tableOutput("tbl_ster"),
                 helpText("멸균 판정은 CSF 유리균 N_c < 10 CFU/mL (배양 검출한계) ",
                          "기준이다. 부착·격리 아집단은 살균효과의 35 % 만 받으므로 ",
                          "언제나 늦게 사라진다 — 완전 치료기간이 필요한 이유가 ",
                          "이 두 곡선의 간극이다.")),

        ## ------------------------------------------------------------------ 4
        tabPanel("4 · 손상 flux (곱)",
                 h4("k_kill × N × Y × (1 − E_dex) — 이 모델의 중심 주장"),
                 plotOutput("p_product", height = "420px"),
                 fluidRow(column(6, plotOutput("p_cargo", height = "260px")),
                          column(6, plotOutput("p_burst_zoom", height = "260px"))),
                 helpText("k_kill 은 첫 투여에서 0 → E_max 로 점프하고 N 은 정확히 ",
                          "그때 최대다. 따라서 손상 flux 의 봉우리는 치료의 첫 몇 ",
                          "시간에 있다. 항생제를 늦추면 N 이 지수적으로 커져 있어 ",
                          "봉우리가 커진다 — 지연이 나쁜 이유가 '균이 많아서'가 ",
                          "아니라 '죽일 것이 많아서'다.")),

        ## ------------------------------------------------------------------ 5
        tabPanel("5 · 염증",
                 h4("CSF 사이토카인 · 호중구 · 단백분해 · 산화 부하"),
                 plotOutput("p_cyto", height = "420px"),
                 plotOutput("p_pmn", height = "260px"),
                 helpText("음영 띠는 폐렴균 수막염에서 보고된 CSF 범위다. ",
                          "덱사메타손 효과구획(IDEX)은 비대칭이다 — 켜짐 t½ 1.5 h, ",
                          "꺼짐 t½ 17 h. 그 비대칭이 투여 시점을 중요하게 만든다.")),

        ## ------------------------------------------------------------------ 6
        tabPanel("6 · 장벽과 침투 (두 개의 문)",
                 h4("염증이 열고 스테로이드가 닫는 같은 문"),
                 plotOutput("p_barrier", height = "400px"),
                 h5("부호가 갈리는 곳"),
                 tableOutput("tbl_signflip"),
                 helpText("장벽 Pb 가 열리면 친수성 항생제가 들어오고 동시에 ",
                          "알부민·호중구·물이 들어온다. 덱사메타손은 둘 다 줄인다. ",
                          "세프트리악손은 C/MIC 여유가 200배라 순이득이지만, ",
                          "세팔로스포린 내성균에 대한 반코마이신은 여유가 없어 ",
                          "같은 조작이 멸균을 2배 늦춘다.")),

        ## ------------------------------------------------------------------ 7
        tabPanel("7 · CSF 진단 패널",
                 h4("실제로 측정되는 값들 — 문헌 범위와 대조"),
                 plotOutput("p_panel", height = "520px"),
                 tableOutput("tbl_panel"),
                 helpText("포도당 경쟁: 호중구 4,000/µL 는 20 mg/dL/h 를 먹고 ",
                          "세균 10⁷ CFU/mL 는 2 mg/dL/h 를 먹는다. 열 배 차이다. ",
                          "그래서 CSF 가 배양음성이 된 뒤에도 저포도당이 남는다 — ",
                          "포도당은 조기 멸균 지표가 될 수 없다.")),

        ## ------------------------------------------------------------------ 8
        tabPanel("8 · 두개내압과 관류",
                 h4("CSF 역학 · 부종 · CPP · 허혈"),
                 plotOutput("p_icp", height = "420px"),
                 plotOutput("p_perf", height = "260px"),
                 helpText("항정상태에서 ICP = P_ss + Q_f · R_out 로 환원되는 ",
                          "지수적 순응도 모델(Marmarou PVI 25 mL)이다. ",
                          "삼투요법의 효력은 1/Pb 에 비례한다 — 부종을 만든 ",
                          "장벽 손상이 삼투제도 새게 한다.")),

        ## ------------------------------------------------------------------ 9
        tabPanel("9 · 임상 엔드포인트",
                 h4("사망 · 청력 · 인지 · 국소결손"),
                 fluidRow(column(5, tableOutput("tbl_out")),
                          column(7, plotOutput("p_hazard", height = "300px"))),
                 plotOutput("p_injury", height = "280px"),
                 helpText("위험함수는 급성 생리이탈만 적분하고 구조적 손상은 ",
                          "종료시점에 1회 평가한다. 영구손상을 시간에 곱하면 ",
                          "회복한 환자도 무한한 위험을 적립하게 된다.")),

        ## ------------------------------------------------------------------ 10
        tabPanel("10 · 덱사메타손 타이밍 스윕",
                 h4("방패를 언제 올리는가"),
                 actionButton("run_sweep", "스윕 실행", class = "btn-warning"),
                 plotOutput("p_sweep", height = "420px"),
                 tableOutput("tbl_sweep"),
                 helpText("사이토카인 폭발은 첫 항생제 투여 후 1–4 시간에 정점을 ",
                          "지난다. 전사효과가 켜지는 데 t½ 1.5 시간이 걸리므로, ",
                          "4 시간 늦은 투여는 폭발을 놓친다. 이것이 '항생제 전 ",
                          "또는 함께'라는 권고의 기전적 내용이다.")),

        ## ------------------------------------------------------------------ 11
        tabPanel("11 · 시나리오 비교",
                 h4("26개 치료 시나리오"),
                 actionButton("run_all", "전체 시나리오 실행 (수 분 소요)",
                              class = "btn-warning"),
                 tableOutput("tbl_all"),
                 hr(),
                 h5("가상 코호트 10명 × (DEX 유/무)"),
                 actionButton("run_coh", "코호트 실행", class = "btn-warning"),
                 tableOutput("tbl_coh"),
                 helpText("대조: de Gans & van de Beek NEJM 2002 폐렴균 아군 — ",
                          "사망 34 % → 14 %, 불량결과 52 % → 26 %.")),

        ## ------------------------------------------------------------------ 12
        tabPanel("12 · 방법론과 검증",
                 h4("이 모델이 어떻게 검증되었는가"),
                 verbatimTextOutput("txt_checks"),
                 hr(),
                 h4("수치검증이 잡아낸 결함"),
                 htmlOutput("txt_bugs"))
      )
    )
  )
)

## =============================================================================
##  Server
## =============================================================================
server <- function(input, output, session) {

  dosing <- reactive({
    wt <- input$wt
    e <- NULL
    add <- function(a, b) if (is.null(a)) b else a + b
    if (input$use_cef) {
      e <- add(e, if (input$cef_ci) ev_cef_ci(input$cef_delay)
                  else ev_cef(input$cef_delay))
    }
    if (input$use_van) {
      e <- add(e, ev(time = input$cef_delay, amt = 15 * wt, cmt = 4,
                     rate = 15 * wt, ii = 6, addl = 55))
    }
    if (input$use_rif) {
      e <- add(e, ev_rif(input$cef_delay - input$rif_lead))
    }
    if (input$use_dex) {
      e <- add(e, ev(time = input$cef_delay + input$dex_time,
                     amt = 0.15 * wt, cmt = 9, rate = 0.15 * wt / 0.25,
                     ii = 6, addl = 4 * input$dex_days - 1))
    }
    if (input$use_mann) {
      e <- add(e, ev(time = 1, amt = 0.5 * wt * 1000, cmt = 12,
                     rate = 0.5 * wt * 1000 / 0.5, ii = 6, addl = 11))
    }
    if (input$use_gly) {
      e <- add(e, ev(time = 0, amt = 1.5 * wt * 1000, cmt = 13, ii = 6, addl = 15))
    }
    e
  })

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    withProgress(message = "적분 중...", {
      run_abm(dosing(), end = input$tend, delta = 0.1,
              N0 = 10^input$logN0, NB0 = 10^input$logNb0,
              HOST_DEF = input$host,
              MU_SCALE = 1 + 0.5 * (1 - input$host),
              MIC_cef = as.numeric(input$mic_cef),
              ANTIPYR = input$antipyr, ANTICONV = input$anticonv,
              CSF_DRAIN = input$drain)
    })
  })

  D <- reactive(as.data.frame(sim()))
  E <- reactive(endpoints_abm(sim()))

  long <- function(d, vars) {
    d %>% select(time, all_of(vars)) %>%
      pivot_longer(-time, names_to = "var", values_to = "val")
  }

  band <- function(v) {
    b <- ref_bands[[v]]
    if (is.null(b)) return(NULL)
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = b[1], ymax = b[2],
             fill = "steelblue", alpha = 0.10)
  }

  ## ---------------------------------------------------------------- tab 1
  output$tbl_profile <- renderTable({
    d <- D()[1, ]
    data.frame(
      항목 = c("CSF 균 밀도 (CFU/mL)", "균혈증 (CFU/mL)", "CSF 백혈구 (/µL)",
               "CSF 단백 (mg/dL)", "CSF 포도당 (mg/dL)", "CSF 포도당비",
               "CSF 락테이트 (mmol/L)", "알부민 지수 Q_alb", "장벽 투과성 Pb",
               "ICP (mmHg)", "CPP (mmHg)", "체온 (°C)", "SOFA"),
      발현시 = c(sprintf("%.1e", 10^input$logN0), sprintf("%.1e", 10^input$logNb0),
                 sprintf("%.0f", d$PMN), sprintf("%.0f", d$PROT_CSF),
                 sprintf("%.1f", d$GLC_CSF), sprintf("%.2f", d$GLCR),
                 sprintf("%.1f", d$LAC_CSF), sprintf("%.0f", d$QALB),
                 sprintf("%.1f", d$PB), sprintf("%.1f", d$ICP),
                 sprintf("%.1f", d$CPPo), sprintf("%.1f", d$TEMP),
                 sprintf("%.1f", d$SOFA)))
  }, striped = TRUE)

  output$tbl_endpoints <- renderTable({
    e <- E()
    data.frame(
      엔드포인트 = c("CSF 멸균 (h)", "부착 아집단 제거 (h)", "최고 ICP (mmHg)",
                     "최저 CPP (mmHg)", "청력역치 이동 (dB)", "인지 z",
                     "사망확률 (%)", "급성 위험 적분"),
      값 = c(ifelse(is.na(e$sterile_h), "미달성", sprintf("%.1f", e$sterile_h)),
             ifelse(is.na(e$adh_clear_h), "미달성", sprintf("%.1f", e$adh_clear_h)),
             sprintf("%.1f", e$peak_ICP), sprintf("%.1f", e$min_CPP),
             sprintf("%.1f", e$hear_dB), sprintf("%.2f", e$cog_z),
             sprintf("%.1f", 100 * e$p_death), sprintf("%.3f", e$haz_acute)))
  }, striped = TRUE)

  output$p_timeline <- renderPlot({
    ev <- dosing()
    if (is.null(ev)) return(NULL)
    df <- as.data.frame(ev@data)
    nm <- c("1" = "세프트리악손", "4" = "반코마이신", "7" = "리팜핀",
            "9" = "덱사메타손", "12" = "만니톨", "13" = "글리세롤")
    df$drug <- nm[as.character(df$cmt)]
    df <- df %>% mutate(n = ifelse(is.na(addl), 1, addl + 1),
                        ii2 = ifelse(is.na(ii), 0, ii))
    ex <- do.call(rbind, lapply(seq_len(nrow(df)), function(i)
      data.frame(drug = df$drug[i],
                 time = df$time[i] + (0:(df$n[i] - 1)) * df$ii2[i])))
    ex <- ex[ex$time <= min(input$tend, 96), ]
    ggplot(ex, aes(time, drug, colour = drug)) +
      geom_point(size = 3, shape = 25, fill = "white") +
      labs(x = "시간 (h)", y = NULL, title = "첫 96시간의 투여 시점") +
      theme_abm + theme(legend.position = "none")
  })

  ## ---------------------------------------------------------------- tab 2
  output$p_pk <- renderPlot({
    d <- D()
    d$CEF_PL <- d$CEF_C / 8
    d$VAN_PL <- d$VAN_C / 20
    d$RIF_PL <- d$RIF_C / 50
    d$DEX_PL <- d$DEX_C / 70
    p1 <- long(d, c("CEF_PL", "VAN_PL", "RIF_PL", "DEX_PL")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.7) +
      scale_y_log10() + labs(x = NULL, y = "혈장 총농도 (mg/L)", colour = NULL,
                             title = "혈장") + theme_abm
    p2 <- long(d, c("CEF_CSF", "VAN_CSF", "RIF_CSF", "DEX_CSF")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.7) +
      scale_y_log10() + labs(x = "시간 (h)", y = "CSF 농도 (mg/L)", colour = NULL,
                             title = "CSF — 여기가 세균이 있는 곳") + theme_abm
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  })

  output$p_pen <- renderPlot({
    long(D(), c("PEN_CEF", "PEN_VAN")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.8) +
      labs(x = "시간 (h)", y = "CSF / 총 혈장농도", colour = NULL) +
      theme_abm
  })

  output$tbl_pk <- renderTable({
    e <- E()
    mic <- as.numeric(input$mic_cef)
    data.frame(
      지표 = c("CEF CSF AUC (mg·h/L)", "VAN CSF AUC (mg·h/L)",
               sprintf("CEF T>4×MIC (h)  [4×MIC = %.2f]", 4 * mic),
               "VAN T>4×MIC (h)  [4 mg/L]", "최고 Pb"),
      값 = c(sprintf("%.0f", D()$AUC_CEF[nrow(D())]),
             sprintf("%.0f", e$AUC_van), sprintf("%.1f", e$T_cef_4MIC),
             sprintf("%.1f", e$T_van_4MIC), sprintf("%.1f", e$peak_Pb)))
  }, striped = TRUE)

  ## ---------------------------------------------------------------- tab 3
  output$p_bact <- renderPlot({
    d <- D()
    d$log_NC <- log10(pmax(d$NC, 1e-2))
    d$log_NADH <- log10(pmax(d$NADH, 1e-2))
    d$log_NB <- log10(pmax(d$NB, 1e-2))
    long(d, c("log_NC", "log_NADH", "log_NB")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
      annotate("text", x = Inf, y = 1.3, hjust = 1.05, size = 3,
               label = "배양 검출한계 10 CFU/mL") +
      labs(x = "시간 (h)", y = "log₁₀ CFU/mL", colour = NULL) + theme_abm
  })

  output$p_kill <- renderPlot({
    long(D(), c("KCEF", "KVAN", "KRIF", "KTOT")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.7) +
      coord_cartesian(xlim = c(0, min(72, input$tend))) +
      labs(x = "시간 (h)", y = "살균속도 (1/h)", colour = NULL,
           title = "Emax 살균속도 — C/MIC 비에 대한 함수") + theme_abm
  })

  output$tbl_ster <- renderTable({
    e <- E()
    data.frame(
      항목 = c("CSF 유리균 멸균", "부착 아집단 제거", "간극"),
      시간_h = c(ifelse(is.na(e$sterile_h), "미달성", sprintf("%.1f", e$sterile_h)),
                 ifelse(is.na(e$adh_clear_h), "미달성", sprintf("%.1f", e$adh_clear_h)),
                 ifelse(is.na(e$adh_clear_h) | is.na(e$sterile_h), "—",
                        sprintf("%.1f", e$adh_clear_h - e$sterile_h))))
  })

  ## ---------------------------------------------------------------- tab 4
  output$p_product <- renderPlot({
    d <- D()
    d$factor_kkill <- d$KTOT / max(d$KTOT, 1e-9)
    d$factor_N <- d$NC / max(d$NC, 1e-9)
    d$factor_shield <- 1 - d$IDEX
    d$product <- d$LYSFLUX / max(d$LYSFLUX, 1e-9)
    long(d, c("factor_kkill", "factor_N", "factor_shield", "product")) %>%
      mutate(var = factor(var,
        levels = c("factor_kkill", "factor_N", "factor_shield", "product"),
        labels = c("k_kill (정규화)", "N (정규화)", "1 − E_dex (방패)",
                   "손상 flux (정규화)"))) %>%
      ggplot(aes(time, val, colour = var)) +
      geom_line(linewidth = 0.9) +
      coord_cartesian(xlim = c(0, min(48, input$tend))) +
      labs(x = "시간 (h)", y = "정규화 값", colour = NULL,
           title = "곱의 세 인자와 그 곱") + theme_abm
  })

  output$p_cargo <- renderPlot({
    long(D(), c("CW", "PLY")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.8) +
      coord_cartesian(xlim = c(0, min(96, input$tend))) +
      labs(x = "시간 (h)", y = "농도", colour = NULL,
           title = "방출된 화물 — 세포벽과 용소") + theme_abm
  })

  output$p_burst_zoom <- renderPlot({
    d <- D()
    ggplot(d, aes(time, LYSFLUX)) +
      geom_area(fill = "firebrick", alpha = 0.35) +
      geom_line(colour = "firebrick", linewidth = 0.8) +
      coord_cartesian(xlim = c(0, 12)) +
      labs(x = "시간 (h)", y = "용해 flux (CWU/mL/h)",
           title = "첫 12시간 — 봉우리는 여기 있다") + theme_abm
  })

  ## ---------------------------------------------------------------- tab 5
  output$p_cyto <- renderPlot({
    d <- D()
    vars <- c("TNF", "IL1", "IL6", "IL10", "CXCL8", "MMP9")
    long(d, vars) %>%
      ggplot(aes(time, val)) + geom_line(colour = "firebrick", linewidth = 0.8) +
      facet_wrap(~var, scales = "free_y", ncol = 3) +
      labs(x = "시간 (h)", y = "CSF 농도 (pg/mL, MMP9 은 ng/mL)") + theme_abm
  })

  output$p_pmn <- renderPlot({
    d <- D()
    p1 <- ggplot(d, aes(time, PMN)) + band("PMN") +
      geom_line(colour = "steelblue4", linewidth = 0.9) +
      labs(x = "시간 (h)", y = "CSF 호중구 (/µL)",
           title = "호중구 (음영 = 문헌 범위)") + theme_abm
    p2 <- long(d, c("MG", "ROS", "IDEX")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.8) +
      labs(x = "시간 (h)", y = NULL, colour = NULL,
           title = "대식세포 활성 · 산화부하 · 스테로이드 효과") + theme_abm
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  })

  ## ---------------------------------------------------------------- tab 6
  output$p_barrier <- renderPlot({
    d <- D()
    p1 <- ggplot(d, aes(time, PB)) +
      geom_line(colour = "darkcyan", linewidth = 1) +
      labs(x = NULL, y = "장벽 투과성 Pb", title = "문") + theme_abm
    p2 <- ggplot(d, aes(time, QALB)) + band("QALB") +
      geom_line(colour = "darkcyan", linewidth = 0.9) +
      labs(x = NULL, y = "Q_alb (×10⁻³)", title = "측정되는 문 — 알부민 지수") +
      theme_abm
    p3 <- ggplot(d, aes(time, VAN_CSF)) +
      geom_line(colour = "purple4", linewidth = 0.9) +
      geom_hline(yintercept = 4, linetype = 2) +
      annotate("text", x = Inf, y = 4.4, hjust = 1.05, size = 3,
               label = "4×MIC (살균 역치)") +
      labs(x = "시간 (h)", y = "반코마이신 CSF (mg/L)",
           title = "문이 닫히면 여기가 먼저 마른다") + theme_abm
    p4 <- ggplot(d, aes(time, VBR)) +
      geom_line(colour = "sienna", linewidth = 0.9) +
      labs(x = "시간 (h)", y = "뇌수분 초과 (mL)",
           title = "문이 닫히면 이것도 줄어든다") + theme_abm
    gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  })

  output$tbl_signflip <- renderTable({
    withProgress(message = "네 팔 비교 중...", {
      as.data.frame(claim3_signflip()) %>%
        transmute(팔 = arm,
                  멸균_h = round(sterile_h, 1),
                  VAN_CSF_AUC = round(AUC_van, 0),
                  최고_Pb = round(peak_Pb, 1),
                  최고_ICP = round(peak_ICP, 1),
                  청력_dB = round(hear_dB, 1),
                  사망_pct = round(100 * p_death, 1))
    })
  }, striped = TRUE)

  ## ---------------------------------------------------------------- tab 7
  output$p_panel <- renderPlot({
    d <- D()
    mk <- function(v, ylab, ttl) {
      g <- ggplot(d, aes(time, .data[[v]]))
      b <- band(v); if (!is.null(b)) g <- g + b
      g + geom_line(colour = "grey20", linewidth = 0.9) +
        labs(x = NULL, y = ylab, title = ttl) + theme_abm
    }
    gridExtra::grid.arrange(
      mk("PMN", "/µL", "CSF 백혈구"),
      mk("PROT_CSF", "mg/dL", "CSF 단백"),
      mk("GLC_CSF", "mg/dL", "CSF 포도당"),
      mk("LAC_CSF", "mmol/L", "CSF 락테이트"),
      mk("QALB", "×10⁻³", "알부민 지수"),
      ggplot(d, aes(time, GLCR)) +
        geom_hline(yintercept = 0.4, linetype = 2, colour = "firebrick") +
        geom_line(colour = "grey20", linewidth = 0.9) +
        labs(x = NULL, y = "비", title = "CSF/혈장 포도당비 (진단 역치 0.4)") +
        theme_abm,
      ncol = 3)
  })

  output$tbl_panel <- renderTable({
    e <- E()
    data.frame(
      지표 = c("최고 백혈구 (/µL)", "최고 단백 (mg/dL)", "최저 포도당 (mg/dL)",
               "최고 락테이트 (mmol/L)", "최고 Q_alb"),
      모델 = c(sprintf("%.0f", e$peak_PMN), sprintf("%.0f", e$peak_Prot),
               sprintf("%.1f", e$min_Glc), sprintf("%.1f", e$peak_Lac),
               sprintf("%.0f", e$peak_Qalb)),
      문헌범위 = c("1,000–5,000", "100–500", "<40 (비 <0.4)", ">3.5 (흔히 6–12)",
                   "30–100+"))
  }, striped = TRUE)

  ## ---------------------------------------------------------------- tab 8
  output$p_icp <- renderPlot({
    d <- D()
    p1 <- ggplot(d, aes(time, ICP)) + band("ICP") +
      geom_line(colour = "darkred", linewidth = 1) +
      geom_hline(yintercept = 25, linetype = 2) +
      labs(x = NULL, y = "ICP (mmHg)", title = "두개내압") + theme_abm
    p2 <- ggplot(d, aes(time, R_OUT)) +
      geom_line(colour = "chocolate4", linewidth = 0.9) +
      labs(x = NULL, y = "mmHg/(mL/h)", title = "CSF 유출저항") + theme_abm
    p3 <- ggplot(d, aes(time, VBR)) +
      geom_line(colour = "sienna", linewidth = 0.9) +
      labs(x = "시간 (h)", y = "mL", title = "뇌수분 초과") + theme_abm
    p4 <- long(d, c("MAP", "CPPo")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 60, linetype = 2, colour = "firebrick") +
      labs(x = "시간 (h)", y = "mmHg", colour = NULL,
           title = "MAP 과 CPP (점선 = CPP 목표 60)") + theme_abm
    gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  })

  output$p_perf <- renderPlot({
    d <- D()
    p1 <- ggplot(d, aes(time, CBFo)) +
      geom_hline(yintercept = 27.5, linetype = 2, colour = "firebrick") +
      geom_line(colour = "navy", linewidth = 0.9) +
      annotate("text", x = Inf, y = 29, hjust = 1.05, size = 3,
               label = "허혈 역치") +
      labs(x = "시간 (h)", y = "CBF (mL/100 g/min)", title = "뇌혈류") + theme_abm
    p2 <- long(d, c("AUTOR", "ISCHo")) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.9) +
      labs(x = "시간 (h)", y = NULL, colour = NULL,
           title = "자동조절 보존도와 허혈지수") + theme_abm
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  })

  ## ---------------------------------------------------------------- tab 9
  output$tbl_out <- renderTable({
    e <- E()
    d <- D()[nrow(D()), ]
    data.frame(
      엔드포인트 = c("사망확률 (%)", "청력역치 이동 (dB)", "어떤 난청 (>25 dB)",
                     "중증 난청 (>60 dB)", "인지 z", "국소결손 확률",
                     "피질 생존분율", "치상핵 생존분율", "유모세포 생존분율",
                     "미로 골화"),
      값 = c(sprintf("%.1f", 100 * e$p_death), sprintf("%.1f", e$hear_dB),
             ifelse(e$hear_dB > 25, "예", "아니오"),
             ifelse(e$hear_dB > 60, "예", "아니오"),
             sprintf("%.2f", e$cog_z), sprintf("%.2f", d$FOCAL),
             sprintf("%.3f", d$NCORT), sprintf("%.3f", d$NDG),
             sprintf("%.3f", d$HC), sprintf("%.3f", d$OSS)))
  }, striped = TRUE)

  output$p_hazard <- renderPlot({
    d <- D()[nrow(D()), ]
    p <- as.list(param(mod))
    df <- data.frame(
      성분 = c("ICP", "CPP", "SOFA", "허혈", "균혈증", "CSF 감염", "기저", "구조손상"),
      기여 = c(p$h_icp * d$I_ICP, p$h_cpp * d$I_CPP, p$h_sofa * d$I_SOFA,
               p$h_isch * d$I_ISCH, p$h_bact * d$I_BACT, p$h_ncsf * d$I_NCSF,
               p$h0 * d$time, p$h_cort_final * (1 - d$NCORT)))
    ggplot(df, aes(reorder(성분, 기여), 기여)) +
      geom_col(fill = "firebrick", alpha = 0.8) + coord_flip() +
      labs(x = NULL, y = "누적 위험 기여", title = "무엇이 실제로 죽이는가") +
      theme_abm
  })

  output$p_injury <- renderPlot({
    long(D(), c("NCORT", "NDG", "HC", "OSS")) %>%
      mutate(var = factor(var, levels = c("NCORT", "NDG", "HC", "OSS"),
                          labels = c("피질 뉴런", "해마 치상핵", "와우 유모세포",
                                     "미로 골화"))) %>%
      ggplot(aes(time, val, colour = var)) + geom_line(linewidth = 0.9) +
      labs(x = "시간 (h)", y = "생존분율 (골화는 진행분율)", colour = NULL) +
      theme_abm
  })

  ## ---------------------------------------------------------------- tab 10
  sweep <- eventReactive(input$run_sweep, {
    withProgress(message = "타이밍 스윕...", claim2_timing())
  })
  output$p_sweep <- renderPlot({
    s <- sweep()
    s2 <- s %>% pivot_longer(-dex_time_h)
    ggplot(s2, aes(dex_time_h, value)) +
      geom_line(colour = "purple4", linewidth = 0.9) +
      geom_point(colour = "purple4") +
      geom_vline(xintercept = 0, linetype = 2) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "덱사메타손 투여 시점 (항생제 기준, h)", y = NULL) + theme_abm
  })
  output$tbl_sweep <- renderTable(sweep(), striped = TRUE, digits = 3)

  ## ---------------------------------------------------------------- tab 11
  allsc <- eventReactive(input$run_all, {
    withProgress(message = "26 시나리오...", run_all_scenarios())
  })
  output$tbl_all <- renderTable({
    allsc() %>%
      transmute(id, 시나리오 = label,
                멸균_h = round(sterile_h, 1), 부착제거_h = round(adh_clear_h, 1),
                최고ICP = round(peak_ICP, 1), 최저CPP = round(min_CPP, 1),
                최저포도당 = round(min_Glc, 1), 최고Lac = round(peak_Lac, 1),
                청력dB = round(hear_dB, 1), 인지z = round(cog_z, 2),
                사망pct = round(100 * p_death, 1))
  }, striped = TRUE)

  coh <- eventReactive(input$run_coh, {
    withProgress(message = "코호트 20 회 적분...", run_cohort())
  })
  output$tbl_coh <- renderTable({
    c0 <- coh()
    rbind(
      c0 %>% transmute(환자 = as.character(id), 균량 = round(logN0, 1),
                       지연h = delay_h, 숙주방어 = host_def,
                       사망_DEX무 = round(100 * death_no_dex, 1),
                       사망_DEX유 = round(100 * death_dex, 1),
                       청력_DEX무 = round(dB_no_dex, 1),
                       청력_DEX유 = round(dB_dex, 1)),
      data.frame(환자 = "평균", 균량 = NA, 지연h = NA, 숙주방어 = NA,
                 사망_DEX무 = round(100 * mean(c0$death_no_dex), 1),
                 사망_DEX유 = round(100 * mean(c0$death_dex), 1),
                 청력_DEX무 = round(mean(c0$dB_no_dex), 1),
                 청력_DEX유 = round(mean(c0$dB_dex), 1)))
  }, striped = TRUE)

  ## ---------------------------------------------------------------- tab 12
  output$txt_checks <- renderPrint(structural_checks())

  output$txt_bugs <- renderUI({
    HTML(paste0(
      "<p>이 모델의 방정식은 R 런타임이 없는 환경에서 작성되었기 때문에, ",
      "먼저 의존성 없는 Python RK4 (<code>abm_reference_python.py</code>) 로 ",
      "독립 구현하고 실제로 적분해서 검증했다. 그 과정에서 <b>32개의 실제 결함</b>",
      "이 드러났다. 전체 목록은 그 파일의 머리말 [F1]–[F32] 에 있고, ",
      "실행 로그 전문은 <code>abm_reference_output.txt</code> 다. 구조를 바꾼 ",
      "것들만 옮기면:</p><ul>",
      "<li><b>F1</b> 포도당 소비항에 기질 의존성이 없어 CSF 포도당이 0 이 된 뒤에도 ",
      "소비가 계속되었다 → 락테이트 243 mmol/L (실측 6–12). 락테이트는 포도당에서만 ",
      "나온다.</li>",
      "<li><b>F14</b> 혈중 세균에 수용력이 없어 N_b 가 2×10¹⁰ CFU/mL 로 발산하고 ",
      "그것이 CSF 를 역오염시켰다.</li>",
      "<li><b>F15</b> 삼투압 환산에 분자량 182 대신 18.2 를 써서 만니톨 0.5 g/kg 이 ",
      "혈장 삼투압을 113 mOsm/kg 올렸다 (실제 11.3). 10배 단위오류.</li>",
      "<li><b>F16</b> 위험함수가 영구손상을 14일 내내 적분해서 회복한 환자도 ",
      "사망확률 99 % 가 되었고, 26개 시나리오 중 22개가 구분되지 않았다.</li>",
      "<li><b>F17</b> 자동조절 공식이 자동조절을 <i>해로운</i> 것으로 만들었다 — ",
      "CPP 38 에서 자동조절이 남은 환자가 압력수동 환자보다 더 허혈했다.</li>",
      "<li><b>F28</b> ICP·CPP 위험항이 선형이면 'CPP 45 를 100 h' 와 'CPP 5 를 10 h' ",
      "가 같은 위험이 된다. 경증의 지속적 이탈이 위험함수를 지배했다.</li>",
      "<li><b>F29</b> 세균 소멸 하한을 도함수 블록 안에 쓰면 불연속 스위치가 되어 ",
      "적분기를 멈춘다. mrgsolve 판은 상태를 건드리지 않고 보고 임계값으로만 ",
      "멸균을 판정한다.</li>",
      "<li><b>F30</b> 위험함수 계수를 눈으로 맞추는 대신 성분별 적분을 상태변수로 ",
      "추가해, 임상시험 목표치(폐렴균 사망 34 %/14 %)에 대해 <i>산술로</i> 풀었다. ",
      "해는 모든 계수에 대한 단일 배율 0.30 이었다.</li>",
      "</ul>"))
  })
}

shinyApp(ui, server)
