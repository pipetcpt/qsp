# =============================================================================
#  X-ALD QSP — Shiny 대시보드
#  X-linked Adrenoleukodystrophy · interactive QSP explorer
#
#  실행 (Run):
#     setwd("x-linked-adrenoleukodystrophy")
#     shiny::runApp("xald_shiny_app.R")
#
#  필요 패키지: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
#
#  이 앱의 설계 의도
#  -----------------------------------------------------------------------------
#  탭 2 ("두 변수")가 이 앱의 존재 이유입니다. 사용자가 어떤 조합을 돌려도
#  화면 왼쪽(측정되는 혈장 C26:0)과 오른쪽(측정되지 않는 뇌염증 스위치)이
#  **서로 따로 움직이는 것**을 보게 됩니다. 로렌조 오일은 왼쪽만 움직이고,
#  유전자치료는 오른쪽만 움직입니다. 그 비대응이 이 질환의 임상 현실입니다.
# =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

source("xald_mrgsolve_model.R")

MOD <- xald_mod()
assign(".XALD_MOD", MOD, envir = globalenv())

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom")

# -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("X-연관 부신백질형성장애 (X-ALD) — QSP 탐색 도구"),
  tags$p(tags$em(paste(
    "하나의 병변(ABCD1 소실), 두 개의 독립 변수 —",
    "측정 가능한 VLCFA 축적(주변부에서 도달)과",
    "측정 불가능한 이중안정 뇌염증 스위치(CNS 내부에서만 도달)."))),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① 환자 (유전형)"),
      sliderInput("mutres", "잔존 ALDP 기능 (MUTRES)", 0, 1, 0.02, step = 0.01),
      sliderInput("fmos", "ALDP 결핍 세포 분율 (여성 보유자 FMOS)", 0, 1, 0, step = 0.05),
      sliderInput("susc", "뇌형 감수성 (SUSCTV)", 0.05, 1.6, 0.46, step = 0.01),
      helpText(HTML(paste0(
        "임계값 <b>SUSC* &asymp; 0.41</b> (모델의 출력, 이분법으로 탐색).<br>",
        "그 아래 = 순수 AMN, 그 위 = 뇌형. 임계값에 가까울수록 발병이 <b>늦습니다</b>",
        " (critical slowing down)."))),
      hr(),
      h4("② 치료"),
      checkboxInput("lo",   "로렌조 오일 (경구, 말초 ELOVL1 경쟁적 저해)", FALSE),
      sliderInput("lo_start", "로렌조 오일 시작 (세)", 0.5, 20, 1.5, step = 0.5),
      checkboxInput("diet", "VLCFA 제한식", FALSE),
      checkboxInput("gt",   "유전자치료 / HSCT", FALSE),
      sliderInput("gt_loes", "이식 시점 (Loes 점수)", 1, 20, 2, step = 1),
      sliderInput("vcn", "벡터 카피수 VCN (0 = 동종 이식)", 0, 2, 0.8, step = 0.1),
      sliderInput("bu", "부설판 노출 배수", 0.5, 1.5, 1.0, step = 0.02),
      checkboxInput("hc",   "하이드로코르티손 보충", FALSE),
      checkboxInput("leri", "레리글리타존 (CNS 침투 PPARγ)", FALSE),
      checkboxInput("sob",  "소베티롬 계열 (CNS 침투 대사 표적)", FALSE),
      checkboxInput("ster", "고용량 스테로이드", FALSE),
      checkboxInput("anti", "항산화 3제 + DMF", FALSE),
      hr(),
      sliderInput("years", "관찰 기간 (년)", 5, 60, 20, step = 5),
      actionButton("run", "시뮬레이션 실행", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. 환자 요약",
                 br(), verbatimTextOutput("summary"),
                 h5("표현형 판정"), verbatimTextOutput("phenotype")),
        tabPanel("2. 두 변수 (핵심)",
                 br(), plotOutput("twovar", height = "460px"),
                 helpText(paste("왼쪽은 측정되는 것, 오른쪽은 측정되지 않는 것입니다.",
                                "치료를 켜고 두 패널이 함께 움직이는지 확인하십시오."))),
        tabPanel("3. VLCFA 대사 · 약물 노출",
                 br(), plotOutput("metab", height = "440px")),
        tabPanel("4. 뇌 스위치 · 이득",
                 br(), plotOutput("switchplot", height = "440px"),
                 helpText("LOOPGAIN > 1 이면 염증이 점화 입력 없이도 스스로 유지됩니다.")),
        tabPanel("5. MRI · 신경 엔드포인트",
                 br(), plotOutput("neuro", height = "440px")),
        tabPanel("6. 부신 축",
                 br(), plotOutput("adrenal", height = "440px")),
        tabPanel("7. 척수형 (AMN)",
                 br(), plotOutput("cord", height = "440px")),
        tabPanel("8. 이식 · 유전자치료",
                 br(), plotOutput("transplant", height = "440px"),
                 helpText(paste("미세아교세포 치환은 느립니다 (KMGREP).",
                                "그 느림이 이식 후 12-18개월 지체 창의 원인입니다."))),
        tabPanel("9. 시나리오 비교",
                 br(),
                 selectInput("scens", "비교할 시나리오", multiple = TRUE,
                             choices = names(xald_scenarios()),
                             selected = c("s10_elicel_early",
                                          "s11_elicel_early_control",
                                          "s07_lo_ccald")),
                 selectInput("scvar", "표시할 변수",
                             choices = c("LOES","NFS","C26","C26LYSO","EDSS",
                                         "MGLCORR","CORTISOL","CHIMER"),
                             selected = "LOES"),
                 actionButton("runsc", "시나리오 실행"),
                 plotOutput("scenplot", height = "420px")),
        tabPanel("10. 집단 · 표현형 분포",
                 br(),
                 numericInput("npop", "개체 수", 150, min = 20, max = 600, step = 10),
                 actionButton("runpop", "집단 시뮬레이션"),
                 plotOutput("popplot", height = "380px"),
                 verbatimTextOutput("poptext")),
        tabPanel("11. 검증 앵커",
                 br(), actionButton("runval", "검증 스위트 실행 (A1-A23)"),
                 DTOutput("valtab"))
      )
    )
  )
)

# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  build <- reactive({
    p <- list(MUTRES = input$mutres, FMOS = input$fmos, SUSCTV = input$susc,
              VCNIN = input$vcn,
              DIETSC = if (input$diet) 0.35 else 1.0,
              STEROIDX = if (input$ster) 1.0 else 0.0,
              NACD = if (input$anti) 1.0 else 0.0,
              DMFD = if (input$anti) 1.0 else 0.0)
    end <- input$years * YR
    ev_list <- list(); txday <- NULL
    if (input$lo)   ev_list <- c(ev_list, list(.lo_regimen(input$lo_start*YR, end)))
    if (input$hc)   ev_list <- c(ev_list, list(.hc_regimen(7.5*YR, end)))
    if (input$leri) ev_list <- c(ev_list, list(.leri_regimen(2*YR, end)))
    if (input$sob)  ev_list <- c(ev_list, list(.sob_regimen(2*YR, end)))
    if (input$gt) {
      # 이식 시점은 사용자가 고른 Loes 값에 도달하는 날 — 즉 짝지은 자연사의 출력
      p0  <- list(MUTRES = input$mutres, FMOS = input$fmos, SUSCTV = input$susc)
      nat <- MOD %>% param(p0) %>% init(xald_birth(MOD, p0)) %>% zero_re() %>%
        mrgsim(end = end, delta = 10, maxsteps = 1e6) %>% as.data.frame()
      i <- which(nat$LOES >= input$gt_loes)
      if (length(i)) {
        ev_list <- c(ev_list, list(.hsct_regimen(nat$time[i[1]],
                                                 vcn = input$vcn,
                                                 bu_scale = input$bu)))
        txday <- nat$time[i[1]]
      }
    }
    e <- if (length(ev_list)) do.call(c, ev_list) else NULL
    list(param = p, ev = e, end = end, txday = txday)
  })

  simdata <- eventReactive(input$run, {
    b <- build()
    init0 <- xald_birth(MOD, b$param)
    m <- MOD %>% param(b$param) %>% init(init0) %>% zero_re()
    d <- if (is.null(b$ev)) {
      m %>% mrgsim(end = b$end, delta = 10, maxsteps = 1e6)
    } else {
      m %>% mrgsim(events = b$ev, end = b$end, delta = 10, maxsteps = 1e6)
    }
    d <- as.data.frame(d)
    d$age <- d$time / YR
    attr(d, "txday") <- b$txday
    d
  }, ignoreNULL = FALSE)

  # ---- 1. summary ----
  output$summary <- renderPrint({
    d <- simdata()
    last <- d[nrow(d), ]
    cat(sprintf("관찰 기간            : %.0f년\n", max(d$age)))
    cat(sprintf("혈장 C26:0           : %.3f umol/L  (정상 0.35, 상승 %.2f배)\n",
                last$C26, last$C26FOLD))
    cat(sprintf("C26:0-lysoPC         : %.3f umol/L\n", last$C26LYSO))
    cat(sprintf("뇌 VLCFA 투과 게이트  : %.3f  (환자에서 포화 -> 표현형 판별 불가)\n",
                last$VLCGATE))
    cat(sprintf("최고 Loes / NFS      : %.1f / %.1f\n", max(d$LOES), max(d$NFS)))
    cat(sprintf("최고 EDSS            : %.2f\n", max(d$EDSS)))
    cat(sprintf("코르티솔 (최종)       : %.0f nmol/L,  ACTH %.0f pg/mL\n",
                last$CORTISOL, last$ACTHOUT))
    cat(sprintf("부신 예비능           : %.3f  (부신부족 역치 0.36)\n", last$ADRRESV))
    cat(sprintf("단핵구 키메리즘        : %.3f,  교정 미세아교세포 %.3f\n",
                last$CHIMER, last$MGLCORR))
    if (max(d$PMDS) > 0)
      cat(sprintf("MDS/AML 누적 확률      : %.3f\n", max(d$PMDS)))
    if (!is.null(attr(d, "txday")))
      cat(sprintf("이식 시점             : %.2f세\n", attr(d, "txday")/YR))
  })

  output$phenotype <- renderPrint({
    d <- simdata()
    fired <- max(d$MGP) > 0.15
    cat(if (input$mutres >= 0.9) "비보유자 / 정상\n"
        else if (fired) "뇌형 (cerebral ALD) — 스위치 점화됨\n"
        else if (input$fmos > 0) "여성 보유자 표현형 (척수형 중심)\n"
        else "순수 AMN (척수형) — 스위치 미점화\n")
    if (fired) {
      i <- which(d$LOES >= 1)
      if (length(i)) cat(sprintf("뇌형 발병 (Loes >= 1) : %.2f세\n", d$age[i[1]]))
      j <- which(d$LOES >= 15)
      if (length(i) && length(j))
        cat(sprintf("Loes 1 -> 15 소요     : %.2f년\n", (d$age[j[1]] - d$age[i[1]])))
    }
    k <- which(d$ADRINSUF >= 1)
    if (length(k)) cat(sprintf("부신부족 발생         : %.2f세\n", d$age[k[1]]))
    cat(sprintf("주요 기능장애 (MFD)   : %s\n", if (max(d$MFD) > 0) "발생" else "없음"))
  })

  # ---- 2. the two variables ----
  output$twovar <- renderPlot({
    d <- simdata()
    a <- ggplot(d, aes(age)) +
      geom_line(aes(y = C26, colour = "혈장 C26:0"), linewidth = 1.1) +
      geom_line(aes(y = C26LYSO*2, colour = "C26:0-lysoPC (x2)"), linewidth = 1) +
      geom_hline(yintercept = 0.35, linetype = 2, colour = "grey40") +
      annotate("text", x = max(d$age)*0.6, y = 0.40, label = "정상 C26:0",
               size = 3.2, colour = "grey30") +
      scale_colour_manual(values = c("혈장 C26:0" = "#c0392b",
                                     "C26:0-lysoPC (x2)" = "#e67e22")) +
      labs(title = "변수 (1) 측정되는 것 — 주변부에서 도달 가능",
           x = "연령 (세)", y = "umol/L", colour = NULL) + THEME
    b <- ggplot(d, aes(age)) +
      geom_line(aes(y = MGP, colour = "활성 미세아교세포 (스위치)"), linewidth = 1.2) +
      geom_line(aes(y = DEMYELIN, colour = "탈수초 분율"), linewidth = 1) +
      geom_line(aes(y = GADENH, colour = "가돌리늄 조영증강"), linewidth = 0.9,
                linetype = 3) +
      scale_colour_manual(values = c("활성 미세아교세포 (스위치)" = "#8e44ad",
                                     "탈수초 분율" = "#2471a3",
                                     "가돌리늄 조영증강" = "#16a085")) +
      labs(title = "변수 (2) 측정되지 않는 것 — CNS 내부에서만 도달 가능",
           x = "연령 (세)", y = "분율", colour = NULL) + THEME
    gridExtra_arrange(a, b)
  })

  # ---- 3. metabolism ----
  output$metab <- renderPlot({
    d <- simdata()
    long <- d %>%
      select(age, `혈장 C26:0` = C26, `뇌 백질 VLCFA` = CBR,
             `척수 VLCFA` = CSC, `부신피질 VLCFA` = CADR,
             `에루크산 (umol/L /100)` = ERUCIC) %>%
      mutate(`에루크산 (umol/L /100)` = `에루크산 (umol/L /100)`/100) %>%
      pivot_longer(-age)
    ggplot(long, aes(age, value, colour = name)) +
      geom_line(linewidth = 1.05) +
      labs(title = "조직별 VLCFA와 약물 노출 (정상 = 1.0 로 정규화된 조직 풀)",
           subtitle = paste("뇌·척수 VLCFA의 약 88%는 국소 합성 —",
                            "에루크산은 BBB를 넘지 못하므로 여기에 닿지 못합니다"),
           x = "연령 (세)", y = "상대값", colour = NULL) + THEME
  })

  # ---- 4. switch ----
  output$switchplot <- renderPlot({
    d <- simdata()
    a <- ggplot(d, aes(age)) +
      geom_line(aes(y = MGP), colour = "#8e44ad", linewidth = 1.2) +
      geom_hline(yintercept = 0.15, linetype = 2, colour = "grey40") +
      labs(title = "활성 미세아교세포 분율 (이중안정 상태변수)",
           x = "연령 (세)", y = "MGP") + THEME
    b <- ggplot(d, aes(age)) +
      geom_line(aes(y = LOOPGAIN), colour = "#c0392b", linewidth = 1.2) +
      geom_hline(yintercept = 1, linetype = 2) +
      labs(title = "자기증폭 루프 이득 (1 = 되돌릴 수 없는 경계)",
           x = "연령 (세)", y = "LOOPGAIN") + THEME
    gridExtra_arrange(a, b)
  })

  # ---- 5. neuro endpoints ----
  output$neuro <- renderPlot({
    d <- simdata()
    long <- d %>% select(age, `Loes (0-34)` = LOES, `NFS (0-25)` = NFS) %>%
      pivot_longer(-age)
    ggplot(long, aes(age, value, colour = name)) +
      geom_line(linewidth = 1.15) +
      geom_hline(yintercept = 9, linetype = 3, colour = "#e67e22") +
      annotate("text", x = min(d$age) + 1, y = 10.2, hjust = 0,
               label = "Loes 9 = 이식 치료 창의 상한", size = 3.2,
               colour = "#e67e22") +
      labs(title = "영상·신경 기능 점수 (한 방향으로만 누적)",
           x = "연령 (세)", y = "점수", colour = NULL) + THEME
  })

  # ---- 6. adrenal ----
  output$adrenal <- renderPlot({
    d <- simdata()
    a <- ggplot(d, aes(age)) +
      geom_line(aes(y = CORTISOL, colour = "코르티솔 (nmol/L)"), linewidth = 1.1) +
      geom_line(aes(y = ACTHOUT, colour = "ACTH (pg/mL)"), linewidth = 1.1) +
      scale_colour_manual(values = c("코르티솔 (nmol/L)" = "#2471a3",
                                     "ACTH (pg/mL)" = "#c0392b")) +
      labs(title = "부신피질 축", x = "연령 (세)", y = NULL, colour = NULL) + THEME
    b <- ggplot(d, aes(age)) +
      geom_line(aes(y = ADRRESV, colour = "스테로이드 생성 예비능"), linewidth = 1.1) +
      geom_line(aes(y = POTASSIUM/10, colour = "혈청 K+ / 10"), linewidth = 1.1) +
      geom_hline(yintercept = 0.36, linetype = 2, colour = "grey40") +
      scale_colour_manual(values = c("스테로이드 생성 예비능" = "#8e44ad",
                                     "혈청 K+ / 10" = "#e67e22")) +
      labs(title = "예비능과 전해질 (점선 = 부신부족 역치)",
           x = "연령 (세)", y = NULL, colour = NULL) + THEME
    gridExtra_arrange(a, b)
  })

  # ---- 7. cord ----
  output$cord <- renderPlot({
    d <- simdata()
    long <- d %>% select(age, `EDSS` = EDSS, `척수 손상 분율` = CORDLOSS,
                         `6분 보행 (m/100)` = WALK6MIN) %>%
      mutate(`6분 보행 (m/100)` = `6분 보행 (m/100)`/100) %>%
      pivot_longer(-age)
    ggplot(long, aes(age, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      labs(title = "척수 축삭병증 (AMN) — 스위치와 무관하게 진행",
           subtitle = "가장 긴 축삭이 먼저 (LENCST > LENDC), CNS 축삭 재생은 0",
           x = "연령 (세)", y = NULL, colour = NULL) + THEME
  })

  # ---- 8. transplant ----
  output$transplant <- renderPlot({
    d <- simdata()
    long <- d %>% select(age, `단핵구 키메리즘` = CHIMER,
                         `교정 미세아교세포` = MGLCORR,
                         `전체 미세아교세포` = MGLTOT,
                         `호중구` = NEUTRO,
                         `MDS/AML 누적확률` = PMDS) %>%
      pivot_longer(-age)
    p <- ggplot(long, aes(age, value, colour = name)) +
      geom_line(linewidth = 1.05) +
      labs(title = "이식 후 재구성 — 빠른 단핵구, 느린 미세아교세포",
           subtitle = "두 곡선 사이의 간격이 곧 치료 지체 창입니다",
           x = "연령 (세)", y = "분율 / 확률", colour = NULL) + THEME
    if (!is.null(attr(d, "txday")))
      p <- p + geom_vline(xintercept = attr(d, "txday")/YR, linetype = 2)
    p
  })

  # ---- 9. scenarios ----
  scendata <- eventReactive(input$runsc, {
    req(input$scens)
    bind_rows(lapply(input$scens, function(nm) {
      d <- xald_run(nm, MOD, delta = 20)
      data.frame(age = d$time/YR, value = d[[input$scvar]], scenario = nm)
    }))
  })
  output$scenplot <- renderPlot({
    ggplot(scendata(), aes(age, value, colour = scenario)) +
      geom_line(linewidth = 1.1) +
      labs(title = paste0("시나리오 비교 — ", input$scvar),
           x = "연령 (세)", y = input$scvar, colour = NULL) + THEME
  })

  # ---- 10. population ----
  popdata <- eventReactive(input$runpop, {
    xald_population(MOD, n = input$npop)
  })
  output$popplot <- renderPlot({
    p <- popdata()
    ggplot(p, aes(SUSC, maxLOES, colour = factor(cerebral))) +
      geom_point(alpha = 0.75, size = 2) +
      geom_vline(xintercept = 0.41, linetype = 2) +
      scale_colour_manual(values = c("0" = "#2471a3", "1" = "#c0392b"),
                          labels = c("척수형 (AMN)", "뇌형"), name = NULL) +
      labs(title = "감수성 하나의 개체간 변이가 두 개의 병을 만든다",
           subtitle = "점선 = 이분법으로 찾은 임계값. 혈장 C26:0는 모든 점에서 동일합니다.",
           x = "감수성 SUSC", y = "최고 Loes 점수") + THEME
  })
  output$poptext <- renderPrint({
    p <- popdata()
    cat(sprintf("뇌형 전환 비율 : %.1f%%   (문헌 35-40%%)\n", 100*mean(p$cerebral)))
    cat(sprintf("부신부족 비율  : %.1f%%   (문헌 약 80%% — 모델은 과대예측)\n",
                100*mean(p$adrinsuf)))
    o <- p$onset[!is.na(p$onset)]
    if (length(o))
      cat(sprintf("뇌형 발병 연령 : 중간 %.1f세 (범위 %.1f-%.1f)  (문헌 4-8세)\n",
                  median(o), min(o), max(o)))
  })

  # ---- 11. validation ----
  valdata <- eventReactive(input$runval, { xald_validate(MOD, verbose = FALSE) })
  output$valtab <- renderDT({
    datatable(valdata(), rownames = FALSE,
              options = list(pageLength = 25, dom = "t")) %>%
      formatStyle("pct", color = styleInterval(c(-20, 20),
                                              c("#c0392b", "#1a7f37", "#c0392b")))
  })
}

# Two-panel helper without taking a hard dependency on gridExtra
gridExtra_arrange <- function(p1, p2) {
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  } else {
    p1   # gridExtra 가 없으면 첫 패널만 표시
  }
}

shinyApp(ui, server)
