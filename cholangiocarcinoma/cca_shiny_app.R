# =============================================================================
#  cca_shiny_app.R
#  Cholangiocarcinoma QSP model — interactive dashboard
# =============================================================================
#
#  13 tabs.  The app is organised around the three claims the model makes,
#  not around the organ systems:
#
#    Tab 4  THE DOSE GATE      — the only tab that matters if claim I is right.
#                                It shows the prescribed dose, the four gates,
#                                and what was actually delivered.  Move the
#                                drainage slider and watch the delivered dose
#                                change without touching a single antitumour
#                                parameter.
#    Tab 3  CLONES             — T_S / T_P / T_R plotted separately, because a
#                                total-volume plot hides the entire mechanism.
#    Tab 11 TWO HAZARDS        — h_tumour and h_biliary on the same axes, and
#                                the fraction of the total hazard that is
#                                biliary.  If that fraction is large the
#                                patient's problem is a stent, not a drug.
#    Tab 13 FALSIFIERS         — one switch each.  Turning the gate off must
#                                make drainage timing stop mattering.
#
#  Run:  shiny::runApp("cca_shiny_app.R")
#  Requires: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
# =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

# The model definition lives in cca_mrgsolve_model.R; sourcing it compiles the
# model and exposes `mod`, `stent_events()` and the scenario drivers.
source("cca_mrgsolve_model.R")

PAL <- c(Sensitive = "#4527A0", Persister = "#00838F", Resistant = "#C62828",
         Total = "#212121", Metastatic = "#EF6C00")

thm <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

# -----------------------------------------------------------------------------
#  UI
# -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("담관암 QSP 모델 · Cholangiocarcinoma QSP model"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("<b>전달된 용량은 입력이 아니라 출력이다</b> · ",
              "<b>내성은 유도되는 것이 아니라 선택된다</b> · ",
              "<b>생존은 서로 다른 시계를 가진 두 개의 위험이다</b>")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 · Patient"),
      sliderInput("tt0", "진단 시 종양 부피 Tumour volume (cm³)",
                  20, 400, 100, step = 10),
      sliderInput("fhilar", "간문부 위치 분율 f_hilar (0 = 간내, 1 = 간문부)",
                  0.02, 0.97, 0.72, step = 0.02),
      checkboxInput("fgfr2", "FGFR2 융합 (FGFR2 fusion)", FALSE),
      checkboxInput("idh1", "IDH1 R132 변이", FALSE),
      checkboxInput("immeng", "면역 관여 종양 (immune-engaged)", FALSE),
      sliderInput("lam0", "종양 증식 속도 LAM0 (1/d)",
                  0.004, 0.025, 0.010, step = 0.001),

      hr(), h4("담도 배액 · Biliary drainage"),
      selectInput("stent", "스텐트", c("금속 SEMS" = "sems",
                                       "플라스틱 plastic" = "plastic",
                                       "없음 none" = "none")),
      sliderInput("stent_delay", "배액 지연 (일) — 항종양 파라미터는 그대로",
                  0, 120, 0, step = 10),

      hr(), h4("전신 치료 · Systemic therapy"),
      checkboxInput("gemcis", "젬시타빈/시스플라틴 (1000 / 25 mg/m², d1·d8 q21)",
                    TRUE),
      sliderInput("ncycle", "주기 수 cycles", 0, 12, 8),
      checkboxInput("durva", "더발루맙 1500 mg", FALSE),
      selectInput("fgi", "FGFR 억제제",
                  c("없음" = "none", "페미가티닙 14/7" = "pem",
                    "페미가티닙 연속" = "pemcont", "푸티바티닙 연속" = "fut")),
      checkboxInput("ivo", "아이보시데닙 500 mg", FALSE),
      checkboxInput("cape", "보조 카페시타빈 6개월 (BILCAP)", FALSE),
      sliderInput("k2l", "2차 항암 강도 K_2L (1/d)", 0, 0.03, 0.010, step = 0.002),

      hr(), h4("반증 스위치 · Falsifiers"),
      checkboxInput("f1", "F1 · 용량 게이트 끄기 (GATE_ON = 0)", FALSE),
      checkboxInput("f2", "F2 · 사전 내성 클론 제거 (MU_RES = 0)", FALSE),
      checkboxInput("f4", "F4 · 기질 투과 장벽 제거 (FPEN_MIN = 1)", FALSE),

      hr(),
      sliderInput("tend", "시뮬레이션 기간 (일)", 180, 2400, 1100, step = 60)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",

        # ---- 1 -------------------------------------------------------------
        tabPanel("1 · 환자 요약",
          h4("입력된 환자와 처방"),
          verbatimTextOutput("summary_txt"),
          h4("모델이 주장하는 세 가지"),
          tags$ol(
            tags$li(HTML("<b>전달 용량은 출력이다.</b> 빌리루빈·호중구·크레아티닌 청소율·수행상태가 ",
                         "투여 시점에 실제 투여량을 결정한다. 종양 → 폐색 → 빌리루빈 → 용량 보류 → ",
                         "종양은 닫힌 양성 되먹임 고리이고, 배액은 그 고리를 끊는 유일한 간선이다.")),
            tags$li(HTML("<b>내성은 선택된다.</b> FGFR2 kinase-domain 변이 클론은 t = 0에 ",
                         "MU_RES·ln(N) 빈도로 이미 존재한다. '내성 발생 시점' 파라미터는 없다.")),
            tags$li(HTML("<b>생존은 두 개의 위험이다.</b> h_tumour는 느리고 부피에 비례하며, ",
                         "h_biliary는 빠르고 재발성이며 배액에 의해 결정된다."))
          ),
          plotOutput("p_overview", height = 380)
        ),

        # ---- 2 -------------------------------------------------------------
        tabPanel("2 · PK",
          h4("혈장 및 세포내 약물동태"),
          plotOutput("p_pk", height = 620),
          tags$p(style = "color:#666",
                 "젬시타빈의 혈장 반감기는 약 6분이므로 효과를 구동하는 것은 혈장 농도가 아니라 ",
                 "세포내 dFdCTP이다. 백금 역시 유리 백금이 아니라 Pt-DNA 부가체가 효과 구동 변수다.")
        ),

        # ---- 3 -------------------------------------------------------------
        tabPanel("3 · 클론 (핵심)",
          h4("민감 클론 · 약물내성 지속세포 · FGFR2 변이 클론"),
          plotOutput("p_clone", height = 420),
          plotOutput("p_clonefrac", height = 260),
          tags$p(style = "color:#666",
                 "총 부피만 그리면 기전이 전부 사라진다. 치료 중 부피가 평평해 보이는 구간은 ",
                 "대부분 T_S가 죽고 T_P로 전환되는 구간이며, 치료 종료 후의 재증식 속도는 ",
                 "T_P의 증식 속도(GP × LAM0)로 결정된다.")
        ),

        # ---- 4 -------------------------------------------------------------
        tabPanel("4 · 용량 게이트 ★",
          h4("규칙 1-4 · 처방된 용량과 실제 전달된 용량"),
          fluidRow(
            column(6, plotOutput("p_gate", height = 380)),
            column(6, plotOutput("p_rdi", height = 380))
          ),
          h4("전달 상대 용량 강도 (RDI)"),
          verbatimTextOutput("rdi_txt"),
          tags$p(style = "color:#666",
                 "배액 지연 슬라이더만 움직여 보라. 항종양 파라미터는 하나도 바뀌지 않는데 ",
                 "RDI가 떨어지고 생존이 짧아진다. F1 스위치를 켜면 그 효과가 사라져야 한다.")
        ),

        # ---- 5 -------------------------------------------------------------
        tabPanel("5 · 담도 폐색 · 배액",
          plotOutput("p_biliary", height = 620),
          tags$p(style = "color:#666",
                 "스텐트 개통성(PATN)은 생물막·슬러지·종양 내증식으로 감쇠한다. ",
                 "금속 스텐트 t½ ≈ 240일, 플라스틱 t½ ≈ 90일.")
        ),

        # ---- 6 -------------------------------------------------------------
        tabPanel("6 · 간 예비능 · ALBI",
          plotOutput("p_liver", height = 560),
          DTOutput("t_albi")
        ),

        # ---- 7 -------------------------------------------------------------
        tabPanel("7 · 혈액학적 독성",
          plotOutput("p_heme", height = 560),
          tags$p(style = "color:#666",
                 "Friberg 구조: 증식 풀 → 3개의 전이 구획 → 순환 호중구, MTT ≈ 125시간, ",
                 "되먹임 지수 γ = 0.16. 호중구 최저치는 규칙 2를 통해 8일째 용량을 잃게 만든다.")
        ),

        # ---- 8 -------------------------------------------------------------
        tabPanel("8 · 면역",
          plotOutput("p_immune", height = 560),
          tags$p(style = "color:#666",
                 "PI_IMMUNE은 용량이 아니라 혼합 지표다. 면역 관여 종양에서만 CD8 효과세포가 ",
                 "증식하며, 이것이 중앙값을 거의 움직이지 않으면서 24개월 꼬리를 벌리는 이유다.")
        ),

        # ---- 9 -------------------------------------------------------------
        tabPanel("9 · 바이오마커",
          plotOutput("p_bio", height = 620),
          tags$p(style = "color:#666",
                 "CA 19-9는 담즙정체만으로도 상승한다 — 황달 상태에서는 종양 지표로 신뢰할 수 없다. ",
                 "인산은 FGFR 표적 결합의 무료 약력학 지표이고, 2-HG는 IDH1 표적 결합 지표다.")
        ),

        # ---- 10 ------------------------------------------------------------
        tabPanel("10 · RECIST · 임상 종말점",
          plotOutput("p_recist", height = 420),
          verbatimTextOutput("endpoint_txt")
        ),

        # ---- 11 ------------------------------------------------------------
        tabPanel("11 · 두 개의 위험 ★",
          plotOutput("p_hazard", height = 420),
          plotOutput("p_surv", height = 300),
          tags$p(style = "color:#666",
                 "위험의 담도 분율이 높다면 이 환자에게 필요한 것은 약이 아니라 스텐트다. ",
                 "이 그림이 모델의 임상적 결론이다.")
        ),

        # ---- 12 ------------------------------------------------------------
        tabPanel("12 · 시나리오 비교",
          h4("표준 시나리오 일괄 실행"),
          actionButton("run_scen", "시나리오 실행 (수 초 소요)"),
          br(), br(),
          plotOutput("p_scen", height = 420),
          DTOutput("t_scen")
        ),

        # ---- 13 ------------------------------------------------------------
        tabPanel("13 · 반증 시험 ★",
          h4("파라미터 하나씩, 나머지는 재적합하지 않음"),
          tableOutput("t_fals"),
          tags$ul(
            tags$li(HTML("<b>F1 GATE_ON = 0</b> — 빌리루빈이 용량을 제한하지 않는다. ",
                         "예측: 배액 지연의 생존 손실이 거의 사라져야 한다.")),
            tags$li(HTML("<b>F2 MU_RES = 0</b> — 사전 내성 클론이 없다. ",
                         "결과(보고됨): 관찰 기간 안에서는 거의 아무것도 바뀌지 않는다. ",
                         "이것은 모델의 <i>음성 결과</i>이며 숨기지 않았다.")),
            tags$li(HTML("<b>F3 PI_IMMUNE = 1</b> — 모든 종양이 면역 관여. ",
                         "예측: 더발루맙이 <i>중앙값</i>을 움직여야 하는데, TOPAZ-1은 아니라고 말한다.")),
            tags$li(HTML("<b>F4 FPEN_MIN = 1</b> — 기질 투과 장벽 없음. ",
                         "예측: 젬시타빈/시스플라틴 반응률이 과녁을 크게 넘어야 한다."))
          )
        )
      )
    )
  )
)

# -----------------------------------------------------------------------------
#  SERVER
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  par_list <- reactive({
    p <- list(
      TT0 = input$tt0, FHILAR = input$fhilar, LAM0 = input$lam0,
      FGFR2 = as.numeric(input$fgfr2), IDH1 = as.numeric(input$idh1),
      IMMENG = as.numeric(input$immeng),
      GEM_MGM2 = if (input$gemcis) 1000 else 0,
      CIS_MGM2 = if (input$gemcis) 25 else 0,
      NCYCLE   = if (input$gemcis) input$ncycle else 0,
      DUR_MG   = if (input$durva) 1500 else 0,
      IVO_MG   = if (input$ivo) 500 else 0,
      CAP_MGM2 = if (input$cape) 1250 else 0,
      CAP_DAYS = if (input$cape) 180 else 0,
      K_2L     = input$k2l,
      KOCC     = if (input$stent == "plastic") 0.0077 else 0.0029,
      GATE_ON  = if (input$f1) 0 else 1,
      MU_RES   = if (input$f2) 0 else 1e-6,
      FPEN_MIN = if (input$f4) 1.0 else 0.30
    )
    p <- c(p, switch(input$fgi,
      none    = list(FGI_MG = 0),
      pem     = list(FGI_MG = 13.5, FGI_ON = 14, FGI_OFF = 7, COVAL = 0, RHO = 0.05),
      pemcont = list(FGI_MG = 13.5, FGI_ON = 1,  FGI_OFF = 0, COVAL = 0, RHO = 0.05),
      fut     = list(FGI_MG = 20,   FGI_ON = 1,  FGI_OFF = 0, COVAL = 1, RHO = 0.45)))
    p
  })

  sim <- reactive({
    m <- param(mod, par_list())
    if (input$stent_delay > 0) m <- init(m, PATN = 0)
    ev <- stent_events(input$stent, input$stent_delay, input$tend)
    out <- if (is.null(ev)) mrgsim(m, end = input$tend, delta = 1)
           else             mrgsim(m, events = ev, end = input$tend, delta = 1)
    as.data.frame(out)
  })

  long <- function(d, cols, key = "var")
    d %>% select(time, all_of(cols)) %>% pivot_longer(-time, names_to = key)

  # ---- 1 --------------------------------------------------------------------
  output$summary_txt <- renderText({
    p <- par_list()
    paste0(
      "종양 ", input$tt0, " cm³ · f_hilar ", input$fhilar,
      " · LAM0 ", input$lam0, " /d\n",
      "유전형: ", if (input$fgfr2) "FGFR2 fusion " else "",
                  if (input$idh1) "IDH1 R132 " else "",
                  if (!input$fgfr2 && !input$idh1) "actionable driver 없음" else "", "\n",
      "면역 관여: ", if (input$immeng) "예 (mixture member)" else "아니오", "\n",
      "배액: ", input$stent, " · 지연 ", input$stent_delay, "일\n",
      "처방: gem/cis ", p$NCYCLE, " 주기 · durvalumab ", p$DUR_MG, " mg · ",
      "FGFRi ", p$FGI_MG, " mg · ivosidenib ", p$IVO_MG, " mg\n",
      "게이트: ", if (p$GATE_ON > 0) "ON" else "OFF (F1)",
      " · 사전 내성 클론: ", if (p$MU_RES > 0) "있음" else "없음 (F2)",
      " · 기질 장벽: ", if (p$FPEN_MIN < 1) "있음" else "없음 (F4)")
  })

  output$p_overview <- renderPlot({
    d <- sim()
    long(d, c("TTOT", "BILI", "GATE", "SURV")) %>%
      mutate(var = factor(var, c("TTOT", "BILI", "GATE", "SURV"),
                          c("종양 부피 (cm³)", "빌리루빈 (mg/dL)",
                            "화학요법 게이트 (0-1)", "생존확률 S(t)"))) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 1, colour = "#37474F") +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "시간 (일)", y = NULL) + thm
  })

  # ---- 2 --------------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim() %>% mutate(Cgem = GEM_C / 22, Ccis = CIS_C / 18,
                          Cdur = DUR_C / 5.6, Cfgi = FGI_C / 235,
                          Civo = IVO_C / 180)
    long(d, c("Cgem", "DFDCTP", "Ccis", "PTDNA", "Cdur", "Cfgi", "Civo", "COV")) %>%
      mutate(var = factor(var,
        c("Cgem", "DFDCTP", "Ccis", "PTDNA", "Cdur", "Cfgi", "Civo", "COV"),
        c("젬시타빈 혈장 (mg/L)", "세포내 dFdCTP (AU)",
          "유리 백금 (mg/L)", "Pt-DNA 부가체 (AU)",
          "더발루맙 (mg/L)", "FGFR 억제제 (mg/L)",
          "아이보시데닙 (mg/L)", "공유결합 점유율"))) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#0277BD") +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      labs(x = "시간 (일)", y = NULL) + thm
  })

  # ---- 3 --------------------------------------------------------------------
  output$p_clone <- renderPlot({
    d <- sim()
    data.frame(time = d$time,
               Sensitive = d$TS, Persister = d$TP,
               Resistant = d$TR, Metastatic = d$TMET, Total = d$TTOT) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL, name = NULL) +
      scale_y_continuous(trans = "log1p") +
      labs(x = "시간 (일)", y = "부피 (cm³, log1p)",
           title = "클론별 종양 부피") + thm
  })

  output$p_clonefrac <- renderPlot({
    d <- sim()
    data.frame(time = d$time, Persister = d$PERSFR, Resistant = d$RESFR) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(x = "시간 (일)", y = "종양 내 분율",
           title = "지속세포 분율과 내성 클론 분율") + thm
  })

  # ---- 4 --------------------------------------------------------------------
  output$p_gate <- renderPlot({
    d <- sim()
    long(d, c("GBILI", "GANC", "GCRCL", "GPS", "GATE")) %>%
      mutate(var = factor(var, c("GBILI", "GANC", "GCRCL", "GPS", "GATE"),
        c("규칙1 빌리루빈", "규칙2 호중구", "규칙3 신기능(시스플라틴)",
          "규칙4 수행상태", "결합 게이트"))) %>%
      ggplot(aes(time, value, colour = var)) + geom_line(linewidth = 0.9) +
      scale_colour_brewer(palette = "Set1", name = NULL) +
      ylim(0, 1.02) + labs(x = "시간 (일)", y = "전달 분율",
                           title = "규칙 1-4") + thm
  })

  output$p_rdi <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, RDI_PC)) +
      geom_line(linewidth = 1, colour = "#C62828") +
      geom_hline(yintercept = 100, linetype = 2, colour = "#999") +
      ylim(0, 105) +
      labs(x = "시간 (일)", y = "누적 전달 상대 용량 강도 (%)",
           title = "처방된 것이 아니라 실제로 들어간 것") + thm
  })

  output$rdi_txt <- renderText({
    d <- sim(); r <- tail(d$RDI_PC, 1)
    paste0("최종 전달 RDI = ", sprintf("%.1f", r), "%\n",
           "빌리루빈 최고치 = ", sprintf("%.1f", max(d$BILI)), " mg/dL\n",
           "호중구 최저치 = ", sprintf("%.2f", min(d$ANC)), " ×10⁹/L\n",
           "담관염 부담 최고치 = ", sprintf("%.2f", max(d$CHOLI)))
  })

  # ---- 5 --------------------------------------------------------------------
  output$p_biliary <- renderPlot({
    d <- sim()
    long(d, c("OBSR", "PATN", "BILI", "ALP", "CHOLI")) %>%
      mutate(var = factor(var, c("OBSR", "PATN", "BILI", "ALP", "CHOLI"),
        c("담도 폐색 분율", "스텐트 개통성", "총 빌리루빈 (mg/dL)",
          "ALP (U/L)", "담관염 부담"))) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#B71C1C", linewidth = 0.9) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      labs(x = "시간 (일)", y = NULL) + thm
  })

  # ---- 6 --------------------------------------------------------------------
  output$p_liver <- renderPlot({
    d <- sim()
    long(d, c("FLR", "ALB", "ALBI", "ALBIG")) %>%
      mutate(var = factor(var, c("FLR", "ALB", "ALBI", "ALBIG"),
        c("기능적 간 예비능 분율", "알부민 (g/dL)", "ALBI 점수", "ALBI 등급"))) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#558B2F", linewidth = 0.9) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      labs(x = "시간 (일)", y = NULL) + thm
  })

  output$t_albi <- renderDT({
    d <- sim()
    idx <- seq(1, nrow(d), by = max(1, floor(nrow(d) / 20)))
    datatable(d[idx, c("time", "BILI", "ALB", "ALBI", "ALBIG", "FLR", "PS")] %>%
                mutate(across(-time, ~round(., 2))),
              options = list(pageLength = 10), rownames = FALSE)
  })

  # ---- 7 --------------------------------------------------------------------
  output$p_heme <- renderPlot({
    d <- sim()
    long(d, c("ANC", "PLT", "CRCL", "NEURO")) %>%
      mutate(var = factor(var, c("ANC", "PLT", "CRCL", "NEURO"),
        c("호중구 (×10⁹/L)", "혈소판 (×10⁹/L)",
          "크레아티닌 청소율 (mL/min)", "누적 신경병증"))) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#6A1B9A", linewidth = 0.9) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      labs(x = "시간 (일)", y = NULL) + thm
  })

  # ---- 8 --------------------------------------------------------------------
  output$p_immune <- renderPlot({
    d <- sim()
    long(d, c("TEFF", "SUPP", "IL6", "IRAE")) %>%
      mutate(var = factor(var, c("TEFF", "SUPP", "IL6", "IRAE"),
        c("종양내 CD8 효과세포", "억제 구획 (Treg/MDSC/TAM)",
          "IL-6", "면역 관련 이상반응 부담"))) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#AD1457", linewidth = 0.9) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      labs(x = "시간 (일)", y = NULL) + thm
  })

  # ---- 9 --------------------------------------------------------------------
  output$p_bio <- renderPlot({
    d <- sim()
    long(d, c("CA199", "CTDNA", "PHOS", "HG2", "SLD")) %>%
      mutate(var = factor(var, c("CA199", "CTDNA", "PHOS", "HG2", "SLD"),
        c("CA 19-9 (U/mL) — 담즙정체에 교란됨", "ctDNA VAF (AU)",
          "혈청 인산 (mg/dL) — FGFR 표적 결합", "혈장 2-HG — IDH1 표적 결합",
          "RECIST 상대 직경합"))) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#00695C", linewidth = 0.9) +
      facet_wrap(~var, scales = "free_y", ncol = 2) +
      labs(x = "시간 (일)", y = NULL) + thm
  })

  # ---- 10 -------------------------------------------------------------------
  output$p_recist <- renderPlot({
    d <- sim()
    nad <- cummin(d$SLD)
    ggplot(d, aes(time, SLD)) +
      geom_line(linewidth = 1) +
      geom_line(aes(y = nad), linetype = 3, colour = "#0277BD") +
      geom_line(aes(y = 1.2 * nad), linetype = 2, colour = "#C62828") +
      geom_hline(yintercept = 0.70, linetype = 2, colour = "#2E7D32") +
      annotate("text", x = max(d$time) * 0.85, y = 0.72, label = "부분관해 임계",
               colour = "#2E7D32", size = 3.5) +
      labs(x = "시간 (일)", y = "상대 직경합",
           title = "RECIST 1.1 · 최저치 대비 20% 증가가 진행") + thm
  })

  output$endpoint_txt <- renderText({
    d <- sim(); nad <- cummin(d$SLD)
    prog <- which(d$SLD >= 1.2 * nad & d$SLD > 1.05 * min(nad))
    tprog <- if (length(prog)) d$time[prog[1]] else NA
    paste0("최저 상대 직경합 = ", sprintf("%.3f", min(d$SLD)),
           "  (부분관해 = ≤ 0.70)\n",
           "진행 시점 = ", ifelse(is.na(tprog), "관찰 기간 내 없음",
                                  paste0(round(tprog), " 일 (", round(tprog/30.44, 1), " 개월)")), "\n",
           "최종 생존확률 S(t) = ", sprintf("%.3f", tail(d$SURV, 1)), "\n",
           "12개월 S = ", sprintf("%.3f", approx(d$time, d$SURV, 365)$y),
           " · 24개월 S = ", sprintf("%.3f", approx(d$time, d$SURV, 730)$y))
  })

  # ---- 11 -------------------------------------------------------------------
  output$p_hazard <- renderPlot({
    d <- sim()
    data.frame(time = d$time, `h_tumour` = d$HAZT, `h_biliary` = d$HAZB,
               check.names = FALSE) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = c(h_tumour = "#37474F", h_biliary = "#C62828"),
                          name = NULL) +
      labs(x = "시간 (일)", y = "순간 위험 (1/일)",
           title = "서로 다른 시계를 가진 두 개의 위험") + thm
  })

  output$p_surv <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, SURV)) + geom_line(linewidth = 1) +
      geom_line(aes(y = HAZFR), linetype = 2, colour = "#C62828") +
      geom_vline(xintercept = c(365, 730), linetype = 3, colour = "#999") +
      ylim(0, 1) +
      labs(x = "시간 (일)", y = "S(t) (실선) · 담도 위험 분율 (점선)") + thm
  })

  # ---- 12 -------------------------------------------------------------------
  scen <- eventReactive(input$run_scen, {
    base <- list(TT0 = input$tt0, LAM0 = input$lam0, K_2L = input$k2l)
    defs <- list(
      list(n = "S1 배액만",            p = c(base, list(FHILAR = .72)), st = "sems", dl = 0),
      list(n = "S2 gem/cis",           p = c(base, list(FHILAR = .72, GEM_MGM2 = 1000,
                                                        CIS_MGM2 = 25, NCYCLE = 8)), st = "sems", dl = 0),
      list(n = "S3 +durvalumab",       p = c(base, list(FHILAR = .72, GEM_MGM2 = 1000,
                                                        CIS_MGM2 = 25, NCYCLE = 8,
                                                        DUR_MG = 1500, IMMENG = 1)), st = "sems", dl = 0),
      list(n = "S4 배액 60일 지연",    p = c(base, list(FHILAR = .72, GEM_MGM2 = 1000,
                                                        CIS_MGM2 = 25, NCYCLE = 8,
                                                        DUR_MG = 1500)), st = "sems", dl = 60),
      list(n = "S5 플라스틱 스텐트",   p = c(base, list(FHILAR = .72, GEM_MGM2 = 1000,
                                                        CIS_MGM2 = 25, NCYCLE = 8,
                                                        KOCC = 0.0077)), st = "plastic", dl = 0),
      list(n = "S6 페미가티닙",        p = c(base, list(FHILAR = .16, FGFR2 = 1, FGI_MG = 13.5)),
           st = "sems", dl = 0),
      list(n = "S7 푸티바티닙",        p = c(base, list(FHILAR = .16, FGFR2 = 1, FGI_MG = 20,
                                                        FGI_ON = 1, FGI_OFF = 0, COVAL = 1,
                                                        RHO = 0.45)), st = "sems", dl = 0),
      list(n = "S8 아이보시데닙",      p = c(base, list(FHILAR = .16, IDH1 = 1, IVO_MG = 500)),
           st = "sems", dl = 0))
    bind_rows(lapply(defs, function(s) {
      m <- param(mod, s$p)
      if (s$dl > 0) m <- init(m, PATN = 0)
      o <- as.data.frame(mrgsim(m, events = stent_events(s$st, s$dl, input$tend),
                                end = input$tend, delta = 7))
      o$scenario <- s$n
      o
    }))
  })

  output$p_scen <- renderPlot({
    ggplot(scen(), aes(time, SURV, colour = scenario)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = c(365, 730), linetype = 3, colour = "#aaa") +
      labs(x = "시간 (일)", y = "S(t)", colour = NULL,
           title = "시나리오별 생존 곡선 (단일 대표 환자)") + thm
  })

  output$t_scen <- renderDT({
    scen() %>% group_by(scenario) %>%
      summarise(`최저 SLD` = round(min(SLD), 3),
                `12개월 S` = round(approx(time, SURV, 365)$y, 3),
                `24개월 S` = round(approx(time, SURV, 730)$y, 3),
                `최종 RDI %` = round(tail(RDI_PC, 1), 1),
                `최고 빌리루빈` = round(max(BILI), 1),
                `담도 위험 분율(최대)` = round(max(HAZFR), 2), .groups = "drop") %>%
      datatable(options = list(pageLength = 10), rownames = FALSE)
  })

  # ---- 13 -------------------------------------------------------------------
  output$t_fals <- renderTable({
    data.frame(
      스위치 = c("F1 GATE_ON = 0", "F2 MU_RES = 0", "F3 PI_IMMUNE = 1",
                 "F4 FPEN_MIN = 1"),
      예측 = c("배액 지연의 생존 손실이 사라진다",
               "FGFR 억제제 PFS가 길어져야 한다",
               "더발루맙이 중앙 생존을 움직여야 한다",
               "gem/cis 반응률이 과녁을 넘는다"),
      `가상시험 결과` = c(
        "통과 — 지연군 mOS 10.2 → 13.1개월, RDI 68% → 100%",
        "음성 결과(보고됨) — 관찰 기간 내 변화 없음",
        "통과 — 중앙 생존이 시야 밖으로, 24개월 생존 76.9%",
        "통과 — 반응률 31.7% → 86.7%"),
      check.names = FALSE)
  }, striped = TRUE, bordered = TRUE, width = "100%")
}

shinyApp(ui, server)
