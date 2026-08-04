## =====================================================================
##  CONTROLLED OVARIAN STIMULATION — Shiny dashboard
##  조절 난소 자극 QSP 모델 인터랙티브 대시보드
##
##  Run:  R -e "shiny::runApp('cos_shiny_app.R', port = 8080)"
##  Needs: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##         and cos_mrgsolve_model.R in the same directory.
##
##  The app is organised around the model's claim, not around its outputs:
##  tab 5 ("트리거") is the point of the whole thing — it shows the SAME
##  ligand pool being read by three kernels, and lets the user watch the
##  VEGF area collapse by an order of magnitude while the maturation edge
##  stays intact.  Tab 7 shows what that does to the patient.
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("cos_mrgsolve_model.R")     # defines mod, cos_protocol, run_cos, ...

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10))

PAL <- c("#2f6fb5", "#a52a2a", "#3d7a2f", "#b5761f", "#6b3fa0",
         "#22707f", "#95552a", "#96231f")

## ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("조절 난소 자극 (Controlled Ovarian Stimulation) — QSP 대시보드"),
  p(style = "color:#555;",
    strong("이 모델의 명제: "),
    "수확량과 OHSS는 같은 상태변수(성장 과립세포 질량)의 두 가지 읽기이므로 ",
    "FSH 용량으로는 분리되지 않는다. 분리는 오직 트리거의 ",
    strong("모양"), "에서 온다 — 성숙은 신호의 ", strong("앞머리"),
    "만 필요하고 OHSS는 ", strong("면적"), "에 비례한다."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      sliderInput("afc", "AFC — 항문난포수 (2–9 mm)", 3, 36, 12, step = 1),
      sliderInput("t50", "FSH 임계값 중앙값 T50 (IU/L)", 6, 14, 9, step = 0.5),
      sliderInput("tone", "GnRH 박동 긴장도 (PCOS 1.7–1.8)", 0.8, 2.0, 1.0,
                  step = 0.1),
      sliderInput("age", "나이 (년)", 24, 45, 32, step = 1),
      hr(),
      h4("자극 (Stimulation)"),
      selectInput("gonad", "고나도트로핀",
                  c("rFSH 매일 (daily rFSH)" = "rfsh",
                    "corifollitropin alfa 150 µg 단회" = "cori",
                    "hp-hMG (LH 활성 포함)" = "hmg")),
      sliderInput("dose", "1일 용량 (IU)", 50, 450, 150, step = 25),
      sliderInput("fshstart", "자극 시작 (주기 일)", 2, 5, 2, step = 1),
      checkboxInput("coast", "coasting (FSH 중단 후 3일)", FALSE),
      conditionalPanel("input.coast == true",
        sliderInput("coastday", "FSH 중단일 (주기 일)", 6, 12, 8, step = 1)),
      hr(),
      h4("서지 예방 (Surge prevention)"),
      selectInput("supp", "방식",
                  c("GnRH 길항제 (antagonist)" = "ant",
                    "프로게스틴 프라이밍 (PPOS)" = "ppos",
                    "없음 (none)" = "none")),
      conditionalPanel("input.supp == 'ant'",
        sliderInput("antstart", "길항제 시작 (주기 일)", 1, 9, 6, step = 1)),
      hr(),
      h4("트리거 (Trigger)"),
      selectInput("trig", "종류",
                  c("hCG 10 000 IU" = "hcg10", "hCG 5 000 IU" = "hcg5",
                    "hCG 2 500 IU" = "hcg25", "hCG 1 500 IU" = "hcg15",
                    "GnRH 작용제 0.2 mg" = "ago",
                    "이중 트리거 (dual)" = "dual")),
      sliderInput("opu", "트리거→채취 간격 (시간)", 26, 46, 36, step = 1),
      hr(),
      h4("이식 · 보조 (Transfer & adjuvants)"),
      radioButtons("transfer", "이식",
                   c("신선 이식 (fresh)" = "fresh",
                     "전동결 (freeze-all)" = "freeze"), inline = TRUE),
      sliderInput("p4sup", "황체기 프로게스테론 (mg/일)", 0, 800, 600,
                  step = 100),
      checkboxInput("cab", "카베르골린 0.5 mg/일", FALSE),
      checkboxInput("letro", "황체기 레트로졸 5 mg/일", FALSE),
      sliderInput("iv", "지지 수액 (L/일)", 0, 3, 0, step = 0.5),
      hr(),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary",
                   style = "width:100%")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        ## ---------------------------------------------------------- 1
        tabPanel("1. 환자 프로파일",
          fluidRow(
            column(6, plotOutput("p_thresh", height = 320)),
            column(6, plotOutput("p_reserve", height = 320))),
          hr(),
          h4("이 환자에 대해 모델이 계산한 것 (입력이 아니라 출력)"),
          DTOutput("t_patient"),
          helpText("난포 수·자극 일수·수확 난자 수는 어느 것도 입력되지 ",
                   "않습니다. AFC와 임계값 분포만 주어지고 나머지는 적분됩니다.")),
        ## ---------------------------------------------------------- 2
        tabPanel("2. 약동학 (PK)",
          plotOutput("p_pk", height = 620),
          helpText("rFSH: F 0.80 · V/F 26 L · 다회투여 t½ 42 h → 150 IU/일에서 ",
                   "혈중 11.9 IU/L. 길항제: Cmax 11.4 ng/mL, t½ 13 h. ",
                   "hCG 10 000 IU: 최고 220 IU/L, t½ 33 h. ",
                   "작용제(트립토렐린): t½ 3 h.")),
        ## ---------------------------------------------------------- 3
        tabPanel("3. 난포 성장 (folliculometry)",
          plotOutput("p_foll", height = 400),
          fluidRow(
            column(6, plotOutput("p_gran", height = 300)),
            column(6, plotOutput("p_count", height = 300))),
          helpText("각 선은 임계값 분위 슬롯 하나(= AFC/10 난포)입니다. ",
                   "낮은 임계값 슬롯이 먼저 자라고, FSH가 임계값 아래로 ",
                   "떨어진 슬롯은 폐쇄(atresia)됩니다.")),
        ## ---------------------------------------------------------- 4
        tabPanel("4. 내분비 (endocrine)",
          plotOutput("p_endo", height = 640),
          helpText("자연 주기에서는 E2·인히빈이 FSH를 7.8 → 2.6 IU/L로 ",
                   "끌어내려 한 개만 배란합니다. 외인성 FSH는 이 하강을 ",
                   "없애므로 다수 난포가 자랍니다 — 같은 방정식입니다.")),
        ## ---------------------------------------------------------- 5
        tabPanel("5. ★ 트리거: 하나의 리간드, 세 개의 커널",
          fluidRow(
            column(7, plotOutput("p_kern", height = 380)),
            column(5, plotOutput("p_area", height = 380))),
          hr(),
          DTOutput("t_trig"),
          helpText(HTML(paste0(
            "<b>지지 커널</b> K = 2 IU/L (기저 박동으로 충분) · ",
            "<b>성숙 커널</b> K = 25 IU/L, n = 6 (서지의 <i>앞머리</i>만 필요, ",
            "이후 34시간 자율 시계) · <b>VEGF 커널</b> K = 40 IU/L, n = 3 ",
            "(<i>면적</i>에 비례). hCG는 앞머리와 면적을 모두 주고, ",
            "작용제 트리거는 앞머리만 줍니다.")))),
        ## ---------------------------------------------------------- 6
        tabPanel("6. 난자 · 배아 · 임상결과",
          fluidRow(
            column(7, plotOutput("p_chain", height = 380)),
            column(5, plotOutput("p_clbr", height = 380))),
          hr(), DTOutput("t_end"),
          helpText("2PN = MII × 0.72 · 배아낭 = 2PN × (0.55 − 0.011·(나이−33)) ",
                   "· 정상배수체 분율 = 0.85/(1+exp((나이−38.2)/4)) ",
                   "· CLBR = 1 − (1 − 0.52)^정상배수체수")),
        ## ---------------------------------------------------------- 7
        tabPanel("7. OHSS",
          plotOutput("p_ohss", height = 620),
          fluidRow(column(12, DTOutput("t_ohss"))),
          helpText("Golan 등급은 파라미터가 아니라 체액 상태에서 계산됩니다: ",
                   "Hct ≥ 45% 또는 복수 ≥ 1.5 L = 중증, Hct ≥ 55% = 최중증. ",
                   "late OHSS는 임신 hCG가 유일한 원인이므로 전동결이 없앱니다.")),
        ## ---------------------------------------------------------- 8
        tabPanel("8. 시나리오 비교",
          h4("현재 환자에서 트리거·이식 전략만 바꾼 6개 arm"),
          actionButton("runcmp", "비교 실행", class = "btn-success"),
          br(), br(),
          plotOutput("p_cmp", height = 380),
          DTOutput("t_cmp")),
        ## ---------------------------------------------------------- 9
        tabPanel("9. 용량 스윕",
          h4("FSH 용량은 수확량과 OHSS를 분리하는가?"),
          actionButton("runsw", "스윕 실행 (75–450 IU)", class = "btn-success"),
          br(), br(),
          plotOutput("p_sweep", height = 420),
          DTOutput("t_sweep"),
          helpText("모델 결과(AFC 12): 3배 용량 → 난자 +15%, 복수 +125%. ",
                   "난자 수확은 코호트 크기에서 포화하지만 VEGF 면적은 ",
                   "포화하지 않습니다.")),
        ## --------------------------------------------------------- 10
        tabPanel("10. 모니터링 차트",
          h4("초음파·호르몬 모니터링 (임상 차트 형식)"),
          DTOutput("t_chart"),
          helpText("트리거 기준: 17 mm 이상 난포 3개. 이 기준이 자극 일수를 ",
                   "결정하므로 '자극 일수'는 입력이 아니라 예측값입니다.")),
        ## --------------------------------------------------------- 11
        tabPanel("11. 가정과 한계",
          h4("이 모델이 재현하는 것"),
          tags$ul(
            tags$li("무자극 주기의 단일 배란 — 난포 수 파라미터 없이"),
            tags$li("자극 일수 11–12일, 난자 8.0 (AFC 12) / 18.5 (AFC 32)"),
            tags$li("트리거일 P4가 난포 수를 세는 계수기가 되는 현상"),
            tags$li("작용제 트리거: 성숙률 동일, 복수 43배 감소"),
            tags$li("34–38시간 채취 창 (두 시계 사이의 간격)"),
            tags$li("길항제 없으면 자극 6–8일에 조기 LH 서지 → 주기 취소")),
          h4("모델이 명시적으로 알지 못하는 것"),
          tags$ul(
            tags$li("몇 개의 난포가 자라는지, 몇 개의 난자가 나오는지"),
            tags$li("hCG가 작용제보다 '위험하다'는 사실"),
            tags$li("OHSS 중증도 척도·위험점수"),
            tags$li("PCOS가 고반응군이라는 사실")),
          h4("알려진 한계 (검증에서 드러난 것)"),
          tags$ul(
            tags$li("자연 주기 우성난포가 16.9 mm에서 멈춤 (관찰 20–22 mm)"),
            tags$li("작용제 트리거의 황체기 결핍을 과소예측 (−41%, 임상은 더 큼)"),
            tags$li("카베르골린 효과가 메타분석(RR 0.38)보다 보수적 (−17%)"),
            tags$li("황체기 레트로졸은 E2를 79% 낮추지만 복수는 4% 미만 변화 ",
                    "— E2를 표지자로, VEGF를 매개자로 본 구조의 반증 가능한 예측")),
          h4("면책"),
          p(style = "color:#a52a2a;",
            "교육·연구용 반정량 모델입니다. 실제 임상 의사결정·처방·규제 ",
            "제출에 사용해서는 안 됩니다."))
      )
    )
  )
)

## ---------------------------------------------------------------------
build_protocol <- function(input) {
  trig <- switch(input$trig, hcg10 = "hcg", hcg5 = "hcg", hcg25 = "hcg",
                 hcg15 = "hcg", ago = "ago", dual = "dual")
  hd <- switch(input$trig, hcg10 = 10000, hcg5 = 5000, hcg25 = 2500,
               hcg15 = 1500, 10000)
  cos_protocol(
    fsh_dose  = if (input$gonad == "cori") input$dose else input$dose,
    fsh_start = input$fshstart,
    fsh_stop  = if (isTRUE(input$coast)) input$coastday else NA,
    hmg_lh    = if (input$gonad == "hmg") 75 else 0,
    cori      = if (input$gonad == "cori") 150 else 0,
    ant_start = if (input$supp == "ant") input$antstart else 0,
    trigger   = trig, hcg_dose = hd,
    opu_delay = input$opu / 24,
    cab = input$cab, letro = input$letro,
    luteal_p4 = input$p4sup,
    fresh = (input$transfer == "fresh"),
    ppos = (input$supp == "ppos"),
    iv_fluid = input$iv,
    name = "current")
}

server <- function(input, output, session) {

  patient <- reactive(list(AFC = input$afc, T50 = input$t50,
                           TONE = input$tone, AGE = input$age))

  sim <- eventReactive(input$go, {
    withProgress(message = "적분 중 (65 ODE)...", {
      run_cos(mod, patient(), build_protocol(input), tend = 32)
    })
  }, ignoreNULL = FALSE)

  ends <- reactive(cos_endpoints(sim()))

  ## ---- 1 patient ----------------------------------------------------
  output$p_thresh <- renderPlot({
    z <- c(-1.6449, -1.0364, -0.6745, -0.3853, -0.1257,
           0.1257, 0.3853, 0.6745, 1.0364, 1.6449)
    th <- input$t50 * exp(0.45 * z)
    amh <- 0.0838 * (0.55 * input$afc + input$afc * (1 - (5/8)^4/(1 + (5/8)^4))) / 0.6
    th <- th * (1 + 0.15 * (amh / 2.5 - 1))
    d <- data.frame(slot = factor(1:10), th = th, m = input$afc / 10)
    ggplot(d, aes(slot, th)) +
      geom_col(aes(fill = th), width = 0.75) +
      geom_hline(yintercept = 11.9 + 0.3, linetype = 2, colour = "#a52a2a") +
      annotate("text", x = 8.5, y = 12.6, label = "150 IU/일에서의 총 FSH",
               colour = "#a52a2a", size = 3.4) +
      scale_fill_gradient(low = "#c8e6b8", high = "#2f6fb5", guide = "none") +
      labs(title = "FSH 임계값 분포 (10 슬롯)",
           subtitle = sprintf("슬롯당 %.1f 난포 · AMH %.2f ng/mL (계산값)",
                              input$afc / 10, amh),
           x = "슬롯 (민감 → 둔감)", y = "임계값 (IU/L)") + THEME
  })

  output$p_reserve <- renderPlot({
    o <- sim()
    d <- o %>% select(time, AMH, INHB, NF11, NF17) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(title = "난소 예비능 표지자와 난포 수",
           subtitle = "AMH는 자극 중 감소한다 (작아진 소난포 풀)",
           x = "주기 일 (d)", y = NULL) + THEME
  })

  output$t_patient <- renderDT({
    e <- ends()
    datatable(data.frame(
      항목 = c("자극 일수 (예측)", "트리거일 E2 (pg/mL)",
               "트리거일 P4 (ng/mL)", "11 mm 이상 난포", "17 mm 이상 난포",
               "채취 난자", "MII 난자", "성숙률", "배아낭", "정상배수체",
               "누적 생존출생률 (%)", "OHSS 등급 (0–4)"),
      값 = round(c(e$stim_days, e$E2_trig, e$P4_trig, e$foll11, e$foll17,
                   e$oocytes, e$MII, e$MII_frac, e$blast, e$euploid,
                   100 * e$CLBR, e$OHSS), 2)),
      options = list(dom = "t", pageLength = 12), rownames = FALSE)
  })

  ## ---- 2 PK ---------------------------------------------------------
  output$p_pk <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time,
        `외인성 FSH (IU/L)` = FSHTOTAL,
        `길항제 (ng/mL)` = CANTOUT,
        `hCG (IU/L)` = CHCGOUT,
        `작용제 (ng/mL)` = CAGOOUT,
        `LH (IU/L)` = LH,
        `LH 등가 리간드 (IU/L)` = LHEQOUT) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      geom_vline(xintercept = o$ttrig[1], linetype = 3) +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(title = "약동학과 리간드 풀",
           subtitle = "점선 = 트리거 시각", x = "주기 일 (d)", y = NULL) + THEME
  })

  ## ---- 3 follicles --------------------------------------------------
  output$p_foll <- renderPlot({
    o <- sim()
    d <- o %>% select(time, D1:D10) %>% pivot_longer(-time)
    d$name <- factor(d$name, levels = paste0("D", 1:10))
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      geom_hline(yintercept = c(11, 17), linetype = 2, colour = "grey45") +
      geom_vline(xintercept = c(o$ttrig[1], o$topu[1]), linetype = 3) +
      scale_colour_viridis_d(option = "D", end = 0.92, name = "슬롯") +
      labs(title = "난포 직경 (초음파 계측에 대응)",
           subtitle = "점선 수평선 = 11 mm(채취 가능) / 17 mm(트리거 기준)",
           x = "주기 일 (d)", y = "직경 (mm)") + THEME
  })

  output$p_gran <- renderPlot({
    o <- sim()
    d <- o %>% select(time, G1:G10) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, group = name, colour = name)) +
      geom_line(linewidth = 0.7, show.legend = FALSE) +
      scale_colour_viridis_d(option = "D", end = 0.92) +
      labs(title = "과립세포 질량 G_i (성숙난포 = 1)",
           subtitle = "이 질량이 E2·P4·난자·VEGF를 동시에 만든다",
           x = "주기 일 (d)", y = NULL) + THEME
  })

  output$p_count <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time, `MG (총 과립세포 질량)` = MGOUT,
                         `난소 부피 (mL)` = OVOL) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(title = "MG — 모델의 중심 상태변수", x = "주기 일 (d)", y = NULL) +
      THEME
  })

  ## ---- 4 endocrine --------------------------------------------------
  output$p_endo <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time,
        `E2 (pg/mL)` = E2, `P4 (ng/mL)` = P4, `LH (IU/L)` = LH,
        `내인성 FSH (IU/L)` = FSHE, `인히빈 B (pg/mL)` = INHB,
        `AMH (ng/mL)` = AMH, `황체 질량` = CL, `서지 준비도 RS` = RS) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      geom_vline(xintercept = o$ttrig[1], linetype = 3) +
      scale_colour_manual(values = rep(PAL, 2), guide = "none") +
      labs(title = "내분비 프로파일", x = "주기 일 (d)", y = NULL) + THEME
  })

  ## ---- 5 the three kernels -----------------------------------------
  kern <- reactive({
    o <- sim()
    o %>% transmute(time, LHEQ = LHEQOUT,
      `지지 (K 2)` = LHEQOUT / (LHEQOUT + 2),
      `성숙 (K 25, n 6)` = (LHEQOUT / 25)^6 / (1 + (LHEQOUT / 25)^6),
      `VEGF (K 40, n 3)` = (LHEQOUT / 40)^3 / (1 + (LHEQOUT / 40)^3))
  })

  output$p_kern <- renderPlot({
    o <- sim(); d <- kern() %>% select(-LHEQ) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1.0) +
      geom_vline(xintercept = c(o$ttrig[1], o$topu[1]), linetype = 3) +
      coord_cartesian(xlim = c(o$ttrig[1] - 1, o$ttrig[1] + 9)) +
      scale_colour_manual(values = c("#2f6fb5", "#a52a2a", "#96231f"),
                          name = NULL) +
      labs(title = "같은 리간드, 세 개의 커널 점유율",
           subtitle = "성숙은 앞머리만, VEGF는 면적을 가져간다",
           x = "주기 일 (d)", y = "점유율 (0–1)") + THEME
  })

  output$p_area <- renderPlot({
    o <- sim(); tr <- o$ttrig[1]
    a <- data.frame(
      kernel = c("성숙 면적 (d)", "VEGF 면적 (d)"),
      value = c(tail(o$EXPLH, 1) - approx(o$time, o$EXPLH, tr)$y,
                tail(o$EXPV, 1) - approx(o$time, o$EXPV, tr)$y))
    ggplot(a, aes(kernel, value, fill = kernel)) +
      geom_col(width = 0.55, show.legend = FALSE) +
      geom_text(aes(label = sprintf("%.2f d", value)), vjust = -0.4) +
      scale_fill_manual(values = c("#a52a2a", "#96231f")) +
      labs(title = "트리거 이후 누적 노출",
           subtitle = "hCG 10 000 IU: 6.2 / 5.2 · 작용제 0.2 mg: 0.5 / 0.4",
           x = NULL, y = "적분 면적 (일)") + THEME
  })

  output$t_trig <- renderDT({
    o <- sim(); tr <- o$ttrig[1]; e <- ends()
    post <- o$time >= tr
    datatable(data.frame(
      항목 = c("리간드 최고치 (IU/L)", "성숙 커널 면적 (d)",
               "VEGF 커널 면적 (d)", "성숙률 (MII/난자)",
               "최고 Hct (%)", "최대 복수 (L)", "황체기 7일 P4 (ng/mL)"),
      값 = round(c(max(o$LHEQOUT[post]), e$AUC_meiosis, e$AUC_vegf,
                   e$MII_frac, e$HCT_max, e$ASC_max, e$P4_lut7), 2)),
      options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- 6 embryology -------------------------------------------------
  output$p_chain <- renderPlot({
    e <- ends()
    d <- data.frame(
      step = factor(c("채취 난자", "MII", "2PN", "배아낭", "정상배수체"),
                    levels = c("채취 난자", "MII", "2PN", "배아낭",
                               "정상배수체")),
      n = c(e$oocytes, e$MII, e$MII * 0.72, e$blast, e$euploid))
    ggplot(d, aes(step, n, fill = step)) +
      geom_col(width = 0.68, show.legend = FALSE) +
      geom_text(aes(label = sprintf("%.1f", n)), vjust = -0.4) +
      scale_fill_manual(values = PAL) +
      labs(title = "난자 → 배아 사슬 (각 단계는 곱셈)",
           subtitle = sprintf("나이 %d세 · 정상배수체 분율 %.2f",
                              input$age, 0.85/(1+exp((input$age-38.2)/4))),
           x = NULL, y = "개수") + THEME
  })

  output$p_clbr <- renderPlot({
    e <- ends(); n <- seq(0, 8, 0.1)
    d <- data.frame(n = n, clbr = 100 * (1 - 0.48^n))
    ggplot(d, aes(n, clbr)) + geom_line(linewidth = 1) +
      geom_point(data = data.frame(n = e$euploid, clbr = 100 * e$CLBR),
                 colour = "#a52a2a", size = 3.5) +
      labs(title = "누적 생존출생률의 포화",
           subtitle = sprintf("이 환자: 정상배수체 %.2f → CLBR %.1f%%",
                              e$euploid, 100 * e$CLBR),
           x = "정상배수체 배아낭 수", y = "CLBR (%)") + THEME
  })

  output$t_end <- renderDT({
    e <- ends()
    datatable(as.data.frame(lapply(e, function(x)
      if (is.numeric(x)) round(x, 2) else x)),
      options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  ## ---- 7 OHSS -------------------------------------------------------
  output$p_ohss <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time,
        `VEGF 신호` = VEGF, `혈관 투과성 PERM` = PERM,
        `복수 (L)` = ASC, `혈장 용적 (L)` = VP,
        `헤마토크릿 (%)` = HCT, `임신 hCG (IU/L)` = PHCG,
        `난소 부피 (mL)` = OVOL, `Golan 등급` = OHSSG) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      geom_vline(xintercept = c(o$ttrig[1], o$topu[1]), linetype = 3) +
      scale_colour_manual(values = rep(c("#96231f", "#a52a2a", "#b5761f",
                                         "#2f6fb5"), 2), guide = "none") +
      labs(title = "OHSS — 같은 질량의 '면적' 읽기",
           subtitle = "이른 OHSS는 외인성 hCG, 늦은 OHSS는 임신 hCG",
           x = "주기 일 (d)", y = NULL) + THEME
  })

  output$t_ohss <- renderDT({
    o <- sim(); e <- ends()
    g <- c("없음", "경증", "중등증", "중증", "최중증")[e$OHSS + 1]
    datatable(data.frame(
      항목 = c("Golan 등급", "최고 Hct (%)", "최대 복수 (L)",
               "최저 혈장 용적 (L)", "최고 난소 부피 (mL)",
               "VEGF 커널 면적 (d)", "이른 OHSS 최고일 (채취 후 일)",
               "늦은 OHSS 발생"),
      값 = c(g, round(e$HCT_max, 1), round(e$ASC_max, 2),
             round(min(o$VP), 2), round(max(o$OVOL), 0),
             round(e$AUC_vegf, 2),
             round(o$time[which.max(o$ASC)] - o$topu[1], 1),
             ifelse(max(o$PHCG) > 5, "예 (신선 이식 + 임신)", "아니오"))),
      options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- 8 scenario comparison ---------------------------------------
  cmp <- eventReactive(input$runcmp, {
    base <- build_protocol(input)
    arms <- list(
      "hCG 10000 + 신선" = modifyList(base, list(trigger = "hcg",
          hcg_dose = 10000, fresh = TRUE, name = "hCG10-fresh")),
      "hCG 10000 + 전동결" = modifyList(base, list(trigger = "hcg",
          hcg_dose = 10000, fresh = FALSE, name = "hCG10-freeze")),
      "hCG 5000" = modifyList(base, list(trigger = "hcg", hcg_dose = 5000,
          fresh = FALSE, name = "hCG5")),
      "작용제 트리거 + 전동결" = modifyList(base, list(trigger = "ago",
          fresh = FALSE, name = "agonist")),
      "이중 트리거" = modifyList(base, list(trigger = "dual", fresh = FALSE,
          name = "dual")),
      "hCG + 카베르골린" = modifyList(base, list(trigger = "hcg",
          cab = TRUE, fresh = FALSE, name = "hCG-cab")))
    withProgress(message = "6개 arm 적분 중...", {
      bind_rows(lapply(names(arms), function(k) {
        e <- cos_endpoints(run_cos(mod, patient(), arms[[k]], tend = 32))
        e$arm <- k; e }))
    })
  })

  output$p_cmp <- renderPlot({
    d <- cmp() %>% select(arm, oocytes, MII, ASC_max, HCT_max, AUC_vegf,
                          CLBR) %>%
      mutate(CLBR = 100 * CLBR) %>% pivot_longer(-arm)
    ggplot(d, aes(reorder(arm, value), value, fill = name)) +
      geom_col(show.legend = FALSE) + coord_flip() +
      facet_wrap(~name, scales = "free_x", ncol = 3) +
      scale_fill_manual(values = PAL) +
      labs(title = "트리거·이식 전략 비교 (같은 환자, 같은 자극)",
           x = NULL, y = NULL) + THEME
  })

  output$t_cmp <- renderDT({
    datatable(cmp() %>% select(arm, stim_days, E2_trig, P4_trig, oocytes, MII,
                               MII_frac, AUC_vegf, HCT_max, ASC_max, OHSS,
                               euploid, CLBR) %>%
              mutate(across(where(is.numeric), ~round(.x, 2))),
              options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  ## ---- 9 dose sweep -------------------------------------------------
  sw <- eventReactive(input$runsw, {
    base <- build_protocol(input)
    withProgress(message = "용량 스윕 적분 중...", {
      bind_rows(lapply(c(75, 112, 150, 225, 300, 450), function(dd) {
        e <- cos_endpoints(run_cos(mod, patient(),
               modifyList(base, list(fsh_dose = dd,
                                     name = paste0(dd, " IU"))), tend = 32))
        e$dose <- dd; e }))
    })
  })

  output$p_sweep <- renderPlot({
    d <- sw()
    s <- max(d$oocytes) / max(d$ASC_max + 1e-6)
    ggplot(d, aes(dose)) +
      geom_line(aes(y = oocytes, colour = "채취 난자"), linewidth = 1.1) +
      geom_point(aes(y = oocytes, colour = "채취 난자"), size = 2) +
      geom_line(aes(y = ASC_max * s, colour = "최대 복수 (L)"),
                linewidth = 1.1) +
      geom_point(aes(y = ASC_max * s, colour = "최대 복수 (L)"), size = 2) +
      scale_y_continuous(name = "채취 난자 (개)",
        sec.axis = sec_axis(~./s, name = "최대 복수 (L)")) +
      scale_colour_manual(values = c("채취 난자" = "#3d7a2f",
                                     "최대 복수 (L)" = "#96231f"),
                          name = NULL) +
      labs(title = "FSH 용량은 수확량과 위험을 분리하지 못한다",
           subtitle = "난자는 코호트 크기에서 포화, 복수는 포화하지 않는다",
           x = "1일 rFSH 용량 (IU)") + THEME
  })

  output$t_sweep <- renderDT({
    datatable(sw() %>% select(dose, stim_days, E2_trig, oocytes, MII, euploid,
                              AUC_vegf, HCT_max, ASC_max, OHSS, CLBR) %>%
              mutate(across(where(is.numeric), ~round(.x, 2))),
              options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  ## ---- 10 clinic chart ---------------------------------------------
  output$t_chart <- renderDT({
    o <- sim()
    d <- o %>% filter(abs(time - round(time)) < 1e-6, time >= 1) %>%
      transmute(`주기 일` = time,
                `E2 (pg/mL)` = round(E2),
                `P4 (ng/mL)` = round(P4, 2),
                `LH (IU/L)` = round(LH, 2),
                `FSH 총 (IU/L)` = round(FSHTOTAL, 1),
                `≥11 mm` = round(NF11, 1),
                `≥14 mm` = round(NF14, 1),
                `≥17 mm` = round(NF17, 1),
                `최대 직경 (mm)` = round(pmax(D1, D2, D3, D4, D5), 1),
                `Hct (%)` = round(HCT, 1),
                `복수 (L)` = round(ASC, 2))
    datatable(d, options = list(pageLength = 20, scrollX = TRUE),
              rownames = FALSE)
  })
}

shinyApp(ui, server)
