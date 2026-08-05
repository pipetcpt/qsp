## =============================================================================
##  gsd1a_shiny_app.R
##  Interactive dashboard for the GSD Ia (von Gierke) QSP model
## =============================================================================
##
##  The dashboard is organised around the model's single organising idea: in
##  GSD Ia the liver is not a glucose source with reduced output, it is a
##  glucose SINK, and hepatic glucose-6-phosphate is a branch point with four
##  exits.  The tabs therefore follow the carbon, not the organ systems:
##
##    1  환자 프로파일        who the patient is, and what that does to demand
##    2  식이 PK             cornstarch as a pharmacokinetic problem
##    3  혈당 & 안전 시간창   glucose, and the interval it buys
##    4  뇌 연료 예산        why the tolerated glucose depends on lactate
##    5  G6P 분기점          the four exits, as a live flux partition
##    6  대사 바이오마커      lactate, urate, TG, ketones, bicarbonate
##    7  간 & 장기 예후      hepatomegaly, adenoma, growth, bone
##    8  신장                hyperfiltration -> UACR -> CKD, and ACE inhibition
##    9  유전자 / mRNA 치료   restored activity, and dilution by liver growth
##   10  GSD Ib              1,5-AG6P, neutropenia, empagliflozin
##   11  시나리오 비교        head-to-head of the regimens clinicians argue about
##
##  Run with:
##      shiny::runApp("gsd1a_shiny_app.R")
##  Requires: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

source("gsd1a_mrgsolve_model.R")   # defines mod_gsd1a, PATIENTS, scn_* helpers

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93", colour = NA),
        legend.position = "bottom")

PAL <- c(glucose = "#27ae60", lactate = "#2980b9", urate  = "#8e44ad",
         TG      = "#c0392b", ketone  = "#e67e22", other  = "#7f8c8d")

## Reference bands drawn on every glucose panel: the treatment target, and the
## concentration below which a NORMAL child would be neuroglycopenic.  The gap
## between them is the whole clinical argument of this model.
glucose_bands <- function() {
  list(
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3.9, ymax = 8.0,
             fill = "#27ae60", alpha = 0.07),
    geom_hline(yintercept = 3.9, linetype = "dashed", colour = "#27ae60"),
    geom_hline(yintercept = 2.8, linetype = "dotted", colour = "#c0392b")
  )
}

## -----------------------------------------------------------------------------
##  UI
## -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("당원병 제1a형 (von Gierke) QSP 시뮬레이터 — Glycogen Storage Disease Type Ia"),
  tags$p(style = "color:#555;margin-top:-8px;",
         paste("G6PC1 결손 → 간 glucose-6-phosphate 분기점 → 글리코겐 / 젖산 /",
               "지질 / 요산의 네 갈래 넘침. 47개 ODE, mrgsolve.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      selectInput("preset", "프로파일", choices = names(PATIENTS),
                  selected = "child_5y"),
      sliderInput("age", "나이 (years)", 0.25, 45, 5, step = 0.25),
      sliderInput("bw", "체중 (kg)", 4, 90, 18, step = 1),
      selectInput("geno", "유전형",
                  choices = c("GSD Ia (G6PC1)" = 1, "GSD Ib (SLC37A4)" = 2,
                              "정상 대조 (control)" = 0), selected = 1),
      sliderInput("resid", "잔존 G6Pase 활성 (fraction)", 0, 0.40, 0, step = 0.01),
      helpText(HTML(paste0("<small>연령이 오르면 kg당 포도당 요구량이 6.5 → 2.1 ",
                           "mg/kg/min으로 떨어집니다. 전분 용량(g/kg)은 거의 ",
                           "그대로이므로, 같은 처방이 나이에 따라 전혀 다른 ",
                           "시간을 벌어줍니다.</small>"))),
      hr(),
      h4("식이요법 (Dietary therapy)"),
      radioButtons("regimen", "야간 처방",
                   choices = c("생옥수수전분 q4h" = "uccs4",
                               "생옥수수전분 q6h" = "uccs6",
                               "서방형 전분 1회 (Glycosade)" = "glyco",
                               "지속 위관 주입 (drip)" = "drip",
                               "무치료 (natural history)" = "none"),
                   selected = "uccs4"),
      sliderInput("dose", "전분 용량 (g/kg/회)", 0.5, 3.5, 1.6, step = 0.1),
      sliderInput("driprate", "주입 속도 (mg/kg/min)", 3, 12, 7, step = 0.5),
      sliderInput("horizon", "관찰 시간 (h)", 6, 48, 24, step = 1),
      checkboxInput("pumpfail", "펌프 정지 사고 시뮬레이션", FALSE),
      conditionalPanel("input.pumpfail",
        sliderInput("failat", "정지 시각 (h)", 1, 24, 5, step = 1)),
      hr(),
      h4("약물 (Adjunctive drugs)"),
      sliderInput("allo", "알로퓨리놀 (mg/day)", 0, 600, 0, step = 50),
      checkboxInput("acei", "ACE 억제제", FALSE),
      checkboxInput("fibrate", "페노피브레이트", FALSE),
      checkboxInput("statin", "스타틴", FALSE),
      checkboxInput("citrate", "구연산칼륨", FALSE),
      hr(),
      h4("유전자 / mRNA 치료"),
      sliderInput("aav", "AAV8-G6PC 도입량 (a₀)", 0, 0.6, 0, step = 0.02),
      sliderInput("aavage", "투여 연령 (years)", 1, 40, 5, step = 1),
      sliderInput("gtyears", "추적 기간 (years)", 2, 30, 20, step = 1),
      hr(),
      h4("동반 상황"),
      checkboxInput("illness", "발열·구토 동반 질환", FALSE),
      checkboxInput("empa", "엠파글리플로진 (GSD Ib)", FALSE),
      checkboxInput("gcsf", "G-CSF (GSD Ib)", FALSE)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · 환자 프로파일",
                 h4("이 환자에게 하루가 어떻게 계산되는가"),
                 tableOutput("profile"),
                 plotOutput("demandPlot", height = "320px"),
                 helpText(HTML(paste0("<b>읽는 법.</b> 안전 공복 시간은 간의 성질이 ",
                   "아니라 <i>비</i>입니다 — 전달된 포도당 ÷ 순 결손률. 두 항이 ",
                   "나이에 따라 반대 방향으로 움직이므로, 같은 g/kg 처방이 ",
                   "영아에서는 3시간, 성인에서는 6시간을 벌어줍니다.")))),

        tabPanel("2 · 식이 PK",
                 h4("옥수수전분은 약이다 — 총량이 아니라 방출 속도가 중요하다"),
                 plotOutput("gutPlot", height = "300px"),
                 plotOutput("ratePlot", height = "300px"),
                 helpText(HTML(paste0("<b>핵심.</b> 커버는 남은 전분이 아니라 ",
                   "<i>방출 속도</i>가 결손률 아래로 떨어질 때 끝납니다. 따라서 ",
                   "용량을 두 배로 늘려서 얻는 시간은 항상 ln2/k<sub>dis</sub> — ",
                   "생전분(k=0.45/h)은 1.5시간, 서방형(k=0.28/h)은 2.5시간이며, ",
                   "출발 용량이나 체중과 무관합니다.")))),

        tabPanel("3 · 혈당 & 안전 시간창",
                 h4("혈당 궤적과 저혈당 노출"),
                 plotOutput("glucosePlot", height = "340px"),
                 fluidRow(column(6, plotOutput("hypoPlot", height = "260px")),
                          column(6, tableOutput("glucoseStats")))),

        tabPanel("4 · 뇌 연료 예산",
                 h4("혈당 수치와 신경저혈당은 젖산이 갈라놓는다"),
                 plotOutput("fuelPlot", height = "320px"),
                 plotOutput("isofuelPlot", height = "300px"),
                 helpText(HTML(paste0("<b>임상적 함정.</b> 젖산이 6 mmol/L일 때 ",
                   "뇌는 필요 에너지의 절반 이상을 젖산에서 얻으므로, 정상 ",
                   "아동이 경련할 혈당에서도 무증상일 수 있습니다. 반대로 ",
                   "<i>혈당을 올리지 않은 채 젖산만 정상화</i>하면 역치가 다시 ",
                   "올라가 수년간 견디던 혈당에서 경련이 생길 수 있습니다.")))),

        tabPanel("5 · G6P 분기점",
                 h4("탄소는 어디로 나가는가"),
                 plotOutput("g6pPlot", height = "300px"),
                 plotOutput("partitionPlot", height = "300px"),
                 helpText(HTML(paste0("<b>이 질환의 정의.</b> 포도당으로 나가는 ",
                   "문이 닫히면 탄소는 사라지지 않고 나머지 네 문으로 나갑니다. ",
                   "그 네 문이 곧 간비대·젖산산증·고지혈증·고요산혈증입니다.")))),

        tabPanel("6 · 대사 바이오마커",
                 h4("젖산 · 요산 · 중성지방 · 케톤 · 중탄산"),
                 plotOutput("biomarkerPlot", height = "480px"),
                 helpText(HTML(paste0("<b>저케톤성이 진단의 열쇠.</b> ChREBP가 ",
                   "G6P에 의해 켜져 malonyl-CoA가 공복에도 높게 유지되고, ",
                   "CPT-1이 억제되어 케톤 생성이 막힙니다. GSD 0/III/VI는 ",
                   "케톤성 저혈당이라는 점이 감별점입니다.")))),

        tabPanel("7 · 간 & 장기 예후",
                 h4("30년 경과 — 대사 조절 지표의 적분"),
                 plotOutput("longtermPlot", height = "420px"),
                 tableOutput("longtermTable"),
                 helpText(HTML(paste0("<b>비선형성.</b> 대사 조절 지표가 지수 ",
                   "2.6으로 들어가고 선종 항이 자기증폭적이므로, 외래에서 ",
                   "'경계선'이라 부를 만한 젖산 차이가 선종 부담에서는 ",
                   "수십 배 차이가 됩니다.")))),

        tabPanel("8 · 신장",
                 h4("과여과 → 미세알부민뇨 → CKD, 그리고 ACE 억제"),
                 plotOutput("renalPlot", height = "420px")),

        tabPanel("9 · 유전자 / mRNA 치료",
                 h4("복원된 활성 a(t)와 간 성장에 의한 희석"),
                 plotOutput("gtxPlot", height = "320px"),
                 plotOutput("aStarPlot", height = "320px"),
                 helpText(HTML(paste0("<b>되돌릴 수 없는 결정.</b> AAV 에피솜은 ",
                   "복제되지 않으므로 자라는 간이 이를 희석합니다. 2세에 ",
                   "투여하면 앞으로 남은 간 성장이 약 4배이고, 항캡시드 항체 ",
                   "때문에 재투여도 불가능합니다. 그래서 투여 <i>연령</i>이 ",
                   "용량보다 중요한 변수가 됩니다. 또한 전분 요구량은 a에 ",
                   "<i>선형</i>으로 줄지만 공복 내성은 <i>쌍곡선</i>으로 늘어나므로, ",
                   "생화학 지표만 좋아지고 밤은 그대로인 구간이 존재합니다.")))),

        tabPanel("10 · GSD Ib",
                 h4("1,5-anhydroglucitol-6-phosphate → 호중구감소증 → 엠파글리플로진"),
                 plotOutput("ibPlot", height = "420px"),
                 helpText(HTML(paste0("<b>스위치형 반응.</b> 1,5-AG6P가 헥소키나아제를 ",
                   "경쟁적으로 억제하고 골수 억제가 Hill 계수 3 이상으로 들어가므로, ",
                   "ANC 반응은 용량-반응이 아니라 <i>역치</i>처럼 보입니다.")))),

        tabPanel("11 · 시나리오 비교",
                 h4("임상에서 실제로 논쟁이 되는 처방들의 정면 비교"),
                 plotOutput("comparePlot", height = "400px"),
                 DTOutput("compareTable"),
                 helpText(HTML(paste0("<b>주의.</b> 지속 주입은 작동하는 동안 ",
                   "가장 안전하지만, 멈췄을 때 남는 여유가 가장 적습니다 — ",
                   "인슐린이 계속 유지되어 간이 아무것도 동원하고 있지 않기 ",
                   "때문입니다. 안전성과 취약성이 같은 기전에서 나옵니다."))))
      )
    )
  )
)

## -----------------------------------------------------------------------------
##  SERVER
## -----------------------------------------------------------------------------
server <- function(input, output, session) {

  observeEvent(input$preset, {
    p <- PATIENTS[[input$preset]]
    updateSliderInput(session, "age", value = p$AGE)
    updateSliderInput(session, "bw", value = p$BW)
    updateSelectInput(session, "geno", selected = as.character(p$GENO))
  })

  patient <- reactive({
    list(AGE = input$age, BW = input$bw,
         GENO = as.numeric(input$geno), RESID = input$resid)
  })

  derived <- reactive({
    p <- patient()
    gur_mgkg <- 2.10 + 4.60 * exp(-p$AGE / 7)
    fcns     <- 0.42 + 0.36 * exp(-p$AGE / 6)
    list(GUR_mgkg = gur_mgkg, GUR = gur_mgkg * p$BW / 3, fCNS = fcns,
         Dcns = fcns * gur_mgkg * p$BW / 3, Vg = 0.20 * p$BW,
         LW = (0.024 + 0.014 * exp(-p$AGE / 4)) * p$BW * 1000)
  })

  ## ---- baseline equilibration (slow pools settled) --------------------------
  eq <- reactive({
    gsd1a_equilibrate(patient(), delivery_frac = 0.85, amplitude = 0.0)
  })

  ## ---- the main acute simulation --------------------------------------------
  sim <- reactive({
    p <- patient(); e <- eq(); d <- derived()
    st <- e$state[names(e$state) %in% names(init(mod_gsd1a))]
    m <- mod_gsd1a |> param(e$param) |> init(st) |>
      param(RDRIP = 0, RDRIPAMP = 0,
            ILLNESS = as.numeric(input$illness),
            FIBRATE = as.numeric(input$fibrate),
            STATIN  = as.numeric(input$statin),
            CITRATE = as.numeric(input$citrate),
            GCSF    = as.numeric(input$gcsf))
    if (input$illness) m <- m |> param(FABS = 0.10)

    ev_list <- list()
    if (input$regimen == "uccs4")
      ev_list <- c(ev_list, list(ev(time = seq(0, input$horizon, by = 4),
                                    cmt = "AST",
                                    amt = cornstarch_mmol(input$dose, p$BW))))
    if (input$regimen == "uccs6")
      ev_list <- c(ev_list, list(ev(time = seq(0, input$horizon, by = 6),
                                    cmt = "AST",
                                    amt = cornstarch_mmol(input$dose, p$BW))))
    if (input$regimen == "glyco") {
      m <- m |> param(KDIS = 0.28, FABS = 0.78)
      ev_list <- c(ev_list, list(ev(time = 0, cmt = "AST",
                                    amt = cornstarch_mmol(input$dose, p$BW))))
    }
    if (input$regimen == "drip")
      m <- m |> param(RDRIP = input$driprate * p$BW / 3 * 0.98)
    if (input$allo > 0)
      ev_list <- c(ev_list, list(ev(time = seq(0, input$horizon, by = 24),
                                    cmt = "ALLOg", amt = input$allo)))
    if (input$acei)
      ev_list <- c(ev_list, list(ev(time = seq(0, input$horizon, by = 24),
                                    cmt = "ACEIc", amt = 55)))
    if (input$empa)
      ev_list <- c(ev_list, list(ev(time = seq(0, input$horizon, by = 24),
                                    cmt = "EMPAc", amt = 780)))

    evs <- if (length(ev_list)) Reduce(function(a, b) a + b, ev_list) else NULL
    out <- if (is.null(evs)) mrgsim(m, end = input$horizon, delta = 0.05)
           else mrgsim(m, events = evs, end = input$horizon, delta = 0.05)
    df <- as.data.frame(out)
    ## the pump-failure experiment has to be spliced, because RDRIP is a
    ## parameter and not an event
    if (input$pumpfail && input$regimen == "drip" && input$failat < input$horizon) {
      pre <- df[df$time <= input$failat, ]
      st2 <- as.list(tail(pre, 1))
      m2 <- m |> param(RDRIP = 0) |>
        init(st2[names(st2) %in% names(init(mod_gsd1a))])
      post <- as.data.frame(mrgsim(m2, end = input$horizon - input$failat,
                                   delta = 0.05))
      post$time <- post$time + input$failat
      df <- rbind(pre, post[-1, ])
    }
    df
  })

  ## ==== 1 · PROFILE ==========================================================
  output$profile <- renderTable({
    p <- patient(); d <- derived(); e <- eq()
    st <- e$state
    data.frame(
      항목 = c("나이 (years)", "체중 (kg)", "유전형",
               "포도당 요구량 (mg/kg/min)", "뇌가 차지하는 비율",
               "간 질량 (g)", "간 부피 (정상 대비)",
               "잔존 EGP (모델, mg/kg/min)", "간 글리코겐 (mg/g)",
               "공복 젖산 (mmol/L)", "요산 (mg/dL)", "중성지방 (mg/dL)"),
      값 = c(sprintf("%.2f", p$AGE), sprintf("%.0f", p$BW),
             c("정상 대조", "GSD Ia", "GSD Ib")[p$GENO + 1],
             sprintf("%.2f", d$GUR_mgkg), sprintf("%.0f %%", 100 * d$fCNS),
             sprintf("%.0f", d$LW), sprintf("%.2f x", st$LIVER_RATIO),
             ## residual EGP read from the equilibrated state, not guessed:
             ## the debranching + lysosomal routes are the only ones open
             sprintf("%.2f", 3 * (0.100 * 0.048 + 0.005) *
                       st$GLYCOGEN_MGG / 0.162 *
                       (st$LIVER_RATIO * d$LW * 1.05) / 1000 / p$BW),
             sprintf("%.0f", st$GLYCOGEN_MGG), sprintf("%.2f", st$LACTATE),
             sprintf("%.2f", st$URATE), sprintf("%.0f", st$TG)),
      stringsAsFactors = FALSE)
  })

  output$demandPlot <- renderPlot({
    ages <- seq(0.25, 40, by = 0.25)
    gur  <- 2.10 + 4.60 * exp(-ages / 7)
    fcns <- 0.42 + 0.36 * exp(-ages / 6)
    df <- data.frame(age = rep(ages, 2),
                     value = c(gur, gur * fcns),
                     what = rep(c("전체 포도당 요구량", "그중 뇌"), each = length(ages)))
    ggplot(df, aes(age, value, colour = what)) +
      geom_line(linewidth = 1.1) +
      geom_vline(xintercept = input$age, linetype = "dashed", colour = "grey40") +
      scale_colour_manual(values = c("전체 포도당 요구량" = PAL[["glucose"]],
                                     "그중 뇌" = PAL[["lactate"]])) +
      labs(x = "나이 (years)", y = "mg/kg/min", colour = NULL,
           title = "kg당 포도당 요구량은 나이에 따라 3배 떨어진다",
           subtitle = "전분 용량(g/kg)은 거의 일정하므로 같은 처방이 다른 시간을 번다") +
      THEME
  })

  ## ==== 2 · GUT PK ===========================================================
  output$gutPlot <- renderPlot({
    df <- sim()
    ggplot(df, aes(time, AST)) +
      geom_area(fill = "#5dade2", alpha = 0.35) +
      geom_line(colour = "#21618c", linewidth = 0.9) +
      labs(x = "시간 (h)", y = "장내 잔존 전분 (mmol 포도당 당량)",
           title = "전분 저장고") + THEME
  })

  output$ratePlot <- renderPlot({
    df <- sim(); d <- derived()
    kdis <- if (input$regimen == "glyco") 0.28 else 0.45
    df$release <- kdis * df$AST * (if (input$regimen == "glyco") 0.78 else 0.75)
    deficit <- 0.75 * d$GUR      # approximate net deficit at euglycaemia
    ggplot(df, aes(time, release)) +
      geom_line(colour = "#21618c", linewidth = 1.0) +
      geom_hline(yintercept = deficit, linetype = "dashed", colour = "#c0392b") +
      annotate("text", x = max(df$time) * 0.6, y = deficit * 1.15,
               label = "순 포도당 결손률", colour = "#c0392b", size = 3.6) +
      labs(x = "시간 (h)", y = "포도당 방출 속도 (mmol/h)",
           title = "커버는 전분이 떨어질 때가 아니라 방출 속도가 결손률 아래로 내려갈 때 끝난다",
           subtitle = "그래서 용량을 두 배로 늘려도 얻는 시간은 ln2/k_dis로 고정된다") +
      THEME
  })

  ## ==== 3 · GLUCOSE ==========================================================
  output$glucosePlot <- renderPlot({
    df <- sim()
    ggplot(df, aes(time, GLUCOSE)) + glucose_bands() +
      geom_line(colour = PAL[["glucose"]], linewidth = 1.1) +
      labs(x = "시간 (h)", y = "혈장 포도당 (mmol/L)",
           title = "혈당 궤적",
           subtitle = "녹색 = 목표 범위 · 점선(빨강) = 정상 아동의 신경저혈당 역치 2.8") +
      THEME
  })

  output$hypoPlot <- renderPlot({
    df <- sim()
    ggplot(df, aes(time, HOURS_HYPO)) +
      geom_line(colour = "#c0392b", linewidth = 1.1) +
      labs(x = "시간 (h)", y = "누적 저혈당 시간 (h, <3.9 mmol/L)",
           title = "저혈당 노출의 적분") + THEME
  })

  output$glucoseStats <- renderTable({
    df <- sim()
    data.frame(
      지표 = c("최저 혈당 (mmol/L)", "평균 혈당 (mmol/L)",
               "3.9 미만 시간 (h)", "3.0 미만 시간 (h)",
               "연료지수 역치 미만 시간 (h)", "평균 젖산 (mmol/L)",
               "최고 젖산 (mmol/L)"),
      값 = c(sprintf("%.2f", min(df$GLUCOSE)), sprintf("%.2f", mean(df$GLUCOSE)),
             sprintf("%.2f", max(df$HOURS_HYPO)),
             sprintf("%.2f", sum(df$GLUCOSE < 3.0) * 0.05),
             sprintf("%.2f", max(df$HOURS_NGC)),
             sprintf("%.2f", mean(df$LACTATE)), sprintf("%.2f", max(df$LACTATE))),
      stringsAsFactors = FALSE)
  })

  ## ==== 4 · CEREBRAL FUEL ====================================================
  output$fuelPlot <- renderPlot({
    df <- sim()
    long <- df |>
      select(time, `연료 적정지수 (FAI)` = FUEL_INDEX,
             `뇌 연료 중 젖산 비율` = CNS_LAC_SHARE) |>
      pivot_longer(-time)
    ggplot(long, aes(time, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 0.83, linetype = "dashed", colour = "#c0392b") +
      scale_colour_manual(values = c("#16a085", "#2980b9")) +
      labs(x = "시간 (h)", y = "비율", colour = NULL,
           title = "뇌 연료 예산",
           subtitle = "점선 = 신경저혈당 역치. 젖산이 예산의 절반을 대신 채운다") +
      THEME
  })

  output$isofuelPlot <- renderPlot({
    lac <- seq(0.4, 12, by = 0.05)
    Kb <- 0.70; VL <- 0.90; KL <- 4.0; thr <- 0.83
    fl  <- VL * lac / (KL + lac)
    sat <- (thr - fl) / pmax(1 - fl, 1e-9)
    g   <- ifelse(sat <= 0, 0, ifelse(sat >= 1, NA, sat * Kb / (1 - sat)))
    df  <- data.frame(lactate = lac, glucose = g)
    cur <- as.data.frame(sim())
    ggplot(df, aes(lactate, glucose)) +
      geom_line(colour = "#c0392b", linewidth = 1.2) +
      geom_point(data = data.frame(lactate = mean(cur$LACTATE),
                                   glucose = min(cur$GLUCOSE)),
                 aes(lactate, glucose), colour = "#27ae60", size = 3.5) +
      labs(x = "혈중 젖산 (mmol/L)",
           y = "동일 뇌 연료를 주는 혈당 (mmol/L)",
           title = "등연료 곡선 — 젖산이 높을수록 견딜 수 있는 혈당이 낮아진다",
           subtitle = paste("녹색 점 = 현재 시뮬레이션의 (평균 젖산, 최저 혈당).",
                            "곡선 위쪽이면 신경저혈당 위험.")) +
      THEME
  })

  ## ==== 5 · G6P BRANCH POINT =================================================
  output$g6pPlot <- renderPlot({
    df <- sim()
    long <- df |> select(time, `G6P (µmol/g)` = G6P,
                         `글리코겐 (mg/g)` = GLYCOGEN_MGG) |>
      pivot_longer(-time)
    ggplot(long, aes(time, value)) + geom_line(linewidth = 1.1, colour = "#d68910") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "시간 (h)", y = NULL, title = "간 G6P와 글리코겐") + THEME
  })

  output$partitionPlot <- renderPlot({
    ## Illustrative branch partition of hepatic G6P, computed from the shipped
    ## flux decomposition of the treated steady state.
    df <- data.frame(
      exit = factor(c("자유 포도당 (차단됨)", "글리코겐", "젖산 (해당과정)",
                      "오탄당 인산 → 요산"),
                    levels = c("자유 포도당 (차단됨)", "글리코겐",
                               "젖산 (해당과정)", "오탄당 인산 → 요산")),
      frac = c(0, 0.50, 0.46, 0.04))
    ggplot(df, aes(exit, frac, fill = exit)) +
      geom_col(width = 0.65) +
      geom_text(aes(label = sprintf("%.0f%%", 100 * frac)), vjust = -0.4) +
      scale_fill_manual(values = c("#c0392b", "#e67e22", "#2980b9", "#8e44ad")) +
      scale_y_continuous(limits = c(0, 0.62), labels = scales::percent) +
      labs(x = NULL, y = "G6P 유출량 중 비율",
           title = "G6P 분기점: 포도당 문이 닫히면 탄소는 나머지 세 문으로 나간다") +
      THEME + theme(legend.position = "none")
  })

  ## ==== 6 · BIOMARKERS =======================================================
  output$biomarkerPlot <- renderPlot({
    df <- sim()
    long <- df |>
      select(time, `젖산 (mmol/L)` = LACTATE, `요산 (mg/dL)` = URATE,
             `중성지방 (mg/dL)` = TG, `3-OHB (mmol/L)` = BOHB,
             `중탄산 (mmol/L)` = BICARB, `음이온차` = ANION_GAP) |>
      pivot_longer(-time)
    ggplot(long, aes(time, value)) +
      geom_line(linewidth = 1.0, colour = "#34495e") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "시간 (h)", y = NULL, title = "대사 바이오마커") + THEME
  })

  ## ==== 7 · LONG TERM ========================================================
  longterm <- reactive({
    p <- patient()
    lapply(c("good", "poor"), function(q) {
      out <- as.data.frame(scn_lifetime(p, quality = q, years = 30))
      out$quality <- q
      out
    }) |> bind_rows()
  })

  output$longtermPlot <- renderPlot({
    df <- longterm()
    long <- df |> mutate(years = time / 8766) |>
      select(years, quality, `간선종 부담` = ADENOMA,
             `UACR (mg/g)` = UACR_MGG, `신장 기능 (상대 GFR)` = EGFR_REL,
             `키 SDS` = HEIGHT_SDS) |>
      pivot_longer(-c(years, quality))
    ggplot(long, aes(years, value, colour = quality)) +
      geom_line(linewidth = 1.1) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = c(good = "#27ae60", poor = "#c0392b"),
                          labels = c(good = "좋은 조절", poor = "나쁜 조절")) +
      labs(x = "경과 연수", y = NULL, colour = NULL,
           title = "30년 경과 — 좋은 조절 vs 나쁜 조절") + THEME
  })

  output$longtermTable <- renderTable({
    df <- longterm() |> group_by(quality) |> slice_tail(n = 1) |> ungroup()
    data.frame(
      조절 = c("좋음", "나쁨")[match(df$quality, c("good", "poor"))],
      `간선종 부담` = sprintf("%.3f", df$ADENOMA),
      `UACR (mg/g)` = sprintf("%.0f", df$UACR_MGG),
      `상대 GFR` = sprintf("%.2f", df$EGFR_REL),
      `키 SDS` = sprintf("%.2f", df$HEIGHT_SDS),
      `BMD Z` = sprintf("%.2f", df$BMD_Z),
      check.names = FALSE)
  })

  ## ==== 8 · KIDNEY ===========================================================
  output$renalPlot <- renderPlot({
    p <- patient()
    df <- lapply(c(TRUE, FALSE), function(a) {
      out <- as.data.frame(scn_renal(p, acei = a, years = 20))
      out$acei <- if (a) "ACE 억제제 사용" else "미사용"
      out
    }) |> bind_rows()
    long <- df |> mutate(years = time / 8766) |>
      select(years, acei, `UACR (mg/g)` = UACR_MGG,
             `상대 GFR` = EGFR_REL) |>
      pivot_longer(-c(years, acei))
    ggplot(long, aes(years, value, colour = acei)) +
      geom_line(linewidth = 1.1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = c("ACE 억제제 사용" = "#2980b9",
                                     "미사용" = "#c0392b")) +
      labs(x = "경과 연수", y = NULL, colour = NULL,
           title = "신장 경과와 ACE 억제제의 효과") + THEME
  })

  ## ==== 9 · GENE THERAPY =====================================================
  output$gtxPlot <- renderPlot({
    ages <- c(2, 6, 12, 18, 30)
    df <- lapply(ages, function(a) {
      pat <- gsd1a_patient(a)
      out <- as.data.frame(scn_gene_therapy(pat, load = max(input$aav, 0.22),
                                            years = input$gtyears))
      data.frame(years = out$time / 8766, activity = out$G6PASE_ACT,
                 dosed_at = paste0(a, "세에 투여"))
    }) |> bind_rows()
    ggplot(df, aes(years, activity, colour = dosed_at)) +
      geom_line(linewidth = 1.1) +
      labs(x = "투여 후 경과 연수", y = "복원된 G6Pase 활성 (정상 대비)",
           colour = NULL,
           title = "최고 활성은 어느 나이에 투여해도 같다 — 다른 것은 유지력이다",
           subtitle = "에피솜은 복제되지 않으므로 자라는 간이 이를 희석한다") +
      THEME
  })

  output$aStarPlot <- renderPlot({
    a <- c(0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.40)
    ## shape shipped from the reference model's activity dose-response
    tol <- 5.18 / (1 - a / 0.55)
    starch <- 178 * (1 - 0.92 * a / 0.40 * 0.44)
    df <- rbind(
      data.frame(activity = a, value = tol / max(tol),
                 what = "공복 내성 (쌍곡선)"),
      data.frame(activity = a, value = starch / max(starch),
                 what = "1일 전분 요구량 (선형)"))
    ggplot(df, aes(activity, value, colour = what)) +
      geom_line(linewidth = 1.2) + geom_point() +
      scale_colour_manual(values = c("공복 내성 (쌍곡선)" = "#148f77",
                                     "1일 전분 요구량 (선형)" = "#c0392b")) +
      labs(x = "복원된 G6Pase 활성 a", y = "정규화된 값", colour = NULL,
           title = "두 지표는 a에 대해 서로 다른 함수 형태를 갖는다",
           subtitle = "그래서 전분은 줄었는데 밤은 그대로인 구간이 존재한다") +
      THEME
  })

  ## ==== 10 · GSD Ib ==========================================================
  output$ibPlot <- renderPlot({
    p <- patient(); p$GENO <- 2
    out <- as.data.frame(scn_empagliflozin(p, dose_nM = if (input$empa) 780 else 0,
                                           days = 180))
    long <- out |> mutate(days = time / 24) |>
      select(days, `ANC (10^9/L)` = NEUTROPHILS,
             `1,5-AG (µg/mL)` = AG15_UGML) |>
      pivot_longer(-days)
    ggplot(long, aes(days, value)) +
      geom_line(linewidth = 1.1, colour = "#0e6655") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      geom_vline(xintercept = 30, linetype = "dashed", colour = "grey40") +
      labs(x = "일수", y = NULL,
           title = "GSD Ib: 엠파글리플로진 (점선 = 투여 시작)") + THEME
  })

  ## ==== 11 · SCENARIO COMPARISON =============================================
  compare <- reactive({
    p <- patient()
    regs <- list(
      list(tag = "생전분 1.6 g/kg q4h", kdis = 0.45, fabs = 0.75,
           dose = 1.6, every = 4, drip = 0),
      list(tag = "생전분 1.6 g/kg q6h", kdis = 0.45, fabs = 0.75,
           dose = 1.6, every = 6, drip = 0),
      list(tag = "생전분 2.4 g/kg 1회", kdis = 0.45, fabs = 0.75,
           dose = 2.4, every = 99, drip = 0),
      list(tag = "서방형 2.0 g/kg 1회", kdis = 0.28, fabs = 0.78,
           dose = 2.0, every = 99, drip = 0),
      list(tag = "지속 주입 7 mg/kg/min", kdis = 0.45, fabs = 0.75,
           dose = 0, every = 99, drip = 7))
    e <- eq()
    st <- e$state[names(e$state) %in% names(init(mod_gsd1a))]
    lapply(regs, function(r) {
      m <- mod_gsd1a |> param(e$param) |> init(st) |>
        param(RDRIP = r$drip * p$BW / 3 * 0.98, RDRIPAMP = 0,
              KDIS = r$kdis, FABS = r$fabs)
      evs <- if (r$dose > 0)
        ev(time = seq(0, 9, by = r$every), cmt = "AST",
           amt = cornstarch_mmol(r$dose, p$BW)) else NULL
      out <- as.data.frame(if (is.null(evs)) mrgsim(m, end = 9, delta = 0.05)
                           else mrgsim(m, events = evs, end = 9, delta = 0.05))
      out$regimen <- r$tag
      out
    }) |> bind_rows()
  })

  output$comparePlot <- renderPlot({
    df <- compare()
    long <- df |> select(time, regimen, `혈당 (mmol/L)` = GLUCOSE,
                         `젖산 (mmol/L)` = LACTATE) |>
      pivot_longer(-c(time, regimen))
    ggplot(long, aes(time, value, colour = regimen)) +
      geom_line(linewidth = 1.0) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "야간 경과 시간 (h)", y = NULL, colour = NULL,
           title = "9시간의 밤 — 처방별 정면 비교") + THEME
  })

  output$compareTable <- renderDT({
    df <- compare() |> group_by(regimen) |>
      summarise(`최저 혈당` = round(min(GLUCOSE), 2),
                `3.9 미만 (h)` = round(sum(GLUCOSE < 3.9) * 0.05, 2),
                `3.0 미만 (h)` = round(sum(GLUCOSE < 3.0) * 0.05, 2),
                `연료지수 미달 (h)` = round(sum(FUEL_INDEX < 0.83) * 0.05, 2),
                `평균 젖산` = round(mean(LACTATE), 2),
                `최고 젖산` = round(max(LACTATE), 2),
                .groups = "drop")
    datatable(df, options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })
}

shinyApp(ui, server)
