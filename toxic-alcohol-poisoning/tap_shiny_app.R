## ===========================================================================
##  tap_shiny_app.R
##  Toxic alcohol poisoning (methanol / ethylene glycol) — interactive QSP
##  dashboard on top of tap_mrgsolve_model.R
##
##  독성 알코올 중독 (메탄올 · 에틸렌글리콜) QSP 대시보드
##
##  Run:
##      shiny::runApp("tap_shiny_app.R")
##  (tap_mrgsolve_model.R must sit in the same directory)
##
##  ---------------------------------------------------------------------------
##  WHAT THIS APP IS FOR
##  ---------------------------------------------------------------------------
##  Not "look at the curves".  Each tab is built to make ONE argument from the
##  model checkable by hand:
##
##   1  환자·노출        who they are, what they drank, and when they arrived
##   2  모체 알코올 PK   the alcohol is an osmole with a reservoir, not a poison
##   3  산성 대사물      the poison is manufactured, and it lags
##   4  두 간극 시계     OG falls, AG rises, the RATIO reads the clock
##   5  산-염기          bicarbonate is titrated, and the space expands
##   6  이온 트래핑      why pH is in the DELIVERY term, drawn as f_HA
##   7  ADH 차단 산술    how much flux each antidote actually leaves running
##   8  신장·결정        the only threshold in the model, and its rate-dependence
##   9  해독제 PK·투석   the margin C/K that decides whether dialysis un-blocks
##  10  임상 엔드포인트  logMAR, putaminal injury, P(death), P(blindness)
##  11  시나리오 비교    25 pre-built scenarios, side by side
##  12  바이오마커       everything a clinician would chart, in one grid
##
##  DISCLAIMER: educational / research model.  Not for clinical use.
## ===========================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("tap_mrgsolve_model.R", local = TRUE)

MOD <- tap_model()
SCN <- tap_scenarios()

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#eef2f7", colour = NA),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom")

PAL <- c("#2b6ca3", "#b03a3a", "#2e7d5b", "#7a3b8f", "#c98b00", "#555555")

## ---------------------------------------------------------------------------
##  helpers
## ---------------------------------------------------------------------------
long_plot <- function(d, vars, labs = NULL, ylab = "", ncol = 2, logy = FALSE) {
  keep <- intersect(vars, names(d))
  dd <- d %>% select(all_of(c("time", keep))) %>%
    pivot_longer(-time, names_to = "v", values_to = "y")
  if (!is.null(labs)) dd$v <- factor(dd$v, levels = keep, labels = labs[keep])
  p <- ggplot(dd, aes(time, y)) +
    geom_line(colour = PAL[1], linewidth = 0.8) +
    facet_wrap(~ v, scales = "free_y", ncol = ncol) +
    labs(x = "시간 since ingestion (h)", y = ylab) + THEME
  if (logy) p <- p + scale_y_log10()
  p
}

hd_bands <- function(p, s) {
  pr <- s$param
  for (w in list(c(pr$HD1ON, pr$HD1OFF), c(pr$HD2ON, pr$HD2OFF))) {
    if (w[1] < 1e5) p <- p + annotate("rect", xmin = w[1], xmax = w[2],
                                      ymin = -Inf, ymax = Inf,
                                      fill = "#2b6ca3", alpha = 0.10)
  }
  if (pr$CRON < 1e5) p <- p + annotate("rect", xmin = pr$CRON, xmax = pr$CROFF,
                                       ymin = -Inf, ymax = Inf,
                                       fill = "#2e7d5b", alpha = 0.08)
  p
}

## Build a scenario from the sidebar controls.
ui_scen <- function(i) {
  hd <- list(); if (i$use_hd) hd <- list(c(i$hd_start, i$hd_start + i$hd_len))
  scen(
    label = "custom", tend = i$tend, wt = i$wt, gfr = i$gfr,
    meoh_gkg = if (i$agent %in% c("methanol", "both")) i$meoh_dose else 0,
    eg_mL    = if (i$agent %in% c("ethylene glycol", "both")) i$eg_dose else 0,
    etoh_gkg = i$etoh_coing,
    fom_start = if (i$antidote == "fomepizole") i$ant_start else NA,
    fom_boost = i$fom_boost,
    eth_start = if (i$antidote == "ethanol") i$ant_start else NA,
    eth_rate_mgkgh = i$eth_rate, eth_hd_mult = i$eth_hd_mult,
    hd = hd, crrt = if (i$use_crrt) c(i$cr_start, i$cr_start + i$cr_len) else NULL,
    bic = if (i$use_bic) list(start = i$bic_start, stop = i$bic_start + i$bic_len,
                              rate = i$bic_rate) else NULL,
    fol = if (i$use_fol) c(i$ant_start, i$tend) else NULL,
    thi = if (i$use_thi) c(i$ant_start, i$tend) else NULL,
    pyr = if (i$use_thi) c(i$ant_start, i$tend) else NULL,
    ca  = if (i$use_ca) list(start = i$ant_start, stop = i$tend, rate = 3) else NULL,
    noacid = if (i$ketone) 1 else 0)
}

## ===========================================================================
##  UI
## ===========================================================================
ui <- fluidPage(
  titlePanel("독성 알코올 중독 QSP 시뮬레이터 · Toxic Alcohol Poisoning (Methanol / Ethylene Glycol)"),
  tags$p(style = "color:#666;margin-top:-10px",
         HTML("모체 알코올은 삼투질이고 <b>독은 대사물</b>입니다. 이 앱의 모든 탭은 ",
              "<i>injury = ∫(ADH flux)dt</i> 라는 하나의 구조에서 따라 나오는 ",
              "주장 하나씩을 손으로 확인할 수 있게 만들어져 있습니다. ",
              "<b>교육·연구 목적 모델이며 임상 사용 불가.</b>")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      sliderInput("wt", "체중 (kg)", 40, 130, 70, 1),
      sliderInput("gfr", "기저 GFR (정상의 배수)", 0.15, 1.2, 1.0, 0.05),
      sliderInput("tend", "관찰 기간 (h)", 24, 168, 96, 12),

      h4("노출 (Exposure)"),
      selectInput("agent", "섭취 물질",
                  c("methanol", "ethylene glycol", "both"), "methanol"),
      conditionalPanel("input.agent != 'ethylene glycol'",
        sliderInput("meoh_dose", "메탄올 (g/kg)", 0, 2.5, 0.7, 0.05)),
      conditionalPanel("input.agent != 'methanol'",
        sliderInput("eg_dose", "에틸렌글리콜 순수 부피 (mL)", 0, 250, 90, 5)),
      sliderInput("etoh_coing", "동반 섭취 에탄올 (g/kg)", 0, 2, 0, 0.1),
      checkboxInput("ketone", "대조군: 대사물이 케톤 (이소프로판올 등가)", FALSE),

      h4("해독제 (Antidote)"),
      selectInput("antidote", "해독제", c("none", "fomepizole", "ethanol"),
                  "fomepizole"),
      sliderInput("ant_start", "투여 시작 시각 (h)", 0, 36, 8, 1),
      conditionalPanel("input.antidote == 'fomepizole'",
        checkboxInput("fom_boost", "투석 중 q4h 재투여", TRUE)),
      conditionalPanel("input.antidote == 'ethanol'",
        sliderInput("eth_rate", "에탄올 유지속도 (mg/kg/h)", 40, 200, 100, 10),
        sliderInput("eth_hd_mult", "투석 중 속도 배수", 1, 4, 2.5, 0.5)),

      h4("체외 제거 (Extracorporeal)"),
      checkboxInput("use_hd", "간헐적 혈액투석 (IHD)", TRUE),
      conditionalPanel("input.use_hd",
        sliderInput("hd_start", "IHD 시작 (h)", 0, 48, 9, 1),
        sliderInput("hd_len", "IHD 지속 (h)", 2, 12, 6, 1)),
      checkboxInput("use_crrt", "CRRT (40 mL/min)", FALSE),
      conditionalPanel("input.use_crrt",
        sliderInput("cr_start", "CRRT 시작 (h)", 0, 48, 9, 1),
        sliderInput("cr_len", "CRRT 지속 (h)", 6, 72, 30, 2)),

      h4("보조 치료 (Adjuncts)"),
      checkboxInput("use_bic", "중탄산나트륨 정주", TRUE),
      conditionalPanel("input.use_bic",
        sliderInput("bic_rate", "중탄산 속도 (mmol/h)", 0, 80, 30, 5),
        sliderInput("bic_start", "시작 (h)", 0, 36, 8, 1),
        sliderInput("bic_len", "지속 (h)", 4, 48, 18, 2)),
      checkboxInput("use_fol", "폴린산 (leucovorin)", FALSE),
      checkboxInput("use_thi", "티아민 + 피리독신", FALSE),
      checkboxInput("use_ca", "칼슘 보충 (경험적)", FALSE),
      hr(),
      helpText(HTML("음영 = 투석 중. <br>",
                    "실행이 느리면 '관찰 기간'을 줄이세요."))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---------------------------------------------------------------- 1
        tabPanel("1 · 환자·노출",
          br(),
          fluidRow(column(6, h4("섭취량과 초기 분포"), tableOutput("t_expose")),
                   column(6, h4("모델이 예상하는 임상 경과"), tableOutput("t_course"))),
          hr(),
          HTML("<p><b>이 탭에서 확인할 것.</b> 같은 <i>총 섭취량</i>이라도
                <i>도착 시각</i>이 결과를 지배합니다. 아래 표의 '이미 산화된 용량
                비율'은 해독제가 더 이상 막을 것이 없는 정도를 뜻하며, 늦게 온
                환자에게 필요한 것은 해독제가 아니라 투석이라는 결론이
                이 한 칸에서 나옵니다.</p>"),
          plotOutput("p_overview", height = "420px")),

        ## ---------------------------------------------------------------- 2
        tabPanel("2 · 모체 알코올 PK",
          br(), plotOutput("p_parent", height = "520px"),
          hr(),
          HTML("<p><b>ADH를 막으면 독이 사라지는 것이 아니라 저수지가 됩니다.</b>
                메탄올은 신장 배출로가 없어 차단 시 반감기가 40시간대로 늘어나고,
                에틸렌글리콜은 20–30%가 소변으로 그대로 나가므로 차단만으로
                끝납니다 — 신기능이 남아 있는 동안에만. 왼쪽 아래의 말초 구획이
                투석 후 반동(rebound)의 출처입니다.</p>"),
          tableOutput("t_halflife")),

        ## ---------------------------------------------------------------- 3
        tabPanel("3 · 산성 대사물",
          br(), plotOutput("p_acidmet", height = "520px"),
          hr(),
          HTML("<p><b>독은 만들어지는 것이므로 늦게 옵니다.</b> 포름산(또는
                글리콜산)의 정점은 모체 알코올의 정점보다 몇 시간 뒤에 오고,
                동반 섭취 에탄올은 그 간격을 더 벌립니다 — 12–24시간의 무증상
                잠복기는 별개의 현상이 아니라 이 지연의 출력입니다.
                엽산 풀(THF)이 소모되면 포름산 제거 능력 자체가 떨어집니다.</p>")),

        ## ---------------------------------------------------------------- 4
        tabPanel("4 · 두 간극 시계",
          br(), plotOutput("p_gaps", height = "420px"),
          hr(),
          HTML("<p><b>같은 반응을 두 번 읽는 것입니다.</b> 삼투압 간극(OG)은
                모체 알코올이 사라진 만큼 떨어지고, 음이온 간극(AG)은 산이 생긴
                만큼 올라갑니다. 따라서 <b>OG/ΔAG 비는 중증도가 아니라
                경과 시간</b>을 읽습니다. 두 값의 <b>합</b>은 섭취량의 근사
                보존량이지만 정확히 보존되지는 않습니다: OG는 알코올의 분포용적
                (~42 L)으로, AG는 중탄산 공간(HCO3에 따라 36→70 L 이상)으로
                나뉘기 때문입니다. 그 편차의 방향과 크기는 5번 탭에 있습니다.</p>"),
          tableOutput("t_clock")),

        ## ---------------------------------------------------------------- 5
        tabPanel("5 · 산-염기",
          br(), plotOutput("p_ab", height = "520px"),
          hr(),
          HTML("<p><b>중탄산은 붓는 것이 아니라 적정하는 것입니다.</b> 이 모델의
                중탄산 정주는 pH 목표에 따라 스스로 꺼집니다. 또한 외인성
                NaHCO&#8323;는 나트륨도 함께 올리므로 <b>음이온 간극은 중탄산
                치료로 거의 변하지 않습니다</b> — 알칼리는 pH를 고치면서 진단
                단서를 지우지 않습니다. 중탄산 공간이 산증이 깊어질수록 커지는
                것(0.40 + 2.6/[HCO3] L/kg)이 HCO3 하강이 감속하는 이유입니다.</p>"),
          tableOutput("t_space")),

        ## ---------------------------------------------------------------- 6
        tabPanel("6 · 이온 트래핑",
          br(),
          fluidRow(column(7, plotOutput("p_trap", height = "400px")),
                   column(5, h4("f_HA(pH) = 1/(1+10^(pH−pKa))"),
                          tableOutput("t_trap"))),
          hr(),
          HTML("<p><b>pH가 예후 인자인 이유는 pH가 전달항 안에 있기 때문입니다.</b>
                막을 통과하는 것은 비이온화 포름산뿐이고, 그 분율은 pH 7.40 →
                6.80에서 약 4배가 됩니다. 같은 혈중 포름산에서 CNS로 가는
                <i>유입 속도상수</i>가 4배, <i>평형 CNS/혈장 비</i>가 1.8배가
                되므로 뇌 부하는 두세 배 차이가 납니다. 중탄산은 이 두 항을
                동시에 되돌리고, 게다가 알칼리뇨에서 이온화된 포름산이 재흡수되지
                않게 하여 신 청소율까지 올립니다 — 하나의 식, 세 개의 항.</p>"),
          plotOutput("p_cns", height = "360px")),

        ## ---------------------------------------------------------------- 7
        tabPanel("7 · ADH 차단 산술",
          br(),
          h4("경쟁적 차단은 자기 기질에 의해 희석됩니다"),
          tableOutput("t_block"),
          hr(),
          HTML("<p><b>‘포메피졸은 에탄올의 1000배’는 침대 옆에서는 틀린
                말입니다.</b> 경쟁적 억제인수는 (1 + ΣS/Km + I/Ki)/(1 + ΣS/Km)
                이므로 기질이 많을수록 작아집니다. 메탄올 100 mg/dL에서 목표
                농도의 에탄올은 여전히 산 생성 플럭스의 약 18%를 남기고,
                포메피졸 10 µg/mL는 0.6%를 남깁니다 — 정직한 비는 약 30배입니다.
                메탄올 400 mg/dL에서는 포메피졸도 2.4%를 남깁니다.</p>"),
          plotOutput("p_flux", height = "380px")),

        ## ---------------------------------------------------------------- 8
        tabPanel("8 · 신장·결정",
          br(), plotOutput("p_kidney", height = "520px"),
          hr(),
          HTML("<p><b>이 모델의 유일한 역치.</b> 옥살산칼슘은 용해도곱을 넘을
                때만 침전하므로, <b>같은 총량을 천천히 전달하면 손상이
                줄어듭니다</b> — 포름산에는 없는 성질입니다. 그리고 침착된
                결정 1 mmol은 칼슘 1 mmol을 데려가므로, 이온화 칼슘은 결정
                부하의 <b>화학량론적 그림자</b>이고 QT는 그 그림자를
                읽습니다.</p>"),
          tableOutput("t_ss")),

        ## ---------------------------------------------------------------- 9
        tabPanel("9 · 해독제 PK · 투석",
          br(), plotOutput("p_ant", height = "420px"),
          hr(),
          h4("투석이 해독제의 차단을 실제로 풀 수 있는가 — 여유(margin)의 문제"),
          tableOutput("t_margin"),
          HTML("<p>중요한 것은 투석이 해독제를 얼마나 빨리 제거하는가가 아니라
                <b>치료농도와 억제상수 사이의 여유 C/K</b>입니다. 포메피졸은
                약 800배, 에탄올은 약 22배입니다. 6시간 세션은 포메피졸의 여유를
                쓸 수 없고(연속 16시간이 필요), 에탄올의 여유는 대부분
                씁니다. 즉 <b>‘투석 중 해독제 증량’은 포메피졸에서는 안전
                여유이고 에탄올에서는 정량적 필수</b>입니다 — 흔히 가르쳐지는
                것과 반대 방향입니다. 위 그래프에서 <code>FOMKIRAT</code>
                (= C/Ki)를 직접 확인하십시오.</p>")),

        ## ---------------------------------------------------------------- 10
        tabPanel("10 · 임상 엔드포인트",
          br(), plotOutput("p_end", height = "520px"),
          hr(), tableOutput("t_end"),
          HTML("<p><b>시력 손실은 적분이 아니라 역치입니다.</b> 망막은 기저핵과
                같은 이유로(산화적 인산화 실패) 망가지므로, 이 모델의 시신경
                손상 구동력은 <i>유리체 포름산</i>이 아니라 <i>망막 ATP
                결손</i>입니다. 선형 구동력으로 썼을 때는 낮은 농도에 오래
                노출된 모든 환자가 실명했고 그것은 임상과 맞지 않았습니다.</p>")),

        ## ---------------------------------------------------------------- 11
        tabPanel("11 · 시나리오 비교",
          br(),
          checkboxGroupInput("cmp", "비교할 시나리오",
                             choices = names(SCN),
                             selected = c("M1_untreated", "M2_fom_early",
                                          "M3_fom_late", "M4_fom_hd", "M9_full"),
                             inline = TRUE),
          selectInput("cmp_var", "표시할 변수",
                      c("pHart", "HCO3o", "AG", "OG", "OGoverAG", "MEOH", "EG",
                        "FORMmM", "GLYCmM", "FCNS", "FVIT", "FOMug", "ETOH",
                        "LOGMAR", "PDEATH", "PBLIND", "GFRpct", "CAIONo",
                        "OXALuM", "CAOXKID", "FLUXM", "FOMKIRAT", "LACTGAP"),
                      "pHart"),
          plotOutput("p_cmp", height = "420px"),
          hr(), h4("25개 시나리오 요약"), tableOutput("t_all")),

        ## ---------------------------------------------------------------- 12
        tabPanel("12 · 바이오마커 대시보드",
          br(), plotOutput("p_bio", height = "760px"),
          hr(),
          HTML("<p><b>젖산 간극(lactate gap)</b>은 공짜 단서입니다. 글리콜산이
                일부 혈액가스 젖산 전극(lactate oxidase)에 교차반응하므로
                POC 젖산이 검사실 젖산보다 훨씬 높게 나옵니다. 이 모델은 그
                인공물을 <code>LACTPOC</code>로 계산해 둡니다 — 에틸렌글리콜을
                의심할 근거이지 젖산산증을 치료할 근거가 아닙니다.</p>"))
      )
    )
  )
)

## ===========================================================================
##  SERVER
## ===========================================================================
server <- function(input, output, session) {

  cur <- reactive({
    s <- ui_scen(input)
    list(s = s, d = tap_run(MOD, s, dt = 0.1))
  })

  ## ------------------------------------------------------------------- 1 --
  output$t_expose <- renderTable({
    i <- input
    data.frame(
      항목 = c("체중 (kg)", "기저 GFR (배수)", "메탄올 (mmol)",
               "에틸렌글리콜 (mmol)", "동반 에탄올 (mmol)",
               "메탄올 최대 예상 산 등가 (mmol)", "에틸렌글리콜 산 등가 (mmol)"),
      값 = c(i$wt, i$gfr,
             round(if (i$agent != "ethylene glycol") mmol_meoh(i$meoh_dose, i$wt) else 0),
             round(if (i$agent != "methanol") mmol_eg(i$eg_dose) else 0),
             round(mmol_etoh(i$etoh_coing, i$wt)),
             round(if (i$agent != "ethylene glycol") mmol_meoh(i$meoh_dose, i$wt) else 0),
             round(1.25 * if (i$agent != "methanol") mmol_eg(i$eg_dose) else 0)))
  }, digits = 2)

  output$t_course <- renderTable({
    d <- cur()$d; fin <- d[nrow(d), ]
    dose <- max(1e-9, mmol_meoh(input$meoh_dose, input$wt) +
                      mmol_eg(if (input$agent != "methanol") input$eg_dose else 0))
    data.frame(
      지표 = c("pH 최저", "HCO3 최저 (mM)", "AG 최고", "OG 최고",
               "포름산 최고 (mM)", "글리콜산 최고 (mM)",
               "이미 산화된 용량 비율 (%)", "GFR 최저 (%)",
               "최종 logMAR", "P(사망)", "P(영구 시각손상)"),
      값 = c(round(min(d$pHart), 3), round(min(d$HCO3o), 1),
             round(max(d$AG), 1), round(max(d$OG), 1),
             round(max(d$FORMmM), 2), round(max(d$GLYCmM), 2),
             round(100 * (fin$MEOHOX + fin$EGOX) / dose, 1),
             round(min(d$GFRpct)), round(fin$LOGMAR, 2),
             round(fin$PDEATH, 3), round(fin$PBLIND, 3)))
  }, digits = 3)

  output$p_overview <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("MEOH", "EG", "FORMmM", "GLYCmM", "pHart", "AG"),
                       labs = c(MEOH = "메탄올 (mg/dL)", EG = "EG (mg/dL)",
                                FORMmM = "포름산 (mM)", GLYCmM = "글리콜산 (mM)",
                                pHart = "동맥 pH", AG = "음이온 간극"),
                       ncol = 3), x$s)
  })

  ## ------------------------------------------------------------------- 2 --
  output$p_parent <- renderPlot({
    x <- cur(); d <- x$d
    d$MEOH_periph <- d$A_MEOH2 / (0.25 * input$wt) * 32.042 / 10
    d$EG_periph   <- d$A_EG2   / (0.27 * input$wt) * 62.068 / 10
    hd_bands(long_plot(d, c("MEOH", "MEOH_periph", "EG", "EG_periph",
                            "ETOH", "OG"),
                       labs = c(MEOH = "메탄올 중앙 (mg/dL)",
                                MEOH_periph = "메탄올 말초 (mg/dL) — 반동의 출처",
                                EG = "EG 중앙 (mg/dL)",
                                EG_periph = "EG 말초 (mg/dL)",
                                ETOH = "에탄올 (mg/dL)",
                                OG = "삼투압 간극 (mOsm)"),
                       ncol = 2), x$s)
  })

  output$t_halflife <- renderTable({
    d <- cur()$d
    f <- function(v) {
      i0 <- which.max(v); m <- seq_along(v) > i0 + 40 & v > 2
      if (sum(m) < 10) return(NA_real_)
      -log(2) / coef(lm(log(v[m]) ~ d$time[m]))[2]
    }
    data.frame(
      물질 = c("메탄올", "에틸렌글리콜"),
      `정점 (mg/dL)` = c(round(max(d$MEOH), 1), round(max(d$EG), 1)),
      `관찰된 소실 반감기 (h)` = round(c(f(d$MEOH), f(d$EG)), 1),
      `<20 mg/dL 도달 (h)` = c(
        suppressWarnings(min(d$time[d$MEOH < 20 & d$time > d$time[which.max(d$MEOH)]])),
        suppressWarnings(min(d$time[d$EG   < 20 & d$time > d$time[which.max(d$EG)]]))),
      check.names = FALSE)
  }, digits = 1)

  ## ------------------------------------------------------------------- 3 --
  output$p_acidmet <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("FORMmM", "FORMmgL", "FCNS", "FVIT",
                              "GLYCmM", "THFo"),
                       labs = c(FORMmM = "혈장 포름산 (mM)",
                                FORMmgL = "혈장 포름산 (mg/L)",
                                FCNS = "CNS 포름산 (mM) — 손상시키는 양",
                                FVIT = "유리체·망막 포름산 (mM)",
                                GLYCmM = "혈장 글리콜산 (mM)",
                                THFo = "간 THF 풀 (상대값)"),
                       ncol = 2), x$s)
  })

  ## ------------------------------------------------------------------- 4 --
  output$p_gaps <- renderPlot({
    x <- cur(); d <- x$d
    dd <- d %>% select(time, OG, DAG, GAPSUM, OGoverAG) %>%
      pivot_longer(-time, names_to = "v", values_to = "y")
    dd$v <- factor(dd$v, levels = c("OG", "DAG", "GAPSUM", "OGoverAG"),
                   labels = c("삼투압 간극 OG (mOsm)",
                              "ΔAG (mEq) — 같은 플럭스, 반대 부호",
                              "OG + ΔAG — 근사 보존량 (섭취량 추정)",
                              "OG/ΔAG — 시계 (경과 시간)"))
    p <- ggplot(dd, aes(time, y)) +
      geom_line(colour = PAL[2], linewidth = 0.9) +
      facet_wrap(~ v, scales = "free_y", ncol = 2) +
      labs(x = "시간 (h)", y = NULL) + THEME
    hd_bands(p, x$s)
  })

  output$t_clock <- renderTable({
    d <- cur()$d
    tt <- c(0.5, 1, 2, 4, 6, 8, 12, 16, 20, 24, 36, 48)
    tt <- tt[tt <= max(d$time)]
    i <- vapply(tt, function(t) which.min(abs(d$time - t)), 1L)
    out <- d[i, c("time", "MEOH", "EG", "FORMmM", "GLYCmM", "OG", "DAG",
                  "OGoverAG", "GAPSUM", "pHart")]
    names(out) <- c("t (h)", "MeOH", "EG", "포름산", "글리콜산",
                    "OG", "ΔAG", "OG/ΔAG (시계)", "합", "pH")
    out
  }, digits = 2)

  ## ------------------------------------------------------------------- 5 --
  output$p_ab <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("pHart", "HCO3o", "PACO2o", "AG", "NAo", "LACTo"),
                       labs = c(pHart = "동맥 pH", HCO3o = "HCO3 (mM)",
                                PACO2o = "PaCO2 (mmHg)", AG = "음이온 간극",
                                NAo = "나트륨 (mM) — 알칼리 부하의 대가",
                                LACTo = "젖산 (mM)"),
                       ncol = 3), x$s)
  })

  output$t_space <- renderTable({
    h <- c(24, 20, 16, 12, 10, 8, 6, 4)
    data.frame(`HCO3 (mM)` = h,
               `중탄산 공간 (L)` = round((0.40 + 2.6 / h) * input$wt, 1),
               `공간 / Vd(메탄올)` = round((0.40 + 2.6 / h) * input$wt /
                                            (0.60 * input$wt), 3),
               check.names = FALSE)
  }, digits = 2)

  ## ------------------------------------------------------------------- 6 --
  output$p_trap <- renderPlot({
    tb <- tap_trapping_table()
    dd <- tb %>% select(pH, rel_entry, eq_CNS_over_plasma) %>%
      pivot_longer(-pH, names_to = "v", values_to = "y")
    dd$v <- factor(dd$v, levels = c("rel_entry", "eq_CNS_over_plasma"),
                   labels = c("CNS 유입 속도상수 (pH 7.40 = 1)",
                              "평형 CNS/혈장 포름산 비"))
    ggplot(dd, aes(pH, y, colour = v)) +
      geom_line(linewidth = 1) + geom_point() +
      scale_x_reverse() + scale_colour_manual(values = PAL[1:2], name = NULL) +
      labs(x = "동맥 pH (오른쪽이 산성)", y = "상대값") + THEME
  })

  output$t_trap <- renderTable({
    tb <- tap_trapping_table()
    names(tb) <- c("pH", "f_HA 혈장", "유입 상대속도", "f_HA 뇌", "평형 CNS/혈장")
    tb
  }, digits = 4)

  output$p_cns <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("fHAratio", "FCNS", "FVIT", "ATPo"),
                       labs = c(fHAratio = "f_HA / f_HA(7.40) — 전달 배수",
                                FCNS = "CNS 포름산 (mM)",
                                FVIT = "유리체 포름산 (mM)",
                                ATPo = "CNS ATP (1 = 정상)"),
                       ncol = 4), x$s)
  })

  ## ------------------------------------------------------------------- 7 --
  output$t_block <- renderTable({
    tb <- tap_blockade_table(MOD)
    names(tb) <- c("메탄올 (mg/dL)", "에탄올 (mg/dL)", "포메피졸 (µg/mL)",
                   "억제 인수", "남은 플럭스 (%)", "플럭스 (mmol/h)")
    tb
  }, digits = 3)

  output$p_flux <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("FLUXM", "FLUXE", "INHFAC", "FLUXLEFT"),
                       labs = c(FLUXM = "메탄올 산화 플럭스 (mmol/h)",
                                FLUXE = "EG 산화 플럭스 (mmol/h)",
                                INHFAC = "억제 인수 (배)",
                                FLUXLEFT = "남은 플럭스 (%)"),
                       ncol = 2), x$s)
  })

  ## ------------------------------------------------------------------- 8 --
  output$p_kidney <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("OXALuM", "SSTUB", "CAOXKID", "PTINJo",
                              "GFRpct", "CAIONo", "QTCms", "SSPLAS"),
                       labs = c(OXALuM = "혈장 옥살산 (µM)",
                                SSTUB = "세뇨관 과포화도 (Ksp 배수)",
                                CAOXKID = "신 CaOx 침착 (mmol, 누적)",
                                PTINJo = "근위세뇨관 손상 (0–1)",
                                GFRpct = "GFR (기저의 %)",
                                CAIONo = "이온화 칼슘 (mM)",
                                QTCms = "QTc (ms)",
                                SSPLAS = "혈장 과포화도"),
                       ncol = 4), x$s)
  })

  output$t_ss <- renderTable({
    ox <- c(2, 5, 10, 20, 40, 80, 160)
    ksp <- 2.32e-3
    data.frame(`옥살산 (µM)` = ox,
               `혈장 SS` = round(1.10 * ox / 1000 / ksp, 2),
               `세뇨관 SS (×50)` = round(1.10 * ox / 1000 / ksp * 50, 0),
               판정 = ifelse(1.10 * ox / 1000 / ksp > 4, "세뇨관 + 조직",
                       ifelse(1.10 * ox / 1000 / ksp * 50 > 8, "세뇨관", "없음")),
               check.names = FALSE)
  }, digits = 2)

  ## ------------------------------------------------------------------- 9 --
  output$p_ant <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("FOMug", "FOMKIRAT", "ETOH", "FLUXM"),
                       labs = c(FOMug = "포메피졸 (µg/mL)",
                                FOMKIRAT = "포메피졸 C/Ki (여유)",
                                ETOH = "에탄올 (mg/dL)",
                                FLUXM = "산 생성 플럭스 (mmol/h)"),
                       ncol = 2, logy = FALSE), x$s)
  })

  output$t_margin <- renderTable({
    tb <- tap_margin_table(MOD)
    names(tb) <- c("해독제", "치료농도 (mM)", "억제상수 (mM)", "투석 CL (L/h)",
                   "중앙 용적 (L)", "여유 C/K (배)", "투석 k (1/h)",
                   "여유를 다 쓰는 데 필요한 투석 시간 (h)")
    tb
  }, digits = 4)

  ## ------------------------------------------------------------------ 10 --
  output$p_end <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d, c("LOGMAR", "OPTICo", "PUTo", "ATPo",
                              "COMA", "PDEATH", "PBLIND", "GFRpct"),
                       labs = c(LOGMAR = "시력 logMAR (0 = 20/20)",
                                OPTICo = "시신경 손상 지수",
                                PUTo = "기저핵(putamen) 손상 지수",
                                ATPo = "CNS ATP",
                                COMA = "의식 저하 지수",
                                PDEATH = "누적 P(사망)",
                                PBLIND = "누적 P(영구 시각손상)",
                                GFRpct = "GFR (%)"),
                       ncol = 4), x$s)
  })

  output$t_end <- renderTable({
    d <- cur()$d; fin <- d[nrow(d), ]
    data.frame(
      엔드포인트 = c("최종 logMAR", "logMAR 해석", "기저핵 손상 지수",
                     "P(사망)", "P(영구 시각손상)", "GFR 최저 (%)",
                     "생존 불가 pH 도달", "관찰 종료 시각 (h)"),
      값 = c(sprintf("%.2f", fin$LOGMAR),
             if (fin$LOGMAR < 0.1) "정상에 가까움" else
             if (fin$LOGMAR < 0.5) "경도 저하" else
             if (fin$LOGMAR < 1.0) "중등도 저하" else
             if (fin$LOGMAR < 1.6) "법적 실명 범위" else "광각 소실 수준",
             sprintf("%.3f", fin$PUTo), sprintf("%.3f", fin$PDEATH),
             sprintf("%.3f", fin$PBLIND), sprintf("%.0f", min(d$GFRpct)),
             if (any(d$FATAL > 0)) "예" else "아니오",
             sprintf("%.1f", fin$time)))
  })

  ## ------------------------------------------------------------------ 11 --
  cmp_data <- reactive({
    req(length(input$cmp) > 0)
    do.call(rbind, lapply(input$cmp, function(k) {
      d <- tap_run(MOD, SCN[[k]], dt = 0.2)
      d$scenario <- k
      d
    }))
  })

  output$p_cmp <- renderPlot({
    d <- cmp_data(); v <- input$cmp_var
    ggplot(d, aes(time, .data[[v]], colour = scenario)) +
      geom_line(linewidth = 0.85) +
      labs(x = "시간 (h)", y = v, colour = NULL) +
      scale_colour_manual(values = rep(PAL, 5)) + THEME
  })

  output$t_all <- renderTable({
    a <- tap_all(MOD, dt = 0.25)
    a$label <- NULL
    a
  }, digits = 3)

  ## ------------------------------------------------------------------ 12 --
  output$p_bio <- renderPlot({
    x <- cur()
    hd_bands(long_plot(x$d,
      c("MEOH", "EG", "FORMmM", "GLYCmM", "OXALuM", "FOMug", "ETOH",
        "pHart", "HCO3o", "PACO2o", "AG", "OG", "OGoverAG", "NAo",
        "CAIONo", "QTCms", "LACTo", "LACTPOC", "LACTGAP", "pHur",
        "GFRpct", "CAOXKID", "ATPo", "LOGMAR"),
      labs = c(MEOH = "메탄올 mg/dL", EG = "EG mg/dL",
               FORMmM = "포름산 mM", GLYCmM = "글리콜산 mM",
               OXALuM = "옥살산 µM", FOMug = "포메피졸 µg/mL",
               ETOH = "에탄올 mg/dL", pHart = "pH", HCO3o = "HCO3",
               PACO2o = "PaCO2", AG = "AG", OG = "OG",
               OGoverAG = "OG/ΔAG 시계", NAo = "Na",
               CAIONo = "이온화 Ca", QTCms = "QTc",
               LACTo = "젖산 (검사실)", LACTPOC = "젖산 (POC, 간섭 포함)",
               LACTGAP = "젖산 간극", pHur = "소변 pH",
               GFRpct = "GFR %", CAOXKID = "신 CaOx mmol",
               ATPo = "CNS ATP", LOGMAR = "logMAR"),
      ncol = 4), x$s)
  })
}

shinyApp(ui, server)
