## =====================================================================
##  icp_shiny_app.R
##  Intrahepatic Cholestasis of Pregnancy (ICP) — QSP 대시보드
##  10 tabs over icp_mrgsolve_model.R
##
##  DESIGN NOTE — why this app is laid out the way it is
##  ---------------------------------------------------
##  The app is built to make ONE thing hard to miss: the bile acid axis
##  and the itch axis are different systems, and the assay the clinic
##  uses reports neither of them cleanly.  So:
##    * Tab 2 always plots the assay value AND the endogenous fraction on
##      the same axes, because the gap between them is the drug effect
##      nobody attributes to the drug.
##    * Tab 4 (fetal) and Tab 6 (itch) can be driven from the same run
##      and deliberately show opposite responses to the same regimen.
##    * Tab 7 refuses to report a single "risk": it decomposes the
##      PITCHES composite so you can see which component moved.
##    * Tab 8 plots the delivery decision as two curves crossing, not as
##      a threshold, because that is what it is.
##
##  USAGE
##      install.packages(c("shiny","mrgsolve","dplyr","tidyr","ggplot2",
##                         "DT","patchwork"))
##      shiny::runApp("icp_shiny_app.R")
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

MODFILE <- "icp_mrgsolve_model"
mod <- mread_cache(MODFILE, ".")

MW_UDCA <- 392.57
ga2d <- function(ga) (ga - 20) * 7

## ---- susceptibility vectors -----------------------------------------
GENOTYPES <- list(
  "정상 임신 (wild type)"            = list(GBSEP = 1.000, GMDR3 = 1.000, GSULT = 1.00),
  "ICP · TBA <40 (경증, ~16)"        = list(GBSEP = 0.795, GMDR3 = 0.695, GSULT = 1.72),
  "ICP · TBA <40 (상단, ~26)"        = list(GBSEP = 0.735, GMDR3 = 0.620, GSULT = 1.95),
  "ICP · TBA 40-99 (중증, ~58)"      = list(GBSEP = 0.650, GMDR3 = 0.500, GSULT = 2.40),
  "ICP · TBA >=100 (최중증, ~138)"   = list(GBSEP = 0.520, GMDR3 = 0.350, GSULT = 3.10),
  "ABCB11 V444A 동형접합"            = list(GBSEP = 0.550, GMDR3 = 0.950, GSULT = 2.20),
  "ABCB4 중증 변이 (MDR3 30%)"       = list(GBSEP = 0.920, GMDR3 = 0.300, GSULT = 2.20),
  "ATP8B1 변이"                      = list(GFIC1 = 0.680, GMDR3 = 0.800, GSULT = 2.40),
  "NR1H4(FXR) 기능저하"              = list(GFXR = 0.550, GBSEP = 0.650, GMDR3 = 0.500, GSULT = 2.40),
  "태반 수송능 저하 (GP 0.6)"        = list(GBSEP = 0.650, GMDR3 = 0.500, GSULT = 2.40, GP = 0.60)
)

## ---- dosing regimen builders ----------------------------------------
build_events <- function(udca_mg, udca_from, rif_mg, rif_from, chol_g, chol_from,
                         same_mg, ibat_umol, ntx_mg, ah_mg, bet_ga, del_ga) {
  e <- NULL
  add <- function(x) if (is.null(e)) x else c(e, x)
  if (udca_mg > 0 && udca_from < del_ga) {
    n <- max(1, round((del_ga - udca_from) * 14))
    e <- add(ev(amt = udca_mg / 2 / MW_UDCA * 1000, cmt = "UDDEP",
                ii = 0.5, addl = n - 1, time = ga2d(udca_from)))
  }
  if (rif_mg > 0 && rif_from < del_ga) {
    n <- max(1, round((del_ga - rif_from) * 14))
    e <- add(ev(amt = rif_mg / 2, cmt = "RIFDEP", ii = 0.5, addl = n - 1,
                time = ga2d(rif_from)))
  }
  if (chol_g > 0 && chol_from < del_ga) {
    n <- max(1, round((del_ga - chol_from) * 28))
    e <- add(ev(amt = chol_g / 4, cmt = "CHOLL", ii = 0.25, addl = n - 1,
                time = ga2d(chol_from)))
  }
  if (same_mg > 0) {
    n <- max(1, round(del_ga - 32))
    e <- add(ev(amt = same_mg, cmt = "SAMC", ii = 1, addl = n - 1,
                time = ga2d(32)))
  }
  if (ibat_umol > 0) {
    n <- max(1, round(del_ga - 32))
    e <- add(ev(amt = ibat_umol, cmt = "ODEV", ii = 1, addl = n - 1,
                time = ga2d(32)))
  }
  if (ntx_mg > 0) {
    n <- max(1, round(del_ga - 32))
    e <- add(ev(amt = ntx_mg, cmt = "NTXC", ii = 1, addl = n - 1,
                time = ga2d(32)))
  }
  if (ah_mg > 0) {
    n <- max(1, round((del_ga - 32) * 3))
    e <- add(ev(amt = ah_mg / 3, cmt = "AHC", ii = 1/3, addl = n - 1,
                time = ga2d(32)))
  }
  if (!is.na(bet_ga) && bet_ga > 0) {
    e <- add(ev(amt = 12, cmt = "BETC", ii = 1, addl = 1, time = ga2d(bet_ga)))
  }
  e
}

run_case <- function(gen, ev_obj, del_ga, twin = 1, extra = list()) {
  p <- c(gen, extra, list(TDEL = ga2d(del_ga), TWIN = twin))
  m <- mod %>% param(p)
  out <- if (is.null(ev_obj)) {
    m %>% mrgsim(end = ga2d(del_ga) + 21, delta = 0.5)
  } else {
    m %>% mrgsim(events = ev_obj, end = ga2d(del_ga) + 21, delta = 0.5)
  }
  as_tibble(out) %>% mutate(GA = 20 + time / 7,
                            phase = ifelse(time < ga2d(del_ga),
                                           "임신 중", "분만 후"))
}

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, colour = "grey35"))

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel("임신성 간내 담즙정체 (ICP) — QSP 대시보드"),
  p(style = "color:#555;",
    strong("교육·연구 목적의 모델입니다."),
    " 임상 의사결정에 사용하지 마십시오. 두 개의 인과 사슬(담즙산 축 · 가려움 축)이",
    " 서로 다른 약에 반응한다는 것이 이 모델의 핵심 주장이며, 탭 4 와 탭 6 을",
    " 같은 처방으로 나란히 보면 그 분리가 드러납니다."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 프로파일"),
      selectInput("gen", "유전 소인 (susceptibility vector)",
                  choices = names(GENOTYPES), selected = names(GENOTYPES)[4]),
      checkboxInput("twin", "쌍태 임신 (성호르몬 부하 ×1.55)", FALSE),
      sliderInput("del_ga", "분만 주수 (GA, weeks)", 34, 41, 39, step = 0.5),
      numericInput("bet_ga", "베타메타손 투여 주수 (0 = 없음)", 0, 0, 40, 0.5),
      hr(),
      h4("담즙산 축 치료"),
      sliderInput("udca_mg", "UDCA (mg/일, 2회 분할)", 0, 2000, 0, step = 250),
      sliderInput("udca_from", "UDCA 시작 주수", 26, 38, 30, step = 0.5),
      sliderInput("chol_g", "콜레스티라민 (g/일, 4회 분할)", 0, 24, 0, step = 4),
      sliderInput("chol_from", "콜레스티라민 시작 주수", 26, 38, 32, step = 0.5),
      sliderInput("ibat", "IBAT 억제제 (가설, 장내 µmol/일)", 0, 8, 0, step = 1),
      sliderInput("same_mg", "SAMe (mg/일 IV)", 0, 1600, 0, step = 200),
      checkboxInput("vitk", "비타민 K 10 mg/일 보충 (34주부터)", FALSE),
      hr(),
      h4("가려움 축 치료"),
      sliderInput("rif_mg", "리팜피신 (mg/일, 2회 분할)", 0, 600, 0, step = 150),
      sliderInput("rif_from", "리팜피신 시작 주수", 26, 38, 32, step = 0.5),
      sliderInput("ntx_mg", "날트렉손 (mg/일)", 0, 50, 0, step = 25),
      sliderInput("ah_mg", "항히스타민 (mg/일)", 0, 24, 0, step = 4),
      hr(),
      h4("비교군 (탭 9)"),
      selectInput("cmp_gen", "비교군 유전 소인", choices = names(GENOTYPES),
                  selected = names(GENOTYPES)[4]),
      sliderInput("cmp_udca", "비교군 UDCA (mg/일)", 0, 2000, 1000, step = 250),
      sliderInput("cmp_rif", "비교군 리팜피신 (mg/일)", 0, 600, 0, step = 150),
      sliderInput("cmp_del", "비교군 분만 주수", 34, 41, 39, step = 0.5)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---- 1 -------------------------------------------------------
        tabPanel(
          "1. 환자 프로파일",
          br(),
          fluidRow(column(6, plotOutput("p_steroid", height = 300)),
                   column(6, plotOutput("p_inhib", height = 300))),
          hr(),
          h4("현재 설정에서의 요약"),
          DTOutput("t_summary"),
          br(),
          div(style = "background:#f7f7fb;padding:10px;border-radius:6px;",
              strong("읽는 방법."),
              " 왼쪽 두 그래프에 이 모델의 시간 구조 전체가 있습니다.",
              " 태반 성호르몬은 만삭까지 지수적으로 오르고, 그 대사체(E2-17G 와",
              " 프로게스테론 설페이트)가 BSEP 를 억제합니다. 그래서 같은 유전형이",
              " 20주에는 조용하고 제3삼분기에 진단 문턱을 넘으며, 분만하면 며칠 안에",
              " 사라집니다 — 어느 쪽 경계에서도 파라미터를 바꾸지 않았습니다.")
        ),

        ## ---- 2 -------------------------------------------------------
        tabPanel(
          "2. 담즙산 검사값 vs 실제",
          br(),
          plotOutput("p_tba", height = 340),
          hr(),
          plotOutput("p_species", height = 300),
          br(),
          div(style = "background:#fff4e5;padding:10px;border-radius:6px;",
              strong("이 탭이 존재하는 이유."),
              " 임상에서 쓰는 총담즙산 검사는 3α-수산화스테로이드 탈수소효소",
              " 방식이므로 UDCA 와 하이오콜산까지 함께 셉니다. 위 그래프의 실선",
              "(검사값)과 점선(내인성 분율)이 벌어지는 폭이 곧 아무도 약에",
              " 귀속시키지 않는 약효입니다. 아래 그래프는 종별 조성으로, UDCA",
              " 투여 시 리토콜산(LCA, 세포독성 가중치 1.0)이 어떻게 변하는지를",
              " 보여줍니다 — 장내 세균의 7β-탈수산화 산물입니다.")
        ),

        ## ---- 3 -------------------------------------------------------
        tabPanel(
          "3. 간 기능 · 약물 PK",
          br(),
          fluidRow(column(6, plotOutput("p_alt", height = 300)),
                   column(6, plotOutput("p_pk", height = 300))),
          hr(),
          fluidRow(column(6, plotOutput("p_transport", height = 300)),
                   column(6, plotOutput("p_cyp", height = 300))),
          br(),
          div(style = "background:#f7f7fb;padding:10px;border-radius:6px;",
              strong("주의."),
              " 오른쪽 아래의 CYP3A4 유도는 리팜피신의 치료 기전(6α-수산화 →",
              " 하이오콜산, SULT2A1 유도 → 신배설)이면서 동시에 UGT1A1 을 같이",
              " 유도해 BSEP 억제제인 E2-17G 를 늘립니다. 모델은 이 상충을",
              " 감추지 않고 그대로 계산합니다.")
        ),

        ## ---- 4 -------------------------------------------------------
        tabPanel(
          "4. 태아 노출과 심근",
          br(),
          fluidRow(column(6, plotOutput("p_fetal", height = 300)),
                   column(6, plotOutput("p_wbar", height = 300))),
          hr(),
          fluidRow(column(6, plotOutput("p_heart", height = 300)),
                   column(6, plotOutput("p_arri", height = 300))),
          br(),
          div(style = "background:#ffecec;padding:10px;border-radius:6px;",
              strong("문턱이 실제로 있는 곳."),
              " 세 개의 층화된 임신에서 모체 총담즙산은 8.7배 차이나는데 부정맥",
              " 지수는 37배 차이납니다. 사산 위험을 모체 담즙산에 쓰면 지수 2.95,",
              " 태아 소수성 부하에 쓰면 2.55, 부정맥 지수에 쓰면 1.55 가 필요합니다.",
              " 즉 100 µmol/L 문턱의 급경사는 태반 수송이나 용량-반응이 아니라",
              " 코넥신-43 탈결합의 문턱 성질에 있습니다. 이 모델의 지도는 반대를",
              " 예상하고 그려졌고, 적합이 그것을 정정했습니다.")
        ),

        ## ---- 5 -------------------------------------------------------
        tabPanel(
          "5. 태반 수송 포화",
          br(),
          plotOutput("p_placenta", height = 340),
          hr(),
          fluidRow(column(6, plotOutput("p_ratio", height = 300)),
                   column(6, plotOutput("p_hyp", height = 300))),
          br(),
          div(style = "background:#f7f7fb;padding:10px;border-radius:6px;",
              strong("태반 수송이 설명하는 것과 설명하지 않는 것."),
              " 확산성 모체→태아 유입과 태아 자체 합성이 수송체 용량 V_P 를",
              " 소진하면 제대:모체 비가 내려가다가 다시 올라갑니다(왼쪽 아래).",
              " 그것이 태아 노출의 절대 수준을 정합니다. 다만 절제 실험",
              "(icp_calibration.py I 절)은 이 포화를 없애도 위험 곡선의 굽음이",
              " 사라지지 않는다는 것을 보입니다 — 굽음은 탭 4 의 심근에 있습니다.")
        ),

        ## ---- 6 -------------------------------------------------------
        tabPanel(
          "6. 가려움 축 (별개의 계)",
          br(),
          fluidRow(column(6, plotOutput("p_vas", height = 300)),
                   column(6, plotOutput("p_atx", height = 300))),
          hr(),
          plotOutput("p_axes", height = 320),
          br(),
          div(style = "background:#ecf0ff;padding:10px;border-radius:6px;",
              strong("왜 UDCA 는 가려움에 거의 듣지 않는가."),
              " 오토탁신은 주로 프로게스테론 설페이트와 에스트라디올로 구동되고",
              " 담즙산에는 약하게, 그리고 포화되게만 반응합니다. 이 구조는 네 가지",
              " 관찰을 동시에 만족시키기 위한 것입니다 — 혈중 오토탁신은 UDCA 로",
              " 내려가지 않고 리팜피신으로는 내려간다, 가려움 정도는 총담즙산과",
              " 잘 맞지 않는다, 상당수 환자에서 가려움이 생화학 이상보다 몇 주",
              " 앞선다, 무증상 고담즙산혈증이 존재한다. 담즙산으로 가려움을 구동하는",
              " 모델은 이 중 어느 것도 재현하지 못합니다.")
        ),

        ## ---- 7 -------------------------------------------------------
        tabPanel(
          "7. 임상 엔드포인트 · PITCHES 분해",
          br(),
          fluidRow(column(6, plotOutput("p_haz", height = 300)),
                   column(6, plotOutput("p_comp", height = 300))),
          hr(),
          h4("복합 엔드포인트의 구성 (분만 시점 기준)"),
          DTOutput("t_endpoints"),
          br(),
          div(style = "background:#fff4e5;padding:10px;border-radius:6px;",
              strong("PITCHES 가 음성이었던 산수."),
              " 복합 엔드포인트의 세 성분 중 둘(37주 미만 분만, 신생아집중치료실",
              " 입원)은 분만 주수와 분만 결정의 함수이고 담즙산의 함수가 아닙니다.",
              " 사산은 위약군 복합의 1% 미만을 차지하므로, 사산을 완전히 없애도",
              " 복합은 0.2%p 밖에 움직이지 않습니다. 604명으로는 자기 기전을",
              " 검출할 수 없었고, 그래서 이 결과는 개별환자자료 메타분석의 사산",
              " 신호와 모순되지 않습니다 — 두 연구는 서로 다른 것을 재고 있습니다.")
        ),

        ## ---- 8 -------------------------------------------------------
        tabPanel(
          "8. 분만 시기 최적화",
          br(),
          sliderInput("wmorb", "신생아 이환 : 사산 효용비 (WMORB)",
                      0.005, 0.15, 0.055, step = 0.005, width = "60%"),
          plotOutput("p_deliv", height = 380),
          hr(),
          DTOutput("t_deliv"),
          br(),
          div(style = "background:#ffecec;padding:10px;border-radius:6px;",
              strong("이 탭은 지침과 부분적으로만 일치합니다."),
              " ≥100 µmol/L 구간에서는 모델이 36–37주를 최적으로 계산하여 지침과",
              " 같은 방향을 가리키고, 치료하면 최적점이 39주로 밀립니다. 그러나",
              " 40–99 구간에서는 39–40주가 나오며 지침의 37–38주를 재현하지",
              " 못합니다. 이유는 산수입니다: 그 구간의 사산 초과 위험은 0.28% 대",
              " 0.13% 로 0.15%p 이고, 어떤 합리적 가중치로도 2–3주 조기분만의",
              " 신생아 비용을 넘지 못합니다. 위 슬라이더를 0.01 아래로 내려 보면",
              " 37주가 최적이 되지만 그러면 <40 구간도 함께 당겨져서 어느 지침도",
              " 권하지 않는 결론이 됩니다. 모델이 무엇을 놓쳤을 가능성이 큰 쪽은",
              " 담즙산의 개인내 변동입니다 — 70 µmol/L 로 측정된 사람이 방문 사이에",
              " 100 을 넘는 날을 보낼 수 있고, 이 모델은 매끄러운 궤적을 돌립니다.")
        ),

        ## ---- 9 -------------------------------------------------------
        tabPanel(
          "9. 시나리오 비교",
          br(),
          plotOutput("p_cmp", height = 480),
          hr(),
          DTOutput("t_cmp"),
          br(),
          div(style = "background:#f7f7fb;padding:10px;border-radius:6px;",
              strong("권하는 비교."),
              " 주 설정을 '리팜피신 300 mg 2회', 비교군을 'UDCA 1000 mg' 으로",
              " 두고 TBA 패널과 VAS 패널을 나란히 보십시오. 두 약의 순위가 두",
              " 패널에서 뒤집힙니다. 그것이 이 모델이 말하려는 전부입니다.")
        ),

        ## ---- 10 ------------------------------------------------------
        tabPanel(
          "10. 바이오마커 · 안전성",
          br(),
          fluidRow(column(6, plotOutput("p_vitk", height = 300)),
                   column(6, plotOutput("p_mec", height = 300))),
          hr(),
          fluidRow(column(6, plotOutput("p_uterus", height = 300)),
                   column(6, plotOutput("p_surf", height = 300))),
          br(),
          div(style = "background:#eefaf6;padding:10px;border-radius:6px;",
              strong("가격이 붙은 상충."),
              " 콜레스티라민은 담즙산을 내리면서 동시에 장내 담즙산 농도를 낮춰",
              " 미셀 형성을 줄이므로 비타민 K 상태를 악화시킵니다(왼쪽 위 INR).",
              " 리팜피신은 풀 자체를 줄여 같은 결과를 냅니다. 두 약 모두 산후 출혈",
              " 위험을 담즙산과 교환하며, 모델은 그 교환비를 계산합니다. UDCA 와",
              " IBAT 억제제는 이 대가를 치르지 않습니다.")
        )
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  gen_main <- reactive(GENOTYPES[[input$gen]])

  ev_main <- reactive({
    build_events(input$udca_mg, input$udca_from, input$rif_mg, input$rif_from,
                 input$chol_g, input$chol_from, input$same_mg, input$ibat,
                 input$ntx_mg, input$ah_mg,
                 if (input$bet_ga > 0) input$bet_ga else NA, input$del_ga)
  })

  sim <- reactive({
    extra <- if (input$vitk) list(VKDOSE = 1.6, VKSTART = ga2d(34)) else list()
    run_case(gen_main(), ev_main(), input$del_ga,
             twin = if (input$twin) 1.55 else 1, extra = extra)
  })

  sim_cmp <- reactive({
    e <- build_events(input$cmp_udca, 30, input$cmp_rif, 32, 0, 32, 0, 0, 0, 0,
                      NA, input$cmp_del)
    run_case(GENOTYPES[[input$cmp_gen]], e, input$cmp_del,
             twin = if (input$twin) 1.55 else 1)
  })

  at_delivery <- function(d, del_ga) {
    d %>% filter(GA <= del_ga + 1e-6) %>% slice_tail(n = 1)
  }

  vline_del <- function(del_ga)
    geom_vline(xintercept = del_ga, linetype = "dotted", colour = "grey40")

  ## ---- tab 1 --------------------------------------------------------
  output$p_steroid <- renderPlot({
    d <- sim()
    d %>% select(GA, `에스트라디올 (nmol/L)` = E2,
                 `프로게스테론 설페이트 (nmol/L)` = P4S,
                 `E2-17G (nmol/L ×10)` = E2G) %>%
      mutate(`E2-17G (nmol/L ×10)` = `E2-17G (nmol/L ×10)` * 10) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) +
      geom_line(linewidth = 0.9) + vline_del(input$del_ga) +
      labs(title = "태반 성호르몬과 그 BSEP 억제 대사체",
           subtitle = "만삭까지 지수적 상승, 분만 후 급속 소실 (점선 = 분만)",
           x = "재태 주수", y = "농도", colour = NULL) + THEME
  })

  output$p_inhib <- renderPlot({
    d <- sim()
    d %>% mutate(`총담즙산 검사값` = TBA,
                 `진단 문턱 10` = 10, `중증 40` = 40, `최중증 100` = 100) %>%
      select(GA, `총담즙산 검사값`, `진단 문턱 10`, `중증 40`, `최중증 100`) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name,
                 linetype = name == "총담즙산 검사값")) +
      geom_line(linewidth = 0.9) + vline_del(input$del_ga) +
      scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"),
                            guide = "none") +
      labs(title = "같은 유전형이 언제 진단 기준을 넘는가",
           subtitle = "발현 시기와 분만 후 소실 모두 파라미터 변경 없이 나온다",
           x = "재태 주수", y = "µmol/L", colour = NULL) + THEME
  })

  output$t_summary <- renderDT({
    d <- at_delivery(sim(), input$del_ga)
    post <- sim() %>% filter(GA > input$del_ga, TBA < 10) %>% slice_head(n = 1)
    onset <- sim() %>% filter(GA <= input$del_ga, TBA >= 10) %>% slice_head(n = 1)
    tibble(
      항목 = c("분만 시 총담즙산 검사값 (µmol/L)",
               "분만 시 내인성 담즙산 (µmol/L)",
               "혈중 UDCA (µmol/L)", "하이오콜산 (µmol/L)",
               "제대혈 총담즙산 (µmol/L)", "태아 소수성 부하 FCL",
               "태아 풀 소수성도 w̄", "부정맥 지수 ARRI",
               "ALT (U/L)", "가려움 VAS (0-10)", "오토탁신 (상대값)", "INR",
               "진단 문턱(10) 통과 주수", "분만 후 10 µmol/L 미만 도달"),
      값 = c(sprintf("%.1f", d$TBA), sprintf("%.1f", d$TBA_ENDO),
             sprintf("%.1f", d$UDCA_P), sprintf("%.1f", d$HCA_P),
             sprintf("%.1f", d$CORD), sprintf("%.2f", d$FCLo),
             sprintf("%.3f", d$WBAR), sprintf("%.4f", d$ARRIo),
             sprintf("%.0f", d$ALT), sprintf("%.1f", d$VAS),
             sprintf("%.2f", d$ATX), sprintf("%.2f", d$INR),
             if (nrow(onset)) sprintf("%.1f 주", onset$GA) else "넘지 않음",
             if (nrow(post)) sprintf("%.0f 일", (post$GA - input$del_ga) * 7)
             else "21일 내 미도달")
    )
  }, options = list(dom = "t", pageLength = 20), rownames = FALSE)

  ## ---- tab 2 --------------------------------------------------------
  output$p_tba <- renderPlot({
    d <- sim()
    d %>% select(GA, `검사값 (UDCA·HCA 포함)` = TBA,
                 `내인성 분율` = TBA_ENDO, `혈중 UDCA` = UDCA_P,
                 `하이오콜산` = HCA_P) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name, linetype = name)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = c(10, 40, 100), linetype = "dotted",
                 colour = "grey60") +
      vline_del(input$del_ga) +
      scale_linetype_manual(values = c("검사값 (UDCA·HCA 포함)" = "solid",
                                       "내인성 분율" = "longdash",
                                       "혈중 UDCA" = "dotted",
                                       "하이오콜산" = "dotdash")) +
      labs(title = "임상이 재는 값과 임상이 재고 있다고 믿는 값",
           subtitle = "두 선의 간격 = 아무도 약에 귀속시키지 않는 약효",
           x = "재태 주수", y = "µmol/L", colour = NULL, linetype = NULL) + THEME
  })

  output$p_species <- renderPlot({
    d <- sim()
    d %>% mutate(CA = SPCA / 18, CDCA = SPCD / 18, DCA = SPDC / 18,
                 LCA = SPLC / 18, UDCA = SPUD / 18) %>%
      select(GA, CA, CDCA, DCA, LCA, UDCA) %>%
      pivot_longer(-GA) %>%
      mutate(name = factor(name, levels = c("CA", "CDCA", "DCA", "LCA", "UDCA"))) %>%
      ggplot(aes(GA, value, fill = name)) +
      geom_area(position = "stack", alpha = 0.85) + vline_del(input$del_ga) +
      scale_fill_brewer(palette = "RdYlBu", direction = -1) +
      labs(title = "모체 혈중 담즙산 종별 조성",
           subtitle = "세포독성 가중치 CA 0.20 · CDCA 0.60 · DCA 0.72 · LCA 1.00 · UDCA 0.02",
           x = "재태 주수", y = "µmol/L", fill = NULL) + THEME
  })

  ## ---- tab 3 --------------------------------------------------------
  output$p_alt <- renderPlot({
    d <- sim()
    d %>% select(GA, `ALT (U/L)` = ALT,
                 `산화 스트레스 ×100` = ROS) %>%
      mutate(`산화 스트레스 ×100` = `산화 스트레스 ×100` * 100) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      vline_del(input$del_ga) +
      labs(title = "간세포 손상", subtitle = "담즙 담즙산:PC 비 상승이 ABCB4 변이의 기전",
           x = "재태 주수", y = NULL, colour = NULL) + THEME
  })

  output$p_pk <- renderPlot({
    d <- sim()
    d %>% mutate(`리팜피신 (mg/L)` = RIFC / 45,
                 `UDCA 혈중 (µmol/L)` = UDCA_P,
                 `SAMe (mg/L)` = SAMC / 22) %>%
      select(GA, `리팜피신 (mg/L)`, `UDCA 혈중 (µmol/L)`, `SAMe (mg/L)`) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.8) +
      vline_del(input$del_ga) +
      labs(title = "약물 노출", subtitle = "투여 중단은 분만 시점",
           x = "재태 주수", y = NULL, colour = NULL) + THEME
  })

  output$p_transport <- renderPlot({
    d <- sim()
    d %>% select(GA, `기저측 배출 MRP3/4/OSTαβ` = MRP4,
                 `SULT2A1 유도` = SULT, `CYP7A1` = CYP7A1,
                 `FGF19` = FGF19) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      vline_del(input$del_ga) +
      labs(title = "적응 반응", subtitle = "탈출 밸브는 기저에서 잠잠하고 담즙정체에서 크게 유도된다",
           x = "재태 주수", y = "상대값", colour = NULL) + THEME
  })

  output$p_cyp <- renderPlot({
    d <- sim()
    d %>% select(GA, `CYP3A4 유도` = CYP3A, `E2-17G (nmol/L)` = E2G) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      vline_del(input$del_ga) +
      labs(title = "리팜피신의 상충하는 두 효과",
           subtitle = "CYP3A4 유도(치료)와 UGT1A1 공동유도 → E2-17G 증가(역효과)",
           x = "재태 주수", y = NULL, colour = NULL) + THEME
  })

  ## ---- tab 4 --------------------------------------------------------
  output$p_fetal <- renderPlot({
    d <- sim()
    d %>% select(GA, `모체 검사값` = TBA, `제대혈 총담즙산` = CORD,
                 `태아 소수성 부하 FCL` = FCLo) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 1) +
      vline_del(input$del_ga) +
      labs(title = "모체 · 태아 노출", subtitle = "태아 부하는 총량이 아니라 소수성 가중합",
           x = "재태 주수", y = "µmol/L", colour = NULL) + THEME
  })

  output$p_wbar <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    ggplot(d, aes(GA, WBAR)) + geom_line(linewidth = 1, colour = "#6a1b9a") +
      labs(title = "태아 풀의 평균 소수성도 w̄",
           subtitle = "소수성 종의 확산 투과도가 더 높아 풀이 커질 뿐 아니라 나빠진다",
           x = "재태 주수", y = "FCL / 제대혈 총담즙산") + THEME
  })

  output$p_heart <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    d %>% select(GA, `코넥신-43 결합 GJ` = GJ, `칼슘 과부하 지수` = FCA) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 1) +
      labs(title = "태아 심근 두 상태", subtitle = "GJ 는 문턱형(Hill 1.6), Ca 는 완만한 포화형",
           x = "재태 주수", y = "상대값 (정상 = 1)", colour = NULL) + THEME
  })

  output$p_arri <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    ggplot(d, aes(GA, ARRIo)) +
      geom_line(linewidth = 1.1, colour = "#b71c1c") +
      labs(title = "부정맥 지수 ARRI = (1 − GJ) × Ca",
           subtitle = "사산 위험 h(t) = HSB0 + HSBSC · ARRI^1.55 · (1 + 0.55·저산소)",
           x = "재태 주수", y = "ARRI") + THEME
  })

  ## ---- tab 5 --------------------------------------------------------
  output$p_placenta <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    d %>% mutate(fCA = FPCA, fCD = FPCD, fDC = FPDC, fLC = FPLC) %>%
      select(GA, `태아 CA` = fCA, `태아 CDCA` = fCD, `태아 DCA` = fDC,
             `태아 LCA` = fLC) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(title = "태아 구획 담즙산 종별 (µmol)",
           subtitle = "확산은 양방향이므로 태아 농도는 모체 농도 약간 위에서 포화한다",
           x = "재태 주수", y = "µmol", colour = NULL) + THEME
  })

  output$p_ratio <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga, TBA > 0.1)
    ggplot(d, aes(TBA, CORD / TBA)) +
      geom_path(linewidth = 1, colour = "#d84315") +
      geom_vline(xintercept = 100, linetype = "dotted") +
      labs(title = "제대혈 : 모체 비", subtitle = "내려가다 다시 올라가는 지점이 수송체 포화",
           x = "모체 총담즙산 (µmol/L)", y = "제대 / 모체") + THEME
  })

  output$p_hyp <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    d %>% select(GA, `융모판 혈관수축` = VASO, `영양막 손상` = TROPH,
                 `태아 저산소 지수` = HYP) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(title = "태반 혈관과 저산소", x = "재태 주수", y = "상대값",
           colour = NULL) + THEME
  })

  ## ---- tab 6 --------------------------------------------------------
  output$p_vas <- renderPlot({
    d <- sim()
    ggplot(d, aes(GA, VAS)) + geom_line(linewidth = 1.1, colour = "#283593") +
      vline_del(input$del_ga) + ylim(0, 10) +
      labs(title = "가려움 VAS (0-10 cm)",
           subtitle = "정상 임신 ~1.3 · 경증 ICP ~4 · 중증 ~5.6 · 최중증 ~7",
           x = "재태 주수", y = "VAS") + THEME
  })

  output$p_atx <- renderPlot({
    d <- sim()
    d %>% select(GA, `오토탁신 ATX` = ATX, `LPA` = LPA,
                 `중추 가려움 상태` = ITCHC) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      vline_del(input$del_ga) +
      labs(title = "가려움 축의 내부", subtitle = "리팜피신만 ATX 를 내린다 (PXR 의존)",
           x = "재태 주수", y = "상대값", colour = NULL) + THEME
  })

  output$p_axes <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    rng <- function(x) if (diff(range(x)) < 1e-9) rep(0.5, length(x)) else
      (x - min(x)) / diff(range(x))
    d %>% mutate(`담즙산 축 (정규화)` = rng(TBA_ENDO),
                 `태아 부하 FCL (정규화)` = rng(FCLo),
                 `가려움 축 VAS (정규화)` = rng(VAS)) %>%
      select(GA, `담즙산 축 (정규화)`, `태아 부하 FCL (정규화)`,
             `가려움 축 VAS (정규화)`) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 1.1) +
      labs(title = "두 축을 같은 축척에 겹쳐 보기",
           subtitle = "처방을 바꿀 때 두 선이 함께 움직이는지 따로 움직이는지를 보라",
           x = "재태 주수", y = "각자의 범위로 정규화", colour = NULL) + THEME
  })

  ## ---- tab 7 --------------------------------------------------------
  output$p_haz <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    d %>% select(GA, `누적 사산 위험 %` = SBRISK,
                 `누적 자연조산 위험 %` = PTBRISK,
                 `누적 태변 위험 %` = MECRISK) %>%
      mutate(across(-GA, ~ .x * 100)) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 1) +
      labs(title = "누적 위험", x = "재태 주수", y = "%", colour = NULL) + THEME
  })

  output$p_comp <- renderPlot({
    d <- at_delivery(sim(), input$del_ga)
    ptb_any <- if (input$del_ga < 37) 1 else d$PTBRISK
    tibble(성분 = c("사산", "37주 미만 분만", "신생아집중치료실 ≥4h"),
           확률 = c(d$SBRISK, ptb_any, d$NICURISK) * 100) %>%
      ggplot(aes(reorder(성분, 확률), 확률, fill = 성분)) +
      geom_col(width = 0.6, show.legend = FALSE) +
      geom_text(aes(label = sprintf("%.2f%%", 확률)), hjust = -0.1, size = 4) +
      coord_flip(clip = "off") +
      scale_fill_manual(values = c("사산" = "#b71c1c",
                                   "37주 미만 분만" = "#f9a825",
                                   "신생아집중치료실 ≥4h" = "#1565c0")) +
      labs(title = "복합 엔드포인트의 구성",
           subtitle = "사산은 복합의 1% 미만 — 없애도 복합은 거의 안 움직인다",
           x = NULL, y = "%") + THEME +
      theme(plot.margin = margin(5, 40, 5, 5))
  })

  output$t_endpoints <- renderDT({
    d <- at_delivery(sim(), input$del_ga)
    ptb_any <- if (input$del_ga < 37) 1 else d$PTBRISK
    comp <- 1 - (1 - d$SBRISK) * (1 - ptb_any) * (1 - d$NICURISK)
    tibble(
      엔드포인트 = c("사산 (누적)", "자연 조산 (누적)", "의인성 조산",
                     "태변 착색 양수", "신생아 호흡곤란", "신생아집중치료실 ≥4h",
                     "복합 엔드포인트 (PITCHES 정의)",
                     "복합 중 사산이 차지하는 비율"),
      `확률 %` = sprintf("%.3f", 100 * c(d$SBRISK, d$PTBRISK,
                                          if (input$del_ga < 37) 1 else 0,
                                          d$MECRISK, d$RDSRISK, d$NICURISK,
                                          comp, d$SBRISK / comp))
    )
  }, options = list(dom = "t", pageLength = 10), rownames = FALSE)

  ## ---- tab 8 --------------------------------------------------------
  deliv_curve <- reactive({
    gen <- gen_main()
    gas <- seq(35, 40, by = 0.5)
    res <- lapply(gas, function(g) {
      e <- build_events(input$udca_mg, input$udca_from, input$rif_mg,
                        input$rif_from, input$chol_g, input$chol_from,
                        input$same_mg, input$ibat, input$ntx_mg, input$ah_mg,
                        if (input$bet_ga > 0) input$bet_ga else NA, g)
      d <- run_case(gen, e, g, twin = if (input$twin) 1.55 else 1)
      a <- at_delivery(d, g)
      tibble(GA = g, 사산 = a$SBRISK,
             이환 = a$NICURISK + a$RDSRISK)
    })
    bind_rows(res) %>% mutate(총손실 = 사산 + input$wmorb * 이환)
  })

  output$p_deliv <- renderPlot({
    d <- deliv_curve()
    opt <- d$GA[which.min(d$총손실)]
    d %>% mutate(`사산 위험` = 사산 * 1000,
                 `신생아 이환 비용 (가중)` = input$wmorb * 이환 * 1000,
                 `합계` = 총손실 * 1000) %>%
      select(GA, `사산 위험`, `신생아 이환 비용 (가중)`, `합계`) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name, linewidth = name == "합계")) +
      geom_line() +
      geom_vline(xintercept = opt, linetype = "dashed", colour = "#b71c1c") +
      annotate("text", x = opt, y = Inf, vjust = 1.5, hjust = -0.1,
               label = sprintf("최적 %.1f 주", opt), colour = "#b71c1c") +
      scale_linewidth_manual(values = c("TRUE" = 1.4, "FALSE" = 0.9),
                             guide = "none") +
      labs(title = "분만 시기는 문턱이 아니라 두 곡선의 교차다",
           subtitle = "왼쪽으로 갈수록 사산은 줄고 신생아 비용은 급격히 늘어난다",
           x = "분만 재태 주수", y = "milli-expected-loss", colour = NULL) + THEME
  })

  output$t_deliv <- renderDT({
    deliv_curve() %>%
      transmute(`분만 주수` = GA,
                `사산 %` = sprintf("%.3f", 100 * 사산),
                `이환 %` = sprintf("%.1f", 100 * 이환),
                `총 손실 (×1000)` = sprintf("%.2f", 1000 * 총손실))
  }, options = list(dom = "t", pageLength = 12), rownames = FALSE)

  ## ---- tab 9 --------------------------------------------------------
  output$p_cmp <- renderPlot({
    a <- sim() %>% mutate(arm = "주 설정")
    b <- sim_cmp() %>% mutate(arm = "비교군")
    bind_rows(a, b) %>%
      select(GA, arm, `총담즙산 검사값` = TBA, `내인성 담즙산` = TBA_ENDO,
             `태아 소수성 부하 FCL` = FCLo, `가려움 VAS` = VAS,
             `ALT` = ALT, `INR` = INR) %>%
      pivot_longer(c(-GA, -arm)) %>%
      mutate(name = factor(name, levels = c("총담즙산 검사값", "내인성 담즙산",
                                            "태아 소수성 부하 FCL",
                                            "가려움 VAS", "ALT", "INR"))) %>%
      ggplot(aes(GA, value, colour = arm)) + geom_line(linewidth = 1) +
      facet_wrap(~ name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = c("주 설정" = "#1565c0",
                                     "비교군" = "#c62828")) +
      labs(title = "두 시나리오 비교",
           subtitle = "담즙산 패널과 가려움 패널에서 순위가 뒤집히는지 확인하라",
           x = "재태 주수", y = NULL, colour = NULL) + THEME
  })

  output$t_cmp <- renderDT({
    a <- at_delivery(sim(), input$del_ga)
    b <- at_delivery(sim_cmp(), input$cmp_del)
    f <- function(x) sprintf("%.2f", x)
    tibble(
      항목 = c("분만 주수", "총담즙산 검사값", "내인성 담즙산",
               "제대혈 총담즙산", "태아 소수성 부하 FCL", "부정맥 지수 ARRI",
               "가려움 VAS", "오토탁신", "ALT", "INR",
               "누적 사산 %", "신생아집중치료실 %"),
      `주 설정` = c(f(input$del_ga), f(a$TBA), f(a$TBA_ENDO), f(a$CORD),
                    f(a$FCLo), sprintf("%.4f", a$ARRIo), f(a$VAS), f(a$ATX),
                    sprintf("%.0f", a$ALT), f(a$INR),
                    sprintf("%.3f", 100 * a$SBRISK),
                    sprintf("%.1f", 100 * a$NICURISK)),
      `비교군` = c(f(input$cmp_del), f(b$TBA), f(b$TBA_ENDO), f(b$CORD),
                   f(b$FCLo), sprintf("%.4f", b$ARRIo), f(b$VAS), f(b$ATX),
                   sprintf("%.0f", b$ALT), f(b$INR),
                   sprintf("%.3f", 100 * b$SBRISK),
                   sprintf("%.1f", 100 * b$NICURISK))
    )
  }, options = list(dom = "t", pageLength = 15), rownames = FALSE)

  ## ---- tab 10 -------------------------------------------------------
  output$p_vitk <- renderPlot({
    d <- sim()
    d %>% select(GA, `비타민 K 상태` = VITK, `INR` = INR,
                 `미셀 형성능 대리 (PIVKA)` = PIVKA) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      vline_del(input$del_ga) +
      labs(title = "지용성 비타민과 응고",
           subtitle = "콜레스티라민과 리팜피신은 담즙산을 비타민 K 와 교환한다",
           x = "재태 주수", y = NULL, colour = NULL) + THEME
  })

  output$p_mec <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    d %>% select(GA, `태변 담즙산 축적 (µmol/100)` = MEC,
                 `누적 태변 위험 %` = MECRISK) %>%
      mutate(`태변 담즙산 축적 (µmol/100)` = `태변 담즙산 축적 (µmol/100)` / 100,
             `누적 태변 위험 %` = `누적 태변 위험 %` * 100) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(title = "태변", subtitle = "만삭 근처에서만 위험이 누적된다 (장운동 성숙)",
           x = "재태 주수", y = NULL, colour = NULL) + THEME
  })

  output$p_uterus <- renderPlot({
    d <- sim() %>% filter(GA <= input$del_ga)
    d %>% select(GA, `옥시토신 수용체 밀도` = OTR, `자궁 활동 지수` = UTA) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(title = "자궁 근육", subtitle = "담즙산이 OTR 을 올려 자연조산 위험을 만든다",
           x = "재태 주수", y = "상대값", colour = NULL) + THEME
  })

  output$p_surf <- renderPlot({
    gas <- seq(34, 41, by = 0.25)
    d <- at_delivery(sim(), input$del_ga)
    tibble(GA = gas,
           `호흡곤란 (스테로이드 없음)` =
             (0.0028 + 1 / (1 + exp(1.15 * (gas - 33.6)))) * 100,
           `호흡곤란 (현재 SURF)` =
             (0.0028 + 1 / (1 + exp(1.15 * (gas - 33.6)))) /
             (1 + 1.55 * max(0, d$SURF - 1)) * 100) %>%
      pivot_longer(-GA) %>%
      ggplot(aes(GA, value, colour = name)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = input$del_ga, linetype = "dotted") +
      labs(title = "의인성 조산의 비용",
           subtitle = "이 곡선의 급경사가 분만 시기 최적화의 절반이다",
           x = "분만 재태 주수", y = "%", colour = NULL) + THEME
  })
}

shinyApp(ui, server)
