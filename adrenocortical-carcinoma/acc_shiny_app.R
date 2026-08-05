## ===========================================================================
##  acc_shiny_app.R
##  Adrenocortical Carcinoma (ACC) QSP dashboard
##  부신피질암 QSP 대시보드
##
##      shiny::runApp("acc_shiny_app.R")
##
##  The app is organised around the model thesis: plasma mitotane is the
##  OVERFLOW of a slowly filling lipid depot, and that single state variable
##  is read three ways that pull in opposite directions. Every tab is one of
##  those reads, or the clock that governs them.
##
##  Tabs
##    1  환자 프로파일    Patient profile & the depot (Vss, tau, Css)
##    2  미토테인 PK      Mitotane PK: plasma vs depot, TDM window
##    3  스테로이드       Steroidogenesis, adrenolysis, free vs total cortisol
##    4  유도 & DDI       CYP3A4 induction and its four victims
##    5  종양 · 엔드포인트 Tumour, RECIST, response
##    6  시나리오 비교     Scenario comparison (any subset of the 18)
##    7  체성분 실험      Body-composition experiment (result A)
##    8  가상 집단        Virtual population (result B)
##    9  면역 · IGF       Immune / checkpoint and the IGF escape route
##   10  안전성           Safety: CNS threshold, liver, marrow, heart, kidney
##   11  문헌 대조        Calibration targets vs model output
## ===========================================================================

## mrgsolve is loaded first, then shiny: mrgsolve exports req() and filter(),
## and the Shiny helpers must win in this file.
library(mrgsolve)
suppressMessages(library(dplyr))
library(shiny)
library(ggplot2)
library(tidyr)
library(gridExtra)

source("acc_mrgsolve_model.R")

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        strip.text = element_text(face = "bold"))

## Therapeutic-window band. annotate() (not geom_rect + aes) so the single
## rectangle is not recycled once per data row.
WINDOW <- list(
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 14, ymax = 20,
           fill = "#c8e6c9", alpha = 0.35),
  geom_hline(yintercept = 14, linetype = "dashed", colour = "#2e7d32"),
  geom_hline(yintercept = 20, linetype = "dashed", colour = "#c62828")
)

## ---------------------------------------------------------------------------
ui <- navbarPage(
  title = "ACC QSP — 부신피질암",
  header = tags$style(HTML("
    .well { background: #fbfbf8; }
    .note { background:#fffde7; border-left:4px solid #c0a020;
            padding:8px 12px; margin:8px 0; font-size:12px; }
    .warn { background:#fdecec; border-left:4px solid #c62828;
            padding:8px 12px; margin:8px 0; font-size:12px; }
  ")),

  ## ---------------- 1. patient profile ------------------------------------
  tabPanel("1 환자 프로파일",
    sidebarLayout(
      sidebarPanel(width = 4,
        h4("체성분 (Body composition)"),
        sliderInput("fat",  "지방량 Adipose mass (kg)", 5, 60, 22, 1),
        sliderInput("bsa",  "BSA (m²)", 1.4, 2.4, 1.85, 0.05),
        checkboxInput("secretor", "코르티솔 분비 종양 (cortisol-secreting)", TRUE),
        hr(),
        h4("미토테인 투여 (Mitotane regimen)"),
        selectInput("regimen", "시작 전략 Starting strategy",
                    c("High-dose start (6→4→3 g/d)" = "high",
                      "Low-dose ramp (1→2→3 g/d)"   = "low",
                      "Flat 3 g/d"                  = "flat3",
                      "Flat 4 g/d"                  = "flat4",
                      "None"                        = "none")),
        checkboxInput("meal", "고지방 식사와 함께 복용 (fatty meal, F 0.55)", TRUE),
        sliderInput("clmult", "청소율 배수 Clearance multiplier", 0.4, 2.5, 1.0, 0.05),
        hr(),
        sliderInput("tend", "시뮬레이션 기간 (days)", 120, 900, 480, 30)
      ),
      mainPanel(width = 8,
        div(class = "note", HTML("<b>이 탭의 요점.</b> Css = D·F/CL 는 <i>어디에</i>
          도달하는지를 정하고, τ = Vss/CL 는 <i>언제</i> 도달하는지를 정합니다.
          지방량은 Css 를 전혀 바꾸지 않고 τ 만 바꿉니다.<br/>
          <b>The point.</b> Css sets WHERE the drug ends up; tau sets WHEN it
          gets there. Adipose mass moves only tau.")),
        fluidRow(
          column(6, h5("저장고 구성 (Depot composition)"), plotOutput("p_depot", height = 240)),
          column(6, h5("두 인자 (The two factors)"), tableOutput("t_depot"))
        ),
        h5("혈장 미토테인 (Plasma mitotane)"), plotOutput("p_cp1", height = 300)
      )
    )),

  ## ---------------- 2. mitotane PK ----------------------------------------
  tabPanel("2 미토테인 PK",
    fluidRow(
      column(6, h5("혈장 vs 지방 저장고"), plotOutput("p_pkboth", height = 300)),
      column(6, h5("체내 총량 중 저장고 비율"), plotOutput("p_depfrac", height = 300))),
    fluidRow(
      column(6, h5("치료 범위 내 누적 시간"), plotOutput("p_tiw", height = 280)),
      column(6, h5("자기유도에 의한 청소율 변화"), plotOutput("p_clmit", height = 280))),
    div(class = "warn", HTML("<b>중단 후에도 계속 투여된다.</b> 저장고가 비워지는 데
      수개월이 걸리므로, 처방을 멈춘 뒤에도 CYP3A4 유도와 부신 기능저하가
      남아 있습니다. 미토테인에는 washout 이 없습니다."))),

  ## ---------------- 3. steroidogenesis -----------------------------------
  tabPanel("3 스테로이드",
    sidebarLayout(
      sidebarPanel(width = 3,
        sliderInput("hcdose", "하이드로코르티손 (mg/day)", 0, 80, 20, 5),
        checkboxInput("cbgon", "CBG 상승 반영 (CBG rise on)", TRUE),
        checkboxInput("indhc", "코르티솔 청소율 유도 반영 (induction on)", TRUE),
        div(class = "note", HTML("두 스위치를 따로 끄면 <b>환자의 오차</b>(유도)와
          <b>검사의 오차</b>(CBG)를 분리할 수 있습니다."))),
      mainPanel(width = 9,
        fluidRow(
          column(6, h5("총 vs 유리 코르티솔"), plotOutput("p_cort", height = 280)),
          column(6, h5("총/유리 비 — 검사 함정"), plotOutput("p_ratio", height = 280))),
        fluidRow(
          column(6, h5("정상 부신 질량 · 종양 분비능"), plotOutput("p_adrn", height = 260)),
          column(6, h5("ACTH · DHEAS · 17-OHP · 알도스테론"), plotOutput("p_ster", height = 260))),
        tableOutput("t_cort")))),

  ## ---------------- 4. induction / DDI -----------------------------------
  tabPanel("4 유도 & DDI",
    fluidRow(
      column(6, h5("CYP3A4 활성 (배수)"), plotOutput("p_enz", height = 300)),
      column(6, h5("에토포시드 청소율 배수"), plotOutput("p_indfold", height = 300))),
    fluidRow(
      column(12, h5("반사실 실험: 피해약물 유도 ON vs OFF"),
             plotOutput("p_counter", height = 320))),
    div(class = "note", HTML("<b>결과 D.</b> 유도는 용량이 아니라 <i>저장고</i>를
      따라갑니다. 그래서 에토포시드 노출이 가장 낮은 시점이 미토테인이 드디어
      치료 농도에 도달한 시점입니다 — 병용요법의 두 축이 하나의 상태변수를 통해
      서로 상쇄됩니다."))),

  ## ---------------- 5. tumour / endpoints --------------------------------
  tabPanel("5 종양 · 엔드포인트",
    sidebarLayout(
      sidebarPanel(width = 3,
        checkboxGroupInput("addrug", "병용 요법 (add-on therapy)",
          c("EDP 항암화학 (q28d)"    = "edp",
            "스트렙토조토신"          = "sz",
            "펨브롤리주맙"            = "pem",
            "린시티닙"                = "lin"),
          selected = "edp"),
        sliderInput("cycles", "화학요법 주기 수", 2, 12, 8, 1),
        sliderInput("tum0", "초기 종양 부피 (mL)", 20, 800, 120, 20)),
      mainPanel(width = 9,
        fluidRow(
          column(6, h5("종양 부피 (민감 · 내성 클론)"), plotOutput("p_tum", height = 290)),
          column(6, h5("RECIST 유사 직경 합"), plotOutput("p_sld", height = 290))),
        h5("엔드포인트 요약"), tableOutput("t_end"),
        div(class = "note", HTML("RECIST 진행은 <b>최저점(nadir)</b> 기준이므로,
          더 깊게 반응한 환자는 +20%/+5 mm 기준선이 더 낮아져 오히려 PFS 가
          짧게 보일 수 있습니다. 기저 직경 회복 시점을 함께 보십시오."))))),

  ## ---------------- 6. scenario comparison -------------------------------
  tabPanel("6 시나리오 비교",
    sidebarLayout(
      sidebarPanel(width = 3,
        checkboxGroupInput("scensel", "시나리오 선택 (18개 중)",
          setNames(names(scenarios), vapply(scenarios, `[[`, "", "label")),
          selected = c("s01", "s03", "s04")),
        selectInput("scenvar", "표시 변수",
          c("Plasma mitotane" = "CMITo", "Tumour volume" = "TUMTOTo",
            "Free cortisol" = "CFREEo", "Total cortisol" = "CORT",
            "CYP3A4 activity" = "ENZ", "Adrenal mass" = "ADRN",
            "Neutrophils" = "CIRCN", "Effector T cells" = "TEFF",
            "Adipose mass" = "FATKG", "RECIST SLD" = "SLDo"))),
      mainPanel(width = 9, plotOutput("p_scen", height = 520),
                tableOutput("t_scen")))),

  ## ---------------- 7. body composition experiment -----------------------
  tabPanel("7 체성분 실험",
    div(class = "note", HTML("<b>결과 A.</b> 동일한 6/4/3 g 요법을 지방량만 바꿔
      투여합니다. 도착점(Css)은 같고 도달 시간(τ)만 달라집니다 — 그리고 20 mg/L
      를 넘겨 신경독성 위험에 놓이는 쪽은 <b>마른</b> 환자입니다.")),
    fluidRow(column(8, plotOutput("p_fatexp", height = 380)),
             column(4, tableOutput("t_fatexp"))),
    h5("도달 시간 vs 지방량 (연속)"), plotOutput("p_fatcurve", height = 300)),

  ## ---------------- 8. virtual population --------------------------------
  tabPanel("8 가상 집단",
    sidebarLayout(
      sidebarPanel(width = 3,
        sliderInput("npop", "가상 환자 수 (per arm)", 20, 300, 80, 20),
        actionButton("runpop", "실행 (run)", class = "btn-primary"),
        div(class = "note", HTML("보고된 겉보기 분포용적 변동(CV ~81.5%)을
          지방량과 청소율에 배분합니다. 표본이 클수록 오래 걸립니다."))),
      mainPanel(width = 9,
        h5("90일 시점 혈장 농도 분포"), plotOutput("p_pop", height = 300),
        h5("분산 분해 — 무엇이 결과를 지배하는가"), tableOutput("t_pop"),
        htmlOutput("popnote")))),

  ## ---------------- 9. immune / IGF --------------------------------------
  tabPanel("9 면역 · IGF",
    fluidRow(
      column(6, h5("효과 T세포 · Treg (코르티솔 의존)"), plotOutput("p_imm", height = 290)),
      column(6, h5("PD-1 점유율 vs 반응"), plotOutput("p_pd1", height = 290))),
    fluidRow(
      column(6, h5("IGF/AKT 신호 · 인슐린 구제"), plotOutput("p_igf", height = 290)),
      column(6, h5("혈당 · 인슐린 (린시티닙)"), plotOutput("p_glu", height = 290))),
    div(class = "warn", HTML("<b>결과 F.</b> 분비형 ACC 는 면역관문 억제제를 만나기
      전에 이미 <i>스스로 스테로이드 전처치</i>를 해 둔 상태입니다. 표적 점유율은
      99% 이상이지만 풀어줄 효과세포가 남아 있지 않습니다."))),

  ## ---------------- 10. safety -------------------------------------------
  tabPanel("10 안전성",
    fluidRow(
      column(6, h5("CNS 손상 점수 · 20 mg/L 역치"), plotOutput("p_ntox", height = 280)),
      column(6, h5("ALT · 유리 T4"), plotOutput("p_liver", height = 280))),
    fluidRow(
      column(6, h5("호중구 (Friberg)"), plotOutput("p_anc", height = 280)),
      column(6, h5("LVEF · eGFR"), plotOutput("p_organ", height = 280))),
    tableOutput("t_safety")),

  ## ---------------- 11. calibration --------------------------------------
  tabPanel("11 문헌 대조",
    h4("보정 목표 대비 모델 출력 (Calibration targets vs model)"),
    tableOutput("t_calib"),
    div(class = "note", HTML("모든 인용은 <code>acc_references.md</code> 참조.
      에토포시드 유도 계수는 미토테인–에토포시드 상호작용 연구가 아니라
      CYP3A4 유도 크기와 에토포시드의 CYP3A4 대사 분율에서 <b>추정</b>한
      값입니다 (모델 파일에 ASSUMED 로 표시)."))),

  ## ---------------- about -------------------------------------------------
  tabPanel("ℹ",
    div(style = "max-width:820px",
      h4("모델 개요"),
      HTML("<p>52개 ODE 구획, 18개 치료 시나리오, 가상 집단 레이어.
        하나의 상태변수(지방 내 미토테인)를 세 가지로 읽습니다:
        <b>① 부신파괴/효소차단(이익)</b>, <b>② CYP3A4 유도(해악)</b>,
        <b>③ 20 mg/L 초과 신경독성(해악)</b>.</p>
        <p><b>교육 및 연구 목적의 모델입니다.</b> 임상 의사결정, 처방,
        규제 제출에 사용하지 마십시오.</p>"))))

## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  sched_of <- function(reg) switch(reg,
    high  = sched_highdose,
    low   = sched_lowdose,
    flat3 = data.frame(g = 3, days = 900),
    flat4 = data.frame(g = 4, days = 900),
    none  = data.frame(g = 0, days = 900))

  mit_ev <- reactive(ev_mit_sched(sched_of(input$regimen)))

  base_par <- reactive(list(
    FATKG0 = input$fat, BSA = input$bsa,
    SECRETOR = as.numeric(input$secretor),
    MEALFAT = as.numeric(input$meal),
    CLMIT0 = 48 * input$clmult))

  ## core simulation used by tabs 1, 2, 4, 10
  sim <- reactive({
    do.call(run, c(list(events = mit_ev(), end = input$tend, delta = 1),
                   base_par()))
  })

  ## steroid simulation (tab 3) -- fine grid, tid hydrocortisone
  sim_ster <- reactive({
    ev_ <- comb(mit_ev(), ev_hc(input$hcdose, input$tend))
    p <- base_par()
    if (!input$cbgon) p$EMAXCBG <- 0
    if (!input$indhc) p$FMHC <- 0
    do.call(run, c(list(events = ev_, end = min(input$tend, 400), delta = 0.05), p))
  })

  ## tumour simulation (tab 5)
  sim_tum <- reactive({
    parts <- list(mit_ev())
    if ("edp" %in% input$addrug) parts <- c(parts, list(ev_edp(input$bsa, input$cycles)))
    if ("sz"  %in% input$addrug) parts <- c(parts, list(ev_sz(input$cycles)))
    if ("pem" %in% input$addrug) parts <- c(parts, list(ev_pem(16)))
    if ("lin" %in% input$addrug) parts <- c(parts, list(ev_lin(min(input$tend, 300))))
    parts <- c(parts, list(ev_hc(50, input$tend)))
    p <- c(base_par(), list(TUM0 = input$tum0))
    do.call(run, c(list(events = do.call(comb, parts),
                        end = input$tend, delta = 0.25), p))
  })

  ## ---- tab 1 -------------------------------------------------------------
  output$p_depot <- renderPlot({
    d <- data.frame(
      part = factor(c("central (blood + lean)", "non-adipose tissue", "ADIPOSE DEPOT"),
                    levels = c("central (blood + lean)", "non-adipose tissue", "ADIPOSE DEPOT")),
      L = c(400, 1100, 203 * input$fat))
    ggplot(d, aes("", L, fill = part)) +
      geom_col(width = 0.6) +
      scale_fill_manual(values = c("#cfe0f5", "#9dc3e6", "#1f4e79")) +
      coord_flip() + labs(x = NULL, y = "Volume (L)", fill = NULL) +
      ggtitle(sprintf("Vss = %.0f L", 1500 + 203 * input$fat)) + THEME
  })

  output$t_depot <- renderTable({
    a <- analytic_depot(3, input$fat, F_ = if (input$meal) 0.55 else 0.35,
                        CL = 48 * input$clmult * 2.05)
    data.frame(
      Quantity = c("Vss (L)", "Css = D·F/CL (mg/L)", "tau = Vss/CL (days)",
                   "terminal t½ (days)", "time to 14 mg/L on 3 g/d (days)"),
      Value = c(sprintf("%.0f", a$Vss), sprintf("%.1f", a$Css),
                sprintf("%.1f", a$tau_days), sprintf("%.1f", a$t_half_days),
                ifelse(is.na(a$time_to_14), "never reaches 14",
                       sprintf("%.0f", a$time_to_14))))
  }, striped = TRUE)

  output$p_cp1 <- renderPlot({
    s <- sim()
    ggplot(s, aes(time, CMITo)) + WINDOW +
      geom_line(linewidth = 1.1, colour = "#1f4e79") +
      labs(x = "Day", y = "Plasma mitotane (mg/L)") +
      ggtitle("14–20 mg/L therapeutic window shaded") + THEME
  })

  ## ---- tab 2 -------------------------------------------------------------
  output$p_pkboth <- renderPlot({
    s <- sim()
    d <- s %>% transmute(time, plasma = CMITo,
                         depot = MITA / (1100 + 203 * input$fat)) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = c(plasma = "#1f4e79", depot = "#8a6d1f")) +
      labs(x = "Day", y = "mg/L", colour = NULL) + THEME
  })
  output$p_depfrac <- renderPlot({
    ggplot(sim(), aes(time, 100 * DEPFRAC)) +
      geom_line(linewidth = 1, colour = "#8a6d1f") + ylim(0, 100) +
      labs(x = "Day", y = "% of body burden in depot") + THEME
  })
  output$p_tiw <- renderPlot({
    ggplot(sim(), aes(time, TIW)) + geom_line(linewidth = 1, colour = "#2e7d32") +
      labs(x = "Day", y = "Cumulative days in 14–20 mg/L") + THEME
  })
  output$p_clmit <- renderPlot({
    s <- sim()
    ggplot(s, aes(time, 48 * input$clmult * (1 + 0.35 * (ENZ - 1)))) +
      geom_line(linewidth = 1, colour = "#5e35b1") +
      labs(x = "Day", y = "Mitotane CL (L/day)") +
      ggtitle("Autoinduction: Css drifts DOWN on a fixed dose") + THEME
  })

  ## ---- tab 3 -------------------------------------------------------------
  output$p_cort <- renderPlot({
    s <- sim_ster()
    d <- s %>% transmute(time, `total cortisol` = CORT,
                         `free cortisol x20` = CFREEo * 20) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.7) +
      labs(x = "Day", y = "ug/dL (free scaled x20)", colour = NULL) + THEME
  })
  output$p_ratio <- renderPlot({
    ggplot(sim_ster(), aes(time, CORT / CFREEo)) +
      geom_line(linewidth = 0.8, colour = "#c62828") +
      geom_hline(yintercept = 16.9, linetype = "dashed") +
      labs(x = "Day", y = "total / free ratio") +
      ggtitle("Dashed = ratio without mitotane. The assay drifts, not the patient.") + THEME
  })
  output$p_adrn <- renderPlot({
    d <- sim_ster() %>% transmute(time, `normal adrenal mass` = ADRN,
                                  `tumour secretory capacity` = STCAP) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = "relative", colour = NULL) + THEME
  })
  output$p_ster <- renderPlot({
    d <- sim_ster() %>% transmute(time, ACTH, DHEAS, OHP17, ALDO) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(x = "Day", y = "relative to normal", colour = NULL) + THEME
  })
  output$t_cort <- renderTable({
    s <- sim_ster(); t1 <- max(s$time)
    w <- s$time >= t1 - 60
    data.frame(
      Metric = c("mean total cortisol (ug/dL)", "mean FREE cortisol (ug/dL)",
                 "total/free ratio", "CBG (fold)", "adrenal mass remaining (%)",
                 "cumulative adrenal-crisis hazard"),
      Value = c(sprintf("%.2f", mean(s$CORT[w])), sprintf("%.3f", mean(s$CFREEo[w])),
                sprintf("%.1f", mean(s$CORT[w]) / mean(s$CFREEo[w])),
                sprintf("%.2f", mean(s$CBG[w])),
                sprintf("%.1f", 100 * s$ADRN[nrow(s)]),
                sprintf("%.2f", s$AIHAZ[nrow(s)])),
      `Model normal` = c("4.99", "0.225", "16.9", "1.00", "100", "-"),
      check.names = FALSE)
  }, striped = TRUE)

  ## ---- tab 4 -------------------------------------------------------------
  output$p_enz <- renderPlot({
    ggplot(sim(), aes(time, ENZ)) + geom_line(linewidth = 1, colour = "#5e35b1") +
      geom_hline(yintercept = 1, linetype = "dotted") +
      labs(x = "Day", y = "CYP3A4 activity (fold of baseline)") + THEME
  })
  output$p_indfold <- renderPlot({
    ggplot(sim(), aes(time, INDFOLD)) + geom_line(linewidth = 1, colour = "#c62828") +
      geom_hline(yintercept = 1, linetype = "dotted") +
      labs(x = "Day", y = "Etoposide clearance (fold)") + THEME
  })
  output$p_counter <- renderPlot({
    ev_ <- comb(mit_ev(), ev_edp(input$bsa, 8), ev_hc(50, 400))
    d <- bind_rows(lapply(c(1, 0), function(iv) {
      p <- c(base_par(), list(INDVIC = iv))
      x <- do.call(run, c(list(events = ev_, end = 300, delta = 0.1), p))
      x %>% transmute(time, CETOo, TUMTOTo,
                      arm = ifelse(iv == 1, "victim induction ON (reality)",
                                            "induction OFF (counterfactual)"))
    }))
    p1 <- ggplot(d, aes(time, CETOo, colour = arm)) + geom_line(linewidth = 0.6) +
      coord_cartesian(xlim = c(140, 175)) +
      scale_colour_manual(values = c("#c62828", "#2e7d32")) +
      labs(x = "Day", y = "Etoposide (mg/L)", colour = NULL) +
      ggtitle("Cycle 6 etoposide exposure") + THEME
    p2 <- ggplot(d, aes(time, TUMTOTo, colour = arm)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = c("#c62828", "#2e7d32")) +
      labs(x = "Day", y = "Tumour (mL)", colour = NULL) +
      ggtitle("Consequence for tumour burden") + THEME
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  })

  ## ---- tab 5 -------------------------------------------------------------
  output$p_tum <- renderPlot({
    d <- sim_tum() %>% transmute(time, sensitive = TUMS, resistant = TUMR,
                                 total = TUMTOTo) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = c(sensitive = "#8b3a2c", resistant = "#e8a89c",
                                     total = "#1f4e79")) +
      labs(x = "Day", y = "Volume (mL)", colour = NULL) + THEME
  })
  output$p_sld <- renderPlot({
    s <- sim_tum(); b <- s$SLDo[1]
    ggplot(s, aes(time, SLDo)) + geom_line(linewidth = 1, colour = "#1f4e79") +
      geom_hline(yintercept = b, linetype = "dotted") +
      geom_hline(yintercept = b * 0.7, linetype = "dashed", colour = "#2e7d32") +
      labs(x = "Day", y = "Sum of diameters (mm)") +
      ggtitle("Dashed green = −30% (PR threshold)") + THEME
  })
  output$t_end <- renderTable({
    s <- sim_tum(); b <- s$SLDo[1]
    nad <- min(s$SLDo); tnad <- s$time[which.min(s$SLDo)]
    prog <- which(s$SLDo >= nad * 1.2 & s$SLDo >= nad + 5 & s$time > tnad)
    back <- which(s$SLDo >= b & s$time > tnad)
    data.frame(
      Endpoint = c("baseline SLD (mm)", "nadir SLD (mm)", "best change (%)",
                   "RECIST category", "PFS by nadir rule (d)",
                   "time to regain baseline (d)", "resistant clone at end (%)"),
      Value = c(sprintf("%.1f", b), sprintf("%.1f", nad),
                sprintf("%+.1f", 100 * (nad - b) / b),
                if (100 * (nad - b) / b <= -30) "PR" else
                  if (max(s$SLDo) >= b * 1.2 + 5) "PD" else "SD",
                if (length(prog)) sprintf("%.0f", s$time[prog[1]]) else "not reached",
                if (length(back)) sprintf("%.0f", s$time[back[1]]) else "not reached",
                sprintf("%.1f", 100 * s$TUMR[nrow(s)] / s$TUMTOTo[nrow(s)])))
  }, striped = TRUE)

  ## ---- tab 6 -------------------------------------------------------------
  scen_data <- reactive({
    shiny::req(input$scensel)
    bind_rows(lapply(input$scensel, function(nm) {
      sc <- scenarios[[nm]]
      x <- run_scenario(sc)
      x$scenario <- sc$label
      x
    }))
  })
  output$p_scen <- renderPlot({
    d <- scen_data()
    ggplot(d, aes(time, .data[[input$scenvar]], colour = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Day", y = input$scenvar, colour = NULL) +
      guides(colour = guide_legend(ncol = 1)) + THEME
  })
  output$t_scen <- renderTable({
    scen_data() %>% group_by(scenario) %>%
      summarise(`Cp d180` = round(CMITo[which.min(abs(time - 180))], 1),
                `tumour d180 (mL)` = round(TUMTOTo[which.min(abs(time - 180))]),
                `mean free cortisol` = round(mean(CFREEo[time > 100]), 3),
                `days in window` = round(max(TIW)), .groups = "drop")
  }, striped = TRUE)

  ## ---- tab 7 -------------------------------------------------------------
  fat_exp <- reactive({
    bind_rows(lapply(c(10, 22, 45), function(f) {
      x <- run(ev_mit_sched(sched_highdose), end = 540, delta = 2,
               FATKG0 = f, SECRETOR = 0)
      x$fat <- paste0(f, " kg fat")
      x
    }))
  })
  output$p_fatexp <- renderPlot({
    ggplot(fat_exp(), aes(time, CMITo, colour = fat)) + WINDOW +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c("#2e7d32", "#1f4e79", "#8a6d1f")) +
      labs(x = "Day", y = "Plasma mitotane (mg/L)", colour = NULL) +
      ggtitle("Same 6/4/3 g regimen; only adipose mass differs") + THEME
  })
  output$t_fatexp <- renderTable({
    fat_exp() %>% group_by(fat) %>%
      summarise(`Vss (L)` = round(VSSo[1]),
                `t→14 (d)` = ifelse(any(CMITo >= 14),
                                    as.character(time[which(CMITo >= 14)[1]]), "never"),
                `peak` = round(max(CMITo), 1),
                `Cp d540` = round(CMITo[which.max(time)], 1),
                `days >20` = round(sum(ABOVE20) * 2), .groups = "drop")
  }, striped = TRUE)
  output$p_fatcurve <- renderPlot({
    fats <- seq(6, 55, by = 4)
    d <- bind_rows(lapply(fats, function(f) {
      x <- run(ev_mit_sched(sched_highdose), end = 540, delta = 2,
               FATKG0 = f, SECRETOR = 0)
      data.frame(fat = f,
                 ttt = ifelse(any(x$CMITo >= 14), x$time[which(x$CMITo >= 14)[1]], NA),
                 css = x$CMITo[which.max(x$time)])
    }))
    d2 <- pivot_longer(d, -fat)
    ggplot(d2, aes(fat, value)) + geom_line(linewidth = 1, colour = "#1f4e79") +
      geom_point() + facet_wrap(~name, scales = "free_y",
        labeller = as_labeller(c(ttt = "time to 14 mg/L (days) — MOVES",
                                 css = "plateau Cp (mg/L) — DOES NOT"))) +
      labs(x = "Adipose mass (kg)", y = NULL) + THEME
  })

  ## ---- tab 8 -------------------------------------------------------------
  pop <- eventReactive(input$runpop, {
    withProgress(message = "Simulating virtual population...",
                 vpop_time_to_target(n = input$npop))
  })
  output$p_pop <- renderPlot({
    d <- pop()
    ggplot(d, aes(cp90, fill = regimen)) +
      geom_histogram(bins = 26, position = "identity", alpha = 0.6) +
      geom_vline(xintercept = 14, linetype = "dashed", colour = "#2e7d32") +
      scale_fill_manual(values = c(high = "#1f4e79", low = "#c0a020")) +
      labs(x = "Plasma mitotane at day 90 (mg/L)", y = "patients", fill = NULL) + THEME
  })
  output$t_pop <- renderTable({
    d <- pop()
    fit <- lm(cp90 ~ regimen + fat + clmult, data = d)
    ss <- anova(fit); tot <- sum(ss[["Sum Sq"]])
    data.frame(Term = rownames(ss),
               `Variance explained (%)` = round(100 * ss[["Sum Sq"]] / tot, 1),
               check.names = FALSE)
  }, striped = TRUE)
  output$popnote <- renderUI({
    d <- pop(); zh <- d[d$regimen == "high", ]
    q <- quantile(zh$fat, c(.1, .9))
    HTML(sprintf("<div class='note'>Regimen contrast on median Cp90:
      <b>%.1f mg/L</b>. Adipose P10→P90 contrast within one arm:
      <b>%.1f mg/L</b>. Correlation of Cp90 with adipose mass: <b>r = %.2f</b>.
      <br/>Both levers matter; only one of them was randomised.</div>",
      median(d$cp90[d$regimen == "high"]) - median(d$cp90[d$regimen == "low"]),
      median(zh$cp90[zh$fat <= q[1]]) - median(zh$cp90[zh$fat >= q[2]]),
      cor(zh$fat, zh$cp90)))
  })

  ## ---- tab 9 -------------------------------------------------------------
  imm_data <- reactive({
    bind_rows(lapply(list(
      list("pembro, non-secreting", ev_pem(16), 0),
      list("pembro, secreting",     ev_pem(16), 1),
      list("pembro + mitotane",     comb(ev_pem(16), ev_mit_sched(sched_highdose),
                                         ev_hc(50, 400)), 1)),
      function(z) {
        x <- run(z[[2]], end = 360, delta = 1, SECRETOR = z[[3]])
        x$arm <- z[[1]]; x
      }))
  })
  output$p_imm <- renderPlot({
    d <- imm_data() %>% select(time, arm, TEFF, TREG) %>% pivot_longer(c(TEFF, TREG))
    ggplot(d, aes(time, value, colour = arm, linetype = name)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Day", y = "relative activity", colour = NULL, linetype = NULL) + THEME
  })
  output$p_pd1 <- renderPlot({
    d <- imm_data()
    ggplot(d, aes(time, TUMTOTo, colour = arm)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = "Tumour (mL)", colour = NULL) +
      ggtitle(sprintf("PD-1 occupancy stays at %.1f%% in every arm",
                      100 * mean(d$OCCPD1o[d$time > 30]))) + THEME
  })
  output$p_igf <- renderPlot({
    d <- bind_rows(
      run(ev_lin(300), end = 300, delta = 1, SECRETOR = 0) %>%
        mutate(arm = "linsitinib (IR-A spared)"),
      run(ev_lin(300), end = 300, delta = 1, SECRETOR = 0, IC50RA = 0.25) %>%
        mutate(arm = "counterfactual: IR-A blocked too"),
      run(NULL, end = 300, delta = 1, SECRETOR = 0) %>% mutate(arm = "no drug"))
    ggplot(d, aes(time, IGFSIG, colour = arm)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = "IGF/AKT/mTOR signal (rel)", colour = NULL) +
      ggtitle("Target engaged, pathway restored") + THEME
  })
  output$p_glu <- renderPlot({
    d <- run(ev_lin(300), end = 300, delta = 1, SECRETOR = 0) %>%
      transmute(time, glucose = GLU, `insulin x40` = INS * 40) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = "mg/dL (insulin scaled)", colour = NULL) + THEME
  })

  ## ---- tab 10 ------------------------------------------------------------
  sim_safe <- reactive({
    do.call(run, c(list(events = comb(mit_ev(), ev_edp(input$bsa, 8), ev_hc(50, 400)),
                        end = min(input$tend, 400), delta = 0.25), base_par()))
  })
  output$p_ntox <- renderPlot({
    s <- sim_safe()
    ggplot(s, aes(time)) +
      geom_line(aes(y = NTOX, colour = "CNS injury score"), linewidth = 1) +
      geom_line(aes(y = CMITo / 40, colour = "plasma mitotane / 40"), linewidth = 0.8) +
      geom_hline(yintercept = 0.5, linetype = "dotted") +
      labs(x = "Day", y = NULL, colour = NULL) + THEME
  })
  output$p_liver <- renderPlot({
    d <- sim_safe() %>% transmute(time, `ALT (U/L)` = ALT, `free T4 x30` = FT4 * 30) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = NULL, colour = NULL) + THEME
  })
  output$p_anc <- renderPlot({
    ggplot(sim_safe(), aes(time, CIRCN)) + geom_line(linewidth = 0.7, colour = "#c62828") +
      geom_hline(yintercept = 1.5, linetype = "dashed") +
      geom_hline(yintercept = 0.5, linetype = "dotted") +
      labs(x = "Day", y = "ANC (10^9/L)") +
      ggtitle("Dashed 1.5 = grade 2; dotted 0.5 = grade 4") + THEME
  })
  output$p_organ <- renderPlot({
    d <- sim_safe() %>% transmute(time, `LVEF (%)` = LVEF,
                                  `eGFR (mL/min)` = GFR) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = NULL, colour = NULL) + THEME
  })
  output$t_safety <- renderTable({
    s <- sim_safe()
    data.frame(
      Readout = c("peak plasma mitotane (mg/L)", "days above 20 mg/L",
                  "days in 14–20 window", "peak CNS injury score",
                  "peak ALT (U/L)", "free T4 nadir (ng/dL)",
                  "ANC nadir (10^9/L)", "LVEF at end (%)",
                  "cumulative anthracycline (mg/m²)", "eGFR at end",
                  "adrenal cortex remaining (%)", "BMD at end (rel)"),
      Value = c(sprintf("%.1f", max(s$CMITo)), sprintf("%.0f", sum(s$ABOVE20) * 0.25),
                sprintf("%.0f", sum(s$INWIN) * 0.25), sprintf("%.3f", max(s$NTOX)),
                sprintf("%.0f", max(s$ALT)), sprintf("%.2f", min(s$FT4)),
                sprintf("%.2f", min(s$CIRCN)), sprintf("%.1f", s$LVEF[nrow(s)]),
                sprintf("%.0f", max(s$DOXCUM)), sprintf("%.0f", s$GFR[nrow(s)]),
                sprintf("%.1f", 100 * s$ADRN[nrow(s)]),
                sprintf("%.3f", s$BMD[nrow(s)])))
  }, striped = TRUE)

  ## ---- tab 11 ------------------------------------------------------------
  output$t_calib <- renderTable({
    s <- run(ev_mit_sched(sched_highdose), end = 480, delta = 1, SECRETOR = 0)
    e <- run(comb(ev_mit_sched(sched_highdose), ev_edp(1.85, 10), ev_hc(50, 480)),
             end = 300, delta = 1)
    data.frame(
      Target = c("Apparent Vss (L)", "Terminal t½ (days)",
                 "Therapeutic window (mg/L)", "Etoposide CL fold-increase",
                 "CBG fold rise", "SHBG fold rise",
                 "Hydrocortisone requirement", "Free cortisol fraction (normal)",
                 "Adrenal cortex loss on therapy"),
      Literature = c("~6,086 (CV 81.5%)", "18–160 (median ~50–60)", "14–20",
                     "~2 (inferred)", "~2", "~2", "roughly doubles",
                     "~4–5%", "near-complete"),
      Model = c(sprintf("%.0f", s$VSSo[1]),
                sprintf("%.0f–%.0f", log(2) * s$VSSo[1] / (48 * 2.05),
                        log(2) * s$VSSo[1] / 48),
                "14–20 (hard-coded readout)",
                sprintf("%.2f", max(e$INDFOLD)),
                sprintf("%.2f", max(s$CBG)), sprintf("%.2f", max(s$SHBG)),
                "20 → 50 mg/day", sprintf("%.1f%%", 100 * s$FFRAC[1]),
                sprintf("%.1f%% remaining", 100 * min(s$ADRN))),
      Source = c("Arshad 2018", "Corso 2021", "Terzolo 2013", "Kroiss 2011 + f_m",
                 "Nader 2006", "Nader 2006", "Chortis 2013", "textbook",
                 "Daffara 2008"),
      check.names = FALSE)
  }, striped = TRUE)
}

shinyApp(ui, server)
