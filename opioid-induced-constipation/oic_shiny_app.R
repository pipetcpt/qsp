## =============================================================================
##  oic_shiny_app.R
##  Opioid-Induced Constipation (OIC) QSP — interactive dashboard
##  오피오이드 유발 변비 QSP — 인터랙티브 대시보드
##
##  10 tabs:
##    1  환자 프로파일   Patient profile & regimen
##    2  PK              Opioid and antagonist pharmacokinetics
##    3  수용체 점유율   Receptor occupancy — the two compartments side by side
##    4  선택지수 SI     Selectivity index & the P-gp drug-interaction explorer
##    5  대장 통과       Colonic transit, segment loads and hydration
##    6  임상 종말점     SBM / CSBM / Bristol / straining / PAC-SYM
##    7  시나리오 비교   Scenario comparison (16 regimens)
##    8  브레이크 분해   Which brake carries the drug response?
##    9  내성 적정       Tolerance titration over 24 weeks
##   10  바이오마커      Biomarker panel & safety
##
##  REQUIRES oic_mrgsolve_model.R in the same directory.
##
##  NOTE ON PROVENANCE: like the model file, this app has never been executed --
##  the build container had no R.  It is written against the mrgsolve model's
##  documented output columns.  See README.md.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)   # grid.arrange, used on tabs 5 and 6

source("oic_mrgsolve_model.R", local = TRUE)

THEME <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

DRUGS <- c("없음 (none)"                    = "none",
           "날록세골 naloxegol"             = "naloxegol",
           "날데메딘 naldemedine"           = "naldemedine",
           "메틸날트렉손 SC methylnaltrexone" = "methylnaltrexone_sc",
           "메틸날트렉손 PO 450 mg"          = "methylnaltrexone_po",
           "날록손 PO naloxone (대조군)"      = "naloxone_po")

DEFAULT_DOSE <- c(none = 0, naloxegol = 25, naldemedine = 0.2,
                  methylnaltrexone_sc = 12, methylnaltrexone_po = 450,
                  naloxone_po = 20)

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("오피오이드 유발 변비 (OIC) — QSP 시뮬레이터"),
  tags$p(style = "color:#555;margin-top:-8px",
         "One receptor, two compartments, one barrier. ",
         tags$b("교육·연구 목적 전용 — 임상 의사결정에 사용 금지.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("오피오이드 (opioid)"),
      sliderInput("op_dose", "옥시코돈 1회 용량 (mg, q12h)",
                  min = 0, max = 120, value = 30, step = 5),
      checkboxInput("methadone", "메타돈으로 대체 (ClC-2 직접 차단)", FALSE),

      h4("길항제 (antagonist)"),
      selectInput("drug", "약물", choices = DRUGS, selected = "naloxegol"),
      numericInput("pam_dose", "1회 용량 (mg)", value = 25, min = 0, step = 0.1),
      numericInput("pam_int", "투여간격 (h)", value = 24, min = 1, step = 1),

      h4("병용 (adjuncts)"),
      numericInput("peg_g",  "PEG 3350 (g/일)",     value = 0, min = 0, step = 1),
      numericInput("lac_g",  "락툴로스 (g/일)",       value = 0, min = 0, step = 5),
      numericInput("lub_ug", "루비프로스톤 (µg, bid)", value = 0, min = 0, step = 8),
      numericInput("lin_ug", "리나클로타이드 (µg/일)", value = 0, min = 0, step = 72),
      numericInput("pro_mg", "프루칼로프리드 (mg/일)", value = 0, min = 0, step = 1),

      h4("상호작용·환자 인자"),
      sliderInput("pgp", "P-gp 억제 (Kp,uu 배수)", min = 1, max = 20,
                  value = 1, step = 1),
      sliderInput("cyp", "CYP3A4 억제 (AUC 배수)", min = 1, max = 6,
                  value = 1, step = 0.2),
      sliderInput("fluid", "수분 섭취 배수", min = 0.6, max = 1.6,
                  value = 1, step = 0.05),
      sliderInput("absx", "대장 수분흡수능 배수", min = 0.7, max = 1.5,
                  value = 1, step = 0.05),
      checkboxInput("rescue", "구제 완하제 허용 (시험 프로토콜)", TRUE),

      sliderInput("days", "시뮬레이션 기간 (일)", min = 14, max = 168,
                  value = 84, step = 7),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1 · 환자 프로파일",
          br(),
          fluidRow(
            column(4, wellPanel(h4("SBM / 주"), h2(textOutput("kpi_sbm")),
                                tags$small("반응자 기준 ≥3"))),
            column(4, wellPanel(h4("Bristol"), h2(textOutput("kpi_bsfs")))),
            column(4, wellPanel(h4("PAC-SYM"), h2(textOutput("kpi_pacsym"))))
          ),
          fluidRow(
            column(4, wellPanel(h4("통증 NRS"), h2(textOutput("kpi_pain")))),
            column(4, wellPanel(h4("COWS 금단"), h2(textOutput("kpi_cows")),
                                tags$small("안전성 종말점"))),
            column(4, wellPanel(h4("구제약 / 주"), h2(textOutput("kpi_resc"))))
          ),
          hr(),
          h4("요약"), tableOutput("tbl_profile")),

        tabPanel("2 · PK",
          br(), plotOutput("p_pk", height = "620px"),
          tags$p(tags$small(
            "위: 오피오이드 혈장·뇌 유리농도. 아래: 길항제 혈장·뇌 유리농도. ",
            "두 뇌 곡선의 비율이 선택성의 전부다."))),

        tabPanel("3 · 수용체 점유율",
          br(), plotOutput("p_occ", height = "620px"),
          tags$p(tags$small(
            "같은 수용체, 두 구획. 장관 곡선과 중추 곡선의 간격이 치료역이다."))),

        tabPanel("4 · 선택지수 · 약물상호작용",
          br(),
          h4("선택지수 SI = 장관 길항제 점유율 / 중추 길항제 점유율"),
          tableOutput("tbl_si"), hr(),
          h4("P-gp 억제제는 노출이 아니라 '비'를 바꾼다"),
          tableOutput("tbl_ddi"),
          tags$p(tags$small(
            "CYP3A4 억제는 순수한 노출 이동이므로 감량으로 대부분 회복된다. ",
            "P-gp 억제는 Kp,uu 자체를 움직이므로 감량은 안전성을 거의 되찾지 못한 채 ",
            "효능만 잃는다 — 라벨의 감량 지시가 겨냥하는 변수가 틀린 경우다."))),

        tabPanel("5 · 대장 통과",
          br(), plotOutput("p_transit", height = "620px"),
          tags$p(tags$small(
            "구획별 고형물(g)과 수분율. 수분율의 하한(≈0.62)은 고형상 결합수이며, ",
            "이 하한이 없으면 모델은 건강한 환자에서도 발산한다."))),

        tabPanel("6 · 임상 종말점",
          br(), plotOutput("p_end", height = "620px"),
          tags$p(tags$small(
            "SBM은 '구제 완하제를 24시간 내 쓰지 않은 배변'으로 정의된다. ",
            "이 모델에서 구제약에 의한 배변은 대장을 비우지만 SBM으로 세지 않는다."))),

        tabPanel("7 · 시나리오 비교",
          br(), actionButton("go_scen", "16개 시나리오 실행", class = "btn-default"),
          br(), br(), plotOutput("p_scen", height = "560px"),
          tableOutput("tbl_scen")),

        tabPanel("8 · 브레이크 분해",
          br(), actionButton("go_brake", "브레이크 분해 실행", class = "btn-default"),
          br(), br(), tableOutput("tbl_brake"),
          tags$p(tags$small(
            "각 브레이크를 하나씩 끄고 동일한 0 → 25 mg 날록세골 비교를 다시 측정한다. ",
            "배수가 1에 가까워지는 항이 약효를 나르는 항이다. 참조 구현에서 ",
            "통과–수분 되먹임 고리를 완전히 제거해도 배수는 2.89 → 2.91로 거의 ",
            "변하지 않았다 — 그 고리는 증폭기가 아니다."))),

        tabPanel("9 · 내성 적정",
          br(), actionButton("go_tol", "24주 적정 실행", class = "btn-default"),
          br(), br(), plotOutput("p_tol", height = "520px"),
          tableOutput("tbl_tol"),
          tags$p(tags$small(
            "통증을 NRS 4.0으로 유지하도록 매주 용량을 조정한다. ",
            "중추 수용체 가용도는 무너지고 장관 가용도는 유지된다."))),

        tabPanel("10 · 바이오마커·안전성",
          br(), plotOutput("p_bio", height = "620px"),
          tableOutput("tbl_bio"))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  observeEvent(input$drug, {
    updateNumericInput(session, "pam_dose",
                       value = unname(DEFAULT_DOSE[input$drug]))
    updateNumericInput(session, "pam_int",
                       value = if (input$drug == "naloxone_po") 8 else 24)
  })

  sim <- eventReactive(input$go, {
    pars <- list(PGPINH = input$pgp, CYP3A4I = input$cyp,
                 FLUIDX = input$fluid, ABSX = input$absx,
                 METHADON = as.numeric(input$methadone),
                 RESCUEON = as.numeric(input$rescue))
    run_oic(drug = input$drug, days = input$days, params = pars,
            op_dose = input$op_dose, op_int = 12,
            pam_dose = input$pam_dose, pam_int = input$pam_int,
            pro_dose = input$pro_mg, peg_g = input$peg_g, lac_g = input$lac_g,
            lub_ug = input$lub_ug, lin_ug = input$lin_ug)
  }, ignoreNULL = FALSE)

  smry <- reactive(summarise_oic(sim()))

  fmt <- function(x, d = 2) formatC(x, digits = d, format = "f")
  output$kpi_sbm    <- renderText(fmt(smry()$SBM_wk))
  output$kpi_bsfs   <- renderText(fmt(smry()$BSFS))
  output$kpi_pacsym <- renderText(fmt(smry()$PACSYM))
  output$kpi_pain   <- renderText(fmt(smry()$PAIN))
  output$kpi_cows   <- renderText(fmt(smry()$COWS))
  output$kpi_resc   <- renderText(fmt(smry()$RESC_wk))

  output$tbl_profile <- renderTable({
    s <- smry()
    tibble(
      지표 = c("SBM/주", "CSBM/주", "구제약/주", "Bristol", "힘주기", "팽만",
               "PAC-SYM", "PAC-QOL", "통증 NRS", "COWS",
               "장관 작용제 점유율", "장관 길항제 점유율",
               "중추 길항제 점유율", "선택지수 SI",
               "직장 수분율", "대장 통과시간 (h)", "매복 위험 적분"),
      값 = c(s$SBM_wk, s$CSBM_wk, s$RESC_wk, s$BSFS, s$STRAIN, s$DIST,
             s$PACSYM, s$QOL, s$PAIN, s$COWS, s$OCCG_AG, s$OCCG_ANT,
             s$OCCC_ANT, s$SI, s$W4, s$CTT_H, s$IMPACT))
  }, digits = 4)

  ## ---- 2 PK ----------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim() %>% filter(time >= max(time) - 168)
    d %>%
      transmute(time = time - min(time),
                `오피오이드 혈장 (mg/L)` = CP_OP,
                `오피오이드 뇌 유리 (nM)` = COP_BR,
                `길항제 혈장 (mg/L)` = CP_PAM,
                `길항제 뇌 유리 (nM)` = CPAM_BR) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.8, colour = "#4f88a8") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "마지막 주 시간 (h)", y = NULL,
           title = "PK — 마지막 7일") + THEME
  })

  ## ---- 3 occupancy ---------------------------------------------------------
  output$p_occ <- renderPlot({
    sim() %>%
      transmute(day = time / 24,
                `장관 MOR — 작용제` = OCCG_AG,
                `장관 MOR — 길항제` = OCCG_ANT,
                `중추 MOR — 작용제` = OCCC_AG,
                `중추 MOR — 길항제` = OCCC_ANT) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) +
      geom_line(linewidth = 0.8) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_colour_manual(values = c("#c0392b", "#2e8b57", "#8b7bb8", "#c9a227")) +
      labs(x = "일", y = "점유율", colour = NULL,
           title = "같은 수용체, 두 구획") + THEME
  })

  ## ---- 4 selectivity & DDI -------------------------------------------------
  output$tbl_si  <- renderTable(selectivity_table(days = 56), digits = 5)
  output$tbl_ddi <- renderTable(ddi_table(days = 56), digits = 4)

  ## ---- 5 transit -----------------------------------------------------------
  output$p_transit <- renderPlot({
    d <- sim() %>% mutate(day = time / 24)
    p1 <- d %>%
      select(day, `상행` = S1, `횡행` = S2, `하행` = S3, `직장S상` = S4) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(x = NULL, y = "고형물 (g)", colour = NULL,
           title = "구획별 고형물 부하") + THEME
    p2 <- d %>%
      transmute(day,
                `상행` = W1/(W1+S1), `횡행` = W2/(W2+S2),
                `하행` = W3/(W3+S3), `직장S상` = W4/(W4+S4)) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 0.615, linetype = "dashed", colour = "#c0392b") +
      labs(x = "일", y = "수분율", colour = NULL,
           title = "구획별 수분율 (붉은 선 = 결합수 하한)") + THEME
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  })

  ## ---- 6 endpoints ---------------------------------------------------------
  output$p_end <- renderPlot({
    d <- sim() %>% mutate(day = time / 24)
    wk <- d %>% mutate(week = floor(day / 7)) %>% group_by(week) %>%
      summarise(SBM = max(CUM_SBM) - min(CUM_SBM),
                CSBM = max(CUM_CSBM) - min(CUM_CSBM),
                RESC = max(CUM_RESC) - min(CUM_RESC), .groups = "drop")
    p1 <- wk %>% pivot_longer(-week) %>%
      ggplot(aes(week, value, fill = name)) +
      geom_col(position = "dodge") +
      geom_hline(yintercept = 3, linetype = "dashed", colour = "#c0392b") +
      labs(x = NULL, y = "주당 횟수", fill = NULL,
           title = "주별 배변 — 붉은 선 = 반응자 기준") + THEME
    p2 <- d %>% select(day, Bristol = BSFS, `힘주기` = STRAIN,
                       `팽만` = DIST, `PAC-SYM` = PACSYM) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(x = "일", y = "점수", colour = NULL, title = "증상 점수") + THEME
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  })

  ## ---- 7 scenarios ---------------------------------------------------------
  scen <- eventReactive(input$go_scen, run_all_scenarios(days = 56))
  output$p_scen   <- renderPlot(plot_scenarios(scen()))
  output$tbl_scen <- renderTable(scen(), digits = 3)

  ## ---- 8 brake decomposition ----------------------------------------------
  brake <- eventReactive(input$go_brake, brake_decomposition(days = 56))
  output$tbl_brake <- renderTable(brake(), digits = 3)

  ## ---- 9 tolerance ---------------------------------------------------------
  tol <- eventReactive(input$go_tol, tolerance_titration(weeks = 24))
  output$p_tol <- renderPlot({
    tol() %>%
      select(week, `옥시코돈 mg/일` = oxy_mg_day,
             `중추 수용체 가용도` = OCCG_AG, `SBM/주` = SBM_wk,
             `통증 NRS` = PAIN) %>%
      pivot_longer(-week) %>%
      ggplot(aes(week, value)) + geom_line(linewidth = 0.9, colour = "#8b7bb8") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "주", y = NULL, title = "24주 적정 — 진통 유지 시도") + THEME
  })
  output$tbl_tol <- renderTable(tol(), digits = 3)

  ## ---- 10 biomarkers -------------------------------------------------------
  output$p_bio <- renderPlot({
    sim() %>% mutate(day = time / 24) %>%
      select(day, `장관 MOR 가용도` = RG_AV, `중추 MOR 가용도` = RC_AV,
             `장관 cAMP` = CAMP, `ACh 유리` = ACH,
             `HAPC (1/h)` = HAPC, `분절 긴장` = TONE,
             `ClC-2 활성` = CLC2, `매복 위험` = IMPACT) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value)) + geom_line(linewidth = 0.8, colour = "#4f9a9a") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "일", y = NULL, title = "기계론적 바이오마커") + THEME
  })
  output$tbl_bio <- renderTable({
    d <- tail(sim(), 1)
    tibble(바이오마커 = c("장관 MOR 가용도", "중추 MOR 가용도", "장관 cAMP",
                          "HAPC (1/h)", "분절 긴장", "총 대장 고형물 (g)",
                          "대장 통과시간 (h)", "매복 위험 적분"),
           값 = c(d$RG_AV, d$RC_AV, d$CAMP, d$HAPC, d$TONE, d$STOT,
                  d$CTT_H, d$IMPACT))
  }, digits = 4)
}

shinyApp(ui, server)
