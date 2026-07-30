# =============================================================================
# lch_shiny_app.R
# Interactive dashboard for the Langerhans Cell Histiocytosis QSP model
# -----------------------------------------------------------------------------
# Run with:  shiny::runApp("lch_shiny_app.R")
# Requires:  shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
#
# The app sources lch_mrgsolve_model.R for the model, phenotypes, regimen
# builders and the run_scenario() driver, so there is exactly one definition of
# the model in this directory.
#
# Eleven tabs:
#   1  환자 프로파일        Patient / cell-of-origin profile
#   2  약물 PK              Drug pharmacokinetics
#   3  MAPK 신호            pERK, cyclin D1, BCL2A1, SASP
#   4  세포 구획 · 병변     Precursor reservoir and lesional burden
#   5  사이토카인 · 분비체  Secretome
#   6  바이오마커           cfDNA, sCD163, CRP, ferritin
#   7  임상 엔드포인트      DAS, response, reactivation
#   8  영구 후유증          Irreversible sequelae (CDI, ND, biliary, lung, bone)
#   9  독성                 ANC (Friberg), skin, LVEF
#  10  시나리오 비교        Scenario comparison
#  11  구조적 선택 검증     Falsification / kill-switch panel
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("lch_mrgsolve_model.R", local = TRUE)

THEME <- theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey92"),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12))

PHENO_LABEL <- c(
  SSb      = "SS-b 단일계통 골 병변 (tissue-restricted origin, FSR 0.10)",
  MSROneg  = "MS RO-negative 다장기 (위험장기 비침범)",
  MSROpos  = "MS RO-positive 다장기 (HSC 기원, 위험장기 침범)",
  CNSrisk  = "두개안면 CNS-risk + 뇌하수체 침범",
  PLCH     = "성인 폐 LCH (흡연자)"
)

SCEN_LABEL <- setNames(
  vapply(scenarios, function(s) s$label, character(1)),
  names(scenarios)
)

long_of <- function(df, cols) {
  df %>% select(time, all_of(cols)) %>% pivot_longer(-time)
}

# =============================================================================
# UI
# =============================================================================
ui <- fluidPage(
  titlePanel("랑게르한스 세포 조직증식증 (LCH) — QSP 모델 대시보드"),
  p(tags$em(paste(
    "62개 ODE · 기원세포 파티션 · ERK 구동 분비체 · 세포정지성 MAPK 억제와",
    "중단 후 재발 · 시간적분으로서의 영구 후유증"))),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("scen", "시나리오 (Scenario)",
                  choices = names(SCEN_LABEL), selected = "S5_vem_continuous"),
      helpText(textOutput("scen_desc")),
      hr(),
      h4("환자 · 기원세포"),
      selectInput("geno", "드라이버 유전형",
                  choices = c("BRAF V600E" = 1, "MAP2K1 / ARAF / 음성" = 2),
                  selected = 1),
      sliderInput("wt", "체중 (kg)", 5, 80, 12, step = 1),
      sliderInput("fsr", "FSR — 기원세포 자기재생능", 0, 1, 1.0, step = 0.05),
      sliderInput("precm0", "PRECM0 — 초기 변이 전구세포 저장고",
                  0, 1.5, 0.45, step = 0.05),
      sliderInput("thr", "THR — 위험장기 seeding 분율", 0, 0.7, 0.34,
                  step = 0.02),
      sliderInput("thp", "THP — 뇌하수체 seeding 분율", 0, 0.5, 0.08,
                  step = 0.02),
      checkboxInput("smoke", "흡연 노출 (폐 LCH 구동)", FALSE),
      hr(),
      h4("약력학 · 구조 스위치"),
      sliderInput("ic50vem", "IC50 vemurafenib (유리, mg/L)",
                  0.005, 0.08, 0.018, step = 0.001),
      sliderInput("kprol", "KPROL — 병변 내 증식률 (1/일)",
                  0, 0.06, 0.022, step = 0.002),
      helpText(HTML(paste(
        "KPROL이 <b>0.0288</b>(최소 면역청소율)을 넘으면 병변이 전구세포",
        "공급 없이도 자립하게 되며, 낮은 Ki-67(≈9%)과 모순됩니다."))),
      sliderInput("kill", "SL_MAPKI_KILL — 반증 스위치",
                  0, 0.4, 0.0, step = 0.05),
      helpText(HTML(paste(
        "0이 아니면 MAPK 억제제가 <b>세포살상</b>이 되어 중단 후 재발이",
        "사라집니다 — 문헌과 모순."))),
      sliderInput("fnapo", "FN_APO — 니치의 세포자멸 보호",
                  0, 1, 0.25, step = 0.05),
      sliderInput("pniche", "PNICHE — 클론 소멸 문턱",
                  0.0005, 0.02, 0.003, step = 0.0005),
      sliderInput("kavp", "KAVP — AVP 뉴런 소실률 (1/일)",
                  0, 0.04, 0.010, step = 0.002),
      hr(),
      actionButton("go", "다시 계산 (Simulate)", class = "btn-primary"),
      br(), br(),
      helpText("슬라이더는 표현형 기본값을 덮어씁니다. 시나리오를 바꾸면",
               "해당 표현형 값으로 초기화하려면 페이지를 새로고침하세요.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. 환자 프로파일",
                 h4("기원세포 파티션 — 이 모델의 첫 번째 구조적 선택"),
                 plotOutput("p_partition", height = "300px"),
                 h4("진단 시점 상태 (run-in 종료 시)"),
                 DTOutput("t_baseline")),
        tabPanel("2. 약물 PK",
                 plotOutput("p_pk", height = "620px"),
                 helpText("PK 파라미터는 70 kg 기준값이며 체중에 대해",
                          "allometric 스케일(CL^0.75, V^1.0)로 보정됩니다.")),
        tabPanel("3. MAPK 신호",
                 plotOutput("p_mapk", height = "600px"),
                 helpText(HTML(paste(
                   "pERK가 기저의 <b>20% 미만</b>(=80% 초과 억제)으로 내려가야",
                   "cyclin D1이 급격히 떨어집니다 — Hill 계수 2가 만드는",
                   "문턱입니다 (Bollag 2010, PMID 20823850)。")))),
        tabPanel("4. 세포 구획 · 병변",
                 plotOutput("p_cells", height = "620px"),
                 helpText(HTML(paste(
                   "핵심 대비: 병변 질량은 사라지는데 <b>골수 전구세포 저장고</b>는",
                   "남습니다. 남은 저장고가 중단 후 재발의 원천입니다.")))),
        tabPanel("5. 사이토카인 · 분비체",
                 plotOutput("p_cyto", height = "620px")),
        tabPanel("6. 바이오마커",
                 plotOutput("p_bio", height = "560px"),
                 helpText("cfDNA는 변이세포의 사멸 flux를 따르므로, 질량이",
                          "남아 있어도 정지성 약물에서는 완만한 플래토를",
                          "만듭니다.")),
        tabPanel("7. 임상 엔드포인트",
                 plotOutput("p_das", height = "420px"),
                 h4("반응 · 재활성화 요약"),
                 DTOutput("t_response")),
        tabPanel("8. 영구 후유증",
                 plotOutput("p_seq", height = "620px"),
                 helpText(HTML(paste(
                   "다섯 개 풀(AVPN·ANTPIT·NEUR·BILF·LUNGC)은 <b>단조",
                   "(monotone)</b>입니다 — 어떤 치료도 되돌리지 못합니다.",
                   "따라서 TTET(활성 질환 일수)이 약효보다 후유증을 더 잘",
                   "설명합니다.")))),
        tabPanel("9. 독성",
                 plotOutput("p_tox", height = "620px")),
        tabPanel("10. 시나리오 비교",
                 checkboxGroupInput(
                   "cmp", "비교할 시나리오",
                   choices = names(SCEN_LABEL),
                   selected = c("S5_vem_continuous", "S6_vem_stop",
                                "S4_cladarac_upfront",
                                "S8_bridge_consolidate"),
                   inline = TRUE),
                 actionButton("go_cmp", "비교 실행", class = "btn-primary"),
                 plotOutput("p_cmp", height = "560px"),
                 DTOutput("t_cmp")),
        tabPanel("11. 구조적 선택 검증",
                 h4("반증 패널 (Falsification panel)"),
                 p(HTML(paste(
                   "각 구조적 선택에는 <b>kill switch</b>가 있습니다.",
                   "스위치를 켰을 때 모델이 문헌과 모순되는 거동을 보여야",
                   "합니다 — 그것이 그 선택이 실제로 일을 하고 있다는 증거입니다."))),
                 actionButton("go_kill", "kill switch 3종 실행",
                              class = "btn-primary"),
                 plotOutput("p_kill", height = "620px"),
                 verbatimTextOutput("t_kill"))
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  output$scen_desc <- renderText(SCEN_LABEL[[input$scen]])

  overrides <- reactive({
    list(GENO = as.numeric(input$geno), WT = input$wt, FSR = input$fsr,
         PRECM0 = input$precm0, THR = input$thr, THP = input$thp,
         SMOKE = as.numeric(input$smoke), IC50_VEM = input$ic50vem,
         KPROL = input$kprol, SL_MAPKI_KILL = input$kill,
         FN_APO = input$fnapo, PNICHE = input$pniche, KAVP = input$kavp)
  })

  sim <- eventReactive(input$go, {
    withProgress(message = "적분 중 (62 ODE)...", value = 0.5, {
      run_scenario(input$scen, overrides())
    })
  }, ignoreNULL = FALSE)

  # ---- 1. patient profile -------------------------------------------------
  output$p_partition <- renderPlot({
    p <- modifyList(pheno[[scenarios[[input$scen]]$kind]], overrides())
    th <- c(Bone = p$THB, Skin = p$THS, `Risk organs` = p$THR,
            Pituitary = p$THP, CNS = p$THC, Lung = p$THL)
    th[is.na(th)] <- 0
    data.frame(site = names(th), theta = as.numeric(th)) %>%
      ggplot(aes(reorder(site, theta), theta, fill = site)) +
      geom_col(show.legend = FALSE) +
      coord_flip() +
      labs(title = paste0("Seeding partition  (FSR = ", p$FSR,
                          ",  PRECM0 = ", p$PRECM0, ")"),
           x = NULL, y = "theta (fraction of circulating progeny)") +
      THEME
  })

  output$t_baseline <- renderDT({
    d <- sim()
    b <- d[1, ]
    data.frame(
      quantity = c("DAS at diagnosis", "Total lesional burden (units)",
                   "Risk-organ burden", "Pituitary burden", "CNS burden",
                   "cfDNA signal", "AVP neuron pool", "Neuron pool",
                   "Biliary fibrosis", "ANC (1e9/L)"),
      value = round(c(b$DAS, b$LTOT_units, b$LRO, b$LPIT, b$LCNS, b$CFDNA,
                      b$AVPN, b$NEUR, b$BILF, b$ANC), 4)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # ---- 2. PK --------------------------------------------------------------
  output$p_pk <- renderPlot({
    long_of(sim(), c("CVEM_mgL", "CDAB_mgL", "CTRA_ngmL", "CVBL_mgL",
                     "CPRE_mgL", "ARATP", "CLATP", "VEMI")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.5) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(title = "약물 농도 · 세포내 활성대사체 · CYP3A4 자가유도",
           x = "진단 후 일수", y = NULL) + THEME
  })

  # ---- 3. MAPK ------------------------------------------------------------
  output$p_mapk <- renderPlot({
    d <- sim()
    long_of(d, c("pERK_pct", "CCND", "BCL", "SASP", "AUCERK")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      geom_hline(data = data.frame(name = "pERK_pct", y = 20),
                 aes(yintercept = y), linetype = 2, colour = "red") +
      labs(title = "MAPK 신호 출력 (붉은 선 = 80% 억제 문턱)",
           x = "진단 후 일수", y = NULL) + THEME
  })

  # ---- 4. cells -----------------------------------------------------------
  output$p_cells <- renderPlot({
    d <- sim()
    long_of(d, c("PRECM", "CIRC", "LBONE", "LSKIN", "LRO", "LPIT", "LCNS",
                 "LLUNG", "TREG", "OCL")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "전구세포 저장고와 장기별 병변 부하",
           x = "진단 후 일수", y = "burden units") + THEME
  })

  # ---- 5. secretome -------------------------------------------------------
  output$p_cyto <- renderPlot({
    long_of(sim(), c("IL1B", "TNFA", "IL6", "OSM", "MMP9", "RANKL")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(title = "ERK/SASP 구동 분비체", x = "진단 후 일수", y = NULL) +
      THEME
  })

  # ---- 6. biomarkers ------------------------------------------------------
  output$p_bio <- renderPlot({
    d <- sim()
    lod <- as.numeric(param(mod)$CF_LOD)
    long_of(d, c("CFDNA", "SCD163", "CRP", "FERR")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      geom_hline(data = data.frame(name = "CFDNA", y = lod),
                 aes(yintercept = y), linetype = 2, colour = "red") +
      labs(title = "바이오마커 (붉은 선 = ddPCR 검출한계)",
           x = "진단 후 일수", y = NULL) + THEME
  })

  # ---- 7. endpoints -------------------------------------------------------
  output$p_das <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, DAS)) +
      geom_line(linewidth = 0.7) +
      geom_hline(yintercept = as.numeric(param(mod)$DAS_ACTIVE),
                 linetype = 2, colour = "red") +
      labs(title = "Disease Activity Score (붉은 선 = 활성 질환 문턱)",
           x = "진단 후 일수", y = "DAS") + THEME
  })

  output$t_response <- renderDT({
    d <- sim()
    d0 <- d$DAS[1]
    at_t <- function(tt) d[which.min(abs(d$time - tt)), ]
    w6 <- at_t(42); w12 <- at_t(84); end <- d[nrow(d), ]
    resp <- function(x) {
      if (x$DAS < 0.1 * d0) "NAD / better"
      else if (x$DAS < d0) "intermediate"
      else "worse"
    }
    data.frame(
      timepoint = c("baseline", "week 6", "week 12", "end of follow-up"),
      DAS = round(c(d0, w6$DAS, w12$DAS, end$DAS), 2),
      cfDNA = round(c(d$CFDNA[1], w6$CFDNA, w12$CFDNA, end$CFDNA), 4),
      reservoir = round(c(d$PRECM[1], w6$PRECM, w12$PRECM, end$PRECM), 4),
      category = c("-", resp(w6), resp(w12), resp(end)),
      TTET_days = round(c(0, w6$TTET, w12$TTET, end$TTET), 0)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # ---- 8. sequelae --------------------------------------------------------
  output$p_seq <- renderPlot({
    d <- sim()
    long_of(d, c("AVPN", "ANTPIT", "NEUR", "BILF", "LUNGC", "BVOL",
                 "TTET", "CUMDAS")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(title = "비가역 장기 풀과 활성 질환 누적",
           x = "진단 후 일수", y = NULL) + THEME
  })

  # ---- 9. toxicity -------------------------------------------------------
  output$p_tox <- renderPlot({
    d <- sim()
    long_of(d, c("ANC", "PROL", "SKTOX", "LVEF")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      geom_hline(data = data.frame(name = "ANC", y = 0.5),
                 aes(yintercept = y), linetype = 2, colour = "red") +
      labs(title = "골수억제(Friberg) · 역설적 피부독성 · LVEF",
           x = "진단 후 일수", y = NULL) + THEME
  })

  # ---- 10. scenario comparison -------------------------------------------
  cmp <- eventReactive(input$go_cmp, {
    req(length(input$cmp) > 0)
    withProgress(message = "시나리오 비교 중...", value = 0.4, {
      bind_rows(lapply(input$cmp, function(s) run_scenario(s, overrides())))
    })
  })

  output$p_cmp <- renderPlot({
    cmp() %>%
      select(time, scenario, DAS, CFDNA, PRECM, AVPN, NEUR, BVOL) %>%
      pivot_longer(-c(time, scenario)) %>%
      ggplot(aes(time, value, colour = scenario)) +
      geom_line(linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "시나리오 비교", x = "진단 후 일수", y = NULL) + THEME
  })

  output$t_cmp <- renderDT({
    cmp() %>% group_by(scenario) %>% slice_tail(n = 1) %>%
      transmute(scenario,
                DAS = round(DAS, 2), cfDNA = round(CFDNA, 4),
                reservoir = round(PRECM, 5),
                AVPN = round(AVPN, 3), NEUR = round(NEUR, 3),
                BILF = round(BILF, 2), LUNGC = round(LUNGC, 2),
                CDI, ND, TTET = round(TTET, 0)) %>%
      datatable(rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })

  # ---- 11. falsification panel -------------------------------------------
  kill <- eventReactive(input$go_kill, {
    withProgress(message = "kill switch 3종 실행 중...", value = 0.3, {
      list(
        c1_ss   = run_scenario("S1_observation"),
        c1_swap = run_scenario("S1_observation", pheno$MSROpos),
        c2_base = run_scenario("S6_vem_stop"),
        c2_kill = run_scenario("S6_vem_stop", list(SL_MAPKI_KILL = 0.25)),
        c3_late = run_scenario("S9_delayed_dx"),
        c3_early = run_scenario("S9b_early_dx")
      )
    })
  })

  output$p_kill <- renderPlot({
    k <- kill()
    bind_rows(
      k$c1_ss    %>% transmute(time, y = LRO,   panel = "① 파티션: SS-b 기원기술자",         arm = "baseline"),
      k$c1_swap  %>% transmute(time, y = LRO,   panel = "① 파티션: MS RO+ 기원기술자로 교체", arm = "switch"),
      k$c2_base  %>% transmute(time, y = CFDNA, panel = "② 세포정지성: SL_MAPKI_KILL = 0",     arm = "baseline"),
      k$c2_kill  %>% transmute(time, y = CFDNA, panel = "② 세포정지성: SL_MAPKI_KILL = 0.25",  arm = "switch"),
      k$c3_early %>% transmute(time, y = AVPN,  panel = "③ 시간적분: 14일에 치료 시작",        arm = "baseline"),
      k$c3_late  %>% transmute(time, y = AVPN,  panel = "③ 시간적분: 180일 지연",              arm = "switch")
    ) %>%
      ggplot(aes(time, y, colour = arm)) + geom_line(linewidth = 0.7) +
      facet_wrap(~panel, scales = "free_y", ncol = 2) +
      labs(title = "세 가지 구조적 선택의 kill switch",
           x = "진단 후 일수", y = NULL) + THEME
  })

  output$t_kill <- renderPrint({
    k <- kill()
    last <- function(d, v) d[[v]][nrow(d)]
    cat("① 기원세포 파티션만 교체 (동일한 lesional rate constants)\n")
    cat(sprintf("   위험장기 부하 d730:  SS-b = %.4f   ->  MS RO+ 기술자 = %.4f\n",
                last(k$c1_ss, "LRO"), last(k$c1_swap, "LRO")))
    cat("\n② MAPK 억제제에 살상항을 부여 (SL_MAPKI_KILL = 0.25)\n")
    cat(sprintf("   중단 후 DAS 최대치:  0 = %.2f   ->  0.25 = %.2f  (재발 소멸)\n",
                max(k$c2_base$DAS[k$c2_base$time > 365]),
                max(k$c2_kill$DAS[k$c2_kill$time > 365])))
    cat(sprintf("   저장고 d365:         0 = %.4f   ->  0.25 = %.6f\n",
                k$c2_base$PRECM[which.min(abs(k$c2_base$time - 365))],
                k$c2_kill$PRECM[which.min(abs(k$c2_kill$time - 365))]))
    cat("\n③ 치료 시작 시점만 변경 (14일 vs 180일)\n")
    cat(sprintf("   AVP 뉴런 풀 d730:    14일 = %.3f  ->  180일 = %.3f\n",
                last(k$c3_early, "AVPN"), last(k$c3_late, "AVPN")))
    cat(sprintf("   중추성 요붕증:       14일 = %d      ->  180일 = %d\n",
                as.integer(last(k$c3_early, "CDI")),
                as.integer(last(k$c3_late, "CDI"))))
    cat(sprintf("   활성 질환 일수 TTET: 14일 = %.0f     ->  180일 = %.0f\n",
                last(k$c3_early, "TTET"), last(k$c3_late, "TTET")))
  })
}

shinyApp(ui, server)
