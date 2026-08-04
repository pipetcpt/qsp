## =============================================================================
##  Niemann-Pick disease type C (NPC) — QSP Shiny dashboard
##  니만-피크병 C형 · 인터랙티브 대시보드
##
##  9 tabs:
##    1. 환자 프로파일     patient profile: genotype x developmental vulnerability
##    2. 약물 PK           drug PK for all four agents
##    3. 리소좀 · 지질     lysosomal biology (the mechanism)
##    4. 소뇌 · 예비능     cerebellum and the damage-reserve integral
##    5. 임상 엔드포인트   clinical endpoints (SARA, NPCCSS 17/5/4, saccades)
##    6. 시나리오 비교     scenario comparison (18 regimens)
##    7. 바이오마커        biomarkers, and the compartment claim
##    8. 시험 설계 실험    trial-design experiment (the design claim)
##    9. 보정 대조표       model vs published targets, including the misses
##
##  Run:  shiny::runApp("npc_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

MODEL_FILE <- "npc_mrgsolve_model.R"
YR <- 365.25

mod_global <- mread_cache("npc", MODEL_FILE)

## ---------------------------------------------------------------------------
## genotype and archetype tables (kept in sync with the .R model's helpers)
## ---------------------------------------------------------------------------
npc_genotypes <- list(
  `WT (정상)`                  = list(f_null = 0.00, theta0 = 1.000, f_npc2 = 1.00),
  `I1061T/I1061T (청소년형)`   = list(f_null = 0.00, theta0 = 0.055, f_npc2 = 1.00),
  `I1061T/null (후기영아형)`   = list(f_null = 0.50, theta0 = 0.055, f_npc2 = 1.00),
  `null/null (주산기형)`       = list(f_null = 1.00, theta0 = 0.055, f_npc2 = 1.00),
  `mild/mild (성인형)`         = list(f_null = 0.00, theta0 = 0.322, f_npc2 = 1.00),
  `NPC2 질환`                  = list(f_null = 0.00, theta0 = 1.000, f_npc2 = 0.03)
)

DRUGS <- c("miglustat", "arimoclomol", "levacetylleucine", "IT adrabetadex")

make_events <- function(drugs, start_yr, dur_yr,
                        mig_mg = 200, ari_mg = 124, nal_g = 4,
                        cd_mg = 900, cd_weeks = 2) {
  t0 <- start_yr * YR
  dur <- dur_yr * YR
  es <- list()
  if ("miglustat" %in% drugs)
    es <- c(es, list(ev(time = t0, amt = mig_mg, cmt = "MIG_GUT",
                        ii = 1/3, addl = ceiling(dur * 3) - 1)))
  if ("arimoclomol" %in% drugs)
    es <- c(es, list(ev(time = t0, amt = ari_mg, cmt = "ARI_GUT",
                        ii = 1/3, addl = ceiling(dur * 3) - 1)))
  if ("levacetylleucine" %in% drugs)
    es <- c(es, list(ev(time = t0, amt = nal_g * 1000 / 3, cmt = "NAL_GUT",
                        ii = 1/3, addl = ceiling(dur * 3) - 1)))
  if ("IT adrabetadex" %in% drugs)
    es <- c(es, list(ev(time = t0, amt = cd_mg, cmt = "CD_CSF",
                        ii = 7 * cd_weeks,
                        addl = max(0, floor(dur / (7 * cd_weeks)) - 1))))
  if (!length(es)) return(ev(time = 0, amt = 0, cmt = "MIG_GUT"))
  Reduce(`+`, es)
}

run_sim <- function(geno, v_dev, liver_risk, years, events, delta = 7, ...) {
  g <- npc_genotypes[[geno]]
  m <- param(mod_global, f_null = g$f_null, theta0 = g$theta0,
             f_npc2 = g$f_npc2, v_dev = v_dev, liver_risk = liver_risk, ...)
  m %>% mrgsim(events = events, end = years * YR, delta = delta,
               hmax = 0.5, atol = 1e-8, rtol = 1e-8) %>% as_tibble()
}

at_age <- function(out, col, age) approx(out$AGE_YR, out[[col]], xout = age, rule = 2)$y

theme_npc <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 10, colour = "grey30"))
}

longify <- function(out, cols) {
  out %>% select(AGE_YR, all_of(cols)) %>%
    pivot_longer(-AGE_YR, names_to = "variable", values_to = "value")
}

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("니만-피크병 C형 (NPC) — QSP 대시보드 / Niemann-Pick type C QSP dashboard"),
  tags$p(style = "color:#555; margin-top:-8px;",
         HTML("혈장 표지자는 <b>내장</b> 저장을 읽고 병은 <b>소뇌</b>에 있다 · ",
              "기능 = 1 − D<sub>rev</sub> − D<sub>irr</sub> · ",
              "각 약은 자기 기전만 볼 수 있는 설계로 시험되었다")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① 환자 (Patient)"),
      selectInput("geno", "유전형 (genotype)",
                  choices = names(npc_genotypes),
                  selected = "I1061T/I1061T (청소년형)"),
      sliderInput("v_dev", "발달 취약성 v_dev", 0, 3, 0.6, step = 0.1),
      helpText(HTML("<small>잔여 NPC1 활성만으로 발병 연령형이 정해지지 않습니다. ",
                    "저장 표현형이 포화하기 때문에 동일 유전형 형제간에도 ",
                    "10년 차이가 납니다.</small>")),
      checkboxInput("liver", "주산기 간질환 위험 (liver_risk)", FALSE),
      numericInput("years", "추적 기간 (년)", 30, min = 2, max = 70),
      hr(),
      h4("② 치료 (Therapy)"),
      checkboxGroupInput("drugs", "약물", choices = DRUGS, selected = character(0)),
      sliderInput("start", "치료 시작 연령 (년)", 0, 30, 13, step = 0.5),
      sliderInput("dur", "치료 기간 (년)", 0.25, 30, 3, step = 0.25),
      hr(),
      h4("③ 용량 (Dose)"),
      numericInput("mig_mg", "미글루스타트 mg/회 (tid)", 200, min = 0, max = 400, step = 25),
      numericInput("ari_mg", "아리모클로몰 mg/회 (tid)", 124, min = 0, max = 372, step = 31),
      numericInput("nal_g",  "레바세틸류신 g/일", 4, min = 0, max = 6, step = 0.5),
      numericInput("cd_mg",  "아드라베타덱스 mg (IT)", 900, min = 0, max = 1800, step = 100),
      numericInput("cd_wk",  "IT 투여 간격 (주)", 2, min = 1, max = 8, step = 1),
      hr(),
      h4("④ 기전 스위치 (Mechanism)"),
      checkboxInput("nal_sym", "레바세틸류신 대증 성분 ON", TRUE),
      checkboxInput("nal_dm",  "레바세틸류신 질병조절 성분 ON", TRUE),
      helpText(HTML("<small>두 스위치를 따로 끄면 8번 탭의 설계 실험을 ",
                    "손으로 재현할 수 있습니다.</small>"))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. 환자 프로파일",
                 br(), fluidRow(column(6, plotOutput("p_profile", height = 300)),
                                column(6, plotOutput("p_func", height = 300))),
                 hr(), h5("요약 지표"), tableOutput("t_profile"),
                 helpText(HTML("f_NPC1 = 정상 대비 리소좀 NPC1 기능. ",
                               "NPC1과 NPC2는 하나의 flux에 <b>직렬</b>로 작용하므로 ",
                               "기능 분율이 곱해집니다."))),
        tabPanel("2. 약물 PK",
                 br(), plotOutput("p_pk", height = 420),
                 hr(), plotOutput("p_pd_target", height = 300),
                 helpText(HTML("미글루스타트의 UGCG 억제는 <b>내장에서 CNS의 약 2배</b>입니다 ",
                               "(뇌:혈장 분배 0.45). 이는 가정이 아니라 예측이며, ",
                               "미글루스타트의 가장 뚜렷한 효과가 내장·구음연하 쪽인 이유입니다."))),
        tabPanel("3. 리소좀 · 지질",
                 br(), plotOutput("p_lipid", height = 380),
                 hr(), plotOutput("p_lyso", height = 320)),
        tabPanel("4. 소뇌 · 예비능",
                 br(), plotOutput("p_dam", height = 320),
                 hr(), plotOutput("p_pc", height = 340),
                 helpText(HTML("<b>DAM은 스트레스의 적분이며 수선항이 없습니다.</b> ",
                               "발병 지연은 스트레스에 <i>반비례</i>하고, 관문이 열린 뒤의 ",
                               "사멸 속도는 <i>포화</i>합니다 — 그래서 진행 기울기가 발병 ",
                               "연령과 거의 무관해집니다 (PMID 19415691)."))),
        tabPanel("5. 임상 엔드포인트",
                 br(), plotOutput("p_scales", height = 380),
                 hr(), fluidRow(column(6, plotOutput("p_sacc", height = 300)),
                                column(6, plotOutput("p_surv", height = 300)))),
        tabPanel("6. 시나리오 비교",
                 br(),
                 helpText("아래 표의 모든 시나리오는 좌측 패널의 환자에게 적용됩니다."),
                 actionButton("run_scen", "18개 시나리오 실행", class = "btn-primary"),
                 br(), br(), DTOutput("t_scen"),
                 hr(), plotOutput("p_scen", height = 420)),
        tabPanel("7. 바이오마커",
                 br(), plotOutput("p_bio", height = 380),
                 hr(), plotOutput("p_compart", height = 320),
                 helpText(HTML("<b>구획 주장.</b> C-triol 생성은 포화합니다. 환자에서는 ",
                               "이미 자기 천장의 약 82%에 있어 더 이상의 저장을 반영하지 ",
                               "못합니다. 그래서 발표된 triol–NPCCSS5 상관이 ρ = 0.265에 ",
                               "불과합니다 (PMID 33228797) — 잡음 때문이 아니라 ",
                               "<i>다른 구획</i>을 재고 있기 때문입니다."))),
        tabPanel("8. 시험 설계 실험",
                 br(),
                 helpText(HTML("동일한 12개월 효능을 갖도록 맞춘 <b>순수 대증약</b>과 ",
                               "<b>순수 질병조절약</b>을 두 개의 발표된 설계에 각각 통과시킵니다.")),
                 numericInput("dsg_weeks", "설계 A: 관찰 기간 (주)", 12, min = 4, max = 104),
                 actionButton("run_dsg", "설계 실험 실행", class = "btn-primary"),
                 br(), br(), tableOutput("t_design"),
                 hr(), plotOutput("p_design", height = 340),
                 helpText(HTML("12주 창은 D<sub>rev</sub>만 봅니다. 12개월 기저치-대비-변화는 ",
                               "상수 편차와 기울기 변화를 구분하지 못합니다. ",
                               "두 성분을 가르는 유일한 설계는 <b>무작위 중단(withdrawal)</b>입니다."))),
        tabPanel("9. 보정 대조표",
                 br(), h4("모델 vs 발표값"), DTOutput("t_calib"),
                 hr(),
                 h4("맞지 않은 항목 (조정하지 않고 보고)"),
                 htmlOutput("misses"))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  ev_now <- reactive({
    make_events(input$drugs, input$start, input$dur,
                input$mig_mg, input$ari_mg, input$nal_g, input$cd_mg, input$cd_wk)
  })

  sim <- reactive({
    run_sim(input$geno, input$v_dev, as.numeric(input$liver),
            input$years, ev_now(), delta = 7,
            nal_sym_on = as.numeric(input$nal_sym),
            nal_dm_on  = as.numeric(input$nal_dm))
  })

  sim_untr <- reactive({
    run_sim(input$geno, input$v_dev, as.numeric(input$liver),
            input$years, ev(time = 0, amt = 0, cmt = "MIG_GUT"), delta = 7)
  })

  ## ---- tab 1 --------------------------------------------------------------
  output$p_profile <- renderPlot({
    o <- sim()
    longify(o, c("f_NPC1_out", "f_eg", "theta")) %>%
      ggplot(aes(AGE_YR, value, colour = variable)) +
      geom_line(linewidth = 0.9) +
      labs(title = "NPC1 접힘 수율과 배출 능력",
           subtitle = "theta = ER 접힘 수율 · f_eg = NPC1 x NPC2 (직렬)",
           x = "연령 (년)", y = "정상 대비 분율", colour = NULL) +
      theme_npc()
  })

  output$p_func <- renderPlot({
    o <- sim()
    longify(o, c("D_REV", "D_IRR", "FUNC")) %>%
      ggplot(aes(AGE_YR, value, colour = variable)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = c(D_REV = "#e67e22", D_IRR = "#c0392b",
                                     FUNC = "#1a7f37")) +
      labs(title = "기능 = 1 − D_rev − D_irr",
           subtitle = "D_rev는 되돌아오고 D_irr은 되돌아오지 않는다",
           x = "연령 (년)", y = NULL, colour = NULL) +
      theme_npc()
  })

  output$t_profile <- renderTable({
    o <- sim(); u <- sim_untr()
    ages <- c(5, 10, 15, 20)
    ages <- ages[ages <= input$years]
    data.frame(
      `연령` = ages,
      `f_NPC1` = round(sapply(ages, function(a) at_age(o, "f_NPC1_out", a)), 4),
      `CHOL_C (x정상)` = round(sapply(ages, function(a) at_age(o, "CHOL_C_fold", a)), 1),
      `DAM/예비능` = round(sapply(ages, function(a) at_age(o, "DAM", a)) / 5707.9, 2),
      `NPCCSS5` = round(sapply(ages, function(a) at_age(o, "NPCCSS5", a)), 2),
      `SARA` = round(sapply(ages, function(a) at_age(o, "SARA", a)), 2),
      `무치료 NPCCSS5` = round(sapply(ages, function(a) at_age(u, "NPCCSS5", a)), 2),
      check.names = FALSE)
  })

  ## ---- tab 2 --------------------------------------------------------------
  output$p_pk <- renderPlot({
    o <- sim() %>% filter(AGE_YR >= input$start, AGE_YR <= input$start + 0.12)
    if (!nrow(o)) o <- sim()
    longify(o, c("C_MIG", "C_ARI", "C_NAL", "CD_CSF_CONC")) %>%
      ggplot(aes(AGE_YR, value)) +
      geom_line(linewidth = 0.7, colour = "#8e44ad") +
      facet_wrap(~variable, scales = "free_y", ncol = 2) +
      labs(title = "약물 농도 (치료 개시 직후 6주)",
           subtitle = "C_MIG uM · C_ARI ng/mL · C_NAL mg/L · CD_CSF mg/L",
           x = "연령 (년)", y = NULL) +
      theme_npc()
  })

  output$p_pd_target <- renderPlot({
    o <- sim()
    longify(o, c("I_mig_v", "I_mig_c", "HSP70", "E_nal_sym", "E_nal_dm")) %>%
      ggplot(aes(AGE_YR, value, colour = variable)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~variable, scales = "free_y", ncol = 3) +
      labs(title = "표적 수준의 약력학",
           subtitle = "I_mig_v = 내장 UGCG 억제 · I_mig_c = CNS UGCG 억제",
           x = "연령 (년)", y = NULL) +
      theme_npc() + theme(legend.position = "none")
  })

  ## ---- tab 3 --------------------------------------------------------------
  output$p_lipid <- renderPlot({
    o <- sim()
    longify(o, c("CHOL_V_fold", "CHOL_C_fold", "GSL_V", "GSL_C", "SPH")) %>%
      ggplot(aes(AGE_YR, value, colour = variable)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~variable, scales = "free_y", ncol = 3) +
      labs(title = "지질 저장 — 내장은 며칠, 뇌는 수년의 시간상수",
           subtitle = "tau_cns = 25 : 정상상태는 같고 도달 속도만 25배 느리다",
           x = "연령 (년)", y = NULL) +
      theme_npc() + theme(legend.position = "none")
  })

  output$p_lyso <- renderPlot({
    o <- sim()
    longify(o, c("CA_LY", "HYD", "AUTOPH", "MITO", "ROS", "GLUC_METAB")) %>%
      ggplot(aes(AGE_YR, value, colour = variable)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~variable, scales = "free_y", ncol = 3) +
      labs(title = "리소좀 · 대사 기능",
           subtitle = "스핑고신 → 산성 Ca²⁺ 저장 차단 → 융합·가수분해 저하 (PMID 18953351)",
           x = "연령 (년)", y = NULL) +
      theme_npc() + theme(legend.position = "none")
  })

  ## ---- tab 4 --------------------------------------------------------------
  output$p_dam <- renderPlot({
    o <- sim()
    o %>% mutate(reserve_used = DAM / 5707.9) %>%
      select(AGE_YR, reserve_used, gate, stress) %>%
      pivot_longer(-AGE_YR) %>%
      ggplot(aes(AGE_YR, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "손상 적분과 예비능 관문",
           subtitle = "reserve_used = 1 을 지나면 푸르킨예 사멸이 시작된다 (되돌릴 수 없음)",
           x = "연령 (년)", y = NULL) +
      theme_npc() + theme(legend.position = "none")
  })

  output$p_pc <- renderPlot({
    o <- sim()
    longify(o, c("PC", "PC_S", "PC_LOST", "INFL", "SYN", "CBL")) %>%
      ggplot(aes(AGE_YR, value, colour = variable)) +
      geom_line(linewidth = 0.9) +
      labs(title = "푸르킨예 세포 풀 · 신경염증 · 소뇌 부피",
           subtitle = "PC + PC_S + PC_LOST = 1 (보존량)",
           x = "연령 (년)", y = "분율 / 지수", colour = NULL) +
      theme_npc()
  })

  ## ---- tab 5 --------------------------------------------------------------
  output$p_scales <- renderPlot({
    o <- sim(); u <- sim_untr()
    bind_rows(mutate(longify(o, c("SARA","NPCCSS5","NPCCSS4","NPCCSS17")), arm = "치료"),
              mutate(longify(u, c("SARA","NPCCSS5","NPCCSS4","NPCCSS17")), arm = "무치료")) %>%
      ggplot(aes(AGE_YR, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~variable, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = c(`치료` = "#8e44ad", `무치료` = "grey40")) +
      labs(title = "임상 척도", x = "연령 (년)", y = NULL, colour = NULL) +
      theme_npc()
  })

  output$p_sacc <- renderPlot({
    o <- sim(); u <- sim_untr()
    bind_rows(mutate(longify(o, c("SACCADE","SWALLOW","HEARING_dB")), arm = "치료"),
              mutate(longify(u, c("SACCADE","SWALLOW","HEARING_dB")), arm = "무치료")) %>%
      ggplot(aes(AGE_YR, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~variable, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = c(`치료` = "#8e44ad", `무치료` = "grey40")) +
      labs(title = "사카드 속도 · 연하 · 청력역치 이동",
           subtitle = "청력역치는 사이클로덱스트린 이독성 (효능과 같은 기전)",
           x = "연령 (년)", y = NULL, colour = NULL) +
      theme_npc()
  })

  output$p_surv <- renderPlot({
    o <- sim(); u <- sim_untr()
    bind_rows(mutate(select(o, AGE_YR, SURV), arm = "치료"),
              mutate(select(u, AGE_YR, SURV), arm = "무치료")) %>%
      ggplot(aes(AGE_YR, SURV, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
      scale_colour_manual(values = c(`치료` = "#8e44ad", `무치료` = "grey40")) +
      labs(title = "생존 확률", subtitle = "위험은 연하기능의 제곱에 비례 (PMID 23039766)",
           x = "연령 (년)", y = "S(t)", colour = NULL) +
      theme_npc()
  })

  ## ---- tab 6 --------------------------------------------------------------
  scen <- eventReactive(input$run_scen, {
    defs <- list(
      list(id = "S01", lab = "무치료",                    d = character(0), s = 13, u = 3),
      list(id = "S06", lab = "미글루스타트",              d = "miglustat",  s = 13, u = 3),
      list(id = "S07", lab = "레바세틸류신",              d = "levacetylleucine", s = 13, u = 3),
      list(id = "S08", lab = "아리모클로몰",              d = "arimoclomol", s = 13, u = 3),
      list(id = "S09", lab = "아리모클로몰+미글루스타트", d = c("arimoclomol","miglustat"), s = 13, u = 3),
      list(id = "S10", lab = "레바세틸류신+미글루스타트", d = c("levacetylleucine","miglustat"), s = 13, u = 3),
      list(id = "S11", lab = "3제 병용",                  d = c("levacetylleucine","miglustat","arimoclomol"), s = 13, u = 3),
      list(id = "S12", lab = "IT 사이클로덱스트린",       d = "IT adrabetadex", s = 13, u = 3),
      list(id = "S13", lab = "IT CD + 미글루스타트",      d = c("IT adrabetadex","miglustat"), s = 13, u = 3),
      list(id = "S14", lab = "3제, 2세부터",              d = c("levacetylleucine","miglustat","arimoclomol"), s = 2,  u = 18),
      list(id = "S15", lab = "3제, 8세부터",              d = c("levacetylleucine","miglustat","arimoclomol"), s = 8,  u = 12),
      list(id = "S16", lab = "3제, 12세부터",             d = c("levacetylleucine","miglustat","arimoclomol"), s = 12, u = 8)
    )
    withProgress(message = "시나리오 실행 중...", {
      purrr::map_dfr(seq_along(defs), function(i) {
        s <- defs[[i]]; incProgress(1 / length(defs))
        o <- run_sim(input$geno, input$v_dev, as.numeric(input$liver), 20,
                     make_events(s$d, s$s, s$u, input$mig_mg, input$ari_mg,
                                 input$nal_g, input$cd_mg, input$cd_wk), delta = 14)
        mutate(o, scenario = s$id, label = s$lab)
      })
    })
  })

  output$t_scen <- renderDT({
    scen() %>% group_by(scenario, label) %>%
      summarise(`NPCCSS5@20` = round(at_age(cur_data_all(), "NPCCSS5", 20), 2),
                `SARA@20` = round(at_age(cur_data_all(), "SARA", 20), 2),
                `D_rev@20` = round(at_age(cur_data_all(), "D_REV", 20), 3),
                `D_irr@20` = round(at_age(cur_data_all(), "D_IRR", 20), 3),
                `PC_lost@20` = round(at_age(cur_data_all(), "PC_LOST", 20), 3),
                `청력 dB@20` = round(at_age(cur_data_all(), "HEARING_dB", 20), 1),
                .groups = "drop") %>%
      datatable(options = list(pageLength = 12, dom = "t"), rownames = FALSE)
  })

  output$p_scen <- renderPlot({
    ggplot(scen(), aes(AGE_YR, NPCCSS5, colour = label)) +
      geom_line(linewidth = 0.8) +
      labs(title = "시나리오별 5영역 NPCCSS",
           subtitle = "조기 개시가 사는 것은 D_rev이고, 이미 잃은 것은 D_irr이다",
           x = "연령 (년)", y = "5-domain NPCCSS", colour = NULL) +
      theme_npc()
  })

  ## ---- tab 8 --------------------------------------------------------------
  design <- eventReactive(input$run_dsg, {
    wk <- input$dsg_weeks
    horizon <- 13 * YR + wk * 7
    run1 <- function(sym, dm, drug_on) {
      e <- if (drug_on) make_events("levacetylleucine", 13, wk * 7 / YR,
                                    nal_g = input$nal_g)
           else ev(time = 0, amt = 0, cmt = "MIG_GUT")
      run_sim(input$geno, input$v_dev, as.numeric(input$liver),
              horizon / YR, e, delta = 3.5,
              nal_sym_on = sym, nal_dm_on = dm)
    }
    pbo <- run1(1, 1, FALSE); sym <- run1(1, 0, TRUE); dmo <- run1(0, 1, TRUE)
    tab <- tibble(
      arm = c("위약", "순수 대증", "순수 질병조절"),
      dSARA = c(at_age(pbo,"SARA",horizon/YR) - at_age(pbo,"SARA",13),
                at_age(sym,"SARA",horizon/YR) - at_age(sym,"SARA",13),
                at_age(dmo,"SARA",horizon/YR) - at_age(dmo,"SARA",13)),
      dN5   = c(at_age(pbo,"NPCCSS5",horizon/YR) - at_age(pbo,"NPCCSS5",13),
                at_age(sym,"NPCCSS5",horizon/YR) - at_age(sym,"NPCCSS5",13),
                at_age(dmo,"NPCCSS5",horizon/YR) - at_age(dmo,"NPCCSS5",13))
    ) %>% mutate(`SARA 위약대비` = dSARA - dSARA[1],
                 `NPCCSS5 위약대비` = dN5 - dN5[1])
    list(tab = tab, traj = bind_rows(mutate(pbo, arm = "위약"),
                                     mutate(sym, arm = "순수 대증"),
                                     mutate(dmo, arm = "순수 질병조절")))
  })

  output$t_design <- renderTable({ design()$tab }, digits = 3)

  output$p_design <- renderPlot({
    design()$traj %>%
      ggplot(aes(AGE_YR, SARA, colour = arm)) +
      geom_line(linewidth = 0.9) +
      labs(title = "같은 12개월 효능, 전혀 다른 단기 SARA 궤적",
           subtitle = "짧은 창은 대증 성분만 본다",
           x = "연령 (년)", y = "SARA", colour = NULL) +
      theme_npc()
  })

  ## ---- tab 9 --------------------------------------------------------------
  output$t_calib <- renderDT({
    oj <- run_sim("I1061T/I1061T (청소년형)", 0.6, 0,
                  35, ev(time = 0, amt = 0, cmt = "MIG_GUT"), delta = 7)
    ow <- run_sim("WT (정상)", 0.6, 0, 45,
                  ev(time = 0, amt = 0, cmt = "MIG_GUT"), delta = 30)
    tibble(
      `목표` = c("T1 혈장 triol, 환자 (ng/mL)", "T2 혈장 triol, 대조 (ng/mL)",
                 "T3 5영역 NPCCSS 기울기 (점/년)", "T4 17영역 기울기 (점/년)",
                 "T11 시험 진입 SARA"),
      `역할` = c("보정", "보정", "보정", "검증", "보정"),
      `발표값` = c(88.31, 5.97, 1.50, 2.80, 15.91),
      `모델` = round(c(at_age(oj, "TRIOL", 10), at_age(ow, "TRIOL", 10),
                       at_age(oj,"NPCCSS5",13) - at_age(oj,"NPCCSS5",12),
                       at_age(oj,"NPCCSS17",13) - at_age(oj,"NPCCSS17",12),
                       at_age(oj, "SARA", 13)), 3)
    ) %>% mutate(`상대오차 %` = round(100 * (`모델` - `발표값`) / `발표값`, 1)) %>%
      datatable(options = list(dom = "t"), rownames = FALSE)
  })

  output$misses <- renderUI({
    HTML(paste0(
      "<ul>",
      "<li><b>아리모클로몰 효과 크기.</b> 발표된 12개월 5영역 차이는 −1.40 ",
      "(95% CI −2.76 ~ −0.03)입니다. 모델은 접힘 수율을 <i>완전히</i> 복구시켜도 ",
      "진행의 약 48%만 제거하며, 13세 진입 환자에서 −0.7점 수준에 머무릅니다. ",
      "모델값은 발표된 신뢰구간 <i>안</i>에 있지만 점추정치보다 훨씬 작습니다. ",
      "모델의 해석: 예비능 적분은 되돌릴 수 없으므로 <b>예비능을 이미 다 쓴 뒤에 ",
      "원발 결함을 고쳐도 진행 속도를 정상으로 되돌릴 수 없습니다.</b> ",
      "더 이른 시점에 진입한 환자에서는 효과가 더 큽니다 (4번 탭 참조).</li>",
      "<li><b>레바세틸류신 장기 연장.</b> 발표된 기저치 대비 SARA 변화는 12개월 ",
      "−1.88, 18개월 −1.64입니다. 모델은 12주 차이(−1.28)를 정확히 재현하지만, ",
      "그 크기의 대증 효과라면 기저치 대비 변화가 12개월에는 이미 0을 넘어야 ",
      "합니다. 발표된 값들이 함의하는 무치료 SARA 진행은 연 0.22점인데, ",
      "모델이 NPCCSS 목표로부터 <i>도출</i>하는 값은 연 1.5점입니다 — ",
      "약 7배 불일치입니다. NPC의 무치료 SARA 기울기는 발표된 값이 없어 ",
      "이 불일치는 미해결로 남겨 둡니다.</li>",
      "<li><b>후기영아형과 청소년형의 분리.</b> 저장 표현형이 포화하기 때문에 ",
      "잔여 NPC1 활성 2배 차이가 정상상태 부하를 거의 바꾸지 못합니다. ",
      "두 연령형은 발달 취약성 v_dev로 <i>따로</i> 지정되며, 유전형에서 ",
      "도출되지 않습니다.</li>",
      "</ul>"))
  })
}

shinyApp(ui, server)
