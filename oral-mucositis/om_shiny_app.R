# =============================================================================
#  Oral Mucositis (OM) QSP — Shiny dashboard
#  구강점막염 정량적 시스템 약리학 대시보드
# =============================================================================
#
#  The app is organised around the model's ONE claim: onset belongs to the
#  INSULT term and duration to the REGENERATION term, and the two belong to
#  different drugs.  Tabs 8-11 are not plots of the model, they are the
#  EXPERIMENTS that test that claim — the same experiments that
#  om_analysis.py runs offline.
#
#  Run:
#    install.packages(c("shiny","mrgsolve","dplyr","tidyr","ggplot2",
#                       "DT","tibble"))
#    shiny::runApp("om_shiny_app.R")
#
#  NOTE: no R toolchain was available where this file was written, so it has
#  not been executed.  It mirrors om_python_reference.py, which was.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(DT)

MODEL_FILE <- "om_mrgsolve_model.R"
mod <- mread_cache("om", MODEL_FILE)

# -----------------------------------------------------------------------------
# schedule builders (mirrors the block at the bottom of the mrgsolve model)
# -----------------------------------------------------------------------------
FX_DUR <- 10 / 1440          # a 2 Gy fraction is delivered over 10 minutes

blank_cov <- function(times) {
  tibble(ID = 1, time = times, evid = 0, cmt = 0, amt = 0, rate = 0)
}

rt_records <- function(nfx = 35, dpf = 2, per_week = 5, t0 = 0,
                       gap = NULL) {
  day <- t0; given <- 0; starts <- c()
  while (given < nfx && (day - t0) < 250) {
    in_gap <- !is.null(gap) && (day - t0) >= gap[1] && (day - t0) < gap[2]
    if ((round(day - t0) %% 7) < per_week && !in_gap) {
      starts <- c(starts, day); given <- given + 1
    }
    day <- day + 1
  }
  bind_rows(lapply(starts, function(s)
    tibble(ID = 1, time = c(s, s + FX_DUR), evid = 0, cmt = 0, amt = 0,
           rate = 0, RRATE = c(dpf / FX_DUR, 0), DPF = dpf)))
}

ev_inf <- function(times, amt, cmt, dur_d) {
  as_tibble(as.data.frame(
    ev(ID = 1, time = times, amt = amt, cmt = cmt, rate = amt / dur_d)))
}

build_schedule <- function(inp) {
  recs <- list()

  if (inp$regimen == "HDM") {
    amt <- inp$mel_dose * inp$bsa * 1000
    recs[[length(recs) + 1]] <- ev_inf(0, amt, "A_mel_c", 0.5 / 24)

  } else if (inp$regimen == "TBI") {
    dur <- FX_DUR
    st <- c(0, 0.35, 1, 1.35, 2, 2.35, 3, 3.35)
    recs[[length(recs) + 1]] <- bind_rows(lapply(st, function(s)
      tibble(ID = 1, time = c(s, s + dur), evid = 0, cmt = 0, amt = 0,
             rate = 0, RRATE = c(1.5 / dur, 0), DPF = 1.5)))
    amt <- inp$cy_equiv * inp$bsa * 1000
    recs[[length(recs) + 1]] <- ev_inf(c(4, 6), amt / 2, "A_mel_c", 0.5 / 24)

  } else if (inp$regimen == "chemoRT") {
    gap <- if (inp$gap_len > 0) c(inp$gap_start, inp$gap_start + inp$gap_len)
           else NULL
    recs[[length(recs) + 1]] <- rt_records(inp$nfx, inp$dpf, inp$fx_week,
                                           gap = gap)
    if (inp$cis_dose > 0) {
      amt <- inp$cis_dose * inp$bsa * 1000
      recs[[length(recs) + 1]] <- ev_inf(c(0, 21, 42), amt, "A_cis_c", 2 / 24)
    }

  } else if (inp$regimen == "5FU_bolus") {
    amt <- inp$fu_dose * inp$bsa * 1000
    recs[[length(recs) + 1]] <- ev_inf(0:4, amt, "A_5fu_c", 5 / 1440)

  } else if (inp$regimen == "5FU_CI") {
    amt <- inp$fu_total * inp$bsa * 1000
    recs[[length(recs) + 1]] <- ev_inf(0, amt, "A_5fu_c", inp$fu_days)
  }

  # ---- interventions --------------------------------------------------------
  if (inp$cryo_h > 0) {
    recs[[length(recs) + 1]] <- tibble(ID = 1, time = c(0, inp$cryo_h / 24),
                                       evid = 0, cmt = 0, amt = 0, rate = 0,
                                       CRYO = c(1, 0))
  }
  if (inp$palifermin) {
    amt <- 60 * inp$wt * 1000
    d0 <- inp$pal_gap
    days <- c(-2, -1, 0) - d0 + inp$pal_anchor
    days2 <- inp$pal_post + c(0, 1, 2)
    recs[[length(recs) + 1]] <-
      as_tibble(as.data.frame(ev(ID = 1, time = c(days, days2), amt = amt,
                                 cmt = "A_pal_c")))
  }
  if (inp$pbm) {
    d <- seq(0, inp$pbm_days, by = 1)
    recs[[length(recs) + 1]] <- bind_rows(lapply(d, function(x)
      tibble(ID = 1, time = c(x, x + 5 / 1440), evid = 0, cmt = 0, amt = 0,
             rate = 0, PBMR = c(inp$pbm_fluence / (5 / 1440), 0))))
  }
  if (inp$benzydamine) {
    recs[[length(recs) + 1]] <- tibble(ID = 1, time = c(0, inp$t_end),
                                       evid = 0, cmt = 0, amt = 0, rate = 0,
                                       BZDL = c(1, 0))
  }
  if (inp$glutamine) {
    recs[[length(recs) + 1]] <- tibble(ID = 1, time = c(0, inp$t_end),
                                       evid = 0, cmt = 0, amt = 0, rate = 0,
                                       GLNL = c(1, 0))
  }

  d <- bind_rows(recs)
  for (nm in c("RRATE", "CRYO", "PBMR", "BZDL", "GLNL")) {
    if (!nm %in% names(d)) d[[nm]] <- 0
    d[[nm]][is.na(d[[nm]])] <- 0
  }
  if (!"DPF" %in% names(d)) d$DPF <- 2
  d$DPF[is.na(d$DPF)] <- 2
  for (nm in c("evid", "cmt", "amt", "rate")) {
    if (!nm %in% names(d)) d[[nm]] <- 0
    d[[nm]][is.na(d[[nm]])] <- 0
  }
  # LOCF the covariates so they persist between records
  d <- d %>% arrange(time)
  d
}

run_sim <- function(inp, pset = list(), sched = NULL, t_end = NULL) {
  s <- if (is.null(sched)) build_schedule(inp) else sched
  te <- if (is.null(t_end)) inp$t_end else t_end
  m <- mod
  pset$sens <- (pset$sens %||% 1) * inp$sens
  if (length(pset)) m <- param(m, pset)
  m %>% data_set(s) %>%
    mrgsim(end = te, delta = 0.05, recsort = 3) %>% as_tibble()
}

`%||%` <- function(a, b) if (is.null(a)) b else a

summarise_run <- function(out, sev = 3) {
  dt <- 0.05
  sm <- out$WHO >= sev
  tibble(
    onset_sev = if (any(sm)) min(out$time[sm]) else NA_real_,
    dur_sev   = sum(sm) * dt,
    peak_area = max(out$ULC),
    peak_who  = max(out$WHO),
    area_auc  = sum(out$ULC) * dt,
    pain_auc  = sum(out$VAS) * dt,
    opioid_pk = max(out$MEDMG),
    anc_nadir = min(out$ANC))
}

# =============================================================================
# UI
# =============================================================================
ui <- fluidPage(
  titlePanel(HTML(paste0(
    "<b>구강점막염 QSP 대시보드</b> — Oral Mucositis QSP<br>",
    "<span style='font-size:14px;color:#555'>",
    "하나의 비대칭: 모든 세포독성은 상피 벨트의 <b>증식 구획</b>을 때리고 ",
    "<b>분화된 장벽</b>은 건드리지 않는다 → <b>발병 시점은 손상 항</b>이, ",
    "<b>지속 기간은 재생 항</b>이 결정한다",
    "</span>"))),
  tags$style(HTML("
    .well { background:#fafbfc; }
    .nav-tabs > li > a { padding:6px 10px; font-size:13px; }
    .clockbox { border-left:4px solid #b9770e; padding:6px 12px;
                background:#fef9e7; margin-bottom:8px; font-size:13px; }
  ")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① 요법 (Regimen)"),
      selectInput("regimen", NULL,
                  c("고용량 멜팔란 (HDM, 자가조혈모세포이식)" = "HDM",
                    "TBI + VP16/Cy 전처치" = "TBI",
                    "두경부 항암방사선 (70 Gy + 시스플라틴)" = "chemoRT",
                    "5-FU 볼루스 d1-5" = "5FU_bolus",
                    "5-FU 96시간 지속주입" = "5FU_CI"),
                  selected = "HDM"),
      conditionalPanel("input.regimen == 'HDM'",
        sliderInput("mel_dose", "멜팔란 (mg/m²)", 100, 200, 200, 10)),
      conditionalPanel("input.regimen == 'TBI'",
        sliderInput("cy_equiv", "알킬화제 등가용량 (mg/m² 멜팔란 상당)",
                    20, 300, 120, 10)),
      conditionalPanel("input.regimen == 'chemoRT'",
        sliderInput("nfx", "분할 횟수", 20, 70, 35, 1),
        sliderInput("dpf", "분할당 선량 (Gy)", 1.0, 2.5, 2.0, 0.1),
        sliderInput("fx_week", "주당 분할", 5, 10, 5, 1),
        sliderInput("cis_dose", "시스플라틴 (mg/m², 0 = RT 단독)",
                    0, 100, 100, 10),
        sliderInput("gap_start", "치료 중단 시작일", 0, 45, 28, 1),
        sliderInput("gap_len", "치료 중단 기간 (일)", 0, 21, 0, 1)),
      conditionalPanel("input.regimen == '5FU_bolus'",
        sliderInput("fu_dose", "5-FU 볼루스 (mg/m²/일)", 200, 600, 425, 25)),
      conditionalPanel("input.regimen == '5FU_CI'",
        sliderInput("fu_total", "5-FU 총량 (mg/m²)", 1000, 6000, 4000, 250),
        sliderInput("fu_days", "주입 기간 (일)", 1, 7, 4, 1)),

      hr(),
      h4("② 손상 축 개입 (INSULT arm)"),
      div(class = "clockbox",
          "시계 1 — 손상: 발병 시점·발생률을 움직인다"),
      sliderInput("cryo_h", "구강 한랭요법 (시간)", 0, 8, 0, 0.5),

      hr(),
      h4("③ 재생 축 개입 (REGENERATION arm)"),
      div(class = "clockbox",
          "시계 2 — 재생: 지속 기간을 움직인다"),
      checkboxInput("palifermin", "팔리페르민 60 µg/kg/일 × 3+3", FALSE),
      conditionalPanel("input.palifermin",
        sliderInput("pal_gap",
                    "마지막 전처치 전 투여와 세포독성 사이 간격 (일)",
                    -1, 3, 1, 0.5),
        sliderInput("pal_anchor", "세포독성 투여 시점 (일)", 0, 8, 0, 1),
        sliderInput("pal_post", "이식 후 재투여 시작일", 1, 14, 3, 1)),
      checkboxInput("pbm", "광생체조절 (매일)", FALSE),
      conditionalPanel("input.pbm",
        sliderInput("pbm_fluence", "조사량 (J/cm²)", 1, 12, 6, 1),
        sliderInput("pbm_days", "시행 기간 (일)", 5, 40, 20, 1)),
      checkboxInput("glutamine", "경구 글루타민", FALSE),

      hr(),
      h4("④ 항염 · 환자 (Anti-inflammatory / host)"),
      checkboxInput("benzydamine", "벤지다민 가글", FALSE),
      sliderInput("sens", "개인 감수성 배수 (sens)", 0.5, 2.0, 1.0, 0.05),
      sliderInput("bsa", "체표면적 (m²)", 1.4, 2.2, 1.8, 0.05),
      sliderInput("wt", "체중 (kg)", 45, 110, 75, 1),
      sliderInput("t_end", "시뮬레이션 기간 (일)", 20, 130, 45, 5),
      helpText(HTML("모든 파라미터의 출처와 보정 내역은 ",
                    "<code>calib.log</code> 및 ",
                    "<code>om_reference_output.txt</code> 참조."))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",

        tabPanel("1. 환자·요법",
          h4("현재 설정된 투여/조사 스케줄"),
          p(HTML("방사선 분할은 <b>10분 창</b>으로 전달됩니다. 적분기가 ",
                 "이 창을 건너뛰면 실제 전달 선량이 처방이 아니라 ",
                 "적분기의 스텝 배치에 좌우되므로, 모든 분할 경계는 ",
                 "데이터셋 레코드로 명시됩니다.")),
          DTOutput("tbl_sched"),
          h4("요약 지표"), DTOutput("tbl_summary")),

        tabPanel("2. 약동학 (PK)",
          h4("혈장 농도"), plotOutput("p_pk", height = 300),
          h4("점막 구획 농도 — 한랭요법이 작용하는 유일한 지점"),
          p(HTML("점막은 관류제한 구획으로, 평형 시간상수가 약 2분입니다. ",
                 "따라서 <b>혈류를 줄이는 것만으로는</b> 노출을 크게 ",
                 "낮출 수 없습니다 — 자세한 분해는 탭 9.")),
          plotOutput("p_muc", height = 300)),

        tabPanel("3. 상피 벨트",
          h4("S → P1 → P2 → P3 → D — 모델의 척추"),
          p(HTML("세포독성은 <b>증식 구획(S, P)</b>만 때립니다. ",
                 "분화된 장벽 D에는 세포독성 항이 없습니다. ",
                 "그래서 잠복기가 존재하며, 그 길이는 신호전달의 ",
                 "수수께끼가 아니라 <b>1/k_shed</b>입니다.")),
          plotOutput("p_belt", height = 380),
          plotOutput("p_freg", height = 240)),

        tabPanel("4. 궤양 · 중증도",
          plotOutput("p_ulcer", height = 340),
          plotOutput("p_scales", height = 300),
          p(HTML("<b>WHO 등급은 순서형이며 4등급에서 포화</b>합니다. ",
                 "OMAS/면적은 포화하지 않습니다 — 탭 12 참조."))),

        tabPanel("5. 염증 캐스케이드",
          h4("Sonis 2·3기 — 하나의 빠른 신호 루프"),
          p(HTML("다섯 '단계'는 순차적 사건이 아니라, 느린 조직 변수(벨트) ",
                 "위에 얹힌 <b>하나의 빠른 신호 루프</b>입니다. ",
                 "TNF → NF-κB 양성 되먹임의 루프 이득은 반드시 1보다 ",
                 "작아야 하며, 그렇지 않으면 모델이 발산합니다.")),
          plotOutput("p_infl", height = 340),
          h4("궤양면의 미생물 집락과 세라마이드"),
          plotOutput("p_micro", height = 260)),

        tabPanel("6. 통증 · 진통",
          plotOutput("p_pain", height = 320),
          h4("마약성 진통제 요구량 (VAS ≤ 4 목표 적정)"),
          plotOutput("p_opioid", height = 260)),

        tabPanel("7. 골수억제 · 감염",
          plotOutput("p_anc", height = 300),
          h4("궤양은 침입 관문 — 위험도 ∝ 면적 × 세균부하 / ANC"),
          p(HTML("점막 곡선과 ANC 곡선은 <b>같은</b> 멜팔란 노출로 ",
                 "구동되지만 시간상수가 완전히 다릅니다. 따라서 위험도는 ",
                 "서로 다른 시점에 정점을 찍는 두 곡선의 곱입니다.")),
          plotOutput("p_inf", height = 260)),

        tabPanel("8. 두 개의 시계",
          h4("모델의 중심 주장을 2×2로 검정한다"),
          p(HTML("손상 항(<code>sens</code>)과 재생 항(<code>lamS</code>)을 ",
                 "각각 ±40% 흔들고, <b>발병 시점</b>과 <b>지속 기간</b>이 ",
                 "어떻게 움직이는지 본다. 주장이 옳다면 각 항은 자기 ",
                 "평가변수만 크게 움직여야 한다.")),
          actionButton("run_clocks", "실행", class = "btn-primary"),
          br(), br(),
          DTOutput("tbl_clocks"),
          plotOutput("p_clocks", height = 320)),

        tabPanel("9. 한랭요법 기준",
          h4("얼음은 입 안에 있는 동안에만 작용한다"),
          p(HTML("따라서 이득 = (1 − f_cryo) × (얼음 창 안에 들어오는 손상 ",
                 "AUC 비율)이고, 반감기 t½ 약물을 투여 직후 T 동안 ",
                 "냉각하면 그 비율은 정확히 <b>1 − 2^(−T/t½)</b>입니다. ",
                 "가이드라인이 보고하는 '어떤 요법에 듣는가'의 패턴이 ",
                 "여기서 <b>유도</b>됩니다.")),
          actionButton("run_cryo", "실행", class = "btn-primary"),
          br(), br(),
          plotOutput("p_cryo", height = 330),
          DTOutput("tbl_cryo")),

        tabPanel("10. 팔리페르민 스케줄링",
          h4("재생을 켜면 표적도 커진다"),
          p(HTML("KGF는 증식 속도를 올립니다 — 그것이 목적입니다. ",
                 "그러나 <b>커진 증식 풀은 세포주기 특이적 세포독성에 ",
                 "더 큰 표적</b>입니다. 동시에 투여하면 이득이 줄고 ",
                 "역전될 수 있으며, 이것이 라벨의 24시간 분리 요구를 ",
                 "낳습니다. 아래 스윕은 그 요구를 <b>재현</b>합니다.")),
          actionButton("run_pal", "실행", class = "btn-primary"),
          br(), br(),
          plotOutput("p_pal", height = 330),
          DTOutput("tbl_pal")),

        tabPanel("11. 분할조사 · 치료중단",
          h4("점막의 회복과 종양의 재증식은 같은 시간을 두고 경쟁한다"),
          actionButton("run_frac", "실행", class = "btn-primary"),
          br(), br(),
          DTOutput("tbl_frac"),
          plotOutput("p_frac", height = 300)),

        tabPanel("12. 척도 포화 · 검정 민감도",
          h4("WHO 척도는 질환이 가장 나쁜 곳에서 검정력을 잃는다"),
          p(HTML("동일한 개입을 WHO 등급과 궤양 면적으로 각각 측정하고, ",
                 "기저 중증도를 <code>sens</code>로만 바꿔 봅니다. ",
                 "WHO 열이 0으로 수축하는데 면적 열은 그렇지 않다면, ",
                 "바뀐 것은 약이 아니라 <b>계측기</b>입니다.")),
          actionButton("run_sat", "실행", class = "btn-primary"),
          br(), br(),
          DTOutput("tbl_sat"),
          plotOutput("p_sat", height = 300)),

        tabPanel("13. 시나리오 비교",
          h4("17개 시나리오 일괄 실행"),
          actionButton("run_scen", "실행", class = "btn-primary"),
          br(), br(),
          DTOutput("tbl_scen"),
          plotOutput("p_scen", height = 340)),

        tabPanel("14. 가상 집단",
          h4("개인차는 구조가 허용하는 세 곳에만 둔다"),
          p(HTML("세포독성 감수성 · 재생 능력 · 장벽 역치. ",
                 "나머지는 공유합니다.")),
          sliderInput("npop", "집단 크기", 20, 400, 100, 20),
          actionButton("run_pop", "실행", class = "btn-primary"),
          br(), br(),
          DTOutput("tbl_pop"),
          plotOutput("p_pop", height = 320)),

        tabPanel("15. 문서",
          h4("이 모델이 쓴 숫자"),
          verbatimTextOutput("txt_prov"))
      )
    )
  )
)

# =============================================================================
# server
# =============================================================================
server <- function(input, output, session) {

  sim <- reactive({
    req(input$regimen)
    run_sim(input)
  })

  output$tbl_sched <- renderDT({
    datatable(build_schedule(input), options = list(pageLength = 12,
                                                    scrollX = TRUE))
  })

  output$tbl_summary <- renderDT({
    datatable(summarise_run(sim()) %>%
                mutate(across(everything(), ~ round(.x, 3))),
              options = list(dom = "t"))
  })

  # ---- 2. PK ---------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim() %>% select(time, CMEL, C5FU, CCIS, CPAL, CMOR) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.7) +
      facet_wrap(~ name, scales = "free_y", nrow = 1) +
      labs(x = "일 (day)", y = "농도") + theme_bw()
  })

  output$p_muc <- renderPlot({
    d <- sim() %>% select(time, Cm_mel, Cm_5fu, Cm_cis, Cm_mtx) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(colour = "#148f77",
                                            linewidth = 0.8) +
      facet_wrap(~ name, scales = "free_y", nrow = 1) +
      labs(x = "일 (day)", y = "점막 농도") + theme_bw()
  })

  # ---- 3. belt -------------------------------------------------------------
  output$p_belt <- renderPlot({
    d <- sim() %>%
      transmute(time, `S 클론원성 기저세포` = S, `Sd 치명 손상 세포` = Sd,
                `P1+P2+P3 증폭구획` = P1 + P2 + P3,
                `D 분화 장벽` = BARRIER) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 0.34, linetype = "dashed",
                 colour = "#c0392b") +
      annotate("text", x = 0, y = 0.36, hjust = 0, size = 3,
               colour = "#c0392b", label = "D_crit — 궤양 역치") +
      labs(x = "일 (day)", y = "정규화 값", colour = NULL) +
      theme_bw() + theme(legend.position = "top")
  })

  output$p_freg <- renderPlot({
    d <- sim() %>% transmute(time, KGFe, PBMe, GLNe) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(x = "일 (day)", y = "재생 조절 풀", colour = NULL) + theme_bw()
  })

  # ---- 4. ulcer ------------------------------------------------------------
  output$p_ulcer <- renderPlot({
    ggplot(sim(), aes(time, ULC)) +
      geom_area(fill = "#e59866", alpha = 0.6) +
      geom_line(linewidth = 0.8) +
      labs(x = "일 (day)", y = "궤양 면적 분율 A_ulc") + theme_bw()
  })

  output$p_scales <- renderPlot({
    d <- sim() %>% transmute(time, `WHO 등급 (0-4)` = WHO,
                             `OMAS (0-5)` = OMAS,
                             `궤양 면적 ×5` = 5 * ULC) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_step(data = ~ filter(.x, name == "WHO 등급 (0-4)"),
                linewidth = 0.9) +
      geom_line(data = ~ filter(.x, name != "WHO 등급 (0-4)"),
                linewidth = 0.8) +
      labs(x = "일 (day)", y = NULL, colour = NULL) +
      theme_bw() + theme(legend.position = "top")
  })

  # ---- 5. inflammation -----------------------------------------------------
  output$p_infl <- renderPlot({
    d <- sim() %>% select(time, ROS, NFkB, TNF, IL1b, IL6) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      facet_wrap(~ name, scales = "free_y") +
      labs(x = "일 (day)", y = NULL) + theme_bw() +
      theme(legend.position = "none")
  })

  output$p_micro <- renderPlot({
    d <- sim() %>% select(time, CER, MB) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(x = "일 (day)", y = NULL, colour = NULL) + theme_bw()
  })

  # ---- 6. pain -------------------------------------------------------------
  output$p_pain <- renderPlot({
    d <- sim() %>% transmute(time, `진통 전 VAS` = VAS_RAW,
                             `진통 후 VAS` = VAS) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 4, linetype = "dashed") +
      labs(x = "일 (day)", y = "VAS (0-10)", colour = NULL) +
      theme_bw() + theme(legend.position = "top")
  })

  output$p_opioid <- renderPlot({
    ggplot(sim(), aes(time, MEDMG)) + geom_area(fill = "#f0b27a",
                                                alpha = 0.7) +
      labs(x = "일 (day)", y = "정맥 모르핀 등가 (mg/일)") + theme_bw()
  })

  # ---- 7. marrow -----------------------------------------------------------
  output$p_anc <- renderPlot({
    ggplot(sim(), aes(time, ANC)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 0.5, linetype = "dashed",
                 colour = "#c0392b") +
      labs(x = "일 (day)", y = "ANC (×10⁹/L)") + theme_bw()
  })

  output$p_inf <- renderPlot({
    d <- sim() %>% transmute(time, `궤양 면적` = ULC,
                             `세균 부하` = MB,
                             `누적 균혈증 위험` = cInf) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(x = "일 (day)", y = NULL, colour = NULL) +
      theme_bw() + theme(legend.position = "top")
  })

  # ---- 8. two clocks -------------------------------------------------------
  clocks <- eventReactive(input$run_clocks, {
    base <- summarise_run(run_sim(input))
    rows <- list(mutate(base, perturbation = "기준 (baseline)"))
    for (r in list(list("손상 ×0.6 (INSULT)", "sens", 0.6),
                   list("손상 ×1.4 (INSULT)", "sens", 1.4),
                   list("재생 ×0.6 (REGEN)", "lamS", 0.6),
                   list("재생 ×1.4 (REGEN)", "lamS", 1.4))) {
      pv <- as.numeric(param(mod)[[r[[2]]]]) * r[[3]]
      p <- setNames(list(pv), r[[2]])
      rows[[length(rows) + 1]] <- summarise_run(run_sim(input, pset = p)) %>%
        mutate(perturbation = r[[1]])
    }
    bind_rows(rows) %>% relocate(perturbation)
  })
  output$tbl_clocks <- renderDT(datatable(
    clocks() %>% mutate(across(where(is.numeric), ~ round(.x, 3))),
    options = list(dom = "t")))
  output$p_clocks <- renderPlot({
    d <- clocks() %>% select(perturbation, onset_sev, dur_sev) %>%
      pivot_longer(-perturbation)
    ggplot(d, aes(perturbation, value, fill = name)) +
      geom_col(position = "dodge") + coord_flip() +
      labs(x = NULL, y = "일 (day)", fill = NULL) + theme_bw()
  })

  # ---- 9. cryotherapy ------------------------------------------------------
  cryo <- eventReactive(input$run_cryo, {
    th <- 1.2 / 24
    hrs <- c(0, 0.25, 0.5, 1, 2, 4, 6, 8)
    bind_rows(lapply(hrs, function(h) {
      i2 <- modifyList(reactiveValuesToList(input), list(cryo_h = h))
      o <- run_sim(i2)
      s <- summarise_run(o)
      tibble(ice_h = h, dur_sev = s$dur_sev, peak_area = s$peak_area,
             auc_muc = max(o$cAUCmuc),
             analytic_cover = 1 - 2 ^ (-(h / 24) / th))
    })) %>% mutate(auc_ratio = auc_muc / auc_muc[1])
  })
  output$p_cryo <- renderPlot({
    d <- cryo()
    ggplot(d, aes(ice_h)) +
      geom_line(aes(y = 1 - auc_ratio, colour = "모델 (AUC 감소분)"),
                linewidth = 1) +
      geom_line(aes(y = analytic_cover * (1 - 0.22),
                    colour = "해석해 (1−2^(−T/t½))×(1−f)"),
                linetype = "dashed", linewidth = 1) +
      labs(x = "얼음 유지 시간 (h)", y = "손상 감소 분율",
           colour = NULL) + theme_bw() + theme(legend.position = "top")
  })
  output$tbl_cryo <- renderDT(datatable(
    cryo() %>% mutate(across(where(is.numeric), ~ round(.x, 4))),
    options = list(dom = "t")))

  # ---- 10. palifermin ------------------------------------------------------
  pal <- eventReactive(input$run_pal, {
    gaps <- c(-1, -0.5, 0, 0.5, 1, 2, 3)
    ref <- summarise_run(run_sim(modifyList(reactiveValuesToList(input),
                                            list(palifermin = FALSE))))
    bind_rows(lapply(gaps, function(g) {
      i2 <- modifyList(reactiveValuesToList(input),
                       list(palifermin = TRUE, pal_gap = g))
      s <- summarise_run(run_sim(i2))
      tibble(gap_d = g, dur_sev = s$dur_sev, peak_area = s$peak_area,
             delta_pct = 100 * (s$dur_sev - ref$dur_sev) / ref$dur_sev)
    })) %>% mutate(no_pal_dur = ref$dur_sev)
  })
  output$p_pal <- renderPlot({
    ggplot(pal(), aes(gap_d, dur_sev)) +
      geom_line(linewidth = 1, colour = "#117864") +
      geom_point(size = 2) +
      geom_hline(aes(yintercept = no_pal_dur), linetype = "dashed",
                 colour = "#c0392b") +
      annotate("text", x = -1, y = Inf, vjust = 1.4, hjust = 0, size = 3.4,
               colour = "#c0392b", label = "팔리페르민 없음") +
      labs(x = "마지막 사전투여 ↔ 세포독성 간격 (일)",
           y = "중증 점막염 지속 기간 (일)") + theme_bw()
  })
  output$tbl_pal <- renderDT(datatable(
    pal() %>% mutate(across(where(is.numeric), ~ round(.x, 3))),
    options = list(dom = "t")))

  # ---- 11. fractionation ---------------------------------------------------
  frac <- eventReactive(input$run_frac, {
    arms <- list(
      list("70 Gy/35 fx/7주 (표준)", 35, 2.0, 5, 0),
      list("81.6 Gy/68 fx b.i.d.", 68, 1.2, 10, 0),
      list("70 Gy/35 fx/6주 (가속)", 35, 2.0, 6, 0),
      list("70 Gy + 7일 중단", 35, 2.0, 5, 7))
    bind_rows(lapply(arms, function(a) {
      i2 <- modifyList(reactiveValuesToList(input),
                       list(regimen = "chemoRT", nfx = a[[2]], dpf = a[[3]],
                            fx_week = a[[4]], gap_len = a[[5]],
                            t_end = 120))
      o <- run_sim(i2); s <- summarise_run(o)
      tibble(arm = a[[1]], onset = s$onset_sev, dur = s$dur_sev,
             peak_area = s$peak_area, BED_tumour = max(o$cBEDt),
             BED_mucosa = max(o$cBEDm))
    }))
  })
  output$tbl_frac <- renderDT(datatable(
    frac() %>% mutate(across(where(is.numeric), ~ round(.x, 2))),
    options = list(dom = "t")))
  output$p_frac <- renderPlot({
    ggplot(frac(), aes(BED_tumour, dur, label = arm)) +
      geom_point(size = 3, colour = "#1f618d") +
      geom_text(vjust = -0.9, size = 3.4) +
      labs(x = "종양 BED (Gy10, 재증식 보정)",
           y = "중증 점막염 지속 기간 (일)") + theme_bw()
  })

  # ---- 12. WHO saturation --------------------------------------------------
  sat <- eventReactive(input$run_sat, {
    bind_rows(lapply(c(0.7, 1.0, 1.3, 1.7), function(sf) {
      i0 <- modifyList(reactiveValuesToList(input),
                       list(sens = sf, cryo_h = 0))
      i1 <- modifyList(reactiveValuesToList(input),
                       list(sens = sf, cryo_h = 6))
      a <- summarise_run(run_sim(i0)); b <- summarise_run(run_sim(i1))
      tibble(sens = sf,
             area_ctrl = a$area_auc, area_trt = b$area_auc,
             d_area_pct = 100 * (b$area_auc - a$area_auc) / a$area_auc,
             who_ctrl = a$dur_sev, who_trt = b$dur_sev,
             d_who_pct = 100 * (b$dur_sev - a$dur_sev) /
               pmax(a$dur_sev, 1e-9))
    }))
  })
  output$tbl_sat <- renderDT(datatable(
    sat() %>% mutate(across(where(is.numeric), ~ round(.x, 2))),
    options = list(dom = "t")))
  output$p_sat <- renderPlot({
    d <- sat() %>% select(sens, `면적 기준 효과 %` = d_area_pct,
                          `WHO 기준 효과 %` = d_who_pct) %>%
      pivot_longer(-sens)
    ggplot(d, aes(sens, value, colour = name)) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      labs(x = "기저 중증도 (sens)", y = "측정된 효과 크기 (%)",
           colour = NULL) + theme_bw() + theme(legend.position = "top")
  })

  # ---- 13. scenarios -------------------------------------------------------
  scen <- eventReactive(input$run_scen, {
    S <- list(
      list("HDM 200, 예방 없음", list(regimen = "HDM", mel_dose = 200)),
      list("HDM 140 (감량)", list(regimen = "HDM", mel_dose = 140)),
      list("HDM 200 + 한랭 30분", list(regimen = "HDM", cryo_h = 0.5)),
      list("HDM 200 + 한랭 6시간", list(regimen = "HDM", cryo_h = 6)),
      list("HDM 200 + 팔리페르민(분리)",
           list(regimen = "HDM", palifermin = TRUE, pal_gap = 1)),
      list("HDM 200 + 팔리페르민(동시)",
           list(regimen = "HDM", palifermin = TRUE, pal_gap = -1)),
      list("HDM 200 + 광생체조절", list(regimen = "HDM", pbm = TRUE)),
      list("HDM 200 + 글루타민", list(regimen = "HDM", glutamine = TRUE)),
      list("HDM 200 + 한랭 + 팔리페르민",
           list(regimen = "HDM", cryo_h = 6, palifermin = TRUE)),
      list("TBI-VP16-Cy (위약)", list(regimen = "TBI")),
      list("TBI-VP16-Cy + 팔리페르민",
           list(regimen = "TBI", palifermin = TRUE)),
      list("두경부 70 Gy + 시스플라틴",
           list(regimen = "chemoRT", t_end = 120)),
      list("두경부 70 Gy 단독",
           list(regimen = "chemoRT", cis_dose = 0, t_end = 120)),
      list("과분할 81.6 Gy/68 fx",
           list(regimen = "chemoRT", nfx = 68, dpf = 1.2, fx_week = 10,
                cis_dose = 0, t_end = 120)),
      list("두경부 항암방사선 + 벤지다민",
           list(regimen = "chemoRT", benzydamine = TRUE, t_end = 120)),
      list("5-FU 볼루스 d1-5", list(regimen = "5FU_bolus", t_end = 40)),
      list("5-FU 96시간 지속주입", list(regimen = "5FU_CI", t_end = 40)))
    bind_rows(lapply(seq_along(S), function(i) {
      i2 <- modifyList(reactiveValuesToList(input), S[[i]][[2]])
      summarise_run(run_sim(i2)) %>% mutate(no = i, scenario = S[[i]][[1]])
    })) %>% relocate(no, scenario)
  })
  output$tbl_scen <- renderDT(datatable(
    scen() %>% mutate(across(where(is.numeric), ~ round(.x, 3))),
    options = list(pageLength = 17, dom = "t")))
  output$p_scen <- renderPlot({
    ggplot(scen(), aes(onset_sev, dur_sev, label = scenario)) +
      geom_point(size = 3, colour = "#b9770e") +
      geom_text(vjust = -0.8, size = 3, check_overlap = TRUE) +
      labs(x = "중증 점막염 발병일 (시계 1: 손상)",
           y = "중증 점막염 지속 기간 (시계 2: 재생)") + theme_bw()
  })

  # ---- 14. virtual population ----------------------------------------------
  pop <- eventReactive(input$run_pop, {
    n <- input$npop
    set.seed(20260806)
    idat <- tibble(
      ID    = 1:n,
      sens  = input$sens * exp(rnorm(n, 0, 0.38) - 0.5 * 0.38 ^ 2),
      lamS  = as.numeric(param(mod)[["lamS"]]) *
                exp(rnorm(n, 0, 0.30) - 0.5 * 0.30 ^ 2),
      Dcrit = 0.34 * exp(rnorm(n, 0, 0.07) - 0.5 * 0.07 ^ 2))
    s <- build_schedule(input)
    s <- bind_rows(lapply(1:n, function(i) mutate(s, ID = i)))
    o <- mod %>% idata_set(idat) %>% data_set(s) %>%
      mrgsim(end = input$t_end, delta = 0.05, recsort = 3) %>% as_tibble()
    o %>% group_by(ID) %>%
      summarise(dur_sev = sum(WHO >= 3) * 0.05,
                peak_who = max(WHO), peak_area = max(ULC),
                .groups = "drop")
  })
  output$tbl_pop <- renderDT({
    d <- pop()
    datatable(tibble(
      `중증(≥3) 발생률` = mean(d$dur_sev > 0),
      `4등급 발생률`    = mean(d$peak_who >= 4),
      `중앙 지속기간(일)` = median(d$dur_sev[d$dur_sev > 0]),
      `평균 최대 면적`   = mean(d$peak_area)) %>%
        mutate(across(everything(), ~ round(.x, 3))),
      options = list(dom = "t"))
  })
  output$p_pop <- renderPlot({
    ggplot(pop(), aes(dur_sev)) +
      geom_histogram(bins = 25, fill = "#5b5b77", colour = "white") +
      labs(x = "중증 점막염 지속 기간 (일)", y = "환자 수") + theme_bw()
  })

  # ---- 15. provenance ------------------------------------------------------
  output$txt_prov <- renderText({
    paste(readLines("om_reference_output.txt", warn = FALSE)[1:60],
          collapse = "\n")
  })
}

shinyApp(ui, server)
