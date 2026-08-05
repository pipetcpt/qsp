## =============================================================================
##  op_shiny_app.R
##  Acute organophosphorus insecticide self-poisoning — QSP dashboard
##  급성 유기인계 살충제 중독 — QSP 대시보드
##
##  Run:
##      library(shiny); library(mrgsolve); library(dplyr); library(tidyr)
##      library(ggplot2); library(DT)
##      shiny::runApp("op_shiny_app.R")
##
##  The app is organised around ONE question: for this patient, at this moment,
##  is the oxime capable of doing anything at all?  That question has a number
##  attached to it —
##
##       OMEGA = k_r_max / (k_i * C_oxon),  free-enzyme ceiling = OMEGA/(1+OMEGA)
##
##  — and every tab is either an input to that number, a consequence of it, or
##  a reminder that the ventilator does not care about it.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

MODEL_FILE <- "op_mrgsolve_model.R"
mod <- mread_cache(file = MODEL_FILE, model = "op")

## ---------------------------------------------------------------------------
## Compound and oxime libraries: they live in the model file's $ENV block, so
## the app and the model can never drift apart.
## ---------------------------------------------------------------------------
MODENV <- tryCatch(mrgsolve::env_get(mod),
                   error = function(e) as.list(mod@envir))
OP_LIBRARY    <- MODENV$OP_LIBRARY
OXIME_LIBRARY <- MODENV$OXIME_LIBRARY
stopifnot(!is.null(OP_LIBRARY), !is.null(OXIME_LIBRARY))

op_params <- function(compound = "chlorpyrifos", oxime = NULL, volume_ml = 50,
                      pon1 = "RR", wt = 60) {
  L <- OP_LIBRARY[[compound]]
  p <- L[c("MW_OP", "KI", "KI_NMJ", "KI_NTE", "T12AGE", "T12SPO", "T12AGNTE",
           "KA_GUT", "FBIO", "V1TH", "VFTH", "QFTH", "VMAXBIO", "KMBIO",
           "CLTHOTH", "CLOXON", "VOXON", "QOXON", "VTOXON", "SOLVFR", "SOLVCV")]
  p$CONC_GL <- L$CONC_GL
  p$DOSE_ML <- volume_ml
  p$PON1SC  <- if (identical(pon1, "QQ")) L$PON1_QQ else 1
  p$WT      <- wt
  if (is.null(oxime) || identical(oxime, "none")) {
    p$OXTYPE <- 0; p$OXRATE <- 0
  } else {
    O <- OXIME_LIBRARY[[oxime]]
    k <- if (identical(O$oxime, "obi")) L$obi else L$pam
    p$KRMAX <- k$KRMAX; p$KDOX <- k$KDOX
    p$OXV1 <- O$V1; p$OXV2 <- O$V2; p$OXCL <- O$CL; p$OXQ <- O$Q
    p$OXMW <- O$MW
    p$OXTYPE <- 1
    p$OXRATE <- if (!is.null(O$inf_mg_kg_h)) O$inf_mg_kg_h * wt * 1000 / O$MW
                else if (!is.null(O$inf_mg_24h)) O$inf_mg_24h / 24 * 1000 / O$MW
                else 0
  }
  p
}

op_events <- function(oxime = NULL, wt = 60, diazepam = TRUE, magnesium = FALSE,
                      scavenger_mg = 0, oxime_start = 1.5, oxime_dur = 168) {
  e <- list()
  if (!is.null(oxime) && !identical(oxime, "none")) {
    O <- OXIME_LIBRARY[[oxime]]
    if (!is.null(O$load_mg_kg))
      e <- c(e, list(ev(time = oxime_start, amt = O$load_mg_kg * wt * 1000 / O$MW,
                        cmt = "A_pam_c")))
    if (!is.null(O$load_mg))
      e <- c(e, list(ev(time = oxime_start, amt = O$load_mg * 1000 / O$MW,
                        cmt = "A_pam_c")))
    if (!is.null(O$bolus_mg))
      e <- c(e, list(ev(time = oxime_start, amt = O$bolus_mg * 1000 / O$MW,
                        cmt = "A_pam_c", ii = O$interval_h,
                        addl = floor(oxime_dur / O$interval_h))))
  }
  if (diazepam)
    e <- c(e, list(ev(time = 0.6, amt = 10 * 1000 / 284.7, cmt = "A_dz",
                      ii = 6, addl = 2)))
  if (magnesium)
    e <- c(e, list(ev(time = 1.0, amt = 16.2, cmt = "A_mg")))
  if (scavenger_mg > 0)
    e <- c(e, list(ev(time = 1.0, amt = scavenger_mg / 1000 / 340000 * 1e6 * 4,
                      cmt = "SCAV")))
  if (!length(e)) return(ev(time = 0, amt = 0, cmt = "A_dz"))
  Reduce(`+`, e)
}

run_case <- function(compound, oxime, volume_ml, pon1, wt, oxime_start,
                     oxime_dur, atr_mode, ventilator, oxygen, icu,
                     glyco, magnesium, scavenger_mg, charcoal_t, tmax = 336) {
  p <- op_params(compound, oxime, volume_ml, pon1, wt)
  p$OXSTART <- oxime_start
  p$OXDUR   <- oxime_dur
  p$ATRMODE <- atr_mode
  p$VENTAVL <- as.numeric(ventilator)
  p$O2AVL   <- as.numeric(oxygen)
  p$ICUSUCT <- as.numeric(icu)
  p$BBBF    <- if (isTRUE(glyco)) 0.02 else 0.35
  p$CHARC_T <- if (is.null(charcoal_t) || charcoal_t < 0) -1 else charcoal_t
  e <- op_events(oxime, wt, TRUE, magnesium, scavenger_mg, oxime_start, oxime_dur)
  mod %>% param(p) %>% ev(e) %>%
    mrgsim(end = tmax, delta = 0.25) %>% as_tibble()
}

## ---------------------------------------------------------------------------
## Closed-form helpers — no ODE needed, these ARE the model's arithmetic
## ---------------------------------------------------------------------------
omega_fn   <- function(krmax, ki, C) krmax / (ki * C)
ceiling_fn <- function(krmax, ki, C) { w <- omega_fn(krmax, ki, C); w / (1 + w) }
achieved_fn <- function(krmax, KD, ki, ks, C, X) {
  kr <- krmax * X / (KD + X); (kr + ks) / (kr + ks + ki * C)
}
Ccrit_fn   <- function(krmax, ki, target = 0.30) (1 - target) / target * krmax / ki
Xreq_fn    <- function(krmax, KD, ki, ks, C, target = 0.30) {
  need <- target / (1 - target) * ki * C - ks
  if (need <= 0) return(0)
  if (need >= krmax) return(NA_real_)      # arithmetically impossible
  KD * need / (krmax - need)
}

THEME <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

## ===========================================================================
## UI
## ===========================================================================
ui <- fluidPage(
  titlePanel("급성 유기인계 살충제 중독 QSP 대시보드 · Acute OP Insecticide Self-Poisoning"),
  tags$p(style = "color:#666;margin-top:-8px;",
         "51-state mrgsolve model · the oxime sufficiency number Ω = k_r_max/(k_i·C_oxon) decides everything the oxime can do"),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① 노출 (Exposure)"),
      selectInput("compound", "화합물 (compound)",
                  choices = c("Chlorpyrifos (diethyl)"   = "chlorpyrifos",
                              "Parathion (diethyl)"      = "parathion",
                              "Dimethoate (dimethyl)"    = "dimethoate",
                              "Fenthion (dimethyl, lipophilic)" = "fenthion"),
                  selected = "chlorpyrifos"),
      sliderInput("volume", "음독량 (mL of concentrate)", 1, 250, 50, step = 1),
      sliderInput("wt", "체중 (kg)", 35, 100, 60, step = 1),
      radioButtons("pon1", "PON1 Q192R 유전형", inline = TRUE,
                   choices = c("R192R (fast)" = "RR", "Q192Q (slow)" = "QQ")),

      hr(), h4("② 해독제 (Antidotes)"),
      selectInput("oxime", "옥심 (oxime)",
                  choices = c("없음 (none)"                       = "none",
                              "Pralidoxime, WHO infusion"         = "pam_who",
                              "Pralidoxime 1 g q6h bolus"         = "pam_bolus",
                              "Obidoxime 250 mg + 750 mg/24h"     = "obidoxime"),
                  selected = "pam_who"),
      sliderInput("oxstart", "옥심 시작 시각 (h)", 0.25, 24, 1.5, step = 0.25),
      sliderInput("oxdur", "옥심 투여 기간 (h)", 12, 336, 168, step = 12),
      radioButtons("atr", "아트로핀 적정 (atropine)", inline = FALSE,
                   choices = c("빠른 배가 프로토콜 (rapid doubling)" = 1,
                               "관행적 적정 (slow, ad hoc)"          = 2,
                               "없음 (none)"                          = 0),
                   selected = 1),
      checkboxInput("glyco", "글리코피롤레이트로 대체 (BBB 비투과)", FALSE),
      checkboxInput("mg", "황산마그네슘 4 g", FALSE),
      sliderInput("scav", "혈장 BChE 생포집제 (mg)", 0, 5000, 0, step = 250),
      sliderInput("charcoal", "활성탄 투여 시각 (h; -1 = 미투여)", -1, 12, -1, step = 0.5),

      hr(), h4("③ 지지요법 (Supportive care)"),
      checkboxInput("vent", "인공환기 가용 (ventilator available)", TRUE),
      checkboxInput("o2", "산소 투여 가능", TRUE),
      checkboxInput("icu", "ICU 기도 간호 (suctioning)", TRUE),

      hr(),
      actionButton("go", "시뮬레이션 실행 (Simulate)", class = "btn-primary",
                   width = "100%"),
      tags$p(style = "font-size:11px;color:#888;margin-top:10px;",
             "교육·연구 목적 전용. 임상 의사결정에 사용하지 마십시오.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## ------------------------------------------------------------------
        tabPanel("① 환자 프로파일",
          br(),
          fluidRow(
            column(3, wellPanel(h5("음독 몰수"), textOutput("kpi_dose"))),
            column(3, wellPanel(h5("옥손 평형 농도"), textOutput("kpi_oxon"))),
            column(3, wellPanel(h5("Ω (24 h)"), textOutput("kpi_omega"))),
            column(3, wellPanel(h5("옥심 천장"), textOutput("kpi_ceiling")))
          ),
          fluidRow(
            column(3, wellPanel(h5("AChE 최저치"), textOutput("kpi_nadir"))),
            column(3, wellPanel(h5("누적 아트로핀"), textOutput("kpi_atr"))),
            column(3, wellPanel(h5("인공환기 시간"), textOutput("kpi_vent"))),
            column(3, wellPanel(h5("14일 사망 확률"), textOutput("kpi_mort")))
          ),
          hr(),
          h4("이 환자에서 옥심이 할 수 있는 일 (what the oxime can do here)"),
          verbatimTextOutput("verdict"),
          hr(),
          plotOutput("plot_profile", height = "420px")
        ),

        ## ------------------------------------------------------------------
        tabPanel("② 독성동태 (PK/TK)",
          br(),
          plotOutput("plot_tk", height = "700px"),
          hr(),
          tags$p("모(母)화합물(thion)은 AChE를 억제하지 않습니다. 억제하는 것은 CYP가 만든 옥손이며, ",
                 "그 생성은 포화합니다. 따라서 대량 음독에서 옥손 고원 농도는 Vmax/CL_oxon 으로 ",
                 "수렴하고 용량과 무관해집니다 — 옥심의 천장이 대량 음독에서 상수가 되는 이유입니다.")
        ),

        ## ------------------------------------------------------------------
        tabPanel("③ 에스터라제 스위치",
          br(),
          plotOutput("plot_esterase", height = "760px"),
          hr(),
          h5("세 상태의 물질수지 (three-state mass balance)"),
          tags$pre("E  --k_i*C_oxon-->  EP  --k_a-->  EP_aged   (되돌릴 수 없음)\nE  <--k_s+k_r[X]--  EP\n\n노화 손실은 사건이 아니라 적분입니다:  aged(T) = 1 - exp(-k_a * ∫ f_EP dt)")
        ),

        ## ------------------------------------------------------------------
        tabPanel("④ 콜린성 증후군",
          br(),
          plotOutput("plot_signs", height = "760px"),
          hr(),
          tags$p("아트로핀 필요량은 (1-E)가 아니라 1/E 로 커집니다. 시냅스 ACh 의 정상상태가 ",
                 "(1+leak)/(E+leak) 라는 쌍곡선이기 때문이며, 경쟁적 길항이므로 필요 농도는 ACh 에 ",
                 "선형입니다. AChE 2%인 환자가 20%인 환자보다 아트로핀을 몇 배 더 쓰는 것은 ",
                 "중증도의 '정도' 문제가 아니라 산술입니다.")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑤ 호흡 · 임상 종말점",
          br(),
          plotOutput("plot_resp", height = "760px"),
          hr(),
          tags$p("환기량은 '뇌가 요구하는 양'과 '근육이 낼 수 있는 양'의 최솟값입니다. ",
                 "안정 시 환기에는 최대 흡기력의 약 28%만 필요하므로(신경근 안전계수), ",
                 "신경근 차단은 중추 억제보다 훨씬 큰 여유를 가집니다. 이 비대칭이 ",
                 "'조기 사망은 중추성'이라는 임상 관찰의 모델 표현입니다.")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑥ 시나리오 비교",
          br(),
          fluidRow(
            column(6, checkboxGroupInput("scen", "비교할 시나리오",
              choices = c("지지요법만 (supportive only)"       = "sup",
                          "2-PAM WHO 주입"                      = "pam_who",
                          "2-PAM 1 g q6h 볼루스"                = "pam_bolus",
                          "오비독심"                             = "obidoxime",
                          "환기기 없음 (2-PAM)"                  = "novent",
                          "느린 아트로핀 적정 (2-PAM)"           = "slowatr",
                          "글리코피롤레이트 (2-PAM)"             = "glyco",
                          "활성탄 1시간 (2-PAM)"                 = "char1",
                          "PON1 Q192Q (2-PAM)"                   = "ponqq"),
              selected = c("sup", "pam_who", "novent", "slowatr"))),
            column(6, radioButtons("scen_y", "표시할 변수",
              choices = c("RBC AChE (%)"        = "AChE_RBC",
                          "누적 사망확률 (%)"    = "MORTALITY",
                          "환기 가능도 VCAP"     = "VCAP",
                          "누적 아트로핀 (mg)"   = "ATRCUM",
                          "POP 중증도 점수"      = "POP",
                          "옥손 (nM)"            = "C_oxon"),
              selected = "AChE_RBC"))
          ),
          plotOutput("plot_scen", height = "440px"),
          hr(),
          DTOutput("tbl_scen")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑦ 옥심 천장 탐색기 (Ω)",
          br(),
          fluidRow(
            column(4, sliderInput("ex_oxon", "혈장 옥손 (nM, log 축)",
                                  min = 0, max = 4, value = 2.1, step = 0.05,
                                  pre = "10^")),
            column(4, sliderInput("ex_target", "목표 유리효소 분율", 0.05, 0.60,
                                  0.30, step = 0.05)),
            column(4, sliderInput("ex_X", "옥심 농도 (µM)", 0, 400, 130, step = 5))
          ),
          plotOutput("plot_omega", height = "420px"),
          hr(),
          DTOutput("tbl_omega"),
          hr(),
          tags$p(strong("읽는 법: "),
            "곡선이 목표선 위로 올라오지 못하는 구간에서는 옥심의 '용량 부족'이 아니라 ",
            "옥심이라는 약물 종류로는 도달이 불가능합니다. 어떤 무작위 시험도 그 구간의 ",
            "환자에서는 이득을 보일 수 없습니다.")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑧ 용량-반응 스윕",
          br(),
          actionButton("go_sweep", "스윕 실행 (8 doses × 2 arms)", class = "btn-info"),
          br(), br(),
          plotOutput("plot_sweep", height = "440px"),
          DTOutput("tbl_sweep")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑨ 가상 임상시험",
          br(),
          fluidRow(
            column(3, numericInput("n_trial", "환자 수", 40, min = 10, max = 200, step = 10)),
            column(3, sliderInput("trial_med", "음독량 중앙값 (mL)", 5, 80, 28, step = 1)),
            column(3, sliderInput("trial_cpf", "diethyl OP 비율", 0, 1, 0.55, step = 0.05)),
            column(3, br(), actionButton("go_trial", "시험 실행", class = "btn-warning"))
          ),
          plotOutput("plot_trial", height = "400px"),
          hr(),
          verbatimTextOutput("txt_trial"),
          tags$p("실사용 무작위 시험은 도착한 환자를 그대로 등록합니다. 천장 위의 환자들이 ",
                 "천장 아래의 실제 효과를 희석하면, 옥심이 무력하다고 가정하지 않고도 ",
                 "음성 결과가 나옵니다.")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑩ 바이오마커 · 지연 증후군",
          br(),
          plotOutput("plot_biomarker", height = "700px"),
          hr(),
          tags$p("적혈구 AChE는 재합성되지 않습니다(적혈구 수명으로만 교체). 근육·뇌 AChE는 ",
                 "약 5일 반감기로 재합성됩니다. 회복기에 RBC AChE 와 임상 상태가 어긋나는 것은 ",
                 "측정 오류가 아니라 두 구획의 회전율이 다르기 때문입니다.")
        ),

        ## ------------------------------------------------------------------
        tabPanel("⑪ 문서",
          br(),
          uiOutput("docs")
        )
      )
    )
  )
)

## ---------------------------------------------------------------------------
## small helper so the app does not hard-depend on patchwork/gridExtra
## ---------------------------------------------------------------------------
gridExtra_arrange <- function(...) {
  ps <- list(...)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(ps, ncol = 1)
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(grobs = ps, ncol = 1)
  } else {
    ps[[1]]
  }
}

## ===========================================================================
## SERVER
## ===========================================================================
server <- function(input, output, session) {

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    withProgress(message = "51-state ODE 적분 중...", value = 0.5, {
      run_case(input$compound,
               if (identical(input$oxime, "none")) NULL else input$oxime,
               input$volume, input$pon1, input$wt, input$oxstart, input$oxdur,
               as.numeric(input$atr), input$vent, input$o2, input$icu,
               input$glyco, input$mg, input$scav, input$charcoal)
    })
  })

  cmp <- reactive(OP_LIBRARY[[input$compound]])

  ## ---------------------------------------------------------------- KPIs
  output$kpi_dose <- renderText({
    L <- cmp()
    mmol <- input$volume * L$CONC_GL / 1000 / L$MW_OP * 1000
    sprintf("%.1f mmol (%.1f g)", mmol, input$volume * L$CONC_GL / 1000)
  })
  output$kpi_oxon <- renderText({
    d <- sim(); sprintf("%.0f nM (24 h)", approx(d$time, d$C_oxon, 24)$y)
  })
  output$kpi_omega <- renderText({
    d <- sim(); sprintf("%.2f", approx(d$time, d$OMEGA, 24)$y)
  })
  output$kpi_ceiling <- renderText({
    d <- sim(); sprintf("%.0f %% 유리효소", 100 * approx(d$time, d$E_CEIL, 24)$y)
  })
  output$kpi_nadir <- renderText({ sprintf("%.2f %%", min(sim()$AChE_RBC)) })
  output$kpi_atr   <- renderText({ sprintf("%.0f mg / 14일", max(sim()$ATRCUM)) })
  output$kpi_vent  <- renderText({ sprintf("%.0f h", max(sim()$VTIME)) })
  output$kpi_mort  <- renderText({ sprintf("%.1f %%", max(sim()$MORTALITY)) })

  output$verdict <- renderText({
    d <- sim(); L <- cmp()
    C  <- approx(d$time, d$C_oxon, 24)$y
    krmax <- if (identical(input$oxime, "obidoxime")) L$obi$KRMAX else L$pam$KRMAX
    KD    <- if (identical(input$oxime, "obidoxime")) L$obi$KDOX  else L$pam$KDOX
    ks <- log(2) / L$T12SPO
    ceil <- ceiling_fn(krmax, L$KI, C)
    Xr30 <- Xreq_fn(krmax, KD, L$KI, ks, C, 0.30)
    Xclin <- if (identical(input$oxime, "obidoxime")) 20 else 130
    ach  <- achieved_fn(krmax, KD, L$KI, ks, C, Xclin)
    phi  <- (log(2) / L$T12AGE) /
      (log(2) / L$T12AGE + ks + krmax * Xclin / (KD + Xclin))
    paste0(
      sprintf("24시간 시점 옥손 %.0f nM,  k_i·C_oxon = %.1f /h,  k_r_max = %.0f /h\n",
              C, L$KI * C, krmax),
      sprintf("Ω = %.2f  →  어떤 용량으로도 도달 가능한 최대 유리효소 = %.0f%%\n",
              omega_fn(krmax, L$KI, C), 100 * ceil),
      sprintf("임상적으로 도달 가능한 농도(%.0f µM)에서 실제 유리효소 = %.0f%%\n",
              Xclin, 100 * ach),
      if (is.na(Xr30))
        "→ 유리효소 30%는 이 옥심으로는 산술적으로 불가능합니다 (k_r_max 천장 아래).\n"
      else if (Xr30 > Xclin)
        sprintf("→ 30%%에 도달하려면 %.0f µM 이 필요하나 임상 도달 농도는 %.0f µM 입니다.\n",
                Xr30, Xclin)
      else
        sprintf("→ 30%% 도달에 필요한 농도는 %.0f µM 로, 현행 용법으로 도달 가능합니다.\n",
                Xr30),
      sprintf("노화 래칫 φ = %.3f (억제 사건 1회가 비가역으로 끝날 확률)\n", phi),
      sprintf("14일 예측 사망확률 %.1f%%,  인공환기 %.0f h,  누적 아트로핀 %.0f mg",
              max(d$MORTALITY), max(d$VTIME), max(d$ATRCUM))
    )
  })

  ## -------------------------------------------------------------- ① profile
  output$plot_profile <- renderPlot({
    d <- sim()
    d %>% select(time, `RBC AChE (%)` = AChE_RBC, `Ω` = OMEGA,
                 `환기 가능도 VCAP` = VCAP, `사망확률 (%)` = MORTALITY) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 1.0, colour = "#1f6fb2") +
      facet_wrap(~name, scales = "free_y") +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "시간 (h)", y = NULL,
           title = "핵심 궤적 (core trajectories)") + THEME
  })

  ## ------------------------------------------------------------------ ② TK
  output$plot_tk <- renderPlot({
    d <- sim()
    bind_rows(
      tibble(time = d$time, value = d$C_thion,  panel = "모화합물 thion (µM)"),
      tibble(time = d$time, value = d$C_oxon,   panel = "혈장 옥손 (nM)"),
      tibble(time = d$time, value = d$C_oxonT,  panel = "조직 옥손 (nM)"),
      tibble(time = d$time, value = d$C_oxime_, panel = "옥심 (µM)"),
      tibble(time = d$time, value = d$Ce_atr,   panel = "아트로핀 효과부위 (nM)"),
      tibble(time = d$time, value = d$kinh_,    panel = "억제 속도 k_i·C_oxon (1/h)"),
      tibble(time = d$time, value = d$kr_,      panel = "재활성화 속도 k_r[X] (1/h)"),
      tibble(time = d$time, value = d$OMEGA,    panel = "Ω = k_r_max/(k_i·C_oxon)")
    ) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.9, colour = "#12507b") +
      facet_wrap(~panel, scales = "free_y", ncol = 2) +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "시간 (h)", y = NULL,
           title = "독성동태와 두 속도상수 (toxicokinetics and the two rate constants)") +
      THEME
  })

  ## ----------------------------------------------------------- ③ esterase
  output$plot_esterase <- renderPlot({
    d <- sim()
    p1 <- d %>% select(time, `유리 E` = E_rbc, `인산화 EP` = EP_rbc,
                       `노화 EP-aged` = EA_rbc) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, 100 * value, fill = name)) +
      geom_area(alpha = 0.85) +
      scale_fill_manual(values = c("유리 E" = "#2e9e5b", "인산화 EP" = "#f0a202",
                                   "노화 EP-aged" = "#c0392b")) +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "시간 (h)", y = "적혈구 AChE 분율 (%)", fill = NULL,
           title = "적혈구 AChE 의 세 상태 — 붉은 영역은 되돌릴 수 없습니다") + THEME
    p2 <- d %>% select(time, `RBC AChE` = E_rbc, `NMJ AChE` = E_nmj,
                       `뇌 AChE` = E_cns, `혈장 BChE` = B_free) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, 100 * value, colour = name)) +
      geom_line(linewidth = 1.0) +
      geom_hline(yintercept = 30, linetype = 2, colour = "#888") +
      annotate("text", x = 250, y = 33, label = "증상 발현 문턱 ~30%",
               size = 3.4, colour = "#666") +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "시간 (h)", y = "기준선 대비 (%)", colour = NULL,
           title = "구획별 에스터라제 — 옥심은 말초에는 닿고 뇌에는 거의 닿지 않습니다") +
      THEME
    p3 <- d %>% select(time, `Ω` = OMEGA, `천장 E_CEIL` = E_CEIL,
                       `준정상상태 E_QSS` = E_QSS, `노화 래칫 φ` = PHI) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#8e44ad") +
      facet_wrap(~name, scales = "free_y", nrow = 1) +
      scale_x_continuous(breaks = c(0, 72, 168, 336)) +
      labs(x = "시간 (h)", y = NULL, title = "모델을 지배하는 두 무차원수") + THEME
    gridExtra_arrange(p1, p2, p3)
  })

  ## -------------------------------------------------------------- ④ signs
  output$plot_signs <- renderPlot({
    d <- sim()
    bind_rows(
      tibble(time = d$time, value = d$ACh_m,   panel = "무스카린성 시냅스 ACh (x 정상)"),
      tibble(time = d$time, value = d$MUSX,    panel = "무스카린성 신호 초과 (0-1)"),
      tibble(time = d$time, value = d$NICX,    panel = "니코틴성 신호 초과 (0-1)"),
      tibble(time = d$time, value = d$CNSX,    panel = "중추 무스카린성 초과 (0-1)"),
      tibble(time = d$time, value = d$SEC,     panel = "기도 분비물 (mL)"),
      tibble(time = d$time, value = d$HR,      panel = "심박수 (bpm)"),
      tibble(time = d$time, value = d$MAP,     panel = "평균동맥압 (mmHg)"),
      tibble(time = d$time, value = d$ATRCUM,  panel = "누적 아트로핀 (mg)"),
      tibble(time = d$time, value = d$Rm_down, panel = "수용체 하향조절 (내성)"),
      tibble(time = d$time, value = d$POP,     panel = "POP 중증도 점수 (0-11)")
    ) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#b7791f") +
      facet_wrap(~panel, scales = "free_y", ncol = 2) +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "시간 (h)", y = NULL,
           title = "콜린성 증후군과 폐쇄루프 아트로핀") + THEME
  })

  ## --------------------------------------------------------------- ⑤ resp
  output$plot_resp <- renderPlot({
    d <- sim()
    bind_rows(
      tibble(time = d$time, value = d$RESPD,     panel = "중추 호흡 구동 RESPD"),
      tibble(time = d$time, value = d$MSTR,      panel = "호흡근력 MSTR"),
      tibble(time = d$time, value = d$VCAP,      panel = "자발 환기 가능도 VCAP"),
      tibble(time = d$time, value = d$VENT_ON,   panel = "인공환기 가동 (0/1)"),
      tibble(time = d$time, value = d$PACO2,     panel = "PaCO₂ (mmHg)"),
      tibble(time = d$time, value = d$PAO2_,     panel = "PaO₂ (mmHg)"),
      tibble(time = d$time, value = d$LUNG,      panel = "흡인성 폐손상 (0-1)"),
      tibble(time = d$time, value = d$MORTALITY, panel = "누적 사망확률 (%)")
    ) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#a03e6f") +
      facet_wrap(~panel, scales = "free_y", ncol = 2) +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "시간 (h)", y = NULL,
           title = "호흡부전 — 이 병의 최종 공통 경로") + THEME
  })

  ## ----------------------------------------------------------- ⑥ scenarios
  scen_runs <- reactive({
    base <- list(compound = input$compound, volume_ml = input$volume,
                 pon1 = input$pon1, wt = input$wt, oxime_start = input$oxstart,
                 oxime_dur = input$oxdur, atr_mode = 1, ventilator = TRUE,
                 oxygen = TRUE, icu = TRUE, glyco = FALSE, magnesium = FALSE,
                 scavenger_mg = 0, charcoal_t = -1)
    defs <- list(
      sup       = modifyList(base, list(oxime = NULL)),
      pam_who   = modifyList(base, list(oxime = "pam_who")),
      pam_bolus = modifyList(base, list(oxime = "pam_bolus")),
      obidoxime = modifyList(base, list(oxime = "obidoxime")),
      novent    = modifyList(base, list(oxime = "pam_who", ventilator = FALSE,
                                        icu = FALSE)),
      slowatr   = modifyList(base, list(oxime = "pam_who", atr_mode = 2)),
      glyco     = modifyList(base, list(oxime = "pam_who", glyco = TRUE)),
      char1     = modifyList(base, list(oxime = "pam_who", charcoal_t = 1)),
      ponqq     = modifyList(base, list(oxime = "pam_who", pon1 = "QQ"))
    )
    sel <- intersect(input$scen, names(defs))
    withProgress(message = "시나리오 시뮬레이션", value = 0, {
      out <- lapply(seq_along(sel), function(i) {
        incProgress(1 / length(sel), detail = sel[i])
        a <- defs[[sel[i]]]
        do.call(run_case, a) %>% mutate(scenario = sel[i])
      })
    })
    bind_rows(out)
  })

  output$plot_scen <- renderPlot({
    d <- scen_runs()
    req(nrow(d) > 0)
    d$y <- d[[input$scen_y]]
    ggplot(d, aes(time, y, colour = scenario)) +
      geom_line(linewidth = 1.0) +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "시간 (h)", y = input$scen_y, colour = NULL,
           title = "시나리오 비교 (scenario comparison)") + THEME
  })

  output$tbl_scen <- renderDT({
    d <- scen_runs(); req(nrow(d) > 0)
    d %>% group_by(scenario) %>%
      summarise(`AChE 최저 (%)` = round(min(AChE_RBC), 2),
                `AChE 72h (%)`  = round(approx(time, AChE_RBC, 72)$y, 1),
                `노화 72h (%)`  = round(approx(time, AGED_RBC, 72)$y, 1),
                `Ω 24h`         = round(approx(time, OMEGA, 24)$y, 2),
                `아트로핀 (mg)` = round(max(ATRCUM)),
                `환기 (h)`      = round(max(VTIME)),
                `IMS 지수 (%)`  = round(max(IMS_IDX)),
                `POP 최대`      = round(max(POP), 1),
                `사망확률 (%)`  = round(max(MORTALITY), 1), .groups = "drop") %>%
      datatable(options = list(dom = "t", pageLength = 12), rownames = FALSE)
  })

  ## --------------------------------------------------------- ⑦ Ω explorer
  omega_tbl <- reactive({
    C <- 10^input$ex_oxon
    rows <- lapply(names(OP_LIBRARY), function(nm) {
      L <- OP_LIBRARY[[nm]]; ks <- log(2) / L$T12SPO
      lapply(c("pam", "obi"), function(ox) {
        k <- if (ox == "obi") L$obi else L$pam
        Xr <- Xreq_fn(k$KRMAX, k$KDOX, L$KI, ks, C, input$ex_target)
        data.frame(
          compound = nm, oxime = ifelse(ox == "obi", "obidoxime", "2-PAM"),
          k_i = L$KI, k_r_max = k$KRMAX,
          Omega = round(omega_fn(k$KRMAX, L$KI, C), 3),
          ceiling_pct = round(100 * ceiling_fn(k$KRMAX, L$KI, C), 1),
          achieved_pct = round(100 * achieved_fn(k$KRMAX, k$KDOX, L$KI, ks, C,
                                                 input$ex_X), 1),
          C_crit_nM = round(Ccrit_fn(k$KRMAX, L$KI, input$ex_target)),
          X_required_uM = ifelse(is.na(Xr), NA, round(Xr)),
          verdict = ifelse(is.na(Xr), "산술적으로 불가능",
                    ifelse(Xr <= input$ex_X, "현재 농도로 도달", "더 높은 농도 필요")))
      }) %>% bind_rows()
    }) %>% bind_rows()
    rows
  })

  output$plot_omega <- renderPlot({
    Cg <- 10^seq(0, 4, length.out = 240)
    dat <- lapply(names(OP_LIBRARY), function(nm) {
      L <- OP_LIBRARY[[nm]]; ks <- log(2) / L$T12SPO
      bind_rows(
        data.frame(compound = nm, oxime = "2-PAM", C = Cg,
                   ceiling = ceiling_fn(L$pam$KRMAX, L$KI, Cg),
                   achieved = achieved_fn(L$pam$KRMAX, L$pam$KDOX, L$KI, ks,
                                          Cg, input$ex_X)),
        data.frame(compound = nm, oxime = "obidoxime", C = Cg,
                   ceiling = ceiling_fn(L$obi$KRMAX, L$KI, Cg),
                   achieved = achieved_fn(L$obi$KRMAX, L$obi$KDOX, L$KI, ks,
                                          Cg, input$ex_X)))
    }) %>% bind_rows() %>%
      pivot_longer(c(ceiling, achieved), names_to = "which", values_to = "E")
    ggplot(dat, aes(C, 100 * E, colour = oxime, linetype = which)) +
      geom_line(linewidth = 0.95) +
      geom_hline(yintercept = 100 * input$ex_target, linetype = 3,
                 colour = "#c0392b", linewidth = 0.8) +
      geom_vline(xintercept = 10^input$ex_oxon, linetype = 2, colour = "#555") +
      facet_wrap(~compound, nrow = 1) +
      scale_x_log10() +
      scale_linetype_manual(values = c(ceiling = 1, achieved = 2),
                            labels = c(achieved = "현재 옥심 농도에서 실제",
                                       ceiling = "무한 용량에서의 천장")) +
      labs(x = "혈장 옥손 농도 (nM, log)", y = "유리 AChE (%)",
           colour = NULL, linetype = NULL,
           title = "옥심 천장 — 실선 아래로 목표선이 지나가면 그 환자에게 옥심은 약이 아닙니다") +
      THEME
  })

  output$tbl_omega <- renderDT({
    datatable(omega_tbl(), options = list(dom = "t", pageLength = 10),
              rownames = FALSE)
  })

  ## ------------------------------------------------------------- ⑧ sweep
  sweep_res <- eventReactive(input$go_sweep, {
    vols <- c(2, 5, 10, 20, 35, 50, 100, 200)
    withProgress(message = "용량-반응 스윕", value = 0, {
      lapply(vols, function(v) {
        incProgress(1 / length(vols), detail = paste0(v, " mL"))
        a <- run_case(input$compound, NULL, v, input$pon1, input$wt, 1.5, 168,
                      1, TRUE, TRUE, TRUE, FALSE, FALSE, 0, -1)
        b <- run_case(input$compound, "pam_who", v, input$pon1, input$wt, 1.5,
                      168, 1, TRUE, TRUE, TRUE, FALSE, FALSE, 0, -1)
        tibble(volume_ml = v,
               oxon_24h  = approx(b$time, b$C_oxon, 24)$y,
               Omega     = approx(b$time, b$OMEGA, 24)$y,
               ceiling   = 100 * approx(b$time, b$E_CEIL, 24)$y,
               AChE72_supportive = approx(a$time, a$AChE_RBC, 72)$y,
               AChE72_oxime      = approx(b$time, b$AChE_RBC, 72)$y,
               mort_supportive   = max(a$MORTALITY),
               mort_oxime        = max(b$MORTALITY))
      }) %>% bind_rows()
    })
  })

  output$plot_sweep <- renderPlot({
    d <- sweep_res()
    d %>% select(volume_ml, `지지요법` = AChE72_supportive,
                 `2-PAM` = AChE72_oxime, `천장` = ceiling) %>%
      pivot_longer(-volume_ml) %>%
      ggplot(aes(volume_ml, value, colour = name)) +
      geom_line(linewidth = 1.0) + geom_point(size = 2) +
      geom_hline(yintercept = 30, linetype = 2, colour = "#888") +
      scale_x_log10() +
      labs(x = "음독량 (mL, log)", y = "72시간 RBC AChE (%)", colour = NULL,
           title = "옥심 효과는 용량이 커질수록 사라진다 — 약 때문이 아니라 옥손 때문에") +
      THEME
  })

  output$tbl_sweep <- renderDT({
    datatable(sweep_res() %>% mutate(across(where(is.numeric), ~round(.x, 2))),
              options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  ## ------------------------------------------------------------- ⑨ trial
  trial_res <- eventReactive(input$go_trial, {
    set.seed(20260805)
    n <- input$n_trial
    withProgress(message = "가상 시험", value = 0, {
      lapply(seq_len(n), function(i) {
        incProgress(1 / n, detail = paste0(i, "/", n))
        ml  <- min(250, max(1, exp(rnorm(1, log(input$trial_med), 0.85))))
        cmpn <- if (runif(1) < input$trial_cpf) "chlorpyrifos" else "dimethoate"
        pon <- if (runif(1) < 0.25) "QQ" else "RR"
        a <- run_case(cmpn, NULL, ml, pon, 60, 1.5, 168, 1, TRUE, TRUE, TRUE,
                      FALSE, FALSE, 0, -1)
        b <- run_case(cmpn, "pam_who", ml, pon, 60, 1.5, 168, 1, TRUE, TRUE,
                      TRUE, FALSE, FALSE, 0, -1)
        tibble(id = i, op = cmpn, ml = ml, pon1 = pon,
               mort_placebo = max(a$MORTALITY), mort_oxime = max(b$MORTALITY),
               vent_placebo = max(a$VTIME),     vent_oxime = max(b$VTIME))
      }) %>% bind_rows()
    })
  })

  output$plot_trial <- renderPlot({
    d <- trial_res()
    d %>% mutate(delta = mort_placebo - mort_oxime) %>%
      ggplot(aes(ml, delta, colour = op)) +
      geom_hline(yintercept = 0, colour = "#999") +
      geom_point(size = 2.4, alpha = 0.85) +
      geom_smooth(se = FALSE, method = "loess", formula = y ~ x,
                  linewidth = 0.8) +
      scale_x_log10() +
      labs(x = "음독량 (mL, log)", y = "옥심에 의한 사망확률 감소 (%p)",
           colour = NULL,
           title = "개별 환자에서의 옥심 효과 — 0 위와 아래가 한 시험 안에 섞인다") +
      THEME
  })

  output$txt_trial <- renderText({
    d <- trial_res()
    strat <- d %>% mutate(s = ifelse(ml <= 15, "<=15 mL", ">15 mL")) %>%
      group_by(op, s) %>%
      summarise(n = n(), pl = mean(mort_placebo), ox = mean(mort_oxime),
                .groups = "drop")
    paste0(
      sprintf("n = %d  |  전체 위약군 사망 %.1f%%  옥심군 %.1f%%  RR = %.2f\n",
              nrow(d), mean(d$mort_placebo), mean(d$mort_oxime),
              mean(d$mort_oxime) / mean(d$mort_placebo)),
      sprintf("인공환기 시간 위약군 %.0f h  옥심군 %.0f h\n",
              mean(d$vent_placebo), mean(d$vent_oxime)),
      paste(apply(strat, 1, function(r)
        sprintf("  %-14s %-8s n=%-4s 위약 %5.1f%%  옥심 %5.1f%%",
                r[["op"]], r[["s"]], r[["n"]],
                as.numeric(r[["pl"]]), as.numeric(r[["ox"]]))), collapse = "\n")
    )
  })

  ## --------------------------------------------------------- ⑩ biomarkers
  output$plot_biomarker <- renderPlot({
    d <- sim()
    bind_rows(
      tibble(time = d$time, value = d$AChE_RBC, panel = "적혈구 AChE (%)"),
      tibble(time = d$time, value = d$BChE_PL,  panel = "혈장 BChE (%)"),
      tibble(time = d$time, value = d$AGED_RBC, panel = "비가역 노화 분율 (%)"),
      tibble(time = d$time, value = 100 * d$E_nmj, panel = "NMJ AChE (%)"),
      tibble(time = d$time, value = d$IMS_IDX,  panel = "중간증후군 지수 (nAChR 차단 %)"),
      tibble(time = d$time, value = 100 * d$N_aged, panel = "노화 NTE (%) — OPIDN 전구"),
      tibble(time = d$time, value = 100 * d$OPIDN, panel = "OPIDN 점수 (%)"),
      tibble(time = d$time, value = d$POP, panel = "POP 중증도 점수")
    ) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#2e7d5b") +
      facet_wrap(~panel, scales = "free_y", ncol = 2) +
      scale_x_continuous(breaks = c(0, 24, 72, 168, 336)) +
      labs(x = "시간 (h)", y = NULL,
           title = "바이오마커와 지연 증후군") + THEME
  })

  ## ---------------------------------------------------------------- ⑪ docs
  output$docs <- renderUI({
    if (file.exists("README.md")) includeMarkdown("README.md")
    else tags$p("README.md 를 찾을 수 없습니다.")
  })
}

shinyApp(ui, server)
