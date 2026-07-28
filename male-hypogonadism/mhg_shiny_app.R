## =============================================================================
##  Male Hypogonadism (남성 성선기능저하증) — QSP Shiny dashboard
##  mhg_shiny_app.R
##
##  Run:
##    setwd("male-hypogonadism"); shiny::runApp("mhg_shiny_app.R")
##
##  Requires mhg_mrgsolve_model.R in the same directory (it is sourced below,
##  which compiles the mrgsolve model once at start-up).
##
##  The dashboard is organised around the three claims the model exists to
##  make, rather than around the compartment list:
##    - the measured variable is not the acting variable   (tabs 2, 3)
##    - the waveform, not just the dose, sets the harm     (tabs 4, 5)
##    - serum T and intratesticular T are different things (tab 6)
##  The remaining tabs carry the organ-level endpoints and the trial ledger.
##
##  DISCLAIMER: research and education only. Not for clinical use.
## =============================================================================

library(shiny)
library(mrgsolve)

source("mhg_mrgsolve_model.R", local = FALSE)

## ---------------------------------------------------------------- helpers ----
`%||%` <- function(a, b) if (is.null(a)) b else a

PALETTE <- c("#2c6fb5", "#1e8449", "#c0392b", "#8e44ad",
             "#d68910", "#16a085", "#7f8c8d", "#c2185b")

REGIMENS <- list(
  "None (untreated)"                 = "none",
  "Gel 1.62% 81 mg/day"              = "gel",
  "SC auto-injector 75 mg weekly"    = "sc",
  "IM cypionate 100 mg weekly"       = "im_wk",
  "IM cypionate 200 mg q2wk"         = "im_q2wk",
  "IM undecanoate 1000 mg q12wk"     = "tu_im",
  "Oral undecanoate 237 mg BID"      = "oral",
  "Pellets 750 mg q4mo"              = "pellet",
  "Nasal gel 11 mg TID"              = "nasal",
  "hCG 1500 IU 3x/week"              = "hcg",
  "Gel + hCG 500 IU EOD"             = "gel_hcg",
  "Clomiphene 25 mg daily"           = "clomi"
)

PATIENTS <- list(
  "Primary (organic testicular failure)" = "organic",
  "Functional (obesity-driven)"          = "functional",
  "Elderly (age-related, high SHBG)"     = "elderly",
  "Secondary (pituitary/hypothalamic)"   = "secondary",
  "Klinefelter 47,XXY"                   = "klinefelter",
  "Opioid-induced androgen deficiency"   = "opioid"
)

pt_pars <- function(key) switch(key,
  organic     = PT_ORGANIC,
  functional  = PT_FUNCTIONAL,
  elderly     = PT_ELDERLY,
  secondary   = PT_SECONDARY,
  klinefelter = PT_KLINEFELTER,
  opioid      = PT_OPIOID,
  PT_ORGANIC)

make_events <- function(key, weeks) switch(key,
  none    = NULL,
  gel     = tst_gel(81, weeks = weeks),
  sc      = tst_sc(75, weeks = weeks),
  im_wk   = tst_im(100, 7, weeks = weeks),
  im_q2wk = tst_im(200, 14, weeks = weeks),
  tu_im   = tst_tu_im(1000, 84, weeks = weeks),
  oral    = tst_oral(237, weeks = weeks),
  pellet  = tst_pellet(750, 120, weeks = weeks),
  nasal   = tst_nasal(11, weeks = weeks),
  hcg     = hcg(1500, 7/3, weeks = weeks),
  gel_hcg = c(tst_gel(81, weeks = weeks), hcg(500, 2, weeks = weeks)),
  clomi   = clomiphene(25, weeks = weeks),
  NULL)

lineplot <- function(dl, var, ylab, main, hlines = NULL, cols = PALETTE) {
  rng <- range(unlist(lapply(dl, function(d) d[[var]])), na.rm = TRUE)
  if (!is.null(hlines)) rng <- range(c(rng, hlines))
  if (diff(rng) == 0) rng <- rng + c(-1, 1)
  plot(NA, xlim = range(dl[[1]]$day), ylim = rng, xlab = "days from start",
       ylab = ylab, main = main, las = 1)
  grid(col = "grey90")
  if (!is.null(hlines)) abline(h = hlines, lty = 2, col = "grey45")
  for (i in seq_along(dl))
    lines(dl[[i]]$day, dl[[i]][[var]], lwd = 2.2, col = cols[(i - 1) %% length(cols) + 1])
  if (length(dl) > 1)
    legend("topright", names(dl), col = cols[seq_along(dl)], lwd = 2.2,
           bty = "n", cex = 0.78)
}

## -------------------------------------------------------------------- UI ----
ui <- fluidPage(
  titlePanel("남성 성선기능저하증 QSP 대시보드 — Male Hypogonadism QSP Dashboard"),
  tags$p(style = "color:#555; margin-top:-8px;",
         "HPG axis · Vermeulen binding equilibrium with regulated SHBG · ",
         "intratesticular testosterone · erythropoiesis · bone · body composition. ",
         tags$b("Research and education only — not for clinical use.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("patient", "환자 유형 (patient archetype)",
                  choices = PATIENTS, selected = "organic"),
      selectInput("reg1", "요법 A (regimen A)", choices = REGIMENS, selected = "gel"),
      selectInput("reg2", "요법 B (regimen B)", choices = REGIMENS, selected = "im_q2wk"),
      selectInput("reg3", "요법 C (regimen C)", choices = REGIMENS, selected = "none"),
      sliderInput("days", "시뮬레이션 기간 (days)", 90, 1095, 365, step = 30),
      hr(),
      h5("생리 파라미터 (physiology)"),
      sliderInput("age", "나이 (years)", 20, 85, 50, step = 1),
      sliderInput("shbg0", "기저 SHBG (nmol/L)", 8, 90, 35, step = 1),
      sliderInput("ins", "인슐린 저항성 (1 = lean)", 0.6, 4.0, 1.0, step = 0.1),
      sliderInput("fat0", "기저 지방량 (kg)", 10, 60, 26, step = 1),
      hr(),
      h5("모델 가정 (exposed assumptions)"),
      sliderInput("ec50h", "EC50_HEPC — 적혈구 반응 (pg/mL)",
                  100, 700, 300, step = 25),
      helpText(style = "font-size:11px;",
               "이 값이 생리적 범위보다 높아야 적혈구 반응이 치료 구간에서 ",
               "볼록해지고 첨두가 지배적이 된다. 보정된 값이며 측정값이 아니다. ",
               "(This single calibrated parameter is what makes waveform matter ",
               "for haematocrit — slide it down and the effect disappears.)"),
      sliderInput("aiflag", "아나스트로졸 병용 (1 mg/day, 0 = 없음)",
                  0, 1, 0, step = 1),
      hr(),
      actionButton("run", "시뮬레이션 실행 (Run)", class = "btn-primary",
                   style = "width:100%;")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("① 환자 프로파일",
          br(),
          fluidRow(column(6, verbatimTextOutput("profile")),
                   column(6, verbatimTextOutput("dx"))),
          hr(), plotOutput("overview", height = "620px")),
        tabPanel("② 결합 평형 (SHBG)",
          br(),
          tags$p("측정되는 총 테스토스테론과 작용하는 유리 테스토스테론 사이의 ",
                 "관계는 SHBG가 지배하는 비선형 질량작용 평형이다. SHBG 완충이 ",
                 "포화하기 때문에 유리 분율은 총 T가 올라갈수록 함께 올라간다 ",
                 "— 즉 FT는 TT의 볼록함수이다."),
          plotOutput("nomogram", height = "420px"),
          hr(), verbatimTextOutput("nomo_txt")),
        tabPanel("③ 진단 역치의 프레임 의존성",
          br(),
          tags$p("단일 총 T 역치(300 ng/dL)는 단일 유리 T 역치(65 pg/mL)가 아니다. ",
                 "SHBG에 따라 같은 역치가 237 ng/dL에서 635 ng/dL까지 움직인다."),
          plotOutput("framing", height = "400px"),
          hr(), tableOutput("frame_tbl")),
        tabPanel("④ PK 파형",
          br(),
          tags$p("같은 주간 용량이라도 제형마다 파형이 다르다. 유리 T 곡선에서 ",
                 "첨두가 총 T보다 더 크게 증폭되는 것이 볼록 결합의 직접적 결과다."),
          plotOutput("pk_tt", height = "300px"),
          plotOutput("pk_ft", height = "300px"),
          hr(), tableOutput("pk_tbl")),
        tabPanel("⑤ 적혈구증가증 (볼록성 원장)",
          br(),
          tags$p("제형 간 적혈구증가증 격차를 분해한다. (a) 볼록 결합은 실재하고 ",
                 "진폭에 비례하지만 작다. (b) 볼록 적혈구 반응은 합성되지 않는다 ",
                 "— 초생리적 첨두가 Hill의 오목한 상단으로 넘어가기 때문이며, ",
                 "이는 모델이 자기 가설을 반증한 지점이다. 용량을 맞추면 파형 효과는 ",
                 "Hct 약 0.1점에 불과하고, 문헌의 격차는 대부분 용량 차이와 ",
                 "고정 역치(Hct>54%)를 집단 분포에 적용할 때 생기는 증폭에서 나온다. ",
                 "즉 약물이 아니라 판정 규칙의 성질이다."),
          plotOutput("hct_plot", height = "340px"),
          hr(), tableOutput("convex_tbl"),
          hr(), h5("EC50_HEPC 민감도"), tableOutput("sens_tbl")),
        tabPanel("⑥ 고환내 T · 생식능",
          br(),
          tags$p("혈청 T 정상화와 고환내 T(ITT) 유지는 서로 다른 문제다. ",
                 "외인성 T는 LH를 끄고 ITT를 정상의 몇 퍼센트로 붕괴시킨다."),
          fluidRow(column(6, plotOutput("itt_plot", height = "320px")),
                   column(6, plotOutput("sperm_plot", height = "320px"))),
          hr(), tableOutput("itt_tbl")),
        tabPanel("⑦ 골 · 체성분",
          br(),
          tags$p("골과 지방에 대한 효과는 대부분 테스토스테론이 아니라 방향화 ",
                 "산물인 에스트라디올을 경유한다 (Finkelstein NEJM 2013). ",
                 "사이드바에서 아나스트로졸을 켜면 그 대가가 보인다."),
          fluidRow(column(6, plotOutput("bmd_plot", height = "300px")),
                   column(6, plotOutput("bc_plot", height = "300px"))),
          hr(), tableOutput("fink_tbl")),
        tabPanel("⑧ 바이오마커",
          br(),
          fluidRow(column(6, plotOutput("bm_lhfsh", height = "290px")),
                   column(6, plotOutput("bm_e2dht", height = "290px"))),
          fluidRow(column(6, plotOutput("bm_shbg", height = "290px")),
                   column(6, plotOutput("bm_psa", height = "290px")))),
        tabPanel("⑨ 시나리오 비교",
          br(), tableOutput("compare_tbl"),
          hr(), plotOutput("compare_plot", height = "420px")),
        tabPanel("⑩ 임상시험 대조표",
          br(),
          tags$p("모델 출력과 발표된 임상시험 종점의 대조. 마지막 항목은 모델이 ",
                 "재현하지 못하는 결과이며, 의도적으로 남겨두었다."),
          verbatimTextOutput("trial_txt"))
      )
    )
  )
)

## ---------------------------------------------------------------- server ----
server <- function(input, output, session) {

  sim <- eventReactive(input$run, {
    weeks <- ceiling(input$days / 7)
    base <- pt_pars(input$patient)
    base$AGE <- input$age
    base$SHBG0 <- input$shbg0
    base$INS_REL <- input$ins
    base$FAT0 <- input$fat0
    base$EC50_HEPC <- input$ec50h

    keys <- c(input$reg1, input$reg2, input$reg3)
    labs <- names(REGIMENS)[match(keys, unlist(REGIMENS))]
    keep <- !duplicated(keys)
    keys <- keys[keep]; labs <- labs[keep]

    out <- lapply(seq_along(keys), function(i) {
      e <- make_events(keys[i], weeks)
      if (input$aiflag > 0.5)
        e <- if (is.null(e)) anastrozole(1, weeks = weeks) else c(e, anastrozole(1, weeks = weeks))
      run_arm(base, e, days = input$days, label = labs[i])
    })
    names(out) <- labs
    list(arms = out, pars = base, days = input$days)
  }, ignoreNULL = FALSE)

  ## ---- tab 1 -----------------------------------------------------------
  output$profile <- renderPrint({
    s <- sim(); p <- s$pars
    cat("환자 프로파일 (patient profile)\n")
    cat("--------------------------------\n")
    cat(sprintf("age                 : %d years\n", as.integer(p$AGE)))
    cat(sprintf("baseline SHBG       : %.0f nmol/L\n", p$SHBG0))
    cat(sprintf("insulin resistance  : %.1f x lean referent\n", p$INS_REL))
    cat(sprintf("fat mass            : %.0f kg\n", p$FAT0))
    cat(sprintf("Leydig loss         : %.0f%%\n", 100 * (p$LEY_DMG %||% 0)))
    cat(sprintf("Sertoli loss        : %.0f%%\n", 100 * (p$SER_DMG %||% 0)))
    cat(sprintf("pituitary capacity  : %.2f\n", p$PITF %||% 1))
    cat(sprintf("EC50_HEPC (assumed) : %.0f pg/mL\n", p$EC50_HEPC))
  })

  output$dx <- renderPrint({
    s <- sim(); d <- s$arms[[1]]
    tt <- .at(d, 0, "TT_out"); ft <- .at(d, 0, "FT_out"); sh <- .at(d, 0, "SHBG")
    cat("치료 전 상태 (pre-treatment, day 0)\n")
    cat("--------------------------------\n")
    cat(sprintf("total T        : %6.0f ng/dL   (%s)\n", tt,
                if (tt < 300) "below the 300 cut-off" else "above the 300 cut-off"))
    cat(sprintf("free T         : %6.1f pg/mL   (%s)\n", ft,
                if (ft < 65) "below the 65 cut-off" else "above the 65 cut-off"))
    cat(sprintf("free fraction  : %6.2f %%\n", .at(d, 0, "FREEPCT")))
    cat(sprintf("SHBG           : %6.1f nmol/L\n", sh))
    cat(sprintf("LH / FSH       : %5.2f / %5.2f IU/L\n",
                .at(d, 0, "LH"), .at(d, 0, "FSH")))
    cat(sprintf("ITT            : %6.1f %% of normal\n", .at(d, 0, "ITT_PCT")))
    cat(sprintf("haematocrit    : %6.1f %%\n", .at(d, 0, "HCT_out")))
    if ((tt < 300) != (ft < 65))
      cat("\n>> 총 T와 유리 T의 판정이 엇갈립니다 — SHBG가 원인입니다.\n",
          ">> TT and FT disagree here. That is the binding protein, not the gonad.\n",
          sep = "")
  })

  output$overview <- renderPlot({
    s <- sim(); dl <- s$arms
    op <- par(mfrow = c(3, 3), mar = c(4, 4.2, 2.4, 1)); on.exit(par(op))
    lineplot(dl, "TT_out", "ng/dL", "Total testosterone", c(300, 1000))
    lineplot(dl, "FT_out", "pg/mL", "Free testosterone", 65)
    lineplot(dl, "SHBG", "nmol/L", "SHBG")
    lineplot(dl, "LH", "IU/L", "LH")
    lineplot(dl, "ITT_PCT", "% of normal", "Intratesticular T", 25)
    lineplot(dl, "E2", "pg/mL", "Estradiol")
    lineplot(dl, "HCT_out", "%", "Haematocrit", 54)
    lineplot(dl, "SPERM_out", "million/mL", "Sperm concentration", 16)
    lineplot(dl, "BMD_TRAB", "% change", "Trabecular vBMD")
  })

  ## ---- tab 2: nomogram --------------------------------------------------
  ftcalc <- function(tt, s) {
    KS <- 1e9; KA <- 3.6e4; MW <- 288.4
    Tm <- tt * 1e-8 / MW; Sm <- s * 1e-9
    N <- 1 + KA * (4.3 * 10 / 66500); a <- N * KS; b <- N + KS * Sm - KS * Tm
    (-b + sqrt(b^2 + 4 * a * Tm)) / (2 * a) * MW * 1e9
  }

  output$nomogram <- renderPlot({
    tt <- seq(50, 2000, 10)
    shbgs <- c(15, 25, 35, 55, 80)
    op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
    plot(NA, xlim = range(tt), ylim = c(0, max(ftcalc(2000, 15))),
         xlab = "total testosterone (ng/dL)", ylab = "free testosterone (pg/mL)",
         main = "Free T vs total T — convex, and SHBG sets the slope", las = 1)
    grid(col = "grey90")
    for (i in seq_along(shbgs))
      lines(tt, ftcalc(tt, shbgs[i]), lwd = 2.4, col = PALETTE[i])
    abline(h = 65, lty = 2); abline(v = 300, lty = 2)
    legend("topleft", sprintf("SHBG %d", shbgs), col = PALETTE[seq_along(shbgs)],
           lwd = 2.4, bty = "n", cex = 0.85)
    ff <- sapply(shbgs, function(s) 100 * ftcalc(tt, s) / (tt * 0.03467 * 288.4))
    plot(NA, xlim = range(tt), ylim = range(ff), xlab = "total testosterone (ng/dL)",
         ylab = "free fraction (%)",
         main = "Free FRACTION rises with total T\n(the SHBG buffer saturates)", las = 1)
    grid(col = "grey90")
    for (i in seq_along(shbgs)) lines(tt, ff[, i], lwd = 2.4, col = PALETTE[i])
  })

  output$nomo_txt <- renderPrint({
    d2 <- diff(diff(sapply(seq(100, 2000, 50), ftcalc, s = 35)))
    cat("At SHBG 35 nmol/L, albumin 4.3 g/dL:\n\n")
    for (tt in c(100, 300, 600, 1000, 1500, 2000))
      cat(sprintf("  TT %5d ng/dL  ->  FT %6.1f pg/mL   free fraction %.2f %%\n",
                  tt, ftcalc(tt, 35), 100 * ftcalc(tt, 35) / (tt * 0.03467 * 288.4)))
    cat(sprintf("\nConvexity: every second difference positive? %s (min %.4f)\n",
                all(d2 > 0), min(d2)))
    cat("\nConsequence (Jensen): a regimen with large peak-to-trough swings delivers\n")
    cat("MORE time-averaged free T than a flat one at the same mean total T.\n")
    cat("Quantified on tab 5 — the effect is real but modest at the binding step.\n")
  })

  ## ---- tab 3: diagnostic framing ---------------------------------------
  output$framing <- renderPlot({
    shbgs <- seq(8, 95, 1)
    tteq <- sapply(shbgs, function(s)
      uniroot(function(x) ftcalc(x, s) - 65, c(20, 4000))$root)
    plot(shbgs, tteq, type = "l", lwd = 3, col = "#c0392b",
         xlab = "SHBG (nmol/L)",
         ylab = "total T equivalent to free T = 65 pg/mL (ng/dL)",
         main = "The same free-T threshold is a moving total-T threshold", las = 1)
    grid(col = "grey90")
    abline(h = 300, lty = 2, lwd = 2, col = "#2c6fb5")
    text(70, 315, "fixed TT cut-off = 300 ng/dL", col = "#2c6fb5", cex = 0.9)
    polygon(c(shbgs, rev(shbgs)), c(tteq, rep(300, length(shbgs))),
            col = adjustcolor("#c0392b", 0.10), border = NA)
    x0 <- shbgs[which.min(abs(tteq - 300))]
    abline(v = x0, lty = 3)
    text(x0, max(tteq) * 0.95, sprintf("agreement only at SHBG %.0f", x0),
         pos = 4, cex = 0.9)
    text(20, 200, "TT cut-off\nOVER-diagnoses", cex = 0.95, col = "#7f8c8d")
    text(80, 430, "TT cut-off\nUNDER-diagnoses", cex = 0.95, col = "#7f8c8d")
  })

  output$frame_tbl <- renderTable({
    shbgs <- c(15, 20, 25, 30, 35, 45, 55, 70, 90)
    tteq <- sapply(shbgs, function(s)
      uniroot(function(x) ftcalc(x, s) - 65, c(20, 4000))$root)
    data.frame(
      `SHBG (nmol/L)` = shbgs,
      `Free T at TT=300 (pg/mL)` = round(ftcalc(300, shbgs), 1),
      `TT giving FT=65 (ng/dL)` = round(tteq, 0),
      Verdict = ifelse(tteq < 300, "TT cut-off over-diagnoses",
                       "TT cut-off under-diagnoses"),
      check.names = FALSE)
  })

  ## ---- tab 4: PK waveforms ---------------------------------------------
  output$pk_tt <- renderPlot({
    s <- sim()
    lineplot(s$arms, "TT_out", "ng/dL",
             "Total testosterone — waveform by formulation", c(300, 1000))
  })
  output$pk_ft <- renderPlot({
    s <- sim()
    lineplot(s$arms, "FT_out", "pg/mL",
             "Free testosterone — peaks are amplified by convex binding", 65)
  })
  output$pk_tbl <- renderTable({
    s <- sim(); dd <- s$days; w1 <- max(0, dd - 90)
    do.call(rbind, lapply(names(s$arms), function(nm) {
      d <- s$arms[[nm]]
      data.frame(Regimen = nm,
        `Cavg TT` = round(.mean_over(d, "TT_out", w1, dd), 0),
        `Cmax TT` = round(max(d$TT_out[d$day >= w1]), 0),
        `Cmin TT` = round(min(d$TT_out[d$day >= w1]), 0),
        `peak:trough` = round(max(d$TT_out[d$day >= w1]) /
                                max(min(d$TT_out[d$day >= w1]), 1), 2),
        `Cavg FT` = round(.mean_over(d, "FT_out", w1, dd), 1),
        `SHBG` = round(.at(d, dd, "SHBG"), 1),
        check.names = FALSE)
    }))
  })

  ## ---- tab 5: erythrocytosis -------------------------------------------
  output$hct_plot <- renderPlot({
    s <- sim()
    op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1)); on.exit(par(op))
    lineplot(s$arms, "HCT_out", "%", "Haematocrit", 54)
    lineplot(s$arms, "PSTOP", "% of population",
             "Population fraction over the 54% stopping rule")
  })

  output$convex_tbl <- renderTable({
    s <- sim(); dd <- s$days; w1 <- max(0, dd - 90)
    ec <- s$pars$EC50_HEPC; hh <- as.numeric(param(mod)$HILL_H)
    tab <- do.call(rbind, lapply(names(s$arms), function(nm) {
      d <- s$arms[[nm]]
      mTT <- .mean_over(d, "TT_out", w1, dd)
      mFT <- .mean_over(d, "FT_out", w1, dd)
      mS  <- .mean_over(d, "SHBG", w1, dd)
      mDR <- .mean_over(d, "DRIVE_out", w1, dd)
      dofm <- mFT^hh / (ec^hh + mFT^hh)
      data.frame(Regimen = nm,
        `(a) convex binding %` = round(100 * (mFT / ftcalc(mTT, mS) - 1), 1),
        `(b) convex response %` = round(100 * (mDR / max(dofm, 1e-9) - 1), 1),
        `mean drive` = round(mDR, 4),
        `dHct` = round(.delta(d, dd, "HCT_out"), 2),
        `(c) % over 54` = round(.at(d, dd, "PSTOP"), 2),
        check.names = FALSE)
    }))
    tab
  })

  output$sens_tbl <- renderTable({
    s <- sim()
    base <- s$pars
    wk <- ceiling(s$days / 7)
    do.call(rbind, lapply(c(150, 225, 300, 450, 600), function(ec) {
      p <- base; p$EC50_HEPC <- ec
      g <- run_arm(p, tst_gel(81, weeks = wk), s$days)
      i <- run_arm(p, tst_im(200, 14, weeks = wk), s$days)
      data.frame(`EC50_HEPC (pg/mL)` = ec,
                 `dHct gel` = round(.delta(g, s$days, "HCT_out"), 2),
                 `dHct IM q2wk` = round(.delta(i, s$days, "HCT_out"), 2),
                 `gap` = round(.at(i, s$days, "HCT_out") - .at(g, s$days, "HCT_out"), 2),
                 `incidence ratio` = round(.at(i, s$days, "PSTOP") /
                                             max(.at(g, s$days, "PSTOP"), 1e-6), 2),
                 check.names = FALSE)
    }))
  })

  ## ---- tab 6: ITT / fertility ------------------------------------------
  output$itt_plot <- renderPlot({
    s <- sim()
    lineplot(s$arms, "ITT_PCT", "% of normal",
             "Intratesticular testosterone", 25)
  })
  output$sperm_plot <- renderPlot({
    s <- sim()
    lineplot(s$arms, "SPERM_out", "million/mL", "Sperm concentration", 16)
  })
  output$itt_tbl <- renderTable({
    s <- sim(); dd <- s$days
    do.call(rbind, lapply(names(s$arms), function(nm) {
      d <- s$arms[[nm]]
      data.frame(Regimen = nm,
        `serum TT (ng/dL)` = round(.mean_over(d, "TT_out", max(0, dd - 30), dd), 0),
        `LH (IU/L)` = round(.at(d, dd, "LH"), 2),
        `ITT (% normal)` = round(.at(d, dd, "ITT_PCT"), 1),
        `sperm (M/mL)` = round(.at(d, dd, "SPERM_out"), 1),
        `fertility` = ifelse(.at(d, dd, "SPERM_out") >= 16, "preserved",
                      ifelse(.at(d, dd, "SPERM_out") >= 1, "impaired", "azoospermic")),
        check.names = FALSE)
    }))
  })

  ## ---- tab 7: bone / body composition ----------------------------------
  output$bmd_plot <- renderPlot({
    s <- sim()
    lineplot(s$arms, "BMD_TRAB", "% change from start", "Trabecular vBMD", 0)
  })
  output$bc_plot <- renderPlot({
    s <- sim()
    op <- par(mfrow = c(1, 2), mar = c(4, 4.2, 3, 0.6)); on.exit(par(op))
    lineplot(s$arms, "LEAN", "kg", "Lean mass")
    lineplot(s$arms, "FAT", "kg", "Fat mass")
  })
  output$fink_tbl <- renderTable({
    ## a compact Finkelstein-style dissociation at fixed dose
    base <- list(AGE = 30, PITF = 0.02, FAT0 = 22, LEAN0 = 62)
    days <- 112; wk <- ceiling(days / 7)
    rows <- list()
    for (g in c(1.25, 2.5, 5, 10)) for (a in c(FALSE, TRUE)) {
      e <- tst_gel(g * 10, weeks = wk)
      if (a) e <- c(e, anastrozole(1, weeks = wk))
      d <- run_arm(base, e, days)
      rows[[length(rows) + 1]] <- data.frame(
        `gel (g/day, 1%)` = g, `anastrozole` = ifelse(a, "yes", "no"),
        `TT (ng/dL)` = round(.at(d, days, "TT_out"), 0),
        `E2 (pg/mL)` = round(.at(d, days, "E2"), 1),
        `d lean (kg)` = round(.delta(d, days, "LEAN"), 2),
        `d fat (kg)` = round(.delta(d, days, "FAT"), 2),
        check.names = FALSE)
    }
    do.call(rbind, rows)
  })

  ## ---- tab 8: biomarkers ------------------------------------------------
  output$bm_lhfsh <- renderPlot({
    s <- sim()
    op <- par(mfrow = c(1, 2), mar = c(4, 4.2, 3, 0.6)); on.exit(par(op))
    lineplot(s$arms, "LH", "IU/L", "LH")
    lineplot(s$arms, "FSH", "IU/L", "FSH")
  })
  output$bm_e2dht <- renderPlot({
    s <- sim()
    op <- par(mfrow = c(1, 2), mar = c(4, 4.2, 3, 0.6)); on.exit(par(op))
    lineplot(s$arms, "E2", "pg/mL", "Estradiol")
    lineplot(s$arms, "DHT", "ng/dL", "DHT")
  })
  output$bm_shbg <- renderPlot({
    s <- sim()
    op <- par(mfrow = c(1, 2), mar = c(4, 4.2, 3, 0.6)); on.exit(par(op))
    lineplot(s$arms, "SHBG", "nmol/L", "SHBG")
    lineplot(s$arms, "FREEPCT", "%", "Free fraction")
  })
  output$bm_psa <- renderPlot({
    s <- sim()
    op <- par(mfrow = c(1, 2), mar = c(4, 4.2, 3, 0.6)); on.exit(par(op))
    lineplot(s$arms, "PSA", "ng/mL", "PSA (saturation model)")
    lineplot(s$arms, "INHB", "pg/mL", "Inhibin B")
  })

  ## ---- tab 9: comparison ------------------------------------------------
  output$compare_tbl <- renderTable({
    s <- sim(); dd <- s$days; w1 <- max(0, dd - 90)
    do.call(rbind, lapply(names(s$arms), function(nm) {
      d <- s$arms[[nm]]
      data.frame(Regimen = nm,
        `Cavg TT` = round(.mean_over(d, "TT_out", w1, dd), 0),
        `Cavg FT` = round(.mean_over(d, "FT_out", w1, dd), 1),
        `SHBG` = round(.at(d, dd, "SHBG"), 1),
        `E2` = round(.at(d, dd, "E2"), 1),
        `dHct` = round(.delta(d, dd, "HCT_out"), 2),
        `dPSA` = round(.delta(d, dd, "PSA"), 2),
        `d lean (kg)` = round(.delta(d, dd, "LEAN"), 2),
        `d fat (kg)` = round(.delta(d, dd, "FAT"), 2),
        `vBMD (%)` = round(.delta(d, dd, "BMD_TRAB"), 2),
        `ITT (%)` = round(.at(d, dd, "ITT_PCT"), 1),
        `libido` = round(.at(d, dd, "LIBIDO"), 2),
        check.names = FALSE)
    }))
  })
  output$compare_plot <- renderPlot({
    s <- sim()
    op <- par(mfrow = c(2, 3), mar = c(4, 4.2, 2.6, 1)); on.exit(par(op))
    lineplot(s$arms, "TT_out", "ng/dL", "Total T", c(300, 1000))
    lineplot(s$arms, "FT_out", "pg/mL", "Free T", 65)
    lineplot(s$arms, "HCT_out", "%", "Haematocrit", 54)
    lineplot(s$arms, "BMD_TRAB", "%", "Trabecular vBMD", 0)
    lineplot(s$arms, "LEAN", "kg", "Lean mass")
    lineplot(s$arms, "LIBIDO", "0-10", "Sexual desire")
  })

  ## ---- tab 10: trial ledger ---------------------------------------------
  output$trial_txt <- renderPrint({
    MHG_trial_ledger()
    invisible(NULL)
  })
}

shinyApp(ui, server)
