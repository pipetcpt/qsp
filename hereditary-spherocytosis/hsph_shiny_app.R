## =====================================================================
##  Hereditary Spherocytosis QSP — Shiny dashboard
##  유전성 구상적혈구증 QSP 대시보드 (14 tabs)
##
##  Front end for hsph_mrgsolve_model.R.  Like the model file, this has NOT
##  been executed (no R toolchain in the build environment); it mirrors the
##  Python reference implementation, which has.
##
##  run:  shiny::runApp("hsph_shiny_app.R")
## =====================================================================
library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

MOD <- mread_cache("hsph_mrgsolve_model", ".")

## ---------------------------------------------------------------- helpers
## the geometry kernel, in R, so the geometry tab does not need the ODEs
vsph  <- function(A) A^1.5 / (6 * sqrt(pi))
dcrit <- function(A, V) {
  s  <- pmin(pmax(V / vsph(A), 1e-9), 1)
  2 * sqrt(A / pi) * cos(acos(-s) / 3 - 2 * pi / 3)
}
ofpoint <- function(A, V, mchc) {          # 50% lysis, %NaCl
  vs <- 0.75 * mchc * V / 100 + 0.04 * V
  vm <- vsph(A)
  ifelse(vm > vs, 290 * (V - vs) / (vm - vs) * (0.9 / 308), 0.9)
}

GENO <- list(
  "정상 (normal)"                              = c(fdef = 0.00, f_b3ves = 1.00),
  "보인자 (HS trait)"                          = c(fdef = 0.12, f_b3ves = 1.00),
  "경증 (mild HS)"                             = c(fdef = 0.20, f_b3ves = 1.00),
  "중등증 ANK1/SPTB (moderate)"                = c(fdef = 0.30, f_b3ves = 1.00),
  "중등증 SLC4A1 band 3 (moderate)"            = c(fdef = 0.30, f_b3ves = 0.15),
  "중증 SPTA1 (severe, recessive)"             = c(fdef = 0.45, f_b3ves = 1.00)
)

run_one <- function(par, end = 1400, delta = 2, burn = 1400) {
  m <- MOD %>% param(par)
  m %>% mrgsim(end = end, delta = delta, hmax = 2) %>% as_tibble()
}

## run to steady state, then apply an intervention at t = 0 of a second run
run_intervention <- function(base_par, iv_par, pre = 1400, post = 1460,
                             delta = 2) {
  s0 <- MOD %>% param(base_par) %>% mrgsim(end = pre, delta = pre) %>%
    as_tibble() %>% tail(1)
  ini <- MOD %>% param(base_par) %>%
    mrgsim(end = pre, delta = pre, obsonly = FALSE)
  MOD %>% param(c(base_par, iv_par)) %>%
    init(as.list(ini@data[nrow(ini@data), names(init(MOD))])) %>%
    mrgsim(end = post, delta = delta, hmax = 2) %>% as_tibble()
}

thm <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        legend.position = "bottom")

## ------------------------------------------------------------------- UI
ui <- fluidPage(
  titlePanel("유전성 구상적혈구증 QSP 모델 · Hereditary Spherocytosis"),
  tags$p(tags$em(paste(
    "적혈구는 막 면적 S와 부피 V, 두 숫자다. 이 두 숫자가 최소 원통 직경",
    "D_c를 정하고, 비장 정맥동 벽이 볼 수 있는 것은 D_c뿐이다.",
    "면적 손실 → D_c 상승 → 비장 삭(cord) 체류 연장 → 면적 손실 —",
    "이 되먹임 고리가 이 질환이다."))),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("geno", "유전자형 / 중증도 (genotype)",
                  names(GENO), selected = "중등증 ANK1/SPTB (moderate)"),
      sliderInput("fdef", "수직 연결 결손 fdef = 1 − 스펙트린 함량",
                  0, 0.6, 0.30, 0.01),
      sliderInput("fb3", "소포의 band 3 함량 f_b3ves",
                  0, 1, 1.00, 0.05,
                  post = "  (1 = 스펙트린/안키린형, 0.15 = band 3형)"),
      hr(),
      radioButtons("spleen", "비장 (spleen)",
                   c("정상 (intact)" = "1",
                     "부분절제 20% 잔존" = "0.2",
                     "부분절제 10% 잔존" = "0.1",
                     "전절제 (total splenectomy)" = "0"), selected = "1"),
      checkboxInput("regrow", "부분절제 후 잔존 비장 재성장", TRUE),
      hr(),
      sliderInput("mita", "미타피밧 mitapivat (mg BID)", 0, 100, 0, 5),
      selectInput("ugt", "UGT1A1",
                  c("정상 (wild type)" = "1",
                    "길버트 증후군 UGT1A1*28/*28" = "0.287")),
      checkboxInput("folate", "엽산 보충 (folate 1–5 mg/day)", TRUE),
      hr(),
      checkboxInput("parvo", "파르보바이러스 B19 재생불량 위기", FALSE),
      numericInput("parvo_t", "  B19 발생 시점 (day)", 30, 0, 300),
      hr(),
      checkboxInput("tx", "정기 수혈 (2 units q4w)", FALSE),
      checkboxInput("chel", "철 킬레이션 (deferasirox)", FALSE),
      hr(),
      sliderInput("tend", "시뮬레이션 기간 (days)", 200, 3650, 730, 30),
      actionButton("go", "실행 (run)", class = "btn-primary")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · 환자 프로파일",
                 h4("정상 상태 요약 (steady state)"),
                 tableOutput("profile"),
                 h4("중증도 스펙트럼 — fdef 하나로 결정된다"),
                 plotOutput("spectrum", height = 340)),
        tabPanel("2 · 적혈구 기하학",
                 h4("D_c = 2√(S/π)·cos(arccos(−V/V_sph)/3 − 2π/3)"),
                 tags$p("적합(fitting)된 것이 하나도 없는 닫힌 식이다."),
                 plotOutput("geomap", height = 380),
                 h4("세포 나이에 따른 (S, V) 궤적"),
                 plotOutput("geotraj", height = 300)),
        tabPanel("3 · 비장 필터",
                 h4("D_c 하나만 보는 필터"),
                 plotOutput("filter", height = 380),
                 tableOutput("filtertab")),
        tabPanel("4 · 혈액검사 (CBC)",
                 plotOutput("cbc", height = 620)),
        tabPanel("5 · 용혈 지표",
                 plotOutput("hemol", height = 560)),
        tabPanel("6 · 청소 기전 분해",
                 tags$p(paste("정상 적혈구는 표지(labelling)되어 죽고, HS",
                              "적혈구는 모양 때문에 죽는다. 같은 장기,",
                              "다른 기전.")),
                 plotOutput("clearance", height = 420),
                 h4("막 손실이 어디서 일어나는가"),
                 tableOutput("stripping")),
        tabPanel("7 · 비장절제",
                 h4("빈혈은 고치고 세포는 고치지 못한다"),
                 plotOutput("splx", height = 430),
                 tableOutput("splxtab")),
        tabPanel("8 · 유전자형 × 비장절제",
                 tags$p(paste("이 모델에서 두 유전자형을 가르는 파라미터는",
                              "하나다: band 3가 소포에 실려 나가는가.")),
                 plotOutput("genoplot", height = 400),
                 tableOutput("genotab")),
        tabPanel("9 · 진단검사",
                 h4("삼투취약성 곡선 (osmotic fragility)"),
                 plotOutput("offcurve", height = 340),
                 h4("EMA 결합 · 삼투구배 ektacytometry"),
                 plotOutput("emaplot", height = 300)),
        tabPanel("10 · 위기 (crises)",
                 h4("파르보바이러스 B19: ΔHb ≈ Hb × 정지일수 / 수명"),
                 plotOutput("crisis", height = 420),
                 tableOutput("crisistab")),
        tabPanel("11 · 약물 (mitapivat)",
                 plotOutput("drug", height = 560),
                 tableOutput("drugtab")),
        tabPanel("12 · 담석 · 철",
                 plotOutput("stone", height = 480),
                 tableOutput("stonetab")),
        tabPanel("13 · 시나리오 비교",
                 plotOutput("scen", height = 620),
                 tableOutput("scentab")),
        tabPanel("14 · 가상 환자군",
                 numericInput("npop", "환자 수", 100, 20, 500, 20),
                 plotOutput("pop", height = 480),
                 tableOutput("poptab"))
      )
    )
  )
)

## --------------------------------------------------------------- server
server <- function(input, output, session) {

  observeEvent(input$geno, {
    g <- GENO[[input$geno]]
    updateSliderInput(session, "fdef", value = unname(g["fdef"]))
    updateSliderInput(session, "fb3",  value = unname(g["f_b3ves"]))
  })

  par_now <- reactive({
    sf <- as.numeric(input$spleen)
    c(fdef      = input$fdef,
      f_b3ves   = input$fb3,
      spl_frac  = sf,
      k_regrow  = if (input$regrow && sf > 0 && sf < 1) 1 / 700 else 0,
      dose_m    = input$mita,
      ugt_f     = as.numeric(input$ugt),
      fol_ok    = if (input$folate) 1.0 else 0.35,
      parvo_t   = if (input$parvo) input$parvo_t else -1,
      tx_start  = if (input$tx) 0 else -1,
      k_chel    = if (input$chel) 0.0016 else 0)
  })

  sim <- eventReactive(input$go, {
    run_one(par_now(), end = input$tend)
  }, ignoreNULL = FALSE)

  ss <- reactive(tail(sim(), 1))

  ## ---- 1 profile -------------------------------------------------------
  output$profile <- renderTable({
    s <- ss()
    tibble(
      지표 = c("혈색소 Hb (g/dL)", "적혈구 RBC (10¹²/L)", "MCV (fL)",
               "MCHC (g/dL)", "망상적혈구 (%)", "적혈구 수명 (일)",
               "총빌리루빈 (mg/dL)", "비장 부피 (mL)",
               "막 면적 S (µm²)", "최소 원통 직경 D_c (µm)",
               "EMA 결합비", "50% 용혈 (%NaCl)",
               "적혈구 결합 IgG (분자/세포)"),
      값 = c(s$Hb, s$RBC, s$MCV, s$MCHC, s$RETpct, s$LIFE, s$TBIL,
             s$SPLEEN, s$AREA, s$DCRIT, s$EMA, s$MCF, s$IGGC)
    )
  }, digits = 3)

  output$spectrum <- renderPlot({
    fs <- seq(0, 0.5, by = 0.05)
    d <- lapply(fs, function(f) {
      p <- par_now(); p["fdef"] <- f
      tail(run_one(p, end = 1200, delta = 1200), 1) %>% mutate(fdef = f)
    }) %>% bind_rows()
    d %>% select(fdef, Hb, RETpct, LIFE, TBIL, AREA, DCRIT) %>%
      pivot_longer(-fdef) %>%
      ggplot(aes(fdef, value)) + geom_line(linewidth = 1, colour = "#c62828") +
      geom_point() + facet_wrap(~name, scales = "free_y") +
      labs(x = "fdef = 1 − 스펙트린 함량", y = NULL) + thm
  })

  ## ---- 2 geometry ------------------------------------------------------
  output$geomap <- renderPlot({
    g <- expand.grid(S = seq(90, 150, 1), V = seq(70, 105, 1)) %>%
      mutate(Dc = dcrit(S, V)) %>% filter(V < vsph(S))
    pts <- tibble(
      lab = c("정상 young", "정상 old", "경증 HS", "중등증 HS", "중증 HS"),
      S = c(140, 125.3, 128, 118, 106),
      V = c(94, 86.1, 89, 86, 80))
    ggplot(g, aes(S, V)) +
      geom_raster(aes(fill = Dc)) +
      geom_contour(aes(z = Dc), colour = "white", breaks = seq(2.6, 4.2, .2)) +
      scale_fill_viridis_c(option = "magma", direction = -1) +
      geom_point(data = pts, size = 3, colour = "cyan") +
      geom_text(data = pts, aes(label = lab), nudge_y = 2.4,
                colour = "cyan", size = 3.4) +
      labs(x = "막 면적 S (µm²)", y = "부피 V (µm³)",
           fill = "D_c (µm)",
           title = "면적을 잃되 부피를 잃지 않으면 D_c는 급히 오른다") + thm
  })

  output$geotraj <- renderPlot({
    s <- sim()
    s %>% select(time, AREA, VOLc, DCRIT, SPHER) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (일)", y = NULL) + thm
  })

  ## ---- 3 splenic filter ------------------------------------------------
  output$filter <- renderPlot({
    p  <- as.list(param(MOD))
    Dc <- seq(2.5, 4.2, 0.01)
    ps <- 1 / (1 + exp(-(Dc - p$D50) / p$wD))
    tc <- p$tau0 * exp((Dc - p$Dc_ref) / p$w_esc)
    Rc <- pmin(p$R_MAX, p$f_pass0 * ps * tc)
    tibble(Dc,
           `p_slow (삭으로 가는 비율)` = ps,
           `τ_cord (일)` = tc,
           `R_cord (삭 체류 시간 분율)` = Rc,
           `막 손실 배수 (1+41·R)` = 1 + (p$cordamp - 1) * Rc) %>%
      pivot_longer(-Dc) %>%
      ggplot(aes(Dc, value)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = c(2.88, 3.27), linetype = 2,
                 colour = c("#2e7d32", "#c62828")) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "최소 원통 직경 D_c (µm)", y = NULL,
           caption = "녹색 = 정상 적혈구, 빨강 = 중등증 HS") + thm
  })

  output$filtertab <- renderTable({
    s <- ss()
    tibble(`삭 체류 분율 R_cord` = s$RCORD,
           `비장에서 일어나는 막 손실 (%)` =
             100 * 41 * s$RCORD / (1 + 41 * s$RCORD),
           `비장 부피 (mL)` = s$SPLEEN)
  }, digits = 4)

  ## ---- 4 CBC -----------------------------------------------------------
  output$cbc <- renderPlot({
    sim() %>% select(time, Hb, Hct, RBC, MCV, MCHC, RETpct, RETabs, LIFE) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#1565c0") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "시간 (일)", y = NULL) + thm
  })

  ## ---- 5 haemolysis ----------------------------------------------------
  output$hemol <- renderPlot({
    sim() %>% select(time, TBIL, BILU, BILC, LDH, HPT, BRPROD, SPLEEN) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#6a1b9a") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "시간 (일)", y = NULL) + thm
  })

  ## ---- 6 clearance decomposition ---------------------------------------
  output$clearance <- renderPlot({
    sim() %>% select(time, FGEOM, FOPS, FSEN, FLIV, FLYS) %>%
      pivot_longer(-time, names_to = "mech") %>%
      mutate(mech = recode(mech,
                           FGEOM = "기하학적 (비장 삭)", FOPS = "옵소닌 (IgG)",
                           FSEN = "노화 표지", FLIV = "간 (기하학적)",
                           FLYS = "혈관내 용혈")) %>%
      ggplot(aes(time, value, fill = mech)) +
      geom_area() + scale_fill_brewer(palette = "Set2") +
      labs(x = "시간 (일)", y = "파괴의 분율", fill = NULL) + thm
  })

  output$stripping <- renderTable({
    lapply(names(GENO), function(g) {
      p <- par_now(); p["fdef"] <- GENO[[g]]["fdef"]
      p["f_b3ves"] <- GENO[[g]]["f_b3ves"]
      s <- tail(run_one(p, end = 1200, delta = 1200), 1)
      tibble(유전자형 = g, R_cord = s$RCORD,
             `비장에서의 막 손실 (%)` = 100 * 41 * s$RCORD /
               (1 + 41 * s$RCORD),
             `기하학적 청소` = s$FGEOM, `옵소닌 청소` = s$FOPS,
             `노화 청소` = s$FSEN)
    }) %>% bind_rows()
  }, digits = 3)

  ## ---- 7 splenectomy ---------------------------------------------------
  output$splx <- renderPlot({
    base <- par_now(); base["spl_frac"] <- 1
    arms <- list(`수술 없음` = c(spl_frac = 1),
                 `부분절제 20%` = c(spl_frac = 0.2, k_regrow = 1 / 700),
                 `부분절제 10%` = c(spl_frac = 0.1, k_regrow = 1 / 700),
                 `전절제` = c(spl_frac = 0))
    lapply(names(arms), function(a)
      run_intervention(base, arms[[a]], post = 1825, delta = 5) %>%
        mutate(arm = a)) %>% bind_rows() %>%
      select(time, arm, Hb, RETpct, TBIL, SPLEEN, AREA, LIFE) %>%
      pivot_longer(c(-time, -arm)) %>%
      ggplot(aes(time / 365, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "수술 후 경과 (년)", y = NULL, colour = NULL) + thm
  })

  output$splxtab <- renderTable({
    base <- par_now(); base["spl_frac"] <- 1
    pre  <- tail(run_one(base, end = 1400, delta = 1400), 1)
    post <- tail(run_intervention(base, c(spl_frac = 0), post = 1400,
                                  delta = 1400), 1)
    tibble(지표 = c("Hb", "망상적혈구 %", "수명 (일)", "총빌리루빈",
                    "막 면적 S", "EMA 결합비", "50% 용혈 %NaCl", "MCHC"),
           수술전 = c(pre$Hb, pre$RETpct, pre$LIFE, pre$TBIL, pre$AREA,
                      pre$EMA, pre$MCF, pre$MCHC),
           수술후 = c(post$Hb, post$RETpct, post$LIFE, post$TBIL, post$AREA,
                      post$EMA, post$MCF, post$MCHC))
  }, digits = 3)

  ## ---- 8 genotype x splenectomy ----------------------------------------
  output$genoplot <- renderPlot({
    lapply(c(1.0, 0.15), function(fb) {
      p <- par_now(); p["f_b3ves"] <- fb; p["spl_frac"] <- 1
      lapply(c(1, 0), function(sf) {
        q <- p; q["spl_frac"] <- sf
        tail(run_one(q, end = 1400, delta = 1400), 1) %>%
          mutate(genotype = ifelse(fb > 0.5, "스펙트린/안키린", "band 3"),
                 spleen = ifelse(sf > 0.5, "비장 있음", "비장절제"))
      }) %>% bind_rows()
    }) %>% bind_rows() %>%
      select(genotype, spleen, Hb, LIFE, IGGC, FOPS, FGEOM) %>%
      pivot_longer(c(-genotype, -spleen)) %>%
      ggplot(aes(spleen, value, fill = genotype)) +
      geom_col(position = "dodge") +
      facet_wrap(~name, scales = "free_y") +
      scale_fill_manual(values = c("#1565c0", "#c62828")) +
      labs(x = NULL, y = NULL, fill = NULL) + thm
  })

  output$genotab <- renderTable({
    lapply(c(1.0, 0.15), function(fb) {
      p <- par_now(); p["f_b3ves"] <- fb
      p["spl_frac"] <- 1; a <- tail(run_one(p, end = 1400, delta = 1400), 1)
      p["spl_frac"] <- 0; b <- tail(run_one(p, end = 1400, delta = 1400), 1)
      tibble(유전자형 = ifelse(fb > .5, "스펙트린/안키린", "band 3"),
             `IgG/세포 (비장 있음)` = a$IGGC,
             `IgG/세포 (절제 후)` = b$IGGC,
             `Hb 전` = a$Hb, `Hb 후` = b$Hb, `ΔHb` = b$Hb - a$Hb)
    }) %>% bind_rows()
  }, digits = 2)

  ## ---- 9 diagnostics ---------------------------------------------------
  output$offcurve <- renderPlot({
    s <- ss()
    nacl <- seq(0.1, 0.85, 0.005)
    cells <- tibble(
      grp = c("정상", "환자"),
      S = c(132, s$AREA), V = c(89, s$VOLc),
      MCHC = c(33.7, s$MCHC))
    d <- lapply(seq_len(nrow(cells)), function(i) {
      m <- ofpoint(cells$S[i], cells$V[i], cells$MCHC[i])
      tibble(grp = cells$grp[i], nacl,
             lysis = 100 / (1 + exp((nacl - m) / 0.045)))
    }) %>% bind_rows()
    ggplot(d, aes(nacl, lysis, colour = grp)) + geom_line(linewidth = 1.1) +
      scale_x_reverse() +
      geom_vline(xintercept = c(0.43, 0.60), linetype = 3) +
      labs(x = "NaCl (%)", y = "용혈 (%)", colour = NULL,
           title = "곡선의 위치는 V_sph = S^1.5/(6√π) 하나로 결정된다") + thm
  })

  output$emaplot <- renderPlot({
    lapply(names(GENO), function(g) {
      p <- par_now(); p["fdef"] <- GENO[[g]]["fdef"]
      p["f_b3ves"] <- GENO[[g]]["f_b3ves"]
      tail(run_one(p, end = 1200, delta = 1200), 1) %>% mutate(grp = g)
    }) %>% bind_rows() %>%
      mutate(`EMA 감소 (%)` = 100 * (1 - EMA)) %>%
      select(grp, `EMA 감소 (%)`, MCF, MCHC) %>%
      pivot_longer(-grp) %>%
      ggplot(aes(grp, value)) + geom_col(fill = "#37474f") +
      facet_wrap(~name, scales = "free_y") +
      geom_hline(data = tibble(name = "EMA 감소 (%)", y = 16),
                 aes(yintercept = y), linetype = 2, colour = "red") +
      coord_flip() + labs(x = NULL, y = NULL) + thm
  })

  ## ---- 10 crises -------------------------------------------------------
  output$crisis <- renderPlot({
    lapply(names(GENO)[c(1, 3, 4, 6)], function(g) {
      p <- par_now(); p["fdef"] <- GENO[[g]]["fdef"]
      p["f_b3ves"] <- GENO[[g]]["f_b3ves"]; p["parvo_t"] <- -1
      run_intervention(p, c(parvo_t = 10, parvo_dur = 8),
                       post = 90, delta = 0.5) %>% mutate(grp = g)
    }) %>% bind_rows() %>%
      select(time, grp, Hb, RETpct, TBIL) %>% pivot_longer(c(-time, -grp)) %>%
      ggplot(aes(time, value, colour = grp)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "감염 후 일수", y = NULL, colour = NULL) + thm
  })

  output$crisistab <- renderTable({
    lapply(names(GENO)[c(1, 3, 4, 6)], function(g) {
      p <- par_now(); p["fdef"] <- GENO[[g]]["fdef"]
      p["f_b3ves"] <- GENO[[g]]["f_b3ves"]; p["parvo_t"] <- -1
      b <- tail(run_one(p, end = 1400, delta = 1400), 1)
      d <- run_intervention(p, c(parvo_t = 10, parvo_dur = 8),
                            post = 90, delta = 0.5)
      tibble(유전자형 = g, `기저 Hb` = b$Hb, `수명 (일)` = b$LIFE,
             `예측 ΔHb = Hb·8/수명` = b$Hb * 8 / b$LIFE,
             `모델 최저 Hb` = min(d$Hb),
             `실제 ΔHb` = b$Hb - min(d$Hb))
    }) %>% bind_rows()
  }, digits = 2)

  ## ---- 11 drug ---------------------------------------------------------
  output$drug <- renderPlot({
    lapply(c(0, 5, 20, 50, 100), function(d) {
      p <- par_now(); p["dose_m"] <- d
      tail(run_one(p, end = 900, delta = 900), 1) %>% mutate(dose = d)
    }) %>% bind_rows() %>%
      select(dose, MITA, Hb, RETpct, TBIL, AREA, LIFE, DCRIT) %>%
      pivot_longer(-dose) %>%
      ggplot(aes(dose, value)) + geom_line(linewidth = 1, colour = "#01579b") +
      geom_point() + facet_wrap(~name, scales = "free_y") +
      labs(x = "미타피밧 (mg BID)", y = NULL) + thm
  })

  output$drugtab <- renderTable({
    base <- par_now()
    arms <- list(`기저` = c(),
                 `미타피밧 100 mg` = c(dose_m = 100),
                 `비장절제` = c(spl_frac = 0),
                 `병용` = c(dose_m = 100, spl_frac = 0))
    lapply(names(arms), function(a) {
      p <- base; if (length(arms[[a]])) p[names(arms[[a]])] <- arms[[a]]
      s <- tail(run_one(p, end = 1200, delta = 1200), 1)
      tibble(요법 = a, Hb = s$Hb, `망상 %` = s$RETpct, `수명` = s$LIFE,
             `ATP` = NA_real_, `총빌리루빈` = s$TBIL)
    }) %>% bind_rows()
  }, digits = 2)

  ## ---- 12 stones and iron ----------------------------------------------
  output$stone <- renderPlot({
    lapply(c("1", "0.287"), function(u) {
      p <- par_now(); p["ugt_f"] <- as.numeric(u)
      run_one(p, end = 40 * 365, delta = 90) %>%
        mutate(grp = ifelse(u == "1", "UGT1A1 정상", "길버트 증후군"))
    }) %>% bind_rows() %>%
      select(time, grp, STONEPCT, TBIL, FELIV = TBIL) %>%
      pivot_longer(c(-time, -grp)) %>%
      ggplot(aes(time / 365, value, colour = grp)) +
      geom_line(linewidth = 1) + facet_wrap(~name, scales = "free_y") +
      labs(x = "나이 (년)", y = NULL, colour = NULL) + thm
  })

  output$stonetab <- renderTable({
    s <- ss()
    tibble(`담석 지수 (%)` = s$STONEPCT, `총빌리루빈` = s$TBIL,
           `빌리루빈 생성 (mg/day)` = s$BRPROD)
  }, digits = 2)

  ## ---- 13 scenarios ----------------------------------------------------
  output$scen <- renderPlot({
    base <- par_now()
    arms <- list(
      `1 자연 경과` = c(),
      `2 전절제` = c(spl_frac = 0),
      `3 부분절제 20%` = c(spl_frac = 0.2, k_regrow = 1 / 700),
      `4 미타피밧 100 mg` = c(dose_m = 100),
      `5 절제 + 미타피밧` = c(spl_frac = 0, dose_m = 100),
      `6 정기 수혈` = c(tx_start = 0),
      `7 엽산 결핍` = c(fol_ok = 0.30),
      `8 길버트 동반` = c(ugt_f = 0.287))
    lapply(names(arms), function(a) {
      p <- base; if (length(arms[[a]])) p[names(arms[[a]])] <- arms[[a]]
      run_one(p, end = input$tend, delta = 5) %>% mutate(arm = a)
    }) %>% bind_rows() %>%
      select(time, arm, Hb, RETpct, TBIL, LIFE, AREA, SPLEEN) %>%
      pivot_longer(c(-time, -arm)) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (일)", y = NULL, colour = NULL) + thm
  })

  output$scentab <- renderTable({
    base <- par_now()
    arms <- list(`자연 경과` = c(), `전절제` = c(spl_frac = 0),
                 `미타피밧` = c(dose_m = 100),
                 `절제 + 미타피밧` = c(spl_frac = 0, dose_m = 100))
    lapply(names(arms), function(a) {
      p <- base; if (length(arms[[a]])) p[names(arms[[a]])] <- arms[[a]]
      s <- tail(run_one(p, end = 1200, delta = 1200), 1)
      tibble(시나리오 = a, Hb = s$Hb, `망상 %` = s$RETpct,
             `수명` = s$LIFE, `빌리루빈` = s$TBIL, `비장 mL` = s$SPLEEN)
    }) %>% bind_rows()
  }, digits = 2)

  ## ---- 14 virtual population -------------------------------------------
  pop <- eventReactive(input$go, {
    set.seed(20260806)
    n <- input$npop
    tibble(id = seq_len(n),
           fdef = pmin(pmax(rnorm(n, 0.28, 0.09), 0.05), 0.55),
           fb3  = ifelse(runif(n) < 0.22, 0.15, 1.00),
           ugt  = ifelse(runif(n) < 0.12, 0.287, 1.00)) %>%
      rowwise() %>%
      mutate(res = list(tail(run_one(c(par_now(),
                                       fdef = fdef, f_b3ves = fb3,
                                       ugt_f = ugt),
                                     end = 1100, delta = 1100), 1))) %>%
      tidyr::unnest(res) %>% ungroup()
  })

  output$pop <- renderPlot({
    pop() %>% select(Hb, RETpct, TBIL, LIFE, MCHC, EMA) %>%
      mutate(`EMA 감소 %` = 100 * (1 - EMA)) %>% select(-EMA) %>%
      pivot_longer(everything()) %>%
      ggplot(aes(value)) + geom_histogram(bins = 25, fill = "#1565c0") +
      facet_wrap(~name, scales = "free") + labs(x = NULL, y = "환자 수") + thm
  })

  output$poptab <- renderTable({
    d <- pop()
    tibble(지표 = c("Hb", "망상 %", "빌리루빈", "수명", "EMA 감소 %"),
           p10 = c(quantile(d$Hb, .1), quantile(d$RETpct, .1),
                   quantile(d$TBIL, .1), quantile(d$LIFE, .1),
                   quantile(100 * (1 - d$EMA), .1)),
           중앙값 = c(median(d$Hb), median(d$RETpct), median(d$TBIL),
                      median(d$LIFE), median(100 * (1 - d$EMA))),
           p90 = c(quantile(d$Hb, .9), quantile(d$RETpct, .9),
                   quantile(d$TBIL, .9), quantile(d$LIFE, .9),
                   quantile(100 * (1 - d$EMA), .9)))
  }, digits = 2)
}

shinyApp(ui, server)
