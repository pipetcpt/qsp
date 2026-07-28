## =============================================================================
## Tardive Dyskinesia QSP — Shiny dashboard
## 지연성 운동이상증 QSP 모델 인터랙티브 대시보드 (9 탭)
##
##   shiny::runApp("td_shiny_app.R")
##
## The model itself lives in td_mrgsolve_model.R; this file only loads it
## (TD_SKIP_RUN prevents the batch analyses from running on source) and wires
## it to a UI.
##
## Tabs
##   1 환자·처방 프로파일   patient, offending drug, adjuncts
##   2 PK · 수용체 점유     concentrations, striatal D2 and VMAT2 occupancy
##   3 가소성 상태          RUP / ROS / SDAM with the latch threshold drawn
##   4 기저핵 회로          IND-GPe-STN-GPi-THAL and the EXC decomposition
##   5 AIMS · 임상 엔드포인트 AIMS, psychosis, parkinsonism, mood, QTc, adherence
##   6 시나리오 비교        up to four strategies on one axis
##   7 가역성 창            withdrawal-timing scan -> point of no return
##   8 반대 레버 프론티어   AIMS vs psychosis / parkinsonism trade-off
##   9 보정 기준·문헌 정박  calibration anchors table
## =============================================================================

library(shiny)
TD_SKIP_RUN <- TRUE
source("td_mrgsolve_model.R")

PAL <- c(aims = "#c0392b", rup = "#8e44ad", sdam = "#d35400",
         ros = "#16a085", occ = "#2471a3", occv = "#117864",
         psy = "#7f8c8d", park = "#b7950b", depr = "#5b2c6f",
         a = "#2471a3", b = "#c0392b", c = "#117864", d = "#b7950b")

lineplot <- function(x, ys, cols, labs, ylab, main, ylim = NULL, hlines = NULL,
                     vlines = NULL, lwd = 2.2) {
  if (is.null(ylim)) ylim <- range(0, unlist(ys), na.rm = TRUE)
  plot(x, ys[[1]], type = "n", ylim = ylim, xlab = "일 (day)", ylab = ylab,
       main = main, las = 1, bty = "l")
  grid(col = "gray90")
  if (!is.null(hlines))
    for (h in seq_along(hlines))
      abline(h = hlines[[h]]$y, lty = 3, col = hlines[[h]]$col)
  if (!is.null(vlines)) for (v in vlines) abline(v = v, lty = 2, col = "gray55")
  for (i in seq_along(ys)) lines(x, ys[[i]], col = cols[i], lwd = lwd)
  legend("topleft", legend = labs, col = cols, lwd = lwd, bty = "n",
         cex = 0.85)
}

## -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("지연성 운동이상증 (Tardive Dyskinesia) QSP 모델 대시보드"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Host)"),
      sliderInput("age", "연령 (year)", 20, 85, 58, step = 1),
      checkboxInput("fga", "1세대 항정신병제 수준 산화 부담 (FGA)", TRUE),
      checkboxInput("dm", "당뇨·대사증후군", FALSE),
      checkboxInput("estro", "에스트로겐 충만 (신경보호)", FALSE),
      sliderInput("gen", "유전 위험 배수 (GEN_RISK)", 0.8, 1.5, 1.0,
                  step = 0.05),
      hr(),
      h4("원인 약물 (Offending drug)"),
      sliderInput("dose1", "초기 용량 (risperidone-eq mg/일)", 0, 16, 8,
                  step = 0.5),
      sliderInput("dose2", "변경 후 용량 (mg/일)", 0, 16, 8, step = 0.5),
      sliderInput("tsw", "용량 변경 시점 (day)", 0, 1825, 730, step = 5),
      sliderInput("flai", "LAI 로 투여되는 비율", 0, 1, 0, step = 0.1),
      checkboxInput("esc", "임상적 증량 정책 (PSYCH>0.20 → 증량)", FALSE),
      hr(),
      h4("TD 치료 (TD-directed therapy)"),
      selectInput("vmat", "VMAT2 억제제",
                  c("없음" = "none", "발베나진 (valbenazine)" = "val",
                    "듀테트라베나진 (deutetrabenazine)" = "dtb")),
      sliderInput("vdose", "VMAT2 억제제 용량 (mg/일)", 0, 120, 80, step = 4),
      sliderInput("tvmat", "VMAT2 억제제 시작 (day)", 0, 1825, 730, step = 5),
      selectInput("cyp", "CYP2D6 표현형",
                  c("정상 (1.0x)" = "1.0", "UM (1.6x)" = "1.6",
                    "IM (0.7x)" = "0.7", "PM (0.5x)" = "0.5",
                    "PM + CYP3A4 억제 (0.35x)" = "0.35")),
      sliderInput("clz", "클로자핀 전환 용량 (mg/일, 0=없음)", 0, 600, 0,
                  step = 25),
      hr(),
      h4("보조 (Adjuncts)"),
      checkboxInput("gkb", "은행엽 EGb761 240 mg (항산화)", FALSE),
      sliderInput("ama", "아만타딘 (mg/일)", 0, 400, 0, step = 50),
      sliderInput("ach", "벤즈트로핀 (mg/일)", 0, 4, 0, step = 0.5),
      checkboxInput("bont", "보툴리눔 독소 (국소)", FALSE),
      hr(),
      sliderInput("days", "시뮬레이션 기간 (day)", 365, 3650, 1825,
                  step = 65)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("1 환자·처방", br(), verbatimTextOutput("profile"),
                 plotOutput("p_summary", height = "340px")),
        tabPanel("2 PK · 점유", br(), plotOutput("p_pk", height = "300px"),
                 plotOutput("p_occ", height = "300px")),
        tabPanel("3 가소성 상태", br(),
                 plotOutput("p_plast", height = "340px"),
                 verbatimTextOutput("latchtxt")),
        tabPanel("4 기저핵 회로", br(), plotOutput("p_bg", height = "300px"),
                 plotOutput("p_exc", height = "300px")),
        tabPanel("5 AIMS · 엔드포인트", br(),
                 plotOutput("p_aims", height = "300px"),
                 plotOutput("p_ep", height = "300px")),
        tabPanel("6 시나리오 비교", br(),
                 plotOutput("p_scen", height = "420px"),
                 tableOutput("t_scen")),
        tabPanel("7 가역성 창", br(),
                 sliderInput("wgrid", "중단 시점 격자 (day)", 60, 1460,
                             c(60, 730), step = 10, width = "60%"),
                 plotOutput("p_window", height = "380px"),
                 tableOutput("t_window")),
        tabPanel("8 반대 레버", br(), plotOutput("p_front", height = "420px"),
                 tableOutput("t_front")),
        tabPanel("9 보정 기준", br(), tableOutput("t_anchor"),
                 verbatimTextOutput("disclaim"))
      )
    )
  )
)

## -----------------------------------------------------------------------------
server <- function(input, output, session) {

  pset <- reactive({
    p <- list(
      AGE = input$age, RISK_FGA = if (input$fga) 1.6 else 1.0,
      DM_RISK = as.numeric(input$dm), ESTROGEN = as.numeric(input$estro),
      GEN_RISK = input$gen,
      DOSE_AP = input$dose1, DOSE_AP2 = input$dose2, TSW_AP = input$tsw,
      F_LAI = input$flai, ESC_ON = as.numeric(input$esc),
      CYP2D6 = as.numeric(input$cyp),
      DOSE_AMA = input$ama, TSTART_AMA = if (input$ama > 0) 0 else 1e6,
      DOSE_ACH = input$ach, TSTART_ACH = if (input$ach > 0) 0 else 1e6,
      GKB_ON = as.numeric(input$gkb),
      TSTART_GKB = if (input$gkb) 0 else 1e6,
      BONT_ON = as.numeric(input$bont),
      TSTART_BONT = if (input$bont) input$tvmat else 1e6)
    if (input$vmat == "val") {
      p$DOSE_VAL <- input$vdose; p$TSTART_VAL <- input$tvmat
    } else if (input$vmat == "dtb") {
      p$DOSE_DTB <- input$vdose; p$TSTART_DTB <- input$tvmat
    }
    if (input$clz > 0) {
      p$DOSE_CLZ <- input$clz; p$TSTART_CLZ <- input$tsw
    }
    p
  })

  run <- reactive(do.call(sim, c(list(days = input$days), pset())))

  ## ---- tab 1 -------------------------------------------------------------
  output$profile <- renderText({
    d <- run(); n <- nrow(d)
    lat <- d$time[which(d$SDAM > 0.49)[1]]
    paste0(
      sprintf("연령 %.0f세 · %s · 유전위험 %.2f\n", input$age,
              if (input$fga) "FGA 수준 산화 부담" else "SGA 수준", input$gen),
      sprintf("D2 점유율: 시작 %.3f → 종료 %.3f | VMAT2 점유율 종료 %.3f\n",
              at_day(d, "OCCo", min(60, input$days)),
              d$OCCo[n], d$OCCVo[n]),
      sprintf("AIMS: 1년 %.2f · 2년 %.2f · 종료 %.2f (최대 %.2f)\n",
              at_day(d, "AIMS", 365), at_day(d, "AIMS", 730), d$AIMS[n],
              max(d$AIMS)),
      sprintf("RUP %.3f · ROS %.3f · SDAM %.3f (잠금 문턱 0.50-0.55)\n",
              d$RUP[n], d$ROS[n], d$SDAM[n]),
      sprintf("잠금 통과 시점: %s\n",
              ifelse(is.na(lat), "통과하지 않음 (가역 상태 유지)",
                     paste0(round(lat), "일"))),
      sprintf("PSYCH %.3f · PARK %.3f · DEPR %.3f · QTc %.2f ms · ADHER %.3f\n",
              d$PSYCH[n], d$PARK[n], d$DEPR[n], d$QTC[n], d$ADHER[n]),
      sprintf("누적 점유-일수 %.0f · 가소성 구동-일수 %.0f\n",
              d$CUMOCC[n],
              sum(diff(d$time) * (head(d$PLAST_DRIVE, -1) +
                                    tail(d$PLAST_DRIVE, -1)) / 2)))
  })

  output$p_summary <- renderPlot({
    d <- run()
    lineplot(d$time, list(d$AIMS, 26 * d$SDAM, 10 * d$RUP),
             c(PAL["aims"], PAL["sdam"], PAL["rup"]),
             c("AIMS (0-26)", "SDAM x 26", "RUP x 10"),
             "지표", "질병 궤적 요약 (Disease trajectory)",
             vlines = c(input$tsw, input$tvmat))
  })

  ## ---- tab 2 -------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- run()
    lineplot(d$time, list(d$C_APo, d$C_NBIo, d$C_HTBo, d$C_CLZo / 10),
             c(PAL["occ"], PAL["c"], PAL["b"], PAL["psy"]),
             c("항정신병제 활성분획 (ng/mL)", "NBI-98782 (ng/mL)",
               "(a+b)-HTBZ (ng/mL)", "클로자핀 /10 (ng/mL)"),
             "농도 (ng/mL)", "약물 농도 (Pharmacokinetics)",
             vlines = c(input$tsw, input$tvmat))
  })
  output$p_occ <- renderPlot({
    d <- run()
    lineplot(d$time, list(d$OCCo, d$OCCVo, d$PLAST_DRIVE),
             c(PAL["occ"], PAL["occv"], PAL["rup"]),
             c("선조체 D2 점유율", "VMAT2 점유율",
               "가소성 구동 Hill(OCC;0.70,4)"),
             "분율 (-)", "수용체 점유율과 가소성 구동", ylim = c(0, 1),
             hlines = list(list(y = 0.65, col = "gray40"),
                           list(y = 0.80, col = "gray40")),
             vlines = c(input$tsw, input$tvmat))
    text(input$days * 0.02, 0.67, "치료 문턱 65%", cex = 0.75, adj = 0)
    text(input$days * 0.02, 0.82, "EPS 문턱 80%", cex = 0.75, adj = 0)
  })

  ## ---- tab 3 -------------------------------------------------------------
  output$p_plast <- renderPlot({
    d <- run()
    lineplot(d$time, list(d$RUP, d$SDAM, d$ROS / 4),
             c(PAL["rup"], PAL["sdam"], PAL["ros"]),
             c("RUP (D2 초민감성)", "SDAM (구조 손상)", "ROS / 4"),
             "상태 (-)", "두 개의 기억: 가역 RUP vs 준-비가역 SDAM",
             hlines = list(list(y = 0.52, col = PAL["sdam"])),
             vlines = c(input$tsw, input$tvmat))
    text(input$days * 0.02, 0.55, "잠금 문턱 SDAM* ~ 0.52", cex = 0.8,
         adj = 0, col = PAL["sdam"])
  })
  output$latchtxt <- renderText({
    d <- run(); n <- nrow(d)
    paste0(
      "SDAM 방정식은 이중안정입니다: 0 (수리 우세) · 문턱 0.50-0.55 · ",
      "상위 0.85-0.90 (자기유지).\n",
      sprintf("현재 SDAM = %.3f → %s\n", d$SDAM[n],
              ifelse(d$SDAM[n] > 0.55, "잠금 상태 (원인약물 중단으로도 완전 회복 불가)",
                     ifelse(d$SDAM[n] > 0.30,
                            "문턱 접근 중 (지금이 개입 시점)",
                            "가역 영역"))),
      sprintf("RUP 은 구동이 사라지면 τ = %.0f일로 감쇠합니다 (SDAM 이 높으면 더 느립니다).\n",
              300 / (1 + 0.8 * d$SDAM[n])))
  })

  ## ---- tab 4 -------------------------------------------------------------
  output$p_bg <- renderPlot({
    d <- run()
    lineplot(d$time, list(d$IND, d$GPE, d$STN, d$GPI, d$THAL),
             c(PAL["a"], PAL["b"], PAL["c"], PAL["d"], PAL["aims"]),
             c("간접로 MSN (IND)", "GPe", "STN", "GPi/SNr", "시상 THAL"),
             "상대 활성 (-)", "기저핵-시상 루프",
             hlines = list(list(y = 1, col = "gray40")),
             vlines = c(input$tsw, input$tvmat))
  })
  output$p_exc <- renderPlot({
    d <- run()
    da_n <- d$DA_No
    stim <- da_n * (1 + 2.2 * d$RUP) * (1 - 0.55 * d$OCCo)
    lineplot(d$time, list(da_n, stim, pmax(stim - 1, 0) + 1.2 * d$SDAM,
                          1.2 * d$SDAM),
             c(PAL["occv"], PAL["rup"], PAL["aims"], PAL["sdam"]),
             c("시냅스 도파민 (기저=1)", "D2STIM", "EXC (총 과잉자극)",
               "EXC 중 구조 손상 기여분"),
             "지표 (-)", "과잉자극의 분해: 가면 효과와 구조적 하한",
             vlines = c(input$tsw, input$tvmat))
  })

  ## ---- tab 5 -------------------------------------------------------------
  output$p_aims <- renderPlot({
    d <- run()
    lineplot(d$time, list(d$AIMS), PAL["aims"], "AIMS 운동증 총점 (0-26)",
             "AIMS", "AIMS 궤적 (측정 지연 τ = 7 d 포함)",
             vlines = c(input$tsw, input$tvmat))
  })
  output$p_ep <- renderPlot({
    d <- run()
    lineplot(d$time, list(d$PSYCH, d$PARK, d$DEPR, d$ADHER, d$QTC / 10),
             c(PAL["psy"], PAL["park"], PAL["depr"], PAL["c"], PAL["b"]),
             c("정신증 부담 PSYCH", "파킨슨증 PARK", "기분 부담 DEPR",
               "순응도 ADHER", "QTc/10 (ms)"),
             "지표 (-)", "동반 엔드포인트 (경쟁하는 위해)",
             vlines = c(input$tsw, input$tvmat))
  })

  ## ---- tab 6 -------------------------------------------------------------
  scen_runs <- reactive({
    base <- pset(); base$DOSE_AP2 <- base$DOSE_AP; base$TSW_AP <- 1e6
    base$DOSE_VAL <- NULL; base$DOSE_DTB <- NULL; base$DOSE_CLZ <- NULL
    tsw <- input$tsw
    arms <- list(
      "계속 유지 (continue)" = list(),
      "50% 감량" = list(DOSE_AP2 = input$dose1 / 2, TSW_AP = tsw),
      "완전 중단" = list(DOSE_AP2 = 0, TSW_AP = tsw),
      "클로자핀 전환" = list(DOSE_AP2 = 0, TSW_AP = tsw, DOSE_CLZ = 350,
                        TSTART_CLZ = tsw),
      "발베나진 80 mg 추가" = list(DOSE_VAL = 80, TSTART_VAL = tsw),
      "클로자핀 + 발베나진" = list(DOSE_AP2 = 0, TSW_AP = tsw,
                            DOSE_CLZ = 350, TSTART_CLZ = tsw,
                            DOSE_VAL = 80, TSTART_VAL = tsw))
    lapply(arms, function(a)
      do.call(sim, c(list(days = input$days), modifyList(base, a))))
  })

  output$p_scen <- renderPlot({
    rs <- scen_runs()
    cols <- c("#7f8c8d", "#2471a3", "#c0392b", "#117864", "#b7950b",
              "#8e44ad")
    lineplot(rs[[1]]$time, lapply(rs, function(d) d$AIMS), cols, names(rs),
             "AIMS (0-26)", "전략별 AIMS 궤적 (같은 환자, 같은 시점)",
             vlines = input$tsw)
  })
  output$t_scen <- renderTable({
    rs <- scen_runs()
    data.frame(
      전략 = names(rs),
      `AIMS 변경시점` = sapply(rs, function(d) at_day(d, "AIMS", input$tsw)),
      `AIMS +6주` = sapply(rs, function(d) at_day(d, "AIMS", input$tsw + 42)),
      `AIMS 종료` = sapply(rs, function(d) tail(d$AIMS, 1)),
      `PSYCH 종료` = sapply(rs, function(d) tail(d$PSYCH, 1)),
      `PARK 종료` = sapply(rs, function(d) tail(d$PARK, 1)),
      `RUP 종료` = sapply(rs, function(d) tail(d$RUP, 1)),
      `SDAM 종료` = sapply(rs, function(d) tail(d$SDAM, 1)),
      check.names = FALSE)
  }, digits = 3)

  ## ---- tab 7 -------------------------------------------------------------
  window_scan <- reactive({
    g <- round(seq(input$wgrid[1], input$wgrid[2], length.out = 12))
    p <- pset(); p$DOSE_AP2 <- 0
    res <- lapply(g, function(s) {
      p$TSW_AP <- s
      d <- do.call(sim, c(list(days = s + 2190), p))
      post <- d[d$time >= s, ]
      data.frame(stop = s, sdam = at_day(d, "SDAM", s),
                 aims_stop = at_day(d, "AIMS", s),
                 peak = max(post$AIMS),
                 peak_day = post$time[which.max(post$AIMS)] - s,
                 aims6 = at_day(d, "AIMS", s + 2190))
    })
    do.call(rbind, res)
  })

  output$p_window <- renderPlot({
    w <- window_scan()
    plot(w$stop, w$aims6, type = "b", pch = 19, col = PAL["aims"], lwd = 2,
         xlab = "원인 약물 중단 시점 (day)",
         ylab = "중단 6년 후 AIMS", las = 1, bty = "l",
         main = "가역성 창: 언제까지 중단이 효과가 있는가",
         ylim = c(0, max(w$aims6, w$peak)))
    grid(col = "gray90")
    lines(w$stop, w$peak, type = "b", pch = 1, col = PAL["park"], lwd = 2)
    lines(w$stop, w$aims_stop, type = "b", pch = 2, col = PAL["occ"],
          lwd = 1.5)
    abline(h = 2, lty = 3)
    legend("topleft", c("중단 6년 후 AIMS (잔존)", "중단 후 최고 AIMS (금단 악화)",
                        "중단 시점 AIMS"),
           col = c(PAL["aims"], PAL["park"], PAL["occ"]), pch = c(19, 1, 2),
           lwd = 2, bty = "n", cex = 0.85)
  })
  output$t_window <- renderTable({
    w <- window_scan()
    w$outcome <- ifelse(w$aims6 > 2, "지속성 TD", "회복")
    names(w) <- c("중단일", "중단시 SDAM", "중단시 AIMS", "최고 AIMS",
                  "최고 도달일", "+6년 AIMS", "결과")
    w
  }, digits = 3)

  ## ---- tab 8 -------------------------------------------------------------
  frontier <- reactive({
    base <- pset(); base$DOSE_AP2 <- base$DOSE_AP; base$TSW_AP <- 1e6
    base$DOSE_VAL <- NULL; base$DOSE_DTB <- NULL; base$DOSE_CLZ <- NULL
    tsw <- input$tsw; readout <- min(tsw + 365, input$days)
    arms <- c(
      lapply(c(0, 0.25, 0.5, 0.75, 1), function(f)
        list(nm = sprintf("D2 차단 %.0f%% 유지", 100 * (1 - f)),
             lever = "원인약물 감량",
             p = list(DOSE_AP2 = input$dose1 * (1 - f), TSW_AP = tsw))),
      lapply(c(20, 40, 80, 120), function(x)
        list(nm = sprintf("발베나진 %d mg", x), lever = "도파민 공급 차단",
             p = list(DOSE_VAL = x, TSTART_VAL = tsw))),
      list(list(nm = "클로자핀 350 전환", lever = "낮은 점유 전환",
                p = list(DOSE_AP2 = 0, TSW_AP = tsw, DOSE_CLZ = 350,
                         TSTART_CLZ = tsw))))
    do.call(rbind, lapply(arms, function(a) {
      d <- do.call(sim, c(list(days = readout + 5), modifyList(base, a$p)))
      data.frame(strategy = a$nm, lever = a$lever,
                 AIMS = at_day(d, "AIMS", readout),
                 PSYCH = at_day(d, "PSYCH", readout),
                 PARK = at_day(d, "PARK", readout),
                 DEPR = at_day(d, "DEPR", readout),
                 ADHER = at_day(d, "ADHER", readout),
                 RUP = at_day(d, "RUP", readout))
    }))
  })

  output$p_front <- renderPlot({
    f <- frontier()
    cols <- c("원인약물 감량" = PAL[["occ"]],
              "도파민 공급 차단" = PAL[["occv"]],
              "낮은 점유 전환" = PAL[["aims"]])
    plot(f$AIMS, f$PSYCH, pch = 19, cex = 1.6, col = cols[f$lever],
         xlab = "AIMS (낮을수록 좋음)", ylab = "정신증 부담 PSYCH",
         main = "반대편 레버: 같은 AIMS 감소가 서로 다른 대가를 지불한다",
         las = 1, bty = "l", xlim = range(f$AIMS) + c(-1, 1),
         ylim = c(0, max(f$PSYCH) + 0.1))
    grid(col = "gray90")
    text(f$AIMS, f$PSYCH, f$strategy, pos = 3, cex = 0.78)
    points(f$AIMS, f$PARK, pch = 1, cex = 1.4, col = cols[f$lever])
    legend("topright", c(names(cols), "속 빈 원 = 파킨슨증 PARK"),
           col = c(cols, "black"), pch = c(19, 19, 19, 1), bty = "n",
           cex = 0.85)
  })
  output$t_front <- renderTable(frontier(), digits = 3)

  ## ---- tab 9 -------------------------------------------------------------
  output$t_anchor <- renderTable({
    data.frame(
      `보정 대상` = c(
        "리스페리돈 활성분획 PK", "선조체 D2 점유 EC50",
        "치료/EPS 점유 문턱", "클로자핀 D2 점유",
        "발베나진 80 mg 노출", "듀테트라베나진 36 mg 노출",
        "KINECT-3 (발베나진 80 mg, 6주)",
        "ARM-TD/AIM-TD (듀테트라베나진 36 mg)",
        "연간 TD 발생률", "금단 유발성 운동증 시점",
        "중단 후 지속성", "은행엽 EGb761 240 mg"),
      `문헌 값` = c(
        "CL/F ~ 5 L/h, t1/2 ~ 20 h, 4 mg → ~33 ng/mL",
        "~12 ng/mL (4 mg ~ 73%, 8 mg ~ 86%)",
        "65% / 78-80%", "350 mg → ~390 ng/mL, D2 ~ 20-40%",
        "NBI-98782 Cave ~ 24 ng/mL (PM 약 2배)",
        "(a+b)-HTBZ ~ 15 ng/mL, QTc ~ 4 ms",
        "AIMS -3.2 (기저 ~10)", "AIMS ~ -3.0",
        "FGA 5-6%/년, SGA 3-4%/년, 55세 이상 25-30%/년",
        "감량 후 4-6주 내 최고",
        "약 1/3 회복, 나머지 수년 지속", "AIMS 개선 (RCT)"),
      `모델 파라미터 / 출력` = c(
        "CL_AP 120 L/일, V2 100 L", "EC50_D2 = 12, HILL_D2 = 1.25",
        "OCC50_P = 0.58, OCC50_PK = 0.78",
        "CL_CLZ 900 L/일, EC50_D2_CLZ = 900 → OCC 0.28",
        "FM_VAL 0.30, CL_NBI 1000 L/일, CYP2D6 x0.5",
        "FM_DTB 0.50, CL_HTB 1200 L/일, QT_HTB 0.25",
        "모델: -3.51 (-33%), 기저 10.59",
        "모델: -3.09 (-29%) at 6주",
        "RISK_FGA 1.6, RISKMOD = 1+0.015(age-40)",
        "모델 최고 시점 = 중단 후 30일",
        "SDAM 이중안정 (문턱 0.52), KOUT_S = 3e-4",
        "ANTIOX_MAX 0.85 → 잠금 400일 → 780일"),
      check.names = FALSE)
  })
  output$disclaim <- renderText(paste(
    "이 대시보드의 모든 수치는 교육·연구 목적의 QSP 모델 출력입니다.",
    "\n임상 의사결정·처방·규제 제출에 사용할 수 없습니다.",
    "\n문헌 목록은 td_references.md, 독립 수치 검증은",
    "td_reference_check.py 를 참조하십시오."))
}

shinyApp(ui, server)
