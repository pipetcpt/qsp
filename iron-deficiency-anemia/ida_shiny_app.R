## =====================================================================
##  ida_shiny_app.R
##  Iron Deficiency Anaemia QSP model — interactive dashboard
##  철결핍성 빈혈 QSP 모델 — 인터랙티브 대시보드
## ---------------------------------------------------------------------
##  Eight tabs:
##    1. Patient profile      — set up the patient and read the baseline
##    2. Oral iron PK         — one dose: the gate opening and closing
##    3. The dosing interval  — fractional vs total absorption, the
##                              central trade-off of the model
##    4. Haemoglobin response — 12-week trajectories and milestones
##    5. Iron distribution    — where the iron actually goes
##    6. Biomarkers           — ferritin, TSAT, sTfR, CHr, hepcidin
##    7. Safety: phosphate    — the carboxymaltose FGF23 penalty
##    8. Scenario comparison  — any two regimens side by side
##
##  Run with:  shiny::runApp("ida_shiny_app.R")
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ---------------------------------------------------------------------
##  the model: sourced from the companion file so there is ONE definition
## ---------------------------------------------------------------------
model_file <- "ida_mrgsolve_model.R"
if (!file.exists(model_file))
  stop("ida_mrgsolve_model.R must sit next to this app.")
src <- readLines(model_file)
eval(parse(text = paste(src[1:max(grep("^ida <- mcode", src))], collapse = "\n")))

WK <- 7 * 24

theme_ida <- function() {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          legend.position = "bottom", legend.title = element_blank(),
          strip.text = element_text(face = "bold"))
}
PAL <- c("#c62828", "#1565c0", "#2e7d32", "#6a1b9a", "#ef6c00",
         "#00838f", "#ad1457", "#4e342e")

oral_ev <- function(dose_mg, interval_h, n_dose, start_h = 0)
  ev(time = start_h, amt = dose_mg, cmt = "A_LUM",
     ii = interval_h, addl = n_dose - 1)
iv_ev <- function(dose_mg, times_h = 0)
  Reduce(`+`, lapply(times_h, function(t)
    ev(time = t, amt = dose_mg, cmt = "A_COL") +
    ev(time = t, amt = dose_mg, cmt = "CUM_IV")))

## build the event object for a named regimen over `weeks` weeks
build_ev <- function(route, dose, interval_h, weeks, iv_n) {
  if (route == "none") return(NULL)
  if (route == "oral") {
    n <- max(1, floor(weeks * 24 * 7 / interval_h))
    oral_ev(dose, interval_h, n)
  } else {
    iv_ev(dose, WK * (seq_len(iv_n) - 1))
  }
}

simulate_case <- function(pars, events, end_h, delta = 12, hmax = 0.5) {
  m <- if (length(pars)) do.call(param, c(list(ida), pars)) else ida
  out <- if (is.null(events)) mrgsim(m, end = end_h, delta = delta, hmax = hmax)
         else mrgsim(m, events = events, end = end_h, delta = delta, hmax = hmax)
  as_tibble(out)
}

## the patient's own baseline: re-equilibrate the model at their parameters
equilibrate <- function(pars, days = 730) {
  m <- if (length(pars)) do.call(param, c(list(ida), pars)) else ida
  tail(as_tibble(mrgsim(m, end = days * 24, delta = 24, hmax = 1)), 1)
}

at_h <- function(d, h) d[which.min(abs(d$time - h)), ]

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel("철결핍성 빈혈 QSP 모델 · Iron Deficiency Anaemia QSP Dashboard"),
  tags$p(style = "color:#555;margin-top:-8px;",
         "37-compartment mrgsolve model · absorption = luminal iron x DMT1 capacity x ",
         "enterocyte export capacity, the third factor set by the hepcidin the ",
         "previous dose generated."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      sliderInput("vbleed", "지속 출혈 Blood loss (mL/day)",
                  min = 0.2, max = 12, value = 6.375, step = 0.125),
      helpText("6.375 mL/day ~ 178 mL/cycle. 0.6 = replete reference."),
      sliderInput("diet", "식이 철 Dietary iron (mg/day)",
                  min = 4, max = 30, value = 16, step = 1),
      sliderInput("il6", "IL-6 input (pg/mL/h) — 염증",
                  min = 0, max = 4, value = 0, step = 0.25),
      helpText("IL-6 steady state = 20 x this value."),
      sliderInput("tmprss6", "TMPRSS6 활성 (1 = normal, <1 = IRIDA)",
                  min = 0.15, max = 1, value = 1, step = 0.05),
      hr(),
      h4("치료 (Therapy)"),
      radioButtons("route", "경로 Route",
                   c("없음 none" = "none", "경구 oral" = "oral",
                     "정맥 IV" = "iv"), selected = "oral"),
      conditionalPanel(
        "input.route == 'oral'",
        sliderInput("odose", "1회 용량 Elemental iron per dose (mg)",
                    min = 15, max = 200, value = 65, step = 5),
        radioButtons("ointerval", "투여 간격 Interval",
                     c("q12h (BID)" = "12", "q24h (daily)" = "24",
                       "q48h (alternate day)" = "48", "q72h" = "72"),
                     selected = "24"),
        checkboxInput("adhere", "GI 부작용에 의한 순응도 저하 반영", TRUE)
      ),
      conditionalPanel(
        "input.route == 'iv'",
        sliderInput("ivdose", "1회 용량 IV dose (mg)",
                    min = 100, max = 1500, value = 1000, step = 50),
        sliderInput("ivn", "투여 횟수 (weekly)", min = 1, max = 5, value = 1),
        radioButtons("ivprod", "제형 Product",
                     c("Ferric carboxymaltose (FCM)" = "1",
                       "Ferric derisomaltose" = "0",
                       "Iron sucrose" = "0.10"), selected = "1"),
        helpText("Only carboxymaltose inhibits FGF23 cleavage.")
      ),
      hr(),
      sliderInput("weeks", "관찰 기간 Follow-up (weeks)",
                  min = 4, max = 24, value = 12, step = 1),
      checkboxInput("reeq", "환자 파라미터로 기저상태 재평형 (느림, ~20 s)", FALSE),
      helpText("Off: start from the calibrated IDA reference state.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",
        tabPanel("1 · 환자 프로파일",
                 br(), h4("기저 상태 Baseline state"),
                 tableOutput("baseline_tbl"),
                 helpText("Reference ranges (adult women): Hb 12.0-15.5 g/dL · MCH 27-33 pg ·",
                          "TSAT 20-45 % · ferritin 15-150 ng/mL · TIBC 250-400 ug/dL ·",
                          "hepcidin 1-20 ng/mL · sTfR 0.8-1.8 mg/L · CHr > 28 pg."),
                 hr(), h4("철 수지 Iron balance at baseline"),
                 plotOutput("balance_plot", height = "300px"),
                 helpText("A negative balance is what drives the disease: the model's IDA",
                          "state is the equilibrium where the falling haemoglobin mass has",
                          "reduced the bleeding iron loss enough to match absorption.")),
        tabPanel("2 · 경구 철 PK",
                 br(), h4("단일 경구 용량 — the gate opening and closing"),
                 plotOutput("single_plot", height = "620px"),
                 tableOutput("single_tbl"),
                 helpText("Serum iron peaks first, hepcidin follows it by ~5 h, and the loss",
                          "of enterocyte export capacity outlasts both because it recovers",
                          "on a 1-2 day clock (KSYN_FPE) rather than hepcidin's 2.5 h.")),
        tabPanel("3 · 투여 간격 (핵심)",
                 br(), h4("The central trade-off: efficiency versus delivery"),
                 plotOutput("interval_plot", height = "400px"),
                 tableOutput("interval_tbl"),
                 helpText("Fractional absorption and total absorbed iron are ordered in",
                          "OPPOSITE directions by dosing frequency. Alternate-day dosing is",
                          "an efficiency and tolerability optimum, not a delivery optimum."),
                 hr(), h4("Refractory window — probe dose after a conditioning dose"),
                 plotOutput("probe_plot", height = "300px")),
        tabPanel("4 · 혈액학 반응",
                 br(), h4("헤모글로빈 및 적혈구 지표"),
                 plotOutput("hb_plot", height = "620px"),
                 tableOutput("hb_tbl")),
        tabPanel("5 · 철 분포",
                 br(), h4("Where the iron goes"),
                 plotOutput("dist_plot", height = "440px"),
                 helpText("Marrow supply is dominated by macrophage recycling of senescent",
                          "red cells; the gut contributes a few mg/day at best. This is why",
                          "intravenous dose escalation buys stores rather than speed."),
                 hr(), h4("조직 철 vs 헤모글로빈 회복 (fatigue / restless legs)"),
                 plotOutput("tissue_plot", height = "320px")),
        tabPanel("6 · 바이오마커",
                 br(), h4("진단 지표 Diagnostic read-outs"),
                 plotOutput("biomarker_plot", height = "620px"),
                 helpText("Under inflammation ferritin rises while TSAT stays low: the",
                          "ferritin/TSAT dissociation of functional iron deficiency.",
                          "sTfR and CHr remain interpretable when ferritin does not.")),
        tabPanel("7 · 안전성: 인 (phosphate)",
                 br(), h4("FGF23 — 인 — 비타민 D 축"),
                 plotOutput("phos_plot", height = "560px"),
                 tableOutput("phos_tbl"),
                 helpText("Carboxymaltose inhibits FGF23 cleavage; derisomaltose does not.",
                          "Both give the same haemoglobin response, so the efficacy endpoint",
                          "cannot distinguish them — only the mechanism can.")),
        tabPanel("8 · 시나리오 비교",
                 br(), h4("두 처방 직접 비교 Compare two regimens"),
                 fluidRow(
                   column(6, selectInput("cmpA", "A",
                     c("oral 65 mg daily", "oral 65 mg BID", "oral 130 mg daily",
                       "oral 130 mg alternate", "oral 65 mg alternate",
                       "IV FCM 1000 mg", "IV FCM 750 mg x2",
                       "IV derisomaltose 1000 mg", "IV sucrose 200 mg x5",
                       "no treatment"), selected = "oral 65 mg daily")),
                   column(6, selectInput("cmpB", "B",
                     c("oral 65 mg daily", "oral 65 mg BID", "oral 130 mg daily",
                       "oral 130 mg alternate", "oral 65 mg alternate",
                       "IV FCM 1000 mg", "IV FCM 750 mg x2",
                       "IV derisomaltose 1000 mg", "IV sucrose 200 mg x5",
                       "no treatment"), selected = "IV FCM 1000 mg"))),
                 plotOutput("cmp_plot", height = "560px"),
                 tableOutput("cmp_tbl"))
      )
    )
  )
)

## =====================================================================
##  server
## =====================================================================
server <- function(input, output, session) {

  patient_pars <- reactive({
    p <- list(VBLEED = input$vbleed, DIET = input$diet,
              IL6_IN = input$il6, TMPRSS6 = input$tmprss6)
    if (identical(input$route, "iv")) p$FCM_FGF <- as.numeric(input$ivprod)
    if (identical(input$route, "oral") && !isTRUE(input$adhere)) p$EMAX_ADH <- 0
    p
  })

  ## the patient's baseline: either the shipped IDA state or re-equilibrated
  base_state <- reactive({
    if (isTRUE(input$reeq)) {
      withProgress(message = "Equilibrating the patient's baseline...", {
        equilibrate(patient_pars())
      })
    } else simulate_case(patient_pars(), NULL, 1, delta = 1)[1, ]
  })

  main_sim <- reactive({
    ev_obj <- build_ev(input$route,
                       if (input$route == "oral") input$odose else input$ivdose,
                       as.numeric(input$ointerval), input$weeks, input$ivn)
    simulate_case(patient_pars(), ev_obj, input$weeks * WK, delta = 6)
  })

  ## ---- 1. patient profile -------------------------------------------
  output$baseline_tbl <- renderTable({
    b <- base_state()
    tibble(
      metric = c("Haemoglobin (g/dL)", "MCH (pg)", "RBC (1e12/L)",
               "Serum iron (ug/dL)", "TIBC (ug/dL)", "TSAT (%)",
               "Hepcidin (ng/mL)", "Ferritin (ng/mL)", "Storage iron (mg)",
               "Tissue iron (mg)", "Reticulocytes (1e9/L)", "CHr (pg)",
               "sTfR (mg/L)", "EPO (mIU/mL)", "Absorption (mg/day)",
               "Marrow iron use (mg/day)", "Enterocyte export capacity",
               "Total body iron (mg)"),
      value = c(b$HB, b$MCH, b$RBCC, b$SI, b$TIBC, b$TSAT, b$HEP, b$FERR,
             b$STORES, b$A_TISS, b$RET_ABS, b$CHR, b$STFR, b$EPO,
             b$ABS_DAY, b$FE_MARROW * 24, b$FPN_ENT, b$BODY_FE)
    )
  }, digits = 3)

  output$balance_plot <- renderPlot({
    b <- base_state()
    absorbed <- b$ABS_DAY
    basal <- 0.0417 * 24
    bleed <- 3.47 * (input$vbleed / (1000 * 24 * 3.9)) * (b$HRBC + b$HRET) * 24
    d <- tibble(term = factor(c("absorbed", "epithelial loss", "blood loss", "net"),
                              levels = c("absorbed", "epithelial loss",
                                         "blood loss", "net")),
                mg = c(absorbed, -basal, -bleed, absorbed - basal - bleed))
    ggplot(d, aes(term, mg, fill = mg > 0)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = 0) +
      geom_text(aes(label = sprintf("%+.2f", mg)), vjust = ifelse(d$mg > 0, -0.4, 1.3)) +
      scale_fill_manual(values = c("TRUE" = "#2e7d32", "FALSE" = "#c62828"),
                        guide = "none") +
      labs(title = "Daily iron balance at the current baseline",
           y = "mg iron / day", x = NULL) + theme_ida()
  })

  ## ---- 2. single oral dose ------------------------------------------
  single <- reactive({
    p <- patient_pars(); p$EMAX_ADH <- 0
    dose <- if (input$route == "oral") input$odose else 60
    with_dose <- simulate_case(p, oral_ev(dose, 24, 1), 72, delta = 0.25, hmax = 0.02)
    no_dose   <- simulate_case(p, NULL, 72, delta = 0.25, hmax = 0.02)
    list(d = with_dose, b = no_dose, dose = dose)
  })

  output$single_plot <- renderPlot({
    s <- single()
    d <- s$d %>%
      transmute(time,
                `Serum iron (ug/dL)` = SI,
                `TSAT (%)` = TSAT,
                `Hepcidin (ng/mL)` = HEP,
                `Enterocyte export capacity` = FPN_ENT,
                `Enterocyte ferritin iron (mg)` = A_EFT,
                `Cumulative absorbed (mg)` = CUM_ABS - CUM_ABS[1]) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) +
      geom_line(colour = "#1565c0", linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      geom_vline(xintercept = c(24, 48), linetype = 3, colour = "#999999") +
      labs(title = sprintf("Single %g mg elemental oral dose", s$dose),
           x = "hours", y = NULL) + theme_ida()
  })

  output$single_tbl <- renderTable({
    s <- single(); d <- s$d
    fia <- 100 * (tail(d$CUM_ABS, 1) - tail(s$b$CUM_ABS, 1)) / s$dose
    tibble(
      metric = c("Peak serum iron (ug/dL)", "Time to peak (h)",
               "Peak hepcidin (ng/mL)", "Hepcidin fold-rise",
               "Hepcidin at 24 h (fold)", "Export capacity nadir (% of baseline)",
               "Export capacity at 24 h (%)", "Export capacity at 48 h (%)",
               "Fractional iron absorption (%)", "Absorbed iron (mg)"),
      value = c(max(d$SI), d$time[which.max(d$SI)], max(d$HEP),
             max(d$HEP) / d$HEP[1], at_h(d, 24)$HEP / d$HEP[1],
             100 * min(d$FPN_ENT) / d$FPN_ENT[1],
             100 * at_h(d, 24)$FPN_ENT / d$FPN_ENT[1],
             100 * at_h(d, 48)$FPN_ENT / d$FPN_ENT[1],
             fia, fia * s$dose / 100))
  }, digits = 3)

  ## ---- 3. dosing interval ------------------------------------------
  interval_data <- reactive({
    p <- patient_pars(); p$EMAX_ADH <- 0
    dose <- if (input$route == "oral") input$odose else 65
    grid <- tibble(interval = c(12, 24, 48, 72)) %>%
      mutate(n = floor(14 * 24 / interval))
    withProgress(message = "Simulating 14-day regimens...", {
      grid$absorbed <- vapply(seq_len(nrow(grid)), function(i) {
        a <- simulate_case(p, oral_ev(dose, grid$interval[i], grid$n[i]),
                           15 * 24, delta = 24, hmax = 0.05)
        b <- simulate_case(p, NULL, 15 * 24, delta = 24, hmax = 0.05)
        tail(a$CUM_ABS, 1) - tail(b$CUM_ABS, 1)
      }, numeric(1))
      ## the same regimens WITH the tolerability penalty
      p2 <- patient_pars(); p2$EMAX_ADH <- 0.45
      grid$absorbed_adh <- vapply(seq_len(nrow(grid)), function(i) {
        a <- simulate_case(p2, oral_ev(dose, grid$interval[i], grid$n[i]),
                           15 * 24, delta = 24, hmax = 0.05)
        b <- simulate_case(p2, NULL, 15 * 24, delta = 24, hmax = 0.05)
        tail(a$CUM_ABS, 1) - tail(b$CUM_ABS, 1)
      }, numeric(1))
    })
    grid %>% mutate(total = dose * n, FIA = 100 * absorbed / total,
                    per_day = absorbed / 14, per_dose = absorbed / n,
                    FIA_adh = 100 * absorbed_adh / total,
                    dose = dose)
  })

  output$interval_plot <- renderPlot({
    g <- interval_data()
    d <- bind_rows(
      tibble(interval = g$interval, metric = "fractional absorption (%)", value = g$FIA),
      tibble(interval = g$interval, metric = "absorbed mg / day", value = g$per_day),
      tibble(interval = g$interval, metric = "absorbed mg / dose", value = g$per_dose))
    ggplot(d, aes(factor(interval), value, group = metric, colour = metric)) +
      geom_line(linewidth = 1) + geom_point(size = 3) +
      facet_wrap(~metric, scales = "free_y") +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(title = sprintf("14 days of %g mg elemental iron per dose", g$dose[1]),
           subtitle = "efficiency rises with the interval; delivery falls",
           x = "dosing interval (h)", y = NULL) + theme_ida()
  })

  output$interval_tbl <- renderTable({
    interval_data() %>%
      transmute(`interval (h)` = interval, doses = n, `total swallowed (mg)` = total,
                `absorbed (mg)` = absorbed, `FIA (%)` = FIA,
                `mg / day` = per_day, `mg / dose` = per_dose,
                `FIA with non-adherence (%)` = FIA_adh)
  }, digits = 3)

  output$probe_plot <- renderPlot({
    p <- patient_pars(); p$EMAX_ADH <- 0
    gaps <- c(4, 8, 12, 24, 36, 48, 72)
    withProgress(message = "Probe-dose experiment...", {
      ref <- {
        a <- simulate_case(p, oral_ev(60, 24, 1), 72, delta = 24, hmax = 0.05)
        b <- simulate_case(p, NULL, 72, delta = 24, hmax = 0.05)
        tail(a$CUM_ABS, 1) - tail(b$CUM_ABS, 1)
      }
      vals <- vapply(gaps, function(g) {
        both <- simulate_case(p, ev(time = 0, amt = 60, cmt = "A_LUM") +
                                 ev(time = g, amt = 60, cmt = "A_LUM"),
                              g + 72, delta = 24, hmax = 0.05)
        cond <- simulate_case(p, ev(time = 0, amt = 60, cmt = "A_LUM"),
                              g + 72, delta = 24, hmax = 0.05)
        tail(both$CUM_ABS, 1) - tail(cond$CUM_ABS, 1)
      }, numeric(1))
    })
    ggplot(tibble(gap = gaps, pct = 100 * vals / ref), aes(gap, pct)) +
      geom_hline(yintercept = 100, linetype = 2, colour = "#999999") +
      geom_line(colour = "#6a1b9a", linewidth = 1) +
      geom_point(size = 3, colour = "#6a1b9a") +
      scale_x_continuous(breaks = gaps) + ylim(60, 105) +
      labs(title = "Absorption of an identical probe dose after a conditioning dose",
           subtitle = "100 % = the same dose on a rested gut",
           x = "hours since the conditioning dose",
           y = "% of rested-gut absorption") + theme_ida()
  })

  ## ---- 4. haematology ----------------------------------------------
  output$hb_plot <- renderPlot({
    d <- main_sim() %>%
      transmute(week = time / 168,
                `Haemoglobin (g/dL)` = HB, `MCH (pg)` = MCH,
                `Reticulocytes (1e9/L)` = RET_ABS,
                `Reticulocyte Hb content CHr (pg)` = CHR,
                `EPO (mIU/mL)` = EPO, `TSAT (%)` = TSAT) %>%
      pivot_longer(-week)
    ggplot(d, aes(week, value)) +
      geom_line(colour = "#c62828", linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(title = "Haematological response", x = "weeks", y = NULL) + theme_ida()
  })

  output$hb_tbl <- renderTable({
    d <- main_sim(); hb0 <- d$HB[1]
    tt <- function(thr) {
      i <- which(d$HB - hb0 >= thr)
      if (length(i)) d$time[i[1]] / 24 else NA_real_
    }
    tibble(metric = c("Baseline Hb (g/dL)", "Hb at end (g/dL)", "dHb (g/dL)",
                    "Days to +1 g/dL", "Days to +2 g/dL",
                    "Reticulocyte peak (1e9/L)", "Day of reticulocyte peak",
                    "Final MCH (pg)", "Final TSAT (%)", "Final ferritin (ng/mL)",
                    "Storage iron at end (mg)", "Iron absorbed (mg)"),
           value = c(hb0, tail(d$HB, 1), tail(d$HB, 1) - hb0, tt(1), tt(2),
                  max(d$RET_ABS), d$time[which.max(d$RET_ABS)] / 24,
                  tail(d$MCH, 1), tail(d$TSAT, 1), tail(d$FERR, 1),
                  tail(d$STORES, 1), tail(d$CUM_ABS, 1) - d$CUM_ABS[1]))
  }, digits = 3)

  ## ---- 5. iron distribution ----------------------------------------
  output$dist_plot <- renderPlot({
    d <- main_sim() %>%
      transmute(week = time / 168,
                `red cell Hb iron` = 3.47 * (HRBC + HRET),
                `macrophage stores` = A_RES, `hepatocyte stores` = A_LIV,
                `non-erythroid tissue` = A_TISS,
                `plasma + colloid` = A_TF + A_NTBI + A_COL) %>%
      pivot_longer(-week)
    ggplot(d, aes(week, value, fill = name)) +
      geom_area(alpha = 0.85) +
      scale_fill_manual(values = PAL) +
      labs(title = "Body iron distribution over time",
           x = "weeks", y = "mg iron") + theme_ida()
  })

  output$tissue_plot <- renderPlot({
    d <- main_sim()
    z <- d[1, ]
    dd <- tibble(week = d$time / 168,
                 `haemoglobin repaired (%)` =
                   100 * (d$HB - z$HB) / max(13.86 - z$HB, 0.01),
                 `tissue iron repaired (%)` =
                   100 * (d$A_TISS - z$A_TISS) / max(386 - z$A_TISS, 0.01)) %>%
      pivot_longer(-week)
    ggplot(dd, aes(week, value, colour = name)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 100, linetype = 2, colour = "#999999") +
      scale_colour_manual(values = c("#c62828", "#4e342e")) +
      labs(title = "Non-erythroid tissue iron lags haemoglobin",
           subtitle = "the oral route repairs the blood before the tissues",
           x = "weeks", y = "% of the deficit repaired") + theme_ida()
  })

  ## ---- 6. biomarkers -----------------------------------------------
  output$biomarker_plot <- renderPlot({
    d <- main_sim() %>%
      transmute(week = time / 168,
                `Ferritin (ng/mL)` = FERR, `TSAT (%)` = TSAT,
                `sTfR (mg/L)` = STFR, `Hepcidin (ng/mL)` = HEP,
                `TIBC (ug/dL)` = TIBC, `Erythroferrone (ng/mL)` = ERFE) %>%
      pivot_longer(-week)
    ggplot(d, aes(week, value)) +
      geom_line(colour = "#00838f", linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(title = "Diagnostic biomarkers", x = "weeks", y = NULL) + theme_ida()
  })

  ## ---- 7. phosphate safety -----------------------------------------
  phos_sims <- reactive({
    dose <- if (input$route == "iv") input$ivdose else 1000
    n <- if (input$route == "iv") input$ivn else 1
    prods <- list(`carboxymaltose (FCM)` = 1, `derisomaltose` = 0,
                  `iron sucrose` = 0.10)
    withProgress(message = "Simulating the FGF23 axis...", {
      bind_rows(lapply(names(prods), function(nm) {
        p <- patient_pars(); p$FCM_FGF <- prods[[nm]]
        simulate_case(p, iv_ev(dose, WK * (seq_len(n) - 1)), 8 * WK, delta = 6) %>%
          mutate(product = nm)
      }))
    })
  })

  output$phos_plot <- renderPlot({
    d <- phos_sims() %>%
      transmute(day = time / 24, product,
                `Intact FGF23 (pg/mL)` = FGF23,
                `Serum phosphate (mg/dL)` = PHOS,
                `1,25(OH)2D (pg/mL)` = CTRIOL, `PTH (pg/mL)` = PTH) %>%
      pivot_longer(-c(day, product))
    ggplot(d, aes(day, value, colour = product)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = PAL) +
      labs(title = "FGF23 - phosphate - vitamin D axis after intravenous iron",
           x = "days", y = NULL) + theme_ida()
  })

  output$phos_tbl <- renderTable({
    phos_sims() %>% group_by(product) %>%
      summarise(`iFGF23 peak (pg/mL)` = max(FGF23),
                `day of peak` = time[which.max(FGF23)] / 24,
                `phosphate nadir (mg/dL)` = min(PHOS),
                `day of nadir` = time[which.min(PHOS)] / 24,
                `days < 2.0` = sum(PHOS < 2.0) * 6 / 24,
                `days < 2.5` = sum(PHOS < 2.5) * 6 / 24,
                `1,25D nadir (pg/mL)` = min(CTRIOL),
                `PTH peak (pg/mL)` = max(PTH), .groups = "drop")
  }, digits = 3)

  ## ---- 8. scenario comparison --------------------------------------
  REG <- list(
    "no treatment"             = list(ev = NULL, p = list()),
    "oral 65 mg daily"         = list(ev = oral_ev(65, 24, 84), p = list()),
    "oral 65 mg BID"           = list(ev = oral_ev(65, 12, 168), p = list()),
    "oral 130 mg daily"        = list(ev = oral_ev(130, 24, 84), p = list()),
    "oral 130 mg alternate"    = list(ev = oral_ev(130, 48, 42), p = list()),
    "oral 65 mg alternate"     = list(ev = oral_ev(65, 48, 42), p = list()),
    "IV FCM 1000 mg"           = list(ev = iv_ev(1000, 0), p = list(FCM_FGF = 1)),
    "IV FCM 750 mg x2"         = list(ev = iv_ev(750, c(0, WK)), p = list(FCM_FGF = 1)),
    "IV derisomaltose 1000 mg" = list(ev = iv_ev(1000, 0), p = list(FCM_FGF = 0)),
    "IV sucrose 200 mg x5"     = list(ev = iv_ev(200, WK * (0:4)), p = list(FCM_FGF = 0.10))
  )

  cmp_data <- reactive({
    pick <- c(input$cmpA, input$cmpB)
    bind_rows(lapply(pick, function(nm) {
      p <- patient_pars(); p[names(REG[[nm]]$p)] <- REG[[nm]]$p
      simulate_case(p, REG[[nm]]$ev, input$weeks * WK, delta = 12) %>%
        mutate(regimen = nm)
    }))
  })

  output$cmp_plot <- renderPlot({
    d <- cmp_data() %>%
      transmute(week = time / 168, regimen,
                `Haemoglobin (g/dL)` = HB, `TSAT (%)` = TSAT,
                `Ferritin (ng/mL)` = FERR, `Storage iron (mg)` = STORES,
                `Serum phosphate (mg/dL)` = PHOS,
                `Cumulative absorbed (mg)` = CUM_ABS - CUM_ABS[1]) %>%
      pivot_longer(-c(week, regimen))
    ggplot(d, aes(week, value, colour = regimen)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = PAL[c(2, 1)]) +
      labs(title = "Head-to-head comparison", x = "weeks", y = NULL) + theme_ida()
  })

  output$cmp_tbl <- renderTable({
    cmp_data() %>% group_by(regimen) %>%
      summarise(`Hb start` = first(HB), `Hb end` = last(HB),
                dHb = last(HB) - first(HB),
                `TSAT end (%)` = last(TSAT), `ferritin end` = last(FERR),
                `stores end (mg)` = last(STORES),
                `retic peak (1e9/L)` = max(RET_ABS),
                `phosphate nadir` = min(PHOS),
                `MCH end (pg)` = last(MCH), .groups = "drop")
  }, digits = 3)
}

shinyApp(ui, server)
