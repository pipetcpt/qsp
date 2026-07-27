## =====================================================================
##  Cancer Anorexia-Cachexia Syndrome (CACS) — QSP Shiny dashboard
##  cacs_shiny_app.R
## ---------------------------------------------------------------------
##  Front end for cacs_mrgsolve_model.R.
##
##  The dashboard is organised around the model's central claim: that mass
##  and function are separate endpoints driven by separate arms, and that
##  a drug's position in the mechanism determines which one it can move.
##  Every tab is there to make one part of that argument inspectable.
##
##    1  환자 프로파일   Patient       tumour, host, staging, the driver
##    2  약동학          PK            all eleven agents, one panel
##    3  ARM A 식욕      Intake        GDF-15 -> GFRAL -> melanocortin -> kcal
##    4  ARM B 이화      Catabolism    STAT3/SMAD -> FoxO -> UPS, and REE
##    5  AXIS C 근질     Quality       the axis mass-only drugs never touch
##    6  체성분          Composition   muscle, fat, water, and what DXA sees
##    7  임상 엔드포인트 Endpoints     weight, grip, FAACT, ECOG, staging
##    8  시나리오 비교   Compare       run several arms side by side
##    9  바이오마커      Biomarkers    CRP, albumin, mGPS, IGF-1, cytokines
##   10  해리 & 이력현상 Dissociation  the LBM-vs-function grid, and timing
##
##  Run:
##    setwd("cancer-cachexia"); shiny::runApp("cacs_shiny_app.R")
##  Requires: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
## =====================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(DT)
})

## ---- load the model without executing its demonstration block ---------
model_env <- new.env()
local({
  src <- readLines("cacs_mrgsolve_model.R")
  cut <- grep("^if \\(identical\\(environment\\(\\), globalenv\\(\\)\\)\\)", src)
  if (length(cut)) src <- src[seq_len(cut[1] - 1)]
  eval(parse(text = paste(src, collapse = "\n")), envir = model_env)
})
mod        <- model_env$mod
CMT        <- model_env$CMT
pat_nsclc  <- model_env$pat_nsclc
pat_panc   <- model_env$pat_panc
pat_trial  <- model_env$pat_trial
pat_gdfhi  <- model_env$pat_gdfhi
pat_early  <- model_env$pat_early
pat_refr   <- model_env$pat_refr

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

PAL <- c("#2f5f9e", "#b04a4a", "#2f7a3f", "#e07b1f", "#7d3cb0",
         "#c2568b", "#3a7f86", "#8a6b2f", "#5c6670", "#a03c6a")

DRUGS <- list(
  anamorelin    = list(label = "Anamorelin 100 mg PO qd",        cmt = "ANAD", dose = 100,  F1 = 0.35, ii = 1),
  megestrol     = list(label = "Megestrol acetate 800 mg/d",     cmt = "MEGD", dose = 800,  F1 = 0.30, ii = 1),
  dexamethasone = list(label = "Dexamethasone 4 mg/d",           cmt = "DEXD", dose = 4,    F1 = 0.80, ii = 1),
  olanzapine    = list(label = "Olanzapine 2.5 mg qhs",          cmt = "OLZD", dose = 2.5,  F1 = 0.60, ii = 1),
  enobosarm     = list(label = "Enobosarm 3 mg PO qd",           cmt = "ENOD", dose = 3,    F1 = 0.60, ii = 1),
  espindolol    = list(label = "Espindolol 10 mg PO bid",        cmt = "ESPD", dose = 10,   F1 = 0.85, ii = 0.5),
  epa           = list(label = "EPA 2 g/d",                      cmt = "EPAD", dose = 2000, F1 = 0.90, ii = 1),
  celecoxib     = list(label = "Celecoxib 200 mg bid",           cmt = "CELD", dose = 200,  F1 = 0.40, ii = 0.5)
)

build_events <- function(sel, days, pons_mg, tcz_on, bim_on) {
  parts <- list()
  for (nm in sel) {
    d <- DRUGS[[nm]]
    if (is.null(d)) next
    parts[[length(parts) + 1]] <- ev(amt = d$dose * d$F1, cmt = CMT(d$cmt),
                                     ii = d$ii, addl = ceiling(days / d$ii) - 1, time = 0)
  }
  if (pons_mg > 0)
    parts[[length(parts) + 1]] <- ev(amt = pons_mg / 147000 * 1e6 * 0.65,
                                     cmt = CMT("POND"), ii = 28,
                                     addl = floor(days / 28), time = 0)
  if (isTRUE(tcz_on))
    parts[[length(parts) + 1]] <- ev(amt = 162 * 0.80, cmt = CMT("TCZD"),
                                     ii = 7, addl = floor(days / 7), time = 0)
  if (isTRUE(bim_on))
    parts[[length(parts) + 1]] <- ev(amt = 10 * 70, cmt = CMT("BIMC"),
                                     ii = 28, addl = floor(days / 28), time = 0)
  if (!length(parts)) return(NULL)
  Reduce(`+`, parts)
}

simulate <- function(patient, pars, events, days) {
  m <- param(mod, modifyList(patient, pars))
  tg <- c(seq(0, 14, 0.25), seq(15, days, 1))
  if (is.null(events)) mrgsim_df(m, end = days, add = tg, delta = 1)
  else mrgsim_df(m, events = events, end = days, add = tg, delta = 1)
}

at_day <- function(d, t) d[which.min(abs(d$time - t)), ]

long <- function(d, vars, labels = vars) {
  out <- d[, c("time", vars), drop = FALSE]
  out <- pivot_longer(out, -time, names_to = "var", values_to = "value")
  out$var <- factor(out$var, levels = vars, labels = labels)
  out
}

facet_plot <- function(d, vars, labels = vars, ncol = 3, ylab = NULL) {
  ggplot(long(d, vars, labels), aes(time, value)) +
    geom_line(linewidth = 0.8, colour = PAL[1]) +
    facet_wrap(~var, scales = "free_y", ncol = ncol) +
    labs(x = "Day", y = ylab) + THEME
}

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel("Cancer Anorexia-Cachexia Syndrome — QSP 시뮬레이터 (CACS)"),
  tags$p(style = "color:#555;margin-top:-10px",
         HTML("<b>ARM A</b> 식욕 (GDF-15 &rarr; GFRAL &rarr; 멜라노코르틴) &middot; ",
              "<b>ARM B</b> 이화 + 과대사 (STAT3/SMAD &rarr; FoxO) &middot; ",
              "<b>AXIS C</b> 근육의 질 (force per kg) &nbsp;&mdash;&nbsp; ",
              "GRIP = MASS &times; QUALITY")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("patient", "환자 유형 (Patient archetype)",
                  c("진행성 비소세포폐암 (NSCLC)"       = "nsclc",
                    "췌장암 - 고카켁신 (pancreatic)"     = "panc",
                    "전카켁시아 조기 (precachexia)"      = "early",
                    "불응성 카켁시아 (refractory)"       = "refr",
                    "등록시험 집단 (trial-like)"         = "trial",
                    "GDF-15 고농도 집단 (enriched)"      = "gdfhi",
                    "건강 대조 (no tumour)"              = "healthy"),
                  selected = "nsclc"),
      sliderInput("days", "시뮬레이션 기간 (days)", 28, 540, 168, step = 14),
      hr(),
      h5("종양 (Tumour)"),
      sliderInput("tvol", "초기 종양 부하 (g)", 0, 400, 25, step = 5),
      sliderInput("kgrow", "Gompertz 성장률 (1/d)", 0.001, 0.020, 0.006, step = 0.001),
      sliderInput("fcach", "조직형 카켁신 유발도", 0.3, 2.0, 1.0, step = 0.1),
      sliderInput("cxsens", "숙주 감수성 CXSENS", 0.5, 1.8, 1.0, step = 0.05),
      hr(),
      h5("항암치료 (Anticancer therapy)"),
      sliderInput("onceff", "종양 사멸률 (1/d)", 0, 0.08, 0, step = 0.005),
      sliderInput("oncstart", "치료 시작일", 0, 400, 14, step = 7),
      sliderInput("myotox", "직접 근독성 (0-1)", 0, 1, 0, step = 0.05),
      hr(),
      h5("약물 (Drugs)"),
      checkboxGroupInput("drugs", NULL,
                         choices = setNames(names(DRUGS), sapply(DRUGS, `[[`, "label"))),
      selectInput("pons", "Ponsegromab SC q4w",
                  c("없음" = "0", "100 mg" = "100", "200 mg" = "200", "400 mg" = "400"),
                  selected = "0"),
      checkboxInput("tcz", "Tocilizumab 162 mg SC weekly", FALSE),
      checkboxInput("bim", "Bimagrumab 10 mg/kg IV q4w", FALSE),
      hr(),
      h5("비약물 (Non-pharmacologic)"),
      sliderInput("ons", "ONS 추가 열량 (kcal/d)", 0, 900, 0, step = 100),
      sliderInput("prot", "단백 처방 (g/kg/d)", 0.6, 2.0, 1.0, step = 0.1),
      sliderInput("exres", "저항운동 (0-1)", 0, 1, 0, step = 0.1),
      sliderInput("exaer", "유산소운동 (0-1)", 0, 1, 0, step = 0.1),
      checkboxInput("leu", "류신/HMB 강화 볼루스", FALSE),
      hr(),
      h5("증상 부담 (Symptom burden)"),
      sliderInput("gibar", "기계적/점막 장벽 GIBAR", 0, 0.6, 0, step = 0.05)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. 환자 프로파일",
                 br(), htmlOutput("summary_box"), hr(),
                 plotOutput("p_driver", height = "420px"),
                 hr(), DTOutput("t_stage")),
        tabPanel("2. 약동학 (PK)",
                 br(), plotOutput("p_pk", height = "560px"),
                 helpText("모든 농도는 선택한 약물에 대해서만 표시됩니다. ",
                          "단클론항체 청소율은 알부민 의존적 FcRn 재순환을 통해 ",
                          "저알부민혈증에서 증가합니다 — 가장 아픈 환자가 가장 적은 노출을 받습니다.")),
        tabPanel("3. ARM A — 식욕/섭취",
                 br(), plotOutput("p_armA", height = "560px"),
                 helpText("렙틴 하강은 정상이라면 AgRP를 켜야 하지만, 중추 IL-1b/PGE2가 ",
                          "그 보상을 차단합니다. 이것이 굶주림과 카켁시아의 차이입니다.")),
        tabPanel("4. ARM B — 이화/과대사",
                 br(), plotOutput("p_armB", height = "560px"),
                 helpText("REE/pREE > 1.10 이 과대사입니다. 섭취가 줄어드는 동안 ",
                          "소비가 늘어나는 것이 이 질환의 정의적 역설입니다.")),
        tabPanel("5. AXIS C — 근육의 질",
                 br(), plotOutput("p_axisC", height = "520px"),
                 hr(), htmlOutput("qual_note")),
        tabPanel("6. 체성분",
                 br(), plotOutput("p_body", height = "520px"),
                 hr(), htmlOutput("dxa_note")),
        tabPanel("7. 임상 엔드포인트",
                 br(), plotOutput("p_end", height = "560px"),
                 hr(), DTOutput("t_end")),
        tabPanel("8. 시나리오 비교",
                 br(),
                 checkboxGroupInput("cmp", "비교할 시나리오",
                                    choices = c("무치료 자연경과" = "none",
                                                "영양 단독" = "nut",
                                                "영양 + 저항운동" = "nutex",
                                                "Anamorelin" = "ana",
                                                "Ponsegromab 400 mg" = "pons",
                                                "Megestrol" = "meg",
                                                "Enobosarm" = "eno",
                                                "Tocilizumab" = "tcz",
                                                "Espindolol" = "esp",
                                                "다중모드 (multimodal)" = "multi"),
                                    selected = c("none", "nut", "ana", "pons", "multi"),
                                    inline = TRUE),
                 plotOutput("p_cmp", height = "560px"),
                 hr(), DTOutput("t_cmp")),
        tabPanel("9. 바이오마커",
                 br(), plotOutput("p_bio", height = "560px"),
                 hr(), htmlOutput("tcz_note")),
        tabPanel("10. 해리 & 이력현상",
                 br(), h4("제지방량과 기능의 해리 (the LBM-function dissociation)"),
                 plotOutput("p_dissoc", height = "380px"),
                 DTOutput("t_dissoc"),
                 hr(), h4("돌아올 수 없는 지점 (the point of no return)"),
                 helpText("동일한 완치적 항암치료를 점점 늦게 시작했을 때 ",
                          "540일째 회복되는 근육량."),
                 plotOutput("p_hyst", height = "340px"))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  base_patient <- reactive({
    switch(input$patient,
           nsclc = pat_nsclc, panc = pat_panc, early = pat_early,
           refr = pat_refr, trial = pat_trial, gdfhi = pat_gdfhi,
           healthy = list(TVOL0 = 0, KGROW = 0.006, FCACHEX = 1, CXSENS = 1))
  })

  observeEvent(input$patient, {
    p <- base_patient()
    updateSliderInput(session, "tvol",   value = p$TVOL0)
    updateSliderInput(session, "kgrow",  value = p$KGROW)
    updateSliderInput(session, "fcach",  value = p$FCACHEX)
    updateSliderInput(session, "cxsens", value = p$CXSENS)
  })

  pars <- reactive({
    list(TVOL0 = input$tvol, KGROW = input$kgrow, FCACHEX = input$fcach,
         CXSENS = input$cxsens, ONCEFF = input$onceff,
         ONCSTART = if (input$onceff > 0) input$oncstart else 1e6,
         MYOTOX = input$myotox, ONSKCAL = input$ons, PROTTG = input$prot,
         EXRES = input$exres, EXAER = input$exaer,
         LEUBOL = as.numeric(isTRUE(input$leu)), GIBAR = input$gibar)
  })

  sim <- reactive({
    ev0 <- build_events(input$drugs, input$days,
                        as.numeric(input$pons), input$tcz, input$bim)
    simulate(list(), pars(), ev0, input$days)
  })

  ## ---- 1. patient profile --------------------------------------------
  output$summary_box <- renderUI({
    d <- sim(); a <- at_day(d, 0); b <- at_day(d, input$days)
    card <- function(lab, val, sub = "")
      sprintf("<div style='display:inline-block;min-width:150px;margin:4px 10px 4px 0;
               padding:8px 12px;border:1px solid #dcdcdc;border-radius:8px'>
               <div style='font-size:11px;color:#666'>%s</div>
               <div style='font-size:19px;font-weight:600'>%s</div>
               <div style='font-size:11px;color:#888'>%s</div></div>", lab, val, sub)
    stage_lab <- c("없음", "전카켁시아", "카켁시아", "불응성")[b$STAGE + 1]
    HTML(paste0(
      card("체중 (kg)", sprintf("%.1f", b$BW), sprintf("기저 %.1f", a$BW)),
      card("6개월 체중감소", sprintf("%.1f %%", b$PCTWL), sprintf("BMI %.1f", b$BMI)),
      card("제지방량 LBM (kg)", sprintf("%.1f", b$LBM), sprintf("골격근 %.1f kg", b$MUSC)),
      card("악력 (kg)", sprintf("%.1f", b$GRIP), sprintf("근질 %.2f", b$MYOQ)),
      card("섭취 (kcal/d)", sprintf("%.0f", b$INTK), sprintf("TEE %.0f", b$TEE)),
      card("REE/pREE", sprintf("%.2f", b$RQREE),
           if (b$RQREE > 1.10) "과대사 (hypermetabolic)" else "정상범위"),
      card("Fearon 병기", stage_lab, sprintf("Martin grade %d", as.integer(b$MGRADE))),
      card("ECOG", sprintf("%.1f", b$ECOG), sprintf("FAACT %.0f/48", b$FAACT)),
      card("생존확률", sprintf("%.0f %%", 100 * b$SURV), sprintf("%d일 시점", input$days))
    ))
  })

  output$p_driver <- renderPlot({
    facet_plot(sim(),
      c("TUMOR", "GDF15", "IL6", "ACTA", "CRP", "ALB"),
      c("종양 부하 (g)", "총 GDF-15 (ng/mL)", "IL-6 (pg/mL)",
        "액티빈 A (pg/mL)", "CRP (mg/L)", "알부민 (g/L)"))
  })

  output$t_stage <- renderDT({
    d <- sim()
    ts <- unique(pmin(c(0, 28, 56, 84, 112, 168, 252, 365, 540), input$days))
    tab <- do.call(rbind, lapply(ts, function(t) {
      b <- at_day(d, t)
      data.frame(Day = t, BW = round(b$BW, 1), `%WL` = round(b$PCTWL, 1),
                 BMI = round(b$BMI, 1), SMI = round(b$SMI, 1),
                 Sarcopenia = ifelse(b$SARCO > 0.5, "yes", "no"),
                 Fearon = c("none", "pre", "cachexia", "refractory")[b$STAGE + 1],
                 Martin = as.integer(b$MGRADE), mGPS = as.integer(b$MGPS),
                 ECOG = round(b$ECOG, 2), Survival = round(b$SURV, 3),
                 check.names = FALSE)
    }))
    datatable(tab, rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  ## ---- 2. PK -----------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim()
    v <- c(); l <- c()
    add <- function(cond, var, lab) { if (cond) { v <<- c(v, var); l <<- c(l, lab) } }
    add("anamorelin"    %in% input$drugs, "CANA", "Anamorelin (ng/mL)")
    add("megestrol"     %in% input$drugs, "CMEG", "Megestrol (ng/mL)")
    add("dexamethasone" %in% input$drugs, "CDEX", "Dexamethasone (ng/mL)")
    add("olanzapine"    %in% input$drugs, "COLZ", "Olanzapine (ng/mL)")
    add("enobosarm"     %in% input$drugs, "CENO", "Enobosarm (ng/mL)")
    add("espindolol"    %in% input$drugs, "CESP", "Espindolol (ng/mL)")
    add(as.numeric(input$pons) > 0,       "CPON", "Ponsegromab (nM)")
    add(isTRUE(input$tcz),                "CTCZ", "Tocilizumab (mg/L)")
    add(isTRUE(input$bim),                "CBIM", "Bimagrumab (mg/L)")
    if (!length(v))
      return(ggplot() + annotate("text", 0, 0, label = "약물을 선택하세요", size = 6) +
               theme_void())
    facet_plot(d, v, l, ncol = 3, ylab = "Concentration")
  })

  ## ---- 3. ARM A ---------------------------------------------------------
  output$p_armA <- renderPlot({
    facet_plot(sim(),
      c("GDFFREE", "BSSIG", "NAUS", "AGRP", "POMC", "ANOR", "LEPT", "GHRL", "INTK"),
      c("유리 GDF-15 (ng/mL)", "뇌간 혐오 신호", "오심 (0-10)",
        "AgRP/NPY 긴장도", "POMC 긴장도", "통합 식욕부진 구동",
        "렙틴 (ng/mL)", "그렐린 (pg/mL)", "에너지 섭취 (kcal/d)"))
  })

  ## ---- 4. ARM B ---------------------------------------------------------
  output$p_armB <- renderPlot({
    facet_plot(sim(),
      c("STAT3", "SMAD", "FOXO", "UPS", "AUTOP", "ARES",
        "IGF1", "CORT", "RQREE"),
      c("근육 STAT3", "SMAD2/3", "FoxO 활성",
        "프로테아좀 flux", "오토파지 flux", "동화 저항 지수",
        "IGF-1 (ng/mL)", "코르티솔 (ug/dL)", "REE / 예측REE"))
  })

  ## ---- 5. AXIS C --------------------------------------------------------
  output$p_axisC <- renderPlot({
    d <- sim()
    facet_plot(d, c("PGC1", "MITO", "ROS", "MYOQ", "SATC", "GRIP"),
               c("PGC-1alpha", "미토콘드리아 함량", "산화 스트레스 (ROS)",
                 "근육의 질 MYOQ (force/kg)", "위성세포 풀 SATC",
                 "악력 (kg) = 질 x 수축성 질량"), ncol = 3)
  })

  output$qual_note <- renderUI({
    d <- sim(); b <- at_day(d, input$days)
    HTML(sprintf(
      "<div style='background:#f5f0fa;padding:12px;border-radius:8px'>
       <b>측정되는 질량과 당기는 질량 (mass DXA sees vs mass that pulls)</b><br/>
       골격근 <b>%.2f kg</b> 중 수축성 질량은 <b>%.2f kg</b>, GH 유래 제지방 수분은
       <b>%.2f kg</b>입니다. DXA가 제지방량으로 세는 것 가운데
       <b>%.2f kg</b>은 힘을 내지 않습니다.<br/>
       근육의 질 MYOQ = <b>%.3f</b> — 이 축은 염증이 내리고 부하(저항운동)만이 올립니다.
       순수 동화 약물은 이 값을 건드리지 못하며, 그래서 제지방량 지표를 맞히고도
       기능 지표를 놓칩니다.</div>",
      b$MUSC, b$MUSC - (b$NONCON - b$LNW), b$LNW, b$NONCON, b$MYOQ))
  })

  ## ---- 6. body composition ---------------------------------------------
  output$p_body <- renderPlot({
    d <- sim()
    ggplot(long(d, c("MUSC", "FATM", "OLBM", "LIVM", "LNW", "EDEM"),
                c("골격근", "지방", "기타 제지방", "간", "GH 제지방 수분", "부종")),
           aes(time, value, fill = var)) +
      geom_area(alpha = 0.85) +
      scale_fill_manual(values = PAL) +
      labs(x = "Day", y = "kg", fill = NULL,
           title = "체성분 구성 — 쌓아 올린 각 층이 체중계의 숫자를 만듭니다") +
      THEME
  })

  output$dxa_note <- renderUI({
    d <- sim(); a <- at_day(d, 0); b <- at_day(d, input$days)
    HTML(sprintf(
      "<div style='background:#f0f4f8;padding:12px;border-radius:8px'>
       체중 변화 <b>%+.2f kg</b> = 근육 <b>%+.2f</b> + 지방 <b>%+.2f</b> +
       기타 제지방 <b>%+.2f</b> + 간 <b>%+.2f</b> + 제지방 수분 <b>%+.2f</b> +
       부종 <b>%+.2f</b>.<br/>
       저알부민혈증과 프로게스틴은 체액을 늘려 <b>실제 조직 손실을 체중계에서 감춥니다</b>.
       체중이 안정적이라는 것이 조직이 안정적이라는 뜻은 아닙니다.</div>",
      b$BW - a$BW, b$MUSC - a$MUSC, b$FATM - a$FATM, b$OLBM - a$OLBM,
      b$LIVM - a$LIVM, b$LNW - a$LNW, b$EDEM - a$EDEM))
  })

  ## ---- 7. endpoints ------------------------------------------------------
  output$p_end <- renderPlot({
    facet_plot(sim(),
      c("BW", "PCTWL", "LBM", "GRIP", "SMI", "FAACT", "ECOG", "MGRADE", "SURV"),
      c("체중 (kg)", "6개월 체중감소 (%)", "제지방량 (kg)", "악력 (kg)",
        "L3 SMI (cm2/m2)", "FAACT A/CS (0-48)", "ECOG PS", "Martin grade",
        "생존확률"))
  })

  output$t_end <- renderDT({
    d <- sim(); a <- at_day(d, 0)
    ts <- unique(pmin(c(28, 56, 84, 112, 168, 252, 365), input$days))
    tab <- do.call(rbind, lapply(ts, function(t) {
      b <- at_day(d, t)
      data.frame(Day = t,
                 `dBW` = round(b$BW - a$BW, 2),
                 `dLBM` = round(b$LBM - a$LBM, 2),
                 `dMuscle` = round(b$MUSC - a$MUSC, 2),
                 `dFat` = round(b$FATM - a$FATM, 2),
                 `dGrip` = round(b$GRIP - a$GRIP, 2),
                 `dFAACT` = round(b$FAACT - a$FAACT, 1),
                 Intake = round(b$INTK), TEE = round(b$TEE),
                 `REE ratio` = round(b$RQREE, 3),
                 check.names = FALSE)
    }))
    datatable(tab, rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  ## ---- 8. scenario comparison ---------------------------------------------
  cmp_defs <- list(
    none  = list(lab = "무치료",             pars = list(), ev = function(dy) NULL),
    nut   = list(lab = "영양 단독",          pars = list(ONSKCAL = 600, PROTTG = 1.5),
                 ev = function(dy) NULL),
    nutex = list(lab = "영양 + 저항운동",    pars = list(ONSKCAL = 600, PROTTG = 1.5,
                                                        EXRES = 1, EXAER = 0.6, LEUBOL = 1),
                 ev = function(dy) NULL),
    ana   = list(lab = "Anamorelin",         pars = list(),
                 ev = function(dy) build_events("anamorelin", dy, 0, FALSE, FALSE)),
    pons  = list(lab = "Ponsegromab 400 mg", pars = list(),
                 ev = function(dy) build_events(character(0), dy, 400, FALSE, FALSE)),
    meg   = list(lab = "Megestrol",          pars = list(),
                 ev = function(dy) build_events("megestrol", dy, 0, FALSE, FALSE)),
    eno   = list(lab = "Enobosarm",          pars = list(),
                 ev = function(dy) build_events("enobosarm", dy, 0, FALSE, FALSE)),
    tcz   = list(lab = "Tocilizumab",        pars = list(),
                 ev = function(dy) build_events(character(0), dy, 0, TRUE, FALSE)),
    esp   = list(lab = "Espindolol",         pars = list(),
                 ev = function(dy) build_events("espindolol", dy, 0, FALSE, FALSE)),
    multi = list(lab = "다중모드",           pars = list(ONSKCAL = 600, PROTTG = 1.5,
                                                        EXRES = 1, EXAER = 0.8, LEUBOL = 1),
                 ev = function(dy) build_events("epa", dy, 400, FALSE, FALSE))
  )

  cmp_runs <- reactive({
    req(length(input$cmp) > 0)
    base <- pars()
    do.call(rbind, lapply(input$cmp, function(k) {
      cd <- cmp_defs[[k]]
      o <- simulate(list(), modifyList(base, cd$pars), cd$ev(input$days), input$days)
      o$arm <- cd$lab
      o
    }))
  })

  output$p_cmp <- renderPlot({
    d <- cmp_runs()
    v <- c("BW", "LBM", "MUSC", "GRIP", "INTK", "FAACT", "MYOQ", "ECOG", "SURV")
    l <- c("체중 (kg)", "제지방량 (kg)", "골격근 (kg)", "악력 (kg)",
           "섭취 (kcal/d)", "FAACT", "근질 MYOQ", "ECOG", "생존확률")
    dd <- pivot_longer(d[, c("time", "arm", v)], -c(time, arm),
                       names_to = "var", values_to = "value")
    dd$var <- factor(dd$var, levels = v, labels = l)
    ggplot(dd, aes(time, value, colour = arm)) +
      geom_line(linewidth = 0.85) +
      facet_wrap(~var, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = PAL) +
      labs(x = "Day", y = NULL, colour = NULL) + THEME
  })

  output$t_cmp <- renderDT({
    d <- cmp_runs()
    tab <- do.call(rbind, lapply(split(d, d$arm), function(x) {
      a <- at_day(x, 0); b <- at_day(x, input$days)
      data.frame(Arm = x$arm[1],
                 dBW = round(b$BW - a$BW, 2), dLBM = round(b$LBM - a$LBM, 2),
                 dMuscle = round(b$MUSC - a$MUSC, 2), dFat = round(b$FATM - a$FATM, 2),
                 dGrip = round(b$GRIP - a$GRIP, 2), MYOQ = round(b$MYOQ, 3),
                 Intake = round(b$INTK), `REE ratio` = round(b$RQREE, 3),
                 CRP = round(b$CRP, 1), ECOG = round(b$ECOG, 2),
                 Survival = round(b$SURV, 3), check.names = FALSE)
    }))
    datatable(tab, rownames = FALSE, options = list(dom = "t", pageLength = 12))
  })

  ## ---- 9. biomarkers -------------------------------------------------------
  output$p_bio <- renderPlot({
    facet_plot(sim(),
      c("CRP", "ALB", "MGPS", "IL6", "TNFA", "MSTN", "IGF1", "TEST", "GDFTOT"),
      c("CRP (mg/L)", "알부민 (g/L)", "mGPS (0-2)", "IL-6 (pg/mL)",
        "TNF-alpha (pg/mL)", "마이오스타틴 (ng/mL) - 역설적 하강",
        "IGF-1 (ng/mL)", "테스토스테론 (ng/dL)", "총 GDF-15 (ng/mL)"))
  })

  output$tcz_note <- renderUI({
    HTML("<div style='background:#fdf6ec;padding:12px;border-radius:8px'>
          <b>두 가지 함정 (two traps in this panel)</b><br/>
          &bull; <b>토실리주맙에서 혈장 IL-6는 올라갑니다.</b> IL-6 수용체를 막으면
          수용체 매개 청소가 사라져 리간드가 축적됩니다. IL-6 상승은 치료 실패가
          아니라 표적 결합의 증거이며, 실제로 판단해야 할 것은 CRP입니다.<br/>
          &bull; <b>마이오스타틴은 카켁시아에서 내려갑니다.</b> 근육 자체가 생산원이기
          때문입니다. 근위축의 원인이 아니라 결과이므로 바이오마커로 오독하면 안 됩니다.
          이화를 끄는 ActRIIB 리간드는 액티빈 A입니다.</div>")
  })

  ## ---- 10. dissociation and hysteresis --------------------------------------
  dissoc <- reactive({
    base <- pars()
    arms <- list(
      list(k = "Anamorelin",   p = list(), e = function(dy) build_events("anamorelin", dy, 0, FALSE, FALSE)),
      list(k = "Enobosarm",    p = list(), e = function(dy) build_events("enobosarm", dy, 0, FALSE, FALSE)),
      list(k = "Bimagrumab",   p = list(), e = function(dy) build_events(character(0), dy, 0, FALSE, TRUE)),
      list(k = "Ponsegromab",  p = list(), e = function(dy) build_events(character(0), dy, 400, FALSE, FALSE)),
      list(k = "Megestrol",    p = list(), e = function(dy) build_events("megestrol", dy, 0, FALSE, FALSE)),
      list(k = "Espindolol",   p = list(), e = function(dy) build_events("espindolol", dy, 0, FALSE, FALSE)),
      list(k = "저항운동",     p = list(EXRES = 1, EXAER = 0.6), e = function(dy) NULL),
      list(k = "다중모드",     p = list(ONSKCAL = 600, PROTTG = 1.5, EXRES = 1, EXAER = 0.8),
           e = function(dy) build_events("epa", dy, 400, FALSE, FALSE))
    )
    ref <- simulate(list(), base, NULL, 84)
    r0 <- at_day(ref, 0); r8 <- at_day(ref, 84)
    do.call(rbind, lapply(arms, function(a) {
      o <- simulate(list(), modifyList(base, a$p), a$e(84), 84)
      x0 <- at_day(o, 0); x8 <- at_day(o, 84)
      dl <- (x8$LBM - x0$LBM) - (r8$LBM - r0$LBM)
      dg <- (x8$GRIP - x0$GRIP) - (r8$GRIP - r0$GRIP)
      data.frame(Intervention = a$k, dLBM = round(dl, 2), dGrip = round(dg, 2),
                 `Grip per kg LBM` = ifelse(abs(dl) < 0.05, NA, round(dg / dl, 2)),
                 `% of pure-mass expectation` =
                   ifelse(abs(dl) < 0.05, NA, round(100 * (dg / dl) / (40 / 25))),
                 check.names = FALSE)
    }))
  })

  output$p_dissoc <- renderPlot({
    d <- dissoc()
    ggplot(d, aes(dLBM, dGrip, label = Intervention)) +
      geom_abline(slope = 40 / 25, intercept = 0, linetype = "dashed", colour = "#999") +
      annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.6, size = 3.4,
               colour = "#777", label = "점선 = 늘어난 질량이 모두 정상 근육일 때의 기대선") +
      geom_hline(yintercept = 0, colour = "#ccc") +
      geom_vline(xintercept = 0, colour = "#ccc") +
      geom_point(size = 3.4, colour = PAL[1]) +
      geom_text(vjust = -0.9, size = 3.6) +
      labs(x = "12주 제지방량 변화, 무치료 대비 (kg)",
           y = "12주 악력 변화, 무치료 대비 (kg)") +
      expand_limits(x = c(-1.5, 4), y = c(-2, 9)) + THEME
  })

  output$t_dissoc <- renderDT(
    datatable(dissoc(), rownames = FALSE, options = list(dom = "t", pageLength = 10)))

  output$p_hyst <- renderPlot({
    base <- pars()
    delays <- c(0, 60, 120, 180, 240, 300)
    res <- do.call(rbind, lapply(delays, function(dl) {
      o <- simulate(list(), modifyList(base, list(ONCEFF = 0.045, ONCSTART = dl)),
                    NULL, 540)
      b <- o[nrow(o), ]
      data.frame(start = dl, muscle = b$MUSC, satc = b$SATC, grip = b$GRIP,
                 pct = 100 * b$MUSC / o$MUSC[1])
    }))
    ggplot(res, aes(start, pct)) +
      geom_line(linewidth = 1, colour = PAL[2]) +
      geom_point(size = 3, colour = PAL[2]) +
      geom_text(aes(label = sprintf("SATC %.2f", satc)), vjust = -1.1, size = 3.3,
                colour = "#666") +
      labs(x = "완치적 항암치료 시작일 (day)",
           y = "540일째 근육량 (발병 전 대비 %)",
           title = "같은 치료, 다른 시점 — 위성세포 풀이 고갈되면 되돌릴 수 없습니다") +
      expand_limits(y = c(60, 105)) + THEME
  })
}

shinyApp(ui, server)
