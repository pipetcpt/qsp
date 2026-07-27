## ============================================================================
##  Acne vulgaris (ACN) — QSP Interactive Dashboard
##  ============================================================================
##  Front-end for acn_mrgsolve_model.R.
##
##  Run:
##      shiny::runApp("acn_shiny_app.R")
##  or from this directory:
##      R -e 'shiny::runApp(".", launch.browser = TRUE)'
##
##  Ten tabs:
##      1. 환자 프로파일     — phenotype, androgen/metabolic set-point, exposome
##      2. 치료 설계         — build a regimen (topical / systemic / hormonal)
##      3. PK                — every drug concentration on one time axis
##      4. 피지선·안드로겐   — SGM, SER, AR signal, SHBG/FAI, IGF-1, mTORC1
##      5. 미생물·항생제 내성— C. acnes planktonic/biofilm and the resistant
##                             fraction; the stewardship argument in one plot
##      6. 염증 캐스케이드   — TLR2 -> IL-1beta -> IL-8 -> neutrophils -> MMP
##      7. 병변 수·엔드포인트— the lesion-transit chain and IGA
##      8. 시나리오 비교     — the 17 prebuilt scenarios side by side
##      9. 안전성            — TG, ALT, K+, mucocutaneous, cumulative mg/kg
##     10. 흉터·색소침착     — PIH and atrophic scarring, and what prevents them
##
##  DISCLAIMER: educational / research tool. Not for clinical use.
## ============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(DT)
})

source("acn_mrgsolve_model.R", local = TRUE, chdir = TRUE)

MOD <- ACN_build()

## ---------------------------------------------------------------------------
##  Plot furniture
## ---------------------------------------------------------------------------
theme_acn <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text      = element_text(face = "bold"),
          legend.position = "bottom",
          legend.title    = element_blank(),
          plot.title      = element_text(face = "bold"))
}

PAL <- c("#C0392B", "#2874A6", "#1E8449", "#B9770E", "#7D3C98",
         "#117A65", "#CA6F1E", "#5D6D7E", "#A93226", "#1F618D")

long_plot <- function(d, vars, labels = vars, ylab = "", title = "",
                      hline = NULL, logy = FALSE) {
  dd <- d[, c("week", vars), drop = FALSE]
  names(dd) <- c("week", labels)
  dd <- tidyr::pivot_longer(dd, -week, names_to = "v", values_to = "y")
  dd$v <- factor(dd$v, levels = labels)
  p <- ggplot(dd, aes(week, y, colour = v)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = rep(PAL, 4)) +
    labs(x = "주 (weeks)", y = ylab, title = title) +
    theme_acn()
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = 2,
                                           colour = "grey40")
  if (logy) p <- p + scale_y_log10()
  p
}

## ---------------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("심상성 여드름 QSP 대시보드 — Acne vulgaris QSP dashboard"),
  tags$p(style = "color:#7B241C;font-weight:bold;",
         "교육·연구용 모델입니다. 임상 의사결정에 사용하지 마십시오. / ",
         "Educational model. Not for clinical use."),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("1. 환자 표현형"),
      selectInput("pheno", "표현형 (phenotype)",
                  choices = c("경증 면포성 (mild comedonal)"      = "mild_comedonal",
                              "중등증 (moderate)"                 = "moderate",
                              "중등증 남성 (moderate male)"       = "moderate_male",
                              "중증 결절성 (severe nodular)"      = "severe_nodular",
                              "성인 여성 (adult female)"          = "adult_female",
                              "PCOS 동반"                         = "pcos",
                              "짙은 피부색 (skin of colour)"      = "skin_of_colour"),
                  selected = "moderate"),
      sliderInput("SEVX",   "체질적 피지 구동 SEVX", 1.0, 2.4, 1.55, 0.05),
      sliderInput("NODPROP","결절 성향 NODPROP",     0.1, 3.0, 1.00, 0.1),
      sliderInput("SKINPIG","색소침착 성향 (Fitzpatrick)", 0.3, 2.4, 1.0, 0.1),

      h4("2. 노출 (Exposome)"),
      sliderInput("GLYLOAD","혈당부하 (0=저GI, 1=고GI)", 0, 1, 0.5, 0.05),
      sliderInput("DAIRY",  "유제품/웨이 섭취",          0, 1, 0.3, 0.05),
      sliderInput("STRESS", "심리적 스트레스",           0, 1, 0.3, 0.05),
      sliderInput("UVX",    "UV/산화 노출",              0, 1, 0.2, 0.05),
      sliderInput("PHYLIA", "C. acnes IA1 우세도",       0, 1, 0.75, 0.05),

      h4("3. 국소 치료"),
      selectInput("ret", "국소 레티노이드",
                  choices = c("없음" = "none", names(ACN_RETPOT)),
                  selected = "adapalene_0.1"),
      checkboxInput("bpo",   "벤조일퍼옥사이드 2.5% QD", FALSE),
      checkboxInput("cli",   "클린다마이신 1% BID",      FALSE),
      checkboxInput("aze",   "아젤라산 15% BID",         FALSE),
      checkboxInput("clasco","클라스코테론 1% BID",      FALSE),
      checkboxInput("dap",   "답손 7.5% QD",             FALSE),
      sliderInput("ADHERE",  "순응도 (adherence)", 0.2, 1.0, 1.0, 0.05),

      h4("4. 전신 치료"),
      selectInput("tet", "테트라사이클린계",
                  choices = c("없음"                        = "none",
                              "독시사이클린 100 mg QD"      = "doxy100",
                              "독시사이클린 40 mg MR (아항균)" = "sub40",
                              "미노사이클린 100 mg QD"      = "mino100",
                              "사레사이클린 1.5 mg/kg QD"   = "sare"),
                  selected = "none"),
      numericInput("tet_wk", "항생제 투여 기간 (주)", 12, 0, 52, 1),
      sliderInput("iso",  "이소트레티노인 (mg/kg/일)", 0, 1.2, 0, 0.05),
      numericInput("iso_wk", "이소트레티노인 기간 (주)", 24, 0, 60, 1),
      checkboxInput("food", "고지방식과 함께 복용", TRUE),
      checkboxInput("lidose", "Lidose/미분화 제형", FALSE),
      sliderInput("spiro","스피로놀락톤 (mg/일)", 0, 200, 0, 25),
      sliderInput("coc",  "COC 에티닐에스트라디올 (µg)", 0, 35, 0, 5),

      h4("5. 시뮬레이션"),
      numericInput("weeks", "시뮬레이션 기간 (주)", 52, 4, 160, 4),
      actionButton("go", "실행 (Run)", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1. 환자 프로파일",
          h4("기저 정상상태 (untreated steady state)"),
          p("모든 시뮬레이션은 이 환자의 자체 정상상태에서 출발합니다. ",
            "병변 수는 입력값이 아니라 안드로겐·대사·노출 설정에서 ",
            "창발(emergent)하는 값입니다."),
          DT::dataTableOutput("tbl_base"),
          hr(),
          h4("네 개의 병인 기둥 (four pathogenic pillars)"),
          plotOutput("p_pillars", height = "320px"),
          hr(),
          verbatimTextOutput("txt_pheno")),

        tabPanel("2. 치료 설계",
          h4("선택한 요법"),
          verbatimTextOutput("txt_regimen"),
          hr(),
          h4("가이드라인 점검 (guideline sanity checks)"),
          htmlOutput("txt_flags")),

        tabPanel("3. PK",
          plotOutput("p_pk_sys", height = "300px"),
          plotOutput("p_pk_top", height = "300px"),
          plotOutput("p_pk_eff", height = "280px")),

        tabPanel("4. 피지선·안드로겐",
          plotOutput("p_seb",  height = "300px"),
          plotOutput("p_horm", height = "300px"),
          plotOutput("p_met",  height = "260px")),

        tabPanel("5. 미생물·내성",
          p(strong("이 탭이 항생제 청지기(stewardship) 논증입니다."),
            " 항생제 단독요법에서는 내성 분획(RESF)이 상승하며 균 부하가 ",
            "다시 올라오고 병변이 반등합니다. BPO를 병용하면 내성 분획이 ",
            "제거되어 반등이 사라집니다."),
          plotOutput("p_micro", height = "320px"),
          plotOutput("p_res",   height = "300px")),

        tabPanel("6. 염증 캐스케이드",
          plotOutput("p_infl", height = "320px"),
          plotOutput("p_ker",  height = "300px")),

        tabPanel("7. 병변 수·엔드포인트",
          plotOutput("p_les",  height = "340px"),
          plotOutput("p_iga",  height = "260px"),
          hr(),
          h4("엔드포인트 요약"),
          DT::dataTableOutput("tbl_end")),

        tabPanel("8. 시나리오 비교",
          checkboxGroupInput("scn", "시나리오 선택", inline = TRUE,
            choices = setNames(1:17, paste0("S", 1:17)),
            selected = c(2, 3, 5, 8, 9, 11)),
          numericInput("scn_wk", "평가 시점 (주)", 12, 1, 80, 1),
          actionButton("go_scn", "시나리오 실행", class = "btn-primary"),
          hr(),
          DT::dataTableOutput("tbl_scn"),
          plotOutput("p_scn", height = "360px")),

        tabPanel("9. 안전성",
          plotOutput("p_saf", height = "320px"),
          plotOutput("p_cum", height = "260px"),
          hr(),
          htmlOutput("txt_saf")),

        tabPanel("10. 흉터·색소침착",
          plotOutput("p_scar", height = "320px"),
          hr(),
          p("위축성 흉터는 결절 부하 x MMP 활성의 시간적분으로 누적되며 ",
            "모델 안에서 되돌릴 수 없습니다. 흉터를 줄이는 유일한 방법은 ",
            "결절이 존재하는 시간을 줄이는 것 — 즉 조기의 충분한 치료입니다."),
          DT::dataTableOutput("tbl_scar")),

        tabPanel("11. 문헌·도움말",
          h4("모델 구조"),
          verbatimTextOutput("txt_help"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
##  Server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  ## ---- assemble the parameter overrides and the event table ---------------
  regimen <- reactive({
    ev  <- NULL
    lab <- character(0)
    prm <- list(ADHERE = input$ADHERE)
    wk  <- input$weeks

    if (input$ret != "none") {
      ev <- c(ev, acn_retinoid(input$ret, days = wk * 7))
      lab <- c(lab, paste("국소 레티노이드:", input$ret))
    }
    if (isTRUE(input$bpo)) {
      ev <- c(ev, acn_bpo(2.5, days = wk * 7)); lab <- c(lab, "BPO 2.5% QD")
    }
    if (isTRUE(input$cli)) {
      ev <- c(ev, acn_clindamycin(1.0, days = wk * 7)); lab <- c(lab, "클린다마이신 1% BID")
    }
    if (isTRUE(input$aze)) {
      ev <- c(ev, acn_azelaic(15, days = wk * 7)); lab <- c(lab, "아젤라산 15% BID")
    }
    if (isTRUE(input$clasco)) {
      ev <- c(ev, acn_clascoterone(1.0, days = wk * 7)); lab <- c(lab, "클라스코테론 1% BID")
    }
    if (isTRUE(input$dap)) {
      ev <- c(ev, acn_dapsone(7.5, days = wk * 7)); lab <- c(lab, "답손 7.5% QD")
    }

    if (input$tet != "none" && input$tet_wk > 0) {
      spec <- switch(input$tet,
        doxy100 = list(dose = 100, p = list()),
        sub40   = list(dose = 40,  p = list(SUBANTI = 1)),
        mino100 = list(dose = 100, p = list(CLTET = 1.6, VTET = 80)),
        sare    = list(dose = round(1.5 * 62), p = list(NARROW = 1, CLTET = 1.75, VTET = 55)))
      ev  <- c(ev, acn_tetracycline(spec$dose, days = input$tet_wk * 7))
      prm <- utils::modifyList(prm, spec$p)
      lab <- c(lab, sprintf("%s x %d주", input$tet, input$tet_wk))
    }

    if (input$iso > 0 && input$iso_wk > 0) {
      ev  <- c(ev, acn_isotretinoin(input$iso, wt = 62, days = input$iso_wk * 7))
      prm <- utils::modifyList(prm, list(FOOD = as.numeric(isTRUE(input$food)),
                                         LIDOSE = as.numeric(isTRUE(input$lidose))))
      lab <- c(lab, sprintf("이소트레티노인 %.2f mg/kg/일 x %d주 (누적 %.0f mg/kg)",
                            input$iso, input$iso_wk, input$iso * input$iso_wk * 7))
    }
    if (input$spiro > 0) {
      ev <- c(ev, acn_spironolactone(input$spiro, days = wk * 7))
      lab <- c(lab, sprintf("스피로놀락톤 %d mg/일", input$spiro))
    }
    if (input$coc > 0) {
      ev <- c(ev, acn_coc(input$coc, cycles = ceiling(wk / 4)))
      lab <- c(lab, sprintf("COC EE %d µg", input$coc))
    }
    if (!length(lab)) lab <- "치료 없음 (자연경과)"
    list(ev = ev, param = prm, label = lab)
  })

  base <- eventReactive(input$go, {
    extra <- list(SEVX = input$SEVX, NODPROP = input$NODPROP,
                  SKINPIG = input$SKINPIG, GLYLOAD = input$GLYLOAD,
                  DAIRY = input$DAIRY, STRESS = input$STRESS,
                  UVX = input$UVX, PHYLIA = input$PHYLIA)
    acn_baseline(MOD, input$pheno, extra = extra)
  }, ignoreNULL = FALSE)

  sim <- eventReactive(input$go, {
    b <- base()
    r <- regimen()
    p <- utils::modifyList(b$param,
                           list(PBOMAX = as.numeric(mrgsolve::param(MOD)$PBOMAX)))
    p <- utils::modifyList(p, r$param)
    m <- mrgsolve::update(MOD, param = p, init = b$init)
    d <- if (is.null(r$ev)) {
      mrgsolve::mrgsim_df(m, end = input$weeks * 168, delta = 12, hmax = 4)
    } else {
      mrgsolve::mrgsim_df(m, events = r$ev, end = input$weeks * 168,
                          delta = 12, hmax = 4)
    }
    d$week <- d$time / 168
    d
  }, ignoreNULL = FALSE)

  ## ---- 1. patient profile -------------------------------------------------
  output$tbl_base <- DT::renderDataTable({
    s <- base()$summary
    data.frame(
      지표 = c("염증성 병변 수", "비염증성 병변 수", "총 병변 수", "IGA (0-4)",
               "미세면포 저장고 MC", "피지 분비율 SER (µg/cm²/min)",
               "C. acnes 부하 (상대)", "내성 분획", "각질화 지수 KER",
               "AR 신호 ARS", "유리 안드로겐 지수 FAI", "IGF-1 (ng/mL)",
               "PIH 지수", "흉터 지수"),
      값 = round(c(s$INFLAM, s$NONINF, s$TOTLES, s$IGA, s$MC, s$SERO,
                   s$CACNT, s$RESF, s$KER, s$ARS, s$FAIO, s$IGF1,
                   s$PIH, s$SCAR), 2))
  }, options = list(dom = "t", pageLength = 20), rownames = FALSE)

  output$p_pillars <- renderPlot({
    s <- base()$summary
    d <- data.frame(
      pillar = factor(c("① 과피지\n(SER/기준)", "② 과각화\n(KER)",
                        "③ C. acnes\n(부하/기준)", "④ 염증\n(IL-1β)"),
                      levels = c("① 과피지\n(SER/기준)", "② 과각화\n(KER)",
                                 "③ C. acnes\n(부하/기준)", "④ 염증\n(IL-1β)")),
      value  = c(s$SERO / 0.9, s$KER, s$CACNT / 1.5, s$IL1B))
    ggplot(d, aes(pillar, value, fill = pillar)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = 1, linetype = 2) +
      scale_fill_manual(values = PAL[c(4, 3, 6, 1)], guide = "none") +
      labs(x = NULL, y = "정상 대비 배수 (fold of healthy reference)",
           title = "네 개의 병인 기둥 — 1.0 = 건강한 모낭피지선") +
      theme_acn()
  })

  output$txt_pheno <- renderPrint({
    p <- base()$param
    cat("표현형 파라미터 (phenotype parameters)\n")
    str(p[order(names(p))], give.attr = FALSE)
  })

  ## ---- 2. regimen ---------------------------------------------------------
  output$txt_regimen <- renderPrint({ cat(paste(regimen()$label, collapse = "\n")) })

  output$txt_flags <- renderUI({
    r  <- regimen(); msg <- character(0)
    abx_top <- isTRUE(input$cli)
    abx_sys <- input$tet != "none" && input$tet_wk > 0 && input$tet != "sub40"
    has_bpo <- isTRUE(input$bpo)
    has_ret <- input$ret != "none"

    if ((abx_top || abx_sys) && !has_bpo)
      msg <- c(msg, paste("⚠️ 항생제를 BPO 없이 사용하고 있습니다. 내성 분획이",
                          "상승하고 효과가 소실됩니다 (탭 5 참조). 가이드라인은",
                          "항생제 사용 시 BPO 또는 레티노이드 병용을 요구합니다."))
    if (abx_sys && input$tet_wk > 16)
      msg <- c(msg, "⚠️ 전신 항생제 16주 초과. 가이드라인 권고는 3–4개월 이내입니다.")
    if ((abx_top || abx_sys) && !has_ret)
      msg <- c(msg, paste("⚠️ 레티노이드 없이 항염 요법만 사용 중입니다.",
                          "미세면포 저장고가 남아 중단 후 재발합니다."))
    if (input$iso > 0 && input$iso * input$iso_wk * 7 < 120)
      msg <- c(msg, sprintf(paste("⚠️ 예상 누적 이소트레티노인 용량 %.0f mg/kg —",
                                  "목표 120–150 mg/kg 미만이며 재발 위험이 높습니다."),
                            input$iso * input$iso_wk * 7))
    if (input$iso > 0 && abx_sys)
      msg <- c(msg, "⛔ 테트라사이클린 + 이소트레티노인 병용은 특발성 두개내압상승 위험으로 금기입니다.")
    if (input$iso > 0)
      msg <- c(msg, "⛔ 가임 여성: 최기형성 — 두 가지 피임법과 임신검사 프로그램이 필수입니다.")
    if (!length(msg)) msg <- "✅ 특이 경고 없음."
    HTML(paste0("<ul>", paste0("<li>", msg, "</li>", collapse = ""), "</ul>"))
  })

  ## ---- 3. PK --------------------------------------------------------------
  output$p_pk_sys <- renderPlot({
    long_plot(sim(), c("CTETO", "CISOO", "COXOO", "CSPIO"),
              c("테트라사이클린", "이소트레티노인", "4-옥소 대사체", "칸레논"),
              "농도 (mg/L)", "전신 약물 농도")
  })
  output$p_pk_top <- renderPlot({
    long_plot(sim(), c("BPOS", "RETS", "CLIS", "AZES", "CLAS", "DAPS"),
              c("BPO", "레티노이드", "클린다마이신", "아젤라산",
                "클라스코테론", "답손"),
              "모낭 내 상대 농도", "국소 약물 모낭 저장고")
  })
  output$p_pk_eff <- renderPlot({
    long_plot(sim(), c("ISOEFFT", "AIEFFT", "CEEO"),
              c("이소트레티노인 효과 (0-1)", "테트라사이클린 항염 효과 (0-1)",
                "EE 농도 (pg/mL)"),
              "", "약력학적 효과 지표")
  })

  ## ---- 4. sebaceous / hormonal -------------------------------------------
  output$p_seb <- renderPlot({
    long_plot(sim(), c("SER", "SGM", "LIP", "DURAB"),
              c("피지 분비율 SER", "피지선 질량 SGM", "지질생성 구동 LIP",
                "영구 선 축소 DURAB"),
              "", "피지선 — 이소트레티노인의 표적", hline = 1)
  })
  output$p_horm <- renderPlot({
    long_plot(sim(), c("TT", "SHBG", "FAIO"),
              c("총 테스토스테론 (ng/dL)", "SHBG (nmol/L)", "유리 안드로겐 지수"),
              "", "안드로겐 축 — COC/스피로놀락톤/클라스코테론의 표적")
  })
  output$p_met <- renderPlot({
    long_plot(sim(), c("IGF1", "INS", "ARS"),
              c("IGF-1 (ng/mL)", "인슐린 (µU/mL)", "AR 신호 (상대)"),
              "", "대사·수용체 신호")
  })

  ## ---- 5. microbiology ----------------------------------------------------
  output$p_micro <- renderPlot({
    long_plot(sim(), c("CAP", "CAB", "CACNT", "FFA", "PORP"),
              c("부유성 C. acnes", "바이오필름", "총 부하",
                "유리지방산", "포르피린"),
              "상대값 (건강 기준 = 1)", "C. acnes 생태")
  })
  output$p_res <- renderPlot({
    d <- sim()
    ggplot(d, aes(week)) +
      geom_line(aes(y = RESF, colour = "내성 분획 RESF"), linewidth = 1.1) +
      geom_line(aes(y = CACNT / max(1e-6, max(d$CACNT)),
                    colour = "총 균 부하 (정규화)"), linewidth = 0.9) +
      geom_line(aes(y = INFLAM / max(1e-6, max(d$INFLAM)),
                    colour = "염증성 병변 (정규화)"), linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(1, 3, 2)]) +
      labs(x = "주 (weeks)", y = "0–1", title = "내성 선택압과 병변 반등") +
      theme_acn()
  })

  ## ---- 6. inflammation ----------------------------------------------------
  output$p_infl <- renderPlot({
    long_plot(sim(), c("TLR", "IL1B", "IL8", "TNF", "IL17", "NEU", "MMP"),
              c("TLR2 신호", "IL-1β", "IL-8", "TNF-α", "IL-17A",
                "호중구", "MMP"),
              "상대값", "선천 → 효과기 염증 캐스케이드", hline = 1)
  })
  output$p_ker <- renderPlot({
    long_plot(sim(), c("KER", "IL1A", "LA", "SQOX"),
              c("각질화 지수 KER", "IL-1α", "리놀레산 (희석)",
                "스쿠알렌 과산화물"),
              "상대값", "누두부 과각화 — 레티노이드의 표적", hline = 1)
  })

  ## ---- 7. lesions ---------------------------------------------------------
  output$p_les <- renderPlot({
    long_plot(sim(), c("MC", "CC", "OC", "PAP", "PUS", "NOD"),
              c("미세면포 (MC, 임상 미검출)", "폐쇄면포", "개방면포",
                "구진", "농포", "결절"),
              "병변 수", "병변 전이 사슬 — MC가 재발의 저장고")
  })
  output$p_iga <- renderPlot({
    long_plot(sim(), c("INFLAM", "NONINF", "IGA"),
              c("염증성 병변 수", "비염증성 병변 수", "IGA (0-4)"),
              "", "임상 엔드포인트")
  })
  output$tbl_end <- DT::renderDataTable({
    d <- sim(); b <- d[1, ]; e <- d[nrow(d), ]
    data.frame(
      지표 = c("염증성 병변", "비염증성 병변", "총 병변", "IGA",
               "피지 분비율", "C. acnes 총부하", "내성 분획",
               "PIH 지수", "흉터 지수"),
      기저 = round(c(b$INFLAM, b$NONINF, b$TOTLES, b$IGA, b$SERO,
                     b$CACNT, b$RESF, b$PIH, b$SCAR), 2),
      최종 = round(c(e$INFLAM, e$NONINF, e$TOTLES, e$IGA, e$SERO,
                     e$CACNT, e$RESF, e$PIH, e$SCAR), 2),
      변화율_pct = round(100 * (c(e$INFLAM, e$NONINF, e$TOTLES, e$IGA, e$SERO,
                                  e$CACNT, e$RESF, e$PIH, e$SCAR) -
                                c(b$INFLAM, b$NONINF, b$TOTLES, b$IGA, b$SERO,
                                  b$CACNT, b$RESF, b$PIH, b$SCAR)) /
                        pmax(1e-6, c(b$INFLAM, b$NONINF, b$TOTLES, b$IGA, b$SERO,
                                     b$CACNT, b$RESF, b$PIH, b$SCAR)), 1))
  }, options = list(dom = "t", pageLength = 12), rownames = FALSE)

  ## ---- 8. scenarios -------------------------------------------------------
  scn_sim <- eventReactive(input$go_scn, {
    ids <- as.integer(input$scn)
    if (!length(ids)) return(NULL)
    ACN_simulate(MOD, which = ids, delta = 24)
  })
  output$tbl_scn <- DT::renderDataTable({
    s <- scn_sim(); if (is.null(s)) return(NULL)
    ACN_summary(s, at_week = input$scn_wk)
  }, options = list(scrollX = TRUE, pageLength = 20), rownames = FALSE)
  output$p_scn <- renderPlot({
    s <- scn_sim(); if (is.null(s)) return(NULL)
    ggplot(s, aes(week, INFLAM, colour = factor(sid))) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = rep(PAL, 3),
                          labels = unique(paste0("S", s$sid, " ", s$scenario))) +
      labs(x = "주 (weeks)", y = "염증성 병변 수",
           title = "시나리오별 염증성 병변 경과") +
      theme_acn() + theme(legend.text = element_text(size = 8))
  })

  ## ---- 9. safety ----------------------------------------------------------
  output$p_saf <- renderPlot({
    long_plot(sim(), c("TG", "ALT", "KSER", "MUCO"),
              c("중성지방 (mg/dL)", "ALT (U/L)", "혈청 칼륨 (mEq/L)",
                "점막피부 독성 (0-10)"),
              "", "안전성 추적")
  })
  output$p_cum <- renderPlot({
    d <- sim()
    ggplot(d, aes(week, CUMISO)) +
      geom_line(linewidth = 1.1, colour = PAL[1]) +
      geom_hline(yintercept = c(120, 150), linetype = 2, colour = "grey35") +
      annotate("text", x = max(d$week) * 0.05, y = 135,
               label = "목표 누적 120–150 mg/kg", hjust = 0, size = 3.5) +
      labs(x = "주 (weeks)", y = "누적 이소트레티노인 (mg/kg)",
           title = "누적 용량 — 재발률을 결정하는 변수") +
      theme_acn()
  })
  output$txt_saf <- renderUI({
    d <- sim(); e <- d[nrow(d), ]
    HTML(sprintf(paste0(
      "<ul><li>최고 중성지방 <b>%.0f mg/dL</b> (기저 %.0f)</li>",
      "<li>최고 ALT <b>%.0f U/L</b></li>",
      "<li>최고 혈청 칼륨 <b>%.2f mEq/L</b></li>",
      "<li>최고 점막피부 독성 점수 <b>%.1f / 10</b> — 구순염은 사실상 전례에서 발생하며 ",
      "복약 확인의 대리지표로 쓰입니다.</li>",
      "<li>총 누적 이소트레티노인 <b>%.0f mg/kg</b></li></ul>"),
      max(d$TG), d$TG[1], max(d$ALT), max(d$KSER), max(d$MUCO), max(d$CUMISO)))
  })

  ## ---- 10. scarring -------------------------------------------------------
  output$p_scar <- renderPlot({
    long_plot(sim(), c("NOD", "MMP", "SCAR", "PIH"),
              c("결절 수", "MMP 활성", "위축성 흉터 지수 (비가역)",
                "PIH 지수"),
              "", "흉터와 색소침착의 축적")
  })
  output$tbl_scar <- DT::renderDataTable({
    d <- sim(); b <- d[1, ]; e <- d[nrow(d), ]
    data.frame(
      지표 = c("누적 결절-시간 (nodule-weeks)", "최종 흉터 지수",
               "최종 PIH 지수", "PIH 반감기 추정 (주)"),
      값 = round(c(sum(d$NOD) * (d$week[2] - d$week[1]),
                   e$SCAR, e$PIH, log(2) / (1.2e-3 * 168)), 2))
  }, options = list(dom = "t"), rownames = FALSE)

  ## ---- 11. help -----------------------------------------------------------
  output$txt_help <- renderPrint({
    cat("Acne vulgaris QSP model\n")
    cat("  - 55 ODE compartments, 213 parameters, 17 scenarios\n")
    cat("  - mechanistic map : acn_qsp_model.dot / .svg / .png\n")
    cat("  - model code      : acn_mrgsolve_model.R\n")
    cat("  - references      : acn_references.md\n\n")
    cat("Reading the model in one sentence:\n")
    cat("  every visible lesion descends from an invisible microcomedone (MC),\n")
    cat("  so anti-inflammatory therapy empties the downstream compartments fast\n")
    cat("  and rebounds, while retinoids drain the reservoir slowly and hold.\n\n")
    cat("Compartments:\n")
    print(names(mrgsolve::init(MOD)))
  })
}

shinyApp(ui, server)
