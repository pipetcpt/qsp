# =====================================================================
#  Bile Acid Diarrhoea / Malabsorption (BAM) — QSP Shiny dashboard
# ---------------------------------------------------------------------
#  The whole app is organised around ONE control panel idea: the disease
#  is written as two numbers, and everything else is a consequence.
#
#      phi    surviving ileal ASBT absorptive capacity  (slider)
#      kappa  ileal FGF19 sensor gain                   (slider)
#
#  Tab 2 ("The central experiment") lets the user FREEZE the FGF19 ->
#  CYP7A1 loop and watch the diarrhoea disappear while the malabsorption
#  stays exactly where it was.  That is the point of the whole model.
#
#  Run:  shiny::runApp("bam_shiny_app.R")
#  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
#  The model is defined in bam_mrgsolve_model.R and sourced here.
# =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("bam_mrgsolve_model.R")   # defines mod, ev_*, bam_readout, sim_bam

THEME <- theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93"),
        legend.position = "bottom")

PAL <- c(healthy = "#4f9d69", patient = "#b3272d", treated = "#1f6fb2",
         frozen  = "#8e6fc4", live = "#cc8400")

# ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel(HTML(paste0(
    "<b>담즙산 설사 / 담즙산 흡수장애 QSP 대시보드</b>",
    "<br><span style='font-size:14px;color:#555'>",
    "Bile Acid Diarrhoea — the loop measures what is ABSORBED, ",
    "the symptom is made by what is SPILLED</span>"))),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① 두 개의 병변축"),
      sliderInput("phi", HTML("<b>&phi;</b> — 잔존 회장 ASBT 흡수능"),
                  min = 0.02, max = 1.0, value = 0.20, step = 0.01),
      sliderInput("kappa", HTML("<b>&kappa;</b> — 회장 FGF19 센서 이득"),
                  min = 0.05, max = 1.2, value = 1.00, step = 0.05),
      helpText(HTML(paste0(
        "<small>Type 1 = &phi;&darr; · Type 2 = &kappa;&darr; · ",
        "Type 3 = 담낭절제 · Type 4 = 약물</small>"))),
      hr(),
      h4("② 구조적 수식자"),
      checkboxInput("chole", "담낭절제 (cholecystectomy)", FALSE),
      sliderInput("entmass", "회장 상피 기능 질량 (ENTMASS)",
                  min = 0.3, max = 1.0, value = 1.0, step = 0.05),
      sliderInput("bai", HTML("미생물 7&alpha;-탈수산화능 (BAI, relative)"),
                  min = 0.1, max = 2.5, value = 1.0, step = 0.1),
      sliderInput("aTR", "통과-흡수 양성되먹임 이득 (aTR)",
                  min = 0.4, max = 3.6, value = 1.55, step = 0.05),
      hr(),
      h4("③ 치료"),
      selectInput("seq", "담즙산 결합수지 (class A)",
                  c("없음" = "0", "colesevelam 1.875 g BID" = "1875_2",
                    "colesevelam 3.75 g BID" = "3750_2",
                    "colestyramine 4 g TID" = "4000_3")),
      selectInput("fxr", "FXR 작용제 (class B)",
                  c("없음" = "0", "obeticholic acid 25 mg QD" = "oca25",
                    "obeticholic acid 10 mg QD" = "oca10",
                    "tropifexor 90 µg QD" = "tro90")),
      selectInput("elo", "ASBT 억제제 (class C)",
                  c("없음" = "0", "elobixibat 5 mg QD" = "5",
                    "elobixibat 10 mg QD" = "10",
                    "elobixibat 15 mg QD" = "15")),
      selectInput("tra", "통과 조절제 (class D)",
                  c("없음" = "0", "loperamide 2 mg BID" = "lop2",
                    "loperamide 4 mg BID" = "lop4",
                    "ondansetron 4 mg BID" = "ond4")),
      checkboxInput("rif", "rifaximin 550 mg TID (class E)", FALSE),
      hr(),
      numericInput("days", "시뮬레이션 기간 (일)", 56, min = 21, max = 120),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary btn-block"),
      hr(),
      helpText(HTML(paste0(
        "<small>모든 파라미터의 출처와 보정 근거는 ",
        "<code>bam_mrgsolve_model.R</code> 하단 PROVENANCE 블록에, ",
        "검증 수치는 <code>bam_verification_output.txt</code>에 있습니다.",
        "</small>")))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        # ---------------------------------------------------------
        tabPanel(
          "1 · 환자 프로파일",
          br(),
          fluidRow(
            column(4, wellPanel(h4("진단 패널"), tableOutput("dxpanel"))),
            column(4, wellPanel(h4("증상"), tableOutput("sxpanel"))),
            column(4, wellPanel(h4("장기계 영향"), tableOutput("orgpanel")))
          ),
          h4("아형 분류 — 두 검사가 두 병변을 따로 읽는다"),
          plotOutput("quadrant", height = "420px"),
          helpText(HTML(paste0(
            "가로축 SeHCAT은 <b>&phi;</b>만, 세로축 C4는 <b>S = f(&phi;,&kappa;)</b>를 ",
            "읽습니다. 좌상 = 전형적 1형, 우상 = <b>SeHCAT 음성 담즙산 설사</b>(2형), ",
            "좌하 = 보상된 흡수장애(무증상), 우하 = 정상.")))
        ),

        # ---------------------------------------------------------
        tabPanel(
          "2 · 중심 실험 (되먹임 얼리기)",
          br(),
          h4("같은 흡수장애, 되먹임만 켜고 끄기"),
          helpText(HTML(paste0(
            "왼쪽은 FGF19&rarr;CYP7A1 되먹임이 살아 있는 생리적 환자, ",
            "오른쪽은 되먹임을 정상 작동점에 <b>동결</b>시킨 같은 환자입니다. ",
            "흡수장애(SeHCAT)는 두 경우 모두 똑같이 무너지지만, ",
            "설사는 되먹임이 살아 있는 쪽에서만 생깁니다."))),
          plotOutput("frozen", height = "460px"),
          br(),
          DT::dataTableOutput("frozentab")
        ),

        # ---------------------------------------------------------
        tabPanel(
          "3 · 담즙산 PK — 장간순환",
          br(),
          fluidRow(
            column(6, plotOutput("pkpool", height = "330px")),
            column(6, plotOutput("pkduo", height = "330px"))
          ),
          fluidRow(
            column(6, plotOutput("pkile", height = "330px")),
            column(6, plotOutput("pkserum", height = "330px"))
          )
        ),

        # ---------------------------------------------------------
        tabPanel(
          "4 · 센서 축 (FXR–FGF19–CYP7A1)",
          br(),
          fluidRow(
            column(6, plotOutput("fgf", height = "330px")),
            column(6, plotOutput("cyp", height = "330px"))
          ),
          fluidRow(
            column(6, plotOutput("c4plot", height = "330px")),
            column(6, plotOutput("loopgain", height = "330px"))
          ),
          helpText(HTML(paste0(
            "우하단은 닫힌 형태의 정상상태 해 ",
            "<code>S = Smax / (1 + (a&middot;&kappa;&middot;x&middot;S)^h)</code>, ",
            "x = f/(1-f). 되먹임이 로그적으로 뻣뻣하기 때문에 재순환비가 ",
            "65배 떨어져도 합성은 4배 남짓만 오릅니다.")))
        ),

        # ---------------------------------------------------------
        tabPanel(
          "5 · 대장 — 3 mM 역치",
          br(),
          fluidRow(
            column(7, plotOutput("colspecies", height = "360px")),
            column(5, plotOutput("threshold", height = "360px"))
          ),
          helpText(HTML(paste0(
            "분비 구동력 D<sub>s</sub>는 <b>자유</b> 담즙산의 효력 가중합입니다 ",
            "(CA 0.15 / CDCA 1.00 / DCA 1.15 / LCA 0.30). ",
            "역치 3 mM은 Mekjian 1971의 사람 대장 관류 실험에서 직접 왔습니다 ",
            "(PMID 4938344)."))),
          plotOutput("transit", height = "320px")
        ),

        # ---------------------------------------------------------
        tabPanel(
          "6 · 임상 엔드포인트",
          br(),
          fluidRow(
            column(6, plotOutput("stoolplot", height = "330px")),
            column(6, plotOutput("bmplot", height = "330px"))
          ),
          fluidRow(
            column(6, plotOutput("fatplot", height = "330px")),
            column(6, plotOutput("ldlplot", height = "330px"))
          )
        ),

        # ---------------------------------------------------------
        tabPanel(
          "7 · 치료 계열 비교",
          br(),
          h4("같은 환자, 다섯 계열의 진입점"),
          plotOutput("drugclass", height = "440px"),
          br(),
          DT::dataTableOutput("drugtab"),
          helpText(HTML(paste0(
            "<b>결합수지는 자기가 치료하려는 구동력을 스스로 올립니다</b> ",
            "— 관강 담즙산을 묶으면 센서 리간드도 사라져 C4가 오릅니다. ",
            "FXR 작용제는 그 신호를 되돌려주므로 병용이 상승적입니다.")))
        ),

        # ---------------------------------------------------------
        tabPanel(
          "8 · 100 cm 규칙",
          br(),
          h4("절제 길이 → 표현형 전환은 유도된 결과입니다"),
          helpText(HTML(paste0(
            "입력은 세 가지뿐입니다: ASBT 밀도의 원위 편중(감쇠거리 60 cm), ",
            "CYP7A1 보상 천장(6배), 유효 임계미셀농도(1.5 mM). ",
            "절제 결과 데이터는 하나도 사용하지 않았습니다."))),
          plotOutput("rule100", height = "480px"),
          br(),
          DT::dataTableOutput("ruletab")
        ),

        # ---------------------------------------------------------
        tabPanel(
          "9 · 바이오마커 & 미생물총",
          br(),
          fluidRow(
            column(6, plotOutput("baiplot", height = "340px")),
            column(6, plotOutput("biomk", height = "340px"))
          ),
          helpText(HTML(paste0(
            "미생물총은 대장에 도달하는 담즙산의 <b>양</b>을 거의 바꾸지 않고 ",
            "<b>효력</b>을 바꿉니다. 그래서 대변 담즙산이 같아도 증상이 ",
            "3배까지 달라질 수 있습니다 — 검사와 증상이 어긋나는 두 번째 이유."))),
          plotOutput("sehcatcurve", height = "330px")
        ),

        # ---------------------------------------------------------
        tabPanel(
          "10 · 시나리오 라이브러리",
          br(),
          h4("14개 사전 정의 시나리오"),
          actionButton("runlib", "라이브러리 전체 실행 (수 분 소요)",
                       class = "btn-warning"),
          br(), br(),
          DT::dataTableOutput("libtab"),
          br(),
          helpText(HTML(paste0(
            "이 표의 수치는 <code>bam_verify_python.py</code>가 독립적으로 ",
            "재현합니다 (<code>bam_verification_output.txt</code>).")))
        )
      )
    )
  )
)

# ---------------------------------------------------------------------
server <- function(input, output, session) {

  build_ev <- function() {
    e <- NULL
    add <- function(a, b) if (is.null(a)) b else a + b
    d <- input$days
    if (input$seq != "0") {
      pr <- strsplit(input$seq, "_")[[1]]
      e <- add(e, ev_seq(as.numeric(pr[1]), as.numeric(pr[2]), d))
    }
    if (input$fxr == "oca25") e <- add(e, ev_oca(25, d))
    if (input$fxr == "oca10") e <- add(e, ev_oca(10, d))
    if (input$fxr == "tro90") e <- add(e, ev_tro(0.090, d))
    if (input$elo != "0")     e <- add(e, ev_elo(as.numeric(input$elo), d))
    if (input$tra == "lop2")  e <- add(e, ev_lop(2, d))
    if (input$tra == "lop4")  e <- add(e, ev_lop(4, d))
    if (input$tra == "ond4")  e <- add(e, ev_ond(4, d))
    if (isTRUE(input$rif))    e <- add(e, ev_rif(550, d))
    e
  }

  pars <- reactive({
    p <- list(phi = input$phi, kappa = input$kappa,
              ENTMASS = input$entmass, aTR = input$aTR,
              kbai = 0.16 * input$bai)
    if (isTRUE(input$chole)) p$GBfun <- 0
    p
  })

  # ---- main simulations, recomputed only on button press -----------
  simset <- eventReactive(input$go, {
    withProgress(message = "시뮬레이션 중...", value = 0, {
      p <- pars(); e <- build_ev()
      incProgress(0.15, detail = "healthy control")
      h <- do.call(sim_bam, c(list(days = input$days, events = NULL),
                              list(phi = 1, kappa = 1)))
      incProgress(0.3, detail = "patient, untreated")
      u <- do.call(sim_bam, c(list(days = input$days, events = NULL), p))
      incProgress(0.3, detail = "patient, treated")
      tr <- if (is.null(e)) u else
        do.call(sim_bam, c(list(days = input$days, events = e), p))
      incProgress(0.25, detail = "feedback frozen")
      fz <- do.call(sim_bam, c(list(days = input$days, events = NULL),
                               c(p, list(FEEDBACK = 0))))
      list(healthy = h, untreated = u, treated = tr, frozen = fz)
    })
  }, ignoreNULL = FALSE)

  ro <- reactive({
    s <- simset()
    list(healthy = bam_readout(s$healthy), untreated = bam_readout(s$untreated),
         treated = bam_readout(s$treated), frozen = bam_readout(s$frozen))
  })

  last_day <- function(out) {
    d <- as.data.frame(out); te <- max(d$time)
    dplyr::filter(d, time >= te - 48) %>% mutate(tod = time - (te - 48))
  }

  # ================= TAB 1 =========================================
  output$dxpanel <- renderTable({
    r <- ro()
    data.frame(
      항목 = c("SeHCAT 7일 잔류 (%)", "공복 혈청 C4 (ng/mL)",
               "공복 FGF19 (pg/mL)", "대변 담즙산 (µmol/day)",
               "담즙산 풀 (g)", "회장 1회 보존율 f"),
      환자 = c(sprintf("%.1f", r$treated$SeHCAT_pct),
               sprintf("%.1f", r$treated$C4),
               sprintf("%.0f", r$treated$FGF19),
               sprintf("%.0f", r$treated$fecBA_umold),
               sprintf("%.2f", r$treated$pool_g),
               sprintf("%.4f", 1 - r$treated$kturn_day / 4.5)),
      정상 = c(sprintf("%.1f", r$healthy$SeHCAT_pct),
               sprintf("%.1f", r$healthy$C4),
               sprintf("%.0f", r$healthy$FGF19),
               sprintf("%.0f", r$healthy$fecBA_umold),
               sprintf("%.2f", r$healthy$pool_g), "~0.97"))
  }, striped = TRUE, width = "100%")

  output$sxpanel <- renderTable({
    r <- ro()
    data.frame(
      항목 = c("대변 수분 (mL/day)", "배변 횟수 (회/day)",
               "Bristol 점수", "대장 분비구동력 (mM-eq)",
               "대장 통과 (상대)"),
      환자 = c(sprintf("%.0f", r$treated$stool_mL),
               sprintf("%.2f", r$treated$BM_per_day),
               sprintf("%.1f", r$treated$bristol),
               sprintf("%.2f", r$treated$Ds),
               sprintf("%.2f", r$treated$TRANS)),
      정상 = c(sprintf("%.0f", r$healthy$stool_mL),
               sprintf("%.2f", r$healthy$BM_per_day),
               sprintf("%.1f", r$healthy$bristol),
               sprintf("%.2f", r$healthy$Ds), "1.00"))
  }, striped = TRUE, width = "100%")

  output$orgpanel <- renderTable({
    r <- ro()
    data.frame(
      항목 = c("대변 지방 (g/day)", "십이지장 최고 [BA] (mM)",
               "LDL-C (mg/dL)", "요중 수산 (mg/day)"),
      환자 = c(sprintf("%.1f", r$treated$fecfat_g),
               sprintf("%.1f", r$treated$Cduo_peak),
               sprintf("%.0f", r$treated$LDLC),
               sprintf("%.0f", r$treated$UOX)),
      정상 = c(sprintf("%.1f", r$healthy$fecfat_g),
               sprintf("%.1f", r$healthy$Cduo_peak),
               sprintf("%.0f", r$healthy$LDLC),
               sprintf("%.0f", r$healthy$UOX)))
  }, striped = TRUE, width = "100%")

  output$quadrant <- renderPlot({
    r <- ro()
    pt <- data.frame(x = r$treated$SeHCAT_pct, y = r$treated$C4,
                     lab = "현재 환자")
    hc <- data.frame(x = r$healthy$SeHCAT_pct, y = r$healthy$C4,
                     lab = "정상 대조")
    ggplot() +
      annotate("rect", xmin = 0, xmax = 15, ymin = 48, ymax = Inf,
               fill = "#ffd9b3", alpha = .5) +
      annotate("rect", xmin = 15, xmax = Inf, ymin = 48, ymax = Inf,
               fill = "#d9c9f0", alpha = .5) +
      annotate("rect", xmin = 0, xmax = 15, ymin = 0, ymax = 48,
               fill = "#e0f0d8", alpha = .5) +
      annotate("text", x = 7, y = 90, label = "1형\n(회장 병변)", size = 5) +
      annotate("text", x = 40, y = 90,
               label = "2형 — SeHCAT 음성\n담즙산 설사", size = 5) +
      annotate("text", x = 7, y = 20, label = "보상된 흡수장애\n(증상 적음)", size = 4.5) +
      annotate("text", x = 40, y = 20, label = "정상", size = 5) +
      geom_vline(xintercept = 15, linetype = 2) +
      geom_hline(yintercept = 48, linetype = 2) +
      geom_point(data = hc, aes(x, y), size = 5, colour = PAL["healthy"]) +
      geom_point(data = pt, aes(x, y), size = 6, colour = PAL["patient"]) +
      geom_text(data = pt, aes(x, y, label = lab), vjust = -1.4, size = 4.5) +
      scale_x_continuous("SeHCAT 7일 잔류 (%)  →  φ를 읽는다",
                         limits = c(0, 60)) +
      scale_y_continuous("공복 혈청 C4 (ng/mL)  →  S = f(φ, κ)를 읽는다",
                         limits = c(0, 120)) + THEME
  })

  # ================= TAB 2 =========================================
  output$frozen <- renderPlot({
    r <- ro()
    d <- bind_rows(
      data.frame(cond = "되먹임 LIVE", metric = "SeHCAT (%)",
                 v = r$untreated$SeHCAT_pct),
      data.frame(cond = "되먹임 FROZEN", metric = "SeHCAT (%)",
                 v = r$frozen$SeHCAT_pct),
      data.frame(cond = "되먹임 LIVE", metric = "대장 담즙산 부하 (µmol/d)",
                 v = r$untreated$fecBA_umold),
      data.frame(cond = "되먹임 FROZEN", metric = "대장 담즙산 부하 (µmol/d)",
                 v = r$frozen$fecBA_umold),
      data.frame(cond = "되먹임 LIVE", metric = "혈청 C4 (ng/mL)",
                 v = r$untreated$C4),
      data.frame(cond = "되먹임 FROZEN", metric = "혈청 C4 (ng/mL)",
                 v = r$frozen$C4),
      data.frame(cond = "되먹임 LIVE", metric = "대변 수분 (mL/d)",
                 v = r$untreated$stool_mL),
      data.frame(cond = "되먹임 FROZEN", metric = "대변 수분 (mL/d)",
                 v = r$frozen$stool_mL))
    ggplot(d, aes(cond, v, fill = cond)) +
      geom_col(width = .6) +
      geom_text(aes(label = sprintf("%.1f", v)), vjust = -.4, size = 4.2) +
      facet_wrap(~metric, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = c("되먹임 LIVE" = unname(PAL["live"]),
                                   "되먹임 FROZEN" = unname(PAL["frozen"]))) +
      labs(x = NULL, y = NULL, fill = NULL) +
      expand_limits(y = 0) + THEME
  })

  output$frozentab <- DT::renderDataTable({
    r <- ro()
    tab <- bind_rows(
      cbind(조건 = "되먹임 LIVE", r$untreated),
      cbind(조건 = "되먹임 FROZEN", r$frozen),
      cbind(조건 = "정상 대조", r$healthy))
    DT::datatable(tab, options = list(dom = "t", scrollX = TRUE),
                  rownames = FALSE) %>%
      DT::formatRound(2:ncol(tab), 2)
  })

  # ================= TAB 3 =========================================
  plot_ts <- function(var, ylab, title) {
    renderPlot({
      s <- simset()
      d <- bind_rows(
        cbind(grp = "정상", last_day(s$healthy)),
        cbind(grp = "환자(무치료)", last_day(s$untreated)),
        cbind(grp = "환자(치료)", last_day(s$treated)))
      ggplot(d, aes(tod, .data[[var]], colour = grp)) +
        geom_line(linewidth = .9) +
        scale_colour_manual(values = c("정상" = unname(PAL["healthy"]),
                                       "환자(무치료)" = unname(PAL["patient"]),
                                       "환자(치료)" = unname(PAL["treated"]))) +
        scale_x_continuous("마지막 48시간 (h)", breaks = seq(0, 48, 8)) +
        labs(y = ylab, title = title, colour = NULL) + THEME
    })
  }
  output$pkpool   <- plot_ts("POOLtot", "µmol", "순환 담즙산 풀")
  output$pkduo    <- plot_ts("CDUO", "mM", "십이지장 담즙산 농도 (CMC 1.5 mM)")
  output$pkile    <- plot_ts("ILE1", "µmol", "말단 회장 관강 (ASBT 구획 1)")
  output$pkserum  <- plot_ts("SERUMBA", "µM", "전신 혈청 담즙산")
  output$fgf      <- plot_ts("FGF19", "pg/mL", "혈장 FGF19 (식후 90-120분 정점)")
  output$cyp      <- plot_ts("CYP7A1", "상대값", "간 CYP7A1 (천장 = 6배)")
  output$c4plot   <- plot_ts("C4", "ng/mL", "혈청 C4 (컷오프 48)")
  output$stoolplot<- plot_ts("WCOL", "mL", "대장 관강 수분")
  output$fatplot  <- plot_ts("FATEFF", "분율", "지방 흡수 효율")
  output$ldlplot  <- plot_ts("LDLC", "mg/dL", "LDL-C")
  output$transit  <- plot_ts("TRANS", "상대값", "대장 통과 속도 (양성 되먹임)")

  output$bmplot <- renderPlot({
    r <- ro()
    d <- data.frame(
      grp = factor(c("정상", "환자(무치료)", "환자(치료)"),
                   levels = c("정상", "환자(무치료)", "환자(치료)")),
      bm  = c(r$healthy$BM_per_day, r$untreated$BM_per_day, r$treated$BM_per_day),
      br  = c(r$healthy$bristol, r$untreated$bristol, r$treated$bristol))
    ggplot(d, aes(grp, bm, fill = grp)) + geom_col(width = .6) +
      geom_text(aes(label = sprintf("%.2f 회/일\nBristol %.1f", bm, br)),
                vjust = -.3, size = 4) +
      scale_fill_manual(values = unname(PAL[c("healthy", "patient", "treated")])) +
      labs(x = NULL, y = "배변 횟수 (회/일)", title = "임상 엔드포인트") +
      expand_limits(y = 0) + guides(fill = "none") + THEME
  })

  # ================= TAB 4: analytic loop ==========================
  output$loopgain <- renderPlot({
    S0 <- 700; Smax <- 6 * S0; h <- 1.5; f0 <- 0.975
    x0 <- f0 / (1 - f0)
    a <- uniroot(function(a) S0 * (1 + (a * x0 * S0)^h) - Smax,
                 c(1e-14, 1))$root
    solveS <- function(f, k)
      uniroot(function(S) S * (1 + (a * k * f / (1 - f) * S)^h) - Smax,
              c(1e-6, Smax))$root
    ff <- seq(0.55, 0.99, length.out = 120)
    d <- bind_rows(lapply(c(1, 0.45, 0.25), function(k)
      data.frame(f = ff, kappa = factor(k),
                 S = sapply(ff, solveS, k = k))))
    ggplot(d, aes(f, S, colour = kappa)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = Smax, linetype = 2) +
      annotate("text", x = 0.6, y = Smax * .96,
               label = "보상 천장 Smax = 6 × 기저", size = 4) +
      scale_x_continuous("회장 1회 보존율 f") +
      scale_y_continuous("정상상태 합성 = 대장 부하 (µmol/day)") +
      labs(colour = "κ", title = "닫힌 형태 정상상태 해") + THEME
  })

  # ================= TAB 5 =========================================
  output$colspecies <- renderPlot({
    s <- simset()
    d <- last_day(s$treated) %>%
      select(tod, CCA, CCDCA, CDCA, CLCA) %>%
      pivot_longer(-tod, names_to = "종", values_to = "µmol") %>%
      mutate(종 = recode(종, CCA = "cholate (w 0.15)",
                         CCDCA = "chenodeoxycholate (w 1.00)",
                         CDCA = "deoxycholate (w 1.15)",
                         CLCA = "lithocholate (w 0.30)"))
    ggplot(d, aes(tod, µmol, fill = 종)) + geom_area(alpha = .85) +
      scale_x_continuous("마지막 48시간 (h)", breaks = seq(0, 48, 8)) +
      labs(title = "대장 담즙산 조성 — 효력은 종에 따라 8배 차이", fill = NULL) +
      THEME
  })

  output$threshold <- renderPlot({
    s <- simset()
    d <- bind_rows(
      cbind(grp = "정상", last_day(s$healthy)),
      cbind(grp = "환자(치료)", last_day(s$treated)))
    ggplot(d, aes(tod, DS, colour = grp)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = 3, linetype = 2, colour = "#b3272d") +
      annotate("text", x = 24, y = 3.25,
               label = "분비 역치 3 mM (Mekjian 1971)", size = 4,
               colour = "#b3272d") +
      scale_colour_manual(values = unname(PAL[c("healthy", "treated")])) +
      scale_x_continuous("마지막 48시간 (h)", breaks = seq(0, 48, 8)) +
      labs(y = "Ds (mM CDCA 등가)", colour = NULL,
           title = "대장 분비 구동력") + THEME
  })

  # ================= TAB 7: drug classes ===========================
  drugcmp <- eventReactive(input$go, {
    withProgress(message = "치료 계열 비교 중...", value = 0, {
      base <- list(phi = input$phi, kappa = input$kappa)
      sc <- list(
        "무치료"                     = NULL,
        "A 결합수지 (colesevelam)"   = ev_seq(1875, 2, 56),
        "B FXR 작용제 (OCA 25 mg)"   = ev_oca(25, 56),
        "B FXR 작용제 (tropifexor)"  = ev_tro(0.090, 56),
        "C ASBT 억제제 (elobixibat)" = ev_elo(10, 56),
        "D 통과 (loperamide)"        = ev_lop(2, 56),
        "D 통과 (ondansetron)"       = ev_ond(4, 56),
        "E 미생물 (rifaximin)"       = ev_rif(550, 56),
        "A + B 병용"                 = ev_seq(1875, 2, 56) + ev_oca(25, 56))
      out <- lapply(names(sc), function(nm) {
        incProgress(1 / length(sc), detail = nm)
        o <- do.call(sim_bam, c(list(days = 56, events = sc[[nm]]), base))
        cbind(치료 = nm, bam_readout(o))
      })
      bind_rows(out)
    })
  }, ignoreNULL = FALSE)

  output$drugclass <- renderPlot({
    d <- drugcmp()
    base <- d$stool_mL[d$치료 == "무치료"]
    d <- d %>% mutate(변화 = 100 * (stool_mL - base) / base,
                      C4변화 = 100 * (C4 - C4[치료 == "무치료"]) /
                        C4[치료 == "무치료"])
    dd <- d %>% select(치료, `대변 수분 변화 (%)` = 변화,
                       `혈청 C4 변화 (%)` = C4변화) %>%
      pivot_longer(-치료)
    ggplot(dd, aes(reorder(치료, value), value, fill = value > 0)) +
      geom_col() + coord_flip() + facet_wrap(~name, scales = "free_x") +
      geom_hline(yintercept = 0) +
      scale_fill_manual(values = c("TRUE" = "#b3272d", "FALSE" = "#1f6fb2")) +
      labs(x = NULL, y = "기저 대비 변화 (%)") +
      guides(fill = "none") + THEME
  })

  output$drugtab <- DT::renderDataTable({
    d <- drugcmp()
    DT::datatable(d, options = list(dom = "t", scrollX = TRUE, pageLength = 12),
                  rownames = FALSE) %>% DT::formatRound(2:ncol(d), 2)
  })

  # ================= TAB 8: the 100 cm rule ========================
  rule <- eventReactive(input$go, {
    withProgress(message = "절제 길이 스윕...", value = 0, {
      Ls <- c(0, 20, 40, 60, 80, 100, 120, 150, 180)
      bind_rows(lapply(Ls, function(L) {
        incProgress(1 / length(Ls), detail = paste0(L, " cm"))
        o <- sim_bam(days = 56, events = NULL, phi = phi_from_resection(L))
        cbind(절제_cm = L, phi = round(phi_from_resection(L), 3), bam_readout(o))
      }))
    })
  }, ignoreNULL = FALSE)

  output$rule100 <- renderPlot({
    d <- rule()
    dd <- d %>% select(절제_cm, `대변 지방 (g/day)` = fecfat_g,
                       `대변 수분 (mL/day)` = stool_mL,
                       `담즙산 풀 (g)` = pool_g,
                       `십이지장 최고 [BA] (mM)` = Cduo_peak) %>%
      pivot_longer(-절제_cm)
    hl <- data.frame(name = c("대변 지방 (g/day)", "십이지장 최고 [BA] (mM)"),
                     y = c(7, 1.5))
    ggplot(dd, aes(절제_cm, value)) +
      geom_line(linewidth = 1, colour = PAL["patient"]) +
      geom_point(size = 2.2) +
      geom_hline(data = hl, aes(yintercept = y), linetype = 2,
                 colour = "#b3272d") +
      geom_vline(xintercept = 100, linetype = 3, linewidth = .9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "절제된 말단 회장 길이 (cm)", y = NULL,
           title = "점선 = 임상 교과서의 100 cm 기준") + THEME
  })

  output$ruletab <- DT::renderDataTable({
    d <- rule()
    DT::datatable(d, options = list(dom = "t", scrollX = TRUE, pageLength = 10),
                  rownames = FALSE) %>% DT::formatRound(3:ncol(d), 2)
  })

  # ================= TAB 9 =========================================
  baisweep <- eventReactive(input$go, {
    withProgress(message = "미생물총 스윕...", value = 0, {
      bs <- c(0.2, 0.5, 1.0, 1.5, 2.0)
      bind_rows(lapply(bs, function(b) {
        incProgress(1 / length(bs), detail = paste0("BAI ", b))
        o <- sim_bam(days = 40, events = NULL, phi = input$phi,
                     kappa = input$kappa, kbai = 0.16 * b)
        cbind(BAI = b, bam_readout(o))
      }))
    })
  }, ignoreNULL = FALSE)

  output$baiplot <- renderPlot({
    d <- baisweep() %>%
      mutate(`대변 담즙산 (상대)` = fecBA_umold / fecBA_umold[BAI == 1],
             `대변 수분 (상대)` = stool_mL / stool_mL[BAI == 1]) %>%
      select(BAI, `대변 담즙산 (상대)`, `대변 수분 (상대)`) %>%
      pivot_longer(-BAI)
    ggplot(d, aes(BAI, value, colour = name)) +
      geom_line(linewidth = 1) + geom_point(size = 2.5) +
      geom_hline(yintercept = 1, linetype = 2) +
      labs(x = "7α-탈수산화능 (BAI, 상대)", y = "BAI = 1 대비 배수",
           colour = NULL, title = "부하는 그대로, 효력만 바뀐다") + THEME
  })

  output$biomk <- renderPlot({
    d <- baisweep()
    ggplot(d, aes(Ds, stool_mL)) +
      geom_path(linewidth = 1, colour = PAL["treated"]) +
      geom_point(aes(size = BAI), colour = PAL["patient"]) +
      geom_vline(xintercept = 3, linetype = 2, colour = "#b3272d") +
      labs(x = "대장 분비 구동력 Ds (mM-eq)", y = "대변 수분 (mL/day)",
           title = "역치 통과가 증상을 만든다") + THEME
  })

  output$sehcatcurve <- renderPlot({
    d <- data.frame(k = seq(0.02, 1.2, length.out = 300)) %>%
      mutate(R7 = 100 * exp(-7 * k))
    ggplot(d, aes(k, R7)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = c(5, 10, 15), linetype = 2,
                 colour = c("#b3272d", "#cc8400", "#4f9d69")) +
      annotate("text", x = 1.0, y = c(8, 13, 18),
               label = c("5 % 중증", "10 % 중등도", "15 % 정상 하한"),
               size = 4) +
      scale_x_continuous("풀 회전율 k = 대변 손실 / 풀 (/day)") +
      scale_y_continuous("SeHCAT 7일 잔류 (%)") +
      labs(title = "SeHCAT은 압축된 척도다 — 임상 등급 전체가 좁은 k 구간에 몰려 있다") +
      THEME
  })

  # ================= TAB 10 ========================================
  lib <- eventReactive(input$runlib, {
    withProgress(message = "시나리오 라이브러리 실행 중...", value = 0.1, {
      bam_scenarios()
    })
  })
  output$libtab <- DT::renderDataTable({
    d <- lib()
    DT::datatable(d, options = list(scrollX = TRUE, pageLength = 15),
                  rownames = FALSE) %>% DT::formatRound(2:ncol(d), 2)
  })
}

shinyApp(ui, server)
