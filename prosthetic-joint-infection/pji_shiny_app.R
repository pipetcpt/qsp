## =============================================================================
##  pji_shiny_app.R — 인공관절 감염 (PJI) QSP 대시보드
##  Prosthetic Joint Infection — interactive QSP dashboard
##
##  실행 (run):
##      # 같은 디렉토리에서
##      shiny::runApp("pji_shiny_app.R")
##  요구 패키지: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
##
##  이 앱은 pji_mrgsolve_model.R 의 모델을 그대로 불러와 13개 탭으로 보여줍니다.
##  탭 3(골농도/MBEC 비)과 탭 7(변이 공급)이 이 모델의 논지입니다 — 나머지 탭은
##  그 두 장의 그림이 왜 임상 지침의 모양을 결정하는지를 보여주는 뒷받침입니다.
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- 모델 로드 ---------------------------------------------------------------
## pji_mrgsolve_model.R 은 모델 객체 pji_mod 와 시나리오 목록 scen 을 만든다.
if (!exists("pji_mod")) {
  src <- file.path(dirname(sys.frame(1)$ofile %||% "."), "pji_mrgsolve_model.R")
  if (!file.exists(src)) src <- "pji_mrgsolve_model.R"
  source(src, local = FALSE)
}

DAY <- 24
THEME <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        plot.title       = element_text(face = "bold", size = 14),
        strip.text       = element_text(face = "bold"))

DRUGCOL <- c(VAN = "#2E6DA4", RIF = "#C0392B", LVX = "#27AE60",
             DAP = "#8E44AD", LZD = "#E67E22", CFZ = "#16A085",
             LOCAL = "#7F8C8D")
POPCOL  <- c("유주균 (planktonic)"       = "#3498DB",
             "바이오필름 (biofilm)"      = "#E67E22",
             "지속균 (persister)"        = "#C0392B",
             "SCV"                       = "#8E44AD",
             "세포내 (intracellular)"    = "#16A085",
             "rpoB 내성"                 = "#000000",
             "gyrA 내성"                 = "#7F8C8D")

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("인공관절 감염 (PJI) QSP 대시보드 — Prosthetic Joint Infection"),
  tags$p(style = "color:#555;margin-top:-10px;",
         "Implant-associated Staphylococcus aureus osteomyelitis · 53-ODE QSP model · 교육/연구용"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("① 환자·병소"),
      radioButtons("organism", "균주", c("MRSA" = 0, "MSSA" = 1), selected = 0, inline = TRUE),
      sliderInput("imsup", "숙주 면역 계수 (1 = 정상)", 0.2, 1.2, 1.0, step = 0.1),
      numericInput("inoc", "접종 균량 (CFU)", 100, min = 1, max = 1e8),
      checkboxInput("implant", "임플란트 존재 (foreign body)", TRUE),
      sliderInput("tdx", "진단·수술 시각 (일)", 7, 240, 21, step = 7),

      hr(), h4("② 수술 전략"),
      selectInput("surg", "1차 수술",
                  c("없음 (항생제 단독)"            = "none",
                    "DAIR (변연절제·임플란트 유지)" = "dair",
                    "1단계 교환"                    = "one",
                    "2단계 교환 (+스페이서)"        = "two")),
      sliderInput("logk", "수술 균 제거량 (log10)", 0, 7, 2.5, step = 0.5),
      conditionalPanel("input.surg == 'two'",
        sliderInput("spacer", "스페이서 용출 가능 반코마이신 (mg)", 0, 3000, 1900, step = 100),
        sliderInput("reimp", "재치환 시점 (수술 후 주)", 4, 16, 8, step = 1)),

      hr(), h4("③ 항생제"),
      checkboxGroupInput("drugs", "요법",
        c("반코마이신 IV 1 g q12h"      = "VAN",
          "리팜피신 PO 450 mg q12h"     = "RIF",
          "레보플록사신 PO 750 mg qd"   = "LVX",
          "답토마이신 IV 8 mg/kg qd"    = "DAP",
          "리네졸리드 PO 600 mg q12h"   = "LZD",
          "세파졸린 IV 2 g q8h (MSSA)"  = "CFZ"),
        selected = c("LVX", "RIF")),
      sliderInput("wk_backbone", "주 약제 투여기간 (주)", 2, 26, 6, step = 1),
      sliderInput("wk_rif",      "리팜피신 투여기간 (주)", 0, 26, 12, step = 1),

      hr(), h4("④ 시뮬레이션"),
      sliderInput("tend", "관찰 기간 (일)", 90, 730, 365, step = 30),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary btn-block"),
      br(),
      tags$small(style = "color:#777;",
        "요약: 전신 항생제로 도달 가능한 유리 골농도는 리팜피신을 제외하면 ",
        "MBEC의 1~2%에 불과합니다. 수술은 그 사실을 우회하는 것이 아니라, ",
        "잔존 균량 N을 낮춰 rpoB 변이체가 존재할 확률 1−e^(−μN)을 떨어뜨려 ",
        "리팜피신을 쓸 수 있게 만드는 단계입니다.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        tabPanel("1. 개요",
          br(), fluidRow(
            column(3, wellPanel(h5("최종 균량 (log10 CFU)"), h3(textOutput("kpi_burden")))),
            column(3, wellPanel(h5("치료 성공 확률"),        h3(textOutput("kpi_cure")))),
            column(3, wellPanel(h5("P(사전 rpoB 변이체)"),   h3(textOutput("kpi_prpob")))),
            column(3, wellPanel(h5("2018 ICM 점수"),         h3(textOutput("kpi_icm"))))),
          plotOutput("p_overview", height = "460px"),
          tags$hr(), htmlOutput("narrative")),

        tabPanel("2. 항생제 PK (혈장 vs 골)",
          br(), plotOutput("p_pk", height = "560px"),
          tags$small("점선 = 혈장 총농도, 실선 = 유리 골농도. ",
                     "두 선의 간격이 곧 (골 침투율 x 유리 분율)입니다.")),

        tabPanel("3. ★ 골농도 / MBEC 비",
          br(), h4("이 탭이 이 모델의 논지입니다"),
          plotOutput("p_ratio", height = "380px"),
          br(), DT::dataTableOutput("t_ratio"),
          tags$small("비율 1을 넘지 못하면 그 약제는 바이오필름을 제거할 수 없습니다. ",
                     "전신 투여로 1에 근접하는 것은 리팜피신뿐이고, 1을 넘기는 것은 ",
                     "항생제 함유 시멘트의 국소 농도뿐입니다.")),

        tabPanel("4. 세균 아집단",
          br(), plotOutput("p_pops", height = "520px"),
          checkboxInput("logscale", "log10 축", TRUE)),

        tabPanel("5. 바이오필름·성숙",
          br(), plotOutput("p_biofilm", height = "520px"),
          tags$small("EPS 성숙도는 바이오필름 내약 배수(MBEC/MIC)를 실제로 곱해 들어가는 ",
                     "상태변수입니다. DAIR가 '3주 안에' 성립하는 이유가 여기 있습니다.")),

        tabPanel("6. 국소 용출 (스페이서)",
          br(), plotOutput("p_local", height = "480px"),
          tags$small("항생제 함유 PMMA는 초기 수일간만 MBEC를 넘습니다. ",
                     "그 뒤로는 전신 요법과 같은 영역으로 떨어집니다.")),

        tabPanel("7. ★ 내성 출현 / 변이 공급",
          br(), plotOutput("p_res", height = "380px"),
          br(), plotOutput("p_mutsupply", height = "300px"),
          tags$small("아래 그림은 시뮬레이션이 아니라 산수입니다: ",
                     "P = 1 − exp(−μN), μ = 1e−8. 변연절제 전후의 N을 표시했습니다.")),

        tabPanel("8. 숙주 면역·염증",
          br(), plotOutput("p_immune", height = "560px")),

        tabPanel("9. 골용해·임플란트 이완",
          br(), plotOutput("p_bone", height = "560px")),

        tabPanel("10. 바이오마커·엔드포인트",
          br(), plotOutput("p_biomarker", height = "520px"),
          tags$small("점선 = 2018 ICM 진단 문턱 (CRP 10 mg/L, 관절액 WBC 3,000/uL, ",
                     "다형핵구 80%, alpha-defensin 1.0 S/CO)")),

        tabPanel("11. 시나리오 비교",
          br(), actionButton("runall", "21개 표준 시나리오 실행", class = "btn-success"),
          br(), br(), DT::dataTableOutput("t_scen"),
          br(), plotOutput("p_scen", height = "460px")),

        tabPanel("12. 독성·안전성",
          br(), plotOutput("p_tox", height = "520px"),
          tags$small("반코마이신 AUC24 > 600 mg·h/L에서 신독성, ",
                     "리네졸리드 장기 투여에서 혈소판 감소, 리팜피신에서 ALT 상승.")),

        tabPanel("13. 파라미터·문헌",
          br(), DT::dataTableOutput("t_par"),
          br(), htmlOutput("reflink"))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  ## ---- 요법 → mrgsolve 이벤트 -----------------------------------------------
  build <- reactive({
    tdx <- input$tdx*DAY
    e <- NULL
    add <- function(e, x) if (is.null(e)) x else e + x
    if ("VAN" %in% input$drugs) e <- add(e, ev_van(tdx, input$wk_backbone))
    if ("LVX" %in% input$drugs) e <- add(e, ev_lvx(tdx, input$wk_backbone))
    if ("DAP" %in% input$drugs) e <- add(e, ev_dap(tdx, input$wk_backbone))
    if ("LZD" %in% input$drugs) e <- add(e, ev_lzd(tdx, input$wk_backbone))
    if ("CFZ" %in% input$drugs) e <- add(e, ev_cfz(tdx, input$wk_backbone))
    if ("RIF" %in% input$drugs && input$wk_rif > 0) e <- add(e, ev_rif(tdx, input$wk_rif))
    if (is.null(e)) e <- ev(time = 0, cmt = "NP", amt = 0)

    p <- list(MSSA = as.numeric(input$organism),
              IMSUP = input$imsup,
              IMPL = as.numeric(input$implant),
              INOC0 = input$inoc)
    if (input$surg != "none") {
      p$TSURG1 <- tdx
      p$LOGK1  <- input$logk
      p$EXCH1  <- ifelse(input$surg == "dair", 0, 1)
    }
    if (input$surg == "two") {
      p$TSURG2 <- tdx + input$reimp*7*DAY
      p$LOGK2  <- 1.5
      p$EXCH2  <- 1
      p$SPCF0  <- 0.63*input$spacer
      p$SPCS0  <- 0.37*input$spacer
    }
    list(ev = e, par = p, tdx = tdx)
  })

  sim <- eventReactive(input$go, {
    b <- build()
    m <- mrgsolve::param(pji_mod, b$par)
    out <- mrgsolve::mrgsim(m, events = b$ev, end = input$tend*DAY, delta = 4) %>%
      as.data.frame()
    out$day <- out$time/DAY
    out
  }, ignoreNULL = FALSE)

  ## ---- KPI -------------------------------------------------------------------
  output$kpi_burden <- renderText(sprintf("%.2f", tail(sim()$LOGNTOT, 1)))
  output$kpi_cure   <- renderText(sprintf("%.1f%%", 100*tail(sim()$PCURE, 1)))
  output$kpi_prpob  <- renderText({
    d <- sim(); sprintf("%.3f", 1 - exp(-1e-8*min(10^d$LOGNTOT))) })
  output$kpi_icm    <- renderText(sprintf("%.0f / 10", max(sim()$ICMSC)))

  output$narrative <- renderUI({
    d <- sim(); b <- build()
    nadir <- min(10^d$LOGNTOT); pr <- 1 - exp(-1e-8*nadir)
    cure  <- tail(d$PCURE, 1)
    HTML(sprintf(
      "<div style='background:#F7F9FB;padding:12px;border-left:4px solid #2E6DA4'>
       <b>모델이 말하는 것</b><br>
       최저 균량 <b>%.2e CFU</b>에서 rpoB 변이체가 이미 존재할 확률은 <b>%.1f%%</b>입니다.
       달성한 최대 (유리 골농도 / MBEC) 비는 <b>%.3f</b>였고,
       최종 치료 성공 확률은 <b>%.1f%%</b>입니다.<br>%s</div>",
      nadir, 100*pr, max(d$RATMAX), 100*cure,
      if (!("RIF" %in% input$drugs))
        "<span style='color:#C0392B'>리팜피신이 빠져 있습니다 — 어떤 전신 약제도 MBEC의 2%를 넘지 못하므로 바이오필름은 제거되지 않습니다.</span>"
      else if (length(setdiff(input$drugs, "RIF")) == 0)
        "<span style='color:#C0392B'>리팜피신 단독요법입니다 — 위 확률만큼의 환자에서 rpoB 변이체가 선택되어 실패합니다.</span>"
      else if (input$surg == "none")
        "<span style='color:#C0392B'>변연절제가 없습니다 — N ~ 1e10에서 변이체 존재 확률은 사실상 1이며, 리팜피신은 시작 시점에 이미 무력합니다.</span>"
      else "<span style='color:#2E7D32'>변연절제 + 리팜피신 병용 — 이 모델이 성립한다고 보는 유일한 조합입니다.</span>"))
  })

  ## ---- 1. 개요 ---------------------------------------------------------------
  output$p_overview <- renderPlot({
    d <- sim()
    dd <- d %>% select(day, `균량 log10 CFU` = LOGNTOT, `CRP (mg/L)` = CRP,
                       `관절액 WBC (/uL)` = SYNWBC, `골량 (% 기저)` = BONEV,
                       `이완 지수` = LOOSEN, `치료 성공 확률` = PCURE) %>%
      pivot_longer(-day)
    ggplot(dd, aes(day, value)) +
      geom_line(linewidth = 0.9, colour = "#2E6DA4") +
      geom_vline(xintercept = input$tdx, linetype = 2, colour = "#C0392B") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "일 (day)", y = NULL,
           title = "PJI 경과 요약 (빨간 점선 = 수술/치료 시작)") + THEME
  })

  ## ---- 2. PK -----------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim()
    pl <- d %>% select(day, VAN = CPVAN, RIF = CPRIF, LVX = CPLVX,
                       DAP = CPDAP, LZD = CPLZD, CFZ = CPCFZ) %>%
      pivot_longer(-day, names_to = "drug", values_to = "conc") %>% mutate(site = "혈장 (총)")
    bn <- d %>% select(day, VAN = VANB, RIF = RIFB, LVX = LVXB,
                       DAP = DAPB, LZD = LZDB, CFZ = CFZB) %>%
      pivot_longer(-day, names_to = "drug", values_to = "conc") %>% mutate(site = "골 (유리)")
    keep <- union(input$drugs, "RIF")
    ggplot(bind_rows(pl, bn) %>% filter(drug %in% keep, conc > 1e-4),
           aes(day, conc, colour = drug, linetype = site)) +
      geom_line(linewidth = 0.8) + scale_y_log10() +
      scale_colour_manual(values = DRUGCOL) +
      scale_linetype_manual(values = c("혈장 (총)" = 2, "골 (유리)" = 1)) +
      labs(x = "일", y = "농도 (mg/L, log)", colour = "약제", linetype = "부위",
           title = "혈장 총농도 대 유리 골농도") + THEME
  })

  ## ---- 3. 골농도 / MBEC ------------------------------------------------------
  output$p_ratio <- renderPlot({
    d <- sim()
    r <- d %>% select(day, VAN = RATVAN, RIF = RATRIF, LVX = RATLVX,
                      DAP = RATDAP, LZD = RATLZD, CFZ = RATCFZ) %>%
      pivot_longer(-day, names_to = "drug", values_to = "ratio") %>%
      filter(drug %in% union(input$drugs, "RIF"), ratio > 1e-6)
    ggplot(r, aes(day, ratio, colour = drug)) +
      geom_hline(yintercept = 1, linetype = 2, colour = "#C0392B", linewidth = 0.9) +
      annotate("text", x = Inf, y = 1.35, hjust = 1.05, label = "MBEC (= 1)",
               colour = "#C0392B", size = 4) +
      geom_line(linewidth = 0.9) + scale_y_log10() +
      scale_colour_manual(values = DRUGCOL) +
      labs(x = "일", y = "C_bone,free / MBEC (log)", colour = "약제",
           title = "달성 가능한 유리 골농도 / 바이오필름 제거 농도") + THEME
  })

  output$t_ratio <- DT::renderDataTable({
    p <- as.list(mrgsolve::param(pji_mod))
    reg <- data.frame(
      약제  = c("반코마이신 1 g q12h", "리팜피신 450 mg q12h (유도 후)",
                "레보플록사신 750 mg qd", "답토마이신 8 mg/kg qd",
                "리네졸리드 600 mg q12h", "세파졸린 2 g q8h"),
      AUC24 = c(2*1000/p$CLVAN, 2*450*p$FRIF/(2*p$CLRIF0), 750*p$FLVX/p$CLLVX,
                560/p$CLDAP, 2*600*p$FLZD/p$CLLZD0, 3*2000/p$CLCFZ),
      f_bone= c(p$PENVAN*p$FUVAN, p$PENRIF*p$FURIF, p$PENLVX*p$FULVX,
                p$PENDAP*p$FUDAP, p$PENLZD*p$FULZD, p$PENCFZ*p$FUCFZ),
      MIC   = c(p$MICVAN, p$MICRIF, p$MICLVX, p$MICDAP, p$MICLZD, p$MICCFZ),
      MBEC  = c(p$MBCVAN, p$MBCRIF, p$MBCLVX, p$MBCDAP, p$MBCLZD, p$MBCCFZ))
    reg$Cbone_free <- round(reg$AUC24/24*reg$f_bone, 3)
    reg$`비 / MIC`  <- round(reg$Cbone_free/reg$MIC, 1)
    reg$`비 / MBEC` <- round(reg$Cbone_free/reg$MBEC, 4)
    reg$AUC24 <- round(reg$AUC24, 0); reg$f_bone <- round(reg$f_bone, 3)
    DT::datatable(reg, rownames = FALSE, options = list(dom = "t", pageLength = 10)) %>%
      DT::formatStyle("비 / MBEC",
        backgroundColor = DT::styleInterval(c(0.02, 0.10), c("#FDECEC", "#FDF3E7", "#EAF6F0")))
  })

  ## ---- 4. 아집단 -------------------------------------------------------------
  output$p_pops <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(day,
      `유주균 (planktonic)`    = NP,  `바이오필름 (biofilm)` = NB,
      `지속균 (persister)`     = NPER,`SCV`                  = NSCV,
      `세포내 (intracellular)` = NIC,
      `rpoB 내성` = RP + RB,           `gyrA 내성` = QP + QB) %>%
      pivot_longer(-day, names_to = "pop", values_to = "n") %>%
      mutate(n = pmax(n, 1e-3))
    g <- ggplot(dd, aes(day, n, colour = pop)) +
      geom_line(linewidth = 0.85) +
      geom_vline(xintercept = input$tdx, linetype = 2, colour = "#555") +
      scale_colour_manual(values = POPCOL) +
      labs(x = "일", y = "CFU", colour = NULL,
           title = "세균 아집단 동태 (점선 = 수술)") + THEME
    if (input$logscale) g <- g + scale_y_log10(limits = c(1e-2, 1e11))
    g
  })

  ## ---- 5. 바이오필름 ---------------------------------------------------------
  output$p_biofilm <- renderPlot({
    d <- sim()
    dd <- d %>% select(day, `EPS 기질 (0-1)` = EPS, `EPS 성숙도` = EPSMAT,
                       `agr 신호 AIP` = AIP, `부골 SEQ` = SEQ,
                       `MDSC` = MDSC, `바이오필름 균량 log10` = LOGNB) %>%
      pivot_longer(-day)
    ggplot(dd, aes(day, value)) + geom_line(linewidth = 0.9, colour = "#4E9E7F") +
      geom_vline(xintercept = input$tdx, linetype = 2, colour = "#C0392B") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "일", y = NULL, title = "바이오필름 기질 성숙과 그 부산물") + THEME
  })

  ## ---- 6. 국소 용출 ----------------------------------------------------------
  output$p_local <- renderPlot({
    d <- sim(); p <- as.list(mrgsolve::param(pji_mod))
    ggplot(d, aes(day, pmax(LOCJ, 1e-3))) +
      geom_line(linewidth = 1, colour = DRUGCOL[["LOCAL"]]) +
      geom_hline(yintercept = p$MBCVAN, linetype = 2, colour = "#C0392B") +
      annotate("text", x = Inf, y = p$MBCVAN*1.6, hjust = 1.05,
               label = "반코마이신 MBEC 512 mg/L", colour = "#C0392B") +
      geom_hline(yintercept = p$MICVAN, linetype = 3, colour = "#2E7D32") +
      scale_y_log10() +
      labs(x = "일", y = "관절액 항생제 농도 (mg/L, log)",
           title = "항생제 함유 시멘트 스페이서의 국소 용출") + THEME
  })

  ## ---- 7. 내성 --------------------------------------------------------------
  output$p_res <- renderPlot({
    d <- sim()
    dd <- d %>% transmute(day,
      `감수성 총합` = pmax(NP + NB + NPER + NSCV + NIC, 1e-3),
      `rpoB 변이`   = pmax(RP + RB, 1e-3),
      `gyrA 변이`   = pmax(QP + QB, 1e-3)) %>%
      pivot_longer(-day, names_to = "clone", values_to = "n")
    ggplot(dd, aes(day, n, colour = clone)) +
      geom_line(linewidth = 0.95) +
      geom_hline(yintercept = 1, linetype = 3, colour = "#888") +
      geom_vline(xintercept = input$tdx, linetype = 2, colour = "#555") +
      scale_y_log10(limits = c(1e-2, 1e11)) +
      scale_colour_manual(values = c("감수성 총합" = "#2E6DA4",
                                     "rpoB 변이" = "#C0392B", "gyrA 변이" = "#7F8C8D")) +
      labs(x = "일", y = "CFU (log)", colour = NULL,
           title = "클론별 동태 — 내성은 선택되는 것이지 만들어지는 것이 아니다") + THEME
  })

  output$p_mutsupply <- renderPlot({
    d <- sim()
    nad <- min(10^d$LOGNTOT); pk <- max(10^d$LOGNTOT)
    N <- 10^seq(3, 11, by = 0.05)
    df <- data.frame(N = N, P = 1 - exp(-1e-8*N))
    ggplot(df, aes(N, P)) + geom_line(linewidth = 1, colour = "#B08040") +
      geom_vline(xintercept = pk,  linetype = 2, colour = "#C0392B") +
      geom_vline(xintercept = max(nad, 1), linetype = 2, colour = "#2E7D32") +
      annotate("text", x = pk, y = 0.15, hjust = 1.1,
               label = sprintf("최대 부하 %.1e", pk), colour = "#C0392B") +
      annotate("text", x = max(nad, 1), y = 0.85, hjust = -0.05,
               label = sprintf("최저 부하 %.1e", nad), colour = "#2E7D32") +
      scale_x_log10() +
      labs(x = "균량 N (CFU, log)", y = "P(사전 rpoB 변이체) = 1 − e^(−μN)",
           title = "μ = 1e−8 — 수술이 실제로 바꾸는 것") + THEME
  })

  ## ---- 8. 면역 --------------------------------------------------------------
  output$p_immune <- renderPlot({
    d <- sim()
    dd <- d %>% select(day, `호중구 (/uL)` = PMN, `단핵구 (/uL)` = MONO,
                       `MDSC` = MDSC, `대식세포` = MAC,
                       `IL-1beta` = IL1B, `TNF-alpha` = TNFA,
                       `IL-6` = IL6, `IL-10` = IL10) %>% pivot_longer(-day)
    ggplot(dd, aes(day, value)) + geom_line(linewidth = 0.85, colour = "#5E86BF") +
      geom_vline(xintercept = input$tdx, linetype = 2, colour = "#C0392B") +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "일", y = NULL, title = "선천면역 세포와 사이토카인") + THEME
  })

  ## ---- 9. 골 ---------------------------------------------------------------
  output$p_bone <- renderPlot({
    d <- sim()
    dd <- d %>% select(day, `RANKL` = RANKL, `OPG` = OPG,
                       `골모세포` = OBL, `파골세포` = OCL,
                       `임플란트 주위 골량 (%)` = BONEV, `이완 지수 (0-100)` = LOOSEN) %>%
      pivot_longer(-day)
    ggplot(dd, aes(day, value)) + geom_line(linewidth = 0.9, colour = "#A98A5F") +
      geom_vline(xintercept = input$tdx, linetype = 2, colour = "#C0392B") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "일", y = NULL, title = "RANKL/OPG 축과 임플란트 주위 골용해") + THEME
  })

  ## ---- 10. 바이오마커 -------------------------------------------------------
  output$p_biomarker <- renderPlot({
    d <- sim()
    dd <- d %>% select(day, `CRP (mg/L)` = CRP, `ESR (mm/h)` = ESR,
                       `관절액 WBC (/uL)` = SYNWBC, `다형핵구 (%)` = PMNPCT,
                       `alpha-defensin (S/CO)` = ADEF, `ICM 점수` = ICMSC) %>%
      pivot_longer(-day)
    thr <- data.frame(name = c("CRP (mg/L)", "ESR (mm/h)", "관절액 WBC (/uL)",
                               "다형핵구 (%)", "alpha-defensin (S/CO)", "ICM 점수"),
                      y = c(10, 30, 3000, 80, 1.0, 6))
    ggplot(dd, aes(day, value)) + geom_line(linewidth = 0.9, colour = "#7592AE") +
      geom_hline(data = thr, aes(yintercept = y), linetype = 2, colour = "#C0392B") +
      geom_vline(xintercept = input$tdx, linetype = 3, colour = "#555") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "일", y = NULL, title = "진단 바이오마커와 2018 ICM 점수") + THEME
  })

  ## ---- 11. 시나리오 비교 ----------------------------------------------------
  allout <- eventReactive(input$runall, {
    withProgress(message = "21개 시나리오 실행 중...", value = 0, {
      res <- lapply(seq_along(scen), function(i) {
        incProgress(1/length(scen), detail = names(scen)[i])
        d <- run_scen(pji_mod, scen[[i]], end = input$tend*DAY)
        d$scenario <- names(scen)[i]; d
      })
      bind_rows(res)
    })
  })

  output$t_scen <- DT::renderDataTable({
    s <- summarise_scen(allout())
    DT::datatable(s, rownames = FALSE,
                  options = list(pageLength = 21, scrollX = TRUE)) %>%
      DT::formatStyle("P_cure",
        backgroundColor = DT::styleInterval(c(0.05, 0.80), c("#FDECEC", "#FDF3E7", "#EAF6F0")))
  })

  output$p_scen <- renderPlot({
    d <- allout()
    ggplot(d, aes(time/DAY, LOGNTOT, colour = scenario)) +
      geom_line(linewidth = 0.7, show.legend = FALSE) +
      facet_wrap(~scenario, ncol = 4) +
      labs(x = "일", y = "총 균량 (log10 CFU)",
           title = "21개 표준 시나리오의 균량 경과") + THEME +
      theme(strip.text = element_text(size = 8))
  })

  ## ---- 12. 독성 -------------------------------------------------------------
  output$p_tox <- renderPlot({
    d <- sim()
    dd <- d %>% select(day, `혈청 크레아티닌 (mg/dL)` = SCR,
                       `VAN AUC24 (mg*h/L)` = AUCF,
                       `혈소판 (1e9/L)` = PLT, `ALT (U/L)` = ALT) %>% pivot_longer(-day)
    thr <- data.frame(name = c("혈청 크레아티닌 (mg/dL)", "VAN AUC24 (mg*h/L)",
                               "혈소판 (1e9/L)", "ALT (U/L)"),
                      y = c(1.35, 600, 100, 75))
    ggplot(dd, aes(day, value)) + geom_line(linewidth = 0.9, colour = "#CB8B62") +
      geom_hline(data = thr, aes(yintercept = y), linetype = 2, colour = "#C0392B") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "일", y = NULL, title = "약물 독성 지표 (점선 = 경고 문턱)") + THEME
  })

  ## ---- 13. 파라미터 ---------------------------------------------------------
  output$t_par <- DT::renderDataTable({
    p <- as.list(mrgsolve::param(pji_mod))
    DT::datatable(data.frame(parameter = names(p), value = unlist(p)),
                  rownames = FALSE, options = list(pageLength = 25, scrollX = TRUE))
  })

  output$reflink <- renderUI(HTML(
    "<p>파라미터 출처와 검증 기준점은 <code>pji_references.md</code>(68편, PMID 전수 대조)와
     <code>pji_mrgsolve_model.R</code> 머리말의 VERIFICATION 절을 참조하십시오.
     기계론적 지도는 <code>pji_qsp_model.svg</code>입니다.</p>
     <p style='color:#C0392B'><b>면책</b>: 교육·연구 목적의 반정량 모델입니다.
     임상 의사결정·처방·규제 제출에 사용하지 마십시오.</p>"))
}

shinyApp(ui, server)
