## ===========================================================================
##  cmv_shiny_app.R
##  Interactive dashboard for the CMV-in-transplant QSP model
##  이식 후 CMV 감염 QSP 모델 — Shiny 대시보드
##
##  DESIGN PRINCIPLE
##  ---------------------------------------------------------------------------
##  The model has exactly two thresholds and one clock, and every tab in this app
##  exists to show a quantity RELATIVE TO ONE OF THEM.  Nothing is plotted as a
##  bare trajectory when it can be plotted as a margin:
##
##    e*  = r0/(r0+DELI) = 0.667   the fraction of virion production a regimen
##                                 must remove for the exponent to change sign
##    E8* = r0/KE8       = 4.81    the CMV-specific CD8 count at which control
##                                 needs no drug at all -- the real endpoint
##    clock  2^(dt/1.2 d)          the fold-rise between PCR draws, which turns a
##                                 nominal trigger into an effective one
##
##  So the PK tab draws a horizontal line at e* and shades the region below it;
##  the immunology tab draws E8* and shades below it; and the strategy tab
##  displays the effective trigger next to the nominal one.  A user who never
##  reads the documentation should still leave knowing which side of each line
##  their regimen sits on.
##
##  TABS
##    1  환자 · 이식 프로파일   patient / transplant profile and risk stratification
##    2  약물 PK              antiviral PK, GCV-TP activation, free fractions
##    3  바이러스 역학          the two DNA species, log10 DNAemia, doubling time
##    4  두 역치 (e* / E8*)     THE CORE TAB: both margins on one screen
##    5  면역 재구성            CD8 effector/memory, CD4, NK, antibody, ISI
##    6  내성                  strain competition, selection window, cross-resistance
##    7  독성 (골수 · 신장)      Friberg ANC chain, eGFR, magnesium, dose log
##    8  전략 비교              prophylaxis vs pre-emptive, monitoring interval
##    9  임상 엔드포인트         disease / rejection hazards, cost, PCR count
##   10  바이오마커 패널         everything a clinic would actually order
##   11  민감도 · 훈련용량 곡선   parameter sweeps incl. the training-dose U-curve
##   12  모델 문서 · 한계        equations, calibration ledger, stated omissions
##
##  RUN
##    shiny::runApp("cmv_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  and cmv_mrgsolve_model.R in the same directory (the driver functions in its
##  trailing comment block must be sourced -- see load_engine() below).
## ===========================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

## ---------------------------------------------------------------------------
##  ENGINE
## ---------------------------------------------------------------------------
MODEL_FILE <- "cmv_mrgsolve_model.R"

load_engine <- function() {
  mod <- mread_cache("cmv", MODEL_FILE)
  ## the R driver (cmv_sim / cmv_scenario / cmv_summary / init_cmv) lives in the
  ## trailing /* ... */ block of the model file so that mread() ignores it.
  src <- readLines(MODEL_FILE)
  a <- grep("^/\\*$", src)
  b <- grep("^\\*/$", src)
  if (length(a) && length(b)) eval(parse(text = paste(
    src[(max(a) + 1):(max(b) - 1)], collapse = "\n")), envir = globalenv())
  mod
}
MOD <- load_engine()
P0  <- as.list(param(MOD))

LN2     <- log(2)
r0_of   <- function(p) LN2 / p$DBL0
EPSTAR  <- function(p) r0_of(p) / (r0_of(p) + LN2 / p$THALFTX)
E8STAR  <- function(p) r0_of(p) / p$KE8
FOLDRISE <- function(dt, p) 2 ^ (dt / p$DBL0)

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "grey92"))

## algebraic effect calculators, identical to $GLOBAL in the model ------------
epol_r <- function(gtp, fk = 1, fosf = 0, ctp = 0, p = P0) {
  x <- (gtp * fk) / p$EC50GTP
  eg <- if (x > 0) x^p$HGTP / (1 + x^p$HGTP) else 0
  ef <- if (fosf > 0) fosf / (fosf + p$EC50FOS) else 0
  ec <- if (ctp > 0) ctp / (ctp + p$EC50CPP) else 0
  1 - (1 - eg) * (1 - ef) * (1 - ec)
}
epack_r <- function(ltvf = 0, mbvf = 0, resl = 1, resm = 1, p = P0) {
  el <- if (ltvf > 0) ltvf / (ltvf + p$EC50LTV * resl) else 0
  em <- if (mbvf > 0) mbvf / (mbvf + p$EC50MBV * resm) else 0
  1 - (1 - el) * (1 - em)
}
gcv_ss <- function(mg_per_day, crcl, p = P0) {
  cl <- (p$CLGCV * (p$FRENGCV * crcl / 100 + (1 - p$FRENGCV))) * 24
  cavg <- mg_per_day / p$MWVGCV * 1000 * p$FVGCV / cl
  c(cavg = cavg, gtp = p$KPHOS * cavg / p$KDEGTP)
}

## ---------------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel(HTML(paste0(
    "<b>이식 후 거대세포바이러스(CMV) 감염 — QSP 대시보드</b>",
    "<br><span style='font-size:14px'>Cytomegalovirus in transplant ",
    "recipients · 48-ODE model · two thresholds and one clock</span>"))),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("환자 (patient)"),
      selectInput("sero", "공여자/수여자 혈청상태 (D/R serostatus)",
                  c("D+/R−  (highest risk, primary infection)" = "D+/R-",
                    "D+/R+" = "D+/R+", "D−/R+" = "D-/R+",
                    "D−/R−  (lowest risk)" = "D-/R-"), "D+/R-"),
      selectInput("induction", "유도요법 (induction)",
                  c("rabbit ATG" = "ATG", "alemtuzumab" = "alemtuzumab",
                    "basiliximab / none" = "basiliximab"), "ATG"),
      sliderInput("crcl", "이식신 기능 CrCl (mL/min)", 15, 120, 95, 5),
      checkboxInput("mpa", "mycophenolate", TRUE),
      checkboxInput("mtor", "mTOR inhibitor 전환 (protective)", FALSE),

      hr(), h4("전략 (strategy)"),
      radioButtons("strategy", NULL,
                   c("universal prophylaxis" = "prophylaxis",
                     "pre-emptive therapy" = "preemptive",
                     "none (natural history)" = "none"), "prophylaxis"),
      conditionalPanel(
        "input.strategy == 'prophylaxis'",
        selectInput("proph_drug", "예방약 (prophylactic agent)",
                    c("valganciclovir 900 mg od" = "VGCV",
                      "letermovir 480 mg od" = "LTV"), "VGCV"),
        sliderInput("proph_days", "예방기간 (days)", 0, 365, 200, 10),
        sliderInput("adhere", "복약순응도 / 유효용량 비율 (adherence)",
                    0.2, 1.0, 1.0, 0.05),
        checkboxInput("post_monitor", "예방 종료 후 PCR 감시 (surveillance)", TRUE),
        checkboxInput("tac_cut_ltv",
                      "letermovir 시작 시 tacrolimus 용량 절반 감량", TRUE)),
      sliderInput("monitor_int", "PCR 감시 간격 Δt (days)", 2, 21, 14, 0.5),
      numericInput("trigger", "치료 개시 역치 (IU/mL)", 1000, 137, 100000, 100),

      hr(), h4("치료 · 구제 (treatment / salvage)"),
      selectInput("tx_drug", "1차 치료 (first-line)",
                  c("valganciclovir" = "VGCV", "ganciclovir IV" = "VGCV",
                    "maribavir" = "MBV", "foscarnet" = "FOS"), "VGCV"),
      selectInput("tx_drug2", "불응성 시 전환 (on refractory CMV)",
                  c("none (continue)" = "VGCV", "maribavir" = "MBV",
                    "foscarnet" = "FOS", "cidofovir" = "CDV",
                    "GCV + maribavir (antagonistic!)" = "GCV+MBV"), "MBV"),
      sliderInput("resist_seed",
                  "치료 시작 시 UL97 소수변이 비율 (minority variant)",
                  0, 0.5, 0, 0.05),
      sliderInput("tx_min_days", "2차 예방 최소 기간 (secondary prophylaxis, d)",
                  0, 90, 0, 7),
      checkboxInput("renal_adjust", "신기능에 따른 용량 조절 (renal dose bands)",
                    TRUE),
      checkboxInput("gcsf", "호중구감소 시 G-CSF", TRUE),

      hr(),
      sliderInput("tend", "관찰기간 (days)", 100, 730, 365, 5),
      actionButton("go", "시뮬레이션 실행 (run)", class = "btn-primary"),
      br(), br(),
      helpText(HTML(paste0(
        "<b>e* = 0.667</b> — 이 값을 넘지 못하는 처방은 라벨과 무관하게 ",
        "돌파감염을 낸다.<br><b>E8* = 4.81 /µL</b> — 약 없이 조절이 되는 ",
        "CMV 특이 CD8 수. 어떤 약도 이 값을 옮기지 못한다.")))
    ),

    mainPanel(
      width = 9,
      fluidRow(
        column(3, wellPanel(strong("e (총 억제율)"), br(),
                            textOutput("kpi_eps", inline = TRUE))),
        column(3, wellPanel(strong("e − e* (여유)"), br(),
                            textOutput("kpi_margin", inline = TRUE))),
        column(3, wellPanel(strong("유효 역치 (effective trigger)"), br(),
                            textOutput("kpi_trigger", inline = TRUE))),
        column(3, wellPanel(strong("E8 − E8* (면역 여유)"), br(),
                            textOutput("kpi_imm", inline = TRUE)))
      ),
      tabsetPanel(
        id = "tabs", type = "tabs",

        tabPanel("1 · 환자 프로파일",
                 h4("위험 계층화 (risk stratification)"),
                 DTOutput("tbl_profile"),
                 h4("초기조건 (initial conditions actually used)"),
                 verbatimTextOutput("txt_init"),
                 helpText(paste(
                   "D+/R− is the highest-risk group in SOLID-ORGAN transplant",
                   "because the recipient has no memory pool at all (EM8 = 0).",
                   "After HSCT the ordering reverses: the RECIPIENT-seropositive",
                   "patient is the one at risk, because the reservoir is the",
                   "host's and the new immune system is the donor's."))),

        tabPanel("2 · 약물 PK",
                 h4("혈장 농도 (plasma concentrations)"),
                 plotOutput("p_pk", height = 300),
                 h4("간시클로버 활성화 단계 (the pUL97-dependent activation step)"),
                 plotOutput("p_gtp", height = 280),
                 helpText(HTML(paste0(
                   "GCV-TP는 <b>바이러스</b> 키나아제 pUL97이 만든다. 따라서 ",
                   "maribavir(pUL97 억제제)는 ganciclovir의 <b>효과항</b>이 아니라 ",
                   "<b>활성화항</b>에 들어가고, 두 약은 구조적으로 길항한다."))),
                 h4("자유약물 분율 (free fraction matters: 99% / 98% bound)"),
                 DTOutput("tbl_free")),

        tabPanel("3 · 바이러스 역학",
                 h4("측정 DNAemia는 두 종류의 합이다"),
                 plotOutput("p_species", height = 330),
                 helpText(HTML(paste0(
                   "Vv = 캡시드에 포장된 virion DNA(유일한 감염성 종). ",
                   "Vl = 죽는 감염세포에서 나오는 유리 DNA(PCR이 보는 것의 약 ",
                   "67%). 중합효소 억제제는 둘 다 줄이고, terminase/kinase ",
                   "억제제는 Vv만 줄이며 concatemer가 절단되지 않아 세포당 DNA를 ",
                   "오히려 늘린다 ⇒ <b>letermovir가 완벽하게 작동해도 혈장 CMV ",
                   "DNA는 거의 떨어지지 않는다.</b>"))),
                 h4("배가시간과 감소반감기 (the two measured rates)"),
                 plotOutput("p_rates", height = 260)),

        tabPanel("4 · 두 역치 (e* / E8*)",
                 h4("역치 1 — 약물: e 가 e* 위에 있는가"),
                 plotOutput("p_eps", height = 300),
                 h4("역치 2 — 면역: E8 가 E8* 에 도달했는가"),
                 plotOutput("p_e8", height = 300),
                 h4("모든 처방을 e* 하나로 채점한 표"),
                 DTOutput("tbl_eps_table")),

        tabPanel("5 · 면역 재구성",
                 h4("CMV 특이 T세포 · NK · 항체"),
                 plotOutput("p_imm", height = 330),
                 h4("면역억제 총량 ISI 와 그 구성"),
                 plotOutput("p_isi", height = 280),
                 helpText(paste(
                   "The expansion term is proportional to ANTIGEN, which is the",
                   "infected-cell pool the antiviral removed. This coupling is",
                   "the whole of late-onset CMV: clearing threshold 1 delays",
                   "reaching threshold 2."))),

        tabPanel("6 · 내성",
                 h4("균주 경쟁 (strain competition in virion DNA)"),
                 plotOutput("p_strain", height = 300),
                 h4("선택압 창 (days inside the selection window)"),
                 plotOutput("p_selwin", height = 250),
                 h4("교차내성 표 (cross-resistance)"),
                 DTOutput("tbl_xres"),
                 helpText(HTML(paste0(
                   "내성은 <b>세기</b>가 아니라 <b>모양</b>이다. UL56 C325Y ",
                   "한 개의 돌연변이는 letermovir의 e를 0.969 → 0.010 (효과의 ",
                   "1.1%만 남음)으로 만들고, UL97 M460V는 ganciclovir의 e를 ",
                   "0.890 → 0.112 (12.6% 남음)로 만든다. 둘 다 '한 번의 ",
                   "돌연변이'지만 증량으로 만회할 여지가 있는 것은 후자뿐이며, ",
                   "그마저도 필요한 GCV-TP가 허가 치료용량의 2.0배다.")))),

        tabPanel("7 · 독성 (골수 · 신장)",
                 h4("Friberg 5-구획 호중구 사슬"),
                 plotOutput("p_anc", height = 300),
                 h4("신기능 · 마그네슘"),
                 plotOutput("p_kidney", height = 280),
                 h4("용량조절 기록 (what the ANC rule actually did)"),
                 verbatimTextOutput("txt_log"),
                 h4("CrCl → 노출 → 효과 · 독성 (the loop, tabulated)"),
                 DTOutput("tbl_renal")),

        tabPanel("8 · 전략 비교",
                 h4("감시간격은 용량이다 (monitoring interval as a dose)"),
                 DTOutput("tbl_clock"),
                 h4("전략별 비교 (run the standard arms side by side)"),
                 actionButton("run_arms", "표준 16개 시나리오 실행",
                              class = "btn-default"),
                 br(), br(),
                 DTOutput("tbl_arms"),
                 plotOutput("p_arms", height = 320)),

        tabPanel("9 · 임상 엔드포인트",
                 h4("누적 위험 (cumulative hazards)"),
                 plotOutput("p_hz", height = 300),
                 h4("요약 (endpoint summary for this run)"),
                 DTOutput("tbl_end"),
                 helpText(paste(
                   "P(disease) and P(rejection) are hazard-function",
                   "probabilities for ONE deterministic virtual patient. They",
                   "are not trial incidence rates and must not be quoted as",
                   "such."))),

        tabPanel("10 · 바이오마커 패널",
                 h4("임상에서 실제로 주문하는 항목들"),
                 plotOutput("p_biomk", height = 520),
                 DTOutput("tbl_biomk")),

        tabPanel("11 · 민감도 · 훈련용량 곡선",
                 h4("훈련용량 U-곡선 (the training-dose U-curve)"),
                 helpText(HTML(paste0(
                   "완전 억제가 최적인가? E8 증식이 항원 의존적이라면 완전히 ",
                   "억제하는 예방은 지속적인 것을 아무것도 남기지 않는다. ",
                   "모델은 200일 예방에서 <b>e* 아래</b>인 순응도 0.40 부근에 ",
                   "최소값을 낸다. <b>처방 권고로 읽지 말 것</b> — 단일 결정론적 ",
                   "환자이며 최소값의 깊이는 KAG에 크게 좌우된다."))),
                 actionButton("run_ucurve", "U-곡선 계산 (slow)",
                              class = "btn-default"),
                 br(), br(),
                 plotOutput("p_ucurve", height = 340),
                 h4("일변량 민감도 (one-at-a-time sensitivity)"),
                 selectInput("sens_par", "파라미터",
                             c("DBL0", "THALFTX", "KE8", "KAG", "KISE", "RHO8",
                               "EC50GTP", "FK_A", "RESBLTV", "AMPPK", "GAM",
                               "EMAXMYE", "KDIS", "LAMT"), "KAG"),
                 plotOutput("p_sens", height = 320)),

        tabPanel("12 · 모델 문서 · 한계",
                 h4("조직 논지 (organising thesis)"),
                 verbatimTextOutput("txt_thesis"),
                 h4("보정 원장 (calibration ledger: what was fitted, to what)"),
                 DTOutput("tbl_ledger"),
                 h4("의도적으로 넣지 않은 것 (stated omissions)"),
                 verbatimTextOutput("txt_gaps"),
                 h4("실행 중 발견되어 수정된 결함 (defects found by running it)"),
                 verbatimTextOutput("txt_defects"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
##  SERVER
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  ## ---- run the engine -----------------------------------------------------
  R <- eventReactive(input$go, ignoreNULL = FALSE, {
    withProgress(message = "integrating 48 ODEs ...", {
      cmv_sim(MOD,
              strategy     = input$strategy,
              proph_drug   = input$proph_drug %||% "VGCV",
              proph_days   = input$proph_days %||% 200,
              monitor_int  = input$monitor_int,
              trigger      = input$trigger,
              tx_drug      = input$tx_drug,
              tx_drug2     = input$tx_drug2,
              tx_min_days  = input$tx_min_days,
              post_monitor = input$post_monitor %||% TRUE,
              renal_adjust = input$renal_adjust,
              tac_cut_ltv  = input$tac_cut_ltv %||% TRUE,
              adhere       = input$adhere %||% 1,
              resist_seed  = input$resist_seed,
              gcsf_policy  = input$gcsf,
              mpa = as.numeric(input$mpa), mtor = as.numeric(input$mtor),
              sero = input$sero, induction = input$induction,
              crcl = input$crcl, tend = input$tend)
    })
  })
  `%||%` <- function(a, b) if (is.null(a)) b else a
  S <- reactive(R()$sim)

  ## ---- KPI strip ----------------------------------------------------------
  output$kpi_eps <- renderText({
    e <- max(S()$EPS_TOT, na.rm = TRUE); sprintf("%.3f (peak on-treatment)", e) })
  output$kpi_margin <- renderText({
    m <- max(S()$MARGIN, na.rm = TRUE)
    sprintf("%+.3f  %s", m, if (m > 0) "CONTROL" else "BREAKTHROUGH") })
  output$kpi_trigger <- renderText({
    f <- FOLDRISE(input$monitor_int, P0)
    sprintf("%.0f IU/mL  (= %.0f × %.1f-fold)", input$trigger * f,
            input$trigger, f) })
  output$kpi_imm <- renderText({
    s <- S(); sprintf("%+.2f /µL at day %.0f", tail(s$IMMMARG, 1),
                      tail(s$time, 1)) })

  ## ---- tab 1: profile -----------------------------------------------------
  output$tbl_profile <- renderDT({
    datatable(data.frame(
      serostatus = c("D+/R−", "D+/R+", "D−/R+", "D−/R−"),
      memory_pool_EM8 = c(0, 20, 22, 0),
      infection_type = c("primary (graft-borne)", "reactivation / superinfection",
                         "reactivation", "none unless transfused"),
      SOT_risk = c("highest", "intermediate", "intermediate", "lowest"),
      HSCT_risk = c("intermediate", "highest (recipient reservoir)",
                    "highest", "lowest"),
      check.names = FALSE), options = list(dom = "t"), rownames = FALSE)
  })
  output$txt_init <- renderPrint({
    str(init_cmv(MOD, input$sero, input$induction, input$crcl)[
      c("LAT", "IW", "E8", "EM8", "E4", "AB", "LYM", "ATGC", "GFR")])
  })

  ## ---- tab 2: PK ----------------------------------------------------------
  output$p_pk <- renderPlot({
    s <- S()
    d <- s %>% select(time, GCV, LTV, MBV, FOS) %>%
      pivot_longer(-time, names_to = "drug", values_to = "conc") %>%
      filter(conc > 1e-6)
    ggplot(d, aes(time, conc, colour = drug)) + geom_line(linewidth = 0.8) +
      scale_y_log10() + facet_wrap(~drug, scales = "free_y") +
      labs(x = "days after transplant", y = "plasma concentration (µM, log)") +
      THEME + theme(legend.position = "none")
  })
  output$p_gtp <- renderPlot({
    s <- S()
    d <- s %>% select(time, `plasma GCV (µM)` = GCV,
                      `GCV-TP (µM-eq)` = GTP) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(x = "days", y = NULL,
           title = paste("GCV-TP t½ 18 h; EC50 =", P0$EC50GTP,
                         "µM-eq, Hill 2 — a UL97 mutant divides GCV-TP by 8")) +
      THEME
  })
  output$tbl_free <- renderDT({
    datatable(data.frame(
      drug = c("ganciclovir", "letermovir", "maribavir", "foscarnet"),
      bound_pct = c("~2", "99", "98", "~17"),
      free_fraction = c(0.98, P0$FULTV, P0$FUMBV, 1.0),
      free_Cavg_uM = c(round(gcv_ss(900, input$crcl)["cavg"], 2),
                       round(480 / P0$MWLTV * 1000 * P0$FLTV /
                               (P0$CLLTV * 24) * P0$FULTV, 4),
                       round(800 / P0$MWMBV * 1000 * P0$FMBV /
                               (P0$CLMBV * 24) * P0$FUMBV, 3),
                       round(90 * 70 * 2 / 126 * 1000 / (P0$CLFOS * 24), 0)),
      EC50_uM = c(NA, P0$EC50LTV, P0$EC50MBV, P0$EC50FOS),
      check.names = FALSE), options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- tab 3: viral species ----------------------------------------------
  output$p_species <- renderPlot({
    s <- S()
    d <- s %>% mutate(Vv = VVW + VVA + VVB, Vl = VLY,
                      measured = Vv + Vl) %>%
      select(time, Vv, Vl, measured) %>% pivot_longer(-time) %>%
      mutate(value = pmax(value, 1))
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 137, linetype = 3) +
      annotate("text", x = 5, y = 165, hjust = 0, size = 3,
               label = "assay LLOQ 137 IU/mL") +
      geom_hline(yintercept = input$trigger, linetype = 2, colour = "red") +
      annotate("text", x = 5, y = input$trigger * 1.3, hjust = 0, size = 3,
               colour = "red", label = "treatment trigger") +
      scale_y_log10() +
      labs(x = "days", y = "IU/mL (log10)",
           colour = NULL, title = "Vv is infectious; Vl is what the PCR mostly sees") +
      THEME
  })
  output$p_rates <- renderPlot({
    s <- S()
    d <- s %>% mutate(l = log10(pmax(DNAEMIA, 1))) %>%
      mutate(slope = c(NA, diff(l)) / c(NA, diff(time))) %>%
      filter(is.finite(slope))
    ggplot(d, aes(time, slope * log(10))) + geom_line(linewidth = 0.7) +
      geom_hline(yintercept = c(LN2 / P0$DBL0, -LN2 / P0$THALFTX),
                 linetype = 2, colour = c("red", "blue")) +
      annotate("text", x = max(d$time) * 0.6, y = LN2 / P0$DBL0 + 0.05,
               size = 3, colour = "red",
               label = "r0 = +0.578/d  (doubling 1.2 d)") +
      annotate("text", x = max(d$time) * 0.6, y = -LN2 / P0$THALFTX - 0.05,
               size = 3, colour = "blue",
               label = "−DELI = −0.289/d  (decline t½ 2.4 d)") +
      labs(x = "days", y = "instantaneous exponent of DNAemia (1/d)") + THEME
  })

  ## ---- tab 4: THE TWO THRESHOLDS -----------------------------------------
  output$p_eps <- renderPlot({
    s <- S(); es <- EPSTAR(P0)
    ggplot(s, aes(time, EPS_TOT)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = es,
               fill = "red", alpha = 0.10) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = es, linetype = 2, colour = "red") +
      annotate("text", x = max(s$time) * 0.02, y = es + 0.03, hjust = 0,
               colour = "red", size = 3.6,
               label = sprintf("e* = %.3f — below this line the exponent is POSITIVE",
                               es)) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = "days", y = "e = total fraction of virion production removed") +
      THEME
  })
  output$p_e8 <- renderPlot({
    s <- S(); e8 <- E8STAR(P0)
    ggplot(s, aes(time, E8EFFC)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = e8,
               fill = "red", alpha = 0.10) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = e8, linetype = 2, colour = "red") +
      annotate("text", x = max(s$time) * 0.02, y = e8 * 1.15, hjust = 0,
               colour = "red", size = 3.6,
               label = sprintf("E8* = %.2f /µL — control with NO drug", e8)) +
      scale_y_log10() +
      labs(x = "days",
           y = "effective CMV-specific CD8 (E8 + 0.5·EM8, cells/µL, log)") +
      THEME
  })
  output$tbl_eps_table <- renderDT({
    es <- EPSTAR(P0)
    g1 <- gcv_ss(900, input$crcl); g2 <- gcv_ss(1800, input$crcl)
    ltvf <- 480 / P0$MWLTV * 1000 * P0$FLTV / (P0$CLLTV * 24) * P0$FULTV
    mbvf <- 800 / P0$MWMBV * 1000 * P0$FMBV / (P0$CLMBV * 24) * P0$FUMBV
    fosc <- 90 * 70 * 2 / 126 * 1000 / (P0$CLFOS * 24)
    rows <- list(
      c("valGCV 900 od / WT",        epol_r(g1["gtp"])),
      c("valGCV 900 od / UL97 mut",  epol_r(g1["gtp"], P0$FK_A)),
      c("valGCV 900 BID / WT",       epol_r(g2["gtp"])),
      c("valGCV 900 BID / UL97 mut", epol_r(g2["gtp"], P0$FK_A)),
      c("letermovir 480 od / WT",    epack_r(ltvf)),
      c("letermovir 480 od / UL56",  epack_r(ltvf, resl = P0$RESBLTV)),
      c("letermovir 240 od (+CsA)",  epack_r(ltvf / 2)),
      c("maribavir 400 BID",         epack_r(mbvf = mbvf)),
      c("foscarnet 90 q12h",         epol_r(0, 1, fosc)),
      c("GCV 900 BID + MBV 400 BID",
        1 - (1 - epol_r(g2["gtp"] / (1 + mbvf / P0$KIMBVK))) *
          (1 - epack_r(mbvf = mbvf))))
    d <- data.frame(regimen = sapply(rows, `[`, 1),
                    e = round(as.numeric(sapply(rows, `[`, 2)), 4))
    d$margin <- round(d$e - es, 3)
    d$verdict <- ifelse(d$margin > 0, "CONTROL", "BREAKTHROUGH")
    datatable(d, options = list(dom = "t", pageLength = 20), rownames = FALSE) %>%
      formatStyle("verdict", backgroundColor =
                    styleEqual(c("CONTROL", "BREAKTHROUGH"),
                               c("#e6f4ea", "#fdecea")))
  })

  ## ---- tab 5: immunity ---------------------------------------------------
  output$p_imm <- renderPlot({
    s <- S()
    d <- s %>% select(time, E8, EM8, E4, NKA, AB) %>% pivot_longer(-time) %>%
      mutate(value = pmax(value, 1e-3))
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      geom_hline(yintercept = E8STAR(P0), linetype = 2, colour = "red") +
      scale_y_log10() + labs(x = "days", y = "cells/µL or relative units (log)",
                             colour = NULL) + THEME
  })
  output$p_isi <- renderPlot({
    s <- S()
    d <- s %>% transmute(time,
                         tacrolimus = P0$WTAC * TAC / 8,
                         lymphodepletion = P0$WLYM * pmax(0, 1 - LYM / P0$LYM0),
                         steroid = P0$WSTER * STER / 20,
                         mycophenolate = P0$WMPA * as.numeric(input$mpa)) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, fill = name)) + geom_area() +
      labs(x = "days", y = "contribution to ISI", fill = NULL) + THEME
  })

  ## ---- tab 6: resistance -------------------------------------------------
  output$p_strain <- renderPlot({
    s <- S()
    d <- s %>% select(time, `wild type` = VVW, `UL97 mutant` = VVA,
                      `UL56 mutant` = VVB) %>%
      pivot_longer(-time) %>% mutate(value = pmax(value, 1e-3))
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      scale_y_log10() + labs(x = "days", y = "virion DNA (IU/mL, log)",
                             colour = NULL) + THEME
  })
  output$p_selwin <- renderPlot({
    s <- S()
    d <- s %>% select(time, `UL97 window` = SELWINA, `UL56 window` = SELWINB) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "days", y = "cumulative days inside the selection window",
           colour = NULL,
           subtitle = paste("selection = wild type held (e > e*) while the",
                            "mutant's own R_eff still exceeds 1")) + THEME
  })
  output$tbl_xres <- renderDT({
    datatable(data.frame(
      mutation = c("UL97 M460V/A594V", "UL54 V715M/A809V", "UL56 C325Y",
                   "UL97 T409M/H411Y/C480F"),
      mechanism = c("loss of GCV phosphorylation (kinase)",
                    "polymerase active site",
                    "terminase subunit", "maribavir binding site (kinase)"),
      ganciclovir = c("resistant (8×)", "resistant", "susceptible",
                      "susceptible"),
      letermovir = c("susceptible", "susceptible", "resistant (>3000×)",
                     "susceptible"),
      maribavir = c("susceptible", "susceptible", "susceptible", "resistant"),
      foscarnet = c("susceptible", "variable", "susceptible", "susceptible"),
      cidofovir = c("susceptible", "resistant", "susceptible", "susceptible"),
      modelled_as_ODE_strain = c("yes (IA)", "no — map/table only", "yes (IB)",
                                 "no — map/table only"),
      check.names = FALSE), options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- tab 7: toxicity ---------------------------------------------------
  output$p_anc <- renderPlot({
    s <- S()
    ggplot(s, aes(time, ANC)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.5,
               fill = "red", alpha = 0.12) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1.0,
               fill = "orange", alpha = 0.12) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = c(0.5, 1.0, 1.5), linetype = 3) +
      labs(x = "days", y = "ANC (10⁹/L)",
           subtitle = "<1.0 → 50% dose cut + G-CSF · <0.5 → hold · >1.5 ×2 → resume") +
      THEME
  })
  output$p_kidney <- renderPlot({
    s <- S()
    d <- s %>% select(time, `eGFR (mL/min/1.73)` = GFR,
                      `magnesium (mg/dL)` = MG,
                      `tubular injury (a.u.)` = TUBI) %>% pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "days", y = NULL) + THEME
  })
  output$txt_log <- renderPrint({ cat(R()$meta$log, sep = "\n") })
  output$tbl_renal <- renderDT({
    bands <- data.frame(crcl = c(100, 80, 60, 50, 40, 30, 20))
    bands$valGCV <- c("900 od", "900 od", "900 od", "450 od", "450 od",
                      "450 q2d", "450 2×/wk")
    mgd <- c(900, 900, 900, 450, 450, 225, 128.6)
    out <- do.call(rbind, lapply(seq_len(nrow(bands)), function(i) {
      g <- gcv_ss(mgd[i], bands$crcl[i])
      ed <- min(0.95, P0$EMAXMYE * g["cavg"] / (P0$EC50MYE + g["cavg"]) + P0$EMPA)
      data.frame(CrCl = bands$crcl[i], dose = bands$valGCV[i],
                 GCV_Cavg_uM = round(g["cavg"], 2),
                 GCV_TP = round(g["gtp"], 1),
                 e = round(epol_r(g["gtp"]), 3),
                 margin = round(epol_r(g["gtp"]) - EPSTAR(P0), 3),
                 ANC_ss = round(P0$CIRC0 * (1 - ed)^(1 / P0$GAM), 2))
    }))
    datatable(out, options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- tab 8: strategy ---------------------------------------------------
  output$tbl_clock <- renderDT({
    dt <- c(2, 3.5, 7, 10, 14, 21)
    datatable(data.frame(
      interval_days = dt,
      fold_rise_between_draws = round(FOLDRISE(dt, P0), 1),
      effective_trigger_IU_mL = round(input$trigger * FOLDRISE(dt, P0)),
      check.names = FALSE), options = list(dom = "t"), rownames = FALSE)
  })
  ARMS <- eventReactive(input$run_arms, {
    withProgress(message = "running the standard arms ...", {
      ids <- names(CMV_SCENARIOS)
      do.call(rbind, lapply(seq_along(ids), function(i) {
        incProgress(1 / length(ids), detail = ids[i])
        r <- cmv_scenario(MOD, ids[i], tend = input$tend)
        cbind(arm = ids[i], cmv_summary(r))
      }))
    })
  })
  output$tbl_arms <- renderDT({
    datatable(ARMS() %>% mutate(across(where(is.numeric), ~round(.x, 3))),
              options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
  })
  output$p_arms <- renderPlot({
    a <- ARMS()
    ggplot(a, aes(reorder(arm, P_disease), P_disease)) +
      geom_col(fill = "grey40") + coord_flip() +
      labs(x = NULL, y = "P(CMV end-organ disease) — one virtual patient") +
      THEME
  })

  ## ---- tab 9: endpoints --------------------------------------------------
  output$p_hz <- renderPlot({
    s <- S()
    d <- s %>% select(time, `P(CMV disease)` = PDIS,
                      `P(acute rejection)` = PREJ) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "days", y = "cumulative probability", colour = NULL) + THEME
  })
  output$tbl_end <- renderDT({
    datatable(cmv_summary(R()) %>% mutate(across(where(is.numeric), ~round(.x, 3))),
              options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  ## ---- tab 10: biomarker panel -------------------------------------------
  output$p_biomk <- renderPlot({
    s <- S()
    d <- s %>% transmute(time,
                         `CMV DNA (log10 IU/mL)` = LOG10VL,
                         `CMV-specific CD8 (/µL)` = E8EFFC,
                         `ANC (10⁹/L)` = ANC,
                         `eGFR` = GFR,
                         `Mg (mg/dL)` = MG,
                         `tacrolimus (ng/mL)` = TAC,
                         `lymphocytes (10³/µL)` = LYM,
                         `e − e*` = MARGIN,
                         `tissue load (log10 IU/g)` = log10(pmax(VT, 1))) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "days after transplant", y = NULL) + THEME
  })
  output$tbl_biomk <- renderDT({
    s <- S()
    idx <- sapply(c(7, 14, 30, 60, 100, 200, 365),
                  function(d) which.min(abs(s$time - d)))
    d <- s[idx, c("time", "LOG10VL", "E8EFFC", "MARGIN", "ANC", "GFR", "MG",
                  "TAC", "ISIC", "PDIS")]
    datatable(d %>% mutate(across(where(is.numeric), ~round(.x, 3))),
              options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  ## ---- tab 11: sweeps ----------------------------------------------------
  UC <- eventReactive(input$run_ucurve, {
    withProgress(message = "sweeping adherence × duration ...", {
      grid <- expand.grid(adhere = c(1, .7, .5, .45, .4, .35, .3, .25, .2),
                          proph_days = c(100, 200))
      do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
        incProgress(1 / nrow(grid))
        r <- cmv_sim(MOD, strategy = "prophylaxis", proph_drug = "VGCV",
                     proph_days = grid$proph_days[i], post_monitor = FALSE,
                     adhere = grid$adhere[i], sero = input$sero,
                     induction = input$induction, crcl = input$crcl,
                     tend = 400)
        g <- gcv_ss(900 * grid$adhere[i], input$crcl)
        data.frame(adhere = grid$adhere[i], proph_days = grid$proph_days[i],
                   e = epol_r(g["gtp"]),
                   P_disease = 1 - exp(-r$state$HZD),
                   selwin = r$state$MUTA)
      }))
    })
  })
  output$p_ucurve <- renderPlot({
    u <- UC()
    ggplot(u, aes(e, P_disease, colour = factor(proph_days))) +
      geom_line(linewidth = 0.9) + geom_point() +
      geom_vline(xintercept = EPSTAR(P0), linetype = 2, colour = "red") +
      annotate("text", x = EPSTAR(P0), y = max(u$P_disease), hjust = -0.05,
               colour = "red", size = 3.5, label = "e* = 0.667") +
      labs(x = "e achieved during prophylaxis",
           y = "P(CMV end-organ disease)", colour = "prophylaxis days",
           subtitle = paste("the minimum sits BELOW e* — partial suppression",
                            "leaves the antigen that trains E8")) + THEME
  })
  output$p_sens <- renderPlot({
    pn <- input$sens_par
    base <- P0[[pn]]
    mult <- c(0.5, 0.75, 1, 1.5, 2)
    d <- do.call(rbind, lapply(mult, function(m) {
      pp <- P0; pp[[pn]] <- base * m
      mm <- MOD %>% param(setNames(list(base * m), pn))
      r <- cmv_sim(mm, strategy = input$strategy,
                   proph_drug = input$proph_drug %||% "VGCV",
                   proph_days = input$proph_days %||% 200,
                   monitor_int = input$monitor_int, trigger = input$trigger,
                   sero = input$sero, induction = input$induction,
                   crcl = input$crcl, tend = min(input$tend, 365))
      data.frame(mult = m, value = base * m,
                 peak_log10VL = max(r$sim$LOG10VL),
                 P_disease = 1 - exp(-r$state$HZD),
                 E8_end = tail(r$sim$E8EFFC, 1))
    }))
    d2 <- d %>% pivot_longer(c(peak_log10VL, P_disease, E8_end))
    ggplot(d2, aes(mult, value)) + geom_line() + geom_point() +
      facet_wrap(~name, scales = "free_y") +
      labs(x = paste("multiple of", pn, "( base =", signif(base, 4), ")"),
           y = NULL) + THEME
  })

  ## ---- tab 12: documentation ---------------------------------------------
  output$txt_thesis <- renderText({ paste(
"r = KPROD*(1-e_pol)*(1-e_pack) - DELI - KE8*E8eff - KENK*NKA",
"",
"Two measured clinical numbers fix KPROD:",
"    untreated DNAemia doubling time  1.2 d  ->  r0   = 0.5776 /d",
"    on-treatment decline half-life   2.4 d  ->  DELI = 0.2888 /d",
"    KPROD = r0 + DELI                       =  0.8664 /d",
"",
"THRESHOLD 1  e*  = r0/KPROD = 0.6667",
"   potency, resistance, renal dose banding and drug interaction are all the",
"   same question asked of this one number.",
"THRESHOLD 2  E8* = r0/KE8   = 4.81 CMV-specific CD8 per uL",
"   control with NO drug.  The endpoint of the illness.  No drug moves it.",
"THE CLOCK    2^(dt/1.2 d) fold-rise between PCR draws",
"   q7d turns a nominal 1000 IU/mL trigger into a 57,000 IU/mL trigger.",
"THE COUPLING d(E8)/dt is proportional to ANTIGEN, i.e. to the infected cells",
"   the drug just removed -- so clearing threshold 1 delays reaching",
"   threshold 2.  That single sentence is late-onset CMV.", sep = "\n") })
  output$tbl_ledger <- renderDT({
    datatable(data.frame(
      quantity = c("KPROD, DELI", "LAMT = 2.00", "KDIS = 0.0300",
                   "GAM / EMAXMYE / EMPA", "EC50 values", "all PK parameters",
                   "AMPPK = 0.50", "KAG = 0.05", "costs"),
      how_it_was_set = c(
        "doubling time 1.2 d and on-treatment decline half-life 2.4 d",
        "ONE number fitted to the observed peak DNAemia of untreated primary D+/R- infection (10^4.83)",
        "ONE number fitted so untreated D+/R- gives P(CMV disease) = 0.60 at one year",
        "re-fitted against CHRONIC-exposure ANC; the published Friberg gamma of 0.17 puts a patient on mycophenolate alone at ANC 0.83",
        "published in-vitro EC50, corrected to free drug where protein binding is high",
        "published population PK (F, V, CL, t1/2, renal fraction)",
        "ASSUMED — no direct human measurement; swept 0-1 in the reference output",
        "WEAKEST parameter in the model; sets the depth of the training-dose U-curve",
        "illustrative list prices, order of magnitude only"),
      check.names = FALSE),
      options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })
  output$txt_gaps <- renderText({ paste(
"* UL54 polymerase mutants: in the map and the cross-resistance table, NOT a",
"  separate ODE strain.",
"* No inter-individual variability. Every probability shown is a hazard for ONE",
"  deterministic virtual patient and must not be read as a trial incidence.",
"* Rejection and eGFR are reduced-form hazard surrogates, not a mechanistic",
"  alloimmune model.",
"* Adoptive CMV-specific T-cell therapy has an input (ACT) but no dosing",
"  schedule; vaccines appear in the map only.",
"* Whole-blood vs plasma assay conversion and the 1-2 log10 inter-laboratory",
"  spread in IU/mL are not modelled.",
"* HSCT-specific engraftment kinetics and GVHD are not modelled; the HSCT risk",
"  reversal is represented only through the initial memory pool.", sep = "\n") })
  output$txt_defects <- renderText({ paste(
"Found by integrating the system, not by inspecting it:",
"",
"1. The tacrolimus drug-interaction was written as a boolean 'if letermovir",
"   present'. Solver round-off of ~1e-12 in the letermovir state latched it on",
"   permanently and trebled the tacrolimus trough (9 -> 32 ng/mL) in arms that",
"   never received letermovir. Now a concentration-dependent Emax term.",
"",
"2. The resistant-strain compartments have an UNSTABLE ZERO: whenever a mutant's",
"   exponent is positive, 1e-16 of solver round-off grows at 0.4-0.5/d and is",
"   macroscopic in ~70 days, so every arm 'developed' resistance out of floating",
"   point. Fixed with an extinction floor at one infected cell body-wide.",
"",
"3. A continuous deterministic mutation flux makes resistance CERTAIN in every",
"   arm by day 60-90, against ~5% observed on valganciclovir prophylaxis,",
"   because it allows 1e-9 of an infected cell to exist and grow. Replaced by",
"   an integral of DAYS SPENT IN THE SELECTION WINDOW.",
"",
"4. Friberg gamma = 0.17 is a transient-chemotherapy fit. Under CHRONIC dosing",
"   the steady state is ANC = CIRC0*(1-Edrug)^(1/gamma), which put a patient on",
"   mycophenolate alone at ANC 0.83 and fired the neutropenia rule in every arm.",
"",
"5. GFR0 was fixed at 55 while the initial condition was CrCl/1.1, so every run",
"   silently decayed from 86 to 55 and attributed it to CMV nephropathy.",
"",
"6. The sanctuary compartment grew without bound in EVERY arm including the",
"   untreated one, because its clearance never exceeded its growth rate.",
"",
"7. Monitoring and disease costs were never charged, so the arm that ordered",
"   twice as many PCRs appeared to be the cheaper one.",
"",
"8. The refractory-CMV switch fired at a fixed day rather than on the guideline",
"   criterion, i.e. before the resistant strain had done anything -- which made",
"   all four salvage arms numerically identical.",
sep = "\n") })
}

shinyApp(ui, server)
