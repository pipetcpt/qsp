###############################################################################
##  Kidney Transplant Rejection & Immunosuppression — QSP interactive dashboard
##  ---------------------------------------------------------------------------
##  ktx_shiny_app.R      companion to ktx_mrgsolve_model.R / ktx_qsp_model.dot
##
##  The point of this app is to let you sit where the transplant clinician sits:
##  you cannot see the alloimmune response, you can only see a trough level, a
##  creatinine, a DSA titre and a viral load — and you have exactly one dial,
##  net immunosuppression.  Turn it down and the graft rejects; turn it up and
##  the virus wins.  Tab 8 makes that trade-off explicit by sweeping the dial.
##
##  Ten tabs:
##    1  Donor & recipient      the fixed risk set, and the resulting net IS
##    2  Drug exposure & TDM    troughs vs protocol targets, MPA AUC, biologics
##    3  Cellular alloimmunity  naive/blast/effector/memory/Treg, graft infiltrate
##    4  DSA, complement, NK    the humoral arm and why rituximab under-performs
##    5  Histology (Banff)      t, i, g, ptc, cg, ci, ah, C4d over 5 years
##    6  Graft function         eGFR, creatinine, proteinuria, dd-cfDNA
##    7  Infection & safety     BK, CMV, leucocytes, glucose, cumulative IS
##    8  The U-curve            sweep net immunosuppression, see both failure modes
##    9  Scenario comparison    stack any number of saved runs side by side
##   10  Parameters & notes     the full annotated parameter set
##
##  Run:  Rscript -e 'shiny::runApp("ktx_shiny_app.R")'
###############################################################################

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## --------------------------------------------------------------------------
## Load the model definition straight out of ktx_mrgsolve_model.R so that the
## app and the batch script can never drift apart.
## --------------------------------------------------------------------------
MODEL_FILE <- "ktx_mrgsolve_model.R"
if (!file.exists(MODEL_FILE)) stop("ktx_mrgsolve_model.R must sit next to this app")
src <- readLines(MODEL_FILE, warn = FALSE)
i1  <- grep("^code <- '", src)[1]
i2  <- grep("^mod <- mcode", src)[1]
eval(parse(text = paste(src[i1:(i2 - 1)], collapse = "\n")))
mod <- mcode("ktx_shiny", code, soloc = tempdir())

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

PAL <- c("#2d5b8c", "#9c2d57", "#1f6b52", "#8a5a20", "#553f85",
         "#7d3269", "#25566f", "#6a5f1e", "#9b3535", "#46632c")

yrs <- function(d) d / 365.25

## a horizontal reference band helper
band <- function(lo, hi, fill = "#4caf50", alpha = 0.10) {
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi,
           fill = fill, alpha = alpha)
}

###############################################################################
## UI
###############################################################################
ui <- fluidPage(
  titlePanel("신장이식 거부반응 QSP 모델 · Kidney Transplant Rejection QSP Model"),
  tags$p(style = "color:#666;margin-top:-8px",
         "59-ODE mechanistic model · alloimmunity vs infection, balanced on one dial: net immunosuppression"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      tags$h4("공여자 / Donor"),
      sliderInput("KDPI", "KDPI (%)", 0, 100, 40, 5),
      sliderInput("CIT", "냉허혈시간 Cold ischaemia (h)", 1, 30, 14, 1),
      checkboxInput("DCD", "DCD 공여자 (donation after circulatory death)", FALSE),

      tags$h4("면역학적 위험 / Immunological risk"),
      sliderInput("HLAMM", "HLA A/B/DR 부적합 (0-6)", 0, 6, 3, 1),
      sliderInput("EPLET", "HLA-DR/DQ 에플렛 부적합 부하", 0, 35, 12, 1),
      numericInput("PREDSA", "사전형성 DSA (MFI)", 0, min = 0, max = 20000, step = 500),
      sliderInput("FCD28N", "CD28-null 기억 T세포 분율 (belatacept 저항성)",
                  0, 1, 0.45, 0.05),
      sliderInput("TMEM0", "이형 교차반응 기억 T세포 풀 (heterologous memory)",
                  0, 1.5, 0.70, 0.05),

      tags$h4("유지요법 / Maintenance"),
      selectInput("regimen", "기본 골격",
                  c("TAC + MPA + steroid"      = "tac",
                    "CsA + MPA + steroid"      = "csa",
                    "Belatacept + MPA + steroid (CNI-free)" = "bela",
                    "TAC(감량) + everolimus"   = "evr")),
      selectInput("induction", "유도요법", c("Basiliximab" = "bas", "rATG" = "atg", "없음" = "none")),
      sliderInput("TACTGT2", "타크로리무스 목표 트로프, 3-12개월 (ng/mL)", 3, 12, 7, 0.5),
      checkboxInput("TDM", "TDM(농도 기반 용량조절) 시행", TRUE),
      checkboxInput("CYP3A5", "CYP3A5 발현형 (*1 보유, 청소율 1.9배)", FALSE),
      sliderInput("MMFDOSE", "MMF 용량 (mg/day)", 0, 3000, 2000, 250),
      numericInput("PREDSTOP", "스테로이드 중단일 (없으면 큰 값)", 100000, min = 5, step = 1),

      tags$h4("복약 순응도 / Adherence"),
      sliderInput("ADHFINAL", "비순응 기간 중 복용 비율", 0, 1, 1, 0.05),
      numericInput("ADHSTART", "비순응 시작일", 100000, min = 0, step = 30),
      numericInput("ADHEND", "순응 회복일", 100000, min = 0, step = 30),
      sliderInput("IPVAMP", "타크로리무스 노출 변동성 IPV (진폭)", 0, 0.8, 0, 0.05),

      tags$h4("감염 / Infection"),
      checkboxInput("BKSUSC", "BK 바이러스 재활성화 소인", FALSE),
      checkboxInput("BKSCREEN", "BK 선별검사 + 선제적 면역억제 감량 (폐루프)", TRUE),
      checkboxInput("CMVDR", "CMV D+/R- 불일치", FALSE),
      numericInput("VGCSTOP", "발간시클로버 예방 종료일", 90, min = 0, step = 30),

      tags$h4("거부반응 치료 / Rejection therapy"),
      checkboxInput("MPON", "메틸프레드니솔론 충격요법", FALSE),
      numericInput("MPSTART", "  충격요법 시작일", 95, min = 0, step = 5),
      checkboxInput("ATG2ON", "ATG 구제요법 (스테로이드 저항성)", FALSE),
      numericInput("ATG2START", "  ATG 구제 시작일", 105, min = 0, step = 5),
      checkboxInput("PLEXON", "혈장교환 (PLEX)", FALSE),
      checkboxInput("IVIGON", "IVIG 2 g/kg", FALSE),
      checkboxInput("RTXON", "리툭시맙 (항CD20)", FALSE),
      checkboxInput("TCZON", "토실리주맙/클라자키주맙 (항IL-6)", FALSE),
      checkboxInput("FZBON", "펠자르타맙 (항CD38)", FALSE),
      numericInput("ABMRDAY", "ABMR 치료 시작일", 730, min = 0, step = 30),

      tags$hr(),
      actionButton("save", "이 실행을 비교탭에 저장", class = "btn-primary"),
      actionButton("clear", "비교 목록 비우기")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 공여자·수혜자", plotOutput("p_profile", height = "300px"),
                 tableOutput("t_profile"), verbatimTextOutput("txt_profile")),
        tabPanel("2 약물 노출·TDM", plotOutput("p_pk", height = "620px")),
        tabPanel("3 세포성 동종면역", plotOutput("p_cell", height = "620px")),
        tabPanel("4 DSA·보체·NK", plotOutput("p_dsa", height = "620px")),
        tabPanel("5 조직학 (Banff)", plotOutput("p_banff", height = "620px")),
        tabPanel("6 이식신 기능", plotOutput("p_func", height = "620px")),
        tabPanel("7 감염·안전성", plotOutput("p_inf", height = "620px")),
        tabPanel("8 U-곡선 (면역억제 균형)",
                 tags$p("면역억제 강도를 훑어(sweep) 거부반응 부담과 감염 부담을 동시에 계산합니다. 수 초 걸립니다."),
                 actionButton("sweep", "U-곡선 계산", class = "btn-warning"),
                 plotOutput("p_ucurve", height = "540px")),
        tabPanel("9 시나리오 비교", plotOutput("p_cmp", height = "620px"),
                 tableOutput("t_cmp")),
        tabPanel("10 파라미터·주석",
                 tags$p("모델의 전체 파라미터 (ktx_mrgsolve_model.R의 $PARAM 블록)."),
                 tableOutput("t_par"))
      )
    )
  )
)

###############################################################################
## SERVER
###############################################################################
server <- function(input, output, session) {

  saved <- reactiveVal(list())

  ## a tiny stacker so the app has no gridExtra dependency
  gridExtra_like <- function(plots, heights) {
    n <- length(plots)
    grid::grid.newpage()
    tot <- sum(heights)
    y0 <- 1
    for (i in seq_len(n)) {
      h <- heights[i] / tot
      print(plots[[i]], vp = grid::viewport(x = 0.5, y = y0 - h / 2,
                                            width = 1, height = h))
      y0 <- y0 - h
    }
  }


  pars <- reactive({
    p <- list(
      KDPI = input$KDPI, CIT = input$CIT, DCD = as.numeric(input$DCD),
      HLAMM = input$HLAMM, EPLET = input$EPLET, PREDSA = input$PREDSA,
      FCD28N = input$FCD28N, TMEM0 = input$TMEM0,
      TACTGT2 = input$TACTGT2, TACTGT3 = max(4, input$TACTGT2 - 1),
      TDM = as.numeric(input$TDM), CYP3A5 = as.numeric(input$CYP3A5),
      MMFDOSE = input$MMFDOSE, MMFON = as.numeric(input$MMFDOSE > 0),
      PREDSTOP = input$PREDSTOP,
      ADHFINAL = input$ADHFINAL, ADHSTART = input$ADHSTART, ADHEND = input$ADHEND,
      IPVAMP = input$IPVAMP, PIPV = 21,
      BKSUSC = as.numeric(input$BKSUSC), BKSCREEN = as.numeric(input$BKSCREEN),
      CMVDR = as.numeric(input$CMVDR), CMVSUSC = as.numeric(input$CMVDR),
      VGCSTOP = input$VGCSTOP,
      MPON = as.numeric(input$MPON), MPSTART = input$MPSTART,
      ATG2ON = as.numeric(input$ATG2ON), ATG2START = input$ATG2START,
      PLEXON = as.numeric(input$PLEXON), PLEXSTART = input$ABMRDAY,
      IVIGON = as.numeric(input$IVIGON), IVIGSTART = input$ABMRDAY + 10,
      RTXON  = as.numeric(input$RTXON),  RTXDAY   = input$ABMRDAY + 12,
      TCZON  = as.numeric(input$TCZON),  TCZSTART = input$ABMRDAY, TCZN = 18,
      FZBON  = as.numeric(input$FZBON),  FZBSTART = input$ABMRDAY
    )
    # maintenance backbone
    p$TACON <- 0; p$CSAON <- 0; p$BELAON <- 0; p$EVRON <- 0
    if (input$regimen == "tac")  p$TACON <- 1
    if (input$regimen == "csa")  p$CSAON <- 1
    if (input$regimen == "bela") { p$BELAON <- 1; p$TDM <- 0 }
    if (input$regimen == "evr") {
      p$TACON <- 1; p$EVRON <- 1; p$MMFON <- 0
      p$TACTGT1 <- 5; p$TACTGT2 <- 4; p$TACTGT3 <- 4
    }
    # induction
    p$BASON <- as.numeric(input$induction == "bas")
    p$ATGON <- as.numeric(input$induction == "atg")
    p
  })

  sim <- reactive({
    mod %>% param(pars()) %>%
      mrgsim(end = 1825, delta = 2, hmax = 3, atol = 1e-8, rtol = 1e-5) %>%
      as_tibble()
  })

  long <- function(d, vars, labs) {
    d %>% select(time, all_of(vars)) %>%
      rename_with(~labs, all_of(vars)) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(name, levels = labs))
  }

  gg <- function(d, ncol = 2, ylab = NULL) {
    ggplot(d, aes(yrs(time), value)) +
      geom_line(linewidth = 0.7, colour = PAL[1]) +
      facet_wrap(~name, scales = "free_y", ncol = ncol) +
      labs(x = "이식 후 경과 (년)", y = ylab) + THEME
  }

  ## ---- 1  profile ---------------------------------------------------------
  output$p_profile <- renderPlot({
    d <- sim()
    long(d, c("NISr", "CNINHr", "MPAINHr", "STINHr", "COSBLKr", "ATGDr"),
         c("Net immunosuppression (NIS)", "Calcineurin inhibition",
           "IMPDH inhibition", "Glucocorticoid effect",
           "Costimulation blockade", "Lymphodepletion")) %>%
      gg(ncol = 3, ylab = "effect (0-1)")
  })

  output$t_profile <- renderTable({
    d <- sim()
    idx <- function(day) which.min(abs(d$time - day))
    data.frame(
      시점 = c("3개월", "12개월", "3년", "5년"),
      eGFR = round(d$EGFRr[c(idx(90), idx(365), idx(1095), idx(1825))], 1),
      `혈청 Cr (mg/dL)` = round(d$SCRDL[c(idx(90), idx(365), idx(1095), idx(1825))], 2),
      `TAC C0` = round(d$TACC0r[c(idx(90), idx(365), idx(1095), idx(1825))], 1),
      `DSA (MFI)` = round(d$DSA[c(idx(90), idx(365), idx(1095), idx(1825))], 0),
      `Banff t` = round(d$BANFF_T[c(idx(90), idx(365), idx(1095), idx(1825))], 2),
      `g+ptc` = round(d$MVISUM[c(idx(90), idx(365), idx(1095), idx(1825))], 2),
      `이식신 생존확률` = round(d$GS[c(idx(90), idx(365), idx(1095), idx(1825))], 3),
      check.names = FALSE
    )
  })

  output$txt_profile <- renderPrint({
    d <- sim()
    cat("지연이식신기능(DGF) 예상:", ifelse(d$DGF[1] > 0.5, "예", "아니오"), "\n")
    cat("최저 eGFR:", round(min(d$EGFRr), 1), "mL/min/1.73 m2  (day",
        round(d$time[which.min(d$EGFRr)]), ")\n")
    cat("최고 Banff t:", round(max(d$BANFF_T), 2),
        " / 최고 g+ptc:", round(max(d$MVISUM), 2), "\n")
    cat("최고 DSA:", round(max(d$DSA)), "MFI\n")
    cat("최고 BK 바이러스 부하:", round(max(d$BKLOG), 2), "log10 copies/mL\n")
    cat("5년 eGFR 기울기:",
        round((d$EGFRr[which.min(abs(d$time - 1825))] -
               d$EGFRr[which.min(abs(d$time - 365))]) / 4, 2),
        "mL/min/1.73 m2/년\n")
  })

  ## ---- 2  PK --------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim()
    p1 <- ggplot(d, aes(yrs(time), TACC0r)) + band(5, 10) +
      geom_line(colour = PAL[1], linewidth = 0.7) +
      labs(x = NULL, y = "타크로리무스 트로프 (ng/mL)",
           title = "TAC C0 (녹색 = 5-10 ng/mL 목표역)") + THEME
    p2 <- ggplot(d, aes(yrs(time), MPAAUC)) + band(30, 60) +
      geom_line(colour = PAL[3], linewidth = 0.7) +
      labs(x = NULL, y = "MPA AUC0-12 (mg*h/L)", title = "MPA 노출 (녹색 = 30-60)") + THEME
    p3 <- long(d, c("CPREDr", "EVRC0", "CSAC0r"),
               c("Prednisolone-eq (ng/mL)", "Everolimus C0 (ng/mL)",
                 "Cyclosporine C0 (ng/mL)")) %>% gg(ncol = 3)
    p4 <- long(d, c("CATGr", "CBELAr", "CRTXr", "CTCZr", "CFZBr"),
               c("rATG (ug/mL)", "Belatacept (ug/mL)", "Rituximab (ug/mL)",
                 "anti-IL-6 (ug/mL)", "anti-CD38 (ug/mL)")) %>% gg(ncol = 5)
    gridExtra_like(list(p1, p2, p3, p4), c(1, 1, 1, 1))
  })

  ## ---- 3  cellular --------------------------------------------------------
  output$p_cell <- renderPlot({
    d <- sim()
    long(d, c("TN", "TACT", "TEFF", "TMEM", "TREG", "TINF", "IL2", "IFNG", "AG"),
         c("Naive alloreactive T (TN)", "Activated blasts (TACT)",
           "Effector T (TEFF)", "Memory T (TMEM)", "Treg", "Graft-infiltrating T (TINF)",
           "IL-2", "IFN-gamma", "Antigen presentation (AG)")) %>%
      gg(ncol = 3, ylab = "relative units")
  })

  ## ---- 4  humoral ---------------------------------------------------------
  output$p_dsa <- renderPlot({
    d <- sim()
    p1 <- ggplot(d, aes(yrs(time), DSA)) +
      geom_hline(yintercept = 1000, linetype = 2, colour = "#9c2d57") +
      geom_line(colour = PAL[2], linewidth = 0.8) +
      labs(x = NULL, y = "DSA (MFI)",
           title = "공여자특이항체 (점선 = 통상적 검출 임계 1000 MFI)") + THEME
    p2 <- long(d, c("BN", "BGC", "PB", "LLPC", "C4D", "NKACT", "ENDO", "MVI"),
               c("Naive/memory B (BN)", "Germinal-centre B (BGC)",
                 "Plasmablasts (PB)", "Long-lived plasma cells (LLPC)",
                 "C4d deposition", "Activated NK", "Endothelial activation",
                 "Microvascular inflammation")) %>% gg(ncol = 4)
    gridExtra_like(list(p1, p2), c(1, 2))
  })

  ## ---- 5  Banff -----------------------------------------------------------
  output$p_banff <- renderPlot({
    d <- sim()
    dd <- long(d, c("BANFF_T", "BANFF_I", "BANFF_G", "BANFF_PTC",
                    "BANFF_CG", "BANFF_CI", "BANFF_AH", "BANFF_C4D"),
               c("t (tubulitis)", "i (interstitial)", "g (glomerulitis)",
                 "ptc (capillaritis)", "cg (transplant glomerulopathy)",
                 "ci (fibrosis)", "ah (arteriolar hyalinosis)", "C4d"))
    ggplot(dd, aes(yrs(time), value)) +
      geom_hline(yintercept = 1, linetype = 3, colour = "grey60") +
      geom_line(linewidth = 0.7, colour = PAL[4]) +
      facet_wrap(~name, ncol = 4) +
      coord_cartesian(ylim = c(0, 3)) +
      labs(x = "이식 후 경과 (년)", y = "Banff 점수 (0-3)",
           title = "Banff 병변 점수 추이") + THEME
  })

  ## ---- 6  function --------------------------------------------------------
  output$p_func <- renderPlot({
    d <- sim()
    p1 <- ggplot(d, aes(yrs(time), EGFRr)) +
      geom_hline(yintercept = 15, linetype = 2, colour = "#9b3535") +
      geom_line(colour = PAL[1], linewidth = 0.8) +
      labs(x = NULL, y = "eGFR (mL/min/1.73 m2)",
           title = "이식신 기능 (점선 = 투석 재개 임계 15)") + THEME
    p2 <- long(d, c("SCRDL", "PROT", "CFDNA", "GS"),
               c("혈청 크레아티닌 (mg/dL)", "단백뇨 UPCR (g/g)",
                 "dd-cfDNA (%)", "사망 검열 이식신 생존확률")) %>% gg(ncol = 4)
    gridExtra_like(list(p1, p2), c(1, 1))
  })

  ## ---- 7  infection -------------------------------------------------------
  output$p_inf <- renderPlot({
    d <- sim()
    p1 <- ggplot(d, aes(yrs(time), BKLOG)) +
      geom_hline(yintercept = 4, linetype = 2, colour = "#9b3535") +
      geom_line(colour = PAL[3], linewidth = 0.8) +
      labs(x = NULL, y = "BK (log10 copies/mL)",
           title = "BK 폴리오마바이러스 (점선 = 10^4, 선제적 감량 임계)") + THEME
    p2 <- long(d, c("CMVLOG", "WBC", "GLU", "NISAUC", "BKVAN", "NISr"),
               c("CMV (log10 IU/mL)", "백혈구 (10^9/L)", "공복혈당 지표 (mmol/L)",
                 "누적 면역억제 노출 (NIS-days)", "BK 신병증 부담",
                 "순 면역억제 (NIS)")) %>% gg(ncol = 3)
    gridExtra_like(list(p1, p2), c(1, 1.4))
  })

  ## ---- 8  the U-curve -----------------------------------------------------
  ucurve <- eventReactive(input$sweep, {
    base <- pars()
    base$BKSUSC <- 1        # a patient who CAN reactivate BK - otherwise there
    base$BKSCREEN <- 0      # is no other arm to the trade-off
    tg <- seq(2, 14, by = 1)
    withProgress(message = "면역억제 강도 훑는 중...", value = 0, {
      do.call(rbind, lapply(seq_along(tg), function(i) {
        incProgress(1 / length(tg))
        p <- base; p$TACTGT1 <- tg[i] + 2; p$TACTGT2 <- tg[i]; p$TACTGT3 <- tg[i]
        s <- mod %>% param(p) %>%
          mrgsim(end = 1825, delta = 5, hmax = 3, atol = 1e-8, rtol = 1e-5) %>%
          as_tibble()
        data.frame(
          target   = tg[i],
          NIS      = mean(s$NISr),
          rejection = max(s$BANFF_T) + max(s$MVISUM) / 2,
          infection = max(s$BKLOG),
          eGFR5     = s$EGFRr[which.min(abs(s$time - 1825))],
          survival  = s$GS[which.min(abs(s$time - 1825))]
        )
      }))
    })
  })

  output$p_ucurve <- renderPlot({
    u <- ucurve()
    dd <- u %>%
      select(NIS, `거부반응 부담 (Banff t + g+ptc/2)` = rejection,
             `감염 부담 (peak BK log10)` = infection,
             `5년 eGFR` = eGFR5, `5년 이식신 생존확률` = survival) %>%
      pivot_longer(-NIS)
    ggplot(dd, aes(NIS, value)) +
      geom_line(linewidth = 0.9, colour = PAL[2]) +
      geom_point(size = 2, colour = PAL[2]) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "평균 순 면역억제 (NIS)", y = NULL,
           title = "면역억제 강도의 U자형 상충관계",
           subtitle = "왼쪽으로 가면 거부반응, 오른쪽으로 가면 바이러스 — 최적점은 그 사이에 있다") +
      THEME
  })

  ## ---- 9  comparison ------------------------------------------------------
  observeEvent(input$save, {
    s <- sim()
    s$run <- paste0("run ", length(saved()) + 1, ": ",
                    input$regimen, "/", input$induction,
                    " TGT", input$TACTGT2,
                    if (input$ADHFINAL < 1) paste0(" ADH", input$ADHFINAL) else "",
                    if (input$BKSUSC) " BK" else "")
    saved(c(saved(), list(s)))
  })
  observeEvent(input$clear, saved(list()))

  output$p_cmp <- renderPlot({
    if (!length(saved())) return(NULL)
    d <- bind_rows(saved())
    d %>%
      select(time, run, eGFR = EGFRr, DSA, `Banff t` = BANFF_T,
             `g+ptc` = MVISUM, `BK log10` = BKLOG, NIS = NISr) %>%
      pivot_longer(-c(time, run)) %>%
      ggplot(aes(yrs(time), value, colour = run)) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = rep(PAL, 3)) +
      labs(x = "이식 후 경과 (년)", y = NULL) + THEME +
      guides(colour = guide_legend(ncol = 2))
  })

  output$t_cmp <- renderTable({
    if (!length(saved())) return(NULL)
    do.call(rbind, lapply(saved(), function(s) {
      i <- function(day) which.min(abs(s$time - day))
      data.frame(
        run = s$run[1],
        `eGFR 1년` = round(s$EGFRr[i(365)], 1),
        `eGFR 5년` = round(s$EGFRr[i(1825)], 1),
        `기울기/년` = round((s$EGFRr[i(1825)] - s$EGFRr[i(365)]) / 4, 2),
        `최고 Banff t` = round(max(s$BANFF_T), 2),
        `DSA 5년` = round(s$DSA[i(1825)]),
        `cg 5년` = round(s$BANFF_CG[i(1825)], 2),
        `BK 최고` = round(max(s$BKLOG), 2),
        `생존확률 5년` = round(s$GS[i(1825)], 3),
        check.names = FALSE
      )
    }))
  })

  ## ---- 10  parameters -----------------------------------------------------
  output$t_par <- renderTable({
    b <- src[(grep("^\\$PARAM", src)[1] + 1):(grep("^\\$CMT", src)[1] - 1)]
    b <- b[grepl(":", b)]
    parts <- strsplit(b, "\\s*:\\s*")
    data.frame(
      parameter   = vapply(parts, `[`, "", 1),
      value       = vapply(parts, `[`, "", 2),
      description = vapply(parts, function(x) paste(x[-c(1, 2)], collapse = " : "), "")
    )
  })
}

shinyApp(ui, server)
