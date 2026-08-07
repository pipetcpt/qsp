## =====================================================================
##  ods_shiny_app.R
##  Osmotic Demyelination Syndrome (ODS) — interactive QSP dashboard
##  삼투성 탈수초 증후군 — 인터랙티브 대시보드 (11 tabs)
##
##  The app is organised around one question that the model answers and a
##  bedside chart cannot:  WHO IS SETTING THE CORRECTION RATE?
##
##  Tab 1  환자 프로파일          Patient & phenotype
##  Tab 2  전해질·체액 PK         Sodium, water, urine
##  Tab 3  신장 — 속도를 정하는 곳 The kidney: AVP, U_osm, EFWC
##  Tab 4  뇌 삼투질 (PD 핵심)     Brain osmolytes and OMEGA
##  Tab 5  손상 캐스케이드         Astrocyte -> microglia -> myelin
##  Tab 6  임상 종말점            Deficit, MRI, the biphasic course
##  Tab 7  중심 실험 (AVP 스위치)  The counterfactual
##  Tab 8  교정 속도 지도         The FOSM x rate safety map
##  Tab 9  재하강 마감시한         The relowering deadline
##  Tab 10 시나리오 비교          Scenario comparison
##  Tab 11 바이오마커 & 모델 설명  Biomarkers and what is fitted
##
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  Run:  shiny::runApp("ods_shiny_app.R")
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source_model <- function() {
  e <- new.env()
  sys.source("ods_mrgsolve_model.R", envir = e)
  e
}
E   <- source_model()
mod <- E$mod

THEME <- theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey93"),
        legend.position = "bottom", panel.grid.minor = element_blank())

## ---------------------------------------------------------------------
##  Phenotype presets.  Each one is a different KIDNEY at the same sodium.
## ---------------------------------------------------------------------
PHENO <- list(
  "SIADH (정상 혈량)"          = list(kind = "siadh",     days = 21,
    note = "AVP가 자율적으로 분비된다. 요삼투압이 높게 고정되어 있고, 요 (Na+K)/혈청 Na 비가 1을 넘으므로 수분제한만으로는 교정되지 않는다."),
  "저혈량성 (구토·이뇨제)"      = list(kind = "hypovol",   days = 7,
    note = "★ ODS가 가장 잘 생기는 표현형. AVP가 용적 자극으로 올라가 있고, 생리식염수로 용적을 채우면 그 자극이 사라지면서 신장이 스스로 나트륨을 올린다."),
  "티아지드 유발"              = list(kind = "thiazide",  days = 21,
    note = "희석 분절(NCC)이 차단되어 최소 요삼투압이 올라가 있다. 약을 끊는 순간 희석능이 돌아온다."),
  "맥주 다음증 / tea-and-toast" = list(kind = "potomania", days = 21,
    note = "용질 섭취가 <250 mOsm/d라 자유수분 배설 능력 자체가 낮다. 식사를 복구하면 그것만으로 과교정이 일어난다."),
  "부신부전"                   = list(kind = "adrenal",   days = 21,
    note = "코르티솔 결핍으로 AVP가 억제되지 않는다. 하이드로코르티손을 주는 순간 수분이뇨가 열린다."),
  "급성 (<12 h, 운동·MDMA)"     = list(kind = "acute",     days = 8/24,
    note = "★ 반대 방향의 병. 유기 삼투질이 아직 빠져나가지 않아 뇌가 부어 있다. 위험은 탈출(herniation)이지 탈수초가 아니다.")
)

build_patient <- function(kind, days, target = 110) {
  if (kind == "acute") {
    w <- uniroot(function(x) {
      o <- mod %>% param(AVPSIADH = 6, ACUTE = 1, WLOAD = x, TWLEND = days) %>%
             mrgsim(end = days, delta = 1/96)
      tail(o$SODIUM, 1) - target
    }, c(1, 60), tol = 1e-4)$root
    p <- list(AVPSIADH = 6, ACUTE = 0)
    o <- mod %>% param(AVPSIADH = 6, ACUTE = 1, WLOAD = w, TWLEND = days) %>%
           mrgsim(end = days, delta = 1/96)
    return(list(param = p, state = as.numeric(tail(o, 1)[, names(init(mod))]), build = o))
  }
  w <- E$solve_win(mod, kind, target, days)
  r <- E$make_chronic(mod, kind, win = w, days = days)
  list(param = r$param, state = r$state, build = r$out)
}

start_at <- function(m, st)
  do.call(init, c(list(m), as.list(setNames(st, names(init(mod))))))

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel("삼투성 탈수초 증후군 (ODS) — QSP 대시보드"),
  tags$p(style = "color:#555;margin-top:-8px",
    HTML("<b>Ω = ORG_set(Osm_eff) − ORG</b> — 유기 삼투질 결핍. ",
         "정상 뇌에서 0, <i>적응된</i> 저나트륨혈증 뇌에서도 0, ",
         "혈장 긴장도가 전사보다 빨리 움직일 때만 양이 된다.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("pheno", "표현형 (phenotype)", names(PHENO)),
      sliderInput("na0", "제시 시점 혈청 [Na] (mmol/L)", 100, 128, 110, step = 1),
      hr(),
      h5("교정 처방 (prescription)"),
      radioButtons("mode", NULL,
        c("속도 지정 3% NaCl 적정" = "ctrl",
          "0.9% 생리식염수만"       = "ns",
          "수분 제한만"             = "fr",
          "톨밥탄 15 mg/일"         = "tlv",
          "경구 요소 30 g/일"       = "urea"), selected = "ctrl"),
      conditionalPanel("input.mode == 'ctrl'",
        sliderInput("rate", "목표 교정 속도 (mmol/L/24 h)", 2, 24, 6, step = 1),
        sliderInput("cap",  "교정 종료 [Na]", 120, 145, 135, step = 1)),
      conditionalPanel("input.mode == 'ns'",
        sliderInput("r09", "0.9% NaCl (L/일)", 0, 4, 2, step = 0.5)),
      conditionalPanel("input.mode == 'fr'",
        sliderInput("win", "수분 섭취 (L/일)", 0.3, 2.5, 1.0, step = 0.1)),
      hr(),
      h5("속도 제어 (rate control)"),
      checkboxInput("ddavp", "선제적 데스모프레신 클램프 2 µg q8h", FALSE),
      checkboxInput("rescue", "과교정 시 재하강 (D5W + DDAVP)", FALSE),
      conditionalPanel("input.rescue",
        sliderInput("trescue", "재하강 시작 (시간)", 6, 96, 24, step = 6),
        sliderInput("nares",   "재하강 목표 [Na]",   110, 128, 118, step = 1)),
      hr(),
      h5("위험 수식인자 (risk modifiers)"),
      sliderInput("fosm", "유기 삼투질 수송 용량 FOSM", 0.35, 1.0, 1.0, step = 0.05),
      helpText(HTML("<small>1.00 정상 · 0.55 알코올중독/영양불량 · 0.70 간질환</small>")),
      sliderInput("fnut", "별아교세포 에너지 예비력 FNUT", 0.6, 1.0, 1.0, step = 0.05),
      sliderInput("kdef", "총 체내 칼륨 결핍 (mmol)", 0, 700, 0, step = 50),
      sliderInput("rkcl", "KCl 보충 (mmol/일)", 0, 200, 0, step = 20),
      hr(),
      h5("보조 치료 (adjuncts, 전임상 근거)"),
      checkboxInput("dex",  "덱사메타손 16 mg/일", FALSE),
      checkboxInput("mino", "미노사이클린 200 mg/일", FALSE),
      hr(),
      sliderInput("tend", "시뮬레이션 기간 (일)", 7, 90, 45, step = 1),
      actionButton("go", "실행 (Run)", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("① 환자 프로파일",
          br(), htmlOutput("phenoNote"),
          h4("모델이 스스로 만들어 낸 적응된 뇌"),
          p(HTML("아래 값은 <b>입력이 아니라 출력</b>이다. 저나트륨혈증을 시뮬레이션으로 발생시키면 뇌가 스스로 삼투질을 버린다.")),
          DTOutput("tblPatient"), br(),
          plotOutput("plotBuild", height = "340px")),

        tabPanel("② 전해질·체액",
          br(), plotOutput("plotNa", height = "300px"),
          plotOutput("plotFluid", height = "300px"),
          h5("24시간 최대 상승폭"), verbatimTextOutput("txtRate")),

        tabPanel("③ 신장 — 속도를 정하는 곳",
          br(),
          p(HTML("<b>EFWC = V<sub>urine</sub> × (1 − (U<sub>Na</sub>+U<sub>K</sub>)/S<sub>Na</sub>)</b>. ",
                 "처방이 아니라 이 값이 교정 속도를 정한다.")),
          plotOutput("plotKidney", height = "420px"),
          h5("Furst 비 (U_Na+U_K)/S_Na — 1을 넘으면 수분제한은 실패한다"),
          plotOutput("plotFurst", height = "220px")),

        tabPanel("④ 뇌 삼투질 (Ω)",
          br(),
          p(HTML("채널로 나가고 <b>전사를 거쳐</b> 들어온다. 그 비대칭이 병이다.")),
          plotOutput("plotOsmo", height = "360px"),
          plotOutput("plotOmega", height = "300px")),

        tabPanel("⑤ 손상 캐스케이드",
          br(),
          p(HTML("순서가 중요하다: <b>별아교세포 → 미세아교세포·BBB → 희소돌기아교세포 → 수초</b>.")),
          plotOutput("plotInjury", height = "480px")),

        tabPanel("⑥ 임상 종말점",
          br(),
          p(HTML("<b>이상성 경과</b>: 나트륨이 정상으로 돌아오고 환자가 좋아진 뒤에 악화가 온다. ",
                 "그리고 MRI는 그보다 더 늦다.")),
          plotOutput("plotClin", height = "380px"),
          DTOutput("tblTiming")),

        tabPanel("⑦ 중심 실험 (AVP 스위치)",
          br(),
          p(HTML("같은 환자, 같은 수액, 고장성 식염수는 <b>양쪽 모두 0 mL</b>. ",
                 "한쪽만 AVP가 생리적으로 반응한다.")),
          actionButton("goCF", "반사실 실험 실행", class = "btn-warning"),
          br(), br(), plotOutput("plotCF", height = "480px"),
          verbatimTextOutput("txtCF")),

        tabPanel("⑧ 교정 속도 지도",
          br(),
          p(HTML("역치 Ω* = 8 mOsm/kg는 <b>모든 환자에게 같은 하나의 숫자</b>다. ",
                 "위험인자는 오직 FOSM(운반체 용량)에만 작용한다.")),
          actionButton("goMap", "안전 지도 계산", class = "btn-warning"),
          br(), br(), plotOutput("plotMap", height = "440px"),
          DTOutput("tblMap")),

        tabPanel("⑨ 재하강 마감시한",
          br(),
          p(HTML("과교정이 일어난 뒤 재하강(D5W + DDAVP + 수분제한)을 언제까지 하면 되는가?")),
          actionButton("goResc", "마감시한 계산", class = "btn-warning"),
          br(), br(), plotOutput("plotResc", height = "400px"),
          DTOutput("tblResc")),

        tabPanel("⑩ 시나리오 비교",
          br(),
          checkboxGroupInput("cmp", "비교할 처방", inline = TRUE,
            choices = c("3% NaCl +6/일", "3% NaCl +12/일", "3% NaCl +20/일",
                        "0.9% NaCl 2 L/일", "수분제한 1 L/일",
                        "톨밥탄 15 mg", "DDAVP 클램프 + 3% NaCl"),
            selected = c("3% NaCl +6/일", "3% NaCl +12/일", "0.9% NaCl 2 L/일")),
          actionButton("goCmp", "비교 실행", class = "btn-warning"),
          br(), br(), plotOutput("plotCmp", height = "460px"),
          DTOutput("tblCmp")),

        tabPanel("⑪ 바이오마커 & 모델 설명",
          br(),
          h4("측정 가능한 것과 측정 불가능한 것"),
          plotOutput("plotBio", height = "340px"),
          hr(),
          h4("무엇을 맞췄고 무엇을 예측했는가"),
          htmlOutput("txtFit"))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  patient <- eventReactive(input$go, {
    ph <- PHENO[[input$pheno]]
    withProgress(message = "적응된 뇌를 시뮬레이션으로 만드는 중...", {
      build_patient(ph$kind, ph$days, input$na0)
    })
  }, ignoreNULL = FALSE)

  rx_param <- reactive({
    p <- list(FOSM = input$fosm, FNUT = input$fnut, RKCL = input$rkcl,
              DEXON = as.numeric(input$dex), MINOON = as.numeric(input$mino))
    switch(input$mode,
      ctrl = c(p, list(CTRLON = 1, RATETGT = input$rate, NASTART = input$na0,
                       NACAP = input$cap, WIN = 1.0)),
      ns   = c(p, list(R09 = input$r09)),
      fr   = c(p, list(WIN = input$win)),
      tlv  = p,
      urea = c(p, list(UREADOSE = 30, WIN = 1.0)))
  })

  sim <- eventReactive(input$go, {
    pt <- patient()
    st <- pt$state
    st[which(names(init(mod)) == "KE")] <- st[which(names(init(mod)) == "KE")] - input$kdef
    pp <- c(pt$param, rx_param())
    if (PHENO[[input$pheno]]$kind == "hypovol")
      pp <- c(pp, list(NALOSS = 0, KLOSS = 0, WLOSS = 0, TLOSSEND = 0))
    ev <- NULL
    if (input$ddavp) ev <- E$ddavp_q8(0, 6)
    if (input$mode == "tlv") ev <- c(ev, E$tolvaptan(7))
    m <- start_at(mod, st) %>% param(pp)
    if (!input$rescue) return(as_tibble(m %>% mrgsim(events = ev, end = input$tend, delta = 0.25)))
    ## the rescue is a change of PRESCRIPTION, so it is run as a second phase
    th <- input$trescue / 24
    a  <- m %>% mrgsim(events = ev, end = th, delta = 0.05)
    st2 <- as.numeric(tail(a, 1)[, names(init(mod))])
    pr <- c(pp, list(RESCUE = 1, TRESCUE = th, NARES = input$nares,
                     DURRES = 1.0, WIN = 0.5, R09 = 0))
    b <- start_at(mod, st2) %>% param(pr) %>%
      mrgsim(events = E$ddavp_q8(th, th + 2.5), end = input$tend, delta = 0.25)
    bind_rows(as_tibble(a), as_tibble(b) %>% filter(time > 0) %>% mutate(time = time))
  }, ignoreNULL = FALSE)

  gg <- function(d, vars, labs_) {
    d %>% select(time, all_of(vars)) %>% pivot_longer(-time) %>%
      mutate(name = factor(name, levels = vars, labels = labs_)) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.85, colour = "#2b6cb0") +
      facet_wrap(~name, scales = "free_y") + labs(x = "일 (days)", y = NULL) + THEME
  }

  ## --- Tab 1 -----------------------------------------------------------
  output$phenoNote <- renderUI({
    HTML(paste0("<div style='background:#f0f6ff;padding:10px;border-left:4px solid #2b6cb0'>",
                PHENO[[input$pheno]]$note, "</div>"))
  })
  output$tblPatient <- renderDT({
    b <- as_tibble(patient()$build); e <- tail(b, 1)
    datatable(tibble(
      항목 = c("혈청 [Na] (mmol/L)", "유효 삼투압 (mOsm/kg)", "뇌 수분 (mL/100 g, 정상 80.00)",
               "총 유기 삼투질 (mOsm/kg, 정상 48.0)", "myo-inositol (정상 7.0)",
               "무기 이온 (정상 227)", "Ω (유기 삼투질 결핍)", "혈장 AVP (pg/mL)",
               "요삼투압 (mOsm/kg)", "요량 (L/일)", "Furst 비 (U_Na+U_K)/S_Na"),
      값 = round(c(e$SODIUM, e$TONICITY, e$BRAINH2O, e$ORGTOT, e$MYOINOS,
                   e$ORGTOT * 0 + (e$TONICITY - 10 - e$ORGTOT), e$OMEGAc,
                   e$AVPpg, e$UOSMOL, e$UVOL, e$FURST), 2)),
      rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })
  output$plotBuild <- renderPlot({
    gg(as_tibble(patient()$build),
       c("SODIUM","BRAINH2O","ORGTOT","MYOINOS"),
       c("혈청 [Na]","뇌 수분 (mL/100 g)","총 유기 삼투질","myo-inositol")) +
      ggtitle("저나트륨혈증이 만들어지는 동안 뇌가 삼투질을 버린다")
  })

  ## --- Tab 2 -----------------------------------------------------------
  output$plotNa <- renderPlot(
    gg(sim(), c("SODIUM","TONICITY","POTASSIUM","BUNmgdL"),
       c("혈청 [Na] (mmol/L)","유효 삼투압","혈청 [K]","BUN (mg/dL)")))
  output$plotFluid <- renderPlot(
    gg(sim(), c("UVOL","UOSMOL","EFWCL","BRAINH2O"),
       c("요량 (L/일)","요삼투압","자유수분 청소율 EFWC (L/일)","뇌 수분")))
  output$txtRate <- renderText({
    d <- sim(); r <- E$na_rise_24(d$time, d$SODIUM)
    lim <- if (input$fosm < 0.8) 8 else 12
    sprintf("어떤 24시간 창에서든 최대 상승폭 = %.1f mmol/L   (이 환자의 모델 한계 ≈ %d)\n%s",
            r, lim, if (r > lim) "*** 한계 초과 ***" else "한계 이내")
  })

  ## --- Tab 3 -----------------------------------------------------------
  output$plotKidney <- renderPlot(
    gg(sim(), c("AVPpg","UOSMOL","UVOL","EFWCL","UNAKc","SODIUM"),
       c("혈장 AVP (pg/mL)","요삼투압 (mOsm/kg)","요량 (L/일)",
         "EFWC (L/일)","요 (Na+K) (mmol/L)","혈청 [Na]")))
  output$plotFurst <- renderPlot({
    ggplot(sim(), aes(time, FURST)) + geom_line(linewidth = 0.9, colour = "#b0416a") +
      geom_hline(yintercept = 1, linetype = 2) +
      labs(x = "일", y = "(U_Na+U_K)/S_Na") + THEME
  })

  ## --- Tab 4 -----------------------------------------------------------
  output$plotOsmo <- renderPlot({
    d <- as_tibble(sim())
    d %>% select(time, ORGTOT, ORGSETC, MYOINOS) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name, ORGTOT = "실제 유기 삼투질 (ORG)",
                           ORGSETC = "현재 긴장도가 요구하는 값 (ORG_set)",
                           MYOINOS = "myo-inositol")) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.95) +
      labs(x = "일", y = "mOsm/kg 뇌수분", colour = NULL,
           title = "따라잡지 못하는 것이 병이다") + THEME
  })
  output$plotOmega <- renderPlot({
    d <- as_tibble(sim())
    ggplot(d, aes(time)) +
      geom_ribbon(aes(ymin = THRESH, ymax = pmax(STRESSc, THRESH)),
                  fill = "#f4b8b8", alpha = 0.7) +
      geom_line(aes(y = STRESSc), linewidth = 1, colour = "#b0416a") +
      geom_line(aes(y = OMEGAc), linewidth = 0.8, colour = "#2f7a53") +
      geom_hline(aes(yintercept = THRESH), linetype = 2) +
      labs(x = "일", y = "mOsm/kg",
           title = "초록 = Ω(유기 삼투질 결핍) · 빨강 = 총 삼투 스트레스 · 점선 = 역치 8.0") + THEME
  })

  ## --- Tab 5 -----------------------------------------------------------
  output$plotInjury <- renderPlot(
    gg(sim(), c("STRESSc","INJRATE","BRAINH2O","SHRINKc","MRIPOS","DEFICIT"),
       c("삼투 스트레스","손상 속도","뇌 수분","뇌 위축 (%)","MRI 신호","임상 결손")))

  ## --- Tab 6 -----------------------------------------------------------
  output$plotClin <- renderPlot({
    d <- as_tibble(sim())
    d %>% select(time, SODIUM, DEFICIT, MRIPOS) %>% pivot_longer(-time) %>%
      mutate(name = recode(name, SODIUM = "혈청 [Na]",
                           DEFICIT = "임상 결손 (0-100)", MRIPOS = "MRI 신호")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.95, colour = "#2b6cb0") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "일", y = NULL,
           title = "나트륨이 먼저 정상이 되고, 환자는 그 다음에 나빠진다") + THEME
  })
  output$tblTiming <- renderDT({
    d <- as_tibble(sim())
    f <- function(v, th) { i <- which(v > th)[1]; if (is.na(i)) NA else d$time[i] }
    datatable(tibble(
      사건 = c("혈청 [Na] 정상화 (>135)", "증상 발현 (결손 >10)",
               "MRI 양성 (신호 >0.25)", "최대 결손"),
      `발생일` = c(f(d$SODIUM, 135), f(d$DEFICIT, 10), f(d$MRIPOS, 0.25),
                   d$time[which.max(d$DEFICIT)])),
      rownames = FALSE, options = list(dom = "t"))
  })

  ## --- Tab 7 : the counterfactual ---------------------------------------
  cf <- eventReactive(input$goCF, {
    pt <- build_patient("hypovol", 7, input$na0)
    base <- c(pt$param, list(NALOSS = 0, KLOSS = 0, WLOSS = 0, TLOSSEND = 0, R09 = 2))
    iavp <- which(names(init(mod)) == "AVP")
    A <- start_at(mod, pt$state) %>% param(base) %>% mrgsim(end = 30, delta = 0.1)
    B <- start_at(mod, pt$state) %>%
      param(c(base, list(AVPFREEZE = 1, TFREEZE = 0, AVPFRZV = pt$state[iavp]))) %>%
      mrgsim(end = 30, delta = 0.1)
    bind_rows(as_tibble(A) %>% mutate(arm = "AVP가 반응함 (생리적)"),
              as_tibble(B) %>% mutate(arm = "AVP 고정 (반사실)"))
  })
  output$plotCF <- renderPlot({
    cf() %>% filter(time <= 10) %>%
      select(time, arm, SODIUM, AVPpg, UOSMOL, EFWCL, STRESSc, DEFICIT) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = recode(name, SODIUM = "혈청 [Na]", AVPpg = "혈장 AVP",
                           UOSMOL = "요삼투압", EFWCL = "EFWC (L/일)",
                           STRESSc = "삼투 스트레스", DEFICIT = "임상 결손")) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c("#b0416a", "#2b6cb0")) +
      labs(x = "일", y = NULL, colour = NULL,
           title = "고장성 식염수는 양쪽 모두 0 mL. 차이는 전부 신장이 만든다.") + THEME
  })
  output$txtCF <- renderText({
    d <- cf()
    s <- d %>% group_by(arm) %>%
      summarise(rise24 = E$na_rise_24(time, SODIUM), str = max(STRESSc),
                def = max(DEFICIT), mye = min(MYE), .groups = "drop")
    paste(apply(s, 1, function(r)
      sprintf("%-24s 24h 상승 %5.1f   최대 스트레스 %5.1f   결손 %5.1f   수초 최저 %5.3f",
              r[1], as.numeric(r[2]), as.numeric(r[3]), as.numeric(r[4]), as.numeric(r[5]))),
      collapse = "\n")
  })

  ## --- Tab 8 : the safety map -------------------------------------------
  smap <- eventReactive(input$goMap, {
    pt <- build_patient("siadh", 21, input$na0)
    grid <- expand.grid(fosm = c(1.00, 0.80, 0.65, 0.55, 0.45),
                        rate = c(4, 6, 8, 10, 12, 14, 16, 20))
    withProgress(message = "안전 지도 계산 중...", {
      grid$stress <- mapply(function(f, r) {
        o <- start_at(mod, pt$state) %>%
          param(c(pt$param, list(CTRLON = 1, RATETGT = r, NASTART = input$na0,
                                 NACAP = 140, WIN = 1, FOSM = f))) %>%
          mrgsim(end = 6, delta = 0.1)
        max(o$STRESSc)
      }, grid$fosm, grid$rate)
    })
    grid
  })
  output$plotMap <- renderPlot({
    ggplot(smap(), aes(factor(rate), factor(fosm), fill = stress)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = sprintf("%.1f", stress)), size = 3.4) +
      scale_fill_gradient2(low = "#cfe9da", mid = "#f7e7c2", high = "#e79a9a",
                           midpoint = 8, name = "최대 스트레스") +
      labs(x = "목표 교정 속도 (mmol/L/24 h)",
           y = "FOSM (유기 삼투질 수송 용량)",
           title = "역치 8.0을 넘는 칸이 위험하다 — 지침의 두 숫자가 여기서 나온다") + THEME
  })
  output$tblMap <- renderDT({
    smap() %>% mutate(안전 = ifelse(stress < 8, "안전", "역치 초과")) %>%
      datatable(rownames = FALSE, options = list(pageLength = 10))
  })

  ## --- Tab 9 : the relowering deadline ------------------------------------
  resc <- eventReactive(input$goResc, {
    pt <- build_patient("hypovol", 7, input$na0)
    pA <- c(pt$param, list(NALOSS = 0, KLOSS = 0, WLOSS = 0, TLOSSEND = 0, R09 = 2))
    Ts <- c(8, 12, 18, 24, 36, 48, 72)
    withProgress(message = "마감시한 계산 중...", {
      out <- lapply(Ts, function(T) {
        th <- T / 24
        a <- start_at(mod, pt$state) %>% param(pA) %>% mrgsim(end = th, delta = 0.05)
        st2 <- as.numeric(tail(a, 1)[, names(init(mod))])
        b <- start_at(mod, st2) %>%
          param(c(pA, list(RESCUE = 1, TRESCUE = th, NARES = 118,
                           DURRES = 1, WIN = 0.5, R09 = 0))) %>%
          mrgsim(events = E$ddavp_q8(th, th + 2.5), end = 90, delta = 0.25)
        tibble(T = T, deficit = max(b$DEFICIT), astro = min(b$AST))
      })
      none <- start_at(mod, pt$state) %>% param(pA) %>% mrgsim(end = 90, delta = 0.25)
      bind_rows(bind_rows(out),
                tibble(T = NA, deficit = max(none$DEFICIT), astro = min(none$AST)))
    })
  })
  output$plotResc <- renderPlot({
    d <- resc(); n <- d$deficit[is.na(d$T)]
    ggplot(filter(d, !is.na(T)), aes(T, deficit)) +
      geom_hline(yintercept = n, linetype = 2, colour = "#b0416a") +
      geom_line(linewidth = 1, colour = "#2b6cb0") + geom_point(size = 3) +
      annotate("text", x = 60, y = n, vjust = -0.6, colour = "#b0416a",
               label = "재하강 없음") +
      labs(x = "재하강 시작 시각 (시간)", y = "최대 임상 결손",
           title = "마감시한은 별아교세포 사멸의 시간상수가 정한다") + THEME
  })
  output$tblResc <- renderDT(datatable(resc(), rownames = FALSE,
                                       options = list(dom = "t")))

  ## --- Tab 10 : scenario comparison -----------------------------------------
  cmp <- eventReactive(input$goCmp, {
    pt <- patient(); base <- pt$param
    if (PHENO[[input$pheno]]$kind == "hypovol")
      base <- c(base, list(NALOSS = 0, KLOSS = 0, WLOSS = 0, TLOSSEND = 0))
    spec <- list(
      "3% NaCl +6/일"        = list(p = list(CTRLON=1, RATETGT=6,  NASTART=input$na0, NACAP=135, WIN=1), e = NULL),
      "3% NaCl +12/일"       = list(p = list(CTRLON=1, RATETGT=12, NASTART=input$na0, NACAP=135, WIN=1), e = NULL),
      "3% NaCl +20/일"       = list(p = list(CTRLON=1, RATETGT=20, NASTART=input$na0, NACAP=135, WIN=1), e = NULL),
      "0.9% NaCl 2 L/일"     = list(p = list(R09 = 2), e = NULL),
      "수분제한 1 L/일"       = list(p = list(WIN = 1.0), e = NULL),
      "톨밥탄 15 mg"          = list(p = list(), e = E$tolvaptan(7)),
      "DDAVP 클램프 + 3% NaCl"= list(p = list(CTRLON=1, RATETGT=6, NASTART=input$na0, NACAP=135, WIN=1),
                                     e = E$ddavp_q8(0, 6)))
    withProgress(message = "시나리오 비교 중...", {
      bind_rows(lapply(input$cmp, function(k) {
        s <- spec[[k]]
        as_tibble(start_at(mod, pt$state) %>%
          param(c(base, s$p, list(FOSM = input$fosm, FNUT = input$fnut))) %>%
          mrgsim(events = s$e, end = input$tend, delta = 0.25)) %>% mutate(arm = k)
      }))
    })
  })
  output$plotCmp <- renderPlot({
    cmp() %>% select(time, arm, SODIUM, EFWCL, OMEGAc, STRESSc, MYE, DEFICIT) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = recode(name, SODIUM = "혈청 [Na]", EFWCL = "EFWC (L/일)",
                           OMEGAc = "Ω", STRESSc = "삼투 스트레스",
                           MYE = "수초 함량", DEFICIT = "임상 결손")) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "일", y = NULL, colour = NULL) + THEME
  })
  output$tblCmp <- renderDT({
    cmp() %>% group_by(arm) %>%
      summarise(`24h 상승` = round(E$na_rise_24(time, SODIUM), 1),
                `[Na] 48h` = round(approx(time, SODIUM, 2)$y, 1),
                `최대 Ω` = round(max(OMEGAc), 1),
                `최대 스트레스` = round(max(STRESSc), 1),
                `수초 최저` = round(min(MYE), 3),
                `최대 결손` = round(max(DEFICIT), 1), .groups = "drop") %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## --- Tab 11 -----------------------------------------------------------------
  output$plotBio <- renderPlot(
    gg(sim(), c("MYOINOS","BRAINH2O","MRIPOS","DEFICIT"),
       c("¹H-MRS myo-inositol (측정 가능)", "뇌 수분 (측정 어려움)",
         "MRI 신호 (늦다)", "임상 결손")))
  output$txtFit <- renderUI(HTML("
<table class='table table-sm'>
<tr><th>정상 생리에 맞춘 것</th><td>수분·용질 정상상태(건강한 모델이 정확한 정상상태가 되도록 WIN을 풀었다: max|dy/dt| = 2.5e-11), Edelman 회귀식, 요 농축 범위, 혈장 요소, AVP 삼투 역치와 기울기</td></tr>
<tr><th>적응 데이터에 맞춘 것</th><td>삼투반응계수 β<sub>i</sub> ([Na] 110 적응 시 총 유기 삼투질 −40%, myo-inositol −66%), 유출/유입 시간상수 (myo-inositol 회복 ~5–6일)</td></tr>
<tr><th>손상 용량–반응에 맞춘 것</th><td>KINJ · KAST · KOLI · KDEM 네 개. 역치 Ω* = 8은 <b>정상위험 한계가 10–12에 오도록</b> 한 번만 보정했다.</td></tr>
<tr><th style='color:#2f7a53'>맞추지 않은 것 (= 예측)</th><td>고위험 한계 6–8 (FOSM 0.55의 결과), 급성/만성 비대칭, 신장이 만드는 자율적 과교정, DDAVP 클램프 효과 크기, 재하강 마감시한, Furst 비 규칙, 요삼투압 &gt; 308일 때 생리식염수의 역효과, MRI 지연</td></tr>
<tr><th style='color:#b0416a'>가장 노출된 가정</th><td>요소는 삼투압 손상 자체(Ω)를 줄이지 못하고 BBB·미세아교세포 팔에서만 작용한다. 교정 속도를 일치시킨 실험에서 요소가 별아교세포 사멸 초기 지표를 줄이면 이 모델은 틀렸다.</td></tr>
<tr><th>알려진 한계</th><td>R이 빌드 컨테이너에 없어 이 모델은 <b>식 검증은 되었으나 컴파일 검증은 되지 않았다</b>. 모든 ODE는 Python/scipy로 독립 재구현하여 적분했다(ods_verify_python.py). Ω는 살아 있는 사람 뇌에서 측정된 적이 없다.</td></tr>
</table>"))
}

shinyApp(ui, server)
