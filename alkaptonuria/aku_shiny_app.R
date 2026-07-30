## ============================================================================
##  Alkaptonuria (AKU) QSP — Shiny dashboard
##  12 tabs.  The app is organised around the three balances that the model is
##  built on, because those are what a user needs to be able to interrogate:
##
##    Balance 1  the flux is conserved; nitisinone changes its EXIT, not its
##               size, so DOSE sets HGA and DIET sets tyrosine
##    Balance 2  the toxic branch is ~1.5e-5 of the flux and is one-way, so
##               benefit is bounded by the integral not yet accrued
##    Balance 3  HGA clearance sits at the renal-plasma-flow ceiling, so the
##               trial endpoint (urinary HGA) and the causal quantity (serum
##               HGA) diverge on treatment
##
##  Run with:  shiny::runApp("aku_shiny_app.R")
##  Requires:  shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  The model file must sit alongside this one.
## ============================================================================

suppressPackageStartupMessages({
  library(shiny); library(mrgsolve); library(dplyr); library(tidyr)
  library(ggplot2); library(DT)
})

source("aku_mrgsolve_model.R", local = TRUE)

AKU_MOD <- aku_model()

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10),
        legend.position = "bottom")

PAL <- c("#e63946", "#457b9d", "#2a9d8f", "#e9c46a", "#7b2cbf",
         "#fb8500", "#8d99ae", "#0288d1")

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel(
    div(
      h3("알캅톤뇨증 QSP 대시보드 · Alkaptonuria QSP Dashboard"),
      p(style = "color:#555;font-size:13px;margin-top:-6px;",
        HTML("HGD 결손 · ochronosis as a one-way integral · nitisinone as a ",
             "<b>flux-redirection</b> agent. ",
             "질병은 보존되는 플럭스의 ",
             "1.5e-5에 해당하는 비가역 분기의 ",
             "<b>적분</b>이다."))
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      sliderInput("bmi", "BMI", 18, 40, 25, 0.5),
      sliderInput("occup", "직업적 부하 (occupational load)", 0.8, 1.8, 1.0, 0.1),
      radioButtons("sexm", "성별 (sex)", c("남성 male" = 1, "여성 female" = 0),
                   selected = 1, inline = TRUE),
      sliderInput("resact", "잔존 HGD 활성 % (residual HGD)", 0, 5, 0, 0.5),
      sliderInput("ckdx", "신기능 저하 배수 (renal decline x)", 1, 4, 1, 0.5),

      hr(), h4("치료 (Treatment)"),
      sliderInput("dose", "니티시논 (mg/day)", 0, 20, 10, 1),
      sliderInput("start", "개시 연령 (age at initiation, y)", 2, 65, 25, 1),
      sliderInput("stop", "중단 연령 (age at stop, y; 70 = 계속)", 5, 70, 70, 1),
      sliderInput("adh", "복약순응도 (adherence)", 0.4, 1, 1, 0.05),

      hr(), h4("식이 (Diet) — 티로신의 유일한 조절 지점"),
      sliderInput("prot", "단백 섭취 (g/day)", 35, 120, 70, 1),
      sliderInput("aasupp", "Phe/Tyr-무함유 보충 (g/day)", 0, 40, 0, 1),

      hr(), h4("보조 치료 (Adjuncts)"),
      checkboxInput("nsaid", "NSAID", FALSE),
      checkboxInput("physio", "물리치료·체중관리 (physiotherapy)", FALSE),
      checkboxInput("alkali", "요 알칼리화·수분 (urine alkalinisation)", FALSE),
      sliderInput("casc", "혈중 아스코르브산 (umol/L)", 30, 150, 50, 5),
      checkboxInput("ideal", "IDEAL 비교약: HGA 차단 + 티로신 상승 없음", FALSE),

      hr(),
      sliderInput("endage", "시뮬레이션 종료 연령", 40, 85, 70, 1),
      actionButton("go", "실행 (Run)", class = "btn-primary btn-block"),
      helpText(style = "font-size:11px;",
               "70년 적분은 수 초가 걸립니다. 반대편 대조군(무치료)은 항상 함께 계산됩니다.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        tabPanel("① 환자 프로파일",
          br(), fluidRow(column(12, uiOutput("kpis"))),
          br(), plotOutput("p_overview", height = "420px"),
          br(), h5("현재 설정 요약"), DTOutput("t_profile")),

        tabPanel("② 니티시논 PK",
          br(), plotOutput("p_pk", height = "300px"),
          plotOutput("p_pkss", height = "260px"),
          helpText(HTML("반감기가 약 54시간이므로 <b>하루 1회 투여와 2회 분할 투여는 ",
                        "거의 구분되지 않는다</b>. 반감기가 투여 간격보다 훨씬 길기 때문이다.")),
          br(), DTOutput("t_pk")),

        tabPanel("③ 경로 플럭스 (Balance 1)",
          br(), h5("들어온 것은 반드시 나간다 — 출구 구성의 변화"),
          plotOutput("p_flux", height = "330px"),
          plotOutput("p_exit", height = "300px"),
          br(), DTOutput("t_mass")),

        tabPanel("④ HGA: 요 vs 혈청 (Balance 3)",
          br(), h5("시험 평가변수와 인과량의 괴리"),
          plotOutput("p_hga", height = "330px"),
          plotOutput("p_gap", height = "300px"),
          helpText(HTML("연골을 물들이는 것은 요중 HGA가 아니라 <b>혈장 HGA</b>다. ",
                        "치료 중 HPPA와 HPLA가 14배 이상 올라가 같은 유기음이온 분비 경로를 ",
                        "경쟁적으로 점유하므로 HGA 클리어런스가 오히려 감소하고, 그 결과 ",
                        "혈장 농도는 소변 배설량보다 훨씬 덜 내려간다."))),

        tabPanel("⑤ 티로시네미아 (용량 vs 식이)",
          br(), h5("2x2: 용량은 HGA를, 식이가 티로신을 정한다"),
          plotOutput("p_tyr", height = "320px"),
          plotOutput("p_2x2", height = "330px"),
          br(), DTOutput("t_tyr")),

        tabPanel("⑥ 색소 적분 (Balance 2)",
          br(), plotOutput("p_pig", height = "330px"),
          plotOutput("p_pigdist", height = "300px"),
          helpText(HTML("관절 연골·디스크·대동맥판의 색소에는 <b>소실항이 없다</b>. ",
                        "그 콜라겐이 교체되지 않기 때문이다. 피부·공막·귀 연골에만 ",
                        "콜라겐 교체 속도에 해당하는 소실항을 두었으며, 보고된 ",
                        "색소 역전이 이 조직들에 국한된다는 것이 모델의 예측이다."))),

        tabPanel("⑦ 관절·척추",
          br(), plotOutput("p_joint", height = "330px"),
          plotOutput("p_spine", height = "300px"),
          br(), DTOutput("t_struct")),

        tabPanel("⑧ 심장·신장·기타 장기",
          br(), plotOutput("p_valve", height = "300px"),
          plotOutput("p_organ", height = "320px")),

        tabPanel("⑨ 임상 엔드포인트·AKUSSI",
          br(), plotOutput("p_akussi", height = "340px"),
          plotOutput("p_hazard", height = "300px"),
          br(), DTOutput("t_end")),

        tabPanel("⑩ 안전성: 각막병증",
          br(), plotOutput("p_kerato", height = "330px"),
          helpText(HTML("관리 임계값: 티로신 700 umol/L 초과 시 단백 제한 강화, ",
                        "900 umol/L 초과 시 Phe/Tyr-무함유 아미노산 보충 추가.")),
          plotOutput("p_safety", height = "290px")),

        tabPanel("⑪ 개시 연령과 헤드룸",
          br(), h5("같은 약, 같은 용량 — 남은 적분만 다르다"),
          actionButton("go_hr", "헤드룸 스캔 실행 (13개 개시 연령)",
                       class = "btn-warning"),
          br(), br(), plotOutput("p_headroom", height = "340px"),
          plotOutput("p_headroom2", height = "300px"),
          br(), DTOutput("t_headroom")),

        tabPanel("⑫ 시나리오 비교·검증",
          br(), h5("사전 정의된 24개 시나리오와 held-out 검증"),
          actionButton("go_scn", "시나리오 세트 실행", class = "btn-warning"),
          actionButton("go_val", "held-out 검증 실행", class = "btn-info"),
          br(), br(), plotOutput("p_scn", height = "340px"),
          br(), h5("held-out 검증 (보정에 쓰이지 않은 값)"),
          DTOutput("t_val"),
          br(), h5("시나리오 종료 시점 요약"), DTOutput("t_scn"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  pars <- reactive({
    list(BMI = input$bmi, OCCUP = input$occup, SEXM = as.numeric(input$sexm),
         RESACT = input$resact/100, CKDX = input$ckdx,
         PROT = input$prot, AASUPP = input$aasupp,
         FBIO = input$adh, CASC = input$casc,
         NSAIDON = as.numeric(input$nsaid), PHYSIO = as.numeric(input$physio),
         ALKALI = as.numeric(input$alkali), IDEALDRG = as.numeric(input$ideal))
  })

  ## treated arm and its matched untreated control -- always both, so that no
  ## conclusion on this dashboard can come from an unmatched comparison
  sims <- eventReactive(input$go, {
    withProgress(message = "적분 중 (integrating 70 years)...", value = 0.2, {
      tx <- sim_aku(AKU_MOD, dose_mg = input$dose, start_age = input$start,
                    stop_age = if (input$stop >= 70) Inf else input$stop,
                    end_age = input$endage, pars = pars(), label = "treated")
      incProgress(0.5)
      ct <- sim_aku(AKU_MOD, dose_mg = 0, start_age = input$start,
                    end_age = input$endage, pars = pars(), label = "control")
      list(tx = tx, ct = ct, both = bind_rows(tx, ct))
    })
  }, ignoreNULL = FALSE)

  fin <- reactive(tail(sims()$tx, 1))
  finc <- reactive(tail(sims()$ct, 1))

  ## ---------------------------------------------------------------- ① profile
  output$kpis <- renderUI({
    f <- fin(); c0 <- finc()
    box <- function(lab, val, sub = "", col = "#457b9d")
      div(style = paste0("display:inline-block;width:15.5%;margin:3px;padding:9px;",
                         "border-left:4px solid ", col, ";background:#fafafa;"),
          div(style = "font-size:10.5px;color:#666;", lab),
          div(style = "font-size:19px;font-weight:bold;", val),
          div(style = "font-size:10px;color:#888;", sub))
    tagList(
      box("혈청 HGA", sprintf("%.2f", f$CHGAo), "umol/L (무치료 대비)", "#e63946"),
      box("요중 HGA", sprintf("%.0f", f$UHGA24), "umol/day", "#c1121f"),
      box("혈장 티로신", sprintf("%.0f", f$CTYRo), "umol/L", "#e9c46a"),
      box("cAKUSSI", sprintf("%.1f", f$CAKUSSI), sprintf("무치료 %.1f", c0$CAKUSSI), "#212121"),
      box("색소 절감", sprintf("%.0f%%", 100*f$PIGSPARE), "counterfactual 대비", "#7b2cbf"),
      box("관절치환 확률", sprintf("%.0f%%", 100*f$PJR), sprintf("무치료 %.0f%%", 100*c0$PJR), "#2a9d8f")
    )
  })

  output$p_overview <- renderPlot({
    d <- sims()$both %>%
      select(AGEYo, scenario, CHGAo, CTYRo, CART, DISCH, CAKUSSI, PAINVAS) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(CHGAo = "혈청 HGA (umol/L)", CTYRo = "혈장 Tyr (umol/L)",
             CART = "온전한 연골 (분율)", DISCH = "디스크 높이 (분율)",
             CAKUSSI = "cAKUSSI", PAINVAS = "통증 VAS")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) +
      geom_vline(xintercept = input$start, linetype = 3, colour = "grey40") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "치료군과 짝지어진 무치료 대조군",
           subtitle = paste0("점선 = 치료 개시 (", input$start, "세). ",
                             "무치료군은 모든 다른 설정이 동일합니다."),
           x = "연령 (years)", y = NULL, colour = NULL) + THEME
  })

  output$t_profile <- renderDT({
    p <- pars()
    datatable(data.frame(parameter = names(p), value = unlist(p)),
              options = list(dom = "t", pageLength = 20), rownames = FALSE)
  })

  ## --------------------------------------------------------------------- ② PK
  output$p_pk <- renderPlot({
    st <- input$start*365.25
    mo <- aku_init(AKU_MOD) %>% param(pars())
    e  <- ev(time = st, amt = input$dose*NT_MG_TO_UMOL, ii = 1, addl = 60)
    o  <- mo %>% ev(e) %>%
      mrgsim(start = st, end = st + 62, delta = 0.05,
             atol = 1e-12, rtol = 1e-10, maxsteps = 2e6) %>% as_tibble()
    ggplot(o, aes((TIME-st), CNTo)) + geom_line(colour = PAL[2], linewidth = 0.8) +
      labs(title = "니티시논 혈장 농도 — 60일 축적",
           subtitle = "t1/2 약 54시간이므로 정상상태 도달에 1-2주가 걸린다",
           x = "치료 개시 후 일수", y = "니티시논 (umol/L)") + THEME
  })

  output$p_pkss <- renderPlot({
    d <- bind_rows(lapply(c(1, 2, 4, 8, 10, 20), function(dd) {
      b <- probe_biochem(AKU_MOD, dd, 0.5, input$start, pars())
      tibble(dose = dd, CNT = b$CNTo, uHGA = b$UHGA24, sHGA = b$CHGAo, Tyr = b$CTYRo)
    }))
    ggplot(d, aes(dose, CNT)) +
      geom_line(colour = PAL[2]) + geom_point(size = 2.6, colour = PAL[2]) +
      labs(title = "용량–정상상태 농도는 선형",
           subtitle = "PK는 선형인데도 뒤따르는 HGA 억제는 그렇지 않다 (다음 탭)",
           x = "니티시논 (mg/day)", y = "정상상태 니티시논 (umol/L)") + THEME
  })

  output$t_pk <- renderDT({
    d <- bind_rows(lapply(c(0, 1, 2, 4, 8, 10, 20), function(dd) {
      b <- probe_biochem(AKU_MOD, dd, 1, input$start, pars())
      tibble(`용량 mg/day` = dd, `니티시논 umol/L` = round(b$CNTo, 2),
             `요중 HGA umol/day` = round(b$UHGA24),
             `혈청 HGA umol/L` = round(b$CHGAo, 2),
             `혈장 Tyr umol/L` = round(b$CTYRo))
    }))
    datatable(d, options = list(dom = "t"), rownames = FALSE)
  })

  ## ------------------------------------------------------------------- ③ flux
  output$p_flux <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, FLUXIN, EXITSUM) %>%
      pivot_longer(-c(AGEYo, scenario))
    ggplot(d, aes(AGEYo, value/1000, colour = name, linetype = scenario)) +
      geom_line(linewidth = 0.85) +
      scale_colour_manual(values = c(FLUXIN = "#2a9d8f", EXITSUM = "#e63946"),
                          labels = c(EXITSUM = "측정된 출구 합", FLUXIN = "식이 Phe+Tyr 유입")) +
      labs(title = "질량보존 점검: 유입과 출구 합",
           subtitle = "니티시논은 이 총량을 바꾸지 않는다. 출구의 구성만 바꾼다.",
           x = "연령 (years)", y = "mmol/day", colour = NULL, linetype = NULL) + THEME
  })

  output$p_exit <- renderPlot({
    d <- sims()$both %>%
      select(AGEYo, scenario, UHGA24, UTYR24, UHPP24, UHPL24, UCON24) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(UHGA24 = "요중 HGA", UTYR24 = "요중 티로신", UHPP24 = "요중 HPPA",
             UHPL24 = "요중 HPLA", UCON24 = "티로신 포합체")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value/1000, fill = name)) +
      geom_area(position = "stack") + facet_wrap(~scenario) +
      scale_fill_manual(values = PAL[1:5]) +
      labs(title = "출구 구성 — 니티시논이 실제로 하는 일",
           subtitle = "HGA 출구를 닫으면 같은 양이 HPPA·HPLA·티로신·포합체로 나간다",
           x = "연령 (years)", y = "mmol/day", fill = NULL) + THEME
  })

  output$t_mass <- renderDT({
    f <- fin(); c0 <- finc()
    d <- tibble(
      `항목` = c("식이 Phe+Tyr 유입", "요중 HGA", "요중 티로신", "요중 HPPA",
               "요중 HPLA", "티로신 포합체", "출구 합", "출구/유입"),
      `치료` = c(f$FLUXIN, f$UHGA24, f$UTYR24, f$UHPP24, f$UHPL24, f$UCON24,
               f$EXITSUM, f$EXITSUM/f$FLUXIN),
      `무치료` = c(c0$FLUXIN, c0$UHGA24, c0$UTYR24, c0$UHPP24, c0$UHPL24,
                 c0$UCON24, c0$EXITSUM, c0$EXITSUM/c0$FLUXIN))
    d$`치료` <- ifelse(d$`항목` == "출구/유입", round(d$`치료`, 3), round(d$`치료`))
    d$`무치료` <- ifelse(d$`항목` == "출구/유입", round(d$`무치료`, 3), round(d$`무치료`))
    datatable(d, options = list(dom = "t"), rownames = FALSE)
  })

  ## ------------------------------------------------------------------- ④ HGA
  output$p_hga <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, CHGAo, UHGA24)
    ggplot(d) +
      geom_line(aes(AGEYo, CHGAo, colour = scenario), linewidth = 0.9) +
      geom_line(aes(AGEYo, UHGA24/1000, colour = scenario), linetype = 2) +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "혈청 HGA (실선, umol/L) vs 요중 HGA (점선, mmol/day)",
           subtitle = paste("두 곡선의 감소폭이 다르다는 것이 이 모델의 핵심 주장이다.",
                            "색소 침착 속도는 실선을 따른다."),
           x = "연령 (years)", y = NULL, colour = NULL) + THEME
  })

  output$p_gap <- renderPlot({
    g <- analysis_endpoint_gap(AKU_MOD, pars())
    d <- g %>% select(dose_mg, uHGA_pct_drop, sHGA_pct_drop) %>%
      pivot_longer(-dose_mg)
    ggplot(d, aes(dose_mg, value, colour = name)) +
      geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
      scale_colour_manual(values = c(uHGA_pct_drop = "#c1121f", sHGA_pct_drop = "#457b9d"),
                          labels = c(sHGA_pct_drop = "혈청 HGA 감소 %",
                                     uHGA_pct_drop = "요중 HGA 감소 % (시험 평가변수)")) +
      labs(title = "평가변수는 인과량보다 좋아 보인다",
           subtitle = "요중 감소율이 혈청 감소율을 체계적으로 앞선다",
           x = "니티시논 (mg/day)", y = "무치료 대비 감소 (%)", colour = NULL) + THEME
  })

  ## ------------------------------------------------------------------- ⑤ Tyr
  output$p_tyr <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, CTYRo, CPHEo, CHPPo, CHPLAo) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(CTYRo = "티로신", CPHEo = "페닐알라닌", CHPPo = "HPPA", CHPLAo = "HPLA")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = name, linetype = scenario)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = c(700, 900), linetype = 3, colour = "grey50") +
      scale_colour_manual(values = PAL[c(4,3,2,5)]) +
      labs(title = "상류 대사물의 역압 (back-pressure)",
           subtitle = "점선 = 티로신 700 / 900 umol/L 관리 임계값",
           x = "연령 (years)", y = "umol/L", colour = NULL, linetype = NULL) + THEME
  })

  output$p_2x2 <- renderPlot({
    d <- analysis_dose_vs_diet(AKU_MOD, pars())
    ggplot(d, aes(dose_mg, sTYR, colour = factor(protein_g))) +
      geom_line(linewidth = 0.9) + geom_point(size = 2) +
      geom_hline(yintercept = 700, linetype = 3) +
      scale_colour_manual(values = PAL[1:5], name = "단백 g/day") +
      labs(title = "티로신은 용량이 아니라 식이가 정한다",
           subtitle = paste("가로축을 20배 움직여도 곡선은 거의 평평하고,",
                            "곡선 사이의 간격(식이)이 전부다"),
           x = "니티시논 (mg/day)", y = "혈장 티로신 (umol/L)") + THEME
  })

  output$t_tyr <- renderDT({
    d <- analysis_dose_vs_diet(AKU_MOD, pars()) %>%
      mutate(across(c(sHGA, sTYR, uHGA24, pigment_rate), ~round(.x, 3)))
    datatable(d, options = list(pageLength = 12), rownames = FALSE)
  })

  ## ---------------------------------------------------------------- ⑥ pigment
  output$p_pig <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, PIGTOT, CFPIG, PDCo, BRCo) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(PIGTOT = "총 색소 (umol)", CFPIG = "무치료 counterfactual (연골, umol)",
             PDCo = "연골 색소 밀도", BRCo = "연골 취성화 (0-1)")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "색소는 적분이고, 치료는 미래의 적분만 건드린다",
           subtitle = "취성화는 색소 밀도의 급한 Hill 함수이므로 3번째 10년대에 갑자기 나타난다",
           x = "연령 (years)", y = NULL, colour = NULL) + THEME
  })

  output$p_pigdist <- renderPlot({
    f <- fin()
    d <- tibble(depot = c("관절 연골", "디스크", "대동맥판", "힘줄",
                          "공막", "귀", "피부", "전립선·신장·기타"),
                umol = c(f$PCART, f$PDISC, f$PVALV, f$PTEND,
                         f$PSCL, f$PEAR, f$PSKIN, f$POTH),
                turnover = c("없음", "없음", "없음", "없음",
                             "있음", "있음", "있음", "없음"))
    ggplot(d, aes(reorder(depot, umol), umol, fill = turnover)) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = c("없음" = "#212529", "있음" = "#457b9d"),
                        name = "콜라겐 교체") +
      labs(title = paste0("색소 침착 분포 (", round(f$AGEYo), "세 시점)"),
           subtitle = "교체가 있는 조직에서만 색소가 옅어질 수 있다 — 관절 연골에는 소실항이 없다",
           x = NULL, y = "umol HGA-equivalents") + THEME
  })

  ## ------------------------------------------------------------------ ⑦ joint
  output$p_joint <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, CART, FRAGC, SYN, SUBCH, OSTEO) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(CART = "온전한 연골", FRAGC = "관절내 색소 파편", SYN = "활막염",
             SUBCH = "연골하 경화", OSTEO = "골극")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "증폭 루프: 파편 → 활막염 → MMP → 더 많은 연골 소실",
           subtitle = "가속이 파라미터가 아니라 이 루프에서 나온다",
           x = "연령 (years)", y = NULL, colour = NULL) + THEME
  })

  output$p_spine <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, DISCH, DCALC, ANKY) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(DISCH = "디스크 높이", DCALC = "디스크 석회화", ANKY = "척추 강직 지수")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "척추 — 색소의 가장 이른 방사선학적 흔적",
           x = "연령 (years)", y = NULL, colour = NULL) + THEME
  })

  output$t_struct <- renderDT({
    f <- fin(); c0 <- finc()
    datatable(tibble(
      `지표` = c("온전한 연골", "활막염", "디스크 높이", "디스크 석회화",
               "척추 강직", "힘줄 온전성", "골밀도"),
      `치료` = round(c(f$CART, f$SYN, f$DISCH, f$DCALC, f$ANKY, f$TENDI, f$BMD), 3),
      `무치료` = round(c(c0$CART, c0$SYN, c0$DISCH, c0$DCALC, c0$ANKY, c0$TENDI, c0$BMD), 3)),
      options = list(dom = "t"), rownames = FALSE)
  })

  ## ------------------------------------------------------------------ ⑧ organ
  output$p_valve <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, PMAXS, AVA, LVMI) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(PMAXS = "최대 압력차 (mmHg)", AVA = "판막 면적 (cm2)", LVMI = "LV 질량 초과")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "대동맥판 — 의도적으로 가역 채널에서 제외된 경로",
           subtitle = "SONIA 2에서 니티시논은 판막 진행을 늦추지 못했다 (p=0.53)",
           x = "연령 (years)", y = NULL, colour = NULL) + THEME
  })

  output$p_organ <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, STONE, PSTONE, RFUN, HEAR, SKINP, TENDI) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(STONE = "신결석 부담", PSTONE = "전립선 결석", RFUN = "상대 신혈류량",
             HEAR = "청력 역치 이동 (dB)", SKINP = "가시 색소", TENDI = "힘줄 온전성")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "장기별 결과 — 신기능 저하는 혈장 HGA를 올려 되먹임한다",
           x = "연령 (years)", y = NULL, colour = NULL) + THEME
  })

  ## ----------------------------------------------------------------- ⑨ AKUSSI
  output$p_akussi <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, CAKUSSI, AKJ, AKS, AKC, PAINVAS) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(CAKUSSI = "cAKUSSI (합)", AKJ = "관절 도메인", AKS = "척추 도메인",
             AKC = "임상 도메인", PAINVAS = "통증 VAS")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y", ncol = 3) +
      geom_vline(xintercept = input$start, linetype = 3, colour = "grey40") +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "cAKUSSI와 그 구성 도메인",
           subtitle = "절대값은 임상 척도와 직접 비교할 수 없다 — 기울기만 비교하라",
           x = "연령 (years)", y = NULL, colour = NULL) + THEME
  })

  output$p_hazard <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, PJR, PAVR, PRUP, PKER) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(PJR = "관절치환", PAVR = "대동맥판치환", PRUP = "힘줄 파열",
             PKER = "각막병증")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = name, linetype = scenario)) +
      geom_line(linewidth = 0.85) +
      scale_colour_manual(values = PAL[c(1,2,3,6)]) +
      labs(title = "누적 사건 확률",
           subtitle = "각막병증만 치료군에서 올라간다 — 나머지는 모두 내려간다",
           x = "연령 (years)", y = "누적 확률", colour = NULL, linetype = NULL) + THEME
  })

  output$t_end <- renderDT({
    f <- fin(); c0 <- finc()
    datatable(tibble(
      `엔드포인트` = c("cAKUSSI", "통증 VAS", "관절치환 확률", "판막치환 확률",
                     "힘줄파열 확률", "각막병증 확률", "색소 절감 분율",
                     "회피된 HGA 노출 (umol/L*day)"),
      `치료` = round(c(f$CAKUSSI, f$PAINVAS, f$PJR, f$PAVR, f$PRUP, f$PKER,
                     f$PIGSPARE, f$AVOIDI), 3),
      `무치료` = round(c(c0$CAKUSSI, c0$PAINVAS, c0$PJR, c0$PAVR, c0$PRUP,
                       c0$PKER, 0, 0), 3)),
      options = list(dom = "t"), rownames = FALSE)
  })

  ## ---------------------------------------------------------------- ⑩ kerato
  output$p_kerato <- renderPlot({
    d <- sims()$both %>% select(AGEYo, scenario, CTYRo, CORTYR, PKER) %>%
      pivot_longer(-c(AGEYo, scenario))
    lab <- c(CTYRo = "혈장 티로신 (umol/L)", CORTYR = "각막 결정 부하",
             PKER = "누적 각막병증 확률")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(AGEYo, value, colour = scenario)) +
      geom_line(linewidth = 0.85) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c(control = "#8d99ae", treated = "#e63946")) +
      labs(title = "각막병증은 용량의 문제가 아니라 식이 관리의 문제다",
           subtitle = "왼쪽 단백 섭취 슬라이더를 내리면 세 패널이 모두 내려간다",
           x = "연령 (years)", y = NULL, colour = NULL) + THEME
  })

  output$p_safety <- renderPlot({
    grid <- expand.grid(dose = c(2, 10, 20), prot = c(42, 56, 70, 84, 105))
    d <- bind_rows(Map(function(dd, pp) {
      b <- probe_biochem(AKU_MOD, dd, 3, input$start,
                         modifyList(pars(), list(PROT = pp)))
      tibble(dose = dd, protein = pp, Tyr = b$CTYRo, crystal = b$CORTYR)
    }, grid$dose, grid$prot))
    ggplot(d, aes(protein, Tyr, colour = factor(dose))) +
      geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
      geom_hline(yintercept = c(700, 900), linetype = 3) +
      scale_colour_manual(values = PAL[1:3], name = "mg/day") +
      labs(title = "단백 섭취 대 티로신, 세 가지 용량에서",
           subtitle = "세 선이 거의 겹친다 — 이것이 이 모델의 가장 검증 가능한 예측이다",
           x = "식이 단백 (g/day)", y = "혈장 티로신 (umol/L)") + THEME
  })

  ## --------------------------------------------------------------- ⑪ headroom
  HR <- eventReactive(input$go_hr, {
    withProgress(message = "13개 개시 연령 스캔 중...", value = 0.3,
                 analysis_headroom(AKU_MOD, pars()))
  })

  output$p_headroom <- renderPlot({
    d <- HR() %>% select(init_age, spared_frac, cart_70, p_joint_repl) %>%
      pivot_longer(-init_age)
    lab <- c(spared_frac = "절감된 색소 분율", cart_70 = "70세 온전 연골",
             p_joint_repl = "70세 관절치환 확률")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(init_age, value, colour = name)) +
      geom_line(linewidth = 1) + geom_point(size = 2.4) +
      scale_colour_manual(values = PAL[c(5,3,1)]) +
      labs(title = "동일한 약·동일한 용량 — 남은 적분만 다르다",
           subtitle = "모든 곡선의 기울기가 곧 '기다림의 비용'이다",
           x = "치료 개시 연령 (years)", y = NULL, colour = NULL) + THEME
  })

  output$p_headroom2 <- renderPlot({
    ggplot(HR(), aes(init_age, cAKUSSI_70)) +
      geom_line(colour = PAL[1], linewidth = 1) + geom_point(size = 2.4, colour = PAL[1]) +
      labs(title = "70세 cAKUSSI 대 개시 연령",
           subtitle = "가로축이 오른쪽으로 갈수록 무치료 궤적에 수렴한다",
           x = "치료 개시 연령 (years)", y = "70세 cAKUSSI") + THEME
  })

  output$t_headroom <- renderDT({
    datatable(HR() %>% mutate(across(where(is.numeric), ~round(.x, 3))),
              options = list(pageLength = 13), rownames = FALSE)
  })

  ## -------------------------------------------------------------- ⑫ scenarios
  SC <- eventReactive(input$go_scn, {
    withProgress(message = "24개 시나리오 실행 중...", value = 0.2,
                 run_scenarios(AKU_MOD))
  })
  VAL <- eventReactive(input$go_val, {
    withProgress(message = "held-out 검증 중...", value = 0.3,
                 validate_aku(AKU_MOD))
  })

  output$p_scn <- renderPlot({
    d <- SC() %>% filter(id %in% c("S01","S09","S10","S11","S12","S13","S23","S24"))
    ggplot(d, aes(AGEYo, CAKUSSI, colour = id)) +
      geom_line(linewidth = 0.9) +
      labs(title = "개시 연령·중단·이상적 비교약 시나리오",
           subtitle = "S01 무치료 · S09-S13 개시 연령 5/15/25/40/55 · S23 45세 중단 · S24 IDEAL",
           x = "연령 (years)", y = "cAKUSSI", colour = NULL) + THEME
  })

  output$t_val <- renderDT({
    datatable(VAL() %>% mutate(across(where(is.numeric), ~signif(.x, 3))),
              options = list(pageLength = 21), rownames = FALSE)
  })

  output$t_scn <- renderDT({
    d <- SC() %>% group_by(id, scenario) %>% slice_tail(n = 1) %>% ungroup() %>%
      select(id, scenario, AGEYo, CTYRo, CHGAo, UHGA24, CAKUSSI, PIGTOT, PJR, PKER) %>%
      mutate(across(where(is.numeric), ~signif(.x, 3)))
    datatable(d, options = list(pageLength = 24, scrollX = TRUE), rownames = FALSE)
  })
}

shinyApp(ui, server)
