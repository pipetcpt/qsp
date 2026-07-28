# =============================================================================
# Drug-Induced Liver Injury (DILI) — QSP simulator (Shiny)
# =============================================================================
#
# A dashboard over dili_mrgsolve_model.R. The app is deliberately organised
# around the model's thesis rather than around its compartments:
#
#   Tab 2  the drug arrives at some RATE
#   Tab 3  the liver disposes of it at some RATE, and that rate has a ceiling
#   Tab 4  when production outruns disposal, a POSITIVE FEEDBACK latches
#   Tab 5  the latch converts stressed hepatocytes into lost mass
#   Tab 6  ALT reports the rate, bilirubin reports the reserve — Hy's Law
#   Tab 7  the bile-acid and immune arms reach the same death node by other routes
#   Tab 8  the antidote window is where the trajectory crosses the separatrix
#   Tab 9  side-by-side scenarios
#
# Run:   R -e "shiny::runApp('dili_shiny_app.R', port=8080)"
# Needs: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
# =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

source("dili_mrgsolve_model.R", local = TRUE)   # defines mod, sim_dili, ...

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93", colour = NA),
        legend.position = "bottom")

PAL <- c("#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD",
         "#8C564B", "#17BECF", "#7F7F7F")

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("약물유발 간손상 (DILI) QSP 시뮬레이터 · Drug-Induced Liver Injury"),
  tags$p(style = "color:#8B0000; font-weight:600; margin-top:-8px;",
         "핵심 명제: 간손상은 용량(dose)이 아니라 속도(rate)의 문제이며, ",
         "JNK–Sab 정귀환 고리 때문에 시스템은 쌍안정(bistable)이다. ",
         "따라서 같은 총 노출량이 무증상일 수도, 치명적일 수도 있다."),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("① 약물 (Drug)"),
      selectInput("archetype", "약물 원형 (archetype)",
                  choices = c("A · 생체활성화형 (APAP)"      = "A",
                              "B · BSEP 억제형 (담즙정체)"   = "B",
                              "C · 특이체질 면역형 (HLA)"    = "C"),
                  selected = "A"),
      numericInput("dose", "1회 용량 (mg/kg)", 350, min = 0, max = 1000, step = 10),
      selectInput("regimen", "투여 방식",
                  choices = c("단회 급성 섭취"        = "single",
                              "6시간 간격 x 7일"      = "q6h7d",
                              "12시간 간격 x 28일"    = "q12h28d",
                              "24시간 간격 x 56일"    = "q24h56d",
                              "지정 시간에 걸쳐 분할" = "spread"),
                  selected = "single"),
      conditionalPanel("input.regimen == 'spread'",
        sliderInput("spread_h", "총 섭취 시간 (h)", 1, 72, 12, step = 1)),

      hr(),
      h4("② 숙주 (Host)"),
      sliderInput("fcyp", "CYP2E1 유도 배수 (만성 음주)", 1, 3, 1, step = 0.1),
      sliderInput("fgsh", "GSH 풀 배수 (금식·영양불량)", 0.4, 1.2, 1, step = 0.05),
      sliderInput("cysbase", "시스테인 공급 set-point", 1.0, 2.0, 1.65, step = 0.05),
      sliderInput("agef", "재생능 (1 = 젊음, 0.5 = 고령)", 0.4, 1.2, 1, step = 0.05),

      hr(),
      h4("③ 약물 부가 책임 (Liabilities)"),
      sliderInput("logki", "BSEP 억제 Ki (log10 µM; 9 = 없음)",
                  -1, 9, 9, step = 0.25),
      sliderInput("func", "탈공역 / β-산화 차단 부담", 0, 0.6, 0, step = 0.05),
      checkboxInput("hla", "HLA 위험 대립유전자 보유", FALSE),
      checkboxInput("ici", "면역관문억제제 병용 (관용 제거)", FALSE),

      hr(),
      h4("④ 치료 (Treatment)"),
      checkboxInput("use_nac", "N-아세틸시스테인 (IV, Prescott)", FALSE),
      conditionalPanel("input.use_nac",
        sliderInput("nac_t", "NAC 개시 시각 (h)", 0, 48, 8, step = 1)),
      checkboxInput("use_ster", "코르티코스테로이드", FALSE),
      conditionalPanel("input.use_ster",
        sliderInput("ster_t", "스테로이드 개시 (일)", 1, 60, 42, step = 1)),

      hr(),
      sliderInput("tend", "시뮬레이션 기간 (h)", 168, 2400, 336, step = 24),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        # ---------------------------------------------------------------- 1
        tabPanel("① 환자 프로파일 / 판정",
          br(),
          fluidRow(
            column(4, wellPanel(h4("결과 판정"), htmlOutput("verdict"))),
            column(4, wellPanel(h4("주요 지표"), tableOutput("keytab"))),
            column(4, wellPanel(h4("이 실행에서 실제로 일어난 일"),
                                htmlOutput("narrative")))
          ),
          hr(),
          h4("생성 속도 대 처리 속도 (production flux vs disposal flux)"),
          p("이 그림이 모델의 전부입니다. 빨간 선(반응성 대사체 생성)이 ",
            "초록 선(GSH 재합성 능력)을 넘어서는 면적이 곧 손상입니다."),
          plotOutput("fluxplot", height = 320)
        ),

        # ---------------------------------------------------------------- 2
        tabPanel("② 약물동태 (PK)",
          br(),
          fluidRow(
            column(6, plotOutput("pk_conc", height = 300)),
            column(6, plotOutput("pk_nomogram", height = 300))
          ),
          hr(),
          h4("대사 경로 배분 — 생체활성화 분율은 용량에 따라 스스로 상승한다"),
          p("황산포합(PAPS)과 글루쿠론산포합(UDPGA)은 보조인자 고갈로 포화되지만 ",
            "CYP 경로는 Km이 높아 거의 선형입니다. 그래서 고용량에서 반응성 ",
            "대사체로 가는 분율이 자동으로 커집니다 — 별도의 스위치는 없습니다."),
          plotOutput("pk_partition", height = 300)
        ),

        # ---------------------------------------------------------------- 3
        tabPanel("③ 해독 능력 (Thiol defence)",
          br(),
          fluidRow(
            column(6, plotOutput("gsh_plot", height = 300)),
            column(6, plotOutput("rm_plot", height = 300))
          ),
          hr(),
          h4("GSH가 게이트다"),
          p("반응성 대사체는 GSH 포합과 단백 공유결합 사이에서 '경쟁 속도'로만 ",
            "나뉩니다. GSH 6.5 mM에서는 약 0.6%만 단백질에 결합하지만 ",
            "0.1 mM에서는 약 25%가 결합합니다. 임계값 파라미터는 없습니다."),
          plotOutput("gate_plot", height = 300)
        ),

        # ---------------------------------------------------------------- 4
        tabPanel("④ 미토콘드리아 · JNK 고리",
          br(),
          fluidRow(
            column(6, plotOutput("mito_plot", height = 300)),
            column(6, plotOutput("jnk_plot", height = 300))
          ),
          hr(),
          h4("쌍안정성과 분리선 (separatrix)"),
          p("아래는 JNK–Sab 고리만 떼어낸 것입니다. 산화환원 여력 ",
            "g = (0.30 + 0.70·GSH/GSH₀)·NRF2 가 떨어지면 OFF 상태가 사라집니다. ",
            "붉은 점선이 분리선이고, 시뮬레이션 궤적이 그 위로 올라가면 손상은 ",
            "돌이킬 수 없습니다."),
          plotOutput("bistab_plot", height = 340)
        ),

        # ---------------------------------------------------------------- 5
        tabPanel("⑤ 간세포 질량 · 재생",
          br(),
          plotOutput("mass_plot", height = 320),
          hr(),
          fluidRow(
            column(6, plotOutput("damp_plot", height = 280)),
            column(6, plotOutput("regen_plot", height = 280))
          )
        ),

        # ---------------------------------------------------------------- 6
        tabPanel("⑥ 임상 검사 · Hy's Law",
          br(),
          fluidRow(
            column(6, plotOutput("lab_enz", height = 300)),
            column(6, plotOutput("lab_fn", height = 300))
          ),
          hr(),
          h4("Hy's Law = 속도(rate) × 예비력(reserve)의 결합 검정"),
          p("ALT는 '죽는 속도'의 저역통과 필터이고, 빌리루빈은 '남은 예비력'이 ",
            "결정합니다(청소율이 생존 질량에 비례). 두 팔을 동시에 요구한다는 ",
            "것은 간이 빨리 죽고 있으면서 동시에 여력이 바닥났음을 요구하는 ",
            "것이며, 그래서 두 팔의 결합만이 사망을 예측합니다."),
          plotOutput("hylaw_plot", height = 340)
        ),

        # ---------------------------------------------------------------- 7
        tabPanel("⑦ 담즙산 · 면역",
          br(),
          fluidRow(
            column(6, plotOutput("bile_plot", height = 300)),
            column(6, plotOutput("rratio_plot", height = 300))
          ),
          hr(),
          fluidRow(
            column(6, plotOutput("innate_plot", height = 280)),
            column(6, plotOutput("adaptive_plot", height = 280))
          )
        ),

        # ---------------------------------------------------------------- 8
        tabPanel("⑧ 해독제 시간창",
          br(),
          p("NAC 개시 시각을 훑습니다. '8시간 규칙'은 모델 안에 상수로 들어 ",
            "있지 않습니다 — 궤적이 분리선을 넘는 시각으로 계산되어 나오며, ",
            "따라서 용량과 숙주 상태에 따라 이동합니다."),
          actionButton("run_nac", "NAC 시간창 계산 (약 40회 시뮬레이션)",
                       class = "btn-warning"),
          br(), br(),
          plotOutput("nacwin_plot", height = 380),
          DTOutput("nacwin_tab")
        ),

        # ---------------------------------------------------------------- 9
        tabPanel("⑨ 시나리오 비교",
          br(),
          checkboxGroupInput("scen_pick", "비교할 시나리오",
            choices = c("S2 150 mg/kg 무치료"            = "S2",
                        "S3 250 mg/kg 무치료"            = "S3",
                        "S4 350 mg/kg 무치료"            = "S4",
                        "S5 350 mg/kg + NAC 8h"          = "S5",
                        "S6 350 mg/kg + NAC 16h"         = "S6",
                        "S8 150 mg/kg 음주·금식 숙주"    = "S8",
                        "S10 BSEP 억제제 28일"           = "S10",
                        "S11 HLA+ 특이체질 56일"         = "S11",
                        "S12 + 면역관문억제제"           = "S12",
                        "S13 + 스테로이드"               = "S13"),
            selected = c("S3", "S4", "S5", "S8"), inline = TRUE),
          actionButton("run_scen", "선택 시나리오 실행", class = "btn-warning"),
          br(), br(),
          plotOutput("scen_plot", height = 420),
          DTOutput("scen_tab")
        ),

        # --------------------------------------------------------------- 10
        tabPanel("⑩ 속도 대 용량",
          br(),
          p("총 용량을 고정하고 도달 속도만 바꿉니다. 간손상이 용량의 문제라면 ",
            "아래 막대는 모두 같은 높이여야 합니다."),
          actionButton("run_rate", "속도 실험 실행", class = "btn-warning"),
          br(), br(),
          plotOutput("rate_plot", height = 360),
          DTOutput("rate_tab")
        )
      )
    )
  )
)

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  arche_params <- reactive({
    switch(input$archetype,
      A = list(),
      B = DRUG_B,
      C = DRUG_C)
  })

  dose_times <- reactive({
    switch(input$regimen,
      single   = 0,
      q6h7d    = seq(0, 162, by = 6),
      q12h28d  = seq(0, 660, by = 12),
      q24h56d  = seq(0, 1320, by = 24),
      spread   = seq(0, input$spread_h, length.out = max(2, input$spread_h)))
  })

  run <- eventReactive(input$go, {
    n  <- length(dose_times())
    ki <- 10^input$logki
    args <- c(
      list(dose_mgkg  = input$dose * n,
           dose_times = dose_times(),
           nac_start  = if (input$use_nac) input$nac_t else NULL,
           tend       = input$tend,
           FCYP    = input$fcyp,
           FGSH    = input$fgsh,
           CYSBASE = input$cysbase,
           AGEF    = input$agef,
           FUNC    = input$func,
           HLA     = as.numeric(input$hla),
           ICI     = as.numeric(input$ici),
           STER    = as.numeric(input$use_ster),
           STER_T0 = if (input$use_ster) input$ster_t * 24 else 1e9),
      arche_params())
    # explicit BSEP slider overrides the archetype default
    if (input$logki < 9) args$KI_BSEP <- ki
    as.data.frame(do.call(sim_dili, args))
  }, ignoreNULL = FALSE)

  smry <- reactive({
    d <- run()
    list(alt = max(d$ALT), ast = max(d$AST), alp = max(d$ALP),
         tbil = max(d$TBIL), inr = max(d$INR),
         R = max(d$ALT)/40 / (max(d$ALP)/120),
         hy = any(d$HYLAW > 0),
         gsh = min(d$GSHPCT), jnk = max(d$JNK), atp = min(d$ATP),
         lost = max(d$LOST), ba = max(d$BAP), mir = max(d$MIR),
         gmin = min(d$REDOXG), fb = max(d$FBIOACT))
  })

  # ----------------------------------------------------------------- tab 1
  output$verdict <- renderUI({
    s <- smry()
    patt <- if (s$R >= 5) "간세포형 (hepatocellular)" else
            if (s$R <= 2) "담즙정체형 (cholestatic)" else "혼합형 (mixed)"
    sev <- if (s$lost < 0.01) list("적응 (adaptation) — 임상 손상 없음", "#1E7B34")
      else if (s$lost < 0.10) list("경증 손상 — 회복 예상", "#B8860B")
      else if (s$lost < 0.40) list("중등도 손상", "#E07B00")
      else if (s$lost < 0.70) list("중증 손상", "#B03A2E")
      else list("전격성 간부전 범위", "#7D1128")
    HTML(sprintf(
      "<div style='font-size:16px;font-weight:700;color:%s'>%s</div>
       <hr style='margin:8px 0'>
       <b>손상 패턴</b>: %s (R = %.2f)<br>
       <b>Hy's Law</b>: <span style='color:%s;font-weight:700'>%s</span><br>
       <b>소실 간세포 질량</b>: %.1f%%<br>
       <b>JNK 고리</b>: %s (peak %.3f)<br>
       <b>최저 산화환원 여력 g</b>: %.2f",
      sev[[2]], sev[[1]], patt, s$R,
      if (s$hy) "#B71C1C" else "#1E7B34", if (s$hy) "충족 (YES)" else "미충족",
      100*s$lost,
      if (s$jnk > 0.3) "점화됨 (latched ON)" else "점화되지 않음", s$jnk,
      s$gmin))
  })

  output$keytab <- renderTable({
    s <- smry()
    data.frame(
      지표 = c("최고 ALT (U/L)", "최고 AST (U/L)", "최고 ALP (U/L)",
               "최고 총빌리루빈 (mg/dL)", "최고 INR", "최저 GSH (% 기저)",
               "최저 ATP (배수)", "최고 혈청 담즙산 (µM)", "최고 miR-122 (배)",
               "최대 생체활성화 분율"),
      값 = c(sprintf("%.0f", s$alt), sprintf("%.0f", s$ast),
             sprintf("%.0f", s$alp), sprintf("%.2f", s$tbil),
             sprintf("%.2f", s$inr), sprintf("%.1f", s$gsh),
             sprintf("%.3f", s$atp), sprintf("%.1f", s$ba),
             sprintf("%.0f", s$mir), sprintf("%.1f%%", 100*s$fb)))
  }, striped = TRUE, width = "100%")

  output$narrative <- renderUI({
    s <- smry(); d <- run()
    steps <- c()
    steps <- c(steps, sprintf("최대 %.0f%%의 간대사가 반응성 대사체로 갔습니다.",
                              100*s$fb))
    steps <- c(steps, sprintf("GSH는 기저의 %.0f%%까지 떨어졌습니다.", s$gsh))
    steps <- c(steps, if (s$jnk > 0.3)
      sprintf("산화환원 여력이 g = %.2f 까지 내려가 JNK 고리가 점화되었습니다.", s$gmin)
      else "산화환원 여력이 유지되어 JNK 고리는 OFF 상태에 머물렀습니다.")
    steps <- c(steps, if (s$ba > 20)
      sprintf("담즙산이 혈청 %.0f µM 까지 올라 담관 손상이 발생했습니다.", s$ba)
      else "담즙산 항상성은 유지되었습니다.")
    steps <- c(steps, sprintf("최종적으로 간세포 질량의 %.1f%%가 소실되었습니다.",
                              100*s$lost))
    HTML(paste0("<ol style='padding-left:18px;margin-bottom:0'>",
                paste0("<li>", steps, "</li>", collapse = ""), "</ol>"))
  })

  output$fluxplot <- renderPlot({
    d <- run()
    p <- as.list(param(mod))
    # reconstruct the two fluxes from captured states
    prod <- p$VMAX_CYP * d$CU/(p$KM_CYP + d$CU) * input$fcyp / p$VLIV / 1000
    disp <- rep(p$VMAX_GSH * input$fgsh, nrow(d)) * (d$NRF2)
    df <- data.frame(time = d$time,
                     `반응성 대사체 생성 (mM/h)` = prod,
                     `GSH 재합성 능력 (mM/h)`    = disp,
                     check.names = FALSE) %>%
      pivot_longer(-time)
    ggplot(df, aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c("#1E7B34", "#B03A2E")) +
      labs(x = "시간 (h)", y = "flux (mM/h)", colour = NULL) + THEME
  })

  # ----------------------------------------------------------------- tab 2
  output$pk_conc <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_line(aes(y = CP, colour = "혈장 (plasma)"), linewidth = 1) +
      geom_line(aes(y = CU, colour = "간 내 유리 (unbound liver)"),
                linewidth = 1) +
      scale_colour_manual(values = PAL[1:2]) +
      labs(title = "약물 농도", x = "시간 (h)", y = "µM", colour = NULL) + THEME
  })

  output$pk_nomogram <- renderPlot({
    d <- run()
    t <- seq(4, 24, by = 0.5)
    nomo <- data.frame(time = t, y = 150 * exp(-log(2)/4 * (t - 4)))
    ggplot() +
      geom_line(data = nomo, aes(time, y), colour = "#B03A2E",
                linetype = 2, linewidth = 1) +
      geom_line(data = subset(d, time <= 24), aes(time, CP_UGML),
                colour = PAL[1], linewidth = 1) +
      annotate("text", x = 16, y = 60, colour = "#B03A2E",
               label = "Rumack-Matthew 치료선\n(4 h, 150 µg/mL)") +
      scale_y_log10() +
      labs(title = "Rumack-Matthew 노모그램 (단회 급성 섭취에만 유효)",
           x = "섭취 후 시간 (h)", y = "혈장 농도 (µg/mL, log)") + THEME
  })

  output$pk_partition <- renderPlot({
    d <- run(); p <- as.list(param(mod))
    ap <- arche_params()
    vs <- if (!is.null(ap$VMAX_SULT)) ap$VMAX_SULT else p$VMAX_SULT
    vu <- if (!is.null(ap$VMAX_UGT))  ap$VMAX_UGT  else p$VMAX_UGT
    vc <- if (!is.null(ap$VMAX_CYP))  ap$VMAX_CYP  else p$VMAX_CYP
    df <- data.frame(
      time = d$time,
      `황산포합 (SULT)`      = vs*d$CU/(p$KM_SULT+d$CU)*d$PAPS,
      `글루쿠론산포합 (UGT)` = vu*d$CU/(p$KM_UGT +d$CU)*d$UDPGA,
      `CYP 생체활성화`       = vc*d$CU/(p$KM_CYP +d$CU)*input$fcyp,
      check.names = FALSE) %>% pivot_longer(-time)
    ggplot(df, aes(time, value, fill = name)) +
      geom_area(position = "fill") +
      scale_fill_manual(values = c("#B03A2E", PAL[3], PAL[4])) +
      scale_y_continuous(labels = scales::percent) +
      labs(x = "시간 (h)", y = "대사 flux 분율", fill = NULL) + THEME
  })

  # ----------------------------------------------------------------- tab 3
  output$gsh_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_hline(yintercept = 30, linetype = 3, colour = "grey40") +
      geom_line(aes(y = GSHPCT, colour = "GSH (% 기저)"), linewidth = 1) +
      geom_line(aes(y = 100*CYS, colour = "시스테인 가용성 (%)"), linewidth = 1) +
      scale_colour_manual(values = PAL[c(3,4)]) +
      labs(title = "글루타티온과 그 율속 기질", x = "시간 (h)",
           y = "% of baseline", colour = NULL) + THEME
  })

  output$rm_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_line(aes(y = RM*1000, colour = "반응성 대사체 (nM)"), linewidth = 1) +
      geom_line(aes(y = ADD, colour = "미토 단백 부가체 (au)"), linewidth = 1) +
      scale_colour_manual(values = c("#7D1128", "#D62728")) +
      labs(title = "반응성 대사체와 공유결합 부가체",
           x = "시간 (h)", y = NULL, colour = NULL) + THEME
  })

  output$gate_plot <- renderPlot({
    p <- as.list(param(mod))
    g <- 10^seq(-2, 1, length.out = 300)
    frac <- p$KBIND / (p$KBIND + p$KGST*g + p$KNQO)
    d <- run()
    ggplot(data.frame(g, frac), aes(g, 100*frac)) +
      geom_line(linewidth = 1.1, colour = "#B03A2E") +
      geom_vline(xintercept = min(d$GSH), linetype = 2, colour = PAL[1]) +
      annotate("text", x = min(d$GSH), y = 40, hjust = -0.1, colour = PAL[1],
               label = sprintf("이 실행의 최저 GSH = %.2f mM", min(d$GSH))) +
      scale_x_log10() +
      labs(title = "GSH 게이트: 단백질에 공유결합하는 반응성 대사체의 분율",
           x = "간 GSH (mM, log)", y = "단백 결합 분율 (%)") + THEME
  })

  # ----------------------------------------------------------------- tab 4
  output$mito_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_hline(yintercept = 0.45, linetype = 3, colour = "grey40") +
      geom_line(aes(y = MITO, colour = "미토콘드리아 기능"), linewidth = 1) +
      geom_line(aes(y = ATP,  colour = "ATP"), linewidth = 1) +
      geom_line(aes(y = ROS/5, colour = "ROS / 5"), linewidth = 1) +
      scale_colour_manual(values = PAL[c(1,3,2)]) +
      labs(title = "생체에너지학 (점선 = MPT ATP 게이트 반개방점)",
           x = "시간 (h)", y = "배수 (fold)", colour = NULL) + THEME
  })

  output$jnk_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_line(aes(y = JNK,    colour = "활성 JNK"), linewidth = 1.1) +
      geom_line(aes(y = REDOXG, colour = "산화환원 여력 g"), linewidth = 1) +
      geom_line(aes(y = NRF2/3, colour = "NRF2 / 3"), linewidth = 1) +
      scale_colour_manual(values = c("#B03A2E", PAL[1], PAL[3])) +
      labs(title = "JNK 고리와 그것을 억제하는 여력",
           x = "시간 (h)", y = NULL, colour = NULL) + THEME
  })

  output$bistab_plot <- renderPlot({
    bb <- analysis_bistability(seq(0.6, 1.6, by = 0.02))
    d  <- run()
    ggplot(bb, aes(g)) +
      geom_line(aes(y = on_state, colour = "안정 ON 상태"), linewidth = 1.2,
                na.rm = TRUE) +
      geom_line(aes(y = separatrix, colour = "분리선 (불안정)"),
                linewidth = 1.2, linetype = 2, na.rm = TRUE) +
      geom_point(data = data.frame(g = d$REDOXG, J = d$JNK),
                 aes(g, J), colour = "grey35", size = 0.35, alpha = 0.5) +
      scale_colour_manual(values = c("#B03A2E", "#7D1128")) +
      labs(title = "JNK–Sab 고리의 분기 다이어그램과 이번 실행의 궤적(회색)",
           x = "산화환원 여력 g", y = "활성 JNK", colour = NULL) + THEME
  })

  # ----------------------------------------------------------------- tab 5
  output$mass_plot <- renderPlot({
    d <- run()
    df <- data.frame(time = d$time,
                     `생존 간세포 (HEP)`     = d$HEP,
                     `스트레스 간세포 (HEPS)` = d$HEPS,
                     `누적 괴사 (NECR)`      = d$NECR,
                     check.names = FALSE) %>% pivot_longer(-time)
    ggplot(df, aes(time, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c(PAL[3], PAL[4], "#B03A2E")) +
      labs(title = "간세포 생명주기", x = "시간 (h)",
           y = "정상 간 대비 분율", colour = NULL) + THEME
  })

  output$damp_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_line(aes(y = DAMP, colour = "DAMP"), linewidth = 1) +
      geom_line(aes(y = KC,   colour = "활성 쿠퍼세포"), linewidth = 1) +
      scale_colour_manual(values = PAL[c(2,5)]) +
      labs(title = "무균성 염증의 개시", x = "시간 (h)", y = NULL,
           colour = NULL) + THEME
  })

  output$regen_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_line(aes(y = TNF,  colour = "TNF-α (재생 프라이밍 + 손상)"),
                linewidth = 1) +
      geom_line(aes(y = IL10, colour = "IL-10 (항염)"), linewidth = 1) +
      scale_colour_manual(values = c("#D62728", PAL[3])) +
      labs(title = "TNF의 양날: 손상 확대이자 재생 신호",
           x = "시간 (h)", y = "au", colour = NULL) + THEME
  })

  # ----------------------------------------------------------------- tab 6
  output$lab_enz <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_hline(yintercept = 120, linetype = 3, colour = "#B71C1C") +
      geom_line(aes(y = ALT, colour = "ALT"), linewidth = 1.1) +
      geom_line(aes(y = AST, colour = "AST"), linewidth = 1) +
      geom_line(aes(y = ALP, colour = "ALP"), linewidth = 1) +
      geom_line(aes(y = MIR*10, colour = "miR-122 × 10"), linewidth = 1) +
      scale_y_log10() +
      scale_colour_manual(values = PAL[c(2,4,5,1)]) +
      labs(title = "효소 및 기전 바이오마커 (점선 = ALT 3×ULN)",
           x = "시간 (h)", y = "U/L 또는 배수 (log)", colour = NULL) + THEME
  })

  output$lab_fn <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_hline(yintercept = 2.4, linetype = 3, colour = "#B71C1C") +
      geom_line(aes(y = TBIL, colour = "총 빌리루빈 (mg/dL)"), linewidth = 1.1) +
      geom_line(aes(y = INR,  colour = "INR"), linewidth = 1.1) +
      scale_colour_manual(values = c("#B8860B", "#7D1128")) +
      labs(title = "합성·배설 기능 (점선 = TBIL 2×ULN)",
           x = "시간 (h)", y = NULL, colour = NULL) + THEME
  })

  output$hylaw_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(ALT, TBIL)) +
      annotate("rect", xmin = 120, xmax = Inf, ymin = 2.4, ymax = Inf,
               fill = "#B71C1C", alpha = 0.10) +
      annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
               colour = "#B71C1C", fontface = 2, label = "Hy's Law 구역") +
      geom_path(aes(colour = time), linewidth = 1.1) +
      geom_vline(xintercept = 120, linetype = 3) +
      geom_hline(yintercept = 2.4, linetype = 3) +
      scale_x_log10() +
      scale_colour_viridis_c(name = "시간 (h)") +
      labs(title = "ALT–빌리루빈 궤적 (속도 축 × 예비력 축)",
           x = "ALT (U/L, log)", y = "총 빌리루빈 (mg/dL)") + THEME
  })

  # ----------------------------------------------------------------- tab 7
  output$bile_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_hline(yintercept = 100, linetype = 3, colour = "grey40") +
      geom_line(aes(y = BAH, colour = "간세포 내 담즙산 (µM)"), linewidth = 1.1) +
      geom_line(aes(y = BAP, colour = "혈청 총 담즙산 (µM)"), linewidth = 1.1) +
      scale_colour_manual(values = c("#B8860B", PAL[1])) +
      labs(title = "담즙산 항상성 (점선 = 세포독성 임계 100 µM)",
           x = "시간 (h)", y = "µM", colour = NULL) + THEME
  })

  output$rratio_plot <- renderPlot({
    d <- subset(run(), ALT > 0 & ALP > 0)
    ggplot(d, aes(time, RRATIO)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 5, ymax = Inf,
               fill = PAL[2], alpha = 0.08) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 2,
               fill = PAL[4], alpha = 0.10) +
      geom_line(linewidth = 1.1, colour = "grey20") +
      scale_y_log10() +
      labs(title = "R-비 (위 = 간세포형, 아래 = 담즙정체형) — 계산으로 나온 결과",
           x = "시간 (h)", y = "R = (ALT/ULN)/(ALP/ULN), log") + THEME
  })

  output$innate_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_line(aes(y = KC,   colour = "쿠퍼세포"), linewidth = 1) +
      geom_line(aes(y = TNF,  colour = "TNF-α"), linewidth = 1) +
      geom_line(aes(y = IL10, colour = "IL-10"), linewidth = 1) +
      scale_colour_manual(values = PAL[c(5,2,3)]) +
      labs(title = "선천면역", x = "시간 (h)", y = "au", colour = NULL) + THEME
  })

  output$adaptive_plot <- renderPlot({
    d <- run()
    ggplot(d, aes(time)) +
      geom_line(aes(y = TCELL, colour = "약물특이 효과 T세포"), linewidth = 1.1) +
      geom_line(aes(y = TREG,  colour = "관용 (Treg/PD-1)"), linewidth = 1.1) +
      scale_colour_manual(values = c("#B03A2E", PAL[3])) +
      labs(title = "적응면역 — 특이체질 DILI를 결정하는 것은 관용이다",
           x = "시간 (h)", y = "분율", colour = NULL) + THEME
  })

  # ----------------------------------------------------------------- tab 8
  nacwin <- eventReactive(input$run_nac, {
    withProgress(message = "NAC 시간창 계산 중...", {
      analysis_nac_window()
    })
  })

  output$nacwin_plot <- renderPlot({
    d <- nacwin()
    d$x <- ifelse(is.na(d$nac_start), 60, d$nac_start)
    ggplot(d, aes(x, 100*lost_mass, colour = host)) +
      geom_line(linewidth = 1.1) + geom_point(size = 2) +
      scale_colour_manual(values = PAL[1:3]) +
      scale_x_continuous(breaks = c(0,8,16,24,32,48,60),
                         labels = c("0","8","16","24","32","48","무치료")) +
      labs(title = "NAC 개시 시각에 따른 간세포 소실 — 절벽의 위치는 이동한다",
           x = "NAC 개시 시각 (h)", y = "소실 간세포 질량 (%)",
           colour = NULL) + THEME
  })

  output$nacwin_tab <- renderDT({
    datatable(nacwin(), options = list(pageLength = 15, scrollX = TRUE)) %>%
      formatRound(columns = which(sapply(nacwin(), is.numeric)), digits = 3)
  })

  # ----------------------------------------------------------------- tab 9
  scen_res <- eventReactive(input$run_scen, {
    sc <- scenarios()
    pick <- input$scen_pick
    withProgress(message = "시나리오 실행 중...", {
      do.call(rbind, lapply(pick, function(k) {
        r <- summarise_run(do.call(sim_dili, sc[[k]]$args))
        cbind(id = k, label = sc[[k]]$label, r)
      }))
    })
  })

  output$scen_plot <- renderPlot({
    d <- scen_res()
    df <- d %>%
      select(label, `최고 ALT` = peak_ALT, `최고 TBIL x1000` = peak_TBIL,
             `소실 질량 %` = lost_mass) %>%
      mutate(`최고 TBIL x1000` = `최고 TBIL x1000`*1000,
             `소실 질량 %` = `소실 질량 %`*100) %>%
      pivot_longer(-label)
    ggplot(df, aes(reorder(label, value), value, fill = name)) +
      geom_col(position = "dodge") + coord_flip() +
      scale_fill_manual(values = PAL[c(2,4,1)]) +
      labs(x = NULL, y = NULL, fill = NULL) + THEME
  })

  output$scen_tab <- renderDT({
    datatable(scen_res(), options = list(pageLength = 12, scrollX = TRUE)) %>%
      formatRound(columns = which(sapply(scen_res(), is.numeric)), digits = 3)
  })

  # ---------------------------------------------------------------- tab 10
  rate_res <- eventReactive(input$run_rate, {
    withProgress(message = "속도 실험 실행 중...", {
      analysis_rate_vs_dose(dose = input$dose)
    })
  })

  output$rate_plot <- renderPlot({
    d <- rate_res()
    ggplot(d, aes(factor(spread_h), 100*lost_mass)) +
      geom_col(fill = "#B03A2E") +
      geom_text(aes(label = sprintf("%.1f%%", 100*lost_mass)), vjust = -0.4) +
      labs(title = sprintf("총 %g mg/kg 를 서로 다른 속도로 투여했을 때",
                           input$dose),
           x = "총 섭취 시간 (h; 0 = 단회 볼루스)",
           y = "소실 간세포 질량 (%)") + THEME
  })

  output$rate_tab <- renderDT({
    datatable(rate_res(), options = list(pageLength = 10, scrollX = TRUE)) %>%
      formatRound(columns = which(sapply(rate_res(), is.numeric)), digits = 3)
  })
}

shinyApp(ui, server)
