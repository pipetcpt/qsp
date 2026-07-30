## =============================================================================
##  nb_shiny_app.R — High-Risk Neuroblastoma QSP explorer
##
##  Nine tabs, organised so that the model's two structural claims are the first
##  things a user can interrogate:
##
##    1  환자·요법 (Patient & regimen)      who is being treated, and with what
##    2  약물동태 (Pharmacokinetics)        four drug classes, four PK shapes
##    3  두 Fc 팔 (The two Fc arms)         the core: ADCC n=1 vs CDC n=2
##    4  전달 장벽 (Delivery barriers)      permeability vs dose, side by side
##    5  종양 부담 (Tumour burden & PD)     four pools, four drug sensitivities
##    6  임상 엔드포인트 (Endpoints)         relapse, pain, ANC/PLT, hearing
##    7  시나리오 비교 (Scenarios)           the nine derived comparisons
##    8  MIBG 선량 (Radiopharmaceutical)     transporter competition + dosimetry
##    9  가상 코호트 (Virtual cohort / EFS)  ANBL0032 and HR-NBL1 reproduction
##
##  Run:  shiny::runApp("nb_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  The model itself lives in nb_mrgsolve_model.R and is sourced below.
## =============================================================================

suppressPackageStartupMessages({
  library(shiny); library(mrgsolve); library(dplyr); library(tidyr)
  library(ggplot2); library(DT)
})

source("nb_mrgsolve_model.R")

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

## a colour per target compartment — used consistently in every tab so that
## "marrow", "nerve" and "solid tumour" always mean the same colour
CMP <- c("solid tumour" = "#C0392B", "bone marrow" = "#1E8449",
         "DRG / nerve" = "#CA6F1E", "plasma" = "#2874A6")

# =============================================================================
ui <- fluidPage(
  titlePanel("고위험 신경교세포종(신경블라스토마) QSP 모델 — High-Risk Neuroblastoma"),
  tags$p(style = "color:#555;margin-top:-10px",
         HTML("두 개의 Fc 이펙터 팔(ADCC ∝ 결합IgG¹ · CDC ∝ 결합IgG²)과 세 개의
               투과성(고형종양 : 골수 : 배근신경절 ≈ 1 : 15 : 200)이 이 모델의 전부입니다.
               <b>The whole model is two Fc arms and three permeabilities.</b>")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      sliderInput("BSA", "체표면적 BSA (m²)", 0.4, 1.6, 0.8, 0.05),
      sliderInput("WT", "체중 (kg)", 8, 60, 20, 1),
      sliderInput("tp0", "초기 고형 부담 (g)", 1, 400, 100, 1),
      sliderInput("tm0", "초기 골수 MRD (10⁹ cells)", 0, 30, 5, 0.5),
      checkboxInput("MYCN", "MYCN 증폭 (amplified)", TRUE),
      selectInput("FCG", "FcγR 유전형 (FCGR3A/2A)",
                  c("저친화 F/F (1.00)" = 1.00, "이형접합 V/F (1.35)" = 1.35,
                    "고친화 V/V (1.80)" = 1.80), selected = 1.35),
      sliderInput("GD2", "GD2 밀도 (10⁶ /cell)", 1, 20, 8, 0.5),
      sliderInput("FMES", "MES(항원저발현) 분율", 0, 0.6, 0.05, 0.01),

      hr(), h4("요법 (Regimen)"),
      checkboxInput("immuno", "항-GD2 항체 (dinutuximab)", TRUE),
      selectInput("fc", "Fc 변이체",
                  c("ch14.18 (C1q 정상)" = 1.00,
                    "hu14.18K322A (C1q ×0.1)" = 0.10), selected = 1.00),
      sliderInput("abdose", "항체 용량 (mg/m²/일)", 0.5, 70, 17.5, 0.5),
      sliderInput("abhrs", "주입 시간 (h)", 0.5, 24, 10, 0.5),
      checkboxInput("use_gm", "GM-CSF (cycles 1,3,5)", TRUE),
      checkboxInput("use_il2", "IL-2 (cycles 2,4 — COG)", TRUE),
      checkboxInput("use_ra", "Isotretinoin", TRUE),
      selectInput("frel", "Isotretinoin 투여방법",
                  c("캡슐 통째로 복용 (F 1.00)" = 1.00,
                    "캡슐 개방·음식에 섞음 (F 0.60)" = 0.60), selected = 1.00),
      checkboxInput("ra_conc", "Isotretinoin을 유도요법과 동시 투여", FALSE),
      selectInput("hdct", "고용량 공고요법", c("CEM", "BuMel")),
      checkboxInput("tandem", "Tandem ASCT (ANBL0532)", FALSE),
      sliderInput("ests", "Sodium thiosulfate 내이 보호", 0, 0.75, 0, 0.05),
      sliderInput("ftumsts", "  …종양 보호(부작용) 가정", 0, 0.30, 0, 0.05),

      hr(), h4("전달 (Delivery)"),
      sliderInput("psgtu", "고형종양 투과성 PSG_TU ×", 0.1, 25, 1, 0.1),
      sliderInput("vifp", "IFP 반감 부피 VIFP (L)", 0.005, 0.30, 0.033, 0.005),

      hr(),
      sliderInput("tend", "시뮬레이션 기간 (일)", 400, 1500, 1200, 50),
      actionButton("go", "실행 (Run)", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ------------------------------------------------- 1. patient & regimen
        tabPanel("1 · 환자·요법",
          h4("적용된 요법 일정 (dosing schedule as built)"),
          DTOutput("schedule"),
          h4("파생 수량 (derived quantities at these settings)"),
          tableOutput("derived"),
          tags$small(HTML(
            "GD2 항원량은 <b>g당 몰수</b>로 계산됩니다: 8×10⁶ copies/cell →
             13.3 nmol/g. 한 코스의 dinutuximab은 373 nmol이므로 항체는
             희소자원이 아니며, 결정하는 것은 <b>혈관 투과성</b>입니다."))),

        ## ------------------------------------------------------------- 2. PK
        tabPanel("2 · 약물동태 (PK)",
          fluidRow(
            column(6, plotOutput("pk_ab", height = 260)),
            column(6, plotOutput("pk_ra", height = 260))),
          fluidRow(
            column(6, plotOutput("pk_ct", height = 260)),
            column(6, plotOutput("pk_cyt", height = 260))),
          tags$small("항체는 혈장·고형종양간질·골수간질·배근신경절간질 네 곳의
                      농도로 표시됩니다. 세 조직의 차이는 투과성뿐입니다.")),

        ## --------------------------------------------------- 3. the two Fc arms
        tabPanel("3 · 두 Fc 팔",
          fluidRow(
            column(6, plotOutput("bpc", height = 280)),
            column(6, plotOutput("arms", height = 280))),
          h4("용량-반응 기하 (dose-response geometry)"),
          plotOutput("ti_curve", height = 300),
          tags$small(HTML(
            "ADCC는 결합 IgG/세포에 <b>1차</b>로 비례하고(하나의 Fc가 하나의
             FcγR을 결합), CDC는 <b>2차</b>로 비례합니다(C1q는 인접한 <b>두</b>
             개의 Fc를 걸쳐야 함). 두 팔이 서로 다른 구획에 있고 그 구획의
             투과성이 다르기 때문에, 최적 용량은 <b>병변의 위치</b>에 따라
             달라집니다."))),

        ## --------------------------------------------- 4. delivery vs dose
        tabPanel("4 · 전달 장벽",
          plotOutput("perm_vs_dose", height = 330),
          h4("주입된 항체의 운명 (fate of the injected antibody)"),
          tableOutput("fate"),
          tags$small("모델에서 용량과 투과성은 곱으로만 등장하므로 고형종양 팔에
                      대해 완전히 교환 가능합니다. 그러나 투과성을 올리는 것은
                      통증 구획을 건드리지 않습니다.")),

        ## ------------------------------------------------------ 5. tumour PD
        tabPanel("5 · 종양 부담·PD",
          plotOutput("burden", height = 300),
          fluidRow(
            column(6, plotOutput("pools", height = 260)),
            column(6, plotOutput("killrates", height = 260))),
          tags$small("네 개의 종양 풀은 서로 다른 약물감수성을 가집니다:
                      TP(세포독성+ADCC+방사선), TQ(G0 — 세포독성 저항),
                      TD(분화 — 세포독성 무반응), TM(골수 — ADCC가 실제로
                      작동하는 곳).")),

        ## ------------------------------------------------------- 6. endpoints
        tabPanel("6 · 임상 엔드포인트",
          fluidRow(
            column(4, wellPanel(h4("재발 시점"), textOutput("t_relapse"))),
            column(4, wellPanel(h4("최저 부담 (nadir)"), textOutput("t_nadir"))),
            column(4, wellPanel(h4("INRC 반응"), textOutput("t_resp")))),
          fluidRow(
            column(6, plotOutput("pain", height = 260)),
            column(6, plotOutput("heme", height = 260))),
          fluidRow(
            column(6, plotOutput("oto", height = 240)),
            column(6, tableOutput("tox_table")))),

        ## ------------------------------------------------------- 7. scenarios
        tabPanel("7 · 시나리오 비교",
          checkboxGroupInput("scen", "비교할 축 (axes to run)",
            choices = c("항체 용량 (Fc 두 팔)" = "dose",
                        "투과성 vs 용량" = "perm",
                        "K322A Fc 변이" = "k322a",
                        "종양 부담 의존성" = "burden",
                        "IL-2 / GM-CSF" = "cyto",
                        "Isotretinoin 스케줄" = "retinoid",
                        "공고요법 강도" = "consol",
                        "STS 이독성 상충" = "sts"),
            selected = c("dose", "k322a"), inline = TRUE),
          actionButton("runscen", "시나리오 실행", class = "btn-warning"),
          hr(), uiOutput("scen_out")),

        ## ------------------------------------------------------------ 8. MIBG
        tabPanel("8 · MIBG 선량",
          fluidRow(
            column(4, sliderInput("mibg_mci", "¹³¹I-MIBG (mCi/kg)", 0, 24, 18, 1)),
            column(4, selectInput("mibg_sa", "비방사능 (specific activity)",
                     c("무담체 NCA (1.1e5 MBq/µmol)" = 1.10e5,
                       "담체첨가 ×10 (1.1e4)" = 1.10e4,
                       "담체첨가 ×100 (1.1e3)" = 1.10e3))),
            column(4, selectInput("netx", "NET 차단 약물 (중단 실패)",
                     c("없음 (1.00)" = 1.00, "삼환계 항우울제 (0.55)" = 0.55,
                       "labetalol (0.35)" = 0.35)))),
          fluidRow(
            column(4, sliderInput("kexcr", "전신 생물학적 제거율 (1/d)",
                                  0.08, 0.40, 0.19, 0.01)),
            column(4, checkboxInput("ki", "요오드화칼륨(KI) 갑상선 차단", TRUE)),
            column(4, br(), actionButton("runmibg", "MIBG 실행",
                                         class = "btn-info"))),
          hr(),
          fluidRow(
            column(6, plotOutput("mibg_dose", height = 280)),
            column(6, plotOutput("mibg_act", height = 280))),
          tableOutput("mibg_tab"),
          tags$small(HTML(
            "전신선량 계수는 <b>보정된</b> 값(0.20 mGy/MBq)이므로 예측이 아닙니다.
             예측되는 것은, 생물학적 제거율이 변할 때 2 Gy 전신선량 목표를 맞추기
             위해 필요한 <b>투여 활성도</b>가 어떻게 변하는지입니다."))),

        ## ---------------------------------------------------------- 9. cohort
        tabPanel("9 · 가상 코호트 (EFS)",
          fluidRow(
            column(3, numericInput("ncoh", "환자 수", 20, 5, 200, 5)),
            column(3, numericInput("seed", "seed", 20260730, 1, 99999999, 1)),
            column(3, checkboxGroupInput("coh_arms", "비교 arm",
                     c("Isotretinoin 단독" = "ctl", "COG 전체" = "cog",
                       "IL-2 제외 (SIOPEN)" = "noil2"),
                     selected = c("ctl", "cog"))),
            column(3, br(), actionButton("runcoh", "코호트 실행",
                                         class = "btn-danger"))),
          hr(),
          plotOutput("efs_curve", height = 340),
          tableOutput("efs_tab"),
          tags$small("보정 목표: ANBL0032 2년 EFS 66% (면역요법) vs 46% (대조).
                      작은 코호트에서는 신뢰구간이 넓으므로 표에 함께 표시합니다."))
      )
    )
  )
)

# =============================================================================
server <- function(input, output, session) {

  ## ------------------------------------------------ parameters from the sidebar
  pset <- reactive({
    list(BSA = input$BSA, WT = input$WT,
         MYCN = as.numeric(input$MYCN),
         FCG = as.numeric(input$FCG),
         GD2DENS = input$GD2 * 1e6,
         FMES = input$FMES,
         C1QEFF = as.numeric(input$fc),
         PSG_TU = 0.020 * input$psgtu,
         VIFP = input$vifp,
         ESTS = input$ests, FTUMSTS = input$ftumsts,
         IMMUNO = as.numeric(input$immuno),
         USE_RA = as.numeric(input$use_ra))
  })

  regimen <- reactive({
    nb_regimen(immuno = input$immuno, use_il2 = input$use_il2,
               use_gm = input$use_gm, use_ra = input$use_ra,
               hdct = input$hdct, tandem = input$tandem,
               ra_concurrent = input$ra_conc, ab_mgm2 = input$abdose,
               ab_hours = input$abhrs, BSA = input$BSA,
               FREL = as.numeric(input$frel))
  })

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    withProgress(message = "적분 중 (integrating)…", {
      nb_run(e = regimen(), tp0 = input$tp0, tm0 = input$tm0,
             end = input$tend, delta = 0.25, param = pset())
    })
  })

  ## ------------------------------------------------------------ 1. schedule
  output$schedule <- renderDT({
    regimen() %>% as_tibble() %>%
      mutate(cmt = as.character(cmt)) %>%
      arrange(time) %>% head(200) %>%
      datatable(options = list(pageLength = 8, dom = "tp"), rownames = FALSE)
  })

  output$derived <- renderTable({
    p <- pset()
    ab_nmol_day <- input$abdose * input$BSA / 1000 / MW_AB * 1e9
    ag_per_g <- 1e9 * p$GD2DENS / 6.022e14
    tibble(
      quantity = c("dinutuximab, nmol/일", "dinutuximab, nmol/코스(4일)",
                   "GD2 항원, nmol/g 종양",
                   "고형종양 총 항원 (nmol)",
                   "PS 곱 — 고형종양 : 골수 : 신경 (L/d)",
                   "ADCC50 / CDC50 (분자/세포)"),
      value = c(sprintf("%.1f", ab_nmol_day),
                sprintf("%.1f", 4 * ab_nmol_day),
                sprintf("%.2f", ag_per_g),
                sprintf("%.0f", ag_per_g * 1.25 * input$tp0),
                sprintf("%.4f : %.3f : %.3f",
                        p$PSG_TU * 1.25 * input$tp0 / 1000, 0.300 * 0.300,
                        4.000 * 0.0020),
                "2.0e4 / 5.0e4"))
  })

  ## ------------------------------------------------------------------- 2. PK
  output$pk_ab <- renderPlot({
    o <- sim()
    o %>% transmute(time,
                    `plasma` = CPn,
                    `solid tumour` = ABT / (0.35 * pmax(VTU, 1e-9)),
                    `bone marrow` = ABM / (0.35 * 0.300),
                    `DRG / nerve` = ABN / (0.35 * 0.0020)) %>%
      pivot_longer(-time) %>%
      filter(time >= 195, time <= 380) %>%
      ggplot(aes(time, pmax(value, 1e-4), colour = name)) +
      geom_line(linewidth = 0.6) + scale_y_log10() +
      scale_colour_manual(values = CMP, name = NULL) +
      labs(title = "항-GD2 항체: 네 구획의 유리 농도",
           subtitle = "면역요법 구간만 표시 · 차이는 투과성뿐",
           x = "일 (day)", y = "유리 항체 (nM, log)") + THEME
  })

  output$pk_ra <- renderPlot({
    sim() %>% filter(time >= 195) %>%
      ggplot(aes(time, CRAo)) + geom_line(colour = "#B9770E", linewidth = 0.6) +
      geom_hline(yintercept = 2, linetype = 2, colour = "grey40") +
      annotate("text", x = Inf, y = 2, label = "2 µM 목표", hjust = 1.05,
               vjust = -0.5, size = 3, colour = "grey30") +
      labs(title = "Isotretinoin (13-cis-RA)",
           subtitle = "CYP26A1 자가유도로 코스 내에서 노출이 스스로 감소",
           x = "일", y = "혈장 농도 (µM)") + THEME
  })

  output$pk_ct <- renderPlot({
    sim() %>% ggplot(aes(time, CCTo)) +
      geom_line(colour = "#5D6D7E", linewidth = 0.5) +
      labs(title = "복합 세포독성 약물 (composite cytotoxic)",
           subtitle = "유도 6주기 → 고용량 공고요법 → ASCT", x = "일",
           y = "농도 (mg/L)") + THEME
  })

  output$pk_cyt <- renderPlot({
    sim() %>% filter(time >= 195) %>%
      transmute(time, `NK (blood)` = NKB, Treg = TREG,
                `ANC ×50` = ANC * 50) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.6) +
      labs(title = "이펙터 공급 (effector supply)",
           subtitle = "IL-2는 Treg를 더 낮은 EC50에서 확장시킨다 (CD25 고친화)",
           x = "일", y = "cells/µL", colour = NULL) + THEME
  })

  ## --------------------------------------------------------- 3. the Fc arms
  output$bpc <- renderPlot({
    sim() %>% filter(time >= 195) %>%
      transmute(time, `solid tumour` = BPCTo, `bone marrow` = BPCMo,
                `DRG / nerve` = BPCNo) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, pmax(value, 1), colour = name)) +
      geom_line(linewidth = 0.7) + scale_y_log10() +
      geom_hline(yintercept = 2e4, linetype = 2, colour = "#1E8449") +
      geom_hline(yintercept = 5e4, linetype = 3, colour = "#CA6F1E") +
      scale_colour_manual(values = CMP, name = NULL) +
      labs(title = "결합 IgG / 세포 — 실제 PD 구동변수",
           subtitle = "파선 = ADCC50 (2×10⁴) · 점선 = CDC50 (5×10⁴)",
           x = "일", y = "분자/세포 (log)") + THEME
  })

  output$arms <- renderPlot({
    o <- sim() %>% filter(time >= 195)
    a <- 0.90 * (o$BPCTo / (o$BPCTo + 2e4))
    m <- 0.90 * 1.25 * (o$BPCMo / (o$BPCMo + 2e4))
    r <- (o$BPCNo / 5e4)^2
    tibble(time = o$time, `ADCC · solid tumour` = a, `ADCC · bone marrow` = m,
           `CDC · nociceptor` = r / (1 + r) * o$CPL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.7) +
      labs(title = "두 Fc 팔의 시간 경과", x = "일",
           subtitle = "CDC는 보체 소모(CPL)로 코스 내에서 감쇠한다",
           y = "상대 구동력", colour = NULL) + THEME
  })

  output$ti_curve <- renderPlot({
    ## algebraic, not simulated — instantaneous geometry at plasma steady state
    d <- 10^seq(log10(0.05), log10(300), length.out = 120)
    scale <- input$abdose
    bT <- 1.42e4 * d / scale; bM <- 3.01e6 * d / scale
    bMs <- 3.73e6 * bM / (bM + 3.73e6)            # marrow saturates on antigen
    bN <- 2.00e5 * d / scale; bNs <- 4.2e5 * bN / (bN + 4.2e5)
    adcc_t <- bT / (bT + 2e4); adcc_m <- bMs / (bMs + 2e4)
    r <- (bNs / 5e4)^2; cdc <- as.numeric(input$fc) * r / (1 + r)
    tibble(dose = d, `ADCC · solid tumour` = adcc_t,
           `ADCC · bone marrow` = adcc_m, `CDC · nociceptor` = cdc) %>%
      pivot_longer(-dose) %>%
      ggplot(aes(dose, value, colour = name)) + geom_line(linewidth = 0.8) +
      geom_vline(xintercept = 17.5, linetype = 2, colour = "grey40") +
      annotate("text", x = 17.5, y = 0.05, label = "임상 용량", angle = 90,
               vjust = -0.4, size = 3, colour = "grey30") +
      scale_x_log10() +
      labs(title = "용량-반응 기하: 골수 팔은 이미 포화, 고형종양 팔은 아직 선형",
           subtitle = "따라서 최적 용량은 병변의 위치에 따라 다르다",
           x = "dinutuximab (mg/m²/일, log)", y = "구동력 (0-1)",
           colour = NULL) + THEME
  })

  ## ------------------------------------------------ 4. permeability vs dose
  output$perm_vs_dose <- renderPlot({
    f <- c(0.1, 0.25, 0.5, 1, 2, 5, 10, 25)
    ## the model enters dose and permeability as a product, so one curve serves
    tibble(fold = f,
           `투과성 ×N (통증 불변)` = 1.42e4 * f,
           `용량 ×N (통증 증가)` = 1.42e4 * f) %>%
      pivot_longer(-fold) %>%
      mutate(drive = value / (value + 2e4)) %>%
      ggplot(aes(fold, drive, colour = name, linetype = name)) +
      geom_line(linewidth = 0.9) + geom_point(size = 1.8) + scale_x_log10() +
      labs(title = "고형종양 ADCC: 투과성과 용량은 완전히 교환 가능",
           subtitle = "곱으로만 등장하기 때문. 그러나 투과성은 통증 구획을 건드리지 않는다",
           x = "배수 (fold)", y = "ADCC 구동력", colour = NULL,
           linetype = NULL) + THEME
  })

  output$fate <- renderTable({
    o <- sim(); dt <- diff(o$time)
    intg <- function(x) sum(dt * head(x, -1))
    dose <- 4 * 5 * input$abdose * input$BSA / 1000 / MW_AB * 1e9
    tibble(compartment = c("주입 총량 (nmol)", "고형종양에 결합 (nmol)",
                           "골수에 결합 (nmol)", "배근신경절에 결합 (nmol)"),
           amount = c(dose, 0.030 * intg(o$BNDT), 0.030 * intg(o$BNDM),
                      0.030 * intg(o$BNDN)),
           `% of dose` = c(100, 100 * 0.030 * intg(o$BNDT) / dose,
                           100 * 0.030 * intg(o$BNDM) / dose,
                           100 * 0.030 * intg(o$BNDN) / dose))
  }, digits = 6)

  ## ---------------------------------------------------------- 5. tumour PD
  output$burden <- renderPlot({
    o <- sim()
    ggplot(o, aes(time, pmax(BURDEN, 1e-9))) +
      geom_line(linewidth = 0.7, colour = "#922B21") + scale_y_log10() +
      geom_hline(yintercept = 0.5, linetype = 2, colour = "grey40") +
      annotate("rect", xmin = 200, xmax = 368, ymin = 1e-9, ymax = Inf,
               alpha = 0.06, fill = "#1E8449") +
      annotate("text", x = 284, y = 1e-8, label = "면역요법 구간", size = 3,
               colour = "#145A32") +
      labs(title = "총 종양 부담", subtitle = "파선 = 재발 검출 임계 (0.5 g)",
           x = "일", y = "10⁹ cells (log)") + THEME
  })

  output$pools <- renderPlot({
    sim() %>% transmute(time, TP, TQ, TD, TM) %>% pivot_longer(-time) %>%
      ggplot(aes(time, pmax(value, 1e-9), colour = name)) +
      geom_line(linewidth = 0.6) + scale_y_log10() +
      labs(title = "네 개의 종양 풀", x = "일", y = "10⁹ cells (log)",
           colour = NULL) + THEME
  })

  output$killrates <- renderPlot({
    o <- sim() %>% filter(time >= 195)
    tibble(time = o$time,
           `ADCC (marrow)` = 0.90 * 1.25 * o$BPCMo / (o$BPCMo + 2e4),
           `ADCC (solid)` = 0.90 * o$BPCTo / (o$BPCTo + 2e4),
           `cytotoxic` = 1.95 * o$CCTo / (o$CCTo + 1.10)) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.6) +
      labs(title = "사멸률 비교 (1/일)", x = "일", y = "kill rate (1/d)",
           colour = NULL) + THEME
  })

  ## ---------------------------------------------------------- 6. endpoints
  output$t_relapse <- renderText({
    t <- nb_relapse(sim())
    if (is.infinite(t)) "재발 없음 (모델상 완전 소멸)" else sprintf("%.0f 일", t)
  })
  output$t_nadir <- renderText(sprintf("%.3g × 10⁹ cells", min(sim()$BURDEN)))
  output$t_resp <- renderText({
    o <- sim(); b0 <- o$BURDEN[1]; bn <- min(o$BURDEN)
    r <- 1 - bn / b0
    if (bn < 1e-6) "CR (모델 검출한계 이하)" else
      if (r > 0.9) "VGPR / PR" else if (r > 0.3) "MR" else "SD / PD"
  })

  output$pain <- renderPlot({
    sim() %>% filter(time >= 195) %>%
      ggplot(aes(time, PAIN)) +
      geom_line(colour = "#CA6F1E", linewidth = 0.7) +
      labs(title = "통증(이질통) 강도 — 보체 매개",
           subtitle = "코스 내 감쇠는 보체 풀 소모에서 나온다",
           x = "일", y = "상대 강도") + THEME
  })

  output$heme <- renderPlot({
    sim() %>% transmute(time, `ANC (10⁹/L)` = ANC,
                        `PLT /100 (10⁹/L)` = PLT / 100) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.6) +
      geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
      labs(title = "골수억제 (Friberg 구조)", x = "일", y = NULL,
           colour = NULL) + THEME
  })

  output$oto <- renderPlot({
    sim() %>% ggplot(aes(time, OTO)) +
      geom_line(colour = "#943126", linewidth = 0.6) +
      labs(title = "누적 이독성 (백금 누적 노출)", x = "일",
           y = "상대 단위") + THEME
  })

  output$tox_table <- renderTable({
    o <- sim()
    tibble(endpoint = c("ANC 최저", "PLT 최저", "통증 최고", "보체 최저 (CPL)",
                        "누적 이독성", "누적 백금 AUC"),
           value = c(min(o$ANC), min(o$PLT), max(o$PAIN), min(o$CPL),
                     max(o$OTO), max(o$PLATAUC)))
  }, digits = 3)

  ## ---------------------------------------------------------- 7. scenarios
  scen_res <- eventReactive(input$runscen, {
    out <- list()
    withProgress(message = "시나리오 실행 중…", {
      if ("dose" %in% input$scen)     out$dose <- nb_scen_dose()
      if ("perm" %in% input$scen)     out$perm <- nb_scen_perm()
      if ("k322a" %in% input$scen)    out$k322a <- nb_scen_k322a()
      if ("burden" %in% input$scen)   out$burden <- nb_scen_burden()
      if ("cyto" %in% input$scen)     out$cyto <- nb_scen_cytokine()
      if ("retinoid" %in% input$scen) out$retinoid <- nb_scen_retinoid()
      if ("sts" %in% input$scen)      out$sts <- nb_scen_sts()
    })
    out
  })

  output$scen_out <- renderUI({
    r <- scen_res()
    if (!length(r)) return(helpText("축을 고르고 '시나리오 실행'을 누르세요."))
    do.call(tagList, lapply(names(r), function(nm) {
      tagList(h4(nm), renderTable(r[[nm]], digits = 4)())
    }))
  })

  ## --------------------------------------------------------------- 8. MIBG
  mibg <- eventReactive(input$runmibg, {
    MBq <- input$mibg_mci * 37 * input$WT
    withProgress(message = "MIBG 적분 중…", {
      nb_run(e = ev_mibg(0, MBq, as.numeric(input$mibg_sa)),
             tp0 = input$tp0, tm0 = input$tm0, end = 60, delta = 0.05,
             param = c(pset(), list(NETX = as.numeric(input$netx),
                                    KEXCR = input$kexcr,
                                    KIBLOCK = if (input$ki) 0.10 else 1.00)))
    })
  })

  output$mibg_dose <- renderPlot({
    mibg() %>% transmute(time, `전신 (whole body)` = DWB,
                         `종양 (tumour)` = DTU) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.7) +
      geom_hline(yintercept = 2, linetype = 2, colour = "grey40") +
      labs(title = "누적 흡수선량", subtitle = "파선 = 2 Gy 전신선량 (조혈모세포 구제 임계)",
           x = "일", y = "Gy", colour = NULL) + THEME
  })

  output$mibg_act <- renderPlot({
    mibg() %>% transmute(time, `혈액 (blood)` = MBC, `종양 (tumour)` = MBT,
                         `갑상선 (thyroid)` = THY) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, pmax(value, 1e-2), colour = name)) +
      geom_line(linewidth = 0.7) + scale_y_log10() +
      labs(title = "잔류 활성도", x = "일", y = "MBq (log)", colour = NULL) +
      THEME
  })

  output$mibg_tab <- renderTable({
    o <- mibg(); MBq <- input$mibg_mci * 37 * input$WT
    tibble(quantity = c("투여 활성도 (MBq)", "전신선량 (Gy)",
                        "mGy / MBq", "종양선량 (Gy)", "ANC 최저",
                        "갑상선 최고 활성도 (MBq)"),
           value = c(MBq, max(o$DWB), 1000 * max(o$DWB) / MBq, max(o$DTU),
                     min(o$ANC), max(o$THY)))
  }, digits = 3)

  ## -------------------------------------------------------------- 9. cohort
  coh <- eventReactive(input$runcoh, {
    co <- nb_cohort(input$ncoh, input$seed)
    res <- list()
    withProgress(message = "코호트 실행 중 (수 분 소요)…", {
      if ("ctl" %in% input$coh_arms)
        res[["Isotretinoin 단독"]] <- nb_efs(co, immuno = FALSE)
      if ("cog" %in% input$coh_arms)
        res[["COG 전체 (Ab+GM+IL2)"]] <- nb_efs(co, immuno = TRUE)
      if ("noil2" %in% input$coh_arms)
        res[["IL-2 제외 (SIOPEN)"]] <- nb_efs(co, immuno = TRUE, use_il2 = FALSE)
    })
    res
  })

  output$efs_curve <- renderPlot({
    r <- coh()
    grid <- seq(0, 1100, 10)
    purrr::imap_dfr(r, function(t, nm)
      tibble(arm = nm, time = grid,
             efs = sapply(grid, function(g) mean(t > g)))) %>%
      ggplot(aes(time, efs, colour = arm)) + geom_step(linewidth = 0.9) +
      geom_vline(xintercept = 730, linetype = 2, colour = "grey40") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = "무사건 생존 (event-free survival)",
           subtitle = "파선 = 2년 (ANBL0032 종점)", x = "일",
           y = "EFS", colour = NULL) + THEME
  })

  output$efs_tab <- renderTable({
    r <- coh()
    purrr::imap_dfr(r, function(t, nm) {
      p2 <- mean(t > 730); n <- length(t)
      se <- sqrt(p2 * (1 - p2) / n)
      tibble(arm = nm, n = n,
             `1년 EFS` = mean(t > 365), `2년 EFS` = p2,
             `2년 95% CI` = sprintf("%.2f – %.2f", max(0, p2 - 1.96 * se),
                                    min(1, p2 + 1.96 * se)),
             `3년 EFS` = mean(t > 1095))
    })
  }, digits = 3)
}

shinyApp(ui, server)
