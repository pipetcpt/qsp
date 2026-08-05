## =============================================================================
##  mpm_shiny_app.R
##  Malignant pleural mesothelioma QSP model -- interactive dashboard
##
##  Run with:
##      shiny::runApp("mpm_shiny_app.R")
##  (the app sources mpm_mrgsolve_model.R from the same directory)
##
##  The app is organised around the model's one structural commitment: the
##  tumour is a RIND whose only degree of freedom is thickness, and every panel
##  shows a DIFFERENT FUNCTIONAL of that same sheet.  Tab 2 is the conversion
##  table between them; tab 3 is the reason they diverge under treatment.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

options(mpm.no_demo = TRUE)   # do not run the model file's own demo block
source("mpm_mrgsolve_model.R", local = TRUE)

theme_set(theme_bw(base_size = 12))

SCEN_CHOICES <- setNames(names(SCENARIOS),
                         vapply(SCENARIOS, function(s) s$label, character(1)))

fp <- function(h, lam) ifelse(h <= 0, 1, (lam / h) * (1 - exp(-h / lam)))


## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("악성 흉막 중피종 QSP 모델 · Malignant Pleural Mesothelioma"),
  tags$p(style = "color:#555;margin-top:-8px;",
         HTML("종양을 <b>덩어리가 아니라 두께 하나만 자유도로 갖는 껍질(rind)</b>로 기술하고, ",
              "임상 지표들을 모두 그 껍질의 서로 다른 범함수로 계산합니다. ",
              "<i>Tumour written as a rind, not a mass: every endpoint is a different ",
              "functional of one thickness.</i>")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 / Patient"),
      sliderInput("emt", "조직형 축 EMT (0 상피양 · 0.5 이상성 · 1 육종양)",
                  0, 1, 0.25, step = 0.05),
      sliderInput("hpar", "진단 시 벽측 흉막 두께 h_par (mm)", 1, 20, 6, step = 0.5),
      sliderInput("hvis", "진단 시 장측 흉막 두께 h_vis (mm)", 1, 20, 4, step = 0.5),
      sliderInput("hfis", "엽간열 / 종격 두께 h_fis (mm)", 1, 20, 5, step = 0.5),
      sliderInput("plv0", "진단 시 흉수량 (mL)", 0, 3000, 628, step = 50),
      sliderInput("fvcpred", "예측 FVC (L)", 2.0, 5.0, 3.3, step = 0.1),
      hr(),
      h4("치료 / Regimen"),
      selectInput("scen", "시나리오", choices = SCEN_CHOICES, selected = "pemcis"),
      selectInput("scen2", "비교 시나리오", choices = SCEN_CHOICES, selected = "nivoipi"),
      sliderInput("tend", "시뮬레이션 기간 (일)", 180, 1500, 1100, step = 30),
      checkboxInput("folate", "엽산 + 비타민 B12 보충", TRUE),
      hr(),
      h4("탐색 / Explore"),
      sliderInput("phi_x", "탐색용 콜라겐 분율 phi", 0, 0.9, 0.25, step = 0.05),
      helpText(HTML("<small>phi는 세 곳에만 들어갑니다: 두께(측정 편향) · ",
                    "침투길이(전달) · 장측 탄성(폐 포획).</small>"))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## ------------------------------------------------------------------
        tabPanel("1. 환자 프로파일",
          br(),
          fluidRow(
            column(6, h4("기저 상태 / Baseline state"), tableOutput("tbl_baseline")),
            column(6, h4("조직형이 바꾸는 다섯 가지"), tableOutput("tbl_emt"))),
          hr(),
          plotOutput("plt_baseline", height = "330px"),
          helpText(HTML("조직형 축 x는 <b>증식 · 화학요법 저항 · 콜라겐 침착 · PD-L1 · ",
                        "T세포 침윤</b> 다섯 가지를 동시에 움직입니다. ",
                        "앞의 세 개는 화학요법에 불리하고 뒤의 두 개는 면역관문 억제에 ",
                        "유리하므로, 두 치료의 비교 부호가 x의 함수로 뒤집힙니다 (탭 9)."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("2. 껍질 기하학",
          br(),
          h4("두께 - 부피 - mRECIST 변환"),
          tableOutput("tbl_geom"),
          hr(),
          fluidRow(
            column(6, plotOutput("plt_geom", height = "320px")),
            column(6,
              h4("같은 '부분 관해'가 뜻하는 서로 다른 세포 살상"),
              verbatimTextOutput("txt_geom"),
              helpText(HTML("껍질은 dV/dh = A 로 <b>일정</b>합니다. 따라서 두께 30% 감소는 ",
                            "부피 30% 감소입니다. 구형 종양에서는 직경 30% 감소가 부피 ",
                            "66% 감소이므로, 같은 단어가 2.2배 다른 세포 살상을 뜻합니다."))))
        ),

        ## ------------------------------------------------------------------
        tabPanel("3. 약물 침투",
          br(),
          h4("깊이 평균 노출  f = (lambda/h)(1 - exp(-h/lambda))"),
          plotOutput("plt_pen", height = "360px"),
          hr(),
          tableOutput("tbl_pen"),
          helpText(HTML("전신 투여는 껍질의 <b>혈관이 있는 바닥</b>에서, 흉강내 투여는 ",
                        "<b>자유 표면</b>에서 들어옵니다. 방향은 반대지만 식은 같고, ",
                        "두 침투 길이 모두 수 mm입니다. 그래서 흉강내 치료는 ",
                        "잔존 두께가 침투 길이보다 얇아진 뒤에만 해석 가능합니다."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("4. 약동학 (PK)",
          br(),
          plotOutput("plt_pk_plasma", height = "300px"),
          plotOutput("plt_pk_tumour", height = "300px"),
          helpText(HTML("항체는 혈장에서 PD-1을 포화시키고도 (RO > 0.9) ",
                        "종양 간질 농도는 혈장의 수 %에 그칩니다 - 결합부위 장벽. ",
                        "penetration이 제한하는 것은 <b>수용체 점유가 아니라 T세포 침윤</b>입니다."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("5. 종양 반응",
          br(),
          plotOutput("plt_resp", height = "400px"),
          hr(),
          fluidRow(
            column(6, h4("측정이 숨기는 살상"), tableOutput("tbl_hidden")),
            column(6, h4("반응 중 콜라겐 분율"), plotOutput("plt_phi", height = "260px"))),
          helpText(HTML("살아있는 세포는 크게 줄지만 두께는 그만큼 줄지 않습니다. ",
                        "차이는 <b>남아있는 콜라겐과 아직 청소되지 않은 괴사물</b>입니다."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("6. 흉막강 · 호흡역학",
          br(),
          plotOutput("plt_pleura", height = "400px"),
          hr(),
          fluidRow(
            column(6, h4("호흡곤란의 두 원인"), tableOutput("tbl_dyspnoea")),
            column(6, h4("흉막유착술 성공 조건"), verbatimTextOutput("txt_pleurodesis"))),
          helpText(HTML("흉수는 <b>제거 가능한</b> 원인이고 폐 포획은 <b>제거 불가능한</b> ",
                        "원인입니다. 장측 껍질이 뻣뻣하면 폐가 펴지지 않으므로 유착술이 ",
                        "실패합니다 - 모델에서는 유착 형성 속도에 apposition 항으로 들어갑니다."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("7. 면역 미세환경",
          br(),
          plotOutput("plt_immune", height = "400px"),
          helpText(HTML("PD-L1은 EMT 축을 따라 0.20 -> 0.60으로 올라가고 IFN-gamma로 ",
                        "추가 유도됩니다. 항CTLA-4는 <b>림프절</b>에서 (침투 장벽 없음), ",
                        "항PD-1은 <b>종양 안에서</b> (침투 장벽 있음) 작용합니다."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("8. 시나리오 비교",
          br(),
          h4("20개 시나리오 종점 표"),
          tableOutput("tbl_scen"),
          helpText(HTML("medOS는 누적 위험이 ln2에 도달하는 시점입니다. ",
                        "관해율은 개인 궤적이 아니라 집단량이므로 여기서는 최선 mRECIST ",
                        "변화만 보고합니다 (집단 ORR은 mpm_calibration.py 참조)."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("9. 조직형 축 교차",
          br(),
          plotOutput("plt_cross", height = "380px"),
          hr(),
          tableOutput("tbl_cross"),
          helpText(HTML("이 곡선은 <b>예측</b>입니다. 보정에 쓴 것은 CheckMate 743의 ",
                        "전체 집단 한 점뿐이고, 상피양 / 비상피양 분리는 EMT 축이 ",
                        "만들어낸 결과입니다 (관측: 상피양 HR 0.86, 비상피양 HR 0.46)."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("10. 독성",
          br(),
          plotOutput("plt_tox", height = "420px"),
          hr(),
          tableOutput("tbl_tox"),
          helpText(HTML("시스플라틴 신독성이 eGFR을 떨어뜨리면 페메트렉시드 청소율이 ",
                        "함께 떨어져 노출이 올라갑니다 - 병용요법 두 약이 <b>신장을 통해 ",
                        "서로를 증폭</b>하는 양성 되먹임입니다."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("11. 바이오마커",
          br(),
          plotOutput("plt_bio", height = "360px"),
          hr(),
          tableOutput("tbl_bio"),
          helpText(HTML("SMRP는 살아있는 세포 표면에서 <b>떨어져 나오는</b> 산물이므로 ",
                        "약이 실제로 바꾸는 양을 따라가고, 두께보다 먼저 움직입니다. ",
                        "신기능이 떨어지면 청소율이 떨어져 값이 올라가는 교란이 있습니다."))
        ),

        ## ------------------------------------------------------------------
        tabPanel("12. 생존",
          br(),
          plotOutput("plt_surv", height = "380px"),
          hr(),
          tableOutput("tbl_surv"),
          helpText(HTML("위험함수는 종양 부피 · ECOG · FVC · 조직형에 걸려 있습니다. ",
                        "부피 계수가 작고 FVC/ECOG 계수가 크기 때문에, 큰 영상 반응이 ",
                        "작은 생존 이득으로만 이어지는 이 질환의 특징이 재현됩니다."))
        )
      )
    )
  )
)


## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  base_param <- reactive({
    list(EMT = input$emt,
         HPAR0 = input$hpar / 10, HVIS0 = input$hvis / 10, HFIS0 = input$hfis / 10,
         PLV0 = input$plv0, FVC_PRED = input$fvcpred,
         FOLATE = as.numeric(input$folate))
  })

  simA <- reactive({
    run_scenario(input$scen, emt = input$emt, end = input$tend,
                 extra = base_param())
  })
  simB <- reactive({
    run_scenario(input$scen2, emt = input$emt, end = input$tend,
                 extra = base_param())
  })
  simBoth <- reactive({
    a <- simA(); b <- simB()
    a$arm <- SCENARIOS[[input$scen]]$label
    b$arm <- SCENARIOS[[input$scen2]]$label
    rbind(a, b)
  })

  ## ---------------- 1. patient profile -------------------------------------
  output$tbl_baseline <- renderTable({
    o <- simA()[1, ]
    data.frame(
      항목 = c("벽측 두께 (mm)", "장측 두께 (mm)", "엽간열 두께 (mm)",
               "mRECIST 합 (mm)", "종양 부피 (cm3)", "콜라겐 분율 phi",
               "부가 탄성 (cmH2O/L)", "흉수 (mL)", "FVC (L)", "FVC (% 예측)",
               "ECOG (연속)"),
      값 = c(round(o$h_par_mm, 2), round(o$h_vis_mm, 2), round(o$h_fis_mm, 2),
             round(o$mRECIST, 1), round(o$Vtumour, 0), round(o$phi, 3),
             round(o$Eadd, 1), round(o$plv_out, 0), round(o$FVC, 2),
             round(o$FVCpct, 0), round(o$ECOG, 2)))
  }, striped = TRUE)

  output$tbl_emt <- renderTable({
    x <- input$emt
    data.frame(
      기전 = c("증식률 kg", "화학요법 살상", "콜라겐 침착 kcol",
               "PD-L1 발현", "T세포 침윤"),
      배수 = c(sprintf("x %.2f", 1 + 0.90 * x),
               sprintf("x %.2f", 1 - 0.55 * x),
               sprintf("x %.2f", 1 + 3.20 * x),
               sprintf("%.2f", 0.20 * (1 + 2.0 * x)),
               sprintf("x %.2f", 1 + 1.00 * x)),
      방향 = c("화학요법에 불리", "화학요법에 불리", "화학요법에 불리",
               "면역관문에 유리", "면역관문에 유리"))
  }, striped = TRUE)

  output$plt_baseline <- renderPlot({
    xs <- seq(0, 1, 0.02)
    d <- rbind(
      data.frame(x = xs, v = 1 + 0.90 * xs, what = "증식률 (chemo 불리)"),
      data.frame(x = xs, v = 1 - 0.55 * xs, what = "화학요법 살상 (불리)"),
      data.frame(x = xs, v = (1 + 3.20 * xs) / 4.2, what = "콜라겐 (정규화, 불리)"),
      data.frame(x = xs, v = 0.20 * (1 + 2.0 * xs) / 0.2 / 3, what = "PD-L1 (정규화, 유리)"),
      data.frame(x = xs, v = (1 + 1.0 * xs) / 2, what = "T세포 침윤 (정규화, 유리)"))
    ggplot(d, aes(x, v, colour = what)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = input$emt, linetype = 2) +
      labs(x = "EMT 축 x", y = "상대값", colour = NULL,
           title = "하나의 축이 다섯 개를 동시에 움직인다")
  })

  ## ---------------- 2. geometry --------------------------------------------
  output$tbl_geom <- renderTable({
    A_par <- 1000 * 0.45; A_vis <- 900 * 0.35; A_fis <- 400 * 0.50
    hs <- c(1, 2, 3, 5, 8, 12, 20)
    data.frame(
      `두께 h (mm)` = hs,
      `벽측 부피 (mL)` = round(A_par * hs / 10, 0),
      `세 잎 합계 (mL)` = round((A_par + A_vis + A_fis) * hs / 10, 0),
      `mRECIST 기여 (mm)` = round(4 * hs, 0),
      check.names = FALSE)
  }, striped = TRUE)

  output$plt_geom <- renderPlot({
    f <- seq(0.3, 1, 0.01)
    d <- rbind(data.frame(f = f, kill = 1 - f, shape = "RIND  (V = A h)"),
               data.frame(f = f, kill = 1 - f^3, shape = "SPHERE (V = 4/3 pi r^3)"))
    ggplot(d, aes(100 * (1 - f), 100 * kill, colour = shape)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = 30, linetype = 2) +
      labs(x = "측정된 선형 감소 (%)", y = "함축된 부피 살상 (%)", colour = NULL,
           title = "같은 30%가 두 배 다른 살상을 뜻한다")
  })

  output$txt_geom <- renderPrint({
    A_par <- 1000 * 0.45
    cat(sprintf("dV/dh (벽측 잎)          : %.0f mL / mm\n", A_par / 10))
    cat(sprintf("껍질  30%% 두께 감소      : %.1f%% 부피 살상\n", 30))
    cat(sprintf("구형  30%% 직경 감소      : %.1f%% 부피 살상\n", 100 * (1 - 0.7^3)))
    cat(sprintf("비율                     : %.2f 배\n", (1 - 0.7^3) / 0.30))
  })

  ## ---------------- 3. penetration -----------------------------------------
  output$plt_pen <- renderPlot({
    phi <- input$phi_x
    hs <- seq(0.05, 2.5, 0.01)
    mk <- function(l0, lab) data.frame(h = hs * 10,
                                       f = fp(hs, l0 * (1 - phi)^1.5), route = lab)
    d <- rbind(mk(0.25, "소분자 (전신)"), mk(0.060, "IgG 항체 (전신)"),
               mk(0.050, "T세포 침윤"), mk(0.35, "흉강내 투여 (자유 표면)"))
    ggplot(d, aes(h, f, colour = route)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = input$hpar, linetype = 2) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "껍질 두께 h (mm)", y = "깊이 평균 노출 분율 f", colour = NULL,
           title = sprintf("phi = %.2f 에서의 침투", phi),
           subtitle = "점선 = 현재 환자의 벽측 두께")
  })

  output$tbl_pen <- renderTable({
    penetration_table(phis = c(input$phi_x))
  }, striped = TRUE)

  ## ---------------- 4. PK ---------------------------------------------------
  output$plt_pk_plasma <- renderPlot({
    simBoth() %>%
      select(time, arm, Ccis, Cpem, Cniv, Cipi, Cbev) %>%
      pivot_longer(-c(time, arm)) %>% filter(value > 1e-6) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line() +
      facet_wrap(~name, scales = "free_y") +
      scale_y_log10() + labs(x = "일", y = "혈장 농도 (mg/L)", colour = NULL,
                             title = "혈장 PK")
  })

  output$plt_pk_tumour <- renderPlot({
    simBoth() %>%
      select(time, arm, CIST, PEMT, NIVT, f_sm, f_ab, f_Tcell) %>%
      pivot_longer(-c(time, arm)) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line() +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "일", y = NULL, colour = NULL,
           title = "종양 내 농도와 침투 분율 (f 는 두께가 줄면 올라간다)")
  })

  ## ---------------- 5. tumour response --------------------------------------
  output$plt_resp <- renderPlot({
    simBoth() %>%
      select(time, arm, mRECIST, Vtumour, Tviable, h_par_mm, h_vis_mm, phi) %>%
      pivot_longer(-c(time, arm)) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "일", y = NULL, colour = NULL,
           title = "같은 껍질의 서로 다른 범함수")
  })

  output$tbl_hidden <- renderTable({
    do.call(rbind, lapply(c(input$scen, input$scen2), function(k) {
      o <- run_scenario(k, emt = input$emt, end = min(input$tend, 700),
                        extra = base_param())
      b <- o[1, ]
      i <- which.min(o$mRECIST)
      data.frame(시나리오 = k,
                 `최선 mRECIST (%)` = round(100 * (min(o$mRECIST) - b$mRECIST) / b$mRECIST, 1),
                 `그 시점 생존세포 (%)` = round(100 * (o$Tviable[i] - b$Tviable) / b$Tviable, 1),
                 `숨겨진 살상 (%p)` = round(
                   100 * (o$Tviable[i] - b$Tviable) / b$Tviable -
                   100 * (min(o$mRECIST) - b$mRECIST) / b$mRECIST, 1),
                 check.names = FALSE)
    }))
  }, striped = TRUE)

  output$plt_phi <- renderPlot({
    ggplot(simBoth(), aes(time, phi, colour = arm)) + geom_line(linewidth = 1) +
      labs(x = "일", y = "콜라겐 분율 phi", colour = NULL,
           title = "반응할수록 껍질은 더 섬유화된다")
  })

  ## ---------------- 6. pleural space and mechanics --------------------------
  output$plt_pleura <- renderPlot({
    simBoth() %>%
      select(time, arm, plv_out, Eadd, FVC, vexp_out, stomablk, SYMPH) %>%
      pivot_longer(-c(time, arm)) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "일", y = NULL, colour = NULL,
           title = "흉수 · 부가 탄성 · 폐 팽창 가능 분율")
  })

  output$tbl_dyspnoea <- renderTable({
    o <- simA()
    idx <- c(1, which(o$time == 90)[1], which(o$time == 180)[1],
             which(o$time == 360)[1])
    idx <- idx[!is.na(idx)]
    data.frame(
      `일` = o$time[idx],
      `흉수 (mL)` = round(o$plv_out[idx], 0),
      `부가 탄성` = round(o$Eadd[idx], 1),
      `FVC (L)` = round(o$FVC[idx], 2),
      `FVC 흉수 없을 때 (L)` = round(o$FVC[idx] /
                                      (1 - 0.55 * o$plv_out[idx] / 3500), 2),
      check.names = FALSE)
  }, striped = TRUE)

  output$txt_pleurodesis <- renderPrint({
    o <- simA()[1, ]
    cat(sprintf("장측 두께      : %.2f mm\n", o$h_vis_mm))
    cat(sprintf("장측 phi       : %.3f\n", o$phi_vis))
    cat(sprintf("부가 탄성      : %.1f cmH2O/L\n", o$Eadd))
    cat(sprintf("포획 역치      : 14.5 cmH2O/L (흉막 압력계 문헌)\n"))
    cat(sprintf("폐 팽창 가능   : %.2f (유착술 형성 속도에 곱해짐)\n", o$vexp_out))
    cat(if (o$Eadd >= 14.5) "-> 폐가 포획되어 있으므로 유착술은 실패할 가능성이 높다\n"
        else "-> 폐가 펴질 수 있으므로 유착술이 성립할 수 있다\n")
  })

  ## ---------------- 7. immune -----------------------------------------------
  output$plt_immune <- renderPlot({
    simBoth() %>%
      select(time, arm, TEFF, TREG, PDL1, IFNG, RO_PD1, RO_CTLA4, f_Tcell, TAM) %>%
      pivot_longer(-c(time, arm)) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "일", y = NULL, colour = NULL,
           title = "면역 구획: 점유는 포화되어도 침윤은 두께에 갇힌다")
  })

  ## ---------------- 8. scenario comparison ----------------------------------
  output$tbl_scen <- renderTable({
    endpoint_table(emt = input$emt)
  }, striped = TRUE, digits = 2)

  ## ---------------- 9. crossover --------------------------------------------
  cross <- reactive({ crossover_curve(seq(0, 1, 0.1)) })

  output$plt_cross <- renderPlot({
    cc <- cross()
    ggplot(cc, aes(EMT)) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_line(aes(y = delta), linewidth = 1.1, colour = "#b3564a") +
      geom_point(aes(y = delta), size = 2) +
      labs(x = "EMT 축 x (0 상피양 - 1 육종양)",
           y = "중앙 생존 차이  IO - 화학요법 (개월)",
           title = "치료 비교의 부호가 조직형 축을 따라 뒤집힌다",
           subtitle = "보정에 쓴 것은 이 곡선 위의 한 점뿐이다")
  })

  output$tbl_cross <- renderTable({
    cross() %>% mutate(across(where(is.numeric), ~round(.x, 2)))
  }, striped = TRUE)

  ## ---------------- 10. toxicity --------------------------------------------
  output$plt_tox <- renderPlot({
    simBoth() %>%
      select(time, arm, ANC, GFR, IRAE, NEURO, PTK, CACHEX) %>%
      pivot_longer(-c(time, arm)) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "일", y = NULL, colour = NULL, title = "독성")
  })

  output$tbl_tox <- renderTable({
    do.call(rbind, lapply(c(input$scen, input$scen2), function(k) {
      o <- run_scenario(k, emt = input$emt, end = input$tend, extra = base_param())
      data.frame(시나리오 = k,
                 `ANC 최저 (10^9/L)` = round(min(o$ANC), 2),
                 `ANC 최저일` = o$time[which.min(o$ANC)],
                 `eGFR 최저` = round(min(o$GFR), 0),
                 `irAE 최고` = round(max(o$IRAE), 2),
                 `신경병증 지수` = round(max(o$NEURO), 2),
                 check.names = FALSE)
    }))
  }, striped = TRUE)

  ## ---------------- 11. biomarker -------------------------------------------
  output$plt_bio <- renderPlot({
    d <- simBoth() %>% group_by(arm) %>%
      mutate(SMRP_rel = 100 * (SMRP - SMRP[10]) / SMRP[10],
             mR_rel   = 100 * (mRECIST - mRECIST[1]) / mRECIST[1],
             viab_rel = 100 * (Tviable - Tviable[1]) / Tviable[1]) %>% ungroup()
    d %>% select(time, arm, SMRP_rel, mR_rel, viab_rel) %>%
      pivot_longer(-c(time, arm)) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~arm) +
      labs(x = "일", y = "기저 대비 변화 (%)", colour = NULL,
           title = "SMRP는 두께가 아니라 살아있는 세포를 따라간다")
  })

  output$tbl_bio <- renderTable({
    o <- simA()
    idx <- vapply(c(0, 30, 60, 90, 180, 270, 360),
                  function(d) which(o$time == d)[1], numeric(1))
    idx <- idx[!is.na(idx)]
    b <- o[max(idx[1], 10), ]
    data.frame(`일` = o$time[idx],
               `SMRP (nM)` = round(o$SMRP[idx], 2),
               `생존세포 (%)` = round(100 * (o$Tviable[idx] - o$Tviable[1]) / o$Tviable[1], 1),
               `mRECIST (%)` = round(100 * (o$mRECIST[idx] - o$mRECIST[1]) / o$mRECIST[1], 1),
               `eGFR` = round(o$GFR[idx], 0),
               check.names = FALSE)
  }, striped = TRUE)

  ## ---------------- 12. survival --------------------------------------------
  output$plt_surv <- renderPlot({
    ggplot(simBoth(), aes(time / 30.44, Surv, colour = arm)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0.5, linetype = 2) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "개월", y = "생존 확률 exp(-누적위험)", colour = NULL,
           title = "모형 생존 곡선")
  })

  output$tbl_surv <- renderTable({
    do.call(rbind, lapply(c(input$scen, input$scen2), function(k) {
      o <- run_scenario(k, emt = input$emt, end = max(input$tend, 800),
                        extra = base_param())
      data.frame(시나리오 = k,
                 `중앙 OS (개월)` = round(median_os(o), 1),
                 `PFS (개월)` = round(pfs_months(o), 1),
                 `1년 생존` = round(o$Surv[which(o$time == 365)[1]], 3),
                 `2년 생존` = round(o$Surv[which(o$time == 730)[1]], 3),
                 check.names = FALSE)
    }))
  }, striped = TRUE)
}

shinyApp(ui, server)
