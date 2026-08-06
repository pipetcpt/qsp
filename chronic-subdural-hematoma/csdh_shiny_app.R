##############################################################################
## csdh_shiny_app.R
## 만성 경막하 혈종 (cSDH) QSP 인터랙티브 대시보드
## Chronic Subdural Haematoma — QSP interactive dashboard
##
## The app is organised around the model's one structural claim: the haematoma
## volume is the FIXED POINT of
##        dV/dt = J_exudation + J_rebleed - J_absorption
## so the first thing every tab shows is which of those three terms the user's
## chosen therapy is actually moving.  Tab 3 (플럭스 균형) is the tab to read
## first; everything else is downstream of it.
##
## Run:
##   shiny::runApp("csdh_shiny_app.R")
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##############################################################################

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("csdh_mrgsolve_model.R")   # defines mod, csdh_sim, ev_* , CSDH_* helpers

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#EEF2F7", colour = NA),
        legend.position = "bottom")

PAL <- c("#2C5F8D", "#8B2635", "#276749", "#9C5B1E", "#553C7B",
         "#B8860B", "#C9548A", "#333333")

##############################################################################
## UI
##############################################################################

ui <- fluidPage(
  titlePanel("만성 경막하 혈종 (cSDH) QSP 시뮬레이터 — Chronic Subdural Haematoma"),
  tags$p(style = "color:#8B0000;font-weight:600;margin-top:-8px;",
         "dV/dt = J_exudation + J_rebleed − J_absorption  ·  ",
         "혈종은 흡수에 실패한 응고물이 아니라 '분비하는 기관'이다"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("① 환자 (Patient)"),
      sliderInput("age", "연령 Age (yr)", 45, 95, 76, 1),
      sliderInput("vres", "두개내 예비공간 V_RES (mL) — 위축",
                  0, 60, 30, 1),
      helpText(tags$small(
        "V_RES는 모델에 두 번, 반대 부호로 들어간다: 압력을 완충해 더 큰 ",
        "혈종을 견디게 하고(이롭다), 동시에 뇌 재팽창을 늦춰 배액 후 공간을 ",
        "남긴다(해롭다).")),
      sliderInput("vpres", "발현 시 부피 V_pres (mL)", 40, 140, 78, 2),

      h4("② 시술 (Procedure)"),
      selectInput("surg", "수술",
                  c("없음 (보존적)" = "none",
                    "천공배액 (배액관 없음)" = "bh",
                    "천공배액 + 배액관 48h" = "bhd",
                    "개두술 + 막 제거" = "crani")),
      sliderInput("wash", "세척 정도 wash (낮을수록 깨끗)", 0.05, 1.0, 0.15, 0.05),
      helpText(tags$small("wash는 공동 내 '농도'를 낮춘다 — 헴 스위치를 끄는 조작.")),
      checkboxInput("mmae", "MMA 색전술 (MMAE)", FALSE),
      sliderInput("mmae_t", "MMAE 시점 (일)", 0, 60, 0, 1),

      h4("③ 약물 (Drugs)"),
      checkboxInput("dex", "덱사메타손 16 mg/d × 2주 감량", FALSE),
      checkboxInput("atv", "아토르바스타틴 20 mg × 8주", FALSE),
      checkboxInput("txa", "트라넥삼산 750 mg/d", FALSE),
      checkboxInput("doac", "아픽사반 재개", FALSE),
      sliderInput("doac_t", "DOAC 재개 시점 (일)", 0, 90, 30, 1),

      h4("④ 시뮬레이션"),
      sliderInput("tend", "관찰 기간 (일)", 60, 365, 180, 5),
      actionButton("go", "실행 (Run)", class = "btn-primary btn-block"),
      hr(),
      helpText(tags$small(
        "모든 곡선은 50-ODE 모델의 출력이며, 파라미터는 Santarius 2009 · ",
        "ATOCH 2018 · Dex-CSDH 2020 · EMBOLISE 2024에 맞추어 보정되었다. ",
        "교육·연구용이며 임상 의사결정에 사용할 수 없다."))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ------------------------------------------------------------
        tabPanel("1. 환자 프로파일",
          br(),
          fluidRow(
            column(6, h4("발현까지의 자연사 (run-in)"),
                   plotOutput("p_runin", height = "300px")),
            column(6, h4("요약"), tableOutput("t_summary"))),
          hr(),
          h4("이 환자에서 V_RES의 두 부호"),
          plotOutput("p_tworole", height = "300px"),
          helpText("왼쪽 축(발현 시 정중선 전위)과 오른쪽 축(재수술 확률)이 ",
                   "위축이 커질 때 반대 방향으로 움직인다.")
        ),

        ## ------------------------------------------------------------
        tabPanel("2. 약물 PK",
          br(), h4("혈장 및 혈종 내 농도"),
          plotOutput("p_pk", height = "420px"),
          hr(),
          h4("TXA: 혈장 t½ 3시간 vs 작용부위 t½ ≈ 6일"),
          plotOutput("p_txa", height = "300px"),
          helpText("TXA는 확산이 아니라 유출액(J_ex)에 실려 공동으로 들어가고 ",
                   "흡수(J_abs)로 씻겨 나간다. 따라서 부위 반감기는 V/J_abs로 ",
                   "결정되며, 혈장 약동학과 무관하게 1–2주가 걸린다.")
        ),

        ## ------------------------------------------------------------
        tabPanel("3. 플럭스 균형 ★",
          br(),
          h4("dV/dt = J_ex + J_rb − J_abs — 세 항의 시간 경과"),
          plotOutput("p_flux", height = "360px"),
          hr(),
          fluidRow(
            column(6, h4("순 변화율 dV/dt"), plotOutput("p_net", height="280px")),
            column(6, h4("Starling 인자"), plotOutput("p_starling", height="280px"))),
          helpText("치료가 어느 항을 건드리는지가 이 모델의 전부다. 수술은 STOCK을, ",
                   "색전술은 SOURCE를, 스타틴은 Lp와 σ를 동시에 움직인다.")
        ),

        ## ------------------------------------------------------------
        tabPanel("4. 혈종 부피·영상",
          br(),
          fluidRow(
            column(6, h4("부피 (mL)"), plotOutput("p_vol", height="290px")),
            column(6, h4("최대 두께 (mm)"), plotOutput("p_thick", height="290px"))),
          hr(),
          fluidRow(
            column(6, h4("CT 밀도 (HU)"), plotOutput("p_hu", height="290px")),
            column(6, h4("막 면적·개방성"), plotOutput("p_area", height="290px")))
        ),

        ## ------------------------------------------------------------
        tabPanel("5. 신생막·신생혈관",
          br(),
          fluidRow(
            column(6, h4("혈관 집단"), plotOutput("p_vessel", height="290px")),
            column(6, h4("성숙도와 누출 인자"), plotOutput("p_mature", height="290px"))),
          hr(),
          h4("사이토카인·프로테아제"),
          plotOutput("p_cyto", height = "300px")
        ),

        ## ------------------------------------------------------------
        tabPanel("6. 헴 스위치 (쌍안정성) ★",
          br(),
          h4("헴 구동 haem_drive = hb²/(1+hb²)"),
          plotOutput("p_switch", height = "300px"),
          hr(),
          h4("세척 정도를 바꾸면 같은 수술이 치유와 재발로 갈라진다"),
          plotOutput("p_bistable", height = "320px"),
          helpText("Hill 지수를 1로 낮추면 이 그림의 모든 곡선이 하나로 겹치고, ",
                   "모델은 어떤 수술로도 낫지 않는 단일 고정점만 갖는다.")
        ),

        ## ------------------------------------------------------------
        tabPanel("7. 섬유소용해 루프",
          br(),
          fluidRow(
            column(6, h4("t-PA · 플라스민"), plotOutput("p_lysis", height="290px")),
            column(6, h4("FDP와 지혈 능력 H"), plotOutput("p_fdp", height="290px"))),
          hr(),
          h4("증폭 인자 1/(1−(1−H)) — 루프가 재출혈을 몇 배로 키우는가"),
          plotOutput("p_gain", height = "280px")
        ),

        ## ------------------------------------------------------------
        tabPanel("8. 뇌 역학",
          br(),
          fluidRow(
            column(6, h4("ICP (양방향)"), plotOutput("p_icp", height="290px")),
            column(6, h4("정중선 전위 (mm)"), plotOutput("p_mls", height="290px"))),
          hr(),
          fluidRow(
            column(6, h4("압박 결손과 유효 예비"),
                   plotOutput("p_comp", height="290px")),
            column(6, h4("잔존 공간 = 재발의 기질"),
                   plotOutput("p_space", height="290px")))
        ),

        ## ------------------------------------------------------------
        tabPanel("9. 임상 엔드포인트",
          br(),
          fluidRow(
            column(6, h4("재수술 확률"), plotOutput("p_reop", height="290px")),
            column(6, h4("양호한 결과 확률 (mRS 0–3)"),
                   plotOutput("p_fav", height="290px"))),
          hr(),
          fluidRow(
            column(6, h4("증상·인지"), plotOutput("p_symp", height="290px")),
            column(6, h4("안전성: 혈당·감염·근병증·혈전"),
                   plotOutput("p_safety", height="290px")))
        ),

        ## ------------------------------------------------------------
        tabPanel("10. 시나리오 비교",
          br(),
          h4("16개 표준 시나리오"),
          checkboxGroupInput("scn_pick", NULL,
            choices = c("보존적"="S01","천공배액"="S02","천공+배액관"="S03",
                        "MMAE 단독"="S05","수술+MMAE"="S06","덱사메타손"="S08",
                        "수술+덱사"="S09","아토르바스타틴"="S10",
                        "TXA"="S12","수술+TXA"="S13","전부"="S16"),
            selected = c("S01","S03","S06","S10"), inline = TRUE),
          plotOutput("p_scn", height = "380px"),
          hr(), h4("엔드포인트 표"), DTOutput("t_scn")
        ),

        ## ------------------------------------------------------------
        tabPanel("11. MMAE 시간 경과 ★",
          br(),
          h4("왜 30일 엔드포인트는 음성이고 180일은 양성인가"),
          plotOutput("p_mmae", height = "340px"),
          hr(),
          h4("절대 위험 감소의 시간 경과"),
          plotOutput("p_mmae_ard", height = "300px"),
          helpText("색전술은 STOCK이 아니라 SOURCE를 없앤다. 기존 신생혈관은 ",
                   "τ ≈ 9–18일로 감쇠하므로 J_ex의 감소는 적분으로 나타나고, ",
                   "따라서 조기 엔드포인트에서는 차이가 작을 수밖에 없다. ",
                   "EMBOLISE(180일, 양성)와 MAGIC-MT(90일, 음성)의 불일치는 ",
                   "이 구조에서 유도되며 가정으로 넣은 것이 아니다.")
        ),

        ## ------------------------------------------------------------
        tabPanel("12. 문헌 대조",
          br(), h4("모델 vs 발표된 시험 수치"),
          DTOutput("t_ledger"),
          hr(),
          helpText("이 표는 모델의 성적표이며, 어긋나는 행은 숨기지 않는다. ",
                   "README.md의 '가장 노출된 예측' 절을 함께 볼 것.")
        )
      )
    )
  )
)

##############################################################################
## SERVER
##############################################################################

server <- function(input, output, session) {

  ## ---- build the requested scenario --------------------------------
  sim <- eventReactive(input$go, {
    doses <- NULL
    if (input$dex)  doses <- c(doses, ev_dex())
    if (input$atv)  doses <- c(doses, ev_atv())
    if (input$txa)  doses <- c(doses, ev_txa(days = input$tend))
    if (input$doac) doses <- c(doses, ev_doac(t0 = input$doac_t,
                                              days = input$tend - input$doac_t))
    surg <- switch(input$surg, none = NA, bh = 0, bhd = 0, crani = 0)
    dd   <- if (identical(input$surg, "bhd")) 2 else 0
    pars <- list(V_RES = input$vres, AGE = input$age)
    if (identical(input$surg, "crani")) {
      pars$K_FUSEX <- 0.20     # membranectomy removes extent
      pars$K_SEPT  <- 0.002
    }
    withProgress(message = "50-ODE 시스템 적분 중...", value = 0.4, {
      csdh_sim("user", tend = input$tend, surgery = surg,
               drain_days = dd, wash = input$wash,
               embolise_at = if (input$mmae) input$mmae_t else NA,
               doses = doses, pars = pars, V_pres = input$vpres)
    })
  }, ignoreNULL = FALSE)

  runin <- eventReactive(input$go, {
    m <- param(mod, list(V_RES = input$vres))
    m %>% init(csdh_fresh()) %>% mrgsim(end = 260, delta = 0.5) %>% as_tibble()
  }, ignoreNULL = FALSE)

  lineplot <- function(d, cols, labs, ylab, hline = NULL) {
    dd <- d %>% select(time, all_of(cols)) %>%
      pivot_longer(-time, names_to = "v", values_to = "y") %>%
      mutate(v = factor(v, levels = cols, labels = labs))
    g <- ggplot(dd, aes(time, y, colour = v)) +
      geom_line(linewidth = 0.85) +
      scale_colour_manual(values = PAL[seq_along(cols)], name = NULL) +
      labs(x = "시간 (일)", y = ylab) + THEME
    if (!is.null(hline))
      g <- g + geom_hline(yintercept = hline, linetype = 2, colour = "#8B0000")
    g
  }

  ## ---- 1. patient profile ------------------------------------------
  output$p_runin <- renderPlot({
    d <- runin()
    lineplot(d, c("VHEM"), c("혈종 부피"), "mL",
             hline = input$vpres) +
      annotate("text", x = 5, y = input$vpres * 1.06,
               label = "발현 문턱", hjust = 0, size = 3.2, colour = "#8B0000")
  })

  output$t_summary <- renderTable({
    d <- sim(); r <- runin()
    i <- which(r$VHEM >= input$vpres)[1]
    data.frame(
      항목 = c("발현까지 잠재기 (일)", "발현 시 두께 (mm)",
               "발현 시 정중선 전위 (mm)", "관찰 종료 시 부피 (mL)",
               "관찰 종료 시 두께 (mm)", "재수술 확률",
               "양호한 결과 확률", "최대 혈당 (mmol/L)"),
      값 = c(sprintf("%.1f", if (is.na(i)) NA else r$time[i]),
             sprintf("%.1f", d$dmax_[1]), sprintf("%.2f", d$MLS_[1]),
             sprintf("%.1f", tail(d$VHEM, 1)), sprintf("%.1f", tail(d$dmax_, 1)),
             sprintf("%.3f", tail(d$P_REOP, 1)),
             sprintf("%.3f", tail(d$P_FAV, 1)), sprintf("%.1f", max(d$GLU))))
  }, striped = TRUE)

  output$p_tworole <- renderPlot({
    vs <- c(5, 12, 20, 30, 40, 50, 60)
    rows <- lapply(vs, function(vr) {
      p <- list(V_RES = vr)
      d <- csdh_sim("s", 180, surgery = 0, drain_days = 2, pars = p)
      data.frame(V_RES = vr, MLS = d$MLS_[1],
                 reop = tail(d$P_REOP, 1))
    })
    dd <- bind_rows(rows)
    ggplot(dd, aes(V_RES)) +
      geom_line(aes(y = MLS, colour = "발현 시 MLS (mm) — 위축이 지켜준다"),
                linewidth = 1) +
      geom_point(aes(y = MLS, colour = "발현 시 MLS (mm) — 위축이 지켜준다")) +
      geom_line(aes(y = reop * 25, colour = "재수술 확률 ×25 — 위축이 해친다"),
                linewidth = 1) +
      geom_point(aes(y = reop * 25, colour = "재수술 확률 ×25 — 위축이 해친다")) +
      scale_y_continuous("발현 시 정중선 전위 (mm)",
                         sec.axis = sec_axis(~./25, name = "재수술 확률")) +
      scale_colour_manual(values = c("#2C5F8D", "#8B2635"), name = NULL) +
      labs(x = "두개내 예비공간 V_RES (mL)") + THEME
  })

  ## ---- 2. PK --------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim()
    p <- as.list(param(mod))
    dd <- d %>% transmute(time,
      `덱사메타손 혈장 (µg/L)` = DEXC / p$DEX_V1 * 1000,
      `아토르바스타틴 혈장 (µg/L)` = ATVC / p$ATV_V1 * 1000,
      `TXA 혈장 (mg/L)` = TXAC / p$TXA_V1,
      `TXA 공동 (mg/L)` = TXAH / pmax(VHEM, 1e-4) * 1000,
      `아픽사반 혈장 (ng/mL)` = DOACC / p$DOAC_V1 * 1000,
      `아픽사반 공동 (ng/mL)` = DOACH / pmax(VHEM, 1e-4) * 1e6) %>%
      pivot_longer(-time, names_to = "v", values_to = "y")
    ggplot(dd, aes(time, y, colour = v)) + geom_line(linewidth = 0.8) +
      facet_wrap(~v, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = rep(PAL, 3), guide = "none") +
      labs(x = "시간 (일)", y = NULL) + THEME
  })

  output$p_txa <- renderPlot({
    d <- csdh_sim("txa", 120, surgery = 0, drain_days = 2,
                  doses = ev_txa(days = 120),
                  pars = list(V_RES = input$vres))
    p <- as.list(param(mod))
    dd <- d %>% transmute(time, 혈장 = TXAC / p$TXA_V1,
                          공동 = TXAH / pmax(VHEM, 1e-4) * 1000) %>%
      pivot_longer(-time, names_to = "v", values_to = "y")
    ggplot(dd, aes(time, y, colour = v)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = c("#2C5F8D", "#8B2635"), name = NULL) +
      labs(x = "시간 (일)", y = "TXA 농도 (mg/L)") + THEME
  })

  ## ---- 3. flux ledger ----------------------------------------------
  output$p_flux <- renderPlot({
    lineplot(sim(), c("JEX", "JRB", "JABS"),
             c("J_ex 유출", "J_rb 재출혈", "J_abs 흡수"), "mL/day")
  })
  output$p_net <- renderPlot({
    d <- sim() %>% mutate(net = JEX + JRB - JABS)
    ggplot(d, aes(time, net)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "#8B0000") +
      geom_line(linewidth = 0.9, colour = "#2C5F8D") +
      labs(x = "시간 (일)", y = "dV/dt (mL/day)") + THEME
  })
  output$p_starling <- renderPlot({
    lineplot(sim(), c("DPNET", "LPEFF", "SIGEFF"),
             c("ΔP 순여과압 (mmHg)", "Lp 유효 전도도", "σ 반사계수"),
             "값")
  })

  ## ---- 4. volume / imaging -----------------------------------------
  output$p_vol   <- renderPlot(lineplot(sim(), "VHEM", "부피", "mL"))
  output$p_thick <- renderPlot(lineplot(sim(), "dmax_", "최대 두께",
                                        "mm", hline = 10))
  output$p_hu    <- renderPlot(lineplot(sim(), "HU_", "CT 밀도", "HU"))
  output$p_area  <- renderPlot(lineplot(sim(), c("A_mem_", "f_pat_"),
                                        c("막 면적 (cm²)", "개방성 f_pat"), "값"))

  ## ---- 5. membrane / vessels ---------------------------------------
  output$p_vessel <- renderPlot(lineplot(sim(), c("NCAP", "NMAT", "NMAC"),
      c("미성숙 혈관 N_CAP", "성숙 혈관 N_MAT", "대식세포 N_MAC"), "밀도 지수"))
  output$p_mature <- renderPlot(lineplot(sim(), c("SIGEFF", "LPEFF", "NPC"),
      c("σ 반사계수", "Lp 전도도", "혈관주위세포 피복"), "값"))
  output$p_cyto <- renderPlot({
    d <- sim() %>% transmute(time, `VEGF (pg/mL)` = CVEGF,
      `Ang-2 (ng/mL)` = CANG2, `MMP-9 (ng/mL)` = CMMP9,
      `IL-6 (pg/mL)` = CIL6, `IL-8 (pg/mL)` = CIL8) %>%
      pivot_longer(-time, names_to = "v", values_to = "y")
    ggplot(d, aes(time, y, colour = v)) + geom_line(linewidth = 0.8) +
      facet_wrap(~v, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = rep(PAL, 2), guide = "none") +
      labs(x = "시간 (일)", y = NULL) + THEME
  })

  ## ---- 6. the switch ------------------------------------------------
  output$p_switch <- renderPlot({
    lineplot(sim(), c("HDRIVE", "CHB_"),
             c("헴 구동 (0–1)", "공동 헤모글로빈 (g/dL)"), "값")
  })
  output$p_bistable <- renderPlot({
    ws <- c(0.05, 0.15, 0.30, 0.55, 1.00)
    dd <- bind_rows(lapply(ws, function(w) {
      d <- csdh_sim("w", 180, surgery = 0, drain_days = 2, wash = w,
                    pars = list(V_RES = input$vres))
      data.frame(time = d$time, V = d$VHEM, wash = sprintf("wash %.2f", w))
    }))
    ggplot(dd, aes(time, V, colour = wash)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 0, linetype = 3) +
      scale_colour_manual(values = PAL[seq_along(ws)], name = NULL) +
      labs(x = "시간 (일)", y = "혈종 부피 (mL)") + THEME
  })

  ## ---- 7. fibrinolysis ----------------------------------------------
  output$p_lysis <- renderPlot(lineplot(sim(), c("CTPA_", "CPLS", "CPAI"),
      c("t-PA (ng/mL)", "플라스민 활성", "PAI-1 (ng/mL)"), "값"))
  output$p_fdp <- renderPlot(lineplot(sim(), c("CFDP_", "HCOMP"),
      c("FDP (µg/mL)", "지혈 능력 H"), "값"))
  output$p_gain <- renderPlot({
    d <- sim() %>% mutate(gain = 1 / pmax(HCOMP, 1e-3))
    ggplot(d, aes(time, gain)) + geom_line(linewidth = 0.9, colour = "#553C7B") +
      geom_hline(yintercept = 1, linetype = 2, colour = "#8B0000") +
      labs(x = "시간 (일)", y = "1/H  (재출혈 증폭 배수)") + THEME
  })

  ## ---- 8. mechanics -------------------------------------------------
  output$p_icp <- renderPlot({
    p <- as.list(param(mod))
    lineplot(sim(), "ICP_", "ICP", "mmHg", hline = p$ICP0)
  })
  output$p_mls  <- renderPlot(lineplot(sim(), "MLS_", "정중선 전위", "mm"))
  output$p_comp <- renderPlot(lineplot(sim(), c("VCOMP", "R_"),
      c("압박 결손 V_comp (mL)", "유효 예비 R (mL)"), "mL"))
  output$p_space <- renderPlot({
    d <- sim() %>% mutate(space = pmax(R_ - VHEM, 0))
    ggplot(d, aes(time, space)) +
      geom_area(fill = "#FFD8A8", colour = "#B8860B") +
      labs(x = "시간 (일)", y = "잔존 경막하 공간 (mL)") + THEME
  })

  ## ---- 9. endpoints --------------------------------------------------
  output$p_reop <- renderPlot(lineplot(sim(), "P_REOP", "재수술 확률", "확률"))
  output$p_fav  <- renderPlot(lineplot(sim(), "P_FAV", "양호한 결과", "확률"))
  output$p_symp <- renderPlot(lineplot(sim(), c("SSYMP", "SCOG", "NNEUR"),
      c("증상 부담", "인지 결손", "신경 온전성"), "값"))
  output$p_safety <- renderPlot(lineplot(sim(), c("GLU", "XINF", "XMYO", "P_THR"),
      c("혈당 (mmol/L)", "감염 부담", "근병증 부담", "혈전 확률"), "값"))

  ## ---- 10. scenarios --------------------------------------------------
  scn_data <- reactive({
    keys <- names(CSDH_SCENARIOS)
    sel  <- keys[substr(keys, 1, 3) %in% input$scn_pick]
    if (!length(sel)) return(NULL)
    bind_rows(lapply(sel, function(k) {
      d <- CSDH_SCENARIOS[[k]]()
      data.frame(time = d$time, V = d$VHEM, dmax = d$dmax_,
                 reop = d$P_REOP, fav = d$P_FAV,
                 scenario = attr(d, "name"))
    }))
  })

  output$p_scn <- renderPlot({
    dd <- scn_data(); req(dd)
    ggplot(dd, aes(time, V, colour = scenario)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = rep(PAL, 3), name = NULL) +
      labs(x = "시간 (일)", y = "혈종 부피 (mL)") + THEME +
      guides(colour = guide_legend(ncol = 2))
  })

  output$t_scn <- renderDT({
    dd <- scn_data(); req(dd)
    dd %>% group_by(scenario) %>%
      summarise(`V 90일` = round(approx(time, V, 90, rule=2)$y, 1),
                `V 180일` = round(approx(time, V, 180, rule=2)$y, 1),
                `두께 180일 (mm)` = round(approx(time, dmax, 180, rule=2)$y, 1),
                `재수술 확률` = round(max(reop), 3),
                `양호한 결과` = round(min(fav), 3), .groups = "drop") %>%
      datatable(options = list(dom = "t", pageLength = 20), rownames = FALSE)
  })

  ## ---- 11. MMAE time course --------------------------------------------
  mmae_pair <- reactive({
    p <- list(V_RES = input$vres)
    list(a = csdh_sim("수술만", 180, surgery = 0, drain_days = 2, pars = p),
         b = csdh_sim("수술+MMAE", 180, surgery = 0, drain_days = 2,
                      embolise_at = 0, pars = p))
  })

  output$p_mmae <- renderPlot({
    z <- mmae_pair()
    dd <- bind_rows(
      data.frame(time = z$a$time, V = z$a$VHEM, NCAP = z$a$NCAP, arm = "수술만"),
      data.frame(time = z$b$time, V = z$b$VHEM, NCAP = z$b$NCAP, arm = "수술+MMAE")) %>%
      pivot_longer(c(V, NCAP), names_to = "v", values_to = "y") %>%
      mutate(v = factor(v, c("V", "NCAP"),
                        c("혈종 부피 (mL)", "미성숙 혈관 N_CAP")))
    ggplot(dd, aes(time, y, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~v, scales = "free_y") +
      geom_vline(xintercept = c(30, 90, 180), linetype = 3, colour = "#666666") +
      scale_colour_manual(values = c("#8B2635", "#276749"), name = NULL) +
      labs(x = "시간 (일)", y = NULL) + THEME
  })

  output$p_mmae_ard <- renderPlot({
    z <- mmae_pair()
    tt <- seq(7, 180, by = 1)
    ard <- approx(z$a$time, z$a$P_REOP, tt, rule=2)$y -
           approx(z$b$time, z$b$P_REOP, tt, rule=2)$y
    ggplot(data.frame(time = tt, ard = ard), aes(time, ard)) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_line(linewidth = 1, colour = "#276749") +
      geom_vline(xintercept = c(30, 90, 180), linetype = 3, colour = "#666666") +
      annotate("text", x = c(30, 90, 180), y = max(ard) * 0.15,
               label = c("30일", "90일\n(MAGIC-MT)", "180일\n(EMBOLISE)"),
               size = 3.1, hjust = -0.05) +
      labs(x = "시간 (일)", y = "절대 위험 감소 (재수술)") + THEME
  })

  ## ---- 12. ledger -------------------------------------------------------
  output$t_ledger <- renderDT({
    withProgress(message = "문헌 대조 계산 중...", value = 0.5, {
      rows <- CSDH_trial_ledger()
    })
    data.frame(엔드포인트 = sapply(rows, `[[`, 1),
               모델 = round(as.numeric(sapply(rows, `[[`, 2)), 3),
               발표값 = sapply(rows, `[[`, 3)) %>%
      datatable(options = list(dom = "t", pageLength = 20), rownames = FALSE)
  })
}

shinyApp(ui, server)
