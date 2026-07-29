# =============================================================================
#  ami_shiny_app.R
#  Acute Myocardial Infarction (STEMI) — QSP interactive dashboard
# =============================================================================
#
#  A front end for ami_mrgsolve_model.R.  The app is organised around the
#  model's thesis rather than around its variable list, so each tab answers one
#  question:
#
#    1  환자 · 병변      Who is this patient, and what did the artery do?
#    2  약동학 (PK)      What is in the blood, and when?
#    3  괴사 파면        Race 1 — does the wavefront beat the clock?
#    4  재관류 손상      Race 2 — acid washout versus the oxidant burst
#    5  미세혈관 폐색    Is the tissue reperfused, or only the artery?
#    6  염증 · 흉터      Does a scar get built before the wall thins?
#    7  재형성 분기      The bifurcation — converge or diverge?
#    8  임상 종말점      What would we have measured?
#    9  시나리오 비교    Head-to-head strategy comparison
#   10  바이오마커       Troponin, CK-MB, CRP, BNP, fibrinogen
#   11  민감도 · 구조    Cut one edge; see what the model was standing on
#
#  RUN
#    install.packages(c("shiny","mrgsolve","dplyr","tidyr","ggplot2",
#                       "DT","bslib","bsicons","scales"))
#    shiny::runApp("ami_shiny_app.R")
#
#  The app never displays a number the model did not compute.  Where a quantity
#  is a reporting convenience rather than a simulated state (the haemorrhage
#  index on tab 8), it is labelled as such in the UI, not only in the code.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)
library(bslib)
library(bsicons)
library(scales)

MODEL_FILE <- "ami_mrgsolve_model.R"
NOT_GIVEN  <- 1e6          # the model's "this intervention was withheld" value

mod <- mread_cache("ami", MODEL_FILE)

# -----------------------------------------------------------------------------
#  helpers
# -----------------------------------------------------------------------------
theme_ami <- function() {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(colour = "grey35", size = 11))
}

# palette: cool for flow/treatment, warm for injury, green for protection
PAL <- c(L1 = "#8B0000", L2 = "#C0392B", L3 = "#E67E22",
         L4 = "#7D9F35", L5 = "#1E8449",
         isch = "#C0392B", rep = "#7B241C", mvo = "#B8860B",
         edv = "#1A5276", ef = "#0E6251", scar = "#145A32")

#' Turn UI inputs into a param() list.  Withheld interventions become 1e6
#' rather than NA, because the model reads them as times inside a tanh switch.
build_params <- function(input) {
  gv <- function(flag, value) if (isTRUE(flag)) value else NOT_GIVEN
  list(
    AAR      = input$aar,
    COLL     = input$coll,
    PRECOND  = as.numeric(input$precond),
    AGE      = input$age,
    SBP0     = input$sbp0,
    HR0      = input$hr0,
    T_OCC    = 0,
    T_PCI    = gv(input$do_pci,  input$t_pci  / 60),
    T_LYT    = gv(input$do_lyt,  input$t_lyt  / 60),
    TNK_DOSE = input$tnk_dose,
    T_ASA    = gv(input$do_asa,  input$t_anti / 60),
    T_TIC    = gv(input$do_tic,  input$t_anti / 60),
    T_HEP    = gv(input$do_hep,  input$t_anti / 60),
    T_MET_IV = gv(input$do_ivbb, input$t_ivbb / 60),
    T_CSA    = gv(input$do_csa,  input$t_csa  / 60),
    T_MET_PO = gv(input$do_bb,   input$t_gdmt),
    T_RAM    = gv(input$do_acei, input$t_gdmt),
    T_EPL    = gv(input$do_mra,  input$t_gdmt + 24),
    T_EMP    = gv(input$do_sglt, input$t_gdmt),
    T_COLC   = gv(input$do_colc, input$t_colc),
    T_CAN    = gv(input$do_can,  input$t_can),
    ARNI     = as.numeric(input$arni)
  )
}

#' Simulate.  `delta` is deliberately fine early: the mPTP transient on tab 4
#' lives inside the first thirty minutes of reflow and a coarse grid hides it.
simulate <- function(pars, end_h, over = list()) {
  p <- modifyList(pars, over)
  fine <- mod %>% param(p) %>% mrgsim(end = min(end_h, 12), delta = 0.005) %>% as_tibble()
  if (end_h <= 12) return(fine)
  coarse <- mod %>% param(p) %>% mrgsim(end = end_h, delta = 0.25) %>% as_tibble()
  bind_rows(fine, filter(coarse, time > 12)) %>% arrange(time)
}

at_time <- function(df, h) df %>% slice(which.min(abs(time - h)))

fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)

# -----------------------------------------------------------------------------
#  UI
# -----------------------------------------------------------------------------
ui <- page_sidebar(
  title = "급성 심근경색 (STEMI) QSP — 두 개의 경주와 하나의 분기",
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  sidebar = sidebar(
    width = 360,
    accordion(
      open = c("환자 · 병변", "재관류 전략"),

      accordion_panel(
        "환자 · 병변",
        sliderInput("aar", "위험 영역 AAR (%LV)", 8, 55, 35, step = 1,
                    post = " %"),
        sliderInput("coll", "부수혈류 COLL (정상 혈류 대비)",
                    0.01, 0.35, 0.10, step = 0.01),
        helpText("부수혈류가 파면의 '시계 속도'를 정합니다. ",
                 "기초 대사 요구(0.20)를 넘는 층은 무운동이 되어도 생존합니다."),
        sliderInput("age", "연령 (년)", 35, 92, 62, step = 1),
        sliderInput("sbp0", "기저 수축기압 (mmHg)", 90, 180, 120, step = 5),
        sliderInput("hr0", "기저 심박수 (1/min)", 50, 110, 72, step = 2),
        checkboxInput("precond", "경색전 협심증 (허혈 전조)", FALSE)
      ),

      accordion_panel(
        "재관류 전략",
        checkboxInput("do_pci", "일차적/구조적 PCI", TRUE),
        sliderInput("t_pci", "PCI 시각 (증상 발현 후 분)", 20, 720, 90, step = 5),
        checkboxInput("do_lyt", "섬유소용해 (테넥테플라제)", FALSE),
        sliderInput("t_lyt", "용해제 볼루스 시각 (분)", 20, 360, 45, step = 5),
        sliderInput("tnk_dose", "테넥테플라제 용량 (mg)", 10, 50, 40, step = 5),
        helpText("두 전략을 동시에 켜면 약물–침습 병용 전략(STREAM형)이 됩니다.")
      ),

      accordion_panel(
        "항혈소판 · 항응고",
        sliderInput("t_anti", "투여 시각 (분)", 5, 360, 60, step = 5),
        checkboxInput("do_asa", "아스피린", TRUE),
        checkboxInput("do_tic", "P2Y12 억제제 (티카그렐러)", TRUE),
        checkboxInput("do_hep", "항응고 (헤파린)", TRUE)
      ),

      accordion_panel(
        "심근보호 (창은 분 단위)",
        checkboxInput("do_ivbb", "정맥 메토프롤롤", FALSE),
        sliderInput("t_ivbb", "정맥 베타차단제 시각 (분)", 10, 360, 40, step = 5),
        checkboxInput("do_csa", "mPTP 억제제 (사이클로스포린)", FALSE),
        sliderInput("t_csa", "CsA 시각 (분)", 10, 360, 85, step = 5),
        helpText("CsA는 재관류 직전에만 의미가 있습니다 — 모델에 시간창 ",
                 "파라미터는 없고, 공극이 몇 분 안에 열리기 때문입니다.")
      ),

      accordion_panel(
        "장기 약물 (GDMT)",
        sliderInput("t_gdmt", "GDMT 시작 (시간)", 6, 168, 24, step = 6),
        checkboxInput("do_acei", "ACE 억제제 (라미프릴)", TRUE),
        checkboxInput("arni", "ARNI로 대체 (사쿠비트릴/발사르탄)", FALSE),
        checkboxInput("do_bb", "경구 베타차단제", TRUE),
        checkboxInput("do_mra", "MR 길항제 (에플레레논)", TRUE),
        checkboxInput("do_sglt", "SGLT2 억제제 (엠파글리플로진)", TRUE),
        checkboxInput("do_colc", "콜히친", FALSE),
        sliderInput("t_colc", "콜히친 시작 (시간)", 6, 336, 24, step = 6),
        checkboxInput("do_can", "항IL-1β (카나키누맙)", FALSE),
        sliderInput("t_can", "항IL-1β 시각 (시간)", 1, 336, 24, step = 1)
      ),

      accordion_panel(
        "시뮬레이션",
        radioButtons("horizon", "시간 지평",
                     c("48시간 (급성)" = 48, "30일" = 720,
                       "180일" = 4320, "1년" = 8760),
                     selected = 4320),
        actionButton("run", "실행", class = "btn-primary w-100")
      )
    )
  ),

  navset_card_tab(
    id = "tabs",

    # ---- 1 ------------------------------------------------------------------
    nav_panel(
      "1 · 환자 · 병변",
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("층별 부수혈류 대 기초 대사 요구"),
             plotOutput("p_supply", height = 330),
             card_footer("수평선이 기초 요구입니다. 그 아래의 층만 죽습니다 — ",
                         "이 한 장이 모델의 생존 판정 전부입니다.")),
        card(card_header("심외막 개통 (PAT)"), plotOutput("p_patency", height = 330))
      ),
      card(card_header("이 실행의 요약"), DTOutput("t_summary"))
    ),

    # ---- 2 ------------------------------------------------------------------
    nav_panel(
      "2 · 약동학",
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("혈전용해 · 항혈소판"), plotOutput("p_pk_acute", height = 330)),
        card(card_header("장기 약물 효과 점유율"), plotOutput("p_pk_chronic", height = 330))
      ),
      card(card_header("섬유소용해계 — 플라스민 · 플라스미노겐 · 피브리노겐"),
           plotOutput("p_lysis", height = 300),
           card_footer("피브리노겐 저점이 전신 출혈 비용의 대리지표입니다."))
    ),

    # ---- 3 ------------------------------------------------------------------
    nav_panel(
      "3 · 괴사 파면 (Race 1)",
      layout_columns(
        col_widths = c(7, 5),
        card(card_header("층별 괴사 — 심내막에서 심외막으로"),
             plotOutput("p_wavefront", height = 380),
             card_footer("층 간의 유일한 차이는 부수혈류입니다. ",
                         "진행 순서를 모델에 알려 준 곳은 없습니다.")),
        card(card_header("에너지 충전도 E"), plotOutput("p_energy", height = 380))
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("경색 크기와 구제 지수"), plotOutput("p_is", height = 300)),
        card(card_header("재관류 시각 스윕 — 시간–근육 곡선"),
             plotOutput("p_salvage_curve", height = 300),
             card_footer("현재 설정의 나머지를 고정하고 PCI 시각만 바꾼 결과."))
      )
    ),

    # ---- 4 ------------------------------------------------------------------
    nav_panel(
      "4 · 재관류 손상 (Race 2)",
      card(card_header("재관류 직후 몇 분 — 산(H⁺)·숙신산·ROS·Ca²⁺·mPTP"),
           plotOutput("p_reperf", height = 400),
           card_footer("산이 세척되면서 공극의 관문이 열리고, 같은 순간 ",
                       "숙신산이 연소되어 산화제가 폭발합니다. ",
                       "이것이 pH 역설이며, 시간창 파라미터는 없습니다.")),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("허혈성 대 재관류성 괴사"), plotOutput("p_split", height = 320)),
        card(card_header("산성 관문과 공극 개방"), plotOutput("p_gate", height = 320))
      )
    ),

    # ---- 5 ------------------------------------------------------------------
    nav_panel(
      "5 · 미세혈관 폐색",
      layout_columns(
        col_widths = c(7, 5),
        card(card_header("MVO와 유효 층 혈류 — 열린 동맥, 닫힌 조직"),
             plotOutput("p_mvo", height = 380),
             card_footer("MVO는 혈류를 줄이고, 줄어든 혈류가 MVO를 키웁니다. ",
                         "no-reflow는 정귀환입니다.")),
        card(card_header("MVO 대 경색 크기"), plotOutput("p_mvo_is", height = 380))
      )
    ),

    # ---- 6 ------------------------------------------------------------------
    nav_panel(
      "6 · 염증 · 흉터",
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("세포 — 호중구 · M1 · M2"), plotOutput("p_cells", height = 330)),
        card(card_header("사이토카인 — IL-1β · IL-6 · TGF-β · CRP"),
             plotOutput("p_cyto", height = 330))
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("흉터 강도 대 MMP 활성 — 경쟁"),
             plotOutput("p_scar", height = 330),
             card_footer("흉터가 늦으면 벽이 얇아집니다. ",
                         "항염을 너무 강하게 하면 이 경쟁에서 지게 됩니다.")),
        card(card_header("경색벽 박화와 두께"), plotOutput("p_thin", height = 330))
      )
    ),

    # ---- 7 ------------------------------------------------------------------
    nav_panel(
      "7 · 재형성 분기",
      card(card_header("Laplace 고리 — 벽응력 · 용적 · 질량"),
           plotOutput("p_remodel", height = 380),
           card_footer("확장은 정귀환, 구심성 비후는 부귀환입니다. ",
                       "수렴이냐 발산이냐는 용량이 아니라 분기입니다.")),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("EDV · ESV · EF 궤적"), plotOutput("p_volumes", height = 330)),
        card(card_header("위험 영역 스윕 — 분리선은 어디인가"),
             plotOutput("p_bifurc", height = 330),
             card_footer("AAR만 바꾸어 1년 EDV 성장률을 계산한 것. ",
                         "임계 경색 크기라는 파라미터는 모델에 없습니다."))
      )
    ),

    # ---- 8 ------------------------------------------------------------------
    nav_panel(
      "8 · 임상 종말점",
      layout_columns(
        col_widths = c(4, 4, 4),
        value_box(title = "경색 크기 (%LV)", value = textOutput("vb_is"),
                  showcase = bsicons::bs_icon("heart-pulse")),
        value_box(title = "구제 지수 (%)", value = textOutput("vb_salv"),
                  showcase = bsicons::bs_icon("shield-check")),
        value_box(title = "MVO (%LV)", value = textOutput("vb_mvo"),
                  showcase = bsicons::bs_icon("droplet"))
      ),
      layout_columns(
        col_widths = c(4, 4, 4),
        value_box(title = "EF, 마지막 시점 (%)", value = textOutput("vb_ef"),
                  showcase = bsicons::bs_icon("activity")),
        value_box(title = "EDV, 마지막 시점 (mL)", value = textOutput("vb_edv"),
                  showcase = bsicons::bs_icon("arrows-expand")),
        value_box(title = "흉터 강도 (0–1)", value = textOutput("vb_scar"),
                  showcase = bsicons::bs_icon("bandaid"))
      ),
      card(card_header("혈역학과 신경호르몬"), plotOutput("p_neuro", height = 330)),
      card(card_header("보고용 지표 (검증된 위험 모형 아님)"),
           DTOutput("t_report"),
           card_footer(strong("주의: "), "두개내 출혈 지수는 플라스민 노출과 ",
                       "연령만으로 만든 보고용 수치이며 절대 위험도로 읽어서는 ",
                       "안 됩니다. 모델은 사건 발생률을 계산하지 않습니다."))
    ),

    # ---- 9 ------------------------------------------------------------------
    nav_panel(
      "9 · 시나리오 비교",
      card(card_header("전략 비교 — 같은 환자, 다른 결정"),
           checkboxGroupInput(
             "scen", NULL, inline = TRUE,
             choices = c("재관류 없음" = "none",
                         "PCI 90분" = "pci90",
                         "PCI 240분" = "pci240",
                         "병원전 용해 45분" = "lyt45",
                         "용해 45분 + PCI 4시간" = "pharminv",
                         "PCI 90분 + 전체 GDMT" = "pci90gdmt",
                         "PCI 90분 + GDMT + 심근보호" = "full"),
             selected = c("none", "pci90", "pci90gdmt", "full")),
           plotOutput("p_scen", height = 420)),
      card(card_header("시나리오 종말점 표"), DTOutput("t_scen"))
    ),

    # ---- 10 -----------------------------------------------------------------
    nav_panel(
      "10 · 바이오마커",
      card(card_header("cTnI와 CK-MB — 세척 현상"),
           plotOutput("p_markers", height = 360),
           card_footer("재관류된 동맥은 조기·고봉 곡선을, 폐색된 동맥은 ",
                       "더 큰 실제 경색에서 늦고 낮은 곡선을 냅니다. ",
                       "그래서 효소로 경색을 재려면 곡선 아래 면적이 필요합니다.")),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("염증 · 벽응력 표지자"), plotOutput("p_bio2", height = 320)),
        card(card_header("정점 대 AUC 대 실제 경색 크기"), DTOutput("t_markers"))
      )
    ),

    # ---- 11 -----------------------------------------------------------------
    nav_panel(
      "11 · 민감도 · 구조",
      card(card_header("한 번에 하나의 기전만 끊기"),
           helpText("각 행은 파라미터 하나를 0 또는 무한대로 보내고 같은 ",
                    "시나리오를 다시 돌린 결과입니다. 모델이 무엇을 딛고 ",
                    "서 있었는지 보여 줍니다."),
           DTOutput("t_sens")),
      card(card_header("구조적 절단의 효과"), plotOutput("p_sens", height = 380))
    )
  )
)

# -----------------------------------------------------------------------------
#  server
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  pars <- eventReactive(input$run, build_params(input), ignoreNULL = FALSE)
  horizon <- eventReactive(input$run, as.numeric(input$horizon), ignoreNULL = FALSE)

  sim <- reactive({
    p <- pars(); p$AAR <- p$AAR / 100
    simulate(p, horizon())
  })

  acute <- reactive(filter(sim(), time <= 48))
  last_row <- reactive(sim() %>% slice_tail(n = 1))

  # ---- 1 --------------------------------------------------------------------
  output$p_supply <- renderPlot({
    p <- pars()
    gc <- c(0.05, 0.35, 0.75, 1.35, 2.50)
    d <- tibble(layer = factor(paste0("L", 1:5), levels = paste0("L", 1:5)),
                collateral = p$COLL * gc)
    ggplot(d, aes(layer, collateral, fill = layer)) +
      geom_col(width = .65) +
      geom_hline(yintercept = 0.20, linetype = "dashed",
                 colour = "#8B0000", linewidth = 1) +
      annotate("text", x = 4.2, y = 0.215, label = "기초 대사 요구 0.20",
               colour = "#8B0000", size = 4) +
      scale_fill_manual(values = PAL[1:5], guide = "none") +
      labs(x = NULL, y = "부수혈류 (정상 대비)",
           title = "죽는 층과 사는 층",
           subtitle = "점선 아래 = 기초 요구를 못 채움 = 사망 경로") +
      theme_ami()
  })

  output$p_patency <- renderPlot({
    ggplot(acute(), aes(time * 60, PATENCY)) +
      geom_line(colour = PAL[["edv"]], linewidth = 1.1) +
      scale_x_continuous(limits = c(0, 480)) +
      labs(x = "증상 발현 후 (분)", y = "심외막 개통 (%)",
           title = "동맥은 언제 열렸는가") +
      theme_ami()
  })

  output$t_summary <- renderDT({
    a <- at_time(sim(), 48); z <- last_row()
    tibble(
      항목 = c("위험 영역 (%LV)", "부수혈류", "경색 크기 48h (%LV)",
               "구제 지수 (%)", "허혈성 성분 (%LV)", "재관류성 성분 (%LV)",
               "재관류 손상 분율 (%)", "MVO (%LV)",
               "EF 48h (%)", "EF 마지막 (%)", "EDV 마지막 (mL)",
               "흉터 강도", "cTnI 정점"),
      값 = c(fmt(input$aar, 0), fmt(pars()$COLL, 3), fmt(a$IS), fmt(a$SALV, 1),
             fmt(a$IS_ISCH), fmt(a$IS_REP), fmt(a$RPFRAC, 1), fmt(a$MVO_LV),
             fmt(a$EF, 1), fmt(z$EF, 1), fmt(z$EDV, 1), fmt(z$SCAR, 3),
             fmt(max(sim()$CTNI), 1))
    ) %>% datatable(rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  # ---- 2 --------------------------------------------------------------------
  output$p_pk_acute <- renderPlot({
    acute() %>%
      transmute(time, `테넥테플라제 (mg/L)` = TNK1 / 4.2,
                `티카그렐러 (mg/L)` = TICc / 88,
                `메토프롤롤 (mg/L)` = METc / 250) %>%
      pivot_longer(-time) %>% filter(value > 1e-9) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = unname(PAL[c(1, 4, 5)]), guide = "none") +
      labs(x = "시간 (h)", y = "혈중 농도") + theme_ami()
  })

  output$p_pk_chronic <- renderPlot({
    sim() %>%
      transmute(time = time / 24,
                ACEi = RAMc, MRA = EPL, SGLT2i = EMP,
                `콜히친` = COLC, `항IL-1β` = CAN) %>%
      pivot_longer(-time) %>% group_by(name) %>%
      filter(max(value) > 1e-9) %>% mutate(value = value / max(value)) %>%
      ungroup() %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      labs(x = "시간 (일)", y = "효과구획 (최대 대비)", colour = NULL,
           title = "장기 약물 노출") + theme_ami()
  })

  output$p_lysis <- renderPlot({
    acute() %>% filter(time <= 12) %>%
      transmute(time, `플라스민` = PLN, `플라스미노겐` = PLG,
                `피브리노겐` = FIBX, `PAI-1` = PAI) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      labs(x = "시간 (h)", y = "정상 대비 / 활성", colour = NULL) + theme_ami()
  })

  # ---- 3 --------------------------------------------------------------------
  output$p_wavefront <- renderPlot({
    acute() %>%
      select(time, NECRO1:NECRO5) %>%
      pivot_longer(-time, names_to = "layer") %>%
      mutate(layer = factor(sub("NECRO", "L", layer),
                            levels = paste0("L", 1:5))) %>%
      ggplot(aes(time, value, colour = layer)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = PAL[1:5],
                          labels = c("L1 심내막하", "L2", "L3 중간벽",
                                     "L4", "L5 심외막하")) +
      labs(x = "시간 (h)", y = "층 내 괴사 (%)", colour = NULL,
           title = "괴사 파면", subtitle = "층 간 차이는 부수혈류 하나뿐") +
      theme_ami()
  })

  output$p_energy <- renderPlot({
    acute() %>% filter(time <= 12) %>%
      select(time, E1, E2, E3, E4, E5) %>%
      pivot_longer(-time, names_to = "layer") %>%
      mutate(layer = factor(sub("E", "L", layer), levels = paste0("L", 1:5))) %>%
      ggplot(aes(time, value, colour = layer)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = PAL[1:5], guide = "none") +
      labs(x = "시간 (h)", y = "에너지 충전도 E",
           title = "에너지는 기초 요구를 못 낼 때만 떨어진다") + theme_ami()
  })

  output$p_is <- renderPlot({
    acute() %>%
      ggplot(aes(time)) +
      geom_line(aes(y = IS, colour = "경색 크기 (%LV)"), linewidth = 1.2) +
      geom_line(aes(y = SALV * input$aar / 100, colour = "구제 (%LV)"),
                linewidth = 1) +
      scale_colour_manual(values = c(unname(PAL[["isch"]]), unname(PAL[["scar"]]))) +
      labs(x = "시간 (h)", y = "%LV", colour = NULL) + theme_ami()
  })

  salvage_curve <- reactive({
    p <- pars(); p$AAR <- p$AAR / 100
    tp <- c(15, 30, 45, 60, 90, 120, 180, 240, 360, 480)
    lapply(tp, function(m) {
      s <- simulate(p, 48, over = list(T_PCI = m / 60, T_ASA = 0.05,
                                       T_TIC = 0.05, T_HEP = 0.05,
                                       T_LYT = NOT_GIVEN)) %>% slice_tail(n = 1)
      tibble(delay = m, IS = s$IS, SALV = s$SALV, RI = s$RPFRAC, MVO = s$MVO_LV)
    }) %>% bind_rows()
  })

  output$p_salvage_curve <- renderPlot({
    salvage_curve() %>%
      ggplot(aes(delay, IS)) +
      geom_line(colour = PAL[["isch"]], linewidth = 1.2) +
      geom_point(size = 2.4, colour = PAL[["isch"]]) +
      labs(x = "재관류 지연 (분)", y = "최종 경색 크기 (%LV)",
           title = "시간–근육 곡선은 직선이 아니다") + theme_ami()
  })

  # ---- 4 --------------------------------------------------------------------
  output$p_reperf <- renderPlot({
    trep <- if (isTRUE(input$do_pci)) input$t_pci / 60 else input$t_lyt / 60
    acute() %>%
      filter(time >= max(trep - 0.25, 0), time <= trep + 2) %>%
      transmute(min = (time - trep) * 60,
                `H⁺ (L1)` = H1, `숙신산 (L1)` = SU1, `ROS` = ROS,
                `Ca²⁺ (L1)` = C1, `mPTP (L1)` = P1 * 10) %>%
      pivot_longer(-min) %>%
      ggplot(aes(min, value, colour = name)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_line(linewidth = 1.1) +
      annotate("text", x = 2, y = Inf, vjust = 1.5, hjust = 0,
               label = "재관류", colour = "grey30") +
      labs(x = "재관류 시점 기준 (분)", y = "값 (mPTP는 ×10)", colour = NULL,
           title = "재관류 후 첫 몇 분",
           subtitle = "산이 씻겨 나가는 순간 공극의 관문이 열린다") +
      theme_ami()
  })

  output$p_split <- renderPlot({
    acute() %>% select(time, `허혈성` = IS_ISCH, `재관류성` = IS_REP) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, fill = name)) +
      geom_area(alpha = .85) +
      scale_fill_manual(values = unname(PAL[c("isch", "rep")])) +
      labs(x = "시간 (h)", y = "%LV", fill = NULL,
           title = "두 기전의 분해 — 가정이 아니라 결과") + theme_ami()
  })

  output$p_gate <- renderPlot({
    acute() %>% filter(time <= 12) %>%
      transmute(time, `산성 관문 1/(1+H/0.3)` = 1 / (1 + H1 / 0.30),
                `mPTP 개방 (L1)` = P1) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c(unname(PAL[["scar"]]), unname(PAL[["rep"]]))) +
      labs(x = "시간 (h)", y = "0–1", colour = NULL) + theme_ami()
  })

  # ---- 5 --------------------------------------------------------------------
  output$p_mvo <- renderPlot({
    acute() %>% transmute(time, `MVO (%LV)` = MVO_LV,
                          `심외막 개통 (%)` = PATENCY) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c(unname(PAL[["edv"]]), unname(PAL[["mvo"]]))) +
      labs(x = "시간 (h)", y = NULL, colour = NULL,
           title = "동맥은 열렸다. 조직은?") + theme_ami()
  })

  output$p_mvo_is <- renderPlot({
    salvage_curve() %>%
      ggplot(aes(IS, MVO, colour = delay)) +
      geom_path(linewidth = 1) + geom_point(size = 3) +
      scale_colour_viridis_c(option = "plasma") +
      labs(x = "경색 크기 (%LV)", y = "MVO (%LV)", colour = "지연 (분)",
           title = "MVO는 경색 크기의 함수가 아니다") + theme_ami()
  })

  # ---- 6 --------------------------------------------------------------------
  output$p_cells <- renderPlot({
    sim() %>% filter(time <= 24 * 21) %>%
      transmute(day = time / 24, `호중구` = NEU, `M1` = M1, `M2` = M2) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1.1) +
      labs(x = "일", y = "상대 수준", colour = NULL,
           title = "염증에서 수복으로") + theme_ami()
  })

  output$p_cyto <- renderPlot({
    sim() %>% filter(time <= 24 * 21) %>%
      transmute(day = time / 24, `IL-1β` = IL1, `IL-6` = IL6,
                `TGF-β` = TGF, `CRP` = CRPX) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1.1) +
      labs(x = "일", y = "상대 수준", colour = NULL) + theme_ami()
  })

  output$p_scar <- renderPlot({
    sim() %>% filter(time <= 24 * 60) %>%
      transmute(day = time / 24, `흉터 강도` = SCAR, `MMP 활성` = MMP) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1.2) +
      scale_colour_manual(values = c(unname(PAL[["mvo"]]), unname(PAL[["scar"]]))) +
      labs(x = "일", y = NULL, colour = NULL,
           title = "흉터를 짓는 속도 대 허무는 속도") + theme_ami()
  })

  output$p_thin <- renderPlot({
    sim() %>% transmute(day = time / 24, `박화 (%)` = THINPCT,
                        `경색부 두께 (cm)` = H_INF * 10) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1.1) +
      labs(x = "일", y = "값 (두께는 mm)", colour = NULL) + theme_ami()
  })

  # ---- 7 --------------------------------------------------------------------
  output$p_remodel <- renderPlot({
    sim() %>% transmute(day = time / 24,
                        `벽응력 (원격)` = WS / first(WS),
                        `벽응력 (경색부)` = WS_INF / first(WS),
                        `EDV / EDV₀` = EDV / first(EDV),
                        `질량 / 질량₀` = MASS / first(MASS)) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) +
      geom_hline(yintercept = 1, linetype = "dotted") +
      geom_line(linewidth = 1.1) +
      labs(x = "일", y = "기저 대비", colour = NULL,
           title = "Laplace 고리",
           subtitle = "확장은 응력을 올리고, 비후는 내린다") + theme_ami()
  })

  output$p_volumes <- renderPlot({
    sim() %>% transmute(day = time / 24, EDV, ESV, EF) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1.1) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = unname(PAL[c("edv", "ef", "isch")]),
                          guide = "none") +
      labs(x = "일", y = NULL) + theme_ami()
  })

  bifurc <- reactive({
    p <- pars()
    aars <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50)
    lapply(aars, function(a) {
      s <- simulate(p, 24 * 365, over = list(AAR = a))
      e30 <- at_time(s, 24 * 30)$EDV
      z <- slice_tail(s, n = 1)
      tibble(AAR = a * 100, IS = z$IS,
             growth = 100 * (z$EDV - e30) / e30, EF = z$EF)
    }) %>% bind_rows()
  })

  output$p_bifurc <- renderPlot({
    bifurc() %>%
      ggplot(aes(IS, growth)) +
      geom_hline(yintercept = c(4, 12), linetype = c("dotted", "dashed"),
                 colour = c("grey50", "#8B0000")) +
      geom_line(colour = PAL[["edv"]], linewidth = 1.2) +
      geom_point(size = 2.6, colour = PAL[["edv"]]) +
      annotate("text", x = Inf, y = 12, hjust = 1.1, vjust = -0.5,
               label = "발산 판정선", colour = "#8B0000", size = 3.6) +
      labs(x = "경색 크기 (%LV)", y = "30일→1년 EDV 성장 (%)",
           title = "분리선은 계산 결과다") + theme_ami()
  })

  # ---- 8 --------------------------------------------------------------------
  output$vb_is   <- renderText(fmt(at_time(sim(), 48)$IS))
  output$vb_salv <- renderText(fmt(at_time(sim(), 48)$SALV, 1))
  output$vb_mvo  <- renderText(fmt(max(sim()$MVO_LV)))
  output$vb_ef   <- renderText(fmt(last_row()$EF, 1))
  output$vb_edv  <- renderText(fmt(last_row()$EDV, 1))
  output$vb_scar <- renderText(fmt(last_row()$SCAR, 3))

  output$p_neuro <- renderPlot({
    sim() %>% transmute(day = time / 24, NE, `Ang II` = ANG, `알도스테론` = ALD,
                        `NP (BNP)` = NTPROBNP, `용적` = VOL) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 1, linetype = "dotted") +
      labs(x = "일", y = "정상 대비", colour = NULL) + theme_ami()
  })

  output$t_report <- renderDT({
    s <- sim()
    pln_auc <- sum(diff(s$time) * head(s$PLN, -1))
    ich <- 0.9 * pln_auc * (1 + 0.055 * (input$age - 50))
    tibble(
      지표 = c("플라스민 노출 (AUC)", "피브리노겐 저점 (정상 대비)",
               "두개내 출혈 지수 (보고용)", "cTnI AUC", "CRP 정점"),
      값 = c(fmt(pln_auc, 4), fmt(min(s$FIBX), 3), fmt(ich, 3),
             fmt(sum(diff(s$time) * head(s$CTNI, -1)), 0), fmt(max(s$CRPX), 2))
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # ---- 9 --------------------------------------------------------------------
  SCEN <- list(
    none      = list(lab = "재관류 없음",
                     ov = list(T_PCI = NOT_GIVEN, T_LYT = NOT_GIVEN)),
    pci90     = list(lab = "PCI 90분",
                     ov = list(T_PCI = 1.5, T_LYT = NOT_GIVEN, T_ASA = 1,
                               T_TIC = 1, T_HEP = 1)),
    pci240    = list(lab = "PCI 240분",
                     ov = list(T_PCI = 4, T_LYT = NOT_GIVEN, T_ASA = 3.5,
                               T_TIC = 3.5, T_HEP = 3.5)),
    lyt45     = list(lab = "병원전 용해 45분",
                     ov = list(T_PCI = NOT_GIVEN, T_LYT = 0.75, T_ASA = 0.65,
                               T_TIC = 0.65, T_HEP = 0.65)),
    pharminv  = list(lab = "용해 45분 + PCI 4시간",
                     ov = list(T_PCI = 4, T_LYT = 0.75, T_ASA = 0.65,
                               T_TIC = 0.65, T_HEP = 0.65)),
    pci90gdmt = list(lab = "PCI 90분 + 전체 GDMT",
                     ov = list(T_PCI = 1.5, T_LYT = NOT_GIVEN, T_ASA = 1,
                               T_TIC = 1, T_HEP = 1, T_RAM = 24,
                               T_MET_PO = 24, T_EPL = 48, T_EMP = 24)),
    full      = list(lab = "PCI 90분 + GDMT + 심근보호",
                     ov = list(T_PCI = 1.5, T_LYT = NOT_GIVEN, T_ASA = 1,
                               T_TIC = 1, T_HEP = 1, T_MET_IV = 0.6,
                               T_CSA = 1.45, T_COLC = 24, T_RAM = 24,
                               T_MET_PO = 24, T_EPL = 48, T_EMP = 24))
  )

  scen_runs <- reactive({
    req(length(input$scen) > 0)
    p <- pars(); p$AAR <- p$AAR / 100
    lapply(input$scen, function(k) {
      simulate(p, horizon(), over = SCEN[[k]]$ov) %>%
        mutate(scenario = SCEN[[k]]$lab)
    }) %>% bind_rows()
  })

  output$p_scen <- renderPlot({
    scen_runs() %>%
      transmute(day = time / 24, scenario,
                `경색 크기 (%LV)` = IS, `MVO (%LV)` = MVO_LV,
                `EF (%)` = EF, `EDV (mL)` = EDV) %>%
      pivot_longer(-c(day, scenario)) %>%
      ggplot(aes(day, value, colour = scenario)) +
      geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "일", y = NULL, colour = NULL) + theme_ami()
  })

  output$t_scen <- renderDT({
    scen_runs() %>% group_by(scenario) %>%
      summarise(`경색 크기 (%LV)` = round(IS[which.min(abs(time - 48))], 2),
                `구제 지수 (%)` = round(SALV[which.min(abs(time - 48))], 1),
                `MVO (%LV)` = round(max(MVO_LV), 2),
                `cTnI 정점` = round(max(CTNI), 1),
                `EF 마지막 (%)` = round(last(EF), 1),
                `EDV 마지막 (mL)` = round(last(EDV), 1),
                `흉터 강도` = round(last(SCAR), 3),
                .groups = "drop") %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # ---- 10 -----------------------------------------------------------------
  output$p_markers <- renderPlot({
    p <- pars(); p$AAR <- p$AAR / 100
    variants <- list(
      list(lab = "PCI 60분", ov = list(T_PCI = 1, T_LYT = NOT_GIVEN)),
      list(lab = "PCI 240분", ov = list(T_PCI = 4, T_LYT = NOT_GIVEN)),
      list(lab = "재관류 없음", ov = list(T_PCI = NOT_GIVEN, T_LYT = NOT_GIVEN))
    )
    lapply(variants, function(v)
      simulate(p, 120, over = v$ov) %>% mutate(scenario = v$lab)) %>%
      bind_rows() %>%
      transmute(time, scenario, `cTnI` = CTNI, `CK-MB` = CKMBP) %>%
      pivot_longer(-c(time, scenario)) %>%
      ggplot(aes(time, value, colour = scenario)) +
      geom_line(linewidth = 1.1) + facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (h)", y = "혈중 농도 (임의 단위)", colour = NULL,
           title = "세척 현상",
           subtitle = "재관류된 동맥 = 조기·고봉. 폐색된 동맥 = 늦고 낮되 경색은 더 크다") +
      theme_ami()
  })

  output$p_bio2 <- renderPlot({
    sim() %>% transmute(day = time / 24, `CRP` = CRPX,
                        `NP (BNP)` = NTPROBNP, `피브리노겐` = FIBX) %>%
      pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 1.1) +
      labs(x = "일", y = "정상 대비", colour = NULL) + theme_ami()
  })

  output$t_markers <- renderDT({
    p <- pars(); p$AAR <- p$AAR / 100
    rows <- list(
      list(lab = "PCI 60분", ov = list(T_PCI = 1, T_LYT = NOT_GIVEN)),
      list(lab = "PCI 120분", ov = list(T_PCI = 2, T_LYT = NOT_GIVEN)),
      list(lab = "PCI 240분", ov = list(T_PCI = 4, T_LYT = NOT_GIVEN)),
      list(lab = "재관류 없음", ov = list(T_PCI = NOT_GIVEN, T_LYT = NOT_GIVEN))
    )
    lapply(rows, function(v) {
      s <- simulate(p, 120, over = v$ov)
      auc <- sum(diff(s$time) * head(s$CTNI, -1))
      tibble(시나리오 = v$lab,
             `TnI 정점` = round(max(s$CTNI), 1),
             `정점 시각 (h)` = round(s$time[which.max(s$CTNI)], 2),
             `TnI AUC` = round(auc, 0),
             `실제 경색 (%LV)` = round(last(s$IS), 2),
             `AUC / 경색` = round(auc / max(last(s$IS), 1e-6), 0))
    }) %>% bind_rows() %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # ---- 11 -----------------------------------------------------------------
  CUTS <- list(
    list(lab = "(온전한 모델)", ov = list()),
    list(lab = "mPTP 없음 (KP_ON = 0)", ov = list(KP_ON = 0)),
    list(lab = "재관류 Ca 급증 없음 (KCA_RP = 0)", ov = list(KCA_RP = 0)),
    list(lab = "산성 관문 제거 (HP50 = 1e6)", ov = list(HP50 = 1e6)),
    list(lab = "숙신산 없음 (KSUC = 0)", ov = list(KSUC = 0)),
    list(lab = "산화제 폭발 없음 (KROS = 0)", ov = list(KROS = 0)),
    list(lab = "부수혈류 기울기 평탄화", ov = list(GC1 = 1, GC2 = 1, GC3 = 1,
                                                  GC4 = 1, GC5 = 1)),
    list(lab = "부수혈류 없음 (COLL = 0)", ov = list(COLL = 0)),
    list(lab = "MVO 귀환 없음", ov = list(GM1 = 0, GM2 = 0, GM3 = 0,
                                          GM4 = 0, GM5 = 0)),
    list(lab = "확장 고리 없음 (KDIL = 0)", ov = list(KDIL = 0)),
    list(lab = "비후 제동 없음 (KHYP = 0)", ov = list(KHYP = 0)),
    list(lab = "혐기 예비 없음 (F_ANAERO = 0)", ov = list(F_ANAERO = 0)),
    list(lab = "허혈 전조 숙주 (PRECOND = 1)", ov = list(PRECOND = 1))
  )

  sens <- reactive({
    p <- pars(); p$AAR <- p$AAR / 100
    lapply(CUTS, function(c) {
      s <- simulate(p, max(horizon(), 24 * 180), over = c$ov)
      a <- at_time(s, 48); z <- slice_tail(s, n = 1)
      tibble(perturbation = c$lab, IS = a$IS, IS_REP = a$IS_REP,
             MVO = max(s$MVO_LV), EDV = z$EDV, EF = z$EF)
    }) %>% bind_rows()
  })

  output$t_sens <- renderDT({
    sens() %>%
      transmute(`절단` = perturbation,
                `경색 크기 (%LV)` = round(IS, 2),
                `재관류 성분 (%LV)` = round(IS_REP, 2),
                `MVO (%LV)` = round(MVO, 2),
                `EDV 180일 (mL)` = round(EDV, 1),
                `EF 180일 (%)` = round(EF, 1)) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 15))
  })

  output$p_sens <- renderPlot({
    base <- sens() %>% filter(perturbation == "(온전한 모델)")
    sens() %>% filter(perturbation != "(온전한 모델)") %>%
      mutate(dIS = IS - base$IS, dEDV = EDV - base$EDV) %>%
      pivot_longer(c(dIS, dEDV)) %>%
      mutate(name = recode(name, dIS = "경색 크기 변화 (%LV)",
                           dEDV = "180일 EDV 변화 (mL)")) %>%
      ggplot(aes(reorder(perturbation, value), value, fill = value > 0)) +
      geom_col() + coord_flip() +
      facet_wrap(~name, scales = "free_x") +
      scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#1E8449"),
                        guide = "none") +
      labs(x = NULL, y = "온전한 모델 대비 변화",
           title = "모델은 무엇을 딛고 서 있었는가") + theme_ami()
  })
}

shinyApp(ui, server)
