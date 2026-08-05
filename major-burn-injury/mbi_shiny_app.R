## =====================================================================
##  mbi_shiny_app.R
##  MAJOR THERMAL BURN INJURY — QSP interactive dashboard
##  =====================================================================
##
##  The app is organised around the model's three statements rather than
##  around organ systems, because the statements are the model:
##
##   (1) THE PRODUCT      dVP = V_infused x f_ret(Pi_p), and Pi_p is
##                        CONVEX in protein, so f_ret decays as you work.
##   (2) THE CONTROLLER   UO-titrated resuscitation is a closed loop whose
##                        gain IS f_ret.  Fluid creep is its fixed point.
##   (3) THE TWO CLOCKS   leak tau ~ 9.5 h vs hypermetabolism tau ~ 104 d,
##                        and the slow clock is driven by A_open(t).
##
##  Tabs
##   1  Patient & protocol      — set the case and the resuscitation plan
##   2  The convex function     — COP vs protein, and where you are on it
##   3  Fluid & Starling        — compartments, fluxes, retention fraction
##   4  The controller          — rate, urine output, loop gain, creep
##   5  Resuscitation morbidity — IAP, ACS, lung water, weight gain
##   6  Hypermetabolism         — REE, catecholamines, core temperature
##   7  Body composition        — lean mass, the S-B difference, bone, fat
##   8  Wound & closure         — A_open, grafts, donor sites, scar
##   9  Inflammation & infection— cytokines, HLA-DR, colony count, sepsis
##  10  Pharmacology & PK       — propranolol, oxandrolone, vancomycin/ARC
##  11  Endpoints & validation  — the calibration table and rBaux
##  12  Scenario comparison     — any subset of the 22 scenarios overlaid
##
##  Run:
##      source("mbi_mrgsolve_model.R")   # defines mod, SCEN, run_scenario
##      shiny::runApp("mbi_shiny_app.R")
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

if (!exists("mod") || !exists("SCEN")) {
  source("mbi_mrgsolve_model.R")
}

THEME <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom", legend.title = element_blank())

LP <- function(C) 2.1*C + 0.16*C^2 + 0.009*C^3

TARGETS <- tibble::tribble(
  ~metric,                                   ~lit,                     ~lo,   ~hi,
  "24-h volume (mL/kg/%TBSA)",               "5.2-6.7",                5.2,   6.7,
  "in:out ratio vs Parkland",                "1.2-1.6",                1.2,   1.6,
  "plasma volume nadir (%)",                 "60-80",                 60,    80,
  "plasma COP at 24 h (mmHg)",               "10-16",                 10,    16,
  "peak weight gain (%)",                    "+15 to +30",            15,    30,
  "peak REE (% predicted)",                  "120-180",              120,   180,
  "lean mass at 14 d, no drug (%)",          "about -9",             -13,    -5,
  "24-h volume (mL/kg)",                     "IAH above ~250",       150,   250
)

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel("중증 화상 QSP 모델 — Major Thermal Burn Injury"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("한 번의 손상 · 두 개의 시계 · 그리고 <b>일하는 동안 이득이 줄어드는 제어기</b>")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("환자 (Patient)"),
      sliderInput("WT",   "체중 Weight (kg)",      30, 150, 80,  step = 5),
      sliderInput("AGE",  "나이 Age (years)",       1,  95, 35,  step = 1),
      sliderInput("TBSA", "화상 범위 %TBSA",        5,  95, 45,  step = 5),
      checkboxInput("INH", "흡입 손상 Inhalation injury", FALSE),

      hr(),
      h4("소생술 (Resuscitation)"),
      selectInput("formula", "시작 공식",
                  c("Parkland 4 mL/kg/%"        = 4,
                    "Modified Brooke 2 mL/kg/%" = 2,
                    "ISBI start-low 2 mL/kg/%"  = 2.0001,
                    "소생술 없음 (historical)"   = 0), selected = 4),
      checkboxInput("titrate", "소변량에 맞춘 적정 (closed loop)", TRUE),
      sliderInput("uotgt", "소변량 목표 (mL/kg/h)", 0.3, 1.5, 0.5, step = 0.1),
      sliderInput("gctrl", "제어기 이득 gain",      0.05, 0.9, 0.32, step = 0.01),
      sliderInput("rmax",  "투여속도 상한 (× 공식)", 1.5, 6.0, 3.0, step = 0.5),
      sliderInput("kopi",  "오피오이드에 의한 소변량 감소", 0, 0.6, 0, step = 0.05),

      hr(),
      h4("콜로이드 · 항산화 (Starling terms)"),
      sliderInput("fcol",  "5% 알부민 분율", 0, 0.6, 0, step = 0.05),
      sliderInput("albstart", "알부민 시작 (h)", 0, 24, 8, step = 1),
      checkboxInput("vitc", "고용량 아스코르브산 66 mg/kg/h × 24 h", FALSE),

      hr(),
      h4("수술 (Surgery)"),
      sliderInput("texc", "첫 가피절제 시점 (day)", 1, 21, 5, step = 1),
      sliderInput("excint", "수술 간격 (h)", 24, 168, 72, step = 24),
      sliderInput("topical", "국소 항균제 효과 (/h)", 0.000, 0.045, 0.022, step = 0.002),

      hr(),
      h4("약물 (Drugs)"),
      sliderInput("prop", "프로프라놀롤 (mg/kg/day)", 0, 6, 0, step = 0.5),
      sliderInput("oxa",  "옥산드롤론 (mg BID)",       0, 20, 0, step = 5),
      sliderInput("tamb", "실내 온도 (°C)",           20, 34, 31, step = 1),

      hr(),
      sliderInput("end", "관찰 기간 (days)", 3, 90, 60, step = 1),
      actionButton("go", "실행 Run", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1 · 환자·프로토콜",
                 br(), htmlOutput("summary_box"), br(),
                 h4("입원 시 요약"), tableOutput("patient_tbl"),
                 h4("이 설정이 무엇을 바꾸는가"),
                 tags$ul(
                   tags$li(HTML("<b>공식과 목표 소변량</b>은 제어기의 <i>설정점</i>만 바꿉니다. 실제 투여량은 누출이 결정합니다.")),
                   tags$li(HTML("<b>알부민 분율</b>은 부피를 더하는 것이 아니라 <b>f_ret를 복원</b>합니다.")),
                   tags$li(HTML("<b>아스코르브산</b>은 Kf 병변의 산화 성분에만 작용합니다.")),
                   tags$li(HTML("<b>가피절제 시점</b>은 느린 시계의 <b>구동원 자체</b>를 제거합니다.")),
                   tags$li(HTML("<b>프로프라놀롤</b>은 구동원이 아니라 <b>전달자</b>를 차단합니다 — 둘은 상가적이지 않습니다."))
                 )),

        tabPanel("2 · 볼록 함수",
                 br(), plotOutput("p_convex", height = 420),
                 HTML("<p>Landis–Pappenheimer: <code>Pi = 2.1C + 0.16C² + 0.009C³</code>.
                       단백질이 절반으로 희석되면 교질삼투압은 <b>62 %</b> 사라집니다.
                       빨간 점은 현재 시뮬레이션이 24시간 동안 이 곡선 위에서 어디로
                       내려가는지를 보여줍니다. 곡선이 볼록하기 때문에 <b>매 리터가
                       직전 리터보다 비쌉니다</b> — 이것이 수액 크립의 수학적 엔진입니다.</p>"),
                 plotOutput("p_fret", height = 300)),

        tabPanel("3 · 체액·Starling",
                 br(), plotOutput("p_compart", height = 340),
                 plotOutput("p_flux", height = 320),
                 plotOutput("p_pressures", height = 300)),

        tabPanel("4 · 제어기",
                 br(), plotOutput("p_ctrl", height = 400),
                 htmlOutput("creep_box"),
                 plotOutput("p_ctrl2", height = 320)),

        tabPanel("5 · 소생술 합병증",
                 br(), plotOutput("p_iap", height = 340),
                 htmlOutput("acs_box"),
                 plotOutput("p_wt", height = 320)),

        tabPanel("6 · 대사항진",
                 br(), plotOutput("p_ree", height = 360),
                 htmlOutput("ree_box"),
                 plotOutput("p_neuro", height = 320)),

        tabPanel("7 · 체성분",
                 br(), plotOutput("p_lbm", height = 360),
                 HTML("<p>순 균형은 <b>큰 두 값의 차</b>입니다. 어느 한쪽이 20 % 움직이면
                       부호가 바뀝니다. 그래서 프로프라놀롤(분해 ↓, 합성 ↑)과
                       옥산드롤론(합성 ↑)이 서로 다른 항을 건드리면서도 같은 방향으로
                       움직입니다.</p>"),
                 plotOutput("p_body", height = 320)),

        tabPanel("8 · 창상·폐쇄",
                 br(), plotOutput("p_wound", height = 380),
                 htmlOutput("closure_box"),
                 plotOutput("p_scar", height = 300)),

        tabPanel("9 · 염증·감염",
                 br(), plotOutput("p_cyto", height = 340),
                 plotOutput("p_infect", height = 340)),

        tabPanel("10 · 약리·PK",
                 br(), plotOutput("p_pk", height = 340),
                 htmlOutput("pk_box"),
                 plotOutput("p_arc", height = 300)),

        tabPanel("11 · 엔드포인트·검증",
                 br(), h4("보정 목표 대비 현재 시뮬레이션"),
                 tableOutput("valid_tbl"),
                 h4("개정 Baux 점수 (외부 검증자)"),
                 htmlOutput("rbaux_box"),
                 plotOutput("p_rbaux", height = 340)),

        tabPanel("12 · 시나리오 비교",
                 br(),
                 checkboxGroupInput("scens", "비교할 시나리오",
                                    choices = names(SCEN),
                                    selected = c("parkland_fixed", "parkland_titrated",
                                                 "creep_uncapped", "albumin_0h", "vitc"),
                                    inline = TRUE),
                 actionButton("gocmp", "비교 실행", class = "btn-primary"),
                 br(), br(),
                 tableOutput("cmp_tbl"),
                 plotOutput("p_cmp1", height = 340),
                 plotOutput("p_cmp2", height = 340))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  sim <- eventReactive(input$go, {
    fm <- as.numeric(input$formula)
    par <- list(
      WT = input$WT, AGE = input$AGE, TBSA = input$TBSA,
      INH = as.numeric(input$INH),
      FORMULA = if (abs(fm - 2.0001) < 1e-6) 2 else fm,
      TITRATE = as.numeric(input$titrate),
      UOTGT = input$uotgt, GCTRL = input$gctrl, RMAX = input$rmax,
      KOPI = input$kopi,
      FCOL = input$fcol,
      ALBSTART = if (input$fcol > 0) input$albstart else 1e9,
      TEXC = input$texc*24, EXCINT = input$excint,
      KTOP = input$topical, TAMB = input$tamb)

    ev <- surg_events(WT = input$WT, TBSA = input$TBSA,
                      texc = input$texc*24, excint = input$excint)
    dv <- drug_events(WT = input$WT,
                      prop_mgkgday = input$prop, prop_start = 120,
                      prop_n = ceiling(input$end*24/6),
                      oxa_mg = input$oxa, oxa_start = 120,
                      oxa_n = ceiling(input$end*24/12),
                      vitc_mgkgh = if (input$vitc) 66 else 0)
    ev <- if (is.null(ev)) dv else if (is.null(dv)) ev else rbind(ev, dv)

    m <- mod %>% param(par)
    out <- if (is.null(ev)) {
      m %>% mrgsim(end = input$end*24, delta = 1, hmax = 0.05)
    } else {
      m %>% ev(as.ev(ev)) %>% mrgsim(end = input$end*24, delta = 1, hmax = 0.05)
    }
    as_tibble(out)
  }, ignoreNULL = FALSE)

  at24 <- function(d) d[which.min(abs(d$time - 24)), ]

  ## ---------------- 1. patient ----------------
  output$patient_tbl <- renderTable({
    d <- sim(); a <- at24(d)
    tibble::tibble(
      항목 = c("체중", "나이", "%TBSA", "흡입 손상", "개정 Baux 점수",
               "rBaux 예측 사망률", "24시간 투여량", "in:out 비",
               "혈장량 최저", "24시간 COP", "최대 REE", "창상 폐쇄", "모델 사망률"),
      값 = c(sprintf("%.0f kg", d$WT[1]), sprintf("%.0f y", d$AGE[1]),
             sprintf("%.0f %%", d$TBSA[1]), ifelse(d$INH[1] > 0, "있음", "없음"),
             sprintf("%.0f", a$RBAUX), sprintf("%.1f %%", a$RBMORT),
             sprintf("%.2f mL/kg/%%  (%.0f mL/kg)", a$MLKGP, a$MLKG),
             sprintf("%.2f", a$INOUT),
             sprintf("%.0f %% of baseline", min(d$VPPCT[d$time <= 24])),
             sprintf("%.1f mmHg", a$COPt),
             sprintf("%.0f %% of predicted", max(d$REEPCT)),
             { z <- which(d$AOPEN <= 0.05*d$TBSA[1])
               if (length(z)) sprintf("%.0f 일", d$time[z[1]]/24) else "미폐쇄" },
             sprintf("%.1f %%", tail(d$MORT, 1))))
  })

  output$summary_box <- renderUI({
    d <- sim(); a <- at24(d)
    col <- if (a$INOUT > 1.6) "#c0392b" else if (a$INOUT > 1.3) "#b9770e" else "#1e8449"
    HTML(sprintf(
      "<div style='padding:12px;border-left:6px solid %s;background:#fbfbfb'>
       <b>24시간 투여량 %.2f mL/kg/%%TBSA</b> (파클랜드 처방의 <b>%.2f배</b>) ·
       혈장량 최저 <b>%.0f %%</b> · 24시간 교질삼투압 <b>%.1f mmHg</b> ·
       최대 복강내압 <b>%.1f mmHg</b> · 최대 REE <b>%.0f %%</b> ·
       60일 사망률 <b>%.1f %%</b> (rBaux 예측 %.1f %%)</div>",
      col, a$MLKGP, a$INOUT, min(d$VPPCT[d$time <= 24]), a$COPt,
      max(d$IAPt), max(d$REEPCT), tail(d$MORT, 1), a$RBMORT))
  })

  ## ---------------- 2. the convex function ----------------
  output$p_convex <- renderPlot({
    d <- sim()
    cc <- seq(0.5, 8, by = 0.05)
    traj <- d %>% filter(time <= 24)
    ggplot() +
      geom_line(aes(cc, LP(cc)), linewidth = 1.1, colour = "#2471a3") +
      geom_path(data = traj, aes(CPt, COPt), colour = "#c0392b", linewidth = 1.2) +
      geom_point(data = traj[c(1, nrow(traj)), ], aes(CPt, COPt),
                 colour = "#c0392b", size = 3) +
      geom_segment(aes(x = 7, xend = 3.5, y = LP(7), yend = LP(7)),
                   linetype = 3, colour = "grey40") +
      geom_segment(aes(x = 3.5, xend = 3.5, y = LP(7), yend = LP(3.5)),
                   linetype = 3, colour = "grey40") +
      annotate("text", x = 5.2, y = LP(7) + 1.5,
               label = "protein -50 %  ->  COP -62 %", colour = "grey30") +
      labs(x = "총 단백질 Total protein (g/dL)",
           y = "교질삼투압 COP (mmHg)",
           title = "볼록성이 만드는 비용: 매 리터가 직전 리터보다 비싸다") + THEME
  })

  output$p_fret <- renderPlot({
    d <- sim() %>% filter(time <= 48)
    ggplot(d, aes(time)) +
      geom_line(aes(y = COPt, colour = "혈장 COP (mmHg)"), linewidth = 1) +
      geom_line(aes(y = CPt*4, colour = "총 단백질 ×4 (g/dL)"), linewidth = 1) +
      labs(x = "시간 (h)", y = "", title = "교질삼투압의 붕괴 — 이것이 곧 f_ret의 붕괴") + THEME
  })

  ## ---------------- 3. fluid ----------------
  output$p_compart <- renderPlot({
    d <- sim() %>% filter(time <= 168)
    d %>% select(time, VP, VIB, VIU, VASC, EVLW) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name, VP = "혈장", VIB = "간질(화상부)",
                           VIU = "간질(비화상부)", VASC = "복강", EVLW = "폐수분")) %>%
      ggplot(aes(time, value/1000, fill = name)) +
      geom_area(alpha = 0.85) +
      labs(x = "시간 (h)", y = "부피 (L)", title = "리터는 어디로 갔는가") + THEME
  })

  output$p_flux <- renderPlot({
    d <- sim() %>% filter(time <= 72)
    ggplot(d, aes(time)) +
      geom_line(aes(y = VPPCT, colour = "혈장량 (% 기저)"), linewidth = 1) +
      geom_hline(yintercept = 100, linetype = 3) +
      labs(x = "시간 (h)", y = "%", title = "혈장량 — 투여한 것이 아니라 유지된 것") + THEME
  })

  output$p_pressures <- renderPlot({
    d <- sim() %>% filter(time <= 48)
    ggplot(d, aes(time)) +
      geom_line(aes(y = COPt, colour = "혈장 COP"), linewidth = 1) +
      geom_line(aes(y = IAPt, colour = "복강내압"), linewidth = 1) +
      geom_hline(yintercept = c(12, 20), linetype = c(3, 2), colour = "#c0392b") +
      annotate("text", x = 2, y = 21, label = "ACS", colour = "#c0392b", hjust = 0) +
      labs(x = "시간 (h)", y = "mmHg", title = "압력들") + THEME
  })

  ## ---------------- 4. controller ----------------
  output$p_ctrl <- renderPlot({
    d <- sim() %>% filter(time <= 48)
    tgt <- input$uotgt*input$WT
    ggplot(d, aes(time)) +
      geom_line(aes(y = RSTATE*100, colour = "투여속도 배수 ×100"), linewidth = 1.1) +
      geom_line(aes(y = UOWIN, colour = "측정 소변량 (mL/h)"), linewidth = 1) +
      geom_line(aes(y = COPt*8, colour = "혈장 COP ×8"), linewidth = 1) +
      geom_hline(yintercept = tgt, linetype = 2, colour = "grey30") +
      annotate("text", x = 44, y = tgt + 8, label = "설정점", colour = "grey30") +
      labs(x = "시간 (h)", y = "",
           title = "제어기는 자기가 파괴하는 이득 위에서 더 세게 민다") + THEME
  })

  output$creep_box <- renderUI({
    d <- sim(); a <- at24(d)
    msg <- if (a$INOUT > 1.6)
      "<b style='color:#c0392b'>수액 크립</b> — 제어기가 상한에 걸렸습니다."
    else if (a$INOUT > 1.3)
      "<b style='color:#b9770e'>경계 범위</b> — 문헌의 in:out 1.2–1.6 안에 있습니다."
    else "<b style='color:#1e8449'>처방량 근처</b>."
    HTML(sprintf("<div style='padding:10px;background:#f7f9fa'>%s 최종 in:out <b>%.2f</b>,
                  총 <b>%.0f mL/kg</b>. 250 mL/kg를 넘으면 복강내 고혈압 위험이 급증합니다
                  (모델 최대 복강내압 %.1f mmHg).</div>",
                 msg, a$INOUT, a$MLKG, max(d$IAPt)))
  })

  output$p_ctrl2 <- renderPlot({
    d <- sim() %>% filter(time <= 72)
    ggplot(d, aes(time)) +
      geom_line(aes(y = UOC/1000, colour = "누적 소변 (L)"), linewidth = 1) +
      geom_line(aes(y = FIN/1000, colour = "누적 투여 (L)"), linewidth = 1) +
      labs(x = "시간 (h)", y = "L", title = "들어간 것과 나온 것") + THEME
  })

  ## ---------------- 5. morbidity ----------------
  output$p_iap <- renderPlot({
    d <- sim() %>% filter(time <= 168)
    ggplot(d, aes(time, IAPt)) + geom_line(linewidth = 1.1, colour = "#c0392b") +
      geom_hline(yintercept = c(12, 20), linetype = c(3, 2)) +
      annotate("text", x = 5, y = 12.8, label = "복강내 고혈압 (IAH)", hjust = 0, size = 3.5) +
      annotate("text", x = 5, y = 20.8, label = "복부구획증후군 (ACS)", hjust = 0, size = 3.5) +
      labs(x = "시간 (h)", y = "복강내압 (mmHg)",
           title = "이 그래프의 내용은 전부 치료가 만든 것이다") + THEME
  })

  output$acs_box <- renderUI({
    d <- sim()
    HTML(sprintf("<div style='padding:10px;background:#fdf2f0'>최대 복강내압
      <b>%.1f mmHg</b> · 최대 폐수분 <b>%.0f mL</b> · 최대 체중 증가 <b>%.0f %%</b>.
      복강내압은 신관류를 깎아 소변량을 줄이고, 제어기는 그것을 저혈량으로 읽어
      <b>더 많이 투여합니다</b> — 이것이 두 번째, 양성 되먹임 고리입니다.</div>",
      max(d$IAPt), max(d$EVLW), max(d$WTGAIN)))
  })

  output$p_wt <- renderPlot({
    d <- sim() %>% filter(time <= 336)
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = WTGAIN, colour = "체액성 체중 증가 (%)"), linewidth = 1.1) +
      geom_line(aes(y = OEDEMA, colour = "세포외 부종 (L)"), linewidth = 1) +
      labs(x = "일 (days)", y = "", title = "부종의 형성과 동원") + THEME
  })

  ## ---------------- 6. hypermetabolism ----------------
  output$p_ree <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = REEPCT, colour = "REE (% 예측치)"), linewidth = 1.2) +
      geom_line(aes(y = AOPENT + 100, colour = "개방 창상 면적 + 100 (%TBSA)"),
                linewidth = 1, linetype = 2) +
      geom_hline(yintercept = 100, linetype = 3) +
      labs(x = "일 (days)", y = "",
           title = "느린 시계: 구동원은 입원 시 %TBSA가 아니라 오늘의 개방 면적") + THEME
  })

  output$ree_box <- renderUI({
    d <- sim()
    z <- which(d$AOPEN <= 0.05*d$TBSA[1])
    cl <- if (length(z)) d$time[z[1]]/24 else NA
    HTML(sprintf("<div style='padding:10px;background:#fefaf0'>최대 REE
      <b>%.0f %%</b> (%s일째) · 관찰 종료 시 <b>%.0f %%</b> · 창상 폐쇄 <b>%s</b>.
      대사항진은 창상이 닫힌 <b>뒤에도</b> 남습니다 (모델 감쇠 시간상수 약 104일) —
      이것이 이 증후군에서 가장 잘 기록된 특징입니다.</div>",
      max(d$REEPCT), round(d$time[which.max(d$REEPCT)]/24),
      tail(d$REEPCT, 1), ifelse(is.na(cl), "미폐쇄", sprintf("%.0f 일", cl))))
  })

  output$p_neuro <- renderPlot({
    d <- sim() %>% filter(time <= 720)
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = CAT, colour = "카테콜아민 (×정상)"), linewidth = 1) +
      geom_line(aes(y = COR/4, colour = "코르티솔 /4 (µg/dL)"), linewidth = 1) +
      geom_line(aes(y = TCORE - 33, colour = "심부체온 −33 (°C)"), linewidth = 1) +
      labs(x = "일 (days)", y = "", title = "전달자 층") + THEME
  })

  ## ---------------- 7. body composition ----------------
  output$p_lbm <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24, LBMPCT)) +
      geom_line(linewidth = 1.2, colour = "#7d3c98") +
      geom_hline(yintercept = 0, linetype = 3) +
      geom_hline(yintercept = -10, linetype = 2, colour = "#c0392b") +
      annotate("text", x = 2, y = -11.5, hjust = 0, size = 3.5, colour = "#c0392b",
               label = "−10 %: 창상 치유 장애가 시작되는 지점") +
      labs(x = "일 (days)", y = "제지방량 변화 (%)",
           title = "순 균형은 큰 두 값의 차 — 그래서 파라미터에 극도로 민감하다") + THEME
  })

  output$p_body <- renderPlot({
    d <- sim() %>% filter(time <= 1440)
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = LBM, colour = "제지방량 (kg)"), linewidth = 1) +
      geom_line(aes(y = FATM, colour = "지방량 (kg)"), linewidth = 1) +
      geom_line(aes(y = BMC*20, colour = "골무기질량 ×20"), linewidth = 1) +
      geom_line(aes(y = HFAT/40, colour = "간 지방 /40 (g)"), linewidth = 1) +
      labs(x = "일 (days)", y = "", title = "체성분과 간") + THEME
  })

  ## ---------------- 8. wound ----------------
  output$p_wound <- renderPlot({
    d <- sim()
    d %>% select(time, AOPEN, AGRF, AHEAL, ADON) %>%
      pivot_longer(-time) %>%
      mutate(name = recode(name, AOPEN = "개방 창상", AGRF = "이식 후 생착 대기",
                           AHEAL = "치유 완료", ADON = "공여부(개방)")) %>%
      ggplot(aes(time/24, value, fill = name)) + geom_area(alpha = 0.85) +
      labs(x = "일 (days)", y = "%TBSA",
           title = "창상 폐쇄 — 공여부가 대형 화상의 율속 자원이다") + THEME
  })

  output$closure_box <- renderUI({
    d <- sim()
    z <- which(d$AOPEN <= 0.05*d$TBSA[1])
    pool <- (100 - d$TBSA[1])*0.60*3
    HTML(sprintf("<div style='padding:10px;background:#f1f9f4'>사용 가능한 공여부
      (미화상 피부의 60 %%, 메시 1:3) = <b>%.0f %%TBSA</b> 상당 · 필요 면적
      <b>%.0f %%TBSA</b> · 폐쇄까지 <b>%s</b>. 필요 면적이 공여부를 넘으면 공여부가
      다시 나을 때까지(약 7일) 기다려야 하며, 그래서 폐쇄 시간은 %%TBSA에 대해
      심하게 비선형입니다.</div>",
      pool, d$TBSA[1],
      ifelse(length(z), sprintf("%.0f 일", d$time[z[1]]/24), "관찰 기간 내 미폐쇄")))
  })

  output$p_scar <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24, SCARC)) + geom_line(linewidth = 1.1, colour = "#c0392b") +
      labs(x = "일 (days)", y = "반흔 콜라겐 (a.u.)",
           title = "비대성 반흔 — 위험은 폐쇄까지 걸린 시간을 따라간다") + THEME
  })

  ## ---------------- 9. inflammation / infection ----------------
  output$p_cyto <- renderPlot({
    d <- sim() %>% filter(time <= 720)
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = IL6, colour = "IL-6 (pg/mL)"), linewidth = 1) +
      geom_line(aes(y = TNFA*4, colour = "TNF-α ×4"), linewidth = 1) +
      geom_line(aes(y = CRP/2, colour = "CRP /2 (mg/L)"), linewidth = 1) +
      geom_line(aes(y = HLADR*2, colour = "단핵구 HLA-DR ×2 (%)"), linewidth = 1) +
      labs(x = "일 (days)", y = "",
           title = "IL-6는 알부민 합성을 억제한다 — 알부민은 음성 급성기 단백") + THEME
  })

  output$p_infect <- renderPlot({
    d <- sim() %>% filter(time <= 1440)
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = BWD, colour = "창상 균수 (log₁₀ CFU/g)"), linewidth = 1.1) +
      geom_line(aes(y = BSYS, colour = "전신 세균 부하 (a.u.)"), linewidth = 1.1) +
      geom_hline(yintercept = 5, linetype = 2, colour = "#c0392b") +
      annotate("text", x = 1, y = 5.3, hjust = 0, size = 3.5, colour = "#c0392b",
               label = "10⁵ CFU/g — 침습성 감염의 고전적 정량 역치") +
      labs(x = "일 (days)", y = "", title = "밀도는 밀도, 침습은 면적") + THEME
  })

  ## ---------------- 10. pharmacology ----------------
  output$p_pk <- renderPlot({
    d <- sim() %>% filter(time <= 720)
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = CPRt, colour = "프로프라놀롤 (ng/mL)"), linewidth = 1) +
      geom_line(aes(y = CVANt, colour = "반코마이신 (mg/L)"), linewidth = 1) +
      labs(x = "일 (days)", y = "농도", title = "약물 농도") + THEME
  })

  output$pk_box <- renderUI({
    d <- sim()
    i <- which.min(abs(d$time - 288))
    cp <- d$CPt[i]
    HTML(sprintf("<div style='padding:10px;background:#eefaf8'>12일째 총 단백질
      <b>%.2f g/dL</b>. 90 %% 결합 약물의 유리분율은 0.10 → <b>%.3f</b>로 오릅니다.
      저추출 약물이라면 <code>CL = fu × CL_int</code> 이므로 정상상태의 <b>유리</b>
      농도는 변하지 않고 <b>총</b> 농도만 떨어집니다 — 총 농도 기준 TDM은
      노출 부족을 <b>과잉 진단</b>합니다. 반대로 증가된 신클리어런스(×%.2f)는
      유리 농도 자체를 낮추므로 TDM이 <b>과소 진단</b>합니다. 화상 환자의 두
      대표적 PK 변화는 정작 중요한 양에 <b>반대 부호</b>로 작용합니다.</div>",
      cp, min(0.10*(7/max(cp, 0.5)), 1), 1 + 0.55))
  })

  output$p_arc <- renderPlot({
    d <- sim() %>% filter(time <= 720)
    ggplot(d, aes(time/24)) +
      geom_line(aes(y = SCR, colour = "혈청 크레아티닌 (mg/dL)"), linewidth = 1) +
      geom_line(aes(y = LAC, colour = "젖산 (mmol/L)"), linewidth = 1) +
      geom_line(aes(y = GLC/50, colour = "혈당 /50 (mg/dL)"), linewidth = 1) +
      labs(x = "일 (days)", y = "", title = "신기능·관류·혈당") + THEME
  })

  ## ---------------- 11. validation ----------------
  output$valid_tbl <- renderTable({
    d <- sim(); a <- at24(d)
    vals <- c(a$MLKGP, a$INOUT, min(d$VPPCT[d$time <= 24]), a$COPt,
              max(d$WTGAIN), max(d$REEPCT),
              d$LBMPCT[which.min(abs(d$time - 336))], a$MLKG)
    TARGETS %>% mutate(
      모델값 = sprintf("%.2f", vals),
      판정 = ifelse(vals >= lo & vals <= hi, "○ 범위 내", "△ 범위 밖")) %>%
      select(지표 = metric, 문헌값 = lit, 모델값, 판정)
  })

  output$rbaux_box <- renderUI({
    d <- sim(); a <- at24(d)
    HTML(sprintf("<div style='padding:10px;background:#f7f9fa'>개정 Baux 점수
      = 나이 + %%TBSA + 17×흡입 = <b>%.0f</b> → 로지스틱 예측 사망률
      <b>%.1f %%</b>. 모델의 기전적 사망률은 <b>%.1f %%</b>입니다.
      <b>rBaux는 모델의 입력이 아닙니다</b> — 위험은 개방 창상의 면적과
      개방 기간, 나이, 흡입 손상, 패혈증, 구획증후군, ARDS, 제지방량 손실에서
      계산되며 rBaux는 그 결과를 채점하는 데만 쓰입니다.</div>",
      a$RBAUX, a$RBMORT, tail(d$MORT, 1)))
  })

  output$p_rbaux <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/24, MORT)) + geom_line(linewidth = 1.2, colour = "#34495e") +
      geom_hline(aes(yintercept = RBMORT), linetype = 2, colour = "#c0392b") +
      annotate("text", x = 1, y = d$RBMORT[1] + 2, hjust = 0, size = 3.5,
               colour = "#c0392b", label = "rBaux 예측치") +
      labs(x = "일 (days)", y = "누적 사망률 (%)",
           title = "기전적으로 계산된 사망률 vs rBaux") + THEME
  })

  ## ---------------- 12. scenario comparison ----------------
  cmp <- eventReactive(input$gocmp, {
    req(length(input$scens) > 0)
    bind_rows(lapply(input$scens, function(s)
      run_scenario(s, WT = input$WT, AGE = input$AGE, TBSA = input$TBSA,
                   INH = as.numeric(input$INH), end = input$end*24)))
  })

  output$cmp_tbl <- renderTable({
    cmp() %>% group_by(scenario) %>% group_modify(~summarise_run(.x)) %>%
      ungroup() %>% select(-scenario1) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
  })

  output$p_cmp1 <- renderPlot({
    cmp() %>% filter(time <= 48) %>%
      ggplot(aes(time, MLKGP, colour = scenario)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = c(4, 5.2, 6.7), linetype = c(2, 3, 3)) +
      labs(x = "시간 (h)", y = "누적 mL/kg/%TBSA",
           title = "24시간 수액량 — 점선은 파클랜드 처방(4)과 관찰 범위(5.2–6.7)") + THEME
  })

  output$p_cmp2 <- renderPlot({
    cmp() %>% ggplot(aes(time/24, REEPCT, colour = scenario)) +
      geom_line(linewidth = 1) +
      labs(x = "일 (days)", y = "REE (% 예측치)",
           title = "느린 시계") + THEME
  })
}

shinyApp(ui, server)
