# =====================================================================
# Short Bowel Syndrome with Chronic Intestinal Failure (SBS-IF)
# 단장증후군 · 만성 장부전 — Shiny 대시보드 (10 tabs)
# ---------------------------------------------------------------------
# The dashboard is organised around the model's central claim: PARENTERAL
# SUPPORT VOLUME IS AN OUTPUT, NOT AN INPUT. So the first thing the user
# does is set an ANATOMY, and the app then SOLVES for the prescription
# that closes the fluid and energy balances. Nothing here lets you type in
# a PN volume, because in this model you cannot.
#
#   Tab 1   환자 프로파일        anatomy, diet, drinking behaviour -> phenotype
#   Tab 2   두 보존 방정식        the fluid and energy budgets, term by term
#   Tab 3   물의 역설            per-stream jejunal driving force and ORS sweep
#   Tab 4   PK                  GLP-2 analogue concentration and occupancy
#   Tab 5   점막 적응            villus, transporters, contact time, citrulline
#   Tab 6   임상 엔드포인트       PN volume, infusion nights, 20% responder
#   Tab 7   시나리오 비교         drugs and adjuncts head to head
#   Tab 8   바이오마커            citrulline, Mg, bicarbonate, micronutrients
#   Tab 9   장기 합병증           IFALD, catheter/access, kidney, bone over 2 y
#   Tab 10  위닝 프로토콜         the weaning algorithm, and the placebo artefact
#
# RUN:  Rscript -e 'shiny::runApp("sbs_shiny_app.R", port = 8080)'
#       (the mrgsolve model is compiled once at start-up; allow ~30 s)
#
# NOTE: educational / research tool. Not validated for clinical
# decision-making, prescribing or regulatory submission.
# =====================================================================

library(shiny)

# ---------------------------------------------------------------------
# Load the model without triggering its scenario/diagnostic suite
# ---------------------------------------------------------------------
Sys.setenv(SBS_NORUN = "1")
MODEL_FILE <- "sbs_mrgsolve_model.R"
if (!file.exists(MODEL_FILE)) {
  # when launched from another working directory, look next to this file
  here <- tryCatch(dirname(normalizePath(sys.frames()[[1]]$ofile)),
                   error = function(e) ".")
  MODEL_FILE <- file.path(here, "sbs_mrgsolve_model.R")
}
suppressMessages(source(MODEL_FILE))   # also defines %||%, mp(), burnin(), sim()

PAL <- c(base = "#3f6d99", drug = "#a01050", alt = "#4a7a4a",
         warn = "#a04040", gold = "#a08000", grey = "#777777",
         purple = "#6a3a8a", orange = "#a04a10")

# ---------------------------------------------------------------------
# Small plotting helpers (base graphics only, no external dependencies)
# ---------------------------------------------------------------------
tsplot <- function(dl, col, ylab, main, cols = NULL, lty = NULL, hline = NULL,
                   legpos = "topright", xlab = "일 (days)") {
  dl <- dl[!vapply(dl, is.null, logical(1))]
  if (!length(dl)) { plot.new(); return(invisible()) }
  rng <- range(unlist(lapply(dl, function(d) d[[col]])), na.rm = TRUE)
  if (!is.finite(rng[1])) rng <- c(0, 1)
  if (diff(rng) < 1e-9) rng <- rng + c(-1, 1) * max(1e-3, abs(rng[1]) * 0.1)
  if (is.null(cols)) cols <- unname(PAL)[seq_along(dl)]
  if (is.null(lty))  lty  <- rep(1, length(dl))
  xr <- range(unlist(lapply(dl, function(d) d$time)))
  par(mar = c(4.2, 4.6, 3.0, 1.2))
  plot(NA, xlim = xr, ylim = rng, xlab = xlab, ylab = ylab, main = main,
       las = 1, cex.main = 1.02, font.main = 1)
  grid(col = "#e6e6e6", lty = 1)
  if (!is.null(hline)) abline(h = hline, col = "#bbbbbb", lty = 2)
  for (i in seq_along(dl))
    lines(dl[[i]]$time, dl[[i]][[col]], col = cols[i], lwd = 2.1, lty = lty[i])
  if (length(dl) > 1 || !is.null(names(dl)))
    legend(legpos, legend = names(dl), col = cols, lty = lty, lwd = 2.1,
           bty = "n", cex = 0.86)
}

barplot2 <- function(v, main, ylab, cols = NULL, horiz = FALSE, digits = 2) {
  par(mar = if (horiz) c(4.2, 11, 3, 1.2) else c(7.5, 4.6, 3, 1.2))
  if (is.null(cols)) cols <- rep(PAL["base"], length(v))
  b <- barplot(v, main = main, ylab = if (horiz) "" else ylab, col = cols,
               border = NA, las = if (horiz) 1 else 2, horiz = horiz,
               cex.names = 0.82, cex.main = 1.02, font.main = 1)
  if (!horiz) abline(h = 0, col = "#888888")
  text(if (horiz) v else b, if (horiz) b else v,
       labels = formatC(v, format = "f", digits = digits),
       pos = if (horiz) 4 else ifelse(v >= 0, 3, 1), cex = 0.76, xpd = TRUE)
  invisible(b)
}

waterfall <- function(labels, values, main, ylab) {
  par(mar = c(8.5, 4.8, 3, 1.2))
  cum <- c(0, cumsum(values))
  ylim <- range(c(cum, 0)) * 1.12
  plot(NA, xlim = c(0.4, length(values) + 0.6), ylim = ylim, xaxt = "n",
       xlab = "", ylab = ylab, main = main, las = 1, cex.main = 1.02,
       font.main = 1)
  grid(col = "#ececec", lty = 1); abline(h = 0, col = "#666666")
  for (i in seq_along(values)) {
    rect(i - 0.36, cum[i], i + 0.36, cum[i + 1], border = NA,
         col = if (values[i] >= 0) "#5b9bd5" else "#d9776a")
    text(i, max(cum[i], cum[i + 1]),
         sprintf("%+.2f", values[i]), pos = 3, cex = 0.72, xpd = TRUE)
  }
  lines(seq_along(values) + 0.36, cum[-1], type = "s", col = "#888888", lty = 3)
  axis(1, at = seq_along(values), labels = labels, las = 2, cex.axis = 0.78)
  invisible(cum)
}

kv_table <- function(m) {
  renderTable({ m }, rownames = FALSE, striped = TRUE, spacing = "xs",
              digits = 3)
}

# ---------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif; }
    h4 { margin-top: 4px; }
    .thesis { background:#fffbe8; border-left:5px solid #a08000;
              padding:9px 13px; margin-bottom:11px; font-size:13.5px; }
    .caution { background:#fdf0f0; border-left:5px solid #a04040;
               padding:8px 12px; margin:9px 0; font-size:12.8px; }
    .note { color:#555; font-size:12.6px; margin:6px 0 12px 0; }
    .well { background:#f7f9fb; }
  "))),

  titlePanel("단장증후군 · 만성 장부전 (SBS-IF) QSP 대시보드 — Short Bowel Syndrome with Chronic Intestinal Failure"),

  div(class = "thesis",
      strong("이 모델의 전제: 비경구영양 용적은 입력이 아니라 출력이다."),
      " 왼쪽에서 설정하는 것은 ",
      em("해부·식이·음수 행동"), "이며, 앱은 수분·나트륨·에너지 보존 방정식을 닫는 ",
      "처방을 ", strong("계산"), "합니다. PN 용적을 직접 입력하는 칸은 없습니다 — ",
      "이 모델에서는 그것이 가능하지 않기 때문입니다.",
      br(),
      tags$span(style = "color:#555",
        "The thing you set is the anatomy; the prescription is SOLVED from the ",
        "fluid and energy budgets. There is deliberately no input box for PN volume.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① 잔여 해부 (Anatomy = parameter)"),
      selectInput("anat", "Messing 해부 유형",
        c("Type 2 공장-대장 (jejunocolic)"      = "jejunocolic",
          "Type 1 말단 공장루 (end-jejunostomy)" = "end_jejunostomy",
          "Type 3 공장-회장-대장 (jejuno-ileo-colic)" = "jejunoileocolic",
          "사용자 지정 (custom)"                  = "custom"),
        selected = "jejunocolic"),
      sliderInput("sbl", "잔여 소장 길이 SBL (cm)", 20, 400, 80, 5),
      sliderInput("ileum", "잔여 회장 (cm) — B12·담즙산·L세포", 0, 200, 0, 5),
      sliderInput("colon", "연속성 대장 비율 COLONFRAC", 0, 1, 0.5, 0.05),
      checkboxInput("icv", "회맹판 보존 (ileocaecal valve)", FALSE),
      sliderInput("mucqual", "잔여 점막의 질 (1 = 정상)", 0.5, 1, 1, 0.05),

      h4("② 식이와 음수 (diet & drinking)"),
      sliderInput("ors", "ORS 비율 (between-meal fluid)", 0, 0.95, 0.40, 0.05),
      sliderInput("drink", "식사 사이 음수량 (L/d)", 0.5, 5, 2.5, 0.1),
      sliderInput("kcal", "경구 열량 기준 (kcal/d)", 1200, 4000, 2400, 100),
      sliderInput("entfrac", "경장 섭취 배율 (0 = 완전 금식)", 0, 1.2, 1.0, 0.02),
      sliderInput("dietna", "식이 나트륨 (mmol/d)", 60, 300, 160, 10),

      h4("③ 약물 (drug)"),
      selectInput("drug", "GLP-2 요법",
        c("없음 (none)" = "none",
          "테두글루타이드 0.05 mg/kg/일 SC" = "ted",
          "아프라글루타이드 5 mg 주 1회"     = "apra",
          "글레파글루타이드 10 mg 주 1회"    = "glep",
          "글레파글루타이드 10 mg 주 2회"    = "glep2",
          "천연 GLP-2 등몰 1일 1회"          = "nat1",
          "천연 GLP-2 등몰 1일 2회"          = "nat2"),
        selected = "ted"),
      checkboxGroupInput("adj", "보조요법 (adjuncts)",
        c("로페라마이드 16 mg/일" = "lop",
          "고용량 PPI"            = "ppi",
          "옥트레오타이드 100 µg tid" = "oct",
          "소마트로핀 0.1 mg/kg/일 x 4주" = "gh",
          "콜레스티라민 4 g tid" = "chol",
          "리팍시민 550 mg bid x 14일" = "rifax",
          "구연산칼륨·저옥살산 식이" = "cit",
          "타우롤리딘 카테터 잠금" = "lock")),
      selectInput("ile", "PN 지질 유탁액 (lipid emulsion)",
        c("대두유 (soybean, phytosterol ~350)" = "soy",
          "SMOF 복합 (composite)"              = "smof",
          "어유 기반 (fish-oil, phytosterol 0)" = "fish",
          "지질 최소화 (lipid minimisation)"    = "min"),
        selected = "soy"),

      h4("④ 위닝 프로토콜 (weaning protocol)"),
      checkboxInput("weanon", "요량 기반 위닝 알고리즘 작동", TRUE),
      sliderInput("urtrig", "PN 감량 트리거 (요량 / 기저 요량)", 1.02, 1.40, 1.10, 0.01),
      sliderInput("ktaper", "최대 감량 속도 (1/일)", 0.001, 0.020, 0.006, 0.001),
      sliderInput("days", "시뮬레이션 기간 (일)", 60, 1095, 168, 7),
      actionButton("go", "시뮬레이션 실행 (run)", class = "btn-primary",
                   style = "width:100%; margin-top:8px;"),
      div(class = "note",
          "해부를 바꾸면 모델은 먼저 정상상태(burn-in)를 다시 찾은 뒤 ",
          "시나리오를 실행합니다. 그래서 모든 비교가 자기 일관된 환자에서 출발합니다.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        # ---------------- TAB 1 -------------------------------------
        tabPanel("1 · 환자 프로파일",
          h4("설정된 해부가 만들어 내는 표현형 (the phenotype the anatomy produces)"),
          div(class = "note",
            "아래 값은 어느 것도 입력이 아닙니다. 잔여 장 길이·대장·회맹판·식이·음수 ",
            "행동만 주어졌고, 나머지는 모두 모델이 정상상태에서 계산한 결과입니다."),
          fluidRow(column(6, tableOutput("t1_pheno")),
                   column(6, tableOutput("t1_absorp"))),
          hr(),
          fluidRow(column(6, plotOutput("t1_anat", height = 300)),
                   column(6, plotOutput("t1_thresh", height = 300))),
          div(class = "caution",
            strong("해부는 상태가 아니라 파라미터입니다."),
            " Messing의 영구 장부전 역치(말단 공장루 ~115 cm, 공장-대장 ~60 cm, ",
            "공장-회장-대장 ~35 cm)가 세 배 가까이 다른 이유는 대장이 하루 3-4 L의 ",
            "수분과 최대 ~1000 kcal을 회수하기 때문입니다.")
        ),

        # ---------------- TAB 2 -------------------------------------
        tabPanel("2 · 두 보존 방정식",
          h4("수분·나트륨·에너지 예산 (the budgets whose residual the central line fills)"),
          fluidRow(column(6, plotOutput("t2_fluid", height = 340)),
                   column(6, plotOutput("t2_energy", height = 340))),
          hr(),
          fluidRow(column(7, tableOutput("t2_terms")),
                   column(5, plotOutput("t2_na", height = 300))),
          div(class = "note",
            "왼쪽 폭포 그림의 마지막 막대가 0에서 벗어난 만큼이 곧 정맥으로 채워야 ",
            "하는 결손입니다. 위닝 알고리즘이 하는 일은 이 결손을 다시 계산하는 것이 ",
            "아니라, 요량이라는 관측 가능한 신호를 통해 그것을 추적하는 것입니다.")
        ),

        # ---------------- TAB 3 -------------------------------------
        tabPanel("3 · 물의 역설",
          h4("공장 Na⁺ 구동력의 부호 변화 (where the sign changes)"),
          div(class = "note",
            "식사와 함께 들어온 유체(stream M)와 식사 사이에 마신 유체(stream D)는 ",
            "나트륨 농도와 포도당 함량이 전혀 다릅니다. 두 흐름에 각각 구동력을 ",
            "계산한 뒤에야 부피 비율로 평균합니다. 먼저 평균해 버리면 이 역설은 사라집니다."),
          fluidRow(column(6, plotOutput("t3_drive", height = 340)),
                   column(6, plotOutput("t3_output", height = 340))),
          hr(),
          fluidRow(column(7, tableOutput("t3_tab")),
                   column(5, plotOutput("t3_pn", height = 320))),
          div(class = "caution",
            strong("갈증이 나는 말단 공장루 환자가 맹물 2 L를 마시면 나트륨과 용적을 ",
            "잃습니다."), " 같은 2 L를 Na⁺ 90 mmol/L의 포도당-식염수로 주면 둘 다 ",
            "얻습니다. 이것은 경고문이 아니라 flux 방정식의 부호 변화입니다.")
        ),

        # ---------------- TAB 4 -------------------------------------
        tabPanel("4 · PK",
          h4("GLP-2 유사체 약동학과 수용체 점유율"),
          fluidRow(column(6, plotOutput("t4_conc", height = 330)),
                   column(6, plotOutput("t4_occ", height = 330))),
          hr(),
          fluidRow(column(6, plotOutput("t4_zoom", height = 320)),
                   column(6, tableOutput("t4_tab"))),
          div(class = "caution",
            strong("점막 반응의 시간 상수(주 단위)가 약물의 반감기(시간 단위)보다 ",
            "훨씬 길다는 것이 핵심입니다."), " 그래서 t½ 2시간의 약물을 하루 한 번 ",
            "투여해도 작동합니다. 반대로 천연 GLP-2는 DPP-4가 Ala2-Gly3에서 절단해 ",
            "t½가 약 7분이므로 등몰 1일 1회로는 점유율을 유지할 수 없습니다.")
        ),

        # ---------------- TAB 5 -------------------------------------
        tabPanel("5 · 점막 적응",
          h4("구조적·기능적 적응 — 영양소로 게이팅되는 곱셈 구조"),
          fluidRow(column(6, plotOutput("t5_vill", height = 330)),
                   column(6, plotOutput("t5_trans", height = 330))),
          hr(),
          fluidRow(column(6, plotOutput("t5_gate", height = 330)),
                   column(6, plotOutput("t5_contact", height = 330))),
          div(class = "note",
            "TROPHIC = (관내 영양소) x (1 + GLP2R 점유 항) x (1 + IGF-1 항). 곱이므로 ",
            "경장 섭취가 0이면 영양 신호도 0이 되고 영양성(trophic) 팔은 완전히 ",
            "멈춥니다. 다만 GLP-2의 운동성(ENS) 팔은 관내 영양소를 필요로 하지 ",
            "않으므로 약효가 완전히 사라지지는 않습니다 — 모델이 스스로 반증한 ",
            "설계 의도이며, 검증 가능한 구조적 예측입니다.")
        ),

        # ---------------- TAB 6 -------------------------------------
        tabPanel("6 · 임상 엔드포인트",
          h4("PN 용적 · 주당 주입 야간 수 · ≥20% 반응자"),
          fluidRow(column(6, plotOutput("t6_pn", height = 330)),
                   column(6, plotOutput("t6_nights", height = 330))),
          hr(),
          fluidRow(column(5, tableOutput("t6_tab")),
                   column(7, plotOutput("t6_out", height = 320))),
          div(class = "note",
            "STEPS 1차 종료점은 20주와 24주에 모두 PN 용적이 20% 이상 감소한 ",
            "환자입니다. 이 앱에서 그 값은 회귀식이 아니라 시험 프로토콜을 ",
            "시뮬레이션한 결과로 나옵니다. 리터가 아니라 ",
            strong("야간(night)"), "이 줄어야 카테터-일이 줄어든다는 점도 함께 보십시오.")
        ),

        # ---------------- TAB 7 -------------------------------------
        tabPanel("7 · 시나리오 비교",
          h4("약물·보조요법 정면 비교 (matched anatomy, matched protocol)"),
          checkboxGroupInput("cmp", "비교할 시나리오",
            c("무투약 (no drug)" = "none",
              "테두글루타이드" = "ted",
              "아프라글루타이드 주 1회" = "apra",
              "글레파글루타이드 주 1회" = "glep",
              "글레파글루타이드 주 2회" = "glep2",
              "천연 GLP-2 등몰 1일 1회" = "nat1",
              "ORS 최적화 (0.9)" = "ors",
              "로페라마이드" = "lop",
              "고용량 PPI" = "ppi",
              "옥트레오타이드" = "oct",
              "테두글루타이드 + ORS + 로페라마이드" = "combo"),
            selected = c("none", "ted", "ors", "combo"), inline = TRUE),
          fluidRow(column(7, plotOutput("t7_pn", height = 350)),
                   column(5, plotOutput("t7_bar", height = 350))),
          hr(),
          tableOutput("t7_tab"),
          div(class = "caution",
            strong("옥트레오타이드는 지금의 배출량을 사고 나중의 적응을 팝니다."),
            " 분비 용적을 줄이는 대신 IGF-1을 억제하고 췌장 외분비를 억제해 ",
            "지방 흡수를 떨어뜨립니다. 시간에 따라 이득의 부호가 달라지는지 ",
            "직접 확인해 보십시오.")
        ),

        # ---------------- TAB 8 -------------------------------------
        tabPanel("8 · 바이오마커",
          h4("시트룰린 · 마그네슘 · 산-염기 · 미량영양소"),
          fluidRow(column(6, plotOutput("t8_citr", height = 320)),
                   column(6, plotOutput("t8_mg", height = 320))),
          hr(),
          fluidRow(column(6, plotOutput("t8_acid", height = 320)),
                   column(6, plotOutput("t8_micro", height = 320))),
          div(class = "caution",
            strong("Mg 결핍은 PTH 분비 자체를 손상시키고 동시에 표적장기 저항을 ",
            "일으킵니다."), " 그래서 마그네슘을 채우기 전까지 저칼슘혈증은 칼슘 ",
            "투여에 반응하지 않습니다. 이것은 모델에서 주장된 것이 아니라 두 개의 ",
            "Mg 의존 항에서 생성됩니다. 시트룰린 <20 µmol/L는 장세포 질량의 직접 ",
            "지표로 취급됩니다.")
        ),

        # ---------------- TAB 9 -------------------------------------
        tabPanel("9 · 장기 합병증",
          h4("2년 경과: IFALD · 카테터와 정맥 접근로 · 신장 · 골"),
          div(class = "note",
            "현대 SBS-IF에서 사망을 결정하는 것은 남은 장이 아니라 치료의 배관입니다. ",
            "장은 좋아지는데 간과 정맥이 나빠지는 경쟁 위험 구조가 실제 의사결정 문제입니다."),
          fluidRow(column(6, plotOutput("t9_liver", height = 320)),
                   column(6, plotOutput("t9_cvc", height = 320))),
          hr(),
          fluidRow(column(6, plotOutput("t9_renal", height = 320)),
                   column(6, plotOutput("t9_bone", height = 320))),
          tableOutput("t9_tab")
        ),

        # ---------------- TAB 10 ------------------------------------
        tabPanel("10 · 위닝 프로토콜",
          h4("요량 기반 감량 알고리즘, 그리고 위약 반응이라는 인공물"),
          div(class = "note",
            "STEPS에서 PN 감량을 결정한 것은 연구자가 아니라 프로토콜이었습니다: ",
            "24시간 요량이 기저치보다 10% 이상 오르면 PN을 줄이고, 체중이나 전해질이 ",
            "떨어지면 보류합니다. 그 알고리즘을 실제 폐쇄 루프 제어기로 구현하면 ",
            "두 가지가 따라 나옵니다 — 감량은 스스로 멈추고(자기 제한), 흡수 예비력이 ",
            "남아 있는 환자는 약과 무관하게 프로토콜에 의해 위닝됩니다."),
          fluidRow(column(6, plotOutput("t10_loop", height = 330)),
                   column(6, plotOutput("t10_placebo", height = 330))),
          hr(),
          fluidRow(column(6, plotOutput("t10_trig", height = 330)),
                   column(6, tableOutput("t10_tab"))),
          div(class = "caution",
            strong("SBS 시험의 위약군 반응 중 상당 부분은 위약 효과가 아니라 ",
            "잘 짜인 진료의 효과입니다."), " 시험 등록 시 양쪽 팔 모두 경구 재수화와 ",
            "식이가 표준화되고, 요량 트리거는 그 개선을 감량으로 번역합니다.")
        )
      ),
      hr(),
      div(class = "note",
        strong("면책:"), " 교육·연구용 QSP 모델입니다. 임상 의사결정·처방·규제 제출에 ",
        "사용해서는 안 됩니다. 흡수·유체 계수는 발표된 표현형을 재현하도록 보정된 ",
        "집단 파라미터이며 개별 환자에게 적용할 수 없습니다. ",
        "근거는 ", code("sbs_references.md"), " (236편, PMID 전수 검증), 구조는 ",
        code("sbs_qsp_model.svg"), " (23 클러스터 · 276 노드) 참조.")
    )
  )
)

# ---------------------------------------------------------------------
# SERVER
# ---------------------------------------------------------------------
server <- function(input, output, session) {

  observeEvent(input$anat, {
    a <- input$anat
    if (a == "end_jejunostomy") {
      updateSliderInput(session, "sbl", value = 90)
      updateSliderInput(session, "ileum", value = 0)
      updateSliderInput(session, "colon", value = 0)
      updateCheckboxInput(session, "icv", value = FALSE)
    } else if (a == "jejunocolic") {
      updateSliderInput(session, "sbl", value = 80)
      updateSliderInput(session, "ileum", value = 0)
      updateSliderInput(session, "colon", value = 0.5)
      updateCheckboxInput(session, "icv", value = FALSE)
    } else if (a == "jejunoileocolic") {
      updateSliderInput(session, "sbl", value = 70)
      updateSliderInput(session, "ileum", value = 25)
      updateSliderInput(session, "colon", value = 1)
      updateCheckboxInput(session, "icv", value = TRUE)
    }
  }, ignoreInit = TRUE)

  # ---- parameter set from the UI -------------------------------------
  basepars <- reactive({
    ilepar <- switch(input$ile,
      soy  = list(PHYTOCONC = 350, OMEGA3 = 0, FLIPIDPN = 1,   ILEDOSE = 0.90),
      smof = list(PHYTOCONC = 100, OMEGA3 = 1, FLIPIDPN = 1,   ILEDOSE = 0.90),
      fish = list(PHYTOCONC = 0,   OMEGA3 = 1, FLIPIDPN = 1,   ILEDOSE = 0.90),
      min  = list(PHYTOCONC = 350, OMEGA3 = 0, FLIPIDPN = 1,   ILEDOSE = 0.14))
    adj <- input$adj %||% character(0)
    c(list(SBL = input$sbl, ILEUMLEN = min(input$ileum, input$sbl),
           COLONFRAC = input$colon, ICV = as.numeric(isTRUE(input$icv)),
           MUCQUAL = input$mucqual, ORSFRAC = input$ors,
           DRINKVOL = input$drink, KCALORAL = input$kcal,
           ENTFRAC = input$entfrac, DIETNA = input$dietna,
           PPIDOSE = if ("ppi" %in% adj) 2 else 0,
           OXCITRATE = as.numeric("cit" %in% adj),
           LOCKFAC = if ("lock" %in% adj) 0.25 else 1,
           WEANON = as.numeric(isTRUE(input$weanon)),
           URTRIG = input$urtrig, KTAPER = input$ktaper),
      ilepar)
  })

  # burn-in depends only on the physiology, not on the drug
  bstate <- reactive({
    p <- basepars()
    p$WEANON <- 1   # the run-in must be allowed to find the prescription
    withProgress(message = "정상상태 탐색 (burn-in)...", value = 0.4,
                 burnin(p))
  })

  runpars <- reactive({
    st <- bstate(); p <- basepars()
    modifyList(p, list(URBASE = st[["URINESM"]],
                       WTTARGET = st[["LEAN"]] + st[["FATM"]]))
  })

  wt0 <- reactive({ st <- bstate(); st[["LEAN"]] + st[["FATM"]] })

  drug_ev <- function(which, days, wt) {
    switch(which,
      none  = NULL,
      ted   = ev_ted(wt, days),
      apra  = ev_apra(days),
      glep  = ev_glep(days, ii = 7),
      glep2 = ev_glep(days, ii = 3.5),
      nat1  = ev_native(wt, days, 1),
      nat2  = ev_native(wt, days, 2),
      NULL)
  }
  adj_ev <- function(adj, days, wt) {
    e <- NULL
    add <- function(x) if (is.null(e)) x else c(e, x)
    if ("lop"   %in% adj) e <- add(ev_lop(days))
    if ("oct"   %in% adj) e <- add(ev_oct(days))
    if ("gh"    %in% adj) e <- add(ev_gh(wt, min(28, days)))
    if ("chol"  %in% adj) e <- add(ev_chol(days))
    if ("rifax" %in% adj) e <- add(ev_rifax(min(14, days)))
    e
  }

  simrun <- eventReactive(input$go, {
    st <- bstate(); p <- runpars(); d <- input$days; w <- wt0()
    adj <- input$adj %||% character(0)
    ev1 <- drug_ev(input$drug, d, w); ev2 <- adj_ev(adj, d, w)
    ee <- if (is.null(ev1)) ev2 else if (is.null(ev2)) ev1 else c(ev1, ev2)
    withProgress(message = "시뮬레이션 실행...", value = 0.6, {
      list(drug = sim(st, p, ee, end = d),
           ref  = sim(st, p, NULL, end = d),
           st = st, p = p, wt = w, days = d)
    })
  }, ignoreNULL = FALSE)

  # =================== TAB 1 ==========================================
  output$t1_pheno <- kv_table({
    r <- simrun(); d <- r$ref
    data.frame(
      `지표 (metric)` = c("처방 PN 용적 (L/주)", "주당 주입 야간 수",
        "처방 PN 열량 (kcal/일)", "변·장루 배출량 (L/일)",
        "장루 유출액 [Na⁺] (mmol/L)", "24시간 요량 (L/일)",
        "요 나트륨 (mmol/일)", "체중 (kg)", "혈장 시트룰린 (µmol/L)",
        "과식 배율 (hyperphagia)", "점막 접촉 시간 (상대값)"),
      `값 (value)` = round(c(d$PNVOLWK[1], d$PNNIGHTS[1], d$PNKCALD[1],
        d$OUT_LD[1], d$STOMANA[1], d$URINE_LD[1], d$UNA[1], d$WT[1],
        d$CITR[1], d$HYPERPH[1], d$CONTC[1]), 2),
      check.names = FALSE)
  })

  output$t1_absorp <- kv_table({
    r <- simrun(); d <- r$ref
    data.frame(
      `흡수 (absorption)` = c("탄수화물 흡수율 (%)", "단백질 흡수율 (%)",
        "지방 흡수율 (%)", "흡수 열량 (kcal/일)", "대장 SCFA 회수 (kcal/일)",
        "총 에너지 소비 (kcal/일)", "미흡수 지방 (g/일)",
        "공장 순 수분 flux (L/일)", "대장 회수 수분 (L/일)",
        "담즙산 풀 (g)", "미셀 기능 (상대값)"),
      `값 (value)` = round(c(100*d$FCHOC[1], 100*d$FPROC[1], 100*d$FFATC[1],
        d$ABSKC[1], d$SCFAKC[1], d$TEEC[1], d$FATMALC[1], d$VJEJC[1],
        d$VCOLC[1], d$BAPOOLC[1], d$MICELC[1]), 2),
      check.names = FALSE)
  })

  output$t1_anat <- renderPlot({
    r <- simrun(); p <- r$p
    v <- c(`잔여 소장\nSBL (cm)` = p$SBL,
           `그중 회장\nileum (cm)` = p$ILEUMLEN,
           `대장 %\ncolon` = 100*p$COLONFRAC,
           `회맹판\nICV` = 100*p$ICV)
    barplot2(v, "설정된 해부 (the anatomy you set)", "cm 또는 %",
             cols = c(PAL["base"], PAL["alt"], PAL["purple"], PAL["gold"]),
             digits = 0)
  })

  output$t1_thresh <- renderPlot({
    st <- bstate(); p <- basepars()
    ls <- seq(30, 250, by = 20)
    pn <- vapply(ls, function(L) {
      q <- modifyList(p, list(SBL = L, ILEUMLEN = min(p$ILEUMLEN, L)))
      s <- burnin(q, days = 500); s[["PNVOL"]] * 7
    }, numeric(1))
    par(mar = c(4.4, 4.8, 3, 1.2))
    plot(ls, pn, type = "b", pch = 19, col = PAL["base"], lwd = 2.1, las = 1,
         xlab = "잔여 소장 길이 (cm)", ylab = "정상상태 PN 용적 (L/주)",
         main = "길이 대 PN 의존도 (현재 대장/회맹판 조건에서)",
         cex.main = 1.0, font.main = 1)
    grid(col = "#ececec", lty = 1)
    abline(v = p$SBL, col = PAL["drug"], lty = 2)
    abline(h = 0.35 * 7, col = "#999999", lty = 3)
    legend("topright", c("현재 환자", "장 자율성 근사 기준"),
           col = c(PAL["drug"], "#999999"), lty = c(2, 3), bty = "n", cex = 0.85)
  })

  # =================== TAB 2 ==========================================
  output$t2_fluid <- renderPlot({
    r <- simrun(); d <- r$ref; i <- 1
    lab <- c("경구 섭취\noral", "PN", "배출\noutput", "불감손실\ninsensible", "요량\nurine")
    val <- c(d$VORALC[i], d$PNVOLWK[i]/7, -d$OUT_LD[i], -0.70, -d$URINE_LD[i])
    waterfall(lab, val, "수분 보존 방정식 (L/일) — 마지막 잔차가 0이면 균형",
              "L/일")
  })

  output$t2_energy <- renderPlot({
    r <- simrun(); d <- r$ref; i <- 1
    lab <- c("흡수 열량\nabsorbed", "대장 SCFA\nsalvage", "PN 열량\nPN kcal",
             "총 소비\n-TEE")
    val <- c(d$ABSKC[i], d$SCFAKC[i], d$PNKCALD[i], -d$TEEC[i])
    waterfall(lab, val, "에너지 보존 방정식 (kcal/일)", "kcal/일")
  })

  output$t2_terms <- kv_table({
    r <- simrun(); d <- r$ref; i <- 1
    data.frame(
      `항 (term)` = c("공장 도달 용적 VDEL (L/일)", "관내 [Na⁺] (mmol/L)",
        "공장 흡수 분율 FRACJ", "공장 순 flux (L/일)", "대장 유입 (L/일)",
        "대장 회수 (L/일)", "배출량 (L/일)", "PN 용적 (L/일)",
        "요량 (L/일)", "체수분 편차 (L)"),
      `값` = round(c(d$VORALC[i] + 5.2, d$CNALUM[i], d$FRACJC[i], d$VJEJC[i],
        d$VORALC[i] + 5.2 - d$VJEJC[i], d$VCOLC[i], d$OUT_LD[i],
        d$PNVOLWK[i]/7, d$URINE_LD[i], d$TBWC[i] - 36), 3),
      check.names = FALSE)
  })

  output$t2_na <- renderPlot({
    r <- simrun()
    tsplot(list(`요 나트륨 (mmol/일)` = r$ref, `투약군` = r$drug),
           "UNA", "요 나트륨 (mmol/일)",
           "나트륨 결핍의 가장 이른 지표: 요 Na⁺ < 10",
           cols = c(PAL["grey"], PAL["drug"]), hline = 10)
  })

  # =================== TAB 3 ==========================================
  ors_sweep <- reactive({
    st <- bstate(); p <- runpars()
    of <- seq(0, 0.95, by = 0.0950)
    do.call(rbind, lapply(of, function(x) {
      d <- sim(st, modifyList(p, list(ORSFRAC = x, WEANON = 0)), end = 45)
      j <- nrow(d)
      data.frame(ORS = x, CNAD = d$CNAD[j], DRIVEM = d$DRIVEM[j],
                 DRIVED = d$DRIVED[j], DRIVE = d$DRIVEC[j],
                 FRACJ = d$FRACJC[j], OUT = d$OUT_LD[j], VJEJ = d$VJEJC[j])
    }))
  })

  output$t3_drive <- renderPlot({
    s <- ors_sweep()
    par(mar = c(4.4, 4.8, 3.2, 1.2))
    plot(s$ORS, s$DRIVED, type = "b", pch = 19, lwd = 2.2, col = PAL["warn"],
         ylim = range(c(s$DRIVED, s$DRIVEM, 0)), las = 1,
         xlab = "식사 사이 유체의 ORS 비율",
         ylab = "구동력 [Na⁺]lum − CEQ (mmol/L)",
         main = "흐름별 구동력: 부호가 바뀌는 곳은 오직 stream D",
         cex.main = 1.0, font.main = 1)
    grid(col = "#ececec", lty = 1); abline(h = 0, lwd = 2, col = "#333333")
    lines(s$ORS, s$DRIVEM, type = "b", pch = 17, lwd = 2.2, col = PAL["alt"])
    lines(s$ORS, s$DRIVE,  type = "b", pch = 15, lwd = 2.2, col = PAL["base"])
    text(0.06, 6, "흡수 absorbs", cex = 0.8, col = "#333333", adj = 0)
    text(0.06, -8, "분비 SECRETES", cex = 0.8, col = PAL["warn"], adj = 0)
    legend("bottomright", c("stream D (식사 사이 음수)", "stream M (식사 + 분비)",
                            "부피 가중 평균"),
           col = c(PAL["warn"], PAL["alt"], PAL["base"]), pch = c(19, 17, 15),
           lwd = 2.2, bty = "n", cex = 0.84)
  })

  output$t3_output <- renderPlot({
    s <- ors_sweep()
    par(mar = c(4.4, 4.8, 3.2, 4.6))
    plot(s$ORS, s$OUT, type = "b", pch = 19, lwd = 2.2, col = PAL["orange"],
         las = 1, xlab = "ORS 비율", ylab = "배출량 (L/일)",
         main = "같은 음수량, 조성만 바꿀 때의 배출량과 공장 flux",
         cex.main = 1.0, font.main = 1)
    grid(col = "#ececec", lty = 1)
    par(new = TRUE)
    plot(s$ORS, s$VJEJ, type = "b", pch = 17, lwd = 2.2, col = PAL["base"],
         axes = FALSE, xlab = "", ylab = "")
    axis(4, las = 1); mtext("공장 순 수분 flux (L/일)", side = 4, line = 3.2,
                            cex = 0.9)
    legend("topright", c("배출량 (좌)", "공장 flux (우)"),
           col = c(PAL["orange"], PAL["base"]), pch = c(19, 17), lwd = 2.2,
           bty = "n", cex = 0.84)
  })

  output$t3_tab <- kv_table({
    s <- ors_sweep()
    data.frame(
      `ORS 비율` = round(s$ORS, 2),
      `음수 [Na⁺]` = round(s$CNAD, 1),
      `stream D 구동력` = round(s$DRIVED, 1),
      `방향` = ifelse(s$DRIVED < 0, "분비 (secretes)", "흡수 (absorbs)"),
      `FRACJ` = round(s$FRACJ, 3),
      `배출량 L/일` = round(s$OUT, 2),
      check.names = FALSE)
  })

  output$t3_pn <- renderPlot({
    st <- bstate(); p <- runpars(); d <- max(84, input$days)
    a <- sim(st, modifyList(p, list(ORSFRAC = 0)),    end = d)
    b <- sim(st, modifyList(p, list(ORSFRAC = 0.45)), end = d)
    cc<- sim(st, modifyList(p, list(ORSFRAC = 0.90)), end = d)
    tsplot(list(`맹물만 (ORS 0)` = a, `ORS 0.45` = b, `ORS 0.90` = cc),
           "PNVOLWK", "PN 용적 (L/주)",
           "경구 재수화만 바꿀 때의 PN 궤적",
           cols = c(PAL["warn"], PAL["gold"], PAL["alt"]))
  })

  # =================== TAB 4 ==========================================
  pk_all <- reactive({
    st <- bstate(); p <- runpars(); w <- wt0(); d <- min(60, input$days)
    m <- init(param(mod, p), as.list(st))
    f <- function(ee) as.data.frame(mrgsim(m, ee, end = d, delta = 0.02))
    list(`테두글루타이드 od` = f(ev_ted(w, d)),
         `아프라글루타이드 qw` = f(ev_apra(d)),
         `글레파글루타이드 qw` = f(ev_glep(d, ii = 7)),
         `천연 GLP-2 등몰 od` = f(ev_native(w, d, 1)))
  })

  output$t4_conc <- renderPlot({
    l <- pk_all()
    dl <- list(`테두글루타이드` = l[[1]], `아프라글루타이드` = l[[2]],
               `글레파글루타이드` = l[[3]])
    tsplot(dl, "CPTED", "테두글루타이드 농도 (mg/L)",
           "SC 투여 후 농도 (테두글루타이드 축)",
           cols = c(PAL["drug"], PAL["base"], PAL["alt"]))
  })

  output$t4_occ <- renderPlot({
    l <- pk_all()
    tsplot(l, "OCCGC", "GLP2R 점유율 (0-1)",
           "수용체 점유율 — 효과부위 구획이 시간 규모를 통합",
           cols = c(PAL["drug"], PAL["base"], PAL["alt"], PAL["warn"]))
  })

  output$t4_zoom <- renderPlot({
    l <- pk_all()
    dl <- lapply(l, function(d) d[d$time <= 3, ])
    tsplot(dl, "OCCGC", "GLP2R 점유율",
           "첫 3일 확대: 천연 GLP-2는 t½ 7분이라 점유율이 유지되지 않는다",
           cols = c(PAL["drug"], PAL["base"], PAL["alt"], PAL["warn"]))
  })

  output$t4_tab <- kv_table({
    l <- pk_all()
    do.call(rbind, lapply(names(l), function(n) {
      d <- l[[n]]; tail30 <- d[d$time >= max(d$time) - 7, ]
      data.frame(`요법` = n,
        `평균 점유율` = round(mean(d$OCCGC), 3),
        `최고 점유율` = round(max(tail30$OCCGC), 3),
        `최저 점유율` = round(min(tail30$OCCGC), 3),
        `최저/최고` = round(min(tail30$OCCGC)/max(max(tail30$OCCGC), 1e-9), 3),
        check.names = FALSE)
    }))
  })

  # =================== TAB 5 ==========================================
  output$t5_vill <- renderPlot({
    r <- simrun()
    tsplot(list(`무투약` = r$ref, `투약군` = r$drug), "VILLC",
           "유효 융모 높이 (상대값)", "구조적 적응",
           cols = c(PAL["grey"], PAL["drug"]))
  })
  output$t5_trans <- renderPlot({
    r <- simrun()
    tsplot(list(`무투약` = r$ref, `투약군` = r$drug), "MUCRELC",
           "점막 질량 지수 MUCREL", "점막 질량 (시트룰린의 결정 인자)",
           cols = c(PAL["grey"], PAL["drug"]))
  })
  output$t5_gate <- renderPlot({
    st <- bstate(); p <- runpars(); w <- wt0(); d <- input$days
    out <- lapply(c(1.0, 0.5, 0.15, 0.02), function(ef) {
      s <- burnin(modifyList(p, list(ENTFRAC = ef)), days = 600)
      q <- modifyList(p, list(ENTFRAC = ef, URBASE = s[["URINESM"]],
                              WTTARGET = s[["LEAN"]] + s[["FATM"]]))
      sim(s, q, ev_ted(s[["LEAN"]] + s[["FATM"]], d), end = d)
    })
    names(out) <- c("경장 100%", "경장 50%", "경장 15%", "완전 금식 2%")
    tsplot(out, "TROPHC", "영양 신호 TROPHIC",
           "영양 게이팅: 곱이므로 금식하면 영양성 팔이 멈춘다",
           cols = c(PAL["alt"], PAL["gold"], PAL["orange"], PAL["warn"]))
  })
  output$t5_contact <- renderPlot({
    r <- simrun()
    tsplot(list(`무투약` = r$ref, `투약군` = r$drug), "CONTC",
           "점막 접촉 시간 (상대값)",
           "회장 브레이크를 통한 접촉 시간 — 모든 흡수 분율의 곱셈 인자",
           cols = c(PAL["grey"], PAL["drug"]))
  })

  # =================== TAB 6 ==========================================
  output$t6_pn <- renderPlot({
    r <- simrun()
    tsplot(list(`무투약 (protocol only)` = r$ref, `선택한 요법` = r$drug),
           "PNVOLWK", "PN 용적 (L/주)",
           "PN 용적 — 알고리즘이 만들어 낸 궤적",
           cols = c(PAL["grey"], PAL["drug"]),
           hline = 0.8 * r$ref$PNVOLWK[1])
  })
  output$t6_nights <- renderPlot({
    r <- simrun()
    tsplot(list(`무투약` = r$ref, `선택한 요법` = r$drug), "PNNIGHTS",
           "주당 주입 야간 수", "리터가 아니라 야간이 줄어야 카테터-일이 준다",
           cols = c(PAL["grey"], PAL["drug"]))
  })
  output$t6_tab <- kv_table({
    r <- simrun(); b0 <- r$ref$PNVOLWK[1]; n <- nrow(r$drug)
    wk20 <- which.min(abs(r$drug$time - 140)); wk24 <- which.min(abs(r$drug$time - 168))
    red <- function(d, i) 100 * (b0 - d$PNVOLWK[i]) / b0
    data.frame(
      `엔드포인트` = c("기저 PN 용적 (L/주)", "최종 PN 용적 (L/주)",
        "PN 용적 변화 (L/주)", "PN 감소율 (%)",
        "20주 감소율 (%)", "24주 감소율 (%)",
        "≥20% 반응자 (20주 AND 24주)", "주당 야간 수 변화",
        "누적 카테터-일", "연간 CRBSI"),
      `투약군` = c(round(b0, 2), round(r$drug$PNVOLWK[n], 2),
        round(r$drug$PNVOLWK[n] - b0, 2), round(red(r$drug, n), 1),
        round(red(r$drug, wk20), 1), round(red(r$drug, wk24), 1),
        ifelse(resp20(r$drug, b0) == 1, "예 (YES)", "아니오"),
        round(r$drug$PNNIGHTS[n] - r$drug$PNNIGHTS[1], 2),
        round(r$drug$CRBSICUM[n] * 0 + r$drug$time[n] * r$drug$PNNIGHTS[n]/7, 0),
        round(r$drug$CRBSIYR[n], 3)),
      `무투약` = c(round(b0, 2), round(r$ref$PNVOLWK[n], 2),
        round(r$ref$PNVOLWK[n] - b0, 2), round(red(r$ref, n), 1),
        round(red(r$ref, wk20), 1), round(red(r$ref, wk24), 1),
        ifelse(resp20(r$ref, b0) == 1, "예 (YES)", "아니오"),
        round(r$ref$PNNIGHTS[n] - r$ref$PNNIGHTS[1], 2),
        round(r$ref$time[n] * r$ref$PNNIGHTS[n]/7, 0),
        round(r$ref$CRBSIYR[n], 3)),
      check.names = FALSE)
  })
  output$t6_out <- renderPlot({
    r <- simrun()
    tsplot(list(`무투약` = r$ref, `선택한 요법` = r$drug), "OUT_LD",
           "변·장루 배출량 (L/일)", "배출량 — 측정 가능한 양",
           cols = c(PAL["grey"], PAL["drug"]))
  })

  # =================== TAB 7 ==========================================
  cmp_runs <- eventReactive(list(input$go, input$cmp), {
    st <- bstate(); p <- runpars(); w <- wt0(); d <- input$days
    sel <- input$cmp %||% "none"
    withProgress(message = "시나리오 비교 계산...", value = 0.5, {
      out <- lapply(sel, function(k) {
        switch(k,
          none  = sim(st, p, NULL, end = d),
          ted   = sim(st, p, ev_ted(w, d), end = d),
          apra  = sim(st, p, ev_apra(d), end = d),
          glep  = sim(st, p, ev_glep(d, ii = 7), end = d),
          glep2 = sim(st, p, ev_glep(d, ii = 3.5), end = d),
          nat1  = sim(st, p, ev_native(w, d, 1), end = d),
          ors   = sim(st, modifyList(p, list(ORSFRAC = 0.9)), NULL, end = d),
          lop   = sim(st, p, ev_lop(d), end = d),
          ppi   = sim(st, modifyList(p, list(PPIDOSE = 2)), NULL, end = d),
          oct   = sim(st, p, ev_oct(d), end = d),
          combo = sim(st, modifyList(p, list(ORSFRAC = 0.9)),
                      c(ev_ted(w, d), ev_lop(d)), end = d))
      })
      names(out) <- sel; out
    })
  }, ignoreNULL = FALSE)

  LBL <- c(none = "무투약", ted = "테두글루타이드", apra = "아프라글루타이드 qw",
           glep = "글레파글루타이드 qw", glep2 = "글레파글루타이드 2x/wk",
           nat1 = "천연 GLP-2 od", ors = "ORS 0.9", lop = "로페라마이드",
           ppi = "고용량 PPI", oct = "옥트레오타이드",
           combo = "테두 + ORS + 로페라마이드")

  output$t7_pn <- renderPlot({
    l <- cmp_runs(); if (!length(l)) return(invisible())
    names(l) <- LBL[names(l)]
    tsplot(l, "PNVOLWK", "PN 용적 (L/주)", "시나리오별 PN 궤적")
  })
  output$t7_bar <- renderPlot({
    l <- cmp_runs(); if (!length(l)) return(invisible())
    v <- vapply(l, function(d) d$PNVOLWK[nrow(d)] - d$PNVOLWK[1], numeric(1))
    names(v) <- LBL[names(l)]
    barplot2(v, "최종 PN 용적 변화 (L/주)", "L/주",
             cols = ifelse(v < 0, PAL["alt"], PAL["warn"]), horiz = TRUE)
  })
  output$t7_tab <- kv_table({
    l <- cmp_runs(); if (!length(l)) return(data.frame())
    do.call(rbind, lapply(names(l), function(k) {
      d <- l[[k]]; n <- nrow(d)
      data.frame(`시나리오` = LBL[[k]],
        `PN L/주` = round(d$PNVOLWK[n], 2),
        `변화` = round(d$PNVOLWK[n] - d$PNVOLWK[1], 2),
        `야간/주` = round(d$PNNIGHTS[n], 2),
        `배출 L/일` = round(d$OUT_LD[n], 2),
        `지방 흡수 %` = round(100*d$FFATC[n], 1),
        `융모` = round(d$VILLC[n], 3),
        `시트룰린` = round(d$CITR[n], 1),
        `체중 kg` = round(d$WT[n], 1),
        `≥20% 반응` = ifelse(resp20(d, d$PNVOLWK[1]) == 1, "예", "-"),
        check.names = FALSE)
    }))
  })

  # =================== TAB 8 ==========================================
  output$t8_citr <- renderPlot({
    r <- simrun()
    tsplot(list(`무투약` = r$ref, `투약군` = r$drug), "CITR",
           "혈장 시트룰린 (µmol/L)",
           "시트룰린 — 장세포 질량의 지표 (<20 = 영구 장부전 시사)",
           cols = c(PAL["grey"], PAL["drug"]), hline = 20)
  })
  output$t8_mg <- renderPlot({
    r <- simrun()
    tsplot(list(`혈청 Mg (무투약)` = r$ref, `혈청 Mg (투약군)` = r$drug),
           "SMGC", "혈청 마그네슘 (mmol/L)",
           "마그네슘 — 배출량에 비례해 소실되고 골에서 완충된다",
           cols = c(PAL["grey"], PAL["drug"]), hline = 0.70)
  })
  output$t8_acid <- renderPlot({
    r <- simrun()
    tsplot(list(`중탄산 (무투약)` = r$ref, `중탄산 (투약군)` = r$drug),
           "HCO3C", "혈청 중탄산 (mmol/L)",
           "만성 대사성 산증 — 골 소실과 결석의 공통 경로",
           cols = c(PAL["grey"], PAL["drug"]), hline = 22)
  })
  output$t8_micro <- renderPlot({
    r <- simrun(); d <- r$ref; n <- nrow(d)
    v <- c(`비타민 D` = d$VITDC[n], `아연 (÷1400)` = d$ZNC[n]/1400,
           `B12 (÷3000)` = d$B12C[n]/3000, `필수지방산 (÷450)` = d$EFAC[n]/450,
           `PTH (÷45)` = d$PTHC[n]/45, `골밀도` = d$BMDC[n])
    barplot2(v, "미량영양소·골 상태 (1.0 = 기준)", "기준 대비 비율",
             cols = ifelse(v < 0.7, PAL["warn"], PAL["base"]))
  })

  # =================== TAB 9 ==========================================
  long_runs <- eventReactive(input$go, {
    st <- bstate(); p <- runpars(); w <- wt0()
    withProgress(message = "2년 궤적 계산...", value = 0.5, {
      list(
        `대두유 ILE` = sim(st, modifyList(p, list(PHYTOCONC = 350, OMEGA3 = 0)),
                          NULL, end = 730),
        `SMOF ILE` = sim(st, modifyList(p, list(PHYTOCONC = 100, OMEGA3 = 1)),
                          NULL, end = 730),
        `어유 ILE` = sim(st, modifyList(p, list(PHYTOCONC = 0, OMEGA3 = 1)),
                          NULL, end = 730),
        `대두유 + 테두글루타이드` = sim(st, modifyList(p, list(PHYTOCONC = 350,
                          OMEGA3 = 0)), ev_ted(w, 730), end = 730),
        `대두유 + 타우롤리딘 잠금` = sim(st, modifyList(p, list(PHYTOCONC = 350,
                          OMEGA3 = 0, LOCKFAC = 0.25, TECHFAC = 0.7)),
                          NULL, end = 730))
    })
  }, ignoreNULL = FALSE)

  output$t9_liver <- renderPlot({
    tsplot(long_runs(), "BILIC", "총 빌리루빈 (mg/dL)",
           "IFALD: 지질 유탁액의 조성이 곧 용량", hline = 2)
  })
  output$t9_cvc <- renderPlot({
    tsplot(long_runs(), "VEINSC", "사용 가능 중심정맥 수",
           "접근로는 소모성 자원이다")
  })
  output$t9_renal <- renderPlot({
    tsplot(long_runs(), "GFRC", "eGFR (mL/min/1.73m²)",
           "신기능 — 만성 탈수와 옥살산 부담의 합")
  })
  output$t9_bone <- renderPlot({
    tsplot(long_runs(), "UOXC", "요 옥살산 (mg/일)",
           "장성 고옥살산뇨 — 대장이 있어야 발생한다", hline = 40)
  })
  output$t9_tab <- kv_table({
    l <- long_runs()
    do.call(rbind, lapply(names(l), function(k) {
      d <- l[[k]]; n <- nrow(d)
      data.frame(`2년 시점` = k,
        `PN L/주` = round(d$PNVOLWK[n], 2),
        `빌리루빈` = round(d$BILIC[n], 2),
        `섬유화 단계` = round(d$FIBC[n], 2),
        `누적 CRBSI` = round(d$CRBSICUM[n], 2),
        `남은 정맥` = round(d$VEINSC[n], 2),
        `eGFR` = round(d$GFRC[n], 1),
        `결석 부담` = round(d$STONEC[n], 3),
        `골밀도` = round(d$BMDC[n], 3),
        `생존 %` = round(d$SURVP[n], 1),
        check.names = FALSE)
    }))
  })

  # =================== TAB 10 =========================================
  output$t10_loop <- renderPlot({
    r <- simrun(); d <- r$drug
    par(mar = c(4.4, 4.8, 3.2, 4.8))
    plot(d$time, d$URSM / r$p$URBASE, type = "l", lwd = 2.2, col = PAL["base"],
         las = 1, xlab = "일 (days)", ylab = "요량 / 기저 요량",
         main = "폐쇄 루프: 요량이 트리거를 넘으면 감량, 넘지 않으면 정지",
         cex.main = 1.0, font.main = 1)
    grid(col = "#ececec", lty = 1)
    abline(h = input$urtrig, col = PAL["warn"], lty = 2, lwd = 2)
    abline(h = 1, col = "#999999", lty = 3)
    par(new = TRUE)
    plot(d$time, d$PNVOLWK, type = "l", lwd = 2.2, col = PAL["drug"],
         axes = FALSE, xlab = "", ylab = "")
    axis(4, las = 1); mtext("PN 용적 (L/주)", side = 4, line = 3.3, cex = 0.9)
    legend("topright", c("요량 비 (좌)", "감량 트리거", "PN 용적 (우)"),
           col = c(PAL["base"], PAL["warn"], PAL["drug"]),
           lty = c(1, 2, 1), lwd = 2.2, bty = "n", cex = 0.84)
  })

  output$t10_placebo <- renderPlot({
    p <- basepars(); w <- NULL
    pre <- modifyList(p, list(ORSFRAC = 0.25))
    s <- burnin(pre, days = 700); wi <- s[["LEAN"]] + s[["FATM"]]
    q <- modifyList(p, list(URBASE = s[["URINESM"]], WTTARGET = wi))
    d <- max(168, input$days)
    a <- sim(s, modifyList(q, list(ORSFRAC = 0.60)), NULL, end = d)
    b <- sim(s, modifyList(q, list(ORSFRAC = 0.60)), ev_ted(wi, d), end = d)
    cc<- sim(s, modifyList(q, list(ORSFRAC = 0.25, WEANON = 0)), NULL, end = d)
    tsplot(list(`위약 + 프로토콜화된 진료` = a, `테두글루타이드 + 프로토콜` = b,
                `프로토콜 없음, 약도 없음` = cc),
           "PNVOLWK", "PN 용적 (L/주)",
           "위약군 반응의 출처: 위약이 아니라 프로토콜",
           cols = c(PAL["gold"], PAL["drug"], PAL["grey"]))
  })

  output$t10_trig <- renderPlot({
    st <- bstate(); p <- runpars(); w <- wt0(); d <- max(168, input$days)
    tr <- c(1.03, 1.05, 1.10, 1.20, 1.30)
    v <- vapply(tr, function(x) {
      dd <- sim(st, modifyList(p, list(URTRIG = x)), ev_ted(w, d), end = d)
      dd$PNVOLWK[nrow(dd)] - dd$PNVOLWK[1]
    }, numeric(1))
    names(v) <- sprintf("트리거 %.2f", tr)
    barplot2(v, "감량 트리거 역치가 측정된 약효를 바꾼다", "PN 변화 (L/주)",
             cols = PAL["purple"])
  })

  output$t10_tab <- kv_table({
    r <- simrun(); d <- r$drug; n <- nrow(d)
    data.frame(
      `제어기 상태` = c("기저 요량 URBASE (L/일)", "현재 평활 요량 (L/일)",
        "요량 비", "감량 트리거", "재증량 역치", "최대 감량 속도 (1/일)",
        "체중 가드레일 (kg)", "현재 체중 (kg)", "탈수 지표", "갈증 반응"),
      `값` = round(c(r$p$URBASE, d$URSM[n], d$URSM[n]/r$p$URBASE,
        r$p$URTRIG, 0.85, r$p$KTAPER, 0.95 * r$p$WTTARGET, d$WT[n],
        d$DEHYDC[n], d$THIRSTC[n]), 3),
      check.names = FALSE)
  })
}

shinyApp(ui, server)
