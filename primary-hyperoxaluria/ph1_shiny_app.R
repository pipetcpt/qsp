# =============================================================================
#  Primary Hyperoxaluria Type 1 (PH1) — Shiny dashboard for the QSP model
#  ph1_shiny_app.R   ·   companion to ph1_mrgsolve_model.R / ph1_qsp_model.dot
#
#  Run:  Rscript -e 'shiny::runApp("ph1_shiny_app.R", port = 8080)'
#
#  The app loads the model source out of ph1_mrgsolve_model.R (the `code`
#  string), so the dashboard and the batch script can never drift apart.
#
#  NINE TABS
#    1  환자·유전형        patient, AGXT genotype, starting kidney damage
#    2  약물 PK·표적관여    siRNA PK, RISC loading, mRNA/enzyme knockdown, PLP
#    3  대사 플럭스        the glyoxylate futile cycle, drawn as fluxes
#    4  요 생화학·과포화    urine oxalate/glycolate/citrate/volume and CaOx RSS
#    5  혈장·전신 옥살로시스 plasma oxalate, bone capacitor, organ deposition
#    6  신기능·진행        eGFR, nephrocalcinosis, fibrosis, nephron mass
#    7  임상 엔드포인트     trial-style % change, stone events, ESKD-free survival
#    8  시나리오 비교       the 28 built-in scenarios side by side
#    9  투석·이식 탐색      dialysis intensity and transplant strategy explorer
# =============================================================================

library(shiny)
library(mrgsolve)   # NOTE: mrgsolve also exports req(); because it is attached
                    # after shiny, shiny::req must be called qualified below.
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# --- load the model definition from the batch script -------------------------
model_file <- "ph1_mrgsolve_model.R"
if (!file.exists(model_file)) {
  stop("ph1_mrgsolve_model.R must sit next to this app (the model code is read from it).")
}
src <- readLines(model_file, warn = FALSE)
i0 <- grep("^code <- '", src)[1]
i1 <- grep("^'$", src)
i1 <- i1[i1 > i0][1]
code <- paste(c(sub("^code <- '", "", src[i0]), src[(i0 + 1):(i1 - 1)]), collapse = "\n")
mod <- mcode_cache("ph1_app", code)

THEME <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "#eef2f7"),
        plot.title = element_text(face = "bold", size = 13))

# Korean display labels are kept in STRINGS, never used as R identifiers, so
# the app parses in any locale (a C-locale Rscript cannot parse non-ASCII
# syntactic names).  ASCII column names are mapped to labels at plot time.
LBL <- c(
  AGTact = "AGT 활성", GOact = "GO 활성", LDHact = "LDH 활성",
  LUMC = "루마시란 혈장 (mg)", LUMR = "루마시란 간 RISC",
  NEDC = "네도시란 혈장 (mg)", NEDR = "네도시란 간 RISC",
  HAO1M = "HAO1 mRNA", GOP = "GO 단백", LDHAM = "LDHA mRNA", LDHP = "LDH 단백",
  Uox24 = "요 옥살산 (mmol/일)", Uglc24 = "요 글리콜산 (mmol/일)",
  Pox = "혈장 옥살산 (µmol/L)", Pglc = "혈장 글리콜산 (µmol/L)",
  Pglx = "혈장 글리옥실산 (µmol/L)", RSSr = "CaOx 상대과포화도",
  Uvol24 = "요량 (L/일)", Ucit24 = "요 구연산 (mmol/일)",
  UoxCr = "요 옥살산:크레아티닌",
  BoneOx = "뼈 저장 (mmol)", SoftOx = "연조직 (mmol)",
  RETINA = "망막", CARDIO = "심장", BONEDIS = "골질환", NEURO = "신경",
  eGFR = "eGFR (mL/min/1.73m²)", NCgrade = "신석회화 등급 (0–3)",
  FIBgrade = "섬유화 등급 (0–3)", StoneEv = "결석 사건 (누적)",
  SurvESKD = "신부전 무발생 생존", load = "신단위당 부하 (정규화)"
)
relab <- function(d) { d$name <- ifelse(is.na(LBL[d$name]), d$name, LBL[d$name]); d }
long <- function(d, cols) {
  relab(pivot_longer(d[, c("time", cols)], -time))
}
named <- function(df, nms) { names(df) <- nms; df }

GENOTYPES <- list(
  "G170R/G170R — 미스타게팅, 피리독신 반응성"        = list(AGTACT0 = 0.02,  B6RESP = 1),
  "F152I — 접힘결함, 부분 반응성"                   = list(AGTACT0 = 0.03,  B6RESP = 0.7),
  "I244T — 부분 반응성"                            = list(AGTACT0 = 0.04,  B6RESP = 0.5),
  "c.33dupC/c.33dupC — 절단형, 비반응성"            = list(AGTACT0 = 0.005, B6RESP = 0),
  "복합 이형접합 (중간 활성)"                        = list(AGTACT0 = 0.08,  B6RESP = 0.4),
  "정상 AGT (건강 대조)"                            = list(AGTACT0 = 1.0,   B6RESP = 0)
)

CKD_PRESET <- list(
  "CKD1–2 (조기 진단)"        = list(NM0 = 0.95, NC0 = 0.00, FIB0 = 0.00, OXABONE0 = 0),
  "CKD3a"                    = list(NM0 = 0.55, NC0 = 0.25, FIB0 = 0.30, OXABONE0 = 40000),
  "CKD4 (지연 진단)"          = list(NM0 = 0.30, NC0 = 0.45, FIB0 = 0.50, OXABONE0 = 150000),
  "ESKD (투석 중)"            = list(NM0 = 0.04, NC0 = 0.80, FIB0 = 0.80, OXABONE0 = 350000)
)

# The 28 scenarios of the batch script, reproduced so tab 8 matches it exactly.
CKD4  <- list(NM0 = 0.30, NC0 = 0.45, FIB0 = 0.50, OXABONE0 = 150000)
ESKDp <- list(NM0 = 0.04, NC0 = 0.80, FIB0 = 0.80, OXABONE0 = 350000, HDON = 1)
CONS  <- list(UVTGT = 3.0, KCITON = 1)
SCEN <- list(
  "01 건강 대조"                            = list(AGTACT0 = 1.0, NM0 = 1.0),
  "02 PH1 미치료 (전형적)"                   = list(),
  "03 PH1 미치료 (null 유전형)"              = list(AGTACT0 = 0.005, B6RESP = 0),
  "04 영아형 옥살로시스"                     = list(AGTACT0 = 0.005, B6RESP = 0, BSA = 0.45,
                                                  WT = 8, GFR0 = 30, PRECUR = 1.6, NM0 = 0.88),
  "05 수분섭취 3 L 단독"                     = list(UVTGT = 3.0),
  "06 수분섭취 + 구연산칼륨"                  = CONS,
  "07 수분섭취 순응도 50%"                   = c(CONS, list(ADHFL = 0.5)),
  "08 피리독신 8 mg/kg (반응형)"             = list(B6ON = 1, B6DOSE = 8),
  "09 피리독신 20 mg/kg (반응형)"            = list(B6ON = 1, B6DOSE = 20),
  "10 피리독신 20 mg/kg (비반응형)"          = list(B6ON = 1, B6DOSE = 20, AGTACT0 = 0.005, B6RESP = 0),
  "11 루마시란"                             = list(LUMON = 1),
  "12 루마시란 + 보존적 치료"                 = c(CONS, list(LUMON = 1)),
  "13 루마시란 영아 용법 (ILLUMINATE-B)"      = list(LUMON = 1, LUMLOAD = 6, LUMDOSE = 3,
                                                  LUMTAU = 30, LUMNL = 3, BSA = 0.55, WT = 10,
                                                  GFR0 = 60, PRECUR = 1.4, NM0 = 0.92),
  "14 루마시란 진행 CKD (ILLUMINATE-C A)"    = c(CKD4, list(LUMON = 1)),
  "15 루마시란 투석 중 (ILLUMINATE-C B)"     = c(ESKDp, list(LUMON = 1)),
  "16 네도시란"                             = list(NEDON = 1),
  "17 네도시란 + 보존적 치료"                 = c(CONS, list(NEDON = 1)),
  "18 이중 RNAi (루마시란+네도시란)"          = list(LUMON = 1, NEDON = 1),
  "19 스티리펜톨 단독"                       = list(STPON = 1),
  "20 경구 옥살산 분해요법"                   = list(OXDON = 1),
  "21 CKD4 지연진단, 기질감소 없음"           = CKD4,
  "22 CKD4 + 루마시란 + 보존적"              = c(CKD4, CONS, list(LUMON = 1)),
  "23 ESKD 통상 주3회 HD"                   = ESKDp,
  "24 ESKD 강화 매일 HD + 야간 PD"           = c(ESKDp, list(HDSESS = 6, PDCL = 6)),
  "25 ESKD → 간·신 동시이식"                 = c(ESKDp, list(TXDAY = 365, KTXON = 1, LTXON = 1)),
  "26 ESKD → 단독 신이식 (RNAi 없음)"        = c(ESKDp, list(TXDAY = 365, KTXON = 1, LTXON = 0)),
  "27 ESKD → 단독 신이식 + 루마시란"          = c(ESKDp, list(TXDAY = 365, KTXON = 1,
                                                           LTXON = 0, LUMON = 1)),
  "28 IL-1 차단 부가 (탐색적)"               = c(CONS, list(ANAON = 1))
)

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel("원발성 과옥살산뇨증 1형 (PH1) QSP 시뮬레이터 · Primary Hyperoxaluria Type 1"),
  tags$p(style = "color:#555;margin-top:-8px",
         "간 AGT 결손 → 글리옥실산 무익회로 → 옥살산 → 결정 신병증 → 전신 옥살로시스. ",
         tags$b("교육·연구 목적 모델이며 임상 의사결정에 사용할 수 없습니다.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 · 유전형"),
      selectInput("geno", "AGXT 유전형", names(GENOTYPES), selected = names(GENOTYPES)[1]),
      selectInput("ckd", "진단 시점 신장 상태", names(CKD_PRESET), selected = names(CKD_PRESET)[1]),
      sliderInput("wt", "체중 (kg)", 5, 100, 70, 1),
      sliderInput("bsa", "체표면적 (m²)", 0.3, 2.2, 1.73, 0.01),
      sliderInput("gfr0", "완전 신단위에서의 GFR (mL/min/1.73 m²)", 25, 130, 105, 5),
      sliderInput("precur", "전구체 공급 배수 (콜라겐 전환)", 0.7, 2.0, 1.0, 0.05),
      hr(),
      h4("기질 감소 요법 (RNAi)"),
      checkboxInput("lum", "루마시란 (HAO1 siRNA)", FALSE),
      conditionalPanel("input.lum",
        sliderInput("lumdose", "유지 용량 (mg/kg SC)", 0.5, 6, 3, 0.5),
        sliderInput("lumtau", "유지 투여간격 (일)", 30, 180, 90, 30)),
      checkboxInput("ned", "네도시란 (LDHA siRNA)", FALSE),
      conditionalPanel("input.ned",
        sliderInput("neddose", "용량 (mg/kg SC)", 0.5, 6, 3.5, 0.5),
        sliderInput("nedtau", "투여간격 (일)", 15, 60, 30, 15)),
      hr(),
      h4("보존적 · 기타 치료"),
      sliderInput("uvtgt", "목표 요량 (L/1.73 m²/일)", 1.0, 4.0, 1.5, 0.1),
      sliderInput("adh", "수분섭취 순응도", 0.3, 1.0, 1.0, 0.05),
      checkboxInput("kcit", "구연산칼륨", FALSE),
      checkboxInput("b6", "피리독신 (비타민 B6)", FALSE),
      conditionalPanel("input.b6",
        sliderInput("b6dose", "용량 (mg/kg/일)", 2, 30, 8, 1)),
      checkboxInput("stp", "스티리펜톨 (LDH-5 억제)", FALSE),
      checkboxInput("oxd", "경구 옥살산 분해 (렐록살리아제 / O. formigenes)", FALSE),
      checkboxInput("ana", "IL-1 차단 (탐색적)", FALSE),
      hr(),
      h4("신대체요법 · 이식"),
      checkboxInput("hd", "혈액투석 강제 시작", FALSE),
      sliderInput("hdsess", "주당 투석 횟수", 0, 7, 3, 1),
      sliderInput("hdhr", "회당 시간 (h)", 2, 8, 4, 0.5),
      sliderInput("pdcl", "추가 복막투석 청소율 (L/일)", 0, 10, 0, 1),
      selectInput("tx", "이식 전략",
                  c("없음", "단독 신이식", "단독 신이식 + RNAi 유지", "간·신 동시이식")),
      conditionalPanel("input.tx != '없음'",
        sliderInput("txday", "이식 시점 (일)", 200, 3650, 365, 50)),
      hr(),
      sliderInput("years", "시뮬레이션 기간 (년)", 2, 25, 15, 1),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary", width = "100%")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · 환자·유전형",
                 br(), h4("설정 요약"), tableOutput("profile"),
                 h4("이 유전형에서 무엇이 결정되는가"),
                 plotOutput("p_geno", height = "420px"),
                 tags$p(style = "color:#555",
                        "잔존 AGT 활성 하나가 글리옥실산의 운명을 결정합니다. ",
                        "피리독신은 미스타게팅 변이(G170R/F152I)에서만 AGT를 되살립니다.")),
        tabPanel("2 · 약물 PK·표적관여",
                 br(), plotOutput("p_pk", height = "330px"),
                 plotOutput("p_kd", height = "330px"),
                 tags$p(style = "color:#555",
                        "siRNA 혈장 농도는 수 시간 내 사라지지만 간 내 RISC 적재 풀은 ",
                        "수 주간 유지되므로 분기 투여가 가능합니다. 효과 발현·소실의 지연은 ",
                        "mRNA가 아니라 단백질 반감기(GO ≈ 7일)가 결정합니다.")),
        tabPanel("3 · 대사 플럭스 (무익회로)",
                 br(), plotOutput("p_flux", height = "360px"),
                 h4("탄소 수지 — 어디로 나가는가"), tableOutput("t_flux"),
                 tags$p(style = "color:#555",
                        "AGT가 사라지면 글리신으로 가던 탄소가 GRHPR↔GO 회로를 돌며 ",
                        "LDH에 반복 노출됩니다. GO를 막으면(루마시란) 탄소는 글리콜산으로, ",
                        "LDHA를 막으면(네도시란) 글리옥실산으로 빠져나갑니다.")),
        tabPanel("4 · 요 생화학·과포화",
                 br(), plotOutput("p_urine", height = "380px"),
                 plotOutput("p_rss", height = "300px"),
                 tags$p(style = "color:#555",
                        "결석을 만드는 것은 요 옥살산 자체가 아니라 상대과포화도(RSS)입니다. ",
                        "수분섭취와 구연산은 옥살산 생성을 전혀 바꾸지 않고 RSS만 낮춥니다.")),
        tabPanel("5 · 혈장·전신 옥살로시스",
                 br(), plotOutput("p_pox", height = "340px"),
                 plotOutput("p_sysox", height = "340px"),
                 tags$p(style = "color:#555",
                        "GFR이 떨어지면 혈장 옥살산은 선형이 아니라 쌍곡선으로 상승합니다. ",
                        "30 µmol/L을 넘으면 뼈가 저장고(capacitor) 역할을 하며, 이식 후에는 ",
                        "반대로 수개월–수년간 방출합니다.")),
        tabPanel("6 · 신기능·진행",
                 br(), plotOutput("p_kidney", height = "400px"),
                 plotOutput("p_loop2", height = "300px"),
                 tags$p(style = "color:#555",
                        "신단위가 줄면 남은 신단위 하나가 감당하는 옥살산 부하가 늘어납니다. ",
                        "요 옥살산이 감소해도 신단위당 부하는 증가하는 것이 가속 기전입니다.")),
        tabPanel("7 · 임상 엔드포인트",
                 br(), h4("시험 형식 요약 (기저 = 180일, 6개월 = 360일)"),
                 tableOutput("t_end"),
                 plotOutput("p_end", height = "360px")),
        tabPanel("8 · 시나리오 비교",
                 br(),
                 checkboxGroupInput("scen", "비교할 시나리오", names(SCEN),
                                    selected = names(SCEN)[c(2, 6, 8, 11, 16, 18)], inline = TRUE),
                 actionButton("goscen", "시나리오 실행", class = "btn-primary"),
                 br(), br(), plotOutput("p_scen", height = "520px"),
                 h4("6개월 시점 비교"), tableOutput("t_scen")),
        tabPanel("9 · 투석·이식 탐색",
                 br(), h4("투석 강도 스윕 — 혈장 옥살산을 30 µmol/L 아래로 내릴 수 있는가"),
                 plotOutput("p_hd", height = "360px"),
                 h4("이식 전략 비교"), plotOutput("p_tx", height = "360px"),
                 tags$p(style = "color:#555",
                        "통상 주 3회 투석의 옥살산 제거량은 생성량보다 작습니다. ",
                        "간 이식만이 생성 자체를 정상화하며, 단독 신이식은 대사결손이 ",
                        "남아 이식신에 옥살로시스가 재발합니다."))
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  pars <- reactive({
    p <- c(GENOTYPES[[input$geno]], CKD_PRESET[[input$ckd]])
    p$WT <- input$wt; p$BSA <- input$bsa; p$GFR0 <- input$gfr0; p$PRECUR <- input$precur
    p$UVTGT <- input$uvtgt; p$ADHFL <- input$adh
    p$LUMON <- as.numeric(input$lum); p$NEDON <- as.numeric(input$ned)
    if (input$lum) { p$LUMDOSE <- input$lumdose; p$LUMTAU <- input$lumtau }
    if (input$ned) { p$NEDDOSE <- input$neddose; p$NEDTAU <- input$nedtau }
    p$KCITON <- as.numeric(input$kcit)
    p$B6ON <- as.numeric(input$b6); if (input$b6) p$B6DOSE <- input$b6dose
    p$STPON <- as.numeric(input$stp); p$OXDON <- as.numeric(input$oxd)
    p$ANAON <- as.numeric(input$ana)
    if (input$hd) p$HDON <- 1
    p$HDSESS <- input$hdsess; p$HDHR <- input$hdhr; p$PDCL <- input$pdcl
    if (input$tx != "없음") {
      p$TXDAY <- input$txday; p$KTXON <- 1
      p$LTXON <- as.numeric(input$tx == "간·신 동시이식")
      if (input$tx == "단독 신이식 + RNAi 유지") p$LUMON <- 1
    }
    p
  })

  sim <- eventReactive(input$go, {
    withProgress(message = "시뮬레이션 실행 중 ...", {
      param(mod, pars()) %>%
        mrgsim(end = input$years * 365, delta = 1, hmax = 2,
               atol = 1e-8, rtol = 1e-6) %>%
        as_tibble()
    })
  }, ignoreNULL = FALSE)

  yr <- function(d) d$time / 365.25

  # ---- tab 1 ---------------------------------------------------------------
  output$profile <- renderTable({
    p <- pars()
    named(data.frame(
      item = c("AGXT 유전형", "잔존 AGT 활성", "피리독신 반응성", "진단 시 신장 상태",
               "체중 / BSA", "기준 GFR", "기질감소 요법", "보존적 치료",
               "신대체요법", "이식"),
      value = c(input$geno,
             sprintf("%.1f%% of normal", 100 * p$AGTACT0),
             sprintf("%.0f%%", 100 * p$B6RESP),
             input$ckd,
             sprintf("%.0f kg / %.2f m²", p$WT, p$BSA),
             sprintf("%.0f mL/min/1.73 m²", p$GFR0),
             paste(c(if (input$lum) "루마시란", if (input$ned) "네도시란",
                     if (input$stp) "스티리펜톨", if (input$b6) "피리독신") %||% "없음",
                   collapse = " + "),
             sprintf("요량 목표 %.1f L, 순응도 %.0f%%%s", input$uvtgt, 100 * input$adh,
                     if (input$kcit) ", 구연산칼륨" else ""),
             if (input$hd || input$hdsess > 0) sprintf("주 %d회 × %.1f h (+PD %d L/일)",
                                                       input$hdsess, input$hdhr, input$pdcl) else "없음",
             input$tx)), c("항목", "값"))
  })

  output$p_geno <- renderPlot({
    d <- sim()
    long(d, c("AGTact", "GOact", "LDHact")) %>%
      ggplot(aes(time / 365.25, value, colour = name)) +
      geom_line(linewidth = 0.9) + ylim(0, 1.05) +
      labs(x = "년", y = "정상 대비 활성 (분율)", colour = NULL,
           title = "간 효소 활성 — 질병과 약물이 함께 결정한다") + THEME
  })

  # ---- tab 2 ---------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim()
    long(d, c("LUMC", "LUMR", "NEDC", "NEDR")) %>%
      ggplot(aes(time / 365.25, value)) + geom_line(linewidth = 0.7, colour = "#1f6fb2") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "년", y = NULL, title = "siRNA 약동학 — 혈장은 빠르게, 간 RISC는 느리게") + THEME
  })
  output$p_kd <- renderPlot({
    d <- sim()
    long(d, c("HAO1M", "GOact", "LDHAM", "LDHact", "AGTact")) %>%
      ggplot(aes(time / 365.25, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(x = "년", y = "분율", colour = NULL,
           title = "표적 관여 — mRNA는 즉시, 단백질은 지연") + THEME
  })

  # ---- tab 3 ---------------------------------------------------------------
  output$p_flux <- renderPlot({
    d <- sim()
    long(d, c("Uox24", "Uglc24", "Pglc", "Pglx")) %>%
      ggplot(aes(time / 365.25, value)) + geom_line(linewidth = 0.8, colour = "#b03030") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "년", y = NULL,
           title = "무익회로의 지문 — 어느 대사체가 오르는지가 작용점을 알려준다") + THEME
  })
  output$t_flux <- renderTable({
    d <- sim(); i <- which.min(abs(d$time - 360)); j <- which.min(abs(d$time - 180))
    b <- c(d$Uox24[j], d$Uglc24[j], d$Pox[j], d$Pglc[j], d$Pglx[j])
    m <- c(d$Uox24[i], d$Uglc24[i], d$Pox[i], d$Pglc[i], d$Pglx[i])
    named(data.frame(
      metric = c("요 옥살산 (mmol/1.73m²/일)", "요 글리콜산 (mmol/1.73m²/일)",
                 "혈장 옥살산 (µmol/L)", "혈장 글리콜산 (µmol/L)",
                 "혈장 글리옥실산 (µmol/L)"),
      base = round(b, 3), month6 = round(m, 3),
      pct = round(100 * (m / b - 1), 1)),
      c("지표", "기저 (180일)", "6개월 (360일)", "변화율 (%)"))
  })

  # ---- tab 4 ---------------------------------------------------------------
  output$p_urine <- renderPlot({
    d <- sim()
    long(d, c("Uox24", "Uvol24", "Ucit24", "UoxCr")) %>%
      ggplot(aes(time / 365.25, value)) + geom_line(linewidth = 0.8, colour = "#2b7a3d") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "년", y = NULL, title = "요 생화학") + THEME
  })
  output$p_rss <- renderPlot({
    d <- sim()
    ggplot(d, aes(time / 365.25, RSSr)) +
      geom_hline(yintercept = 4, linetype = 2, colour = "#c0392b") +
      annotate("text", x = max(yr(d)) * 0.75, y = 4.4, label = "준안정 한계 (자발 핵형성)",
               colour = "#c0392b", size = 3.5) +
      geom_line(linewidth = 1, colour = "#c0392b") +
      labs(x = "년", y = "CaOx 상대과포화도",
           title = "과포화도 — 이 선 아래로 내려가면 새 결석이 생기지 않는다") + THEME
  })

  # ---- tab 5 ---------------------------------------------------------------
  output$p_pox <- renderPlot({
    d <- sim()
    ggplot(d, aes(time / 365.25, Pox)) +
      geom_hline(yintercept = 30, linetype = 2, colour = "#6a1b9a") +
      annotate("text", x = max(yr(d)) * 0.7, y = 32.5, label = "전신 옥살로시스 침착 역치",
               colour = "#6a1b9a", size = 3.5) +
      geom_line(linewidth = 1, colour = "#6a1b9a") +
      labs(x = "년", y = "혈장 옥살산 (µmol/L)",
           title = "혈장 옥살산 — GFR이 떨어지면 쌍곡선으로 상승한다") + THEME
  })
  output$p_sysox <- renderPlot({
    d <- sim()
    long(d, c("BoneOx", "SoftOx", "RETINA", "CARDIO", "BONEDIS", "NEURO")) %>%
      ggplot(aes(time / 365.25, value)) + geom_line(linewidth = 0.8, colour = "#6a1b9a") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "년", y = NULL, title = "뼈 저장고와 장기별 침착") + THEME
  })

  # ---- tab 6 ---------------------------------------------------------------
  output$p_kidney <- renderPlot({
    d <- sim()
    long(d, c("eGFR", "NCgrade", "FIBgrade", "StoneEv")) %>%
      ggplot(aes(time / 365.25, value)) + geom_line(linewidth = 0.9, colour = "#1f6fb2") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "년", y = NULL, title = "신기능과 병변 진행") + THEME
  })
  output$p_loop2 <- renderPlot({
    d <- sim() %>% mutate(load = Uox24 / pmax(eGFR / 105, 0.02))
    long(d, c("Uox24", "load")) %>%
      ggplot(aes(time / 365.25, value, colour = name)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "년", y = NULL, colour = NULL,
           title = "가속 기전 — 요 옥살산은 내려가는데 신단위당 부하는 올라간다") +
      THEME + theme(legend.position = "none")
  })

  # ---- tab 7 ---------------------------------------------------------------
  output$t_end <- renderTable({
    d <- sim()
    a <- function(day, col) d[[col]][which.min(abs(d$time - day))]
    ey <- input$years * 365
    named(data.frame(
      endpoint = c("요 옥살산 기저 (mmol/1.73m²/일)", "요 옥살산 6개월",
                     "요 옥살산 변화율 (%)", "혈장 옥살산 기저 (µmol/L)",
                     "혈장 옥살산 6개월", "혈장 옥살산 변화율 (%)",
                     "CaOx RSS 6개월", "eGFR 1년", "eGFR 최종",
                     "신석회화 등급 최종", "누적 결석 사건", "전신 옥살로시스 지표",
                     "말기신질환 도달", "신부전 무발생 생존확률"),
      value = c(sprintf("%.3f", a(180, "Uox24")), sprintf("%.3f", a(360, "Uox24")),
             sprintf("%+.1f", 100 * (a(360, "Uox24") / a(180, "Uox24") - 1)),
             sprintf("%.1f", a(180, "Pox")), sprintf("%.1f", a(360, "Pox")),
             sprintf("%+.1f", 100 * (a(360, "Pox") / a(180, "Pox") - 1)),
             sprintf("%.2f", a(360, "RSSr")), sprintf("%.1f", a(365, "eGFR")),
             sprintf("%.1f", a(ey, "eGFR")), sprintf("%.2f", a(ey, "NCgrade")),
             sprintf("%.1f", a(ey, "StoneEv")), sprintf("%.3f", a(ey, "SYSOX")),
             if (any(d$eGFR < 15)) sprintf("%.1f년", min(d$time[d$eGFR < 15]) / 365.25) else "미도달",
             sprintf("%.3f", a(ey, "SurvESKD")))),
      c("엔드포인트", "값"))
  })
  output$p_end <- renderPlot({
    d <- sim()
    long(d, c("Uox24", "Pox", "eGFR", "SurvESKD")) %>%
      ggplot(aes(time / 365.25, value)) + geom_line(linewidth = 0.9) +
      geom_vline(xintercept = 180 / 365.25, linetype = 3) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "년", y = NULL,
           title = "임상 엔드포인트 (점선 = 치료 시작)") + THEME
  })

  # ---- tab 8 ---------------------------------------------------------------
  scen_sim <- eventReactive(input$goscen, {
    shiny::req(length(input$scen) > 0)
    withProgress(message = "시나리오 실행 중 ...", {
      bind_rows(lapply(input$scen, function(nm) {
        p <- SCEN[[nm]]
        m <- if (length(p)) param(mod, p) else mod
        out <- m %>% mrgsim(end = min(input$years, 15) * 365, delta = 2,
                            hmax = 2, atol = 1e-8, rtol = 1e-6) %>% as_tibble()
        out$scenario <- nm
        out
      }))
    })
  })
  output$p_scen <- renderPlot({
    d <- scen_sim()
    relab(pivot_longer(d[, c("time", "scenario", "Uox24", "Pox", "Pglc",
                             "RSSr", "eGFR", "NCgrade")], -c(time, scenario))) %>%
      ggplot(aes(time / 365.25, value, colour = scenario)) +
      geom_line(linewidth = 0.7) + facet_wrap(~name, scales = "free_y") +
      labs(x = "년", y = NULL, colour = NULL, title = "시나리오 비교") +
      THEME + guides(colour = guide_legend(ncol = 2))
  })
  output$t_scen <- renderTable({
    d <- scen_sim()
    tb <- d %>% group_by(scenario) %>%
      summarise(
        uox_b   = round(Uox24[which.min(abs(time - 180))], 3),
        uox_m6  = round(Uox24[which.min(abs(time - 360))], 3),
        pct     = round(100 * (Uox24[which.min(abs(time - 360))] /
                                Uox24[which.min(abs(time - 180))] - 1), 1),
        uglc_m6 = round(Uglc24[which.min(abs(time - 360))], 3),
        pox_m6  = round(Pox[which.min(abs(time - 360))], 1),
        rss_m6  = round(RSSr[which.min(abs(time - 360))], 2),
        egfr_y5 = round(eGFR[which.min(abs(time - 1825))], 1),
        egfr_end = round(last(eGFR), 1),
        stones  = round(last(StoneEv), 1),
        .groups = "drop")
    named(as.data.frame(tb),
          c("시나리오", "요Ox 기저", "요Ox 6개월", "변화 %", "요Glc 6개월",
            "혈장Ox 6개월", "RSS 6개월", "eGFR 5년", "eGFR 최종", "결석 사건"))
  })

  # ---- tab 9 ---------------------------------------------------------------
  output$p_hd <- renderPlot({
    base <- c(ESKDp, list(WT = input$wt, BSA = input$bsa))
    grid <- expand.grid(sess = c(3, 4, 5, 6, 7), pd = c(0, 4, 8))
    d <- bind_rows(lapply(seq_len(nrow(grid)), function(k) {
      out <- param(mod, c(base, list(HDSESS = grid$sess[k], PDCL = grid$pd[k]))) %>%
        mrgsim(end = 730, delta = 5, hmax = 2, atol = 1e-8, rtol = 1e-6) %>% as_tibble()
      out$sess <- grid$sess[k]; out$pd <- paste0("PD ", grid$pd[k], " L/일"); out
    }))
    d %>% filter(time > 300) %>% group_by(sess, pd) %>%
      summarise(Pox = mean(Pox), .groups = "drop") %>%
      ggplot(aes(sess, Pox, colour = pd)) +
      geom_hline(yintercept = 30, linetype = 2, colour = "#6a1b9a") +
      geom_line(linewidth = 1) + geom_point(size = 2.5) +
      labs(x = "주당 혈액투석 횟수", y = "정상상태 혈장 옥살산 (µmol/L)", colour = NULL,
           title = "투석 강도로 옥살로시스 역치를 넘어설 수 있는가") + THEME
  })
  output$p_tx <- renderPlot({
    strat <- list(
      "이식 없음 (투석 유지)" = ESKDp,
      "단독 신이식" = c(ESKDp, list(TXDAY = 365, KTXON = 1, LTXON = 0)),
      "단독 신이식 + 루마시란" = c(ESKDp, list(TXDAY = 365, KTXON = 1, LTXON = 0, LUMON = 1)),
      "간·신 동시이식" = c(ESKDp, list(TXDAY = 365, KTXON = 1, LTXON = 1)))
    d <- bind_rows(lapply(names(strat), function(nm) {
      out <- param(mod, strat[[nm]]) %>%
        mrgsim(end = 2920, delta = 5, hmax = 2, atol = 1e-8, rtol = 1e-6) %>% as_tibble()
      out$strategy <- nm; out
    }))
    relab(pivot_longer(d[, c("time", "strategy", "eGFR", "Pox", "Uox24",
                             "NCgrade")], -c(time, strategy))) %>%
      ggplot(aes(time / 365.25, value, colour = strategy)) +
      geom_vline(xintercept = 1, linetype = 3) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      labs(x = "년", y = NULL, colour = NULL,
           title = "이식 전략 (점선 = 이식 시점) — 간을 바꿔야 생성이 멈춘다") + THEME
  })
}

`%||%` <- function(a, b) if (length(a) == 0) b else a

shinyApp(ui, server)
