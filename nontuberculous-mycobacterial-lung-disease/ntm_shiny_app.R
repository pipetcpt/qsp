## =============================================================================
##  ntm_shiny_app.R — Interactive dashboard for the MAC-PD QSP model
##  Nontuberculous Mycobacterial Lung Disease (Mycobacterium avium complex)
##
##  Run with:   shiny::runApp("ntm_shiny_app.R")
##  Requires:   shiny, mrgsolve, dplyr, tidyr, ggplot2, DT, bslib
##              plus ntm_mrgsolve_model.R in the same directory
##
##  11 tabs:
##    1  환자 프로파일          Patient profile & host phenotype
##    2  약물 PK                 Drug PK across plasma / ELF / intracellular
##    3  세균 니치               Where the bacteria actually are (4 niches)
##    4  임상 엔드포인트         Culture conversion, QOL-B, BMI, cavity
##    5  내성 진화               rrl mutant selection and the per-niche gate
##    6  시나리오 비교           12 regimens side by side
##    7  바이오마커              Immune / inflammatory / airway readouts
##    8  독성 · 안전성           Ototoxicity, optic, renal, hepatic, QTc
##    9  pH 역설 탐색기          The one input, two consequences, explored
##   10  민감도 분석             Local sensitivity on the conversion endpoint
##   11  모델 정보               Structure, assumptions, references
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

## ---------------------------------------------------------------------------
## Load the model. ntm_mrgsolve_model.R defines `ntm_code` and compiles `mod`;
## we source it in a sandbox environment so the scenario simulations at the
## bottom of that file do not run every time the app starts.
## ---------------------------------------------------------------------------
build_model <- function() {
  src  <- readLines("ntm_mrgsolve_model.R")
  stop_at <- grep("^##  COMPARTMENT NUMBERS FOR DOSING", src)[1]
  if (is.na(stop_at)) stop_at <- length(src)
  e <- new.env()
  eval(parse(text = paste(src[1:(stop_at - 2)], collapse = "\n")), envir = e)
  e$mod
}
mod <- build_model()

CMT <- setNames(seq_along(mrgsolve::cmt(mod)), mrgsolve::cmt(mod))

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#EEF2F7", colour = NA),
        legend.position  = "bottom")

PAL <- c("#0277BD", "#C62828", "#2E7D32", "#EF6C00",
         "#5E35B1", "#00838F", "#AD1457", "#455A64")

## ---------------------------------------------------------------------------
## Regimen builder — turns UI switches into an mrgsolve event object
## ---------------------------------------------------------------------------
build_events <- function(tx_days, macro, macro_freq, emb, rif, cfz,
                         amk_iv, amk_iv_wks, alis, wt) {
  evs <- list()
  add <- function(cmt, amt, ii, until, rate = 0) {
    if (amt <= 0 || until <= 0) return(NULL)
    ev(amt = amt, cmt = cmt, ii = ii, addl = max(0, floor(until / ii)),
       time = 0, rate = rate)
  }
  ii_m <- if (macro_freq == "TIW") 7/3 else if (macro_freq == "BID") 0.5 else 1
  if (macro != "none")
    evs <- c(evs, list(add(CMT["MGUT"],
                           if (macro == "CLR") 500 else if (ii_m > 1) 500 else 250,
                           ii_m, tx_days)))
  if (emb)  evs <- c(evs, list(add(CMT["EGUT"],
                                   if (ii_m > 1) 25 * wt else 15 * wt, ii_m, tx_days)))
  if (rif)  evs <- c(evs, list(add(CMT["FGUT"], 600, ii_m, tx_days)))
  if (cfz)  evs <- c(evs, list(add(CMT["CGUT"], 100, 1, tx_days)))
  if (amk_iv) evs <- c(evs, list(add(CMT["KCEN"], 15 * wt, 7/3,
                                     amk_iv_wks * 7, rate = 15 * wt * 24)))
  if (alis) evs <- c(evs, list(add(CMT["KLIP"], 590, 1, tx_days)))
  evs <- Filter(Negate(is.null), evs)
  if (!length(evs)) return(NULL)
  Reduce(c, evs)
}

run_sim <- function(pars, events, end = 540) {
  s <- mod %>% param(pars)
  out <- if (is.null(events)) mrgsim(s, end = end, delta = 1)
         else mrgsim(s, events = events, end = end, delta = 1)
  as_tibble(out)
}

## ---------------------------------------------------------------------------
## Pre-defined comparison scenarios (tab 6)
## ---------------------------------------------------------------------------
SCENARIOS <- list(
  "1. 관찰 (watchful waiting)" =
    list(p = list(CAVFLAG = 0), e = NULL),
  "2. AZM+EMB+RIF 주3회 (결절기관지확장형)" =
    list(p = list(CAVFLAG = 0, MACTYPE = 0),
         e = list(tx = 365, macro = "AZM", freq = "TIW", emb = TRUE, rif = TRUE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = FALSE)),
  "3. AZM+EMB+RIF 매일 (공동형)" =
    list(p = list(CAVFLAG = 1, MACTYPE = 0),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = TRUE, rif = TRUE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = FALSE)),
  "4. + IV 아미카신 3개월 (공동형)" =
    list(p = list(CAVFLAG = 1),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = TRUE, rif = TRUE,
                  cfz = FALSE, iv = TRUE, wks = 13, alis = FALSE)),
  "5. + ALIS 흡입 590 mg/일 (CONVERT)" =
    list(p = list(CAVFLAG = 1),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = TRUE, rif = TRUE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = TRUE)),
  "6. 아지트로마이신 단독요법 (금기)" =
    list(p = list(CAVFLAG = 0),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = FALSE, rif = FALSE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = FALSE)),
  "7. CLR+EMB+RIF (CYP3A 상호작용)" =
    list(p = list(CAVFLAG = 1, MACTYPE = 1),
         e = list(tx = 365, macro = "CLR", freq = "BID", emb = TRUE, rif = TRUE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = FALSE)),
  "8. AZM+EMB+CFZ+ALIS (리팜핀 배제)" =
    list(p = list(CAVFLAG = 1),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = TRUE, rif = FALSE,
                  cfz = TRUE, iv = FALSE, wks = 0, alis = TRUE)),
  "9. 8번 + 기도청결요법 + 영양지원" =
    list(p = list(CAVFLAG = 1, ACT = 1, NUTR = 1),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = TRUE, rif = FALSE,
                  cfz = TRUE, iv = FALSE, wks = 0, alis = TRUE)),
  "10. 항IFN-γ 자가항체 숙주 + ALIS" =
    list(p = list(CAVFLAG = 1, IFNCAP = 0.15),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = TRUE, rif = TRUE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = TRUE)),
  "11. 공동형 AZM 단독요법 (숙주 정상)" =
    list(p = list(CAVFLAG = 1),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = FALSE, rif = FALSE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = FALSE)),
  "12. 공동형 단독요법 + 항IFN-γ (내성 sweep)" =
    list(p = list(CAVFLAG = 1, IFNCAP = 0.15),
         e = list(tx = 365, macro = "AZM", freq = "QD", emb = FALSE, rif = FALSE,
                  cfz = FALSE, iv = FALSE, wks = 0, alis = FALSE))
)

scenario_events <- function(s, wt = 52) {
  if (is.null(s$e)) return(NULL)
  with(s$e, build_events(tx, macro, freq, emb, rif, cfz, iv, wks, alis, wt))
}

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: -apple-system, 'Segoe UI', 'Noto Sans KR', sans-serif; }
    .banner { background:linear-gradient(90deg,#0277BD,#5E35B1);
              color:white; padding:14px 18px; border-radius:8px;
              margin-bottom:14px; }
    .banner h3 { margin:0 0 4px 0; font-weight:700; }
    .banner p  { margin:0; font-size:12.5px; opacity:.93; }
    .keybox { background:#FFF8E1; border-left:5px solid #F9A825;
              padding:10px 14px; border-radius:5px; margin:10px 0;
              font-size:13px; }
    .well { background:#F7F9FB; }
  "))),

  div(class = "banner",
      h3("비결핵 항산균 폐질환 (MAC-PD) QSP 모델 대시보드"),
      p("Nontuberculous Mycobacterial Lung Disease · Mycobacterium avium complex — ",
        "47-ODE niche-resolved PK/PD model  |  ",
        "핵심 구조: 세균은 4개의 물리적 니치에 나뉘어 있고, 식포 pH 5.2 하나가 ",
        "마크로라이드의 축적(151배)과 활성소실(12.6배)을 동시에 만든다.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("환자 표현형"),
      radioButtons("cav", "질환 표현형",
                   c("결절기관지확장형 (nodular-bronchiectatic)" = 0,
                     "섬유공동형 (fibrocavitary)" = 1),
                   selected = 1),
      sliderInput("wt", "체중 (kg)", 38, 80, 52, 1),
      sliderInput("ifncap", "Th1/IFN-γ 축 기능 (1 = 정상, ↓ = 항IFN-γ 항체·MSMD)",
                  0.05, 1, 1, 0.05),
      sliderInput("mcc0", "기저 점액섬모청소능 MCC₀", 0.20, 0.95, 0.75, 0.05),

      hr(), h4("약물 요법"),
      selectInput("macro", "마크로라이드",
                  c("아지트로마이신 (AZM)" = "AZM",
                    "클래리트로마이신 (CLR)" = "CLR",
                    "없음" = "none"), selected = "AZM"),
      selectInput("freq", "투여 빈도",
                  c("매일 (QD)" = "QD", "주 3회 (TIW)" = "TIW",
                    "1일 2회 (BID, CLR)" = "BID"), selected = "QD"),
      checkboxInput("emb", "에탐부톨 (EMB)", TRUE),
      checkboxInput("rif", "리팜핀 (RIF)", TRUE),
      checkboxInput("cfz", "클로파지민 (CFZ) 100 mg/일", FALSE),
      checkboxInput("alis", "ALIS 흡입 590 mg/일 (리포솜 아미카신)", FALSE),
      checkboxInput("ivamk", "IV 아미카신 15 mg/kg 주3회", FALSE),
      conditionalPanel("input.ivamk == true",
                       sliderInput("ivwks", "IV 아미카신 투여기간 (주)", 2, 26, 12, 1)),
      sliderInput("tx", "총 치료기간 (일)", 90, 540, 365, 30),

      hr(), h4("비약물 중재"),
      checkboxInput("act", "기도청결요법 (HFCWO + 고장성 식염수)", FALSE),
      checkboxInput("nutr", "영양 지원", FALSE),

      hr(),
      sliderInput("phag", "식포 pH (구조적 핵심 파라미터)", 4.5, 7.4, 5.2, 0.1),
      helpText("이 슬라이더 하나가 마크로라이드의 세포내 축적비와 MIC 상승을 ",
               "동시에 바꿉니다. 9번 탭에서 그 결과를 직접 확인하세요."),

      hr(),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---- TAB 1 ---------------------------------------------------------
        tabPanel("1. 환자 프로파일",
          br(),
          div(class = "keybox",
              strong("숙주 표현형이 왜 중요한가 — "),
              "MAC-PD는 병원체의 독력보다 숙주의 기도청소능과 체형에 의해 규정됩니다. ",
              "마른 체형·흉곽 이형성·CFTR 이형접합·섬모 기능 변이가 MCC를 낮추고, ",
              "낮아진 MCC는 세균 정체 → 바이오필름 → 기도 손상 → 다시 MCC 저하로 ",
              "닫히는 Cole 악순환의 입구가 됩니다."),
          fluidRow(
            column(6, plotOutput("p_profile", height = 320)),
            column(6, plotOutput("p_mcc", height = 320))),
          hr(),
          h4("설정 요약"), DTOutput("t_profile")),

        ## ---- TAB 2 ---------------------------------------------------------
        tabPanel("2. 약물 PK",
          br(),
          div(class = "keybox",
              strong("한 약물, 네 개의 농도. "),
              "혈장 농도는 이 질환에서 거의 아무것도 예측하지 못합니다. ",
              "아래 그림에서 같은 약물의 혈장 / ELF / 세포내 농도가 자릿수 단위로 ",
              "다르다는 점, 그리고 IV 아미카신의 세포내 농도가 정확히 0이라는 점을 ",
              "확인하세요."),
          fluidRow(
            column(6, plotOutput("p_pk_macro", height = 300)),
            column(6, plotOutput("p_pk_amk",  height = 300))),
          fluidRow(
            column(6, plotOutput("p_pk_other", height = 300)),
            column(6, plotOutput("p_enz",      height = 300))),
          hr(), h4("정상상태 노출 요약 (day 180–365)"), DTOutput("t_pk")),

        ## ---- TAB 3 ---------------------------------------------------------
        tabPanel("3. 세균 니치",
          br(),
          div(class = "keybox",
              strong("공동형 대 결절기관지확장형은 다른 모델이 아니라 같은 모델의 ",
                     "다른 초기값입니다. "),
              "CAVFLAG는 오직 B_C(건락 구획)의 초기값 하나만 바꿉니다. ",
              "경구요법이 B_C에 도달하지 못하기 때문에(마크로라이드 침투 0.15, ",
              "아미카신 0.02), 두 표현형의 배양음전율 격차가 산술적으로 발생합니다."),
          plotOutput("p_niche", height = 380),
          fluidRow(
            column(6, plotOutput("p_sputum", height = 300)),
            column(6, plotOutput("p_nichebar", height = 300)))),

        ## ---- TAB 4 ---------------------------------------------------------
        tabPanel("4. 임상 엔드포인트",
          br(),
          fluidRow(
            column(6, plotOutput("p_conv",  height = 300)),
            column(6, plotOutput("p_qolb",  height = 300))),
          fluidRow(
            column(6, plotOutput("p_cav",   height = 300)),
            column(6, plotOutput("p_bmi",   height = 300))),
          hr(), h4("주요 엔드포인트"), DTOutput("t_end")),

        ## ---- TAB 5 ---------------------------------------------------------
        tabPanel("5. 내성 진화",
          br(),
          div(class = "keybox",
              strong("내성은 스위치가 아니라 나눗셈입니다. "),
              "Φ = (마크로라이드 살균속도) / (총 살균속도). 어떤 니치에서 Φ가 1에 ",
              "가까워지면 그 니치에서 요법은 '기능적 마크로라이드 단독요법'이며, ",
              "MAC은 rrl 유전자가 단일 사본이므로 한 번의 점돌연변이로 완전 내성이 ",
              "됩니다. 식포 안에서는 동반약물이 도달하지 않으므로 Φ_I가 구조적으로 ",
              "높습니다 — ALIS를 추가하면 그 분모가 처음으로 채워집니다. ",
              br(), br(),
              "다만 검증 결과 Φ가 1에 가까워지는 것만으로는 내성이 자리잡지 않습니다. ",
              "내성균이 실제로 우세해지려면 ", strong("두 조건이 동시에"), " 필요합니다 — ",
              "① 선택압(Φ→1)과 ② 내성균의 순증식률이 양수가 되는 숙주 조건",
              "(높은 균량 + 손상된 IFN-γ/점액섬모 청소). 12번 시나리오에서만 ",
              "내성이 100%까지 sweep 합니다."),
          fluidRow(
            column(6, plotOutput("p_res",  height = 320)),
            column(6, plotOutput("p_phi",  height = 320))),
          plotOutput("p_respop", height = 300)),

        ## ---- TAB 6 ---------------------------------------------------------
        tabPanel("6. 시나리오 비교",
          br(),
          checkboxGroupInput("scn", "비교할 요법 선택",
                             choices = names(SCENARIOS),
                             selected = names(SCENARIOS)[c(1, 2, 3, 5, 6)],
                             inline = FALSE),
          actionButton("go_scn", "시나리오 비교 실행", class = "btn-primary"),
          br(), br(),
          plotOutput("p_scn_sput", height = 340),
          fluidRow(
            column(6, plotOutput("p_scn_res", height = 300)),
            column(6, plotOutput("p_scn_oto", height = 300))),
          hr(), h4("시나리오 요약표"), DTOutput("t_scn")),

        ## ---- TAB 7 ---------------------------------------------------------
        tabPanel("7. 바이오마커",
          br(),
          fluidRow(
            column(6, plotOutput("p_immune", height = 300)),
            column(6, plotOutput("p_protease", height = 300))),
          fluidRow(
            column(6, plotOutput("p_airway", height = 300)),
            column(6, plotOutput("p_symp",   height = 300)))),

        ## ---- TAB 8 ---------------------------------------------------------
        tabPanel("8. 독성 · 안전성",
          br(),
          div(class = "keybox",
              strong("이독성은 폐가 아니라 혈장을 따라갑니다. "),
              "외림프 아미카신은 KCEN(혈장)에서 유입되고, 효능은 KELF/KMAC(폐)에서 ",
              "나옵니다. 두 경로가 분리되어 있기 때문에 ALIS는 '폐는 높게, 귀는 낮게'가 ",
              "가정이 아니라 결과가 됩니다. 아래 lung:ear 비를 확인하세요."),
          fluidRow(
            column(6, plotOutput("p_oto",  height = 300)),
            column(6, plotOutput("p_lungear", height = 300))),
          fluidRow(
            column(6, plotOutput("p_tox_other", height = 300)),
            column(6, plotOutput("p_qtc", height = 300))),
          hr(), h4("최대 독성 지표"), DTOutput("t_tox")),

        ## ---- TAB 9 ---------------------------------------------------------
        tabPanel("9. pH 역설 탐색기",
          br(),
          div(class = "keybox",
              strong("입력은 하나, 결과는 둘. "),
              "식포 pH를 바꾸면 (i) Henderson-Hasselbalch 이온 포획비 R_trap 과 ",
              "(ii) 마크로라이드 MIC 상승배수가 동시에 움직입니다. 두 곡선이 교차하는 ",
              "지점 아래에서는 '세포 내로 더 많이 들어갈수록 덜 듣는' 구간이 됩니다."),
          fluidRow(
            column(7, plotOutput("p_paradox", height = 380)),
            column(5, br(), DTOutput("t_paradox"))),
          hr(),
          sliderInput("pka", "마크로라이드 pKa", 7.5, 9.8, 8.7, 0.1, width = "60%"),
          sliderInput("gm2", "pH 단위당 MIC 상승 (log10, 마크로라이드 GM)",
                      0.1, 1.2, 0.5, 0.05, width = "60%"),
          plotOutput("p_paradox_net", height = 300)),

        ## ---- TAB 10 --------------------------------------------------------
        tabPanel("10. 민감도 분석",
          br(),
          helpText("현재 설정된 요법에 대해 각 파라미터를 ±30% 변화시켰을 때 ",
                   "12개월 시점 객담 부하(log10 CFU/mL)가 얼마나 움직이는지 봅니다."),
          actionButton("go_sens", "민감도 분석 실행 (수 초 소요)",
                       class = "btn-warning"),
          br(), br(),
          plotOutput("p_sens", height = 460),
          DTOutput("t_sens")),

        ## ---- TAB 11 --------------------------------------------------------
        tabPanel("11. 모델 정보",
          br(),
          h4("구조"),
          tags$ul(
            tags$li(strong("46개 ODE"), " — 6개 약물 PK(경구 마크로라이드·EMB·RIF·CFZ, ",
                    "IV 아미카신, 흡입 리포솜 아미카신) + CYP3A 유도 + 6개 세균 구획 ",
                    "+ 숙주 면역/기도/구조 + 5개 독성 구획"),
            tags$li(strong("4개 물리적 니치"), " — B_E(기도내강 유리균) · B_B(바이오필름) · ",
                    "B_I(대식세포 식포) · B_C(건락/공동벽). 각 약물은 니치별로 서로 다른 ",
                    "유효농도를 가진다."),
            tags$li(strong("3개 내성 구획"), " — rrl A2058G 변이체가 세포외/바이오필름/",
                    "세포내에서 각각 복제 비례로 생성되고, 니치별 살균비에 의해 ",
                    "선택된다. 초기값은 MU × 해당 구획 균량이므로 균량이 낮으면 ",
                    "변이체가 애초에 존재하지 않는다."),
            tags$li(strong("12개 시나리오"), " — 관찰, 지침 3제(주3회/매일), IV 아미카신, ",
                    "ALIS, 단독요법, CYP3A 상호작용, 리팜핀 배제, 비약물 중재, ",
                    "면역저하 숙주, 내성 sweep")),
          h4("유도되는 값 (파라미터가 아님)"),
          tags$ul(
            tags$li("R_trap = (1+10^(pKa−pH_in))/(1+10^(pKa−pH_out)) ≈ 151× — ",
                    "아지트로마이신의 대식세포 축적"),
            tags$li("MIC 상승 = 10^(GM×(7.4−pH_in)) ≈ 12.6× — 같은 pH에서의 활성 소실"),
            tags$li("순 효력 이득 = 151 / 12.6 ≈ 12× — '100–1000배 폐 이행'의 실제 가치"),
            tags$li("Φ_I = 식포 내 마크로라이드 살균 분율 — 내성 선택 압력의 크기")),
          h4("한계"),
          tags$ul(
            tags$li("결정론적 단일 환자 모델입니다. 개체간 변이는 포함되지 않으므로 ",
                    "배양음전은 확률이 아니라 역치 통과로 표현됩니다. 확률적 소멸은 ",
                    "1 CFU 미만에서 증식을 멈추는 소멸 바닥(EXFLOOR)으로 근사했습니다."),
            tags$li("재감염과 재발을 구분하지 않습니다(유전형 분석이 필요한 구분)."),
            tags$li("파라미터는 문헌 기반 근사치이며 환자 데이터에 적합되지 ",
                    "않았습니다. 임상 의사결정에 사용할 수 없습니다.")),
          h4("참고문헌"),
          p("전체 목록은 ", code("ntm_references.md"), " 참조 (76개 항목, 12개 섹션)."),
          hr(),
          p(em("Claude Code Routine · QSP Disease Model Library · 2026-08-01")))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  cur_params <- reactive({
    list(CAVFLAG = as.numeric(input$cav),
         MACTYPE = ifelse(input$macro == "CLR", 1, 0),
         WTBL    = input$wt,
         IFNCAP  = input$ifncap,
         MCC0    = input$mcc0,
         PHPHAG  = input$phag,
         ACT     = as.numeric(input$act),
         NUTR    = as.numeric(input$nutr))
  })

  sim <- eventReactive(input$go, {
    ev_obj <- build_events(input$tx, input$macro, input$freq, input$emb,
                           input$rif, input$cfz, input$ivamk,
                           if (is.null(input$ivwks)) 0 else input$ivwks,
                           input$alis, input$wt)
    run_sim(cur_params(), ev_obj)
  }, ignoreNULL = FALSE)

  ## ---- helpers -------------------------------------------------------------
  lp <- function(d, ycol, ylab, title, sub = NULL, logy = FALSE, col = PAL[1]) {
    g <- ggplot(d, aes(time, .data[[ycol]])) +
      geom_line(linewidth = 0.9, colour = col) +
      labs(x = "시간 (일)", y = ylab, title = title, subtitle = sub) + THEME
    if (logy) g <- g + scale_y_log10()
    g
  }

  ## ================= TAB 1 =================================================
  output$p_profile <- renderPlot({
    d <- sim()
    d %>% select(time, WT, BMI) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = PAL) +
      labs(x = "시간 (일)", y = NULL,
           title = "체중 / BMI — TNF 구동 소모 고리",
           subtitle = "저체중은 위험인자이자 결과입니다 (양방향 화살표)") +
      THEME + theme(legend.position = "none")
  })

  output$p_mcc <- renderPlot({
    d <- sim()
    d %>% select(time, MCC, MUC, BRO) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = PAL[c(3, 4, 2)]) +
      labs(x = "시간 (일)", y = NULL,
           title = "기도 방어와 구조 손상 (Cole 악순환)",
           subtitle = "MCC ↓ → 점액 ↑ → 기관지확장증 ↑ → 다시 MCC ↓") +
      THEME + theme(legend.position = "none")
  })

  output$t_profile <- renderDT({
    p <- cur_params()
    tibble(항목 = c("질환 표현형", "체중", "Th1/IFN-γ 기능", "기저 MCC",
                    "식포 pH", "마크로라이드", "투여빈도", "EMB", "RIF",
                    "CFZ", "ALIS", "IV 아미카신", "기도청결", "영양지원",
                    "치료기간"),
           값 = c(ifelse(input$cav == 1, "섬유공동형", "결절기관지확장형"),
                  paste(input$wt, "kg"), input$ifncap, input$mcc0, input$phag,
                  input$macro, input$freq,
                  ifelse(input$emb, "예", "아니오"),
                  ifelse(input$rif, "예", "아니오"),
                  ifelse(input$cfz, "예", "아니오"),
                  ifelse(input$alis, "예", "아니오"),
                  ifelse(input$ivamk, paste0("예 (", input$ivwks, "주)"), "아니오"),
                  ifelse(input$act, "예", "아니오"),
                  ifelse(input$nutr, "예", "아니오"),
                  paste(input$tx, "일"))) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  ## ================= TAB 2 =================================================
  output$p_pk_macro <- renderPlot({
    sim() %>% select(time, `혈장` = CMPo, `ELF` = CMEo, `세포내(대식세포)` = CMIo) %>%
      pivot_longer(-time) %>% filter(value > 0) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) + scale_y_log10() +
      scale_colour_manual(values = PAL[c(1, 3, 5)]) +
      labs(x = "시간 (일)", y = "농도 (mg/L, log)",
           title = "마크로라이드: 같은 약물, 세 자릿수 차이",
           subtitle = "세포내 농도는 pKa와 식포 pH로부터 유도됩니다", colour = NULL) +
      THEME
  })

  output$p_pk_amk <- renderPlot({
    sim() %>% select(time, `혈장` = CKPo, `ELF(유리)` = CKEo,
                     `세포내(리포솜 경유)` = CKIo, `외림프` = KPERI) %>%
      pivot_longer(-time) %>%
      mutate(value = pmax(value, 1e-4)) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) + scale_y_log10() +
      scale_colour_manual(values = PAL[c(1, 3, 5, 2)]) +
      labs(x = "시간 (일)", y = "농도 (mg/L, log)",
           title = "아미카신: 같은 분자, 다른 주소",
           subtitle = "IV만 투여하면 '세포내' 곡선은 바닥에 붙어 있습니다",
           colour = NULL) + THEME
  })

  output$p_pk_other <- renderPlot({
    sim() %>% select(time, `에탐부톨` = CEPo, `리팜핀` = CFPo,
                     `클로파지민` = CCPo) %>%
      pivot_longer(-time) %>% mutate(value = pmax(value, 1e-4)) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) + scale_y_log10() +
      scale_colour_manual(values = PAL[c(3, 4, 8)]) +
      labs(x = "시간 (일)", y = "혈장 농도 (mg/L, log)",
           title = "동반 약물 혈장 농도", colour = NULL) + THEME
  })

  output$p_enz <- renderPlot({
    lp(sim(), "ENZ", "상대 CYP3A4 효소량 (기저 = 1)",
       "리팜핀에 의한 CYP3A 유도",
       "클래리트로마이신을 선택하면 이 곡선이 그 약의 노출을 그대로 깎습니다",
       col = PAL[4])
  })

  output$t_pk <- renderDT({
    d <- sim() %>% filter(time >= 180, time <= 365)
    tibble(
      지표 = c("마크로라이드 혈장 (mg/L)", "마크로라이드 ELF (mg/L)",
               "마크로라이드 세포내 (mg/L)", "세포내:혈장 비",
               "아미카신 ELF (mg/L)", "아미카신 세포내 (mg/L)",
               "아미카신 혈장 (mg/L)", "외림프 아미카신 (mg/L)",
               "세포외 C/MIC", "세포내 C/MIC", "CYP3A 유도배수"),
      평균 = round(c(mean(d$CMPo), mean(d$CMEo), mean(d$CMIo),
                     mean(d$CMIo) / max(mean(d$CMPo), 1e-9),
                     mean(d$CKEo), mean(d$CKIo), mean(d$CKPo), mean(d$KPERI),
                     mean(d$CME_MIC), mean(d$CMI_MIC), mean(d$ENZ)), 3)) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 15))
  })

  ## ================= TAB 3 =================================================
  output$p_niche <- renderPlot({
    sim() %>% select(time, `B_E 기도내강` = BE, `B_B 바이오필름` = BB,
                     `B_I 세포내` = BI, `B_C 건락/공동` = BC) %>%
      pivot_longer(-time, names_to = "niche", values_to = "cfu") %>%
      mutate(cfu = pmax(cfu, 1)) %>%
      ggplot(aes(time, log10(cfu), colour = niche)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c("#FFB74D", "#FF8A65", "#FF7043", "#D84315")) +
      labs(x = "시간 (일)", y = "log10 CFU", colour = NULL,
           title = "세균은 어디에 있는가 — 4개 니치의 시간경과",
           subtitle = "경구요법이 끝난 뒤 다시 올라오는 곡선이 재발의 저장고입니다") +
      THEME
  })

  output$p_sputum <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, LOGSPUT)) +
      geom_line(linewidth = 1, colour = PAL[1]) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = PAL[2]) +
      annotate("text", x = max(d$time) * 0.7, y = 1.35,
               label = "배양 음성 역치", size = 3.4, colour = PAL[2]) +
      labs(x = "시간 (일)", y = "객담 log10 CFU/mL",
           title = "객담 균량과 배양음전") + THEME
  })

  output$p_nichebar <- renderPlot({
    d <- sim()
    tp <- c(0, 90, 180, 365, max(d$time))
    d %>% filter(time %in% tp) %>%
      select(time, BE, BB, BI, BC) %>%
      pivot_longer(-time, names_to = "niche", values_to = "cfu") %>%
      mutate(cfu = pmax(cfu, 1), time = factor(time)) %>%
      ggplot(aes(time, log10(cfu), fill = niche)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c("#FFB74D", "#FF8A65", "#FF7043", "#D84315")) +
      labs(x = "시점 (일)", y = "log10 CFU", fill = NULL,
           title = "시점별 니치 구성") + THEME
  })

  ## ================= TAB 4 =================================================
  output$p_conv <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, TNEG)) +
      geom_line(linewidth = 1, colour = PAL[3]) +
      labs(x = "시간 (일)", y = "누적 배양음성 일수",
           title = "배양음전 유지 기간",
           subtitle = "지침 목표: 음전 후 12개월(365일) 유지") + THEME
  })
  output$p_qolb <- renderPlot(lp(sim(), "QOLB", "QOL-B 호흡영역 (0–100)",
                                 "삶의 질", "MCID ≈ 10점", col = PAL[5]))
  output$p_cav  <- renderPlot(lp(sim(), "CAV", "공동 부피 (cm³)",
                                 "공동 크기", col = PAL[2]))
  output$p_bmi  <- renderPlot(lp(sim(), "BMI", "BMI (kg/m²)",
                                 "체질량지수", col = PAL[6]))

  output$t_end <- renderDT({
    d <- sim()
    at <- function(t, col) d[[col]][which.min(abs(d$time - t))]
    tibble(
      엔드포인트 = c("객담 log10 CFU/mL (기저)", "객담 (6개월)", "객담 (12개월)",
                     "객담 (18개월, 치료 종료 후)",
                     "6개월 내 배양음전", "12개월 내 배양음전",
                     "누적 음성 일수", "내성 분율 (12개월)",
                     "공동 부피 (12개월, cm³)", "QOL-B (12개월)", "BMI (12개월)"),
      값 = c(round(at(0, "LOGSPUT"), 2), round(at(180, "LOGSPUT"), 2),
             round(at(365, "LOGSPUT"), 2), round(at(540, "LOGSPUT"), 2),
             ifelse(any(d$CULTNEG[d$time <= 180] == 1), "예", "아니오"),
             ifelse(any(d$CULTNEG[d$time <= 365] == 1), "예", "아니오"),
             round(max(d$TNEG), 0), signif(at(365, "RESFRAC"), 3),
             round(at(365, "CAV"), 1), round(at(365, "QOLB"), 1),
             round(at(365, "BMI"), 2))) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 15))
  })

  ## ================= TAB 5 =================================================
  output$p_res <- renderPlot({
    sim() %>% mutate(RESFRAC = pmax(RESFRAC, 1e-12)) %>%
      ggplot(aes(time, RESFRAC)) +
      geom_line(linewidth = 1, colour = PAL[2]) + scale_y_log10() +
      labs(x = "시간 (일)", y = "내성 분율 (log)",
           title = "rrl A2058G 내성균 분율",
           subtitle = "MAC은 rrl 단일 사본 — 한 번의 점돌연변이가 완전 내성") + THEME
  })

  output$p_phi <- renderPlot({
    sim() %>% select(time, `Φ_E 세포외` = PHI_E, `Φ_I 세포내` = PHI_I) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey40") +
      scale_colour_manual(values = PAL[c(1, 2)]) + ylim(0, 1) +
      labs(x = "시간 (일)", y = "Φ = 마크로라이드 살균 분율", colour = NULL,
           title = "니치별 내성 선택 게이트",
           subtitle = "점선(0.9) 위 = 사실상 마크로라이드 단독요법 구간") + THEME
  })

  output$p_respop <- renderPlot({
    sim() %>% select(time, `감수성 세포외` = BE, `내성 세포외` = RE,
                     `감수성 바이오필름` = BB, `내성 바이오필름` = RB,
                     `감수성 세포내` = BI, `내성 세포내` = RI) %>%
      pivot_longer(-time) %>% mutate(value = pmax(value, 1e-3)) %>%
      ggplot(aes(time, log10(value), colour = name)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(1, 2, 3, 4, 6, 7)]) +
      labs(x = "시간 (일)", y = "log10 CFU", colour = NULL,
           title = "감수성 대 내성 집단 (세포외 / 바이오필름 / 세포내)") + THEME
  })

  ## ================= TAB 6 =================================================
  scn_sim <- eventReactive(input$go_scn, {
    req(length(input$scn) > 0)
    withProgress(message = "시나리오 시뮬레이션 중...", {
      lapply(input$scn, function(nm) {
        s <- SCENARIOS[[nm]]
        p <- modifyList(list(WTBL = input$wt, PHPHAG = input$phag), s$p)
        incProgress(1 / length(input$scn), detail = nm)
        run_sim(p, scenario_events(s, input$wt)) %>% mutate(scenario = nm)
      }) %>% bind_rows()
    })
  })

  output$p_scn_sput <- renderPlot({
    ggplot(scn_sim(), aes(time, LOGSPUT, colour = scenario)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      labs(x = "시간 (일)", y = "객담 log10 CFU/mL", colour = NULL,
           title = "요법별 객담 균량") +
      THEME + theme(legend.text = element_text(size = 8))
  })

  output$p_scn_res <- renderPlot({
    scn_sim() %>% mutate(RESFRAC = pmax(RESFRAC, 1e-12)) %>%
      ggplot(aes(time, RESFRAC, colour = scenario)) +
      geom_line(linewidth = 0.9) + scale_y_log10() +
      labs(x = "시간 (일)", y = "내성 분율", colour = NULL,
           title = "요법별 내성 발현") +
      THEME + theme(legend.position = "none")
  })

  output$p_scn_oto <- renderPlot({
    ggplot(scn_sim(), aes(time, HEARDB, colour = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(x = "시간 (일)", y = "청력역치 상승 (dB)", colour = NULL,
           title = "요법별 이독성") +
      THEME + theme(legend.position = "none")
  })

  output$t_scn <- renderDT({
    scn_sim() %>% group_by(scenario) %>%
      summarise(`객담 6개월` = round(LOGSPUT[which.min(abs(time - 180))], 2),
                `객담 12개월` = round(LOGSPUT[which.min(abs(time - 365))], 2),
                `객담 18개월` = round(LOGSPUT[which.min(abs(time - 540))], 2),
                `6개월 음전` = ifelse(any(CULTNEG[time <= 180] == 1), "○", "×"),
                `12개월 음전` = ifelse(any(CULTNEG[time <= 365] == 1), "○", "×"),
                `음성 일수` = round(max(TNEG)),
                `내성분율(12M)` = signif(RESFRAC[which.min(abs(time - 365))], 3),
                `공동(12M)` = round(CAV[which.min(abs(time - 365))], 1),
                `QOL-B(12M)` = round(QOLB[which.min(abs(time - 365))], 1),
                `청력(dB)` = round(max(HEARDB), 1),
                `QTc(ms)` = round(max(QTCMS), 1),
                .groups = "drop") %>%
      datatable(rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })

  ## ================= TAB 7 =================================================
  output$p_immune <- renderPlot({
    sim() %>% select(time, `대식세포` = MPH, `IFN-γ` = IFNG, `TNF-α` = TNF) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(1, 3, 2)]) +
      labs(x = "시간 (일)", y = "상대 단위", colour = NULL,
           title = "면역 축") + THEME
  })
  output$p_protease <- renderPlot({
    sim() %>% select(time, `호중구` = NEU, `MMP-1/9` = MMP) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(4, 7)]) +
      labs(x = "시간 (일)", y = "상대 단위", colour = NULL,
           title = "호중구–단백분해효소 축",
           subtitle = "마크로라이드의 항염 작용이 균이 죽기 전에 여기를 먼저 낮춥니다") +
      THEME
  })
  output$p_airway <- renderPlot({
    sim() %>% select(time, `점액 부하` = MUC, `MCC` = MCC,
                     `기관지확장증 점수` = BRO) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = PAL[c(6, 3, 2)]) +
      labs(x = "시간 (일)", y = NULL, title = "기도 상태") +
      THEME + theme(legend.position = "none")
  })
  output$p_symp <- renderPlot(lp(sim(), "SYM", "증상 점수 (0–10)",
                                 "증상 (기침·객담·피로·발열)", col = PAL[7]))

  ## ================= TAB 8 =================================================
  output$p_oto <- renderPlot(lp(sim(), "HEARDB", "청력역치 상승 (dB)",
                                "이독성 — 누적 외림프 노출 구동", col = PAL[2]))
  output$p_lungear <- renderPlot({
    sim() %>% mutate(LUNGEAR = pmax(LUNGEAR, 1e-3)) %>%
      ggplot(aes(time, LUNGEAR)) +
      geom_line(linewidth = 1, colour = PAL[5]) + scale_y_log10() +
      labs(x = "시간 (일)", y = "ELF : 외림프 아미카신 비 (log)",
           title = "폐 : 귀 비율 — ALIS의 존재 이유",
           subtitle = "IV에서는 이 비가 낮고, ALIS에서는 자릿수 단위로 올라갑니다") +
      THEME
  })
  output$p_tox_other <- renderPlot({
    sim() %>% select(time, `시신경(EMB)` = OPT, `신세뇨관(아미카신)` = NEP,
                     `간(RIF)` = HEP) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(3, 5, 4)]) +
      labs(x = "시간 (일)", y = "손상 지표", colour = NULL,
           title = "기타 장기 독성") + THEME
  })
  output$p_qtc <- renderPlot(lp(sim(), "QTCMS", "QTc 연장 (ms)",
                                "심전도 — 마크로라이드 hERG 차단",
                                "클래리트로마이신 > 아지트로마이신", col = PAL[8]))

  output$t_tox <- renderDT({
    d <- sim()
    tibble(독성 = c("청력역치 상승 (dB)", "시신경 손상 지표 (0–1)",
                    "신세뇨관 손상 지표 (0–1)", "간손상 지표", "최대 QTc 연장 (ms)",
                    "최대 외림프 아미카신 (mg/L)", "최대 혈장 아미카신 (mg/L)"),
           최댓값 = round(c(max(d$HEARDB), max(d$OPT), max(d$NEP), max(d$HEP),
                            max(d$QTCMS), max(d$KPERI), max(d$CKPo)), 3)) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## ================= TAB 9 =================================================
  paradox_df <- reactive({
    pH <- seq(4.5, 7.4, by = 0.05)
    tibble(pH = pH,
           R_trap = (1 + 10^(input$pka - pH)) / (1 + 10^(input$pka - 7.4)),
           MIC_fold = 10^(input$gm2 * (7.4 - pH))) %>%
      mutate(net = R_trap / MIC_fold)
  })

  output$p_paradox <- renderPlot({
    paradox_df() %>%
      select(pH, `R_trap (축적)` = R_trap, `MIC 상승배수 (활성소실)` = MIC_fold) %>%
      pivot_longer(-pH) %>%
      ggplot(aes(pH, value, colour = name)) +
      geom_line(linewidth = 1.1) + scale_y_log10() +
      geom_vline(xintercept = input$phag, linetype = "dashed", colour = "grey30") +
      annotate("text", x = input$phag, y = 1.2, hjust = -0.1, size = 3.5,
               label = paste0("현재 식포 pH = ", input$phag)) +
      scale_colour_manual(values = c(PAL[1], PAL[2])) +
      labs(x = "식포 pH", y = "배수 (log)", colour = NULL,
           title = "하나의 입력, 두 개의 결과",
           subtitle = "pH를 내리면 더 많이 들어가고(파랑), 동시에 덜 듣습니다(빨강)") +
      THEME
  })

  output$p_paradox_net <- renderPlot({
    paradox_df() %>%
      ggplot(aes(pH, net)) +
      geom_line(linewidth = 1.2, colour = PAL[5]) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = PAL[2]) +
      geom_vline(xintercept = input$phag, linetype = "dashed", colour = "grey30") +
      scale_y_log10() +
      labs(x = "식포 pH", y = "순 효력 이득 = R_trap / MIC 상승배수 (log)",
           title = "실제로 남는 것",
           subtitle = "빨간 선(=1) 아래로 내려가면 세포내 이행이 오히려 손해입니다") +
      THEME
  })

  output$t_paradox <- renderDT({
    d <- paradox_df()
    row <- d[which.min(abs(d$pH - input$phag)), ]
    tibble(항목 = c("식포 pH", "마크로라이드 pKa", "이온포획 축적비 R_trap",
                    "MIC 상승배수", "순 효력 이득",
                    "세포외 대비 실질 이득"),
           값 = c(input$phag, input$pka, round(row$R_trap, 1),
                  round(row$MIC_fold, 2), round(row$net, 1),
                  paste0(round(row$net, 1), " 배"))) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## ================= TAB 10 ================================================
  sens_res <- eventReactive(input$go_sens, {
    ps <- c("PHPHAG", "GM", "GK", "PCSM", "PCSK", "PBFK", "PBFL", "KUPT",
            "TOLCS", "TOLBF", "MU", "IFNCAP", "MCC0", "EMAXPERM", "RMELF")
    ev_obj <- build_events(input$tx, input$macro, input$freq, input$emb,
                           input$rif, input$cfz, input$ivamk,
                           if (is.null(input$ivwks)) 0 else input$ivwks,
                           input$alis, input$wt)
    base_pars <- cur_params()
    base_out <- run_sim(base_pars, ev_obj, end = 365)
    base_val <- base_out$LOGSPUT[nrow(base_out)]
    withProgress(message = "민감도 분석 중...", {
      lapply(ps, function(p) {
        incProgress(1 / length(ps), detail = p)
        b <- as.numeric(param(mod)[[p]])
        if (p %in% names(base_pars)) b <- as.numeric(base_pars[[p]])
        vapply(c(0.7, 1.3), function(f) {
          pp <- base_pars; pp[[p]] <- b * f
          o <- run_sim(pp, ev_obj, end = 365)
          o$LOGSPUT[nrow(o)]
        }, numeric(1)) -> v
        tibble(parameter = p, low = v[1], high = v[2], base = base_val,
               swing = abs(v[2] - v[1]))
      }) %>% bind_rows() %>% arrange(desc(swing))
    })
  })

  output$p_sens <- renderPlot({
    d <- sens_res()
    d %>% mutate(parameter = factor(parameter, levels = rev(d$parameter))) %>%
      ggplot() +
      geom_segment(aes(y = parameter, yend = parameter, x = low, xend = high),
                   linewidth = 5, colour = PAL[1], alpha = 0.55) +
      geom_point(aes(y = parameter, x = base), colour = PAL[2], size = 2.6) +
      labs(x = "12개월 객담 log10 CFU/mL", y = NULL,
           title = "국소 민감도 — 파라미터 ±30%",
           subtitle = "빨간 점 = 기준값. 막대가 길수록 결과를 크게 움직이는 파라미터") +
      THEME
  })

  output$t_sens <- renderDT({
    sens_res() %>% mutate(across(where(is.numeric), ~round(.x, 3))) %>%
      datatable(rownames = FALSE, options = list(dom = "tp", pageLength = 15))
  })
}

shinyApp(ui, server)
