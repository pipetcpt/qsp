## =====================================================================
##  CVID QSP MODEL — SHINY DASHBOARD
##  보통 가변 면역결핍증 (Common Variable Immunodeficiency)
##  ---------------------------------------------------------------------
##  10 tabs, organised around the two arms of the disease:
##
##    ARM 1  antibody deficiency  — tabs 2, 3, 4  (replaceable)
##    ARM 2  immune dysregulation — tabs 6, 7     (not replaceable)
##    the layer neither arm reverses — tab 5      (irreversible)
##
##  The controls are deliberately arranged so that the two questions a
##  CVID clinician actually asks are one slider apart:
##      "how much immunoglobulin, given how?"        (ARM 1 panel)
##      "and what do I do about the dysregulation?"  (ARM 2 panel)
##
##  Run:  shiny::runApp("cvid_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## =====================================================================

suppressPackageStartupMessages({
  library(shiny); library(mrgsolve); library(dplyr)
  library(tidyr); library(ggplot2); library(DT)
})

## ---------------------------------------------------------------------
##  Load the model. The model file ends in an invisible list, and its
##  pre-runs solve the healthy and CVID baselines, so sourcing it gives
##  us the compiled model plus both sets of initial conditions.
## ---------------------------------------------------------------------
MODEL_FILE <- "cvid_mrgsolve_model.R"
message("Sourcing ", MODEL_FILE, " (compiles the model and solves baselines)...")
BUILD <- source(MODEL_FILE, local = new.env())$value
mod        <- BUILD$mod
CVID_SS    <- BUILD$cvid
HEALTHY_SS <- BUILD$healthy
CVID_GENO  <- list(HEALTHY = 0, DAMAGEON = 1, FCSR = 0.06, FPCD = 0.60,
                   FTACI = 0.35, FTFH = 0.55)

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10))

PAL <- c("#1565C0", "#C62828", "#2E7D32", "#6A1B9A", "#EF6C00", "#00838F")

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel(HTML(paste0(
    "<b>보통 가변 면역결핍증 (CVID) — QSP 시뮬레이터</b>",
    "<br><span style='font-size:13px;color:#555'>",
    "ARM 1 항체결핍(보충 가능) &nbsp;∥&nbsp; ARM 2 면역조절이상(보충 불가) ",
    "&nbsp;→&nbsp; 비가역적 구조 손상</span>"))),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("환자 (Patient)"),
      sliderInput("wt", "체중 Body weight (kg)", 40, 120, 70, 5),
      sliderInput("delay", "진단 지연 Diagnostic delay (년)", 0, 15, 2, 1),
      sliderInput("years", "시뮬레이션 기간 (년)", 5, 30, 20, 5),
      helpText(HTML("<small>지연 기간 동안 래칫은 <b>제동 없이</b> 돌아갑니다.
                     이 슬라이더가 모델에서 가장 강력한 변수입니다.</small>")),

      hr(), h4("ARM 1 · 면역글로불린 보충"),
      radioButtons("route", "경로 Route",
                   c("IVIG (q3-4주)" = "IV",
                     "SCIG (주 1-2회)" = "SC",
                     "fSCIG (히알루로니다제, q4주)" = "FSC"),
                   selected = "IV"),
      sliderInput("dose", "용량 Dose (mg/kg per 4 weeks, IV-equivalent)",
                  200, 1500, 500, 50),
      conditionalPanel("input.route == 'IV'",
        radioButtons("tau_iv", "투여 간격", c("q3주" = 21, "q4주" = 28), 28)),
      conditionalPanel("input.route == 'SC'",
        radioButtons("tau_sc", "투여 빈도",
                     c("주 1회" = 7, "주 2회" = 3.5, "매일" = 1), 7)),
      checkboxInput("scadj", "SC 용량 보정계수 1.37 적용 (AUC 매칭)", TRUE),

      hr(), h4("ARM 2 · 면역조절이상"),
      sliderInput("dysgeno", "조절이상 부하 Dysregulation load", 0, 1, 0, 0.05),
      selectInput("pheno", "우세 표현형 (Chapel 표현형)",
                  c("없음 (합병증 없는 CVID)" = "none",
                    "GLILD / 육아종" = "glild",
                    "자가면역 세포감소증 (ITP/AIHA)" = "cytopenia",
                    "장병증 / 단백소실" = "enteropathy",
                    "다클론성 림프증식" = "lproc")),
      selectInput("geno", "유전자형 Genotype",
                  c("미확인 (다인성)" = "poly",
                    "CTLA4 반접합부족증 (CHAI)" = "ctla4",
                    "LRBA 결핍" = "lrba",
                    "APDS (PIK3CD GOF)" = "apds")),

      hr(), h4("치료 Therapies"),
      checkboxGroupInput("tx", NULL,
        c("리툭시맙 Rituximab" = "rtx",
          "아자티오프린/MMF" = "aza",
          "프레드니손 Prednisone" = "gc",
          "아바타셉트 Abatacept" = "aba",
          "시롤리무스 Sirolimus" = "siro",
          "레니올리십 Leniolisib (APDS)" = "lenio",
          "JAK 억제제" = "jak",
          "엘트롬보팍 Eltrombopag" = "elt",
          "비장절제 Splenectomy" = "splenec",
          "아지트로마이신 예방" = "azm",
          "기도 청결 물리치료" = "physio")),
      sliderInput("tstart", "ARM 2 치료 시작 (년)", 0, 20, 5, 1),

      hr(),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary btn-block"),
      checkboxInput("cmp", "표준 요법(IVIG 500 q4w)과 비교", TRUE)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---- 1 -----------------------------------------------------
        tabPanel("1 · 환자 프로파일",
          br(), fluidRow(
            column(4, wellPanel(h5("면역 표현형 (Immune phenotype)"),
                                tableOutput("t_pheno"))),
            column(4, wellPanel(h5("20년 시점 결과 (Outcomes)"),
                                tableOutput("t_outcome"))),
            column(4, wellPanel(h5("모델 구조"),
                                tableOutput("t_struct")))),
          hr(), h5("EUROclass / Freiburg 위치"),
          plotOutput("p_euro", height = "300px"),
          hr(),
          HTML("<div style='background:#FFFDE7;border-left:4px solid #F9A825;
                padding:10px'><b>이 모델의 전제</b><br>
                CVID는 하나의 검사실 정의를 공유하는 두 개의 질환이며,
                우리가 투여하는 치료는 그중 <b>하나만</b> 치료합니다.
                <b>ARM 1</b>(항체결핍 → 감염)은 혈청 IgG 농도로 정량적으로
                예측되고 보충요법이 거의 해결합니다.
                <b>ARM 2</b>(면역조절이상 → 자가면역·GLILD·장병증·림프종)는
                T세포·관용 장애이며 보충요법이 사실상 아무 영향을 주지
                않습니다. 그리고 사람을 죽이는 것은 후자입니다
                (Resnick 2012: 비감염 합병증 사망 RR ≈ 11배).</div>")),

        ## ---- 2 -----------------------------------------------------
        tabPanel("2 · IgG PK",
          br(), plotOutput("p_pk", height = "330px"),
          hr(), fluidRow(
            column(6, plotOutput("p_pk_zoom", height = "300px")),
            column(6, br(), tableOutput("t_pk"),
                   HTML("<small><b>같은 월 용량, 다른 trough.</b>
                    trough는 월 투여량만으로 정해지지 않고 <i>어떻게 나누어
                    투여했는지</i>로 정해집니다. IVIG q4주는 600-900 mg/dL의
                    톱니파를 만들고, 같은 AUC의 주 1회 SCIG는 거의 평탄하므로
                    <b>같은 밀리그램으로 더 높은 trough</b>를 얻습니다.
                    </small>")))),

        ## ---- 3 -----------------------------------------------------
        tabPanel("3 · 노출-반응 (Exposure-Response)",
          br(), fluidRow(
            column(6, plotOutput("p_er", height = "340px")),
            column(6, plotOutput("p_er2", height = "340px"))),
          hr(),
          HTML("<div style='background:#E3F2FD;border-left:4px solid #1565C0;
                padding:10px'><b>Orange 2010의 기울기는 맞고, 그 따름정리는
                틀렸습니다.</b><br>
                17개 IVIG 연구의 메타회귀: 폐렴 발생률은 IgG trough 100 mg/dL
                상승마다 약 <b>27% 감소</b>합니다. 이 기울기는 견고합니다.
                함께 인용되는 '1400 mg/dL에서 폐렴 0'은 <i>비율에 직선을
                맞춘</i> 인공물입니다. 이 모델은 기울기는 유지하되 Hill 형태로
                구현하여, 방어가 <b>무균에 도달하지 않으면서</b> 개선되고
                800-1000 mg/dL 부근에 수확 체감 무릎이 생깁니다.
                직선형은 '1400까지 올리고 잊어라'라고 말하고, Hill형은
                '그 위에서는 사는 것이 거의 없다'고 말합니다.</div>"),
          hr(), h5("볼록성(Jensen) 벌점 — 프로파일 모양만으로 생기는 초과 위험"),
          plotOutput("p_jensen", height = "300px")),

        ## ---- 4 -----------------------------------------------------
        tabPanel("4 · 감염 (Infection)",
          br(), plotOutput("p_inf", height = "330px"),
          hr(), fluidRow(
            column(7, plotOutput("p_divergence", height = "300px")),
            column(5, br(), HTML("<b>폐렴과 부비동염의 용량-반응은 다릅니다.</b>
              <br><small>혈청 IgG는 기도 표면으로 잘 투과되지 않고,
              분비형 IgA는 <b>어떤 제품으로도 보충되지 않습니다</b>.
              따라서 점막 감염률에는 IgG와 <b>무관한 바닥</b>이 존재하며,
              보충 중인 환자도 연 2-3회의 부비동염·기관지염을 겪습니다.
              결과: 재발 <b>폐렴</b>에 대한 증량은 가치가 있고,
              지속 <b>부비동염</b>에 대한 증량은 대체로 무의미합니다.
              그 바닥에 작용하는 것은 아지트로마이신입니다.</small>")))),

        ## ---- 5 -----------------------------------------------------
        tabPanel("5 · 비가역 손상 (The Ratchet)",
          br(), plotOutput("p_ratchet", height = "340px"),
          hr(), fluidRow(
            column(6, plotOutput("p_delay", height = "320px")),
            column(6, plotOutput("p_lung", height = "320px"))),
          hr(),
          HTML("<div style='background:#FFEBEE;border-left:4px solid #B71C1C;
                padding:10px'><b>래칫 (The ratchet).</b><br>
                집락화 → IL-8 → 호중구 엘라스타제 → 기도벽 파괴 →
                <b>기관지확장증(비가역)</b> → 청소율 저하 → 점액 정체 →
                <i>더 많은 집락화</i>. 비가역 상태변수를 내부에 품은
                양성 되먹임 고리입니다.<br>
                동일한 정상상태 trough, 동일한 모든 것을 가진 두 환자는
                한 명이 1년째, 다른 한 명이 7년째 진단되었다는 이유만으로
                <b>영구적으로 다른 FEV1</b>을 갖습니다. 보충요법은 래칫을
                멈추지만 되돌리지는 못합니다.<br>
                따라서 CVID에서 가장 가치 있는 개입은 더 좋은 면역글로불린이
                아니라 <b>더 이른 혈청 IgG 측정</b>입니다.</div>")),

        ## ---- 6 -----------------------------------------------------
        tabPanel("6 · ARM 2 · 면역조절이상",
          br(), plotOutput("p_arm2", height = "340px"),
          hr(), fluidRow(
            column(6, plotOutput("p_cyto", height = "300px")),
            column(6, plotOutput("p_glild", height = "300px"))),
          hr(),
          HTML("<div style='background:#F3E5F5;border-left:4px solid #4A148C;
                padding:10px'><b>구조적 영(null).</b> GLILD 방정식에는
                IgG 항이 <b>직접적으로 존재하지 않습니다</b>. 모델은
                IVIG 300 → 1000 mg/kg 전 범위에서 GLILD 활성도를 8% 남짓
                움직이는데, 이조차 직접 효과가 아니라
                감염 → 만성 활성화 → CD21low B → 림프 집합체의 <b>간접</b>
                경로입니다. 같은 모델에서 리툭시맙+아자티오프린은 80%를
                움직입니다(약 10배 차이). 보충요법을 최대로 올려도 GLILD는
                거의 그대로라는 뜻이며, 이는 파라미터가 아니라
                <b>구조</b>에서 나옵니다.</div>")),

        ## ---- 7 -----------------------------------------------------
        tabPanel("7 · 면역억제의 대가",
          br(), plotOutput("p_isb", height = "330px"),
          hr(), fluidRow(
            column(6, plotOutput("p_tradeoff", height = "310px")),
            column(6, br(), tableOutput("t_tradeoff"),
              HTML("<small><b>CVID에서 리툭시맙은 상대적으로 저렴합니다.</b>
              다른 질환에서 리툭시맙의 주된 체액성 대가는 이차성
              감마글로불린저하증입니다. CVID에서는 그 대가가 <b>이미
              지불되었고 이미 보충되고 있습니다</b>. 모델에서 동일한
              B세포 고갈(-98%)이 정상 숙주에서는 IgG를 -11% 떨어뜨리지만
              보충 중인 CVID에서는 -3%만 떨어뜨립니다. 리툭시맙의 진짜
              대가는 체액성이 아니라 <b>세포성</b>입니다.<br><br>
              반대 방향의 대가도 있습니다. ARM 2 치료는 모두 순
              면역억제를 올리고, 이는 곧바로 ARM 1 감염 위험으로
              되돌아옵니다.</small>")))),

        ## ---- 8 -----------------------------------------------------
        tabPanel("8 · 시나리오 비교",
          br(), fluidRow(
            column(3, selectInput("cmp_set", "비교 세트",
              c("용량 (300-1000 mg/kg)" = "dose",
                "경로·간격 (AUC 매칭)" = "route",
                "진단 지연 (1/4/7/15년)" = "delay",
                "불응성 ITP 4개 전략" = "itp",
                "GLILD 4개 전략" = "glild"))),
            column(3, selectInput("cmp_y", "종점",
              c("혈청 IgG" = "IGG", "폐렴률" = "PNEU_YR",
                "부비동염률" = "SINO_YR", "기관지확장증(Reiff)" = "BE",
                "FEV1 %pred" = "FEV1", "DLCO %pred" = "DLCO",
                "GLILD 활성도" = "GLILD", "혈소판" = "PLT",
                "ARM 2 부담" = "ARM2", "삶의 질" = "QOL"))),
            column(6, br(), helpText("각 세트는 실행 시 새로 시뮬레이션됩니다."))),
          plotOutput("p_cmp", height = "380px"),
          hr(), DTOutput("t_cmp")),

        ## ---- 9 -----------------------------------------------------
        tabPanel("9 · 바이오마커",
          br(), plotOutput("p_bio", height = "360px"),
          hr(), fluidRow(
            column(6, plotOutput("p_bcell", height = "300px")),
            column(6, plotOutput("p_baff", height = "300px"))),
          hr(), HTML("<small><b>보충요법이 가려버리는 지표와 가리지 못하는
            지표.</b> 혈청 IgG는 투여 후 해석 불가가 됩니다. 그러나
            <b>IgA·IgM</b>, <b>전환 기억 B %</b>, <b>CD21low %</b>,
            <b>유리 경쇄</b>는 외인성 IgG의 영향을 받지 않으므로 잔존
            내인성 B세포 기능과 조절이상 궤적을 계속 읽을 수 있습니다.
            가용성 <b>BAFF</b>는 B세포 sink가 줄어들기 때문에 CVID에서
            상승해 있고, 리툭시맙 후에는 더 올라갑니다.</small>")),

        ## ---- 10 ----------------------------------------------------
        tabPanel("10 · 검증 & 진단",
          br(), h5("공개 문헌 앵커 대조 (Validation anchors)"),
          DTOutput("t_anchor"),
          hr(), h5("모델 진단 (Diagnostics — 34개, 모두 실패할 수 있는 검사)"),
          DTOutput("t_diag"),
          hr(), HTML("<small>진단 D7은 원래 <i>정확한 0</i>을 기대하고
            작성되었고, 모델이 그 기대를 <b>반증</b>했습니다. GLILD 방정식에
            IgG 항은 없지만 감염 → 만성 활성화 경로를 통한 8% 남짓의
            간접 민감도가 남습니다. 생물학적으로 실재하는 경로이므로
            제거하지 않고 <b>보고</b>합니다. D15(레니올리십)는 임의의
            림프절 축소 항을 제거하고 PI3Kδ 억제만으로 남긴 결과 관찰치
            -39%에 대해 -21%로 <b>과소예측</b>합니다. 이 역시 맞추지 않고
            그대로 보고합니다.</small>"))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  ## ---- assemble parameters from the UI ------------------------------
  pars_of <- function() {
    p <- CVID_GENO
    p$WT <- input$wt
    p$VC <- 0.055 * input$wt * 10      # dL
    p$VP <- 0.045 * input$wt * 10
    p$DYSGENO <- input$dysgeno

    ## Chapel phenotypes are largely mutually exclusive: emphasise one arm
    p$KINAA <- 0.006; p$KINLA <- 0.0020; p$KINGUT <- 0.0080
    if (input$pheno == "glild")      p$KINLA  <- 0.0030
    if (input$pheno == "cytopenia")  p$KINAA  <- 0.0300
    if (input$pheno == "enteropathy")p$KINGUT <- 0.0120
    if (input$pheno == "lproc")    { p$KINSPL <- 0.0045; p$KINLAD <- 0.0100 }

    if (input$geno == "ctla4") p$CTLA4G  <- 0.5
    if (input$geno == "lrba")  p$CTLA4G  <- 0.3
    if (input$geno == "apds")  p$PI3KGOF <- 0.8

    t0 <- input$tstart * 365
    if ("aza"     %in% input$tx) { p$AZAON <- 1; p$AZAT0 <- t0 }
    if ("jak"     %in% input$tx) { p$JAKON <- 1; p$JAKT0 <- t0 }
    if ("azm"     %in% input$tx) { p$AZMON <- 1; p$AZMT0 <- input$delay*365 }
    if ("physio"  %in% input$tx)   p$PHYSIOON <- 1
    if ("splenec" %in% input$tx) { p$SPLENEC <- 1; p$TSPLENEC <- t0 }
    if (input$route == "FSC")      p$HYAL <- 1
    p
  }

  ## ---- assemble the event table -------------------------------------
  ev_of <- function(dose = input$dose, route = input$route) {
    yrs <- input$years; wt <- input$wt; d0 <- input$delay * 365
    tau <- switch(route, IV = as.numeric(input$tau_iv),
                         SC = as.numeric(input$tau_sc), FSC = 28)
    ## dose is specified as IV-equivalent mg/kg per 4 weeks
    per_dose <- dose * tau / 28
    if (route == "SC"  && isTRUE(input$scadj)) per_dose <- per_dose / 0.73
    if (route == "FSC" && isTRUE(input$scadj)) per_dose <- per_dose / 0.93
    cmtn <- if (route == "IV") "IGG_C" else "IGG_SC"
    nn <- max(1, floor((yrs*365 - d0)/tau) + 1)
    e <- as.data.frame(ev(amt = per_dose*wt, cmt = cmtn, time = d0,
                          ii = tau, addl = nn - 1, evid = 1))
    t0 <- input$tstart * 365; add <- list()
    if ("rtx" %in% input$tx)
      add[[length(add)+1]] <- as.data.frame(ev(
        amt = 700, cmt = "RTX_C", time = t0, ii = 7, addl = 3))
      if ("rtx" %in% input$tx)
        add[[length(add)+1]] <- as.data.frame(ev(
          amt = 700, cmt = "RTX_C", time = t0 + 182, ii = 182,
          addl = max(0, floor((yrs*365 - t0)/182) - 1)))
    if ("gc" %in% input$tx)
      add[[length(add)+1]] <- as.data.frame(ev(
        amt = 20, cmt = "PRED_C", time = t0, ii = 1,
        addl = max(0, yrs*365 - t0)))
    if ("aba" %in% input$tx)
      add[[length(add)+1]] <- as.data.frame(ev(
        amt = 10*wt, cmt = "ABA_C", time = t0, ii = 28,
        addl = max(0, ceiling((yrs*365 - t0)/28))))
    if ("siro" %in% input$tx)
      add[[length(add)+1]] <- as.data.frame(ev(
        amt = 2000, cmt = "SIRO_C", time = t0, ii = 1,
        addl = max(0, yrs*365 - t0)))
    if ("lenio" %in% input$tx)
      add[[length(add)+1]] <- as.data.frame(ev(
        amt = 70000, cmt = "LENIO_C", time = t0, ii = 0.5,
        addl = max(0, (yrs*365 - t0)*2)))
    if ("elt" %in% input$tx)
      add[[length(add)+1]] <- as.data.frame(ev(
        amt = 50, cmt = "ELT_C", time = t0, ii = 1,
        addl = max(0, yrs*365 - t0)))
    out <- bind_rows(c(list(e), add)); out$ID <- 1
    out %>% arrange(time)
  }

  sim1 <- function(pars, events, yrs, delta = 7) {
    mod %>% param(pars) %>% init(CVID_SS) %>%
      mrgsim(data = events, end = yrs*365, delta = delta,
             atol = 1e-8, rtol = 1e-6) %>% as_tibble()
  }

  R <- eventReactive(input$go, ignoreNULL = FALSE, {
    withProgress(message = "시뮬레이션 중...", value = 0.3, {
      main <- sim1(pars_of(), ev_of(), input$years) %>% mutate(arm = "선택 요법")
      incProgress(0.4)
      ref <- if (isTRUE(input$cmp)) {
        p <- CVID_GENO; p$WT <- input$wt
        p$VC <- 0.055*input$wt*10; p$VP <- 0.045*input$wt*10
        p$DYSGENO <- input$dysgeno
        sim1(p, ev_of(500, "IV"), input$years) %>%
          mutate(arm = "표준 IVIG 500 q4w")
      } else NULL
      list(main = main, both = bind_rows(main, ref))
    })
  })

  yr <- function(d) d$time / 365

  ## ================= TAB 1 ==========================================
  output$t_pheno <- renderTable({
    d <- R()$main; l <- tail(d, 1)
    tibble(지표 = c("혈청 IgG (mg/dL)", "전환 기억 B (% of B)",
                    "CD21low B (% of B)", "나이브 B (% of B)",
                    "총 B세포 (/uL)", "가용성 BAFF (x정상)"),
           값 = sprintf("%.1f", c(l$IGG, l$SMBPCT, l$CD21PCT,
                                 l$NAIVEPCT, l$BTOTAL, l$BAFF)))
  }, colnames = TRUE)

  output$t_outcome <- renderTable({
    d <- R()$main; l <- tail(d, 1)
    w <- d %>% filter(time >= max(time) - 28)
    tibble(지표 = c("IgG trough (mg/dL)", "폐렴 (/년)", "부비동염 (/년)",
                    "기관지확장증 (Reiff 0-18)", "FEV1 (%pred)",
                    "DLCO (%pred)", "누적 폐렴", "삶의 질"),
           값 = sprintf("%.2f", c(min(w$IGG), mean(w$PNEU_YR),
                                 mean(w$SINO_YR), l$BE, l$FEV1, l$DLCO,
                                 l$PNEU_CUM, l$QOL)))
  })

  output$t_struct <- renderTable({
    tibble(항목 = c("ODE 구획", "파라미터", "시나리오", "진단 검사",
                    "지도 노드", "참고문헌"),
           값 = c("74", "224", "30", "34 (34 PASS)", "234", "118"))
  })

  output$p_euro <- renderPlot({
    l <- tail(R()$main, 1)
    grid <- expand.grid(smb = seq(0, 12, 0.25), c21 = seq(0, 30, 0.5)) %>%
      mutate(cls = case_when(
        smb < 2 & c21 > 10 ~ "smB− CD21low↑ (Freiburg Ia)\n비장비대·육아종·자가면역",
        smb < 2            ~ "smB− CD21low정상 (Freiburg Ib)",
        c21 > 10           ~ "smB+ CD21low↑",
        TRUE               ~ "smB+ CD21low정상 (경증 표현형)"))
    ggplot(grid, aes(smb, c21, fill = cls)) + geom_raster(alpha = 0.5) +
      geom_vline(xintercept = 2, linetype = 2) +
      geom_hline(yintercept = 10, linetype = 2) +
      annotate("point", x = min(l$SMBPCT, 12), y = min(l$CD21PCT, 30),
               size = 6, shape = 21, fill = "#C62828", colour = "black",
               stroke = 1.2) +
      annotate("text", x = min(l$SMBPCT, 12) + 0.4, y = min(l$CD21PCT, 30) + 1.4,
               label = "이 환자", fontface = "bold", size = 4) +
      scale_fill_manual(values = c("#FFCDD2", "#FFE0B2", "#C5CAE9", "#C8E6C9")) +
      labs(x = "전환 기억 B세포 (% of B)", y = "CD21low B세포 (% of B)",
           title = "EUROclass / Freiburg 분류 위치",
           subtitle = "smB− 그리고 CD21low>10%가 조절이상 표현형을 예측합니다") +
      THEME
  })

  ## ================= TAB 2 ==========================================
  output$p_pk <- renderPlot({
    d <- R()$both
    ggplot(d, aes(yr(d), IGG, colour = arm)) + geom_line(linewidth = 0.5) +
      geom_hline(yintercept = c(500, 700, 1000), linetype = 3,
                 colour = "grey45") +
      annotate("text", x = 0.3, y = c(520, 720, 1020), hjust = 0, size = 3,
               colour = "grey35", label = c("500", "700", "1000 mg/dL")) +
      scale_colour_manual(values = PAL) +
      labs(x = "시간 (년)", y = "혈청 IgG (mg/dL)",
           title = "면역글로불린 보충 PK — 전체 경과",
           subtitle = "정상상태까지 약 5 반감기(3-6개월). 지연 기간에는 내인성 IgG만 존재합니다") +
      THEME
  })

  output$p_pk_zoom <- renderPlot({
    d <- R()$both; te <- max(d$time)
    dz <- d %>% filter(time >= te - 84)
    ggplot(dz, aes((time - (te - 84)), IGG, colour = arm)) +
      geom_line(linewidth = 0.8) +
      scale_colour_manual(values = PAL) +
      labs(x = "마지막 12주 (일)", y = "혈청 IgG (mg/dL)",
           title = "정상상태 프로파일 (마지막 12주)",
           subtitle = "톱니파의 진폭이 곧 볼록성 벌점의 크기입니다") + THEME
  })

  output$t_pk <- renderTable({
    d <- R()$both; te <- max(d$time)
    d %>% filter(time >= te - 28) %>% group_by(arm) %>%
      summarise(`평균` = round(mean(IGG)), `trough` = round(min(IGG)),
                `peak` = round(max(IGG)),
                `변동폭` = round(max(IGG) - min(IGG)),
                `CL (dL/일)` = round(mean(CLTOT), 2),
                `폐렴/년` = round(mean(PNEU_YR), 4), .groups = "drop")
  })

  ## ================= TAB 3 ==========================================
  ER <- reactive({
    pp <- as.list(param(mod))
    cc <- seq(50, 1600, 10)
    brd <- pp$BRDPOOL*cc/(cc + pp$BRDC50) + pp$BRDENDOG*(3.2/pp$MEMREF)
    ops <- cc*(0.5 + 0.5*brd/(pp$BRDPOOL + pp$BRDENDOG))
    iga <- 0.24
    mp  <- ops*(0.35 + 0.65*iga)
    tibble(IgG = cc,
           pneu = pp$RMAXPNEU/(1 + (ops/pp$C50PNEU)^pp$HILLPNEU),
           sino = pp$RSINOFLR*(1 + 0.9*(1 - iga)) +
                  pp$RMAXSINO/(1 + (mp/pp$C50SINO)^pp$HILLSINO),
           lin  = pmax(0, 0.113*(1 - 0.0027*(cc - 500))))
  })

  output$p_er <- renderPlot({
    e <- ER(); l <- tail(R()$main, 1)
    ggplot(e, aes(IgG)) +
      geom_line(aes(y = pneu, colour = "Hill형 (이 모델)"), linewidth = 1) +
      geom_line(aes(y = lin, colour = "직선형 (Orange 외삽)"),
                linewidth = 0.8, linetype = 2) +
      geom_point(data = tibble(IgG = 500, y = 0.113), aes(y = y),
                 size = 3, colour = "black") +
      annotate("text", x = 560, y = 0.125, hjust = 0, size = 3.2,
               label = "앵커: 0.113/년 @ 500 mg/dL") +
      geom_vline(xintercept = l$IGG, colour = "#C62828", linetype = 3) +
      geom_vline(xintercept = 1400, colour = "grey60", linetype = 3) +
      annotate("text", x = 1395, y = 0.28, hjust = 1, size = 3.2,
               colour = "grey40", label = "직선형이 0에 닿는 지점\n(인공물)") +
      scale_colour_manual(values = c("#1565C0", "#C62828")) +
      labs(x = "혈청 IgG (mg/dL)", y = "폐렴 (사건/년)",
           title = "ARM 1 노출-반응",
           subtitle = "두 형태는 500-900 구간에서 거의 겹치고, 그 밖에서 전혀 다른 말을 합니다") +
      THEME
  })

  output$p_er2 <- renderPlot({
    e <- ER() %>% mutate(slope = c(NA, diff(log(pneu))/diff(IgG))) %>%
      mutate(pct100 = (exp(slope*100) - 1)*100)
    ggplot(e, aes(IgG, pct100)) + geom_line(linewidth = 1, colour = "#1565C0") +
      geom_hline(yintercept = -27, linetype = 2, colour = "#C62828") +
      annotate("text", x = 1450, y = -25, hjust = 1, size = 3.4,
               colour = "#C62828", label = "Orange 2010: −27% / 100 mg/dL") +
      labs(x = "혈청 IgG (mg/dL)", y = "100 mg/dL 증량당 폐렴률 변화 (%)",
           title = "기울기는 농도에 따라 변합니다",
           subtitle = "800-1000 mg/dL 부근이 수확 체감의 무릎입니다") + THEME
  })

  output$p_jensen <- renderPlot({
    e <- ER(); m <- 700; a <- 350
    seg <- tibble(x = c(m - a, m + a),
                  y = approx(e$IgG, e$pneu, c(m - a, m + a))$y)
    flat <- approx(e$IgG, e$pneu, m)$y
    ggplot(e %>% filter(IgG > 250, IgG < 1200), aes(IgG, pneu)) +
      geom_line(linewidth = 1, colour = "#1565C0") +
      geom_line(data = seg, aes(x, y), linetype = 2, colour = "#C62828") +
      geom_point(data = seg, aes(x, y), colour = "#C62828", size = 2.5) +
      annotate("point", x = m, y = flat, colour = "#2E7D32", size = 3.5) +
      annotate("point", x = m, y = mean(seg$y), colour = "#C62828", size = 3.5) +
      annotate("segment", x = m, xend = m, y = flat, yend = mean(seg$y),
               arrow = arrow(length = unit(2, "mm"), ends = "both"),
               colour = "black") +
      annotate("text", x = m + 25, y = (flat + mean(seg$y))/2, hjust = 0,
               size = 3.6, fontface = "bold",
               label = sprintf("볼록성 벌점\n%+.0f%%",
                               100*(mean(seg$y)/flat - 1))) +
      labs(x = "혈청 IgG (mg/dL)", y = "폐렴 (사건/년)",
           title = "왜 평탄한 프로파일이 같은 AUC에서 더 좋은가 — Jensen 부등식",
           subtitle = "곡선이 볼록하므로 동일 평균에서 변동하는 프로파일의 평균 위험이 더 높습니다") +
      THEME
  })

  ## ================= TAB 4 ==========================================
  output$p_inf <- renderPlot({
    d <- R()$both %>%
      select(time, arm, `폐렴 (/년)` = PNEU_YR,
             `부비동·기관지 감염 (/년)` = SINO_YR,
             `누적 폐렴` = PNEU_CUM, `누적 침습감염` = INVAS_CUM) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time/365, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) +
      labs(x = "시간 (년)", y = NULL, title = "감염 부담 (ARM 1)",
           subtitle = "지연 기간의 계단은 보충 시작 시점을 표시합니다") + THEME
  })

  output$p_divergence <- renderPlot({
    e <- ER() %>% mutate(pneu_rel = pneu/pneu[which.min(abs(IgG - 400))],
                         sino_rel = sino/sino[which.min(abs(IgG - 400))]) %>%
      filter(IgG >= 300, IgG <= 1500) %>%
      select(IgG, `폐렴 (상대)` = pneu_rel, `부비동염 (상대)` = sino_rel) %>%
      pivot_longer(-IgG)
    ggplot(e, aes(IgG, value, colour = name)) + geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c("#C62828", "#00838F")) +
      scale_y_continuous(limits = c(0, 1.35)) +
      labs(x = "혈청 IgG (mg/dL)", y = "IgG 400 mg/dL 대비 상대 발생률",
           title = "증량은 폐렴을 없애지만 부비동염은 거의 못 건드립니다",
           subtitle = "점막 바닥은 분비형 IgA 부재에서 오고, 어떤 용량으로도 도달할 수 없습니다") +
      THEME
  })

  ## ================= TAB 5 ==========================================
  output$p_ratchet <- renderPlot({
    d <- R()$main %>%
      select(time, `기도 집락화` = COLON, `기도 염증` = AIRINF,
             `엘라스타제 활성` = NEACT,
             `기관지확장증 (Reiff)` = BE,
             `연간 진행 (Reiff/년)` = BEFLUX, `점액 부하` = MUCUS) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/365, value)) +
      geom_line(linewidth = 0.8, colour = "#B71C1C") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (년)", y = NULL,
           title = "래칫의 각 단계",
           subtitle = "적절한 trough에 도달하면 진행률이 0으로 떨어지고, 축적된 손상은 남습니다") +
      THEME
  })

  output$p_delay <- renderPlot({
    withProgress(message = "지연 시나리오...", {
      p <- pars_of()
      dd <- lapply(c(1, 4, 7, 15), function(k) {
        d0 <- k*365; tau <- 28
        nn <- max(1, floor((input$years*365 - d0)/tau) + 1)
        e <- as.data.frame(ev(amt = 500*input$wt, cmt = "IGG_C", time = d0,
                              ii = tau, addl = nn - 1, evid = 1)); e$ID <- 1
        sim1(p, e, input$years, delta = 28) %>%
          mutate(lab = sprintf("지연 %2d년", k))
      }) %>% bind_rows()
    })
    ggplot(dd, aes(time/365, FEV1, colour = lab)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "시간 (년)", y = "FEV1 (% 예측치)",
           title = "동일한 치료, 다른 시작 시점",
           subtitle = "trough는 모두 같습니다. 차이는 전부 비가역적입니다") + THEME
  })

  output$p_lung <- renderPlot({
    d <- R()$both %>%
      select(time, arm, `FEV1 %pred` = FEV1, `FVC %pred` = FVC,
             `DLCO %pred` = DLCO) %>% pivot_longer(-c(time, arm))
    ggplot(d, aes(time/365, value, colour = arm, linetype = name)) +
      geom_line(linewidth = 0.8) + scale_colour_manual(values = PAL) +
      labs(x = "시간 (년)", y = "% 예측치", title = "폐기능",
           subtitle = "DLCO는 기관지확장증보다 GLILD/섬유화에 더 민감합니다") +
      THEME + theme(legend.box = "vertical")
  })

  ## ================= TAB 6 ==========================================
  output$p_arm2 <- renderPlot({
    d <- R()$both %>%
      select(time, arm, `GLILD 활성도` = GLILD, `간질 섬유화` = FIB,
             `비장 비대` = SPLEEN, `림프절 비대` = LAD,
             `ARM 2 총부담` = ARM2, `생존 확률` = SURV) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time/365, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) +
      labs(x = "시간 (년)", y = NULL,
           title = "ARM 2 — 보충요법이 건드리지 못하는 질환",
           subtitle = "두 곡선이 겹친다면, 그것이 바로 구조적 영(null)입니다") +
      THEME
  })

  output$p_cyto <- renderPlot({
    d <- R()$main %>%
      select(time, `혈소판 (10^9/L)` = PLT, `헤모글로빈 (g/dL)` = HGB,
             `알부민 (g/dL)` = ALB) %>% pivot_longer(-time)
    ggplot(d, aes(time/365, value)) + geom_line(linewidth = 0.9,
                                                colour = "#AD1457") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "시간 (년)", y = NULL, title = "세포감소증과 단백 소실") + THEME
  })

  output$p_glild <- renderPlot({
    d <- R()$main %>%
      select(time, `림프 집합체` = LYMPHAGG, `육아종` = GRAN,
             `장병증 활성도` = GUT, `단백소실 (PLE)` = PLE) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/365, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) +
      labs(x = "시간 (년)", y = NULL,
           title = "GLILD 구성요소와 장병증",
           subtitle = "PLE가 켜지면 IgG 청소율이 올라가 같은 용량의 trough가 무너집니다") +
      THEME
  })

  ## ================= TAB 7 ==========================================
  output$p_isb <- renderPlot({
    d <- R()$main %>%
      select(time, `순 면역억제 부담` = ISBURDEN,
             `감수성 배수` = SUSC, `폐렴 (/년)` = PNEU_YR,
             `IgG 청소율 (dL/일)` = CLTOT) %>% pivot_longer(-time)
    ggplot(d, aes(time/365, value)) + geom_line(linewidth = 0.9,
                                                colour = "#4527A0") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (년)", y = NULL,
           title = "ARM 2 치료가 ARM 1 질환으로 되돌아오는 경로",
           subtitle = "순 면역억제는 감수성을 올리고, 감수성은 감염률을 올립니다") +
      THEME
  })

  output$p_tradeoff <- renderPlot({
    withProgress(message = "리툭시맙 비대칭 계산...", {
      rtxev <- as.data.frame(ev(amt = 700, cmt = "RTX_C", time = 0,
                                ii = 7, addl = 3)); rtxev$ID <- 1
      nrm <- mod %>% param(HEALTHY = 1, DAMAGEON = 0) %>% init(HEALTHY_SS) %>%
        mrgsim(data = rtxev, end = 3*365, delta = 7) %>% as_tibble() %>%
        mutate(arm = "정상 숙주")
      sc <- as.data.frame(ev(amt = 500/4/0.73*input$wt, cmt = "IGG_SC",
                             time = 0, ii = 7, addl = 155))
      cv <- mod %>% param(CVID_GENO) %>% param(DAMAGEON = 0) %>%
        init(CVID_SS) %>%
        mrgsim(data = (bind_rows(rtxev, sc) %>% mutate(ID = 1) %>%
                       arrange(time)),
               end = 3*365, delta = 7) %>% as_tibble() %>%
        mutate(arm = "CVID (보충 중)")
      dd <- bind_rows(nrm, cv) %>% group_by(arm) %>%
        mutate(IgG_rel = 100*IGG/IGG[which.min(abs(time - 180))],
               B_rel = 100*BTOTAL/BTOTAL[1]) %>% ungroup()
    })
    ggplot(dd, aes(time/365)) +
      geom_line(aes(y = B_rel, colour = arm, linetype = "총 B세포"),
                linewidth = 0.9) +
      geom_line(aes(y = IgG_rel, colour = arm, linetype = "혈청 IgG"),
                linewidth = 0.9) +
      scale_colour_manual(values = c("#C62828", "#1565C0")) +
      labs(x = "리툭시맙 후 시간 (년)", y = "기저치 대비 (%)",
           title = "리툭시맙: 같은 고갈, 다른 대가",
           subtitle = "체액성 대가는 CVID에서 이미 지불되고 이미 보충되어 있습니다") +
      THEME + theme(legend.box = "vertical")
  })

  output$t_tradeoff <- renderTable({
    d <- R()$main; w <- d %>% filter(time >= max(time) - 28)
    tibble(항목 = c("순 면역억제 부담", "감수성 배수",
                    "누적 면역억제 노출", "누적 침습감염",
                    "ARM 2 부담", "모델 생존확률"),
           값 = sprintf("%.3f", c(mean(w$ISBURDEN), mean(w$SUSC),
                                 tail(d$ISCUM,1), tail(d$INVAS_CUM,1),
                                 tail(d$ARM2,1), tail(d$SURV,1))))
  })

  ## ================= TAB 8 ==========================================
  CMP <- reactive({
    p <- pars_of(); yrs <- input$years; wt <- input$wt
    mkiv <- function(dose, tau = 28, d0 = input$delay*365) {
      nn <- max(1, floor((yrs*365 - d0)/tau) + 1)
      e <- as.data.frame(ev(amt = dose*tau/28*wt, cmt = "IGG_C", time = d0,
                            ii = tau, addl = nn - 1, evid = 1)); e$ID <- 1; e
    }
    mksc <- function(dose, tau) {
      d0 <- input$delay*365
      nn <- max(1, floor((yrs*365 - d0)/tau) + 1)
      e <- as.data.frame(ev(amt = dose*tau/28/0.73*wt, cmt = "IGG_SC",
                            time = d0, ii = tau, addl = nn - 1, evid = 1))
      e$ID <- 1; e
    }
    withProgress(message = "비교 세트 시뮬레이션...", {
      switch(input$cmp_set,
        dose = lapply(c(300, 400, 500, 600, 800, 1000), function(k)
          sim1(p, mkiv(k), yrs, 14) %>% mutate(lab = paste0(k, " mg/kg q4w"))),
        route = list(
          sim1(p, mkiv(500, 28), yrs, 2) %>% mutate(lab = "IVIG q4주"),
          sim1(p, mkiv(500, 21), yrs, 2) %>% mutate(lab = "IVIG q3주"),
          sim1(p, mksc(500, 7),  yrs, 2) %>% mutate(lab = "SCIG 주1회"),
          sim1(p, mksc(500, 3.5),yrs, 2) %>% mutate(lab = "SCIG 주2회"),
          sim1(modifyList(p, list(HYAL = 1)), mksc(500, 28), yrs, 2) %>%
            mutate(lab = "fSCIG q4주")),
        delay = lapply(c(1, 4, 7, 15), function(k)
          sim1(p, mkiv(500, 28, k*365), yrs, 14) %>%
            mutate(lab = sprintf("지연 %2d년", k))),
        itp = { pi <- modifyList(p, list(DYSGENO = 0.6, KINAA = 0.030,
                                         KINLA = 0.004))
          t0 <- input$tstart*365; base <- mkiv(500)
          add <- function(x) (bind_rows(base, as.data.frame(x)) %>%
                              mutate(ID = 1) %>% arrange(time))
          list(
            sim1(pi, base, yrs, 14) %>% mutate(lab = "IgG 보충만"),
            sim1(pi, add(ev(amt = 20, cmt = "PRED_C", time = t0, ii = 1,
                            addl = yrs*365 - t0)), yrs, 14) %>%
              mutate(lab = "프레드니손"),
            sim1(pi, add(ev(amt = 700, cmt = "RTX_C", time = t0, ii = 365,
                            addl = 14)), yrs, 14) %>%
              mutate(lab = "리툭시맙"),
            sim1(modifyList(pi, list(SPLENEC = 1, TSPLENEC = t0)), base,
                 yrs, 14) %>% mutate(lab = "비장절제"),
            sim1(pi, add(ev(amt = 50, cmt = "ELT_C", time = t0, ii = 1,
                            addl = yrs*365 - t0)), yrs, 14) %>%
              mutate(lab = "엘트롬보팍")) },
        glild = { pg <- modifyList(p, list(DYSGENO = 0.7, KINLA = 0.003))
          t0 <- input$tstart*365; base <- mkiv(500)
          add <- function(x) (bind_rows(base, as.data.frame(x)) %>%
                              mutate(ID = 1) %>% arrange(time))
          list(
            sim1(pg, base, yrs, 14) %>% mutate(lab = "IgG 보충만"),
            sim1(pg, add(ev(amt = 20, cmt = "PRED_C", time = t0, ii = 1,
                            addl = yrs*365 - t0)), yrs, 14) %>%
              mutate(lab = "프레드니손"),
            sim1(modifyList(pg, list(AZAON = 1, AZAT0 = t0)),
                 add(ev(amt = 700, cmt = "RTX_C", time = t0, ii = 182,
                        addl = 30)), yrs, 14) %>%
              mutate(lab = "리툭시맙+AZA"),
            sim1(modifyList(pg, list(JAKON = 1, JAKT0 = t0)), base,
                 yrs, 14) %>% mutate(lab = "JAK 억제제")) }
      ) %>% bind_rows()
    })
  })

  output$p_cmp <- renderPlot({
    d <- CMP(); yv <- input$cmp_y
    ggplot(d, aes(time/365, .data[[yv]], colour = lab)) +
      geom_line(linewidth = 0.85) +
      scale_colour_manual(values = rep(PAL, 3)) +
      labs(x = "시간 (년)", y = yv,
           title = paste("시나리오 비교 —", yv)) + THEME
  })

  output$t_cmp <- renderDT({
    d <- CMP(); te <- max(d$time)
    d %>% group_by(lab) %>%
      summarise(`IgG trough` = round(min(IGG[time >= te - 28])),
                `변동폭` = round(diff(range(IGG[time >= te - 28]))),
                `폐렴/년` = round(mean(PNEU_YR[time >= te - 28]), 4),
                `부비동염/년` = round(mean(SINO_YR[time >= te - 28]), 2),
                `Reiff` = round(last(BE), 2), `FEV1` = round(last(FEV1), 1),
                `DLCO` = round(last(DLCO), 1),
                `GLILD` = round(last(GLILD), 3),
                `혈소판` = round(last(PLT)),
                `ARM2` = round(last(ARM2), 2),
                `QoL` = round(last(QOL), 3),
                `생존` = round(last(SURV), 3), .groups = "drop") %>%
      datatable(options = list(dom = "t", pageLength = 12), rownames = FALSE)
  })

  ## ================= TAB 9 ==========================================
  output$p_bio <- renderPlot({
    d <- R()$main %>%
      select(time, `혈청 IgG (mg/dL)` = IGG,
             `유효 옵소닌 (mg/dL 등가)` = OPSONIN,
             `점막 방어 (mg/dL 등가)` = MUCPROT,
             `가용성 BAFF (x정상)` = BAFF,
             `전환 기억 B (%)` = SMBPCT, `CD21low B (%)` = CD21PCT) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/365, value)) + geom_line(linewidth = 0.85,
                                                colour = "#00695C") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "시간 (년)", y = NULL, title = "바이오마커 궤적",
           subtitle = "옵소닌과 점막 방어의 간격이 곧 IgA 공백입니다") + THEME
  })

  output$p_bcell <- renderPlot({
    d <- R()$main %>%
      select(time, `나이브 B` = NAIVEPCT, `전환 기억 B` = SMBPCT,
             `CD21low B` = CD21PCT) %>% pivot_longer(-time)
    ggplot(d, aes(time/365, value, fill = name)) +
      geom_area(position = "stack", alpha = 0.85) +
      scale_fill_manual(values = c("#90CAF9", "#EF9A9A", "#A5D6A7")) +
      labs(x = "시간 (년)", y = "B세포 아형 (% of B)",
           title = "B세포 구획 구성") + THEME
  })

  output$p_baff <- renderPlot({
    d <- R()$main %>%
      select(time, `BAFF` = BAFF, `총 B세포 (/uL)` = BTOTAL) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/365, value)) + geom_line(linewidth = 0.9,
                                                colour = "#00838F") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "시간 (년)", y = NULL,
           title = "BAFF와 B세포 sink",
           subtitle = "sink가 줄면 BAFF가 오릅니다 — 리툭시맙 후 2-4배") + THEME
  })

  ## ================= TAB 10 =========================================
  output$t_anchor <- renderDT({
    datatable(BUILD$anchors, options = list(dom = "t", pageLength = 20),
              rownames = FALSE)
  })
  output$t_diag <- renderDT({
    BUILD$diagnostics %>%
      mutate(pass = ifelse(pass, "PASS", "FAIL")) %>%
      datatable(options = list(dom = "tp", pageLength = 15), rownames = FALSE)
  })
}

shinyApp(ui, server)
