##############################################################################
# Fabry Disease QSP — Shiny Interactive Dashboard
# 파브리병 정량적 시스템 약리학 모델 — Shiny 대시보드
# 8 Tabs: 환자 프로파일 · PK/효소 · Gb3 동역학 · 신장 기능 · 심장 기능 ·
#         치료 시나리오 비교 · 바이오마커 패널 · 가상 환자 집단
##############################################################################

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

# ── Colour palette ────────────────────────────────────────────────────────
scenario_colors <- c(
  "자연경과"           = "#616161",
  "아갈시다제 베타"    = "#1565C0",
  "아갈시다제 알파"    = "#0288D1",
  "미갈라스타트"       = "#2E7D32",
  "페구니알시다제 알파" = "#7B1FA2",
  "ERT+루세라스탓"    = "#E65100"
)

##############################################################################
# CORE SIMULATION FUNCTION (simplified analytical / semi-ODE)
##############################################################################
simulate_fabry <- function(
    phenotype      = "classic",   # "classic" | "late"
    is_amenable    = FALSE,       # migalastat amenable mutation
    treatment      = "none",      # "none"|"agaB"|"agaA"|"migalastat"|"pegu"|"combo"
    agaB_dose      = 1.0,         # mg/kg
    agaA_dose      = 0.2,         # mg/kg
    mig_dose       = 150,         # mg QOD
    luc_dose       = 1000,        # mg TID
    sim_years      = 5,
    age_start      = 30,
    sex            = "male",
    eGFR0          = 90
) {
  t <- seq(0, sim_years * 365, by = 7)  # weekly steps
  n <- length(t)

  # Baseline enzyme activity
  E0 <- if (phenotype == "classic") 0.05 else 2.5   # nmol/h/mg

  # Treatment-specific enzyme restoration
  E_max_treat <- switch(treatment,
    "none"       = 0,
    "agaB"       = 70 * agaB_dose,   # ~70 for 1 mg/kg
    "agaA"       = 15 * agaA_dose / 0.2, # ~15 for 0.2 mg/kg
    "migalastat" = if (is_amenable) 6.0 else 0,
    "pegu"       = 80 * agaB_dose,   # similar to agaB but longer duration
    "combo"      = 70 * agaB_dose,   # ERT component
    0
  )
  E_total <- E0 + E_max_treat

  # GCS inhibition (lucerastat component)
  GCS_inhib_frac <- if (treatment %in% c("combo")) 0.40 else 0.0

  # Enzyme activity over time (approach to steady-state, weeks)
  t_ss <- 4 * 7  # 4 weeks to reach ~63% of SS enzyme effect
  E_t  <- E_total - (E_total - E0) * exp(-t / t_ss)

  # Lyso-Gb3 dynamics
  LGB3_0 <- if (phenotype == "classic") 12.0 else 4.5   # μg/L
  LGB3_ss_treat <- LGB3_0 * (1 - min(0.85, (E_t[n] - E0) / 40)) *
                   (1 - GCS_inhib_frac * 0.35)
  LGB3_t <- LGB3_0 - (LGB3_0 - pmax(LGB3_ss_treat, 1.5)) *
            (1 - exp(-t / (90 * 7)))  # 90 weeks to reach new SS

  # Gb3 plasma
  Gb3_0 <- if (phenotype == "classic") 0.8 else 0.4
  Gb3_ss <- Gb3_0 * (1 - min(0.75, (E_total - E0) / 50))
  Gb3_t  <- Gb3_0 - (Gb3_0 - pmax(Gb3_ss, 0.05)) * (1 - exp(-t / (78 * 7)))

  # eGFR decline
  # Natural history: ~ -3 to -6 mL/min/1.73m²/yr (classic male)
  eGFR_annual_loss_base <- if (sex == "male" & phenotype == "classic") -4.5
                           else if (sex == "female") -1.5 else -2.0

  eGFR_treat_benefit <- switch(treatment,
    "agaB"  = 0.65,  # 65% reduction in GFR loss (Banikazemi 2007)
    "agaA"  = 0.50,
    "pegu"  = 0.70,
    "migalastat" = if (is_amenable) 0.60 else 0,
    "combo" = 0.75,
    "none"  = 0
  )
  eGFR_annual_loss <- eGFR_annual_loss_base * (1 - eGFR_treat_benefit)
  eGFR_t <- pmax(0, eGFR0 + eGFR_annual_loss * t / 365)

  # LVMi
  LVMi0 <- if (phenotype == "classic" & sex == "male") 148 else 120
  LVMi_rate <- if (phenotype == "classic") 3.0 else 1.5  # g/m²/yr natural
  LVMi_benefit <- switch(treatment,
    "agaB"  = 0.80, "agaA" = 0.65, "pegu" = 0.85,
    "migalastat" = if (is_amenable) 0.60 else 0,
    "combo" = 0.90, "none" = 0
  )
  LVMi_t <- LVMi0 + (LVMi_rate * (1 - LVMi_benefit)) * t / 365

  # Neuropathic Pain (BPI-SF)
  Pain0 <- 6.5
  Pain_benefit <- switch(treatment,
    "agaB"  = 0.35, "agaA" = 0.30, "pegu" = 0.40,
    "migalastat" = if (is_amenable) 0.45 else 0,
    "combo" = 0.45, "none" = 0
  )
  Pain_t <- Pain0 * (1 - Pain_benefit * (1 - exp(-t / (52 * 7))))

  # UPCR
  UPCR0 <- if (phenotype == "classic") 350 else 150
  UPCR_benefit <- switch(treatment,
    "agaB" = 0.40, "agaA" = 0.30, "pegu" = 0.45,
    "migalastat" = if (is_amenable) 0.50 else 0,
    "combo" = 0.55, "none" = 0
  )
  UPCR_t <- UPCR0 * (1 - UPCR_benefit * (1 - exp(-t / (24 * 7))))

  data.frame(
    time     = t,
    Year     = t / 365,
    LysoGb3  = pmax(0.5, LGB3_t),
    Gb3_plasma = pmax(0.02, Gb3_t),
    eGFR     = eGFR_t,
    LVMi     = LVMi_t,
    Pain     = pmax(1, Pain_t),
    UPCR     = pmax(30, UPCR_t),
    EnzymeActivity = E_t,
    Inflammation = 1.8 - 0.9 * (E_t - E0) / (E_total - E0 + 0.1)
  )
}

##############################################################################
# UI
##############################################################################
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "파브리병 QSP 대시보드"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("1. 환자 프로파일",      tabName = "patient",   icon = icon("user")),
      menuItem("2. PK / 효소 활성",     tabName = "pk",        icon = icon("flask")),
      menuItem("3. Gb3 동역학",         tabName = "gb3",       icon = icon("dna")),
      menuItem("4. 신장 기능",          tabName = "renal",     icon = icon("filter")),
      menuItem("5. 심장 기능",          tabName = "cardiac",   icon = icon("heart")),
      menuItem("6. 시나리오 비교",      tabName = "scenarios", icon = icon("chart-line")),
      menuItem("7. 바이오마커 패널",    tabName = "biomarkers",icon = icon("vial")),
      menuItem("8. 가상 환자 집단",     tabName = "population",icon = icon("users"))
    )
  ),

  dashboardBody(
    tabItems(

      ##########################################################################
      # TAB 1: PATIENT PROFILE
      ##########################################################################
      tabItem("patient",
        fluidRow(
          box(title = "파브리병 환자 설정", width = 4, status = "primary", solidHeader = TRUE,
            selectInput("phenotype", "표현형 (Phenotype)",
                        choices = c("고전형 (Classic)" = "classic",
                                    "후기 발현형 (Late-onset)" = "late")),
            radioButtons("sex", "성별", choices = c("남성" = "male", "여성" = "female"), inline = TRUE),
            sliderInput("age_start", "진단 연령", min = 5, max = 70, value = 25, step = 1),
            sliderInput("eGFR0", "기저 eGFR (mL/min/1.73m²)", min = 30, max = 120, value = 90),
            sliderInput("sim_years", "시뮬레이션 기간 (년)", min = 1, max = 15, value = 5)
          ),
          box(title = "돌연변이 및 치료 적합성", width = 4, status = "warning", solidHeader = TRUE,
            checkboxInput("is_amenable", "미갈라스타트 적합 변이 보유 (~40%)", value = FALSE),
            helpText("미갈라스타트(Galafold)는 HEK293 세포 검정으로 적합성이 확인된\n변이에서만 효과적입니다."),
            selectInput("treatment", "치료 방법 선택",
              choices = c(
                "치료 없음 (자연경과)" = "none",
                "아갈시다제 베타 (Fabrazyme)" = "agaB",
                "아갈시다제 알파 (Replagal)" = "agaA",
                "미갈라스타트 (Galafold)" = "migalastat",
                "페구니알시다제 알파 (Elfabrio)" = "pegu",
                "ERT + 루세라스탓 (병용)" = "combo"
              ))
          ),
          box(title = "용량 조절", width = 4, status = "info", solidHeader = TRUE,
            sliderInput("agaB_dose", "아갈시다제 베타 (mg/kg)", min = 0.3, max = 2.0, value = 1.0, step = 0.1),
            sliderInput("agaA_dose", "아갈시다제 알파 (mg/kg)", min = 0.1, max = 0.6, value = 0.2, step = 0.05),
            sliderInput("mig_dose",  "미갈라스타트 (mg, QOD)", min = 50, max = 300, value = 150, step = 25),
            sliderInput("luc_dose",  "루세라스탓 (mg, TID)",   min = 250, max = 1500, value = 1000, step = 250)
          )
        ),
        fluidRow(
          box(title = "파브리병 개요 — 치료 표적 경로", width = 12, status = "success",
            p(strong("파브리병(Fabry disease)"), "은 X-연관 열성 유전 질환으로,",
              strong("GLA 유전자 변이"), "에 의해 α-갈락토시다제 A(α-Gal A) 효소가 결핍됩니다."),
            p("결핍된 효소로 인해 글리코스핑고지질(Gb3 및 lyso-Gb3)이 신장, 심장, CNS, 피부 등에
              점진적으로 축적되어 다장기 손상을 유발합니다."),
            tags$ul(
              tags$li(strong("효소대체요법(ERT):"), " 아갈시다제 베타(Fabrazyme) 1 mg/kg Q2W,
                      아갈시다제 알파(Replagal) 0.2 mg/kg Q2W,
                      페구니알시다제 알파(Elfabrio) 1 mg/kg Q4W"),
              tags$li(strong("샤페론 요법:"), " 미갈라스타트(Galafold) 150 mg QOD —
                      적합 변이(~40%)에만 유효, 경구 복용 가능"),
              tags$li(strong("기질감소요법(SRT):"), " 루세라스탓(Lucerastat) 1000 mg TID —
                      GCS 억제로 Gb3 전구체 감소 (MODIFY 임상시험)")
            )
          )
        )
      ),

      ##########################################################################
      # TAB 2: PK / ENZYME ACTIVITY
      ##########################################################################
      tabItem("pk",
        fluidRow(
          box(title = "α-갈락토시다제 A 효소 활성 경과", width = 6, status = "primary", solidHeader = TRUE,
            plotOutput("plot_enzyme")
          ),
          box(title = "치료별 효소 활성 복원 메커니즘", width = 6, status = "info",
            h4("ERT (효소대체요법)"),
            p("IV 주입 후 M6P 수용체를 통해 리소솜으로 전달됩니다. 아갈시다제 베타는
              t½ ~45분의 짧은 혈장 반감기를 가지며, 페구니알시다제 알파(PEGylation)는
              t½ ~80시간으로 연장됩니다."),
            h4("미갈라스타트 (샤페론)"),
            p("소분자 샤페론으로 변이 α-Gal A 단백질의 잘못된 접힘을 방지하고
              리소솜 전달을 개선합니다. 적합 변이(GalaxIES 검정)에서만 유효합니다."),
            h4("루세라스탓 (SRT)"),
            p("글루코실세라미드 합성효소(GCS/UGCG) 억제를 통해 Gb3 전구체 합성을 줄입니다.
              MODIFY 임상시험에서 Gb3 감소를 확인하였습니다.")
          )
        ),
        fluidRow(
          box(title = "치료 PK 파라미터 요약", width = 12,
            DTOutput("pk_param_table")
          )
        )
      ),

      ##########################################################################
      # TAB 3: Gb3 DYNAMICS
      ##########################################################################
      tabItem("gb3",
        fluidRow(
          box(title = "혈장 Lyso-Gb3 (민감한 치료 반응 바이오마커)", width = 6,
              status = "danger", solidHeader = TRUE,
            plotOutput("plot_lysoGb3")
          ),
          box(title = "혈장 Gb3 경과", width = 6, status = "warning", solidHeader = TRUE,
            plotOutput("plot_plasma_Gb3")
          )
        ),
        fluidRow(
          box(title = "Lyso-Gb3 임상 의의", width = 12,
            p(strong("혈장 Lyso-Gb3"), "는 현재 가장 민감한 파브리병 바이오마커로,
              기저 질환 활성도 및 치료 반응을 반영합니다."),
            p("정상 참고 범위: 남성 <2 μg/L, 여성 <2 μg/L (일부 기관 기준 <1.4 μg/L)"),
            p("고전형 남성 미치료 시: 평균 30–80 μg/L"),
            tags$ul(
              tags$li("Aerts JM et al. (2008) PNAS: Lyso-Gb3 발견 및 임상 의의"),
              tags$li("Smid BE et al. (2014) Mol Genet Metab: ERT 치료 후 급격한 감소"),
              tags$li("Nowak A et al. (2017) Kidney Int: eGFR 기울기와 상관")
            )
          )
        )
      ),

      ##########################################################################
      # TAB 4: RENAL FUNCTION
      ##########################################################################
      tabItem("renal",
        fluidRow(
          box(title = "eGFR 경과 (신장 기능 보존)", width = 6,
              status = "danger", solidHeader = TRUE,
            plotOutput("plot_eGFR")
          ),
          box(title = "UPCR 경과 (단백뇨)", width = 6, status = "warning", solidHeader = TRUE,
            plotOutput("plot_UPCR")
          )
        ),
        fluidRow(
          box(title = "파브리 신장병 병기", width = 12,
            p(strong("자연 경과:"), "치료 없이 고전형 남성에서 eGFR이 연간 약 -3~-12 mL/min/1.73m²
              감소합니다(Warnock 2012)."),
            p(strong("ERT 신장 보호:"), "Banikazemi et al. (2007) Annals of Internal Medicine에서
              복합 신장 사건 발생률을 61% 감소시킴."),
            p(strong("치료 목표:"), "eGFR 기울기를 정상 노화에 해당하는 -1 mL/min/1.73m²/yr 수준으로 개선."),
            h5("CKD 병기 기준 (KDIGO 2012):"),
            tags$ul(
              tags$li("G1: ≥90 mL/min/1.73m²"),
              tags$li("G2: 60–89"),
              tags$li("G3a: 45–59"),
              tags$li("G3b: 30–44"),
              tags$li("G4: 15–29"),
              tags$li("G5 (ESRD): <15")
            )
          )
        )
      ),

      ##########################################################################
      # TAB 5: CARDIAC FUNCTION
      ##########################################################################
      tabItem("cardiac",
        fluidRow(
          box(title = "LVMi 경과 (심장 비후 역전)", width = 7,
              status = "primary", solidHeader = TRUE,
            plotOutput("plot_LVMi")
          ),
          box(title = "신경병성 통증 (BPI-SF)", width = 5,
              status = "warning", solidHeader = TRUE,
            plotOutput("plot_Pain")
          )
        ),
        fluidRow(
          box(title = "파브리 심장병증 — 핵심 사항", width = 12,
            p(strong("파브리 심장병증"), "은 좌심실 비후(LVH), 이완기 기능 장애, 전도 이상,
              부정맥 및 심근 섬유화(LGE)를 포함합니다."),
            p("남성 정상 LVMi 상한값: 115 g/m², 여성: 95 g/m²"),
            tags$ul(
              tags$li("Weidemann F et al. (2003) Circulation: ERT 치료 후 LVMi 감소"),
              tags$li("Desnick RJ (2010) J Inherit Metab Dis: 심장 Gb3 ERT 반응성"),
              tags$li("Linhart A (2020) Eur Heart J: 파브리 심장병증 자연 경과")
            )
          )
        )
      ),

      ##########################################################################
      # TAB 6: SCENARIO COMPARISON
      ##########################################################################
      tabItem("scenarios",
        fluidRow(
          box(title = "6가지 치료 시나리오 비교 — Lyso-Gb3", width = 6,
              status = "primary", solidHeader = TRUE,
            plotOutput("plot_scen_lgb3")
          ),
          box(title = "6가지 치료 시나리오 비교 — eGFR", width = 6,
              status = "danger", solidHeader = TRUE,
            plotOutput("plot_scen_eGFR")
          )
        ),
        fluidRow(
          box(title = "시나리오 비교 — LVMi", width = 6,
              status = "warning", solidHeader = TRUE,
            plotOutput("plot_scen_LVMi")
          ),
          box(title = "치료 시나리오별 5년 예상 결과", width = 6,
            DTOutput("outcome_table")
          )
        )
      ),

      ##########################################################################
      # TAB 7: BIOMARKER PANEL
      ##########################################################################
      tabItem("biomarkers",
        fluidRow(
          box(title = "바이오마커 대시보드", width = 12, status = "success",
            fluidRow(
              valueBoxOutput("vb_lysoGb3"),
              valueBoxOutput("vb_eGFR"),
              valueBoxOutput("vb_LVMi")
            ),
            fluidRow(
              valueBoxOutput("vb_UPCR"),
              valueBoxOutput("vb_pain"),
              valueBoxOutput("vb_enzyme")
            )
          )
        ),
        fluidRow(
          box(title = "바이오마커 해석 가이드", width = 12,
            DTOutput("biomarker_guide_table")
          )
        )
      ),

      ##########################################################################
      # TAB 8: VIRTUAL POPULATION
      ##########################################################################
      tabItem("population",
        fluidRow(
          box(title = "가상 환자 집단 설정 (Monte Carlo)", width = 4,
              status = "info", solidHeader = TRUE,
            sliderInput("n_patients", "시뮬레이션 환자 수", min = 20, max = 200, value = 50),
            sliderInput("cv_enzyme", "효소 활성 CV% (개인간 변이)", min = 10, max = 80, value = 40),
            sliderInput("cv_eGFR", "기저 eGFR CV%", min = 10, max = 50, value = 25),
            selectInput("pop_treatment", "치료법",
              choices = c("none", "agaB", "agaA", "migalastat", "pegu", "combo")),
            actionButton("run_pop", "집단 시뮬레이션 실행", class = "btn-primary")
          ),
          box(title = "eGFR 집단 분포 (5년 예측)", width = 8, status = "success",
            plotOutput("plot_pop_eGFR")
          )
        ),
        fluidRow(
          box(title = "Lyso-Gb3 반응 분포", width = 6,
            plotOutput("plot_pop_lgb3")
          ),
          box(title = "ERT 반응자 vs 비반응자 비율", width = 6,
            plotOutput("plot_responder")
          )
        )
      )
    )
  )
)

##############################################################################
# SERVER
##############################################################################
server <- function(input, output, session) {

  # ── Reactive simulation (single patient) ─────────────────────────────
  sim_data <- reactive({
    simulate_fabry(
      phenotype   = input$phenotype,
      is_amenable = input$is_amenable,
      treatment   = input$treatment,
      agaB_dose   = input$agaB_dose,
      agaA_dose   = input$agaA_dose,
      mig_dose    = input$mig_dose,
      luc_dose    = input$luc_dose,
      sim_years   = input$sim_years,
      age_start   = input$age_start,
      sex         = input$sex,
      eGFR0       = input$eGFR0
    )
  })

  # ── All 6 scenarios for comparison ───────────────────────────────────
  all_scen_data <- reactive({
    treatments <- c("none", "agaB", "agaA", "migalastat", "pegu", "combo")
    labels     <- c("자연경과", "아갈시다제 베타", "아갈시다제 알파",
                    "미갈라스타트", "페구니알시다제 알파", "ERT+루세라스탓")
    purrr::map2_dfr(treatments, labels, function(trt, lbl) {
      simulate_fabry(
        phenotype   = input$phenotype,
        is_amenable = TRUE,  # assume amenable for migalastat comparison
        treatment   = trt,
        sim_years   = input$sim_years,
        eGFR0       = input$eGFR0
      ) %>% mutate(Scenario = lbl)
    })
  })

  # ── TAB 2: Enzyme Activity ────────────────────────────────────────────
  output$plot_enzyme <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(Year, EnzymeActivity)) +
      geom_line(color = "#1565C0", size = 1.3) +
      geom_hline(yintercept = 8.0, linetype = "dashed", color = "gray50") +
      annotate("text", x = max(d$Year) * 0.8, y = 8.5,
               label = "정상 범위 (≥8 nmol/h/mg)", size = 3.5) +
      labs(title = "α-갈락토시다제 A 활성 경과",
           x = "시간 (년)", y = "α-Gal A (nmol/h/mg)") +
      theme_bw(base_size = 13)
  })

  output$pk_param_table <- renderDT({
    data.frame(
      약물 = c("아갈시다제 베타 (Fabrazyme)", "아갈시다제 알파 (Replagal)",
               "페구니알시다제 알파 (Elfabrio)", "미갈라스타트 (Galafold)", "루세라스탓"),
      용법 = c("1 mg/kg IV Q2W", "0.2 mg/kg IV Q2W", "1 mg/kg IV Q4W",
               "150 mg PO QOD", "1000 mg PO TID"),
      `t½` = c("~45분", "~45–110분", "~80시간", "~3.5시간", "~8시간"),
      `F(%)` = c("100 (IV)", "100 (IV)", "100 (IV)", "75%", "65%"),
      작용기전 = c("M6P 수용체 → 리소솜 전달", "M6P 수용체 → 리소솜 전달",
                   "PEGylation → 연장 t½", "α-Gal A 샤페론 (적합 변이)", "GCS 억제 (SRT)"),
      `EC50/IC50` = c("약 1.5 ng/mL (Km)", "약 1.5 ng/mL (Km)", "약 0.8 ng/mL", "0.25 μg/mL", "0.18 μg/mL"),
      `Emax_Gb3감소(%)` = c("~80%", "~70%", "~85%", "~60% (적합)", "~40%"),
      임상시험 = c("FABRY-001 (Eng 2001 NEJM)", "Schiffmann 2001 Ann Intern Med",
                   "BRIGHT (Schiffmann 2021 JAMA)", "ATTRACT (Germain 2016 NEJM)",
                   "MODIFY (Lenders 2022 Lancet DE)")
    ) %>% datatable(options = list(pageLength = 5, scrollX = TRUE), rownames = FALSE)
  })

  # ── TAB 3: Gb3 Dynamics ───────────────────────────────────────────────
  output$plot_lysoGb3 <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(Year, LysoGb3)) +
      geom_line(color = "#C62828", size = 1.3) +
      geom_hline(yintercept = 2.0, linetype = "dashed", color = "forestgreen") +
      annotate("text", x = max(d$Year) * 0.75, y = 2.5,
               label = "정상 상한값 (<2 μg/L)", color = "forestgreen", size = 3.5) +
      labs(title = "혈장 Lyso-Gb3 경과",
           x = "시간 (년)", y = "Lyso-Gb3 (μg/L)") +
      theme_bw(base_size = 13)
  })

  output$plot_plasma_Gb3 <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(Year, Gb3_plasma)) +
      geom_line(color = "#E65100", size = 1.3) +
      geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray50") +
      labs(title = "혈장 Gb3 (트리헥소실세라미드)",
           x = "시간 (년)", y = "Gb3 (μg/mL)") +
      theme_bw(base_size = 13)
  })

  # ── TAB 4: Renal ──────────────────────────────────────────────────────
  output$plot_eGFR <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(Year, eGFR)) +
      geom_line(color = "#AD1457", size = 1.3) +
      geom_hline(yintercept = 60, linetype = "dashed", color = "red") +
      annotate("text", x = max(d$Year) * 0.7, y = 63, label = "CKD G3 경계 (60)",
               color = "red", size = 3.5) +
      labs(title = "eGFR 경과 (CKD-EPI)",
           x = "시간 (년)", y = "eGFR (mL/min/1.73m²)") +
      ylim(0, 120) + theme_bw(base_size = 13)
  })

  output$plot_UPCR <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(Year, UPCR)) +
      geom_line(color = "#6A1B9A", size = 1.3) +
      geom_hline(yintercept = 300, linetype = "dashed", color = "orange") +
      labs(title = "UPCR (단백뇨)",
           x = "시간 (년)", y = "UPCR (mg/g Cr)") +
      theme_bw(base_size = 13)
  })

  # ── TAB 5: Cardiac ────────────────────────────────────────────────────
  output$plot_LVMi <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(Year, LVMi)) +
      geom_line(color = "#1565C0", size = 1.3) +
      geom_hline(yintercept = 115, linetype = "dashed", color = "blue") +
      annotate("text", x = max(d$Year) * 0.7, y = 118, label = "남성 정상 상한값 (115 g/m²)",
               color = "blue", size = 3.5) +
      labs(title = "LVMi (좌심실 질량 지수)",
           x = "시간 (년)", y = "LVMi (g/m²)") +
      theme_bw(base_size = 13)
  })

  output$plot_Pain <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(Year, Pain)) +
      geom_line(color = "#FF6F00", size = 1.3) +
      ylim(0, 10) +
      labs(title = "신경병성 통증 (BPI-SF)",
           x = "시간 (년)", y = "통증 점수 (0–10)") +
      theme_bw(base_size = 13)
  })

  # ── TAB 6: Scenario Comparison ────────────────────────────────────────
  output$plot_scen_lgb3 <- renderPlot({
    d <- all_scen_data()
    ggplot(d, aes(Year, LysoGb3, color = Scenario)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = scenario_colors) +
      geom_hline(yintercept = 2.0, linetype = "dashed") +
      labs(title = "Lyso-Gb3 시나리오 비교",
           x = "시간 (년)", y = "Lyso-Gb3 (μg/L)") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom", legend.title = element_blank())
  })

  output$plot_scen_eGFR <- renderPlot({
    d <- all_scen_data()
    ggplot(d, aes(Year, eGFR, color = Scenario)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = scenario_colors) +
      geom_hline(yintercept = 60, linetype = "dashed", color = "red") +
      labs(title = "eGFR 시나리오 비교",
           x = "시간 (년)", y = "eGFR (mL/min/1.73m²)") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom", legend.title = element_blank())
  })

  output$plot_scen_LVMi <- renderPlot({
    d <- all_scen_data()
    ggplot(d, aes(Year, LVMi, color = Scenario)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = scenario_colors) +
      labs(title = "LVMi 시나리오 비교",
           x = "시간 (년)", y = "LVMi (g/m²)") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom", legend.title = element_blank())
  })

  output$outcome_table <- renderDT({
    d <- all_scen_data()
    last_yr <- max(d$Year)
    d %>% filter(abs(Year - last_yr) < 0.15) %>%
      group_by(Scenario) %>%
      summarise(
        `Lyso-Gb3 (μg/L)` = round(mean(LysoGb3), 1),
        `eGFR`             = round(mean(eGFR), 1),
        `LVMi (g/m²)`      = round(mean(LVMi), 1),
        `통증 BPI`          = round(mean(Pain), 1),
        .groups = "drop"
      ) %>%
      datatable(rownames = FALSE,
                options = list(pageLength = 8, dom = "t"))
  })

  # ── TAB 7: Biomarker Valuebox ─────────────────────────────────────────
  output$vb_lysoGb3 <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(round(d$LysoGb3, 1), "Lyso-Gb3 (μg/L)", icon = icon("vial"),
             color = if (d$LysoGb3 < 2) "green" else if (d$LysoGb3 < 10) "yellow" else "red")
  })
  output$vb_eGFR <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(round(d$eGFR, 1), "eGFR (mL/min/1.73m²)", icon = icon("filter"),
             color = if (d$eGFR >= 60) "green" else if (d$eGFR >= 30) "yellow" else "red")
  })
  output$vb_LVMi <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(round(d$LVMi, 1), "LVMi (g/m²)", icon = icon("heart"),
             color = if (d$LVMi <= 115) "green" else if (d$LVMi <= 140) "yellow" else "red")
  })
  output$vb_UPCR <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(round(d$UPCR, 0), "UPCR (mg/g)", icon = icon("tint"),
             color = if (d$UPCR < 150) "green" else if (d$UPCR < 500) "yellow" else "red")
  })
  output$vb_pain <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(round(d$Pain, 1), "통증 BPI-SF (0–10)", icon = icon("bolt"),
             color = if (d$Pain < 3) "green" else if (d$Pain < 6) "yellow" else "red")
  })
  output$vb_enzyme <- renderValueBox({
    d <- tail(sim_data(), 1)
    valueBox(round(d$EnzymeActivity, 2), "α-Gal A (nmol/h/mg)", icon = icon("flask"),
             color = if (d$EnzymeActivity >= 5) "green" else if (d$EnzymeActivity >= 1) "yellow" else "red")
  })

  output$biomarker_guide_table <- renderDT({
    data.frame(
      바이오마커 = c("혈장 Lyso-Gb3", "혈장 Gb3", "백혈구 α-Gal A 활성",
                     "소변 Gb3", "eGFR (CKD-EPI)", "LVMi", "UPCR"),
      단위 = c("μg/L", "μg/mL", "nmol/h/mg", "nmol/mg Cr",
               "mL/min/1.73m²", "g/m²", "mg/g Cr"),
      정상범위 = c("<2.0", "<0.1", "≥8 (남)", "낮음", "≥60 (G1-2)", "≤115 (남)", "<150"),
      임상의의 = c("치료 반응 민감 지표 (ERT/샤페론)", "장기 Gb3 축적 간접 지표",
                   "진단 기준 (남성)", "신장 Gb3 침착 직접 측정", "신기능 보존 목표",
                   "심장 비후 역전 목표", "사구체 단백질 여과 장벽"),
      치료반응성 = c("빠름 (수주)", "중간 (수개월)", "중간", "빠름", "느림 (1–2년)",
                     "느림 (1–2년)", "중간 (6개월)")
    ) %>% datatable(rownames = FALSE, options = list(pageLength = 10, dom = "t"))
  })

  # ── TAB 8: Population ─────────────────────────────────────────────────
  pop_data <- eventReactive(input$run_pop, {
    n <- input$n_patients
    cv_e  <- input$cv_enzyme / 100
    cv_gfr <- input$cv_eGFR / 100
    set.seed(42)

    purrr::map_dfr(seq_len(n), function(i) {
      this_eGFR <- rnorm(1, input$eGFR0, input$eGFR0 * cv_gfr) %>% pmax(20)
      simulate_fabry(
        phenotype   = input$phenotype,
        is_amenable = TRUE,
        treatment   = input$pop_treatment,
        sim_years   = input$sim_years,
        eGFR0       = this_eGFR
      ) %>% mutate(PatientID = i)
    })
  })

  output$plot_pop_eGFR <- renderPlot({
    req(pop_data())
    d <- pop_data()
    ggplot(d, aes(Year, eGFR, group = PatientID)) +
      geom_line(alpha = 0.25, color = "#1565C0") +
      stat_summary(aes(group = 1), fun = median, geom = "line",
                   color = "navy", size = 1.5, linetype = "solid") +
      geom_hline(yintercept = 60, linetype = "dashed", color = "red") +
      labs(title = "가상 환자 집단 eGFR 분포 (median = 진한 선)",
           x = "시간 (년)", y = "eGFR (mL/min/1.73m²)") +
      theme_bw(base_size = 13)
  })

  output$plot_pop_lgb3 <- renderPlot({
    req(pop_data())
    d <- pop_data()
    ggplot(d, aes(Year, LysoGb3, group = PatientID)) +
      geom_line(alpha = 0.25, color = "#C62828") +
      stat_summary(aes(group = 1), fun = median, geom = "line",
                   color = "darkred", size = 1.5) +
      geom_hline(yintercept = 2.0, linetype = "dashed") +
      labs(title = "가상 환자 집단 Lyso-Gb3 분포",
           x = "시간 (년)", y = "Lyso-Gb3 (μg/L)") +
      theme_bw(base_size = 13)
  })

  output$plot_responder <- renderPlot({
    req(pop_data())
    d <- pop_data()
    resp <- d %>%
      filter(abs(Year - max(Year)) < 0.2) %>%
      group_by(PatientID) %>%
      summarise(LysoGb3_end = mean(LysoGb3)) %>%
      mutate(Response = ifelse(LysoGb3_end < 5.0, "반응자 (<5 μg/L)", "비반응자"))
    ggplot(resp, aes(Response, fill = Response)) +
      geom_bar() +
      scale_fill_manual(values = c("반응자 (<5 μg/L)" = "#2E7D32", "비반응자" = "#C62828")) +
      labs(title = "치료 종료 시점 반응자 비율 (Lyso-Gb3 <5 μg/L 기준)",
           x = "", y = "환자 수") +
      theme_bw(base_size = 13) + theme(legend.position = "none")
  })
}

##############################################################################
# RUN APP
##############################################################################
shinyApp(ui, server)
