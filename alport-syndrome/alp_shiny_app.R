## =============================================================================
## Alport Syndrome (AS) QSP Shiny App — skeleton
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## Companion model: alp_mrgsolve_model.R
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread("alp_mrgsolve_model.R")

SCENARIOS <- c(
  "1. 자연경과 - XLAS 남성"                          = "natural_xlas_m",
  "2. 자연경과 - ARAS (상염색체 열성)"                = "natural_aras",
  "3. 자연경과 - ADAS/이형접합 thin-BM"               = "natural_adas",
  "4. Ramipril (ACEi), 단백뇨 발생 후 시작"           = "ramipril_late",
  "5. Ramipril 조기(무증상기) 시작 (EARLY PRO-TECT)"  = "ramipril_early",
  "6. Losartan (ARB) 단독"                            = "losartan",
  "7. Sparsentan (dual ETA/AT1) 단독"                 = "sparsentan",
  "8. Bardoxolone methyl (Nrf2 활성화제)"             = "bardoxolone",
  "9. Lademirsen (RG-012, anti-miR-21 ASO)"           = "lademirsen",
  "10. 병용: RAAS 최대차단 + 다파글리플로진 + Sparsentan" = "combo_max"
)

GENOTYPES <- c(
  "ADAS - 이형접합 thin-BM (경증)"           = "0.4",
  "XLAS 여성 - 모자이크 (경증-중등도)"        = "0.55",
  "XLAS 남성 / ARAS (표준 중증도)"           = "1.0",
  "ARAS - truncating variant (최중증)"       = "1.15"
)

ui <- navbarPage(
  title = "Alport Syndrome (AS) QSP Explorer",

  ## ---- Tab 1: 환자 프로파일 ----
  tabPanel("환자 프로파일",
    sidebarLayout(
      sidebarPanel(
        h4("환자/유전형 파라미터"),
        selectInput("genotype", "유전형 (중증도)", choices = GENOTYPES, selected = "1.0"),
        numericInput("age_yr", "연령 (세)", value = 10, min = 0, max = 60),
        numericInput("horizon_yr", "시뮬레이션 기간 (년)", value = 25, min = 1, max = 50),
        hr(),
        actionButton("simulate", "시뮬레이션 실행", class = "btn-primary")
      ),
      mainPanel(
        h4("환자 요약"),
        verbatimTextOutput("patient_summary"),
        plotOutput("baseline_gbm_plot")
      )
    )
  ),

  ## ---- Tab 2: PK (약물 농도) ----
  tabPanel("PK (약동학)",
    sidebarLayout(
      sidebarPanel(
        selectInput("pk_drug", "약물 선택",
          choices = c("Ramipril/Ramiprilat", "Losartan/E-3174", "Sparsentan",
                      "Bardoxolone methyl", "Lademirsen", "Dapagliflozin", "Finerenone")),
        numericInput("dose", "용량 (mg 또는 상대단위)", value = 5),
        numericInput("interval_h", "투여 간격 (h)", value = 24)
      ),
      mainPanel(plotOutput("pk_plot"))
    )
  ),

  ## ---- Tab 3: PD 주요지표 (신장) ----
  tabPanel("PD 주요지표 (신장)",
    sidebarLayout(
      sidebarPanel(
        selectInput("scenario_pd", "치료 시나리오", choices = SCENARIOS)
      ),
      mainPanel(
        plotOutput("gbm_podocyte_plot"),
        plotOutput("fibrosis_hemodynamics_plot")
      )
    )
  ),

  ## ---- Tab 4: 임상 엔드포인트 ----
  tabPanel("임상 엔드포인트",
    sidebarLayout(
      sidebarPanel(
        selectInput("scenario_clin", "치료 시나리오", choices = SCENARIOS)
      ),
      mainPanel(
        plotOutput("egfr_trajectory_plot"),
        plotOutput("uacr_plot"),
        verbatimTextOutput("esrd_time_summary")
      )
    )
  ),

  ## ---- Tab 5: 시나리오 비교 ----
  tabPanel("시나리오 비교",
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("compare_scenarios", "비교할 시나리오 (최대 5개)",
          choices = SCENARIOS, selected = c("natural_xlas_m", "ramipril_late", "ramipril_early", "combo_max"))
      ),
      mainPanel(
        plotOutput("compare_egfr_plot"),
        DTOutput("compare_table")
      )
    )
  ),

  ## ---- Tab 6: 바이오마커 ----
  tabPanel("바이오마커",
    sidebarLayout(
      sidebarPanel(
        selectInput("biomarker", "바이오마커",
          choices = c("UACR", "GBM 구조 지표", "miR-21 활성도", "섬유화(IFTA) 지표"))
      ),
      mainPanel(plotOutput("biomarker_plot"))
    )
  ),

  ## ---- Tab 7: 이신장외 병변 (와우·안구) ----
  tabPanel("이신장외 병변 (와우·안구)",
    sidebarLayout(
      sidebarPanel(
        selectInput("scenario_extra", "치료 시나리오", choices = SCENARIOS)
      ),
      mainPanel(
        plotOutput("hearing_loss_plot"),
        plotOutput("ocular_score_plot")
      )
    )
  ),

  ## ---- Tab 8: 참고문헌 ----
  tabPanel("참고문헌",
    mainPanel(
      includeMarkdown("alp_references.md")
    )
  )
)

server <- function(input, output, session) {

  sim_data <- eventReactive(input$simulate, {
    mod %>%
      param(SEVERITY = as.numeric(input$genotype)) %>%
      mrgsim(end = 24 * 365 * input$horizon_yr, delta = 24 * 7) %>%
      as_tibble()
  })

  output$patient_summary <- renderPrint({
    cat("Genotype severity:", input$genotype, "\n")
    cat("Age (yr):", input$age_yr, "\n")
    cat("Simulation horizon (yr):", input$horizon_yr, "\n")
  })

  output$baseline_gbm_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, GBM_INTEG)) + geom_line(color = "#2874A6") +
      labs(x = "Years", y = "GBM structural integrity", title = "GBM Integrity — Natural History")
  })

  output$pk_plot <- renderPlot({
    e <- ev(amt = input$dose, cmt = "RAM_GUT", ii = input$interval_h, addl = 30)
    out <- mod %>% ev(e) %>% mrgsim(end = 24 * 40, delta = 1) %>% as_tibble()
    ggplot(out, aes(time, RAM_CENT)) + geom_line(color = "#B03A2E") +
      labs(x = "Time (h)", y = "Plasma amount", title = paste(input$pk_drug, "PK profile"))
  })

  output$gbm_podocyte_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS)) +
      geom_line(aes(y = GBM_INTEG, color = "GBM integrity")) +
      geom_line(aes(y = PODO_FRAC, color = "Podocyte fraction")) +
      labs(x = "Years", y = "Fraction", color = "", title = "GBM & Podocyte Trajectory")
  })

  output$fibrosis_hemodynamics_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, FIBROSIS)) + geom_line(color = "#943126") +
      labs(x = "Years", y = "Fibrosis (IFTA) index", title = "Fibrotic Progression")
  })

  output$egfr_trajectory_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, EGFR)) + geom_line(color = "#212F3C") +
      geom_hline(yintercept = 15, linetype = "dashed", color = "red") +
      labs(x = "Years", y = "eGFR (mL/min/1.73m2)", title = "eGFR Trajectory (ESRD threshold = 15)")
  })

  output$uacr_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, UACR)) + geom_line(color = "#6E2C00") +
      labs(x = "Years", y = "UACR (mg/g)", title = "Proteinuria Progression")
  })

  output$esrd_time_summary <- renderPrint({
    df <- sim_data()
    esrd_row <- df[df$EGFR <= 15, ]
    if (nrow(esrd_row) > 0) {
      cat("Time to ESRD (eGFR<=15):", round(esrd_row$FX_YEARS[1], 1), "years\n")
    } else {
      cat("ESRD not reached within simulation horizon\n")
    }
  })

  output$compare_egfr_plot <- renderPlot({
    plot.new()
    title("Scenario comparison placeholder — loop mod %>% param(...) %>% ev(...) per scenario")
  })

  output$compare_table <- renderDT({
    datatable(data.frame(Scenario = names(SCENARIOS), Time_to_ESRD_yr = NA))
  })

  output$biomarker_plot <- renderPlot({
    df <- sim_data()
    col <- switch(input$biomarker,
      "UACR" = "UACR", "GBM 구조 지표" = "GBM_INTEG",
      "miR-21 활성도" = "MIR21", "섬유화(IFTA) 지표" = "FIBROSIS")
    ggplot(df, aes(FX_YEARS, .data[[col]])) + geom_line(color = "#1B4F72") +
      labs(x = "Years", y = col, title = paste(input$biomarker, "Trajectory"))
  })

  output$hearing_loss_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, HEARING_LOSS)) + geom_line(color = "#117A65") +
      labs(x = "Years", y = "Hearing threshold shift (dB HL)", title = "Progressive Sensorineural Hearing Loss")
  })

  output$ocular_score_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_YEARS, OCULAR_SCORE)) + geom_line(color = "#CA6F1E") +
      labs(x = "Years", y = "Ocular severity index", title = "Lenticonus / Retinopathy Progression")
  })
}

shinyApp(ui, server)
