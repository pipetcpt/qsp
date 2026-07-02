## =============================================================================
## PTSD QSP Shiny App — skeleton
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## Companion model: ptsd_mrgsolve_model.R
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread("ptsd_mrgsolve_model.R")

SCENARIOS <- c(
  "1. 자연경과 - 중등도 트라우마 (미치료)"                       = "natural_moderate",
  "2. 자연경과 - 중증/반복 트라우마 + FKBP5 + 해리형 (최악예후)"  = "natural_severe_dissoc",
  "3. Sertraline 100-200mg/day"                                  = "sertraline",
  "4. Paroxetine 20-50mg/day"                                    = "paroxetine",
  "5. Prazosin 1-15mg qhs (악몽 부가요법)"                        = "prazosin_addon",
  "6. 외상중심 심리치료 주1회 x12 (PE/CPT/EMDR)"                  = "trauma_psychotherapy",
  "7. Ketamine 0.5mg/kg IV x6 (2주간)"                            = "ketamine",
  "8. MDMA-보조 심리치료 (3세션 + 12회 준비/통합치료)"             = "mdma_therapy",
  "9. 고회복탄력성 + 경증 트라우마 (자연관해)"                     = "high_resilience",
  "10. 병용: SSRI + 주간 심리치료 + Prazosin"                     = "combination_soc"
)

TRAUMA_TYPES <- c(
  "경증/단일 트라우마 (SEVERITY=0.6)" = "0.6",
  "표준/중등도 (전투, 폭력) (SEVERITY=1.0)" = "1.0",
  "중증/반복 (아동학대+재외상) (SEVERITY=1.5)" = "1.5"
)

ui <- navbarPage(
  title = "PTSD QSP Explorer",

  ## ---- Tab 1: 환자 프로파일 ----
  tabPanel("환자 프로파일",
    sidebarLayout(
      sidebarPanel(
        h4("환자/위험요인 파라미터"),
        selectInput("trauma_severity", "트라우마 중증도", choices = TRAUMA_TYPES, selected = "1.0"),
        checkboxInput("fkbp5_risk", "FKBP5 위험 대립유전자 x 아동기 트라우마", value = FALSE),
        checkboxInput("dissociation", "해리 하위유형(주변외상성 해리)", value = FALSE),
        numericInput("resilience", "회복탄력성 지수 (0.5=취약 ~ 1.5=회복력 높음)", value = 1.0, min = 0.5, max = 1.5, step = 0.1),
        numericInput("age_yr", "연령 (세)", value = 32, min = 18, max = 80),
        numericInput("horizon_wk", "시뮬레이션 기간 (주)", value = 52, min = 4, max = 260),
        hr(),
        actionButton("simulate", "시뮬레이션 실행", class = "btn-primary")
      ),
      mainPanel(
        h4("환자 요약"),
        verbatimTextOutput("patient_summary"),
        plotOutput("baseline_caps_plot")
      )
    )
  ),

  ## ---- Tab 2: PK (약동학) ----
  tabPanel("PK (약동학)",
    sidebarLayout(
      sidebarPanel(
        selectInput("pk_drug", "약물 선택",
          choices = c("Sertraline", "Paroxetine", "Prazosin", "Ketamine", "MDMA")),
        numericInput("dose", "용량 (mg)", value = 100),
        numericInput("interval_h", "투여 간격 (h)", value = 24)
      ),
      mainPanel(plotOutput("pk_plot"))
    )
  ),

  ## ---- Tab 3: PD 주요지표 (공포회로/HPA) ----
  tabPanel("PD 주요지표 (공포회로/HPA)",
    sidebarLayout(
      sidebarPanel(
        selectInput("scenario_pd", "치료 시나리오", choices = SCENARIOS)
      ),
      mainPanel(
        plotOutput("circuit_plot"),
        plotOutput("cortisol_ne_plot")
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
        plotOutput("caps5_trajectory_plot"),
        plotOutput("symptom_clusters_plot"),
        verbatimTextOutput("remission_summary")
      )
    )
  ),

  ## ---- Tab 5: 시나리오 비교 ----
  tabPanel("시나리오 비교",
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("compare_scenarios", "비교할 시나리오 (최대 5개)",
          choices = SCENARIOS, selected = c("natural_moderate", "sertraline", "trauma_psychotherapy", "mdma_therapy", "combination_soc"))
      ),
      mainPanel(
        plotOutput("compare_caps5_plot"),
        DTOutput("compare_table")
      )
    )
  ),

  ## ---- Tab 6: 바이오마커 ----
  tabPanel("바이오마커",
    sidebarLayout(
      sidebarPanel(
        selectInput("biomarker", "바이오마커",
          choices = c("코르티솔", "청반/NE 긴장도", "편도체 반응성 지수", "vmPFC 긴장도", "소거기억 강도", "수면장애 지수"))
      ),
      mainPanel(plotOutput("biomarker_plot"))
    )
  ),

  ## ---- Tab 7: 수면/야간증상 ----
  tabPanel("수면/야간증상",
    sidebarLayout(
      sidebarPanel(
        selectInput("scenario_sleep", "치료 시나리오", choices = SCENARIOS)
      ),
      mainPanel(
        plotOutput("sleep_plot")
      )
    )
  ),

  ## ---- Tab 8: 참고문헌 ----
  tabPanel("참고문헌",
    mainPanel(
      includeMarkdown("ptsd_references.md")
    )
  )
)

server <- function(input, output, session) {

  sim_data <- eventReactive(input$simulate, {
    mod %>%
      param(TRAUMA_SEVERITY = as.numeric(input$trauma_severity),
            FKBP5_RISK = as.numeric(input$fkbp5_risk),
            DISSOCIATION = as.numeric(input$dissociation),
            RESILIENCE = input$resilience) %>%
      mrgsim(end = 24 * 7 * input$horizon_wk, delta = 24) %>%
      as_tibble()
  })

  output$patient_summary <- renderPrint({
    cat("Trauma severity:", input$trauma_severity, "\n")
    cat("FKBP5 risk allele x childhood trauma:", input$fkbp5_risk, "\n")
    cat("Dissociative subtype:", input$dissociation, "\n")
    cat("Resilience index:", input$resilience, "\n")
    cat("Age (yr):", input$age_yr, "\n")
    cat("Simulation horizon (wk):", input$horizon_wk, "\n")
  })

  output$baseline_caps_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_WEEKS, CAPS5)) + geom_line(color = "#2874A6") +
      labs(x = "Weeks", y = "CAPS-5-like total score", title = "Symptom Severity — Natural History")
  })

  output$pk_plot <- renderPlot({
    cmt <- switch(input$pk_drug,
      "Sertraline" = "SERT_GUT", "Paroxetine" = "PAROX_GUT", "Prazosin" = "PRAZ_GUT",
      "Ketamine" = "KET_CENT", "MDMA" = "MDMA_GUT")
    e <- ev(amt = input$dose, cmt = cmt, ii = input$interval_h, addl = 30)
    out <- mod %>% ev(e) %>% mrgsim(end = 24 * 14, delta = 1) %>% as_tibble()
    ycol <- switch(input$pk_drug,
      "Sertraline" = "SERT_CENT", "Paroxetine" = "PAROX_CENT", "Prazosin" = "PRAZ_CENT",
      "Ketamine" = "KET_CENT", "MDMA" = "MDMA_CENT")
    ggplot(out, aes(time, .data[[ycol]])) + geom_line(color = "#B03A2E") +
      labs(x = "Time (h)", y = "Plasma amount", title = paste(input$pk_drug, "PK profile"))
  })

  output$circuit_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_WEEKS)) +
      geom_line(aes(y = AMYG_REACT, color = "Amygdala reactivity")) +
      geom_line(aes(y = VMPFC_TONE, color = "vmPFC inhibitory tone")) +
      geom_line(aes(y = EXT_MEM, color = "Extinction-memory strength")) +
      geom_line(aes(y = FEAR_MEM, color = "Fear-memory strength")) +
      labs(x = "Weeks", y = "Index (0-1)", color = "", title = "Fear-Extinction Circuit Trajectory")
  })

  output$cortisol_ne_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_WEEKS)) +
      geom_line(aes(y = CORTISOL, color = "Serum cortisol (ug/dL)")) +
      geom_line(aes(y = NE_TONE * 20, color = "NE/LC tone (x20 scaled)")) +
      labs(x = "Weeks", y = "Value", color = "", title = "HPA Axis & Noradrenergic Tone")
  })

  output$caps5_trajectory_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_WEEKS, CAPS5)) + geom_line(color = "#212F3C") +
      geom_hline(yintercept = 20, linetype = "dashed", color = "darkgreen") +
      labs(x = "Weeks", y = "CAPS-5-like total score", title = "CAPS-5 Trajectory (remission threshold <=20)")
  })

  output$symptom_clusters_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_WEEKS)) +
      geom_line(aes(y = INTRUSION, color = "B: Intrusion")) +
      geom_line(aes(y = AVOIDANCE, color = "C: Avoidance")) +
      geom_line(aes(y = NEGCOG, color = "D: Negative cognition/mood")) +
      geom_line(aes(y = HYPERAROUSE, color = "E: Hyperarousal")) +
      labs(x = "Weeks", y = "Cluster severity (0-40)", color = "", title = "DSM-5 Symptom Cluster Trajectories")
  })

  output$remission_summary <- renderPrint({
    df <- sim_data()
    rem_row <- df[df$CAPS5 <= 20, ]
    if (nrow(rem_row) > 0) {
      cat("Time to remission (CAPS-5<=20):", round(rem_row$FX_WEEKS[1], 1), "weeks\n")
    } else {
      cat("Remission not reached within simulation horizon\n")
    }
  })

  output$compare_caps5_plot <- renderPlot({
    plot.new()
    title("Scenario comparison placeholder — loop mod %>% param(...) %>% ev(...) per scenario")
  })

  output$compare_table <- renderDT({
    datatable(data.frame(Scenario = names(SCENARIOS), Time_to_Remission_wk = NA))
  })

  output$biomarker_plot <- renderPlot({
    df <- sim_data()
    col <- switch(input$biomarker,
      "코르티솔" = "CORTISOL", "청반/NE 긴장도" = "NE_TONE",
      "편도체 반응성 지수" = "AMYG_REACT", "vmPFC 긴장도" = "VMPFC_TONE",
      "소거기억 강도" = "EXT_MEM", "수면장애 지수" = "SLEEP_DIST")
    ggplot(df, aes(FX_WEEKS, .data[[col]])) + geom_line(color = "#1B4F72") +
      labs(x = "Weeks", y = col, title = paste(input$biomarker, "Trajectory"))
  })

  output$sleep_plot <- renderPlot({
    df <- sim_data()
    ggplot(df, aes(FX_WEEKS, SLEEP_DIST)) + geom_line(color = "#78281F") +
      labs(x = "Weeks", y = "Sleep disturbance index (0-1)", title = "Nightmare/Sleep Disturbance Trajectory")
  })
}

shinyApp(ui, server)
