## =============================================================================
## Trigeminal Neuralgia (TN) QSP Shiny App — skeleton
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT, plotly (optional)
## Companion model: tn_mrgsolve_model.R
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread(here_or_local <- "tn_mrgsolve_model.R")

SCENARIOS <- c(
  "1. Untreated natural history"                 = "untreated",
  "2. Carbamazepine monotherapy (titrated)"       = "cbz_mono",
  "3. Oxcarbazepine monotherapy"                  = "oxc_mono",
  "4. CBZ + Baclofen combination"                 = "cbz_baclofen",
  "5. MVD (microvascular decompression), day 14"  = "mvd_post",
  "6. CBZ-intolerant -> Gabapentin + Pregabalin"  = "gbp_pgb_switch",
  "7. Percutaneous RF rhizotomy, day 30"          = "rf_rhizotomy"
)

ui <- navbarPage(
  title = "Trigeminal Neuralgia (TN) QSP Explorer",

  ## ---- Tab 1: Patient Profile ----
  tabPanel("환자 프로파일",
    sidebarLayout(
      sidebarPanel(
        h4("환자/질환 파라미터"),
        sliderInput("nvc_severity", "신경혈관압박(NVC) 중증도 (0-2)", 0, 2, 1.2, step = 0.1),
        selectInput("branch", "침범 분지", choices = c("V2 상악", "V3 하악", "V2+V3 복합", "V1 안신경(드묾)")),
        checkboxInput("secondary_ms", "이차성 TN (다발성경화증 동반)", FALSE),
        hr(),
        actionButton("simulate", "시뮬레이션 실행", class = "btn-primary")
      ),
      mainPanel(
        h4("환자 요약"),
        verbatimTextOutput("patient_summary"),
        plotOutput("baseline_severity_plot")
      )
    )
  ),

  ## ---- Tab 2: PK ----
  tabPanel("PK (약동학)",
    sidebarLayout(
      sidebarPanel(
        selectInput("pk_scenario", "시나리오", choices = SCENARIOS),
        checkboxGroupInput("pk_drugs", "표시 약물",
          choices = c("Carbamazepine" = "CBZ_conc", "CBZ-epoxide" = "CBZ_epox_c",
                      "Oxcarbazepine-MHD" = "OXC_mhd_c", "Baclofen" = "BAC_conc",
                      "Gabapentin" = "GBP_conc", "Pregabalin" = "PGB_conc"),
          selected = c("CBZ_conc", "CBZ_epox_c"))
      ),
      mainPanel(plotOutput("pk_plot", height = 450), DTOutput("pk_table"))
    )
  ),

  ## ---- Tab 3: 경로 PD (Nav channelopathy / central sensitization) ----
  tabPanel("경로 PD",
    sidebarLayout(
      sidebarPanel(selectInput("pd_scenario", "시나리오", choices = SCENARIOS)),
      mainPanel(
        plotOutput("navblock_plot", height = 300),
        plotOutput("centsens_plot", height = 300)
      )
    )
  ),

  ## ---- Tab 4: 임상 엔드포인트 ----
  tabPanel("임상 엔드포인트",
    sidebarLayout(
      sidebarPanel(selectInput("clin_scenario", "시나리오", choices = SCENARIOS)),
      mainPanel(
        plotOutput("paroxysm_plot", height = 300),
        plotOutput("pain_nrs_plot", height = 300),
        verbatimTextOutput("bni_summary")
      )
    )
  ),

  ## ---- Tab 5: 시나리오 비교 ----
  tabPanel("시나리오 비교",
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("compare_scenarios", "비교할 시나리오",
          choices = SCENARIOS, selected = c("untreated", "cbz_mono", "mvd_post"))
      ),
      mainPanel(
        plotOutput("compare_paroxysm_plot", height = 350),
        plotOutput("compare_pain_plot", height = 350),
        DTOutput("compare_summary_table")
      )
    )
  ),

  ## ---- Tab 6: 바이오마커/영상 ----
  tabPanel("바이오마커",
    mainPanel(
      h4("MRI FIESTA/CISS 신경혈관압박 소견 (개념적 지표)"),
      plotOutput("nav_upreg_plot", height = 300),
      p("실제 앱에서는 압박 중증도와 Nav 채널 발현 지수를 영상 바이오마커와 연계 표시.")
    )
  ),

  ## ---- Tab 7: 안전성 ----
  tabPanel("안전성",
    sidebarLayout(
      sidebarPanel(selectInput("safety_scenario", "시나리오", choices = SCENARIOS)),
      mainPanel(
        plotOutput("sodium_plot", height = 300),
        plotOutput("sedation_plot", height = 300),
        helpText("저나트륨혈증 경고 기준: Na+ < 135 mEq/L (OXC > CBZ 위험).")
      )
    )
  ),

  ## ---- Tab 8: 참고문헌 ----
  tabPanel("참고문헌",
    mainPanel(
      h4("주요 문헌"),
      includeMarkdown("tn_references.md")
    )
  )
)

server <- function(input, output, session) {

  run_sim <- function(scenario_key, end_time = 4320) {
    mod2 <- mod %>% param(NVC_SEVERITY = input$nvc_severity)

    dosing <- switch(scenario_key,
      untreated      = ev(time = 0, amt = 0, cmt = "CBZ_GUT"),
      cbz_mono       = ev(time = seq(0, end_time, 12), amt = 200, cmt = "CBZ_GUT", ii = 12, addl = 359),
      oxc_mono       = ev(time = seq(0, end_time, 12), amt = 300, cmt = "OXC_GUT", ii = 12, addl = 359),
      cbz_baclofen   = seq(ev(time = seq(0, end_time, 12), amt = 200, cmt = "CBZ_GUT", ii = 12, addl = 359),
                            ev(time = seq(0, end_time, 8), amt = 10, cmt = "BAC_GUT", ii = 8, addl = 539)),
      mvd_post       = ev(time = seq(0, 312, 12), amt = 200, cmt = "CBZ_GUT", ii = 12, addl = 25),
      gbp_pgb_switch = seq(ev(time = seq(0, end_time, 8), amt = 300, cmt = "GBP_GUT", ii = 8, addl = 539),
                            ev(time = seq(0, end_time, 12), amt = 75, cmt = "PGB_GUT", ii = 12, addl = 359)),
      rf_rhizotomy   = ev(time = 0, amt = 0, cmt = "CBZ_GUT")
    )

    if (scenario_key == "mvd_post") mod2 <- mod2 %>% param(MVD_ON = 1, MVD_TIME = 336)
    if (scenario_key == "rf_rhizotomy") mod2 <- mod2 %>% param(RF_ON = 1, RF_TIME = 720)

    mod2 %>% ev(dosing) %>% mrgsim(end = end_time) %>% as_tibble()
  }

  output$patient_summary <- renderText({
    paste0("NVC 중증도: ", input$nvc_severity, "\n분지: ", input$branch,
           "\n이차성(MS): ", ifelse(input$secondary_ms, "예", "아니오"))
  })

  output$baseline_severity_plot <- renderPlot({
    df <- run_sim("untreated")
    ggplot(df, aes(time/24, PAROX)) + geom_line(color = "firebrick") +
      labs(x = "일(day)", y = "발작 빈도 (회/일)", title = "무치료 자연경과") + theme_minimal()
  })

  output$pk_plot <- renderPlot({
    df <- run_sim(input$pk_scenario)
    df_long <- df %>% select(time, all_of(input$pk_drugs)) %>%
      pivot_longer(-time, names_to = "analyte", values_to = "conc")
    ggplot(df_long, aes(time/24, conc, color = analyte)) + geom_line() +
      labs(x = "일(day)", y = "농도", color = "") + theme_minimal()
  })
  output$pk_table <- renderDT({ run_sim(input$pk_scenario) %>% select(time, all_of(input$pk_drugs)) %>% head(50) })

  output$navblock_plot <- renderPlot({
    df <- run_sim(input$pd_scenario)
    ggplot(df, aes(time/24, navblock_tot)) + geom_line(color = "steelblue") +
      labs(x = "일", y = "Nav 채널 차단 분율", title = "약물 유도 Nav 차단") + theme_minimal()
  })
  output$centsens_plot <- renderPlot({
    df <- run_sim(input$pd_scenario)
    ggplot(df, aes(time/24, CENTSENS)) + geom_line(color = "purple") +
      labs(x = "일", y = "중추감작 지수") + theme_minimal()
  })

  output$paroxysm_plot <- renderPlot({
    df <- run_sim(input$clin_scenario)
    ggplot(df, aes(time/24, PAROX)) + geom_line(color = "firebrick") +
      labs(x = "일", y = "발작 빈도 (회/일)") + theme_minimal()
  })
  output$pain_nrs_plot <- renderPlot({
    df <- run_sim(input$clin_scenario)
    ggplot(df, aes(time/24, PAIN)) + geom_line(color = "darkorange") +
      labs(x = "일", y = "통증 NRS (0-10)") + theme_minimal()
  })
  output$bni_summary <- renderText({
    df <- run_sim(input$clin_scenario)
    last_pain <- tail(df$PAIN, 1)
    bni <- cut(last_pain, breaks = c(-1, 0.5, 3, 6, 8, 11), labels = c("I", "II", "III", "IV", "V"))
    paste0("최종 시점 BNI Pain Scale 근사: ", bni)
  })

  output$compare_paroxysm_plot <- renderPlot({
    dfs <- lapply(input$compare_scenarios, function(s) run_sim(s) %>% mutate(scenario = s))
    df_all <- bind_rows(dfs)
    ggplot(df_all, aes(time/24, PAROX, color = scenario)) + geom_line() +
      labs(x = "일", y = "발작 빈도 (회/일)", color = "시나리오") + theme_minimal()
  })
  output$compare_pain_plot <- renderPlot({
    dfs <- lapply(input$compare_scenarios, function(s) run_sim(s) %>% mutate(scenario = s))
    df_all <- bind_rows(dfs)
    ggplot(df_all, aes(time/24, PAIN, color = scenario)) + geom_line() +
      labs(x = "일", y = "통증 NRS", color = "시나리오") + theme_minimal()
  })
  output$compare_summary_table <- renderDT({
    dfs <- lapply(input$compare_scenarios, function(s) {
      d <- run_sim(s); tibble(scenario = s, final_paroxysm = tail(d$PAROX, 1),
                               final_pain = tail(d$PAIN, 1), final_Na = tail(d$NA_PLASMA, 1))
    })
    bind_rows(dfs)
  })

  output$nav_upreg_plot <- renderPlot({
    df <- run_sim("untreated")
    ggplot(df, aes(time/24, NAV_UPREG)) + geom_line(color = "brown") +
      labs(x = "일", y = "Nav 채널 발현/이소성흥분성 지수") + theme_minimal()
  })

  output$sodium_plot <- renderPlot({
    df <- run_sim(input$safety_scenario)
    ggplot(df, aes(time/24, NA_PLASMA)) + geom_line(color = "darkgreen") +
      geom_hline(yintercept = 135, linetype = "dashed", color = "red") +
      labs(x = "일", y = "혈장 Na+ (mEq/L)") + theme_minimal()
  })
  output$sedation_plot <- renderPlot({
    df <- run_sim(input$safety_scenario)
    ggplot(df, aes(time/24, SEDATION)) + geom_line(color = "gray30") +
      labs(x = "일", y = "졸림/운동실조 점수 (0-10)") + theme_minimal()
  })
}

shinyApp(ui, server)
