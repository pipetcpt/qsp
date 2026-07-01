## =============================================================================
## Osteogenesis Imperfecta (OI) QSP Shiny App — skeleton
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## Companion model: oi_mrgsolve_model.R
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread("oi_mrgsolve_model.R")

SCENARIOS <- c(
  "1. 자연경과 (미치료, type III/IV)"          = "untreated",
  "2. IV Pamidronate 주기요법 (Glorieux 1998)" = "pamidronate",
  "3. IV Zoledronic acid 연 2회"                = "zoledronic",
  "4. Denosumab SC (severe recessive OI)"       = "denosumab",
  "5. Teriparatide SC (성인 type I 한정)"       = "teriparatide",
  "6. Setrusumab (항-스클레로스틴, 투자적)"     = "setrusumab",
  "7. Fresolimumab (항-TGF-b, 투자적)"          = "fresolimumab"
)

ui <- navbarPage(
  title = "Osteogenesis Imperfecta (OI) QSP Explorer",

  ## ---- Tab 1: 환자 프로파일 ----
  tabPanel("환자 프로파일",
    sidebarLayout(
      sidebarPanel(
        h4("환자/질환 파라미터"),
        selectInput("oi_type", "Sillence 분류",
          choices = c("Type I (경증)" = "0.4", "Type III/IV (중등도-중증)" = "1.0",
                      "Type V-VIII (희귀 열성형)" = "1.4")),
        numericInput("age_yr", "연령 (세)", value = 8, min = 0, max = 80),
        numericInput("weight_kg", "체중 (kg)", value = 25, min = 3, max = 100),
        checkboxInput("growth_plate_open", "성장판 개방 (소아)", TRUE),
        hr(),
        actionButton("simulate", "시뮬레이션 실행", class = "btn-primary")
      ),
      mainPanel(
        h4("환자 요약"),
        verbatimTextOutput("patient_summary"),
        plotOutput("baseline_bmd_plot")
      )
    )
  ),

  ## ---- Tab 2: PK ----
  tabPanel("PK (약동학)",
    sidebarLayout(
      sidebarPanel(
        selectInput("pk_scenario", "시나리오", choices = SCENARIOS),
        checkboxGroupInput("pk_drugs", "표시 약물",
          choices = c("Pamidronate" = "PAM_Cp", "Zoledronic acid" = "ZOL_Cp",
                      "Denosumab" = "DMAB_Cp", "Teriparatide" = "TPTD_Cp",
                      "Setrusumab" = "SETRU_Cp", "Fresolimumab" = "FRESO_Cp"),
          selected = c("PAM_Cp", "ZOL_Cp"))
      ),
      mainPanel(plotOutput("pk_plot", height = 450), DTOutput("pk_table"))
    )
  ),

  ## ---- Tab 3: PD 주요지표 (골재형성 경로) ----
  tabPanel("PD 주요지표",
    sidebarLayout(
      sidebarPanel(selectInput("pd_scenario", "시나리오", choices = SCENARIOS)),
      mainPanel(
        plotOutput("sost_rankl_plot", height = 300),
        plotOutput("ob_oc_plot", height = 300)
      )
    )
  ),

  ## ---- Tab 4: 임상 엔드포인트 ----
  tabPanel("임상 엔드포인트",
    sidebarLayout(
      sidebarPanel(selectInput("clin_scenario", "시나리오", choices = SCENARIOS)),
      mainPanel(
        plotOutput("bmd_plot", height = 300),
        plotOutput("fracture_plot", height = 300),
        verbatimTextOutput("clin_summary")
      )
    )
  ),

  ## ---- Tab 5: 시나리오 비교 ----
  tabPanel("시나리오 비교",
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("compare_scenarios", "비교할 시나리오",
          choices = SCENARIOS, selected = c("untreated", "pamidronate", "zoledronic", "setrusumab"))
      ),
      mainPanel(
        plotOutput("compare_bmd_plot", height = 350),
        plotOutput("compare_fracture_plot", height = 350),
        DTOutput("compare_summary_table")
      )
    )
  ),

  ## ---- Tab 6: 바이오마커 (골대사표지자) ----
  tabPanel("바이오마커",
    sidebarLayout(
      sidebarPanel(selectInput("bio_scenario", "시나리오", choices = SCENARIOS)),
      mainPanel(
        plotOutput("p1np_ctx_plot", height = 350),
        helpText("P1NP: 골형성 표지자, CTX: 골흡수 표지자. OI는 전형적으로 고회전율(high-turnover) 양상.")
      )
    )
  ),

  ## ---- Tab 7: 성장/골격 ----
  tabPanel("성장/골격",
    sidebarLayout(
      sidebarPanel(selectInput("growth_scenario", "시나리오", choices = SCENARIOS)),
      mainPanel(
        plotOutput("height_plot", height = 300),
        helpText("Teriparatide는 성장판 개방 소아에서는 골육종 위험으로 금기(투여 시 경고 표시).")
      )
    )
  ),

  ## ---- Tab 8: 참고문헌 ----
  tabPanel("참고문헌",
    mainPanel(
      h4("주요 문헌"),
      includeMarkdown("oi_references.md")
    )
  )
)

server <- function(input, output, session) {

  run_sim <- function(scenario_key, end_time = 24*365*4) {
    mod2 <- mod %>% param(SEVERITY = as.numeric(input$oi_type))

    dosing <- switch(scenario_key,
      untreated    = ev(time = 0, amt = 0, cmt = "PAM_CENT"),
      pamidronate  = ev(time = seq(0, end_time, 24*90), amt = input$weight_kg*1,
                         cmt = "PAM_CENT", rate = input$weight_kg*1/2, ii = 24*90, addl = 15),
      zoledronic   = ev(time = seq(0, end_time, 24*182), amt = input$weight_kg*0.05,
                         cmt = "ZOL_CENT", rate = input$weight_kg*0.05/1, ii = 24*182, addl = 7),
      denosumab    = ev(time = seq(0, end_time, 24*90), amt = input$weight_kg*1,
                         cmt = "DMAB_DEPOT", ii = 24*90, addl = 15),
      teriparatide = ev(time = seq(0, end_time, 24), amt = 20, cmt = "TPTD_DEPOT",
                         ii = 24, addl = 1459),
      setrusumab   = ev(time = seq(0, end_time, 24*28), amt = input$weight_kg*20,
                         cmt = "SETRU_CENT", rate = input$weight_kg*20/1, ii = 24*28, addl = 51),
      fresolimumab = ev(time = seq(0, end_time, 24*28), amt = input$weight_kg*1,
                         cmt = "FRESO_CENT", rate = input$weight_kg*1/1, ii = 24*28, addl = 51)
    )

    mod2 %>% ev(dosing) %>% mrgsim(end = end_time, delta = 24) %>% as_tibble()
  }

  output$patient_summary <- renderText({
    paste0("Sillence 분류 중증도: ", input$oi_type, "\n연령: ", input$age_yr, "세\n",
           "체중: ", input$weight_kg, "kg\n성장판 개방: ", ifelse(input$growth_plate_open, "예", "아니오"))
  })

  output$baseline_bmd_plot <- renderPlot({
    df <- run_sim("untreated")
    ggplot(df, aes(time/365, BMC)) + geom_line(color = "firebrick") +
      labs(x = "연(year)", y = "BMC 지수 (기저 대비 변화)", title = "무치료 자연경과") + theme_minimal()
  })

  output$pk_plot <- renderPlot({
    df <- run_sim(input$pk_scenario)
    df_long <- df %>% select(time, all_of(input$pk_drugs)) %>%
      pivot_longer(-time, names_to = "analyte", values_to = "conc")
    ggplot(df_long, aes(time/24, conc, color = analyte)) + geom_line() +
      labs(x = "일(day)", y = "농도 (mg/L 또는 ug/L)", color = "") + theme_minimal()
  })
  output$pk_table <- renderDT({ run_sim(input$pk_scenario) %>% select(time, all_of(input$pk_drugs)) %>% head(50) })

  output$sost_rankl_plot <- renderPlot({
    df <- run_sim(input$pd_scenario)
    ggplot(df, aes(time/365)) +
      geom_line(aes(y = SOST_eff, color = "Sclerostin (유효)")) +
      geom_line(aes(y = RANKL_free, color = "RANKL (유리)")) +
      labs(x = "연", y = "상대 수준", color = "") + theme_minimal()
  })
  output$ob_oc_plot <- renderPlot({
    df <- run_sim(input$pd_scenario)
    ggplot(df, aes(time/365)) +
      geom_line(aes(y = OB, color = "조골세포 활성")) +
      geom_line(aes(y = OC, color = "파골세포 활성")) +
      labs(x = "연", y = "활성 지수", color = "") + theme_minimal()
  })

  output$bmd_plot <- renderPlot({
    df <- run_sim(input$clin_scenario)
    ggplot(df, aes(time/365, BMC)) + geom_line(color = "steelblue") +
      labs(x = "연", y = "BMC 지수 (Z-score 변화)") + theme_minimal()
  })
  output$fracture_plot <- renderPlot({
    df <- run_sim(input$clin_scenario)
    ggplot(df, aes(time/365, FX_CUM)) + geom_line(color = "darkorange") +
      labs(x = "연", y = "누적 골절 수") + theme_minimal()
  })
  output$clin_summary <- renderText({
    df <- run_sim(input$clin_scenario)
    paste0("4년 시점 BMC 변화: ", round(tail(df$BMC, 1), 2),
           "\n4년 시점 누적 골절: ", round(tail(df$FX_CUM, 1), 2))
  })

  output$compare_bmd_plot <- renderPlot({
    dfs <- lapply(input$compare_scenarios, function(s) run_sim(s) %>% mutate(scenario = s))
    df_all <- bind_rows(dfs)
    ggplot(df_all, aes(time/365, BMC, color = scenario)) + geom_line() +
      labs(x = "연", y = "BMC 지수", color = "시나리오") + theme_minimal()
  })
  output$compare_fracture_plot <- renderPlot({
    dfs <- lapply(input$compare_scenarios, function(s) run_sim(s) %>% mutate(scenario = s))
    df_all <- bind_rows(dfs)
    ggplot(df_all, aes(time/365, FX_CUM, color = scenario)) + geom_line() +
      labs(x = "연", y = "누적 골절 수", color = "시나리오") + theme_minimal()
  })
  output$compare_summary_table <- renderDT({
    dfs <- lapply(input$compare_scenarios, function(s) {
      d <- run_sim(s); tibble(scenario = s, final_BMC = tail(d$BMC, 1), final_fracture = tail(d$FX_CUM, 1))
    })
    bind_rows(dfs)
  })

  output$p1np_ctx_plot <- renderPlot({
    df <- run_sim(input$bio_scenario)
    ggplot(df, aes(time/365)) +
      geom_line(aes(y = P1NP, color = "P1NP (형성)")) +
      geom_line(aes(y = CTX*50, color = "CTX x50 (흡수)")) +
      labs(x = "연", y = "표지자 수준", color = "") + theme_minimal()
  })

  output$height_plot <- renderPlot({
    df <- run_sim(input$growth_scenario)
    ggplot(df, aes(time/365, HEIGHT_Z)) + geom_line(color = "darkgreen") +
      labs(x = "연", y = "신장 Z-score") + theme_minimal()
  })
}

shinyApp(ui, server)
