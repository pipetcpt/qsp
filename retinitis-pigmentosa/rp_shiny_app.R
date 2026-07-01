## =============================================================================
## Retinitis Pigmentosa (RP) QSP Shiny App — skeleton
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## Companion model: rp_mrgsolve_model.R
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

mod <- mread("rp_mrgsolve_model.R")

SCENARIOS <- c(
  "1. 자연경과 - RHO-adRP"                       = "natural_rho",
  "2. 자연경과 - RPGR-XLRP"                      = "natural_xlrp",
  "3. 자연경과 - RPE65-LCA/EOSRD"                = "natural_rpe65",
  "4. Voretigene neparvovec (Luxturna, RPE65)"   = "voretigene",
  "5. 투자적 RPGR 유전자치료 (XLRP)"             = "rpgr_gt",
  "6. MCO-010 광유전학 치료 (말기, 유전형 무관)" = "mco010",
  "7. CNTF 이식형 서방출 임플란트"                = "cntf",
  "8. N-아세틸시스테인 (경구 항산화제)"          = "nac",
  "9. 비타민A palmitate + DHA (경구)"            = "vita",
  "10. Voretigene neparvovec + NAC 병용"         = "combo"
)

GENOTYPES <- c(
  "USH2A (경도-중등도, 서서히 진행)" = "0.75",
  "RHO-상염색체우성 RP (중등도)"     = "1.0",
  "RPGR-X연관 RP (빠른 진행)"        = "1.6",
  "RPE65-LCA/조기발현중증망막이영양증" = "2.5"
)

ui <- navbarPage(
  title = "Retinitis Pigmentosa (RP) QSP Explorer",

  ## ---- Tab 1: 환자 프로파일 ----
  tabPanel("환자 프로파일",
    sidebarLayout(
      sidebarPanel(
        h4("환자/유전형 파라미터"),
        selectInput("genotype", "유전형 (중증도)", choices = GENOTYPES, selected = "1.0"),
        numericInput("age_yr", "연령 (세)", value = 20, min = 0, max = 90),
        numericInput("horizon_yr", "시뮬레이션 기간 (년)", value = 15, min = 1, max = 40),
        hr(),
        actionButton("simulate", "시뮬레이션 실행", class = "btn-primary")
      ),
      mainPanel(
        h4("환자 요약"),
        verbatimTextOutput("patient_summary"),
        plotOutput("baseline_rodcone_plot")
      )
    )
  ),

  ## ---- Tab 2: PK (유전자벡터 발현 · 약물 농도) ----
  tabPanel("PK (약동학)",
    sidebarLayout(
      sidebarPanel(
        selectInput("pk_scenario", "시나리오", choices = SCENARIOS),
        checkboxGroupInput("pk_analytes", "표시 항목",
          choices = c("RPE65 발현(GT65_EXPR)" = "GT65_EXPR", "RPGR 발현(GTRPGR_EXPR)" = "GTRPGR_EXPR",
                      "MCO 옵신 발현(OP_EXPR)" = "OP_EXPR", "CNTF 조직농도" = "CNTF_DEV",
                      "NAC 혈장농도" = "NAC_CENT", "비타민A 저장고" = "VITA_STORE"),
          selected = c("GT65_EXPR", "CNTF_DEV"))
      ),
      mainPanel(plotOutput("pk_plot", height = 450), DTOutput("pk_table"))
    )
  ),

  ## ---- Tab 3: PD 주요지표 (광수용체 생존/산화스트레스) ----
  tabPanel("PD 주요지표",
    sidebarLayout(
      sidebarPanel(selectInput("pd_scenario", "시나리오", choices = SCENARIOS)),
      mainPanel(
        plotOutput("rodcone_plot", height = 300),
        plotOutput("ros_microglia_plot", height = 300)
      )
    )
  ),

  ## ---- Tab 4: 임상 엔드포인트 ----
  tabPanel("임상 엔드포인트",
    sidebarLayout(
      sidebarPanel(selectInput("clin_scenario", "시나리오", choices = SCENARIOS)),
      mainPanel(
        plotOutput("erg_plot", height = 280),
        plotOutput("vf_bcva_plot", height = 280),
        plotOutput("fst_mlmt_plot", height = 280),
        verbatimTextOutput("clin_summary")
      )
    )
  ),

  ## ---- Tab 5: 시나리오 비교 ----
  tabPanel("시나리오 비교",
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("compare_scenarios", "비교할 시나리오",
          choices = SCENARIOS, selected = c("natural_rpe65", "voretigene", "mco010", "cntf"))
      ),
      mainPanel(
        plotOutput("compare_mlmt_plot", height = 350),
        plotOutput("compare_fst_plot", height = 350),
        DTOutput("compare_summary_table")
      )
    )
  ),

  ## ---- Tab 6: 바이오마커 (RGC/황반부종) ----
  tabPanel("바이오마커",
    sidebarLayout(
      sidebarPanel(selectInput("bio_scenario", "시나리오", choices = SCENARIOS),
                   checkboxInput("cai_on", "CAI 치료 (도르졸라마이드/아세타졸아마이드) 병용", FALSE)),
      mainPanel(
        plotOutput("rgc_plot", height = 300),
        plotOutput("cme_plot", height = 300),
        helpText("RGC_FRAC: 망막신경절세포 생존율(광유전학/보철 치료 반응성 결정). CME_CST: 낭포황반부종 중심망막두께.")
      )
    )
  ),

  ## ---- Tab 7: 유전자치료/광유전학 반응 ----
  tabPanel("유전자치료 · 광유전학",
    sidebarLayout(
      sidebarPanel(selectInput("gt_scenario", "시나리오", choices = SCENARIOS,
                                selected = "voretigene")),
      mainPanel(
        plotOutput("visualcycle_plot", height = 300),
        plotOutput("bypass_plot", height = 300),
        helpText("Voretigene neparvovec: 시각회로(visual cycle) 회복 지표. MCO-010: 광수용체 생존과 무관한 우회(bypass) 광감작 신호(RGC 생존율에 비례).")
      )
    )
  ),

  ## ---- Tab 8: 참고문헌 ----
  tabPanel("참고문헌",
    mainPanel(
      h4("주요 문헌"),
      includeMarkdown("rp_references.md")
    )
  )
)

server <- function(input, output, session) {

  run_sim <- function(scenario_key, end_time = NULL) {
    if (is.null(end_time)) end_time <- 24*365*input$horizon_yr

    mod2 <- mod %>% param(SEVERITY = as.numeric(input$genotype))

    is_rpe65_geno <- as.numeric(input$genotype) == 2.5
    is_xlrp_geno  <- as.numeric(input$genotype) == 1.6

    mod2 <- switch(scenario_key,
      natural_rho    = mod2 %>% param(IS_RPE65 = 0, IS_XLRP = 0),
      natural_xlrp   = mod2 %>% param(IS_RPE65 = 0, IS_XLRP = 1),
      natural_rpe65  = mod2 %>% param(IS_RPE65 = 1, IS_XLRP = 0),
      voretigene     = mod2 %>% param(IS_RPE65 = 1, IS_XLRP = 0),
      rpgr_gt        = mod2 %>% param(IS_RPE65 = 0, IS_XLRP = 1),
      mco010         = mod2 %>% param(IS_RPE65 = as.numeric(is_rpe65_geno), IS_XLRP = as.numeric(is_xlrp_geno)),
      cntf           = mod2 %>% param(CNTF_PROD = 0.050),
      nac            = mod2,
      vita           = mod2 %>% param(VITA_KIN = 0.0020),
      combo          = mod2 %>% param(IS_RPE65 = 1, IS_XLRP = 0)
    )
    mod2 <- mod2 %>% param(CAI_ON = as.numeric(input$cai_on))

    dosing <- switch(scenario_key,
      natural_rho    = ev(time = 0, amt = 0, cmt = "GT65_VG"),
      natural_xlrp   = ev(time = 0, amt = 0, cmt = "GT65_VG"),
      natural_rpe65  = ev(time = 0, amt = 0, cmt = "GT65_VG"),
      voretigene     = ev(time = 0, amt = 100, cmt = "GT65_VG"),
      rpgr_gt        = ev(time = 0, amt = 100, cmt = "GTRPGR_VG"),
      mco010         = ev(time = 0, amt = 100, cmt = "OP_VG"),
      cntf           = ev(time = 0, amt = 0, cmt = "GT65_VG"),
      nac            = ev(time = seq(0, end_time, 8), amt = 600, cmt = "NAC_GUT", ii = 8, addl = length(seq(0, end_time, 8))),
      vita           = ev(time = 0, amt = 0, cmt = "GT65_VG"),
      combo          = ev(time = 0, amt = 100, cmt = "GT65_VG")
    )

    out <- mod2 %>% ev(dosing) %>% mrgsim(end = end_time, delta = 24) %>% as_tibble()
    out
  }

  output$patient_summary <- renderText({
    paste0("유전형/중증도 계수: ", input$genotype, "\n연령: ", input$age_yr, "세\n",
           "시뮬레이션 기간: ", input$horizon_yr, "년")
  })

  output$baseline_rodcone_plot <- renderPlot({
    df <- run_sim("natural_rho")
    ggplot(df, aes(time/(24*365))) +
      geom_line(aes(y = ROD_FRAC, color = "간상세포 생존율")) +
      geom_line(aes(y = CONE_FRAC, color = "원추세포 생존율")) +
      labs(x = "연(year)", y = "생존 분율", color = "", title = "자연경과 (미치료)") + theme_minimal()
  })

  output$pk_plot <- renderPlot({
    df <- run_sim(input$pk_scenario)
    df_long <- df %>% select(time, all_of(input$pk_analytes)) %>%
      pivot_longer(-time, names_to = "analyte", values_to = "value")
    ggplot(df_long, aes(time/(24*365), value, color = analyte)) + geom_line() +
      labs(x = "연(year)", y = "발현/농도 (상대 단위)", color = "") + theme_minimal()
  })
  output$pk_table <- renderDT({ run_sim(input$pk_scenario) %>% select(time, all_of(input$pk_analytes)) %>% head(50) })

  output$rodcone_plot <- renderPlot({
    df <- run_sim(input$pd_scenario)
    ggplot(df, aes(time/(24*365))) +
      geom_line(aes(y = ROD_FRAC, color = "간상세포 생존율")) +
      geom_line(aes(y = CONE_FRAC, color = "원추세포 생존율")) +
      labs(x = "연", y = "생존 분율", color = "") + theme_minimal()
  })
  output$ros_microglia_plot <- renderPlot({
    df <- run_sim(input$pd_scenario)
    ggplot(df, aes(time/(24*365))) +
      geom_line(aes(y = ROS, color = "산화스트레스(ROS)")) +
      geom_line(aes(y = MICROGLIA, color = "미세아교세포 활성")) +
      labs(x = "연", y = "상대 지수", color = "") + theme_minimal()
  })

  output$erg_plot <- renderPlot({
    df <- run_sim(input$clin_scenario)
    ggplot(df, aes(time/(24*365))) +
      geom_line(aes(y = ERG_ROD, color = "ERG 간상세포 b파")) +
      geom_line(aes(y = ERG_CONE, color = "ERG 원추세포 진폭")) +
      labs(x = "연", y = "진폭 (uV)", color = "") + theme_minimal()
  })
  output$vf_bcva_plot <- renderPlot({
    df <- run_sim(input$clin_scenario)
    ggplot(df, aes(time/(24*365))) +
      geom_line(aes(y = VF_AREA/10, color = "시야 면적/10 (deg2)")) +
      geom_line(aes(y = BCVA*500, color = "BCVA x500 (logMAR)")) +
      labs(x = "연", y = "값 (스케일 조정)", color = "") + theme_minimal()
  })
  output$fst_mlmt_plot <- renderPlot({
    df <- run_sim(input$clin_scenario)
    ggplot(df, aes(time/(24*365))) +
      geom_line(aes(y = FST, color = "FST (dB, 낮을수록 좋음)")) +
      geom_line(aes(y = MLMT*10, color = "MLMT x10")) +
      labs(x = "연", y = "값", color = "") + theme_minimal()
  })
  output$clin_summary <- renderText({
    df <- run_sim(input$clin_scenario)
    paste0(input$horizon_yr, "년 시점 ERG 간상세포 진폭: ", round(tail(df$ERG_ROD, 1), 1), " uV\n",
           input$horizon_yr, "년 시점 시야 면적: ", round(tail(df$VF_AREA, 1), 0), " deg2\n",
           input$horizon_yr, "년 시점 BCVA: ", round(tail(df$BCVA, 1), 2), " logMAR\n",
           input$horizon_yr, "년 시점 MLMT: ", round(tail(df$MLMT, 1), 2))
  })

  output$compare_mlmt_plot <- renderPlot({
    dfs <- lapply(input$compare_scenarios, function(s) run_sim(s) %>% mutate(scenario = s))
    df_all <- bind_rows(dfs)
    ggplot(df_all, aes(time/(24*365), MLMT, color = scenario)) + geom_line() +
      labs(x = "연", y = "MLMT score", color = "시나리오") + theme_minimal()
  })
  output$compare_fst_plot <- renderPlot({
    dfs <- lapply(input$compare_scenarios, function(s) run_sim(s) %>% mutate(scenario = s))
    df_all <- bind_rows(dfs)
    ggplot(df_all, aes(time/(24*365), FST, color = scenario)) + geom_line() +
      labs(x = "연", y = "FST (dB)", color = "시나리오") + theme_minimal()
  })
  output$compare_summary_table <- renderDT({
    dfs <- lapply(input$compare_scenarios, function(s) {
      d <- run_sim(s)
      tibble(scenario = s, final_ROD = tail(d$ROD_FRAC, 1), final_CONE = tail(d$CONE_FRAC, 1),
             final_MLMT = tail(d$MLMT, 1), final_FST = tail(d$FST, 1))
    })
    bind_rows(dfs)
  })

  output$rgc_plot <- renderPlot({
    df <- run_sim(input$bio_scenario)
    ggplot(df, aes(time/(24*365), RGC_FRAC)) + geom_line(color = "darkgreen") +
      labs(x = "연", y = "RGC 생존 분율") + theme_minimal()
  })
  output$cme_plot <- renderPlot({
    df <- run_sim(input$bio_scenario)
    ggplot(df, aes(time/(24*365), CME_CST)) + geom_line(color = "darkorange") +
      labs(x = "연", y = "중심망막두께 (um)") + theme_minimal()
  })

  output$visualcycle_plot <- renderPlot({
    df <- run_sim(input$gt_scenario)
    ggplot(df, aes(time/(24*365), visual_cycle)) + geom_line(color = "steelblue") +
      labs(x = "연", y = "시각회로 플럭스 (상대치)", title = "RPE65 시각회로 회복") + theme_minimal()
  })
  output$bypass_plot <- renderPlot({
    df <- run_sim(input$gt_scenario)
    ggplot(df, aes(time/(24*365), op_bypass)) + geom_line(color = "purple") +
      labs(x = "연", y = "광유전학 우회 신호 (상대치)", title = "MCO-010 우회 경로 (RGC 생존율에 비례)") + theme_minimal()
  })
}

shinyApp(ui, server)
