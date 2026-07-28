## =============================================================================
## sah_shiny_app.R
## Aneurysmal Subarachnoid Haemorrhage -> Delayed Cerebral Ischaemia (aSAH-DCI)
## Interactive QSP dashboard  ·  동맥류성 지주막하출혈 후 지연성 뇌허혈 대시보드
##
## 8 tabs:
##   1  환자 프로파일        Patient profile & the shared reserve
##   2  약물 PK              Drug pharmacokinetics
##   3  혈액-헤모글로빈 시계  The haemoglobin clock (why day 4-10 exists)
##   4  4개 소비자           The four consumers of one reserve
##   5  관류와 산소           Perfusion, oxygen delivery & extraction
##   6  임상 엔드포인트       Clinical endpoints
##   7  시나리오 비교         Scenario comparison (15 arms)
##   8  바이오마커 · 모니터링  Biomarkers and bedside monitoring
##
## The dashboard is deliberately built around the model's central claim: the
## reserve gauge on tab 1 and the stacked consumer plot on tab 4 are the two
## figures that make "the angiogram is not the endpoint" visible.
##
## Run:  shiny::runApp("sah_shiny_app.R")
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("sah_mrgsolve_model.R")   # provides sah_mod, sah_population, scenarios,
                                 # sah_simulate(), sah_endpoints(), helpers

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#eef2f5", colour = NA),
        legend.position = "bottom")

ARM_CHOICES <- setNames(names(sah_scenario_defs),
                        vapply(sah_scenario_defs, function(x) x$label, ""))

## =============================================================================
## UI
## =============================================================================
ui <- fluidPage(
  titlePanel(HTML(paste0(
    "<b>aSAH &#8594; DCI QSP Dashboard</b> &nbsp;",
    "<span style='font-size:14px;color:#555'>동맥류성 지주막하출혈 후 지연성 뇌허혈 ",
    "— 하나의 세동맥 예비능, 네 명의 소비자</span>"))),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      sliderInput("age", "Age (years)", 25, 88, 56, 1),
      selectInput("wfns", "WFNS grade", choices = 1:5, selected = 3),
      selectInput("mfi", "modified Fisher grade", choices = 1:4, selected = 3),
      checkboxInput("hp22", "Hp2-2 haptoglobin genotype", FALSE),
      sliderInput("hgb", "Haemoglobin (g/dL)", 7, 17, 12.6, 0.1),
      sliderInput("map0", "Baseline MAP (mmHg)", 70, 132, 96, 1),
      sliderInput("aregb", "Autoregulatory gain (1 = intact)", 0.10, 0.99, 0.86, 0.01),
      hr(),
      h4("표현형 (Phenotype multipliers)"),
      sliderInput("spsn", "Large-artery vasoreactivity", 0.3, 3.0, 1.0, 0.05),
      sliderInput("rmsn", "Microvascular reactivity", 0.3, 3.0, 1.0, 0.05),
      sliderInput("thrp", "Microthrombosis propensity", 0.2, 4.0, 1.0, 0.05),
      sliderInput("vrisk", "Tissue volume at risk (mL)", 10, 150, 42, 1),
      hr(),
      h4("치료 (Treatment)"),
      checkboxInput("nim", "Oral nimodipine 60 mg q4h", TRUE),
      checkboxInput("nimiv", "IV nimodipine 2 mg/h instead", FALSE),
      sliderInput("clz", "Clazosentan (mg/h)", 0, 15, 0, 1),
      checkboxInput("cil", "Cilostazol 100 mg bid", FALSE),
      checkboxInput("stat", "Simvastatin 40 mg", FALSE),
      checkboxInput("drain", "Early lumbar / cisternal drainage", FALSE),
      checkboxInput("nic", "Intrathecal nicardipine implant", FALSE),
      checkboxInput("ket", "Ketamine (SD suppression)", FALSE),
      sliderInput("mil", "Milrinone (ug/kg/min, d4-14)", 0, 1.5, 0, 0.1),
      sliderInput("maph", "Induced hypertension (+mmHg, d4-14)", 0, 30, 0, 5),
      sliderInput("hgbtgt", "Transfusion threshold (g/dL, 0 = off)", 0, 12, 0, 0.5),
      hr(),
      h4("가상 임상시험 (Virtual trial)"),
      sliderInput("npop", "Population size", 100, 900, 300, 50),
      actionButton("run_pop", "Run population", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## ---------------------------------------------------------- tab 1
        tabPanel(
          "1. 환자 프로파일",
          br(),
          fluidRow(
            column(7, h4("공유 예비능 게이지 (the shared reserve)"),
                   plotOutput("reserve_gauge", height = "300px")),
            column(5, h4("요약"), tableOutput("profile_tbl"))
          ),
          hr(),
          h4("핵심 상태변수 (key states over 21 days)"),
          plotOutput("profile_states", height = "420px"),
          helpText(HTML(paste0(
            "<b>How to read this tab.</b> The gauge shows the arteriolar reserve ",
            "(RA0 &#8594; RAmin) and how much of it each consumer has taken. ",
            "Ischaemia begins where the stacked demand crosses the top of the gauge — ",
            "not where the angiogram looks worst.")))
        ),

        ## ---------------------------------------------------------- tab 2
        tabPanel(
          "2. 약물 PK",
          br(),
          plotOutput("pk_plot", height = "480px"),
          hr(),
          h4("정상상태 노출과 수용체 점유 (steady-state exposure and effect)"),
          tableOutput("pk_tbl"),
          helpText(HTML(paste0(
            "Nimodipine's oral bioavailability is ~13% because of CYP3A4 first pass, ",
            "which is why the IV route reaches the same exposure with a much larger ",
            "MAP penalty — the same drug, a different harm channel.")))
        ),

        ## ---------------------------------------------------------- tab 3
        tabPanel(
          "3. 헤모글로빈 시계",
          br(),
          plotOutput("clock_plot", height = "460px"),
          hr(),
          fluidRow(
            column(6, h4("전이 사슬 (the transit chain)"),
                   plotOutput("clock_chain", height = "300px")),
            column(6, h4("소거 용량 (clearance capacity)"),
                   plotOutput("clock_clear", height = "300px"))
          ),
          helpText(HTML(paste0(
            "<b>There is no day-4 switch in this model.</b> The window is the ",
            "convolution of a transit (clot &#8594; erythrocytes &#8594; free oxyHb) ",
            "with an inducible sink (HO-1, t&#189; 1.9 d). Heavier clot moves the peak ",
            "earlier and higher; Hp2-2 makes the sink weaker; drainage truncates the ",
            "INPUT rather than blocking an effector.")))
        ),

        ## ---------------------------------------------------------- tab 4
        tabPanel(
          "4. 네 개의 소비자",
          br(),
          h4("예비능 소비 분해 (stacked reserve consumption)"),
          plotOutput("consumers_stack", height = "420px"),
          hr(),
          fluidRow(
            column(6, h4("소비자별 지분 (share of total demand)"),
                   plotOutput("consumers_share", height = "320px")),
            column(6, h4("허혈 여유 (ischaemic margin = CBFmax / CBFrequired)"),
                   plotOutput("margin_plot", height = "320px"))
          ),
          helpText(HTML(paste0(
            "<b>The redundancy theorem, visually.</b> Drag the clazosentan slider up: ",
            "the orange (large-artery) band shrinks dramatically and the stack barely ",
            "drops below the reserve line, because the other three consumers simply ",
            "take up the slack. That is the whole CONSCIOUS-2 result in one figure.")))
        ),

        ## ---------------------------------------------------------- tab 5
        tabPanel(
          "5. 관류와 산소",
          br(),
          plotOutput("perf_plot", height = "500px"),
          hr(),
          fluidRow(
            column(6, h4("저항 성분 (resistance components)"),
                   plotOutput("resist_plot", height = "300px")),
            column(6, h4("추출 한계 (CTH degrades OEFmax)"),
                   plotOutput("oef_plot", height = "300px"))
          ),
          helpText(HTML(paste0(
            "CaO2 (anaemia, transfusion) and CMRO2 (fever, spreading depolarisation) ",
            "enter the SAME node as vasospasm. A patient at Hb 8 g/dL is consuming ",
            "the reserve just as surely as one with a 40% calibre loss.")))
        ),

        ## ---------------------------------------------------------- tab 6
        tabPanel(
          "6. 임상 엔드포인트",
          br(),
          fluidRow(
            column(6, h4("이 환자의 엔드포인트"), tableOutput("endpoint_tbl")),
            column(6, h4("결과 확률 분해 (outcome logit contributions)"),
                   plotOutput("outcome_waterfall", height = "320px"))
          ),
          hr(),
          h4("가상 집단 엔드포인트 (virtual population, current treatment)"),
          plotOutput("pop_endpoints", height = "380px"),
          DTOutput("pop_tbl")
        ),

        ## ---------------------------------------------------------- tab 7
        tabPanel(
          "7. 시나리오 비교",
          br(),
          checkboxGroupInput(
            "arms", "Arms to compare",
            choices = ARM_CHOICES,
            selected = c("S1", "S2", "S4", "S6", "S8"), inline = TRUE),
          actionButton("run_arms", "Run selected arms", class = "btn-primary"),
          hr(),
          h4("발생률과 위험비 (incidence and risk ratios vs S2)"),
          DTOutput("arm_tbl"),
          br(),
          plotOutput("arm_plot", height = "420px"),
          helpText(HTML(paste0(
            "<b>The diagnostic figure of this whole library entry.</b> Plot ",
            "RR(angiographic vasospasm) against RR(poor outcome) across arms. If ",
            "vasospasm were the mechanism, the points would lie on the diagonal. ",
            "They do not: clazosentan sits at the far left (angiography won) and at ",
            "1.0 vertically (endpoint unmoved).")))
        ),

        ## ---------------------------------------------------------- tab 8
        tabPanel(
          "8. 바이오마커 · 모니터링",
          br(),
          plotOutput("biomarker_plot", height = "460px"),
          hr(),
          fluidRow(
            column(6, h4("진단 성능 (TCD vs PbtO2 as DCI tests)"),
                   tableOutput("test_tbl")),
            column(6, h4("선행 시간 (lead times in DCI patients)"),
                   tableOutput("lead_tbl"))
          ),
          helpText(HTML(paste0(
            "TCD velocity is a measurement of the ONE consumer angiography can see, ",
            "and it is confounded by flow: when CBF falls, velocity falls too, so a ",
            "severely spastic, severely hypoperfused vessel can read deceptively low. ",
            "PbtO2 measures the shared node itself, which is why it dominates as a test.")))
        )
      )
    )
  )
)

## =============================================================================
## SERVER
## =============================================================================
server <- function(input, output, session) {

  ## ---- one-subject parameter set --------------------------------------------
  subject <- reactive({
    data.frame(
      ID = 1, AGE = input$age, WFNS = as.numeric(input$wfns),
      MFI = as.numeric(input$mfi), HP22 = as.numeric(input$hp22),
      MAP0 = input$map0, HGB0 = input$hgb, THRP = input$thrp,
      SPSN = input$spsn, RMSN = input$rmsn, AREGB = input$aregb,
      EBIB = min(max(0.08 + 0.155 * (as.numeric(input$wfns) - 1), 0), 1),
      RA0i = 0.85, RAMINi = 0.20, CLF = 1, VRISKi = input$vrisk,
      EVDF = as.numeric(as.numeric(input$mfi) %in% c(2, 4))
    )
  })

  treatment <- reactive({
    sc <- list(label = "custom", TSTART = 0.5, TSTOP = 21)
    if (input$nim && !input$nimiv) sc$NIMPO <- 60
    if (input$nimiv) sc$NIMIV <- 2
    if (input$clz > 0) sc$CLZR <- input$clz
    if (input$cil) sc$CILD <- 100
    if (input$stat) sc$STATON <- 1
    if (input$drain) { sc$DRAIN <- 2.20; sc$DRAINHB <- 1.85 }
    if (input$nic) sc$NICD <- 4
    if (input$ket) sc$KETON <- 1
    if (input$mil > 0) sc$MILR <- input$mil
    if (input$maph > 0) sc$MAPH <- input$maph
    if (input$hgbtgt > 0) sc$HGBTGT <- input$hgbtgt
    sc
  })

  sim1 <- reactive(sah_simulate(subject(), treatment(), end = 21, delta = 0.05))

  ## ---- tab 1 ---------------------------------------------------------------
  output$reserve_gauge <- renderPlot({
    s <- sim1()
    d <- s %>%
      select(time, DLARGE, DMICRO, DTHROMB, DCPPd, DSDd, RESERVE) %>%
      pivot_longer(DLARGE:DSDd, names_to = "consumer", values_to = "demand") %>%
      mutate(consumer = recode(consumer,
        DLARGE = "1 large-artery spasm", DMICRO = "2 microvascular tone",
        DTHROMB = "3 microthrombosis", DCPPd = "4a CPP loss",
        DSDd = "4b SD inverse coupling"))
    ggplot(d, aes(time, demand, fill = consumer)) +
      geom_area(colour = "white", linewidth = 0.2) +
      geom_hline(aes(yintercept = RESERVE), linetype = "dashed",
                 linewidth = 1.1, colour = "#b8860b") +
      annotate("text", x = 1, y = max(d$RESERVE) * 1.06, hjust = 0,
               label = "vasodilatory reserve", colour = "#b8860b", size = 4) +
      scale_fill_manual(values = c("#e8871a", "#2a9d8f", "#c1121f",
                                   "#457b9d", "#5f0f40")) +
      labs(x = "day", y = "reserve demand (resistance units)", fill = NULL) +
      THEME
  })

  output$profile_tbl <- renderTable({
    s <- sim1(); ep <- sah_endpoints(s, subject())
    data.frame(
      metric = c("peak calibre loss (%)", "peak TCD (cm/s)", "CBF nadir",
                 "PbtO2 nadir (mmHg)", "min ischaemic margin",
                 "infarct volume (mL)", "P(mRS 4-6)", "Na nadir (mmol/L)"),
      value = c(sprintf("%.0f", 100 * ep$spt_max), sprintf("%.0f", ep$tcd_max),
                sprintf("%.1f", min(s$CBFo)), sprintf("%.1f", ep$pbt_min),
                sprintf("%.2f", min(s$MARGIN)), sprintf("%.1f", ep$infvol),
                sprintf("%.1f%%", 100 * ep$ppoor), sprintf("%.1f", ep$na_min)))
  })

  output$profile_states <- renderPlot({
    s <- sim1()
    d <- s %>% select(time, OXYHB, ET1, SPTOT, RMIC, MTHR, AREG, ISCHo, INFVOL) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.9, colour = "#1d3557") +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "day", y = NULL) + THEME
  })

  ## ---- tab 2 ---------------------------------------------------------------
  output$pk_plot <- renderPlot({
    s <- sim1()
    d <- s %>% select(time, NIMC, CLAZ, CILO, MILRC, NICA, KETA) %>%
      pivot_longer(-time) %>% filter(value > 1e-9)
    if (!nrow(d)) return(ggplot() + annotate("text", 1, 1, label = "no drug selected") +
                         theme_void())
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = "concentration (ng/mL)", colour = NULL) +
      THEME + theme(legend.position = "none")
  })

  output$pk_tbl <- renderTable({
    s <- sim1(); w <- s$time > 5 & s$time < 13
    data.frame(
      drug = c("nimodipine", "clazosentan", "cilostazol", "milrinone",
               "IT nicardipine", "ketamine"),
      Css_ng_mL = c(mean(s$NIMC[w]), mean(s$CLAZ[w]), mean(s$CILO[w]),
                    mean(s$MILRC[w]), mean(s$NICA[w]), mean(s$KETA[w])),
      effect_fraction = c(
        mean(s$NIMC[w]) / (mean(s$NIMC[w]) + 30),
        mean(s$CLAZ[w]) / (mean(s$CLAZ[w]) + 400),
        mean(s$CILO[w]) / (mean(s$CILO[w]) + 600),
        mean(s$MILRC[w]) / (mean(s$MILRC[w]) + 150),
        mean(s$NICA[w]) / (mean(s$NICA[w]) + 2000),
        mean(s$KETA[w]) / (mean(s$KETA[w]) + 1000)))
  }, digits = 3)

  ## ---- tab 3 ---------------------------------------------------------------
  output$clock_plot <- renderPlot({
    s <- sim1()
    d <- s %>% select(time, CLOT, RBCL, OXYHB, HEME, HO1, BOX, ET1, SPTOT) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.9, colour = "#6a1b9a") +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      geom_vline(xintercept = c(4, 10), linetype = "dotted", colour = "grey40") +
      labs(x = "day", y = NULL,
           subtitle = "dotted lines = the conventional 'day 4-10' window, which the model never encodes") +
      THEME
  })

  output$clock_chain <- renderPlot({
    s <- sim1()
    d <- s %>% select(time, CLOT, RBCL, OXYHB) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1.1) +
      labs(x = "day", y = "relative amount", colour = NULL) + THEME
  })

  output$clock_clear <- renderPlot({
    s <- sim1()
    d <- s %>% select(time, HP, HO1) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1.1) +
      labs(x = "day", y = "capacity (relative)", colour = NULL,
           subtitle = "HO-1 is INDUCIBLE — its lag is why free oxyHb persists into week 2") +
      THEME
  })

  ## ---- tab 4 ---------------------------------------------------------------
  output$consumers_stack <- renderPlot({
    s <- sim1()
    d <- s %>% select(time, DLARGE, DMICRO, DTHROMB, DCPPd, DSDd, RESERVE) %>%
      pivot_longer(DLARGE:DSDd, names_to = "consumer", values_to = "demand")
    ggplot(d, aes(time, demand, fill = consumer)) +
      geom_area(colour = "white", linewidth = 0.2) +
      geom_hline(aes(yintercept = RESERVE), linetype = "dashed", linewidth = 1.1,
                 colour = "#b8860b") +
      scale_fill_brewer(palette = "Set2") +
      labs(x = "day", y = "reserve demand", fill = NULL) + THEME
  })

  output$consumers_share <- renderPlot({
    s <- sim1()
    d <- s %>% select(time, DLARGE, DMICRO, DTHROMB, DCPPd, DSDd) %>%
      pivot_longer(-time, names_to = "consumer", values_to = "demand") %>%
      group_by(time) %>% mutate(share = demand / sum(demand)) %>% ungroup()
    ggplot(d, aes(time, share, fill = consumer)) +
      geom_area(position = "fill", colour = "white", linewidth = 0.2) +
      scale_fill_brewer(palette = "Set2") +
      scale_y_continuous(labels = scales::percent) +
      labs(x = "day", y = "share of total demand", fill = NULL) + THEME
  })

  output$margin_plot <- renderPlot({
    s <- sim1()
    ggplot(s, aes(time, MARGIN)) + geom_line(linewidth = 1.1, colour = "#1d3557") +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "#c1121f") +
      annotate("text", x = 1, y = 1.03, hjust = 0, label = "ischaemic threshold",
               colour = "#c1121f", size = 3.6) +
      labs(x = "day", y = "CBFmax / CBFrequired") + THEME
  })

  ## ---- tab 5 ---------------------------------------------------------------
  output$perf_plot <- renderPlot({
    s <- sim1()
    d <- s %>% select(time, CPPo, CBFo, ISCHo, PBTO2, OEFo, TCD, ICP, MAP) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.9, colour = "#01579b") +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "day", y = NULL) + THEME
  })

  output$resist_plot <- renderPlot({
    s <- sim1()
    d <- s %>% select(time, RLo, RAo, RMo) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, fill = name)) +
      geom_area(colour = "white", linewidth = 0.2) +
      scale_fill_manual(values = c(RLo = "#e8871a", RAo = "#ffd400", RMo = "#2a9d8f")) +
      labs(x = "day", y = "resistance (mmHg per mL/100g/min)", fill = NULL) + THEME
  })

  output$oef_plot <- renderPlot({
    s <- sim1()
    d <- s %>% select(time, CTHo, OEFMAXo, OEFo) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1.0) +
      labs(x = "day", y = NULL, colour = NULL,
           subtitle = "when OEF reaches OEFmax there is no extraction reserve left") +
      THEME
  })

  ## ---- tab 6 ---------------------------------------------------------------
  output$endpoint_tbl <- renderTable({
    ep <- sah_endpoints(sim1(), subject())
    data.frame(
      endpoint = c("any angiographic vasospasm (>=25%)",
                   "moderate-severe (>=33%)", "severe (>=50%)",
                   "DCI", "infarct volume (mL)", "P(mRS 4-6)"),
      value = c(ifelse(ep$ang_any > 0.5, "yes", "no"),
                ifelse(ep$ang_ms > 0.5, "yes", "no"),
                ifelse(ep$ang_sev > 0.5, "yes", "no"),
                ifelse(ep$dci > 0.5, "yes", "no"),
                sprintf("%.1f", ep$infvol), sprintf("%.1f%%", 100 * ep$ppoor)))
  })

  output$outcome_waterfall <- renderPlot({
    ep <- sah_endpoints(sim1(), subject())
    b <- c(EBI = 3.15 * ep$EBIB, infarct = 1.08 * log1p(ep$infvol),
           ICP = 1.05 * min(max((ep$icp_max - 15) / 15, 0), 2),
           age = 0.036 * (ep$AGE - 56),
           harm = 1.55 * ep$harm, SD = 0.85 * log1p(ep$sdcum) / 3,
           hypoNa = 0.45 * (ep$na_min < 130))
    d <- data.frame(term = names(b), logit = as.numeric(b))
    ggplot(d, aes(reorder(term, logit), logit)) +
      geom_col(fill = "#457b9d") + coord_flip() +
      labs(x = NULL, y = "contribution to outcome logit") + THEME
  })

  pop_res <- eventReactive(input$run_pop, {
    pop <- sah_population(input$npop)
    sah_endpoints(sah_simulate(pop, treatment()), pop)
  })

  output$pop_endpoints <- renderPlot({
    ep <- pop_res()
    d <- data.frame(
      endpoint = c("any vasospasm", "moderate-severe", "severe", "DCI",
                   "poor outcome"),
      pct = 100 * c(mean(ep$ang_any), mean(ep$ang_ms), mean(ep$ang_sev),
                    mean(ep$dci), mean(ep$ppoor)))
    ggplot(d, aes(reorder(endpoint, pct), pct)) +
      geom_col(fill = "#2a9d8f") + coord_flip() +
      geom_text(aes(label = sprintf("%.1f%%", pct)), hjust = -0.1) +
      expand_limits(y = max(d$pct) * 1.15) +
      labs(x = NULL, y = "% of virtual population") + THEME
  })

  output$pop_tbl <- renderDT({
    ep <- pop_res()
    datatable(ep %>% select(ID, MFI, HP22, AREGB, spt_max, tcd_max, pbt_min,
                            dci, infvol, ppoor) %>%
                mutate(across(where(is.numeric), ~round(.x, 3))),
              options = list(pageLength = 8, scrollX = TRUE), rownames = FALSE)
  })

  ## ---- tab 7 ---------------------------------------------------------------
  arm_res <- eventReactive(input$run_arms, {
    pop <- sah_population(input$npop)
    arms <- input$arms
    if (!"S2" %in% arms) arms <- c("S2", arms)
    eps <- lapply(arms, function(k) sah_endpoints(sah_simulate(pop, sah_scenario_defs[[k]]), pop))
    names(eps) <- arms
    ref <- eps[["S2"]]
    do.call(rbind, lapply(arms, function(k) {
      e <- eps[[k]]
      data.frame(arm = k, label = sah_scenario_defs[[k]]$label,
                 ang_ms = 100 * mean(e$ang_ms), dci = 100 * mean(e$dci),
                 poor = 100 * mean(e$ppoor),
                 rr_ang_ms = mean(e$ang_ms) / mean(ref$ang_ms),
                 rr_dci = mean(e$dci) / mean(ref$dci),
                 rr_poor = mean(e$ppoor) / mean(ref$ppoor))
    }))
  })

  output$arm_tbl <- renderDT({
    datatable(arm_res() %>% mutate(across(where(is.numeric), ~round(.x, 3))),
              options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
  })

  output$arm_plot <- renderPlot({
    d <- arm_res()
    ggplot(d, aes(rr_ang_ms, rr_poor, label = arm)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
      geom_hline(yintercept = 1, linewidth = 0.3) +
      geom_point(aes(size = dci), colour = "#c1121f", alpha = 0.8) +
      geom_text(vjust = -1.1, size = 4) +
      scale_size_continuous(range = c(3, 10)) +
      labs(x = "RR of moderate-severe angiographic vasospasm",
           y = "RR of poor outcome (mRS 4-6)", size = "DCI %",
           subtitle = paste("If vasospasm were the mechanism every arm would sit on the",
                            "dashed diagonal. The distance from it is the redundancy.")) +
      THEME
  })

  ## ---- tab 8 ---------------------------------------------------------------
  output$biomarker_plot <- renderPlot({
    s <- sim1()
    d <- s %>% select(time, TCD, LINDEG, PBTO2, OXYHB, ET1, INFL, SDCUM, NAS) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.9, colour = "#4e342e") +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "day", y = NULL) + THEME
  })

  output$test_tbl <- renderTable({
    ep <- pop_res()
    tests <- function(x, thr, above = TRUE) {
      pos <- if (above) x > thr else x < thr
      c(sensitivity = 100 * sum(pos & ep$dci > 0.5) / max(sum(ep$dci > 0.5), 1),
        specificity = 100 * sum(!pos & ep$dci < 0.5) / max(sum(ep$dci < 0.5), 1))
    }
    rbind(
      data.frame(test = "TCD > 120 cm/s",  t(tests(ep$tcd_max, 120))),
      data.frame(test = "TCD > 160 cm/s",  t(tests(ep$tcd_max, 160))),
      data.frame(test = "TCD > 200 cm/s",  t(tests(ep$tcd_max, 200))),
      data.frame(test = "PbtO2 < 20 mmHg", t(tests(ep$pbt_min, 20, FALSE))),
      data.frame(test = "PbtO2 < 15 mmHg", t(tests(ep$pbt_min, 15, FALSE))))
  }, digits = 1)

  output$lead_tbl <- renderTable({
    pop <- sah_population(min(input$npop, 300))
    sim <- sah_simulate(pop, treatment())
    ep  <- sah_endpoints(sim, pop)
    ids <- ep$ID[ep$dci > 0.5]
    if (!length(ids)) return(data.frame(note = "no DCI in this population"))
    s <- sim[sim$ID %in% ids, ]
    f <- function(col, thr, above, lab) {
      v <- vapply(split(s, s$ID), function(d) {
        x <- if (identical(col, "dINF")) c(0, diff(d$INFVOL)) else d[[col]]
        ok <- (if (above) x > thr else x < thr) & d$time > 1
        if (any(ok)) d$time[which(ok)[1]] else NA_real_
      }, numeric(1))
      v <- v[!is.na(v)]
      data.frame(signal = lab, median_day = median(v), n = length(v))
    }
    rbind(f("TCD", 120, TRUE, "TCD > 120 cm/s"),
          f("PBTO2", 20, FALSE, "PbtO2 < 20 mmHg"),
          f("dINF", 0.05, TRUE, "infarct starts"))
  }, digits = 1)
}

shinyApp(ui, server)
