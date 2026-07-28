## ============================================================================
## Anorexia Nervosa (AN) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 9 tabs: Patient profile · Energy balance & body composition · Neuroendocrine
##         axes · Drug PK/PD · Refeeding safety · Bone health · Clinical
##         endpoints · Scenario comparison · References
## Dependencies: shiny, shinydashboard, mrgsolve, dplyr, tidyr, ggplot2, DT
##
## Run with:  shiny::runApp("an_shiny_app.R")
##            (an_mrgsolve_model.R must sit in the same directory)
## ----------------------------------------------------------------------------
## EDUCATIONAL / RESEARCH USE ONLY. Not a refeeding protocol and not a
## substitute for specialist eating-disorder medical care.
## ============================================================================
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DT)
})

## ---------- Lazy model load -------------------------------------------------
get_model <- function() {
  if (!exists(".AN_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.AN_MOD)) {
    assign(".AN_MOD",
           mread_cache(model = "an", file = "an_mrgsolve_model.R", project = "."),
           envir = .GlobalEnv)
  }
  .GlobalEnv$.AN_MOD
}

## ---------- Scenario definitions -------------------------------------------
SCENARIOS <- list(
  "1. Untreated natural history" = list(
    par = list(EI_START = 0, EI_ESCAL = 0, INPATIENT = 0),
    tx = 0, olz = FALSE, flx = FALSE, e2 = FALSE),
  "2. Inpatient higher-calorie refeeding (2000 kcal start)" = list(
    par = list(EI_START = 2000, EI_ESCAL = 200, EI_PLAN_MAX = 3400, INPATIENT = 1),
    tx = 1.4, olz = FALSE, flx = FALSE, e2 = FALSE),
  "3. Conservative slow refeeding (1000 kcal start)" = list(
    par = list(EI_START = 1000, EI_ESCAL = 60, EI_PLAN_MAX = 3000, INPATIENT = 1),
    tx = 1.4, olz = FALSE, flx = FALSE, e2 = FALSE),
  "4. Steep escalation, no phosphate prophylaxis" = list(
    par = list(EI_START = 2400, EI_ESCAL = 400, EI_PLAN_MAX = 3600, INPATIENT = 1,
               ILLNESS_DUR = 4, BASE_BMI = 12.5),
    tx = 1.4, olz = FALSE, flx = FALSE, e2 = FALSE),
  "5. Steep escalation + PO4/Mg/thiamine prophylaxis" = list(
    par = list(EI_START = 2400, EI_ESCAL = 400, EI_PLAN_MAX = 3600, INPATIENT = 1,
               ILLNESS_DUR = 4, BASE_BMI = 12.5,
               PHOS_SUP = 0.9, MG_SUP = 0.25, THIA_SUP = 0.10),
    tx = 1.4, olz = FALSE, flx = FALSE, e2 = FALSE),
  "6. Outpatient FBT (adolescent)" = list(
    par = list(EI_START = 2200, EI_ESCAL = 40, ADOLESCENT = 1, FBT = 1),
    tx = 1.6, olz = FALSE, flx = FALSE, e2 = FALSE),
  "7. Outpatient CBT-E (adult)" = list(
    par = list(EI_START = 2200, EI_ESCAL = 40, ADOLESCENT = 0, FBT = 0),
    tx = 1.6, olz = FALSE, flx = FALSE, e2 = FALSE),
  "8. Olanzapine 10 mg + outpatient therapy" = list(
    par = list(EI_START = 2200, EI_ESCAL = 40),
    tx = 1.6, olz = TRUE, flx = FALSE, e2 = FALSE),
  "9. Fluoxetine 60 mg (nutritional gate demo)" = list(
    par = list(EI_START = 2200, EI_ESCAL = 40),
    tx = 1.6, olz = FALSE, flx = TRUE, e2 = FALSE),
  "10. Transdermal E2 + discharge at day 90" = list(
    par = list(EI_START = 2000, EI_ESCAL = 200, INPATIENT = 1,
               DISCHARGE_DAY = 90, POST_DC_FAC = 0.72),
    tx = 1.4, olz = FALSE, flx = FALSE, e2 = TRUE)
)

build_events <- function(spec, horizon) {
  e <- NULL
  add <- function(a, b) if (is.null(a)) b else a + b
  n_wk <- max(floor(horizon / 7), 1)
  if (spec$tx > 0) {
    n_tx <- if (!is.null(spec$par$DISCHARGE_DAY) && spec$par$DISCHARGE_DAY > 0)
      max(floor(min(spec$par$DISCHARGE_DAY, horizon) / 7), 1) else n_wk
    e <- add(e, ev(amt = spec$tx, cmt = "THERAPY", ii = 7, addl = n_tx - 1))
  }
  if (isTRUE(spec$olz)) e <- add(e, ev(amt = 1.0, cmt = "OLZ_GUT", ii = 1, addl = horizon - 1))
  if (isTRUE(spec$flx)) e <- add(e, ev(amt = 1.0, cmt = "FLX_CENT", ii = 1, addl = horizon - 1))
  if (isTRUE(spec$e2))  e <- add(e, ev(amt = 1.0, cmt = "E2_PATCH", ii = 3.5, addl = floor(horizon / 3.5) - 1))
  if (is.null(e)) e <- ev(amt = 0, cmt = "THERAPY")
  e
}

run_sim <- function(scen_name, horizon, user_par, extra = list()) {
  spec <- SCENARIOS[[scen_name]]
  mod  <- get_model()
  par  <- modifyList(modifyList(spec$par, user_par), extra)
  mod %>% param(par) %>%
    mrgsim(events = build_events(spec, horizon), end = horizon, delta = 0.25) %>%
    as_tibble() %>% mutate(scenario = scen_name)
}

theme_an <- function() {
  theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          panel.grid.minor = element_blank())
}

lineplot <- function(df, y, ylab, hline = NULL, hlab = NULL) {
  p <- ggplot(df, aes(x = time, y = .data[[y]])) +
    geom_line(linewidth = 0.9, colour = "#2c6fbb") +
    labs(x = "Time (days)", y = ylab) + theme_an()
  if (!is.null(hline))
    p <- p + geom_hline(yintercept = hline, linetype = 2, colour = "#c0392b")
  p
}

multiplot <- function(df, ys, labels, ylab) {
  d <- df %>% select(time, all_of(ys)) %>%
    pivot_longer(-time, names_to = "var", values_to = "val") %>%
    mutate(var = factor(var, levels = ys, labels = labels))
  ggplot(d, aes(time, val, colour = var)) + geom_line(linewidth = 0.9) +
    labs(x = "Time (days)", y = ylab) + theme_an()
}

## ---------- UI --------------------------------------------------------------
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "Anorexia Nervosa QSP", titleWidth = 300),
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("1. Patient profile",        tabName = "patient", icon = icon("user")),
      menuItem("2. Energy & body comp",     tabName = "energy",  icon = icon("weight-scale")),
      menuItem("3. Neuroendocrine axes",    tabName = "endo",    icon = icon("dna")),
      menuItem("4. Drug PK/PD",             tabName = "pk",      icon = icon("pills")),
      menuItem("5. Refeeding safety",       tabName = "refeed",  icon = icon("triangle-exclamation")),
      menuItem("6. Bone health",            tabName = "bone",    icon = icon("bone")),
      menuItem("7. Clinical endpoints",     tabName = "clin",    icon = icon("notes-medical")),
      menuItem("8. Scenario comparison",    tabName = "compare", icon = icon("chart-line")),
      menuItem("9. References",             tabName = "refs",    icon = icon("book"))
    ),
    hr(),
    selectInput("scenario", "Scenario:", names(SCENARIOS),
                selected = "2. Inpatient higher-calorie refeeding (2000 kcal start)"),
    sliderInput("horizon", "Horizon (days):", 30, 730, 365, step = 5),
    sliderInput("base_bmi", "Presenting BMI (kg/m2):", 10, 18, 14, step = 0.25),
    sliderInput("median_bmi", "Median BMI for age/height:", 18, 23, 20.5, step = 0.1),
    sliderInput("ht", "Height (m):", 1.40, 1.90, 1.65, step = 0.01),
    sliderInput("sev", "ED cognition severity (0-1):", 0, 1, 0.80, step = 0.05),
    sliderInput("dur", "Illness duration (years):", 0, 15, 1.5, step = 0.5),
    checkboxInput("purge", "Binge/purge subtype", FALSE),
    checkboxInput("adolescent", "Adolescent", TRUE),
    checkboxInput("smoker", "Smoker (CYP1A2 induction)", FALSE),
    h5("Repletion / add-ons", style = "padding-left:15px;color:#eee"),
    sliderInput("phos_sup", "Phosphate repletion (mg/dL/day):", 0, 2, 0, step = 0.1),
    sliderInput("k_sup", "Potassium repletion (mEq/L/day):", 0, 1.5, 0, step = 0.1),
    sliderInput("thia_sup", "Thiamine repletion (1/day):", 0, 0.3, 0, step = 0.02),
    actionButton("run", "Run simulation", icon = icon("play"),
                 style = "color:#fff;background:#6f42c1;margin:10px 15px;")
  ),
  dashboardBody(
    tabItems(
      ## ---- 1. Patient profile ----
      tabItem("patient",
        fluidRow(
          valueBoxOutput("vb_bmi", width = 3),
          valueBoxOutput("vb_pmbmi", width = 3),
          valueBoxOutput("vb_hr", width = 3),
          valueBoxOutput("vb_phos", width = 3)
        ),
        fluidRow(
          box(width = 7, title = "Baseline patient state", status = "primary",
              solidHeader = TRUE, DTOutput("patient_table")),
          box(width = 5, title = "How this model is wired", status = "info",
              solidHeader = TRUE,
              p(strong("Anorexia nervosa"), "is modeled as an eating-disorder",
                "cognition drive that suppresses intake below requirement. The",
                "resulting energy deficit and depleted body stores are combined",
                "into a single starvation index that drives every downstream",
                "axis: leptin, T3, cortisol, GnRH/estradiol, IGF-1, bone",
                "remodeling and heart rate."),
              p(strong("Two clinically important feedbacks are explicit:"),
                tags$ul(
                  tags$li("Starvation perpetuates the cognition (Minnesota effect) — ",
                          "so weight restoration is a prerequisite for, not a ",
                          "consequence of, psychological recovery."),
                  tags$li("SSRI efficacy is gated by nutritional state, which is ",
                          "why fluoxetine does nothing in the underweight patient."))),
              p(strong("Safety:"), "the refeeding insulin surge is driven by intake",
                "ABOVE the tissue adaptation state, so escalation rate — not",
                "absolute calories — sets the hypophosphatemia risk."))
        )
      ),

      ## ---- 2. Energy balance & body composition ----
      tabItem("energy",
        fluidRow(
          box(width = 6, title = "Prescribed plan vs actual intake vs expenditure",
              status = "primary", solidHeader = TRUE, plotOutput("p_intake", height = 300)),
          box(width = 6, title = "Net energy balance (kcal/day)",
              status = "primary", solidHeader = TRUE, plotOutput("p_eb", height = 300))
        ),
        fluidRow(
          box(width = 6, title = "BMI and % median BMI", status = "success",
              solidHeader = TRUE, plotOutput("p_bmi", height = 300)),
          box(width = 6, title = "Fat mass / fat-free mass (kg)", status = "success",
              solidHeader = TRUE, plotOutput("p_comp", height = 300))
        ),
        fluidRow(
          box(width = 6, title = "Adaptive thermogenesis & compulsive exercise",
              status = "warning", solidHeader = TRUE, plotOutput("p_adapt", height = 300)),
          box(width = 6, title = "ED cognition drive (EDE-Q global proxy)",
              status = "warning", solidHeader = TRUE, plotOutput("p_drive", height = 300))
        )
      ),

      ## ---- 3. Neuroendocrine axes ----
      tabItem("endo",
        fluidRow(
          box(width = 6, title = "Leptin & ghrelin", status = "primary",
              solidHeader = TRUE, plotOutput("p_lepghr", height = 300)),
          box(width = 6, title = "Free T3 (low-T3 syndrome) & cortisol",
              status = "primary", solidHeader = TRUE, plotOutput("p_t3cort", height = 300))
        ),
        fluidRow(
          box(width = 6, title = "GnRH pulse index & estradiol", status = "success",
              solidHeader = TRUE, plotOutput("p_hpg", height = 300)),
          box(width = 6, title = "IGF-1 (GH resistance)", status = "success",
              solidHeader = TRUE, plotOutput("p_igf", height = 300))
        ),
        fluidRow(
          box(width = 12, title = "Menses resumption probability", status = "info",
              solidHeader = TRUE, plotOutput("p_menses", height = 260),
              p(em("Leptin gate at 1.85 ng/mL with a slow pulse-generator recovery ",
                   "reproduces the 3-6 month lag between weight restoration and ",
                   "return of menses.")))
        )
      ),

      ## ---- 4. Drug PK/PD ----
      tabItem("pk",
        fluidRow(
          box(width = 6, title = "Olanzapine exposure & appetite effect",
              status = "primary", solidHeader = TRUE, plotOutput("p_olz", height = 300)),
          box(width = 6, title = "Fluoxetine + norfluoxetine exposure",
              status = "primary", solidHeader = TRUE, plotOutput("p_flx", height = 300))
        ),
        fluidRow(
          box(width = 12, title = "SSRI nutritional gate (why fluoxetine fails at low weight)",
              status = "danger", solidHeader = TRUE, plotOutput("p_gate", height = 300),
              p(em("The gate is a steep Hill function of % median BMI (50% at 90%mBMI). ",
                   "Below weight restoration the same plasma exposure produces ",
                   "essentially no effect on eating-disorder cognition — the model's ",
                   "encoding of Attia 1998 and Walsh 2006.")))
        ),
        fluidRow(
          box(width = 12, title = "Therapy engagement effect (FBT / CBT-E)",
              status = "success", solidHeader = TRUE, plotOutput("p_tx", height = 260))
        )
      ),

      ## ---- 5. Refeeding safety ----
      tabItem("refeed",
        fluidRow(
          valueBoxOutput("vb_nadir", width = 4),
          valueBoxOutput("vb_rsrisk", width = 4),
          valueBoxOutput("vb_qtc", width = 4)
        ),
        fluidRow(
          box(width = 6, title = "Serum phosphate (dashed = 2.5 mg/dL threshold)",
              status = "danger", solidHeader = TRUE, plotOutput("p_phos", height = 300)),
          box(width = 6, title = "Potassium & magnesium", status = "danger",
              solidHeader = TRUE, plotOutput("p_kmg", height = 300))
        ),
        fluidRow(
          box(width = 6, title = "Refeeding insulin-surge index & occult depletion",
              status = "warning", solidHeader = TRUE, plotOutput("p_surge", height = 300)),
          box(width = 6, title = "Thiamine stores & composite refeeding risk",
              status = "warning", solidHeader = TRUE, plotOutput("p_thia", height = 300))
        ),
        fluidRow(
          box(width = 12, title = "Heart rate and QTc", status = "info",
              solidHeader = TRUE, plotOutput("p_cardiac", height = 280))
        )
      ),

      ## ---- 6. Bone health ----
      tabItem("bone",
        fluidRow(
          box(width = 6, title = "Bone turnover markers (P1NP formation / CTX resorption)",
              status = "primary", solidHeader = TRUE, plotOutput("p_markers", height = 320)),
          box(width = 6, title = "BMD Z-score trajectory", status = "primary",
              solidHeader = TRUE, plotOutput("p_bmd", height = 320))
        ),
        fluidRow(
          box(width = 12, title = "Bone-directed therapy comparison", status = "success",
              solidHeader = TRUE, plotOutput("p_bonedrug", height = 320),
              p(em("Compares no bone therapy, transdermal estradiol, oral contraceptive ",
                   "(hepatic first pass suppresses IGF-1 and therefore formation), ",
                   "teriparatide and risedronate on the currently selected patient.")))
        )
      ),

      ## ---- 7. Clinical endpoints ----
      tabItem("clin",
        fluidRow(
          box(width = 12, title = "Endpoint summary (day 0 / 30 / 90 / 180 / end)",
              status = "primary", solidHeader = TRUE, DTOutput("endpoint_table"))
        ),
        fluidRow(
          box(width = 6, title = "Weight restoration & remission flags",
              status = "success", solidHeader = TRUE, plotOutput("p_flags", height = 300)),
          box(width = 6, title = "Medical instability indicator", status = "danger",
              solidHeader = TRUE, plotOutput("p_unstable", height = 300))
        )
      ),

      ## ---- 8. Scenario comparison ----
      tabItem("compare",
        fluidRow(
          box(width = 12, title = "All 10 scenarios on the current patient profile",
              status = "primary", solidHeader = TRUE,
              actionButton("runall", "Run all scenarios", icon = icon("layer-group"),
                           style = "color:#fff;background:#6f42c1;"),
              br(), br(),
              plotOutput("p_cmp_bmi", height = 330))
        ),
        fluidRow(
          box(width = 6, title = "Phosphate", status = "danger", solidHeader = TRUE,
              plotOutput("p_cmp_phos", height = 300)),
          box(width = 6, title = "BMD Z-score", status = "info", solidHeader = TRUE,
              plotOutput("p_cmp_bmd", height = 300))
        ),
        fluidRow(
          box(width = 6, title = "ED cognition drive", status = "warning",
              solidHeader = TRUE, plotOutput("p_cmp_drive", height = 300)),
          box(width = 6, title = "Menses probability", status = "success",
              solidHeader = TRUE, plotOutput("p_cmp_menses", height = 300))
        ),
        fluidRow(
          box(width = 12, title = "Scenario endpoint table", status = "primary",
              solidHeader = TRUE, DTOutput("cmp_table"))
        )
      ),

      ## ---- 9. References ----
      tabItem("refs",
        box(width = 12, title = "References", status = "primary", solidHeader = TRUE,
            p("The full annotated bibliography (35+ PubMed-linked references) lives in ",
              code("an_references.md"), " next to this app."),
            tags$ul(
              tags$li("Keys A et al. The Biology of Human Starvation. 1950 — semi-starvation reproduces AN cognition."),
              tags$li("Hall KD. Lancet 2011;378:826 — energy-partition modeling of body-weight change."),
              tags$li("Misra M, Klibanski A. Lancet Diabetes Endocrinol 2014;2:581 — skeletal effects of AN."),
              tags$li("Misra M et al. J Bone Miner Res 2011;26:2430 — transdermal estradiol improves BMD."),
              tags$li("Fazeli PK et al. J Clin Endocrinol Metab 2014;99:1322 — teriparatide in adult AN."),
              tags$li("Garber AK et al. J Adolesc Health 2021;68:1 — higher-calorie refeeding (StRONG trial)."),
              tags$li("Attia E et al. Am J Psychiatry 2019;176:449 — olanzapine RCT in outpatient AN."),
              tags$li("Walsh BT et al. JAMA 2006;295:2605 — fluoxetine after weight restoration."),
              tags$li("Lock J et al. Arch Gen Psychiatry 2010;67:1025 — FBT vs individual therapy."),
              tags$li("Arcelus J et al. Arch Gen Psychiatry 2011;68:724 — mortality meta-analysis (SMR 5.86).")),
            hr(),
            p(strong("Disclaimer:"), "Educational and research use only. This model is",
              "not validated for clinical decision-making, and nothing in it is a",
              "refeeding protocol. Refeeding a severely malnourished patient requires",
              "specialist inpatient care with electrolyte monitoring."))
      )
    )
  )
)

## ---------- Server ----------------------------------------------------------
server <- function(input, output, session) {

  user_par <- reactive({
    list(HT = input$ht, BASE_BMI = input$base_bmi, MEDIAN_BMI = input$median_bmi,
         SEV = input$sev, ILLNESS_DUR = input$dur,
         PURGE = as.numeric(input$purge), ADOLESCENT = as.numeric(input$adolescent),
         SMOKER = as.numeric(input$smoker),
         PHOS_SUP = input$phos_sup, K_SUP = input$k_sup, THIA_SUP = input$thia_sup)
  })

  sim <- eventReactive(input$run, {
    run_sim(input$scenario, input$horizon, user_par())
  }, ignoreNULL = FALSE)

  all_sims <- eventReactive(input$runall, {
    bind_rows(lapply(names(SCENARIOS), function(s)
      run_sim(s, input$horizon, user_par())))
  })

  ## ---- value boxes ----
  output$vb_bmi <- renderValueBox({
    d <- sim(); v <- tail(d$BMI, 1)
    valueBox(sprintf("%.1f", v), "Final BMI (kg/m2)", icon = icon("weight-scale"),
             color = if (v >= 19) "green" else if (v >= 16) "yellow" else "red")
  })
  output$vb_pmbmi <- renderValueBox({
    d <- sim(); v <- tail(d$Percent_median_BMI, 1)
    valueBox(sprintf("%.0f%%", v), "Final % median BMI", icon = icon("percent"),
             color = if (v >= 95) "green" else if (v >= 85) "yellow" else "red")
  })
  output$vb_hr <- renderValueBox({
    d <- sim(); v <- min(d$Heart_rate)
    valueBox(sprintf("%.0f", v), "Lowest heart rate (bpm)", icon = icon("heart-pulse"),
             color = if (v >= 55) "green" else if (v >= 45) "yellow" else "red")
  })
  output$vb_phos <- renderValueBox({
    d <- sim(); v <- min(d$Phosphate)
    valueBox(sprintf("%.2f", v), "Phosphate nadir (mg/dL)", icon = icon("flask"),
             color = if (v >= 3.0) "green" else if (v >= 2.5) "yellow" else "red")
  })
  output$vb_nadir <- renderValueBox({
    d <- sim(); v <- min(d$Phosphate); t <- d$time[which.min(d$Phosphate)]
    valueBox(sprintf("%.2f @ d%.0f", v, t), "Phosphate nadir (mg/dL)",
             icon = icon("arrow-down"), color = if (v >= 2.5) "green" else "red")
  })
  output$vb_rsrisk <- renderValueBox({
    d <- sim(); v <- max(d$Refeeding_risk)
    valueBox(sprintf("%.2f", v), "Peak refeeding-risk index", icon = icon("triangle-exclamation"),
             color = if (v < 0.3) "green" else if (v < 0.6) "yellow" else "red")
  })
  output$vb_qtc <- renderValueBox({
    d <- sim(); v <- max(d$QTc_ms)
    valueBox(sprintf("%.0f", v), "Peak QTc (ms)", icon = icon("wave-square"),
             color = if (v < 450) "green" else if (v < 480) "yellow" else "red")
  })

  ## ---- 1. patient table ----
  output$patient_table <- renderDT({
    d <- sim(); d0 <- d[1, ]
    datatable(data.frame(
      Variable = c("BMI (kg/m2)", "% median BMI", "Fat mass (kg)", "Fat-free mass (kg)",
                   "Energy intake (kcal/d)", "Energy expenditure (kcal/d)",
                   "Leptin (ng/mL)", "Free T3 (pg/mL)", "Cortisol (ug/dL)",
                   "Estradiol (pg/mL)", "IGF-1 (ng/mL)", "Heart rate (bpm)",
                   "QTc (ms)", "Phosphate (mg/dL)", "Occult PO4 depletion index",
                   "BMD Z-score", "ED drive (0-1)"),
      Baseline = round(c(d0$BMI, d0$Percent_median_BMI, d0$Fat_mass_kg, d0$Fat_free_mass_kg,
                         d0$Energy_intake, d0$Energy_expenditure, d0$Leptin, d0$Free_T3,
                         d0$Cortisol, d0$Estradiol, d0$IGF1_level, d0$Heart_rate,
                         d0$QTc_ms, d0$Phosphate, d0$Phos_depletion_idx,
                         d0$BMD_Zscore, d0$ED_drive), 2)),
      options = list(dom = "t", pageLength = 20), rownames = FALSE)
  })

  ## ---- 2. energy ----
  output$p_intake <- renderPlot(
    multiplot(sim(), c("Prescribed_plan", "Energy_intake", "Energy_expenditure"),
              c("Prescribed plan", "Actual intake", "Total expenditure"), "kcal/day"))
  output$p_eb <- renderPlot(
    lineplot(sim(), "Energy_balance", "Net energy balance (kcal/day)", hline = 0))
  output$p_bmi <- renderPlot(
    multiplot(sim(), c("BMI", "Percent_median_BMI"),
              c("BMI (kg/m2)", "% median BMI"), "value"))
  output$p_comp <- renderPlot(
    multiplot(sim(), c("Fat_mass_kg", "Fat_free_mass_kg", "Body_weight_kg"),
              c("Fat mass", "Fat-free mass", "Body weight"), "kg"))
  output$p_adapt <- renderPlot(
    multiplot(sim(), c("Adaptive_thermogen", "Exercise_kcal"),
              c("Adaptive thermogenesis (fraction)", "Exercise (kcal/day)"), "value"))
  output$p_drive <- renderPlot(
    multiplot(sim(), c("ED_drive", "EDEQ_global"),
              c("ED drive (0-1)", "EDE-Q global (0-6)"), "score"))

  ## ---- 3. endocrine ----
  output$p_lepghr <- renderPlot(
    multiplot(sim(), c("Leptin", "Ghrelin"),
              c("Leptin (ng/mL)", "Ghrelin (pg/mL)"), "concentration"))
  output$p_t3cort <- renderPlot(
    multiplot(sim(), c("Free_T3", "Cortisol"),
              c("Free T3 (pg/mL)", "Cortisol (ug/dL)"), "concentration"))
  output$p_hpg <- renderPlot(
    multiplot(sim(), c("GnRH_index", "Estradiol"),
              c("GnRH pulse index (0-1)", "Estradiol (pg/mL)"), "value"))
  output$p_igf <- renderPlot(lineplot(sim(), "IGF1_level", "IGF-1 (ng/mL)"))
  output$p_menses <- renderPlot(
    lineplot(sim(), "Menses_probability", "P(menses resumed)", hline = 0.5))

  ## ---- 4. PK/PD ----
  output$p_olz <- renderPlot(lineplot(sim(), "Olanzapine_exposure", "Olanzapine exposure (au)"))
  output$p_flx <- renderPlot(lineplot(sim(), "Fluoxetine_total", "Fluoxetine + norfluoxetine (au)"))
  output$p_gate <- renderPlot(
    multiplot(sim(), c("SSRI_gate", "ED_drive"),
              c("Nutritional gate (0-1)", "ED drive (0-1)"), "value"))
  output$p_tx <- renderPlot(lineplot(sim(), "Therapy_effect", "Therapy effect (fractional)"))

  ## ---- 5. refeeding ----
  output$p_phos <- renderPlot(lineplot(sim(), "Phosphate", "Phosphate (mg/dL)", hline = 2.5))
  output$p_kmg <- renderPlot(
    multiplot(sim(), c("Potassium", "Magnesium"),
              c("Potassium (mEq/L)", "Magnesium (mg/dL)"), "concentration"))
  output$p_surge <- renderPlot(
    multiplot(sim(), c("Refeed_surge", "Phos_depletion_idx"),
              c("Insulin-surge index", "Occult PO4 depletion"), "index"))
  output$p_thia <- renderPlot(
    multiplot(sim(), c("Thiamine", "Refeeding_risk"),
              c("Thiamine stores (0-1)", "Composite refeeding risk"), "index"))
  output$p_cardiac <- renderPlot(
    multiplot(sim(), c("Heart_rate", "QTc_ms"),
              c("Heart rate (bpm)", "QTc (ms)"), "value"))

  ## ---- 6. bone ----
  output$p_markers <- renderPlot(
    multiplot(sim(), c("P1NP_marker", "CTX_marker"),
              c("P1NP (ug/L)", "CTX (ng/mL) x1"), "marker"))
  output$p_bmd <- renderPlot(lineplot(sim(), "BMD_Zscore", "BMD Z-score", hline = -2))
  output$p_bonedrug <- renderPlot({
    h <- input$horizon
    base_spec <- SCENARIOS[[input$scenario]]
    mod <- get_model()
    arms <- list(
      "No bone therapy"      = list(ev = NULL, par = list()),
      "Transdermal E2"       = list(ev = ev(amt = 1, cmt = "E2_PATCH", ii = 3.5,
                                            addl = max(floor(h / 3.5) - 1, 0)), par = list()),
      "Oral contraceptive"   = list(ev = NULL, par = list(OCP_ORAL = 1)),
      "Teriparatide 20 ug qd"= list(ev = ev(amt = 1, cmt = "TPTD", ii = 1,
                                            addl = h - 1), par = list()),
      "Risedronate weekly"   = list(ev = ev(amt = 1, cmt = "BIS", ii = 7,
                                            addl = max(floor(h / 7) - 1, 0)), par = list())
    )
    base_ev <- build_events(base_spec, h)
    out <- bind_rows(lapply(names(arms), function(a) {
      e <- if (is.null(arms[[a]]$ev)) base_ev else base_ev + arms[[a]]$ev
      p <- modifyList(modifyList(base_spec$par, user_par()), arms[[a]]$par)
      mod %>% param(p) %>% mrgsim(events = e, end = h, delta = 1) %>%
        as_tibble() %>% mutate(arm = a)
    }))
    ggplot(out, aes(time, BMD_Zscore, colour = arm)) + geom_line(linewidth = 0.9) +
      labs(x = "Time (days)", y = "BMD Z-score") + theme_an()
  })

  ## ---- 7. endpoints ----
  output$endpoint_table <- renderDT({
    d <- sim()
    pick <- function(t) d[which.min(abs(d$time - t)), ]
    ts <- c(0, 30, 90, 180, max(d$time))
    rows <- lapply(ts, pick)
    tab <- data.frame(
      Day = round(ts),
      BMI = sapply(rows, function(r) round(r$BMI, 1)),
      `Pct median BMI` = sapply(rows, function(r) round(r$Percent_median_BMI, 0)),
      `EDE-Q` = sapply(rows, function(r) round(r$EDEQ_global, 2)),
      `Leptin` = sapply(rows, function(r) round(r$Leptin, 2)),
      `Free T3` = sapply(rows, function(r) round(r$Free_T3, 2)),
      `Estradiol` = sapply(rows, function(r) round(r$Estradiol, 1)),
      `IGF-1` = sapply(rows, function(r) round(r$IGF1_level, 0)),
      `BMD Z` = sapply(rows, function(r) round(r$BMD_Zscore, 2)),
      `HR` = sapply(rows, function(r) round(r$Heart_rate, 0)),
      `QTc` = sapply(rows, function(r) round(r$QTc_ms, 0)),
      `PO4` = sapply(rows, function(r) round(r$Phosphate, 2)),
      `P(menses)` = sapply(rows, function(r) round(r$Menses_probability, 2)),
      check.names = FALSE)
    datatable(tab, options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })
  output$p_flags <- renderPlot(
    multiplot(sim(), c("Weight_restored", "Remission"),
              c("Weight restored (>=95% mBMI)", "Remission (+ drive < 0.30)"), "flag"))
  output$p_unstable <- renderPlot(
    lineplot(sim(), "Medically_unstable", "Medically unstable (0/1)"))

  ## ---- 8. comparison ----
  cmp_plot <- function(y, ylab, hline = NULL) {
    d <- all_sims()
    p <- ggplot(d, aes(time, .data[[y]], colour = scenario)) +
      geom_line(linewidth = 0.8) + labs(x = "Time (days)", y = ylab) +
      theme_an() + theme(legend.text = element_text(size = 8))
    if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = 2)
    p
  }
  output$p_cmp_bmi    <- renderPlot(cmp_plot("BMI", "BMI (kg/m2)"))
  output$p_cmp_phos   <- renderPlot(cmp_plot("Phosphate", "Phosphate (mg/dL)", 2.5))
  output$p_cmp_bmd    <- renderPlot(cmp_plot("BMD_Zscore", "BMD Z-score"))
  output$p_cmp_drive  <- renderPlot(cmp_plot("ED_drive", "ED drive (0-1)"))
  output$p_cmp_menses <- renderPlot(cmp_plot("Menses_probability", "P(menses)"))
  output$cmp_table <- renderDT({
    d <- all_sims()
    tab <- d %>% group_by(scenario) %>%
      summarise(`Final BMI` = round(last(BMI), 1),
                `Final %mBMI` = round(last(Percent_median_BMI), 0),
                `PO4 nadir` = round(min(Phosphate), 2),
                `Peak RS risk` = round(max(Refeeding_risk), 2),
                `Lowest HR` = round(min(Heart_rate), 0),
                `Peak QTc` = round(max(QTc_ms), 0),
                `Final BMD Z` = round(last(BMD_Zscore), 2),
                `Final EDE-Q` = round(last(EDEQ_global), 2),
                `P(menses) end` = round(last(Menses_probability), 2),
                .groups = "drop")
    datatable(tab, options = list(dom = "t", pageLength = 12), rownames = FALSE)
  })
}

shinyApp(ui, server)
