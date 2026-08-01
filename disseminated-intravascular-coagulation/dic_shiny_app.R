## =====================================================================
##  dic_shiny_app.R
##  Disseminated Intravascular Coagulation (DIC) — QSP explorer
##  파종성 혈관내 응고 · 정량적 시스템 약리학 대시보드
##
##  Ten tabs, organised around the model's structural claim: DIC is one
##  consumption process read through TWO independently-set clocks.
##
##    1  Patient & trigger      — set the aetiology, i.e. set the clocks
##    2  Clock 1 · acute phase  — why septic DIC does not bleed
##    3  Clock 2 · fibrinolysis — why APL DIC does
##    4  Coagulation PD         — thrombin, factors, the consumption clock
##    5  Fibrin: rate vs stock  — the distinction drugs are blind to
##    6  Drug PK                — every agent's exposure on one axis
##    7  The three products     — heparin/AT, APC, haemostatic failure
##    8  Diagnostic scores      — ISTH · SIC · JAAM as OUTPUTS
##    9  Scenario comparison    — the 16 therapeutic arms side by side
##   10  Trial reconciliation   — where the model agrees and disagrees
##
##  Run with:  shiny::runApp("dic_shiny_app.R")
##  Requires:  shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  Sources:   dic_mrgsolve_model.R (must sit in the same directory)
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("dic_mrgsolve_model.R", local = TRUE)

PAL <- c(sepsis = "#C62828", apl = "#1565C0", obstetric = "#6A1B9A",
         sepsis_noDIC = "#2E7D32", healthy = "#546E7A", treated = "#EF6C00")

theme_dic <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}

lineplot <- function(df, ycol, ylab, hlines = NULL, logy = FALSE) {
  p <- ggplot(df, aes(time / 24, .data[[ycol]], colour = arm)) +
    geom_line(linewidth = 0.9) +
    labs(x = "Days from presentation", y = ylab, colour = NULL) +
    theme_dic()
  if (!is.null(hlines))
    p <- p + geom_hline(yintercept = hlines, linetype = "dashed",
                        colour = "grey45", linewidth = 0.3)
  if (logy) p <- p + scale_y_log10()
  p
}

## =====================================================================
##  UI
## =====================================================================
ui <- navbarPage(
  title = "DIC QSP — 파종성 혈관내 응고",
  theme = NULL,
  header = tags$style(HTML(
    ".well{background:#F7F8FA;} .nav-tabs{font-weight:600;}
     .noteblock{background:#FFF9E6;border-left:4px solid #F5A623;
                padding:10px 14px;margin:8px 0;font-size:13px;}")),

  ## ------------------------------------------------------------------
  tabPanel(
    "1 · 환자 & 유발인자",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h4("Aetiology — this is where the two clocks are set"),
        radioButtons("aet", "기저 질환 (aetiology)",
          choices = c("Septic DIC"                = "sepsis",
                      "Same sepsis, DIC-resistant"= "sepsis_noDIC",
                      "APL (acute promyelocytic leukaemia)" = "apl",
                      "Obstetric / major trauma"  = "obstetric",
                      "Healthy control"           = "healthy"),
          selected = "sepsis"),
        hr(),
        h5("CLOCK 1 — acute-phase sign"),
        sliderInput("SIL6", "IL-6 gain (SIL6): sepsis 1.0, APL 0.10",
                    0, 1.5, 1.0, 0.05),
        sliderInput("APOSF", "Max IL-6-driven rise in fibrinogen synthesis (APOSF)",
                    0, 6, 3, 0.25),
        h5("CLOCK 2 — PAI-1 : t-PA"),
        sliderInput("SPAI", "PAI-1 induction gain (SPAI): sepsis 1.0, APL 0.15",
                    0, 1.5, 1.0, 0.05),
        sliderInput("KANX", "Annexin-A2 amplification of t-PA on blasts (KANX)",
                    0, 6, 3, 0.25),
        hr(),
        h5("Host & severity"),
        sliderInput("PATH0", "Initial pathogen / DAMP load", 0, 2, 1, 0.05),
        sliderInput("SRCCTL", "Source control (multiplier on clearance)",
                    0.5, 4, 1, 0.1),
        sliderInput("BLAST0", "APL blast burden at diagnosis (x10^9/L)",
                    0, 100, 30, 5),
        sliderInput("MSUP", "Marrow output multiplier (0.20 in APL)",
                    0.05, 1.2, 1.0, 0.05),
        numericInput("WTkg", "Body weight (kg)", 70, 35, 150, 1),
        numericInput("endh", "Simulation horizon (h)", 672, 48, 1344, 24)
      ),
      mainPanel(
        width = 8,
        div(class = "noteblock",
            strong("Read this tab as the model's central claim."),
            " Thrombotic and haemorrhagic DIC are not two diseases with two",
            " models. They are the same 49 equations at two settings of",
            " SIL6/APOSF (clock 1) and SPAI/KANX (clock 2). Everything on the",
            " other nine tabs follows arithmetically from what you set here."),
        plotOutput("p_trigger", height = "320px"),
        plotOutput("p_endothel", height = "300px"),
        tableOutput("t_profile")
      )
    )
  ),

  ## ------------------------------------------------------------------
  tabPanel(
    "2 · Clock 1 · 급성기 반응",
    fluidRow(
      column(12, div(class = "noteblock",
        strong("Why septic DIC does not bleed."),
        " Fibrinogen synthesis is multiplied by (1 + APOSF x APF); antithrombin,",
        " protein C and free protein S are DIVIDED by (1 + ANEG x APF). One",
        " IL-6 signal, two opposite signs. At the SAME consumption rate the",
        " simulated fibrinogen steady state is 281 mg/dL in sepsis and 87 mg/dL",
        " in APL — so the septic patient has unopposed thrombin WITH substrate",
        " (microthrombosis) and the APL patient has neither (haemorrhage)."))),
    fluidRow(
      column(6, plotOutput("p_clock1_fib", height = "300px")),
      column(6, plotOutput("p_clock1_anticoag", height = "300px"))),
    fluidRow(
      column(6, plotOutput("p_halflife", height = "320px")),
      column(6, tableOutput("t_clock1")))
  ),

  ## ------------------------------------------------------------------
  tabPanel(
    "3 · Clock 2 · 섬유소용해",
    fluidRow(
      column(12, div(class = "noteblock",
        strong("Why APL DIC does bleed."),
        " free t-PA = tPA_nM x (FLOC + (1-FLOC)/(1 + PAI1_nM/KPAI)) x annexin.",
        " In simulated sepsis total t-PA rises 4-fold and free t-PA still FALLS",
        " to 0.66x of healthy, because PAI-1 goes to a 42:1 molar excess.",
        " In APL the ratio is 3.9:1 and free t-PA is 13x the septic value."))),
    fluidRow(
      column(6, plotOutput("p_clock2_pai", height = "300px")),
      column(6, plotOutput("p_clock2_free", height = "300px"))),
    fluidRow(
      column(6, plotOutput("p_plasmin", height = "300px")),
      column(6, plotOutput("p_a2ap", height = "300px")))
  ),

  ## ------------------------------------------------------------------
  tabPanel(
    "4 · 응고 PD",
    fluidRow(
      column(6, plotOutput("p_thrombin", height = "300px")),
      column(6, plotOutput("p_kcons", height = "300px"))),
    fluidRow(
      column(6, plotOutput("p_factors", height = "320px")),
      column(6, plotOutput("p_platelet", height = "320px")))
  ),

  ## ------------------------------------------------------------------
  tabPanel(
    "5 · 피브린: 속도 vs 재고",
    fluidRow(
      column(12, div(class = "noteblock",
        strong("The distinction most DIC drugs are blind to."),
        " Every anticoagulant acts on the deposition RATE. Organ failure is",
        " driven by the deposited STOCK. The two are different state",
        " variables with different time constants, so a drug's effect on the",
        " first does not translate one-for-one into an effect on the second."))),
    fluidRow(
      column(6, plotOutput("p_rate_stock", height = "320px")),
      column(6, plotOutput("p_organ", height = "320px"))),
    fluidRow(column(12, plotOutput("p_ddimer", height = "280px")))
  ),

  ## ------------------------------------------------------------------
  tabPanel(
    "6 · 약물 PK",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        h4("Dosing"),
        checkboxInput("useUFH", "UFH infusion", FALSE),
        sliderInput("doseUFH", "U/kg/h", 0, 30, 18, 1),
        checkboxInput("useENOX", "Enoxaparin SC q12h", FALSE),
        sliderInput("doseENOX", "mg/kg", 0, 1.5, 1, 0.1),
        checkboxInput("useATC", "Antithrombin concentrate", FALSE),
        sliderInput("doseATC", "loading IU/kg", 0, 150, 86, 2),
        checkboxInput("useRTM", "Thrombomodulin alfa (ART-123)", FALSE),
        sliderInput("doseRTM", "mg/kg/day", 0, 0.12, 0.06, 0.01),
        checkboxInput("useAPC", "Drotrecogin alfa", FALSE),
        sliderInput("doseAPC", "ug/kg/h", 0, 48, 24, 2),
        sliderInput("startAPC", "start time (h)", 0, 120, 6, 3),
        checkboxInput("useTXA", "Tranexamic acid", FALSE),
        checkboxInput("useARG", "Argatroban", FALSE),
        checkboxInput("useATRA", "ATRA (APL only)", FALSE),
        checkboxInput("useFIBC", "Fibrinogen concentrate 4 g q12h", FALSE),
        checkboxInput("usePLTX", "Platelets 1 unit/day", FALSE)
      ),
      mainPanel(
        width = 9,
        plotOutput("p_pk", height = "380px"),
        plotOutput("p_pkpd", height = "320px"))
    )
  ),

  ## ------------------------------------------------------------------
  tabPanel(
    "7 · 세 개의 곱",
    fluidRow(
      column(12, div(class = "noteblock",
        strong("Three places where the pharmacology is a PRODUCT, and the"),
        strong(" clinical consequence follows from the arithmetic.")))),
    fluidRow(
      column(6, h5("(1) Heparin is a catalyst, and AT saturates it"),
                plotOutput("p_hepceiling", height = "300px")),
      column(6, h5("(2) APC generation = kPCA x THR x TM x EPCR x PC"),
                plotOutput("p_apcproduct", height = "300px"))),
    fluidRow(
      column(6, h5("(3) Haemostatic failure: the worst term dominates"),
                plotOutput("p_bleedprod", height = "300px")),
      column(6, tableOutput("t_products")))
  ),

  ## ------------------------------------------------------------------
  tabPanel(
    "8 · 진단 점수",
    fluidRow(
      column(12, div(class = "noteblock",
        strong("The scores are OUTPUTS of the state vector, not inputs."),
        " The ISTH score awards fibrinogen a point only below 100 mg/dL — and",
        " in sepsis fibrinogen never gets there, for the reason on tab 2.",
        " Watch the 24-hour mark in a septic run: ISTH reads 2 (not overt)",
        " while SIC already reads 4. That gap is why SIC was proposed, and the",
        " model reproduces it from the equations rather than being told."))),
    fluidRow(
      column(7, plotOutput("p_scores", height = "340px")),
      column(5, plotOutput("p_score_parts", height = "340px"))),
    fluidRow(column(12, DTOutput("t_labs")))
  ),

  ## ------------------------------------------------------------------
  tabPanel(
    "9 · 시나리오 비교",
    sidebarLayout(
      sidebarPanel(width = 3,
        checkboxGroupInput("arms", "Arms to compare", choices = names(SCEN),
          selected = c("C_septic_DIC", "D_UFH", "F_AT_alone",
                       "H_thrombomodulin", "I_rhAPC_early", "L_TXA_sepsis")),
        selectInput("cmpvar", "Variable",
          choices = c("Platelets" = "PLTc", "Fibrinogen" = "FIBc", "INR" = "INR",
                      "D-dimer" = "DDc", "Antithrombin %" = "ATpc",
                      "Protein C %" = "PCpc", "Thrombin (nM)" = "THR",
                      "Plasmin (nM)" = "PLN", "Fibrin stock" = "FDEP",
                      "SOFA" = "SOFA", "ISTH" = "ISTH", "SIC" = "SIC",
                      "Bleeding index" = "BIDX", "28-d mortality %" = "MORT"),
          selected = "FDEP"),
        actionButton("runall", "Run all arms", class = "btn-primary")),
      mainPanel(width = 9,
        plotOutput("p_compare", height = "380px"),
        DTOutput("t_compare"))
    )
  ),

  ## ------------------------------------------------------------------
  tabPanel(
    "10 · 임상시험 대조",
    fluidRow(column(12,
      h4("What the model gets right, and what it does not"),
      div(class = "noteblock",
        "A model that only reproduced the trials it was built from would be",
        " worthless. The table below states the disagreements as plainly as",
        " the agreements. The heparin rows are the model's most exposed",
        " prediction and should be read as a falsifiable hypothesis."),
      DTOutput("t_trials"),
      br(),
      h4("Did 'given too late' explain PROWESS-SHOCK?"),
      p("A start-time sweep of drotrecogin alfa. If the 'the fibrin stock",
        " already exists' story were sufficient, benefit would fall",
        " monotonically with start time. In this model it does not: benefit",
        " PEAKS at 24-48 h, tracking peak fibrin deposition, and only",
        " collapses beyond ~72 h. The model refutes that prior hypothesis."),
      plotOutput("p_timing", height = "320px")))
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  ## --- build the run from the sidebar ---------------------------------
  sim <- reactive({
    WTl <- input$WTkg
    par <- list(SIL6 = input$SIL6, SPAI = input$SPAI, APOSF = input$APOSF,
                KANX = input$KANX, SRCCTL = input$SRCCTL, MSUP = input$MSUP)
    A   <- AETIOLOGY[[input$aet]]
    ini <- A$init
    if (!is.null(ini$PATH))  ini$PATH  <- input$PATH0
    if (!is.null(ini$BLAST)) ini$BLAST <- input$BLAST0

    evs <- NULL; infs <- NULL
    if (input$useUFH)
      infs <- inf_c(infs, infusion("RATE_UFH", input$doseUFH * WTl, 0, 336))
    if (input$useENOX)
      evs <- ev_c(evs, bolus("ENXD", input$doseENOX * WTl * 100, seq(0, 335, 12)))
    if (input$useATC)
      evs <- ev_c(evs, bolus("AT", 0.014 * input$doseATC, 0),
                       bolus("AT", 0.014 * input$doseATC / 2, seq(12, 84, 12)))
    if (input$useRTM)
      evs <- ev_c(evs, bolus("RTM", input$doseRTM * WTl, seq(0, 120, 24)))
    if (input$useAPC)
      infs <- inf_c(infs, infusion("RATE_APC", input$doseAPC * WTl,
                                   input$startAPC, input$startAPC + 96))
    if (input$useTXA) {
      evs  <- ev_c(evs, bolus("TXAC", 1000, 0))
      infs <- inf_c(infs, infusion("RATE_TXA", 125, 0, 96))
    }
    if (input$useARG)  infs <- inf_c(infs, infusion("RATE_ARG", 8.4, 0, 336))
    if (input$useATRA) evs <- ev_c(evs, bolus("ATRAD", 40, seq(0, input$endh, 12)))
    if (input$useFIBC) evs <- ev_c(evs, bolus("FIB", 100, seq(0, 72, 12)))
    if (input$usePLTX) evs <- ev_c(evs, bolus("PLT", 40, seq(0, 144, 24)))

    AET_LOCAL <- AETIOLOGY
    AET_LOCAL[[input$aet]]$init <- ini
    old <- AETIOLOGY; AETIOLOGY <<- AET_LOCAL
    on.exit(AETIOLOGY <<- old, add = TRUE)

    out <- run_dic(aet = input$aet, events = evs, infusions = infs,
                   par = par, end = input$endh)
    out$arm <- "current run"
    out
  })

  ## --- reference untreated run of the same aetiology -------------------
  ref <- reactive({
    r <- run_dic(aet = input$aet, par = list(SIL6 = input$SIL6,
                 SPAI = input$SPAI, APOSF = input$APOSF, KANX = input$KANX),
                 end = input$endh)
    r$arm <- "untreated reference"; r
  })
  both <- reactive(bind_rows(ref(), sim()))

  ## ---- tab 1 ----------------------------------------------------------
  output$p_trigger <- renderPlot({
    d <- sim() %>% select(time, PATH, BLAST, IL6, TNF, HIST, TF) %>%
      pivot_longer(-time)
    ggplot(d, aes(time / 24, value)) + geom_line(colour = PAL[["sepsis"]], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", nrow = 2) +
      labs(x = "Days", y = NULL, title = "Trigger and innate amplification") +
      theme_dic()
  })
  output$p_endothel <- renderPlot({
    d <- sim() %>% select(time, TM, EPCR, GLX) %>% pivot_longer(-time)
    ggplot(d, aes(time / 24, 100 * value, colour = name)) + geom_line(linewidth = 1) +
      labs(x = "Days", y = "% of normal", colour = NULL,
           title = "Endothelial anticoagulant surface — the brakes are on the target") +
      theme_dic()
  })
  output$t_profile <- renderTable({
    s <- sim()
    idx <- sapply(c(24, 48, 96, 168), function(h) which.min(abs(s$time - h)))
    data.frame(Hour = s$time[idx],
               Platelets = round(s$PLTc[idx]), Fibrinogen = round(s$FIBc[idx]),
               INR = round(s$INR[idx], 2), `D-dimer` = round(s$DDc[idx], 1),
               `AT %` = round(s$ATpc[idx]), `PC %` = round(s$PCpc[idx]),
               ISTH = s$ISTH[idx], SIC = s$SIC[idx], SOFA = round(s$SOFA[idx], 1),
               check.names = FALSE)
  })

  ## ---- tab 2 : clock 1 -------------------------------------------------
  output$p_clock1_fib <- renderPlot(
    lineplot(both(), "FIBc", "Fibrinogen (mg/dL)", hlines = c(100, 150)) +
      ggtitle("Fibrinogen — a POSITIVE acute-phase reactant"))
  output$p_clock1_anticoag <- renderPlot({
    d <- both() %>% select(time, arm, AT = ATpc, `Protein C` = PCpc) %>%
      pivot_longer(c(AT, `Protein C`))
    ggplot(d, aes(time / 24, value, colour = arm, linetype = name)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Days", y = "% of normal", colour = NULL, linetype = NULL,
           title = "AT and protein C — NEGATIVE acute-phase reactants") + theme_dic()
  })
  output$p_halflife <- renderPlot({
    d <- data.frame(
      protein = c("Protein C", "Factor VII", "Factor VIII", "Factor V",
                  "Factor X", "Antithrombin", "Prothrombin", "Fibrinogen", "Platelets"),
      thalf   = c(6, 5, 10, 15, 40, 65, 65, 100, 200))
    d$protein <- factor(d$protein, levels = d$protein[order(d$thalf)])
    ggplot(d, aes(protein, thalf)) +
      geom_col(fill = "#4DB6AC") + coord_flip() +
      labs(x = NULL, y = "Plasma half-life (h)",
           title = "Turnover ordering: protein C reaches its floor in ~22 h,\nfibrinogen takes ~5 days") +
      theme_dic()
  })
  output$t_clock1 <- renderTable({
    s <- sim(); i <- which.min(abs(s$time - 48))
    data.frame(
      Quantity = c("Acute-phase factor APF", "Fibrinogen synthesis multiplier",
                   "Fibrinogen (mg/dL)", "Antithrombin (%)", "Protein C (%)"),
      `At 48 h` = c(round(s$APFo[i], 3), round(s$FIBSYN[i], 2),
                    round(s$FIBc[i]), round(s$ATpc[i]), round(s$PCpc[i])),
      check.names = FALSE)
  })

  ## ---- tab 3 : clock 2 -------------------------------------------------
  output$p_clock2_pai <- renderPlot({
    d <- both() %>% select(time, arm, `t-PA (ng/mL)` = TPA, `PAI-1 (ng/mL)` = PAI1) %>%
      pivot_longer(c(-time, -arm))
    ggplot(d, aes(time / 24, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Days", y = NULL, colour = NULL,
           title = "t-PA rises — and PAI-1 rises far more") + theme_dic()
  })
  output$p_clock2_free <- renderPlot(
    lineplot(both(), "FREETPA", "Free t-PA (nM)") +
      ggtitle("Free t-PA: what actually reaches plasminogen"))
  output$p_plasmin <- renderPlot(
    lineplot(both(), "PLN", "Plasmin (nM)") + ggtitle("Plasmin"))
  output$p_a2ap <- renderPlot(
    lineplot(both(), "A2pc", "alpha-2-antiplasmin (%)", hlines = 50) +
      ggtitle("alpha-2-antiplasmin — systemic fibrinogenolysis begins only\nonce this is exhausted"))

  ## ---- tab 4 -----------------------------------------------------------
  output$p_thrombin <- renderPlot(
    lineplot(both(), "THR", "Thrombin (nM)") + ggtitle("Thrombin — the hub variable"))
  output$p_kcons <- renderPlot(
    lineplot(both(), "SFM", "Soluble fibrin monomer (nM)") +
      ggtitle("Soluble fibrin monomer — the earliest marker"))
  output$p_factors <- renderPlot({
    d <- sim() %>% select(time, FII, FV, FVII, FX, TFPI, PS) %>% pivot_longer(-time)
    ggplot(d, aes(time / 24, 100 * value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "Days", y = "% of normal", colour = NULL,
           title = "One consumption clock, different half-lives") + theme_dic()
  })
  output$p_platelet <- renderPlot(
    lineplot(both(), "PLTc", "Platelets (x10^9/L)", hlines = c(50, 100, 150)) +
      ggtitle("Platelets"))

  ## ---- tab 5 -----------------------------------------------------------
  output$p_rate_stock <- renderPlot(
    lineplot(both(), "FDEP", "Deposited microvascular fibrin (AU)") +
      ggtitle("The STOCK — what organ failure is actually driven by"))
  output$p_organ <- renderPlot({
    d <- sim() %>% select(time, Kidney = OKID, Liver = OLIV,
                          Lung = OLUN, Brain = OCNS) %>% pivot_longer(-time)
    ggplot(d, aes(time / 24, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "Days", y = "Injury index (0-1)", colour = NULL,
           title = "Organ injury — liver injury feeds back on synthesis") + theme_dic()
  })
  output$p_ddimer <- renderPlot(
    lineplot(both(), "DDc", "D-dimer (mg/L FEU)", hlines = c(0.5, 5)) +
      ggtitle("D-dimer requires BOTH cross-linked fibrin AND plasmin — a product"))

  ## ---- tab 6 -----------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim() %>%
      transmute(time,
        `UFH (anti-Xa IU/mL)`   = HEPC / (3.5 * 1000),
        `Enoxaparin (IU/mL)`    = ENXC / (5 * 1000),
        `rTM (ug/mL)`           = RTM / 3.5,
        `Drotrecogin (nM)`      = APCX / 10 * 0.017857,
        `Tranexamic acid (ug/mL)` = TXAC / 12,
        `Argatroban (ug/mL)`    = ARGA / 12,
        `ATRA (mg/L)`           = ATRAC / 100) %>%
      pivot_longer(-time) %>% filter(value > 1e-9)
    if (!nrow(d)) return(ggplot() + labs(title = "No drug selected") + theme_dic())
    ggplot(d, aes(time / 24, value)) + geom_line(colour = PAL[["treated"]], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Days", y = NULL, title = "Drug exposure") + theme_dic()
  })
  output$p_pkpd <- renderPlot({
    d <- both() %>% select(time, arm, Thrombin = THR, `Fibrin stock` = FDEP,
                           `Bleeding index` = BIDX) %>% pivot_longer(c(-time, -arm))
    ggplot(d, aes(time / 24, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "Days", y = NULL, colour = NULL, title = "PD response") + theme_dic()
  })

  ## ---- tab 7 : the three products --------------------------------------
  output$p_hepceiling <- renderPlot({
    P <- as.list(param(mod))
    at <- seq(0.05, 2.4, 0.01)
    cap <- function(chep) {
      amp <- 1 + P$EMAXHI * chep / (P$EC50H + chep)
      atc <- at / (P$KMAT + at) * (P$KMAT + 1)
      P$KATII * (at + (amp - 1) * atc)
    }
    d <- rbind(data.frame(AT = 100 * at, rate = cap(0),   hep = "no heparin"),
               data.frame(AT = 100 * at, rate = cap(0.3), hep = "UFH 0.3 IU/mL"),
               data.frame(AT = 100 * at, rate = cap(0.7), hep = "UFH 0.7 IU/mL"))
    ggplot(d, aes(AT, rate, colour = hep)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = c(40, 100), linetype = "dashed", colour = "grey40") +
      labs(x = "Antithrombin activity (%)", y = "Thrombin inactivation rate (1/h)",
           colour = NULL,
           title = "Heparin is a CATALYST: the AT term saturates.\nAt AT 40% the ceiling is 72% of normal; AT 220% buys only +26%.") +
      theme_dic()
  })
  output$p_apcproduct <- renderPlot({
    d <- sim() %>% select(time, TM, EPCR, PC) %>%
      mutate(Product = TM * EPCR * PC) %>% pivot_longer(-time)
    ggplot(d, aes(time / 24, value, colour = name,
                  linewidth = name == "Product")) +
      geom_line() + scale_linewidth_manual(values = c(0.7, 1.6), guide = "none") +
      labs(x = "Days", y = "Fraction of normal", colour = NULL,
           title = "Three fractions multiply. At 48 h in septic DIC:\n0.318 x 0.559 x 0.333 = 0.059 of normal.") +
      theme_dic()
  })
  output$p_bleedprod <- renderPlot({
    s <- sim(); P <- as.list(param(mod))
    d <- data.frame(time = s$time,
      Platelets  = 1 / (1 + (s$PLTc / P$KBPLT)^P$GBP),
      Fibrinogen = 1 / (1 + (s$FIBc / P$KBFIB)^P$GBF),
      `PT/factors` = { x <- pmax(s$PT - 12, 0) / P$KBPT; x^P$GBT / (1 + x^P$GBT) },
      Lysis      = s$PLN / (P$KBLYS + s$PLN),
      check.names = FALSE) %>% pivot_longer(-time)
    ggplot(d, aes(time / 24, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "Days", y = "Component hazard", colour = NULL,
           title = "The worst component dominates the product —\nso correct the worst variable first") +
      theme_dic()
  })
  output$t_products <- renderTable({
    s <- sim(); i <- which.min(abs(s$time - 48))
    data.frame(
      Product = c("TM x EPCR x PC", "APC per unit thrombin",
                  "PAI-1 : t-PA (molar)", "Free t-PA (nM)"),
      `At 48 h` = c(round(s$APCPROD[i], 3), signif(s$APCTHR[i], 3),
                    round(s$PAITPA[i], 1), signif(s$FREETPA[i], 3)),
      check.names = FALSE)
  })

  ## ---- tab 8 -----------------------------------------------------------
  output$p_scores <- renderPlot({
    d <- sim() %>% select(time, ISTH, SIC, SOFA) %>% pivot_longer(-time)
    ggplot(d, aes(time / 24, value, colour = name)) + geom_step(linewidth = 0.9) +
      geom_hline(yintercept = 5, linetype = "dashed", colour = "#C62828") +
      geom_hline(yintercept = 4, linetype = "dotted", colour = "#1565C0") +
      labs(x = "Days", y = "Score", colour = NULL,
           title = "ISTH >= 5 = overt DIC (red);  SIC >= 4 (blue)") + theme_dic()
  })
  output$p_score_parts <- renderPlot({
    s <- sim()
    d <- data.frame(time = s$time,
      Platelets  = ifelse(s$PLTc < 50, 2, ifelse(s$PLTc < 100, 1, 0)),
      `D-dimer`  = ifelse(s$DDc > 5, 3, ifelse(s$DDc > 1, 2, 0)),
      PT         = ifelse(s$PT - 12 > 6, 2, ifelse(s$PT - 12 > 3, 1, 0)),
      Fibrinogen = ifelse(s$FIBc < 100, 1, 0),
      check.names = FALSE) %>% pivot_longer(-time)
    ggplot(d, aes(time / 24, value, fill = name)) +
      geom_area(position = "stack") +
      labs(x = "Days", y = "ISTH points", fill = NULL,
           title = "Where the ISTH points come from —\nnote the fibrinogen band is empty in sepsis") +
      theme_dic()
  })
  output$t_labs <- renderDT({
    s <- sim()
    idx <- seq(1, nrow(s), by = max(1, round(12 / (s$time[2] - s$time[1]))))
    datatable(s[idx, c("time", "PLTc", "FIBc", "PT", "INR", "DDc", "ATpc",
                       "PCpc", "A2pc", "THR", "PLN", "FDEP", "SOFA",
                       "ISTH", "SIC", "BIDX", "MORT")] %>%
                mutate(across(where(is.numeric), ~ round(.x, 2))),
              options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)
  })

  ## ---- tab 9 -----------------------------------------------------------
  allarms <- eventReactive(input$runall, {
    withProgress(message = "Running scenarios", value = 0, {
      n <- length(input$arms)
      bind_rows(lapply(seq_along(input$arms), function(k) {
        incProgress(1 / n, detail = input$arms[k])
        run_scenario(input$arms[k])
      }))
    })
  })
  output$p_compare <- renderPlot({
    d <- allarms(); req(d)
    ggplot(d, aes(time / 24, .data[[input$cmpvar]], colour = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(x = "Days", y = input$cmpvar, colour = NULL,
           title = "Therapeutic scenario comparison") + theme_dic()
  })
  output$t_compare <- renderDT({
    d <- allarms(); req(d)
    datatable(summarise_dic(d, hours = c(48, 672)),
              options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
  })

  ## ---- tab 10 ----------------------------------------------------------
  output$t_trials <- renderDT({
    datatable(data.frame(
      Intervention = c("Drotrecogin alfa (rhAPC)", "Antithrombin, no heparin",
                       "Antithrombin + heparin", "Thrombomodulin alfa",
                       "ATRA in APL", "Antifibrinolytic + ATRA in APL",
                       "UFH in septic DIC"),
      Trial = c("PROWESS 2001", "KyberSept 2001 (subgroup)", "KyberSept 2001",
                "SCARLET 2019", "APL trials since 1988", "Avvisati 1989 / ELN guidance",
                "meta-analyses"),
      `Trial result` = c("24.7% vs 30.8% (-6.1 pts)", "37.8% vs 43.6% (-5.8 pts)",
                         "38.9% vs 38.7% (0), bleeding 22.0% vs 12.8%",
                         "26.8% vs 29.4% (-2.6, p=0.32)",
                         "early haemorrhagic death cut by ~2/3 to 3/4",
                         "no benefit; thrombosis reported with ATRA",
                         "at best a few points in DIC subgroups"),
      `Model result` = c("-5.9 pts", "-7.3 pts",
                         "-5.3 pts increment, higher bleeding burden",
                         "+0.5 pts (essentially neutral)",
                         "32.8% -> 7.5%", "7.5% -> 21.8% (harm)", "-15.0 pts"),
      Verdict = c("agrees", "agrees", "partly agrees (direction, not magnitude)",
                  "agrees", "agrees", "agrees",
                  "DISAGREES — the model's most exposed prediction"),
      check.names = FALSE),
      options = list(dom = "t", scrollX = TRUE), rownames = FALSE) %>%
      formatStyle("Verdict", target = "row",
                  backgroundColor = styleEqual(
                    "DISAGREES — the model's most exposed prediction", "#FFE0E0"))
  })
  output$p_timing <- renderPlot({
    d <- data.frame(start_h = c(0, 6, 12, 24, 30, 48, 72, 96),
                    MORT = c(37.8, 37.4, 36.9, 36.1, 35.8, 35.9, 37.3, 39.1))
    ggplot(d, aes(start_h, MORT)) +
      geom_line(linewidth = 1, colour = PAL[["treated"]]) + geom_point(size = 2.5) +
      geom_hline(yintercept = 43.3, linetype = "dashed", colour = "#C62828") +
      annotate("text", x = 70, y = 43.8, label = "no drug: 43.3%",
               colour = "#C62828", size = 3.5) +
      labs(x = "Drotrecogin alfa start time (h)", y = "28-day mortality (%)",
           title = "Benefit peaks at 24-48 h, tracking peak fibrin deposition") +
      theme_dic()
  })
}

shinyApp(ui, server)
