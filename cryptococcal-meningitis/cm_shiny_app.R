## =====================================================================
##  HIV-associated Cryptococcal Meningitis -- QSP dashboard
##  Companion Shiny application for cm_mrgsolve_model.R
## ---------------------------------------------------------------------
##  The application is organised around the model's central claim: the
##  disease runs on TWO CLOCKS, and every trial endpoint measures the
##  fast one.  Tab 3 puts them on the same axes so the gap is visible.
##
##  Tabs
##    1  Patient & regimen     -- build a virtual patient and a regimen
##    2  Drug exposure (PK)    -- plasma, CSF, CNS and renal compartments
##    3  TWO CLOCKS            -- viable burden vs capsular antigen
##    4  Pressure & perfusion  -- the Davson volume budget, ICP, LPs
##    5  Immunity & IRIS       -- CD4, CSF WBC, cytokines, IRIS drive
##    6  Clinical endpoints    -- EFA, sterility, mortality, disability
##    7  Scenario comparison   -- all 16 regimens side by side
##    8  Safety               -- GFR, K, Hb, ANC, ALT and the 5FC trap
##    9  Calibration          -- model vs published trials, with residuals
##
##  Run with:  shiny::runApp("cm_shiny_app.R")
##  Requires:  shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("cm_mrgsolve_model.R", local = TRUE)

THEME <- theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey93", colour = NA),
        legend.position = "bottom", panel.grid.minor = element_blank())

CLOCK_COLS <- c("viable burden (clock 1)" = "#2e7d32",
                "capsular antigen (clock 2)" = "#e65100")

## ---------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("HIV-associated Cryptococcal Meningitis — QSP dashboard"),
  tags$p(style = "color:#555;margin-top:-8px",
         tags$em(paste("Two clocks: viable yeast clear in days, capsular",
                       "antigen in weeks. Early fungicidal activity sees",
                       "only the first."))),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Virtual patient"),
      selectInput("phenotype", "Preset phenotype",
                  choices = c("median trial participant" = "median",
                              "upper-quartile burden (high pressure)" = "high_burden",
                              "paucicellular, CD4 12 (COAT high-risk)" = "pauci",
                              "partially preserved immunity" = "immune"),
                  selected = "median"),
      sliderInput("logCFU", "Presenting CSF burden (log10 CFU/mL)",
                  3.0, 7.0, 5.0, step = 0.1),
      sliderInput("cd4", "CD4 count (cells/uL)", 5, 350, 25, step = 5),
      sliderInput("wt", "Body weight (kg)", 35, 100, 60, step = 1),

      hr(), h4("Induction regimen"),
      selectInput("scenario", "Preset regimen",
                  choices = setNames(names(cm_scenarios),
                                     vapply(cm_scenarios,
                                            function(x) x$label, "")),
                  selected = "ambition"),
      checkboxInput("custom", "Override the preset below", FALSE),
      conditionalPanel(
        "input.custom == true",
        radioButtons("ambform", "Amphotericin",
                     c("none" = "none",
                       "deoxycholate, daily" = "d",
                       "liposomal, single 10 mg/kg" = "l"), "l"),
        conditionalPanel("input.ambform == 'd'",
                         sliderInput("ambdose", "AmB-d dose (mg/kg/day)",
                                     0.5, 1.2, 1.0, step = 0.05),
                         sliderInput("ambdays", "days of AmB-d", 1, 28, 7)),
        checkboxInput("use5fc", "flucytosine 100 mg/kg/day", TRUE),
        conditionalPanel("input.use5fc == true",
                         sliderInput("fcdays", "days of flucytosine", 1, 28, 14)),
        sliderInput("flu1", "fluconazole induction (mg/day)", 0, 1200, 1200,
                    step = 100),
        sliderInput("flu2", "fluconazole consolidation (mg/day)", 0, 800, 800,
                    step = 100)
      ),

      hr(), h4("Adjuncts"),
      sliderInput("lpn", "Therapeutic lumbar punctures (n)", 0, 14, 0),
      conditionalPanel("input.lpn > 0",
                       sliderInput("lpevery", "interval between LPs (days)",
                                   1, 4, 2, step = 1)),
      sliderInput("artday", "ART start (day; 200 = never)", 5, 200, 35, step = 1),
      checkboxInput("dex", "adjunctive dexamethasone (CryptoDex taper)", FALSE),
      checkboxInput("sert", "sertraline 400 mg/day (ASTRO-CM)", FALSE),
      checkboxInput("ifng", "interferon gamma 100 ug SC x6", FALSE),

      hr(),
      sliderInput("days", "Simulation horizon (days)", 28, 182, 70, step = 7),
      actionButton("go", "Simulate", class = "btn-primary"),
      tags$hr(),
      tags$small(style = "color:#777",
                 "Educational research model. Not for clinical use.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## ---- 1 ------------------------------------------------------
        tabPanel(
          "1 · Patient & regimen",
          h4("Presenting state (the model's own quasi-equilibrium)"),
          tableOutput("baseTab"),
          h4("Regimen as the model sees it"),
          verbatimTextOutput("regTxt"),
          h4("What the presenting burden sets"),
          plotOutput("basePlot", height = 320),
          tags$p(tags$small(paste(
            "Presenting fungal burden is not just a prognostic marker in this",
            "model: it sets the antigen pool, which sets outflow resistance,",
            "which sets opening pressure. Burden and pressure are therefore",
            "correlated at baseline but decouple during therapy, which is the",
            "point of tab 3.")))
        ),

        ## ---- 2 ------------------------------------------------------
        tabPanel(
          "2 · Drug exposure (PK)",
          fluidRow(
            column(6, plotOutput("pkAmb", height = 300)),
            column(6, plotOutput("pkOral", height = 300))),
          fluidRow(
            column(6, plotOutput("pkBrain", height = 300)),
            column(6, plotOutput("pkKid", height = 300))),
          h4("Exposure summary"),
          tableOutput("pkTab"),
          tags$p(tags$small(paste(
            "The comparison that matters here is CNS AUC against renal",
            "cortical AUC. A single 10 mg/kg liposomal dose matches seven",
            "daily 1 mg/kg deoxycholate doses on the first and roughly",
            "halves the second.")))
        ),

        ## ---- 3 ------------------------------------------------------
        tabPanel(
          "3 · TWO CLOCKS",
          h4("Viable burden and capsular antigen, normalised to presentation"),
          plotOutput("clocks", height = 380),
          fluidRow(
            column(6,
                   h4("Clock 1 — what a CSF culture sees"),
                   plotOutput("clock1", height = 280)),
            column(6,
                   h4("Clock 2 — what no antifungal touches"),
                   plotOutput("clock2", height = 280))),
          h4("The gap, quantified"),
          tableOutput("clockTab"),
          tags$p(tags$small(paste(
            "Time to a negative culture is a drug property. Time for the",
            "antigen pool to halve is a property of its own clearance, about",
            "13 days, and is nearly the same whatever the regimen. Only bulk",
            "removal by lumbar puncture shortens it.")))
        ),

        ## ---- 4 ------------------------------------------------------
        tabPanel(
          "4 · Pressure & perfusion",
          h4("Intracranial pressure is a residual, not a state"),
          plotOutput("icpPlot", height = 320),
          fluidRow(
            column(6, plotOutput("routPlot", height = 280)),
            column(6, plotOutput("volPlot", height = 280))),
          h4("Pressure burden and drainage"),
          tableOutput("icpTab"),
          tags$p(tags$small(paste(
            "ICP = Pss + Pel(Vex + oedema) and absorption = Pel / Rout.",
            "Nothing in the model sets ICP directly, so pressure can only be",
            "changed by altering formation, resistance or volume. Antigen",
            "sets the resistance; the needle removes both volume and",
            "antigen.")))
        ),

        ## ---- 5 ------------------------------------------------------
        tabPanel(
          "5 · Immunity & IRIS",
          fluidRow(
            column(6, plotOutput("immPlot", height = 300)),
            column(6, plotOutput("cytPlot", height = 300))),
          h4("IRIS is the product of a rate and a stock"),
          plotOutput("irisPlot", height = 300),
          h4("ART timing sweep"),
          plotOutput("artSweep", height = 300),
          tags$p(tags$small(paste(
            "IRIS drive = (rate of immune recovery) x (antigen present when",
            "recovery starts). Deferring ART does not make the immune",
            "response smaller; it makes the antigen it meets smaller. That",
            "is why the CSF-leucocyte-poor patient, who clears antigen most",
            "slowly, is the one most harmed by early ART.")))
        ),

        ## ---- 6 ------------------------------------------------------
        tabPanel(
          "6 · Clinical endpoints",
          fluidRow(
            column(6, plotOutput("mortPlot", height = 320)),
            column(6, plotOutput("hazPlot", height = 320))),
          h4("Which clock is killing this patient?"),
          plotOutput("hazShare", height = 300),
          h4("Endpoint summary"),
          tableOutput("endTab")
        ),

        ## ---- 7 ------------------------------------------------------
        tabPanel(
          "7 · Scenario comparison",
          h4("All 16 regimens on the current virtual patient"),
          checkboxGroupInput("cmpSel", NULL, inline = TRUE,
                             choices = names(cm_scenarios),
                             selected = c("flu1200", "oral_acta", "ambd_5fc",
                                          "ambd_5fc_1wk", "ambition",
                                          "ambition_lp")),
          plotOutput("cmpBurden", height = 300),
          plotOutput("cmpICP", height = 280),
          plotOutput("cmpMort", height = 280),
          DTOutput("cmpTab")
        ),

        ## ---- 8 ------------------------------------------------------
        tabPanel(
          "8 · Safety",
          fluidRow(
            column(6, plotOutput("safeRenal", height = 300)),
            column(6, plotOutput("safeHaem", height = 300))),
          h4("The flucytosine clearance trap"),
          plotOutput("trapPlot", height = 300),
          tableOutput("safeTab"),
          tags$p(tags$small(paste(
            "Flucytosine is renally cleared and amphotericin is",
            "nephrotoxic, so the same flucytosine dose reaches a higher",
            "plasma concentration, and suppresses the marrow more, when it",
            "is given beside amphotericin than beside fluconazole. The",
            "interaction is a shared clearance organ, not additive marrow",
            "toxicity.")))
        ),

        ## ---- 9 ------------------------------------------------------
        tabPanel(
          "9 · Calibration",
          h4("Model against published trials"),
          DTOutput("calTab"),
          fluidRow(
            column(6, plotOutput("calEFA", height = 320)),
            column(6, plotOutput("calMort", height = 320))),
          h4("Known residuals, stated rather than hidden"),
          tags$ul(
            tags$li(paste("10-week mortality: mean absolute error ~2.1",
                          "percentage points across 9 regimens from 5",
                          "trials.")),
            tags$li(paste("2-week mortality is over-predicted by 3-4 points.",
                          "Cause identified: the burden hazard is",
                          "instantaneous rather than integrated over the",
                          "preceding days.")),
            tags$li(paste("Worst EFA residual: the 1-week AmB-d + 5FC arm,",
                          "-0.36 modelled against -0.42 published.")),
            tags$li(paste("The partner-drug hazard ratio is 0.85 against",
                          "ACTA's 0.62, so ergosterol antagonism explains",
                          "the direction but only about half the magnitude.")),
            tags$li(paste("The ART-timing curve has no optimum, because the",
                          "model contains no competing risk from leaving HIV",
                          "untreated. Do not read it beyond ~8 weeks."))
          )
        )
      )
    )
  )
)

## ---------------------------------------------------------------------
##  Server
## ---------------------------------------------------------------------
server <- function(input, output, session) {

  PRESET <- list(median = c(5.0, 25), high_burden = c(5.9, 25),
                 pauci = c(5.3, 12), immune = c(4.3, 150))

  observeEvent(input$phenotype, {
    v <- PRESET[[input$phenotype]]
    updateSliderInput(session, "logCFU", value = v[1])
    updateSliderInput(session, "cd4", value = v[2])
  })

  ## --- initial state: relax the slow variables at the chosen burden ---
  init_state <- reactive({
    cm_burnin(logCFU = input$logCFU, CD4 = input$cd4)
  })

  ## --- assemble the parameter list -----------------------------------
  par_list <- reactive({
    p <- cm_scenarios[[input$scenario]]$par
    if (isTRUE(input$custom)) {
      p <- list(KSUPP = 0.55)
      if (input$ambform == "d") {
        p$AMBD_MGKG <- input$ambdose; p$AMBD_DAYS <- input$ambdays
      } else if (input$ambform == "l") {
        p$AMBL_MGKG <- 10
      }
      if (isTRUE(input$use5fc)) {
        p$FC_MGKG <- 100; p$FC_DAYS <- input$fcdays
      }
      p$FLU1_MG <- input$flu1; p$FLU1_DAYS <- 14
      p$FLU2_MG <- input$flu2
    }
    p$WT      <- input$wt
    p$ART_DAY <- if (input$artday >= 200) 1e6 else input$artday
    p$LP_N    <- input$lpn
    p$LP_EVERY <- input$lpevery
    if (isTRUE(input$dex))  { p$DEX_MGKG <- 0.30; p$DEX_DAYS <- 42 }
    if (isTRUE(input$sert)) { p$SERT_MG <- 400;  p$SERT_DAYS <- 14 }
    if (isTRUE(input$ifng)) { p$IFN_UG <- 100;   p$IFN_N <- 6 }
    p
  })

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    y0 <- init_state()
    keep <- intersect(names(y0), names(init(mod)))
    p <- par_list(); p$INITMODE <- 0; p$CLAMPF <- 0
    mod %>%
      init(as.list(y0[keep])) %>%
      param(p) %>%
      mrgsim(end = input$days, delta = 0.25, hmax = 0.02) %>%
      as_tibble()
  })

  ## ---------------- tab 1 --------------------------------------------
  output$baseTab <- renderTable({
    y <- init_state()
    tibble(
      quantity = c("CSF burden (log10 CFU/mL)", "CSF GXM (ug/mL)",
                   "CSF CrAg titre (1:x)", "opening pressure (mmH2O)",
                   "outflow resistance (cmH2O/(mL/min))",
                   "CSF leucocytes (/uL)", "CSF IFN-gamma (pg/mL)",
                   "CD4 (/uL)", "haemoglobin (g/dL)"),
      value = c(sprintf("%.2f", log10(y$Fe + y$Fres + y$Ft + y$Fi)),
                sprintf("%.0f", y$GXM), sprintf("%.0f", 160 * y$GXM),
                sprintf("%.0f", y$ICP_o),
                sprintf("%.1f", y$Rout * 1440),
                sprintf("%.1f", y$WBC), sprintf("%.0f", y$IFNG),
                sprintf("%.0f", y$CD4), sprintf("%.1f", y$Hb)))
  }, colnames = FALSE)

  output$regTxt <- renderPrint({
    str(par_list(), give.attr = FALSE)
  })

  output$basePlot <- renderPlot({
    grid <- seq(3.5, 6.5, by = 0.25)
    d <- bind_rows(lapply(grid, function(g) {
      y <- cm_burnin(logCFU = g, CD4 = input$cd4)
      tibble(logCFU = g, GXM = y$GXM, ICP = y$ICP_o, WBC = y$WBC)
    }))
    d %>% pivot_longer(-logCFU) %>%
      ggplot(aes(logCFU, value)) +
      geom_line(linewidth = 1, colour = "#00695c") + geom_point() +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "presenting burden (log10 CFU/mL)", y = NULL,
           title = "Burden sets antigen, antigen sets pressure") + THEME
  })

  ## ---------------- tab 2 --------------------------------------------
  output$pkAmb <- renderPlot({
    sim() %>% transmute(time, `AmB-d plasma` = Ad / 30,
                        `L-AmB plasma` = Al / 5) %>%
      pivot_longer(-time) %>% filter(value > 1e-4) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_y_log10() +
      labs(x = "day", y = "ug/mL", colour = NULL,
           title = "Amphotericin plasma concentration") + THEME
  })
  output$pkOral <- renderPlot({
    sim() %>% transmute(time, `5FC plasma` = Cfc_o, `5FC CSF` = FCcsf,
                        `fluconazole plasma` = Cfl_o,
                        `fluconazole CSF` = FLcsf) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      labs(x = "day", y = "ug/mL", colour = NULL,
           title = "Oral partner drugs") + THEME
  })
  output$pkBrain <- renderPlot({
    sim() %>% ggplot(aes(time, Cbr_o)) +
      geom_line(linewidth = 1.1, colour = "#00695c") +
      geom_hline(yintercept = 0.15, linetype = 2, colour = "grey40") +
      annotate("text", x = Inf, y = 0.15, hjust = 1.05, vjust = -0.5,
               size = 3, label = "effect threshold") +
      labs(x = "day", y = "ug/mL",
           title = "CNS effect-site amphotericin") + THEME
  })
  output$pkKid <- renderPlot({
    sim() %>% ggplot(aes(time, Akid / 0.30)) +
      geom_line(linewidth = 1.1, colour = "#795548") +
      labs(x = "day", y = "ug/mL",
           title = "Renal cortical amphotericin (the toxicity driver)") + THEME
  })
  output$pkTab <- renderTable({
    s <- sim()
    tibble(quantity = c("CNS AmB AUC 0-14 (ug.d/mL)",
                        "renal AmB AUC 0-14 (ug.d/mL)",
                        "days with CNS AmB > 0.15 ug/mL",
                        "peak plasma 5FC (ug/mL)",
                        "CSF 5FC AUC 0-14 (ug.d/mL)",
                        "mean CSF fluconazole (ug/mL)",
                        "ergosterol nadir (fraction of normal)",
                        "peak effective AmB EC50 (ug/mL)"),
           value = c(sprintf("%.2f", approx(s$time, s$ABR, 14)$y),
                     sprintf("%.1f", sum(diff(s$time) *
                       head(s$Akid / 0.30, -1)) ),
                     sprintf("%.1f", sum(diff(s$time)[head(s$Cbr_o, -1) > .15])),
                     sprintf("%.1f", max(s$Cfc_o)),
                     sprintf("%.0f", approx(s$time, s$AFC, 14)$y),
                     sprintf("%.1f", mean(s$FLcsf[s$time <= 14])),
                     sprintf("%.3f", min(s$ERG)),
                     sprintf("%.2f", max(s$EC50Aeff))))
  }, colnames = FALSE)

  ## ---------------- tab 3 --------------------------------------------
  output$clocks <- renderPlot({
    s <- sim()
    b0 <- s$logCFU_o[1]; g0v <- s$GXM[1]
    tibble(time = s$time,
           `viable burden (clock 1)` = 10^(s$logCFU_o - b0),
           `capsular antigen (clock 2)` = s$GXM / g0v) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.3) +
      geom_hline(yintercept = 0.5, linetype = 3, colour = "grey40") +
      scale_y_log10(limits = c(1e-6, 3)) +
      scale_colour_manual(values = CLOCK_COLS) +
      labs(x = "day", y = "fraction of presenting value", colour = NULL,
           title = paste("The same patient, two clocks:",
                         "each normalised to its presenting value")) + THEME
  })
  output$clock1 <- renderPlot({
    sim() %>% transmute(time, extracellular = Fe, resistant = Fres,
                        persister = Ft, intracellular = Fi,
                        parenchymal = Fp) %>%
      pivot_longer(-time) %>% mutate(value = pmax(value, 1e-4)) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_y_log10() +
      labs(x = "day", y = "CFU/mL", colour = NULL) + THEME
  })
  output$clock2 <- renderPlot({
    sim() %>% transmute(time, `CSF GXM (ug/mL)` = GXM,
                        `outflow resistance (x normal)` = Rout / 0.005556) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 1, colour = "#e65100") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "day", y = NULL) + THEME
  })
  output$clockTab <- renderTable({
    s <- sim()
    st <- s$time[which(s$Sterile == 1)[1]]
    ghalf <- s$time[which(s$GXM <= 0.5 * s$GXM[1])[1]]
    ic <- s$time[which(s$ICP_o < 250)[1]]
    tibble(quantity = c("EFA over days 0-14 (log10 CFU/mL/day)",
                        "day of first negative CSF culture",
                        "CSF GXM on that day (ug/mL)",
                        "GXM as % of presenting value on that day",
                        "day CSF GXM has halved",
                        "day ICP falls below 250 mmH2O",
                        "lag: antigen half-life minus sterility (days)"),
           value = c(sprintf("%.3f", cm_efa(s)),
                     ifelse(is.na(st), "not reached", sprintf("%.1f", st)),
                     ifelse(is.na(st), "-",
                            sprintf("%.1f", approx(s$time, s$GXM, st)$y)),
                     ifelse(is.na(st), "-",
                            sprintf("%.0f%%", 100 * approx(s$time, s$GXM, st)$y /
                                      s$GXM[1])),
                     ifelse(is.na(ghalf), "not reached", sprintf("%.1f", ghalf)),
                     ifelse(is.na(ic), "not reached", sprintf("%.1f", ic)),
                     ifelse(is.na(st) || is.na(ghalf), "-",
                            sprintf("%.1f", ghalf - st))))
  }, colnames = FALSE)

  ## ---------------- tab 4 --------------------------------------------
  output$icpPlot <- renderPlot({
    sim() %>% ggplot(aes(time, ICP_o)) +
      geom_hline(yintercept = 250, linetype = 2, colour = "#c62828") +
      geom_line(linewidth = 1.2, colour = "#1565c0") +
      labs(x = "day", y = "mmH2O",
           title = "CSF opening pressure (dashed line = 250 mmH2O)") + THEME
  })
  output$routPlot <- renderPlot({
    sim() %>% ggplot(aes(GXM, Rout * 1440)) +
      geom_path(linewidth = 1, colour = "#e65100",
                arrow = arrow(length = unit(0.15, "cm"))) +
      labs(x = "CSF GXM (ug/mL)", y = "outflow resistance (cmH2O/(mL/min))",
           title = "Resistance follows antigen, with hysteresis") + THEME
  })
  output$volPlot <- renderPlot({
    sim() %>% transmute(time, `excess CSF volume (mL)` = Vex,
                        `oedema (mL-equivalent)` = EDEMA) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      labs(x = "day", y = "mL", colour = NULL,
           title = "The volume budget behind the pressure") + THEME
  })
  output$icpTab <- renderTable({
    s <- sim()
    tibble(quantity = c("peak ICP (mmH2O)", "day of peak ICP",
                        "pressure-time integral above 250 (mmH2O.day)",
                        "days spent above 250 mmH2O",
                        "minimum cerebral perfusion pressure (mmHg)",
                        "CSF volume drained (mL)",
                        "GXM removed by drainage (mg)"),
           value = c(sprintf("%.0f", max(s$ICP_o)),
                     sprintf("%.1f", s$time[which.max(s$ICP_o)]),
                     sprintf("%.0f", tail(s$AICP, 1)),
                     sprintf("%.1f", sum(diff(s$time)[head(s$ICP_o, -1) > 250])),
                     sprintf("%.0f", min(s$CPP_o)),
                     sprintf("%.0f", tail(s$CLPV, 1)),
                     sprintf("%.2f", tail(s$CLPG, 1) / 1000)))
  }, colnames = FALSE)

  ## ---------------- tab 5 --------------------------------------------
  output$immPlot <- renderPlot({
    sim() %>% transmute(time, `CD4 (cells/uL)` = CD4,
                        `CSF leucocytes (/uL)` = WBC,
                        `HIV RNA (log10)` = VL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 1, colour = "#6a1b9a") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "day", y = NULL) + THEME
  })
  output$cytPlot <- renderPlot({
    sim() %>% transmute(time, `IFN-gamma` = IFNG, `TNF/IL-6 index` = PROIN,
                        `IL-10` = IL10, `macrophage index x100` = 100 * MAC) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      labs(x = "day", y = "pg/mL (or index)", colour = NULL,
           title = "CSF cytokine network") + THEME
  })
  output$irisPlot <- renderPlot({
    s <- sim()
    tibble(time = s$time,
           `IRIS activity` = s$IRISa,
           `antigen stock (GXM/200)` = s$GXM / 200,
           `recovery rate (immcomp - lag) x10` =
             10 * pmax(0, s$CD4 / (s$CD4 + 110) - s$IMML)) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1.1) +
      labs(x = "day", y = NULL, colour = NULL,
           title = "IRIS = rate x stock. Neither factor alone predicts it.") +
      THEME
  })
  output$artSweep <- renderPlot({
    y0 <- init_state(); keep <- intersect(names(y0), names(init(mod)))
    d <- bind_rows(lapply(c(7, 10, 14, 21, 28, 35, 42, 56), function(a) {
      s <- mod %>% init(as.list(y0[keep])) %>%
        param(modifyList(par_list(),
                         list(ART_DAY = a, INITMODE = 0, CLAMPF = 0))) %>%
        mrgsim(end = 182, delta = 1, hmax = 0.05) %>% as_tibble()
      tibble(art_day = a, mortality = tail(s$Mortality, 1),
             iris_peak = max(s$IRISa),
             gxm_at_art = approx(s$time, s$GXM, a)$y)
    }))
    d %>% pivot_longer(-art_day) %>%
      ggplot(aes(art_day, value)) +
      geom_line(linewidth = 1.1, colour = "#ad1457") + geom_point() +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day of ART initiation", y = NULL,
           title = "26-week mortality, IRIS drive and the antigen ART meets") +
      THEME
  })

  ## ---------------- tab 6 --------------------------------------------
  output$mortPlot <- renderPlot({
    sim() %>% ggplot(aes(time, 100 * Mortality)) +
      geom_line(linewidth = 1.2, colour = "#b71c1c") +
      geom_vline(xintercept = c(14, 70), linetype = 3) +
      labs(x = "day", y = "cumulative mortality (%)",
           title = "Cumulative mortality (dashed: 2 and 10 weeks)") + THEME
  })
  output$hazPlot <- renderPlot({
    s <- sim()
    tibble(time = head(s$time, -1),
           hazard = diff(s$HAZ) / diff(s$time)) %>%
      ggplot(aes(time, hazard)) + geom_line(linewidth = 1, colour = "#e65100") +
      labs(x = "day", y = "instantaneous hazard (1/day)",
           title = "Hazard rate: front-loaded, as the trials show") + THEME
  })
  output$hazShare <- renderPlot({
    s <- sim(); p <- as.list(param(mod))
    haz <- tibble(
      burden      = p$hB * pmax(0, s$logCFU_o - 3),
      unsterile   = ifelse(s$Sterile == 0, p$hSTER, 0),
      pressure    = p$hICP * pmax(0, s$ICP_o - 250) / 250,
      perfusion   = p$hCPP * pmax(0, 60 - s$CPP_o) / 20,
      injury      = p$hAMS * s$NEUR * 10,
      IRIS        = p$hIRIS * s$IRISa,
      anaemia     = p$hAN * pmax(0, p$Hb0 - s$Hb),
      neutropenia = p$hNEUT * pmax(0, p$ANCthr - s$ANC),
      hypokalaemia= p$hK * pmax(0, p$Kthr - s$Kser),
      renal       = p$hGFR * pmax(0, p$GFRthr - s$GFR) / 30,
      background  = p$h0)
    tot <- vapply(haz, function(v) sum(head(v, -1) * diff(s$time)), 0)
    tibble(term = names(tot), share = 100 * tot / sum(tot)) %>%
      filter(share > 0.05) %>%
      ggplot(aes(reorder(term, share), share)) +
      geom_col(fill = "#37474f") + coord_flip() +
      labs(x = NULL, y = "% of cumulative hazard",
           title = "Separable hazard terms: which clock killed the patient") +
      THEME
  })
  output$endTab <- renderTable({
    s <- sim()
    st <- s$time[which(s$Sterile == 1)[1]]
    tibble(endpoint = c("EFA days 0-14 (log10 CFU/mL/day)",
                        "burden at day 14 (log10 CFU/mL)",
                        "day of first negative CSF culture",
                        "2-week mortality (%)", "10-week mortality (%)",
                        "mortality at end of horizon (%)",
                        "peak reversible injury index",
                        "permanent disability index",
                        "peak IRIS activity"),
           value = c(sprintf("%.3f", cm_efa(s)),
                     sprintf("%.2f", approx(s$time, s$logCFU_o, 14)$y),
                     ifelse(is.na(st), "not reached", sprintf("%.1f", st)),
                     sprintf("%.1f", 100 * approx(s$time, s$Mortality, 14)$y),
                     sprintf("%.1f", 100 * approx(s$time, s$Mortality,
                                                  min(70, max(s$time)))$y),
                     sprintf("%.1f", 100 * tail(s$Mortality, 1)),
                     sprintf("%.3f", max(s$NEUR)),
                     sprintf("%.3f", tail(s$DIS, 1)),
                     sprintf("%.3f", max(s$IRISa))))
  }, colnames = FALSE)

  ## ---------------- tab 7 --------------------------------------------
  cmp <- reactive({
    y0 <- init_state(); keep <- intersect(names(y0), names(init(mod)))
    bind_rows(lapply(input$cmpSel, function(nm) {
      p <- modifyList(cm_scenarios[[nm]]$par,
                      list(WT = input$wt, INITMODE = 0, CLAMPF = 0,
                           ART_DAY = if (input$artday >= 200) 1e6 else input$artday))
      mod %>% init(as.list(y0[keep])) %>% param(p) %>%
        mrgsim(end = input$days, delta = 0.5, hmax = 0.02) %>%
        as_tibble() %>% mutate(scenario = nm)
    }))
  })
  output$cmpBurden <- renderPlot({
    cmp() %>% ggplot(aes(time, logCFU_o, colour = scenario)) +
      geom_line(linewidth = 1) +
      labs(x = "day", y = "log10 CFU/mL", colour = NULL,
           title = "Clock 1 by regimen") + THEME
  })
  output$cmpICP <- renderPlot({
    cmp() %>% ggplot(aes(time, ICP_o, colour = scenario)) +
      geom_hline(yintercept = 250, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 1) +
      labs(x = "day", y = "mmH2O", colour = NULL,
           title = "Pressure by regimen: nearly the same curve") + THEME
  })
  output$cmpMort <- renderPlot({
    cmp() %>% ggplot(aes(time, 100 * Mortality, colour = scenario)) +
      geom_line(linewidth = 1) +
      labs(x = "day", y = "%", colour = NULL,
           title = "Cumulative mortality by regimen") + THEME
  })
  output$cmpTab <- renderDT({
    cmp() %>% group_by(scenario) %>%
      summarise(EFA = round(cm_efa(cur_data_all()), 3),
                d14_logCFU = round(approx(time, logCFU_o, 14)$y, 2),
                sterile_day = { i <- which(Sterile == 1)
                                if (length(i)) round(time[i[1]], 1) else NA },
                ICP_peak = round(max(ICP_o)),
                AUC_ICP_250 = round(max(AICP)),
                m2 = round(100 * approx(time, Mortality, 14)$y, 1),
                m10 = round(100 * approx(time, Mortality,
                                         min(70, max(time)))$y, 1),
                Hb_nadir = round(min(Hb), 2),
                ANC_nadir = round(min(ANC), 2),
                .groups = "drop") %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  ## ---------------- tab 8 --------------------------------------------
  output$safeRenal <- renderPlot({
    sim() %>% transmute(time, `GFR (mL/min)` = GFR,
                        `potassium (mmol/L)` = Kser) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 1, colour = "#795548") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "day", y = NULL) + THEME
  })
  output$safeHaem <- renderPlot({
    sim() %>% transmute(time, `haemoglobin (g/dL)` = Hb,
                        `neutrophils (1e9/L)` = ANC, `ALT (U/L)` = ALT) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 1, colour = "#5d4037") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "day", y = NULL) + THEME
  })
  output$trapPlot <- renderPlot({
    s <- sim()
    ggplot(s, aes(GFR, Cfc_o)) +
      geom_path(colour = "#00838f", linewidth = 1,
                arrow = arrow(length = unit(0.15, "cm"))) +
      labs(x = "GFR (mL/min)", y = "plasma flucytosine (ug/mL)",
           title = paste("Amphotericin lowers the organ that clears",
                         "flucytosine")) + THEME
  })
  output$safeTab <- renderTable({
    s <- sim()
    tibble(quantity = c("GFR nadir (mL/min)", "potassium nadir (mmol/L)",
                        "haemoglobin nadir (g/dL)",
                        "neutrophil nadir (1e9/L)", "peak ALT (U/L)",
                        "peak plasma 5FC (ug/mL)",
                        "days with 5FC above 100 ug/mL"),
           value = c(sprintf("%.0f", min(s$GFR)), sprintf("%.2f", min(s$Kser)),
                     sprintf("%.2f", min(s$Hb)), sprintf("%.2f", min(s$ANC)),
                     sprintf("%.0f", max(s$ALT)),
                     sprintf("%.1f", max(s$Cfc_o)),
                     sprintf("%.1f",
                             sum(diff(s$time)[head(s$Cfc_o, -1) > 100]))))
  }, colnames = FALSE)

  ## ---------------- tab 9 --------------------------------------------
  calib <- reactive({
    withProgress(message = "Running all calibration scenarios",
                 cm_calibration(days = 70))
  })
  output$calTab <- renderDT({
    calib() %>%
      transmute(scenario, source,
                EFA_model = round(efa_model, 3), EFA_published = efa_pub,
                EFA_residual = round(efa_model - efa_pub, 3),
                m10_model = round(100 * m10_model, 1),
                m10_published = round(100 * m10_pub, 1),
                m10_residual = round(100 * (m10_model - m10_pub), 1),
                sterile_day = round(sterile_day, 1)) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })
  output$calEFA <- renderPlot({
    calib() %>% ggplot(aes(efa_pub, efa_model, label = scenario)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      geom_point(size = 3, colour = "#2e7d32") +
      geom_text(hjust = -0.1, size = 3) +
      expand_limits(x = c(-0.6, 0), y = c(-0.6, 0)) +
      labs(x = "published EFA", y = "model EFA",
           title = "Early fungicidal activity") + THEME
  })
  output$calMort <- renderPlot({
    calib() %>% filter(!is.na(m10_pub)) %>%
      ggplot(aes(100 * m10_pub, 100 * m10_model, label = scenario)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      geom_point(size = 3, colour = "#b71c1c") +
      geom_text(hjust = -0.1, size = 3) +
      expand_limits(x = c(20, 65), y = c(20, 65)) +
      labs(x = "published 10-week mortality (%)",
           y = "model 10-week mortality (%)",
           title = "10-week mortality") + THEME
  })
}

shinyApp(ui, server)
