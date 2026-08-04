## =============================================================================
##  cmd_shiny_app.R
##  Coronary Microvascular Dysfunction (CMD / ANOCA-INOCA) — QSP dashboard
## =============================================================================
##
##  The app is built around one idea, and every tab is a different way of
##  refusing to collapse it:
##
##      CFR = v_hyper / v_rest is a RATIO.
##      The same low value comes from a raised DENOMINATOR (resting flow too
##      high — the functional endotype) or a floored NUMERATOR (hyperaemic flow
##      capped by structure — the structural endotype), and the two need
##      opposite drugs.
##
##  Tabs
##  ----
##   1. Patient & endotype      — set the physiology, see the endotype emerge
##   2. Coronary function test  — CFR, both resistances, IMR, MRR, ACh reserve
##   3. The ratio, opened up    — numerator and denominator plotted separately
##   4. Heart rate enters twice — demand vs the diastolic perfusion window
##   5. Transmural perfusion    — subendocardial deficit through the day
##   6. Drug PK / PD            — 14 agents, concentrations and effect handles
##   7. Endpoints               — SAQ, angina rate, Bruce time, hazards
##   8. Scenario comparison     — any two regimens side by side, 24 presets
##   9. Trial reproductions     — CorMicA, RWISE, PRIZE, WARRIOR
##  10. Virtual population      — prevalence, endotype split, response spread
##
##  Requires cmd_mrgsolve_model.R in the same directory.
##  ⚠ Educational / research tool.  Not for clinical use.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("cmd_mrgsolve_model.R")

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#eef2f6", colour = NA),
        legend.position = "bottom")

PAL <- c(control = "#2b8cbe", functional = "#e6550d",
         structural = "#756bb1", vasospastic = "#c51b8a",
         noncardiac = "#31a354")

## -----------------------------------------------------------------------------
##  UI
## -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Coronary Microvascular Dysfunction (ANOCA / INOCA) — QSP dashboard"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("CFR is a ratio: the same low number is a raised resting",
              "denominator (<b>functional</b>) or a floored hyperaemic numerator",
              "(<b>structural</b>). Supply is bought only in diastole, demand is",
              "spent only in systole &mdash; so heart rate enters the oxygen",
              "balance twice, with the same sign.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("endotype", "Endotype preset",
                  choices = names(cmd_endotypes), selected = "functional"),
      hr(),
      h5("Physiology (overrides the preset)"),
      sliderInput("rmin_f", "Structural floor RMIN_F (minimal resistance)",
                  min = 1.0, max = 2.0, value = 1.0, step = 0.02),
      sliderInput("auto_off", "Controller offset AUTO_OFF (functional lesion)",
                  min = 0.0, max = 1.5, value = 0.92, step = 0.02),
      sliderInput("infl", "Inflammatory drive INFL", min = 1.0, max = 3.0,
                  value = 1.4, step = 0.05),
      sliderInput("rock_d", "Rho-kinase drive", min = 1.0, max = 4.0,
                  value = 1.3, step = 0.05),
      sliderInput("spasm", "Spasm susceptibility", min = 0, max = 1,
                  value = 0.10, step = 0.05),
      sliderInput("sens_bas", "Primary central sensitisation",
                  min = 0, max = 2, value = 0.75, step = 0.05),
      checkboxInput("geno", "rs9349379-G carrier (raised ET-1)", TRUE),
      hr(),
      h5("Therapy"),
      selectInput("regimen", "Regimen",
                  choices = c("untreated", "ivabradine", "bisoprolol",
                              "nebivolol", "ranolazine", "amlodipine",
                              "diltiazem", "nicorandil", "zibotentan",
                              "zibotentan_cf", "zib_nobp", "zib_nofluid",
                              "fasudil", "sildenafil", "aminophylline",
                              "imipramine", "trimetazidine", "rehab",
                              "imt_warrior", "iva_ran", "cormica_func",
                              "cormica_struct", "cormica_vaso",
                              "cormica_noncard"),
                  selected = "untreated"),
      selectInput("regimen2", "Comparator regimen (tab 8)",
                  choices = c("untreated", "ivabradine", "ranolazine",
                              "amlodipine", "nicorandil", "zibotentan",
                              "aminophylline", "imipramine", "rehab",
                              "cormica_func"),
                  selected = "ivabradine"),
      sliderInput("days", "Follow-up (days)", min = 14, max = 730,
                  value = 168, step = 14),
      hr(),
      actionButton("run", "Run", class = "btn-primary"),
      tags$p(style = "font-size:11px;color:#777;margin-top:12px",
             "Calibrated against Rahman 2019 (endotypes), Kelshiker 2022",
             "(hazard slope), CorMicA, RWISE, PRIZE and WARRIOR. All numbers",
             "in the reference output are produced by cmd_reference_model.py.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. Patient & endotype",
                 h4("Emergent endotype call"),
                 verbatimTextOutput("call"),
                 h4("What the sliders did to the physiology"),
                 tableOutput("phys"),
                 plotOutput("statePlot", height = "420px")),
        tabPanel("2. Coronary function test",
                 h4("The invasive study, as performed"),
                 tableOutput("cft"),
                 plotOutput("cftPlot", height = "380px"),
                 tags$p(style = "color:#555",
                        HTML("Adenosine abolishes the metabolic and myogenic",
                             "terms and 88% of constrictor tone (calibrated);",
                             "the residue is why an endothelin or Rho-kinase",
                             "antagonist can move a <i>hyperaemic</i>",
                             "measurement at all."))),
        tabPanel("3. The ratio, opened up",
                 h4("Numerator and denominator, never merged"),
                 plotOutput("ratioPlot", height = "460px"),
                 tableOutput("ratioTab"),
                 tags$p(style = "color:#555",
                        HTML("Two patients with the same CFR can differ by",
                             "~57% in absolute maximal flow. A ratio cannot see",
                             "that; an absolute flow can. Note also that",
                             "microvascular resistance reserve does NOT rescue",
                             "this &mdash; a drug lowers pressure in both states,",
                             "so the correction cancels and MRR tracks CFR."))),
        tabPanel("4. Heart rate enters twice",
                 h4("Demand and the diastolic perfusion window"),
                 sliderInput("hrwl", "Fixed workload for the sweep",
                             min = 1.0, max = 4.0, value = 3.2, step = 0.1),
                 plotOutput("hrPlot", height = "440px"),
                 tableOutput("hrTab"),
                 tags$p(style = "color:#555",
                        HTML("43% of the fall in subendocardial deficit from",
                             "68 &rarr; 55 bpm is the diastolic time window and",
                             "57% is lower demand, so a pure rate-slowing drug",
                             "has a supply-side action and is not",
                             "interchangeable with removing demand another way."))),
        tabPanel("5. Transmural perfusion",
                 h4("Subendocardial deficit across the day"),
                 plotOutput("transPlot", height = "460px"),
                 tags$p(style = "color:#555",
                        HTML("The workload profile is a real diurnal covariate:",
                             "15 h at rest, 6 h light, 2.5 h moderate, 0.5 h",
                             "stair-climbing. Exercise filling pressure is",
                             "subtracted from subendocardial driving pressure,",
                             "which is where the CMD/HFpEF overlap lives."))),
        tabPanel("6. Drug PK / PD",
                 h4("Concentrations and the handles they pull"),
                 plotOutput("pkPlot", height = "340px"),
                 plotOutput("pdPlot", height = "340px")),
        tabPanel("7. Endpoints",
                 h4("What the patient and the trial actually see"),
                 plotOutput("epPlot", height = "460px"),
                 tableOutput("epTab")),
        tabPanel("8. Scenario comparison",
                 h4("Two regimens, same patient"),
                 plotOutput("cmpPlot", height = "460px"),
                 tableOutput("cmpTab"),
                 h4("All 24 preset scenarios"),
                 tableOutput("allTab")),
        tabPanel("9. Trial reproductions",
                 h4("CorMicA — stratified vs unguided care"),
                 tableOutput("cormicaTab"),
                 verbatimTextOutput("cormicaTxt"),
                 h4("RWISE — a real subgroup effect, diluted away"),
                 tableOutput("rwiseTab"),
                 verbatimTextOutput("rwiseTxt"),
                 h4("WARRIOR — contamination and timescale"),
                 tableOutput("warriorTab"),
                 verbatimTextOutput("warriorTxt")),
        tabPanel("10. Virtual population",
                 h4("Prevalence and the endotype split"),
                 numericInput("npop", "Population size", 120, min = 20,
                              max = 500, step = 20),
                 actionButton("runpop", "Sample population"),
                 plotOutput("popPlot", height = "420px"),
                 tableOutput("popTab"),
                 tags$p(style = "color:#555",
                        HTML("Risk factors are sampled first and the structural",
                             "floor and controller offset drawn conditionally on",
                             "them; independent draws produce patients with",
                             "severe remodelling and no risk factor at all.")))
      )
    )
  )
)

## -----------------------------------------------------------------------------
##  SERVER
## -----------------------------------------------------------------------------
server <- function(input, output, session) {

  ## keep the sliders in step with the preset
  observeEvent(input$endotype, {
    p <- cmd_endotypes[[input$endotype]]
    updateSliderInput(session, "rmin_f", value = p$RMIN_F)
    updateSliderInput(session, "auto_off", value = p$AUTO_OFF)
    updateSliderInput(session, "infl", value = p$INFL)
    updateSliderInput(session, "rock_d", value = p$ROCK_D)
    updateSliderInput(session, "spasm", value = p$SPASM)
    updateSliderInput(session, "sens_bas", value = p$SENS_BAS)
    updateCheckboxInput(session, "geno", value = p$GENO > 0.5)
  })

  overrides <- reactive({
    list(RMIN_F = input$rmin_f, AUTO_OFF = input$auto_off, INFL = input$infl,
         ROCK_D = input$rock_d, SPASM = input$spasm,
         SENS_BAS = input$sens_bas, GENO = as.numeric(input$geno),
         KSBP = if (input$rmin_f > 1.15) 37 else 19)
  })

  sim <- eventReactive(input$run, {
    cmd_simulate(input$endotype, input$regimen, days = input$days,
                 delta = 1, extra = overrides())
  }, ignoreNULL = FALSE)

  sim2 <- eventReactive(input$run, {
    cmd_simulate(input$endotype, input$regimen2, days = input$days,
                 delta = 1, extra = overrides())
  }, ignoreNULL = FALSE)

  last <- reactive({ s <- sim(); s[nrow(s), ] })

  ## ---- tab 1 --------------------------------------------------------------
  output$call <- renderText({
    r <- last()
    lab <- if (r$CFR >= 2.5) {
      if (r$ACHFR < 0.9) "vasospastic angina (normal CFR, positive ACh test)"
      else "no coronary microvascular dysfunction"
    } else if (r$MR_HYP < 2.5) {
      "CMD, FUNCTIONAL endotype: normal minimal resistance, raised resting flow"
    } else {
      "CMD, STRUCTURAL endotype: raised minimal resistance"
    }
    paste0(lab, "\n",
           sprintf("CFR %.2f   resting MR %.2f   hyperaemic MR %.2f   IMR %.1f U\n",
                   r$CFR, r$MR_RST, r$MR_HYP, r$IMR),
           sprintf("resting MBF %.2f   hyperaemic MBF %.2f mL/min/g   ",
                   r$MBF_RST, r$MBF_HYP),
           sprintf("resting O2 extraction %.2f\n", r$E_RST),
           sprintf("SAQ %.1f   angina %.1f/week   hyperaemic endo/epi %.2f",
                   r$SAQ, r$ANG, r$ENDOEPI))
  })

  output$phys <- renderTable({
    r <- last()
    data.frame(
      quantity = c("minimal resistance RMIN x structure", "resting arteriolar tone",
                   "hyperaemic subendo/subepi ratio", "LV end-diastolic pressure",
                   "24 h subendocardial deficit", "24 h afferent drive",
                   "central sensitisation", "sympathetic tone"),
      value = c(sprintf("%.2f mmHg/(cm/s)", r$MR_HYP),
                sprintf("%.3f", 0.5 * (r$TONE_E + r$TONE_P)),
                sprintf("%.3f", r$ENDOEPI),
                sprintf("%.1f mmHg", r$LVEDP),
                sprintf("%.3f mL O2/min/100 g", r$BI),
                sprintf("%.3f", r$BN),
                sprintf("%.2f", r$SENS),
                sprintf("%.2f", r$SYMP)),
      stringsAsFactors = FALSE)
  }, rownames = FALSE)

  output$statePlot <- renderPlot({
    s <- sim()
    d <- s %>%
      select(day, NO, ROS, ET1, ROCK, ML, CAPD, LVEDP, CAD) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(colour = "#2b8cbe", linewidth = 0.7) +
      facet_wrap(~ name, scales = "free_y", ncol = 4) +
      labs(x = "day", y = NULL, title = "Mechanistic states") + THEME
  })

  ## ---- tab 2 --------------------------------------------------------------
  output$cft <- renderTable({
    r <- last()
    data.frame(
      index = c("CFR", "resting MR", "hyperaemic MR", "IMR", "MRR",
                "ACh flow reserve", "hyperaemic endo/epi", "DPTI", "SPTI",
                "SEVR"),
      value = c(sprintf("%.2f", r$CFR), sprintf("%.2f", r$MR_RST),
                sprintf("%.2f", r$MR_HYP), sprintf("%.1f U", r$IMR),
                sprintf("%.2f", r$MRR), sprintf("%.2f", r$ACHFR),
                sprintf("%.2f", r$ENDOEPI), sprintf("%.1f", r$DPTI),
                sprintf("%.1f", r$SPTI), sprintf("%.2f", r$SEVR)),
      threshold = c("< 2.0 or < 2.5", "—", ">= 2.5 = structural", ">= 25",
                    "< 3.0", "< 0.9 abnormal", "< 1.0 = ischaemia", "—", "—",
                    "—"),
      stringsAsFactors = FALSE)
  }, rownames = FALSE)

  output$cftPlot <- renderPlot({
    s <- sim()
    d <- s %>% select(day, CFR, MR_RST, MR_HYP, IMR) %>% pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(linewidth = 0.8, colour = "#756bb1") +
      facet_wrap(~ name, scales = "free_y") +
      labs(x = "day", y = NULL, title = "Indices over time") + THEME
  })

  ## ---- tab 3 --------------------------------------------------------------
  output$ratioPlot <- renderPlot({
    tab <- do.call(rbind, lapply(names(cmd_endotypes), function(k) {
      s <- cmd_simulate(k, "untreated", days = 2, delta = 1)
      r <- s[nrow(s), ]
      data.frame(endotype = k, MBF_rest = r$MBF_RST, MBF_hyp = r$MBF_HYP,
                 CFR = r$CFR, stringsAsFactors = FALSE)
    }))
    d <- tab %>% pivot_longer(c(MBF_rest, MBF_hyp))
    ggplot(d, aes(endotype, value, fill = name)) +
      geom_col(position = "dodge") +
      geom_text(data = tab, aes(endotype, MBF_hyp + 0.18,
                                label = sprintf("CFR %.2f", CFR)),
                inherit.aes = FALSE, size = 4) +
      scale_fill_manual(values = c(MBF_rest = "#9ecae1", MBF_hyp = "#08519c"),
                        name = NULL) +
      labs(x = NULL, y = "myocardial blood flow (mL/min/g)",
           title = "The same ratio, different absolute flows") + THEME
  })

  output$ratioTab <- renderTable({
    do.call(rbind, lapply(names(cmd_endotypes), function(k) {
      s <- cmd_simulate(k, "untreated", days = 2, delta = 1)
      r <- s[nrow(s), ]
      data.frame(endotype = k, CFR = round(r$CFR, 2),
                 MR_rest = round(r$MR_RST, 2), MR_hyp = round(r$MR_HYP, 2),
                 MRR = round(r$MRR, 2), IMR = round(r$IMR, 1),
                 MBF_rest = round(r$MBF_RST, 2), MBF_hyp = round(r$MBF_HYP, 2),
                 E_rest = round(r$E_RST, 3), stringsAsFactors = FALSE)
    }))
  }, rownames = FALSE)

  ## ---- tab 4 --------------------------------------------------------------
  hrsweep <- reactive({
    hrs <- c(50, 55, 60, 68, 75, 85, 95, 110)
    do.call(rbind, lapply(hrs, function(hr) {
      par <- c(cmd_endotypes[[input$endotype]], overrides(),
               list(WL_FIX = input$hrwl, HR0 = hr))
      s <- param(cmd_mod, par) %>%
        mrgsim(end = 48, delta = 2, hmax = 0.05) %>% as.data.frame()
      r <- s[nrow(s), ]
      data.frame(HR = hr, DPTI = r$DPTI, SEVR = r$SEVR, deficit = r$BI,
                 afferent = r$BN, endo_epi = r$ENDOEPI,
                 stringsAsFactors = FALSE)
    }))
  })

  output$hrPlot <- renderPlot({
    d <- hrsweep() %>% pivot_longer(-HR)
    ggplot(d, aes(HR, value)) + geom_line(linewidth = 0.9, colour = "#e6550d") +
      geom_point() + facet_wrap(~ name, scales = "free_y") +
      labs(x = "resting heart rate (bpm)", y = NULL,
           title = sprintf("Heart-rate sweep at fixed workload %.1f", input$hrwl)) +
      THEME
  })
  output$hrTab <- renderTable({ hrsweep() %>% mutate(across(-HR, ~round(., 3))) },
                              rownames = FALSE)

  ## ---- tab 5 --------------------------------------------------------------
  output$transPlot <- renderPlot({
    s <- sim()
    last3 <- s %>% filter(day >= max(day) - 3)
    d <- last3 %>% select(day, ENDOEPI, DPTI, SEVR, BI, BN) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(colour = "#31a354", linewidth = 0.7) +
      facet_wrap(~ name, scales = "free_y", ncol = 2) +
      labs(x = "day", y = NULL,
           title = "Last three days: the diurnal workload profile at work") +
      THEME
  })

  ## ---- tab 6 --------------------------------------------------------------
  output$pkPlot <- renderPlot({
    s <- sim()
    cols <- grep("^C_|^M_", names(s), value = TRUE)
    d <- s %>% select(day, all_of(cols)) %>% pivot_longer(-day) %>%
      group_by(name) %>% filter(max(value) > 1e-9) %>% ungroup()
    if (!nrow(d)) return(ggplot() + labs(title = "no drug on board") + THEME)
    ggplot(d, aes(day, value, colour = name)) + geom_line() +
      facet_wrap(~ name, scales = "free_y") +
      labs(x = "day", y = "concentration (mg/L)",
           title = "Pharmacokinetics") + THEME + theme(legend.position = "none")
  })

  output$pdPlot <- renderPlot({
    s <- sim()
    d <- s %>% select(day, HR_RST, MR_HYP, LVEDP, CAD, NO, ROCK) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(colour = "#c51b8a", linewidth = 0.7) +
      facet_wrap(~ name, scales = "free_y") +
      labs(x = "day", y = NULL, title = "The handles the drugs pull") + THEME
  })

  ## ---- tab 7 --------------------------------------------------------------
  output$epPlot <- renderPlot({
    s <- sim()
    d <- s %>% select(day, SAQ, ANG, SENS, TNI, BNP, CHHOSP) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(colour = "#08519c", linewidth = 0.8) +
      facet_wrap(~ name, scales = "free_y") +
      labs(x = "day", y = NULL, title = "Clinical endpoints") + THEME
  })

  output$epTab <- renderTable({
    r <- last()
    data.frame(endpoint = c("SAQ summary", "angina episodes/week",
                            "central sensitisation", "hs-cTnI (ng/L)",
                            "NT-proBNP (pg/mL)",
                            "cumulative hospitalisation hazard",
                            "cumulative mortality hazard"),
               value = round(c(r$SAQ, r$ANG, r$SENS, r$TNI, r$BNP,
                               r$CHHOSP, r$CHMORT), 4),
               stringsAsFactors = FALSE)
  }, rownames = FALSE)

  ## ---- tab 8 --------------------------------------------------------------
  output$cmpPlot <- renderPlot({
    d <- rbind(sim() %>% mutate(arm = input$regimen),
               sim2() %>% mutate(arm = input$regimen2)) %>%
      select(day, arm, SAQ, ANG, CFR, MR_HYP, BI, BN) %>%
      pivot_longer(c(SAQ, ANG, CFR, MR_HYP, BI, BN))
    ggplot(d, aes(day, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~ name, scales = "free_y") +
      scale_colour_manual(values = c("#e6550d", "#2b8cbe"), name = NULL) +
      labs(x = "day", y = NULL, title = "Head to head") + THEME
  })

  output$cmpTab <- renderTable({
    a <- last(); b <- sim2()[nrow(sim2()), ]
    data.frame(quantity = c("CFR", "hyperaemic MR", "24 h deficit",
                            "24 h afferent drive", "angina/week", "SAQ"),
               arm1 = round(c(a$CFR, a$MR_HYP, a$BI, a$BN, a$ANG, a$SAQ), 3),
               arm2 = round(c(b$CFR, b$MR_HYP, b$BI, b$BN, b$ANG, b$SAQ), 3),
               difference = round(c(b$CFR - a$CFR, b$MR_HYP - a$MR_HYP,
                                    b$BI - a$BI, b$BN - a$BN, b$ANG - a$ANG,
                                    b$SAQ - a$SAQ), 3),
               stringsAsFactors = FALSE)
  }, rownames = FALSE)

  output$allTab <- renderTable({
    cmd_run_all(days = 168) %>%
      mutate(across(where(is.numeric), ~round(., 3)))
  }, rownames = FALSE)

  ## ---- tab 9 --------------------------------------------------------------
  cormica <- reactive(cmd_cormica())
  output$cormicaTab <- renderTable({
    cormica() %>% mutate(across(where(is.numeric), ~round(., 2)))
  }, rownames = FALSE)
  output$cormicaTxt <- renderText({
    sprintf(paste("population mean dSAQ %+.2f U (CorMicA observed +11.70 U,",
                  "95%% CI 5.0-18.4).\nThe stratified arm wins for a mechanical",
                  "reason: the same prescription is right for one endotype and",
                  "wrong\nfor another, so an unstratified arm averages a benefit",
                  "with a harm."),
            attr(cormica(), "population_dSAQ"))
  })

  rwise <- reactive(cmd_rwise())
  output$rwiseTab <- renderTable({
    rwise() %>% mutate(across(where(is.numeric), ~round(., 3)))
  }, rownames = FALSE)
  output$rwiseTxt <- renderText({
    sprintf(paste("CFR<2.5 stratum %+.2f U; whole randomised cohort %+.2f U.",
                  "\nNothing about the drug changed between those two numbers -",
                  "only the prevalence of the mechanism it targets."),
            attr(rwise(), "dSAQ_CMD_stratum"),
            attr(rwise(), "dSAQ_whole_cohort"))
  })

  warrior <- reactive(cmd_warrior())
  output$warriorTab <- renderTable({
    warrior()$observable %>% mutate(HR_observed = round(HR_observed, 3))
  }, rownames = FALSE)
  output$warriorTxt <- renderText({
    w <- warrior()
    sprintf(paste("model true hazard ratio %.3f; CFR gain %+.3f over 2.5 y but",
                  "SAQ gain only %+.2f U.\nWARRIOR observed HR 1.13 (0.94-1.37)",
                  "with a contamination-adjusted estimate of 0.74.\nA structural",
                  "target scored on a symptom endpoint is a mismatch of",
                  "instruments."),
            w$hr_true, w$dCFR, w$dSAQ)
  })

  ## ---- tab 10 -------------------------------------------------------------
  pop <- eventReactive(input$runpop, {
    cmd_population(n = input$npop, days = 84)
  })

  output$popPlot <- renderPlot({
    p <- pop()
    p$stratum <- ifelse(p$CFR >= 2.5, "CFR >= 2.5",
                        ifelse(p$MR_hyp < 2.5, "functional", "structural"))
    ggplot(p, aes(MR_hyp, MR_rest, colour = stratum)) +
      geom_point(size = 2, alpha = 0.8) +
      geom_vline(xintercept = 2.5, linetype = 2) +
      labs(x = "hyperaemic MR (mmHg/(cm/s))",
           y = "resting MR (mmHg/(cm/s))",
           title = "The plane the endotype call is made in") + THEME
  })

  output$popTab <- renderTable({
    p <- pop()
    cmd <- p[p$CFR < 2.5, ]
    data.frame(quantity = c("n", "CFR < 2.5", "of those, functional",
                            "of those, structural", "CFR < 2.0", "IMR >= 25"),
               model = c(nrow(p),
                         sprintf("%.0f%%", 100 * nrow(cmd) / nrow(p)),
                         sprintf("%.0f%%", 100 * sum(cmd$MR_hyp < 2.5) /
                                   max(1, nrow(cmd))),
                         sprintf("%.0f%%", 100 * sum(cmd$MR_hyp >= 2.5) /
                                   max(1, nrow(cmd))),
                         sprintf("%.0f%%", 100 * mean(p$CFR < 2.0)),
                         sprintf("%.0f%%", 100 * mean(p$IMR >= 25))),
               observed = c("—", "53% (Rahman 2019)", "62%", "38%", "—", "—"),
               stringsAsFactors = FALSE)
  }, rownames = FALSE)
}

shinyApp(ui, server)
