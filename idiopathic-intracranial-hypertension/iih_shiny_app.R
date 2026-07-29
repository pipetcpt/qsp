# =============================================================================
# Idiopathic Intracranial Hypertension -- QSP explorer (Shiny)
# =============================================================================
#
# This dashboard is deliberately built around ONE idea: intracranial pressure is
# the fixed point of a loop, so the app never shows a pressure trace without
# also showing (a) the loop gain that produced it and (b) the patient's own
# untreated counterfactual.  A number of the tabs exist specifically to make a
# clinically misleading measurement look misleading:
#
#   Tab 3  splits the OCT RNFL number into the two states it is a sum of.
#   Tab 5  shows that pulsatile tinnitus tracks the venous gradient and not the
#          pressure, so a drug and a stent that lower ICP equally do not
#          abolish it equally.
#   Tab 6  shows that headache carries a sensitisation state with its own
#          hysteresis, so normalising pressure has a ceiling.
#   Tab 8  shows the lumbar puncture emptying and refilling in half an hour
#          while the symptom states move for weeks.
#
# Run with:   shiny::runApp("iih_shiny_app.R")
# Requires:   shiny, mrgsolve, dplyr, tidyr, ggplot2, DT, patchwork
# and the model file in the same directory:
source("iih_mrgsolve_model.R")

library(shiny)
library(tidyr)
library(DT)
library(patchwork)

CM <- 1.35951
THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        strip.background = element_rect(fill = "#eef2f6", colour = NA),
        plot.title = element_text(face = "bold", size = 12))

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Idiopathic intracranial hypertension - QSP explorer"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("ICP = P<sub>cv</sub> + Q<sub>v</sub>(R<sub>ts0</sub> +
               R<sub>tsf</sub> + R<sub>tsc</sub>(ICP)) +
               R<sub>out</sub>F<sub>form</sub> &mdash; a fixed point, not an
               accumulation. Loop gain g, amplification 1/(1&minus;g).")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Patient"),
      sliderInput("bmi", "BMI at presentation", 24, 52, 38, step = 0.5),
      sliderInput("pcrit", HTML("Sinus wall support PCRIT (mmHg)<br/>
                                 <small>lower = more collapsible = higher
                                 loop gain</small>"),
                  8.0, 16.0, 10.95, step = 0.05),
      sliderInput("pstiff", HTML("Collapse sharpness PSTIFF (mmHg)<br/>
                                  <small>&lt;1.5 is the bistable,
                                  potentially fulminant regime</small>"),
                  0.6, 6.0, 4.14, step = 0.02),
      sliderInput("runin", "Untreated run-in before therapy (days)",
                  0, 365, 120, step = 5),
      sliderInput("analg", "Analgesic exposure (0-1)", 0, 1, 0.25, step = 0.05),

      hr(), h4("Drugs"),
      sliderInput("acz", "Acetazolamide (mg/day)", 0, 4000, 0, step = 250),
      sliderInput("tpm", "Topiramate (mg/day)", 0, 400, 0, step = 25),
      sliderInput("ex",  "Exenatide (ug/day)", 0, 40, 0, step = 5),
      sliderInput("sem", "Semaglutide (mg/week)", 0, 2.4, 0, step = 0.2),
      sliderInput("fur", "Furosemide (mg/day)", 0, 160, 0, step = 20),
      sliderInput("azd", "AZD4017 (mg/day)", 0, 800, 0, step = 100),
      sliderInput("pred", "Prednisolone (mg/day)", 0, 60, 0, step = 5),
      checkboxInput("citrate", "Potassium citrate cover", FALSE),
      sliderInput("tstop", "Stop all drugs on day", 0, 730, 730, step = 10),

      hr(), h4("Weight and procedures"),
      sliderInput("diet", "Weight programme target (fraction)",
                  0, 0.30, 0, step = 0.01),
      sliderInput("bar", "Bariatric surgery target (fraction)",
                  0, 0.35, 0, step = 0.01),
      numericInput("tbar", "  day of surgery", 1e9),
      numericInput("tstent", "Venous sinus stent on day (1e9 = never)", 1e9),
      numericInput("tshunt", "CSF shunt on day", 1e9),
      numericInput("tonsf", "ONSF on day", 1e9),
      numericInput("tlp", "Lumbar puncture on day", 1e9),
      sliderInput("days", "Simulation horizon (days)", 90, 1095, 365,
                  step = 15),
      actionButton("go", "Simulate", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1. Patient and the loop",
                 br(), fluidRow(column(6, verbatimTextOutput("fixpoint")),
                                column(6, plotOutput("gcurve", height = 300))),
                 h5("Where this patient's pressure comes from, term by term"),
                 plotOutput("decomp", height = 320),
                 DTOutput("decomptab")),

        tabPanel("2. Pressure and the venous gradient",
                 br(), plotOutput("icp", height = 380),
                 h5("The gradient is the amplifier, and it is measurable"),
                 plotOutput("grad", height = 300)),

        tabPanel("3. Optic disc: the OCT trap",
                 br(),
                 tags$p(HTML("<b>RNFL(measured) = surviving axons + prelaminar
                              swelling.</b> Swelling needs living axons to
                              swell, so a recovering disc and a dying disc pass
                              through the same measured thickness. The dashed
                              line is what OCT reports; the solid lines are the
                              two states it is a sum of.")),
                 plotOutput("onh", height = 400),
                 plotOutput("tlpd", height = 300)),

        tabPanel("4. Visual endpoints",
                 br(), plotOutput("vision", height = 400),
                 h5("Injury is an integral of supra-threshold gradient, so"),
                 h5("time-to-treatment competes with dose"),
                 DTOutput("delaytab")),

        tabPanel("5. Gradient-driven symptoms",
                 br(),
                 tags$p(HTML("Pulsatile tinnitus reads <b>Q<sub>v</sub>R<sub>ts</sub></b>,
                              not ICP. Diplopia reads the pressure peak. They
                              therefore dissociate under treatments that lower
                              pressure by the same amount.")),
                 plotOutput("gradsym", height = 380),
                 DTOutput("dissoctab")),

        tabPanel("6. Headache",
                 br(),
                 tags$p(HTML("Central sensitisation has its own hysteresis and
                              medication overuse blocks its decay, so the SAME
                              normalised pressure leaves a different headache
                              burden.")),
                 plotOutput("headache", height = 400),
                 DTOutput("hatab")),

        tabPanel("7. Drug exposure, tolerability, adherence",
                 br(), plotOutput("pk", height = 380),
                 h5("The dose-response is bent by adherence, not by the target"),
                 plotOutput("dr", height = 320),
                 DTOutput("drtab")),

        tabPanel("8. Lumbar puncture",
                 br(),
                 tags$p(HTML("&tau; = C&middot;R<sub>out</sub>. The pressure is
                              back within a cmH<sub>2</sub>O in under an hour.
                              Whatever relief a diagnostic LP gives for days is
                              not a pressure effect.")),
                 plotOutput("lpfast", height = 320),
                 plotOutput("lpslow", height = 320)),

        tabPanel("9. Scenario comparison",
                 br(), selectInput("scen", "Scenarios", multiple = TRUE,
                                   choices = names(SCENARIOS),
                                   selected = c("s01_untreated", "s03_acz_std",
                                                "s07_diet", "s09_exenatide",
                                                "s13_stent")),
                 plotOutput("scenplot", height = 460),
                 DTOutput("scentab")),

        tabPanel("10. Susceptibility and folds",
                 br(),
                 tags$p(HTML("Only <b>PCRIT</b> and <b>PSTIFF</b> change down
                              these tables. Nothing pharmacological does. The
                              responder / non-responder dichotomy and the
                              fulminant phenotype both live here.")),
                 plotOutput("spectrum", height = 360),
                 DTOutput("spectab"),
                 h5("Continuation: how many stable pressures does this patient have?"),
                 plotOutput("contplot", height = 320)),

        tabPanel("11. Calibration and biases",
                 br(), h5("Simulated versus published"),
                 DTOutput("caltab"),
                 h5("Known biases, kept rather than fitted"),
                 verbatimTextOutput("biases"))
      )
    )
  )
)

# -----------------------------------------------------------------------------
# server
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  pars <- reactive(list(PCRIT = input$pcrit, PSTIFF = input$pstiff))

  patient <- reactive({
    iih_patient(bmi = input$bmi, pars = pars())
  })

  regimen_pars <- reactive({
    list(DOSACZ = input$acz, DOSTPM = input$tpm, DOSEX = input$ex,
         DOSSEM = input$sem, DOSFUR = input$fur, DOSAZD = input$azd,
         DOSPRED = input$pred, CITRATE = as.numeric(input$citrate),
         FDIET = input$diet, TDIET = 0, FBAR = input$bar, TBAR = input$tbar,
         TSTENT = input$tstent, TSHUNT = input$tshunt, TONSF = input$tonsf,
         TLP = input$tlp, ANALG = input$analg, TSTOP = input$tstop)
  })

  sim <- eventReactive(input$go, {
    m <- patient() %>% param(ANALG = input$analg)
    treated <- do.call(iih_run, c(list(m, days = input$days, delta = 1,
                                       runin = input$runin),
                                  regimen_pars()))
    ctrl <- iih_run(m, days = input$days, delta = 1, runin = input$runin,
                    ANALG = input$analg)
    bind_rows(as.data.frame(treated) %>% mutate(arm = "treated"),
              as.data.frame(ctrl)    %>% mutate(arm = "untreated counterfactual"))
  }, ignoreNULL = FALSE)

  # ---- tab 1 ---------------------------------------------------------------
  output$fixpoint <- renderPrint({
    s <- attr(patient(), "steady")
    cat(sprintf(
      "Untreated fixed point\n---------------------\n"))
    cat(sprintf("stable pressures found      : %d\n", s$nstable))
    cat(sprintf("opening pressure            : %.1f cmH2O  (%.2f mmHg)\n",
                s$icp_cmH2O, s$icp))
    cat(sprintf("central venous pressure     : %.2f mmHg\n", s$pcv))
    cat(sprintf("transverse sinus gradient   : %.2f mmHg\n", s$grad))
    cat(sprintf("CSF outflow resistance      : %.2f mmHg.min/mL\n", s$rout))
    cat(sprintf("CSF formation               : %.3f mL/min (%.0f mL/day)\n",
                s$fform, s$fform * 1440))
    cat(sprintf("R_out * F_form              : %.2f mmHg\n", s$rout * s$fform))
    cat(sprintf("\nLOOP GAIN g                 : %.3f\n", s$g))
    cat(sprintf("amplification 1/(1-g)       : %.2f\n", s$amplification))
    cat(sprintf("  -> 1 mmHg delivered through R_out*F_form moves ICP %.2f mmHg\n",
                s$amplification))
    cat(sprintf("  -> 1 mmHg delivered through P_cv moves ICP 1.00 mmHg\n"))
    hm <- with(as.list(param(iih)),
               FSEG * QV0 * RTSCMAX / (4 * input$pstiff))
    cat(sprintf("\nmax collapse slope h'max    : %.3f  (%s)\n", hm,
                ifelse(hm > 1, "BISTABLE REGIME - a fold exists",
                       "monostable - one pressure at every parameter")))
  })

  output$gcurve <- renderPlot({
    d <- iih_decompose(seq(24, 50, by = 1), pars = pars())
    ggplot(d, aes(BMI, g)) +
      geom_line(colour = "#c0392b", linewidth = 1) +
      geom_hline(yintercept = 1, linetype = 2) +
      geom_vline(xintercept = input$bmi, linetype = 3) +
      labs(title = "loop gain against habitus",
           y = "g = dP_sss/dICP", x = "BMI") + THEME
  })

  output$decomp <- renderPlot({
    d <- iih_decompose(seq(24, 50, by = 1), pars = pars()) %>%
      transmute(BMI,
                `central venous P_cv` = PCV * CM,
                `venous gradient Q_v R_ts` = venous_gradient * CM,
                `CSF term R_out F_form` = Rout_Fform * CM) %>%
      pivot_longer(-BMI)
    ggplot(d, aes(BMI, value, fill = name)) +
      geom_area() +
      geom_hline(yintercept = 25, linetype = 2) +
      annotate("text", x = 25, y = 26.5, hjust = 0, size = 3.2,
               label = "diagnostic threshold 25 cmH2O") +
      geom_vline(xintercept = input$bmi, linetype = 3) +
      scale_fill_manual(values = c("#b9d2ea", "#f4c6b6", "#cdbdea")) +
      labs(title = "the three terms of the Davson relation, stacked",
           y = "cmH2O", x = "BMI") + THEME
  })

  output$decomptab <- renderDT(
    datatable(iih_decompose(seq(24, 50, by = 2), pars = pars()) %>%
                mutate(across(where(is.numeric), ~round(.x, 3))),
              options = list(pageLength = 8, dom = "tp"), rownames = FALSE))

  # ---- tab 2 ---------------------------------------------------------------
  output$icp <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, ICPCM, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 25, linetype = 2) +
      scale_colour_manual(values = c("treated" = "#1f77b4",
                                     "untreated counterfactual" = "#999999")) +
      labs(title = "opening pressure against its own counterfactual",
           x = "day from start of therapy", y = "cmH2O") + THEME
  })

  output$grad <- renderPlot({
    d <- sim() %>%
      select(time, arm, `sinus gradient (mmHg)` = GRADmm,
             `sagittal sinus pressure (mmHg)` = PSSSmm,
             `central venous pressure (mmHg)` = PCV) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c("treated" = "#1f77b4",
                                     "untreated counterfactual" = "#999999")) +
      labs(x = "day", y = NULL) + THEME
  })

  # ---- tab 3 ---------------------------------------------------------------
  output$onh <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_line(aes(y = RNFLum, colour = "OCT RNFL (measured)"),
                linetype = 2, linewidth = 1.1) +
      geom_line(aes(y = AXON, colour = "surviving axons"), linewidth = 0.9) +
      geom_line(aes(y = RNFLS, colour = "prelaminar swelling"),
                linewidth = 0.9) +
      facet_wrap(~arm) +
      scale_colour_manual(values = c("OCT RNFL (measured)" = "#2c3e50",
                                     "surviving axons" = "#2e7d32",
                                     "prelaminar swelling" = "#c0392b")) +
      labs(title = "one measured number, two states, opposite meanings",
           x = "day", y = "micrometre") + THEME
  })

  output$tlpd <- renderPlot({
    d <- sim() %>%
      select(time, arm, `ICP (mmHg)` = ICP,
             `sheath pressure PONS (mmHg)` = PONS,
             `IOP (mmHg)` = IOP,
             `translaminar gradient (mmHg)` = TLPDmm) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c("treated" = "#1f77b4",
                                     "untreated counterfactual" = "#999999")) +
      labs(x = "day", y = NULL,
           subtitle = "the disc reads PONS - IOP; a carbonic anhydrase inhibitor moves both") +
      THEME
  })

  # ---- tab 4 ---------------------------------------------------------------
  output$vision <- renderPlot({
    d <- sim() %>%
      select(time, arm, `mean deviation (dB)` = MD,
             `Frisen grade` = FRISEN,
             `cumulative exposure (mmHg.day)` = CUMEXP,
             `diplopia burden` = DIPL) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) +
      geom_line(linewidth = 0.9) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c("treated" = "#1f77b4",
                                     "untreated counterfactual" = "#999999")) +
      labs(x = "day", y = NULL) + THEME
  })

  output$delaytab <- renderDT(
    datatable(iih_delay_vs_dose(bmi = input$bmi) %>%
                mutate(across(where(is.numeric), ~round(.x, 2))),
              options = list(pageLength = 10, dom = "t"), rownames = FALSE))

  # ---- tab 5 ---------------------------------------------------------------
  output$gradsym <- renderPlot({
    d <- sim() %>%
      select(time, arm, `pulsatile tinnitus (0-10)` = TINN,
             `sinus gradient (mmHg)` = GRADmm,
             `opening pressure (cmH2O)` = ICPCM,
             `diplopia burden` = DIPL) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) +
      geom_line(linewidth = 0.9) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c("treated" = "#1f77b4",
                                     "untreated counterfactual" = "#999999")) +
      labs(x = "day", y = NULL) + THEME
  })

  output$dissoctab <- renderDT({
    m <- patient()
    b <- iih_present(m, input$runin)
    rows <- list("acetazolamide 2 g/day" = list(DOSACZ = 2000),
                 "weight programme -15 %" = list(FDIET = 0.15, TDIET = 0),
                 "venous sinus stent" = list(TSTENT = 0),
                 "CSF shunt" = list(TSHUNT = 0),
                 "ONSF" = list(TONSF = 0))
    d <- bind_rows(lapply(names(rows), function(nm) {
      e <- tail(as.data.frame(do.call(iih_run,
                c(list(m, days = 180, runin = input$runin), rows[[nm]]))), 1)
      data.frame(intervention = nm,
                 dICP_cmH2O = round(e$ICPCM - b$ICPCM, 2),
                 dGradient_mmHg = round(e$GRADmm - b$GRADmm, 2),
                 dTLPD_mmHg = round(e$TLPDmm - b$TLPDmm, 2),
                 dTinnitus = round(e$TINN - b$TINN, 2),
                 dHeadacheDays = round(e$MHD - b$MHD, 2),
                 dMD_dB = round(e$MD - b$MD, 2))
    }))
    datatable(d, options = list(dom = "t"), rownames = FALSE)
  })

  # ---- tab 6 ---------------------------------------------------------------
  output$headache <- renderPlot({
    d <- sim() %>%
      select(time, arm, `monthly headache days` = MHD,
             `central sensitisation` = STG,
             `medication overuse` = MOH, `CGRP tone` = CGRP) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) +
      geom_line(linewidth = 0.9) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c("treated" = "#1f77b4",
                                     "untreated counterfactual" = "#999999")) +
      labs(x = "day", y = NULL) + THEME
  })

  output$hatab <- renderDT(
    datatable(iih_headache(bmi = input$bmi) %>%
                mutate(across(where(is.numeric), ~round(.x, 3))),
              options = list(pageLength = 8, dom = "t"), rownames = FALSE))

  # ---- tab 7 ---------------------------------------------------------------
  output$pk <- renderPlot({
    d <- sim() %>%
      select(time, arm, `acetazolamide (mg/L)` = CACZmgL,
             `topiramate (mg/L)` = CTPMmgL, `exenatide (ng/mL)` = CEXngmL,
             `bicarbonate (mEq/L)` = HCO3, `paraesthesia (0-1)` = PARES,
             `adherence` = ADH) %>%
      pivot_longer(-c(time, arm)) %>% filter(arm == "treated")
    ggplot(d, aes(time, value)) +
      geom_line(colour = "#1f77b4", linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL,
           subtitle = "prescribed dose enters the ODE multiplied by the adherence state") +
      THEME
  })

  output$dr <- renderPlot({
    d <- iih_dose_response(bmi = input$bmi)
    p1 <- ggplot(d, aes(dose_mg_day, dICP)) +
      geom_line(colour = "#1f77b4") + geom_point() +
      labs(x = "acetazolamide (mg/day)", y = "dICP (cmH2O)",
           title = "pressure") + THEME
    p2 <- ggplot(d, aes(dose_mg_day, adherence)) +
      geom_line(colour = "#c0392b") + geom_point() +
      labs(x = "acetazolamide (mg/day)", y = "adherence",
           title = "why it flattens") + THEME
    p3 <- ggplot(d, aes(dose_mg_day, TLPD)) +
      geom_line(colour = "#2e7d32") + geom_point() +
      labs(x = "acetazolamide (mg/day)", y = "translaminar gradient (mmHg)",
           title = "at the disc") + THEME
    p1 | p2 | p3
  })

  output$drtab <- renderDT(
    datatable(iih_dose_response(bmi = input$bmi) %>%
                mutate(across(where(is.numeric), ~round(.x, 3))),
              options = list(pageLength = 10, dom = "t"), rownames = FALSE))

  # ---- tab 8 ---------------------------------------------------------------
  lp <- reactive(iih_lp(bmi = input$bmi))

  output$lpfast <- renderPlot({
    ggplot(lp()$minutes, aes(min, ICP_cmH2O)) +
      geom_line(colour = "#7c2f56", linewidth = 1) +
      labs(title = "25 mL removed at minute 0",
           x = "minutes", y = "cmH2O") + THEME
  })

  output$lpslow <- renderPlot({
    d <- lp()$days %>% pivot_longer(-day)
    ggplot(d, aes(day, value)) + geom_line(colour = "#7c2f56", linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "the same puncture over 30 days",
           x = "day", y = NULL) + THEME
  })

  # ---- tab 9 ---------------------------------------------------------------
  scen <- eventReactive(input$go, {
    base <- patient()
    bind_rows(lapply(input$scen, function(nm) {
      sc <- SCENARIOS[[nm]]
      as.data.frame(do.call(iih_run, c(list(base, days = input$days,
                                            runin = input$runin), sc$pars))) %>%
        mutate(scenario = nm, label = sc$label)
    }))
  }, ignoreNULL = FALSE)

  output$scenplot <- renderPlot({
    d <- scen() %>%
      select(time, label, `opening pressure (cmH2O)` = ICPCM,
             `OCT RNFL (um)` = RNFLum, `mean deviation (dB)` = MD,
             `monthly headache days` = MHD,
             `sinus gradient (mmHg)` = GRADmm,
             `pulsatile tinnitus` = TINN) %>%
      pivot_longer(-c(time, label))
    ggplot(d, aes(time, value, colour = label)) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL) + THEME +
      guides(colour = guide_legend(ncol = 2))
  })

  output$scentab <- renderDT({
    d <- scen() %>% group_by(label) %>% slice_tail(n = 1) %>% ungroup() %>%
      transmute(label, ICP_cmH2O = round(ICPCM, 1),
                gradient = round(GRADmm, 2), RNFL = round(RNFLum, 1),
                Frisen = round(FRISEN, 2), MD = round(MD, 2),
                headache_days = round(MHD, 1), tinnitus = round(TINN, 2),
                BMI = round(BMIout, 1), adherence = round(ADH, 2))
    datatable(d, options = list(dom = "t"), rownames = FALSE)
  })

  # ---- tab 10 --------------------------------------------------------------
  output$spectrum <- renderPlot({
    d <- iih_responder_spectrum()
    p1 <- ggplot(d, aes(g, -dICP_cmH2O)) + geom_line(colour = "#1f77b4") +
      geom_point() +
      labs(x = "loop gain g", y = "ICP fall (cmH2O)",
           title = "identical 30 % secretion cut, different sinus mechanics") +
      THEME
    p2 <- ggplot(d, aes(PCRIT, ICP_cmH2O)) + geom_line(colour = "#c0392b") +
      geom_hline(yintercept = 25, linetype = 2) +
      labs(x = "PCRIT (mmHg)", y = "untreated ICP (cmH2O)",
           title = "and different untreated pressures") + THEME
    p1 | p2
  })

  output$spectab <- renderDT(
    datatable(iih_responder_spectrum() %>%
                mutate(across(where(is.numeric), ~round(.x, 3))),
              options = list(pageLength = 8, dom = "tp"), rownames = FALSE))

  output$contplot <- renderPlot({
    d <- iih_continuation(pcrit = input$pcrit)
    ggplot(d, aes(K_mmHg)) +
      geom_line(aes(y = lowest_cmH2O, colour = "lowest stable pressure")) +
      geom_line(aes(y = highest_cmH2O, colour = "highest stable pressure")) +
      geom_point(aes(y = highest_cmH2O, size = n_stable > 1),
                 colour = "#c0392b") +
      scale_size_manual(values = c(`FALSE` = 0, `TRUE` = 2),
                        guide = "none") +
      labs(x = "K = R_out * F_form (mmHg)", y = "cmH2O",
           title = "continuation in the CSF-side parameter; red = bistable") +
      THEME
  })

  # ---- tab 11 --------------------------------------------------------------
  output$caltab <- renderDT({
    m <- iih_patient(bmi = 38)
    b <- iih_present(m, 120)
    arms <- list(
      list("Sinclair 2010 diet -15.7 %, 3 mo", list(FDIET = 0.157, TDIET = 0),
           90, -8.4),
      list("Mitchell 2023 exenatide, 12 wk", list(DOSEX = 20), 84, -5.6),
      list("Markey 2020 AZD4017, 12 wk", list(DOSAZD = 800), 84, 0.3),
      list("IIHTT placebo arm, 6 mo", list(FDIET = 0.033, TDIET = 0), 180, NA),
      list("IIHTT acetazolamide arm, 6 mo",
           list(DOSACZ = 2000, FDIET = 0.070, TDIET = 0), 180, NA),
      list("IIH:WT bariatric, 12 mo", list(FBAR = 0.21, TBAR = 0), 365, NA))
    d <- bind_rows(lapply(arms, function(a) {
      e <- tail(as.data.frame(do.call(iih_run,
                c(list(m, days = a[[3]]), a[[2]]))), 1)
      data.frame(arm = a[[1]], day = a[[3]],
                 simulated_dICP_cmH2O = round(e$ICPCM - b$ICPCM, 2),
                 published_dICP_cmH2O = a[[4]],
                 simulated_dMD_dB = round(e$MD - b$MD, 2),
                 simulated_dFrisen = round(e$FRISEN - b$FRISEN, 2))
    }))
    datatable(d, options = list(dom = "t"), rownames = FALSE)
  })

  output$biases <- renderPrint({
    cat(
"1. The between-arm ICP difference in the IIHTT reconstruction is about 1.4x\n",
"   larger than published (-6.1 vs -4.4 cmH2O).  The published placebo arm\n",
"   fell nearly 10 cmH2O on a 3 % weight loss, which no mechanism here\n",
"   reproduces, and the 6-month opening pressure was a single measurement in\n",
"   a subset.  Not tuned away.\n\n",
"2. The simulated presenting transverse sinus gradient (8-10 mmHg) sits at\n",
"   the low end of the 10-25 mmHg reported in stented cohorts, which are\n",
"   selected for severe stenosis.\n\n",
"3. Body weight is UNAMPLIFIED in this model only because the collapse is\n",
"   taken to be driven by the mid-segment sinus pressure (FSEG = 0.5).  At\n",
"   FSEG = 0 (collapse driven at the jugular end) weight becomes amplified\n",
"   as well.  This is the largest structural uncertainty in the model and\n",
"   the whole FSEG sensitivity is printed in ANALYSIS 1 of the Python\n",
"   reference rather than hidden.\n\n",
"4. Dosing is continuous rather than event-based so that the tolerability ->\n",
"   adherence -> exposure loop can close inside the ODE system.  For\n",
"   acetazolamide the daily-average concentration that drives the\n",
"   pharmacology differs from 500 mg TID by < 4 % because the secretory\n",
"   state has a 3.6 h time constant.\n\n",
"5. At the trial-calibrated parameters the loop gain never exceeds ~0.6, so\n",
"   the model is MONOSTABLE: there is one pressure at every parameter set,\n",
"   and a finite course of a drug cannot be curative.  The fold requires a\n",
"   sinus roughly three times more compliant (PSTIFF < 1.5).  That is a\n",
"   prediction about which patients could ever stop treatment, and it was\n",
"   not what the model was built expecting.\n")
  })
}

shinyApp(ui, server)
