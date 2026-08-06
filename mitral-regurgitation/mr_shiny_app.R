## =============================================================================
##  Mitral Regurgitation -- interactive QSP dashboard
##  ---------------------------------------------------------------------------
##  The app is organised around the model's organising idea rather than around
##  its state variables: a regurgitant volume means nothing until it has been
##  divided by something, and each tab is one of the denominators.
##
##    Tab  1  Patient          build a patient, see the beat solved
##    Tab  2  Denominator 1    atrial compliance: acute versus chronic
##    Tab  3  Denominator 2    regurgitant fraction and the low-flow trap
##    Tab  4  Denominator 3    proportionality, COAPT versus MITRA-FR
##    Tab  5  Denominator 4    ejection fraction and the operative threshold
##    Tab  6  Denominator 5    PISA geometry and the severity threshold
##    Tab  7  Pharmacokinetics drug exposure and where each agent enters
##    Tab  8  The vortex       coaptation reserve, leaflets, loop gain
##    Tab  9  Right heart      pulmonary vascular ratchet, the second barrier
##    Tab 10  Endpoints        symptoms, hospitalisation, survival
##    Tab 11  Scenarios        side-by-side comparison of interventions
##    Tab 12  Biomarkers       BNP, eGFR, neurohormonal state
##
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2
##  Run with: shiny::runApp("mr_shiny_app.R")
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

mod <- mread_cache("mr_qsp", "mr_mrgsolve_model.R")

## ---- the two trial-matched virtual patients, from mr_reference_output.txt ---
PATIENTS <- list(
  "Healthy" = list(),
  "Acute severe MR (papillary rupture)" = list(EROApri = 0.60),
  "Mild chronic MR" = list(EROApri = 0.10),
  "Compensated severe primary MR" = list(EROApri = 0.42, V0d = 47.0,
                                         LVmass = 205.0, V0la = 52.0),
  "COAPT-like (disproportionate)" = list(
    Ees = 0.9680, V0d = 80.8715, LVmass = 266.53, Fib = 0.0448,
    V0la = 81.274, Fibla = 0.0886, AFb = 0.0260, Ann = 7.6352,
    Aleaf = 10.2408, Rpul = 0.0930, Rpulfix = 0.0133, Eesrv = 0.4620,
    Vtot = 5299.49),
  "MITRA-FR-like (proportionate)" = list(
    Ees = 1.9503, V0d = 141.1263, LVmass = 319.87, Fib = 0.0376,
    V0la = 33.026, Fibla = 0.0499, Ann = 8.2989, Aleaf = 12.5540,
    Vtot = 4907.68),
  "Atrial functional MR (AF, normal LV)" = list(AFb = 0.90, V0la = 190,
                                                Ann = 9.5)
)

GDMT <- list(rate_fur = 40, rate_bb = 100, rate_sac = 194,
             rate_val = 206, rate_mra = 25, rate_sg = 10)

theme_mr <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"),
          legend.position = "bottom")
}

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("Mitral regurgitation — one regurgitant volume, five denominators"),
  tags$p(style = "color:#555;max-width:1000px;",
         paste("A regurgitant orifice produces one number that echocardiography",
               "reports. That number is meaningless until it is divided by",
               "something: atrial compliance decides congestion, total stroke",
               "volume decides the grade, ventricular volume decides whether the",
               "valve or the ventricle is the disease, the afterload the leak",
               "removes hides contractility, and PISA geometry biases the",
               "measurement itself.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("patient", "Patient", names(PATIENTS),
                  selected = "COAPT-like (disproportionate)"),
      sliderInput("tend", "Follow-up (days)", 30, 3650, 730, step = 30),
      checkboxInput("gdmt", "Guideline-directed medical therapy", TRUE),
      hr(),
      h5("Lesion"),
      sliderInput("eropri", "Degenerative orifice EROA (cm2)", 0, 1.0, 0, 0.01),
      sliderInput("kprog", "Degenerative progression (1/day)", 0, 0.002, 0, 0.0001),
      hr(),
      h5("Device / surgery"),
      sliderInput("teerf", "Orifice abolished by device (fraction)", 0, 1, 0, 0.01),
      sliderInput("teerrmv", "Mitral resistance added (mmHg.s/mL)",
                  0, 0.08, 0, 0.005),
      checkboxInput("ring", "Annuloplasty ring (clamp the annulus)", FALSE),
      sliderInput("ringann", "Ring annular area (cm2)", 4, 10, 5.5, 0.1),
      checkboxInput("crt", "Cardiac resynchronisation", FALSE),
      hr(),
      h5("Drugs"),
      sliderInput("fur", "Furosemide (mg/day)", 0, 240, 40, 10),
      sliderInput("arni", "Sacubitril/valsartan (mg/day, total)", 0, 400, 0, 20),
      sliderInput("bb", "Beta-blocker (mg/day)", 0, 200, 0, 10),
      sliderInput("mra", "MRA (mg/day)", 0, 50, 0, 5),
      sliderInput("sg", "SGLT2 inhibitor (mg/day)", 0, 10, 0, 1),
      sliderInput("snp", "Nitroprusside (effect units)", 0, 8, 0, 0.5),
      sliderInput("dob", "Dobutamine (effect units)", 0, 8, 0, 0.5),
      hr(),
      h5("Physiology overrides"),
      sliderInput("hrfix", "Clamp heart rate (0 = off)", 0, 130, 0, 5),
      sliderInput("kpisa", "PISA inflation for a crescent", 0, 1.0, 0.45, 0.05)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 Patient",
                 h4("The beat, solved"),
                 tableOutput("beat"),
                 plotOutput("p_overview", height = "460px")),
        tabPanel("2 Denominator 1 · atrium",
                 h4("Congestion is the regurgitant volume divided by the operating atrial compliance"),
                 tags$p("A virgin atrium sits high on a steep pressure-volume curve, so a regurgitant volume it has never seen produces a giant v wave. A chronically dilated atrium has shifted that curve right and flattened it, so the identical volume is nearly silent. Atrial fibrosis then stiffens it back, which is how a long-compensated patient decompensates with an unchanged valve."),
                 plotOutput("p_den1", height = "420px"),
                 h4("The same regurgitant volume, swept across atrial compliance"),
                 plotOutput("p_den1b", height = "300px")),
        tabPanel("3 Denominator 2 · fraction",
                 h4("Regurgitant fraction rises when the ventricle fails, with the orifice unchanged"),
                 tags$p("Regurgitant fraction is the leak divided by total stroke volume. In a low-output patient the same regurgitant volume is a larger fraction, so the reported severity climbs as the ventricle deteriorates — the grade is partly a statement about the ventricle."),
                 plotOutput("p_den2", height = "380px"),
                 plotOutput("p_den2b", height = "320px")),
        tabPanel("4 Denominator 3 · proportionality",
                 h4("Valve or ventricle? The COAPT / MITRA-FR axis"),
                 tags$p("The same device, abolishing the same fraction of the orifice, is applied to both trial-matched patients. Any difference comes only from how large the orifice is relative to the ventricle behind it."),
                 plotOutput("p_den3", height = "400px"),
                 h4("Device efficacy sweep"),
                 tags$p("In the disproportionate patient the endpoint improves monotonically. In the proportionate patient a perfect valve barely moves it."),
                 plotOutput("p_den3b", height = "320px")),
        tabPanel("5 Denominator 4 · ejection fraction",
                 h4("Why the operative threshold is 60% and not 50%"),
                 tags$p("Ejection fraction is computed against a stroke volume that contains the leak, and the leak is ejected into a low-pressure sink. Abolish the orifice with contractility untouched and EF falls while forward output rises: the two facts are the same fact seen through different denominators."),
                 tableOutput("t_den4"),
                 plotOutput("p_den4", height = "400px")),
        tabPanel("6 Denominator 5 · measurement",
                 h4("PISA assumes a hemisphere; the functional orifice is a crescent"),
                 tags$p("The severity threshold of 0.4 cm2 was calibrated on round degenerative orifices and is then applied to crescentic functional ones. Slide the inflation factor and watch the reported grade move without the physiology changing at all."),
                 plotOutput("p_den5", height = "360px"),
                 h4("Self-consistency audit of the reported trial baselines"),
                 tags$p("Take the reported LVEDV, LVEF and PISA-EROA at face value and ask what cardiac index the patient must have been living at."),
                 tableOutput("t_den5")),
        tabPanel("7 Pharmacokinetics",
                 h4("Exposure, and where each agent enters the physiology"),
                 plotOutput("p_pk", height = "420px"),
                 tableOutput("t_pk")),
        tabPanel("8 The vortex",
                 h4("MR begets MR: coaptation reserve, leaflet supply, loop gain"),
                 tags$p("The orifice in secondary MR is not a hole but a geometric deficit: regurgitation exists only when the coaptation area the geometry demands exceeds the leaflet area available. Leaflet growth is rate-limited and only partial, which is why the speed of dilation matters and not only its final size."),
                 plotOutput("p_vortex", height = "440px"),
                 plotOutput("p_vortex2", height = "300px")),
        tabPanel("9 Right heart",
                 h4("The second barrier: an irreversible pulmonary vascular ratchet"),
                 tags$p("A fraction of every rise in pulmonary vascular resistance never regresses. An identical procedure therefore buys progressively less the longer it is deferred: timing is part of the intervention, not a scheduling detail."),
                 plotOutput("p_rv", height = "440px")),
        tabPanel("10 Endpoints",
                 h4("Symptoms, hospitalisation, survival"),
                 tags$p("Two numbers in the whole model are fitted to outcome data, both from the COAPT control arm. Everything else here is prediction."),
                 plotOutput("p_end", height = "440px"),
                 tableOutput("t_end")),
        tabPanel("11 Scenarios",
                 h4("Which arm of the feedback loop does each procedure cut?"),
                 tags$p("An annuloplasty ring clamps the annulus but leaves the tethering arm intact, and the ventricle keeps dilating underneath it. Replacement holds, at the cost of a gradient."),
                 plotOutput("p_scen", height = "520px"),
                 tableOutput("t_scen")),
        tabPanel("12 Biomarkers",
                 h4("BNP, renal function and the neurohormonal state"),
                 plotOutput("p_bio", height = "440px"))
      )
    )
  )
)

## =============================================================================
##  server
## =============================================================================
server <- function(input, output, session) {

  base_init <- reactive({
    iv <- PATIENTS[[input$patient]]
    if (input$eropri > 0) iv$EROApri <- input$eropri
    if (input$teerf > 0)  iv$TEERf <- input$teerf
    if (input$teerrmv > 0) iv$TEERrmv <- input$teerrmv
    if (isTRUE(input$ring)) { iv$RingF <- 1; iv$Ann <- input$ringann }
    iv
  })

  base_param <- reactive({
    p <- list(k_pri_prog = input$kprog, crt = as.numeric(input$crt),
              k_pisa_sec = input$kpisa, snp_inf = input$snp, dob_inf = input$dob,
              rate_fur = input$fur, rate_bb = input$bb, rate_mra = input$mra,
              rate_sg = input$sg,
              rate_sac = input$arni * 0.485, rate_val = input$arni * 0.515)
    if (isTRUE(input$gdmt)) p <- modifyList(GDMT, p[lengths(p) > 0])
    if (input$hrfix > 0) p$HR_fix <- input$hrfix
    p
  })

  sim <- function(init = list(), param = list(), tend = NULL, delta = NULL) {
    tend <- if (is.null(tend)) input$tend else tend
    delta <- if (is.null(delta)) max(tend / 400, 0.05) else delta
    mod %>%
      init(modifyList(base_init(), init)) %>%
      param(modifyList(base_param(), param)) %>%
      mrgsim(end = tend, delta = delta) %>%
      as_tibble()
  }

  run_main <- reactive(sim())

  ## ---- tab 1 -------------------------------------------------------------
  output$beat <- renderTable({
    d <- run_main(); a <- d[1, ]; b <- d[nrow(d), ]
    tibble(
      quantity = c("EROA, true (cm2)", "EROA, PISA-reported (cm2)",
                   "regurgitant volume (mL)", "regurgitant fraction",
                   "LV end-diastolic volume (mL)", "LVEF", "forward EF",
                   "cardiac index (L/min/m2)", "mean LA pressure (mmHg)",
                   "v wave (mmHg)", "operating LA compliance (mL/mmHg)",
                   "effective wedge (mmHg)", "mean PA pressure (mmHg)",
                   "CVP (mmHg)", "mitral gradient (mmHg)",
                   "coaptation reserve (cm2)"),
      day_0 = c(a$EROA, a$EROApisa, a$RVol, a$RF, a$EDV, a$EF, a$EFfwd, a$CI,
                a$LAP, a$vwave, a$ClaOp, a$Ppcw, a$PAP, a$CVP, a$dPmv,
                a$CoaptResv),
      final = c(b$EROA, b$EROApisa, b$RVol, b$RF, b$EDV, b$EF, b$EFfwd, b$CI,
                b$LAP, b$vwave, b$ClaOp, b$Ppcw, b$PAP, b$CVP, b$dPmv,
                b$CoaptResv))
  }, digits = 3)

  output$p_overview <- renderPlot({
    run_main() %>%
      select(time, EROA, RVol, RF, EDV, EF, CI, Ppcw, LAvol) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(name, c("EROA", "RVol", "RF", "EDV", "EF", "CI",
                                   "Ppcw", "LAvol"))) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.8, colour = "#22558A") +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "days", y = NULL) + theme_mr()
  })

  ## ---- tab 2: denominator 1 ---------------------------------------------
  output$p_den1 <- renderPlot({
    run_main() %>%
      select(time, RVol, ClaOp, vwave, LAP, Ppcw, LAvol) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.8, colour = "#1F6B60") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "days", y = NULL,
           title = "The numerator, the denominator, and what comes out") +
      theme_mr()
  })

  output$p_den1b <- renderPlot({
    ## sweep the atrium from virgin to chronically dilated, orifice fixed
    grid <- lapply(seq(32, 210, by = 8), function(v) {
      d <- sim(init = list(V0la = v, EROApri = max(input$eropri, 0.45)),
               tend = 0.5, delta = 0.5)
      tibble(V0la = v, ClaOp = d$ClaOp[1], vwave = d$vwave[1],
             Ppcw = d$Ppcw[1], RVol = d$RVol[1])
    }) %>% bind_rows()
    grid %>% pivot_longer(-V0la) %>%
      ggplot(aes(V0la, value)) +
      geom_line(linewidth = 0.9, colour = "#1F6B60") +
      facet_wrap(~name, scales = "free_y", nrow = 1) +
      labs(x = "LA unstressed volume (mL) — virgin to chronically dilated",
           y = NULL,
           title = "Hold the leak fixed and move only the atrium") +
      theme_mr()
  })

  ## ---- tab 3: denominator 2 ---------------------------------------------
  output$p_den2 <- renderPlot({
    run_main() %>%
      select(time, RVol, SVtot, SVfwd, RF) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.8, colour = "#8E2B3F") +
      facet_wrap(~name, scales = "free_y", nrow = 1) +
      labs(x = "days", y = NULL,
           title = "The leak, the total, and the ratio of the two") + theme_mr()
  })

  output$p_den2b <- renderPlot({
    ## the low-flow trap: drop contractility, hold the orifice fixed
    grid <- lapply(seq(0.6, 3.0, by = 0.15), function(e) {
      d <- sim(init = list(Ees = e, EROApri = max(input$eropri, 0.35)),
               tend = 0.5, delta = 0.5)
      tibble(Ees = e, RVol = d$RVol[1], SVtot = d$SVtot[1], RF = d$RF[1],
             CI = d$CI[1])
    }) %>% bind_rows()
    grid %>% pivot_longer(-Ees) %>%
      ggplot(aes(Ees, value)) +
      geom_line(linewidth = 0.9, colour = "#8E2B3F") +
      facet_wrap(~name, scales = "free_y", nrow = 1) +
      labs(x = "contractility E_es (mmHg/mL)", y = NULL,
           title = "The low-flow trap: severity grade rises as the ventricle fails, orifice unchanged") +
      theme_mr()
  })

  ## ---- tab 4: denominator 3 ---------------------------------------------
  output$p_den3 <- renderPlot({
    both <- lapply(c("COAPT-like (disproportionate)",
                     "MITRA-FR-like (proportionate)"), function(nm) {
      lapply(c(0, input$teerf %||% 0.68), function(tf) {
        d <- mod %>% init(modifyList(PATIENTS[[nm]],
                                     list(TEERf = tf, TEERrmv = 0.030))) %>%
          param(modifyList(GDMT, base_param())) %>%
          mrgsim(end = input$tend, delta = input$tend / 300) %>% as_tibble()
        d$patient <- nm
        d$arm <- ifelse(tf > 0, "device", "control")
        d
      }) %>% bind_rows()
    }) %>% bind_rows()
    both %>% select(time, patient, arm, RVol, RF, CI, Ppcw, HFH, SurvP) %>%
      pivot_longer(c(-time, -patient, -arm)) %>%
      ggplot(aes(time, value, colour = arm, linetype = patient)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = c(control = "#8E2B3F", device = "#22558A")) +
      labs(x = "days", y = NULL, colour = NULL, linetype = NULL,
           title = "Same device, two denominators") + theme_mr()
  })

  output$p_den3b <- renderPlot({
    grid <- lapply(c("COAPT-like (disproportionate)",
                     "MITRA-FR-like (proportionate)"), function(nm) {
      lapply(seq(0, 1, by = 0.1), function(tf) {
        d <- mod %>% init(modifyList(PATIENTS[[nm]],
                                     list(TEERf = tf, TEERrmv = 0.030))) %>%
          param(modifyList(GDMT, base_param())) %>%
          mrgsim(end = 365, delta = 5) %>% as_tibble()
        tibble(patient = nm, TEERf = tf, HFH = tail(d$HFH, 1),
               death = 1 - tail(d$SurvP, 1), Ppcw = tail(d$Ppcw, 1),
               CI = tail(d$CI, 1))
      }) %>% bind_rows()
    }) %>% bind_rows()
    grid %>% pivot_longer(c(-patient, -TEERf)) %>%
      ggplot(aes(TEERf, value, colour = patient)) +
      geom_line(linewidth = 1) + geom_point(size = 1.4) +
      facet_wrap(~name, scales = "free_y", nrow = 1) +
      scale_colour_manual(values = c("#22558A", "#B08A20")) +
      labs(x = "fraction of the orifice abolished", y = NULL, colour = NULL,
           title = "A perfect valve does not rescue a proportionate patient") +
      theme_mr()
  })

  ## ---- tab 5: denominator 4 ---------------------------------------------
  den4 <- reactive({
    pre <- sim(tend = 0.5, delta = 0.5)
    post <- sim(init = list(EROApri = 0, TEERf = 1), tend = 0.5, delta = 0.5)
    list(pre = pre[1, ], post = post[1, ])
  })

  output$t_den4 <- renderTable({
    d <- den4()
    tibble(quantity = c("EROA (cm2)", "regurgitant volume (mL)",
                        "total stroke volume (mL)", "forward stroke volume (mL)",
                        "LVEF", "forward EF", "end-systolic pressure (mmHg)",
                        "end-systolic wall stress (mmHg)",
                        "cardiac index (L/min/m2)", "effective wedge (mmHg)"),
           pre_op = c(d$pre$EROA, d$pre$RVol, d$pre$SVtot, d$pre$SVfwd,
                      d$pre$EF, d$pre$EFfwd, d$pre$Pes, d$pre$sigES,
                      d$pre$CI, d$pre$Ppcw),
           post_op = c(d$post$EROA, d$post$RVol, d$post$SVtot, d$post$SVfwd,
                       d$post$EF, d$post$EFfwd, d$post$Pes, d$post$sigES,
                       d$post$CI, d$post$Ppcw))
  }, digits = 3)

  output$p_den4 <- renderPlot({
    grid <- lapply(seq(0.7, 5.0, by = 0.2), function(e) {
      a <- sim(init = list(Ees = e), tend = 0.5, delta = 0.5)
      b <- sim(init = list(Ees = e, EROApri = 0, TEERf = 1),
               tend = 0.5, delta = 0.5)
      tibble(Ees = e, pre_EF = a$EF[1], post_EF = b$EF[1],
             pre_CI = a$CI[1], post_CI = b$CI[1])
    }) %>% bind_rows()
    ggplot(grid, aes(pre_EF)) +
      geom_line(aes(y = post_EF), linewidth = 1, colour = "#8E2B3F") +
      geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey50") +
      geom_hline(yintercept = 0.50, linetype = 2, colour = "#22558A") +
      geom_vline(xintercept = 0.60, linetype = 2, colour = "#B08A20") +
      annotate("text", x = 0.605, y = 0.20, hjust = 0, size = 3.4,
               colour = "#B08A20", label = "guideline: operate at EF 0.60") +
      annotate("text", x = min(grid$pre_EF), y = 0.515, hjust = 0, size = 3.4,
               colour = "#22558A", label = "post-operative EF 0.50") +
      labs(x = "PRE-operative ejection fraction",
           y = "POST-operative ejection fraction",
           title = "The threshold is derived, not conventional",
           subtitle = paste("Where the red curve crosses the blue line is the",
                            "pre-operative EF at which the operation leaves a",
                            "failing ventricle")) +
      theme_mr()
  })

  ## ---- tab 6: denominator 5 ---------------------------------------------
  output$p_den5 <- renderPlot({
    grid <- lapply(seq(0, 1, by = 0.05), function(k) {
      d <- sim(param = list(k_pisa_sec = k), tend = 0.5, delta = 0.5)
      tibble(k_pisa_sec = k, EROA_true = d$EROA[1],
             EROA_reported = d$EROApisa[1])
    }) %>% bind_rows()
    grid %>% pivot_longer(-k_pisa_sec) %>%
      ggplot(aes(k_pisa_sec, value, colour = name)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0.40, linetype = 2, colour = "#B08A20") +
      annotate("text", x = 0, y = 0.415, hjust = 0, size = 3.4,
               colour = "#B08A20",
               label = "0.4 cm2: the 'severe' threshold, calibrated on round orifices") +
      scale_colour_manual(values = c(EROA_true = "#22558A",
                                     EROA_reported = "#8E2B3F")) +
      labs(x = "PISA inflation factor for a crescentic orifice",
           y = "orifice area (cm2)", colour = NULL,
           title = "The same valve, graded differently") + theme_mr()
  })

  output$t_den5 <- renderTable({
    ## reported baselines taken at face value
    trials <- tibble(
      trial = c("COAPT", "MITRA-FR"),
      EDV = c(192.7, 256.5), EF = c(0.313, 0.333), EROA_pisa = c(0.41, 0.31),
      HR = c(76.2, 70.3), Pes = c(93, 93), LAP = c(20, 15))
    out <- lapply(seq_len(nrow(trials)), function(i) {
      tr <- trials[i, ]
      lapply(c(1.0, 1.2, 1.45, 1.8), function(k) {
        ero <- tr$EROA_pisa / k
        Tr_ <- (0.42 - 0.0016 * tr$HR) + 0.070
        rv <- ero * Tr_ * 40 * sqrt(max(0.90 * tr$Pes - tr$LAP, 0.5))
        svf <- tr$EDV * tr$EF - rv
        tibble(trial = tr$trial, k_PISA = k, EROA_true = ero,
               RVol = rv, SV_forward = svf,
               CI = svf * tr$HR / 1000 / 1.90)
      }) %>% bind_rows()
    }) %>% bind_rows()
    out
  }, digits = 3)

  ## ---- tab 7: PK ---------------------------------------------------------
  output$p_pk <- renderPlot({
    sim(tend = min(input$tend, 60)) %>%
      select(time, Ac_fur, Ac_sac, Ac_val, Ac_bb, Ac_mra, Ac_sg,
             Ce_snp, Ce_dob) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.8, colour = "#1F6B62") +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "days", y = "amount (mg) or effect units", NULL,
           title = "Exposure") + theme_mr()
  })

  output$t_pk <- renderTable({
    tibble(agent = c("Furosemide", "Sacubitril/valsartan", "Beta-blocker",
                     "MRA", "SGLT2 inhibitor", "Nitroprusside", "Dobutamine"),
           enters_the_model_at = c(
             "natriuresis, with braking and a hypovolaemia floor; activates RAAS",
             "neprilysin: venous capacitance, natriuresis, antifibrotic; AT1: SVR, fibrosis",
             "heart rate down; contractility down acutely, RECOVERED chronically",
             "antifibrotic in both ventricle and atrium, mild natriuresis",
             "lowers the blood-volume SETPOINT by ~6%, mildly antifibrotic",
             "systemic resistance down 42% — redistributes flow forward at once",
             "contractility and rate up"))
  })

  ## ---- tab 8: the vortex -------------------------------------------------
  output$p_vortex <- renderPlot({
    run_main() %>%
      select(time, Areq, Aleaf, CoaptResv, EROA, CDepth, Ann, EDV, RF) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.8, colour = "#8E3F18") +
      geom_hline(data = ~filter(.x, name == "CoaptResv"),
                 aes(yintercept = 0), linetype = 2, colour = "grey40") +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "days", y = NULL,
           title = "Demand, supply, and the reserve between them — MR begins when the reserve crosses zero") +
      theme_mr()
  })

  output$p_vortex2 <- renderPlot({
    grid <- lapply(seq(30, 220, by = 10), function(v) {
      d <- sim(init = list(V0d = v), tend = 0.5, delta = 0.5)
      tibble(V0d = v, EDV = d$EDV[1], EROA = d$EROA[1],
             CoaptResv = d$CoaptResv[1], RF = d$RF[1])
    }) %>% bind_rows()
    grid %>% pivot_longer(c(-V0d, -EDV)) %>%
      ggplot(aes(EDV, value)) +
      geom_line(linewidth = 0.9, colour = "#8E3F18") +
      facet_wrap(~name, scales = "free_y", nrow = 1) +
      labs(x = "LV end-diastolic volume (mL)", y = NULL,
           title = "Sweep chamber size and watch the orifice open") + theme_mr()
  })

  ## ---- tab 9: right heart -----------------------------------------------
  output$p_rv <- renderPlot({
    waits <- c(0, 365, 730, 1460, 2190)
    out <- lapply(waits, function(w) {
      d <- mod %>% init(base_init()) %>% param(base_param()) %>%
        mrgsim(end = w + 365, delta = 15) %>% as_tibble()
      d$wait <- paste0(round(w / 365, 1), " y")
      d
    }) %>% bind_rows()
    out %>% select(time, wait, PVRWU, Rpulfix, PAP, CVP, Ppcw, CI) %>%
      pivot_longer(c(-time, -wait)) %>%
      ggplot(aes(time, value, colour = wait)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "days", y = NULL, colour = "deferral",
           title = "The irreversible component only goes up") + theme_mr()
  })

  ## ---- tab 10: endpoints -------------------------------------------------
  output$p_end <- renderPlot({
    run_main() %>%
      select(time, NYHAi, HFH, SurvP, Ppcw, CI, Oedema) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.8, colour = "#41415C") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "days", y = NULL) + theme_mr()
  })

  output$t_end <- renderTable({
    d <- run_main(); b <- d[nrow(d), ]
    tibble(endpoint = c("symptom index (1-4.4)",
                        "expected HF hospitalisations",
                        "annualised HF hospitalisation rate",
                        "survival probability", "days with alveolar oedema",
                        "eGFR (mL/min/1.73m2)", "BNP (x normal)"),
           value = c(b$NYHAi, b$HFH, b$HFH / (max(b$time, 1) / 365), b$SurvP,
                     sum(d$Oedema) * mean(diff(d$time)), b$eGFR, b$BNP))
  }, digits = 3)

  ## ---- tab 11: scenarios -------------------------------------------------
  scen <- reactive({
    defs <- list(
      "medical therapy only" = list(),
      "TEER (edge-to-edge)" = list(TEERf = 0.68, TEERrmv = 0.030),
      "annuloplasty ring" = list(RingF = 1, Ann = 5.50),
      "ring + leaflet repair" = list(RingF = 1, Ann = 5.50, TEERf = 0.45),
      "replacement, chordal-sparing" = list(TEERf = 0.97, TEERrmv = 0.012),
      "replacement, chordae divided" = list(TEERf = 0.97, TEERrmv = 0.012)
    )
    lapply(names(defs), function(nm) {
      iv <- modifyList(base_init(), defs[[nm]])
      if (nm == "replacement, chordae divided") {
        e0 <- if (!is.null(iv$Ees)) iv$Ees else 2.90
        iv$Ees <- e0 * 0.90
      }
      d <- mod %>% init(iv) %>% param(base_param()) %>%
        mrgsim(end = input$tend, delta = input$tend / 300) %>% as_tibble()
      d$scenario <- nm
      d
    }) %>% bind_rows()
  })

  output$p_scen <- renderPlot({
    scen() %>% select(time, scenario, EROA, RF, EDV, Ppcw, CI, dPmv,
                      HFH, SurvP, NYHAi) %>%
      pivot_longer(c(-time, -scenario)) %>%
      ggplot(aes(time, value, colour = scenario)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "days", y = NULL, colour = NULL) + theme_mr()
  })

  output$t_scen <- renderTable({
    scen() %>% group_by(scenario) %>% slice_tail(n = 1) %>% ungroup() %>%
      transmute(scenario, EROA, RF, EDV, Ppcw, CI, mitral_gradient = dPmv,
                HFH, survival = SurvP)
  }, digits = 3)

  ## ---- tab 12: biomarkers ------------------------------------------------
  output$p_bio <- renderPlot({
    run_main() %>%
      select(time, BNP, eGFR, NE, Ang, Ald, Vtot, Fib, Fibla, AFb) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) +
      geom_line(linewidth = 0.8, colour = "#84245C") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "days", y = NULL) + theme_mr()
  })
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

shinyApp(ui, server)
