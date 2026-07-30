## =============================================================================
##  lipo_shiny_app.R
##  Lipodystrophy syndromes QSP model — interactive dashboard
##
##  The app is organised around the model's one structural claim: adipose tissue
##  does two jobs, so the disease has two deficits, and each drug touches only
##  one of them.  Tab 3 (the buffer) and tab 4 (the leptin transducer) are the
##  two halves; every other tab is a consequence of one of them.
##
##  Run:
##    setwd("lipodystrophy-syndromes")
##    shiny::runApp("lipo_shiny_app.R")
##
##  Requires: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
##  (falls back to base graphics if ggplot2 is unavailable)
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

source("lipo_mrgsolve_model.R")

PHENO_CHOICES <- setNames(names(phenotypes),
                          vapply(phenotypes, function(p) p$label, character(1)))

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10),
        legend.position = "bottom")

## Two-arm colour convention used everywhere: grey = untreated, colour = treated
ARMCOL <- c("untreated" = "#8c8c8c", "treated" = "#1f6feb",
            "comparator" = "#e07b39")

## -----------------------------------------------------------------------------
##  UI
## -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Lipodystrophy syndromes — QSP simulator"),
  tags$p(style = "color:#555;margin-top:-8px",
         paste("One lesion, two deficits: adipose tissue is both a finite",
               "triglyceride buffer and the endocrine reporter of its own size.",
               "Educational / research model — not for clinical use.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("pheno", "Phenotype (a capacity distribution, not a label)",
                  choices = PHENO_CHOICES, selected = "cgl"),
      sliderInput("fc1", "Femorogluteal capacity (x normal)", 0, 2, 0.02, 0.01),
      sliderInput("fc2", "Upper-trunk/facial capacity (x normal)", 0, 2, 0.05, 0.01),
      sliderInput("fcv", "Visceral capacity (x normal)", 0, 2, 0.05, 0.01),
      sliderInput("pread", "Preadipocyte reserve (gates PPARg)", 0, 1, 0.03, 0.01),
      helpText("Changing the phenotype selector resets these three ceilings."),

      hr(), h4("Metreleptin"),
      checkboxInput("met_on", "Metreleptin", TRUE),
      sliderInput("met_dose", "Dose (mg/kg/d)", 0, 0.20, 0.06, 0.01),
      sliderInput("met_start", "Start (year)", 0, 10, 2, 0.5),
      checkboxInput("met_stop", "Withdraw part-way", FALSE),
      sliderInput("met_stopy", "Withdraw at (year)", 1, 15, 6, 0.5),
      checkboxInput("neut", "Neutralising anti-drug antibodies", FALSE),
      checkboxInput("clamp", "Clamp food intake (pair-fed)", FALSE),

      hr(), h4("Other therapy"),
      checkboxInput("pio_on", "Pioglitazone 45 mg", FALSE),
      checkboxInput("vol_on", "Volanesorsen 300 mg q1w (APOC3 ASO)", FALSE),
      checkboxInput("evi_on", "Evinacumab 900 mg q4w (anti-ANGPTL3)", FALSE),
      checkboxInput("fib_on", "Fenofibrate 145 mg", FALSE),
      checkboxInput("om3_on", "Omega-3 4 g/d", FALSE),
      checkboxInput("glp_on", "GLP-1 RA 1 mg q1w", FALSE),
      checkboxInput("met2_on", "Metformin 2 g/d", FALSE),
      checkboxInput("sgl_on", "SGLT2 inhibitor 10 mg", FALSE),
      sliderInput("ins_dose", "Insulin (U/d)", 0, 400, 0, 25),
      sliderInput("dietfat", "Dietary fat (% kcal)", 10, 40, 35, 1),

      hr(),
      sliderInput("years", "Simulate (years)", 4, 25, 12, 1),
      actionButton("go", "Simulate", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient profile", br(),
                 fluidRow(column(6, tableOutput("profile")),
                          column(6, plotOutput("depots", height = "300px"))),
                 hr(), htmlOutput("profile_note")),
        tabPanel("2 · Metreleptin PK", br(),
                 plotOutput("pk", height = "300px"),
                 plotOutput("pk24", height = "280px"),
                 htmlOutput("pk_note")),
        tabPanel("3 · The buffer & overflow", br(),
                 plotOutput("flux", height = "330px"),
                 plotOutput("ectopic", height = "300px"),
                 htmlOutput("flux_note")),
        tabPanel("4 · Leptin transducer", br(),
                 plotOutput("transducer", height = "330px"),
                 plotOutput("signal", height = "280px"),
                 htmlOutput("trans_note")),
        tabPanel("5 · Glycaemia", br(),
                 plotOutput("glyc", height = "380px"),
                 tableOutput("glyc_tab")),
        tabPanel("6 · Lipids & pancreatitis", br(),
                 plotOutput("lipid", height = "360px"),
                 tableOutput("panc_tab"),
                 htmlOutput("lipid_note")),
        tabPanel("7 · Liver", br(),
                 plotOutput("liver", height = "380px"),
                 htmlOutput("liver_note")),
        tabPanel("8 · Kidney, endocrine, body composition", br(),
                 plotOutput("other", height = "400px")),
        tabPanel("9 · Clinical endpoints", br(),
                 h4("Change from the day treatment started"),
                 tableOutput("endpoints"),
                 plotOutput("endplot", height = "300px")),
        tabPanel("10 · Scenario comparison", br(),
                 selectInput("scn_pick", "Scenarios", multiple = TRUE,
                             choices = names(scenarios),
                             selected = c("s1", "s6", "s16", "s18")),
                 actionButton("run_scn", "Run selected", class = "btn-primary"),
                 br(), br(), DTOutput("scn_tab"), plotOutput("scn_plot", height = "340px")),
        tabPanel("11 · Falsification tests", br(),
                 helpText(paste("Each table below is an attempt to break one",
                                "structural claim of the model. Test 2 also",
                                "records a structure that WAS broken and",
                                "discarded during development.")),
                 selectInput("inf_pick", "Test",
                             choices = c("1 occupancy x deficit" = "i1",
                                         "2 pair-fed clamp" = "i2",
                                         "3 TG/hepatic decoupling" = "i3",
                                         "4 PPARg needs a substrate" = "i4",
                                         "5 diet equivalence" = "i5",
                                         "6 treatment window" = "i6",
                                         "7 neutralising ADA" = "i7")),
                 actionButton("run_inf", "Run test", class = "btn-primary"),
                 br(), br(), DTOutput("inf_tab"), htmlOutput("inf_note")),
        tabPanel("12 · Biomarkers", br(),
                 plotOutput("bio", height = "420px"),
                 htmlOutput("bio_note"))
      )
    )
  )
)

## -----------------------------------------------------------------------------
##  SERVER
## -----------------------------------------------------------------------------
server <- function(input, output, session) {

  ## Selecting a phenotype rewrites the capacity sliders, because in this model
  ## the phenotype IS the capacity distribution.
  observeEvent(input$pheno, {
    p <- phenotypes[[input$pheno]]
    updateSliderInput(session, "fc1", value = p$FC1)
    updateSliderInput(session, "fc2", value = p$FC2)
    updateSliderInput(session, "fcv", value = p$FCV)
    updateSliderInput(session, "pread", value = p$PREAD)
  })

  base <- reactive({
    p <- pheno_param(phenotypes[[input$pheno]])
    p$FC1 <- input$fc1; p$FC2 <- input$fc2; p$FCV <- input$fcv
    p$PREAD <- input$pread
    baseline_pheno(p)
  })

  sim <- eventReactive(input$go, {
    b   <- base()
    end <- input$years*365
    t0  <- input$met_start*365
    dur <- if (input$met_stop) max(1, input$met_stopy*365 - t0) else end - t0

    ev_list <- list()
    if (input$met_on && input$met_dose > 0 && dur > 0)
      ev_list <- c(ev_list, list(met_sc(input$met_dose, 60, start = t0, dur = dur)))
    if (input$pio_on) ev_list <- c(ev_list, list(pio_oral(45, t0, end - t0)))
    if (input$vol_on) ev_list <- c(ev_list, list(vol_sc(300, t0, end - t0)))
    if (input$evi_on) ev_list <- c(ev_list, list(evi_iv(900, t0, ceiling((end - t0)/28))))
    if (input$fib_on) ev_list <- c(ev_list, list(fen_oral(145, t0, end - t0)))
    if (input$glp_on) ev_list <- c(ev_list, list(glp_sc(1, t0, end - t0)))
    if (input$met2_on) ev_list <- c(ev_list, list(metf_oral(1000, t0, end - t0)))
    if (input$sgl_on) ev_list <- c(ev_list, list(sgl_oral(10, t0, end - t0)))

    ps <- list()
    if (input$neut) ps$NEUT <- 1
    if (input$clamp) ps$EICLAMP <- b$init$EIS
    if (input$om3_on) { ps$OM3DOSE <- 1; ps$OM3T <- t0 }
    if (input$ins_dose > 0) { ps$INSDOSE <- input$ins_dose; ps$INST <- t0 }
    if (input$dietfat != 35) { ps$DIETT <- t0; ps$DIETFAT <- input$dietfat/100 }

    m <- zero_re(b$mod)
    if (length(ps)) m <- param(m, ps)

    trt <- if (length(ev_list))
      mrgsim_df(m, events = do.call(c, ev_list), end = end, delta = 1)
    else mrgsim_df(m, end = end, delta = 1)
    unt <- mrgsim_df(zero_re(b$mod), end = end, delta = 1)

    trt$arm <- "treated"; unt$arm <- "untreated"
    list(trt = trt, unt = unt, both = rbind(unt, trt), t0 = t0, end = end,
         base = b, mod = m, ev = ev_list)
  }, ignoreNULL = FALSE)

  yr <- function(d) d$time/365

  two_arm <- function(col, ylab, title, sub = NULL, logy = FALSE) {
    s <- sim(); d <- s$both
    d$y <- d[[col]]
    p <- ggplot(d, aes(time/365, y, colour = arm)) +
      geom_vline(xintercept = s$t0/365, linetype = 3, colour = "grey45") +
      geom_line(linewidth = 0.8) +
      scale_colour_manual(values = ARMCOL, name = NULL) +
      labs(x = "Years", y = ylab, title = title, subtitle = sub) + THEME
    if (logy) p <- p + scale_y_log10()
    p
  }

  ## ---- 1 patient profile -------------------------------------------------
  output$profile <- renderTable({
    s <- sim(); a <- s$unt[1, ]
    data.frame(
      Quantity = c("Total storage capacity (kg)", "Total fat mass (kg)",
                   "Subcutaneous fill (fraction)", "Leptin (ng/mL)",
                   "Leptin signal S(L)", "SIGNAL DEFICIT (1-S)",
                   "Energy intake (kcal/d)", "Influx J_in (g/d)",
                   "Adipose storage flux (g/d)", "OVERFLOW J_ov (g/d)",
                   "Hepatic fat fraction (%)", "Triglyceride (mg/dL)",
                   "Glucose (mg/dL)", "HbA1c (%)", "Insulin (uU/mL)",
                   "Insulin resistance index", "Adiponectin (ug/mL)",
                   "ALT (U/L)", "Fibrosis stage", "Proteinuria (g/d)"),
      Value = c(round(a$t_CAP, 2), round(a$t_FAT, 2), round(a$t_fill, 2),
                round(a$t_LTOT, 2), round(a$t_SD, 3), round(a$t_DEF, 3),
                round(a$EIS), round(a$t_JIN, 1), round(a$t_JST, 1),
                round(a$t_JOV, 1), round(a$t_HFF, 1), round(a$PTG),
                round(a$GLU), round(a$A1C, 2), round(a$INS, 1),
                round(a$IRX, 2), round(a$ADPN, 1), round(a$ALT),
                round(a$FIB, 2), round(a$PROT, 2)))
  }, striped = TRUE)

  output$depots <- renderPlot({
    s <- sim(); a <- s$unt[1, ]
    cap <- c(femorogluteal = 12*input$fc1*a$CAPM, `upper trunk/face` = 5*input$fc2*a$CAPM,
             visceral = 3*input$fcv)
    act <- c(a$ASC1, a$ASC2, a$AVIS)
    d <- data.frame(depot = rep(names(cap), 2),
                    kg = c(cap - act, act),
                    what = rep(c("headroom", "stored"), each = 3))
    ggplot(d, aes(depot, kg, fill = what)) +
      geom_col(width = 0.6) +
      scale_fill_manual(values = c(headroom = "#dfe6ee", stored = "#c9a227"), name = NULL) +
      labs(title = "Depot capacity and how full it is", y = "kg TG", x = NULL) +
      THEME + theme(axis.text.x = element_text(angle = 20, hjust = 1))
  })

  output$profile_note <- renderUI(HTML(paste0(
    "<p style='color:#444'>The capacity ceilings on the left are the ONLY ",
    "phenotype inputs. Leptin, hepatic fat, triglyceride and HbA1c shown here ",
    "are outputs of running those ceilings to steady state, not settings. ",
    "A depot that is nearly full has almost no headroom left, which is why the ",
    "storage flux collapses long before the fat mass reaches zero.</p>")))

  ## ---- 2 metreleptin PK --------------------------------------------------
  output$pk <- renderPlot({
    two_arm("t_LTOT", "Leptin, assay-equivalent (ng/mL)",
            "Total leptin: endogenous + metreleptin",
            "Daily dosing of a ~4-hour half-life drug: every point is near-trough")
  })

  output$pk24 <- renderPlot({
    s <- sim()
    if (!input$met_on || input$met_dose == 0) return(NULL)
    d0 <- max(s$t0 + 365, s$t0 + 30)
    out <- mrgsim_df(s$mod, events = met_sc(input$met_dose, 60, start = s$t0,
                                           dur = s$end - s$t0),
                     end = d0 + 1, delta = 0.02, recsort = 3)
    w <- out[out$time >= d0 & out$time <= d0 + 1, ]
    w$hour <- (w$time - d0)*24
    ggplot(w, aes(hour)) +
      geom_line(aes(y = t_LTOT, colour = "total"), linewidth = 0.9) +
      geom_line(aes(y = t_LEFF, colour = "receptor-available"), linewidth = 0.9) +
      geom_hline(yintercept = 4, linetype = 2, colour = "#b03060") +
      annotate("text", x = 1, y = 4.6, label = "transducer EC50", size = 3,
               colour = "#b03060", hjust = 0) +
      scale_colour_manual(values = c(total = "#1f6feb",
                                     `receptor-available` = "#2e8b57"), name = NULL) +
      labs(title = "One dosing interval at steady state",
           subtitle = sprintf("Mean 24-h signal S = %.3f; %.0f%% of the day above EC50",
                              mean(w$t_SD), 100*mean(w$t_LEFF > 4)),
           x = "Hours after dose", y = "ng/mL") + THEME
  })

  output$pk_note <- renderUI(HTML(paste0(
    "<p style='color:#444'>Target engagement is easy: at the approved dose the ",
    "receptor sees leptin above the transducer EC50 for essentially the whole ",
    "day in every phenotype. That is precisely why occupancy cannot explain ",
    "the difference in outcome between generalised and partial lipodystrophy ",
    "(tab 11, test 1). Neutralising antibodies split the two curves: total ",
    "leptin stays high while receptor-available leptin collapses.</p>")))

  ## ---- 3 buffer and overflow --------------------------------------------
  output$flux <- renderPlot({
    s <- sim()
    d <- s$both %>%
      select(time, arm, `influx J_in` = t_JIN, `stored` = t_JST,
             `oxidised` = t_OXW, `OVERFLOW J_ov` = t_JOV) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time/365, value, colour = name, linetype = arm)) +
      geom_vline(xintercept = s$t0/365, linetype = 3, colour = "grey45") +
      geom_line(linewidth = 0.75) +
      scale_colour_manual(values = c("influx J_in" = "#d2691e", stored = "#6b8e23",
                                     oxidised = "#4a90b8",
                                     "OVERFLOW J_ov" = "#b03060"), name = NULL) +
      scale_linetype_manual(values = c(untreated = 2, treated = 1), name = NULL) +
      labs(x = "Years", y = "g TG / day",
           title = "The buffer balance",
           subtitle = "J_ov = J_in - storage - oxidation; everything downstream follows this line") +
      THEME
  })

  output$ectopic <- renderPlot({
    s <- sim()
    d <- s$both %>%
      select(time, arm, `liver (g)` = LTG, `muscle (g)` = IMCL) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time/365, value, colour = arm)) +
      geom_vline(xintercept = s$t0/365, linetype = 3, colour = "grey45") +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = ARMCOL, name = NULL) +
      labs(x = "Years", y = NULL, title = "Ectopic pools fed by the overflow") + THEME
  })

  output$flux_note <- renderUI(HTML(paste0(
    "<p style='color:#444'>Storage flux is capacity-limited: it is the product ",
    "of adipocyte number (capacity), remaining headroom, substrate saturation ",
    "and insulin action. Note that insulin ACTION is nearly preserved in ",
    "insulin-resistant patients because the adipocyte sees a 4-5x higher ",
    "insulin concentration — so capacity, not signalling, is what limits ",
    "storage. Only two interventions move the overflow line: anything that ",
    "lowers J_in (metreleptin via appetite, a GLP-1 RA, a low-fat diet) and ",
    "anything that raises capacity (a PPARg agonist, if preadipocytes ",
    "remain).</p>")))

  ## ---- 4 leptin transducer ----------------------------------------------
  output$transducer <- renderPlot({
    s <- sim(); a <- s$unt[1, ]
    L <- seq(0.05, 60, length.out = 400)
    S <- L^2/(L^2 + 4^2)
    d <- data.frame(L, S)
    pts <- data.frame(
      L = c(a$t_LTOT, s$trt[nrow(s$trt), ]$t_LTOT),
      S = c(a$t_SD, s$trt[nrow(s$trt), ]$t_SD),
      what = c("baseline", "on treatment (trough)"))
    ggplot(d, aes(L, S)) +
      geom_line(linewidth = 0.9, colour = "#2e8b57") +
      geom_vline(xintercept = 4, linetype = 2, colour = "#b03060") +
      geom_point(data = pts, aes(L, S, shape = what), size = 3.5) +
      geom_segment(data = pts[1, ], aes(x = L, xend = L, y = S, yend = 1),
                   arrow = arrow(length = unit(0.15, "cm")), colour = "#b03060") +
      annotate("text", x = pts$L[1]*1.15, y = min(1, pts$S[1] + 0.5),
               label = sprintf("deficit = %.2f\n(the achievable gain)", 1 - pts$S[1]),
               hjust = 0, size = 3.4, colour = "#b03060") +
      scale_x_log10() +
      labs(x = "Total leptin (ng/mL, log scale)", y = "Hypothalamic signal S(L)",
           title = "The saturating transducer, and where this patient sits on it",
           subtitle = "The drug moves the patient rightwards; the DEFICIT is the vertical distance left to travel") +
      THEME
  })

  output$signal <- renderPlot({
    s <- sim()
    d <- s$both %>% select(time, arm, signal = t_SD, deficit = t_DEF,
                           occupancy = t_OCC) %>% pivot_longer(-c(time, arm))
    ggplot(d, aes(time/365, value, colour = name, linetype = arm)) +
      geom_vline(xintercept = s$t0/365, linetype = 3, colour = "grey45") +
      geom_line(linewidth = 0.75) +
      scale_colour_manual(values = c(signal = "#2e8b57", deficit = "#b03060",
                                     occupancy = "#1f6feb"), name = NULL) +
      scale_linetype_manual(values = c(untreated = 2, treated = 1), name = NULL) +
      labs(x = "Years", y = NULL,
           title = "Occupancy is a property of the dose; the deficit is a property of the patient") +
      THEME
  })

  output$trans_note <- renderUI(HTML(paste0(
    "<p style='color:#444'>Effect = occupancy x deficit. Because S(L) ",
    "saturates near 10 ng/mL, a patient whose endogenous leptin is already 8-12 ",
    "ng/mL has almost no vertical distance left to travel, however much drug ",
    "is given. Doubling the dose in that patient raises occupancy and leaves ",
    "the effect where it was (tab 11, test 1). The same saturation, ",
    "approached from above and combined with the SOCS3 brake, is why leptin ",
    "does nothing in common obesity.</p>")))

  ## ---- 5 glycaemia -------------------------------------------------------
  output$glyc <- renderPlot({
    s <- sim()
    d <- s$both %>% select(time, arm, `HbA1c (%)` = A1C, `glucose (mg/dL)` = GLU,
                           `insulin (uU/mL)` = t_INST, `beta-cell function` = BCF,
                           `insulin resistance` = IRX) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time/365, value, colour = arm)) +
      geom_vline(xintercept = s$t0/365, linetype = 3, colour = "grey45") +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = ARMCOL, name = NULL) +
      labs(x = "Years", y = NULL, title = "Glycaemic trajectory") + THEME
  })

  output$glyc_tab <- renderTable({
    s <- sim()
    i <- function(d, t) d[which.min(abs(d$time - t)), ]
    tt <- c(s$t0, s$t0 + 182, s$t0 + 365, s$end)
    do.call(rbind, lapply(tt, function(t) data.frame(
      year = round(t/365, 1),
      HbA1c_untreated = round(i(s$unt, t)$A1C, 2),
      HbA1c_treated   = round(i(s$trt, t)$A1C, 2),
      insulin_treated = round(i(s$trt, t)$t_INST, 1),
      IR_treated      = round(i(s$trt, t)$IRX, 2))))
  }, striped = TRUE)

  ## ---- 6 lipids ----------------------------------------------------------
  output$lipid <- renderPlot({
    s <- sim()
    d <- s$both
    ggplot(d, aes(time/365, PTG, colour = arm)) +
      geom_hline(yintercept = 1000, linetype = 2, colour = "#b03060") +
      annotate("text", x = 0.2, y = 1150, label = "pancreatitis threshold",
               size = 3.2, colour = "#b03060", hjust = 0) +
      geom_vline(xintercept = s$t0/365, linetype = 3, colour = "grey45") +
      geom_line(linewidth = 0.85) +
      scale_colour_manual(values = ARMCOL, name = NULL) +
      scale_y_log10() +
      labs(x = "Years", y = "Plasma triglyceride (mg/dL, log)",
           title = "Plasma triglyceride is tier 5 of the overflow cascade",
           subtitle = "What matters clinically is time above the threshold, not the mean") +
      THEME
  })

  output$panc_tab <- renderTable({
    s <- sim()
    i <- function(d, t) d[which.min(abs(d$time - t)), ]
    data.frame(
      arm = c("untreated", "treated"),
      TG_end = c(round(i(s$unt, s$end)$PTG), round(i(s$trt, s$end)$PTG)),
      days_above_1000 = c(round(i(s$unt, s$end)$TAT - i(s$unt, s$t0)$TAT),
                          round(i(s$trt, s$end)$TAT - i(s$trt, s$t0)$TAT)),
      cumulative_hazard = c(round(i(s$unt, s$end)$PANCH - i(s$unt, s$t0)$PANCH, 3),
                            round(i(s$trt, s$end)$PANCH - i(s$trt, s$t0)$PANCH, 3)))
  }, striped = TRUE)

  output$lipid_note <- renderUI(HTML(paste0(
    "<p style='color:#444'>An APOC3 or ANGPTL3-directed agent clears this ",
    "compartment without touching influx or capacity, so it wins on ",
    "time-above-threshold while leaving hepatic fat and HbA1c almost exactly ",
    "where they were. That dissociation is the model's sharpest testable ",
    "prediction (tab 11, test 3): if hepatic fat fell alongside triglyceride, ",
    "the tier structure would be wrong.</p>")))

  ## ---- 7 liver -----------------------------------------------------------
  output$liver <- renderPlot({
    s <- sim()
    d <- s$both %>% select(time, arm, `hepatic fat (%)` = t_HFF, `ALT (U/L)` = ALT,
                           `inflammation` = INFL, `fibrosis stage` = FIB) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time/365, value, colour = arm)) +
      geom_vline(xintercept = s$t0/365, linetype = 3, colour = "grey45") +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = ARMCOL, name = NULL) +
      labs(x = "Years", y = NULL,
           title = "The liver, and the only state in this model that does not fully reverse") +
      THEME
  })

  output$liver_note <- renderUI(HTML(paste0(
    "<p style='color:#444'>Hepatic fat, ALT and inflammation relax back within ",
    "months of removing the driver. Fibrosis regresses 67 times more slowly ",
    "than it forms, so its stage at the moment treatment starts sets the ",
    "CEILING of what can be recovered — which is the whole argument for early ",
    "initiation (tab 11, test 6). Nothing here is bistable: withdraw the drug ",
    "and every other state returns to where it began.</p>")))

  ## ---- 8 other organs ----------------------------------------------------
  output$other <- renderPlot({
    s <- sim()
    d <- s$both %>% select(time, arm, `proteinuria (g/d)` = PROT,
                           `eGFR` = EGFR, `androgen index` = ANDX,
                           `body weight (kg)` = t_BW, `fat mass (kg)` = t_FAT,
                           `adiponectin (ug/mL)` = ADPN) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time/365, value, colour = arm)) +
      geom_vline(xintercept = s$t0/365, linetype = 3, colour = "grey45") +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = ARMCOL, name = NULL) +
      labs(x = "Years", y = NULL, title = "Organ and body-composition consequences") + THEME
  })

  ## ---- 9 endpoints -------------------------------------------------------
  output$endpoints <- renderTable({
    s <- sim()
    i <- function(d, t) d[which.min(abs(d$time - t)), ]
    a <- i(s$trt, s$t0); b1 <- i(s$trt, s$t0 + 365); b2 <- i(s$trt, s$end)
    u2 <- i(s$unt, s$end)
    data.frame(
      Endpoint = c("HbA1c (%)", "Triglyceride (%)", "Hepatic fat (%rel)",
                   "Leptin (ng/mL)", "Energy intake (kcal/d)", "Fat mass (kg)",
                   "ALT (U/L)", "Fibrosis stage", "Proteinuria (g/d)"),
      `At 1 year` = c(round(b1$A1C - a$A1C, 2),
                      round(100*(b1$PTG - a$PTG)/a$PTG, 1),
                      round(100*(b1$t_HFF - a$t_HFF)/a$t_HFF, 1),
                      round(b1$t_LTOT - a$t_LTOT, 2),
                      round(b1$EIS - a$EIS),
                      round(b1$t_FAT - a$t_FAT, 2),
                      round(b1$ALT - a$ALT), round(b1$FIB - a$FIB, 2),
                      round(b1$PROT - a$PROT, 2)),
      `At end` = c(round(b2$A1C - a$A1C, 2),
                   round(100*(b2$PTG - a$PTG)/a$PTG, 1),
                   round(100*(b2$t_HFF - a$t_HFF)/a$t_HFF, 1),
                   round(b2$t_LTOT - a$t_LTOT, 2),
                   round(b2$EIS - a$EIS), round(b2$t_FAT - a$t_FAT, 2),
                   round(b2$ALT - a$ALT), round(b2$FIB - a$FIB, 2),
                   round(b2$PROT - a$PROT, 2)),
      `Untreated at end` = c(round(u2$A1C - a$A1C, 2),
                             round(100*(u2$PTG - a$PTG)/a$PTG, 1),
                             round(100*(u2$t_HFF - a$t_HFF)/a$t_HFF, 1),
                             round(u2$t_LTOT - a$t_LTOT, 2),
                             round(u2$EIS - a$EIS), round(u2$t_FAT - a$t_FAT, 2),
                             round(u2$ALT - a$ALT), round(u2$FIB - a$FIB, 2),
                             round(u2$PROT - a$PROT, 2)),
      check.names = FALSE)
  }, striped = TRUE)

  output$endplot <- renderPlot({
    s <- sim()
    d <- s$both %>% select(time, arm, HbA1c = A1C, `hepatic fat %` = t_HFF) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time/365, value, colour = arm)) +
      geom_vline(xintercept = s$t0/365, linetype = 3, colour = "grey45") +
      geom_line(linewidth = 0.9) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = ARMCOL, name = NULL) +
      labs(x = "Years", y = NULL, title = "Primary endpoints") + THEME
  })

  ## ---- 10 scenario comparison -------------------------------------------
  scn_run <- eventReactive(input$run_scn, {
    picks <- input$scn_pick
    if (!length(picks)) return(NULL)
    lapply(picks, function(k) scenarios[[k]]())
  })

  output$scn_tab <- renderDT({
    d <- scn_run(); if (is.null(d)) return(NULL)
    datatable(summarise_all(d), options = list(pageLength = 12, scrollX = TRUE),
              rownames = FALSE)
  })

  output$scn_plot <- renderPlot({
    d <- scn_run(); if (is.null(d)) return(NULL)
    all <- do.call(rbind, lapply(d, function(x)
      x[, c("time", "scenario", "A1C", "PTG", "t_HFF")]))
    dl <- all %>% rename(HbA1c = A1C, TG = PTG, `hepatic fat %` = t_HFF) %>%
      pivot_longer(-c(time, scenario))
    ggplot(dl, aes(time/365, value, colour = scenario)) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      labs(x = "Years", y = NULL, title = "Scenario comparison") +
      THEME + theme(legend.text = element_text(size = 7))
  })

  ## ---- 11 falsification tests -------------------------------------------
  inf_run <- eventReactive(input$run_inf, {
    switch(input$inf_pick,
           i1 = inference1_dose_vs_deficit(),
           i2 = inference2_pairfed(),
           i3 = inference3_decoupling(),
           i4 = inference4_substrate(),
           i5 = inference5_diet_equivalence(),
           i6 = inference6_window(),
           i7 = inference7_ada())
  })

  output$inf_tab <- renderDT({
    d <- inf_run(); if (is.null(d)) return(NULL)
    datatable(d, options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  output$inf_note <- renderUI({
    txt <- switch(input$inf_pick,
      i1 = paste("Occupancy rises across the four arms while the effect does",
                 "not. Doubling the dose in the near-normal-leptin patient",
                 "raises 24-hour occupancy and leaves HbA1c and hepatic fat",
                 "where they were. The deficit, not the dose, is the constraint."),
      i2 = paste("With food intake clamped at the untreated value, 57% of the",
                 "hepatic-fat benefit and 35% of the HbA1c benefit survive, so",
                 "the direct hepatic arm of leptin is real but smaller than the",
                 "appetite arm. An earlier version of the model, in which leptin",
                 "inhibited VLDL secretion directly, predicted the clamp would",
                 "make steatosis 71% WORSE; that structure was discarded because",
                 "suppressing the largest hepatic disposal route at unchanged",
                 "substrate delivery cannot lower the pool."),
      i3 = paste("Volanesorsen removes 75-80% of plasma triglyceride and moves",
                 "hepatic fat by 0%. Metreleptin moves both. The ratio column",
                 "is the discriminating statistic."),
      i4 = paste("Identical pioglitazone exposure in both arms. Capacity rises",
                 "only where preadipocytes exist. Nothing in the model says",
                 "'TZDs do not work in generalised lipodystrophy'."),
      i5 = paste("The diet is titrated, not guessed, so that influx matches the",
                 "drug arm exactly. Hepatic fat then matches almost exactly,",
                 "while triglyceride and HbA1c do not: the residual is the",
                 "leptin signal, which the diet arm actually LOSES as fat mass",
                 "falls further."),
      i6 = paste("Same endpoint hepatic fat and HbA1c whenever treatment",
                 "starts; different fibrosis. Start time sets the ceiling of",
                 "recovery, not its rate."),
      i7 = paste("Total leptin, which is what a clinical assay measures, stays",
                 "high while receptor-available leptin collapses. The result",
                 "looks exactly like non-adherence and is not."))
    HTML(paste0("<p style='color:#444'>", txt, "</p>"))
  })

  ## ---- 12 biomarkers -----------------------------------------------------
  output$bio <- renderPlot({
    s <- sim()
    d <- s$both %>% select(time, arm, `leptin (ng/mL)` = t_LTOT,
                           `adiponectin (ug/mL)` = ADPN, `NEFA (mmol/L)` = NEFA,
                           `APOC3 (rel)` = APOC3, `ANGPTL3 (rel)` = ANG3,
                           `capacity multiplier` = CAPM,
                           `neutralising titre` = NAB,
                           `islet lipid (rel)` = PLIP) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time/365, value, colour = arm)) +
      geom_vline(xintercept = s$t0/365, linetype = 3, colour = "grey45") +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = ARMCOL, name = NULL) +
      labs(x = "Years", y = NULL, title = "Mechanistic biomarkers") + THEME
  })

  output$bio_note <- renderUI(HTML(paste0(
    "<p style='color:#444'>Two of these are PD markers that a trial would ",
    "actually follow (leptin, APOC3 knockdown), two are mechanism read-outs ",
    "with no available therapy (adiponectin, islet lipid), and the ",
    "neutralising titre is the one that dissociates exposure from effect.</p>")))
}

shinyApp(ui, server)
