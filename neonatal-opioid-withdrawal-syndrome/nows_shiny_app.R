## =============================================================================
##  Neonatal opioid withdrawal syndrome (NOWS) — interactive QSP dashboard
##
##  The app is organised around the model's one claim: the quantity that drives
##  the disease (GAP = A - ITONE) is not the quantity that is scored, and not
##  the quantity that is dosed.  Every tab is a different way of putting those
##  three quantities next to each other.
##
##  Run:  shiny::runApp("nows_shiny_app.R")
##  Requires: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT, gridExtra
##  NOT FOR CLINICAL USE.
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

## The model definition lives in nows_mrgsolve_model.R.  Sourcing it builds
## `mod`, `sim_nows()`, `nows_endpoints()` and `scenarios`.
source("nows_mrgsolve_model.R")

PAL <- c(gap = "#c2410c", tone = "#1d4ed8", adapt = "#7c2d12",
         score = "#b45309", dose = "#0f766e", care = "#15803d",
         esc = "#6d28d9", ref = "#94a3b8")

theme_nows <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "#475569", size = 10),
          legend.position = "bottom")
}

## -----------------------------------------------------------------------------
## UI
## -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Neonatal opioid withdrawal syndrome — QSP simulator"),
  tags$p(style = "color:#475569;margin-top:-10px;",
         "GAP = A − ITONE.  A noradrenergic disease, treated with a µ-agonist, ",
         "scored with an instrument that measures neither.  Educational model — not for clinical use."),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Exposure"),
      selectInput("MDRUG", "Maternal opioid",
                  c("Methadone" = 1, "Buprenorphine" = 2,
                    "Short-acting (heroin/oxycodone)" = 3, "None" = 0),
                  selected = 1),
      sliderInput("MDOSE", "Maternal daily dose (mg/day)", 2, 220, 90, step = 2),
      sliderInput("GA", "Gestational age at birth (wk)", 30, 41.5, 39, step = 0.5),
      sliderInput("WT0", "Birth weight (kg)", 1.2, 4.6, 3.30, step = 0.05),
      checkboxInput("FBZD", "Benzodiazepine co-exposure", FALSE),
      sliderInput("FNIC", "Maternal cigarettes/day", 0, 30, 0, step = 1),

      hr(), h4("Care and assessment"),
      sliderInput("CARE", "Non-pharmacologic care bundle (0–1)", 0, 1, 0.55, step = 0.05),
      checkboxInput("BF", "Breastfeeding", FALSE),
      radioButtons("ESCMODE", "Assessment driving the decision",
                   c("Finnegan score" = 0, "Eat–Sleep–Console" = 1), selected = 0),
      sliderInput("THRSTART", "Finnegan treatment threshold", 5, 14, 8, step = 0.5),

      hr(), h4("Pharmacotherapy"),
      selectInput("TRTDRUG", "Treatment opioid",
                  c("Oral morphine" = 1, "Oral methadone" = 2,
                    "Sublingual buprenorphine" = 3), selected = 1),
      sliderInput("DINIT", "Starting dose (mg/kg/day morphine-eq)", 0.10, 0.80, 0.32, step = 0.02),
      sliderInput("FTARGET", "Score the protocol titrates toward", 3, 8, 6, step = 0.2),
      selectInput("TAPMODE", "Weaning rule",
                  c("−10%/day of CURRENT dose" = 0,
                    "−10%/day of STABILISATION dose" = 1,
                    "−10%/day of STANDARD dose" = 2,
                    "Track A (unobservable oracle)" = 3), selected = 1),
      sliderInput("TAPFRAC", "Daily wean fraction", 0.02, 0.40, 0.10, step = 0.01),
      sliderInput("CDOSE", "Clonidine (µg/kg/day)", 0, 15, 0, step = 1),
      sliderInput("PDOSE", "Phenobarbital maintenance (mg/kg/day)", 0, 8, 0, step = 0.5),
      sliderInput("PLOAD", "Phenobarbital load (mg/kg)", 0, 25, 0, step = 1),
      sliderInput("TMAX", "Simulation length (days)", 15, 60, 45, step = 5)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## ---------------------------------------------------------------
        tabPanel(
          "1 · Infant profile",
          br(),
          fluidRow(column(12, verbatimTextOutput("profile"))),
          plotOutput("p_profile", height = "430px"),
          helpText("The set-point A at birth, the fraction of it that is the ",
                   "durable (in-utero) pool, and the treatment the protocol ",
                   "ends up choosing for this infant.")
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "2 · Pharmacokinetics",
          br(),
          plotOutput("p_pk", height = "460px"),
          helpText("The in-utero opioid washing out and the treatment opioid ",
                   "coming in.  The washout is the clock that sets onset: the ",
                   "maternal DOSE only shifts this curve sideways.")
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "3 · The gap",
          br(),
          plotOutput("p_gap", height = "480px"),
          helpText("A (what the brain expects), ITONE (what it is getting), and ",
                   "the difference.  Below zero is over-sedation, which no ",
                   "Finnegan item scores.")
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "4 · Score vs gap",
          br(),
          plotOutput("p_score", height = "460px"),
          helpText("The same infant read two ways.  Divergence between the two ",
                   "panels is the model's central claim, and it is largest when ",
                   "phenobarbital is on board or the environment is arousing.")
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "5 · Eat–Sleep–Console",
          br(),
          plotOutput("p_esc", height = "460px"),
          helpText("The three functional items as continuous states.  Feeding ",
                   "is the only one that falls in BOTH directions — it is the ",
                   "over-treatment alarm the score does not have.")
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "6 · Dose and receptor",
          br(),
          plotOutput("p_dose", height = "460px"),
          helpText("Dose, µ occupancy and receptor availability.  Note how far ",
                   "occupancy moves for a given dose change at the top of the ",
                   "Emax curve versus near EC50.")
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "7 · Weaning rules",
          br(),
          plotOutput("p_wean", height = "420px"),
          DTOutput("t_wean"),
          helpText("Four rules run on the SAME infant.  They are not better and ",
                   "worse; they are points on one trade-off curve between ",
                   "treatment days and cumulative withdrawal burden.")
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "8 · Trade-off frontier",
          br(),
          plotOutput("p_front", height = "480px"),
          helpText("Sweeping the suppression target traces the exchange rate ",
                   "between days of treatment and ∫GAP⁺dt.  Suppressing very ",
                   "hard buys withdrawal burden back at the cost of days AND ",
                   "of over-sedation.")
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "9 · Maternal dose",
          br(),
          plotOutput("p_mdose", height = "460px"),
          DTOutput("t_mdose"),
          helpText("Maternal methadone dose changes WHEN withdrawal starts, ",
                   "not HOW BAD it is, because the adaptation map is saturated ",
                   "at every clinical dose.")
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "10 · Preterm experiment",
          br(),
          plotOutput("p_preterm", height = "460px"),
          DTOutput("t_preterm"),
          helpText("Preterm NOWS is milder for two reasons that cohort data ",
                   "always confound.  The third curve forces a 34-week infant ",
                   "to carry a term-level adaptation pool, isolating the ",
                   "pharmacokinetic half.")
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "11 · Adjuncts",
          br(),
          plotOutput("p_adj", height = "440px"),
          DTOutput("t_adj"),
          helpText("Clonidine fills the same GIRK conductance without ",
                   "re-driving the adaptation pool.  Phenobarbital lowers the ",
                   "score without lowering the gap — which is useful only when ",
                   "part of the score is not opioid withdrawal.")
        ),

        ## ---------------------------------------------------------------
        tabPanel(
          "12 · Endpoints",
          br(),
          DTOutput("t_end"),
          br(),
          plotOutput("p_growth", height = "360px")
        )
      )
    )
  )
)

## -----------------------------------------------------------------------------
## SERVER
## -----------------------------------------------------------------------------
server <- function(input, output, session) {

  base_par <- reactive({
    list(MDRUG = as.numeric(input$MDRUG), MDOSE = input$MDOSE,
         GA = input$GA, WT0 = input$WT0,
         FBZD = as.numeric(input$FBZD), FNIC = input$FNIC,
         CARE = input$CARE, BF = as.numeric(input$BF),
         ESCMODE = as.numeric(input$ESCMODE), THRSTART = input$THRSTART,
         TRTDRUG = as.numeric(input$TRTDRUG), DINIT = input$DINIT,
         FTARGET = input$FTARGET, TAPMODE = as.numeric(input$TAPMODE),
         TAPFRAC = input$TAPFRAC, CDOSE = input$CDOSE,
         PDOSE = input$PDOSE, PLOAD = input$PLOAD)
  })

  runp <- function(extra = list()) {
    p <- modifyList(base_par(), extra)
    do.call(sim_nows, c(list(tmax = input$TMAX * 24), p)) |>
      mutate(day = time / 24)
  }

  out <- reactive(runp())

  ## ---- 1. profile ----------------------------------------------------------
  output$profile <- renderPrint({
    o <- out(); e <- nows_endpoints(o)
    cat(sprintf(
      "Set-point at birth  A0 = %.3f   (durable pool %.3f, re-inducible %.3f)\n",
      e$A0, o$AD[1], o$AT[1]))
    cat(sprintf("Peak Finnegan %.1f on day %.1f   ·  first crossing of 8: day %s\n",
                e$peak_finnegan, e$day_of_peak,
                ifelse(is.na(e$onset_day), "never", sprintf("%.1f", e$onset_day))))
    cat(sprintf("Pharmacotherapy: %s   ·  peak dose %.2f mg/kg/day  ·  %.1f days\n",
                ifelse(e$treated, "yes", "no"), e$max_dose, e$treat_days))
    cat(sprintf("Cumulative morphine equivalents %.2f mg/kg  ·  weight gain %.0f g\n",
                e$cum_morphine, e$weight_gain_g))
    cat(sprintf("Withdrawal burden ∫GAP+dt = %.1f   ·  over-sedation ∫GAP−dt = %.1f\n",
                e$auc_gap, e$auc_sedation))
  })

  output$p_profile <- renderPlot({
    o <- out()
    d <- o |> select(day, A = A_out, ITONE = ITONE_o, Finnegan = FNAS_S,
                     Dose = DOSE) |>
      pivot_longer(-day)
    ggplot(d, aes(day, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = c(A = PAL["adapt"], ITONE = PAL["tone"],
                                     Finnegan = PAL["score"], Dose = PAL["dose"]),
                          guide = "none") +
      labs(x = "postnatal day", y = NULL,
           title = "The four quantities this model keeps apart") +
      theme_nows()
  })

  ## ---- 2. PK ---------------------------------------------------------------
  output$p_pk <- renderPlot({
    o <- out()
    d <- o |> select(day, `in-utero methadone` = CD_o,
                     `in-utero buprenorphine` = CB_o,
                     `treatment morphine` = CM_o, `M6G` = CG_o,
                     `clonidine` = CC_o) |>
      pivot_longer(-day) |> filter(value > 1e-4)
    ggplot(d, aes(day, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      scale_y_log10() +
      labs(x = "postnatal day", y = "plasma concentration (ng/mL, log)",
           colour = NULL,
           title = "The washout that opens the gap, and what is given to close it") +
      theme_nows()
  })

  ## ---- 3. gap --------------------------------------------------------------
  output$p_gap <- renderPlot({
    o <- out()
    ggplot(o, aes(day)) +
      geom_hline(yintercept = 0, colour = "grey50", linetype = 2) +
      geom_ribbon(aes(ymin = pmin(GAP_o, 0), ymax = 0), fill = "#93c5fd", alpha = 0.5) +
      geom_ribbon(aes(ymin = 0, ymax = pmax(GAP_o, 0)), fill = "#fdba74", alpha = 0.6) +
      geom_line(aes(y = A_out, colour = "A (set-point)"), linewidth = 0.9) +
      geom_line(aes(y = ITONE_o, colour = "ITONE (inhibitory tone)"), linewidth = 0.9) +
      geom_line(aes(y = GAP_o, colour = "GAP"), linewidth = 1.1) +
      scale_colour_manual(values = c("A (set-point)" = PAL[["adapt"]],
                                     "ITONE (inhibitory tone)" = PAL[["tone"]],
                                     "GAP" = PAL[["gap"]])) +
      labs(x = "postnatal day", y = "fraction of maximal GIRK inhibition",
           colour = NULL,
           title = "The gap",
           subtitle = "orange = withdrawal drive · blue = over-sedation, which nothing scores") +
      theme_nows()
  })

  ## ---- 4. score vs gap -----------------------------------------------------
  output$p_score <- renderPlot({
    o <- out()
    a <- ggplot(o, aes(day, FNAS_S)) +
      geom_hline(yintercept = 8, linetype = 2, colour = "grey50") +
      geom_line(colour = PAL[["score"]], linewidth = 1) +
      labs(x = NULL, y = "smoothed Finnegan score",
           title = "What is measured") + theme_nows()
    b <- ggplot(o, aes(day, GAP_o)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
      geom_line(colour = PAL[["gap"]], linewidth = 1) +
      labs(x = "postnatal day", y = "GAP",
           title = "What is happening") + theme_nows()
    gridExtra::grid.arrange(a, b, ncol = 1)
  })

  ## ---- 5. ESC --------------------------------------------------------------
  output$p_esc <- renderPlot({
    o <- out()
    d <- o |> select(day, Eat = EAT, Sleep = SLP, Console = CONS) |>
      pivot_longer(-day)
    ggplot(d, aes(day, value, colour = name)) +
      geom_hline(yintercept = 0.55, linetype = 2, colour = "grey50") +
      geom_line(linewidth = 1) +
      ylim(0, 1) +
      labs(x = "postnatal day", y = "function preserved (0–1)", colour = NULL,
           title = "Eat · Sleep · Console",
           subtitle = "dashed line = the functional-failure threshold") +
      theme_nows()
  })

  ## ---- 6. dose and receptor ------------------------------------------------
  output$p_dose <- renderPlot({
    o <- out()
    d <- o |> select(day, `dose (mg/kg/day)` = DOSE,
                     `µ effect EMU` = EMU_o,
                     `α2 effect EA2` = EA2_o,
                     `receptor availability` = RMU) |>
      pivot_longer(-day)
    ggplot(d, aes(day, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "postnatal day", y = NULL,
           title = "Dose, occupancy, and the receptor pool") +
      theme_nows() + theme(legend.position = "none")
  })

  ## ---- 7. weaning rules ----------------------------------------------------
  wean_runs <- reactive({
    modes <- c("−10%/day of current" = 0, "−10%/day of stabilisation" = 1,
               "−10%/day of standard" = 2, "A-tracking oracle" = 3)
    bind_rows(lapply(names(modes), function(nm) {
      o <- runp(list(TAPMODE = unname(modes[nm])))
      o$rule <- nm; o
    }))
  })

  output$p_wean <- renderPlot({
    d <- wean_runs()
    ggplot(d, aes(day, DOSE, colour = rule)) +
      geom_line(linewidth = 0.95) +
      labs(x = "postnatal day", y = "dose (mg/kg/day morphine-eq)", colour = NULL,
           title = "The same infant, four weaning rules") +
      theme_nows()
  })

  output$t_wean <- renderDT({
    d <- wean_runs()
    r <- d |> group_by(rule) |> group_modify(~ nows_endpoints(.x)) |> ungroup() |>
      transmute(rule, treat_days = round(treat_days, 1),
                stop_day = round(stop_day, 1),
                cum_morphine = round(cum_morphine, 2),
                withdrawal_burden = round(auc_gap, 1),
                over_sedation = round(auc_sedation, 1),
                weight_gain_g = round(weight_gain_g))
    datatable(r, rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- 8. frontier ---------------------------------------------------------
  output$p_front <- renderPlot({
    tgt <- seq(3, 7.8, by = 0.6)
    d <- bind_rows(lapply(tgt, function(g) {
      o <- runp(list(FTARGET = g))
      cbind(target = g, nows_endpoints(o))
    }))
    ggplot(d, aes(auc_gap, treat_days)) +
      geom_path(colour = PAL[["ref"]]) +
      geom_point(aes(size = auc_sedation, colour = target)) +
      geom_text(aes(label = sprintf("%.1f", target)), vjust = -1.1, size = 3) +
      scale_colour_gradient(low = PAL[["gap"]], high = PAL[["tone"]]) +
      labs(x = "cumulative withdrawal burden  ∫GAP⁺dt",
           y = "days of pharmacotherapy",
           size = "over-sedation", colour = "score target",
           title = "The exchange rate nobody sets explicitly",
           subtitle = "each point is a suppression target; point size is cumulative over-sedation") +
      theme_nows()
  })

  ## ---- 9. maternal dose ----------------------------------------------------
  mdose_runs <- reactive({
    doses <- c(40, 70, 100, 130, 160)
    bind_rows(lapply(doses, function(md) {
      o <- runp(list(MDOSE = md, MDRUG = 1)); o$mdose <- md; o
    }))
  })

  output$p_mdose <- renderPlot({
    d <- mdose_runs()
    ggplot(d, aes(day, FNAS_S, colour = factor(mdose))) +
      geom_hline(yintercept = 8, linetype = 2, colour = "grey50") +
      geom_line(linewidth = 0.9) +
      labs(x = "postnatal day", y = "smoothed Finnegan score",
           colour = "maternal methadone (mg/day)",
           title = "Maternal dose moves the curve sideways, not upwards") +
      theme_nows()
  })

  output$t_mdose <- renderDT({
    d <- mdose_runs()
    r <- d |> group_by(mdose) |> group_modify(~ nows_endpoints(.x)) |> ungroup() |>
      transmute(`maternal mg/day` = mdose, A0 = round(A0, 3),
                `peak score` = round(peak_finnegan, 1),
                `onset (day)` = round(onset_day, 2),
                `treatment days` = round(treat_days, 1))
    datatable(r, rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- 10. preterm ---------------------------------------------------------
  preterm_runs <- reactive({
    cfg <- list(
      list(nm = "term 39 wk", p = list(GA = 39, WT0 = 3.30)),
      list(nm = "preterm 34 wk (own adaptation)", p = list(GA = 34, WT0 = 2.10)),
      list(nm = "preterm 34 wk with TERM adaptation",
           p = list(GA = 34, WT0 = 2.10, ADFORCE = 0.578)))
    bind_rows(lapply(cfg, function(c) { o <- runp(c$p); o$arm <- c$nm; o }))
  })

  output$p_preterm <- renderPlot({
    d <- preterm_runs()
    ggplot(d, aes(day, FNAS_S, colour = arm)) +
      geom_hline(yintercept = 8, linetype = 2, colour = "grey50") +
      geom_line(linewidth = 0.95) +
      labs(x = "postnatal day", y = "smoothed Finnegan score", colour = NULL,
           title = "Two reasons preterm NOWS is milder, separated") +
      theme_nows()
  })

  output$t_preterm <- renderDT({
    d <- preterm_runs()
    r <- d |> group_by(arm) |> group_modify(~ nows_endpoints(.x)) |> ungroup() |>
      transmute(arm, A0 = round(A0, 3), `peak score` = round(peak_finnegan, 1),
                `treatment days` = round(treat_days, 1),
                `peak dose` = round(max_dose, 2),
                `cum morphine mg/kg` = round(cum_morphine, 2))
    datatable(r, rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- 11. adjuncts --------------------------------------------------------
  adj_runs <- reactive({
    cfg <- list(
      list(nm = "morphine alone", p = list(CDOSE = 0, PDOSE = 0, PLOAD = 0)),
      list(nm = "+ clonidine 6 µg/kg/day", p = list(CDOSE = 6)),
      list(nm = "+ clonidine 12 µg/kg/day", p = list(CDOSE = 12)),
      list(nm = "+ phenobarbital", p = list(PDOSE = 5, PLOAD = 20)))
    bind_rows(lapply(cfg, function(c) { o <- runp(c$p); o$arm <- c$nm; o }))
  })

  output$p_adj <- renderPlot({
    d <- adj_runs() |>
      select(day, arm, `Finnegan score` = FNAS_S, `GAP` = GAP_o) |>
      pivot_longer(c(`Finnegan score`, GAP))
    ggplot(d, aes(day, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "postnatal day", y = NULL, colour = NULL,
           title = "Which adjunct moves the score, and which moves the gap") +
      theme_nows()
  })

  output$t_adj <- renderDT({
    d <- adj_runs()
    r <- d |> group_by(arm) |> group_modify(~ nows_endpoints(.x)) |> ungroup() |>
      transmute(arm, `treatment days` = round(treat_days, 1),
                `cum morphine mg/kg` = round(cum_morphine, 2),
                `withdrawal burden` = round(auc_gap, 1),
                `weight gain (g)` = round(weight_gain_g))
    datatable(r, rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- 12. endpoints -------------------------------------------------------
  output$t_end <- renderDT({
    e <- nows_endpoints(out())
    r <- data.frame(endpoint = names(e), value = unlist(lapply(e, function(x)
      if (is.numeric(x)) round(x, 3) else as.character(x))))
    datatable(r, rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  output$p_growth <- renderPlot({
    o <- out()
    ggplot(o, aes(day, WT * 1000)) +
      geom_line(colour = PAL[["care"]], linewidth = 1) +
      labs(x = "postnatal day", y = "weight (g)",
           title = "Growth — the endpoint that reads the gap from the other side") +
      theme_nows()
  })
}

shinyApp(ui, server)
