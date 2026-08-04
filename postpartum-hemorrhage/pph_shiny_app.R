## =====================================================================
##  PRIMARY POSTPARTUM HAEMORRHAGE — interactive QSP dashboard (Shiny)
##  Front-end for pph_mrgsolve_model.R
##
##  The app is organised around the model's thesis rather than around its
##  compartments: the left panel is "what kind of uterus is this, and what did
##  you do, and when", and the tabs answer, in order,
##
##    1  Patient & bleeding      — the flux, the store, and the clock
##    2  Uterine mechanics       — tone, the Poiseuille valve, the ceiling
##    3  Drug PK / receptor      — concentrations and the OTR that vanishes
##    4  Haemostasis             — the shear gate, clot, maturation, lysis
##    5  Resuscitation & milieu  — products, calcium, temperature, acid-base
##    6  Clinical endpoints      — loss, shock index, oxygen debt, organ risk
##    7  Scenario comparison     — any number of saved arms side by side
##    8  Severity sweep          — the escalation ladder and the knife edge
##    9  Population              — loss distribution and the refractory tail
##
##  NOT FOR CLINICAL USE. Educational / research QSP model.
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

## The model is defined in pph_mrgsolve_model.R; sourcing it compiles `pph` and
## defines the dosing/parameter builders.  The scenario functions there are only
## defined, not run, so sourcing is cheap.
if (!exists("pph")) source("pph_mrgsolve_model.R")

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#EEEEEE", colour = NA),
        legend.position = "bottom")

## ---------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Primary postpartum haemorrhage — QSP simulator"),
  p(em(paste("A flux problem against a store that only looks large.",
             "Term uterine blood flow 750 mL/min; the whole volume expansion",
             "of pregnancy is 2.8 minutes of it."))),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("The uterus"),
      sliderInput("ATN", "Atony severity (ATN)", 0, 1, 0.74, step = 0.01),
      helpText("0 normal · 0.6 moderate · 0.74 severe · 0.88 refractory · 0.94 catastrophic"),
      sliderInput("OTR0", "Oxytocin receptor availability at delivery",
                  0.1, 1, 1.0, step = 0.01),
      helpText("1.0 = agonist-naive; ~0.4 after a long augmented labour"),
      sliderInput("UBF0", "Uterine blood flow (mL/min)", 400, 1000, 750, step = 10),

      h4("The other three Ts"),
      sliderInput("KTR", "Trauma: laceration leak at normal BP (mL/min)",
                  0, 400, 0, step = 10),
      numericInput("TREPAIR", "Laceration repaired at (min, blank = never)",
                   value = NA, min = 0, max = 240),
      sliderInput("RET0", "Tissue: retained fraction", 0, 0.6, 0, step = 0.05),
      numericInput("TMANUAL", "Manual removal at (min)", value = NA),
      sliderInput("FIB0", "Thrombin: starting fibrinogen (g/L)",
                  1.0, 6.0, 4.5, step = 0.1),

      h4("Uterotonics (time of administration, min)"),
      checkboxInput("useoxy", "Oxytocin 5 IU IV + 40 IU / 4 h", TRUE),
      numericInput("toxy", NULL, value = 5, min = 0, max = 240),
      sliderInput("oxymult", "Oxytocin dose multiplier", 0.5, 4, 1, step = 0.5),
      checkboxInput("usecbt", "Carbetocin 100 ug IV", FALSE),
      numericInput("tcbt", NULL, value = 5),
      checkboxInput("useerg", "Ergometrine 0.5 mg IM", FALSE),
      numericInput("terg", NULL, value = 15),
      checkboxInput("usepgf", "Carboprost 250 ug IM (x2, 15 min apart)", FALSE),
      numericInput("tpgf", NULL, value = 20),
      checkboxInput("usemso", "Misoprostol 800 ug SL", FALSE),
      numericInput("tmso", NULL, value = 20),
      checkboxInput("usemass", "Uterine massage / bimanual compression", TRUE),
      numericInput("tmass", NULL, value = 4),

      h4("Antifibrinolytic"),
      checkboxInput("usetxa", "Tranexamic acid 1 g IV (repeat at +30 min)", TRUE),
      sliderInput("ttxa", "TXA given at (min)", 0, 180, 10, step = 5),

      h4("Mechanical / surgical"),
      checkboxInput("usebal", "Intrauterine balloon tamponade", FALSE),
      sliderInput("tbal", "Balloon at (min)", 0, 180, 30, step = 5),
      checkboxInput("useaoc", "Aortic compression (bridge)", FALSE),
      sliderInput("taoc", "Aortic compression from (min)", 0, 60, 8, step = 1),
      sliderInput("taoce", "... until (min)", 5, 120, 35, step = 5),
      checkboxInput("useual", "Uterine artery ligation", FALSE),
      numericInput("tual", NULL, value = 50),
      checkboxInput("usehyst", "Hysterectomy", FALSE),
      numericInput("thyst", NULL, value = 60),

      h4("Resuscitation"),
      selectInput("resus", NULL,
                  c("none", "goal-directed (1:1 + Fg + Ca + warmer)",
                    "crystalloid-first (RBC only, no warmer)"),
                  selected = "goal-directed (1:1 + Fg + Ca + warmer)"),
      checkboxInput("warmer", "Fluid warmer", TRUE),
      checkboxInput("calcium", "Calcium replacement", TRUE),
      sliderInput("maxrbc", "Red-cell units available", 0, 20, 10, step = 1),

      hr(),
      numericInput("tend", "Simulate to (min)", value = 240, min = 30, max = 480),
      actionButton("save", "Save this arm for comparison", class = "btn-primary"),
      actionButton("clear", "Clear saved arms")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient & bleeding",
                 fluidRow(column(12, htmlOutput("verdict"))),
                 plotOutput("p_bleed", height = "600px"),
                 tableOutput("t_milestones")),
        tabPanel("2 · Uterine mechanics", plotOutput("p_uterus", height = "700px"),
                 htmlOutput("note_uterus")),
        tabPanel("3 · Drug PK & receptor", plotOutput("p_pk", height = "700px"),
                 htmlOutput("note_pk")),
        tabPanel("4 · Haemostasis", plotOutput("p_clot", height = "700px"),
                 htmlOutput("note_clot")),
        tabPanel("5 · Resuscitation & milieu", plotOutput("p_resus", height = "700px"),
                 tableOutput("t_products")),
        tabPanel("6 · Clinical endpoints", plotOutput("p_end", height = "620px"),
                 tableOutput("t_end")),
        tabPanel("7 · Scenario comparison", plotOutput("p_cmp", height = "520px"),
                 tableOutput("t_cmp")),
        tabPanel("8 · Severity sweep", plotOutput("p_sweep", height = "560px"),
                 htmlOutput("note_sweep")),
        tabPanel("9 · Population", plotOutput("p_pop", height = "520px"),
                 tableOutput("t_pop"))
      )
    )
  )
)

## ---------------------------------------------------------------------
##  SERVER
## ---------------------------------------------------------------------
server <- function(input, output, session) {

  na_never <- function(x) if (is.null(x) || is.na(x)) 1e9 else x

  build <- reactive({
    pars <- list(
      ATN = input$ATN, OTR0 = input$OTR0, UBF0 = input$UBF0,
      KTR = input$KTR, RET0 = input$RET0, FIB0 = input$FIB0,
      TREPAIR = na_never(input$TREPAIR), TMANUAL = na_never(input$TMANUAL),
      TMASS = if (input$usemass) input$tmass else 1e9,
      TMASSE = if (input$usemass) 90 else 1e9,
      TOXYI = if (input$useoxy) input$toxy else 1e9,
      TOXYE = if (input$useoxy) input$toxy + 240 else 1e9,
      ROXYI = 0.1667 * input$oxymult,
      TBAL = if (input$usebal) input$tbal else 1e9,
      TAOC = if (input$useaoc) input$taoc else 1e9,
      TAOCE = if (input$useaoc) input$taoce else 1e9,
      TUAL = if (input$useual) na_never(input$tual) else 1e9,
      THYST = if (input$usehyst) na_never(input$thyst) else 1e9,
      MAXRBC = input$maxrbc,
      WARMER = as.numeric(input$warmer)
    )

    pars <- switch(input$resus,
      "none" = c(pars, list(TRESUS = 1e9)),
      "goal-directed (1:1 + Fg + Ca + warmer)" =
        c(pars, list(TRESUS = 12, RCRYST = 12, TCRYSTE = 40, HBTRIG = 7,
                     MTPLOSS = 1500, QSTOP = 25, RBCRATE = 0.2, FFPRATIO = 1,
                     PLTTRIG = 75, PLTRATE = 0.05, FIBTRIG = 2, RFGC = 0.4,
                     RCA = if (input$calcium) 0.6 else 0, CATRIG = 1.0)),
      "crystalloid-first (RBC only, no warmer)" =
        c(pars, list(TRESUS = 12, RCRYST = 70, TCRYSTE = 90, HBTRIG = 7,
                     MTPLOSS = 1500, QSTOP = 25, RBCRATE = 0.2, FFPRATIO = 0,
                     FIBTRIG = 0, RCA = if (input$calcium) 0.6 else 0,
                     CATRIG = 1.0))
    )

    e <- ev(time = 0, amt = 0, cmt = "OXY")
    if (input$useoxy)  e <- e + ev(time = input$toxy, amt = 5 * input$oxymult, cmt = "OXY")
    if (input$usecbt)  e <- e + ev(time = na_never(input$tcbt), amt = 100, cmt = "CBT")
    if (input$useerg)  e <- e + ev(time = na_never(input$terg), amt = 0.5, cmt = "ERD")
    if (input$usepgf)  e <- e + ev(time = na_never(input$tpgf), amt = 250, cmt = "PGD") +
                                ev(time = na_never(input$tpgf) + 15, amt = 250, cmt = "PGD")
    if (input$usemso)  e <- e + ev(time = na_never(input$tmso), amt = 800, cmt = "MSD")
    if (input$usetxa)  e <- e + ev(time = input$ttxa, amt = 1000, cmt = "TX1") +
                                ev(time = input$ttxa + 30, amt = 1000, cmt = "TX1")
    list(pars = pars, ev = e)
  })

  sim <- reactive({
    b <- build()
    zero_re(param(pph, b$pars)) %>%
      mrgsim(events = b$ev, end = input$tend, delta = 0.5) %>%
      as_tibble()
  })

  ## ---- saved arms -------------------------------------------------
  saved <- reactiveVal(list())
  observeEvent(input$save, {
    lab <- sprintf("ATN %.2f · OTR %.2f%s%s%s", input$ATN, input$OTR0,
                   if (input$usetxa) sprintf(" · TXA %g", input$ttxa) else " · no TXA",
                   if (input$usebal) sprintf(" · balloon %g", input$tbal) else "",
                   if (input$resus == "none") " · no resus" else "")
    s <- saved(); s[[lab]] <- sim(); saved(s)
  })
  observeEvent(input$clear, saved(list()))

  ## ---- helpers ----------------------------------------------------
  endpoints <- function(d) {
    dead <- d %>% filter(Lethal > 0) %>% slice(1)
    tibble(
      `Blood loss (mL)` = round(max(d$LOSS)),
      `Outcome` = if (nrow(dead)) sprintf("death at %.0f min", dead$time[1])
                  else if (tail(d$Qbleed, 1) < 25) "haemostasis achieved"
                  else "still bleeding at end of simulation",
      `Hb nadir (g/dL)` = round(min(d$Hb), 1),
      `Fibrinogen nadir (g/L)` = round(min(d$Fib), 2),
      `MAP nadir (mmHg)` = round(min(d$MAP)),
      `Peak shock index` = round(max(d$ShockIx), 2),
      `Tone at end` = round(tail(d$TONE, 1), 2),
      `RBC units` = round(max(d$RBCU), 1),
      `FFP units` = round(max(d$FFPU), 1),
      `Fibrinogen given (g)` = round(max(d$FGCU), 1),
      `iCa nadir (mmol/L)` = round(min(d$CAI), 2),
      `Temperature nadir (C)` = round(min(d$TMP), 1),
      `O2 debt (mL)` = round(max(d$ODEBT)),
      `min with MAP<50` = round(max(d$TSEV)),
      `min with MAP<65` = round(max(d$TAKI)),
      `min with Fib<2` = round(max(d$XSEV))
    )
  }

  long <- function(d, cols) {
    d %>% select(time, all_of(names(cols))) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(cols[name], levels = unname(cols)))
  }

  ## ---- tab 1 ------------------------------------------------------
  output$verdict <- renderUI({
    e <- endpoints(sim())
    col <- if (grepl("death", e$Outcome)) "#B71C1C"
           else if (grepl("still", e$Outcome)) "#E65100" else "#1B5E20"
    HTML(sprintf(
      "<div style='padding:10px;border-left:6px solid %s;background:#FAFAFA'>
       <b style='font-size:16px'>%s &mdash; %s mL</b><br/>
       Haemoglobin nadir %s g/dL, fibrinogen nadir %s g/L, MAP nadir %s mmHg,
       accumulated oxygen debt %s mL O<sub>2</sub> (lethal 8400).</div>",
      col, e$Outcome, format(e$`Blood loss (mL)`, big.mark = ","),
      e$`Hb nadir (g/dL)`, e$`Fibrinogen nadir (g/L)`, e$`MAP nadir (mmHg)`,
      format(e$`O2 debt (mL)`, big.mark = ",")))
  })

  output$p_bleed <- renderPlot({
    cols <- c(Qbleed = "leak (mL/min)", LOSS = "cumulative loss (mL)",
              MAP = "MAP (mmHg)", HR = "heart rate (/min)",
              Hb = "haemoglobin (g/dL)", DO2 = "oxygen delivery (mL/min)")
    ggplot(long(sim(), cols), aes(time, value)) +
      geom_line(linewidth = 0.9, colour = "#B71C1C") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "minutes after delivery of the placenta", y = NULL,
           title = "The flux, and what it does to the store") + THEME
  })

  output$t_milestones <- renderTable({
    d <- sim()
    m <- c(500, 1000, 1500, 2500, 3150)
    tibble(`Milestone (mL)` = m,
           `Reached at (min)` = sapply(m, function(x) {
             i <- which(d$LOSS >= x)[1]
             if (is.na(i)) "not reached" else sprintf("%.1f", d$time[i])
           }))
  })

  ## ---- tab 2 ------------------------------------------------------
  output$p_uterus <- renderPlot({
    cols <- c(TONE = "myometrial tone", ToneTgt = "tone target (drugs+milieu)",
              Drive = "total contractile drive", Patency = "patency (1-tone)^4",
              Quterus = "uterine leak (mL/min)", Uflow = "uterine perfusion (mL/min)")
    ggplot(long(sim(), cols), aes(time, value)) +
      geom_line(linewidth = 0.9, colour = "#880E4F") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "minutes", y = NULL,
           title = "The living ligature: tone, the Poiseuille valve, and the leak it sets") +
      THEME
  })
  output$note_uterus <- renderUI(HTML(
    "<p><em>Patency goes as the fourth power of (1&minus;tone), so the band between
     a firm and a boggy uterus is a few percent of contractile capacity. Tone at
     0.91 leaks 0.05 mL/min; at 0.50, 47 mL/min; at 0.20, 307 mL/min. Note that
     the drive can be high while the tone target stays low &mdash; that is the
     ceiling (CAP) imposed by atony itself.</em></p>"))

  ## ---- tab 3 ------------------------------------------------------
  output$p_pk <- renderPlot({
    cols <- c(Oxytocin = "oxytocin (IU/L)", Carbopr = "carboprost (ug/L)",
              Ergomet = "ergometrine (ug/L)", Misopr = "misoprostol (ug/L)",
              TXA = "tranexamic acid (mg/L)", OTR = "available OTR fraction")
    ggplot(long(sim(), cols), aes(time, value)) +
      geom_line(linewidth = 0.9, colour = "#01579B") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "minutes", y = NULL,
           title = "Concentrations, and the receptor that disappears while you use it") +
      THEME
  })
  output$note_pk <- renderUI(HTML(
    "<p><em>Watch the OTR panel. A standard infusion takes receptor availability
     from 1.00 to about 0.09 over two hours, which is why raising the oxytocin
     dose recovers little and crossing to an FP, &alpha;<sub>1</sub> or EP
     receptor agonist recovers a lot. Set the receptor slider to 0.38 (an
     augmented labour) and compare the dose multiplier against adding
     carboprost.</em></p>"))

  ## ---- tab 4 ------------------------------------------------------
  output$p_clot <- renderPlot({
    cols <- c(ShearGate = "shear gate 1/(1+(Q/80)^2)", CLT = "placental-bed clot",
              MAT = "thrombus organisation", Seal = "effective seal",
              Plasmin = "plasmin (baseline 1)", TXAprot = "clot protection by TXA")
    ggplot(long(sim(), cols), aes(time, value)) +
      geom_line(linewidth = 0.9, colour = "#004D40") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "minutes", y = NULL,
           title = "Coagulation is locked out until the mechanics work") + THEME
  })
  output$note_clot <- renderUI(HTML(
    "<p><em>The shear gate is the lock: at 168 mL/min it sits at 0.19 and the
     placental-bed clot plateaus around 0.22, so no transfusion strategy can
     produce haemostasis. Reduce the leak &mdash; with tone or with a balloon &mdash;
     and the gate opens, the clot forms, and only then does it begin to organise
     into a durable seal. That organisation window is what every mechanical
     manoeuvre is really buying, and what TXA protects.</em></p>"))

  ## ---- tab 5 ------------------------------------------------------
  output$p_resus <- renderPlot({
    cols <- c(Fib = "fibrinogen (g/L)", Plt = "platelets (10^9/L)",
              Factors = "pooled factor activity", CAI = "ionised calcium (mmol/L)",
              TMP = "core temperature (C)", pH = "arterial pH")
    d <- long(sim(), cols)
    thr <- tibble(name = factor(c("fibrinogen (g/L)", "platelets (10^9/L)",
                                  "pooled factor activity",
                                  "ionised calcium (mmol/L)",
                                  "core temperature (C)", "arterial pH"),
                                levels = unname(cols)),
                  y = c(2.0, 50, 0.30, 1.0, 36, 7.20))
    ggplot(d, aes(time, value)) +
      geom_hline(data = thr, aes(yintercept = y), linetype = 2, colour = "#B71C1C") +
      geom_line(linewidth = 0.9, colour = "#0D47A1") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "minutes", y = NULL,
           title = "What the resuscitation costs (dashed = critical value)") + THEME
  })
  output$t_products <- renderTable({
    d <- sim()
    tibble(Product = c("red cells (units)", "FFP (units)",
                       "fibrinogen concentrate (g)", "platelets (units)",
                       "calcium (mmol)", "crystalloid still intravascular (mL)"),
           Given = c(round(max(d$RBCU), 1), round(max(d$FFPU), 1),
                     round(max(d$FGCU), 1), round(max(d$PLTU), 1),
                     round(max(d$CAU), 1), round(max(d$CIV))))
  })

  ## ---- tab 6 ------------------------------------------------------
  output$p_end <- renderPlot({
    cols <- c(LOSS = "cumulative loss (mL)", ShockIx = "shock index",
              ODEBT = "accumulated oxygen debt (mL O2)",
              TSEV = "minutes with MAP < 50 (Sheehan risk)",
              TAKI = "minutes with MAP < 65 (AKI risk)",
              UrineRate = "urine output (mL/min)")
    ggplot(long(sim(), cols), aes(time, value)) +
      geom_line(linewidth = 0.9, colour = "#37474F") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "minutes", y = NULL, title = "Endpoints, most of them integrals of time") +
      THEME
  })
  output$t_end <- renderTable({
    endpoints(sim()) %>%
      mutate(across(everything(), as.character)) %>%
      pivot_longer(everything(), names_to = "Endpoint", values_to = "Value")
  })

  ## ---- tab 7 ------------------------------------------------------
  output$p_cmp <- renderPlot({
    s <- saved()
    if (!length(s)) {
      return(ggplot() + annotate("text", 0, 0, label =
        "Set up an arm on the left and press 'Save this arm for comparison'.") +
        theme_void())
    }
    d <- bind_rows(lapply(names(s), function(n)
      s[[n]] %>% select(time, LOSS, Qbleed, MAP, Fib) %>% mutate(arm = n)))
    d %>% pivot_longer(c(LOSS, Qbleed, MAP, Fib)) %>%
      mutate(name = recode(name, LOSS = "cumulative loss (mL)",
                           Qbleed = "leak (mL/min)", MAP = "MAP (mmHg)",
                           Fib = "fibrinogen (g/L)")) %>%
      ggplot(aes(time, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "minutes", y = NULL, colour = NULL) + THEME
  })
  output$t_cmp <- renderTable({
    s <- saved()
    if (!length(s)) return(NULL)
    bind_rows(lapply(names(s), function(n) endpoints(s[[n]]) %>% mutate(arm = n))) %>%
      select(arm, everything()) %>% mutate(across(everything(), as.character))
  })

  ## ---- tab 8 ------------------------------------------------------
  output$p_sweep <- renderPlot({
    b <- build()
    grid <- seq(0.5, 0.95, by = 0.025)
    res <- bind_rows(lapply(grid, function(a) {
      d <- zero_re(param(pph, modifyList(b$pars, list(ATN = a)))) %>%
        mrgsim(events = b$ev, end = 240, delta = 2) %>% as_tibble()
      dead <- d %>% filter(Lethal > 0) %>% slice(1)
      tibble(ATN = a, loss = max(d$LOSS),
             outcome = if (nrow(dead)) "death" else
                       if (tail(d$Qbleed, 1) < 25) "haemostasis" else "still bleeding")
    }))
    ggplot(res, aes(ATN, loss, colour = outcome)) +
      geom_line(colour = "grey60") + geom_point(size = 3) +
      scale_colour_manual(values = c(death = "#B71C1C", haemostasis = "#1B5E20",
                                     `still bleeding` = "#E65100")) +
      geom_hline(yintercept = c(1000, 1500), linetype = 2, colour = "grey40") +
      labs(x = "atony severity (ATN)", y = "cumulative blood loss (mL)",
           colour = NULL,
           title = "The knife edge: this treatment plan swept across severity") +
      THEME
  })
  output$note_sweep <- renderUI(HTML(
    "<p><em>Because patency is quartic in tone, outcome is close to a step
     function of contractile capacity. With no treatment the step sits between
     ATN 0.60 (2218 mL, survives) and ATN 0.63 (death at 70 min). Each therapy
     you add on the left moves the step to the right &mdash; that displacement,
     not the change in average blood loss, is what the therapy is worth.</em></p>"))

  ## ---- tab 9 ------------------------------------------------------
  ## only run the 300-subject simulation while the population tab is open
  pop <- reactive({
    req(input$tabs == "9 · Population")
    b <- build()
    set.seed(20260804)
    n <- 300
    idata <- tibble(ID = seq_len(n),
                    ATN = plogis(qlogis(0.55) + rnorm(n, 0, 0.9)),
                    OTR0 = plogis(qlogis(0.75) + rnorm(n, 0, 0.8)))
    zero_re(param(pph, b$pars)) %>%
      mrgsim(idata = idata, events = b$ev, end = 240, delta = 4) %>%
      as_tibble() %>% group_by(ID) %>%
      summarise(loss = max(LOSS), lethal = max(Lethal),
                atn = first(AtonySeverity), .groups = "drop")
  })

  output$p_pop <- renderPlot({
    d <- pop()
    ggplot(d, aes(loss)) +
      geom_histogram(bins = 40, fill = "#5E35B1", colour = "white") +
      geom_vline(xintercept = c(1000, 1500, 2500), linetype = 2, colour = "#B71C1C") +
      labs(x = "cumulative blood loss (mL)", y = "women",
           title = paste("300 women on the same plan: the distribution is",
                         "bimodal, not bell-shaped")) + THEME
  })
  output$t_pop <- renderTable({
    d <- pop()
    tibble(`median loss (mL)` = round(median(d$loss)),
           `PPH >=1000 mL` = sprintf("%.1f%%", 100 * mean(d$loss >= 1000)),
           `severe >=1500 mL` = sprintf("%.1f%%", 100 * mean(d$loss >= 1500)),
           `massive >=2500 mL` = sprintf("%.1f%%", 100 * mean(d$loss >= 2500)),
           `lethal` = sprintf("%.1f%%", 100 * mean(d$lethal > 0)))
  })
}

shinyApp(ui, server)
