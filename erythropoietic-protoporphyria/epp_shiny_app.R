## =============================================================================
##  ERYTHROPOIETIC PROTOPORPHYRIA (EPP) / X-LINKED PROTOPORPHYRIA (XLP)
##  Interactive QSP dashboard  —  Shiny front end for epp_mrgsolve_model.R
## =============================================================================
##
##  Run with:
##      shiny::runApp("epp_shiny_app.R")
##  (the file must sit next to epp_mrgsolve_model.R, from which it sources the
##   compiled mrgsolve model and the scenario builders)
##
##  NINE TABS
##  ---------
##   1  Patient & genotype      the two lesions (FECH LoF / ALAS2 GoF), MC1R
##                              status, iron status, behaviour
##   2  Heme pathway            flux through the eight enzymes and where it
##                              overflows; the derived ~35% threshold
##   3  Drug PK                 four disposition topologies on one axis, and
##                              the afamelanotide PK/PD hysteresis
##   4  Photobiology            the Soret band, epidermal transmittance,
##                              and tolerance time as a ratio of products
##   5  Clinical endpoints      tolerance time, pain-free sun hours, pain NRS
##   6  Scenario comparison     all nine therapeutic scenarios side by side
##                              plus the additive-vs-multiplicative test
##   7  Hepatic bistability     continuation of the cholestasis feedback loop
##                              and the location of the two folds
##   8  Control loop            the retreat-delay sweep: prodrome vs reaction
##   9  Diagnostics             the Zn-PP signature that separates EPP, XLP,
##                              iron deficiency and lead
##
##  EDUCATIONAL / RESEARCH USE ONLY — not validated for clinical decisions.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

## ---- load the model and helpers --------------------------------------------
## `epp_mrgsolve_model.R` only auto-runs when sourced as a script
## (sys.nframe() == 0L), so sourcing it here just gives us `mod` and helpers.
source("epp_mrgsolve_model.R", local = FALSE)

PAL <- c("#b03030", "#2060a0", "#207040", "#8040a0",
         "#c07030", "#606060", "#20808a", "#a02060", "#808020")

theme_epp <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#f0f0f0"),
          legend.position  = "bottom")
}

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel(paste("Erythropoietic Protoporphyria / X-Linked Protoporphyria —",
                   "QSP dashboard")),
  tags$p(style = "color:#555;margin-top:-10px;",
         paste("Protoporphyrin IX is both the product of the pathway and the",
               "substrate of the broken enzyme. Every intervention therefore",
               "acts on a ratio, and the sign of its effect is a property of",
               "which arm is rate-limiting — not of the drug.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Genotype"),
      sliderInput("FRES", "Residual ferrochelatase activity (% of normal)",
                  min = 2, max = 100, value = 15, step = 1),
      helpText(paste("Derived threshold is ~34.5%: normal chelation capacity",
                     "is 2.9x the normal ALAS2 flux, so below 1/2.9 the",
                     "terminal step can no longer keep up.")),
      sliderInput("GOF", "ALAS2 gain-of-function (XLP; 1 = none)",
                  min = 1, max = 6, value = 1, step = 0.1),
      sliderInput("FMC1R", "MC1R functionality (0.3 = red-hair variant)",
                  min = 0.1, max = 1, value = 1, step = 0.05),

      h4("Iron status"),
      sliderInput("FEDOSE", "Iron supply relative to baseline",
                  min = -0.6, max = 2, value = 0, step = 0.05),
      helpText(paste("Iron enters at three points with opposite signs, so the",
                     "net effect is computed: strongly helpful in XLP,",
                     "essentially inert in EPP.")),

      h4("Therapy"),
      checkboxInput("useAfa", "Afamelanotide 16 mg implant q60d", FALSE),
      checkboxInput("useDer", "Dersimelagon orally once daily", FALSE),
      conditionalPanel("input.useDer",
        sliderInput("derDose", "Dersimelagon dose (mg)", 25, 600, 300, 25)),
      checkboxInput("useBit", "Bitopertin orally once daily", FALSE),
      conditionalPanel("input.useBit",
        sliderInput("bitDose", "Bitopertin dose (mg)", 5, 120, 60, 5)),
      checkboxInput("useBC", "Beta-carotene 180 mg daily", FALSE),
      checkboxInput("useChol", "Cholestyramine 16 g daily", FALSE),
      checkboxInput("useTx", "Regular RBC transfusion", FALSE),

      h4("Environment & behaviour"),
      sliderInput("FCLOUD", "Cloud transmission at 400-410 nm",
                  0.2, 1, 1, 0.05),
      sliderInput("FCLOTH", "Clothing / opaque screen transmission",
                  0.05, 1, 1, 0.05),
      checkboxInput("glass", "Indoors behind window glass (transmits 405 nm)",
                    FALSE),
      sliderInput("TAUDELmin", "Retreat delay after the prodrome (min)",
                  0.5, 90, 2, 0.5),
      sliderInput("days", "Simulation length (days)", 30, 365, 180, 5),

      hr(),
      actionButton("run", "Run simulation", class = "btn-primary"),
      tags$p(style = "font-size:11px;color:#777;margin-top:10px;",
             paste("Educational / research model. Not validated for clinical",
                   "decision-making."))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## -------------------------------------------------- 1
        tabPanel(
          "1. Patient & genotype",
          br(),
          fluidRow(
            column(6, h4("Predicted steady-state phenotype"),
                   tableOutput("tbl_profile")),
            column(6, h4("Where this patient sits on the FECH curve"),
                   plotOutput("plt_fresCurve", height = "320px"))),
          hr(),
          htmlOutput("txt_profile")
        ),

        ## -------------------------------------------------- 2
        tabPanel(
          "2. Heme pathway",
          br(),
          plotOutput("plt_pathway", height = "340px"),
          hr(),
          h4("Flux balance at the terminal step"),
          tableOutput("tbl_flux"),
          helpText(paste("Upstream flux is INTACT in this disease. What fails",
                         "is disposal, so the pathway keeps delivering",
                         "substrate to a closed door and the excess leaves the",
                         "erythron as overflow."))
        ),

        ## -------------------------------------------------- 3
        tabPanel(
          "3. Drug PK / PD",
          br(),
          plotOutput("plt_pk", height = "300px"),
          hr(),
          h4("Afamelanotide: a five-day drug with a sixty-day effect"),
          plotOutput("plt_hysteresis", height = "300px"),
          helpText(paste("Plasma half-life ~30 min and the implant is spent by",
                         "about day 5, yet protection goes on rising for",
                         "another two to three weeks and peaks around day 20.",
                         "The dose interval is set by melanin turnover (~35 d),",
                         "not by the pharmacokinetics."))
        ),

        ## -------------------------------------------------- 4
        tabPanel(
          "4. Photobiology",
          br(),
          fluidRow(
            column(6, h4("Action spectrum and what blocks it"),
                   plotOutput("plt_spectrum", height = "300px")),
            column(6, h4("Tolerance time is a ratio of products"),
                   plotOutput("plt_ttolSurface", height = "300px"))),
          hr(),
          h4("Environmental modifiers"),
          tableOutput("tbl_env"),
          helpText(paste("SPF sunscreen appears at transmittance 1.00 on",
                         "purpose: ultraviolet filters do not attenuate the",
                         "400-410 nm Soret band that drives this disease."))
        ),

        ## -------------------------------------------------- 5
        tabPanel(
          "5. Clinical endpoints",
          br(),
          plotOutput("plt_endpoints", height = "380px"),
          hr(),
          fluidRow(
            column(6, h4("Season summary"), tableOutput("tbl_endpoints")),
            column(6, h4("Pain course after one bad exposure"),
                   plotOutput("plt_painCourse", height = "260px")))
        ),

        ## -------------------------------------------------- 6
        tabPanel(
          "6. Scenario comparison",
          br(),
          plotOutput("plt_scenarios", height = "400px"),
          hr(),
          h4("Nine scenarios, 180-day season means"),
          tableOutput("tbl_scenarios"),
          hr(),
          h4("Do the two axes add, or multiply?"),
          tableOutput("tbl_bliss"),
          helpText(paste("Photodynamic dose is a PRODUCT of irradiance,",
                         "epidermal transmittance and photosensitiser",
                         "concentration. Shielding and source-reduction act on",
                         "different factors of that product, so they compose",
                         "multiplicatively — supra-additively on any linear",
                         "scale."))
        ),

        ## -------------------------------------------------- 7
        tabPanel(
          "7. Hepatic bistability",
          br(),
          plotOutput("plt_continuation", height = "380px"),
          hr(),
          fluidRow(
            column(6, h4("Steady states at this genotype"),
                   tableOutput("tbl_branches")),
            column(6, h4("Rescue from the cholestatic branch"),
                   tableOutput("tbl_rescue"))),
          helpText(paste("Liver transplantation is deliberately absent from the",
                         "rescue table: it changes none of these parameters,",
                         "because the source is the marrow. That is exactly",
                         "why protoporphyric disease recurs in the graft and",
                         "why marrow transplantation is the curative option."))
        ),

        ## -------------------------------------------------- 8
        tabPanel(
          "8. The control loop",
          br(),
          plotOutput("plt_delay", height = "340px"),
          hr(),
          tableOutput("tbl_delay"),
          helpText(paste("The prodrome is the sensor, walking into shade is the",
                         "actuator, and human reaction time is the loop delay.",
                         "A shielding drug slows the rate at which dose",
                         "accrues, which makes the same delay non-limiting —",
                         "a benefit invisible to an 'hours in sunlight'",
                         "endpoint."))
        ),

        ## -------------------------------------------------- 9
        tabPanel(
          "9. Diagnostics",
          br(),
          h4("The zinc-protoporphyrin signature"),
          plotOutput("plt_diag", height = "340px"),
          hr(),
          tableOutput("tbl_diag"),
          helpText(paste("Ferrochelatase is also the enzyme that inserts zinc.",
                         "EPP breaks it, so the substrate piles up and the",
                         "metal cannot go in either. XLP leaves it intact and",
                         "merely floods it. Iron deficiency and lead raise",
                         "Zn-PP with a normal metal-free PPIX. Three different",
                         "positions in one plane, from one set of equations."))
        )
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  ## ---- assemble parameters from the sidebar --------------------------------
  parset <- reactive({
    list(FRES   = input$FRES / 100,
         GOF    = input$GOF,
         FMC1R  = input$FMC1R,
         FEDOSE = input$FEDOSE,
         FCLOUD = input$FCLOUD,
         FCLOTH = input$FCLOTH,
         FGLASS = if (isTRUE(input$glass)) 0.90 else 1.00,
         TAUDEL = input$TAUDELmin / 60,
         CHOLEFF= if (isTRUE(input$useChol)) 0.85 else 0.0,
         ERYSUP = if (isTRUE(input$useTx))   0.40 else 1.0)
  })

  events <- reactive({
    d <- input$days
    e <- NULL
    add <- function(a, b) if (is.null(a)) b else c(a, b)
    if (isTRUE(input$useAfa))
      e <- add(e, ev_afamelanotide(n = max(1, floor(d / 60)), 60, 16))
    if (isTRUE(input$useDer)) e <- add(e, ev_dersimelagon(d, input$derDose))
    if (isTRUE(input$useBit)) e <- add(e, ev_bitopertin(d, input$bitDose))
    if (isTRUE(input$useBC))  e <- add(e, ev_betacarotene(d, 180))
    e
  })

  ## ---- steady state, then the simulation -----------------------------------
  ss <- eventReactive(input$run, {
    withProgress(message = "Equilibrating the slow pools...", value = 0.3,
                 epp_steady(mod, parset()))
  }, ignoreNULL = FALSE)

  sim <- eventReactive(input$run, {
    withProgress(message = "Simulating...", value = 0.6, {
      m <- init_from(mod, ss()) %>% param(parset()) %>% param(EXPOSE = 1)
      e <- events()
      out <- if (is.null(e))
        m %>% mrgsim(end = input$days * 24, delta = 6, hmax = 0.05)
      else
        m %>% mrgsim(events = e, end = input$days * 24, delta = 6, hmax = 0.05)
      as.data.frame(out)
    })
  }, ignoreNULL = FALSE)

  ## =========================================================== TAB 1
  output$tbl_profile <- renderTable({
    s <- ss()
    data.frame(
      Quantity = c("Erythrocyte metal-free PPIX (umol/L)",
                   "Erythrocyte Zn-protoporphyrin (umol/L)",
                   "Total erythrocyte porphyrin (umol/L)",
                   "Metal-free fraction (%)",
                   "Plasma PPIX (umol/L)",
                   "Dermal (photoactive) PPIX (umol/L)",
                   "Hepatic PPIX (umol/L)",
                   "Cholestasis index (0-1)",
                   "ALT (U/L)", "Total bilirubin (mg/dL)",
                   "Tolerance time to prodrome (min)",
                   "Pain-free sun (h/day, censored)"),
      Value = c(s$PRBC, s$ZRBC, s$TOTEP, s$MFFRAC, s$PPL, s$PSK, s$PLIV,
                s$CHOL, s$ALT, s$TBIL, s$TTOLMIN, s$SUNHRD))
  }, digits = 3)

  output$plt_fresCurve <- renderPlot({
    grid <- seq(0.04, 1.0, by = 0.02)
    d <- do.call(rbind, lapply(grid, function(f) {
      x <- epp_steady(mod, modifyList(parset(), list(FRES = f)), days = 900)
      data.frame(FRES = 100 * f, PPIX = x$PRBC, TTOL = x$TTOLMIN)
    }))
    ggplot(d, aes(FRES, PPIX)) +
      geom_line(linewidth = 1, colour = PAL[1]) +
      geom_vline(xintercept = 34.5, linetype = "dashed", colour = PAL[6]) +
      annotate("text", x = 36, y = max(d$PPIX) * 0.9, hjust = 0, size = 3.4,
               label = "derived threshold 34.5%\n(= 1 / capacity ratio 2.9)") +
      geom_point(data = data.frame(FRES = input$FRES,
                                   PPIX = ss()$PRBC),
                 size = 3, colour = PAL[2]) +
      scale_x_reverse() +
      labs(x = "residual ferrochelatase activity (%)  [decreasing ->]",
           y = "erythrocyte PPIX (umol/L)") +
      theme_epp()
  })

  output$txt_profile <- renderUI({
    s <- ss()
    HTML(sprintf(paste0(
      "<p>At <b>%.0f%%</b> residual ferrochelatase activity the model predicts ",
      "an erythrocyte protoporphyrin of <b>%.1f umol/L</b>, of which ",
      "<b>%.0f%%</b> is metal-free, and a tolerance time of <b>%.0f minutes</b> ",
      "of direct summer sun before the prodrome.</p>",
      "<p>The threshold marked on the curve is not an input. Normal ",
      "chelation capacity in this model is 2.9 times the normal ALAS2 flux, so ",
      "the terminal step stops keeping up below 1/2.9 = 34.5%% residual ",
      "activity. That is why a null FECH allele is silent in a parent and ",
      "causal in a child who also inherited the common hypomorphic ",
      "IVS3-48T&gt;C allele in trans: the disease needs the ratio, not the ",
      "allele.</p>"),
      input$FRES, s$PRBC, s$MFFRAC, s$TTOLMIN))
  })

  ## =========================================================== TAB 2
  output$plt_pathway <- renderPlot({
    d <- sim()
    d %>% select(time, ALA, PBG, UPG, CPG, PPG, PPIXE) %>%
      pivot_longer(-time) %>%
      mutate(name = factor(name,
                           levels = c("ALA","PBG","UPG","CPG","PPG","PPIXE"),
                           labels = c("ALA","porphobilinogen",
                                      "uroporphyrinogen III",
                                      "coproporphyrinogen III",
                                      "protoporphyrinogen IX",
                                      "PROTOPORPHYRIN IX"))) %>%
      ggplot(aes(time / 24, value, colour = name)) +
      geom_line(linewidth = 0.8) +
      scale_colour_manual(values = PAL) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = "umol/L (erythron)", colour = NULL,
           title = "Intermediates are unremarkable; the terminal pool is not") +
      theme_epp() + theme(legend.position = "none")
  })

  output$tbl_flux <- renderTable({
    s <- ss(); p <- parset()
    fe <- s$FE
    fisc  <- (fe / (20 + fe)) / (20 / (20 + 20))
    ffech <- p$FRES * fisc
    cap   <- 435 * ffech * (fe / (10 + fe))
    sat   <- s$PPIXE / (1 + s$PPIXE)
    vfe   <- cap * sat
    jesc  <- 0.672 * s$PPIXE
    data.frame(
      Quantity = c("ALAS2 flux into the pathway (umol/L/h)",
                   "Maximum chelation CAPACITY (umol/L/h)",
                   "capacity / flux ratio",
                   "actual heme synthesis (umol/L/h)",
                   "OVERFLOW out of the erythron (umol/L/h)",
                   "fraction of flux that overflows (%)"),
      Value = c(vfe + jesc, cap, cap / (vfe + jesc), vfe, jesc,
                100 * jesc / (vfe + jesc)))
  }, digits = 3)

  ## =========================================================== TAB 3
  output$plt_pk <- renderPlot({
    d <- sim()
    d %>% select(time,
                 `afamelanotide (ug/L)` = AC,
                 `dersimelagon (ug/L)`  = DC,
                 `bitopertin (ug/L)`    = BC,
                 `beta-carotene (ug/L)` = CC) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 24, value, colour = name)) +
      geom_line(linewidth = 0.8) +
      scale_colour_manual(values = PAL) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = "plasma concentration", colour = NULL,
           title = "Four drugs, four disposition topologies") +
      theme_epp() + theme(legend.position = "none")
  })

  output$plt_hysteresis <- renderPlot({
    o <- init_from(mod, ss()) %>% param(parset()) %>% param(EXPOSE = 0) %>%
      mrgsim(events = ev_afamelanotide(3, 60, 16), end = 240 * 24, delta = 2) %>%
      as.data.frame()
    o %>% select(time, `plasma (ug/L)` = AC, `melanin OD` = MELIDX,
                 `tolerance (min)` = TTOLMIN) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 24, value)) +
      geom_line(linewidth = 0.8, colour = PAL[4]) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "day", y = NULL) + theme_epp()
  })

  ## =========================================================== TAB 4
  output$plt_spectrum <- renderPlot({
    lam <- seq(290, 700, by = 2)
    soret <- exp(-((lam - 405)^2) / (2 * 9^2))
    q1 <- 0.06 * exp(-((lam - 506)^2) / (2 * 10^2))
    q2 <- 0.04 * exp(-((lam - 541)^2) / (2 * 10^2))
    q3 <- 0.03 * exp(-((lam - 576)^2) / (2 * 10^2))
    q4 <- 0.05 * exp(-((lam - 630)^2) / (2 * 12^2))
    d <- data.frame(lam = lam, PPIX = soret + q1 + q2 + q3 + q4,
                    glass = as.numeric(lam > 340) * 0.9,
                    spf   = 1 - 0.98 * as.numeric(lam < 400),
                    melanin = exp(-(700 - lam) / 260))
    d %>% pivot_longer(-lam) %>%
      mutate(name = recode(name,
                           PPIX = "PPIX absorption (Soret + Q bands)",
                           glass = "window-glass transmission",
                           spf = "SPF 50 sunscreen transmission",
                           melanin = "eumelanin transmission")) %>%
      ggplot(aes(lam, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      annotate("rect", xmin = 400, xmax = 410, ymin = 0, ymax = 1.05,
               alpha = 0.15, fill = PAL[1]) +
      annotate("text", x = 405, y = 1.10, label = "400-410 nm", size = 3.2) +
      scale_colour_manual(values = PAL) +
      labs(x = "wavelength (nm)", y = "relative", colour = NULL,
           title = "Sunscreen stops at 400 nm; the disease starts at 400 nm") +
      theme_epp() + theme(legend.position = "bottom",
                          legend.text = element_text(size = 8))
  })

  output$plt_ttolSurface <- renderPlot({
    g <- expand.grid(mel = seq(1, 3, by = 0.05),
                     ppix = seq(0.2, 2.5, by = 0.05))
    g$ttol <- 60 / (8.20 * exp(-0.693 * g$mel) * g$ppix)
    s <- ss()
    ggplot(g, aes(ppix, mel, z = ttol)) +
      geom_raster(aes(fill = pmin(ttol, 120))) +
      geom_contour(breaks = c(10, 20, 30, 60, 120), colour = "white") +
      geom_point(data = data.frame(ppix = s$PSK, mel = s$MELIDX),
                 aes(ppix, mel), inherit.aes = FALSE, size = 3,
                 colour = "white") +
      scale_fill_viridis_c(name = "tolerance\n(min)") +
      labs(x = "dermal PPIX (umol/L)  [the source axis]",
           y = "melanin OD  [the shielding axis]",
           title = "Contours are hyperbolae, not lines") +
      theme_epp()
  })

  output$tbl_env <- renderTable({
    s <- ss()
    conds <- data.frame(
      condition = c("open summer sun", "behind ordinary window glass",
                    "overcast sky", "SPF 50 chemical sunscreen",
                    "iron-oxide tinted sunscreen",
                    "long sleeves + wide-brim hat", "full opaque cover"),
      transmitted = c(1.00, 0.90, 0.50, 1.00, 0.35, 0.15, 0.05))
    conds$tolerance_min <- s$TTOLMIN / conds$transmitted
    conds$sun_h_day <- 8 * (1 - exp(-3 * (conds$tolerance_min / 60) / 8))
    conds
  }, digits = 2)

  ## =========================================================== TAB 5
  output$plt_endpoints <- renderPlot({
    sim() %>%
      select(time, `erythrocyte PPIX (umol/L)` = PRBC,
             `melanin OD` = MELIDX,
             `tolerance time (min)` = TTOLMIN,
             `pain-free sun (h/day)` = SUNHRD,
             `pain NRS` = PAIN,
             `cumulative sun (h)` = SUNCUM) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time / 24, value)) +
      geom_line(linewidth = 0.8, colour = PAL[1]) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL) + theme_epp()
  })

  output$tbl_endpoints <- renderTable({
    d <- sim()
    data.frame(
      Endpoint = c("mean tolerance time (min)",
                   "mean pain-free sun (h/day)",
                   paste0("season pain-free sun (h over ", input$days, " d)"),
                   "final erythrocyte PPIX (umol/L)",
                   "mean melanin OD",
                   "peak pain NRS",
                   "phototoxic reaction-days"),
      Value = c(mean(d$TTOLMIN), mean(d$SUNHRD), tail(d$SUNCUM, 1),
                tail(d$PRBC, 1), mean(d$MELIDX), max(d$PAIN),
                tail(d$RXNCUM, 1)))
  }, digits = 2)

  output$plt_painCourse <- renderPlot({
    m <- init_from(mod, ss()) %>% param(parset())
    one <- m %>% param(EXPOSE = 1) %>%
      mrgsim(end = 24, delta = 0.25, hmax = 0.01) %>% as.data.frame()
    st <- tail(one, 1)
    rec <- init_from(mod, st) %>% param(parset()) %>% param(EXPOSE = 0) %>%
      mrgsim(end = 96, delta = 1, hmax = 0.05) %>% as.data.frame()
    rec$hours <- rec$time
    ggplot(rec, aes(hours, PAIN)) +
      geom_line(linewidth = 1, colour = PAL[1]) +
      labs(x = "hours after the exposure", y = "pain NRS 0-10",
           title = "A three-microsecond species, a three-day illness",
           subtitle = paste("singlet oxygen lives ~3 us; the persistence is in",
                            "the slow inflammatory states it switched on")) +
      ylim(0, 10) + theme_epp()
  })

  ## =========================================================== TAB 6
  scen <- eventReactive(input$run, {
    withProgress(message = "Running nine scenarios...", value = 0.5, {
      bind_rows(
        scenario_untreated(mod),
        scenario_afamelanotide(mod),
        scenario_dersimelagon(mod, 300),
        scenario_bitopertin(mod, 60),
        scenario_combination(mod),
        scenario_xlp(mod))
    })
  }, ignoreNULL = FALSE)

  output$plt_scenarios <- renderPlot({
    scen() %>%
      select(time, scenario, `erythrocyte PPIX (umol/L)` = PRBC,
             `Zn-PP (umol/L)` = ZRBC, `melanin OD` = MELIDX,
             `tolerance time (min)` = TTOLMIN) %>%
      pivot_longer(-c(time, scenario)) %>%
      ggplot(aes(time / 24, value, colour = scenario)) +
      geom_line(linewidth = 0.8) +
      scale_colour_manual(values = PAL) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL, colour = NULL,
           title = "The two axes are orthogonal in mechanism",
           subtitle = paste("MC1R agonists move melanin and never PPIX;",
                            "bitopertin moves PPIX and never melanin")) +
      theme_epp()
  })

  output$tbl_scenarios <- renderTable({
    scen() %>% group_by(scenario) %>%
      summarise(`RBC PPIX` = last(PRBC), `Zn-PP` = last(ZRBC),
                `metal-free %` = last(MFFRAC), melanin = mean(MELIDX),
                `tolerance min` = mean(TTOLMIN),
                `sun h/day` = mean(SUNHRD),
                `season h` = mean(SUNHRD) * 180, .groups = "drop")
  }, digits = 2)

  output$tbl_bliss <- renderTable({ combination_index(mod) }, digits = 3)

  ## =========================================================== TAB 7
  output$plt_continuation <- renderPlot({
    d <- hepatic_continuation(mod, seq(0.30, 0.05, by = -0.01))
    d %>% select(FRES, `hepatic PPIX (umol/L)` = PLIV,
                 `cholestasis index` = CHOL, `ALT (U/L)` = ALT,
                 `bilirubin (mg/dL)` = TBIL) %>%
      pivot_longer(-FRES) %>%
      ggplot(aes(100 * FRES, value)) +
      geom_line(linewidth = 1, colour = PAL[1]) +
      geom_vline(xintercept = c(16.1, 7.4), linetype = "dashed",
                 colour = PAL[6]) +
      facet_wrap(~name, scales = "free_y") +
      scale_x_reverse() +
      labs(x = "residual ferrochelatase activity (%)  [decreasing ->]",
           y = NULL,
           title = "Protoporphyric hepatopathy is a saddle-node, not an accumulation",
           subtitle = paste("dashed lines: the two folds. Between them the",
                            "patient is bistable — biochemically stable for",
                            "decades, one transient away from an irreversible",
                            "jump.")) +
      theme_epp()
  })

  output$tbl_branches <- renderTable({
    s <- ss()
    data.frame(Quantity = c("hepatic PPIX (umol/L)", "cholestasis index",
                            "bile flow (fraction of normal)",
                            "ALT (U/L)", "bilirubin (mg/dL)"),
               Value = c(s$PLIV, s$CHOL, s$BILEFLOW, s$ALT, s$TBIL))
  }, digits = 3)

  output$tbl_rescue <- renderTable({
    scenario_hepatopathy_rescue(mod) %>%
      select(arm, `hepatic PPIX` = PLIV, cholestasis = CHOL,
             ALT, bilirubin = TBIL, `RBC PPIX` = PRBC)
  }, digits = 2)

  ## =========================================================== TAB 8
  delaySweep <- eventReactive(input$run, {
    scenario_control_loop(mod, c(0.5, 1, 2, 5, 10, 20, 40, 90))
  }, ignoreNULL = FALSE)

  output$plt_delay <- renderPlot({
    delaySweep() %>%
      select(delay_min, untreated = peak_pain_untreated,
             afamelanotide = peak_pain_afa) %>%
      pivot_longer(-delay_min) %>%
      ggplot(aes(delay_min, value, colour = name)) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      scale_x_log10() +
      scale_colour_manual(values = PAL[1:2]) +
      labs(x = "retreat delay after the prodrome (min, log scale)",
           y = "peak pain NRS 0-10", colour = NULL,
           title = "The reaction is all-or-nothing, and shielding moves the step",
           subtitle = paste("same sun, same PPIX; only the loop delay differs.",
                            "Between 20 and 40 minutes the mast-cell and",
                            "complement arm becomes self-amplifying.")) +
      theme_epp()
  })

  output$tbl_delay <- renderTable({ delaySweep() }, digits = 2)

  ## =========================================================== TAB 9
  output$plt_diag <- renderPlot({
    cases <- list(
      "normal"                 = list(FRES = 1.00, GOF = 1.0, FEDOSE = 0),
      "silent FECH carrier"    = list(FRES = 0.50, GOF = 1.0, FEDOSE = 0),
      "EPP (FECH 25%)"         = list(FRES = 0.25, GOF = 1.0, FEDOSE = 0),
      "EPP (FECH 15%)"         = list(FRES = 0.15, GOF = 1.0, FEDOSE = 0),
      "EPP (FECH 8%)"          = list(FRES = 0.08, GOF = 1.0, FEDOSE = 0),
      "XLP (ALAS2 GoF)"        = list(FRES = 1.00, GOF = 3.5, FEDOSE = 0),
      "iron deficiency"        = list(FRES = 1.00, GOF = 1.0, FEDOSE = -0.5))
    d <- do.call(rbind, lapply(names(cases), function(n) {
      x <- epp_steady(mod, cases[[n]], days = 1200)
      data.frame(case = n, PPIX = x$PRBC, ZNPP = x$ZRBC,
                 MFFRAC = x$MFFRAC, TOTEP = x$TOTEP)
    }))
    ggplot(d, aes(PPIX, ZNPP, colour = case)) +
      geom_point(size = 4) +
      geom_text(aes(label = case), size = 3.4, vjust = -1.1,
                show.legend = FALSE) +
      scale_x_log10() + scale_y_log10() +
      scale_colour_manual(values = PAL) +
      labs(x = "metal-free erythrocyte PPIX (umol/L, log)",
           y = "zinc protoporphyrin (umol/L, log)",
           title = "One plane separates four different diseases",
           subtitle = paste("EPP moves right (substrate up, metal blocked);",
                            "XLP moves up and right; iron deficiency moves up",
                            "only.")) +
      theme_epp() + theme(legend.position = "none")
  })

  output$tbl_diag <- renderTable({
    cases <- list(
      "normal"              = list(FRES = 1.00, GOF = 1.0, FEDOSE = 0),
      "EPP (FECH 15%)"      = list(FRES = 0.15, GOF = 1.0, FEDOSE = 0),
      "XLP (ALAS2 GoF)"     = list(FRES = 1.00, GOF = 3.5, FEDOSE = 0),
      "iron deficiency"     = list(FRES = 1.00, GOF = 1.0, FEDOSE = -0.5))
    do.call(rbind, lapply(names(cases), function(n) {
      x <- epp_steady(mod, cases[[n]], days = 1200)
      data.frame(case = n,
                 `metal-free PPIX` = x$PRBC, `Zn-PP` = x$ZRBC,
                 `total EP` = x$TOTEP, `metal-free %` = x$MFFRAC,
                 `tolerance min` = x$TTOLMIN, check.names = FALSE)
    }))
  }, digits = 2)
}

shinyApp(ui, server)
