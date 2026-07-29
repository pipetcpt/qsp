## ===========================================================================
##  ABPA QSP — Shiny dashboard
##  ---------------------------------------------------------------------------
##  10 tabs.  The layout is deliberately opinionated in three places, because
##  the model's three main results are all results about MEASUREMENT and a
##  dashboard that lets you plot the wrong thing would undo them:
##
##   1. Total IgE is NEVER plotted alone.  The IgE tab always shows total and
##      free on the same panel with a warning band, because they move in
##      opposite directions on omalizumab and the ABPA response criterion is
##      written on the one that goes the wrong way.
##   2. Fungal burden is NEVER shown without its PARTITION.  Every fungal panel
##      shows the sanctuary share alongside the total, because a falling total
##      with a rising sanctuary share is a different clinical situation from a
##      falling total with a falling share.
##   3. No steroid arm is shown without its currency.  The corticosteroid tab
##      carries cortisol, BMD, HbA1c and cumulative prednisolone-equivalent on
##      the same screen as the efficacy curves.
##
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  Run:      shiny::runApp("abpa_shiny_app.R")
## ===========================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

MOD <- mread_cache("abpa_mrgsolve_model", ".")

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom")

PAL <- c(untreated = "#7f7f7f", prednisolone = "#2e7d32", itraconazole = "#c62828",
         voriconazole = "#ad1457", `neb L-AmB` = "#6a1b9a", omalizumab = "#1565c0",
         mepolizumab = "#00838f", benralizumab = "#00695c", dupilumab = "#ef6c00",
         `dupi + itra` = "#4527a0", `bud + itra` = "#795548",
         combination = "#000000")

## ---------------------------------------------------------------------------
##  dosing helpers (mirror abpa_reference_implementation.py exactly)
## ---------------------------------------------------------------------------
isham_taper <- function(days, wt) {
  d0 <- 0.5 * wt; out <- NULL
  for (day in 0:(days - 1)) {
    dose <- if (day < 14) d0
            else if (day < 70) { if ((day - 14) %% 2 == 0) d0 else 0 }
            else if (day < 154) max(d0 - 5 * ((day - 70) %/% 14 + 1), 0)
            else 0
    if (dose > 0) out <- c(out, ev(time = day, amt = dose, cmt = "APRD"))
  }
  out
}

build_events <- function(inp) {
  d <- inp$days; e <- NULL
  if (inp$itra > 0)  e <- c(e, ev(time = 0, amt = inp$itra * 2, cmt = "AITR",
                                  ii = 1, addl = d - 1))
  if (inp$vori > 0)  e <- c(e, ev(time = 0, amt = inp$vori * 2, cmt = "AVOR",
                                  ii = 1, addl = d - 1))
  if (inp$amb > 0)   e <- c(e, ev(time = 0, amt = inp$amb, cmt = "AMB",
                                  ii = 2, addl = floor(d / 2) - 1))
  if (inp$oma > 0)   e <- c(e, ev(time = 0, amt = inp$oma, cmt = "OMAD",
                                  ii = inp$oma_ii, addl = floor(d / inp$oma_ii) - 1))
  if (inp$mepo > 0)  e <- c(e, ev(time = 0, amt = inp$mepo, cmt = "MEPD",
                                  ii = 28, addl = floor(d / 28) - 1))
  if (inp$benra > 0) e <- c(e, ev(time = 0, amt = inp$benra, cmt = "BEND",
                                  ii = 28, addl = 2),
                              ev(time = 84, amt = inp$benra, cmt = "BEND",
                                 ii = 56, addl = max(floor((d - 84) / 56) - 1, 0)))
  if (inp$dupi > 0)  e <- c(e, ev(time = 0, amt = inp$dupi, cmt = "DUPD",
                                  ii = 14, addl = floor(d / 14) - 1))
  if (inp$teze > 0)  e <- c(e, ev(time = 0, amt = inp$teze, cmt = "TEZD",
                                  ii = 28, addl = floor(d / 28) - 1))
  if (inp$bud > 0)   e <- c(e, ev(time = 0, amt = inp$bud, cmt = "ABUD",
                                  ii = 1, addl = d - 1))
  if (inp$steroid == "ISHAM taper")  e <- c(e, isham_taper(d, inp$wt))
  if (inp$steroid == "maintenance")  e <- c(e, ev(time = 0, amt = inp$pred_mg,
                                                  cmt = "APRD", ii = 1, addl = d - 1))
  if (inp$steroid == "methylpred")   e <- c(e, ev(time = 0, amt = inp$pred_mg * 0.8,
                                                  cmt = "AMPD", ii = 1, addl = d - 1))
  if (is.null(e)) e <- ev(time = 0, amt = 0, cmt = "AITR")
  e
}

run_arm <- function(inp, overrides = list()) {
  m <- MOD
  if (length(overrides)) m <- do.call(param, c(list(m), overrides))
  m <- param(m, f_pen = inp$f_pen, gp = inp$gp, cyp2c19 = inp$cyp2c19,
             kelCX = log(2) / inp$cplx_t12)
  init0 <- list()
  if (inp$severity > 0) init0 <- list(BRON = inp$severity, PLUG = 1 + 0.06 * inp$severity)
  if (length(init0)) m <- do.call(init, c(list(m), init0))
  if (inp$lavage) m <- init(m, FPLG = 0.02, PLUG = 0.05)
  m %>% ev(build_events(inp)) %>%
    mrgsim(end = inp$days, delta = 1, atol = 1e-8, rtol = 1e-6) %>% as_tibble()
}

PRESETS <- list(
  "untreated"        = list(itra = 0,   steroid = "none"),
  "prednisolone"     = list(itra = 0,   steroid = "ISHAM taper"),
  "itraconazole"     = list(itra = 200, steroid = "none"),
  "voriconazole"     = list(vori = 200, steroid = "none"),
  "neb L-AmB"        = list(amb = 10,   steroid = "none"),
  "omalizumab"       = list(oma = 375,  steroid = "none"),
  "mepolizumab"      = list(mepo = 100, steroid = "none"),
  "benralizumab"     = list(benra = 30, steroid = "none"),
  "dupilumab"        = list(dupi = 300, steroid = "none"),
  "dupi + itra"      = list(dupi = 300, itra = 200, steroid = "none"),
  "bud + itra"       = list(bud = 1.6,  itra = 200, steroid = "none")
)

BLANK <- list(days = 364, wt = 70, itra = 0, vori = 0, amb = 0, oma = 0,
              oma_ii = 14, mepo = 0, benra = 0, dupi = 0, teze = 0, bud = 0,
              pred_mg = 10, steroid = "none", f_pen = 0.10, gp = 0.35,
              cyp2c19 = 1.0, cplx_t12 = 8, severity = 0, lavage = FALSE)

## ===========================================================================
##  UI
## ===========================================================================
ui <- fluidPage(
  titlePanel("ABPA QSP — the sanctuary partition, the IgE trap and the CYP3A4 confound"),
  tags$head(tags$style(HTML("
    .warn {background:#fff3cd; border-left:5px solid #c9a227; padding:9px 12px;
           margin-bottom:10px; font-size:13px;}
    .bad  {background:#fdecea; border-left:5px solid #c0392b; padding:9px 12px;
           margin-bottom:10px; font-size:13px;}
    .ok   {background:#eaf4ea; border-left:5px solid #2e7d32; padding:9px 12px;
           margin-bottom:10px; font-size:13px;}
  "))),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("days", "Horizon (days)", 84, 1825, 364, step = 28),
      sliderInput("wt", "Weight (kg)", 40, 130, 70, step = 5),
      sliderInput("severity", "Baseline bronchiectasis score (0 = model default 3)",
                  0, 18, 0, step = 1),
      checkboxInput("lavage", "Bronchoscopic plug removal at day 0", FALSE),

      h4("Antifungal"),
      sliderInput("itra", "Itraconazole (mg BID)", 0, 400, 0, step = 50),
      sliderInput("vori", "Voriconazole (mg BID)", 0, 400, 0, step = 50),
      selectInput("cyp2c19", "CYP2C19 phenotype",
                  c("normal metaboliser" = 1.0, "poor metaboliser" = 0.40,
                    "ultrarapid" = 1.6), selected = 1.0),
      sliderInput("amb", "Nebulised L-AmB (mg q48h)", 0, 25, 0, step = 5),

      h4("Corticosteroid"),
      selectInput("steroid", "Systemic steroid",
                  c("none", "ISHAM taper", "maintenance", "methylpred")),
      sliderInput("pred_mg", "Maintenance dose (mg/d prednisolone-equiv)",
                  0, 40, 10, step = 2.5),
      sliderInput("bud", "Inhaled budesonide (mg/d)", 0, 3.2, 0, step = 0.4),

      h4("Biologic"),
      sliderInput("oma", "Omalizumab (mg per dose)", 0, 1500, 0, step = 75),
      selectInput("oma_ii", "Omalizumab interval (days)", c(14, 28), selected = 14),
      sliderInput("mepo", "Mepolizumab (mg q4w)", 0, 300, 0, step = 100),
      sliderInput("benra", "Benralizumab (mg)", 0, 60, 0, step = 30),
      sliderInput("dupi", "Dupilumab (mg q2w)", 0, 600, 0, step = 100),
      sliderInput("teze", "Tezepelumab (mg q4w)", 0, 420, 0, step = 70),

      h4("Structural parameters (the model's argument lives here)"),
      sliderInput("f_pen", "f_pen — azole penetration into a plug", 0.01, 1.0, 0.10,
                  step = 0.01),
      sliderInput("gp", "g_p — intra-plug fungal growth (NEVER MEASURED)",
                  0.05, 0.80, 0.35, step = 0.01),
      sliderInput("cplx_t12", "omalizumab:IgE complex half-life (d)", 3, 26, 8,
                  step = 1),
      hr(),
      checkboxGroupInput("compare", "Comparator arms",
                         names(PRESETS),
                         selected = c("untreated", "itraconazole", "dupi + itra"))
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("① Patient profile",
          div(class = "warn", HTML("<b>Read this before using the other tabs.</b>
            This model exists to make three measurement problems visible, so the
            dashboard refuses to display certain things in isolation: total IgE is
            always paired with free IgE, fungal burden is always paired with its
            sanctuary share, and no steroid arm is shown without its toxicity
            currency. Those pairings are the point, not decoration.")),
          fluidRow(column(6, plotOutput("p_profile", height = 340)),
                   column(6, plotOutput("p_endpoints", height = 340))),
          h4("State of the patient at the end of the horizon"),
          DTOutput("t_profile")),

        tabPanel("② Sanctuary partition",
          div(class = "bad", HTML("<b>The core of the model.</b> FLUM sees unbound
            plasma drug; FPLG sees <code>f_pen</code> × that. A falling total burden
            with a <i>rising</i> sanctuary share is not the same clinical situation as
            a falling total with a falling share — the second is a regimen that will
            hold, the first is one that will relapse.")),
          fluidRow(column(6, plotOutput("p_partition", height = 320)),
                   column(6, plotOutput("p_sanctfrac", height = 320))),
          fluidRow(column(12, plotOutput("p_threshold", height = 330))),
          div(class = "warn", HTML("The dashed line is the luminal kill rate the
            current azole regimen actually delivers. Where the solid E* curve rises
            above it, the regimen is <b>suppressive by construction</b> — no dose
            increase within the class ceiling changes that, only something that
            raises <code>k_out</code> does."))),

        tabPanel("③ Antifungal PK",
          div(class = "warn", HTML("Itraconazole clearance is
            <b>autoinhibiting</b>, so exposure is supra-proportional to dose;
            voriconazole clearance is saturable and CYP2C19-dependent. Both are
            reasons dose-proportional intuition fails here.")),
          fluidRow(column(6, plotOutput("p_azole_pk", height = 320)),
                   column(6, plotOutput("p_unbound", height = 320))),
          fluidRow(column(12, plotOutput("p_cyp", height = 300)))),

        tabPanel("④ IgE — the trap",
          div(class = "bad", HTML("<b>Total IgE and free IgE move in OPPOSITE
            directions on omalizumab.</b> Total rises roughly three-fold because the
            omalizumab:IgE complex is cleared more slowly than free IgE; free falls by
            ~96%. The ABPA response criterion (\"total IgE falls 35–50%\") is written
            on the curve that goes the wrong way, and will therefore score a working
            drug as a failure. Move the complex half-life slider: the size of the rise
            changes, the direction never does.")),
          fluidRow(column(6, plotOutput("p_ige", height = 330)),
                   column(6, plotOutput("p_fcer", height = 330))),
          fluidRow(column(12, DTOutput("t_ige_criterion")))),

        tabPanel("⑤ Omalizumab dosing",
          div(class = "warn", HTML("Omalizumab neutralises IgE stoichiometrically, so
            the requirement is a <b>molar flux</b> condition: drug delivered per day
            must exceed IgE produced per day. That is why the approved table is
            proportional to IgE × weight, and why it stops where the maximum dose
            stops matching the flux. ABPA's diagnostic threshold sits above that
            point.")),
          fluidRow(column(6, plotOutput("p_flux", height = 340)),
                   column(6, plotOutput("p_freeige_dose", height = 340))),
          DTOutput("t_dose")),

        tabPanel("⑥ Type-2 axis & mucus",
          fluidRow(column(6, plotOutput("p_cyto", height = 320)),
                   column(6, plotOutput("p_eos", height = 320))),
          fluidRow(column(6, plotOutput("p_plug", height = 320)),
                   column(6, plotOutput("p_kout", height = 320))),
          div(class = "warn", HTML("<code>k_out</code>, the plug clearance rate, is
            the hinge of the whole model: it is a <b>type-2</b> parameter, and it sets
            the <b>antifungal's</b> eradication threshold."))),

        tabPanel("⑦ Corticosteroid & its currency",
          div(class = "bad", HTML("Every panel on this tab is one screen on purpose.
            A steroid arm shown without cortisol, bone density, HbA1c and cumulative
            dose is a half-reported experiment.")),
          fluidRow(column(6, plotOutput("p_cs", height = 300)),
                   column(6, plotOutput("p_cort", height = 300))),
          fluidRow(column(4, plotOutput("p_bmd", height = 280)),
                   column(4, plotOutput("p_a1c", height = 280)),
                   column(4, plotOutput("p_cumo", height = 280)))),

        tabPanel("⑧ The CYP3A4 confound",
          div(class = "bad", HTML("Itraconazole inhibits CYP3A4 at the gut wall and in
            the liver, so adding it to a corticosteroid raises that steroid's exposure.
            Part of what reads as antifungal efficacy is therefore a drug interaction —
            and its size depends on <b>which steroid the protocol specified</b>:
            ~1.1× for prednisolone, ~2.6× for methylprednisolone, ~4× for inhaled
            budesonide.")),
          fluidRow(column(12, plotOutput("p_ddi", height = 340))),
          h4("Decomposition of the current arm"),
          DTOutput("t_ddi"),
          div(class = "warn", HTML("The <i>DDI off</i> row is a counterfactual, not a
            regimen: it is the same schedule with <code>DDI_OFF = 1</code>. The gap
            between it and the full arm is the part of the benefit that is
            pharmacokinetic."))),

        tabPanel("⑨ Scenario comparison",
          fluidRow(column(12, plotOutput("p_compare", height = 460))),
          h4("End-of-horizon comparison, every arm priced in its own currency"),
          DTOutput("t_compare")),

        tabPanel("⑩ Structure & irreversibility",
          div(class = "bad", HTML("Bronchiectasis is the only endpoint in this model
            that cannot be given back. It is an integrator: exacerbations and
            persistent plugs add to it and nothing subtracts. Every other curve on
            this dashboard can be recovered by a good enough regimen; this one
            cannot.")),
          fluidRow(column(6, plotOutput("p_bron", height = 330)),
                   column(6, plotOutput("p_fev1", height = 330))),
          fluidRow(column(12, plotOutput("p_hazard", height = 300))),
          div(class = "warn", HTML("Caveat the model cannot resolve: whether ABPA
            <i>causes</i> bronchiectasis or co-occurs with it is not settled in the
            literature (see §⑪ of the references). The integrator assumes causation,
            so every number on this tab inherits that assumption."))
        )
      )
    )
  )
)

## ===========================================================================
##  server
## ===========================================================================
server <- function(input, output, session) {

  inp <- reactive({
    modifyList(BLANK, list(
      days = input$days, wt = input$wt, itra = input$itra, vori = input$vori,
      amb = input$amb, oma = input$oma, oma_ii = as.numeric(input$oma_ii),
      mepo = input$mepo, benra = input$benra, dupi = input$dupi,
      teze = input$teze, bud = input$bud, pred_mg = input$pred_mg,
      steroid = input$steroid, f_pen = input$f_pen, gp = input$gp,
      cyp2c19 = as.numeric(input$cyp2c19), cplx_t12 = input$cplx_t12,
      severity = input$severity, lavage = input$lavage))
  })

  sim  <- reactive(run_arm(inp()))
  base <- reactive(run_arm(modifyList(inp(), list(
    itra = 0, vori = 0, amb = 0, oma = 0, mepo = 0, benra = 0, dupi = 0,
    teze = 0, bud = 0, steroid = "none", lavage = FALSE))))

  comps <- reactive({
    req(length(input$compare) > 0)
    bind_rows(lapply(input$compare, function(nm) {
      run_arm(modifyList(inp(), modifyList(
        list(itra = 0, vori = 0, amb = 0, oma = 0, mepo = 0, benra = 0,
             dupi = 0, teze = 0, bud = 0, steroid = "none"),
        PRESETS[[nm]]))) %>% mutate(arm = nm)
    }))
  })

  lp <- function(d, cols, ylab, title, sub = NULL) {
    d %>% select(time, all_of(cols)) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 0.85) +
      labs(x = "day", y = ylab, title = title, subtitle = sub, colour = NULL) +
      THEME
  }

  ## ---- ① profile ----------------------------------------------------------
  output$p_profile <- renderPlot(
    lp(sim(), c("IGE_TOTAL", "EOSB"), "IU/mL  |  cells/µL",
       "Serology", "total IgE and blood eosinophils") +
      scale_y_log10())

  output$p_endpoints <- renderPlot(
    lp(sim(), c("FEV1", "PLUG_SCORE", "BRON"), "value",
       "Endpoints", "FEV1 % pred · mucus-plug score /18 · bronchiectasis /18"))

  output$t_profile <- renderDT({
    s <- tail(sim(), 1)
    tibble(quantity = c("total IgE (IU/mL)", "FREE IgE (IU/mL)",
                        "free IgE (ng/mL)", "blood eosinophils (/µL)",
                        "mucus-plug score /18", "total fungal burden",
                        "sanctuary share", "FEV1 (% pred)",
                        "bronchiectasis /18", "expected exacerbations",
                        "cortisol (µg/dL)", "BMD (fraction)",
                        "cumulative prednisolone-equiv (mg)"),
           value = signif(c(s$IGE_TOTAL, s$IGE_FREE, s$IGE_FREE_NG, s$EOSB,
                            s$PLUG_SCORE, s$FTOT, s$SANCT_FRAC, s$FEV1, s$BRON,
                            s$CHAZ, s$CORT, s$BMD, s$CUMO), 4)) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  ## ---- ② sanctuary --------------------------------------------------------
  output$p_partition <- renderPlot(
    lp(sim(), c("FLUM", "FPLG", "FTOT"), "relative burden",
       "The partition",
       "FLUM sees full unbound drug; FPLG sees f_pen × that"))

  output$p_sanctfrac <- renderPlot(
    ggplot(sim(), aes(time, SANCT_FRAC)) +
      geom_line(linewidth = 0.9, colour = "#a01a45") +
      geom_hline(yintercept = 0.5, linetype = 3) +
      labs(x = "day", y = "FPLG / (FLUM + FPLG)",
           title = "Sanctuary share of the surviving organisms",
           subtitle = "a falling total with a rising share is a regimen that will relapse") +
      ylim(0, 1) + THEME)

  output$p_threshold <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_line(aes(y = E_FLOOR, colour = "E* floor = (g_p − k_out)/f_pen"),
                linewidth = 1) +
      geom_line(aes(y = E_LUM_CAP, colour = "delivered luminal kill rate"),
                linewidth = 1, linetype = 2) +
      geom_hline(aes(yintercept = 1.55), linetype = 3) +
      annotate("text", x = max(d$time) * 0.02, y = 1.60, hjust = 0, size = 3.4,
               label = "azole class ceiling (Emax)") +
      labs(x = "day", y = "kill rate (1/d)", colour = NULL,
           title = "Does the regimen clear the sanctuary, or only suppress it?",
           subtitle = "where the floor sits above the delivered rate, no azole dose within the class ceiling suffices") +
      THEME
  })

  ## ---- ③ antifungal PK ---------------------------------------------------
  output$p_azole_pk <- renderPlot(
    lp(sim(), c("ITRA", "OHIT", "VORI", "AMB"), "mg/L",
       "Antifungal concentrations",
       "itraconazole clearance is autoinhibiting; voriconazole is saturable"))

  output$p_unbound <- renderPlot({
    p <- as.list(param(MOD))
    sim() %>% transmute(time,
                        `unbound itraconazole` = p$fu_itra * ITRA,
                        `unbound OH-itra` = p$fu_oh * OHIT,
                        `unbound voriconazole` = p$fu_vori * VORI) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.85) +
      geom_hline(yintercept = p$EC50_af, linetype = 3) +
      annotate("text", x = 0, y = p$EC50_af * 1.35, hjust = 0, size = 3.3,
               label = "EC50 for fungal killing") +
      labs(x = "day", y = "mg/L unbound", colour = NULL,
           title = "Unbound concentration is what the fungus sees",
           subtitle = "itraconazole is 99.8% protein-bound: total concentration is a poor guide") +
      THEME
  })

  output$p_cyp <- renderPlot(
    ggplot(sim(), aes(time, I3A4_CAP)) +
      geom_line(linewidth = 0.9, colour = "#1f5a27") +
      geom_hline(yintercept = 1, linetype = 3) +
      ylim(0, 1.05) +
      labs(x = "day", y = "fractional CYP3A4 activity remaining",
           title = "The drug interaction, as one number",
           subtitle = "any co-administered CYP3A4 substrate has its exposure divided by roughly this") +
      THEME)

  ## ---- ④ IgE trap --------------------------------------------------------
  output$p_ige <- renderPlot({
    d <- sim(); b <- tail(base(), 1)
    d %>% select(time, IGE_TOTAL, IGE_FREE) %>% pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0.65 * d$IGE_TOTAL[1], linetype = 2,
                 colour = "#c0392b") +
      annotate("text", x = 0, y = 0.65 * d$IGE_TOTAL[1] * 0.86, hjust = 0,
               size = 3.4, colour = "#c0392b",
               label = "ABPA response criterion: total IgE −35%") +
      scale_y_log10() +
      scale_colour_manual(values = c(IGE_TOTAL = "#c0392b", IGE_FREE = "#1565c0")) +
      labs(x = "day", y = "IU/mL (log)", colour = NULL,
           title = "The two IgE curves go opposite ways",
           subtitle = "the criterion is drawn on the red one") +
      THEME
  })

  output$p_fcer <- renderPlot(
    lp(sim(), c("FCER_OCC", "FCER", "EFFECTOR"), "relative",
       "FcεRI occupancy, density and effector activation",
       "occupancy barely moves; the benefit is in the slow fall of receptor density"))

  output$t_ige_criterion <- renderDT({
    d <- sim(); s <- tail(d, 1)
    tibble(readout = c("total IgE, day 0 (IU/mL)", "total IgE, end (IU/mL)",
                       "total IgE fold-change", "free IgE, day 0 (IU/mL)",
                       "free IgE, end (IU/mL)", "free IgE suppression (%)",
                       "free IgE, end (ng/mL)  [target < 25]",
                       "effector activation drop (%)",
                       "response by TOTAL-IgE criterion",
                       "response by FREE-IgE criterion"),
           value = c(signif(d$IGE_TOTAL[1], 4), signif(s$IGE_TOTAL, 4),
                     signif(s$IGE_TOTAL / d$IGE_TOTAL[1], 3),
                     signif(d$IGE_FREE[1], 4), signif(s$IGE_FREE, 4),
                     signif(100 * (1 - s$IGE_FREE / d$IGE_FREE[1]), 3),
                     signif(s$IGE_FREE_NG, 4),
                     signif(100 * (1 - s$EFFECTOR / d$EFFECTOR[1]), 3),
                     ifelse(s$IGE_TOTAL <= 0.65 * d$IGE_TOTAL[1],
                            "met", "NOT MET"),
                     ifelse(s$IGE_FREE <= 0.50 * d$IGE_FREE[1],
                            "responder", "non-responder"))) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 12))
  })

  ## ---- ⑤ omalizumab dosing ----------------------------------------------
  output$p_flux <- renderPlot({
    p <- as.list(param(MOD))
    OMA_MG_TO_NM <- 1.2293; IGE_IU_TO_NM <- 0.01263
    ige0 <- seq(250, 20000, by = 250)
    need <- p$kdegE * ige0 * IGE_IU_TO_NM / p$eta_oma        # nM/d of omalizumab
    supply <- function(mg, ii) mg * p$F_oma * OMA_MG_TO_NM / ii
    tibble(ige0,
           `IgE production (as omalizumab needed)` = need,
           `375 mg q2w delivers` = supply(375, 14),
           `current dose delivers` = supply(max(input$oma, 1e-9),
                                            as.numeric(input$oma_ii))) %>%
      pivot_longer(-ige0) %>%
      ggplot(aes(ige0, value, colour = name)) + geom_line(linewidth = 0.9) +
      geom_vline(xintercept = 1500, linetype = 2) +
      annotate("text", x = 1650, y = max(need) * 0.9, hjust = 0, size = 3.4,
               label = "label cap: total IgE 1500 IU/mL") +
      geom_vline(xintercept = 1000, linetype = 3, colour = "#c0392b") +
      annotate("text", x = 1050, y = max(need) * 0.6, hjust = 0, size = 3.4,
               colour = "#c0392b", label = "ABPA diagnostic threshold") +
      labs(x = "baseline total IgE (IU/mL)", y = "nM/day", colour = NULL,
           title = "The dosing table is a molar flux-matching rule",
           subtitle = "where the supply line falls below the demand line, free IgE cannot be suppressed") +
      THEME
  })

  output$p_freeige_dose <- renderPlot({
    doses <- c(75, 150, 225, 300, 375, 525, 750, 1050, 1500)
    d <- bind_rows(lapply(doses, function(mg) {
      tail(run_arm(modifyList(inp(), list(oma = mg))), 1) %>%
        transmute(dose = mg, free_ng = IGE_FREE_NG, total = IGE_TOTAL)
    }))
    ggplot(d, aes(dose, free_ng)) +
      geom_line(linewidth = 0.9, colour = "#1565c0") + geom_point() +
      geom_hline(yintercept = 25, linetype = 2, colour = "#c0392b") +
      annotate("text", x = min(doses), y = 30, hjust = 0, size = 3.4,
               colour = "#c0392b", label = "classical free-IgE target 25 ng/mL") +
      geom_vline(xintercept = 375, linetype = 3) +
      scale_y_log10() +
      labs(x = "omalizumab dose per administration (mg)",
           y = "free IgE at end of horizon (ng/mL, log)",
           title = "Dose–response for the species that matters",
           subtitle = "the dotted line is the maximum dose on the approved table") +
      THEME
  })

  output$t_dose <- renderDT({
    doses <- c(150, 225, 300, 375, 525, 750, 1050, 1500)
    bind_rows(lapply(doses, function(mg) {
      s <- tail(run_arm(modifyList(inp(), list(oma = mg))), 1)
      tibble(`dose (mg)` = mg,
             `total IgE (IU/mL)` = signif(s$IGE_TOTAL, 4),
             `free IgE (IU/mL)` = signif(s$IGE_FREE, 4),
             `free IgE (ng/mL)` = signif(s$IGE_FREE_NG, 4),
             `target met` = ifelse(s$IGE_FREE_NG < 25, "yes", "no"),
             `on approved table` = ifelse(mg <= 375, "yes", "NO"))
    })) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- ⑥ type-2 axis ----------------------------------------------------
  output$p_cyto <- renderPlot(
    lp(sim(), c("AG", "TH2", "IL13", "IL5"), "relative",
       "Antigen and the type-2 amplifier",
       "antigen is released by BOTH fungal compartments"))
  output$p_eos <- renderPlot(
    lp(sim(), c("EOSB", "EOSA", "EPX"), "cells/µL  |  relative",
       "Eosinophils and granule chemistry",
       "EPX cross-links mucin, which is how eosinophils slow plug clearance") +
      scale_y_log10())
  output$p_plug <- renderPlot(
    lp(sim(), c("PLUG", "PLUG_SCORE"), "relative  |  score /18",
       "The plug", "the compartment that hides the organism"))
  output$p_kout <- renderPlot(
    ggplot(sim(), aes(time, KOUT_PLUG)) +
      geom_line(linewidth = 0.9, colour = "#5b569b") +
      geom_hline(aes(yintercept = as.list(param(MOD))$gp), linetype = 2,
                 colour = "#a01a45") +
      annotate("text", x = 0, y = as.list(param(MOD))$gp * 1.06, hjust = 0,
               size = 3.4, colour = "#a01a45",
               label = "g_p — above this line the sanctuary self-sustains") +
      labs(x = "day", y = "k_out (1/d)",
           title = "Plug clearance — a type-2 parameter that sets an antifungal threshold",
           subtitle = "the single most consequential quantity in the model") +
      THEME)

  ## ---- ⑦ steroid currency -----------------------------------------------
  output$p_cs <- renderPlot(
    lp(sim(), c("PRED", "MPRD", "BUD"), "mg/L", "Corticosteroid concentrations",
       "note the effect of any CYP3A4 inhibitor on these curves"))
  output$p_cort <- renderPlot(
    ggplot(sim(), aes(time, CORT)) + geom_line(linewidth = 0.9, colour = "#8b5b5b") +
      geom_hline(yintercept = 5, linetype = 2, colour = "#c0392b") +
      annotate("text", x = 0, y = 5.6, hjust = 0, size = 3.4, colour = "#c0392b",
               label = "adrenal insufficiency threshold") +
      labs(x = "day", y = "morning cortisol (µg/dL)",
           title = "HPA axis", subtitle = "the currency the efficacy is paid in") +
      THEME)
  output$p_bmd <- renderPlot(
    lp(sim(), "BMD", "fraction of baseline", "Bone mineral density"))
  output$p_a1c <- renderPlot(lp(sim(), "HBA1C", "%", "HbA1c"))
  output$p_cumo <- renderPlot(
    lp(sim(), "CUMO", "mg", "Cumulative prednisolone-equivalent"))

  ## ---- ⑧ the DDI --------------------------------------------------------
  ddi_arms <- reactive({
    i <- inp()
    bind_rows(
      run_arm(i) %>% mutate(arm = "as specified"),
      run_arm(i, overrides = list(DDI_OFF = 1)) %>% mutate(arm = "CYP3A4 interaction OFF"),
      run_arm(i, overrides = list(AF_OFF = 1)) %>% mutate(arm = "antifungal OFF (DDI only)"))
  })

  output$p_ddi <- renderPlot(
    ddi_arms() %>% select(time, arm, FEV1, PLUG_SCORE, FTOT, AUCCS) %>%
      pivot_longer(c(FEV1, PLUG_SCORE, FTOT, AUCCS)) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.85) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL, colour = NULL,
           title = "Counterfactuals: which part of the benefit is pharmacokinetic?",
           subtitle = "AUCCS is the GR exposure actually delivered — the honest denominator") +
      THEME)

  output$t_ddi <- renderDT({
    d <- ddi_arms() %>% group_by(arm) %>% slice_tail(n = 1) %>% ungroup()
    ref <- d$FEV1[d$arm == "as specified"]
    d %>% transmute(arm,
                    `AUC of GR effect` = signif(AUCCS, 4),
                    `total fungus` = signif(FTOT, 4),
                    `plug score` = signif(PLUG_SCORE, 3),
                    `FEV1` = signif(FEV1, 4),
                    `FEV1 vs full arm` = signif(FEV1 - ref, 3)) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- ⑨ comparison ------------------------------------------------------
  output$p_compare <- renderPlot(
    comps() %>% select(time, arm, IGE_TOTAL, IGE_FREE, EOSB, PLUG_SCORE,
                       FTOT, SANCT_FRAC, FEV1, BRON, HAZ_YR) %>%
      pivot_longer(-c(time, arm)) %>%
      mutate(name = factor(name, levels = c("IGE_TOTAL", "IGE_FREE", "EOSB",
                                            "PLUG_SCORE", "FTOT", "SANCT_FRAC",
                                            "FEV1", "BRON", "HAZ_YR"))) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = PAL, na.value = "grey40") +
      labs(x = "day", y = NULL, colour = NULL,
           title = "Every arm, every currency",
           subtitle = "total IgE and free IgE are adjacent panels on purpose") +
      THEME)

  output$t_compare <- renderDT({
    comps() %>% group_by(arm) %>% slice_tail(n = 1) %>% ungroup() %>%
      transmute(arm,
                `total IgE` = round(IGE_TOTAL),
                `free IgE` = round(IGE_FREE, 1),
                `eos /µL` = round(EOSB),
                `plug /18` = round(PLUG_SCORE, 2),
                `fungus` = signif(FTOT, 3),
                `sanctuary %` = round(100 * SANCT_FRAC, 1),
                `FEV1` = round(FEV1, 1),
                `bronchiect.` = round(BRON, 3),
                `exac (exp.)` = round(CHAZ, 3),
                `OCS mg` = round(CUMO),
                `cortisol` = round(CORT, 1),
                `clears sanctuary?` = ifelse(SANCT_OK > 0.5, "yes", "no")) %>%
      datatable(rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })

  ## ---- ⑩ structure -------------------------------------------------------
  output$p_bron <- renderPlot(
    comps() %>%
      ggplot(aes(time, BRON, colour = arm)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL, na.value = "grey40") +
      labs(x = "day", y = "bronchiectasis score /18", colour = NULL,
           title = "The one curve that never comes down",
           subtitle = "an integrator: exacerbations and persistent plugs add, nothing subtracts") +
      THEME)
  output$p_fev1 <- renderPlot(
    comps() %>% ggplot(aes(time, FEV1, colour = arm)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL, na.value = "grey40") +
      labs(x = "day", y = "FEV1 (% predicted)", colour = NULL,
           title = "FEV1 = reversible (plug, AHR) + fixed (bronchiectasis)") +
      THEME)
  output$p_hazard <- renderPlot(
    comps() %>% ggplot(aes(time, HAZ_YR, colour = arm)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL, na.value = "grey40") +
      labs(x = "day", y = "exacerbations per year", colour = NULL,
           title = "Instantaneous exacerbation hazard",
           subtitle = "each event feeds the bronchiectasis integrator") +
      THEME)
}

shinyApp(ui, server)
