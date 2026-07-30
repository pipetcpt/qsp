## =============================================================================
##  tts_shiny_app.R
##  Takotsubo syndrome QSP model — interactive dashboard
##
##  The app is organised around the model's one structural claim: the apex is
##  not a weak segment, it is the segment with the most beta2 receptors, so it
##  crosses the Gs-to-Gi threshold first.  Tab 3 (the switch) is therefore the
##  centre of the app; tabs 4-9 are consequences of it, tab 10 is the routing
##  claim about inotropes, and tab 11 is the set of experiments that would
##  falsify the whole thing.
##
##  Every slider that matters is a PHYSICAL input (receptor density, beta2
##  fraction, innervation, wall thickness, surge amplitude, adrenaline fraction,
##  oestradiol, septal geometry) rather than an outcome.  You cannot set the
##  ejection fraction in this app; you can only set the things that produce it.
##
##  Run:
##    setwd("takotsubo-syndrome")
##    shiny::runApp("tts_shiny_app.R")
##
##  Requires: shiny, mrgsolve  (ggplot2 / DT used if available, base graphics
##  and plain tables otherwise, so the app runs on a minimal install)
## =============================================================================

library(shiny)
library(mrgsolve)

HAS_GG <- requireNamespace("ggplot2", quietly = TRUE)
HAS_DT <- requireNamespace("DT",      quietly = TRUE)
if (HAS_GG) library(ggplot2)

source("tts_mrgsolve_model.R")

## -----------------------------------------------------------------------------
##  Presentation helpers
## -----------------------------------------------------------------------------
PHENO_CHOICES <- setNames(names(phenotypes),
                          vapply(phenotypes, function(p) p$label, character(1)))

SEGCOL  <- c(Apex = "#c2185b", Mid = "#e07b39", Base = "#1d6fb8")
ARMCOL  <- c(untreated = "#8c8c8c", treated = "#1f6feb", comparator = "#e07b39")

THEME <- if (HAS_GG) {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey35", size = 10),
          legend.position = "bottom")
} else NULL

## Plot a set of columns against time, in whichever graphics system is present.
tsplot <- function(df, cols, labels = cols, cols_colour = NULL,
                   title = "", subtitle = "", ylab = "", xunit = c("d", "h"),
                   ymin = NA, hline = NULL) {
  xunit <- match.arg(xunit)
  x <- if (xunit == "d") df$time/24 else df$time
  xlab <- if (xunit == "d") "days from trigger" else "hours from trigger"
  if (HAS_GG) {
    long <- do.call(rbind, lapply(seq_along(cols), function(i) data.frame(
      x = x, y = df[[cols[i]]], series = labels[i])))
    long$series <- factor(long$series, levels = labels)
    p <- ggplot(long, aes(x, y, colour = series)) +
      geom_line(linewidth = 0.85) +
      labs(title = title, subtitle = subtitle, x = xlab, y = ylab,
           colour = NULL) + THEME
    if (!is.null(cols_colour))
      p <- p + scale_colour_manual(values = setNames(cols_colour, labels))
    if (!is.null(hline))
      p <- p + geom_hline(yintercept = hline, linetype = 2, colour = "grey45")
    if (!is.na(ymin)) p <- p + expand_limits(y = ymin)
    return(p)
  }
  ys <- unlist(df[cols]); rng <- range(c(ys, ymin, hline), na.rm = TRUE)
  plot(x, df[[cols[1]]], type = "l", ylim = rng, xlab = xlab, ylab = ylab,
       main = title, col = if (is.null(cols_colour)) 1 else cols_colour[1],
       lwd = 2)
  if (length(cols) > 1)
    for (i in 2:length(cols))
      lines(x, df[[cols[i]]],
            col = if (is.null(cols_colour)) i else cols_colour[i], lwd = 2)
  if (!is.null(hline)) abline(h = hline, lty = 2, col = "grey45")
  legend("topright", legend = labels, bty = "n", lwd = 2,
         col = if (is.null(cols_colour)) seq_along(cols) else cols_colour)
  invisible(NULL)
}

tbl_out <- function(df) if (HAS_DT) DT::datatable(
  df, rownames = FALSE, options = list(dom = "t", pageLength = 40,
                                       scrollX = TRUE)) else df

renderTbl <- function(expr) if (HAS_DT) DT::renderDataTable(expr) else
  renderTable(expr)
tblUI <- function(id) if (HAS_DT) DT::dataTableOutput(id) else tableOutput(id)

## -----------------------------------------------------------------------------
##  UI
## -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Takotsubo syndrome — QSP model (64 ODEs)"),
  tags$p(style = "color:#555;margin-top:-8px",
         paste("Apical ballooning as a THRESHOLD, not a location: one beta2-AR",
               "Gs-to-Gi switch function evaluated at three receptor densities.",
               "Educational / research use only.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("pheno", "Phenotype (trigger + susceptibility)",
                  choices = PHENO_CHOICES, selected = "emo"),
      tags$hr(),
      tags$b("Trigger"),
      sliderInput("amp", "Surge amplitude (nmol/L/h)",
                  40, 900, 250, step = 10),
      sliderInput("frace", "Fraction routed to ADRENALINE",
                  0, 1, 0.75, step = 0.01),
      sliderInput("tau", "Surge decay constant (h)", 0.5, 72, 1.6, step = 0.5),
      tags$hr(),
      tags$b("Susceptibility"),
      sliderInput("e2", "Oestradiol (pg/mL)", 5, 150, 15, step = 5),
      sliderInput("sept", "Septal geometry SEPT (sets whether LVOTO is possible)",
                  0, 1.5, 0.35, step = 0.05),
      tags$hr(),
      tags$b("The receptor field — the ONLY inter-segment difference"),
      sliderInput("rho_ap", "Apical beta-AR density (base = 1.0)",
                  0.6, 2.0, 1.40, step = 0.05),
      sliderInput("fb2_ap", "Apical beta2 fraction", 0.10, 0.70, 0.42,
                  step = 0.01),
      sliderInput("fb2_bs", "Basal beta2 fraction", 0.10, 0.70, 0.24,
                  step = 0.01),
      sliderInput("inn_ap", "Apical innervation (MIBG gradient)",
                  0.3, 1.5, 0.62, step = 0.02),
      tags$hr(),
      tags$b("Treatment"),
      selectInput("drug", "Acute arm",
                  choices = c("none", "dobutamine", "adrenaline",
                              "noradrenaline", "milrinone", "levosimendan",
                              "esmolol", "phenylephrine", "IABP", "Impella"),
                  selected = "none"),
      sliderInput("dose", "Acute arm intensity (fraction of the reference dose)",
                  0, 2, 1, step = 0.1),
      checkboxInput("acei", "Ramipril 5 mg daily from day 1", FALSE),
      checkboxInput("bb",   "Metoprolol 50 mg bd from day 1", FALSE),
      checkboxInput("doac", "Apixaban 5 mg bd from day 1", FALSE),
      checkboxInput("ptx",  "Pertussis toxin (abolish Gi — in silico only)",
                    FALSE),
      sliderInput("qtdrug", "QT-prolonging comedication burden", 0, 1, 0,
                  step = 0.1),
      tags$hr(),
      sliderInput("end", "Simulate to (days)", 7, 365, 90, step = 7),
      actionButton("go", "Run", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        tabPanel("1 · Patient",
          h4("The derived resting state"),
          p("Nothing on this tab is typed in. The three segments are given a",
            "receptor density, a beta2 fraction and an innervation weight; the",
            "resting Gs coupling, cAMP, PKA and contractility are then solved",
            "from those. Note how close the three resting Gs values are: the",
            "higher apical receptor density is offset by the lower apical",
            "innervation, so",
            tags$b("at rest the apex is not a weak segment."),
            "It only becomes one when a blood-borne beta2 agonist arrives —",
            "which is why this disease needs adrenaline rather than",
            "sympathetic nerve traffic."),
          tblUI("t_rest"),
          br(),
          h4("Trigger profile actually applied"),
          plotOutput("p_surge", height = 260)
        ),

        tabPanel("2 · Catecholamines",
          h4("Plasma and interstitial catecholamines"),
          p("Adrenaline is blood-borne and reaches all segments equally.",
            "Noradrenaline is released at nerve terminals, so it follows the",
            "innervation gradient — and it is 25-fold weaker at beta2.",
            "That asymmetry is the whole reason the two agonists produce",
            "different diseases."),
          fluidRow(column(6, plotOutput("p_cat", height = 300)),
                   column(6, plotOutput("p_occ", height = 300))),
          h4("Reflex limb"),
          p("Once forward flow falls, the baroreflex sustains interstitial",
            "noradrenaline for days after the primary surge has cleared. That",
            "is what keeps the BASE hyperkinetic, and therefore what keeps an",
            "LVOT gradient open."),
          plotOutput("p_refl", height = 260)
        ),

        tabPanel("3 · The switch",
          h4("Gs → Gi trafficking: one function, three densities"),
          p("The phosphorylated (Gi-coupled) beta2 fraction PHOS is computed by",
            "the same equation in all three segments. The regional difference",
            "comes from the beta2 DENSITY that multiplies it, and from the",
            "PKA feedback that the denser segment generates for itself."),
          fluidRow(column(6, plotOutput("p_phos", height = 300)),
                   column(6, plotOutput("p_gi",   height = 300))),
          fluidRow(column(6, plotOutput("p_camp", height = 300)),
                   column(6, plotOutput("p_qss",  height = 300))),
          h4("Where does the threshold sit?"),
          p("A steady adrenaline infusion is ramped and the concentration at",
            "which each segment's Gi coupling exceeds 20% of its beta2 pool is",
            "recorded. The ordering is an output, not an input."),
          actionButton("scan", "Run the threshold scan"),
          tblUI("t_scan")
        ),

        tabPanel("4 · Segmental mechanics",
          h4("Where the ballooning comes from"),
          p("Apical shortening, mid shortening and basal shortening, plus the",
            "ballooning index (basal minus apical). No node in the model is",
            "told to be akinetic."),
          fluidRow(column(7, plotOutput("p_sh",   height = 320)),
                   column(5, plotOutput("p_ball", height = 320))),
          h4("Contractility and its three penalties"),
          plotOutput("p_cont", height = 300),
          tblUI("t_mech")
        ),

        tabPanel("5 · LVEF & endpoints",
          fluidRow(column(7, plotOutput("p_ef",   height = 330)),
                   column(5, plotOutput("p_wmsi", height = 330))),
          h4("Summary of this arm"),
          tblUI("t_end")
        ),

        tabPanel("6 · Haemodynamics & LVOT",
          fluidRow(column(6, plotOutput("p_hemo", height = 300)),
                   column(6, plotOutput("p_grad", height = 300))),
          p(tags$b("Read the two panels together."),
            "A gradient is the product of a patient property (septal geometry,",
            "the SEPT slider) and a physiological state (basal",
            "hypercontractility). Neither alone produces obstruction, which is",
            "why the same surge obstructs one patient and not another."),
          plotOutput("p_cong", height = 300)
        ),

        tabPanel("7 · Biomarkers",
          fluidRow(column(6, plotOutput("p_tni", height = 300)),
                   column(6, plotOutput("p_bnp", height = 300))),
          h4("The discordance"),
          p("Troponin reports lost membrane integrity; NT-proBNP reports",
            "end-diastolic wall stress. In this model the primary lesion is a",
            "SIGNALLING state of surviving myocytes, so the first is small",
            "while the second is large. The ratio is the model's diagnostic",
            "prediction, not a fitted output — it falls out of the two",
            "different physical drivers."),
          plotOutput("p_ratio", height = 300),
          fluidRow(column(6, plotOutput("p_infl", height = 280)))
        ),

        tabPanel("8 · Electrophysiology",
          fluidRow(column(6, plotOutput("p_qtc",  height = 300)),
                   column(6, plotOutput("p_tdp",  height = 300))),
          p("QTc is a slow remodelling state driven by tissue oedema AND by the",
            "apex-to-base dispersion of contractile state. Because the driver",
            "enters through a 35-hour lag, the QTc peak arrives on day 2-3",
            "rather than on arrival — which is what is observed."),
          plotOutput("p_k", height = 260)
        ),

        tabPanel("9 · Thrombus",
          fluidRow(column(6, plotOutput("p_stasis", height = 300)),
                   column(6, plotOutput("p_thr",    height = 300))),
          p("Stasis is not an independent risk factor here. It is computed from",
            "the apical akinesia and the cavity size, so a patient whose apex",
            "recovers quickly is predicted to need less anticoagulation."),
          actionButton("thrgo", "Compare anticoagulation arms"),
          tblUI("t_thr")
        ),

        tabPanel("10 · Inotrope routing",
          h4("Matched contractile drive, not matched dose"),
          p("Dobutamine and levosimendan are titrated to add the SAME",
            "mass-weighted contractile drive at hour 8. Comparing them at",
            "matched dose compares potency; comparing them at matched drive",
            "isolates the route. One raises cAMP (and therefore the",
            "phosphorylation that created the disease); the other sensitises",
            "the myofilament and leaves cAMP alone."),
          actionButton("mdgo", "Run the matched-drive comparison"),
          tblUI("t_md"),
          br(),
          h4("Every cAMP-route agent, side by side"),
          actionButton("inogo", "Run the inotrope panel"),
          tblUI("t_ino")
        ),

        tabPanel("11 · Falsification",
          h4("Three experiments that would break the model"),
          tags$ol(
            tags$li(tags$b("A pure noradrenaline surge must NOT balloon the apex."),
                    " NE is ~25-fold weaker at beta2, so it should not recruit Gi",
                    " anywhere. If the model balloons on NE alone, the agonist-",
                    "identity claim is wrong."),
            tags$li(tags$b("Pertussis toxin must abolish it."),
                    " If blocking Gi leaves the apical akinesis intact, the",
                    " akinesis was not coming from the switch."),
            tags$li(tags$b("Milrinone must reproduce dobutamine's harm."),
                    " Milrinone raises cAMP DOWNSTREAM of the receptor. If the",
                    " harm were a receptor-occupancy artefact, milrinone would",
                    " be safe. The model says it is not.")),
          actionButton("fgo", "Run all three"),
          tblUI("t_fals")
        ),

        tabPanel("12 · Causal decomposition",
          h4("Four counterfactual integrators"),
          p("cAMP has a ~1-minute half-life, so its quasi-steady-state value is",
            "available in closed form; that makes exact counterfactuals cheap.",
            "Each line below is the ejection fraction the same patient would",
            "have if ONE mechanism were removed. The shares are therefore",
            "measured rather than asserted — and they are not additive,",
            "because each counterfactual removes one mechanism alone."),
          plotOutput("p_cf", height = 320),
          actionButton("hdgo", "Tabulate the shares across phenotypes"),
          tblUI("t_hd")
        ),

        tabPanel("13 · Population & recurrence",
          h4("A virtual cohort"),
          p("Receptor density, septal geometry, surge amplitude, adrenaline",
            "fraction and oestradiol are sampled; every clinical readout is",
            "then a consequence. The spread of LVOT gradients is produced by",
            "the septal-geometry distribution alone."),
          sliderInput("npop", "Cohort size", 25, 400, 120, step = 25),
          actionButton("popgo", "Run the cohort"),
          plotOutput("p_pop", height = 320),
          tblUI("t_pop"),
          br(),
          h4("One-year recurrence: ACE inhibition vs beta-blockade"),
          p("The recurrence hazard is driven by residual receptor",
            "phosphorylation and reduced by ACE inhibition. Beta-blockade",
            "enters a different term and barely moves it — which is the",
            "registry observation, arrived at by routing rather than by rule."),
          actionButton("recgo", "Run the one-year arms"),
          tblUI("t_rec")
        ),

        tabPanel("14 · Calibration",
          h4("Eight parameters, eight published anchors"),
          p("Nelder-Mead against the acute apical phenotype. The oedema and",
            "stunning kinetics and the oedema-to-QTc slope are held FIXED:",
            "left free they trade against each other without limit (their",
            "product is what the QTc anchor sees), which buys a better",
            "objective and an uninterpretable model."),
          actionButton("calgo", "Show residuals at the shipped parameter values"),
          verbatimTextOutput("t_cal"),
          br(),
          h4("What the model does NOT reproduce"),
          tags$ul(
            tags$li("An isolated mid-ventricular ring cannot be produced by any",
                    " monotonic apex>mid>base gradient. The model gives apical",
                    " or apical-plus-mid patterns."),
            tags$li("Reverse (basal) takotsubo is not derived. It is reproduced",
                    " only by INVERTING the receptor gradient, which the model",
                    " states as a prediction about reverse TTS rather than as a",
                    " result."),
            tags$li("The oedema state was calibrated to the QTc time course",
                    " (days 2-3), which is faster than CMR T2 normalisation",
                    " (weeks to months). The model therefore under-represents",
                    " the persistence of tissue oedema."))
        )
      )
    )
  )
)

## -----------------------------------------------------------------------------
##  SERVER
## -----------------------------------------------------------------------------
DRUG_REF <- list(
  dobutamine    = list(par = "RATE_DOB", ref = 105,  stop = 52),
  adrenaline    = list(par = "RATE_EPI", ref = 5,    stop = 28),
  noradrenaline = list(par = "RATE_NEP", ref = 8,    stop = 52),
  milrinone     = list(par = "RATE_MIL", ref = 3.0,  stop = 52),
  levosimendan  = list(par = "RATE_LEV", ref = 0.42, stop = 28),
  esmolol       = list(par = "RATE_ESM", ref = 420,  stop = 52),
  phenylephrine = list(par = "RATE_PHE", ref = 30,   stop = 52),
  IABP          = list(par = "IABP",     ref = 1,    stop = 72),
  Impella       = list(par = "IMPELLA",  ref = 1,    stop = 72)
)

server <- function(input, output, session) {

  ## Any phenotype selection resets the physical sliders to that phenotype
  observeEvent(input$pheno, {
    p <- phenotypes[[input$pheno]]
    updateSliderInput(session, "amp",   value = p$AMP_TOT)
    updateSliderInput(session, "frace", value = p$FRAC_E)
    updateSliderInput(session, "tau",   value = p$TAU_SUR)
    updateSliderInput(session, "e2",    value = p$E2)
    updateSliderInput(session, "sept",  value = p$SEPT)
    updateSliderInput(session, "rho_ap", value = if (is.null(p$RHO_AP)) 1.40 else p$RHO_AP)
    updateSliderInput(session, "fb2_ap", value = if (is.null(p$FB2_AP)) 0.42 else p$FB2_AP)
    updateSliderInput(session, "fb2_bs", value = if (is.null(p$FB2_BS)) 0.24 else p$FB2_BS)
    updateSliderInput(session, "inn_ap", value = if (is.null(p$INN_AP)) 0.62 else p$INN_AP)
  })

  user_par <- reactive(list(
    AMP_TOT = input$amp, FRAC_E = input$frace, TAU_SUR = input$tau,
    E2 = input$e2, SEPT = input$sept, RHO_AP = input$rho_ap,
    FB2_AP = input$fb2_ap, FB2_BS = input$fb2_bs, INN_AP = input$inn_ap,
    PTX = as.numeric(input$ptx), QTDRUG = input$qtdrug))

  ## Baseline is re-derived from the sliders, not edited
  base_mod <- reactive({
    p <- user_par(); p0 <- p; p0$AMP_TOT <- 0
    m <- param(mod, p0)
    o <- mrgsim_df(zero_re(m), end = 336, delta = 336, recsort = 3)
    st <- as.list(o[nrow(o), intersect(names(o), names(as.list(init(mod))))])
    list(mod = init(param(mod, p), st), init = st)
  })

  oral_ev <- reactive({
    evs <- list()
    if (input$acei) evs[[length(evs)+1]] <- oral_events$ramipril
    if (input$bb)   evs[[length(evs)+1]] <- oral_events$metoprolol
    if (input$doac) evs[[length(evs)+1]] <- oral_events$apixaban
    if (!length(evs)) return(NULL)
    Reduce(function(a, b) rbind(as.data.frame(a), as.data.frame(b)), evs)
  })

  acute_changes <- reactive({
    if (identical(input$drug, "none") || input$dose == 0) return(list())
    d <- DRUG_REF[[input$drug]]
    v <- d$ref * input$dose
    list(list(t = 4,      par = setNames(list(v), d$par)),
         list(t = d$stop, par = setNames(list(0), d$par)))
  })

  run <- eventReactive(input$go, {
    b   <- base_mod()
    end <- input$end * 24
    ev  <- oral_ev()
    ch  <- acute_changes()
    withProgress(message = "integrating", value = 0.4, {
      if (is.null(ev)) {
        run_arm(b, ch, end = end, delta = max(0.25, end/4000))
      } else if (!length(ch)) {
        run_oral(b, ev, end = end, delta = max(0.25, end/4000))
      } else {
        ## both: apply the acute changes segment-wise with the oral events on
        m <- zero_re(b$mod)
        ts <- vapply(ch, function(x) x$t, numeric(1))
        segs <- c(0, sort(ts), end); st <- b$init; outs <- list()
        for (i in seq_len(length(segs) - 1)) {
          if (i > 1) m <- param(m, ch[[i - 1]]$par)
          e2 <- as.data.frame(ev)
          e2 <- e2[e2$time >= segs[i] & e2$time < segs[i + 1], , drop = FALSE]
          if (nrow(e2)) e2$time <- e2$time - segs[i]
          m2 <- init(m, st)
          o  <- if (nrow(e2))
                  mrgsim_df(m2, events = e2, end = segs[i+1] - segs[i],
                            delta = max(0.25, end/4000), recsort = 3)
                else
                  mrgsim_df(m2, end = segs[i+1] - segs[i],
                            delta = max(0.25, end/4000), recsort = 3)
          o$time <- o$time + segs[i]
          st <- as.list(o[nrow(o), intersect(names(o),
                                             names(as.list(init(mod))))])
          outs[[i]] <- if (i == 1) o else o[-1, ]
        }
        do.call(rbind, outs)
      }
    })
  }, ignoreNULL = FALSE)

  ## Untreated comparator, same patient
  run_ref <- eventReactive(input$go, {
    b <- base_mod()
    mrgsim_df(zero_re(b$mod), end = input$end*24,
              delta = max(0.25, input$end*24/4000), recsort = 3)
  }, ignoreNULL = FALSE)

  ## ---------------------------------------------------------------- tab 1
  output$t_rest <- renderTbl({
    b <- base_mod(); o <- mrgsim_df(zero_re(param(b$mod, AMP_TOT = 0)),
                                    end = 0, delta = 1, recsort = 3)
    r <- o[1, ]
    tbl_out(data.frame(
      Segment  = c("Apex", "Mid", "Base"),
      mass_frac = c(0.28, 0.34, 0.38),
      betaAR_density = c(input$rho_ap, 1.15, 1.00),
      beta2_fraction = c(input$fb2_ap, 0.32, input$fb2_bs),
      innervation = c(input$inn_ap, 0.90, 1.20),
      wall_thickness = c(0.36, 0.72, 1.00),
      resting_Gs = round(c(r$GS_AP, NA, r$GS_BS), 4),
      resting_beta2_occ = round(c(r$O2_AP, NA, r$O2_BS), 4),
      resting_shortening = round(c(r$SH_AP, r$SH_MD, r$SH_BS), 3)))
  })

  output$p_surge <- renderPlot({
    t <- seq(0, 96, 0.25)
    g <- 1/(1 + input$e2/30)
    s <- input$amp*g*exp(-t/input$tau)
    df <- data.frame(time = t, adrenaline = input$frace*s,
                     noradrenaline = (1 - input$frace)*s)
    tsplot(df, c("adrenaline", "noradrenaline"),
           c("adrenaline secretion", "noradrenaline release"),
           c("#c0392b", "#e07b39"), xunit = "h",
           title = "Secretion input actually applied",
           subtitle = sprintf("oestradiol gain %.2f (E2 = %.0f pg/mL)", g, input$e2),
           ylab = "nmol/L/h")
  })

  ## ---------------------------------------------------------------- tab 2
  output$p_cat <- renderPlot(tsplot(run(), c("EPI", "NE", "NEI"),
    c("plasma adrenaline", "plasma noradrenaline", "interstitial NE"),
    c("#c0392b", "#e07b39", "#8a4a10"), xunit = "h",
    title = "Catecholamines", ylab = "nmol/L"))

  output$p_occ <- renderPlot(tsplot(run(), c("O2_AP", "O2_BS"),
    c("apical beta2 occupancy", "basal beta2 occupancy"),
    c(SEGCOL[["Apex"]], SEGCOL[["Base"]]), xunit = "h",
    title = "beta2 occupancy", ylab = "fraction", hline = 0.30,
    subtitle = "dashed line = the Gi trafficking threshold THR_SW"))

  output$p_refl <- renderPlot(tsplot(run(), c("NEI", "CO"),
    c("interstitial NE (nmol/L)", "cardiac output (L/min)"),
    c("#8a4a10", "#1d6fb8"),
    title = "The reflex limb", ylab = "",
    subtitle = "interstitial NE tracks the forward-flow deficit, not the clock"))

  ## ---------------------------------------------------------------- tab 3
  output$p_phos <- renderPlot(tsplot(run(),
    c("PHOS_AP", "PHOS_MD", "PHOS_BS"), c("Apex", "Mid", "Base"),
    unname(SEGCOL), title = "Gi-coupled beta2 fraction (PHOS)",
    ylab = "fraction", subtitle = "same equation in all three segments"))
  output$p_gi <- renderPlot(tsplot(run(), c("GI_AP", "GI_BS"),
    c("Apex", "Base"), unname(SEGCOL[c(1, 3)]),
    title = "Gi coupling = density x occupancy x PHOS",
    ylab = "arbitrary units",
    subtitle = "the density factor is where the region enters"))
  output$p_camp <- renderPlot(tsplot(run(),
    c("CAMP_AP", "CAMP_MD", "CAMP_BS"), c("Apex", "Mid", "Base"),
    unname(SEGCOL), title = "Segmental cAMP", ylab = "normalised", hline = 1))
  output$p_qss <- renderPlot(tsplot(run(), c("CAMP_AP", "CQ_AP"),
    c("integrated cAMP", "quasi-steady-state cAMP"),
    c("#c2185b", "#333333"),
    title = "Is the QSS approximation good?", ylab = "normalised",
    subtitle = "the counterfactuals use the QSS form; the two must coincide"))

  scan_res <- eventReactive(input$scan, withProgress(
    message = "ramping adrenaline", value = 0.5, threshold_scan()))
  output$t_scan <- renderTbl(tbl_out(scan_res()))

  ## ---------------------------------------------------------------- tab 4
  output$p_sh <- renderPlot(tsplot(run(), c("SH_AP", "SH_MD", "SH_BS"),
    c("Apex", "Mid", "Base"), unname(SEGCOL),
    title = "Segmental fractional shortening", ylab = "fraction",
    hline = 0, subtitle = "below zero = systolic bulging"))
  output$p_ball <- renderPlot(tsplot(run(), c("BALL"), "ballooning index",
    "#2a6099", title = "Ballooning index", ylab = "basal minus apical",
    hline = 0))
  output$p_cont <- renderPlot(tsplot(run(),
    c("CONT_AP", "CONT_BS", "CASENS"),
    c("apical contractility", "basal contractility",
      "myofilament Ca sensitisation"),
    c(SEGCOL[["Apex"]], SEGCOL[["Base"]], "#2e8b57"),
    title = "Contractility and the drug that bypasses cAMP", ylab = "-",
    hline = 1))
  output$t_mech <- renderTbl({
    o <- run(); i <- which.min(o$LVEF)
    tbl_out(data.frame(
      quantity = c("time of LVEF nadir (h)", "apical shortening at nadir",
                   "basal shortening at nadir", "ballooning index at nadir",
                   "variant code (+1 apical / 0 none / -1 reverse)",
                   "wall motion score index at nadir"),
      value = c(round(o$time[i], 1), round(o$SH_AP[i], 3),
                round(o$SH_BS[i], 3), round(o$BALL[i], 3),
                o$VARIANT[i], round(o$WMSI[i], 2))))
  })

  ## ---------------------------------------------------------------- tab 5
  output$p_ef <- renderPlot({
    o <- run(); r <- run_ref()
    df <- data.frame(time = o$time, treated = o$LVEF,
                     untreated = r$LVEF[seq_len(nrow(o))])
    tsplot(df, c("untreated", "treated"), c("untreated", "this arm"),
           unname(ARMCOL[c(1, 2)]), title = "LVEF", ylab = "%",
           hline = 50, subtitle = "dashed line = 50%")
  })
  output$p_wmsi <- renderPlot(tsplot(run(), c("WMSI", "GLS"),
    c("wall motion score index", "global longitudinal strain (%)"),
    c("#2a6099", "#4b8b3b"), title = "Regional function", ylab = ""))
  output$t_end <- renderTbl({
    o <- run(); o$scenario <- "this arm"; o$label <- "interactive"
    tbl_out(t(summarise_arm(o)))
  })

  ## ---------------------------------------------------------------- tab 6
  output$p_hemo <- renderPlot(tsplot(run(), c("CO", "MAP", "SV"),
    c("cardiac output (L/min)", "mean arterial pressure (mmHg)",
      "stroke volume (mL)"), c("#1d6fb8", "#3b6fa8", "#7aa8d8"),
    title = "Forward flow", ylab = ""))
  output$p_grad <- renderPlot(tsplot(run(), c("GRAD", "FFLOSS"),
    c("LVOT gradient (mmHg)", "fractional forward-flow loss"),
    c("#b5651d", "#c0392b"), title = "Dynamic obstruction", ylab = ""))
  output$p_cong <- renderPlot(tsplot(run(), c("PCWP"),
    "pulmonary capillary wedge pressure", "#1f4e79",
    title = "Congestion", ylab = "mmHg", hline = 18))

  ## ---------------------------------------------------------------- tab 7
  output$p_tni <- renderPlot(tsplot(run(), "TNI", "hs-cTnI",
    "#b03a4b", title = "Troponin", ylab = "ng/mL"))
  output$p_bnp <- renderPlot(tsplot(run(), "NTBNP", "NT-proBNP",
    "#1f4e79", title = "NT-proBNP", ylab = "pg/mL"))
  output$p_ratio <- renderPlot(tsplot(run(), "BNPTNI",
    "NT-proBNP : troponin", "#333333",
    title = "The discordance ratio", ylab = "pg/mL per ng/mL"))
  output$p_infl <- renderPlot(tsplot(run(), c("IL6", "CRP"),
    c("IL-6 (pg/mL)", "CRP (mg/L)"), c("#a94064", "#d17ba0"),
    title = "Inflammation", ylab = ""))

  ## ---------------------------------------------------------------- tab 8
  output$p_qtc <- renderPlot(tsplot(run(), "QTC", "QTc", "#5a5aa8",
    title = "QTc", ylab = "ms", hline = 500,
    subtitle = "dashed line = 500 ms"))
  output$p_tdp <- renderPlot(tsplot(run(), "PTDP",
    "cumulative probability of torsade", "#c0392b",
    title = "Torsade hazard", ylab = "probability"))
  output$p_k <- renderPlot(tsplot(run(), "OED_AP", "apical oedema",
    "#fadadd", title = "The slow driver behind the QTc", ylab = "fraction"))

  ## ---------------------------------------------------------------- tab 9
  output$p_stasis <- renderPlot(tsplot(run(), "STASIS", "apical stasis index",
    "#a0522d", title = "Stasis is computed, not assumed", ylab = "-"))
  output$p_thr <- renderPlot(tsplot(run(), c("THR", "PEMB"),
    c("LV thrombus (arbitrary mass)", "cumulative embolic probability"),
    c("#a0522d", "#c0392b"), title = "Thrombus and embolism", ylab = ""))
  thr_res <- eventReactive(input$thrgo, withProgress(
    message = "running arms", value = 0.5, thrombus_arms()))
  output$t_thr <- renderTbl(tbl_out(thr_res()))

  ## --------------------------------------------------------------- tab 10
  md_res <- eventReactive(input$mdgo, withProgress(
    message = "titrating to matched drive", value = 0.4, {
      m <- matched_drive()
      rbind(cbind(agent = "dobutamine", rate = round(m$rate_dob, 1), m$dob),
            cbind(agent = "levosimendan", rate = round(m$rate_lev, 3), m$lev))
    }))
  output$t_md <- renderTbl(tbl_out(md_res()))

  ino_res <- eventReactive(input$inogo, withProgress(
    message = "running the inotrope panel", value = 0.4, {
      ids <- c("S02", "S08", "S09", "S10", "S11", "S12", "S13")
      do.call(rbind, lapply(ids, function(i) {
        o <- run_scenario(i, end = 720, delta = 1)
        summarise_arm(o)[, c("scenario", "label", "LVEF_nadir", "BALL_peak",
                             "GRAD_peak", "CO_min", "TNI_peak", "SHR_GI")]
      }))
    }))
  output$t_ino <- renderTbl(tbl_out(ino_res()))

  ## --------------------------------------------------------------- tab 11
  f_res <- eventReactive(input$fgo, withProgress(
    message = "running the falsification set", value = 0.4, falsification()))
  output$t_fals <- renderTbl(tbl_out(f_res()))

  ## --------------------------------------------------------------- tab 12
  output$p_cf <- renderPlot({
    o <- run()
    df <- data.frame(time = o$time, actual = 100*o$LVEF/100,
                     no_Gi = 100*o$EF_NOGI, no_oedema = 100*o$EF_NOED,
                     normal_energetics = 100*o$EF_NATP)
    df$actual <- o$LVEF
    tsplot(df, c("actual", "no_Gi", "no_oedema", "normal_energetics"),
           c("actual", "if Gi were abolished", "if oedema were zero",
             "if energetics were normal"),
           c("#333333", "#c2185b", "#a94064", "#a1743b"),
           title = "Counterfactual ejection fractions", ylab = "%")
  })
  hd_res <- eventReactive(input$hdgo, withProgress(
    message = "running the headroom set", value = 0.4, headroom()))
  output$t_hd <- renderTbl(tbl_out(hd_res()))

  ## --------------------------------------------------------------- tab 13
  pop_res <- eventReactive(input$popgo, withProgress(
    message = "running the cohort", value = 0.4, population(n = input$npop)))
  output$p_pop <- renderPlot({
    p <- pop_res()
    if (HAS_GG) {
      ggplot(p, aes(GRAD_peak, LVEF_nadir)) +
        geom_point(alpha = 0.6, colour = "#1f6feb") +
        labs(title = "Peak LVOT gradient vs LVEF nadir in the cohort",
             subtitle = "the gradient spread comes from septal geometry alone",
             x = "peak LVOT gradient (mmHg)", y = "LVEF nadir (%)") + THEME
    } else {
      plot(p$GRAD_peak, p$LVEF_nadir, pch = 16, col = "#1f6feb",
           xlab = "peak LVOT gradient (mmHg)", ylab = "LVEF nadir (%)",
           main = "Cohort")
    }
  })
  output$t_pop <- renderTbl({
    p <- pop_res()
    nf <- attr(p, "n_failed")
    q <- function(x) round(quantile(x, c(0.05, 0.5, 0.95), na.rm = TRUE), 2)
    tbl_out(data.frame(
      variable = c("LVEF nadir (%)", "LVEF day 30 (%)", "peak gradient (mmHg)",
                   "ballooning index", "peak troponin (ng/mL)", "peak QTc (ms)"),
      note = c(sprintf("%d retained, %d dropped as outside the domain of validity",
                       nrow(p), if (is.null(nf)) 0L else nf), "", "", "", "", ""),
      p5  = c(q(p$LVEF_nadir)[1], q(p$LVEF_d30)[1], q(p$GRAD_peak)[1],
              q(p$BALL_peak)[1], q(p$TNI_peak)[1], q(p$QTC_peak)[1]),
      median = c(q(p$LVEF_nadir)[2], q(p$LVEF_d30)[2], q(p$GRAD_peak)[2],
                 q(p$BALL_peak)[2], q(p$TNI_peak)[2], q(p$QTC_peak)[2]),
      p95 = c(q(p$LVEF_nadir)[3], q(p$LVEF_d30)[3], q(p$GRAD_peak)[3],
              q(p$BALL_peak)[3], q(p$TNI_peak)[3], q(p$QTC_peak)[3])))
  })
  rec_res <- eventReactive(input$recgo, withProgress(
    message = "running one year", value = 0.4, recurrence_1y()))
  output$t_rec <- renderTbl(tbl_out(rec_res()))

  ## --------------------------------------------------------------- tab 14
  output$t_cal <- renderPrint({
    input$calgo
    objective(log(fit_pars), verbose = TRUE)
  })
}

shinyApp(ui, server)
