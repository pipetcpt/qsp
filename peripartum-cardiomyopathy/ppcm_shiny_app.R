##  Peripartum Cardiomyopathy (PPCM) — Shiny Dashboard
##  ============================================================================
##  Companion interface for ppcm_mrgsolve_model.R (30-ODE QSP model).
##
##  Run with:
##      shiny::runApp("ppcm_shiny_app.R")
##  or place beside the model file and:
##      library(shiny); library(mrgsolve); library(dplyr); library(tidyr)
##      library(ggplot2); library(purrr)
##      shiny::runApp()
##
##  Eight tabs:
##    1. Patient profile      — susceptibility, delivery timing, presentation
##    2. Drug PK / exposure   — bromocriptine plasma + all HF-agent exposures
##    3. Prolactin axis       — 23K/16K prolactin, cathepsin D, ROS, lactation
##    4. Microvascular & miR  — sFlt-1, free VEGF, capillary density, miR-146a
##    5. LV function          — LVEF, LVEDV/LVEDD, stroke volume, cardiac output
##    6. Clinical endpoints   — NYHA, 6MWT, NT-proBNP, recovery flag, MACE hazard
##    7. Scenario comparison  — all nine prebuilt regimens side by side
##    8. Biomarker panel      — tracked biomarkers on a common normalised scale
##
##  The app deliberately surfaces the model's two hard pharmacological gates:
##  when the patient is still ANTEPARTUM, the ACEi/ARB and bromocriptine effect
##  traces are pinned at zero regardless of the dosing you request, and the
##  plots are annotated with the delivery day.
##
##  Disclaimer: research / education / hypothesis generation only.
##  ============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

MODEL_FILE <- "ppcm_mrgsolve_model.R"
mod <- mread_cache(MODEL_FILE)

## ---------------------------------------------------------------------------
## Dosing builders (mirror the helpers documented in the model file)
## ---------------------------------------------------------------------------

ppcm_brc <- function(regimen = c("none", "1wk", "8wk"), t0 = 0) {
  regimen <- match.arg(regimen)
  if (regimen == "none") return(NULL)
  if (regimen == "1wk") {
    return(ev(amt = 1, cmt = "BRC_G", time = t0, ii = 0.5, addl = 13))
  }
  ev(amt = 1, cmt = "BRC_G", time = t0,      ii = 0.5, addl = 27) +
    ev(amt = 1, cmt = "BRC_G", time = t0 + 14, ii = 1.0, addl = 41)
}

ppcm_hf <- function(end = 365, t0 = 0, bb = TRUE, rasi = TRUE, mra = TRUE,
                    diu = TRUE, sg = FALSE, vaso = FALSE, ac = FALSE,
                    dig = FALSE) {
  n <- max(1, floor(end - t0))
  e <- NULL
  add <- function(e, on, cmt, ii = 1, addl = n - 1) {
    if (!on) return(e)
    d <- ev(amt = 1, cmt = cmt, time = t0, ii = ii, addl = addl)
    if (is.null(e)) d else e + d
  }
  e <- add(e, bb,   "BB");   e <- add(e, rasi, "RASI")
  e <- add(e, mra,  "MRA");  e <- add(e, diu,  "DIU")
  e <- add(e, sg,   "SG");   e <- add(e, ac,   "AC")
  e <- add(e, dig,  "DIG")
  e <- add(e, vaso, "VD", ii = 1 / 3, addl = 3 * n - 1)
  e
}

## Nine prebuilt scenarios, used by the comparison tab
SCENARIOS <- list(
  "1. No therapy (natural history)" =
    list(par = list(BF = 1, T_DELIVERY = -7, SEV = 0.72), ev = NULL),
  "2. Standard HF therapy only" =
    list(par = list(BF = 1, T_DELIVERY = -7, SEV = 0.72), hf = list()),
  "3. Standard + bromocriptine 1 week + LMWH" =
    list(par = list(BF = 1, T_DELIVERY = -7, SEV = 0.72),
         hf = list(ac = TRUE), brc = "1wk"),
  "4. BOARD: standard + bromocriptine 8 weeks + LMWH" =
    list(par = list(BF = 1, T_DELIVERY = -7, SEV = 0.72),
         hf = list(ac = TRUE), brc = "8wk"),
  "5. BOARD + SGLT2 inhibitor" =
    list(par = list(BF = 1, T_DELIVERY = -7, SEV = 0.72),
         hf = list(ac = TRUE, sg = TRUE), brc = "8wk"),
  "6. TTN truncating-variant carrier on BOARD" =
    list(par = list(BF = 1, T_DELIVERY = -7, SEV = 0.72, GEN_TTN = 1),
         hf = list(ac = TRUE), brc = "8wk"),
  "7. Antepartum preeclamptic (ACEi gated off until day 14)" =
    list(par = list(BF = 1, T_DELIVERY = 14, SEV = 0.72, PREECL = 1),
         hf = list(mra = FALSE, vaso = TRUE, ac = TRUE), brc = "8wk",
         brc_t0 = 14),
  "8. Severe, BOARD WITHOUT anticoagulation" =
    list(par = list(BF = 1, T_DELIVERY = -7, SEV = 0.90),
         hf = list(ac = FALSE), brc = "8wk"),
  "9. Severe, BOARD WITH anticoagulation" =
    list(par = list(BF = 1, T_DELIVERY = -7, SEV = 0.90),
         hf = list(ac = TRUE), brc = "8wk")
)

run_scenario <- function(spec, end_d = 365) {
  m <- mod
  if (length(spec$par)) m <- do.call(param, c(list(m), spec$par))
  e <- NULL
  if (!is.null(spec$hf)) {
    e <- do.call(ppcm_hf, c(list(end = end_d), spec$hf))
  }
  if (!is.null(spec$brc)) {
    b <- ppcm_brc(spec$brc, t0 = if (is.null(spec$brc_t0)) 0 else spec$brc_t0)
    e <- if (is.null(e)) b else e + b
  }
  if (is.null(e)) e <- ev(amt = 0, cmt = "BB", time = 0)
  mrgsim(m, e, end = end_d, delta = 0.5) %>% as_tibble()
}

## ---------------------------------------------------------------------------
## Plot helpers
## ---------------------------------------------------------------------------

THEME <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.background = element_rect(fill = "grey92"),
        plot.title = element_text(face = "bold"))

## vertical marker for the delivery event
delivery_marker <- function(td) {
  if (is.null(td) || td <= 0) return(NULL)
  list(
    geom_vline(xintercept = td, linetype = "dashed", colour = "firebrick"),
    annotate("text", x = td, y = Inf, label = " delivery", hjust = 0,
             vjust = 1.5, size = 3.2, colour = "firebrick")
  )
}

long_plot <- function(df, vars, labels, title, ylab, td = NULL, logy = FALSE) {
  d <- df %>%
    select(time, all_of(vars)) %>%
    pivot_longer(-time, names_to = "var", values_to = "value") %>%
    mutate(var = factor(var, levels = vars, labels = labels))
  p <- ggplot(d, aes(time, value, colour = var)) +
    geom_line(linewidth = 0.8) +
    labs(title = title, x = "Days since diagnosis", y = ylab) +
    THEME
  if (logy) p <- p + scale_y_log10()
  if (!is.null(td)) p <- p + delivery_marker(td)
  p
}

facet_plot <- function(df, vars, labels, title, td = NULL) {
  d <- df %>%
    select(time, all_of(vars)) %>%
    pivot_longer(-time, names_to = "var", values_to = "value") %>%
    mutate(var = factor(var, levels = vars, labels = labels))
  p <- ggplot(d, aes(time, value)) +
    geom_line(linewidth = 0.8, colour = "steelblue4") +
    facet_wrap(~var, scales = "free_y") +
    labs(title = title, x = "Days since diagnosis", y = NULL) +
    THEME + theme(legend.position = "none")
  if (!is.null(td)) p <- p + delivery_marker(td)
  p
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Peripartum Cardiomyopathy — QSP / PK-PD Explorer"),
  tags$p(tags$em(
    "16-kDa prolactin (vasoinhibin) / sFlt-1 two-hit anti-angiogenic model with ",
    "exosomal miR-146a suppression of cardiomyocyte survival signalling. ",
    "Research and teaching use only — not for clinical decision-making."
  )),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("SEV", "Presentation severity (0-1)", 0.2, 1.0, 0.72, 0.02),
      helpText("0.72 corresponds to a presenting LVEF of about 29%."),
      sliderInput("T_DELIVERY", "Delivery day (relative to diagnosis)",
                  -30, 30, -7, 1),
      helpText(tags$b("Negative = already delivered."),
               " Positive means an antepartum diagnosis: ACEi/ARB and",
               " bromocriptine effects are gated OFF until this day."),
      checkboxInput("GEN_TTN", "TTN (or other) truncating variant carrier", FALSE),
      checkboxInput("PREECL", "Preeclampsia / hypertensive disorder", FALSE),
      checkboxInput("BF", "Breastfeeding at presentation", TRUE),
      sliderInput("LOSS_STAT3", "STAT3 / MnSOD antioxidant deficit", 0, 0.8, 0.45, 0.05),

      h4("Bromocriptine"),
      radioButtons("brc", NULL,
                   c("None" = "none",
                     "1 week (2.5 mg BID x 7 d)" = "1wk",
                     "8 weeks (BID x 2 wk, then OD x 6 wk)" = "8wk"),
                   selected = "8wk"),

      h4("Heart-failure therapy"),
      checkboxInput("bb",   "Beta-blocker", TRUE),
      checkboxInput("rasi", "ACE inhibitor / ARB / ARNI", TRUE),
      checkboxInput("mra",  "MRA (spironolactone / eplerenone)", TRUE),
      checkboxInput("diu",  "Loop diuretic", TRUE),
      checkboxInput("sg",   "SGLT2 inhibitor", FALSE),
      checkboxInput("vaso", "Hydralazine + nitrate (antepartum-safe)", FALSE),
      checkboxInput("ac",   "Anticoagulation (LMWH / warfarin)", TRUE),
      checkboxInput("dig",  "Digoxin", FALSE),

      h4("Simulation"),
      sliderInput("end_d", "Horizon (days)", 90, 730, 365, 5),
      actionButton("go", "Simulate", class = "btn-primary"),
      hr(),
      helpText(tags$b("Note:"), " bromocriptine suppresses lactation, and in this",
               " model involution is irreversible once complete — that is what",
               " separates the 1-week from the 8-week regimen, together with how",
               " much of the early high-oxidative-stress window each covers.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## ------------------------------------------------------- 1
        tabPanel(
          "1. Patient profile",
          br(),
          fluidRow(
            column(6, h4("Presentation state"), tableOutput("tbl_presentation")),
            column(6, h4("Active gates and therapy"), tableOutput("tbl_gates"))
          ),
          hr(),
          h4("Susceptibility and load drivers"),
          plotOutput("p_profile", height = "380px"),
          helpText("Load index is the peripartum haemodynamic multiplier: elevated",
                   " through late gestation, spiking after delivery (uterine",
                   " autotransfusion) and decaying with a 10-day time constant.")
        ),

        ## ------------------------------------------------------- 2
        tabPanel(
          "2. Drug PK / exposure",
          br(),
          h4("Bromocriptine plasma exposure and D2-mediated effect"),
          plotOutput("p_brc", height = "320px"),
          hr(),
          h4("Heart-failure agent exposures"),
          plotOutput("p_hfpk", height = "340px"),
          helpText("Exposures are in dose-proportional arbitrary units. Only",
                   " bromocriptine has a structured depot plus two-compartment",
                   " plasma model; the other agents are single exposure",
                   " compartments driving Emax/EC50 effects.")
        ),

        ## ------------------------------------------------------- 3
        tabPanel(
          "3. Prolactin axis",
          br(),
          h4("The PPCM switch: prolactin substrate x activated cathepsin D"),
          plotOutput("p_prl", height = "420px"),
          hr(),
          h4("Lactation state"),
          plotOutput("p_lact", height = "260px"),
          helpText("16-kDa prolactin is generated ONLY when cathepsin D is",
                   " activated above its physiological baseline. A normal",
                   " lactating woman therefore has abundant 23-kDa prolactin",
                   " and no vasoinhibin — which is why breastfeeding alone is",
                   " not pathogenic in this model.")
        ),

        ## ------------------------------------------------------- 4
        tabPanel(
          "4. Microvascular & miR-146a",
          br(),
          h4("Two-hit anti-angiogenic milieu and capillary density"),
          plotOutput("p_angio", height = "420px"),
          hr(),
          h4("Exosomal miR-146a and cardiomyocyte survival signalling"),
          plotOutput("p_mir", height = "300px"),
          helpText("sFlt-1 is produced only while the placenta is in situ, so it",
                   " disappears after delivery; the damage it does to capillary",
                   " density persists. Survival signalling is the composite of",
                   " Erbb4, Nras and Notch1.")
        ),

        ## ------------------------------------------------------- 5
        tabPanel(
          "5. LV function",
          br(),
          h4("Left ventricular ejection fraction"),
          plotOutput("p_lvef", height = "340px"),
          hr(),
          h4("Chamber geometry and output"),
          plotOutput("p_hemo", height = "380px"),
          helpText("The dashed line marks the LVEF >= 50% recovery threshold used",
                   " in the trial literature.")
        ),

        ## ------------------------------------------------------- 6
        tabPanel(
          "6. Clinical endpoints",
          br(),
          fluidRow(
            column(6, h4("Endpoint summary"), tableOutput("tbl_endpoints")),
            column(6, h4("Tissue outcome"), tableOutput("tbl_tissue"))
          ),
          hr(),
          plotOutput("p_endpoints", height = "420px"),
          helpText("MACE hazard is an illustrative algebraic surrogate, not a",
                   " fitted survival model. LV thrombus propensity rises with",
                   " intracavitary stasis (LVEF < 35%) and while bromocriptine",
                   " is on board; anticoagulation acts on it and on nothing else.")
        ),

        ## ------------------------------------------------------- 7
        tabPanel(
          "7. Scenario comparison",
          br(),
          h4("Nine prebuilt regimens"),
          actionButton("go_sc", "Run all scenarios", class = "btn-primary"),
          helpText("Runs the nine scenarios documented in the model file. The",
                   " sidebar controls do not affect this tab."),
          br(),
          plotOutput("p_scen_lvef", height = "460px"),
          hr(),
          h4("LVEF at 30 / 90 / 180 / 365 days"),
          tableOutput("tbl_scen"),
          hr(),
          h4("Irreversible scar and LV thrombus by scenario"),
          plotOutput("p_scen_scar", height = "360px")
        ),

        ## ------------------------------------------------------- 8
        tabPanel(
          "8. Biomarker panel",
          br(),
          h4("Tracked biomarkers"),
          plotOutput("p_biomarkers", height = "520px"),
          hr(),
          h4("NT-proBNP (log scale)"),
          plotOutput("p_bnp", height = "300px"),
          helpText("Serum miR-146a is the one biomarker in this panel that is",
                   " reported to be specific to PPCM rather than to heart failure",
                   " in general, and to fall with bromocriptine.")
        )
      )
    )
  )
)

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------

server <- function(input, output, session) {

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    td <- input$T_DELIVERY
    m <- mod %>% param(
      SEV        = input$SEV,
      T_DELIVERY = td,
      GEN_TTN    = as.numeric(input$GEN_TTN),
      PREECL     = as.numeric(input$PREECL),
      BF         = as.numeric(input$BF),
      LOSS_STAT3 = input$LOSS_STAT3
    )
    # bromocriptine cannot start before delivery
    brc_t0 <- max(0, td)
    e <- ppcm_hf(end = input$end_d, bb = input$bb, rasi = input$rasi,
                 mra = input$mra, diu = input$diu, sg = input$sg,
                 vaso = input$vaso, ac = input$ac, dig = input$dig)
    b <- ppcm_brc(input$brc, t0 = brc_t0)
    if (!is.null(b)) e <- if (is.null(e)) b else e + b
    if (is.null(e)) e <- ev(amt = 0, cmt = "BB", time = 0)
    mrgsim(m, e, end = input$end_d, delta = 0.5) %>% as_tibble()
  })

  td_val <- reactive(input$T_DELIVERY)

  at_day <- function(df, d, col) {
    if (max(df$time) < d) return(NA_real_)
    df[[col]][which.min(abs(df$time - d))]
  }

  ## ------------------------------------------------ tab 1
  output$tbl_presentation <- renderTable({
    d <- sim()
    tibble(
      Quantity = c("LVEF at presentation (%)", "LVEDV (mL)", "LVEDD (cm)",
                   "NT-proBNP (pg/mL)", "NYHA class", "Cardiac output (L/min)",
                   "23-kDa prolactin (ng/mL)", "16-kDa prolactin (au)",
                   "Capillary density (fraction)"),
      Value = c(sprintf("%.1f", d$LVEF[1]), sprintf("%.0f", d$LVEDV[1]),
                sprintf("%.1f", d$LVEDD_cm[1]), sprintf("%.0f", d$NTproBNP[1]),
                sprintf("%.1f", d$NYHA[1]), sprintf("%.2f", d$Cardiac_output[1]),
                sprintf("%.0f", d$Prolactin_23kDa[1]),
                sprintf("%.2f", d$Prolactin_16kDa[1]),
                sprintf("%.2f", d$Capillary_density[1]))
    )
  }, digits = 2)

  output$tbl_gates <- renderTable({
    td <- td_val()
    ante <- td > 0
    tibble(
      Item = c("Delivery day", "Antepartum at t = 0",
               "ACEi/ARB/ARNI effect", "Bromocriptine effect",
               "Bromocriptine regimen", "Anticoagulation",
               "Antepartum-safe vasodilator"),
      Status = c(
        sprintf("%+d", td),
        if (ante) "YES" else "no (postpartum)",
        if (ante && input$rasi) "DOSED BUT GATED OFF until delivery"
        else if (input$rasi) "active" else "not prescribed",
        if (ante && input$brc != "none") sprintf("starts on day %d", td)
        else if (input$brc != "none") "active" else "not prescribed",
        switch(input$brc, none = "none", `1wk` = "1 week", `8wk` = "8 weeks"),
        if (input$ac) "on" else tags$b("OFF"),
        if (input$vaso) "on" else "off"
      )
    )
  }, sanitize.text.function = identity)

  output$p_profile <- renderPlot({
    facet_plot(sim(),
               c("Load_index", "Myocardial_ROS", "Cathepsin_D", "Inflammation",
                 "Neurohormonal", "Congestion"),
               c("Peripartum load index", "Myocardial ROS", "Active cathepsin D",
                 "Inflammatory burden", "Neurohormonal activation",
                 "Congestion / volume excess"),
               "Susceptibility and load drivers", td = td_val())
  })

  ## ------------------------------------------------ tab 2
  output$p_brc <- renderPlot({
    long_plot(sim(), c("BRC_C", "BRC_P", "Bromocriptine_effect"),
              c("Bromocriptine central (au)", "Peripheral (au)",
                "Fractional prolactin suppression"),
              "Bromocriptine PK and D2-mediated effect", "au / fraction",
              td = td_val())
  })

  output$p_hfpk <- renderPlot({
    facet_plot(sim(), c("BB", "RASI", "MRA", "DIU", "SG", "VD", "AC", "DIG"),
               c("Beta-blocker", "ACEi/ARB/ARNI", "MRA", "Loop diuretic",
                 "SGLT2 inhibitor", "Hydralazine+nitrate", "Anticoagulant",
                 "Digoxin"),
               "Heart-failure agent exposures (au)", td = td_val())
  })

  ## ------------------------------------------------ tab 3
  output$p_prl <- renderPlot({
    facet_plot(sim(),
               c("Prolactin_23kDa", "Cathepsin_D", "Prolactin_16kDa",
                 "Myocardial_ROS"),
               c("23-kDa prolactin (ng/mL)", "Active cathepsin D (au)",
                 "16-kDa prolactin / vasoinhibin (au)", "Myocardial ROS (au)"),
               "Prolactin cleavage axis", td = td_val())
  })

  output$p_lact <- renderPlot({
    long_plot(sim(), c("Lactation", "Bromocriptine_effect"),
              c("Lactation state (0-1)", "Bromocriptine effect"),
              "Lactation involution is irreversible once complete", "fraction",
              td = td_val())
  })

  ## ------------------------------------------------ tab 4
  output$p_angio <- renderPlot({
    facet_plot(sim(),
               c("sFlt1_level", "Free_VEGF", "Capillary_density",
                 "Irreversible_scar"),
               c("Placental sFlt-1 (au)", "Free VEGF (au)",
                 "Capillary density (fraction of normal)",
                 "Irreversible scar (fraction)"),
               "Anti-angiogenic milieu and its structural consequence",
               td = td_val())
  })

  output$p_mir <- renderPlot({
    long_plot(sim(), c("miR146a_load", "Survival_signalling",
                       "Viable_myocyte_mass"),
              c("Cardiomyocyte miR-146a (au)",
                "Erbb4/Nras/Notch1 signalling capacity",
                "Viable contractile myocyte mass"),
              "miR-146a exosome axis", "au / fraction", td = td_val())
  })

  ## ------------------------------------------------ tab 5
  output$p_lvef <- renderPlot({
    d <- sim()
    p <- ggplot(d, aes(time, LVEF)) +
      geom_line(linewidth = 1.1, colour = "steelblue4") +
      geom_hline(yintercept = 50, linetype = "dashed", colour = "darkgreen") +
      annotate("text", x = max(d$time), y = 50, label = "recovery threshold ",
               hjust = 1, vjust = -0.6, size = 3.4, colour = "darkgreen") +
      geom_hline(yintercept = 35, linetype = "dotted", colour = "firebrick") +
      annotate("text", x = max(d$time), y = 35,
               label = "thrombus / WCD threshold ", hjust = 1, vjust = -0.6,
               size = 3.4, colour = "firebrick") +
      labs(title = "LVEF trajectory", x = "Days since diagnosis",
           y = "LVEF (%)") +
      coord_cartesian(ylim = c(0, 72)) + THEME
    p + delivery_marker(td_val())
  })

  output$p_hemo <- renderPlot({
    facet_plot(sim(),
               c("LVEDV", "LVEDD_cm", "Stroke_volume", "Cardiac_output",
                 "Heart_rate", "Wall_stress"),
               c("LV end-diastolic volume (mL)", "LV end-diastolic diameter (cm)",
                 "Stroke volume (mL)", "Cardiac output (L/min)",
                 "Heart rate (bpm)", "Normalised LV wall stress"),
               "Chamber geometry and haemodynamics", td = td_val())
  })

  ## ------------------------------------------------ tab 6
  output$tbl_endpoints <- renderTable({
    d <- sim()
    days <- c(30, 90, 180, 365)
    days <- days[days <= max(d$time)]
    tibble(
      Day = days,
      `LVEF (%)` = sapply(days, at_day, df = d, col = "LVEF"),
      `NYHA` = sapply(days, at_day, df = d, col = "NYHA"),
      `6MWT (m)` = sapply(days, at_day, df = d, col = "SixMWT_m"),
      `NT-proBNP` = sapply(days, at_day, df = d, col = "NTproBNP"),
      `Recovered` = ifelse(sapply(days, at_day, df = d,
                                  col = "Recovered_flag") > 0.5, "yes", "no")
    )
  }, digits = 1)

  output$tbl_tissue <- renderTable({
    d <- sim()
    tibble(
      Quantity = c("Peak LV thrombus propensity (au)",
                   "Irreversible scar at end (fraction)",
                   "Fibrosis at end (au)",
                   "Capillary density at end (fraction)",
                   "Lowest LVEF reached (%)",
                   "LVEF at end of horizon (%)",
                   "Peak MACE hazard (per year, surrogate)"),
      Value = c(sprintf("%.2f", max(d$LV_thrombus)),
                sprintf("%.2f", tail(d$Irreversible_scar, 1)),
                sprintf("%.2f", tail(d$Fibrosis, 1)),
                sprintf("%.2f", tail(d$Capillary_density, 1)),
                sprintf("%.1f", min(d$LVEF)),
                sprintf("%.1f", tail(d$LVEF, 1)),
                sprintf("%.2f", max(d$MACE_hazard_yr)))
    )
  })

  output$p_endpoints <- renderPlot({
    facet_plot(sim(),
               c("NYHA", "SixMWT_m", "NTproBNP", "LV_thrombus",
                 "MACE_hazard_yr", "Recovered_flag"),
               c("NYHA class", "6-minute walk distance (m)",
                 "NT-proBNP (pg/mL)", "LV thrombus propensity (au)",
                 "MACE hazard (per year, surrogate)",
                 "LVEF >= 50% (recovery flag)"),
               "Clinical endpoints", td = td_val())
  })

  ## ------------------------------------------------ tab 7
  scen <- eventReactive(input$go_sc, {
    end_d <- input$end_d
    imap_dfr(SCENARIOS, function(spec, lab) {
      run_scenario(spec, end_d) %>% mutate(scenario = lab)
    })
  })

  output$p_scen_lvef <- renderPlot({
    d <- scen()
    ggplot(d, aes(time, LVEF, colour = scenario)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 50, linetype = "dashed", colour = "grey30") +
      labs(title = "LVEF by scenario", x = "Days since diagnosis",
           y = "LVEF (%)") +
      coord_cartesian(ylim = c(0, 72)) +
      guides(colour = guide_legend(ncol = 2)) + THEME
  })

  output$tbl_scen <- renderTable({
    d <- scen()
    d %>%
      group_by(scenario) %>%
      summarise(
        `d30`  = LVEF[which.min(abs(time - 30))],
        `d90`  = LVEF[which.min(abs(time - 90))],
        `d180` = LVEF[which.min(abs(time - 180))],
        `d365` = if (max(time) >= 365) LVEF[which.min(abs(time - 365))] else NA_real_,
        `peak thrombus` = max(LV_thrombus),
        `scar at end` = tail(Irreversible_scar, 1),
        .groups = "drop"
      )
  }, digits = 2)

  output$p_scen_scar <- renderPlot({
    d <- scen() %>%
      select(time, scenario, Irreversible_scar, LV_thrombus) %>%
      pivot_longer(c(Irreversible_scar, LV_thrombus),
                   names_to = "var", values_to = "value") %>%
      mutate(var = recode(var,
                          Irreversible_scar = "Irreversible scar (fraction)",
                          LV_thrombus = "LV thrombus propensity (au)"))
    ggplot(d, aes(time, value, colour = scenario)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Days since diagnosis", y = NULL,
           title = "Scenarios 8 and 9 differ ONLY in thrombus burden") +
      guides(colour = guide_legend(ncol = 2)) + THEME
  })

  ## ------------------------------------------------ tab 8
  output$p_biomarkers <- renderPlot({
    facet_plot(sim(),
               c("Prolactin_16kDa", "miR146a_load", "sFlt1_level", "Free_VEGF",
                 "Myocardial_ROS", "Inflammation", "Fibrosis", "Neurohormonal",
                 "Capillary_density"),
               c("16-kDa prolactin (au)", "miR-146a (au; PPCM-specific)",
                 "sFlt-1 (au)", "Free VEGF (au)", "Myocardial ROS (au)",
                 "Inflammatory burden (au)", "Fibrosis (au)",
                 "Neurohormonal activation (au)", "Capillary density"),
               "Biomarker panel", td = td_val())
  })

  output$p_bnp <- renderPlot({
    long_plot(sim(), "NTproBNP", "NT-proBNP",
              "NT-proBNP (log scale)", "pg/mL", td = td_val(), logy = TRUE)
  })
}

shinyApp(ui, server)
