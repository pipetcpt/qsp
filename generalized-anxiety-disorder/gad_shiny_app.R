# =============================================================================
# Generalized Anxiety Disorder QSP model -- Shiny dashboard
# =============================================================================
# Runs gad_mrgsolve_model.R interactively.  The organising idea of the whole
# app is the one equation the model is built around:
#
#       Phi = (E_amy * S_glu) / (C_pfc * (1 + e_ex*EXPECT) * I_gaba)
#
# so the app is arranged so that a user can watch each of the four factors
# move on its own clock, and then watch what the READOUT does to them.
#
# Tabs
#   1  Patient & regimen     -- who is being treated and with what
#   2  Pharmacokinetics      -- concentrations and receptor occupancy
#   3  The four factors      -- E_amy, S_glu, C_pfc, I_gaba and Phi
#   4  The four clocks       -- normalised onset curves on one axis
#   5  Clinical endpoints    -- HAM-A total, psychic and somatic, MADRS
#   6  Trial simulator       -- placebo arm, visits, measured vs pharmacological
#   7  Scenario comparison   -- up to six regimens side by side
#   8  Dose-response         -- the occupancy hyperbola and what it implies
#   9  Benzodiazepine        -- tolerance, dependence and taper rebound
#  10  Adverse effects       -- each burden with its own tolerance clock
#  11  Assay sensitivity     -- sweep site expectancy, watch the delta vanish
#  12  Relapse prevention    -- randomised withdrawal
#
# Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
# Run with:  shiny::runApp("gad_shiny_app.R")
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

MOD <- mread_cache("gad", "gad_mrgsolve_model.R")

PAL <- c("#3E8F86", "#C46A78", "#4B87B8", "#C08A3E", "#5E9E4E",
         "#8E6FB5", "#B85A8E", "#8F8B6A")

theme_gad <- function() {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "#555555"))
}

DRUGS <- c("none", "escitalopram", "venlafaxine ER", "duloxetine",
           "pregabalin", "lorazepam", "alprazolam", "buspirone",
           "quetiapine XR")

DEFAULT_DOSE <- c("none" = 0, "escitalopram" = 10, "venlafaxine ER" = 150,
                  "duloxetine" = 60, "pregabalin" = 300, "lorazepam" = 3,
                  "alprazolam" = 1.5, "buspirone" = 45, "quetiapine XR" = 150)

# ---- build a dosing data set -------------------------------------------------
drug_events <- function(drug, mg, tstart = 0, tstop = 84, id = 1) {
  if (drug == "none" || mg <= 0) return(NULL)
  spec <- switch(drug,
    "escitalopram"  = list(cmt = 1,  per_day = 1),
    "venlafaxine ER" = list(cmt = 5, per_day = 1),
    "duloxetine"    = list(cmt = 8,  per_day = 1),
    "pregabalin"    = list(cmt = 10, per_day = 2),
    "lorazepam"     = list(cmt = 12, per_day = 3),
    "alprazolam"    = list(cmt = 12, per_day = 3),
    "buspirone"     = list(cmt = 15, per_day = 3),
    "quetiapine XR" = list(cmt = 18, per_day = 1))
  ii <- 1 / spec$per_day
  n  <- max(1, floor((tstop - tstart) / ii))
  data.frame(ID = id, time = tstart, cmt = spec$cmt,
             amt = mg / spec$per_day * 1000, ii = ii, addl = n - 1, evid = 1)
}

visit_events <- function(visits, id = 1) {
  if (length(visits) == 0) return(NULL)
  data.frame(ID = id, time = visits, cmt = 50, amt = 1, ii = 0, addl = 0, evid = 1)
}

VISITS_8WK  <- c(0, 7, 14, 28, 42, 56)
VISITS_12WK <- c(VISITS_8WK, 70, 84)

bzd_params <- function(drug) {
  if (drug == "alprazolam")
    list(BZD_POT = 1.53, F_bzd = 0.88, V_bzd = 65, CL_bzd = 90)
  else
    list(BZD_POT = 1.00, F_bzd = 0.90, V_bzd = 90, CL_bzd = 108)
}

run_one <- function(drug, mg, dis, crcl, fm2d6, adherence, cbt,
                    tstop, visits, fluct_on, drug_stop = Inf, id = 1) {
  ev <- rbind(drug_events(drug, mg * adherence, 0, min(drug_stop, tstop), id),
              visit_events(visits, id))
  if (is.null(ev)) ev <- data.frame(ID = id, time = 0, cmt = 1, amt = 0,
                                    ii = 0, addl = 0, evid = 1)
  p <- c(list(DIS = dis, crcl = crcl, fm2d6 = fm2d6, CBT = as.numeric(cbt)),
         bzd_params(drug))
  m <- MOD
  if (!fluct_on) p$FLUCT0 <- 0
  m <- param(m, p)
  out <- mrgsim_d(m, ev, end = tstop, delta = 0.25) %>% as_tibble()
  out$arm <- sprintf("%s %g", drug, mg)
  out
}

# =============================================================================
ui <- navbarPage(
  title = "GAD QSP  ·  Phi = (E_amy × S_glu) / (C_pfc × I_gaba)",
  header = tags$style(HTML(".well{background:#F7F9FA;} .nav-tabs{font-size:13px;}")),

  # ---- 1 -------------------------------------------------------------------
  tabPanel("1 · Patient & regimen",
    sidebarLayout(
      sidebarPanel(width = 4,
        h4("Patient"),
        sliderInput("dis", "Disease severity DIS (1.0 = typical trial patient)",
                    0.4, 1.7, 1.0, 0.05),
        sliderInput("crcl", "Creatinine clearance (mL/min)", 15, 130, 100, 5),
        selectInput("cyp2d6", "CYP2D6 phenotype",
                    c("extensive (fm 0.70)" = 0.70,
                      "poor (fm 0.15)" = 0.15,
                      "ultrarapid (fm 0.85)" = 0.85), selected = 0.70),
        sliderInput("adh", "Adherence", 0.3, 1.0, 1.0, 0.05),
        checkboxInput("cbt", "Concurrent CBT (12 weekly sessions)", FALSE),
        hr(),
        h4("Regimen"),
        selectInput("drug", "Drug", DRUGS, selected = "escitalopram"),
        uiOutput("dose_ui"),
        sliderInput("tstop", "Simulation length (days)", 28, 560, 84, 7),
        numericInput("dstop", "Stop drug at day (Inf = never)", Inf),
        hr(),
        h4("Trial context"),
        checkboxInput("visits", "Study visits (weeks 0,1,2,4,6,8,10,12)", TRUE),
        checkboxInput("fluct", "Enrolment at a symptom peak (regression to the mean)", TRUE),
        helpText("Switch both off to see the PHARMACOLOGICAL clock alone. ",
                 "The difference between the two is the fifth clock, and it is ",
                 "usually larger than the difference between drugs.")
      ),
      mainPanel(width = 8,
        h4("The one equation"),
        verbatimTextOutput("eqn"),
        plotOutput("p_overview", height = "460px"),
        DTOutput("t_summary"))
    )),

  # ---- 2 -------------------------------------------------------------------
  tabPanel("2 · Pharmacokinetics",
    fluidRow(column(12,
      h4("Plasma concentrations"), plotOutput("p_pk", height = "320px"),
      h4("Receptor occupancy"), plotOutput("p_occ", height = "320px"),
      helpText("Occupancy anchors: SERT EC50 5 ng/mL for escitalopram (80% at ",
               "10 mg); NET EC50 260 ng/mL for venlafaxine total active moiety ",
               "(Arakawa 2019 measured 8-61% across 75-300 mg); ",
               "benzodiazepine-site EC50 96 ng/mL lorazepam-equivalent ",
               "(Atack 2007) — which is why anxiolytic doses sit near 20% ",
               "occupancy, not near saturation.")))),

  # ---- 3 -------------------------------------------------------------------
  tabPanel("3 · The four factors",
    fluidRow(column(12,
      plotOutput("p_factors", height = "420px"),
      h4("Phi and its numerator / denominator"),
      plotOutput("p_phi", height = "300px"),
      helpText("E_amy and S_glu are the numerator; C_pfc and I_gaba are the ",
               "denominator. A drug that halves the numerator and a drug that ",
               "doubles the denominator do the same thing to Phi — but not the ",
               "same thing to HAM-A, because the readout saturates.")))),

  # ---- 4 -------------------------------------------------------------------
  tabPanel("4 · The four clocks",
    sidebarLayout(
      sidebarPanel(width = 3,
        helpText("Each curve is scaled to its own eventual effect, so this axis ",
                 "shows SPEED only. Nothing here is fitted: the ordering is a ",
                 "consequence of where each drug acts."),
        checkboxGroupInput("clock_drugs", "Include",
          c("escitalopram 10", "pregabalin 300", "lorazepam 3", "quetiapine XR 150"),
          selected = c("escitalopram 10", "pregabalin 300", "lorazepam 3"))),
      mainPanel(width = 9,
        plotOutput("p_clocks", height = "420px"),
        h4("Escitalopram, step by step down the cascade"),
        plotOutput("p_cascade", height = "320px"),
        helpText("SERT occupancy is essentially complete within a day. The ",
                 "score is not, because extracellular 5-HT can only rise as the ",
                 "somatodendritic 5-HT1A autoreceptor desensitises, and the ",
                 "plasticity that raises C_pfc then has its own time constant. ",
                 "Two slow steps in series is what makes the onset SIGMOID.")))),

  # ---- 5 -------------------------------------------------------------------
  tabPanel("5 · Clinical endpoints",
    fluidRow(column(12,
      plotOutput("p_hama", height = "360px"),
      fluidRow(column(6, plotOutput("p_clusters", height = "300px")),
               column(6, plotOutput("p_other", height = "300px"))),
      DTOutput("t_endpoints")))),

  # ---- 6 -------------------------------------------------------------------
  tabPanel("6 · Trial simulator",
    sidebarLayout(
      sidebarPanel(width = 3,
        sliderInput("t_eex", "Site expectancy coefficient e_ex", 0, 1.6, 0.55, 0.05),
        sliderInput("t_fluct", "Enrolment peak (HAM-A points)", 0, 12, 6, 0.5),
        selectInput("t_drug", "Active arm", DRUGS[-1], selected = "escitalopram"),
        uiOutput("t_dose_ui")),
      mainPanel(width = 9,
        plotOutput("p_trial", height = "380px"),
        h4("Measured delta vs pharmacological change"),
        DTOutput("t_trial"),
        helpText("The placebo arm reaches roughly half of its eight-week effect ",
                 "in the first week. The fifth clock is the fastest in the ",
                 "system, which is exactly why measured onset differences ",
                 "between drugs are so much smaller than pharmacological ones.")))),

  # ---- 7 -------------------------------------------------------------------
  tabPanel("7 · Scenario comparison",
    sidebarLayout(
      sidebarPanel(width = 3,
        lapply(1:6, function(i) {
          tagList(
            selectInput(paste0("s", i, "_drug"), paste("Arm", i), DRUGS,
                        selected = c("none", "escitalopram", "pregabalin",
                                     "venlafaxine ER", "lorazepam",
                                     "quetiapine XR")[i]),
            numericInput(paste0("s", i, "_mg"), NULL,
                         c(0, 10, 300, 150, 3, 150)[i]))
        }),
        sliderInput("s_tstop", "Days", 28, 168, 84, 7)),
      mainPanel(width = 9,
        plotOutput("p_scen", height = "420px"),
        DTOutput("t_scen")))),

  # ---- 8 -------------------------------------------------------------------
  tabPanel("8 · Dose-response",
    sidebarLayout(
      sidebarPanel(width = 3,
        selectInput("dr_drug", "Drug",
                    c("escitalopram", "venlafaxine ER", "pregabalin",
                      "quetiapine XR", "duloxetine"), "escitalopram"),
        textInput("dr_doses", "Doses (comma separated, mg/d)", "5,10,20,30,40")),
      mainPanel(width = 9,
        plotOutput("p_dr", height = "380px"),
        DTOutput("t_dr"),
        helpText("For escitalopram the curve is flat because C/(C+5) is already ",
                 "0.80 at a steady-state 20 ng/mL. That is a derivation, not a ",
                 "fit. For venlafaxine the SERT arm is equally flat but the NET ",
                 "arm is not, which is where the SNRI dose-response comes from.")))),

  # ---- 9 -------------------------------------------------------------------
  tabPanel("9 · Benzodiazepine",
    sidebarLayout(
      sidebarPanel(width = 3,
        selectInput("bz_drug", "Agent", c("lorazepam", "alprazolam"), "lorazepam"),
        numericInput("bz_mg", "Daily dose (mg)", 3),
        numericInput("bz_dur", "Duration (days)", 84),
        selectInput("bz_stop", "Discontinuation",
                    c("abrupt", "4-week taper", "continue"), "abrupt")),
      mainPanel(width = 9,
        plotOutput("p_bz", height = "360px"),
        h4("The two receptor pools and the one adaptation state"),
        plotOutput("p_bzpool", height = "300px"),
        helpText("R_a1 (sedation) tolerates with a time constant of days; R_a2 ",
                 "(anxiolysis) with one of weeks. DEPEND is a single adaptation ",
                 "state that subtracts a fifth of the effect while the drug is ",
                 "on board and ALL of it when the drug is removed — the same ",
                 "parameter, opposite sign. That is a falsifiable coupling.")))),

  # ---- 10 ------------------------------------------------------------------
  tabPanel("10 · Adverse effects",
    fluidRow(column(12,
      plotOutput("p_ae", height = "400px"),
      plotOutput("p_aeburden", height = "260px"),
      helpText("Nausea tolerates in about a week, sedation in a fortnight, ",
               "sexual dysfunction does not tolerate at all. That ordering, ",
               "not the peak magnitude, is what determines whether a patient ",
               "is still taking the drug at month six.")))),

  # ---- 11 ------------------------------------------------------------------
  tabPanel("11 · Assay sensitivity",
    fluidRow(column(12,
      h4("Sweep the site expectancy coefficient; the pharmacology never changes"),
      plotOutput("p_assay", height = "400px"),
      DTOutput("t_assay"),
      helpText("Expectancy raises C_pfc, C_pfc is in the DENOMINATOR of Phi, ",
               "and HAM-A is a saturating function of Phi. A placebo-responsive ",
               "site is therefore already on the flat part of the readout when ",
               "the drug arrives. Assay sensitivity is a property of the ",
               "instrument, derived here rather than assumed.")))),

  # ---- 12 ------------------------------------------------------------------
  tabPanel("12 · Relapse prevention",
    sidebarLayout(
      sidebarPanel(width = 3,
        numericInput("rp_n", "Virtual subjects", 200),
        numericInput("rp_ol", "Open-label weeks", 12),
        numericInput("rp_mg", "Escitalopram dose (mg)", 20),
        numericInput("rp_fu", "Follow-up weeks", 76),
        numericInput("rp_seed", "Seed", 20260806),
        actionButton("rp_go", "Run", class = "btn-primary")),
      mainPanel(width = 9,
        plotOutput("p_relapse", height = "400px"),
        DTOutput("t_relapse"),
        helpText("Allgulander 2006 (PMID 16316482) reported 56% relapse on ",
                 "placebo and 19% on escitalopram over 24-76 weeks after a ",
                 "12-week open-label lead-in. Nothing in the acute-phase ",
                 "calibration knew about those numbers.")))))

# =============================================================================
server <- function(input, output, session) {

  output$dose_ui <- renderUI({
    numericInput("mg", "Total daily dose (mg)",
                 value = unname(DEFAULT_DOSE[input$drug]))
  })
  output$t_dose_ui <- renderUI({
    numericInput("t_mg", "Total daily dose (mg)",
                 value = unname(DEFAULT_DOSE[input$t_drug]))
  })

  visits_vec <- reactive(if (isTRUE(input$visits)) VISITS_12WK else numeric(0))

  sim <- reactive({
    req(input$mg)
    run_one(input$drug, input$mg, input$dis, input$crcl,
            as.numeric(input$cyp2d6), input$adh, input$cbt,
            input$tstop, visits_vec(), input$fluct,
            drug_stop = ifelse(is.finite(input$dstop), input$dstop, Inf))
  })

  output$eqn <- renderText({
    d <- sim()
    e <- tail(d, 1)
    sprintf(paste0("Phi = (E_amy x S_glu) / (C_pfc(1+e_ex*EXPECT) x I_gaba)\n",
                   "    = (%.3f x %.3f) / (%.3f x %.3f) = %.3f",
                   "   [normalised to healthy: %.3f]\n",
                   "HAM-A = %.1f  (psychic %.1f + somatic %.1f)"),
            e$EAMY, e$SGLU, e$CPFCE, e$IGABA, e$PHI, e$PHI_N,
            e$HAMA, e$HAMA_PSY, e$HAMA_SOM)
  })

  output$p_overview <- renderPlot({
    d <- sim()
    d %>% select(time, HAMA, HAMA_PSY, HAMA_SOM, PHI_N, MADRS) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) +
      labs(x = "days", y = NULL, title = "Overview",
           subtitle = "one gain, and what it is read out as") +
      theme_gad() + theme(legend.position = "none")
  })

  output$t_summary <- renderDT({
    d <- sim()
    wk <- c(0, 7, 14, 28, 56, 84)
    wk <- wk[wk <= max(d$time)]
    idx <- sapply(wk, function(t) which.min(abs(d$time - t)))
    tibble(day = wk,
           HAMA = round(d$HAMA[idx], 2),
           change = round(d$HAMA[idx] - d$HAMA[idx[1]], 2),
           psychic = round(d$HAMA_PSY[idx], 2),
           somatic = round(d$HAMA_SOM[idx], 2),
           Phi_n = round(d$PHI_N[idx], 3),
           E_amy = round(d$EAMY[idx], 3),
           S_glu = round(d$SGLU[idx], 3),
           C_pfc = round(d$CPFCE[idx], 3),
           I_gaba = round(d$IGABA[idx], 3)) %>%
      datatable(options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  output$p_pk <- renderPlot({
    d <- sim()
    d %>% select(time, escitalopram = Cesc_o, venlafaxine = Cven_o,
                 ODV = Codv_o, duloxetine = Cdlx_o, pregabalin = Cpgb_o,
                 benzodiazepine = Cbzd_o, quetiapine = Cqtp_o) %>%
      pivot_longer(-time) %>% filter(value > 1e-6) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) +
      labs(x = "days", y = "ng/mL") + theme_gad() +
      theme(legend.position = "none")
  })

  output$p_occ <- renderPlot({
    d <- sim()
    d %>% select(time, SERT = OCC_SERT, NET = OCC_NET, `alpha2delta` = OCC_A2D,
                 `BZ site` = OCC_BZ, H1 = OCC_H1) %>%
      pivot_longer(-time) %>% filter(value > 1e-6) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) + ylim(0, 1) +
      labs(x = "days", y = "fractional occupancy", colour = NULL) + theme_gad()
  })

  output$p_factors <- renderPlot({
    d <- sim()
    d %>% select(time, `E_amy (numerator)` = EAMY, `S_glu (numerator)` = SGLU,
                 `C_pfc (denominator)` = CPFCE, `I_gaba (denominator)` = IGABA) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1.1) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) +
      labs(x = "days", y = NULL,
           title = "The four factors, each on its own clock") +
      theme_gad() + theme(legend.position = "none")
  })

  output$p_phi <- renderPlot({
    d <- sim()
    d %>% mutate(numerator = EAMY * SGLU, denominator = CPFCE * IGABA) %>%
      select(time, numerator, denominator, Phi_n = PHI_N) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1.1) +
      scale_colour_manual(values = PAL) +
      labs(x = "days", y = NULL, colour = NULL) + theme_gad()
  })

  clocks <- reactive({
    specs <- list("escitalopram 10" = c("escitalopram", 10),
                  "pregabalin 300" = c("pregabalin", 300),
                  "lorazepam 3" = c("lorazepam", 3),
                  "quetiapine XR 150" = c("quetiapine XR", 150))
    sel <- input$clock_drugs
    bind_rows(lapply(sel, function(nm) {
      s <- specs[[nm]]
      run_one(s[1], as.numeric(s[2]), input$dis, input$crcl, 0.70, 1, FALSE,
              84, numeric(0), FALSE) %>% mutate(arm = nm)
    }))
  })

  output$p_clocks <- renderPlot({
    clocks() %>% group_by(arm) %>%
      mutate(frac = (HAMA - first(HAMA)) / (last(HAMA) - first(HAMA))) %>%
      ungroup() %>%
      ggplot(aes(time, frac, colour = arm)) + geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 0.5, linetype = 2, colour = "#999999") +
      scale_colour_manual(values = PAL) + xlim(0, 56) +
      labs(x = "days", y = "fraction of the drug's own eventual effect",
           colour = NULL, title = "Speed, with magnitude divided out",
           subtitle = "no trial contact: this is the pharmacological clock alone") +
      theme_gad()
  })

  output$p_cascade <- renderPlot({
    d <- run_one("escitalopram", 10, input$dis, input$crcl, 0.70, 1, FALSE,
                 84, numeric(0), FALSE)
    d %>% transmute(time,
                    `SERT occupancy` = OCC_SERT / max(OCC_SERT),
                    `autoreceptor desensitisation` = 1 - (AUTO - min(AUTO)) / (1 - min(AUTO)),
                    `extracellular 5-HT` = (SHT - first(SHT)) / (last(SHT) - first(SHT)),
                    `C_pfc` = (CPFCE - first(CPFCE)) / (last(CPFCE) - first(CPFCE)),
                    `HAM-A` = (HAMA - first(HAMA)) / (last(HAMA) - first(HAMA))) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1.1) +
      scale_colour_manual(values = PAL) +
      labs(x = "days", y = "fraction of eventual change", colour = NULL) +
      theme_gad()
  })

  output$p_hama <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, HAMA)) +
      geom_line(linewidth = 1.2, colour = PAL[1]) +
      geom_hline(yintercept = 7, linetype = 2, colour = "#5E9E4E") +
      annotate("text", x = max(d$time) * 0.8, y = 7.8, label = "remission (HAM-A 7)",
               colour = "#5E9E4E", size = 3.5) +
      geom_hline(yintercept = 15, linetype = 2, colour = "#C46A78") +
      annotate("text", x = max(d$time) * 0.8, y = 15.8, label = "relapse threshold (15)",
               colour = "#C46A78", size = 3.5) +
      labs(x = "days", y = "HAM-A total") + theme_gad()
  })

  output$p_clusters <- renderPlot({
    d <- sim()
    d %>% select(time, psychic = HAMA_PSY, somatic = HAMA_SOM) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1.1) +
      scale_colour_manual(values = PAL[c(2, 3)]) +
      labs(x = "days", y = "cluster score", colour = NULL,
           title = "Psychic and somatic clusters move at different speeds") +
      theme_gad()
  })

  output$p_other <- renderPlot({
    d <- sim()
    d %>% select(time, MADRS, `sleep deficit x5` = SLEEPD) %>%
      mutate(`sleep deficit x5` = `sleep deficit x5` * 5) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1.1) +
      scale_colour_manual(values = PAL[c(4, 6)]) +
      labs(x = "days", y = NULL, colour = NULL) + theme_gad()
  })

  output$t_endpoints <- renderDT({
    d <- sim()
    e <- tail(d, 1)
    b <- head(d, 1)
    tibble(endpoint = c("HAM-A total", "HAM-A psychic", "HAM-A somatic",
                        "MADRS", "response (>=50% reduction)",
                        "remission (HAM-A <= 7)", "weight change (kg)"),
           baseline = c(b$HAMA, b$HAMA_PSY, b$HAMA_SOM, b$MADRS, NA, NA, 0),
           final = c(e$HAMA, e$HAMA_PSY, e$HAMA_SOM, e$MADRS,
                     as.numeric(e$HAMA <= 0.5 * b$HAMA),
                     as.numeric(e$HAMA <= 7), e$WT)) %>%
      mutate(across(where(is.numeric), ~round(.x, 2))) %>%
      datatable(options = list(dom = "t"), rownames = FALSE)
  })

  trial <- reactive({
    req(input$t_mg)
    m <- param(MOD, e_ex = input$t_eex, FLUCT0 = input$t_fluct, DIS = input$dis)
    mk <- function(drug, mg, id) {
      ev <- rbind(drug_events(drug, mg, 0, 84, id), visit_events(VISITS_12WK, id))
      if (is.null(ev)) ev <- data.frame(ID = id, time = 0, cmt = 1, amt = 0,
                                        ii = 0, addl = 0, evid = 1)
      mm <- param(m, bzd_params(drug))
      mrgsim_d(mm, ev, end = 84, delta = 1) %>% as_tibble() %>%
        mutate(arm = if (drug == "none") "placebo" else sprintf("%s %g", drug, mg))
    }
    pharm <- function(drug, mg, id) {
      ev <- drug_events(drug, mg, 0, 84, id)
      if (is.null(ev)) ev <- data.frame(ID = id, time = 0, cmt = 1, amt = 0,
                                        ii = 0, addl = 0, evid = 1)
      mm <- param(m, c(bzd_params(drug), list(FLUCT0 = 0)))
      mrgsim_d(mm, ev, end = 84, delta = 1) %>% as_tibble() %>%
        mutate(arm = "pharmacological only")
    }
    bind_rows(mk("none", 0, 1), mk(input$t_drug, input$t_mg, 2),
              pharm(input$t_drug, input$t_mg, 3))
  })

  output$p_trial <- renderPlot({
    trial() %>% group_by(arm) %>% mutate(chg = HAMA - first(HAMA)) %>% ungroup() %>%
      ggplot(aes(time, chg, colour = arm)) + geom_line(linewidth = 1.2) +
      scale_colour_manual(values = PAL) +
      labs(x = "days", y = "HAM-A change from baseline", colour = NULL) +
      theme_gad()
  })

  output$t_trial <- renderDT({
    d <- trial()
    wk <- c(7, 14, 28, 56, 84)
    f <- function(a, t) {
      s <- d %>% filter(arm == a)
      s$HAMA[which.min(abs(s$time - t))] - s$HAMA[1]
    }
    arms <- unique(d$arm)
    act <- setdiff(arms, c("placebo", "pharmacological only"))
    tibble(day = wk,
           placebo = sapply(wk, f, a = "placebo"),
           active = sapply(wk, f, a = act),
           `measured delta` = active - placebo,
           `pharmacological only` = sapply(wk, f, a = "pharmacological only")) %>%
      mutate(across(where(is.numeric), ~round(.x, 2))) %>%
      datatable(options = list(dom = "t"), rownames = FALSE)
  })

  scen <- reactive({
    bind_rows(lapply(1:6, function(i) {
      dg <- input[[paste0("s", i, "_drug")]]
      mg <- input[[paste0("s", i, "_mg")]]
      req(!is.null(dg), !is.null(mg))
      run_one(dg, mg, input$dis, input$crcl, as.numeric(input$cyp2d6),
              input$adh, input$cbt, input$s_tstop, VISITS_12WK, TRUE, id = i) %>%
        mutate(arm = if (dg == "none") "placebo" else sprintf("%s %g mg", dg, mg))
    }))
  })

  output$p_scen <- renderPlot({
    scen() %>% group_by(arm) %>% mutate(chg = HAMA - first(HAMA)) %>% ungroup() %>%
      ggplot(aes(time, chg, colour = arm)) + geom_line(linewidth = 1.1) +
      scale_colour_manual(values = PAL) +
      labs(x = "days", y = "HAM-A change from baseline", colour = NULL) +
      theme_gad()
  })

  output$t_scen <- renderDT({
    scen() %>% group_by(arm) %>%
      summarise(baseline = round(first(HAMA), 2),
                `week 1` = round(HAMA[which.min(abs(time - 7))] - first(HAMA), 2),
                `week 4` = round(HAMA[which.min(abs(time - 28))] - first(HAMA), 2),
                `week 8` = round(HAMA[which.min(abs(time - 56))] - first(HAMA), 2),
                `final Phi_n` = round(last(PHI_N), 3),
                `AE burden` = round(last(AE_BURDEN), 3), .groups = "drop") %>%
      datatable(options = list(dom = "t"), rownames = FALSE)
  })

  dr <- reactive({
    doses <- as.numeric(strsplit(input$dr_doses, ",")[[1]])
    doses <- doses[!is.na(doses)]
    bind_rows(lapply(seq_along(doses), function(i) {
      run_one(input$dr_drug, doses[i], input$dis, input$crcl,
              as.numeric(input$cyp2d6), 1, FALSE, 56, VISITS_8WK, TRUE, id = i) %>%
        mutate(dose = doses[i])
    }))
  })

  output$p_dr <- renderPlot({
    s <- dr() %>% group_by(dose) %>%
      summarise(occ = max(c(last(OCC_SERT), last(OCC_A2D), last(OCC_H1))),
                chg = last(HAMA) - first(HAMA), .groups = "drop")
    p1 <- ggplot(s, aes(dose, occ)) + geom_line(colour = PAL[1], linewidth = 1.2) +
      geom_point(size = 3, colour = PAL[1]) + ylim(0, 1) +
      labs(x = "mg/d", y = "occupancy at the primary target") + theme_gad()
    p2 <- ggplot(s, aes(dose, chg)) + geom_line(colour = PAL[2], linewidth = 1.2) +
      geom_point(size = 3, colour = PAL[2]) +
      labs(x = "mg/d", y = "week-8 HAM-A change") + theme_gad()
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  })

  output$t_dr <- renderDT({
    dr() %>% group_by(dose) %>%
      summarise(Css = round(max(c(last(Cesc_o), last(Cven_o) + last(Codv_o),
                                  last(Cpgb_o), last(Cqtp_o), last(Cdlx_o))), 1),
                OCC_SERT = round(last(OCC_SERT), 3),
                OCC_NET = round(last(OCC_NET), 3),
                OCC_A2D = round(last(OCC_A2D), 3),
                `week-8 change` = round(last(HAMA) - first(HAMA), 2),
                `AE burden` = round(last(AE_BURDEN), 3), .groups = "drop") %>%
      datatable(options = list(dom = "t"), rownames = FALSE)
  })

  bz <- reactive({
    dur <- input$bz_dur
    tstop <- dur + 84
    id <- 1
    ev <- switch(input$bz_stop,
      "abrupt" = drug_events(input$bz_drug, input$bz_mg, 0, dur, id),
      "continue" = drug_events(input$bz_drug, input$bz_mg, 0, tstop, id),
      "4-week taper" = rbind(
        drug_events(input$bz_drug, input$bz_mg, 0, dur, id),
        drug_events(input$bz_drug, input$bz_mg * 0.75, dur, dur + 7, id),
        drug_events(input$bz_drug, input$bz_mg * 0.50, dur + 7, dur + 14, id),
        drug_events(input$bz_drug, input$bz_mg * 0.25, dur + 14, dur + 21, id),
        drug_events(input$bz_drug, input$bz_mg * 0.125, dur + 21, dur + 28, id)))
    m <- param(MOD, c(bzd_params(input$bz_drug),
                      list(DIS = input$dis, FLUCT0 = 0)))
    mrgsim_d(m, ev, end = tstop, delta = 0.5) %>% as_tibble()
  })

  output$p_bz <- renderPlot({
    d <- bz()
    base <- d$HAMA[1]
    ggplot(d, aes(time, HAMA)) + geom_line(linewidth = 1.2, colour = PAL[3]) +
      geom_hline(yintercept = base, linetype = 2, colour = "#999999") +
      annotate("text", x = max(d$time) * 0.05, y = base + 0.8,
               label = "untreated baseline", hjust = 0, size = 3.5,
               colour = "#777777") +
      geom_vline(xintercept = input$bz_dur, linetype = 3, colour = "#C46A78") +
      labs(x = "days", y = "HAM-A total",
           title = "One adaptation state gives tolerance on drug and rebound off it") +
      theme_gad()
  })

  output$p_bzpool <- renderPlot({
    bz() %>% select(time, `R_a1 (sedation pool)` = RA1,
                    `R_a2 (anxiolysis pool)` = RA2,
                    `DEPEND (adaptation)` = DEPEND,
                    `BZ occupancy` = OCC_BZ) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1.1) +
      scale_colour_manual(values = PAL) +
      labs(x = "days", y = NULL, colour = NULL) + theme_gad()
  })

  output$p_ae <- renderPlot({
    d <- sim()
    d %>% select(time, nausea = NAUSEA, sedation = SEDATION,
                 dizziness = DIZZINESS, activation = ACTIVATION,
                 `sexual dysfunction` = SEXD, `weight (kg)` = WT) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1.1) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) +
      labs(x = "days", y = NULL,
           title = "Each adverse effect has its own tolerance clock") +
      theme_gad() + theme(legend.position = "none")
  })

  output$p_aeburden <- renderPlot({
    d <- sim()
    d %>% select(time, `AE burden` = AE_BURDEN,
                 `probability still in trial` = SURVIVAL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1.2) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL[c(2, 5)]) +
      labs(x = "days", y = NULL) + theme_gad() +
      theme(legend.position = "none")
  })

  assay <- reactive({
    exs <- seq(0, 1.5, by = 0.15)
    bind_rows(lapply(seq_along(exs), function(i) {
      m <- param(MOD, e_ex = exs[i], DIS = input$dis)
      run <- function(drug, mg, id) {
        ev <- rbind(drug_events(drug, mg, 0, 56, id), visit_events(VISITS_8WK, id))
        if (is.null(ev)) ev <- data.frame(ID = id, time = 0, cmt = 1, amt = 0,
                                          ii = 0, addl = 0, evid = 1)
        mrgsim_d(param(m, bzd_params(drug)), ev, end = 56, delta = 7) %>% as_tibble()
      }
      p <- run("none", 0, 1)
      a <- run("escitalopram", 10, 2)
      g <- run("pregabalin", 300, 3)
      tibble(e_ex = exs[i],
             placebo = last(p$HAMA) - first(p$HAMA),
             escitalopram = last(a$HAMA) - first(a$HAMA),
             pregabalin = last(g$HAMA) - first(g$HAMA)) %>%
        mutate(`delta escitalopram` = escitalopram - placebo,
               `delta pregabalin` = pregabalin - placebo)
    }))
  })

  output$p_assay <- renderPlot({
    assay() %>% select(e_ex, `placebo change` = placebo,
                       `delta escitalopram`, `delta pregabalin`) %>%
      pivot_longer(-e_ex) %>%
      ggplot(aes(e_ex, value, colour = name)) +
      geom_line(linewidth = 1.2) + geom_point(size = 2) +
      scale_colour_manual(values = PAL) +
      labs(x = "site expectancy coefficient e_ex",
           y = "week-8 HAM-A (points)", colour = NULL,
           title = "Identical pharmacology, vanishing measured effect") +
      theme_gad()
  })

  output$t_assay <- renderDT({
    assay() %>% mutate(across(where(is.numeric), ~round(.x, 2))) %>%
      datatable(options = list(dom = "t", pageLength = 12), rownames = FALSE)
  })

  relapse <- eventReactive(input$rp_go, {
    set.seed(input$rp_seed)
    n <- input$rp_n
    dis <- pmin(pmax(rnorm(n, 1.0, 0.22), 0.45), 1.70)
    ol_days <- input$rp_ol * 7
    fu_days <- input$rp_fu * 7
    ev <- do.call(rbind, lapply(1:n, function(i)
      rbind(drug_events("escitalopram", input$rp_mg, 0, ol_days, i),
            visit_events(seq(0, ol_days, by = 14), i))))
    idata <- data.frame(ID = 1:n, DIS = dis)
    o1 <- mrgsim_d(MOD, ev, idata = idata, end = ol_days, delta = ol_days) %>%
      as_tibble() %>% filter(time == ol_days)
    hm12 <- o1$HAMA + rnorm(n, 0, 2.2)
    resp <- which(hm12 <= 10)
    if (length(resp) < 10) return(NULL)
    cont <- resp[seq(1, length(resp), by = 2)]
    pbo <- setdiff(resp, cont)
    mk2 <- function(ids, on) {
      do.call(rbind, lapply(ids, function(i) {
        e <- rbind(if (on) drug_events("escitalopram", input$rp_mg, 0, fu_days, i) else NULL,
                   visit_events(seq(0, fu_days, by = 28), i))
        if (is.null(e)) e <- data.frame(ID = i, time = 0, cmt = 1, amt = 0,
                                        ii = 0, addl = 0, evid = 1)
        e
      }))
    }
    # continue the phase-1 state by re-running with the phase-1 lead-in prefixed
    run2 <- function(ids, on) {
      ev2 <- do.call(rbind, lapply(ids, function(i)
        rbind(drug_events("escitalopram", input$rp_mg, 0, ol_days, i),
              visit_events(seq(0, ol_days, by = 14), i),
              if (on) drug_events("escitalopram", input$rp_mg, ol_days,
                                  ol_days + fu_days, i) else NULL,
              visit_events(ol_days + seq(0, fu_days, by = 28), i))))
      mrgsim_d(MOD, ev2, idata = idata[idata$ID %in% ids, , drop = FALSE],
               end = ol_days + fu_days, delta = 28) %>% as_tibble() %>%
        filter(time >= ol_days)
    }
    a <- run2(cont, TRUE) %>% mutate(arm = "escitalopram continued")
    b <- run2(pbo, FALSE) %>% mutate(arm = "placebo")
    bind_rows(a, b) %>%
      mutate(HAMA_obs = HAMA + rnorm(n(), 0, 2.2))
  })

  output$p_relapse <- renderPlot({
    d <- relapse()
    req(d)
    d %>% group_by(arm, ID) %>% arrange(time) %>%
      mutate(rel = cummax(as.integer(HAMA_obs >= 15))) %>%
      group_by(arm, time) %>%
      summarise(free = 100 * (1 - mean(rel)), .groups = "drop") %>%
      ggplot(aes((time - min(time)) / 7, free, colour = arm)) +
      geom_step(linewidth = 1.2) +
      scale_colour_manual(values = PAL[c(1, 2)]) + ylim(0, 100) +
      labs(x = "weeks after randomisation", y = "% relapse-free", colour = NULL,
           title = "Randomised withdrawal",
           subtitle = "reported: 19% relapse on escitalopram, 56% on placebo") +
      theme_gad()
  })

  output$t_relapse <- renderDT({
    d <- relapse()
    req(d)
    d %>% group_by(arm, ID) %>%
      summarise(relapsed = any(HAMA_obs >= 15), .groups = "drop") %>%
      group_by(arm) %>%
      summarise(n = n(), `relapse %` = round(100 * mean(relapsed), 1),
                .groups = "drop") %>%
      datatable(options = list(dom = "t"), rownames = FALSE)
  })
}

shinyApp(ui, server)
