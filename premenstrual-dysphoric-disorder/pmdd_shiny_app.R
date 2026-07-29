##  Premenstrual Dysphoric Disorder (PMDD) — Shiny dashboard for the QSP model
##  ============================================================================
##  Run with:
##      library(shiny); runApp("pmdd_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  and pmdd_mrgsolve_model.R in the same directory.
##
##  The app is organised around the model's central claim: PMDD is a GAIN
##  disorder, so the two controls that matter most are KP (where the peak of the
##  inverted-U neurosteroid response sits) and the two plasticity gains — NOT
##  any hormone level. Tab 4 draws the inverted U with the patient's own luteal
##  trajectory laid on top of it, and tab 10 runs the five structural
##  predictions side by side with their falsified (level-detector) counterparts.
##  ============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

MODEL_FILE <- "pmdd_mrgsolve_model.R"
TCYC <- 28

mod <- mread_cache("pmdd", MODEL_FILE)

## ---------------------------------------------------------------- helpers ---
cyclic_times <- function(tend, every = 1, lo = 1, hi = 29, tcyc = TCYC) {
  t <- seq(0, tend - every, by = every)
  cd <- (t %% tcyc) + 1
  t[cd >= lo & cd < hi]
}

dose_df <- function(cmt, amt, times) {
  if (!length(times)) return(NULL)
  data.frame(ID = 1, time = times, cmt = cmt, amt = amt, evid = 1)
}

## Every treatment the app can switch on, expressed as (param, dosing) pairs.
REGIMENS <- list(
  "None (untreated)"                          = "none",
  "Sertraline 50 mg — luteal only (cd 15-28)" = "sert_lut",
  "Sertraline 50 mg — continuous"             = "sert_cont",
  "Sertraline 50 mg — symptom onset (cd 22-)" = "sert_onset",
  "Fluoxetine 20 mg — continuous"             = "flx",
  "Drospirenone/EE — 24/4"                    = "drsp244",
  "Drospirenone/EE — 21/7"                    = "drsp217",
  "Drospirenone/EE — continuous"              = "drspcont",
  "Leuprolide 3.75 mg q28d"                   = "leup",
  "Leuprolide + E2/P4 add-back (from d112)"   = "addback",
  "Sepranolone 10 mg SC q48h (luteal)"        = "sepra",
  "Dutasteride 2.5 mg/d"                      = "duta25",
  "Dutasteride 0.5 mg/d"                      = "duta05",
  "Neurosteroid analogue — sedative dose"     = "allohigh",
  "Neurosteroid analogue — sub-sedative dose" = "allolow",
  "Alprazolam 0.25 mg TID (luteal)"           = "alpz",
  "Ulipristal acetate 5 mg/d"                 = "upa",
  "Calcium carbonate 1200 mg/d"               = "calcium",
  "Spironolactone 100 mg/d (luteal)"          = "spiro",
  "CBT / mindfulness"                         = "cbt",
  "Sertraline (luteal) + drospirenone 24/4"   = "combo",
  "Bilateral oophorectomy + transdermal E2"   = "ooph"
)

build_run <- function(key, tend) {
  p <- list(); d <- NULL
  switch(key,
    none       = {},
    sert_lut   = d <- dose_df("SERT_D", 50, cyclic_times(tend, 1, 15, 29)),
    sert_cont  = d <- dose_df("SERT_D", 50, seq(0, tend - 1)),
    sert_onset = d <- dose_df("SERT_D", 50, cyclic_times(tend, 1, 22, 29)),
    flx        = d <- dose_df("FLX_C", 20 * 1000 * 0.9 / 2500, seq(0, tend - 1)),
    drsp244    = d <- dose_df("DRSP_D", 3, cyclic_times(tend, 1, 1, 25)),
    drsp217    = d <- dose_df("DRSP_D", 3, cyclic_times(tend, 1, 1, 22)),
    drspcont   = d <- dose_df("DRSP_D", 3, seq(0, tend - 1)),
    leup       = d <- dose_df("LEUP_D", 3.75, seq(0, tend - 1, by = 28)),
    addback    = {
      p <- list(ADDBACK_E2 = 210, ADDBACK_P4 = 20, ADDBACK_T0 = 112)
      d <- dose_df("LEUP_D", 3.75, seq(0, tend - 1, by = 28))
    },
    sepra      = d <- dose_df("SEPRA_D", 10, cyclic_times(tend, 2, 15, 29)),
    duta25     = d <- dose_df("DUTA_D", 2.5, seq(0, tend - 1)),
    duta05     = d <- dose_df("DUTA_D", 0.5, seq(0, tend - 1)),
    allohigh   = d <- dose_df("ZUR", 20, cyclic_times(tend, 1, 22, 29)),
    allolow    = d <- dose_df("ZUR", 3, cyclic_times(tend, 1, 22, 29)),
    alpz       = d <- dose_df("ALPZ_D", 0.25, cyclic_times(tend, 1/3, 15, 29)),
    upa        = d <- dose_df("UPA_C", 5 * 1000 * 0.8 / 500, seq(0, tend - 1)),
    calcium    = p <- list(CA_SUPP = 0.18),
    spiro      = p <- list(SPIRO = 1.6),
    cbt        = p <- list(CBT = 0.18),
    combo      = d <- rbind(dose_df("SERT_D", 50, cyclic_times(tend, 1, 15, 29)),
                            dose_df("DRSP_D", 3, cyclic_times(tend, 1, 1, 25))),
    ooph       = p <- list(OOPH = 1, ADDBACK_E2 = 170)
  )
  list(param = p, dose = d)
}

run_sim <- function(key, phen, tend, extra = list()) {
  r <- build_run(key, tend)
  m <- param(mod, c(r$param, phen, extra))
  s <- if (is.null(r$dose)) mrgsim(m, end = tend, delta = 0.25)
       else mrgsim(m, data = r$dose, end = tend, delta = 0.25)
  as.data.frame(s)
}

## Trial-style endpoints on one cycle -----------------------------------------
endpoint <- function(d, cyc, tcyc = TCYC) {
  w <- d %>% filter(time >= (cyc - 1) * tcyc, time < cyc * tcyc)
  fol <- w %>% filter(cycle_day >= 5, cycle_day < 11)
  lut <- w %>% filter(cycle_day >= 21, cycle_day < 29)
  list(
    fol       = mean(fol$DRSP_11),
    lut       = mean(lut$DRSP_11),
    peak      = max(lut$DRSP_11),
    peak_day  = lut$cycle_day[which.max(lut$DRSP_11)],
    core_fol  = mean(fol$DRSP_CORE),
    core_lut  = mean(lut$DRSP_CORE),
    core_inc  = 100 * (mean(lut$DRSP_CORE) - mean(fol$DRSP_CORE)) / mean(fol$DRSP_CORE),
    allo_peak = max(w$ALLOB),
    allo_day  = w$cycle_day[which.max(w$ALLOB)],
    hf        = mean(w$HF),
    bmd       = tail(w$BMD, 1)
  )
}

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

## ------------------------------------------------------------------- UI -----
ui <- fluidPage(
  titlePanel("Premenstrual Dysphoric Disorder — QSP model explorer"),
  p(HTML("A 60-compartment QSP model in which PMDD is a <b>gain</b> disorder: ",
         "the neurosteroid response is <b>non-monotonic</b> (inverted-U, peak at ",
         "<code>KP</code>) and its gain is set by a <b>slow, symmetric change ",
         "detector</b> in GABA-A subunit composition. Hormone levels are ",
         "identical in cases and controls — only <code>KP</code> and two gains ",
         "differ.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Sensitivity phenotype"),
      sliderInput("KP", "KP — neurosteroid load at the response peak (nM)",
                  min = 1.5, max = 16, value = 5.0, step = 0.5),
      sliderInput("GAMMA_A4", "alpha4 gain (withdrawal sensitisation)",
                  min = 0, max = 2.5, value = 1.5, step = 0.05),
      sliderInput("GAMMA_D", "delta gain (negative = tolerance)",
                  min = -1, max = 0.5, value = -0.35, step = 0.05),
      actionButton("preset_pmdd", "PMDD preset"),
      actionButton("preset_ctrl", "Control preset"),
      hr(),
      h4("Physiology"),
      sliderInput("TCYC", "cycle length (d)", 24, 34, 28, step = 1),
      sliderInput("STRESS", "background stress load", 0, 1, 0.30, step = 0.05),
      sliderInput("K5AR", "SRD5A1 activity (x nominal)", 0.4, 2.0, 1.0, step = 0.05),
      hr(),
      h4("Treatment"),
      selectInput("regA", "Regimen A", choices = REGIMENS, selected = "none"),
      selectInput("regB", "Regimen B (comparator)", choices = REGIMENS,
                  selected = "sert_lut"),
      sliderInput("ncyc", "cycles simulated", 4, 12, 6, step = 1),
      sliderInput("evalcyc", "cycle to report", 2, 12, 5, step = 1),
      hr(),
      checkboxInput("falsify",
                    "Falsification mode (level detector, no plasticity)", FALSE),
      helpText("Falsification mode sets MONOTONE = 1, RATE_OFF = 1 and both",
               "plasticity gains to 0. Every prediction in tab 10 should",
               "collapse or invert.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient profile", br(),
                 fluidRow(column(6, DTOutput("tbl_profile")),
                          column(6, plotOutput("plt_profile", height = 320))),
                 br(), verbatimTextOutput("txt_profile")),
        tabPanel("2 · Cycle engine", br(),
                 plotOutput("plt_cycle", height = 560)),
        tabPanel("3 · Neurosteroids", br(),
                 plotOutput("plt_ns", height = 520),
                 helpText("Plasma and brain allopregnanolone, 5-alpha-DHP and",
                          "SRD5A1 activity. Dutasteride acts on SRD5A1;",
                          "SSRIs raise ALLO through the 3-alpha-HSD shift.")),
        tabPanel("4 · The inverted U", br(),
                 plotOutput("plt_invu", height = 520),
                 helpText("The unimodal transduction curve for the current",
                          "phenotype, with the simulated luteal trajectory laid",
                          "over it. If the operating range straddles the peak,",
                          "raising and lowering the load can BOTH help.")),
        tabPanel("5 · Receptor plasticity", br(),
                 plotOutput("plt_plast", height = 520)),
        tabPanel("6 · Drug PK", br(),
                 plotOutput("plt_pk", height = 520)),
        tabPanel("7 · DRSP domains", br(),
                 plotOutput("plt_domains", height = 620)),
        tabPanel("8 · Clinical endpoints", br(),
                 DTOutput("tbl_endpoints"), br(),
                 plotOutput("plt_endpoint", height = 380)),
        tabPanel("9 · Scenario comparison", br(),
                 plotOutput("plt_compare", height = 420),
                 br(), DTOutput("tbl_compare")),
        tabPanel("10 · Predictions [P1]-[P5]", br(),
                 DTOutput("tbl_pred"), br(),
                 plotOutput("plt_pred", height = 380),
                 helpText("Each row is computed live from the model. Switch on",
                          "falsification mode in the sidebar and watch the",
                          "signs change.")),
        tabPanel("11 · Safety / organ systems", br(),
                 plotOutput("plt_safety", height = 520),
                 helpText("Hypoestrogenic cost of ovarian suppression: hot",
                          "flushes and lumbar BMD. Add-back protects the",
                          "skeleton; oophorectomy plus estradiol does too.")),
        tabPanel("12 · Model notes", br(), htmlOutput("txt_notes"))
      )
    )
  )
)

## ---------------------------------------------------------------- server ----
server <- function(input, output, session) {

  observeEvent(input$preset_pmdd, {
    updateSliderInput(session, "KP", value = 5.0)
    updateSliderInput(session, "GAMMA_A4", value = 1.5)
    updateSliderInput(session, "GAMMA_D", value = -0.35)
  })
  observeEvent(input$preset_ctrl, {
    updateSliderInput(session, "KP", value = 13.0)
    updateSliderInput(session, "GAMMA_A4", value = 0.25)
    updateSliderInput(session, "GAMMA_D", value = -0.05)
  })

  phen <- reactive({
    p <- list(KP = input$KP, GAMMA_A4 = input$GAMMA_A4, GAMMA_D = input$GAMMA_D,
              TCYC = input$TCYC, STRESS = input$STRESS,
              K5AR = 3.0 * input$K5AR)
    if (isTRUE(input$falsify))
      p <- c(p, list(MONOTONE = 1, RATE_OFF = 1, GAMMA_D = 0, GAMMA_A4 = 0))
    p
  })

  tend <- reactive(input$ncyc * input$TCYC)
  cyc  <- reactive(min(input$evalcyc, input$ncyc))

  simA <- reactive(run_sim(input$regA, phen(), tend()))
  simB <- reactive(run_sim(input$regB, phen(), tend()))
  simRef <- reactive(run_sim("none", phen(), tend()))

  ## one evaluation cycle, both arms, long format
  cyc_long <- reactive({
    lo <- (cyc() - 1) * input$TCYC; hi <- cyc() * input$TCYC
    bind_rows(
      simA() %>% filter(time >= lo, time < hi) %>% mutate(arm = "A"),
      simB() %>% filter(time >= lo, time < hi) %>% mutate(arm = "B")
    )
  })

  ## ---- 1. profile ----------------------------------------------------------
  output$tbl_profile <- renderDT({
    e <- endpoint(simRef(), cyc(), input$TCYC)
    data.frame(
      Quantity = c("KP (nM)", "alpha4 gain", "delta gain", "cycle length (d)",
                   "E2 peak (pg/mL)", "P4 peak (ng/mL)",
                   "brain ALLO follicular (nM)", "brain ALLO luteal peak (nM)",
                   "ALLO peak (cycle day)", "DRSP-11 follicular",
                   "DRSP-11 luteal peak", "symptom peak (cycle day)",
                   "core luteal increase (%)", "meets DSM-5 >=30% criterion"),
      Value = c(sprintf("%.1f", input$KP), sprintf("%.2f", input$GAMMA_A4),
                sprintf("%.2f", input$GAMMA_D), input$TCYC,
                sprintf("%.0f", max(simRef()$E2)),
                sprintf("%.1f", max(simRef()$P4)),
                sprintf("%.2f", min(simRef()$ALLOB)),
                sprintf("%.2f", e$allo_peak), sprintf("%.1f", e$allo_day),
                sprintf("%.1f", e$fol), sprintf("%.1f", e$peak),
                sprintf("%.1f", e$peak_day), sprintf("%+.0f", e$core_inc),
                ifelse(e$core_inc >= 30, "YES", "no")),
      check.names = FALSE)
  }, options = list(dom = "t", pageLength = 20), rownames = FALSE)

  output$plt_profile <- renderPlot({
    lo <- (cyc() - 1) * input$TCYC; hi <- cyc() * input$TCYC
    d <- simRef() %>% filter(time >= lo, time < hi)
    ggplot(d, aes(cycle_day)) +
      geom_line(aes(y = DRSP_11, colour = "DRSP-11"), linewidth = 1) +
      geom_line(aes(y = 11 + 55 * ALLOB / max(ALLOB), colour = "brain ALLO (scaled)"),
                linewidth = 0.8, linetype = 2) +
      scale_colour_manual(values = c("DRSP-11" = "#c92a2a",
                                     "brain ALLO (scaled)" = "#1c7ed6")) +
      labs(x = "cycle day", y = "DRSP-11 (11-66)", colour = NULL,
           title = "Untreated cycle for this phenotype",
           subtitle = "the symptom peak should sit days AFTER the ALLO peak") +
      THEME
  })

  output$txt_profile <- renderText({
    e <- endpoint(simRef(), cyc(), input$TCYC)
    sprintf(paste("Symptom peak lags the allopregnanolone peak by %.1f days.",
                  "\nCore affective luteal increase %+.0f%% (DSM-5 needs >=30%%).",
                  "\nMode: %s"),
            e$peak_day - e$allo_day, e$core_inc,
            ifelse(isTRUE(input$falsify),
                   "FALSIFICATION (level detector)", "base model"))
  })

  ## ---- 2. cycle engine -----------------------------------------------------
  output$plt_cycle <- renderPlot({
    lo <- (cyc() - 1) * input$TCYC; hi <- cyc() * input$TCYC
    d <- cyc_long() %>%
      select(cycle_day, arm, FSH, LH, FOLL, CL, E2, P4) %>%
      pivot_longer(-c(cycle_day, arm))
    d$name <- factor(d$name, levels = c("FSH", "LH", "FOLL", "CL", "E2", "P4"))
    ggplot(d, aes(cycle_day, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "cycle day", y = NULL, colour = "arm",
           title = "HPO cycle engine",
           subtitle = "A = first regimen, B = comparator") + THEME
  })

  ## ---- 3. neurosteroids ----------------------------------------------------
  output$plt_ns <- renderPlot({
    d <- cyc_long() %>%
      select(cycle_day, arm, ALLOP, ALLOB, DHP, SRD5A1_activity, NS_load, ISOB) %>%
      pivot_longer(-c(cycle_day, arm))
    ggplot(d, aes(cycle_day, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "cycle day", y = NULL, colour = "arm",
           title = "Neurosteroidogenesis and effective load") + THEME
  })

  ## ---- 4. the inverted U ---------------------------------------------------
  output$plt_invu <- renderPlot({
    L <- seq(0, 20, by = 0.05)
    kp <- input$KP
    curve <- data.frame(L = L, drive = (L / kp) * exp(1 - L / kp))
    d <- cyc_long() %>% filter(arm == "A")
    traj <- data.frame(L = d$NS_load,
                       drive = (d$NS_load / kp) * exp(1 - d$NS_load / kp),
                       cycle_day = d$cycle_day)
    ggplot(curve, aes(L, drive)) +
      geom_line(linewidth = 1.1, colour = "#495057") +
      geom_vline(xintercept = kp, linetype = 3) +
      geom_point(data = traj, aes(colour = cycle_day), size = 2, alpha = 0.75) +
      scale_colour_viridis_c(option = "plasma") +
      annotate("text", x = kp, y = 1.05, hjust = -0.1, size = 3.4,
               label = sprintf("KP = %.1f nM (flat top: dDrive/dL = 0)", kp)) +
      labs(x = "effective neurosteroid load L (nM)",
           y = "position term of NS_DRIVE",
           colour = "cycle day",
           title = "Non-monotonic transduction with the patient's trajectory",
           subtitle = paste("ascending limb: more ALLO is worse ·",
                            "descending limb: more ALLO is better")) + THEME
  })

  ## ---- 5. plasticity -------------------------------------------------------
  output$plt_plast <- renderPlot({
    d <- cyc_long() %>%
      select(cycle_day, arm, DELTA, ALPHA4, NS_sensitivity, Change_signal,
             NS_drive, Affective_effect) %>%
      pivot_longer(-c(cycle_day, arm))
    ggplot(d, aes(cycle_day, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "cycle day", y = NULL, colour = "arm",
           title = "GABA-A composition, the change signal, and the gain",
           subtitle = paste("alpha4 follows the RATE of change; delta follows",
                            "the exposure history and is tolerance-like")) + THEME
  })

  ## ---- 6. PK ---------------------------------------------------------------
  output$plt_pk <- renderPlot({
    d <- bind_rows(simA() %>% mutate(arm = "A"), simB() %>% mutate(arm = "B")) %>%
      select(time, arm, SERT_C, DMS_C, FLX_C, NFLX_C, DRSP_C, LEUP_C,
             SEPRA_C, DUTA_C, ALPZ_C, UPA_C, ZUR) %>%
      pivot_longer(-c(time, arm)) %>%
      group_by(name) %>% filter(max(value) > 1e-6) %>% ungroup()
    if (!nrow(d))
      return(ggplot() + annotate("text", 0, 0, label = "no drug in either arm") +
               theme_void())
    ggplot(d, aes(time, value, colour = arm)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "time (d)", y = "concentration (ng/mL; ZUR in nM)", colour = "arm",
           title = "Drug pharmacokinetics") + THEME
  })

  ## ---- 7. DRSP domains ----------------------------------------------------
  output$plt_domains <- renderPlot({
    d <- cyc_long() %>%
      select(cycle_day, arm, S_IRR, S_DEP, S_ANX, S_LAB, S_ANH, S_COG,
             S_FAT, S_APP, S_SLP, S_OVR, S_PHY) %>%
      pivot_longer(-c(cycle_day, arm))
    ggplot(d, aes(cycle_day, value, colour = arm)) +
      geom_line(linewidth = 0.85) +
      facet_wrap(~name, ncol = 3) + ylim(1, 6) +
      labs(x = "cycle day", y = "domain score (1-6)", colour = "arm",
           title = "The eleven DRSP domains",
           subtitle = paste("core affective domains: irritability, depressed",
                            "mood, anxiety, lability, anhedonia")) + THEME
  })

  ## ---- 8. endpoints -------------------------------------------------------
  output$tbl_endpoints <- renderDT({
    eA <- endpoint(simA(), cyc(), input$TCYC)
    eB <- endpoint(simB(), cyc(), input$TCYC)
    eR <- endpoint(simRef(), cyc(), input$TCYC)
    red <- function(e) 100 * (1 - (e$lut - 11) / (eR$lut - 11))
    data.frame(
      Endpoint = c("DRSP-11 follicular", "DRSP-11 luteal mean",
                   "DRSP-11 luteal peak", "symptom peak (cycle day)",
                   "core luteal increase (%)",
                   "luteal burden reduction vs untreated (%)",
                   "hot flushes / day", "lumbar BMD change (%)"),
      A = c(sprintf("%.1f", eA$fol), sprintf("%.1f", eA$lut),
            sprintf("%.1f", eA$peak), sprintf("%.1f", eA$peak_day),
            sprintf("%+.0f", eA$core_inc), sprintf("%+.1f", red(eA)),
            sprintf("%.2f", eA$hf), sprintf("%+.2f", eA$bmd)),
      B = c(sprintf("%.1f", eB$fol), sprintf("%.1f", eB$lut),
            sprintf("%.1f", eB$peak), sprintf("%.1f", eB$peak_day),
            sprintf("%+.0f", eB$core_inc), sprintf("%+.1f", red(eB)),
            sprintf("%.2f", eB$hf), sprintf("%+.2f", eB$bmd)),
      check.names = FALSE)
  }, options = list(dom = "t"), rownames = FALSE)

  output$plt_endpoint <- renderPlot({
    d <- cyc_long()
    ggplot(d, aes(cycle_day, DRSP_11, colour = arm)) +
      geom_line(linewidth = 1) +
      geom_ribbon(data = data.frame(x = c(21, 28)),
                  aes(x = x, ymin = -Inf, ymax = Inf), inherit.aes = FALSE,
                  fill = "#f1f3f5", alpha = 0.5) +
      geom_line(linewidth = 1) +
      labs(x = "cycle day", y = "DRSP-11 (11-66)", colour = "arm",
           title = "Daily symptom trajectory over the reported cycle",
           subtitle = "shaded band = the luteal window used for the endpoint") +
      THEME
  })

  ## ---- 9. scenario comparison ---------------------------------------------
  COMPARE_SET <- c("none", "sert_lut", "sert_cont", "drsp244", "drsp217",
                   "leup", "sepra", "duta25", "duta05", "allohigh", "allolow",
                   "combo")
  compare <- reactive({
    lo <- (cyc() - 1) * input$TCYC; hi <- cyc() * input$TCYC
    ref <- endpoint(simRef(), cyc(), input$TCYC)
    res <- lapply(COMPARE_SET, function(k) {
      s <- run_sim(k, phen(), tend())
      e <- endpoint(s, cyc(), input$TCYC)
      data.frame(regimen = names(REGIMENS)[match(k, unlist(REGIMENS))],
                 key = k, luteal = e$lut, peak = e$peak,
                 reduction = 100 * (1 - (e$lut - 11) / (ref$lut - 11)),
                 peak_day = e$peak_day)
    })
    bind_rows(res)
  })

  output$plt_compare <- renderPlot({
    d <- compare() %>% filter(key != "none") %>%
      mutate(regimen = reorder(regimen, reduction))
    ggplot(d, aes(reduction, regimen,
                  fill = reduction > 0)) +
      geom_col(width = 0.7) +
      scale_fill_manual(values = c("TRUE" = "#2f9e44", "FALSE" = "#c92a2a"),
                        guide = "none") +
      labs(x = "% reduction of the luteal DRSP burden above floor", y = NULL,
           title = "Every regimen, one endpoint",
           subtitle = "negative = symptoms made worse") + THEME
  })

  output$tbl_compare <- renderDT({
    compare() %>% select(-key) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
  }, options = list(dom = "tp", pageLength = 12), rownames = FALSE)

  ## ---- 10. the five predictions -------------------------------------------
  output$tbl_pred <- renderDT({
    ph <- phen()
    ref <- endpoint(run_sim("none", ph, tend()), cyc(), input$TCYC)
    red <- function(k, t = tend(), c = cyc()) {
      e <- endpoint(run_sim(k, ph, t), c, input$TCYC)
      100 * (1 - (e$lut - 11) / (ref$lut - 11))
    }
    lag <- ref$peak_day - ref$allo_day
    d25 <- red("duta25", 10 * input$TCYC, min(9, 10))
    d05 <- red("duta05", 10 * input$TCYC, min(9, 10))
    ah  <- red("allohigh"); al <- red("allolow")
    r244 <- red("drsp244"); r217 <- red("drsp217")
    data.frame(
      Prediction = c("[P1] symptom peak lags the ALLO peak",
                     "[P2] dutasteride 2.5 mg vs 0.5 mg",
                     "[P3] neurosteroid: sedative vs sub-sedative dose",
                     "[P4] 24/4 vs 21/7 hormone-free interval",
                     "[P5] see the add-back regimen on tab 8"),
      Model = c(sprintf("%+.1f d", lag),
                sprintf("%+.1f%% vs %+.1f%%", d25, d05),
                sprintf("%+.1f%% vs %+.1f%%", ah, al),
                sprintf("%+.1f%% vs %+.1f%% (gap %.1f pt)", r244, r217, r244 - r217),
                "flare then remission under constant hormones"),
      Expected = c("3-6 d (level detector: ~1 d)",
                   "2.5 mg works, 0.5 mg does not",
                   "opposite SIGNS from one drug",
                   "24/4 better by >=8 points",
                   "transient, not sustained"),
      check.names = FALSE)
  }, options = list(dom = "t"), rownames = FALSE)

  output$plt_pred <- renderPlot({
    ph <- phen()
    hi <- run_sim("allohigh", ph, tend()); lo <- run_sim("allolow", ph, tend())
    nn <- run_sim("none", ph, tend())
    a <- (cyc() - 1) * input$TCYC; b <- cyc() * input$TCYC
    d <- bind_rows(
      nn %>% filter(time >= a, time < b) %>% mutate(arm = "untreated"),
      hi %>% filter(time >= a, time < b) %>% mutate(arm = "sedative dose"),
      lo %>% filter(time >= a, time < b) %>% mutate(arm = "sub-sedative dose"))
    ggplot(d, aes(cycle_day, DRSP_11, colour = arm)) +
      geom_line(linewidth = 1) +
      labs(x = "cycle day", y = "DRSP-11", colour = NULL,
           title = "[P3] one drug, two doses, opposite directions",
           subtitle = paste("both doses of the same neurosteroid analogue,",
                            "given on cycle days 22-28")) + THEME
  })

  ## ---- 11. safety ---------------------------------------------------------
  output$plt_safety <- renderPlot({
    d <- bind_rows(simA() %>% mutate(arm = "A"), simB() %>% mutate(arm = "B")) %>%
      select(time, arm, HF, BMD, E2_total, Hypoestrogenism, ENDO) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "time (d)", y = NULL, colour = "arm",
           title = "Organ-system and safety readouts") + THEME
  })

  ## ---- 12. notes ----------------------------------------------------------
  output$txt_notes <- renderUI({
    HTML(paste0(
      "<h4>What this model commits to</h4>",
      "<p>PMDD is modelled as a <b>gain</b> disorder. Luteal E2, P4 and ",
      "allopregnanolone are identical in the PMDD and control presets — only ",
      "<code>KP</code> (where the inverted-U response peaks) and the two ",
      "plasticity gains differ. Set the control preset and watch the same ",
      "hormones stop producing symptoms.</p>",
      "<h4>Structure</h4><ul>",
      "<li>Non-monotonic transduction: <code>NS_DRIVE = SENS · (L/KP) · ",
      "exp(1 − L/KP)</code>, peaking at <code>L = KP</code>.</li>",
      "<li>A slow, symmetric change detector: delta follows the exposure ",
      "history (tolerance), alpha4 follows the rate of change of neurosteroid ",
      "and progestogen tone (sensitisation).</li>",
      "<li>60 ODE compartments; days as the time unit; DRSP-11 as the ",
      "endpoint (11 = no symptoms, 66 = maximum).</li></ul>",
      "<h4>Two fitted parameters</h4>",
      "<p><code>KI_ISO</code> was fitted to the sepranolone effect size and ",
      "<code>W_SHT_AMY</code> to the sertraline effect size. Everything else ",
      "comes from published PK or from the physiological ranges asserted in ",
      "the twin's check table.</p>",
      "<h4>Files</h4><ul>",
      "<li><code>pmdd_mrgsolve_model.R</code> — the model</li>",
      "<li><code>pmdd_python_twin.py</code> — a dependency-free numerical twin ",
      "that asserts 38 checks</li>",
      "<li><code>pmdd_qsp_model.dot/.svg/.png</code> — the mechanistic map</li>",
      "<li><code>pmdd_references.md</code> — 59 PubMed-verified references</li>",
      "</ul>",
      "<p><b>Educational and research use only.</b> Not validated for clinical ",
      "decision-making.</p>"))
  })
}

shinyApp(ui, server)
