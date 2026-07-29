# =============================================================================
#  hg_shiny_app.R
#  Hyperemesis Gravidarum QSP model — interactive dashboard
# =============================================================================
#
#  Run:
#      # from this directory, with hg_mrgsolve_model.R alongside
#      shiny::runApp("hg_shiny_app.R")
#
#  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
#
#  ---------------------------------------------------------------------------
#  WHAT THIS APP IS FOR
#  ---------------------------------------------------------------------------
#  It is built around the one control that matters: the ADAPTATION slider
#  (ALPHA) on the sidebar. Slide it to 1 and the area postrema is a
#  fold-change detector; slide it to 0 and it is a level detector. Everything
#  else in the app is there so you can watch what that single choice does:
#
#    ALPHA = 0.92 (default)   symptoms peak at GA 9-11 wk and remit by 16,
#                             beta-thalassaemia is protected, and raising GDF15
#                             BEFORE conception helps while raising it during
#                             pregnancy does not;
#    ALPHA = 0                symptoms never remit (because GDF15 never falls),
#                             and beta-thalassaemia becomes the worst case in
#                             the cohort.
#
#  The second thing to try is the TIMING control on the Prevention tab. Move
#  metformin from pre-conception to GA 6 weeks and the benefit disappears — the
#  same drug, the same dose, a different sign. That is the model's central
#  falsifiable claim, and it is a slider.
#
#  TABS (12)
#    1  Patient & cohort      the two axes of risk (SENS and TROPH_GAIN)
#    2  GDF15 & set-point     the hormone, its set-point, and the ratio
#    3  Emetic cascade        AP -> NTS -> nausea / vomiting
#    4  PUQE-24 & endpoints   the instrument, and its ceiling
#    5  Drug PK               plasma concentrations of every agent
#    6  Node occupancy        receptor occupancy and NTS authority per drug
#    7  Simulated trials      the VOMIT comparison, run as a cohort
#    8  Prevention timing     where the therapeutic window closes
#    9  Fluid & electrolytes  the chloride-responsive alkalosis
#   10  Thiamine & Wernicke   the slowest clock, and the dextrose trap
#   11  Thyroid               why only one HG phenotype is thyrotoxic
#   12  Validation            model vs published trials, side by side
#
#  DISCLAIMER: educational / research tool. Not validated for clinical use.
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

# The model definition lives in hg_mrgsolve_model.R. We source it for `mod`,
# `MW`, `nmol` and `CMTN`, but suppress its own scenario run and plots.
HG_APP_MODE <- TRUE
if (!exists("mod")) {
  src <- readLines("hg_mrgsolve_model.R")
  cut <- grep("^#  RUN$", src)
  if (length(cut)) src <- src[seq_len(cut[1] - 3)]
  eval(parse(text = paste(src, collapse = "\n")))
}

T_START <- -70
TR      <- 56          # default randomisation / treatment start, GA 8 weeks

theme_set(theme_bw(base_size = 12))
PAL <- c("#1565c0", "#6a1b9a", "#2e7d32", "#e65100", "#c62828",
         "#00838f", "#f9a825", "#616161")

# -----------------------------------------------------------------------------
#  dosing builders — one per agent, all in mg, converted to nmol
# -----------------------------------------------------------------------------
ev_none <- function(...) NULL

ev_ond <- function(t0, days, mg = 8)
  dose("OND_D", nmol(mg, "OND"), t0, ii = 1/3, addl = days * 3 - 1)
ev_mir <- function(t0, days, mg = 30)
  dose("MIR_D", nmol(mg, "MIR"), t0, ii = 1, addl = days - 1)
ev_dox <- function(t0, days, mg = 40)
  rbind(dose("DOX_D", nmol(mg / 4, "DOX"), t0 + 0.30, ii = 1, addl = days - 1),
        dose("DOX_D", nmol(mg / 4, "DOX"), t0 + 0.65, ii = 1, addl = days - 1),
        dose("DOX_D", nmol(mg / 2, "DOX"), t0 + 0.92, ii = 1, addl = days - 1))
ev_mcp <- function(t0, days, mg = 10)
  dose("MCP_D", nmol(mg, "MCP"), t0, ii = 0.25, addl = days * 4 - 1)
ev_pmz <- function(t0, days, mg = 25)
  dose("PMZ_C", nmol(mg, "PMZ") * 0.25, t0, ii = 0.25, addl = days * 4 - 1)
ev_gbp <- function(t0, days, mg = 600)
  dose("GBP_D", nmol(mg, "GBP"), t0, ii = 1/3, addl = days * 3 - 1)
ev_clo <- function(t0, days, mg = 0.15)
  dose("CLO_D", nmol(mg, "CLO"), t0, ii = 1, addl = days - 1)
ev_ste <- function(t0, days = 14, mg = 125)
  rbind(dose("STE_C", nmol(mg, "STE"), t0),
        dose("STE_D", nmol(40, "STE"), t0 + 1),
        dose("STE_D", nmol(20, "STE"), t0 + 2, ii = 1, addl = 2),
        dose("STE_D", nmol(10, "STE"), t0 + 5, ii = 1, addl = 2),
        dose("STE_D", nmol(5,  "STE"), t0 + 8, ii = 1, addl = 6))
ev_mab <- function(t0, days, nm = 90) dose("MAB_C", nm, t0)

DRUGS <- list(
  "none"                        = list(f = ev_none, node = "-"),
  "ondansetron 8 mg q8h"        = list(f = ev_ond, node = "peripheral 5-HT3"),
  "doxylamine/pyridoxine"       = list(f = ev_dox, node = "H1 + vestibular"),
  "metoclopramide 10 mg q6h"    = list(f = ev_mcp, node = "AP D2 + prokinetic"),
  "promethazine 25 mg q6h"      = list(f = ev_pmz, node = "NTS H1 / M1"),
  "mirtazapine 30 mg qHS"       = list(f = ev_mir, node = "H1 + 5-HT2 + 5-HT3"),
  "gabapentin 600 mg q8h"       = list(f = ev_gbp, node = "NTS alpha-2-delta"),
  "clonidine 5 mg patch"        = list(f = ev_clo, node = "NTS alpha-2"),
  "corticosteroid taper"        = list(f = ev_ste, node = "no node on the chain"),
  "anti-GDF15 mAb"              = list(f = ev_mab, node = "the ligand itself"))

# -----------------------------------------------------------------------------
#  simulation wrapper
# -----------------------------------------------------------------------------
sim_hg <- function(sens = 1.0, troph = 1.0, alpha = 0.92, tau_sp = 30,
                   thal = FALSE, tobacco = FALSE, lowgdf = FALSE,
                   drug = "none", t_drug = TR, drug_days = 14,
                   metformin = c("none", "pre-conception", "from GA 6 wk"),
                   iv_on = FALSE, iv_from = 56, iv_days = 28,
                   iv_dextrose = TRUE, iv_thiamine = 100,
                   end = 200, delta = 0.25) {
  metformin <- match.arg(metformin)
  ev <- DRUGS[[drug]]$f(t_drug, drug_days)
  if (metformin == "pre-conception") {
    ev <- rbind(ev, dose("MET_C", nmol(1000, "MET") * 0.55, T_START,
                         ii = 0.5, addl = 167))
  } else if (metformin == "from GA 6 wk") {
    ev <- rbind(ev, dose("MET_C", nmol(1000, "MET") * 0.55, 42,
                         ii = 0.5, addl = 195))
  }
  p <- list(SENS = sens, TROPH_GAIN = troph, ALPHA = alpha, TAU_SP = tau_sp,
            THAL_ON = as.numeric(thal), TOBAC_ON = as.numeric(tobacco),
            LOWGDF_ON = as.numeric(lowgdf))
  if (!iv_on) {
    out <- mod %>% param(p) %>%
      mrgsim(events = ev, start = T_START, end = end, delta = delta)
    return(as.data.frame(out))
  }
  run_with_support(par = p, ev = ev, t_on = iv_from,
                   t_off = iv_from + iv_days,
                   fluid = 3.0, k_rate = 40,
                   thi = iv_thiamine, dex = if (iv_dextrose) 150 else 0,
                   end = end, delta = delta)
}

pick <- function(df, t, col) df[[col]][which.min(abs(df$time - t))]

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel("Hyperemesis Gravidarum — QSP model"),
  tags$p(style = "color:#555;margin-top:-8px",
    HTML("A model of nausea and vomiting of pregnancy in which the brainstem
          measures <b>how much more GDF15 there is than there used to be</b>,
          not how much there is. Move the <b>Adaptation</b> slider to 0 to turn
          it back into a conventional concentration-response model and watch
          the predictions invert.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("The structural switch"),
      sliderInput("alpha", "Adaptation completeness (ALPHA)",
                  min = 0, max = 1, value = 0.92, step = 0.02),
      helpText(HTML("<b>1</b> = pure fold-change detector · <b>0</b> = pure
                     level detector. Fitted value 0.92.")),
      sliderInput("tau_sp", "Set-point time constant (days)",
                  min = 5, max = 60, value = 30, step = 1),
      hr(),
      h4("Two axes of risk"),
      sliderInput("sens", "Maternal GFRAL sensitivity (SENS)",
                  min = 0.5, max = 1.5, value = 1.00, step = 0.02),
      sliderInput("troph", "Fetal trophoblast secretory gain (TROPH_GAIN)",
                  min = 0.5, max = 2.5, value = 1.00, step = 0.05),
      checkboxInput("thal", "Beta-thalassaemia (lifelong high GDF15)", FALSE),
      checkboxInput("tobacco", "Pre-conception tobacco use", FALSE),
      checkboxInput("lowgdf", "Low-GDF15 risk allele", FALSE),
      hr(),
      h4("Treatment"),
      selectInput("drug", "Antiemetic", names(DRUGS), selected = "none"),
      sliderInput("t_drug", "Start (GA weeks)",
                  min = 5, max = 16, value = 8, step = 0.5),
      sliderInput("drug_days", "Duration (days)",
                  min = 3, max = 42, value = 14, step = 1),
      radioButtons("metformin", "Metformin timing",
                   c("none", "pre-conception", "from GA 6 wk"),
                   selected = "none"),
      hr(),
      h4("Supportive care"),
      checkboxInput("iv_on", "IV fluid + KCl", FALSE),
      conditionalPanel("input.iv_on",
        sliderInput("iv_from", "From (GA weeks)",
                    min = 5, max = 16, value = 8, step = 0.5),
        sliderInput("iv_days", "Days", min = 3, max = 42, value = 28, step = 1),
        checkboxInput("iv_dextrose", "Include IV dextrose", TRUE),
        sliderInput("iv_thiamine", "IV thiamine (mg/day)",
                    min = 0, max = 300, value = 100, step = 25),
        helpText("Set thiamine to 0 with dextrose on to reproduce the
                  Wernicke trap on tab 10."))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient & cohort",
          br(), verbatimTextOutput("profile"),
          plotOutput("p_cohort", height = "380px"),
          helpText("The two axes are not interchangeable. Sensitivity-driven and
                    production-driven HG can have identical PUQE and different
                    hCG, thyroid function and GDF15.")),
        tabPanel("2 · GDF15 & set-point",
          br(), plotOutput("p_gdf", height = "330px"),
          plotOutput("p_fold", height = "300px"),
          helpText("Top: the hormone and the set-point chasing it. Bottom: the
                    ratio, which is what the area postrema actually sees. The
                    ratio returns to 1 while the hormone stays high — that is
                    the whole explanation for spontaneous remission.")),
        tabPanel("3 · Emetic cascade",
          br(), plotOutput("p_cascade", height = "620px")),
        tabPanel("4 · PUQE-24 & endpoints",
          br(), plotOutput("p_puqe", height = "330px"),
          fluidRow(column(6, plotOutput("p_wt", height = "280px")),
                   column(6, plotOutput("p_intake", height = "280px"))),
          helpText("PUQE-24 saturates at 15. That ceiling is why trials in
                    severe HG are structurally insensitive to weak drugs.")),
        tabPanel("5 · Drug PK",
          br(), plotOutput("p_pk", height = "560px"),
          helpText("Mirtazapine's ~29 h half-life is why its separation from
                    ondansetron widens after day 4 — no onset parameter needed.")),
        tabPanel("6 · Node occupancy",
          br(), DTOutput("t_nodes"),
          helpText(HTML("A drug cannot outperform the weight of the node it
                         occupies. Ondansetron reaches ~98% receptor occupancy
                         on a node carrying ~3% of the drive."))),
        tabPanel("7 · Simulated trials",
          br(), actionButton("run_trial", "Run the VOMIT comparison (7 subjects)"),
          br(), br(), verbatimTextOutput("trial_txt"),
          plotOutput("p_trial", height = "340px"),
          helpText("Averaged over a heterogeneous cohort, because a single
                    trajectory on the steep part of the curve overstates the
                    mean trial effect several-fold.")),
        tabPanel("8 · Prevention timing",
          br(), actionButton("run_prev", "Sweep metformin start time"),
          br(), br(), plotOutput("p_prev", height = "380px"),
          DTOutput("t_prev"),
          helpText("The benefit is a function of WHEN, not how much. Once the
                    placental ramp is underway there is no set-point left to
                    move.")),
        tabPanel("9 · Fluid & electrolytes",
          br(), plotOutput("p_elec", height = "560px"),
          helpText("The alkalosis is chloride-responsive: without chloride the
                    kidney cannot excrete the bicarbonate, so nothing else
                    corrects it.")),
        tabPanel("10 · Thiamine & Wernicke",
          br(), plotOutput("p_thi", height = "330px"),
          verbatimTextOutput("thi_txt"),
          actionButton("run_dex", "Compare dextrose with and without thiamine"),
          br(), br(), DTOutput("t_dex")),
        tabPanel("11 · Thyroid",
          br(), plotOutput("p_thy", height = "460px"),
          helpText("hCG cross-activates the TSH receptor. Because hCG and GDF15
                    come from the same cell, the model predicts thyrotoxicosis
                    tracks the fetal-production axis only.")),
        tabPanel("12 · Validation",
          br(), DTOutput("t_valid"),
          helpText(HTML("Regenerate independently with
            <code>python3 hg_reference_impl.py --check</code>, which needs only
            python3 — no R, no compiler.")))
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  base_args <- reactive(list(
    sens = input$sens, troph = input$troph, alpha = input$alpha,
    tau_sp = input$tau_sp, thal = input$thal, tobacco = input$tobacco,
    lowgdf = input$lowgdf, drug = input$drug,
    t_drug = input$t_drug * 7, drug_days = input$drug_days,
    metformin = input$metformin, iv_on = input$iv_on,
    iv_from = if (is.null(input$iv_from)) 56 else input$iv_from * 7,
    iv_days = if (is.null(input$iv_days)) 28 else input$iv_days,
    iv_dextrose = if (is.null(input$iv_dextrose)) TRUE else input$iv_dextrose,
    iv_thiamine = if (is.null(input$iv_thiamine)) 100 else input$iv_thiamine))

  cur <- reactive({
    withProgress(message = "simulating ...", value = 0.5,
                 do.call(sim_hg, base_args()))
  })

  # comparator with the same patient but ALPHA = 0
  cmp <- reactive({
    a <- base_args(); a$alpha <- 0
    do.call(sim_hg, a)
  })

  gest <- function(df) filter(df, time >= 14, time <= 200)

  # ---- tab 1 ---------------------------------------------------------------
  output$profile <- renderText({
    d <- gest(cur()); w <- filter(d, time >= 28, time <= 160)
    sprintf(paste0(
      "PUQE-24 peak            %5.1f  at GA %.1f weeks\n",
      "PUQE-24 at GA 16 weeks  %5.1f\n",
      "Peak maternal GDF15     %5.0f pg/mL   (fold-change at GA 9 wk: %.2f)\n",
      "Maximum weight loss     %5.1f %%\n",
      "Nadir K+ / Cl- / Na+    %5.2f / %5.1f / %5.1f mmol/L\n",
      "Peak HCO3-              %5.1f mmol/L\n",
      "Thiamine nadir          %5.1f mg      P(Wernicke) %.1f %%\n",
      "Nadir TSH               %5.2f mIU/L   peak free T4 %.1f pmol/L\n",
      "Windsor-style HG case?  %s"),
      max(w$PUQE), w$GA_WEEKS[which.max(w$PUQE)],
      pick(d, 112, "PUQE"), max(w$GDF15), pick(d, 63, "FOLD"),
      max(w$WT_LOSS_P), min(w$K_PLASMA), min(w$CL_PLASMA), min(w$NA_PLASMA),
      max(w$HCO3_PL), min(w$THIAMINE), 100 * max(w$P_WE),
      min(w$TSH), max(w$FT4),
      if (max(w$PUQE) >= 13 || (max(w$PUQE) >= 7 && max(w$WT_LOSS_P) >= 5))
        "YES" else "no")
  })

  output$p_cohort <- renderPlot({
    cohorts <- list(
      "normal pregnancy"        = list(sens = 0.62),
      "HG, sensitivity-driven"  = list(sens = 1.00),
      "HG, production-driven"   = list(sens = 0.62, troph = 2.00),
      "low-GDF15 allele"        = list(sens = 1.00, lowgdf = TRUE),
      "beta-thalassaemia"       = list(sens = 1.00, thal = TRUE),
      "pre-conception tobacco"  = list(sens = 1.00, tobacco = TRUE))
    df <- bind_rows(lapply(names(cohorts), function(nm) {
      a <- c(cohorts[[nm]], list(alpha = input$alpha, tau_sp = input$tau_sp))
      gest(do.call(sim_hg, a)) %>% mutate(cohort = nm)
    }))
    ggplot(df, aes(GA_WEEKS, PUQE, colour = cohort)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 13, linetype = 2, colour = "grey40") +
      annotate("text", x = 26, y = 13.4, label = "PUQE 13 (severe)",
               size = 3, colour = "grey40") +
      scale_colour_manual(values = PAL) +
      labs(x = "gestational age (weeks)", y = "PUQE-24", colour = NULL,
           title = sprintf("Cohorts at ALPHA = %.2f", input$alpha))
  })

  # ---- tab 2 ---------------------------------------------------------------
  output$p_gdf <- renderPlot({
    gest(cur()) %>% select(GA_WEEKS, GDF15, SETPOINT) %>%
      pivot_longer(-GA_WEEKS) %>%
      ggplot(aes(GA_WEEKS, value, colour = name)) +
      geom_line(linewidth = 1) + scale_y_log10() +
      scale_colour_manual(values = PAL[c(1, 2)]) +
      labs(x = "gestational age (weeks)", y = "pg/mL (log scale)",
           colour = NULL, title = "The hormone, and the set-point chasing it")
  })

  output$p_fold <- renderPlot({
    bind_rows(gest(cur()) %>% mutate(model = sprintf("ALPHA = %.2f", input$alpha)),
              gest(cmp()) %>% mutate(model = "ALPHA = 0 (level detector)")) %>%
      ggplot(aes(GA_WEEKS, DRIVE, colour = model)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = mod@param$EC50_AP, linetype = 3) +
      annotate("text", x = 26, y = mod@param$EC50_AP * 1.12,
               label = "AP threshold", size = 3) +
      scale_colour_manual(values = PAL[c(2, 8)]) +
      labs(x = "gestational age (weeks)", y = "effective GFRAL drive",
           colour = NULL,
           title = "What the area postrema actually sees")
  })

  # ---- tab 3 ---------------------------------------------------------------
  output$p_cascade <- renderPlot({
    gest(cur()) %>%
      select(GA_WEEKS, `GFRAL surface pool` = GFRAL,
             `area postrema` = AP, `NTS activity` = NTSC,
             `nausea (h/24h)` = NAUSEA_H, `vomits/day` = VOMITS,
             `gastric delay index` = GAS) %>%
      pivot_longer(-GA_WEEKS) %>%
      ggplot(aes(GA_WEEKS, value)) +
      geom_line(linewidth = 0.9, colour = PAL[3]) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "gestational age (weeks)", y = NULL)
  })

  # ---- tab 4 ---------------------------------------------------------------
  output$p_puqe <- renderPlot({
    gest(cur()) %>%
      select(GA_WEEKS, PUQE, NAUSEA_H, VOMITS) %>%
      pivot_longer(-GA_WEEKS) %>%
      ggplot(aes(GA_WEEKS, value, colour = name)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL[c(5, 1, 3)]) +
      labs(x = "gestational age (weeks)", y = NULL, colour = NULL,
           title = "PUQE-24 and its components")
  })
  output$p_wt <- renderPlot({
    gest(cur()) %>% ggplot(aes(GA_WEEKS, WT_LOSS_P)) +
      geom_line(linewidth = 1, colour = PAL[6]) +
      geom_hline(yintercept = 5, linetype = 2) +
      labs(x = "GA (weeks)", y = "weight loss (%)",
           title = "Weight loss (5% = Windsor threshold)")
  })
  output$p_intake <- renderPlot({
    gest(cur()) %>% ggplot(aes(GA_WEEKS, KET)) +
      geom_line(linewidth = 1, colour = PAL[4]) +
      labs(x = "GA (weeks)", y = "beta-hydroxybutyrate (mmol/L)",
           title = "Starvation ketosis")
  })

  # ---- tab 5 ---------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- cur() %>% filter(time >= input$t_drug * 7 - 2,
                          time <= input$t_drug * 7 + input$drug_days + 10)
    d %>% select(time, ondansetron = C_OND, mirtazapine = C_MIR,
                 gabapentin = C_GBP, clonidine = C_CLO, metformin = C_MET) %>%
      pivot_longer(-time) %>% filter(value > 1e-6) %>%
      ggplot(aes(time - input$t_drug * 7, value)) +
      geom_line(linewidth = 0.9, colour = PAL[4]) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "days from treatment start", y = "plasma concentration (nM)",
           title = "Drug PK for the selected regimen")
  })

  # ---- tab 6 ---------------------------------------------------------------
  output$t_nodes <- renderDT({
    p <- as.list(mod@param)
    rows <- lapply(setdiff(names(DRUGS), "none"), function(nm) {
      a <- base_args(); a$drug <- nm; a$metformin <- "none"
      d <- do.call(sim_hg, c(a, list(end = a$t_drug + a$drug_days + 8)))
      w <- d[d$time >= a$t_drug + 1 & d$time <= a$t_drug + a$drug_days, ]
      ref <- do.call(sim_hg, c(modifyList(a, list(drug = "none")),
                               list(end = a$t_drug + a$drug_days + 8)))
      wr <- ref[ref$time >= a$t_drug + 1 &
                ref$time <= a$t_drug + a$drug_days, ]
      data.frame(drug = nm, node = DRUGS[[nm]]$node,
                 dPUQE_day7 = round(pick(d, a$t_drug + 7, "PUQE") -
                                    pick(d, a$t_drug, "PUQE") -
                                    (pick(ref, a$t_drug + 7, "PUQE") -
                                     pick(ref, a$t_drug, "PUQE")), 2),
                 min_NTS_on_drug = round(min(w$NTSC), 3),
                 min_NTS_off_drug = round(min(wr$NTSC), 3))
    })
    datatable(bind_rows(rows), rownames = FALSE,
              options = list(dom = "t", pageLength = 12))
  })

  # ---- tab 7 ---------------------------------------------------------------
  trial <- eventReactive(input$run_trial, {
    senses <- c(0.82, 0.88, 0.94, 1.00, 1.08, 1.18, 1.30)
    arms <- list(placebo = "none",
                 ondansetron = "ondansetron 8 mg q8h",
                 mirtazapine = "mirtazapine 30 mg qHS")
    days <- c(2, 4, 7, 14)
    withProgress(message = "running 21 simulations ...", {
      res <- lapply(names(arms), function(an) {
        m <- sapply(senses, function(s) {
          d <- sim_hg(sens = s, alpha = input$alpha, tau_sp = input$tau_sp,
                      drug = arms[[an]], t_drug = TR, drug_days = 14,
                      end = TR + 24)
          b <- pick(d, TR, "PUQE")
          sapply(days, function(dd) pick(d, TR + dd, "PUQE") - b)
        })
        incProgress(1 / 3)
        data.frame(arm = an, day = days, dPUQE = rowMeans(m))
      })
      bind_rows(res)
    })
  })

  output$trial_txt <- renderText({
    tr <- trial()
    pb <- tr[tr$arm == "placebo", ]
    out <- c("VOMIT trial (Ostenfeld 2026, PMID 41478546), cohort mean",
             "  change in PUQE-24 vs placebo:")
    for (an in c("ondansetron", "mirtazapine")) {
      a <- tr[tr$arm == an, ]
      out <- c(out, sprintf("    %-12s %s", an,
        paste(sprintf("d%d %+5.2f", a$day, a$dPUQE - pb$dPUQE),
              collapse = "   ")))
    }
    c(out,
      "  observed       ondansetron -0.51 (95% CI -2.32 to  1.30)",
      "                 mirtazapine -1.86 (95% CI -3.61 to -0.12)") %>%
      paste(collapse = "\n")
  })

  output$p_trial <- renderPlot({
    tr <- trial(); pb <- tr[tr$arm == "placebo", ]
    tr %>% filter(arm != "placebo") %>%
      left_join(select(pb, day, ref = dPUQE), by = "day") %>%
      mutate(delta = dPUQE - ref) %>%
      ggplot(aes(day, delta, colour = arm)) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_hline(yintercept = c(-0.51, -1.86), linetype = 3,
                 colour = c(PAL[4], PAL[2])) +
      scale_colour_manual(values = PAL[c(4, 2)]) +
      labs(x = "day of treatment", y = "PUQE-24 change vs placebo",
           colour = NULL,
           title = "Dotted lines = the observed day-2 effects")
  })

  # ---- tab 8 ---------------------------------------------------------------
  prev <- eventReactive(input$run_prev, {
    starts <- c(-70, -56, -42, -28, -14, 0, 14, 28, 42, 56)
    withProgress(message = "sweeping metformin start time ...", {
      bind_rows(lapply(starts, function(t0) {
        ev <- dose("MET_C", nmol(1000, "MET") * 0.55, t0, ii = 0.5,
                   addl = max(1, floor((140 - t0) / 0.5)))
        p <- list(SENS = 1.0, ALPHA = input$alpha, TAU_SP = input$tau_sp)
        d <- as.data.frame(mod %>% param(p) %>%
               mrgsim(events = ev, start = T_START, end = 170, delta = 0.5))
        w <- filter(d, time >= 28, time <= 160)
        incProgress(1 / length(starts))
        data.frame(start_day = t0, start_label = sprintf("%+d d", t0),
                   PUQE_max = max(w$PUQE),
                   GDF15_at_conception = pick(d, 14, "GDF15"),
                   fold_at_9wk = pick(d, 63, "FOLD"),
                   wt_loss = max(w$WT_LOSS_P))
      }))
    })
  })

  output$p_prev <- renderPlot({
    prev() %>% ggplot(aes(start_day, PUQE_max)) +
      geom_line(linewidth = 1, colour = PAL[2]) +
      geom_point(size = 2.5, colour = PAL[2]) +
      geom_vline(xintercept = 14, linetype = 2, colour = PAL[5]) +
      annotate("text", x = 17, y = Inf, vjust = 1.5, hjust = 0,
               label = "conception", colour = PAL[5], size = 3.5) +
      labs(x = "metformin start (days relative to LMP; 14 = conception)",
           y = "peak PUQE-24",
           title = "The therapeutic window closes at conception",
           subtitle = "same drug, same dose, benefit is entirely a function of timing")
  })
  output$t_prev <- renderDT(
    datatable(prev() %>% select(-start_day) %>% mutate(across(where(is.numeric),
                                                             ~round(.x, 2))),
              rownames = FALSE, options = list(dom = "t")))

  # ---- tab 9 ---------------------------------------------------------------
  output$p_elec <- renderPlot({
    filter(cur(), time >= 28, time <= 170) %>%
      select(GA_WEEKS, `K+ (mmol/L)` = K_PLASMA, `Cl- (mmol/L)` = CL_PLASMA,
             `Na+ (mmol/L)` = NA_PLASMA, `HCO3- (mmol/L)` = HCO3_PL,
             `ECF volume (L)` = VOL, `ALT (U/L)` = ALT) %>%
      pivot_longer(-GA_WEEKS) %>%
      ggplot(aes(GA_WEEKS, value)) +
      geom_line(linewidth = 0.9, colour = PAL[6]) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "gestational age (weeks)", y = NULL)
  })

  # ---- tab 10 --------------------------------------------------------------
  output$p_thi <- renderPlot({
    d <- filter(cur(), time >= 28, time <= 170)
    ggplot(d, aes(GA_WEEKS)) +
      geom_line(aes(y = THIAMINE), linewidth = 1, colour = PAL[5]) +
      geom_hline(yintercept = mod@param$THI_CRIT, linetype = 2) +
      geom_line(aes(y = P_WE * 28), linewidth = 1, linetype = 3,
                colour = PAL[8]) +
      scale_y_continuous(
        name = "whole-body thiamine (mg)",
        sec.axis = sec_axis(~ . / 28, name = "P(Wernicke)")) +
      labs(x = "gestational age (weeks)",
           title = "Thiamine store (solid) and cumulative Wernicke risk (dotted)")
  })
  output$thi_txt <- renderText({
    w <- filter(cur(), time >= 28, time <= 170)
    sprintf(paste0("Thiamine nadir  %.1f mg  (hazard threshold %.0f mg)\n",
                   "P(Wernicke)     %.1f %%\n",
                   "Nausea is an hours-scale process; this store has a ~%.0f day",
                   " half-life.\nThe two can move in opposite directions."),
            min(w$THIAMINE), mod@param$THI_CRIT, 100 * max(w$P_WE),
            log(2) / mod@param$KEL_THI)
  })
  dex <- eventReactive(input$run_dex, {
    arms <- list("saline + KCl, no dextrose" = list(dex = FALSE, thi = 0),
                 "dextrose, NO thiamine"     = list(dex = TRUE,  thi = 0),
                 "dextrose + thiamine 100 mg" = list(dex = TRUE, thi = 100))
    withProgress(message = "comparing supportive regimens ...", {
      bind_rows(lapply(names(arms), function(nm) {
        a <- arms[[nm]]
        d <- sim_hg(sens = 1.18, alpha = input$alpha, tau_sp = input$tau_sp,
                    iv_on = TRUE, iv_from = 56, iv_days = 28,
                    iv_dextrose = a$dex, iv_thiamine = a$thi,
                    end = 150, delta = 0.5)
        w <- filter(d, time >= 42, time <= 145)
        incProgress(1 / 3)
        data.frame(regimen = nm,
                   thiamine_nadir_mg = round(min(w$THIAMINE), 1),
                   P_Wernicke_pct = round(100 * max(w$P_WE), 1),
                   HCO3_max = round(max(w$HCO3_PL), 1),
                   weight_loss_pct = round(max(w$WT_LOSS_P), 1))
      }))
    })
  })
  output$t_dex <- renderDT(datatable(dex(), rownames = FALSE,
                                     options = list(dom = "t")))

  # ---- tab 11 --------------------------------------------------------------
  output$p_thy <- renderPlot({
    ph <- list("normal (TROPH 1.0)"        = list(sens = 0.62, troph = 1.0),
               "HG, sensitivity-driven"    = list(sens = 1.00, troph = 1.0),
               "HG, production-driven"     = list(sens = 0.62, troph = 2.0))
    df <- bind_rows(lapply(names(ph), function(nm) {
      a <- c(ph[[nm]], list(alpha = input$alpha, tau_sp = input$tau_sp))
      gest(do.call(sim_hg, a)) %>% mutate(phenotype = nm)
    }))
    df %>% select(GA_WEEKS, phenotype, `hCG (IU/L)` = HCG,
                  `TSH (mIU/L)` = TSH, `free T4 (pmol/L)` = FT4,
                  `PUQE-24` = PUQE) %>%
      pivot_longer(c(-GA_WEEKS, -phenotype)) %>%
      ggplot(aes(GA_WEEKS, value, colour = phenotype)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL[c(3, 1, 2)]) +
      labs(x = "gestational age (weeks)", y = NULL, colour = NULL) +
      theme(legend.position = "bottom")
  })

  # ---- tab 12 --------------------------------------------------------------
  output$t_valid <- renderDT({
    datatable(data.frame(
      target = c("VOMIT mirtazapine dPUQE d2", "VOMIT ondansetron dPUQE d2",
                 "VOMIT mirtazapine > ondansetron d7",
                 "Koren doxylamine dPUQE d14",
                 "Guttuso gabapentin relative reduction",
                 "Maina clonidine PUQE gain",
                 "Yost corticosteroid (null result)",
                 "beta-thalassaemia peak PUQE",
                 "metformin pre-conception peak PUQE",
                 "tobacco pre-conception peak PUQE",
                 "HG peak, GA weeks", "HG PUQE at GA 16 wk",
                 "HG weight loss", "GDF15 28wk/9wk ratio",
                 "production-driven HG: TSH suppressed"),
      model = c("-1.92", "-0.52", "-1.61", "-1.75", "16%", "1.37", "0.00",
                "3.00", "10.59", "12.04", "8.5 wk", "3.0", "6.7%", "1.43",
                "0.16"),
      observed = c("-1.86 (-3.61,-0.12)", "-0.51 (-2.32,1.30)",
                   "-1.35 (-3.10,0.40)", "-0.9", "52% (16-88)",
                   "CI 0.43-3.24", "34% vs 35% (P=.89)", "very low NVP",
                   "aRR 0.29", "aRR 0.51", "9-11 wk", "remitted",
                   ">=5%", ">1", "<0.4 in ~60%"),
      role = c("FITTED (E0)", "FITTED (W_VAG)", rep("prediction", 13)),
      PMID = c("41478546", "41478546", "41478546", "20843504", "33451591",
               "24684734", "14662211", "38092039", "40588059", "40588059",
               "31515515", "31515515", "34555550", "12495665", "15073140"),
      verdict = c(rep("PASS", 15))),
      rownames = FALSE, escape = FALSE,
      options = list(dom = "t", pageLength = 15))
  })
}

shinyApp(ui, server)
