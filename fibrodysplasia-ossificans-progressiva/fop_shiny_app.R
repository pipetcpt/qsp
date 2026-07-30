# =====================================================================
# Fibrodysplasia Ossificans Progressiva (FOP) — Shiny dashboard
#   Front end for fop_mrgsolve_model.R (42-ODE QSP model).
#   Run:  shiny::runApp("fop_shiny_app.R")
#
#   The tabs are organised around the four axes of the model rather than
#   around organ systems, because in FOP the interesting quantities are
#   all comparisons: rate against stock, age against potency, efficacy
#   against physeal toxicity, and endpoint against mechanism.
#
#   11 tabs:
#     1  Patient & regimen      — set up the virtual patient
#     2  Drug PK                — all seven agents
#     3  Signalling & lesion    — activin A, pSMAD1/5/8, the relay
#     4  HO burden: RATE vs STOCK (Axis 2)
#     5  Clinical endpoints     — CAJIS, FVC, survival
#     6  Age at start (Axis 3)
#     7  Dose vs physeal toxicity (Axis 4)
#     8  Anti-flare ceiling (Axis 1)
#     9  Mechanism-class ceilings
#    10  Imaging endpoints      — WBCT vs NaF PET dilution
#    11  Trial design           — variability and external-control bias
# =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

source("fop_mrgsolve_model.R")

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#eceff1"),
        legend.position = "bottom")

PAL <- c(untreated = "#455a64", treated = "#c62828", alt = "#00695c",
         third = "#ef6c00", fourth = "#4527a0")

# ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Fibrodysplasia Ossificans Progressiva — QSP model explorer"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Virtual patient"),
      sliderInput("age0", "Age at treatment start (yr)", 4, 45, 15, step = 1),
      sliderInput("wt", "Body weight (kg)", 10, 90, 55, step = 5),
      radioButtons("sex", "Sex (sets physeal fusion age)",
                   c("Female (14 yr)" = 0, "Male (16 yr)" = 1), inline = TRUE),
      sliderInput("years", "Simulation length (yr)", 0.5, 40, 3, step = 0.5),
      hr(),
      h4("Regimen"),
      selectInput("drug", "Agent",
                  c("none",
                    "palovarotene (RARgamma, oral)",
                    "garetosmab (anti-activin A, IV q4w)",
                    "ALK2 kinase inhibitor (oral)",
                    "saracatinib (SRC/ALK2, oral)",
                    "sirolimus (mTORC1, oral)",
                    "IL-1 blockade (anakinra, SC)")),
      numericInput("dose", "Daily dose (mg; mg/kg for the mAb)", 5, 0, 1000),
      checkboxInput("flare_dose", "Add MOVE flare-up dosing (20 mg x 4 wk, then 10 mg x 8 wk)", FALSE),
      hr(),
      h4("Flare-up events"),
      checkboxInput("do_flare", "Simulate a discrete flare-up", TRUE),
      sliderInput("flare_at", "Flare at month", 0, 24, 6, step = 1),
      sliderInput("flare_mag", "Flare magnitude (1 = typical, 6 = IM injection)",
                  1, 10, 1, step = 1),
      checkboxInput("pred", "Prednisone burst 2 mg/kg/d x 4 d at the flare", FALSE),
      hr(),
      h4("Structural switches"),
      sliderInput("phi", "Flare-attributable fraction of new HO (PHI_FL)",
                  0.05, 0.9, 0.295, step = 0.005),
      sliderInput("leak", "Ligand-independent receptor leak (F_LEAK)",
                  0, 0.6, 0.10, step = 0.01),
      sliderInput("eps", "Anti-ACVR1 agonist conversion (EPS_AGON)",
                  0, 1, 0, step = 0.05),
      sliderInput("rsel", "Physeal selectivity R_SEL = EC50(PPC)/EC50(HO)",
                  0.5, 40, 3.698, step = 0.1)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 Patient & regimen",
                 h4("Starting burden and what the model assumes about it"),
                 verbatimTextOutput("setup"),
                 h4("Calibration against published anchors"),
                 tableOutput("calib")),
        tabPanel("2 Drug PK", plotOutput("pk", height = 520),
                 helpText("Concentrations in ng/mL for small molecules and mg/L for garetosmab.")),
        tabPanel("3 Signalling & lesion", plotOutput("sig", height = 560),
                 helpText(paste("The receptor step is deliberately LINEAR in activin A so that the",
                                "episodic and smouldering routes stay exactly separable."))),
        tabPanel("4 HO: rate vs stock", plotOutput("burden", height = 480),
                 tableOutput("burden_tab"),
                 helpText(paste("Axis 2: the trial endpoint is the slope, the disability is the area.",
                                "A large reduction in the first is a small reduction in the second."))),
        tabPanel("5 Clinical endpoints", plotOutput("clin", height = 560),
                 helpText(paste("CAJIS is on its true 0-30 scale but is NOT numerically calibrated",
                                "to a published cohort mean; read it as ordinal."))),
        tabPanel("6 Age at start (Axis 3)", plotOutput("agestart", height = 420),
                 tableOutput("agestart_tab")),
        tabPanel("7 Dose vs physis (Axis 4)", plotOutput("tox", height = 460),
                 tableOutput("tox_tab"),
                 helpText(paste("Efficacy and premature physeal closure share one exposure axis.",
                                "Move the R_SEL slider to see what selectivity a next-generation",
                                "retinoid would need."))),
        tabPanel("8 Anti-flare ceiling (Axis 1)", plotOutput("ceiling", height = 420),
                 tableOutput("ceiling_tab")),
        tabPanel("9 Class ceilings", tableOutput("classceil"),
                 plotOutput("classplot", height = 400)),
        tabPanel("10 Imaging endpoints", plotOutput("imaging", height = 480),
                 tableOutput("imaging_tab"),
                 helpText(paste("Total NaF lesion activity mixes new lesions with remodelling of the",
                                "mature stock. That is the dilution that sank the LUMINA-1 primary."))),
        tabPanel("11 Trial design", h4("Between-subject variability"),
                 tableOutput("trialvar"),
                 h4("External-control bias in a single-arm design"),
                 tableOutput("extctrl"),
                 plotOutput("extplot", height = 340))
      )
    )
  )
)

# ---------------------------------------------------------------------
server <- function(input, output, session) {

  dur <- reactive(max(30, round(yr(input$years))))

  pars <- reactive({
    # PHI_FL is set by rebalancing the episodic activin A source
    p    <- as.list(param(fop))
    sigsm <- p$KSIG_A * (p$KIN_ASM / p$KOUT_ACTA) / (1 - input$leak)
    kin_afl <- input$phi / (1 - input$phi) * sigsm * p$KOUT_ACTA / p$KSIG_A
    list(SEX = as.numeric(input$sex), WT = input$wt, F_LEAK = input$leak,
         KIN_AFL = kin_afl, EPS_AGON = input$eps, R_SEL = input$rsel)
  })

  events <- reactive({
    d <- dur(); e <- NULL
    cmt <- switch(input$drug,
                  "palovarotene (RARgamma, oral)"        = CMT_PAL,
                  "garetosmab (anti-activin A, IV q4w)"  = CMT_GAR,
                  "ALK2 kinase inhibitor (oral)"         = CMT_AK,
                  "saracatinib (SRC/ALK2, oral)"         = CMT_SAR,
                  "sirolimus (mTORC1, oral)"             = CMT_SIR,
                  "IL-1 blockade (anakinra, SC)"         = CMT_IL1, NULL)
    if (!is.null(cmt) && input$dose > 0) {
      e <- if (identical(cmt, CMT_GAR))
        ev(amt = input$dose * input$wt, ii = 28, addl = floor(d / 28), cmt = cmt)
      else oral_daily(input$dose, 0, d, cmt)
    }
    ft <- round(yr(input$flare_at / 12))
    if (input$flare_dose && identical(cmt, CMT_PAL))
      e <- c(e, oral_daily(20, ft, 28, CMT_PAL), oral_daily(10, ft + 28, 56, CMT_PAL))
    if (input$do_flare) {
      e <- c(e, flare_events(ft, amt = input$flare_mag))
      if (input$pred) e <- c(e, oral_daily(2 * input$wt, ft, 4, CMT_PRE))
    }
    e
  })

  base_run <- reactive({
    ev0 <- if (input$do_flare)
      flare_events(round(yr(input$flare_at / 12)), amt = input$flare_mag) else NULL
    do.call(sim_pt, c(list(input$age0, dur(), max(1, round(dur() / 300)), ev0), pars()))
  })

  trt_run <- reactive({
    do.call(sim_pt, c(list(input$age0, dur(), max(1, round(dur() / 300)), events()),
                      pars()))
  })

  both <- reactive({
    bind_rows(mutate(base_run(), arm = "untreated"),
              mutate(trt_run(),  arm = "treated"))
  })

  # ---- 1. setup ----
  output$setup <- renderPrint({
    cat("Age at start                    :", input$age0, "yr\n")
    cat("Heterotopic bone already formed :", round(ho_at_age(input$age0), 0), "mL\n")
    cat("  (from the natural-history trajectory; a scenario that ignores this\n")
    cat("   gets every percent-of-total-volume statement wrong)\n")
    cat("Ossifiable territory remaining  :",
        round(100 * (1 - ho_at_age(input$age0) / 2000), 1), "%\n")
    cat("Growth plates                   :",
        ifelse(input$age0 < 14 + 2 * as.numeric(input$sex), "OPEN - physeal hazard applies",
               "fused - no physeal hazard"), "\n")
    cat("Flare-attributable share of HO  :", input$phi, "\n")
    cat("Anti-flare therapy ceiling      :", round(100 * input$phi, 1), "%\n")
  })
  output$calib <- renderTable(calibration_report())

  # ---- 2. PK ----
  output$pk <- renderPlot({
    d <- trt_run() %>%
      select(time, C_PALO, C_GARE, C_ALK2I, C_SIRO) %>%
      pivot_longer(-time)
    ggplot(d, aes(time / 30.4, value)) + geom_line(colour = PAL["treated"], linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "months", y = "concentration", title = "Drug exposure") + THEME
  })

  # ---- 3. signalling ----
  output$sig <- renderPlot({
    d <- both() %>%
      select(time, arm, ACTA_FREE, SIG_TOTAL, FRAC_FLARE, LESIONS) %>%
      pivot_longer(c(-time, -arm))
    ggplot(d, aes(time / 30.4, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") + scale_colour_manual(values = PAL) +
      labs(x = "months", y = NULL, title = "Ligand, signal, route split, lesion count") + THEME
  })

  # ---- 4. rate vs stock ----
  output$burden <- renderPlot({
    d <- both() %>% select(time, arm, HO_TOTAL, NEW_HO) %>%
      pivot_longer(c(-time, -arm)) %>%
      mutate(name = recode(name, HO_TOTAL = "TOTAL stock (mL) - the disability",
                           NEW_HO = "NEW HO since t0 (mL) - the trial endpoint"))
    ggplot(d, aes(time / 30.4, value, colour = arm)) + geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") + scale_colour_manual(values = PAL) +
      labs(x = "months", y = "mL") + THEME
  })
  output$burden_tab <- renderTable({
    b <- base_run(); t <- trt_run(); d <- dur()
    ru <- ann_new_ho(b, 0, d); rt <- ann_new_ho(t, 0, d)
    data.frame(quantity = c("annualised new HO (mL/yr)", "TOTAL stock at end (mL)"),
               untreated = round(c(ru, tail(b$HO_TOTAL, 1)), 2),
               treated   = round(c(rt, tail(t$HO_TOTAL, 1)), 2),
               reduction_pct = round(100 * (1 - c(rt / ru,
                                 tail(t$HO_TOTAL, 1) / tail(b$HO_TOTAL, 1))), 1))
  })

  # ---- 5. clinical ----
  output$clin <- renderPlot({
    d <- both() %>% mutate(age = input$age0 + time / 365.25) %>%
      select(age, arm, CAJIS, FVC_PCT, SURV, HO_THOR) %>%
      pivot_longer(c(-age, -arm))
    ggplot(d, aes(age, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") + scale_colour_manual(values = PAL) +
      labs(x = "age (yr)", y = NULL) + THEME
  })

  # ---- 6. age at start ----
  output$agestart <- renderPlot({
    d <- sc05_age_at_start()
    ggplot(d, aes(factor(start_age), pct_burden_averted_at_40)) +
      geom_col(fill = PAL["alt"]) +
      geom_text(aes(label = paste0(pct_burden_averted_at_40, "%")), vjust = -0.4) +
      labs(x = "age at treatment start (yr)",
           y = "% of the age-40 burden averted",
           title = "Same drug, same efficacy - the start age decides the outcome") + THEME
  })
  output$agestart_tab <- renderTable(sc05_age_at_start())

  # ---- 7. dose vs physis ----
  tox <- reactive({
    d <- dur()
    doses <- c(0, 0.5, 1, 2, 3, 5, 7.5, 10)
    base <- ann_new_ho(do.call(sim_pt, c(list(input$age0, d, 30.4), pars())), 0, d)
    do.call(rbind, lapply(doses, function(x) {
      o <- do.call(sim_pt, c(list(input$age0, d, 30.4,
             if (x == 0) NULL else oral_daily(x, 0, d, CMT_PAL)), pars()))
      data.frame(dose_mg = x,
                 HO_suppression_pct = round(100 * (1 - ann_new_ho(o, 0, d) / base), 1),
                 PPC_probability_pct = round(100 * tail(o$PPC_PROB, 1), 1),
                 height_Z_loss = round(-tail(o$HEIGHT_Z, 1), 2))
    }))
  })
  output$tox <- renderPlot({
    ggplot(tox(), aes(HO_suppression_pct, PPC_probability_pct)) +
      geom_path(colour = PAL["treated"], linewidth = 1) +
      geom_point(size = 2.5, colour = PAL["treated"]) +
      geom_text(aes(label = paste0(dose_mg, " mg")), hjust = -0.2, size = 3.2) +
      geom_hline(yintercept = 10, linetype = 2, colour = "grey40") +
      annotate("text", x = 5, y = 12, label = "10% physeal hazard", hjust = 0, size = 3.2) +
      labs(x = "% suppression of new HO", y = "premature physeal closure probability (%)",
           title = "One exposure axis, two effects: the trade-off cannot be dosed away") + THEME
  })
  output$tox_tab <- renderTable(tox())

  # ---- 8. anti-flare ceiling ----
  output$ceiling <- renderPlot({
    d <- dur()
    fr <- seq(0, 1, 0.1)
    base <- ann_new_ho(do.call(sim_pt, c(list(input$age0, d, 30.4), pars())), 0, d)
    p0 <- pars()
    res <- do.call(rbind, lapply(fr, function(f) {
      p <- p0; p$KIN_AFL <- p0$KIN_AFL * (1 - f)
      o <- do.call(sim_pt, c(list(input$age0, d, 30.4, NULL), p))
      data.frame(flare_suppression = 100 * f,
                 HO_reduction = 100 * (1 - ann_new_ho(o, 0, d) / base))
    }))
    ggplot(res, aes(flare_suppression, HO_reduction)) +
      geom_line(linewidth = 1, colour = PAL["fourth"]) +
      geom_hline(yintercept = 100 * input$phi, linetype = 2) +
      geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey50") +
      annotate("text", x = 5, y = 100 * input$phi + 3,
               label = paste0("ceiling = PHI_FL = ", round(100 * input$phi, 1), "%"),
               hjust = 0, size = 3.5) +
      labs(x = "% of flare-ups abolished", y = "% reduction in new HO",
           title = "Axis 1: perfect flare control is not disease control") + THEME
  })
  output$ceiling_tab <- renderTable(sc07_antiflare_ceiling(age0 = 15))

  # ---- 9. class ceilings ----
  output$classceil <- renderTable(class_ceilings())
  output$classplot <- renderPlot({
    d <- class_ceilings()
    ggplot(d, aes(reorder(mechanism_class, analytic_ceiling_pct), analytic_ceiling_pct)) +
      geom_col(fill = PAL["fourth"]) + coord_flip() +
      geom_text(aes(label = paste0(analytic_ceiling_pct, "%")), hjust = -0.15) +
      ylim(0, 110) + labs(x = NULL, y = "maximal achievable reduction in new HO (%)",
           title = "Where each mechanism class caps out, independent of potency") + THEME
  })

  # ---- 10. imaging ----
  output$imaging <- renderPlot({
    d <- both() %>% select(time, arm, NAF_TOTAL, NAF_NEW, NAF_OLD, NAF_FRAC_NEW) %>%
      pivot_longer(c(-time, -arm))
    ggplot(d, aes(time / 30.4, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") + scale_colour_manual(values = PAL) +
      labs(x = "months", y = "NaF PET activity (a.u.)") + THEME
  })
  output$imaging_tab <- renderTable(sc06_lumina1())

  # ---- 11. trial design ----
  output$trialvar <- renderTable(sc14_trial_variability(60))
  output$extctrl  <- renderTable(sc15_external_control_bias(age0 = input$age0))
  output$extplot  <- renderPlot({
    d <- sc15_external_control_bias(age0 = input$age0)
    ggplot(d, aes(trial_vs_NHS_intrinsic_activity)) +
      geom_line(aes(y = apparent_reduction_vs_NHS_pct, colour = "apparent vs external control"),
                linewidth = 1) +
      geom_line(aes(y = true_drug_effect_pct, colour = "true drug effect"), linewidth = 1) +
      geom_line(aes(y = apparent_reduction_with_NO_drug_pct,
                    colour = "apparent WITHOUT any drug"), linewidth = 1, linetype = 2) +
      scale_colour_manual(values = c("apparent vs external control" = PAL[["treated"]],
                                     "true drug effect" = PAL[["alt"]],
                                     "apparent WITHOUT any drug" = PAL[["third"]])) +
      labs(x = "intrinsic activity of the trial cohort relative to the external control",
           y = "% reduction in new HO", colour = NULL,
           title = "Cohort selection alone can manufacture the whole effect size") + THEME
  })
}

shinyApp(ui, server)
