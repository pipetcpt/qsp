##############################################################################
## flu_shiny_app.R
## Influenza A QSP model — interactive dashboard
## ===========================================================================
##
## The app is organised around the model's central claim rather than around
## its compartments.  Tab 2 is the point: it draws the RESIDUAL AUC R(t), the
## exact upper bound on what an antiviral started at t can remove, and shows
## how much of that bound the chosen regimen actually takes.  Everything else
## is either the machinery that produces R(t) or a consequence of it.
##
##   1  Patient & regimen      who is being simulated and what they receive
##   2  The bound              R(t), achieved reduction, and the two headrooms
##   3  Viral kinetics         URT / LRT titres, wild type vs resistant
##   4  Epithelium & immunity  target pool, interferon, CD8, antibody, IgA
##   5  Drug PK/PD             concentrations and the fractional blockade
##   6  Clinical endpoints     symptom score, TTAS rule, fever, SpO2
##   7  Scenario comparison    the shipped arms side by side
##   8  Resistance             competitive release over dose time x potency
##   9  Operator decomposition which term of the loop each therapy touches
##  10  Trial ledger           model vs published, and the five discrepancies
##
## Run:
##   R -e "shiny::runApp('flu_shiny_app.R', port = 8080)"
##
## Requires flu_mrgsolve_model.R in the same directory.
##############################################################################

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

source("flu_mrgsolve_model.R", local = TRUE)

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10),
        legend.position = "bottom")

PAL <- c(placebo = "#8a8a8a", oseltamivir = "#1f77b4", baloxavir = "#7b3fa0",
         peramivir = "#17becf", favipiravir = "#8b5a2b", mAb = "#2ca02c",
         combination = "#d62728", "wild type" = "#c0392b",
         resistant = "#d97706")

REGIMENS <- c("none (placebo)", "oseltamivir 75 mg BID x5d",
              "baloxavir 40 mg single", "peramivir 600 mg IV",
              "favipiravir 1800/800 BID x5d", "anti-HA mAb IV",
              "baloxavir + oseltamivir")

build_doses <- function(regimen, t_rx, steroid) {
  d <- switch(regimen,
    "none (placebo)"                = NULL,
    "oseltamivir 75 mg BID x5d"     = rx_oseltamivir(t_rx),
    "baloxavir 40 mg single"        = rx_baloxavir(t_rx),
    "peramivir 600 mg IV"           = rx_peramivir(t_rx),
    "favipiravir 1800/800 BID x5d"  = rx_favipiravir(t_rx),
    "anti-HA mAb IV"                = rx_mab(t_rx),
    "baloxavir + oseltamivir"       = rbind(rx_baloxavir(t_rx),
                                            rx_oseltamivir(t_rx)))
  if (steroid) d <- rbind(d, rx_dexamethasone(t_rx))
  d
}

##############################################################################
## UI
##############################################################################
ui <- fluidPage(
  titlePanel("Influenza A — QSP model: the bound, the operators, and the release"),
  tags$p(style = "color:#555;margin-top:-8px;",
         "An antiviral cannot subtract viral load that has already happened. ",
         tags$b("R(t) = the residual viral AUC"), " is the exact upper bound on ",
         "what a drug started at t can remove — no potency crosses it."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("host", "Host phenotype",
                  c("adult, prior exposure", "adult, naive",
                    "child (naive, larger epithelium)",
                    "elderly (weak CD8)", "immunocompromised",
                    "vaccinated (mucosal IgA)"),
                  selected = "adult, naive"),
      sliderInput("t0u", "URT target-cell pool (×10⁸)", 2, 8, 4, 0.5),
      sliderInput("fifn", "Interferon capacity (×)", 0.2, 2, 1, 0.05),
      sliderInput("fctl", "CD8 capacity (×)", 0.05, 2, 1, 0.05),
      sliderInput("fab",  "Humoral capacity (×)", 0.05, 2, 1, 0.05),
      sliderInput("iga0", "Pre-existing mucosal IgA", 0, 1.2, 0, 0.05),

      hr(), h4("Regimen"),
      selectInput("regimen", "Antiviral", REGIMENS,
                  selected = "baloxavir 40 mg single"),
      sliderInput("delay", "First dose (h after symptom onset)",
                  -24, 120, 24, 6),
      checkboxInput("steroid", "add dexamethasone 6 mg daily × 5 d", FALSE),

      hr(), h4("Resistance"),
      selectInput("resprof", "Resistance profile",
                  c("none (drug-susceptible)", "PA/I38T (baloxavir, ×50)",
                    "NA/H275Y (oseltamivir, ×300)"),
                  selected = "none (drug-susceptible)"),
      sliderInput("cost", "Replicative fitness cost of the mutant",
                  0, 0.6, 0.18, 0.02),

      hr(), h4("Pharmacodynamics"),
      sliderInput("emaxbx", "Baloxavir Emax", 0.5, 0.999999, 0.9999,
                  step = 0.0001),
      sliderInput("hillbx", "Baloxavir Hill slope", 0.5, 3, 2, 0.1),
      sliderInput("emaxnai", "NAI Emax", 0.5, 0.999, 0.995, 0.001),
      sliderInput("wvir", "Symptom drive tracking titre (WVIR)", 0, 1, 0.6, 0.05),

      hr(),
      sliderInput("tmax", "Simulation horizon (days)", 14, 40, 30, 2),
      helpText("Every default is the value calibrated in ",
               tags$code("flu_reference_check.py"), ".")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient & regimen",
                 br(), verbatimTextOutput("summary_txt"),
                 h5("Calibration against untreated natural history"),
                 tableOutput("calib_tbl")),
        tabPanel("2 · The bound",
                 br(),
                 plotOutput("bound_plot", height = "330px"),
                 plotOutput("headroom_plot", height = "280px"),
                 tableOutput("bound_tbl")),
        tabPanel("3 · Viral kinetics",
                 br(),
                 plotOutput("viral_plot", height = "340px"),
                 plotOutput("strain_plot", height = "280px")),
        tabPanel("4 · Epithelium & immunity",
                 br(),
                 plotOutput("epi_plot", height = "300px"),
                 plotOutput("imm_plot", height = "320px")),
        tabPanel("5 · Drug PK/PD",
                 br(),
                 plotOutput("pk_plot", height = "300px"),
                 plotOutput("pd_plot", height = "300px"),
                 tableOutput("pd_tbl")),
        tabPanel("6 · Clinical endpoints",
                 br(),
                 plotOutput("sym_plot", height = "320px"),
                 plotOutput("clin_plot", height = "300px"),
                 tableOutput("endpoint_tbl")),
        tabPanel("7 · Scenario comparison",
                 br(), tableOutput("scen_tbl"),
                 plotOutput("scen_plot", height = "340px")),
        tabPanel("8 · Resistance",
                 br(),
                 plotOutput("release_plot", height = "340px"),
                 tableOutput("release_tbl"),
                 tags$p(style = "color:#a33;",
                        tags$b("Deterministic caveat. "),
                        "A lineage at 10⁻⁴ is 'present' in every simulated ",
                        "patient. Real I38T emergence is a stochastic ",
                        "establishment event in ~10% of treated adults. ",
                        "Read this as selection pressure, not viral load.")),
        tabPanel("9 · Operators",
                 br(), tableOutput("oper_tbl"),
                 plotOutput("oper_plot", height = "320px")),
        tabPanel("10 · Trial ledger",
                 br(), tableOutput("ledger_tbl"),
                 htmlOutput("discrepancies"))
      )
    )
  )
)

##############################################################################
## SERVER
##############################################################################
server <- function(input, output, session) {

  host_overrides <- reactive({
    switch(input$host,
      "adult, prior exposure"             = list(CTL0 = 6.0, IGA0 = 0.30),
      "adult, naive"                      = list(),
      "child (naive, larger epithelium)"  = list(T0U = 6e8, FAB = 0.65),
      "elderly (weak CD8)"                = list(FCTL = 0.35, FAB = 0.6),
      "immunocompromised"                 = list(FCTL = 0.12, FAB = 0.15,
                                                 FIFN = 0.55),
      "vaccinated (mucosal IgA)"          = list(IGA0 = 0.55, CTL0 = 6.0))
  })

  res_overrides <- reactive({
    switch(input$resprof,
      "none (drug-susceptible)"       = list(RF_BX = 1, RF_NAI = 1),
      "PA/I38T (baloxavir, ×50)"      = list(RF_BX = 50, RF_NAI = 1),
      "NA/H275Y (oseltamivir, ×300)"  = list(RF_BX = 1, RF_NAI = 300))
  })

  overrides <- reactive({
    ov <- modifyList(host_overrides(), res_overrides())
    ov <- modifyList(ov, list(
      T0U = input$t0u * 1e8, FIFN = input$fifn, FCTL = input$fctl,
      FAB = input$fab, IGA0 = input$iga0, COST = input$cost,
      EMAX_BX = input$emaxbx, HILL_BX = input$hillbx,
      EMAX_NAI = input$emaxnai, EMAX_PR = input$emaxnai, WVIR = input$wvir))
    ## host presets set CTL0/T0U/etc.; explicit sliders win where they overlap
    if (!is.null(host_overrides()$CTL0)) ov$CTL0 <- host_overrides()$CTL0
    ov
  })

  base_run <- reactive({
    do.call(sim, c(list(doses = NULL, tmax = input$tmax), overrides()))
  })
  onset <- reactive(symptom_onset(base_run()))
  t_rx  <- reactive(max(0, onset() + input$delay / 24))

  treated <- reactive({
    d <- build_doses(input$regimen, t_rx(), input$steroid)
    do.call(sim, c(list(doses = d, tmax = input$tmax), overrides()))
  })

  ## ---- tab 1 -------------------------------------------------------------
  output$summary_txt <- renderPrint({
    o <- base_run(); tr <- treated(); r <- summarise_run(tr, t_rx())
    cat("Host                    :", input$host, "\n")
    cat("Regimen                 :", input$regimen,
        if (input$steroid) "+ dexamethasone" else "", "\n")
    cat("Resistance profile      :", input$resprof, "\n\n")
    cat(sprintf("Symptom onset           : %.1f h post-infection\n", onset() * 24))
    cat(sprintf("First dose              : %.1f h p.i. (%+d h vs onset)\n",
                t_rx() * 24, input$delay))
    cat(sprintf("Untreated peak titre    : %.2f log10 at %.1f h\n",
                max(o$LOGV), o$time[which.max(o$LOGV)] * 24))
    cat(sprintf("Target pool at dosing   : %.2f %% of T0\n",
                o$TFRAC[which.min(abs(o$time - t_rx()))] * 100))
    cat("\n--- treated arm ---\n")
    cat(sprintf("Time to alleviation     : %s h\n", format(round(r$ttas_h, 1))))
    cat(sprintf("Shedding cessation      : %s h after dose\n",
                format(round(r$shed_h, 1))))
    cat(sprintf("24-h titre change       : %.2f log10\n", r$drop24))
    cat(sprintf("Viral AUC               : %.2f log10.d  (untreated %.2f)\n",
                r$auc_log, viral_auc(o)))
    cat(sprintf("Epithelium lost         : %.1f %% of T0\n", r$epi_lost * 100))
    cat(sprintf("Peak resistant titre    : %.2f log10\n", r$mut_peak))
  })

  output$calib_tbl <- renderTable(FLU_baseline(), digits = 2)

  ## ---- tab 2: the bound --------------------------------------------------
  bound_data <- reactive({
    o <- base_run(); R <- residual_auc(o)
    data.frame(time_h = o$time * 24, R = R, pct = R / R[1] * 100,
               logv = o$LOGV)
  })

  output$bound_plot <- renderPlot({
    d <- bound_data(); ons <- onset() * 24
    ggplot(d, aes(time_h, pct)) +
      geom_area(fill = "#ffd0d0", alpha = 0.6) +
      geom_line(colour = "#c0392b", linewidth = 1) +
      geom_vline(xintercept = ons, linetype = 2, colour = "#555") +
      geom_vline(xintercept = t_rx() * 24, linetype = 1, colour = "#1f77b4",
                 linewidth = 0.9) +
      annotate("text", x = ons, y = 96, label = " symptom onset",
               hjust = 0, size = 3.4, colour = "#555") +
      annotate("text", x = t_rx() * 24, y = 85, label = " first dose",
               hjust = 0, size = 3.4, colour = "#1f77b4") +
      labs(title = "R(t) — the residual viral AUC still ahead of the patient",
           subtitle = paste0("Everything any antiviral started at t could ever ",
                             "remove lies underneath this curve."),
           x = "hours post-infection", y = "% of total viral AUC remaining") +
      THEME
  })

  output$headroom_plot <- renderPlot({
    o <- base_run(); R <- residual_auc(o); total <- R[1]
    ons <- onset()
    hs  <- c(0, 6, 12, 24, 36, 48, 72, 96)
    d <- do.call(rbind, lapply(hs, function(h) {
      trx <- ons + h / 24; i <- which.min(abs(o$time - trx))
      tr  <- do.call(sim, c(list(doses = build_doses(input$regimen, trx,
                                                     input$steroid),
                                 tmax = input$tmax), overrides()))
      data.frame(h = h, bound = R[i], achieved = total - viral_auc(tr))
    }))
    d$left <- d$bound - d$achieved
    dl <- pivot_longer(d[, c("h", "achieved", "left")], -h)
    ggplot(dl, aes(factor(h), value, fill = name)) +
      geom_col(width = 0.7) +
      scale_fill_manual(values = c(achieved = "#2ca02c", left = "#e0e0e0"),
                        labels = c("taken by the drug", "left on the table"),
                        name = NULL) +
      labs(title = "The bound falls faster than any drug can close on it",
           subtitle = "bar height = R(t); green = the reduction actually achieved",
           x = "hours after symptom onset at which the drug is given",
           y = "log10 TCID50/mL · day") +
      THEME
  })

  output$bound_tbl <- renderTable({
    o <- base_run(); R <- residual_auc(o); total <- R[1]; ons <- onset()
    do.call(rbind, lapply(c(0, 12, 24, 48, 72), function(h) {
      trx <- ons + h / 24; i <- which.min(abs(o$time - trx))
      tr  <- do.call(sim, c(list(doses = build_doses(input$regimen, trx,
                                                     input$steroid),
                                 tmax = input$tmax), overrides()))
      ach <- total - viral_auc(tr)
      data.frame(`h after onset` = h,
                 `T remaining %` = o$TFRAC[i] * 100,
                 `bound R` = R[i], `% of total` = R[i] / total * 100,
                 achieved = ach, `efficiency %` = ach / R[i] * 100,
                 check.names = FALSE)
    }))
  }, digits = 2)

  ## ---- tab 3: viral kinetics ---------------------------------------------
  output$viral_plot <- renderPlot({
    o <- base_run(); tr <- treated()
    d <- rbind(
      data.frame(t = o$time * 24,  y = o$LOGV,  arm = "untreated", site = "URT"),
      data.frame(t = tr$time * 24, y = tr$LOGV, arm = "treated",   site = "URT"),
      data.frame(t = o$time * 24,  y = o$LOGVL, arm = "untreated", site = "LRT"),
      data.frame(t = tr$time * 24, y = tr$LOGVL, arm = "treated",  site = "LRT"))
    ggplot(d, aes(t, y, colour = arm, linetype = site)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = LOD, linetype = 3, colour = "#999") +
      geom_vline(xintercept = t_rx() * 24, colour = "#1f77b4", linetype = 2) +
      scale_colour_manual(values = c(untreated = "#8a8a8a", treated = "#7b3fa0"),
                          name = NULL) +
      coord_cartesian(ylim = c(-2, 8)) +
      labs(title = "Viral titre, upper and lower respiratory tract",
           subtitle = "dotted line = assay floor; blue dashed = first dose",
           x = "hours post-infection", y = "log10 TCID50/mL") + THEME
  })

  output$strain_plot <- renderPlot({
    tr <- treated()
    d <- rbind(data.frame(t = tr$time * 24, y = tr$LOGVW, strain = "wild type"),
               data.frame(t = tr$time * 24, y = tr$LOGVM, strain = "resistant"))
    ggplot(d, aes(t, y, colour = strain)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = LOD, linetype = 3, colour = "#999") +
      geom_vline(xintercept = t_rx() * 24, colour = "#1f77b4", linetype = 2) +
      scale_colour_manual(values = PAL[c("wild type", "resistant")], name = NULL) +
      coord_cartesian(ylim = c(-6, 8)) +
      labs(title = "Wild type versus resistant subpopulation",
           subtitle = paste("The mutant is never seeded — it arises from the",
                            "mutation term and grows when its competitor is removed."),
           x = "hours post-infection", y = "log10 TCID50/mL") + THEME
  })

  ## ---- tab 4: epithelium and immunity ------------------------------------
  output$epi_plot <- renderPlot({
    o <- base_run(); tr <- treated()
    d <- rbind(
      data.frame(t = o$time * 24,  y = o$TFRAC * 100,  arm = "untreated"),
      data.frame(t = tr$time * 24, y = tr$TFRAC * 100, arm = "treated"))
    ggplot(d, aes(t, y, colour = arm)) + geom_line(linewidth = 1) +
      geom_vline(xintercept = t_rx() * 24, colour = "#1f77b4", linetype = 2) +
      scale_colour_manual(values = c(untreated = "#8a8a8a", treated = "#7b3fa0"),
                          name = NULL) +
      labs(title = "The susceptible epithelium — the resource being spent",
           subtitle = "the level at the dose line is what the drug still has to protect",
           x = "hours post-infection", y = "% of T0 still susceptible") + THEME
  })

  output$imm_plot <- renderPlot({
    tr <- treated()
    d <- rbind(
      data.frame(t = tr$time * 24, y = tr$FU / max(tr$FU + 1e-12),  v = "interferon (URT)"),
      data.frame(t = tr$time * 24, y = tr$CTL / max(tr$CTL),        v = "CD8 effectors"),
      data.frame(t = tr$time * 24, y = tr$AB  / max(tr$AB + 1e-12), v = "serum IgG"),
      data.frame(t = tr$time * 24, y = tr$IGA / max(tr$IGA + 1e-12), v = "mucosal IgA"),
      data.frame(t = tr$time * 24, y = tr$IL6 / max(tr$IL6),        v = "IL-6"))
    ggplot(d, aes(t, y, colour = v)) + geom_line(linewidth = 0.9) +
      labs(title = "Immune effectors, each scaled to its own maximum",
           subtitle = "the ordering in time is the point, not the heights",
           x = "hours post-infection", y = "fraction of own peak",
           colour = NULL) + THEME
  })

  ## ---- tab 5: PK/PD ------------------------------------------------------
  output$pk_plot <- renderPlot({
    tr <- treated()
    d <- rbind(
      data.frame(t = tr$time * 24, y = tr$COCE,   drug = "oseltamivir carboxylate (ELF)"),
      data.frame(t = tr$time * 24, y = tr$CBXE,   drug = "baloxavir acid, free (ELF)"),
      data.frame(t = tr$time * 24, y = tr$CBXTOT, drug = "baloxavir acid, total plasma"))
    d <- d[d$y > 1e-4, ]
    if (!nrow(d)) return(NULL)
    ggplot(d, aes(t, y, colour = drug)) + geom_line(linewidth = 0.9) +
      scale_y_log10() +
      labs(title = "Drug concentrations at the effect site",
           x = "hours post-infection", y = "ng/mL (log scale)", colour = NULL) +
      THEME
  })

  output$pd_plot <- renderPlot({
    tr <- treated(); p <- as.list(param(flu))
    e_bx <- input$emaxbx * tr$CBXE^input$hillbx /
            (tr$CBXE^input$hillbx + p$EC50_BX^input$hillbx + 1e-30)
    e_na <- input$emaxnai * tr$COCE / (tr$COCE + p$EC50_NAI + 1e-30)
    d <- rbind(data.frame(t = tr$time * 24, y = e_bx * 100,
                          op = "baloxavir: transcription + production"),
               data.frame(t = tr$time * 24, y = e_na * 100,
                          op = "NAI: release"))
    ggplot(d, aes(t, y, colour = op)) + geom_line(linewidth = 0.9) +
      coord_cartesian(ylim = c(0, 100)) +
      geom_vline(xintercept = t_rx() * 24, colour = "#1f77b4", linetype = 2) +
      labs(title = "Fractional blockade — what the concentrations buy",
           subtitle = paste("With a Michaelis (Hill = 1) slope the residual",
                            "production fraction cannot fall below EC50/(C+EC50),",
                            "\nand the model cannot reproduce the published 24-h fall",
                            "at any Emax. Move the Hill slider to see this."),
           x = "hours post-infection", y = "% blockade", colour = NULL) + THEME
  })

  output$pd_tbl <- renderTable({
    x <- FLU_pd_calibration()
    data.frame(quantity = names(x), value = unlist(x), row.names = NULL)
  }, digits = 4)

  ## ---- tab 6: clinical ---------------------------------------------------
  output$sym_plot <- renderPlot({
    o <- base_run(); tr <- treated(); p <- as.list(param(flu))
    d <- rbind(data.frame(t = o$time * 24,  y = o$SYM,  arm = "untreated"),
               data.frame(t = tr$time * 24, y = tr$SYM, arm = "treated"))
    ggplot(d, aes(t, y, colour = arm)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = p$SALLEV, linetype = 2, colour = "#2ca02c") +
      geom_vline(xintercept = t_rx() * 24, colour = "#1f77b4", linetype = 2) +
      annotate("text", x = 0, y = p$SALLEV + 0.6, hjust = 0, size = 3.4,
               colour = "#2ca02c",
               label = "alleviation ceiling: all 7 symptoms mild/absent") +
      scale_colour_manual(values = c(untreated = "#8a8a8a", treated = "#7b3fa0"),
                          name = NULL) +
      labs(title = "Composite symptom score (0–21) and the CAPSTONE rule",
           subtitle = "TTAS requires the score to stay under the ceiling for 21.5 h",
           x = "hours post-infection", y = "symptom score") + THEME
  })

  output$clin_plot <- renderPlot({
    tr <- treated()
    d <- rbind(data.frame(t = tr$time * 24, y = tr$TEMPC, v = "temperature (°C)"),
               data.frame(t = tr$time * 24, y = tr$SPO2,  v = "SpO2 (%)"),
               data.frame(t = tr$time * 24, y = tr$BAC * 10,
                          v = "bacterial burden (log10 CFU × 10)"))
    ggplot(d, aes(t, y, colour = v)) + geom_line(linewidth = 0.9) +
      facet_wrap(~v, scales = "free_y", ncol = 3) +
      labs(x = "hours post-infection", y = NULL, colour = NULL,
           title = "Fever, oxygenation and the bacterial sequel") +
      THEME + theme(legend.position = "none")
  })

  output$endpoint_tbl <- renderTable({
    r <- summarise_run(treated(), t_rx())
    data.frame(endpoint = c("TTAS (h)", "shedding cessation (h)",
                            "24-h titre change (log10)", "viral AUC (log10.d)",
                            "fever duration (h)", "peak symptom score",
                            "peak IL-6 (pg/mL)", "epithelium lost (%)",
                            "minimum SpO2 (%)", "peak bacterial load (log10)"),
               value = c(r$ttas_h, r$shed_h, r$drop24, r$auc_log, r$fever_h,
                         r$peak_sym, r$peak_il6, r$epi_lost * 100,
                         r$spo2_min, r$bac_peak))
  }, digits = 2)

  ## ---- tab 7: scenarios --------------------------------------------------
  output$scen_tbl <- renderTable(FLU_run_scenarios(), digits = 2)

  output$scen_plot <- renderPlot({
    ons <- onset(); trx <- ons + RX_DELAY_H / 24
    arms <- list(placebo = NULL, oseltamivir = rx_oseltamivir(trx),
                 baloxavir = rx_baloxavir(trx),
                 combination = rbind(rx_baloxavir(trx), rx_oseltamivir(trx)))
    d <- do.call(rbind, lapply(names(arms), function(nm) {
      o <- do.call(sim, c(list(doses = arms[[nm]], tmax = input$tmax),
                          overrides()))
      data.frame(t = o$time * 24, y = o$LOGV, arm = nm)
    }))
    ggplot(d, aes(t, y, colour = arm)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = LOD, linetype = 3, colour = "#999") +
      geom_vline(xintercept = trx * 24, colour = "#555", linetype = 2) +
      scale_colour_manual(values = PAL[names(arms)], name = NULL) +
      coord_cartesian(ylim = c(-2, 8)) +
      labs(title = "The shipped arms, all randomised at the same clock time",
           x = "hours post-infection", y = "log10 TCID50/mL") + THEME
  })

  ## ---- tab 8: resistance -------------------------------------------------
  release <- reactive({
    ons <- onset()
    hs  <- c(0, 12, 24, 36, 48, 72)
    ems <- c(0, 0.9, 0.99, 0.999, 0.9999, 0.999999)
    do.call(rbind, lapply(hs, function(h) {
      trx <- ons + h / 24
      do.call(rbind, lapply(ems, function(em) {
        o <- do.call(sim, c(list(doses = rx_baloxavir(trx), tmax = input$tmax,
                                 EMAX_BX = em, RF_BX = 50, RF_NAI = 1,
                                 COST = input$cost),
                            host_overrides()))
        data.frame(h = h, Emax = em, mut_peak = max(o$LOGVM[o$time >= trx]))
      }))
    }))
  })

  output$release_plot <- renderPlot({
    d <- release()
    ggplot(d, aes(factor(Emax), factor(h), fill = mut_peak)) +
      geom_tile(colour = "white") +
      geom_text(aes(label = sprintf("%.1f", mut_peak)), size = 3.2) +
      scale_fill_gradient(low = "#fff5e8", high = "#d97706",
                          name = "peak mutant\nlog10") +
      labs(title = "Competitive release over dose time × potency",
           subtitle = paste("Across a row: selection is not monotone in potency.",
                            "Down a column: release needs a field to release —",
                            "\ndosing late selects little, because the wild type",
                            "has already eaten the epithelium."),
           x = "baloxavir Emax", y = "hours after symptom onset at dosing") +
      THEME
  })

  output$release_tbl <- renderTable(FLU_host_comparison(), digits = 2)

  ## ---- tab 9: operators --------------------------------------------------
  output$oper_tbl <- renderTable(FLU_operator_decomposition(), digits = 2)

  output$oper_plot <- renderPlot({
    d <- FLU_operator_decomposition()
    d$operator <- factor(d$operator, levels = d$operator[order(d$AUC)])
    ggplot(d, aes(operator, AUC)) +
      geom_col(fill = "#7b3fa0", alpha = 0.85, width = 0.7) +
      coord_flip() +
      labs(title = "Viral AUC under each operator alone at 95%",
           subtitle = paste("Infected-cell death is the strongest operator in",
                            "the model — and no licensed antiviral uses it."),
           x = NULL, y = "viral AUC (log10.d)") + THEME
  })

  ## ---- tab 10: ledger ----------------------------------------------------
  output$ledger_tbl <- renderTable(FLU_trial_ledger(), digits = 2)

  output$discrepancies <- renderUI({
    HTML(paste0(
      "<h4>The five discrepancies this model reports rather than hides</h4><ol>",
      "<li><b>The Hill slope is forced by the data.</b> With a Michaelis ",
      "(h = 1) concentration–response the 24-hour titre fall saturates near ",
      "−3.3 log₁₀ at <i>any</i> E<sub>max</sub>; CAPSTONE-1 reports −4.8. The ",
      "trial datum is evidence about the <i>shape</i> of the response, not ",
      "only its potency. The calibrated in-vivo EC50 is also 13× below the ",
      "in-vitro value on a free-drug basis.</li>",
      "<li><b>Late therapy in the immunocompromised.</b> The model says it ",
      "keeps most of its value. That matches practice; no randomised trial ",
      "has tested it, so the prediction is unfalsified, not validated.</li>",
      "<li><b>The two drugs should not be equal on symptoms.</b> Any model in ",
      "which symptoms track virus predicts a larger benefit for baloxavir ",
      "than for oseltamivir. CAPSTONE-1 found 53.7 vs 53.8 h despite a 2-log ",
      "difference in day-2 titre.</li>",
      "<li><b>Resistance here is deterministic.</b> Real I38T emergence is a ",
      "stochastic establishment event in ~10% of treated adults; the model ",
      "gives an ensemble expectation that happens to nobody.</li>",
      "<li><b>The host ordering comes out backwards.</b> The model makes prior ",
      "immunity <i>increase</i> competitive release, because anything that ",
      "slows the wild type leaves a larger field. That is the reverse of the ",
      "children-vs-adults comparison, and it is the clearest sign that a ",
      "well-mixed epithelium is the wrong idealisation for this question.</li>",
      "</ol><p style='color:#666'>Educational and research use only. Not ",
      "validated for clinical decision-making.</p>"))
  })
}

shinyApp(ui, server)
