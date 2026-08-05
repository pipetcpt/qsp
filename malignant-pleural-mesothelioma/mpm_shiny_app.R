# =============================================================================
#  mpm_shiny_app.R
#  Interactive dashboard for the malignant pleural mesothelioma QSP model
#  (mpm_mrgsolve_model.R, 51 ODEs).
#
#  The app is organised around the one structural claim the model makes: the
#  tumour is a SHEET, so burden is AREA x THICKNESS and every clinical quantity
#  is downstream of which of the two factors a treatment moves.  Ten tabs:
#
#    1  Patient & tumour     — set the rind (area, thickness), histology,
#                              genotype, renal function; see the sheet drawn
#    2  Geometry             — A(t), h(t), V(t) side by side; the growth front
#    3  Delivery geometry    — penetration depth, exposed fraction, sanctuary,
#                              and the systemic/intrapleural complementarity
#    4  Drug PK              — pemetrexed (plasma + polyglutamates), platinum,
#                              mAbs with receptor occupancy, ADI/arginine
#    5  mRECIST vs truth     — the measured thickness sum against the viable
#                              cell count, and the fibrotic fraction between
#    6  Pleural space        — effusion, symphysis, VEGF, and the non-monotone
#                              relationship between effusion and tumour area
#    7  Toxicity             — ANC, CrCl, the nephrotoxicity feedback loop,
#                              homocysteine, irAE
#    8  Biomarker            — serum mesothelin with and without the renal
#                              confounder
#    9  Endpoints & survival — FVC, ECOG, hazard, survival, cohort medians
#   10  Scenario comparison  — any subset of the 24 built-in regimens, plus
#                              the trial calibration table
#
#  Run:  shiny::runApp("mpm_shiny_app.R")
# =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

source("mpm_mrgsolve_model.R", local = TRUE)

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

PAL <- c("#2C6E9B", "#C0392B", "#4F9C68", "#C9A227", "#8A76C4",
         "#C25C93", "#C2803F", "#4A91B5", "#7A8794", "#5B8C5A")

# -----------------------------------------------------------------------------
#  UI
# -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Malignant pleural mesothelioma — QSP model (area x thickness)"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("hist", "Histology",
                  c("epithelioid", "biphasic", "sarcomatoid"), "epithelioid"),
      sliderInput("h0", "Rind thickness at baseline (cm)",
                  0.15, 3.00, 0.80, 0.05),
      sliderInput("a0", "Involved pleural area (cm2)  [ceiling 1300]",
                  100, 1250, 450, 25),
      sliderInput("mfrac", "Fraction of the rind that is already matrix",
                  0.0, 0.6, 0.28, 0.02),
      sliderInput("crcl", "Creatinine clearance (mL/min)", 30, 130, 95, 5),

      h4("Genotype"),
      checkboxInput("ass1", "ASS1 promoter methylated (arginine auxotroph)", TRUE),
      checkboxInput("mtap", "CDKN2A/MTAP co-deleted (PRMT5i eligible)", TRUE),
      sliderInput("ercc1", "ERCC1 repair capacity (1 = normal)", 0.3, 2.5, 1.0, 0.1),

      h4("Regimen"),
      checkboxGroupInput(
        "drugs", NULL,
        c("Cisplatin" = "cis", "Pemetrexed" = "pem",
          "Folate + B12 supplementation" = "supp",
          "Bevacizumab 15 mg/kg q3w" = "bev",
          "Nivolumab 3 mg/kg q2w" = "nivo",
          "Ipilimumab 1 mg/kg q6w" = "ipi",
          "ADI-PEG20 36 mg/m2 weekly" = "adi",
          "PRMT5 inhibitor (continuous)" = "prmt5",
          "Intrapleural cisplatin 100 mg/m2" = "ip",
          "Talc pleurodesis (day 5)" = "talc",
          "Extended P/D (day 70)" = "epd",
          "Hemithoracic RT (day 110-145)" = "rt"),
        selected = c("cis", "pem", "supp")),
      sliderInput("ncyc", "Chemotherapy cycles", 0, 12, 6, 1),
      sliderInput("iodur", "Immunotherapy duration (days)", 0, 730, 730, 30),
      checkboxInput("adaptive", "Apply protocol dose delays / reductions", TRUE),
      sliderInput("tmax", "Horizon (days)", 200, 1250, 900, 50),
      actionButton("go", "Simulate", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient & tumour",
                 fluidRow(column(6, plotOutput("sheetPlot", height = 320)),
                          column(6, tableOutput("baseTbl"))),
                 helpText(HTML(
                   "The rectangle is the pleural surface of one hemithorax ",
                   "(~1300 cm<sup>2</sup>); the shaded part is involved, and the bar ",
                   "beneath it is the rind in cross-section, drawn to scale. ",
                   "Everything else in this app follows from these two numbers."))),

        tabPanel("2 · Geometry",
                 plotOutput("geomPlot", height = 520),
                 helpText(HTML(
                   "<b>Area and thickness are regulated separately.</b> The ",
                   "growth front is a perimeter, so dA/dt scales with ",
                   "2&radic;(&pi;A) and the volume doubling time lengthens on ",
                   "its own, with no resistance mechanism."))),

        tabPanel("3 · Delivery geometry",
                 fluidRow(column(7, plotOutput("penPlot", height = 340)),
                          column(5, plotOutput("sanctPlot", height = 340))),
                 plotOutput("routePlot", height = 260),
                 helpText(HTML(
                   "<b>The rind is fed from one face.</b> f_exp = (L/h)(1 - ",
                   "e<sup>-h/L</sup>). Intrapleural drug enters from the other ",
                   "face, so the two routes cover complementary depths. ",
                   "At the advancing margin (~0.6 mm) f_exp &rarr; 1, which is ",
                   "why chemotherapy stops spread long before it dents bulk."))),

        tabPanel("4 · Drug PK",
                 plotOutput("pkPlot", height = 560),
                 helpText(HTML(
                   "Pemetrexed polyglutamates are the pharmacologically active ",
                   "species and are retained for days; the plasma curve is over ",
                   "in hours. Receptor occupancy is shown for both checkpoints."))),

        tabPanel("5 · mRECIST vs truth",
                 plotOutput("recistPlot", height = 420),
                 fluidRow(column(6, plotOutput("fibPlot", height = 300)),
                          column(6, tableOutput("respTbl"))),
                 helpText(HTML(
                   "<b>mRECIST sums six thicknesses.</b> It sees h and is blind ",
                   "to A, and h contains the fibro-necrotic matrix left by the ",
                   "cells the drug already killed. A sheet that spreads thins ",
                   "itself, so the scan can improve while the tumour grows."))),

        tabPanel("6 · Pleural space",
                 plotOutput("pleuraPlot", height = 420),
                 plotOutput("effAreaPlot", height = 280),
                 helpText(HTML(
                   "Formation rises with involved area; stomatal drainage falls ",
                   "with it; the pleural space has a finite volume. The product ",
                   "is non-monotone &mdash; effusion peaks mid-disease and the ",
                   "late hemithorax goes dry."))),

        tabPanel("7 · Toxicity",
                 plotOutput("toxPlot", height = 520),
                 helpText(HTML(
                   "<b>The nephrotoxicity loop.</b> Cisplatin lowers CrCl, ",
                   "pemetrexed clearance is proportional to CrCl, so exposure ",
                   "rises and the next nadir is deeper. Folate raises the ",
                   "marrow EC50 far more than the tumour EC50."))),

        tabPanel("8 · Biomarker",
                 plotOutput("smrpPlot", height = 420),
                 helpText(HTML(
                   "Soluble mesothelin is a 40 kDa fragment cleared by ",
                   "glomerular filtration. The dashed line is the same tumour ",
                   "with renal function held at baseline: the gap is the part of ",
                   "the biomarker signal that is kidney, not cancer."))),

        tabPanel("9 · Endpoints & survival",
                 plotOutput("endPlot", height = 400),
                 plotOutput("survPlot", height = 320),
                 tableOutput("survTbl")),

        tabPanel("10 · Scenario comparison",
                 checkboxGroupInput("scens", "Regimens", names(SCEN),
                                    selected = names(SCEN)[c(1, 2, 4, 6, 7, 11)],
                                    inline = TRUE),
                 actionButton("goscen", "Run selected", class = "btn-primary"),
                 plotOutput("scenPlot", height = 420),
                 tableOutput("scenTbl"),
                 h4("Trial calibration"),
                 tableOutput("calTbl"),
                 helpText(HTML(
                   "Control arms differ by more than four months across these ",
                   "trials, so the model is anchored on three absolute medians ",
                   "and otherwise calibrated to <b>within-trial hazard ",
                   "ratios</b>.")))
      )
    )
  )
)

# -----------------------------------------------------------------------------
#  SERVER
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  build <- eventReactive(input$go, {
    d <- input$drugs
    pars <- list(H0_RIND = input$h0, A0_RIND = input$a0, MFRAC0 = input$mfrac,
                 CRCL0 = input$crcl, ERCC1F = input$ercc1,
                 ASS1LOSS = as.numeric(input$ass1))
    ev_all <- NULL
    add <- function(e) if (is.null(ev_all)) e else ev_all + e

    if (("cis" %in% d || "pem" %in% d) && input$ncyc > 0)
      ev_all <- add(chemo_events(input$ncyc, pem = "pem" %in% d,
                                 plat = if ("cis" %in% d) "cis" else "none"))
    if ("bev" %in% d)  ev_all <- add(bev_events(input$tmax))
    if ("nivo" %in% d || "ipi" %in% d)
      ev_all <- add(io_events(input$iodur, nivo = "nivo" %in% d,
                              ipi = "ipi" %in% d))
    if ("adi" %in% d)  ev_all <- add(adi_events(min(input$tmax, 730)))
    if ("ip"  %in% d)  ev_all <- add(ip_events(max(input$ncyc, 1)))

    wins <- NULL
    if ("talc"  %in% d) wins <- rbind(wins, data.frame(name = "TALCON", t0 = 5,
                                                       t1 = 7, value = 1))
    if ("rt"    %in% d) wins <- rbind(wins, data.frame(name = "RTRATE", t0 = 110,
                                                       t1 = 145, value = 0.055))
    if ("prmt5" %in% d && input$mtap)
      wins <- rbind(wins, data.frame(name = "PRMT5", t0 = 0, t1 = input$tmax,
                                     value = 0.0075))
    surg <- if ("epd" %in% d)
      list(time = 70, h_res = 0.020, mort = 0.040, fvc_loss = 0.14,
           hazx = 1.30, lbm_loss = 0.10) else NULL

    sim <- sim_scenario(input$hist, ev_all, pars, surgery = surg,
                        windows = wins, tmax = input$tmax,
                        supplemented = "supp" %in% d)

    # counterfactual with renal function pinned, for the biomarker tab
    sim_kid <- sim_scenario(input$hist, ev_all,
                            modifyList(pars, list(KNEPH = 0, KNEPHIRR = 0)),
                            surgery = surg, windows = wins, tmax = input$tmax,
                            supplemented = "supp" %in% d)
    list(sim = sim, kid = sim_kid)
  }, ignoreNULL = FALSE)

  S <- reactive(build()$sim)

  # ---- tab 1 ---------------------------------------------------------------
  output$sheetPlot <- renderPlot({
    a <- input$a0; h <- input$h0; Spl <- 1300
    df <- data.frame(x = c(0, 1), y = c(0, 1))
    ggplot() +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0.45, ymax = 1,
               fill = "grey92", colour = "grey60") +
      annotate("rect", xmin = 0, xmax = a / Spl, ymin = 0.45, ymax = 1,
               fill = PAL[3], alpha = 0.65) +
      annotate("text", x = 0.5, y = 1.06,
               label = sprintf("pleural surface 1300 cm2 — involved %.0f cm2 (%.0f %%)",
                               a, 100 * a / Spl)) +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0.10,
               ymax = 0.10 + 0.28 * min(h, 3) / 3, fill = PAL[2], alpha = 0.7) +
      annotate("text", x = 0.5, y = 0.03,
               label = sprintf("rind thickness %.2f cm  ->  volume %.0f cm3",
                               h, a * h)) +
      coord_cartesian(xlim = c(0, 1), ylim = c(-0.02, 1.12)) +
      theme_void()
  })

  output$baseTbl <- renderTable({
    s <- S()[1, ]
    data.frame(
      quantity = c("involved area (cm2)", "measured thickness (cm)",
                   "true volume (cm3)", "viable tumour (cm3)",
                   "matrix fraction of the rind (%)",
                   "mRECIST 6-site sum (mm)", "penetration depth L_p (mm)",
                   "exposed fraction f_exp", "geometric sanctuary",
                   "effusion (mL)", "serum mesothelin (nM)",
                   "FVC (% predicted)", "ECOG"),
      value = signif(c(s$A, s$hmeas, s$VTUM, s$N, 100 * s$FIBFR, s$MRECIST,
                       s$LPMM, s$FEXP, s$SANCT, s$VEFF, s$SMRP, s$FVCP,
                       s$ECOG), 3))
  })

  # ---- tab 2 ---------------------------------------------------------------
  output$geomPlot <- renderPlot({
    s <- S()
    d <- s %>% transmute(time,
                         `involved area A (cm2)` = A,
                         `measured thickness h (cm)` = hmeas,
                         `viable thickness (cm)` = hvia,
                         `burden V = A x h (cm3)` = VTUM,
                         `viable tumour N (cm3)` = N,
                         `growth front perimeter (cm)` = 2 * sqrt(pi * A)) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(colour = PAL[1], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") + labs(x = "day", y = NULL) + THEME
  })

  # ---- tab 3 ---------------------------------------------------------------
  output$penPlot <- renderPlot({
    s <- S()
    d <- s %>% transmute(time, `L_p (mm)` = LPMM,
                         `exposed fraction` = FEXP,
                         `sanctuary fraction` = SANCT) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL) + labs(x = "day", y = NULL,
                                               colour = NULL) + THEME
  })

  output$sanctPlot <- renderPlot({
    hs <- seq(0.1, 3.2, by = 0.05)
    fp <- function(L, h) (L / h) * (1 - exp(-h / L))
    d <- bind_rows(
      data.frame(h = hs, f = fp(0.075, hs), route = "systemic small molecule"),
      data.frame(h = hs, f = fp(0.085, hs), route = "intrapleural"),
      data.frame(h = hs, f = fp(0.020, hs), route = "150 kDa IgG"),
      data.frame(h = hs, f = fp(0.030, hs), route = "T cells"))
    ggplot(d, aes(h, f, colour = route)) + geom_line(linewidth = 0.9) +
      geom_vline(xintercept = input$h0, linetype = 2, colour = "grey40") +
      scale_colour_manual(values = PAL) +
      labs(x = "rind thickness h (cm)", y = "exposed fraction f_exp",
           colour = NULL) + THEME
  })

  output$routePlot <- renderPlot({
    hs <- seq(0.1, 3.2, by = 0.05)
    fp <- function(L, h) (L / h) * (1 - exp(-h / L))
    fs <- fp(0.075, hs); fi <- fp(0.085, hs)
    d <- data.frame(h = hs,
                    systemic = fs,
                    `intrapleural adds` = pmax(0, pmin(fi, 1 - fs)),
                    sanctuary = pmax(0, 1 - fs - fi), check.names = FALSE) %>%
      pivot_longer(-h)
    ggplot(d, aes(h, value, fill = name)) + geom_area() +
      scale_fill_manual(values = c(PAL[1], PAL[4], "grey85")) +
      labs(x = "rind thickness h (cm)", y = "fraction of the rind",
           fill = NULL,
           title = "The two routes attack opposite faces, so their coverage adds") +
      THEME
  })

  # ---- tab 4 ---------------------------------------------------------------
  output$pkPlot <- renderPlot({
    s <- S()
    d <- s %>% transmute(time,
                         `pemetrexed plasma (mg/L)` = CPEM,
                         `tumour polyglutamates (mg/L)` = PEM_TP,
                         `marrow polyglutamates (mg/L)` = PEM_MP,
                         `free platinum (mg/L)` = CCIS,
                         `Pt-DNA adducts` = ADD,
                         `bevacizumab (mg/L)` = CBEV,
                         `PD-1 occupancy` = ROPD1,
                         `CTLA-4 occupancy` = ROCT4,
                         `plasma arginine (uM)` = ARG) %>% pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(colour = PAL[4], linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") + labs(x = "day", y = NULL) + THEME
  })

  # ---- tab 5 ---------------------------------------------------------------
  output$recistPlot <- renderPlot({
    s <- S()
    d <- data.frame(
      time = s$time,
      `mRECIST thickness sum` = 100 * (s$MRECIST - s$MRECIST[1]) / s$MRECIST[1],
      `viable tumour cells` = 100 * (s$N - s$N[1]) / s$N[1],
      `true volume A x h` = 100 * (s$VTUM - s$VTUM[1]) / s$VTUM[1],
      `involved area` = 100 * (s$A - s$A[1]) / s$A[1],
      check.names = FALSE) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = c(-30, 20), linetype = 2, colour = "grey55") +
      annotate("text", x = Inf, y = -30, hjust = 1.05, vjust = -0.4,
               label = "PR threshold", size = 3, colour = "grey40") +
      annotate("text", x = Inf, y = 20, hjust = 1.05, vjust = -0.4,
               label = "PD threshold", size = 3, colour = "grey40") +
      scale_colour_manual(values = PAL) +
      labs(x = "day", y = "% change from baseline", colour = NULL) + THEME
  })

  output$fibPlot <- renderPlot({
    s <- S()
    ggplot(data.frame(time = s$time, f = 100 * s$FIBFR),
           aes(time, f)) + geom_area(fill = PAL[7], alpha = 0.6) +
      labs(x = "day", y = "% of the MEASURED rind that is matrix, not tumour") +
      THEME
  })

  output$respTbl <- renderTable({
    s <- S(); n <- min(300, nrow(s))
    best_m <- min(100 * (s$MRECIST[1:n] - s$MRECIST[1]) / s$MRECIST[1])
    i <- which.min(s$MRECIST[1:n])
    data.frame(
      quantity = c("best mRECIST change (%)", "viable-cell change at that day (%)",
                   "gap (percentage points)", "mRECIST category",
                   "day of best response"),
      value = c(sprintf("%.1f", best_m),
                sprintf("%.1f", 100 * (s$N[i] - s$N[1]) / s$N[1]),
                sprintf("%.1f", 100 * (s$N[i] - s$N[1]) / s$N[1] - best_m),
                if (best_m <= -30) "PR" else if (best_m < 20) "SD" else "PD",
                sprintf("%d", s$time[i])))
  })

  # ---- tab 6 ---------------------------------------------------------------
  output$pleuraPlot <- renderPlot({
    s <- S()
    d <- s %>% transmute(time, `effusion (mL)` = VEFF, `symphysis` = PSY,
                         `free VEGF (pg/mL)` = VEGFF,
                         `interstitial pressure (mmHg)` = IFPo) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(colour = PAL[8], linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") + labs(x = "day", y = NULL) + THEME
  })

  output$effAreaPlot <- renderPlot({
    s <- S()
    ggplot(data.frame(A = s$A, V = s$VEFF, t = s$time),
           aes(A, V, colour = t)) + geom_path(linewidth = 1) +
      scale_colour_viridis_c(name = "day") +
      labs(x = "involved pleural area (cm2)", y = "effusion volume (mL)",
           title = "Effusion is non-monotone in tumour area") + THEME
  })

  # ---- tab 7 ---------------------------------------------------------------
  output$toxPlot <- renderPlot({
    s <- S()
    d <- s %>% transmute(time, `ANC (10^9/L)` = CIRC,
                         `creatinine clearance (mL/min)` = CRCL,
                         `homocysteine (uM)` = HCY,
                         `plasma folate (nM)` = FOL,
                         `cumulative pemetrexed AUC` = AUCP,
                         `irAE intensity` = IRAE) %>% pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(colour = PAL[2], linewidth = 0.85) +
      geom_hline(data = data.frame(name = "ANC (10^9/L)", y = 0.5),
                 aes(yintercept = y), linetype = 2, colour = "grey50") +
      facet_wrap(~name, scales = "free_y") + labs(x = "day", y = NULL) + THEME
  })

  # ---- tab 8 ---------------------------------------------------------------
  output$smrpPlot <- renderPlot({
    b <- build()
    d <- data.frame(time = b$sim$time,
                    `observed SMRP` = b$sim$SMRP,
                    `SMRP if CrCl were held at baseline` = b$kid$SMRP,
                    `viable tumour (scaled)` = b$sim$N / b$sim$N[1] * b$sim$SMRP[1],
                    check.names = FALSE) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name, linetype = name)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(2, 1, 3)]) +
      scale_linetype_manual(values = c(1, 2, 3)) +
      labs(x = "day", y = "serum mesothelin (nM)", colour = NULL,
           linetype = NULL) + THEME
  })

  # ---- tab 9 ---------------------------------------------------------------
  output$endPlot <- renderPlot({
    s <- S()
    d <- s %>% transmute(time, `FVC (% predicted)` = FVCP, `ECOG` = ECOG,
                         `pain index` = PAINi, `lean body mass (kg)` = LBM,
                         `IL-6 (pg/mL)` = IL6, `invasion depth (cm)` = Z) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(colour = PAL[6], linewidth = 0.85) +
      facet_wrap(~name, scales = "free_y") + labs(x = "day", y = NULL) + THEME
  })

  output$survPlot <- renderPlot({
    s <- S()
    ggplot(data.frame(t = s$time / 30.44, S = s$SURV), aes(t, S)) +
      geom_line(linewidth = 1, colour = PAL[1]) +
      geom_hline(yintercept = 0.5, linetype = 2, colour = "grey50") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "months", y = "survival probability") + THEME
  })

  output$survTbl <- renderTable({
    cc <- cohort_curve(function(h) {
      isolate({
        d <- input$drugs
        pars <- list(H0_RIND = input$h0, A0_RIND = input$a0,
                     MFRAC0 = input$mfrac, CRCL0 = input$crcl,
                     ERCC1F = input$ercc1, ASS1LOSS = as.numeric(input$ass1))
        ev_all <- if (input$ncyc > 0 && ("cis" %in% d || "pem" %in% d))
          chemo_events(input$ncyc, pem = "pem" %in% d) else NULL
        sim_scenario(h, ev_all, pars, tmax = 1250,
                     supplemented = "supp" %in% d)
      })
    })
    data.frame(readout = c("median OS, this patient (mo)",
                           "median OS, trial-mix cohort (mo)"),
               value = c(round(median_from_S(S()$time, S()$SURV) / 30.44, 1),
                         round(cc$median_mo, 1)))
  })

  # ---- tab 10 --------------------------------------------------------------
  scenres <- eventReactive(input$goscen, {
    lapply(input$scens, function(nm) {
      s <- SCEN[[nm]](input$hist); s$scenario <- nm; s
    }) %>% bind_rows()
  })

  output$scenPlot <- renderPlot({
    d <- scenres()
    dd <- d %>% group_by(scenario) %>%
      mutate(mR = 100 * (MRECIST - first(MRECIST)) / first(MRECIST)) %>%
      ungroup() %>%
      select(time, scenario, `survival` = SURV, `mRECIST % change` = mR,
             `involved area (cm2)` = A, `viable tumour (cm3)` = N) %>%
      pivot_longer(-c(time, scenario))
    ggplot(dd, aes(time, value, colour = scenario)) +
      geom_line(linewidth = 0.8) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = rep(PAL, 3)) +
      labs(x = "day", y = NULL, colour = NULL) + THEME +
      theme(legend.text = element_text(size = 8))
  })

  output$scenTbl <- renderTable({
    d <- scenres()
    d %>% group_by(scenario) %>%
      summarise(`median OS (mo)` = round(median_from_S(time, SURV) / 30.44, 1),
                `best mRECIST %` = round(min(100 * (MRECIST - first(MRECIST)) /
                                               first(MRECIST)), 1),
                `viable nadir %` = round(min(100 * (N - first(N)) / first(N)), 1),
                `ANC nadir` = round(min(CIRC), 2),
                `peak effusion (mL)` = round(max(VEFF)),
                `FVC at 6 mo` = round(FVCP[which.min(abs(time - 183))], 1),
                .groups = "drop")
  })

  output$calTbl <- renderTable({
    data.frame(
      comparison = c("best supportive care (absolute)",
                     "cisplatin alone (absolute)",
                     "cisplatin + pemetrexed (absolute)",
                     "cis+pem vs cisplatin", "+ bevacizumab",
                     "nivolumab+ipilimumab vs chemo",
                     "  epithelioid", "  non-epithelioid",
                     "pembrolizumab + chemo vs chemo",
                     "nivolumab 2L vs placebo",
                     "pegargiminase + chemo (non-epithelioid)",
                     "extended P/D + chemo vs chemo",
                     "extrapleural pneumonectomy vs none"),
      trial = c("MS01", "EMPHACIS", "EMPHACIS", "EMPHACIS", "MAPS",
                "CheckMate 743", "CheckMate 743", "CheckMate 743",
                "IND.227/KN-483", "CONFIRM", "ATOMIC-Meso", "MARS2", "MARS"),
      observed = c("7.6 mo", "9.3 mo", "12.1 mo", "HR 0.77", "HR 0.77",
                   "HR 0.74", "HR 0.86", "HR 0.46", "HR 0.79", "HR 0.69",
                   "HR 0.71", "HR 1.28", "HR 1.90"))
  })
}

shinyApp(ui, server)
