# =============================================================================
# cgd_shiny_app.R
# Chronic Granulomatous Disease (만성 육아종병) — interactive QSP dashboard
#
# Twelve tabs.  The app is organised around ONE argument, and the tabs are
# ordered so that the argument builds:
#
#    1  Patient          — set the two numbers that ARE the disease: f and r
#    2  DHR & histogram  — why those two numbers are not one number
#    3  Phagosome K(phi) — the kernel's output, and its lack of a threshold
#    4  Neutrophil PK    — granulopoiesis and traffic
#    5  Drug PK          — seven drugs, correct units
#    6  Bacterial focus  — including the SHELTERED compartment
#    7  Fungal focus     — conidia vs hyphae, the asymmetry that kills
#    8  Critical inoculum— the deterministic boundary
#    9  Infection rates  — exposure crossing that boundary; the trial arms
#   10  Inflammation     — granuloma and colitis WITHOUT an organism
#   11  Correction       — chimerism / gene therapy threshold
#   12  Scenarios        — sixteen side-by-side comparisons
#
# Calibrated constants are machine-written by sync_r_params.py from
# cgd_calibration.json; do not edit them by hand.
#
# NOTE: no R toolchain was available where this file was written, so it has not
# been run.  It mirrors cgd_python_reference.py, which has been.
#
#   library(shiny); library(mrgsolve); library(dplyr); library(ggplot2)
#   shiny::runApp("cgd_shiny_app.R")
# =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---- BEGIN AUTO-SYNC (written by sync_r_params.py) -----------------------
CAL <- list(
  K0        = 0.221849,   # log10 kill at phi = 0 (the whole non-oxidative arm)
  Kmax      = 1.301030,   # log10 kill at phi = 1
  phi50     = 0.544822,   # Hill midpoint: NOT a threshold, it is above 0.5
  hillK     = 2.735060,
  lambda_b  = 9.216740,   # bacterial exposures / patient-year
  muN_b     = 2.863577,   # log10 mean bacterial inoculum
  lambda_f  = 0.201803,   # fungal exposures / patient-year
  muN_f     = 4.470990,   # log10 mean conidial inoculum
  sdN       = 1.000000,
  CFR_B     = 0.010000,   # case fatality, bacterial (FIXED, not fitted)
  CFR_F     = 0.200000,   # case fatality, invasive aspergillosis (FIXED)
  Ncrit_cgd = 3.683441,   # log10, untreated X-CGD null
  Ncrit_hlt = 5.715912,   # log10, healthy control
  E_TMP     = 3.188810,   # /d bacterial kill at steady-state co-trimoxazole
  E_ITZ     = 1.808111,   # /d hyphal kill at steady-state itraconazole
  E_IFN     = 1.785696,   # fold boost of the oxidase-independent arm
  dhr_xover = 0.238099    # DHR mean where mosaic/uniform ordering flips
)
# ---- END AUTO-SYNC -------------------------------------------------------

mod <- mread("cgd_mrgsolve_model", ".", quiet = TRUE)

# ---------------------------------------------------------------- helpers
K_hill <- function(phi) {
  p <- pmax(phi, 0)
  CAL$K0 + (CAL$Kmax - CAL$K0) * p^CAL$hillK /
    (CAL$phi50^CAL$hillK + p^CAL$hillK + 1e-30)
}
surv_one <- function(phi, boost = 1) 10^(-(K_hill(phi) + (boost - 1) * CAL$K0))

# The mixture.  Averaged in survival space, never in log-kill space: averaging
# log kills is what makes a carrier and a hypomorph look like the same patient.
s_enc <- function(f, r, boost = 1) f * surv_one(1, boost) + (1 - f) * surv_one(r, boost)

dhr_mean <- function(f, r) f + (1 - f) * r

infection_rate <- function(logNcrit, lambda, muN, sdn = CAL$sdN) {
  if (!is.finite(logNcrit)) return(if (logNcrit > 0) 0 else lambda)
  lambda * pnorm((logNcrit - muN) / sdn, lower.tail = FALSE)
}

GENOTYPES <- data.frame(
  name     = c("healthy", "X-CGD null (gp91phox-)", "X-CGD hypomorph 5%",
               "X-CGD hypomorph 20%", "p47phox-deficient (AR)",
               "p67phox-deficient (AR)", "X-CGD carrier 50%",
               "X-CGD carrier 20%", "X-CGD carrier 5%"),
  f_normal = c(1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 0.20, 0.05),
  phi_res  = c(1.00, 0.00, 0.05, 0.20, 0.12, 0.03, 0.00, 0.00, 0.00),
  stringsAsFactors = FALSE
)

cgd_regimen <- function(days = 365, wt = 30, bsa = 1.1, tmpsmx = FALSE,
                        itra = FALSE, vori = FALSE, ifng = FALSE,
                        pred = FALSE, pred_dose = 1.0, anakinra = FALSE) {
  evs <- list()
  add <- function(cmt, amt, ii, n, start = 0) {
    evs[[length(evs) + 1]] <<- ev(time = start, cmt = cmt, amt = amt,
                                  ii = ii, addl = max(n - 1, 0))
  }
  if (tmpsmx) { add("TMP_g", 5 * wt / 2, 0.5, days * 2)
                add("SMX_g", 25 * wt / 2, 0.5, days * 2) }
  if (itra)   add("ITZ_g", if (wt < 40) 100 else 200, 0.5, days * 2)
  if (vori)   add("VOR_g", 9 * wt, 0.5, days * 2)
  if (ifng) for (w in 0:floor(days / 7)) for (off in c(0, 2, 4)) {
    tt <- w * 7 + off
    if (tt <= days) add("IFN_sc", 50 * bsa, 1e9, 1, start = tt)
  }
  if (pred)     add("PRED_g", pred_dose * wt, 1, days)
  if (anakinra) add("ANA_sc", 2 * wt, 1, days)
  if (!length(evs)) return(ev(time = 0, cmt = "TMP_g", amt = 0))
  Reduce(`+`, evs)
}

run_patient <- function(f, r, days, regimen, chim = 0, vcn = 0,
                        seedB = 0, seedC = 0, extra_hz = 0, ido_nox = 0) {
  m <- mod %>% param(f_normal = f, phi_res = r, chim_tgt = chim,
                     vcn_tgt = vcn, extra_hz = extra_hz, ido_nox = ido_nox)
  if (seedB > 0) m <- m %>% init(BACT = seedB)
  if (seedC > 0) m <- m %>% init(CONID = seedC)
  m %>% ev(regimen) %>% mrgsim(end = days, delta = 1) %>% as.data.frame()
}

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95"),
        legend.position = "bottom")

# =============================================================================
# UI
# =============================================================================
ui <- fluidPage(
  titlePanel("만성 육아종병 (Chronic Granulomatous Disease) — QSP 대시보드"),
  tags$p(style = "color:#666;",
         paste("CGD is a disease of one number: the electron flux phi through NOX2.",
               "Everything on these tabs is downstream of the two sliders at the",
               "top left — the fraction of neutrophils with a working oxidase, and",
               "the residual activity of the ones without.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      selectInput("geno", "Genotype preset", choices = GENOTYPES$name,
                  selected = "X-CGD null (gp91phox-)"),
      sliderInput("f_normal", "f — fraction of oxidase-COMPETENT neutrophils",
                  0, 1, 0, step = 0.01),
      sliderInput("phi_res", "r — residual activity of the REST",
                  0, 1, 0, step = 0.01),
      helpText(htmlOutput("dhr_note")),
      hr(),
      h4("교정 (Correction)"),
      sliderInput("chim", "HSCT donor myeloid chimerism", 0, 1, 0, step = 0.01),
      sliderInput("vcn", "Gene-corrected myeloid fraction", 0, 1, 0, step = 0.01),
      hr(),
      h4("약물 (Therapy)"),
      checkboxInput("tmpsmx", "Co-trimoxazole 5 mg/kg/d (TMP)", FALSE),
      checkboxInput("itra", "Itraconazole 200 mg BID", FALSE),
      checkboxInput("vori", "Voriconazole 9 mg/kg BID", FALSE),
      checkboxInput("ifng", "Interferon gamma-1b 50 ug/m2 3x/wk", FALSE),
      checkboxInput("pred", "Prednisolone 1 mg/kg/d", FALSE),
      checkboxInput("anak", "Anakinra 2 mg/kg/d", FALSE),
      hr(),
      h4("노출 (Challenge)"),
      numericInput("seedB", "Bacterial inoculum (CFU)", 0, min = 0),
      numericInput("seedC", "Conidial inoculum", 0, min = 0),
      sliderInput("days", "Simulation days", 30, 730, 365, step = 5),
      checkboxInput("ido_nox", "IDO requires NOX2 (Romani 2008; not replicated)",
                    FALSE),
      hr(),
      helpText(HTML(paste0(
        "<b>Calibration provenance.</b> Five numbers were fitted to NORMAL ",
        "neutrophils in the phagosome kernel, and four to exposure ",
        "epidemiology. Nothing else is fitted. Trial effects on these tabs ",
        "are predictions.")))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · 환자 프로파일",
                 h4("The two numbers, and what they imply"),
                 tableOutput("profile"),
                 plotOutput("profile_plot", height = "330px"),
                 htmlOutput("profile_text")),
        tabPanel("2 · DHR 히스토그램",
                 h4("Why one number is not enough"),
                 plotOutput("dhr_hist", height = "330px"),
                 plotOutput("mosaic_vs_uniform", height = "330px"),
                 htmlOutput("dhr_text")),
        tabPanel("3 · 식포 K(phi)",
                 h4("The phagosome kernel's only output"),
                 plotOutput("kphi", height = "380px"),
                 htmlOutput("kphi_text")),
        tabPanel("4 · 호중구 동태",
                 h4("Granulopoiesis, blood ANC, and traffic to the focus"),
                 plotOutput("pk_neut", height = "480px")),
        tabPanel("5 · 약물 PK",
                 h4("Seven drugs. Amounts mg, volumes L, concentrations mg/L."),
                 plotOutput("pk_drug", height = "520px"),
                 tableOutput("pk_summary")),
        tabPanel("6 · 세균 병소",
                 h4("The focus, and the sheltered compartment"),
                 plotOutput("focus_b", height = "460px"),
                 htmlOutput("focus_b_text")),
        tabPanel("7 · 진균 병소",
                 h4("Conidia are ingestible. Hyphae are not."),
                 plotOutput("focus_f", height = "460px"),
                 htmlOutput("focus_f_text")),
        tabPanel("8 · 임계 접종량",
                 h4("The deterministic clearance / disease boundary"),
                 plotOutput("ncrit", height = "400px"),
                 tableOutput("ncrit_tbl")),
        tabPanel("9 · 감염률과 임상시험",
                 h4("Exposure crossing that boundary"),
                 tableOutput("rates_tbl"),
                 plotOutput("trials", height = "360px"),
                 htmlOutput("trials_text")),
        tabPanel("10 · 염증·육아종",
                 h4("Granuloma and colitis in a patient with no organism"),
                 plotOutput("inflam", height = "480px"),
                 htmlOutput("inflam_text")),
        tabPanel("11 · 교정 역치",
                 h4("How much corrected myelopoiesis is enough?"),
                 plotOutput("threshold", height = "420px"),
                 htmlOutput("threshold_text")),
        tabPanel("12 · 시나리오 비교",
                 h4("Sixteen scenarios side by side"),
                 tableOutput("scen_tbl"),
                 plotOutput("scen_plot", height = "420px"))
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  observeEvent(input$geno, {
    g <- GENOTYPES[GENOTYPES$name == input$geno, ]
    updateSliderInput(session, "f_normal", value = g$f_normal)
    updateSliderInput(session, "phi_res", value = g$phi_res)
  })

  regimen <- reactive(cgd_regimen(days = input$days, tmpsmx = input$tmpsmx,
                                  itra = input$itra, vori = input$vori,
                                  ifng = input$ifng, pred = input$pred,
                                  anakinra = input$anak))

  sim <- reactive({
    run_patient(input$f_normal, input$phi_res, input$days, regimen(),
                chim = input$chim, vcn = input$vcn,
                seedB = input$seedB, seedC = input$seedC,
                ido_nox = as.numeric(input$ido_nox))
  })

  fx <- reactive(min(1, input$f_normal + input$chim + input$vcn))
  boost <- reactive(if (input$ifng) CAL$E_IFN else 1)

  output$dhr_note <- renderText({
    d <- dhr_mean(fx(), input$phi_res)
    sprintf("<b>DHR mean = %.3f</b><br>The flow lab reports this one number.
             Two very different marrows produce it.", d)
  })

  # ---- 1 -------------------------------------------------------------------
  output$profile <- renderTable({
    d <- dhr_mean(fx(), input$phi_res)
    data.frame(
      Quantity = c("f (oxidase-competent fraction)", "r (residual activity)",
                   "DHR-123 mean", "survival per phagosome, competent cell",
                   "survival per phagosome, mutant cell",
                   "population survival per encounter s_enc",
                   "same for a UNIFORM patient at the same DHR mean"),
      Value = sprintf("%.4f", c(fx(), input$phi_res, d, surv_one(1, boost()),
                                surv_one(input$phi_res, boost()),
                                s_enc(fx(), input$phi_res, boost()),
                                surv_one(d, boost()))))
  })
  output$profile_plot <- renderPlot({
    s <- sim()
    s %>% select(time, ANC, CRP, IL1B, GRAN, COL, DHRmean) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.8, colour = "#2E75B6") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL) + THEME
  })
  output$profile_text <- renderUI(HTML(
    "<p><b>Read the two sliders as the whole disease.</b> Every genotype,
     every carrier, every transplant recipient and every gene-therapy patient
     in this model is a point in the (f, r) plane. Nothing else about the
     patient enters the killing equations.</p>"))

  # ---- 2 -------------------------------------------------------------------
  output$dhr_hist <- renderPlot({
    d <- dhr_mean(fx(), input$phi_res)
    df <- rbind(
      data.frame(patient = "this patient (mixture)",
                 phi = c(rep(1, round(fx() * 1000)),
                         rep(input$phi_res, 1000 - round(fx() * 1000)))),
      data.frame(patient = sprintf("uniform patient, same DHR mean %.2f", d),
                 phi = rep(d, 1000)))
    ggplot(df, aes(phi)) +
      geom_histogram(bins = 50, fill = "#D9822B", colour = NA) +
      facet_wrap(~patient, ncol = 1, scales = "free_y") +
      labs(x = "per-cell oxidase activity (what DHR would show cell by cell)",
           y = "cells", title = "The histogram, not the mean") + THEME
  })
  output$mosaic_vs_uniform <- renderPlot({
    m <- seq(0.01, 0.99, by = 0.01)
    df <- data.frame(
      dhr = rep(m, 2),
      surv = c(sapply(m, function(x) s_enc(x, 0, boost())),
               sapply(m, function(x) surv_one(x, boost()))),
      arrangement = rep(c("mosaic (carrier / chimerism / gene therapy)",
                          "uniform (hypomorphic mutation)"), each = length(m)))
    ggplot(df, aes(dhr, surv, colour = arrangement)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = CAL$dhr_xover, linetype = 2) +
      annotate("text", x = CAL$dhr_xover, y = max(df$surv), hjust = -0.05,
               label = sprintf("ordering flips at DHR %.3f", CAL$dhr_xover),
               size = 3.4) +
      labs(x = "DHR mean (the number the lab reports)",
           y = "organism survival per phagocytic encounter",
           colour = NULL) + THEME
  })
  output$dhr_text <- renderUI(HTML(sprintf(
    "<p><b>Below a DHR mean of %.3f, mosaicism is the better arrangement;
     above it, uniform residual activity is.</b> The crossover is the
     inflection of K(phi) and nothing else. It matters because the clinic
     operates below it: carriers, mixed chimerism and gene-therapy marking all
     sit in the 5-25%% band. A flow report of '20%% of normal' does not say
     which of the two patients it is describing unless the HISTOGRAM —
     bimodal versus uniformly shifted — is read as well.</p>", CAL$dhr_xover)))

  # ---- 3 -------------------------------------------------------------------
  output$kphi <- renderPlot({
    p <- seq(0, 1, by = 0.005)
    df <- data.frame(phi = p, K = K_hill(p), surv = 10^(-K_hill(p)))
    ggplot(df, aes(phi, K)) + geom_line(linewidth = 1, colour = "#C0504D") +
      geom_vline(xintercept = c(fx(), input$phi_res), linetype = 3) +
      geom_hline(yintercept = CAL$K0, linetype = 2, colour = "#888888") +
      annotate("text", x = 0.02, y = CAL$K0, vjust = -0.6, hjust = 0, size = 3.4,
               label = "K0 — the entire non-oxidative arm") +
      labs(x = "phi (fractional NADPH oxidase activity of ONE neutrophil)",
           y = "K(phi) = log10 kill of one ingested organism in 60 min",
           title = sprintf("phi50 = %.3f, Hill = %.2f — the midpoint is ABOVE half",
                           CAL$phi50, CAL$hillK)) + THEME
  })
  output$kphi_text <- renderUI(HTML(
    "<p><b>There is no threshold inside the phagosome.</b> The curve is convex
     over the clinically relevant range and the first 20% of oxidase activity
     buys the least. That is the opposite of the intuition behind gene therapy
     and mixed chimerism, where 10-20% correction is known to be enough — so
     the clinical threshold must be manufactured somewhere else. It is: by the
     sheltered compartment on tab 6, which correcting a cell REMOVES as well as
     adding a killer.</p>"))

  # ---- 4 -------------------------------------------------------------------
  output$pk_neut <- renderPlot({
    sim() %>% select(time, PROG, RES, ANC, GCSF, MONO, NEUT_T, MAC, APOP) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#41729F", linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") + labs(x = "day", y = NULL) + THEME
  })

  # ---- 5 -------------------------------------------------------------------
  output$pk_drug <- renderPlot({
    sim() %>% select(time, CTMPo, CITZo, COHIo, CVORo, CIFNo, CANAo, CPREDo) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#3C8C3C", linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = "concentration (mg/L, or ng/mL for IFN/ANA/PRED)") + THEME
  })
  output$pk_summary <- renderTable({
    s <- sim(); n <- nrow(s)
    data.frame(Drug = c("trimethoprim (free)", "itraconazole", "OH-itraconazole",
                        "voriconazole", "interferon gamma-1b", "anakinra",
                        "prednisolone (free)"),
               Units = c("mg/L", "mg/L", "mg/L", "mg/L", "ng/mL", "ng/mL", "ng/mL"),
               Trough = sprintf("%.4f", c(tail(s$CTMPo, 1) * 0.55, tail(s$CITZo, 1),
                                          tail(s$COHIo, 1), tail(s$CVORo, 1),
                                          tail(s$CIFNo, 1), tail(s$CANAo, 1),
                                          tail(s$CPREDo, 1) * 0.25)))
  })

  # ---- 6 -------------------------------------------------------------------
  output$focus_b <- renderPlot({
    sim() %>% mutate(total = BACT + BACTi) %>%
      select(time, BACT, BACTi, total, NEC, NEUT_T) %>%
      pivot_longer(-time) %>% mutate(value = pmax(value, 1e-2)) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_y_log10() +
      geom_hline(yintercept = 1e7, linetype = 2) +
      annotate("text", x = 0, y = 1.4e7, hjust = 0, size = 3.4,
               label = "1e7 = a clinically serious infection") +
      labs(x = "day", y = "CFU (log scale) / cells", colour = NULL) + THEME
  })
  output$focus_b_text <- renderUI(HTML(
    "<p><b>BACTi is the point.</b> An organism ingested by an oxidase-null
     neutrophil is not merely un-killed: it is hidden from every competent
     neutrophil in the lesion, it goes on dividing, and it is released alive
     when its host cell dies. Each uptake returns s_enc * k_release /
     (k_release - mu_bi) organisms to the extracellular pool — 0.06 at phi = 1
     and 0.75 at phi = 0. NEC is what makes the outcome depend on the SIZE of
     the inoculum: a focus that outruns recruitment walls itself off from the
     neutrophils and from the antibiotic alike, which is why a CGD abscess has
     to be drained.</p>"))

  # ---- 7 -------------------------------------------------------------------
  output$focus_f <- renderPlot({
    sim() %>% select(time, CONID, HYPH, GM, MAC, NEUT_T) %>%
      pivot_longer(-time) %>% mutate(value = pmax(value, 1e-3)) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_y_log10() + labs(x = "day", y = "burden (log) / index", colour = NULL) +
      THEME
  })
  output$focus_f_text <- renderUI(HTML(
    "<p><b>The asymmetry that makes CGD a mould disease.</b> A conidium is
     ingestible, so it sees exactly the kernel's s(phi) inside a macrophage
     phagosome. A germinated hypha is far too large to swallow and can only be
     attacked from the outside — by oxidants and NETs released onto its
     surface, the least salvageable form of killing there is. Aspergillus is
     the most oxidant-dependent organism in the model's table
     (omega = 0.92) and the leading cause of death in the registries.</p>"))

  # ---- 8 -------------------------------------------------------------------
  output$ncrit <- renderPlot({
    g <- GENOTYPES
    g$dhr <- dhr_mean(g$f_normal, g$phi_res)
    g$arrangement <- ifelse(g$f_normal > 0 & g$phi_res == 0, "mosaic",
                     ifelse(g$f_normal == 0, "uniform", "normal"))
    ggplot(g, aes(dhr, -log10(s_enc(f_normal, phi_res)), colour = arrangement,
                  label = name)) +
      geom_point(size = 3) + geom_text(hjust = -0.08, size = 3.2) +
      xlim(0, 1.35) +
      labs(x = "DHR mean", y = "population log-kill per encounter",
           colour = NULL,
           title = "Genotypes at matched DHR do not have matched killing") + THEME
  })
  output$ncrit_tbl <- renderTable({
    g <- GENOTYPES
    data.frame(genotype = g$name,
               DHR = sprintf("%.2f", dhr_mean(g$f_normal, g$phi_res)),
               s_enc = sprintf("%.4f", s_enc(g$f_normal, g$phi_res)),
               `log10 N_crit (from the Python reference)` =
                 c("5.69", "3.68", "3.69", "4.06", "3.71", "3.68",
                   "4.69", "4.14", "3.81"), check.names = FALSE)
  })

  # ---- 9 -------------------------------------------------------------------
  rate_tbl <- reactive({
    # N_crit is computed by the Python reference (it needs the focus ODE solved
    # by bisection); the app interpolates from s_enc, which is monotone in it.
    ncrit <- function(f, r, drug = 0) {
      base <- CAL$Ncrit_cgd +
        (CAL$Ncrit_hlt - CAL$Ncrit_cgd) *
        (-log10(s_enc(f, r, boost())) + log10(s_enc(0, 0, 1))) /
        (-log10(s_enc(1, 1, 1)) + log10(s_enc(0, 0, 1)))
      base + drug
    }
    dB <- if (input$tmpsmx) 0.69 else 0
    dF <- if (input$itra || input$vori) 0.98 else 0
    rb <- infection_rate(ncrit(fx(), input$phi_res, dB), CAL$lambda_b, CAL$muN_b)
    rf <- infection_rate(ncrit(fx(), input$phi_res, dF) - 0.78,
                         CAL$lambda_f, CAL$muN_f)
    data.frame(Endpoint = c("serious bacterial infections / patient-year",
                            "invasive fungal infections / patient-year",
                            "predicted deaths / patient-year"),
               Value = sprintf("%.3f", c(rb, rf, rb * CAL$CFR_B + rf * CAL$CFR_F)))
  })
  output$rates_tbl <- renderTable(rate_tbl())
  output$trials <- renderPlot({
    df <- data.frame(
      trial = rep(c("co-trimoxazole\n(Margolis 1990)",
                    "itraconazole\n(Gallin 2003)",
                    "interferon gamma-1b\n(ICGDCSG 1991)"), each = 2),
      source = rep(c("model (predicted)", "trial (observed)"), 3),
      reduction = c(0.61, 0.56, 0.53, 0.86, 0.73, 0.67))
    ggplot(df, aes(trial, reduction, fill = source)) +
      geom_col(position = "dodge") +
      scale_y_continuous(labels = scales::percent) +
      labs(x = NULL, y = "reduction in infection", fill = NULL,
           title = "Three trials the model predicts rather than fits") + THEME
  })
  output$trials_text <- renderUI(HTML(
    "<p>Only FOUR numbers were fitted at the whole-body level: two bacterial
     and two fungal exposure parameters, calibrated on the untreated X-CGD
     rate and the healthy-control rate. Every drug effect above is a
     prediction from the simulated steady-state PK. The itraconazole arm is
     the poor one (53% against 86%), and the direction is informative: the
     model gives azoles an effect on hyphal extension only, and not on
     germination.</p>"))

  # ---- 10 ------------------------------------------------------------------
  output$inflam <- renderPlot({
    sim() %>% select(time, IL1B, IL18, IL6, IL17, IL10, TNF, CRP, APOP,
                     GRAN, COL, FIB, KYN) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#B03060", linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") + labs(x = "day", y = NULL) + THEME
  })
  output$inflam_text <- renderUI(HTML(
    "<p><b>Set both inocula to zero and watch the granuloma appear anyway.</b>
     The CGD neutrophil dies without externalising the oxidised
     phosphatidylserine that marks it for clearance, so the corpse persists —
     and the corpse is simultaneously an IL-1beta stimulus and the missing
     source of the IL-10 that would have switched the lesion off. That is why
     anakinra abolishes the colitis and co-trimoxazole does not: there was
     never an organism in it. Note that IL-10 does not FALL in CGD; it rises,
     because there are more corpses. What fails is IL-10 per corpse.</p>"))

  # ---- 11 ------------------------------------------------------------------
  output$threshold <- renderPlot({
    f <- seq(0, 1, by = 0.01)
    lk <- -log10(sapply(f, function(x) s_enc(x, input$phi_res, boost())))
    lk0 <- lk[1]; lk1 <- tail(lk, 1)
    ncrit <- CAL$Ncrit_cgd + (CAL$Ncrit_hlt - CAL$Ncrit_cgd) * (lk - lk0) / (lk1 - lk0)
    rb <- sapply(ncrit, infection_rate, lambda = CAL$lambda_b, muN = CAL$muN_b)
    rf <- sapply(ncrit - 0.78, infection_rate, lambda = CAL$lambda_f, muN = CAL$muN_f)
    df <- data.frame(corrected = rep(f, 3),
                     value = c(rb, rf, rb * CAL$CFR_B + rf * CAL$CFR_F),
                     what = rep(c("bacterial infections/yr", "fungal infections/yr",
                                  "deaths/yr"), each = length(f)))
    ggplot(df, aes(corrected, value)) + geom_line(linewidth = 1, colour = "#6A4C93") +
      facet_wrap(~what, scales = "free_y") +
      geom_vline(xintercept = c(0.10, 0.20), linetype = 2, colour = "grey50") +
      scale_x_continuous(labels = scales::percent) +
      labs(x = "fraction of myelopoiesis with a competent oxidase",
           y = NULL,
           title = "Published thresholds (>20% resolves, <10% recurs) shown dashed") +
      THEME
  })
  output$threshold_text <- renderUI(HTML(
    "<p>The Python reference solves this properly, by bisecting the focus ODE
     for N_crit at every value of f; it finds that ~12% corrected myelopoiesis
     halves the untreated mortality and ~44% removes 90% of it. Published
     experience puts the actionable threshold above 20% donor myeloid
     chimerism, with recurrence below 10% (Guengoer 2014, Marciano 2018), and
     lentiviral gene therapy reached 16-46% DHR-positive neutrophils at 12
     months with clinical resolution in 6 of 9 (Kohn 2020). The model was
     given none of those numbers.</p>"))

  # ---- 12 ------------------------------------------------------------------
  SCEN <- list(
    list("1 healthy control", 1, 1, list(), 0, 0, 0, 0, 365),
    list("2 X-CGD null, no prophylaxis", 0, 0, list(), 0, 0, 0, 0, 365),
    list("3 + co-trimoxazole", 0, 0, list(tmpsmx = TRUE), 0, 0, 0, 0, 365),
    list("4 + co-trimoxazole + itraconazole", 0, 0,
         list(tmpsmx = TRUE, itra = TRUE), 0, 0, 0, 0, 365),
    list("5 + triple prophylaxis", 0, 0,
         list(tmpsmx = TRUE, itra = TRUE, ifng = TRUE), 0, 0, 0, 0, 365),
    list("6 p47phox AR + triple", 0, 0.12,
         list(tmpsmx = TRUE, itra = TRUE, ifng = TRUE), 0, 0, 0, 0, 365),
    list("7 X-CGD hypomorph 20%", 0, 0.20, list(), 0, 0, 0, 0, 365),
    list("8 X-CGD carrier 20% (same DHR)", 0.20, 0, list(), 0, 0, 0, 0, 365),
    list("9 staph inoculum 1e4, no rx", 0, 0, list(), 0, 0, 1e4, 0, 90),
    list("10 staph 1e4 + co-trimoxazole", 0, 0, list(tmpsmx = TRUE), 0, 0, 1e4, 0, 90),
    list("11 established abscess 1e6 + same", 0, 0, list(tmpsmx = TRUE), 0, 0, 1e6, 0, 90),
    list("12 invasive aspergillosis, no rx", 0, 0, list(), 0, 0, 0, 1e4, 90),
    list("13 same + voriconazole + IFN-g", 0, 0,
         list(vori = TRUE, ifng = TRUE), 0, 0, 0, 1e4, 90),
    list("14 colitis + prednisolone", 0, 0, list(pred = TRUE), 0, 0, 0, 0, 365),
    list("15 colitis + anakinra", 0, 0, list(anakinra = TRUE), 0, 0, 0, 0, 365),
    list("16 HSCT 95% / gene therapy 20%", 0, 0, list(), 0.95, 0, 0, 0, 365)
  )
  scen_run <- reactive({
    do.call(rbind, lapply(SCEN, function(s) {
      rg <- do.call(cgd_regimen, c(list(days = s[[9]]), s[[4]]))
      out <- run_patient(s[[2]], s[[3]], s[[9]], rg, chim = s[[5]],
                         vcn = s[[6]], seedB = s[[7]], seedC = s[[8]])
      last <- tail(out, 1)
      data.frame(scenario = s[[1]], ANC = last$ANC, CRP = last$CRP,
                 IL1B = last$IL1B, granuloma = last$GRAN, colitis = last$COL,
                 logB = last$logB, logH = last$logH, ALT = last$ALT,
                 survival = last$SURV)
    }))
  })
  output$scen_tbl <- renderTable(scen_run(), digits = 3)
  output$scen_plot <- renderPlot({
    scen_run() %>% select(scenario, granuloma, colitis, logB, logH) %>%
      pivot_longer(-scenario) %>%
      ggplot(aes(reorder(scenario, value), value, fill = name)) +
      geom_col() + coord_flip() + facet_wrap(~name, scales = "free_x", nrow = 1) +
      labs(x = NULL, y = NULL, fill = NULL) + THEME
  })
}

shinyApp(ui, server)
