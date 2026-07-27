## =====================================================================
##  Charcot-Marie-Tooth disease type 1A (CMT1A) — QSP Shiny dashboard
## ---------------------------------------------------------------------
##  Run with:   shiny::runApp("cmt1a_shiny_app.R")
##  Requires:   shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
##
##  The app drives the model in cmt1a_mrgsolve_model.R. Every simulation
##  starts from a healthy newborn nerve and integrates copy number over
##  the life course, so moving the copy-number control does not switch
##  between stored phenotypes - it re-generates one.
##
##  Ten tabs:
##    1  Patient           — genotype, stature, modifiers, life course
##    2  Gene dosage       — mRNA, folding, aggregates, the U-curve
##    3  Schwann cell      — c-Jun/EGR2, differentiation, myelin, NCV
##    4  Axon margins      — supply/demand by length class (the geometry)
##    5  Motor unit reserve— MUNE, unit size, the silent decade
##    6  Clinical          — CMTNS-R, ONLS, walking, strength
##    7  Biomarkers        — NfL, MRI fat fraction, CMAP, SNAP, NCV
##    8  Therapeutic window— the knockdown sweep and its far wall
##    9  Drug PK/PD        — ascorbate, PXT3003, oligonucleotide, HDAC6i
##   10  Trial designer    — endpoint choice, required N, treatment window
## =====================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

source("cmt1a_mrgsolve_model.R")

YRD <- 365.25
theme_cmt <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom")

PAL <- c("CMT1A (3 copies)" = "#c0392b", "Healthy (2 copies)" = "#2c7fb8",
         "HNPP (1 copy)"    = "#7a4fa5", "Treated" = "#218c5a",
         "foot"  = "#c0392b", "leg" = "#e08214", "hand" = "#2c7fb8",
         "proximal" = "#218c5a")

## ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Charcot-Marie-Tooth 1A — Quantitative Systems Pharmacology"),
  tags$p(style = "color:#666;margin-top:-10px;",
         "A gene-dosage model. The disease parameter is an integer; the phenotype is integrated, not assumed."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Genotype and constitution"),
      sliderInput("cn", "PMP22 copy number", min = 1, max = 4, value = 3, step = 1),
      helpText(HTML("1 = HNPP &middot; 2 = normal &middot; 3 = CMT1A &middot; 4 = homozygous duplication")),
      sliderInput("htf", "Stature factor (axon length)", 0.85, 1.20, 1.00, 0.01),
      sliderInput("modf", "Aggregate-clearance modifier", 0.6, 1.8, 1.0, 0.05),
      hr(),
      h4("Intervention"),
      selectInput("drug", "Therapy",
                  c("None" = "none",
                    "PMP22 oligonucleotide (SC q4w)" = "oli",
                    "Ascorbic acid (oral)"           = "aa",
                    "PXT3003"                        = "pxt",
                    "HDAC6 inhibitor"                = "hdi",
                    "AAV1.NT-3"                      = "nt3",
                    "Exercise programme"             = "ex",
                    "Vincristine (4 weekly doses)"   = "vcr")),
      conditionalPanel("input.drug == 'oli'",
                       sliderInput("olidose", "Dose per injection (mg)", 0, 120, 20, 5)),
      conditionalPanel("input.drug == 'aa'",
                       sliderInput("aadose", "Daily dose (mg)", 0, 6000, 4000, 250)),
      conditionalPanel("input.drug == 'pxt'",
                       sliderInput("pxtlev", "PLEO-CMT dose level", 1, 2, 2, 1)),
      conditionalPanel("input.drug == 'hdi'",
                       sliderInput("hdidose", "Daily dose (mg)", 0, 90, 30, 5)),
      sliderInput("tstart", "Treatment start (age, years)", 0, 60, 30, 1),
      sliderInput("tstop",  "Treatment stop (age, years)",  1, 80, 80, 1),
      hr(),
      checkboxInput("preg", "Simulate a pregnancy at age 30", FALSE),
      sliderInput("maxage", "Follow to age (years)", 20, 80, 70, 5),
      hr(),
      actionButton("go", "Simulate", class = "btn-primary btn-block"),
      tags$small(style = "color:#888;",
                 "A life-course run takes a few seconds. Comparators (healthy, untreated CMT1A) are always shown.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 Patient",        br(), verbatimTextOutput("summary"), br(), DTOutput("snap")),
        tabPanel("2 Gene dosage",    br(), plotOutput("p_dose", height = 620)),
        tabPanel("3 Schwann cell",   br(), plotOutput("p_sc",   height = 620)),
        tabPanel("4 Axon margins",   br(), plotOutput("p_marg", height = 620),
                 tags$p(style="color:#666;",
                        "Margin = delivered support / terminal demand. The classes differ only in length.")),
        tabPanel("5 Reserve",        br(), plotOutput("p_res",  height = 620),
                 tags$p(style="color:#666;",
                        "Motor unit number falls from the first decade; measured strength does not follow until the reserve is spent.")),
        tabPanel("6 Clinical",       br(), plotOutput("p_clin", height = 620)),
        tabPanel("7 Biomarkers",     br(), plotOutput("p_bio",  height = 620)),
        tabPanel("8 Window",         br(), plotOutput("p_win",  height = 560),
                 tags$p(style="color:#666;",
                        "The near wall is aggregate burden; the far wall is tomacula. They are different pathways, which is why the far wall is invisible on CMTNS.")),
        tabPanel("9 Drug PK/PD",     br(), plotOutput("p_pk",   height = 620)),
        tabPanel("10 Trial designer", br(), DTOutput("tbl_trial"), br(), plotOutput("p_trial", height = 380))
      )
    )
  )
)

## ---------------------------------------------------------------------
build_run <- function(input, cn = NULL, drug = NULL) {
  cn   <- if (is.null(cn)) input$cn else cn
  drug <- if (is.null(drug)) input$drug else drug
  p <- list(CN = cn, HTF = input$htf, MODF = input$modf)
  e <- NULL
  t0 <- input$tstart * YRD
  dur <- max(1, (input$tstop - input$tstart)) * YRD
  if (drug == "oli" && input$olidose > 0) {
    p$OLI_ON <- 1
    e <- ev(time = t0, amt = input$olidose, cmt = "OLISC",
            ii = 28, addl = max(0, round(dur / 28) - 1))
  } else if (drug == "aa" && input$aadose > 0) {
    p$AA_ON <- 1; p$AA_DTOT <- input$aadose
    e <- ev(time = t0, amt = (input$aadose / 3) * 5.678, cmt = "AAD",
            ii = 1/3, addl = max(0, round(dur * 3) - 1))
  } else if (drug == "pxt") {
    p$PXT_ON <- 1
    e <- pxt_events(input$pxtlev, t0, dur)
  } else if (drug == "hdi" && input$hdidose > 0) {
    p$HDI_ON <- 1
    e <- ev(time = t0, amt = input$hdidose, cmt = "HDID",
            ii = 1, addl = max(0, round(dur) - 1))
  } else if (drug == "nt3") {
    p$NT3_ON <- 1
  } else if (drug == "ex") {
    p$KEXER <- 1
  } else if (drug == "vcr") {
    p$VCR_ON <- 1
    e <- ev(time = t0, amt = 2000, cmt = "VCRC", ii = 7, addl = 3)
  }
  if (isTRUE(input$preg)) { p$PREG_ON <- 1; p$PREG_T0 <- 30 * YRD; p$SEXF <- 1 }
  m <- param(mod, p)
  out <- if (is.null(e)) mrgsim(m, end = input$maxage * YRD, delta = 20, hmax = 5)
         else mrgsim(m, events = e, end = input$maxage * YRD, delta = 20, hmax = 5)
  as.data.frame(out)
}

server <- function(input, output, session) {

  sims <- eventReactive(input$go, {
    withProgress(message = "Integrating the life course...", value = 0, {
      incProgress(0.15)
      pat <- build_run(input)
      incProgress(0.45)
      unt <- build_run(input, cn = input$cn, drug = "none")
      incProgress(0.75)
      hea <- build_run(input, cn = 2, drug = "none")
      incProgress(1)
      list(pat = pat, unt = unt, hea = hea)
    })
  }, ignoreNULL = FALSE)

  longify <- function(s, vars) {
    bind_rows(
      s$hea %>% select(AGE_Y, all_of(vars)) %>% mutate(arm = "Healthy (2 copies)"),
      s$unt %>% select(AGE_Y, all_of(vars)) %>% mutate(arm = "Untreated"),
      s$pat %>% select(AGE_Y, all_of(vars)) %>% mutate(arm = "Simulated patient")
    ) %>% pivot_longer(all_of(vars), names_to = "var", values_to = "val")
  }

  panel <- function(s, vars, labs, title) {
    d <- longify(s, vars)
    d$var <- factor(d$var, levels = vars, labels = labs)
    ggplot(d, aes(AGE_Y, val, colour = arm)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~var, scales = "free_y") +
      scale_colour_manual(values = c("Healthy (2 copies)" = "#2c7fb8",
                                     "Untreated" = "#c0392b",
                                     "Simulated patient" = "#218c5a")) +
      labs(title = title, x = "Age (years)", y = NULL, colour = NULL) + theme_cmt
  }

  ## ---- 1 patient ----
  output$summary <- renderPrint({
    s <- sims(); a <- min(input$maxage, 40)
    g <- function(d, v) d[[v]][which.min(abs(d$AGE_Y - a))]
    cat(sprintf("PMP22 copy number %d   stature factor %.2f   clearance modifier %.2f\n",
                input$cn, input$htf, input$modf))
    cat(sprintf("Therapy: %s (age %d to %d)\n\n", input$drug, input$tstart, input$tstop))
    cat(sprintf("At age %d:\n", a))
    cat(sprintf("  PMP22 mRNA      %6.3f      membrane stoichiometry gm  %6.3f\n",
                g(s$pat,"MRNA_REL"), g(s$pat,"GM_OUT")))
    cat(sprintf("  aggregate load  %6.3f      tomacula burden            %6.3f\n",
                g(s$pat,"PAGG_OUT"), g(s$pat,"TOM_OUT")))
    cat(sprintf("  motor NCV       %6.1f m/s  MUNE (leg)                 %6.1f %%\n",
                g(s$pat,"NCV_OUT"), g(s$pat,"MUNE_PCT")))
    cat(sprintf("  CMTNS-R         %6.2f      ONLS                       %6.2f\n",
                g(s$pat,"CMTNS"), g(s$pat,"ONLS")))
    cat(sprintf("  dorsiflexion    %6.1f %%    calf fat fraction          %6.1f %%\n",
                g(s$pat,"DORSI"), g(s$pat,"FFL")))
    cat(sprintf("  plasma NfL      %6.1f      pressure-palsy hazard      %6.2f /y\n",
                g(s$pat,"NFL"), g(s$pat,"PPHAZ_Y")))
    cat(sprintf("\n  untreated comparator CMTNS-R %.2f   healthy comparator %.2f\n",
                g(s$unt,"CMTNS"), g(s$hea,"CMTNS")))
  })

  output$snap <- renderDT({
    s <- sims()
    ages <- c(5, 10, 20, 30, 40, 50, 60, 70)
    ages <- ages[ages <= input$maxage]
    v <- c("NCV_OUT","MUNE_PCT","CMAP_ULN","SNAP_SUR","CMTNS","ONLS","DORSI","GRIP","FFL","NFL","T10MWT")
    tab <- do.call(rbind, lapply(ages, function(a) {
      r <- s$pat[which.min(abs(s$pat$AGE_Y - a)), v]
      data.frame(Age = a, round(r, 2), row.names = NULL)
    }))
    datatable(tab, options = list(dom = "t", pageLength = 20),
              caption = "Simulated patient trajectory")
  })

  ## ---- 2 gene dosage ----
  output$p_dose <- renderPlot({
    panel(sims(), c("MRNA_REL","GM_OUT","PAGG_OUT","TOM_OUT","KD_TOT","SCD_OUT"),
          c("PMP22 mRNA (x normal)","membrane stoichiometry gm",
            "aggregate burden","tomacula burden",
            "fractional PMP22 suppression","Schwann-cell differentiation"),
          "Gene dosage and its two downstream arms")
  })

  ## ---- 3 Schwann cell ----
  output$p_sc <- renderPlot({
    panel(sims(), c("SCD_OUT","DBG_LAC","DBG_NDI","NCV_OUT","CMAP_ULN","DBG_MITO"),
          c("myelinating differentiation","glial lactate support",
            "nodal disorganisation","motor NCV (m/s)",
            "ulnar CMAP (mV)","mitochondrial health"),
          "Schwann cell state, myelin and conduction (note the flat NCV after childhood)")
  })

  ## ---- 4 margins ----
  output$p_marg <- renderPlot({
    s <- sims()
    d <- s$pat %>%
      select(AGE_Y, foot = MARGIN_F, leg = MARGIN_L, hand = MARGIN_H, proximal = MARGIN_P) %>%
      pivot_longer(-AGE_Y, names_to = "class", values_to = "margin")
    d2 <- s$pat %>%
      select(AGE_Y, foot = AXF, leg = AXL, hand = AXH, proximal = AXP) %>%
      pivot_longer(-AGE_Y, names_to = "class", values_to = "axons")
    g1 <- ggplot(d, aes(AGE_Y, margin, colour = class)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
      scale_colour_manual(values = PAL) +
      labs(title = "Supply / demand margin by axon length class",
           subtitle = "dashed line = the threshold; classes differ only in length",
           x = "Age (years)", y = "margin", colour = NULL) + theme_cmt
    g2 <- ggplot(d2, aes(AGE_Y, 100 * axons, colour = class)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PAL) +
      labs(title = "Surviving motor axons", x = "Age (years)", y = "% of birth", colour = NULL) +
      theme_cmt
    gridExtra::grid.arrange(g1, g2, ncol = 1)
  })

  ## ---- 5 reserve ----
  output$p_res <- renderPlot({
    s <- sims()
    d <- s$pat %>% transmute(AGE_Y,
      `motor unit number (%)` = MUNE_PCT,
      `mean unit size (x normal)` = RSL * 20,
      `dorsiflexion strength (%)` = DORSI,
      `calf fat fraction (%)` = FFL) %>%
      pivot_longer(-AGE_Y)
    ggplot(d, aes(AGE_Y, value, colour = name)) + geom_line(linewidth = 1) +
      labs(title = "The silent decade: reserve is consumed long before strength moves",
           subtitle = "unit size is plotted x20 to share the axis",
           x = "Age (years)", y = NULL, colour = NULL) + theme_cmt
  })

  ## ---- 6 clinical ----
  output$p_clin <- renderPlot({
    panel(sims(), c("CMTNS","ONLS","DORSI","GRIP","T10MWT","D6MWT"),
          c("CMTNS-R (0-36)","ONLS (0-12)","ankle dorsiflexion (% pred)",
            "grip (% pred)","10-metre walk (s)","6-minute walk (m)"),
          "Clinical endpoints")
  })

  ## ---- 7 biomarkers ----
  output$p_bio <- renderPlot({
    panel(sims(), c("NFL","FFL","CMAP_ULN","SNAP_SUR","MUNE_PCT","NCV_OUT"),
          c("plasma NfL (pg/mL)","calf fat fraction (%)","ulnar CMAP (mV)",
            "sural SNAP (uV)","MUNE (%)","motor NCV (m/s)"),
          "Biomarkers, ranked left-to-right by how much they actually move")
  })

  ## ---- 8 therapeutic window ----
  win <- reactive({
    res <- lapply(seq(0, 0.80, by = 0.05), function(kd) {
      cnv <- if (kd >= 0.85) 1e6 else 3.0 * kd / (0.85 - kd)
      m <- param(mod, CN = 3, OLI_ON = 1, KNV_OUT = 0, KPL_NV = 0,
                 HTF = input$htf, MODF = input$modf)
      m <- init(m, OLINV = cnv * 0.60)
      o <- as.data.frame(mrgsim(m, end = 50 * YRD, delta = 200, hmax = 5))
      r <- o[nrow(o), ]
      data.frame(kd = 100 * kd, mRNA = r$MRNA_REL, gm = r$GM_OUT,
                 PAGG = r$PAGG_OUT, TOM = r$TOM_OUT,
                 CMTNS = r$CMTNS, PPHAZ = r$PPHAZ_Y)
    })
    bind_rows(res)
  })

  output$p_win <- renderPlot({
    w <- win()
    d <- w %>% select(kd, `CMTNS-R at 50 y` = CMTNS,
                      `pressure palsies / year` = PPHAZ,
                      `aggregate burden` = PAGG,
                      `PMP22 mRNA (x normal)` = mRNA) %>%
      pivot_longer(-kd)
    ggplot(d, aes(kd, value)) + geom_line(linewidth = 1, colour = "#c0392b") +
      geom_vline(xintercept = 33.3, linetype = 2, colour = "#218c5a") +
      geom_vline(xintercept = 66.7, linetype = 2, colour = "#7a4fa5") +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "The therapeutic window has a far wall, and it is not on the disability scale",
           subtitle = "green = 33% knockdown restores 2-copy mRNA; purple = 67% reaches the 1-copy (HNPP) level",
           x = "PMP22 knockdown (%)", y = NULL) + theme_cmt
  })

  ## ---- 9 PK ----
  output$p_pk <- renderPlot({
    s <- sims()
    v <- switch(input$drug,
      aa  = c("AA_PLASMA","KD_TOT"),
      oli = c("OLI_NRV","KD_TOT"),
      pxt = c("BAC_CONC","KD_TOT"),
      hdi = c("VEL_TR","MARGIN_F"),
      vcr = c("VEL_TR","MARGIN_F"),
      c("KD_TOT","MRNA_REL"))
    d <- s$pat %>% select(AGE_Y, all_of(v)) %>% pivot_longer(-AGE_Y)
    ggplot(d, aes(AGE_Y, value)) + geom_line(linewidth = 0.9, colour = "#2c7fb8") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(title = paste("Pharmacokinetics / pharmacodynamics:", input$drug),
           x = "Age (years)", y = NULL) + theme_cmt
  })

  ## ---- 10 trial designer ----
  output$tbl_trial <- renderDT({
    s <- sims()
    a0 <- max(input$tstart, 1); a1 <- a0 + 2
    gv <- function(d, a, v) d[[v]][which.min(abs(d$AGE_Y - a))]
    eps <- list(
      list("CMTNS-R (points)", "CMTNS", 1.8),
      list("ONLS (points)", "ONLS", 0.9),
      list("calf fat fraction (%)", "FFL", 1.6),
      list("plasma NfL (pg/mL)", "NFL", 6.0),
      list("ulnar CMAP (mV)", "CMAP_ULN", 0.9),
      list("motor NCV (m/s)", "NCV_OUT", 2.5))
    rows <- lapply(eps, function(e) {
      du <- gv(s$unt, a1, e[[2]]) - gv(s$unt, a0, e[[2]])
      dt <- gv(s$pat, a1, e[[2]]) - gv(s$pat, a0, e[[2]])
      n <- if (abs(dt - du) < 1e-9) Inf else ceiling(2 * (1.96 + 0.84)^2 * e[[3]]^2 / (dt - du)^2)
      data.frame(Endpoint = e[[1]], `Assumed SD` = e[[3]],
                 `24-mo placebo` = round(du, 4), `24-mo treated` = round(dt, 4),
                 `N per arm` = ifelse(is.finite(n) && n < 1e7, format(n, big.mark = ","), ">10,000,000"),
                 check.names = FALSE)
    })
    datatable(bind_rows(rows), options = list(dom = "t"),
              caption = "Sample size for 80% power, alpha 0.05, two-year trial starting at the treatment-start age")
  })

  output$p_trial <- renderPlot({
    s <- sims()
    d <- bind_rows(
      s$unt %>% transmute(AGE_Y, CMTNS, arm = "Untreated"),
      s$pat %>% transmute(AGE_Y, CMTNS, arm = "Simulated patient"))
    ggplot(d, aes(AGE_Y, CMTNS, colour = arm)) + geom_line(linewidth = 1) +
      annotate("rect", xmin = input$tstart, xmax = min(input$tstop, input$maxage),
               ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "#218c5a") +
      scale_colour_manual(values = c("Untreated" = "#c0392b", "Simulated patient" = "#218c5a")) +
      labs(title = "Treatment window", x = "Age (years)", y = "CMTNS-R", colour = NULL) +
      theme_cmt
  })
}

shinyApp(ui, server)
