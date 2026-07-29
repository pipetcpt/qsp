# =============================================================================
#  Dravet Syndrome QSP model — Shiny dashboard
# =============================================================================
#
#  Run with:
#      Rscript -e 'shiny::runApp("dravet_shiny_app.R", port = 8080)'
#
#  Requires: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
#  The ODE model itself lives in dravet_mrgsolve_model.R and is sourced below,
#  so the app and the command-line model can never drift apart.
#
#  Twelve tabs, each answering one question:
#     1  Patient & regimen ...... who is this, and what are they taking
#     2  Drug exposure .......... PK of every moiety, including the metabolites
#     3  Drug interactions ...... the CYP layer and the norclobazam multiplier
#     4  Target biology ......... SCN1A -> Nav1.1 -> interneuron firing capacity
#     5  E:I balance ............ where inhibition is lost and what restores it
#     6  Seizure endpoints ...... MCSF, the trial-style read-out, status
#     7  Route decomposition .... how much of the effect is the drug-drug
#                                 interaction, and how much is the drug
#     8  Sodium-channel paradox . the same drug in two hosts
#     9  Fever challenge ........ thermal susceptibility
#    10  Therapeutic index ...... seizure reduction bought per unit of sedation
#    11  Safety ................. somnolence, weight, transaminases, valve
#    12  Long-horizon outcomes .. development, status, SUDEP over years
# =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

source("dravet_mrgsolve_model.R", local = TRUE)

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  strip.background = element_rect(fill = "#eef2f7", colour = NA),
                  plot.title = element_text(face = "bold")))
PAL <- c("#2c6fbb", "#c0392b", "#2e7d32", "#b8860b", "#6a3d9a", "#00838f",
         "#8d6e63")

GENO_CHOICES <- c("Ultrarapid (UM)" = 1.60, "Normal (NM)" = 1.00,
                  "Intermediate (IM)" = 0.55, "Poor (PM)" = 0.20)

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  titlePanel("Dravet Syndrome — Quantitative Systems Pharmacology Explorer"),
  tags$p(style = "color:#555;margin-top:-8px;",
         "SCN1A haploinsufficiency → Nav1.1 loss in GABAergic interneurons",
         "→ threshold failure of fast-spiking inhibition → seizures.",
         tags$b("Educational model — not for clinical use.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Host"),
      sliderInput("allele", "Functional SCN1A allele dose", 0.30, 1.00, 0.50,
                  step = 0.05),
      helpText("0.5 = Dravet haploinsufficiency. 1.0 = healthy control.",
               "This one number is the disease; every drug effect below is",
               "computed from the same equations regardless of its value."),
      selectInput("geno", "CYP2C19 genotype", GENO_CHOICES, selected = 1.00),
      sliderInput("age0", "Age at start (years)", 0.5, 12, 1, step = 0.5),

      h4("Background therapy"),
      sliderInput("clb", "Clobazam (mg/kg/day)", 0, 1.5, 0.5, step = 0.05),
      sliderInput("vpa", "Valproate (mg/kg/day)", 0, 60, 30, step = 5),

      h4("Add-on (from day 42)"),
      sliderInput("stp", "Stiripentol (mg/kg/day)", 0, 100, 0, step = 5),
      sliderInput("cbd", "Cannabidiol (mg/kg/day)", 0, 25, 0, step = 1),
      sliderInput("ffa", "Fenfluramine (mg/kg/day)", 0, 0.7, 0, step = 0.05),
      sliderInput("aso", "Zorevunersen (mg intrathecal q90d)", 0, 90, 0,
                  step = 10),
      sliderInput("ltg", "Lamotrigine (mg/kg/day)", 0, 10, 0, step = 0.5),
      tags$div(style = "color:#b03030;font-size:11px;",
               "Lamotrigine is contraindicated in Dravet syndrome. It is",
               "included precisely so the model can be asked to reproduce",
               "that fact rather than be told it."),
      checkboxInput("keto", "Ketogenic diet", FALSE),

      h4("Challenge"),
      checkboxInput("fever", "Febrile illness on day 100", FALSE),
      sliderInput("tend", "Simulate to day", 140, 1825, 140, step = 20),
      hr(),
      h4("Decomposition switches"),
      checkboxInput("pk_on", "Allow CYP inhibition (PK route)", TRUE),
      checkboxInput("pd_on", "Allow direct target engagement (PD route)", TRUE),
      helpText("Turn one off to see how much of the benefit travels down",
               "each path into the same clinical endpoint.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",
        tabPanel("1 Patient & regimen",
                 br(), DTOutput("tbl_regimen"),
                 br(), verbatimTextOutput("txt_summary")),
        tabPanel("2 Drug exposure", br(), plotOutput("p_pk", height = 620)),
        tabPanel("3 Drug interactions",
                 br(), plotOutput("p_ddi", height = 340),
                 br(), DTOutput("tbl_ddi")),
        tabPanel("4 Target biology", br(), plotOutput("p_target", height = 560),
                 br(), htmlOutput("txt_target")),
        tabPanel("5 E:I balance", br(), plotOutput("p_ei", height = 560)),
        tabPanel("6 Seizure endpoints",
                 br(), plotOutput("p_sz", height = 400),
                 br(), DTOutput("tbl_endpoints")),
        tabPanel("7 Route decomposition",
                 br(), h4("How much of the effect is the drug-drug interaction?"),
                 DTOutput("tbl_decomp"), br(),
                 plotOutput("p_headroom", height = 320), br(),
                 htmlOutput("txt_decomp")),
        tabPanel("8 Sodium-channel paradox",
                 br(), plotOutput("p_flip", height = 420),
                 br(), DTOutput("tbl_flip"), br(), htmlOutput("txt_flip")),
        tabPanel("9 Fever challenge", br(), plotOutput("p_fever", height = 560)),
        tabPanel("10 Therapeutic index",
                 br(), plotOutput("p_ti", height = 420),
                 br(), DTOutput("tbl_ti")),
        tabPanel("11 Safety", br(), plotOutput("p_safety", height = 620)),
        tabPanel("12 Long-horizon outcomes",
                 br(), plotOutput("p_long", height = 560),
                 br(), DTOutput("tbl_long"))
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  pset <- reactive({
    list(ALLELE = input$allele, E2C19_GENO = as.numeric(input$geno),
         AGE0 = input$age0, KETO = as.numeric(input$keto),
         PK_ROUTE = as.numeric(input$pk_on), PD_ROUTE = as.numeric(input$pd_on),
         BG_CLB = input$clb, BG_VPA = input$vpa)
  })

  regimen <- reactive({
    tend <- input$tend
    ev <- bind_ev(
      background_reg(input$clb, input$vpa, tend),
      if (input$stp > 0) titrated("STP_G", input$stp, 3, tend),
      if (input$cbd > 0) titrated("CBD_G", input$cbd, 2, tend),
      if (input$ffa > 0) titrated("FFA_G", input$ffa, 2, tend),
      if (input$ltg > 0) dose_seq("NVB_G", input$ltg, ADDON_START, tend, 2),
      if (input$aso > 0)
        data.frame(ID = 1,
                   time = seq(ADDON_START, tend, by = 90),
                   amt = input$aso, cmt = "ASO_CSF", evid = 1),
      if (input$fever)
        data.frame(ID = 1, time = 100, amt = 1.0, cmt = "INFECT", evid = 1)
    )
    ev
  })

  sim <- reactive({
    m <- param(mod, pset())
    as.data.frame(mrgsim_d(m, regimen(), end = input$tend, delta = 0.05,
                           atol = 1e-9, rtol = 1e-9, hmax = 0.1))
  })

  # --- helper: run an arbitrary variant of the current patient --------------
  variant <- function(over_param = list(), over_dose = NULL, tend = NULL) {
    tend <- if (is.null(tend)) input$tend else tend
    ev <- if (is.null(over_dose)) regimen() else over_dose
    m <- param(mod, modifyList(pset(), over_param))
    as.data.frame(mrgsim_d(m, ev, end = tend, delta = 0.05,
                           atol = 1e-9, rtol = 1e-9, hmax = 0.1))
  }

  mcsf_pair <- function(d) {
    cum <- function(tt) approx(d$time, d$BURD, xout = tt, ties = "ordered")$y
    b <- (cum(BASE_END) - cum(BASE_START)) / (BASE_END - BASE_START) * 30.4375
    hi <- min(max(d$time), TRT_END)
    lo <- min(TRT_START, hi - 1)
    t2 <- (cum(hi) - cum(lo)) / (hi - lo) * 30.4375
    c(baseline = b, treated = t2, reduction = 100 * (b - t2) / b)
  }

  # ------------------------------------------------------ 1 patient & regimen
  output$tbl_regimen <- renderDT({
    r <- regimen()
    if (!nrow(r)) return(datatable(data.frame(note = "no doses")))
    r %>% group_by(cmt) %>%
      summarise(`first dose (day)` = min(time), `last dose (day)` = max(time),
                `n doses` = n(), `mg/kg per dose` = round(max(amt), 4),
                .groups = "drop") %>%
      rename(compartment = cmt) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  output$txt_summary <- renderPrint({
    d <- sim(); p <- mcsf_pair(d)
    tl <- tail(d, 1)
    cat(sprintf("Host: SCN1A allele dose %.2f (%s)\n", input$allele,
                if (input$allele >= 0.99) "healthy control" else "Dravet"))
    cat(sprintf("CYP2C19 activity %.2f\n\n", as.numeric(input$geno)))
    cat(sprintf("Baseline MCSF (days %d-%d)      : %8.2f /month\n",
                BASE_START, BASE_END, p[["baseline"]]))
    cat(sprintf("On-treatment MCSF              : %8.2f /month\n", p[["treated"]]))
    cat(sprintf("Reduction                      : %8.1f %%\n\n", p[["reduction"]]))
    cat(sprintf("Norclobazam                    : %8.3f mg/L\n", tl$C_NCLB))
    cat(sprintf("N-CLB : CLB ratio              : %8.2f\n", tl$RATIO_NCLB))
    cat(sprintf("Interneuron firing capacity    : %8.4f\n", tl$CAP_INT))
    cat(sprintf("E:I ratio                      : %8.4f\n", tl$EI))
    cat(sprintf("Somnolence (0-1)               : %8.3f\n", tl$SOMN))
    cat(sprintf("Status epilepticus / year      : %8.2f\n", tl$SE_YR))
    cat(sprintf("Developmental quotient         : %8.1f\n", tl$DQ))
  })

  # ------------------------------------------------------------ 2 drug exposure
  output$p_pk <- renderPlot({
    d <- sim()
    long <- d %>%
      select(time, Clobazam = C_CLB, Norclobazam = C_NCLB,
             Stiripentol = C_STP, Cannabidiol = C_CBD, `7-OH-CBD` = C_7OH,
             Fenfluramine = C_FFA, Norfenfluramine = C_NOR,
             Valproate = C_VPA, Lamotrigine = C_NVB) %>%
      pivot_longer(-time) %>%
      group_by(name) %>% filter(max(value) > 1e-9) %>% ungroup()
    if (!nrow(long)) return(NULL)
    ggplot(long, aes(time, value, colour = name)) +
      geom_line(linewidth = 0.6, show.legend = FALSE) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = rep(PAL, 3)) +
      labs(title = "Plasma concentrations of every modelled moiety",
           subtitle = paste("Metabolites are shown separately because in this",
                            "disease they are not a detail: norclobazam is the",
                            "main active species."),
           x = "Day", y = "mg/L")
  })

  # ------------------------------------------------------- 3 drug interactions
  output$p_ddi <- renderPlot({
    d <- sim() %>%
      select(time, CYP2C19 = E2C19, CYP3A4 = E3A4, `CYP1A2/2B6/2D6` = E1A2,
             UGT1A4 = EUGT) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = 0.8) + scale_colour_manual(values = PAL) +
      coord_cartesian(ylim = c(0, NA)) +
      labs(title = "Enzyme activity — where the drug-drug interactions live",
           subtitle = "Relative to the patient's own genotype-determined baseline",
           x = "Day", y = "Relative activity", colour = NULL)
  })

  output$tbl_ddi <- renderDT({
    ref <- variant(over_dose = background_reg(input$clb, input$vpa, input$tend))
    cur <- sim()
    tl <- function(x) mean(tail(x, 400))
    data.frame(
      moiety = c("Clobazam", "Norclobazam", "N-CLB : CLB ratio",
                 "CYP2C19 activity"),
      `background only` = round(c(tl(ref$C_CLB), tl(ref$C_NCLB),
                                  tl(ref$RATIO_NCLB), tl(ref$E2C19)), 3),
      `current regimen` = round(c(tl(cur$C_CLB), tl(cur$C_NCLB),
                                  tl(cur$RATIO_NCLB), tl(cur$E2C19)), 3),
      `fold change` = round(c(tl(cur$C_CLB) / tl(ref$C_CLB),
                              tl(cur$C_NCLB) / tl(ref$C_NCLB),
                              tl(cur$RATIO_NCLB) / tl(ref$RATIO_NCLB),
                              tl(cur$E2C19) / tl(ref$E2C19)), 2),
      check.names = FALSE) %>%
      datatable(rownames = FALSE, options = list(dom = "t"),
                caption = paste("Observed for comparison: cannabidiol raises",
                                "norclobazam ~6.0x and clobazam ~1.6x",
                                "(Geffrey 2015, PMID 26114620); stiripentol",
                                "raises norclobazam ~2.5x."))
  })

  # ---------------------------------------------------------- 4 target biology
  output$p_target <- renderPlot({
    d <- sim() %>%
      select(time, `Nav1.1 function` = NAV_INT,
             `Interneuron firing capacity` = CAP_INT,
             `Pyramidal excitability` = CAP_PYR,
             `GABA-A receptor availability` = RGABA,
             `Core temperature` = TCORE, `Brain ASO` = ASO_BR) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = 0.8, show.legend = FALSE) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = rep(PAL, 2)) +
      labs(title = "From transcript to firing capacity",
           x = "Day", y = NULL)
  })

  output$txt_target <- renderUI({
    ceil <- input$allele / (1 - 0.33)
    HTML(sprintf(
      "<p><b>The antisense ceiling is arithmetic, not pharmacological.</b>
       Skipping the poison exon recovers the transcripts the <i>intact</i>
       allele wastes and nothing more. With an allele dose of %.2f the ceiling
       is <b>%.3f</b> against a healthy value of 1.000, so this approach
       cannot reach a normal phenotype however well it works.</p>",
      input$allele, ceil))
  })

  # ------------------------------------------------------------- 5 E:I balance
  output$p_ei <- renderPlot({
    d <- sim() %>%
      select(time, `Inhibitory drive` = INH, `E:I ratio` = EI,
             `Benzodiazepine-site gain` = GAIN_PAM,
             `Stiripentol-site gain` = GAIN_STP,
             `PAM occupancy` = PAM) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = 0.8, show.legend = FALSE) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = rep(PAL, 2)) +
      labs(title = "Excitation-inhibition balance and what modulates it",
           subtitle = paste("The benzodiazepine site and the stiripentol site",
                            "are separate saturable terms, because they are",
                            "separate binding sites."),
           x = "Day", y = NULL)
  })

  # -------------------------------------------------------- 6 seizure endpoints
  output$p_sz <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_line(aes(y = MCSF, colour = "instantaneous"), linewidth = 0.5,
                alpha = 0.55) +
      geom_line(aes(y = MCSF_F, colour = "28-day filtered"), linewidth = 1.0) +
      annotate("rect", xmin = BASE_START, xmax = BASE_END, ymin = -Inf,
               ymax = Inf, alpha = 0.10, fill = PAL[1]) +
      annotate("rect", xmin = TRT_START, xmax = min(TRT_END, input$tend),
               ymin = -Inf, ymax = Inf, alpha = 0.10, fill = PAL[3]) +
      scale_colour_manual(values = c("instantaneous" = PAL[7],
                                     "28-day filtered" = PAL[1])) +
      labs(title = "Monthly convulsive seizure frequency",
           subtitle = paste("Blue band = baseline observation window;",
                            "green band = maintenance window. The trial-style",
                            "read-out compares the two."),
           x = "Day", y = "Seizures / month", colour = NULL)
  })

  output$tbl_endpoints <- renderDT({
    d <- sim(); p <- mcsf_pair(d); tl <- tail(d, 1)
    data.frame(
      endpoint = c("Baseline MCSF (/month)", "On-treatment MCSF (/month)",
                   "Reduction (%)", "Status epilepticus (/year)",
                   "Cumulative convulsive seizures",
                   "Cumulative SE episodes", "SUDEP risk to date (%)"),
      value = round(c(p[["baseline"]], p[["treated"]], p[["reduction"]],
                      tl$SE_YR, tl$BURD, tl$SECNT, tl$SUDEP_CUM), 2)) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # ---------------------------------------------------- 7 route decomposition
  decomp_tbl <- reactive({
    both <- mcsf_pair(variant())[["reduction"]]
    pk <- mcsf_pair(variant(list(PD_ROUTE = 0)))[["reduction"]]
    pdo <- mcsf_pair(variant(list(PK_ROUTE = 0)))[["reduction"]]
    data.frame(
      route = c("Both routes (as prescribed)",
                "PK route only (CYP inhibition, direct target deleted)",
                "PD route only (direct target, CYP inhibition deleted)"),
      `reduction (%)` = round(c(both, pk, pdo), 1),
      `share of total (%)` = c(100,
                               round(100 * pk / both, 1),
                               round(100 * pdo / both, 1)),
      check.names = FALSE)
  })

  output$tbl_decomp <- renderDT(
    datatable(decomp_tbl(), rownames = FALSE, options = list(dom = "t")))

  output$p_headroom <- renderPlot({
    xs <- c(1, 1.5, 2, 3, 4, 6, 8, 10)
    base <- as.list(param(mod))$CL_NCLB
    res <- vapply(xs, function(x) {
      d <- variant(list(CL_NCLB = base / x))
      c(mcsf_pair(d)[["baseline"]], mean(tail(d$C_NCLB, 400)))
    }, numeric(2))
    df <- data.frame(fold = xs, mcsf = res[1, ], nclb = res[2, ])
    ggplot(df, aes(nclb * 1000, mcsf)) +
      geom_line(colour = PAL[2], linewidth = 1) +
      geom_point(colour = PAL[2], size = 2) +
      geom_vline(xintercept = 1103, linetype = 2, colour = PAL[1]) +
      annotate("text", x = 1103, y = max(df$mcsf), hjust = -0.05, size = 3.4,
               colour = PAL[1],
               label = "1103 ng/mL: Hashi 2015 ≥" ) +
      labs(title = "Headroom test: force norclobazam up and change nothing else",
           subtitle = paste("If the benzodiazepine site is already saturated at",
                            "the concentration a standard clobazam dose",
                            "reaches, no interaction can buy much."),
           x = "Norclobazam (ng/mL)", y = "MCSF (/month)")
  })

  output$txt_decomp <- renderUI({
    t <- decomp_tbl()
    HTML(sprintf(
      "<p>With the current regimen the CYP-inhibition route accounts for
       <b>%.1f%%</b> of the total effect and direct target engagement for
       <b>%.1f%%</b>. The reason the interaction buys so little is the curve
       above: Hashi 2015 found norclobazam near 1100&nbsp;ng/mL already
       associated with &ge;90%% seizure control, and a standard clobazam dose
       in this model reaches about 1600&nbsp;ng/mL &mdash; so the site is past
       its useful range <i>before</i> any interacting drug is added. Note also
       that the curve turns back up at high concentrations, because sustained
       occupancy recruits receptor tolerance.</p>",
      t[2, 3], t[3, 3]))
  })

  # ------------------------------------------------- 8 sodium-channel paradox
  flip_data <- reactive({
    dose <- if (input$ltg > 0) input$ltg else 5
    lapply(c(Dravet = 0.5, `Healthy control` = 1.0), function(a) {
      ev <- bind_ev(background_reg(input$clb, input$vpa, 140),
                    dose_seq("NVB_G", dose, ADDON_START, 140, 2))
      variant(list(ALLELE = a), over_dose = ev, tend = 140)
    })
  })

  output$p_flip <- renderPlot({
    fd <- flip_data()
    df <- bind_rows(lapply(names(fd), function(n)
      transform(fd[[n]][, c("time", "MCSF", "CAP_INT")], host = n)))
    ggplot(df, aes(time, MCSF, colour = host)) +
      geom_line(linewidth = 0.9) +
      geom_vline(xintercept = ADDON_START, linetype = 2) +
      scale_y_log10() +
      scale_colour_manual(values = c(Dravet = PAL[2],
                                     `Healthy control` = PAL[3])) +
      labs(title = "The same sodium-channel blocker in two hosts",
           subtitle = paste("Identical drug, dose and equations. Only the",
                            "SCN1A allele dose differs. Log scale."),
           x = "Day", y = "MCSF / month (log)", colour = NULL)
  })

  output$tbl_flip <- renderDT({
    fd <- flip_data()
    rows <- lapply(names(fd), function(n) {
      p <- mcsf_pair(fd[[n]])
      data.frame(host = n, `baseline MCSF` = round(p[["baseline"]], 2),
                 `on drug` = round(p[["treated"]], 2),
                 `change (%)` = round(-p[["reduction"]], 1),
                 `interneuron capacity` = round(tail(fd[[n]]$CAP_INT, 1), 4),
                 check.names = FALSE)
    })
    datatable(bind_rows(rows), rownames = FALSE, options = list(dom = "t"))
  })

  output$txt_flip <- renderUI(HTML(
    "<p><b>This is the model's falsification test.</b> Nothing in the code
     asks whether the patient has Dravet syndrome before computing a drug
     effect. The sign flip comes from one asymmetry: sustaining 100-200&nbsp;Hz
     firing costs a fast-spiking interneuron far more available sodium
     conductance than a pyramidal cell needs for 5-20&nbsp;Hz, so interneuron
     firing capacity is a steep threshold function (Hill&nbsp;6) while
     pyramidal excitability is shallow (Hill&nbsp;1.5). A healthy brain sits
     far above that threshold and tolerates a partial block; a Dravet brain
     sits just above it and falls off. Set the allele slider to 1.0 and the
     same drug becomes an anticonvulsant.</p>"))

  # ----------------------------------------------------------- 9 fever
  output$p_fever <- renderPlot({
    ev <- bind_ev(regimen(),
                  data.frame(ID = 1, time = 100, amt = 1.0, cmt = "INFECT",
                             evid = 1))
    hosts <- c(Dravet = 0.5, `Healthy control` = 1.0)
    df <- bind_rows(lapply(names(hosts), function(n) {
      d <- variant(list(ALLELE = hosts[[n]]), over_dose = ev, tend = 140)
      transform(d[d$time > 90 & d$time < 120,
                  c("time", "TCORE", "CAP_INT", "MCSF")], host = n)
    }))
    long <- pivot_longer(df, c(TCORE, CAP_INT, MCSF))
    long$name <- factor(long$name, c("TCORE", "CAP_INT", "MCSF"),
                        c("Core temperature (°C)",
                          "Interneuron firing capacity",
                          "MCSF-equivalent rate"))
    ggplot(long, aes(time, value, colour = host)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      scale_colour_manual(values = c(Dravet = PAL[2],
                                     `Healthy control` = PAL[3])) +
      labs(title = "Febrile challenge on day 100",
           subtitle = paste("Thermal sensitivity is mutant-selective, so the",
                            "same fever costs the two hosts very different",
                            "amounts of inhibition."),
           x = "Day", y = NULL, colour = NULL)
  })

  # ------------------------------------------------------ 10 therapeutic index
  ti_tbl <- reactive({
    tend <- 140
    arms <- list(
      list(lab = "Background only", ev = background_reg(input$clb, input$vpa, tend)),
      list(lab = "+ stiripentol 50",
           ev = bind_ev(background_reg(input$clb, input$vpa, tend),
                        titrated("STP_G", 50, 3, tend))),
      list(lab = "+ cannabidiol 20",
           ev = bind_ev(background_reg(input$clb, input$vpa, tend),
                        titrated("CBD_G", 20, 2, tend))),
      list(lab = "+ fenfluramine 0.7",
           ev = bind_ev(background_reg(input$clb, input$vpa, tend),
                        titrated("FFA_G", 0.7, 2, tend))),
      list(lab = "+ zorevunersen 60",
           ev = bind_ev(background_reg(input$clb, input$vpa, tend),
                        data.frame(ID = 1, time = ADDON_START, amt = 60,
                                   cmt = "ASO_CSF", evid = 1))),
      list(lab = "Clobazam doubled",
           ev = background_reg(min(input$clb * 2, 2), input$vpa, tend))
    )
    ref <- NULL
    rows <- lapply(arms, function(a) {
      d <- variant(over_dose = a$ev, tend = tend)
      p <- mcsf_pair(d); som <- tail(d$SOMN, 1)
      if (is.null(ref)) ref <<- c(p[["treated"]], som)
      gain <- 100 * (ref[1] - p[["treated"]]) / ref[1]
      dsom <- som - ref[2]
      data.frame(arm = a$lab, MCSF = round(p[["treated"]], 2),
                 `gain vs reference (pp)` = round(gain, 1),
                 somnolence = round(som, 3),
                 `added sedation` = round(dsom, 3),
                 `pp per unit sedation` =
                   ifelse(abs(dsom) < 1e-4, NA, round(gain / dsom, 0)),
                 check.names = FALSE)
    })
    bind_rows(rows)
  })

  output$tbl_ti <- renderDT(datatable(
    ti_tbl(), rownames = FALSE, options = list(dom = "t"),
    caption = paste("Blank in the last column means the arm reduced seizures",
                    "with no added sedation at all, because it does not touch",
                    "the benzodiazepine site.")))

  output$p_ti <- renderPlot({
    df <- ti_tbl()
    ggplot(df, aes(somnolence, MCSF, label = arm)) +
      geom_point(size = 3.2, colour = PAL[1]) +
      geom_text(hjust = -0.08, size = 3.4) +
      expand_limits(x = max(df$somnolence) * 1.5) +
      labs(title = "Seizure control against sedation",
           subtitle = paste("Down and to the left is better. Arms that work",
                            "through norclobazam move right as they move",
                            "down; arms that do not, move straight down."),
           x = "Somnolence (0-1)", y = "MCSF / month")
  })

  # ------------------------------------------------------------- 11 safety
  output$p_safety <- renderPlot({
    d <- sim() %>%
      select(time, `Somnolence (0-1)` = SOMN, `Weight z-change` = WGT,
             `ALT (U/L)` = ALT, `Valve index` = VLV,
             `Cumulative SUDEP risk (%)` = SUDEP_CUM,
             `Status epilepticus / year` = SE_YR) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = 0.8, show.legend = FALSE) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = rep(PAL, 2)) +
      labs(title = "Treatment-emergent safety read-outs",
           subtitle = paste("Somnolence tracks norclobazam, not the drug that",
                            "was added; the transaminase signal needs",
                            "valproate as well as cannabidiol."),
           x = "Day", y = NULL)
  })

  # ------------------------------------------------- 12 long-horizon outcomes
  output$p_long <- renderPlot({
    tend <- max(input$tend, 1095)
    ev <- bind_ev(regimen(),
                  data.frame(ID = 1, time = seq(150, tend, by = 91),
                             amt = 0.9, cmt = "INFECT", evid = 1))
    d <- variant(over_dose = ev, tend = tend) %>%
      select(time, `Developmental quotient` = DQ,
             `Cumulative seizures` = BURD,
             `Cumulative SE episodes` = SECNT,
             `Cumulative SUDEP risk (%)` = SUDEP_CUM) %>%
      pivot_longer(-time)
    ggplot(d, aes(time / 365.25, value, colour = name)) +
      geom_line(linewidth = 0.9, show.legend = FALSE) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = rep(PAL, 2)) +
      labs(title = "Long-horizon outcomes with a febrile illness every quarter",
           subtitle = paste("Developmental loss is driven by seizure burden",
                            "AND by sedation, which is why a regimen can be",
                            "worse for development than the seizures it",
                            "prevents."),
           x = "Age from start (years)", y = NULL)
  })

  output$tbl_long <- renderDT({
    tend <- max(input$tend, 1095)
    ev <- bind_ev(regimen(),
                  data.frame(ID = 1, time = seq(150, tend, by = 91),
                             amt = 0.9, cmt = "INFECT", evid = 1))
    tl <- tail(variant(over_dose = ev, tend = tend), 1)
    data.frame(
      outcome = c("Developmental quotient", "Cumulative convulsive seizures",
                  "Cumulative status epilepticus episodes",
                  "Cumulative SUDEP risk (%)", "Weight z-change",
                  "ALT (U/L)", "Valve index"),
      value = round(c(tl$DQ, tl$BURD, tl$SECNT, tl$SUDEP_CUM, tl$WGT,
                      tl$ALT, tl$VLV), 2)) %>%
      datatable(rownames = FALSE, options = list(dom = "t"),
                caption = sprintf("At day %.0f (%.1f years from start).",
                                  tend, tend / 365.25))
  })
}

shinyApp(ui, server)
