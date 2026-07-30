## ============================================================================
##  Chagas disease / Chronic Chagas Cardiomyopathy — Shiny dashboard
##  14 tabs over the 68-ODE mrgsolve model in chg_mrgsolve_model.R
## ----------------------------------------------------------------------------
##  The app is organised around the four things the model is FOR, not around
##  the compartments it happens to contain:
##
##    (1) tabs 3-5  the reservoir: why blood qPCR and cure are different
##                  endpoints, and why a more potent drug can lose
##    (2) tabs 6-9  the myocardium: how a parasite that is barely there
##                  produces a cardiomyopathy that certainly is
##    (3) tabs 10-12 the trials: BENDITA, STOP-CHAGAS and BENEFIT reproduced,
##                  and the parasite-attributable fraction that explains all three
##    (4) tab 13    the therapeutic index: 2 weeks versus 8 weeks
##
##  Run with:  shiny::runApp("chg_shiny_app.R")
## ============================================================================

suppressPackageStartupMessages({
  library(shiny); library(mrgsolve); library(ggplot2); library(dplyr); library(tidyr)
})

source("chg_mrgsolve_model.R")
MOD  <- build_chg()
YEAR <- 365.25

theme_chg <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey35", size = 10),
          legend.position = "bottom")
}
PAL <- c("#1F5C99", "#B03030", "#2E8B57", "#CC7700", "#7A4FA3", "#008B8B", "#999999")

## ---------------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Chagas disease / chronic Chagas cardiomyopathy — QSP dashboard"),
  tags$p(style = "color:#555;margin-top:-8px",
         "68-ODE mrgsolve model. Three time-scales: killing (days), the dormant reservoir (months), ",
         "myocardial remodelling (decades). Every trial anchor is in chg_references.md."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("susc", "Susceptibility (SUSC)", 0.05, 1.5, 0.55, 0.05),
      sliderInput("age0", "Age at infection (years)", 1, 55, 30, 1),
      sliderInput("horizon", "Simulation horizon (years)", 2, 60, 40, 1),
      selectInput("strain", "Parasite strain (DTU)",
                  c("TcII / TcV / TcVI  (Southern Cone)" = 1,
                    "TcI  (Colombia, Central America)"   = 8), selected = 1),
      sliderInput("gitrop", "Digestive tropism", 0.2, 4, 1, 0.2),

      h4("Trypanocidal therapy"),
      selectInput("drug", "Agent",
                  c("none", "benznidazole", "nifurtimox", "posaconazole",
                    "fosravuconazole", "fexinidazole",
                    "benznidazole + posaconazole",
                    "hypothetical sterilising agent")),
      sliderInput("tstart", "Given after N years of infection", 0, 40, 25, 1),
      sliderInput("dose", "Daily dose (mg)", 50, 900, 300, 50),
      sliderInput("ddays", "Duration (days)", 3, 365, 60, 1),
      checkboxInput("weekly", "Weekly instead of daily dosing", FALSE),
      checkboxInput("sterile", "Follow the STERILISED branch of the cure lottery", TRUE),

      h4("Cardiovascular therapy"),
      checkboxInput("acei", "ACE inhibitor / ARB", FALSE),
      checkboxInput("mra",  "MRA (spironolactone)", FALSE),
      checkboxInput("carv", "Carvedilol", FALSE),
      checkboxInput("amio", "Amiodarone", FALSE),
      checkboxInput("icd",  "ICD", FALSE),
      checkboxInput("oac",  "Anticoagulation", FALSE),
      checkboxInput("afib", "Atrial fibrillation present", FALSE),
      hr(),
      helpText("Comparator in every plot is the same patient, never treated.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("1 · Patient profile",  plotOutput("p_profile", height = 560), tableOutput("t_profile")),
        tabPanel("2 · Drug PK",          plotOutput("p_pk", height = 560), verbatimTextOutput("v_pk")),
        tabPanel("3 · Parasite dynamics",plotOutput("p_para", height = 560), htmlOutput("h_para")),
        tabPanel("4 · Reservoir & cure", plotOutput("p_cure", height = 520), htmlOutput("h_cure")),
        tabPanel("5 · qPCR vs cure",     plotOutput("p_pcr", height = 520), htmlOutput("h_pcr")),
        tabPanel("6 · Immunity",         plotOutput("p_imm", height = 560)),
        tabPanel("7 · Myocardial remodelling", plotOutput("p_myo", height = 560)),
        tabPanel("8 · Arrhythmia substrate",   plotOutput("p_arr", height = 560), htmlOutput("h_arr")),
        tabPanel("9 · Clinical endpoints",     plotOutput("p_end", height = 560), tableOutput("t_end")),
        tabPanel("10 · Trial reproduction",    plotOutput("p_trial", height = 560), tableOutput("t_trial")),
        tabPanel("11 · Attributable fraction", plotOutput("p_paf", height = 520), htmlOutput("h_paf")),
        tabPanel("12 · Age at treatment",      plotOutput("p_age", height = 560), tableOutput("t_age")),
        tabPanel("13 · Safety & index",        plotOutput("p_safe", height = 520), htmlOutput("h_safe")),
        tabPanel("14 · Digestive form",        plotOutput("p_gi", height = 520))
      )
    )
  )
)

## ---------------------------------------------------------------------------
##  server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  base_par <- reactive({
    list(SUSC = input$susc, AGE0 = input$age0, STRAIN = as.numeric(input$strain),
         GITROP = input$gitrop,
         ACEI_ON = as.numeric(input$acei), MRA_ON = as.numeric(input$mra),
         ICD_ON = as.numeric(input$icd), OAC_ON = as.numeric(input$oac),
         AFIB = as.numeric(input$afib))
  })

  ev_treat <- reactive({
    t0 <- input$tstart * YEAR
    ii <- if (input$weekly) 7 else 1
    switch(input$drug,
      "none"                        = NULL,
      "benznidazole"                = ev_bzn(t0, input$dose, input$ddays, ii),
      "nifurtimox"                  = ev_nfx(t0, input$dose/3, input$ddays),
      "posaconazole"                = ev_azl(t0, input$dose/2, input$ddays),
      "fosravuconazole"             = ev_azl(t0, input$dose, input$ddays, ii = 1),
      "fexinidazole"                = ev_fex(t0, input$dose*4, min(input$ddays, 10)),
      "benznidazole + posaconazole" = c(ev_bzn(t0, input$dose, input$ddays, ii),
                                        ev_azl(t0, 400, input$ddays)),
      "hypothetical sterilising agent" = ev_bzn(t0, input$dose, input$ddays, ii))
  })

  extra_par <- reactive({
    p <- list(ONDRUG = as.numeric(input$drug != "none"),
              STERILE = as.numeric(input$sterile))
    if (input$drug == "fosravuconazole") p$AZLTYPE <- 1
    if (input$drug == "hypothetical sterilising agent") p$EMAXCD_N <- 1.35
    p
  })

  ## treated and matched untreated trajectories -------------------------------
  sim_pair <- reactive({
    end <- input$horizon * YEAR
    dlt <- max(3, end/1500)
    cvev <- NULL
    if (input$carv) cvev <- c(cvev, ev_carv(input$tstart*YEAR, 25, end))
    if (input$amio) cvev <- c(cvev, ev_amio(input$tstart*YEAR, 200, end))
    evt <- ev_treat()
    evb <- if (is.null(evt)) cvev else if (is.null(cvev)) evt else c(evt, cvev)
    trt <- chg_run(MOD, evb, end = end, delta = dlt,
                   param = c(base_par(), extra_par()))
    unt <- chg_run(MOD, NULL, end = end, delta = dlt, param = base_par())
    trt$arm <- "treated"; unt$arm <- "untreated"
    list(trt = trt, unt = unt, both = rbind(trt, unt))
  })

  ## fine-grained run for anything that needs the parasite nadir --------------
  sim_fine <- reactive({
    end <- min(input$horizon, max(3, input$tstart + 3)) * YEAR
    evt <- ev_treat()
    chg_run(MOD, evt, end = end, delta = 1, param = c(base_par(), extra_par()))
  })

  cure_prob <- reactive({
    d <- sim_fine()
    if (input$drug == "none") return(0)
    p_cure(MOD, d, from = input$tstart*YEAR)
  })

  lineplot <- function(d, vars, labs, title, sub, ylab, logy = FALSE) {
    dd <- d %>% select(AGE, arm, all_of(vars)) %>%
      pivot_longer(all_of(vars), names_to = "v", values_to = "y")
    dd$v <- factor(dd$v, levels = vars, labels = labs)
    g <- ggplot(dd, aes(AGE, y, colour = v, linetype = arm)) +
      geom_line(linewidth = 0.75) +
      scale_colour_manual(values = PAL, name = NULL) +
      scale_linetype_manual(values = c(treated = "solid", untreated = "22"), name = NULL) +
      labs(title = title, subtitle = sub, x = "age (years)", y = ylab) + theme_chg()
    if (logy) g <- g + scale_y_log10()
    g
  }

  ## 1 patient profile --------------------------------------------------------
  output$p_profile <- renderPlot({
    d <- sim_pair()$both
    lineplot(d, c("EFout","LGE","EDV","BNP","RASSI","AGE"),
             c("EF (%)","LGE (g)","LVEDV (mL)","NT-proBNP (rel)","Rassi score","age"),
             "Patient profile over the whole natural history",
             "dashed = the same patient never treated", "value") +
      facet_wrap(~v, scales = "free_y")
  })
  output$t_profile <- renderTable({
    p <- sim_pair(); i <- nrow(p$trt)
    data.frame(
      quantity  = c("EF (%)","LGE mass (g)","LVEDV (mL)","collagen fraction",
                    "scar heterogeneity","sympathetic denervation",
                    "conduction reserve","apical aneurysm index",
                    "cumulative mortality","Rassi score"),
      treated   = round(c(p$trt$EFout[i], p$trt$LGE[i], p$trt$EDV[i], p$trt$COL[i],
                          p$trt$SCARH[i], p$trt$SYMPD[i], p$trt$COND[i], p$trt$APEX[i],
                          p$trt$MORT[i], p$trt$RASSI[i]), 3),
      untreated = round(c(p$unt$EFout[i], p$unt$LGE[i], p$unt$EDV[i], p$unt$COL[i],
                          p$unt$SCARH[i], p$unt$SYMPD[i], p$unt$COND[i], p$unt$APEX[i],
                          p$unt$MORT[i], p$unt$RASSI[i]), 3))
  })

  ## 2 PK ---------------------------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim_fine()
    d$t <- (d$time - input$tstart*YEAR)
    d <- d[d$t >= -2 & d$t <= input$ddays + 60, ]
    dd <- data.frame(t = rep(d$t, 4),
                     conc = c(d$CEN_B/32, d$CEN_A/300, d$CEN_N/100, d$CEN_M/60),
                     drug = rep(c("benznidazole","azole","nifurtimox","amiodarone"),
                                each = nrow(d)))
    dd <- dd[dd$conc > 1e-6, ]
    ggplot(dd, aes(t, conc, colour = drug)) + geom_line(linewidth = 0.8) +
      scale_y_log10() + scale_colour_manual(values = PAL, name = NULL) +
      labs(title = "Plasma concentration-time profiles",
           subtitle = "day 0 = first dose; note the azole's much longer terminal phase",
           x = "days from first dose", y = "concentration (mg/L, log scale)") + theme_chg()
  })
  output$v_pk <- renderPrint({
    d <- sim_fine(); w <- d[d$time >= input$tstart*YEAR &
                            d$time <= input$tstart*YEAR + input$ddays, ]
    cat(sprintf("mean benznidazole concentration on treatment : %.2f mg/L\n", mean(w$CEN_B)/32))
    cat(sprintf("mean azole concentration on treatment        : %.3f mg/L\n", mean(w$CEN_A)/300))
    cat(sprintf("cumulative benznidazole AUC                  : %.0f mg*day/L\n", max(w$BZNCUM)))
    cat(sprintf("mean static (replication block)              : %.3f\n", mean(w$STAT)))
    cat(sprintf("mean cidal rate, replicating forms           : %.3f /day\n", mean(w$CIDR)))
    cat(sprintf("mean cidal rate, DORMANT forms               : %.4f /day\n", mean(w$CIDD)))
  })

  ## 3 parasite ---------------------------------------------------------------
  output$p_para <- renderPlot({
    d <- sim_fine(); d$t <- d$time - input$tstart*YEAR
    dd <- data.frame(t = rep(d$t, 3),
                     n = pmax(1e-6, c(d$PTOTo - d$PDORMo, d$PDORMo, d$PBLD*5000)),
                     pool = rep(c("replicating amastigotes","DORMANT amastigotes",
                                  "blood trypomastigotes (x blood volume)"), each = nrow(d)))
    ggplot(dd, aes(t, n, colour = pool)) + geom_line(linewidth = 0.85) +
      scale_y_log10() + scale_colour_manual(values = PAL[c(1,2,3)], name = NULL) +
      geom_vline(xintercept = c(0, input$ddays), linetype = "33", colour = "grey40") +
      labs(title = "Parasite compartments through and after treatment",
           subtitle = "vertical lines bracket the treatment course",
           x = "days from first dose", y = "organism equivalents (log scale)") + theme_chg()
  })
  output$h_para <- renderUI(HTML(paste0(
    "<p style='color:#444'>The replicating pool falls fast under any agent. The dormant pool is the ",
    "one that decides the outcome, and only a <b>cidal</b> agent touches it. ",
    "Minimum total burden reached: <b>", signif(min(sim_fine()$PTOTo), 3),
    "</b> organism equivalents.</p>")))

  ## 4 reservoir & cure -------------------------------------------------------
  output$p_cure <- renderPlot({
    durs <- c(3, 7, 14, 21, 28, 42, 56, 90, 180, 365)
    res <- lapply(durs, function(dy) {
      ev <- switch(input$drug,
                   "posaconazole"    = ev_azl(input$tstart*YEAR, input$dose/2, dy),
                   "fosravuconazole" = ev_azl(input$tstart*YEAR, input$dose, dy, ii = 1),
                   "nifurtimox"      = ev_nfx(input$tstart*YEAR, input$dose/3, dy),
                   ev_bzn(input$tstart*YEAR, input$dose, dy,
                          if (input$weekly) 7 else 1))
      d <- chg_run(MOD, ev, end = (input$tstart + 2.5)*YEAR, delta = 1,
                   param = c(base_par(), extra_par()))
      pc <- p_cure(MOD, d, from = input$tstart*YEAR)
      data.frame(days = dy, cure = pc, reported = p_reported(pc))
    })
    r <- do.call(rbind, res)
    ggplot(r, aes(days, reported)) +
      geom_line(colour = PAL[1], linewidth = 0.9) +
      geom_point(colour = PAL[1], size = 2) +
      geom_hline(yintercept = 0.89, linetype = "33", colour = PAL[2]) +
      annotate("text", x = 200, y = 0.905, label = "BENDITA 8-week arm (89%)",
               colour = PAL[2], size = 3.4) +
      scale_x_log10(breaks = durs) + ylim(0, 1) +
      labs(title = "Cure saturates with duration — the BENDITA result as a curve",
           subtitle = "reported sustained parasitological clearance vs treatment duration",
           x = "treatment duration (days, log scale)",
           y = "sustained clearance a trial would report") + theme_chg()
  })
  output$h_cure <- renderUI(HTML(sprintf(
    "<p style='color:#444'>Probability of sterile cure for the selected regimen: <b>%.3f</b>; ",
    cure_prob()) %>% paste0(sprintf(
    "a trial would report <b>%.3f</b> after allowing for assay sensitivity and the per-visit ",
    p_reported(cure_prob())),
    "false-negative rate that also produces the placebo 'conversion' rate. ",
    "For a purely static agent this curve is flat and near zero at every duration.</p>")))

  ## 5 qPCR vs cure -----------------------------------------------------------
  output$p_pcr <- renderPlot({
    t0 <- input$tstart*YEAR
    arms <- list(benznidazole = ev_bzn(t0, 400, 60), posaconazole = ev_azl(t0, 400, 60))
    dd <- do.call(rbind, lapply(names(arms), function(a) {
      d <- chg_run(MOD, arms[[a]], end = t0 + 1.5*YEAR, delta = 1,
                   param = c(base_par(), list(ONDRUG = 1, STERILE = 0)))
      data.frame(t = d$time - t0, PBLD = pmax(1e-9, d$PBLD), arm = a)
    }))
    ggplot(dd, aes(t, PBLD, colour = arm)) + geom_line(linewidth = 0.85) +
      geom_hline(yintercept = param(MOD)$LOD_PCR, linetype = "33", colour = "grey30") +
      annotate("text", x = 350, y = param(MOD)$LOD_PCR*1.6, label = "qPCR limit of detection",
               size = 3.3, colour = "grey30") +
      scale_y_log10() + scale_colour_manual(values = PAL[c(1,2)], name = NULL) +
      labs(title = "The endpoint the trial measured is not the endpoint that mattered",
           subtitle = "STOP-CHAGAS reproduced: both arms qPCR-negative at day 30; only one stays negative",
           x = "days from first dose", y = "blood trypomastigotes (par-eq/mL, log)") + theme_chg()
  })
  output$h_pcr <- renderUI(HTML(
    "<p style='color:#444'>Blood trypomastigotes are a <b>flux</b>, not a reservoir: they are the ",
    "output of amastigote replication. A drug that merely arrests replication empties the blood ",
    "without killing anything, so qPCR reports success while the reservoir is intact. ",
    "Observed in STOP-CHAGAS: posaconazole 93.3% negative at day 30 (better than benznidazole's 89.7%) ",
    "and 13.3% sustained at day 180 (against benznidazole 86.7%).</p>"))

  ## 6 immunity ---------------------------------------------------------------
  output$p_imm <- renderPlot({
    d <- sim_pair()$both
    lineplot(d, c("TH1","TREG","MPHI","TNFA","TGFB","AAB"),
             c("TH1","Treg","macrophage","TNF-alpha","TGF-beta","autoantibody"),
             "Immune and autoimmune state",
             "the autoantibody curve is the one that does not come back down after cure",
             "relative units") + facet_wrap(~v, scales = "free_y")
  })

  ## 7 myocardium -------------------------------------------------------------
  output$p_myo <- renderPlot({
    d <- sim_pair()$both
    lineplot(d, c("CMYO","COL","MFB","MVD","EDV","WSo"),
             c("viable myocytes","collagen fraction","myofibroblasts",
               "microvascular dysfunction","LVEDV (mL)","wall stress (mmHg)"),
             "Myocardial remodelling",
             "collagen and wall stress form the loop that does not need a parasite",
             "value") + facet_wrap(~v, scales = "free_y")
  })

  ## 8 arrhythmia -------------------------------------------------------------
  output$p_arr <- renderPlot({
    d <- sim_pair()$both
    lineplot(d, c("SCARH","SYMPD","PSYMD","COND","SANF","APEX"),
             c("scar heterogeneity","sympathetic denervation","vagal denervation",
               "conduction reserve","sinus node reserve","apical aneurysm"),
             "Arrhythmic substrate and conduction system",
             "none of these six curves has a recovery term: all are irreversible",
             "index (0-1)") + facet_wrap(~v, scales = "free_y")
  })
  output$h_arr <- renderUI(HTML(
    "<p style='color:#444'>Sudden death in this model is driven by <b>scar heterogeneity and ",
    "sympathetic denervation</b>, not by how much fibrosis there is. That is why arrhythmia ",
    "severity in Chagas tracks the extent of the MIBG defect rather than the extent of the ",
    "perfusion defect, and it is the reason amiodarone and an ICD move a different endpoint ",
    "from the one an ACE inhibitor moves.</p>"))

  ## 9 endpoints --------------------------------------------------------------
  output$p_end <- renderPlot({
    p <- sim_pair()
    dd <- rbind(
      data.frame(AGE = p$trt$AGE, cif = p$trt$CIFSCD, mode = "sudden death",   arm = "treated"),
      data.frame(AGE = p$trt$AGE, cif = p$trt$CIFHF,  mode = "pump failure",   arm = "treated"),
      data.frame(AGE = p$trt$AGE, cif = p$trt$CIFSTK, mode = "stroke",         arm = "treated"),
      data.frame(AGE = p$unt$AGE, cif = p$unt$CIFSCD, mode = "sudden death",   arm = "untreated"),
      data.frame(AGE = p$unt$AGE, cif = p$unt$CIFHF,  mode = "pump failure",   arm = "untreated"),
      data.frame(AGE = p$unt$AGE, cif = p$unt$CIFSTK, mode = "stroke",         arm = "untreated"))
    ggplot(dd, aes(AGE, cif, colour = mode, linetype = arm)) + geom_line(linewidth = 0.85) +
      scale_colour_manual(values = PAL, name = NULL) +
      scale_linetype_manual(values = c(treated = "solid", untreated = "22"), name = NULL) +
      labs(title = "Cumulative incidence by mode of death (true competing risks)",
           subtitle = "these are CIFs, not 1 - exp(-H): suppressing one mode raises the others",
           x = "age (years)", y = "cumulative incidence") + theme_chg()
  })
  output$t_end <- renderTable({
    p <- sim_pair(); i <- nrow(p$trt)
    data.frame(endpoint  = c("sudden cardiac death","pump-failure death","fatal stroke",
                             "all-cause mortality","pacemaker / ICD","new heart failure",
                             "sustained VT","BENEFIT composite"),
               treated   = round(c(p$trt$CIFSCD[i], p$trt$CIFHF[i], p$trt$CIFSTK[i], p$trt$MORT[i],
                                   p$trt$CIFPPM[i], p$trt$CIFNHF[i], p$trt$CIFVT[i], p$trt$CIFCMP[i]), 3),
               untreated = round(c(p$unt$CIFSCD[i], p$unt$CIFHF[i], p$unt$CIFSTK[i], p$unt$MORT[i],
                                   p$unt$CIFPPM[i], p$unt$CIFNHF[i], p$unt$CIFVT[i], p$unt$CIFCMP[i]), 3))
  })

  ## 10 trial reproduction ----------------------------------------------------
  trialtab <- reactive({
    bend <- function(dose, days, ii = 1) {
      d <- chg_run(MOD, ev_bzn(0, dose, days, ii), end = 2*YEAR, delta = 1,
                   param = list(SUSC = 0.55, ONDRUG = 1))
      p_reported(p_cure(MOD, d))
    }
    dB <- chg_run(MOD, ev_bzn(0, 400, 60), end = 2*YEAR, delta = 1,
                  param = list(SUSC = 0.55, ONDRUG = 1))
    dA <- chg_run(MOD, ev_azl(0, 400, 60), end = 2*YEAR, delta = 1,
                  param = list(SUSC = 0.55, ONDRUG = 1))
    data.frame(
      trial = c(rep("BENDITA", 5), rep("STOP-CHAGAS", 2)),
      arm = c("placebo","BZN 300 x 8 wk","BZN 300 x 4 wk","BZN 300 x 2 wk","BZN 150 x 4 wk",
              "BZN 400 x 60 d","posaconazole x 60 d"),
      observed = c(0.03, 0.89, 0.89, 0.83, 0.83, 0.867, 0.133),
      predicted = round(c(p_reported(0), bend(300,56), bend(300,28), bend(300,14), bend(150,28),
                          p_reported(p_cure(MOD, dB), q_false = 0.10),
                          p_reported(p_cure(MOD, dA), q_false = 0.10)), 3))
  })
  output$p_trial <- renderPlot({
    r <- trialtab()
    ggplot(r, aes(observed, predicted, colour = trial, label = arm)) +
      geom_abline(slope = 1, linetype = "33", colour = "grey50") +
      geom_point(size = 3.2) + geom_text(vjust = -1, size = 3.1, show.legend = FALSE) +
      xlim(0, 1.05) + ylim(0, 1.05) +
      scale_colour_manual(values = PAL[c(1,2)], name = NULL) +
      labs(title = "Observed versus predicted, two independent randomised trials",
           subtitle = "points on the dashed line are exact; the model was not refitted per trial",
           x = "observed sustained parasitological clearance",
           y = "model prediction") + theme_chg()
  })
  output$t_trial <- renderTable(trialtab())

  ## 11 attributable fraction -------------------------------------------------
  output$p_paf <- renderPlot({
    d <- chg_run(MOD, NULL, end = 45*YEAR, delta = 30.4, param = base_par())
    ggplot(d, aes(AGE, PAF)) + geom_line(colour = PAL[1], linewidth = 1) +
      geom_vline(xintercept = input$age0 + input$tstart, linetype = "33", colour = PAL[2]) +
      annotate("text", x = input$age0 + input$tstart, y = 0.9,
               label = "your treatment time", hjust = -0.05, colour = PAL[2], size = 3.4) +
      ylim(0, 1) +
      labs(title = "Parasite-attributable fraction of the myocardial injury rate",
           subtitle = "there is NO age term anywhere in the model: this decay is the growth of the denominator",
           x = "age (years)", y = "PAF") + theme_chg()
  })
  output$h_paf <- renderUI({
    d <- chg_run(MOD, NULL, end = 45*YEAR, delta = 30.4, param = base_par())
    tt <- input$tstart*YEAR
    paf <- d$PAF[which.min(abs(d$time - tt))]
    pw1 <- chg_power(1 - 1.00*paf, 0.29); pw2 <- chg_power(1 - 0.20*paf, 0.29)
    HTML(sprintf(paste0("<p style='color:#444'>At your chosen treatment time the PAF is <b>%.3f</b>. ",
      "If a trypanocide sterilised <b>every</b> patient the largest achievable hazard ratio would be ",
      "<b>%.3f</b>, needing about <b>%s</b> patients at a 29%% event rate. At the sterilisation ",
      "fraction BENEFIT actually achieved (~20%%) the achievable hazard ratio is <b>%.3f</b> and the ",
      "trial would need about <b>%s</b>. BENEFIT randomised 2,854.</p>"),
      paf, 1-paf, format(round(pw1$n), big.mark = ","),
      1-0.2*paf, format(round(pw2$n), big.mark = ",")))
  })

  ## 12 age at treatment ------------------------------------------------------
  agetab <- reactive({
    times <- c(2, 6, 10, 14, 18, 22, 25, 30)
    do.call(rbind, lapply(times, function(ty) {
      d <- chg_run(MOD, ev_bzn(ty*YEAR, 300, 56), end = 50*YEAR, delta = 91.3,
                   param = list(SUSC = input$susc, AGE0 = input$age0,
                                ONDRUG = 1, STERILE = 1))
      i <- nrow(d)
      data.frame(years_of_infection = ty, EF_at_end = round(d$EFout[i], 1),
                 LGE = round(d$LGE[i], 1), mortality_50y = round(d$MORT[i], 3))
    }))
  })
  output$p_age <- renderPlot({
    r <- agetab()
    u <- chg_run(MOD, NULL, end = 50*YEAR, delta = 91.3,
                 param = list(SUSC = input$susc, AGE0 = input$age0))
    ggplot(r, aes(years_of_infection, mortality_50y)) +
      geom_line(colour = PAL[1], linewidth = 0.9) + geom_point(size = 2.4, colour = PAL[1]) +
      geom_hline(yintercept = u$MORT[nrow(u)], linetype = "33", colour = PAL[2]) +
      annotate("text", x = 22, y = u$MORT[nrow(u)] - 0.03, label = "never treated",
               colour = PAL[2], size = 3.4) +
      labs(title = "The same drug, the same dose, the same disease — only the timing differs",
           subtitle = "50-year cumulative mortality in the sterilised branch, by years of infection at treatment",
           x = "years of infection when benznidazole was given",
           y = "50-year cumulative mortality") + theme_chg()
  })
  output$t_age <- renderTable(agetab())

  ## 13 safety ----------------------------------------------------------------
  output$p_safe <- renderPlot({
    durs <- c(7, 14, 21, 28, 42, 56, 80)
    r <- do.call(rbind, lapply(durs, function(dy) {
      d <- chg_run(MOD, ev_bzn(0, 300, dy), end = 1.5*YEAR, delta = 1,
                   param = list(SUSC = 0.55, ONDRUG = 1))
      i <- nrow(d)
      disc <- 1 - exp(-(d$HRASH[i] + d$HNEU[i]))
      pc <- p_cure(MOD, d)
      data.frame(days = dy, value = c(p_reported(pc), disc),
                 what = c("sustained clearance", "permanent discontinuation"))
    }))
    ggplot(r, aes(days, value, colour = what)) + geom_line(linewidth = 0.9) +
      geom_point(size = 2.2) + scale_colour_manual(values = PAL[c(1,2)], name = NULL) +
      ylim(0, 1) +
      labs(title = "Therapeutic index of duration",
           subtitle = "efficacy saturates by two weeks; toxicity does not saturate at all",
           x = "treatment duration (days)", y = "probability") + theme_chg()
  })
  output$h_safe <- renderUI(HTML(
    "<p style='color:#444'>BENDITA reported <b>zero</b> treatment discontinuations in the 2-week arm ",
    "and 7% overall, while STOP-CHAGAS' 60-day benznidazole arms reached 31.7%. In this model the ",
    "cutaneous hazard peaks in the third week and the neuropathic hazard is a function of ",
    "<i>cumulative</i> exposure, so the two curves above separate for a structural reason rather ",
    "than a fitted one.</p>"))

  ## 14 digestive form --------------------------------------------------------
  output$p_gi <- renderPlot({
    d <- sim_pair()$both
    lineplot(d, c("ENSN","ESOD","COLD"),
             c("myenteric neurons","oesophageal dilatation","colonic dilatation"),
             "Digestive form",
             "megaviscera require >50% enteric neuron loss, and neurons do not regenerate",
             "index") + facet_wrap(~v, scales = "free_y")
  })
}

shinyApp(ui, server)
