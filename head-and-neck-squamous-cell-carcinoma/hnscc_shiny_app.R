## =============================================================================
##  hnscc_shiny_app.R
##  Head and Neck Squamous Cell Carcinoma — interactive QSP dashboard
##
##  Eleven tabs.  The design principle is that the app should make the model's
##  ONE claim checkable by hand rather than merely displayed:
##
##      Delta(log kill) = (agent effect on factor i) x (HEADROOM in factor i)
##
##  So the "Resistance headroom" tab plots the four measured headrooms side by
##  side with the log kill each agent actually bought, and the "Matched
##  contrasts" tab runs pairs of arms that differ in exactly one thing.  Every
##  other tab is a conventional read-out.
##
##  Run:   shiny::runApp("hnscc_shiny_app.R")
##  Needs: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT (DT optional)
##
##  Licence: see repository LICENSE.  EDUCATIONAL / RESEARCH USE ONLY.
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)

source("hnscc_mrgsolve_model.R")   # builds `mod` and the dosing helpers

HAS_DT <- requireNamespace("DT", quietly = TRUE)

theme_set(theme_minimal(base_size = 12) +
          theme(panel.grid.minor = element_blank(),
                strip.background = element_rect(fill = "#eef2f7", colour = NA),
                legend.position  = "bottom"))

PAL <- c("#2b6cb0", "#c05621", "#2f855a", "#805ad5", "#b83280",
         "#4a5568", "#d69e2e", "#319795")

at <- function(x, col, t) x[[col]][which.min(abs(x$time - t))]

## -----------------------------------------------------------------------------
##  UI
## -----------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Head & Neck Squamous Cell Carcinoma — QSP model explorer"),
  tags$p(style = "color:#666;margin-top:-8px",
         "72 ODEs. Four resistance factors (hypoxia x repair x repopulation x ",
         "immune escape). Every headroom shown is measured by a counterfactual ",
         "integrator running alongside the real one, not assumed."),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Tumour & patient"),
      radioButtons("hpv", "Aetiology", inline = TRUE,
                   choices = c("HPV-negative" = 0, "HPV-positive" = 1),
                   selected = 0),
      sliderInput("v0",  "Gross tumour volume (cm3)", 5, 90, 32, step = 1),
      sliderInput("py",  "Tobacco pack-years", 0, 60, 35, step = 5),
      checkboxInput("smokrt", "Continues smoking during radiotherapy", FALSE),
      sliderInput("cps", "PD-L1 combined positive score (CPS)", 0, 100, 10, step = 1),
      sliderInput("hla", "HLA class-I / B2M loss (fraction of clones)", 0, 0.6, 0.15, step = 0.01),
      sliderInput("hbg", "Haemoglobin (g/dL)", 8, 16, 13.5, step = 0.1),
      sliderInput("crcl","Creatinine clearance (mL/min)", 30, 130, 100, step = 5),
      sliderInput("age", "Age (years)", 30, 90, 60, step = 1),
      sliderInput("ecog","ECOG performance status", 0, 3, 1, step = 1),
      sliderInput("bsa", "Body surface area (m2)", 1.3, 2.4, 1.8, step = 0.05),

      hr(),
      h4("Radiotherapy"),
      sliderInput("gy",  "Dose per fraction (Gy)", 1.0, 2.5, 2.0, step = 0.1),
      sliderInput("nfx", "Number of fractions", 20, 70, 35, step = 1),
      radioButtons("perday", "Fractions per day", inline = TRUE,
                   choices = c("1" = 1, "2 (hyperfractionated)" = 2), selected = 1),
      sliderInput("gapafter", "Unplanned break after fraction (0 = none)", 0, 34, 0, step = 1),
      sliderInput("gaplen",   "Length of break (days)", 0, 21, 0, step = 1),
      sliderInput("fpar", "Fraction of dose to parotid (0.35 IMRT, 0.80 3D-CRT)",
                  0.1, 0.95, 0.35, step = 0.05),

      hr(),
      h4("Systemic therapy"),
      checkboxGroupInput("drugs", NULL,
        choices = c("Cisplatin 100 mg/m2 q3w x3"   = "cis3",
                    "Cisplatin 40 mg/m2 weekly x7" = "cisw",
                    "Carboplatin AUC 5 q3w x3"     = "carbo",
                    "Cetuximab 400 -> 250 mg/m2/wk"= "cetux",
                    "Pembrolizumab 200 mg q3w"     = "pembro",
                    "5-FU 1000 mg/m2/d x4 q3w"     = "fu",
                    "Docetaxel 75 mg/m2 q3w"       = "doce",
                    "Nimorazole (hypoxic sensitiser)" = "nimo"),
        selected = "cis3"),

      hr(),
      sliderInput("end", "Simulation horizon (days)", 200, 1825, 730, step = 5),
      actionButton("go", "Simulate", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient profile",     tableOutput("profile"),
                                            plotOutput("phenoPlot", height = 330)),
        tabPanel("2 · Drug PK",             plotOutput("pkPlot",    height = 620)),
        tabPanel("3 · Target occupancy",    plotOutput("occPlot",   height = 480),
                                            tableOutput("occTab")),
        tabPanel("4 · Tumour populations",  plotOutput("popPlot",   height = 620)),
        tabPanel("5 · Radiobiology",        plotOutput("rbPlot",    height = 620)),
        tabPanel("6 · Resistance headroom", plotOutput("hdPlot",    height = 420),
                                            tableOutput("hdTab")),
        tabPanel("7 · Immune microenvironment", plotOutput("immPlot", height = 620)),
        tabPanel("8 · Clinical endpoints",  plotOutput("epPlot",    height = 480),
                                            tableOutput("epTab")),
        tabPanel("9 · Toxicity",            plotOutput("toxPlot",   height = 700)),
        tabPanel("10 · Matched contrasts",  plotOutput("cmpPlot",   height = 430),
                                            if (HAS_DT) DT::dataTableOutput("cmpTab")
                                            else tableOutput("cmpTab")),
        tabPanel("11 · Biomarkers",         plotOutput("bmPlot",    height = 620))
      )
    )
  )
)

## -----------------------------------------------------------------------------
##  SERVER
## -----------------------------------------------------------------------------
server <- function(input, output, session) {

  covars <- reactive({
    list(HPV = as.numeric(input$hpv), PY = input$py,
         SMOKRT = as.numeric(input$smokrt), CPS = input$cps,
         HLALOSS = input$hla, HBG = input$hbg, CRCL0 = input$crcl,
         AGE = input$age, ECOG = input$ecog, BSA = input$bsa,
         V0 = input$v0, FPAR = input$fpar,
         CARBOF = as.numeric("carbo" %in% input$drugs),
         NIMOF  = as.numeric("nimo"  %in% input$drugs))
  })

  ## Build the event table from the sidebar
  build_ev <- function(drugs, bsa) {
    ga <- if (input$gapafter > 0 && input$gaplen > 0) input$gapafter else NA
    ev <- rt_course(input$gy, input$nfx, 0, as.numeric(input$perday),
                    gapafter = ga, gaplen = input$gaplen)
    if ("cis3"   %in% drugs) ev <- c(ev, cis_q3w(100, bsa, 3, 0))
    if ("cisw"   %in% drugs) ev <- c(ev, cis_wk(40, bsa, 7, 0))
    if ("carbo"  %in% drugs) ev <- c(ev, carbo_q3w(5, input$crcl, 3, 0))
    if ("cetux"  %in% drugs) ev <- c(ev, cetux(bsa, 8, -7))
    if ("pembro" %in% drugs) ev <- c(ev, pembro(200, 17, 0))
    if ("fu"     %in% drugs) ev <- c(ev, fu_ci(1000, bsa, 4, 3, 0))
    if ("doce"   %in% drugs) ev <- c(ev, doce_q3w(75, bsa, 3, 0))
    ev
  }

  sim <- eventReactive(input$go, {
    m  <- zero_re(param(mod, covars()))
    ev <- build_ev(input$drugs, input$bsa)
    mrgsim_df(m, events = ev, end = input$end, delta = 0.5)
  }, ignoreNULL = FALSE)

  long <- function(d, cols, labels = cols) {
    d %>% select(time, all_of(cols)) %>%
      pivot_longer(-time, names_to = "var", values_to = "value") %>%
      mutate(var = factor(var, levels = cols, labels = labels))
  }

  fx_line <- function() {
    ga <- if (input$gapafter > 0 && input$gaplen > 0) input$gapafter else NA
    tt <- fx_times(input$nfx, 0, as.numeric(input$perday), ga, input$gaplen)
    geom_vline(xintercept = range(tt), linetype = 3, colour = "#999")
  }

  ## ---------------------------------------------------------------- 1 profile
  output$profile <- renderTable({
    d <- sim()
    hpvlab <- if (input$hpv == "1") "HPV16 / p16-positive" else "HPV-negative (tobacco/alcohol)"
    data.frame(
      Quantity = c("Aetiology", "Baseline volume (cm3)", "Baseline clonogens",
                   "Baseline tumour pO2 (mmHg)", "Baseline hypoxic fraction",
                   "Baseline clonogen alpha (/Gy)", "Prescribed dose (Gy)",
                   "Overall treatment time (days)", "Systemic therapy"),
      Value = c(hpvlab,
                sprintf("%.1f", d$VOL[1]),
                format(signif(d$CSCN[1], 3), big.mark = ","),
                sprintf("%.1f", d$PO2OUT[1]),
                sprintf("%.3f", d$HFOUT[1]),
                sprintf("%.3f", d$ALPHAO[1]),
                sprintf("%.1f", input$gy * input$nfx),
                sprintf("%.0f", diff(range(fx_times(input$nfx, 0,
                          as.numeric(input$perday),
                          if (input$gapafter > 0 && input$gaplen > 0) input$gapafter else NA,
                          input$gaplen)))),
                if (length(input$drugs)) paste(input$drugs, collapse = ", ") else "none"))
  })

  output$phenoPlot <- renderPlot({
    d <- sim()
    ggplot(long(d, c("VOL", "PO2OUT", "HFOUT", "ALPHAO"),
                c("Volume (cm3)", "pO2 (mmHg)", "Hypoxic fraction",
                  "Clonogen alpha (/Gy)")),
           aes(time, value)) +
      geom_line(colour = PAL[1], linewidth = 0.8) + fx_line() +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Days", y = NULL,
           title = "Tumour geometry and the two quantities that set radiosensitivity")
  })

  ## --------------------------------------------------------------------- 2 PK
  output$pkPlot <- renderPlot({
    d <- sim() %>% filter(time <= min(input$end, 160))
    ggplot(long(d, c("CFREEO", "CCTXP", "CPEMP", "CPTOT"),
                c("Free platinum (mg/L)", "Cetuximab (mg/L)",
                  "Anti-PD-1 (mg/L)", "Cumulative platinum (mg/m2)")),
           aes(time, value)) +
      geom_line(colour = PAL[2], linewidth = 0.7) +
      facet_wrap(~var, scales = "free_y", ncol = 1) +
      labs(x = "Days", y = NULL,
           title = "Plasma exposure",
           subtitle = "Free platinum has a 25-minute half-life, so its peaks are narrow spikes on this scale")
  })

  ## ------------------------------------------------------------- 3 occupancy
  output$occPlot <- renderPlot({
    d <- sim() %>% filter(time <= min(input$end, 160))
    ggplot(long(d, c("OCCEGO", "OCCPDO", "NHEJO", "MUREPO"),
                c("EGFR occupancy", "PD-1 occupancy",
                  "NHEJ capacity (norm.)", "Sublethal repair rate (/d)")),
           aes(time, value)) +
      geom_line(colour = PAL[3], linewidth = 0.8) +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Days", y = NULL,
           title = "Target engagement and what it does to repair",
           subtitle = "Cetuximab lowers NHEJ capacity through nuclear EGFR; platinum adducts lower it further")
  })

  output$occTab <- renderTable({
    d <- sim()
    data.frame(
      Metric = c("Peak EGFR occupancy", "Peak PD-1 occupancy",
                 "Minimum NHEJ capacity", "Minimum sublethal repair rate (/d)",
                 "Peak platinum adduct load"),
      Value = sprintf("%.3f", c(max(d$OCCEGO), max(d$OCCPDO),
                                min(d$NHEJO), min(d$MUREPO), max(d$PTDNA))))
  })

  ## ---------------------------------------------------------- 4 populations
  output$popPlot <- renderPlot({
    d <- sim()
    p1 <- ggplot(long(d, c("CSCN", "CSCHF"),
                      c("Surviving clonogens", "Hypoxic fraction of clonogens")),
                 aes(time, value)) +
      geom_line(colour = PAL[1], linewidth = 0.8) + fx_line() +
      facet_wrap(~var, scales = "free_y") +
      scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1)) +
      labs(x = "Days", y = NULL,
           title = "Clonogens: the only cells that decide locoregional control",
           subtitle = "The hypoxic share falls as the tumour shrinks -- that is reoxygenation, and it is an output")
    p2 <- ggplot(long(d, c("TUMOX", "TUMHY", "TUMDM", "TUMNC"),
                      c("Oxic proliferating", "Hypoxic quiescent",
                        "Lethally damaged", "Necrotic")),
                 aes(time, value, colour = var)) +
      geom_line(linewidth = 0.8) + scale_colour_manual(values = PAL) +
      labs(x = "Days", y = "Cells", colour = NULL, title = "Bulk tumour compartments")
    gridExtra_stack(p1, p2)
  })

  ## --------------------------------------------------------- 5 radiobiology
  output$rbPlot <- renderPlot({
    d <- sim() %>% filter(time <= min(input$end, 200))
    ggplot(long(d, c("OERBAR", "RHOO", "HDPOP", "TOTLK"),
                c("Population-average OER", "Self-renewal probability rho",
                  "Log-kill handed back by repopulation (fraction)",
                  "Net log kill (nats)")),
           aes(time, value)) +
      geom_line(colour = PAL[4], linewidth = 0.8) + fx_line() +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Days", y = NULL,
           title = "The five R's, as they actually behave in this course",
           subtitle = "rho and HDPOP rise about three weeks in: that is the accelerated-repopulation kick-off, set by KDAMOUT")
  })

  ## ------------------------------------------------------------- 6 headroom
  output$hdPlot <- renderPlot({
    d <- sim()
    hd <- long(d, c("HDHYP", "HDREPR", "HDPOP", "HDIMM"),
               c("Hypoxia", "Intact repair", "Repopulation", "Immune escape")) %>%
      filter(time > 5, time <= min(input$end, 200))
    ggplot(hd, aes(time, value, colour = var)) +
      geom_line(linewidth = 1) + scale_colour_manual(values = PAL[c(1,2,3,4)]) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(x = "Days", y = "Fraction of achievable log kill lost", colour = NULL,
           title = "The four resistance headrooms, measured",
           subtitle = "Each is 1 - (actual log kill) / (log kill with that one factor switched off)")
  })

  output$hdTab <- renderTable({
    d <- sim(); t50 <- min(50, max(d$time))
    data.frame(
      Factor = c("Hypoxia", "Intact repair", "Repopulation", "Immune escape",
                 "Composite resistance R"),
      `Headroom at day 50` = c(sprintf("%.1f %%", 100 * at(d, "HDHYP", t50)),
                               sprintf("%.1f %%", 100 * at(d, "HDREPR", t50)),
                               sprintf("%.1f %%", 100 * at(d, "HDPOP", t50)),
                               sprintf("%.1f %%", 100 * at(d, "HDIMM", t50)),
                               sprintf("%.2f x",  at(d, "RESIST", t50))),
      Interpretation = c("removable by nimorazole / carbogen / transfusion",
                         "removable by platinum (alpha up, NHEJ down)",
                         "removable by shortening treatment time, partly by cetuximab",
                         "removable by PD-1 blockade",
                         "product of the four"),
      check.names = FALSE)
  })

  ## --------------------------------------------------------------- 7 immune
  output$immPlot <- renderPlot({
    d <- sim()
    ggplot(long(d, c("CD8T", "TEX", "TREG", "MDSC", "IFNG", "PDL1"),
                c("Effector CD8", "Exhausted CD8", "Treg", "MDSC",
                  "IFN-gamma", "PD-L1")),
           aes(time, value)) +
      geom_line(colour = PAL[5], linewidth = 0.8) + fx_line() +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Days", y = "Normalised", title = "Tumour immune microenvironment",
           subtitle = "PD-L1 sits at CPS/20 plus IFN-gamma-driven adaptive induction")
  })

  ## ------------------------------------------------------------ 8 endpoints
  output$epPlot <- renderPlot({
    d <- sim()
    ggplot(long(d, c("RECST", "LRC", "OS3", "VOL"),
                c("Change in longest diameter (%)",
                  "Locoregional control probability",
                  "3-year overall survival surrogate", "Volume (cm3)")),
           aes(time, value)) +
      geom_line(colour = PAL[1], linewidth = 0.9) + fx_line() +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Days", y = NULL, title = "Clinical endpoints")
  })

  output$epTab <- renderTable({
    d <- sim(); tt <- min(730, max(d$time))
    data.frame(
      Endpoint = c("Peak net log10 kill", "Clonogen nadir",
                   "Locoregional control", "3-year OS surrogate",
                   "Best RECIST change (%)", "Volume at end of horizon (cm3)"),
      Value = c(sprintf("%.2f", max(d$LOGK10)),
                signif(min(d$CSCN), 3),
                sprintf("%.3f", at(d, "LRC", tt)),
                sprintf("%.3f", at(d, "OS3", tt)),
                sprintf("%+.1f", min(d$RECST)),
                sprintf("%.1f", d$VOL[nrow(d)])))
  })

  ## ------------------------------------------------------------- 9 toxicity
  output$toxPlot <- renderPlot({
    d <- sim()
    ggplot(long(d, c("MUCGR", "SALFLO", "XEROG", "DYSPHG", "HLOSS",
                     "EGFRR", "ANCNAD", "RASHG", "PCTWL"),
                c("Mucositis grade", "Salivary flow (mL/min)", "Xerostomia grade",
                  "Dysphagia grade", "Hearing shift (dB)", "eGFR (mL/min)",
                  "ANC (10^9/L)", "Rash grade", "Weight loss (%)")),
           aes(time, value)) +
      geom_line(colour = PAL[6], linewidth = 0.8) + fx_line() +
      facet_wrap(~var, scales = "free_y", ncol = 3) +
      labs(x = "Days", y = NULL,
           title = "Normal-tissue and systemic toxicity",
           subtitle = "Mucositis is the competing clock: it is what forces the treatment breaks that feed repopulation")
  })

  ## ----------------------------------------------------- 10 matched contrasts
  ##  Each pair differs in exactly one thing, so the difference is attributable.
  contrasts <- reactive({
    bsa <- input$bsa
    RTx <- rt_course(2, 35); CISx <- cis_q3w(100, bsa, 3, 0); CTXx <- cetux(bsa, 8, -7)
    arms <- list(
      list("HPV- RT alone",        HPVNEG, RTx),
      list("HPV- + cisplatin",     HPVNEG, c(RTx, CISx)),
      list("HPV- + cetuximab",     HPVNEG, c(RTx, CTXx)),
      list("HPV- + both",          HPVNEG, c(RTx, CISx, CTXx)),
      list("HPV- + nimorazole",    modifyList(HPVNEG, list(NIMOF = 1)), RTx),
      list("HPV+ RT alone",        HPVPOS, RTx),
      list("HPV+ + cisplatin",     HPVPOS, c(RTx, CISx)),
      list("HPV+ + cetuximab",     HPVPOS, c(RTx, CTXx)),
      list("HPV+ + nimorazole",    modifyList(HPVPOS, list(NIMOF = 1)), RTx)
    )
    do.call(rbind, lapply(arms, function(a) {
      o <- mrgsim_df(zero_re(param(mod, a[[2]])), events = a[[3]], end = 400, delta = 2)
      data.frame(arm = a[[1]], logkill = max(o$LOGK10), LRC = max(o$LRC),
                 HD_hyp = at(o, "HDHYP", 50), HD_rep = at(o, "HDREPR", 50),
                 HD_pop = at(o, "HDPOP", 50), HD_imm = at(o, "HDIMM", 50),
                 mucositis = max(o$MUCGR), stringsAsFactors = FALSE)
    }))
  })

  output$cmpPlot <- renderPlot({
    cc <- contrasts()
    cc$grp <- ifelse(grepl("^HPV\\+", cc$arm), "HPV-positive", "HPV-negative")
    ggplot(cc, aes(reorder(arm, logkill), logkill, fill = grp)) +
      geom_col(width = 0.65) +
      geom_text(aes(label = sprintf("%.2f", logkill)), hjust = -0.15, size = 3.4) +
      coord_flip() + scale_fill_manual(values = PAL[c(1, 3)]) +
      expand_limits(y = max(cc$logkill) * 1.12) +
      labs(x = NULL, y = "Net log10 kill delivered to the clonogen pool",
           fill = NULL,
           title = "Matched arms: identical 70 Gy, one thing changed at a time",
           subtitle = paste("Cetuximab buys less in HPV-positive disease because the EGFR and",
                            "repopulation headroom it targets is smaller there"))
  })

  cmp_table <- reactive({
    cc <- contrasts()
    cc[] <- lapply(cc, function(x) if (is.numeric(x)) round(x, 3) else x)
    cc
  })
  if (HAS_DT) {
    output$cmpTab <- DT::renderDataTable(cmp_table(), options = list(dom = "t", pageLength = 12))
  } else {
    output$cmpTab <- renderTable(cmp_table())
  }

  ## ---------------------------------------------------------- 11 biomarkers
  output$bmPlot <- renderPlot({
    d <- sim()
    ggplot(long(d, c("CTD", "IL6", "PDL1", "MGSER", "TSHX", "PTDNA"),
                c("HPV ctDNA (copies/mL)", "IL-6 (norm.)", "PD-L1 (norm.)",
                  "Serum magnesium (mmol/L)", "TSH (mIU/L)",
                  "Platinum-DNA adducts (norm.)")),
           aes(time, value)) +
      geom_line(colour = PAL[8], linewidth = 0.8) + fx_line() +
      facet_wrap(~var, scales = "free_y") +
      labs(x = "Days", y = NULL, title = "Circulating and tissue biomarkers",
           subtitle = "ctDNA is released only by HPV-positive tumours and clears with a 2.6-hour half-life")
  })
}

## A tiny vertical stacker so the app does not need gridExtra
gridExtra_stack <- function(p1, p2) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(p1, p2, ncol = 1)
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  } else {
    p1
  }
}

shinyApp(ui, server)
