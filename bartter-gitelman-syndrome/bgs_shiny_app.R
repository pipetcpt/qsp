## =============================================================================
##  bgs_shiny_app.R — Bartter & Gitelman Syndrome QSP Dashboard
##
##  An interactive front-end for bgs_mrgsolve_model.R. The app is organised
##  around one question: what changes when you move the SAME transport lesion
##  from the thick ascending limb to the distal convoluted tubule?
##
##  10 tabs:
##    1  Patient & genotype       — set FTAL / FDCT and see the lesion position
##    2  Nephron segment flux     — the Na+ waterfall, segment by segment
##    3  Drug PK & target engagement
##    4  Serum electrolytes & acid-base
##    5  Urinary excretion & the calcium sign flip
##    6  RAAS / macula densa / PGE2 axis
##    7  Clinical endpoints (growth, QTc, symptoms, eGFR)
##    8  Scenario comparison
##    9  Diagnostic test simulator (HCTZ / furosemide challenge)
##   10  Virtual population (genotype sweep)
##
##  Run with:
##    setwd("bartter-gitelman-syndrome"); shiny::runApp("bgs_shiny_app.R")
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)
library(ggrepel)

source("bgs_mrgsolve_model.R")

PAL <- c(TAL = "#2E6FA7", DCT = "#B4553F", NORM = "#4C8C5A",
         DRUG = "#7A5AA8", WARN = "#C0392B", GREY = "#7A7A7A")

theme_bgs <- function() {
  theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(colour = "grey35", size = 11),
          legend.position = "bottom",
          strip.text = element_text(face = "bold"))
}

## Reference ranges used for shading throughout the app
REF <- list(K = c(3.5, 5.1), MG = c(0.70, 1.05), HCO3 = c(22, 28),
            CL = c(98, 107), QTC = c(350, 450), UCACR = c(0.15, 0.45))

GENO <- bgs_genotypes()

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  title = "Bartter & Gitelman Syndrome — QSP Dashboard",
  tags$head(tags$style(HTML("
    .bgs-title { font-size: 25px; font-weight: 700; margin-bottom: 2px; }
    .bgs-sub   { color:#666; font-size: 13px; margin-bottom: 14px; }
    .keybox    { background:#F4F7FB; border-left:4px solid #2E6FA7;
                 padding:9px 13px; margin:9px 0; font-size:12.5px; }
    .warnbox   { background:#FDF3F2; border-left:4px solid #C0392B;
                 padding:9px 13px; margin:9px 0; font-size:12.5px; }
    .vbox      { background:#FFF; border:1px solid #E3E3E3; border-radius:6px;
                 padding:9px 11px; text-align:center; }
    .vnum      { font-size:21px; font-weight:700; }
    .vlab      { font-size:11px; color:#777; }
  "))),

  div(class = "bgs-title", "Bartter & Gitelman Syndrome — QSP Dashboard"),
  div(class = "bgs-sub",
      paste("Salt-losing tubulopathies. One transport lesion, two anatomical",
            "positions: FTAL (thick ascending limb) vs FDCT (distal",
            "convoluted tubule). 38-ODE mrgsolve model.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      selectInput("geno", "Genotype preset",
                  choices = GENO$genotype, selected = "Gitelman"),

      sliderInput("ftal", "FTAL — TAL NaCl capacity", 0.02, 1, 1, step = 0.01),
      sliderInput("fdct", "FDCT — DCT NCC capacity", 0.02, 1, 0.15, step = 0.01),
      helpText(HTML(paste0("<b>FTAL &lt; 1</b> = Bartter-type (upstream of the ",
                           "macula densa).<br><b>FDCT &lt; 1</b> = Gitelman-",
                           "type (downstream)."))),

      hr(),
      h5("Patient"),
      numericInput("bw", "Body weight (kg)", 60, 4, 120, step = 1),
      numericInput("htcm", "Height (cm)", 168, 50, 200, step = 1),
      numericInput("gfr0", "Baseline eGFR (mL/min/1.73m²)", 100, 15, 140,
                   step = 5),
      checkboxInput("peds", "Growing child (activate growth endpoint)", FALSE),

      hr(),
      h5("Diet"),
      sliderInput("nain", "Na intake (mmol/day)", 40, 400, 150, step = 10),
      sliderInput("kin", "K intake (mmol/day)", 20, 150, 70, step = 5),
      sliderInput("mgin", "Mg intake (mmol/day)", 4, 30, 12, step = 1),

      hr(),
      h5("Therapy"),
      sliderInput("indo", "Indomethacin (mg/kg/day, TID)", 0, 5, 0, step = 0.25),
      sliderInput("cele", "Celecoxib (mg/kg/day, BID)", 0, 12, 0, step = 0.5),
      sliderInput("amil", "Amiloride (mg/day, BID)", 0, 30, 0, step = 1),
      sliderInput("spiro", "Spironolactone (mg/day)", 0, 300, 0, step = 25),
      sliderInput("ena", "Enalapril (mg/day)", 0, 20, 0, step = 1),
      sliderInput("kcl", "Oral KCl (mmol/day, TID)", 0, 200, 0, step = 10),
      sliderInput("mgo", "Oral Mg (mmol/day)", 0, 60, 0, step = 2),
      sliderInput("mgn", "Mg doses per day", 1, 6, 4, step = 1),
      sliderInput("nacl", "NaCl supplement (mmol/day)", 0, 300, 0, step = 20),
      sliderInput("adher", "Adherence", 0.2, 1, 1, step = 0.05),
      checkboxInput("rhgh", "Recombinant GH", FALSE),

      hr(),
      sliderInput("days", "Simulation horizon (days)", 90, 1825, 730,
                  step = 30),
      actionButton("go", "Simulate", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      fluidRow(
        column(2, div(class = "vbox", div(class = "vnum", textOutput("vK")),
                      div(class = "vlab", "serum K (mmol/L)"))),
        column(2, div(class = "vbox", div(class = "vnum", textOutput("vMG")),
                      div(class = "vlab", "serum Mg (mmol/L)"))),
        column(2, div(class = "vbox", div(class = "vnum", textOutput("vHC")),
                      div(class = "vlab", "HCO3 (mmol/L)"))),
        column(2, div(class = "vbox", div(class = "vnum", textOutput("vCA")),
                      div(class = "vlab", "urine Ca/Cr"))),
        column(2, div(class = "vbox", div(class = "vnum", textOutput("vPG")),
                      div(class = "vlab", "urine PGE2 (xULN)"))),
        column(2, div(class = "vbox", div(class = "vnum", textOutput("vQT")),
                      div(class = "vlab", "QTc (ms)")))
      ),
      br(),

      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---- 1 -------------------------------------------------------------
        tabPanel(
          "1 · Patient & genotype",
          br(),
          uiOutput("lesionBox"),
          fluidRow(column(6, plotOutput("pLesion", height = 300)),
                   column(6, plotOutput("pCapacity", height = 300))),
          h5("Genotype library"),
          DTOutput("tGeno")
        ),

        ## ---- 2 -------------------------------------------------------------
        tabPanel(
          "2 · Nephron segment flux",
          br(),
          div(class = "keybox", HTML(paste(
            "The sodium waterfall. Each bar is what LEAVES the segment.",
            "A TAL lesion and a DCT lesion both raise distal delivery — but",
            "only the TAL lesion changes what the macula densa, which sits at",
            "the END of the TAL, is able to sense."))),
          plotOutput("pWaterfall", height = 330),
          fluidRow(column(6, plotOutput("pSegFrac", height = 290)),
                   column(6, plotOutput("pDeliv", height = 290)))
        ),

        ## ---- 3 -------------------------------------------------------------
        tabPanel(
          "3 · Drug PK & target engagement",
          br(),
          fluidRow(column(6, plotOutput("pPK", height = 320)),
                   column(6, plotOutput("pTE", height = 320))),
          div(class = "keybox", HTML(paste(
            "Target engagement is identical in both syndromes — the DRUG does",
            "the same thing. What differs is whether the engaged target sits",
            "on a loop the patient is actually running."))),
          plotOutput("pPKzoom", height = 280)
        ),

        ## ---- 4 -------------------------------------------------------------
        tabPanel(
          "4 · Serum electrolytes & acid-base",
          br(),
          plotOutput("pLytes", height = 430),
          div(class = "keybox", HTML(paste(
            "Alkalosis corrects only when CHLORIDE is repleted: chloride",
            "depletion raises the renal bicarbonate threshold. Try KCl vs",
            "an equal K load without chloride — the model separates them."))),
          fluidRow(column(6, plotOutput("pKcompart", height = 280)),
                   column(6, plotOutput("pBic", height = 280)))
        ),

        ## ---- 5 -------------------------------------------------------------
        tabPanel(
          "5 · Urine & the calcium sign flip",
          br(),
          div(class = "keybox", HTML(paste(
            "Urinary Ca/Cr changes SIGN with lesion position. TAL lesions",
            "collapse the lumen-positive potential (hypercalciuria →",
            "nephrocalcinosis); DCT lesions leave it intact and volume",
            "contraction drives proximal Ca retention (hypocalciuria)."))),
          fluidRow(column(6, plotOutput("pUCa", height = 320)),
                   column(6, plotOutput("pPD", height = 320))),
          plotOutput("pUrine", height = 330)
        ),

        ## ---- 6 -------------------------------------------------------------
        tabPanel(
          "6 · RAAS / macula densa / PGE2",
          br(),
          div(class = "keybox", HTML(paste(
            "The macula densa reads luminal NaCl through its OWN NKCC2, which",
            "carries the same lesion. So a TAL lesion blinds it even though",
            "luminal NaCl is HIGH — exactly as furosemide does. That is why",
            "urinary PGE2 rises in Bartter and not in Gitelman."))),
          fluidRow(column(6, plotOutput("pMD", height = 320)),
                   column(6, plotOutput("pRAAS", height = 320))),
          plotOutput("pPGE", height = 300)
        ),

        ## ---- 7 -------------------------------------------------------------
        tabPanel(
          "7 · Clinical endpoints",
          br(),
          fluidRow(column(6, plotOutput("pGrowth", height = 300)),
                   column(6, plotOutput("pQTC", height = 300))),
          fluidRow(column(6, plotOutput("pSympt", height = 300)),
                   column(6, plotOutput("pKidney", height = 300))),
          uiOutput("safetyBox")
        ),

        ## ---- 8 -------------------------------------------------------------
        tabPanel(
          "8 · Scenario comparison",
          br(),
          checkboxGroupInput(
            "scen", "Scenarios to run",
            choices = setNames(names(bgs_scenarios()),
                               vapply(bgs_scenarios(), `[[`, "", "label")),
            selected = c("S2", "S3", "S5", "S6"), inline = FALSE),
          actionButton("goScen", "Run selected scenarios",
                       class = "btn-primary"),
          br(), br(),
          plotOutput("pScen", height = 480),
          h5("Endpoints at end of horizon"),
          DTOutput("tScen")
        ),

        ## ---- 9 -------------------------------------------------------------
        tabPanel(
          "9 · Diagnostic test simulator",
          br(),
          div(class = "keybox", HTML(paste(
            "A thiazide or furosemide challenge is an in-vivo target-engagement",
            "assay. The response is FLAT when the transporter the drug blocks",
            "is already offline — thiazide in Gitelman, furosemide in",
            "Bartter I/II."))),
          fluidRow(
            column(4, radioButtons("dxdrug", "Challenge",
                                   c("Hydrochlorothiazide" = "hctz",
                                     "Furosemide" = "furo"), "hctz")),
            column(4, sliderInput("dxblock", "Fractional block", 0.3, 0.95,
                                  0.90, step = 0.05)),
            column(4, actionButton("goDx", "Run challenge",
                                   class = "btn-primary"))
          ),
          plotOutput("pDx", height = 400),
          DTOutput("tDx")
        ),

        ## ---- 10 ------------------------------------------------------------
        tabPanel(
          "10 · Virtual population",
          br(),
          actionButton("goVpop", "Sweep the genotype library",
                       class = "btn-primary"),
          br(), br(),
          plotOutput("pVpop", height = 400),
          DTOutput("tVpop"),
          div(class = "warnbox", HTML(paste(
            "<b>Educational / research model only.</b> Not validated for",
            "clinical decision-making, prescribing, or regulatory use.")))
        )
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  ## --- genotype preset fills the two capacity sliders -----------------------
  observeEvent(input$geno, {
    g <- GENO[GENO$genotype == input$geno, ]
    updateSliderInput(session, "ftal", value = g$FTAL)
    updateSliderInput(session, "fdct", value = g$FDCT)
    updateCheckboxInput(session, "peds", value = as.logical(g$PEDS))
  })

  pars <- reactive({
    g <- GENO[GENO$genotype == input$geno, ]
    list(FTAL = input$ftal, FDCT = input$fdct,
         FROMKCD = g$FROMKCD, PEDS = as.numeric(input$peds),
         BW = input$bw, HTCM = input$htcm, GFR0 = input$gfr0,
         NAIN = input$nain, KIN = input$kin, MGIN = input$mgin,
         ADHER0 = input$adher, RHGH = if (input$rhgh) 0.6 else 0)
  })

  regimen <- reactive({
    d <- input$days
    ev_list <- list()
    if (input$indo  > 0) ev_list <- c(ev_list, list(
      dose_indomethacin(input$indo, input$bw, 3, 0, d)))
    if (input$cele  > 0) ev_list <- c(ev_list, list(
      dose_celecoxib(input$cele, input$bw, 2, 0, d)))
    if (input$amil  > 0) ev_list <- c(ev_list, list(
      dose_amiloride(input$amil, 2, 0, d)))
    if (input$spiro > 0) ev_list <- c(ev_list, list(
      dose_spironolactone(input$spiro, 1, 0, d)))
    if (input$ena   > 0) ev_list <- c(ev_list, list(
      dose_enalapril(input$ena, 0, d)))
    if (input$kcl   > 0) ev_list <- c(ev_list, list(
      dose_kcl(input$kcl, 3, 0, d)))
    if (input$mgo   > 0) ev_list <- c(ev_list, list(
      dose_mg(input$mgo, input$mgn, 0, d)))
    if (input$nacl  > 0) ev_list <- c(ev_list, list(
      dose_nacl(input$nacl, 3, 0, d)))
    if (!length(ev_list)) return(NULL)
    Reduce(c, ev_list)
  })

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    withProgress(message = "Burning in to the genotype steady state…", {
      m <- bgs_init_at(bgs, pars(), days = 600)
      incProgress(0.6, detail = "Simulating therapy…")
      rg <- regimen()
      out <- if (is.null(rg)) {
        mrgsim(m, end = input$days, delta = 1, hmax = 0.25)
      } else {
        mrgsim(m, events = rg, end = input$days, delta = 1, hmax = 0.25)
      }
      as.data.frame(out)
    })
  })

  last <- reactive({ d <- sim(); d[nrow(d), ] })

  ## --- value boxes ----------------------------------------------------------
  fmt <- function(x, n = 2) formatC(x, digits = n, format = "f")
  output$vK  <- renderText(fmt(last()$CK))
  output$vMG <- renderText(fmt(last()$CMG))
  output$vHC <- renderText(fmt(last()$CHC, 1))
  output$vCA <- renderText(fmt(last()$UCACR, 3))
  output$vPG <- renderText(fmt(last()$UPGE2, 1))
  output$vQT <- renderText(fmt(last()$QTC, 0))

  ## ========================================================================
  ##  TAB 1 — lesion position
  ## ========================================================================
  output$lesionBox <- renderUI({
    ftal <- input$ftal; fdct <- input$fdct
    msg <- if (ftal < 0.7 && fdct < 0.7) {
      paste0("Mixed lesion (FTAL ", ftal, ", FDCT ", fdct, "): the ",
             "CLCNKB/BSND pattern — ClC-Kb is expressed in BOTH segments, ",
             "which is exactly why type III overlaps Gitelman.")
    } else if (ftal < 0.7) {
      paste0("Lesion UPSTREAM of the macula densa (FTAL ", ftal, "): expect a ",
             "PGE2 storm, hypercalciuria, polyuria, growth failure — and a ",
             "real response to COX inhibition.")
    } else if (fdct < 0.7) {
      paste0("Lesion DOWNSTREAM of the macula densa (FDCT ", fdct, "): expect ",
             "NORMAL urinary PGE2, hypocalciuria, severe hypomagnesaemia — ",
             "and little benefit from COX inhibition.")
    } else {
      "Both capacities near normal — control physiology."
    }
    div(class = "keybox", HTML(paste0("<b>Lesion position:</b> ", msg)))
  })

  output$pLesion <- renderPlot({
    d <- data.frame(
      seg = factor(c("PT", "TAL", "MD", "DCT", "ASDN"),
                   levels = c("PT", "TAL", "MD", "DCT", "ASDN")),
      x = 1:5,
      cap = c(1, input$ftal, NA, input$fdct, 1))
    ggplot(d, aes(x, 1)) +
      geom_segment(aes(x = 0.6, xend = 5.4, y = 1, yend = 1),
                   linewidth = 3, colour = "grey85") +
      geom_point(aes(colour = cap, size = ifelse(is.na(cap), 5, 12))) +
      geom_text(aes(label = seg), vjust = -2.4, fontface = "bold") +
      geom_text(aes(label = ifelse(is.na(cap), "sensor",
                                   paste0(round(cap * 100), "%"))),
                vjust = 3.4, size = 4) +
      scale_colour_gradient(low = PAL[["WARN"]], high = PAL[["NORM"]],
                            limits = c(0, 1), na.value = "#E8A33D",
                            name = "capacity") +
      scale_size_identity() +
      coord_cartesian(ylim = c(0.6, 1.4)) +
      labs(title = "Where the lesion sits",
           subtitle = "MD = macula densa, at the END of the TAL") +
      theme_bgs() +
      theme(axis.text = element_blank(), axis.title = element_blank(),
            panel.grid = element_blank())
  })

  output$pCapacity <- renderPlot({
    d <- sim()
    dd <- d %>% select(time, ATAL, ANCC, PDREL) %>%
      pivot_longer(-time)
    lab <- c(ATAL = "effective TAL activity", ANCC = "effective NCC activity",
             PDREL = "TAL lumen-positive PD")
    ggplot(dd, aes(time, value, colour = name)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = c(ATAL = PAL[["TAL"]], ANCC = PAL[["DCT"]],
                                     PDREL = PAL[["GREY"]]),
                          labels = lab, name = NULL) +
      labs(title = "Effective transport activity", x = "day",
           y = "relative to healthy") +
      theme_bgs()
  })

  output$tGeno <- renderDT(datatable(GENO, rownames = FALSE,
                                     options = list(dom = "t", pageLength = 12)))

  ## ========================================================================
  ##  TAB 2 — segment flux waterfall
  ## ========================================================================
  output$pWaterfall <- renderPlot({
    r <- last()
    d <- data.frame(
      stage = factor(c("filtered", "left PT", "left TAL (= at MD)",
                       "left DCT (= to ASDN)", "urine"),
                     levels = c("filtered", "left PT", "left TAL (= at MD)",
                                "left DCT (= to ASDN)", "urine")),
      load = c(r$FILTNA, r$FILTNA - r$RNAPT, r$LDCT, r$LASDN, r$UNA))
    ggplot(d, aes(stage, load)) +
      geom_col(fill = PAL[["TAL"]], width = 0.62) +
      geom_text(aes(label = round(load)), vjust = -0.5, size = 4.2) +
      scale_y_log10() +
      labs(title = "Sodium waterfall along the nephron",
           subtitle = "mmol/day leaving each segment (log scale)",
           x = NULL, y = "Na+ load (mmol/day)") +
      theme_bgs()
  })

  output$pSegFrac <- renderPlot({
    d <- sim()
    dd <- d %>%
      transmute(time,
                PT = RNAPT / FILTNA, TAL = RNATAL / FILTNA,
                DCT = RNADCT / FILTNA, ASDN = RNACD / FILTNA,
                urine = UNA / FILTNA) %>%
      pivot_longer(-time)
    ggplot(dd, aes(time, value, fill = name)) +
      geom_area() +
      scale_fill_brewer(palette = "Blues", direction = -1, name = NULL) +
      labs(title = "Fractional Na+ handling", x = "day",
           y = "fraction of filtered load") +
      theme_bgs()
  })

  output$pDeliv <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_line(aes(y = DELIVREL, colour = "distal Na delivery"),
                linewidth = 0.9) +
      geom_line(aes(y = LUMREL, colour = "luminal NaCl at MD"),
                linewidth = 0.9) +
      geom_line(aes(y = MDSENSE, colour = "SENSED NaCl at MD"),
                linewidth = 1.2) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
      scale_colour_manual(values = c("distal Na delivery" = PAL[["DCT"]],
                                     "luminal NaCl at MD" = PAL[["GREY"]],
                                     "SENSED NaCl at MD" = PAL[["WARN"]]),
                          name = NULL) +
      labs(title = "Delivered vs sensed",
           subtitle = paste("Luminal NaCl at the macula densa can be HIGH",
                            "while the sensed signal is LOW"),
           x = "day", y = "relative to healthy") +
      theme_bgs()
  })

  ## ========================================================================
  ##  TAB 3 — PK / target engagement
  ## ========================================================================
  output$pPK <- renderPlot({
    d <- sim()
    dd <- d %>% select(time, CIND, CCEL, CAMI, CCAN, CACE) %>%
      pivot_longer(-time) %>% filter(value > 1e-9)
    if (!nrow(dd)) return(ggplot() + labs(title = "No drug on board") +
                            theme_bgs())
    ggplot(dd, aes(time, value, colour = name)) +
      geom_line(linewidth = 0.7) +
      scale_y_log10() +
      scale_colour_brewer(palette = "Dark2", name = NULL,
                          labels = c(CACE = "enalaprilat", CAMI = "amiloride",
                                     CCAN = "canrenone", CCEL = "celecoxib",
                                     CIND = "indomethacin")) +
      labs(title = "Plasma concentrations", x = "day", y = "mg/L (log)") +
      theme_bgs()
  })

  output$pTE <- renderPlot({
    d <- sim()
    dd <- d %>% select(time, INH2, INH1, AMIB, MRB) %>% pivot_longer(-time)
    ggplot(dd, aes(time, 100 * value, colour = name)) +
      geom_line(linewidth = 0.9) +
      scale_colour_manual(values = c(INH2 = PAL[["WARN"]], INH1 = PAL[["GREY"]],
                                     AMIB = PAL[["DRUG"]], MRB = PAL[["DCT"]]),
                          labels = c(AMIB = "ENaC block", INH1 = "COX-1 (gastric)",
                                     INH2 = "renal COX-2", MRB = "MR block"),
                          name = NULL) +
      ylim(0, 100) +
      labs(title = "Target engagement", x = "day", y = "% inhibition") +
      theme_bgs()
  })

  output$pPKzoom <- renderPlot({
    d <- sim()
    w <- d %>% filter(time >= max(time) - 3)
    dd <- w %>% select(time, CIND, CCEL, CAMI) %>% pivot_longer(-time) %>%
      filter(value > 1e-9)
    if (!nrow(dd)) return(ggplot() + labs(title = "No drug on board") +
                            theme_bgs())
    ggplot(dd, aes(time, value, colour = name)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_brewer(palette = "Dark2", guide = "none") +
      labs(title = "Steady-state dosing interval (last 3 days)",
           x = "day", y = "mg/L") +
      theme_bgs()
  })

  ## ========================================================================
  ##  TAB 4 — serum electrolytes & acid-base
  ## ========================================================================
  output$pLytes <- renderPlot({
    d <- sim()
    dd <- d %>%
      select(time, `serum K` = CK, `serum Mg` = CMG, `serum Cl` = CCL,
             `serum HCO3` = CHC, `serum Na` = CNA) %>%
      pivot_longer(-time)
    rng <- data.frame(
      name = c("serum K", "serum Mg", "serum Cl", "serum HCO3", "serum Na"),
      lo = c(REF$K[1], REF$MG[1], REF$CL[1], REF$HCO3[1], 135),
      hi = c(REF$K[2], REF$MG[2], REF$CL[2], REF$HCO3[2], 145))
    ggplot(dd, aes(time, value)) +
      geom_rect(data = rng, inherit.aes = FALSE,
                aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi),
                fill = PAL[["NORM"]], alpha = 0.12) +
      geom_line(linewidth = 0.9, colour = PAL[["TAL"]]) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "Serum electrolytes and acid-base",
           subtitle = "shaded band = reference range",
           x = "day", y = NULL) +
      theme_bgs()
  })

  output$pKcompart <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_line(aes(y = CK, colour = "extracellular"), linewidth = 1) +
      geom_line(aes(y = CKIC / 29.3, colour = "intracellular / 29.3"),
                linewidth = 1) +
      scale_colour_manual(values = c(extracellular = PAL[["TAL"]],
                                     `intracellular / 29.3` = PAL[["DCT"]]),
                          name = NULL) +
      labs(title = "Potassium: the hidden intracellular deficit",
           subtitle = "serum K understates total-body depletion",
           x = "day", y = "mmol/L") +
      theme_bgs()
  })

  output$pBic <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_line(aes(y = CHC, colour = "serum HCO3"), linewidth = 1) +
      geom_line(aes(y = CCL - 70, colour = "serum Cl − 70"), linewidth = 1) +
      scale_colour_manual(values = c(`serum HCO3` = PAL[["WARN"]],
                                     `serum Cl − 70` = PAL[["TAL"]]),
                          name = NULL) +
      labs(title = "Chloride depletion drives the alkalosis",
           subtitle = "HCO3 falls only as Cl is repleted",
           x = "day", y = "mmol/L") +
      theme_bgs()
  })

  ## ========================================================================
  ##  TAB 5 — urine and the calcium sign flip
  ## ========================================================================
  output$pUCa <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, UCACR)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = REF$UCACR[1],
               ymax = REF$UCACR[2], fill = PAL[["NORM"]], alpha = 0.12) +
      geom_hline(yintercept = 0.07, linetype = 2, colour = PAL[["DCT"]]) +
      geom_hline(yintercept = 0.20, linetype = 2, colour = PAL[["TAL"]]) +
      geom_line(linewidth = 1.1, colour = PAL[["WARN"]]) +
      annotate("text", x = Inf, y = 0.07, hjust = 1.05, vjust = -0.5,
               label = "Gitelman < 0.07", size = 3.4, colour = PAL[["DCT"]]) +
      annotate("text", x = Inf, y = 0.20, hjust = 1.05, vjust = -0.5,
               label = "Bartter > 0.20", size = 3.4, colour = PAL[["TAL"]]) +
      labs(title = "Urinary Ca/creatinine — the discriminator",
           x = "day", y = "mmol/mmol") +
      theme_bgs()
  })

  output$pPD <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_line(aes(y = PDREL, colour = "TAL lumen-positive PD"),
                linewidth = 1) +
      geom_line(aes(y = UCA / max(d$UCA), colour = "urinary Ca (scaled)"),
                linewidth = 1) +
      scale_colour_manual(values = c(`TAL lumen-positive PD` = PAL[["TAL"]],
                                     `urinary Ca (scaled)` = PAL[["WARN"]]),
                          name = NULL) +
      labs(title = "Lose the potential, lose the calcium",
           x = "day", y = "relative") +
      theme_bgs()
  })

  output$pUrine <- renderPlot({
    d <- sim()
    dd <- d %>%
      select(time, `U Na (mmol/d)` = UNA, `U K (mmol/d)` = UK,
             `U Cl (mmol/d)` = UCL, `U Ca (mmol/d)` = UCA,
             `U Mg (mmol/d)` = UMG, `FE-Mg (%)` = FEMG,
             `urine volume (L/d)` = UVOL,
             `urine osm (mosm/kg)` = UOSMA) %>%
      pivot_longer(-time)
    ggplot(dd, aes(time, value)) +
      geom_line(linewidth = 0.85, colour = PAL[["DCT"]]) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(title = "Urinary excretion", x = "day", y = NULL) +
      theme_bgs()
  })

  ## ========================================================================
  ##  TAB 6 — macula densa / RAAS / PGE2
  ## ========================================================================
  output$pMD <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
      geom_line(aes(y = LUMREL, colour = "luminal NaCl (delivered)"),
                linewidth = 1) +
      geom_line(aes(y = MDSENSE, colour = "SENSED NaCl"), linewidth = 1.3) +
      scale_colour_manual(values = c(`luminal NaCl (delivered)` = PAL[["GREY"]],
                                     `SENSED NaCl` = PAL[["WARN"]]),
                          name = NULL) +
      labs(title = "The positional sensor",
           subtitle = "the gap between these two lines IS Bartter syndrome",
           x = "day", y = "relative to healthy") +
      theme_bgs()
  })

  output$pRAAS <- renderPlot({
    d <- sim()
    dd <- d %>%
      transmute(time, PRA = PRA, `Ang II` = ANGREL, aldosterone = ALDREL) %>%
      pivot_longer(-time)
    ggplot(dd, aes(time, value, colour = name)) +
      geom_line(linewidth = 0.95) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey55") +
      scale_colour_manual(values = c(PRA = PAL[["WARN"]],
                                     `Ang II` = PAL[["TAL"]],
                                     aldosterone = PAL[["DCT"]]), name = NULL) +
      labs(title = "Secondary hyperreninaemic hyperaldosteronism",
           subtitle = "…with a normal blood pressure",
           x = "day", y = "x baseline") +
      theme_bgs()
  })

  output$pPGE <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, UPGE2)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1.5,
               fill = PAL[["NORM"]], alpha = 0.12) +
      geom_line(linewidth = 1.2, colour = PAL[["WARN"]]) +
      labs(title = "Urinary PGE2 — the hyperprostaglandin-E phenotype",
           subtitle = paste("Elevated only when the lesion is in the TAL.",
                            "This is the loop COX inhibitors act on."),
           x = "day", y = "x upper limit of normal") +
      theme_bgs()
  })

  ## ========================================================================
  ##  TAB 7 — clinical endpoints
  ## ========================================================================
  output$pGrowth <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_line(aes(y = HTSD, colour = "height SDS"), linewidth = 1.1) +
      geom_line(aes(y = IGF1 * 2 - 2, colour = "IGF-1 (rescaled)"),
                linewidth = 0.9) +
      geom_hline(yintercept = -2, linetype = 2, colour = PAL[["WARN"]]) +
      scale_colour_manual(values = c(`height SDS` = PAL[["NORM"]],
                                     `IGF-1 (rescaled)` = PAL[["GREY"]]),
                          name = NULL) +
      labs(title = "Growth", subtitle = "dashed line = −2 SDS",
           x = "day", y = "SDS") +
      theme_bgs()
  })

  output$pQTC <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, QTC)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 450,
               fill = PAL[["NORM"]], alpha = 0.12) +
      geom_line(linewidth = 1.1, colour = PAL[["WARN"]]) +
      geom_hline(yintercept = 500, linetype = 2, colour = PAL[["WARN"]]) +
      labs(title = "QTc", subtitle = "driven jointly by K and Mg",
           x = "day", y = "ms") +
      theme_bgs()
  })

  output$pSympt <- renderPlot({
    d <- sim()
    dd <- d %>% select(time, cramps = CRAMP, fatigue = FATIG) %>%
      pivot_longer(-time)
    ggplot(dd, aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c(cramps = PAL[["DCT"]],
                                     fatigue = PAL[["GREY"]]), name = NULL) +
      ylim(0, 10) +
      labs(title = "Symptom burden", x = "day", y = "score (0-10)") +
      theme_bgs()
  })

  output$pKidney <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_line(aes(y = EGFR, colour = "eGFR"), linewidth = 1.1) +
      geom_line(aes(y = 30 * NEPH + 60, colour = "nephrocalcinosis (rescaled)"),
                linewidth = 0.9) +
      scale_colour_manual(values = c(eGFR = PAL[["TAL"]],
                                     `nephrocalcinosis (rescaled)` =
                                       PAL[["WARN"]]), name = NULL) +
      labs(title = "Kidney outcome",
           subtitle = "nephrocalcinosis grade rescaled: 60 + 30 x grade",
           x = "day", y = "mL/min/1.73m²") +
      theme_bgs()
  })

  output$safetyBox <- renderUI({
    r <- last()
    msgs <- character()
    if (r$ULCER > 0.30)
      msgs <- c(msgs, sprintf(paste("NSAID gastropathy risk index %.2f —",
                                    "the commonest reason therapy is stopped."),
                              r$ULCER))
    if (r$AKIFLAG > 0.5)
      msgs <- c(msgs, "PG-dependent GFR loss > 20%: acute kidney injury risk.")
    if (r$CK > 5.3)
      msgs <- c(msgs, sprintf(paste("Iatrogenic HYPERKALAEMIA (K %.2f) —",
                                    "stacked K-sparing agents + KCl."), r$CK))
    if (r$QTC > 500)
      msgs <- c(msgs, sprintf("QTc %.0f ms: torsade risk.", r$QTC))
    if (!length(msgs))
      return(div(class = "keybox", "No safety flags at the current settings."))
    div(class = "warnbox", HTML(paste0("<b>Safety flags</b><ul><li>",
                                       paste(msgs, collapse = "</li><li>"),
                                       "</li></ul>")))
  })

  ## ========================================================================
  ##  TAB 8 — scenario comparison
  ## ========================================================================
  scenData <- eventReactive(input$goScen, {
    req(length(input$scen) > 0)
    withProgress(message = "Running scenarios…", {
      out <- lapply(seq_along(input$scen), function(i) {
        incProgress(1 / length(input$scen), detail = input$scen[i])
        bgs_run_scenario(input$scen[i], days = input$days, delta = 7)
      })
      bind_rows(out)
    })
  })

  output$pScen <- renderPlot({
    d <- scenData()
    dd <- d %>%
      select(time, scenario, `serum K` = CK, `serum Mg` = CMG,
             `serum HCO3` = CHC, `urine Ca/Cr` = UCACR,
             `urine PGE2` = UPGE2, `QTc` = QTC) %>%
      pivot_longer(c(-time, -scenario))
    ggplot(dd, aes(time, value, colour = scenario)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_colour_brewer(palette = "Dark2", name = NULL) +
      labs(title = "Scenario comparison", x = "day", y = NULL) +
      theme_bgs() +
      theme(legend.position = "bottom", legend.direction = "vertical")
  })

  output$tScen <- renderDT({
    d <- scenData()
    tab <- d %>% group_by(scenario) %>% filter(time == max(time)) %>%
      ungroup() %>%
      transmute(scenario,
                K = round(CK, 2), Mg = round(CMG, 2), HCO3 = round(CHC, 1),
                Cl = round(CCL, 1), `uCa/Cr` = round(UCACR, 3),
                `FE-Mg %` = round(FEMG, 1), `uPGE2` = round(UPGE2, 2),
                `Uvol L/d` = round(UVOL, 2), QTc = round(QTC),
                cramps = round(CRAMP, 1), `Ht SDS` = round(HTSD, 2),
                nephrocalc = round(NEPH, 2), eGFR = round(EGFR, 1),
                `ulcer risk` = round(ULCER, 2))
    datatable(tab, rownames = FALSE,
              options = list(dom = "t", scrollX = TRUE, pageLength = 20))
  })

  ## ========================================================================
  ##  TAB 9 — diagnostic challenge
  ## ========================================================================
  dxData <- eventReactive(input$goDx, {
    withProgress(message = "Running the challenge in every genotype…", {
      blk <- input$dxblock
      p0 <- list(DXT0 = 2, DXDUR = 0.25)
      if (input$dxdrug == "hctz") p0$DXHCTZ <- blk else p0$DXFURO <- blk
      out <- lapply(seq_len(nrow(GENO)), function(i) {
        incProgress(1 / nrow(GENO), detail = GENO$genotype[i])
        pg <- list(FTAL = GENO$FTAL[i], FDCT = GENO$FDCT[i],
                   FROMKCD = GENO$FROMKCD[i], PEDS = GENO$PEDS[i])
        m <- bgs_init_at(bgs, c(pg, p0), days = 600)
        df <- as.data.frame(mrgsim(m, end = 5, delta = 0.05, hmax = 0.01))
        df$genotype <- GENO$genotype[i]
        df
      })
      bind_rows(out)
    })
  })

  output$pDx <- renderPlot({
    d <- dxData()
    dd <- d %>% select(time, genotype, `urinary Cl (mmol/d)` = UCL,
                       `urinary Na (mmol/d)` = UNA,
                       `urine volume (L/d)` = UVOL) %>%
      pivot_longer(c(-time, -genotype))
    ggplot(dd, aes(time, value, colour = genotype)) +
      geom_vline(xintercept = 2, linetype = 2, colour = "grey45") +
      geom_line(linewidth = 0.85) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(title = paste0("Response to a ",
                          ifelse(input$dxdrug == "hctz", "thiazide",
                                 "furosemide"), " challenge"),
           subtitle = paste("A FLAT response means the transporter the drug",
                            "blocks is already offline"),
           x = "day", y = NULL) +
      theme_bgs() +
      theme(legend.direction = "vertical")
  })

  output$tDx <- renderDT({
    d <- dxData()
    tab <- d %>% group_by(genotype) %>%
      summarise(base_UCl = round(UCL[which.min(abs(time - 1.9))], 1),
                peak_UCl = round(max(UCL[time >= 2 & time <= 2.5]), 1),
                .groups = "drop") %>%
      mutate(delta = round(peak_UCl - base_UCl, 1),
             `% rise` = round(100 * delta / base_UCl, 1),
             interpretation = ifelse(`% rise` < 15,
                                     "FLAT — target already offline",
                                     "responsive"))
    datatable(tab, rownames = FALSE,
              options = list(dom = "t", pageLength = 12))
  })

  ## ========================================================================
  ##  TAB 10 — virtual population
  ## ========================================================================
  vpopData <- eventReactive(input$goVpop, {
    withProgress(message = "Sweeping the genotype library…", {
      bgs_vpop(days = 600)
    })
  })

  output$pVpop <- renderPlot({
    d <- vpopData()
    ggplot(d, aes(uCaCr, uPGE2, colour = Mg, size = 4 - K)) +
      geom_vline(xintercept = c(0.07, 0.20), linetype = 2,
                 colour = "grey55") +
      geom_hline(yintercept = 1.5, linetype = 2, colour = "grey55") +
      geom_point(alpha = 0.9) +
      geom_text_repel(aes(label = genotype), size = 3.4,
                               show.legend = FALSE) +
      scale_x_log10() +
      scale_colour_gradient(low = PAL[["WARN"]], high = PAL[["NORM"]],
                            name = "serum Mg") +
      scale_size_continuous(name = "K deficit") +
      labs(title = "The two-dimensional phenotype space",
           subtitle = paste("Genotypes separate along urinary Ca/Cr (x) and",
                            "urinary PGE2 (y) with no genotype-specific",
                            "parameters"),
           x = "urinary Ca/creatinine (mmol/mmol, log)",
           y = "urinary PGE2 (x ULN)") +
      theme_bgs()
  })

  output$tVpop <- renderDT(datatable(vpopData(), rownames = FALSE,
                                     options = list(dom = "t", scrollX = TRUE,
                                                    pageLength = 12)))
}

shinyApp(ui, server)
