## =====================================================================
##  pws_shiny_app.R
##  Prader-Willi Syndrome QSP model — interactive dashboard
##  프래더-윌리 증후군 QSP 모델 — 인터랙티브 대시보드
##
##  10 tabs:
##    1  Patient & environment   — the lesion and the four food knobs
##    2  The lesion              — PC1/3 branches, Eq. A / Eq. B, live
##    3  Drug PK                 — apparent day-scale exposures
##    4  Satiety integrator      — the five arms and the harmonic mean
##    5  Bistability             — the phase line, its fixed points, the cusp
##    6  Body composition        — energy balance, lean target, fat buffer
##    7  Growth & somatotropic   — height, IGF-1 SDS, bone age
##    8  Airway                  — THE TWO CLOCKS and the GH window
##    9  Endpoints               — HQ-CT, %fat, AHI, HbA1c, BMD, Cobb
##   10  Scenario comparison     — the drug panel, biomarker vs endpoint
##
##  DESIGN NOTE
##  -----------
##  Most disease dashboards plot endpoints.  This one plots the MECHANISM
##  next to the endpoint on purpose, because the whole content of this
##  model is that the two are nearly orthogonal: tab 4 shows which satiety
##  arm a drug reaches, tab 5 shows whether the patient can leave the
##  food-seeking state at all, and tab 10 shows that the agents which move
##  the ghrelin biomarker most move HQ-CT least.  A dashboard that showed
##  only tab 9 would let a user conclude the opposite of what the model
##  says.
##
##  RUN
##    library(shiny); library(mrgsolve); library(dplyr); library(tidyr)
##    library(ggplot2)
##    shiny::runApp("pws_shiny_app.R")
##
##  The model file `pws_mrgsolve_model.R` must sit in the same directory.
##  Every default here reproduces a row of `pws_reference_output.txt`.
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

MODFILE <- "pws_mrgsolve_model.R"
YEARS <- function(y) y * 365.25

mod <- mread_cache("pws", MODFILE)

## ---- the four food environments, exactly as in the reference model ----
ENVS <- list(
  "non-PWS control (self-regulated)" = list(AVAIL = 1.40, CUE = 1.00, TITR = 0),
  "PWS, unmanaged (free access)"     = list(AVAIL = 1.40, CUE = 1.00, TITR = 1e-9),
  "PWS, normal portions, no titration" = list(AVAIL = 1.00, CUE = 0.55, TITR = 1e-9),
  "PWS, standard management"         = list(AVAIL = 1.00, CUE = 0.28, TITR = 1),
  "PWS, strict environment"          = list(AVAIL = 0.90, CUE = 0.10, TITR = 1)
)

## ---- the drug panel ---------------------------------------------------
DRUGS <- list(
  "none"                       = NULL,
  "somatropin 0.035 mg/kg/d"   = list(cmt = "AGHD", amt = NA, ii = 1,   mgkg = 0.035),
  "long-acting GH weekly"      = list(cmt = "ALGD", amt = NA, ii = 7,   mgkg = 0.245),
  "carbetocin 3.2 mg TID"      = list(cmt = "ACBD", amt = 3.2,  ii = 1/3),
  "carbetocin 9.6 mg TID"      = list(cmt = "ACBD", amt = 9.6,  ii = 1/3),
  "DCCR 5.1 mg/kg/d"           = list(cmt = "ADZD", amt = NA, ii = 1,   mgkg = 5.1),
  "setmelanotide 3 mg/d"       = list(cmt = "ASMD", amt = 3.0,  ii = 1),
  "semaglutide 2.4 mg weekly"  = list(cmt = "ASGD", amt = 2.4,  ii = 7),
  "octreotide LAR 30 mg q28d"  = list(cmt = "AOCD", amt = 30,   ii = 28),
  "livoletide 60 ug/kg/d"      = list(cmt = "ALVD", amt = NA, ii = 1,   mgkg = 0.060),
  "testosterone 100 mg q14d"   = list(cmt = "ATSD", amt = 100,  ii = 14),
  "metformin 1500 mg/d"        = list(cmt = "AMFD", amt = 1500, ii = 1)
)

## a body weight only used to turn mg/kg into mg for the event table;
## the model itself always normalizes concentrations by the live weight
make_ev <- function(name, start_y, stop_y, wt_kg) {
  d <- DRUGS[[name]]
  if (is.null(d)) return(NULL)
  amt <- if (is.na(d$amt)) d$mgkg * wt_kg else d$amt
  n <- max(1, floor((YEARS(stop_y) - YEARS(start_y)) / d$ii))
  ev(time = YEARS(start_y), amt = amt, cmt = d$cmt, ii = d$ii, addl = n - 1)
}

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom")

## =====================================================================
##  UI
## =====================================================================
ui <- fluidPage(
  titlePanel("Prader-Willi Syndrome — QSP model dashboard"),
  p(tags$i(paste("one convertase, five branches, a harmonic-mean satiety",
                 "integrator, and a bistable food-seeking state"))),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("les", "Genotype",
                  c("PWS (paternal 15q11-q13 silent)" = 1,
                    "non-PWS control" = 0), selected = 1),
      selectInput("subtype", "Subtype (behaviour modifier)",
                  c("deletion (FSUB 1.15)" = 1.15,
                    "mUPD (FSUB 0.90)" = 0.90,
                    "neutral (1.00)" = 1.00), selected = 1.15),
      sliderInput("dpc13", "PC1/3 activity loss (DPC13)",
                  0.0, 0.85, 0.60, step = 0.05),
      sliderInput("foxtn", "Surviving PVN oxytocin neurons (FOXTN)",
                  0.4, 1.0, 0.75, step = 0.05),
      sliderInput("fosa", "Upper-airway vulnerability (FOSA)",
                  0.6, 3.0, 1.0, step = 0.1),
      hr(),
      h4("Food environment"),
      selectInput("env", "Preset", names(ENVS),
                  selected = "PWS, standard management"),
      checkboxInput("envman", "override the preset", FALSE),
      conditionalPanel(
        "input.envman == true",
        sliderInput("avail", "portions offered / requirement (AVAIL)",
                    0.5, 1.4, 1.0, step = 0.02),
        sliderInput("cue", "food cue & access exposure (CUE)",
                    0.05, 1.0, 0.28, step = 0.01),
        sliderInput("titr", "caregiver titration to the growth chart (TITR)",
                    0, 1, 1, step = 0.1),
        sliderInput("bwtgtr", "weight-for-height target held (BWTGTR)",
                    0.9, 1.6, 1.15, step = 0.05)
      ),
      hr(),
      h4("Therapy"),
      selectInput("drug1", "Drug 1", names(DRUGS),
                  selected = "somatropin 0.035 mg/kg/d"),
      sliderInput("d1start", "start (years)", 0.5, 20, 1, step = 0.5),
      selectInput("drug2", "Drug 2 (add-on)", names(DRUGS),
                  selected = "none"),
      sliderInput("d2start", "start (years)", 0.5, 25, 12, step = 0.5),
      hr(),
      sliderInput("tend", "simulate to age (years)", 2, 30, 25, step = 1),
      numericInput("wt", "weight used for mg/kg dosing (kg)", 30, 3, 150),
      actionButton("go", "Run", class = "btn-primary"),
      hr(),
      p(tags$small(paste("Every default reproduces a row of",
                         "pws_reference_output.txt.  Educational and",
                         "research use only — not for clinical decisions.")))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient & environment",
                 h4("what the four environment knobs actually do"),
                 tableOutput("envtab"),
                 plotOutput("envplot", height = 340),
                 htmlOutput("envnote")),
        tabPanel("2 · The lesion",
                 h4("PC1/3 branches — Eq. A and Eq. B, evaluated live"),
                 tableOutput("branchtab"),
                 plotOutput("branchplot", height = 360),
                 htmlOutput("branchnote")),
        tabPanel("3 · Drug PK",
                 plotOutput("pkplot", height = 500),
                 htmlOutput("pknote")),
        tabPanel("4 · Satiety integrator",
                 h4("the five arms and the harmonic mean"),
                 plotOutput("armplot", height = 330),
                 plotOutput("satplot", height = 260),
                 htmlOutput("satnote")),
        tabPanel("5 · Bistability",
                 h4("the phase line of dSEEK/dt, and its fixed points"),
                 plotOutput("phaseplot", height = 340),
                 tableOutput("fptab"),
                 plotOutput("cuspplot", height = 300),
                 htmlOutput("bistnote")),
        tabPanel("6 · Body composition",
                 plotOutput("compplot", height = 340),
                 plotOutput("energyplot", height = 300),
                 htmlOutput("compnote")),
        tabPanel("7 · Growth & somatotropic",
                 plotOutput("growplot", height = 460),
                 htmlOutput("grownote")),
        tabPanel("8 · Airway (two clocks)",
                 h4("the transient the model does NOT fit"),
                 plotOutput("airplot", height = 460),
                 tableOutput("airtab"),
                 htmlOutput("airnote")),
        tabPanel("9 · Endpoints",
                 plotOutput("endplot", height = 520),
                 tableOutput("endtab")),
        tabPanel("10 · Scenario comparison",
                 h4("biomarker versus endpoint, across the whole panel"),
                 actionButton("runpanel", "Run the panel (8-13 weeks at age 12)",
                              class = "btn-warning"),
                 tableOutput("paneltab"),
                 plotOutput("orthoplot", height = 380),
                 htmlOutput("panelnote"))
      )
    )
  )
)

## =====================================================================
##  SERVER
## =====================================================================
server <- function(input, output, session) {

  envpar <- reactive({
    if (isTRUE(input$envman)) {
      list(AVAIL = input$avail, CUE = input$cue, TITR = input$titr,
           BWTGTR = input$bwtgtr)
    } else ENVS[[input$env]]
  })

  parlist <- reactive({
    c(list(LES = as.numeric(input$les), FSUB = as.numeric(input$subtype),
           DPC13 = input$dpc13, FOXTN = input$foxtn, FOSA = input$fosa),
      envpar())
  })

  sim <- eventReactive(input$go, {
    e1 <- make_ev(input$drug1, input$d1start, input$tend, input$wt)
    e2 <- make_ev(input$drug2, input$d2start, input$tend, input$wt)
    evs <- Filter(Negate(is.null), list(e1, e2))
    m <- mod %>% param(parlist())
    if (length(evs)) m <- m %>% ev(do.call(c, evs))
    m %>% mrgsim(end = YEARS(input$tend), delta = 7, hmax = 0.125) %>%
      as_tibble()
  }, ignoreNULL = FALSE)

  ## ---- tab 1 ---------------------------------------------------------
  output$envtab <- renderTable({
    p <- envpar()
    data.frame(
      knob = c("AVAIL", "CUE", "TITR", "BWTGTR"),
      value = c(p$AVAIL, p$CUE, p$TITR,
                if (is.null(p$BWTGTR)) 1.15 else p$BWTGTR),
      meaning = c("portions offered, as a multiple of the age-normative requirement",
                  "external food-cue and access exposure driving food-seeking",
                  "0 = portions ignore the growth chart; 1 = titrated to a weight target",
                  "the weight-for-height the family or clinic is actually holding"))
  }, digits = 3)

  output$envplot <- renderPlot({
    d <- sim()
    d %>% select(time, EI_OUT, TEE_OUT, EIFRAC, TEEFRAC) %>%
      mutate(age = time / 365.25) %>%
      pivot_longer(c(EI_OUT, TEE_OUT)) %>%
      ggplot(aes(age, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "age (y)", y = "kcal/d", colour = NULL,
           title = "energy intake versus expenditure",
           subtitle = paste("where these two curves cross is where weight",
                            "stops changing")) + THEME
  })

  output$envnote <- renderUI(HTML(paste0(
    "<p><b>Why intake is anchored to what is OFFERED, not to what is spent.</b> ",
    "A PWS child in Miller phase 2a is fed age-normative portions by parents ",
    "who do not yet restrict.  Because expenditure is low — low lean mass ",
    "times low activity, about 0.72 of predicted — a <i>normal</i> intake is ",
    "a surplus.  That is why weight gain <b>precedes</b> hyperphagia in this ",
    "model without any appetite parameter changing.  Set TITR = 0 and watch ",
    "weight-for-age SDS turn positive years before HQ-CT moves.</p>")))

  ## ---- tab 2 ---------------------------------------------------------
  output$branchtab <- renderTable({
    pc <- 1 - input$dpc13 * as.numeric(input$les)
    eps <- c("pro-oxytocin -> OXT" = 6.00, "POMC -> alpha-MSH" = 1.50,
             "pro-GHRH -> GHRH" = 1.00, "proinsulin -> insulin" = 0.15,
             "proghrelin -> acyl-ghrelin" = 0.10)
    data.frame(branch = names(eps), eps = as.numeric(eps),
               `product loss %` = 100 * (1 - (1 + eps) / (1 + eps / pc)),
               `precursor x` = (1 + eps) / (pc + eps),
               `precursor:product` = ((1 + eps) / (pc + eps)) /
                 ((1 + eps) / (1 + eps / pc)),
               check.names = FALSE)
  }, digits = 4)

  output$branchplot <- renderPlot({
    pcs <- seq(1, 0.15, by = -0.01)
    eps <- c(`pro-oxytocin (6.0)` = 6.0, `POMC (1.5)` = 1.5,
             `pro-GHRH (1.0)` = 1.0, `proinsulin (0.15)` = 0.15,
             `proghrelin (0.10)` = 0.10)
    d <- do.call(rbind, lapply(names(eps), function(n) {
      e <- eps[[n]]
      rbind(data.frame(pc = pcs, branch = n, quantity = "product loss",
                       value = 1 - (1 + e) / (1 + e / pcs)),
            data.frame(pc = pcs, branch = n, quantity = "precursor accumulation",
                       value = (1 + e) / (pcs + e) - 1))
    }))
    ggplot(d, aes(pc, value, colour = branch)) + geom_line(linewidth = 0.9) +
      facet_wrap(~quantity, scales = "free_y") +
      geom_vline(xintercept = 1 - input$dpc13 * as.numeric(input$les),
                 linetype = 2) +
      scale_x_reverse() +
      labs(x = "PC1/3 activity", y = NULL, colour = NULL,
           title = "the two panels run in OPPOSITE directions",
           subtitle = paste("the branches that lose product are exactly the",
                            "ones that do NOT accumulate precursor")) + THEME
  })

  output$branchnote <- renderUI(HTML(paste0(
    "<p><b>An exact identity, and the measurement theory it forces.</b><br>",
    "precursor + product = S &nbsp;&rarr;&nbsp; measures <i>synthesis</i>, ",
    "exactly blind to PC1/3.<br>",
    "precursor / product = 1/PC13 &nbsp;&rarr;&nbsp; measures the ",
    "<i>convertase</i>, identical for all five branches, blind to eps.<br>",
    "product alone &nbsp;&rarr;&nbsp; the only quantity carrying branch ",
    "information, and the hardest to assay in a living patient.</p>",
    "<p>The last column of the table above is the same number in every row, ",
    "by algebra.  A trial using a cross-reacting oxytocin assay for target ",
    "engagement is therefore measuring a quantity the model says ",
    "<b>cannot respond</b>.</p>")))

  ## ---- tab 3 ---------------------------------------------------------
  output$pkplot <- renderPlot({
    d <- sim() %>% mutate(age = time / 365.25) %>%
      select(age, C_GH, C_CB, C_DZ, C_SM, C_SG, C_OC, C_LV, C_TS, C_MF) %>%
      pivot_longer(-age) %>% filter(value > 1e-9)
    if (!nrow(d)) return(ggplot() + labs(title = "no drug given") + THEME)
    ggplot(d, aes(age, value)) + geom_line(linewidth = 0.7, colour = "#00695c") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "age (y)", y = "concentration",
           title = "apparent day-scale exposures",
           subtitle = paste("calibrated to reproduce the correct 24-h AUC,",
                            "not the true Cmax")) + THEME
  })

  output$pknote <- renderUI(HTML(paste0(
    "<p>The model runs on a <b>day scale</b>, so every drug is represented by ",
    "an apparent exposure whose 24-h AUC is right.  For somatropin this ",
    "matters: FPOTGH = 0.63 is the pulsatility correction, because pulsatile ",
    "pituitary GH is more IGF-1-efficient per unit AUC than a flat exposure. ",
    "Without it a flattened surrogate over-delivers IGF-1 by about 3 SDS.</p>")))

  ## ---- tab 4 ---------------------------------------------------------
  output$armplot <- renderPlot({
    sim() %>% mutate(age = time / 365.25) %>%
      select(age, `x1 oxytocin relay (w .44)` = ARM_OXT,
             `x2 vagal (w .16)` = ARM_VAGAL, `x3 PYY (w .13)` = ARM_PYY,
             `x4 GLP-1 (w .14)` = ARM_GLP1,
             `x5 leptin/insulin (w .13)` = ARM_LEPINS) %>%
      pivot_longer(-age) %>%
      ggplot(aes(age, value, colour = name)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 1, linetype = 2) +
      labs(x = "age (y)", y = "arm value (control = 1)", colour = NULL,
           title = "the five satiety arms") + THEME
  })

  output$satplot <- renderPlot({
    sim() %>% mutate(age = time / 365.25) %>%
      select(age, SAT_OUT, RELAY_INT) %>% pivot_longer(-age) %>%
      ggplot(aes(age, value, colour = name)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 1, linetype = 2) +
      labs(x = "age (y)", y = NULL, colour = NULL,
           title = "the harmonic mean, and relay integrity",
           subtitle = paste("relay integrity also gates leptin's",
                            "anorexigenic signal")) + THEME
  })

  output$satnote <- renderUI(HTML(paste0(
    "<p><b>Why an intact arm cannot rescue a broken one.</b> ",
    "SAT = 1/&Sigma;(w<sub>i</sub>/x<sub>i</sub>), and ",
    "dSAT/dx<sub>i</sub> = SAT&sup2;&middot;w<sub>i</sub>/x<sub>i</sub>&sup2; — ",
    "sensitivity grows <i>quadratically</i> as an arm gets worse, so a harmonic ",
    "mean concentrates all available leverage on its weakest term.  Normalizing ",
    "all four non-oxytocin arms moves SAT from 0.46 to only 0.49 of a ",
    "possible 1.00.</p>",
    "<p><b>And why MC4R agonism is bounded.</b> alpha-MSH is not its own arm: ",
    "it is a saturating input gain on the oxytocin arm, because MC4R satiety ",
    "signalling is relayed through PVN oxytocin neurons.  relay(&infin;) = ",
    "KREL+1 = 1.20, and PWS's block is <i>below</i> MC4R.  Select ",
    "setmelanotide and then carbetocin and compare x1: the same receptor-level ",
    "logic that makes setmelanotide work in POMC and LEPR deficiency makes it ",
    "fail here.</p>")))

  ## ---- tab 5 ---------------------------------------------------------
  seek_rhs <- function(S, DRVe, SAT, G, kson = 0.120, kshalf = 0.450,
                       nself = 4, ksoff = 2.20) {
    sh <- S^nself
    kson * DRVe * (1 - S) + G * sh / (kshalf^nself + sh) * (1 - S) -
      ksoff * S * SAT
  }
  fps <- function(DRVe, SAT, G, n = 4001) {
    S <- seq(0, 1, length.out = n); f <- seek_rhs(S, DRVe, SAT, G)
    idx <- which(diff(sign(f)) != 0)
    out <- vapply(idx, function(i) {
      uniroot(function(s) seek_rhs(s, DRVe, SAT, G),
              c(S[i], S[i + 1]))$root }, numeric(1))
    data.frame(S = out,
               stability = ifelse(seek_rhs(out + 1e-4, DRVe, SAT, G) < 0,
                                  "stable", "SADDLE"))
  }
  last_state <- reactive({
    d <- sim(); tail(d, 1)
  })

  output$phaseplot <- renderPlot({
    s <- last_state()
    S <- seq(0, 1, length.out = 601)
    d <- data.frame(S = S, dS = seek_rhs(S, s$DRIVE_EFF, s$SAT_OUT, s$SELF_GAIN))
    ggplot(d, aes(S, dS)) + geom_line(linewidth = 0.9, colour = "#c62828") +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = s$SEEK, linetype = 2) +
      labs(x = "SEEK", y = "dSEEK/dt",
           title = "the phase line at the end of the simulation",
           subtitle = sprintf(paste("DRVe = %.3f, SAT = %.3f, self-gain = %.3f;",
                                    "dashed line = where the patient is"),
                              s$DRIVE_EFF, s$SAT_OUT, s$SELF_GAIN)) + THEME
  })

  output$fptab <- renderTable({
    s <- last_state()
    f <- fps(s$DRIVE_EFF, s$SAT_OUT, s$SELF_GAIN)
    f$`patient is here` <- ifelse(abs(f$S - s$SEEK) < 0.03, "<<<", "")
    f
  }, digits = 4)

  output$cuspplot <- renderPlot({
    s <- last_state()
    cues <- seq(0.05, 1.2, by = 0.01)
    d <- do.call(rbind, lapply(cues, function(c) {
      f <- fps(s$DRIVE_INT * c, s$SAT_OUT, s$SELF_GAIN)
      data.frame(cue = c, f)
    }))
    ggplot(d, aes(cue, S, colour = stability)) + geom_point(size = 0.7) +
      labs(x = "food-cue and access exposure (CUE)", y = "fixed points of SEEK",
           colour = NULL,
           title = "the bifurcation diagram",
           subtitle = paste("at high CUE the LOW branch is annihilated —",
                            "no dose can hold a state that does not exist")) +
      THEME
  })

  output$bistnote <- renderUI(HTML(paste0(
    "<p><b>Food security does not lower the drive; it restores the existence ",
    "of a low state.</b> At the model's PWS inputs the saddle-node sits at ",
    "CUE* &asymp; 0.58: above it the vector field has only the food-seeking ",
    "state, so pharmacology has nothing to stabilize.  This is why ",
    "DESTINY-PWS and CARE-PWS both required a stable food-secure environment ",
    "to enrol — and why their effect sizes must not be extrapolated to ",
    "families who do not have one.</p>",
    "<p><b>A trial-design consequence.</b> Because the state is bistable, the ",
    "HQ-CT response should be <b>bimodal</b> (latched versus not) rather than ",
    "shifted, so a mean difference of one or two points may be the wrong ",
    "summary statistic.  And leaving a latched state needs an oxytocin-arm ",
    "gain of about 2.4&times;, where carbetocin's own analytic ceiling is ",
    "1.49&times; — the gap is efficacy, not potency.</p>")))

  ## ---- tab 6 ---------------------------------------------------------
  output$compplot <- renderPlot({
    sim() %>% mutate(age = time / 365.25) %>%
      select(age, `fat mass` = FM, `fat-free mass` = LFM,
             `lean target` = LEAN_TGT, `weight target` = WT_TGT) %>%
      pivot_longer(-age) %>%
      ggplot(aes(age, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "age (y)", y = "kg", colour = NULL,
           title = "lean mass follows a developmental target; FAT IS THE BUFFER",
           subtitle = paste("energy is exactly conserved, so PWS's lean deficit",
                            "is derived and not imposed")) + THEME
  })

  output$energyplot <- renderPlot({
    sim() %>% mutate(age = time / 365.25) %>%
      select(age, PBF, TEEFRAC, EIFRAC, MEALLOAD) %>% pivot_longer(-age) %>%
      ggplot(aes(age, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "age (y)", y = NULL, colour = NULL) +
      THEME + theme(legend.position = "none")
  })

  output$compnote <- renderUI(HTML(paste0(
    "<p>With TITR = 1 the caregiver titrates portions to a weight-for-height ",
    "target, so <b>body weight is pinned and composition is decided by the ",
    "lean-mass target</b>.  That is the mechanism by which growth hormone ",
    "lowers %fat while barely moving BMI — the finding the GH trials in PWS ",
    "keep reporting.  Turn somatropin on and off and watch %fat move while ",
    "the weight target does not.</p>")))

  ## ---- tab 7 ---------------------------------------------------------
  output$growplot <- renderPlot({
    sim() %>% mutate(age = time / 365.25) %>%
      select(age, `height (cm)` = HT, `height SDS` = HTSDS,
             `IGF-1 SDS` = IGFSDS, `GH signal (rel)` = GHREL_OUT,
             `growth velocity (cm/y)` = GV_OUT, `bone age (y)` = BA) %>%
      pivot_longer(-age) %>%
      ggplot(aes(age, value)) + geom_line(linewidth = 0.9, colour = "#1565c0") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "age (y)", y = NULL,
           title = "growth and the somatotropic axis") + THEME
  })

  output$grownote <- renderUI(HTML(paste0(
    "<p>Growth velocity has a large <b>GH-independent floor</b> (GVFL = 0.72), ",
    "which is why untreated PWS loses about 20 cm of adult height and not 60. ",
    "The pubertal component is gated by age-normalized sex steroid, so PWS's ",
    "hypogonadism blunts the spurt specifically.  Epiphyses also close later ",
    "in PWS, which is why late growth hormone still buys height and why ",
    "sex-steroid replacement shortens that window as it works.</p>")))

  ## ---- tab 8 ---------------------------------------------------------
  output$airplot <- renderPlot({
    sim() %>% mutate(age = time / 365.25) %>%
      select(age, `AHI (events/h)` = AHI, `lymphoid volume` = LYMPH,
             `respiratory muscle` = RMS, `%fat` = PBF) %>%
      pivot_longer(-age) %>%
      ggplot(aes(age, value)) + geom_line(linewidth = 0.9, colour = "#00838f") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "age (y)", y = NULL,
           title = "the two clocks, and what they do to the airway",
           subtitle = paste("tau_lymphoid ~ 20 d  <  tau_muscle 75 d  <",
                            "tau_fat (months)")) + THEME
  })

  output$airtab <- renderTable({
    d <- sim() %>% mutate(age = time / 365.25)
    st <- input$d1start
    d %>% filter(age >= st - 0.02, age <= st + 2.0) %>%
      mutate(week = round((age - st) * 52.18)) %>%
      filter(week %in% c(0, 2, 4, 6, 8, 12, 26, 52, 104)) %>%
      group_by(week) %>% slice(1) %>% ungroup() %>%
      select(week, LYMPH, RMS, PBF, AHI)
  }, digits = 3)

  output$airnote <- renderUI(HTML(paste0(
    "<p><b>The window is emergent, not fitted.</b> Growth hormone reaches the ",
    "upper airway twice: fast, through IGF-1 driving lymphoid and ",
    "adenotonsillar hypertrophy (tau ~20 d, with a 100-d involution), and ",
    "slowly, through lean mass raising respiratory muscle strength and fat ",
    "mass falling.  Because the fast clock leads, AHI rises for 4-8 weeks ",
    "and then falls below baseline.  Airway vulnerability (FOSA) multiplies ",
    "only the obstructive terms, so it amplifies the <i>window</i> and not the ",
    "plateau: 11.2 &rarr; 22.1 &rarr; 4.2 events/h at FOSA = 1.9.</p>",
    "<p>This is the mechanistic content of <i>polysomnography before, and 6-8 ",
    "weeks after, starting growth hormone</i>, and of not starting during a ",
    "respiratory infection.  Set FOSA to 1.9 and re-run.</p>")))

  ## ---- tab 9 ---------------------------------------------------------
  output$endplot <- renderPlot({
    sim() %>% mutate(age = time / 365.25) %>%
      select(age, `HQ-CT` = HQ, `%fat` = PBF, BMI, `height SDS` = HTSDS,
             `IGF-1 SDS` = IGFSDS, `AHI` = AHI, `HbA1c (%)` = HBA1C,
             `BMD (g/cm2)` = BMD, `Cobb (deg)` = COBB,
             `behaviour` = BEH, `SEEK` = SEEK,
             `adiponectin` = ADPN) %>%
      pivot_longer(-age) %>%
      ggplot(aes(age, value)) + geom_line(linewidth = 0.9, colour = "#ad1457") +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      labs(x = "age (y)", y = NULL, title = "clinical endpoints") + THEME
  })

  output$endtab <- renderTable({
    tail(sim(), 1) %>%
      transmute(age = time / 365.25, BMI, `%fat` = PBF, `HQ-CT` = HQ,
                SEEK, AHI, `height` = HT, `height SDS` = HTSDS,
                `IGF-1 SDS` = IGFSDS, HbA1c = HBA1C, BMD, Cobb = COBB)
  }, digits = 3)

  ## ---- tab 10 --------------------------------------------------------
  panel <- eventReactive(input$runpanel, {
    base <- c(list(LES = 1, FSUB = as.numeric(input$subtype),
                   DPC13 = input$dpc13, FOXTN = input$foxtn),
              ENVS[["PWS, standard management"]])
    gh <- ev(time = YEARS(1), amt = 0.035 * input$wt, cmt = "AGHD",
             ii = 1, addl = 8760)
    arms <- list(
      "placebo"                  = NULL,
      "carbetocin 3.2 mg TID"    = ev(time = YEARS(12), amt = 3.2, cmt = "ACBD",
                                      ii = 1/3, addl = 167),
      "carbetocin 9.6 mg TID"    = ev(time = YEARS(12), amt = 9.6, cmt = "ACBD",
                                      ii = 1/3, addl = 167),
      "DCCR 5.1 mg/kg"           = ev(time = YEARS(12), amt = 5.1 * input$wt,
                                      cmt = "ADZD", ii = 1, addl = 90),
      "setmelanotide 3 mg"       = ev(time = YEARS(12), amt = 3.0, cmt = "ASMD",
                                      ii = 1, addl = 83),
      "livoletide 60 ug/kg"      = ev(time = YEARS(12), amt = 0.060 * input$wt,
                                      cmt = "ALVD", ii = 1, addl = 83),
      "octreotide LAR 30 mg"     = ev(time = YEARS(12), amt = 30, cmt = "AOCD",
                                      ii = 28, addl = 2),
      "semaglutide 2.4 mg QW"    = ev(time = YEARS(12), amt = 2.4, cmt = "ASGD",
                                      ii = 7, addl = 12),
      "metformin 1500 mg"        = ev(time = YEARS(12), amt = 1500, cmt = "AMFD",
                                      ii = 1, addl = 90))
    res <- lapply(names(arms), function(nm) {
      evs <- if (is.null(arms[[nm]])) gh else c(gh, arms[[nm]])
      out <- mod %>% param(base) %>% ev(evs) %>%
        mrgsim(end = YEARS(12) + 91, delta = 7, hmax = 0.125) %>% as_tibble()
      b <- out %>% filter(abs(time - YEARS(12)) == min(abs(time - YEARS(12)))) %>%
        slice(1)
      e <- tail(out, 1)
      tibble(arm = nm, dHQ = e$HQ - b$HQ,
             dAG = 100 * (e$AG / b$AG - 1),
             dAGUAG = 100 * (e$AGUAGEFF / b$AGUAGEFF - 1),
             dIGF1 = 100 * (e$IGF1 / b$IGF1 - 1),
             dPBF = e$PBF - b$PBF, dA1C = e$HBA1C - b$HBA1C,
             SEEKend = e$SEEK)
    })
    bind_rows(res) %>%
      mutate(dHQ_vs_placebo = dHQ - dHQ[arm == "placebo"])
  })

  output$paneltab <- renderTable({ panel() }, digits = 3)

  output$orthoplot <- renderPlot({
    d <- panel() %>% filter(arm != "placebo")
    ggplot(d, aes(dAG, dHQ_vs_placebo, label = arm)) +
      geom_point(size = 3, colour = "#c62828") +
      geom_text(hjust = -0.08, size = 3.4) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_vline(xintercept = 0, linetype = 2) +
      expand_limits(x = c(-110, 60)) +
      labs(x = "change in acyl-ghrelin (%)",
           y = "placebo-corrected change in HQ-CT",
           title = "the biomarker and the endpoint are nearly ORTHOGONAL",
           subtitle = paste("the agents that move ghrelin most move HQ-CT",
                            "least; octreotide is the extreme case")) + THEME
  })

  output$panelnote <- renderUI(HTML(paste0(
    "<p>Three of these arms were designed around ghrelin and all three fail ",
    "the endpoint while succeeding on the biomarker.  In the model that is one ",
    "number: the GHS-R1a occupancy is <b>already saturated</b> at PWS ",
    "concentrations (K = 300 pg/mL, PWS &asymp; 870), so a two-fold ",
    "hyperghrelinaemia buys only +27% of receptor arm and abolishing it ",
    "entirely can return at most 2.7% of the drive.</p>",
    "<p><b>Honest labelling.</b> The ghrelin-arm weight, the OXTR Emax and the ",
    "K-ATP Emax are <i>calibrated</i> to reported trial outcomes.  What is ",
    "<i>derived</i> is the structure that makes those outcomes mutually ",
    "consistent: the Miller phase sequence, the airway window, the inverted ",
    "carbetocin dose-response, this orthogonality, the bimodality prediction, ",
    "and the 2.4&times; oxytocin-arm gain a next-generation agonist would have ",
    "to reach.</p>")))
}

shinyApp(ui, server)
