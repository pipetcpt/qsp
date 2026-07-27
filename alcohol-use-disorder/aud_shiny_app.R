###############################################################################
##  Alcohol Use Disorder (AUD) — QSP interactive dashboard
##  ---------------------------------------------------------------------------
##  aud_shiny_app.R      companion to aud_mrgsolve_model.R / aud_qsp_model.dot
##
##  The point of this app is to let you drive the CLOSED LOOP by hand:
##  you set the environment, the host and the regimen, and the app shows you
##  the drinking the model then generates — not the other way round.
##
##  Nine tabs:
##    1  Patient & environment   who they are, where they drink, what they took
##    2  Ethanol PK              BAC, acetaldehyde, first pass, metabolic tolerance
##    3  Reinforcement arm       beta-endorphin -> MOR -> dopamine -> sensitisation
##    4  Negative-affect arm     GABA/NMDA/GLT-1/dynorphin/CRF -> craving
##    5  Drinking & endpoints    drinks/day, %HDD, PDA, DPDD, WHO risk level
##    6  Withdrawal & detox      CIWA-Ar, excitation index, seizure risk, benzos
##    7  Biomarkers              PEth, CDT, GGT, MCV, AST/ALT
##    8  Organ systems           liver, thiamine, blood pressure, cardiomyopathy
##    9  Scenario comparison     stack any number of saved runs side by side
##
##  Run:  Rscript -e 'shiny::runApp("aud_shiny_app.R")'
###############################################################################

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

## --------------------------------------------------------------------------
## Load the model definition straight out of aud_mrgsolve_model.R so that the
## app and the batch script can never drift apart.
## --------------------------------------------------------------------------
MODEL_FILE <- "aud_mrgsolve_model.R"
src  <- readLines(MODEL_FILE, warn = FALSE)
i1   <- grep("^code <- '", src)[1]
i2   <- grep("^mod <- mcode_cache", src)[1]
eval(parse(text = paste(src[i1:(i2 - 1)], collapse = "\n")))
mod  <- mcode_cache("aud_shiny", code, soloc = tempdir())

DAY <- 24

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"))

PAL <- c("#2b6cb0", "#c05621", "#2f855a", "#805ad5", "#b83280",
         "#276749", "#975a16", "#822727")

## --------------------------------------------------------------------------
## Endpoint mapping (identical to the batch script)
## --------------------------------------------------------------------------
Z0 <- -1.35; Z1 <- 1.60; Z2 <- 0.60; Z3 <- 0.12; NBK <- 2.4

ztnb_ge <- function(thresh, mu, size) {
  p0 <- dnbinom(0, mu = mu, size = size)
  pl <- pnbinom(thresh - 1, mu = mu, size = size)
  max(0, (1 - pl)/(1 - p0))
}

endpoints <- function(dsum, ref, abstg, hdd_thr = 5) {
  mu  <- mean(dsum$DPD); neg <- mean(dsum$NEGAFF); crv <- mean(dsum$CRAVE)
  dneg <- if (is.null(ref)) 0 else (ref$neg - neg)/max(ref$neg, 1e-6)
  dcrv <- if (is.null(ref)) 0 else (ref$crv - crv)/max(ref$crv, 1e-6)
  pi_ab <- min(max(plogis(Z0 + Z1*abstg + Z2*dneg + Z3*dcrv), 0.001), 0.999)
  dpdd  <- mu/(1 - pi_ab)
  phdd  <- 100*(1 - pi_ab)*ztnb_ge(hdd_thr, dpdd, NBK)
  list(mu = mu, neg = neg, crv = crv, PDA = 100*pi_ab, DPDD = dpdd,
       HDD = phdd, TAC = mu*14,
       WHO = cut(mu*14, c(-1, 1, 40, 60, 100, 1e6),
                 labels = c("abstinent", "low", "medium", "high", "very high")))
}

daily <- function(d) {
  d %>% mutate(day = floor(time/DAY)) %>% group_by(day) %>%
    summarise(GPD = mean(GPD), DPD = mean(DPD_ACT), BACmax = max(BACO),
              CIWAmax = max(CIWA), CRAVE = mean(CRAVE), NEGAFF = mean(NEGAFF),
              PETH = mean(PETH), GGT = mean(GGT), CDT = mean(CDT),
              MCV = mean(MCV), ALT = mean(ALT), AST = mean(AST),
              SBP = mean(SBP), STEAT = mean(STEAT), THIA = mean(THIA),
              FIB = mean(FIB), CMYO = mean(CMYO), .groups = "drop")
}

## --------------------------------------------------------------------------
## Pre-computed disease-history states (built once at start-up)
## --------------------------------------------------------------------------
message("building AUD history states (this takes ~1 min)...")
.burn <- function(env, vuln) {
  o <- mod %>% param(ENVDRIVE = env, VULN = vuln) %>%
    mrgsim(end = (3*365 + 21/24)*DAY, delta = 6, hmax = 0.5,
           maxsteps = 500000) %>% as_tibble()
  s <- as.list(o[nrow(o), names(init(mod))])
  s$ETOHCUM <- 0; s$BZDCUM <- 0
  s
}
HIST <- list(
  `healthy (never a problem drinker)` = as.list(init(mod)),
  `moderate drinker (~2 drinks/day)`  = .burn(1.6, 1.00),
  `AUD, trial-entry severity`         = .burn(3.4, 1.15),
  `AUD, severe`                       = .burn(4.6, 1.35))
message("done.")

DRUGS <- c("none", "naltrexone 50 mg PO", "XR-naltrexone 380 mg IM q4wk",
           "acamprosate 666 mg tid", "naltrexone + acamprosate",
           "disulfiram 250 mg (supervised)", "disulfiram 250 mg (unsupervised)",
           "topiramate 300 mg/day", "gabapentin 1800 mg/day",
           "nalmefene 18 mg as-needed", "baclofen 180 mg/day",
           "ondansetron 4 ug/kg bid", "semaglutide 1.0 mg weekly")

build_ev <- function(drug, dur_d, bw) {
  n <- function(ii, d = dur_d) max(0, floor(d*DAY/ii) - 1)
  switch(drug,
    "none" = NULL,
    "naltrexone 50 mg PO" =
      ev(time = 0, amt = 50, cmt = "NTX_GUT", ii = 24, addl = n(24)),
    "XR-naltrexone 380 mg IM q4wk" =
      ev(time = 0, amt = 380, cmt = "NTX_DEP", ii = 28*DAY, addl = n(28*DAY)),
    "acamprosate 666 mg tid" =
      ev(time = 0, amt = 666, cmt = "ACP_GUT", ii = 8, addl = n(8)),
    "naltrexone + acamprosate" =
      ev(time = 0, amt = 50, cmt = "NTX_GUT", ii = 24, addl = n(24)) +
      ev(time = 0, amt = 666, cmt = "ACP_GUT", ii = 8, addl = n(8)),
    "disulfiram 250 mg (supervised)" =
      ev(time = 0, amt = 250, cmt = "DSF_GUT", ii = 24, addl = n(24)),
    "disulfiram 250 mg (unsupervised)" =
      ev(time = 0, amt = 250, cmt = "DSF_GUT", ii = 96, addl = n(96)),
    "topiramate 300 mg/day" =
      ev(time = 0, amt = 50, cmt = "TPM_GUT", ii = 24, addl = 13) +
      ev(time = 14*DAY, amt = 150, cmt = "TPM_GUT", ii = 24, addl = 13) +
      ev(time = 28*DAY, amt = 300, cmt = "TPM_GUT", ii = 24, addl = n(24, dur_d - 28)),
    "gabapentin 1800 mg/day" =
      ev(time = 0, amt = 600, cmt = "GBP_GUT", ii = 8, addl = n(8)),
    "nalmefene 18 mg as-needed" =
      ev(time = 14, amt = 18, cmt = "NMF_GUT", ii = 24, addl = n(24)),
    "baclofen 180 mg/day" =
      ev(time = 0, amt = 60, cmt = "BACL_GUT", ii = 8, addl = n(8)),
    "ondansetron 4 ug/kg bid" =
      ev(time = 0, amt = 0.004*bw, cmt = "OND_GUT", ii = 12, addl = n(12)),
    "semaglutide 1.0 mg weekly" =
      ev(time = 0, amt = 0.25, cmt = "SEM_SC", ii = 168, addl = 3) +
      ev(time = 4*168, amt = 0.5, cmt = "SEM_SC", ii = 168, addl = 3) +
      ev(time = 8*168, amt = 1.0, cmt = "SEM_SC", ii = 168, addl = n(168, dur_d - 56)))
}

drug_par <- function(drug) {
  if (drug == "disulfiram 250 mg (supervised)") list(SUPERV = 1) else list()
}

###############################################################################
## UI
###############################################################################

ui <- fluidPage(
  titlePanel("Alcohol Use Disorder — QSP dashboard (closed-loop drinking model)"),
  tags$p(style = "color:#555;margin-top:-8px;",
         HTML("In this model the <b>dose is an output</b>: you set the host, the
               environment and the regimen, and the simulation generates the
               drinking. 72 ODE compartments · educational use only.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("hist", "Drinking history (starting state)",
                  choices = names(HIST), selected = "AUD, trial-entry severity"),
      sliderInput("bw", "Body weight (kg)", 45, 130, 75, 1),
      radioButtons("sex", "Sex", c("male" = 0, "female" = 1), inline = TRUE),
      sliderInput("fed", "Fed state (0 fasted – 1 fed)", 0, 1, 0.5, 0.05),
      h4("Pharmacogenetics"),
      sliderInput("aldhf", "ALDH2 activity (1 = wild-type, 0.2 = *2 het)",
                  0.02, 1.2, 1.0, 0.02),
      sliderInput("adhf", "ADH1B activity (2.5 = ADH1B*2)", 0.5, 3.0, 1.0, 0.1),
      sliderInput("oprmf", "OPRM1 naltrexone potency (1.6 = 118G carrier)",
                  0.6, 2.2, 1.0, 0.1),
      sliderInput("httf", "5-HTTLPR ondansetron responsiveness",
                  0.5, 2.5, 1.0, 0.1),
      h4("Environment & goal"),
      sliderInput("env", "Environmental drinking pressure", 0, 7, 3.4, 0.1),
      sliderInput("vuln", "Neuro-vulnerability", 0.6, 2.0, 1.15, 0.05),
      sliderInput("abstg", "Abstinence-goal strength (psychosocial)",
                  0, 1, 0.55, 0.05),
      checkboxInput("dry", "Enforce abstinence (inpatient / detox)", FALSE),
      h4("Regimen"),
      selectInput("drug", "Medication", choices = DRUGS,
                  selected = "naltrexone 50 mg PO"),
      checkboxInput("bzd", "Symptom-triggered benzodiazepine protocol", FALSE),
      checkboxInput("thia", "Thiamine replacement", FALSE),
      checkboxInput("praz", "Prazosin (alpha1 blocker)", FALSE),
      sliderInput("dur", "Simulation length (days)", 14, 365, 180, 7),
      actionButton("go", "Simulate", class = "btn-primary"),
      tags$hr(),
      textInput("label", "Label for comparison", "run 1"),
      actionButton("save", "Save to comparison"),
      actionButton("clear", "Clear comparison")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 · Patient & environment",
                 br(), verbatimTextOutput("summary"),
                 plotOutput("p_profile", height = 380)),
        tabPanel("2 · Ethanol PK",
                 br(), plotOutput("p_bac", height = 300),
                 plotOutput("p_metab", height = 300)),
        tabPanel("3 · Reinforcement arm",
                 br(), plotOutput("p_reward", height = 340),
                 plotOutput("p_sens", height = 280)),
        tabPanel("4 · Negative-affect arm",
                 br(), plotOutput("p_allo", height = 340),
                 plotOutput("p_neg", height = 280)),
        tabPanel("5 · Drinking & endpoints",
                 br(), plotOutput("p_drink", height = 300),
                 tableOutput("t_end"),
                 plotOutput("p_who", height = 240)),
        tabPanel("6 · Withdrawal & detox",
                 br(), plotOutput("p_ciwa", height = 340),
                 plotOutput("p_bzd", height = 260)),
        tabPanel("7 · Biomarkers",
                 br(), plotOutput("p_bm", height = 460),
                 tableOutput("t_bm")),
        tabPanel("8 · Organ systems",
                 br(), plotOutput("p_organ", height = 460)),
        tabPanel("9 · Scenario comparison",
                 br(), tableOutput("t_cmp"),
                 plotOutput("p_cmp", height = 380))
      )
    )
  )
)

###############################################################################
## SERVER
###############################################################################

server <- function(input, output, session) {

  store <- reactiveVal(list())

  simulate <- eventReactive(input$go, {
    ini <- HIST[[input$hist]]
    par <- c(list(BW = input$bw, SEXF = as.numeric(input$sex), FED = input$fed,
                  ALDHF = input$aldhf, ADHF = input$adhf, OPRMF = input$oprmf,
                  HTTF = input$httf, ENVDRIVE = input$env, VULN = input$vuln,
                  ABSTG = input$abstg,
                  DRINKSW = if (input$dry) 0 else 1,
                  BZDON = if (input$bzd) 1 else 0,
                  THIADOSE = if (input$thia) 0.02 else 0,
                  PRAZON = if (input$praz) 1 else 0),
             drug_par(input$drug))
    e <- build_ev(input$drug, input$dur, input$bw)
    m <- mod %>% init(ini) %>% param(par)
    out <- if (is.null(e)) {
      m %>% mrgsim(end = input$dur*DAY, delta = 0.25, hmax = 0.5,
                   maxsteps = 500000)
    } else {
      m %>% mrgsim(events = e, end = input$dur*DAY, delta = 0.25, hmax = 0.5,
                   maxsteps = 500000)
    }
    out <- as_tibble(out)
    list(raw = out, day = daily(out), lab = input$label,
         abstg = input$abstg, drug = input$drug)
  })

  ## reference (no-medication) run for the endpoint mapping ------------------
  reference <- reactive({
    s <- simulate()
    ini <- HIST[[input$hist]]
    par <- list(BW = input$bw, SEXF = as.numeric(input$sex), FED = input$fed,
                ALDHF = input$aldhf, ADHF = input$adhf,
                ENVDRIVE = input$env, VULN = input$vuln, ABSTG = input$abstg,
                DRINKSW = if (input$dry) 0 else 1)
    o <- mod %>% init(ini) %>% param(par) %>%
      mrgsim(end = input$dur*DAY, delta = 1, hmax = 0.5, maxsteps = 500000) %>%
      as_tibble() %>% daily()
    w <- o %>% filter(day >= min(30, max(o$day)/3))
    list(neg = mean(w$NEGAFF), crv = mean(w$CRAVE), mu = mean(w$DPD))
  })

  wnd <- function(d) d %>% filter(day >= min(30, max(d$day)/3))

  ## ---- tab 1 --------------------------------------------------------------
  output$summary <- renderPrint({
    s <- simulate(); r <- reference()
    e <- endpoints(wnd(s$day), r, s$abstg)
    e0 <- endpoints(wnd(s$day), NULL, s$abstg)
    last <- tail(s$day, 1)
    cat(sprintf("history          : %s\n", input$hist))
    cat(sprintf("medication       : %s\n", s$drug))
    cat(sprintf("consumption      : %.1f drinks/day  (%.0f g ethanol/day)  WHO risk level: %s\n",
                e$mu, e$TAC, as.character(e$WHO)))
    cat(sprintf("endpoints        : %%HDD %.1f   %%PDA %.1f   DPDD %.2f\n",
                e$HDD, e$PDA, e$DPDD))
    cat(sprintf("versus no medication in the same patient: dTAC %+.0f g/day\n",
                (e$mu - r$mu)*14))
    cat(sprintf("craving (PACS)   : %.2f      negative affect: %.2f\n", e$crv, e$neg))
    cat(sprintf("peak CIWA-Ar     : %.1f      max modelled seizure risk: %.1f %%\n",
                max(s$raw$CIWA), max(s$raw$SEIZP)))
    cat(sprintf("biomarkers (end) : PEth %.0f ng/mL  CDT %.2f %%  GGT %.0f U/L  MCV %.1f fL\n",
                last$PETH, last$CDT, last$GGT, last$MCV))
    cat(sprintf("liver            : AST %.0f / ALT %.0f (ratio %.2f)  steatosis %.0f %%  fibrosis %.2f\n",
                last$AST, last$ALT, last$AST/last$ALT, 100*last$STEAT, last$FIB))
    cat(sprintf("MOR occupancy    : %.1f %%    total benzodiazepine %.0f mg diazepam-eq\n",
                mean(s$raw$MOROCC), max(s$raw$BZDCUM)))
  })

  output$p_profile <- renderPlot({
    s <- simulate()
    d <- s$day %>% select(day, `drinks/day` = DPD, `craving (PACS)` = CRAVE,
                          `negative affect` = NEGAFF, `peak daily CIWA` = CIWAmax) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value)) +
      geom_line(colour = PAL[1], linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "The closed loop, day by day", x = "day", y = NULL) + THEME
  })

  ## ---- tab 2 --------------------------------------------------------------
  output$p_bac <- renderPlot({
    s <- simulate()
    last3 <- s$raw %>% filter(time >= max(time) - 72)
    ggplot(last3, aes(time/24, BACO)) +
      geom_hline(yintercept = 0.08, linetype = 2, colour = "#b83280") +
      annotate("text", x = -Inf, y = 0.083, hjust = -0.1, size = 3.2,
               colour = "#b83280", label = "binge threshold 0.08 g/dL") +
      geom_line(colour = PAL[1], linewidth = 0.7) +
      labs(title = "Blood alcohol concentration — last three simulated days",
           x = "day", y = "BAC (g/dL)") + THEME
  })

  output$p_metab <- renderPlot({
    s <- simulate()
    last3 <- s$raw %>% filter(time >= max(time) - 72) %>%
      select(time, `acetaldehyde (uM)` = CACD, `CYP2E1 (fold)` = CYP2E1,
             `NADH/NAD+ (rel)` = NADH, `lactate (mM)` = LACT) %>%
      pivot_longer(-time)
    ggplot(last3, aes(time/24, value)) +
      geom_line(colour = PAL[2], linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Oxidative metabolism and redox", x = "day", y = NULL) + THEME
  })

  ## ---- tab 3 --------------------------------------------------------------
  output$p_reward <- renderPlot({
    s <- simulate()
    last3 <- s$raw %>% filter(time >= max(time) - 72) %>%
      select(time, `beta-endorphin` = BEND, `MOR activation` = MORACT,
             `MOR occupancy by antagonist (%)` = MOROCC,
             `accumbal dopamine` = DA_NAC) %>%
      pivot_longer(-time)
    ggplot(last3, aes(time/24, value)) +
      geom_line(colour = PAL[3], linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Positive-reinforcement arm (the naltrexone / nalmefene target)",
           x = "day", y = NULL) + THEME
  })

  output$p_sens <- renderPlot({
    s <- simulate()
    d <- s$raw %>% filter(time %% 6 < 0.3) %>%
      select(time, `incentive salience (CUE)` = CUE, `habit strength` = HABIT,
             `reward set-point shift` = ALLO, `prefrontal control` = PFC) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/24, value)) +
      geom_line(colour = PAL[4], linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Incentive sensitisation, habit capture, control erosion",
           x = "day", y = NULL) + THEME
  })

  ## ---- tab 4 --------------------------------------------------------------
  output$p_allo <- renderPlot({
    s <- simulate()
    d <- s$raw %>% filter(time %% 6 < 0.3) %>%
      select(time, `GABA-A subunit index` = GABAA_SUB, `NMDA upregulation` = NMDA_UP,
             `GLT-1 capacity` = GLT1, `extracellular glutamate` = GLU_NAC,
             `dynorphin tone` = DYN, `CeA CRF tone` = CRF_CEA) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/24, value)) +
      geom_line(colour = PAL[5], linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Allostatic neuroadaptation (the acamprosate / gabapentin target)",
           x = "day", y = NULL) + THEME
  })

  output$p_neg <- renderPlot({
    s <- simulate()
    d <- s$raw %>% filter(time %% 6 < 0.3) %>%
      select(time, `negative affect` = NEGAFF, `sleep disruption` = SLEEPD,
             `NPY` = NPY_CEA, `cortisol` = CORT) %>%
      pivot_longer(-time)
    ggplot(d, aes(time/24, value)) +
      geom_line(colour = PAL[6], linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Negative affect, sleep and the HPA axis", x = "day", y = NULL) +
      THEME
  })

  ## ---- tab 5 --------------------------------------------------------------
  output$p_drink <- renderPlot({
    s <- simulate()
    ggplot(s$day, aes(day)) +
      geom_col(aes(y = GPD), fill = PAL[1], alpha = 0.75) +
      geom_hline(yintercept = c(40, 60, 100), linetype = 2,
                 colour = c("#2f855a", "#c05621", "#822727")) +
      labs(title = "Generated consumption and WHO risk-level bands",
           subtitle = "dashed lines: WHO medium (40) / high (60) / very high (100) g per day",
           x = "day", y = "g ethanol per day") + THEME
  })

  output$t_end <- renderTable({
    s <- simulate(); r <- reference()
    e <- endpoints(wnd(s$day), r, s$abstg)
    data.frame(
      endpoint = c("drinks per day (mean)", "g ethanol per day",
                   "% heavy drinking days", "% days abstinent",
                   "drinks per drinking day", "WHO risk level",
                   "delta TAC vs no medication (g/day)"),
      value = c(sprintf("%.2f", e$mu), sprintf("%.0f", e$TAC),
                sprintf("%.1f", e$HDD), sprintf("%.1f", e$PDA),
                sprintf("%.2f", e$DPDD), as.character(e$WHO),
                sprintf("%+.0f", (e$mu - r$mu)*14)))
  })

  output$p_who <- renderPlot({
    s <- simulate()
    d <- s$day %>% mutate(lvl = cut(GPD, c(-1, 1, 40, 60, 100, 1e6),
      labels = c("abstinent", "low", "medium", "high", "very high")))
    ggplot(d, aes(day, 1, fill = lvl)) + geom_tile() +
      scale_fill_manual(values = c("#2f855a", "#68a691", "#e9c46a",
                                   "#e76f51", "#9b2226"), drop = FALSE) +
      labs(title = "WHO risk drinking level over time", x = "day",
           y = NULL, fill = NULL) +
      THEME + theme(axis.text.y = element_blank())
  })

  ## ---- tab 6 --------------------------------------------------------------
  output$p_ciwa <- renderPlot({
    s <- simulate()
    d <- s$raw %>% filter(time <= 21*24) %>%
      select(time, `CIWA-Ar` = CIWA, `excitation index` = EXC,
             `GABA tone` = GABA_TONE, `seizure risk (%)` = SEIZP) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) +
      geom_line(colour = PAL[8], linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Withdrawal (first 21 days)", x = "hours", y = NULL) + THEME
  })

  output$p_bzd <- renderPlot({
    s <- simulate()
    d <- s$raw %>% filter(time <= 21*24) %>%
      select(time, `benzodiazepine (mg/L diazepam-eq)` = CBZD,
             `cumulative dose (mg)` = BZDCUM, `kindling index` = KINDLE) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) +
      geom_line(colour = PAL[7], linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Symptom-triggered benzodiazepine exposure and kindling",
           x = "hours", y = NULL) + THEME
  })

  ## ---- tab 7 --------------------------------------------------------------
  output$p_bm <- renderPlot({
    s <- simulate()
    d <- s$day %>%
      select(day, `PEth (ng/mL)` = PETH, `CDT (%)` = CDT, `GGT (U/L)` = GGT,
             `MCV (fL)` = MCV, `AST (U/L)` = AST, `ALT (U/L)` = ALT) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value)) +
      geom_line(colour = PAL[1], linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Objective consumption biomarkers", x = "day", y = NULL) + THEME
  })

  output$t_bm <- renderTable({
    s <- simulate(); a <- s$day[1, ]; b <- tail(s$day, 1)
    data.frame(marker = c("PEth (ng/mL)", "CDT (%)", "GGT (U/L)", "MCV (fL)",
                          "AST (U/L)", "ALT (U/L)", "AST/ALT"),
               start = sprintf("%.2f", c(a$PETH, a$CDT, a$GGT, a$MCV, a$AST,
                                         a$ALT, a$AST/a$ALT)),
               end   = sprintf("%.2f", c(b$PETH, b$CDT, b$GGT, b$MCV, b$AST,
                                         b$ALT, b$AST/b$ALT)),
               `reference / cut-off` = c("< 20 abstinent; > 200 heavy",
                                         "< 1.7", "< 55", "80-100", "< 40",
                                         "< 45", "> 2 suggests alcohol"),
               check.names = FALSE)
  })

  ## ---- tab 8 --------------------------------------------------------------
  output$p_organ <- renderPlot({
    s <- simulate()
    d <- s$day %>%
      select(day, `steatosis (fraction)` = STEAT, `fibrosis stage` = FIB,
             `thiamine (relative)` = THIA, `systolic BP (mmHg)` = SBP,
             `cardiomyopathy index` = CMYO) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value)) +
      geom_line(colour = PAL[6], linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "End-organ trajectories", x = "day", y = NULL) + THEME
  })

  ## ---- tab 9 --------------------------------------------------------------
  observeEvent(input$save, {
    s <- simulate(); r <- reference()
    e <- endpoints(wnd(s$day), r, s$abstg)
    st <- store()
    st[[length(st) + 1]] <- list(
      lab = s$lab, drug = s$drug, day = s$day,
      row = data.frame(label = s$lab, medication = s$drug,
                       `g/day` = round(e$TAC), DPDD = round(e$DPDD, 2),
                       `%HDD` = round(e$HDD, 1), `%PDA` = round(e$PDA, 1),
                       craving = round(e$crv, 2),
                       `peak CIWA` = round(max(s$raw$CIWA), 1),
                       `PEth end` = round(tail(s$day$PETH, 1)),
                       check.names = FALSE))
    store(st)
  })
  observeEvent(input$clear, store(list()))

  output$t_cmp <- renderTable({
    st <- store(); if (!length(st)) return(NULL)
    bind_rows(lapply(st, `[[`, "row"))
  })

  output$p_cmp <- renderPlot({
    st <- store(); if (!length(st)) return(NULL)
    d <- bind_rows(lapply(st, function(x)
      x$day %>% transmute(day, GPD, CRAVE, PETH, label = x$lab))) %>%
      pivot_longer(c(GPD, CRAVE, PETH))
    ggplot(d, aes(day, value, colour = label)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y",
                 labeller = as_labeller(c(GPD = "g ethanol/day",
                                          CRAVE = "craving (PACS)",
                                          PETH = "PEth (ng/mL)"))) +
      scale_colour_manual(values = rep(PAL, 4)) +
      labs(title = "Saved scenarios", x = "day", y = NULL, colour = NULL) + THEME
  })
}

shinyApp(ui, server)
