## =============================================================================
##  Carbon Monoxide Poisoning -- interactive QSP dashboard
##  Front-end for co_mrgsolve_model.R (45 ODEs)
## =============================================================================
##
##  The dashboard is organised around the model's central claim: there are TWO
##  CO occupancies with two different clocks, and the clinical instruments only
##  see the fast one.  Tab 3 shows the fast pool, tab 4 the slow pool, and tab 8
##  shows what the monitors report -- the three together are the argument.
##
##  Tabs
##  ----
##   1. Patient & exposure    subject, source, enclosure, therapy timing
##   2. Gas exchange & CFK     half-life, transfer resistance, the pressure floor
##   3. Blood CO (fast pool)   COHb, oxygen content, the allosteric left shift
##   4. Tissue CO (slow pool)  myoglobin, cytochrome c oxidase, the two clocks
##   5. Energetics & injury    brain, watershed, myocardium, necrosis, troponin
##   6. Oxidative cascade      xanthine oxidase, ROS, neutrophils, peroxidation
##   7. DNS switch             adduct burden vs threshold, clone, demyelination
##   8. Monitoring model       SpO2 vs SaO2, the saturation gap, lactate
##   9. Scenario comparison    up to 6 arms side by side
##  10. HBO timing             benefit as a function of time to first session
##  11. Virtual population     DNS incidence as a threshold distribution
##  12. Fetal & special        fetus, anaemia, smoker, paediatric, elderly
##
##  Run with:  shiny::runApp("co_shiny_app.R")
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

MODEL_FILE <- "co_mrgsolve_model.R"
mod <- mread(MODEL_FILE)

MINUTE <- 1
HOUR   <- 60
DAY    <- 1440

## exposure scale calibrated by root-finding on peak COHb (see README)
EXPOSURES <- list(
  "Mild  (COHb ~10%)"      = list(ppm =  538, dur = 60),
  "Moderate (COHb ~25%)"   = list(ppm = 1409, dur = 60),
  "Severe (COHb ~40%)"     = list(ppm = 2298, dur = 60),
  "Critical (COHb ~55%)"   = list(ppm = 4198, dur = 45),
  "Custom"                 = list(ppm = 1500, dur = 60)
)

THERAPIES <- c("None (room air)", "Nasal cannula 6 L/min",
               "Non-rebreather mask", "Intubated FiO2 1.0",
               "Intubated + hyperventilation", "Carbogen 95/5",
               "HBO x1 at 3.0 ATA", "HBO x3 at 3.0 ATA (Weaver)",
               "HBO x3 at 2.0 ATA", "HBO x3 delayed 20 h")

theme_co <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}
PAL <- c("#2E6F9E", "#C2603A", "#4E9440", "#8E6FB5", "#C4566A", "#B79A16")

## ---------------------------------------------------------------- simulation
therapy_params <- function(label, t_arrive, o2_hours) {
  p <- list(to2_start = 1e9, to2_stop = 1e9, FiO2_trt = 0.85, VAtrt = 4200,
            FiCO2 = 0, thbo_start = 1e9, ATA = 1, nhbo = 1, hbo_gap = 480)
  if (label == "None (room air)") return(p)
  p$to2_start <- t_arrive
  p$to2_stop  <- t_arrive + o2_hours*HOUR
  if (label == "Nasal cannula 6 L/min")        p$FiO2_trt <- 0.44
  if (label == "Intubated FiO2 1.0")           p$FiO2_trt <- 1.00
  if (label == "Intubated + hyperventilation") { p$FiO2_trt <- 1.00; p$VAtrt <- 8000 }
  if (label == "Carbogen 95/5")                { p$FiO2_trt <- 0.95; p$VAtrt <- 8000
                                                 p$FiCO2 <- 0.05 }
  if (grepl("^HBO", label)) {
    p$to2_stop   <- t_arrive + 48*HOUR
    p$thbo_start <- t_arrive + 30
    p$ATA        <- if (grepl("2.0 ATA", label)) 2.0 else 3.0
    p$nhbo       <- if (grepl("x1", label)) 1 else 3
    if (grepl("delayed", label)) p$thbo_start <- t_arrive + 20*HOUR
  }
  p
}

simulate_arm <- function(ppm, dur, therapy, t_arrive, o2_hours, Hb, WT,
                         cn_rate = 0, ohcbl = 0, nac = 0, allo = 0,
                         hypotherm = 0, days = 45) {
  pp <- c(list(ppm_fix = ppm, texp = dur, Hb = Hb, WT = WT,
               CN_rate = cn_rate, hypotherm = hypotherm),
          therapy_params(therapy, t_arrive, o2_hours))
  m <- mod %>% param(pp)
  ev_list <- list()
  if (ohcbl > 0) ev_list[[length(ev_list)+1]] <-
      ev(time = t_arrive + 5, amt = ohcbl/1355*1e6/18, cmt = "OHCbl")
  if (nac   > 0) ev_list[[length(ev_list)+1]] <-
      ev(time = t_arrive, amt = nac*WT, cmt = "NACc")
  if (allo  > 0) ev_list[[length(ev_list)+1]] <-
      ev(time = t_arrive, amt = allo, cmt = "Oxy")
  if (length(ev_list)) m <- m %>% ev(do.call(c, ev_list))
  m %>% mrgsim(end = days*DAY, delta = 5, hmax = 5) %>% as_tibble()
}

## ---------------------------------------------------------------------- UI
ui <- fluidPage(
  titlePanel("Carbon Monoxide Poisoning — QSP dashboard"),
  tags$p(tags$em(paste(
    "Two occupancies, two clocks: carboxyhaemoglobin is fast and measured;",
    "cytochrome c oxidase is slow and unmeasured. Educational model only —",
    "not for clinical use."))),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Exposure"),
      selectInput("exposure", "Severity preset", names(EXPOSURES),
                  selected = "Severe (COHb ~40%)"),
      conditionalPanel("input.exposure == 'Custom'",
        sliderInput("ppm", "Ambient CO (ppm)", 100, 10000, 1500, 50),
        sliderInput("dur", "Exposure duration (min)", 5, 480, 60, 5)),
      sliderInput("cn_rate", "Cyanide co-exposure (µM/min, fire smoke)",
                  0, 10, 0, 0.5),
      hr(),
      h4("Subject"),
      sliderInput("Hb", "Haemoglobin (g/dL)", 6, 18, 15, 0.5),
      sliderInput("WT", "Weight (kg)", 10, 120, 70, 1),
      hr(),
      h4("Treatment"),
      selectInput("therapy", "Primary therapy", THERAPIES,
                  selected = "Non-rebreather mask"),
      sliderInput("t_arrive", "Time to treatment (min from exposure start)",
                  30, 1440, 90, 10),
      sliderInput("o2_hours", "Normobaric oxygen duration (h)", 0, 48, 6, 1),
      hr(),
      h4("Adjuncts"),
      sliderInput("ohcbl", "Hydroxocobalamin (g)", 0, 10, 0, 1),
      sliderInput("nac", "N-acetylcysteine (mg/kg)", 0, 300, 0, 25),
      sliderInput("allo", "Oxypurinol (mg/L initial)", 0, 40, 0, 5),
      checkboxInput("hypotherm", "Targeted temperature management", FALSE),
      hr(),
      sliderInput("days", "Follow-up (days)", 1, 90, 45, 1),
      sliderInput("zoom", "Acute-phase zoom (h)", 2, 48, 12, 1)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        tabPanel("1. Patient & exposure",
          h4("Exposure and treatment as configured"),
          verbatimTextOutput("summary_txt"),
          plotOutput("exposure_plot", height = "300px"),
          h4("Derived acute indices"), tableOutput("acute_tbl")),

        tabPanel("2. Gas exchange (CFK)",
          h4("Coburn–Forster–Kane: what actually sets the half-life"),
          tags$p(paste(
            "The transfer resistance is B = 1/DLco + PL/VA. At rest about 81%",
            "of B is ventilation and only 19% membrane diffusion, so the",
            "half-life is a ventilation property.")),
          plotOutput("cfk_plot", height = "320px"),
          h4("Half-life by delivery device, and the hyperbaric floor"),
          tableOutput("halflife_tbl"),
          tags$p(tags$strong("The floor: "), paste(
            "raising chamber pressure raises the driving force PcO2 and the",
            "resistance PL/VA in exact proportion, so PcO2/B tends to VA and",
            "the rate constant tends to a ceiling VA/(M[O2Hb]Vb). The",
            "half-life cannot fall below ln2*M*[O2Hb]*Vb/VA = 31.3 min at any",
            "pressure. The observed ~20 min at 2.5–3 ATA is below that",
            "floor, so CFK is falsified in the chamber — reported here as a",
            "miss, not fitted away.")),
          plotOutput("floor_plot", height = "300px")),

        tabPanel("3. Blood CO (fast pool)",
          h4("Carboxyhaemoglobin and the two hits of one ligand"),
          plotOutput("cohb_plot", height = "300px"),
          h4("Oxygen content, dissolved oxygen, and the allosteric left shift"),
          plotOutput("o2_plot", height = "320px"),
          h4("CO against the anaemia with the same arterial oxygen content"),
          tableOutput("anaemia_tbl")),

        tabPanel("4. Tissue CO (slow pool)",
          h4("The pool nobody measures"),
          plotOutput("tissue_plot", height = "330px"),
          h4("Cytochrome c oxidase: CO binds only the reduced a3 haem, so oxygen protects"),
          plotOutput("cco_plot", height = "300px"),
          h4("The two clocks"), tableOutput("clocks_tbl")),

        tabPanel("5. Energetics & injury",
          h4("Brain, watershed and myocardial energy charge"),
          plotOutput("energy_plot", height = "330px"),
          h4("Organ injury markers"),
          plotOutput("organ_plot", height = "330px")),

        tabPanel("6. Oxidative cascade",
          h4("The burst that the treatment creates"),
          tags$p(paste(
            "Reactive oxygen production requires BOTH a prior energy failure",
            "(to convert xanthine dehydrogenase to the oxidase) AND restored",
            "oxygen as co-substrate. The burst is therefore generated by",
            "reoxygenation, which is to say by the treatment.")),
          plotOutput("ox_plot", height = "340px"),
          h4("Innate amplification and glutathione"),
          plotOutput("innate_plot", height = "320px")),

        tabPanel("7. DNS switch",
          h4("Delayed neurological sequelae as a bistable switch"),
          verbatimTextOutput("dns_txt"),
          plotOutput("dns_plot", height = "340px"),
          h4("Adduct burden against the computed threshold"),
          plotOutput("thresh_plot", height = "300px")),

        tabPanel("8. Monitoring model",
          h4("What the monitors report against what is true"),
          tags$p(paste(
            "COHb absorbs at 660 nm almost exactly as oxyhaemoglobin does and",
            "is invisible at 940 nm, so a two-wavelength oximeter attributes",
            "the CO-occupied fraction to oxygenated haemoglobin. SpO2 reports",
            "(O2Hb + COHb): the more CO the patient carries, the more",
            "reassuring the number becomes.")),
          plotOutput("monitor_plot", height = "330px"),
          h4("Saturation gap across the COHb range"),
          plotOutput("gap_plot", height = "300px"),
          h4("Lactate — the only bedside discriminator between CO and cyanide"),
          plotOutput("lac_plot", height = "280px")),

        tabPanel("9. Scenario comparison",
          h4("Treatment arms side by side"),
          checkboxGroupInput("arms", "Arms to compare", THERAPIES,
            selected = c("None (room air)", "Non-rebreather mask",
                         "HBO x3 at 3.0 ATA (Weaver)", "HBO x3 delayed 20 h"),
            inline = TRUE),
          plotOutput("compare_plot", height = "480px"),
          tableOutput("compare_tbl")),

        tabPanel("10. HBO timing",
          h4("The window is the adduct, not the carboxyhaemoglobin"),
          tags$p(paste(
            "Benefit is near-maximal within the first hour, has halved by",
            "6–8 h and is gone by 24 h — while COHb is already below 5%",
            "in every arm beyond about 6 h. Weaver 2002 delivered the first",
            "session at a median of roughly 4 h and was positive;",
            "Scheinkestel 1999 treated after ICU stabilisation, frequently",
            "beyond 12 h, and was null.")),
          sliderInput("hbo_grid", "Latest first-session time to test (h)",
                      6, 48, 36, 2),
          plotOutput("timing_plot", height = "400px"),
          tableOutput("timing_tbl")),

        tabPanel("11. Virtual population",
          h4("The clinical incidence curve is a distribution of thresholds"),
          tags$p(paste(
            "Within one patient the response is a step. The graded incidence",
            "seen in cohorts is therefore between-patient variation in where",
            "that step sits, not a graded dose-response within anyone.")),
          sliderInput("npop", "Population size", 20, 400, 120, 20),
          actionButton("runpop", "Simulate population", class = "btn-primary"),
          plotOutput("pop_plot", height = "380px"),
          tableOutput("pop_tbl")),

        tabPanel("12. Fetal & special",
          h4("The fetal compartment sets the treatment duration"),
          tags$p(paste(
            "Fetal haemoglobin binds CO about 1.8x more avidly and the fetal",
            "compartment equilibrates slowly, so the fetus lags, peaks higher",
            "and clears far more slowly. Maternal COHb reaches 5% at 204 min",
            "on a mask while the fetal compartment does not until 755 min —",
            "which is where the bedside 'five times' rule comes from.")),
          plotOutput("fetal_plot", height = "340px"),
          h4("Haemoglobin sweep: the same COHb percentage is not the same disease"),
          plotOutput("hb_plot", height = "320px"))
      )
    )
  )
)

## ------------------------------------------------------------------- server
server <- function(input, output, session) {

  expo <- reactive({
    e <- EXPOSURES[[input$exposure]]
    if (input$exposure == "Custom") list(ppm = input$ppm, dur = input$dur) else e
  })

  sim <- reactive({
    e <- expo()
    simulate_arm(e$ppm, e$dur, input$therapy, input$t_arrive, input$o2_hours,
                 input$Hb, input$WT, cn_rate = input$cn_rate,
                 ohcbl = input$ohcbl, nac = input$nac, allo = input$allo,
                 hypotherm = as.numeric(input$hypotherm), days = input$days)
  })

  zoom <- reactive(function(d) filter(d, time <= input$zoom*HOUR))

  ## ---- tab 1 -------------------------------------------------------------
  output$summary_txt <- renderText({
    e <- expo(); d <- sim()
    tp <- therapy_params(input$therapy, input$t_arrive, input$o2_hours)
    paste0(
      sprintf("Exposure    : %.0f ppm for %.0f min  (Haber dose %.0f ppm-min)\n",
              e$ppm, e$dur, e$ppm*e$dur),
      sprintf("Subject     : Hb %.1f g/dL, %.0f kg\n", input$Hb, input$WT),
      sprintf("Therapy     : %s, starting %.0f min after exposure began\n",
              input$therapy, input$t_arrive),
      if (tp$thbo_start < 1e8)
        sprintf("Hyperbaric  : %.0f session(s) at %.1f ATA, first at %.0f min\n",
                tp$nhbo, tp$ATA, tp$thbo_start) else "",
      sprintf("Peak COHb   : %.1f%%   Peak fetal COHb : %.1f%%\n",
              max(d$COHb_pct), max(d$FetCO_pct)),
      sprintf("Peak CcO inhibition (brain) : %.1f%%\n", 100*max(d$CcOb)),
      sprintf("Peak MBP adduct burden      : %.3f  (threshold 0.4425)\n",
              max(d$MBPad)),
      sprintf("Demyelination at day %.0f     : %.3f  -> DNS %s\n",
              input$days, tail(d$Demy, 1),
              ifelse(tail(d$Demy, 1) > 0.05, "YES", "no")),
      sprintf("Cognitive score at day %.0f  : %.1f%% of baseline",
              input$days, 100*tail(d$Cog, 1)))
  })

  output$exposure_plot <- renderPlot({
    e <- expo(); d <- zoom()(sim())
    d %>% transmute(time = time/HOUR,
                    `Ambient CO (ppm/100)` = ifelse(time*HOUR < e$dur, e$ppm/100, 0),
                    `COHb (%)` = COHb_pct,
                    `Fetal COHb (%)` = FetCO_pct) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.05) +
      scale_colour_manual(values = PAL) +
      labs(x = "hours from start of exposure", y = NULL, colour = NULL,
           title = "Exposure and the fast pool") + theme_co()
  })

  output$acute_tbl <- renderTable({
    d <- sim()
    tibble(
      Index = c("Peak COHb (%)", "Peak SpO2 displayed (%)",
                "Minimum true SaO2 (%)", "Peak saturation gap (points)",
                "Minimum CaO2 (mL/dL)", "Minimum tissue PO2 (mmHg)",
                "Peak lactate (mM)", "Peak troponin I (ng/mL)",
                "Peak creatine kinase (U/L)", "Peak ICP (mmHg)"),
      Value = c(sprintf("%.1f", max(d$COHb_pct)),
                sprintf("%.1f", max(d$SpO2_displayed)),
                sprintf("%.1f", min(d$SaO2_pct)),
                sprintf("%.1f", max(d$sat_gap)),
                sprintf("%.2f", min(d$CaO2)),
                sprintf("%.1f", min(d$PtO2_mmHg)),
                sprintf("%.1f", max(d$Lac)),
                sprintf("%.1f", max(d$TnI)),
                sprintf("%.0f", max(d$CK)),
                sprintf("%.1f", max(d$ICP))))
  })

  ## ---- tab 2 -------------------------------------------------------------
  cfk_halflife <- function(PcO2, VA = 4200, ATA = 1, Hb = 15, Vb = 5500,
                           DL = 25, M = 245, F0 = 0.30) {
    PL <- 760*ATA - 47; cap <- 1.34*Hb/100
    B  <- 1/DL + PL/VA
    ## closed-form on the linearised equation is enough for display
    O2 <- cap*(1 - F0/2)
    k  <- PcO2/(M*O2*Vb*B)
    log(2)/k
  }

  output$cfk_plot <- renderPlot({
    grid <- expand_grid(VA = seq(2500, 12000, 250),
                        dev = c("room air", "non-rebreather", "FiO2 1.0")) %>%
      mutate(PcO2 = c("room air" = 100, "non-rebreather" = 432,
                      "FiO2 1.0" = 560)[dev],
             t12 = mapply(cfk_halflife, PcO2, VA))
    ggplot(grid, aes(VA/1000, t12, colour = dev)) +
      geom_line(linewidth = 1.05) +
      geom_vline(xintercept = 4.2, linetype = 3) +
      scale_colour_manual(values = PAL) + scale_y_log10() +
      labs(x = "alveolar ventilation (L/min)", y = "COHb half-life (min, log)",
           colour = NULL, title = "The half-life is a ventilation property") +
      theme_co()
  })

  output$halflife_tbl <- renderTable({
    tibble(Delivery = c("Room air", "Nasal cannula 6 L/min", "Non-rebreather mask",
                        "Intubated FiO2 1.0", "HBO 2.0 ATA", "HBO 2.4 ATA",
                        "HBO 3.0 ATA", "Analytic floor (any pressure)"),
           `Model t1/2 (min)` = c("312.9", "124.1", "71.7", "55.2",
                                  "40.5", "39.5", "38.6", "31.3"),
           `Observed (min)`   = c("320", "—", "74", "—",
                                  "~20–25", "~20", "~20", "—"))
  })

  output$floor_plot <- renderPlot({
    ata <- seq(1, 12, 0.1)
    tibble(ata,
           t12 = mapply(function(a) cfk_halflife(760*a - 47 - 50, ATA = a), ata)) %>%
      ggplot(aes(ata)) +
      geom_line(aes(y = t12, colour = "CFK prediction"), linewidth = 1.1) +
      geom_hline(aes(yintercept = 31.3, colour = "analytic floor 31.3 min"),
                 linetype = 2, linewidth = 0.9) +
      geom_hline(aes(yintercept = 20, colour = "observed at 2.5–3 ATA"),
                 linetype = 3, linewidth = 0.9) +
      scale_colour_manual(values = c("CFK prediction" = PAL[1],
                                     "analytic floor 31.3 min" = PAL[2],
                                     "observed at 2.5–3 ATA" = PAL[3])) +
      labs(x = "chamber pressure (ATA)", y = "COHb half-life (min)",
           colour = NULL,
           title = "Pressure buys almost nothing, and the data sit below the floor") +
      theme_co()
  })

  ## ---- tab 3 -------------------------------------------------------------
  output$cohb_plot <- renderPlot({
    zoom()(sim()) %>%
      transmute(time = time/HOUR, `COHb (%)` = COHb_pct,
                `true SaO2 (%)` = SaO2_pct,
                `SpO2 displayed (%)` = SpO2_displayed) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.05) + scale_colour_manual(values = PAL) +
      labs(x = "hours", y = "%", colour = NULL,
           title = "The fast pool, and what the monitor says about it") + theme_co()
  })

  output$o2_plot <- renderPlot({
    zoom()(sim()) %>%
      transmute(time = time/HOUR, `CaO2 (mL/dL)` = CaO2,
                `dissolved O2 (mL/dL)` = O2_diss,
                `P50 effective (mmHg)` = P50_eff,
                `tissue PO2 (mmHg)` = PtO2_mmHg) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.05) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(x = "hours", y = NULL,
           title = "Content, dissolved oxygen and the left shift") + theme_co()
  })

  output$anaemia_tbl <- renderTable({
    P50 <- 26.8; n <- 2.7
    inv <- function(s, p50) p50*(s/(1-s))^(1/n)
    Ca  <- function(Hb, F) 1.34*Hb*0.97*(1-F) + 0.003*95
    pv  <- function(Hb, F) {
      cap <- 1.34*Hb*(1-F); Cv <- Ca(Hb, F) - 5
      inv(min(max((Cv - 0.12)/cap, 1e-6), 1-1e-6), P50*(1 - 0.75*F))
    }
    Fv <- c(0, .1, .2, .3, .4, .5)
    Hbeq <- sapply(Fv, function(F)
      uniroot(function(h) Ca(h, 0) - Ca(15, F), c(0.5, 25))$root)
    tibble(`COHb (%)` = sprintf("%.0f", 100*Fv),
           `CaO2 (mL/dL)` = sprintf("%.2f", sapply(Fv, function(F) Ca(15, F))),
           `equivalent anaemia (g/dL)` = sprintf("%.2f", Hbeq),
           `tissue PO2 with CO` = sprintf("%.1f", sapply(Fv, function(F) pv(15, F))),
           `tissue PO2 of that anaemia` = sprintf("%.1f", sapply(Hbeq, function(h) pv(h, 0))),
           `deficit (mmHg)` = sprintf("%.1f",
              sapply(seq_along(Fv), function(i) pv(Hbeq[i], 0) - pv(15, Fv[i]))))
  })

  ## ---- tab 4 -------------------------------------------------------------
  output$tissue_plot <- renderPlot({
    d <- filter(sim(), time <= min(input$days*DAY, 3*DAY))
    d %>% transmute(time = time/HOUR,
                    `brain tissue CO` = Pbr/max(Pbr),
                    `myocardial CO` = Pht/max(Pht),
                    `skeletal muscle CO` = Pmu/max(Pmu),
                    `COHb` = COHb_pct/max(COHb_pct)) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.05) + scale_colour_manual(values = PAL) +
      labs(x = "hours", y = "fraction of own peak", colour = NULL,
           title = "Normalised to each peak: the reservoirs empty on different clocks") +
      theme_co()
  })

  output$cco_plot <- renderPlot({
    d <- filter(sim(), time <= min(input$days*DAY, 10*DAY))
    d %>% transmute(time = time/HOUR,
                    `brain CcO inhibited (%)` = 100*CcOb,
                    `myocardial CcO inhibited (%)` = 100*CcOh,
                    `total incl. cyanide (%)` = CcO_total_pct,
                    `cardiac MbCO (%)` = 100*MbHt) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.05) + scale_colour_manual(values = PAL) +
      labs(x = "hours", y = "%", colour = NULL,
           title = "The slow pool") + theme_co()
  })

  output$clocks_tbl <- renderTable({
    d <- sim()
    tfrac <- function(v, fr) {
      i <- which.max(v); j <- which(v[i:length(v)] < fr*max(v))
      if (!length(j)) NA_real_ else d$time[i + j[1] - 1]
    }
    tibble(Pool = c("COHb (measured)", "brain tissue CO", "cardiac MbCO",
                    "skeletal muscle CO", "brain cytochrome c oxidase"),
           `peak` = sprintf("%.3g", c(max(d$COHb_pct), max(d$Pbr), max(d$MbHt),
                                      max(d$Pmu), max(d$CcOb))),
           `t(50% of peak), min` = sprintf("%.0f",
              c(tfrac(d$COHb_pct,.5), tfrac(d$Pbr,.5), tfrac(d$MbHt,.5),
                tfrac(d$Pmu,.5), tfrac(d$CcOb,.5))),
           `t(10% of peak), min` = sprintf("%.0f",
              c(tfrac(d$COHb_pct,.1), tfrac(d$Pbr,.1), tfrac(d$MbHt,.1),
                tfrac(d$Pmu,.1), tfrac(d$CcOb,.1))))
  })

  ## ---- tab 5 -------------------------------------------------------------
  output$energy_plot <- renderPlot({
    filter(sim(), time <= min(input$days*DAY, 2*DAY)) %>%
      transmute(time = time/HOUR, `brain` = ATPb, `watershed (GP)` = ATPgp,
                `myocardium` = ATPh) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.05) +
      geom_hline(yintercept = 0.42, linetype = 2, colour = "grey40") +
      annotate("text", x = Inf, y = 0.42, hjust = 1.05, vjust = -0.5,
               label = "necrosis threshold", size = 3, colour = "grey30") +
      scale_colour_manual(values = PAL) +
      labs(x = "hours", y = "energy charge (fraction of baseline)", colour = NULL,
           title = "The watershed is not handicapped at baseline — it is steeper") +
      theme_co()
  })

  output$organ_plot <- renderPlot({
    sim() %>% transmute(time = time/DAY, `necrosis` = Necb, `oedema` = Edema,
                        `troponin I (ng/mL)` = TnI, `EF` = EF,
                        `CK (kU/L)` = CK/1000, `creatinine (mg/dL)` = Cr) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 1.0, colour = PAL[2]) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "days", y = NULL, title = "Organ injury") + theme_co()
  })

  ## ---- tab 6 -------------------------------------------------------------
  output$ox_plot <- renderPlot({
    filter(sim(), time <= min(input$days*DAY, 3*DAY)) %>%
      transmute(time = time/HOUR, `xanthine oxidase` = XO, `ROS` = ROS,
                `NO / peroxynitrite` = NOx, `lipid peroxidation` = LPO) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.05) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(x = "hours", y = NULL,
           title = "The oxidative burst needs prior failure AND restored oxygen") +
      theme_co()
  })

  output$innate_plot <- renderPlot({
    filter(sim(), time <= min(input$days*DAY, 5*DAY)) %>%
      transmute(time = time/HOUR, `adherent neutrophils` = Neut,
                `myeloperoxidase` = MPO, `glutathione` = GSH,
                `microglia` = Micro) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.05) + facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL, guide = "none") +
      labs(x = "hours", y = NULL,
           title = "Innate amplification; hyperbaric oxygen blocks the first step directly") +
      theme_co()
  })

  ## ---- tab 7 -------------------------------------------------------------
  output$dns_txt <- renderText({
    d <- sim(); Mc <- 0.4425
    paste0(
      "The autoreactive clone expands if and only if  Tprol*H(MBPad) > Tdeath.\n",
      "With Tdeath/Tprol = 0.4750, theta = 0.450 and n = 6 this inverts to\n",
      "   MBPad_crit = theta*((Tdeath/Tprol)/(1 - Tdeath/Tprol))^(1/n) = 0.4425\n",
      "Demyelination liberates further antigen (epitope spreading), so above the\n",
      "threshold the loop LATCHES.  The weeks-long latency and the all-or-none\n",
      "character are consequences of that structure, not assumptions.\n\n",
      sprintf("This patient: peak MBPad = %.3f  =  %.2f x threshold  ->  %s\n",
              max(d$MBPad), max(d$MBPad)/Mc,
              ifelse(max(d$MBPad) > Mc, "clone expands", "clone contracts")),
      sprintf("Demyelination at day %.0f = %.3f;  cognitive score %.1f%% of baseline",
              input$days, tail(d$Demy, 1), 100*tail(d$Cog, 1)))
  })

  output$dns_plot <- renderPlot({
    sim() %>% transmute(time = time/DAY, `MBP adduct` = MBPad,
                        `autoreactive clone` = Tcell,
                        `demyelination` = Demy,
                        `cognitive score` = Cog) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.1) + scale_colour_manual(values = PAL) +
      labs(x = "days", y = NULL, colour = NULL,
           title = "The lucid interval is the time the switch takes to climb") +
      theme_co()
  })

  output$thresh_plot <- renderPlot({
    d <- sim()
    ggplot(d, aes(time/DAY, MBPad)) +
      geom_line(linewidth = 1.1, colour = PAL[1]) +
      geom_hline(yintercept = 0.4425, linetype = 2, colour = PAL[2], linewidth = 1) +
      annotate("text", x = Inf, y = 0.4425, hjust = 1.05, vjust = -0.6,
               label = "MBPad_crit = 0.4425", colour = PAL[2], size = 3.5) +
      labs(x = "days", y = "adduct burden",
           title = "Time spent above the tolerance threshold is what decides DNS") +
      theme_co()
  })

  ## ---- tab 8 -------------------------------------------------------------
  output$monitor_plot <- renderPlot({
    zoom()(sim()) %>%
      transmute(time = time/HOUR, `SpO2 displayed` = SpO2_displayed,
                `true SaO2` = SaO2_pct, `saturation gap` = sat_gap,
                `PaO2/10 (mmHg)` = PaO2_mmHg/10) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.05) + scale_colour_manual(values = PAL) +
      labs(x = "hours", y = NULL, colour = NULL,
           title = "The sicker the patient, the more reassuring the number") +
      theme_co()
  })

  output$gap_plot <- renderPlot({
    e660 <- c(O2 = .081, HH = .845, CO = .083)
    e940 <- c(O2 = .290, HH = .170, CO = .010)
    f <- function(FCO) {
      fO2 <- (1-FCO)*.97; fHH <- (1-FCO)*.03
      R <- (fO2*e660["O2"] + fHH*e660["HH"] + FCO*e660["CO"]) /
           (fO2*e940["O2"] + fHH*e940["HH"] + FCO*e940["CO"])
      s <- (e660["HH"] - R*e940["HH"]) /
           ((e660["HH"]-e660["O2"]) - R*(e940["HH"]-e940["O2"]))
      c(sp = unname(s), tru = fO2)
    }
    FCO <- seq(0, .7, .01)
    m <- t(sapply(FCO, f))
    tibble(COHb = 100*FCO, `SpO2 displayed` = 100*m[, "sp"],
           `true SaO2` = 100*m[, "tru"], `gap` = 100*(m[,"sp"] - m[,"tru"])) %>%
      pivot_longer(-COHb) %>%
      ggplot(aes(COHb, value, colour = name)) +
      geom_line(linewidth = 1.1) + scale_colour_manual(values = PAL) +
      labs(x = "true COHb (%)", y = "%", colour = NULL,
           title = "SpO2 tracks (O2Hb + COHb) almost exactly") + theme_co()
  })

  output$lac_plot <- renderPlot({
    filter(sim(), time <= min(input$days*DAY, 2*DAY)) %>%
      transmute(time = time/HOUR, `lactate (mM)` = Lac,
                `blood cyanide (µM)` = CNbl,
                `CcO inhibited total (%)` = CcO_total_pct/10) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) +
      geom_line(linewidth = 1.05) + scale_colour_manual(values = PAL) +
      labs(x = "hours", y = NULL, colour = NULL,
           title = "A lactate that will not fall on 100% oxygen means cyanide") +
      theme_co()
  })

  ## ---- tab 9 -------------------------------------------------------------
  comparison <- reactive({
    e <- expo()
    req(length(input$arms) > 0)
    lapply(input$arms, function(a) {
      simulate_arm(e$ppm, e$dur, a, input$t_arrive, input$o2_hours,
                   input$Hb, input$WT, cn_rate = input$cn_rate,
                   days = input$days) %>% mutate(arm = a)
    }) %>% bind_rows()
  })

  output$compare_plot <- renderPlot({
    comparison() %>%
      transmute(arm, time = time/HOUR, `COHb (%)` = COHb_pct,
                `brain CcO inhibited (%)` = 100*CcOb,
                `lipid peroxidation` = LPO, `MBP adduct` = MBPad,
                `demyelination` = Demy, `cognitive score` = Cog) %>%
      pivot_longer(-c(arm, time)) %>%
      ggplot(aes(time, value, colour = arm)) +
      geom_line(linewidth = 1.0) +
      facet_wrap(~name, scales = "free", ncol = 3) +
      scale_x_log10() + scale_colour_manual(values = rep(PAL, 3)) +
      labs(x = "hours (log)", y = NULL, colour = NULL,
           title = "Treatment arms") + theme_co()
  })

  output$compare_tbl <- renderTable({
    comparison() %>% group_by(arm) %>%
      summarise(`peak COHb (%)` = sprintf("%.1f", max(COHb_pct)),
                `peak CcO (%)`  = sprintf("%.1f", 100*max(CcOb)),
                `peak MBPad`    = sprintf("%.3f", max(MBPad)),
                `demyelination` = sprintf("%.3f", last(Demy)),
                `necrosis`      = sprintf("%.3f", last(Necb)),
                `cognitive (%)` = sprintf("%.1f", 100*last(Cog)),
                `DNS`           = ifelse(last(Demy) > 0.05, "YES", "–"),
                .groups = "drop")
  })

  ## ---- tab 10 ------------------------------------------------------------
  timing <- reactive({
    e <- expo()
    hrs <- c(1, 1.5, 2, 3, 4, 6, 8, 12, 18, 24, 36)
    hrs <- hrs[hrs <= input$hbo_grid]
    base <- simulate_arm(e$ppm, e$dur, "Non-rebreather mask", input$t_arrive,
                         48, input$Hb, input$WT, days = input$days)
    ref <- tail(base$Cog, 1)
    lapply(hrs, function(h) {
      pp <- c(list(ppm_fix = e$ppm, texp = e$dur, Hb = input$Hb),
              list(to2_start = input$t_arrive, to2_stop = input$t_arrive + 48*HOUR,
                   thbo_start = h*HOUR, ATA = 3, nhbo = 3, hbo_gap = 480))
      d <- mod %>% param(pp) %>% mrgsim(end = input$days*DAY, delta = 15, hmax = 5) %>%
        as_tibble()
      tibble(hbo_h = h, cohb_then = approx(d$time, d$COHb_pct, h*HOUR)$y,
             mbp = max(d$MBPad), demy = tail(d$Demy, 1), cog = tail(d$Cog, 1),
             recovered = 100*(tail(d$Cog, 1) - ref)/max(1e-9, 1 - ref))
    }) %>% bind_rows() %>% mutate(ref_cog = ref)
  })

  output$timing_plot <- renderPlot({
    timing() %>%
      select(hbo_h, `COHb at treatment (%)` = cohb_then,
             `peak MBP adduct` = mbp, `cognitive score` = cog,
             `% of loss recovered` = recovered) %>%
      pivot_longer(-hbo_h) %>%
      ggplot(aes(hbo_h, value)) +
      geom_line(linewidth = 1.1, colour = PAL[1]) +
      geom_point(size = 1.8, colour = PAL[1]) +
      facet_wrap(~name, scales = "free_y") +
      geom_vline(xintercept = 6, linetype = 3) +
      labs(x = "time of first hyperbaric session (h)", y = NULL,
           title = "Benefit decays with the adduct, not with the COHb") + theme_co()
  })

  output$timing_tbl <- renderTable({
    timing() %>% transmute(`HBO start (h)` = sprintf("%.1f", hbo_h),
                           `COHb then (%)` = sprintf("%.1f", cohb_then),
                           `peak MBPad` = sprintf("%.3f", mbp),
                           `demyelination` = sprintf("%.3f", demy),
                           `cognitive` = sprintf("%.3f", cog),
                           `% recovered` = sprintf("%.1f", recovered))
  })

  ## ---- tab 11 ------------------------------------------------------------
  pop <- eventReactive(input$runpop, {
    e <- expo(); n <- input$npop
    set.seed(20260806)
    withProgress(message = "Simulating virtual population", value = 0, {
      lapply(seq_len(n), function(i) {
        incProgress(1/n)
        pp <- c(list(ppm_fix = e$ppm, texp = e$dur,
                     Hb = max(8, min(18, rnorm(1, 14.2, 1.7))),
                     VA = max(3000, min(5600, rnorm(1, 4200, 500))),
                     kAd = 0.00050*exp(rnorm(1, 0, .30)),
                     Tprol = 0.00040*exp(rnorm(1, 0, .22)),
                     kon_cco = 0.10*exp(rnorm(1, 0, .25)),
                     CBF0 = max(32, min(68, rnorm(1, 50, 7))),
                     kXOon = 0.028*exp(rnorm(1, 0, .30))),
                therapy_params(input$therapy, input$t_arrive, input$o2_hours))
        d <- mod %>% param(pp) %>% mrgsim(end = 45*DAY, delta = 60, hmax = 5) %>%
          as_tibble()
        tibble(id = i, cohb = max(d$COHb_pct), mbp = max(d$MBPad),
               demy = tail(d$Demy, 1), cog = tail(d$Cog, 1))
      }) %>% bind_rows()
    })
  })

  output$pop_plot <- renderPlot({
    p <- pop()
    ggplot(p, aes(mbp, cog, colour = demy > 0.05)) +
      geom_point(size = 2.2, alpha = 0.8) +
      geom_vline(xintercept = 0.4425, linetype = 2, colour = "grey30") +
      annotate("text", x = 0.4425, y = Inf, vjust = 1.5, hjust = -0.05,
               label = "MBPad_crit", size = 3.4) +
      scale_colour_manual(values = c("FALSE" = PAL[1], "TRUE" = PAL[2]),
                          labels = c("no DNS", "DNS"), name = NULL) +
      labs(x = "peak MBP adduct burden", y = "cognitive score at 45 days",
           title = "One threshold, many patients: the incidence curve is its distribution") +
      theme_co()
  })

  output$pop_tbl <- renderTable({
    p <- pop()
    tibble(Statistic = c("n", "mean peak COHb (%)", "DNS incidence (%)",
                         "mean cognitive score", "mean score, DNS patients",
                         "mean score, non-DNS patients"),
           Value = c(sprintf("%d", nrow(p)),
                     sprintf("%.1f", mean(p$cohb)),
                     sprintf("%.1f", 100*mean(p$demy > 0.05)),
                     sprintf("%.3f", mean(p$cog)),
                     sprintf("%.3f", ifelse(any(p$demy > .05),
                                            mean(p$cog[p$demy > .05]), NA)),
                     sprintf("%.3f", ifelse(any(p$demy <= .05),
                                            mean(p$cog[p$demy <= .05]), NA))))
  })

  ## ---- tab 12 ------------------------------------------------------------
  output$fetal_plot <- renderPlot({
    e <- expo()
    lapply(c("None (room air)", "Non-rebreather mask", "Intubated FiO2 1.0"),
           function(a) simulate_arm(e$ppm, e$dur, a, input$t_arrive, 48,
                                    input$Hb, input$WT, days = 3) %>%
             mutate(arm = a)) %>% bind_rows() %>%
      select(arm, time, `maternal COHb (%)` = COHb_pct,
             `fetal COHb (%)` = FetCO_pct) %>%
      pivot_longer(-c(arm, time)) %>%
      ggplot(aes(time/HOUR, value, colour = name)) +
      geom_line(linewidth = 1.05) + facet_wrap(~arm) +
      geom_hline(yintercept = 5, linetype = 2, colour = "grey40") +
      scale_colour_manual(values = PAL) +
      labs(x = "hours", y = "%", colour = NULL,
           title = "Treat to the fetal clock, not to the maternal COHb") + theme_co()
  })

  output$hb_plot <- renderPlot({
    e <- expo()
    lapply(c(8, 10, 12, 15, 17), function(h)
      simulate_arm(e$ppm, e$dur, input$therapy, input$t_arrive, input$o2_hours,
                   h, input$WT, days = input$days) %>% mutate(Hb = h)) %>%
      bind_rows() %>% group_by(Hb) %>%
      summarise(`peak COHb (%)` = max(COHb_pct),
                `min CaO2 (mL/dL)` = min(CaO2),
                `min tissue PO2 (mmHg)` = min(PtO2_mmHg),
                `peak MBP adduct` = max(MBPad),
                `cognitive score` = last(Cog), .groups = "drop") %>%
      pivot_longer(-Hb) %>%
      ggplot(aes(Hb, value)) + geom_line(linewidth = 1.1, colour = PAL[4]) +
      geom_point(size = 2, colour = PAL[4]) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "haemoglobin (g/dL)", y = NULL,
           title = "The same exposure is a different disease at a different haemoglobin") +
      theme_co()
  })
}

shinyApp(ui, server)
