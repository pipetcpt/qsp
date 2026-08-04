## =====================================================================
##  MALIGNANT HYPERTHERMIA — QSP Shiny dashboard
##  Front-end for mh_mrgsolve_model.R (39 ODEs, 18 scenarios)
##
##  The app is deliberately organised around the model's central claim:
##  ONE driver (myoplasmic Ca2+) read through TWO stores of very different
##  size. Tab 3 shows the two stores side by side on the same time axis so
##  that the reason EtCO2 precedes fever is visible rather than asserted,
##  and tab 6 lets the user reproduce the two "normal monitor, dead
##  patient" arms by hand.
##
##  Run with:   shiny::runApp("mh_shiny_app.R")
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("mh_mrgsolve_model.R", local = TRUE)   # defines mod, GENO, AGENT, scen, ...

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

PAL <- c("#1565C0", "#C62828", "#2E7D32", "#EF6C00", "#6A1B9A",
         "#00838F", "#AD1457", "#455A64")

# ---------------------------------------------------------------------
#  Simulation driver used by every tab
# ---------------------------------------------------------------------
simulate_mh <- function(genotype, agent, mac_mult, sux, wt,
                        t_recognise, dantrolene_delay, dantrolene_mgkg,
                        n_bolus, maintenance,
                        do_stop_agent, do_charcoal, vent_mult, cool_W, surf_W,
                        ca_o = 1, nac = 0, end = 480) {

  ag  <- AGENT[[agent]]
  td  <- t_recognise + dantrolene_delay
  pp  <- c(GENO[[genotype]],
           list(WT = wt, POTENCY = ag$POTENCY, VAP = ag$MAC * mac_mult,
                CA_O = ca_o, NAC = nac,
                TVAPOFF = if (do_stop_agent) t_recognise else 1e9,
                TCHARC  = if (do_charcoal)   t_recognise else 1e9,
                TVENT   = if (vent_mult > 1) t_recognise else 1e9,
                VENTMLT = vent_mult,
                TCOOL   = if (cool_W + surf_W > 0) t_recognise else 1e9,
                COOL_W  = cool_W, SURF_W = surf_W))

  e <- ev(time = 0, amt = 0, cmt = "DAN1")
  if (sux)  e <- c(e, ev(time = 2, amt = 1.0 * wt, cmt = "SCHA"))
  if (dantrolene_mgkg > 0 && n_bolus > 0) {
    e <- c(e, ev(time = td, amt = dantrolene_mgkg * wt, cmt = "DAN1",
                 ii = 6, addl = n_bolus - 1))
    if (maintenance)
      e <- c(e, ev(time = td + 60, amt = 1.0 * wt, cmt = "DAN1",
                   ii = 360, addl = 7))
  }
  param(mod, pp) %>% mrgsim(events = e, end = end, delta = 0.5) %>% as_tibble()
}

first_lethal <- function(o) {
  i <- which(o$LETHAL > 0)
  if (length(i)) o$time[i[1]] else NA_real_
}

kpi_box <- function(label, value, sub = "", colour = "#1565C0") {
  div(style = paste0("flex:1;min-width:150px;border-left:5px solid ", colour,
                     ";background:#FFFFFF;padding:10px 14px;margin:6px;",
                     "border-radius:6px;box-shadow:0 1px 3px rgba(0,0,0,.12);"),
      div(style = "font-size:11px;color:#5A6270;text-transform:uppercase;
                   letter-spacing:.04em;", label),
      div(style = paste0("font-size:24px;font-weight:700;color:", colour, ";"), value),
      div(style = "font-size:11px;color:#78838F;", sub))
}

long_plot <- function(o, vars, labels, title, subtitle = NULL, vline = NULL) {
  d <- o %>% select(time, all_of(vars)) %>%
    pivot_longer(-time, names_to = "v", values_to = "y") %>%
    mutate(v = factor(v, levels = vars, labels = labels))
  p <- ggplot(d, aes(time, y, colour = v)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~v, scales = "free_y") +
    scale_colour_manual(values = rep(PAL, 6), guide = "none") +
    labs(title = title, subtitle = subtitle, x = "time (min)", y = NULL) + THEME
  if (!is.null(vline))
    p <- p + geom_vline(xintercept = vline, linetype = 2, colour = "#B71C1C")
  p
}

# =====================================================================
#  UI
# =====================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body{background:#F6F7F9;}
    .kpirow{display:flex;flex-wrap:wrap;margin-bottom:8px;}
    .note{background:#FFF8E1;border-left:4px solid #FF6F00;padding:10px 14px;
          border-radius:4px;font-size:13px;margin:8px 0;}
    .thesis{background:#E3F2FD;border-left:4px solid #0D47A1;padding:10px 14px;
          border-radius:4px;font-size:13px;margin:8px 0;}
    .warn{background:#FFEBEE;border-left:4px solid #B71C1C;padding:10px 14px;
          border-radius:4px;font-size:13px;margin:8px 0;}
  "))),

  titlePanel("Malignant Hyperthermia — QSP dashboard"),
  div(class = "thesis",
      strong("One driver, two stores. "),
      "Myoplasmic Ca", tags$sup("2+"), " drives everything. The body's CO",
      tags$sub("2"), " store is 60 mL/mmHg and its core heat store is 158 kJ/°C; ",
      "divided by the same excess flux those give 0.15 min per mmHg of EtCO",
      tags$sub("2"), " and 17 min per °C of core temperature. That ratio — not ",
      "convention — is why capnography is the first alarm and why the disease ",
      "is named after the last thing that happens."),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("genotype", "Genotype",
                  choices = c("Non-susceptible (MHN)"        = "MHN",
                              "RYR1 low-penetrance"          = "MHS_low",
                              "CACNA1S p.Arg1086His"         = "MHS_CACNA1S",
                              "RYR1 high-penetrance (R614C)" = "MHS_high",
                              "RYR1-related myopathy"        = "RYR1_myopathy"),
                  selected = "MHS_high"),
      sliderInput("wt", "Weight (kg)", 20, 140, 70, step = 5),

      h4("Anaesthetic"),
      selectInput("agent", "Volatile agent",
                  choices = names(AGENT), selected = "sevoflurane"),
      sliderInput("mac", "MAC multiple", 0, 2, 1, step = 0.1),
      checkboxInput("sux", "Succinylcholine 1 mg/kg at induction", TRUE),

      h4("Recognition and response"),
      sliderInput("trec", "Event recognised at (min)", 5, 150, 45, step = 5),
      checkboxInput("stopagent", "Close the vaporiser + 10 L/min O₂", TRUE),
      checkboxInput("charcoal", "Activated-charcoal filters", TRUE),
      sliderInput("vent", "Minute-ventilation multiplier", 1, 4, 2, step = 0.5),
      sliderInput("coolW", "Core (lavage/intravascular) cooling (W)", 0, 400, 200, step = 25),
      sliderInput("surfW", "Surface cooling (W)", 0, 300, 100, step = 25),

      h4("Dantrolene"),
      sliderInput("dandelay", "Delay from recognition to first dose (min)",
                  0, 120, 2, step = 1),
      helpText("2 min ≈ Ryanodex (250 mg / 5 mL). 18 min ≈ Dantrium/Revonto ",
               "(9 vials of 20 mg, each in 60 mL of sterile water)."),
      sliderInput("danmgkg", "Bolus size (mg/kg)", 0, 5, 2.5, step = 0.5),
      sliderInput("nbolus", "Number of boluses (q6 min)", 0, 8, 4, step = 1),
      checkboxInput("maint", "1 mg/kg q6h maintenance for 48 h", TRUE),

      h4("Mechanism probes"),
      sliderInput("cao", "Extracellular Ca²⁺ (fraction of normal)", 0, 1, 1, step = 0.1),
      sliderInput("nac", "Antioxidant block of RyR1 oxidation", 0, 1, 0, step = 0.1),
      sliderInput("end", "Simulation length (min)", 120, 900, 480, step = 60)
    ),

    mainPanel(
      width = 9,
      uiOutput("kpis"),
      tabsetPanel(
        type = "tabs",

        tabPanel("1 · Patient & event",
          br(), plotOutput("p_overview", height = "560px"),
          div(class = "note", textOutput("txt_overview"))),

        tabPanel("2 · Trigger PK",
          br(),
          div(class = "thesis",
              "MAC is a CNS potency scale. RyR1 reads molar concentration, and ",
              "1 MAC spans 0.75% (halothane) to 6.0% (desflurane) — an 8-fold ",
              "spread. The well-perfused muscle compartment is the trigger ",
              "compartment; its ~2 h time constant is why onset is delayed and ",
              "variable, and why the weakest trigger gives the latest, most ",
              "easily-missed presentation."),
          plotOutput("p_pk", height = "440px"),
          h4("Circuit washout"), tableOutput("t_washout"),
          div(class = "note",
              "Flushing the machine cannot get the circuit below what the ",
              "patient is exhaling. Charcoal filters scrub the exhaled agent ",
              "too, which is why they have an asymptote that the flush does not.")),

        tabPanel("3 · The two stores",
          br(),
          div(class = "thesis",
              "The point of this tab: both curves are driven by the SAME event. ",
              "Their separation is entirely a store-size effect."),
          plotOutput("p_stores", height = "420px"),
          tableOutput("t_capacitance"),
          div(class = "warn", htmlOutput("txt_capacitance"))),

        tabPanel("4 · Ca²⁺ and energetics",
          br(), plotOutput("p_ca", height = "560px"),
          div(class = "note",
              "At steady state SR release and SERCA uptake cancel, so the ",
              "sustained myoplasmic Ca²⁺ level is set by the sarcolemmal ",
              "balance (store-operated entry in, PMCA/NCX out). Drag the ",
              "extracellular-Ca²⁺ slider to see the plateau collapse. ",
              "Note also that phosphocreatine is spent before ATP falls at all ",
              "— and that once ATP does fall, the resulting rigor costs no ",
              "further ATP and cannot be reversed.")),

        tabPanel("5 · Clinical endpoints",
          br(), plotOutput("p_endpoints", height = "560px"),
          h4("MHAUS Clinical Grading Scale (computed as an output, not an input)"),
          tableOutput("t_cgs")),

        tabPanel("6 · Masking",
          br(),
          div(class = "warn",
              strong("The uncomfortable result. "),
              "Hyperventilation acts on the CO", tags$sub("2"), " elimination term ",
              "and cooling acts on the heat-loss term. Neither touches RyR1 open ",
              "probability. Each therefore normalises the monitor it acts on while ",
              "the event continues underneath."),
          plotOutput("p_masking", height = "480px"),
          tableOutput("t_masking")),

        tabPanel("7 · Scenario comparison",
          br(),
          checkboxGroupInput("scensel", "Scenarios", choices = names(scen),
                             selected = names(scen)[c(2, 3, 6, 7, 10, 11)],
                             inline = TRUE),
          plotOutput("p_scen", height = "480px"),
          tableOutput("t_scen")),

        tabPanel("8 · Time to dantrolene",
          br(),
          div(class = "thesis",
              "Larach 2010 (n = 286 North American MH Registry events) found that ",
              "each 30 min of delay multiplied the odds of a complication by 1.61. ",
              "The model reproduces the direction and the monotonicity but is ",
              "steeper than the registry — see the README for why that ",
              "discrepancy is informative rather than fatal."),
          plotOutput("p_delay", height = "420px"), tableOutput("t_delay")),

        tabPanel("9 · Contracture test",
          br(),
          div(class = "thesis",
              "The IVCT/CHCT measures exactly the parameter the genotype moves ",
              "(EC50 for channel opening), which is why a functional assay still ",
              "outperforms sequencing. European MH Group: MHS if ≥ 0.2 g of ",
              "threshold contracture at ≤ 2% halothane or ≤ 2 mM caffeine."),
          plotOutput("p_ivct", height = "380px"),
          tableOutput("t_ivct"),
          h4("Extracellular-Ca²⁺ dependence"), tableOutput("t_ivct_ca")),

        tabPanel("10 · Biomarkers & organs",
          br(), plotOutput("p_biomarker", height = "560px"),
          div(class = "note",
              "Myoglobin has a 2.5 h plasma half-life and CK a 36 h one: ",
              "myoglobin is early and brief, so a normal value at 12 h excludes ",
              "nothing. Destroying 1% of muscle mass moves serum K⁺ by about ",
              "2.4 mmol/L before any acidotic transcellular shift is counted.")),

        tabPanel("11 · Model card",
          br(), htmlOutput("modelcard"))
      )
    )
  )
)

# =====================================================================
#  SERVER
# =====================================================================
server <- function(input, output, session) {

  sim <- reactive({
    simulate_mh(input$genotype, input$agent, input$mac, input$sux, input$wt,
                input$trec, input$dandelay, input$danmgkg, input$nbolus,
                input$maint, input$stopagent, input$charcoal, input$vent,
                input$coolW, input$surfW, input$cao, input$nac, input$end)
  })

  # ---------------- KPI row -------------------------------------------
  output$kpis <- renderUI({
    o <- sim(); tl <- first_lethal(o)
    i <- if (is.na(tl)) nrow(o) else which(o$LETHAL > 0)[1]
    w <- seq_len(i)
    div(class = "kpirow",
      kpi_box("peak EtCO₂", sprintf("%.0f", max(o$ETCO2[w])), "mmHg", "#1565C0"),
      kpi_box("peak core temp", sprintf("%.1f", max(o$TCORE[w])), "°C", "#C62828"),
      kpi_box("peak CK", format(round(o$CK[i]), big.mark = ","), "U/L", "#EF6C00"),
      kpi_box("peak K⁺", sprintf("%.1f", max(o$KTOT[w])), "mmol/L", "#6A1B9A"),
      kpi_box("lowest pH", sprintf("%.2f", min(o$PHART[w])), "arterial", "#00838F"),
      kpi_box("muscle destroyed", sprintf("%.2f%%", o$DESTRPC[i]), "of mass", "#AD1457"),
      kpi_box("CGS rank", sprintf("%.0f", max(o$cgsRank[w])), "1-6 (Larach)", "#455A64"),
      kpi_box("outcome", if (is.na(tl)) "survives" else sprintf("%.0f min", tl),
              if (is.na(tl)) "within the run" else "to lethal physiology",
              if (is.na(tl)) "#2E7D32" else "#B71C1C")
    )
  })

  # ---------------- 1 overview -----------------------------------------
  output$p_overview <- renderPlot({
    long_plot(sim(),
      c("ETCO2","TCORE","HR","PHART","LACS","CK","KTOT","FORCE","Po"),
      c("EtCO₂ (mmHg)","core temp (°C)","heart rate (/min)",
        "arterial pH","lactate (mmol/L)","CK (U/L)","K⁺ (mmol/L)",
        "muscle activation","RyR1 open probability"),
      "The event", vline = input$trec)
  })
  output$txt_overview <- renderText({
    o <- sim(); tl <- first_lethal(o)
    sprintf(paste("Red dashed line = the moment the event is recognised (%d min).",
                  "RyR1 open probability peaks at %.5f (resting %.5f).",
                  "Redox sensitisation of the channel reaches %.2f, and it is that",
                  "term — not the agent — that keeps the event running after",
                  "the vaporiser is closed. %s"),
            input$trec, max(o$Po), min(o$Po), max(o$SENS),
            if (is.na(tl)) "The patient survives this run."
            else sprintf("Lethal physiology at %.0f min.", tl))
  })

  # ---------------- 2 PK -------------------------------------------------
  output$p_pk <- renderPlot({
    long_plot(sim(), c("CIRC","ALV","MUSF","MUSS","DANE","SCHA"),
      c("circuit agent (vol%)","alveolar agent (vol%)",
        "muscle, well-perfused (vol%) — TRIGGER","muscle, bulk (vol%)",
        "dantrolene effect site (mg/L)","succinylcholine (mg)"),
      "Trigger and drug pharmacokinetics", vline = input$trec)
  })
  output$t_washout <- renderTable({ washout() }, digits = 1)

  # ---------------- 3 the two stores -------------------------------------
  output$p_stores <- renderPlot({
    o <- sim()
    b <- o %>% transmute(time,
                `EtCO2 (fraction of its own excursion)` =
                  (ETCO2 - first(ETCO2)) / max(1e-9, max(ETCO2) - first(ETCO2)),
                `core temp (fraction of its own excursion)` =
                  (TCORE - first(TCORE)) / max(1e-9, max(TCORE) - first(TCORE))) %>%
      pivot_longer(-time, names_to = "v", values_to = "y")
    ggplot(b, aes(time, y, colour = v)) + geom_line(linewidth = 1.1) +
      scale_colour_manual(values = c("#1565C0", "#C62828"), name = NULL) +
      geom_vline(xintercept = input$trec, linetype = 2, colour = "#B71C1C") +
      labs(title = "The same event, normalised, read through two stores",
           subtitle = paste("CO2 store 60 mL/mmHg vs core heat store 158 kJ/degC.",
                            "Nothing else differs."),
           x = "time (min)", y = "fraction of own excursion") + THEME
  })
  output$t_capacitance <- renderTable({
    ce <- capacitance_experiment()
    tibble(quantity = names(ce), value = unlist(ce))
  }, digits = 2)
  output$txt_capacitance <- renderUI({
    ce <- capacitance_experiment()
    HTML(sprintf(paste0(
      "<b>%.0f mmHg of EtCO<sub>2</sub> per °C.</b> One degree of fever costs ",
      "the same elapsed time as %.0f mmHg of EtCO<sub>2</sub>. Since the ",
      "ventilation-limited EtCO<sub>2</sub> plateau is only about %.0f mmHg, the ",
      "entire capnograph scale is exhausted before the thermometer moves a ",
      "single degree. Waiting for fever costs roughly %.0f min of unopposed ",
      "rhabdomyolysis."),
      ce$`mmHg of EtCO2 per degC`, ce$`mmHg of EtCO2 per degC`,
      ce$`EtCO2 plateau (mmHg)`, ce$`core 38.5 degC (min)`))
  })

  # ---------------- 4 calcium / energetics -------------------------------
  output$p_ca <- renderPlot({
    long_plot(sim(),
      c("CAM","CASR","CAMITO","MITOEFF","ATP","PCR","ATPuse","Q_MUS","SENS"),
      c("myoplasmic Ca²⁺ (µM) — THE DRIVER","SR Ca²⁺ (µM)",
        "mitochondrial Ca²⁺ (µM)","mitochondrial coupling efficiency",
        "ATP (mM)","phosphocreatine (mM)","ATP turnover (µM/min)",
        "muscle perfusion (L/min)","RyR1 oxidation (SENS)"),
      "Calcium handling and energetics", vline = input$trec)
  })

  # ---------------- 5 endpoints -------------------------------------------
  output$p_endpoints <- renderPlot({
    long_plot(sim(),
      c("ETCO2","TCORE","cgs","DESTRPC","SCR","GFRF","FIB","PLT","NEURO"),
      c("EtCO₂ (mmHg)","core temp (°C)","CGS raw score",
        "muscle destroyed (%)","creatinine (mg/dL)","GFR (fraction)",
        "fibrinogen (mg/dL)","platelets (10⁹/L)","cerebral injury index"),
      "Clinical endpoints", vline = input$trec)
  })
  output$t_cgs <- renderTable({
    o <- sim()
    tibble(`peak raw score` = max(o$cgs), `peak rank` = max(o$cgsRank),
           interpretation = c("1 almost never / 2 unlikely / 3 somewhat less than likely /
                               4 somewhat greater than likely / 5 very likely /
                               6 almost certain")) },
    digits = 0)

  # ---------------- 6 masking ---------------------------------------------
  masking <- reactive({
    base <- function(...) simulate_mh(input$genotype, input$agent, input$mac,
                                      input$sux, input$wt, input$trec, ...,
                                      ca_o = 1, nac = 0, end = input$end)
    list(
      untreated = base(0, 0, 0, FALSE, FALSE, FALSE, 1, 0, 0),
      vent_only = base(0, 0, 0, FALSE, FALSE, FALSE, 3, 0, 0),
      cool_only = base(0, 0, 0, FALSE, FALSE, FALSE, 1, 200, 100),
      full      = base(2, 2.5, 4, TRUE, TRUE, TRUE, 2, 200, 100))
  })
  output$p_masking <- renderPlot({
    m <- masking()
    d <- bind_rows(lapply(names(m), function(n)
      m[[n]] %>% transmute(time, arm = n, EtCO2 = ETCO2, core = TCORE, CK = CK))) %>%
      pivot_longer(c(EtCO2, core, CK), names_to = "v", values_to = "y")
    ggplot(d, aes(time, y, colour = arm)) + geom_line(linewidth = 1) +
      facet_wrap(~v, scales = "free_y") +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(title = "Two arms produce a normal monitor and a dead patient",
           subtitle = paste("hyperventilation normalises EtCO2; cooling normalises",
                            "temperature; neither touches RyR1"),
           x = "time (min)", y = NULL) + THEME
  })
  output$t_masking <- renderTable({
    m <- masking()
    do.call(rbind, lapply(names(m), function(n) {
      o <- m[[n]]; tl <- first_lethal(o)
      i <- if (is.na(tl)) nrow(o) else which(o$LETHAL > 0)[1]; w <- seq_len(i)
      tibble(arm = n, `peak EtCO2` = max(o$ETCO2[w]), `peak core` = max(o$TCORE[w]),
             CK = o$CK[i], `peak K+` = max(o$KTOT[w]),
             outcome = if (is.na(tl)) "survives" else sprintf("dies at %.0f min", tl))
    }))
  }, digits = 1)

  # ---------------- 7 scenarios --------------------------------------------
  scen_res <- reactive({
    req(input$scensel)
    do.call(rbind, lapply(input$scensel, summarise_scenario, end = input$end))
  })
  output$p_scen <- renderPlot({
    req(input$scensel)
    d <- do.call(rbind, lapply(input$scensel, run_scenario, end = input$end)) %>%
      select(time, scenario, ETCO2, TCORE, CK, KTOT) %>%
      pivot_longer(c(ETCO2, TCORE, CK, KTOT), names_to = "v", values_to = "y")
    ggplot(d, aes(time, y, colour = scenario)) + geom_line(linewidth = 0.9) +
      facet_wrap(~v, scales = "free_y") +
      scale_colour_manual(values = rep(PAL, 4), name = NULL) +
      labs(title = "Scenario comparison", x = "time (min)", y = NULL) + THEME
  })
  output$t_scen <- renderTable({ scen_res() }, digits = 2)

  # ---------------- 8 delay --------------------------------------------------
  dc <- reactive(delay_curve())
  output$p_delay <- renderPlot({
    d <- dc() %>% pivot_longer(-delay_min, names_to = "v", values_to = "y")
    ggplot(d, aes(delay_min, y, colour = v)) +
      geom_line(linewidth = 1) + geom_point() +
      facet_wrap(~v, scales = "free_y") +
      scale_colour_manual(values = rep(PAL, 3), guide = "none") +
      labs(title = "The price of delay",
           x = "minutes from recognition to first dantrolene dose", y = NULL) + THEME
  })
  output$t_delay <- renderTable({ dc() }, digits = 2)

  # ---------------- 9 contracture test ---------------------------------------
  output$p_ivct <- renderPlot({
    hal <- seq(0, 3, 0.1)
    d <- do.call(rbind, lapply(names(GENO), function(g) {
      b <- ivct(g)
      tibble(genotype = g, halothane = hal,
             g_force = sapply(hal, function(h) ivct(g, halothane = h) - b))
    }))
    ggplot(d, aes(halothane, g_force, colour = genotype)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0.2, linetype = 2, colour = "#B71C1C") +
      geom_vline(xintercept = 2.0, linetype = 3, colour = "#455A64") +
      scale_colour_manual(values = PAL, name = NULL) +
      labs(title = "Simulated halothane contracture test",
           subtitle = "dashed = 0.2 g EMHG threshold; dotted = 2% halothane",
           x = "bath halothane (vol%)", y = "threshold contracture (g)") + THEME
  })
  output$t_ivct    <- renderTable({ ivct_panel() }, digits = 2)
  output$t_ivct_ca <- renderTable({ ivct_ca_dependence() }, digits = 2)

  # ---------------- 10 biomarkers ----------------------------------------------
  output$p_biomarker <- renderPlot({
    long_plot(sim(),
      c("CK","MB","KTOT","PHOS","LACS","BE","MBT","DDIM","TDOSE"),
      c("CK (U/L)","myoglobin (ng/mL)","K⁺ (mmol/L)","phosphate (mmol/L)",
        "lactate (mmol/L)","base excess (mEq/L)","tubular cast burden",
        "D-dimer (ng/mL FEU)","thermal dose (°C·min > 40)"),
      "Biomarkers and organ injury", vline = input$trec)
  })

  # ---------------- 11 model card -----------------------------------------------
  output$modelcard <- renderUI(HTML('
  <h4>What this model is</h4>
  <p>A 39-compartment ODE model of malignant hyperthermia in a 70 kg adult,
  spanning volatile-agent pharmacokinetics, RyR1 gating, sarcoplasmic and
  sarcolemmal calcium handling, muscle energetics, whole-body CO<sub>2</sub> and
  heat balance, rhabdomyolysis and its systemic consequences, and dantrolene
  PK/PD. Time unit: minutes.</p>

  <h4>What is fitted</h4>
  <ol>
    <li>EC50 and resting open probability per genotype, set so the simulated
        contracture test classifies MHN as MHN and MHS as MHS at European MH
        Group thresholds.</li>
    <li>Store-operated Ca<sup>2+</sup> entry gain per genotype, set so sustained
        myoplasmic Ca<sup>2+</sup> reaches about 1.05 &micro;M.</li>
    <li>The muscle injury gain, set against observed CK ranges.</li>
    <li>Dantrolene PK, from Flewellen 1983.</li>
    <li>Relative RyR1 potency per vol% of the five volatile agents, set so that
        onset ranks halothane &lt; enflurane &lt; isoflurane &lt; sevoflurane
        &lt; desflurane at 1 MAC.</li>
  </ol>

  <h4>What is derived and therefore falsifiable</h4>
  <ul>
    <li>The 111 mmHg-per-&deg;C capacitance ratio and the EtCO<sub>2</sub>-first
        ordering of signs.</li>
    <li>The failure of agent withdrawal alone.</li>
    <li>Both masking arms.</li>
    <li>RQ rising above 1 (bicarbonate buffering of lactate liberates
        CO<sub>2</sub> that never came from O<sub>2</sub>).</li>
    <li>The extracellular-Ca<sup>2+</sup> dependence of the contracture.</li>
    <li>Non-depolarising blockers having exactly zero effect on rigidity.</li>
  </ul>

  <h4>Where the model disagrees with the textbook, and says so</h4>
  <ul>
    <li>A core temperature rise of 1&ndash;2 &deg;C every 5 min requires
        500&ndash;1000 W of net heat storage. At the model&rsquo;s aerobic
        ceiling that is reachable only transiently, only in the core
        compartment, and only at near-total muscle recruitment. A fulminant
        event runs at 3&ndash;5 &deg;C/h here.</li>
    <li>Succinylcholine advances sustained onset by only about 3 min, which is
        weaker than its epidemiological association with fulminant MH. This is
        the model&rsquo;s least well-supported structural choice.</li>
    <li>The complication rate rises more steeply with time-to-dantrolene than
        Larach 2010&rsquo;s odds ratio of 1.61 per 30 min.</li>
  </ul>

  <h4>Not for clinical use</h4>
  <p>Educational and research model only. It has not been validated against
  patient-level data and must not be used for any clinical decision. In a real
  suspected MH event, follow the MHAUS / EMHG protocol and telephone your
  national MH hotline.</p>'))
}

shinyApp(ui, server)
