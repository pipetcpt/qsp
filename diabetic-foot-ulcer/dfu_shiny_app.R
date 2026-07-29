## ============================================================================
##  Diabetic Foot Ulcer (DFU) QSP — Shiny Dashboard
##  ---------------------------------------------------------------------------
##  A 12-tab interactive front end to dfu_mrgsolve_model.R.
##
##  The layout deliberately mirrors the six structural commitments of the model,
##  so that each tab answers one clinical question that the model can actually
##  settle:
##
##     Tab 1  Patient        who is this foot? (neuropathy, perfusion, deformity)
##     Tab 2  Offloading     is the device on the foot? (efficacy x adherence)
##     Tab 3  Wound geometry the perimeter law and the PAR4 decision rule
##     Tab 4  Oxygen gate    is anything anabolic even possible here?
##     Tab 5  Inflammation   the M1->M2 switch that never happens
##     Tab 6  Proteases      MMP-9/TIMP-1 and what it does to applied drug
##     Tab 7  Repair         angiogenesis, fibroblasts, senescence, granulation
##     Tab 8  Infection      bioburden, biofilm tolerance, osteomyelitis
##     Tab 9  Therapy PK     every course, its window and its tissue exposure
##     Tab 10 Endpoints      PAR4, closure, time-to-heal, amputation
##     Tab 11 Remission      the part after closure, where the disease lives
##     Tab 12 Compare        scenario library, side by side
##
##  Run with:   shiny::runApp("dfu_shiny_app.R")
##  Requires:   shiny, mrgsolve, dplyr, tidyr, ggplot2, DT, purrr
## ============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)
library(purrr)

mod <- mread_cache("dfu", "dfu_mrgsolve_model.R")

## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------
debride_ev <- function(period = 14, until = 140, amt = 1) {
  if (is.null(period) || period <= 0) return(NULL)
  ev(time = seq(0, until, by = period), amt = amt, cmt = 11)
}

run_one <- function(pars, deb_period = 14, end = 540, delta = 0.5) {
  m <- mod %>% param(pars)
  e <- debride_ev(deb_period)
  out <- if (is.null(e)) mrgsim(m, end = end, delta = delta)
         else            mrgsim(m, events = e, end = end, delta = delta)
  as.data.frame(out)
}

## Time to closure, PAR4 and the closed-form perimeter-law prediction
endpoints_of <- function(d, area0, close = 0.05) {
  a28 <- d$AREA[which.min(abs(d$time - 28))]
  par4 <- 1 - a28 / area0
  tclose <- if (any(d$AREA <= close)) d$time[which.max(d$AREA <= close)] else NA_real_
  list(
    par4      = par4,
    a28       = a28,
    a84       = d$AREA[which.min(abs(d$time - 84))],
    tclose    = tclose,
    tpred     = if (par4 > 0 && par4 < 1) 28 / (1 - sqrt(1 - par4)) else NA_real_,
    closed12  = !is.na(tclose) && tclose <= 84,
    rec12     = if (!is.na(tclose) && tclose + 365 <= max(d$time)) {
                  h2 <- d$P_REC[which.min(abs(d$time - (tclose + 365)))]
                  h1 <- d$P_REC[which.min(abs(d$time - tclose))]
                  h2 - h1
                } else NA_real_,
    amp365    = d$P_AMP[which.min(abs(d$time - 365))]
  )
}

theme_dfu <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")
}

OFFLOAD_DEVICES <- list(
  "None"                            = c(eff = 0.00, adh = 1.00),
  "Half-shoe / therapeutic shoe"    = c(eff = 0.40, adh = 0.55),
  "Removable cast walker (RCW)"     = c(eff = 0.87, adh = 0.28),
  "Instant TCC (rendered irremov.)" = c(eff = 0.87, adh = 1.00),
  "Total contact cast (TCC)"        = c(eff = 0.85, adh = 1.00)
)

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Diabetic Foot Ulcer — Quantitative Systems Pharmacology Dashboard"),
  p(em(paste("Research / education only. Closure is modelled as a perimeter",
             "process; offloading as device efficacy x adherence; topical growth",
             "factor as a race against MMP-9."))),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Patient"),
      sliderInput("AREA0", "Wound area at presentation (cm²)",
                  0.2, 12, 2.0, step = 0.1),
      sliderInput("DEPTH0", "Wound depth (mm)", 0.5, 10, 3.0, step = 0.5),
      sliderInput("HBA1C0", "HbA1c at presentation (%)", 6, 14, 9.0, step = 0.1),
      sliderInput("PERF0", "Effective wound perfusion (0–1)",
                  0.15, 1.0, 0.85, step = 0.01),
      helpText("≈0.85 neuropathic, ≈0.55 neuro-ischaemic, ≤0.35 CLTI"),
      sliderInput("DEFORMITY", "Foot deformity (0–1)", 0, 1, 0.45, step = 0.05),
      sliderInput("BACT0", "Bioburden at presentation (log₁₀ CFU/g)",
                  3, 9, 5.2, step = 0.1),
      checkboxInput("PROBE_BONE", "Probe-to-bone positive", FALSE),

      hr(), h4("Offloading"),
      selectInput("device", "Device", names(OFFLOAD_DEVICES),
                  selected = "Removable cast walker (RCW)"),
      sliderInput("ADHERENCE", "Adherence — fraction of daily activity worn",
                  0, 1, 0.28, step = 0.01),
      helpText(strong("Armstrong 2003:"), "removable walkers are worn during",
               "28% of daily activity. Drag this to 1.0 to see the whole",
               "TCC-vs-walker gap disappear."),

      hr(), h4("Wound bed"),
      selectInput("deb", "Sharp debridement",
                  c("Weekly" = 7, "Every 2 weeks" = 14, "Every 4 weeks" = 28,
                    "None" = 0), selected = 14),

      hr(), h4("Topical / device therapy"),
      checkboxInput("bec",  "Becaplermin 0.01% gel (rhPDGF-BB)", FALSE),
      checkboxInput("nosf", "TLC-NOSF sucrose octasulfate dressing", FALSE),
      checkboxInput("esm",  "Esmolol 14% topical gel", FALSE),
      checkboxInput("ctp",  "Cellular / tissue product", FALSE),
      checkboxInput("npwt", "Negative-pressure wound therapy", FALSE),
      checkboxInput("hbot", "Hyperbaric oxygen (40 sessions)", FALSE),

      hr(), h4("Infection / vascular"),
      checkboxInput("abx",  "Systemic antibiotic", FALSE),
      sliderInput("T_ABX_END", "Antibiotic course (days)", 0, 84, 14, step = 7),
      checkboxInput("rif",  "Rifampicin adjunct (bone / biofilm)", FALSE),
      checkboxInput("revasc", "Revascularise", FALSE),
      sliderInput("T_REVASC", "Revascularisation day", 0, 90, 14, step = 1),

      hr(), h4("Systemic / remission care"),
      sliderInput("HBA1C_TGT", "HbA1c target (%)", 6, 14, 9.0, step = 0.1),
      checkboxInput("FOOTWEAR", "Therapeutic footwear after closure", FALSE),
      checkboxInput("TEMPMON",  "Home skin-temperature monitoring", FALSE),
      checkboxInput("EDU",      "Structured education / self-inspection", FALSE),
      checkboxInput("SURVEIL",  "Podiatric surveillance", FALSE),

      hr(),
      sliderInput("end", "Simulation horizon (days)", 120, 1825, 540, step = 30),
      actionButton("go", "Simulate", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        ## ---- 1 ------------------------------------------------------------
        tabPanel(
          "1 · Patient profile",
          h4("The three axes that make a foot ulcerate"),
          plotOutput("p_profile", height = 300),
          hr(),
          h4("Slow axes over the full horizon"),
          plotOutput("p_slow", height = 300),
          helpText(paste("Neuropathy relaxes with a ~4.6-year time constant at",
                         "HbA1c 7.0. Whatever you do to glycaemia today changes",
                         "the NEXT ulcer, not this one."))
        ),

        ## ---- 2 ------------------------------------------------------------
        tabPanel(
          "2 · Offloading",
          h4("Effective offloading = device efficacy × adherence"),
          fluidRow(
            column(6, plotOutput("p_offload_bar", height = 320)),
            column(6, plotOutput("p_offload_area", height = 320))
          ),
          hr(),
          h4("Adherence sweep for the current device"),
          plotOutput("p_adh_sweep", height = 300),
          DTOutput("t_offload")
        ),

        ## ---- 3 ------------------------------------------------------------
        tabPanel(
          "3 · Wound geometry & the perimeter law",
          h4("Area, radius and perimeter"),
          plotOutput("p_geom", height = 320),
          helpText(paste("Under dA/dt = −k·P the RADIUS falls linearly. A",
                         "straight line in the middle panel means the model is",
                         "obeying the perimeter law and k is constant.")),
          hr(),
          h4("PAR₄ and what it predicts"),
          fluidRow(
            column(6, plotOutput("p_par4", height = 300)),
            column(6, DTOutput("t_par4"))
          ),
          helpText(HTML(paste0(
            "Given PAR<sub>4</sub>, the perimeter law gives a size-free forecast ",
            "<b>t<sub>heal</sub> = 28 / (1 − √(1 − PAR<sub>4</sub>))</b>. ",
            "It over-predicts whenever k is still rising during weeks 1–4, ",
            "which is why the 50% rule is a conservative screen.")))
        ),

        ## ---- 4 ------------------------------------------------------------
        tabPanel(
          "4 · Oxygen gate",
          h4("Wound pO₂ and the single Hill gate that multiplies every anabolic term"),
          plotOutput("p_o2", height = 340),
          hr(),
          fluidRow(
            column(6, plotOutput("p_hif", height = 300)),
            column(6, helpText(HTML(paste0(
              "<b>The HIF paradox.</b> The diabetic wound is hypoxic, so the ",
              "hypoxia signal is high — but methylglyoxal adducts on HIF-1α/p300 ",
              "block transactivation, so the response the wound is entitled to ",
              "never arrives. <br><br>",
              "<b>The gate.</b> OG = pO₂³/(pO₂³ + 26³). At TcPO₂ 20 mmHg the gate ",
              "is ~0.3; at 50 mmHg it is ~0.86. Below the gate, collagen ",
              "hydroxylation, angiogenesis, the M1→M2 switch and the neutrophil ",
              "oxidative burst all fail together — which is why no biologic ",
              "works on an ischaemic wound and why revascularisation comes first."))))
          )
        ),

        ## ---- 5 ------------------------------------------------------------
        tabPanel(
          "5 · Inflammation",
          h4("The M1 → M2 switch that never happens"),
          plotOutput("p_inflam", height = 340),
          hr(),
          h4("IL-1β self-amplification"),
          plotOutput("p_il1", height = 280),
          helpText(paste("IL-1β sustains M1 polarisation, which sustains IL-1β.",
                         "The switch to a reparative macrophage phenotype is",
                         "gated by IL-1β, oxygen and hyperglycaemia together —",
                         "fixing only one of the three does not release it."))
        ),

        ## ---- 6 ------------------------------------------------------------
        tabPanel(
          "6 · Proteases & drug survival",
          h4("MMP-9 / TIMP-1 — the healing predictor"),
          plotOutput("p_prot", height = 320),
          hr(),
          h4("What the protease environment does to applied becaplermin"),
          plotOutput("p_bec", height = 320),
          helpText(HTML(paste0(
            "Becaplermin is cleared at (K<sub>out</sub> + K<sub>deg</sub>·MMP9/TIMP1). ",
            "Turn on the TLC-NOSF dressing and watch wound becaplermin rise ",
            "without changing the dose. At the receptor the interaction is ",
            "superadditive (~1.45× the sum of the single-agent gains in PDGF ",
            "drive); at the whole wound it comes out almost exactly ",
            "multiplicative on healing rate, because MMP-9 excess damages the ",
            "wound through several parallel channels the dressing fixes anyway.")))
        ),

        ## ---- 7 ------------------------------------------------------------
        tabPanel(
          "7 · Repair tissue",
          h4("Angiogenesis"),
          plotOutput("p_angio", height = 300),
          hr(),
          h4("Fibroblasts, senescence and matrix"),
          plotOutput("p_fib", height = 300),
          hr(),
          h4("Granulation fill, depth and keratinocyte migratory competence"),
          plotOutput("p_gran", height = 300)
        ),

        ## ---- 8 ------------------------------------------------------------
        tabPanel(
          "8 · Infection & biofilm",
          h4("Bioburden, biofilm and the fraction of the antibiotic that reaches it"),
          plotOutput("p_inf", height = 340),
          helpText(paste("The lower panel is the percentage of the antibiotic's",
                         "potential kill that survives biofilm tolerance. Sharp",
                         "debridement opens a window; the biofilm closes it",
                         "within 2–3 days.")),
          hr(),
          h4("Osteomyelitis and the exposed-bone nidus"),
          plotOutput("p_osteo", height = 300),
          helpText(paste("Antibiotic reduces the bone burden; only resection",
                         "removes the substrate the regrowth term feeds on."))
        ),

        ## ---- 9 ------------------------------------------------------------
        tabPanel(
          "9 · Therapy PK",
          h4("Topical and device exposures"),
          plotOutput("p_pk_topical", height = 320),
          hr(),
          h4("Systemic antibiotic: central, wound tissue, bone"),
          plotOutput("p_pk_abx", height = 320),
          hr(),
          h4("Cumulative treatment burden and antibiotic-days"),
          plotOutput("p_burden", height = 250)
        ),

        ## ---- 10 -----------------------------------------------------------
        tabPanel(
          "10 · Clinical endpoints",
          fluidRow(
            column(4, wellPanel(h4("PAR₄"), h2(textOutput("e_par4")))),
            column(4, wellPanel(h4("Time to closure"), h2(textOutput("e_tclose")))),
            column(4, wellPanel(h4("Closed by 12 weeks"), h2(textOutput("e_c12"))))
          ),
          fluidRow(
            column(4, wellPanel(h4("Recurrence @ 12 mo"), h2(textOutput("e_rec")))),
            column(4, wellPanel(h4("Amputation risk @ 1 y"), h2(textOutput("e_amp")))),
            column(4, wellPanel(h4("Antibiotic-days"), h2(textOutput("e_abxd"))))
          ),
          hr(),
          plotOutput("p_endpoints", height = 340),
          DTOutput("t_endpoints")
        ),

        ## ---- 11 -----------------------------------------------------------
        tabPanel(
          "11 · Remission (after closure)",
          h4("Closure is remission, not cure"),
          plotOutput("p_remission", height = 340),
          helpText(HTML(paste0(
            "Closing the wound resets AREA. It does <b>not</b> reset neuropathy, ",
            "deformity or plantar pressure, and it leaves scar at ~80% of normal ",
            "tensile strength. Everything on this tab happens at the point where ",
            "most wound-care trials stop measuring."))),
          hr(),
          h4("What each element of remission care is worth"),
          plotOutput("p_remcare", height = 320),
          DTOutput("t_remcare")
        ),

        ## ---- 12 -----------------------------------------------------------
        tabPanel(
          "12 · Scenario comparison",
          h4("Prebuilt scenario library"),
          checkboxGroupInput(
            "scen", NULL, inline = TRUE,
            choices = c("no offloading", "standard care (RCW 28%)",
                        "total contact cast", "RCW at 100% adherence",
                        "+ becaplermin", "+ TLC-NOSF", "+ both",
                        "ischaemic", "revascularised d14",
                        "optimal bundle", "optimal bundle + remission care"),
            selected = c("standard care (RCW 28%)", "total contact cast",
                         "RCW at 100% adherence", "optimal bundle + remission care")
          ),
          actionButton("go2", "Run comparison", class = "btn-primary"),
          hr(),
          plotOutput("p_compare", height = 400),
          DTOutput("t_compare")
        )
      )
    )
  )
)

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  ## keep the adherence slider in step with the chosen device
  observeEvent(input$device, {
    updateSliderInput(session, "ADHERENCE",
                      value = unname(OFFLOAD_DEVICES[[input$device]]["adh"]))
  })

  current_pars <- reactive({
    dev <- OFFLOAD_DEVICES[[input$device]]
    list(
      AREA0 = input$AREA0, DEPTH0 = input$DEPTH0,
      HBA1C0 = input$HBA1C0, HBA1C_TGT = input$HBA1C_TGT,
      PERF0 = input$PERF0, PERF_TGT = input$PERF0,
      DEFORMITY = input$DEFORMITY, BACT0 = input$BACT0,
      PROBE_BONE = as.numeric(input$PROBE_BONE),
      DEV_EFF = unname(dev["eff"]), ADHERENCE = input$ADHERENCE,
      R_BEC = as.numeric(input$bec), R_NOSF = as.numeric(input$nosf),
      R_ESM = as.numeric(input$esm), R_CTP = as.numeric(input$ctp),
      NPWT = as.numeric(input$npwt), R_OXY = as.numeric(input$hbot),
      R_ABX = as.numeric(input$abx), T_ABX_END = input$T_ABX_END,
      RIF_ON = as.numeric(input$rif),
      T_REVASC = if (input$revasc) input$T_REVASC else 1e6,
      REVASC_GAIN = if (input$revasc) 0.45 else 0,
      FOOTWEAR = as.numeric(input$FOOTWEAR), TEMPMON = as.numeric(input$TEMPMON),
      EDU = as.numeric(input$EDU), SURVEIL = as.numeric(input$SURVEIL)
    )
  })

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    run_one(current_pars(), deb_period = as.numeric(input$deb), end = input$end)
  })

  ep <- reactive(endpoints_of(sim(), input$AREA0))

  long <- function(d, vars) {
    d %>% select(time, all_of(vars)) %>% pivot_longer(-time)
  }

  ## ---- 1 patient ----------------------------------------------------------
  output$p_profile <- renderPlot({
    long(sim(), c("NEURO", "PERF", "STRESS_o", "TCPO2")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#2c6fbb") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL) + theme_dfu()
  })
  output$p_slow <- renderPlot({
    d <- run_one(current_pars(), as.numeric(input$deb), end = max(input$end, 1825))
    long(d, c("HBA1C", "NEURO")) %>%
      ggplot(aes(time / 365, value)) + geom_line(linewidth = 1, colour = "#8a3ffc") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "years", y = NULL) + theme_dfu()
  })

  ## ---- 2 offloading -------------------------------------------------------
  offload_grid <- reactive({
    imap_dfr(OFFLOAD_DEVICES, function(v, nm) {
      p <- modifyList(current_pars(),
                      list(DEV_EFF = unname(v["eff"]), ADHERENCE = unname(v["adh"])))
      d <- run_one(p, as.numeric(input$deb), end = input$end)
      e <- endpoints_of(d, input$AREA0)
      tibble(device = nm, eff = v["eff"], adh = v["adh"],
             effective = v["eff"] * v["adh"],
             t_close = e$tclose, par4 = 100 * e$par4)
    })
  })

  output$p_offload_bar <- renderPlot({
    offload_grid() %>%
      pivot_longer(c(eff, adh, effective)) %>%
      mutate(name = factor(name, c("eff", "adh", "effective"),
                           c("device efficacy", "adherence", "EFFECTIVE"))) %>%
      ggplot(aes(reorder(device, value), value, fill = name)) +
      geom_col(position = "dodge") + coord_flip() +
      scale_fill_manual(values = c("#9ec5e8", "#f5b971", "#2c6fbb")) +
      labs(x = NULL, y = "fraction", fill = NULL) + theme_dfu()
  })

  output$p_offload_area <- renderPlot({
    imap_dfr(OFFLOAD_DEVICES, function(v, nm) {
      p <- modifyList(current_pars(),
                      list(DEV_EFF = unname(v["eff"]), ADHERENCE = unname(v["adh"])))
      run_one(p, as.numeric(input$deb), end = min(input$end, 200)) %>%
        transmute(time, AREA, device = nm)
    }) %>%
      ggplot(aes(time, AREA, colour = device)) + geom_line(linewidth = 0.9) +
      labs(x = "day", y = "wound area (cm²)", colour = NULL) + theme_dfu()
  })

  output$p_adh_sweep <- renderPlot({
    dev <- OFFLOAD_DEVICES[[input$device]]
    map_dfr(seq(0, 1, by = 0.1), function(a) {
      p <- modifyList(current_pars(),
                      list(DEV_EFF = unname(dev["eff"]), ADHERENCE = a))
      d <- run_one(p, as.numeric(input$deb), end = input$end)
      tibble(adherence = a, t_close = endpoints_of(d, input$AREA0)$tclose)
    }) %>%
      ggplot(aes(adherence, t_close)) +
      geom_line(linewidth = 1, colour = "#d1495b") + geom_point(size = 2) +
      geom_vline(xintercept = 0.28, linetype = 2) +
      annotate("text", x = 0.28, y = Inf, vjust = 1.5, hjust = -0.05,
               label = "observed RCW wear (Armstrong 2003)", size = 3.2) +
      labs(x = "adherence (fraction of daily activity worn)",
           y = "time to closure (days)") + theme_dfu()
  })

  output$t_offload <- renderDT(datatable(offload_grid(), options = list(dom = "t")))

  ## ---- 3 geometry ---------------------------------------------------------
  output$p_geom <- renderPlot({
    long(sim(), c("AREA", "RADIUS", "PERIM_o")) %>%
      mutate(name = factor(name, c("AREA", "RADIUS", "PERIM_o"),
                           c("area (cm²)", "radius (cm) — should be LINEAR",
                             "perimeter (cm)"))) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#2c6fbb") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL) + theme_dfu()
  })

  par4_grid <- reactive({
    map_dfr(c(0.5, 1, 2, 4, 8), function(a0) {
      d <- run_one(modifyList(current_pars(), list(AREA0 = a0)),
                   as.numeric(input$deb), end = input$end)
      e <- endpoints_of(d, a0)
      tibble(A0 = a0, r0 = sqrt(a0 / pi), PAR4 = 100 * e$par4,
             t_close = e$tclose, t_pred_PAR4 = e$tpred)
    })
  })
  output$p_par4 <- renderPlot({
    par4_grid() %>%
      ggplot(aes(PAR4, t_close)) +
      geom_line(colour = "#2c6fbb", linewidth = 1) +
      geom_point(aes(size = A0)) +
      geom_line(aes(y = t_pred_PAR4), linetype = 2, colour = "#d1495b") +
      geom_vline(xintercept = 50, linetype = 3) +
      labs(x = "PAR₄ (%)", y = "time to closure (days)", size = "A₀ (cm²)",
           caption = "dashed red = 28/(1−√(1−PAR₄)); dotted = Sheehan 50% rule") +
      theme_dfu()
  })
  output$t_par4 <- renderDT(datatable(par4_grid(), options = list(dom = "t")) %>%
                              formatRound(2:5, 2))

  ## ---- 4 oxygen -----------------------------------------------------------
  output$p_o2 <- renderPlot({
    long(sim(), c("TCPO2", "OG_o", "PERF", "VASC")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#0f766e") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL) + theme_dfu()
  })
  output$p_hif <- renderPlot({
    long(sim(), c("VEGF", "SDF1", "EPC")) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "day", y = "normalised", colour = NULL) + theme_dfu()
  })

  ## ---- 5 inflammation -----------------------------------------------------
  output$p_inflam <- renderPlot({
    long(sim(), c("NEUT", "M1", "M2", "ROS")) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "day", y = "normalised wound units", colour = NULL) + theme_dfu()
  })
  output$p_il1 <- renderPlot({
    sim() %>% mutate(`M2/M1` = M2 / pmax(M1, 1e-6)) %>%
      long(c("IL1B", "M2/M1")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#b45309") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL) + theme_dfu()
  })

  ## ---- 6 proteases --------------------------------------------------------
  output$p_prot <- renderPlot({
    long(sim(), c("MMP9", "TIMP1", "PROT_o")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#7c2d12") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL) + theme_dfu()
  })
  output$p_bec <- renderPlot({
    bind_rows(
      run_one(modifyList(current_pars(), list(R_BEC = 1, R_NOSF = 0)),
              as.numeric(input$deb), end = input$end) %>%
        transmute(time, BEC_W, PDGF_TOT_o, arm = "becaplermin alone"),
      run_one(modifyList(current_pars(), list(R_BEC = 1, R_NOSF = 1)),
              as.numeric(input$deb), end = input$end) %>%
        transmute(time, BEC_W, PDGF_TOT_o, arm = "becaplermin + TLC-NOSF")
    ) %>% pivot_longer(c(BEC_W, PDGF_TOT_o)) %>%
      ggplot(aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c("#d1495b", "#2c6fbb")) +
      labs(x = "day", y = NULL, colour = NULL) + theme_dfu()
  })

  ## ---- 7 repair -----------------------------------------------------------
  output$p_angio <- renderPlot({
    long(sim(), c("VEGF", "EPC", "VASC")) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "day", y = NULL, colour = NULL) + theme_dfu()
  })
  output$p_fib <- renderPlot({
    long(sim(), c("FIB", "SEN", "COL")) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "day", y = NULL, colour = NULL) + theme_dfu()
  })
  output$p_gran <- renderPlot({
    long(sim(), c("GRAN", "KERA", "DEPTH")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#166534") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL) + theme_dfu()
  })

  ## ---- 8 infection --------------------------------------------------------
  output$p_inf <- renderPlot({
    long(sim(), c("BACT", "BIOF", "BIOF_TOL")) %>%
      mutate(name = factor(name, c("BACT", "BIOF", "BIOF_TOL"),
                           c("bioburden (log₁₀ CFU/g)", "biofilm biomass",
                             "% of antibiotic effect surviving tolerance"))) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#9a3412") +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "day", y = NULL) + theme_dfu()
  })
  output$p_osteo <- renderPlot({
    long(sim(), c("OSTEO", "BEXP", "ABX_B")) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "day", y = NULL, colour = NULL) + theme_dfu()
  })

  ## ---- 9 PK ---------------------------------------------------------------
  output$p_pk_topical <- renderPlot({
    long(sim(), c("BEC_W", "NOSF", "ESM_W", "CTP", "OXY")) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "day", y = "au", colour = NULL) + theme_dfu()
  })
  output$p_pk_abx <- renderPlot({
    long(sim(), c("ABX_C", "ABX_W", "ABX_B")) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "day", y = "au", colour = NULL) + theme_dfu()
  })
  output$p_burden <- renderPlot({
    long(sim(), c("TOXIDX", "ABXD")) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#525252") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL) + theme_dfu()
  })

  ## ---- 10 endpoints -------------------------------------------------------
  output$e_par4   <- renderText(sprintf("%.1f%%", 100 * ep()$par4))
  output$e_tclose <- renderText(ifelse(is.na(ep()$tclose), "not closed",
                                       sprintf("%.0f d", ep()$tclose)))
  output$e_c12    <- renderText(ifelse(ep()$closed12, "yes", "no"))
  output$e_rec    <- renderText(ifelse(is.na(ep()$rec12), "—",
                                       sprintf("%.1f%%", ep()$rec12)))
  output$e_amp    <- renderText(sprintf("%.1f%%", ep()$amp365))
  output$e_abxd   <- renderText(sprintf("%.0f d", max(sim()$ABXD)))

  output$p_endpoints <- renderPlot({
    sim() %>%
      ggplot(aes(time, AREA)) +
      geom_line(linewidth = 1, colour = "#2c6fbb") +
      geom_hline(yintercept = 0.05, linetype = 2) +
      geom_vline(xintercept = c(28, 84), linetype = 3) +
      annotate("text", x = 28, y = Inf, vjust = 1.4, label = "PAR₄", size = 3.5) +
      annotate("text", x = 84, y = Inf, vjust = 1.4, label = "12 wk", size = 3.5) +
      labs(x = "day", y = "wound area (cm²)") + theme_dfu()
  })
  output$t_endpoints <- renderDT({
    e <- ep()
    datatable(tibble(
      metric = c("PAR₄ (%)", "area at 12 wk (cm²)", "time to closure (d)",
                 "perimeter-law forecast from PAR₄ (d)",
                 "recurrence @ 12 mo (%)", "amputation @ 1 y (%)"),
      value = c(100 * e$par4, e$a84, e$tclose, e$tpred, e$rec12, e$amp365)
    ), options = list(dom = "t")) %>% formatRound("value", 2)
  })

  ## ---- 11 remission -------------------------------------------------------
  output$p_remission <- renderPlot({
    long(sim(), c("SCAR", "NEURO", "STRESS_o", "P_REC")) %>%
      mutate(name = factor(name, c("SCAR", "NEURO", "STRESS_o", "P_REC"),
                           c("scar maturation", "neuropathy (unchanged)",
                             "mechanical stress", "cumulative re-ulceration (%)"))) %>%
      ggplot(aes(time, value)) + geom_line(linewidth = 0.9, colour = "#7e22ce") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL) + theme_dfu()
  })

  remcare <- reactive({
    combos <- list(
      "nothing"                     = list(),
      "+ therapeutic footwear"      = list(FOOTWEAR = 1),
      "+ temperature monitoring"    = list(TEMPMON = 1),
      "+ education + surveillance"  = list(EDU = 1, SURVEIL = 1),
      "+ HbA1c 7.0"                 = list(HBA1C_TGT = 7),
      "everything"                  = list(FOOTWEAR = 1, TEMPMON = 1, EDU = 1,
                                           SURVEIL = 1, HBA1C_TGT = 7)
    )
    imap_dfr(combos, function(extra, nm) {
      d <- run_one(modifyList(current_pars(), extra),
                   as.numeric(input$deb), end = max(input$end, 540))
      e <- endpoints_of(d, input$AREA0)
      tibble(aftercare = nm, t_close = e$tclose, recurrence_12mo = e$rec12)
    })
  })
  output$p_remcare <- renderPlot({
    remcare() %>%
      pivot_longer(c(t_close, recurrence_12mo)) %>%
      ggplot(aes(reorder(aftercare, value), value, fill = name)) +
      geom_col(show.legend = FALSE) + coord_flip() +
      facet_wrap(~name, scales = "free_x") +
      scale_fill_manual(values = c("#7e22ce", "#2c6fbb")) +
      labs(x = NULL, y = NULL) + theme_dfu()
  })
  output$t_remcare <- renderDT(datatable(remcare(), options = list(dom = "t")) %>%
                                 formatRound(2:3, 1))

  ## ---- 12 comparison ------------------------------------------------------
  SCEN <- list(
    "no offloading"                 = list(DEV_EFF = 0.00, ADHERENCE = 1.00),
    "standard care (RCW 28%)"       = list(DEV_EFF = 0.87, ADHERENCE = 0.28),
    "total contact cast"            = list(DEV_EFF = 0.85, ADHERENCE = 1.00),
    "RCW at 100% adherence"         = list(DEV_EFF = 0.87, ADHERENCE = 1.00),
    "+ becaplermin"                 = list(DEV_EFF = 0.87, ADHERENCE = 0.28, R_BEC = 1),
    "+ TLC-NOSF"                    = list(DEV_EFF = 0.87, ADHERENCE = 0.28, R_NOSF = 1),
    "+ both"                        = list(DEV_EFF = 0.87, ADHERENCE = 0.28,
                                           R_BEC = 1, R_NOSF = 1),
    "ischaemic"                     = list(DEV_EFF = 0.87, ADHERENCE = 0.28,
                                           PERF0 = 0.32, PERF_TGT = 0.32),
    "revascularised d14"            = list(DEV_EFF = 0.87, ADHERENCE = 0.28,
                                           PERF0 = 0.32, PERF_TGT = 0.32,
                                           T_REVASC = 14, REVASC_GAIN = 0.45),
    "optimal bundle"                = list(DEV_EFF = 0.85, ADHERENCE = 1.00,
                                           R_BEC = 1, R_NOSF = 1),
    "optimal bundle + remission care" = list(DEV_EFF = 0.85, ADHERENCE = 1.00,
                                             R_BEC = 1, R_NOSF = 1, FOOTWEAR = 1,
                                             TEMPMON = 1, EDU = 1, SURVEIL = 1,
                                             HBA1C_TGT = 7)
  )

  comp <- eventReactive(input$go2, {
    req(length(input$scen) > 0)
    imap_dfr(SCEN[input$scen], function(extra, nm) {
      per <- if (grepl("optimal", nm)) 7 else as.numeric(input$deb)
      run_one(modifyList(current_pars(), extra), per, end = max(input$end, 540)) %>%
        mutate(scenario = nm)
    })
  })

  output$p_compare <- renderPlot({
    comp() %>% filter(time <= 200) %>%
      ggplot(aes(time, AREA, colour = scenario)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 0.05, linetype = 2) +
      labs(x = "day", y = "wound area (cm²)", colour = NULL) + theme_dfu()
  })
  output$t_compare <- renderDT({
    comp() %>% group_by(scenario) %>% group_modify(~{
      e <- endpoints_of(.x, input$AREA0)
      tibble(PAR4 = 100 * e$par4, t_close = e$tclose,
             closed_12wk = e$closed12, recurrence_12mo = e$rec12,
             amputation_1y = e$amp365)
    }) %>% ungroup() %>% datatable(options = list(dom = "t")) %>%
      formatRound(c(2, 3, 5, 6), 1)
  })
}

shinyApp(ui, server)
