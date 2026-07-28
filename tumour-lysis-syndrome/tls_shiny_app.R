##############################################################################
## tls_shiny_app.R
## Tumour Lysis Syndrome QSP — interactive dashboard
## ===========================================================================
## Front end for tls_mrgsolve_model.R.  The tabs are organised around the
## STRUCTURE of the model rather than around organ systems, because the whole
## point is that the disease is a race between a release flux and a clearance
## capacity, closed into a loop through the kidney:
##
##   1  Patient & tumour        who is being treated and how big the load is
##   2  The race                UA_req vs UA_crit — the two threshold curves
##   3  Drug PK                 every exposure the model tracks
##   4  Solute trajectories     urate, K, phosphate, calcium, magnesium
##   5  Kidney & crystals       supersaturation, deposition, GFR, urine flow
##   6  Operator comparison     which term of the loop each therapy touches
##   7  Urine pH optimum        the alkalinisation trade, both sides priced
##   8  Potassium rescue        redistribution vs clearance, total body K
##   9  Clinical endpoints      Cairo-Bishop, hazards, trial ledger
##  10  Parameters & provenance where every number came from
##
## Run:
##   shiny::runApp("tls_shiny_app.R")
##############################################################################

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

source("tls_mrgsolve_model.R")

PAL <- c(urate = "#7e57c2", K = "#00897b", PO4 = "#1e88e5", Ca = "#f9a825",
         Cr = "#c62828", eGFR = "#455a64", crystal = "#5d4037",
         cap = "#2e7d32", xan = "#8bc34a")

theme_tls <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey35", size = 10),
          legend.position = "bottom", legend.title = element_blank(),
          strip.text = element_text(face = "bold"))
}

##############################################################################
## UI
##############################################################################

ui <- fluidPage(
  titlePanel("Tumour Lysis Syndrome — QSP dashboard"),
  tags$p(style = "color:#555;margin-top:-8px",
         tags$em(paste("Release flux vs clearance capacity, closed into a loop",
                       "through the kidney. 48 ODEs. Educational model —",
                       "not for clinical use."))),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient and tumour"),
      sliderInput("N0", "Tumour burden (×10¹² cells)",
                  min = 0.05, max = 10, value = 3, step = 0.05),
      helpText(textOutput("burden_note", inline = TRUE)),
      selectInput("td", "Tumour doubling time",
                  c("Burkitt / B-ALL (30 h)" = 30,
                    "aggressive lymphoma (72 h)" = 72,
                    "AML (200 h)" = 200,
                    "CLL (2000 h)" = 2000), selected = 30),
      sliderInput("neph0", "Starting nephron mass (prior CKD)",
                  min = 0.25, max = 1.0, value = 1.0, step = 0.05),
      checkboxInput("g6pd", "G6PD deficient", FALSE),

      h4("Cytoreduction"),
      radioButtons("trigger", NULL,
                   c("Cytotoxic chemotherapy at t = 0" = "chemo",
                     "Venetoclax (choose a ramp)" = "ven",
                     "Spontaneous lysis only" = "none")),
      conditionalPanel(
        "input.trigger == 'ven'",
        selectInput("ramp", "Venetoclax schedule",
                    c("400 mg from day 1", "200 mg from day 1",
                      "100 mg from day 1", "50 mg from day 1",
                      "2-step 20/50", "5-week label ramp",
                      "8-week slow ramp"),
                    selected = "5-week label ramp")),
      numericInput("pred_days", "Steroid prephase (days before t = 0; 0 = none)",
                   0, min = 0, max = 10, step = 1),

      h4("Prophylaxis and support"),
      selectInput("hydration", "Fluid",
                  c("2 L/day (standard)" = HYD_STD,
                    "3 L/m²/day (guideline high risk)" = HYD_AGGR,
                    "4 L/m²/day (maximal)" = HYD_MAX),
                  selected = HYD_AGGR),
      selectInput("xou", "Urate-directed therapy",
                  c("none", "allopurinol 300 mg", "allopurinol 600 mg",
                    "febuxostat 120 mg", "rasburicase 0.20 mg/kg",
                    "rasburicase 0.20 mg/kg × 5",
                    "rasburicase 0.40 mg/kg × 5"),
                  selected = "rasburicase 0.20 mg/kg"),
      numericInput("lead", "Lead time before t = 0 (h)", 0,
                   min = -0, max = 168, step = 12),
      sliderInput("hco3", "Bicarbonate infusion (mmol/h)",
                  min = 0, max = 200, value = 0, step = 5),
      helpText(textOutput("ph_note", inline = TRUE)),
      checkboxInput("seve", "Phosphate binder (sevelamer q8h)", FALSE),
      checkboxInput("szc", "Potassium binder (SZC 10 g q8h)", FALSE),
      checkboxInput("furo", "Furosemide 40 mg q8h", FALSE),
      sliderInput("ca_inf", "IV calcium (mmol/h)",
                  min = 0, max = 6, value = 0, step = 0.5),
      checkboxInput("hd", "Dialysis, day 2–5", FALSE),
      numericInput("tend", "Simulation horizon (h)", 336,
                   min = 96, max = 1200, step = 24),
      hr(),
      actionButton("go", "Simulate", class = "btn-primary", width = "100%")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        ## ------------------------------------------------------------------
        tabPanel(
          "1 · Patient & tumour",
          br(),
          fluidRow(column(6, plotOutput("p_tumour", height = 300)),
                   column(6, plotOutput("p_release", height = 300))),
          h4("Released load against the pools it lands in"),
          tableOutput("t_load"),
          tags$small(paste(
            "Intracellular contents per 10¹² cells: K⁺ 38.5 mmol (110 mmol/L",
            "cell volume × 0.35 pL/cell), phosphate 60 mmol (~48 of it",
            "nucleic-acid phosphate at 16 pg nucleic acid per cell, MW 330),",
            "purine 24 mmol. The extracellular pools those land in are 69 mmol",
            "of potassium and 19 mmol of phosphate, which is the whole",
            "problem."))
        ),
        ## ------------------------------------------------------------------
        tabPanel(
          "2 · The race",
          br(),
          p(strong("UA_req"), " is the plasma urate at which renal excretion",
            "balances the peak purine release flux. ", strong("UA_crit"),
            " is the plasma urate at which tubular fluid reaches the",
            "metastable limit and urate begins to precipitate. The disease is",
            "the region where UA_req exceeds UA_crit: the kidney cannot",
            "excrete the load without crystallising it."),
          plotOutput("p_race", height = 380),
          fluidRow(column(6, h4("Thresholds for this prescription"),
                          tableOutput("t_thresholds")),
                   column(6, h4("The same race for the other two solutes"),
                          tableOutput("t_hierarchy"))),
          tags$small(paste(
            "Note that the Cairo–Bishop laboratory-TLS urate cut-off of",
            "8.0 mg/dL falls inside the UA_crit range this model computes from",
            "solubility and urine flow alone (6.8 mg/dL at 2 L/day, 10.5 at",
            "4 L/m²/day). The criterion is reproduced as a crystallisation",
            "threshold rather than fitted to."))
        ),
        ## ------------------------------------------------------------------
        tabPanel(
          "3 · Drug PK",
          br(), plotOutput("p_pk", height = 460),
          h4("Exposure summary"), tableOutput("t_pk"),
          tags$small(paste(
            "Oxypurinol clearance scales with GFR, so the drug given to",
            "prevent the kidney injury accumulates because of it. Febuxostat",
            "is hepatically cleared and does not."))
        ),
        ## ------------------------------------------------------------------
        tabPanel(
          "4 · Solutes",
          br(), plotOutput("p_solutes", height = 520),
          fluidRow(column(6, h4("Peaks and nadirs"), tableOutput("t_peaks")),
                   column(6, h4("Purine cascade"), plotOutput("p_purine",
                                                              height = 260)))
        ),
        ## ------------------------------------------------------------------
        tabPanel(
          "5 · Kidney & crystals",
          br(), plotOutput("p_kidney", height = 480),
          fluidRow(column(6, plotOutput("p_ss", height = 300)),
                   column(6, plotOutput("p_injury", height = 300))),
          tags$small(paste(
            "Crystal growth is bounded by delivery: deposition can never",
            "exceed 0.95 × the tubular delivery rate of the solute. Tubular",
            "obstruction also limits urine flow independently of GFR, which is",
            "what makes crystal nephropathy oliguric rather than merely",
            "azotaemic."))
        ),
        ## ------------------------------------------------------------------
        tabPanel(
          "6 · Operator comparison",
          br(),
          p("Every therapy is classified by which term of the loop it",
            "touches. The model's central claim is that the classes are not",
            "interchangeable — no amount of a flux operator empties a pool,",
            "and no amount of a sequestration operator removes phosphate that",
            "came from inside the patient."),
          plotOutput("p_operators", height = 400),
          tableOutput("t_operators"),
          tags$small(paste(
            "POOL = rasburicase, dialysis · FLUX = allopurinol, febuxostat ·",
            "FLUX-SHAPING = venetoclax ramp, steroid prephase · DILUTION =",
            "hydration · SPECIATION = bicarbonate · REDISTRIBUTION =",
            "insulin/glucose, β₂-agonist · SEQUESTRATION = sevelamer, SZC."))
        ),
        ## ------------------------------------------------------------------
        tabPanel(
          "7 · Urine pH optimum",
          br(),
          p("Uric acid solubility rises with pH (pKa 5.75). So does the",
            "calcium-phosphate driving force, because the precipitating",
            "species is HPO₄²⁻ (pKa₂ 6.80), and systemic alkalosis lowers",
            "ionised calcium as well. Both sides are priced here, and the",
            "optimum moves with whichever solute is rate-limiting."),
          fluidRow(column(6, plotOutput("p_ph_curves", height = 320)),
                   column(6, plotOutput("p_ph_outcome", height = 320))),
          tableOutput("t_ph"),
          tags$small(paste(
            "In the reference integration the optimum urine pH is 6.6 at a",
            "3×10¹²-cell burden without rasburicase, 7.0 at 6×10¹² without",
            "rasburicase, and 5.9 — i.e. no alkali at all — once rasburicase",
            "is on board. Routine alkalinisation was abandoned; the model",
            "derives why rather than asserting it."))
        ),
        ## ------------------------------------------------------------------
        tabPanel(
          "8 · Potassium rescue",
          br(),
          p("Insulin/glucose and β₂-agonists move potassium between",
            "compartments. Dialysis, gut binders and (in the volume-replete",
            "patient) diuretics remove it from the body. Total body",
            "exchangeable potassium is plotted alongside the serum",
            "concentration so the difference is visible rather than asserted."),
          plotOutput("p_krescue", height = 400),
          tableOutput("t_krescue")
        ),
        ## ------------------------------------------------------------------
        tabPanel(
          "9 · Clinical endpoints",
          br(),
          fluidRow(column(6, h4("Cairo–Bishop classification"),
                          tableOutput("t_cairo")),
                   column(6, h4("Cumulative hazards"),
                          plotOutput("p_hazard", height = 280))),
          h4("Safety limb: the oxidant load of urate oxidase"),
          plotOutput("p_oxidant", height = 260),
          tags$small(paste(
            "Urate oxidase makes one H₂O₂ per urate destroyed, so the patient",
            "with the largest tumour receives the largest peroxide load: the",
            "toxicity is proportional to the reason for giving the drug. G6PD",
            "deficiency is a contraindication for exactly that reason."))
        ),
        ## ------------------------------------------------------------------
        tabPanel(
          "10 · Parameters & provenance",
          br(),
          h4("Where the numbers came from"),
          tags$ul(
            tags$li(strong("Cell contents: "),
                    "K⁺ 110 mmol/L of cell volume at 0.35 pL per blast;",
                    "nucleic acid 16 pg/cell at MW 330 gives 48 mmol of",
                    "nucleotide and therefore ~24 mmol of purine and ~48 mmol",
                    "of nucleic-acid phosphate per 10¹² cells."),
            tags$li(strong("Urate solubility: "),
                    "S = S_HU·(1 + 10^(pH − 5.75)) with S_HU 0.655 mmol/L,",
                    "giving 130 mg/L at pH 5.0 against a reported ~150."),
            tags$li(strong("Urate handling: "),
                    "basal FE 8% and 700 mg/day excretion; FE rises to ~30% as",
                    "URAT1/GLUT9 reabsorption saturates."),
            tags$li(strong("Rasburicase: "),
                    "V1 8 L, CL 0.30 L/h (t½ 18 h); Km for urate ~25 µmol/L,",
                    "two orders of magnitude below TLS urate, so the reaction",
                    "is zero order and a dose buys a fixed mmol/h."),
            tags$li(strong("Allopurinol/oxypurinol: "),
                    "allopurinol t½ 1.4 h, oxypurinol t½ 23 h with renal",
                    "clearance; XO IC50 1.1 mg/L for oxypurinol."),
            tags$li(strong("Phosphate: "),
                    "TmP/GFR ≈ 0.90 mmol/L with FGF23-like and PTH-driven",
                    "suppression; basal excretion 23 mmol/day."),
            tags$li(strong("Calcium phosphate: "),
                    "systemic threshold set at Ca × PO₄ = 4.84 mmol²/L², which",
                    "is the classical 60 mg²/dL².")),
          h4("What is calibrated rather than derived"),
          tags$ul(
            tags$li("The crystal nucleation rate constants, the obstruction",
                    "half-max crystal mass, and the nephron loss/recovery",
                    "rates are calibrated so that an unprophylaxed",
                    "high-burden Burkitt patient reaches a peak creatinine",
                    "ratio of about 3 with recovery over two weeks. They are",
                    "not measured quantities."),
            tags$li("The interstitial calcium-phosphate deposition rate",
                    "constant is calibrated to give hypocalcaemia at the",
                    "phosphate concentrations at which it is clinically",
                    "observed."),
            tags$li("Hazard functions for arrhythmia, seizure and dialysis are",
                    "shapes chosen to be flat at baseline and steep in the",
                    "pathological range; the absolute probabilities should not",
                    "be quoted.")),
          h4("Full parameter block"),
          tableOutput("t_params")
        )
      )
    )
  )
)

##############################################################################
## SERVER
##############################################################################

server <- function(input, output, session) {

  ramp_table <- function(nm) {
    switch(nm,
      "400 mg from day 1" = data.frame(time = 0, dose = 400),
      "200 mg from day 1" = data.frame(time = 0, dose = 200),
      "100 mg from day 1" = data.frame(time = 0, dose = 100),
      "50 mg from day 1"  = data.frame(time = 0, dose = 50),
      "2-step 20/50"      = data.frame(time = c(0, 168), dose = c(20, 50)),
      "5-week label ramp" = data.frame(time = c(0, 168, 336, 504, 672),
                                       dose = c(20, 50, 100, 200, 400)),
      "8-week slow ramp"  = data.frame(time = c(0, 168, 336, 504, 672, 840),
                                       dose = c(10, 20, 50, 100, 200, 400)))
  }

  args <- reactive({
    a <- list(N0 = input$N0, hydration = as.numeric(input$hydration),
              tend = input$tend, hco3 = input$hco3, ca_inf = input$ca_inf,
              neph0 = input$neph0,
              extra_param = list(TD = as.numeric(input$td),
                                 G6PD = ifelse(input$g6pd, 0, 1)))
    lead <- -abs(input$lead)
    a$chemo <- if (input$trigger == "chemo") c(0, 24, 1.0) else NULL
    if (input$trigger == "ven") a$ven <- ramp_table(input$ramp)
    if (input$pred_days > 0) {
      a$pred_start <- -24 * input$pred_days
      a$pred_days <- input$pred_days
    }
    switch(input$xou,
      "allopurinol 300 mg" = { a$allo_start <- lead },
      "allopurinol 600 mg" = { a$allo_start <- lead; a$allo_dose <- 600 },
      "febuxostat 120 mg"  = { a$febu_start <- lead },
      "rasburicase 0.20 mg/kg" = { a$rasb_start <- 0 },
      "rasburicase 0.20 mg/kg × 5" = { a$rasb_start <- 0; a$rasb_days <- 5 },
      "rasburicase 0.40 mg/kg × 5" = { a$rasb_start <- 0; a$rasb_days <- 5
                                       a$rasb_mgkg <- 0.40 })
    if (input$seve) a$seve_start <- 0
    if (input$szc)  a$szc_start  <- 0
    if (input$furo) a$furo_start <- 0
    if (input$hd)   a$dialysis   <- c(48, 120, 1.2)
    a
  })

  sim <- eventReactive(input$go, { do.call(tls_sim, args()) },
                       ignoreNULL = FALSE)
  smry <- reactive(tls_summary(sim()))

  ## ---------------- helper notes ----------------
  output$burden_note <- renderText({
    j <- peak_release_flux(input$N0, td = as.numeric(input$td))
    sprintf("≈ %.1f g of releasable uric acid; peak purine flux %.2f mmol/h",
            input$N0 * pget("Q_PUR") * 0.168, j)
  })
  output$ph_note <- renderText({
    ph <- pget("PH_BASE") + pget("PH_RISE") *
      hillf(input$hco3 / pget("KEL_HCO3"), pget("KH_HCO3"))
    sprintf("steady-state urine pH ≈ %.2f", ph)
  })

  ## ---------------- tab 1 ----------------
  output$p_tumour <- renderPlot({
    d <- sim()
    ggplot(d, aes(time / 24)) +
      geom_line(aes(y = N_VIA, colour = "viable"), linewidth = 1) +
      geom_line(aes(y = N_LYS, colour = "committed to lysis"), linewidth = 1) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
      scale_colour_manual(values = c(viable = "#37474f",
                                     `committed to lysis` = "#d84315")) +
      labs(title = "Tumour", subtitle = "t = 0 is the start of cytoreduction",
           x = "days", y = "×10¹² cells") + theme_tls()
  })

  output$p_release <- renderPlot({
    d <- sim()
    rel <- pget("KLYS") * d$N_LYS
    df <- data.frame(t = d$time / 24,
                     urate = rel * pget("Q_PUR"),
                     potassium = rel * pget("Q_K"),
                     phosphate = rel * pget("Q_PO4")) %>%
      pivot_longer(-t)
    caps <- data.frame(
      name = c("urate", "potassium", "phosphate"),
      cap = c(fe_ua_of(0.476) * pget("GFR0") * 0.476, k_capacity(6.0), 12.0))
    ggplot(df, aes(t, value, colour = name)) +
      geom_line(linewidth = 1) +
      geom_hline(data = caps, aes(yintercept = cap, colour = name),
                 linetype = 2) +
      scale_colour_manual(values = c(urate = PAL[["urate"]],
                                     potassium = PAL[["K"]],
                                     phosphate = PAL[["PO4"]])) +
      labs(title = "Release flux against clearance capacity",
           subtitle = "dashed = the capacity the intact kidney has for that solute",
           x = "days", y = "mmol/h") + theme_tls()
  })

  output$t_load <- renderTable({
    N0 <- input$N0
    data.frame(
      solute = c("potassium", "phosphate", "purine → urate"),
      `released (mmol)` = c(pget("Q_K"), pget("Q_PO4"), pget("Q_PUR")) * N0,
      `extracellular pool (mmol)` = c(4.10, 1.15, 0.32) * pget("VEC0"),
      `ratio` = (c(pget("Q_K"), pget("Q_PO4"), pget("Q_PUR")) * N0) /
        (c(4.10, 1.15, 0.32) * pget("VEC0")),
      check.names = FALSE)
  }, digits = 1)

  ## ---------------- tab 2 ----------------
  output$p_race <- renderPlot({
    qu <- as.numeric(input$hydration) - pget("INSENS")
    ph <- pget("PH_BASE") + pget("PH_RISE") *
      hillf(input$hco3 / pget("KEL_HCO3"), pget("KH_HCO3"))
    uc <- ua_critical(qu, ph) * 16.81
    b <- seq(0.05, 10, length.out = 60)
    ur <- sapply(b, function(n)
      ua_required(peak_release_flux(n, td = as.numeric(input$td))) * 16.81)
    df <- data.frame(burden = b, UA_req = ur, UA_crit = uc)
    cross <- if (any(ur > uc)) b[which(ur > uc)[1]] else NA
    ggplot(df, aes(burden)) +
      geom_ribbon(aes(ymin = UA_crit, ymax = pmax(UA_req, UA_crit)),
                  fill = "#7e57c2", alpha = 0.15) +
      geom_line(aes(y = UA_req, colour = "UA_req (needed to excrete)"),
                linewidth = 1.1) +
      geom_line(aes(y = UA_crit, colour = "UA_crit (precipitation begins)"),
                linewidth = 1.1) +
      geom_hline(yintercept = 8, linetype = 3) +
      annotate("text", x = max(b) * 0.75, y = 8.8, size = 3.2,
               label = "Cairo–Bishop urate criterion, 8.0 mg/dL") +
      { if (!is.na(cross)) geom_vline(xintercept = cross, linetype = 2,
                                      colour = "#4527a0") } +
      geom_point(data = data.frame(x = input$N0,
                                   y = ua_required(peak_release_flux(
                                     input$N0, td = as.numeric(input$td))) * 16.81),
                 aes(x, y), size = 3, colour = "#4527a0") +
      scale_colour_manual(values = c("#c62828", "#1e88e5")) +
      scale_y_log10() +
      labs(title = "The race, as two threshold concentrations",
           subtitle = paste0("shaded = the region where the load cannot be ",
                             "excreted without crystallising; crossing at ",
                             ifelse(is.na(cross), "—",
                                    sprintf("%.2f×10¹² cells", cross))),
           x = "tumour burden (×10¹² cells)", y = "plasma urate (mg/dL, log)") +
      theme_tls()
  })

  output$t_thresholds <- renderTable({
    do.call(rbind, lapply(list(c(HYD_STD, 5.90, "2 L/day, pH 5.9"),
                               c(HYD_AGGR, 5.90, "3 L/m²/day, pH 5.9"),
                               c(HYD_MAX, 5.90, "4 L/m²/day, pH 5.9"),
                               c(HYD_AGGR, 6.50, "3 L/m²/day, pH 6.5"),
                               c(HYD_AGGR, 7.00, "3 L/m²/day, pH 7.0")),
      function(r) data.frame(
        prescription = r[[3]],
        `UA_crit (mg/dL)` = ua_critical(as.numeric(r[[1]]) - pget("INSENS"),
                                        as.numeric(r[[2]])) * 16.81,
        check.names = FALSE)))
  }, digits = 1)

  output$t_hierarchy <- renderTable({
    h <- data.frame(
      solute = c("urate", "potassium", "phosphate"),
      `peak flux (mmol/h)` = c(
        peak_release_flux(input$N0, q = pget("Q_PUR"), td = as.numeric(input$td)),
        peak_release_flux(input$N0, q = pget("Q_K"), td = as.numeric(input$td)),
        peak_release_flux(input$N0, q = pget("Q_PO4"), td = as.numeric(input$td))),
      `capacity (mmol/h)` = c(fe_ua_of(0.476) * pget("GFR0") * 0.476,
                              k_capacity(6.0), 12.0),
      check.names = FALSE)
    h$`flux / capacity` <- h[[2]] / h[[3]]
    h
  }, digits = 2)

  ## ---------------- tab 3 ----------------
  output$p_pk <- renderPlot({
    d <- sim()
    df <- data.frame(t = d$time,
                     `oxypurinol (mg/L)` = d$C_OXY,
                     `febuxostat (mg/L)` = d$C_FEBU,
                     `rasburicase (mg/L)` = d$C_RASB,
                     `venetoclax (mg/L)` = d$C_VEN,
                     `XO inhibition (fraction)` = d$XOI,
                     `urine pH` = d$urinepH,
                     check.names = FALSE) %>% pivot_longer(-t)
    ggplot(df, aes(t / 24, value)) + geom_line(linewidth = 0.9, colour = "#37474f") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
      labs(title = "Exposures", x = "days", y = NULL) + theme_tls()
  })

  output$t_pk <- renderTable({
    d <- sim()
    data.frame(quantity = c("oxypurinol Cmax", "oxypurinol at day 7",
                            "febuxostat Cmax", "rasburicase Cmax",
                            "peak XO inhibition", "peak urine pH"),
               value = c(max(d$C_OXY), approx(d$time, d$C_OXY, 168)$y,
                         max(d$C_FEBU), max(d$C_RASB), max(d$XOI),
                         max(d$urinepH)))
  }, digits = 3)

  ## ---------------- tab 4 ----------------
  output$p_solutes <- renderPlot({
    d <- sim()
    df <- data.frame(t = d$time,
                     `urate (mg/dL)` = d$UA_mgdl,
                     `potassium (mmol/L)` = d$K,
                     `phosphate (mmol/L)` = d$PO4,
                     `ionised calcium (mmol/L)` = d$CaIon,
                     `magnesium (mmol/L)` = d$Mg,
                     `LDH (U/L)` = d$LDH_UL,
                     check.names = FALSE) %>% pivot_longer(-t)
    thr <- data.frame(name = c("urate (mg/dL)", "potassium (mmol/L)",
                               "phosphate (mmol/L)"),
                      y = c(8.0, 6.0, 1.45))
    ggplot(df, aes(t / 24, value)) +
      geom_line(linewidth = 1, colour = "#1e88e5") +
      geom_hline(data = thr, aes(yintercept = y), linetype = 2,
                 colour = "#c62828") +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "Solute trajectories",
           subtitle = "red dashed = Cairo–Bishop laboratory-TLS thresholds",
           x = "days", y = NULL) + theme_tls()
  })

  output$t_peaks <- renderTable({
    s <- smry()
    data.frame(endpoint = c("peak urate (mg/dL)", "peak K (mmol/L)",
                            "peak phosphate (mmol/L)",
                            "ionised calcium nadir (mmol/L)",
                            "peak creatinine ratio", "eGFR nadir (mL/min)",
                            "peak LDH (U/L)", "tumour lysed (×10¹² cells)"),
               value = c(s$UA_peak, s$K_peak, s$PO4_peak, s$Ca_nadir,
                         s$Cr_ratio, s$eGFR_nadir, s$LDH_peak, s$lysed))
  }, digits = 2)

  output$p_purine <- renderPlot({
    d <- sim()
    df <- data.frame(t = d$time,
                     hypoxanthine = d$HYPOX / pget("VEC0"),
                     xanthine = d$XAN / pget("VEC0"),
                     urate = d$URATE / pget("VEC0"),
                     allantoin = d$ALLANT / pget("VEC0")) %>% pivot_longer(-t)
    ggplot(df, aes(t / 24, value, colour = name)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = pget("S_XAN"), linetype = 3, colour = "#8bc34a") +
      annotate("text", x = Inf, y = pget("S_XAN"), hjust = 1.05, vjust = -0.5,
               size = 3, label = "xanthine solubility") +
      labs(title = "Purine cascade", x = "days", y = "mmol/L") + theme_tls()
  })

  ## ---------------- tab 5 ----------------
  output$p_kidney <- renderPlot({
    d <- sim()
    df <- data.frame(t = d$time,
                     `eGFR (mL/min)` = d$eGFR,
                     `creatinine (mg/dL)` = d$Cr,
                     `urine output (L/day)` = d$QU_Lday,
                     `obstruction (fraction)` = d$OBSTR,
                     `nephron mass` = d$NEPH,
                     `tubular injury` = d$TUBINJ,
                     check.names = FALSE) %>% pivot_longer(-t)
    ggplot(df, aes(t / 24, value)) + geom_line(linewidth = 1, colour = "#455a64") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
      labs(title = "Kidney", x = "days", y = NULL) + theme_tls()
  })

  output$p_ss <- renderPlot({
    d <- sim()
    df <- data.frame(t = d$time, `urate supersaturation` = d$SS_UA,
                     `systemic Ca × PO₄` = d$SS_SYS,
                     check.names = FALSE) %>% pivot_longer(-t)
    ggplot(df, aes(t / 24, value, colour = name)) + geom_line(linewidth = 1) +
      geom_hline(yintercept = pget("SS_TH_UA"), linetype = 2,
                 colour = PAL[["urate"]]) +
      geom_hline(yintercept = 1, linetype = 2, colour = PAL[["cap"]]) +
      scale_colour_manual(values = c(PAL[["cap"]], PAL[["urate"]])) +
      labs(title = "Supersaturation",
           subtitle = "dashed = the nucleation threshold for each species",
           x = "days", y = "SS") + theme_tls()
  })

  output$p_injury <- renderPlot({
    d <- sim()
    df <- data.frame(t = d$time, urate = d$XT_UA, xanthine = d$XT_XAN,
                     `calcium phosphate` = d$XT_CAP,
                     check.names = FALSE) %>% pivot_longer(-t)
    ggplot(df, aes(t / 24, value, fill = name)) +
      geom_area(position = "stack", alpha = 0.85) +
      scale_fill_manual(values = c(`calcium phosphate` = PAL[["cap"]],
                                   urate = PAL[["urate"]],
                                   xanthine = PAL[["xan"]])) +
      labs(title = "Renal crystal burden by species",
           subtitle = "obstruction half-max at 12 mmol",
           x = "days", y = "mmol") + theme_tls()
  })

  ## ---------------- tab 6 ----------------
  op_res <- eventReactive(input$go, {
    N0 <- input$N0
    arms <- list(
      "nothing (2 L/day)"      = list(hydration = HYD_STD),
      "DILUTION 3 L/m²/day"    = list(hydration = HYD_AGGR),
      "DILUTION 4 L/m²/day"    = list(hydration = HYD_MAX),
      "FLUX allopurinol"       = list(allo_start = 0),
      "FLUX febuxostat"        = list(febu_start = 0),
      "POOL rasburicase"       = list(rasb_start = 0),
      "POOL dialysis"          = list(dialysis = c(48, 120, 1.2)),
      "SPECIATION pH 7.5"      = list(hco3 = 45),
      "SEQUESTRATION PO₄ bind" = list(seve_start = 0),
      "FLUX-SHAPE prephase"    = list(pred_start = -120, pred_days = 5))
    do.call(rbind, lapply(names(arms), function(nm) {
      a <- arms[[nm]]; a$N0 <- N0
      a$extra_param <- list(TD = as.numeric(input$td))
      s <- tls_summary(do.call(tls_sim, a), nm)
      s[, c("label", "UA_peak", "K_peak", "PO4_peak", "Ca_nadir", "Cr_ratio",
            "XT_UA", "XT_CAP", "P_rrt")]
    }))
  }, ignoreNULL = FALSE)

  output$p_operators <- renderPlot({
    r <- op_res()
    r$class <- sub(" .*", "", r$label)
    ggplot(r, aes(reorder(label, -Cr_ratio), Cr_ratio, fill = class)) +
      geom_col() + coord_flip() +
      geom_hline(yintercept = 1, linetype = 2) +
      labs(title = "Peak creatinine ratio by operator class",
           subtitle = paste0("tumour burden ", input$N0, "×10¹² cells"),
           x = NULL, y = "peak creatinine / baseline") + theme_tls()
  })
  output$t_operators <- renderTable(op_res(), digits = 3)

  ## ---------------- tab 7 ----------------
  output$p_ph_curves <- renderPlot({
    ph <- seq(5.2, 7.9, by = 0.05)
    f0 <- 1 / (1 + 10^(pget("PKA2_PO4") - pget("PH_BASE")))
    s0 <- pget("S_HU") * (1 + 10^(pget("PH_BASE") - pget("PKA_UA")))
    df <- data.frame(
      pH = ph,
      `urate solubility (× pH 5.9)` =
        pget("S_HU") * (1 + 10^(ph - pget("PKA_UA"))) / s0,
      `HPO₄²⁻ fraction (× pH 5.9)` =
        (1 / (1 + 10^(pget("PKA2_PO4") - ph))) / f0,
      check.names = FALSE) %>% pivot_longer(-pH)
    ggplot(df, aes(pH, value, colour = name)) + geom_line(linewidth = 1.1) +
      scale_y_log10() +
      scale_colour_manual(values = c(PAL[["cap"]], PAL[["urate"]])) +
      labs(title = "Both sides of the pH trade",
           subtitle = "benefit on urate, cost on calcium phosphate",
           x = "urine pH", y = "fold change vs pH 5.9 (log)") + theme_tls()
  })

  ph_res <- eventReactive(input$go, {
    do.call(rbind, lapply(c(0, 16, 40, 90, 200), function(h) {
      a <- args(); a$hco3 <- h
      d <- do.call(tls_sim, a); s <- tls_summary(d)
      data.frame(`urine pH` = tail(d$urinepH, 1), `urate crystal` = s$XT_UA,
                 `CaP crystal` = s$XT_CAP, `ionised Ca nadir` = s$Ca_nadir,
                 `peak Cr ratio` = s$Cr_ratio, `P(RRT)` = s$P_rrt,
                 `P(seizure)` = s$P_sz, check.names = FALSE)
    }))
  }, ignoreNULL = FALSE)

  output$p_ph_outcome <- renderPlot({
    r <- ph_res()
    ggplot(r, aes(`urine pH`, `peak Cr ratio`)) +
      geom_line(linewidth = 1, colour = PAL[["Cr"]]) +
      geom_point(size = 2.5, colour = PAL[["Cr"]]) +
      labs(title = "What the trade costs at this burden",
           subtitle = "the minimum is the optimum urine pH for THIS patient",
           x = "urine pH", y = "peak creatinine ratio") + theme_tls()
  })
  output$t_ph <- renderTable(ph_res(), digits = 3)

  ## ---------------- tab 8 ----------------
  krescue <- eventReactive(input$go, { TLS_potassium_rescue() },
                           ignoreNULL = FALSE)
  output$t_krescue <- renderTable(krescue(), digits = 2)
  output$p_krescue <- renderPlot({
    r <- krescue()
    df <- r %>% select(rescue, K_2h, K_6h, K_12h, K_24h) %>%
      pivot_longer(-rescue) %>%
      mutate(hour = as.numeric(gsub("[^0-9]", "", name)))
    ggplot(df, aes(hour, value, colour = rescue)) +
      geom_line(linewidth = 1) + geom_point(size = 2) +
      geom_hline(yintercept = 6.0, linetype = 2, colour = "#c62828") +
      labs(title = "Serum potassium after each rescue",
           subtitle = paste("the shift arms come back up; see the dTBK column",
                            "for why"),
           x = "hours", y = "K⁺ (mmol/L)") + theme_tls()
  })

  ## ---------------- tab 9 ----------------
  output$t_cairo <- renderTable({
    cb <- cairo_bishop(sim())
    data.frame(criterion = c("urate ≥ 8 mg/dL or +25%",
                             "potassium ≥ 6.0 or +25%",
                             "phosphate ≥ 1.45 or +25%",
                             "calcium ≤ 1.75 or −25%",
                             "criteria met", "laboratory TLS", "clinical TLS"),
               value = c(ifelse(cb$urate, "yes", "no"),
                         ifelse(cb$K, "yes", "no"),
                         ifelse(cb$PO4, "yes", "no"),
                         ifelse(cb$Ca, "yes", "no"),
                         as.character(cb$n),
                         ifelse(cb$LTLS, "YES", "no"),
                         ifelse(cb$CTLS, "YES", "no")))
  })

  output$p_hazard <- renderPlot({
    d <- sim()
    df <- data.frame(t = d$time, arrhythmia = d$P_ARR, seizure = d$P_SZ,
                     `renal replacement` = d$P_RRT,
                     check.names = FALSE) %>% pivot_longer(-t)
    ggplot(df, aes(t / 24, 100 * value, colour = name)) +
      geom_line(linewidth = 1) +
      labs(title = "Cumulative event probability",
           subtitle = "shapes, not calibrated absolute risks",
           x = "days", y = "%") + theme_tls()
  })

  output$p_oxidant <- renderPlot({
    d <- sim()
    df <- data.frame(t = d$time, `cumulative H₂O₂ (mmol)` = d$H2O2,
                     `methaemoglobin (%)` = d$MetHb,
                     `haemoglobin (g/dL)` = d$HB,
                     check.names = FALSE) %>% pivot_longer(-t)
    ggplot(df, aes(t / 24, value)) + geom_line(linewidth = 1, colour = "#c62828") +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Urate oxidase oxidant limb",
           subtitle = ifelse(input$g6pd, "G6PD DEFICIENT",
                             "normal G6PD — tick the box to see the deficient case"),
           x = "days", y = NULL) + theme_tls()
  })

  ## ---------------- tab 10 ----------------
  output$t_params <- renderTable({
    p <- as.list(param(mod))
    data.frame(parameter = names(p), value = unlist(p, use.names = FALSE))
  }, digits = 5)
}

shinyApp(ui, server)
