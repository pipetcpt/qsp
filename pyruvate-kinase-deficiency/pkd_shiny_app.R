## =============================================================================
##  PYRUVATE KINASE DEFICIENCY (PKLR) -- INTERACTIVE QSP DASHBOARD
## =============================================================================
##  Front end for `pkd_mrgsolve_model.R` (80 ODEs, 14 red cell age cohorts).
##
##  THE APP IS BUILT AROUND ONE ARGUMENT, and the tab order follows it:
##
##    a pyruvate kinase lesion lowers ATP AND raises 2,3-BPG, because the enzyme
##    sits at the last ATP-generating step, DOWNSTREAM of the Rapoport-Luebering
##    branch point.  Those two effects hit the patient from opposite sides.  So a
##    PK activator raises haemoglobin BY LOWERING 2,3-BPG -- adding oxygen
##    carrier while subtracting oxygen unloading.
##
##  Tab 4 ("Oxygen transport") is therefore not a secondary readout, it is the
##  point of the app: it puts haemoglobin next to the three quantities that
##  actually describe oxygen transport (equivalent haemoglobin, required cardiac
##  output, tissue PO2) and lets the user watch them disagree.  Every plot that
##  shows haemoglobin also shows equivalent haemoglobin, deliberately, so the
##  user cannot read one without the other.
##
##  Run:  shiny::runApp("pkd_shiny_app.R")
##  Requires: shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
## =============================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

## The model source lives in pkd_mrgsolve_model.R; source it for `mod` and GENO.
## We source in a child environment so the scenario block at the bottom of that
## file does not run on app start-up.
.pkd_env <- new.env()
local({
  src <- readLines("pkd_mrgsolve_model.R")
  cut <- grep("^##  SCENARIOS", src)
  if (length(cut)) src <- src[seq_len(cut[1] - 3)]
  eval(parse(text = paste(src, collapse = "\n")), envir = .pkd_env)
}, envir = .pkd_env)
mod <- get("mod", envir = .pkd_env)

GENO <- tibble::tribble(
  ~label,                        ~ALPHA, ~TAUPK, ~ACTIVAT,
  "very mild / compensated",       0.30,   60,      1.0,
  "mild",                          0.22,   55,      1.0,
  "moderate",                      0.16,   50,      1.0,
  "severe, not transfused",        0.12,   50,      1.0,
  "transfusion dependent",         0.09,   45,      1.0,
  "non-missense / null allele",    0.12,   50,      0.0
)

CAPS <- c("HB", "RETPCT", "DPGTOT", "DPGFREE", "P50", "PVO2T", "COREQ", "HBEQ",
          "EXTB", "SAO2", "BILI", "ATPYNG", "ATPOLD", "PKREL", "SEQFRAC",
          "TSATT", "NCELL")

burn_in <- function(p, days = 800) {
  do.call(param, c(list(mod), p)) %>%
    mrgsim(end = days, delta = days, hmax = 0.5) %>%
    as.data.frame() %>% slice(n())
}

init_from <- function(row) {
  keep <- intersect(names(row), mrgsolve::cmt(mod))
  as.list(row[, keep, drop = FALSE])
}

run_one <- function(p, days, dose = NULL, y0 = NULL) {
  if (is.null(y0)) y0 <- burn_in(p)
  m <- do.call(param, c(list(mod), p)) %>% init(init_from(y0))
  if (!is.null(dose)) m <- m %>% ev(dose)
  m %>% mrgsim(end = days, delta = max(0.25, days / 600), hmax = 0.5) %>%
    as.data.frame()
}

mit_ev <- function(mg, days) {
  if (mg <= 0) return(NULL)
  ev(amt = mg, cmt = "MGUT", ii = 0.5, addl = max(0, days / 0.5 - 1))
}
eta_ev <- function(mg, days) {
  if (mg <= 0) return(NULL)
  ev(amt = mg, cmt = "EGUT", ii = 1, addl = max(0, days - 1))
}
teb_ev <- function(mg, days) {
  if (mg <= 0) return(NULL)
  ev(amt = mg, cmt = "TGUT", ii = 1, addl = max(0, days - 1))
}
dfx_ev <- function(mg, days) {
  if (mg <= 0) return(NULL)
  ev(amt = mg, cmt = "DGUT", ii = 1, addl = max(0, days - 1))
}
tx_ev <- function(units, ii, days) {
  if (units <= 0) return(NULL)
  ev(amt = units * 1.05 * 10 / 30, cmt = "ND1", ii = ii,
     addl = max(0, floor(days / ii) - 1))
}

thm <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("Pyruvate Kinase Deficiency (PKLR) — QSP dashboard"),
  tags$p(style = "color:#555;margin-top:-8px",
    tags$b("One lesion, two opposing consequences: "),
    "pyruvate kinase is the last ATP-generating step of glycolysis and sits ",
    "downstream of the Rapoport–Luebering branch point, so losing it lowers ATP ",
    "(the cell dies sooner) and raises 2,3-BPG (each surviving gram of ",
    "haemoglobin unloads more oxygen). A PK activator therefore raises ",
    "haemoglobin by lowering 2,3-BPG. Watch tab 4."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      selectInput("geno", "Genotype / phenotype", choices = GENO$label,
                  selected = GENO$label[4]),
      sliderInput("alpha", "Residual PKR activity of fresh protein (α)",
                  0.04, 0.60, 0.12, 0.01),
      sliderInput("taupk", "PKR thermolability τ_PK (d)", 15, 170, 50, 5),
      sliderInput("activ", "Activatable protein fraction (0 = null allele)",
                  0, 1, 1, 0.05),
      selectInput("ugt", "UGT1A1 genotype",
                  c("*1/*1 (normal)" = 1.0, "*1/*28 (Gilbert)" = 0.7,
                    "*28/*28" = 0.3), selected = 1.0),
      checkboxInput("splen", "Splenectomised", FALSE),
      checkboxInput("chol", "Cholecystectomy performed", FALSE),
      hr(),
      h4("Therapy"),
      sliderInput("mit", "Mitapivat (mg BID)", 0, 100, 50, 5),
      sliderInput("eta", "Etavopivat (mg OD)", 0, 600, 0, 50),
      sliderInput("teb", "Tebapivat (mg OD)", 0, 1, 0, 0.1),
      sliderInput("dfx", "Deferasirox (mg/d)", 0, 2000, 0, 250),
      sliderInput("txu", "Transfusion (units per episode)", 0, 4, 0, 1),
      sliderInput("txii", "Transfusion interval (d)", 14, 56, 21, 7),
      sliderInput("gtf", "Gene therapy: fraction of marrow output corrected",
                  0, 1, 0, 0.05),
      selectInput("cyp", "Co-medication (CYP3A)",
                  c("none" = 1.0, "strong inducer (rifampicin)" = 3.0,
                    "inhibitor (fluconazole)" = 0.45), selected = 1.0),
      hr(),
      sliderInput("days", "Simulate (days)", 60, 1460, 365, 30),
      actionButton("go", "Run", class = "btn-primary btn-block")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        ## ---- 1 -----------------------------------------------------------
        tabPanel("1 · Patient profile",
          h4("Untreated steady state for this genotype"),
          p("Set before any therapy is applied. The point of this tab is that ",
            "haemoglobin SATURATES as a severity readout: above the compensation ",
            "knee the marrow replaces what is destroyed and Hb is flat, while ",
            "reticulocytes and bilirubin keep rising monotonically."),
          DTOutput("tblProfile"),
          br(),
          plotOutput("plotSpectrum", height = "380px"),
          p(em("Critical red cell lifespan L* = 120 d ÷ (maximum amplification ÷ ",
               "basal amplification) = 120 ÷ 6.7 ≈ 18 d. Left of that line the ",
               "marrow is saturated and Hb tracks lifespan; right of it Hb is flat."))
        ),
        ## ---- 2 -----------------------------------------------------------
        tabPanel("2 · Drug PK",
          h4("Pharmacokinetics of the activators"),
          p("Mitapivat is a CYP3A substrate AND a moderate inducer, so it ",
            "accelerates its own clearance: steady-state exposure is lower than ",
            "a single dose predicts. Only the RED CELL compartment drives effect."),
          plotOutput("plotPK", height = "560px")
        ),
        ## ---- 3 -----------------------------------------------------------
        tabPanel("3 · Red cell metabolism",
          h4("Where the lesion actually bites, along the cell age axis"),
          p("Mutant PK-R is thermolabile and an anucleate cell cannot replace it, ",
            "so the lesion DEEPENS with age. An allosteric activator multiplies ",
            "protein that is present and therefore cannot rescue the oldest ",
            "cohorts; only the slower thermostabilisation limb (which lengthens ",
            "τ_PK itself) reaches them."),
          plotOutput("plotMetab", height = "520px"),
          br(),
          plotOutput("plotStab", height = "300px")
        ),
        ## ---- 4 -----------------------------------------------------------
        tabPanel("4 · Oxygen transport ★",
          h4("Haemoglobin is not a sufficient statistic"),
          p(strong("This is the tab the model was built for."), " Haemoglobin and ",
            "2,3-BPG move in opposite directions under treatment. The three ",
            "curves next to Hb are three closures on the same physiology: ",
            "EQUIVALENT haemoglobin (the Hb a normal-P50 subject would need to ",
            "unload the same oxygen), the cardiac output REQUIRED at a fixed ",
            "tissue PO2, and the tissue PO2 achieved at a fixed cardiac output. ",
            "The model cannot arbitrate between them, so it reports all three."),
          plotOutput("plotO2", height = "540px"),
          br(),
          h4("Break-even: how far may 2,3-BPG fall before the Hb gain is spent?"),
          DTOutput("tblBreakeven"),
          p(em("Read the last column as a falsifiable threshold. At a +1.5 g/dL ",
               "response — the ACTIVATE primary endpoint — the drug may lower ",
               "2,3-BPG by at most this much before the oxygen-transport gain ",
               "becomes a loss."))
        ),
        ## ---- 5 -----------------------------------------------------------
        tabPanel("5 · Haemolysis & clinical endpoints",
          h4("Haemoglobin, reticulocytes, bilirubin, LDH, haptoglobin"),
          plotOutput("plotHem", height = "560px"),
          br(),
          DTOutput("tblEnd")
        ),
        ## ---- 6 -----------------------------------------------------------
        tabPanel("6 · Spleen & sequestration",
          h4("Holding versus killing — and why they cannot be separated"),
          p("The red pulp is modelled as a metabolic stress test: glycolytic ",
            "capacity is reduced, mechanical load raised, and a reticulocyte ",
            "loses the mitochondrial subsidy exactly where its pump load is ",
            "highest. It is also RETAINED there. With reversible pooling the ",
            "fractional rise in reticulocytes after splenectomy must EQUAL the ",
            "fractional rise in haemoglobin — both are the same wasted fraction ",
            "of marrow output. That is a prediction the registry data can refute."),
          plotOutput("plotSpleen", height = "480px"),
          br(),
          DTOutput("tblSplx")
        ),
        ## ---- 7 -----------------------------------------------------------
        tabPanel("7 · Iron & hepcidin",
          h4("The benefit the haemoglobin endpoint cannot see"),
          p("PK deficiency loads iron WITHOUT transfusion: marrow expansion ",
            "raises erythroferrone, which suppresses hepcidin, which lifts the ",
            "brake on duodenal absorption. A PK activator reduces the ",
            "destruction the marrow was compensating for, so erythroferrone ",
            "falls and hepcidin rises. This benefit is LARGEST in patients whose ",
            "haemoglobin barely moves — the compensated ones, whom a ",
            "haemoglobin-response endpoint calls non-responders."),
          plotOutput("plotIron", height = "560px")
        ),
        ## ---- 8 -----------------------------------------------------------
        tabPanel("8 · Scenario comparison",
          h4("Standard regimens, same patient"),
          p("Each row is burned in to the untreated steady state for the ",
            "selected genotype first, then the intervention is applied."),
          actionButton("goCmp", "Run comparison (slow)", class = "btn-warning"),
          br(), br(),
          DTOutput("tblCmp"),
          br(),
          plotOutput("plotCmp", height = "420px")
        ),
        ## ---- 9 -----------------------------------------------------------
        tabPanel("9 · Diagnostics & biomarkers",
          h4("Why the enzyme assay can read normal"),
          p("Assayed red cell PK activity is a cohort-weighted mean and ",
            "reticulocytes carry freshly made enzyme, so the sicker the patient ",
            "the MORE the assay over-reads. This is the mechanism behind the ",
            "standard advice to interpret PK activity relative to another ",
            "age-dependent enzyme (the PK/hexokinase ratio) and to distrust a ",
            "normal absolute value in the presence of reticulocytosis."),
          plotOutput("plotDiag", height = "440px"),
          br(),
          DTOutput("tblDiag")
        ),
        ## ---- 10 ----------------------------------------------------------
        tabPanel("10 · Model & caveats",
          h4("What this model is, and what it cannot do"),
          htmlOutput("caveats")
        )
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  observeEvent(input$geno, {
    g <- GENO[GENO$label == input$geno, ]
    updateSliderInput(session, "alpha", value = g$ALPHA)
    updateSliderInput(session, "taupk", value = g$TAUPK)
    updateSliderInput(session, "activ", value = g$ACTIVAT)
  })

  pset <- reactive({
    list(ALPHA = input$alpha, TAUPK = input$taupk, ACTIVAT = input$activ,
         UGT = as.numeric(input$ugt), SPLEEN = as.numeric(!input$splen),
         CHOL = as.numeric(input$chol), GTFRAC = input$gtf,
         CYPEXT = as.numeric(input$cyp))
  })

  base0 <- reactive({
    p <- pset(); p$GTFRAC <- 0; p$SPLEEN <- 1
    burn_in(p)
  })

  sim <- eventReactive(input$go, {
    withProgress(message = "Solving 80 ODEs (exact glycolysis at every step)…", {
      p <- pset()
      d <- input$days
      doses <- Filter(Negate(is.null), list(
        mit_ev(input$mit, d), eta_ev(input$eta, d), teb_ev(input$teb, d),
        dfx_ev(input$dfx, d), tx_ev(input$txu, input$txii, d)))
      dose <- if (length(doses)) Reduce(c, doses) else NULL
      run_one(p, d, dose = dose, y0 = base0())
    })
  }, ignoreNULL = FALSE)

  ## ---- tab 1 -------------------------------------------------------------
  spectrum <- reactive({
    bind_rows(lapply(seq_len(nrow(GENO)), function(i) {
      g <- GENO[i, ]
      r <- burn_in(list(ALPHA = g$ALPHA, TAUPK = g$TAUPK, ACTIVAT = g$ACTIVAT))
      tibble(label = g$label, alpha = g$ALPHA,
             Hb = r$HB, ret = r$RETPCT, bili = r$BILI, DPG = r$DPGTOT,
             P50 = r$P50, HBEQ = r$HBEQ,
             lifespan = 120 * r$NCELL / max(r$NCELL, 1e-9))
    }))
  })

  output$tblProfile <- renderDT({
    r <- base0()
    tibble(quantity = c("haemoglobin (g/dL)", "EQUIVALENT haemoglobin (g/dL)",
                        "reticulocytes (%)", "total bilirubin (mg/dL)",
                        "2,3-BPG total (mM)", "free 2,3-BPG (mM)",
                        "P50 (mmHg)", "tissue PO2 (mmHg)",
                        "required cardiac output (L/min)",
                        "ATP, youngest cohort (mM)", "ATP, oldest cohort (mM)",
                        "assayed PK activity (rel.)",
                        "reticulocyte mass held in spleen (frac)",
                        "liver iron (mg Fe/g dw)", "serum ferritin (ng/mL)"),
           value = round(c(r$HB, r$HBEQ, r$RETPCT, r$BILI, r$DPGTOT, r$DPGFREE,
                           r$P50, r$PVO2T, r$COREQ, r$ATPYNG, r$ATPOLD,
                           r$PKREL, r$SEQFRAC, r$LIC, r$FERR), 3)) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  output$plotSpectrum <- renderPlot({
    s <- spectrum() %>%
      select(label, alpha, Hb, HBEQ, ret, bili, P50) %>%
      pivot_longer(-c(label, alpha))
    ggplot(s, aes(alpha, value)) +
      geom_line(colour = "#1f6fb2") + geom_point(size = 2) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "residual PKR activity of fresh protein (α)", y = NULL,
           title = "Haemoglobin saturates; reticulocytes and bilirubin do not") +
      thm
  })

  ## ---- tab 2 -------------------------------------------------------------
  output$plotPK <- renderPlot({
    d <- sim()
    d %>% transmute(time,
                    `mitapivat, plasma (mg/L)` = MC / 40,
                    `mitapivat, red cell (mg/L)` = MRBC,
                    `CYP3A4 (relative)` = CYP,
                    `etavopivat, red cell (mg/L)` = ERBC,
                    `tebapivat, red cell (mg/L)` = TRBC,
                    `deferasirox (mg/L)` = DCEN) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) +
      geom_line(colour = "#0f7a52") +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL,
           title = "Only the red cell compartment drives effect; CYP3A4 rises by auto-induction") +
      thm
  })

  ## ---- tab 3 -------------------------------------------------------------
  output$plotMetab <- renderPlot({
    d <- sim()
    d %>% transmute(time,
                    `ATP, youngest cohort (mM)` = ATPYNG,
                    `ATP, oldest cohort (mM)` = ATPOLD,
                    `free 2,3-BPG (mM)` = DPGFREE,
                    `total 2,3-BPG (mM)` = DPGTOT) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#b07d1a", linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL,
           title = "ATP up and 2,3-BPG down: the same drug action, opposite directions") +
      thm
  })

  output$plotStab <- renderPlot({
    d <- sim()
    d %>% transmute(time, `PKR stabilisation factor (× τ_PK)` = STAB) %>%
      ggplot(aes(time, `PKR stabilisation factor (× τ_PK)`)) +
      geom_line(colour = "#7b2d8e", linewidth = 0.9) +
      labs(x = "day", title = paste0("The slow limb: thermostabilisation takes ",
           "months and is the only thing that reaches OLD cohorts")) + thm
  })

  ## ---- tab 4 -------------------------------------------------------------
  output$plotO2 <- renderPlot({
    d <- sim()
    b <- base0()
    d %>% transmute(time,
      `haemoglobin (g/dL)` = HB,
      `EQUIVALENT haemoglobin (g/dL)` = HBEQ,
      `P50 (mmHg)` = P50,
      `O2 extracted per L blood` = EXTB,
      `required cardiac output (L/min)` = COREQ,
      `tissue PO2 (mmHg)` = PVO2T) %>%
      pivot_longer(-time) %>%
      mutate(grp = ifelse(grepl("EQUIVALENT|cardiac|tissue PO2|extracted", name),
                          "oxygen transport", "what the trial measures")) %>%
      ggplot(aes(time, value, colour = grp)) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = c("what the trial measures" = "#c0392b",
                                     "oxygen transport" = "#5b2d8e"), name = NULL) +
      labs(x = "day", y = NULL,
           title = "Haemoglobin rises. Oxygen transport does not follow it.") +
      thm
  })

  output$tblBreakeven <- renderDT({
    b <- base0()
    p <- as.list(param(mod))
    hill <- function(po2, p50) { r <- (po2 / p50)^p$NHILL; r / (1 + r) }
    ext <- function(p50) hill(p$PAO2, p50) - hill(p$PVO2REF, p50)
    e0 <- ext(b$P50)
    out <- lapply(c(1.0, 1.5, 2.0, 3.0), function(dhb) {
      hb1 <- b$HB + dhb
      target <- e0 * b$HB / hb1
      lo <- 10; hi <- 80
      for (k in 1:80) {
        mid <- (lo + hi) / 2
        if (ext(mid) < target) lo <- mid else hi <- mid
      }
      p50s <- (lo + hi) / 2
      dpgs <- p$DPG0 * (p50s / p$P50REF)^(1 / p$NDPG)
      tibble(`ΔHb (g/dL)` = dhb, `Hb reached` = round(hb1, 2),
             `break-even P50` = round(p50s, 2),
             `break-even 2,3-BPG (mM)` = round(dpgs, 3),
             `max 2,3-BPG fall` = sprintf("%.1f%%", 100 * (1 - dpgs / b$DPGTOT)))
    }) %>% bind_rows()
    datatable(out, rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- tab 5 -------------------------------------------------------------
  output$plotHem <- renderPlot({
    sim() %>% transmute(time,
      `haemoglobin (g/dL)` = HB, `reticulocytes (%)` = RETPCT,
      `total bilirubin (mg/dL)` = BILI, `LDH (U/L)` = LDH,
      `haptoglobin (g/L)` = HAP, `gallstone hazard index` = GS,
      `spleen volume (mL)` = SPLV, `transfused units` = TXU) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#8c1c1c", linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL) + thm
  })

  output$tblEnd <- renderDT({
    d <- sim(); b <- base0(); e <- slice(d, n())
    tibble(endpoint = c("haemoglobin (g/dL)", "EQUIVALENT haemoglobin (g/dL)",
                        "reticulocytes (%)", "total bilirubin (mg/dL)",
                        "LDH (U/L)", "liver iron (mg Fe/g dw)",
                        "transfused units", "P50 (mmHg)"),
           baseline = round(c(b$HB, b$HBEQ, b$RETPCT, b$BILI, b$LDH, b$LIC,
                              0, b$P50), 3),
           final = round(c(e$HB, e$HBEQ, e$RETPCT, e$BILI, e$LDH, e$LIC,
                           e$TXU, e$P50), 3),
           change = round(c(e$HB - b$HB, e$HBEQ - b$HBEQ, e$RETPCT - b$RETPCT,
                            e$BILI - b$BILI, e$LDH - b$LDH, e$LIC - b$LIC,
                            e$TXU, e$P50 - b$P50), 3)) %>%
      datatable(rownames = FALSE, options = list(dom = "t")) %>%
      formatStyle("change", color = styleInterval(0, c("#b03030", "#207020")))
  })

  ## ---- tab 6 -------------------------------------------------------------
  output$plotSpleen <- renderPlot({
    sim() %>% transmute(time,
      `sequestered reticulocytes (10^12/L)` = SEQR,
      `fraction of retic mass in spleen` = SEQFRAC,
      `spleen volume (mL)` = SPLV,
      `reticulocytes (%)` = RETPCT) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#4a7c2f", linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") + labs(x = "day", y = NULL) + thm
  })

  output$tblSplx <- renderDT({
    p <- pset()
    p$GTFRAC <- 0
    a <- p; a$SPLEEN <- 1
    b <- p; b$SPLEEN <- 0
    r0 <- burn_in(a)
    r1 <- run_one(b, 540, y0 = r0) %>% slice(n())
    tibble(quantity = c("haemoglobin (g/dL)", "reticulocytes (%)"),
           `spleen present` = round(c(r0$HB, r0$RETPCT), 3),
           `splenectomised` = round(c(r1$HB, r1$RETPCT), 3),
           `fractional change` = sprintf(
             "%+.1f%%", 100 * (c(r1$HB, r1$RETPCT) / c(r0$HB, r0$RETPCT) - 1))) %>%
      datatable(rownames = FALSE, options = list(dom = "t"),
                caption = paste("Model prediction: these two fractional changes",
                                "must be EQUAL if splenic reticulocyte",
                                "destruction is the mechanism."))
  })

  ## ---- tab 7 -------------------------------------------------------------
  output$plotIron <- renderPlot({
    sim() %>% transmute(time,
      `erythroferrone (ng/L)` = ERFE, `hepcidin (ng/L)` = HEP,
      `erythropoietin (IU/L)` = EPO, `transferrin saturation` = TSATT,
      `liver iron (mg Fe/g dw)` = LIC, `serum ferritin (ng/mL)` = FERR,
      `cardiac iron (mg Fe/g dw)` = CARDFE,
      `plasma iron (mg)` = FEP) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value)) + geom_line(colour = "#a06a00", linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "day", y = NULL,
           title = "ERFE down, hepcidin up, liver iron down — the van Beers 2024 signature") +
      thm
  })

  ## ---- tab 8 -------------------------------------------------------------
  cmp <- eventReactive(input$goCmp, {
    withProgress(message = "Running scenario panel…", {
      p <- pset(); p$GTFRAC <- 0; p$SPLEEN <- 1
      y0 <- burn_in(p)
      d <- 365
      scen <- list(
        `untreated`                = list(p = p, dose = NULL),
        `mitapivat 5 mg BID`       = list(p = p, dose = mit_ev(5, d)),
        `mitapivat 20 mg BID`      = list(p = p, dose = mit_ev(20, d)),
        `mitapivat 50 mg BID`      = list(p = p, dose = mit_ev(50, d)),
        `etavopivat 400 mg OD`     = list(p = p, dose = eta_ev(400, d)),
        `splenectomy`              = list(p = modifyList(p, list(SPLEEN = 0)), dose = NULL),
        `splenectomy + mitapivat`  = list(p = modifyList(p, list(SPLEEN = 0)), dose = mit_ev(50, d)),
        `gene therapy 45%`         = list(p = modifyList(p, list(GTFRAC = 0.45)), dose = NULL)
      )
      bind_rows(lapply(names(scen), function(nm) {
        s <- scen[[nm]]
        r <- if (is.null(s$dose) && identical(s$p, p)) y0 else
          slice(run_one(s$p, d, dose = s$dose, y0 = y0), n())
        tibble(scenario = nm, Hb = r$HB, HBEQ = r$HBEQ, ret = r$RETPCT,
               DPG = r$DPGTOT, P50 = r$P50, COREQ = r$COREQ,
               bili = r$BILI, LIC = r$LIC)
      }))
    })
  })

  output$tblCmp <- renderDT({
    cmp() %>% mutate(across(where(is.numeric), ~round(.x, 3))) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 12))
  })

  output$plotCmp <- renderPlot({
    cmp() %>% select(scenario, Hb, HBEQ) %>%
      pivot_longer(-scenario) %>%
      mutate(name = recode(name, Hb = "haemoglobin",
                           HBEQ = "EQUIVALENT haemoglobin (P50-corrected)")) %>%
      ggplot(aes(reorder(scenario, value), value, fill = name)) +
      geom_col(position = "dodge") + coord_flip() +
      scale_fill_manual(values = c("haemoglobin" = "#c0392b",
        "EQUIVALENT haemoglobin (P50-corrected)" = "#5b2d8e"), name = NULL) +
      labs(x = NULL, y = "g/dL",
           title = "The two bars disagree. That disagreement is the model's point.") +
      thm
  })

  ## ---- tab 9 -------------------------------------------------------------
  output$plotDiag <- renderPlot({
    spectrum() %>%
      mutate(assayed = NA_real_) -> s
    d <- bind_rows(lapply(seq_len(nrow(GENO)), function(i) {
      g <- GENO[i, ]
      r <- burn_in(list(ALPHA = g$ALPHA, TAUPK = g$TAUPK, ACTIVAT = g$ACTIVAT))
      tibble(label = g$label, `true α` = g$ALPHA, `assayed activity` = r$PKREL,
             `over-read ratio` = r$PKREL / g$ALPHA, `reticulocytes (%)` = r$RETPCT)
    }))
    d %>% pivot_longer(-c(label, `true α`)) %>%
      ggplot(aes(`true α`, value)) + geom_line(colour = "#444") +
      geom_point(size = 2) + facet_wrap(~name, scales = "free_y") +
      labs(y = NULL, title = paste("The sicker the patient, the more the assay",
                                   "over-reads")) + thm
  })

  output$tblDiag <- renderDT({
    bind_rows(lapply(seq_len(nrow(GENO)), function(i) {
      g <- GENO[i, ]
      r <- burn_in(list(ALPHA = g$ALPHA, TAUPK = g$TAUPK, ACTIVAT = g$ACTIVAT))
      tibble(phenotype = g$label, `true α` = g$ALPHA,
             `assayed (rel.)` = round(r$PKREL, 3),
             `ratio` = round(r$PKREL / g$ALPHA, 3),
             `reticulocytes (%)` = round(r$RETPCT, 1))
    })) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- tab 10 ------------------------------------------------------------
  output$caveats <- renderUI({
    HTML('
<p><b>Structure.</b> 80 ODEs. 14 red cell age cohorts each carrying their own
adenylate pool, a splenic sequestration compartment, 7 cohorts for
gene-therapy-corrected cells, 3 transfused-cell states, marrow (progenitor +
two erythroblast pools), erythropoietin, the erythroferrone–hepcidin–iron axis
with liver/cardiac/macrophage compartments, bilirubin and hepatobiliary
consequences, three PK-R activators, deferasirox, and endocrine/bone/vascular
sequelae. The fast glycolytic subsystem is solved to quasi-steady state exactly
at every derivative evaluation, for 29 distinct cell populations.</p>

<p><b>Four fitted parameters</b>: marrow output (set so a wild-type subject sits
at Hb 15.0 g/dL), the mutant PKR decay constant, and the two haemolysis gains.
Everything else is literature or back-calculated from the normal erythrocyte
operating point. Nineteen healthy-physiology quantities are then
<i>predictions</i> — red cell lifespan, reticulocyte percentage, bilirubin
production, the 25 mg/d iron recycling flux — and are checked in
<code>pkd_reference_output.txt</code> §1 before any patient is simulated.</p>

<p><b>Thirteen defects</b> were found by running the equations rather than by
reading them; each is logged with its fix site in
<code>pkd_reference_output.txt</code> §7. The instructive ones were not coding
slips but places where a physiologically plausible equation gave a
quantitatively impossible answer — a bare phosphoglycerate-kinase-equilibrium
form that predicted a 9% rather than 2-fold 2,3-BPG rise; Michaelis-Menten
bilirubin conjugation that assigned a moderately affected patient a bilirubin of
23 000 mg/dL; aggregate rather than per-cohort mass balance that <i>created</i>
red cells and reported haemoglobins of 54 g/dL.</p>

<p><b>Limitations, stated so they are not mistaken for results.</b></p>
<ul>
<li>pH is not a state. The Bohr effect and the strong pH dependence of both BPGM
activities are folded into fixed constants, so acid-base disturbance cannot be
simulated.</li>
<li>Free versus Mg-bound ADP is not resolved; the PGK equilibrium uses total ADP,
which overstates how far ATP/ADP falls.</li>
<li><b>The 3-PG inhibition constant of the 2,3-BPG phosphatase (30 µM) is the
single most load-bearing parameter and the least well pinned by data.</b> The
oxygen-transport conclusion scales with the size of the 2,3-BPG excursion, so it
inherits that uncertainty. The break-even table on tab 4 exists precisely so the
conclusion can be re-scored against a different 2,3-BPG estimate.</li>
<li>Splenic destruction and sequestration share one compartment with a mean-field
transit; there is no distribution of transit times.</li>
<li>The oxygen module is whole-body and single-tissue, so it cannot represent the
regional differences in extraction where a right-shifted curve actually helps or
hurts most.</li>
<li>Cardiac output responds to haemoglobin, not to tissue PO2, so the model
cannot arbitrate between the three closures on tab 4 — it reports all of them.</li>
</ul>

<p><b>This is an educational and research model.</b> It has not been validated
against individual patient data and must not be used for clinical decisions,
prescribing, or regulatory submission.</p>')
  })
}

shinyApp(ui, server)
