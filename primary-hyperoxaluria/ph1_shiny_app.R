# =====================================================================
#  PRIMARY HYPEROXALURIA (PH1 / PH2 / PH3) — Shiny dashboard
#  ---------------------------------------------------------------------
#  A front end for ph1_mrgsolve_model.R. The organising idea of the model
#  is carried into the UI: the app never lets you set a disease severity.
#  You set a GENOTYPE (residual enzyme activity) and a THERAPY, and every
#  number on every tab is an output of the mass balance
#
#     production = renal excretion + enteric elimination + d(tissue)/dt
#
#  Eleven tabs:
#    1  Mass balance          the balance sheet itself, drawn as a budget
#    2  Genotype              the AGT nonlinearity, live
#    3  Urinary oxalate       the efficacy biomarker while eGFR is intact
#    4  Plasma oxalate        the biomarker past the threshold
#    5  Biomarker inversion   the eGFR at which Uox stops meaning anything
#    6  Kidney                eGFR trajectory, the crystal-injury cascade
#    7  Stones                AP(CaOx) index and event rate (Layer B)
#    8  Systemic oxalosis     bone, retina, heart, nerve, skin, marrow
#    9  Drug PK/PD            siRNA plasma vs RISC vs enzyme vs effect
#   10  Dialysis & transplant the arithmetic of removal, and unloading
#   11  Scenario comparison   any two arms side by side
#
#  RUN:  shiny::runApp("ph1_shiny_app.R")
#  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
# =====================================================================

suppressMessages({
  library(shiny); library(mrgsolve); library(dplyr)
  library(tidyr); library(ggplot2); library(DT)
})

# ---- load the model ------------------------------------------------
# The model file also runs 34 scenarios and 16 diagnostics when sourced
# directly; here we only want the compiled model object, so the driver
# sections are skipped by reading the file up to the driver marker.
MODEL_FILE <- "ph1_mrgsolve_model.R"
src <- readLines(MODEL_FILE, warn = FALSE)
cut <- grep("^#\\s+2\\. SIMULATION INFRASTRUCTURE", src)
if (length(cut)) src <- src[1:(cut[1] - 3)]
eval(parse(text = paste(src, collapse = "\n")), envir = globalenv())

YR <- 365.25
STATES <- names(mrgsolve::init(mod))

GENOTYPES <- list(
  "Healthy (AGT 100%)"                       = list(FAGT = 1.00, B6RESC = 0.000),
  "Heterozygous carrier (AGT 50%)"           = list(FAGT = 0.50, B6RESC = 0.000),
  "PH1 attenuated (AGT 10%)"                 = list(FAGT = 0.10, B6RESC = 0.000),
  "PH1 G170R homozygous (B6-responsive)"     = list(FAGT = 0.005, B6RESC = 0.035),
  "PH1 F152I (B6-responsive)"                = list(FAGT = 0.010, B6RESC = 0.028),
  "PH1 I244T (partially B6-responsive)"      = list(FAGT = 0.005, B6RESC = 0.015),
  "PH1 G82E (PLP site destroyed)"            = list(FAGT = 0.005, B6RESC = 0.000),
  "PH1 null / frameshift (B6-unresponsive)"  = list(FAGT = 0.000, B6RESC = 0.000),
  "PH2 (GRHPR null)"                         = list(FGRHPR = 0),
  "PH3 (HOGA1 null, GRHPR-inhibition only)"  = list(FHOGA1 = 0),
  "PH3 (HOGA1 null + cytosolic aldolase)"    = list(FHOGA1 = 0, HOGCYT = 1.8),
  "Enteric hyperoxaluria (AGT normal)"       = list(OX_DIET = 4200, FABS_OX = 0.30)
)

# ---- one simulation ------------------------------------------------
simulate <- function(g, tx_age, yrs, lum, ned, stp, pn, uvol, kcit, mg, thiaz,
                     oxf, hd, hdh, pd, il1, asc, growth, steroid) {
  pars <- c(GENOTYPES[[g]],
            list(U_VOL_TGT = uvol, KCIT_DOSE = kcit, MG_DOSE = mg,
                 THIAZ = thiaz, OXF_DEG = oxf, PN_MGKG = pn,
                 HD_PERWK = hd, HD_HRS = hdh, PD_ON = as.numeric(pd),
                 IL1_BLOCK = il1, ASC_DOSE = asc, GROWTH = growth,
                 STEROID = steroid))
  t0 <- tx_age * YR
  n  <- max(yrs - tx_age, 0)
  bag <- list()
  if (lum) bag <- c(bag, list(ev(amt = 12900, cmt = "LUM_SC", ii = 30, addl = 2, time = t0),
                              ev(amt = 12900, cmt = "LUM_SC", ii = 90,
                                 addl = ceiling(n * 4) + 1, time = t0 + 90)))
  if (ned) bag <- c(bag, list(ev(amt = 10000, cmt = "NED_SC", ii = 30,
                                 addl = ceiling(n * 12) + 1, time = t0)))
  if (stp) bag <- c(bag, list(ev(amt = 1750, cmt = "STP_GUT", ii = 0.5,
                                 addl = ceiling(n * 730) + 1, time = t0)))
  evs <- if (!length(bag)) NULL else Reduce(function(a, b) c(a, b), bag)
  mod %>% param(pars) %>%
    mrgsim(events = evs, end = yrs * YR, delta = 15, hmax = 0.5) %>% as_tibble()
}

thm <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom")

vline_tx <- function(tx) geom_vline(xintercept = tx, linetype = 3, colour = "grey35")

# =====================================================================
#  UI
# =====================================================================
ui <- fluidPage(
  titlePanel("Primary Hyperoxaluria — QSP dashboard (PH1 / PH2 / PH3)"),
  tags$p(tags$em(paste(
    "Oxalate is a terminal metabolite: production = renal excretion +",
    "enteric elimination + d(tissue burden)/dt.",
    "There is no severity slider in this app — you set a genotype and a",
    "therapy, and everything else is an output of that balance."))),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("g", "Genotype / disease", names(GENOTYPES),
                  selected = "PH1 G170R homozygous (B6-responsive)"),
      sliderInput("yrs", "Follow-up (years of age)", 5, 60, 40, 1),
      sliderInput("tx",  "Age therapy starts (y)",   0, 40,  5, 1),
      tags$hr(), tags$strong("LAYER A — production"),
      checkboxInput("lum", "Lumasiran 3 mg/kg (monthly x3, then quarterly)", FALSE),
      checkboxInput("ned", "Nedosiran 160 mg SC monthly", FALSE),
      checkboxInput("stp", "Stiripentol 50 mg/kg/day (oral, reaches kidney)", FALSE),
      sliderInput("pn", "Pyridoxine (mg/kg/day)", 0, 30, 0, 1),
      tags$hr(), tags$strong("LAYER B — crystallisation"),
      sliderInput("uvol", "Target urine volume (L/day)", 0.8, 4.0, 1.5, 0.1),
      sliderInput("kcit", "Potassium citrate (mmol citrate/day)", 0, 12, 0, 1),
      sliderInput("mg",   "Magnesium (mmol/day)", 0, 12, 0, 1),
      sliderInput("thiaz","Thiazide effect on urinary Ca (fraction)", 0, 0.5, 0, 0.05),
      tags$hr(), tags$strong("LAYER C — removal"),
      sliderInput("oxf", "Enteric oxalate degradation (fraction)", 0, 0.9, 0, 0.05),
      sliderInput("hd",  "Haemodialysis sessions/week", 0, 7, 0, 1),
      sliderInput("hdh", "Hours per HD session", 2, 8, 4, 0.5),
      checkboxInput("pd", "Peritoneal dialysis", FALSE),
      tags$hr(), tags$strong("Experimental / modifiers"),
      sliderInput("il1", "Anti-IL-1 / NLRP3 blockade (fraction)", 0, 0.9, 0, 0.05),
      sliderInput("asc", "Extra oxalate from ascorbate (umol/day)", 0, 600, 0, 25),
      sliderInput("growth", "Collagen-turnover / growth multiplier", 0.5, 2.0, 1.0, 0.1),
      sliderInput("steroid", "Corticosteroid effect on bone resorption", 0, 2, 0, 0.25),
      tags$hr(),
      tags$small(tags$em(paste(
        "Educational / research model only. Not validated for clinical",
        "use, prescribing, or regulatory submission.")))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 Mass balance", br(),
                 fluidRow(column(6, plotOutput("p_bal", height = 330)),
                          column(6, plotOutput("p_burden", height = 330))),
                 h4("The balance sheet at the end of follow-up"),
                 DTOutput("t_bal"),
                 helpText(paste(
                   "Production is set by the genotype and reduced only by Layer A.",
                   "Renal excretion is the product of oxalate clearance and plasma",
                   "oxalate, and that clearance is itself destroyed by the oxalate.",
                   "Whatever the three sinks do not remove is deposited in tissue —",
                   "that residual is the disease."))),
        tabPanel("2 Genotype", br(),
                 plotOutput("p_geno", height = 380),
                 h4("Residual AGT activity vs urinary oxalate"),
                 DTOutput("t_geno"),
                 helpText(paste(
                   "AGT is a high-Vmax, low-Km enzyme, so this curve is steeply",
                   "nonlinear. 50% of the enzyme is indistinguishable from 100%",
                   "(which is why PH1 is recessive and carriers are healthy), while",
                   "moving from 0.5% to 3% — a chaperone-sized gain — removes about",
                   "a third of the oxalate (which is why pyridoxine is worth giving",
                   "to the genotypes it can rescue)."))),
        tabPanel("3 Urinary oxalate", br(),
                 plotOutput("p_uox", height = 320),
                 plotOutput("p_glyc", height = 280),
                 helpText(paste(
                   "Upper: 24-h urinary oxalate against the 0.46 mmol/1.73m2/day",
                   "upper limit of normal. Lower: urinary glycolate. On lumasiran the",
                   "glycolate rise is not a side effect — it is the receipt for HAO1",
                   "knockdown, because the blocked substrate has to go somewhere."))),
        tabPanel("4 Plasma oxalate", br(),
                 plotOutput("p_pox", height = 340),
                 plotOutput("p_gx", height = 260),
                 helpText(paste(
                   "The dashed line is the calcium-oxalate solubility limit of plasma",
                   "(~30 umol/L). Above it, oxalate has nowhere to go but tissue and",
                   "the disease stops being a kidney disease. Lower panel: plasma",
                   "glyoxylate — it rises on an LDHA-directed drug because that drug",
                   "blocks the exit rather than the entrance."))),
        tabPanel("5 Biomarker inversion", br(),
                 plotOutput("p_inv", height = 380),
                 h4("Production held constant; only nephron mass varied"),
                 DTOutput("t_inv"),
                 helpText(paste(
                   "Uox = oxalate clearance x plasma oxalate. Below the deposition",
                   "threshold Uox equals production and is a clean efficacy",
                   "biomarker. Past it, falling GFR pulls Uox DOWN while the disease",
                   "accelerates. A falling Uox in advanced PH is not a therapeutic",
                   "success — which is why ILLUMINATE-A (eGFR >= 30) could use 24-h",
                   "Uox and ILLUMINATE-C had to switch to plasma oxalate."))),
        tabPanel("6 Kidney", br(),
                 plotOutput("p_egfr", height = 320),
                 plotOutput("p_casc", height = 300),
                 helpText(paste(
                   "There is no prescribed eGFR trajectory anywhere in the model.",
                   "Crystals adhere, activate NLRP3, generate IL-1beta and TGF-beta,",
                   "and the resulting fibrosis removes nephrons; fewer nephrons means",
                   "less clearance, higher plasma oxalate and more deposition. The",
                   "gain of that loop rises as eGFR falls, which is what makes the",
                   "terminal decline abrupt."))),
        tabPanel("7 Stones", br(),
                 plotOutput("p_ap", height = 320),
                 plotOutput("p_stone", height = 280),
                 helpText(paste(
                   "The AP(CaOx) index is a Tiselius-style supersaturation surrogate",
                   "computed from oxalate, calcium, citrate, magnesium and volume.",
                   "Raising urine volume and adding citrate change NO mass-balance",
                   "term and still cut supersaturation by half — because the stone",
                   "endpoint is a ratio, not a flux. That is also the limit of",
                   "supportive care: it protects the kidney and never closes the",
                   "balance."))),
        tabPanel("8 Systemic oxalosis", br(),
                 plotOutput("p_ox", height = 340),
                 plotOutput("p_bone", height = 280),
                 helpText(paste(
                   "Organ deposition begins only once plasma oxalate crosses the",
                   "solubility limit, so these curves are flat for years and then",
                   "are not. Bone is modelled as two pools: a surface pool that",
                   "exchanges in months (this is what unloads after a transplant)",
                   "and a deep crystal-incorporated pool that empties over years."))),
        tabPanel("9 Drug PK/PD", br(),
                 fluidRow(column(6, plotOutput("p_pk", height = 300)),
                          column(6, plotOutput("p_risc", height = 300))),
                 plotOutput("p_enz", height = 280),
                 helpText(paste(
                   "A GalNAc-siRNA has a plasma half-life of hours and an effect",
                   "that lasts months. The active species is the RISC-loaded strand",
                   "inside the hepatocyte, which is why quarterly dosing works and",
                   "why the schedule is loading-then-maintenance. Note that renal",
                   "LDHA is never silenced: the asialoglycoprotein receptor is a",
                   "hepatocyte receptor, so delivery — not target — decides which",
                   "disease the drug treats."))),
        tabPanel("10 Dialysis & transplant", br(),
                 h4("Weekly oxalate budget on the current dialysis prescription"),
                 DTOutput("t_dial"),
                 plotOutput("p_dial", height = 320),
                 helpText(paste(
                   "Removal is integrated from the actual intradialytic clearance",
                   "acting on the actual (falling, then rebounding) plasma level, so",
                   "the weekly total is computed rather than assumed. In PH1 it does",
                   "not match weekly production: conventional thrice-weekly",
                   "haemodialysis leaves a shortfall of several mmol every week.",
                   "Cutting production with a Layer A drug is what turns an",
                   "insufficient prescription into a sufficient one."))),
        tabPanel("11 Scenario comparison", br(),
                 fluidRow(
                   column(4, selectInput("cg", "Comparator genotype", names(GENOTYPES),
                                         selected = "PH1 G170R homozygous (B6-responsive)")),
                   column(4, checkboxGroupInput("cdrug", "Comparator therapy",
                                                c("Lumasiran" = "lum", "Nedosiran" = "ned",
                                                  "Stiripentol" = "stp",
                                                  "Hyperhydration + citrate" = "layb",
                                                  "Pyridoxine 10 mg/kg" = "b6"))),
                   column(4, helpText("The left arm is whatever the sidebar says."))),
                 plotOutput("p_cmp", height = 460),
                 DTOutput("t_cmp"))
      )
    )
  )
)

# =====================================================================
#  SERVER
# =====================================================================
server <- function(input, output, session) {

  d <- reactive({
    simulate(input$g, input$tx, input$yrs, input$lum, input$ned, input$stp,
             input$pn, input$uvol, input$kcit, input$mg, input$thiaz,
             input$oxf, input$hd, input$hdh, input$pd, input$il1,
             input$asc, input$growth, input$steroid)
  })

  # ---- 1 mass balance ----------------------------------------------
  output$p_bal <- renderPlot({
    x <- d()
    df <- x %>% transmute(AGE_Y,
                          Production = PROD_TOT,
                          `Renal excretion` = c(0, diff(CUM_UOX) / diff(time)),
                          `Enteric` = c(0, diff(CUM_ENT) / diff(time)),
                          `Dialysis` = c(0, diff(CUM_HD) / diff(time))) %>%
      pivot_longer(-AGE_Y)
    ggplot(df, aes(AGE_Y, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(title = "Mass-balance terms", x = "age (years)",
           y = "umol/day", colour = NULL) + vline_tx(input$tx) + thm
  })
  output$p_burden <- renderPlot({
    x <- d()
    ggplot(x, aes(AGE_Y, TBURDEN)) +
      geom_area(fill = "#ef9a9a", alpha = 0.6) + geom_line(colour = "#b71c1c") +
      labs(title = "Accumulated tissue oxalate (the unclosed residual)",
           x = "age (years)", y = "mmol") + vline_tx(input$tx) + thm
  })
  output$t_bal <- renderDT({
    l <- tail(d(), 1)
    acc <- l$OX_PL + l$OX_ECF + l$OX_BONE + l$OX_BONE_D + l$OX_SOFT + l$OX_KID
    tibble(Term = c("cumulative production", "urinary excretion",
                    "enteric elimination", "dialytic removal",
                    "accumulated in tissue", "closure error (%)"),
           mmol = round(c(l$CUM_PROD, l$CUM_UOX, l$CUM_ENT, l$CUM_HD, acc,
                          100 * (l$CUM_PROD - l$CUM_UOX - l$CUM_ENT - l$CUM_HD - acc) /
                            max(l$CUM_PROD, 1)) / c(1000, 1000, 1000, 1000, 1000, 1), 4)) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # ---- 2 genotype ---------------------------------------------------
  geno_scan <- reactive({
    f <- c(1, 0.5, 0.25, 0.15, 0.10, 0.07, 0.05, 0.03, 0.02, 0.01, 0.005, 0)
    u <- sapply(f, function(z)
      tail(mod %>% param(FAGT = z, PN_MGKG = input$pn) %>%
             mrgsim(end = 3 * YR, delta = 90, hmax = 0.5) %>% as_tibble(), 1)$UOX_MMOL)
    tibble(`AGT (%)` = 100 * f, Uox = round(u, 3),
           `x normal` = round(u / u[1], 2), `x ULN` = round(u / 0.46, 2))
  })
  output$p_geno <- renderPlot({
    g <- geno_scan()
    ggplot(g, aes(`AGT (%)`, Uox)) + geom_line(linewidth = 1) + geom_point(size = 2) +
      geom_hline(yintercept = 0.46, linetype = 2, colour = "#c62828") +
      annotate("text", x = 60, y = 0.50, label = "upper limit of normal",
               colour = "#c62828", size = 3.4) +
      scale_x_log10() +
      labs(title = "The AGT nonlinearity: recessive inheritance and chaperone rescue in one curve",
           x = "residual AGT activity (% of normal, log scale)",
           y = "urinary oxalate (mmol/1.73m2/day)") + thm
  })
  output$t_geno <- renderDT(datatable(geno_scan(), rownames = FALSE,
                                      options = list(dom = "t", pageLength = 12)))

  # ---- 3 urinary oxalate -------------------------------------------
  output$p_uox <- renderPlot({
    ggplot(d(), aes(AGE_Y, UOX_MMOL)) + geom_line(linewidth = 1, colour = "#00695c") +
      geom_hline(yintercept = 0.46, linetype = 2, colour = "#c62828") +
      labs(title = "24-h urinary oxalate", x = "age (years)",
           y = "mmol/1.73m2/day") + vline_tx(input$tx) + thm
  })
  output$p_glyc <- renderPlot({
    ggplot(d(), aes(AGE_Y, UGLYC_MMOL)) + geom_line(linewidth = 1, colour = "#ef6c00") +
      labs(title = "Urinary glycolate — the on-target receipt for HAO1 knockdown",
           x = "age (years)", y = "mmol/day") + vline_tx(input$tx) + thm
  })

  # ---- 4 plasma oxalate --------------------------------------------
  output$p_pox <- renderPlot({
    ggplot(d(), aes(AGE_Y, POX_OUT)) + geom_line(linewidth = 1, colour = "#4a148c") +
      geom_hline(yintercept = 30, linetype = 2, colour = "#c62828") +
      annotate("text", x = min(d()$AGE_Y) + 2, y = 33,
               label = "plasma CaOx solubility limit", colour = "#c62828", size = 3.4,
               hjust = 0) +
      labs(title = "Plasma oxalate and the threshold that defines systemic oxalosis",
           x = "age (years)", y = "umol/L") + vline_tx(input$tx) + thm
  })
  output$p_gx <- renderPlot({
    ggplot(d(), aes(AGE_Y, PGX_OUT)) + geom_line(linewidth = 1, colour = "#1565c0") +
      labs(title = "Plasma glyoxylate — rises when the EXIT rather than the entrance is blocked",
           x = "age (years)", y = "umol/L") + vline_tx(input$tx) + thm
  })

  # ---- 5 inversion --------------------------------------------------
  inv <- reactive({
    st <- as.numeric(tail(simulate(input$g, 99, 2, FALSE, FALSE, FALSE, input$pn,
                                   input$uvol, input$kcit, input$mg, input$thiaz,
                                   input$oxf, 0, 4, FALSE, 0, input$asc,
                                   input$growth, 0), 1)[STATES])
    names(st) <- STATES
    st[c("OX_BONE", "OX_BONE_D", "OX_SOFT", "OX_KID", "FIBROSIS", "CRYST",
         "NLRP3", "IL1B", "TGFB")] <- 0
    nf <- c(1, 0.7, 0.5, 0.35, 0.25, 0.18, 0.12, 0.08, 0.05, 0.03, 0.02, 0.014)
    do.call(rbind, lapply(nf, function(n) {
      st2 <- st; st2["NEPHRON"] <- n
      r <- tail(mod %>% param(c(GENOTYPES[[input$g]], list(K_LOSS = 0, INIT_ON = 0))) %>%
                  init(as.list(st2)) %>%
                  mrgsim(end = 400, delta = 20, hmax = 0.5) %>% as_tibble(), 1)
      tibble(eGFR = round(r$EGFR, 1), Pox = round(r$POX_OUT, 2),
             Uox = round(r$UOX_MMOL, 3))
    })) %>% mutate(`Uox % of eGFR125` = round(100 * Uox / Uox[1], 1))
  })
  output$p_inv <- renderPlot({
    i <- inv()
    ggplot(i, aes(eGFR)) +
      geom_line(aes(y = 100 * Uox / i$Uox[1], colour = "Uox (% of eGFR 125)"),
                linewidth = 1) +
      geom_line(aes(y = Pox, colour = "plasma oxalate (umol/L)"), linewidth = 1) +
      geom_hline(yintercept = 30, linetype = 2, colour = "grey40") +
      scale_x_reverse() +
      labs(title = "The inversion: production is CONSTANT across this whole plot",
           x = "eGFR (mL/min/1.73m2, falling to the right)", y = NULL,
           colour = NULL) + thm
  })
  output$t_inv <- renderDT(datatable(inv(), rownames = FALSE,
                                     options = list(dom = "t", pageLength = 12)))

  # ---- 6 kidney -----------------------------------------------------
  output$p_egfr <- renderPlot({
    ggplot(d(), aes(AGE_Y, EGFR)) + geom_line(linewidth = 1, colour = "#0d47a1") +
      geom_hline(yintercept = c(60, 30, 15), linetype = 3, colour = "grey45") +
      annotate("text", x = min(d()$AGE_Y), y = c(63, 33, 18),
               label = c("CKD3", "CKD4", "ESKD"), hjust = 0, size = 3.2,
               colour = "grey35") +
      labs(title = "eGFR — an output of the crystal-injury loop, never prescribed",
           x = "age (years)", y = "mL/min/1.73m2") + vline_tx(input$tx) + thm
  })
  output$p_casc <- renderPlot({
    x <- d() %>% select(AGE_Y, CRYST, NLRP3, IL1B, TGFB, FIBROSIS, TUB_INJ) %>%
      pivot_longer(-AGE_Y)
    ggplot(x, aes(AGE_Y, value, colour = name)) + geom_line(linewidth = 0.85) +
      labs(title = "Crystal -> NLRP3 -> IL-1beta -> TGF-beta -> fibrosis",
           x = "age (years)", y = "relative", colour = NULL) +
      vline_tx(input$tx) + thm
  })

  # ---- 7 stones -----------------------------------------------------
  output$p_ap <- renderPlot({
    ggplot(d(), aes(AGE_Y, AP_CAOX)) + geom_line(linewidth = 1, colour = "#f57f17") +
      geom_hline(yintercept = 1, linetype = 2, colour = "#c62828") +
      labs(title = "AP(CaOx) index — a RATIO, which is why Layer B works without changing a flux",
           x = "age (years)", y = "AP(CaOx) index") + vline_tx(input$tx) + thm
  })
  output$p_stone <- renderPlot({
    ggplot(d(), aes(AGE_Y, STONE_RATE)) + geom_line(linewidth = 1, colour = "#6a1b9a") +
      labs(title = "Symptomatic stone event rate", x = "age (years)",
           y = "events/year") + vline_tx(input$tx) + thm
  })

  # ---- 8 systemic oxalosis -----------------------------------------
  output$p_ox <- renderPlot({
    x <- d() %>% select(AGE_Y, RETINA, MYOCARD, NERVE, SKIN, BONE_DIS, MARROW) %>%
      pivot_longer(-AGE_Y)
    ggplot(x, aes(AGE_Y, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(title = "Organ deposition scores — flat for years, then not",
           x = "age (years)", y = "score", colour = NULL) +
      vline_tx(input$tx) + thm
  })
  output$p_bone <- renderPlot({
    x <- d() %>% transmute(AGE_Y,
                           `bone surface (exchangeable)` = OX_BONE / 1000,
                           `bone deep (crystal-bound)` = OX_BONE_D / 1000,
                           `soft tissue` = OX_SOFT / 1000,
                           `kidney parenchyma` = OX_KID / 1000) %>%
      pivot_longer(-AGE_Y)
    ggplot(x, aes(AGE_Y, value, fill = name)) + geom_area(alpha = 0.75) +
      labs(title = "Where the unexcreted oxalate went", x = "age (years)",
           y = "mmol", fill = NULL) + vline_tx(input$tx) + thm
  })

  # ---- 9 drug PK/PD -------------------------------------------------
  output$p_pk <- renderPlot({
    x <- d() %>% filter(AGE_Y >= input$tx) %>%
      select(AGE_Y, LUM_PL, NED_PL) %>% pivot_longer(-AGE_Y)
    ggplot(x, aes(AGE_Y, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(title = "Plasma siRNA (hours)", x = "age (years)", y = "nmol",
           colour = NULL) + thm
  })
  output$p_risc <- renderPlot({
    x <- d() %>% filter(AGE_Y >= input$tx) %>%
      select(AGE_Y, LUM_RISC, NED_RISC) %>% pivot_longer(-AGE_Y)
    ggplot(x, aes(AGE_Y, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(title = "RISC-loaded siRNA (months) — the active species",
           x = "age (years)", y = "nmol", colour = NULL) + thm
  })
  output$p_enz <- renderPlot({
    x <- d() %>% transmute(AGE_Y, `GO protein (HAO1)` = FGO_OUT,
                           `hepatic LDHA` = FLDH_OUT, `renal LDHA` = LDHA_K) %>%
      pivot_longer(-AGE_Y)
    ggplot(x, aes(AGE_Y, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(title = "Enzyme levels — renal LDHA is never reached by a GalNAc-siRNA",
           x = "age (years)", y = "fraction of normal", colour = NULL) +
      vline_tx(input$tx) + thm
  })

  # ---- 10 dialysis --------------------------------------------------
  output$t_dial <- renderDT({
    x <- d(); n <- nrow(x); i <- which(x$time >= max(x$time) - 7)[1]
    g <- function(col) (x[[col]][n] - x[[col]][i]) / 1000
    pr <- g("CUM_PROD")
    tibble(Term = c("production", "dialytic removal", "urinary excretion",
                    "enteric elimination", "SHORTFALL deposited in tissue",
                    "plasma oxalate (umol/L)"),
           `mmol/week` = round(c(pr, g("CUM_HD"), g("CUM_UOX"), g("CUM_ENT"),
                                 pr - g("CUM_HD") - g("CUM_UOX") - g("CUM_ENT"),
                                 tail(x$POX_OUT, 1)), 2)) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })
  output$p_dial <- renderPlot({
    x <- d() %>% transmute(AGE_Y, production = PROD_TOT,
                           removal = c(0, diff(CUM_HD + CUM_UOX + CUM_ENT) / diff(time))) %>%
      pivot_longer(-AGE_Y)
    ggplot(x, aes(AGE_Y, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(title = "Production vs total removal (the gap is the disease)",
           x = "age (years)", y = "umol/day", colour = NULL) + thm
  })

  # ---- 11 comparison ------------------------------------------------
  cmp <- reactive({
    cd <- input$cdrug
    simulate(input$cg, input$tx, input$yrs,
             "lum" %in% cd, "ned" %in% cd, "stp" %in% cd,
             if ("b6" %in% cd) 10 else 0,
             if ("layb" %in% cd) 3.0 else 1.5,
             if ("layb" %in% cd) 4 else 0,
             0, 0, 0, input$hd, input$hdh, input$pd, 0, 0, 1, 0)
  })
  output$p_cmp <- renderPlot({
    a <- d() %>% mutate(arm = "sidebar arm")
    b <- cmp() %>% mutate(arm = "comparator")
    x <- bind_rows(a, b) %>%
      select(AGE_Y, arm, UOX_MMOL, POX_OUT, EGFR, AP_CAOX, TBURDEN, OXALOSIS) %>%
      pivot_longer(-c(AGE_Y, arm))
    ggplot(x, aes(AGE_Y, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "age (years)", y = NULL, colour = NULL) + thm
  })
  output$t_cmp <- renderDT({
    f <- function(x, nm) tail(x, 1) %>%
      transmute(arm = nm, Uox = round(UOX_MMOL, 3), Pox = round(POX_OUT, 2),
                eGFR = round(EGFR, 1), AP = round(AP_CAOX, 2),
                `stones/y` = round(STONE_RATE, 2),
                `tissue mmol` = round(TBURDEN, 1),
                oxalosis = round(OXALOSIS, 3))
    bind_rows(f(d(), "sidebar arm"), f(cmp(), "comparator")) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })
}

shinyApp(ui, server)
