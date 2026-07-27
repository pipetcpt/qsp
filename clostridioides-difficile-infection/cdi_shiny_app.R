# =============================================================================
#  Clostridioides difficile Infection — QSP Shiny dashboard
# =============================================================================
#
#  Interactive front end for cdi_mrgsolve_model.R (61 ODE compartments).
#
#  The app is organised around the question the model exists to answer: not
#  "does this drug kill C. difficile" (they all do) but "what ecological state
#  does the patient come out of therapy in, and does that state permit the
#  spore reservoir to re-fire".  Every tab is a different view of that race.
#
#  Tabs
#  ----
#   1. Patient & regimen ...... covariates, the index antibiotic, the exposure,
#                               and the regimen builder; shows the resulting
#                               clinical course and the endpoint summary card
#   2. Ecology ................ the six guilds, Shannon diversity, and the
#                               anaerobe-niche occupancy that IS colonization
#                               resistance
#   3. Bile acids & niche ..... conjugated / free / secondary bile acids, the
#                               germination and outgrowth modifiers they produce,
#                               plus the sialic acid and Stickland pools
#   4. Pathogen & toxin ....... spore / vegetative / mucosal / reservoir burdens,
#                               PaLoc induction, luminal and mucosal toxin
#   5. Mucosa & immunity ...... colonocytes, stem pool, tight junctions, mucus,
#                               permeability, cytokines, neutrophils, membranes
#   6. Clinical endpoints ..... stools/day, WBC, creatinine, albumin, IDSA
#                               severity banding, cure / relapse markers
#   7. Drug exposure .......... faecal and plasma concentrations for every agent,
#                               target engagement, and collateral guild kill
#   8. Regimen comparison ..... all built-in scenarios side by side, with the
#                               end-of-therapy ecology table and the mapped
#                               8-week recurrence probability
#   9. Model & references ..... structure, calibration anchors, caveats
#
#  Run:  shiny::runApp("cdi_shiny_app.R")
#        (expects cdi_mrgsolve_model.R in the same directory)
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

# ---- load the model without triggering its own scenario report ---------------
# CDI_SOURCE_ONLY suppresses the report block at the bottom of the model file,
# so sourcing it here costs nothing but the compile.
CDI_SOURCE_ONLY <- TRUE
find_model <- function() {
  cands <- c("cdi_mrgsolve_model.R",
             file.path(getwd(), "cdi_mrgsolve_model.R"),
             file.path(dirname(sys.frame(1)$ofile %||% "."), "cdi_mrgsolve_model.R"))
  hit <- cands[file.exists(cands)]
  if (!length(hit)) stop("cdi_mrgsolve_model.R not found next to the app")
  hit[1]
}
# local = TRUE keeps the model objects in the app environment (reachable from
# the server function by lexical scoping) AND lets the sourced file see
# CDI_SOURCE_ONLY, which is what suppresses its 18-scenario report.
source(find_model(), local = TRUE)

# ---- shared plot styling ----------------------------------------------------
PAL <- c(
  SBA = "#1e7b3c", BUT = "#4aa564", BAC = "#8a7418", BIF = "#146b6b",
  ENT = "#c0392b", ENC = "#e67e22",
  TCA = "#d9a377", CA = "#c9ad55", CDCA = "#8a7418", DCA = "#1e7b3c", LCA = "#4aa564",
  veg = "#9b2222", muc = "#c0392b", spl = "#d98f8f", spb = "#5b2a8d",
  a   = "#e67e22", b = "#9b2222",
  drug1 = "#1b4f8f", drug2 = "#2a7fbf", drug3 = "#6b2185", drug4 = "#146b6b"
)

theme_cdi <- function() {
  theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 10, colour = "grey35"),
          legend.position = "bottom", legend.title = element_blank(),
          strip.text = element_text(face = "bold", size = 10))
}

# Therapy band.  On a log10 axis an -Inf/Inf rect makes scales::log_breaks
# take log(-Inf), so pass a finite y range when the panel is log-scaled.
rx_band <- function(t0, t1, ylim = NULL) {
  if (is.na(t0) || is.na(t1)) return(NULL)
  ymin <- if (is.null(ylim)) -Inf else ylim[1]
  ymax <- if (is.null(ylim)) Inf  else ylim[2]
  annotate("rect", xmin = t0, xmax = t1, ymin = ymin, ymax = ymax,
           fill = "#1b4f8f", alpha = 0.06)
}

long_plot <- function(d, vars, labels, title, sub, ylab, cols = NULL,
                      logy = FALSE, hline = NULL, rx = NULL) {
  dd <- d %>% select(time, all_of(vars)) %>%
    pivot_longer(-time, names_to = "v", values_to = "y") %>%
    mutate(v = factor(v, levels = vars, labels = labels))
  p <- ggplot(dd, aes(time, y, colour = v)) +
    (if (!is.null(rx)) rx_band(rx[1], rx[2]) else NULL) +
    geom_line(linewidth = 0.85) +
    labs(title = title, subtitle = sub, x = "Day", y = ylab) +
    theme_cdi()
  if (!is.null(cols)) p <- p + scale_colour_manual(values = unname(cols))
  # the stiff solver can leave states a hair below zero; floor before log10 so
  # the axis transform never produces NaN
  if (logy) {
    dd$y <- pmax(dd$y, 1e-12)
    yl <- range(dd$y, finite = TRUE)
    p <- ggplot(dd, aes(time, y, colour = v)) +
      (if (!is.null(rx)) rx_band(rx[1], rx[2], yl) else NULL) +
      geom_line(linewidth = 0.85) +
      labs(title = title, subtitle = sub, x = "Day", y = ylab) +
      theme_cdi() + scale_y_log10()
    if (!is.null(cols)) p <- p + scale_colour_manual(values = unname(cols))
  }
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = 2,
                                           colour = "grey45", linewidth = 0.4)
  p
}

# =============================================================================
#  UI
# =============================================================================
ui <- fluidPage(
  title = "C. difficile infection — QSP dashboard",
  tags$head(tags$style(HTML("
    body { background:#fbfbfa; }
    h4 { margin-top:4px; }
    .wellish { background:#fff; border:1px solid #e6e6e2; border-radius:8px;
               padding:12px 14px; margin-bottom:12px; }
    .kpi { display:inline-block; min-width:132px; margin:4px 10px 4px 0;
           padding:8px 10px; border-radius:8px; background:#fff;
           border:1px solid #e6e6e2; }
    .kpi .lab { font-size:10.5px; color:#666; text-transform:uppercase;
                letter-spacing:.04em; }
    .kpi .val { font-size:19px; font-weight:700; }
    .sev0 { color:#1e7b3c } .sev1 { color:#b8860b } .sev2 { color:#a52020 }
    .note { font-size:11.5px; color:#555; }
  "))),

  titlePanel(div(
    h3(HTML("<i>Clostridioides difficile</i> infection — QSP dashboard"),
       style = "margin-bottom:2px"),
    div(class = "note",
        "61-compartment mrgsolve model: microbial ecology × bile-acid",
        "germination switch × PaLoc toxin output × mucosal barrier.",
        "Shaded band = days on anti-C. difficile therapy.")
  )),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(class = "wellish",
        h4("Patient"),
        sliderInput("wt", "Weight (kg)", 40, 140, 70, 5),
        sliderInput("igg0", "Baseline anti-TcdB IgG (1 = median)", 0.1, 4.0, 1.0, 0.1),
        sliderInput("immcomp", "Humoral response capacity", 0.1, 1.6, 1.0, 0.05),
        sliderInput("fiber", "Fermentable fibre intake", 0.2, 1.4, 1.0, 0.1),
        checkboxInput("rt027", "Ribotype 027 / NAP1 (high toxin, high sporulation)", FALSE),
        checkboxInput("vre", "Vancomycin-resistant Enterococcus", FALSE),
        checkboxInput("ppi", "On a proton-pump inhibitor", FALSE)
      ),
      div(class = "wellish",
        h4("Index (precipitating) antibiotic"),
        sliderInput("abx_days", "Duration (days from day 0)", 0, 45, 7, 1),
        sliderInput("abx_amt", "Relative dysbiotic potency", 0.2, 2.0, 1.0, 0.1),
        sliderInput("inoc", "Ingested spore inoculum (×10⁶ CFU/g)",
                    0, 0.6, 0.05, 0.01)
      ),
      div(class = "wellish",
        h4("Anti-C. difficile therapy"),
        sliderInput("tstart", "Day therapy starts", 5, 20, 9, 1),
        selectInput("abx", "Antibacterial",
          c("none", "vancomycin 125 mg q6h", "vancomycin taper/pulse",
            "fidaxomicin 200 mg q12h", "fidaxomicin extended-pulsed (EXTEND)",
            "metronidazole 500 mg q8h", "ridinilazole 200 mg q12h"),
          selected = "vancomycin 125 mg q6h"),
        sliderInput("rxdays", "Course length (days)", 5, 20, 10, 1),
        checkboxInput("rifax", "Rifaximin 400 mg q8h × 14 d chaser", FALSE),
        checkboxInput("tige", "IV tigecycline salvage", FALSE)
      ),
      div(class = "wellish",
        h4("Recurrence prevention"),
        checkboxInput("bez", "Bezlotoxumab 10 mg/kg IV, single dose", FALSE),
        selectInput("lbp", "Microbiome restoration",
          c("none", "FMT", "SER-109 (oral Firmicutes spores ×3 d)",
            "RBX2660 (rectal suspension)"), selected = "none"),
        sliderInput("lbp_day", "Days after last antibacterial dose", 0.5, 10, 1.5, 0.5)
      ),
      div(class = "wellish",
        sliderInput("tend", "Follow-up horizon (days)", 40, 180, 90, 10),
        actionButton("run", "Simulate", class = "btn-primary btn-block")
      )
    ),

    mainPanel(
      width = 9,
      uiOutput("kpis"),
      tabsetPanel(
        id = "tabs", type = "tabs",
        tabPanel("1 · Course",
          plotOutput("p_course", height = 320),
          plotOutput("p_course2", height = 300),
          div(class = "note", textOutput("course_txt"))),
        tabPanel("2 · Ecology",
          plotOutput("p_guild", height = 320),
          plotOutput("p_div", height = 280),
          div(class = "note",
              "MB_SBA is the bai⁺ 7α-dehydroxylating Clostridia guild.",
              "Its relative abundance at the end of therapy is the single best",
              "predictor of what happens next, because it alone makes the",
              "secondary bile acids that suppress germination and outgrowth.")),
        tabPanel("3 · Bile acids & niche",
          plotOutput("p_ba", height = 320),
          plotOutput("p_niche", height = 300),
          div(class = "note",
              "Vancomycin spares Bacteroidetes, so deconjugation continues and",
              "the conjugated pool rises only modestly — the permissive change",
              "is the disappearance of deoxycholate and lithocholate, plus a",
              "nutrient niche that Bacteroidetes actively open by liberating",
              "sialic acid from mucin.")),
        tabPanel("4 · Pathogen & toxin",
          plotOutput("p_cd", height = 320),
          plotOutput("p_tox", height = 300),
          div(class = "note",
              "Toxin lags the bacterial peak: CodY/CcpA repress the PaLoc while",
              "Stickland amino acids are abundant, so induction rises as the",
              "population exhausts its own niche.")),
        tabPanel("5 · Mucosa & immunity",
          plotOutput("p_epi", height = 320),
          plotOutput("p_imm", height = 300)),
        tabPanel("6 · Clinical endpoints",
          plotOutput("p_clin", height = 340),
          plotOutput("p_sev", height = 260),
          div(class = "note",
              "IDSA/SHEA severe CDI = WBC ≥ 15 × 10⁹/L or creatinine ≥ 1.5 mg/dL.")),
        tabPanel("7 · Drug exposure",
          plotOutput("p_pk", height = 330),
          plotOutput("p_engage", height = 300),
          div(class = "note",
              "Note metronidazole: it reaches the lumen almost entirely by",
              "secretion across inflamed mucosa, so its faecal concentration",
              "falls as the patient improves.")),
        tabPanel("8 · Regimen comparison",
          div(class = "note",
              "Runs every built-in scenario on the CURRENT patient covariates."),
          actionButton("runall", "Run all scenarios", class = "btn-default"),
          plotOutput("p_cmp", height = 380),
          tableOutput("t_cmp")),
        tabPanel("9 · Model & references",
          htmlOutput("about"))
      )
    )
  )
)

# =============================================================================
#  SERVER
# =============================================================================
server <- function(input, output, session) {

  # ---- assemble the event table from the UI ---------------------------------
  build <- reactive({
    ts <- input$tstart
    ev_list <- list()
    if (input$abx_days > 0)
      ev_list[[length(ev_list) + 1]] <-
        ev(time = 0, amt = input$abx_amt, cmt = "ABX_D", ii = 1,
           addl = input$abx_days - 1)
    if (input$inoc > 0)
      ev_list[[length(ev_list) + 1]] <-
        ev(time = 2, amt = input$inoc * (1 + 1.5 * as.numeric(input$ppi)),
           cmt = "CD_SPORE_L")

    rx0 <- rx1 <- NA_real_
    a <- input$abx
    if (a == "vancomycin 125 mg q6h") {
      ev_list[[length(ev_list) + 1]] <- ev_vanco(ts, input$rxdays)
      rx0 <- ts; rx1 <- ts + input$rxdays
    } else if (a == "vancomycin taper/pulse") {
      ev_list[[length(ev_list) + 1]] <- ev_vanco_taper(ts)
      rx0 <- ts; rx1 <- ts + 55
    } else if (a == "fidaxomicin 200 mg q12h") {
      ev_list[[length(ev_list) + 1]] <- ev_fidaxo(ts, input$rxdays)
      rx0 <- ts; rx1 <- ts + input$rxdays
    } else if (a == "fidaxomicin extended-pulsed (EXTEND)") {
      ev_list[[length(ev_list) + 1]] <- ev_fidaxo_extend(ts)
      rx0 <- ts; rx1 <- ts + 25
    } else if (a == "metronidazole 500 mg q8h") {
      ev_list[[length(ev_list) + 1]] <- ev_metro(ts, input$rxdays)
      rx0 <- ts; rx1 <- ts + input$rxdays
    } else if (a == "ridinilazole 200 mg q12h") {
      ev_list[[length(ev_list) + 1]] <- ev_ridini(ts, input$rxdays)
      rx0 <- ts; rx1 <- ts + input$rxdays
    }
    rxend <- if (is.na(rx1)) ts else rx1

    if (isTRUE(input$rifax))
      ev_list[[length(ev_list) + 1]] <- ev_rifax(rxend, 14)
    if (isTRUE(input$bez))
      ev_list[[length(ev_list) + 1]] <- ev_bezlo(ts + 1, input$wt)

    p <- LBP_NONE
    if (input$lbp == "FMT") {
      p <- LBP_FMT
      ev_list[[length(ev_list) + 1]] <- ev_lbp(rxend + input$lbp_day)
    } else if (input$lbp == "SER-109 (oral Firmicutes spores ×3 d)") {
      p <- LBP_SER
      ev_list[[length(ev_list) + 1]] <-
        ev_lbp(rxend + input$lbp_day, amt = 0.40, n = 3)
    } else if (input$lbp == "RBX2660 (rectal suspension)") {
      p <- LBP_RBX
      ev_list[[length(ev_list) + 1]] <- ev_lbp(rxend + input$lbp_day, amt = 1.0)
    }

    if (!length(ev_list)) ev_list[[1]] <- ev(time = 0, amt = 0, cmt = "LBP")
    e <- Reduce(c, ev_list)

    pars <- c(p, list(
      WT = input$wt, IGG0 = input$igg0, IMMCOMP = input$immcomp,
      FIBER = input$fiber, RT027 = as.numeric(input$rt027),
      VRE = as.numeric(input$vre), PPI = as.numeric(input$ppi),
      TIGE = as.numeric(input$tige)))

    list(ev = e, par = pars, rx = c(rx0, rx1), ts = ts, rxend = rxend)
  })

  sim <- eventReactive(input$run, {
    b <- build()
    d <- mrgsim_e(param(mod, b$par), b$ev, end = input$tend, delta = 0.25,
                  atol = 1e-8, rtol = 1e-6, maxsteps = 200000) %>% as_tibble()
    list(d = d, b = b, e = endpoints(d, tstart = b$ts, t_rx_end = rx_end(b$ev)))
  }, ignoreNULL = FALSE)

  RX <- reactive(sim()$b$rx)

  # ---- KPI strip ------------------------------------------------------------
  output$kpis <- renderUI({
    e <- sim()$e
    sevcls <- if (e$peak_WBC >= 20 || e$max_CRE >= 2) "sev2" else
              if (e$peak_WBC >= 15 || e$max_CRE >= 1.5) "sev1" else "sev0"
    k <- function(lab, val, cls = "") div(class = "kpi",
      div(class = "lab", lab), div(class = paste("val", cls), val))
    # 8-week recurrence probability, using the mapping calibrated in the model
    p <- rri_to_prob(e$RRI_rx, -2.236, 1.994)
    div(
      k("Peak stools/day", sprintf("%.1f", e$peak_stool)),
      k("Peak WBC", sprintf("%.1f", e$peak_WBC), sevcls),
      k("Max creatinine", sprintf("%.2f", e$max_CRE), sevcls),
      k("Nadir albumin", sprintf("%.2f", e$min_ALB)),
      k("Time to resolution",
        if (is.na(e$ttrod)) "not reached" else sprintf("%.0f d", e$ttrod)),
      k("Relapse in model",
        if (e$recurred) sprintf("day %.0f", e$recur_day) else "none"),
      k("MB_SBA at end of rx", sprintf("%.1f%% of normal", 100 * e$sba_frac_rx)),
      k("Deoxycholate at end of rx", sprintf("%.0f µM", e$dca_rx)),
      k("Mapped 8-wk recurrence", sprintf("%.0f%%", 100 * p),
        if (p > 0.22) "sev2" else if (p > 0.14) "sev1" else "sev0")
    )
  })

  # ---- 1. course ------------------------------------------------------------
  output$p_course <- renderPlot({
    d <- sim()$d
    long_plot(d, c("STOOL", "WBC", "CRE"),
      c("Unformed stools/day", "WBC (10⁹/L)", "Creatinine (mg/dL)"),
      "Clinical course", "dashed lines: cure threshold (3 stools/day) and IDSA severity cut-offs",
      "value", rx = RX()) +
      geom_hline(yintercept = c(3, 15, 1.5), linetype = 3, colour = "grey55",
                 linewidth = 0.35)
  })
  output$p_course2 <- renderPlot({
    d <- sim()$d
    long_plot(d, c("LGVEG", "LGSPB", "TOXLo"),
      c("log10 luminal C. difficile (CFU/g)",
        "log10 mucosal spore reservoir (CFU/g)",
        "Mucosal toxin load (ng/mL-equiv)"),
      "What drives the course",
      "the reservoir outlives the drug; the guild decides whether it can re-fire",
      "value", rx = RX())
  })
  output$course_txt <- renderText({
    e <- sim()$e
    paste0(
      "Symptomatic: ", e$symptomatic,
      if (!is.na(e$onset_day)) paste0(" (onset day ", e$onset_day, ")") else "",
      ". Days with ≥3 unformed stools: ", e$sick_days,
      ". Days meeting IDSA severity: ", e$severe_days,
      ". Nadir viable colonocyte fraction: ", sprintf("%.2f", e$min_EPI),
      ". Recurrence risk index at assessment (day ", e$assess_day, "): ",
      sprintf("%.3f", e$RRI_rx), ".")
  })

  # ---- 2. ecology -----------------------------------------------------------
  output$p_guild <- renderPlot({
    d <- sim()$d
    long_plot(d, c("MB_SBA", "MB_BUT", "MB_BAC", "MB_BIF", "MB_ENT", "MB_ENC"),
      c("bai⁺ Clostridia (MB_SBA)", "Butyrogens (MB_BUT)",
        "Bacteroidetes (MB_BAC)", "Bifidobacterium (MB_BIF)",
        "Enterobacteriaceae (MB_ENT)", "Enterococcus (MB_ENC)"),
      "Microbiota guilds", "relative abundance, log scale",
      "relative abundance", cols = PAL[1:6], logy = TRUE, rx = RX())
  })
  output$p_div <- renderPlot({
    d <- sim()$d
    long_plot(d, c("SHAN", "fANAo"),
      c("Shannon diversity (modelled guilds)",
        "Anaerobe-niche occupancy (fraction of healthy)"),
      "Colonization resistance", "the quantity that sets the C. difficile carrying capacity",
      "value", rx = RX())
  })

  # ---- 3. bile acids & niche ----------------------------------------------
  output$p_ba <- renderPlot({
    d <- sim()$d
    long_plot(d, c("BA_TCA", "BA_CA", "BA_CDCA", "BA_DCA", "BA_LCA"),
      c("Conjugated primary (germinant)", "Cholate", "Chenodeoxycholate (inhibitor)",
        "Deoxycholate", "Lithocholate"),
      "Colonic bile acids", "µM in colonic water",
      "µM", cols = PAL[c("TCA","CA","CDCA","DCA","LCA")], rx = RX())
  })
  output$p_niche <- renderPlot({
    d <- sim()$d
    long_plot(d, c("NUT_SIA", "NUT_AA", "SCFA_BUT"),
      c("Free sialic acid + succinate (mM)", "Stickland amino acids (mM)",
        "Butyrate (mM)"),
      "Luminal nutrient niche",
      "the niche opens because the fermenters that normally consume it are gone",
      "mM", rx = RX())
  })

  # ---- 4. pathogen & toxin -----------------------------------------------
  output$p_cd <- renderPlot({
    d <- sim()$d
    long_plot(d, c("LGSPL", "LGVEG", "LGMUC", "LGSPB"),
      c("Luminal spores", "Luminal vegetative", "Mucosa-adherent",
        "Mucosal spore reservoir"),
      "C. difficile compartments", "log10 CFU/g faeces",
      "log10 CFU/g", cols = PAL[c("spl","veg","muc","spb")], rx = RX())
  })
  output$p_tox <- renderPlot({
    d <- sim()$d
    long_plot(d, c("TCDA", "TCDB", "TCDA_MUC", "TCDB_MUC", "TOX_CPLX"),
      c("Luminal TcdA", "Luminal TcdB", "Mucosal TcdA", "Mucosal TcdB",
        "TcdB·bezlotoxumab complex"),
      "Toxins", "ng/mL; the complexed pool is what a neutralising antibody removes",
      "ng/mL", rx = RX())
  })

  # ---- 5. mucosa & immunity ---------------------------------------------
  output$p_epi <- renderPlot({
    d <- sim()$d
    long_plot(d, c("EPI", "EPI_SC", "EPI_TJ", "EPI_MUCUS", "EPI_PERM"),
      c("Viable colonocytes", "Lgr5⁺ stem/progenitor pool",
        "Tight-junction integrity", "MUC2 mucus layer",
        "Paracellular permeability index"),
      "Colonic epithelium",
      "TcdB blocks Frizzled, so the stem pool that has to repair the damage is itself a target",
      "fraction of normal / index", rx = RX())
  })
  output$p_imm <- renderPlot({
    d <- sim()$d
    long_plot(d, c("IM_IL8", "IM_IL1B", "IM_TNF", "IM_NEUT", "IM_IL22",
                   "IM_PSM", "AB_IGG"),
      c("IL-8/CXCL8", "IL-1β (pyrin/NLRP3)", "TNF-α/IL-6",
        "Mucosal neutrophils", "IL-22 (protective)", "Pseudomembranes",
        "Anti-TcdB IgG"),
      "Immune response", "normalised units (0 at health)", "units", rx = RX())
  })

  # ---- 6. clinical endpoints ---------------------------------------------
  output$p_clin <- renderPlot({
    d <- sim()$d
    long_plot(d, c("STOOL", "WBC", "ALB", "CRE"),
      c("Unformed stools/day", "WBC (10⁹/L)", "Albumin (g/dL)",
        "Creatinine (mg/dL)"),
      "Clinical endpoints",
      "albumin falls through permeability-driven protein loss, creatinine through volume depletion",
      "value", rx = RX()) + facet_wrap(~v, scales = "free_y")
  })
  output$p_sev <- renderPlot({
    d <- sim()$d
    long_plot(d, c("SEVERE", "FULMIN", "SYMPT"),
      c("IDSA/SHEA severe", "Fulminant", "Symptomatic (≥3 stools/day)"),
      "Severity banding over time", "1 = criterion met", "flag", rx = RX())
  })

  # ---- 7. drug exposure ---------------------------------------------------
  output$p_pk <- renderPlot({
    d <- sim()$d
    long_plot(d, c("CVANo", "CFDXo", "COPo", "CMTZFo", "CRIDo", "CRFXo"),
      c("Vancomycin (faecal)", "Fidaxomicin (faecal)", "OP-1118 (faecal)",
        "Metronidazole (faecal)", "Ridinilazole (faecal)", "Rifaximin (faecal)"),
      "Colonic drug exposure", "µg/g faeces", "µg/g", rx = RX())
  })
  output$p_engage <- renderPlot({
    d <- sim()$d
    long_plot(d, c("CMTZSo", "CBEZo"),
      c("Metronidazole, plasma (mg/L)", "Bezlotoxumab, plasma (mg/L)"),
      "Systemic exposure", "the two agents in the model that act systemically",
      "mg/L", rx = RX()) + facet_wrap(~v, scales = "free_y")
  })

  # ---- 8. regimen comparison ---------------------------------------------
  cmp <- eventReactive(input$runall, {
    b <- build()
    covs <- b$par[intersect(names(b$par),
      c("WT","IGG0","IMMCOMP","FIBER","RT027","VRE","PPI","TIGE"))]
    withProgress(message = "Running scenarios", value = 0, {
      out <- lapply(seq_along(scen), function(i) {
        nm <- names(scen)[i]
        incProgress(1 / length(scen), detail = nm)
        s <- scen[[nm]]
        m <- param(mod, c(s$par, covs))
        d <- mrgsim_e(m, s$ev, end = input$tend, delta = 0.25,
                      atol = 1e-8, rtol = 1e-6, maxsteps = 200000) %>% as_tibble()
        e <- endpoints(d, t_rx_end = rx_end(s$ev)); e$scenario <- nm
        list(d = d %>% mutate(scenario = nm), e = e)
      })
    })
    list(d = bind_rows(lapply(out, `[[`, "d")),
         e = bind_rows(lapply(out, `[[`, "e")))
  })

  output$p_cmp <- renderPlot({
    cd <- cmp()$d
    ggplot(cd, aes(time, STOOL, colour = scenario)) +
      geom_hline(yintercept = 3, linetype = 2, colour = "grey55") +
      geom_line(linewidth = 0.7) +
      labs(title = "Unformed stools per day, all built-in scenarios",
           subtitle = "dashed line = clinical-cure threshold; a second crossing is a modelled relapse",
           x = "Day", y = "Unformed stools / day") +
      theme_cdi() + theme(legend.position = "right", legend.text = element_text(size = 8))
  })

  output$t_cmp <- renderTable({
    e <- cmp()$e
    lg <- try(fit_rri_logistic(e), silent = TRUE)
    b0 <- if (inherits(lg, "try-error")) -2.236 else lg$b0
    b1 <- if (inherits(lg, "try-error")) 1.994 else lg$b1
    e %>% transmute(
      Scenario      = scenario,
      `Peak stools` = round(peak_stool, 1),
      `Peak WBC`    = round(peak_WBC, 1),
      `TTROD (d)`   = ifelse(is.na(ttrod), "—", sprintf("%.0f", ttrod)),
      `Relapse`     = ifelse(recurred, sprintf("day %.0f", recur_day), "none"),
      `MB_SBA end rx (% normal)` = round(100 * sba_frac_rx, 1),
      `DCA end rx (µM)`     = round(dca_rx),
      `Reservoir (log10 CFU/g)`  = round(spb_rx, 2),
      RRI           = round(RRI_rx, 3),
      `8-wk recurrence` = sprintf("%.0f%%", 100 * rri_to_prob(RRI_rx, b0, b1))
    )
  }, striped = TRUE, spacing = "xs")

  # ---- 9. about ------------------------------------------------------------
  output$about <- renderUI(HTML('
<h4>What this model is</h4>
<p>A 61-compartment mrgsolve QSP model of <i>Clostridioides difficile</i>
infection built around the recurrence loop rather than around the organism.
Killing <i>C. difficile</i> is the easy part: oral vancomycin sterilises the
stool within days and cures roughly four patients in five &mdash; and then one in
four relapses, because the drug that cleared the organism also held down the
microbial guild that provides colonization resistance while a drug-proof spore
reservoir waited.</p>

<h4>Three coupled layers</h4>
<ul>
<li><b>Ecology</b> &mdash; six guilds under drug-specific kill, and the two
currencies through which they suppress <i>C. difficile</i>: secondary bile acids
(bai&#8314; 7&alpha;-dehydroxylation) and the luminal nutrient niche. Crucially the
nutrient-competition guild set <i>excludes</i> Bacteroidetes, which liberate
sialic acid rather than consuming it &mdash; which is why vancomycin, a
Bacteroides-sparing drug, still leaves the niche wide open.</li>
<li><b>Pathogen</b> &mdash; spore &rarr; CspC-mediated germination &rarr; vegetative
outgrowth &rarr; sporulation, with a mucosa/biofilm spore reservoir that is the
actual seed of relapse.</li>
<li><b>Host</b> &mdash; PaLoc-regulated toxin output (CodY/CcpA nutrient
repression makes toxin lag the bacterial peak), Rho-GTPase glucosylation,
Frizzled blockade of stem-cell renewal, pyrin-inflammasome IL-1&beta;,
neutrophilic pseudomembranes, and the endpoints trials actually measure.</li>
</ul>

<h4>Why the arms differ in the model</h4>
<p>Two therapeutic quantities are structurally separate:
<code>KILLCD</code> (how fast the organism dies) and the collateral kill of
<code>MB_SBA</code> (how hard the restoring guild is hit). Fidaxomicin and
ridinilazole score well on both and do not relapse; vancomycin and metronidazole
cure the episode and buy the relapse; FMT, SER-109 and RBX2660 abort it by
re-seeding the guild; bezlotoxumab acts on neither but neutralises the toxin
through the vulnerable window.</p>

<h4>Calibration</h4>
<p>The drug-free baseline is a <i>solved</i> steady state &mdash; $MAIN inverts the
bile-acid cascade algebraically, including the closed-form Michaelis-Menten root
for cholate and chenodeoxycholate &mdash; and drifts by 0.000% over 90 simulated
days. Faecal and plasma exposures reproduce published values (vancomycin
~900 &micro;g/g; fidaxomicin + OP-1118 ~900 &micro;g/g; metronidazole ~10 &micro;g/g faecal and
~14 mg/L plasma; bezlotoxumab C<sub>max</sub> ~230 mg/L). Time to resolution of
diarrhoea is 3&ndash;4 days on therapy against a median of 2&ndash;4 days in trials.
The mechanistic recurrence risk index maps onto observed 8-week recurrence
across six independent phase 3 anchors (Louie 2011, Cornely 2012, Johnson 2014,
Wilcox 2017, van Nood 2013, Feuerstadt 2022) with R&sup2; &asymp; 0.73.</p>

<h4>Limits</h4>
<p>Semi-quantitative and deterministic. A deterministic trajectory either
relapses or it does not; a real cohort splits, which is why the arm-level
recurrence rates come from the explicit, auditable RRI&rarr;probability mapping
rather than from counting relapses in single runs. Parameters are
literature-anchored order-of-magnitude estimates, not a fitted population model.
The FMT and SER-109 anchors come from recurrent-CDI populations with a higher
baseline hazard than the first-episode arms, so they are conservative.</p>

<p><b>Not for clinical decision-making, prescribing, or regulatory
submission.</b> Full source list with 129 PubMed-verified citations in
<code>cdi_references.md</code>.</p>'))
}

shinyApp(ui, server)
