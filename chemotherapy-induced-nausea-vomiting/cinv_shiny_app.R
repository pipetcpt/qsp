## =============================================================================
##  cinv_shiny_app.R
##  Chemotherapy-Induced Nausea and Vomiting (CINV) — interactive QSP dashboard
##
##  The app is built around ONE claim and every tab exists to let a user try to
##  break it:
##
##      Every antiemetic subtracts ONE ADDITIVE TERM from the emetic drive, and
##      every endpoint is exp(-integral of a threshold hazard).  Additivity is
##      therefore exact ONLY in the drive, and is destroyed on the way to the
##      endpoint by the saturating hazard - which is where all reported
##      "synergy" between antiemetics actually comes from.
##
##  Tab 4 ("Additive terms") shows the sum decomposed, so you can watch a drug
##  delete its own term and leave the others untouched.  Tab 9 ("log-additivity")
##  tests the claim numerically.  Tab 10 shows the corollary that the absolute
##  benefit of the fourth agent is NON-MONOTONIC in baseline risk.
##
##  RUN:  shiny::runApp("cinv_shiny_app.R")
##        (expects cinv_mrgsolve_model.R in the same directory)
## =============================================================================

library(shiny)
library(mrgsolve)
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))

source("cinv_mrgsolve_model.R", local = FALSE)

theme_set(theme_bw(base_size = 12) +
          theme(panel.grid.minor = element_blank(),
                strip.background = element_rect(fill = "grey92", colour = NA),
                legend.position = "bottom"))

PAL <- c("#2c5aa0", "#c0392b", "#1c7a4a", "#8a5cc0", "#c07a10",
         "#1c6a7a", "#a0306a", "#5a5a5a")

DRUGS <- c("Ondansetron 8 mg IV" = "ond", "Ondansetron 32 mg IV" = "ond32",
           "Granisetron 1 mg IV" = "gra", "Granisetron ER-SC 10 mg" = "graER",
           "Palonosetron 0.25 mg IV" = "pal",
           "Aprepitant 125/80/80 PO" = "apr", "Fosaprepitant 150 mg IV" = "fos",
           "Netupitant 300 mg PO" = "net", "Rolapitant 180 mg PO" = "rol",
           "Dexamethasone (full 20 mg d1)" = "dexF",
           "Dexamethasone (NK1-corrected 12 mg d1)" = "dexC",
           "Olanzapine 10 mg" = "ola10", "Olanzapine 5 mg" = "ola5",
           "Metoclopramide 20 mg" = "mcp", "Lorazepam 1 mg" = "lor",
           "Dronabinol 5 mg" = "dro")

build_events <- function(sel) {
  e <- list()
  add <- function(x) e[[length(e) + 1]] <<- x
  if ("ond"   %in% sel) add(dose_ond_iv(8))
  if ("ond32" %in% sel) add(dose_ond_iv(32))
  if ("gra"   %in% sel) add(dose_gra_iv(1))
  if ("graER" %in% sel) add(dose_graER(10))
  if ("pal"   %in% sel) add(dose_pal_iv(0.25))
  if ("apr"   %in% sel) add(APR_STD())
  if ("fos"   %in% sel) add(dose_fos(150))
  if ("net"   %in% sel) add(dose_net(300))
  if ("rol"   %in% sel) add(dose_rol(180))
  if ("dexF"  %in% sel) add(DEX_STD())
  if ("dexC"  %in% sel) add(DEX_NK1())
  if ("ola10" %in% sel) add(dose_ola(10))
  if ("ola5"  %in% sel) add(dose_ola(5))
  if ("mcp"   %in% sel) add(dose_mcp(20, c(-0.5, 6, 12, 24, 36, 48, 60, 72)))
  if ("lor"   %in% sel) add(dose_lor(1, c(-1, 12, 24, 36, 48)))
  if ("dro"   %in% sel) add(dose_dro(5, c(-1, 6, 12, 24, 36, 48, 60, 72)))
  if (!length(e)) return(NULL)
  do.call(rbind, e)
}

simulate_custom <- function(sel, chemo, covars, end = 120, cycles = 1,
                            par_over = list()) {
  em <- EMETO[[chemo]]
  pars <- c(list(EMETO_P = em[["P"]], EMETO_C = em[["C"]], CISRATE = em[["rate"]]),
            covars, par_over)
  m2 <- do.call(param, c(list(mod), pars))
  ev <- build_events(sel)
  if (cycles > 1 && !is.null(ev)) {
    ev <- do.call(rbind, lapply(0:(cycles - 1), function(k) {
      x <- ev; x$time <- x$time + k * 504; x }))
  }
  tstart <- 0
  if (!is.null(ev)) {
    ev <- ev[order(ev$time), , drop = FALSE]
    ev$ID <- 1L
    ev <- ev[, c("ID", "time", "cmt", "amt", "evid")]
    tstart <- min(0, min(ev$time))
  }
  out <- if (is.null(ev)) {
    mrgsim(m2, start = tstart, end = end, delta = 0.25, atol = 1e-9, rtol = 1e-8)
  } else {
    mrgsim(m2, data = ev, start = tstart, end = end, delta = 0.25,
           atol = 1e-9, rtol = 1e-8, recsort = 3)
  }
  as.data.frame(out)
}

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("CINV QSP dashboard — an additive drive, a threshold hazard, an exponential endpoint"),
  tags$head(tags$style(HTML("
    .thesis {background:#f3eefb; border-left:5px solid #5a1c9a; padding:9px 13px;
             margin-bottom:11px; font-size:13px;}
    .warnbox{background:#fdf0e8; border-left:5px solid #a04a00; padding:8px 12px;
             font-size:12.5px;}
    .kpi {background:#eef4fb; border:1px solid #bcd8f2; border-radius:6px;
          padding:7px 10px; margin:3px; display:inline-block; min-width:132px;}
    .kpi b{font-size:17px; color:#1a3d6b;}
  "))),
  div(class = "thesis", HTML(
    "<b>The claim this dashboard exists to test.</b> Each antiemetic removes exactly
     one additive term <i>w<sub>j</sub>A<sub>j</sub></i> from the emetic drive; the
     endpoint is CR = E<sub>S</sub>[exp(&minus;S&middot;&int;&lambda;)]. Additivity is
     exact only in the <i>drive</i>; the saturating hazard destroys it before the
     endpoint. Add and remove agents in the sidebar and watch the decomposition in
     <b>Additive terms</b>, then check the arithmetic in <b>log-additivity</b>.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("chemo", "Chemotherapy regimen",
                  choices = c("Cisplatin (HEC)" = "cisplatin",
                              "Anthracycline + cyclophosphamide (HEC)" = "AC",
                              "Carboplatin AUC 5 (MEC)" = "carboplatin",
                              "Oxaliplatin (MEC)" = "oxaliplatin",
                              "Other moderately emetogenic" = "MEC",
                              "Low emetogenic" = "LEC"),
                  selected = "cisplatin"),
      checkboxGroupInput("drugs", "Antiemetic regimen",
                         choices = DRUGS, selected = c("ond", "dexF")),
      hr(),
      h5("Patient"),
      radioButtons("sex", "Sex", c("Female" = 1, "Male" = 0), inline = TRUE),
      sliderInput("age", "Age (years)", 20, 85, 52, step = 1),
      checkboxInput("alcohol", "Chronic heavy alcohol history (protective)", FALSE),
      checkboxInput("motion", "Motion sickness / pregnancy emesis history", FALSE),
      checkboxInput("prior", "CINV in a previous cycle", FALSE),
      selectInput("cyp2d6", "CYP2D6 phenotype",
                  c("Normal metaboliser" = 1, "Poor metaboliser" = 0.1,
                    "Ultrarapid metaboliser" = 2.5), selected = 1),
      hr(),
      sliderInput("cycles", "Cycles to simulate", 1, 6, 1, step = 1),
      sliderInput("omega", "Population susceptibility SD (omega)",
                  0.2, 1.6, 0.86, step = 0.02),
      helpText("omega is what makes CR(0-120) exceed CR(0-24)xCR(24-120).")
    ),

    mainPanel(
      width = 9,
      uiOutput("kpis"),
      tabsetPanel(
        id = "tabs", type = "tabs",
        tabPanel("1 · Patient & endpoints", br(),
                 h4("Complete response and no-nausea, by window"),
                 plotOutput("p_endpoints", height = "300px"),
                 h4("Endpoint table"), tableOutput("t_endpoints"),
                 div(class = "warnbox", textOutput("risknote"))),
        tabPanel("2 · Antiemetic PK", br(),
                 h4("Plasma concentrations of the selected agents"),
                 plotOutput("p_pk", height = "360px"),
                 helpText("Note the four-order-of-magnitude spread in half-life:
                           ondansetron 4 h, palonosetron 40 h, netupitant 88 h,
                           rolapitant 180 h. A single rolapitant dose covers the
                           whole 120 h window; ondansetron does not reach day 2.")),
        tabPanel("3 · Receptor occupancy", br(),
                 h4("5-HT3 and NK1 receptor pools — who is bound to what"),
                 plotOutput("p_occ", height = "400px"),
                 helpText("The 5-HT3 pool is modelled KINETICALLY because
                           palonosetron's off-rate (0.10 /h) and internalisation
                           are the whole point; equilibrium occupancy cannot
                           represent them.")),
        tabPanel("4 · Additive terms  ★", br(),
                 div(class = "thesis", HTML("<b>This is the thesis, drawn.</b>
                   The emetic drive is a SUM. Each agent deletes its own term and
                   leaves the others exactly where they were. The area postrema
                   term (WAP&middot;FG) is downstream of no receptor, so nothing in
                   the formulary can remove it — that is why no regimen reaches
                   100% control.")),
                 plotOutput("p_terms", height = "340px"),
                 h4("Drive, its tonic floor, and the excess that drives the hazard"),
                 plotOutput("p_drive", height = "260px")),
        tabPanel("5 · Acute vs delayed", br(),
                 h4("Two terms, two time courses — the phases are an OUTPUT"),
                 plotOutput("p_phase", height = "380px"),
                 helpText("Peripheral 5-HT is a BURST (sigmoidal in mucosal
                           platinum, peak 4-6 h). Substance P is a transcriptional
                           RAMP peaking near 40 h. A 5-HT3 antagonist deletes the
                           first and cannot touch the second.")),
        tabPanel("6 · Nausea vs vomiting", br(),
                 h4("Two different weighted sums, two different endpoints"),
                 plotOutput("p_nausea", height = "360px"),
                 h4("Cortical-sum composition at peak"),
                 plotOutput("p_ctxbar", height = "260px"),
                 helpText("Olanzapine is the only agent that subtracts four
                           cortical terms (D2, 5-HT2A, H1, M1) at once.")),
        tabPanel("7 · Biomarkers", br(),
                 h4("Observable biomarkers"),
                 plotOutput("p_bio", height = "420px"),
                 helpText("Urinary 5-HIAA is the acute arm made visible; plasma
                           vasopressin and the gastric slow-wave coupling index
                           track the nausea arm.")),
        tabPanel("8 · Dexamethasone DDI  ★", br(),
                 div(class = "thesis", HTML("The 20 mg &rarr; 12 mg dexamethasone
                   dose correction is not a safety footnote bolted on: it is a
                   COMPUTED consequence of CYP3A4 inhibition by aprepitant or
                   netupitant. Rolapitant inhibits CYP2D6 and not CYP3A4, so it
                   needs no correction — the class splits.")),
                 plotOutput("p_ddi", height = "320px"),
                 h4("Dexamethasone exposure ratio vs no NK1 antagonist"),
                 tableOutput("t_ddi")),
        tabPanel("9 · log-additivity  ★", br(),
                 div(class = "thesis", HTML("The two bars would be equal if the
                   saturating hazard preserved additivity on its way to the
                   endpoint. It does not - in this model the combination recovers
                   only about two thirds of the sum of the single-agent effects.
                   Additivity is exact in the DRIVE and nowhere after it.")),
                 plotOutput("p_logadd", height = "340px"),
                 tableOutput("t_logadd")),
        tabPanel("10 · Who benefits from the 4th agent  ★", br(),
                 div(class = "thesis", HTML("Corollary of the same algebra:
                   CR is MULTIPLIED by exp(&Delta;I), so the ABSOLUTE gain from
                   adding olanzapine is a non-monotonic function of baseline risk.
                   It is small in patients already controlled and small in patients
                   whose drive is overwhelming.")),
                 plotOutput("p_fourth", height = "380px"),
                 tableOutput("t_fourth")),
        tabPanel("11 · Scenario comparison", br(),
                 selectInput("scens", "Scenarios", choices = names(SCEN),
                             multiple = TRUE,
                             selected = c("S03_ond_dex", "S04_ond_dex_apr",
                                          "S06_pal_dex_apr", "S12_quadruplet")),
                 actionButton("go_scen", "Run selected scenarios"),
                 br(), br(), plotOutput("p_scen", height = "340px"),
                 tableOutput("t_scen")),
        tabPanel("12 · Anticipatory carryover", br(),
                 h4("Cycle-to-cycle conditioning (set 'Cycles' above to 6)"),
                 plotOutput("p_antic", height = "380px"),
                 helpText("Extinction is far slower than acquisition, so a badly
                           controlled first cycle raises the conditioned term for
                           the rest of the course. No receptor antagonist reaches
                           this term; only GABA-A and behavioural therapy do.")),
        tabPanel("13 · Safety", br(),
                 h4("Off-target and safety readouts"),
                 plotOutput("p_safety", height = "420px"),
                 helpText("The QTc panel reproduces the dose dependence that led
                           to the 2012 withdrawal of the 32 mg intravenous
                           ondansetron dose.")),
        tabPanel("14 · Systemic consequences", br(),
                 h4("Fluid, potassium, creatinine, dose intensity, quality of life"),
                 plotOutput("p_sys", height = "420px"),
                 helpText("The terminal endpoint is relative dose intensity: the
                           mechanism by which uncontrolled CINV reaches survival.")),
        tabPanel("15 · Model / provenance", br(),
                 h4("Fitted vs held-out"), tableOutput("t_prov"),
                 h4("Tonic self-test (must be ~0 — not a calibrated constant)"), verbatimTextOutput("v_tonic"),
                 h4("Full parameter set"), tableOutput("t_par"))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  covars <- reactive(list(
    FEMALE = as.numeric(input$sex), AGE = input$age,
    ALCOHOL = as.numeric(input$alcohol), MOTION = as.numeric(input$motion),
    PRIORCINV = as.numeric(input$prior), PH2D6 = as.numeric(input$cyp2d6)))

  sim <- reactive({
    end <- if (input$cycles > 1) input$cycles * 504 else 120
    simulate_custom(input$drugs, input$chemo, covars(), end = end,
                    cycles = input$cycles)
  })

  ep <- reactive(endpoints(sim(), omega = input$omega))

  output$kpis <- renderUI({
    e <- ep()
    fmt <- function(lab, v, s = "%.1f%%") div(class = "kpi", lab, br(),
                                              tags$b(sprintf(s, v)))
    tagList(
      fmt("CR acute (0-24 h)", 100 * e$CR_acute),
      fmt("CR delayed (24-120 h)", 100 * e$CR_delayed),
      fmt("CR overall (0-120 h)", 100 * e$CR_overall),
      fmt("No nausea (0-120 h)", 100 * e$NN_overall),
      fmt("Mean episodes 0-24 h", e$emesis_24, "%.2f"),
      fmt("Peak nausea (0-10)", e$peak_nausea, "%.1f"),
      fmt("Peak dQTcF (ms)", e$peak_QTc, "%.1f"))
  })

  output$p_endpoints <- renderPlot({
    e <- ep()
    d <- data.frame(
      window = factor(rep(c("acute 0-24 h", "delayed 24-120 h", "overall 0-120 h"), 2),
                      levels = c("acute 0-24 h", "delayed 24-120 h", "overall 0-120 h")),
      endpoint = rep(c("complete response", "no nausea"), each = 3),
      value = 100 * c(e$CR_acute, e$CR_delayed, e$CR_overall,
                      e$NN_acute, e$NN_delayed, e$NN_overall))
    ggplot(d, aes(window, value, fill = endpoint)) +
      geom_col(position = "dodge", width = 0.7) +
      geom_text(aes(label = sprintf("%.1f", value)),
                position = position_dodge(0.7), vjust = -0.4, size = 3.6) +
      scale_fill_manual(values = PAL[1:2]) +
      labs(x = NULL, y = "% of patients", fill = NULL) + ylim(0, 105)
  })

  output$t_endpoints <- renderTable({
    e <- ep()
    data.frame(Endpoint = c("CR acute", "CR delayed", "CR overall",
                            "No nausea overall", "Mean episodes 0-24 h",
                            "Mean episodes 0-120 h", "Rescue doses 0-120 h",
                            "Peak nausea", "FLIE", "Relative dose intensity",
                            "Hazard integral acute", "Hazard integral delayed"),
               Value = sprintf("%.4g", c(e$CR_acute, e$CR_delayed, e$CR_overall,
                 e$NN_overall, e$emesis_24, e$emesis_120, e$rescue_120,
                 e$peak_nausea, e$FLIE, e$RDI, e$I_acute, e$I_delayed)))
  }, rownames = FALSE)

  output$risknote <- renderText({
    d <- sim()
    s <- d$SUSC_o[1]
    sprintf(paste("Patient susceptibility multiplier = %.2f. This multiplies the",
                  "HAZARD, so it changes CR multiplicatively, not additively.",
                  "CR(0-24) x CR(24-120) = %.3f but the joint CR(0-120) = %.3f —",
                  "the difference is entirely the population variance omega."),
            s, ep()$CR_acute * ep()$CR_delayed, ep()$CR_overall)
  })

  output$p_pk <- renderPlot({
    d <- sim()
    keep <- c(cOND_o = "ondansetron", cPAL_o = "palonosetron",
              cAPR_o = "aprepitant", cDEX_o = "dexamethasone",
              cOLA_o = "olanzapine")
    dl <- d %>% select(time, all_of(names(keep))) %>%
      pivot_longer(-time, names_to = "k", values_to = "conc") %>%
      mutate(drug = keep[k]) %>% filter(conc > 1e-4)
    if (!nrow(dl)) return(ggplot() + labs(title = "no agent selected"))
    ggplot(dl, aes(time, conc, colour = drug)) + geom_line(linewidth = 0.85) +
      scale_y_log10() + scale_colour_manual(values = PAL) +
      labs(x = "time (h)", y = "plasma concentration (nM, log scale)", colour = NULL)
  })

  output$p_occ <- renderPlot({
    d <- sim()
    a <- d %>% select(time, R3S, R3O, R3G, R3P, R3I) %>%
      pivot_longer(-time, names_to = "state", values_to = "frac") %>%
      mutate(pool = "5-HT3 receptor pool")
    b <- d %>% select(time, N1S, N1A, N1N, N1R) %>%
      pivot_longer(-time, names_to = "state", values_to = "frac") %>%
      mutate(pool = "NK1 receptor pool")
    ggplot(rbind(a, b), aes(time, frac, fill = state)) +
      geom_area(position = "stack") + facet_wrap(~pool, ncol = 1) +
      scale_fill_brewer(palette = "Set2") + ylim(0, 1) +
      labs(x = "time (h)", y = "fraction of receptor pool", fill = NULL)
  })

  output$p_terms <- renderPlot({
    d <- sim() %>% select(time, TERM_VAG, TERM_NK1, TERM_D2, TERM_AP) %>%
      pivot_longer(-time, names_to = "term", values_to = "v") %>%
      mutate(term = recode(term,
        TERM_VAG = "vagal / 5-HT3  (WVAG.VAG)",
        TERM_NK1 = "central NK1  (WNKC.N1S)",
        TERM_D2  = "dopamine D2  (WD2.occD2)",
        TERM_AP  = "area postrema, RECEPTOR-INDEPENDENT  (WAP.FG)"))
    ggplot(d, aes(time, v, fill = term)) + geom_area(position = "stack") +
      scale_fill_manual(values = c(PAL[1], PAL[4], PAL[3], PAL[2])) +
      labs(x = "time (h)", y = "additive contribution to the emetic sum",
           fill = NULL) + guides(fill = guide_legend(nrow = 2))
  })

  output$p_drive <- renderPlot({
    d <- sim()
    dl <- data.frame(time = rep(d$time, 2),
                     v = c(d$DRIVE, d$EBS_o),
                     k = rep(c("DRIVE (total sum)", "EXCESS over tonic (drives the hazard)"),
                             each = nrow(d)))
    ggplot(dl, aes(time, v, colour = k)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = as.numeric(param(mod)$THR), linetype = "dashed",
                 colour = "grey40") +
      annotate("text", x = max(d$time) * 0.82, y = as.numeric(param(mod)$THR),
               label = "THR (hazard half-max)", vjust = -0.5, size = 3.3,
               colour = "grey30") +
      scale_colour_manual(values = c(PAL[1], PAL[2])) +
      labs(x = "time (h)", y = "drive (-)", colour = NULL)
  })

  output$p_phase <- renderPlot({
    d <- sim()
    dl <- data.frame(
      time = rep(d$time, 4),
      v = c(d$SHT_G / max(1e-9, max(d$SHT_G)), d$SP_C / max(1e-9, max(d$SP_C)),
            d$FG_o / max(1e-9, max(d$FG_o)), d$Z2 / max(1e-9, max(d$Z2))),
      k = rep(c("gut 5-HT (ACUTE arm)", "central substance P (DELAYED arm)",
                "mucosal injury signal FG", "delayed transit signal Z2"),
              each = nrow(d)))
    ggplot(dl, aes(time, v, colour = k)) + geom_line(linewidth = 0.9) +
      geom_vline(xintercept = 24, linetype = "dotted") +
      annotate("text", x = 24, y = 1.03, label = "24 h", size = 3.2, hjust = -0.15) +
      scale_colour_manual(values = PAL[c(2, 4, 5, 6)]) +
      labs(x = "time (h)", y = "normalised to own maximum", colour = NULL) +
      guides(colour = guide_legend(nrow = 2))
  })

  output$p_nausea <- renderPlot({
    d <- sim()
    dl <- data.frame(time = rep(d$time, 3),
      v = c(d$NAUSEA, 10 * d$lam_o / max(1e-9, max(d$lam_o)), d$EMES),
      k = rep(c("NAUSEA (graded, cortical sum)",
                "emetic hazard (threshold, brainstem sum; scaled)",
                "cumulative expected episodes"), each = nrow(d)))
    ggplot(dl, aes(time, v, colour = k)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = PAL[c(4, 2, 1)]) +
      labs(x = "time (h)", y = "value", colour = NULL) +
      guides(colour = guide_legend(nrow = 3))
  })

  output$p_ctxbar <- renderPlot({
    d <- sim()
    i <- which.max(d$CTX)
    p <- param(mod)
    terms <- c(
      "vagal" = as.numeric(p$NVAG) * d$VAG[i],
      "central NK1" = as.numeric(p$NNK) * d$N1S[i],
      "vasopressin" = as.numeric(p$NAVP) * d$AVP[i] / as.numeric(p$KAVP),
      "gastric dysrhythmia" = as.numeric(p$NSW) * max(0, 1 - d$SWC[i]),
      "anticipatory + anxiety" = as.numeric(p$NANT) * (d$ANTIC[i] + d$ANX[i]),
      "receptor-independent" = as.numeric(p$NAP) * d$FG_o[i])
    dd <- data.frame(term = names(terms), v = as.numeric(terms))
    ggplot(dd, aes(reorder(term, v), v)) +
      geom_col(fill = PAL[4], width = 0.65) + coord_flip() +
      labs(x = NULL, y = "contribution to the cortical sum at its peak")
  })

  output$p_bio <- renderPlot({
    d <- sim()
    dl <- data.frame(time = rep(d$time, 5),
      v = c(d$UHIAA, d$AVP, d$SWC, d$SHT_P, d$GVOL / 100),
      k = rep(c("cumulative urinary 5-HIAA (umol)", "plasma vasopressin (pg/mL)",
                "slow-wave coupling index (-)", "plasma 5-HT (nM)",
                "gastric volume (100 mL)"), each = nrow(d)))
    ggplot(dl, aes(time, v)) + geom_line(colour = PAL[1], linewidth = 0.85) +
      facet_wrap(~k, scales = "free_y", ncol = 2) +
      labs(x = "time (h)", y = NULL)
  })

  ddi <- reactive({
    base <- list(); res <- list()
    arms <- list("dexamethasone alone" = c("dexF"),
                 "+ aprepitant" = c("dexF", "apr"),
                 "+ netupitant" = c("dexF", "net"),
                 "+ rolapitant" = c("dexF", "rol"))
    for (nm in names(arms)) {
      d <- simulate_custom(arms[[nm]], input$chemo, covars(), end = 120)
      auc <- sum(diff(d$time) * head(d$cDEX_o, -1))
      res[[nm]] <- data.frame(arm = nm, AUC = auc,
                              Cmax = max(d$cDEX_o), CYP3A4 = min(d$E3A4))
    }
    r <- do.call(rbind, res)
    r$AUC_ratio <- r$AUC / r$AUC[1]
    r
  })

  output$p_ddi <- renderPlot({
    r <- ddi()
    ggplot(r, aes(reorder(arm, AUC_ratio), AUC_ratio)) +
      geom_col(fill = PAL[2], width = 0.6) + coord_flip() +
      geom_text(aes(label = sprintf("x%.2f", AUC_ratio)), hjust = -0.15, size = 4) +
      ylim(0, max(r$AUC_ratio) * 1.18) +
      labs(x = NULL, y = "dexamethasone AUC ratio vs no NK1 antagonist")
  })
  output$t_ddi <- renderTable({
    r <- ddi()
    data.frame(Arm = r$arm, `AUC (nM h)` = sprintf("%.0f", r$AUC),
               `AUC ratio` = sprintf("%.2f", r$AUC_ratio),
               `Cmax (nM)` = sprintf("%.0f", r$Cmax), check.names = FALSE)
  }, rownames = FALSE)

  logadd <- reactive({
    arms <- list("no prophylaxis" = character(0),
                 "dexamethasone only" = "dexF",
                 "ondansetron only" = "ond",
                 "aprepitant only" = "apr",
                 "olanzapine only" = "ola10",
                 "ond + dex" = c("ond", "dexF"),
                 "ond + dex + apr" = c("ond", "dexC", "apr"),
                 "quadruplet" = c("ond", "dexC", "apr", "ola10"))
    do.call(rbind, lapply(names(arms), function(nm) {
      d <- simulate_custom(arms[[nm]], input$chemo, covars(), end = 120)
      I <- at_time(d, "HINT", 120)
      data.frame(arm = nm, I = I, CR = frailty_CR(I, input$omega))
    }))
  })

  output$p_logadd <- renderPlot({
    r <- logadd()
    I0 <- r$I[r$arm == "no prophylaxis"]
    singles <- r %>% filter(arm %in% c("dexamethasone only", "ondansetron only",
                                       "aprepitant only", "olanzapine only")) %>%
      mutate(dI = I0 - I)
    pred <- I0 - sum(singles$dI)
    obs <- r$I[r$arm == "quadruplet"]
    dd <- data.frame(
      quantity = c("sum of single-agent dI", "dI of the quadruplet"),
      v = c(sum(singles$dI), I0 - obs))
    ggplot(dd, aes(quantity, v)) +
      geom_col(fill = PAL[4], width = 0.5) +
      geom_text(aes(label = sprintf("%.2f", v)), vjust = -0.4, size = 4.4) +
      labs(x = NULL,
           y = "reduction in the hazard integral (dI)",
           subtitle = paste0("Additivity holds in dI (hence in log CR). ",
             "Residual = ", sprintf("%.1f%%", 100 * (1 - (I0 - obs) / sum(singles$dI))),
             " and comes only from the saturating threshold, not from receptor interaction."))
  })

  output$t_logadd <- renderTable({
    r <- logadd()
    I0 <- r$I[r$arm == "no prophylaxis"]
    data.frame(Arm = r$arm,
               `hazard integral I` = sprintf("%.3f", r$I),
               `dI vs none` = sprintf("%.3f", I0 - r$I),
               `CR 0-120 h` = sprintf("%.3f", r$CR),
               `log CR` = sprintf("%.3f", log(pmax(r$CR, 1e-12))),
               check.names = FALSE)
  }, rownames = FALSE)

  output$p_fourth <- renderPlot({
    trip <- simulate_custom(c("pal", "dexC", "apr"), input$chemo, covars(), end = 120)
    quad <- simulate_custom(c("pal", "dexC", "apr", "ola10"), input$chemo,
                            covars(), end = 120)
    I3 <- at_time(trip, "HINT", 120)
    I4 <- at_time(quad, "HINT", 120)
    k <- exp(seq(-2.5, 2.5, length.out = 61))
    dd <- data.frame(
      CR_triplet = sapply(k, function(x) frailty_CR(x * I3, input$omega)),
      abs_gain = sapply(k, function(x) frailty_CR(x * I4, input$omega) -
                                       frailty_CR(x * I3, input$omega)))
    pk <- dd[which.max(dd$abs_gain), ]
    ggplot(dd, aes(100 * CR_triplet, 100 * abs_gain)) +
      geom_line(linewidth = 1, colour = PAL[2]) +
      geom_point(data = pk, aes(100 * CR_triplet, 100 * abs_gain), size = 3) +
      labs(x = "baseline CR on the triplet (%)",
           y = "ABSOLUTE gain in CR from adding olanzapine (percentage points)",
           subtitle = sprintf(paste("Maximal absolute benefit %.1f points, at a",
             "baseline CR of %.0f%%. The benefit vanishes at BOTH ends — which is",
             "why a fourth agent is worth adding in HEC and not in LEC."),
             100 * pk$abs_gain, 100 * pk$CR_triplet))
  })

  output$t_fourth <- renderTable({
    trip <- simulate_custom(c("pal", "dexC", "apr"), input$chemo, covars(), end = 120)
    quad <- simulate_custom(c("pal", "dexC", "apr", "ola10"), input$chemo,
                            covars(), end = 120)
    I3 <- at_time(trip, "HINT", 120)
    I4 <- at_time(quad, "HINT", 120)
    k <- c(0.15, 0.35, 0.7, 1, 1.6, 3, 6)
    data.frame(`susceptibility x` = k,
      `CR triplet` = sprintf("%.3f", sapply(k, function(x) frailty_CR(x * I3, input$omega))),
      `CR quadruplet` = sprintf("%.3f", sapply(k, function(x) frailty_CR(x * I4, input$omega))),
      `absolute gain` = sprintf("%+.3f", sapply(k, function(x)
        frailty_CR(x * I4, input$omega) - frailty_CR(x * I3, input$omega))),
      `CR ratio` = sprintf("%.2f", sapply(k, function(x)
        frailty_CR(x * I4, input$omega) / frailty_CR(x * I3, input$omega))),
      check.names = FALSE)
  }, rownames = FALSE)

  scen_res <- eventReactive(input$go_scen, {
    do.call(rbind, lapply(input$scens, function(n) {
      d <- run_scenario(n)
      cbind(scenario = n, endpoints(d, omega = input$omega))
    }))
  })

  output$p_scen <- renderPlot({
    r <- scen_res()
    dl <- r %>% select(scenario, CR_acute, CR_delayed, CR_overall, NN_overall) %>%
      pivot_longer(-scenario, names_to = "endpoint", values_to = "v")
    ggplot(dl, aes(scenario, 100 * v, fill = endpoint)) +
      geom_col(position = "dodge", width = 0.75) +
      scale_fill_manual(values = PAL[1:4]) + coord_flip() +
      labs(x = NULL, y = "% of patients", fill = NULL)
  })
  output$t_scen <- renderTable({
    r <- scen_res()
    r %>% mutate(across(where(is.numeric), ~sprintf("%.3g", .x)))
  }, rownames = FALSE)

  output$p_antic <- renderPlot({
    d <- sim()
    dl <- data.frame(time = rep(d$time, 3),
      v = c(d$ANTIC, d$ANX, d$NAUSEA / 10),
      k = rep(c("conditioned strength ANTIC", "anxiety ANX", "nausea / 10"),
              each = nrow(d)))
    ggplot(dl, aes(time / 24, v, colour = k)) + geom_line(linewidth = 0.85) +
      scale_colour_manual(values = PAL[c(2, 5, 4)]) +
      labs(x = "time (days)", y = "value", colour = NULL)
  })

  output$p_safety <- renderPlot({
    d <- sim()
    dl <- data.frame(time = rep(d$time, 5),
      v = c(d$QTC, d$GLU, d$SED, d$AKATH, d$CONST),
      k = rep(c("dQTcF (ms)", "plasma glucose (mmol/L)", "sedation (0-10)",
                "akathisia / EPS score", "constipation index"), each = nrow(d)))
    ggplot(dl, aes(time, v)) + geom_line(colour = PAL[2], linewidth = 0.85) +
      facet_wrap(~k, scales = "free_y", ncol = 2) + labs(x = "time (h)", y = NULL)
  })

  output$p_sys <- renderPlot({
    d <- sim()
    dl <- data.frame(time = rep(d$time, 5),
      v = c(d$ECFV, d$KP, d$SCR, 100 * d$RDI, d$FLIE),
      k = rep(c("ECF volume deficit (L)", "plasma potassium (mmol/L)",
                "serum creatinine (mg/dL)", "relative dose intensity (%)",
                "FLIE score"), each = nrow(d)))
    ggplot(dl, aes(time, v)) + geom_line(colour = PAL[3], linewidth = 0.85) +
      facet_wrap(~k, scales = "free_y", ncol = 2) + labs(x = "time (h)", y = NULL)
  })

  output$t_prov <- renderTable({
    data.frame(
      Parameter = c("KRE", "EPPT", "WNKC", "WAP", "THR", "EDEXBS", "OMEGA",
                    "NK50", "NN50", "tonic reference values", "everything else"),
      Status = c(rep("FITTED", 9), "ANALYTIC functions of the parameters (self-test = 0)",
                 "literature prior (PK, Ki, koff, physiology)"),
      `Anchored to` = c(
        rep("Hesketh 2003 JCO both arms x 3 windows + untreated natural history", 7),
        "Navari 2016 NEJM CONTROL arm no-nausea + untreated natural history",
        "Navari 2016 NEJM CONTROL arm no-nausea + untreated natural history",
        "-", "-"),
      check.names = FALSE)
  }, rownames = FALSE)

  output$v_tonic <- renderPrint({
    em <- EMETO[[input$chemo]]
    print(tonic_selftest(mod, c(list(EMETO_P = em[["P"]], EMETO_C = em[["C"]],
                                     CISRATE = em[["rate"]]), covars())))
  })

  output$t_par <- renderTable({
    p <- as.list(param(mod))
    data.frame(Parameter = names(p), Value = sprintf("%.6g", unlist(p)))
  }, rownames = FALSE)
}

shinyApp(ui, server)
