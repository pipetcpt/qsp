## ===========================================================================
##  mmd_shiny_app.R — Moyamoya Disease QSP explorer
## ===========================================================================
##
##  The app is organised around the model's thesis rather than around its
##  compartments: every tab exists to make ONE claim inspectable.
##
##    1. Patient          who this hemisphere is, and where it sits on gS*
##    2. Reserve cliff    the three things that switch together at gS*
##    3. Haemodynamics    the network, the four inlets, pressures and flows
##    4. Oxygen           OEF ceiling, and the DERIVED infarct threshold
##    5. Probes           acetazolamide vs PaCO2 — not the same probe
##    6. Blood pressure   the paradox, and its computable optimum
##    7. Revascularisation direct / indirect / combined, and the graft
##    8. Hyperperfusion   focal, relative, and priced in pre-op ischaemia
##    9. Endpoints        infarct, haemorrhage, cognition, event hazards
##   10. Scenarios        side-by-side comparison of any two arms
##   11. Population       what actually predicts the outcome
##   12. Model            equations, parameters, provenance, validation
##
##  Run:  shiny::runApp("mmd_shiny_app.R")
##  Requires mmd_mrgsolve_model.R in the same directory.
##
##  EDUCATIONAL / RESEARCH USE ONLY — not for clinical decision-making.
## ===========================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

## --- load the model (compiles once) ---------------------------------------
local({
  src <- readLines("mmd_mrgsolve_model.R")
  k <- grep("^if \\(identical\\(environment", src)
  if (length(k) == 0) k <- length(src) + 1
  eval(parse(text = paste(src[seq_len(k - 1)], collapse = "\n")),
       envir = globalenv())
})

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey35"),
        legend.position = "bottom", legend.title = element_blank())

CBF_CRIT_LAB <- "CBF_crit = CMRO2/(CaO2 x OEF_max)"

## ===========================================================================
##  UI
## ===========================================================================
ui <- fluidPage(
  title = "Moyamoya Disease QSP",
  tags$head(tags$style(HTML("
    .hdr {background:#1f3864;color:#fff;padding:10px 14px;margin-bottom:8px;
          border-radius:4px}
    .hdr h3 {margin:0;font-weight:700}
    .hdr p  {margin:3px 0 0 0;font-size:12px;opacity:.9}
    .claim  {background:#fff8e1;border-left:4px solid #f5b942;padding:8px 12px;
             margin-bottom:10px;font-size:13px}
    .warn   {background:#fdecea;border-left:4px solid #c0392b;padding:8px 12px;
             margin-bottom:10px;font-size:13px}
    .kpi    {background:#f4f7fb;border:1px solid #d6e0ef;border-radius:4px;
             padding:8px;text-align:center;margin-bottom:6px}
    .kpi .v {font-size:20px;font-weight:700;color:#1f3864}
    .kpi .l {font-size:10px;color:#555;text-transform:uppercase}
  "))),

  div(class = "hdr",
      h3("Moyamoya Disease — Quantitative Systems Pharmacology"),
      p("Reserve, not flow, is the state variable. 36 ODEs + a three-node ",
        "haemodynamic network whose autoregulator has a floor.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("arch", "Archetype",
                  choices = c("Adult, ischaemic onset" = "adult_isch",
                              "Adult, haemorrhagic onset" = "adult_haem",
                              "Paediatric" = "paediatric",
                              "Sickle-cell moyamoya" = "sickle_cell",
                              "Slow / asymptomatic" = "asymptomatic"),
                  selected = "adult_isch"),
      sliderInput("days", "Follow-up (years)", 1, 15, 10, step = 1),
      hr(),
      h5("Patient biology"),
      sliderInput("ANGIO", "Angiogenic capacity ANGIO", 0.10, 1.10, 0.58,
                  step = 0.02),
      sliderInput("K_SMC", "Lesion growth K_SMC (x1e-4 /d)", 0.4, 9.0, 2.6,
                  step = 0.1),
      selectInput("RNF", "RNF213 R4810K",
                  c("wild type (quasi-moyamoya)" = "0",
                    "heterozygote" = "1", "homozygote" = "1.6"),
                  selected = "1"),
      sliderInput("PCA_INV", "PCA involvement (0-1)", 0, 1, 0, step = 0.25),
      sliderInput("HB", "Haemoglobin (g/dL)", 6, 17, 15, step = 0.5),
      sliderInput("MAP0", "Baseline MAP (mmHg)", 65, 125, 90, step = 1),
      sliderInput("G_LEAK", "Pial coupling g_leak", 0.15, 2.5, 0.80,
                  step = 0.05),
      hr(),
      h5("Intervention"),
      selectInput("surg", "Revascularisation",
                  c("none", "direct", "indirect", "combined")),
      sliderInput("surg_y", "Operate at year", 0.5, 12, 5, step = 0.5),
      checkboxGroupInput("drugs", "Medical therapy",
                         choices = c("Aspirin 100 mg/d" = "ASA",
                                     "Cilostazol 100 mg b.d." = "CILO",
                                     "Atorvastatin 20 mg/d" = "STAT",
                                     "Nifedipine GITS 60 mg/d" = "NIF",
                                     "Minocycline 100 mg b.d." = "MIN",
                                     "Antihypertensive" = "AHT")),
      sliderInput("paco2", "PaCO2 (mmHg)", 20, 60, 40, step = 1),
      hr(),
      actionButton("run", "Run", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1 Patient",
                 br(), fluidRow(
                   column(3, uiOutput("kpi_gs")),
                   column(3, uiOutput("kpi_cvr")),
                   column(3, uiOutput("kpi_cbf")),
                   column(3, uiOutput("kpi_crit"))),
                 div(class = "claim", textOutput("patient_claim")),
                 plotOutput("p_patient", height = "330px"),
                 h5("Where this hemisphere sits"), tableOutput("t_patient")),

        tabPanel("2 Reserve cliff",
                 br(), div(class = "claim",
                           HTML("A <b>critical inlet conductance gS*</b> exists.
                            Above it CBF is defended and dCBF/dMAP = 0; below it
                            the arteriole is on its floor, CBF becomes
                            pressure-passive, and the measured acetazolamide
                            response turns <b>negative</b>. All three switch at
                            the same point — and angiographic stenosis alone
                            cannot locate it, because gS = g_ICA + g_moya +
                            g_PVA + g_bypass.")),
                 plotOutput("p_cliff", height = "420px"),
                 tableOutput("t_cliff")),

        tabPanel("3 Haemodynamics",
                 br(), plotOutput("p_press", height = "280px"),
                 plotOutput("p_inlets", height = "280px"),
                 plotOutput("p_flows", height = "260px")),

        tabPanel("4 Oxygen",
                 br(), div(class = "claim",
                           HTML("The penumbral threshold is <b>not a fitted
                            parameter</b>. It is CMRO2/(CaO2 x OEF_max) —
                            19.7 mL/100g/min at Hb 15, and 37.0 at Hb 8. In
                            sickle-cell moyamoya the same brain infarcts at
                            flows a normal brain tolerates, and transfusion,
                            which touches no vessel, is what moves it.")),
                 plotOutput("p_ox", height = "300px"),
                 plotOutput("p_hb", height = "280px")),

        tabPanel("5 Probes",
                 br(), div(class = "claim",
                           HTML("Acetazolamide acts on the <b>arteriole</b>.
                            PaCO2 acts on the arteriole <b>and on the pial
                            collateral conduits</b>. In a normal brain that
                            distinction is invisible; in moyamoya the conduits
                            are the circulation, so the two probes disagree —
                            and an acetazolamide study and a breath-hold study
                            are not interchangeable.")),
                 plotOutput("p_probe", height = "340px"),
                 tableOutput("t_probe"),
                 div(class = "warn",
                     HTML("<b>The crying child.</b> Hyperventilation lowers
                       demand (protective) and constricts the conduits
                       (harmful). Where the arteriole has range left, the fall
                       in flow is large and safe. Where it does not, a smaller
                       fall is an infarct."))),

        tabPanel("6 Blood pressure",
                 br(), div(class = "claim",
                           HTML("Ischaemic hazard runs through the
                            pressure-passive term; haemorrhagic hazard runs
                            through periventricular wall stress. They have
                            <b>opposite sign in MAP</b>, so an optimum exists —
                            and it is a different number in the two
                            phenotypes.")),
                 plotOutput("p_map", height = "380px"),
                 tableOutput("t_map")),

        tabPanel("7 Revascularisation",
                 br(), div(class = "claim",
                           HTML("A graft is a <b>parallel conductance</b> — the
                            only intervention in this model that changes gS. It
                            treats both phenotypes by one mechanism: raising
                            P_A lifts flow AND collapses the gradient across
                            the fragile periventricular route, which then
                            regresses. The operation retires the vessel that
                            was going to bleed.")),
                 plotOutput("p_surg", height = "300px"),
                 plotOutput("p_graft", height = "280px")),

        tabPanel("8 Hyperperfusion",
                 br(), div(class = "warn",
                           HTML("A lumped territory <b>cannot</b> hyperperfuse.
                            The syndrome is necessarily <b>focal</b> (the cortex
                            the graft is sewn to, at P_F &gt;&gt; P_A) and
                            <b>relative</b> (against the flow that cortex had
                            adapted to). Vasoparalysis sets its height; barrier
                            re-adaptation sets its duration.")),
                 plotOutput("p_hyper", height = "340px"),
                 h5("Pial coupling trades territorial gain against focal surge"),
                 plotOutput("p_leak", height = "300px")),

        tabPanel("9 Endpoints",
                 br(), plotOutput("p_end", height = "330px"),
                 plotOutput("p_haz", height = "300px"),
                 tableOutput("t_end")),

        tabPanel("10 Scenarios",
                 br(), fluidRow(
                   column(4, selectInput("cmpA", "Arm A",
                                         names(SCENARIOS), names(SCENARIOS)[2])),
                   column(4, selectInput("cmpB", "Arm B",
                                         names(SCENARIOS), names(SCENARIOS)[11])),
                   column(4, selectInput("cmpv", "Variable",
                                         c("CBFWS", "CBFA", "CBFF", "CVRINT",
                                           "INFPCT", "PISCH", "PHEM", "gS",
                                           "gPVA", "SIGPVA", "OEFA", "HYPERREL"),
                                         "CBFWS"))),
                 verbatimTextOutput("cmp_what"),
                 plotOutput("p_cmp", height = "380px")),

        tabPanel("11 Population",
                 br(), fluidRow(
                   column(4, sliderInput("npop", "Virtual patients", 40, 400,
                                         150, step = 20)),
                   column(4, sliderInput("popy", "Horizon (years)", 2, 12, 5,
                                         step = 1)),
                   column(4, br(), actionButton("runpop", "Simulate cohort",
                                                class = "btn-primary"))),
                 plotOutput("p_pop", height = "340px"),
                 h5("Correlates of the 5-year outcome"),
                 tableOutput("t_pop")),

        tabPanel("12 Model",
                 br(), h4("The central identity"),
                 verbatimTextOutput("eqn"),
                 h4("Parameters"), tableOutput("t_par"),
                 h4("Provenance and validation"), verbatimTextOutput("prov"))
      )
    )
  )
)

## ===========================================================================
##  SERVER
## ===========================================================================
server <- function(input, output, session) {

  pars <- reactive({
    list(ANGIO = input$ANGIO, K_SMC = input$K_SMC * 1e-4,
         RNF = as.numeric(input$RNF),
         SEC = if (as.numeric(input$RNF) == 0) 1.0 else 0.0,
         PCA_INV = input$PCA_INV, HB = input$HB, MAP0 = input$MAP0,
         G_LEAK = input$G_LEAK, PACO2 = input$paco2)
  })

  drugspec <- reactive({
    d <- list()
    if ("ASA"  %in% input$drugs) d$ASA  <- list(amt = 100, ii = 1,   addl = 9999)
    if ("CILO" %in% input$drugs) d$CILO <- list(amt = 100, ii = 0.5, addl = 19999)
    if ("STAT" %in% input$drugs) d$STAT <- list(amt = 20,  ii = 1,   addl = 9999)
    if ("NIF"  %in% input$drugs) d$NIF  <- list(amt = 60,  ii = 1,   addl = 9999)
    if ("MIN"  %in% input$drugs) d$MIN  <- list(amt = 100, ii = 0.5, addl = 19999)
    if ("AHT"  %in% input$drugs) d$AHT  <- list(amt = 10,  ii = 1,   addl = 9999)
    if (length(d) == 0) NULL else d
  })

  sim <- eventReactive(input$run, {
    surg <- if (input$surg == "none") NULL
            else list(kind = input$surg, day = input$surg_y * 365)
    run_mmd(arch = input$arch, days = input$days * 365, delta = 2,
            drugs = drugspec(), surgery = surg, pars = pars()) %>%
      mutate(year = time / 365)
  }, ignoreNULL = FALSE)

  kpi <- function(v, l) HTML(sprintf(
    "<div class='kpi'><div class='v'>%s</div><div class='l'>%s</div></div>", v, l))
  last <- reactive(tail(sim(), 1))

  output$kpi_gs   <- renderUI(kpi(sprintf("%.2f", last()$gS),
                                  "inlet conductance gS"))
  output$kpi_cvr  <- renderUI(kpi(sprintf("%.0f%%", last()$CVRINT),
                                  "intrinsic reserve"))
  output$kpi_cbf  <- renderUI(kpi(sprintf("%.0f", last()$CBFWS),
                                  "watershed CBF"))
  output$kpi_crit <- renderUI(kpi(sprintf("%.1f", last()$CBFCRIT),
                                  "infarct threshold"))

  output$patient_claim <- renderText({
    d <- last()
    side <- if (d$ARPOS > 0.99) "BELOW gS* — the arteriole is on its dilatory floor, so CBF is pressure-passive and the acetazolamide response is negative."
            else "ABOVE gS* — flow is still defended, and reserve is being spent silently."
    sprintf("This hemisphere is %s  Watershed margin to threshold: %+.1f mL/100g/min.",
            side, d$MARGIN)
  })

  output$p_patient <- renderPlot({
    d <- sim() %>% select(year, CBFA, CBFF, CBFWS, CBFCRIT, CVRINT) %>%
      pivot_longer(-c(year, CBFCRIT))
    ggplot(d, aes(year, value, colour = name)) +
      geom_hline(aes(yintercept = CBFCRIT), linetype = 2, colour = "grey40") +
      geom_line(linewidth = 0.8) +
      annotate("text", x = Inf, y = last()$CBFCRIT, hjust = 1.02, vjust = -0.5,
               label = CBF_CRIT_LAB, size = 3, colour = "grey30") +
      labs(title = "Flow is defended for years while reserve is spent",
           subtitle = "Stenosis is linear in time; conductance is its fourth power",
           x = "year", y = "mL/100g/min  |  reserve %") + THEME
  })

  output$t_patient <- renderTable({
    d <- last()
    data.frame(
      quantity = c("stenosis (diameter)", "g_ICA", "g_moya", "g_PVA",
                   "g_coll", "g_bypass", "gS total", "arteriole position",
                   "P_A (mmHg)", "P_F (mmHg)", "OEF", "dCBF/dMAP",
                   "measured reactivity %"),
      value = sprintf("%.3f", c(d$STEN, d$gICA, d$gMOYA, d$gPVA, d$gCOLL,
                                d$gBYP, d$gS, d$ARPOS, d$PA, d$PF, d$OEFA,
                                d$DCBFDMAP, d$CVRMEAS)))
  }, striped = TRUE)

  ## --- 2. the cliff -------------------------------------------------------
  cliff <- reactive({
    st <- seq(0.02, 0.62, by = 0.02)
    bind_rows(lapply(st, function(x) {
      q <- mod %>% param(ARCHETYPE[[input$arch]]) %>% param(pars()) %>%
        init(init_state(ARCHETYPE[[input$arch]])) %>%
        init(SMC = x / 0.945) %>% param(PROG = 0) %>%
        mrgsim(end = 0, delta = 1) %>% as_tibble()
      q$STENset <- x; q
    }))
  })

  output$p_cliff <- renderPlot({
    d <- cliff() %>%
      select(gS, `CBF_A` = CBFA, `intrinsic reserve %` = CVRINT,
             `measured reactivity %` = CVRMEAS, `dCBF/dMAP` = DCBFDMAP) %>%
      pivot_longer(-gS)
    ggplot(d, aes(gS, value)) +
      geom_hline(yintercept = 0, colour = "grey70") +
      geom_line(linewidth = 0.9, colour = "#1f3864") +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Three things switch at the same inlet conductance",
           subtitle = "gS* ~ 2.20 mL/min/mmHg for this geometry",
           x = "total inlet conductance gS (mL/min/mmHg)", y = NULL) + THEME
  })

  output$t_cliff <- renderTable({
    cliff() %>%
      transmute(gS = round(gS, 3), stenosis = round(STENset, 3),
                CBF_A = round(CBFA, 1), CBF_ws = round(CBFWS, 1),
                arteriole = round(ARPOS, 3), CVR_int = round(CVRINT, 1),
                CVR_meas = round(CVRMEAS, 2),
                dCBF_dMAP = round(DCBFDMAP, 3)) %>%
      filter(row_number() %% 3 == 1)
  }, striped = TRUE)

  ## --- 3. haemodynamics ---------------------------------------------------
  output$p_press <- renderPlot({
    sim() %>% select(year, MAPx, PA, PF, PB) %>% pivot_longer(-year) %>%
      ggplot(aes(year, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(title = "Nodal pressures", x = "year", y = "mmHg") + THEME
  })
  output$p_inlets <- renderPlot({
    sim() %>% select(year, gICA, gMOYA, gPVA, gCOLL, gBYP) %>%
      pivot_longer(-year) %>%
      ggplot(aes(year, value, fill = name)) +
      geom_area(alpha = 0.85, colour = "white", linewidth = 0.2) +
      labs(title = "The four inlets, in parallel",
           subtitle = "the antegrade ICA is replaced, not supplemented",
           x = "year", y = "conductance (mL/min/mmHg)") + THEME
  })
  output$p_flows <- renderPlot({
    sim() %>% select(year, QICA, QMOYA, QPVA, QCOLL, QBYP, QLEAK) %>%
      pivot_longer(-year) %>%
      ggplot(aes(year, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(title = "Flow through each route", x = "year", y = "mL/min") + THEME
  })

  ## --- 4. oxygen ----------------------------------------------------------
  output$p_ox <- renderPlot({
    sim() %>% select(year, OEFA, OEFWS) %>% pivot_longer(-year) %>%
      ggplot(aes(year, value, colour = name)) +
      geom_hline(yintercept = 0.85, linetype = 2, colour = "#c0392b") +
      geom_line(linewidth = 0.9) +
      annotate("text", x = 0, y = 0.86, hjust = 0, size = 3,
               colour = "#c0392b", label = "OEF ceiling — misery perfusion") +
      labs(title = "Oxygen extraction rises first, then runs out",
           x = "year", y = "OEF") + THEME
  })
  output$p_hb <- renderPlot({
    hb <- seq(6, 17, by = 0.25)
    data.frame(HB = hb, CBF_crit = 3.30 / (1.34 * hb * 0.98 / 100 * 0.85)) %>%
      ggplot(aes(HB, CBF_crit)) +
      geom_line(linewidth = 1, colour = "#1f3864") +
      geom_vline(xintercept = input$HB, linetype = 3) +
      geom_hline(yintercept = 19.71, linetype = 2, colour = "grey50") +
      labs(title = "The infarct threshold is set by the blood, not the vessel",
           subtitle = "CBF_crit = CMRO2/(CaO2 x OEF_max) — dotted line = this patient",
           x = "haemoglobin (g/dL)", y = "CBF_crit (mL/100g/min)") + THEME
  })

  ## --- 5. probes ----------------------------------------------------------
  probe_tab <- reactive({
    base <- sim(); dy <- max(base$time)
    st <- as.numeric(tail(base, 1)[names(init_state(ARCHETYPE[[input$arch]]))])
    names(st) <- names(init_state(ARCHETYPE[[input$arch]]))
    one <- function(lbl, extra) {
      q <- mod %>% param(ARCHETYPE[[input$arch]]) %>% param(pars()) %>%
        init(st) %>% param(PROG = 0) %>% param(extra) %>%
        mrgsim(end = 0.05, delta = 0.005) %>% as_tibble()
      j <- which.max(abs(q$CBFA - q$CBFA[1]))
      data.frame(probe = lbl, CBF_A = q$CBFA[j], CBF_ws = q$CBFWS[j],
                 CBF_B = q$CBFB[j])
    }
    b <- one("baseline (PaCO2 40)", list(PACO2 = 40))
    az <- mod %>% param(ARCHETYPE[[input$arch]]) %>% param(pars()) %>%
      init(st) %>% init(AZ_C = 1000) %>% param(PROG = 0, PACO2 = 40) %>%
      mrgsim(end = 0.05, delta = 0.002) %>% as_tibble()
    ja <- which.max(abs(az$CBFA - b$CBF_A))
    bind_rows(
      b,
      data.frame(probe = "acetazolamide 1 g IV", CBF_A = az$CBFA[ja],
                 CBF_ws = az$CBFWS[ja], CBF_B = az$CBFB[ja]),
      one("hypercapnia (PaCO2 50)", list(PACO2 = 50)),
      one("hypocapnia (PaCO2 30)", list(PACO2 = 30)),
      one("hyperventilation (PaCO2 25)", list(PACO2 = 25))
    ) %>% mutate(dA = (CBF_A / CBF_A[1] - 1) * 100,
                 dB = (CBF_B / CBF_B[1] - 1) * 100)
  })

  output$p_probe <- renderPlot({
    probe_tab() %>% select(probe, dA, dB) %>%
      pivot_longer(-probe, names_to = "territory") %>%
      mutate(territory = recode(territory, dA = "affected (MCA)",
                                dB = "donor (PCA)"),
             probe = factor(probe, levels = probe_tab()$probe)) %>%
      ggplot(aes(probe, value, fill = territory)) +
      geom_col(position = "dodge") + geom_hline(yintercept = 0) +
      coord_flip() +
      labs(title = "The donor always answers. The affected territory may answer backwards.",
           subtitle = "negative = STEAL: dilating the donor takes flow from the territory it supplies",
           x = NULL, y = "% change in CBF") + THEME
  })
  output$t_probe <- renderTable({
    probe_tab() %>% mutate(across(where(is.numeric), ~round(., 2)))
  }, striped = TRUE)

  ## --- 6. blood pressure --------------------------------------------------
  map_tab <- reactive({
    base <- sim()
    st <- as.numeric(tail(base, 1)[names(init_state(ARCHETYPE[[input$arch]]))])
    names(st) <- names(init_state(ARCHETYPE[[input$arch]]))
    bind_rows(lapply(seq(65, 125, by = 5), function(mp) {
      q <- mod %>% param(ARCHETYPE[[input$arch]]) %>% param(pars()) %>%
        init(st) %>% param(PROG = 0, MAP_MULT = mp / input$MAP0,
                           MAP_MULT_T = -1e9) %>%
        mrgsim(end = 0, delta = 1) %>% as_tibble()
      sd <- 4.50 * (1 + 0.70 * (q$ARPOS[1] > 0.99))
      ih <- 0.0480 * exp(-max(min((q$CBFWS[1] - q$CBFCRIT[1]) / sd, 25), -6)) +
        0.075 * 0.35 * st[["GMOYA"]]
      hh <- 0.4140 * st[["ANEU"]] * q$SIGPVA[1]^2
      data.frame(MAP = mp, CBF_A = q$CBFA[1], CBF_ws = q$CBFWS[1],
                 sig_pva = q$SIGPVA[1], ischaemic = ih, haemorrhagic = hh,
                 total = ih + hh)
    }))
  })
  output$p_map <- renderPlot({
    d <- map_tab(); opt <- d$MAP[which.min(d$total)]
    d %>% select(MAP, ischaemic, haemorrhagic, total) %>% pivot_longer(-MAP) %>%
      ggplot(aes(MAP, value, colour = name)) +
      geom_vline(xintercept = opt, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 0.9) +
      annotate("text", x = opt, y = Inf, vjust = 1.5, hjust = -0.05, size = 3.4,
               label = sprintf("optimum %d mmHg", opt)) +
      labs(title = "Two hazards, opposite sign in MAP",
           subtitle = "disease state frozen; only the pressure is swept",
           x = "MAP (mmHg)", y = "annual hazard") + THEME
  })
  output$t_map <- renderTable({
    map_tab() %>% mutate(across(where(is.numeric), ~round(., 4)))
  }, striped = TRUE)

  ## --- 7. revascularisation ----------------------------------------------
  surg_cmp <- reactive({
    ops <- c("none", "direct", "indirect", "combined")
    bind_rows(lapply(ops, function(op) {
      s <- if (op == "none") NULL else list(kind = op, day = input$surg_y * 365)
      run_mmd(arch = input$arch, days = input$days * 365, delta = 5,
              drugs = drugspec(), surgery = s, pars = pars()) %>%
        mutate(year = time / 365, arm = op)
    }))
  })
  output$p_surg <- renderPlot({
    surg_cmp() %>% select(year, arm, CBFWS, CBFCRIT) %>%
      ggplot(aes(year, CBFWS, colour = arm)) +
      geom_hline(aes(yintercept = CBFCRIT), linetype = 2, colour = "grey40") +
      geom_vline(xintercept = input$surg_y, linetype = 3) +
      geom_line(linewidth = 0.9) +
      labs(title = "Watershed perfusion by operation",
           subtitle = "the borderzone is what the graft has to reach",
           x = "year", y = "CBF watershed (mL/100g/min)") + THEME
  })
  output$p_graft <- renderPlot({
    surg_cmp() %>% filter(arm != "none") %>%
      select(year, arm, gBYP, gPVA, QPVA) %>% pivot_longer(-c(year, arm)) %>%
      ggplot(aes(year, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "The graft grows; the fragile route retires",
           subtitle = "raising P_A collapses the gradient across the periventricular anastomosis",
           x = "year", y = NULL) + THEME
  })

  ## --- 8. hyperperfusion --------------------------------------------------
  output$p_hyper <- renderPlot({
    ton <- input$surg_y * 365
    d <- run_mmd(arch = input$arch, days = ton + 120, delta = 0.05,
                 surgery = list(kind = "direct", day = ton), pars = pars()) %>%
      filter(time >= ton - 10) %>% mutate(day = time - ton)
    d %>% select(day, `CBF_F (focal)` = CBFF, `adapted set-point` = CBFAD,
                 `CBF_A (territory)` = CBFA) %>% pivot_longer(-day) %>%
      ggplot(aes(day, value, colour = name)) + geom_line(linewidth = 0.9) +
      geom_vline(xintercept = 0, linetype = 3) +
      labs(title = "Hyperperfusion is focal and relative",
           subtitle = "50 mL/100g/min is normal, and it injures a barrier adapted to 18",
           x = "days from surgery", y = "mL/100g/min") + THEME
  })
  output$p_leak <- renderPlot({
    ton <- input$surg_y * 365
    bind_rows(lapply(c(0.25, 0.35, 0.55, 0.80, 1.20, 1.80), function(gl) {
      p <- pars(); p$G_LEAK <- gl
      d <- run_mmd(arch = input$arch, days = ton + 120, delta = 0.2,
                   surgery = list(kind = "direct", day = ton), pars = p)
      m <- d$time >= ton
      pre <- approx(d$time, d$CBFWS, ton - 1)$y
      data.frame(g_leak = gl, territorial_gain = tail(d$CBFWS, 1) - pre,
                 focal_surge = max(d$HYPERREL[m]))
    })) %>% pivot_longer(-g_leak) %>%
      ggplot(aes(g_leak, value, colour = name)) +
      geom_line(linewidth = 1) + geom_point() +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Pial coupling trades territorial gain against focal surge",
           subtitle = "g_leak belongs to the patient, not the operation — which is why hyperperfusion is a minority event",
           x = "pial coupling g_leak (mL/min/mmHg)", y = NULL) + THEME
  })

  ## --- 9. endpoints -------------------------------------------------------
  output$p_end <- renderPlot({
    sim() %>% select(year, `infarct %` = INFPCT, `cognition z` = COG,
                     `oedema` = EDEMA, `microaneurysm` = ANEU) %>%
      pivot_longer(-year) %>%
      ggplot(aes(year, value, colour = name)) + geom_line(linewidth = 0.85) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Tissue endpoints", x = "year", y = NULL) + THEME
  })
  output$p_haz <- renderPlot({
    sim() %>% select(year, `P(ischaemic event)` = PISCH,
                     `P(haemorrhage)` = PHEM) %>% pivot_longer(-year) %>%
      ggplot(aes(year, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_y_continuous(labels = scales::percent) +
      labs(title = "Cumulative event probability",
           subtitle = "ischaemic hazard is exponential in the margin to CBF_crit; haemorrhagic hazard reads periventricular wall stress",
           x = "year", y = NULL) + THEME
  })
  output$t_end <- renderTable({
    d <- last()
    data.frame(endpoint = c("infarct fraction %", "cognition (z)",
                            "P(ischaemic event)", "P(haemorrhage)",
                            "cumulative haematoma (mL)",
                            "days below threshold", "REMOD"),
               value = sprintf("%.3f", c(d$INFPCT, d$COG, d$PISCH, d$PHEM,
                                         d$HEMV, d$TIAB, d$REMOD)))
  }, striped = TRUE)

  ## --- 10. scenarios ------------------------------------------------------
  output$cmp_what <- renderText({
    paste0("A — ", input$cmpA, ": ",
           gsub("\\s+", " ", SCENARIOS[[input$cmpA]]$what), "\n\n",
           "B — ", input$cmpB, ": ",
           gsub("\\s+", " ", SCENARIOS[[input$cmpB]]$what))
  })
  output$p_cmp <- renderPlot({
    d <- bind_rows(run_scenario(input$cmpA), run_scenario(input$cmpB)) %>%
      mutate(year = time / 365)
    ggplot(d, aes(year, .data[[input$cmpv]], colour = scenario)) +
      geom_line(linewidth = 0.9) +
      labs(title = paste("Scenario comparison:", input$cmpv),
           x = "year", y = input$cmpv) + THEME
  })

  ## --- 11. population ----------------------------------------------------
  pop <- eventReactive(input$runpop, {
    withProgress(message = "simulating cohort", {
      virtual_population(n = input$npop, days = input$popy * 365)
    })
  })
  output$p_pop <- renderPlot({
    p <- pop()
    ggplot(p, aes(CVRINT, INFPCT, colour = factor(PCA_INV), size = gS)) +
      geom_point(alpha = 0.7) +
      labs(title = "Reserve predicts tissue loss; the angiogram does not",
           subtitle = "colour = PCA involvement, size = inlet conductance gS",
           x = "intrinsic reserve at horizon (%)",
           y = "infarct fraction (%)") + THEME
  })
  output$t_pop <- renderTable({
    p <- pop()
    xs <- c("RNF", "ANGIO", "K_SMC", "HB", "MAP0", "PCA_INV", "G_LEAK",
            "STEN", "gS", "CBFWS", "CVRINT")
    ys <- c("INFPCT", "PISCH", "PHEM")
    out <- expand.grid(predictor = xs, outcome = ys, stringsAsFactors = FALSE)
    out$r <- round(mapply(function(a, b)
      suppressWarnings(cor(p[[a]], p[[b]])), out$predictor, out$outcome), 3)
    out %>% pivot_wider(names_from = outcome, values_from = r)
  }, striped = TRUE)

  ## --- 12. model ---------------------------------------------------------
  output$eqn <- renderText(paste(
    "THREE-NODE NETWORK (one hemisphere)",
    "",
    "  A: gSA(Pa-P_A) + gc(P_B-P_A) + gl(P_F-P_A) = gA(P_A-Pv)",
    "  F: (gSF+gb)(Pa-P_F) + gl(P_A-P_F)          = gF(P_F-Pv)",
    "  B: gp(Pa-P_B) + gc(P_A-P_B)                = gB(P_B-Pv)",
    "",
    "  gS = g_ICA(1-STEN)^4 + g_moya + g_PVA + g_bypass",
    "  autoregulation picks gA, gF, gB to meet demand, CLAMPED to",
    "     [1/R_art_max_eff, 1/R_art_min]   <-- the whole model is this clamp",
    "",
    "DERIVED, NOT FITTED",
    "  CBF_crit = CMRO2 / (CaO2 * OEF_max)",
    "           = 19.71 mL/100g/min at Hb 15 g/dL",
    "           = 36.96 mL/100g/min at Hb 8 g/dL",
    "  gS*      = 2.201 mL/min/mmHg (reserve = 0, CBF becomes pressure-passive,",
    "             measured reactivity turns negative -- all at the same point)",
    "",
    "WATERSHED",
    "  CBF_ws = CBF_A * [1 - WS_K*(1 - antegrade flow fraction)]",
    "  collateral flow reaches the cortical surface first, so the deep",
    "  borderzone is supplied last -- which is why the infarcts are there",
    "",
    "PERIVENTRICULAR WALL STRESS (the haemorrhage source)",
    "  sigma_PVA = 0.5(MAP + P_A) * dilatation^1.5 / sigma_ref,",
    "  dilatation = (g_PVA/g_ref)^(1/8)   [G ~ N r^4, both N and r grow]",
    "",
    "EVENT HAZARDS",
    "  ischaemic:    h = h0 * exp(-(CBF_ws - CBF_crit)/sd),",
    "                sd widens when the arteriole is on its floor",
    "  haemorrhagic: h = h0 * ANEU * sigma_PVA^2",
    sep = "\n"))

  output$t_par <- renderTable({
    p <- as.list(param(mod))
    data.frame(parameter = names(p), value = unlist(p, use.names = FALSE))
  }, striped = TRUE)

  output$prov <- renderText(paste(
    "CALIBRATION ANCHORS (PMIDs in mmd_references.md)",
    "  normal cortical CBF 50, CMRO2 3.3, OEF 0.33",
    "  large-artery share of cerebrovascular resistance ~25%",
    "  acetazolamide 1 g IV raises normal CBF 30-40%   -> AZ_EMAX 0.45",
    "  CO2 reactivity ~3.5 %/mmHg PaCO2                -> K_CO2 0.035",
    "  symptomatic adult MMD: affected CBF 15-25% below normal, OEF 0.42-0.50",
    "  JAM Trial CONSERVATIVE 5-y rebleeding 31.6%     -> HEM_HAZ0 0.414",
    "     (the bypass HR is then a PREDICTION, not a fit)",
    "  post-bypass hyperperfusion peaks day 2-7, resolves over 2-3 weeks",
    "     -> TAU_ADAPT 22 d and K_REM_OFF 0.05 /d",
    "",
    "VALIDATION",
    "  Every equation is implemented TWICE: once here in mrgsolve/C++ and once",
    "  in mmd_reference_model.py (numpy/scipy LSODA).  The two agree to a",
    "  median of 0.000% and a maximum of 0.15% across 80 paired values",
    "  (mmd_cross_validation.txt).  Run validate() for the anchor checks and",
    "  jam_emulation() for the trial reproduction.",
    "",
    "KNOWN LIMITATION",
    "  A lumped territory cannot hyperperfuse; the focal compartment exists",
    "  because of that, and its coupling g_leak is the least well constrained",
    "  parameter in the model.  It is exposed as a slider for that reason.",
    "",
    "EDUCATIONAL / RESEARCH USE ONLY -- not for clinical decision-making.",
    sep = "\n"))
}

shinyApp(ui, server)
