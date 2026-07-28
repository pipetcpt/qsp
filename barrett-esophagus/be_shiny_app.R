## ============================================================================
## Barrett's Esophagus (BE) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 9 tabs: Patient profile · Drug PK & acid PD · Reflux/injury axis ·
##         Metaplasia & proliferation · Clonal evolution · Clinical endpoints ·
##         Scenario comparison · Ablation & recurrence · Pharmacogenomics &
##         safety   (+ References)
## Dependencies: shiny, shinydashboard, mrgsolve, dplyr, tidyr, ggplot2, DT
## Run:  shiny::runApp("be_shiny_app.R")   (with be_mrgsolve_model.R alongside)
## ----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DT)
})

## ---------- Lazy model load -------------------------------------------------
get_model <- function() {
  if (!exists(".BE_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.BE_MOD)) {
    assign(".BE_MOD", mread_cache("be_mrgsolve_model", project = "."), envir = .GlobalEnv)
  }
  .GlobalEnv$.BE_MOD
}

## ---------- Dosing vocabulary (1.0 au = esomeprazole 20 mg / ASA 300 mg) ----
ppi_low  <- function(d) ev(amt = 1.0, cmt = "PPI_DEPOT", ii = 1,   addl = max(d - 1, 0))
ppi_high <- function(d) ev(amt = 2.0, cmt = "PPI_DEPOT", ii = 0.5, addl = max(2 * d - 1, 0))
asa_300  <- function(d) ev(amt = 1.0, cmt = "ASA_DEPOT", ii = 1,   addl = max(d - 1, 0))
vpz_20   <- function(d) ev(amt = 1.0, cmt = "VPZ_DEPOT", ii = 1,   addl = max(d - 1, 0))
udca_600 <- function(d) ev(amt = 1.0, cmt = "UDCA_CENT", ii = 0.5, addl = max(2 * d - 1, 0))
rfa_x    <- function(start, n, every = 56)
  ev(amt = 1.0, cmt = "ABL_PULSE", time = start, ii = every, addl = max(n - 1, 0))
emr_1    <- function(start) ev(amt = 1.0, cmt = "EMR_PULSE", time = start)

NULL_EV <- ev(amt = 0, cmt = "PPI_DEPOT")

## ---------- Patient archetypes ---------------------------------------------
ARCHETYPES <- list(
  "NDBE, 3 cm (non-dysplastic)" =
    list(BASE_BE_LEN = 3, INIT_FP16 = 0.45, INIT_FLGD = 0.03,
         INIT_FP53 = 0.05, INIT_FHGD = 0.00, INIT_ANEUPL = 0.05),
  "Confirmed LGD, 5 cm" =
    list(BASE_BE_LEN = 5, INIT_FP16 = 0.80, INIT_FLGD = 0.50,
         INIT_FP53 = 0.50, INIT_FHGD = 0.00, INIT_ANEUPL = 0.25),
  "HGD, 6 cm (pre-EET)" =
    list(BASE_BE_LEN = 6, INIT_FP16 = 0.90, INIT_FLGD = 0.60,
         INIT_FP53 = 0.70, INIT_FHGD = 0.45, INIT_ANEUPL = 0.45),
  "Ultra-short segment, 1 cm" =
    list(BASE_BE_LEN = 1, INIT_FP16 = 0.30, INIT_FLGD = 0.01,
         INIT_FP53 = 0.02, INIT_FHGD = 0.00, INIT_ANEUPL = 0.02)
)

CYP_PHENO <- c("UM (*17/*17)" = 1.8, "RM (*1/*17)" = 1.4, "NM (*1/*1)" = 1.0,
               "IM (*1/*2)" = 0.55, "PM (*2/*2)" = 0.25)

SCENARIOS <- c(
  "No therapy (natural history)",
  "Esomeprazole 20 mg OD (AspECT low dose)",
  "Esomeprazole 40 mg BID (AspECT high dose)",
  "Esomeprazole 20 mg OD + aspirin 300 mg",
  "Esomeprazole 40 mg BID + aspirin 300 mg",
  "Vonoprazan 20 mg OD (P-CAB)",
  "High-dose PPI + UDCA",
  "Prompt RFA x3 on high-dose PPI",
  "EMR + RFA x4 on high-dose PPI",
  "Antireflux surgery at year 1"
)

build_events <- function(scenario, horizon_d) {
  d <- horizon_d
  switch(scenario,
    "No therapy (natural history)"                  = NULL_EV,
    "Esomeprazole 20 mg OD (AspECT low dose)"       = ppi_low(d),
    "Esomeprazole 40 mg BID (AspECT high dose)"     = ppi_high(d),
    "Esomeprazole 20 mg OD + aspirin 300 mg"        = ppi_low(d) + asa_300(d),
    "Esomeprazole 40 mg BID + aspirin 300 mg"       = ppi_high(d) + asa_300(d),
    "Vonoprazan 20 mg OD (P-CAB)"                   = vpz_20(d),
    "High-dose PPI + UDCA"                          = ppi_high(d) + udca_600(d),
    "Prompt RFA x3 on high-dose PPI"                 = ppi_high(d) + rfa_x(30, 3),
    "EMR + RFA x4 on high-dose PPI"                 = ppi_high(d) + emr_1(30) + rfa_x(65, 4),
    "Antireflux surgery at year 1"                  = NULL_EV,
    NULL_EV
  )
}

## ---------- Simulation wrapper ---------------------------------------------
run_sim <- function(scenario, horizon_y, inp) {
  mod <- get_model()
  d   <- horizon_y * 365
  arch <- ARCHETYPES[[inp$archetype]]
  par <- c(arch, list(
    BASE_AET    = inp$base_aet,
    BASE_BILE   = inp$base_bile,
    OBESITY     = inp$obesity,
    SMOKE       = as.numeric(inp$smoke),
    HH_SIZE     = inp$hh,
    AGE0        = inp$age,
    CL2C19      = as.numeric(CYP_PHENO[[inp$cyp]]),
    WL_START    = if (inp$weightloss) inp$wl_year * 365 else 1e6,
    WL_MAG      = inp$wl_mag,
    BACLOFEN    = as.numeric(inp$baclofen),
    ALGINATE    = as.numeric(inp$alginate),
    FUNDO       = if (scenario == "Antireflux surgery at year 1") 1 else 0,
    FUNDO_START = 365
  ))
  if (!is.null(inp$be_len_override) && inp$be_len_override > 0)
    par$BASE_BE_LEN <- inp$be_len_override

  mod %>% param(par) %>%
    mrgsim(events = build_events(scenario, d), end = d, delta = 5) %>%
    as_tibble() %>%
    mutate(scenario = scenario, years = time / 365)
}

## ---------- Plot theme ------------------------------------------------------
th <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.text = element_text(face = "bold"))

line_plot <- function(df, yvar, ylab, title, colour = "scenario", logy = FALSE) {
  p <- ggplot(df, aes(x = years, y = .data[[yvar]], colour = .data[[colour]])) +
    geom_line(linewidth = 0.9) +
    labs(x = "Years", y = ylab, title = title) + th
  if (logy) p <- p + scale_y_log10()
  p
}

multi_plot <- function(df, vars, labels, title, colour = "scenario") {
  long <- df %>% select(years, !!colour := all_of(colour), all_of(vars)) %>%
    pivot_longer(all_of(vars), names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, levels = vars, labels = labels))
  ggplot(long, aes(x = years, y = value, colour = .data[[colour]])) +
    geom_line(linewidth = 0.85) +
    facet_wrap(~variable, scales = "free_y") +
    labs(x = "Years", y = NULL, title = title) + th
}

## ============================================================================
## UI
## ============================================================================
ui <- dashboardPage(
  dashboardHeader(title = "Barrett's Esophagus QSP", titleWidth = 320),
  dashboardSidebar(
    width = 320,
    sidebarMenu(
      menuItem("1. Patient profile",        tabName = "profile",  icon = icon("user")),
      menuItem("2. Drug PK & acid PD",      tabName = "pk",       icon = icon("pills")),
      menuItem("3. Reflux & injury axis",   tabName = "injury",   icon = icon("fire")),
      menuItem("4. Metaplasia & growth",    tabName = "meta",     icon = icon("dna")),
      menuItem("5. Clonal evolution",       tabName = "clonal",   icon = icon("sitemap")),
      menuItem("6. Clinical endpoints",     tabName = "endpoint", icon = icon("heart-pulse")),
      menuItem("7. Scenario comparison",    tabName = "compare",  icon = icon("chart-line")),
      menuItem("8. Ablation & recurrence",  tabName = "ablation", icon = icon("bolt")),
      menuItem("9. PGx & safety",           tabName = "safety",   icon = icon("shield-halved")),
      menuItem("References",                tabName = "refs",     icon = icon("book"))
    ),
    hr(),
    selectInput("archetype", "Patient archetype", names(ARCHETYPES)),
    selectInput("scenario", "Primary scenario", SCENARIOS, selected = SCENARIOS[3]),
    sliderInput("horizon", "Horizon (years)", 1, 20, 10, step = 1),
    hr(),
    sliderInput("base_aet",  "Off-therapy AET (% pH<4)", 4, 30, 12, step = 1),
    sliderInput("base_bile", "Off-therapy bile exposure (0-100)", 0, 100, 35, step = 5),
    sliderInput("be_len_override", "Prague M override (cm; 0 = archetype)", 0, 12, 0, step = 1),
    sliderInput("obesity",   "Visceral adiposity (0-1)", 0, 1, 0.5, step = 0.05),
    sliderInput("hh",        "Hiatal hernia (cm)", 0, 8, 3, step = 1),
    sliderInput("age",       "Age (years)", 30, 85, 58, step = 1),
    checkboxInput("smoke", "Current smoker", TRUE),
    selectInput("cyp", "CYP2C19 phenotype", names(CYP_PHENO), selected = "NM (*1/*1)"),
    hr(),
    checkboxInput("baclofen", "Baclofen (TLESR reduction)", FALSE),
    checkboxInput("alginate", "Alginate raft", FALSE),
    checkboxInput("weightloss", "Weight-loss intervention", FALSE),
    conditionalPanel(
      "input.weightloss == true",
      sliderInput("wl_year", "Weight loss starts (year)", 0, 10, 2, step = 1),
      sliderInput("wl_mag", "Visceral fat reduction (fraction)", 0.1, 0.8, 0.5, step = 0.05)
    ),
    hr(),
    checkboxGroupInput("compare_set", "Scenarios to compare (tab 7)",
                       choices = SCENARIOS,
                       selected = SCENARIOS[c(1, 2, 3, 5, 8)])
  ),

  dashboardBody(
    tabItems(
      ## ---- 1. Patient profile ---------------------------------------------
      tabItem(
        "profile",
        fluidRow(
          valueBoxOutput("vb_aet", 3), valueBoxOutput("vb_len", 3),
          valueBoxOutput("vb_haz", 3), valueBoxOutput("vb_cum", 3)
        ),
        fluidRow(
          box(width = 7, title = "Baseline state vs on-therapy state", status = "primary",
              solidHeader = TRUE, DTOutput("tbl_profile")),
          box(width = 5, title = "How to read this model", status = "info",
              solidHeader = TRUE,
              HTML("<p><b>The central tension.</b> Acid suppression is not reflux
              suppression. A PPI removes the H<sup>+</sup> arm of the insult and
              inactivates pepsin above pH 4, but it does not stop the volume
              reflux that delivers bile acids — and deoxycholate is
              <i>more</i> membrane-permeant and a stronger NF-κB → CDX2 inducer
              at neutral pH. Profound acid suppression also raises gastrin 2-4×,
              and gastrin is trophic to Barrett epithelium through CCK2R.</p>
              <p>So in this model high-dose PPI produces a large fall in acid
              exposure and symptoms but only a <b>partial</b> fall in the
              mutational drive that governs progression — which is why AspECT
              found HR 0.73 for high-dose PPI, and why adding aspirin (acting on
              the PGE2/apoptosis arm the PPI cannot reach) was additive
              (HR 0.59).</p>
              <p>Set <code>K_PH_BILE</code> and <code>K_GAS_KI</code> to zero and
              PPI monotherapy becomes curative in the model. That is the model's
              falsifiable content.</p>"))
        )
      ),

      ## ---- 2. Drug PK & acid PD -------------------------------------------
      tabItem(
        "pk",
        fluidRow(
          box(width = 12, title = "PPI plasma exposure and proton-pump occupancy (first 30 days)",
              status = "primary", solidHeader = TRUE, plotOutput("p_pk_short", height = 330))
        ),
        fluidRow(
          box(width = 12, title = "Acid blockade, COX activity, gastrin over the full horizon",
              status = "primary", solidHeader = TRUE, plotOutput("p_pd_long", height = 380))
        )
      ),

      ## ---- 3. Reflux & injury axis ----------------------------------------
      tabItem(
        "injury",
        fluidRow(
          box(width = 12, title = "Acid exposure time, bile exposure and effective bile toxicity",
              status = "danger", solidHeader = TRUE, plotOutput("p_exposure", height = 360),
              footer = "Note that BILE does not fall with acid suppression, and the effective (pH-corrected) bile toxicity RISES.")
        ),
        fluidRow(
          box(width = 12, title = "Injury, inflammation, tissue PGE2",
              status = "danger", solidHeader = TRUE, plotOutput("p_inflam", height = 340))
        )
      ),

      ## ---- 4. Metaplasia & growth -----------------------------------------
      tabItem(
        "meta",
        fluidRow(
          box(width = 12, title = "CDX2 drive, Prague M segment length, Ki-67, mutational drive",
              status = "success", solidHeader = TRUE, plotOutput("p_meta", height = 420))
        ),
        fluidRow(
          box(width = 12, title = "Mutational drive decomposition (injury × proliferation)",
              status = "success", solidHeader = TRUE, plotOutput("p_mut", height = 320))
        )
      ),

      ## ---- 5. Clonal evolution --------------------------------------------
      tabItem(
        "clonal",
        fluidRow(
          box(width = 12, title = "Clone fractions along the CDKN2A → TP53 → HGD cascade",
              status = "warning", solidHeader = TRUE, plotOutput("p_clonal", height = 420))
        ),
        fluidRow(
          box(width = 6, title = "Aneuploidy index", status = "warning", solidHeader = TRUE,
              plotOutput("p_aneu", height = 300)),
          box(width = 6, title = "Invasive EAC burden", status = "warning", solidHeader = TRUE,
              plotOutput("p_eac", height = 300))
        )
      ),

      ## ---- 6. Clinical endpoints ------------------------------------------
      tabItem(
        "endpoint",
        fluidRow(
          box(width = 12, title = "Cumulative incidence: EAC and the AspECT composite (HGD or EAC)",
              status = "primary", solidHeader = TRUE, plotOutput("p_cum", height = 360))
        ),
        fluidRow(
          box(width = 6, title = "Annualized hazards (%/yr)", status = "primary",
              solidHeader = TRUE, plotOutput("p_haz", height = 320)),
          box(width = 6, title = "Symptom score (heartburn / regurgitation, 0-10)",
              status = "primary", solidHeader = TRUE, plotOutput("p_sympt", height = 320))
        )
      ),

      ## ---- 7. Scenario comparison -----------------------------------------
      tabItem(
        "compare",
        fluidRow(
          box(width = 12, title = "Cumulative AspECT composite (HGD or EAC) by scenario",
              status = "primary", solidHeader = TRUE, plotOutput("p_cmp_cum", height = 380))
        ),
        fluidRow(
          box(width = 12, title = "Endpoint table at horizon, with ratio vs low-dose PPI (AspECT reference arm)",
              status = "primary", solidHeader = TRUE, DTOutput("tbl_compare"))
        ),
        fluidRow(
          box(width = 12, title = "Acid exposure, mutational drive and segment length by scenario",
              status = "primary", solidHeader = TRUE, plotOutput("p_cmp_multi", height = 360))
        )
      ),

      ## ---- 8. Ablation & recurrence ---------------------------------------
      tabItem(
        "ablation",
        fluidRow(
          box(width = 4, title = "Ablation course", status = "warning", solidHeader = TRUE,
              sliderInput("abl_start", "First RFA session (year)", 0, 10, 0.1, step = 0.1),
              sliderInput("abl_n", "Number of RFA sessions", 1, 6, 3, step = 1),
              sliderInput("abl_int", "Interval between sessions (weeks)", 4, 16, 8, step = 2),
              checkboxInput("abl_emr", "Precede with EMR of a visible lesion", FALSE),
              radioButtons("abl_reflux", "Post-ablation reflux control",
                           c("High-dose PPI" = "high", "Low-dose PPI" = "low",
                             "Antireflux surgery" = "surg", "None" = "none"),
                           selected = "high")),
          box(width = 8, title = "Segment length, neosquamous cover, dysplastic fractions",
              status = "warning", solidHeader = TRUE, plotOutput("p_abl", height = 430))
        ),
        fluidRow(
          box(width = 12, title = "Recurrence is a reflux-control problem, not an ablation-technique problem",
              status = "info", solidHeader = TRUE, plotOutput("p_abl_recur", height = 340),
              footer = "Post-ablation regrowth in this model is not a fitted recurrence rate: the segment's target length is still set by the patient's own CDX2 drive, so residual/buried glands regrow whenever acid AND bile exposure remain uncontrolled. Same ablation, 10-year risk 4.35% on high-dose PPI vs 14.85% with no acid control vs 4.80% after antireflux surgery. Move the 'First RFA session' slider to see the cost of delay: 3.8% at 3 years if prompt, 19.8% if deferred a year.")
        )
      ),

      ## ---- 9. PGx & safety -------------------------------------------------
      tabItem(
        "safety",
        fluidRow(
          box(width = 12, title = "CYP2C19 phenotype vs acid control on esomeprazole 20 mg OD",
              status = "primary", solidHeader = TRUE, plotOutput("p_cyp", height = 360))
        ),
        fluidRow(
          box(width = 6, title = "Upper-GI bleeding risk (%/yr)", status = "danger",
              solidHeader = TRUE, plotOutput("p_bleed", height = 300)),
          box(width = 6, title = "Gastrin fold-rise (long-term PPI exposure surrogate)",
              status = "danger", solidHeader = TRUE, plotOutput("p_gastrin", height = 300))
        ),
        fluidRow(
          box(width = 12, title = "Benefit vs harm at horizon", status = "info",
              solidHeader = TRUE, DTOutput("tbl_safety"))
        )
      ),

      ## ---- References ------------------------------------------------------
      tabItem(
        "refs",
        box(width = 12, title = "Key calibration anchors", status = "info", solidHeader = TRUE,
            HTML("<ul>
            <li><b>AspECT</b> — Jankowski 2018 <i>Lancet</i> 392:400. 2557 patients,
            2×2 factorial, median 8.9 y. Composite (death/EAC/HGD): high-dose PPI
            HR 0.73; aspirin HR 0.93 (NS alone); high-dose PPI + aspirin vs
            low-dose PPI HR 0.59.</li>
            <li><b>SURF</b> — Phoa 2014 <i>JAMA</i> 311:1209. Confirmed LGD, RFA vs
            surveillance: progression to HGD/EAC 1.5% vs 26.5% at 3 years.</li>
            <li><b>AIM Dysplasia</b> — Shaheen 2009 <i>NEJM</i> 360:2277. RFA:
            CE-D 90.5%, CE-IM 77.4%, stricture 6%.</li>
            <li><b>Hvid-Jensen 2011</b> <i>NEJM</i> 365:1375. NDBE → EAC
            1.2 per 1000 person-years (0.12%/yr).</li>
            <li><b>Duits 2015</b> <i>Gut</i> 64:700. Expert-confirmed LGD ~9%/yr;
            73% of community LGD downstaged.</li>
            <li><b>Weaver 2014</b> <i>Nat Genet</i> 46:837. Mutation ordering:
            CDKN2A early, TP53 at the HGD transition, SMAD4 EAC-restricted.</li>
            <li><b>Jenkins 2007</b> <i>Carcinogenesis</i> 28:136. Deoxycholate is
            genotoxic via ROS <i>at neutral pH</i> — the basis of K_PH_BILE.</li>
            <li><b>Haigh 2003</b> <i>Gastroenterology</i> 124:615. Gastrin drives
            Barrett proliferation through CCK2R — the basis of K_GAS_KI.</li>
            </ul>
            <p>The full annotated list (74 references) is in
            <code>be_references.md</code>.</p>
            <p><b>Disclaimer.</b> Educational and research use only. Not validated
            for clinical decision-making or regulatory submission.</p>"))
      )
    )
  )
)

## ============================================================================
## Server
## ============================================================================
server <- function(input, output, session) {

  inp <- reactive({
    list(archetype = input$archetype, base_aet = input$base_aet,
         base_bile = input$base_bile, obesity = input$obesity, smoke = input$smoke,
         hh = input$hh, age = input$age, cyp = input$cyp,
         weightloss = input$weightloss, wl_year = input$wl_year %||% 2,
         wl_mag = input$wl_mag %||% 0.5, baclofen = input$baclofen,
         alginate = input$alginate, be_len_override = input$be_len_override)
  })
  `%||%` <- function(a, b) if (is.null(a)) b else a

  sim_main <- reactive(run_sim(input$scenario, input$horizon, inp()))
  sim_ref  <- reactive(run_sim("No therapy (natural history)", input$horizon, inp()))

  sim_cmp <- reactive({
    req(length(input$compare_set) > 0)
    bind_rows(lapply(input$compare_set, run_sim, horizon_y = input$horizon, inp = inp()))
  })

  both <- reactive(bind_rows(sim_main(), sim_ref() %>% mutate(scenario = "No therapy (reference)")))

  last <- function(df) df %>% group_by(scenario) %>% slice_tail(n = 1) %>% ungroup()

  ## ---- 1. Profile -----------------------------------------------------------
  output$vb_aet <- renderValueBox({
    v <- last(sim_main())$AET_pct
    valueBox(sprintf("%.1f %%", v), "Acid exposure time on therapy",
             icon = icon("droplet"), color = if (v < 4) "green" else if (v < 8) "yellow" else "red")
  })
  output$vb_len <- renderValueBox({
    v <- last(sim_main())$Prague_M_cm
    valueBox(sprintf("%.1f cm", v), "Prague M at horizon",
             icon = icon("ruler"), color = if (v < 1) "green" else if (v < 3) "yellow" else "red")
  })
  output$vb_haz <- renderValueBox({
    v <- last(sim_main())$Hazard_comp_pct_yr
    valueBox(sprintf("%.2f %%/yr", v), "HGD-or-EAC hazard at horizon",
             icon = icon("triangle-exclamation"),
             color = if (v < 0.3) "green" else if (v < 2) "yellow" else "red")
  })
  output$vb_cum <- renderValueBox({
    v <- last(sim_main())$CumInc_comp_pct
    valueBox(sprintf("%.2f %%", v), sprintf("Cumulative HGD/EAC at %d y", input$horizon),
             icon = icon("chart-area"),
             color = if (v < 1) "green" else if (v < 10) "yellow" else "red")
  })

  output$tbl_profile <- renderDT({
    d <- sim_main()
    f <- d %>% slice(1); l <- d %>% slice_tail(n = 1)
    tibble(
      Variable = c("Acid exposure time (%)", "Acid blockade (fraction)",
                   "Bile exposure (index)", "Effective bile toxicity",
                   "Injury index", "Inflammation index", "Tissue PGE2 (pg/mg)",
                   "Gastrin (fold of baseline)", "CDX2 drive", "Prague M (cm)",
                   "Ki-67 (%)", "Mutational drive (rel. to untreated)",
                   "p16-null fraction", "LGD fraction", "TP53-mutant fraction",
                   "HGD fraction", "Aneuploidy", "Symptom score (0-10)",
                   "HGD/EAC hazard (%/yr)", "Cumulative HGD/EAC (%)"),
      Baseline = round(c(f$AET_pct, f$Acid_blockade, f$Bile_exposure, f$Bile_effective,
                         f$Injury_index, f$Inflammation, f$PGE2_tissue, f$Gastrin_fold,
                         f$CDX2_drive, f$Prague_M_cm, f$Ki67_pct, f$Mutational_drive,
                         f$p16_null_frac, f$LGD_frac, f$TP53_mut_frac, f$HGD_frac,
                         f$Aneuploidy, f$Symptom_score, f$Hazard_comp_pct_yr,
                         f$CumInc_comp_pct), 3),
      `At horizon` = round(c(l$AET_pct, l$Acid_blockade, l$Bile_exposure, l$Bile_effective,
                             l$Injury_index, l$Inflammation, l$PGE2_tissue, l$Gastrin_fold,
                             l$CDX2_drive, l$Prague_M_cm, l$Ki67_pct, l$Mutational_drive,
                             l$p16_null_frac, l$LGD_frac, l$TP53_mut_frac, l$HGD_frac,
                             l$Aneuploidy, l$Symptom_score, l$Hazard_comp_pct_yr,
                             l$CumInc_comp_pct), 3)
    ) %>% datatable(rownames = FALSE, options = list(pageLength = 20, dom = "t"))
  })

  ## ---- 2. PK / PD -----------------------------------------------------------
  output$p_pk_short <- renderPlot({
    mod <- get_model(); d <- 30
    arch <- ARCHETYPES[[input$archetype]]
    bind_rows(
      mod %>% param(c(arch, list(CL2C19 = as.numeric(CYP_PHENO[[input$cyp]])))) %>%
        mrgsim(events = ppi_low(d), end = d, delta = 0.02) %>% as_tibble() %>%
        mutate(scenario = "Esomeprazole 20 mg OD"),
      mod %>% param(c(arch, list(CL2C19 = as.numeric(CYP_PHENO[[input$cyp]])))) %>%
        mrgsim(events = ppi_high(d), end = d, delta = 0.02) %>% as_tibble() %>%
        mutate(scenario = "Esomeprazole 40 mg BID")
    ) %>% mutate(years = time) %>%
      multi_plot(c("PPI_conc", "Pump_active_frac", "Acid_blockade", "AET_pct"),
                 c("PPI plasma (au)", "Active pump fraction",
                   "Acid blockade (0-1)", "AET (% pH<4)"),
                 "Days (not years) — pump occupancy accumulates over ~5 days") +
      labs(x = "Days")
  })

  output$p_pd_long <- renderPlot({
    multi_plot(both(),
               c("Acid_blockade", "COX_activity", "Gastrin_fold", "pH4_holding_pct"),
               c("Acid blockade (0-1)", "Mucosal COX activity",
                 "Gastrin (fold)", "pH>4 holding time (%)"),
               "Pharmacodynamic steady states")
  })

  ## ---- 3. Injury axis -------------------------------------------------------
  output$p_exposure <- renderPlot({
    multi_plot(both(), c("AET_pct", "Bile_exposure", "Bile_effective", "DeMeester_proxy"),
               c("AET (% pH<4)", "Bile exposure (index)",
                 "Effective bile toxicity", "DeMeester score (proxy)"),
               "The two arms of the refluxate behave differently under acid suppression")
  })
  output$p_inflam <- renderPlot({
    multi_plot(both(), c("Injury_index", "Inflammation", "PGE2_tissue", "Ki67_pct"),
               c("Injury / 8-OHdG index", "Inflammation (NF-κB composite)",
                 "Tissue PGE2 (pg/mg)", "Ki-67 (%)"),
               "Injury → inflammation → PGE2 → proliferation")
  })

  ## ---- 4. Metaplasia --------------------------------------------------------
  output$p_meta <- renderPlot({
    multi_plot(both(), c("CDX2_drive", "Prague_M_cm", "Ki67_pct", "Mutational_drive"),
               c("CDX2 drive (0-1)", "Prague M (cm)", "Ki-67 (%)",
                 "Mutational drive (rel.)"),
               "Metaplastic program and its proliferative consequence")
  })
  output$p_mut <- renderPlot({
    line_plot(both(), "Mutational_drive", "MUT (1.0 = untreated baseline)",
              "Mutational drive = (injury / injury_ref) × (Ki-67 / Ki-67_ref)")
  })

  ## ---- 5. Clonal ------------------------------------------------------------
  output$p_clonal <- renderPlot({
    multi_plot(both(), c("p16_null_frac", "LGD_frac", "TP53_mut_frac", "HGD_frac"),
               c("p16-null fraction", "LGD fraction", "TP53-mutant fraction",
                 "HGD fraction"),
               "Stepwise clonal cascade within the metaplastic field")
  })
  output$p_aneu <- renderPlot(line_plot(both(), "Aneuploidy", "Index (0-1)", "Aneuploidy"))
  output$p_eac  <- renderPlot(line_plot(both(), "EAC_burden", "Burden units", "Invasive EAC"))

  ## ---- 6. Endpoints ---------------------------------------------------------
  output$p_cum <- renderPlot({
    multi_plot(both(), c("CumInc_EAC_pct", "CumInc_comp_pct"),
               c("Cumulative EAC (%)", "Cumulative HGD or EAC (%)"),
               "Cumulative incidence")
  })
  output$p_haz <- renderPlot({
    multi_plot(both(), c("Hazard_LGD_pct_yr", "Hazard_HGD_pct_yr",
                         "Hazard_EAC_pct_yr", "Hazard_comp_pct_yr"),
               c("LGD (%/yr)", "HGD (%/yr)", "EAC (%/yr)", "HGD or EAC (%/yr)"),
               "Instantaneous annualized hazards")
  })
  output$p_sympt <- renderPlot(line_plot(both(), "Symptom_score", "Score (0-10)", "Symptoms"))

  ## ---- 7. Comparison --------------------------------------------------------
  output$p_cmp_cum <- renderPlot({
    line_plot(sim_cmp(), "CumInc_comp_pct", "Cumulative HGD or EAC (%)",
              "Scenario comparison — AspECT composite endpoint")
  })
  output$p_cmp_multi <- renderPlot({
    multi_plot(sim_cmp(), c("AET_pct", "Mutational_drive", "Prague_M_cm", "Symptom_score"),
               c("AET (% pH<4)", "Mutational drive", "Prague M (cm)", "Symptom score"),
               "Mechanistic separation between scenarios")
  })
  output$tbl_compare <- renderDT({
    l <- last(sim_cmp())
    ref <- l %>% filter(grepl("^Esomeprazole 20 mg OD \\(AspECT", scenario)) %>%
      pull(CumInc_comp_pct)
    refv <- if (length(ref) == 1 && ref > 0) ref else NA_real_
    l %>% transmute(
      Scenario = scenario,
      `AET (%)` = round(AET_pct, 2),
      `Bile eff.` = round(Bile_effective, 1),
      `MUT` = round(Mutational_drive, 3),
      `Prague M (cm)` = round(Prague_M_cm, 2),
      `LGD frac` = round(LGD_frac, 4),
      `HGD frac` = round(HGD_frac, 4),
      `Cum EAC (%)` = round(CumInc_EAC_pct, 3),
      `Cum HGD/EAC (%)` = round(CumInc_comp_pct, 3),
      `Ratio vs low-dose PPI` = round(CumInc_comp_pct / refv, 3),
      `Symptoms` = round(Symptom_score, 2),
      `Bleed (%/yr)` = round(Bleed_risk_pct_yr, 3)
    ) %>% datatable(rownames = FALSE, options = list(pageLength = 12, scrollX = TRUE))
  })

  ## ---- 8. Ablation ----------------------------------------------------------
  sim_abl <- reactive({
    mod <- get_model(); d <- input$horizon * 365
    arch <- ARCHETYPES[[input$archetype]]
    base_par <- c(arch, list(BASE_AET = input$base_aet, BASE_BILE = input$base_bile,
                             OBESITY = input$obesity, SMOKE = as.numeric(input$smoke),
                             CL2C19 = as.numeric(CYP_PHENO[[input$cyp]])))
    acid_ev <- switch(input$abl_reflux, high = ppi_high(d), low = ppi_low(d),
                      surg = NULL_EV, none = NULL_EV)
    par <- base_par
    if (input$abl_reflux == "surg") par <- c(par, list(FUNDO = 1, FUNDO_START = 1))
    abl <- rfa_x(input$abl_start * 365, input$abl_n, input$abl_int * 7)
    if (input$abl_emr) abl <- abl + emr_1(max(input$abl_start * 365 - 35, 0))
    bind_rows(
      mod %>% param(par) %>% mrgsim(events = acid_ev + abl, end = d, delta = 5) %>%
        as_tibble() %>% mutate(scenario = "Ablation course", years = time / 365),
      mod %>% param(par) %>% mrgsim(events = acid_ev, end = d, delta = 5) %>%
        as_tibble() %>% mutate(scenario = "No ablation", years = time / 365)
    )
  })

  output$p_abl <- renderPlot({
    multi_plot(sim_abl(), c("Prague_M_cm", "Neosquamous_frac", "LGD_frac", "HGD_frac"),
               c("Prague M (cm)", "Neosquamous fraction", "LGD fraction", "HGD fraction"),
               "Eradication and what comes back")
  })
  output$p_abl_recur <- renderPlot({
    multi_plot(sim_abl(), c("Prague_M_cm", "Stricture_burden", "CumInc_comp_pct", "Bile_effective"),
               c("Prague M (cm)", "Stricture burden", "Cumulative HGD/EAC (%)",
                 "Effective bile toxicity"),
               "Durability depends on the residual reflux drive")
  })

  ## ---- 9. PGx & safety ------------------------------------------------------
  sim_cyp <- reactive({
    mod <- get_model(); d <- input$horizon * 365
    arch <- ARCHETYPES[[input$archetype]]
    bind_rows(lapply(names(CYP_PHENO), function(nm) {
      mod %>% param(c(arch, list(CL2C19 = as.numeric(CYP_PHENO[[nm]]),
                                 BASE_AET = input$base_aet, BASE_BILE = input$base_bile,
                                 OBESITY = input$obesity, SMOKE = as.numeric(input$smoke)))) %>%
        mrgsim(events = ppi_low(d), end = d, delta = 5) %>% as_tibble() %>%
        mutate(scenario = nm, years = time / 365)
    }))
  })

  output$p_cyp <- renderPlot({
    multi_plot(sim_cyp(), c("Acid_blockade", "AET_pct", "Symptom_score", "CumInc_comp_pct"),
               c("Acid blockade (0-1)", "AET (% pH<4)", "Symptom score",
                 "Cumulative HGD/EAC (%)"),
               "Standard-dose PPI delivers very different acid control by CYP2C19 phenotype")
  })
  output$p_bleed   <- renderPlot(line_plot(sim_cmp(), "Bleed_risk_pct_yr", "%/yr",
                                           "Upper-GI bleeding risk"))
  output$p_gastrin <- renderPlot(line_plot(sim_cmp(), "Gastrin_fold", "Fold of baseline",
                                           "Serum gastrin"))
  output$tbl_safety <- renderDT({
    last(sim_cmp()) %>% transmute(
      Scenario = scenario,
      `Cum HGD/EAC (%)` = round(CumInc_comp_pct, 3),
      `Bleeding (%/yr)` = round(Bleed_risk_pct_yr, 3),
      `Gastrin (fold)` = round(Gastrin_fold, 2),
      `PPI safety index` = round(PPI_safety_index, 3),
      `Stricture burden` = round(Stricture_burden, 3),
      `Symptoms (0-10)` = round(Symptom_score, 2)
    ) %>% datatable(rownames = FALSE, options = list(pageLength = 12, dom = "t", scrollX = TRUE))
  })
}

shinyApp(ui, server)
