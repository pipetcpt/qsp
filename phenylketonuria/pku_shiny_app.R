## ============================================================================
## Phenylketonuria (PKU) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs: Patient/Genotype · Drug PK · Amino-Acid PD (Plasma/Brain) ·
##         Neurocognitive Endpoints · Scenario comparison · Immunogenicity &
##         Growth Biomarkers · Maternal PKU & Safety · References
## Dependencies: shiny, shinydashboard, mrgsolve, dplyr, tidyr, ggplot2, DT
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

## ---------- Lazy model load ----------
get_model <- function() {
  if (!exists(".PKU_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.PKU_MOD)) {
    assign(".PKU_MOD", mread_cache("pku_mrgsolve_model.R"), envir = .GlobalEnv)
  }
  .GlobalEnv$.PKU_MOD
}

SCENARIO_LIST <- c(
  "1. Untreated classic PKU (natural history)",
  "2. Diet-only, well-controlled classic PKU",
  "3. Diet + sapropterin (BH4-responsive moderate PKU)",
  "4. Sapropterin trial, non-responder (classic PKU)",
  "5. Pegvaliase induction -> maintenance (severe classic PKU)",
  "6. Pegvaliase with heightened immunogenicity",
  "7. Untreated mild hyperphenylalaninemia (MHP)",
  "8. Maternal PKU - poor pre-conception control",
  "9. Maternal PKU - optimized pre-conception control"
)

## ---------- Scenario -> params/events ----------
scenario_setup <- function(scenario, wt, age0) {
  switch(scenario,
    "1. Untreated classic PKU (natural history)" =
      list(par = list(GENOTYPE=0, WT=wt, DIET_PHE_MGKGD=50, TYR_DIET_MGKGD=30, AGE0=age0), ev = NULL),
    "2. Diet-only, well-controlled classic PKU" =
      list(par = list(GENOTYPE=0, WT=wt, DIET_PHE_MGKGD=12, TYR_DIET_MGKGD=70, AGE0=age0), ev = NULL),
    "3. Diet + sapropterin (BH4-responsive moderate PKU)" =
      list(par = list(GENOTYPE=1, BH4_RESPONSIVE=1, WT=wt, DIET_PHE_MGKGD=25, TYR_DIET_MGKGD=50, AGE0=age0),
           ev = ev(amt=150, cmt="GUT_BH4", ii=24, addl=3650)),
    "4. Sapropterin trial, non-responder (classic PKU)" =
      list(par = list(GENOTYPE=0, BH4_RESPONSIVE=0, WT=wt, DIET_PHE_MGKGD=12, TYR_DIET_MGKGD=70, AGE0=age0),
           ev = ev(amt=300, cmt="GUT_BH4", ii=24, addl=3650)),
    "5. Pegvaliase induction -> maintenance (severe classic PKU)" =
      list(par = list(GENOTYPE=0, WT=wt, DIET_PHE_MGKGD=45, TYR_DIET_MGKGD=70, AGE0=age0),
           ev = ev(amt=20, cmt="SC_PEG", ii=24, addl=3650)),
    "6. Pegvaliase with heightened immunogenicity" =
      list(par = list(GENOTYPE=0, WT=wt, DIET_PHE_MGKGD=45, TYR_DIET_MGKGD=70, AGE0=age0, KIN_ADA=0.05),
           ev = ev(amt=20, cmt="SC_PEG", ii=24, addl=3650)),
    "7. Untreated mild hyperphenylalaninemia (MHP)" =
      list(par = list(GENOTYPE=3, WT=wt, DIET_PHE_MGKGD=50, TYR_DIET_MGKGD=30, AGE0=age0), ev = NULL),
    "8. Maternal PKU - poor pre-conception control" =
      list(par = list(GENOTYPE=0, MATERNAL=1, WT=wt, DIET_PHE_MGKGD=45, TYR_DIET_MGKGD=40, AGE0=25), ev = NULL),
    "9. Maternal PKU - optimized pre-conception control" =
      list(par = list(GENOTYPE=0, MATERNAL=1, WT=wt, DIET_PHE_MGKGD=10, TYR_DIET_MGKGD=70, AGE0=25), ev = NULL)
  )
}

run_sim <- function(scenario, horizon_h, wt, age0) {
  mod <- get_model()
  s <- scenario_setup(scenario, wt, age0)
  m <- mod %>% param(s$par)
  ev_obj <- s$ev
  out <- if (is.null(ev_obj)) {
    m %>% mrgsim(end = horizon_h, delta = 24)
  } else {
    m %>% mrgsim(events = ev_obj, end = horizon_h, delta = 24)
  }
  out %>% as_tibble() %>% mutate(scenario = scenario)
}

## ---------- UI ----------
ui <- dashboardPage(
  dashboardHeader(title = "PKU QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient / genotype",       tabName = "patient", icon = icon("user")),
      menuItem("2. Drug PK",                  tabName = "pk",      icon = icon("syringe")),
      menuItem("3. Amino-acid PD",            tabName = "aa",      icon = icon("dna")),
      menuItem("4. Neurocognitive endpoints", tabName = "neuro",   icon = icon("brain")),
      menuItem("5. Scenario comparison",      tabName = "compare", icon = icon("layer-group")),
      menuItem("6. Immunogenicity / growth",  tabName = "safety2", icon = icon("shield-virus")),
      menuItem("7. Maternal PKU & safety",    tabName = "maternal",icon = icon("baby")),
      menuItem("8. References",               tabName = "refs",    icon = icon("book"))
    ),
    hr(),
    selectInput("scenario", "Scenario:", SCENARIO_LIST,
                selected = "2. Diet-only, well-controlled classic PKU"),
    sliderInput("horizon_y", "Simulation horizon (years):", 0.5, 18, 8, step = 0.5),
    sliderInput("wt",   "Body weight (kg):", 3, 90, 8),
    sliderInput("age0",  "Starting age (years):", 0, 40, 0.02, step = 0.1),
    actionButton("run", "Run simulation", icon = icon("play"),
                 style = "color:#fff;background:#1f6feb")
  ),
  dashboardBody(
    tabItems(
      tabItem("patient",
        fluidRow(
          box(width = 12, title = "Patient / genotype summary", status = "primary",
              solidHeader = TRUE, DTOutput("patient_table"), br(),
              p(strong("Genotype legend:"), " Classic PKU (residual PAH <1%) · ",
                "Moderate (1-5%) · Mild (5-10%) · Mild hyperphenylalaninemia ",
                "(MHP, >10%, generally benign, no treatment required). ",
                "BH4-responsive genotypes are typically missense/misfolding ",
                "variants rescuable by sapropterin cofactor/chaperone action; ",
                "null alleles are non-responsive.")
          )
        ),
        fluidRow(
          box(width = 12, title = "About this dashboard", status = "info",
              p("This dashboard runs the mrgsolve QSP model for phenylketonuria ",
                "(PKU), linking dietary phenylalanine intake, saturable hepatic ",
                "PAH/BH4 enzymatic clearance, saturable blood-brain-barrier LAT1 ",
                "amino-acid transport competition, and downstream cerebral ",
                "monoamine synthesis, myelination and neurocognitive endpoints."),
              p("Pick a scenario and patient profile in the left panel, then ",
                "press ", strong("Run simulation"), " to update all plots. ",
                "The ACMG/NIH lifelong target range is Phe 120-360 umol/L.")
          )
        )
      ),

      tabItem("pk",
        fluidRow(
          box(width = 6, title = "Sapropterin plasma concentration", plotOutput("pk_bh4_plot", height = 360)),
          box(width = 6, title = "Pegvaliase plasma concentration", plotOutput("pk_peg_plot", height = 360))
        )
      ),

      tabItem("aa",
        fluidRow(
          box(width = 6, title = "Plasma Phe / Tyr / Trp (umol/L)", plotOutput("plasma_aa_plot", height = 380)),
          box(width = 6, title = "Brain (BBB) Phe / Tyr / Trp indices", plotOutput("brain_aa_plot", height = 380))
        ),
        fluidRow(
          box(width = 12, title = "Plasma Phe vs. ACMG/NIH target band (120-360 umol/L)",
              status = "warning", solidHeader = TRUE, plotOutput("phe_band_plot", height = 320))
        )
      ),

      tabItem("neuro",
        fluidRow(
          box(width = 6, title = "Dopamine & serotonin synthesis indices", plotOutput("monoamine_plot", 340)),
          box(width = 6, title = "Myelin/white-matter integrity index", plotOutput("myelin_plot", 340))
        ),
        fluidRow(
          box(width = 6, title = "Cumulative IQ deficit (points)", plotOutput("iq_plot", 340)),
          box(width = 6, title = "Executive-function deficit index", plotOutput("execfx_plot", 340))
        )
      ),

      tabItem("compare",
        fluidRow(
          box(width = 12, title = "Scenario comparison panel", status = "warning", solidHeader = TRUE,
              p("Runs all nine single-stage built-in scenarios with the current ",
                "weight/age profile; press the button below."),
              actionButton("run_all", "Run all scenarios", icon = icon("rocket"),
                           style = "color:#fff;background:#0f5132"),
              br(), br(), plotOutput("compare_plot", height = 600)
          )
        ),
        fluidRow(
          box(width = 12, title = "Endpoint summary table (final time point)",
              status = "info", DTOutput("endpoint_table"))
        )
      ),

      tabItem("safety2",
        fluidRow(
          box(width = 6, title = "Anti-drug antibody (ADA) titer & tolerization", plotOutput("ada_plot", 340)),
          box(width = 6, title = "Growth Z-score deficit", plotOutput("growth_plot", 340))
        )
      ),

      tabItem("maternal",
        fluidRow(
          box(width = 12, title = "Maternal PKU fetal-risk flags", status = "danger", solidHeader = TRUE,
              p("Active only when the scenario sets ", code("MATERNAL=1"), ". Congenital ",
                "heart disease (CHD) risk rises sharply above periconceptional Phe ",
                "~900-1200 umol/L; microcephaly/IQ-deficit risk with sustained Phe ",
                ">600 umol/L through pregnancy (Rouse 2000 J Pediatr; Koch 2000 Mol Genet Metab)."),
              plotOutput("maternal_plot", height = 340)
          )
        )
      ),

      tabItem("refs",
        fluidRow(
          box(width = 12, title = "Key references", status = "primary", solidHeader = TRUE,
              p("See ", code("pku_references.md"), " in this directory for the full ",
                "curated list (30+ PubMed-linked references)."),
              tags$ul(
                tags$li("Levy HL, et al. 2007 Lancet — sapropterin PKU-001 phase 3 (PMID 17368981)."),
                tags$li("Trefz FK, et al. 2009 Pediatrics — sapropterin PKU-003 pediatric (PMID 19153186)."),
                tags$li("Thomas J, et al. 2018 Mol Genet Metab — pegvaliase PRISM-1 (PMID 29174875)."),
                tags$li("Longo N, et al. 2019 Genet Med — pegvaliase PRISM-2 long-term (PMID 30190611)."),
                tags$li("Vockley J, et al. 2014 Genet Med — ACMG PKU management guideline (PMID 24385074)."),
                tags$li("Koch R, et al. 2003 Pediatrics — Maternal PKU Collaborative Study (PMID 12456904)."),
                tags$li("Waisbren SE, et al. 2007 Mol Genet Metab — treatment-timing IQ meta-analysis (PMID 17532249)."),
                tags$li("Diamond A, et al. 1997 Child Development — prefrontal dopamine-dependent deficits (PMID 9299839).")
              )
          )
        )
      )
    )
  )
)

## ---------- Server ----------
server <- function(input, output, session) {

  results <- reactiveVal(NULL)
  all_results <- reactiveVal(NULL)

  observeEvent(input$run, {
    showNotification("Running mrgsolve simulation…", type = "message", duration = 1)
    results(run_sim(input$scenario, input$horizon_y * 8760, input$wt, input$age0))
  }, ignoreNULL = FALSE)

  observeEvent(input$run_all, {
    showNotification("Running 9 scenarios…", type = "message", duration = 1)
    out <- lapply(SCENARIO_LIST, function(sc) run_sim(sc, input$horizon_y * 8760, input$wt, input$age0))
    all_results(bind_rows(out))
  })

  # --- Patient table ---
  output$patient_table <- renderDT({
    tibble(
      Field = c("Weight (kg)", "Starting age (y)", "Scenario", "Horizon (y)"),
      Value = c(input$wt, input$age0, input$scenario, input$horizon_y)
    ) %>% datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # --- PK ---
  output$pk_bh4_plot <- renderPlot({
    df <- results(); shiny::req(df)
    ggplot(df, aes(time/8760, conc_sapropterin)) +
      geom_line(colour = "#1f6feb", linewidth = 0.8) +
      labs(x = "Age (years)", y = "Sapropterin plasma conc (mg/L)") +
      theme_minimal(base_size = 13)
  })
  output$pk_peg_plot <- renderPlot({
    df <- results(); shiny::req(df)
    ggplot(df, aes(time/8760, conc_pegvaliase)) +
      geom_line(colour = "#cc3344", linewidth = 0.8) +
      labs(x = "Age (years)", y = "Pegvaliase plasma conc (mg/L)") +
      theme_minimal(base_size = 13)
  })

  # --- Amino-acid PD ---
  output$plasma_aa_plot <- renderPlot({
    df <- results(); shiny::req(df)
    df %>% select(time, Plasma_Phe_umolL, Plasma_Tyr_umolL, Plasma_Trp_umolL) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/8760, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        labs(x = "Age (years)", y = "Plasma conc (umol/L)", colour = "") +
        theme_minimal(base_size = 12)
  })
  output$brain_aa_plot <- renderPlot({
    df <- results(); shiny::req(df)
    df %>% select(time, Brain_Phe_idx, Brain_Tyr_idx, Brain_Trp_idx) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time/8760, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        labs(x = "Age (years)", y = "Brain index (umol/L-equivalent)", colour = "") +
        theme_minimal(base_size = 12)
  })
  output$phe_band_plot <- renderPlot({
    df <- results(); shiny::req(df)
    ggplot(df, aes(time/8760, Plasma_Phe_umolL)) +
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=120, ymax=360), fill="#d4edda", alpha=0.03, inherit.aes=FALSE) +
      geom_hline(yintercept = c(120,360), lty = 2, colour = "#2e7d32") +
      geom_line(colour = "#212529", linewidth = 1) +
      labs(x = "Age (years)", y = "Plasma Phe (umol/L)",
           caption = "Green dashed lines = ACMG/NIH lifelong target band (120-360 umol/L)") +
      theme_minimal(base_size = 13)
  })

  # --- Neurocognitive ---
  output$monoamine_plot <- renderPlot({
    df <- results(); shiny::req(df)
    df %>% select(time, Dopamine_idx, Serotonin_idx) %>% pivot_longer(-time) %>%
      ggplot(aes(time/8760, value, colour = name)) +
        geom_line(linewidth = 0.9) + geom_hline(yintercept = 1.0, lty = 2, colour = "grey50") +
        labs(x = "Age (years)", y = "Synthesis-capacity index (1.0=normal)", colour = "") +
        theme_minimal(base_size = 12)
  })
  output$myelin_plot <- renderPlot({
    df <- results(); shiny::req(df)
    ggplot(df, aes(time/8760, Myelin_idx)) +
      geom_line(colour = "#6d4c41", linewidth = 1) +
      scale_y_continuous(limits = c(0,1)) +
      labs(x = "Age (years)", y = "Myelin/white-matter integrity (1.0=normal)") +
      theme_minimal(base_size = 13)
  })
  output$iq_plot <- renderPlot({
    df <- results(); shiny::req(df)
    ggplot(df, aes(time/8760, IQ_deficit_pts)) +
      geom_line(colour = "#9b1c1c", linewidth = 1) +
      labs(x = "Age (years)", y = "Cumulative IQ deficit (points)") +
      theme_minimal(base_size = 13)
  })
  output$execfx_plot <- renderPlot({
    df <- results(); shiny::req(df)
    ggplot(df, aes(time/8760, ExecFunction_deficit)) +
      geom_line(colour = "#e65100", linewidth = 1) +
      labs(x = "Age (years)", y = "Executive-function deficit index (0-100)") +
      theme_minimal(base_size = 13)
  })

  # --- Immunogenicity / growth ---
  output$ada_plot <- renderPlot({
    df <- results(); shiny::req(df)
    df %>% select(time, ADA_titer, Tolerization_frac) %>% pivot_longer(-time) %>%
      ggplot(aes(time/8760, value, colour = name)) +
        geom_line(linewidth = 0.9) +
        facet_wrap(~name, scales = "free_y") +
        labs(x = "Age (years)", y = NULL) +
        theme_minimal(base_size = 12) + theme(legend.position = "none")
  })
  output$growth_plot <- renderPlot({
    df <- results(); shiny::req(df)
    ggplot(df, aes(time/8760, GrowthZ_deficit)) +
      geom_line(colour = "#00695c", linewidth = 1) +
      labs(x = "Age (years)", y = "Growth Z-score deficit (SD units)") +
      theme_minimal(base_size = 13)
  })

  # --- Maternal ---
  output$maternal_plot <- renderPlot({
    df <- results(); shiny::req(df)
    df %>% select(time, Plasma_Phe_umolL, Maternal_CHD_risk_flag, Maternal_Microcephaly_risk_flag) %>%
      mutate(week = time/168) %>%
      ggplot(aes(week, Plasma_Phe_umolL)) +
        geom_hline(yintercept = c(600,900), lty = 2, colour = c("#f9a825","#c62828")) +
        geom_line(colour = "#4a148c", linewidth = 1) +
        labs(x = "Gestational week", y = "Maternal plasma Phe (umol/L)",
             caption = "Orange = microcephaly-risk threshold (600); red = CHD-risk threshold (900)") +
        theme_minimal(base_size = 13)
  })

  # --- Scenario comparison ---
  output$compare_plot <- renderPlot({
    df <- all_results(); shiny::req(df)
    ggplot(df, aes(time/8760, Plasma_Phe_umolL, colour = scenario)) +
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=120, ymax=360), fill="#d4edda", alpha=0.02, inherit.aes=FALSE) +
      geom_line(linewidth = 0.9) +
      labs(x = "Age (years)", y = "Plasma Phe (umol/L)", colour = "Scenario") +
      theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
      guides(colour = guide_legend(ncol = 2))
  })
  output$endpoint_table <- renderDT({
    df <- all_results(); shiny::req(df)
    df %>% group_by(scenario) %>% slice_tail(n = 1) %>%
      select(scenario, Plasma_Phe_umolL, Plasma_Tyr_umolL, IQ_deficit_pts,
             ExecFunction_deficit, Dopamine_idx, ACMG_band_flag) %>%
      mutate(across(where(is.numeric), ~round(.x, 2))) %>%
      datatable(rownames = FALSE)
  })
}

shinyApp(ui, server)
