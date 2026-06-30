# =============================================================================
# Cluster Headache (CH) QSP — Shiny dashboard
#
# 8 tabs:
#   1. Patient & Disease Profile
#   2. Drug PK (sumatriptan/zolmitriptan/verapamil/lithium/topiramate/galcanezumab)
#   3. Pathway PD (hypothalamic drive, CGRP/PACAP, pial tone)
#   4. Clinical Endpoints (attacks/week, responder %, pain-free 15 min)
#   5. Scenario Comparison (S0..S6)
#   6. Biomarkers (CGRP serum, VIP, hazard)
#   7. Safety (verapamil ECG, lithium tox, topiramate cognition)
#   8. References & Notes
#
# Requires:  shiny, mrgsolve, tidyverse, plotly, DT, shinyWidgets
# Driver model: ch_mrgsolve_model.R  (sourced from same directory)
# =============================================================================

if (!exists("ch_mod")) source("ch_mrgsolve_model.R")

library(shiny);  library(plotly); library(DT); library(shinyWidgets)
library(tidyverse); library(mrgsolve)

scenario_labels <- c(
  "S0" = "S0: No treatment",
  "S1" = "S1: O2 + Sumatriptan SC (acute)",
  "S2" = "S2: Verapamil 240 mg BID",
  "S3" = "S3: Verapamil + Lithium",
  "S4" = "S4: Galcanezumab 300 mg q4w",
  "S5" = "S5: Prednisone bridge + Verapamil",
  "S6" = "S6: GON block + Verapamil"
)

# ---- UI ---------------------------------------------------------------------
ui <- navbarPage(
  "Cluster Headache QSP", id = "ch_nav",

  tabPanel("1. Patient",
    sidebarLayout(
      sidebarPanel(width = 4,
        h4("Demographics"),
        numericInput("WT",  "Weight (kg)", 78, min = 40, max = 150),
        numericInput("CrCL","CrCL (mL/min)", 100, min = 15, max = 180),
        radioButtons("SEX","Sex", c("Male"=1,"Female"=0), inline = TRUE),
        radioButtons("SMOKER","Smoking", c("Yes"=1,"No"=0), inline = TRUE),
        radioButtons("CHRONIC","Chronic CH?", c("No (episodic)"=0,"Yes (chronic)"=1), inline = TRUE),
        sliderInput("BASE_HYPO","Baseline hypothalamic drive", 0, 0.3, 0.05, 0.01),
        sliderInput("BOUT_AMP","In-bout amplitude", 0.3, 1.0, 0.95, 0.05),
        actionButton("rerun","Run / Refresh", class="btn-primary")
      ),
      mainPanel(
        h3("Patient & disease summary"),
        verbatimTextOutput("pt_summary"),
        plotlyOutput("circ_plot", height = 320)
      )
    )
  ),

  tabPanel("2. Drug PK",
    sidebarLayout(
      sidebarPanel(width = 3,
        h4("Acute drugs"),
        numericInput("DOSE_SUMA","Sumatriptan SC dose (mg/attack)", 6, 0, 12, 1),
        numericInput("DOSE_ZOL","Zolmitriptan IN dose (mg/attack)", 5, 0, 10, 1),
        h4("Preventive drugs"),
        numericInput("DOSE_VERA","Verapamil per dose (mg, BID)", 240, 0, 480, 40),
        numericInput("DOSE_LI","Lithium per dose (mg, BID)",     300, 0, 600, 100),
        numericInput("DOSE_TOPI","Topiramate per dose (mg, BID)", 50, 0, 200, 25),
        numericInput("DOSE_GALCA","Galcanezumab SC (mg q4w)",    300, 0, 300, 100),
        numericInput("DOSE_PRED","Prednisone (mg/d × 5 d)",       60, 0, 80, 10)
      ),
      mainPanel(
        plotlyOutput("pk_acute", height = 300),
        plotlyOutput("pk_prev",  height = 320),
        plotlyOutput("pk_mab",   height = 300)
      )
    )
  ),

  tabPanel("3. Pathway PD",
    plotlyOutput("hypo_drive", height = 300),
    plotlyOutput("cgrp_pacap", height = 320),
    plotlyOutput("pial_tone",  height = 280)
  ),

  tabPanel("4. Clinical Endpoints",
    fluidRow(
      column(6, plotlyOutput("attacks_week")),
      column(6, plotlyOutput("hazard_ts"))
    ),
    h4("Trial-anchored endpoint table"),
    DTOutput("endpoint_tbl")
  ),

  tabPanel("5. Scenario Comparison",
    checkboxGroupInput("scenarios_sel","Scenarios", choices = names(scenario_labels),
                       selected = c("S0","S1","S2","S4"), inline = TRUE),
    plotlyOutput("scn_attacks", height = 360),
    plotlyOutput("scn_cumat",   height = 320)
  ),

  tabPanel("6. Biomarkers",
    plotlyOutput("biom_cgrp", height = 300),
    plotlyOutput("biom_pacap", height = 300),
    plotlyOutput("biom_pial",  height = 280)
  ),

  tabPanel("7. Safety",
    h4("Verapamil — predicted trough vs PR interval"),
    plotlyOutput("safe_vera", height = 320),
    h4("Lithium — trough vs toxicity band"),
    plotlyOutput("safe_li",   height = 300),
    h4("Galcanezumab — exposure window"),
    plotlyOutput("safe_galca", height = 300)
  ),

  tabPanel("8. References",
    includeMarkdown("ch_references.md")
  )
)

# ---- SERVER -----------------------------------------------------------------
server <- function(input, output, session) {

  run_sim <- eventReactive(input$rerun, ignoreNULL = FALSE, {
    pars <- list(
      WT        = input$WT,
      CrCL      = input$CrCL,
      SEX       = as.numeric(input$SEX),
      SMOKER    = as.numeric(input$SMOKER),
      CHRONIC   = as.numeric(input$CHRONIC),
      BASE_HYPO = input$BASE_HYPO,
      BOUT_AMP  = input$BOUT_AMP
    )

    weeks <- 12
    end_h <- 24*7*weeks
    base_mod <- ch_mod %>% param(pars)

    build_doses <- function(sc) {
      switch(sc,
        "S0" = NULL,
        "S1" = c(ev(time = 24+4, amt = input$DOSE_SUMA, cmt = "SUMA_DEPOT"),
                 ev(time = 24+4, amt = 0.78,            cmt = "O2_EFFECT", evid = 1)),
        "S2" = ev(amt = input$DOSE_VERA, cmt = "VERA_DEPOT", ii = 12, addl = 2*weeks*7 - 1),
        "S3" = c(ev(amt = input$DOSE_VERA, cmt="VERA_DEPOT", ii=12, addl=2*weeks*7-1),
                 ev(amt = input$DOSE_LI,   cmt="LI_DEPOT",   ii=12, addl=2*weeks*7-1)),
        "S4" = ev(amt = input$DOSE_GALCA, cmt = "GALCA_SC", ii = 24*28,
                  addl = floor(weeks/4)),
        "S5" = c(ev(amt = input$DOSE_PRED, cmt="PRED_DEPOT", ii=24, addl=4),
                 ev(amt = 40,               cmt="PRED_DEPOT", ii=24, addl=2, time=24*5),
                 ev(amt = 20,               cmt="PRED_DEPOT", ii=24, addl=2, time=24*8),
                 ev(amt = input$DOSE_VERA,  cmt="VERA_DEPOT", ii=12, addl=2*weeks*7-1)),
        "S6" = c(ev(time = 0, amt = 0.65, cmt = "GON_EFF", evid = 1),
                 ev(amt = input$DOSE_VERA, cmt="VERA_DEPOT", ii=12, addl=2*weeks*7-1))
      )
    }

    res <- map_dfr(names(scenario_labels), function(sc) {
      evx <- build_doses(sc)
      mod <- if (!is.null(evx)) base_mod %>% ev(evx) else base_mod
      out <- mod %>% mrgsim(end = end_h, delta = 2) %>% as_tibble()
      out$scenario <- sc
      out
    })
    res
  })

  output$pt_summary <- renderPrint({
    cat("Cluster headache subtype:",
        ifelse(input$CHRONIC=="1","Chronic","Episodic"), "\n")
    cat("Sex:", ifelse(input$SEX=="1","M","F"), " | Wt:", input$WT,"kg | CrCL:",
        input$CrCL,"mL/min\n")
    cat("Hypothalamic baseline drive:", input$BASE_HYPO, "\n")
  })

  output$circ_plot <- renderPlotly({
    h <- 0:23
    drive <- 1 + 0.35*cos(2*pi*(h-3)/24)
    p <- plot_ly(x=h, y=drive, type="scatter", mode="lines",
                 line=list(color="#C0392B", width=3)) %>%
      layout(title="24-h circadian gating of hypothalamic drive (3 AM peak)",
             xaxis=list(title="Hour of day"), yaxis=list(title="Relative drive"))
    p
  })

  output$pk_acute <- renderPlotly({
    df <- run_sim() %>% filter(scenario=="S1")
    plot_ly(df, x=~time/24) %>%
      add_lines(y=~ConcSuma, name="Sumatriptan SC (ng/mL)") %>%
      layout(title="Acute abortive PK", xaxis=list(title="Days"), yaxis=list(title="Conc (ng/mL)"))
  })

  output$pk_prev <- renderPlotly({
    df <- run_sim() %>% filter(scenario=="S3")
    plot_ly(df, x=~time/24) %>%
      add_lines(y=~ConcVera, name="Verapamil ng/mL") %>%
      add_lines(y=~ConcLi*1000, name="Lithium µEq/L") %>%
      layout(title="Preventive PK (Verapamil & Lithium)",
             xaxis=list(title="Days"), yaxis=list(title="Conc"))
  })

  output$pk_mab <- renderPlotly({
    df <- run_sim() %>% filter(scenario=="S4")
    plot_ly(df, x=~time/24, y=~ConcGalca, type="scatter", mode="lines",
            line=list(color="#138D75", width=2)) %>%
      layout(title="Galcanezumab µg/mL — monthly 300 mg SC",
             xaxis=list(title="Days"), yaxis=list(title="µg/mL"))
  })

  output$hypo_drive <- renderPlotly({
    df <- run_sim() %>% filter(scenario %in% c("S0","S2","S4"))
    plot_ly(df, x=~time/24, y=~HypoDrive, color=~scenario, type="scatter", mode="lines") %>%
      layout(title="Hypothalamic drive vs treatment",
             xaxis=list(title="Days"), yaxis=list(title="HYPO_DRIVE (a.u.)"))
  })

  output$cgrp_pacap <- renderPlotly({
    df <- run_sim() %>% filter(scenario %in% c("S0","S4"))
    p1 <- plot_ly(df, x=~time/24, y=~CGRPtone, color=~scenario, type="scatter", mode="lines") %>%
      layout(title="CGRP tone (a.u.)", xaxis=list(title="Days"))
    p1
  })

  output$pial_tone <- renderPlotly({
    df <- run_sim()
    plot_ly(df, x=~time/24, y=~Pialtone, color=~scenario, type="scatter", mode="lines") %>%
      layout(title="Pial vasodilation tone",
             xaxis=list(title="Days"), yaxis=list(title="PIAL (a.u.)"))
  })

  output$attacks_week <- renderPlotly({
    df <- run_sim() %>% filter(scenario %in% input$scenarios_sel)
    plot_ly(df, x=~time/24, y=~AttacksWeek, color=~scenario, type="scatter", mode="lines") %>%
      layout(title="Attacks per week",
             xaxis=list(title="Days"), yaxis=list(title="Attacks/wk"))
  })

  output$hazard_ts <- renderPlotly({
    df <- run_sim() %>% filter(scenario %in% input$scenarios_sel)
    plot_ly(df, x=~time/24, y=~HazardPerH, color=~scenario, type="scatter", mode="lines") %>%
      layout(title="Instantaneous attack hazard /h",
             xaxis=list(title="Days"), yaxis=list(title="Hazard /h"))
  })

  output$endpoint_tbl <- renderDT({
    df <- run_sim() %>%
      group_by(scenario) %>%
      summarise(mean_atks_wk    = mean(AttacksWeek),
                cum_atks_12wk   = max(CUM_ATTACKS),
                .groups="drop") %>%
      mutate(label = scenario_labels[scenario])
    datatable(df[, c("label","mean_atks_wk","cum_atks_12wk")],
              options = list(pageLength = 7, dom='t'),
              colnames = c("Scenario","Mean attacks/wk","Cumulative attacks (12 wk)")) %>%
      formatRound(2:3, 1)
  })

  output$scn_attacks <- renderPlotly({
    df <- run_sim() %>% filter(scenario %in% input$scenarios_sel)
    plot_ly(df, x=~time/24, y=~AttacksWeek, color=~scenario, type="scatter", mode="lines") %>%
      layout(title="Attacks/week — scenarios", xaxis=list(title="Days"))
  })

  output$scn_cumat <- renderPlotly({
    df <- run_sim() %>% filter(scenario %in% input$scenarios_sel)
    plot_ly(df, x=~time/24, y=~CUM_ATTACKS, color=~scenario, type="scatter", mode="lines") %>%
      layout(title="Cumulative attacks", xaxis=list(title="Days"))
  })

  output$biom_cgrp <- renderPlotly({
    df <- run_sim() %>% filter(scenario %in% input$scenarios_sel)
    plot_ly(df, x=~time/24, y=~CGRPtone, color=~scenario, type="scatter", mode="lines") %>%
      layout(title="Serum CGRP surrogate", xaxis=list(title="Days"))
  })
  output$biom_pacap <- renderPlotly({
    df <- run_sim() %>% filter(scenario %in% input$scenarios_sel)
    plot_ly(df, x=~time/24, y=~PACAPtone, color=~scenario, type="scatter", mode="lines") %>%
      layout(title="PACAP tone", xaxis=list(title="Days"))
  })
  output$biom_pial <- renderPlotly({
    df <- run_sim() %>% filter(scenario %in% input$scenarios_sel)
    plot_ly(df, x=~time/24, y=~Pialtone, color=~scenario, type="scatter", mode="lines") %>%
      layout(title="Pial tone", xaxis=list(title="Days"))
  })

  output$safe_vera <- renderPlotly({
    df <- run_sim() %>% filter(scenario %in% c("S2","S3","S5","S6"))
    plot_ly(df, x=~time/24, y=~ConcVera, color=~scenario, type="scatter", mode="lines") %>%
      layout(title="Verapamil ng/mL  (PR > 200 ms risk if peaks >250 ng/mL)",
             xaxis=list(title="Days"), yaxis=list(title="ng/mL"))
  })

  output$safe_li <- renderPlotly({
    df <- run_sim() %>% filter(scenario == "S3")
    plot_ly(df, x=~time/24, y=~ConcLi, type="scatter", mode="lines",
            line=list(color="#7B3FBF", width=2)) %>%
      layout(title="Lithium trough (target 0.6-0.8 mEq/L; >1.2 toxicity)",
             xaxis=list(title="Days"), yaxis=list(title="mEq/L"),
             shapes = list(
               list(type="rect", fillcolor="#C8E6C9", line=list(width=0), opacity=0.4,
                    x0=0, x1=24*7*12/24, y0=0.6, y1=0.8),
               list(type="rect", fillcolor="#FFCCBC", line=list(width=0), opacity=0.3,
                    x0=0, x1=24*7*12/24, y0=1.2, y1=2.5)
             ))
  })

  output$safe_galca <- renderPlotly({
    df <- run_sim() %>% filter(scenario == "S4")
    plot_ly(df, x=~time/24, y=~ConcGalca, type="scatter", mode="lines",
            line=list(color="#138D75", width=2)) %>%
      layout(title="Galcanezumab µg/mL — chronic exposure",
             xaxis=list(title="Days"), yaxis=list(title="µg/mL"))
  })
}

# ---- Run --------------------------------------------------------------------
# shinyApp(ui = ui, server = server)
