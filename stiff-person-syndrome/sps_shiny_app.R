## =============================================================
## Stiff Person Syndrome (SPS) QSP — Shiny app skeleton
## Tabs (8):
##   1. Patient profile (anthropometrics, severity, comorbidity)
##   2. Anti-GAD65 antibody PK/dynamics (serum & CSF, IVIG, PLEX)
##   3. Drug PK (Diazepam/DMD, Baclofen oral/IT, Gabapentin,
##      Rituximab, Prednisolone)
##   4. CNS GABAergic biomarkers (GAD enzyme, GABA pool, MN exc.)
##   5. Clinical endpoints (HSI, stiffness, spasm freq, falls)
##   6. Treatment scenario comparison (overlay up to 4 regimens)
##   7. Safety / risk panel (BMD loss, infection, BZD tolerance)
##   8. Population sensitivity & VPC-style envelopes
## =============================================================

library(shiny)
library(bslib)
library(plotly)
library(dplyr)
library(tidyr)
library(ggplot2)

source("sps_mrgsolve_model.R")   # loads sps_mod and sps_run()

severity_levels <- c("mild", "moderate", "severe", "crisis")
scenarios <- c("dx_bzd_only",
               "bzd_baclofen_combo",
               "bzd_ivig_q4w",
               "rtx_induction",
               "plex_rescue",
               "intrathecal_baclofen")
scen_labels <- c("Diazepam mono",
                 "BZD + oral baclofen",
                 "BZD + IVIG q4w",
                 "Rituximab induction",
                 "PLEX rescue",
                 "Intrathecal baclofen")
names(scenarios) <- scen_labels

ui <- page_navbar(
  title = "SPS-QSP Simulator",
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  ## ============ TAB 1: Patient profile ====================
  nav_panel("1. Patient profile",
    layout_sidebar(
      sidebar = sidebar(
        h4("Patient"),
        numericInput("wt", "Body weight (kg)", 70, min = 30, max = 150),
        numericInput("age", "Age (years)", 45, min = 18, max = 85),
        selectInput("severity", "Disease severity",
                    choices = severity_levels, selected = "moderate"),
        selectInput("phenotype", "Phenotype",
                    choices = c("Classic SPS", "Stiff-limb (focal)",
                                "PERM (GlyR+)", "Paraneoplastic",
                                "SPS-plus / cerebellar")),
        checkboxGroupInput("comorbid", "Comorbidities",
                           choices = c("Type 1 diabetes", "Hashimoto",
                                       "Vitiligo", "Pernicious anemia",
                                       "Anxiety / agoraphobia")),
        actionButton("apply", "Apply patient", class = "btn-primary")
      ),
      card(card_header("Estimated baseline biomarkers"),
        tableOutput("baseline_tbl"),
        verbatimTextOutput("profile_note"))
    )
  ),

  ## ============ TAB 2: Antibody dynamics ==================
  nav_panel("2. Anti-GAD65 dynamics",
    layout_sidebar(
      sidebar = sidebar(
        h4("Immunomodulation"),
        sliderInput("ivig_g", "IVIG total dose per cycle (g)",
                    0, 200, 140, step = 10),
        sliderInput("ivig_cycles", "Number of cycles", 0, 12, 6),
        sliderInput("ivig_ii", "Cycle interval (weeks)", 2, 8, 4),
        checkboxInput("plex_on", "Add PLEX rescue (5 sessions)", FALSE),
        sliderInput("horizon", "Horizon (days)", 30, 365, 180, step = 30)
      ),
      card(card_header("Serum and CSF anti-GAD65"),
        plotlyOutput("plot_ab"),
        helpText("Serum anti-GAD65 in U/mL; CSF in U/mL (right axis).")),
      card(card_header("B-cell / plasma-cell tracks"),
        plotlyOutput("plot_bcell"))
    )
  ),

  ## ============ TAB 3: Drug PK ============================
  nav_panel("3. Drug PK",
    layout_sidebar(
      sidebar = sidebar(
        h4("Symptomatic regimen"),
        sliderInput("diaz_daily", "Diazepam daily (mg)", 0, 80, 30),
        sliderInput("bac_daily",  "Baclofen oral daily (mg)", 0, 120, 0),
        sliderInput("bac_IT",     "Intrathecal baclofen (mcg/d)", 0, 1000, 0),
        sliderInput("gab_daily",  "Gabapentin daily (mg)", 0, 3600, 0),
        sliderInput("pred_daily", "Prednisolone daily (mg)", 0, 100, 0),
        sliderInput("rtx_dose",   "Rituximab cycle dose (mg)", 0, 2000, 0),
        sliderInput("rtx_cycles", "Rituximab doses (n, q14d)", 0, 4, 2)
      ),
      card(card_header("Plasma & brain concentrations"),
        plotlyOutput("plot_pk", height = "500px"))
    )
  ),

  ## ============ TAB 4: GABAergic biomarkers ===============
  nav_panel("4. GABAergic biomarkers",
    layout_sidebar(
      sidebar = sidebar(
        h4("Trace overlay"),
        checkboxGroupInput("gaba_tracks", "Tracks",
                           choices = c("GAD activity (%)", "GABA pool (%)",
                                       "α-MN excitability (%)",
                                       "Total inhibition (×100)"),
                           selected = c("GAD activity (%)",
                                        "GABA pool (%)",
                                        "α-MN excitability (%)")),
        helpText("Driver: anti-GAD65 CSF → GAD enzyme inhibition")
      ),
      card(card_header("GAD / GABA / α-MN time course"),
        plotlyOutput("plot_gabaergic", height = "500px"))
    )
  ),

  ## ============ TAB 5: Clinical endpoints =================
  nav_panel("5. Clinical endpoints",
    layout_sidebar(
      sidebar = sidebar(
        h4("Endpoints"),
        helpText("Outputs reflect simulated patient response."),
        sliderInput("ep_horizon", "Horizon (days)", 30, 365, 180, step = 30)
      ),
      layout_columns(
        card(card_header("Heightened Sensitivity Index (HSI)"),
             plotlyOutput("plot_hsi")),
        card(card_header("Spasm frequency (events/d)"),
             plotlyOutput("plot_spasm"))
      ),
      layout_columns(
        card(card_header("Stiffness score (0-100)"),
             plotlyOutput("plot_stiff")),
        card(card_header("Falls projection"),
             plotlyOutput("plot_falls"))
      )
    )
  ),

  ## ============ TAB 6: Scenario comparison ================
  nav_panel("6. Scenario comparison",
    layout_sidebar(
      sidebar = sidebar(
        h4("Overlay up to 4 scenarios"),
        checkboxGroupInput("scen_pick", "Scenarios",
                           choices = scenarios,
                           selected = scenarios[c(1, 3, 4)]),
        selectInput("scen_severity", "Severity",
                    choices = severity_levels, selected = "moderate"),
        sliderInput("scen_horizon", "Horizon (days)", 30, 365, 180, step = 30)
      ),
      layout_columns(
        card(card_header("Stiffness trajectory"),
             plotlyOutput("plot_scen_stiff")),
        card(card_header("Anti-GAD65 trajectory"),
             plotlyOutput("plot_scen_ab"))
      ),
      card(card_header("Endpoint summary at horizon"),
           tableOutput("scen_table"))
    )
  ),

  ## ============ TAB 7: Safety / risk panel ================
  nav_panel("7. Safety & risk",
    layout_sidebar(
      sidebar = sidebar(
        h4("Risk look-ups"),
        sliderInput("safety_horizon", "Horizon (months)", 1, 24, 12),
        sliderInput("bzd_high", "High-dose BZD threshold (mg/d)", 20, 80, 50)
      ),
      layout_columns(
        card(card_header("BMD trajectory (lumbar)"),
             plotlyOutput("plot_bmd")),
        card(card_header("Cumulative steroid exposure"),
             plotlyOutput("plot_steroid"))
      ),
      card(card_header("Adverse-event projection"),
        tableOutput("safety_tbl"))
    )
  ),

  ## ============ TAB 8: Population sensitivity =============
  nav_panel("8. Population sensitivity",
    layout_sidebar(
      sidebar = sidebar(
        h4("Virtual cohort"),
        sliderInput("n_pop", "Cohort size (N)", 20, 500, 100, step = 20),
        sliderInput("ab_cv", "Anti-GAD CV%", 0, 60, 30),
        sliderInput("gad_cv", "GAD activity CV%", 0, 60, 25),
        sliderInput("pk_cv", "PK CL CV%", 0, 60, 30),
        actionButton("run_pop", "Simulate cohort", class = "btn-primary")
      ),
      card(card_header("VPC envelope — stiffness"),
        plotlyOutput("plot_vpc")),
      card(card_header("Responder analysis (≥30% HSI drop)"),
        tableOutput("vpc_tbl"))
    )
  )
)

server <- function(input, output, session){

  ## ----- baseline preview ------------------------------------
  output$baseline_tbl <- renderTable({
    sev <- input$severity %||% "moderate"
    vals <- switch(sev,
      mild      = c("Anti-GAD65 (U/mL)" = 4000,  "GAD (%)" = 80, "Stiffness" = 25, "Spasms/d" = 0.5),
      moderate  = c("Anti-GAD65 (U/mL)" = 10000, "GAD (%)" = 55, "Stiffness" = 55, "Spasms/d" = 1.5),
      severe    = c("Anti-GAD65 (U/mL)" = 25000, "GAD (%)" = 35, "Stiffness" = 80, "Spasms/d" = 4),
      crisis    = c("Anti-GAD65 (U/mL)" = 40000, "GAD (%)" = 25, "Stiffness" = 92, "Spasms/d" = 7))
    data.frame(Marker = names(vals), Value = unname(vals))
  })

  output$profile_note <- renderText({
    paste("Phenotype:", input$phenotype,
          "\nComorbidities:", paste(input$comorbid, collapse = ", "))
  })

  ## ----- Tab 2: antibody dynamics ----------------------------
  ab_sim <- reactive({
    sps_run(scenario = if(input$plex_on) "plex_rescue" else "bzd_ivig_q4w",
            severity = input$severity %||% "moderate",
            horizon_d = input$horizon)
  })
  output$plot_ab <- renderPlotly({
    d <- ab_sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~ABAB_S,   name = "Serum anti-GAD65 (U/mL)") %>%
      add_lines(y = ~ABAB_CSF*100, name = "CSF anti-GAD65 (U/mL × 100)",
                line = list(dash = "dash")) %>%
      layout(xaxis = list(title = "Days"),
             yaxis = list(title = "Anti-GAD65"))
  })
  output$plot_bcell <- renderPlotly({
    d <- ab_sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~BNAIVE, name = "Naive B (cells/uL)") %>%
      add_lines(y = ~PBLAST*10, name = "Plasmablasts ×10") %>%
      add_lines(y = ~LLPC*10,   name = "LLPC ×10") %>%
      layout(xaxis = list(title = "Days"),
             yaxis = list(title = "Count / rel U"))
  })

  ## ----- Tab 3: PK ------------------------------------------
  pk_sim <- reactive({
    # custom regimen, run mrgsolve directly
    mod <- sps_mod
    if(!is.null(input$severity)) mod <- sps_active_init(mod, input$severity)
    ev_in <- ev(amt = 0, cmt = "DIAZ_GUT", time = -1)
    if(input$diaz_daily > 0)
      ev_in <- c(ev_in, ev(amt = input$diaz_daily/3, ii = 8,
                            addl = input$horizon*3-1, cmt = "DIAZ_GUT"))
    if(input$bac_daily > 0)
      ev_in <- c(ev_in, ev(amt = input$bac_daily/3, ii = 8,
                            addl = input$horizon*3-1, cmt = "BAC_GUT"))
    if(input$bac_IT > 0)
      ev_in <- c(ev_in, ev(amt = input$bac_IT/1000, ii = 24,
                            addl = input$horizon-1, cmt = "BAC_CSF"))
    if(input$gab_daily > 0)
      ev_in <- c(ev_in, ev(amt = input$gab_daily/3, ii = 8,
                            addl = input$horizon*3-1, cmt = "GAB_GUT"))
    if(input$pred_daily > 0)
      ev_in <- c(ev_in, ev(amt = input$pred_daily, ii = 24,
                            addl = input$horizon-1, cmt = "PRED_GUT"))
    if(input$rtx_dose > 0 && input$rtx_cycles > 0)
      ev_in <- c(ev_in, ev(amt = input$rtx_dose, ii = 24*14,
                            addl = input$rtx_cycles-1, cmt = "RTX_C",
                            rate = input$rtx_dose/4))
    mrgsim(mod, ev = ev_in, end = (input$horizon %||% 90)*24, delta = 0.5) %>%
      as.data.frame()
  })
  output$plot_pk <- renderPlotly({
    d <- pk_sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~diaz_total_ngml, name = "Diazepam ng/mL") %>%
      add_lines(y = ~dmd_ngml,        name = "DMD ng/mL") %>%
      add_lines(y = ~bac_csf_ngml,    name = "Baclofen CSF ng/mL") %>%
      add_lines(y = ~gab_mgL*1000,    name = "Gabapentin ng/mL") %>%
      add_lines(y = ~rtx_ugml*1000,   name = "Rituximab ng/mL") %>%
      layout(xaxis = list(title = "Days"),
             yaxis = list(title = "Conc (ng/mL)", type = "log"))
  })

  ## ----- Tab 4: GABAergic -----------------------------------
  gaba_sim <- reactive({ pk_sim() })
  output$plot_gabaergic <- renderPlotly({
    d <- gaba_sim()
    pl <- plot_ly(d, x = ~time/24)
    if("GAD activity (%)" %in% input$gaba_tracks)
      pl <- add_lines(pl, y = ~gad_pct,  name = "GAD activity (%)")
    if("GABA pool (%)" %in% input$gaba_tracks)
      pl <- add_lines(pl, y = ~gaba_pct, name = "GABA pool (%)")
    if("α-MN excitability (%)" %in% input$gaba_tracks)
      pl <- add_lines(pl, y = ~mn_pct,   name = "α-MN excitability (%)")
    if("Total inhibition (×100)" %in% input$gaba_tracks)
      pl <- add_lines(pl, y = ~inh_total*100, name = "Total inhibition ×100")
    pl %>% layout(xaxis = list(title = "Days"), yaxis = list(title = "% / score"))
  })

  ## ----- Tab 5: clinical -------------------------------------
  ep_sim <- reactive({ pk_sim() })
  output$plot_hsi <- renderPlotly({
    d <- ep_sim();  plot_ly(d, x = ~time/24, y = ~HSI, name = "HSI") %>% add_lines()
  })
  output$plot_spasm <- renderPlotly({
    d <- ep_sim();  plot_ly(d, x = ~time/24, y = ~SPASM, name = "Spasm/d") %>% add_lines()
  })
  output$plot_stiff <- renderPlotly({
    d <- ep_sim();  plot_ly(d, x = ~time/24, y = ~STIFF, name = "Stiffness") %>% add_lines()
  })
  output$plot_falls <- renderPlotly({
    d <- ep_sim()
    d <- d %>% mutate(falls_per_y = pmin(120, 0.02*STIFF*SPASM))
    plot_ly(d, x = ~time/24, y = ~falls_per_y, name = "Falls/y") %>% add_lines()
  })

  ## ----- Tab 6: scenario comparison --------------------------
  scen_sim <- reactive({
    do.call(rbind,
      lapply(input$scen_pick, function(s){
        r <- sps_run(scenario = s,
                     severity = input$scen_severity,
                     horizon_d = input$scen_horizon)
        r$scenario <- names(scenarios)[scenarios == s]
        r
      }))
  })
  output$plot_scen_stiff <- renderPlotly({
    d <- scen_sim()
    plot_ly(d, x = ~time/24, y = ~STIFF, color = ~scenario) %>% add_lines()
  })
  output$plot_scen_ab <- renderPlotly({
    d <- scen_sim()
    plot_ly(d, x = ~time/24, y = ~ABAB_S, color = ~scenario) %>% add_lines()
  })
  output$scen_table <- renderTable({
    d <- scen_sim()
    d %>% group_by(scenario) %>%
      summarize(`Final stiffness`   = round(tail(STIFF, 1), 1),
                `Final spasm/d`     = round(tail(SPASM, 1), 2),
                `Anti-GAD (%base)`  = round(tail(abab_pct, 1), 1),
                `GABA (%base)`      = round(tail(gaba_pct, 1), 1))
  })

  ## ----- Tab 7: safety --------------------------------------
  safety_sim <- reactive({ pk_sim() })
  output$plot_bmd <- renderPlotly({
    d <- safety_sim(); plot_ly(d, x = ~time/24, y = ~BMD) %>% add_lines()
  })
  output$plot_steroid <- renderPlotly({
    d <- safety_sim()
    d <- d %>% mutate(cum_pred = cumsum(pmax(0, PRED_C_C %||% 0))*0.5)
    plot_ly(d, x = ~time/24, y = ~cum_pred) %>% add_lines()
  })
  output$safety_tbl <- renderTable({
    d <- safety_sim()
    data.frame(
      Risk = c("BZD daily > threshold (% time)",
               "Cumulative steroid (g equiv)",
               "BMD loss (%)",
               "Projected infection events/y"),
      Value = c(
        round(100*mean(d$diaz_total_ngml > input$bzd_high*5, na.rm=TRUE), 1),
        round(sum(pmax(0, d$PRED_C_C %||% 0))*0.5/1000, 2),
        round(100*(1 - tail(d$BMD,1)/1.0), 2),
        round(0.1 + 0.5 * (input$rtx_dose > 0), 2)
      )
    )
  })

  ## ----- Tab 8: VPC ------------------------------------------
  vpc_sim <- eventReactive(input$run_pop, {
    n <- input$n_pop
    set.seed(42)
    runs <- lapply(seq_len(n), function(i){
      mod_i <- sps_mod
      ab_mul   <- rlnorm(1, 0, input$ab_cv/100)
      gad_mul  <- rlnorm(1, 0, input$gad_cv/100)
      cl_mul   <- rlnorm(1, 0, input$pk_cv/100)
      mod_i <- param(mod_i,
                     CL_DIAZ = sps_mod$CL_DIAZ * cl_mul,
                     CL_BAC  = sps_mod$CL_BAC  * cl_mul,
                     CL_IGG  = sps_mod$CL_IGG  * cl_mul)
      mod_i <- sps_active_init(mod_i, "moderate")
      mod_i <- init(mod_i,
                    ABAB_S = pmin(80000, max(200, 10000 * ab_mul)),
                    GAD    = pmin(95, max(15, 55 * gad_mul)))
      out <- mrgsim(mod_i,
                    ev = c(ev(amt = 5, ii = 8, addl = 200, cmt = "DIAZ_GUT"),
                           ev(amt = 70, ii = 24*28, addl = 5,
                              cmt = "IVIG_C")),
                    end = 24*180, delta = 6) %>% as.data.frame()
      out$ID <- i
      out
    })
    do.call(rbind, runs)
  })
  output$plot_vpc <- renderPlotly({
    d <- vpc_sim()
    pct <- d %>% group_by(time) %>%
      summarize(p05 = quantile(STIFF, 0.05),
                p50 = quantile(STIFF, 0.50),
                p95 = quantile(STIFF, 0.95))
    plot_ly(pct, x = ~time/24) %>%
      add_ribbons(ymin = ~p05, ymax = ~p95, name = "5–95%",
                  line = list(width = 0),
                  fillcolor = "rgba(58,120,194,0.25)") %>%
      add_lines(y = ~p50, name = "Median") %>%
      layout(xaxis = list(title = "Days"),
             yaxis = list(title = "Stiffness"))
  })
  output$vpc_tbl <- renderTable({
    d <- vpc_sim()
    base <- d %>% filter(time == 0) %>% select(ID, base = STIFF)
    last <- d %>% group_by(ID) %>% summarize(last = tail(STIFF, 1))
    df <- inner_join(base, last, by = "ID")
    data.frame(
      Metric = c("N simulated",
                 "Median ΔStiffness",
                 "% responders (≥30% drop)",
                 "% with crisis (Stiff>90 sustained)"),
      Value  = c(nrow(df),
                 round(median(df$base - df$last), 1),
                 round(100*mean((df$base - df$last)/df$base >= 0.30), 1),
                 round(100*mean(df$last > 90), 1))
    )
  })

}

shinyApp(ui, server)
