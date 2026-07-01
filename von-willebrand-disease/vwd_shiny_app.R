# =====================================================================
# Von Willebrand Disease (VWD) — QSP Shiny Dashboard
# Author : Claude Code Routine (2026-07-01)
# Deps   : shiny, bslib, mrgsolve, dplyr, tidyr, ggplot2, plotly, DT
# Run    : shiny::runApp("vwd_shiny_app.R")
# =====================================================================

library(shiny)
library(bslib)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

source("vwd_mrgsolve_model.R", local = TRUE)   # brings `vwd_mod`

ui <- page_navbar(
  title = "Von Willebrand Disease QSP Explorer",
  theme = bs_theme(bootswatch = "flatly", primary = "#AD1457"),
  sidebar = sidebar(
    width = 340, open = "always",
    h5("Patient profile"),
    selectInput("genotype", "VWD type",
                c("Type 1 (mild, partial quantitative)" = "type1",
                  "Type 2A (multimerization defect)" = "type2a",
                  "Type 2B (GPIba gain-of-function)" = "type2b",
                  "Type 2N (defective FVIII binding)" = "type2n",
                  "Type 3 (severe, total deficiency)" = "type3",
                  "Acquired VWS (aortic stenosis/lymphoproliferative)" = "acquired")),
    sliderInput("wt", "Body weight (kg)", 40, 120, 70, step = 1),
    checkboxInput("pregnant", "3rd-trimester pregnancy (physiologic VWF/FVIII boost)", FALSE),
    hr(),
    h5("Therapy"),
    selectInput("regimen", "Regimen",
                c("Untreated (natural history)",
                  "DDAVP IV single dose (0.3 mcg/kg)",
                  "DDAVP intranasal, repeated (q12h x3d)",
                  "Recombinant VWF (vonicog alfa-like) QD",
                  "Plasma-derived VWF/FVIII concentrate QD",
                  "Tranexamic acid + hormonal therapy (menorrhagia)",
                  "Major surgery: PK-guided PD concentrate q12h")),
    checkboxInput("fluid_restrict", "Post-DDAVP fluid restriction (24h)", FALSE),
    checkboxInput("hormonal_on", "Hormonal therapy active (COC/LNG-IUS)", FALSE),
    sliderInput("days", "Simulation horizon (days)", 1, 30, 7, step = 1),
    actionButton("run", "Run simulation", class = "btn-primary w-100")
  ),
  nav_panel("① Patient & Overview",
    layout_columns(
      col_widths = c(4,4,4),
      value_box(title = "VWD type", value = textOutput("geno_label"), theme = "secondary"),
      value_box(title = "Baseline VWF:RCo", value = textOutput("base_rco"), theme = "secondary"),
      value_box(title = "Regimen", value = textOutput("reg_label"), theme = "primary")
    ),
    card(card_header("Mechanistic map (12-cluster QSP schematic)"), htmlOutput("schematic"))
  ),
  nav_panel("② Drug PK",
    layout_columns(col_widths = c(6,6),
      card(card_header("DDAVP / Tranexamic acid"), plotlyOutput("plt_pk_small")),
      card(card_header("Recombinant VWF / Plasma-derived concentrate (IU/dL equiv.)"), plotlyOutput("plt_pk_factor"))
    )
  ),
  nav_panel("③ Hemostatic biomarkers",
    layout_columns(col_widths = c(6,6),
      card(card_header("VWF:Ag / VWF:RCo / HMWM fraction"), plotlyOutput("plt_vwf")),
      card(card_header("FVIII:C / Platelet count / ADAMTS13 activity"), plotlyOutput("plt_fviii_plt"))
    )
  ),
  nav_panel("④ Clinical endpoints",
    layout_columns(col_widths = c(4,4,4),
      value_box(title = "Bleeding score @ end", value = textOutput("bleed_end"), theme = "warning"),
      value_box(title = "Menstrual loss (mL/cyc)", value = textOutput("mens_end"), theme = "info"),
      value_box(title = "Hemoglobin (g/dL)", value = textOutput("hb_end"), theme = "success")
    ),
    card(card_header("Bleeding score, menstrual & GI blood loss trajectories"), plotlyOutput("plt_clinical"))
  ),
  nav_panel("⑤ Scenario comparison",
    card(card_header("Compare all regimens (VWF:RCo & bleeding score)"), plotlyOutput("plt_compare")),
    card(card_header("Table of endpoints at simulation end"), DTOutput("tbl_compare"))
  ),
  nav_panel("⑥ Biomarkers & diagnostics",
    layout_columns(col_widths = c(6,6),
      card(card_header("VWFpp (acute WPB-release marker)"), plotlyOutput("plt_vwfpp")),
      card(card_header("VWF:RCo/Ag ratio proxy (qualitative-defect discriminator)"), plotlyOutput("plt_ratio"))
    )
  ),
  nav_panel("⑦ Safety",
    layout_columns(col_widths = c(6,6),
      card(card_header("Serum sodium (DDAVP hyponatremia risk)"), plotlyOutput("plt_na")),
      card(card_header("Thrombotic-risk index (factor overcorrection)"), plotlyOutput("plt_thromb"))
    )
  ),
  nav_panel("⑧ References",
    card(card_header("Key references"), htmlOutput("refs"))
  )
)

genotype_params <- function(g, pregnant) {
  base <- switch(g,
    "type1"    = list(VWFAG_BASE = 35, VWFRCO_BASE = 30, FVIII_BASE = 42, HMWM_BASE = 0.85),
    "type2a"   = list(VWFAG_BASE = 45, VWFRCO_BASE = 15, FVIII_BASE = 55, HMWM_BASE = 0.20),
    "type2b"   = list(VWFAG_BASE = 45, VWFRCO_BASE = 22, FVIII_BASE = 55, HMWM_BASE = 0.55,
                       TYPE2B_FLAG = 1, PLT_BASE = 130),
    "type2n"   = list(VWFAG_BASE = 60, VWFRCO_BASE = 55, FVIII_BASE = 12, HMWM_BASE = 0.90,
                       TYPE2N_FLAG = 1),
    "type3"    = list(VWFAG_BASE = 2,  VWFRCO_BASE = 1,  FVIII_BASE = 6,  HMWM_BASE = 0.02),
    "acquired" = list(VWFAG_BASE = 40, VWFRCO_BASE = 18, FVIII_BASE = 45, HMWM_BASE = 0.30,
                       ACQUIRED_FLAG = 1, ADAMTS13_BASE = 0.6)
  )
  if (pregnant) base$PREG_BOOST <- 2.5
  base
}

run_regimen <- function(regimen, days, wt, genotype, pregnant, fluid_restrict, hormonal_on) {
  tmax <- days * 24
  params <- genotype_params(genotype, pregnant)
  params$FLUID_RESTRICT <- as.numeric(fluid_restrict)
  params$HORMONAL_ON    <- as.numeric(hormonal_on)
  mod <- vwd_mod %>% param(params)

  events <- switch(regimen,
    "Untreated (natural history)" = ev(amt = 0),
    "DDAVP IV single dose (0.3 mcg/kg)" =
      ev(cmt = "DDAVP_DEPOT", amt = 0.3 * wt, time = 0),
    "DDAVP intranasal, repeated (q12h x3d)" =
      ev(cmt = "DDAVP_DEPOT", amt = 300, ii = 12, addl = 5, time = 0),
    "Recombinant VWF (vonicog alfa-like) QD" =
      ev(cmt = "RVWF_CP", amt = 50 * wt, ii = 24, addl = floor(tmax/24) - 1, time = 0),
    "Plasma-derived VWF/FVIII concentrate QD" =
      ev(cmt = "PDVWF_CP", amt = 50 * wt, ii = 24, addl = floor(tmax/24) - 1, time = 0),
    "Tranexamic acid + hormonal therapy (menorrhagia)" =
      ev(cmt = "TXA_DEPOT", amt = 1300, ii = 8, addl = floor(tmax/8) - 1, time = 0),
    "Major surgery: PK-guided PD concentrate q12h" =
      ev(cmt = "PDVWF_CP", amt = 60 * wt, ii = 12, addl = floor(tmax/12) - 1, time = 0)
  )

  if (regimen == "Tranexamic acid + hormonal therapy (menorrhagia)") mod <- mod %>% param(HORMONAL_ON = 1)

  mod %>% ev(events) %>% mrgsim(end = tmax, delta = 1) %>% as_tibble()
}

server <- function(input, output, session) {

  sim <- eventReactive(input$run, {
    run_regimen(input$regimen, input$days, input$wt, input$genotype,
                input$pregnant, input$fluid_restrict, input$hormonal_on)
  }, ignoreNULL = FALSE)

  cmp <- eventReactive(input$run, {
    regs <- c("Untreated (natural history)",
              "DDAVP IV single dose (0.3 mcg/kg)",
              "DDAVP intranasal, repeated (q12h x3d)",
              "Recombinant VWF (vonicog alfa-like) QD",
              "Plasma-derived VWF/FVIII concentrate QD",
              "Tranexamic acid + hormonal therapy (menorrhagia)",
              "Major surgery: PK-guided PD concentrate q12h")
    bind_rows(lapply(regs, function(r) {
      run_regimen(r, input$days, input$wt, input$genotype, input$pregnant,
                  input$fluid_restrict, input$hormonal_on) %>% mutate(regimen = r)
    }))
  }, ignoreNULL = FALSE)

  output$geno_label <- renderText({
    switch(input$genotype, type1="Type 1", type2a="Type 2A", type2b="Type 2B",
           type2n="Type 2N", type3="Type 3", acquired="Acquired VWS")
  })
  output$base_rco <- renderText(sprintf("%.0f IU/dL", genotype_params(input$genotype, input$pregnant)$VWFRCO_BASE))
  output$reg_label <- renderText(input$regimen)

  output$schematic <- renderUI({
    tags$div(
      tags$p("VWF gene defects cause quantitative (Type 1/3) or qualitative (Type 2A/2B/2M/2N) ",
             "deficiency, impairing platelet-collagen/GPIba adhesion and/or FVIII stabilization."),
      tags$p("Desmopressin (DDAVP) triggers acute Weibel-Palade body release via a V2-like endothelial ",
             "receptor. Replacement therapy (recombinant or plasma-derived VWF/FVIII) restores multimers ",
             "directly. Tranexamic acid blocks fibrinolysis; hormonal therapy reduces menstrual blood loss."),
      tags$img(src = "vwd_qsp_model.png", style = "max-width:100%;")
    )
  })

  output$plt_pk_small <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~DDAVP_CP, name = "DDAVP Cp (ng/mL equiv)") %>%
      add_lines(y = ~TXA_CP, name = "TXA Cp (mcg/mL)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Conc."))
  })

  output$plt_pk_factor <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~RVWF_CONC, name = "rVWF contribution (IU/dL)") %>%
      add_lines(y = ~PDVWF_CONC, name = "PD concentrate contribution (IU/dL)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "IU/dL equiv"))
  })

  output$plt_vwf <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~VWF_AG_TOTAL, name = "VWF:Ag (IU/dL)") %>%
      add_lines(y = ~VWF_RCO_TOTAL, name = "VWF:RCo (IU/dL)") %>%
      add_lines(y = ~HMWM_EFFECTIVE*100, name = "HMWM fraction (x100)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "IU/dL / a.u."))
  })

  output$plt_fviii_plt <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~FVIII_C_TOTAL, name = "FVIII:C (IU/dL)") %>%
      add_lines(y = ~PLT_COUNT, name = "Platelets (x10^9/L)") %>%
      add_lines(y = ~ADAMTS13_ACT*100, name = "ADAMTS13 activity (x100)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "value"))
  })

  output$bleed_end <- renderText(sprintf("%.1f", tail(sim()$BLEED_SCORE, 1)))
  output$mens_end  <- renderText(sprintf("%.0f", tail(sim()$MENS_LOSS, 1)))
  output$hb_end    <- renderText(sprintf("%.1f", tail(sim()$HB, 1)))

  output$plt_clinical <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~BLEED_SCORE, name = "Bleeding score (a.u.)") %>%
      add_lines(y = ~MENS_LOSS, name = "Menstrual loss (mL/cycle)") %>%
      add_lines(y = ~GI_LOSS, name = "GI loss (mL/day)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "value"))
  })

  output$plt_compare <- renderPlotly({
    d <- cmp()
    p1 <- plot_ly(d, x = ~time/24, y = ~VWF_RCO_TOTAL, color = ~regimen, type = "scatter", mode = "lines",
                  legendgroup = ~regimen, showlegend = TRUE)
    p2 <- plot_ly(d, x = ~time/24, y = ~BLEED_SCORE, color = ~regimen, type = "scatter", mode = "lines",
                  legendgroup = ~regimen, showlegend = FALSE)
    subplot(p1, p2, nrows = 2, shareX = TRUE) %>%
      layout(xaxis = list(title = "Day"))
  })

  output$tbl_compare <- renderDT({
    d <- cmp() %>% group_by(regimen) %>% slice_tail(n = 1) %>%
      transmute(Regimen = regimen,
                `VWF:RCo (IU/dL)` = round(VWF_RCO_TOTAL, 1),
                `FVIII:C (IU/dL)` = round(FVIII_C_TOTAL, 1),
                `Platelets` = round(PLT_COUNT, 0),
                `Bleeding score` = round(BLEED_SCORE, 2),
                `Hemoglobin` = round(HB, 2),
                `Thrombotic risk` = round(THROMB_RISK, 2))
    datatable(d, options = list(pageLength = 10))
  })

  output$plt_vwfpp <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24, y = ~VWFPP, type = "scatter", mode = "lines",
            name = "VWFpp (IU/dL equiv)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "IU/dL equiv"))
  })

  output$plt_ratio <- renderPlotly({
    d <- sim() %>% mutate(rco_ag_ratio = VWF_RCO_TOTAL / (VWF_AG_TOTAL + 1e-6))
    plot_ly(d, x = ~time/24, y = ~rco_ag_ratio, type = "scatter", mode = "lines",
            name = "VWF:RCo/Ag ratio") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "ratio (<0.7 suggests qualitative defect)"))
  })

  output$plt_na <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24, y = ~NA_SERUM, type = "scatter", mode = "lines",
            name = "Serum sodium (mmol/L)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "mmol/L"))
  })

  output$plt_thromb <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24, y = ~THROMB_RISK, type = "scatter", mode = "lines",
            name = "Thrombotic-risk index") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "a.u."))
  })

  output$refs <- renderUI({
    tags$div(
      tags$p(tags$b("Key guidelines & pivotal trials:")),
      tags$ul(
        tags$li("James PD, Connell NT, et al. Blood Adv. 2021 (PMID 33570651) — ASH/ISTH/NHF/WFH 2021 diagnosis guideline."),
        tags$li("Connell NT, Flood VH, et al. Blood Adv. 2021 (PMID 33570647) — ASH/ISTH/NHF/WFH 2021 management guideline."),
        tags$li("Leebeek FW, Eikenboom JC. N Engl J Med. 2016 (PMID 27959741) — comprehensive VWD review."),
        tags$li("Gill JC, et al. Blood. 2015 (PMID 26239086) — recombinant VWF (vonicog alfa) pivotal trial."),
        tags$li("Mannucci PM. Blood. 1997 (PMID 9326215) — DDAVP pharmacology, the first 20 years.")
      ),
      tags$p("See vwd_references.md for the full bibliography.")
    )
  })
}

shinyApp(ui, server)
