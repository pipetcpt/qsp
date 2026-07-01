# =====================================================================
# Mucopolysaccharidosis Type I (MPS I) — QSP Shiny Dashboard
# Author : Claude Code Routine (2026-07-01)
# Deps   : shiny, bslib, mrgsolve, dplyr, tidyr, ggplot2, plotly, DT
# Run    : shiny::runApp("mps1_shiny_app.R")
# =====================================================================

library(shiny)
library(bslib)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

source("mps1_mrgsolve_model.R", local = TRUE)   # brings `mps1_mod`

REGIMENS <- c(
  "Untreated (natural history)",
  "Laronidase ERT 0.58 mg/kg IV weekly",
  "HSCT alone (early transplant)",
  "HSCT alone (delayed transplant)",
  "ERT bridging (~15 wk) then HSCT",
  "ERT, high anti-drug antibody titer",
  "ERT, poor adherence (60% of infusions)",
  "ERT + genistein (investigational SRT)",
  "ERT + AAV9 CNS gene therapy (investigational)",
  "HSCT + long-term low-dose ERT (residual disease)"
)

ui <- page_navbar(
  title = "Mucopolysaccharidosis Type I (MPS I) QSP Explorer",
  theme = bs_theme(bootswatch = "flatly", primary = "#00838F"),
  sidebar = sidebar(
    width = 340, open = "always",
    h5("Patient profile"),
    selectInput("phenotype", "Phenotype", c("Severe (Hurler, CNS+)" = "1", "Attenuated (Hurler-Scheie/Scheie)" = "0")),
    sliderInput("wt", "Body weight (kg)", 5, 40, 13, step = 1),
    sliderInput("age_transplant", "Age at HSCT (months, if applicable)", 3, 48, 12, step = 1),
    hr(),
    h5("Therapy"),
    selectInput("regimen", "Regimen", REGIMENS, selected = REGIMENS[1]),
    sliderInput("years", "Simulation horizon (yr)", 1, 5, 5, step = 1),
    actionButton("run", "Run simulation", class = "btn-primary w-100")
  ),
  nav_panel("① Patient & Overview",
    layout_columns(
      col_widths = c(4,4,4),
      value_box(title = "Phenotype", value = textOutput("pheno_label"), theme = "secondary"),
      value_box(title = "Regimen", value = textOutput("reg_label"), theme = "primary"),
      value_box(title = "Untreated median survival", value = "~6-10 yr (severe, CNS-driven)", theme = "secondary")
    ),
    card(card_header("Mechanistic map (11-cluster QSP schematic)"), htmlOutput("schematic"))
  ),
  nav_panel("② Drug PK",
    layout_columns(col_widths = c(6,6),
      card(card_header("Laronidase plasma vs. tissue-retained enzyme"), plotlyOutput("plt_pk_laro")),
      card(card_header("Anti-drug antibody (ADA) / genistein exposure"), plotlyOutput("plt_pk_other"))
    )
  ),
  nav_panel("③ Enzyme access / GAG burden",
    layout_columns(col_widths = c(6,6),
      card(card_header("Enzyme-access index by compartment (systemic/CNS/cartilage)"), plotlyOutput("plt_enzaccess")),
      card(card_header("GAG burden pools (normalized)"), plotlyOutput("plt_gag"))
    )
  ),
  nav_panel("④ Clinical endpoints",
    layout_columns(col_widths = c(3,3,3,3),
      value_box(title = "Urinary GAG @ end", value = textOutput("ugag_end"), theme = "success"),
      value_box(title = "DQ @ end", value = textOutput("dq_end"), theme = "success"),
      value_box(title = "Liver/spleen index @ end", value = textOutput("livspleen_end"), theme = "info"),
      value_box(title = "Survival probability @ end", value = textOutput("surv_end"), theme = "info")
    ),
    card(card_header("uGAG, liver/spleen, FVC, AHI, joint ROM, DQ trajectories"), plotlyOutput("plt_clinical"))
  ),
  nav_panel("⑤ Scenario comparison",
    card(card_header("Compare all pre-configured regimens (urinary GAG)"), plotlyOutput("plt_compare")),
    card(card_header("Table of endpoints at simulation end"), DTOutput("tbl_compare"))
  ),
  nav_panel("⑥ Biomarkers",
    layout_columns(col_widths = c(6,6),
      card(card_header("Cardiac valve/LV-mass index & corneal clouding"), plotlyOutput("plt_valve_cornea")),
      card(card_header("Height Z-score"), plotlyOutput("plt_height"))
    )
  ),
  nav_panel("⑦ Safety / Survival",
    layout_columns(col_widths = c(6,6),
      card(card_header("Cumulative mortality hazard"), plotlyOutput("plt_hazard")),
      card(card_header("Survival probability"), plotlyOutput("plt_survival"))
    )
  ),
  nav_panel("⑧ References",
    card(card_header("Key references"), htmlOutput("refs"))
  )
)

wt_mgkg <- function(wt, dose_mgkg) dose_mgkg * wt

run_regimen <- function(regimen, years, wt, phenotype, age_transplant) {
  tmax <- years * 8760
  laro_dose <- wt_mgkg(wt, 0.58)
  laro_full   <- ev(cmt = "LARO_CENT", amt = laro_dose, ii = 168, addl = floor(tmax/168) - 1, time = 0, rate = laro_dose/4)
  laro_bridge <- ev(cmt = "LARO_CENT", amt = laro_dose, ii = 168, addl = 14, time = 0, rate = laro_dose/4)
  gen_ev      <- ev(cmt = "GEN_GUT", amt = wt_mgkg(wt, 10), ii = 24, addl = floor(tmax/24) - 1, time = 0)

  pars <- list(PHENOTYPE = as.numeric(phenotype), HSCT_FLAG = 0, AAV9_FLAG = 0,
               AGE_AT_TRANSPLANT_MO = age_transplant, ADHERENCE_ERT = 1.0, ADA_TARGET = 0.30)
  events <- ev(amt = 0)

  if (regimen == "Untreated (natural history)") {
    events <- ev(amt = 0)
  } else if (regimen == "Laronidase ERT 0.58 mg/kg IV weekly") {
    events <- laro_full
  } else if (regimen == "HSCT alone (early transplant)") {
    pars$HSCT_FLAG <- 1; pars$AGE_AT_TRANSPLANT_MO <- 9
  } else if (regimen == "HSCT alone (delayed transplant)") {
    pars$HSCT_FLAG <- 1; pars$AGE_AT_TRANSPLANT_MO <- 30
  } else if (regimen == "ERT bridging (~15 wk) then HSCT") {
    pars$HSCT_FLAG <- 1; events <- laro_bridge
  } else if (regimen == "ERT, high anti-drug antibody titer") {
    events <- laro_full; pars$ADA_TARGET <- 0.85
  } else if (regimen == "ERT, poor adherence (60% of infusions)") {
    events <- laro_full; pars$ADHERENCE_ERT <- 0.6
  } else if (regimen == "ERT + genistein (investigational SRT)") {
    events <- c(laro_full, gen_ev)
  } else if (regimen == "ERT + AAV9 CNS gene therapy (investigational)") {
    events <- laro_full; pars$AAV9_FLAG <- 1
  } else if (regimen == "HSCT + long-term low-dose ERT (residual disease)") {
    events <- laro_full; pars$HSCT_FLAG <- 1
  }

  mod <- do.call(param, c(list(.x = mps1_mod), pars))
  mod %>% ev(events) %>% mrgsim(end = tmax, delta = 24) %>% as_tibble()
}

server <- function(input, output, session) {

  sim <- eventReactive(input$run, {
    run_regimen(input$regimen, input$years, input$wt, input$phenotype, input$age_transplant)
  }, ignoreNULL = FALSE)

  cmp <- eventReactive(input$run, {
    bind_rows(lapply(REGIMENS, function(r) {
      run_regimen(r, input$years, input$wt, input$phenotype, input$age_transplant) %>% mutate(regimen = r)
    }))
  }, ignoreNULL = FALSE)

  output$pheno_label <- renderText(ifelse(input$phenotype == "1", "Severe (Hurler)", "Attenuated (Hurler-Scheie/Scheie)"))
  output$reg_label <- renderText(input$regimen)

  output$schematic <- renderUI({
    tags$div(
      tags$p("IDUA (4p16.3) loss-of-function → alpha-L-iduronidase deficiency → lysosomal ",
             "dermatan/heparan sulfate accumulation → multi-organ storage: dysostosis multiplex, ",
             "valvular/myocardial disease, upper-airway obstruction, hepatosplenomegaly, corneal ",
             "clouding, and (severe Hurler) progressive neurodegeneration behind an intact blood-brain barrier."),
      tags$p("Laronidase (IV ERT) restores enzyme activity in well-perfused visceral tissue via M6PR-mediated ",
             "uptake but cannot cross the BBB. HSCT achieves donor cross-correction of visceral tissue AND ",
             "(uniquely) CNS microglial replacement — the only validated route to halt neurodegeneration, ",
             "provided transplantation occurs early. Investigational: lentiviral HSC gene therapy, AAV9 ",
             "CNS-directed gene therapy, and oral substrate-reduction therapy (genistein)."),
      tags$img(src = "mps1_qsp_model.png", style = "max-width:100%;")
    )
  })

  output$plt_pk_laro <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~LARO_CP, name = "Plasma Cp (ng/mL equiv)") %>%
      add_lines(y = ~LARO_TISSUE, name = "Tissue-retained enzyme (mg)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "value"))
  })

  output$plt_pk_other <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~ADA_INHIB, name = "ADA-mediated inhibition (fraction)") %>%
      add_lines(y = ~GEN_SYNTH_RED, name = "Genistein GAG-synthesis reduction (fraction)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "fraction"))
  })

  output$plt_enzaccess <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~ENZ_ACCESS_SYS, name = "Systemic (visceral)") %>%
      add_lines(y = ~ENZ_ACCESS_CNS, name = "CNS") %>%
      add_lines(y = ~ENZ_ACCESS_CART, name = "Cartilage/bone") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "enzyme-access index (0-1)"))
  })

  output$plt_gag <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~GAG_SYS, name = "Systemic GAG") %>%
      add_lines(y = ~GAG_CNS, name = "CNS GAG") %>%
      add_lines(y = ~GAG_CART, name = "Cartilage GAG") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "normalized burden (1.0 = untreated)"))
  })

  output$ugag_end <- renderText(sprintf("%.0f µg/mg Cr", tail(sim()$UGAG, 1)))
  output$dq_end <- renderText(sprintf("%.0f", tail(sim()$DQ, 1)))
  output$livspleen_end <- renderText(sprintf("%.2fx normal", tail(sim()$LIVSPLEEN, 1)))
  output$surv_end <- renderText(sprintf("%.1f%%", 100*tail(sim()$SURVIVAL, 1)))

  output$plt_clinical <- renderPlotly({
    d <- sim()
    p1 <- plot_ly(d, x = ~time/24, y = ~UGAG, type = "scatter", mode = "lines", name = "uGAG")
    p2 <- plot_ly(d, x = ~time/24, y = ~LIVSPLEEN, type = "scatter", mode = "lines", name = "Liver/spleen index")
    p3 <- plot_ly(d, x = ~time/24) %>% add_lines(y = ~FVC, name = "FVC %pred") %>% add_lines(y = ~AHI, name = "AHI")
    p4 <- plot_ly(d, x = ~time/24) %>% add_lines(y = ~JOINTROM, name = "Joint ROM") %>% add_lines(y = ~DQ, name = "DQ")
    subplot(p1, p2, p3, p4, nrows = 4, shareX = TRUE) %>% layout(xaxis = list(title = "Day"))
  })

  output$plt_compare <- renderPlotly({
    d <- cmp()
    plot_ly(d, x = ~time/24, y = ~UGAG, color = ~regimen, type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Urinary GAG (µg/mg creatinine)"))
  })

  output$tbl_compare <- renderDT({
    d <- cmp() %>% group_by(regimen) %>% slice_tail(n = 1) %>%
      transmute(Regimen = regimen,
                `uGAG (µg/mg Cr)` = round(UGAG, 0),
                `Liver/spleen index` = round(LIVSPLEEN, 2),
                `FVC %pred` = round(FVC, 1),
                `AHI (events/h)` = round(AHI, 1),
                `Joint ROM (deg)` = round(JOINTROM, 0),
                `DQ` = round(DQ, 0),
                `Height Z` = round(HEIGHTZ, 2),
                `Survival` = sprintf("%.1f%%", 100*SURVIVAL))
    datatable(d, options = list(pageLength = 10))
  })

  output$plt_valve_cornea <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~VALVE, name = "Valve/LV-mass index") %>%
      add_lines(y = ~CORNEA, name = "Corneal clouding score") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "value"))
  })

  output$plt_height <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24, y = ~HEIGHTZ, type = "scatter", mode = "lines", name = "Height Z-score") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Z-score"))
  })

  output$plt_hazard <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24, y = ~HAZARD, type = "scatter", mode = "lines", name = "Cumulative hazard") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "cumulative hazard index"))
  })

  output$plt_survival <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24, y = ~SURVIVAL, type = "scatter", mode = "lines", name = "Survival probability") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "probability", range = c(0,1)))
  })

  output$refs <- renderUI({
    tags$div(
      tags$p(tags$b("Pivotal trials & guidelines:")),
      tags$ul(
        tags$li("Wraith JE, et al. J Pediatr 2004 (PMID 15126990) — pivotal placebo-controlled laronidase RCT."),
        tags$li("Clarke LA, et al. Pediatrics 2009 (PMID 19117887) — laronidase long-term open-label extension."),
        tags$li("Kakkis ED, et al. N Engl J Med 2001 (PMID 11172140) — first-in-human laronidase phase 1/2."),
        tags$li("Peters C, et al. Blood 1998 (PMID 9516162) — HSCT outcome registry."),
        tags$li("Aldenhoven M, et al. Blood 2015 (PMID 25624320) — long-term HSCT outcome, transplant-age effect."),
        tags$li("Boelens JJ, et al. Blood 2013 (PMID 23493783) — multicenter HSCT donor-source outcome study."),
        tags$li("Gentner B, et al. N Engl J Med 2021 (PMID 34788506) — lentiviral HSC gene therapy (OTL-203).")
      ),
      tags$p("See mps1_references.md for the full bibliography.")
    )
  })
}

shinyApp(ui, server)
