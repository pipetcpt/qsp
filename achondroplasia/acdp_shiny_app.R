# =====================================================================
# Achondroplasia (ACH) — QSP Shiny Dashboard
# Author : Claude Code Routine (2026-07-01)
# Deps   : shiny, bslib, mrgsolve, dplyr, tidyr, ggplot2, plotly, DT
# Run    : shiny::runApp("acdp_shiny_app.R")
# =====================================================================

library(shiny)
library(bslib)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

source("acdp_mrgsolve_model.R", local = TRUE)   # brings `acdp_mod`

ui <- page_navbar(
  title = "Achondroplasia QSP Explorer",
  theme = bs_theme(bootswatch = "flatly", primary = "#6A1B9A"),
  sidebar = sidebar(
    width = 340, open = "always",
    h5("Patient profile"),
    sliderInput("age", "Age at start (yr)", 0.5, 14, 5, step = 0.5),
    sliderInput("wt", "Body weight (kg)", 5, 45, 15, step = 1),
    sliderInput("height0", "Baseline height (cm)", 55, 130, 85, step = 1),
    sliderInput("heightz0", "Baseline height Z-score", -7, -2, -5, step = 0.5),
    sliderInput("fm0", "Baseline foramen magnum area (mm^2)", 150, 400, 280, step = 10),
    checkboxInput("severe_fms", "Severe foramen magnum stenosis (surgical candidate)", FALSE),
    hr(),
    h5("Therapy"),
    selectInput("regimen", "Regimen",
                c("Untreated (natural history)",
                  "Vosoritide 15 µg/kg SC QD (approved dose)",
                  "Vosoritide 2.5 µg/kg SC QD",
                  "Vosoritide 7.5 µg/kg SC QD",
                  "Vosoritide 30 µg/kg SC QD",
                  "TransCon CNP (navepegritide) SC QW",
                  "Infigratinib PO QD (FGFR1-3 TKI)",
                  "Growth hormone (off-label)",
                  "Vosoritide + FMD surgery",
                  "Vosoritide, 60% adherence")),
    sliderInput("years", "Simulation horizon (yr)", 1, 5, 1, step = 1),
    actionButton("run", "Run simulation", class = "btn-primary w-100")
  ),
  nav_panel("① Patient & Overview",
    layout_columns(
      col_widths = c(4,4,4),
      value_box(title = "Baseline height Z", value = textOutput("bz"), theme = "secondary"),
      value_box(title = "Untreated AGV", value = "3.9 cm/yr", theme = "secondary"),
      value_box(title = "Regimen", value = textOutput("reg_label"), theme = "primary")
    ),
    card(card_header("Mechanistic map (15-cluster QSP schematic)"), htmlOutput("schematic"))
  ),
  nav_panel("② Drug PK",
    layout_columns(col_widths = c(6,6),
      card(card_header("Vosoritide / TransCon CNP (free-CNP moiety)"), plotlyOutput("plt_pk_cnp")),
      card(card_header("Infigratinib"), plotlyOutput("plt_pk_infig"))
    )
  ),
  nav_panel("③ Pathway PD (pERK / cGMP / chondrocyte rescue)",
    layout_columns(col_widths = c(6,6),
      card(card_header("pERK activity vs. cGMP counter-signal"), plotlyOutput("plt_perk")),
      card(card_header("Chondrocyte proliferation index"), plotlyOutput("plt_chondro"))
    )
  ),
  nav_panel("④ Growth endpoints",
    layout_columns(col_widths = c(4,4,4),
      value_box(title = "AGV @ end (cm/yr)", value = textOutput("agv_end"), theme = "success"),
      value_box(title = "Height Z @ end", value = textOutput("hz_end"), theme = "success"),
      value_box(title = "Cumulative height (cm)", value = textOutput("ht_end"), theme = "info")
    ),
    card(card_header("Height Z-score & cumulative height trajectory"), plotlyOutput("plt_growth"))
  ),
  nav_panel("⑤ Structural & biomarker endpoints",
    layout_columns(col_widths = c(6,6),
      card(card_header("Foramen magnum area · Spinal canal Z-score"), plotlyOutput("plt_structural")),
      card(card_header("OSA-AHI · Otitis media rate · BMI-Z"), plotlyOutput("plt_biomarkers"))
    )
  ),
  nav_panel("⑥ Scenario comparison",
    card(card_header("Compare all pre-configured regimens (Height Z-score)"), plotlyOutput("plt_compare")),
    card(card_header("Table of endpoints at simulation end"), DTOutput("tbl_compare"))
  ),
  nav_panel("⑦ Safety",
    layout_columns(col_widths = c(6,6),
      card(card_header("Hemodynamics: MAP · HR"), plotlyOutput("plt_hemo")),
      card(card_header("Serum phosphate (FGFR1 off-target)"), plotlyOutput("plt_phos"))
    )
  ),
  nav_panel("⑧ References",
    card(card_header("Key references"), htmlOutput("refs"))
  )
)

wt_ugkg <- function(wt, dose_ugkg) dose_ugkg * wt

run_regimen <- function(regimen, years, wt, height0, heightz0, fm0, severe_fms) {
  tmax <- years * 8760
  base <- acdp_mod %>% param(HEIGHT0 = height0, HEIGHTZ0 = heightz0, FMAREA0 = fm0)

  events <- switch(regimen,
    "Untreated (natural history)" = ev(amt = 0),
    "Vosoritide 15 µg/kg SC QD (approved dose)" =
      ev(cmt = "VOS_DEPOT", amt = wt_ugkg(wt, 15), ii = 24, addl = floor(tmax/24) - 1, time = 0),
    "Vosoritide 2.5 µg/kg SC QD" =
      ev(cmt = "VOS_DEPOT", amt = wt_ugkg(wt, 2.5), ii = 24, addl = floor(tmax/24) - 1, time = 0),
    "Vosoritide 7.5 µg/kg SC QD" =
      ev(cmt = "VOS_DEPOT", amt = wt_ugkg(wt, 7.5), ii = 24, addl = floor(tmax/24) - 1, time = 0),
    "Vosoritide 30 µg/kg SC QD" =
      ev(cmt = "VOS_DEPOT", amt = wt_ugkg(wt, 30), ii = 24, addl = floor(tmax/24) - 1, time = 0),
    "TransCon CNP (navepegritide) SC QW" =
      ev(cmt = "TCNP_DEPOT", amt = wt_ugkg(wt, 100), ii = 168, addl = floor(tmax/168) - 1, time = 0),
    "Infigratinib PO QD (FGFR1-3 TKI)" =
      ev(cmt = "INFIG_GUT", amt = wt_ugkg(wt, 0.5), ii = 24, addl = floor(tmax/24) - 1, time = 0),
    "Growth hormone (off-label)" = ev(amt = 0),
    "Vosoritide + FMD surgery" =
      ev(cmt = "VOS_DEPOT", amt = wt_ugkg(wt, 15), ii = 24, addl = floor(tmax/24) - 1, time = 0),
    "Vosoritide, 60% adherence" =
      ev(cmt = "VOS_DEPOT", amt = wt_ugkg(wt, 15), ii = 24, addl = floor(tmax/24) - 1, time = 0)
  )

  mod <- base
  if (regimen == "Growth hormone (off-label)") mod <- mod %>% param(GH_ON = 1)
  if (regimen == "Vosoritide, 60% adherence")  mod <- mod %>% param(ADHERENCE = 0.6)
  if (regimen == "Vosoritide + FMD surgery" && severe_fms) mod <- mod %>% param(FMAREA0 = fm0 + 60)

  mod %>% ev(events) %>% mrgsim(end = tmax, delta = 24) %>% as_tibble()
}

server <- function(input, output, session) {

  sim <- eventReactive(input$run, {
    run_regimen(input$regimen, input$years, input$wt, input$height0,
                input$heightz0, input$fm0, input$severe_fms)
  }, ignoreNULL = FALSE)

  cmp <- eventReactive(input$run, {
    regs <- c("Untreated (natural history)",
              "Vosoritide 15 µg/kg SC QD (approved dose)",
              "Vosoritide 2.5 µg/kg SC QD", "Vosoritide 7.5 µg/kg SC QD",
              "Vosoritide 30 µg/kg SC QD", "TransCon CNP (navepegritide) SC QW",
              "Infigratinib PO QD (FGFR1-3 TKI)", "Growth hormone (off-label)",
              "Vosoritide + FMD surgery", "Vosoritide, 60% adherence")
    bind_rows(lapply(regs, function(r) {
      run_regimen(r, input$years, input$wt, input$height0, input$heightz0,
                  input$fm0, input$severe_fms) %>% mutate(regimen = r)
    }))
  }, ignoreNULL = FALSE)

  output$bz <- renderText(sprintf("%.1f", input$heightz0))
  output$reg_label <- renderText(input$regimen)

  output$schematic <- renderUI({
    tags$div(
      tags$p("FGFR3 G380R → constitutive RAS-RAF-MEK-ERK hyperactivation → growth-plate ",
             "chondrocyte proliferation/differentiation suppression → rhizomelic short ",
             "stature, foramen magnum stenosis, spinal stenosis, OSA, otitis media."),
      tags$p("Vosoritide / TransCon CNP act on NPR-B → cGMP → PKGII → inhibitory RAF1 ",
             "phosphorylation, counter-regulating the MAPK cascade. Infigratinib blocks ",
             "the FGFR3 kinase domain directly (with FGFR1 off-target liability)."),
      tags$img(src = "acdp_qsp_model.png", style = "max-width:100%;")
    )
  })

  output$plt_pk_cnp <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~VOS_CP, name = "Vosoritide Cp") %>%
      add_lines(y = ~TCNP_CP, name = "TransCon-released free CNP") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Conc. (ng/mL equiv)"))
  })

  output$plt_pk_infig <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24, y = ~INFIG_CP, type = "scatter", mode = "lines",
            name = "Infigratinib Cp") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Conc. (ng/mL equiv)"))
  })

  output$plt_perk <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~PERK, name = "pERK activity") %>%
      add_lines(y = ~CGMP_SIG, name = "cGMP counter-signal") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "a.u."))
  })

  output$plt_chondro <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24, y = ~CHONDRO, type = "scatter", mode = "lines",
            name = "Chondrocyte proliferation index") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "index (0-1)"))
  })

  output$agv_end <- renderText(sprintf("%.2f", tail(sim()$AGV_CALC, 1)))
  output$hz_end  <- renderText(sprintf("%.2f", tail(sim()$HEIGHTZ, 1)))
  output$ht_end  <- renderText(sprintf("%.1f", tail(sim()$HEIGHT_CM, 1)))

  output$plt_growth <- renderPlotly({
    d <- sim()
    p1 <- plot_ly(d, x = ~time/24, y = ~HEIGHTZ, type = "scatter", mode = "lines", name = "Height Z")
    p2 <- plot_ly(d, x = ~time/24, y = ~HEIGHT_CM, type = "scatter", mode = "lines", name = "Height (cm)")
    subplot(p1, p2, nrows = 2, shareX = TRUE) %>% layout(xaxis = list(title = "Day"))
  })

  output$plt_structural <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~FMAREA, name = "Foramen magnum area (mm^2)") %>%
      add_lines(y = ~SPCANALZ*100, name = "Spinal canal Z (x100)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "value"))
  })

  output$plt_biomarkers <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~AHI, name = "OSA-AHI (events/h)") %>%
      add_lines(y = ~OTITIS, name = "Otitis media (episodes/yr)") %>%
      add_lines(y = ~BMIZ, name = "BMI-Z") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "value"))
  })

  output$plt_compare <- renderPlotly({
    d <- cmp()
    plot_ly(d, x = ~time/24, y = ~HEIGHTZ, color = ~regimen, type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Height Z-score"))
  })

  output$tbl_compare <- renderDT({
    d <- cmp() %>% group_by(regimen) %>% slice_tail(n = 1) %>%
      transmute(Regimen = regimen,
                `AGV (cm/yr)` = round(AGV_CALC, 2),
                `Height Z` = round(HEIGHTZ, 2),
                `Height (cm)` = round(HEIGHT_CM, 1),
                `Foramen magnum (mm^2)` = round(FMAREA, 0),
                `OSA-AHI` = round(AHI, 1),
                `Serum phosphate` = round(PHOS, 2))
    datatable(d, options = list(pageLength = 10))
  })

  output$plt_hemo <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~MAP_BP, name = "MAP (mmHg)") %>%
      add_lines(y = ~HR, name = "HR (bpm)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "value"))
  })

  output$plt_phos <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24, y = ~PHOS, type = "scatter", mode = "lines",
            name = "Serum phosphate (mg/dL)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "mg/dL"))
  })

  output$refs <- renderUI({
    tags$div(
      tags$p(tags$b("Pivotal trials:")),
      tags$ul(
        tags$li("Savarirayan R, et al. NEJM 2019 (PMID 31269546) — Vosoritide phase 2 dose-finding."),
        tags$li("Savarirayan R, et al. Lancet 2020 (PMID 32891212) — Vosoritide phase 3, 52-week RCT (ΔAGV +1.57 cm/yr)."),
        tags$li("Savarirayan R, et al. 2021 (PMID 34341520) — Vosoritide phase 3 open-label extension, 2-year."),
        tags$li("BridgeBio/QED — PROPEL 2 (NEJM 2025, PMID 39555818) and PROPEL 3 (NEJM 2026) phase 2/3 infigratinib; positive, NDA planned Q3-2026."),
        tags$li("Ascendis Pharma — ApproaCH phase 3 (TransCon CNP / navepegritide); FDA-approved Feb-2026 as YUVIWEL.")
      ),
      tags$p("See acdp_references.md for the full 50-item bibliography.")
    )
  })
}

shinyApp(ui, server)
