# =====================================================================
# X-Linked Hypophosphatemia (XLH) — QSP Shiny Dashboard
# Author : Claude Code Routine (2026-07-01)
# Deps   : shiny, bslib, mrgsolve, dplyr, tidyr, ggplot2, plotly, DT
# Run    : shiny::runApp("xlh_shiny_app.R")
# =====================================================================

library(shiny)
library(bslib)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

source("xlh_mrgsolve_model.R", local = TRUE)   # brings `xlh_mod`

ui <- page_navbar(
  title = "X-Linked Hypophosphatemia QSP Explorer",
  theme = bs_theme(bootswatch = "flatly", primary = "#0277BD"),
  sidebar = sidebar(
    width = 340, open = "always",
    h5("Patient profile"),
    selectInput("poptype", "Population", c("Pediatric (5 yr)", "Adult")),
    sliderInput("wt", "Body weight (kg)", 10, 100, 18, step = 1),
    sliderInput("phos0", "Baseline serum phosphate (mg/dL)", 1.2, 3.0, 2.2, step = 0.1),
    sliderInput("rss0", "Baseline Rickets Severity Score (ped)", 0, 10, 5.5, step = 0.5),
    hr(),
    h5("Therapy"),
    selectInput("regimen", "Regimen",
                c("Untreated (natural history)",
                  "Conventional therapy (oral phosphate + calcitriol)",
                  "Burosumab 0.8 mg/kg SC Q2W (pediatric approved)",
                  "Burosumab 2.0 mg/kg SC Q2W (max pediatric)",
                  "Burosumab 1.0 mg/kg SC Q4W (adult approved)",
                  "Switch: conventional → burosumab",
                  "Conventional therapy, poor GI adherence (60%)",
                  "Burosumab, supratherapeutic overcorrection",
                  "Conventional therapy, long-term (tertiary HPT risk)",
                  "Burosumab adult, 60% adherence")),
    sliderInput("years", "Simulation horizon (yr)", 1, 4, 1, step = 1),
    actionButton("run", "Run simulation", class = "btn-primary w-100")
  ),
  nav_panel("① Patient & Overview",
    layout_columns(
      col_widths = c(4,4,4),
      value_box(title = "Baseline serum Pi", value = textOutput("pi0"), theme = "secondary"),
      value_box(title = "Untreated TmP/GFR", value = "1.8 mg/dL", theme = "secondary"),
      value_box(title = "Regimen", value = textOutput("reg_label"), theme = "primary")
    ),
    card(card_header("Mechanistic map (15-cluster QSP schematic)"), htmlOutput("schematic"))
  ),
  nav_panel("② Drug PK",
    layout_columns(col_widths = c(6,6),
      card(card_header("Burosumab concentration-time"), plotlyOutput("plt_pk_buro")),
      card(card_header("Oral phosphate / calcitriol exposure"), plotlyOutput("plt_pk_conv"))
    )
  ),
  nav_panel("③ Pathway PD (FGF23 / NPT2 / TmP-GFR)",
    layout_columns(col_widths = c(6,6),
      card(card_header("FGF23 neutralization & NPT2a/c rescue"), plotlyOutput("plt_fgf23")),
      card(card_header("TmP/GFR & serum phosphate"), plotlyOutput("plt_tmpgfr"))
    )
  ),
  nav_panel("④ Clinical endpoints",
    layout_columns(col_widths = c(4,4,4),
      value_box(title = "RSS @ end", value = textOutput("rss_end"), theme = "success"),
      value_box(title = "Height Z @ end (ped)", value = textOutput("hz_end"), theme = "success"),
      value_box(title = "6MWT Δ (adult, m)", value = textOutput("sixmwt_end"), theme = "info")
    ),
    card(card_header("RSS, growth (ped) & 6MWT/WOMAC (adult) trajectories"), plotlyOutput("plt_clinical"))
  ),
  nav_panel("⑤ Scenario comparison",
    card(card_header("Compare all pre-configured regimens (serum phosphate)"), plotlyOutput("plt_compare")),
    card(card_header("Table of endpoints at simulation end"), DTOutput("tbl_compare"))
  ),
  nav_panel("⑥ Biomarkers",
    layout_columns(col_widths = c(6,6),
      card(card_header("1,25(OH)2D · PTH"), plotlyOutput("plt_vitd_pth")),
      card(card_header("Bone-specific ALP (BSAP)"), plotlyOutput("plt_bsap"))
    )
  ),
  nav_panel("⑦ Safety",
    layout_columns(col_widths = c(6,6),
      card(card_header("Urine Ca/Cr ratio"), plotlyOutput("plt_ucacr")),
      card(card_header("Nephrocalcinosis risk index"), plotlyOutput("plt_nephrocalc"))
    )
  ),
  nav_panel("⑧ References",
    card(card_header("Key references"), htmlOutput("refs"))
  )
)

wt_mgkg <- function(wt, dose_mgkg) dose_mgkg * wt

run_regimen <- function(regimen, years, wt, phos0, rss0) {
  tmax <- years * 8760
  base <- xlh_mod %>% param(PHOS0 = phos0, RSS0 = rss0)

  events <- switch(regimen,
    "Untreated (natural history)" = ev(amt = 0),
    "Conventional therapy (oral phosphate + calcitriol)" =
      c(ev(cmt = "PHOSORAL_GUT", amt = wt_mgkg(wt, 10), ii = 6, addl = floor(tmax/6) - 1, time = 0),
        ev(cmt = "CALC_GUT", amt = wt_mgkg(wt, 0.030), ii = 24, addl = floor(tmax/24) - 1, time = 0)),
    "Burosumab 0.8 mg/kg SC Q2W (pediatric approved)" =
      ev(cmt = "BURO_DEPOT", amt = wt_mgkg(wt, 0.8), ii = 336, addl = floor(tmax/336) - 1, time = 0),
    "Burosumab 2.0 mg/kg SC Q2W (max pediatric)" =
      ev(cmt = "BURO_DEPOT", amt = wt_mgkg(wt, 2.0), ii = 336, addl = floor(tmax/336) - 1, time = 0),
    "Burosumab 1.0 mg/kg SC Q4W (adult approved)" =
      ev(cmt = "BURO_DEPOT", amt = wt_mgkg(wt, 1.0), ii = 672, addl = floor(tmax/672) - 1, time = 0),
    "Switch: conventional → burosumab" =
      ev(cmt = "BURO_DEPOT", amt = wt_mgkg(wt, 0.8), ii = 336, addl = floor(tmax/336) - 1, time = 0),
    "Conventional therapy, poor GI adherence (60%)" =
      ev(cmt = "PHOSORAL_GUT", amt = wt_mgkg(wt, 10), ii = 6, addl = floor(tmax/6) - 1, time = 0),
    "Burosumab, supratherapeutic overcorrection" =
      ev(cmt = "BURO_DEPOT", amt = wt_mgkg(wt, 3.5), ii = 336, addl = floor(tmax/336) - 1, time = 0),
    "Conventional therapy, long-term (tertiary HPT risk)" =
      ev(cmt = "CALC_GUT", amt = wt_mgkg(wt, 0.030), ii = 24, addl = floor(tmax/24) - 1, time = 0),
    "Burosumab adult, 60% adherence" =
      ev(cmt = "BURO_DEPOT", amt = wt_mgkg(wt, 1.0), ii = 672, addl = floor(tmax/672) - 1, time = 0)
  )

  mod <- base
  if (regimen == "Conventional therapy, poor GI adherence (60%)") mod <- mod %>% param(ADHERENCE_CONV = 0.6)
  if (regimen == "Burosumab adult, 60% adherence")                mod <- mod %>% param(ADHERENCE_BURO = 0.6)

  mod %>% ev(events) %>% mrgsim(end = tmax, delta = 24) %>% as_tibble()
}

server <- function(input, output, session) {

  sim <- eventReactive(input$run, {
    run_regimen(input$regimen, input$years, input$wt, input$phos0, input$rss0)
  }, ignoreNULL = FALSE)

  cmp <- eventReactive(input$run, {
    regs <- c("Untreated (natural history)",
              "Conventional therapy (oral phosphate + calcitriol)",
              "Burosumab 0.8 mg/kg SC Q2W (pediatric approved)",
              "Burosumab 2.0 mg/kg SC Q2W (max pediatric)",
              "Burosumab 1.0 mg/kg SC Q4W (adult approved)",
              "Conventional therapy, poor GI adherence (60%)",
              "Burosumab, supratherapeutic overcorrection",
              "Conventional therapy, long-term (tertiary HPT risk)",
              "Burosumab adult, 60% adherence")
    bind_rows(lapply(regs, function(r) {
      run_regimen(r, input$years, input$wt, input$phos0, input$rss0) %>% mutate(regimen = r)
    }))
  }, ignoreNULL = FALSE)

  output$pi0 <- renderText(sprintf("%.1f mg/dL", input$phos0))
  output$reg_label <- renderText(input$regimen)

  output$schematic <- renderUI({
    tags$div(
      tags$p("PHEX loss-of-function → osteocyte FGF23 overproduction → FGF23-FGFR1c/",
             "αKlotho renal signaling → NPT2a/NPT2c internalization → renal phosphate ",
             "wasting (TmP/GFR↓) + suppressed 1α-hydroxylase (calcitriol↓) → chronic ",
             "hypophosphatemia → rickets (pediatric) / osteomalacia (adult), enthesopathy, myopathy."),
      tags$p("Burosumab neutralizes circulating FGF23 directly, restoring NPT2a/c-mediated ",
             "phosphate reabsorption AND de-repressing endogenous calcitriol synthesis. ",
             "Conventional therapy (oral phosphate + calcitriol) bypasses the FGF23 axis via ",
             "direct substrate replacement, producing transient serum Pi spikes and hypercalciuria."),
      tags$img(src = "xlh_qsp_model.png", style = "max-width:100%;")
    )
  })

  output$plt_pk_buro <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24, y = ~BURO_CP, type = "scatter", mode = "lines", name = "Burosumab Cp") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Conc. (ng/mL equiv)"))
  })

  output$plt_pk_conv <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~PHOSORAL_SIG, name = "Oral phosphate exposure signal") %>%
      add_lines(y = ~CALC_CENT, name = "Oral calcitriol Cp (ng/mL equiv)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "a.u."))
  })

  output$plt_fgf23 <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~FGF23_NEUT, name = "FGF23 neutralization (%)") %>%
      add_lines(y = ~NPT2, name = "NPT2a/c activity (0-1)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "fraction"))
  })

  output$plt_tmpgfr <- renderPlotly({
    d <- sim()
    p1 <- plot_ly(d, x = ~time/24, y = ~TMPGFR, type = "scatter", mode = "lines", name = "TmP/GFR")
    p2 <- plot_ly(d, x = ~time/24, y = ~PHOS, type = "scatter", mode = "lines", name = "Serum Pi")
    subplot(p1, p2, nrows = 2, shareX = TRUE) %>% layout(xaxis = list(title = "Day"))
  })

  output$rss_end <- renderText(sprintf("%.2f", tail(sim()$RSS, 1)))
  output$hz_end  <- renderText(sprintf("%.2f", tail(sim()$HEIGHTZ_XLH, 1)))
  output$sixmwt_end <- renderText(sprintf("%+.0f", tail(sim()$SIXMWT, 1) - sim()$SIXMWT[1]))

  output$plt_clinical <- renderPlotly({
    d <- sim()
    p1 <- plot_ly(d, x = ~time/24, y = ~RSS, type = "scatter", mode = "lines", name = "RSS")
    p2 <- plot_ly(d, x = ~time/24, y = ~HEIGHTZ_XLH, type = "scatter", mode = "lines", name = "Height Z")
    p3 <- plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~SIXMWT, name = "6MWT (m)") %>% add_lines(y = ~WOMAC, name = "WOMAC")
    subplot(p1, p2, p3, nrows = 3, shareX = TRUE) %>% layout(xaxis = list(title = "Day"))
  })

  output$plt_compare <- renderPlotly({
    d <- cmp()
    plot_ly(d, x = ~time/24, y = ~PHOS, color = ~regimen, type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Serum phosphate (mg/dL)"))
  })

  output$tbl_compare <- renderDT({
    d <- cmp() %>% group_by(regimen) %>% slice_tail(n = 1) %>%
      transmute(Regimen = regimen,
                `Serum Pi (mg/dL)` = round(PHOS, 2),
                `TmP/GFR` = round(TMPGFR, 2),
                `RSS` = round(RSS, 2),
                `Height Z (ped)` = round(HEIGHTZ_XLH, 2),
                `6MWT (m)` = round(SIXMWT, 0),
                `Nephrocalcinosis risk` = round(NEPHROCALC, 3))
    datatable(d, options = list(pageLength = 10))
  })

  output$plt_vitd_pth <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~CALCITRIOL, name = "1,25(OH)2D (pg/mL)") %>%
      add_lines(y = ~PTH, name = "PTH (pg/mL)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "value"))
  })

  output$plt_bsap <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24, y = ~BSAP, type = "scatter", mode = "lines", name = "BSAP (µg/L)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "µg/L"))
  })

  output$plt_ucacr <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24, y = ~UCACR, type = "scatter", mode = "lines", name = "Urine Ca/Cr") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "mg/mg"))
  })

  output$plt_nephrocalc <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24, y = ~NEPHROCALC, type = "scatter", mode = "lines", name = "Nephrocalcinosis risk") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "risk index (0-1)"))
  })

  output$refs <- renderUI({
    tags$div(
      tags$p(tags$b("Pivotal trials & guidelines:")),
      tags$ul(
        tags$li("Carpenter TO, et al. NEJM 2018 (PMID 29791829) — Burosumab phase 2 pediatric dose-finding."),
        tags$li("Imel EA, et al. Lancet 2019 (PMID 31104833) — CL303 phase 3 RCT, burosumab vs conventional therapy."),
        tags$li("Whyte MP, et al. Lancet Diabetes Endocrinol 2019 (PMID 31104830) — phase 2, ages 1-4 yr."),
        tags$li("Insogna KL, et al. JBMR 2018 (PMID 29947083) — AXLES1 adult phase 3 placebo-controlled RCT."),
        tags$li("Haffner D, et al. Nat Rev Nephrol 2019 (PMID 31068690) — international consensus guideline.")
      ),
      tags$p("See xlh_references.md for the full bibliography.")
    )
  })
}

shinyApp(ui, server)
