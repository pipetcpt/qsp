# =====================================================================
# Neurofibromatosis Type 1 (NF1) — QSP Shiny Dashboard
# Author : Claude Code Routine (2026-07-01)
# Deps   : shiny, bslib, mrgsolve, dplyr, tidyr, ggplot2, plotly, DT
# Run    : shiny::runApp("nf1_shiny_app.R")
# =====================================================================

library(shiny)
library(bslib)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

source("nf1_mrgsolve_model.R", local = TRUE)   # brings `nf1_mod`

ui <- page_navbar(
  title = "Neurofibromatosis Type 1 (NF1-PN) QSP Explorer",
  theme = bs_theme(bootswatch = "flatly", primary = "#6A1B9A"),
  sidebar = sidebar(
    width = 360, open = "always",
    h5("Patient profile"),
    selectInput("poptype", "Population", c("Pediatric (10 yr)", "Adult")),
    sliderInput("bsa", "Body surface area (m2)", 0.6, 2.2, 1.10, step = 0.05),
    sliderInput("pn_vol0", "Baseline target-PN volume (mL index)", 40, 250, 100, step = 5),
    checkboxInput("puberty", "Puberty / pregnancy growth acceleration", FALSE),
    checkboxInput("opg", "Concurrent optic pathway glioma (OPG)", FALSE),
    hr(),
    h5("Therapy"),
    selectInput("regimen", "Regimen",
                c("Untreated (natural history)",
                  "Selumetinib 25 mg/m2 BID (SPRINT, pediatric)",
                  "Selumetinib 20 mg/m2 BID (dose-reduced, AE)",
                  "Mirdametinib 2 mg/m2 BID, 3wk-on/1wk-off (ReNeu, pediatric)",
                  "Mirdametinib 2 mg/m2 BID, 3wk-on/1wk-off (ReNeu, adult)",
                  "Selumetinib, drug holiday then rechallenge",
                  "Selumetinib, poor adherence (60%)",
                  "Trametinib off-label approximation",
                  "Mirdametinib adult, long-term (5 yr)")),
    sliderInput("years", "Simulation horizon (yr)", 1, 5, 2, step = 1),
    actionButton("run", "Run simulation", class = "btn-primary w-100")
  ),
  nav_panel("① Patient & Overview",
    layout_columns(
      col_widths = c(4,4,4),
      value_box(title = "Baseline PN volume", value = textOutput("pnvol0"), theme = "secondary"),
      value_box(title = "Baseline tumor pain (NRS-11)", value = "5.5", theme = "secondary"),
      value_box(title = "Regimen", value = textOutput("reg_label"), theme = "primary")
    ),
    card(card_header("Neurofibromin-RAS-MAPK mechanistic map"), uiOutput("schematic"))
  ),
  nav_panel("② Drug PK",
    card(card_header("Selumetinib / mirdametinib plasma concentration"), plotlyOutput("plt_pk", height = 380))
  ),
  nav_panel("③ Pathway PD (MEK inhibition / pERK / resistance)",
    card(card_header("Fractional MEK inhibition, pERK suppression, adaptive resistance"), plotlyOutput("plt_pd", height = 380))
  ),
  nav_panel("④ Clinical endpoints (PN / OPG / cNF / pain / QoL)",
    layout_columns(
      col_widths = c(4,4,4),
      value_box(title = "PN volume response", value = textOutput("pn_resp_end"), theme = "success"),
      value_box(title = "Pain (NRS-11, end)", value = textOutput("pain_end"), theme = "success"),
      value_box(title = "HRQoL (end)", value = textOutput("qol_end"), theme = "success")
    ),
    plotlyOutput("plt_clinical", height = 480)
  ),
  nav_panel("⑤ Scenario comparison",
    plotlyOutput("plt_compare", height = 380),
    DTOutput("tbl_compare")
  ),
  nav_panel("⑥ Biomarkers (pERK/resistance/growth)",
    plotlyOutput("plt_biomarkers", height = 420)
  ),
  nav_panel("⑦ Safety (LVEF / dermatologic / CPK)",
    plotlyOutput("plt_safety", height = 420)
  ),
  nav_panel("⑧ References",
    card(card_header("Key references"), htmlOutput("refs"))
  )
)

bsa_mg <- function(bsa, dose_mgm2) dose_mgm2 * bsa

run_regimen <- function(regimen, years, bsa, pn_vol0, puberty, opg) {
  tmax <- years * 8760
  base <- nf1_mod %>% param(PN_VOL0 = pn_vol0,
                             PUBERTY_MULT = if (puberty) 2.2 else 1.0,
                             OPG_ONSET = if (opg) 1 else 0)

  events <- switch(regimen,
    "Untreated (natural history)" = ev(amt = 0),
    "Selumetinib 25 mg/m2 BID (SPRINT, pediatric)" =
      ev(cmt = "SEL_GUT", amt = bsa_mg(bsa, 25), ii = 12, addl = floor(tmax/12) - 1, time = 0),
    "Selumetinib 20 mg/m2 BID (dose-reduced, AE)" =
      ev(cmt = "SEL_GUT", amt = bsa_mg(bsa, 20), ii = 12, addl = floor(tmax/12) - 1, time = 0),
    "Mirdametinib 2 mg/m2 BID, 3wk-on/1wk-off (ReNeu, pediatric)" =
      ev(cmt = "MIR_GUT", amt = bsa_mg(bsa, 2), ii = 12, addl = floor(tmax/12) - 1, time = 0),
    "Mirdametinib 2 mg/m2 BID, 3wk-on/1wk-off (ReNeu, adult)" =
      ev(cmt = "MIR_GUT", amt = bsa_mg(bsa, 2), ii = 12, addl = floor(tmax/12) - 1, time = 0),
    "Selumetinib, drug holiday then rechallenge" =
      c(ev(cmt = "SEL_GUT", amt = bsa_mg(bsa, 25), ii = 12, addl = floor((tmax*0.4)/12) - 1, time = 0),
        ev(cmt = "SEL_GUT", amt = bsa_mg(bsa, 25), ii = 12, addl = floor((tmax*0.4)/12) - 1, time = tmax*0.6)),
    "Selumetinib, poor adherence (60%)" =
      ev(cmt = "SEL_GUT", amt = bsa_mg(bsa, 25), ii = 12, addl = floor(tmax/12) - 1, time = 0),
    "Trametinib off-label approximation" =
      ev(cmt = "MIR_GUT", amt = bsa_mg(bsa, 2), ii = 24, addl = floor(tmax/24) - 1, time = 0),
    "Mirdametinib adult, long-term (5 yr)" =
      ev(cmt = "MIR_GUT", amt = bsa_mg(bsa, 2), ii = 12, addl = floor(tmax/12) - 1, time = 0)
  )

  mod <- base
  if (regimen == "Selumetinib, poor adherence (60%)") mod <- mod %>% param(ADHERENCE_SEL = 0.6)
  if (regimen == "Trametinib off-label approximation") mod <- mod %>% param(EC50_MIR = 15, CL_MIR = 8, V_MIR = 110)

  mod %>% ev(events) %>% mrgsim(end = tmax, delta = 24) %>% as_tibble()
}

server <- function(input, output, session) {

  sim <- eventReactive(input$run, {
    run_regimen(input$regimen, input$years, input$bsa, input$pn_vol0, input$puberty, input$opg)
  }, ignoreNULL = FALSE)

  cmp <- eventReactive(input$run, {
    regs <- c("Untreated (natural history)",
              "Selumetinib 25 mg/m2 BID (SPRINT, pediatric)",
              "Selumetinib 20 mg/m2 BID (dose-reduced, AE)",
              "Mirdametinib 2 mg/m2 BID, 3wk-on/1wk-off (ReNeu, pediatric)",
              "Selumetinib, poor adherence (60%)",
              "Trametinib off-label approximation",
              "Mirdametinib adult, long-term (5 yr)")
    bind_rows(lapply(regs, function(r) {
      run_regimen(r, input$years, input$bsa, input$pn_vol0, input$puberty, input$opg) %>% mutate(regimen = r)
    }))
  }, ignoreNULL = FALSE)

  output$pnvol0 <- renderText(sprintf("%.0f mL (index)", input$pn_vol0))
  output$reg_label <- renderText(input$regimen)

  output$schematic <- renderUI({
    tags$div(
      tags$p("NF1 (17q11.2) biallelic loss → neurofibromin RAS-GAP deficiency → ",
             "constitutive RAS-GTP → RAF-MEK-ERK hyperactivation in NF1-null Schwann cells → ",
             "plexiform neurofibroma (PN) growth, cutaneous neurofibroma accumulation, optic ",
             "pathway glioma, and MPNST transformation risk."),
      tags$p("Selumetinib / mirdametinib (oral MEK1/2 inhibitors) suppress pERK, driving ",
             "tumor-growth-inhibition of PN/OPG volume (REiNS ≥20% response threshold), with ",
             "adaptive RTK-feedback resistance on chronic dosing and rebound on discontinuation."),
      tags$img(src = "nf1_qsp_model.png", style = "max-width:100%;")
    )
  })

  output$plt_pk <- renderPlotly({
    d <- sim()
    plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~SEL_CP, name = "Selumetinib Cp (ng/mL)") %>%
      add_lines(y = ~MIR_CP, name = "Mirdametinib Cp (ng/mL)") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Conc. (ng/mL)"))
  })

  output$plt_pd <- renderPlotly({
    d <- sim()
    p1 <- plot_ly(d, x = ~time/24, y = ~TOTAL_INHIB, type = "scatter", mode = "lines", name = "MEK inhibition fraction")
    p2 <- plot_ly(d, x = ~time/24, y = ~PERK_SUPPRESSION, type = "scatter", mode = "lines", name = "pERK suppression")
    p3 <- plot_ly(d, x = ~time/24, y = ~RESIST, type = "scatter", mode = "lines", name = "Adaptive resistance")
    subplot(p1, p2, p3, nrows = 3, shareX = TRUE) %>% layout(xaxis = list(title = "Day"))
  })

  output$pn_resp_end <- renderText(sprintf("%+.1f%%", tail(sim()$PN_RESPONSE_PCT, 1)))
  output$pain_end <- renderText(sprintf("%.1f", tail(sim()$PAIN, 1)))
  output$qol_end <- renderText(sprintf("%.0f", tail(sim()$QOL, 1)))

  output$plt_clinical <- renderPlotly({
    d <- sim()
    p1 <- plot_ly(d, x = ~time/24, y = ~PN_TOTAL, type = "scatter", mode = "lines", name = "PN volume")
    p2 <- plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~OPG_VOL, name = "OPG volume") %>% add_lines(y = ~CNF_BURDEN, name = "cNF burden")
    p3 <- plot_ly(d, x = ~time/24) %>%
      add_lines(y = ~PAIN, name = "Pain (NRS-11)") %>% add_lines(y = ~QOL, name = "HRQoL") %>%
      add_lines(y = ~VISION, name = "Visual acuity index")
    subplot(p1, p2, p3, nrows = 3, shareX = TRUE) %>% layout(xaxis = list(title = "Day"))
  })

  output$plt_compare <- renderPlotly({
    d <- cmp()
    plot_ly(d, x = ~time/24, y = ~PN_RESPONSE_PCT, color = ~regimen, type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "PN volume change (%)"))
  })

  output$tbl_compare <- renderDT({
    d <- cmp() %>% group_by(regimen) %>% slice_tail(n = 1) %>%
      transmute(Regimen = regimen,
                `PN response (%)` = round(PN_RESPONSE_PCT, 1),
                `Pain (NRS-11)` = round(PAIN, 1),
                `HRQoL` = round(QOL, 0),
                `LVEF (%)` = round(LVEF, 1),
                `Dermatologic AE` = round(DERM_AE, 0),
                `CPK (xULN)` = round(CPK_AE, 2))
    datatable(d, options = list(pageLength = 10))
  })

  output$plt_biomarkers <- renderPlotly({
    d <- sim()
    p1 <- plot_ly(d, x = ~time/24, y = ~RESIST, type = "scatter", mode = "lines", name = "Adaptive resistance")
    p2 <- plot_ly(d, x = ~time/24, y = ~GROWTHZ, type = "scatter", mode = "lines", name = "Pediatric growth Z-score")
    subplot(p1, p2, nrows = 2, shareX = TRUE) %>% layout(xaxis = list(title = "Day"))
  })

  output$plt_safety <- renderPlotly({
    d <- sim()
    p1 <- plot_ly(d, x = ~time/24, y = ~LVEF, type = "scatter", mode = "lines", name = "LVEF (%)")
    p2 <- plot_ly(d, x = ~time/24, y = ~DERM_AE, type = "scatter", mode = "lines", name = "Dermatologic AE composite")
    p3 <- plot_ly(d, x = ~time/24, y = ~CPK_AE, type = "scatter", mode = "lines", name = "CPK (xULN)")
    subplot(p1, p2, p3, nrows = 3, shareX = TRUE) %>% layout(xaxis = list(title = "Day"))
  })

  output$refs <- renderUI({
    tags$div(
      tags$p(tags$b("Pivotal trials & guidelines:")),
      tags$ul(
        tags$li("Gross AM, et al. NEJM 2020 (PMID 32187457) — SPRINT phase 2, selumetinib 25 mg/m2 BID pediatric inoperable PN."),
        tags$li("Dombi E, et al. NEJM 2016 (PMID 28029918) — phase 1 dose-finding, 20-30 mg/m2 BID."),
        tags$li("Moertel CL, et al. J Clin Oncol 2025 (PMID 39514826) — ReNeu phase 2b, mirdametinib adults+children."),
        tags$li("Gutmann DH, et al. Nat Rev Dis Primers 2017 (PMID 28230061) — NF1 disease primer."),
        tags$li("Legius E, et al. Genet Med 2021 (PMID 34012067) — revised NF1/Legius diagnostic criteria."),
        tags$li("Dombi E, et al. Neurology 2013 (PMID 24249804) — REiNS volumetric response criteria.")
      ),
      tags$p("See nf1_references.md for the full bibliography.")
    )
  })
}

shinyApp(ui, server)
