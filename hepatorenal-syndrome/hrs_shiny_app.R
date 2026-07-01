# =====================================================================
# Hepatorenal Syndrome (HRS) — QSP Shiny Dashboard
# Author : Claude Code Routine (2026-07-01)
# Deps   : shiny, bslib, mrgsolve, dplyr, ggplot2, plotly, DT
# Run    : shiny::runApp("hrs_shiny_app.R")
# =====================================================================

library(shiny)
library(bslib)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)

source("hrs_mrgsolve_model.R", local = TRUE)   # brings `hrs_mod`

ui <- page_navbar(
  title = "HRS QSP Explorer",
  theme = bs_theme(bootswatch = "flatly", primary = "#0277BD"),
  sidebar = sidebar(
    width = 340, open = "always",
    h5("Patient profile"),
    selectInput("child", "Child-Pugh class", c("B (7-9)", "C (10-15)"), "C (10-15)"),
    sliderInput("meld", "MELD-Na", min = 15, max = 40, value = 26, step = 1),
    sliderInput("base_map", "Baseline MAP (mmHg)", 55, 90, 70),
    sliderInput("base_cr", "Baseline Cr (mg/dL)", 1.0, 5.0, 2.5, step = 0.1),
    sliderInput("base_alb", "Baseline albumin (g/dL)", 1.5, 4.0, 2.8, step = 0.1),
    checkboxGroupInput("precip", "Precipitants",
                       choices = c("SBP" = "SBP", "GI bleed" = "GIB",
                                   "LVP w/o albumin" = "LVP", "NSAID exposure" = "NSAID")),
    hr(),
    h5("Therapy"),
    selectInput("regimen", "Regimen",
                c("Natural history", "Terlipressin bolus + albumin",
                  "Terlipressin CI + albumin", "Norepinephrine + albumin",
                  "Midodrine + Octreotide + albumin",
                  "Albumin only", "TIPS surrogate")),
    sliderInput("terli_dose", "Terlipressin dose per bolus (mg)", 0.5, 2.0, 1.0, step = 0.5),
    sliderInput("terli_int", "Terli dosing interval (h)", 4, 12, 6, step = 2),
    sliderInput("ne_rate", "NE infusion (mg/h)", 0, 3, 0.5, step = 0.1),
    sliderInput("alb_dose", "Daily albumin (g)", 0, 100, 40, step = 5),
    sliderInput("days", "Simulation horizon (d)", 3, 21, 14, step = 1),
    actionButton("run", "Run simulation", class = "btn-primary w-100")
  ),
  nav_panel("① Overview",
    layout_columns(
      col_widths = c(6,6),
      card(card_header("MAP · Renal blood flow"), plotlyOutput("plt_hemo")),
      card(card_header("GFR · Serum creatinine"), plotlyOutput("plt_renal"))
    ),
    card(card_header("Model schematic — 15-cluster QSP map"), htmlOutput("schematic"))
  ),
  nav_panel("② Drug PK",
    layout_columns(col_widths=c(6,6),
      card(card_header("Terlipressin / Lysyl-VP"), plotlyOutput("plt_terli")),
      card(card_header("NE · Midodrine · Octreotide"), plotlyOutput("plt_pk_other"))
    ),
    card(card_header("Albumin concentration"), plotlyOutput("plt_alb"))
  ),
  nav_panel("③ Neurohormonal",
    layout_columns(col_widths=c(6,6),
      card(card_header("RAAS: renin & aldosterone"), plotlyOutput("plt_raas")),
      card(card_header("SNS & AVP"), plotlyOutput("plt_sns"))
    )
  ),
  nav_panel("④ Renal & urinary",
    layout_columns(col_widths=c(6,6),
      card(card_header("Urine Na · output"), plotlyOutput("plt_urine")),
      card(card_header("Serum Na"), plotlyOutput("plt_na"))
    )
  ),
  nav_panel("⑤ Clinical endpoints",
    layout_columns(col_widths=c(4,4,4),
      value_box(title = "HRS reversal (48h)", value = textOutput("reversal"), theme = "primary"),
      value_box(title = "30-day survival", value = textOutput("s30"), theme = "success"),
      value_box(title = "90-day survival", value = textOutput("s90"), theme = "info")
    ),
    card(card_header("MELD trajectory · Survival"), plotlyOutput("plt_endpoint"))
  ),
  nav_panel("⑥ Scenario comparison",
    card(card_header("Compare all pre-configured regimens"), plotlyOutput("plt_compare")),
    card(card_header("Table of steady-state (day 7)"), DTOutput("tbl_compare"))
  ),
  nav_panel("⑦ Safety",
    layout_columns(col_widths=c(6,6),
      card(card_header("Ischemic AUC (Terlipressin)"), plotlyOutput("plt_isch")),
      card(card_header("MAP overshoot / arrhythmia risk"), plotlyOutput("plt_safety"))
    )
  ),
  nav_panel("⑧ References",
    card(card_header("Key references"), htmlOutput("refs"))
  )
)

run_scenario <- function(regimen, days, terli_dose, terli_int, ne_rate, alb_dose, precip,
                         base_map, base_cr, base_alb) {
  mod <- hrs_mod %>% param(BASE_MAP = base_map, BASE_CR = base_cr, BASE_ALB = base_alb,
                           FLAG_SBP = as.numeric("SBP" %in% precip),
                           FLAG_GIB = as.numeric("GIB" %in% precip),
                           FLAG_LVP = as.numeric("LVP" %in% precip),
                           FLAG_NSAID = as.numeric("NSAID" %in% precip))
  tmax <- days*24
  n_addl_terli <- floor(tmax/terli_int) - 1
  n_addl_alb <- days - 1
  ev_alb <- ev(cmt="ALB_C", amt=alb_dose, ii=24, addl=n_addl_alb, time=0)
  events <- switch(regimen,
    "Natural history" = ev(amt=0),
    "Terlipressin bolus + albumin" = c(
      ev(cmt="TERLI_C", amt=terli_dose*1e6, ii=terli_int, addl=n_addl_terli, time=0), ev_alb),
    "Terlipressin CI + albumin" = c(
      ev(cmt="TERLI_C", amt=terli_dose*4*1e6, rate=terli_dose*4*1e6/24,
         ii=24, addl=days-1, time=0), ev_alb),
    "Norepinephrine + albumin" = c(
      ev(cmt="NE_C", amt=ne_rate*1e6, rate=ne_rate*1e6, ii=1, addl=tmax-1, time=0), ev_alb),
    "Midodrine + Octreotide + albumin" = c(
      ev(cmt="GUT_MID", amt=12.5, ii=8, addl=floor(tmax/8)-1, time=0),
      ev(cmt="DEP_OCT", amt=0.2, ii=8, addl=floor(tmax/8)-1, time=0), ev_alb),
    "Albumin only" = ev_alb,
    "TIPS surrogate" = ev(amt=0)
  )
  if (regimen == "TIPS surrogate") {
    mod <- mod %>% param(BASE_RBF = 800, BASE_SVR = 900)
  }
  mod %>% ev(events) %>% mrgsim(end=tmax, delta=0.5) %>% as_tibble()
}

server <- function(input, output, session) {
  sim <- eventReactive(input$run, {
    run_scenario(input$regimen, input$days, input$terli_dose, input$terli_int,
                 input$ne_rate, input$alb_dose, input$precip,
                 input$base_map, input$base_cr, input$base_alb)
  }, ignoreNULL = FALSE)

  cmp <- eventReactive(input$run, {
    regs <- c("Natural history", "Terlipressin bolus + albumin",
              "Terlipressin CI + albumin", "Norepinephrine + albumin",
              "Midodrine + Octreotide + albumin", "Albumin only", "TIPS surrogate")
    bind_rows(lapply(regs, function(r) {
      run_scenario(r, input$days, input$terli_dose, input$terli_int,
                   input$ne_rate, input$alb_dose, input$precip,
                   input$base_map, input$base_cr, input$base_alb) %>%
        mutate(regimen = r)
    }))
  }, ignoreNULL = FALSE)

  ln <- function(df, y, ttl) plot_ly(df, x = ~time/24, y = df[[y]], type="scatter",
                                     mode="lines") %>%
                              layout(title = ttl, xaxis=list(title="Day"),
                                     yaxis=list(title=y))

  output$plt_hemo <- renderPlotly({
    d <- sim()
    plot_ly(d, x=~time/24) %>%
      add_lines(y=~MAP_c, name="MAP (mmHg)") %>%
      add_lines(y=~RBF_c/10, name="RBF (mL/min /10)") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="value"))
  })
  output$plt_renal <- renderPlotly({
    d <- sim()
    plot_ly(d, x=~time/24) %>%
      add_lines(y=~GFR_c, name="GFR (mL/min)") %>%
      add_lines(y=~SCR_c*30, name="sCr × 30 (mg/dL scaled)") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="value"))
  })
  output$plt_terli <- renderPlotly({
    d <- sim()
    plot_ly(d, x=~time/24) %>%
      add_lines(y=~TERLI_C/1e3, name="Terli central (μg)") %>%
      add_lines(y=~LVP_C/1e3, name="Lysyl-VP central (μg)") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="μg or ng/mL"))
  })
  output$plt_pk_other <- renderPlotly({
    d <- sim()
    plot_ly(d, x=~time/24) %>%
      add_lines(y=~CP_NE, name="NE (ng/mL)") %>%
      add_lines(y=~CP_DES, name="Desglymid (ng/mL)") %>%
      add_lines(y=~CP_OCT, name="Octreotide (ng/mL)")
  })
  output$plt_alb <- renderPlotly({ ln(sim(), "ALB_C", "Albumin cpt (g)") })
  output$plt_raas <- renderPlotly({
    d <- sim()
    plot_ly(d, x=~time/24) %>%
      add_lines(y=~PRA_c, name="Renin activity") %>%
      add_lines(y=~ALDO_c/100, name="Aldo (÷100)")
  })
  output$plt_sns <- renderPlotly({
    d <- sim()
    plot_ly(d, x=~time/24) %>%
      add_lines(y=~NE_endo/100, name="Endog NE (÷100)") %>%
      add_lines(y=~AVP_c, name="AVP")
  })
  output$plt_urine <- renderPlotly({
    d <- sim()
    plot_ly(d, x=~time/24) %>%
      add_lines(y=~UNa_c, name="Urine Na (mEq/L)") %>%
      add_lines(y=~UOUT_c/50, name="Urine out ÷50 (mL/d)")
  })
  output$plt_na <- renderPlotly({ ln(sim(), "Na_c", "Serum Na (mEq/L)") })
  output$plt_isch <- renderPlotly({ ln(sim(), "IschAUC", "Ischemic AUC (ng·h/mL)") })
  output$plt_safety <- renderPlotly({
    d <- sim(); plot_ly(d, x=~time/24) %>%
      add_lines(y=~MAP_c, name="MAP") %>%
      add_lines(y=~ifelse(MAP_c>110,1,0)*100, name="Overshoot flag ×100")
  })
  output$plt_endpoint <- renderPlotly({
    d <- sim(); plot_ly(d, x=~time/24) %>%
      add_lines(y=~MELD_pred, name="MELD estimate") %>%
      add_lines(y=~S30*100, name="30-day surv %") %>%
      add_lines(y=~S90*100, name="90-day surv %")
  })
  output$plt_compare <- renderPlotly({
    d <- cmp() %>% group_by(regimen, time) %>% summarise(SCR=mean(SCR_c),.groups="drop")
    plot_ly(d, x=~time/24, y=~SCR, color=~regimen, type="scatter", mode="lines") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="sCr (mg/dL)"))
  })
  output$tbl_compare <- renderDT({
    cmp() %>% group_by(regimen) %>%
      slice_max(time, n=1) %>% ungroup() %>%
      transmute(regimen, day=time/24, MAP=round(MAP_c,1),
                GFR=round(GFR_c,1), sCr=round(SCR_c,2),
                MELD=round(MELD_pred,1), S30=round(S30,2),
                Reversal=HRS_reversal, Isch=IschRisk) %>%
      datatable(options=list(pageLength=10, dom="t"), rownames=FALSE)
  })
  output$reversal <- renderText({
    d <- sim(); if(any(d$HRS_reversal==1)) "Achieved" else "Not achieved"
  })
  output$s30 <- renderText({ sprintf("%.1f%%", tail(sim()$S30,1)*100) })
  output$s90 <- renderText({ sprintf("%.1f%%", tail(sim()$S90,1)*100) })
  output$schematic <- renderText({
    "<div style='text-align:center'><img src='hrs_qsp_model.svg' style='max-width:100%;height:auto'/></div>"
  })
  output$refs <- renderText({
    "<ul>
    <li>Wong F et al. NEJM 2021 (CONFIRM) — terlipressin + albumin for HRS-AKI.</li>
    <li>Angeli P et al. J Hepatol 2019 — ICA-2019 revised HRS diagnostic criteria.</li>
    <li>Sanyal AJ et al. Gastroenterology 2008 (OT-0401) — terli reversal.</li>
    <li>China L et al. NEJM 2021 (ATTIRE) — targeted albumin therapy.</li>
    <li>Caraceni P et al. Lancet 2018 (ANSWER) — long-term albumin in ascites.</li>
    <li>Cavallin M et al. Hepatology 2015 — terli CI vs bolus.</li>
    <li>Boyer TD et al. Gastroenterology 2016 (REVERSE) — terli in HRS-1.</li>
    </ul>"
  })
}

shinyApp(ui, server)
