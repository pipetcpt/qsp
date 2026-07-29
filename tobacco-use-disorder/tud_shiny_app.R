## ===========================================================================
##  Tobacco Use Disorder (Nicotine Dependence) — QSP Shiny dashboard
##  ---------------------------------------------------------------------------
##  File      : tud_shiny_app.R
##  Companion : tud_mrgsolve_model.R  (must be in the same directory)
##  Run       : shiny::runApp("tud_shiny_app.R")
##
##  8 tabs:
##    1. Patient profile      — phenotype, CYP2A6 genotype, dependence severity
##    2. Nicotine PK          — nicotine / cotinine / 3HC, NMR, delivery route
##    3. Receptor pharmacology— a4b2* occupancy, desensitization, upregulation
##    4. Reward & withdrawal  — dopamine tone, set-point, deficit, MNWS, QSU
##    5. Clinical endpoints   — CAR wk 9-12, 7-day PP, CO verification
##    6. Scenario comparison  — all 10 therapy arms side by side
##    7. Biomarkers           — NMR treatment matching, CO, cotinine, PET beta2*
##    8. Safety & tolerability— nausea, sleep, weight, FEV1, persistence
## ===========================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

source("tud_mrgsolve_model.R", local = TRUE)   # provides mod, scenarios, builders

PAL <- c("#B71C1C", "#1565C0", "#2E7D32", "#EF6C00", "#6A1B9A",
         "#00838F", "#AD1457", "#4E342E", "#37474F", "#827717")

theme_tud <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          legend.position  = "bottom",
          plot.title       = element_text(face = "bold", size = 13))
}

quit_line <- function() {
  list(geom_vline(xintercept = 0, linetype = 2, colour = "grey35"),
       annotate("text", x = 0, y = Inf, label = " TQD", hjust = 0, vjust = 1.4,
                size = 3, colour = "grey35"))
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Tobacco Use Disorder (Nicotine Dependence) — QSP Dashboard"),
  tags$p(tags$em(paste(
    "Educational QSP model: nicotine PK + CYP2A6 pharmacogenetics,",
    "a4b2* nicotinic receptor desensitization/upregulation, mesolimbic dopamine,",
    "allostatic withdrawal, craving, and lapse hazard. NOT for clinical use."))),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Therapy"),
      selectInput("scenario", "Scenario", choices = names(scenarios),
                  selected = names(scenarios)[5]),
      checkboxInput("counsel", "Behavioural support / counselling", TRUE),

      hr(), h4("Smoker phenotype"),
      sliderInput("cpd0",    "Cigarettes per day",       5, 60, 20, 1),
      sliderInput("niccig",  "Nicotine absorbed per cigarette (mg)",
                  0.02, 2.0, 1.1, 0.01),
      sliderInput("yrsmoke", "Years of smoking",          1, 60, 25, 1),

      hr(), h4("Pharmacogenetics"),
      sliderInput("f2a6", "CYP2A6 activity (F2A6)", 0.15, 2.2, 1.0, 0.05),
      helpText(HTML(paste0(
        "0.35 = slow (*4/*4, *9), 1.0 = normal, 1.7 = fast.<br>",
        "This single scalar sets NMR, nicotine clearance and therefore how much ",
        "relief a <b>fixed</b> NRT dose delivers."))),
      sliderInput("f2b6", "CYP2B6 activity (bupropion)", 0.3, 1.6, 1.0, 0.05),
      sliderInput("a5",   "CHRNA5 rs16969968 risk alleles", 0, 2, 0, 1),
      checkboxInput("psy", "Psychiatric history", FALSE),

      hr(), h4("Receptor pharmacology"),
      sliderInput("emaxv", "Varenicline intrinsic activity (EMAXV)",
                  0, 1, 0.45, 0.05),
      helpText("0 = pure antagonist, 1 = full agonist. Sweep this on tab 3."),
      sliderInput("emaxr", "Max upregulation (EMAXR)", 0, 1.5, 0.75, 0.05),
      sliderInput("phid",  "Signal loss per unit desensitization (PHID)",
                  0.3, 0.95, 0.78, 0.01),

      hr(), h4("Renal function"),
      sliderInput("rfv", "Varenicline CL multiplier (RFV)", 0.4, 1.2, 1.0, 0.05),
      helpText("CrCl < 30 mL/min: label reduces dose; set RFV ~ 0.5."),

      hr(),
      sliderInput("horizon", "Follow-up (weeks after TQD)", 12, 52, 24, 4),
      actionButton("go", "Run simulation", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1. Patient profile",
          br(),
          fluidRow(column(6, h4("Derived baseline"), tableOutput("profile")),
                   column(6, h4("Steady-state smoker phenotype"),
                          tableOutput("ss_tbl"))),
          hr(),
          h4("Occupancy calibration (Brody 2006)"),
          plotOutput("occ_curve", height = 300),
          helpText(paste("Effective in vivo Kd = 5.4 nM reproduces the reported",
                         "EC50 of 0.87 ng/mL for 50% beta2* occupancy."))
        ),

        tabPanel("2. Nicotine PK",
          br(), plotOutput("pk_plot", height = 430),
          hr(), plotOutput("pk_zoom", height = 280),
          helpText(paste("Lower panel: 72 h around the quit date. Note the",
                         "overnight trough and morning resensitization that",
                         "produce time-to-first-cigarette craving."))
        ),

        tabPanel("3. Receptor pharmacology",
          br(), plotOutput("rec_plot", height = 430),
          hr(), h4("Why PARTIAL agonism — sweep of varenicline intrinsic activity"),
          plotOutput("emaxv_sweep", height = 300),
          tableOutput("emaxv_tbl")
        ),

        tabPanel("4. Reward & withdrawal",
          br(), plotOutput("pd_plot", height = 460),
          hr(),
          helpText(HTML(paste0(
            "<b>DA</b> = dopaminergic supply (fast, minutes). ",
            "<b>SETP</b> = hedonic set-point, a ~14-day moving average of DA. ",
            "<b>DEF = SETP - DA</b> is the withdrawal driver. In a smoker at ",
            "steady state DEF is ~0; withdrawal is the transient between two ",
            "matched states, and its duration is the set-point time constant, ",
            "not a fitted curve.")))
        ),

        tabPanel("5. Clinical endpoints",
          br(), plotOutput("car_plot", height = 380),
          hr(), h4("Endpoint summary (current scenario)"),
          tableOutput("ep_tbl"),
          hr(),
          h4("Validation against trial continuous-abstinence rates"),
          actionButton("go_val", "Run validation (7 arms x 5 CYP2A6 strata)",
                       class = "btn-warning"),
          br(), br(), tableOutput("val_tbl"),
          helpText(paste("Population-averaged over the CYP2A6 distribution.",
                         "Anchors: EAGLES (Anthenelli 2016), ORCA-2 (Rigotti",
                         "2023), Cochrane NRT reviews. Reference RMSE ~2.8",
                         "percentage points.")),
          helpText(paste("CAR = continuous abstinence probability. Simulated",
                         "conditional on remaining abstinent, so PK/PD curves",
                         "describe the abstinent trajectory."))
        ),

        tabPanel("6. Scenario comparison",
          br(),
          checkboxGroupInput("cmp", "Arms to compare", choices = names(scenarios),
                             selected = names(scenarios)[c(2, 3, 4, 5, 7, 8, 9)],
                             inline = FALSE),
          actionButton("go_cmp", "Run comparison", class = "btn-primary"),
          br(), br(),
          plotOutput("cmp_plot", height = 520),
          hr(), tableOutput("cmp_tbl"),
          hr(),
          h4("Varenicline + patch: competitive displacement at a4b2*"),
          actionButton("go_disp", "Run displacement analysis",
                       class = "btn-warning"),
          br(), br(), tableOutput("disp_tbl"),
          helpText(paste("The model predicts NON-additivity: adding a patch",
                         "raises total occupancy only slightly but nearly",
                         "halves varenicline's own occupancy, and it is that",
                         "term which carries the partial-agonist blockade.",
                         "This disagrees with Koegelenberg 2014 and agrees",
                         "with later neutral trials - a falsifiable",
                         "prediction, not a fitted result."))
        ),

        tabPanel("7. Biomarkers",
          br(),
          h4("NMR treatment matching (Lerman 2015 Lancet Respir Med)"),
          actionButton("go_nmr", "Run NMR experiment", class = "btn-primary"),
          br(), br(),
          plotOutput("nmr_plot", height = 340), tableOutput("nmr_tbl"),
          hr(), h4("Exposure & verification biomarkers over time"),
          plotOutput("bm_plot", height = 320)
        ),

        tabPanel("8. Safety & tolerability",
          br(), plotOutput("safety_plot", height = 460),
          hr(), tableOutput("safety_tbl"),
          helpText(paste("Weight gain approaches ~5 kg in sustained abstainers",
                         "(Aubin 2012); FEV1 slope halves after cessation",
                         "(Fletcher-Peto / Lung Health Study)."))
        ),

        tabPanel("About",
          br(),
          h4("Model structure"),
          tags$ul(
            tags$li("31 ODE compartments (nicotine 2-cmt + biophase, cotinine, 3HC, patch and buccal depots, varenicline, cytisinicline, bupropion + hydroxybupropion, receptor desensitization and upregulation, dopamine, set-point, allostatic load, MNWS, QSU, habit, abstinence survival, COHb, weight, FEV1, nausea, tolerance, sleep, persistence)"),
            tags$li("Competitive a4b2* occupancy with intrinsic-activity weighting; non-competitive block by hydroxybupropion"),
            tags$li("Route-resolved phasic reinforcement (cigarette 1.0, oral NRT 0.45, patch 0.05), scaled by reinforcement-episode frequency (cigarettes/day)"),
            tags$li("CYP2A6-driven nicotine titration: intake scales with clearance^0.5, so fast metabolizers smoke more and a fixed NRT dose under-replaces them"),
            tags$li("Lapse hazard as a proportional-hazards function of craving, withdrawal and habit")
          ),
          h4("Key calibration anchors"),
          tags$ul(
            tags$li("Nicotine CL 72 L/h, t1/2 ~1.8 h; 20 cig/day -> 11-17 ng/mL"),
            tags$li("beta2* occupancy EC50 0.87 ng/mL (Brody 2006)"),
            tags$li("Varenicline 1 mg BID -> ~90% occupancy (Lotfipour 2012)"),
            tags$li("beta2* upregulation normalizes over 3-4 weeks (Cosgrove 2009)"),
            tags$li("CAR wk 9-12 EAGLES: varenicline 33.5 / bupropion 22.6 / patch 23.4 / placebo 12.5 %")
          ),
          h4("Disclaimer"),
          tags$p(tags$strong(paste(
            "Educational and research use only. Not validated, not for clinical",
            "decision-making, prescribing, or regulatory submission.")))
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# SERVER
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  user_par <- reactive({
    list(CPD0 = input$cpd0, NICCIG = input$niccig, YRSMOKE = input$yrsmoke,
         F2A6 = input$f2a6, F2B6 = input$f2b6, A5RISK = input$a5,
         PSYHX = as.numeric(input$psy), COUNSEL = as.numeric(input$counsel),
         EMAXV = input$emaxv, EMAXR = input$emaxr, PHID = input$phid,
         RFV = input$rfv)
  })

  sim <- eventReactive(input$go, ignoreNULL = FALSE, {
    s <- scenarios[[input$scenario]]
    s$par <- modifyList(s$par, user_par())
    tq <- s$par$TQD
    end <- if (is.finite(tq) && tq < 1e5) tq + input$horizon*7*24 else
             TQD_H + input$horizon*7*24
    run_scenario(input$scenario, s, end = end, delta = 2) |>
      mutate(day = (time - min(tq, TQD_H))/24)
  })

  # ---- tab 1 -------------------------------------------------------------
  output$profile <- renderTable({
    p <- user_par()
    cl <- 18 + 54*p$F2A6
    data.frame(
      Quantity = c("Cigarettes / day", "Nicotine intake (mg/day)",
                   "Years smoked", "Habit strength (0-1)",
                   "CYP2A6 activity", "Nicotine clearance (L/h)",
                   "Nicotine t1/2 (h)", "Predicted NMR (3HC/cotinine)"),
      Value = c(p$CPD0, round(p$CPD0*p$NICCIG, 1), p$YRSMOKE,
                round(p$YRSMOKE/(p$YRSMOKE + 8), 3), p$F2A6, round(cl, 1),
                round(0.693*182/cl, 2), round(0.45*p$F2A6*3.0/3.5, 3))
    )
  }, digits = 3)

  output$ss_tbl <- renderTable({
    d <- sim(); pre <- d[d$day < 0, ]
    if (!nrow(pre)) return(NULL)
    last <- tail(pre, 24)
    data.frame(
      Variable = c("Plasma nicotine (ng/mL, 24 h mean)", "Cotinine (ng/mL)",
                   "NMR", "a4b2* occupancy", "Desensitized fraction",
                   "Receptor pool RUP", "Dopamine tone DA", "Set-point SETP",
                   "MNWS", "QSU", "Exhaled CO (ppm)"),
      Value = c(mean(last$NIC_NGML), mean(last$COT_NGML), mean(last$NMR_),
                mean(last$OCC_TOT), mean(last$DES), mean(last$RUP),
                mean(last$DA), mean(last$SETP), mean(last$WD), mean(last$QSU),
                mean(last$CO_PPM))
    )
  }, digits = 3)

  output$occ_curve <- renderPlot({
    d <- occupancy_curve()
    ggplot(d, aes(nic_ngml, 100*occupancy)) +
      geom_line(linewidth = 1, colour = PAL[1]) +
      geom_hline(yintercept = 50, linetype = 3) +
      geom_vline(xintercept = 0.87, linetype = 3) +
      annotate("point", x = 0.87, y = 50, size = 3) +
      annotate("text", x = 0.95, y = 42,
               label = "Brody 2006: 0.87 ng/mL = 50%", hjust = 0, size = 3.4) +
      scale_x_log10() +
      labs(x = "Plasma nicotine (ng/mL, log scale)",
           y = expression(beta*"2* occupancy (%)")) +
      theme_tud()
  })

  # ---- tab 2 -------------------------------------------------------------
  output$pk_plot <- renderPlot({
    d <- sim() |>
      select(day, NIC_NGML, COT_NGML, HC_NGML, NMR_, VAR_NGML, CYT_NGML,
             BUP_NGML, OHB_NGML) |>
      pivot_longer(-day) |> filter(!is.na(value))
    ggplot(d, aes(day, value)) +
      geom_line(colour = PAL[2], linewidth = 0.5) +
      facet_wrap(~name, scales = "free_y", ncol = 2) + quit_line() +
      labs(x = "Days from TQD", y = NULL, title = "Nicotine and drug PK") +
      theme_tud()
  })

  output$pk_zoom <- renderPlot({
    d <- sim() |> filter(day >= -2, day <= 3)
    ggplot(d, aes(day, NIC_NGML)) +
      geom_line(colour = PAL[1], linewidth = 0.8) + quit_line() +
      labs(x = "Days from TQD", y = "Plasma nicotine (ng/mL)",
           title = "Diurnal nicotine profile around the quit date") +
      theme_tud()
  })

  # ---- tab 3 -------------------------------------------------------------
  output$rec_plot <- renderPlot({
    d <- sim() |>
      select(day, OCC_NIC, OCC_VAR, OCC_CYT, OCC_TOT, DES, RUP, ACT_) |>
      pivot_longer(-day)
    ggplot(d, aes(day, value)) +
      geom_line(colour = PAL[5], linewidth = 0.55) +
      facet_wrap(~name, scales = "free_y", ncol = 2) + quit_line() +
      labs(x = "Days from TQD", y = NULL,
           title = "a4b2* occupancy, desensitization, receptor pool, activation") +
      theme_tud()
  })

  sweep_res <- eventReactive(input$go, ignoreNULL = FALSE, {
    partial_agonism_sweep()
  })

  output$emaxv_sweep <- renderPlot({
    d <- sweep_res() |>
      select(EMAXV, mnws_peak, lapse_reward, CAR_wk12) |>
      pivot_longer(-EMAXV)
    ggplot(d, aes(EMAXV, value, colour = name)) +
      geom_line(linewidth = 0.9) + geom_point(size = 1.6) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = PAL) +
      labs(x = "Varenicline intrinsic activity (0 = antagonist, 1 = full agonist)",
           y = NULL) +
      theme_tud() + theme(legend.position = "none")
  })
  output$emaxv_tbl <- renderTable(sweep_res(), digits = 3)

  val_res  <- eventReactive(input$go_val,  { validate_against_trials() })
  output$val_tbl <- renderTable(val_res(), digits = 2)

  disp_res <- eventReactive(input$go_disp, { combination_displacement() })
  output$disp_tbl <- renderTable(disp_res(), digits = 3)

  # ---- tab 4 -------------------------------------------------------------
  output$pd_plot <- renderPlot({
    d <- sim() |> select(day, DA, SETP, DEF_, ALLO, WD, QSU, HABIT) |>
      pivot_longer(-day)
    ggplot(d, aes(day, value)) +
      geom_line(colour = PAL[3], linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 2) + quit_line() +
      labs(x = "Days from TQD", y = NULL,
           title = "Dopamine supply vs hedonic set-point; withdrawal and craving") +
      theme_tud()
  })

  # ---- tab 5 -------------------------------------------------------------
  output$car_plot <- renderPlot({
    d <- sim()
    ggplot(d, aes(day, CAR_PCT)) +
      geom_line(linewidth = 1, colour = PAL[1]) + quit_line() +
      geom_vline(xintercept = 84, linetype = 3) +
      annotate("text", x = 84, y = 100, label = " week 12", hjust = 0, size = 3.2) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = "Days from TQD", y = "Continuous abstinence (%)",
           title = "Abstinence survival curve") +
      theme_tud()
  })

  output$ep_tbl <- renderTable({
    summarise_endpoints(sim()) |> select(-scenario) |> t() |> as.data.frame() |>
      setNames("Value") |> tibble::rownames_to_column("Endpoint")
  }, digits = 3)

  # ---- tab 6 -------------------------------------------------------------
  cmp <- eventReactive(input$go_cmp, {
    req(length(input$cmp) > 0)
    end <- TQD_H + input$horizon*7*24
    bind_rows(lapply(input$cmp, function(n) {
      s <- scenarios[[n]]; s$par <- modifyList(s$par, user_par())
      run_scenario(n, s, end = end, delta = 4)
    })) |> mutate(day = (time - TQD_H)/24)
  })

  output$cmp_plot <- renderPlot({
    d <- cmp() |> select(day, scenario, NIC_NGML, OCC_TOT, RUP, DA, WD, QSU,
                         CAR_PCT, WT) |> pivot_longer(-c(day, scenario))
    ggplot(d, aes(day, value, colour = scenario)) +
      geom_line(linewidth = 0.55) +
      facet_wrap(~name, scales = "free_y", ncol = 2) + quit_line() +
      scale_colour_manual(values = PAL) +
      labs(x = "Days from TQD", y = NULL) +
      theme_tud() +
      guides(colour = guide_legend(ncol = 2, title = NULL)) +
      theme(legend.text = element_text(size = 7))
  })
  output$cmp_tbl <- renderTable(summarise_endpoints(cmp()), digits = 2)

  # ---- tab 7 -------------------------------------------------------------
  nmr_res <- eventReactive(input$go_nmr, { nmr_experiment() })

  output$nmr_plot <- renderPlot({
    ggplot(nmr_res(), aes(NMR, CAR_wk12, colour = arm)) +
      geom_line(linewidth = 1) + geom_point(size = 2.6) +
      scale_colour_manual(values = PAL[c(2, 1)]) +
      labs(x = "Nicotine metabolite ratio (3HC / cotinine)",
           y = "CAR week 9-12 (%)",
           title = "One PK parameter (CYP2A6) creates the NMR x treatment interaction",
           subtitle = "Patch relief degrades as metabolism accelerates; varenicline is renally cleared and unaffected") +
      theme_tud()
  })
  output$nmr_tbl <- renderTable(nmr_res(), digits = 3)

  output$bm_plot <- renderPlot({
    d <- sim() |> select(day, COT_NGML, CO_PPM, NMR_, RUP) |> pivot_longer(-day)
    ggplot(d, aes(day, value)) +
      geom_line(colour = PAL[6], linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y") + quit_line() +
      labs(x = "Days from TQD", y = NULL) + theme_tud()
  })

  # ---- tab 8 -------------------------------------------------------------
  output$safety_plot <- renderPlot({
    d <- sim() |> select(day, NAUS, TOL, SLP, WT, FEV_PCT, ADH_PCT) |>
      pivot_longer(-day)
    ggplot(d, aes(day, value)) +
      geom_line(colour = PAL[4], linewidth = 0.6) +
      facet_wrap(~name, scales = "free_y", ncol = 2) + quit_line() +
      labs(x = "Days from TQD", y = NULL,
           title = "Tolerability, persistence and organ-system outcomes") +
      theme_tud()
  })

  output$safety_tbl <- renderTable({
    d <- sim()
    data.frame(
      Endpoint = c("Peak nausea index", "Peak sleep disturbance",
                   "Weight change at end (kg)", "FEV1 at end (% baseline)",
                   "Medication persistence at end (%)"),
      Value = c(max(d$NAUS), max(d$SLP), tail(d$WT, 1), tail(d$FEV_PCT, 1),
                tail(d$ADH_PCT, 1))
    )
  }, digits = 3)
}

shinyApp(ui, server)
