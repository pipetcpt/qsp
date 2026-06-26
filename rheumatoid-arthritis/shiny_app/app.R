## ============================================================
## Rheumatoid Arthritis QSP — Shiny App
## Interactive PK/PD Simulator
## ============================================================
## Dependencies: shiny, shinydashboard, mrgsolve, dplyr,
##               ggplot2, plotly, DT, tidyr, viridis
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)

## ── Compile mrgsolve model (same as ra_model.R) ──────────────
source("../ra_model.R", local = TRUE)
base_mod <- mod_ra  # compiled model from ra_model.R

## ══════════════════════════════════════════════════════════════
## UI
## ══════════════════════════════════════════════════════════════
ui <- dashboardPage(
  skin = "red",

  ## ── Header ─────────────────────────────────────────────────
  dashboardHeader(
    title = span(
      icon("pills"), "RA-QSP Simulator",
      style = "font-size:16px; font-weight:bold;"
    ),
    titleWidth = 260
  ),

  ## ── Sidebar ─────────────────────────────────────────────────
  dashboardSidebar(
    width = 260,

    sidebarMenu(
      id = "sidebar_menu",
      menuItem("PK / PD Simulation", tabName = "simulation",
               icon = icon("chart-line")),
      menuItem("Parameter Sensitivity", tabName = "sensitivity",
               icon = icon("sliders-h")),
      menuItem("Mechanistic Map",     tabName = "map",
               icon = icon("project-diagram")),
      menuItem("About",              tabName = "about",
               icon = icon("info-circle"))
    ),

    hr(),

    ## ── Drug Selection ─────────────────────────────────────
    h5("  Drug Regimen", style = "color:#ddd; margin-left:10px;"),

    # Tocilizumab
    checkboxInput("use_tcz", "Tocilizumab (IL-6Rα mAb)", value = TRUE),
    conditionalPanel(
      condition = "input.use_tcz",
      selectInput("tcz_route", "Route", choices = c("SC", "IV"), selected = "SC"),
      sliderInput("tcz_dose", "TCZ Dose (mg)", 80, 800, 162, step = 10),
      selectInput("tcz_interval", "Interval (days)",
                  choices = c("7 (weekly)"=7, "14 (q2w)"=14, "28 (q4w)"=28),
                  selected = 14)
    ),

    hr(style = "border-color:#444"),

    # Methotrexate
    checkboxInput("use_mtx", "Methotrexate (weekly, oral)", value = TRUE),
    conditionalPanel(
      condition = "input.use_mtx",
      sliderInput("mtx_dose", "MTX Dose (mg/week)", 5, 25, 15, step = 2.5)
    ),

    hr(style = "border-color:#444"),

    # NSAIDs / Glucocorticoids (qualitative switch)
    checkboxInput("use_nsaid", "NSAID (COX-2 inhibitor)", value = FALSE),
    conditionalPanel(
      condition = "input.use_nsaid",
      sliderInput("nsaid_cox_inh", "COX-2 Inhibition (%)", 0, 90, 60, 5)
    ),
    checkboxInput("use_gc", "Glucocorticoid (prednisone)", value = FALSE),
    conditionalPanel(
      condition = "input.use_gc",
      sliderInput("gc_nfkb_inh", "NF-κB Inhibition (%)", 0, 80, 40, 5)
    ),

    hr(style = "border-color:#444"),

    ## ── Patient Parameters ────────────────────────────────
    h5("  Patient Parameters", style = "color:#ddd; margin-left:10px;"),
    sliderInput("patient_wt", "Body weight (kg)", 40, 120, 70, step = 5),
    sliderInput("baseline_das28", "Baseline DAS28-CRP", 3.2, 8.0, 5.8, 0.1),
    sliderInput("baseline_crp",   "Baseline CRP (mg/L)", 5, 80, 20, 1),
    sliderInput("sim_weeks",       "Simulation Duration (weeks)", 12, 52, 24, 4),

    hr(style = "border-color:#444"),

    ## ── Action Button ─────────────────────────────────────
    div(
      actionButton("run_sim", "  Run Simulation",
                   icon = icon("play"),
                   style = "width:90%; background-color:#c0392b; color:white;
                            border:none; margin-left:5%; font-weight:bold;"),
      style = "text-align:center; padding:5px;"
    )
  ),

  ## ── Body ────────────────────────────────────────────────────
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background-color: #f9f9f9; }
        .small-box .icon-large { font-size: 70px; }
        .nav-tabs-custom { background: white; }
        .box { border-top-color: #c0392b; }
      "))
    ),

    tabItems(

      ## ── TAB 1: Simulation ─────────────────────────────────
      tabItem(tabName = "simulation",

        # Value boxes (week 24 summary)
        fluidRow(
          valueBoxOutput("vbox_das28",   width = 3),
          valueBoxOutput("vbox_crp",     width = 3),
          valueBoxOutput("vbox_acr50",   width = 3),
          valueBoxOutput("vbox_eular",   width = 3)
        ),

        # Main plots (row 1)
        fluidRow(
          tabBox(title = "Drug Pharmacokinetics", width = 6,
            tabPanel("Plasma Concentration",
              plotlyOutput("plot_pk", height = 300)),
            tabPanel("Receptor Blockade",
              plotlyOutput("plot_rec_occ", height = 300))
          ),
          tabBox(title = "Key Biomarkers", width = 6,
            tabPanel("CRP",
              plotlyOutput("plot_crp", height = 300)),
            tabPanel("Cytokines",
              plotlyOutput("plot_cytokines", height = 300)),
            tabPanel("RANKL",
              plotlyOutput("plot_rankl", height = 300))
          )
        ),

        # Main plots (row 2)
        fluidRow(
          box(title = "DAS28-CRP Disease Activity Score", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("plot_das28", height = 320)),
          box(title = "ACR20/50/70 Response", width = 6,
              status = "success", solidHeader = TRUE,
              plotlyOutput("plot_acr", height = 320))
        ),

        # Output table
        fluidRow(
          box(title = "Simulation Output (weekly)",
              width = 12, status = "primary", solidHeader = TRUE,
              collapsible = TRUE, collapsed = TRUE,
              DTOutput("table_sim"))
        )
      ),

      ## ── TAB 2: Sensitivity Analysis ──────────────────────
      tabItem(tabName = "sensitivity",

        fluidRow(
          box(title = "Parameter Sensitivity Controls",
              width = 4, status = "warning", solidHeader = TRUE,
              selectInput("sens_param", "Parameter to vary",
                          choices = c(
                            "CL_TCZ (clearance)"    = "CL_TCZ",
                            "V1_TCZ (volume)"       = "V1_TCZ",
                            "F_SC_TCZ (bioavail.)"  = "F_SC_TCZ",
                            "EC50_IL6_CRP"          = "EC50_IL6_CRP",
                            "ksyn_IL6"              = "ksyn_IL6",
                            "ksyn_CRP"              = "ksyn_CRP"
                          )),
              sliderInput("sens_fold", "Fold range (relative to nominal)",
                          0.1, 3.0, c(0.5, 2.0), step = 0.1),
              sliderInput("sens_n", "Number of levels", 3, 9, 5, step = 2),
              actionButton("run_sens", "Run Sensitivity",
                           icon = icon("play"),
                           style = "background-color:#e67e22; color:white;
                                    border:none; font-weight:bold;")
          ),
          box(title = "DAS28-CRP Sensitivity", width = 8,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("plot_sens_das28", height = 380))
        ),

        fluidRow(
          box(title = "CRP Sensitivity", width = 6,
              plotlyOutput("plot_sens_crp", height = 300)),
          box(title = "Peak Conc. vs. Parameter", width = 6,
              plotlyOutput("plot_sens_pk", height = 300))
        )
      ),

      ## ── TAB 3: Mechanistic Map ────────────────────────────
      tabItem(tabName = "map",
        fluidRow(
          box(title = "RA QSP Mechanistic Map", width = 12,
              status = "danger", solidHeader = TRUE,
              p(style = "color:gray;",
                "Full mechanistic network: immune activation, cytokine signaling,
                 JAK-STAT/NF-κB pathways, synovial pathology, bone/cartilage, pain,
                 PK/PD mechanisms for all approved DMARDs and biologics."),
              imageOutput("map_image", height = "700px"),
              br(),
              downloadButton("dl_svg", "Download SVG", class = "btn-info"),
              downloadButton("dl_png", "Download PNG", class = "btn-success"),
              downloadButton("dl_dot", "Download DOT", class = "btn-default")
          )
        )
      ),

      ## ── TAB 4: About ────────────────────────────────────
      tabItem(tabName = "about",
        fluidRow(
          box(title = "About This QSP Model", width = 8,
              status = "primary", solidHeader = TRUE,
              h4("Model Overview"),
              p("This Quantitative Systems Pharmacology (QSP) model simulates
                the pharmacokinetics and pharmacodynamics of disease-modifying
                antirheumatic drugs (DMARDs) in Rheumatoid Arthritis (RA)."),
              h4("Modeled Mechanisms"),
              tags$ul(
                tags$li("Tocilizumab (anti-IL-6Rα): 2-compartment PK + TMDD
                          for soluble IL-6 receptor dynamics"),
                tags$li("Methotrexate: 1-compartment PK, TNF-α inhibition via
                          adenosine/AICAR pathway (Emax model)"),
                tags$li("TNF-α dynamics: RA-elevated baseline, MTX-sensitive production"),
                tags$li("IL-6 dynamics: TNF-driven feedback, production unaffected by IL6R blockade"),
                tags$li("CRP: STAT3/IL-6-driven synthesis; markedly reduced by tocilizumab"),
                tags$li("RANKL: TNF-driven upregulation; reflects bone erosion risk"),
                tags$li("DAS28-CRP: Composite of SJC28, TJC28, CRP, Patient Global"),
                tags$li("ACR20/50/70 & EULAR response thresholds")
              ),
              h4("Key Parameters"),
              tableOutput("table_params"),
              h4("Limitations"),
              tags$ul(
                tags$li("Single-patient deterministic model; does not capture between-subject variability"),
                tags$li("IL-6 production is not directly inhibited by IL-6R blockade in this model
                          (clinically IL-6 rises transiently; modeled implicitly via TMDD)"),
                tags$li("No biologic immunogenicity (ADA) dynamics in current version"),
                tags$li("Pain/TJC28 uses simplified inflammation proxy (not nociceptor model)")
              ),
              h4("References"),
              tags$a("See full references list (references.md)",
                     href = "https://github.com/pipetcpt/qsp/blob/main/rheumatoid-arthritis/references.md",
                     target = "_blank")
          ),
          box(title = "Model Quick Reference", width = 4,
              status = "info", solidHeader = TRUE,
              h5("PK Equations (TCZ)"),
              p(withMathJax(
                "$$\\frac{dC_1}{dt} = F_{SC}k_a D_{depot} - \\frac{CL}{V_1}C_1 - \\frac{Q}{V_1}C_1 + \\frac{Q}{V_2}C_2 - k_{on}C_1 R_{free} + k_{off}RC$$"
              )),
              h5("TMDD (IL-6R)"),
              p(withMathJax(
                "$$\\frac{dR_{free}}{dt} = k_{syn,R} - k_{deg,R}R_{free} - k_{on}C_1 R_{free} + k_{off}RC$$"
              )),
              h5("CRP PD"),
              p(withMathJax(
                "$$\\frac{dCRP}{dt} = k_{syn,CRP}\\frac{IL_6^{eff}}{EC_{50}+IL_6^{eff}} - k_{out,CRP}\\cdot CRP$$"
              )),
              h5("DAS28-CRP"),
              p(withMathJax(
                "$$DAS28 = 0.56\\sqrt{TJC} + 0.28\\sqrt{SJC} + 0.36\\ln(CRP+1) + 0.014 P_{Global} + 0.96$$"
              ))
          )
        )
      )
    )
  )
)


## ══════════════════════════════════════════════════════════════
## SERVER
## ══════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  ## ── Reactive simulation ───────────────────────────────────
  sim_result <- eventReactive(input$run_sim, {
    # Update model parameters based on UI inputs
    m <- param(base_mod,
      WT          = input$patient_wt,
      DAS28_base  = input$baseline_das28,
      CRP_RA0     = input$baseline_crp
    )

    # Apply NSAID effect (reduces COX-2 → lower PGE2 → lower TJC28)
    nsaid_eff <- if (input$use_nsaid) input$nsaid_cox_inh / 100 else 0
    # Apply GC effect (NF-κB suppression → lower cytokine production)
    gc_eff    <- if (input$use_gc) input$gc_nfkb_inh / 100 else 0

    # Scale cytokine synthesis rates (qualitative NSAID/GC effect)
    m <- param(m,
      ksyn_TNF = 1.20 * (1 - 0.3 * gc_eff),
      ksyn_IL6 = 6.0  * (1 - 0.2 * gc_eff)
    )

    # Update initial conditions
    m <- init_ra(m)

    # Build events
    ev_list <- list()
    if (input$use_tcz) {
      interval <- as.numeric(input$tcz_interval)
      n_doses  <- ceiling((input$sim_weeks * 7) / interval) + 1
      tcz_times <- seq(0, by = interval, length.out = n_doses)
      tcz_times <- tcz_times[tcz_times <= input$sim_weeks * 7]
      if (input$tcz_route == "SC") {
        ev_list[[1]] <- ev(cmt = "DEPOT_TCZ", amt = input$tcz_dose, time = tcz_times)
      } else {
        ev_list[[1]] <- ev(cmt = "C1_TCZ",    amt = input$tcz_dose, time = tcz_times,
                           rate = input$tcz_dose * 24)  # 1-hr infusion
      }
    }
    if (input$use_mtx) {
      mtx_times <- seq(0, by = 7, length.out = input$sim_weeks)
      mtx_times <- mtx_times[mtx_times <= input$sim_weeks * 7]
      ev_list[[2]] <- ev(cmt = "GI_MTX", amt = input$mtx_dose, time = mtx_times)
    }

    events <- if (length(ev_list) > 0) do.call(c, ev_list) else ev(time = 0, amt = 0, cmt = 1)

    mrgsim(m, events = events,
           end = input$sim_weeks * 7, delta = 0.5,
           carry_out = "evid") %>%
      as.data.frame() %>%
      filter(evid == 0)
  }, ignoreNULL = FALSE)

  ## ── Value Boxes ──────────────────────────────────────────
  w24_summary <- reactive({
    df <- sim_result()
    end_t <- max(df$time)
    df %>% filter(abs(time - end_t) < 0.6) %>% slice(1)
  })

  output$vbox_das28 <- renderValueBox({
    d <- w24_summary()
    das <- round(d$DAS28_CRP_out, 2)
    color <- if (das < 2.6) "green" else if (das < 3.2) "olive" else if (das < 5.1) "yellow" else "red"
    category <- if (das < 2.6) "Remission" else if (das < 3.2) "Low" else if (das < 5.1) "Moderate" else "High"
    valueBox(
      das, paste("DAS28-CRP (Wk", round(input$sim_weeks), ")\n", category),
      icon = icon("stethoscope"), color = color
    )
  })

  output$vbox_crp <- renderValueBox({
    d <- w24_summary()
    crp <- round(d$CRP_out, 1)
    color <- if (crp < 5) "green" else if (crp < 20) "yellow" else "red"
    valueBox(
      paste(crp, "mg/L"), "CRP",
      icon = icon("vial"), color = color
    )
  })

  output$vbox_acr50 <- renderValueBox({
    d <- w24_summary()
    acr50 <- round(d$pct_improvement, 1)
    achieved <- if (d$ACR50 == 1) "ACR50 Achieved ✓" else if (d$ACR20 == 1) "ACR20 Achieved" else "No ACR Response"
    color <- if (d$ACR70 == 1) "green" else if (d$ACR50 == 1) "olive" else if (d$ACR20 == 1) "yellow" else "red"
    valueBox(
      paste0(acr50, "%"), paste("Improvement\n", achieved),
      icon = icon("chart-bar"), color = color
    )
  })

  output$vbox_eular <- renderValueBox({
    d <- w24_summary()
    eular_txt <- if (d$EULAR_remission == 1) "EULAR Remission ✓" else if (d$EULAR_good == 1) "EULAR Good" else "Moderate/None"
    color <- if (d$EULAR_remission == 1) "green" else if (d$EULAR_good == 1) "olive" else "yellow"
    valueBox(
      round(d$DAS28_CRP_out, 2), eular_txt,
      icon = icon("check-circle"), color = color
    )
  })

  ## ── Plot: PK ──────────────────────────────────────────────
  output$plot_pk <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time, C_TCZ_mgL)) +
      geom_line(color = "#c0392b", linewidth = 1) +
      geom_hline(yintercept = 1, linetype = 2, color = "gray50") +
      scale_x_continuous("Time (days)") +
      scale_y_continuous("TCZ Conc. (mg/L)") +
      theme_bw() +
      labs(title = NULL)
    ggplotly(p) %>% layout(hovermode = "x unified")
  })

  output$plot_rec_occ <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time, fR_blocked * 100)) +
      geom_line(color = "#8e44ad", linewidth = 1) +
      geom_hline(yintercept = 80, linetype = 2, color = "gray50") +
      scale_x_continuous("Time (days)") +
      scale_y_continuous("IL-6R Blockade (%)", limits = c(0, 100)) +
      theme_bw()
    ggplotly(p) %>% layout(hovermode = "x unified")
  })

  ## ── Plot: CRP ─────────────────────────────────────────────
  output$plot_crp <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time, CRP_out)) +
      geom_line(color = "#e74c3c", linewidth = 1) +
      geom_hline(yintercept = 5,  linetype = 2, color = "green4") +
      geom_hline(yintercept = 10, linetype = 3, color = "orange") +
      scale_x_continuous("Time (days)") +
      scale_y_continuous("CRP (mg/L)") +
      theme_bw()
    ggplotly(p) %>% layout(hovermode = "x unified")
  })

  ## ── Plot: Cytokines ───────────────────────────────────────
  output$plot_cytokines <- renderPlotly({
    df <- sim_result() %>%
      select(time, `TNF-α (ng/mL)` = TNFa_out, `IL-6 (pg/mL)` = IL6_out) %>%
      pivot_longer(-time, names_to = "Cytokine", values_to = "Conc")

    # Normalize to baseline for overlay
    p <- ggplot(df, aes(time, Conc, color = Cytokine)) +
      geom_line(linewidth = 1) +
      facet_wrap(~Cytokine, scales = "free_y") +
      scale_x_continuous("Time (days)") +
      scale_y_continuous("Concentration") +
      scale_color_manual(values = c("#e74c3c","#3498db")) +
      theme_bw() + theme(legend.position = "none")
    ggplotly(p)
  })

  ## ── Plot: RANKL ───────────────────────────────────────────
  output$plot_rankl <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time, RANKL_out)) +
      geom_line(color = "#795548", linewidth = 1) +
      scale_x_continuous("Time (days)") +
      scale_y_continuous("RANKL (pmol/L)") +
      theme_bw() +
      labs(subtitle = "Bone erosion driver (TNF-driven)")
    ggplotly(p) %>% layout(hovermode = "x unified")
  })

  ## ── Plot: DAS28 ──────────────────────────────────────────
  output$plot_das28 <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(time, DAS28_CRP_out)) +
      geom_line(color = "#c0392b", linewidth = 1.2) +
      geom_ribbon(aes(ymin = 0, ymax = 2.6), fill = "#27ae60", alpha = 0.10) +
      geom_ribbon(aes(ymin = 2.6, ymax = 3.2), fill = "#f1c40f", alpha = 0.10) +
      geom_ribbon(aes(ymin = 3.2, ymax = 5.1), fill = "#e67e22", alpha = 0.10) +
      geom_ribbon(aes(ymin = 5.1, ymax = 9.0), fill = "#e74c3c", alpha = 0.08) +
      geom_hline(yintercept = 2.6, linetype = 2, color = "#27ae60") +
      geom_hline(yintercept = 3.2, linetype = 2, color = "#f1c40f") +
      geom_hline(yintercept = 5.1, linetype = 2, color = "#e74c3c") +
      scale_x_continuous("Time (days)", breaks = seq(0, 365, 28)) +
      scale_y_continuous("DAS28-CRP", limits = c(0, 8)) +
      theme_bw() +
      labs(subtitle = "Green=Remission, Yellow=Low, Orange=Moderate, Red=High")
    ggplotly(p) %>% layout(hovermode = "x unified")
  })

  ## ── Plot: ACR ─────────────────────────────────────────────
  output$plot_acr <- renderPlotly({
    df <- sim_result() %>%
      select(time, pct_improvement, ACR20, ACR50, ACR70) %>%
      pivot_longer(ACR20:ACR70, names_to = "Threshold", values_to = "Achieved")

    p <- ggplot(df, aes(time, pct_improvement)) +
      geom_line(color = "#2c3e50", linewidth = 1, linetype = "dashed") +
      geom_hline(yintercept = 20, linetype = 2, color = "#27ae60", linewidth = 0.8) +
      geom_hline(yintercept = 50, linetype = 2, color = "#f39c12", linewidth = 0.8) +
      geom_hline(yintercept = 70, linetype = 2, color = "#e74c3c", linewidth = 0.8) +
      annotate("text", x = 0, y = 22, label = "ACR20", size = 3.5, color = "#27ae60", hjust = 0) +
      annotate("text", x = 0, y = 52, label = "ACR50", size = 3.5, color = "#f39c12", hjust = 0) +
      annotate("text", x = 0, y = 72, label = "ACR70", size = 3.5, color = "#e74c3c", hjust = 0) +
      scale_x_continuous("Time (days)", breaks = seq(0, 365, 28)) +
      scale_y_continuous("DAS28 Improvement (%)", limits = c(-5, 100)) +
      theme_bw() +
      labs(subtitle = "Percentage improvement from pre-treatment DAS28")
    ggplotly(p) %>% layout(hovermode = "x unified")
  })

  ## ── Table: Simulation ─────────────────────────────────────
  output$table_sim <- renderDT({
    df <- sim_result() %>%
      filter(time %% 7 < 0.6) %>%   # weekly snapshots
      mutate(Week = round(time / 7)) %>%
      select(Week,
             `TCZ (mg/L)` = C_TCZ_mgL,
             `IL6R Block (%)` = fR_blocked,
             `TNF-α (ng/mL)` = TNFa_out,
             `IL-6 (pg/mL)` = IL6_out,
             `CRP (mg/L)` = CRP_out,
             SJC28 = SJC28_out,
             TJC28 = TJC28_out,
             `DAS28-CRP` = DAS28_CRP_out,
             `Δ DAS28 (%)` = pct_improvement,
             ACR20, ACR50, ACR70,
             EULAR = EULAR_remission,
             `HAQ-DI` = HAQ_DI) %>%
      mutate(across(where(is.numeric), ~ round(.x, 2)),
             `IL6R Block (%)` = round(`IL6R Block (%)` * 100, 1))

    datatable(df,
      options = list(pageLength = 12, scrollX = TRUE,
                     dom = 'Bfrtip', buttons = c('csv','excel')),
      extensions = 'Buttons',
      rownames = FALSE,
      class = "table-striped table-bordered compact"
    ) %>%
      formatStyle("DAS28-CRP",
        backgroundColor = styleInterval(c(2.6, 3.2, 5.1),
          c("#d5f5e3","#a9dfbf","#fdebd0","#f5cba7")))
  })

  ## ── Sensitivity Analysis ─────────────────────────────────
  sens_result <- eventReactive(input$run_sens, {
    param_name <- input$sens_param
    fold_range <- input$sens_fold
    n_levels   <- input$sens_n

    nominal <- param(base_mod)[[param_name]]
    fold_seq <- seq(fold_range[1], fold_range[2], length.out = n_levels)
    param_vals <- nominal * fold_seq

    ev_sens <- make_events(
      tcz_dose_mg  = 162,
      tcz_interval = 14,
      tcz_route    = "SC",
      tcz_only     = FALSE,
      mtx_dose_mg  = 15
    )

    results_list <- lapply(seq_along(param_vals), function(i) {
      m_tmp <- param(base_mod, setNames(list(param_vals[i]), param_name))
      m_tmp <- init_ra(m_tmp)
      df <- mrgsim(m_tmp, events = ev_sens, end = 168, delta = 0.5,
                   carry_out = "evid") %>%
        as.data.frame() %>%
        filter(evid == 0)
      df$param_level <- sprintf("%.2f× (%.3g)", fold_seq[i], param_vals[i])
      df$fold        <- fold_seq[i]
      df
    })
    bind_rows(results_list)
  })

  output$plot_sens_das28 <- renderPlotly({
    df <- sens_result()
    p <- ggplot(df, aes(time, DAS28_CRP_out,
                        color = factor(round(fold, 2)),
                        group = factor(round(fold, 2)))) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 2.6, linetype = 2, color = "green4") +
      geom_hline(yintercept = 5.1, linetype = 2, color = "red3") +
      scale_color_viridis_d(name = "Fold change") +
      scale_x_continuous("Time (days)") +
      scale_y_continuous("DAS28-CRP", limits = c(1, 8)) +
      theme_bw() +
      labs(title = paste("Sensitivity:", input$sens_param))
    ggplotly(p)
  })

  output$plot_sens_crp <- renderPlotly({
    df <- sens_result()
    p <- ggplot(df, aes(time, CRP_out,
                        color = factor(round(fold, 2)))) +
      geom_line(linewidth = 0.9) +
      scale_color_viridis_d(name = "Fold") +
      scale_x_continuous("Time (days)") +
      scale_y_continuous("CRP (mg/L)") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_sens_pk <- renderPlotly({
    df <- sens_result() %>%
      group_by(fold) %>%
      summarise(Cmax = max(C_TCZ_mgL), Ctrough = min(C_TCZ_mgL[time > 14]))
    p <- ggplot(df) +
      geom_line(aes(fold, Cmax, color = "Cmax"), linewidth = 1) +
      geom_line(aes(fold, Ctrough, color = "Ctrough"), linewidth = 1, linetype = 2) +
      scale_x_continuous("Fold change of parameter") +
      scale_y_continuous("TCZ concentration (mg/L)") +
      scale_color_manual(values = c(Cmax = "#c0392b", Ctrough = "#8e44ad"), name = "") +
      theme_bw()
    ggplotly(p)
  })

  ## ── Mechanistic Map ──────────────────────────────────────
  output$map_image <- renderImage({
    list(
      src = normalizePath("../ra_qsp.png"),
      contentType = "image/png",
      width = "100%"
    )
  }, deleteFile = FALSE)

  output$dl_svg <- downloadHandler(
    filename = "ra_qsp_mechanistic_map.svg",
    content  = function(file) file.copy("../ra_qsp.svg", file)
  )
  output$dl_png <- downloadHandler(
    filename = "ra_qsp_mechanistic_map.png",
    content  = function(file) file.copy("../ra_qsp.png", file)
  )
  output$dl_dot <- downloadHandler(
    filename = "ra_qsp.dot",
    content  = function(file) file.copy("../ra_qsp.dot", file)
  )

  ## ── About: Parameter table ───────────────────────────────
  output$table_params <- renderTable({
    data.frame(
      Parameter = c("CL_TCZ","V1_TCZ","Q_TCZ","V2_TCZ",
                    "F_SC_TCZ","ka_SC_TCZ",
                    "kon_TCZ","koff_TCZ","KD",
                    "ksyn_R","kdeg_R",
                    "Emax_MTX","EC50_MTX",
                    "EC50_IL6_CRP","kout_CRP"),
      Value = c(0.224, 3.72, 0.697, 2.91,
                0.80, 0.29,
                0.0272, 0.0096, 0.35,
                0.132, 0.066,
                0.50, 0.05,
                8.0, 0.693),
      Unit = c("L/day","L","L/day","L",
               "—","1/day",
               "1/(nM·day)","1/day","nM",
               "nM/day","1/day",
               "—","mg/L",
               "pg/mL","1/day"),
      Source = c(rep("Gibiansky 2012",4),
                 rep("Levi 2018",2),
                 rep("TMDD fit",3),
                 rep("Simeoni 2019",2),
                 rep("Cronstein 2005",2),
                 rep("Levi 2018",2))
    )
  }, striped = TRUE, bordered = TRUE, hover = TRUE)

}

## ── Run app ──────────────────────────────────────────────────
shinyApp(ui, server)
