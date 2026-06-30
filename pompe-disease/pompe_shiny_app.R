## =============================================================
## Pompe Disease QSP — Shiny dashboard skeleton
##
## Run with:  shiny::runApp("pompe_shiny_app.R")
##
## Requires:  shiny, shinydashboard, mrgsolve, ggplot2, dplyr, DT
## Loads ODE model from pompe_mrgsolve_model.R (sourced lazily).
## =============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)

source("pompe_mrgsolve_model.R", local = FALSE)

scenarios <- c("No treatment"         = "no_tx",
               "Alglucosidase alfa"   = "alglu",
               "Avalglucosidase alfa" = "aval",
               "Cipaglucosidase + Miglustat" = "cipa_mig",
               "AAV9-GAA gene therapy" = "aav_gt",
               "Alglucosidase + ITI (CRIM-)" = "alglu_iti")

ui <- dashboardPage(
  dashboardHeader(title = "Pompe Disease QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient profile",     tabName = "patient",  icon = icon("user")),
      menuItem("Pharmacokinetics",    tabName = "pk",       icon = icon("chart-line")),
      menuItem("PD: tissue & glycogen", tabName = "pd",     icon = icon("dna")),
      menuItem("Clinical endpoints",  tabName = "endpoints",icon = icon("walking")),
      menuItem("Scenario comparison", tabName = "compare",  icon = icon("balance-scale")),
      menuItem("Biomarkers & safety", tabName = "biomark",  icon = icon("vial")),
      menuItem("ADA / immunology",    tabName = "ada",      icon = icon("shield-virus")),
      menuItem("References",          tabName = "refs",     icon = icon("book"))
    )
  ),
  dashboardBody(
    tabItems(
      ## ---- 1. Patient profile ----
      tabItem("patient",
        fluidRow(
          box(width = 4, title = "Patient & disease setup", status = "primary",
              radioButtons("phenotype", "Phenotype",
                           choices = c("LOPD adult" = "LOPD",
                                       "IOPD infantile" = "IOPD")),
              numericInput("wt", "Body weight (kg)", value = 70, min = 3, max = 120),
              numericInput("gaa_base", "Residual GAA activity (frac. of normal)",
                           value = 0.10, min = 0, max = 0.4, step = 0.01),
              checkboxInput("crim_neg", "CRIM-negative (high-ADA risk)", value = FALSE),
              selectInput("scenario", "Treatment scenario", choices = scenarios,
                          selected = "alglu"),
              sliderInput("years", "Simulation horizon (years)", min = 0.5, max = 10,
                          value = 3, step = 0.5),
              actionButton("go", "Run simulation", class = "btn-primary")
          ),
          box(width = 8, title = "Phenotype summary", status = "info",
              htmlOutput("phenotype_summary"))
        )
      ),
      ## ---- 2. PK ----
      tabItem("pk",
        fluidRow(
          box(width = 12, title = "Plasma drug concentrations",
              plotOutput("pk_plot", height = 380))
        )
      ),
      ## ---- 3. PD ----
      tabItem("pd",
        fluidRow(
          box(width = 6, title = "Lysosomal GAA pool (tissue)",
              plotOutput("pd_gaa", height = 320)),
          box(width = 6, title = "Lysosomal glycogen",
              plotOutput("pd_glyc", height = 320))
        )
      ),
      ## ---- 4. Clinical endpoints ----
      tabItem("endpoints",
        fluidRow(
          box(width = 6, title = "FVC upright (% predicted)",
              plotOutput("ep_fvc", height = 280)),
          box(width = 6, title = "6-minute walk distance",
              plotOutput("ep_smwt", height = 280)),
          box(width = 6, title = "LV mass index",
              plotOutput("ep_lvmi", height = 280)),
          box(width = 6, title = "Ventilator failure hazard",
              plotOutput("ep_vent", height = 280))
        )
      ),
      ## ---- 5. Scenario comparison ----
      tabItem("compare",
        fluidRow(
          box(width = 12, title = "Multi-scenario simulation",
              checkboxGroupInput("comp_scen", "Scenarios",
                                 choices = scenarios,
                                 selected = c("no_tx","alglu","aval","cipa_mig")),
              actionButton("go_compare", "Run comparison"),
              plotOutput("compare_plot", height = 460),
              DT::DTOutput("compare_table"))
        )
      ),
      ## ---- 6. Biomarkers ----
      tabItem("biomark",
        fluidRow(
          box(width = 6, title = "Plasma Hex4",
              plotOutput("bm_hex4", height = 260)),
          box(width = 6, title = "Serum CK",
              plotOutput("bm_ck", height = 260)),
          box(width = 6, title = "NT-proBNP",
              plotOutput("bm_ntbnp", height = 260)),
          box(width = 6, title = "LV ejection fraction",
              plotOutput("bm_ef", height = 260))
        )
      ),
      ## ---- 7. ADA / immunology ----
      tabItem("ada",
        fluidRow(
          box(width = 6, title = "Anti-GAA ADA titre",
              plotOutput("ada_titre", height = 280)),
          box(width = 6, title = "Tissue uptake block (%)",
              plotOutput("ada_block", height = 280))
        )
      ),
      ## ---- 8. References ----
      tabItem("refs",
        includeMarkdown("pompe_references.md")
      )
    )
  )
)

server <- function(input, output, session) {

  sim <- eventReactive(input$go, {
    iopd <- input$phenotype == "IOPD"
    out <- pompe_run(input$scenario, iopd = iopd, years = input$years)
    out$WT       <- input$wt
    out$gaa_base <- input$gaa_base
    out$crim_neg <- input$crim_neg
    out
  })

  comp_sim <- eventReactive(input$go_compare, {
    iopd <- input$phenotype == "IOPD"
    lapply(input$comp_scen, function(s) pompe_run(s, iopd = iopd, years = input$years)) |>
      bind_rows()
  })

  output$phenotype_summary <- renderUI({
    HTML(paste0(
      "<p><b>Phenotype:</b> ", input$phenotype,
      " · <b>Body weight:</b> ", input$wt, " kg</p>",
      "<p><b>Residual GAA:</b> ", round(input$gaa_base*100,1), "%",
      " · <b>CRIM:</b> ", ifelse(input$crim_neg, "negative", "positive"),
      " · <b>Scenario:</b> ", names(scenarios)[scenarios == input$scenario], "</p>",
      "<p><i>Press 'Run simulation' to populate the dashboard.</i></p>"
    ))
  })

  ## ---- PK plot ----
  output$pk_plot <- renderPlot({
    df <- sim()
    pk <- df |>
      select(time, Cp_alglu, Cp_aval, Cp_cipa, Cp_mig) |>
      tidyr::pivot_longer(-time, names_to = "drug", values_to = "Cp")
    ggplot(pk, aes(time, Cp, color = drug)) +
      geom_line(linewidth = 1) +
      facet_wrap(~drug, scales = "free_y") +
      labs(x = "Time (days)", y = "Plasma concentration (mg/L)") +
      theme_bw() + theme(legend.position = "none")
  })

  output$pd_gaa <- renderPlot({
    df <- sim()
    long <- df |> select(time, GAA_M, GAA_C, GAA_D) |>
      tidyr::pivot_longer(-time, names_to = "tissue", values_to = "gaa")
    ggplot(long, aes(time, gaa, color = tissue)) +
      geom_line(linewidth = 1) +
      labs(x = "Time (days)", y = "Lysosomal GAA (a.u.)") +
      theme_bw()
  })

  output$pd_glyc <- renderPlot({
    df <- sim()
    long <- df |> select(time, GLYC_M, GLYC_C, GLYC_D) |>
      tidyr::pivot_longer(-time, names_to = "tissue", values_to = "glyc")
    ggplot(long, aes(time, glyc, color = tissue)) +
      geom_line(linewidth = 1) +
      labs(x = "Time (days)", y = "Lysosomal glycogen (a.u.)") +
      theme_bw()
  })

  output$ep_fvc  <- renderPlot({ ggplot(sim(), aes(time, FVC_UP)) + geom_line() + theme_bw() })
  output$ep_smwt <- renderPlot({ ggplot(sim(), aes(time, SMWT))  + geom_line() + theme_bw() })
  output$ep_lvmi <- renderPlot({ ggplot(sim(), aes(time, LVMI))  + geom_line() + theme_bw() })
  output$ep_vent <- renderPlot({ ggplot(sim(), aes(time, VENT_RISK)) + geom_line() + theme_bw() })

  output$compare_plot <- renderPlot({
    df <- comp_sim()
    long <- df |> select(time, scenario, FVC_UP, LVMI, SMWT) |>
      tidyr::pivot_longer(c(FVC_UP, LVMI, SMWT), names_to = "metric", values_to = "value")
    ggplot(long, aes(time, value, color = scenario)) +
      geom_line(linewidth = 1) +
      facet_wrap(~metric, scales = "free_y") + theme_bw()
  })

  output$compare_table <- DT::renderDT({
    df <- comp_sim()
    df |> group_by(scenario) |>
      summarise(last_FVC = round(tail(FVC_UP,1),1),
                last_LVMI = round(tail(LVMI,1),1),
                last_SMWT = round(tail(SMWT,1),0),
                last_HEX4 = round(tail(HEX4,1),2),
                last_ADA  = round(tail(ADA_T,1),1))
  })

  output$bm_hex4  <- renderPlot({ ggplot(sim(), aes(time, HEX4))    + geom_line() + theme_bw() })
  output$bm_ck    <- renderPlot({ ggplot(sim(), aes(time, CK))      + geom_line() + theme_bw() })
  output$bm_ntbnp <- renderPlot({ ggplot(sim(), aes(time, NTproBNP))+ geom_line() + theme_bw() })
  output$bm_ef    <- renderPlot({ ggplot(sim(), aes(time, EF_LV))   + geom_line() + theme_bw() })

  output$ada_titre <- renderPlot({ ggplot(sim(), aes(time, ADA_T))    + geom_line() + theme_bw() })
  output$ada_block <- renderPlot({ ggplot(sim(), aes(time, ada_block*100)) + geom_line() +
      labs(y = "Uptake block (%)") + theme_bw() })
}

shinyApp(ui, server)
