# =====================================================================
# Carcinoid Syndrome QSP — Shiny dashboard
# Date  : 2026-06-30
#
# UI :
#   Tab 1 Patient / Tumor profile (site, grade, Ki-67, SSTR PET, CHD)
#   Tab 2 Drug PK & SSTR occupancy (octreotide/lanreotide/pasireotide)
#   Tab 3 Tryptophan→5-HTP→5-HT→5-HIAA pathway viewer
#   Tab 4 Symptom dynamics (BM/day, flushing/day, bronchospasm)
#   Tab 5 Tumor burden & PRRT/TKI/mTORi growth modulation
#   Tab 6 Carcinoid heart disease progression (5-HT2B → valve fibrosis)
#   Tab 7 Scenario comparator (12 prebuilt regimens)
#   Tab 8 Biomarkers (CgA, NT-proBNP, 24-h urinary 5-HIAA)
# =====================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(plotly)
  library(DT)
  library(dplyr)
  library(tidyr)
})

source("carcsyn_mrgsolve_model.R")

ui <- dashboardPage(
  dashboardHeader(title = "Carcinoid QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient profile",       tabName = "pt",     icon = icon("user")),
      menuItem("Drug PK & SSTR",        tabName = "pk",     icon = icon("syringe")),
      menuItem("Serotonin pathway",     tabName = "ser",    icon = icon("dna")),
      menuItem("Symptoms (BM / flush)", tabName = "symp",   icon = icon("heart-pulse")),
      menuItem("Tumor burden",          tabName = "tumor",  icon = icon("circle-radiation")),
      menuItem("Carcinoid heart dz",    tabName = "chd",    icon = icon("heart")),
      menuItem("Scenarios (n=12)",      tabName = "scen",   icon = icon("flask")),
      menuItem("Biomarkers",            tabName = "bio",    icon = icon("vial"))
    ),
    hr(),
    h5("Patient covariates", style = "padding-left:15px;"),
    sliderInput("age",   "Age (yr)",     30, 90, 63),
    sliderInput("wt",    "Weight (kg)",  40, 130, 72),
    sliderInput("ki67",  "Ki-67 (%)",     1, 25, 6),
    sliderInput("liver", "Liver burden (%)", 5, 90, 35),
    selectInput("site",  "Primary site",
                choices = c("Midgut (ileal)", "Foregut/bronchial",
                            "Hindgut", "Unknown primary")),
    selectInput("ssta",  "SSTR2 PET Krenning",
                choices = c("Grade 4 (strong)" = 4,
                            "Grade 3" = 3,
                            "Grade 2" = 2,
                            "Grade <2 (PRRT contraindicated)" = 1))
  ),
  dashboardBody(
    tabItems(
      tabItem("pt",
        fluidRow(
          box(width=6, title="Patient summary", status="primary", solidHeader=TRUE,
              verbatimTextOutput("ptsummary")),
          box(width=6, title="Tumor grade & SSTR PET", status="info", solidHeader=TRUE,
              plotlyOutput("ptplot"))
        )),
      tabItem("pk",
        fluidRow(
          box(width=6, title="Octreotide / Lanreotide / Pasireotide",
              plotlyOutput("pk_ssa")),
          box(width=6, title="SSTR2/5 fractional occupancy",
              plotlyOutput("pk_occ"))
        )),
      tabItem("ser",
        fluidRow(
          box(width=12, title="Tryptophan → 5-HTP → 5-HT → 5-HIAA",
              plotlyOutput("ser_path"))
        )),
      tabItem("symp",
        fluidRow(
          box(width=6, title="Bowel movements / day", plotlyOutput("bm_plot")),
          box(width=6, title="Flushing episodes / day", plotlyOutput("flush_plot"))
        )),
      tabItem("tumor",
        fluidRow(
          box(width=12, title="Hepatic tumor burden under treatment",
              plotlyOutput("tumor_plot"))
        )),
      tabItem("chd",
        fluidRow(
          box(width=6, title="Valve collagen / Hassan equivalent",
              plotlyOutput("valve_plot")),
          box(width=6, title="NT-proBNP", plotlyOutput("bnp_plot"))
        )),
      tabItem("scen",
        fluidRow(
          box(width=12, status="primary", solidHeader=TRUE,
              title="Scenario library (TELESTAR, CLARINET, NETTER-1/2, RADIANT-4)",
              checkboxGroupInput("scenpick", "Select scenarios to compare",
                choices = c("S1_natural_history",
                            "S2_octreotide_LAR_30mg",
                            "S3_lanreotide_autogel_120mg",
                            "S4_octreotide_plus_telotristat",
                            "S5_pasireotide_LAR",
                            "S6_everolimus_10mg",
                            "S7_sunitinib_proxy_TKI",
                            "S8_PRRT_177Lu_DOTATATE",
                            "S9_IFN_alpha",
                            "S10_HAE_then_octreotide",
                            "S11_carcinoid_crisis_prevention",
                            "S12_quad_therapy"),
                selected = c("S1_natural_history",
                             "S2_octreotide_LAR_30mg",
                             "S4_octreotide_plus_telotristat",
                             "S8_PRRT_177Lu_DOTATATE"),
                inline = TRUE)),
          box(width=12, plotlyOutput("scen_compare"))
        )),
      tabItem("bio",
        fluidRow(
          box(width=4, title="Plasma 5-HT (nmol/L)",  plotlyOutput("bio_ser")),
          box(width=4, title="Urinary 5-HIAA (mg/24h)", plotlyOutput("bio_hiaa")),
          box(width=4, title="Platelet 5-HT (ng/mL)", plotlyOutput("bio_plt"))
        ),
        fluidRow(
          box(width=12, title="Biomarker summary table",
              DTOutput("bio_table"))
        ))
    )
  )
)

server <- function(input, output, session) {

  out_sim <- reactive({
    req(mod)
    scen <- build_scenarios()
    pick <- input$scenpick %||% names(scen)[1:4]
    run_scenarios(mod, scen[pick])
  })

  output$ptsummary <- renderPrint({
    cat("Age           :", input$age, "yr\n")
    cat("Weight        :", input$wt, "kg\n")
    cat("Primary site  :", input$site, "\n")
    cat("Ki-67         :", input$ki67, "%\n")
    cat("Hepatic load  :", input$liver, "%\n")
    cat("SSTR Krenning :", input$ssta, "\n")
    cat("PRRT eligible :", if (as.numeric(input$ssta) >= 3) "YES" else "NO", "\n")
  })

  output$ptplot <- renderPlotly({
    plot_ly(
      x = c("Ki-67 (%)", "Liver replacement (%)", "Krenning"),
      y = c(input$ki67, input$liver, as.numeric(input$ssta)*10),
      type = "bar", marker = list(color = c("#EF5350","#FFA726","#42A5F5"))
    ) %>% layout(yaxis = list(title=""))
  })

  output$pk_ssa <- renderPlotly({
    s <- out_sim(); req(nrow(s)>0)
    plot_ly(s, x = ~time/24, y = ~C_OCT, color = ~scenario, type="scatter", mode="lines",
            name="Octreotide") %>%
      add_trace(y = ~C_LAN, name = "Lanreotide") %>%
      add_trace(y = ~C_PAS, name = "Pasireotide") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Drug (ng/mL)"))
  })

  output$pk_occ <- renderPlotly({
    s <- out_sim(); req(nrow(s)>0)
    plot_ly(s, x=~time/24, y=~SSTR_TOTAL, color=~scenario,
            type="scatter", mode="lines") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="SSTR2/5 occupancy"))
  })

  output$ser_path <- renderPlotly({
    s <- out_sim(); req(nrow(s)>0)
    plot_ly(s, x=~time/24) %>%
      add_lines(y=~SER_P,   name="Plasma 5-HT (nmol/L)") %>%
      add_lines(y=~SER_PLT/100, name="Platelet 5-HT /100 (ng/mL)") %>%
      add_lines(y=~HIAA_U,  name="Urinary 5-HIAA (mg/24h)") %>%
      layout(xaxis=list(title="Day"), yaxis=list(title="Marker"))
  })

  output$bm_plot    <- renderPlotly(plot_ly(out_sim(), x=~time/24, y=~BM,
                                            color=~scenario, type="scatter", mode="lines"))
  output$flush_plot <- renderPlotly(plot_ly(out_sim(), x=~time/24, y=~FLUSH,
                                            color=~scenario, type="scatter", mode="lines"))
  output$tumor_plot <- renderPlotly(plot_ly(out_sim(), x=~time/24, y=~TUMOR,
                                            color=~scenario, type="scatter", mode="lines"))
  output$valve_plot <- renderPlotly(plot_ly(out_sim(), x=~time/24, y=~VALVE,
                                            color=~scenario, type="scatter", mode="lines"))
  output$bnp_plot   <- renderPlotly(plot_ly(out_sim(), x=~time/24, y=~NTproBNP,
                                            color=~scenario, type="scatter", mode="lines"))

  output$scen_compare <- renderPlotly({
    s <- out_sim(); req(nrow(s)>0)
    summ <- s %>% group_by(scenario) %>% summarize(
      delta_BM    = round(mean(BM[time>24*180])    - mean(BM[time<24*30]), 2),
      delta_FLUSH = round(mean(FLUSH[time>24*180]) - mean(FLUSH[time<24*30]), 2),
      delta_TUMOR_pct = round(100*(mean(TUMOR[time>24*180]) -
                                    mean(TUMOR[time<24*30]))/mean(TUMOR[time<24*30]), 1)
    )
    plot_ly(summ, x=~scenario, y=~delta_BM, type="bar", name="ΔBM/day") %>%
      add_trace(y=~delta_FLUSH, name="ΔFlush/day") %>%
      add_trace(y=~delta_TUMOR_pct, name="ΔTumor (%)")
  })

  output$bio_ser  <- renderPlotly(plot_ly(out_sim(), x=~time/24, y=~SER_P,
                                          color=~scenario, type="scatter", mode="lines"))
  output$bio_hiaa <- renderPlotly(plot_ly(out_sim(), x=~time/24, y=~HIAA_U,
                                          color=~scenario, type="scatter", mode="lines"))
  output$bio_plt  <- renderPlotly(plot_ly(out_sim(), x=~time/24, y=~SER_PLT,
                                          color=~scenario, type="scatter", mode="lines"))
  output$bio_table <- renderDT({
    s <- out_sim()
    s %>% group_by(scenario) %>% summarize(
      max_5HT = round(max(SER_P), 2),
      min_5HT = round(min(SER_P), 2),
      mean_5HIAA = round(mean(HIAA_U), 2),
      mean_BM = round(mean(BM), 2),
      mean_FLUSH = round(mean(FLUSH), 2),
      end_TUMOR = round(tail(TUMOR,1), 1),
      end_NTproBNP = round(tail(NTproBNP,1), 1)
    )
  })
}

`%||%` <- function(a,b) if (is.null(a)) b else a

# shinyApp(ui, server)
