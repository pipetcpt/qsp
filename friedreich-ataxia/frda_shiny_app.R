# =============================================================================
# Friedreich Ataxia (FRDA) — Shiny dashboard
# 8 tabs: Patient · FXN/Genetics · Mito-Energetics · CNS · Cardiac · Pancreas ·
#         Therapy comparison · Biomarkers
# =============================================================================
# Run:
#   library(shiny); library(mrgsolve); library(dplyr); library(ggplot2);
#   library(DT); library(tidyr)
#   source("frda_mrgsolve_model.R")
#   shiny::runApp("frda_shiny_app.R")

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

if (!exists("mod")) source("frda_mrgsolve_model.R")

ui <- fluidPage(
  titlePanel("Friedreich Ataxia (FRDA) — QSP Dashboard"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient profile"),
      sliderInput("GAA1", "Short-allele GAA repeats", 100, 1700, 650, 50),
      sliderInput("GAA2", "Long-allele GAA repeats", 100, 1700, 850, 50),
      sliderInput("AAO",  "Age at onset (years)",   3, 40, 12, 1),
      radioButtons("SEX", "Sex", choices = c("F"=1, "M"=0), inline = TRUE),
      hr(),
      h4("Therapy"),
      checkboxInput("USE_OMAV", "Omaveloxolone 150 mg QD", TRUE),
      checkboxInput("USE_IDB",  "Idebenone 450 mg TID",     FALSE),
      checkboxInput("USE_DFP",  "Deferiprone 25 mg/kg/d",   FALSE),
      checkboxInput("USE_NOM",  "Nomlabofusp 50 mg SC QD",  FALSE),
      checkboxInput("USE_AAV",  "AAVrh10-FXN gene Tx (single IV)", FALSE),
      checkboxInput("USE_ACEi", "ACEi (cardiac)",           FALSE),
      hr(),
      sliderInput("SIM_DAYS", "Simulation horizon (days)", 30, 365*3, 365, 30),
      actionButton("RUN", "Run simulation", class = "btn-primary")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("1. Patient",       plotOutput("plot_overview", height = 380),
                                     DTOutput("tab_summary")),
        tabPanel("2. FXN/Genetics",  plotOutput("plot_fxn", height = 320),
                                     plotOutput("plot_fes", height = 280)),
        tabPanel("3. Mito-Energetics", plotOutput("plot_etc", height = 300),
                                       plotOutput("plot_ros", height = 280)),
        tabPanel("4. CNS",           plotOutput("plot_drg", height = 280),
                                     plotOutput("plot_cb",  height = 280),
                                     plotOutput("plot_mfars", height = 280)),
        tabPanel("5. Cardiac",       plotOutput("plot_lvmi", height = 320),
                                     plotOutput("plot_atp", height = 280)),
        tabPanel("6. Pancreas",      plotOutput("plot_bcell", height = 280),
                                     plotOutput("plot_glc", height = 280),
                                     plotOutput("plot_hba1c", height = 280)),
        tabPanel("7. Therapy comparison",
                                     plotOutput("plot_compare", height = 480),
                                     DTOutput("tab_compare")),
        tabPanel("8. Biomarkers / Safety",
                                     plotOutput("plot_alt", height = 280),
                                     plotOutput("plot_pk", height = 320),
                                     DTOutput("tab_biomarker"))
      )
    )
  )
)

server <- function(input, output, session) {

  sim_one <- function(patient, regimen, days) {
    e <- ev(amt = 0, cmt = "OMAV_GUT")
    if (regimen$omav)  e <- e + ev(amt = 150,  cmt = "OMAV_GUT", ii = 24, addl = days-1)
    if (regimen$idb)   e <- e + ev(amt = 450,  cmt = "IDB_GUT",  ii = 8,  addl = days*3-1)
    if (regimen$dfp)   e <- e + ev(amt = 800,  cmt = "DFP_GUT",  ii = 12, addl = days*2-1)
    if (regimen$nom)   e <- e + ev(amt = 50,   cmt = "NOM_SC",   ii = 24, addl = days-1)
    mod %>%
      param(
        GAA1 = patient$GAA1, GAA2 = patient$GAA2,
        AAO  = patient$AAO,  SEX  = as.numeric(patient$SEX),
        AAV_dose_flag = if (regimen$aav) 1 else 0,
        ACEi_eff      = if (regimen$ace) 1 else 0
      ) %>%
      mrgsim(events = e, end = days, delta = 1) %>%
      as.data.frame()
  }

  sim <- eventReactive(input$RUN, {
    patient <- list(GAA1 = input$GAA1, GAA2 = input$GAA2,
                    AAO  = input$AAO,  SEX  = input$SEX)
    regimen <- list(omav = input$USE_OMAV, idb = input$USE_IDB,
                    dfp  = input$USE_DFP,  nom = input$USE_NOM,
                    aav  = input$USE_AAV,  ace = input$USE_ACEi)
    sim_one(patient, regimen, input$SIM_DAYS)
  }, ignoreNULL = FALSE)

  # ----- Tab 1: overview -----
  output$plot_overview <- renderPlot({
    d <- sim()
    df <- d %>%
      select(time, FXN, ETC, mFARS, LVMI, Glucose) %>%
      pivot_longer(-time)
    ggplot(df, aes(time, value)) +
      geom_line(linewidth = 1, color = "steelblue") +
      facet_wrap(~name, scales = "free_y") +
      theme_minimal(base_size = 12) +
      labs(x = "Day", y = NULL, title = "Disease trajectory")
  })

  output$tab_summary <- renderDT({
    d <- sim()
    last <- tail(d, 1)
    data.frame(
      Metric = c("FXN (%WT)", "Fe-S pool (%)", "ETC (%)", "ROS (a.u.)",
                 "DRG (%)", "Cerebellar fn (%)", "LVMI (g/m²)",
                 "β-cell (%)", "Glucose (mmol/L)", "HbA1c (%)",
                 "mFARS", "T25FW (s)"),
      Day_0  = round(c(d$FXN[1], d$FeS[1], d$ETC[1], d$ROS[1], d$DRG_pool[1],
                       d$CB_func[1], d$LVMI[1], d$Bcell[1], d$Glucose[1],
                       d$HbA1c[1], d$mFARS[1], d$T25FW[1]), 2),
      Final  = round(c(last$FXN, last$FeS, last$ETC, last$ROS, last$DRG_pool,
                       last$CB_func, last$LVMI, last$Bcell, last$Glucose,
                       last$HbA1c, last$mFARS, last$T25FW), 2)
    )
  }, options = list(dom = "t", pageLength = 20))

  # ----- Tab 2 -----
  output$plot_fxn <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, FXN)) + geom_line(color = "darkred", linewidth = 1) +
      labs(x = "Day", y = "Frataxin (%WT)", title = "Frataxin trajectory") +
      theme_minimal(base_size = 12) +
      geom_hline(yintercept = 35, linetype = 2, color = "grey") +
      annotate("text", x = 0, y = 37, label = "Carrier ~50%", hjust = 0, size = 3)
  })

  output$plot_fes <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, FeS)) + geom_line(color = "darkorange", linewidth = 1) +
      labs(x = "Day", y = "Fe-S cluster pool (%)", title = "Mitochondrial Fe-S biogenesis") +
      theme_minimal(base_size = 12)
  })

  # ----- Tab 3 -----
  output$plot_etc <- renderPlot({
    d <- sim() %>% select(time, ETC, ATP) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, color = name)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = "Capacity", title = "ETC & ATP") +
      theme_minimal(base_size = 12)
  })

  output$plot_ros <- renderPlot({
    d <- sim() %>% select(time, ROS, AOX) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, color = name)) + geom_line(linewidth = 1) +
      labs(x = "Day", y = "a.u.", title = "ROS vs antioxidant capacity") +
      theme_minimal(base_size = 12)
  })

  # ----- Tab 4 -----
  output$plot_drg <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, DRG_pool)) + geom_line(color = "purple", linewidth = 1) +
      labs(x = "Day", y = "DRG (%)", title = "Dorsal root ganglion survival") +
      theme_minimal(base_size = 12)
  })
  output$plot_cb <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, CB_func)) + geom_line(color = "purple4", linewidth = 1) +
      labs(x = "Day", y = "Cerebellar function (%)", title = "Cerebellar / dentate function") +
      theme_minimal(base_size = 12)
  })
  output$plot_mfars <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, mFARS)) + geom_line(color = "tomato", linewidth = 1.1) +
      labs(x = "Day", y = "mFARS", title = "Clinical ataxia score (mFARS)") +
      theme_minimal(base_size = 12)
  })

  # ----- Tab 5 -----
  output$plot_lvmi <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, LVMI)) + geom_line(color = "firebrick", linewidth = 1) +
      geom_hline(yintercept = 95, linetype = 2, color = "grey40") +
      labs(x = "Day", y = "LVMI (g/m²)", title = "LV mass — cardiomyopathy") +
      theme_minimal(base_size = 12)
  })
  output$plot_atp <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, ATP)) + geom_line(color = "darkgreen", linewidth = 1) +
      labs(x = "Day", y = "ATP (a.u.)", title = "Cardiac energetics") +
      theme_minimal(base_size = 12)
  })

  # ----- Tab 6 -----
  output$plot_bcell <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, Bcell)) + geom_line(color = "tan4", linewidth = 1) +
      labs(x = "Day", y = "β-cell mass (%)", title = "Pancreatic β-cell pool") +
      theme_minimal(base_size = 12)
  })
  output$plot_glc <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, Glucose)) + geom_line(color = "tan2", linewidth = 1) +
      labs(x = "Day", y = "Glucose (mmol/L)", title = "Plasma glucose") +
      theme_minimal(base_size = 12)
  })
  output$plot_hba1c <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, HbA1c)) + geom_line(color = "chocolate", linewidth = 1) +
      geom_hline(yintercept = 6.5, linetype = 2, color = "red") +
      labs(x = "Day", y = "HbA1c (%)", title = "HbA1c") +
      theme_minimal(base_size = 12)
  })

  # ----- Tab 7 -----
  output$plot_compare <- renderPlot({
    patient <- list(GAA1 = input$GAA1, GAA2 = input$GAA2,
                    AAO = input$AAO, SEX = input$SEX)
    regs <- list(
      "Natural history" = list(omav=FALSE, idb=FALSE, dfp=FALSE, nom=FALSE, aav=FALSE, ace=FALSE),
      "Omaveloxolone"   = list(omav=TRUE,  idb=FALSE, dfp=FALSE, nom=FALSE, aav=FALSE, ace=FALSE),
      "Idebenone"       = list(omav=FALSE, idb=TRUE,  dfp=FALSE, nom=FALSE, aav=FALSE, ace=FALSE),
      "Deferiprone"     = list(omav=FALSE, idb=FALSE, dfp=TRUE,  nom=FALSE, aav=FALSE, ace=FALSE),
      "Nomlabofusp"     = list(omav=FALSE, idb=FALSE, dfp=FALSE, nom=TRUE,  aav=FALSE, ace=FALSE),
      "AAV gene Tx"     = list(omav=FALSE, idb=FALSE, dfp=FALSE, nom=FALSE, aav=TRUE,  ace=FALSE),
      "Omav + ACEi"     = list(omav=TRUE,  idb=FALSE, dfp=FALSE, nom=FALSE, aav=FALSE, ace=TRUE)
    )
    out <- lapply(names(regs), function(nm) {
      d <- sim_one(patient, regs[[nm]], input$SIM_DAYS)
      d$regimen <- nm; d
    }) %>% bind_rows()
    df <- out %>% select(time, regimen, mFARS, LVMI, HbA1c, FXN) %>%
      pivot_longer(c(mFARS, LVMI, HbA1c, FXN))
    ggplot(df, aes(time, value, color = regimen)) +
      geom_line(linewidth = 1) +
      facet_wrap(~name, scales = "free_y") +
      theme_minimal(base_size = 12) +
      labs(x = "Day", y = NULL, title = "Therapy comparison")
  })

  output$tab_compare <- renderDT({
    patient <- list(GAA1 = input$GAA1, GAA2 = input$GAA2,
                    AAO = input$AAO, SEX = input$SEX)
    regs <- list(
      "Natural history" = list(omav=FALSE, idb=FALSE, dfp=FALSE, nom=FALSE, aav=FALSE, ace=FALSE),
      "Omaveloxolone"   = list(omav=TRUE,  idb=FALSE, dfp=FALSE, nom=FALSE, aav=FALSE, ace=FALSE),
      "Idebenone"       = list(omav=FALSE, idb=TRUE,  dfp=FALSE, nom=FALSE, aav=FALSE, ace=FALSE),
      "Deferiprone"     = list(omav=FALSE, idb=FALSE, dfp=TRUE,  nom=FALSE, aav=FALSE, ace=FALSE),
      "Nomlabofusp"     = list(omav=FALSE, idb=FALSE, dfp=FALSE, nom=TRUE,  aav=FALSE, ace=FALSE),
      "AAV gene Tx"     = list(omav=FALSE, idb=FALSE, dfp=FALSE, nom=FALSE, aav=TRUE,  ace=FALSE),
      "Omav + ACEi"     = list(omav=TRUE,  idb=FALSE, dfp=FALSE, nom=FALSE, aav=FALSE, ace=TRUE)
    )
    do.call(rbind, lapply(names(regs), function(nm) {
      d <- sim_one(patient, regs[[nm]], input$SIM_DAYS)
      l <- tail(d, 1)
      data.frame(
        Regimen = nm,
        FXN     = round(l$FXN, 1),
        mFARS   = round(l$mFARS, 1),
        LVMI    = round(l$LVMI, 0),
        HbA1c   = round(l$HbA1c, 2),
        DRG_pct = round(l$DRG_pool, 1),
        ALT     = round(l$ALT_OMAV, 0)
      )
    }))
  }, options = list(dom = "t"))

  # ----- Tab 8 -----
  output$plot_alt <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, ALT_OMAV)) + geom_line(color = "olivedrab", linewidth = 1) +
      geom_hline(yintercept = 100, linetype = 2, color = "red") +
      labs(x = "Day", y = "ALT (U/L)", title = "Omaveloxolone class-effect ALT") +
      theme_minimal(base_size = 12)
  })
  output$plot_pk <- renderPlot({
    d <- sim()
    ggplot(d, aes(time, OMAV_CP)) + geom_line(color = "navy", linewidth = 1) +
      labs(x = "Day", y = "Omav plasma (ng/mL)", title = "Omaveloxolone PK") +
      theme_minimal(base_size = 12)
  })
  output$tab_biomarker <- renderDT({
    d <- sim()
    l <- tail(d, 1)
    data.frame(
      Biomarker = c("Plasma FXN (proxy %WT)", "Fe-S pool (%)", "ETC (%)",
                    "ROS (a.u.)", "AOX (a.u.)", "LVMI (g/m²)",
                    "HbA1c (%)", "mFARS", "T25FW (s)", "ALT (U/L)"),
      Value = round(c(l$FXN, l$FeS, l$ETC, l$ROS, l$AOX, l$LVMI,
                      l$HbA1c, l$mFARS, l$T25FW, l$ALT_OMAV), 2),
      Reference = c("WT 100%; carriers ~50%", "≥80% normal", "≥80% normal",
                    "≈0.20 baseline", "≈1.0 baseline", "<95 normal",
                    "<6.5% non-DM", "0 (none) – 93 (max)", "5 s baseline", "<55 U/L")
    )
  }, options = list(dom = "t", pageLength = 20))
}

shinyApp(ui, server)
