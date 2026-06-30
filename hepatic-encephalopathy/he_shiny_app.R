## =====================================================================
## Hepatic Encephalopathy (HE) — QSP Shiny App
## ---------------------------------------------------------------------
## 8 tabs:
##   1. Patient Profile / Phenotype
##   2. Drug PK (Lactulose, Rifaximin, LOLA, BCAA, Benzoate, GPB,
##      Albumin, Flumazenil)
##   3. Gut-Liver-Systemic NH3 kinetics
##   4. Brain NH3 / Glutamine / Astrocyte swelling
##   5. Neuroinflammation (LPS, TNFα, GABA-A tone, Mn)
##   6. Clinical endpoints (West Haven, PHES surrogate, mortality hazard)
##   7. Scenario comparison (9 treatment ladders)
##   8. Biomarkers / Fischer ratio / Trial overlay
## =====================================================================
library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("he_mrgsolve_model.R", local = TRUE, chdir = TRUE)
# mod, all_sc available

ui <- dashboardPage(
  dashboardHeader(title = "Hepatic Encephalopathy QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient Profile",  tabName = "profile",   icon = icon("user")),
      menuItem("2. Drug PK",          tabName = "pk",        icon = icon("pills")),
      menuItem("3. Gut-Liver NH3",    tabName = "ammonia",   icon = icon("flask")),
      menuItem("4. Brain / Astrocyte",tabName = "brain",     icon = icon("brain")),
      menuItem("5. Neuroinflammation",tabName = "infl",      icon = icon("fire")),
      menuItem("6. Clinical Endpoints",tabName = "clinical", icon = icon("heartbeat")),
      menuItem("7. Scenarios",        tabName = "scenario",  icon = icon("layer-group")),
      menuItem("8. Biomarkers",       tabName = "bio",       icon = icon("chart-line"))
    ),
    hr(),
    h4("Patient Phenotype"),
    sliderInput("MELD",   "MELD score",                  6, 40, 20),
    sliderInput("fps",    "Portosystemic shunt fraction",0, 0.95, 0.40, 0.05),
    sliderInput("HEPMASS","Hepatocyte mass (functional)",0.05, 1.0, 0.45, 0.05),
    sliderInput("SARCO",  "Sarcopenia (0=none, 1=severe)",0, 1, 0.30, 0.05),
    sliderInput("PROTEIN","Protein intake (g/kg/d)",      0.5, 2.0, 1.2, 0.1),
    checkboxInput("ALF_FLAG","Acute liver failure (ALF)", FALSE),
    hr(),
    h4("Therapy"),
    sliderInput("LAC_DOSE","Lactulose (g/d)",   0, 90, 30, 5),
    sliderInput("RIF_DOSE","Rifaximin (mg/d)",  0, 1650, 1100, 50),
    checkboxInput("LOLA_on","LOLA 20 g IV/d",  FALSE),
    checkboxInput("BCAA_on","BCAA 12 g TID",    FALSE),
    checkboxInput("BZ_on",  "Na-Benzoate 5 g TID", FALSE),
    checkboxInput("GPB_on", "GPB 6 mL TID",     FALSE),
    checkboxInput("PROB_on","Probiotic / FMT",  FALSE),
    checkboxInput("ALB_on", "Albumin 50 g IV",  FALSE),
    checkboxInput("FLU_on", "Flumazenil 1 mg IV (rescue)", FALSE),
    sliderInput("SIM_DAYS","Simulation days",   3, 90, 30, 1)
  ),
  dashboardBody(
    tabItems(
      tabItem("profile",
        fluidRow(
          box(width=12, title="Patient phenotype summary", status="primary", solidHeader=TRUE,
              verbatimTextOutput("phenotype_summary"))),
        fluidRow(
          box(width=6, title="MELD / Child-Pugh visual", plotOutput("meld_plot", height=320)),
          box(width=6, title="Hepatocyte mass · sarcopenia · shunt", plotOutput("phenotype_bars", height=320)))
      ),
      tabItem("pk",
        fluidRow(
          box(width=6, title="LOLA / BCAA / Benzoate / GPB PK", plotOutput("pk_plot", height=320)),
          box(width=6, title="Albumin / Flumazenil PK", plotOutput("pk_alb_flu", height=320))),
        fluidRow(
          box(width=12, title="Gut pH (lactulose) & urease activity (rifaximin)", plotOutput("pk_gut", height=300)))
      ),
      tabItem("ammonia",
        fluidRow(
          box(width=6, title="Gut → Portal → Systemic NH3", plotOutput("nh3_axis", height=340)),
          box(width=6, title="Muscle GS & hepatic ureagenesis", plotOutput("nh3_clear", height=340))),
        fluidRow(
          box(width=12, title="NH3 timecourse vs trial benchmarks (Kircheis 1997 LOLA −38%)",
              plotOutput("nh3_trial", height=300)))
      ),
      tabItem("brain",
        fluidRow(
          box(width=6, title="Brain NH3 & Glutamine", plotOutput("brain_nh3", height=320)),
          box(width=6, title="Astrocyte swelling index", plotOutput("brain_swell", height=320))),
        fluidRow(
          box(width=12, title="Manganese basal ganglia accumulation",
              plotOutput("brain_mn", height=280)))
      ),
      tabItem("infl",
        fluidRow(
          box(width=6, title="LPS / TNFα", plotOutput("infl_lps", height=320)),
          box(width=6, title="GABA-A PAM tone", plotOutput("infl_gaba", height=320))),
        fluidRow(
          box(width=12, title="NH3 × inflammation synergy (Shawcross 2004)",
              plotOutput("infl_synergy", height=280)))
      ),
      tabItem("clinical",
        fluidRow(
          box(width=6, title="West Haven grade trajectory", plotOutput("wh_traj", height=340)),
          box(width=6, title="Cumulative mortality hazard", plotOutput("hazard", height=340))),
        fluidRow(
          box(width=12, title="Predicted minimal vs overt HE (PHES surrogate)",
              plotOutput("phes", height=280)))
      ),
      tabItem("scenario",
        fluidRow(
          box(width=12, title="9 treatment-ladder scenarios — West Haven, NH3, swelling overlay",
              plotOutput("scenario_panel", height=520))),
        fluidRow(
          box(width=12, title="Scenario comparison table at day 30",
              DT::DTOutput("scenario_table")))
      ),
      tabItem("bio",
        fluidRow(
          box(width=6, title="Fischer ratio (BCAA/AAA)", plotOutput("fischer", height=320)),
          box(width=6, title="Plasma glutamine, brain glutamine", plotOutput("gln_both", height=320))),
        fluidRow(
          box(width=12, title="Trial overlays (Bass 2010, Sharma 2013, Les 2011)",
              plotOutput("trial_overlay", height=300)))
      )
    )
  )
)

server <- function(input, output, session){

  sim_data <- reactive({
    end_t <- input$SIM_DAYS * 24
    p <- list(
      MELD=input$MELD, fps=input$fps, HEPMASS=input$HEPMASS, SARCO=input$SARCO,
      PROTEIN=input$PROTEIN, ALF_FLAG=as.numeric(input$ALF_FLAG),
      LAC_DOSE=input$LAC_DOSE, RIF_DOSE=input$RIF_DOSE
    )
    m <- param(mod, p)
    evlist <- list()
    if(input$LOLA_on) evlist[[length(evlist)+1]] <- ev(time=seq(0,end_t,24), amt=20000, cmt="LOLA_C")
    if(input$BCAA_on) evlist[[length(evlist)+1]] <- ev(time=seq(0,end_t,8),  amt=12,    cmt="BCAA_GUT")
    if(input$BZ_on)   evlist[[length(evlist)+1]] <- ev(time=seq(0,end_t,8),  amt=5000,  cmt="BZ_GUT")
    if(input$GPB_on)  evlist[[length(evlist)+1]] <- ev(time=seq(0,end_t,8),  amt=6000,  cmt="GPB_GUT")
    if(input$PROB_on) evlist[[length(evlist)+1]] <- ev(time=0,                amt=1,    cmt="PROB_E")
    if(input$ALB_on)  evlist[[length(evlist)+1]] <- ev(time=seq(0,end_t,12), amt=50,    cmt="ALB_C")
    if(input$FLU_on)  evlist[[length(evlist)+1]] <- ev(time=24,               amt=1,    cmt="FLU_C")
    events <- if(length(evlist)==0) ev(time=0, amt=0, cmt=1) else do.call(c, evlist)
    out <- mrgsim(m, events=events, end=end_t, delta=2) %>% as_tibble()
    out$day <- out$time/24
    out
  })

  output$phenotype_summary <- renderPrint({
    cat("Patient phenotype:\n")
    cat("  MELD              :", input$MELD, "\n")
    cat("  Shunt fraction    :", input$fps, "\n")
    cat("  Hepatocyte mass   :", input$HEPMASS, "\n")
    cat("  Sarcopenia        :", input$SARCO, "\n")
    cat("  Protein intake    :", input$PROTEIN, "g/kg/d\n")
    cat("  ALF flag          :", input$ALF_FLAG, "\n\n")
    cat("Therapy: Lact", input$LAC_DOSE, "g/d, Rif", input$RIF_DOSE, "mg/d\n")
    cat("  +LOLA:", input$LOLA_on, "  +BCAA:", input$BCAA_on, "\n")
    cat("  +Benzoate:", input$BZ_on, "  +GPB:", input$GPB_on, "\n")
    cat("  +Probiotic:", input$PROB_on, "  +Albumin:", input$ALB_on, "  +Flumazenil:", input$FLU_on, "\n")
  })

  output$meld_plot <- renderPlot({
    d <- data.frame(score=c("MELD", "Child-Pugh (est.)"),
                    value=c(input$MELD, 5 + input$MELD/4))
    ggplot(d, aes(score, value)) +
      geom_col(fill=c("#1976d2","#c62828")) + theme_minimal(14) +
      labs(y="Score", x="")
  })

  output$phenotype_bars <- renderPlot({
    d <- data.frame(feature=c("Hepatocyte mass","Sarcopenia","Shunt frac"),
                    value=c(input$HEPMASS,input$SARCO,input$fps))
    ggplot(d, aes(feature, value, fill=feature)) +
      geom_col() + theme_minimal(14) +
      scale_fill_manual(values=c("#2e7d32","#ef6c00","#5e35b1")) +
      labs(y="0–1") + guides(fill="none")
  })

  output$pk_plot <- renderPlot({
    d <- sim_data() %>% select(day, LOLA_C, BCAA_C, BZ_C, GPB_C) %>%
      pivot_longer(-day) %>% mutate(value = pmax(value,1e-6))
    ggplot(d, aes(day, value, color=name)) + geom_line(linewidth=1) +
      theme_minimal(13) + labs(x="Day", y="Concentration") +
      scale_y_log10()
  })

  output$pk_alb_flu <- renderPlot({
    d <- sim_data() %>% select(day, ALB_C, FLU_C) %>% pivot_longer(-day)
    ggplot(d, aes(day, value, color=name)) + geom_line(linewidth=1) +
      theme_minimal(13) + labs(x="Day", y="Concentration")
  })

  output$pk_gut <- renderPlot({
    d <- sim_data() %>% select(day, pH_lumen, urease_eff) %>% pivot_longer(-day)
    ggplot(d, aes(day, value, color=name)) + geom_line(linewidth=1) +
      facet_wrap(~name, scales="free_y") + theme_minimal(13) +
      labs(x="Day", y="Value")
  })

  output$nh3_axis <- renderPlot({
    d <- sim_data() %>% select(day, GUT_NH3, PORTAL_NH3, SYS_NH3) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value, color=name)) + geom_line(linewidth=1) +
      theme_minimal(13) + labs(x="Day", y="NH3 (µmol/L)")
  })

  output$nh3_clear <- renderPlot({
    d <- sim_data() %>% select(day, hepClear, muscleClear) %>%
      pivot_longer(-day)
    ggplot(d, aes(day, value, color=name)) + geom_line(linewidth=1) +
      theme_minimal(13) + labs(x="Day", y="Clearance rate (µmol/L/h)")
  })

  output$nh3_trial <- renderPlot({
    d <- sim_data()
    benchmark <- 75 * 0.62
    ggplot(d, aes(day, SYS_NH3)) + geom_line(linewidth=1.2, color="#1976d2") +
      geom_hline(yintercept=benchmark, lty=2, color="#c62828") +
      annotate("text", x=max(d$day)*0.7, y=benchmark+5,
               label="LOLA target (−38%, Kircheis 1997)", color="#c62828") +
      theme_minimal(13) + labs(x="Day", y="Systemic NH3 (µmol/L)")
  })

  output$brain_nh3 <- renderPlot({
    d <- sim_data() %>% select(day, BRAIN_NH3, BRAIN_GLN) %>% pivot_longer(-day)
    ggplot(d, aes(day, value, color=name)) + geom_line(linewidth=1) +
      facet_wrap(~name, scales="free_y") + theme_minimal(13)
  })

  output$brain_swell <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, SWELL)) + geom_area(fill="#ec407a", alpha=0.4) +
      geom_line(linewidth=1.2, color="#c2185b") +
      theme_minimal(13) + labs(y="Astrocyte swelling index (0-1)", x="Day")
  })

  output$brain_mn <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, Mn_BRAIN)) + geom_line(linewidth=1.2, color="#6a1b9a") +
      theme_minimal(13) + labs(y="Brain Mn (nmol/g)", x="Day")
  })

  output$infl_lps <- renderPlot({
    d <- sim_data() %>% select(day, LPS, TNFa) %>% pivot_longer(-day)
    ggplot(d, aes(day, value, color=name)) + geom_line(linewidth=1) +
      facet_wrap(~name, scales="free_y") + theme_minimal(13)
  })

  output$infl_gaba <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, GABA_PAM)) +
      geom_line(linewidth=1.2, color="#fbc02d") +
      geom_hline(yintercept=1, lty=2) + theme_minimal(13) +
      labs(y="GABA-A PAM (relative)", x="Day")
  })

  output$infl_synergy <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, infl_amp)) + geom_line(linewidth=1.2, color="#c62828") +
      theme_minimal(13) +
      labs(y="NH3-inflammation synergy factor", x="Day",
           title="Shawcross 2004: LPS + NH3 → severe HE risk")
  })

  output$wh_traj <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, WH)) + geom_line(linewidth=1.4, color="#37474f") +
      geom_hline(yintercept=c(1,2,3), lty=3) +
      annotate("text", x=max(d$day)*0.05, y=c(1.05,2.05,3.05),
               label=c("Grade I","Grade II","Grade III"), hjust=0, size=3) +
      theme_minimal(13) + labs(y="West Haven (continuous)", x="Day") +
      coord_cartesian(ylim=c(0,4))
  })

  output$hazard <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, DEATH_HZ)) + geom_line(linewidth=1.4, color="#b71c1c") +
      theme_minimal(13) + labs(y="Cumulative hazard (a.u.)", x="Day")
  })

  output$phes <- renderPlot({
    d <- sim_data() %>% mutate(PHES = -4*pmin(WH,4)/4)
    ggplot(d, aes(day, PHES)) + geom_line(linewidth=1.2, color="#0097a7") +
      geom_hline(yintercept=-4, lty=2, color="#c62828") +
      annotate("text", x=max(d$day)*0.6, y=-3.7, label="MHE threshold (-4)") +
      theme_minimal(13) + labs(y="PHES surrogate", x="Day")
  })

  output$scenario_panel <- renderPlot({
    d <- all_sc %>% mutate(day = time/24)
    p1 <- ggplot(d, aes(day, WH, color=scenario)) + geom_line(linewidth=0.9) +
      theme_minimal(12) + labs(y="West Haven", x="Day") + theme(legend.position="bottom")
    p1
  })

  output$scenario_table <- DT::renderDT({
    all_sc %>% filter(time==720) %>%
      select(scenario, SYS_NH3, BRAIN_NH3, SWELL, WH, DEATH_HZ) %>%
      mutate(across(where(is.numeric), ~round(.,2))) %>%
      DT::datatable(options=list(pageLength=10))
  })

  output$fischer <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(day, Fischer)) + geom_line(linewidth=1.2, color="#558b2f") +
      geom_hline(yintercept=3.5, lty=2) +
      annotate("text", x=max(d$day)*0.6, y=3.6, label="Normal Fischer ratio ~3.5") +
      theme_minimal(13) + labs(y="BCAA/AAA", x="Day")
  })

  output$gln_both <- renderPlot({
    d <- sim_data() %>% select(day, GLN_PL, BRAIN_GLN) %>% pivot_longer(-day)
    ggplot(d, aes(day, value, color=name)) + geom_line(linewidth=1) +
      facet_wrap(~name, scales="free_y") + theme_minimal(13)
  })

  output$trial_overlay <- renderPlot({
    d <- sim_data()
    benchmarks <- data.frame(
      trial = c("Bass 2010 (RR 0.42)","Sharma 2013 (-24% mort.)","Les 2011 (-44% HE)"),
      y = c(2.5, 2.2, 1.9)
    )
    ggplot(d, aes(day, WH)) + geom_line(linewidth=1.2) +
      geom_hline(data=benchmarks, aes(yintercept=y, color=trial), lty=2) +
      theme_minimal(13) + labs(y="West Haven", x="Day")
  })
}

shinyApp(ui, server)
