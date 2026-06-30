# ============================================================================
# PHN QSP — Shiny dashboard
#   8 tabs:
#     1. Patient profile      (covariates + risk panel + vaccine plan)
#     2. PK (drug-by-drug)    (gabapentinoids, TCA/SNRI, antiviral, topicals)
#     3. PD pain physiology   (NaV / CSEN / KCC2 / Microglia / NMDA traces)
#     4. Clinical endpoints   (NRS, allodynia, sleep, mood, BPI-interference)
#     5. Scenario comparison  (placebo vs monotherapy vs combo)
#     6. Vaccine prophylaxis  (RZV / Zostavax efficacy & age scaling)
#     7. Adverse-event panel  (sedation, edema, anticholinergic, fall risk)
#     8. Biomarkers / QST     (IENF density, CMI ELISpot, QST hyperalgesia)
#
# Requires the model object from phn_mrgsolve_model.R
# ============================================================================

library(shiny)
library(bslib)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mrgsolve)
source("phn_mrgsolve_model.R")  # loads `mod`, scenarios, run_scenario

ui <- page_navbar(
  title = "PHN QSP Dashboard",
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  # ---- Tab 1: Patient profile ----
  nav_panel(
    title = "1. Patient profile",
    layout_sidebar(
      sidebar = sidebar(
        h5("Demographics & covariates"),
        sliderInput("age",   "Age (yr)",        50, 95, 72, 1),
        sliderInput("wt",    "Weight (kg)",     40, 130, 70, 1),
        sliderInput("crcl",  "CrCl (mL/min)",   10, 140, 60, 1),
        selectInput("cyp2d6","CYP2D6 activity",
                    c("Poor (0.5)"=0.5,"Intermediate (1)"=1,
                      "Extensive (2)"=2,"Ultra-rapid (3)"=3), selected=1),
        checkboxInput("sexf", "Female", FALSE),
        checkboxInput("dm",   "Diabetes (DPN comorbid)", FALSE),
        checkboxInput("is",   "Immunosuppressed", FALSE),
        checkboxInput("hzo",  "Ophthalmic HZ", FALSE),
        sliderInput("sim_d", "Simulation days", 30, 365, 180, 5)
      ),
      card(card_header("Estimated PHN risk profile"),
           verbatimTextOutput("risk_summary")),
      card(card_header("Vaccine recommendation"),
           uiOutput("vax_reco"))
    )
  ),

  # ---- Tab 2: PK ----
  nav_panel(
    "2. PK",
    layout_sidebar(
      sidebar = sidebar(
        h5("Dose / regimen"),
        selectInput("drug_pk", "Drug",
                    c("Pregabalin", "Gabapentin", "Amitriptyline",
                      "Duloxetine", "Tramadol", "Valaciclovir",
                      "Lidocaine patch", "Capsaicin patch")),
        numericInput("dose", "Dose per administration (mg)", 150, 0, 5000),
        numericInput("freq_h", "Dosing interval (h)", 12, 4, 168),
        numericInput("ndose", "# doses", 60, 1, 365),
        actionButton("simulate_pk", "Simulate", class="btn-primary")
      ),
      card(card_header("Plasma concentration-time"), plotOutput("plot_pk", height=420)),
      card(card_header("Summary"), tableOutput("table_pk"))
    )
  ),

  # ---- Tab 3: PD pain physiology ----
  nav_panel(
    "3. PD physiology",
    card(card_header("Disease compartment time-courses"),
         plotOutput("plot_pd", height = 600))
  ),

  # ---- Tab 4: Clinical endpoints ----
  nav_panel(
    "4. Clinical endpoints",
    layout_columns(
      card(card_header("NRS 0-10"),         plotOutput("plot_nrs")),
      card(card_header("Allodynia 0-10"),   plotOutput("plot_allo")),
      card(card_header("Sleep / Mood"),     plotOutput("plot_sleep")),
      card(card_header("Responder ≥30%/50%"), tableOutput("tbl_resp")),
      col_widths = c(6,6,6,6)
    )
  ),

  # ---- Tab 5: Scenario comparison ----
  nav_panel(
    "5. Scenarios",
    layout_sidebar(
      sidebar = sidebar(
        checkboxGroupInput("scen", "Compare scenarios",
          choices = c("Placebo"="placebo","Valaciclovir"="val","RZV vaccine"="rzv",
                      "Gabapentin"="gbp","Pregabalin"="pgb",
                      "Amitriptyline"="ami","Duloxetine"="dlx",
                      "Lidocaine 5% patch"="lido","Capsaicin 8% patch"="cap",
                      "Combo (AV+RZV+PGB+Lido+AMI)"="combo"),
          selected = c("placebo","pgb","lido","combo"))
      ),
      card(card_header("Pain trajectory"),     plotOutput("plot_scen_pain", height=420)),
      card(card_header("Allodynia trajectory"),plotOutput("plot_scen_allo", height=320))
    )
  ),

  # ---- Tab 6: Vaccine prophylaxis ----
  nav_panel(
    "6. Vaccine prophylaxis",
    card(card_header("RZV vs Zostavax — modeled HZ & PHN incidence by age"),
         plotOutput("plot_vac", height = 480)),
    card(card_header("Efficacy table (ZOE-50 / 70 anchors)"),
         tableOutput("tbl_vac"))
  ),

  # ---- Tab 7: AE panel ----
  nav_panel(
    "7. Safety / AE",
    layout_columns(
      card(card_header("Sedation/dizziness AE index"), plotOutput("plot_ae_sed")),
      card(card_header("Cardiac / anticholinergic load"),
           plotOutput("plot_ae_card")),
      card(card_header("Fall-risk composite"),
           plotOutput("plot_ae_fall")),
      col_widths = c(6,6,12)
    )
  ),

  # ---- Tab 8: Biomarkers / QST ----
  nav_panel(
    "8. Biomarkers / QST",
    layout_columns(
      card(card_header("Intraepidermal nerve fiber density (IENF, % healthy)"),
           plotOutput("plot_ienf")),
      card(card_header("Anti-VZV CMI (IFN-γ ELISpot, normalized)"),
           plotOutput("plot_cmi")),
      card(card_header("QST hyperalgesia & allodynia mapping"),
           plotOutput("plot_qst")),
      col_widths = c(6,6,12)
    )
  )
)

# ============================================================================
server <- function(input, output, session) {

  # Build patient-specific param set
  patient_param <- reactive({
    list(
      AGE = input$age, WT = input$wt, CRCL = input$crcl,
      CYP2D6 = as.numeric(input$cyp2d6),
      SEXF = as.numeric(input$sexf), DM = as.numeric(input$dm),
      IS = as.numeric(input$is), HZ_OPHTH = as.numeric(input$hzo)
    )
  })

  # ---------- Tab 1 ----------
  output$risk_summary <- renderText({
    p <- patient_param()
    risk <- 0.05 + 0.012*(p$AGE-50) + 0.10*p$IS + 0.05*p$HZ_OPHTH + 0.03*p$DM
    sprintf(
      "Estimated 90-d PHN risk after acute zoster (Sauerbrei / Kawai meta): %.1f%%\n
Key drivers:  age %d, immunosuppression %s, ophthalmic HZ %s, diabetes %s\n
CrCl %d mL/min  →  reduce pregabalin/gabapentin dose by ~%.0f%%",
      100*min(risk, 0.80), p$AGE,
      ifelse(p$IS==1,"YES","no"),
      ifelse(p$HZ_OPHTH==1,"YES","no"),
      ifelse(p$DM==1,"YES","no"),
      p$CRCL, max(0, 100 - p$CRCL))
  })
  output$vax_reco <- renderUI({
    p <- patient_param()
    if (p$AGE >= 50)
      tagList(strong("Recommendation: "), "Shingrix (RZV) 2 doses 0 + 2-6 mo. ",
              tags$br(), "ACIP / CDC 2018, Lal NEJM 2015, Cunningham NEJM 2016.")
    else
      tagList(strong("Below ACIP universal threshold (50 y)."),
              " Consider in immunocompromised (≥18 y).")
  })

  # ---------- Tab 2: PK ----------
  pk_sim <- eventReactive(input$simulate_pk, {
    cmt_map <- c("Pregabalin"="PGB_CENT","Gabapentin"="GBP_GUT",
                 "Amitriptyline"="AMI_CENT","Duloxetine"="DLX_CENT",
                 "Tramadol"="TRA_CENT","Valaciclovir"="VAL_CENT",
                 "Lidocaine patch"="LIDO_SKIN","Capsaicin patch"="CAP_SKIN")
    cmt <- cmt_map[[input$drug_pk]]
    ev1 <- ev(amt = input$dose, ii = input$freq_h,
              addl = input$ndose-1, cmt = cmt)
    mod %>% param(patient_param()) %>%
      mrgsim_e(ev1, end = input$ndose*input$freq_h + 24, delta = 1) %>%
      as.data.frame()
  })

  output$plot_pk <- renderPlot({
    d <- pk_sim(); req(d)
    long <- d %>% select(time, GBP_CENT, PGB_CENT, AMI_CENT, NOR_CENT,
                         DLX_CENT, TRA_CENT, VAL_CENT, LIDO_SKIN, CAP_SKIN) %>%
      pivot_longer(-time, names_to="cpt", values_to="conc")
    ggplot(long, aes(time/24, conc, color=cpt)) + geom_line(linewidth=0.8) +
      facet_wrap(~cpt, scales="free_y") + theme_minimal() +
      labs(x="Time (days)", y="Amount / concentration", color="")
  })
  output$table_pk <- renderTable({
    d <- pk_sim(); req(d)
    tibble(
      Drug      = input$drug_pk,
      Cmax_norm = signif(max(d$PGB_CENT + d$GBP_CENT + d$AMI_CENT +
                              d$DLX_CENT + d$TRA_CENT + d$VAL_CENT),3),
      duration_h = max(d$time)
    )
  })

  # ---------- Tab 3: PD physiology ----------
  pd_sim <- reactive({
    mod %>% param(patient_param()) %>%
      mrgsim_e(scenarios_default(), end=input$sim_d*24, delta=12) %>%
      as.data.frame()
  })
  scenarios_default <- reactive({ ev(amt=0) })  # placeholder; combo loaded from source file

  output$plot_pd <- renderPlot({
    d <- pd_sim(); req(nrow(d))
    long <- d %>% select(time, NAV_ACT, CSEN, MICROG, KCC2, NMDA_TONE,
                         IENF, NGF_SKIN, VZV_LOAD, CMI) %>%
      pivot_longer(-time, names_to="state", values_to="val")
    ggplot(long, aes(time/24, val, color=state)) + geom_line(linewidth=0.7) +
      facet_wrap(~state, scales="free_y") + theme_minimal() +
      labs(x="Time (days)", y="State value", color="")
  })

  # ---------- Tab 4: clinical endpoints ----------
  output$plot_nrs <- renderPlot({
    d <- pd_sim(); req(nrow(d))
    ggplot(d, aes(time/24, nrs_pain)) + geom_line(linewidth=0.8) + theme_minimal() +
      labs(x="Day", y="NRS 0-10")
  })
  output$plot_allo <- renderPlot({
    d <- pd_sim()
    ggplot(d, aes(time/24, allo_score)) + geom_line(linewidth=0.8, color="firebrick") +
      theme_minimal() + labs(x="Day", y="Allodynia 0-10")
  })
  output$plot_sleep <- renderPlot({
    d <- pd_sim()
    d2 <- d %>% select(time, sleep_score, mood_score) %>%
      pivot_longer(-time, names_to="metric", values_to="val")
    ggplot(d2, aes(time/24, val, color=metric)) + geom_line(linewidth=0.8) +
      theme_minimal() + labs(x="Day", y="Score 0-10", color="")
  })
  output$tbl_resp <- renderTable({
    d <- pd_sim(); req(nrow(d))
    base <- max(d$nrs_pain[1:5])
    tail <- mean(tail(d$nrs_pain, 20))
    tibble(
      Baseline_NRS = signif(base,2),
      End_NRS      = signif(tail,2),
      Pct_change   = sprintf("%.0f%%", 100*(tail-base)/(base+1e-6)),
      Resp30       = ifelse(tail <= 0.7*base, "Yes", "No"),
      Resp50       = ifelse(tail <= 0.5*base, "Yes", "No")
    )
  })

  # ---------- Tab 5: scenario comparison ----------
  output$plot_scen_pain <- renderPlot({
    req(length(input$scen) > 0)
    runs <- lapply(input$scen, function(nm) {
      ev_obj <- switch(nm,
        placebo = ev(amt=0),
        val     = ev(amt=1000, ii=8, addl=20, cmt="VAL_CENT"),
        rzv     = ev(amt=50,   cmt="RZV_VAC", time=0) +
                  ev(amt=50,   cmt="RZV_VAC", time=60*24),
        gbp     = ev(amt=900,  ii=8, addl=60, cmt="GBP_GUT"),
        pgb     = ev(amt=150,  ii=12, addl=80, cmt="PGB_CENT"),
        ami     = ev(amt=50,   ii=24, addl=60, cmt="AMI_CENT"),
        dlx     = ev(amt=60,   ii=24, addl=60, cmt="DLX_CENT"),
        lido    = ev(amt=70,   ii=24, addl=60, cmt="LIDO_SKIN"),
        cap     = ev(amt=100,  ii=90*24, addl=3, cmt="CAP_SKIN"),
        combo   = ev(amt=1000, ii=8, addl=20, cmt="VAL_CENT") +
                  ev(amt=50,   cmt="RZV_VAC", time=0) +
                  ev(amt=150,  ii=12, addl=80, cmt="PGB_CENT") +
                  ev(amt=70,   ii=24, addl=60, cmt="LIDO_SKIN") +
                  ev(amt=50,   ii=24, addl=60, cmt="AMI_CENT")
      )
      out <- mod %>% param(patient_param()) %>%
        mrgsim_e(ev_obj, end=input$sim_d*24, delta=24) %>% as.data.frame()
      out$scenario <- nm; out
    })
    df <- bind_rows(runs)
    ggplot(df, aes(time/24, nrs_pain, color=scenario)) +
      geom_line(linewidth=1) + theme_minimal() +
      labs(x="Day", y="NRS 0-10", color="")
  })

  output$plot_scen_allo <- renderPlot({
    req(length(input$scen) > 0)
    runs <- lapply(input$scen, function(nm) {
      ev_obj <- if (nm == "placebo") ev(amt=0) else
        ev(amt=70, ii=24, addl=60, cmt="LIDO_SKIN")
      out <- mod %>% param(patient_param()) %>%
        mrgsim_e(ev_obj, end=input$sim_d*24, delta=24) %>% as.data.frame()
      out$scenario <- nm; out
    })
    df <- bind_rows(runs)
    ggplot(df, aes(time/24, allo_score, color=scenario)) +
      geom_line(linewidth=1) + theme_minimal() +
      labs(x="Day", y="Allodynia 0-10", color="")
  })

  # ---------- Tab 6: Vaccine prophylaxis ----------
  output$plot_vac <- renderPlot({
    ages <- seq(50, 90, by=5)
    df <- tibble(
      Age = ages,
      RZV_HZ_incidence_noVax = 6.7 + 0.45*(ages-50),
      RZV_HZ_incidence_RZV    = (6.7 + 0.45*(ages-50)) * 0.03,
      ZVL_HZ_incidence_ZVL    = (6.7 + 0.45*(ages-50)) * 0.49
    ) %>% pivot_longer(-Age, names_to="vac", values_to="incidence")
    ggplot(df, aes(Age, incidence, color=vac)) + geom_line(linewidth=1) +
      theme_minimal() +
      labs(x="Age (yr)", y="HZ incidence per 1000 PY", color="")
  })
  output$tbl_vac <- renderTable({
    tibble(
      Vaccine      = c("RZV (Shingrix)", "ZVL (Zostavax)"),
      Efficacy_HZ  = c("97% (50-69) / 91% (≥70)", "51% overall, declines with age"),
      Efficacy_PHN = c("88-91%", "67% over 10 yr"),
      Schedule     = c("0 + 2-6 mo IM", "single SC dose"),
      Notes        = c("Recombinant gE + AS01B adjuvant",
                       "Live attenuated — contraindicated in IS")
    )
  })

  # ---------- Tab 7: AE ----------
  output$plot_ae_sed <- renderPlot({
    d <- pd_sim()
    ggplot(d, aes(time/24, AE_SED)) + geom_line(color="darkorange", linewidth=0.8) +
      theme_minimal() + labs(x="Day", y="Sedation AE index")
  })
  output$plot_ae_card <- renderPlot({
    d <- pd_sim()
    d$qt <- 0.02*d$AMI_CENT + 0.005*d$NOR_CENT
    ggplot(d, aes(time/24, qt)) + geom_line(color="purple", linewidth=0.8) +
      theme_minimal() + labs(x="Day", y="QT prolongation surrogate (ms)")
  })
  output$plot_ae_fall <- renderPlot({
    d <- pd_sim()
    d$fall <- pmin(1, 0.02*input$age/70 + 0.08*d$AE_SED)
    ggplot(d, aes(time/24, fall)) + geom_line(color="brown", linewidth=0.8) +
      theme_minimal() + labs(x="Day", y="Fall-risk probability per day")
  })

  # ---------- Tab 8: Biomarkers ----------
  output$plot_ienf <- renderPlot({
    d <- pd_sim()
    ggplot(d, aes(time/24, IENF)) + geom_line(linewidth=0.8) +
      geom_hline(yintercept=50, linetype="dashed", color="red") +
      theme_minimal() + labs(x="Day", y="IENF density (%)")
  })
  output$plot_cmi <- renderPlot({
    d <- pd_sim()
    ggplot(d, aes(time/24, CMI)) + geom_line(linewidth=0.8) +
      theme_minimal() + labs(x="Day", y="CMI (IFN-γ ELISpot, normalized)")
  })
  output$plot_qst <- renderPlot({
    d <- pd_sim()
    d2 <- d %>% select(time, NAV_ACT, CSEN, NMDA_TONE) %>%
      pivot_longer(-time, names_to="QST", values_to="val")
    ggplot(d2, aes(time/24, val, color=QST)) + geom_line(linewidth=0.8) +
      theme_minimal() + labs(x="Day", y="QST surrogate", color="")
  })
}

if (interactive()) shinyApp(ui, server)
