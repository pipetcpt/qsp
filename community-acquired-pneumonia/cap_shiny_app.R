## ============================================================================
## Community-Acquired Pneumonia (CAP) QSP — Shiny Dashboard
## ----------------------------------------------------------------------------
## 8 tabs: Patient profile · Antibiotic PK & ELF exposure · PK/PD indices
##         (fT>MIC, fAUC/MIC) · Pathogen & PAMP dynamics · Host inflammation
##         & biomarkers · Gas exchange & organ failure · Scenario comparison ·
##         References
##
## The dashboard is built around the three claims the underlying model makes:
##   (1) fT>MIC and fAUC/MIC are OUTPUTS, not inputs  -> tab 3
##   (2) beta-lactam killing is a PAMP source, so the first dose transiently
##       makes inflammation worse, and a macrolide blunts that even when the
##       isolate is macrolide-resistant  -> tab 4 (toggle "ermB resistant")
##   (3) corticosteroid has two arms and its net sign depends on whether the
##       antibiotic is working  -> tab 7 (compare scenarios 6 and 7)
##
## Dependencies: shiny, shinydashboard, mrgsolve, dplyr, tidyr, ggplot2, DT
## Run with:  shiny::runApp("cap_shiny_app.R")
##            (cap_mrgsolve_model.R must sit in the same directory)
## ----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DT)
})

## ---------- Lazy model load -------------------------------------------------
get_model <- function() {
  if (!exists(".CAP_MOD", envir = .GlobalEnv) || is.null(.GlobalEnv$.CAP_MOD)) {
    assign(".CAP_MOD", mread_cache("cap_mrgsolve_model", project = "."),
           envir = .GlobalEnv)
  }
  .GlobalEnv$.CAP_MOD
}

## ---------- Regimen catalogue -----------------------------------------------
REGIMENS <- c(
  "None (untreated natural history)",
  "Amoxicillin 1 g PO q8h",
  "Ceftriaxone 2 g IV q24h",
  "Ceftriaxone + azithromycin (guideline combination)",
  "Levofloxacin 750 mg PO/IV q24h",
  "Ceftriaxone + azithromycin + hydrocortisone (severe CAP)",
  "Ceftriaxone 5-day short course"
)

n_daily <- function(n) max(as.integer(n) - 1L, 0L)

build_events <- function(regimen, t0, days_abx, days_steroid) {
  none <- ev(amt = 0, cmt = "CEFC", time = 0)
  cef  <- function(d) ev(amt = 2000, cmt = "CEFC", time = t0, ii = 24, addl = n_daily(d))
  amx  <- function(d) ev(amt = 1000, cmt = "AMXD", time = t0, ii = 8,
                         addl = n_daily(3 * d))
  azi  <- function(d) ev(amt = 500,  cmt = "AZID", time = t0, ii = 24,
                         addl = n_daily(min(d, 5)))
  lvx  <- function(d) ev(amt = 750,  cmt = "LVXD", time = t0, ii = 24, addl = n_daily(d))
  hc   <- function(d) ev(amt = 200,  cmt = "HCC",  time = t0, rate = 200 / 24,
                         ii = 24, addl = n_daily(d))

  switch(regimen,
    "None (untreated natural history)"                         = none,
    "Amoxicillin 1 g PO q8h"                                   = amx(days_abx),
    "Ceftriaxone 2 g IV q24h"                                  = cef(days_abx),
    "Ceftriaxone + azithromycin (guideline combination)"       = cef(days_abx) + azi(days_abx),
    "Levofloxacin 750 mg PO/IV q24h"                           = lvx(days_abx),
    "Ceftriaxone + azithromycin + hydrocortisone (severe CAP)" = cef(days_abx) + azi(days_abx) +
                                                                 hc(days_steroid),
    "Ceftriaxone 5-day short course"                           = cef(5),
    none)
}

## ---------- Simulation driver ------------------------------------------------
run_sim <- function(p) {
  mod <- get_model()
  evs <- build_events(p$regimen, p$t0, p$days_abx, p$days_steroid)

  pars <- list(
    AGE       = p$age,
    WT        = p$wt,
    CRCL      = p$crcl,
    IMMSUP    = p$immsup,
    SEV0      = p$sev0,
    FIO2      = p$fio2,
    VIRAL     = as.numeric(p$viral),
    VACC      = as.numeric(p$vacc),
    VASOPRESS = as.numeric(p$vaso),
    MIC_CEF   = p$mic_cef,
    MIC_AMX   = p$mic_amx,
    MIC_LVX   = p$mic_lvx,
    MIC_AZI   = if (p$erm) 64 else p$mic_azi,
    YLYS_BL   = p$ylys,
    IMAX_GCP  = p$imax_gcp,
    KPERS     = p$kpers
  )

  do.call(param, c(list(mod), pars)) |>
    mrgsim(events = evs, end = p$horizon * 24, delta = 0.25) |>
    as_tibble() |>
    mutate(day = time / 24)
}

## Run every prebuilt scenario for the comparison tab
run_all_scenarios <- function(horizon, t0) {
  base <- list(age = 68, wt = 70, crcl = 90, immsup = 0, fio2 = 0.21,
               viral = FALSE, vacc = FALSE, vaso = FALSE,
               mic_cef = 0.25, mic_amx = 0.03, mic_lvx = 1, mic_azi = 0.12,
               erm = FALSE, ylys = 1.0, imax_gcp = 0.45, kpers = 0.010,
               horizon = horizon, t0 = t0, days_abx = 7, days_steroid = 4)

  spec <- list(
    "1 Untreated"                = modifyList(base, list(regimen = REGIMENS[1], sev0 = 0.55)),
    "2 Amoxicillin (mild)"       = modifyList(base, list(regimen = REGIMENS[2], sev0 = 0.30, age = 52)),
    "3 Ceftriaxone (ward)"       = modifyList(base, list(regimen = REGIMENS[3], sev0 = 0.55)),
    "4 Ceftriaxone+azithro"      = modifyList(base, list(regimen = REGIMENS[4], sev0 = 0.55)),
    "5 Levofloxacin"             = modifyList(base, list(regimen = REGIMENS[5], sev0 = 0.55)),
    "6 Severe + steroid"         = modifyList(base, list(regimen = REGIMENS[6], sev0 = 0.85,
                                                         vaso = TRUE, fio2 = 0.60, age = 70)),
    "7 Severe + steroid, PNSP"   = modifyList(base, list(regimen = REGIMENS[6], sev0 = 0.85,
                                                         vaso = TRUE, fio2 = 0.60, age = 70,
                                                         mic_cef = 4.0, erm = TRUE)),
    "8 Ceftriaxone 5-day"        = modifyList(base, list(regimen = REGIMENS[7], sev0 = 0.55)),
    "9 Delayed 18 h (severe)"    = modifyList(base, list(regimen = REGIMENS[4], sev0 = 0.85,
                                                         vaso = TRUE, fio2 = 0.60, age = 70,
                                                         t0 = 18))
  )
  bind_rows(lapply(names(spec), function(nm) run_sim(spec[[nm]]) |> mutate(scenario = nm)))
}

THEME <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.background = element_rect(fill = "grey92"))

lineplot <- function(df, y, ylab, logy = FALSE, hline = NULL, colour = NULL) {
  g <- ggplot(df, aes(x = day, y = .data[[y]],
                      colour = if (is.null(colour)) NULL else .data[[colour]])) +
    geom_line(linewidth = 0.8) +
    labs(x = "Day", y = ylab) + THEME
  if (!is.null(hline)) g <- g + geom_hline(yintercept = hline, linetype = 2, colour = "grey40")
  if (logy) g <- g + scale_y_log10()
  g
}

## ============================== UI ==========================================
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "CAP QSP Dashboard", titleWidth = 320),
  dashboardSidebar(
    width = 320,
    sidebarMenu(
      menuItem("1 · Patient profile",        tabName = "profile",  icon = icon("user")),
      menuItem("2 · Antibiotic PK / ELF",    tabName = "pk",       icon = icon("prescription-bottle")),
      menuItem("3 · PK/PD indices",          tabName = "pkpd",     icon = icon("bullseye")),
      menuItem("4 · Pathogen & PAMP",        tabName = "bug",      icon = icon("bacterium")),
      menuItem("5 · Inflammation & markers", tabName = "inflam",   icon = icon("fire")),
      menuItem("6 · Gas exchange & organs",  tabName = "organ",    icon = icon("lungs")),
      menuItem("7 · Scenario comparison",    tabName = "compare",  icon = icon("layer-group")),
      menuItem("8 · References",             tabName = "refs",     icon = icon("book"))
    ),
    hr(),
    selectInput("regimen", "Regimen", choices = REGIMENS, selected = REGIMENS[4]),
    sliderInput("sev0",   "Presenting severity (0-1)", 0.1, 1.0, 0.55, 0.05),
    sliderInput("t0",     "Door-to-antibiotic time (h)", 0, 36, 6, 1),
    sliderInput("days_abx", "Antibiotic duration (days)", 3, 14, 7, 1),
    sliderInput("days_steroid", "Hydrocortisone duration (days)", 0, 8, 4, 1),
    sliderInput("horizon", "Simulation horizon (days)", 5, 28, 14, 1),
    actionButton("go", "Run simulation", icon = icon("play"),
                 class = "btn-danger", width = "90%")
  ),
  dashboardBody(
    tabItems(
      ## ---- 1 patient profile -------------------------------------------------
      tabItem(
        "profile",
        fluidRow(
          box(title = "Host", width = 4, status = "primary", solidHeader = TRUE,
              sliderInput("age", "Age (years)", 18, 95, 68, 1),
              sliderInput("wt", "Weight (kg)", 40, 130, 70, 1),
              sliderInput("crcl", "Creatinine clearance (mL/min)", 15, 150, 90, 5),
              sliderInput("immsup", "Immunosuppression (0-1)", 0, 1, 0, 0.1),
              checkboxInput("vacc", "PCV20 / PPSV23 vaccinated", FALSE),
              checkboxInput("viral", "Antecedent influenza co-infection", FALSE)),
          box(title = "Isolate susceptibility (MIC, mg/L)", width = 4,
              status = "warning", solidHeader = TRUE,
              sliderInput("mic_cef", "Ceftriaxone", 0.03, 8, 0.25, 0.01),
              sliderInput("mic_amx", "Amoxicillin / penicillin", 0.01, 8, 0.03, 0.01),
              sliderInput("mic_lvx", "Levofloxacin", 0.25, 16, 1, 0.25),
              sliderInput("mic_azi", "Azithromycin", 0.03, 8, 0.12, 0.03),
              checkboxInput("erm", "ermB macrolide-resistant (MIC forced to 64)", FALSE),
              helpText("Tick ermB and keep the macrolide on: the KILL channel dies, ",
                       "the immunomodulatory channel does not.")),
          box(title = "Support & mechanism switches", width = 4,
              status = "danger", solidHeader = TRUE,
              sliderInput("fio2", "FiO2", 0.21, 1.0, 0.21, 0.01),
              checkboxInput("vaso", "Vasopressor / fluid support", FALSE),
              sliderInput("ylys", "beta-lactam lysis yield YLYS_BL", 0, 1.2, 1.0, 0.05),
              sliderInput("imax_gcp", "Steroid immunosuppressive arm IMAX_GCP", 0, 0.8, 0.45, 0.05),
              sliderInput("kpers", "Persister switch rate KPERS (1/h)", 0, 0.05, 0.010, 0.002),
              helpText("Set YLYS_BL to 0 and the post-first-dose inflammatory surge ",
                       "disappears. Set IMAX_GCP to 0 and steroid can no longer harm."))
        ),
        fluidRow(
          valueBoxOutput("vb_curb", width = 3),
          valueBoxOutput("vb_sofa", width = 3),
          valueBoxOutput("vb_stable", width = 3),
          valueBoxOutput("vb_mort", width = 3)
        ),
        fluidRow(
          box(title = "Vital-sign trajectory", width = 12, status = "primary",
              solidHeader = TRUE, plotOutput("p_vitals", height = 340))
        )
      ),

      ## ---- 2 PK --------------------------------------------------------------
      tabItem(
        "pk",
        fluidRow(
          box(title = "Plasma concentration", width = 6, status = "primary",
              solidHeader = TRUE, plotOutput("p_pk_plasma", height = 320)),
          box(title = "Free epithelial-lining-fluid concentration (the site that matters)",
              width = 6, status = "primary", solidHeader = TRUE,
              plotOutput("p_pk_elf", height = 320))
        ),
        fluidRow(
          box(title = "Azithromycin: why a 0.4 mg/L serum drug works in the lung",
              width = 6, status = "info", solidHeader = TRUE,
              plotOutput("p_azi", height = 300)),
          box(title = "Hydrocortisone plasma vs effect site, and the glucose cost",
              width = 6, status = "info", solidHeader = TRUE,
              plotOutput("p_hc", height = 300))
        )
      ),

      ## ---- 3 PK/PD -----------------------------------------------------------
      tabItem(
        "pkpd",
        fluidRow(
          box(width = 12, status = "warning", solidHeader = TRUE,
              title = "These two indices are model OUTPUTS, not inputs",
              HTML("<p>Nothing in the model is told that beta-lactams are ",
                   "'time-dependent' and quinolones 'concentration-dependent'. ",
                   "Both drugs use the same saturating Emax function of free ELF ",
                   "concentration with EC50 proportional to the isolate's MIC. ",
                   "The only difference is the Hill coefficient (1.5 vs 2.6). ",
                   "fT&gt;MIC and fAUC/MIC below are time-integrals of those ",
                   "curves, accumulated in the TAM and AUCF compartments.</p>"))
        ),
        fluidRow(
          box(title = "fT>MIC (% of elapsed time), active beta-lactam",
              width = 6, status = "primary", solidHeader = TRUE,
              plotOutput("p_ftmic", height = 320)),
          box(title = "Running free ELF AUC / MIC, levofloxacin",
              width = 6, status = "primary", solidHeader = TRUE,
              plotOutput("p_aucmic", height = 320))
        ),
        fluidRow(
          box(title = "Instantaneous kill pressure vs MIC multiples",
              width = 12, status = "info", solidHeader = TRUE,
              plotOutput("p_killcurve", height = 320),
              helpText("Free ELF concentration expressed in MIC multiples, log scale. ",
                       "Where the trace sits relative to 1x MIC is the whole of PK/PD."))
        )
      ),

      ## ---- 4 pathogen --------------------------------------------------------
      tabItem(
        "bug",
        fluidRow(
          box(title = "Bacterial burden: alveolar, persister, blood",
              width = 7, status = "primary", solidHeader = TRUE,
              plotOutput("p_bugs", height = 360)),
          box(title = "Free PAMP pool — the lysis-coupled term",
              width = 5, status = "danger", solidHeader = TRUE,
              plotOutput("p_pamp", height = 360),
              helpText("Watch the first 24 h after the first beta-lactam dose. ",
                       "Killing liberates cell wall and pneumolysin; the PAMP ",
                       "pool spikes before it falls. Slide YLYS_BL to 0 to remove it."))
        ),
        fluidRow(
          box(title = "Persisters (BP) are what decide relapse after a short course",
              width = 12, status = "warning", solidHeader = TRUE,
              plotOutput("p_pers", height = 300))
        )
      ),

      ## ---- 5 inflammation ----------------------------------------------------
      tabItem(
        "inflam",
        fluidRow(
          box(title = "Cytokine network", width = 7, status = "danger",
              solidHeader = TRUE, plotOutput("p_cyt", height = 360)),
          box(title = "Cellular effectors", width = 5, status = "primary",
              solidHeader = TRUE, plotOutput("p_cells", height = 360))
        ),
        fluidRow(
          box(title = "CRP vs procalcitonin — different drivers, different half-lives",
              width = 7, status = "info", solidHeader = TRUE,
              plotOutput("p_bio", height = 320),
              helpText("PCT is driven by the bacterial PAMP pool and suppressed by the ",
                       "antiviral interferon programme; CRP is driven by IL-6. Tick ",
                       "'influenza co-infection' and the two separate.")),
          box(title = "Resolution arm (SPMs) and lymphocyte count", width = 5,
              status = "success", solidHeader = TRUE, plotOutput("p_res", height = 320))
        )
      ),

      ## ---- 6 organ -----------------------------------------------------------
      tabItem(
        "organ",
        fluidRow(
          box(title = "Alveolar oedema, barrier permeability, surfactant",
              width = 6, status = "primary", solidHeader = TRUE,
              plotOutput("p_lung", height = 320)),
          box(title = "Shunt fraction and PaO2/FiO2", width = 6, status = "primary",
              solidHeader = TRUE, plotOutput("p_gas", height = 320))
        ),
        fluidRow(
          box(title = "Haemodynamics, lactate, renal injury", width = 6,
              status = "danger", solidHeader = TRUE, plotOutput("p_sys", height = 320)),
          box(title = "SOFA and cumulative mortality hazard", width = 6,
              status = "danger", solidHeader = TRUE, plotOutput("p_sofa", height = 320))
        )
      ),

      ## ---- 7 comparison ------------------------------------------------------
      tabItem(
        "compare",
        fluidRow(
          box(width = 12, status = "warning", solidHeader = TRUE,
              title = "Nine prebuilt scenarios on one parameter set",
              HTML("<p>Scenarios <b>6</b> and <b>7</b> are the point of the model: ",
                   "identical hydrocortisone regimen, identical host, the only ",
                   "difference is the isolate's MIC. Where the antibiotic sterilises, ",
                   "steroid is protective; where it does not, the same steroid is ",
                   "harmful, because the arm that suppresses opsonophagocytosis is no ",
                   "longer redundant.</p>"),
              actionButton("go_all", "Run all nine scenarios", icon = icon("play"),
                           class = "btn-warning"))
        ),
        fluidRow(
          box(title = "Bacterial burden", width = 6, status = "primary",
              solidHeader = TRUE, plotOutput("c_bug", height = 320)),
          box(title = "PaO2/FiO2", width = 6, status = "primary",
              solidHeader = TRUE, plotOutput("c_pf", height = 320))
        ),
        fluidRow(
          box(title = "Cumulative mortality (%)", width = 6, status = "danger",
              solidHeader = TRUE, plotOutput("c_mort", height = 320)),
          box(title = "PAMP pool", width = 6, status = "danger",
              solidHeader = TRUE, plotOutput("c_pamp", height = 320))
        ),
        fluidRow(
          box(title = "Endpoint summary", width = 12, status = "info",
              solidHeader = TRUE, DTOutput("c_tab"))
        )
      ),

      ## ---- 8 references ------------------------------------------------------
      tabItem(
        "refs",
        fluidRow(
          box(width = 12, status = "primary", solidHeader = TRUE,
              title = "Primary calibration anchors",
              HTML("
<ul>
<li><b>Dequin 2023 (CAPE-COD, NEJM)</b> — hydrocortisone 200 mg/d in severe CAP,
    28-day mortality 6.2% vs 11.9%.</li>
<li><b>Torres 2015 (JAMA)</b> / <b>Blum 2015 (Lancet)</b> — corticosteroid reduces
    treatment failure and shortens time to clinical stability.</li>
<li><b>Halm 1998 (JAMA)</b> — median time to clinical stability ~3 days; the Halm
    criteria are implemented verbatim as <code>Clinically_stable</code>.</li>
<li><b>Craig 1998</b> / <b>Drusano 2004</b> — the PK/PD index framework that this
    model deliberately derives rather than assumes.</li>
<li><b>Nau &amp; Eiffert 2002 (Clin Microbiol Rev)</b> and <b>Spreer 2003 (AAC)</b> —
    antibiotic-class-dependent liberation of pro-inflammatory bacterial
    compounds; non-bacteriolytic agents release far less pneumolysin.</li>
<li><b>Martin-Loeches 2010</b> / <b>Sligl 2014</b> — macrolide combination and
    survival in severe CAP.</li>
<li><b>Uranga 2016 (JAMA Intern Med)</b> / <b>Dinh 2021 (Lancet)</b> — short-course
    therapy.</li>
<li><b>Schuetz 2018 (Lancet Infect Dis)</b> — procalcitonin-guided duration.</li>
<li><b>Kumar 2006 (Crit Care Med)</b> — time-to-effective-antimicrobial geometry.</li>
</ul>
<p>The complete annotated bibliography with PubMed links is in
<code>cap_references.md</code>.</p>")),
          box(width = 12, status = "danger", solidHeader = TRUE, title = "Disclaimer",
              HTML("<p>Educational and research use only. Parameters are literature-
                   anchored approximations, not fitted to patient-level data. This
                   dashboard must not be used for clinical decisions, prescribing, or
                   regulatory submission.</p>"))
        )
      )
    )
  )
)

## ============================== SERVER ======================================
server <- function(input, output, session) {

  params <- reactive({
    list(regimen = input$regimen, sev0 = input$sev0, t0 = input$t0,
         days_abx = input$days_abx, days_steroid = input$days_steroid,
         horizon = input$horizon, age = input$age, wt = input$wt,
         crcl = input$crcl, immsup = input$immsup, fio2 = input$fio2,
         viral = input$viral, vacc = input$vacc, vaso = input$vaso,
         mic_cef = input$mic_cef, mic_amx = input$mic_amx,
         mic_lvx = input$mic_lvx, mic_azi = input$mic_azi, erm = input$erm,
         ylys = input$ylys, imax_gcp = input$imax_gcp, kpers = input$kpers)
  })

  sim <- eventReactive(input$go, ignoreNULL = FALSE, { run_sim(params()) })
  allsim <- eventReactive(input$go_all, { run_all_scenarios(input$horizon, input$t0) })

  long <- function(df, cols) {
    df |> select(day, all_of(cols)) |> pivot_longer(-day, names_to = "series")
  }

  ## ---- value boxes ---------------------------------------------------------
  output$vb_curb <- renderValueBox({
    v <- max(sim()$CURB65)
    valueBox(v, "Peak CURB-65", icon = icon("triangle-exclamation"),
             color = if (v >= 3) "red" else if (v >= 2) "yellow" else "green")
  })
  output$vb_sofa <- renderValueBox({
    v <- max(sim()$SOFA_score)
    valueBox(v, "Peak SOFA", icon = icon("heart-pulse"),
             color = if (v >= 6) "red" else if (v >= 3) "yellow" else "green")
  })
  output$vb_stable <- renderValueBox({
    d <- sim() |> filter(Clinically_stable == 1, day > 0.5)
    v <- if (nrow(d)) sprintf("%.1f d", min(d$day)) else "not reached"
    valueBox(v, "Time to clinical stability (Halm)", icon = icon("clock"), color = "aqua")
  })
  output$vb_mort <- renderValueBox({
    v <- sprintf("%.1f%%", max(sim()$Mortality_pct))
    valueBox(v, "Cumulative mortality", icon = icon("skull"), color = "purple")
  })

  ## ---- tab 1 ---------------------------------------------------------------
  output$p_vitals <- renderPlot({
    long(sim(), c("Temperature", "Heart_rate", "Resp_rate", "MAP_mmHg", "SpO2_pct")) |>
      ggplot(aes(day, value)) + geom_line(linewidth = 0.8, colour = "#b03030") +
      facet_wrap(~series, scales = "free_y", ncol = 3) +
      labs(x = "Day", y = NULL) + THEME
  })

  ## ---- tab 2 ---------------------------------------------------------------
  lineplot_df <- function(d) {
    ggplot(d, aes(day, value, colour = series)) + geom_line(linewidth = 0.8) +
      labs(x = "Day", y = "Concentration (mg/L)") + THEME
  }
  output$p_pk_plasma <- renderPlot({
    long(sim(), c("Ceftriaxone_plasma", "Hydrocortisone_plasma")) |>
      lineplot_df()
  })
  output$p_pk_elf <- renderPlot({
    long(sim(), c("Ceftriaxone_ELF_free", "Amoxicillin_ELF_free",
                  "Azithromycin_ELF", "Levofloxacin_ELF_free")) |>
      ggplot(aes(day, value, colour = series)) + geom_line(linewidth = 0.8) +
      scale_y_log10() +
      geom_hline(yintercept = input$mic_cef, linetype = 2, colour = "grey40") +
      labs(x = "Day", y = "Free ELF concentration (mg/L, log)",
           caption = "dashed line = ceftriaxone MIC") + THEME
  })
  output$p_azi <- renderPlot({
    long(sim(), c("Azithromycin_ELF")) |>
      ggplot(aes(day, value)) + geom_line(linewidth = 0.9, colour = "#2060a0") +
      labs(x = "Day", y = "Azithromycin ELF (mg/L)") + THEME
  })
  output$p_hc <- renderPlot({
    long(sim(), c("Hydrocortisone_plasma", "Glucose")) |>
      ggplot(aes(day, value)) + geom_line(linewidth = 0.9, colour = "#207040") +
      facet_wrap(~series, scales = "free_y") + labs(x = "Day", y = NULL) + THEME
  })

  ## ---- tab 3 ---------------------------------------------------------------
  output$p_ftmic <- renderPlot({
    lineplot(sim(), "fT_above_MIC_pct", "fT>MIC (% of elapsed time)", hline = c(40, 100))
  })
  output$p_aucmic <- renderPlot({
    lineplot(sim(), "fAUC_over_MIC", "cumulative fAUC/MIC", hline = 30)
  })
  output$p_killcurve <- renderPlot({
    d <- sim() |>
      transmute(day,
                Ceftriaxone = Ceftriaxone_ELF_free / input$mic_cef,
                Amoxicillin = Amoxicillin_ELF_free / input$mic_amx,
                Levofloxacin = Levofloxacin_ELF_free / input$mic_lvx,
                Azithromycin = Azithromycin_ELF / (if (input$erm) 64 else input$mic_azi)) |>
      pivot_longer(-day, names_to = "series") |>
      filter(value > 1e-4)
    ggplot(d, aes(day, value, colour = series)) + geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 1, linetype = 2) + scale_y_log10() +
      labs(x = "Day", y = "Free ELF concentration / MIC (log)") + THEME
  })

  ## ---- tab 4 ---------------------------------------------------------------
  output$p_bugs <- renderPlot({
    long(sim(), c("Bacteria_log10", "Persisters_log10", "Bacteraemia_log10")) |>
      ggplot(aes(day, value, colour = series)) + geom_line(linewidth = 0.8) +
      labs(x = "Day", y = "log10 CFU/mL") + THEME
  })
  output$p_pamp <- renderPlot({
    lineplot(sim(), "PAMP_pool", "Free PAMP pool (au)")
  })
  output$p_pers <- renderPlot({
    lineplot(sim(), "Persisters_log10", "Persister burden (log10 CFU/mL)")
  })

  ## ---- tab 5 ---------------------------------------------------------------
  output$p_cyt <- renderPlot({
    long(sim(), c("TNF_alpha", "IL1_beta", "IL6_level", "IL8_level", "IL10_level")) |>
      ggplot(aes(day, value, colour = series)) + geom_line(linewidth = 0.8) +
      labs(x = "Day", y = "pg/mL") + THEME
  })
  output$p_cells <- renderPlot({
    long(sim(), c("Alveolar_macrophage", "Lung_neutrophils", "Blood_neutrophils")) |>
      ggplot(aes(day, value)) + geom_line(linewidth = 0.8, colour = "#404080") +
      facet_wrap(~series, scales = "free_y", ncol = 1) + labs(x = "Day", y = NULL) + THEME
  })
  output$p_bio <- renderPlot({
    long(sim(), c("CRP_level", "PCT_level")) |>
      ggplot(aes(day, value)) + geom_line(linewidth = 0.9, colour = "#806020") +
      facet_wrap(~series, scales = "free_y", ncol = 1) + labs(x = "Day", y = NULL) + THEME
  })
  output$p_res <- renderPlot({
    long(sim(), c("SPM_level", "Lymphocytes")) |>
      ggplot(aes(day, value)) + geom_line(linewidth = 0.9, colour = "#207040") +
      facet_wrap(~series, scales = "free_y", ncol = 1) + labs(x = "Day", y = NULL) + THEME
  })

  ## ---- tab 6 ---------------------------------------------------------------
  output$p_lung <- renderPlot({
    long(sim(), c("Oedema_fraction", "Permeability", "Surfactant_function",
                  "Fibroproliferation")) |>
      ggplot(aes(day, value)) + geom_line(linewidth = 0.8, colour = "#206070") +
      facet_wrap(~series, scales = "free_y") + labs(x = "Day", y = NULL) + THEME
  })
  output$p_gas <- renderPlot({
    long(sim(), c("Shunt_fraction", "PaO2_FiO2", "SpO2_pct")) |>
      ggplot(aes(day, value)) + geom_line(linewidth = 0.8, colour = "#2060a0") +
      facet_wrap(~series, scales = "free_y", ncol = 1) + labs(x = "Day", y = NULL) + THEME
  })
  output$p_sys <- renderPlot({
    long(sim(), c("MAP_mmHg", "Lactate", "Renal_injury")) |>
      ggplot(aes(day, value)) + geom_line(linewidth = 0.8, colour = "#a03030") +
      facet_wrap(~series, scales = "free_y", ncol = 1) + labs(x = "Day", y = NULL) + THEME
  })
  output$p_sofa <- renderPlot({
    long(sim(), c("SOFA_score", "Mortality_pct")) |>
      ggplot(aes(day, value)) + geom_line(linewidth = 0.9, colour = "#602080") +
      facet_wrap(~series, scales = "free_y", ncol = 1) + labs(x = "Day", y = NULL) + THEME
  })

  ## ---- tab 7 ---------------------------------------------------------------
  output$c_bug  <- renderPlot(lineplot(allsim(), "Bacteria_log10", "log10 CFU/mL",
                                       colour = "scenario"))
  output$c_pf   <- renderPlot(lineplot(allsim(), "PaO2_FiO2", "PaO2/FiO2",
                                       hline = c(100, 300), colour = "scenario"))
  output$c_mort <- renderPlot(lineplot(allsim(), "Mortality_pct", "Cumulative mortality (%)",
                                       colour = "scenario"))
  output$c_pamp <- renderPlot(lineplot(allsim(), "PAMP_pool", "PAMP pool (au)",
                                       colour = "scenario"))
  output$c_tab  <- renderDT({
    allsim() |>
      group_by(scenario) |>
      summarise(
        `Peak SOFA`            = max(SOFA_score),
        `Nadir PaO2/FiO2`      = round(min(PaO2_FiO2)),
        `Peak PAMP`            = round(max(PAMP_pool), 2),
        `Peak CRP (mg/L)`      = round(max(CRP_level)),
        `Peak PCT (ng/mL)`     = round(max(PCT_level), 2),
        `Day bacteria <1e2`    = {
          i <- which(Bacteria_log10 < 2 & day > 0.5)
          if (length(i)) round(day[min(i)], 2) else NA_real_
        },
        `Day clinically stable` = {
          i <- which(Clinically_stable == 1 & day > 0.5)
          if (length(i)) round(day[min(i)], 2) else NA_real_
        },
        `Mortality (%)`        = round(max(Mortality_pct), 2),
        .groups = "drop") |>
      datatable(options = list(dom = "t", pageLength = 12), rownames = FALSE)
  })
}

shinyApp(ui, server)
