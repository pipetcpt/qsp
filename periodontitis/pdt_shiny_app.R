##  Periodontitis QSP — Shiny Dashboard
##  ============================================================================
##  Interactive front end for pdt_mrgsolve_model.R.
##
##  Run with:
##      shiny::runApp("pdt_shiny_app.R")
##  or from this directory:
##      R -e 'shiny::runApp("pdt_shiny_app.R", launch.browser = TRUE)'
##
##  The dashboard is organised around the three claims the model is built on,
##  rather than around the compartment list:
##
##    * the disease is a BISTABLE loop  -> Tab 8 (Separatrix) lets you find the
##      tipping point for the current patient and see whether the instrument
##      can actually cross it;
##    * neutrophils are recruited and DISARMED -> Tab 4 plots the killing
##      capacity as a state variable and Tab 9 lets you switch the subversion
##      off entirely;
##    * bone is a RATCHET -> Tab 6 separates what recovers from what does not,
##      and Tab 10 shows what a year of delay costs.
##
##  Nine of the ten tabs are simulation views; the last is the site-population
##  view that explains why individual sites are bimodal and trial means are not.
##  ============================================================================

library(shiny)
library(mrgsolve)

source("pdt_mrgsolve_model.R")
MOD <- pdt_model()

PAL <- c(base = "#2c6fbb", drug = "#c0392b", alt = "#27ae60",
         warn = "#e67e22", grey = "#7f8c8d", purple = "#8e44ad")

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Periodontitis QSP Dashboard — dysbiosis, neutrophil subversion, and the bone ratchet"),
  tags$p(style = "color:#666;margin-top:-8px",
         "Research / teaching model. Not for clinical use."),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Site"),
      sliderInput("bloss0", "Baseline alveolar bone loss (mm)",
                  min = 0.1, max = 8.5, value = 6.0, step = 0.1),
      helpText("Probing depth is not set directly — it is computed from the ",
               "attachment level, the gingival margin and the probe-penetration ",
               "artefact, so it emerges from the site's state."),

      h4("Patient"),
      checkboxInput("smoke", "Heavy smoker", FALSE),
      sliderInput("hbset", "Intrinsic HbA1c without periodontal contribution (%)",
                  min = 5.0, max = 11.0, value = 5.4, step = 0.1),
      sliderInput("ohi", "Supragingival plaque control (0 = none, 1 = excellent)",
                  min = 0, max = 1, value = 0.6, step = 0.05),
      checkboxInput("ladi", "Leukocyte adhesion deficiency-I", FALSE),
      checkboxInput("aa", "A. actinomycetemcomitans co-infection", FALSE),

      h4("Mechanical therapy"),
      radioButtons("mech", NULL,
                   c("None" = "none",
                     "Scaling and root planing" = "srp",
                     "Open flap surgery" = "surg",
                     "Surgery + regeneration (EMD/GTR)" = "regen"),
                   selected = "srp"),
      checkboxInput("spt", "3-monthly supportive periodontal therapy", TRUE),

      h4("Adjuncts"),
      checkboxInput("am",   "Amoxicillin 500 + metronidazole 400 TID x 7 d", FALSE),
      checkboxInput("sdd",  "Sub-antimicrobial doxycycline 20 mg BID x 9 mo", FALSE),
      checkboxInput("mino", "Local minocycline microspheres (0, 3, 6 mo)", FALSE),
      checkboxInput("chx",  "Chlorhexidine 0.12% rinse x 60 d", FALSE),
      checkboxInput("amy",  "AMY-101 (C3 inhibitor) local, weekly x 3", FALSE),
      checkboxInput("rve",  "Resolvin E1 local, weekly x 12", FALSE),
      checkboxInput("atnf", "Anti-TNF 5 mg/kg q8w", FALSE),
      checkboxInput("dmab", "Denosumab 60 mg SC q6mo", FALSE),
      checkboxInput("nsaid","Systemic NSAID (COX inhibition)", FALSE),
      sliderInput("glyc", "Glycaemic therapy effect on HbA1c set point (%)",
                  min = 0, max = 3, value = 0, step = 0.1),

      h4("Mechanistic switches"),
      sliderInput("esub", "C5aR1 x TLR2 subversion strength (E_SUBVERT)",
                  min = 0, max = 6, value = 4, step = 0.25),
      helpText("Set this to zero and the same bacteria meet a competent ",
               "neutrophil. It is the single most consequential number in the model."),

      h4("Horizon"),
      sliderInput("end", "Simulation length (days)", min = 180, max = 3650,
                  value = 730, step = 30),
      actionButton("go", "Simulate", class = "btn-primary", width = "100%")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1. Patient profile",
                 br(),
                 fluidRow(column(6, tableOutput("profile")),
                          column(6, tableOutput("staging"))),
                 h4("Where this site sits relative to its own tipping point"),
                 verbatimTextOutput("verdict")),

        tabPanel("2. Drug PK",
                 br(), plotOutput("pk_plasma", height = "300px"),
                 plotOutput("pk_local", height = "300px")),

        tabPanel("3. Microbiology",
                 br(), plotOutput("micro", height = "340px"),
                 plotOutput("keystone", height = "300px")),

        tabPanel("4. Innate immunity and the subversion",
                 br(), plotOutput("innate", height = "340px"),
                 plotOutput("clearance", height = "300px")),

        tabPanel("5. Mediators and proteases",
                 br(), plotOutput("mediators", height = "340px"),
                 plotOutput("mmp", height = "300px")),

        tabPanel("6. Clinical endpoints",
                 br(), plotOutput("clinical", height = "340px"),
                 h4("What the probing-depth reduction is actually made of"),
                 tableOutput("decomp"),
                 plotOutput("ratchet", height = "300px")),

        tabPanel("7. Systemic spillover",
                 br(), plotOutput("systemic", height = "340px"),
                 plotOutput("pisa", height = "300px")),

        tabPanel("8. Separatrix",
                 br(),
                 helpText("Bisection on the residual dysbiotic biomass ",
                          "immediately after debridement: how far down must the ",
                          "biofilm be pushed for this site not to relapse?"),
                 actionButton("sep_go", "Locate the separatrix for this site"),
                 verbatimTextOutput("sep")),

        tabPanel("9. Scenario comparison",
                 br(),
                 helpText("Eleven prebuilt arms at a deep site. Adjuncts separate ",
                          "only where closed instrumentation cannot cross the ",
                          "separatrix on its own."),
                 actionButton("scn_go", "Run all scenarios"),
                 plotOutput("scenarios", height = "420px"),
                 tableOutput("scn_tab")),

        tabPanel("10. Site population and critical depth",
                 br(),
                 helpText("Individual sites are bimodal; the clinical mean is an ",
                          "average across the separatrix. This tab reproduces the ",
                          "depth bands trials report, and the critical probing ",
                          "depth of Lindhe 1982 falls out of it."),
                 actionButton("pop_go", "Run the site population (slow, ~2 min)"),
                 tableOutput("bands"),
                 plotOutput("critdepth", height = "360px")),

        tabPanel("11. Calibration",
                 br(),
                 helpText("Every target is from the literature. None was used as ",
                          "an objective function."),
                 actionButton("cal_go", "Run calibration report (slow, ~3 min)"),
                 tableOutput("cal"))
      )
    )
  )
)

## ---------------------------------------------------------------------------
## SERVER
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  pars <- reactive({
    list(SMOKE     = as.numeric(input$smoke),
         HB_SET    = input$hbset,
         GLYC_RX   = input$glyc,
         OHI       = input$ohi,
         LADI      = as.numeric(input$ladi),
         AA_ON     = as.numeric(input$aa),
         NSAID_ON  = as.numeric(input$nsaid),
         E_SUBVERT = input$esub)
  })

  patient <- reactive({
    PDT_patient(MOD, bloss0 = input$bloss0, param = pars())
  })

  events <- reactive({
    e <- NULL
    add <- function(x) if (is.null(e)) x else c(e, x)
    if (input$mech == "srp")  e <- add(ev_srp())
    if (input$mech %in% c("surg", "regen")) e <- add(ev_surgery())
    if (input$mech == "regen") e <- add(ev_regen())
    if (input$spt && input$mech != "none")
      e <- add(ev_spt(from = 90, to = max(90, input$end)))
    if (input$am)   e <- add(ev_amox_metro())
    if (input$sdd)  e <- add(ev_sdd())
    if (input$mino) e <- add(ev_minocycline())
    if (input$chx)  e <- add(ev_chx())
    if (input$amy)  e <- add(ev_amy101())
    if (input$rve)  e <- add(ev_rve1())
    if (input$atnf) e <- add(ev_antitnf(n = max(1, ceiling(input$end/56))))
    if (input$dmab) e <- add(ev_denosumab(n = max(1, ceiling(input$end/182))))
    e
  })

  sim <- eventReactive(input$go, {
    pt <- param(patient(), pars())
    e  <- events()
    if (is.null(e)) mrgsim_df(pt, end = input$end, delta = 1)
    else            mrgsim_df(pt, events = e, end = input$end, delta = 1)
  }, ignoreNULL = FALSE)

  ## ---- helpers -------------------------------------------------------------
  ln <- function(d, vars, labs, cols, ylab, main, ylim = NULL, log = "") {
    yl <- if (is.null(ylim)) range(unlist(d[vars]), na.rm = TRUE) else ylim
    plot(d$time, d[[vars[1]]], type = "n", ylim = yl, log = log,
         xlab = "Time (days)", ylab = ylab, main = main, las = 1, bty = "l")
    grid(col = "#eeeeee")
    for (i in seq_along(vars))
      lines(d$time, d[[vars[i]]], col = cols[i], lwd = 2)
    legend("topright", labs, col = cols, lwd = 2, bty = "n", cex = 0.85)
  }

  ## ---- 1. profile ----------------------------------------------------------
  output$profile <- renderTable({
    d <- sim(); i0 <- 1; i1 <- nrow(d)
    data.frame(
      Variable = c("Probing depth (mm)", "Clinical attachment level (mm)",
                   "Alveolar bone loss (mm)", "Gingival margin (mm)",
                   "Bleeding on probing (%)", "Inflammatory tone (health = 1)",
                   "Neutrophil killing capacity (% of health)",
                   "Dysbiotic biomass", "P. gingivalis (% of biofilm)",
                   "GCF flow (uL/min)"),
      Baseline = round(c(d$PPD_o[i0], d$CAL_m[i0], d$BLOSS[i0], d$GM[i0],
                         d$BOPs_o[i0], d$INFLAM_o[i0], 100*d$PKN_o[i0],
                         d$BIOD[i0], 100*d$PG[i0], d$GCF_ul[i0]), 2),
      End      = round(c(d$PPD_o[i1], d$CAL_m[i1], d$BLOSS[i1], d$GM[i1],
                         d$BOPs_o[i1], d$INFLAM_o[i1], 100*d$PKN_o[i1],
                         d$BIOD[i1], 100*d$PG[i1], d$GCF_ul[i1]), 2))
  })

  output$staging <- renderTable({
    d <- sim(); i1 <- nrow(d)
    stg <- function(b) if (b < 1) "Stage I" else if (b < 2) "Stage II" else
                       if (b < 5) "Stage III" else "Stage IV"
    data.frame(
      Item = c("2018 EFP/AAP stage (bone loss)", "Tooth survival (%)",
               "Cumulative MRONJ risk (%)", "Periodontal inflamed surface area (mm2)",
               "CRP (mg/L)", "HbA1c (%)", "Flow-mediated dilation (%)",
               "Achievable closed biofilm disruption (%)",
               "Achievable closed calculus removal (%)"),
      Value = c(stg(d$BLOSS[i1]),
                sprintf("%.1f", d$TSURV_o[i1]), sprintf("%.3f", d$MRONJ_o[i1]),
                sprintf("%.0f", d$PISA_o[i1]),  sprintf("%.2f", d$CRP[i1]),
                sprintf("%.2f", d$HBA1C[i1]),   sprintf("%.2f", d$FMD[i1]),
                sprintf("%.0f", 100*d$SRPBIO_o[1]),
                sprintf("%.0f", 100*d$SRPEFF_o[1])))
  })

  output$verdict <- renderText({
    d <- sim(); i1 <- nrow(d)
    healed <- d$INFLAM_o[i1] < 1.6
    paste0(
      if (healed) "The site fell into the HEALTH basin and stayed there.\n"
      else        "The site RELAPSED into the disease basin.\n",
      sprintf("  Inflammatory tone at end of horizon : %.2f (health = 1)\n", d$INFLAM_o[i1]),
      sprintf("  Bone lost over the horizon          : %.2f mm (never returns)\n",
              d$BLOSS[i1] - d$BLOSS[1]),
      sprintf("  Attachment lost over the horizon    : %+.2f mm\n",
              d$CAL_m[i1] - d$CAL_m[1]),
      sprintf("  Permanent instrumentation scar      : %.2f mm\n", d$SCAR[i1]),
      "\nBone loss is the only irreversible output in this model. Everything ",
      "else you see relax back is a state variable; the bone is an integral.")
  })

  ## ---- 2. PK ---------------------------------------------------------------
  output$pk_plasma <- renderPlot({
    d <- sim()
    d$CDOX <- d$DOXC/100; d$CAMX <- d$AMXC/25
    d$CMTZ <- d$MTZC/50;  d$CTNF <- d$TNFIC/3; d$CDMB <- d$DMABC/3.5
    ln(d, c("CDOX","CAMX","CMTZ","CTNF","CDMB"),
       c("Doxycycline (SDD)","Amoxicillin","Metronidazole","Anti-TNF","Denosumab"),
       c(PAL["drug"],PAL["base"],PAL["alt"],PAL["purple"],PAL["warn"]),
       "Plasma concentration (mg/L)", "Systemic pharmacokinetics")
  })
  output$pk_local <- renderPlot({
    d <- sim()
    ln(d, c("DOXG","MINLOC","AMYLOC","RVELOC","CHXLOC"),
       c("Doxycycline in GCF","Minocycline in pocket","AMY-101 (local)",
         "Resolvin E1 (local)","Chlorhexidine (supragingival)"),
       c(PAL["drug"],PAL["base"],PAL["purple"],PAL["alt"],PAL["grey"]),
       "Local concentration / effect-site units", "Local delivery")
  })

  ## ---- 3. microbiology -----------------------------------------------------
  output$micro <- renderPlot({
    d <- sim()
    ln(d, c("BIOD","BIOS","CALC"),
       c("Dysbiotic biomass","Symbiont biomass","Calculus reservoir"),
       c(PAL["drug"],PAL["alt"],PAL["warn"]),
       "Normalised units", "Biofilm ecology")
  })
  output$keystone <- renderPlot({
    d <- sim()
    d$PGpct <- 100*d$PG
    ln(d, c("PGpct","GING","LPSN"),
       c("P. gingivalis (% of biofilm)","Gingipain activity","Tissue PAMP load"),
       c(PAL["drug"],PAL["purple"],PAL["base"]),
       "Units", "Keystone pathogen: low abundance, high leverage")
  })

  ## ---- 4. innate -----------------------------------------------------------
  output$innate <- renderPlot({
    d <- sim()
    ln(d, c("PMN","PKN_o","C5A"),
       c("Neutrophils recruited","Killing capacity (x health)","Local C5a"),
       c(PAL["base"],PAL["alt"],PAL["drug"]),
       "Normalised units", "Recruited but disarmed")
  })
  output$clearance <- renderPlot({
    d <- sim()
    d$CLEAR <- d$PMN * d$PKN_o
    ln(d, c("PMN","CLEAR"),
       c("Neutrophils present","Effective bacterial clearance (PMN x capacity)"),
       c(PAL["grey"],PAL["drug"]),
       "Normalised units",
       "Why more neutrophils is not more clearance")
  })

  ## ---- 5. mediators --------------------------------------------------------
  output$mediators <- renderPlot({
    d <- sim()
    ln(d, c("IL1B","TNFA","IL17","PGE2","INFLAM_o"),
       c("IL-1beta","TNF-alpha","IL-17A","PGE2","Composite tone"),
       c(PAL["drug"],PAL["base"],PAL["purple"],PAL["warn"],"black"),
       "Normalised units (health = 1)", "Cytokine network")
  })
  output$mmp <- renderPlot({
    d <- sim()
    ln(d, c("MMP8A_o","TIMP","RATIO_o","OCL"),
       c("Active MMP-8","TIMP-1","Free RANKL : OPG","Osteoclast pool"),
       c(PAL["drug"],PAL["alt"],PAL["purple"],PAL["warn"]),
       "Normalised units (health = 1)", "Proteolysis and osteoclastogenesis")
  })

  ## ---- 6. clinical ---------------------------------------------------------
  output$clinical <- renderPlot({
    d <- sim()
    ln(d, c("PPD_o","CAL_m","BLOSS","GM"),
       c("Probing depth (measured)","Attachment level (measured)",
         "Alveolar bone loss","Gingival margin"),
       c(PAL["base"],PAL["drug"],PAL["warn"],PAL["alt"]),
       "mm from the cemento-enamel junction", "Clinical endpoints")
  })
  output$decomp <- renderTable({
    d <- sim(); i0 <- 1; i1 <- nrow(d)
    dPPD <- d$PPD_o[i0] - d$PPD_o[i1]
    dCAL <- d$CAL[i0]   - d$CAL[i1]
    dGM  <- d$GM[i1]    - d$GM[i0]
    dPEN <- (d$PPD_o[i0] - d$PPDT_o[i0]) - (d$PPD_o[i1] - d$PPDT_o[i1])
    data.frame(Component = c("True attachment gain", "Recession",
                             "Reduced probe penetration", "TOTAL"),
               mm = round(c(dCAL, dGM, dPEN, dPPD), 3),
               `Percent of total` = round(100*c(dCAL, dGM, dPEN, dPPD)/dPPD, 1),
               check.names = FALSE)
  })
  output$ratchet <- renderPlot({
    d <- sim()
    plot(d$time, d$BLOSS, type = "l", lwd = 3, col = PAL["warn"],
         xlab = "Time (days)", ylab = "mm", las = 1, bty = "l",
         ylim = range(c(d$BLOSS, d$CAL, d$SCAR)),
         main = "The ratchet: what recovers and what does not")
    grid(col = "#eeeeee")
    lines(d$time, d$CAL,  lwd = 2, col = PAL["drug"])
    lines(d$time, d$SCAR, lwd = 2, col = PAL["grey"], lty = 2)
    lines(d$time, d$DIST_o, lwd = 2, col = PAL["base"], lty = 3)
    legend("topleft",
           c("Alveolar bone loss (monotone)", "Attachment level (partly recovers)",
             "Permanent instrumentation scar", "Plaque front to bone crest (mm)"),
           col = c(PAL["warn"], PAL["drug"], PAL["grey"], PAL["base"]),
           lwd = 2, lty = c(1,1,2,3), bty = "n", cex = 0.85)
  })

  ## ---- 7. systemic ---------------------------------------------------------
  output$systemic <- renderPlot({
    d <- sim()
    ln(d, c("IL6S","CRP","INSR","HBA1C","FMD"),
       c("Systemic IL-6","CRP (mg/L)","Insulin resistance","HbA1c (%)","FMD (%)"),
       c(PAL["base"],PAL["drug"],PAL["warn"],PAL["purple"],PAL["alt"]),
       "Mixed units", "Systemic consequences of the periodontal lesion")
  })
  output$pisa <- renderPlot({
    d <- sim()
    ln(d, c("PISA_o"), "Periodontal inflamed surface area",
       PAL["drug"], "mm2",
       "PISA: the dose the rest of the body receives")
  })

  ## ---- 8. separatrix -------------------------------------------------------
  sep <- eventReactive(input$sep_go, {
    PDT_separatrix(MOD, bloss0 = input$bloss0, param = pars())
  })
  output$sep <- renderPrint({ sep() })

  ## ---- 9. scenarios --------------------------------------------------------
  scn <- eventReactive(input$scn_go, {
    PDT_simulate_scenarios(MOD, end = min(input$end, 730), delta = 5)
  })
  output$scenarios <- renderPlot({
    s <- scn(); labs <- unique(s$scenario)
    cols <- rainbow(length(labs), v = 0.8)
    plot(range(s$time), range(s$PPD_o), type = "n", xlab = "Time (days)",
         ylab = "Probing depth (mm)", las = 1, bty = "l",
         main = "Scenario comparison, deep site")
    grid(col = "#eeeeee")
    for (i in seq_along(labs)) {
      d <- s[s$scenario == labs[i], ]
      lines(d$time, d$PPD_o, col = cols[i], lwd = 2)
    }
    legend("topright", labs, col = cols, lwd = 2, bty = "n", cex = 0.7)
  })
  output$scn_tab <- renderTable({
    s <- scn()
    do.call(rbind, lapply(split(s, s$scenario), function(d) {
      i1 <- nrow(d)
      data.frame(Scenario = d$scenario[1],
                 `PPD baseline` = round(d$PPD_o[1], 2),
                 `PPD end` = round(d$PPD_o[i1], 2),
                 `CAL gain` = round(d$CAL_m[1] - d$CAL_m[i1], 2),
                 `Bone lost` = round(d$BLOSS[i1] - d$BLOSS[1], 2),
                 `BOP %` = round(d$BOPs_o[i1]),
                 check.names = FALSE)
    }))
  })

  ## ---- 10. population ------------------------------------------------------
  pop <- eventReactive(input$pop_go, {
    list(SRP     = PDT_site_population(MOD, "SRP",     param = pars()),
         SURGERY = PDT_site_population(MOD, "SURGERY", param = pars()))
  })
  output$bands <- renderTable({
    p <- pop()
    rbind(cbind(Therapy = "SRP",     PDT_bands(p$SRP)),
          cbind(Therapy = "Surgery", PDT_bands(p$SURGERY)))
  })
  output$critdepth <- renderPlot({
    cd <- PDT_critical_depth(MOD, pops = pop())
    plot(cd$PPD0, cd$CALgain_SRP, type = "l", lwd = 3, col = PAL["base"],
         ylim = range(c(cd$CALgain_SRP, cd$CALgain_SURG)),
         xlab = "Baseline probing depth (mm)",
         ylab = "Attachment gain at 12 months (mm)", las = 1, bty = "l",
         main = "Critical probing depth (emergent, not fitted)")
    grid(col = "#eeeeee")
    lines(cd$PPD0, cd$CALgain_SURG, lwd = 3, col = PAL["drug"])
    abline(h = 0, lty = 2, col = "#999999")
    legend("topleft", c("Scaling and root planing", "Open flap surgery",
                        sprintf("critical PD, SRP = %.2f mm", attr(cd,"critical_PD_SRP")),
                        sprintf("critical PD, surgery = %.2f mm", attr(cd,"critical_PD_SURGERY"))),
           col = c(PAL["base"], PAL["drug"], NA, NA), lwd = c(3,3,NA,NA), bty = "n", cex = 0.85)
  })

  ## ---- 11. calibration -----------------------------------------------------
  cal <- eventReactive(input$cal_go, { PDT_calibration_report(MOD, verbose = FALSE) })
  output$cal <- renderTable({ cal() })
}

shinyApp(ui, server)
