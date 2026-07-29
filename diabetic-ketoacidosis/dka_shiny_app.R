## =============================================================================
##  dka_shiny_app.R
##  Interactive dashboard for the Diabetic Ketoacidosis / HHS QSP model
##  (42 ODEs; see dka_mrgsolve_model.R for the equations and calibration notes)
## =============================================================================
##
##  DESIGN PRINCIPLE OF THIS DASHBOARD
##  ----------------------------------
##  Every screen is built so that the user cannot look at one number in
##  isolation, because the central claim of the model is that the numbers we
##  measure in a hyperglycaemic crisis resolve in the reverse order of their
##  information content.  Concretely:
##
##    * the glucose panel always plots the ketone on the same axis of time, with
##      the resolution thresholds marked, so the lag is impossible to miss;
##    * the acid-base panel shows the strong-ion LEDGER, not just the pH, so it
##      is visible that bicarbonate is a small residual of large terms;
##    * the potassium panel always shows the total-body deficit next to the
##      serum value, because that is the whole point about potassium;
##    * the ketone panel separates BHB from acetoacetate and overlays the
##      "nitroprusside-visible" fraction, so the strip paradox is on screen;
##    * every treatment trajectory is drawn against its own UNTREATED
##      counterfactual from the identical presenting state.
##
##  The presenting state is never entered by the user.  The user sets the
##  PATHOPHYSIOLOGY (residual beta-cell function, water access, illness
##  severity, drugs, body weight, duration of the prodrome) and the model
##  computes what walks through the door.
##
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  Run:      shiny::runApp("dka_shiny_app.R")
##  EDUCATIONAL / RESEARCH USE ONLY.
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

source("dka_mrgsolve_model.R", local = TRUE)   # builds `dka`, patient(), etc.

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey93"),
        legend.position = "bottom")

PAL <- c(treated = "#1f6fb4", untreated = "#b02418", target = "#1d6b3a",
         a = "#1f6fb4", b = "#b02418", c = "#8a6d1a", d = "#5b2d8e",
         e = "#14655f", f = "#a04a12")

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Diabetic Ketoacidosis / Hyperglycaemic Hyperosmolar State — QSP model"),
  tags$p(style = "color:#666;margin-top:-8px;",
         HTML("Two deficits, two disposal architectures, one cause. ",
              "Glucose has two exits (insulin-mediated disposal and a renal ",
              "escape valve whose conductance is the extracellular volume); ",
              "ketoacid has one, and it saturates. <b>Educational use only.</b>")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("1. The patient (pathophysiology)"),
      sliderInput("BW", "Body weight (kg)", 20, 120, 70, 5),
      sliderInput("BETA", "Residual beta-cell function (0 = type 1)",
                  0, 0.6, 0, 0.02),
      sliderInput("WATER", "Access to / drive for oral water",
                  0, 1, 1, 0.05),
      sliderInput("ILL0", "Severity of the precipitating illness",
                  0, 1, 0.55, 0.05),
      sliderInput("PRODROME", "Hours of untreated evolution before arrival",
                  6, 72, 24, 2),
      checkboxInput("SGLT2", "SGLT2 inhibitor on board", FALSE),
      checkboxInput("ETOH", "Ethanol (raises hepatic NADH/NAD+)", FALSE),
      sliderInput("GFRMAX", "Baseline GFR (L/h; 7.2 = 120 mL/min)",
                  1.5, 9, 7.2, 0.3),

      hr(),
      h4("2. The protocol (treatment)"),
      sliderInput("INS", "Insulin infusion (U/kg/h)", 0, 0.4, 0.10, 0.01),
      selectInput("FLUID", "Resuscitation fluid",
                  c("0.9% saline" = "NS", "Plasma-Lyte 148" = "PLASMALYTE",
                    "Ringer's lactate" = "LR", "0.45% saline" = "HALF_NS")),
      sliderInput("FH", "First-hour bolus (mL/kg)", 0, 30, 15, 1),
      sliderInput("MAINT", "Maintenance rate (mL/h)", 50, 900, 250, 25),
      sliderInput("KCL", "KCl in the fluid (mmol/L)", 0, 60, 40, 5),
      sliderInput("BIC", "NaHCO3 bolus over 2 h (mmol)", 0, 300, 0, 25),
      checkboxInput("KRULE", "Apply the ADA rule: hold insulin if K < 3.3",
                    TRUE),
      sliderInput("HRS", "Hours of treatment to simulate", 6, 48, 24, 2),

      hr(),
      actionButton("go", "Simulate", class = "btn-primary btn-block"),
      tags$p(style = "font-size:11px;color:#888;margin-top:8px;",
             "The presenting laboratory values are never entered here. They are",
             "integrated from a healthy steady state with insulin withdrawn,",
             "and whatever that lands on is where treatment starts.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---- 1 ------------------------------------------------------------
        tabPanel(
          "1. Presentation",
          h4("What walked through the door — computed, not entered"),
          fluidRow(column(7, DTOutput("tbl_present")),
                   column(5, plotOutput("plt_leadin", height = "380px"))),
          hr(),
          h5("Severity classification"),
          verbatimTextOutput("txt_severity"),
          tags$p(style = "color:#666;font-size:12px;",
                 "The prodrome is a positive feedback loop: hyperglycaemia",
                 "drives an osmotic diuresis, the diuresis contracts the",
                 "extracellular volume, the fall in GFR closes the renal escape",
                 "valve, and the glucose rises further. Ketosis then closes the",
                 "second loop by stopping the patient drinking.")
        ),

        ## ---- 2 ------------------------------------------------------------
        tabPanel(
          "2. Insulin PK/PD",
          h4("Plasma insulin, its two effect sites, and the four half-maximal constants"),
          plotOutput("plt_pk", height = "300px"),
          plotOutput("plt_pd", height = "300px"),
          tags$p(style = "color:#666;font-size:12px;",
                 "Suppression of lipolysis is half-maximal at 15 uU/mL;",
                 "stimulation of glucose disposal at 60. Everything the model",
                 "says about low-dose insulin follows from that ordering.",
                 "Peripheral actions follow plasma insulin within minutes;",
                 "hepatic glucose output follows over hours.")
        ),

        ## ---- 3 ------------------------------------------------------------
        tabPanel(
          "3. Glucose — and why it misleads",
          h4("Glucose and beta-hydroxybutyrate on the same time axis"),
          plotOutput("plt_glu_bhb", height = "340px"),
          h5("Where the glucose actually goes"),
          plotOutput("plt_gflux", height = "280px"),
          DTOutput("tbl_glucose")
        ),

        ## ---- 4 ------------------------------------------------------------
        tabPanel(
          "4. Ketones and the strip paradox",
          h4("BHB, acetoacetate, acetone, and the nitroprusside-visible fraction"),
          plotOutput("plt_ketones", height = "330px"),
          plotOutput("plt_ratio", height = "260px"),
          tags$p(style = "color:#666;font-size:12px;",
                 "A urine nitroprusside strip detects acetoacetate and acetone",
                 "but not BHB. As hepatic redox normalises during treatment, BHB",
                 "is converted to acetoacetate, so the strip improves far more",
                 "slowly than the disease does — and can appear to worsen.")
        ),

        ## ---- 5 ------------------------------------------------------------
        tabPanel(
          "5. Acid-base ledger",
          h4("Bicarbonate is a residual, not a pool"),
          plotOutput("plt_ledger", height = "330px"),
          plotOutput("plt_ab", height = "300px"),
          DTOutput("tbl_ab"),
          tags$p(style = "color:#666;font-size:12px;",
                 "Each bar decomposes the electroneutrality condition. The",
                 "bicarbonate is the small difference left over, which is why it",
                 "is the least robust of the acid-base variables once chloride",
                 "is being infused, and why the anion gap is the better bedside",
                 "quantity.")
        ),

        ## ---- 6 ------------------------------------------------------------
        tabPanel(
          "6. Electrolytes and the potassium deficit",
          h4("Serum potassium next to total-body potassium"),
          plotOutput("plt_k", height = "320px"),
          h4("Sodium, chloride, and iatrogenic hyperchloraemia"),
          plotOutput("plt_nacl", height = "300px"),
          DTOutput("tbl_lytes")
        ),

        ## ---- 7 ------------------------------------------------------------
        tabPanel(
          "7. Volume, kidney and the escape valve",
          h4("Extracellular volume, GFR, urine flow and renal glucose clearance"),
          plotOutput("plt_vol", height = "330px"),
          h4("The escape valve: renal glucose clearance versus ECF volume"),
          plotOutput("plt_valve", height = "300px")
        ),

        ## ---- 8 ------------------------------------------------------------
        tabPanel(
          "8. Brain and mental status",
          h4("Osmolality, the brain's osmolyte store, brain water and GCS"),
          plotOutput("plt_brain", height = "340px"),
          h5("Fluid-rate arms (the PECARN FLUID question)"),
          plotOutput("plt_fluidrate", height = "280px"),
          tags$p(style = "color:#666;font-size:12px;",
                 "The osmolyte store empties over days while plasma osmolality",
                 "is corrected over hours. In this model that rate mismatch is",
                 "the smaller of two contributions; the larger is an ischaemic",
                 "insult accrued BEFORE treatment starts, which is why a slower",
                 "drip does not protect the brain.")
        ),

        ## ---- 9 ------------------------------------------------------------
        tabPanel(
          "9. Resolution and scenario comparison",
          h4("When does each criterion cross?"),
          DTOutput("tbl_resolution"),
          h4("Insulin dose-response"),
          plotOutput("plt_dose", height = "320px"),
          DTOutput("tbl_dose"),
          h4("Crystalloid comparison"),
          plotOutput("plt_fluidcmp", height = "300px")
        ),

        ## ---- 10 -----------------------------------------------------------
        tabPanel(
          "10. Phenotype atlas",
          h4("One model, several diseases — set by two knobs"),
          DTOutput("tbl_pheno"),
          plotOutput("plt_pheno", height = "360px"),
          tags$p(style = "color:#666;font-size:12px;",
                 "Residual beta-cell function restrains lipolysis (and, through",
                 "portal insulin, the CPT-1 gate) long before it restrains",
                 "glucose. Water access sets how far the glucose can run. Those",
                 "two knobs move the same equations from ketoacidosis to a",
                 "hyperosmolar state; an SGLT2 inhibitor holds the renal valve",
                 "open and produces a euglycaemic presentation.")
        ),

        ## ---- 11 -----------------------------------------------------------
        tabPanel(
          "11. Model card",
          h4("What this model is, and where not to trust it"),
          htmlOutput("txt_card")
        )
      )
    )
  )
)

## ---------------------------------------------------------------------------
## SERVER
## ---------------------------------------------------------------------------
server <- function(input, output, session) {

  pat <- reactive({
    patient(BW = input$BW, BETA = input$BETA, WATER = input$WATER,
            ILL0 = input$ILL0, SGLT2 = as.numeric(input$SGLT2),
            ALCOHOL = if (input$ETOH) 1.2 else 0,
            GFRMAX = input$GFRMAX)
  })

  ## ---- lead-in (the presentation) -------------------------------------
  lead <- eventReactive(input$go, {
    leadin(pat(), input$PRODROME)
  }, ignoreNULL = FALSE)

  init0 <- reactive(start_from(lead()))

  ## ---- treated and untreated counterfactual ---------------------------
  treated <- eventReactive(input$go, {
    run_protocol(pat(), init0(), hours = input$HRS,
                 fluid = input$FLUID,
                 dex_fluid = if (input$FLUID %in% c("PLASMALYTE", "LR"))
                   "D5PL" else "D5HALF",
                 first_hour = input$FH / 1000 * input$BW,
                 maint = input$MAINT / 1000,
                 ins = input$INS, ins_taper = input$INS / 2,
                 kcl = input$KCL, k_rule = input$KRULE,
                 bicarb_mmol = input$BIC)
  }, ignoreNULL = FALSE)

  untreated <- eventReactive(input$go, {
    p <- modifyList(pat(), c(list(RATE_FL = 0, INS_IV = 0, KCL = 0, BICARB = 0),
                            fluid_pars("NS")))
    dka %>% param(p) %>% init(init0()) %>%
      mrgsim(end = input$HRS, delta = 0.25) %>% as_tibble()
  }, ignoreNULL = FALSE)

  both <- reactive({
    bind_rows(mutate(treated(), arm = "treated"),
              mutate(untreated(), arm = "untreated"))
  })

  ## ---- 1. presentation -------------------------------------------------
  output$tbl_present <- renderDT({
    r <- lead()[nrow(lead()), ]
    tibble::tibble(
      quantity = c("glucose (mg/dL)", "pH", "bicarbonate (mmol/L)",
                   "anion gap (mEq/L)", "beta-hydroxybutyrate (mmol/L)",
                   "acetoacetate (mmol/L)", "BHB:AcAc ratio",
                   "sodium, measured (mmol/L)", "sodium, corrected (mmol/L)",
                   "potassium (mmol/L)", "chloride (mmol/L)",
                   "PCO2 (mmHg)", "effective osmolality (mOsm/kg)",
                   "urea nitrogen (mg/dL)",
                   "creatinine, true (mg/dL)",
                   "creatinine, Jaffe assay (mg/dL)",
                   "GFR (mL/min)", "extracellular volume (L)",
                   "total body water deficit (L)",
                   "total body potassium (mmol)",
                   "retained ketoanion = potential bicarbonate (mmol)",
                   "cumulative urine output (L)",
                   "Glasgow Coma Scale"),
      value = round(c(r$GLUmgdl, r$PHa, r$BICARB_mM, r$ANIONGAP, r$BHBmM,
                      r$ACACmM, r$BHB_ACAC, r$NAmM, r$NA_CORR, r$KmM, r$CLmM,
                      r$PCO2, r$OSMeff, r$BUN, r$CREA_TRUE, r$CREA_JAFFE,
                      r$GFRmlmin, r$VECFL,
                      pat()$FTBW * pat()$BW - r$TBWL, r$KTOT,
                      r$POTBICARB, r$UVOL, r$GCS), 2))
  }, options = list(dom = "t", pageLength = 25), rownames = FALSE)

  output$plt_leadin <- renderPlot({
    d <- lead() %>%
      select(time, glucose = GLUmgdl, BHB = BHBmM, HCO3 = BICARB_mM,
             `GFR (mL/min)` = GFRmlmin) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.9, colour = PAL["b"]) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "hours of untreated evolution", y = NULL,
           title = "The prodrome, integrated") + THEME
  })

  output$txt_severity <- renderPrint({
    r <- lead()[nrow(lead()), ]
    sev <- if (r$PHa < 7.00 || r$BICARB_mM < 10) "SEVERE" else
      if (r$PHa < 7.24 || r$BICARB_mM < 15) "MODERATE" else "MILD"
    hhs <- r$OSMeff > 320 && r$BHBmM < 3.0
    cat(sprintf("ADA severity band (by pH and bicarbonate): %s\n", sev))
    cat(sprintf("Meets hyperosmolar (HHS) pattern (osm > 320, BHB < 3.0): %s\n",
                ifelse(hhs, "YES", "no")))
    cat(sprintf("Euglycaemic pattern (glucose < 250 with gap > 16): %s\n",
                ifelse(r$GLUmgdl < 250 && r$ANIONGAP > 16, "YES", "no")))
    cat(sprintf("\nWater deficit %.1f L (%.1f%% of body weight); potassium deficit %.0f mmol (%.1f mmol/kg)\n",
                pat()$FTBW * pat()$BW - r$TBWL,
                100 * (pat()$FTBW * pat()$BW - r$TBWL) / pat()$BW,
                pat()$K0 * pat()$FECF * pat()$BW +
                  pat()$KIC0 * (pat()$FTBW - pat()$FECF) * pat()$BW - r$KTOT,
                (pat()$K0 * pat()$FECF * pat()$BW +
                   pat()$KIC0 * (pat()$FTBW - pat()$FECF) * pat()$BW - r$KTOT) /
                  pat()$BW))
  })

  ## ---- 2. insulin -------------------------------------------------------
  output$plt_pk <- renderPlot({
    d <- treated() %>%
      select(time, `plasma insulin` = INSPuU,
             `peripheral effect site` = INSEFuU,
             `portal (endogenous privilege)` = IPORTuU) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 15, linetype = 2, colour = PAL["target"]) +
      geom_hline(yintercept = 60, linetype = 3, colour = PAL["b"]) +
      annotate("text", x = Inf, y = 15, hjust = 1.05, vjust = -0.5, size = 3,
               label = "lipolysis IC50 = 15", colour = PAL["target"]) +
      annotate("text", x = Inf, y = 60, hjust = 1.05, vjust = -0.5, size = 3,
               label = "disposal EC50 = 60", colour = PAL["b"]) +
      scale_y_log10() + scale_colour_manual(values = unname(PAL[c("a","c","d")])) +
      labs(x = "hours of treatment", y = "uU/mL (log)", colour = NULL) + THEME
  })

  output$plt_pd <- renderPlot({
    d <- treated() %>%
      select(time, `lipolysis suppression` = fLIPsupp,
             `disposal stimulation` = fUPstim, `CPT-1 gate open` = CPT1gate) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = unname(PAL[c("d","a","b")])) +
      labs(x = "hours of treatment", y = "fraction of maximum",
           colour = NULL,
           title = "The antiketogenic arm saturates; the glucose arm does not") +
      THEME
  })

  ## ---- 3. glucose -------------------------------------------------------
  output$plt_glu_bhb <- renderPlot({
    d <- both() %>%
      select(time, arm, `glucose (mg/dL)` = GLUmgdl, `BHB (mmol/L)` = BHBmM) %>%
      pivot_longer(c(-time, -arm))
    thr <- tibble::tibble(name = c("glucose (mg/dL)", "BHB (mmol/L)"),
                          y = c(250, 0.6))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      geom_hline(data = thr, aes(yintercept = y), linetype = 2,
                 colour = PAL["target"]) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = unname(PAL[c("treated", "untreated")])) +
      labs(x = "hours of treatment", y = NULL, colour = NULL,
           title = "Glucose reaches target hours before the ketosis clears") +
      THEME
  })

  output$plt_gflux <- renderPlot({
    d <- treated() %>%
      select(time, `hepatic output` = HGPflux, `total disposal` = UPTflux,
             `urinary loss` = UGLUrate) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = unname(PAL[c("c","a","e")])) +
      labs(x = "hours", y = "mmol/h", colour = NULL,
           title = "A large part of the early glucose fall is renal") + THEME
  })

  output$tbl_glucose <- renderDT({
    treated() %>%
      filter(time %in% c(0, 0.5, 1, 2, 4, 6, 8, 12, 24)) %>%
      transmute(hour = time, glucose = round(GLUmgdl), BHB = round(BHBmM, 2),
                HGP = round(HGPflux), disposal = round(UPTflux),
                urine = round(UGLUrate),
                `renal CL (L/h)` = round(RENAL_CL, 2),
                `GFR (mL/min)` = round(GFRmlmin))
  }, options = list(dom = "t"), rownames = FALSE)

  ## ---- 4. ketones -------------------------------------------------------
  output$plt_ketones <- renderPlot({
    d <- treated() %>%
      transmute(time, BHB = BHBmM, acetoacetate = ACACmM, acetone = ACETmM,
                `nitroprusside-visible (AcAc + acetone)` = ACACmM + ACETmM,
                `total ketone` = KETmM) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = unname(PAL[c("a","b","c","d","e")])) +
      labs(x = "hours", y = "mmol/L", colour = NULL) + THEME
  })

  output$plt_ratio <- renderPlot({
    d <- treated() %>% transmute(time, `BHB:AcAc` = BHB_ACAC,
                                 `hepatic NADH/NAD+ index` = REDOX) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      scale_colour_manual(values = unname(PAL[c("b","d")])) +
      labs(x = "hours", y = NULL, colour = NULL,
           title = "As redox normalises, BHB converts to the species the strip sees") +
      THEME
  })

  ## ---- 5. acid-base -----------------------------------------------------
  output$plt_ledger <- renderPlot({
    d <- treated() %>% filter(time %in% c(0, 2, 4, 8, 12, 24))
    conc <- (pat()$FECF * pat()$BW) / d$VECFL
    led <- bind_rows(
      tibble::tibble(time = d$time, term = "strong ion difference",
                     value = d$SIDmEq),
      tibble::tibble(time = d$time, term = "weak acid (albumin, phosphate)",
                     value = -(pat()$ALB * conc * (0.123 * d$PHa - 0.631) +
                                 d$PHOSmM * (0.309 * d$PHa - 0.469))),
      tibble::tibble(time = d$time, term = "ketoanion", value = -d$KETmM),
      tibble::tibble(time = d$time, term = "lactate", value = -d$LACmM),
      tibble::tibble(time = d$time, term = "other strong anion",
                     value = -pat()$SIGO * conc),
      tibble::tibble(time = d$time, term = "= bicarbonate",
                     value = d$BICARB_mM))
    led$term <- factor(led$term, levels = c("strong ion difference",
      "weak acid (albumin, phosphate)", "ketoanion", "lactate",
      "other strong anion", "= bicarbonate"))
    ggplot(led, aes(factor(time), value, fill = term)) +
      geom_col(position = "dodge") +
      geom_hline(yintercept = 0, colour = "grey40") +
      labs(x = "hours of treatment", y = "mEq/L", fill = NULL,
           title = "The electroneutrality ledger: bicarbonate is what is left over") +
      THEME
  })

  output$plt_ab <- renderPlot({
    d <- both() %>% select(time, arm, pH = PHa, `HCO3 (mmol/L)` = BICARB_mM,
                           `anion gap` = ANIONGAP, `PCO2 (mmHg)` = PCO2) %>%
      pivot_longer(c(-time, -arm))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = unname(PAL[c("treated", "untreated")])) +
      labs(x = "hours", y = NULL, colour = NULL) + THEME
  })

  output$tbl_ab <- renderDT({
    treated() %>% filter(time %in% c(0, 1, 2, 4, 6, 8, 12, 18, 24)) %>%
      transmute(hour = time, pH = round(PHa, 3), HCO3 = round(BICARB_mM, 1),
                PCO2 = round(PCO2), `anion gap` = round(ANIONGAP, 1),
                SID = round(SIDmEq, 1), ketone = round(KETmM, 2),
                Cl = round(CLmM, 1),
                `NH4 excretion (mmol/h)` = round(UNH4rate, 2))
  }, options = list(dom = "t"), rownames = FALSE)

  ## ---- 6. electrolytes --------------------------------------------------
  output$plt_k <- renderPlot({
    tot0 <- pat()$K0 * pat()$FECF * pat()$BW +
      pat()$KIC0 * (pat()$FTBW - pat()$FECF) * pat()$BW
    d <- bind_rows(
      treated() %>% transmute(time, name = "serum K (mmol/L)", value = KmM),
      treated() %>% transmute(time, name = "total-body K deficit (mmol)",
                              value = tot0 - KTOT))
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.9, colour = PAL["a"]) +
      geom_hline(data = tibble::tibble(name = "serum K (mmol/L)", y = 3.3),
                 aes(yintercept = y), linetype = 2, colour = PAL["b"]) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "hours", y = NULL,
           title = "The serum value and the deficit move in opposite directions") +
      THEME
  })

  output$plt_nacl <- renderPlot({
    d <- both() %>% select(time, arm, `Na measured` = NAmM,
                           `Na corrected` = NA_CORR, `Cl` = CLmM,
                           `phosphate` = PHOSmM) %>%
      pivot_longer(c(-time, -arm))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = unname(PAL[c("treated", "untreated")])) +
      labs(x = "hours", y = "mmol/L", colour = NULL) + THEME
  })

  output$tbl_lytes <- renderDT({
    treated() %>% filter(time %in% c(0, 2, 4, 8, 12, 24)) %>%
      transmute(hour = time, Na = round(NAmM, 1),
                `Na corrected` = round(NA_CORR, 1), Cl = round(CLmM, 1),
                K = round(KmM, 2), `total body K` = round(KTOT),
                `urinary K (mmol/h)` = round(UKrate, 1),
                phosphate = round(PHOSmM, 2))
  }, options = list(dom = "t"), rownames = FALSE)

  ## ---- 7. volume / kidney ----------------------------------------------
  output$plt_vol <- renderPlot({
    d <- both() %>% select(time, arm, `ECF volume (L)` = VECFL,
                           `GFR (mL/min)` = GFRmlmin,
                           `urine flow (L/h)` = UVLh,
                           `renal glucose CL (L/h)` = RENAL_CL) %>%
      pivot_longer(c(-time, -arm))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = unname(PAL[c("treated", "untreated")])) +
      labs(x = "hours", y = NULL, colour = NULL) + THEME
  })

  output$plt_valve <- renderPlot({
    waters <- c(1, 0.7, 0.5, 0.35, 0.25, 0.15, 0.05)
    d <- bind_rows(lapply(waters, function(w) {
      s <- leadin(patient(BW = input$BW, BETA = input$BETA, WATER = w,
                          ILL0 = input$ILL0, GFRMAX = input$GFRMAX), 36)
      r <- s[nrow(s), ]
      tibble::tibble(water = w, glucose = r$GLUmgdl, ecf = r$VECFL,
                     gfr = r$GFRmlmin, renal_cl = r$RENAL_CL,
                     osm = r$OSMeff, gcs = r$GCS)
    }))
    ggplot(d, aes(ecf, glucose)) +
      geom_line(colour = PAL["c"], linewidth = 1) +
      geom_point(aes(size = renal_cl, colour = gcs)) +
      scale_colour_gradient(low = PAL["b"], high = PAL["target"]) +
      labs(x = "extracellular volume at 36 h (L)",
           y = "glucose (mg/dL)", size = "renal glucose\nclearance (L/h)",
           colour = "GCS",
           title = "Same insulin deficit; only the drinking capacity differs") +
      THEME
  })

  ## ---- 8. brain ---------------------------------------------------------
  output$plt_brain <- renderPlot({
    d <- both() %>% select(time, arm, `effective osmolality` = OSMeff,
                           `brain osmolyte store (mOsm/kg)` = OSMB,
                           `intracranial pressure (mmHg)` = ICPmmHg,
                           `Glasgow Coma Scale` = GCS) %>%
      pivot_longer(c(-time, -arm))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = unname(PAL[c("treated", "untreated")])) +
      labs(x = "hours", y = NULL, colour = NULL) + THEME
  })

  output$plt_fluidrate <- renderPlot({
    arms <- list(slow = c(10, 125), standard = c(15, 250),
                 fast = c(20, 500), `very fast` = c(30, 750))
    d <- bind_rows(lapply(names(arms), function(n) {
      a <- arms[[n]]
      s <- run_protocol(pat(), init0(), hours = 12,
                        first_hour = a[1] / 1000 * input$BW,
                        maint = a[2] / 1000, ins = input$INS,
                        ins_taper = input$INS / 2, kcl = input$KCL)
      mutate(select(s, time, ICPmmHg, OSMeff, GCS), arm = n)
    }))
    d$arm <- factor(d$arm, levels = names(arms))
    ggplot(d, aes(time, ICPmmHg, colour = arm)) + geom_line(linewidth = 0.9) +
      labs(x = "hours", y = "intracranial pressure (mmHg)", colour = NULL,
           title = "A 6-fold range of fluid rate barely moves the brain") + THEME
  })

  ## ---- 9. resolution / comparison --------------------------------------
  output$tbl_resolution <- renderDT({
    resolution_table(treated()) %>% mutate(hours = round(hours, 2))
  }, options = list(dom = "t"), rownames = FALSE)

  doses <- reactive({
    dd <- c(0.025, 0.05, 0.10, 0.14, 0.20, 0.40)
    lapply(setNames(dd, paste0(dd, " U/kg/h")), function(x)
      run_protocol(pat(), init0(), hours = max(24, input$HRS),
                   fluid = input$FLUID, ins = x, ins_taper = x / 2,
                   kcl = input$KCL, k_rule = input$KRULE))
  })

  output$plt_dose <- renderPlot({
    d <- bind_rows(lapply(names(doses()), function(n)
      mutate(select(doses()[[n]], time, BHBmM, GLUmgdl, KmM), dose = n)))
    d2 <- pivot_longer(d, c(BHBmM, GLUmgdl, KmM))
    ggplot(d2, aes(time, value, colour = dose)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(x = "hours", y = NULL, colour = NULL,
           title = "Ketone curves converge; glucose and potassium do not") + THEME
  })

  output$tbl_dose <- renderDT({
    bind_rows(lapply(names(doses()), function(n) {
      s <- doses()[[n]]
      tibble::tibble(dose = n, `steady-state uU/mL` = round(max(s$INSPuU)),
                     `t BHB<0.6 (h)` = round(crossing_time(s, "BHBmM", 0.6, FALSE), 1),
                     `t glucose<250 (h)` = round(crossing_time(s, "GLUmgdl", 250, FALSE), 1),
                     `K nadir` = round(min(s$KmM), 2),
                     `glucose nadir` = round(min(s$GLUmgdl)))
    }))
  }, options = list(dom = "t"), rownames = FALSE)

  output$plt_fluidcmp <- renderPlot({
    fl <- c(NS = "0.9% saline", PLASMALYTE = "Plasma-Lyte 148", LR = "Ringer's lactate")
    d <- bind_rows(lapply(names(fl), function(f) {
      s <- run_protocol(pat(), init0(), hours = max(24, input$HRS), fluid = f,
                        dex_fluid = if (f == "NS") "D5HALF" else "D5PL",
                        ins = input$INS, ins_taper = input$INS / 2,
                        kcl = input$KCL, k_rule = input$KRULE)
      mutate(select(s, time, BICARB_mM, CLmM, PHa, ANIONGAP), fluid = fl[[f]])
    })) %>% pivot_longer(c(BICARB_mM, CLmM, PHa, ANIONGAP))
    ggplot(d, aes(time, value, colour = fluid)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y") +
      scale_colour_manual(values = unname(PAL[c("a","target","c")])) +
      labs(x = "hours", y = NULL, colour = NULL,
           title = "Chloride load versus base: the two are the same axis") + THEME
  })

  ## ---- 10. phenotypes ---------------------------------------------------
  phenos <- reactive(sc10_phenotypes())

  output$tbl_pheno <- renderDT({
    bind_rows(lapply(names(phenos()), function(n)
      mutate(presentation_row(phenos()[[n]]), phenotype = n, .before = 1))) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
  }, options = list(dom = "t", scrollX = TRUE), rownames = FALSE)

  output$plt_pheno <- renderPlot({
    d <- bind_rows(lapply(names(phenos()), function(n) {
      r <- phenos()[[n]][nrow(phenos()[[n]]), ]
      tibble::tibble(phenotype = n, glucose = r$GLUmgdl, BHB = r$BHBmM,
                     osm = r$OSMeff, pH = r$PHa)
    }))
    ggplot(d, aes(glucose, BHB, colour = osm, size = 7.4 - pH)) +
      geom_point() +
      ggrepel::geom_text_repel(aes(label = phenotype), size = 3,
                               show.legend = FALSE) +
      scale_colour_gradient(low = PAL["a"], high = PAL["b"]) +
      labs(x = "glucose at presentation (mg/dL)",
           y = "beta-hydroxybutyrate (mmol/L)",
           colour = "effective\nosmolality", size = "acidaemia\n(7.4 - pH)",
           title = "The glucose-ketone plane, and where each phenotype lands") +
      THEME
  })

  ## ---- 11. model card ---------------------------------------------------
  output$txt_card <- renderUI({
    HTML('
<h5>What it is</h5>
<p>42 coupled ordinary differential equations covering insulin pharmacokinetics
and two effect sites, hepatic and peripheral glucose metabolism, adipose
lipolysis, hepatic ketogenesis through an explicit CPT-1 gate, saturable ketone
disposal, physicochemical acid-base, renal handling of glucose, ketoanion,
sodium, chloride, potassium and ammonium, extracellular and intracellular volume
with osmotic exchange, counter-regulatory hormones, a brain osmolyte store, and
clinical read-outs.</p>

<h5>Three structural choices worth knowing about</h5>
<ol>
<li><b>Bicarbonate is not a state variable.</b> It is the residual of the
strong-ion difference after weak acids, ketoanion and lactate are subtracted, so
the conversion of an organic acidosis into a hyperchloraemic one requires no
book-keeping term: it happens because excreting a ketoanion with sodium removes a
cation as well as an anion.</li>
<li><b>Endogenous insulin is privileged at the liver.</b> Because secreted
insulin passes the liver first, the CPT-1 gate sees several times the peripheral
concentration. Exogenous insulin does not. This is why residual beta-cell
function prevents ketosis while permitting extreme hyperglycaemia.</li>
<li><b>The presenting state is computed.</b> No laboratory value is entered; the
prodrome is integrated from a healthy steady state with insulin withdrawn.</li>
</ol>

<h5>Where not to trust it</h5>
<ul>
<li><b>Bicarbonate and pH are ill-conditioned here</b>, being the small
difference of large strong-ion terms. A 30% change in the fractional excretion of
chloride alone moves the presenting bicarbonate by several mmol/L while the anion
gap and the BHB barely move. Treat the gap and the ketone as the reliable
outputs. This is a property of the chemistry, not a bug.</li>
<li><b>The cerebral-oedema block is an inference, not a measurement.</b> It
reproduces the reported risk factors and the negative fluid-rate trial, but the
partition between osmotic and ischaemic contributions cannot be validated in
patients.</li>
<li><b>Ketone oxidation capacity</b> is constrained by steady-state
concentrations rather than measured directly, and the split between its saturable
and non-saturable arms is a modelling choice.</li>
<li><b>The portal:peripheral insulin ratio</b> carries the entire DKA-versus-HHS
distinction and is taken from fasting first-pass studies; whether it holds during
a crisis is unknown.</li>
<li>No pharmacogenomics, no cardiac electrophysiology beyond a potassium
read-out, no coagulation, and no infection dynamics: the precipitating illness is
a single decaying scalar.</li>
</ul>

<h5>Not for clinical use</h5>
<p>This is a teaching and hypothesis-generating model. It has not been fitted to
individual patient data and must not be used to guide treatment.</p>')
  })
}

shinyApp(ui, server)
