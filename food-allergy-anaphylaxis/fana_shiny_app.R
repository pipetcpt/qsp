## =====================================================================
##  IgE-mediated food allergy / anaphylaxis — QSP Shiny dashboard
## ---------------------------------------------------------------------
##  Run with:   shiny::runApp("fana_shiny_app.R")
##  Requires:   shiny, mrgsolve, ggplot2, dplyr, tidyr, DT
##
##  The app drives fana_mrgsolve_model.R. Nothing here is a stored
##  phenotype: every patient is GENERATED from sIgE, total IgE, age,
##  weight and comorbidity, and the eliciting dose is RE-COMPUTED by
##  bisection each time those move. That is deliberate — the point of
##  the model is that the threshold is derived, not tabulated.
##
##  Nine tabs:
##    1  Patient          — build the patient, see the surface terms
##    2  The engine       — rho, L, f and the square law made visible
##    3  Challenge        — a DBPCFC in this patient, dose by dose
##    4  Anaphylaxis      — the acute event, mediators to MAP
##    5  Adrenaline       — the timing integral
##    6  Anti-IgE         — the fast arm and the slow arm
##    7  Immunotherapy    — IgG4, anergy, Treg, and withdrawal
##    8  Threshold ladder — what every intervention is worth, in logs
##    9  Cofactors        — the within-patient variance
## =====================================================================

library(shiny)
library(mrgsolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

source("fana_mrgsolve_model.R")

theme_fa <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "#666666", size = 10),
        legend.position = "bottom")

PAL <- c("untreated"    = "#b3352f", "omalizumab" = "#2f8b52",
         "OIT"          = "#b5761a", "combination" = "#5b2f86",
         "dupilumab"    = "#2f6ea8", "quiet day"  = "#2f8b52",
         "exercise + NSAID" = "#b3352f",
         "histamine" = "#c1651b", "PAF" = "#b3352f",
         "cysLT" = "#2f8b52", "tryptase" = "#5b2f86")

MIN <- function(h) h*60

## ---------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("IgE-mediated food allergy and anaphylaxis — QSP dashboard"),
  tags$p(style = "color:#666;margin-top:-8px;",
         HTML("The eliciting dose is the output, not an input. Everything on the left \
               regenerates a patient; the threshold is re-derived by bisection.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Patient"),
      sliderInput("sige",  "Specific IgE (kU_A/L)", 0.5, 300, 40, step = 0.5),
      sliderInput("tige",  "TOTAL serum IgE (IU/mL)", 20, 4000, 300, step = 10),
      sliderInput("agey",  "Age (years)", 1, 60, 10, step = 1),
      sliderInput("wt",    "Weight (kg)", 8, 110, 35, step = 1),
      sliderInput("mcb",   "Mast cell burden (1 = normal, HaT ~2)", 0.5, 3, 1, step = 0.1),
      sliderInput("rel",   "Intrinsic releasability", 0.2, 1.5, 1, step = 0.05),
      checkboxInput("asthma", "Uncontrolled asthma", FALSE),
      checkboxInput("bblock", "Beta-blocker / ACE inhibitor", FALSE),
      sliderInput("pafah", "PAF-acetylhydrolase activity", 0.2, 1.5, 1, step = 0.05),
      hr(),
      h4("Cofactors (act on delivery)"),
      checkboxInput("cex",  "Exercise", FALSE),
      checkboxInput("cns",  "NSAID / aspirin", FALSE),
      checkboxInput("cet",  "Alcohol", FALSE),
      checkboxInput("cinf", "Acute infection / fever", FALSE),
      checkboxInput("cppi", "Proton pump inhibitor", FALSE),
      hr(),
      h4("Acute exposure"),
      sliderInput("dose", "Ingested peanut protein (mg, log scale)",
                  min = 0, max = 4.3, value = 3.3, step = 0.05,
                  pre = "10^"),
      sliderInput("epidelay", "Adrenaline delay (min; 60 = none)", 0, 60, 10, step = 1),
      checkboxInput("supine", "Supine, legs elevated", TRUE),
      checkboxInput("premed", "Cetirizine 10 mg premedication", FALSE),
      hr(),
      actionButton("recompute", "Re-compute eliciting dose",
                   class = "btn-primary", width = "100%"),
      tags$p(style = "font-size:11px;color:#888;margin-top:8px;",
             "Bisection takes a few seconds.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1 Patient",
          br(),
          fluidRow(
            column(6, wellPanel(htmlOutput("patient_card"))),
            column(6, wellPanel(htmlOutput("threshold_card")))
          ),
          plotOutput("p_surface_bar", height = "300px"),
          tags$p(style="color:#666;font-size:12px;",
                 "The three terms are multiplied and then SQUARED. A patient is dangerous \
                  because of the product, not because of any one factor.")
        ),

        tabPanel("2 The engine",
          br(),
          plotOutput("p_squarelaw", height = "340px"),
          plotOutput("p_transfer",  height = "300px"),
          tags$p(style="color:#666;font-size:12px;",
                 "Upper: the eliciting dose against the surface index, on log-log axes. \
                  The slope should be -2 — that is the model's central structural claim. \
                  Lower: the degranulation transfer function, which is what turns a smooth \
                  cross-link curve into a threshold.")
        ),

        tabPanel("3 Challenge",
          br(),
          plotOutput("p_dbpcfc", height = "460px"),
          DTOutput("t_dbpcfc")
        ),

        tabPanel("4 Anaphylaxis",
          br(),
          plotOutput("p_mediators", height = "260px"),
          plotOutput("p_physiology", height = "340px"),
          wellPanel(htmlOutput("acute_card"))
        ),

        tabPanel("5 Adrenaline",
          br(),
          plotOutput("p_epi_pk", height = "260px"),
          plotOutput("p_epi_timing", height = "330px"),
          DTOutput("t_epi")
        ),

        tabPanel("6 Anti-IgE",
          br(),
          fluidRow(
            column(4, sliderInput("oma_wk", "Weeks of omalizumab", 2, 52, 16)),
            column(4, htmlOutput("oma_dose_card")),
            column(4, htmlOutput("oma_free_card"))
          ),
          plotOutput("p_oma", height = "430px"),
          tags$p(style="color:#666;font-size:12px;",
                 "Free IgE is at target within days; the receptor takes months. \
                  Both sit inside the square, so the second arm is worth as much as the first.")
        ),

        tabPanel("7 Immunotherapy",
          br(),
          fluidRow(
            column(4, sliderInput("oit_mg", "Maintenance dose (mg/day)", 0, 600, 300, step = 10)),
            column(4, sliderInput("oit_wk", "Weeks on therapy", 4, 130, 78)),
            column(4, sliderInput("off_wk", "Weeks after stopping", 0, 78, 26))
          ),
          plotOutput("p_oit", height = "460px"),
          tags$p(style="color:#666;font-size:12px;",
                 "IgG4 comes on over months and goes off over months; mast cell anergy \
                  comes on in days and goes off in days. Desensitisation is the product \
                  of the two, which is why a missed week matters and a missed year is fatal \
                  to the effect.")
        ),

        tabPanel("8 Threshold ladder",
          br(),
          actionButton("run_ladder", "Compute the ladder (slow: ~1 min)",
                       class = "btn-warning"),
          br(), br(),
          plotOutput("p_ladder", height = "440px"),
          DTOutput("t_ladder")
        ),

        tabPanel("9 Cofactors",
          br(),
          actionButton("run_cof", "Compute cofactor shifts", class = "btn-warning"),
          br(), br(),
          plotOutput("p_cof", height = "380px"),
          plotOutput("p_cof_event", height = "320px")
        )
      )
    )
  )
)

## ---------------------------------------------------------------------
##  SERVER
## ---------------------------------------------------------------------
server <- function(input, output, session) {

  par_now <- reactive({
    patient(WT = input$wt, AGEY = input$agey,
            SIGE0 = input$sige, TIGE0 = max(input$tige, input$sige*1.2),
            MCBURDEN = input$mcb, RELEASE = input$rel,
            ASTHMA = as.numeric(input$asthma),
            BBLOCK = as.numeric(input$bblock),
            PAFAH = input$pafah,
            COF_EX = as.numeric(input$cex), COF_NSAID = as.numeric(input$cns),
            COF_ETOH = as.numeric(input$cet), COF_INF = as.numeric(input$cinf),
            COF_PPI = as.numeric(input$cppi),
            SUPINE = as.numeric(input$supine))
  })

  dose_mg <- reactive(10^input$dose)

  ## baseline surface terms (no allergen, t = 0)
  surf_now <- reactive({
    mod %>% param(par_now()) %>% mrgsim(end = 1, delta = 1) %>%
      as.data.frame() %>% slice(1)
  })

  ed_now <- eventReactive(input$recompute, {
    withProgress(message = "Bisecting for the eliciting dose", value = 0.5, {
      find_ED(par_now())
    })
  }, ignoreNULL = FALSE)

  ## ---- tab 1 --------------------------------------------------------
  output$patient_card <- renderUI({
    s <- surf_now()
    HTML(sprintf(
      "<b>Generated patient</b><br>
       specific IgE &nbsp; <b>%.1f</b> kU<sub>A</sub>/L<br>
       total IgE &nbsp;&nbsp;&nbsp;&nbsp; <b>%.0f</b> IU/mL<br>
       <b>f = sIgE/total = %.3f</b><br>
       receptor density &rho; = %.3f<br>
       occupancy L = %.3f<br>
       <b>surface index &rho;&middot;L&middot;f = %.4f</b> (%.2f&times; reference)<br>
       baseline IgG4 %.1f &mu;g/mL &nbsp; interception &times;%.2f",
      s$sIgE_IU, s$TotIgE_IU, s$Fspec, s$RHO, s$Locc_occ,
      s$SURFidx, s$SURFrel, s$IgG4_ug, s$Interc))
  })

  output$threshold_card <- renderUI({
    e <- ed_now()
    if (is.na(e)) return(HTML("<b>Eliciting dose</b><br>no reaction up to 30 g — \
                               this patient is not clinically reactive"))
    HTML(sprintf(
      "<b>Eliciting dose (grade &ge;2)</b><br>
       <span style='font-size:30px;color:#b3352f'><b>%s mg</b></span> peanut protein<br>
       &asymp; %.2f g of peanut (~%.1f%% of one kernel)<br><br>
       tolerates a single 600 mg dose: <b>%s</b><br>
       VITAL ED05 reference dose for peanut is 0.2&ndash;2 mg — \
       precautionary labelling protects the population, not this patient",
      formatC(e, format = "f", digits = 1, big.mark = ","),
      e/250, 100*e/250, ifelse(e >= 600, "YES", "no")))
  })

  output$p_surface_bar <- renderPlot({
    s <- surf_now()
    d <- data.frame(
      term = factor(c("receptor density  ρ", "occupancy  L",
                      "specific fraction  f", "product  ρ·L·f",
                      "SQUARED (the engine)"),
                    levels = c("receptor density  ρ", "occupancy  L",
                               "specific fraction  f", "product  ρ·L·f",
                               "SQUARED (the engine)")),
      value = c(s$RHO, s$Locc_occ, s$Fspec, s$SURFidx, s$SURFrel^2 * 0.1209))
    ggplot(d, aes(term, value, fill = term)) +
      geom_col(width = .62, show.legend = FALSE) +
      geom_text(aes(label = sprintf("%.4f", value)), vjust = -0.4, size = 4) +
      scale_fill_manual(values = c("#2f8b80", "#2f8b80", "#7a4fa5",
                                   "#b5761a", "#b3352f")) +
      labs(title = "What sets this patient's threshold",
           subtitle = "the three terms multiply, and then the product is squared",
           x = NULL, y = NULL) +
      expand_limits(y = max(d$value)*1.18) + theme_fa
  })

  ## ---- tab 2 --------------------------------------------------------
  output$p_squarelaw <- renderPlot({
    ss <- c(5, 10, 20, 40, 80, 160)
    d <- do.call(rbind, lapply(ss, function(x) {
      p <- modifyList(par_now(), list(SIGE0 = x))
      s <- mod %>% param(p) %>% mrgsim(end = 1, delta = 1) %>%
        as.data.frame() %>% slice(1)
      data.frame(sIgE = x, surf = s$SURFrel, ED = find_ED(p, tol = 0.05))
    }))
    d <- d[!is.na(d$ED), ]
    fit <- lm(log10(ED) ~ log10(surf), data = d)
    ggplot(d, aes(surf, ED)) +
      geom_point(size = 3.4, colour = "#b3352f") +
      geom_smooth(method = "lm", se = FALSE, colour = "#2f8b52", linewidth = .9) +
      scale_x_log10() + scale_y_log10() +
      labs(title = sprintf("The square law: slope = %.2f (theory -2.00)",
                           coef(fit)[2]),
           subtitle = "eliciting dose against the surface index, sIgE swept 5-160 kU/L",
           x = "surface index ρ·L·f (relative)",
           y = "eliciting dose (mg peanut protein)") + theme_fa
  })

  output$p_transfer <- renderPlot({
    xl <- 10^seq(-4, 0, length.out = 300)
    p  <- as.list(param(mod))
    d  <- data.frame(XL = xl,
                     degran = xl^p$HDG/(p$XL50^p$HDG + xl^p$HDG))
    ggplot(d, aes(XL, degran)) +
      geom_line(colour = "#b5761a", linewidth = 1.1) +
      geom_vline(xintercept = p$XL50, linetype = 2, colour = "#888888") +
      annotate("text", x = p$XL50, y = .95, label = "  XL50", hjust = 0,
               colour = "#666666", size = 3.6) +
      scale_x_log10() +
      labs(title = "The degranulation transfer function",
           subtitle = "a smooth cross-link curve becomes a clinical threshold here",
           x = "cross-link index X", y = "fractional degranulation rate") +
      theme_fa
  })

  ## ---- tab 3 --------------------------------------------------------
  chal <- reactive({
    mod %>% param(par_now()) %>% ev(ev_challenge_practall()) %>%
      mrgsim(end = 10, delta = 0.01) %>% as.data.frame()
  })

  output$p_dbpcfc <- renderPlot({
    d <- chal() %>%
      select(time, `skin score` = URTskin, `GI score` = GISY,
             `FEV1 (% pred)` = FEV1, `MAP (mmHg)` = MAP,
             `severity grade` = SEV, `tryptase (ng/mL)` = Tryptase) %>%
      pivot_longer(-time)
    doses <- c(3, 10, 30, 100, 300, 600, 1000)
    vl <- data.frame(time = (seq_along(doses) - 1)*0.5, lab = doses)
    ggplot(d, aes(time, value)) +
      geom_vline(data = vl, aes(xintercept = time), colour = "#cccccc",
                 linetype = 3) +
      geom_line(colour = "#b3352f", linewidth = .85) +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(title = "PRACTALL escalation challenge in this patient",
           subtitle = "single doses 3 / 10 / 30 / 100 / 300 / 600 / 1000 mg at 30-min intervals",
           x = "time (h)", y = NULL) + theme_fa
  })

  output$t_dbpcfc <- renderDT({
    doses <- c(3, 10, 30, 100, 300, 600, 1000)
    d <- do.call(rbind, lapply(doses, function(z) {
      s <- mod %>% param(par_now()) %>% ev(ev_allergen(z)) %>%
        mrgsim(end = 6, delta = 0.02) %>% as.data.frame()
      data.frame(`dose (mg)` = z,
                 `peak skin` = round(max(s$URTskin), 2),
                 `peak GI` = round(max(s$GISY), 2),
                 `FEV1 nadir` = round(min(s$FEV1), 1),
                 `MAP nadir` = round(min(s$MAP), 1),
                 `peak tryptase` = round(max(s$Tryptase), 1),
                 `grade` = as.integer(max(s$SEV)),
                 check.names = FALSE)
    }))
    datatable(d, rownames = FALSE, options = list(dom = "t", pageLength = 10))
  })

  ## ---- tab 4 --------------------------------------------------------
  acute <- reactive({
    p <- par_now()
    e <- ev_allergen(dose_mg())
    if (input$premed) e <- ev_cetirizine(10, time = -2) + e
    if (input$epidelay < 60) e <- e + ev_epi_im(0.3, time = input$epidelay/60)
    mod %>% param(p) %>% ev(e) %>%
      mrgsim(start = -2, end = 8, delta = 0.004) %>% as.data.frame() %>%
      filter(time >= -0.1)
  })

  output$p_mediators <- renderPlot({
    d <- acute() %>%
      select(time, histamine = HIST, PAF = PAF, cysLT = CYSLT) %>%
      pivot_longer(-time)
    ggplot(d, aes(MIN(time), value, colour = name)) +
      geom_line(linewidth = .95) +
      scale_colour_manual(values = PAL, name = NULL) +
      coord_cartesian(xlim = c(0, 240)) +
      labs(title = sprintf("Mediator waves after %s mg peanut protein",
                           formatC(dose_mg(), format = "f", digits = 0, big.mark = ",")),
           subtitle = "histamine peaks and clears in minutes; cysLT sustains the lesion",
           x = "minutes", y = "pool (1 = full anaphylactic level)") + theme_fa
  })

  output$p_physiology <- renderPlot({
    d <- acute() %>%
      select(time, `MAP (mmHg)` = MAP, `FEV1 (% pred)` = FEV1,
             `plasma volume deficit (%)` = PVdef,
             `tryptase (ng/mL)` = Tryptase,
             `skin score` = URTskin, `severity grade` = SEV) %>%
      mutate(`plasma volume deficit (%)` = 100*`plasma volume deficit (%)`) %>%
      pivot_longer(-time)
    ggplot(d, aes(MIN(time), value)) +
      geom_line(colour = "#b3352f", linewidth = .9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      coord_cartesian(xlim = c(0, 300)) +
      labs(title = "Physiology of the event", x = "minutes", y = NULL,
           subtitle = "MAP holds flat while the reserve is spent, then falls") +
      theme_fa
  })

  output$acute_card <- renderUI({
    s <- acute()
    HTML(sprintf(
      "<b>Event summary</b> &nbsp;|&nbsp; dose %s mg &nbsp;|&nbsp; adrenaline %s<br>
       MAP nadir <b>%.0f</b> mmHg &nbsp; &middot; &nbsp;
       FEV1 nadir <b>%.0f%%</b> &nbsp; &middot; &nbsp;
       peak plasma volume deficit <b>%.1f%%</b> &nbsp; &middot; &nbsp;
       peak tryptase <b>%.1f</b> ng/mL &nbsp; &middot; &nbsp;
       maximum grade <b>%d</b><br>
       harm integral &int;(plasma deficit)dt = <b>%.3f</b> L&middot;h",
      formatC(dose_mg(), format = "f", digits = 0, big.mark = ","),
      ifelse(input$epidelay < 60, paste0("at ", input$epidelay, " min"), "none"),
      min(s$MAP), min(s$FEV1), 100*max(s$PVdef), max(s$Tryptase),
      as.integer(max(s$SEV)), max(s$AUCLEAK)))
  })

  ## ---- tab 5 --------------------------------------------------------
  output$p_epi_pk <- renderPlot({
    im <- mod %>% param(par_now()) %>% ev(ev_epi_im(0.3)) %>%
      mrgsim(end = 1, delta = 0.002) %>% as.data.frame() %>%
      mutate(route = "IM 0.3 mg (thigh)")
    ins <- mod %>% param(par_now()) %>% ev(ev_epi_in(2.0)) %>%
      mrgsim(end = 1, delta = 0.002) %>% as.data.frame() %>%
      mutate(route = "intranasal 2 mg")
    d <- rbind(im, ins)
    ggplot(d, aes(MIN(time), Cepi_ngmL, colour = route)) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c("#b3352f", "#2f6ea8"), name = NULL) +
      labs(title = "Adrenaline pharmacokinetics",
           subtitle = "t½ ~2.5 min: one dose is a short window, which is why a third of episodes need a second",
           x = "minutes", y = "plasma adrenaline (ng/mL)") + theme_fa
  })

  epi_sweep <- reactive({
    do.call(rbind, lapply(c(2, 5, 10, 20, 30, 45, NA), function(dl) {
      s <- sc_accidental(par_now(), dose_mg(), epi_delay_min = dl,
                         supine = as.numeric(input$supine))
      data.frame(delay = ifelse(is.na(dl), 99, dl),
                 label = ifelse(is.na(dl), "none", paste0(dl, " min")),
                 map_nadir = min(s$MAP), fev_nadir = min(s$FEV1),
                 pv_def = 100*max(s$PVdef), harm = max(s$AUCLEAK),
                 grade = max(s$SEV))
    }))
  })

  output$p_epi_timing <- renderPlot({
    d <- epi_sweep() %>%
      select(label, delay, `MAP nadir (mmHg)` = map_nadir,
             `peak PV deficit (%)` = pv_def,
             `harm integral (L*h)` = harm) %>%
      pivot_longer(-c(label, delay))
    ggplot(d, aes(reorder(label, delay), value, group = 1)) +
      geom_line(colour = "#4a6785", linewidth = .9) +
      geom_point(size = 3, colour = "#b3352f") +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "Adrenaline delay is an integral, not a switch",
           subtitle = "the drug does not become less effective; the deficit it must repair becomes larger",
           x = "delay to adrenaline", y = NULL) + theme_fa
  })

  output$t_epi <- renderDT({
    d <- epi_sweep() %>%
      transmute(`delay` = label,
                `MAP nadir` = round(map_nadir, 1),
                `FEV1 nadir` = round(fev_nadir, 1),
                `peak PV deficit (%)` = round(pv_def, 1),
                `harm integral (L*h)` = round(harm, 3),
                `max grade` = as.integer(grade))
    datatable(d, rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- tab 6 --------------------------------------------------------
  oma <- reactive({
    sc_omalizumab(par_now(), weeks = input$oma_wk)
  })

  output$oma_dose_card <- renderUI({
    dt <- omalizumab_dose(input$wt, input$tige)
    HTML(sprintf("<b>Label-derived dose</b><br>%d mg every %d weeks<br>
                  <span style='font-size:11px;color:#888'>the table exists because \
                  the drug must supply enough BINDING SITES for the IgE present</span>",
                 dt$dose_mg, round(dt$interval_h/168)))
  })

  output$oma_free_card <- renderUI({
    d <- oma(); l <- d[nrow(d), ]
    HTML(sprintf("<b>At week %d</b><br>free IgE %.1f ng/mL %s<br>total IgE %.0f IU/mL (%.1f&times;)<br>
                  &rho; %.3f &nbsp; surface %.3f&times; baseline",
                 input$oma_wk, l$FreeIgE_ng,
                 ifelse(l$FreeIgE_ng < 25, "<span style='color:#2f8b52'>(at target)</span>",
                        "<span style='color:#b3352f'>(ABOVE target — expect non-response)</span>"),
                 l$TotIgE_IU, l$TotIgE_IU/d$TotIgE_IU[1],
                 l$RHO, l$SURFrel/d$SURFrel[1]))
  })

  output$p_oma <- renderPlot({
    d <- oma() %>%
      mutate(week = time/168) %>%
      select(week, `free IgE (ng/mL)` = FreeIgE_ng,
             `total IgE (IU/mL)` = TotIgE_IU,
             `FcepsilonRI density ρ` = RHO,
             `surface index (relative)` = SURFrel,
             `omalizumab (µg/mL)` = Oma_ug,
             `predicted ED multiplier` = SURFrel) %>%
      mutate(`predicted ED multiplier` =
               (`surface index (relative)`[1]/`surface index (relative)`)^2) %>%
      pivot_longer(-week)
    ggplot(d, aes(week, value)) +
      geom_line(colour = "#2f8b52", linewidth = .95) +
      geom_vline(xintercept = 16, linetype = 2, colour = "#888888") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "Anti-IgE: the fast arm and the slow arm",
           subtitle = "dashed line = week 16, where the trials place the challenge",
           x = "weeks", y = NULL) + theme_fa
  })

  ## ---- tab 7 --------------------------------------------------------
  oit <- reactive({
    sc_oit(par_now(), maint_mg = max(input$oit_mg, 1),
           build_wk = max(4, round(input$oit_wk/3)),
           maint_wk = input$oit_wk - max(4, round(input$oit_wk/3)),
           off_wk = max(input$off_wk, 1))
  })

  output$p_oit <- renderPlot({
    d <- oit() %>%
      mutate(week = time/168) %>%
      select(week, `sIgG4 (µg/mL)` = IgG4_ug,
             `IgG4 : IgE ratio` = G4E_ratio,
             `allergen interception (×)` = Interc,
             `mast cell anergy` = ANERG,
             `Treg (relative)` = TREG,
             `sIgE (IU/mL)` = sIgE_IU) %>%
      pivot_longer(-week)
    ggplot(d, aes(week, value)) +
      geom_line(colour = "#b5761a", linewidth = .95) +
      geom_vline(xintercept = input$oit_wk, linetype = 2, colour = "#b3352f") +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(title = "Oral immunotherapy: two arms with two time constants",
           subtitle = "red dashed line = therapy stops",
           x = "weeks", y = NULL) + theme_fa
  })

  ## ---- tab 8 --------------------------------------------------------
  ladder <- eventReactive(input$run_ladder, {
    withProgress(message = "Computing the threshold ladder", value = 0, {
      p  <- par_now()
      dt <- omalizumab_dose(p$WT, p$TIGE0)
      out <- list()
      step <- function(nm, ed) { incProgress(1/6); data.frame(arm = nm, ED = ed) }

      ed0 <- find_ED(p, tol = 0.04)
      out[[1]] <- step("untreated", ed0)

      st_oma <- state_at(p, ev_omalizumab(dt$dose_mg, dt$interval_h, 14), 16*168)
      out[[2]] <- step("omalizumab wk16", find_ED(p, init0 = st_oma, tol = 0.04))

      poit <- modifyList(p, list(OITON = 1, OITDOSE = 300))
      st_oit <- state_at(poit, ev_oit_daily(300, 26*7), 26*168)
      out[[3]] <- step("OIT 300 mg/d, 6 mo", find_ED(p, init0 = st_oit, tol = 0.04))

      st_cmb <- state_at(poit, ev_omalizumab(dt$dose_mg, dt$interval_h, 14) +
                           ev_oit_daily(300, 18*7, start = 8*168), 26*168)
      out[[4]] <- step("omalizumab + OIT", find_ED(p, init0 = st_cmb, tol = 0.04))

      st_dup <- state_at(p, ev_dupilumab(300, 2*168, 12), 24*168)
      out[[5]] <- step("dupilumab 24 wk", find_ED(p, init0 = st_dup, tol = 0.04))

      out[[6]] <- step("BTK inhibitor",
                       find_ED(modifyList(p, list(BTKI = 1)), tol = 0.04))
      do.call(rbind, out)
    })
  })

  output$p_ladder <- renderPlot({
    d <- ladder()
    d$logshift <- log10(d$ED/d$ED[d$arm == "untreated"])
    ggplot(d, aes(reorder(arm, ED), ED, fill = arm)) +
      geom_col(width = .6, show.legend = FALSE) +
      geom_hline(yintercept = 600, linetype = 2, colour = "#b3352f") +
      annotate("text", x = .7, y = 700, label = "trial endpoint: 600 mg",
               hjust = 0, colour = "#b3352f", size = 3.6) +
      geom_text(aes(label = sprintf("%s mg  (%+.2f log)",
                                    formatC(ED, format = "f", digits = 0, big.mark = ","),
                                    logshift)),
                hjust = -0.05, size = 3.7) +
      scale_y_log10(expand = expansion(mult = c(0, .45))) +
      coord_flip() +
      scale_fill_manual(values = c("untreated" = "#b3352f",
                                   "omalizumab wk16" = "#2f8b52",
                                   "OIT 300 mg/d, 6 mo" = "#b5761a",
                                   "omalizumab + OIT" = "#5b2f86",
                                   "dupilumab 24 wk" = "#2f6ea8",
                                   "BTK inhibitor" = "#4a6785")) +
      labs(title = "What each intervention is worth, in logs of threshold",
           subtitle = "surface-acting drugs buy two logs per log; allergen-acting drugs buy one",
           x = NULL, y = "eliciting dose (mg peanut protein, log scale)") +
      theme_fa
  })

  output$t_ladder <- renderDT({
    d <- ladder()
    d$`fold vs untreated` <- round(d$ED/d$ED[d$arm == "untreated"], 1)
    d$`log10 shift` <- round(log10(d$ED/d$ED[d$arm == "untreated"]), 2)
    d$ED <- round(d$ED, 1)
    datatable(d, rownames = FALSE, options = list(dom = "t"))
  })

  ## ---- tab 9 --------------------------------------------------------
  cof <- eventReactive(input$run_cof, {
    withProgress(message = "Sweeping cofactors", value = 0, {
      base <- patient(WT = input$wt, AGEY = input$agey, SIGE0 = input$sige,
                      TIGE0 = input$tige, MCBURDEN = input$mcb,
                      RELEASE = input$rel)
      combos <- list("none" = list(),
                     "exercise" = list(COF_EX = 1),
                     "NSAID" = list(COF_NSAID = 1),
                     "alcohol" = list(COF_ETOH = 1),
                     "infection" = list(COF_INF = 1),
                     "PPI" = list(COF_PPI = 1),
                     "exercise + NSAID" = list(COF_EX = 1, COF_NSAID = 1),
                     "exercise + alcohol" = list(COF_EX = 1, COF_ETOH = 1))
      do.call(rbind, lapply(names(combos), function(nm) {
        incProgress(1/length(combos))
        data.frame(cofactor = nm,
                   ED = find_ED(modifyList(base, combos[[nm]]), tol = 0.04))
      }))
    })
  })

  output$p_cof <- renderPlot({
    d <- cof()
    d$shift <- log10(d$ED/d$ED[d$cofactor == "none"])
    ggplot(d, aes(reorder(cofactor, -shift), shift, fill = shift)) +
      geom_col(width = .6, show.legend = FALSE) +
      scale_fill_gradient(low = "#b3352f", high = "#dbe8f6") +
      geom_text(aes(label = sprintf("%s mg", formatC(ED, format = "f", digits = 1))),
                vjust = 1.3, size = 3.6) +
      coord_flip() +
      labs(title = "Cofactors move the threshold; they do not change the antibodies",
           subtitle = "within-patient variation across repeat challenges is ~0.4-0.5 log10 — this is where it comes from",
           x = NULL, y = "log10 shift in eliciting dose") + theme_fa
  })

  output$p_cof_event <- renderPlot({
    d <- sc_cofactor(par_now(), dose_mg()) %>%
      select(time, arm, `MAP (mmHg)` = MAP, `FEV1 (% pred)` = FEV1,
             `skin score` = URTskin, `severity grade` = SEV) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(MIN(time), value, colour = arm)) +
      geom_line(linewidth = .9) +
      scale_colour_manual(values = PAL, name = NULL) +
      facet_wrap(~name, scales = "free_y") +
      coord_cartesian(xlim = c(0, 240)) +
      labs(title = sprintf("The same %s mg, on two different days",
                           formatC(dose_mg(), format = "f", digits = 0, big.mark = ",")),
           x = "minutes", y = NULL) + theme_fa
  })
}

shinyApp(ui, server)
