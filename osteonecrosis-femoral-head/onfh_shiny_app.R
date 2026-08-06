## =============================================================================
##  ONFH QSP model -- Shiny dashboard
##  Osteonecrosis of the Femoral Head
## =============================================================================
##
##  Run with:
##      install.packages(c("shiny", "mrgsolve", "dplyr", "tidyr", "ggplot2", "DT"))
##      shiny::runApp("onfh_shiny_app.R")
##
##  The app is built around ONE question, asked twelve ways:
##
##      why does a femoral head that is fully mineralised, as stiff as it ever
##      was, and carrying a perfectly ordinary load, break?
##
##  The answer the model gives is that it breaks at the CONICAL INTERFACE, and
##  it breaks there BECAUSE IT IS HEALING.  Tab 4 shows the interface losing
##  more than half of its strength by month 4-6 (modulus down to 35% of normal)
##  and getting all of it back by year 3; tab 6 shows that the microdamage clock is running fastest exactly
##  then; tab 9 shows that where the transition sits as a function of lesion
##  geometry is the Kerboull angle, which the model was never shown.
##
##  Tabs
##    1  Patient and lesion       geometry, JIC type, the two areas
##    2  The three clocks         perfusion / repair / fatigue on one time axis
##    3  Perfusion and compartment  the closed-compartment Starling resistor
##    4  The reparative front      the vulnerable window
##    5  Mechanics                 stress and strength for all five probes
##    6  Damage and collapse       which structure fails, and when
##    7  Clinical endpoints        VAS, Harris Hip Score, arthroplasty risk
##    8  Treatment comparison      ten arms on a matched hip
##    9  Staging explorer          collapse map over CNA x location
##   10  Virtual cohort            collapse-free survival by JIC class
##   11  Sensitivity               local sensitivity of 5-year depression
##   12  Model card                what is fitted, what is predicted, misses
## =============================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

MODEL_FILE <- "onfh_mrgsolve_model.R"
mod <- mread(MODEL_FILE)

DAY_PER_MO <- 30.44
TEND <- 1826

## ---- JIC presets (mid-class by construction) --------------------------------
JIC_PRESETS <- list(
  "A  (medial third)"                    = c(alpha = 28, phi = -55),
  "B  (medial two thirds)"               = c(alpha = 34, phi = -34),
  "C1 (beyond two thirds, inside rim)"   = c(alpha = 45, phi = -12),
  "C2 (beyond the acetabular edge)"      = c(alpha = 58, phi =   4)
)

jic_type <- function(theta_L, theta_c = 50) {
  ifelse(theta_L <= -theta_c / 3, "A",
    ifelse(theta_L <= theta_c / 3, "B",
      ifelse(theta_L <= theta_c, "C1", "C2")))
}

## ---- treatment arms ---------------------------------------------------------
ARMS <- list(
  "untreated" = list(par = list(), ev = NULL),
  "alendronate 70 mg weekly x 24 mo" =
    list(par = list(), ev = ev(amt = 16.0, cmt = "ALN_c", ii = 7, addl = 103)),
  "denosumab 60 mg q6m x 24 mo" =
    list(par = list(), ev = ev(amt = 18.75, cmt = "DMB_c", ii = 182.5, addl = 3)),
  "core decompression (single track)" =
    list(par = list(CD_time = 30, CD_ntrack = 1, PWB = 0.55), ev = NULL),
  "multiple drilling + BMAC" =
    list(par = list(CD_time = 30, CD_ntrack = 3, BMAC_dose = 2, PWB = 0.55),
         ev = NULL),
  "teriparatide 20 ug/d x 18 mo" =
    list(par = list(), ev = ev(amt = 20, cmt = "TPT_sc", ii = 1, addl = 547)),
  "alendronate 12 mo -> teriparatide 18 mo" =
    list(par = list(),
         ev = seq(ev(amt = 16.0, cmt = "ALN_c", ii = 7, addl = 51),
                  ev(amt = 20, cmt = "TPT_sc", ii = 1, addl = 547))),
  "iloprost 5-day infusion" =
    list(par = list(), ev = ev(amt = 1.2, cmt = "ILO", ii = 1, addl = 4, time = 14)),
  "structural graft + BMAC" =
    list(par = list(CD_time = 30, CD_ntrack = 3, BMAC_dose = 2,
                    SUPPORT_GRAFT = 1, PWB = 0.50), ev = NULL),
  "protected weight bearing alone" =
    list(par = list(PWB = 0.45), ev = NULL)
)

run_arm <- function(alpha, phi, base_par, arm, tend = TEND) {
  m <- mod %>% param(alpha_deg = alpha, phi_deg = phi)
  if (length(base_par)) m <- do.call(param, c(list(m), base_par))
  if (length(arm$par))  m <- do.call(param, c(list(m), arm$par))
  if (is.null(arm$ev)) {
    out <- m %>% mrgsim(end = tend, delta = 2)
  } else {
    out <- m %>% mrgsim(arm$ev, end = tend, delta = 2)
  }
  as_tibble(out)
}

collapse_month <- function(df, thr = 2) {
  i <- which(df$DEPR >= thr)
  if (!length(i)) return(NA_real_)
  df$time[i[1]] / DAY_PER_MO
}

theme_onfh <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("Osteonecrosis of the Femoral Head — QSP model"),
  tags$p(tags$em(paste(
    "Three clocks run in a closed compartment: perfusion (hours), repair",
    "(months), fatigue (years). The hip collapses because it started to heal."))),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Lesion"),
      selectInput("preset", "JIC preset", names(JIC_PRESETS), selected = names(JIC_PRESETS)[3]),
      sliderInput("alpha", "cone half-angle α (deg)  →  CNA = 4α",
                  10, 78, 45, step = 1),
      sliderInput("phi", "axis tilt φ (deg, + = lateral)", -70, 34, -12, step = 1),
      htmlOutput("geomBadge"),
      hr(),
      h4("Patient"),
      sliderInput("bw", "body weight (kg)", 44, 115, 70, step = 1),
      sliderInput("ncyc", "gait cycles per hip per day", 1000, 15000, 5000, step = 250),
      sliderInput("str", "trabecular strength S₀ (MPa)", 6, 14, 9.5, step = 0.1),
      sliderInput("xf0", "front penetration at presentation (mm)", 0.3, 8, 0.5, step = 0.1),
      hr(),
      h4("Repair biology"),
      sliderInput("vfront", "creeping substitution speed (mm/month)",
                  0.15, 1.8, 0.75, step = 0.05),
      sliderInput("respmax", "maximum interface turnover", 0.25, 0.75, 0.55, step = 0.01),
      sliderInput("kfill", "refill rate multiplier", 0.4, 2.0, 1.0, step = 0.05),
      hr(),
      selectInput("arm", "treatment arm", names(ARMS)),
      helpText("Everything on this page is recomputed from the",
               "49-ODE model; nothing is tabulated.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",
        tabPanel("1 · Patient & lesion",
                 plotOutput("headPlot", height = "380px"),
                 h4("What the geometry alone already decides"),
                 tableOutput("geomTable"),
                 helpText(HTML(paste0(
                   "The loaded cap grows like α² and the conical interface only ",
                   "like α, so <b>A_cap/A_int = 2 tan(α/2)</b> rises without ",
                   "bound. Every newton landing on the necrotic cap has to ",
                   "leave through that cone, because the apex is a point.")))),
        tabPanel("2 · The three clocks", plotOutput("clocks", height = "620px")),
        tabPanel("3 · Perfusion & compartment", plotOutput("perf", height = "620px")),
        tabPanel("4 · The reparative front", plotOutput("front", height = "620px"),
                 htmlOutput("windowText")),
        tabPanel("5 · Mechanics", plotOutput("mech", height = "620px"),
                 tableOutput("mechTable")),
        tabPanel("6 · Damage & collapse", plotOutput("dmg", height = "620px"),
                 htmlOutput("failMode")),
        tabPanel("7 · Clinical endpoints", plotOutput("clin", height = "620px"),
                 tableOutput("clinTable")),
        tabPanel("8 · Treatment comparison",
                 plotOutput("arms", height = "460px"), DTOutput("armTable")),
        tabPanel("9 · Staging explorer",
                 plotOutput("stage", height = "560px"),
                 helpText(paste(
                   "Collapse map over lesion size and position. The white",
                   "contour is the model's collapse boundary; Kerboull's",
                   "published risk bands (<190 low, 190-240 intermediate,",
                   ">240 high) are drawn for comparison. The model was never",
                   "shown them."))),
        tabPanel("10 · Virtual cohort",
                 numericInput("ncoh", "cohort size", 120, 40, 400, 20),
                 actionButton("runcoh", "simulate cohort"),
                 plotOutput("km", height = "480px"), tableOutput("kmTable")),
        tabPanel("11 · Sensitivity", plotOutput("sens", height = "560px")),
        tabPanel("12 · Model card", htmlOutput("card"))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  observeEvent(input$preset, {
    p <- JIC_PRESETS[[input$preset]]
    updateSliderInput(session, "alpha", value = p[["alpha"]])
    updateSliderInput(session, "phi",   value = p[["phi"]])
  })

  base_par <- reactive(list(
    BW = input$bw, Ncyc = input$ncyc, S_tr0 = input$str, XF0 = input$xf0,
    v_front = input$vfront / DAY_PER_MO, RESP_max = input$respmax,
    k_fill = 0.030 * input$kfill
  ))

  sim <- reactive(run_arm(input$alpha, input$phi, base_par(), ARMS[[input$arm]]))
  ctl <- reactive(run_arm(input$alpha, input$phi, base_par(), ARMS[["untreated"]]))

  theta_L <- reactive(input$alpha + input$phi)

  output$geomBadge <- renderUI({
    tl <- theta_L()
    HTML(sprintf(
      "<div style='padding:6px;border-radius:6px;background:#eef3ff'>
       <b>JIC type %s</b> &nbsp;·&nbsp; θ_L = %.0f°
       &nbsp;·&nbsp; Kerboull CNA = %.0f°</div>",
      jic_type(tl), tl, 4 * input$alpha))
  })

  ## ---- 1 head diagram -------------------------------------------------------
  output$headPlot <- renderPlot({
    a <- input$alpha; ph <- input$phi; R <- 22.5; tc <- 50
    circ <- tibble(t = seq(0, 2 * pi, length.out = 400)) %>%
      mutate(x = R * sin(t), y = R * cos(t))
    wedge <- tibble(t = seq((ph - a) * pi / 180, (ph + a) * pi / 180,
                            length.out = 200)) %>%
      mutate(x = R * sin(t), y = R * cos(t)) %>%
      bind_rows(tibble(x = 0, y = 0), .)
    contact <- tibble(t = seq(-tc * pi / 180, tc * pi / 180, length.out = 200)) %>%
      mutate(x = R * 1.06 * sin(t), y = R * 1.06 * cos(t))
    press <- tibble(t = seq(-tc * pi / 180, tc * pi / 180, length.out = 200)) %>%
      mutate(p = cos(t),
             x = (R * 1.10 + 6 * p) * sin(t), y = (R * 1.10 + 6 * p) * cos(t))
    ggplot() +
      geom_polygon(data = wedge, aes(x, y), fill = "#e6b8ae", colour = NA) +
      geom_path(data = circ, aes(x, y), linewidth = 1.1, colour = "#333333") +
      geom_path(data = contact, aes(x, y), linewidth = 2.4, colour = "#6f9dc9",
                alpha = .8) +
      geom_path(data = press, aes(x, y), linewidth = 0.8, colour = "#c98f2e") +
      geom_segment(aes(x = 0, y = 0, xend = 0, yend = R * 1.30),
                   arrow = arrow(length = unit(3, "mm"), ends = "first"),
                   linewidth = 1.0, colour = "#8f2f22") +
      annotate("text", x = 0, y = R * 1.42, label = "hip joint resultant force",
               colour = "#8f2f22", size = 3.6) +
      annotate("text", x = R * 0.95, y = -R * 0.3, label = "lateral", size = 3.4) +
      annotate("text", x = -R * 0.95, y = -R * 0.3, label = "medial", size = 3.4) +
      coord_fixed(xlim = c(-R * 1.5, R * 1.5), ylim = c(-R * 1.2, R * 1.6)) +
      labs(title = "Coronal section: the necrotic cone in the contact field",
           subtitle = paste0("shaded = necrotic sector (2α = ", 2 * input$alpha,
                             "°); blue arc = acetabular contact patch; ",
                             "gold curve = p(θ) = p₀cos θ"),
           x = NULL, y = NULL) +
      theme_void(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))
  })

  output$geomTable <- renderTable({
    d <- sim()[1, ]
    tibble(quantity = c("peak contact pressure p₀ (MPa)",
                        "loaded cap area A_cap (mm²)",
                        "conical interface area A_int (mm²)",
                        "geometric amplification A_cap/A_int = 2 tan(α/2)",
                        "fraction of joint load on the lesion",
                        "lateral boundary θ_L (deg)",
                        "JIC type", "Kerboull CNA (deg)"),
           value = c(sprintf("%.3f", d$P0), sprintf("%.0f", d$A_CAP),
                     sprintf("%.0f", d$A_INT), sprintf("%.3f", d$AMPLIF),
                     sprintf("%.3f", d$F_LOAD), sprintf("%.1f", d$THETA_L),
                     jic_type(d$THETA_L), sprintf("%.0f", d$CNA)))
  }, striped = TRUE)

  ## ---- 2 three clocks -------------------------------------------------------
  output$clocks <- renderPlot({
    d <- sim() %>%
      transmute(month = time / DAY_PER_MO,
                `clock 1 · intraosseous pressure (mmHg)` = PIO,
                `clock 1 · marrow adipocyte fraction` = ADIPO,
                `clock 2 · front penetration XF (mm)` = XF,
                `clock 2 · interface bone volume fraction β` = BETA_EFF,
                `clock 3 · interface stress / strength` = SR_INT,
                `clock 3 · head depression (mm)` = DEPR) %>%
      pivot_longer(-month, names_to = "v", values_to = "y")
    d$v <- factor(d$v, levels = unique(d$v))
    ggplot(d, aes(month, y)) +
      geom_line(linewidth = 1.0, colour = "#2b5f8e") +
      facet_wrap(~ v, scales = "free_y", ncol = 2) +
      labs(title = "The three clocks", x = "months since the lesion formed",
           y = NULL) + theme_onfh
  })

  ## ---- 3 perfusion ----------------------------------------------------------
  output$perf <- renderPlot({
    sim() %>%
      transmute(month = time / DAY_PER_MO,
                `intraosseous pressure (mmHg)` = PIO,
                `marrow adipocyte fraction` = ADIPO,
                `microvascular occlusion` = THROMB,
                `marrow oedema fraction` = EDEMA,
                `endothelial competence` = ENDO,
                `PAI-1 (fold)` = PAI1) %>%
      pivot_longer(-month) %>%
      ggplot(aes(month, value)) + geom_line(linewidth = 1, colour = "#8a5fa8") +
      facet_wrap(~ name, scales = "free_y") +
      labs(title = "Clock 1 — the closed compartment",
           subtitle = paste("Perfusion is a Starling resistor:",
                            "Q = (P_a − max(P_v, P_io)) / R.",
                            "Once P_io exceeds venous pressure the head is",
                            "perfused by the difference, not by the arterial",
                            "pressure."),
           x = "months", y = NULL) + theme_onfh
  })

  ## ---- 4 front --------------------------------------------------------------
  output$front <- renderPlot({
    sim() %>%
      transmute(month = time / DAY_PER_MO,
                `front penetration XF (mm)` = XF,
                `open resorption cavity CAV` = CAV,
                `new bone NB` = NB,
                `mineralisation MINZ` = MINZ,
                `interface bone volume fraction β` = BETA_EFF,
                `interface modulus E_int (MPa)` = E_INT) %>%
      pivot_longer(-month) %>%
      ggplot(aes(month, value)) + geom_line(linewidth = 1, colour = "#3f7a52") +
      facet_wrap(~ name, scales = "free_y") +
      labs(title = "Clock 2 — creeping substitution, and why it weakens the hip",
           subtitle = paste("Osteoclasts excavate before osteoblasts refill,",
                            "and the moving front keeps resetting the active",
                            "zone to mid-cycle, so the weakness persists for as",
                            "long as the front is still crossing."),
           x = "months", y = NULL) + theme_onfh
  })

  output$windowText <- renderUI({
    d <- sim()
    i <- which.min(d$BETA_EFF)
    HTML(sprintf(paste0(
      "<div style='padding:8px;background:#f0f9f2;border-radius:6px'>",
      "Interface strength bottoms out at <b>month %.1f</b> ",
      "(β = %.3f, modulus %.0f MPa = %.0f%% of normal, strength %.2f MPa). ",
      "It is back above 0.95 by month %s. ",
      "<b>The necrotic bone never changes.</b> Everything that moves in this ",
      "panel is the repair.</div>"),
      d$time[i] / DAY_PER_MO, d$BETA_EFF[i], d$E_INT[i],
      100 * d$E_INT[i] / 620, d$S_INT[i],
      { j <- which(d$BETA_EFF > 0.95 & d$time > d$time[i])
        if (length(j)) sprintf("%.0f", d$time[j[1]] / DAY_PER_MO) else "—" }))
  })

  ## ---- 5 mechanics ----------------------------------------------------------
  output$mech <- renderPlot({
    sim() %>%
      transmute(month = time / DAY_PER_MO,
                `interface traction σ_int (MPa)` = SIGMA_INT,
                `interface strength S_int (MPa)` = S_INT,
                `stress / strength at the interface` = SR_INT) %>%
      pivot_longer(-month) %>%
      ggplot(aes(month, value)) + geom_line(linewidth = 1.1, colour = "#a86f14") +
      facet_wrap(~ name, scales = "free_y", ncol = 1) +
      labs(title = "Stress is nearly constant; strength is what moves",
           subtitle = paste("The traction is set once, by geometry and body",
                            "weight. The ratio rises only because the",
                            "denominator falls."),
           x = "months", y = NULL) + theme_onfh
  })

  output$mechTable <- renderTable({
    pre <- sim() %>% filter(time <= 550)
    tibble(probe = c("reparative interface", "untouched necrotic core",
                     "subchondral plate", "acetabular rim",
                     "adjacent living bone (negative control)"),
           `peak stress / strength, first 18 months` = sprintf(
             "%.3f", c(max(pre$SR_INT), max(pre$SR_NEC), max(pre$SR_PL),
                       max(pre$SR_RIM), max(pre$SR_LIV))),
           `damage cleared by remodelling?` =
             c("yes, but the tissue is weak", "NO — no osteocytes to signal it",
               "partly", "only if the rim sits on living bone", "yes"))
  }, striped = TRUE)

  ## ---- 6 damage -------------------------------------------------------------
  output$dmg <- renderPlot({
    sim() %>%
      transmute(month = time / DAY_PER_MO,
                `microdamage at the interface` = pmin(D_int, 1.5),
                `microdamage in necrotic bone` = pmin(D_nec, 1.5),
                `microdamage in living bone` = pmin(D_liv, 1.5),
                `crescent extent` = CRESC,
                `head depression (mm)` = DEPR,
                `stress / strength` = SR_INT) %>%
      pivot_longer(-month) %>%
      ggplot(aes(month, value)) + geom_line(linewidth = 1, colour = "#8f2f22") +
      geom_hline(data = function(z) subset(z, grepl("microdamage", z$name)),
                 aes(yintercept = 1), linetype = 2, colour = "grey40") +
      facet_wrap(~ name, scales = "free_y") +
      labs(title = "Clock 3 — fatigue",
           subtitle = paste("Damage in LIVING bone reaches a steady state,",
                            "because targeted remodelling removes it. Damage in",
                            "NECROTIC bone has no removal term at all — but the",
                            "necrotic bone is strong, so it never gets there",
                            "first."),
           x = "months", y = NULL) + theme_onfh
  })

  output$failMode <- renderUI({
    d <- sim(); m <- collapse_month(d)
    HTML(sprintf(
      "<div style='padding:8px;background:#fdf0ee;border-radius:6px'>%s</div>",
      if (is.na(m)) "This hip does not reach 2 mm of depression within 5 years."
      else sprintf(paste("Collapse (2 mm depression) at <b>month %.1f</b>.",
                         "Depression at 5 years %.2f mm; Harris Hip Score %.0f."),
                   m, tail(d$DEPR, 1), tail(d$HHS, 1))))
  })

  ## ---- 7 clinical -----------------------------------------------------------
  output$clin <- renderPlot({
    bind_rows(mutate(sim(), arm = input$arm), mutate(ctl(), arm = "untreated")) %>%
      transmute(arm, month = time / DAY_PER_MO,
                `pain VAS (0-10)` = VAS,
                `Harris Hip Score` = HHS,
                `head depression (mm)` = DEPR,
                `cumulative P(arthroplasty)` = P_THA) %>%
      pivot_longer(-c(arm, month)) %>%
      ggplot(aes(month, value, colour = arm)) + geom_line(linewidth = 1) +
      facet_wrap(~ name, scales = "free_y") +
      scale_colour_manual(values = c("#2b5f8e", "#c0392b")) +
      labs(title = "Clinical endpoints", x = "months", y = NULL, colour = NULL) +
      theme_onfh
  })

  output$clinTable <- renderTable({
    a <- sim(); b <- ctl()
    tibble(endpoint = c("collapse (month)", "depression at 5 y (mm)",
                        "Harris Hip Score at 5 y", "P(arthroplasty) at 5 y"),
           `selected arm` = c(sprintf("%.1f", collapse_month(a)),
                              sprintf("%.2f", tail(a$DEPR, 1)),
                              sprintf("%.0f", tail(a$HHS, 1)),
                              sprintf("%.3f", tail(a$P_THA, 1))),
           untreated = c(sprintf("%.1f", collapse_month(b)),
                         sprintf("%.2f", tail(b$DEPR, 1)),
                         sprintf("%.0f", tail(b$HHS, 1)),
                         sprintf("%.3f", tail(b$P_THA, 1))))
  }, striped = TRUE)

  ## ---- 8 arms ---------------------------------------------------------------
  armRuns <- reactive({
    lapply(names(ARMS), function(nm)
      mutate(run_arm(input$alpha, input$phi, base_par(), ARMS[[nm]]), arm = nm)) %>%
      bind_rows()
  })

  output$arms <- renderPlot({
    armRuns() %>%
      ggplot(aes(time / DAY_PER_MO, DEPR, colour = arm)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 2, linetype = 2) +
      annotate("text", x = 3, y = 2.25, hjust = 0, size = 3.4,
               label = "radiographic collapse (ARCO IIIA/IIIB)") +
      labs(title = "Ten arms on the same hip",
           x = "months", y = "femoral head depression (mm)", colour = NULL) +
      theme_onfh + theme(legend.position = "right")
  })

  output$armTable <- renderDT({
    ar <- armRuns()
    lapply(unique(ar$arm), function(nm) {
      d <- ar[ar$arm == nm, ]
      tibble(arm = nm,
             `collapse (month)` = round(collapse_month(d), 1),
             `depression 5 y (mm)` = round(last(d$DEPR), 2),
             `minimum interface β` = round(min(d$BETA_EFF), 3),
             `peak stress/strength` = round(max(d$SR_INT[d$time <= 550]), 3),
             `HHS 5 y` = round(last(d$HHS)),
             `P(THA) 5 y` = round(last(d$P_THA), 3))
    }) %>% bind_rows() %>%
      datatable(options = list(pageLength = 12, dom = "t"), rownames = FALSE)
  })

  ## ---- 9 staging map --------------------------------------------------------
  output$stage <- renderPlot({
    cna <- seq(80, 300, by = 20)
    tl  <- seq(-40, 70, by = 10)
    grid <- expand.grid(CNA = cna, theta_L = tl) %>%
      mutate(alpha = CNA / 4, phi = theta_L - alpha) %>%
      filter(phi >= -85, theta_L - 2 * alpha >= -115)
    res <- lapply(seq_len(nrow(grid)), function(i) {
      d <- run_arm(grid$alpha[i], grid$phi[i], base_par(), ARMS[[input$arm]])
      tibble(CNA = grid$CNA[i], theta_L = grid$theta_L[i],
             tcol = collapse_month(d), depr = last(d$DEPR))
    }) %>% bind_rows()
    ggplot(res, aes(CNA, theta_L, fill = pmin(depr, 8))) +
      geom_tile() +
      scale_fill_gradientn(colours = c("#f7fbff", "#c6dbef", "#fdd0a2",
                                       "#fc9272", "#a50f15"),
                           name = "depression\nat 5 y (mm)") +
      geom_vline(xintercept = c(190, 240), linetype = 2, colour = "white") +
      annotate("text", x = 190, y = max(tl), label = "Kerboull 190°",
               angle = 90, vjust = -0.4, hjust = 1, colour = "white", size = 3.2) +
      annotate("text", x = 240, y = max(tl), label = "240°",
               angle = 90, vjust = -0.4, hjust = 1, colour = "white", size = 3.2) +
      geom_hline(yintercept = c(-50/3, 50/3, 50), linetype = 3) +
      labs(title = "Collapse map: lesion size against lesion position",
           subtitle = paste("Horizontal lines are the JIC A/B, B/C1 and C1/C2",
                            "boundaries, which are geometric definitions, not",
                            "model parameters."),
           x = "Kerboull combined necrotic angle (deg)",
           y = "lateral lesion boundary θ_L (deg)") + theme_onfh
  })

  ## ---- 10 cohort ------------------------------------------------------------
  cohort <- eventReactive(input$runcoh, {
    set.seed(20260806)
    n <- input$ncoh
    alpha <- pmin(pmax(rnorm(n, 38, 14), 8), 78)
    phi   <- pmin(pmax(rnorm(n, -14, 20), -70), 34)
    bw    <- pmin(pmax(rnorm(n, 70, 11), 44), 112)
    ncy   <- pmin(pmax(rlnorm(n, log(5000), 0.28), 900), 18000)
    s0    <- pmin(pmax(rnorm(n, 9.5, 0.95), 6), 13)
    lapply(seq_len(n), function(i) {
      bp <- list(BW = bw[i], Ncyc = ncy[i], S_tr0 = s0[i],
                 XF0 = runif(1, 0.3, 0.45 * max(2, 22.5 * sin(alpha[i] * pi / 180))))
      d <- run_arm(alpha[i], phi[i], bp, ARMS[[input$arm]])
      tibble(jic = jic_type(alpha[i] + phi[i]), tcol = collapse_month(d))
    }) %>% bind_rows()
  })

  output$km <- renderPlot({
    ch <- cohort()
    grid <- tibble(month = seq(0, 60, 1))
    km <- ch %>% group_by(jic) %>% group_modify(function(d, k) {
      grid %>% mutate(S = sapply(month, function(m)
        mean(is.na(d$tcol) | d$tcol > m)))
    }) %>% ungroup()
    ggplot(km, aes(month, S, colour = jic)) + geom_step(linewidth = 1.1) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = "Collapse-free survival by JIC class",
           subtitle = paste("Only the C1 rate was used to fit the single",
                            "damage constant; A, B and C2 are predictions."),
           x = "months", y = "collapse-free survival", colour = "JIC") +
      theme_onfh
  })

  output$kmTable <- renderTable({
    cohort() %>% group_by(jic) %>%
      summarise(n = n(),
                `collapse by 2 y` = mean(!is.na(tcol) & tcol <= 24),
                `by 3 y` = mean(!is.na(tcol) & tcol <= 36),
                `by 5 y` = mean(!is.na(tcol) & tcol <= 60), .groups = "drop")
  }, digits = 3, striped = TRUE)

  ## ---- 11 sensitivity -------------------------------------------------------
  output$sens <- renderPlot({
    keys <- c("k_dmg", "m_fat", "RESP_max", "v_front", "k_fill", "k_minz",
              "S_tr0", "E_tr0", "Ncyc", "f_hip", "BW", "k_resp", "theta_c",
              "zeta_max", "chi_rim0", "f_spread", "L_eng", "W_front")
    base <- run_arm(input$alpha, input$phi, base_par(), ARMS[[input$arm]])
    ref <- last(base$DEPR)
    res <- lapply(keys, function(k) {
      v0 <- as.numeric(as.list(param(mod))[[k]])
      lo <- run_arm(input$alpha, input$phi,
                    modifyList(base_par(), setNames(list(v0 * 0.8), k)),
                    ARMS[[input$arm]])
      hi <- run_arm(input$alpha, input$phi,
                    modifyList(base_par(), setNames(list(v0 * 1.2), k)),
                    ARMS[[input$arm]])
      tibble(param = k,
             sens = (last(hi$DEPR) - last(lo$DEPR)) / (0.4 * max(ref, 1e-6)))
    }) %>% bind_rows()
    ggplot(res, aes(reorder(param, abs(sens)), sens)) +
      geom_col(fill = "#2b5f8e") + coord_flip() +
      labs(title = "Local sensitivity of 5-year head depression (±20%)",
           x = NULL, y = "normalised sensitivity") + theme_onfh
  })

  ## ---- 12 model card --------------------------------------------------------
  output$card <- renderUI(HTML("
<h3>What is fitted</h3>
<p>Three numbers. <code>k_dmg</code>, the microdamage rate constant, anchored on
ONE quantity: the five-year collapse rate of untreated JIC C1 hips.
<code>h0_tha</code>, the baseline arthroplasty hazard. <code>k_pain</code>, a
pain scale. Nothing else in the model was tuned against an outcome.</p>

<h3>What is predicted</h3>
<ul>
<li>Collapse rates for JIC types A, B and C2.</li>
<li>The Kerboull combined-necrotic-angle threshold. The model contains no
    number resembling 190°; it contains a cone, a pressure field, and a
    power-law S-N curve.</li>
<li>The <b>timing</b> of the hazard: the interface loses strength at months
    4-12 and gets it back by year 3, so the collapse hazard is not monotone.
    It peaks early and then vanishes — which is why a hip still intact at four
    years almost never collapses afterwards.</li>
<li>The ceiling on bisphosphonate benefit. Alendronate acts on one limb of the
    interface and cannot make bone; the coupling term means it also removes
    part of the stimulus for the formation it is protecting.</li>
<li>The sign flip of core decompression with lesion size: venting and channel
    creation do not scale with the lesion, but the bone the track removes
    does.</li>
<li>The dominance of peak steroid dose over cumulative steroid dose, which
    follows from a single Hill coefficient of 2.5 on the adipogenic drive.</li>
</ul>

<h3>Known misses, stated rather than fitted away</h3>
<ul>
<li>Type B collapses far too rarely &mdash; 2.5% against a reported 15&ndash;25%.
    This is the model's largest quantitative failure. Real JIC typing has
    measurement error and real lesions are not cones; both would move mass from
    B into C1 in the published series.</li>
<li>The Kerboull threshold lands near 150&deg;, not 190&deg;. The gradient
    across his risk bands is right; the intercept is not.</li>
<li>Protected weight bearing alone abolishes collapse here, and clinically it
    does not. Gait cycles enter the damage law linearly and the model has no
    notion of compliance, muscle loading, or the fact that nobody stays
    non-weight-bearing for the 6&ndash;12 months the window lasts.</li>
<li>The two randomised alendronate trials disagree with each other, and the
    model cannot reproduce both. See README.md for what it does say about the
    maximum effect size the mechanism allows.</li>
<li>Damage is tracked at five fixed probes, not as a field. A real crescent
    propagates.</li>
</ul>

<h3>Not for clinical use</h3>
<p>This is a research and teaching model. It has not been validated against
patient-level data and must not inform care.</p>"))
}

shinyApp(ui, server)
