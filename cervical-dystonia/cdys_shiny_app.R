## =====================================================================
##  cdys_shiny_app.R
##  Cervical dystonia QSP model -- interactive dashboard
##
##      library(shiny); library(mrgsolve)
##      shiny::runApp("cdys_shiny_app.R")
##
##  TAB 2 IS THE POINT OF THIS APP.
##  Everything else is conventional PK/PD plotting.  Tab 2 ("The Bound") draws
##  the one line no botulinum toxin can cross at a given injection geometry,
##  and lets the user try to cross it with dose.  It cannot be done: the line
##  contains no drug parameter.
##
##      L(t) >= L_min = (1 - phi) * D_cen ,    phi = rho * sum(w over injected)
##
##  Ten tabs:
##    1  Patient & plan          -- who is being treated and where the needle goes
##    2  THE BOUND               -- phi, the ceiling, and the four levers
##    3  Toxin PK                -- A -> B -> C, and why the drug is gone by day 2
##    4  Chemodenervation PD     -- SNARE, safety factor, sprouting
##    5  Clinical endpoints      -- TWSTRS subscales against trial anchors
##    6  Dose-response           -- the plateau, and what it measures
##    7  Adverse effects         -- spread, dysphagia, neck weakness, dry mouth
##    8  Chronic cycles          -- the central ratchet in the troughs
##    9  Immunogenicity          -- antibody, non-response, serotype rescue
##   10  Scenario comparison     -- all 15 scenarios side by side
## =====================================================================

library(shiny)
library(mrgsolve)

MODEL_FILE <- "cdys_mrgsolve_model.R"
mod <- mread_cache("cdys", MODEL_FILE)

## $ENV holds the scenario list and helpers defined in the model file.
MENV <- if (exists("env_get", asNamespace("mrgsolve"))) mrgsolve::env_get(mod) else mod@envir

## --- fixed reference values, matching the calibrated model -------------
TW0        <- 42.94
RHO_CAL    <- 0.5463
MCID       <- 4.5
W_DEFAULT  <- c(0.20, 0.24, 0.13, 0.08, 0.12, 0.06, 0.10, 0.07)
MUSCLES    <- c("SCM (contra)", "Splenius capitis", "Trapezius",
                "Levator scapulae", "Semispinalis (deep)", "Scalenes",
                "Obliquus cap. inf. (deep)", "Longus colli (unreachable)")
PROD_NAMES <- c("incobotulinumtoxinA", "onabotulinumtoxinA",
                "abobotulinumtoxinA", "daxibotulinumtoxinA",
                "rimabotulinumtoxinB")

## published anchors, for overlay on tab 5
ANCHORS <- data.frame(week = c(4, 12), dTW = c(-10.5, -6.5))

phi_of <- function(pattern, rho) rho * sum(W_DEFAULT[pattern > 0])

tw_floor <- function(phi) {
  L <- pmax(1 - phi, 0)
  sev <- 35 * L^1.6 / (L^1.6 + 0.8347^1.6)
  a <- 0.03 * L^1.5
  pn <- a / (a + 0.03)
  dis <- 30 * 0.8 * (0.55 * sev / 35 + 0.45 * pn)
  c(severity = sev, pain = 20 * pn, disability = dis,
    total = sev + dis + 20 * pn)
}

## ---------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Cervical Dystonia QSP model — botulinum neurotoxin PK/PD"),
  tags$head(tags$style(HTML("
    .keynum {font-size:22px; font-weight:700; color:#8a2be2;}
    .bnd {background:#fffbe6; border-left:4px solid #e0c060; padding:8px;
          margin-bottom:8px;}
    .warn {background:#fff0f0; border-left:4px solid #d68a8a; padding:8px;}
  "))),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Product & dose"),
      selectInput("product", "Product", PROD_NAMES, selected = PROD_NAMES[1]),
      sliderInput("dose_mult", "Dose multiplier (1 = 240 U standard)",
                  0.25, 8, 1, step = 0.25),
      sliderInput("upermL", "Dilution (Units per mL)", 12.5, 200, 50, step = 12.5),

      h4("Target list — Units per muscle"),
      lapply(seq_along(MUSCLES), function(i) {
        numericInput(paste0("d", i), MUSCLES[i],
                     value = c(50, 90, 60, 40, 0, 0, 0, 0)[i],
                     min = 0, max = 400, step = 5)
      }),

      h4("Geometry — the ceiling lives here"),
      sliderInput("rho", "rho: share of a targeted muscle actually reached",
                  0.2, 1.0, RHO_CAL, step = 0.01),
      helpText("Calibrated value 0.55, recovered from the dose-response",
               "plateau. Raising it is what EMG / ultrasound guidance does."),
      checkboxInput("bound", "Show THE BOUND (perfect permanent blockade)",
                    FALSE),

      h4("Schedule"),
      sliderInput("interval", "Re-injection interval (days)", 42, 240, 84,
                  step = 7),
      sliderInput("n_inj", "Number of injections", 1, 12, 1, step = 1),
      sliderInput("end", "Simulate to (days)", 120, 1200, 260, step = 20),

      h4("Adjuncts"),
      sliderInput("thp", "Trihexyphenidyl (mg/d)", 0, 30, 0, step = 2),
      sliderInput("bac", "Baclofen (mg/d)", 0, 80, 0, step = 10),
      sliderInput("clz", "Clonazepam (mg/d)", 0, 4, 0, step = 0.5),
      checkboxInput("dbs", "GPi deep brain stimulation from day 0", FALSE),
      sliderInput("denerv", "Selective peripheral denervation (fraction of unreached drive)",
                  0, 0.8, 0, step = 0.05)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",

        tabPanel("1 · Patient & plan",
          h4("Injection plan"), tableOutput("plan"),
          h4("What the plan implies"), htmlOutput("plansummary"),
          h4("Muscle torque shares"), plotOutput("wplot", height = "260px")),

        tabPanel("2 · THE BOUND",
          div(class = "bnd",
            h4("L(t) ≥ L_min = (1 − φ)·D_cen,  φ = ρ · Σ w(injected)"),
            p(strong("There is no drug parameter on the right-hand side."),
              "Not dose, not potency, not serotype, not dilution, not interval.",
              "φ is a property of the injection plan and of this patient's",
              "anatomy. Try to beat the dashed line with the dose slider.")),
          htmlOutput("boundtext"),
          plotOutput("boundplot", height = "380px"),
          h4("The ceiling as a function of φ"),
          plotOutput("phicurve", height = "300px"),
          h4("Four levers, priced at the current plan"),
          tableOutput("levers")),

        tabPanel("3 · Toxin PK",
          h4("A → B → C: the drug is gone long before the effect peaks"),
          plotOutput("pk", height = "340px"),
          htmlOutput("pktext"),
          h4("Toxin arriving where it was not wanted"),
          plotOutput("pkspread", height = "260px")),

        tabPanel("4 · Chemodenervation PD",
          h4("Intact SNARE fraction, transmission, and sprouting"),
          plotOutput("pd", height = "340px"),
          h4("The safety factor makes E(S) a threshold function"),
          plotOutput("safety", height = "280px"),
          htmlOutput("pdtext")),

        tabPanel("5 · Clinical endpoints",
          h4("TWSTRS subscales"), plotOutput("tw", height = "360px"),
          h4("Change from baseline, against the two fitted anchors"),
          plotOutput("dtw", height = "300px"),
          tableOutput("endtable")),

        tabPanel("6 · Dose-response",
          h4("The plateau is not saturation — it is a measurement of φ"),
          plotOutput("dr", height = "360px"),
          htmlOutput("drtext"),
          h4("Benefit is logarithmic in dose; spread is linear"),
          plotOutput("drlog", height = "300px")),

        tabPanel("7 · Adverse effects",
          h4("Spread compartments and their consequences"),
          plotOutput("ae", height = "340px"),
          div(class = "warn", htmlOutput("aetext")),
          h4("Efficacy and dysphagia risk against dilution"),
          plotOutput("dilution", height = "300px")),

        tabPanel("8 · Chronic cycles",
          h4("Peaks are peripheral and cyclical; TROUGHS are central and ratchet"),
          plotOutput("chronic", height = "360px"),
          h4("Central drive D_cen across cycles"),
          plotOutput("dcen", height = "260px"),
          tableOutput("cyctable")),

        tabPanel("9 · Immunogenicity",
          h4("Antigen is protein load, not Units"),
          plotOutput("imm", height = "340px"),
          htmlOutput("immtext"),
          h4("Cumulative antigen by product at this plan"),
          tableOutput("loadtable")),

        tabPanel("10 · Scenario comparison",
          h4("All 15 predefined scenarios"),
          p("Each row is run from the model file's own scenario list, so this",
            "table is reproducible outside the app."),
          tableOutput("scen"))
      )
    )
  )
)

## ---------------------------------------------------------------------
## SERVER
## ---------------------------------------------------------------------
server <- function(input, output, session) {

  pattern <- reactive({
    p <- vapply(1:8, function(i) {
      v <- input[[paste0("d", i)]]
      if (is.null(v) || is.na(v)) 0 else v
    }, numeric(1))
    p * input$dose_mult
  })

  PRODUCTS <- list(
    incobotulinumtoxinA = list(serotype = "A", unit_scale = 1.00, load = 0.44/100,
                               klc_mult = 1.00, auto_pref = 1.0),
    onabotulinumtoxinA  = list(serotype = "A", unit_scale = 1.00, load = 5.00/100,
                               klc_mult = 1.00, auto_pref = 1.0),
    abobotulinumtoxinA  = list(serotype = "A", unit_scale = 0.34, load = 4.35/500,
                               klc_mult = 1.00, auto_pref = 1.0),
    daxibotulinumtoxinA = list(serotype = "A", unit_scale = 1.00, load = 0.50/100,
                               klc_mult = 0.50, auto_pref = 1.0),
    rimabotulinumtoxinB = list(serotype = "B", unit_scale = 0.03, load = 5.00/5000,
                               klc_mult = 1.45, auto_pref = 4.0)
  )

  ## build and run one simulation for an arbitrary pattern
  runsim <- function(pat, rho_mult = NULL, block = FALSE, product = NULL,
                     n_inj = NULL, interval = NULL, end = NULL, upermL = NULL) {
    pr <- PRODUCTS[[if (is.null(product)) input$product else product]]
    n_inj <- if (is.null(n_inj)) input$n_inj else n_inj
    interval <- if (is.null(interval)) input$interval else interval
    end <- if (is.null(end)) input$end else end
    upermL <- if (is.null(upermL)) input$upermL else upermL
    cmtA <- paste0("A", 1:8, "_")

    m <- param(mod,
      c(setNames(as.list(pat), paste0("D", 1:8)),
        setNames(as.list(pat / upermL), paste0("V", 1:8)),
        list(KLC = 0.02465 * pr$klc_mult, AUTOPREF = pr$auto_pref,
             RHO = input$rho,
             RHOMULT = if (is.null(rho_mult)) 1 else rho_mult,
             BLOCK = if (block) 1 else 0,
             THPDOSE = input$thp, BACDOSE = input$bac, CLZDOSE = input$clz,
             DBSON = if (input$dbs) 0 else -1,
             DENERV = input$denerv)))

    if (all(pat <= 0) || n_inj < 1) return(mrgsim(m, end = end, delta = 0.5))

    times <- (seq_len(n_inj) - 1) * interval
    ev <- do.call(rbind, lapply(times, function(t0) {
      rbind(
        do.call(rbind, lapply(which(pat > 0), function(i)
          data.frame(time = t0, cmt = cmtA[i],
                     amt = pat[i] * pr$unit_scale, evid = 1))),
        data.frame(time = t0,
                   cmt = if (pr$serotype == "A") "AGA" else "AGB",
                   amt = sum(pat) * pr$load, evid = 1))
    }))
    mrgsim(m, data = ev, end = end, delta = 0.5)
  }

  sim  <- reactive(as.data.frame(runsim(pattern())))
  simB <- reactive(as.data.frame(runsim(pattern(), block = TRUE)))

  phi_now <- reactive(phi_of(pattern(), input$rho))

  ## ---------------- tab 1 ----------------
  output$plan <- renderTable({
    pat <- pattern()
    data.frame(Muscle = MUSCLES, `Torque share w` = W_DEFAULT,
               `Units` = pat,
               `Volume (mL)` = round(pat / input$upermL, 2),
               `On the list` = ifelse(pat > 0, "yes", "—"),
               check.names = FALSE)
  }, digits = 2)

  output$plansummary <- renderUI({
    pat <- pattern(); ph <- phi_now()
    fl <- tw_floor(ph)
    HTML(sprintf(
      "<p>Total dose <b>%.0f U</b> in <b>%d</b> muscles, total volume
       <b>%.1f mL</b>.</p>
       <p>Σw over the injected muscles = <b>%.3f</b>; ρ = <b>%.2f</b>;
       so <span class='keynum'>φ = %.3f</span>.</p>
       <p><b>%.0f%% of this patient's dystonic drive is in tissue this plan
       does not affect at all.</b> The asymptotic best possible TWSTRS total
       at this φ is <b>%.1f</b> (from an untreated %.1f), i.e. a ceiling of
       <b>%+.1f</b> points — and that is for permanent, complete, side-effect-free
       blockade, which no injection achieves.</p>",
      sum(pat), sum(pat > 0), sum(pat) / input$upermL,
      sum(W_DEFAULT[pat > 0]), input$rho, ph,
      100 * (1 - ph), fl["total"], TW0, fl["total"] - TW0))
  })

  output$wplot <- renderPlot({
    pat <- pattern()
    cols <- ifelse(pat > 0, "#e0a878", "#bbbbbb")
    bp <- barplot(W_DEFAULT, names.arg = NULL, col = cols, border = NA,
                  ylab = "share of dystonic torque, w",
                  main = "orange = on the target list, grey = untouched")
    text(bp, par("usr")[3] - 0.008, srt = 30, adj = 1, labels = MUSCLES,
         xpd = TRUE, cex = 0.75)
    abline(h = 0)
  })

  ## ---------------- tab 2 : THE BOUND ----------------
  output$boundtext <- renderUI({
    d <- sim(); b <- simB(); ph <- phi_now()
    ach <- min(d$TWSTRS_TOTAL) - TW0
    bnd <- min(b$TWSTRS_TOTAL) - TW0
    HTML(sprintf(
      "<p>At φ = <b>%.3f</b>: this plan reaches <b>%+.2f</b> TWSTRS points at
       its nadir. Perfect permanent blockade of the same motor units reaches
       <b>%+.2f</b>. The whole of what a better or bigger toxin could ever add
       is <span class='keynum'>%+.2f</span> points.</p>",
      ph, ach, bnd, bnd - ach))
  })

  output$boundplot <- renderPlot({
    d <- sim(); b <- simB()
    fl <- tw_floor(phi_now())["total"]
    yl <- range(c(d$TWSTRS_TOTAL, b$TWSTRS_TOTAL, fl, TW0)) + c(-2, 2)
    plot(d$time, d$TWSTRS_TOTAL, type = "l", lwd = 2.5, col = "#2b6cb0",
         ylim = yl, xlab = "day", ylab = "TWSTRS total (0–85)",
         main = "you cannot get below the dashed lines with any dose")
    lines(b$time, b$TWSTRS_TOTAL, lwd = 2, col = "#c05621", lty = 1)
    abline(h = TW0, col = "grey60", lty = 3)
    abline(h = fl, col = "#b7791f", lty = 2, lwd = 2)
    legend("topright", bty = "n", cex = 0.85,
           legend = c("this plan", "perfect blockade (the integrated bound)",
                      "asymptotic floor (1−φ), held forever",
                      "untreated baseline"),
           col = c("#2b6cb0", "#c05621", "#b7791f", "grey60"),
           lty = c(1, 1, 2, 3), lwd = c(2.5, 2, 2, 1))
  })

  output$phicurve <- renderPlot({
    phs <- seq(0, 1, by = 0.01)
    fl <- vapply(phs, function(p) tw_floor(p)["total"], numeric(1))
    plot(phs, fl - TW0, type = "l", lwd = 2.5, col = "#b7791f",
         xlab = expression(phi ~ "= " ~ rho %.% Sigma * w),
         ylab = "best attainable ΔTWSTRS (points)",
         main = "the ceiling is a function of geometry alone")
    ph <- phi_now()
    points(ph, tw_floor(ph)["total"] - TW0, pch = 19, cex = 1.6, col = "#8a2be2")
    text(ph, tw_floor(ph)["total"] - TW0, "  this plan", adj = 0, cex = 0.9)
    abline(v = c(0.36, 0.48), col = "grey70", lty = 3)
    text(0.36, min(fl - TW0) * 0.35, "standard\nsurface plan", cex = 0.75,
         col = "grey30")
    text(0.48, min(fl - TW0) * 0.65, "extended\nguided plan", cex = 0.75,
         col = "grey30")
  })

  output$levers <- renderTable({
    pat <- pattern()
    ach <- min(sim()$TWSTRS_TOTAL) - TW0
    b_now <- min(simB()$TWSTRS_TOTAL) - TW0
    b_rho <- min(as.data.frame(
      runsim(pat, rho_mult = 1 / input$rho, block = TRUE))$TWSTRS_TOTAL) - TW0
    pat_ext <- pat
    pat_ext[c(5, 7)] <- pmax(pat_ext[c(5, 7)], 30)
    b_ext <- min(as.data.frame(
      runsim(pat_ext, rho_mult = 1 / input$rho, block = TRUE))$TWSTRS_TOTAL) - TW0
    data.frame(
      Lever = c("1 · more or better toxin (at this φ)",
                "2 · perfect needle placement (ρ → 1)",
                "3 · add the deep muscles (Σw)",
                "4 · this plan, as it stands"),
      `Ceiling ΔTWSTRS` = c(b_now, b_rho, b_ext, ach),
      `Headroom` = c(b_now - ach, b_rho - b_now, b_ext - b_rho, 0),
      check.names = FALSE)
  }, digits = 2)

  ## ---------------- tab 3 : PK ----------------
  output$pk <- renderPlot({
    d <- sim()
    par(mar = c(4, 4, 3, 4))
    plot(d$time, d$C_SCM, type = "l", lwd = 2.5, col = "#c53030",
         xlab = "day", ylab = "active light chain C (U-eq), SCM",
         main = "the drug leaves in hours; its light chain stays for weeks")
    par(new = TRUE)
    plot(d$time, d$S_SCM, type = "l", lwd = 2, col = "#2b6cb0",
         axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1))
    axis(4); mtext("intact SNARE fraction S", side = 4, line = 2.4)
    legend("right", bty = "n", legend = c("C, active light chain", "S, intact SNARE"),
           col = c("#c53030", "#2b6cb0"), lwd = 2, cex = 0.85)
  })

  output$pktext <- renderUI({
    HTML("<p>Free toxin has a half-life of about <b>6 hours</b> in the muscle;
      the internalising pool about <b>1.5 days</b>; the catalytically active
      light chain about <b>28 days</b> (this is the calibrated parameter); and
      SNAP-25 resynthesis about <b>5 days</b>. <b>Duration of action is not a
      pharmacokinetic property of the injectate.</b> By day 2 essentially none
      of the injected toxin remains, and the clinical effect has not yet
      peaked.</p>")
  })

  output$pkspread <- renderPlot({
    d <- sim()
    plot(d$time, d$DEF_SW, type = "l", lwd = 2.5, col = "#c05621", ylim = c(0, 1),
         xlab = "day", ylab = "transmission deficit",
         main = "swallow and autonomic compartments — nobody aimed here")
    lines(d$time, d$DEF_AU, lwd = 2, col = "#2f855a")
    legend("topright", bty = "n", legend = c("pharyngeal (→ dysphagia)",
           "autonomic (→ dry mouth)"), col = c("#c05621", "#2f855a"),
           lwd = 2, cex = 0.85)
  })

  ## ---------------- tab 4 : PD ----------------
  output$pd <- renderPlot({
    d <- sim()
    plot(d$time, d$S_SCM, type = "l", lwd = 2.5, col = "#2b6cb0", ylim = c(0, 1),
         xlab = "day", ylab = "fraction",
         main = "SNARE, parent-terminal transmission, and sprouts")
    lines(d$time, d$E_SCM, lwd = 2.5, col = "#c53030")
    lines(d$time, d$Q_SCM, lwd = 2, col = "#2f855a", lty = 2)
    legend("right", bty = "n", cex = 0.85,
           legend = c("S, intact SNARE", "E(S), parent transmission",
                      "Q, sprout capacity"),
           col = c("#2b6cb0", "#c53030", "#2f855a"), lty = c(1, 1, 2), lwd = 2)
  })

  output$safety <- renderPlot({
    S <- seq(0, 1, by = 0.005)
    r <- S^3 / (S^3 + 0.20^3)
    E <- (3 * r / (1 + 3 * r)) / (3 / 4)
    plot(S, E, type = "l", lwd = 3, col = "#c53030",
         xlab = "intact SNARE fraction S", ylab = "transmission efficacy E(S)",
         main = "half the substrate can go with no measurable weakness")
    lines(S, r, lwd = 2, col = "grey55", lty = 2)
    abline(v = 0.20, col = "grey70", lty = 3)
    legend("bottomright", bty = "n", cex = 0.85,
           legend = c("E(S) with safety factor SF = 3", "raw release capacity r(S)"),
           col = c("#c53030", "grey55"), lty = c(1, 2), lwd = c(3, 2))
  })

  output$pdtext <- renderUI({
    HTML("<p>Two nonlinearities, both real. SNARE cooperativity makes release
      capacity a Hill function of S; the neuromuscular <b>safety factor</b>
      means release must fall below about 1/3 of normal before any weakness
      appears. Together they make E(S) a <b>threshold</b> function — which is
      why the dose→duration relationship is logarithmic, and why the clinical
      dose-response plateaus.</p>
      <p>Sprouts draw on the same cytosolic SNARE pool the light chain has
      destroyed, so their release is gated by S too. In this model that makes
      sprouting a <b>minor</b> contributor to wearing-off — worth about three
      days of duration — which contradicts the usual account. See A5.</p>")
  })

  ## ---------------- tab 5 : endpoints ----------------
  output$tw <- renderPlot({
    d <- sim()
    plot(d$time, d$TWSTRS_TOTAL, type = "l", lwd = 3, col = "#1a202c",
         ylim = c(0, 50), xlab = "day", ylab = "points",
         main = "TWSTRS total and its three subscales")
    lines(d$time, d$SEV, lwd = 2, col = "#2b6cb0")
    lines(d$time, d$DISAB, lwd = 2, col = "#2f855a")
    lines(d$time, d$TWSTRS_PAIN, lwd = 2, col = "#c53030")
    abline(h = TW0, col = "grey65", lty = 3)
    legend("topright", bty = "n", cex = 0.85,
           legend = c("total (0–85)", "severity (0–35)", "disability (0–30)",
                      "pain (0–20)"),
           col = c("#1a202c", "#2b6cb0", "#2f855a", "#c53030"), lwd = c(3,2,2,2))
  })

  output$dtw <- renderPlot({
    d <- sim()
    plot(d$time / 7, d$TWSTRS_TOTAL - TW0, type = "l", lwd = 3, col = "#1a202c",
         xlab = "week", ylab = "ΔTWSTRS total from baseline",
         main = "filled circles = the two anchors the model was fitted to")
    abline(h = 0, col = "grey65", lty = 3)
    abline(h = -MCID, col = "#b7791f", lty = 2)
    text(max(d$time / 7) * 0.75, -MCID + 0.7,
         sprintf("MCID = %.1f points", MCID), cex = 0.8, col = "#b7791f")
    points(ANCHORS$week, ANCHORS$dTW, pch = 19, cex = 1.8, col = "#8a2be2")
  })

  output$endtable <- renderTable({
    d <- sim()
    at <- function(x) approx(d$time, d$TWSTRS_TOTAL, xout = x, rule = 2)$y - TW0
    g <- TW0 - d$TWSTRS_TOTAL
    i <- which.max(g)
    j <- which(seq_along(g) > i & g <= MCID)
    dur <- if (length(j)) d$time[min(j)] else NA_real_
    data.frame(Quantity = c("ΔTWSTRS week 4 (fitted)",
                            "ΔTWSTRS week 12 (fitted)",
                            "nadir ΔTWSTRS", "day of nadir",
                            sprintf("duration of benefit (gain > %.1f)", MCID),
                            "peak P(dysphagia)", "peak P(neck weakness)"),
               Value = c(at(28), at(84), min(d$TWSTRS_TOTAL) - TW0,
                         d$time[which.min(d$TWSTRS_TOTAL)], dur,
                         max(d$P_DYSPHAGIA), max(d$P_NECKWEAK)))
  }, digits = 2)

  ## ---------------- tab 6 : dose-response ----------------
  dr_data <- reactive({
    pat <- pattern() / input$dose_mult
    mults <- c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16)
    do.call(rbind, lapply(mults, function(f) {
      d <- as.data.frame(runsim(pat * f, n_inj = 1, end = 420))
      g <- TW0 - d$TWSTRS_TOTAL
      i <- which.max(g)
      j <- which(seq_along(g) > i & g <= MCID)
      data.frame(dose = sum(pat * f), mult = f,
                 nadir = min(d$TWSTRS_TOTAL) - TW0,
                 duration = if (length(j)) d$time[min(j)] else max(d$time),
                 p_dys = max(d$P_DYSPHAGIA),
                 peak_csw = max(d$CSW))
    }))
  })

  output$dr <- renderPlot({
    r <- dr_data(); ph <- phi_now()
    par(mar = c(4, 4, 3, 4))
    plot(r$dose, r$nadir, type = "b", pch = 19, lwd = 2.5, col = "#2b6cb0",
         log = "x", xlab = "total dose (U, log scale)",
         ylab = "nadir ΔTWSTRS total",
         main = "benefit plateaus; dysphagia does not")
    abline(h = tw_floor(ph)["total"] - TW0, col = "#b7791f", lty = 2, lwd = 2)
    par(new = TRUE)
    plot(r$dose, r$p_dys, type = "b", pch = 17, lwd = 2, col = "#c53030",
         log = "x", axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1))
    axis(4); mtext("P(dysphagia)", side = 4, line = 2.4)
    legend("bottomleft", bty = "n", cex = 0.85,
           legend = c("nadir ΔTWSTRS", "asymptotic floor at this φ",
                      "P(dysphagia)"),
           col = c("#2b6cb0", "#b7791f", "#c53030"), lty = c(1, 2, 1),
           pch = c(19, NA, 17), lwd = 2)
  })

  output$drtext <- renderUI({
    HTML("<p>A dose-response curve that flattens while adverse effects keep
      climbing is <b>not pharmacological saturation</b> — the toxin has not run
      out of anything. It is the needle running out of disease to reach. Read
      that way, the plateau is a <b>measurement of φ</b>, and fitting it gives
      ρ ≈ 0.55: a standard surface injection affects only about 36% of the
      dystonic drive.</p>")
  })

  output$drlog <- renderPlot({
    r <- dr_data()
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    plot(r$mult, r$duration, type = "b", pch = 19, log = "x", lwd = 2,
         col = "#2f855a", xlab = "dose multiple (log)", ylab = "duration (d)",
         main = "duration ~ log(dose)")
    plot(r$mult, r$peak_csw, type = "b", pch = 17, lwd = 2, col = "#c53030",
         xlab = "dose multiple (linear)", ylab = "peak pharyngeal light chain",
         main = "spread ~ dose (or worse)")
    par(mfrow = c(1, 1))
  })

  ## ---------------- tab 7 : adverse effects ----------------
  output$ae <- renderPlot({
    d <- sim()
    plot(d$time, d$P_DYSPHAGIA, type = "l", lwd = 2.5, col = "#c53030",
         ylim = c(0, 1), xlab = "day", ylab = "probability",
         main = "modelled adverse-effect probabilities")
    lines(d$time, d$P_NECKWEAK, lwd = 2.5, col = "#b7791f")
    lines(d$time, d$P_DRYMOUTH, lwd = 2, col = "#2f855a", lty = 2)
    legend("topright", bty = "n", cex = 0.85,
           legend = c("dysphagia", "neck weakness / head drop", "dry mouth"),
           col = c("#c53030", "#b7791f", "#2f855a"), lty = c(1, 1, 2), lwd = 2)
  })

  output$aetext <- renderUI({
    HTML("<p><b>The muscles that generate the dystonic torque are the muscles
      that hold the head up.</b> The antagonist of the disease is the agonist
      of the posture. That is a geometric conflict, not a dosing one, and it is
      why the bound on tab 2 is not clinically attainable even in principle:
      perfect blockade of the posterior group carries a high probability of head
      drop. No improvement in any toxin resolves it.</p>")
  })

  output$dilution <- renderPlot({
    pat <- pattern()
    ups <- c(200, 100, 50, 25, 12.5)
    r <- do.call(rbind, lapply(ups, function(u) {
      d <- as.data.frame(runsim(pat, n_inj = 1, end = 260, upermL = u))
      data.frame(upermL = u, nadir = min(d$TWSTRS_TOTAL) - TW0,
                 p_dys = max(d$P_DYSPHAGIA))
    }))
    par(mar = c(4, 4, 3, 4))
    plot(r$upermL, r$nadir, type = "b", pch = 19, lwd = 2.5, col = "#2b6cb0",
         log = "x", xlab = "Units per mL (log; more dilute to the left)",
         ylab = "nadir ΔTWSTRS", main = "concentration is a nearly free safety lever")
    par(new = TRUE)
    plot(r$upermL, r$p_dys, type = "b", pch = 17, lwd = 2, col = "#c53030",
         log = "x", axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1))
    axis(4); mtext("P(dysphagia)", side = 4, line = 2.4)
  })

  ## ---------------- tab 8 : chronic ----------------
  output$chronic <- renderPlot({
    n <- max(input$n_inj, 6)
    d <- as.data.frame(runsim(pattern(), n_inj = n, end = n * input$interval))
    plot(d$time, d$TWSTRS_TOTAL, type = "l", lwd = 2, col = "#1a202c",
         xlab = "day", ylab = "TWSTRS total",
         main = "peaks stay put; troughs drift downward")
    abline(h = TW0, col = "grey65", lty = 3)
    tr <- vapply(seq_len(n), function(k) {
      m <- d$time >= (k - 1) * input$interval & d$time <= k * input$interval
      max(d$TWSTRS_TOTAL[m])
    }, numeric(1))
    points(seq_len(n) * input$interval, tr, pch = 19, col = "#c53030", cex = 1.2)
    lines(seq_len(n) * input$interval, tr, col = "#c53030", lty = 2)
    legend("bottomleft", bty = "n", cex = 0.85, legend = "cycle troughs",
           col = "#c53030", pch = 19, lty = 2)
  })

  output$dcen <- renderPlot({
    n <- max(input$n_inj, 6)
    d <- as.data.frame(runsim(pattern(), n_inj = n, end = n * input$interval))
    plot(d$time, d$DCEN, type = "l", lwd = 2.5, col = "#8a2be2",
         xlab = "day", ylab = expression(D[cen]),
         main = "BoNT blocks intrafusal terminals too, so it reaches the centre")
  })

  output$cyctable <- renderTable({
    n <- max(input$n_inj, 6)
    d <- as.data.frame(runsim(pattern(), n_inj = n, end = n * input$interval))
    do.call(rbind, lapply(seq_len(n), function(k) {
      m <- d$time >= (k - 1) * input$interval & d$time <= k * input$interval
      data.frame(cycle = k, nadir = min(d$TWSTRS_TOTAL[m]),
                 trough = max(d$TWSTRS_TOTAL[m]),
                 D_cen = approx(d$time, d$DCEN,
                                xout = k * input$interval, rule = 2)$y,
                 mean_gain = mean(d$GAIN[m]))
    }))
  }, digits = 3)

  ## ---------------- tab 9 : immunogenicity ----------------
  output$imm <- renderPlot({
    n <- max(input$n_inj, 20)
    d <- as.data.frame(runsim(pattern(), n_inj = n, end = n * input$interval))
    par(mar = c(4, 4, 3, 4))
    plot(d$time, d$NABA, type = "l", lwd = 2.5, col = "#c53030",
         xlab = "day", ylab = "neutralising antibody (serotype A)",
         main = "antibody, and the potency it removes")
    lines(d$time, d$NABB, lwd = 2, col = "#2f855a", lty = 2)
    par(new = TRUE)
    plot(d$time, d$GATE_A, type = "l", lwd = 2, col = "#2b6cb0",
         axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1))
    axis(4); mtext("fraction of injected potency surviving", side = 4, line = 2.4)
    legend("left", bty = "n", cex = 0.85,
           legend = c("Nab serotype A", "Nab serotype B", "potency gate A"),
           col = c("#c53030", "#2f855a", "#2b6cb0"), lty = c(1, 2, 1), lwd = 2)
  })

  output$immtext <- renderUI({
    HTML("<p>The antigen is <b>protein load</b>, not Units: 5 ng per 100 U for
      onabotulinumtoxinA against 0.44 ng per 100 U for incobotulinumtoxinA, an
      eleven-fold difference at the same nominal dose. Serotypes A and B share
      no neutralising epitopes, so a switch to B rescues potency — finitely,
      because the B pool immunises independently.</p>
      <p>At the <b>observed</b> antibody rate, however, shorter intervals are
      simply better on this axis, for every product: the median patient never
      loses meaningful potency. The 12-week convention is optimal only for
      patients above roughly the 99th centile of antibody propensity. It is a
      rule calibrated on the tail and applied to everyone — defensible as
      insurance against an absorbing loss, but it is a decision about variance,
      not about the average patient. See A7.</p>")
  })

  output$loadtable <- renderTable({
    pat <- pattern()
    n <- max(input$n_inj, 20)
    do.call(rbind, lapply(names(PRODUCTS), function(nm) {
      pr <- PRODUCTS[[nm]]
      data.frame(Product = nm, Serotype = pr$serotype,
                 `ng per 100 label-U` = pr$load * 100,
                 `cumulative ng over the plan` = sum(pat) * pr$load * n,
                 check.names = FALSE)
    }))
  }, digits = 3)

  ## ---------------- tab 10 : scenarios ----------------
  output$scen <- renderTable({
    nms <- names(MENV$SCENARIOS)
    if (is.null(nms)) nms <- character(0)
    do.call(rbind, lapply(nms, function(nm) {
      out <- try(MENV$sim_scenario(mod, nm, delta = 1), silent = TRUE)
      if (inherits(out, "try-error")) return(NULL)
      d <- as.data.frame(out)
      data.frame(scenario = nm,
                 label = MENV$SCENARIOS[[nm]]$label,
                 phi = phi_of(MENV$SCENARIOS[[nm]]$pattern, RHO_CAL),
                 wk4 = approx(d$time, d$TWSTRS_TOTAL, 28, rule = 2)$y - TW0,
                 nadir = min(d$TWSTRS_TOTAL) - TW0,
                 peak_P_dysphagia = max(d$P_DYSPHAGIA))
    }))
  }, digits = 3)
}

shinyApp(ui, server)
