## ============================================================
##  Obesity QSP — Interactive Shiny Dashboard
##  6 Tabs: Patient Profile · Drug PK · PD Biomarkers ·
##          Clinical Endpoints · Scenario Comparison ·
##          Cardiometabolic Risk
##
##  Drugs modeled:
##   ① Placebo
##   ② Semaglutide 2.4 mg SC QW
##   ③ Tirzepatide 15 mg SC QW
##   ④ Orlistat 120 mg PO TID
##   ⑤ Phentermine/Topiramate 15/92 mg QD
##
##  Author: Claude Code Routine (CCR) — 2026-06-18
## ============================================================

library(shiny)
library(bslib)
library(plotly)
library(DT)
library(dplyr)
library(tidyr)
library(deSolve)

# ── ODE System (pure deSolve; self-contained, no mrgsolve dependency) ──────────
obesity_ode <- function(t, state, pars) {
  with(as.list(c(state, pars)), {

    # Drug concentrations / target engagement
    sema_glp1 <- pmax(0, SEMA_C / (SEMA_C + EC50_sema_glp1))
    tirz_glp1 <- pmax(0, TIRZ_C / (TIRZ_C + EC50_tirz_glp1))
    tirz_gip  <- pmax(0, TIRZ_C / (TIRZ_C + EC50_tirz_gip))

    glp1r_target <- pmin(sema_glp1 + tirz_glp1, 1.0)
    gipr_target  <- pmin(tirz_gip, 1.0)

    orl_eff  <- Imax_orl  * ORL_GUT / (ORL_GUT + IC50_orl)
    cns_eff  <- Emax_cns  * CNS_C   / (CNS_C   + EC50_cns)

    # Food intake reduction
    fi_glp1 <- Emax_food_glp1 * GLP1R_OCC / (GLP1R_OCC + EC50_food_rec)
    fi_gip  <- Emax_food_gip  * GIPR_OCC  / (GIPR_OCC  + EC50_food_rec)
    total_fi <- pmin(fi_glp1 + fi_gip + orl_eff + cns_eff, 0.70)

    kin_food <- kout_food * FOOD_R0

    # Gastric emptying
    gastric_ss <- GASTRIC0 * (1.0 - Emax_gastric * GLP1R_OCC)

    # Ghrelin
    ghrelin_ss <- GHRELIN_R0 * (1.0 - Imax_ghrelin * GLP1R_OCC)

    # Insulin secretion (incretin potentiation)
    inc_eff <- 1.0 +
      Emax_glp1_ins * GLP1R_OCC / (GLP1R_OCC + EC50_ins_sec) +
      Emax_gip_ins  * GIPR_OCC  / (GIPR_OCC  + EC50_ins_sec)

    # Glucose
    EGP_frac <- pmax(1.0 - Imax_ins * INSULIN_P / (INSULIN_P + IC50_ins_HGP), 0.05)

    # Body composition
    kcal_per_h   <- 104.2
    energy_def_h <- kcal_per_h * (1.0 - FOOD_R)

    # Leptin steady-state
    leptin_ss <- k_lep_adip * ADIP

    # Triglyceride steady-state
    trig_ss <- TG0 * (ADIP / A0)^TG_adip_exp

    # HbA1c
    hba1c_ss <- (GLUCOSE_P + 46.7) / 28.7

    # Inflammation
    inflam_ss <- INF0 * (ADIP / A0)^inflam_adip_exp

    # HOMA-IR
    homa_ss <- (INSULIN_P * GLUCOSE_P) / 22.5

    list(c(
      # PK
      dSEMA_GUT  = -ka_sema * SEMA_GUT,
      dSEMA_C    = ka_sema * SEMA_GUT * F_sema / Vd_sema - (CL_sema/Vd_sema)*SEMA_C,
      dTIRZ_GUT  = -ka_tirz * TIRZ_GUT,
      dTIRZ_C    = ka_tirz * TIRZ_GUT * F_tirz / Vd_tirz - (CL_tirz/Vd_tirz)*TIRZ_C,
      dORL_GUT   = -ke_orl * ORL_GUT,
      dCNS_C     = -ke_cns * CNS_C,
      # Receptor occupancy
      dGLP1R_OCC = kout_rec * (glp1r_target - GLP1R_OCC),
      dGIPR_OCC  = kout_rec * (gipr_target  - GIPR_OCC),
      # PD mediators
      dFOOD_R    = kin_food * (1.0 - total_fi) - kout_food * FOOD_R,
      dGASTRIC   = kout_gastric * (gastric_ss - GASTRIC),
      dGHRELIN_R = kout_ghrelin * (ghrelin_ss - GHRELIN_R),
      # Metabolic
      dINSULIN_P = (kout_ins * I0) * (GLUCOSE_P/G0) * inc_eff - kout_ins * INSULIN_P,
      dGLUCOSE_P = kHGP * EGP_frac - kGdisp*(INSULIN_P/I0)*GLUCOSE_P
                   + kout_gluc*(G0*FOOD_R - GLUCOSE_P)*0.1,
      # Composition
      dADIP      = -energy_def_h * adip_frac / calden_fat,
      dBWT_C     = -energy_def_h / calden_fat,
      dLEPTIN_P  = kout_leptin * (leptin_ss - LEPTIN_P),
      # Biomarkers
      dTRIG_P    = kout_trig   * (trig_ss    - TRIG_P),
      dHBA1C_C   = kout_hba1c  * (hba1c_ss   - HBA1C_C),
      dINFLAM_I  = kout_inflam * (inflam_ss  - INFLAM_I),
      dHOMA_IR_C = kout_homa   * (homa_ss    - HOMA_IR_C)
    ))
  })
}

# ── Default parameters ─────────────────────────────────────────────────────────
default_pars <- list(
  ka_sema=0.0177, CL_sema=0.0403, Vd_sema=12.4, F_sema=0.89,
  ka_tirz=0.0187, CL_tirz=0.0500, Vd_tirz=9.20, F_tirz=0.80,
  ka_orl=0.50, ke_orl=0.35, Imax_orl=0.30, IC50_orl=0.55,
  ka_cns=0.15, ke_cns=0.07, Emax_cns=0.40, EC50_cns=0.50,
  EC50_sema_glp1=0.016, EC50_tirz_glp1=0.050, EC50_tirz_gip=0.013,
  kout_rec=0.10,
  Emax_food_glp1=0.42, Emax_food_gip=0.08, EC50_food_rec=0.50,
  kout_food=0.0069, FOOD_R0=1.0,
  kout_gastric=0.10, Emax_gastric=0.25, GASTRIC0=1.0,
  kout_ghrelin=0.06, Imax_ghrelin=0.30, GHRELIN_R0=1.0,
  kout_ins=0.12, Emax_glp1_ins=1.8, Emax_gip_ins=1.2, EC50_ins_sec=0.6,
  kHGP=1.10, kGdisp=0.015, Imax_ins=0.80, IC50_ins_HGP=5.0, kout_gluc=0.05,
  kout_adip=0.00024, calden_fat=7700, adip_frac=0.87, kout_bwt=0.00024,
  kout_leptin=0.30, k_lep_adip=0.714,
  kout_trig=0.008, TG_adip_exp=1.50,
  kout_hba1c=0.00046,
  kout_inflam=0.004, inflam_adip_exp=1.20,
  kout_homa=0.05
)

state_names <- c("SEMA_GUT","SEMA_C","TIRZ_GUT","TIRZ_C","ORL_GUT","CNS_C",
                 "GLP1R_OCC","GIPR_OCC","FOOD_R","GASTRIC","GHRELIN_R",
                 "INSULIN_P","GLUCOSE_P","ADIP","BWT_C","LEPTIN_P",
                 "TRIG_P","HBA1C_C","INFLAM_I","HOMA_IR_C")

# ── Dosing event function ──────────────────────────────────────────────────────
make_events <- function(arm, sim_hours) {
  events <- list()
  if (arm == "sema") {
    # Semaglutide escalation: 0.25→0.5→1.0→1.7→2.4 mg QW
    doses_mg <- c(rep(0.25,2), rep(0.5,2), rep(1.0,4), rep(1.7,4), rep(2.4,60))
    doses_mg <- doses_mg[seq_len(min(length(doses_mg), floor(sim_hours/(7*24))))]
    times    <- seq(0, (length(doses_mg)-1)*7*24, by=7*24)
    nmol     <- doses_mg / 4113.6 * 1e6
    events   <- data.frame(var="SEMA_GUT", time=times, value=nmol, method="add")
  } else if (arm == "tirz") {
    doses_mg <- c(rep(2.5,4), rep(5,4), rep(7.5,4), rep(10,4), rep(12.5,4), rep(15,60))
    doses_mg <- doses_mg[seq_len(min(length(doses_mg), floor(sim_hours/(7*24))))]
    times    <- seq(0, (length(doses_mg)-1)*7*24, by=7*24)
    nmol     <- doses_mg / 4813.5 * 1e6
    events   <- data.frame(var="TIRZ_GUT", time=times, value=nmol, method="add")
  } else if (arm == "orlistat") {
    times    <- seq(0, sim_hours-8, by=8)
    events   <- data.frame(var="ORL_GUT", time=times, value=120, method="add")
  } else if (arm == "cns") {
    times    <- seq(0, sim_hours-24, by=24)
    events   <- data.frame(var="CNS_C", time=times, value=1.0, method="add")
  }
  events
}

# ── Run simulation ─────────────────────────────────────────────────────────────
run_sim <- function(arm, sim_weeks=72, bw0=106, height=1.73, a0=44,
                    i0=18, g0=112, l0=32, tg0=185, hba0=5.8,
                    inf0=1.5, homa0=4.5) {
  sim_hours <- sim_weeks * 7 * 24
  times     <- seq(0, sim_hours, by=24)

  y0 <- c(SEMA_GUT=0, SEMA_C=0, TIRZ_GUT=0, TIRZ_C=0, ORL_GUT=0, CNS_C=0,
           GLP1R_OCC=0, GIPR_OCC=0, FOOD_R=1.0, GASTRIC=1.0, GHRELIN_R=1.0,
           INSULIN_P=i0, GLUCOSE_P=g0,
           ADIP=a0, BWT_C=bw0, LEPTIN_P=l0,
           TRIG_P=tg0, HBA1C_C=hba0, INFLAM_I=inf0, HOMA_IR_C=homa0)

  pars <- c(default_pars, I0=i0, G0=g0, A0=a0, BW0=bw0,
            TG0=tg0, INF0=inf0, L0=l0, HOMA0=homa0)

  ev <- make_events(arm, sim_hours)

  tryCatch({
    out <- ode(y=y0, times=times, func=obesity_ode, parms=pars,
               method="lsoda", events=list(data=ev))
    df <- as.data.frame(out)
    colnames(df)[1] <- "time"
    df$week      <- df$time / (7*24)
    df$bwt_pct   <- (df$BWT_C - bw0) / bw0 * 100
    df$bmi       <- df$BWT_C / height^2
    df$arm       <- arm
    df
  }, error=function(e) { message("ODE error: ", e$message); NULL })
}

# ── Arm metadata ───────────────────────────────────────────────────────────────
arm_meta <- list(
  placebo   = list(label="① Placebo",                    color="#9E9E9E"),
  sema      = list(label="② Semaglutide 2.4 mg QW",     color="#7E57C2"),
  tirz      = list(label="③ Tirzepatide 15 mg QW",      color="#1565C0"),
  orlistat  = list(label="④ Orlistat 120 mg TID",       color="#2E7D32"),
  cns       = list(label="⑤ Phentermine/Topiramate QD", color="#E65100")
)

# ── Shiny UI ───────────────────────────────────────────────────────────────────
ui <- page_navbar(
  title = "Obesity QSP Dashboard",
  theme = bs_theme(bootswatch="flatly", primary="#1565C0"),
  bg    = "#0D47A1", fg="#FFFFFF",

  # ── Tab 1: Patient Profile ──────────────────────────────────────────────────
  nav_panel("① Patient Profile",
    layout_sidebar(
      sidebar=sidebar(
        width=300,
        h5("Patient Characteristics"),
        numericInput("bw0",    "Body Weight (kg)",         value=106, min=60, max=200),
        numericInput("height", "Height (m)",                value=1.73, min=1.4, max=2.0, step=0.01),
        numericInput("age",    "Age (years)",               value=45, min=18, max=80),
        selectInput("sex",     "Sex",                       choices=c("Female","Male")),
        hr(),
        h5("Metabolic Baseline"),
        numericInput("g0",     "Fasting Glucose (mg/dL)",  value=112, min=70, max=300),
        numericInput("i0",     "Fasting Insulin (µU/mL)",  value=18, min=2, max=100),
        numericInput("tg0",    "Triglycerides (mg/dL)",    value=185, min=50, max=600),
        numericInput("hba0",   "HbA1c (%)",                value=5.8, min=4, max=12),
        numericInput("a0",     "Adipose Mass (kg)",         value=44, min=5, max=120),
        hr(),
        h5("Simulation"),
        sliderInput("sim_weeks","Duration (weeks)",         min=12, max=104, value=72, step=4),
        actionButton("run_sim","▶ Run Simulation",
                     class="btn-primary btn-lg w-100")
      ),
      # Main area: patient summary cards
      layout_columns(
        col_widths=c(3,3,3,3),
        value_box("BMI", textOutput("bmi_val"), showcase=bsicons::bs_icon("person-fill"),
                  theme="primary"),
        value_box("Obesity Class", textOutput("obes_class"), showcase=bsicons::bs_icon("clipboard2-pulse"),
                  theme="danger"),
        value_box("HOMA-IR", textOutput("homa_val"), showcase=bsicons::bs_icon("activity"),
                  theme="warning"),
        value_box("Metabolic Risk", textOutput("risk_val"), showcase=bsicons::bs_icon("heart-pulse"),
                  theme="info")
      ),
      br(),
      plotlyOutput("natural_hist", height="380px"),
      br(),
      DT::dataTableOutput("patient_tbl")
    )
  ),

  # ── Tab 2: Drug PK ──────────────────────────────────────────────────────────
  nav_panel("② Drug PK",
    layout_sidebar(
      sidebar=sidebar(
        width=280,
        h5("Drug Selection"),
        selectInput("pk_drug","Drug",
                    choices=c("Semaglutide 2.4mg QW"="sema",
                              "Tirzepatide 15mg QW"="tirz",
                              "Orlistat 120mg TID"="orlistat",
                              "Phentermine/Topiramate"="cns")),
        checkboxInput("pk_log","Log Y-axis",value=FALSE),
        hr(),
        h5("PK Reference"),
        tableOutput("pk_ref_tbl")
      ),
      plotlyOutput("pk_conc_plot", height="350px"),
      br(),
      plotlyOutput("pk_rec_occ_plot", height="300px")
    )
  ),

  # ── Tab 3: PD Biomarkers ────────────────────────────────────────────────────
  nav_panel("③ PD Biomarkers",
    layout_sidebar(
      sidebar=sidebar(
        width=260,
        h5("Display Options"),
        checkboxGroupInput("arms_pd","Arms to show",
                           choices=c("Placebo"="placebo","Semaglutide"="sema",
                                     "Tirzepatide"="tirz","Orlistat"="orlistat",
                                     "CNS Agent"="cns"),
                           selected=c("placebo","sema","tirz")),
        sliderInput("pd_week_range","Week range",min=0,max=104,value=c(0,72),step=4)
      ),
      layout_columns(
        col_widths=c(6,6),
        plotlyOutput("pd_insulin",  height="280px"),
        plotlyOutput("pd_glucose",  height="280px"),
        plotlyOutput("pd_food",     height="280px"),
        plotlyOutput("pd_ghrelin",  height="280px"),
        plotlyOutput("pd_leptin",   height="280px"),
        plotlyOutput("pd_inflam",   height="280px")
      )
    )
  ),

  # ── Tab 4: Clinical Endpoints ────────────────────────────────────────────────
  nav_panel("④ Clinical Endpoints",
    layout_sidebar(
      sidebar=sidebar(
        width=260,
        h5("View Options"),
        checkboxGroupInput("arms_ep","Arms",
                           choices=c("Placebo"="placebo","Semaglutide"="sema",
                                     "Tirzepatide"="tirz","Orlistat"="orlistat",
                                     "CNS Agent"="cns"),
                           selected=c("placebo","sema","tirz","orlistat","cns")),
        sliderInput("ep_week","Summary week",min=4,max=104,value=52,step=4)
      ),
      layout_columns(
        col_widths=c(6,6),
        plotlyOutput("ep_bwt",    height="300px"),
        plotlyOutput("ep_bmi",    height="300px"),
        plotlyOutput("ep_hba1c",  height="300px"),
        plotlyOutput("ep_trig",   height="300px"),
        plotlyOutput("ep_homa",   height="300px"),
        plotlyOutput("ep_waist",  height="300px")
      )
    )
  ),

  # ── Tab 5: Scenario Comparison ───────────────────────────────────────────────
  nav_panel("⑤ Scenario Comparison",
    layout_sidebar(
      sidebar=sidebar(
        width=260,
        h5("Comparison Options"),
        sliderInput("cmp_week","Summary at Week",min=4,max=104,value=52,step=4),
        selectInput("cmp_var","Primary Variable",
                    choices=c("Body Weight %Δ"="bwt_pct",
                              "Fasting Glucose"="GLUCOSE_P",
                              "HbA1c"="HBA1C_C",
                              "Triglycerides"="TRIG_P",
                              "HOMA-IR"="HOMA_IR_C",
                              "Leptin"="LEPTIN_P",
                              "Inflammation"="INFLAM_I")),
        downloadButton("dl_cmp","Download CSV")
      ),
      plotlyOutput("cmp_bar",    height="350px"),
      br(),
      plotlyOutput("cmp_spider", height="350px"),
      br(),
      DT::dataTableOutput("cmp_table")
    )
  ),

  # ── Tab 6: Cardiometabolic Risk ──────────────────────────────────────────────
  nav_panel("⑥ Cardiometabolic Risk",
    layout_sidebar(
      sidebar=sidebar(
        width=280,
        h5("Risk Calculator Inputs"),
        selectInput("risk_drug","Drug Scenario",
                    choices=c("Placebo"="placebo","Semaglutide"="sema",
                              "Tirzepatide"="tirz","Orlistat"="orlistat",
                              "CNS Agent"="cns")),
        sliderInput("risk_week","Evaluate at Week",min=4,max=104,value=52,step=4),
        numericInput("ldl0","LDL-C (mg/dL)",value=130,min=50,max=300),
        numericInput("hdl0","HDL-C (mg/dL)",value=45,min=20,max=100),
        checkboxInput("smoker","Current Smoker",value=FALSE),
        hr(),
        h5("Dose-Response"),
        selectInput("dr_drug","Drug",
                    choices=c("Semaglutide"="sema","Tirzepatide"="tirz")),
        sliderInput("dr_weeks","Duration (wk)",min=12,max=104,value=68,step=4)
      ),
      layout_columns(
        col_widths=c(3,3,3,3),
        value_box("Wt Loss",  textOutput("risk_wt"),   showcase=bsicons::bs_icon("arrow-down"),    theme="success"),
        value_box("ASCVD Δ",  textOutput("risk_ascvd"),showcase=bsicons::bs_icon("heart-pulse"),   theme="danger"),
        value_box("HbA1c Δ",  textOutput("risk_hba"),  showcase=bsicons::bs_icon("droplet-fill"),  theme="warning"),
        value_box("TG Δ",     textOutput("risk_tg"),   showcase=bsicons::bs_icon("activity"),      theme="info")
      ),
      br(),
      layout_columns(
        col_widths=c(6,6),
        plotlyOutput("risk_scatter", height="320px"),
        plotlyOutput("dr_curve",     height="320px")
      ),
      br(),
      DT::dataTableOutput("risk_tbl")
    )
  )
)

# ── Shiny Server ───────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive: run all 5 simulations whenever inputs change
  sims <- eventReactive(input$run_sim, {
    withProgress(message="Running QSP simulations...", {
      params <- list(bw0=input$bw0, height=input$height,
                     g0=input$g0, i0=input$i0, tg0=input$tg0,
                     hba0=input$hba0, a0=input$a0,
                     sim_weeks=input$sim_weeks)

      results <- list()
      arms    <- c("placebo","sema","tirz","orlistat","cns")
      for(i in seq_along(arms)) {
        setProgress(i/5, detail=arm_meta[[arms[i]]]$label)
        df <- do.call(run_sim, c(list(arm=arms[i]), params))
        if(!is.null(df)) results[[arms[i]]] <- df
      }
      results
    })
  }, ignoreNULL=FALSE)

  # Combine sims into one data frame
  all_df <- reactive({
    s <- sims()
    if(length(s)==0) return(NULL)
    bind_rows(lapply(names(s), function(a) {
      s[[a]] %>% mutate(arm_label = arm_meta[[a]]$label,
                        arm_color = arm_meta[[a]]$color)
    }))
  })

  # ── Tab 1: Patient Summary ────────────────────────────────────────────────────
  output$bmi_val    <- renderText(sprintf("%.1f kg/m²", input$bw0/input$height^2))
  output$obes_class <- renderText({
    bmi <- input$bw0/input$height^2
    if(bmi<25)"Normal" else if(bmi<30)"Overweight" else if(bmi<35)"Class I" else if(bmi<40)"Class II" else "Class III"
  })
  output$homa_val   <- renderText(sprintf("%.1f", input$i0*input$g0/22.5))
  output$risk_val   <- renderText({
    homa <- input$i0*input$g0/22.5
    bmi  <- input$bw0/input$height^2
    if(homa>5|bmi>40)"High" else if(homa>2.5|bmi>30)"Moderate" else "Low"
  })

  output$natural_hist <- renderPlotly({
    df <- all_df()
    if(is.null(df)) return(NULL)
    pl <- df %>%
      filter(arm=="placebo") %>%
      plot_ly(x=~week, y=~bwt_pct, type="scatter", mode="lines",
              line=list(color="#9E9E9E", width=2)) %>%
      layout(title="Natural History — Placebo Body Weight (%)",
             xaxis=list(title="Week"),
             yaxis=list(title="Weight Change (%)"),
             paper_bgcolor="#FAFAFA")
    pl
  })

  output$patient_tbl <- DT::renderDataTable({
    data.frame(
      Parameter = c("Body Weight","Height","BMI","Obesity Class",
                    "Fasting Glucose","Fasting Insulin","HOMA-IR",
                    "HbA1c","Triglycerides","Adipose Mass"),
      Value     = c(
        paste0(input$bw0," kg"),
        paste0(input$height," m"),
        sprintf("%.1f kg/m²", input$bw0/input$height^2),
        {bmi<-input$bw0/input$height^2; if(bmi<25)"Normal" else if(bmi<30)"Overweight" else if(bmi<35)"Class I Obese" else if(bmi<40)"Class II Obese" else "Class III Obese"},
        paste0(input$g0," mg/dL"),
        paste0(input$i0," µU/mL"),
        sprintf("%.1f", input$i0*input$g0/22.5),
        paste0(input$hba0," %"),
        paste0(input$tg0," mg/dL"),
        paste0(input$a0," kg")
      )
    )
  }, options=list(dom="t", pageLength=10))

  # ── Tab 2: Drug PK ───────────────────────────────────────────────────────────
  output$pk_conc_plot <- renderPlotly({
    df <- all_df()
    if(is.null(df)) return(NULL)
    arm <- input$pk_drug
    d   <- filter(df, .data$arm==arm)
    conc_col <- switch(arm, sema="SEMA_C", tirz="TIRZ_C",
                             orlistat="ORL_GUT", cns="CNS_C")
    y_lab <- switch(arm, sema="Semaglutide (nmol/L)", tirz="Tirzepatide (nmol/L)",
                    orlistat="Orlistat Gut Lumen (mg)", cns="CNS Agent (rel.)")
    y_val <- d[[conc_col]]
    if(input$pk_log) y_val <- log10(pmax(y_val, 1e-6))
    ylab  <- if(input$pk_log) paste0("log10 ", y_lab) else y_lab

    plot_ly(d, x=~week, y=y_val, type="scatter", mode="lines",
            line=list(color=arm_meta[[arm]]$color, width=2)) %>%
      layout(title=paste("PK:", arm_meta[[arm]]$label),
             xaxis=list(title="Week"), yaxis=list(title=ylab),
             paper_bgcolor="#FAFAFA")
  })

  output$pk_rec_occ_plot <- renderPlotly({
    df <- all_df()
    if(is.null(df)) return(NULL)
    arm <- input$pk_drug
    d   <- filter(df, .data$arm==arm)
    if(!arm %in% c("sema","tirz")) return(NULL)

    plot_ly(d, x=~week, y=~GLP1R_OCC*100, type="scatter", mode="lines",
            name="GLP-1R (%)", line=list(color="#0288D1",width=2)) %>%
      add_trace(y=~GIPR_OCC*100, name="GIPR (%)",
                line=list(color="#F9A825",width=2)) %>%
      layout(title="Receptor Occupancy (%) — Escalation Phase",
             xaxis=list(title="Week"), yaxis=list(title="Occupancy (%)",range=c(0,100)),
             paper_bgcolor="#FAFAFA")
  })

  output$pk_ref_tbl <- renderTable({
    tribble(
      ~Param, ~Sema, ~Tirz,
      "MW (g/mol)", "4114", "4814",
      "t½ (h)",     "168",  "120",
      "F (%)",      "89",   "80",
      "Vd (L)",     "12.4", "9.2",
      "EC50 (nM)",  "0.016","0.050/0.013"
    )
  })

  # ── Tab 3: PD Biomarkers ─────────────────────────────────────────────────────
  pd_data <- reactive({
    df <- all_df()
    if(is.null(df)) return(NULL)
    df %>%
      filter(arm %in% input$arms_pd,
             week >= input$pd_week_range[1],
             week <= input$pd_week_range[2])
  })

  make_pd_plot <- function(var, title, ylab, ref_line=NULL) {
    renderPlotly({
      d <- pd_data()
      if(is.null(d)) return(NULL)
      arms_shown <- unique(d$arm)
      colors <- setNames(sapply(arms_shown, function(a) arm_meta[[a]]$color), arms_shown)
      labels <- setNames(sapply(arms_shown, function(a) arm_meta[[a]]$label), arms_shown)

      p <- plot_ly()
      for(a in arms_shown) {
        da <- filter(d, arm==a)
        p  <- add_trace(p, data=da, x=~week, y=~.data[[var]],
                        type="scatter", mode="lines",
                        name=labels[[a]], line=list(color=colors[[a]], width=1.8))
      }
      if(!is.null(ref_line))
        p <- add_segments(p, x=0, xend=max(d$week,na.rm=TRUE),
                          y=ref_line, yend=ref_line,
                          line=list(dash="dash",color="red",width=1),
                          showlegend=FALSE)
      layout(p, title=title, xaxis=list(title="Week"),
             yaxis=list(title=ylab), showlegend=TRUE,
             paper_bgcolor="#FAFAFA")
    })
  }

  output$pd_insulin  <- make_pd_plot("INSULIN_P", "Plasma Insulin (µU/mL)","Insulin",ref_line=15)
  output$pd_glucose  <- make_pd_plot("GLUCOSE_P", "Fasting Glucose (mg/dL)","Glucose",ref_line=100)
  output$pd_food     <- make_pd_plot("FOOD_R",    "Food Intake (rel.)","Relative Intake",ref_line=0.7)
  output$pd_ghrelin  <- make_pd_plot("GHRELIN_R", "Ghrelin (rel.)","Relative Ghrelin")
  output$pd_leptin   <- make_pd_plot("LEPTIN_P",  "Plasma Leptin (ng/mL)","Leptin")
  output$pd_inflam   <- make_pd_plot("INFLAM_I",  "Inflammation Index","Inflam. Index",ref_line=1.0)

  # ── Tab 4: Clinical Endpoints ─────────────────────────────────────────────────
  ep_data <- reactive({
    df <- all_df()
    if(is.null(df)) return(NULL)
    filter(df, arm %in% input$arms_ep)
  })

  make_ep_plot <- function(var, title, ylab, ref_lines=NULL) {
    renderPlotly({
      d <- ep_data()
      if(is.null(d)) return(NULL)
      arms_shown <- unique(d$arm)
      colors <- setNames(sapply(arms_shown, function(a) arm_meta[[a]]$color), arms_shown)
      labels <- setNames(sapply(arms_shown, function(a) arm_meta[[a]]$label), arms_shown)

      p <- plot_ly()
      for(a in arms_shown) {
        da <- filter(d, arm==a)
        yv <- if(var=="bmi") da$bmi else if(var=="waist") 88.2+1.15*(da$BWT_C-106) else da[[var]]
        p  <- add_trace(p, x=da$week, y=yv, type="scatter", mode="lines",
                        name=labels[[a]], line=list(color=colors[[a]], width=1.8))
      }
      if(!is.null(ref_lines)) {
        for(rl in ref_lines)
          p <- add_segments(p, x=0, xend=max(d$week,na.rm=TRUE),
                            y=rl$y, yend=rl$y,
                            line=list(dash="dash", color=rl$col, width=1),
                            showlegend=FALSE)
      }
      layout(p, title=title, xaxis=list(title="Week"),
             yaxis=list(title=ylab), paper_bgcolor="#FAFAFA")
    })
  }

  output$ep_bwt   <- make_ep_plot("bwt_pct","Body Weight Change (%)","% from Baseline",
                                   list(list(y=-5,col="orange"),list(y=-10,col="red"),list(y=-15,col="darkred")))
  output$ep_bmi   <- make_ep_plot("bmi","BMI (kg/m²)","BMI",
                                   list(list(y=30,col="orange"),list(y=25,col="green")))
  output$ep_hba1c <- make_ep_plot("HBA1C_C","HbA1c (%)","HbA1c (%)",
                                   list(list(y=5.7,col="orange"),list(y=6.5,col="red")))
  output$ep_trig  <- make_ep_plot("TRIG_P","Triglycerides (mg/dL)","TG (mg/dL)",
                                   list(list(y=150,col="orange"),list(y=200,col="red")))
  output$ep_homa  <- make_ep_plot("HOMA_IR_C","HOMA-IR","HOMA-IR",
                                   list(list(y=2.5,col="red")))
  output$ep_waist <- make_ep_plot("waist","Waist Circumference (est. cm)","Waist (cm)",
                                   list(list(y=88,col="orange"),list(y=80,col="green")))

  # ── Tab 5: Scenario Comparison ────────────────────────────────────────────────
  cmp_data <- reactive({
    df <- all_df()
    if(is.null(df)) return(NULL)
    target_wk <- input$cmp_week
    df %>%
      group_by(arm, arm_label, arm_color) %>%
      arrange(week) %>%
      filter(abs(week-target_wk)==min(abs(week-target_wk))) %>%
      slice(1) %>%
      ungroup()
  })

  output$cmp_bar <- renderPlotly({
    d <- cmp_data()
    if(is.null(d)) return(NULL)
    var <- input$cmp_var
    yv  <- d[[var]]
    y0  <- switch(var, bwt_pct=0, GLUCOSE_P=input$g0, HBA1C_C=input$hba0,
                  TRIG_P=input$tg0, HOMA_IR_C=input$i0*input$g0/22.5,
                  LEPTIN_P=32, INFLAM_I=1.5)
    delta <- yv - y0

    plot_ly(d, x=~arm_label, y=delta, type="bar",
            marker=list(color=d$arm_color)) %>%
      layout(title=paste("Change from Baseline —", var, "at Week", input$cmp_week),
             xaxis=list(title=""),
             yaxis=list(title=paste("Δ", var)),
             paper_bgcolor="#FAFAFA")
  })

  output$cmp_spider <- renderPlotly({
    d <- cmp_data()
    if(is.null(d)) return(NULL)

    vars  <- c("bwt_pct","GLUCOSE_P","HBA1C_C","TRIG_P","HOMA_IR_C","LEPTIN_P")
    vlabs <- c("WT%","Glucose","HbA1c","TG","HOMA-IR","Leptin")
    baselines <- c(0, input$g0, input$hba0, input$tg0, input$i0*input$g0/22.5, 32)
    # Normalize: lower is better → invert for radar (lower = closer to center)

    p <- plot_ly(type="scatterpolar", fill="toself")
    for(i in seq_len(nrow(d))) {
      row_vals <- as.numeric(d[i, vars])
      rel_vals <- 100 * (1 - (row_vals - baselines) / (baselines + 1e-9))
      rel_vals <- pmax(pmin(rel_vals, 200), 0)
      p <- add_trace(p, r=c(rel_vals, rel_vals[1]),
                     theta=c(vlabs, vlabs[1]),
                     name=d$arm_label[i],
                     line=list(color=d$arm_color[i]))
    }
    layout(p, polar=list(radialaxis=list(range=c(0,150))),
           title=paste("Relative Improvement at Week", input$cmp_week),
           paper_bgcolor="#FAFAFA")
  })

  output$cmp_table <- DT::renderDataTable({
    d <- cmp_data()
    if(is.null(d)) return(NULL)
    d %>%
      select(arm_label, BWT_C, bwt_pct, bmi, GLUCOSE_P, HBA1C_C,
             TRIG_P, HOMA_IR_C, LEPTIN_P, INFLAM_I) %>%
      mutate(across(where(is.numeric), ~round(.,2))) %>%
      rename(Arm=arm_label, "Wt(kg)"=BWT_C, "Wt%Δ"=bwt_pct,
             BMI=bmi, "Glucose"=GLUCOSE_P, "HbA1c%"=HBA1C_C,
             "TG(mg/dL)"=TRIG_P, "HOMA-IR"=HOMA_IR_C,
             "Leptin"=LEPTIN_P, "Inflam"=INFLAM_I)
  }, options=list(scrollX=TRUE, dom="t"))

  output$dl_cmp <- downloadHandler(
    filename=function() paste0("obesity_qsp_comparison_wk",input$cmp_week,".csv"),
    content=function(f) {
      d <- cmp_data()
      write.csv(d, f, row.names=FALSE)
    }
  )

  # ── Tab 6: Cardiometabolic Risk ───────────────────────────────────────────────
  risk_arm_data <- reactive({
    df <- all_df()
    if(is.null(df)) return(NULL)
    arm <- input$risk_drug
    filter(df, .data$arm==arm)
  })

  # Pooled Cohort Eq. simplified (10-yr ASCVD %)
  ascvd_risk <- function(age, sex, sbp, tc, hdl, smoker, dm) {
    # White female formula (simplified Goff 2014)
    lp <- -29.799 + 4.884*log(age) + 13.540*log(tc) - 3.114*log(hdl) +
           2.019*log(sbp) + 7.574*smoker + 0.661*dm
    100*(1 - 0.9665^exp(lp - 26.114))
  }

  output$risk_wt <- renderText({
    d <- risk_arm_data()
    if(is.null(d)) return("—")
    tw <- input$risk_week
    row <- d[which.min(abs(d$week-tw)),]
    sprintf("%.1f%%", row$bwt_pct)
  })

  output$risk_ascvd <- renderText({
    d <- risk_arm_data()
    if(is.null(d)) return("—")
    tw <- input$risk_week
    row <- d[which.min(abs(d$week-tw)),]
    # ASCVD assumes SBP improves by ~0.5 mmHg per kg lost
    delta_wt <- row$BWT_C - input$bw0
    sbp_now  <- 135 + delta_wt*0.5
    tc       <- input$ldl0 + input$hdl0 + row$TRIG_P/5
    dm       <- if(row$HBA1C_C >= 6.5) 1 else 0
    r0 <- ascvd_risk(input$age, input$sex, 135, input$ldl0+input$hdl0+input$tg0/5,
                     input$hdl0, as.integer(input$smoker), dm)
    r1 <- ascvd_risk(input$age, input$sex, sbp_now, tc,
                     input$hdl0, as.integer(input$smoker), dm)
    sprintf("%.1f → %.1f%%", r0, r1)
  })

  output$risk_hba  <- renderText({
    d <- risk_arm_data()
    if(is.null(d)) return("—")
    tw  <- input$risk_week
    row <- d[which.min(abs(d$week-tw)),]
    sprintf("%.2f → %.2f%%", input$hba0, row$HBA1C_C)
  })

  output$risk_tg <- renderText({
    d <- risk_arm_data()
    if(is.null(d)) return("—")
    tw  <- input$risk_week
    row <- d[which.min(abs(d$week-tw)),]
    sprintf("%.0f → %.0f mg/dL", input$tg0, row$TRIG_P)
  })

  output$risk_scatter <- renderPlotly({
    df <- all_df()
    if(is.null(df)) return(NULL)
    tw <- input$risk_week
    d  <- df %>%
      group_by(arm, arm_label, arm_color) %>%
      arrange(week) %>% filter(abs(week-tw)==min(abs(week-tw))) %>%
      slice(1) %>% ungroup()
    plot_ly(d, x=~bwt_pct, y=~HOMA_IR_C, type="scatter", mode="markers+text",
            text=~arm_label, textposition="top center",
            marker=list(color=d$arm_color, size=16, line=list(color="white",width=2))) %>%
      layout(title=paste("Weight Loss vs HOMA-IR at Week", tw),
             xaxis=list(title="Body Weight Change (%)"),
             yaxis=list(title="HOMA-IR"),
             paper_bgcolor="#FAFAFA")
  })

  output$dr_curve <- renderPlotly({
    # Dose-response: tirzepatide 5 doses, semaglutide 5 doses
    dr_arm <- input$dr_drug
    tw     <- input$dr_weeks
    doses  <- if(dr_arm=="sema") c(0.5, 1.0, 1.7, 2.0, 2.4) else c(5, 7.5, 10, 12.5, 15)
    MW     <- if(dr_arm=="sema") 4113.6 else 4813.5
    cmt    <- if(dr_arm=="sema") "SEMA_GUT" else "TIRZ_GUT"

    results <- lapply(doses, function(d_mg) {
      nmol  <- d_mg / MW * 1e6
      times <- seq(0, (tw-1)*7*24, by=7*24)
      ev_df <- data.frame(var=cmt, time=times, value=nmol, method="add")
      y0 <- c(SEMA_GUT=0, SEMA_C=0, TIRZ_GUT=0, TIRZ_C=0, ORL_GUT=0, CNS_C=0,
               GLP1R_OCC=0, GIPR_OCC=0, FOOD_R=1.0, GASTRIC=1.0, GHRELIN_R=1.0,
               INSULIN_P=input$i0, GLUCOSE_P=input$g0,
               ADIP=input$a0, BWT_C=input$bw0, LEPTIN_P=32,
               TRIG_P=input$tg0, HBA1C_C=input$hba0,
               INFLAM_I=1.5, HOMA_IR_C=input$i0*input$g0/22.5)
      pars <- c(default_pars, I0=input$i0, G0=input$g0, A0=input$a0,
                BW0=input$bw0, TG0=input$tg0, INF0=1.5, L0=32, HOMA0=input$i0*input$g0/22.5)
      tryCatch({
        out <- ode(y=y0, times=c(0, tw*7*24), func=obesity_ode,
                   parms=pars, method="lsoda",
                   events=list(data=ev_df))
        df_o <- as.data.frame(out)
        bwt_pct_end <- (df_o$BWT_C[nrow(df_o)] - input$bw0)/input$bw0*100
        data.frame(dose=d_mg, bwt_pct=bwt_pct_end)
      }, error=function(e) data.frame(dose=d_mg, bwt_pct=NA))
    }) %>% bind_rows()

    col <- arm_meta[[dr_arm]]$color
    plot_ly(results, x=~dose, y=~bwt_pct, type="scatter", mode="lines+markers",
            line=list(color=col, width=2.5),
            marker=list(color=col, size=10)) %>%
      layout(title=paste("Dose-Response:", arm_meta[[dr_arm]]$label, "at Week", tw),
             xaxis=list(title=paste0("Dose (mg)")),
             yaxis=list(title="Body Weight Change (%)"),
             paper_bgcolor="#FAFAFA")
  })

  output$risk_tbl <- DT::renderDataTable({
    df <- all_df()
    if(is.null(df)) return(NULL)
    tw <- input$risk_week
    d  <- df %>%
      group_by(arm, arm_label) %>%
      arrange(week) %>% filter(abs(week-tw)==min(abs(week-tw))) %>%
      slice(1) %>% ungroup() %>%
      mutate(
        wt_loss_pct = round(bwt_pct, 1),
        glucose_chg = round(GLUCOSE_P - input$g0, 1),
        hba1c_chg   = round(HBA1C_C  - input$hba0, 2),
        tg_chg      = round(TRIG_P   - input$tg0, 1),
        homa_chg    = round(HOMA_IR_C - input$i0*input$g0/22.5, 2),
        leptin_chg  = round(LEPTIN_P  - 32, 1)
      ) %>%
      select(arm_label, wt_loss_pct, glucose_chg, hba1c_chg, tg_chg, homa_chg, leptin_chg) %>%
      rename(Arm="arm_label", "Wt%Δ"="wt_loss_pct",
             "ΔGlucose"="glucose_chg", "ΔHbA1c"="hba1c_chg",
             "ΔTG"="tg_chg","ΔHOMA-IR"="homa_chg","ΔLeptin"="leptin_chg")
  }, options=list(dom="t", scrollX=TRUE))
}

# ── Launch ─────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
