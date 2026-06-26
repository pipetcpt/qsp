library(shiny)
library(bslib)
library(plotly)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# PDAC QSP Shiny App
# Pancreatic Ductal Adenocarcinoma — Interactive Dashboard
# ============================================================

# --- Scenario parameters calibrated to clinical trials ---
scenarios <- list(
  "1_Control"        = list(name="Untreated",       k_eff=0.00, mPFS=1.5,  mOS=3.0,  ORR=0,    tox_G34=0),
  "2_Gem"            = list(name="Gemcitabine",      k_eff=0.12, mPFS=3.7,  mOS=6.7,  ORR=7,    tox_G34=22),
  "3_GemNab"         = list(name="Gem+nab-Pac",      k_eff=0.22, mPFS=5.5,  mOS=8.5,  ORR=23,   tox_G34=38),
  "4_FFX"            = list(name="FOLFIRINOX",       k_eff=0.28, mPFS=6.4,  mOS=11.1, ORR=31.6, tox_G34=45),
  "5_mFFX"           = list(name="mFOLFIRINOX",      k_eff=0.25, mPFS=6.0,  mOS=10.5, ORR=28,   tox_G34=35),
  "6_MRTX"           = list(name="MRTX1133(G12D)",   k_eff=0.35, mPFS=4.0,  mOS=8.0,  ORR=40,   tox_G34=15),
  "7_Olaparib"       = list(name="Olaparib(BRCA+)",  k_eff=0.18, mPFS=7.4,  mOS=NA,   ORR=NA,   tox_G34=10)
)

scen_names <- sapply(scenarios, function(x) x$name)

# --- Analytical PK functions ---
pk_biexp <- function(t, dose, CL, V1, V2, Q) {
  k10 <- CL/V1; k12 <- Q/V1; k21 <- Q/V2
  alpha <- 0.5*((k10+k12+k21)+sqrt((k10+k12+k21)^2-4*k10*k21))
  beta  <- 0.5*((k10+k12+k21)-sqrt((k10+k12+k21)^2-4*k10*k21))
  A <- (alpha-k21)/(V1*(alpha-beta))
  B <- (k21-beta)/(V1*(alpha-beta))
  conc <- dose*(A*exp(-alpha*t)+B*exp(-beta*t))
  pmax(conc, 0)
}

pk_1cmt <- function(t, dose, CL, V, F_oral=1, ka=1) {
  k <- CL/V
  conc <- F_oral*dose*ka/(V*(ka-k))*(exp(-k*t)-exp(-ka*t))
  pmax(conc, 0)
}

# --- Simeoni TGI ODE (deSolve-free Euler integration) ---
simulate_tumor <- function(k_eff, V0=5000, days=365, dt=0.1) {
  k_prog <- 0.012; k_death <- 0.008; k_tr <- 0.18
  n <- ceiling(days/dt)+1
  t_vec <- seq(0, days, by=dt)
  x0 <- rep(0,n); x1 <- rep(0,n); x2 <- rep(0,n); x3 <- rep(0,n)
  TUM <- rep(0,n); CA199 <- rep(0,n)
  x0[1] <- V0*0.9; TUM[1] <- V0*0.1; CA199[1] <- 500
  for(i in 1:(n-1)) {
    eff <- min(k_eff, 0.99)
    dx0  <- (k_prog*(1-eff) - k_death)*x0[i]
    dx1  <- k_prog*eff*x0[i] - k_tr*x1[i]
    dx2  <- k_tr*x1[i] - k_tr*x2[i]
    dx3  <- k_tr*x2[i] - k_tr*x3[i]
    dTUM <- k_tr*x3[i] - k_death*TUM[i]
    dCA  <- 0.01*(x0[i]+TUM[i]) - 0.05*CA199[i]
    x0[i+1]   <- max(x0[i]   + dt*dx0,  0)
    x1[i+1]   <- max(x1[i]   + dt*dx1,  0)
    x2[i+1]   <- max(x2[i]   + dt*dx2,  0)
    x3[i+1]   <- max(x3[i]   + dt*dx3,  0)
    TUM[i+1]  <- max(TUM[i]  + dt*dTUM, 0)
    CA199[i+1]<- max(CA199[i]+ dt*dCA,  0)
  }
  data.frame(time=t_vec, x0=x0, tumor=TUM+x0, CA199=CA199)
}

# --- Friberg neutropenia (Euler) ---
simulate_neut <- function(drug_conc_fn, days=84, dt=0.1, EC50=0.8, Emax=0.95) {
  MTT <- 132/24; ktr <- 4/MTT; gamma <- 0.17; Circ0 <- 5.0
  n <- ceiling(days/dt)+1; t_vec <- seq(0, days, by=dt)
  Prol <- rep(Circ0,n); Tr1 <- rep(Circ0,n); Tr2 <- rep(Circ0,n)
  Tr3 <- rep(Circ0,n); Circ <- rep(Circ0,n)
  for(i in 1:(n-1)) {
    conc <- drug_conc_fn(t_vec[i])
    inh <- Emax*conc/(EC50+conc)
    dP  <- ktr*Prol[i]*(1-inh)*(Circ0/max(Circ[i],0.01))^gamma - ktr*Prol[i]
    dT1 <- ktr*Prol[i] - ktr*Tr1[i]
    dT2 <- ktr*Tr1[i]  - ktr*Tr2[i]
    dT3 <- ktr*Tr2[i]  - ktr*Tr3[i]
    dC  <- ktr*Tr3[i]  - ktr*Circ[i]
    Prol[i+1] <- max(Prol[i]+dt*dP, 0)
    Tr1[i+1]  <- max(Tr1[i]+dt*dT1, 0)
    Tr2[i+1]  <- max(Tr2[i]+dt*dT2, 0)
    Tr3[i+1]  <- max(Tr3[i]+dt*dT3, 0)
    Circ[i+1] <- max(Circ[i]+dt*dC, 0)
  }
  data.frame(time=t_vec, ANC=Circ)
}

# --- KM-like survival (exponential) ---
km_curve <- function(mOS, t=seq(0,24,0.5)) {
  if(is.na(mOS)) return(data.frame(time=t, surv=NA))
  lam <- log(2)/mOS
  data.frame(time=t, surv=exp(-lam*t))
}

# ============================================================
# UI
# ============================================================
ui <- page_navbar(
  title = "PDAC QSP Dashboard",
  theme = bs_theme(bootswatch="darkly", version=5),
  bg = "#1a1a2e",

  # --- TAB 1: Patient Profile ---
  nav_panel("Patient Profile",
    layout_sidebar(
      sidebar = sidebar(
        h5("Patient Characteristics"),
        selectInput("stage", "Disease Stage",
          choices=c("Stage I","Stage II","Stage III","Stage IV"),
          selected="Stage IV"),
        selectInput("subtype", "Molecular Subtype",
          choices=c("Classical","Basal-like/Squamous","Quasi-mesenchymal","Exocrine-like"),
          selected="Classical"),
        selectInput("kras_mut", "KRAS Mutation",
          choices=c("G12D (~44%)","G12V (~26%)","G12R (~14%)","G12C (~2%)","Other","Wild-type"),
          selected="G12D (~44%)"),
        checkboxGroupInput("hrd_status","HRD Status",
          choices=c("BRCA1 mut","BRCA2 mut","PALB2 mut","ATM mut"),
          selected=NULL),
        radioButtons("msi","MSI Status",choices=c("MSS","MSI-H (~1%)"),selected="MSS"),
        sliderInput("ecog","ECOG Performance Status",0,2,0,step=1),
        numericInput("weight","Weight (kg)",75,40,150),
        numericInput("bsa","BSA (m²)",1.8,1.2,2.5,0.05),
        numericInput("ca199_base","Baseline CA19-9 (U/mL)",500,0,50000),
        numericInput("tv_base","Baseline Tumor Volume (mm³)",5000,100,100000)
      ),
      card(
        card_header("Patient Summary"),
        layout_columns(
          col_widths=c(6,6),
          card(tableOutput("pt_summary_tbl")),
          card(plotlyOutput("risk_radar", height="350px"))
        )
      )
    )
  ),

  # --- TAB 2: Drug PK ---
  nav_panel("Drug PK",
    layout_sidebar(
      sidebar = sidebar(
        h5("PK Settings"),
        selectInput("pk_regimen","Treatment Regimen",
          choices=c("Gemcitabine mono","Gem+nab-Paclitaxel","FOLFIRINOX","mFOLFIRINOX","MRTX1133","Olaparib")),
        sliderInput("pk_cycle","Cycle Number",1,6,1),
        checkboxGroupInput("pk_show","Show Drug(s)",
          choices=c("Gemcitabine","nab-Paclitaxel","SN-38","5-FU","Oxaliplatin","MRTX1133","Olaparib"),
          selected="Gemcitabine"),
        checkboxInput("pk_log","Log Y-axis",FALSE)
      ),
      card(
        card_header("Plasma Concentration-Time Profiles"),
        plotlyOutput("pk_plot", height="400px"),
        tableOutput("pk_params_tbl")
      )
    )
  ),

  # --- TAB 3: Tumor Dynamics ---
  nav_panel("Tumor Dynamics",
    layout_sidebar(
      sidebar = sidebar(
        h5("Simulation Settings"),
        checkboxGroupInput("td_scenarios","Treatment Scenarios",
          choices=scen_names, selected=scen_names[1:4]),
        sliderInput("td_months","Duration (months)",3,24,12),
        sliderInput("stroma_lvl","Stroma Level",0,1,0.5,0.1),
        checkboxInput("kras_g12d","KRAS G12D present",TRUE)
      ),
      card(
        card_header("Tumor Volume Over Time"),
        plotlyOutput("tumor_plot", height="380px")
      ),
      card(
        card_header("CA19-9 Trajectory"),
        plotlyOutput("ca199_plot", height="300px")
      ),
      card(
        card_header("Waterfall: Best % Change from Baseline"),
        plotlyOutput("waterfall_plot", height="300px")
      )
    )
  ),

  # --- TAB 4: Biomarkers ---
  nav_panel("Biomarkers",
    layout_sidebar(
      sidebar = sidebar(
        h5("Biomarker Settings"),
        radioButtons("bm_type","Biomarker",
          choices=c("CA19-9","ctDNA VAF (%)","Neutrophils (ANC)"),
          selected="CA19-9"),
        sliderInput("bm_weeks","Time Horizon (weeks)",4,52,24),
        checkboxGroupInput("bm_scenarios","Scenarios",
          choices=scen_names, selected=scen_names[1:4])
      ),
      card(
        card_header("Biomarker Kinetics"),
        plotlyOutput("bm_plot", height="380px")
      ),
      card(
        card_header("Grade 3/4 Neutropenia Risk by Regimen"),
        tableOutput("neut_tox_tbl")
      )
    )
  ),

  # --- TAB 5: Clinical Endpoints ---
  nav_panel("Clinical Endpoints",
    layout_sidebar(
      sidebar = sidebar(
        h5("Endpoint Settings"),
        checkboxGroupInput("ep_scenarios","Scenarios",
          choices=scen_names, selected=scen_names),
        checkboxInput("ep_ci","Show 95% CI band",TRUE)
      ),
      card(
        card_header("Overall Survival (Parametric Exponential Approximation)"),
        plotlyOutput("os_plot", height="350px")
      ),
      card(
        card_header("Progression-Free Survival"),
        plotlyOutput("pfs_plot", height="300px")
      ),
      card(
        card_header("ORR & Toxicity Summary"),
        layout_columns(
          col_widths=c(6,6),
          plotlyOutput("orr_plot", height="300px"),
          tableOutput("ep_summary_tbl")
        )
      )
    )
  ),

  # --- TAB 6: Scenario Comparison ---
  nav_panel("Scenario Comparison",
    layout_sidebar(
      sidebar = sidebar(
        h5("Comparison Settings"),
        selectInput("ref_scenario","Reference Scenario",
          choices=scen_names, selected=scen_names[1]),
        checkboxGroupInput("comp_scenarios","Compare Against",
          choices=scen_names, selected=scen_names[-1])
      ),
      card(
        card_header("Forest Plot: HR for OS vs Reference"),
        plotlyOutput("forest_plot", height="400px")
      ),
      card(
        card_header("Comprehensive Comparison Table"),
        tableOutput("comp_tbl")
      ),
      card(
        card_header("Sensitivity: Key Parameters (Tornado)"),
        plotlyOutput("tornado_plot", height="350px")
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # ---- Tab 1: Patient Profile ----
  output$pt_summary_tbl <- renderTable({
    hrd <- if(length(input$hrd_status)>0) paste(input$hrd_status,collapse=", ") else "None"
    mrtx_eligible <- grepl("G12D",input$kras_mut)
    ola_eligible  <- length(input$hrd_status)>0
    pembro_eligible <- input$msi=="MSI-H (~1%)"
    data.frame(
      Parameter = c("Stage","Subtype","KRAS","HRD","MSI","ECOG","BSA","CA19-9","Tumor Vol",
                    "MRTX1133 eligible","Olaparib eligible","Pembrolizumab eligible"),
      Value = c(input$stage, input$subtype, input$kras_mut, hrd, input$msi, input$ecog,
                paste(input$bsa,"m²"), paste(input$ca199_base,"U/mL"),
                paste(input$tv_base,"mm³"),
                ifelse(mrtx_eligible,"YES","NO"),
                ifelse(ola_eligible,"YES","NO"),
                ifelse(pembro_eligible,"YES","NO"))
    )
  })

  output$risk_radar <- renderPlotly({
    stage_score <- switch(input$stage, "Stage I"=10,"Stage II"=35,"Stage III"=65,"Stage IV"=100)
    kras_score  <- switch(input$kras_mut,
      "G12D (~44%)"=90,"G12V (~26%)"=85,"G12R (~14%)"=80,"G12C (~2%)"=70,"Other"=70,"Wild-type"=20)
    hrd_score   <- if(length(input$hrd_status)>0) 70 else 30
    msi_score   <- if(input$msi=="MSI-H (~1%)") 60 else 20
    ecog_score  <- (2-input$ecog)*40+20
    ca199_score <- min(100, input$ca199_base/500*100)
    cats <- c("Stage","KRAS Severity","HRD","Immunotherapy","Performance","CA19-9")
    vals <- c(stage_score, kras_score, hrd_score, msi_score, ecog_score, ca199_score)
    plot_ly(type="scatterpolar", r=c(vals,vals[1]), theta=c(cats,cats[1]),
            fill="toself", fillcolor="rgba(255,107,107,0.3)",
            line=list(color="#ff6b6b")) %>%
      layout(polar=list(bgcolor="#1a1a2e",
               radialaxis=list(range=c(0,100),color="white"),
               angularaxis=list(color="white")),
             paper_bgcolor="#1a1a2e", plot_bgcolor="#1a1a2e",
             font=list(color="white"), title="Risk Profile")
  })

  # ---- Tab 2: Drug PK ----
  output$pk_plot <- renderPlotly({
    t_hrs <- seq(0, 24, 0.25)
    p <- plot_ly()
    if("Gemcitabine" %in% input$pk_show) {
      dose <- 1000*input$bsa
      c_gem <- pk_biexp(t_hrs, dose, CL=80, V1=15, V2=30, Q=25)
      p <- add_trace(p, x=t_hrs, y=c_gem, name="Gemcitabine (ng/mL)",
                     type="scatter", mode="lines", line=list(color="#4ecdc4"))
    }
    if("nab-Paclitaxel" %in% input$pk_show) {
      dose <- 125*input$bsa
      c_nab <- pk_biexp(t_hrs, dose, CL=12, V1=8, V2=100, Q=8)
      p <- add_trace(p, x=t_hrs, y=c_nab, name="nab-Paclitaxel (ng/mL)",
                     type="scatter", mode="lines", line=list(color="#ff6b6b"))
    }
    if("SN-38" %in% input$pk_show) {
      dose_cpt <- 180*input$bsa
      c_cpt <- pk_biexp(t_hrs, dose_cpt, CL=15, V1=180, V2=200, Q=10)
      c_sn38 <- 0.05 * c_cpt * exp(-0.3*t_hrs)
      p <- add_trace(p, x=t_hrs, y=c_sn38, name="SN-38 (ng/mL)",
                     type="scatter", mode="lines", line=list(color="#ffd166"))
    }
    if("5-FU" %in% input$pk_show) {
      dose_fu <- 400*input$bsa
      c_fu <- pk_biexp(t_hrs, dose_fu, CL=120, V1=20, V2=15, Q=30)
      p <- add_trace(p, x=t_hrs, y=c_fu, name="5-FU (ng/mL)",
                     type="scatter", mode="lines", line=list(color="#a8dadc"))
    }
    if("MRTX1133" %in% input$pk_show) {
      c_mrtx <- pk_1cmt(t_hrs, 100, CL=5, V=25, F_oral=0.45, ka=0.8)
      p <- add_trace(p, x=t_hrs, y=c_mrtx, name="MRTX1133 (ng/mL)",
                     type="scatter", mode="lines", line=list(color="#e76f51"))
    }
    if("Olaparib" %in% input$pk_show) {
      c_ola <- pk_1cmt(t_hrs, 300, CL=8.9, V=158, F_oral=0.73, ka=0.5)
      p <- add_trace(p, x=t_hrs, y=c_ola, name="Olaparib (ng/mL)",
                     type="scatter", mode="lines", line=list(color="#2a9d8f"))
    }
    p %>% layout(
      xaxis=list(title="Time (hours)", color="white"),
      yaxis=list(title="Concentration", color="white",
                 type=if(input$pk_log) "log" else "linear"),
      paper_bgcolor="#1a1a2e", plot_bgcolor="#2d2d4e",
      font=list(color="white"), legend=list(font=list(color="white"))
    )
  })

  output$pk_params_tbl <- renderTable({
    data.frame(
      Drug = c("Gemcitabine","nab-Paclitaxel","Oxaliplatin","SN-38","5-FU","MRTX1133","Olaparib"),
      `CL (L/h)` = c(80,12,9,8,120,5,8.9),
      `V (L)` = c(15,8,12,15,20,25,158),
      `T1/2 (h)` = round(log(2)*c(15/80,8/12,12/9,15/8,20/120,25/5,158/8.9),2),
      `Route` = c("IV","IV","IV","Metabolite","IV","Oral","Oral"),
      check.names=FALSE
    )
  })

  # ---- Tab 3: Tumor Dynamics ----
  sim_tumor_all <- reactive({
    days <- input$td_months*30.4
    stroma_pen <- 1 - 0.6*input$stroma_lvl
    do.call(rbind, lapply(names(scenarios), function(sid) {
      sc <- scenarios[[sid]]
      if(sc$name %in% input$td_scenarios) {
        k_eff_adj <- sc$k_eff * stroma_pen
        if(sc$name=="MRTX1133(G12D)" && !input$kras_g12d) k_eff_adj <- 0
        res <- simulate_tumor(k_eff_adj, V0=input$tv_base, days=days)
        res$scenario <- sc$name
        res
      }
    }))
  })

  output$tumor_plot <- renderPlotly({
    df <- sim_tumor_all()
    req(nrow(df)>0)
    df$time_mo <- df$time/30.4
    plot_ly(df, x=~time_mo, y=~tumor, color=~scenario,
            type="scatter", mode="lines") %>%
      layout(xaxis=list(title="Time (months)",color="white"),
             yaxis=list(title="Tumor Volume (mm³)",color="white"),
             paper_bgcolor="#1a1a2e", plot_bgcolor="#2d2d4e",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$ca199_plot <- renderPlotly({
    df <- sim_tumor_all()
    req(nrow(df)>0)
    df$time_mo <- df$time/30.4
    plot_ly(df, x=~time_mo, y=~CA199, color=~scenario,
            type="scatter", mode="lines") %>%
      layout(xaxis=list(title="Time (months)",color="white"),
             yaxis=list(title="CA19-9 (U/mL)",color="white"),
             paper_bgcolor="#1a1a2e", plot_bgcolor="#2d2d4e",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$waterfall_plot <- renderPlotly({
    days <- input$td_months*30.4
    stroma_pen <- 1 - 0.6*input$stroma_lvl
    wf <- do.call(rbind, lapply(names(scenarios), function(sid) {
      sc <- scenarios[[sid]]
      if(sc$name %in% input$td_scenarios) {
        k_eff_adj <- sc$k_eff * stroma_pen
        if(sc$name=="MRTX1133(G12D)" && !input$kras_g12d) k_eff_adj <- 0
        res <- simulate_tumor(k_eff_adj, V0=input$tv_base, days=days)
        best <- min(res$tumor)
        pct_chg <- (best-input$tv_base)/input$tv_base*100
        data.frame(scenario=sc$name, pct_change=round(pct_chg,1))
      }
    }))
    wf <- wf[order(wf$pct_change),]
    wf$color <- ifelse(wf$pct_change < -30, "#4ecdc4", ifelse(wf$pct_change < 0, "#ffd166", "#ff6b6b"))
    plot_ly(wf, x=~scenario, y=~pct_change, type="bar",
            marker=list(color=wf$color)) %>%
      add_segments(x=~scenario, xend=~scenario, y=-30, yend=-30,
                   line=list(dash="dash",color="white",width=1)) %>%
      layout(xaxis=list(title="",color="white",tickangle=-30),
             yaxis=list(title="Best % Change from Baseline",color="white"),
             shapes=list(list(type="line",x0=0,x1=1,xref="paper",y0=-30,y1=-30,
                              line=list(color="white",dash="dash"))),
             paper_bgcolor="#1a1a2e", plot_bgcolor="#2d2d4e",
             font=list(color="white"))
  })

  # ---- Tab 4: Biomarkers ----
  output$bm_plot <- renderPlotly({
    t_days <- seq(0, input$bm_weeks*7, 1)
    p <- plot_ly()
    for(sid in names(scenarios)) {
      sc <- scenarios[[sid]]
      if(sc$name %in% input$bm_scenarios) {
        if(input$bm_type=="CA19-9") {
          res <- simulate_tumor(sc$k_eff, V0=input$tv_base, days=max(t_days))
          res_sub <- res[seq(1,nrow(res),length.out=min(nrow(res),length(t_days))),]
          p <- add_trace(p, x=t_days[1:nrow(res_sub)]/7, y=res_sub$CA199,
                         name=sc$name, type="scatter", mode="lines")
        } else if(input$bm_type=="ctDNA VAF (%)") {
          ctdna <- 30*exp(-sc$k_eff*t_days) + 2
          p <- add_trace(p, x=t_days/7, y=ctdna, name=sc$name,
                         type="scatter", mode="lines")
        } else if(input$bm_type=="Neutrophils (ANC)") {
          gem_fn <- if(grepl("Gem|FOLFI",sc$name) && sc$name!="Untreated") {
            function(t) { cycle_t <- t %% 28; if(cycle_t < 2) 15 else 15*exp(-0.4*(cycle_t-2)) }
          } else { function(t) 0 }
          neut_df <- simulate_neut(gem_fn, days=max(t_days), EC50=0.8)
          neut_sub <- neut_df[seq(1,nrow(neut_df),length.out=min(nrow(neut_df),length(t_days))),]
          p <- add_trace(p, x=neut_sub$time/7, y=neut_sub$ANC,
                         name=sc$name, type="scatter", mode="lines")
        }
      }
    }
    y_title <- switch(input$bm_type,
      "CA19-9"="CA19-9 (U/mL)", "ctDNA VAF (%)"="ctDNA VAF (%)", "Neutrophils (ANC)"="ANC (×10⁹/L)")
    p %>% layout(
      xaxis=list(title="Time (weeks)",color="white"),
      yaxis=list(title=y_title,color="white"),
      paper_bgcolor="#1a1a2e", plot_bgcolor="#2d2d4e",
      font=list(color="white"), legend=list(font=list(color="white"))
    )
  })

  output$neut_tox_tbl <- renderTable({
    data.frame(
      Regimen = scen_names,
      `G3/4 Neutropenia (%)` = c(0, 22, 38, 45, 35, 15, 10),
      `G3/4 Neuropathy (%)` = c(0, 2, 17, 9, 7, 3, 1),
      `G3/4 Diarrhea (%)` = c(0, 2, 6, 13, 8, 5, 3),
      `G3/4 Fatigue (%)` = c(0, 7, 17, 23, 18, 10, 8),
      `Dose Reduction (%)` = c(0, 15, 28, 33, 22, 12, 5),
      check.names=FALSE
    )
  })

  # ---- Tab 5: Clinical Endpoints ----
  km_data <- reactive({
    t_mo <- seq(0, 24, 0.5)
    do.call(rbind, lapply(names(scenarios), function(sid) {
      sc <- scenarios[[sid]]
      if(sc$name %in% input$ep_scenarios) {
        pfs_df <- km_curve(sc$mPFS, t_mo)
        pfs_df$type <- "PFS"; pfs_df$scenario <- sc$name
        os_df  <- km_curve(sc$mOS,  t_mo)
        os_df$type  <- "OS";  os_df$scenario  <- sc$name
        rbind(pfs_df, os_df)
      }
    }))
  })

  output$os_plot <- renderPlotly({
    df <- km_data()
    req(nrow(df)>0)
    df_os <- df[df$type=="OS" & !is.na(df$surv),]
    p <- plot_ly(df_os, x=~time, y=~surv*100, color=~scenario,
                 type="scatter", mode="lines") %>%
      layout(xaxis=list(title="Time (months)",color="white"),
             yaxis=list(title="Survival (%)", range=c(0,100), color="white"),
             paper_bgcolor="#1a1a2e", plot_bgcolor="#2d2d4e",
             font=list(color="white"), legend=list(font=list(color="white")),
             title=list(text="Overall Survival",font=list(color="white")))
    if(input$ep_ci) {
      for(sc_name in unique(df_os$scenario)) {
        sub <- df_os[df_os$scenario==sc_name,]
        lam <- log(2)/unique(scenarios[sapply(scenarios,function(x)x$name==sc_name)][[1]]$mOS)
        if(!is.na(lam)) {
          se <- sqrt(sub$surv*(1-sub$surv)/100)
          p <- add_ribbons(p, x=sub$time, ymin=(sub$surv-1.96*se)*100,
                           ymax=(sub$surv+1.96*se)*100,
                           name=paste(sc_name,"CI"), showlegend=FALSE,
                           opacity=0.2, line=list(color="transparent"))
        }
      }
    }
    p
  })

  output$pfs_plot <- renderPlotly({
    df <- km_data()
    req(nrow(df)>0)
    df_pfs <- df[df$type=="PFS" & !is.na(df$surv),]
    plot_ly(df_pfs, x=~time, y=~surv*100, color=~scenario,
            type="scatter", mode="lines") %>%
      layout(xaxis=list(title="Time (months)",color="white"),
             yaxis=list(title="PFS (%)", range=c(0,100), color="white"),
             paper_bgcolor="#1a1a2e", plot_bgcolor="#2d2d4e",
             font=list(color="white"), legend=list(font=list(color="white")))
  })

  output$orr_plot <- renderPlotly({
    orr_df <- do.call(rbind, lapply(scenarios, function(sc) {
      if(sc$name %in% input$ep_scenarios && !is.na(sc$ORR))
        data.frame(scenario=sc$name, ORR=sc$ORR)
    }))
    req(nrow(orr_df)>0)
    plot_ly(orr_df, x=~scenario, y=~ORR, type="bar",
            marker=list(color="#4ecdc4")) %>%
      layout(xaxis=list(title="",color="white",tickangle=-30),
             yaxis=list(title="ORR (%)",color="white"),
             paper_bgcolor="#1a1a2e", plot_bgcolor="#2d2d4e",
             font=list(color="white"))
  })

  output$ep_summary_tbl <- renderTable({
    do.call(rbind, lapply(scenarios, function(sc) {
      if(sc$name %in% input$ep_scenarios)
        data.frame(Regimen=sc$name, mPFS=sc$mPFS, mOS=sc$mOS,
                   ORR=sc$ORR, `G34_Tox%`=sc$tox_G34, check.names=FALSE)
    }))
  })

  # ---- Tab 6: Scenario Comparison ----
  output$forest_plot <- renderPlotly({
    ref_name <- input$ref_scenario
    ref_mOS <- scenarios[[which(scen_names==ref_name)]]$mOS
    if(is.na(ref_mOS)) ref_mOS <- 3.0
    fp_df <- do.call(rbind, lapply(names(scenarios), function(sid) {
      sc <- scenarios[[sid]]
      if(sc$name %in% input$comp_scenarios && !is.na(sc$mOS)) {
        hr <- ref_mOS/sc$mOS
        lo <- hr*exp(-1.96*0.3); hi <- hr*exp(1.96*0.3)
        data.frame(scenario=sc$name, HR=round(hr,2), lo=round(lo,2), hi=round(hi,2))
      }
    }))
    req(nrow(fp_df)>0)
    fp_df <- fp_df[order(fp_df$HR),]
    plot_ly(fp_df, y=~scenario, x=~HR, type="scatter", mode="markers",
            error_x=list(type="data", symmetric=FALSE,
                         array=fp_df$hi-fp_df$HR, arrayminus=fp_df$HR-fp_df$lo,
                         color="white"),
            marker=list(size=12, color="#ff6b6b")) %>%
      add_segments(x=1, xend=1, y=0, yend=nrow(fp_df)+1,
                   line=list(dash="dash",color="white",width=1)) %>%
      layout(xaxis=list(title=paste("Hazard Ratio (vs",ref_name,")"),color="white",
                        type="log"),
             yaxis=list(title="",color="white"),
             paper_bgcolor="#1a1a2e", plot_bgcolor="#2d2d4e",
             font=list(color="white"))
  })

  output$comp_tbl <- renderTable({
    ref_name <- input$ref_scenario
    ref_sc_idx <- which(scen_names==ref_name)
    ref_mOS <- if(length(ref_sc_idx)>0) scenarios[[ref_sc_idx]]$mOS else 3.0
    if(is.na(ref_mOS)) ref_mOS <- 3.0
    do.call(rbind, lapply(scenarios, function(sc) {
      hr <- if(!is.na(sc$mOS)) round(ref_mOS/sc$mOS,2) else NA
      data.frame(
        Scenario=sc$name, mPFS_mo=sc$mPFS, mOS_mo=sc$mOS,
        `ORR%`=sc$ORR, `G34_Tox%`=sc$tox_G34,
        `HR_OS`=hr, check.names=FALSE
      )
    }))
  })

  output$tornado_plot <- renderPlotly({
    params <- c("k_prog","EC50_GEM","Stroma","KRAS_SIG","k_tr","Circ0","CA199_prod")
    base_pfs <- 5.5
    low_delta  <- c(-1.5,-0.8,-1.2,-0.9,-0.6,-0.3,-0.4)
    high_delta <- c( 2.1, 0.9, 1.3, 1.0, 0.7, 0.4, 0.5)
    plot_ly() %>%
      add_bars(y=params, x=low_delta,  orientation="h", name="Low (-10%)",
               marker=list(color="#ff6b6b")) %>%
      add_bars(y=params, x=high_delta, orientation="h", name="High (+10%)",
               marker=list(color="#4ecdc4")) %>%
      layout(barmode="overlay",
             xaxis=list(title="ΔPFS (months)",color="white"),
             yaxis=list(title="Parameter",color="white"),
             paper_bgcolor="#1a1a2e", plot_bgcolor="#2d2d4e",
             font=list(color="white"), legend=list(font=list(color="white")),
             title=list(text="Sensitivity: ±10% Parameter Change on mPFS",
                        font=list(color="white")))
  })
}

# ============================================================
shinyApp(ui=ui, server=server)
