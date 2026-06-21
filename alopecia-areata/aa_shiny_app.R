## ============================================================
## Alopecia Areata QSP — Interactive Shiny Dashboard
## 6 Tabs: Patient Profile · Drug PK · PD Biomarkers ·
##         Clinical Endpoints · Scenario Comparison · Risk & Biomarker
## ============================================================
## Dependencies: shiny, shinydashboard, plotly, DT, deSolve, dplyr,
##               tidyr, ggplot2
## Run: shiny::runApp("aa_shiny_app.R")
## ============================================================

library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)

## ============================================================
## ODE System (deSolve)
## ============================================================

aa_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    # Active plasma concentration (ng/mL)
    Cp_bari <- if (DRUG_ARM %in% c(2,3)) (ACENT / Vc_bari) * 1000 else 0
    Cp_ritl <- if (DRUG_ARM == 4)        (ACENT / Vc_ritl) * 1000 else 0
    Cp_tofa <- if (DRUG_ARM == 5)        (ACENT / Vc_tofa) * 1000 else 0

    # JAK inhibition
    inhib_j1_bari <- if(DRUG_ARM %in% c(2,3)) Emax_bari*Cp_bari/(IC50_j1_bari+Cp_bari) else 0
    inhib_j2_bari <- if(DRUG_ARM %in% c(2,3)) Emax_bari*Cp_bari/(IC50_j2_bari+Cp_bari) else 0
    inhib_j3_ritl <- if(DRUG_ARM==4) Emax_ritl * JAK3B else 0
    inhib_j1_tofa <- if(DRUG_ARM==5) Emax_tofa*Cp_tofa/(IC50_j1_tofa+Cp_tofa) else 0
    inhib_j3_tofa <- if(DRUG_ARM==5) Emax_tofa*Cp_tofa/(IC50_j3_tofa+Cp_tofa) else 0

    stat1_inhib <- min(inhib_j1_bari + inhib_j1_tofa, 0.95)
    stat5_inhib <- min(max(inhib_j2_bari, inhib_j3_ritl, inhib_j3_tofa), 0.95)
    jak_eff     <- max(1 - (stat1_inhib + stat5_inhib)/2, 0.05)

    # PK ODEs
    if (DRUG_ARM == 2) {
      ka_a=ka_bari; CL_a=CL_bari; Vc_a=Vc_bari; CLd_a=CLd_bari; Vp_a=Vp_bari
    } else if (DRUG_ARM == 3) {
      ka_a=ka_bari; CL_a=CL_bari; Vc_a=Vc_bari; CLd_a=CLd_bari; Vp_a=Vp_bari
    } else if (DRUG_ARM == 4) {
      ka_a=ka_ritl; CL_a=CL_ritl; Vc_a=Vc_ritl; CLd_a=0; Vp_a=0
    } else if (DRUG_ARM == 5) {
      ka_a=ka_tofa; CL_a=CL_tofa; Vc_a=Vc_tofa; CLd_a=0; Vp_a=0
    } else {
      ka_a=0; CL_a=1; Vc_a=1; CLd_a=0; Vp_a=0
    }

    dAGUT  <- -ka_a * AGUT
    dACENT <- ka_a*AGUT - (CL_a/Vc_a)*ACENT - (CLd_a/Vc_a)*ACENT + (CLd_a/Vp_a)*APERI
    dAPERI <- if(Vp_a>0) (CLd_a/Vc_a)*ACENT - (CLd_a/Vp_a)*APERI else 0

    # Covalent JAK3 binding (ritlecitinib)
    Cp_r2 <- if(DRUG_ARM==4) (ACENT/Vc_ritl)*1000 else 0
    dJAK3B <- kbind_jak3*Cp_r2*(1-JAK3B) - ksyn_jak3*JAK3B

    # Immune cells
    dNKG2DL <- kprod_il15*(1.5 + 0.5*IFNG) - kdeg_il15*NKG2DL
    nkc_stim <- IL15C * NKG2DL
    dNKC <- kprol_nk*nkc_stim - kdeath_nk*NKC - stat5_inhib*0.15*NKC
    priming <- kprol_cd8 * IFNG * NKG2DL * (1 - IPIDX)
    dCD8N  <- 0.05 - 0.02*CD8N - priming*CD8N
    dCD8E  <- priming*CD8N*(1 + 0.5*IL15C) - kdeath_cd8e*CD8E - TREG*0.1*CD8E - stat5_inhib*0.20*CD8E
    dTREG  <- kprol_treg*(1 - stat1_inhib*0.3) - kdeath_treg*TREG + 0.01*(treg_ss - TREG)

    # Cytokines
    dIFNG   <- kprod_ifng*(CD8E + 0.4*NKC)*jak_eff + 0.05*DISEASE_ON - kdeg_ifng*IFNG - stat1_inhib*0.4*IFNG
    dIL15C  <- kprod_il15*(1 + 0.3*IFNG)*DISEASE_ON - kdeg_il15*IL15C - stat1_inhib*0.2*IL15C
    dCXCL10 <- kprod_cx10*IFNG - kdeg_cx10*CXCL10

    # JAK/STAT
    dPSTAT1 <- kpSTAT1*IFNG*(1-stat1_inhib) - kdepSTAT1*PSTAT1
    dPSTAT5 <- kpSTAT5*IL15C*(1-stat5_inhib) - kdepSTAT5*PSTAT5

    # Hair follicle
    dIPIDX  <- kip_restore*(1 - IPIDX) - kip_decay*IFNG*IPIDX
    catagen_r <- kcatagen + kcat_ifng*IFNG*(1-IPIDX)
    dANAGEN <- kanagen*(1-ANAGEN) - catagen_r*ANAGEN
    dHAIRDEN<- kdensity*ANAGEN*HAIRDEN*(1-HAIRDEN/100) - kloss*INFLAM*(100-HAIRDEN)/100 +
               0.1*(100-HAIRDEN)*(ANAGEN-0.3)
    dHAIRDEN<- if(HAIRDEN<=0 & dHAIRDEN<0) 0 else if(HAIRDEN>=100 & dHAIRDEN>0) 0 else dHAIRDEN

    # Inflammation & SALT
    dINFLAM <- 0.3*CD8E*IFNG - 0.2*TREG - 0.15*INFLAM - stat1_inhib*0.3*INFLAM
    dSALT   <- kSALT_inc*(100 - HAIRDEN - SALT)
    dDUPIL  <- -CL_dup/Vc_dup * DUPIL

    list(c(dAGUT, dACENT, dAPERI, dCD8N, dCD8E, dTREG, dNKC,
           dIFNG, dIL15C, dCXCL10, dPSTAT1, dPSTAT5,
           dIPIDX, dANAGEN, dHAIRDEN, dSALT, dINFLAM, dNKG2DL, dJAK3B, dDUPIL),
         Cp_bari  = Cp_bari,
         Cp_ritl  = Cp_ritl,
         Cp_tofa  = Cp_tofa,
         stat1_inh_pct = stat1_inhib*100,
         stat5_inh_pct = stat5_inhib*100,
         SALT50 = as.integer(SALT <= 25),
         SALT90 = as.integer(SALT <= 5))
  })
}

## Base parameters
base_params <- list(
  # PK
  ka_bari=1.35, F_bari=0.79, Vc_bari=19.3, Vp_bari=29.4, CLd_bari=3.1, CL_bari=6.2,
  ka_ritl=2.0,  F_ritl=0.70, Vc_ritl=110.0, CL_ritl=66.0,
  ka_tofa=2.5,  F_tofa=0.74, Vc_tofa=29.0, CL_tofa=30.0,
  CL_dup=0.0052, Vc_dup=4.8,
  kbind_jak3=0.08, ksyn_jak3=0.04,
  # PD
  IC50_j1_bari=2.2, IC50_j2_bari=2.1, Emax_bari=0.90,
  IC50_j1_tofa=1.8, IC50_j3_tofa=1.2, Emax_tofa=0.92,
  Emax_ritl=0.85,
  # Immune
  kprol_cd8=0.15, kdeath_cd8e=0.08, kprol_treg=0.02, kdeath_treg=0.04, treg_ss=0.5,
  kprol_nk=0.12, kdeath_nk=0.10,
  kprod_ifng=0.50, kdeg_ifng=0.20, kprod_il15=0.08, kdeg_il15=0.15,
  kprod_cx10=0.30, kdeg_cx10=0.25,
  kpSTAT1=0.60, kdepSTAT1=0.40, kpSTAT5=0.50, kdepSTAT5=0.35,
  # Hair follicle
  kip_decay=0.005, kip_restore=0.002, kanagen=0.010, kcatagen=0.003, kcat_ifng=0.020,
  kdensity=0.003, kloss=0.008,
  # Clinical
  kSALT_inc=0.15,
  DISEASE_ON=1, DRUG_ARM=1
)

## Dosing event generator
make_dose_events <- function(drug_arm, dose_mg, sim_h) {
  if (drug_arm == 1) return(NULL)
  if (drug_arm == 5) {
    # BID: every 12h
    data.frame(var="AGUT", time=sort(c(seq(0,sim_h,24), seq(12,sim_h,24))),
               value=dose_mg/2, method="add")
  } else {
    data.frame(var="AGUT", time=seq(0, sim_h, 24), value=dose_mg, method="add")
  }
}

## Run simulation with dosing
run_sim_shiny <- function(drug_arm, dose_mg, disease_on, sim_weeks,
                           salt_baseline, hairden_baseline, parms_override = list()) {
  p <- modifyList(base_params, c(list(DRUG_ARM=drug_arm, DISEASE_ON=disease_on),
                                  parms_override))

  initial_state <- c(
    AGUT=0, ACENT=0, APERI=0,
    CD8N=1.0, CD8E=0.1, TREG=0.5, NKC=0.8,
    IFNG=1.2, IL15C=1.1, CXCL10=1.5,
    PSTAT1=0.35, PSTAT5=0.25,
    IPIDX=0.30, ANAGEN=0.40, HAIRDEN=hairden_baseline,
    SALT=salt_baseline, INFLAM=0.60, NKG2DL=2.0, JAK3B=0.0, DUPIL=0.0
  )

  sim_h   <- sim_weeks * 168
  times   <- seq(0, sim_h, by=1)
  events  <- make_dose_events(drug_arm, dose_mg, sim_h)

  if (!is.null(events)) {
    ev_list <- lapply(1:nrow(events), function(i) {
      list(var=events$var[i], time=events$time[i], value=events$value[i], method=events$method[i])
    })
    out <- tryCatch(
      ode(y=initial_state, times=times, func=aa_ode, parms=p,
          events=list(data=events), method="lsoda"),
      error=function(e) {
        ode(y=initial_state, times=times, func=aa_ode, parms=p, method="euler", hini=0.1)
      }
    )
  } else {
    out <- ode(y=initial_state, times=times, func=aa_ode, parms=p, method="lsoda")
  }

  df <- as.data.frame(out)
  df$WEEKS <- df$time / 168
  df
}

## ============================================================
## UI
## ============================================================

ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "Alopecia Areata QSP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("① Patient Profile",    tabName="patient",   icon=icon("user")),
      menuItem("② Drug PK",            tabName="pk",        icon=icon("pills")),
      menuItem("③ PD Biomarkers",      tabName="pd",        icon=icon("dna")),
      menuItem("④ Clinical Endpoints", tabName="clinical",  icon=icon("chart-line")),
      menuItem("⑤ Scenario Comparison",tabName="scenario",  icon=icon("layer-group")),
      menuItem("⑥ Risk & Biomarker",   tabName="risk",      icon=icon("microscope"))
    ),
    hr(),
    h5("Global Settings", style="color:#ccc; padding-left:15px"),
    selectInput("drug_arm", "Treatment",
                choices=c("Placebo"=1, "Baricitinib 4 mg QD"=2,
                          "Baricitinib 2 mg QD"=3, "Ritlecitinib 50 mg QD"=4,
                          "Tofacitinib 5 mg BID"=5),
                selected=2),
    sliderInput("sim_weeks", "Simulation (weeks)", 12, 52, 36, step=4),
    sliderInput("salt_init", "Baseline SALT Score", 10, 100, 50, step=5),
    hr(),
    h5("Disease Severity", style="color:#ccc; padding-left:15px"),
    sliderInput("ifng_mult", "IFN-γ Drive (×)", 0.5, 3.0, 1.2, step=0.1),
    sliderInput("treg_ss",   "Treg Set-Point",  0.1, 1.0, 0.5, step=0.05)
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f5f5f5; }
      .box { border-top-color: #7B1FA2; }
    "))),

    tabItems(

      ## ---- Tab 1: Patient Profile ----
      tabItem(tabName="patient",
        fluidRow(
          box(title="Patient Demographics", width=4, status="purple", solidHeader=TRUE,
              numericInput("age", "Age (years)", 35, 10, 80),
              selectInput("sex", "Sex", c("Female"="F","Male"="M")),
              numericInput("bmi", "BMI (kg/m²)", 24, 15, 45),
              selectInput("aa_type",  "AA Subtype",
                          c("Patchy"="patchy","Alopecia Totalis"="at","Alopecia Universalis"="au")),
              selectInput("aa_duration", "Duration",
                          c("<1 year"="short","1–5 years"="mid",">5 years"="long")),
              checkboxGroupInput("comorbid", "Comorbidities",
                                 c("Atopic dermatitis","Thyroid disease","Vitiligo","Anxiety/Depression"))
          ),
          box(title="Disease Severity Overview", width=4, status="purple", solidHeader=TRUE,
              valueBoxOutput("vbox_salt",    width=12),
              valueBoxOutput("vbox_hairden", width=12),
              valueBoxOutput("vbox_ip",      width=12),
              valueBoxOutput("vbox_dx",      width=12)
          ),
          box(title="Natural History Trajectory (Untreated)", width=4, status="purple", solidHeader=TRUE,
              plotlyOutput("p_natural_hist", height=280)
          )
        ),
        fluidRow(
          box(title="AA Classification & Disease Characteristics", width=12, status="purple",
              DTOutput("tbl_classification"))
        )
      ),

      ## ---- Tab 2: Drug PK ----
      tabItem(tabName="pk",
        fluidRow(
          box(title="PK Settings", width=3, status="olive", solidHeader=TRUE,
              sliderInput("dose_mg", "Daily Dose (mg)", 1, 200, 4, step=1),
              selectInput("pk_drug2", "Drug for PK Plot",
                          c("Baricitinib"="bari","Ritlecitinib"="ritl","Tofacitinib"="tofa")),
              numericInput("pk_weeks", "Display (weeks)", 4, 1, 12),
              h5("PK Parameters"),
              tableOutput("tbl_pk_params")
          ),
          box(title="Plasma Concentration vs. Time", width=9, status="olive", solidHeader=TRUE,
              plotlyOutput("p_pk_plot", height=340),
              hr(),
              fluidRow(
                valueBoxOutput("vbox_cmax", width=4),
                valueBoxOutput("vbox_cmin", width=4),
                valueBoxOutput("vbox_auc",  width=4)
              )
          )
        ),
        fluidRow(
          box(title="PK Parameter Reference Table", width=12, status="olive",
              DTOutput("tbl_pk_full"))
        )
      ),

      ## ---- Tab 3: PD Biomarkers ----
      tabItem(tabName="pd",
        fluidRow(
          box(title="p-STAT1 Inhibition (%)", width=6, status="maroon", solidHeader=TRUE,
              plotlyOutput("p_stat1", height=250)),
          box(title="p-STAT5 Inhibition (%)", width=6, status="maroon", solidHeader=TRUE,
              plotlyOutput("p_stat5", height=250))
        ),
        fluidRow(
          box(title="IFN-γ Level (relative)", width=4, status="maroon", solidHeader=TRUE,
              plotlyOutput("p_ifng", height=220)),
          box(title="Serum CXCL10/IP-10 (IFN Signature)", width=4, status="maroon", solidHeader=TRUE,
              plotlyOutput("p_cx10", height=220)),
          box(title="NKG2D Ligand (MICA/MICB)", width=4, status="maroon", solidHeader=TRUE,
              plotlyOutput("p_nkg2dl", height=220))
        ),
        fluidRow(
          box(title="Perifollicular CD8+ T Cells", width=6, status="maroon", solidHeader=TRUE,
              plotlyOutput("p_cd8e", height=220)),
          box(title="Immune Privilege Index", width=6, status="maroon", solidHeader=TRUE,
              plotlyOutput("p_ipidx", height=220))
        )
      ),

      ## ---- Tab 4: Clinical Endpoints ----
      tabItem(tabName="clinical",
        fluidRow(
          box(title="SALT Score Over Time", width=8, status="teal", solidHeader=TRUE,
              plotlyOutput("p_salt_main", height=320),
              hr(),
              fluidRow(
                valueBoxOutput("vbox_salt_w24", width=4),
                valueBoxOutput("vbox_salt50",   width=4),
                valueBoxOutput("vbox_salt90",   width=4)
              )
          ),
          box(title="Hair Density Recovery (%)", width=4, status="teal", solidHeader=TRUE,
              plotlyOutput("p_hairden", height=200),
              hr(),
              plotlyOutput("p_anagen", height=180))
        ),
        fluidRow(
          box(title="OLSS Score (Eyebrow/Lash Regrowth Proxy)", width=6, status="teal", solidHeader=TRUE,
              plotlyOutput("p_olss", height=220)),
          box(title="Inflammation Index Over Time", width=6, status="teal", solidHeader=TRUE,
              plotlyOutput("p_inflam", height=220))
        )
      ),

      ## ---- Tab 5: Scenario Comparison ----
      tabItem(tabName="scenario",
        fluidRow(
          box(title="All 5 Scenarios — SALT Score", width=12, status="blue", solidHeader=TRUE,
              plotlyOutput("p_scenario_salt", height=350))
        ),
        fluidRow(
          box(title="IFN-γ Comparison", width=6, status="blue", solidHeader=TRUE,
              plotlyOutput("p_scenario_ifng", height=250)),
          box(title="Hair Density Comparison", width=6, status="blue", solidHeader=TRUE,
              plotlyOutput("p_scenario_hairden", height=250))
        ),
        fluidRow(
          box(title="Summary Table at Selected Week", width=12, status="blue",
              sliderInput("summary_week", "Timepoint (weeks)", 4, 52, 36, step=4),
              DTOutput("tbl_scenario_summary"))
        )
      ),

      ## ---- Tab 6: Risk & Biomarker ----
      tabItem(tabName="risk",
        fluidRow(
          box(title="Dose–Response: Baricitinib SALT Score", width=6, status="purple", solidHeader=TRUE,
              plotlyOutput("p_dose_resp", height=280)),
          box(title="Biomarker Correlations (CXCL10 vs SALT)", width=6, status="purple", solidHeader=TRUE,
              plotlyOutput("p_biom_corr", height=280))
        ),
        fluidRow(
          box(title="AA Risk Stratification", width=4, status="purple", solidHeader=TRUE,
              h4("Prognosis Calculator"),
              numericInput("salt_prog", "Current SALT", 50, 0, 100),
              numericInput("dur_prog",  "Duration (yrs)", 2, 0, 30),
              selectInput("subtype_prog", "Subtype",
                          c("Patchy"="patchy","Alopecia Totalis"="at","Alopecia Universalis"="au")),
              actionButton("calc_prog", "Calculate Risk", class="btn-warning"),
              hr(),
              htmlOutput("risk_result")
          ),
          box(title="Biomarker Reference Ranges", width=4, status="purple", solidHeader=TRUE,
              DTOutput("tbl_biomarkers")),
          box(title="Weekly Biomarker Summary (downloadable)", width=4, status="purple", solidHeader=TRUE,
              DTOutput("tbl_weekly_biom"),
              downloadButton("dl_sim", "Download CSV"))
        )
      )
    )
  )
)

## ============================================================
## Server
## ============================================================

server <- function(input, output, session) {

  ## ---- Reactive: single-arm simulation ----
  sim_data <- reactive({
    parms_ov <- list(
      kprod_ifng = base_params$kprod_ifng * input$ifng_mult,
      treg_ss    = input$treg_ss
    )
    run_sim_shiny(
      drug_arm       = as.integer(input$drug_arm),
      dose_mg        = input$dose_mg,
      disease_on     = 1,
      sim_weeks      = input$sim_weeks,
      salt_baseline  = input$salt_init,
      hairden_baseline= 100 - input$salt_init,
      parms_override = parms_ov
    )
  })

  ## ---- Reactive: all 5 scenarios ----
  sim_all_arms <- reactive({
    arms_cfg <- list(
      list(arm=1, label="Placebo",              dose=0),
      list(arm=2, label="Baricitinib 4 mg QD",  dose=4),
      list(arm=3, label="Baricitinib 2 mg QD",  dose=2),
      list(arm=4, label="Ritlecitinib 50 mg QD",dose=50),
      list(arm=5, label="Tofacitinib 5 mg BID", dose=5)
    )
    parms_ov <- list(kprod_ifng=base_params$kprod_ifng*input$ifng_mult, treg_ss=input$treg_ss)
    bind_rows(lapply(arms_cfg, function(a) {
      df <- run_sim_shiny(a$arm, a$dose, 1, input$sim_weeks,
                          input$salt_init, 100-input$salt_init, parms_ov)
      df$ARM <- a$label
      df
    }))
  })

  arm_colors <- c("Placebo"="#9E9E9E","Baricitinib 4 mg QD"="#1565C0",
                  "Baricitinib 2 mg QD"="#42A5F5","Ritlecitinib 50 mg QD"="#7B1FA2",
                  "Tofacitinib 5 mg BID"="#2E7D32")

  ## ---- Tab 1: Patient Profile ----
  sim_pbo <- reactive({
    run_sim_shiny(1, 0, 1, input$sim_weeks, input$salt_init, 100-input$salt_init)
  })

  output$vbox_salt <- renderValueBox({
    last <- tail(sim_pbo(), 1)
    valueBox(sprintf("%.0f", last$SALT), "SALT Score (untreated)",
             icon=icon("percentage"), color="red")
  })
  output$vbox_hairden <- renderValueBox({
    last <- tail(sim_pbo(), 1)
    valueBox(sprintf("%.0f%%", last$HAIRDEN), "Hair Density",
             icon=icon("user"), color="purple")
  })
  output$vbox_ip <- renderValueBox({
    last <- tail(sim_pbo(), 1)
    valueBox(sprintf("%.2f", last$IPIDX), "Immune Privilege Index",
             icon=icon("shield-alt"), color="orange")
  })
  output$vbox_dx <- renderValueBox({
    subtype <- switch(input$aa_type, patchy="Patchy AA", at="Alopecia Totalis", au="Alopecia Universalis")
    valueBox(subtype, "Subtype", icon=icon("diagnoses"), color="blue")
  })

  output$p_natural_hist <- renderPlotly({
    df <- sim_pbo()
    plot_ly(df, x=~WEEKS, y=~SALT, type='scatter', mode='lines',
            line=list(color='#E53935', width=2), name="SALT") %>%
      add_lines(y=~HAIRDEN, line=list(color='#43A047', dash='dash'), name="Hair Density") %>%
      layout(title="Untreated Disease Course",
             xaxis=list(title="Weeks"), yaxis=list(title="Score / %"))
  })

  output$tbl_classification <- renderDT({
    tibble(
      "Feature"     = c("SALT Range","Classification","Prognosis","Typical onset","Nail involvement","Comorbidity risk"),
      "Patchy AA"   = c("1–99%","S1-S5","Variable, may resolve","Any age","10–20%","Atopic dermatitis 25%"),
      "Alopecia Totalis" = c("100% scalp","AT","Chronic in 50%","Often childhood","15–30%","Thyroid 14–30%"),
      "Alopecia Universalis" = c("All body hair","AU","Usually chronic","Childhood/young adult","20–50%","Multiple autoimmune")
    ) %>% datatable(options=list(dom='t', pageLength=6), rownames=FALSE)
  })

  ## ---- Tab 2: Drug PK ----
  output$tbl_pk_params <- renderTable({
    tibble(
      Parameter = c("Baricitinib","Ritlecitinib","Tofacitinib"),
      "Dose (mg)" = c("2 or 4 QD","50 QD","5 BID"),
      "F (%)" = c("79","70","74"),
      "t½ (h)" = c("12","~2","3"),
      "Vc (L)" = c("19.3","110","29"),
      "CL (L/h)" = c("6.2","66","30")
    )
  })

  output$p_pk_plot <- renderPlotly({
    df <- sim_data() %>% filter(WEEKS <= input$pk_weeks)
    pk_col <- switch(input$pk_drug2, bari="Cp_bari", ritl="Cp_ritl", tofa="Cp_tofa")
    y_vals <- df[[pk_col]]
    plot_ly(df, x=~time, y=y_vals, type='scatter', mode='lines',
            line=list(color='#00695C', width=2)) %>%
      layout(title=paste(tools::toTitleCase(input$pk_drug2), "Plasma Concentration"),
             xaxis=list(title="Time (hours)"), yaxis=list(title="Cp (ng/mL)"))
  })

  output$vbox_cmax <- renderValueBox({
    df <- sim_data(); pk_col <- switch(input$pk_drug2, bari="Cp_bari", ritl="Cp_ritl", tofa="Cp_tofa")
    valueBox(sprintf("%.1f ng/mL", max(df[[pk_col]], na.rm=TRUE)), "Cmax", icon=icon("arrow-up"), color="green")
  })
  output$vbox_cmin <- renderValueBox({
    df <- sim_data(); pk_col <- switch(input$pk_drug2, bari="Cp_bari", ritl="Cp_ritl", tofa="Cp_tofa")
    ss_vals <- df[[pk_col]][df$WEEKS > (input$sim_weeks-1)]
    valueBox(sprintf("%.1f ng/mL", min(ss_vals, na.rm=TRUE)), "Cmin (SS)", icon=icon("arrow-down"), color="orange")
  })
  output$vbox_auc <- renderValueBox({
    df <- sim_data(); pk_col <- switch(input$pk_drug2, bari="Cp_bari", ritl="Cp_ritl", tofa="Cp_tofa")
    auc_val <- sum(diff(df$time) * (head(df[[pk_col]],-1) + tail(df[[pk_col]],-1))/2, na.rm=TRUE)
    valueBox(sprintf("%.0f ng·h/mL", auc_val), "AUCtotal", icon=icon("chart-area"), color="blue")
  })

  output$tbl_pk_full <- renderDT({
    tibble(
      Drug=c("Baricitinib","Ritlecitinib","Tofacitinib","Dupilumab"),
      Mechanism=c("JAK1/2 reversible","JAK3/TEC covalent","JAK1/3 reversible","IL-4Rα mAb"),
      "Route/Freq"=c("PO QD","PO QD","PO BID","SC Q2W"),
      "Approved dose (AA)"=c("4 mg","50 mg","Off-label 5 mg","Off-label"),
      "F (%)"=c("79","70","74","64"),
      "Tmax (h)"=c("0.7–1.0","1.0–2.0","0.5–1.0","168"),
      "t½ (h)"=c("12","~2","3","528"),
      "CL (L/h)"=c("6.2","66","30","0.005"),
      "Vc (L)"=c("19.3","110","29","4.8"),
      "Primary elim"=c("Renal 75%","CYP3A4","CYP3A4/2C19","Receptor-mediated"),
      "Key trial"=c("BRAVE-AA1/2","ALLEGRO","Liu 2019","Experimental")
    ) %>% datatable(options=list(scrollX=TRUE, pageLength=4), rownames=FALSE)
  })

  ## ---- Tab 3: PD Biomarkers ----
  output$p_stat1 <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~WEEKS, y=~stat1_inh_pct, type='scatter', mode='lines',
            line=list(color='#7B1FA2', width=2)) %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="p-STAT1 Inhibition (%)"))
  })
  output$p_stat5 <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~WEEKS, y=~stat5_inh_pct, type='scatter', mode='lines',
            line=list(color='#1565C0', width=2)) %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="p-STAT5 Inhibition (%)"))
  })
  output$p_ifng <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~WEEKS, y=~IFNG, type='scatter', mode='lines',
            line=list(color='#E65100', width=2)) %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="IFN-γ (rel.)"))
  })
  output$p_cx10 <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~WEEKS, y=~CXCL10, type='scatter', mode='lines',
            line=list(color='#FF6F00', width=2)) %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="CXCL10 (rel.)"))
  })
  output$p_nkg2dl <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~WEEKS, y=~NKG2DL, type='scatter', mode='lines',
            line=list(color='#1A237E', width=2)) %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="NKG2D Ligand (rel.)"))
  })
  output$p_cd8e <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~WEEKS, y=~CD8E, type='scatter', mode='lines',
            line=list(color='#C62828', width=2)) %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="CD8+ Effector (rel.)"))
  })
  output$p_ipidx <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~WEEKS, y=~IPIDX, type='scatter', mode='lines',
            line=list(color='#2E7D32', width=2)) %>%
      add_lines(y=rep(1, nrow(df)), line=list(dash='dash', color='gray'), name="Healthy") %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="IP Index (0–1)"))
  })

  ## ---- Tab 4: Clinical Endpoints ----
  output$p_salt_main <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~WEEKS, y=~SALT, type='scatter', mode='lines',
            line=list(color='#7B1FA2', width=2.5), name="SALT") %>%
      add_lines(y=rep(25, nrow(df)), line=list(dash='dash', color='#1565C0'), name="SALT50 threshold") %>%
      add_lines(y=rep(5,  nrow(df)), line=list(dash='dot',  color='#E65100'), name="SALT90 threshold") %>%
      layout(title="SALT Score Trajectory",
             xaxis=list(title="Weeks"), yaxis=list(title="SALT Score (0=none, 100=complete)"))
  })
  output$vbox_salt_w24 <- renderValueBox({
    df <- sim_data(); w24 <- df[df$WEEKS >= 23.9 & df$WEEKS <= 24.1,]
    if(nrow(w24)>0) val <- sprintf("%.0f", mean(w24$SALT)) else val <- "N/A"
    valueBox(val, "SALT at Week 24", icon=icon("calendar"), color="purple")
  })
  output$vbox_salt50 <- renderValueBox({
    df <- sim_data(); last <- tail(df,1)
    valueBox(if(last$SALT50>0.5)"RESPONDER" else "NON-RESPONDER",
             "SALT50 Status", icon=icon("check"), color=if(last$SALT50>0.5)"green" else "red")
  })
  output$vbox_salt90 <- renderValueBox({
    df <- sim_data(); last <- tail(df,1)
    valueBox(if(last$SALT90>0.5)"RESPONDER" else "NON-RESPONDER",
             "SALT90 Status", icon=icon("star"), color=if(last$SALT90>0.5)"green" else "orange")
  })
  output$p_hairden <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~WEEKS, y=~HAIRDEN, type='scatter', mode='lines',
            line=list(color='#43A047', width=2)) %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="Hair Density (%)", range=c(0,100)))
  })
  output$p_anagen <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~WEEKS, y=~ANAGEN, type='scatter', mode='lines',
            line=list(color='#66BB6A', width=2)) %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="Anagen Fraction", range=c(0,1)))
  })
  output$p_olss <- renderPlotly({
    df <- sim_data()
    olss_val <- ifelse(df$HAIRDEN > 80, 3, ifelse(df$HAIRDEN > 50, 2, 1))
    plot_ly(df, x=~WEEKS, y=olss_val, type='scatter', mode='lines',
            line=list(color='#FF8F00', width=2)) %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="OLSS Score (1=absent, 3=complete)",
                                                    range=c(0,4)))
  })
  output$p_inflam <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~WEEKS, y=~INFLAM, type='scatter', mode='lines',
            line=list(color='#E53935', width=2)) %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="Inflammation Index (0–1)"))
  })

  ## ---- Tab 5: Scenario Comparison ----
  output$p_scenario_salt <- renderPlotly({
    df <- sim_all_arms()
    p <- plot_ly()
    for(arm in unique(df$ARM)) {
      sub <- df[df$ARM==arm,]
      p <- add_lines(p, x=sub$WEEKS, y=sub$SALT, name=arm,
                     line=list(color=arm_colors[arm], width=2.5))
    }
    p %>% add_lines(x=df$WEEKS[df$ARM=="Placebo"], y=rep(25, sum(df$ARM=="Placebo")),
                    line=list(dash='dash',color='gray'), name="SALT50 threshold") %>%
      layout(title="All Scenarios — SALT Score (36 Weeks)",
             xaxis=list(title="Weeks"), yaxis=list(title="SALT Score"))
  })
  output$p_scenario_ifng <- renderPlotly({
    df <- sim_all_arms()
    p <- plot_ly()
    for(arm in unique(df$ARM))
      p <- add_lines(p, x=df$WEEKS[df$ARM==arm], y=df$IFNG[df$ARM==arm],
                     name=arm, line=list(color=arm_colors[arm], width=2))
    p %>% layout(xaxis=list(title="Weeks"), yaxis=list(title="IFN-γ (rel.)"))
  })
  output$p_scenario_hairden <- renderPlotly({
    df <- sim_all_arms()
    p <- plot_ly()
    for(arm in unique(df$ARM))
      p <- add_lines(p, x=df$WEEKS[df$ARM==arm], y=df$HAIRDEN[df$ARM==arm],
                     name=arm, line=list(color=arm_colors[arm], width=2))
    p %>% layout(xaxis=list(title="Weeks"), yaxis=list(title="Hair Density (%)"))
  })
  output$tbl_scenario_summary <- renderDT({
    df <- sim_all_arms()
    wk <- input$summary_week
    tbl <- df %>%
      filter(abs(WEEKS - wk) < 0.6) %>%
      group_by(ARM) %>%
      slice_tail(n=1) %>%
      select(ARM, SALT, HAIRDEN, IFNG, CXCL10, PSTAT1, stat1_inh_pct, SALT50, SALT90) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
    datatable(tbl, options=list(dom='t'), rownames=FALSE)
  })

  ## ---- Tab 6: Risk & Biomarker ----
  output$p_dose_resp <- renderPlotly({
    doses <- c(1, 2, 4, 6, 8)
    parms_ov <- list(kprod_ifng=base_params$kprod_ifng*input$ifng_mult, treg_ss=input$treg_ss)
    results <- lapply(doses, function(d) {
      df <- run_sim_shiny(2, d, 1, input$sim_weeks, input$salt_init, 100-input$salt_init, parms_ov)
      data.frame(dose=d, SALT=tail(df$SALT,1), HAIRDEN=tail(df$HAIRDEN,1))
    })
    dr <- bind_rows(results)
    plot_ly(dr, x=~dose, y=~SALT, type='scatter', mode='lines+markers',
            line=list(color='#1565C0', width=2), marker=list(size=8)) %>%
      layout(title="Baricitinib Dose–Response (SALT at study end)",
             xaxis=list(title="Daily Dose (mg)"), yaxis=list(title="SALT Score"))
  })

  output$p_biom_corr <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x=~CXCL10, y=~SALT, type='scatter', mode='markers',
            marker=list(color=~WEEKS, colorscale='Viridis', size=4, showscale=TRUE),
            text=~paste("Week:", round(WEEKS,1))) %>%
      layout(title="CXCL10 vs. SALT Score",
             xaxis=list(title="CXCL10 (rel.)"), yaxis=list(title="SALT Score"))
  })

  observeEvent(input$calc_prog, {
    salt  <- input$salt_prog
    dur   <- input$dur_prog
    subtp <- input$subtype_prog
    risk_score <- salt/100 * 0.4 + min(dur/10, 1) * 0.3 +
                  ifelse(subtp=="au", 0.3, ifelse(subtp=="at", 0.2, 0.1))
    risk_cat <- if(risk_score > 0.6) "HIGH" else if(risk_score > 0.35) "MODERATE" else "LOW"
    risk_col <- if(risk_score > 0.6) "red" else if(risk_score > 0.35) "orange" else "green"
    output$risk_result <- renderUI({
      tagList(
        h4(sprintf("Risk Category: %s", risk_cat), style=paste0("color:", risk_col)),
        p(sprintf("Risk Score: %.2f", risk_score)),
        p(sprintf("Spontaneous remission likelihood: %s",
                  if(risk_score<0.35) "30-50%" else if(risk_score<0.6) "10-25%" else "<10%")),
        p(sprintf("JAK inhibitor recommendation: %s",
                  if(risk_score>0.35) "INDICATED (SALT ≥20% + significant duration)" else "Consider watchful waiting"))
      )
    })
  })

  output$tbl_biomarkers <- renderDT({
    tibble(
      Biomarker=c("Serum CXCL10/IP-10","CD8+ perifollicular T cells",
                  "p-STAT1 (PBMC)","p-STAT5 (PBMC)","Serum IFN-γ",
                  "NKG2D ligands (MICA)","Foxp3+ Treg density","IL-15","SALT score"),
      "Normal"=c("<200 pg/mL","<5%","<10%","<15%","<1 pg/mL","Low",">5/HPF","<100 pg/mL","0"),
      "Active AA"=c("400-2000 pg/mL","20-50%","25-60%","20-50%","2-10 pg/mL","High","<2/HPF","200-500 pg/mL","10-100"),
      "On JAK inh"=c("↓50-80%","↓40-70%","↓50-90%","↓40-85%","↓40-70%","↓30-50%","↑","↓20-40%","↓SALT50 36%")
    ) %>% datatable(options=list(dom='t', pageLength=9), rownames=FALSE)
  })

  output$tbl_weekly_biom <- renderDT({
    df <- sim_data() %>%
      filter(WEEKS %in% c(0,4,8,12,16,20,24,36)) %>%
      select(WEEKS, SALT, HAIRDEN, IFNG, CXCL10, PSTAT1, CD8E, IPIDX) %>%
      mutate(across(where(is.numeric), ~round(.x, 3)))
    datatable(df, options=list(pageLength=8, dom='t'), rownames=FALSE)
  })

  output$dl_sim <- downloadHandler(
    filename = function() paste0("aa_sim_", Sys.Date(), ".csv"),
    content  = function(file) {
      write.csv(sim_data(), file, row.names=FALSE)
    }
  )
}

shinyApp(ui=ui, server=server)
