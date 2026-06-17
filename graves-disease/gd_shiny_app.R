################################################################################
# Graves' Disease QSP — Interactive Shiny Dashboard
# 6 Tabs: Patient Profile · PK · Thyroid Hormone PD · Cardiovascular & Bone
#         Immune & Ophthalmopathy · Scenario Comparison
################################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(patchwork)

# ─── Inline mrgsolve model (abbreviated for Shiny) ────────────────────────────
graves_code <- '
$PARAM
ksyn_T4=1.12, ksyn_T3=0.18, kD1_T4T3=0.060, kD3_T4rT3=0.015, kD1_rT3=0.40,
kel_T4=0.099, kel_T3=0.693, kel_rT3=1.386, frac_fT4=0.00025, frac_fT3=0.003,
ksyn_TSH=0.50, kel_TSH=1.386, IC50_fT4=14.0, IC50_fT3=4.5, n_TSH=2.0,
fT4_norm=16.0, fT3_norm=5.5, E_TRAb_max=3.5, EC50_TRAb=8.0,
E_TSH_max=1.0, EC50_TSH=1.5,
ksyn_Bcell=0.05, kdeath_Bcell=0.02, kprod_TRAb=0.20, kel_TRAb=0.021,
kgrowth_thy=0.01, kdeath_thy=0.003, RAI_kill=0.10, NIS_RAI=0.50,
kphy_RAI=0.087, kbio_RAI=0.15, kelb_RAI=0.050,
ka_MMI=8.0, kel_MMI=4.0, F_MMI=0.93, IC50_MMI=0.15,
ka_PTU=5.0, kel_PTU=9.24, IC50_PTU=1.0, IC50_PTU_D1=5.0,
ka_PROP=6.0, kel_PROP=4.16, EC50_PROP_HR=0.02,
kBR_base=1.0, kBR_T3stim=0.25, kBR_kel=0.10,
HR_base=72.0, kHR_T3=0.80, HR_max_stim=40.0, HR_EC50_fT3=8.0,
kGO_TRAb=0.08, kGO_kel=0.015, GCS_GO_inh=0.70, TPO_inh_max=0.98,
MMI_Bcell=0.35, RTX_Bcell=0.95,
Thyroid_mass0=1.0, TRAb0=25.0, Bcell0=10.0,
USE_MMI=0, USE_PTU=0, USE_RAI=0, USE_PROP=0, USE_LT4=0, USE_RTX=0, USE_GCS=0

$CMT MMI PTU_C RAI_SERUM RAI_THYR PROP_C TSH_C T4_C T3_C fT4_C fT3_C rT3_C
     TRAb_C Bcell_C BoneResor HR_dev GO_act ThyMass TPO_inh

$ODE
dxdt_MMI=     -kel_MMI*MMI;
dxdt_PTU_C=   -kel_PTU*PTU_C;
dxdt_RAI_SERUM= -NIS_RAI*RAI_SERUM*(ThyMass/Thyroid_mass0)-kelb_RAI*RAI_SERUM;
dxdt_RAI_THYR=  NIS_RAI*RAI_SERUM*(ThyMass/Thyroid_mass0)-(kphy_RAI+kbio_RAI)*RAI_THYR;
dxdt_PROP_C=  -kel_PROP*PROP_C;
double TPO_inh_MMI=(USE_MMI>0.5)?TPO_inh_max*MMI/(IC50_MMI+MMI):0.0;
double TPO_inh_PTU=(USE_PTU>0.5)?TPO_inh_max*PTU_C/(IC50_PTU+PTU_C):0.0;
double TPO_tot=1.0-(1.0-TPO_inh_MMI)*(1.0-TPO_inh_PTU);
dxdt_TPO_inh=(TPO_tot-TPO_inh)*2.0;
double D1f=1.0-(USE_PTU>0.5?0.5*PTU_C/(IC50_PTU_D1+PTU_C):0.0)-(USE_PROP>0.5?0.3*PROP_C/(EC50_PROP_HR+PROP_C):0.0);
if(D1f<0.05)D1f=0.05;
double TRAb_s=E_TRAb_max*TRAb_C/(EC50_TRAb+TRAb_C);
double TSH_s=E_TSH_max*TSH_C/(EC50_TSH+TSH_C);
double TSHR_a=TSH_s+TRAb_s;
double RAI_dmg=(USE_RAI>0.5)?RAI_kill*RAI_THYR:0.0;
dxdt_ThyMass=kgrowth_thy*ThyMass*(TSHR_a-1.0)-kdeath_thy*ThyMass-RAI_dmg*ThyMass;
if(ThyMass<0.05)dxdt_ThyMass=0.0;
double sf=ThyMass*TSHR_a*(1.0-TPO_inh);
double T4s=ksyn_T4*sf; double T3s=ksyn_T3*sf;
dxdt_T4_C=T4s-kel_T4*T4_C-kD1_T4T3*D1f*fT4_C-kD3_T4rT3*fT4_C;
dxdt_T3_C=T3s+kD1_T4T3*D1f*fT4_C-kel_T3*T3_C;
dxdt_fT4_C=frac_fT4*T4s-kel_T4*fT4_C-kD1_T4T3*D1f*fT4_C-kD3_T4rT3*fT4_C;
dxdt_fT3_C=frac_fT3*T3s+kD1_T4T3*D1f*fT4_C*frac_fT4/frac_fT3*0.15-kel_T3*fT3_C;
dxdt_rT3_C=kD3_T4rT3*fT4_C-kD1_rT3*rT3_C;
double fb=(pow(fT4_C/fT4_norm,n_TSH)+pow(fT3_C/fT3_norm,n_TSH))/2.0;
dxdt_TSH_C=ksyn_TSH/fb-kel_TSH*TSH_C;
double Bg=ksyn_Bcell*Bcell_C; double Bd=kdeath_Bcell*Bcell_C;
double mi=(USE_MMI>0.5)?MMI_Bcell*MMI/(IC50_MMI+MMI):0.0;
double ri=(USE_RTX>0.5)?RTX_Bcell:0.0;
dxdt_Bcell_C=Bg-Bd-mi*Bcell_C-ri*Bcell_C;
dxdt_TRAb_C=kprod_TRAb*Bcell_C-kel_TRAb*TRAb_C;
dxdt_BoneResor=kBR_T3stim*(fT3_C/fT3_norm-1.0)-kBR_kel*(BoneResor-1.0);
double fT3e=HR_max_stim*fT3_C/(HR_EC50_fT3+fT3_C);
double ph=(USE_PROP>0.5)?HR_max_stim*PROP_C/(EC50_PROP_HR+PROP_C):0.0;
dxdt_HR_dev=0.5*(fT3e-ph-HR_dev);
double gcs=(USE_GCS>0.5)?GCS_GO_inh:0.0;
dxdt_GO_act=kGO_TRAb*TRAb_C*(1.0-gcs)-kGO_kel*GO_act;

$TABLE
double TSH_obs=TSH_C; double fT4_obs=fT4_C; double fT3_obs=fT3_C;
double TRAb_obs=TRAb_C; double HR_obs=HR_base+HR_dev;
$CAPTURE TSH_obs fT4_obs fT3_obs TRAb_obs HR_obs ThyMass TRAb_C BoneResor GO_act TPO_inh
'

gd_mod <- mrgsolve::mcode("graves_shiny", graves_code, quiet=TRUE)

# ─── Helper ────────────────────────────────────────────────────────────────────
run_sim <- function(mod, params, inits, duration=730, mmi_dose=0, ptu_dose=0,
                    rai_dose=0, prop_dose=0) {
  evs <- NULL
  if(mmi_dose  > 0) evs <- c(evs, list(ev(amt=mmi_dose/100,  ii=0.33, addl=duration*3, cmt=1, time=0)))
  if(ptu_dose  > 0) evs <- c(evs, list(ev(amt=ptu_dose/300,  ii=0.25, addl=duration*4, cmt=2, time=0)))
  if(rai_dose  > 0) evs <- c(evs, list(ev(amt=rai_dose,      cmt=3, time=0)))
  if(prop_dose > 0) evs <- c(evs, list(ev(amt=prop_dose/640, ii=0.17, addl=duration*6, cmt=5, time=0)))

  m <- do.call(param, c(list(mod), params))
  m <- do.call(init,  c(list(m),   inits))

  if(!is.null(evs)) {
    all_ev <- evs[[1]]
    if(length(evs) > 1) for(e in evs[-1]) all_ev <- all_ev + e
    mrgsim(m, events=all_ev, end=duration, delta=1) %>% as.data.frame()
  } else {
    mrgsim(m, end=duration, delta=1) %>% as.data.frame()
  }
}

# ─── UI ────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Graves' Disease QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",       tabName="tab_patient",  icon=icon("user")),
      menuItem("Drug PK",               tabName="tab_pk",       icon=icon("pills")),
      menuItem("Thyroid Hormone PD",    tabName="tab_hormone",  icon=icon("chart-line")),
      menuItem("Cardiovascular & Bone", tabName="tab_cvbone",   icon=icon("heart")),
      menuItem("Immune & Ophthalmo",    tabName="tab_immune",   icon=icon("eye")),
      menuItem("Scenario Comparison",   tabName="tab_scenario", icon=icon("layer-group"))
    ),

    hr(),
    h5("  Disease Settings", style="color:#aaa; padding-left:15px;"),
    sliderInput("TRAb_init",  "Initial TRAb (IU/L)",  min=2, max=60, value=25),
    sliderInput("fT4_init",   "Initial fT4 (pmol/L)", min=22, max=80, value=35),
    sliderInput("sim_dur",    "Simulation (days)",     min=90, max=1095, value=730, step=30),

    hr(),
    h5("  Drug Therapy",  style="color:#aaa; padding-left:15px;"),
    checkboxInput("use_mmi",  "Methimazole (MMI)",   value=FALSE),
    conditionalPanel("input.use_mmi",
      sliderInput("dose_mmi", "MMI Dose (mg/day)", 5, 60, 30, step=5)),
    checkboxInput("use_ptu",  "PTU", value=FALSE),
    conditionalPanel("input.use_ptu",
      sliderInput("dose_ptu", "PTU Dose (mg/day)", 50, 600, 300, step=50)),
    checkboxInput("use_rai",  "Radioiodine (¹³¹I)",  value=FALSE),
    conditionalPanel("input.use_rai",
      sliderInput("dose_rai", "¹³¹I (mCi)", 5, 30, 15)),
    checkboxInput("use_prop", "Propranolol",  value=FALSE),
    conditionalPanel("input.use_prop",
      sliderInput("dose_prop", "Propranolol (mg/day)", 20, 160, 80, step=20)),
    checkboxInput("use_gcs",  "Glucocorticoids (GO)", value=FALSE),
    checkboxInput("use_rtx",  "Rituximab",  value=FALSE),

    hr(),
    actionButton("run_sim", "Run Simulation", class="btn-primary btn-block",
                 icon=icon("play"))
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background:#f4f6f9; }
      .box { border-top-color: #3c8dbc; }
      .value-box .icon { font-size:36px; }
    "))),

    tabItems(
      # ── TAB 1: Patient Profile ──────────────────────────────────────────────
      tabItem("tab_patient",
        fluidRow(
          valueBoxOutput("vbox_TSH",   width=3),
          valueBoxOutput("vbox_fT4",   width=3),
          valueBoxOutput("vbox_TRAb",  width=3),
          valueBoxOutput("vbox_HR",    width=3)
        ),
        fluidRow(
          box(title="Patient Status Summary", width=12, solidHeader=TRUE,
              status="primary",
              fluidRow(
                column(4,
                  h4("Diagnosis"),
                  uiOutput("diag_status"),
                  hr(),
                  h4("Current Therapy"),
                  uiOutput("therapy_list")
                ),
                column(4,
                  h4("Lab Results at Endpoint"),
                  tableOutput("lab_table")
                ),
                column(4,
                  h4("Disease Severity"),
                  plotOutput("radar_plot", height="260px")
                )
              )
          )
        ),
        fluidRow(
          box(title="Thyroid Hormone Time Course", width=12, solidHeader=TRUE,
              plotlyOutput("overview_plot", height="320px"))
        )
      ),

      # ── TAB 2: Drug PK ──────────────────────────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(title="MMI / PTU Plasma Concentration",  width=6, solidHeader=TRUE,
              status="warning", plotlyOutput("pk_plot_mmi", height="300px")),
          box(title="Propranolol Plasma Concentration", width=6, solidHeader=TRUE,
              status="info",    plotlyOutput("pk_plot_prop", height="300px"))
        ),
        fluidRow(
          box(title="Radioiodine ¹³¹I Kinetics", width=6, solidHeader=TRUE,
              status="danger", plotlyOutput("pk_plot_rai", height="300px")),
          box(title="TPO Inhibition by ATD", width=6, solidHeader=TRUE,
              status="primary", plotlyOutput("pk_plot_tpo", height="300px"))
        ),
        fluidRow(
          box(title="Drug PK Parameters", width=12, solidHeader=TRUE,
              DTOutput("pk_param_table"))
        )
      ),

      # ── TAB 3: Thyroid Hormone PD ────────────────────────────────────────────
      tabItem("tab_hormone",
        fluidRow(
          box(title="Serum TSH",             width=4, solidHeader=TRUE,
              status="primary", plotlyOutput("pd_tsh",  height="260px")),
          box(title="Free T4 (fT4)",          width=4, solidHeader=TRUE,
              status="warning", plotlyOutput("pd_fT4",  height="260px")),
          box(title="Free T3 (fT3)",          width=4, solidHeader=TRUE,
              status="warning", plotlyOutput("pd_fT3",  height="260px"))
        ),
        fluidRow(
          box(title="Total T4 & T3",         width=4, solidHeader=TRUE,
              plotlyOutput("pd_tot_T4T3", height="260px")),
          box(title="Reverse T3 (rT3)",      width=4, solidHeader=TRUE,
              plotlyOutput("pd_rT3", height="260px")),
          box(title="Thyroid Mass (Goiter)",  width=4, solidHeader=TRUE,
              status="danger", plotlyOutput("pd_thyroid_mass", height="260px"))
        )
      ),

      # ── TAB 4: Cardiovascular & Bone ─────────────────────────────────────────
      tabItem("tab_cvbone",
        fluidRow(
          box(title="Heart Rate",       width=6, solidHeader=TRUE,
              status="danger", plotlyOutput("cv_hr",  height="300px")),
          box(title="Bone Resorption Index (Osteoclast Activity)", width=6,
              solidHeader=TRUE, status="warning", plotlyOutput("cv_bone", height="300px"))
        ),
        fluidRow(
          box(title="Cumulative BMD Loss Estimate", width=6, solidHeader=TRUE,
              plotlyOutput("cv_bmd", height="280px")),
          box(title="CV & Bone Risk Summary", width=6, solidHeader=TRUE,
              tableOutput("cv_bone_table"))
        )
      ),

      # ── TAB 5: Immune & Ophthalmopathy ───────────────────────────────────────
      tabItem("tab_immune",
        fluidRow(
          box(title="TRAb (Thyroid-Stimulating Antibody)", width=6,
              solidHeader=TRUE, status="danger", plotlyOutput("imm_trab", height="300px")),
          box(title="B-cell Activity",           width=6, solidHeader=TRUE,
              status="warning", plotlyOutput("imm_bcell", height="300px"))
        ),
        fluidRow(
          box(title="Graves Ophthalmopathy (Orbital Activation Score)", width=6,
              solidHeader=TRUE, status="primary", plotlyOutput("imm_go", height="300px")),
          box(title="Remission Probability Over Time", width=6,
              solidHeader=TRUE, plotlyOutput("imm_remission", height="300px"))
        )
      ),

      # ── TAB 6: Scenario Comparison ───────────────────────────────────────────
      tabItem("tab_scenario",
        fluidRow(
          box(title="Scenario Setup", width=4, solidHeader=TRUE, status="primary",
              checkboxGroupInput("scenarios_sel",
                label="Select Scenarios to Compare:",
                choices=c("Untreated",
                          "MMI 30mg/day",
                          "MMI 60mg/day (Block)",
                          "Radioiodine 15 mCi",
                          "MMI + Propranolol",
                          "PTU 300mg/day",
                          "Rituximab + MMI"),
                selected=c("Untreated","MMI 30mg/day","Radioiodine 15 mCi")),
              hr(),
              radioButtons("comp_var", "Compare Biomarker:",
                choices=c("TSH_obs","fT4_obs","fT3_obs","TRAb_obs",
                          "HR_obs","BoneResor","GO_act","ThyMass"),
                selected="fT4_obs"),
              actionButton("run_compare", "Run Comparison",
                           class="btn-success btn-block", icon=icon("sync"))
          ),
          box(title="Scenario Comparison Plot", width=8, solidHeader=TRUE,
              status="info", plotlyOutput("scenario_plot", height="500px"))
        ),
        fluidRow(
          box(title="Endpoint Summary Table", width=12, solidHeader=TRUE,
              DTOutput("scenario_table"))
        )
      )
    )
  )
)

# ─── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Reactive simulation ────────────────────────────────────────────────────
  sim_data <- eventReactive(input$run_sim, {
    params <- list(
      USE_MMI  = as.integer(input$use_mmi),
      USE_PTU  = as.integer(input$use_ptu),
      USE_RAI  = as.integer(input$use_rai),
      USE_PROP = as.integer(input$use_prop),
      USE_GCS  = as.integer(input$use_gcs),
      USE_RTX  = as.integer(input$use_rtx)
    )
    inits <- list(
      TRAb_C  = input$TRAb_init,
      Bcell_C = 10,
      TSH_C   = 0.01,
      fT4_C   = input$fT4_init,
      fT3_C   = input$fT4_init * 0.34,
      T4_C    = input$fT4_init * 6.3,
      T3_C    = 5.0,
      rT3_C   = 0.3,
      ThyMass = 1.5,
      BoneResor = 1.6,
      HR_dev  = 25,
      GO_act  = input$TRAb_init * 1.4
    )
    run_sim(gd_mod, params, inits,
            duration  = input$sim_dur,
            mmi_dose  = ifelse(input$use_mmi,  input$dose_mmi,  0),
            ptu_dose  = ifelse(input$use_ptu,  input$dose_ptu,  0),
            rai_dose  = ifelse(input$use_rai,  input$dose_rai,  0),
            prop_dose = ifelse(input$use_prop, input$dose_prop, 0))
  }, ignoreNULL=FALSE)

  # Initial sim on startup
  observe({
    if(is.null(sim_data())) isolate(shinyjs::click("run_sim"))
  })

  last_row <- reactive({ tail(sim_data(), 1) })

  # ── Value Boxes ─────────────────────────────────────────────────────────────
  output$vbox_TSH <- renderValueBox({
    tsh <- round(last_row()$TSH_obs, 3)
    col <- if(tsh < 0.4) "red" else if(tsh > 4.5) "orange" else "green"
    valueBox(paste(tsh, "mIU/L"), "TSH", icon=icon("thermometer"), color=col)
  })
  output$vbox_fT4 <- renderValueBox({
    ft4 <- round(last_row()$fT4_obs, 1)
    col <- if(ft4 > 22) "red" else if(ft4 < 12) "orange" else "green"
    valueBox(paste(ft4, "pmol/L"), "Free T4", icon=icon("flask"), color=col)
  })
  output$vbox_TRAb <- renderValueBox({
    tr <- round(last_row()$TRAb_obs, 1)
    col <- if(tr > 1.75) "red" else "green"
    valueBox(paste(tr, "IU/L"), "TRAb", icon=icon("virus"), color=col)
  })
  output$vbox_HR <- renderValueBox({
    hr <- round(last_row()$HR_obs, 0)
    col <- if(hr > 100) "red" else if(hr < 60) "orange" else "green"
    valueBox(paste(hr, "bpm"), "Heart Rate", icon=icon("heartbeat"), color=col)
  })

  # ── Overview plot ────────────────────────────────────────────────────────────
  output$overview_plot <- renderPlotly({
    df <- sim_data() %>% mutate(months = TIME/30)
    p <- plot_ly(df, x=~months) %>%
      add_lines(y=~fT4_obs,  name="fT4 (pmol/L)",  line=list(color="orange")) %>%
      add_lines(y=~TSH_obs,  name="TSH (mIU/L)",   line=list(color="blue", dash="dash")) %>%
      add_lines(y=~TRAb_obs, name="TRAb (IU/L)",   line=list(color="red")) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="Concentration"),
             legend=list(orientation="h"),
             shapes=list(
               list(type="rect", x0=0, x1=max(df$months), y0=12, y1=22,
                    fillcolor="rgba(0,200,0,0.08)", line=list(width=0))
             ))
    p
  })

  # ── PK Plots ─────────────────────────────────────────────────────────────────
  output$pk_plot_mmi  <- renderPlotly({
    df <- sim_data() %>% mutate(months=TIME/30)
    plot_ly(df, x=~months, y=~MMI, type="scatter", mode="lines",
            line=list(color="purple")) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="MMI conc (AU)"))
  })
  output$pk_plot_prop <- renderPlotly({
    df <- sim_data() %>% mutate(months=TIME/30)
    plot_ly(df, x=~months, y=~PROP_C, type="scatter", mode="lines",
            line=list(color="#2980b9")) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="Propranolol (AU)"))
  })
  output$pk_plot_rai  <- renderPlotly({
    df <- sim_data() %>% mutate(months=TIME/30)
    plot_ly(df, x=~months) %>%
      add_lines(y=~RAI_SERUM, name="Serum ¹³¹I", line=list(color="orange")) %>%
      add_lines(y=~RAI_THYR,  name="Thyroid ¹³¹I", line=list(color="red")) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="¹³¹I (AU)"))
  })
  output$pk_plot_tpo  <- renderPlotly({
    df <- sim_data() %>% mutate(months=TIME/30, pct=TPO_inh*100)
    plot_ly(df, x=~months, y=~pct, type="scatter", mode="lines",
            line=list(color="darkgreen")) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="TPO Inhibition (%)"))
  })

  # ── PK Parameters Table ───────────────────────────────────────────────────────
  output$pk_param_table <- renderDT({
    data.frame(
      Drug       = c("Methimazole (MMI)","PTU","Radioiodine ¹³¹I","Propranolol","Levothyroxine"),
      Bioavail   = c("93%","75%","~100%","26%","70-80%"),
      Tmax       = c("1h","1-2h","<1h","1-2h","2-4h"),
      Half_life  = c("4-6h","1-2h","8d (phys)","4h","7d"),
      Vd         = c("35L","20L","body","3.9 L/kg","10L"),
      Mechanism  = c("TPO inhibition + immunosuppression",
                     "TPO + D1 inhibition","β-radiation thyroid cell kill",
                     "β1-blocker + D1 inhibition","T4 replacement"),
      stringsAsFactors=FALSE
    ) %>% datatable(options=list(dom="t"), rownames=FALSE)
  })

  # ── Hormone Plots ─────────────────────────────────────────────────────────────
  make_hormone_plot <- function(var, title, reflo=NULL, refhi=NULL, col="steelblue") {
    renderPlotly({
      df <- sim_data() %>% mutate(months=TIME/30)
      p  <- plot_ly(df, x=~months, y=~.data[[var]], type="scatter", mode="lines",
                    line=list(color=col)) %>%
        layout(xaxis=list(title="Time (months)"), yaxis=list(title=title))
      if(!is.null(reflo)) p <- p %>% add_lines(y=reflo, name="Lower ref", line=list(dash="dash", color="grey"))
      if(!is.null(refhi)) p <- p %>% add_lines(y=refhi, name="Upper ref", line=list(dash="dash", color="grey"))
      p
    })
  }

  output$pd_tsh  <- make_hormone_plot("TSH_obs", "TSH (mIU/L)",  0.4, 4.5, "royalblue")
  output$pd_fT4  <- make_hormone_plot("fT4_obs", "fT4 (pmol/L)", 12,  22,  "darkorange")
  output$pd_fT3  <- make_hormone_plot("fT3_obs", "fT3 (pmol/L)", 3.5, 7.5, "goldenrod")
  output$pd_rT3  <- make_hormone_plot("rT3_C",   "rT3 (AU)",     col="grey40")
  output$pd_thyroid_mass <- make_hormone_plot("ThyMass","Thyroid Mass (rel)", 0.8, 1.2, "tomato")

  output$pd_tot_T4T3 <- renderPlotly({
    df <- sim_data() %>% mutate(months=TIME/30)
    plot_ly(df, x=~months) %>%
      add_lines(y=~T4_C, name="Total T4", line=list(color="darkorange")) %>%
      add_lines(y=~T3_C, name="Total T3", line=list(color="goldenrod")) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="Concentration (nmol/L)"))
  })

  # ── CV & Bone Plots ───────────────────────────────────────────────────────────
  output$cv_hr <- renderPlotly({
    df <- sim_data() %>% mutate(months=TIME/30)
    plot_ly(df, x=~months, y=~HR_obs, type="scatter", mode="lines",
            line=list(color="crimson")) %>%
      add_lines(y=100, name="Tachycardia cutoff", line=list(dash="dash",color="grey")) %>%
      add_lines(y=72,  name="Normal HR",           line=list(dash="dot", color="green")) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="Heart Rate (bpm)"))
  })
  output$cv_bone <- renderPlotly({
    df <- sim_data() %>% mutate(months=TIME/30)
    plot_ly(df, x=~months, y=~BoneResor, type="scatter", mode="lines",
            line=list(color="saddlebrown")) %>%
      add_lines(y=1.0, name="Normal", line=list(dash="dash", color="grey")) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="Bone Resorption Index"))
  })
  output$cv_bmd <- renderPlotly({
    df <- sim_data() %>% mutate(months=TIME/30,
                                 bmd_loss = cumsum(pmax(0, BoneResor-1)*0.005/365))
    plot_ly(df, x=~months, y=~bmd_loss, type="scatter", mode="lines",
            line=list(color="peru")) %>%
      layout(xaxis=list(title="Time (months)"),
             yaxis=list(title="Estimated BMD Loss (fraction)"))
  })
  output$cv_bone_table <- renderTable({
    lr <- last_row()
    data.frame(
      Parameter = c("Heart Rate","Tachycardia?","Bone Resorption Index","Fracture Risk"),
      Value     = c(round(lr$HR_obs,0),
                    ifelse(lr$HR_obs>100,"Yes","No"),
                    round(lr$BoneResor,2),
                    ifelse(lr$BoneResor>1.5,"Elevated","Normal"))
    )
  })

  # ── Immune & GO Plots ─────────────────────────────────────────────────────────
  output$imm_trab <- renderPlotly({
    df <- sim_data() %>% mutate(months=TIME/30)
    plot_ly(df, x=~months, y=~TRAb_obs, type="scatter", mode="lines",
            line=list(color="darkred")) %>%
      add_lines(y=1.75, name="Positive cutoff", line=list(dash="dash", color="red")) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="TRAb (IU/L)"))
  })
  output$imm_bcell <- renderPlotly({
    df <- sim_data() %>% mutate(months=TIME/30)
    plot_ly(df, x=~months, y=~Bcell_C, type="scatter", mode="lines",
            line=list(color="purple")) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="Disease B cells (AU)"))
  })
  output$imm_go <- renderPlotly({
    df <- sim_data() %>% mutate(months=TIME/30)
    plot_ly(df, x=~months, y=~GO_act, type="scatter", mode="lines",
            line=list(color="darkblue")) %>%
      add_lines(y=3, name="CAS≥3 (active)", line=list(dash="dash", color="orange")) %>%
      layout(xaxis=list(title="Time (months)"), yaxis=list(title="GO Orbital Score (AU)"))
  })
  output$imm_remission <- renderPlotly({
    df <- sim_data() %>% mutate(months=TIME/30,
                                 in_remission = as.integer(TRAb_obs < 1.75 & TSH_obs > 0.4 & TSH_obs < 4.5))
    plot_ly(df, x=~months, y=~in_remission, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(39,174,96,0.3)",
            line=list(color="#27ae60")) %>%
      layout(xaxis=list(title="Time (months)"),
             yaxis=list(title="Remission (1=Yes, 0=No)", range=c(0,1.5)))
  })

  # ── Scenario Comparison ───────────────────────────────────────────────────────
  scenario_sims <- eventReactive(input$run_compare, {
    dur   <- input$sim_dur
    inits <- list(TRAb_C=25, Bcell_C=10, TSH_C=0.01, fT4_C=35, fT3_C=12,
                  T4_C=220, T3_C=5.0, rT3_C=0.3, ThyMass=1.5,
                  BoneResor=1.6, HR_dev=25, GO_act=35)

    scen_list <- list(
      "Untreated"           = list(params=list(USE_MMI=0,USE_PTU=0,USE_RAI=0),
                                   mmi=0,ptu=0,rai=0,prop=0),
      "MMI 30mg/day"        = list(params=list(USE_MMI=1), mmi=30, ptu=0, rai=0, prop=0),
      "MMI 60mg/day (Block)"= list(params=list(USE_MMI=1), mmi=60, ptu=0, rai=0, prop=0),
      "Radioiodine 15 mCi"  = list(params=list(USE_RAI=1), mmi=0,  ptu=0, rai=15,prop=0),
      "MMI + Propranolol"   = list(params=list(USE_MMI=1,USE_PROP=1), mmi=30, ptu=0, rai=0, prop=80),
      "PTU 300mg/day"       = list(params=list(USE_PTU=1), mmi=0, ptu=300, rai=0, prop=0),
      "Rituximab + MMI"     = list(params=list(USE_MMI=1,USE_RTX=1), mmi=30, ptu=0, rai=0, prop=0)
    )
    sel <- input$scenarios_sel
    lapply(sel, function(nm) {
      s <- scen_list[[nm]]
      run_sim(gd_mod, s$params, inits, dur, s$mmi, s$ptu, s$rai, s$prop) %>%
        mutate(scenario=nm, months=TIME/30)
    }) %>% bind_rows()
  })

  output$scenario_plot <- renderPlotly({
    df  <- scenario_sims()
    var <- input$comp_var
    p   <- plot_ly()
    for(scen in unique(df$scenario)) {
      sub <- df %>% filter(scenario == scen)
      p   <- p %>% add_lines(data=sub, x=~months, y=as.formula(paste0("~",var)),
                              name=scen)
    }
    p %>% layout(xaxis=list(title="Time (months)"),
                 yaxis=list(title=var),
                 legend=list(orientation="h"))
  })

  output$scenario_table <- renderDT({
    df  <- scenario_sims()
    var <- input$comp_var
    df %>%
      group_by(scenario) %>%
      summarise(
        Peak   = round(max(.data[[var]], na.rm=TRUE), 2),
        Nadir  = round(min(.data[[var]], na.rm=TRUE), 2),
        Final  = round(last(.data[[var]]), 2),
        AUC    = round(sum(.data[[var]]) * mean(diff(TIME)), 1)
      ) %>% datatable(options=list(dom="t"), rownames=FALSE)
  })

  # ── Lab Table ────────────────────────────────────────────────────────────────
  output$lab_table <- renderTable({
    lr <- last_row()
    data.frame(
      Biomarker = c("TSH","fT4","fT3","TRAb","Heart Rate","Thyroid Mass","Bone Resorption","GO Score"),
      Value     = c(sprintf("%.3f mIU/L", lr$TSH_obs),
                    sprintf("%.1f pmol/L", lr$fT4_obs),
                    sprintf("%.1f pmol/L", lr$fT3_obs),
                    sprintf("%.1f IU/L",   lr$TRAb_obs),
                    sprintf("%.0f bpm",    lr$HR_obs),
                    sprintf("%.2f (rel)",  lr$ThyMass),
                    sprintf("%.2f (AU)",   lr$BoneResor),
                    sprintf("%.1f (AU)",   lr$GO_act)),
      Normal   = c("0.4-4.5","12-22","3.5-7.5","<1.75","60-100","~1.0","1.0","<5"),
      Status   = c(
        ifelse(lr$TSH_obs<0.4, "LOW", ifelse(lr$TSH_obs>4.5,"HIGH","Normal")),
        ifelse(lr$fT4_obs>22,  "HIGH","Normal"),
        ifelse(lr$fT3_obs>7.5, "HIGH","Normal"),
        ifelse(lr$TRAb_obs>1.75,"POSITIVE","Negative"),
        ifelse(lr$HR_obs>100,   "TACHY","Normal"),
        ifelse(lr$ThyMass>1.2,  "Enlarged","Normal"),
        ifelse(lr$BoneResor>1.3,"Elevated","Normal"),
        ifelse(lr$GO_act>5,     "Active GO","Inactive")
      )
    )
  })

  # ── Diagnosis / Therapy status ────────────────────────────────────────────────
  output$diag_status <- renderUI({
    lr <- last_row()
    status <- if(lr$TSH_obs < 0.1 && lr$fT4_obs > 22 && lr$TRAb_obs > 1.75) {
      tags$div(tags$span("Active Graves' Hyperthyroidism", class="label label-danger",
                         style="font-size:14px"))
    } else if(lr$TSH_obs > 0.4 && lr$TSH_obs < 4.5) {
      tags$div(tags$span("Euthyroid (controlled)", class="label label-success",
                         style="font-size:14px"))
    } else {
      tags$div(tags$span("Borderline / Subclinical", class="label label-warning",
                         style="font-size:14px"))
    }
    status
  })
  output$therapy_list <- renderUI({
    therapies <- c()
    if(input$use_mmi)  therapies <- c(therapies, paste("MMI", input$dose_mmi, "mg/day"))
    if(input$use_ptu)  therapies <- c(therapies, paste("PTU", input$dose_ptu, "mg/day"))
    if(input$use_rai)  therapies <- c(therapies, paste("¹³¹I", input$dose_rai, "mCi"))
    if(input$use_prop) therapies <- c(therapies, paste("Propranolol", input$dose_prop, "mg/day"))
    if(input$use_gcs)  therapies <- c(therapies, "Glucocorticoids")
    if(input$use_rtx)  therapies <- c(therapies, "Rituximab")
    if(length(therapies)==0) therapies <- "None (untreated)"
    tags$ul(lapply(therapies, tags$li))
  })

  # Placeholder radar
  output$radar_plot <- renderPlot({
    lr   <- last_row()
    vars <- c("TRAb","fT4","HR","BoneRes","GO")
    vals <- c(min(1, lr$TRAb_obs/30), min(1, lr$fT4_obs/60),
              min(1,(lr$HR_obs-60)/80), min(1,(lr$BoneResor-1)/2),
              min(1, lr$GO_act/60))
    df   <- data.frame(var=factor(vars, levels=vars), val=vals, norm=0.4)
    ggplot(df, aes(x=var, y=val)) +
      geom_col(fill="#e74c3c", alpha=0.7, width=0.6) +
      geom_hline(yintercept=0.4, linetype="dashed", color="#27ae60") +
      coord_polar() +
      scale_y_continuous(limits=c(0,1)) +
      theme_minimal(base_size=11) +
      labs(title="Disease Activity\n(normalized 0–1)")
  })
}

shinyApp(ui, server)
