## ============================================================
## Gastroparesis QSP Interactive Dashboard — Shiny App
## 7 Tabs: Patient Profile · PK · PD · GI Motility ·
##         Clinical Endpoints · Scenario Comparison · Biomarkers
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(plotly)

# ---- Embedded mini-model (analytic approximation for Shiny) ----
simulate_gp <- function(
    DM_flag   = 1,
    HbA1c_0   = 9.5,
    ICC_0     = 0.40,
    nNOS_0    = 0.35,
    Pyloric_0 = 0.70,
    Antral_0  = 0.38,
    Meal_size = 300,
    Dose_MCP  = 0,     # mg
    Dose_DOM  = 0,     # mg
    Dose_ERY  = 0,     # mg
    Dose_PRU  = 0,     # mg
    Dose_REL  = 0,     # mcg
    Dose_OND  = 0,     # mg
    sim_days  = 7
) {
  dt   <- 0.25          # h step
  tend <- sim_days * 24
  t    <- seq(0, tend, by = dt)
  n    <- length(t)

  # PK parameters
  PK <- list(
    MCP = list(ka=1.2, CL=70, V=110, F=0.75, t_half=5.0),
    DOM = list(ka=0.8, CL=50, V=140, F=0.15, t_half=7.7),
    ERY = list(ka=0.6, CL=35, V=60,  F=0.35, t_half=1.5),
    PRU = list(ka=0.9, CL=18, V=400, F=0.90, t_half=26.0),
    REL = list(ka=0.5, CL=8,  V=20,  F=1.00, t_half=1.7),
    OND = list(ka=1.0, CL=40, V=160, F=0.60, t_half=3.0)
  )

  # Dosing intervals
  ii <- list(MCP=6, DOM=8, ERY=8, PRU=24, REL=12, OND=8)
  doses <- list(
    MCP = Dose_MCP,  DOM = Dose_DOM, ERY = Dose_ERY,
    PRU = Dose_PRU,  REL = Dose_REL, OND = Dose_OND
  )

  # Pre-compute drug concentration profiles (1-compartment oral)
  drug_conc <- function(drug_name) {
    d  <- doses[[drug_name]]
    pk <- PK[[drug_name]]
    iv <- ii[[drug_name]]
    if (d == 0) return(rep(0, n))
    n_doses <- floor(tend / iv) + 1
    dose_times <- seq(0, by = iv, length.out = n_doses)
    C <- rep(0, n)
    ka <- pk$ka; ke <- log(2) / pk$t_half
    Vd <- pk$V; F  <- pk$F
    for (td in dose_times) {
      dt_v <- t - td
      idx  <- which(dt_v >= 0)
      if (length(idx) == 0) next
      C[idx] <- C[idx] + (d * F * ka) / (Vd * (ka - ke)) *
                          (exp(-ke * dt_v[idx]) - exp(-ka * dt_v[idx]))
    }
    pmax(C * 1000 / Vd, 0)   # ng/mL
  }

  C_MCP <- drug_conc("MCP")
  C_DOM <- drug_conc("DOM")
  C_ERY <- drug_conc("ERY")
  C_PRU <- drug_conc("PRU")
  C_REL <- drug_conc("REL")
  C_OND <- drug_conc("OND")

  # PD effects at each time point
  E_D2  <- pmax(pmin(
    1 - (1 - 0.85*C_MCP/(5 + C_MCP)) * (1 - 0.80*C_DOM/(3 + C_DOM)), 1), 0)
  E_5HT4<- pmax(pmin(
    0.50*C_MCP/(50+C_MCP) + 0.92*C_PRU/(3+C_PRU) *
    (1 - 0.50*C_MCP/(50+C_MCP)), 1), 0)
  E_Mot <- 0.75 * C_ERY / (800 + C_ERY)
  E_GHSR<- 0.80 * C_REL / (15  + C_REL)
  E_5HT3<- 0.88 * C_OND / (8   + C_OND)  # antiemetic

  # Disease state (slow ODE, Euler)
  ICC  <- rep(ICC_0,    n)
  nNOS <- rep(nNOS_0,   n)
  k_nNOS_deg <- 0.0003 * DM_flag
  k_nNOS_syn <- 0.0002
  k_ICC_deg  <- 0.0002 * DM_flag
  k_ICC_rec  <- 0.0001
  for (i in 2:n) {
    nNOS[i] <- nNOS[i-1] + dt*(k_nNOS_syn - k_nNOS_deg*nNOS[i-1])
    ICC[i]  <- ICC[i-1]  + dt*(k_ICC_rec*(1-ICC[i-1]) - k_ICC_deg*ICC[i-1])
    nNOS[i] <- pmax(pmin(nNOS[i], 1), 0)
    ICC[i]  <- pmax(pmin(ICC[i],  1), 0)
  }

  # Antral & pyloric (fast, use direct formula)
  Antral <- pmax(pmin(
    Antral_0 + 0.30*E_D2 + 0.40*E_5HT4 + 0.35*E_Mot + 0.30*E_GHSR, 1), 0)
  Pyloric<- pmax(pmin(Pyloric_0 - 0.25*E_5HT4, 1), 0)

  # Meal schedule (every 6h starting t=0)
  meal_times <- seq(0, tend-1, by=6)
  GasVol <- rep(0, n)
  for (tm in meal_times) {
    idx <- which(abs(t - tm) < dt/2)[1]
    if (!is.na(idx)) GasVol[idx] <- GasVol[idx] + 300
  }

  # Propagate gastric emptying
  BG_pen <- DM_flag * (HbA1c_0 - 7.5) * 0.02
  for (i in 2:n) {
    ICC_eff  <- 0.5 + 0.5*ICC[i]
    nNOS_eff <- 0.5 + 0.5*nNOS[i]
    A  <- Antral[i]; P <- Pyloric[i]
    GER <- 0.025 * ICC_eff * nNOS_eff *
           (A^2 / (0.25 + A^2)) / (1 + P/0.4) * (1 - BG_pen)
    GER <- pmax(GER, 0)
    GasVol[i] <- pmax(GasVol[i] + GasVol[i-1] * (1 - GER*dt), 0)
  }

  # GCSI
  Ret   <- pmin(GasVol / 300, 1)
  Nausea_comp <- (1 - E_5HT3) * (1 - 0.3*E_D2)
  GCSI  <- pmin(0.40*5*Ret + 0.30*3*Nausea_comp + 0.30*5*(1-Antral), 5)

  tibble(
    time     = t,
    C_MCP    = C_MCP, C_DOM = C_DOM, C_ERY = C_ERY,
    C_PRU    = C_PRU, C_REL = C_REL, C_OND = C_OND,
    D2_occ   = E_D2  * 100,
    HT4_act  = E_5HT4 * 100,
    Mot_act  = E_Mot  * 100,
    GHSR_act = E_GHSR * 100,
    Antral   = Antral * 100,
    Pyloric  = Pyloric* 100,
    ICC_pct  = ICC    * 100,
    nNOS_pct = nNOS   * 100,
    GasVol   = GasVol,
    Ret4h    = Ret    * 100,
    GCSI     = GCSI
  )
}

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "Gastroparesis QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName="tab_patient",   icon=icon("user-md")),
      menuItem("Pharmacokinetics",   tabName="tab_PK",        icon=icon("flask")),
      menuItem("Pharmacodynamics",   tabName="tab_PD",        icon=icon("pills")),
      menuItem("GI Motility",        tabName="tab_motility",  icon=icon("wave-square")),
      menuItem("Clinical Endpoints", tabName="tab_endpoints", icon=icon("chart-bar")),
      menuItem("Scenario Comparison",tabName="tab_scenarios", icon=icon("balance-scale")),
      menuItem("Biomarkers",         tabName="tab_biomarkers",icon=icon("dna"))
    ),
    hr(),
    h5("Treatment Selection", style="padding-left:15px; color:#aaa"),
    checkboxGroupInput("drugs",
      label = NULL,
      choices = c("Metoclopramide 10mg QID" = "MCP",
                  "Domperidone 10mg TID"    = "DOM",
                  "Erythromycin 250mg TID"  = "ERY",
                  "Prucalopride 2mg QD"     = "PRU",
                  "Relamorelin 100mcg BID"  = "REL",
                  "Ondansetron 8mg TID"     = "OND"),
      selected = character(0)
    ),
    hr(),
    sliderInput("sim_days","Simulation (days)",1,28,7,step=1),
    actionButton("run","Run Simulation", icon=icon("play"),
                 style="margin:10px; background:#2E7D32; color:white; width:85%")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper, .right-side { background-color: #f9f9f9; }
      .box { border-radius: 8px; }
      .info-box { border-radius: 8px; }
    "))),
    tabItems(

      ## TAB 1: Patient Profile
      tabItem("tab_patient",
        fluidRow(
          box(title="Disease Subtype", width=4, status="primary",
            radioButtons("dm_flag","Gastroparesis Type",
              choices=c("Diabetic (DM)"=1, "Idiopathic"=0), selected=1),
            conditionalPanel("input.dm_flag==1",
              sliderInput("hba1c","HbA1c (%)", 5.5, 14, 9.5, 0.5)
            )
          ),
          box(title="Disease Severity", width=4, status="warning",
            sliderInput("icc0","ICC Density (% of normal)",10,90,40,5),
            sliderInput("nnos0","nNOS Activity (% of normal)",10,90,35,5)
          ),
          box(title="Gastric Motor Baseline", width=4, status="danger",
            sliderInput("antral0","Antral Contractility (%)",10,90,38,5),
            sliderInput("pyloric0","Pyloric Tone (0=open, 100=tight)",
                        20,90,70,5),
            sliderInput("meal_size","Meal Volume (mL)",100,600,300,50)
          )
        ),
        fluidRow(
          valueBoxOutput("vb_gcsi_base",  width=3),
          valueBoxOutput("vb_icc_base",   width=3),
          valueBoxOutput("vb_nnos_base",  width=3),
          valueBoxOutput("vb_ger_base",   width=3)
        ),
        fluidRow(
          box(title="Disease Severity Radar", width=6, status="primary",
            plotOutput("radar_disease", height=300)
          ),
          box(title="Clinical Context", width=6, status="info",
            HTML("<h4>Gastroparesis Diagnostic Criteria</h4>
            <ul>
              <li>Gastric retention &gt;10% at 4h on nuclear scintigraphy</li>
              <li>GCSI total score ≥2.0 (scale 0–5)</li>
              <li>Exclusion of mechanical obstruction</li>
            </ul>
            <h4>Disease Severity Classification (Abell et al.)</h4>
            <table border='1' cellpadding='4' width='100%'>
              <tr><th>Grade</th><th>4h Retention</th><th>GCSI</th></tr>
              <tr><td>Mild</td><td>10–20%</td><td>≤2.5</td></tr>
              <tr><td>Moderate</td><td>20–35%</td><td>2.5–3.5</td></tr>
              <tr><td>Severe</td><td>&gt;35%</td><td>&gt;3.5</td></tr>
            </table>
            <br>
            <b>ICC Loss</b>: Normal ICC density ~8,000/cm². DM patients retain
            ~40–60% after 10y. nNOS neurons absent/depleted in ~85% of
            gastroparesis biopsies (Grover et al., Gastroenterology 2011).")
          )
        )
      ),

      ## TAB 2: Pharmacokinetics
      tabItem("tab_PK",
        fluidRow(
          box(title="Plasma Concentration Profiles", width=12, status="primary",
            plotlyOutput("pk_plot", height=400))
        ),
        fluidRow(
          box(title="PK Parameter Summary", width=6, status="info",
            tableOutput("pk_table")
          ),
          box(title="Drug PK Properties", width=6, status="warning",
            HTML("<table border='1' cellpadding='5' width='100%'>
              <tr><th>Drug</th><th>Class</th><th>t½</th><th>F%</th><th>BBB</th></tr>
              <tr><td>Metoclopramide</td><td>D2/5HT4</td><td>5h</td><td>75%</td><td>Yes</td></tr>
              <tr><td>Domperidone</td><td>D2 (periph)</td><td>7.7h</td><td>15%</td><td>No</td></tr>
              <tr><td>Erythromycin</td><td>Motilin agonist</td><td>1.5h</td><td>35%</td><td>Limited</td></tr>
              <tr><td>Prucalopride</td><td>5-HT4 (full)</td><td>26h</td><td>90%</td><td>Limited</td></tr>
              <tr><td>Relamorelin</td><td>Ghrelin/GHSR</td><td>1.7h</td><td>SC</td><td>No</td></tr>
              <tr><td>Ondansetron</td><td>5-HT3 block</td><td>3h</td><td>60%</td><td>Limited</td></tr>
            </table>")
          )
        )
      ),

      ## TAB 3: Pharmacodynamics
      tabItem("tab_PD",
        fluidRow(
          box(title="Receptor Occupancy / Activation", width=12, status="success",
            plotlyOutput("pd_receptor_plot", height=350))
        ),
        fluidRow(
          box(title="Mechanism of Action Summary", width=6,
            HTML("<h4>Prokinetic Mechanisms</h4>
            <ul>
              <li><b>Metoclopramide</b>: D2 block → ↑ ACh release; partial 5-HT4 agonist</li>
              <li><b>Domperidone</b>: Peripheral D2 block only (no CNS EPS)</li>
              <li><b>Erythromycin</b>: Motilin receptor agonist → Phase III MMC</li>
              <li><b>Prucalopride</b>: Full 5-HT4 agonist → ↑ ACh, antral contractions</li>
              <li><b>Relamorelin</b>: GHSR (ghrelin) agonist → antral contractility ↑</li>
            </ul>
            <h4>Antiemetic Mechanisms</h4>
            <ul>
              <li><b>Ondansetron</b>: 5-HT3 block at CTZ → ↓ nausea/vomiting</li>
              <li><b>Metoclopramide</b>: D2+5HT3 block at CTZ → antiemetic</li>
            </ul>")
          ),
          box(title="PD Effect vs Concentration", width=6,
            plotlyOutput("pd_emax_plot", height=300))
        )
      ),

      ## TAB 4: GI Motility
      tabItem("tab_motility",
        fluidRow(
          box(title="Gastric Volume (Content) Over Time", width=6, status="primary",
            plotlyOutput("gasvol_plot", height=300)),
          box(title="Gastric Emptying Rate (%/h)", width=6, status="success",
            plotlyOutput("ger_plot", height=300))
        ),
        fluidRow(
          box(title="Antral Contractility (%)", width=6,
            plotlyOutput("antral_plot", height=280)),
          box(title="Pyloric Tone (%)", width=6,
            plotlyOutput("pyloric_plot", height=280))
        )
      ),

      ## TAB 5: Clinical Endpoints
      tabItem("tab_endpoints",
        fluidRow(
          valueBoxOutput("vb_gcsi",    width=3),
          valueBoxOutput("vb_ret4h",   width=3),
          valueBoxOutput("vb_nausea",  width=3),
          valueBoxOutput("vb_d2",      width=3)
        ),
        fluidRow(
          box(title="GCSI Score Over Time", width=6, status="danger",
            plotlyOutput("gcsi_plot", height=300)),
          box(title="4-Hour Gastric Retention (%)", width=6, status="warning",
            plotlyOutput("ret4h_plot", height=300))
        ),
        fluidRow(
          box(title="Response Thresholds (ROME/ACG Guidelines)", width=12,
            HTML("<table border='1' cellpadding='6' width='100%'>
              <tr><th>Measure</th><th>Normal</th><th>Mild GP</th><th>Moderate GP</th><th>Severe GP</th></tr>
              <tr><td>4h Retention</td><td>&lt;10%</td><td>10-20%</td><td>20-35%</td><td>&gt;35%</td></tr>
              <tr><td>GCSI Score</td><td>&lt;1.0</td><td>1-2.5</td><td>2.5-3.5</td><td>&gt;3.5</td></tr>
              <tr><td>Nausea VAS</td><td>&lt;10mm</td><td>10-30mm</td><td>30-60mm</td><td>&gt;60mm</td></tr>
              <tr><td>GER T½</td><td>60-90min</td><td>90-120min</td><td>120-180min</td><td>&gt;180min</td></tr>
            </table>")
          )
        )
      ),

      ## TAB 6: Scenario Comparison
      tabItem("tab_scenarios",
        fluidRow(
          box(title="Multi-Scenario GCSI Comparison", width=12, status="primary",
            plotlyOutput("scenario_gcsi_plot", height=350))
        ),
        fluidRow(
          box(title="Run All 7 Pre-defined Scenarios",
              width=12, status="info",
              actionButton("run_all_scenarios","Run All Scenarios",
                           icon=icon("sync"),
                           style="background:#1565C0; color:white"),
              br(), br(),
              tableOutput("scenario_table"))
        )
      ),

      ## TAB 7: Biomarkers
      tabItem("tab_biomarkers",
        fluidRow(
          box(title="ICC Density & nNOS Activity — Disease Progression", width=12,
            plotlyOutput("biomarker_plot", height=350))
        ),
        fluidRow(
          box(title="Gastroparesis Biomarker Panel", width=6, status="warning",
            HTML("<h4>Structural Biomarkers</h4>
            <ul>
              <li><b>ICC Density</b> (c-Kit+ on full-thickness biopsy): &lt;50% normal = depleted</li>
              <li><b>nNOS Neurons</b> (myenteric plexus): absent in 85% of GP biopsies</li>
              <li><b>CD206+ M2 Macrophages</b>: protective; reduced in DM-GP</li>
              <li><b>HO-1 Expression</b>: cytoprotective enzyme, impaired in GP</li>
            </ul>
            <h4>Functional Biomarkers</h4>
            <ul>
              <li><b>Gastric Scintigraphy (4h)</b>: Gold-standard &gt;10% = abnormal</li>
              <li><b>SmartPill / Wireless Motility Capsule</b>: pH, pressure, T</li>
              <li><b>13C-Octanoic Acid Breath Test</b>: non-radioactive GES</li>
            </ul>")
          ),
          box(title="Plasma / Serum Biomarkers", width=6, status="info",
            HTML("<h4>Laboratory Markers</h4>
            <ul>
              <li><b>Ghrelin</b>: ↓ in DM-gastroparesis; predicts response to REL</li>
              <li><b>Motilin</b>: Phase III MMC-related; may be ↓ in idiopathic GP</li>
              <li><b>HbA1c</b>: Primary driver of DM-GP severity; target &lt;7%</li>
              <li><b>C-Reactive Protein</b>: Systemic inflammation marker</li>
              <li><b>Nutritional Panel</b>: Albumin, pre-albumin, B12, D3</li>
              <li><b>Anti-Hu antibody</b>: Paraneoplastic GP screening</li>
            </ul>
            <h4>Emerging Biomarkers</h4>
            <ul>
              <li>Serum SCF (Stem Cell Factor) — ICC survival signal</li>
              <li>Plasma substance P — ENS inflammation marker</li>
              <li>Urinary 5-HIAA — 5-HT metabolism</li>
            </ul>")
          )
        )
      )
    )
  )
)

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

  # Reactive simulation
  sim_data <- eventReactive(input$run, {
    doses_on <- input$drugs
    simulate_gp(
      DM_flag   = as.numeric(input$dm_flag),
      HbA1c_0   = if(as.numeric(input$dm_flag)==1) input$hba1c else 5.5,
      ICC_0     = input$icc0 / 100,
      nNOS_0    = input$nnos0 / 100,
      Pyloric_0 = input$pyloric0 / 100,
      Antral_0  = input$antral0 / 100,
      Meal_size = input$meal_size,
      Dose_MCP  = if("MCP" %in% doses_on) 10 else 0,
      Dose_DOM  = if("DOM" %in% doses_on) 10 else 0,
      Dose_ERY  = if("ERY" %in% doses_on) 250 else 0,
      Dose_PRU  = if("PRU" %in% doses_on) 2 else 0,
      Dose_REL  = if("REL" %in% doses_on) 100 else 0,
      Dose_OND  = if("OND" %in% doses_on) 8 else 0,
      sim_days  = input$sim_days
    )
  }, ignoreNULL=FALSE)

  # Value boxes TAB1
  output$vb_gcsi_base <- renderValueBox({
    valueBox(round(3.2, 1), "Baseline GCSI", icon=icon("exclamation-triangle"),
             color=if(input$icc0<40)"red" else "yellow")
  })
  output$vb_icc_base <- renderValueBox({
    valueBox(paste0(input$icc0,"%"), "ICC Density", icon=icon("microscope"),
             color=if(input$icc0<40)"red" else "orange")
  })
  output$vb_nnos_base <- renderValueBox({
    valueBox(paste0(input$nnos0,"%"), "nNOS Activity", icon=icon("bolt"),
             color=if(input$nnos0<35)"red" else "orange")
  })
  output$vb_ger_base <- renderValueBox({
    ger_est <- round(0.025 * (0.5+input$icc0/200) * (0.5+input$nnos0/200) *
                     (input$antral0/100)^2 / (0.25+(input$antral0/100)^2) /
                     (1 + (input$pyloric0/100)/0.4) * 60, 1)
    valueBox(paste0(ger_est,"%/h"), "Est. GER", icon=icon("tachometer-alt"),
             color=if(ger_est<8)"red" else "green")
  })

  # Radar chart (simplified bar plot)
  output$radar_disease <- renderPlot({
    df <- data.frame(
      Measure = c("ICC Density","nNOS Activity","Antral Fxn","Pyloric Open","GER"),
      Value   = c(input$icc0, input$nnos0, input$antral0,
                  100-input$pyloric0,
                  pmax(pmin(0.025*(0.5+input$icc0/200)*(0.5+input$nnos0/200)*
                       (input$antral0/100)^2/(0.25+(input$antral0/100)^2)/
                       (1+(input$pyloric0/100)/0.4)*6000, 100), 0)),
      Normal  = c(100,100,100,70,100)
    )
    df_long <- df %>% pivot_longer(c(Value, Normal), names_to="Type", values_to="Val")
    ggplot(df_long, aes(x=Measure, y=Val, fill=Type)) +
      geom_col(position="dodge", alpha=0.8) +
      scale_fill_manual(values=c(Value="#EF5350",Normal="#42A5F5")) +
      labs(title="Patient vs Normal (%)", y="%", x=NULL, fill=NULL) +
      theme_minimal(base_size=12) + coord_flip()
  })

  # PK Plot
  output$pk_plot <- renderPlotly({
    d <- sim_data()
    if(is.null(d)) return(NULL)
    drugs_sel <- input$drugs
    if(length(drugs_sel)==0) {
      return(plot_ly() %>% layout(title="No drugs selected"))
    }
    fig <- plot_ly()
    cols <- c(MCP="#1565C0",DOM="#AD1457",ERY="#E65100",
              PRU="#2E7D32",REL="#6A1B9A",OND="#00838F")
    names_map <- c(MCP="Metoclopramide",DOM="Domperidone",ERY="Erythromycin",
                   PRU="Prucalopride",REL="Relamorelin",OND="Ondansetron")
    for(dr in drugs_sel) {
      col_nm <- paste0("C_", dr)
      if(col_nm %in% names(d)) {
        fig <- fig %>% add_lines(x=d$time, y=d[[col_nm]],
                                  name=names_map[dr], line=list(color=cols[dr]))
      }
    }
    fig %>% layout(xaxis=list(title="Time (h)"),
                   yaxis=list(title="Concentration (ng/mL)"),
                   legend=list(orientation="h"))
  })

  output$pk_table <- renderTable({
    data.frame(
      Drug     = c("Metoclopramide","Domperidone","Erythromycin",
                   "Prucalopride","Relamorelin","Ondansetron"),
      t_half_h = c(5.0, 7.7, 1.5, 26.0, 1.7, 3.0),
      F_pct    = c(75, 15, 35, 90, "SC", 60),
      Vd_L     = c(110, 140, 60, 400, 20, 160),
      CL_L_h   = c(70, 50, 35, 18, 8, 40)
    )
  }, striped=TRUE, bordered=TRUE)

  # PD Receptor Plot
  output$pd_receptor_plot <- renderPlotly({
    d <- sim_data()
    if(is.null(d)) return(NULL)
    d2 <- d %>% select(time, D2_occ, HT4_act, Mot_act, GHSR_act) %>%
      pivot_longer(-time, names_to="Receptor", values_to="Effect_pct")
    d2$Receptor <- recode(d2$Receptor,
      "D2_occ"  = "D2 Occupancy",
      "HT4_act" = "5-HT4 Activation",
      "Mot_act" = "Motilin Activation",
      "GHSR_act"= "Ghrelin/GHSR"
    )
    fig <- plot_ly(d2, x=~time, y=~Effect_pct, color=~Receptor,
                   type="scatter", mode="lines")
    fig %>% layout(xaxis=list(title="Time (h)"),
                   yaxis=list(title="Effect (%)", range=c(0,100)),
                   legend=list(orientation="h"))
  })

  output$pd_emax_plot <- renderPlotly({
    C_range <- seq(0, 200, length.out=200)
    df <- data.frame(
      C_range = C_range,
      MCP_D2  = 0.85 * C_range / (5 + C_range) * 100,
      PRU_5HT4= 0.92 * C_range / (3 + C_range) * 100,
      OND_5HT3= 0.88 * C_range / (8 + C_range) * 100
    )
    plot_ly(df, x=~C_range) %>%
      add_lines(y=~MCP_D2, name="MCP: D2 Block", line=list(color="#1565C0")) %>%
      add_lines(y=~PRU_5HT4, name="PRU: 5-HT4 Agonism", line=list(color="#2E7D32")) %>%
      add_lines(y=~OND_5HT3, name="OND: 5-HT3 Block", line=list(color="#00838F")) %>%
      layout(xaxis=list(title="Drug Concentration (ng/mL)"),
             yaxis=list(title="Effect (%)"),
             legend=list(orientation="h"))
  })

  # GI Motility plots
  output$gasvol_plot <- renderPlotly({
    d <- sim_data()
    if(is.null(d)) return(NULL)
    plot_ly(d, x=~time, y=~GasVol, type="scatter", mode="lines",
            line=list(color="#E65100")) %>%
      layout(xaxis=list(title="Time (h)"),
             yaxis=list(title="Gastric Volume (mL)"))
  })

  output$ger_plot <- renderPlotly({
    d <- sim_data()
    if(is.null(d)) return(NULL)
    ger_series <- c(0, diff(d$GasVol))
    ger_series[ger_series > 0] <- 0  # only emptying (negative dV)
    ger_pct <- -ger_series / pmax(d$GasVol, 1) * 400  # rough %/h
    plot_ly(x=d$time, y=pmax(ger_pct, 0), type="scatter", mode="lines",
            line=list(color="#2E7D32")) %>%
      layout(xaxis=list(title="Time (h)"),
             yaxis=list(title="GER (% meal/h)"))
  })

  output$antral_plot <- renderPlotly({
    d <- sim_data()
    if(is.null(d)) return(NULL)
    plot_ly(d, x=~time, y=~Antral, type="scatter", mode="lines",
            line=list(color="#1565C0")) %>%
      add_lines(y=rep(100,nrow(d)), name="Normal", line=list(dash="dash", color="grey")) %>%
      layout(xaxis=list(title="Time (h)"),
             yaxis=list(title="Antral Contractility (%)"))
  })

  output$pyloric_plot <- renderPlotly({
    d <- sim_data()
    if(is.null(d)) return(NULL)
    plot_ly(d, x=~time, y=~Pyloric, type="scatter", mode="lines",
            line=list(color="#AD1457")) %>%
      layout(xaxis=list(title="Time (h)"),
             yaxis=list(title="Pyloric Tone (0=open, 100=tight)"))
  })

  # Clinical Endpoints
  d_late <- reactive({
    d <- sim_data()
    if(is.null(d)) return(NULL)
    d %>% filter(time > max(time)*0.8)
  })

  output$vb_gcsi <- renderValueBox({
    d <- d_late()
    v <- if(!is.null(d)) round(mean(d$GCSI),2) else "--"
    col <- if(is.numeric(v) && v < 2.5) "green" else if(is.numeric(v) && v < 3.5) "yellow" else "red"
    valueBox(v, "GCSI Score", icon=icon("stethoscope"), color=col)
  })
  output$vb_ret4h <- renderValueBox({
    d <- d_late()
    v <- if(!is.null(d)) round(mean(d$Ret4h),1) else "--"
    col <- if(is.numeric(v) && v < 10) "green" else if(is.numeric(v) && v < 20) "yellow" else "red"
    valueBox(paste0(v,"%"), "4h Retention", icon=icon("hourglass-half"), color=col)
  })
  output$vb_nausea <- renderValueBox({
    d <- d_late()
    v <- if(!is.null(d)) round(mean(d$GCSI * 0.3 * 100/3), 0) else "--"
    valueBox(paste0(v,"mm"), "Nausea VAS", icon=icon("frown"), color="orange")
  })
  output$vb_d2 <- renderValueBox({
    d <- d_late()
    v <- if(!is.null(d)) round(mean(d$D2_occ),0) else "--"
    valueBox(paste0(v,"%"), "D2 Occupancy", icon=icon("pills"), color="blue")
  })

  output$gcsi_plot <- renderPlotly({
    d <- sim_data()
    if(is.null(d)) return(NULL)
    plot_ly(d, x=~time, y=~GCSI, type="scatter", mode="lines",
            line=list(color="#C62828")) %>%
      add_lines(y=rep(2.0,nrow(d)), name="Moderate threshold",
                line=list(dash="dash", color="orange")) %>%
      layout(xaxis=list(title="Time (h)"),
             yaxis=list(title="GCSI Score", range=c(0,5)))
  })

  output$ret4h_plot <- renderPlotly({
    d <- sim_data()
    if(is.null(d)) return(NULL)
    plot_ly(d, x=~time, y=~Ret4h, type="scatter", mode="lines",
            line=list(color="#E65100")) %>%
      add_lines(y=rep(10,nrow(d)), name="Normal limit",
                line=list(dash="dash", color="green")) %>%
      layout(xaxis=list(title="Time (h)"),
             yaxis=list(title="Gastric Retention (%)"))
  })

  # Scenario Comparison
  all_scenarios <- eventReactive(input$run_all_scenarios, {
    common <- list(
      DM_flag=1, HbA1c_0=9.5, ICC_0=0.40, nNOS_0=0.35,
      Pyloric_0=0.70, Antral_0=0.38, Meal_size=300, sim_days=7
    )
    scenarios <- list(
      list(name="S0: Untreated", doses=list()),
      list(name="S1: Metoclopramide",  doses=list(Dose_MCP=10)),
      list(name="S2: Domperidone",     doses=list(Dose_DOM=10)),
      list(name="S3: Erythromycin",    doses=list(Dose_ERY=250)),
      list(name="S4: Prucalopride",    doses=list(Dose_PRU=2)),
      list(name="S5: Relamorelin",     doses=list(Dose_REL=100)),
      list(name="S6: PRU+Ondansetron", doses=list(Dose_PRU=2, Dose_OND=8))
    )
    bind_rows(lapply(scenarios, function(sc) {
      args <- c(common, sc$doses)
      do.call(simulate_gp, args) %>%
        filter(time > 96) %>%
        summarise(
          Scenario    = sc$name,
          GCSI_mean   = round(mean(GCSI),2),
          Ret4h_mean  = round(mean(Ret4h),1),
          D2_occ_mean = round(mean(D2_occ),1),
          HT4_act_mean= round(mean(HT4_act),1)
        )
    }))
  })

  output$scenario_table <- renderTable({
    d <- all_scenarios()
    if(is.null(d)) return(NULL)
    d %>% rename(
      "GCSI (0-5)"       = GCSI_mean,
      "4h Retention (%)" = Ret4h_mean,
      "D2 Occ. (%)"      = D2_occ_mean,
      "5-HT4 Act. (%)"   = HT4_act_mean
    )
  }, striped=TRUE, bordered=TRUE)

  output$scenario_gcsi_plot <- renderPlotly({
    d <- all_scenarios()
    if(is.null(d)) return(NULL)
    plot_ly(d, x=~Scenario, y=~GCSI_mean, type="bar",
            color=~Scenario, showlegend=FALSE) %>%
      add_lines(y=rep(2.0, nrow(d)), name="Moderate threshold",
                line=list(dash="dash", color="orange"), showlegend=TRUE) %>%
      layout(xaxis=list(title=NULL, tickangle=-30),
             yaxis=list(title="GCSI Score (Day 5-7 mean)", range=c(0,5)))
  })

  # Biomarkers
  output$biomarker_plot <- renderPlotly({
    d <- sim_data()
    if(is.null(d)) return(NULL)
    plot_ly(d, x=~time) %>%
      add_lines(y=~ICC_pct,  name="ICC Density (%)", line=list(color="#1565C0")) %>%
      add_lines(y=~nNOS_pct, name="nNOS Activity (%)", line=list(color="#2E7D32")) %>%
      add_lines(y=rep(50,nrow(d)), name="50% threshold",
                line=list(dash="dash", color="red")) %>%
      layout(xaxis=list(title="Time (h)"),
             yaxis=list(title="% of Normal", range=c(0,100)),
             legend=list(orientation="h"))
  })
}

## ============================================================
## LAUNCH
## ============================================================
shinyApp(ui = ui, server = server)
