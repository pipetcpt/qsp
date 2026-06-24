## =============================================================================
## Hidradenitis Suppurativa (HS) — Shiny QSP Dashboard
## =============================================================================
## Tabs:
##   1. Patient Profile & Inputs
##   2. Drug PK Profiles
##   3. PD — Cytokine & Cell Dynamics
##   4. Clinical Endpoints (AN, IHS4, HiSCR)
##   5. Scenario Comparison
##   6. Biomarker & Virtual Population
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)
library(shinycssloaders)

## ---- mrgsolve model code (embedded) ----------------------------------------
model_code <- '
$PARAM
ADA_ka=0.280, ADA_F=0.640, ADA_Vc=3.20, ADA_Vp=4.10, ADA_CL=0.220, ADA_Q=0.500
ADA_IC50=0.200, ADA_Imax=0.980, ADA_ADA_effect=0.0
SEC_ka=0.300, SEC_F=0.730, SEC_Vc=7.20, SEC_CL=0.190
SEC_IC50=0.080, SEC_Imax=0.970
BIM_ka=0.280, BIM_F=0.740, BIM_Vc=5.40, BIM_CL=0.170
BIM_IC50=0.060, BIM_Imax=0.990
TNF_kin=0.500, TNF_kout=0.140, TNF_ss=3.571
IL17A_kin=0.200, IL17A_kout=0.100, IL17A_ss=2.000
IL6_kin=2.000, IL6_kout=0.500, IL6_amp_TNF=0.300, IL6_amp_IL17=0.200
IL1b_kin=0.100, IL1b_kout=0.200, IL1b_amp_TNF=0.250
IL23_kin=0.080, IL23_kout=0.120, IL23_amp=0.200
Th17_kin=0.050, Th17_kout=0.030, Th17_IL23=0.400, Th17_IL1b=0.300
M1_kin=0.080, M1_kout=0.050, M1_TNF=0.300, M1_IL17=0.200
Neu_kin=0.200, Neu_kout=0.400, Neu_IL8=0.500
AN_kin=0.020, AN_kout=0.010, AN_drive=0.800, AN_base=5.000
Fist_kin=0.002, Fist_kout=0.003
IHS4_w_AN=1.000, IHS4_w_ab=2.000, IHS4_w_fi=3.000
BW=80.0, SEX=0.0, HURLEY=2.0, SMOKE=1.0

$CMT ADA_abs ADA_Cc ADA_Cp SEC_abs SEC_Cc BIM_abs BIM_Cc
     TNF TNF_ADA IL17A IL17A_anti IL6 IL1b IL23
     Th17_idx M1_idx Neu_idx AN Fist IHS4

$MAIN
double ADA_C = ADA_Cc / ADA_Vc;
double SEC_C = SEC_Cc / SEC_Vc;
double BIM_C = BIM_Cc / BIM_Vc;
double I_ADA = ADA_Imax * ADA_C  / (ADA_IC50  + ADA_C  + 1e-12);
double I_SEC = SEC_Imax * SEC_C  / (SEC_IC50  + SEC_C  + 1e-12);
double I_BIM = BIM_Imax * BIM_C  / (BIM_IC50  + BIM_C  + 1e-12);
double I_IL17 = fmax(I_SEC, I_BIM);
double ADA_CL_eff = ADA_CL * (1.0 + ADA_ADA_effect);
double inflam_idx = (TNF/TNF_ss + IL17A/IL17A_ss + M1_idx + Neu_idx) / 4.0;
if(inflam_idx < 0) inflam_idx = 0;
double smoke_mult = 1.0 + 0.30 * SMOKE;

$ODE
dxdt_ADA_abs = -ADA_ka * ADA_abs;
dxdt_ADA_Cc  =  ADA_ka * ADA_F * ADA_abs - (ADA_CL_eff + ADA_Q)/ADA_Vc * ADA_Cc + ADA_Q/ADA_Vp * ADA_Cp;
dxdt_ADA_Cp  =  ADA_Q/ADA_Vc * ADA_Cc - ADA_Q/ADA_Vp * ADA_Cp;
dxdt_SEC_abs = -SEC_ka * SEC_abs;
dxdt_SEC_Cc  =  SEC_ka * SEC_F * SEC_abs - SEC_CL/SEC_Vc * SEC_Cc;
dxdt_BIM_abs = -BIM_ka * BIM_abs;
dxdt_BIM_Cc  =  BIM_ka * BIM_F * BIM_abs - BIM_CL/BIM_Vc * BIM_Cc;
double TNF_prod = TNF_kin * M1_idx * smoke_mult * (1.0 - I_ADA);
dxdt_TNF     = TNF_prod - TNF_kout * TNF;
dxdt_TNF_ADA = I_ADA * TNF_kin * M1_idx * smoke_mult - TNF_kout * TNF_ADA;
double IL17A_prod = IL17A_kin * Th17_idx * (1.0 + IL23/IL23_kin) * (1.0 - I_IL17);
dxdt_IL17A      = IL17A_prod - IL17A_kout * IL17A;
dxdt_IL17A_anti = I_IL17 * IL17A_kin * Th17_idx - IL17A_kout * IL17A_anti;
dxdt_IL6 = IL6_kin + IL6_amp_TNF*TNF/(TNF_ss+TNF) + IL6_amp_IL17*IL17A/(IL17A_ss+IL17A) - IL6_kout*IL6;
dxdt_IL1b = IL1b_kin*(1.0+IL1b_amp_TNF*TNF/TNF_ss)*M1_idx - IL1b_kout*IL1b;
dxdt_IL23 = IL23_kin*(1.0+IL23_amp*M1_idx) - IL23_kout*IL23;
double Th17_expand = Th17_IL23*IL23/(IL23_kin+IL23) + Th17_IL1b*IL1b/(IL1b_kin/IL1b_kout+IL1b);
dxdt_Th17_idx = Th17_kin*(1.0+Th17_expand) - Th17_kout*Th17_idx;
double M1_expand = M1_TNF*TNF/(TNF_ss+TNF) + M1_IL17*IL17A/(IL17A_ss+IL17A);
dxdt_M1_idx = M1_kin*(1.0+M1_expand)*smoke_mult - M1_kout*M1_idx;
dxdt_Neu_idx = Neu_kin*(1.0+Neu_IL8*IL17A/(IL17A_ss+IL17A)) - Neu_kout*Neu_idx;
dxdt_AN = AN_kin*inflam_idx*AN_drive - AN_kout*AN;
dxdt_Fist = Fist_kin*inflam_idx - Fist_kout*Fist;
dxdt_IHS4 = 0;

$TABLE
double C_ADA = ADA_Cc/ADA_Vc;
double C_SEC = SEC_Cc/SEC_Vc;
double C_BIM = BIM_Cc/BIM_Vc;
double I_ADA_t = ADA_Imax*C_ADA/(ADA_IC50+C_ADA+1e-12);
double I_SEC_t = SEC_Imax*C_SEC/(SEC_IC50+C_SEC+1e-12);
double I_BIM_t = BIM_Imax*C_BIM/(BIM_IC50+C_BIM+1e-12);
double IHS4_cur = IHS4_w_AN*AN + IHS4_w_fi*Fist;
double HiSCR = (AN <= AN_base*0.5) ? 1.0 : 0.0;
double VAS_pain = fmin(10.0, AN*0.8);
double DLQI = fmin(30.0, IHS4_cur*0.8 + VAS_pain*0.5);
double AN_pct = (AN - AN_base)/AN_base*100.0;

$CAPTURE C_ADA C_SEC C_BIM I_ADA_t I_SEC_t I_BIM_t
         TNF IL17A IL6 IL1b IL23
         Th17_idx M1_idx Neu_idx
         AN Fist IHS4_cur HiSCR VAS_pain DLQI AN_pct

$INIT
ADA_abs=0, ADA_Cc=0, ADA_Cp=0,
SEC_abs=0, SEC_Cc=0,
BIM_abs=0, BIM_Cc=0,
TNF=3.571, TNF_ADA=0,
IL17A=2.000, IL17A_anti=0,
IL6=8.000, IL1b=0.800, IL23=0.800,
Th17_idx=1.5, M1_idx=1.8, Neu_idx=2.0,
AN=5.0, Fist=2.0, IHS4=9.0
'

## Pre-compile model at startup
mod_global <- mcode("HS_Shiny", model_code, quiet = TRUE)

## Helper: build event table
make_events <- function(drug, ada_dose, sec_dose, bim_dose, duration, ada_ada_effect) {
  evs <- list()
  day_max <- duration * 7

  if ("adalimumab" %in% drug) {
    evs[[length(evs) + 1]] <- ev(amt = 160,      time = 0,  cmt = "ADA_abs")
    evs[[length(evs) + 1]] <- ev(amt = 80,        time = 14, cmt = "ADA_abs")
    evs[[length(evs) + 1]] <- ev(amt = ada_dose,  time = 28, cmt = "ADA_abs",
                                  ii = 14, addl = floor((day_max - 28) / 14))
  }
  if ("secukinumab" %in% drug) {
    for (t in c(0, 7, 14, 21, 28))
      evs[[length(evs) + 1]] <- ev(amt = sec_dose, time = t, cmt = "SEC_abs")
    evs[[length(evs) + 1]] <- ev(amt = sec_dose, time = 56, cmt = "SEC_abs",
                                  ii = 28, addl = floor((day_max - 56) / 28))
  }
  if ("bimekizumab" %in% drug) {
    evs[[length(evs) + 1]] <- ev(amt = bim_dose, time = 0, cmt = "BIM_abs",
                                  ii = 14, addl = floor(day_max / 14))
  }
  if (length(evs) == 0) return(ev(amt = 0, time = 0, cmt = "ADA_abs"))
  Reduce("+", evs)
}

## ---- UI --------------------------------------------------------------------
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "HS QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "tab1", icon = icon("user")),
      menuItem("Drug PK",             tabName = "tab2", icon = icon("pills")),
      menuItem("PD — Cytokines",      tabName = "tab3", icon = icon("flask")),
      menuItem("Clinical Endpoints",  tabName = "tab4", icon = icon("chart-line")),
      menuItem("Scenario Comparison", tabName = "tab5", icon = icon("sliders-h")),
      menuItem("Biomarker & VPop",    tabName = "tab6", icon = icon("dna"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box-header { background-color: #4a235a !important; color: white !important; }
      .small-box { border-radius: 8px; }
    "))),
    tabItems(

      ## ----------------------------------------------------------------
      ## TAB 1: Patient Profile
      ## ----------------------------------------------------------------
      tabItem(
        tabName = "tab1",
        fluidRow(
          box(
            title = "Patient Characteristics", status = "primary", solidHeader = TRUE,
            width = 4,
            selectInput("sex",    "Sex",
                        choices = c("Female" = 0, "Male" = 1), selected = 0),
            sliderInput("bw",     "Body Weight (kg)", min = 40, max = 150, value = 80, step = 5),
            sliderInput("age",    "Age (years)",      min = 18, max = 70,  value = 32, step = 1),
            selectInput("smoke",  "Smoking Status",
                        choices = c("Non-smoker" = 0, "Smoker" = 1), selected = 1),
            selectInput("hurley", "Hurley Stage",
                        choices = c("Stage I (mild)" = 1,
                                    "Stage II (moderate)" = 2,
                                    "Stage III (severe)" = 3),
                        selected = 2),
            sliderInput("an_base","Baseline AN Count", min = 1, max = 30, value = 5)
          ),
          box(
            title = "Disease Background", status = "warning", solidHeader = TRUE,
            width = 4,
            checkboxGroupInput("regions", "Affected Body Regions",
              choices = c("Axillary"        = "axilla",
                          "Inguinal/Groin"  = "groin",
                          "Inframammary"    = "infra",
                          "Buttocks/Perianal" = "butt"),
              selected = c("axilla", "groin")),
            sliderInput("fist_base", "Baseline Fistulae", min = 0, max = 10, value = 2),
            checkboxInput("ibd_comorbidity", "IBD Comorbidity", FALSE),
            checkboxInput("metabolic_synd",  "Metabolic Syndrome", FALSE),
            sliderInput("disease_duration",  "Disease Duration (years)", 1, 20, 5)
          ),
          box(
            title = "Prior Treatment History", status = "danger", solidHeader = TRUE,
            width = 4,
            checkboxGroupInput("prior_tx", "Prior Therapies",
              choices = c("Antibiotics (tetracycline class)"   = "abx",
                          "Adalimumab (failed)"                = "ada_fail",
                          "Hormonal therapy"                   = "horm",
                          "Surgical drainage"                   = "surg"),
              selected = c("abx")),
            sliderInput("prior_biologic_years", "Years on Prior Biologic", 0, 10, 0),
            helpText("Prior adalimumab failure may increase ADA antibody risk.")
          )
        ),
        fluidRow(
          valueBoxOutput("vbox_hurley"),
          valueBoxOutput("vbox_ihs4_0"),
          valueBoxOutput("vbox_dlqi_0")
        ),
        fluidRow(
          box(
            title = "Baseline Disease Characteristics Summary",
            status = "info", solidHeader = TRUE, width = 12,
            DT::dataTableOutput("baseline_table")
          )
        )
      ),

      ## ----------------------------------------------------------------
      ## TAB 2: Drug PK
      ## ----------------------------------------------------------------
      tabItem(
        tabName = "tab2",
        fluidRow(
          box(
            title = "Drug Selection & Dosing", status = "primary", solidHeader = TRUE,
            width = 3,
            checkboxGroupInput("drugs_pk", "Select Drug(s)",
              choices = c("Adalimumab"  = "adalimumab",
                          "Secukinumab" = "secukinumab",
                          "Bimekizumab" = "bimekizumab"),
              selected = "adalimumab"),
            conditionalPanel(
              "input.drugs_pk.includes('adalimumab')",
              sliderInput("ada_dose", "Adalimumab maintenance (mg Q2W)",
                          min = 20, max = 80, value = 40, step = 20),
              checkboxInput("ada_ada", "Anti-drug antibodies (+)", FALSE)
            ),
            conditionalPanel(
              "input.drugs_pk.includes('secukinumab')",
              selectInput("sec_dose", "Secukinumab dose",
                          choices = c("150 mg" = 150, "300 mg" = 300), selected = 300)
            ),
            conditionalPanel(
              "input.drugs_pk.includes('bimekizumab')",
              sliderInput("bim_dose", "Bimekizumab (mg Q2W)",
                          min = 160, max = 480, value = 320, step = 160)
            ),
            sliderInput("sim_weeks_pk", "Simulation Duration (weeks)", 4, 104, 52)
          ),
          box(
            title = "Drug Concentration–Time Profiles",
            status = "success", solidHeader = TRUE, width = 9,
            withSpinner(plotlyOutput("pk_plot", height = "450px"))
          )
        ),
        fluidRow(
          box(
            title = "PK Summary Table (Steady-state estimates)",
            status = "info", solidHeader = TRUE, width = 12,
            DT::dataTableOutput("pk_table")
          )
        )
      ),

      ## ----------------------------------------------------------------
      ## TAB 3: PD — Cytokine & Cell Dynamics
      ## ----------------------------------------------------------------
      tabItem(
        tabName = "tab3",
        fluidRow(
          box(
            title = "Cytokine Dynamics", status = "danger", solidHeader = TRUE,
            width = 6,
            withSpinner(plotlyOutput("cyt_plot", height = "400px"))
          ),
          box(
            title = "Inflammatory Cell Indices", status = "warning", solidHeader = TRUE,
            width = 6,
            withSpinner(plotlyOutput("cell_plot", height = "400px"))
          )
        ),
        fluidRow(
          box(
            title = "Drug Inhibition of Key Targets",
            status = "primary", solidHeader = TRUE, width = 6,
            withSpinner(plotlyOutput("inhib_plot", height = "300px"))
          ),
          box(
            title = "Disease Mediator Summary at Week 16 & 52",
            status = "info", solidHeader = TRUE, width = 6,
            DT::dataTableOutput("pd_table")
          )
        )
      ),

      ## ----------------------------------------------------------------
      ## TAB 4: Clinical Endpoints
      ## ----------------------------------------------------------------
      tabItem(
        tabName = "tab4",
        fluidRow(
          valueBoxOutput("vbox_hiscr16"),
          valueBoxOutput("vbox_ihs4_52"),
          valueBoxOutput("vbox_vas_52")
        ),
        fluidRow(
          box(
            title = "Abscess + Nodule (AN) Count Over Time",
            status = "danger", solidHeader = TRUE, width = 6,
            withSpinner(plotlyOutput("an_plot", height = "350px"))
          ),
          box(
            title = "IHS4 Score & HiSCR Response",
            status = "warning", solidHeader = TRUE, width = 6,
            withSpinner(plotlyOutput("ihs4_plot", height = "350px"))
          )
        ),
        fluidRow(
          box(
            title = "Quality of Life (DLQI) & VAS Pain",
            status = "success", solidHeader = TRUE, width = 6,
            withSpinner(plotlyOutput("qol_plot", height = "300px"))
          ),
          box(
            title = "Fistula/Sinus Tract Score",
            status = "primary", solidHeader = TRUE, width = 6,
            withSpinner(plotlyOutput("fist_plot", height = "300px"))
          )
        )
      ),

      ## ----------------------------------------------------------------
      ## TAB 5: Scenario Comparison
      ## ----------------------------------------------------------------
      tabItem(
        tabName = "tab5",
        fluidRow(
          box(
            title = "Comparator Scenarios", status = "primary", solidHeader = TRUE,
            width = 3,
            checkboxGroupInput("compare_scenarios", "Select Scenarios",
              choices = c(
                "No Treatment"          = "none",
                "Adalimumab 40 mg Q2W"  = "ada",
                "Secukinumab 300 mg"    = "sec",
                "Bimekizumab 320 mg"    = "bim",
                "ADA + SEC Combo"       = "ada_sec"
              ),
              selected = c("none", "ada", "bim")),
            sliderInput("comp_weeks", "Weeks to compare", 16, 104, 52)
          ),
          box(
            title = "AN Count Comparison", status = "success", solidHeader = TRUE, width = 9,
            withSpinner(plotlyOutput("comp_an", height = "350px"))
          )
        ),
        fluidRow(
          box(
            title = "IHS4 Comparison", status = "warning", solidHeader = TRUE, width = 6,
            withSpinner(plotlyOutput("comp_ihs4", height = "300px"))
          ),
          box(
            title = "HiSCR Response Rates by Scenario",
            status = "danger", solidHeader = TRUE, width = 6,
            withSpinner(plotlyOutput("comp_hiscr", height = "300px"))
          )
        ),
        fluidRow(
          box(
            title = "Comparative Efficacy Table",
            status = "info", solidHeader = TRUE, width = 12,
            DT::dataTableOutput("comp_table")
          )
        )
      ),

      ## ----------------------------------------------------------------
      ## TAB 6: Biomarker & Virtual Population
      ## ----------------------------------------------------------------
      tabItem(
        tabName = "tab6",
        fluidRow(
          box(
            title = "Biomarker Dashboard", status = "danger", solidHeader = TRUE,
            width = 6,
            withSpinner(plotlyOutput("bm_spider", height = "350px"))
          ),
          box(
            title = "Th17/Treg Ratio & IL-17A Over Time",
            status = "warning", solidHeader = TRUE, width = 6,
            withSpinner(plotlyOutput("th17_plot", height = "350px"))
          )
        ),
        fluidRow(
          box(
            title = "Virtual Population — HiSCR Distribution",
            status = "primary", solidHeader = TRUE,
            width = 6,
            sliderInput("vpop_n", "Number of virtual patients", 20, 200, 50, step = 10),
            sliderInput("vpop_bw_sd", "Body weight variability (SD, kg)", 5, 30, 15),
            withSpinner(plotlyOutput("vpop_hiscr", height = "300px"))
          ),
          box(
            title = "VPop — Smoking Impact on Outcomes",
            status = "success", solidHeader = TRUE, width = 6,
            withSpinner(plotlyOutput("vpop_smoke", height = "300px"))
          )
        )
      )
    )
  )
)

## ---- SERVER -----------------------------------------------------------------
server <- function(input, output, session) {

  ## Reactive: build patient-specific parameters
  pt_params <- reactive({
    list(
      BW      = as.numeric(input$bw),
      SEX     = as.numeric(input$sex),
      SMOKE   = as.numeric(input$smoke),
      HURLEY  = as.numeric(input$hurley),
      AN_base = as.numeric(input$an_base),
      ADA_ADA_effect = if (isTRUE(input$ada_ada)) 0.5 else 0.0
    )
  })

  ## Reactive: update model with patient params
  pt_mod <- reactive({
    p <- pt_params()
    m <- do.call(param, c(list(mod_global), p))
    init(m, AN = p$AN_base, Fist = input$fist_base)
  })

  ## Reactive: run simulation for single drug panel
  sim_pk <- reactive({
    drugs   <- input$drugs_pk
    ada_d   <- as.numeric(input$ada_dose)
    sec_d   <- as.numeric(input$sec_dose)
    bim_d   <- as.numeric(input$bim_dose)
    ada_ada <- isTRUE(input$ada_ada)
    weeks   <- input$sim_weeks_pk

    ev_tab  <- make_events(drugs, ada_d, sec_d, bim_d, weeks, ada_ada)
    m       <- pt_mod()
    if (ada_ada) m <- param(m, ADA_ADA_effect = 0.5)
    mrgsim(m, events = ev_tab, end = weeks * 7, delta = 1) %>% as.data.frame()
  })

  ## Reactive: scenario comparison
  sim_compare <- reactive({
    weeks <- input$comp_weeks
    scen  <- input$compare_scenarios
    m     <- pt_mod()

    scenario_evs <- list(
      none    = ev(amt = 0, time = 0, cmt = "ADA_abs"),
      ada     = make_events("adalimumab",  40,  300, 320, weeks, FALSE),
      sec     = make_events("secukinumab", 40,  300, 320, weeks, FALSE),
      bim     = make_events("bimekizumab", 40,  300, 320, weeks, FALSE),
      ada_sec = make_events(c("adalimumab","secukinumab"), 40, 300, 320, weeks, FALSE)
    )

    labels <- c(none="No Treatment", ada="Adalimumab 40Q2W",
                sec="Secukinumab 300", bim="Bimekizumab 320",
                ada_sec="ADA+SEC Combo")

    results <- lapply(scen, function(sc) {
      mrgsim(m, events = scenario_evs[[sc]], end = weeks * 7, delta = 1) %>%
        as.data.frame() %>%
        mutate(Scenario = labels[[sc]])
    })
    bind_rows(results)
  })

  ## ------- Tab 1 outputs -------
  output$vbox_hurley <- renderValueBox({
    valueBox(
      paste("Hurley", c("I","II","III")[as.numeric(input$hurley)]),
      "Disease Stage",
      icon = icon("layer-group"), color = "orange"
    )
  })
  output$vbox_ihs4_0 <- renderValueBox({
    ihs4 <- input$an_base * 1 + input$fist_base * 3
    valueBox(round(ihs4, 1), "Baseline IHS4", icon = icon("thermometer-half"), color = "red")
  })
  output$vbox_dlqi_0 <- renderValueBox({
    dlqi <- min(30, (input$an_base * 1 + input$fist_base * 3) * 0.8 + input$an_base * 0.8 * 0.5)
    valueBox(round(dlqi, 1), "Est. Baseline DLQI", icon = icon("smile-beam"), color = "purple")
  })
  output$baseline_table <- DT::renderDataTable({
    data.frame(
      Parameter = c("Age","Sex","Body Weight","Hurley Stage","AN Count",
                    "Fistulae","Smoking","Disease Duration","Comorbidities"),
      Value = c(
        input$age,
        ifelse(input$sex == 0, "Female", "Male"),
        paste(input$bw, "kg"),
        c("I","II","III")[as.numeric(input$hurley)],
        input$an_base,
        input$fist_base,
        ifelse(input$smoke == 1, "Yes", "No"),
        paste(input$disease_duration, "years"),
        paste(c(if(input$ibd_comorbidity) "IBD", if(input$metabolic_synd) "MetSyn"), collapse=", ")
      )
    )
  }, options = list(dom = "t", paging = FALSE))

  ## ------- Tab 2 outputs -------
  output$pk_plot <- renderPlotly({
    df <- sim_pk()
    df_long <- df %>%
      pivot_longer(c(C_ADA, C_SEC, C_BIM), names_to="Drug", values_to="Conc") %>%
      mutate(Drug = recode(Drug, C_ADA="Adalimumab", C_SEC="Secukinumab", C_BIM="Bimekizumab"))

    p <- ggplot(df_long, aes(x=time/7, y=Conc, color=Drug)) +
      geom_line(linewidth=1) +
      labs(x="Time (weeks)", y="Concentration (μg/mL)", color="Drug") +
      scale_color_manual(values=c(Adalimumab="#2980B9", Secukinumab="#27AE60", Bimekizumab="#8E44AD")) +
      theme_bw()
    ggplotly(p) %>% layout(legend=list(orientation="h"))
  })

  output$pk_table <- DT::renderDataTable({
    df <- sim_pk() %>% filter(time == max(time))
    data.frame(
      Drug        = c("Adalimumab","Secukinumab","Bimekizumab"),
      Trough_ugmL = round(c(df$C_ADA[1], df$C_SEC[1], df$C_BIM[1]), 2),
      Inhibition  = paste0(round(c(df$I_ADA_t[1]*100, df$I_SEC_t[1]*100, df$I_BIM_t[1]*100), 1),"%"),
      Target      = c("TNF-α","IL-17A","IL-17A/F")
    )
  }, options = list(dom="t", paging=FALSE))

  ## ------- Tab 3 outputs -------
  output$cyt_plot <- renderPlotly({
    df <- sim_pk() %>%
      pivot_longer(c(TNF,IL17A,IL6,IL1b,IL23), names_to="Cytokine", values_to="Val")
    p <- ggplot(df, aes(x=time/7, y=Val, color=Cytokine)) +
      geom_line(linewidth=0.9) +
      facet_wrap(~Cytokine, scales="free_y", nrow=2) +
      labs(x="Time (weeks)", y="Concentration") +
      theme_bw()
    ggplotly(p)
  })

  output$cell_plot <- renderPlotly({
    df <- sim_pk() %>%
      pivot_longer(c(Th17_idx, M1_idx, Neu_idx), names_to="Cell", values_to="Index") %>%
      mutate(Cell=recode(Cell, Th17_idx="Th17", M1_idx="M1 Macrophage", Neu_idx="Neutrophil"))
    p <- ggplot(df, aes(x=time/7, y=Index, color=Cell)) +
      geom_line(linewidth=1) +
      geom_hline(yintercept=1, linetype="dashed", color="gray60") +
      labs(x="Time (weeks)", y="Index (1=normal)") +
      theme_bw()
    ggplotly(p)
  })

  output$inhib_plot <- renderPlotly({
    df <- sim_pk() %>%
      pivot_longer(c(I_ADA_t, I_SEC_t, I_BIM_t), names_to="Drug", values_to="Inhibition") %>%
      mutate(Drug=recode(Drug, I_ADA_t="ADA→TNF", I_SEC_t="SEC→IL-17A", I_BIM_t="BIM→IL-17A/F"),
             Inhibition = Inhibition * 100)
    p <- ggplot(df, aes(x=time/7, y=Inhibition, color=Drug)) +
      geom_line(linewidth=1) +
      scale_y_continuous(limits=c(0,100)) +
      labs(x="Time (weeks)", y="Target Inhibition (%)") +
      theme_bw()
    ggplotly(p)
  })

  output$pd_table <- DT::renderDataTable({
    df <- sim_pk()
    wks <- c(16, 52) * 7
    df %>%
      filter(time %in% wks) %>%
      mutate(Week = paste0("Wk ", time/7)) %>%
      select(Week, TNF, IL17A, IL6, IL1b, Th17_idx, M1_idx, HiSCR) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
  }, options=list(dom="t", paging=FALSE))

  ## ------- Tab 4 outputs -------
  sim_with_noTx <- reactive({
    m  <- pt_mod()
    ev_none <- ev(amt=0, time=0, cmt="ADA_abs")
    weeks   <- input$sim_weeks_pk
    drugs   <- input$drugs_pk
    ada_d   <- as.numeric(input$ada_dose)
    sec_d   <- as.numeric(input$sec_dose)
    bim_d   <- as.numeric(input$bim_dose)
    ev_tx   <- make_events(drugs, ada_d, sec_d, bim_d, weeks, isTRUE(input$ada_ada))
    tx  <- mrgsim(m, events=ev_tx,   end=weeks*7, delta=1) %>% as.data.frame() %>% mutate(Arm="Treatment")
    ntx <- mrgsim(m, events=ev_none, end=weeks*7, delta=1) %>% as.data.frame() %>% mutate(Arm="No Treatment")
    bind_rows(tx, ntx)
  })

  output$vbox_hiscr16 <- renderValueBox({
    df <- sim_with_noTx() %>% filter(Arm=="Treatment", abs(time-112)<1)
    rate <- if (nrow(df) > 0) round(mean(df$HiSCR)*100,1) else 0
    valueBox(paste0(rate, "%"), "HiSCR @ Wk 16", icon=icon("check-circle"), color="green")
  })
  output$vbox_ihs4_52 <- renderValueBox({
    df <- sim_with_noTx() %>% filter(Arm=="Treatment", abs(time-364)<2)
    val <- if (nrow(df) > 0) round(mean(df$IHS4_cur),1) else 0
    valueBox(val, "IHS4 @ Wk 52", icon=icon("chart-bar"), color="orange")
  })
  output$vbox_vas_52 <- renderValueBox({
    df <- sim_with_noTx() %>% filter(Arm=="Treatment", abs(time-364)<2)
    val <- if (nrow(df) > 0) round(mean(df$VAS_pain),1) else 0
    valueBox(val, "VAS Pain @ Wk 52", icon=icon("heartbeat"), color="red")
  })

  output$an_plot <- renderPlotly({
    df <- sim_with_noTx()
    p <- ggplot(df, aes(x=time/7, y=AN, color=Arm)) +
      geom_line(linewidth=1.1) +
      geom_hline(yintercept=input$an_base*0.5, linetype="dashed", color="darkgreen") +
      annotate("text", x=4, y=input$an_base*0.5+0.3, label="HiSCR cutoff", size=3, color="darkgreen") +
      scale_color_manual(values=c(Treatment="#2980B9", "No Treatment"="#E74C3C")) +
      labs(x="Time (weeks)", y="AN Count", color="") + theme_bw()
    ggplotly(p)
  })
  output$ihs4_plot <- renderPlotly({
    df <- sim_with_noTx()
    p <- ggplot(df, aes(x=time/7, y=IHS4_cur, color=Arm)) +
      geom_line(linewidth=1.1) +
      scale_color_manual(values=c(Treatment="#2980B9", "No Treatment"="#E74C3C")) +
      labs(x="Time (weeks)", y="IHS4 Score", color="") + theme_bw()
    ggplotly(p)
  })
  output$qol_plot <- renderPlotly({
    df <- sim_with_noTx() %>%
      pivot_longer(c(DLQI, VAS_pain), names_to="Score", values_to="Val")
    p <- ggplot(df, aes(x=time/7, y=Val, color=Arm, linetype=Score)) +
      geom_line(linewidth=0.9) +
      facet_wrap(~Score, scales="free_y") +
      scale_color_manual(values=c(Treatment="#2980B9","No Treatment"="#E74C3C")) +
      labs(x="Time (weeks)", y="Score", color="") + theme_bw()
    ggplotly(p)
  })
  output$fist_plot <- renderPlotly({
    df <- sim_with_noTx()
    p <- ggplot(df, aes(x=time/7, y=Fist, color=Arm)) +
      geom_line(linewidth=1.1) +
      scale_color_manual(values=c(Treatment="#2980B9","No Treatment"="#E74C3C")) +
      labs(x="Time (weeks)", y="Fistula Score", color="") + theme_bw()
    ggplotly(p)
  })

  ## ------- Tab 5 outputs -------
  output$comp_an <- renderPlotly({
    df <- sim_compare()
    p <- ggplot(df, aes(x=time/7, y=AN, color=Scenario)) +
      geom_line(linewidth=1) +
      scale_color_brewer(palette="Set1") +
      labs(x="Weeks", y="AN Count", color="") + theme_bw()
    ggplotly(p)
  })
  output$comp_ihs4 <- renderPlotly({
    df <- sim_compare()
    p <- ggplot(df, aes(x=time/7, y=IHS4_cur, color=Scenario)) +
      geom_line(linewidth=1) +
      scale_color_brewer(palette="Set1") +
      labs(x="Weeks", y="IHS4", color="") + theme_bw()
    ggplotly(p)
  })
  output$comp_hiscr <- renderPlotly({
    df <- sim_compare() %>%
      filter(abs(time - input$comp_weeks*7) < 1) %>%
      group_by(Scenario) %>%
      summarise(HiSCR_rate=mean(HiSCR)*100, .groups="drop")
    p <- ggplot(df, aes(x=reorder(Scenario, -HiSCR_rate), y=HiSCR_rate, fill=Scenario)) +
      geom_col() +
      scale_fill_brewer(palette="Set1") +
      labs(x="", y=paste0("HiSCR Rate (%) at Wk ", input$comp_weeks), fill="") +
      theme_bw() + theme(axis.text.x=element_text(angle=30, hjust=1))
    ggplotly(p)
  })
  output$comp_table <- DT::renderDataTable({
    df <- sim_compare()
    wk16 <- df %>% filter(abs(time-112)<1) %>%
      group_by(Scenario) %>%
      summarise(
        HiSCR_Wk16 = paste0(round(mean(HiSCR)*100,1),"%"),
        AN_Wk16    = round(mean(AN),1),
        IHS4_Wk16  = round(mean(IHS4_cur),1),
        .groups    = "drop"
      )
    wk52 <- df %>% filter(abs(time-364)<2) %>%
      group_by(Scenario) %>%
      summarise(
        HiSCR_Wk52 = paste0(round(mean(HiSCR)*100,1),"%"),
        AN_Wk52    = round(mean(AN),1),
        IHS4_Wk52  = round(mean(IHS4_cur),1),
        DLQI_Wk52  = round(mean(DLQI),1),
        .groups    = "drop"
      )
    full_join(wk16, wk52, by="Scenario")
  }, options=list(dom="t", paging=FALSE))

  ## ------- Tab 6 outputs -------
  output$bm_spider <- renderPlotly({
    df <- sim_pk() %>% filter(abs(time-112)<1)
    bm <- data.frame(
      Biomarker = c("TNF-α", "IL-17A", "IL-6", "IL-1β", "Th17 idx", "M1 Macro", "Neutrophil"),
      Baseline  = c(3.571, 2.0, 8.0, 0.8, 1.5, 1.8, 2.0),
      Treatment = c(df$TNF[1], df$IL17A[1], df$IL6[1], df$IL1b[1],
                    df$Th17_idx[1], df$M1_idx[1], df$Neu_idx[1])
    ) %>%
      mutate(Pct_Baseline = Treatment/Baseline*100)
    plot_ly(type="scatterpolar", mode="lines+markers",
            r=c(bm$Pct_Baseline, bm$Pct_Baseline[1]),
            theta=c(bm$Biomarker, bm$Biomarker[1]),
            fill="toself", name="Wk 16 vs Baseline (%)") %>%
      layout(polar=list(radialaxis=list(range=c(0,120))))
  })

  output$th17_plot <- renderPlotly({
    df <- sim_pk() %>%
      mutate(Th17_Treg_ratio = Th17_idx / 0.7)
    p <- ggplot(df, aes(x=time/7)) +
      geom_line(aes(y=IL17A, color="IL-17A (pg/mL)"), linewidth=1) +
      geom_line(aes(y=Th17_Treg_ratio*0.5, color="Th17/Treg ratio"), linetype="dashed", linewidth=1) +
      scale_color_manual(values=c("IL-17A (pg/mL)"="#0000FF","Th17/Treg ratio"="#E74C3C")) +
      labs(x="Time (weeks)", y="Value", color="") + theme_bw()
    ggplotly(p)
  })

  output$vpop_hiscr <- renderPlotly({
    n   <- input$vpop_n
    bw_sd <- input$vpop_bw_sd
    set.seed(42)
    bws    <- rnorm(n, mean=as.numeric(input$bw), sd=bw_sd)
    smokes <- rbinom(n, 1, 0.5)
    drugs  <- input$drugs_pk
    ada_d  <- as.numeric(input$ada_dose)
    sec_d  <- as.numeric(input$sec_dose)
    bim_d  <- as.numeric(input$bim_dose)

    vp_hiscr <- sapply(seq_len(n), function(i) {
      m_vp <- param(mod_global, BW=pmax(40,bws[i]), SMOKE=smokes[i])
      ev_vp <- make_events(drugs, ada_d, sec_d, bim_d, 16, FALSE)
      out   <- mrgsim(m_vp, events=ev_vp, end=112, delta=7) %>% as.data.frame()
      tail(out$HiSCR, 1)
    })

    df_vp <- data.frame(HiSCR=factor(vp_hiscr, levels=c(0,1), labels=c("Non-responder","Responder")))
    p <- ggplot(df_vp, aes(x=HiSCR, fill=HiSCR)) +
      geom_bar(width=0.5) +
      scale_fill_manual(values=c(Responder="#27AE60", "Non-responder"="#E74C3C")) +
      labs(x="HiSCR at Week 16", y="Count", fill="") + theme_bw()
    ggplotly(p)
  })

  output$vpop_smoke <- renderPlotly({
    n <- min(input$vpop_n, 100)
    set.seed(99)
    bws <- rnorm(n, mean=as.numeric(input$bw), sd=15)
    drugs  <- input$drugs_pk
    ada_d  <- as.numeric(input$ada_dose)
    sec_d  <- as.numeric(input$sec_dose)
    bim_d  <- as.numeric(input$bim_dose)

    vp_rows <- lapply(c(0,1), function(sm) {
      sapply(seq_len(n), function(i) {
        m_vp <- param(mod_global, BW=pmax(40,bws[i]), SMOKE=sm)
        ev_vp <- make_events(drugs, ada_d, sec_d, bim_d, 52, FALSE)
        out   <- mrgsim(m_vp, events=ev_vp, end=364, delta=7) %>% as.data.frame()
        data.frame(IHS4_52 = tail(out$IHS4_cur, 1),
                   VAS_52  = tail(out$VAS_pain, 1),
                   Smoking = ifelse(sm==1,"Smoker","Non-smoker"))
      }) %>% t() %>% as.data.frame() %>%
        mutate(across(c(IHS4_52, VAS_52), as.numeric))
    })
    df_smoke <- bind_rows(vp_rows)
    p <- ggplot(df_smoke, aes(x=Smoking, y=IHS4_52, fill=Smoking)) +
      geom_boxplot(alpha=0.7) +
      scale_fill_manual(values=c(Smoker="#E74C3C","Non-smoker"="#3498DB")) +
      labs(x="Smoking Status", y="IHS4 at Week 52", fill="") + theme_bw()
    ggplotly(p)
  })
}

## ---- Launch app -------------------------------------------------------------
shinyApp(ui = ui, server = server)
