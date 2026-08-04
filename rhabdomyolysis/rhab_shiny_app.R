## =====================================================================
##  rhab_shiny_app.R
##  Rhabdomyolysis-Induced AKI — interactive QSP dashboard (10 tabs)
##  횡문근분해증 유발 급성 신손상 — 인터랙티브 QSP 대시보드
## ---------------------------------------------------------------------
##  The app is built around the model's central claim: the nephrotoxic
##  exposure is a PRODUCT of three terms and each therapy moves exactly
##  one of them.  The "Product decomposition" tab therefore shows the
##  three terms side by side rather than only their result, and the
##  "Therapy trade-offs" tab sweeps flow against pH so the diminishing
##  returns are visible instead of asserted.
##
##  Run with:  shiny::runApp("rhab_shiny_app.R")
##  Requires:  shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

## the model and the scenario library live in the model file
source("rhab_mrgsolve_model.R")

THEME <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#ECEFF1", colour = NA),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

## reference bands drawn on the marker panels
BANDS <- list(
  SCr  = c(0.6, 1.3), Kc = c(3.5, 5.0), CAc = c(1.15, 1.30),
  POc  = c(0.8, 1.45), HCO3 = c(22, 26), Clc = c(98, 107), Nac = c(135, 145)
)

## ---------------------------------------------------------------------
##  UI
## ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("횡문근분해증 유발 급성 신손상 QSP 모델 · Rhabdomyolysis-Induced AKI"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("Tubular toxic exposure &asymp; (filtered myoglobin) &divide; ",
              "(urine flow) &times; f<sub>diss</sub>(urine pH). ",
              "Crystalloid moves the denominator, alkali moves the chemistry, ",
              "extracorporeal therapy moves the numerator &mdash; and because ",
              "they multiply, each is worth less once the others succeed.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 · Patient"),
      sliderInput("minj", "위험 근육량 · muscle at risk (kg)",
                  min = 0, max = 14, value = 6, step = 0.5),
      sliderInput("age", "연령 · age (years)", min = 16, max = 90, value = 45),
      radioButtons("sex", "성별 · sex", inline = TRUE,
                   choices = c("남 male" = 0, "여 female" = 1), selected = 0),
      selectInput("cause", "원인 · cause",
                  choices = c("압사/구획 crush or compression" = "crush",
                              "운동성 exertional" = "exert",
                              "약물(스타틴 등) drug/statin" = "drug")),
      conditionalPanel(
        "input.cause == 'crush'",
        sliderInput("entrap", "매몰 시간 · entrapment (h)",
                    min = 0, max = 24, value = 8, step = 1),
        sliderInput("ccomp", "구획 유순도 · compartment compliance (L/mmHg)",
                    min = 0.006, max = 0.040, value = 0.028, step = 0.002),
        checkboxInput("fascio", "근막절개 · fasciotomy", FALSE),
        conditionalPanel("input.fascio",
                         sliderInput("fasct", "근막절개 시점 · at (h)",
                                     min = 0, max = 48, value = 10, step = 1))
      ),
      hr(),
      h4("수액 · Crystalloid  (the FLOW term)"),
      selectInput("fluid", "종류 · type",
                  choices = c("0.9% saline" = "NS", "Ringer lactate" = "LR",
                              "Balanced (Plasma-Lyte-like)" = "PL",
                              "0.45% saline + bicarbonate" = "BIC")),
      sliderInput("flrate", "속도 · rate (mL/h)",
                  min = 0, max = 1200, value = 500, step = 50),
      sliderInput("flon", "시작 시점 · started at (h)",
                  min = 0, max = 48, value = 8, step = 1),
      hr(),
      h4("알칼리화 · Alkali  (the pH term)"),
      sliderInput("bicr", HTML("NaHCO<sub>3</sub> (mmol/h)"),
                  min = 0, max = 40, value = 0, step = 2.5),
      sliderInput("aczr", "Acetazolamide (mg/h)",
                  min = 0, max = 60, value = 0, step = 5),
      hr(),
      h4("체외 순환 · Extracorporeal  (the SOURCE term)"),
      sliderInput("crrtq", "청소율 · clearance (mL/min)",
                  min = 0, max = 60, value = 0, step = 5),
      sliderInput("crrton", "시작 시점 · started at (h)",
                  min = 0, max = 72, value = 24, step = 1),
      hr(),
      h4("보조 요법 · Adjuncts"),
      sliderInput("manr", "Mannitol (g/h)", min = 0, max = 15, value = 0, step = 1),
      sliderInput("furr", "Furosemide (mg/h)", min = 0, max = 30, value = 0, step = 5),
      sliderInput("dantr", "Dantrolene (mg/h)", min = 0, max = 20, value = 0, step = 2),
      checkboxInput("insulin", "Insulin + dextrose (K⁺ shift)", FALSE),
      sliderInput("eaox", "항산화제 heme-arm block (fraction)",
                  min = 0, max = 0.95, value = 0, step = 0.05),
      hr(),
      sliderInput("tend", "관찰 기간 · follow-up (days)",
                  min = 3, max = 90, value = 28, step = 1),
      checkboxInput("cmp", "미치료 대조군 겹쳐 보기 · overlay untreated control", TRUE)
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs", type = "tabs",

        ## ---- 1 ------------------------------------------------------
        tabPanel(
          "1 · 환자 프로파일",
          h4("한눈에 보는 경과 · Course at a glance"),
          fluidRow(column(12, uiOutput("verdict"))),
          plotOutput("p_overview", height = "560px"),
          helpText("붉은 띠는 매몰 기간, 파란 띠는 수액 투여 기간입니다. ",
                   "Red band = entrapment, blue band = crystalloid infusion.")
        ),

        ## ---- 2 ------------------------------------------------------
        tabPanel(
          "2 · 두 개의 시계 (CK vs Mb)",
          h4("One release flux, two half-lives"),
          plotOutput("p_clocks", height = "420px"),
          plotOutput("p_ratio", height = "260px"),
          helpText(HTML("CK (86 kDa, not filtered, t&frac12; ~36 h) integrates the ",
                        "insult; myoglobin (17.8 kDa, filtered, t&frac12; ~2.4 h) ",
                        "tracks its instantaneous rate. The lag between the two ",
                        "peaks is why a CK threshold discriminates poorly, and the ",
                        "CK:myoglobin ratio is a readout of renal myoglobin ",
                        "clearance."))
        ),

        ## ---- 3 ------------------------------------------------------
        tabPanel(
          "3 · 곱의 분해 (product decomposition)",
          h4("The three terms that multiply"),
          plotOutput("p_product", height = "520px"),
          helpText(HTML("Top: filtered load (the numerator). Middle: urine flow ",
                        "(the denominator). Bottom: f<sub>diss</sub>(urine pH) ",
                        "(the chemistry). The resulting instantaneous toxic ",
                        "exposure is the fourth panel. Moving any single term ",
                        "leaves the others in charge."))
        ),

        ## ---- 4 ------------------------------------------------------
        tabPanel(
          "4 · 세뇨관 처리 (Tm 포화)",
          h4("A saturable transporter splits one load into two arms"),
          plotOutput("p_tubular", height = "440px"),
          plotOutput("p_arms", height = "260px"),
          helpText("Below Tm almost everything is endocytosed (arm A, heme) and ",
                   "the urine dipstick can be negative; above Tm the excess ",
                   "spills distally (arm B, casts). Because the per-nephron ",
                   "burden is divided by surviving nephron mass, injury per ",
                   "nephron accelerates as nephrons are lost.")
        ),

        ## ---- 5 ------------------------------------------------------
        tabPanel(
          "5 · 압력 노드 (ΔP)",
          h4("One pressure variable, two opposite gates"),
          plotOutput("p_press", height = "440px"),
          plotOutput("p_gates", height = "240px"),
          helpText(HTML("Necrosis needs ischaemia; washout needs perfusion. ",
                        "While the limb is entrapped the first gate is open and ",
                        "the second is shut &mdash; muscle dies silently. ",
                        "Extrication (or fasciotomy) flips both, which is the ",
                        "origin of the post-extrication potassium surge and of ",
                        "the limb-versus-kidney trade-off."))
        ),

        ## ---- 6 ------------------------------------------------------
        tabPanel(
          "6 · 전해질 · 산염기",
          h4("Potassium feeds back; calcium is biphasic"),
          plotOutput("p_elec", height = "560px"),
          helpText("Renal potassium excretion is written as urine flow times a ",
                   "regulated urine concentration, so anuria stops it ",
                   "arithmetically. Calcium is consumed into an injured-muscle ",
                   "deposit early and released when macrophages clear the ",
                   "scaffold — the rebound hypercalcaemia of the recovery phase.")
        ),

        ## ---- 7 ------------------------------------------------------
        tabPanel(
          "7 · 신장 결과 · 임상 엔드포인트",
          h4("Nephron mass, fibrosis and staging"),
          plotOutput("p_renal", height = "440px"),
          h4("계산된 엔드포인트 · computed endpoints"),
          DT::dataTableOutput("t_endpoints")
        ),

        ## ---- 8 ------------------------------------------------------
        tabPanel(
          "8 · 치료 상충 관계 (flow × pH)",
          h4("Diminishing returns, computed not asserted"),
          actionButton("run_sweep", "스윕 실행 · run the 4 × 3 sweep",
                       class = "btn-primary"),
          helpText("Twelve simulations: urine-driving fluid rate against three ",
                   "alkali doses, at the patient set on the left."),
          plotOutput("p_sweep", height = "380px"),
          DT::dataTableOutput("t_sweep")
        ),

        ## ---- 9 ------------------------------------------------------
        tabPanel(
          "9 · 시나리오 비교 (18개)",
          h4("The scenario library"),
          checkboxGroupInput("scen_pick", NULL, inline = TRUE,
                             choices = names(scenarios),
                             selected = names(scenarios)[c(2, 3, 4, 5, 6)]),
          actionButton("run_scen", "선택 시나리오 실행 · run selection",
                       class = "btn-primary"),
          plotOutput("p_scen", height = "460px"),
          DT::dataTableOutput("t_scen")
        ),

        ## ---- 10 -----------------------------------------------------
        tabPanel(
          "10 · 바이오마커 · 모니터링",
          h4("What you can actually measure, and when"),
          plotOutput("p_bio", height = "480px"),
          helpText("Urine pH and the osmolar gap are the two bedside variables ",
                   "that decide whether the alkali and the mannitol arms of the ",
                   "model are doing what was intended. Dipstick myoglobinuria ",
                   "turns positive only once the reabsorptive maximum is ",
                   "exceeded.")
        ),

        ## ---- about --------------------------------------------------
        tabPanel(
          "ℹ️ 모델 정보",
          h4("Model provenance"),
          tags$ul(
            tags$li("47 ODE compartments; 187 parameters; 18 pre-built scenarios."),
            tags$li(HTML("Every ODE was independently re-implemented in ",
                         "Python/scipy for calibration, and the ",
                         "<code>$ODE</code> block of the model file was then ",
                         "extracted, compiled and integrated separately to ",
                         "confirm the port. Eight defects found that way are ",
                         "documented in the model file.")),
            tags$li(HTML("The healthy control (scenario 1) is an exact steady ",
                         "state: twelve monitored concentrations drift &lt;3% ",
                         "over 14 simulated days.")),
            tags$li(HTML("KDIGO stage and the McMahon score are computed as ",
                         "<i>outputs</i> in <code>$TABLE</code>, not supplied ",
                         "as inputs."))
          ),
          h4("한계 · Limitations"),
          tags$ul(
            tags$li("단일 구획 근육 모델이며 구획압은 하나로 통합되어 있습니다."),
            tags$li("패혈증·창상 감염·절단이 모델에 없으므로 '근막절개 미시행' ",
                    "시나리오를 임상 권고로 읽어서는 안 됩니다."),
            tags$li("사망은 보정된 생존 모델이 아니라 누적 부정맥 위험도로만 ",
                    "표현됩니다."),
            tags$li("교육 및 가설 생성 목적이며 임상 의사결정 도구가 아닙니다.")
          )
        )
      )
    )
  )
)

## ---------------------------------------------------------------------
##  server
## ---------------------------------------------------------------------
server <- function(input, output, session) {

  ## ---- assemble the parameter set from the sidebar ----------------
  pset <- reactive({
    p <- as.list(FLUIDS[[input$fluid]])
    crush <- input$cause == "crush"
    ent <- if (crush) input$entrap else 0
    p$MINJ0   <- input$minj
    p$AGE     <- input$age
    p$SEXF    <- as.numeric(input$sex)
    p$CRUSH   <- as.numeric(crush)
    p$PEXT    <- if (crush) 150 else 0
    p$PEXTOFF <- ent
    p$ORALON  <- ent
    p$CCOMP   <- if (crush) input$ccomp else 0.028
    p$FASCT   <- if (crush && isTRUE(input$fascio)) input$fasct else 1e6
    ## non-mechanical causes enter through the primary-necrosis door
    if (input$cause == "exert") {
      p$PRIMR <- input$minj/20; p$PRIMON <- 0; p$PRIMOFF <- 20
    } else if (input$cause == "drug") {
      p$PRIMR <- input$minj/75; p$PRIMON <- 0; p$PRIMOFF <- 72
    }
    ## crystalloid: full rate to 48 h, half to 120 h, then stop
    p$FLRATE  <- input$flrate/1000
    p$FLON    <- input$flon
    p$FLOFF   <- 48
    p$FLRATE2 <- input$flrate/2000
    p$FLOFF2  <- 120
    ## alkali
    p$BICR <- input$bicr; p$BICON <- input$flon; p$BICOFF <- 72
    p$ACZR <- input$aczr; p$ACZON <- input$flon; p$ACZOFF <- 96
    ## extracorporeal
    p$CRRTQ <- input$crrtq*60/1000; p$CRRTON <- input$crrton; p$CRRTOFF <- 336
    ## adjuncts
    p$MANR <- input$manr; p$MANON <- input$flon; p$MANOFF <- 72
    p$FURR <- input$furr; p$FURON <- input$flon; p$FUROFF <- 72
    p$DANTR <- input$dantr; p$DANTON <- 1; p$DANTOFF <- 24
    if (isTRUE(input$insulin)) { p$INSR <- 1; p$INSON <- ent; p$INSOFF <- ent + 12 }
    p$EAOX <- input$eaox
    if (input$eaox > 0) { p$AOXR <- 1; p$AOXON <- input$flon; p$AOXOFF <- 168 }
    p
  })

  sim <- reactive({
    rhab %>% param(pset()) %>%
      mrgsim(end = input$tend*24, delta = 0.25, hmax = 0.5) %>%
      as_tibble() %>% mutate(arm = "선택 요법 · your regimen")
  })

  ## the same patient with no therapy at all
  ctrl <- reactive({
    p <- pset()
    for (k in c("FLRATE", "FLRATE2", "BICR", "ACZR", "CRRTQ", "MANR", "FURR",
                "DANTR", "INSR", "AOXR", "EAOX")) p[[k]] <- 0
    p$FASCT <- 1e6
    rhab %>% param(p) %>%
      mrgsim(end = input$tend*24, delta = 0.25, hmax = 0.5) %>%
      as_tibble() %>% mutate(arm = "미치료 대조군 · untreated")
  })

  both <- reactive(if (isTRUE(input$cmp)) bind_rows(sim(), ctrl()) else sim())

  ## ---- derived quantities the model does not capture directly -----
  derive <- function(d) {
    pp <- as.list(param(rhab))
    ps <- pset()
    for (n in names(ps)) pp[[n]] <- ps[[n]]
    d %>% mutate(
      day   = time/24,
      VS    = pmax(ECFV, 6) + pmax(MSEQ, 0),
      pHu   = { fa <- pmin(1, (pmax(HCO3 - pp$PHUOFF, 0)^2 /
                                (pmax(HCO3 - pp$PHUOFF, 0)^2 + pp$PHUK^2)) *
                             (1/(1 + 6*pmax(0, 1 - pmax(ECFV, 6)/pp$ECFV0))) *
                             pmin(1, Kc/3.6))
                az <- pp$EACZ*(ACZ/(ACZ + 120))
                fa <- fa + az*(1 - fa)
                pp$PHUMIN + (pp$PHUMAX - pp$PHUMIN)*pmin(fa, 1) },
      fdiss = 1/(1 + 10^((pHu - pp$PHD)/pp$PHDS)),
      GFRe  = pmax(pp$GFRMAX*NEPH*1, 0.004),
      FL    = GFRe*MBc*pp$SIEV,
      Tm    = pp$TMMB*pmax(NEPH, 0.02),
      Reab  = FL/(1 + FL/pmax(Tm, 0.5)),
      Exc   = FL - Reab,
      MAP   = pp$MAP0*(0.42 + 0.58*(pmax(ECFV, 6)/pp$ECFV0)^1.3),
      Pcomp = pp$PCOMP0 + pp$FCOMP*pmax(MSEQ, 0)/
              (pp$CCOMP*ifelse(time >= pp$FASCT, pp$FASCF, 1)) +
              ifelse(time < pp$PEXTOFF, pp$PEXT, 0),
      dP    = MAP - Pcomp,
      gate1 = 1/(1 + exp(pmax(pmin((dP - pp$DPCRIT)/6, 50), -50))),
      gate2 = 1 - gate1,
      osmgap = (MANC/pp$VMAN)*1000/182.17
    )
  }

  dd <- reactive(derive(both()))

  ## ---- 1 : overview ----------------------------------------------
  output$verdict <- renderUI({
    d <- dd() %>% filter(arm == "선택 요법 · your regimen")
    kd <- max(d$KDIGO); scr <- max(d$SCr); k <- max(d$Kc)
    eg <- tail(d$eGFR, 1); ol <- max(d$ANURH); mcm <- tail(d$MCMAHON, 1)
    col <- c("#2E7D32", "#F9A825", "#EF6C00", "#C62828")[kd + 1]
    div(style = paste0("border-left:6px solid ", col,
                       ";background:#FAFAFA;padding:10px 14px;margin-bottom:10px"),
        HTML(sprintf(
          "<b style='color:%s'>KDIGO %d</b> &nbsp;·&nbsp; peak creatinine <b>%.2f mg/dL</b>
           &nbsp;·&nbsp; peak K<sup>+</sup> <b>%.2f mmol/L</b>
           &nbsp;·&nbsp; oligo-anuric <b>%.0f h</b>
           &nbsp;·&nbsp; McMahon <b>%.1f</b>
           &nbsp;·&nbsp; eGFR at day %d <b>%.0f mL/min</b>
           &nbsp;·&nbsp; peak CK <b>%s U/L</b>",
          col, kd, scr, k, ol, mcm, input$tend, eg,
          format(round(max(d$CK)), big.mark = ","))))
  })

  output$p_overview <- renderPlot({
    d <- dd() %>%
      select(day, arm, CK, MBc, SCr, Kc, CAc, HCO3, eGFR, NEPH) %>%
      pivot_longer(-c(day, arm))
    lab <- c(CK = "CK (U/L)", MBc = "myoglobin (mg/L)", SCr = "creatinine (mg/dL)",
             Kc = "K+ (mmol/L)", CAc = "ionised Ca (mmol/L)",
             HCO3 = "HCO3- (mmol/L)", eGFR = "eGFR (mL/min)",
             NEPH = "nephron mass (fraction)")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(day, value, colour = arm)) +
      annotate("rect", xmin = 0, xmax = if (input$cause == "crush") input$entrap/24 else 0,
               ymin = -Inf, ymax = Inf, fill = "#B71C1C", alpha = 0.10) +
      annotate("rect", xmin = input$flon/24, xmax = 2,
               ymin = -Inf, ymax = Inf, fill = "#0D47A1", alpha = 0.06) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      scale_colour_manual(values = c("#C62828", "#1565C0"), name = NULL) +
      labs(x = "day", y = NULL) + THEME
  })

  ## ---- 2 : two clocks --------------------------------------------
  output$p_clocks <- renderPlot({
    d <- dd() %>% filter(arm == "선택 요법 · your regimen", time <= 168)
    pk <- d %>% summarise(tck = time[which.max(CK)], tmb = time[which.max(MBc)])
    ggplot(d, aes(time)) +
      geom_line(aes(y = CK/1000, colour = "CK (kU/L) — integrator"), linewidth = 1) +
      geom_line(aes(y = MBc, colour = "myoglobin (mg/L) — rate meter"), linewidth = 1) +
      geom_vline(xintercept = pk$tmb, linetype = 2, colour = "#1B5E20") +
      geom_vline(xintercept = pk$tck, linetype = 2, colour = "#B71C1C") +
      annotate("segment", x = pk$tmb, xend = pk$tck, y = 0, yend = 0,
               arrow = arrow(ends = "both", length = unit(2, "mm"))) +
      annotate("label", x = (pk$tmb + pk$tck)/2, y = 0,
               label = sprintf("lag %.1f h", pk$tck - pk$tmb), size = 3.4) +
      scale_colour_manual(values = c("#B71C1C", "#1B5E20"), name = NULL) +
      labs(x = "hours", y = NULL, title = "Same release flux, two half-lives") + THEME
  })

  output$p_ratio <- renderPlot({
    d <- dd() %>% filter(time <= 336, time > 6)
    ggplot(d, aes(time, CKMB, colour = arm)) + geom_line(linewidth = 0.8) +
      scale_y_log10() +
      scale_colour_manual(values = c("#C62828", "#1565C0"), name = NULL) +
      labs(x = "hours", y = "CK : myoglobin ratio",
           subtitle = "the ratio rises as renal myoglobin clearance is lost") + THEME
  })

  ## ---- 3 : product decomposition ---------------------------------
  output$p_product <- renderPlot({
    d <- dd() %>% filter(arm == "선택 요법 · your regimen", time <= 168) %>%
      transmute(time,
                `1 filtered load (mg/h)  — numerator` = FL,
                `2 urine flow (mL/h)  — denominator` = 1000*(FL*0 +
                   pmax(c(0, diff(UOUT))/0.25*1000, 0)/1000),
                `3 f_diss(urine pH)  — chemistry` = fdiss,
                `4 instantaneous toxic exposure` = c(0, diff(TOXAUC))/0.25) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(colour = "#4527A0", linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 1) +
      labs(x = "hours", y = NULL) + THEME
  })

  ## ---- 4 : tubular handling --------------------------------------
  output$p_tubular <- renderPlot({
    d <- dd() %>% filter(arm == "선택 요법 · your regimen", time <= 168)
    ggplot(d, aes(time)) +
      geom_ribbon(aes(ymin = 0, ymax = Reab, fill = "reabsorbed (arm A: heme)"),
                  alpha = 0.55) +
      geom_ribbon(aes(ymin = Reab, ymax = FL, fill = "spilled distally (arm B: casts)"),
                  alpha = 0.55) +
      geom_line(aes(y = Tm, linetype = "Tm (falls with nephron mass)"),
                linewidth = 0.8) +
      scale_fill_manual(values = c("#F57F17", "#3E2723"), name = NULL) +
      scale_linetype_manual(values = 2, name = NULL) +
      labs(x = "hours", y = "mg/h",
           title = "One filtered load, split by a saturable transporter") + THEME
  })

  output$p_arms <- renderPlot({
    d <- dd() %>% filter(arm == "선택 요법 · your regimen", time <= 336)
    ggplot(d, aes(time)) +
      geom_line(aes(y = CAST, colour = "cast burden"), linewidth = 0.9) +
      geom_line(aes(y = TUBI, colour = "tubular injury"), linewidth = 0.9) +
      geom_line(aes(y = NEPH, colour = "nephron mass"), linewidth = 0.9) +
      geom_line(aes(y = FIB, colour = "fibrosis"), linewidth = 0.9) +
      scale_colour_manual(values = c("#3E2723", "#3949AB", "#B71C1C", "#00695C"),
                          name = NULL) +
      labs(x = "hours", y = "fraction") + THEME
  })

  ## ---- 5 : pressures ---------------------------------------------
  output$p_press <- renderPlot({
    d <- dd() %>% filter(arm == "선택 요법 · your regimen", time <= 240)
    ggplot(d, aes(time)) +
      geom_hline(yintercept = 30, linetype = 3, colour = "#B71C1C") +
      geom_line(aes(y = MAP, colour = "MAP"), linewidth = 0.9) +
      geom_line(aes(y = Pcomp, colour = "compartment pressure"), linewidth = 0.9) +
      geom_line(aes(y = dP, colour = "delta P = MAP - Pcomp"), linewidth = 1.1) +
      annotate("text", x = max(d$time)*0.85, y = 34, size = 3.2,
               label = "fasciotomy threshold 30 mmHg", colour = "#B71C1C") +
      scale_colour_manual(values = c("#4527A0", "#7E57C2", "#B71C1C"), name = NULL) +
      labs(x = "hours", y = "mmHg") + THEME
  })

  output$p_gates <- renderPlot({
    d <- dd() %>% filter(arm == "선택 요법 · your regimen", time <= 240)
    ggplot(d, aes(time)) +
      geom_line(aes(y = gate1, colour = "gate 1 — ischaemia (makes necrosis)"),
                linewidth = 1) +
      geom_line(aes(y = gate2, colour = "gate 2 — perfusion (makes washout)"),
                linewidth = 1) +
      scale_colour_manual(values = c("#311B92", "#00838F"), name = NULL) +
      labs(x = "hours", y = "gate (0-1)",
           subtitle = "the same ΔP opens one and closes the other") + THEME
  })

  ## ---- 6 : electrolytes ------------------------------------------
  output$p_elec <- renderPlot({
    d <- dd() %>%
      select(day, arm, Kc, CAc, POc, HCO3, Clc, Nac, URc, CADEP) %>%
      pivot_longer(-c(day, arm))
    lab <- c(Kc = "K+ (mmol/L)", CAc = "ionised Ca (mmol/L)",
             POc = "phosphate (mmol/L)", HCO3 = "HCO3- (mmol/L)",
             Clc = "Cl- (mmol/L)", Nac = "Na+ (mmol/L)",
             URc = "urate (mg/dL)", CADEP = "Ca-PO4 deposit (mmol)")
    d$name <- factor(lab[d$name], levels = lab)
    ggplot(d, aes(day, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 4) +
      scale_colour_manual(values = c("#C62828", "#1565C0"), name = NULL) +
      labs(x = "day", y = NULL) + THEME
  })

  ## ---- 7 : renal outcome -----------------------------------------
  output$p_renal <- renderPlot({
    d <- dd() %>% select(day, arm, NEPH, FIB, eGFR, SCr, BUNc, KDIGO) %>%
      pivot_longer(-c(day, arm))
    ggplot(d, aes(day, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = c("#C62828", "#1565C0"), name = NULL) +
      labs(x = "day", y = NULL) + THEME
  })

  output$t_endpoints <- DT::renderDataTable({
    dd() %>% group_by(arm) %>% summarise(
      `peak CK (U/L)` = round(max(CK)),
      `time to CK peak (h)` = round(time[which.max(CK)], 1),
      `peak myoglobin (mg/L)` = round(max(MBc), 1),
      `peak creatinine (mg/dL)` = round(max(SCr), 2),
      `KDIGO` = max(KDIGO),
      `McMahon` = round(tail(MCMAHON, 1), 1),
      `peak K+ (mmol/L)` = round(max(Kc), 2),
      `hours K+ > 6` = round(max(HKH), 1),
      `ionised Ca nadir` = round(min(CAc), 2),
      `oligo-anuric hours` = round(max(ANURH)),
      `urine 0-24 h (mL)` = round(approx(time, UOUT, 24)$y*1000),
      `myoglobin filtered (g)` = round(max(MBFILT)/1000, 2),
      `muscle lysed (kg)` = round(max(MLYS), 2),
      `third space peak (L)` = round(max(MSEQ), 1),
      `fibrosis (fraction)` = round(tail(FIB, 1), 3),
      `final eGFR (mL/min)` = round(tail(eGFR, 1), 1),
      .groups = "drop") %>%
      t() %>% as.data.frame() %>% setNames(NULL) %>%
      DT::datatable(options = list(dom = "t", pageLength = 20), rownames = TRUE)
  })

  ## ---- 8 : flow x pH sweep ---------------------------------------
  sweep <- eventReactive(input$run_sweep, {
    base <- pset()
    grid <- expand.grid(flow = c(100, 250, 500, 1000), bic = c(0, 12.5, 30))
    withProgress(message = "12 simulations…", {
      bind_rows(lapply(seq_len(nrow(grid)), function(i) {
        incProgress(1/nrow(grid))
        p <- base
        p$FLRATE <- grid$flow[i]/1000; p$FLRATE2 <- grid$flow[i]/2000
        p$FLOFF <- 48; p$FLOFF2 <- 120
        p$BICR <- grid$bic[i]; p$BICOFF <- 72
        rhab %>% param(p) %>% mrgsim(end = 1440, delta = 1, hmax = 0.5) %>%
          as_tibble() %>%
          summarise(flow = grid$flow[i], bic = grid$bic[i],
                    TOXAUC = max(TOXAUC), SCr = max(SCr),
                    KDIGO = max(KDIGO), eGFR = tail(eGFR, 1),
                    pHu24 = NA_real_)
      }))
    })
  })

  output$p_sweep <- renderPlot({
    s <- sweep()
    ggplot(s, aes(factor(flow), TOXAUC, fill = factor(bic))) +
      geom_col(position = position_dodge(0.8), width = 0.7) +
      geom_text(aes(label = round(TOXAUC)), position = position_dodge(0.8),
                vjust = -0.4, size = 3) +
      scale_fill_manual(values = c("#B0BEC5", "#66BB6A", "#1B5E20"),
                        name = "NaHCO3 (mmol/h)") +
      labs(x = "crystalloid rate (mL/h)", y = "cumulative toxic exposure",
           title = "The same alkali dose buys less as the flow term succeeds") +
      THEME
  })

  output$t_sweep <- DT::renderDataTable({
    s <- sweep() %>% select(flow, bic, TOXAUC, SCr, KDIGO, eGFR) %>%
      mutate(across(c(TOXAUC, SCr, eGFR), ~round(.x, 2)))
    DT::datatable(s, options = list(dom = "tp", pageLength = 12), rownames = FALSE)
  })

  ## ---- 9 : scenario library --------------------------------------
  scen <- eventReactive(input$run_scen, {
    picks <- input$scen_pick
    validate(need(length(picks) > 0, "시나리오를 선택하세요."))
    withProgress(message = "running scenarios…", {
      bind_rows(lapply(picks, function(nm) {
        incProgress(1/length(picks))
        run_scenario(nm, end = 1440, delta = 1)
      }))
    })
  })

  output$p_scen <- renderPlot({
    d <- scen() %>% mutate(day = time/24) %>%
      select(day, scenario, CK, MBc, SCr, Kc, eGFR, CAc) %>%
      pivot_longer(-c(day, scenario))
    ggplot(d, aes(day, value, colour = scenario)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "day", y = NULL) + THEME +
      theme(legend.text = element_text(size = 8))
  })

  output$t_scen <- DT::renderDataTable({
    s <- summarise_scenario(scen()) %>% mutate(across(where(is.numeric), ~round(.x, 2)))
    DT::datatable(s, options = list(dom = "tp", scrollX = TRUE, pageLength = 20),
                  rownames = FALSE)
  })

  ## ---- 10 : biomarkers / monitoring ------------------------------
  output$p_bio <- renderPlot({
    d <- dd() %>% filter(time <= 336) %>%
      transmute(time, arm,
                `urine pH` = pHu,
                `plasma pH` = pHpl,
                `distal myoglobin spill (mg/h)` = Exc,
                `osmolar gap from mannitol (mOsm)` = osmgap,
                `CK : myoglobin ratio` = CKMB,
                `cumulative arrhythmia hazard` = ARRH) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = c("#C62828", "#1565C0"), name = NULL) +
      labs(x = "hours", y = NULL) + THEME
  })
}

shinyApp(ui, server)
