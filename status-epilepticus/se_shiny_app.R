## =====================================================================
##  se_shiny_app.R
##  Status Epilepticus QSP model — interactive dashboard
##  경련지속상태 QSP 모델 — Shiny 대시보드
##
##  10 tabs, each answering one question the model exists to answer:
##
##   1  환자·병인 (Patient & aetiology)   — set the DRIVE, see where the
##                                          thresholds land for this patient
##   2  치료 타임라인 (Treatment timeline) — build a drug schedule minute by
##                                          minute and watch the clock run
##   3  수용체 시계 (The receptor clock)   — R_SYN out, NR_SYN in, R_EXTRA flat
##   4  약물 PK (Pharmacokinetics)         — 11 drugs, plasma and effect site
##   5  1분의 가격 (The price of a minute) — CREQ_BZD(t) and its asymptote
##   6  교차점 (The crossover)             — benzodiazepine vs ketamine target
##   7  운동-뇌파 해리 (Motor vs EEG)      — what the bedside sees vs the brain
##   8  전신 생리 (Systemic physiology)    — the phase I / phase II switch
##   9  뇌 손상 (Injury & outcome)         — the second integrator
##  10  가상 집단 (Virtual population)     — response RATES, ESETT / RAMPART
##
##  run:  shiny::runApp("se_shiny_app.R")
##  (loads se_mrgsolve_model.R from the same directory)
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ---------------------------------------------------------------------
#  Load the model without executing the script's own reporting section
# ---------------------------------------------------------------------
.src <- readLines("se_mrgsolve_model.R")
.cut <- grep("^#  MAIN$", .src)
if (length(.cut)) .src <- .src[1:(.cut[1] - 2)]
eval(parse(text = paste(.src, collapse = "\n")))

theme_se <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom")

PALETTE <- c("#1E88E5", "#E65100", "#1B5E20", "#880E4F", "#4527A0",
             "#00838F", "#B71C1C", "#546E7A", "#F9A825", "#6A1B9A")

# =====================================================================
#  UI
# =====================================================================
ui <- fluidPage(
  titlePanel("경련지속상태 QSP 모델 · Status Epilepticus — a receptor-trafficking clock"),
  tags$p(style = "color:#555;margin-top:-8px",
         HTML("<b>59 ODEs · 11 drugs · 22 scenarios.</b> One clock drives two receptor pools in opposite directions and leaves a third alone. Every curve below is that arithmetic.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("① 환자 · 병인"),
      selectInput("aetiology", "병인 (aetiology)",
                  choices = c("AED 중단 (withdrawal)"          = 1.15,
                              "대사성 (metabolic)"             = 1.35,
                              "전형적 (typical adult SE)"      = 1.55,
                              "급성 뇌졸중 (acute stroke)"     = 1.80,
                              "CNS 감염 (CNS infection)"       = 2.10,
                              "자가면역 뇌염 (anti-NMDAR)"     = 2.40,
                              "FIRES / NORSE"                  = 2.60,
                              "저산소성 (post-anoxic)"         = 3.20),
                  selected = 1.55),
      checkboxInput("reversible", "병인이 치료 가능 (drive removed at 30 min)", FALSE),
      sliderInput("alb", "혈청 알부민 (g/dL)", 1.8, 5.0, 4.0, 0.1),
      sliderInput("crcl", "크레아티닌 청소율 (mL/min)", 10, 140, 90, 5),
      checkboxInput("anak", "IL-1 차단 (anakinra)", FALSE),
      checkboxInput("gluciv", "포도당 + 티아민 정주", FALSE),
      checkboxInput("cool", "목표체온조절 (TTM)", FALSE),

      hr(),
      h4("② 1차 약제 (first line)"),
      selectInput("bzd", "benzodiazepine",
                  c("없음 (none)" = "none",
                    "lorazepam 4 mg IV"        = "lzp4",
                    "lorazepam 2 mg IV (저용량)" = "lzp2",
                    "midazolam 10 mg IM"       = "mdzim",
                    "midazolam 0.2 mg/kg IV"   = "mdziv",
                    "diazepam 10 mg IV"        = "dzp"),
                  selected = "lzp4"),
      sliderInput("t_bzd", "투여 시각 (min)", 0, 120, 8, 1),
      checkboxInput("bzd_rep", "4분 뒤 2회차 투여", TRUE),

      hr(),
      h4("③ 2차 약제 (second line)"),
      selectInput("second", "agent",
                  c("없음 (none)" = "none",
                    "levetiracetam 60 mg/kg" = "lev",
                    "fosphenytoin 20 mg PE/kg" = "fos",
                    "valproate 40 mg/kg" = "vpa",
                    "phenobarbital 20 mg/kg" = "pb"),
                  selected = "lev"),
      sliderInput("t_second", "투여 시각 (min)", 5, 240, 30, 5),

      hr(),
      h4("④ 3차 약제 (anaesthesia)"),
      checkboxGroupInput("third", NULL,
                         c("midazolam infusion 0.4 mg/kg/h" = "mdzinf",
                           "propofol 2 mg/kg + 6 mg/kg/h"   = "prop",
                           "ketamine 2 mg/kg + 3 mg/kg/h"   = "ket",
                           "allopregnanolone 86 ug/kg/h"    = "allo")),
      sliderInput("t_third", "시작 시각 (min)", 30, 480, 120, 10),

      hr(),
      sliderInput("tend", "시뮬레이션 길이 (min)", 120, 1440, 480, 60),
      helpText(HTML("<small>모든 탭은 같은 하나의 시뮬레이션을 읽습니다. 파라미터를 바꾸면 10개 탭이 동시에 갱신됩니다.</small>"))
    ),

    mainPanel(
      width = 9,
      fluidRow(
        column(12,
          div(style = "background:#F1F8E9;border-left:5px solid #33691E;padding:8px 12px;margin-bottom:10px",
              htmlOutput("verdict")))
      ),
      tabsetPanel(
        id = "tabs",
        tabPanel("① 환자·병인",
                 br(), plotOutput("p_gain", height = 320),
                 br(), tableOutput("t_thresholds"),
                 helpText("G = 흥분성 이득 / 억제성 긴장도. G > 1 이면 발작이 자라고, G < 1 이면 멎습니다. 병인은 분자만 움직이고, 약물은 대부분 분모만 움직입니다.")),
        tabPanel("② 치료 타임라인",
                 br(), plotOutput("p_timeline", height = 430),
                 br(), tableOutput("t_summary")),
        tabPanel("③ 수용체 시계",
                 br(), plotOutput("p_clock", height = 400),
                 br(), plotOutput("p_pools", height = 250),
                 helpText("빨간 곡선(NR_SYN)이 올라가는 동안 파란 곡선(R_SYN × F_BZS)이 내려갑니다. 초록(R_EXTRA)은 거의 움직이지 않습니다 — 늦게까지 듣는 약들이 바로 이 풀을 씁니다.")),
        tabPanel("④ 약물 PK",
                 br(), plotOutput("p_pk", height = 430),
                 br(), plotOutput("p_occ", height = 260)),
        tabPanel("⑤ 1분의 가격",
                 br(), plotOutput("p_creq", height = 380),
                 br(), tableOutput("t_creq"),
                 helpText("CREQ_BZD(t) = 지금 이 순간 G를 1로 되돌리는 데 필요한 벤조디아제핀 효과부위 농도(EC50 배수). 수직 점근선 = 어떤 용량으로도 안 되는 시각.")),
        tabPanel("⑥ 교차점",
                 br(), plotOutput("p_cross", height = 380),
                 br(), tableOutput("t_cross"),
                 helpText("동일한 점유율(0.75)에서 각 약물이 G를 몇 % 깎는지. 표적만 비교하기 위해 용량·역가는 같게 고정했습니다.")),
        tabPanel("⑦ 운동-뇌파 해리",
                 br(), plotOutput("p_diss", height = 400),
                 br(), tableOutput("t_diss"),
                 helpText("VA Cooperative: 같은 lorazepam이 현성 SE에서 64.9%, subtle SE에서 7.7-24.2%. 약이 바뀐 게 아니라 읽는 창이 바뀐 것입니다.")),
        tabPanel("⑧ 전신 생리",
                 br(), plotOutput("p_sys", height = 470),
                 helpText("30분 근처에서 보상기(카테콜아민 상승·CBF 증가)가 비보상기(자동조절 실패·저혈압·저혈당·고체온)로 넘어갑니다.")),
        tabPanel("⑨ 뇌 손상",
                 br(), plotOutput("p_injury", height = 400),
                 br(), tableOutput("t_injury"),
                 helpText("SEIZ와 INJURY는 서로 다른 상태변수입니다. '발작이 멎었다'와 '뇌를 지켰다'는 같은 말이 아닙니다.")),
        tabPanel("⑩ 가상 집단",
                 br(),
                 fluidRow(
                   column(4, numericInput("npop", "집단 크기", 200, 50, 600, 50)),
                   column(4, sliderInput("meddrive", "EDRIVE 중앙값", 1.2, 2.2, 1.55, 0.05)),
                   column(4, sliderInput("sddrive", "EDRIVE 로그 SD", 0.10, 0.60, 0.30, 0.02))),
                 actionButton("runpop", "집단 시뮬레이션 실행", class = "btn-primary"),
                 br(), br(),
                 tableOutput("t_pop"),
                 plotOutput("p_pop", height = 320),
                 helpText("보정 목표: 1차 벤조디아제핀 59-73% (PHTSE / RAMPART), 벤조 실패자에서 2차 약제 45-47% (ESETT)."))
      )
    )
  )
)

# =====================================================================
#  SERVER
# =====================================================================
server <- function(input, output, session) {

  build_events <- reactive({
    e <- NULL
    add <- function(x) if (is.null(e)) x else e + x

    tb <- input$t_bzd
    if (input$bzd != "none") {
      e <- switch(input$bzd,
                  lzp4  = d_lzp(tb, 4),
                  lzp2  = d_lzp(tb, 2),
                  mdzim = d_mdz_im(tb, 10),
                  mdziv = d_mdz_iv(tb, 0.2),
                  dzp   = d_dzp(tb, 10))
      if (input$bzd_rep) {
        e2 <- switch(input$bzd,
                     lzp4  = d_lzp(tb + 4, 4),
                     lzp2  = d_lzp(tb + 4, 2),
                     mdzim = d_mdz_im(tb + 4, 10),
                     mdziv = d_mdz_iv(tb + 4, 0.2),
                     dzp   = d_dzp(tb + 4, 10))
        e <- e + e2
      }
    }
    if (input$second != "none") {
      s2 <- switch(input$second,
                   lev = d_lev(input$t_second, 60),
                   fos = d_fos(input$t_second, 20),
                   vpa = d_vpa(input$t_second, 40),
                   pb  = d_pb(input$t_second, 20))
      e <- add(s2)
    }
    tt <- input$t_third
    dur <- max(30, input$tend - tt)
    if ("mdzinf" %in% input$third)
      e <- add(d_mdz_iv(tt, 0.2) + d_mdz_inf(tt, dur, 0.4))
    if ("prop" %in% input$third)
      e <- add(d_pro(tt, 2) + d_pro_inf(tt, dur, 6))
    if ("ket" %in% input$third)
      e <- add(d_ket(tt, 2) + d_ket_inf(tt, dur, 3))
    if ("allo" %in% input$third)
      e <- add(d_allo_inf(tt, dur, 86))
    if (is.null(e)) e <- ev(time = 0, amt = 0, cmt = 1)
    e
  })

  sim <- reactive({
    p <- list(EDRIVE = as.numeric(input$aetiology),
              ALB = input$alb, CRCL = input$crcl,
              ANAK = as.numeric(input$anak),
              GLUCIV = as.numeric(input$gluciv),
              COOL = as.numeric(input$cool),
              TDRIVE = if (input$reversible) 30 else 1e6,
              DRESID = 0.10)
    param(mod, p) %>%
      mrgsim_df(events = build_events(), end = input$tend, delta = 0.5,
                maxsteps = 2000000)
  })

  # untreated counterfactual with the same aetiology — the reference clock
  sim0 <- reactive({
    param(mod, EDRIVE = as.numeric(input$aetiology)) %>%
      mrgsim_df(events = ev(time = 0, amt = 0, cmt = 1),
                end = max(input$tend, 720), delta = 1, maxsteps = 2000000)
  })

  tstop <- reactive({
    o <- sim(); n5 <- 10
    r <- rle(o$SEIZ < 0.05)
    i <- which(r$values & r$lengths >= n5)
    if (!length(i)) return(NA_real_)
    pos <- if (i[1] == 1) 1 else sum(r$lengths[1:(i[1] - 1)]) + 1
    o$time[pos]
  })

  # ---------------- verdict banner ----------------
  output$verdict <- renderUI({
    o <- sim(); ts <- tstop()
    loss <- max(o$NEURLOSS); burden <- max(o$TSEIZ)
    if (is.na(ts)) {
      HTML(sprintf("<b style='color:#B71C1C'>발작 지속 (not terminated within %.0f min)</b> &nbsp;·&nbsp; 누적 발작 시간 %.0f분 &nbsp;·&nbsp; 해마 신경세포 손실 %.1f%% &nbsp;·&nbsp; 최저 MAP %.0f mmHg",
                   input$tend, burden, loss, min(o$MAP)))
    } else {
      HTML(sprintf("<b style='color:#1B5E20'>발작 종료 @ %.1f분</b> &nbsp;·&nbsp; 누적 발작 시간 %.1f분 &nbsp;·&nbsp; 해마 신경세포 손실 %.1f%% &nbsp;·&nbsp; 최저 MAP %.0f mmHg &nbsp;·&nbsp; 삽관 %s",
                   ts, burden, loss, min(o$MAP),
                   ifelse(any(o$INTUB > 0.5), "필요", "불필요")))
    }
  })

  # ---------------- ① gain ----------------
  output$p_gain <- renderPlot({
    o <- sim()
    d <- o %>% select(time, oEGAIN, oINH_TOT, GNET) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 1) +
      facet_wrap(~ name, scales = "free_y", ncol = 3,
                 labeller = as_labeller(c(oEGAIN = "흥분성 이득 E",
                                          oINH_TOT = "억제성 긴장도 I",
                                          GNET = "G = E / I"))) +
      scale_colour_manual(values = PALETTE) +
      labs(x = "시간 (min)", y = NULL) + theme_se + theme(legend.position = "none")
  })

  output$t_thresholds <- renderTable({
    o0 <- sim0()
    i <- function(tt) which.min(abs(o0$time - tt))
    tv <- c(5, 15, 30, 60, 120, 240)
    data.frame(`시각 (min)` = tv,
               `무치료 G` = round(o0$GNET[sapply(tv, i)], 2),
               `벤조 표적 풀` = round(o0$oTGT_BZD[sapply(tv, i)], 3),
               `NMDA 전도도` = round(o0$oNMDA[sapply(tv, i)], 3),
               `필요 벤조 농도 (EC50배)` =
                 ifelse(o0$CREQ_BZD[sapply(tv, i)] < 0, "불가능",
                        sprintf("%.2f", o0$CREQ_BZD[sapply(tv, i)])),
               check.names = FALSE)
  }, digits = 3)

  # ---------------- ② timeline ----------------
  output$p_timeline <- renderPlot({
    o <- sim()
    d <- o %>% select(time, EEG = EEGOUT, Motor = MOTOROUT,
                      `G (gain)` = GNET, `Benzo pool` = oTGT_BZD,
                      `NMDA cond.` = oNMDA, `Injury (%)` = NEURLOSS) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = 1) +
      facet_wrap(~ name, scales = "free_y", ncol = 2) +
      scale_colour_manual(values = PALETTE) +
      labs(x = "시간 (min)", y = NULL, title = "하나의 시뮬레이션, 여섯 개의 읽는 창") +
      theme_se + theme(legend.position = "none")
  })

  output$t_summary <- renderTable({
    o <- sim()
    data.frame(
      항목 = c("발작 종료 시각 (min)", "누적 발작 시간 (min)",
               "최대 세포내 Cl- (mM)", "60분 R_SYN", "60분 NMDA 전도도",
               "해마 신경세포 손실 (%)", "최고 CK (U/L)", "최저 MAP (mmHg)",
               "최저 혈당 (mmol/L)", "최고 체온 (°C)", "삽관", "뇌전증 유발 부담"),
      값 = c(ifelse(is.na(tstop()), "종료 안 됨", sprintf("%.1f", tstop())),
             sprintf("%.1f", max(o$TSEIZ)),
             sprintf("%.1f", max(o$CLI)),
             sprintf("%.3f", o$RSYN[which.min(abs(o$time - 60))]),
             sprintf("%.3f", o$oNMDA[which.min(abs(o$time - 60))]),
             sprintf("%.1f", max(o$NEURLOSS)),
             sprintf("%.0f", max(o$CK)),
             sprintf("%.0f", min(o$MAP)),
             sprintf("%.2f", min(o$GLUCP)),
             sprintf("%.2f", max(o$TEMP)),
             ifelse(any(o$INTUB > 0.5), "필요", "불필요"),
             sprintf("%.2f", max(o$EPG))))
  })

  # ---------------- ③ receptor clock ----------------
  output$p_clock <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time,
                         `R_SYN × F_BZS (벤조 표적)` = oTGT_BZD,
                         `R_SYN (시냅스 GABA-A)` = RSYN,
                         `R_EXTRA (시냅스외 δ)` = REXTRA,
                         `NR_SYN × Mg 해제 (NMDA)` = oNMDA) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) +
      geom_hline(yintercept = 1, linetype = 3, colour = "grey60") +
      geom_line(linewidth = 1.2) +
      scale_colour_manual(values = c("#1E88E5", "#64B5F6", "#1B5E20", "#E65100")) +
      labs(x = "시간 (min)", y = "기저치 대비 (fraction of baseline)",
           colour = NULL, title = "하나의 시계, 서로 반대 방향의 두 수용체") +
      theme_se
  })

  output$p_pools <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time, `KCC2` = KCC2, `[Cl-]i / 20 mM` = CLI / 20,
                         `Cl 계수 (CLFAC)` = oCLFAC,
                         `신경펩타이드 균형` = NETPEP,
                         `아데노신` = ADO) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PALETTE[c(6, 4, 1, 5, 3)]) +
      labs(x = "시간 (min)", y = NULL, colour = NULL,
           title = "억제를 갉아먹는 나머지 세 경로") + theme_se
  })

  # ---------------- ④ PK ----------------
  output$p_pk <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time,
                         `lorazepam (ng/mL)` = CP_LZP,
                         `midazolam (ng/mL)` = CP_MDZ,
                         `diazepam (ng/mL)` = CP_DZP,
                         `levetiracetam (mg/L)` = CP_LEV,
                         `phenytoin total (mg/L)` = CP_PHT,
                         `valproate total (mg/L)` = CP_VPA,
                         `phenobarbital (mg/L)` = CP_PB,
                         `ketamine (ng/mL)` = CP_KET,
                         `propofol (mg/L)` = CP_PRO,
                         `allopregnanolone (ng/mL)` = CP_ALLO) %>%
      pivot_longer(-time) %>% group_by(name) %>% filter(max(value) > 1e-6) %>% ungroup()
    if (!nrow(d)) return(NULL)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      facet_wrap(~ name, scales = "free_y") +
      scale_colour_manual(values = PALETTE) +
      labs(x = "시간 (min)", y = "혈장 농도") + theme_se +
      theme(legend.position = "none")
  })

  output$p_occ <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time, benzodiazepine = oOCC_BZD, ketamine = oOCC_KET,
                         levetiracetam = oOCC_LEV, phenytoin = oOCC_PHT,
                         valproate = oOCC_VPA, phenobarbital = oOCC_PB,
                         propofol = oOCC_PRO) %>%
      pivot_longer(-time) %>% group_by(name) %>% filter(max(value) > 1e-4) %>% ungroup()
    if (!nrow(d)) return(NULL)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      scale_colour_manual(values = PALETTE) + ylim(0, 1) +
      labs(x = "시간 (min)", y = "부위 점유율 (occupancy)", colour = NULL,
           title = "점유율은 약이 결정하고, 효과는 점유율 × 남아 있는 표적입니다") +
      theme_se
  })

  # ---------------- ⑤ price of a minute ----------------
  output$p_creq <- renderPlot({
    o0 <- sim0()
    d <- o0 %>% filter(time <= 240) %>%
      transmute(time,
                benzodiazepine = ifelse(CREQ_BZD < 0, NA, CREQ_BZD),
                phenobarbital  = ifelse(CREQ_PB  < 0, NA, CREQ_PB)) %>%
      pivot_longer(-time)
    tinf <- o0$time[which(o0$CREQ_BZD < 0)[1]]
    ggplot(d, aes(time, value, colour = name)) +
      geom_line(linewidth = 1.2, na.rm = TRUE) +
      { if (!is.na(tinf)) geom_vline(xintercept = tinf, linetype = 2, colour = "#B71C1C") } +
      { if (!is.na(tinf)) annotate("text", x = tinf, y = Inf, vjust = 1.4, hjust = -0.05,
                                   label = sprintf("어떤 벤조 용량도 불가능: %.0f분", tinf),
                                   colour = "#B71C1C", size = 4) } +
      scale_y_log10() +
      scale_colour_manual(values = c("#1E88E5", "#E65100")) +
      labs(x = "시간 (min)", y = "G=1 을 만드는 데 필요한 효과부위 농도 (EC50 배수, log)",
           colour = NULL, title = "1분의 가격은 일정하지 않습니다") + theme_se
  })

  output$t_creq <- renderTable({
    o0 <- sim0(); tv <- c(5, 10, 20, 30, 45, 60, 120, 240)
    i <- sapply(tv, function(tt) which.min(abs(o0$time - tt)))
    data.frame(`시각 (min)` = tv,
               `R_SYN` = round(o0$RSYN[i], 3),
               `F_BZS` = round(o0$FBZS[i], 3),
               `벤조 표적` = round(o0$oTGT_BZD[i], 3),
               `필요 벤조 (EC50배)` = ifelse(o0$CREQ_BZD[i] < 0, "불가능", sprintf("%.2f", o0$CREQ_BZD[i])),
               `필요 페노바르비탈 (EC50배)` = ifelse(o0$CREQ_PB[i] < 0, "불가능", sprintf("%.2f", o0$CREQ_PB[i])),
               check.names = FALSE)
  })

  # ---------------- ⑥ crossover ----------------
  cross <- reactive({
    o <- sim0(); p <- as.list(param(mod)); occ <- 0.75
    tv <- seq(5, min(720, max(o$time)), by = 5)
    i <- sapply(tv, function(tt) which.min(abs(o$time - tt)))
    inh0 <- o$oINH_TOT[i]; eg0 <- o$oEGAIN[i]; g0 <- eg0 / inh0
    tgt_b <- o$oTGT_BZD[i]
    tgt_p <- p$FSYNPB * o$RSYN[i] + (1 - p$FSYNPB) * o$REXTRA[i]
    blk <- p$EMAXKET * occ * (p$UDMIN + (1 - p$UDMIN) * o$SEIZ[i])
    data.frame(time = tv,
               benzodiazepine = 1 - (eg0 / (inh0 + p$EMAXBZD * occ * tgt_b)) / g0,
               phenobarbital  = 1 - (eg0 / (inh0 + p$EMAXPB  * occ * tgt_p)) / g0,
               ketamine       = blk * o$KETSHARE[i])
  })

  output$p_cross <- renderPlot({
    d <- cross(); dl <- pivot_longer(d, -time)
    xc <- d$time[which(d$ketamine > d$benzodiazepine)[1]]
    ggplot(dl, aes(time, value, colour = name)) + geom_line(linewidth = 1.2) +
      { if (!is.na(xc)) geom_vline(xintercept = xc, linetype = 2, colour = "grey30") } +
      { if (!is.na(xc)) annotate("text", x = xc, y = 0.05, hjust = -0.05,
                                 label = sprintf("교차: %.0f분", xc), size = 4) } +
      scale_colour_manual(values = c("#1E88E5", "#E65100", "#1B5E20")) +
      labs(x = "시간 (min)", y = "동일 점유율(0.75)에서 G 감소분", colour = NULL,
           title = "표적만 비교했을 때 — 하나는 떠나고 하나는 도착합니다") + theme_se
  })

  output$t_cross <- renderTable({
    d <- cross(); tv <- c(5, 30, 60, 120, 240, 480)
    d <- d[sapply(tv, function(tt) which.min(abs(d$time - tt))), ]
    data.frame(`시각 (min)` = d$time,
               `벤조` = round(d$benzodiazepine, 3),
               `페노바르비탈` = round(d$phenobarbital, 3),
               `케타민` = round(d$ketamine, 3), check.names = FALSE)
  })

  # ---------------- ⑦ dissociation ----------------
  output$p_diss <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time, `EEG 부담 (electrographic)` = EEGOUT,
                         `운동 발현 (motor)` = MOTOROUT,
                         `운동 이득 MOTG` = MOTG,
                         `벤조 점유율` = oOCC_BZD) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 0.15, linetype = 3, colour = "grey50") +
      scale_colour_manual(values = c("#880E4F", "#B71C1C", "#546E7A", "#1E88E5")) +
      labs(x = "시간 (min)", y = NULL, colour = NULL,
           title = "침상에서 보이는 것과 뇌가 하고 있는 것") + theme_se
  })

  output$t_diss <- renderTable({
    o <- sim(); tv <- c(15, 30, 60, 90, 120, 180, 240)
    tv <- tv[tv <= max(o$time)]
    i <- sapply(tv, function(tt) which.min(abs(o$time - tt)))
    data.frame(`시각 (min)` = tv,
               `운동 출력` = round(o$MOTOROUT[i], 3),
               `EEG 부담` = round(o$EEGOUT[i], 3),
               `침상 판정` = ifelse(o$MOTOROUT[i] < 0.15, "조절된 것처럼 보임", "발작 지속"),
               `실제` = ifelse(o$EEGOUT[i] < 0.05, "조절됨", "발작 지속"),
               check.names = FALSE)
  })

  # ---------------- ⑧ systemic ----------------
  output$p_sys <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time, `MAP (mmHg)` = MAP, `CBF (relative)` = oCBF,
                         `공급-수요 불일치` = oMISMATCH,
                         `혈당 (mmol/L)` = GLUCP, `젖산 (mmol/L)` = LAC,
                         `체온 (°C)` = TEMP, `CK (U/L)` = CK,
                         `자동조절 AUTO` = AUTO, `호흡 부담` = RESPD) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      facet_wrap(~ name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = PALETTE) +
      labs(x = "시간 (min)", y = NULL) + theme_se + theme(legend.position = "none")
  })

  # ---------------- ⑨ injury ----------------
  output$p_injury <- renderPlot({
    o <- sim()
    d <- o %>% transmute(time, `글루타메이트 (uM)` = GLU, `[Ca2+]i` = CAI,
                         `누적 손상` = INJURY, `해마 신경세포 손실 (%)` = NEURLOSS,
                         `IL-1beta` = IL1B, `뇌전증 유발` = EPG) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      facet_wrap(~ name, scales = "free_y", ncol = 3) +
      scale_colour_manual(values = PALETTE) +
      labs(x = "시간 (min)", y = NULL) + theme_se + theme(legend.position = "none")
  })

  output$t_injury <- renderTable({
    o <- sim()
    data.frame(
      항목 = c("누적 발작 시간 (min)", "최고 [Ca2+]i", "누적 손상 (a.u.)",
               "해마 신경세포 손실 (%)", "최고 IL-1beta", "BBB 투과성 배수",
               "P-gp 유도 배수", "뇌:혈장 분배 지수", "뇌전증 유발 부담"),
      값 = round(c(max(o$TSEIZ), max(o$CAI), max(o$INJURY), max(o$NEURLOSS),
                   max(o$IL1B), max(o$BBBP), max(o$PGP), min(o$BPRATIO),
                   max(o$EPG)), 3))
  })

  # ---------------- ⑩ virtual population ----------------
  popres <- eventReactive(input$runpop, {
    n <- input$npop
    set.seed(20260804)
    id <- data.frame(ID = 1:n, EDRIVE = input$meddrive * exp(rnorm(n, 0, input$sddrive)))
    o1 <- mod %>% mrgsim_df(events = d_lzp(20, 4), idata = id, end = 100, delta = 1,
                            maxsteps = 2000000, obsonly = TRUE)
    f1 <- o1 %>% group_by(ID) %>%
      summarise(ok = any(SEIZ[time >= 20 & time <= 80] < 0.05), .groups = "drop")
    id2 <- id[id$ID %in% f1$ID[!f1$ok], ]
    rate2 <- function(e) {
      if (!nrow(id2)) return(NA_real_)
      o <- mod %>% mrgsim_df(events = e, idata = id2, end = 150, delta = 1,
                             maxsteps = 2000000, obsonly = TRUE)
      100 * mean((o %>% group_by(ID) %>%
                    summarise(ok = any(SEIZ[time >= 65 & time <= 125] < 0.05),
                              .groups = "drop"))$ok)
    }
    im <- mod %>% mrgsim_df(events = d_mdz_im(8, 10), idata = id, end = 100, delta = 1,
                            maxsteps = 2000000, obsonly = TRUE)
    iv <- mod %>% mrgsim_df(events = d_lzp(12.5, 4), idata = id, end = 100, delta = 1,
                            maxsteps = 2000000, obsonly = TRUE)
    rr <- function(o, a, b) 100 * mean((o %>% group_by(ID) %>%
      summarise(ok = any(SEIZ[time >= a & time <= b] < 0.05), .groups = "drop"))$ok)
    list(
      tab = data.frame(
        endpoint = c("1차 벤조디아제핀 (20분)", "ESETT levetiracetam",
                     "ESETT fosphenytoin", "ESETT valproate",
                     "RAMPART IM midazolam", "RAMPART IV lorazepam"),
        model_pct = round(c(100 * mean(f1$ok),
                            rate2(d_lzp(20, 4) + d_lev(65, 60)),
                            rate2(d_lzp(20, 4) + d_fos(65, 20)),
                            rate2(d_lzp(20, 4) + d_vpa(65, 40)),
                            rr(im, 8, 68), rr(iv, 12.5, 72.5)), 1),
        observed_pct = c("59-73", "47", "45", "46", "73.4", "63.4")),
      first = do.call(rbind, lapply(c(5, 10, 20, 30, 45, 60), function(tt) {
        o <- mod %>% mrgsim_df(events = d_lzp(tt, 4), idata = id, end = tt + 70,
                               delta = 1, maxsteps = 2000000, obsonly = TRUE)
        data.frame(line = "1차 (benzodiazepine)", t = tt, pct = rr(o, tt, tt + 60))
      })),
      second = if (nrow(id2)) do.call(rbind, lapply(c(30, 45, 60, 90, 120), function(tt) {
        o <- mod %>% mrgsim_df(events = d_lzp(20, 4) + d_lev(tt, 60), idata = id2,
                               end = tt + 70, delta = 1, maxsteps = 2000000, obsonly = TRUE)
        data.frame(line = "2차 (levetiracetam, 벤조 실패자)", t = tt, pct = rr(o, tt, tt + 60))
      })) else NULL)
  })

  output$t_pop <- renderTable({ popres()$tab })

  output$p_pop <- renderPlot({
    r <- popres(); d <- rbind(r$first, r$second)
    ggplot(d, aes(t, pct, colour = line)) +
      geom_line(linewidth = 1.2) + geom_point(size = 2.5) +
      scale_colour_manual(values = c("#1E88E5", "#1B5E20")) + ylim(0, 100) +
      labs(x = "투여 시각 (min)", y = "반응률 (%)", colour = NULL,
           title = "두 개의 절벽 — 그리고 그 사이의 40분") + theme_se
  })
}

shinyApp(ui, server)
