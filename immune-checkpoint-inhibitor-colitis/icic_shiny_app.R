## ---------------------------------------------------------------------------
## icic_shiny_app.R
## ===========================================================================
## Interactive dashboard for the IMMUNE CHECKPOINT INHIBITOR COLITIS QSP model.
##
## The app is built around the model's central claim, so the controls are
## arranged to let a user break it:
##
##   Phi_col = (a_eff*Teff + a_trm*Trm) / (a_reg*Treg + kappa)
##
##   - move the ANTI-PD-1 dose slider across its whole range and watch the
##     occupancy trace refuse to move (tab 3).  That is the flat dose-toxicity
##     curve, drawn.
##   - move the IPILIMUMAB dose slider and watch the ADCC trace move
##     proportionally while its own occupancy trace stays flat (tab 3).
##   - watch S_eff fall for weeks with a flat stool chart (tab 4).  The dashed
##     line is S* = L/A_max: nothing is reported until the curve crosses it.
##   - choose a rescue agent and read the onset time off tab 5, and the
##     tumour cost off tab 7.  They are not the same ranking.
##
## Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
## Run:  shiny::runApp("icic_shiny_app.R")
## ---------------------------------------------------------------------------

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

MODEL_FILE <- "icic_mrgsolve_model.R"
BW <- 70

mod <- mread_cache("icic", MODEL_FILE)

## Phenotype presets. The MEDIAN patient does not get colitis at any of these
## doses -- only ~7-12% of ipilimumab-treated patients do -- so a dashboard
## that opened on the median would look broken. It opens on the susceptible
## phenotype and says so.
PHENO <- list(
  "Median patient (does NOT develop colitis — this is correct)" =
    list(par = list(ag_scale = 1.00, kappa = 0.15, phi_FCGR = 1.00),
         ini = list(Dv = 1.00, Bu = 1.00, Trm = 0.02)),
  "Susceptible (FCGR3A V/V, thin regulatory floor, prior antibiotics)" =
    list(par = list(ag_scale = 1.45, kappa = 0.11, phi_FCGR = 1.60),
         ini = list(Dv = 0.70, Bu = 0.70^1.3, Trm = 0.06)),
  "Highly susceptible (steroid-refractory phenotype)" =
    list(par = list(ag_scale = 1.75, kappa = 0.095, phi_FCGR = 1.60),
         ini = list(Dv = 0.60, Bu = 0.60^1.3, Trm = 0.12)),
  "Pre-existing microscopic colitis (latent primed pool)" =
    list(par = list(ag_scale = 1.30, kappa = 0.13, phi_FCGR = 1.00),
         ini = list(Dv = 0.85, Bu = 0.85^1.3, Trm = 0.30, Ent = 0.85)),
  "FCGR3A F/F (low-affinity FcgammaRIIIa — inefficient ADCC)" =
    list(par = list(ag_scale = 1.45, kappa = 0.11, phi_FCGR = 0.55),
         ini = list(Dv = 0.70, Bu = 0.70^1.3, Trm = 0.06))
)

## Prednisolone as a continuous mg/day input: a 3 h half-life drug with a
## ~12 h genomic effect half-life is pharmacodynamically indistinguishable
## from a daily dose, and this keeps the taper a smooth covariate.
pred_rate <- function(tau, top, plateau, taper_wk) {
  ifelse(tau < 0, 0,
         ifelse(tau < plateau, top,
                {
                  wk <- (tau - plateau) / 7
                  ifelse(wk >= taper_wk, 0, top * (1 - wk / taper_wk))
                }))
}

build_data <- function(inp) {
  tt <- seq(0, inp$tend, by = 1)
  cov <- data.frame(ID = 1, time = tt, evid = 0, amt = 0, cmt = 0,
                    PRED_RATE = pred_rate(tt - inp$ster_day, inp$ster_mgkg * BW,
                                          inp$ster_plateau, inp$ster_taper),
                    JAK_RATE = ifelse(inp$rescue == "tofacitinib" &
                                        tt >= inp$resc_day &
                                        tt < inp$resc_day + 60, 20, 0),
                    ABX_ON = ifelse(tt < inp$abx_days, 1, 0))

  ev <- data.frame()
  add <- function(times, amt, cmt)
    data.frame(ID = 1, time = times, amt = amt, cmt = cmt, evid = 1)

  if (inp$ipi_mgkg > 0)
    ev <- bind_rows(ev, add(seq(0, by = 21, length.out = inp$ipi_n),
                            inp$ipi_mgkg * BW, 1))
  if (inp$pd1_mgkg > 0)
    ev <- bind_rows(ev, add(seq(0, by = 14, length.out = inp$pd1_n),
                            inp$pd1_mgkg * BW, 3))
  if (inp$rescue == "infliximab")
    ev <- bind_rows(ev, add(inp$resc_day + c(0, 14, 42), 5 * BW, 5))
  if (inp$rescue == "vedolizumab")
    ev <- bind_rows(ev, add(inp$resc_day + c(0, 14, 42), 300, 7))
  if (inp$rescue == "tocilizumab")
    ev <- bind_rows(ev, add(inp$resc_day + c(0, 28), 8 * BW, 9))

  if (nrow(ev)) {
    ev$PRED_RATE <- NA; ev$JAK_RATE <- NA; ev$ABX_ON <- NA
    out <- bind_rows(cov, ev) %>% arrange(time, desc(evid)) %>%
      fill(PRED_RATE, JAK_RATE, ABX_ON, .direction = "downup")
  } else out <- cov
  out
}

run_sim <- function(inp) {
  ph <- PHENO[[inp$pheno]]
  m <- mod %>% param(ph$par) %>% init(ph$ini)
  if (inp$fmt_day > 0) {
    ## FMT is a reset of the antigen/SCFA arm, so it enters as a step change
    ## in Dv rather than as a drug.
    d1 <- build_data(inp) %>% filter(time <= inp$fmt_day)
    r1 <- m %>% data_set(d1) %>% mrgsim(end = inp$fmt_day, delta = 0.5)
    st <- as.list(tail(as_tibble(r1), 1))
    st$Dv <- min(1, st$Dv + 0.60)
    keep <- names(init(m))
    m2 <- m %>% init(st[keep])
    d2 <- build_data(inp) %>% filter(time > inp$fmt_day) %>%
      mutate(time = time - inp$fmt_day)
    r2 <- m2 %>% data_set(d2) %>% mrgsim(end = inp$tend - inp$fmt_day,
                                         delta = 0.5)
    return(bind_rows(as_tibble(r1),
                     as_tibble(r2) %>% mutate(time = time + inp$fmt_day)))
  }
  m %>% data_set(build_data(inp)) %>% mrgsim(end = inp$tend, delta = 0.5) %>%
    as_tibble()
}

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10))

## ===========================================================================
## UI — eight tabs
## ===========================================================================
ui <- fluidPage(
  titlePanel(HTML(paste0(
    "<b>면역관문억제제 유발 대장염 QSP 대시보드</b>",
    "<br><span style='font-size:14px;color:#555'>Immune Checkpoint Inhibitor",
    "-Induced Colitis — a ratio with a saturated numerator and a depletable",
    " denominator</span>"))),
  tags$hr(),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("pheno", "환자 표현형 (patient phenotype)",
                  choices = names(PHENO), selected = names(PHENO)[2]),
      tags$hr(),
      tags$b("면역관문억제제 (checkpoint inhibitor)"),
      sliderInput("ipi_mgkg", "ipilimumab (mg/kg q3w)", 0, 10, 3, step = 0.5),
      sliderInput("ipi_n", "ipilimumab 투여 횟수", 0, 8, 4, step = 1),
      sliderInput("pd1_mgkg", "anti-PD-1 (mg/kg q2w)", 0, 10, 0, step = 0.5),
      sliderInput("pd1_n", "anti-PD-1 투여 횟수", 0, 20, 13, step = 1),
      helpText(HTML(paste0("<i>anti-PD-1 슬라이더를 끝까지 움직여도 3번 탭의",
                           " 점유율 곡선은 거의 변하지 않습니다. 그것이 평탄한",
                           " 용량-독성 관계입니다.</i>"))),
      tags$hr(),
      tags$b("구조 요법 (rescue therapy)"),
      sliderInput("ster_day", "스테로이드 시작일", 0, 180, 53, step = 1),
      sliderInput("ster_mgkg", "prednisolone (mg/kg/day)", 0, 2, 1, step = 0.1),
      sliderInput("ster_plateau", "고정 용량 유지일", 0, 21, 7, step = 1),
      sliderInput("ster_taper", "감량 기간 (주)", 1, 12, 4, step = 1),
      radioButtons("rescue", "추가 구조 약제",
                   c("none", "infliximab", "vedolizumab", "tocilizumab",
                     "tofacitinib"), selected = "none"),
      sliderInput("resc_day", "구조 약제 시작일", 0, 200, 60, step = 1),
      sliderInput("fmt_day", "FMT 시행일 (0 = 시행 안 함)", 0, 200, 0,
                  step = 1),
      tags$hr(),
      sliderInput("abx_days", "치료 전 항생제 노출 (일)", 0, 30, 0, step = 1),
      sliderInput("tend", "관찰 기간 (일)", 90, 365, 210, step = 15),
      tags$hr(),
      helpText(HTML(paste0("<span style='font-size:11px'>교육·연구 목적의 ",
                           "모델입니다. 임상 의사결정에 사용하지 마십시오.",
                           "</span>")))
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        ## ---- 1 -----------------------------------------------------------
        tabPanel(
          "1. 환자 프로파일",
          br(),
          fluidRow(column(12, uiOutput("verdict"))),
          br(),
          fluidRow(column(6, plotOutput("p_ratio", height = 300)),
                   column(6, plotOutput("p_cells", height = 300))),
          br(),
          fluidRow(column(12, DTOutput("t_summary")))
        ),
        ## ---- 2 -----------------------------------------------------------
        tabPanel(
          "2. 약동학 (PK)",
          br(),
          plotOutput("p_pk", height = 360),
          br(),
          plotOutput("p_clifx", height = 300),
          helpText(HTML(paste0("<b>질병이 자기 해독제를 먹는다.</b> ",
                               "중증 대장염의 단백소실장병증이 알부민을 떨어뜨리고,",
                               " 알부민은 단일클론항체 청소율의 가장 큰 공변량입니다",
                               " (CL ∝ ALB<sup>-0.9</sup>). 같은 5 mg/kg가 가장",
                               " 아픈 환자에게 가장 낮은 노출을 전달합니다.")))
        ),
        ## ---- 3 -----------------------------------------------------------
        tabPanel(
          "3. 표적 결합 — 두 개의 다른 노출-반응 형태",
          br(),
          plotOutput("p_occ", height = 340),
          br(),
          fluidRow(column(6, plotOutput("p_doseresp", height = 320)),
                   column(6, plotOutput("p_G", height = 320))),
          helpText(HTML(paste0("PD-1 점유율은 임상 용량 전 범위에서 0.98-0.999",
                               "입니다 — 찾을 용량-반응이 없습니다. ADCC는 포화되지",
                               " 않으므로 ipilimumab의 용량-반응은 전부 <b>분모</b>",
                               "에서 나옵니다.")))
        ),
        ## ---- 4 -----------------------------------------------------------
        tabPanel(
          "4. 흡수 예비능과 검열된 증상",
          br(),
          plotOutput("p_reserve", height = 380),
          br(),
          plotOutput("p_stool", height = 300),
          helpText(HTML(paste0("<b>S* = L/A_max = 0.389.</b> 대장 흡수 기능의 ",
                               "61%가 완전히 정상적인 배변 기록과 함께 소실될 수",
                               " 있습니다. 등급은 늦고 손실이 큰 계기이며, 증상",
                               "-유발 치료는 태생적으로 늦습니다.")))
        ),
        ## ---- 5 -----------------------------------------------------------
        tabPanel(
          "5. 임상 엔드포인트",
          br(),
          plotOutput("p_grade", height = 320),
          br(),
          fluidRow(column(6, plotOutput("p_epi", height = 300)),
                   column(6, plotOutput("p_ulcer", height = 300))),
          br(),
          DTOutput("t_endpoints")
        ),
        ## ---- 6 -----------------------------------------------------------
        tabPanel(
          "6. 바이오마커",
          br(),
          plotOutput("p_biomark", height = 360),
          br(),
          plotOutput("p_cyto", height = 320),
          helpText(HTML(paste0("칼프로텍틴과 내시경은 <b>병변</b>을 읽으므로 ",
                               "검열되지 않습니다. 그래서 증상보다 먼저 움직입니다.")))
        ),
        ## ---- 7 -----------------------------------------------------------
        tabPanel(
          "7. 종양 — 선택성이 통화다",
          br(),
          plotOutput("p_tumor", height = 340),
          br(),
          fluidRow(column(6, plotOutput("p_imm", height = 300)),
                   column(6, plotOutput("p_haz", height = 300))),
          helpText(HTML(paste0("χ (비-장선택성): prednisolone 1.00 · ",
                               "infliximab 0.60 · vedolizumab 0.10. ",
                               "α4β7:MAdCAM-1은 종양이 쓰지 않는 <b>장 주소</b>",
                               "이므로, 같은 대장염 조절을 훨씬 싸게 살 수 있습니다.")))
        ),
        ## ---- 8 -----------------------------------------------------------
        tabPanel(
          "8. 시나리오 비교",
          br(),
          checkboxGroupInput(
            "cmp", "비교할 시나리오",
            choices = c("ipilimumab 1 mg/kg", "ipilimumab 3 mg/kg",
                        "ipilimumab 10 mg/kg", "anti-PD-1 3 mg/kg",
                        "anti-PD-1 10 mg/kg", "ipi 3 + nivo 1"),
            selected = c("ipilimumab 1 mg/kg", "ipilimumab 3 mg/kg",
                         "ipilimumab 10 mg/kg", "anti-PD-1 3 mg/kg"),
            inline = TRUE),
          plotOutput("p_cmp", height = 420),
          br(),
          DTOutput("t_cmp"),
          helpText(HTML(paste0("두 가지를 확인하십시오. (1) anti-PD-1 3 vs 10 ",
                               "mg/kg는 사실상 겹칩니다. (2) 병용은 두 단독의 합이",
                               " 아닙니다 — 등급이 역치 지표이고, ipilimumab이 이미",
                               " 취약한 환자를 그 역치 너머로 밀어놓았기 때문입니다.")))
        )
      )
    )
  )
)

## ===========================================================================
## SERVER
## ===========================================================================
server <- function(input, output, session) {

  sim <- reactive({
    inp <- reactiveValuesToList(input)
    run_sim(inp)
  })

  ## ---- 1. verdict banner --------------------------------------------------
  output$verdict <- renderUI({
    s <- sim()
    pg <- max(s$grade_o); mn <- min(s$S_eff_o)
    d2 <- s$time[which(s$grade_o >= 2)][1]
    dc <- s$time[which(s$Calpro >= 200)][1]
    col <- c("#2e7d32", "#f9a825", "#ef6c00", "#c62828")[pg + 1]
    HTML(sprintf(paste0(
      "<div style='padding:14px;border-left:6px solid %s;background:#fafafa'>",
      "<span style='font-size:20px;font-weight:bold;color:%s'>최고 등급 G%d",
      "</span><br><span style='color:#444'>최저 흡수 기능 S_eff = %.3f ",
      "(역치 S* = %.3f) · 잔여 예비능 %.0f%%<br>",
      "칼프로텍틴 &gt;200 µg/g: %s일 · 등급 2 도달: %s일 · <b>선행 시간 %s일</b>",
      "<br>최저 알부민 %.2f g/dL · 최대 궤양 지수 %.2f · 누적 스테로이드 %.0f mg",
      "</span></div>"),
      col, col, pg, mn, s$S_star_o[1], 100 * max(0, mn - s$S_star_o[1]) /
        s$S_star_o[1],
      ifelse(is.na(dc), "—", format(dc)), ifelse(is.na(d2), "—", format(d2)),
      ifelse(is.na(d2) || is.na(dc), "—", format(d2 - dc)),
      min(s$Alb), max(s$Ulcer), max(s$SterCum)))
  })

  output$p_ratio <- renderPlot({
    s <- sim()
    s %>% select(time, `G = Reg0/Reg (조절 해제)` = G_o,
                 `Ppd1 (PD-1 해제)` = Ppd1_o) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 1, linetype = 3) +
      labs(x = "일", y = "배수 (baseline = 1)",
           title = "비(ratio)의 두 인자",
           subtitle = "분모는 고갈되고 (G가 쌍곡선으로 상승), 분자는 포화된다") +
      THEME
  })

  output$p_cells <- renderPlot({
    s <- sim()
    s %>% select(time, Teff, Trm, Treg, Nclone) %>% pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "일", y = "상대량", title = "세포 구획",
           subtitle = "Trm은 약보다 오래 남는다 — 방아쇠가 아니라 상태다") +
      THEME
  })

  output$t_summary <- renderDT({
    s <- sim()
    data.frame(
      지표 = c("최고 CTCAE 설사 등급", "최저 S_eff", "검열 역치 S*",
               "최고 칼프로텍틴 (µg/g)", "최고 CRP (mg/L)",
               "최저 알부민 (g/dL)", "최대 궤양 지수", "최저 은와(ISC) 지수",
               "d180 Trm", "누적 스테로이드 (mg)", "누적 감염 위험",
               "d180 종양 부담"),
      값 = c(max(s$grade_o), round(min(s$S_eff_o), 3), round(s$S_star_o[1], 3),
             round(max(s$Calpro)), round(max(s$CRP), 1), round(min(s$Alb), 2),
             round(max(s$Ulcer), 2), round(min(s$ISC), 3),
             round(s$Trm[nrow(s)], 3), round(max(s$SterCum)),
             round(max(s$InfHaz), 3), round(s$Tumor[nrow(s)], 1))) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 20))
  })

  ## ---- 2. PK --------------------------------------------------------------
  output$p_pk <- renderPlot({
    s <- sim()
    s %>% select(time, ipilimumab = C_ipi_o, `anti-PD-1` = C_pd1_o,
                 infliximab = C_ifx_o, vedolizumab = C_vdz_o) %>%
      pivot_longer(-time) %>% filter(value > 1e-6) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      labs(x = "일", y = "혈장 농도 (µg/mL)", title = "단일클론항체 약동학") +
      THEME
  })

  output$p_clifx <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_line(aes(y = CL_ifx_o, colour = "infliximab CL (L/일)"),
                linewidth = 0.9) +
      geom_line(aes(y = Alb / 10, colour = "알부민 / 10 (g/dL)"),
                linewidth = 0.9) +
      labs(x = "일", y = NULL,
           title = "질병 중증도가 자기 해독제의 청소율을 올린다",
           subtitle = "CL_ifx = 0.40 · (ALB/4.0)^-0.9 · (1 + 0.006·CRP)") +
      THEME
  })

  ## ---- 3. target engagement ----------------------------------------------
  output$p_occ <- renderPlot({
    s <- sim()
    s %>% select(time, `PD-1 점유율 (포화)` = O_PD1_o,
                 `CTLA-4 점유율 (포화)` = O_CTLA4_o,
                 `ADCC 구동 (미포화)` = A_adcc_o) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 1) +
      ylim(0, 1) +
      labs(x = "일", y = "최대치 대비 분율",
           title = "왜 한 약은 용량-반응이 있고 다른 약은 없는가",
           subtitle = paste0("점유율은 어느 임상 용량에서도 포화 상태다. ",
                             "포화되지 않는 것은 ADCC뿐이다.")) +
      THEME
  })

  output$p_doseresp <- renderPlot({
    ## computed analytically from the same equations, not simulated
    p <- as.list(param(mod))
    d <- seq(0.1, 10, length.out = 120)
    pd1 <- d * BW / (p$CL_pd1 * 14); ipi <- d * BW / (p$CL_ipi * 21)
    tibble(
      dose = rep(d, 2),
      value = c((pd1 * p$BDC * 1e3 / p$MW_pd1) /
                  (pd1 * p$BDC * 1e3 / p$MW_pd1 + p$KD_PD1_nM),
                (ipi * p$BDC) / (ipi * p$BDC + p$EC50_ADCC)),
      what = rep(c("anti-PD-1 점유율", "ipilimumab ADCC 구동"), each = 120)
    ) %>% ggplot(aes(dose, value, colour = what)) +
      geom_line(linewidth = 1) + ylim(0, 1) +
      scale_x_log10(breaks = c(0.1, 0.3, 1, 3, 10)) +
      labs(x = "용량 (mg/kg, 로그)", y = "최대치 대비 분율",
           title = "노출-반응 형태 두 가지",
           subtitle = "하나는 30배 용량 범위 전체에서 평평하다") +
      THEME
  })

  output$p_G <- renderPlot({
    p <- as.list(param(mod))
    treg <- seq(0.02, 1, length.out = 200)
    tibble(treg, G = (p$a_reg * 1 + p$kappa) / (p$a_reg * treg + p$kappa)) %>%
      ggplot(aes(treg, G)) + geom_line(linewidth = 1, colour = "#c0392b") +
      geom_hline(yintercept = (1 + p$kappa) / p$kappa, linetype = 2) +
      annotate("text", x = 0.55, y = (1 + p$kappa) / p$kappa * 0.93, size = 3.4,
               label = sprintf("천장 G_max = Reg0/kappa = %.1f\n(비-Treg 바닥 kappa가 정한다)",
                               (1 + p$kappa) / p$kappa)) +
      labs(x = "잔존 대장 Treg (기저 대비)", y = "조절 해제 계수 G",
           title = "고갈 가능한 분모는 쌍곡선이다",
           subtitle = "모든 점유율 항이 평평해진 뒤에도 G는 계속 오른다") +
      THEME
  })

  ## ---- 4. reserve ---------------------------------------------------------
  output$p_reserve <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_ribbon(aes(ymin = S_star_o, ymax = pmax(S_eff_o, S_star_o)),
                  fill = "#d4f0ef", alpha = 0.65) +
      geom_line(aes(y = S_eff_o), linewidth = 1, colour = "#17807a") +
      geom_hline(aes(yintercept = S_star_o), linetype = 2) +
      annotate("text", x = max(s$time) * 0.02, y = s$S_star_o[1] + 0.04,
               hjust = 0, size = 3.6,
               label = "S* = L/A_max = 0.389 — 이 아래로 내려가야 비로소 배변 기록이 움직인다") +
      ylim(0, 1.02) +
      labs(x = "일", y = "유효 흡수 기능 S_eff",
           title = "첫 증상 전에 소진되는 예비능",
           subtitle = "음영 부분이 남아 있는 예비능. 증상은 0이고 병변은 진행 중이다") +
      THEME
  })

  output$p_stool <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_line(aes(y = Stool_o, colour = "변 수분량 (mL/일)"), linewidth = 0.9) +
      geom_line(aes(y = Calpro, colour = "칼프로텍틴 (µg/g)"), linewidth = 0.9) +
      geom_hline(yintercept = 200, linetype = 3) +
      labs(x = "일", y = NULL,
           title = "검열된 계기와 검열되지 않은 계기",
           subtitle = "칼프로텍틴이 먼저 움직인다. 그것이 조기 개입의 근거다") +
      THEME
  })

  ## ---- 5. endpoints -------------------------------------------------------
  output$p_grade <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_step(aes(y = grade_o), linewidth = 1, colour = "#7e5109") +
      geom_line(aes(y = dfreq_o), linewidth = 0.7, colour = "#b9770e",
                linetype = 2) +
      scale_y_continuous("CTCAE 등급 (실선) · 배변 증가 횟수/일 (점선)",
                         breaks = 0:12) +
      labs(x = "일", title = "임상 엔드포인트",
           subtitle = "G1 <4회 · G2 4-6회 · G3 ≥7회 (기저 대비 증가)") +
      THEME
  })

  output$p_epi <- renderPlot({
    s <- sim()
    s %>% select(time, `대장세포 질량 Ent` = Ent, `밀착연접 TJ` = TJ,
                 `은와 줄기세포 ISC` = ISC, `점액 Muc` = Muc) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      ylim(0, 1.05) +
      labs(x = "일", y = "기저 대비", title = "상피 구획",
           subtitle = "ISC 회복은 20일 시상수 — 스테로이드가 되돌리지 못하는 것") +
      THEME
  })

  output$p_ulcer <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_line(aes(y = Ulcer, colour = "궤양 지수"), linewidth = 0.9) +
      geom_line(aes(y = Alb, colour = "알부민 (g/dL)"), linewidth = 0.9) +
      labs(x = "일", y = NULL, title = "궤양과 단백 소실",
           subtitle = "깊은 궤양은 스테로이드 불응성의 예측 인자다") +
      THEME
  })

  output$t_endpoints <- renderDT({
    s <- sim()
    first_at <- function(v, thr) {
      i <- which(v >= thr)[1]; if (is.na(i)) NA else s$time[i]
    }
    data.frame(
      계기 = c("병변 (S_eff가 기저의 90% 미만)", "칼프로텍틴 >150 µg/g",
               "칼프로텍틴 >200 µg/g", "CTCAE 등급 1", "CTCAE 등급 2",
               "CTCAE 등급 3"),
      최초일 = c(s$time[which(s$S_eff_o < 0.9 * s$S_eff_o[1])[1]],
                 first_at(s$Calpro, 150), first_at(s$Calpro, 200),
                 first_at(s$grade_o, 1), first_at(s$grade_o, 2),
                 first_at(s$grade_o, 3))) %>%
      datatable(rownames = FALSE, options = list(dom = "t"),
                caption = "검열 사다리: 왼쪽에서 오른쪽으로 읽으면 예비능이 소진되는 과정이다")
  })

  ## ---- 6. biomarkers ------------------------------------------------------
  output$p_biomark <- renderPlot({
    s <- sim()
    s %>% select(time, `칼프로텍틴 (µg/g)` = Calpro, `CRP (mg/L)` = CRP,
                 `알부민 ×100 (g/dL)` = Alb) %>%
      mutate(`알부민 ×100 (g/dL)` = `알부민 ×100 (g/dL)` * 100) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      facet_wrap(~name, scales = "free_y", ncol = 3) +
      labs(x = "일", y = NULL, title = "바이오마커") +
      THEME + theme(legend.position = "none")
  })

  output$p_cyto <- renderPlot({
    s <- sim()
    s %>% select(time, `IFN-γ` = IFNg, `TNF-α` = TNFa, `IL-6` = IL6,
                 `IL-15 (Trm 연료)` = IL15, `CXCL10` = CXCL10,
                 `MAdCAM-1` = MADCAM) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.85) +
      labs(x = "일", y = "기저 대비 배수", title = "사이토카인 망",
           subtitle = "IFN-γ가 프로그램을 정하고 TNF-α가 집행한다") +
      THEME
  })

  ## ---- 7. tumour ----------------------------------------------------------
  output$p_tumor <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_line(aes(y = Tumor, colour = "종양 부담"), linewidth = 1) +
      geom_line(aes(y = TumTeff * 200, colour = "종양내 효과기 ×200"),
                linewidth = 0.85) +
      labs(x = "일", y = NULL, title = "구조 요법의 종양 비용",
           subtitle = "대장을 상하게 하는 클론이 종양을 죽이는 클론이다") +
      THEME
  })

  output$p_imm <- renderPlot({
    s <- sim()
    s %>% select(time, `전신 면역억제 (χ 가중)` = Imm_sys_o,
                 `E_pred` = E_pred_o, `E_ifx` = E_ifx_o,
                 `E_vdz (장 선택적)` = E_vdz_o) %>%
      pivot_longer(-time) %>%
      ggplot(aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      ylim(0, 1) +
      labs(x = "일", y = "효과 분율", title = "면역억제의 구성",
           subtitle = "vedolizumab은 χ=0.10으로만 전신 계정에서 인출한다") +
      THEME
  })

  output$p_haz <- renderPlot({
    s <- sim()
    ggplot(s, aes(time)) +
      geom_line(aes(y = InfHaz, colour = "누적 감염 위험"), linewidth = 0.9) +
      geom_line(aes(y = SterCum / 20000, colour = "누적 스테로이드 / 20000"),
                linewidth = 0.9) +
      labs(x = "일", y = NULL, title = "면역억제의 대가") + THEME
  })

  ## ---- 8. scenario comparison --------------------------------------------
  cmp_sim <- reactive({
    inp <- reactiveValuesToList(input)
    specs <- list(
      "ipilimumab 1 mg/kg"  = list(ipi_mgkg = 1,  pd1_mgkg = 0),
      "ipilimumab 3 mg/kg"  = list(ipi_mgkg = 3,  pd1_mgkg = 0),
      "ipilimumab 10 mg/kg" = list(ipi_mgkg = 10, pd1_mgkg = 0),
      "anti-PD-1 3 mg/kg"   = list(ipi_mgkg = 0,  pd1_mgkg = 3),
      "anti-PD-1 10 mg/kg"  = list(ipi_mgkg = 0,  pd1_mgkg = 10),
      "ipi 3 + nivo 1"      = list(ipi_mgkg = 3,  pd1_mgkg = 1)
    )
    bind_rows(lapply(input$cmp, function(nm) {
      i2 <- modifyList(inp, specs[[nm]])
      run_sim(i2) %>% mutate(arm = nm)
    }))
  })

  output$p_cmp <- renderPlot({
    s <- cmp_sim(); req(nrow(s) > 0)
    ggplot(s, aes(time, S_eff_o, colour = arm)) + geom_line(linewidth = 0.9) +
      geom_hline(aes(yintercept = S_star_o), linetype = 2) + ylim(0, 1.02) +
      labs(x = "일", y = "유효 흡수 기능 S_eff",
           title = "시나리오 비교 — 역치를 넘는가, 넘지 않는가",
           subtitle = "anti-PD-1 3 vs 10 mg/kg는 겹친다. 그것이 요점이다") +
      THEME
  })

  output$t_cmp <- renderDT({
    s <- cmp_sim(); req(nrow(s) > 0)
    s %>% group_by(arm) %>%
      summarise(`최고 등급` = max(grade_o),
                `최저 S_eff` = round(min(S_eff_o), 3),
                `등급2 도달일` = suppressWarnings(min(time[grade_o >= 2])),
                `최고 칼프로텍틴` = round(max(Calpro)),
                `최저 Treg` = round(min(Treg), 3),
                `최고 G` = round(max(G_o), 2),
                `최고 ADCC 구동` = round(max(A_adcc_o), 3),
                `d180 종양` = round(last(Tumor), 1), .groups = "drop") %>%
      mutate(`등급2 도달일` = ifelse(is.infinite(`등급2 도달일`), NA,
                                     `등급2 도달일`)) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })
}

shinyApp(ui, server)
