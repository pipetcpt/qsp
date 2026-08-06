## =============================================================================
##  ctcl_shiny_app.R
##  Cutaneous T-cell lymphoma (mycosis fungoides / Sézary syndrome)
##  Interactive dashboard for the QSP model in ctcl_mrgsolve_model.R
##
##  The app is organised around the model's one claim: the clone lives in three
##  boxes, each therapy reaches a different subset of them, and each response
##  criterion reads only one of them.  Tab 3 (compartments) and tab 6 (global
##  response) are therefore the two tabs that matter; the rest support them.
##
##  Run with:  shiny::runApp("ctcl_shiny_app.R")
## =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(mrgsolve)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(DT)
})

source("ctcl_mrgsolve_model.R", local = TRUE)

THEME <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey35", size = 10))

PAL <- c(skin = "#c0392b", blood = "#1f6fb2", node = "#7f8c8d",
         resistant = "#e07a80", transformed = "#7b241c",
         surveillance = "#1f6fb2", drug = "#3d8b40")

REGIMENS <- c("없음 (natural history)", "협대역 UVB + 국소 스테로이드",
              "클로르메틴 젤", "저선량 TSEB 12 Gy", "벡사로텐 300 mg/m²",
              "인터페론 알파 3 MU tiw", "벡사로텐 + 인터페론",
              "보리노스타트 400 mg qd", "로미뎁신 14 mg/m²",
              "모가물리주맙 1 mg/kg", "모가물리주맙 3 mg/kg",
              "브렌툭시맙 베도틴 ×16", "로미뎁신 ×2 → 브렌툭시맙",
              "젬시타빈 ×6", "저용량 메토트렉세이트",
              "저용량 알렘투주맙 SC", "ECP q4w", "ECP + 인터페론",
              "모가물리주맙 + ECP", "항포도상구균 항생제 (d90–100)")

build_ev <- function(name, wt, bsa, dur) {
  switch(name,
    "없음 (natural history)"        = NULL,
    "협대역 UVB + 국소 스테로이드"  = c(ev_nbuvb(min(dur, 180)), ev_steroid(min(dur, 180))),
    "클로르메틴 젤"                 = ev_chlormethine(dur),
    "저선량 TSEB 12 Gy"             = ev_tseb(12, 8),
    "벡사로텐 300 mg/m²"            = ev_bexarotene(bsa, dur),
    "인터페론 알파 3 MU tiw"        = ev_interferon(dur),
    "벡사로텐 + 인터페론"           = c(ev_bexarotene(bsa, dur), ev_interferon(dur)),
    "보리노스타트 400 mg qd"        = ev_vorinostat(dur),
    "로미뎁신 14 mg/m²"             = ev_romidepsin(bsa, ceiling(dur / 28)),
    "모가물리주맙 1 mg/kg"          = ev_mogamulizumab(wt, dur),
    "모가물리주맙 3 mg/kg"          = c(ev(time = 0, amt = 3 * wt, cmt = "MOG1", ii = 7, addl = 3),
                                        ev(time = 28, amt = 3 * wt, cmt = "MOG1",
                                           ii = 14, addl = max(0, floor((dur - 28) / 14)))),
    "브렌툭시맙 베도틴 ×16"         = ev_brentuximab(wt, 16),
    "로미뎁신 ×2 → 브렌툭시맙"      = c(ev_romidepsin(bsa, 2),
                                        ev(time = 56, amt = 1.8 * wt, cmt = "BV1",
                                           ii = 21, addl = 12)),
    "젬시타빈 ×6"                   = ev_gemcitabine(bsa, 6),
    "저용량 메토트렉세이트"         = ev_methotrexate(dur),
    "저용량 알렘투주맙 SC"          = ev_alemtuzumab(min(dur, 180)),
    "ECP q4w"                       = ev_ecp(dur),
    "ECP + 인터페론"                = c(ev_ecp(dur), ev_interferon(dur)),
    "모가물리주맙 + ECP"            = c(ev_mogamulizumab(wt, dur), ev_ecp(dur)),
    "항포도상구균 항생제 (d90–100)" = ev_antibiotic(90, 10))
}

## =============================================================================
##  UI
## =============================================================================
ui <- fluidPage(
  titlePanel("피부 T세포 림프종 (MF / Sézary) QSP 대시보드 — 한 클론, 세 구획, 두 개의 숙주 변수"),
  tags$p(style = "color:#555;margin-top:-8px;",
         "Cutaneous T-cell lymphoma QSP dashboard · mSWAT는 피부만, B-score는 혈액만 읽습니다. ",
         "GLOBAL 반응은 그 둘의 논리곱입니다."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("환자 (Patient)"),
      selectInput("stage", "병기 (TNMB stage)",
                  c("IA", "IB", "IIA", "IIB", "IIB_LCT", "IIIA", "IIIB",
                    "IVA1", "IVA2", "IVB"), selected = "IIB"),
      radioButtons("kind", "표현형 (residency phenotype)",
                   c("MF (CD69 높음 · 피부 상주)" = "MF",
                     "Sézary (CD69 낮음 · 재순환)" = "SS"), selected = "MF"),
      sliderInput("cd69", "CD69 (조직 상주도)", 0.05, 0.97, 0.85, 0.01),
      sliderInput("ccr7", "CCR7 (재순환 프로그램)", 0.02, 0.98, 0.10, 0.01),
      sliderInput("ccr4", "CCR4 양성 분획", 0.10, 1.0, 0.85, 0.05),
      sliderInput("cd30", "CD30 발현", 0.0, 0.95, 0.10, 0.05),
      sliderInput("fres", "약물저항 클론 분획 FRES", 0.01, 0.90, 0.30, 0.01),
      sliderInput("sensf", "약물 감수성 배수", 0.2, 3.0, 1.0, 0.1),
      sliderInput("nkf", "NK 밀도 배수", 0.2, 2.5, 1.0, 0.1),
      sliderInput("deficit", "감시 결손 (질환 공격성)", 0.0, 0.45, 0.25, 0.01),
      numericInput("wt", "체중 (kg)", 75, 40, 140, 1),
      numericInput("bsa", "체표면적 (m²)", 1.85, 1.3, 2.5, 0.05),
      hr(),
      h4("치료 (Regimen)"),
      selectInput("reg", "요법 A", REGIMENS, selected = "모가물리주맙 1 mg/kg"),
      selectInput("reg2", "요법 B (비교)", REGIMENS, selected = "보리노스타트 400 mg qd"),
      sliderInput("dur", "투여 기간 (일)", 84, 730, 546, 7),
      sliderInput("end", "시뮬레이션 기간 (일)", 180, 1095, 730, 30),
      actionButton("go", "시뮬레이션 실행", class = "btn-primary btn-block")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("① 환자 프로파일",
                 br(), htmlOutput("profile"), br(),
                 plotOutput("p_profile", height = "330px"),
                 DT::dataTableOutput("t_profile")),
        tabPanel("② 약동학 (PK)",
                 br(), plotOutput("p_pk", height = "330px"),
                 h5("조직 노출과 수용체 점유율 — 피부에서 부족한 것은 약물이 아니다"),
                 plotOutput("p_occ", height = "280px")),
        tabPanel("③ 구획별 종양 부하",
                 br(), plotOutput("p_clone", height = "380px"),
                 h5("감수성 · 저항성 아클론"),
                 plotOutput("p_subclone", height = "280px")),
        tabPanel("④ 피부 종말점 (mSWAT)",
                 br(), plotOutput("p_mswat", height = "340px"),
                 h5("반점 · 판 · 종양의 서로 다른 소실 시간상수"),
                 plotOutput("p_morph", height = "280px")),
        tabPanel("⑤ 혈액 · 림프절",
                 br(), plotOutput("p_blood", height = "340px"),
                 plotOutput("p_node", height = "250px")),
        tabPanel("⑥ GLOBAL 반응 (Olsen 논리곱)",
                 br(), htmlOutput("respbox"),
                 plotOutput("p_resp", height = "300px"),
                 DT::dataTableOutput("t_resp")),
        tabPanel("⑦ 감시 변수 · 지속성",
                 br(),
                 tags$p("ORR은 '얼마나 죽였나'를 읽고, TTNT는 '누가 남아서 감시하는가'가 정합니다."),
                 plotOutput("p_surv", height = "340px"),
                 plotOutput("p_ttnt", height = "260px")),
        tabPanel("⑧ 소양증 · 삶의 질",
                 br(), plotOutput("p_itch", height = "330px"),
                 h5("종양 부하와 소양증의 분리"),
                 plotOutput("p_itchdecouple", height = "280px")),
        tabPanel("⑨ 혈청 바이오마커",
                 br(), plotOutput("p_biomk", height = "420px")),
        tabPanel("⑩ 장벽 · 포도상구균 · 감염",
                 br(), plotOutput("p_barrier", height = "340px"),
                 plotOutput("p_infect", height = "260px")),
        tabPanel("⑪ 독성 · 안전성",
                 br(), plotOutput("p_tox", height = "420px"),
                 DT::dataTableOutput("t_tox")),
        tabPanel("⑫ 시나리오 비교",
                 br(), plotOutput("p_cmp", height = "380px"),
                 DT::dataTableOutput("t_cmp")),
        tabPanel("⑬ 구조 실험",
                 br(),
                 tags$p(strong("E1"), " — 모가물리주맙의 구획 격차는 노출이 아니라 효과세포 때문이다."),
                 DT::dataTableOutput("t_e1"), br(),
                 tags$p(strong("E6"), " — CD69 하나만 움직여 MF에서 Sézary로."),
                 plotOutput("p_e6", height = "280px"))
      )
    )
  )
)

## =============================================================================
##  SERVER
## =============================================================================
server <- function(input, output, session) {

  observeEvent(input$kind, {
    if (input$kind == "MF") {
      updateSliderInput(session, "cd69", value = 0.85)
      updateSliderInput(session, "ccr7", value = 0.10)
    } else {
      updateSliderInput(session, "cd69", value = 0.20)
      updateSliderInput(session, "ccr7", value = 0.85)
    }
  })

  covs <- reactive(list(CD69 = input$cd69, CCR7 = input$ccr7, CCR4 = input$ccr4,
                        CD30 = input$cd30, FRES = input$fres, SENSF = input$sensf,
                        NKF = input$nkf))

  patient <- eventReactive(input$go, {
    withProgress(message = "환자 평형화 중 (stage = 감시 스위치의 평형점)", value = 0.3, {
      ctcl_patient(input$stage, input$kind, input$deficit, covs())
    })
  }, ignoreNULL = FALSE)

  simA <- eventReactive(input$go, {
    pt <- patient(); e <- build_ev(input$reg, input$wt, input$bsa, input$dur)
    d <- if (is.null(e)) mrgsim(pt$mod, end = input$end, delta = 1)
         else mrgsim(pt$mod, events = e, end = input$end, delta = 1)
    as.data.frame(d)
  }, ignoreNULL = FALSE)

  simB <- eventReactive(input$go, {
    pt <- patient(); e <- build_ev(input$reg2, input$wt, input$bsa, input$dur)
    d <- if (is.null(e)) mrgsim(pt$mod, end = input$end, delta = 1)
         else mrgsim(pt$mod, events = e, end = input$end, delta = 1)
    as.data.frame(d)
  }, ignoreNULL = FALSE)

  both <- reactive(bind_rows(mutate(simA(), arm = input$reg),
                             mutate(simB(), arm = input$reg2)))

  ## ---- ① profile ----------------------------------------------------------
  output$profile <- renderUI({
    pt <- patient(); x <- simA()
    HTML(sprintf(
      "<div style='background:#f6f6f9;padding:10px;border-radius:6px'>
       <b>병기</b> %s · <b>표현형</b> %s ·
       <b>평형 면역적합도 IMMF<sub>eq</sub></b> %.2f → <b>적용값</b> %.2f
       (감시 결손 %.0f%%)<br>
       <b>기저</b> mSWAT %.1f · Sézary %.0f/µL (B%d) · N-score %d ·
       소양증 NRS %.1f · TARC %.0f pg/mL<br>
       <span style='color:#666'>병기는 초기값이 아니라 감시 스위치의 평형점으로 계산되며,
       피부·혈액 분배는 CD69/CCR7로부터 <i>유도</i>됩니다.</span></div>",
      pt$stage, pt$kind, pt$immf_eq, pt$immf, 100 * input$deficit,
      x$MSWAT[1], x$SEZCT[1], as.integer(x$BSCORE[1]), as.integer(x$NSCORE[1]),
      x$PRUR[1], x$TARC[1]))
  })

  output$p_profile <- renderPlot({
    x <- simA()[1, ]
    d <- data.frame(
      compartment = factor(c("피부 감수성", "피부 저항성", "전환 아클론",
                             "혈액 감수성", "혈액 저항성", "림프절", "내장"),
                           levels = c("피부 감수성", "피부 저항성", "전환 아클론",
                                      "혈액 감수성", "혈액 저항성", "림프절", "내장")),
      cells = c(x$NSK, x$NSKR, x$NTR, x$NBL, x$NBLR, x$NLN, x$NVS))
    ggplot(d, aes(compartment, cells, fill = compartment)) +
      geom_col(width = 0.65, show.legend = FALSE) +
      labs(title = "기저 클론 분포 (10⁹ cells)",
           subtitle = "MF/SS의 차이는 방정식이 아니라 CD69·CCR7 값에서 나옵니다",
           x = NULL, y = "10⁹ cells") + THEME
  })

  output$t_profile <- DT::renderDataTable({
    x <- simA()[1, ]
    DT::datatable(data.frame(
      항목 = c("mSWAT", "반점 %BSA", "판 %BSA", "종양 %BSA", "Sézary /µL",
               "B-score", "N-score", "M-score", "T-stage", "소양증 NRS",
               "TARC pg/mL", "sIL-2R U/mL", "LDH U/L", "장벽 무결성", "CD4 /µL"),
      값 = round(c(x$MSWAT, x$APAT, x$APLQ, x$ATUM, x$SEZCT, x$BSCORE, x$NSCORE,
                   x$MSCORE, x$TSKIN, x$PRUR, x$TARC, x$SIL2R, x$LDH, x$BARR,
                   x$CD4N), 2)),
      options = list(dom = "t", pageLength = 20), rownames = FALSE)
  })

  ## ---- ② PK ---------------------------------------------------------------
  output$p_pk <- renderPlot({
    x <- simA()
    d <- x %>% select(time, `모가물리주맙 혈장 (µg/mL)` = CMOGP,
                      `모가물리주맙 피부 (µg/mL)` = CMOGSKO,
                      `BV ADC 혈장 (µg/mL)` = CBVP,
                      `유리 MMAE (ng/mL)` = CMMP) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(colour = PAL[["drug"]], linewidth = 0.7) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "약동학", subtitle = "피부 농도 = 0.157 × 혈장 (Shah & Betts 2013, 적합하지 않은 상수)",
           x = "일 (day)", y = NULL) + THEME
  })

  output$p_occ <- renderPlot({
    x <- simA()
    d <- x %>% select(time, `혈액 CCR4 점유율` = OCCBL, `피부 CCR4 점유율` = OCCSK) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = "CCR4 수용체 점유율", x = "일", y = "점유율", colour = NULL,
           subtitle = "피부 점유율도 거의 1 — 용량을 올려도 피부 반응이 움직이지 않는 이유") + THEME
  })

  ## ---- ③ compartments -----------------------------------------------------
  output$p_clone <- renderPlot({
    d <- both() %>%
      transmute(time, arm, 피부 = NSK + NSKR + NTR, 혈액 = NBL + NBLR,
                림프절 = NLN) %>% pivot_longer(c(피부, 혈액, 림프절))
    ggplot(d, aes(time, value, colour = name, linetype = arm)) +
      geom_line(linewidth = 0.8) + scale_y_log10() +
      scale_colour_manual(values = c(피부 = PAL[["skin"]], 혈액 = PAL[["blood"]],
                                     림프절 = PAL[["node"]])) +
      labs(title = "구획별 클론 부하", subtitle = "같은 약이 구획마다 다른 속도로 듣습니다",
           x = "일", y = "10⁹ cells (log)", colour = NULL, linetype = NULL) + THEME
  })

  output$p_subclone <- renderPlot({
    d <- simA() %>% transmute(time, 감수성 = NSK, 저항성 = NSKR, 전환 = NTR) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, fill = name)) + geom_area(alpha = 0.85) +
      scale_fill_manual(values = c(감수성 = PAL[["skin"]], 저항성 = PAL[["resistant"]],
                                   전환 = PAL[["transformed"]])) +
      labs(title = "피부 클론의 아집단", x = "일", y = "10⁹ cells", fill = NULL,
           subtitle = "저항성 아클론은 약에는 저항하지만 CD8 감시에는 저항하지 않습니다") + THEME
  })

  ## ---- ④ mSWAT ------------------------------------------------------------
  output$p_mswat <- renderPlot({
    d <- both()
    b <- d %>% group_by(arm) %>% slice(1) %>% transmute(arm, base = MSWAT)
    d <- left_join(d, b, by = "arm")
    ggplot(d, aes(time, MSWAT, colour = arm)) + geom_line(linewidth = 0.9) +
      geom_hline(aes(yintercept = 0.5 * base, colour = arm), linetype = "dashed",
                 alpha = 0.5) +
      labs(title = "mSWAT (피부 구획 점수)", x = "일", y = "mSWAT", colour = NULL,
           subtitle = "점선 = 각 팔의 50% 감소선 (피부 부분반응 기준)") + THEME
  })

  output$p_morph <- renderPlot({
    d <- simA() %>% transmute(time, 반점 = APAT, 판 = APLQ, 종양 = ATUM) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(title = "병변 형태별 면적", x = "일", y = "%BSA", colour = NULL,
           subtitle = "τ: 반점 8일 · 판 25일 · 종양 12일 — 반응 확인에 4주가 필요한 이유") + THEME
  })

  ## ---- ⑤ blood / node -----------------------------------------------------
  output$p_blood <- renderPlot({
    ggplot(both(), aes(time, SEZCT, colour = arm)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = c(250, 1000), linetype = "dotted", colour = "grey40") +
      scale_y_log10() +
      labs(title = "혈중 Sézary 세포수", x = "일", y = "cells/µL (log)", colour = NULL,
           subtitle = "점선 = B1 (250) 및 B2 (1000) 경계") + THEME
  })
  output$p_node <- renderPlot({
    ggplot(both(), aes(time, NLN, colour = arm)) + geom_line(linewidth = 0.9) +
      labs(title = "림프절 클론 부하", x = "일", y = "10⁹ cells", colour = NULL) + THEME
  })

  ## ---- ⑥ global response --------------------------------------------------
  respA <- reactive(ctcl_endpoints(simA()))
  respB <- reactive(ctcl_endpoints(simB()))

  output$respbox <- renderUI({
    a <- respA()
    HTML(sprintf(
      "<div style='background:#f3eefa;padding:10px;border-radius:6px'>
       <b>%s</b> — 최고 mSWAT 감소 <b>%.1f%%</b> ·
       피부반응일 %s · 혈액반응일 %s · <b>GLOBAL 반응일 %s</b> ·
       PFS %s일 · TTNT %s일<br>
       <span style='color:#666'>GLOBAL은 침범한 모든 구획에서 동시에 반응이 성립해야 합니다.
       한 구획만 좋아도 GLOBAL은 움직이지 않습니다.</span></div>",
      input$reg, a$best_skin_pct,
      ifelse(is.na(a$skin_response_day), "—", a$skin_response_day),
      ifelse(is.na(a$blood_response_day), "—", a$blood_response_day),
      ifelse(is.na(a$global_response_day), "—", a$global_response_day),
      ifelse(is.na(a$PFS_day), "미도달", a$PFS_day),
      ifelse(is.na(a$TTNT_day), "미도달", a$TTNT_day)))
  })

  output$p_resp <- renderPlot({
    x <- simA()
    b <- c(x$MSWAT[1], x$SEZCT[1], x$NLN[1])
    d <- data.frame(time = x$time,
                    피부 = 100 * (1 - x$MSWAT / b[1]),
                    혈액 = 100 * (1 - x$SEZCT / max(b[2], 1e-9)),
                    림프절 = 100 * (1 - x$NLN / max(b[3], 1e-9))) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 50, linetype = "dashed", colour = "grey40") +
      scale_colour_manual(values = c(피부 = PAL[["skin"]], 혈액 = PAL[["blood"]],
                                     림프절 = PAL[["node"]])) +
      labs(title = "구획별 기저 대비 감소율", x = "일", y = "% 감소", colour = NULL,
           subtitle = "점선 = 50% (부분반응). GLOBAL 반응은 세 선이 모두 위에 있어야 성립") + THEME
  })

  output$t_resp <- DT::renderDataTable({
    DT::datatable(bind_rows(cbind(arm = input$reg, respA()),
                            cbind(arm = input$reg2, respB())),
                  options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  ## ---- ⑦ surveillance -----------------------------------------------------
  output$p_surv <- renderPlot({
    d <- both() %>% group_by(arm) %>%
      transmute(time, `피부 CD8 감시` = SURVSK / first(SURVSK),
                `조절 T세포` = TREG / first(TREG),
                `NK 세포` = NKB / first(NKB),
                `DC 활성` = DCA / first(DCA)) %>% ungroup() %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = name, linetype = arm)) +
      geom_line(linewidth = 0.8) + geom_hline(yintercept = 1, colour = "grey70") +
      labs(title = "두 번째 상태변수 — 감시 (기저 대비)", x = "일", y = "배수",
           colour = NULL, linetype = NULL,
           subtitle = "화학요법은 클론과 감시를 함께 지웁니다. ECP·인터페론은 감시만 올립니다.") + THEME
  })

  output$p_ttnt <- renderPlot({
    d <- both() %>% group_by(arm) %>%
      transmute(time, mSWAT = MSWAT / first(MSWAT) * 100) %>% ungroup()
    ggplot(d, aes(time, mSWAT, colour = arm)) + geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 75, linetype = "dashed", colour = "grey40") +
      labs(title = "다음 치료까지의 기간 (TTNT) 기준선", x = "일",
           y = "기저 대비 mSWAT (%)", colour = NULL,
           subtitle = "점선 75% 재도달 시점 = TTNT") + THEME
  })

  ## ---- ⑧ pruritus ---------------------------------------------------------
  output$p_itch <- renderPlot({
    ggplot(both(), aes(time, PRUR, colour = arm)) + geom_line(linewidth = 0.9) +
      scale_y_continuous(limits = c(0, 10)) +
      labs(title = "소양증 NRS", x = "일", y = "NRS 0–10", colour = NULL,
           subtitle = "IL-31 · Th2 · 장벽 · 초항원의 합. 중추 감작(SENS)이 느린 꼬리를 만듭니다.") + THEME
  })
  output$p_itchdecouple <- renderPlot({
    x <- simA()
    d <- data.frame(time = x$time,
                    `종양 부하 (기저=1)` = (x$NSK + x$NSKR + x$NTR) / (x$NSK[1] + x$NSKR[1] + x$NTR[1]),
                    `소양증 (기저=1)` = x$PRUR / x$PRUR[1],
                    `IL-31 (기저=1)` = x$IL31 / x$IL31[1],
                    check.names = FALSE) %>% pivot_longer(-time)
    ggplot(d, aes(time, value, colour = name)) + geom_line(linewidth = 0.8) +
      labs(title = "소양증은 종양 부하를 그대로 따라가지 않는다", x = "일",
           y = "기저 대비", colour = NULL) + THEME
  })

  ## ---- ⑨ biomarkers -------------------------------------------------------
  output$p_biomk <- renderPlot({
    d <- both() %>% select(time, arm, `TARC (pg/mL)` = TARC, `sIL-2R (U/mL)` = SIL2R,
                           `LDH (U/L)` = LDH, `Th2 톤` = TH2, `IL-31` = IL31,
                           `IFN-γ` = IFNG) %>% pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "혈청 바이오마커와 사이토카인", x = "일", y = NULL, colour = NULL) + THEME
  })

  ## ---- ⑩ barrier ----------------------------------------------------------
  output$p_barrier <- renderPlot({
    d <- both() %>% select(time, arm, `장벽 무결성` = BARR,
                           `황색포도상구균` = SAUR, `초항원` = SAG) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "장벽 – 포도상구균 – 초항원 되먹임", x = "일", y = NULL, colour = NULL,
           subtitle = "항생제는 클론을 죽이지 않고도 이 고리를 끊습니다") + THEME
  })
  output$p_infect <- renderPlot({
    d <- both() %>% select(time, arm, `누적 감염 위험` = HZI,
                           `누적 질환 위험` = HZD, `생존 확률` = SURV) %>%
      pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "위험 적분과 생존", x = "일", y = NULL, colour = NULL) + THEME
  })

  ## ---- ⑪ toxicity ---------------------------------------------------------
  output$p_tox <- renderPlot({
    d <- both() %>% select(time, arm, `호중구 (10⁹/L)` = CIRC, `혈소판 (10⁹/L)` = PLT,
                           `말초신경병증` = PN, `CD4 (/µL)` = CD4N,
                           `유리 T4 (ng/dL)` = FT4, `중성지방 (mg/dL)` = TG,
                           `모가 연관 발진` = MAR) %>% pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.8) +
      facet_wrap(~name, scales = "free_y") +
      labs(title = "독성 지표", x = "일", y = NULL, colour = NULL) + THEME
  })
  output$t_tox <- DT::renderDataTable({
    f <- function(x, nm) data.frame(
      arm = nm, ANC최저 = round(min(x$CIRC), 2), 혈소판최저 = round(min(x$PLT), 0),
      신경병증최대 = round(max(x$PN), 2), CD4최저 = round(min(x$CD4N), 0),
      FT4최저 = round(min(x$FT4), 2), 중성지방최대 = round(max(x$TG), 0),
      발진최대 = round(max(x$MAR), 2))
    DT::datatable(bind_rows(f(simA(), input$reg), f(simB(), input$reg2)),
                  options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- ⑫ comparison -------------------------------------------------------
  output$p_cmp <- renderPlot({
    d <- both() %>% group_by(arm) %>%
      transmute(time, mSWAT = MSWAT / first(MSWAT),
                Sezary = SEZCT / max(first(SEZCT), 1e-9),
                감시 = SURVSK / first(SURVSK), NRS = PRUR / first(PRUR)) %>%
      ungroup() %>% pivot_longer(-c(time, arm))
    ggplot(d, aes(time, value, colour = arm)) + geom_line(linewidth = 0.85) +
      facet_wrap(~name, scales = "free_y") + geom_hline(yintercept = 1, colour = "grey80") +
      labs(title = "요법 A vs B (모두 기저 = 1)", x = "일", y = "기저 대비", colour = NULL) + THEME
  })
  output$t_cmp <- DT::renderDataTable({
    DT::datatable(bind_rows(cbind(arm = input$reg, respA()),
                            cbind(arm = input$reg2, respB())),
                  options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  ## ---- ⑬ structural experiments -------------------------------------------
  output$t_e1 <- DT::renderDataTable({
    DT::datatable(E1_effector_not_exposure(), options = list(dom = "t"), rownames = FALSE)
  })
  output$p_e6 <- renderPlot({
    d <- E6_residency_sweep()
    ggplot(d, aes(CD69, Sezary_per_uL)) +
      geom_line(colour = PAL[["blood"]], linewidth = 1) + geom_point(size = 2) +
      geom_hline(yintercept = c(250, 1000), linetype = "dotted") +
      scale_y_log10() + scale_x_reverse() +
      labs(title = "CD69 하나로 MF에서 Sézary까지",
           subtitle = "같은 병기, 같은 방정식 — 조직상주도만 바꿉니다",
           x = "CD69 (조직 상주도, 오른쪽이 낮음)", y = "Sézary /µL (log)") + THEME
  })
}

shinyApp(ui, server)
