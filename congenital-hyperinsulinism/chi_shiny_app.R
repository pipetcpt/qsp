## =====================================================================
##  Congenital Hyperinsulinism (CHI) — QSP explorer (Shiny)
##  ---------------------------------------------------------------------
##  The app is built around ONE control: the slider for g, the fraction of
##  beta-cell K_ATP conductance that survives.  Every tab is a different
##  way of reading the consequences of that one number.
##
##  Run:  shiny::runApp("chi_shiny_app.R")
##  Requires: shiny, mrgsolve, dplyr, tidyr, ggplot2, DT
##  The model code is sourced from chi_mrgsolve_model.R (same directory).
## =====================================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

## ---------------------------------------------------------------------
## Load the compiled model + helpers from the model file.  We source it in
## a throw-away environment and lift out what we need, so that running the
## model file's own scenarios does not happen here.
## ---------------------------------------------------------------------
.src <- new.env()
local({
  code <- readLines("chi_mrgsolve_model.R")
  ## keep everything up to (but not including) the scenario section
  cut <- grep("^##  SCENARIO 1", code)
  if (length(cut)) code <- code[seq_len(cut[1] - 3)]
  eval(parse(text = paste(code, collapse = "\n")), envir = .src)
})
mod  <- .src$mod
geno <- .src$geno

FEED_MGKGMIN <- 845 * 8 / 1440    # 4.69 mg/kg/min enteral baseline

GENO_LABEL <- c(
  "정상 신생아 (normal)"                              = "normal",
  "K_ATP 우성 (ABCC8 dominant, g≈0.60)"               = "katp_dominant",
  "K_ATP 열성 미만성 (biallelic diffuse, g≈0.02)"      = "katp_recess_diff",
  "K_ATP 국소형 (focal, paternal UPD 11p15)"           = "katp_focal",
  "GDH-HI / HI-HA (GLUD1 활성화)"                     = "gdh_hi",
  "GCK 활성화 (set-point 이동)"                        = "gck_activating",
  "HNF4A (일시적, 후기 MODY)"                          = "hnf4a",
  "SCHAD/HADH 결핍"                                   = "schad",
  "SLC16A1 (운동유발성)"                               = "slc16a1",
  "췌장 아전절제 후 (remnant 5%)"                      = "post_pancreatect"
)

feeds     <- function(days = 4)  ev(amt = 845, cmt = "GGUT", time = 0,
                                   ii = 3, addl = days * 8 - 1)
diazoxide <- function(d, days = 6) ev(amt = d/3, cmt = "DZXg", time = 0,
                                   ii = 8, addl = days * 3 - 1)
octreotide<- function(d, days = 6) ev(amt = d/4, cmt = "OCTs", time = 0,
                                   ii = 6, addl = days * 4 - 1)
ersodetug <- function(d)          ev(amt = d/0.055, cmt = "ERS", time = 0)
nifedipine<- function(d, days=6)  ev(amt = d/1.2, cmt = "NIF", time = 0,
                                   ii = 8, addl = days * 3 - 1)
sirolimus <- function(d, days=10) ev(amt = d, cmt = "SIRg", time = 0,
                                   ii = 24, addl = days - 1)

## ---------------------------------------------------------------------
## one simulation, assembled from the UI state
## ---------------------------------------------------------------------
simulate <- function(gt, g_override = NA, dzx = 0, oct = 0, ers = 0, nif = 0,
                     sir = 0, gcg = 0, remnant = 1, fed = TRUE,
                     closed_loop = TRUE, target = 70, end = 72,
                     fast_from = NA) {
  p <- geno[[gt]]
  if (!is.na(g_override) && gt != "normal") p$g_ab <- g_override
  p$BMASS0 <- remnant
  p$GCG_inf <- gcg
  if (closed_loop) { p$loop <- 1; p$Gtarget <- target }

  e <- NULL
  nfeed <- if (is.na(fast_from)) ceiling(end/24) + 1 else max(1, floor(fast_from/24) + 1)
  if (fed) {
    nf <- if (is.na(fast_from)) nfeed * 8 else max(1, floor(fast_from/3))
    e <- ev(amt = 845, cmt = "GGUT", time = 0, ii = 3, addl = nf - 1)
  }
  add <- function(x) if (is.null(e)) x else c(e, x)
  if (dzx > 0) e <- add(diazoxide(dzx, ceiling(end/24) + 1))
  if (oct > 0) e <- add(octreotide(oct, ceiling(end/24) + 1))
  if (ers > 0) e <- add(ersodetug(ers))
  if (nif > 0) e <- add(nifedipine(nif, ceiling(end/24) + 1))
  if (sir > 0) e <- add(sirolimus(sir, ceiling(end/24) + 1))

  out <- mod %>% param(p)
  if (!is.null(e)) out <- out %>% ev(e)
  out %>% mrgsim(end = end, delta = 0.05) %>% as_tibble()
}

tail_stats <- function(d, frac = 0.25) {
  tmax <- max(d$time); d <- filter(d, time > tmax * (1 - frac))
  tibble(
    glucose_mean = mean(d$GLU), glucose_min = min(d$GLU),
    insulin      = mean(d$INS), cpeptide = mean(d$CPEP),
    GIR          = mean(d$GIRout), total_glc = mean(d$TOTGLC),
    BOHB         = mean(d$BOHB), FFA = mean(d$FFA),
    ammonia      = mean(d$NH3), glycogen = mean(d$GLY),
    fuel_ratio   = mean(d$fuelrat), Vm = mean(d$Vm_ab),
    pct_below70  = 100 * mean(d$GLU < 70), pct_below54 = 100 * mean(d$GLU < 54)
  )
}

THEME <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom", plot.title = element_text(face = "bold"))

## =====================================================================
## UI
## =====================================================================
ui <- fluidPage(
  titlePanel("선천성 고인슐린증 (Congenital Hyperinsulinism) — QSP 탐색기"),
  tags$p(style = "color:#555;margin-top:-8px",
    HTML("이 앱의 모든 결과는 <b>단 하나의 숫자</b>, 즉 잔존 베타세포 K<sub>ATP</sub> 전도도 <b>g</b>의 산술적 결과입니다. ",
         "디아족사이드는 g에 <b>곱해지고</b>, 옥트레오타이드는 전압 분배기에 <b>더해지며</b>, ",
         "글루카곤·포도당·에르소데투그·수술은 g와 <b>무관</b>합니다.")),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("gt", "유전형 (genotype)", choices = GENO_LABEL,
                  selected = "katp_recess_diff"),
      sliderInput("g", HTML("잔존 K<sub>ATP</sub> 전도도 <b>g</b>"),
                  min = 0, max = 1, value = 0.02, step = 0.01),
      helpText(HTML("<small>g≈0–0.05 열성 미만성 · 0.5–0.8 우성 · 1.0 정상<br>",
                    "g&lt;0.1에서는 <b>중증도가 포화</b>됩니다 (막전위가 이미 leak 전위에 고정).</small>")),
      hr(),
      h5("치료 (Therapy)"),
      sliderInput("dzx", "디아족사이드 (mg/kg/day)", 0, 20, 0, step = 2.5),
      sliderInput("oct", "옥트레오타이드 (µg/kg/day)", 0, 40, 0, step = 5),
      sliderInput("gcg", "글루카곤 지속주입 (µg/kg/h)", 0, 20, 0, step = 2.5),
      sliderInput("ers", "에르소데투그 (mg/kg IV)", 0, 12, 0, step = 3),
      sliderInput("nif", "니페디핀 (mg/kg q8h)", 0, 1, 0, step = 0.25),
      sliderInput("sir", "시롤리무스 (mg/kg/day)", 0, 0.12, 0, step = 0.02),
      hr(),
      h5("수술 / 관리"),
      sliderInput("remnant", "베타세포 잔존량 (수술 후)", 0.02, 1.0, 1.0, step = 0.02),
      checkboxInput("loop", "정맥 포도당 자동적정 (closed loop)", TRUE),
      sliderInput("target", "목표 혈당 (mg/dL)", 50, 100, 70, step = 5),
      checkboxInput("fed", "경장수유 q3h", TRUE),
      sliderInput("end", "시뮬레이션 기간 (h)", 24, 168, 72, step = 12)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("① 환자 프로파일",
          br(), fluidRow(column(12, DTOutput("profile"))),
          br(), plotOutput("betacell", height = "330px"),
          helpText(HTML("<small>왼쪽: 전압 분배기. K<sub>ATP</sub> 전도도가 leak을 압도하면 세포는 침묵하고, ",
                        "사라지면 −20 mV에 고정되어 <b>모든 혈당에서</b> 인슐린을 분비합니다.</small>"))),
        tabPanel("② PK (약물 농도)",
          br(), plotOutput("pk", height = "560px"),
          helpText(HTML("<small>디아족사이드 t½ ≈ 20 h (신생아에서 최대 30 h) — 정상상태 도달에 3–5일. ",
                        "옥트레오타이드 t½ ≈ 1.7 h이므로 q6h 투여에서 큰 진폭을 보입니다.</small>"))),
        tabPanel("③ PD 주요지표",
          br(), plotOutput("pd", height = "620px")),
        tabPanel("④ 임상 엔드포인트",
          br(), fluidRow(column(12, DTOutput("endpoints"))),
          br(), plotOutput("gir", height = "330px")),
        tabPanel("⑤ 시나리오 비교",
          br(),
          helpText("동일 유전형에서 각 약제를 단독으로 적용해 정맥 포도당 요구량을 비교합니다."),
          plotOutput("compare", height = "380px"),
          br(), DTOutput("compare_tbl")),
        tabPanel("⑥ 바이오마커",
          br(), plotOutput("biomarker", height = "620px"),
          helpText(HTML("<small>케톤 부재는 CHI의 <b>진단적 핵심</b>입니다. 지방분해 IC₅₀(12) &lt; 케톤생성(15) ",
                        "&lt; 당분해(30) &lt; 당신생(45) µU/mL이므로, 저혈당을 일으킬 만한 인슐린은 ",
                        "이미 케톤을 껐습니다. 암모니아는 GDH-HI에서만 상승하고 디아족사이드에 반응하지 않습니다.</small>"))),
        tabPanel("⑦ 뇌 연료 (핵심)",
          br(), plotOutput("fuel", height = "400px"),
          br(), plotOutput("fuelthresh", height = "330px"),
          helpText(HTML("<small>모델이 계산한 값: 케톤이 0일 때 뇌는 <b>43.7 mg/dL</b>에서 연료가 고갈되지만, ",
                        "BOHB 2 mM에서는 <b>25.3 mg/dL</b>까지 견딥니다. 고인슐린증은 포도당을 빼앗기 전에 ",
                        "이미 <b>18.4 mg/dL의 안전역</b>을 삭제합니다 — 케톤성 저혈당은 40대까지 허용되지만 ",
                        "CHI 목표가 70 mg/dL인 이유입니다.</small>"))),
        tabPanel("⑧ g-축 탐색",
          br(),
          helpText(HTML("모델의 중심 주장: <b>같은 엔드포인트, 같은 용량, 정반대의 결론.</b>")),
          plotOutput("gsweep", height = "420px"),
          br(), DTOutput("gsweep_tbl")),
        tabPanel("⑨ 수술 계획",
          br(), plotOutput("surgery", height = "380px"),
          br(), DTOutput("surgery_tbl"),
          helpText(HTML("<small>수술은 <b>B(질량)</b>를 줄이고 <b>g</b>는 그대로 둡니다. 그래서 하나의 방정식이 ",
                        "지속 저혈당과 수술 후 당뇨를 모두 만듭니다. 잔존량 0.20–0.30이 유일한 창(window)이고, ",
                        "≤0.10은 즉시 당뇨입니다. 게다가 성장이 같은 잔존량을 다시 읽습니다: 3.5 kg에서 충분한 ",
                        "2 % 잔존은 20 kg에서 정상의 0.35 %에 불과합니다.</small>"))),
        tabPanel("⑩ 문서",
          br(),
          h4("모델이 무엇을 보정했고 무엇을 예측했는가"),
          HTML("<p>보정에 쓴 숫자는 <b>9개</b>이고 그중 8개는 <b>정상 신생아 생리</b>입니다:
                혈당 75 mg/dL에서 인슐린 5 µU/mL, 정상 포도당→인슐린 용량반응의 <i>모양</i>
                (75 mg/dL 대비 45에서 0.25배·90에서 1.8배·150에서 12배·250에서 24배),
                신생아 뇌 포도당 소비 4.2 mg/kg/min, 만삭 간 글리코겐 약 2000 mg/kg,
                그리고 인슐린 IC₅₀의 서열입니다. CHI에서 가져온 숫자는
                <b>단 하나</b>(gGIRK=1.2, 옥트레오타이드가 중증 미만성 CHI의 포도당 요구량을
                약 절반으로 줄인다는 관찰)뿐입니다.</p>
                <p>따라서 아래는 모두 <b>예측</b>입니다: 중증 미만성 CHI의 총 포도당 요구량
                15.9 mg/kg/min · 디아족사이드 반응률의 유전형 분리(g=0.02에서 0.1%, g≥0.3에서 100%) ·
                g&lt;0.1에서의 중증도 포화 · 케톤 부재 · 글루카곤 자극검사 양성(+62 mg/dL) ·
                절제 범위의 좁은 창 · GDH-HI의 류신 감수성과 디아족사이드에 반응하지 않는 암모니아 ·
                니페디핀의 실패(기전은 옳지만 EC₅₀에 도달할 수 없음).</p>
                <h4>정직하게 밝히는 긴장</h4>
                <p>예측된 혈장 인슐린(≈83 µU/mL)은 보고된 critical sample 값(10–50 µU/mL)보다
                <b>높습니다</b>. 이것은 선택이 아니라 강제된 산술입니다 — 신생아 인슐린 감수성에서
                16 mg/kg/min의 포도당 요구량을 20 µU/mL로 설명할 수는 없습니다. 이 모델에서 가장
                노출된 파라미터이며, README에 T1–T4로 모두 기재했습니다.</p>
                <p style='color:#a33'><b>면책:</b> 교육·연구 목적의 QSP 모델입니다.
                임상적 의사결정이나 투약에 사용하지 마십시오.</p>")
        )
      )
    )
  )
)

## =====================================================================
## SERVER
## =====================================================================
server <- function(input, output, session) {

  observeEvent(input$gt, {
    gdef <- geno[[input$gt]]$g_ab
    if (!is.null(gdef)) updateSliderInput(session, "g", value = gdef)
    rem <- geno[[input$gt]]$BMASS0
    updateSliderInput(session, "remnant", value = if (is.null(rem)) 1 else rem)
  })

  sim <- reactive({
    simulate(input$gt, g_override = input$g, dzx = input$dzx, oct = input$oct,
             ers = input$ers, nif = input$nif, sir = input$sir, gcg = input$gcg,
             remnant = input$remnant, fed = input$fed,
             closed_loop = input$loop, target = input$target, end = input$end)
  })

  ## ---------------- ① profile ----------------
  output$profile <- renderDT({
    s <- tail_stats(sim())
    datatable(tibble(
      항목 = c("잔존 K_ATP 전도도 g", "베타세포 막전위 (mV)",
               "평균 혈당 (mg/dL)", "최저 혈당 (mg/dL)",
               "혈장 인슐린 (µU/mL)", "C-펩타이드 (ng/mL)",
               "정맥 포도당 (mg/kg/min)", "총 포도당 공급 (mg/kg/min)",
               "간 글리코겐 (mg/kg)", "뇌 연료 공급/요구 비"),
      값 = c(sprintf("%.2f", input$g), sprintf("%.1f", s$Vm),
             sprintf("%.1f", s$glucose_mean), sprintf("%.1f", s$glucose_min),
             sprintf("%.1f", s$insulin), sprintf("%.2f", s$cpeptide),
             sprintf("%.2f", s$GIR), sprintf("%.2f", s$total_glc),
             sprintf("%.0f", s$glycogen), sprintf("%.3f", s$fuel_ratio))
    ), options = list(dom = "t", pageLength = 12), rownames = FALSE)
  })

  output$betacell <- renderPlot({
    gg <- seq(0, 1, 0.01)
    p  <- as.list(param(mod))
    Gs <- c(45, 60, 90, 150)
    d <- expand.grid(g = gg, G = Gs) %>% as_tibble() %>% rowwise() %>%
      mutate({
        Gb <- G/18.016
        Rt <- p$Rbas + (p$Rmax - p$Rbas)*Gb^p$hG/(p$KG^p$hG + Gb^p$hG) +
              p$kAA*(p$AA0/(p$KAA + p$AA0))
        Po <- 1/(1 + (Rt/p$R50)^p$nR)
        GK <- p$gKmax*g*Po
        tibble(Vm = ((GK)*p$EK + p$gleak*p$Eleak)/(GK + p$gleak))
      }) %>% ungroup()
    ggplot(d, aes(g, Vm, colour = factor(G))) +
      geom_hline(yintercept = -50, linetype = "dashed", colour = "grey40") +
      annotate("text", x = 0.8, y = -47, label = "Ca_v 활성화 문턱 (−50 mV)",
               size = 3.4, colour = "grey30") +
      geom_vline(xintercept = input$g, colour = "#d1495b", linewidth = 1) +
      geom_line(linewidth = 1.1) + scale_x_reverse() +
      labs(title = "g 하나가 막전위를 정한다",
           x = "잔존 K_ATP 전도도 g (오른쪽 = 정상)", y = "안정막전위 (mV)",
           colour = "혈당 (mg/dL)") + THEME
  })

  ## ---------------- ② PK ----------------
  output$pk <- renderPlot({
    d <- sim() %>%
      select(time, `디아족사이드 (µg/mL)` = DZX, `옥트레오타이드 (ng/mL)` = OCT,
             `에르소데투그 (µg/mL)` = ERS, `시롤리무스 (mg/L)` = SIR,
             `니페디핀 (mg/L)` = NIF, `글루카곤 (pg/mL)` = GCG) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.9, colour = "#2f6f9f") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "시간 (h)", y = NULL, title = "약동학") + THEME
  })

  ## ---------------- ③ PD ----------------
  output$pd <- renderPlot({
    d <- sim() %>%
      select(time, `혈당 (mg/dL)` = GLU, `CGM 간질액 (mg/dL)` = GLUi,
             `인슐린 (µU/mL)` = INS, `수용체 인슐린 X (µU/mL)` = insact,
             `C-펩타이드 (ng/mL)` = CPEP,
             `베타세포 막전위 (mV)` = Vm_ab) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.9, colour = "#c0504d") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "시간 (h)", y = NULL,
           title = "약력학 — 혈당·인슐린·막전위") + THEME
  })

  ## ---------------- ④ endpoints ----------------
  output$endpoints <- renderDT({
    s <- tail_stats(sim()); d <- sim()
    datatable(tibble(
      엔드포인트 = c("평균 혈당", "최저 혈당", "70 mg/dL 미만 시간 비율",
                     "54 mg/dL 미만 시간 비율", "저혈당 면적 (mg/dL·h)",
                     "뇌 연료 부족 누적시간 (h)", "신경발달 결손 지수",
                     "정맥 포도당 요구량", "총 포도당 공급"),
      값 = c(sprintf("%.1f mg/dL", s$glucose_mean),
             sprintf("%.1f mg/dL", s$glucose_min),
             sprintf("%.1f %%", s$pct_below70), sprintf("%.1f %%", s$pct_below54),
             sprintf("%.0f", max(d$AUCHYPO)), sprintf("%.2f", max(d$TFUEL)),
             sprintf("%.2f", max(d$DEV)),
             sprintf("%.2f mg/kg/min", s$GIR),
             sprintf("%.2f mg/kg/min", s$total_glc))
    ), options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  output$gir <- renderPlot({
    d <- sim()
    ggplot(d, aes(time)) +
      geom_ribbon(aes(ymin = 0, ymax = GIRout), fill = "#93cdf7", alpha = 0.6) +
      geom_line(aes(y = GLU/10), colour = "#c0504d", linewidth = 1) +
      scale_y_continuous(name = "정맥 포도당 (mg/kg/min)",
        sec.axis = sec_axis(~.*10, name = "혈당 (mg/dL)")) +
      labs(x = "시간 (h)", title = "포도당 요구량과 혈당") + THEME
  })

  ## ---------------- ⑤ compare ----------------
  cmp <- reactive({
    arms <- list("무치료" = list(), "디아족사이드 15" = list(dzx = 15),
                 "옥트레오타이드 30" = list(oct = 30),
                 "글루카곤 15" = list(gcg = 15),
                 "에르소데투그 9" = list(ers = 9),
                 "니페디핀 0.5" = list(nif = 0.5),
                 "옥트레오타이드+글루카곤" = list(oct = 30, gcg = 15))
    lapply(names(arms), function(n) {
      a <- arms[[n]]
      d <- do.call(simulate, c(list(gt = input$gt, g_override = input$g,
                                    remnant = input$remnant, fed = input$fed,
                                    closed_loop = TRUE, target = input$target,
                                    end = 72), a))
      tail_stats(d) %>% mutate(arm = n)
    }) %>% bind_rows()
  })

  output$compare <- renderPlot({
    d <- cmp() %>% mutate(arm = factor(arm, levels = arm))
    ggplot(d, aes(arm, GIR, fill = GIR)) +
      geom_col(width = 0.65) + coord_flip() +
      scale_fill_gradient(low = "#4c9f70", high = "#d1495b", guide = "none") +
      labs(x = NULL, y = "정맥 포도당 요구량 (mg/kg/min)",
           title = sprintf("g = %.2f 에서의 치료 반응", input$g)) + THEME
  })

  output$compare_tbl <- renderDT({
    d <- cmp()
    base <- d$GIR[1]
    datatable(d %>% transmute(치료 = arm,
        `정맥 포도당` = round(GIR, 2), `감소율 %` = round(100*(base - GIR)/max(base, 1e-9), 1),
        `평균 혈당` = round(glucose_mean, 1), `인슐린` = round(insulin, 1),
        `BOHB` = round(BOHB, 3)),
      options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  ## ---------------- ⑥ biomarkers ----------------
  output$biomarker <- renderPlot({
    d <- sim() %>%
      select(time, `BOHB (mmol/L)` = BOHB, `유리지방산 (mmol/L)` = FFA,
             `간 글리코겐 (mg/kg)` = GLY, `암모니아 (µmol/L)` = NH3,
             `젖산 (mmol/L)` = LAC, `글루카곤 (pg/mL)` = GCG) %>%
      pivot_longer(-time)
    ggplot(d, aes(time, value)) + geom_line(linewidth = 0.9, colour = "#8a7355") +
      facet_wrap(~name, scales = "free_y", ncol = 2) +
      labs(x = "시간 (h)", y = NULL, title = "바이오마커") + THEME
  })

  ## ---------------- ⑦ cerebral fuel ----------------
  output$fuel <- renderPlot({
    p <- as.list(param(mod)); d <- sim()
    dd <- d %>% mutate(
      glc  = p$UbrMax*GLU/(p$Kbr + GLU)*(1 + p$acbf*pmax(0, p$Gcbf - GLU)/p$Gcbf),
      ket  = p$VketMax*BOHB/(p$Kket + BOHB),
      lact = p$VlacMax*LAC/(p$Klac + LAC)) %>%
      select(time, `포도당` = glc, `케톤` = ket, `젖산` = lact) %>%
      pivot_longer(-time)
    ggplot(dd, aes(time, value, fill = name)) +
      geom_area() +
      geom_hline(yintercept = p$CMRreq, linetype = "dashed", linewidth = 1) +
      annotate("text", x = max(dd$time)*0.1, y = p$CMRreq*1.06,
               label = "뇌 에너지 요구량", size = 3.6) +
      scale_fill_manual(values = c(`포도당` = "#f39c9c", `케톤` = "#cbb99a",
                                   `젖산` = "#d9ecf7")) +
      labs(x = "시간 (h)", y = "뇌 연료 공급 (mg/kg/min, 포도당 등가)",
           fill = NULL, title = "뇌는 합(sum)으로 먹는다") + THEME
  })

  output$fuelthresh <- renderPlot({
    p <- as.list(param(mod))
    f <- function(G, B) {
      p$UbrMax*G/(p$Kbr + G)*(1 + p$acbf*max(0, p$Gcbf - G)/p$Gcbf) +
        p$VketMax*B/(p$Kket + B) + p$VlacMax*1.3/(p$Klac + 1.3)
    }
    d <- tibble(BOHB = seq(0, 5, 0.05)) %>% rowwise() %>%
      mutate(thr = {
        lo <- 1; hi <- 200
        for (i in 1:60) { m <- (lo + hi)/2; if (f(m, BOHB) < p$CMRreq) lo <- m else hi <- m }
        (lo + hi)/2 }) %>% ungroup()
    ggplot(d, aes(BOHB, thr)) +
      geom_line(linewidth = 1.2, colour = "#b5423f") +
      geom_point(data = tibble(BOHB = 0, thr = 43.7), size = 3, colour = "#d1495b") +
      annotate("text", x = 0.9, y = 46,
               label = "케톤 0 (CHI): 43.7 mg/dL", size = 3.8, colour = "#d1495b") +
      annotate("text", x = 2.6, y = 27,
               label = "BOHB 2 mM (케톤성): 25.3 mg/dL", size = 3.8, colour = "#8a7355") +
      labs(x = "혈장 BOHB (mmol/L)",
           y = "뇌 연료가 고갈되는 혈당 (mg/dL)",
           title = "고인슐린증은 포도당을 빼앗기 전에 안전역을 삭제한다") + THEME
  })

  ## ---------------- ⑧ g sweep ----------------
  gsw <- reactive({
    gs <- c(0.6, 0.4, 0.3, 0.2, 0.1, 0.05, 0.02)
    lapply(gs, function(g) {
      b <- simulate(input$gt, g_override = g, closed_loop = TRUE,
                    target = input$target, end = 72)
      dz <- simulate(input$gt, g_override = g, dzx = 15, closed_loop = TRUE,
                     target = input$target, end = 72)
      oc <- simulate(input$gt, g_override = g, oct = 30, closed_loop = TRUE,
                     target = input$target, end = 72)
      b0 <- tail_stats(b)$GIR
      tibble(g = g, GIR0 = b0,
             `디아족사이드` = if (b0 > 0.05) 100*(b0 - tail_stats(dz)$GIR)/b0 else NA,
             `옥트레오타이드` = if (b0 > 0.05) 100*(b0 - tail_stats(oc)$GIR)/b0 else NA)
    }) %>% bind_rows()
  })

  output$gsweep <- renderPlot({
    d <- gsw() %>% pivot_longer(c(`디아족사이드`, `옥트레오타이드`))
    ggplot(d, aes(g, value, colour = name)) +
      geom_line(linewidth = 1.3) + geom_point(size = 2.6) + scale_x_reverse() +
      scale_colour_manual(values = c(`디아족사이드` = "#2f6f9f",
                                     `옥트레오타이드` = "#c0504d")) +
      labs(x = "잔존 K_ATP 전도도 g (오른쪽 = 정상에 가까움)",
           y = "정맥 포도당 요구량 감소 (%)", colour = NULL,
           title = "하나의 방정식, 두 개의 결론",
           subtitle = "디아족사이드는 g에 곱해지고(g→0에서 소멸) 옥트레오타이드는 분배기에 더해진다") +
      THEME
  })

  output$gsweep_tbl <- renderDT({
    datatable(gsw() %>% transmute(g, `무치료 포도당` = round(GIR0, 2),
        `디아족사이드 감소%` = round(`디아족사이드`, 1),
        `옥트레오타이드 감소%` = round(`옥트레오타이드`, 1)),
      options = list(dom = "t"), rownames = FALSE)
  })

  ## ---------------- ⑨ surgery ----------------
  surg <- reactive({
    lapply(c(1.0, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02), function(b) {
      d <- simulate(input$gt, g_override = input$g, remnant = b,
                    closed_loop = TRUE, target = input$target, end = 48)
      tail_stats(d) %>% mutate(remnant = b)
    }) %>% bind_rows()
  })

  output$surgery <- renderPlot({
    d <- surg()
    ggplot(d, aes(remnant)) +
      geom_hline(yintercept = c(70, 126), linetype = "dashed", colour = "grey45") +
      geom_line(aes(y = glucose_mean), linewidth = 1.3, colour = "#c0504d") +
      geom_point(aes(y = glucose_mean), size = 2.6, colour = "#c0504d") +
      geom_col(aes(y = GIR*10), fill = "#93cdf7", alpha = 0.55, width = 0.02) +
      scale_x_log10() +
      scale_y_continuous(name = "평균 혈당 (mg/dL)",
        sec.axis = sec_axis(~./10, name = "정맥 포도당 (mg/kg/min)")) +
      labs(x = "베타세포 잔존량 (log)",
           title = "절제는 B를 줄이고 g는 남긴다 — 그래서 창(window)이 좁다") + THEME
  })

  output$surgery_tbl <- renderDT({
    datatable(surg() %>% transmute(`잔존량` = remnant,
        `정맥 포도당` = round(GIR, 2), `평균 혈당` = round(glucose_mean, 1),
        `인슐린` = round(insulin, 1),
        해석 = ifelse(GIR > 0.3, "여전히 포도당 의존",
                 ifelse(glucose_mean < 126, "정상 혈당", "당뇨"))),
      options = list(dom = "t"), rownames = FALSE)
  })
}

shinyApp(ui, server)
