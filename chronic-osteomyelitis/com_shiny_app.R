## ============================================================================
##  com_shiny_app.R
##  만성 골수염 · 임플란트 연관 골감염 QSP 모델 — 인터랙티브 대시보드
##
##  실행:
##      shiny::runApp("com_shiny_app.R")
##  또는
##      Rscript -e 'shiny::runApp("com_shiny_app.R", port = 8080)'
##
##  이 앱은 하나의 질문에 답하도록 설계되었다:
##
##      "이 환자에게 이 요법은 멸균 가능한가, 아니면 억제만 가능한가?"
##
##  그 답은 슬라이더를 움직여 얻는 그림이 아니라 화면 왼쪽 위에 계속 떠 있는
##  하나의 숫자 RSTER 로 주어진다. RSTER <= 1 이면 기간을 아무리 늘려도
##  멸균은 오지 않는다 — 탭 ⑩ 이 그것을 실험으로 보여준다.
##
##  탭 구성 (10개)
##    ①  환자 · 개입 요약        — 무엇을 바꾸면 무엇이 바뀌는지 한 화면
##    ②  두 개의 벌점            — 이 모델의 논지 (공간 벌점 vs 표현형 벌점)
##    ③  표적부위 약동학          — 혈장 → 골간질 → 격리골 → 세포내
##    ④  세균 아집단             — 5개 풀의 log10 궤적
##    ⑤  살균 동역학 · RSTER      — 약 · 면역 · 성장의 경쟁
##    ⑥  뼈 · 관류 (악순환)      — 병이 스스로 만드는 약동학 장벽
##    ⑦  임상 엔드포인트 · 바이오마커 — CRP vs ESR, 통증, 저수지, 재발확률
##    ⑧  독성 · 안전성           — 기간이 필요한 효능 vs 기간이 만드는 독성
##    ⑨  시나리오 비교           — 11개 표준 시나리오 나란히
##    ⑩  기간 스윕               — 기간이 언제 도움이 되고 언제 무의미한가
##
##  교육·연구용. 임상 의사결정에 사용하지 말 것.
## ============================================================================

library(shiny)
source("com_mrgsolve_model.R")

PAL <- c(plank = "#2563eb", bfa = "#f97316", bfp = "#b91c1c",
         ic   = "#7c3aed", res = "#0f766e", tot = "#111827",
         drugA = "#2563eb", drugB = "#16a34a", loc = "#0891b2",
         bone = "#78716c", seq = "#57534e", perf = "#db2777",
         crp  = "#ca8a04", esr = "#a16207", pain = "#9333ea",
         good = "#16a34a", bad = "#dc2626", grey = "#94a3b8")

panelplot <- function(...) {
  par(mar = c(4.0, 4.4, 2.4, 1.0), mgp = c(2.5, 0.7, 0), las = 1,
      cex.axis = 0.85, cex.lab = 0.95, cex.main = 1.0, bty = "l")
}
vlines <- function(tsurg, tstop) {
  if (is.finite(tsurg) && tsurg < 1e5)
    abline(v = tsurg, col = PAL[["good"]], lty = 2, lwd = 1.4)
  if (is.finite(tstop) && tstop < 1e5)
    abline(v = tstop, col = PAL[["grey"]], lty = 3, lwd = 1.4)
}

## ============================================================================
##  UI
## ============================================================================

ui <- fluidPage(
  titlePanel("만성 골수염 · 임플란트 연관 골감염 — QSP 시뮬레이터"),
  tags$p(style = "color:#475569;margin-top:-8px;",
    tags$b("논지: "),
    "치료 실패는 서로 독립인 두 벌점이 곱해져 생긴다 — ",
    tags$b("공간 벌점"), " (관류 없는 격리골로의 확산: eta = tanh(phi)/phi, 용량으로 밀 수 있다) 과 ",
    tags$b("표현형 벌점"), " (성장속도 의존 살균: mu→0 이면 kill→0, 어떤 용량으로도 밀 수 없다). ",
    "그래서 ", tags$b("소파술은 초기조건"), "을, ", tags$b("리팜핀은 살균천장"), "을, ",
    tags$b("기간은 도달시간"), "만을 바꾼다."),
  hr(),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("① 병소 · 숙주"),
      sliderInput("seq0", "초기 격리골 부피 SEQ (cm3)", 0.2, 20, 2.0, 0.2),
      sliderInput("aimp0", "임플란트 표면적 (cm2) — 0 = 이물질 없음", 0, 80, 30, 5),
      sliderInput("perf0", "초기 관류 분율 PERF", 0.05, 1.0, 0.37, 0.01),
      sliderInput("padsev", "PAD 중증도 (당뇨발)", 0, 0.9, 0, 0.05),
      sliderInput("hostf", "숙주 면역 능력 (1 = 정상)", 0.3, 1.2, 1.0, 0.05),
      sliderInput("egfr0", "eGFR (mL/min/1.73m2)", 20, 120, 90, 5),

      h4("② 항생제"),
      selectInput("drugA", "백본 (성장의존)", names(drug_A_library), "vancomycin"),
      numericInput("doseA", "백본 1회 용량 (mg)", 1250, 0, 6000, 50),
      selectInput("iiA", "백본 투여 간격", c("q8h" = 1/3, "q12h" = 1/2, "q24h" = 1), 1/2),
      checkboxInput("ivA", "정주 투여 (해제하면 경구)", TRUE),
      selectInput("drugB", "파트너 (성장비의존)", names(drug_B_library), "rifampicin"),
      numericInput("doseB", "파트너 1회 용량 (mg)", 450, 0, 1200, 50),
      selectInput("iiB", "파트너 투여 간격", c("q12h" = 1/2, "q24h" = 1), 1/2),
      sliderInput("weeks", "투여 기간 (주)", 1, 52, 6, 1),

      h4("③ 수술 · 국소"),
      checkboxInput("surg", "소파술 시행", TRUE),
      sliderInput("tsurg", "소파술 시점 (일)", 1, 180, 2, 1),
      sliderInput("deblog", "소파에 의한 log10 세균 감소 (P0)", 0, 7, 5.5, 0.5),
      sliderInput("seqf", "격리골 절제 비율", 0, 0.99, 0.95, 0.01),
      checkboxInput("imprem", "임플란트 제거", TRUE),
      sliderInput("flap", "피판/근육 이식에 의한 관류 회복", 0, 0.35, 0.20, 0.05),
      checkboxInput("revasc", "혈관재생술 (PAD 교정)", FALSE),
      numericInput("aloc", "국소 항생제 depot (mg, 0 = 없음)", 0, 0, 6000, 250),

      h4("④ 시뮬레이션"),
      sliderInput("tend", "관찰 기간 (일)", 100, 900, 400, 50),
      actionButton("go", "다시 계산", class = "btn-primary")
    ),

    mainPanel(
      width = 9,
      htmlOutput("verdict"),
      tabsetPanel(
        id = "tabs", type = "tabs",
        tabPanel("① 환자·개입 요약",
                 br(), tableOutput("summ"),
                 plotOutput("ovw", height = "430px"),
                 htmlOutput("interp")),
        tabPanel("② 두 개의 벌점",
                 br(), htmlOutput("penexp"),
                 plotOutput("pen", height = "560px")),
        tabPanel("③ 표적부위 약동학", br(), plotOutput("pk", height = "620px"),
                 htmlOutput("pkexp")),
        tabPanel("④ 세균 아집단", br(), plotOutput("bact", height = "560px"),
                 htmlOutput("bactexp")),
        tabPanel("⑤ 살균 동역학 · RSTER", br(), plotOutput("kill", height = "600px"),
                 htmlOutput("killexp")),
        tabPanel("⑥ 뼈 · 관류 (악순환)", br(), plotOutput("bone", height = "620px"),
                 htmlOutput("boneexp")),
        tabPanel("⑦ 임상 엔드포인트", br(), plotOutput("clin", height = "600px"),
                 htmlOutput("clinexp")),
        tabPanel("⑧ 독성", br(), plotOutput("tox", height = "560px"),
                 htmlOutput("toxexp")),
        tabPanel("⑨ 시나리오 비교", br(),
                 actionButton("runscn", "표준 시나리오 11개 실행 (수십 초)",
                              class = "btn-warning"),
                 br(), br(), tableOutput("scntab"), plotOutput("scnplot", height = "440px")),
        tabPanel("⑩ 기간 스윕", br(),
                 actionButton("runsweep", "기간 스윕 실행 (2·4·6·12·26주 x 리팜핀 유무)",
                              class = "btn-warning"),
                 br(), br(), tableOutput("swtab"), plotOutput("swplot", height = "420px"),
                 htmlOutput("swexp"))
      )
    )
  )
)

## ============================================================================
##  SERVER
## ============================================================================

server <- function(input, output, session) {

  built <- eventReactive(input$go, ignoreNULL = FALSE, {
    days <- input$weeks * 7
    m <- set_drugs(mod, input$drugA, input$drugB)
    tsurg <- if (input$surg) input$tsurg else 1e6

    m <- param(m,
      HOSTF   = input$hostf,
      PADSEV  = input$padsev,
      EGFR0   = input$egfr0,
      TSURG   = tsurg,
      SURGDUR = 0.25,
      DEBLOG  = if (input$surg) input$deblog else 0,
      DEBSEQF = if (input$surg) input$seqf   else 0,
      DEBPUSF = if (input$surg) 0.95 else 0,
      DEBICF  = if (input$surg) min(0.95, input$deblog / 6) else 0,
      IMPREM  = if (input$surg && input$imprem) 1 else 0,
      FLAPB   = if (input$surg) input$flap else 0,
      TREVASC = if (input$revasc) max(1, tsurg + 1) else 1e6,
      REVASCB = if (input$revasc) 0.75 else 0,
      TABX0   = 0,
      TABX1   = days
    )
    m <- update(m, init = list(SEQ = input$seq0, AIMP = input$aimp0,
                               PERF = input$perf0, EGFRC = input$egfr0))

    evs <- list()
    if (input$doseA > 0) {
      iiA <- as.numeric(input$iiA)
      evs[[length(evs) + 1]] <- if (input$ivA)
        ev_iv(input$doseA, iiA, days, "CENA", inf_h = 1.5)
      else ev_po(input$doseA, iiA, days, "GUTA")
    }
    if (input$doseB > 0 && input$drugB != "none")
      evs[[length(evs) + 1]] <- ev_po(input$doseB, as.numeric(input$iiB), days, "GUTB")
    ## 국소 depot 은 소파술 때 삽입하지만, 소파술을 하지 않는 시나리오에서는
    ## 0일차에 넣는다. tsurg 는 소파술이 없을 때 1e6 이므로 그 값을 그대로 쓰면
    ## solver 가 100만일까지 적분하려 들어 멈춘다.
    if (input$aloc > 0) {
      tloc <- if (input$surg) min(input$tsurg, input$tend) else 0
      evs[[length(evs) + 1]] <- ev_bolus(input$aloc, "ALOC", time = tloc)
    }

    e <- if (length(evs)) do.call(ev_bind, evs) else ev_none()
    out <- as.data.frame(run_scn(m, e, end = input$tend, delta = 0.25))
    list(sim = out, tsurg = tsurg, tstop = days, mod = m)
  })

  idx_at <- function(d, t) which.min(abs(d$time - t))

  ## ---- always-visible verdict -------------------------------------------
  output$verdict <- renderUI({
    b <- built(); d <- b$sim
    i7 <- idx_at(d, min(7, b$tstop))
    itx <- idx_at(d, b$tstop)
    rs <- d$RSTER[i7]; res <- d$RESV[itx]; pr <- d$PRLPS[itx]
    col <- if (rs > 1) PAL[["good"]] else PAL[["bad"]]
    msg <- if (rs > 1)
      "멸균 가능 (RSTER &gt; 1) — 남은 문제는 기간뿐이다."
    else
      "멸균 불가 (RSTER &le; 1) — 이 요법은 억제요법이다. 기간을 늘려도 완치는 오지 않는다."
    HTML(sprintf(
      "<div style='border-left:6px solid %s;background:#f8fafc;padding:10px 14px;margin-bottom:10px;'>
       <span style='font-size:1.25em;font-weight:700;color:%s'>RSTER(7일) = %.3f</span>
       &nbsp;&nbsp;<span style='color:%s;font-weight:600'>%s</span><br>
       <span style='color:#475569'>치료 종료 시 저수지 = <b>%.3g CFU</b> ·
       재발 확률 = <b>%.1f%%</b> · 격리골 침투효율 eta(백본) = <b>%.3f</b> ·
       확산거리 L = <b>%.3f cm</b> · 지속체 살균 = <b>%.3f/d</b> (약) + <b>%.3f/d</b> (면역)</span></div>",
      col, col, rs, col, msg, res, 100 * pr,
      d$ETAA[i7], d$LEFFO[i7], d$KILLPS[i7], d$KIMMPS[i7]))
  })

  ## ---- ① summary ---------------------------------------------------------
  output$summ <- renderTable({
    b <- built(); d <- b$sim
    i0 <- 1; i7 <- idx_at(d, min(7, b$tstop)); itx <- idx_at(d, b$tstop); ie <- nrow(d)
    pick <- function(i) c(d$LBTOT[i], d$LBFP[i], d$LBIC[i], d$LBRES[i], d$SEQ[i],
                          d$PERF[i], d$BMLOSS[i], d$CRP[i], d$ESR[i], d$PAIN[i])
    tb <- data.frame(
      metric = c("\uc9c0\ud45c \u2014 \ucd1d \uc138\uae30\uc7ac", "", "", "", "", "", "", "", "", ""),
      d0 = pick(i0), d7 = pick(i7), dtx = pick(itx), dend = pick(ie),
      stringsAsFactors = FALSE)
    tb$metric <- c("\ucd1d \uc138\uad70 log10", "\uc9c0\uc18d\uccb4 log10",
                   "\uc138\ud5b5\ub0b4 log10", "rpoB \ub0b4\uc131 log10",
                   "\uac70\ub9ac\uacf0 SEQ (cm3)", "\uad00\ub958 PERF",
                   "\ud53c\uc9c8\uacf0 \uc18c\uc2e4 (%)",
                   "CRP (mg/L)", "ESR (mm/h)", "\ud1b5\uc99d NRS")
    names(tb) <- c("\uc9c0\ud45c", "0d", "7d",
                   "\uce58\ub8cc\uc885\ub8cc", "\uad00\uc0b0\uc885\ub8cc")
    tb
  }, digits = 3, striped = TRUE, hover = TRUE, width = "100%")

  output$ovw <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(1, 3)); panelplot()
    plot(d$time, d$LBTOT, type = "l", lwd = 2.4, col = PAL[["tot"]], ylim = c(-6.5, 10.5),
         xlab = "일", ylab = "log10 CFU", main = "총 세균 부하")
    lines(d$time, d$LBFP, lwd = 2, col = PAL[["bfp"]])
    abline(h = 0, col = PAL[["grey"]], lty = 3)
    vlines(b$tsurg, b$tstop)
    legend("topright", c("총량", "지속체", "1 CFU"), col = c(PAL[["tot"]], PAL[["bfp"]], PAL[["grey"]]),
           lty = c(1, 1, 3), lwd = 2, bty = "n", cex = 0.85)

    plot(d$time, d$RSTER, type = "l", lwd = 2.4, col = PAL[["good"]],
         xlab = "일", ylab = "RSTER", main = "멸균비 RSTER")
    abline(h = 1, col = PAL[["bad"]], lty = 2, lwd = 1.6)
    vlines(b$tsurg, b$tstop)

    plot(d$time, d$PERF, type = "l", lwd = 2.4, col = PAL[["perf"]], ylim = c(0, 1),
         xlab = "일", ylab = "PERF / eta", main = "관류와 침투효율")
    lines(d$time, d$ETAA, lwd = 2, col = PAL[["drugA"]])
    lines(d$time, d$ETAB, lwd = 2, col = PAL[["drugB"]], lty = 2)
    vlines(b$tsurg, b$tstop)
    legend("bottomright", c("PERF", "eta 백본", "eta 파트너"),
           col = c(PAL[["perf"]], PAL[["drugA"]], PAL[["drugB"]]),
           lty = c(1, 1, 2), lwd = 2, bty = "n", cex = 0.85)
  })

  output$interp <- renderUI(HTML(
    "<p style='color:#475569'>세로 점선: <span style='color:#16a34a'>초록 = 소파술</span>,
     <span style='color:#94a3b8'>회색 = 항생제 중단</span>.
     중단 이후에도 총량이 계속 내려가면 저수지가 실제로 비었다는 뜻이고,
     중단 후 다시 올라가면 살아남은 저수지에서 재증식한 것이다.
     오른쪽 그림의 eta 는 <b>병소 기하가 바뀌면 같은 약의 침투가 바뀐다</b>는 사실을
     보여준다 — 소파술의 효과 중 일부는 세균을 줄인 것이 아니라 확산거리를 줄인 것이다.</p>"))

  ## ---- ② the two penalties ----------------------------------------------
  output$penexp <- renderUI(HTML(
    "<p style='color:#475569'><b>왼쪽(공간 벌점)</b>은 농도 문제다 — 곡선을 따라 오른쪽으로 갈수록
     침투가 나빠지지만, 리<b>ㅁ</b>에서의 농도를 100배 올리면 심부 농도도 100배 오른다.
     <b>오른쪽(표현형 벌점)</b>은 농도 문제가 아니다 — 농도를 포화시켜도(h=1로 두어도)
     지속체 막대의 높이는 오직 Emax_gi 로만 결정된다. 리팜핀만 그 막대를 갖고 있고,
     그것이 이 병에서 리팜핀이 특별한 유일한 이유다.</p>"))

  output$pen <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(2, 2)); panelplot()

    ## (a) Thiele penalty curves
    Ls <- seq(0.005, 0.6, length.out = 200)
    pA <- drug_A_library[[input$drugA]]
    eta <- function(L, k, D) { p <- L * sqrt(k / D); ifelse(p < 1e-6, 1, tanh(p) / p) }
    plot(Ls, eta(Ls, pA$AKLS, pA$ADEFF), type = "l", lwd = 2.6, col = PAL[["drugA"]],
         ylim = c(0, 1), xlab = "확산거리 L (cm)", ylab = "침투효율 eta",
         main = "(a) 공간 벌점: eta = tanh(phi)/phi")
    pB <- drug_B_library[[input$drugB]]
    lines(Ls, eta(Ls, pB$BKLS, pB$BDEFF), lwd = 2.6, col = PAL[["drugB"]], lty = 2)
    Lnow <- d$LEFFO[idx_at(d, min(7, b$tstop))]
    abline(v = Lnow, col = PAL[["bad"]], lty = 3, lwd = 1.6)
    text(Lnow, 0.95, sprintf(" 현재 L = %.3f", Lnow), col = PAL[["bad"]], cex = 0.85, adj = 0)
    legend("topright", c(input$drugA, input$drugB), col = c(PAL[["drugA"]], PAL[["drugB"]]),
           lty = c(1, 2), lwd = 2.4, bty = "n", cex = 0.85)

    ## (b) Emax collapse across pools, at saturating concentration
    pools <- c("부유\n(mu/mumax=0.35)", "필름 변연\n(0.12)", "지속체\n(0)")
    kA <- c(pA$AEMGD * 0.35 + pA$AEMGI, pA$AEMGD * 0.12 + pA$AEMGI, pA$AEMGI)
    kB <- if (input$drugB == "none") c(0, 0, 0) else
      c(pB$BEMGD * 0.35 + pB$BEMGI, pB$BEMGD * 0.12 + pB$BEMGI, pB$BEMGI)
    mm <- rbind(kA, kB)
    rownames(mm) <- c("A", "B")
    barplot(mm, beside = TRUE, names.arg = pools, col = c(PAL[["drugA"]], PAL[["drugB"]]),
            ylab = "포화농도에서의 살균속도 (1/d)",
            main = "(b) 표현형 벌점: Emax 붕괴 (농도 포화 가정)")
    legend("topright", c("백본", "파트너"), fill = c(PAL[["drugA"]], PAL[["drugB"]]),
           bty = "n", cex = 0.85)
    abline(h = 6 * 0.12, col = PAL[["bad"]], lty = 2)
    text(1, 6 * 0.12, " 바이오필름 성장속도", col = PAL[["bad"]], adj = 0, cex = 0.8)

    ## (c) the two penalties multiplied, over time
    plot(d$time, d$ETAA, type = "l", lwd = 2.4, col = PAL[["drugA"]], ylim = c(0, 1),
         xlab = "일", ylab = "무차원", main = "(c) 두 벌점의 시간 변화")
    lines(d$time, d$MUBF / 6, lwd = 2.4, col = PAL[["bfa"]])
    lines(d$time, d$ETAA * d$MUBF / 6, lwd = 3, col = PAL[["bad"]])
    vlines(b$tsurg, b$tstop)
    legend("topright", c("eta (공간)", "mu/mumax (표현형)", "곱 (실효 살균 여력)"),
           col = c(PAL[["drugA"]], PAL[["bfa"]], PAL[["bad"]]), lwd = 2.4, bty = "n", cex = 0.8)

    ## (d) T_ster map
    kwake <- as.numeric(param(b$mod)$KWAKE)
    kills <- seq(0, 1.5, length.out = 200)
    plot(kills, log(1e8) / (kills + kwake), type = "l", lwd = 2.4, col = PAL[["bfp"]],
         log = "y", ylim = c(1, 400), xlab = "지속체 살균속도 KILLPS (1/d)",
         ylab = "T_ster (일, 로그축)", main = "(d) T_ster = ln(P0)/(KILLPS + KWAKE)")
    lines(kills, log(1e5) / (kills + kwake), lwd = 2.4, col = PAL[["bfa"]])
    lines(kills, log(1e2) / (kills + kwake), lwd = 2.4, col = PAL[["good"]])
    abline(h = c(42, 84), col = PAL[["grey"]], lty = 3)
    text(1.2, 45, "6주", col = PAL[["grey"]], cex = 0.8)
    text(1.2, 90, "12주", col = PAL[["grey"]], cex = 0.8)
    kp <- d$KILLPS[idx_at(d, min(7, b$tstop))] + d$KIMMPS[idx_at(d, min(7, b$tstop))]
    abline(v = kp, col = "#111827", lty = 2, lwd = 1.6)
    legend("topright", c("P0 = 1e8 (소파 없음)", "P0 = 1e5", "P0 = 1e2 (근치적 소파)"),
           col = c(PAL[["bfp"]], PAL[["bfa"]], PAL[["good"]]), lwd = 2.4, bty = "n", cex = 0.8)
  })

  ## ---- ③ target-site PK --------------------------------------------------
  output$pk <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(2, 2)); panelplot()
    w <- d$time <= min(max(d$time), b$tstop + 14)

    plot(d$time[w], d$CAOBS[w], type = "l", lwd = 1.6, col = PAL[["drugA"]],
         xlab = "일", ylab = "mg/L", main = "(a) 혈장 농도")
    lines(d$time[w], d$CBOBS[w], lwd = 1.6, col = PAL[["drugB"]])
    legend("topright", c("백본", "파트너"), col = c(PAL[["drugA"]], PAL[["drugB"]]),
           lwd = 2, bty = "n", cex = 0.85)

    ymax <- max(1e-3, d$CBA[w], d$CSA[w], d$CIA[w], na.rm = TRUE)
    plot(d$time[w], pmax(d$CBA[w], 1e-4), type = "l", lwd = 2, col = PAL[["drugA"]],
         log = "y", ylim = c(max(1e-4, ymax / 1e4), ymax * 1.3),
         xlab = "일", ylab = "mg/L (로그축)", main = "(b) 백본: 구획별 농도")
    lines(d$time[w], pmax(d$CSA[w], 1e-4), lwd = 2, col = PAL[["seq"]])
    lines(d$time[w], pmax(d$CIA[w], 1e-4), lwd = 2, col = PAL[["ic"]])
    if (max(d$CLOC) > 0) lines(d$time[w], pmax(d$CLOC[w], 1e-4), lwd = 2, col = PAL[["loc"]], lty = 2)
    legend("topright", c("골간질", "격리골/필름", "세포내", "국소 창상액"),
           col = c(PAL[["drugA"]], PAL[["seq"]], PAL[["ic"]], PAL[["loc"]]),
           lty = c(1, 1, 1, 2), lwd = 2, bty = "n", cex = 0.8)

    ymax2 <- max(1e-4, d$CBB[w], d$CSB[w], d$CIB[w], na.rm = TRUE)
    plot(d$time[w], pmax(d$CBB[w], 1e-5), type = "l", lwd = 2, col = PAL[["drugB"]],
         log = "y", ylim = c(max(1e-5, ymax2 / 1e4), ymax2 * 1.3),
         xlab = "일", ylab = "mg/L (로그축)", main = "(c) 파트너: 구획별 농도")
    lines(d$time[w], pmax(d$CSB[w], 1e-5), lwd = 2, col = PAL[["seq"]])
    lines(d$time[w], pmax(d$CIB[w], 1e-5), lwd = 2, col = PAL[["ic"]])
    legend("topright", c("골간질", "격리골/필름", "세포내"),
           col = c(PAL[["drugB"]], PAL[["seq"]], PAL[["ic"]]), lwd = 2, bty = "n", cex = 0.8)

    plot(d$time[w], d$BSRA[w], type = "l", lwd = 2.2, col = PAL[["bone"]],
         ylim = c(0, max(0.5, max(d$BSRA[w], na.rm = TRUE))),
         xlab = "일", ylab = "비", main = "(d) 겉보기 골:혈청 비 · 유도 상태")
    lines(d$time[w], d$ETAA[w], lwd = 2, col = PAL[["drugA"]], lty = 2)
    vlines(b$tsurg, b$tstop)
    legend("topright", c("골:혈청 비 (관류 의존)", "eta (기하 의존)"),
           col = c(PAL[["bone"]], PAL[["drugA"]]), lty = c(1, 2), lwd = 2, bty = "n", cex = 0.8)
  })

  output$pkexp <- renderUI(HTML(
    "<p style='color:#475569'>(b)/(c) 의 세 선 사이 간격이 이 병의 약동학이다.
     혈장에서 골간질로 한 번, 골간질에서 격리골 심부로 또 한 번 곱셈으로 깎인다.
     국소 depot(파란 점선)은 <b>첫 번째 곱셈만</b> 우회한다 — 두 번째(eta)는 그대로 남는다.
     그래서 탭 ⑨ 의 S7 에서 창상액 농도가 MIC의 수십 배여도 바이오필름은 살아남는다.</p>"))

  ## ---- ④ bacterial pools --------------------------------------------------
  output$bact <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(1, 2)); panelplot()
    plot(d$time, d$LBTOT, type = "l", lwd = 3, col = PAL[["tot"]], ylim = c(-6.5, 10.5),
         xlab = "일", ylab = "log10 CFU", main = "(a) 아집단별 궤적")
    lines(d$time, d$LBPL,  lwd = 2, col = PAL[["plank"]])
    lines(d$time, d$LBFA,  lwd = 2, col = PAL[["bfa"]])
    lines(d$time, d$LBFP,  lwd = 2.4, col = PAL[["bfp"]])
    lines(d$time, d$LBIC,  lwd = 2, col = PAL[["ic"]])
    lines(d$time, d$LBRES, lwd = 2, col = PAL[["res"]], lty = 2)
    abline(h = 0, col = PAL[["grey"]], lty = 3)
    vlines(b$tsurg, b$tstop)
    legend("bottomleft", c("총량", "부유", "필름 변연", "지속체", "세포내", "rpoB 내성"),
           col = c(PAL[["tot"]], PAL[["plank"]], PAL[["bfa"]], PAL[["bfp"]], PAL[["ic"]], PAL[["res"]]),
           lty = c(1, 1, 1, 1, 1, 2), lwd = 2, bty = "n", cex = 0.8, ncol = 2)

    plot(d$time, d$RESV, type = "l", lwd = 2.6, col = PAL[["bfp"]], log = "y",
         ylim = c(max(1e-6, min(d$RESV[d$RESV > 0])), max(1e10, max(d$RESV))),
         xlab = "일", ylab = "저수지 (CFU, 로그축)",
         main = "(b) 재발의 원료: 저수지와 재발확률")
    abline(h = 1, col = PAL[["bad"]], lty = 2, lwd = 1.6)
    par(new = TRUE)
    plot(d$time, d$PRLPS, type = "l", lwd = 2, col = PAL[["grey"]], axes = FALSE,
         xlab = "", ylab = "", ylim = c(0, 1))
    axis(4, col = PAL[["grey"]], col.axis = PAL[["grey"]])
    mtext("재발확률", side = 4, line = 2.2, col = PAL[["grey"]], cex = 0.85, las = 0)
    vlines(b$tsurg, b$tstop)
  })

  output$bactexp <- renderUI(HTML(
    "<p style='color:#475569'>지속체(빨강)가 다른 선들보다 늦게 내려간다면 그 요법에는
     Emax_gi 가 없다는 뜻이다. 그 경우 지속체는 <b>각성(KWAKE)을 통해서만</b> 사라지므로
     필요한 기간이 수백 일 규모가 된다 — 임상적으로는 '완치되지 않는다'와 같은 말이다.
     오른쪽 그림의 빨간 수평선(1 CFU)이 실질적인 완치 경계다: 결정론적 ODE 는 1 CFU 이하를
     표현할 수 없으므로, 그 아래는 재발확률 1-exp(-p*N) 로 넘긴다.</p>"))

  ## ---- ⑤ killing ----------------------------------------------------------
  output$kill <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(2, 2)); panelplot()

    plot(d$time, d$KILLBF, type = "l", lwd = 2.4, col = PAL[["drugA"]],
         ylim = c(0, max(0.5, max(d$KILLBF, d$MUBF, d$KIMMBF, na.rm = TRUE))),
         xlab = "일", ylab = "1/d", main = "(a) 필름 변연: 약 vs 면역 vs 성장")
    lines(d$time, d$KIMMBF, lwd = 2.4, col = PAL[["ic"]])
    lines(d$time, d$MUBF, lwd = 2.4, col = PAL[["bad"]], lty = 2)
    vlines(b$tsurg, b$tstop)
    legend("topright", c("약물 살균", "면역 살균", "성장 (넘어야 하는 선)"),
           col = c(PAL[["drugA"]], PAL[["ic"]], PAL[["bad"]]),
           lty = c(1, 1, 2), lwd = 2.2, bty = "n", cex = 0.8)

    plot(d$time, d$KILLPS, type = "l", lwd = 2.4, col = PAL[["drugA"]],
         ylim = c(0, max(0.2, max(d$KILLPS + d$KIMMPS, na.rm = TRUE))),
         xlab = "일", ylab = "1/d", main = "(b) 지속체: 무엇이 실제로 죽이는가")
    lines(d$time, d$KIMMPS, lwd = 2.4, col = PAL[["ic"]])
    abline(h = as.numeric(param(b$mod)$KWAKE), col = PAL[["grey"]], lty = 3, lwd = 1.6)
    vlines(b$tsurg, b$tstop)
    legend("topright", c("약물", "면역(감염세포 제거 포함)", "각성속도 KWAKE"),
           col = c(PAL[["drugA"]], PAL[["ic"]], PAL[["grey"]]),
           lty = c(1, 1, 3), lwd = 2.2, bty = "n", cex = 0.8)

    plot(d$time, d$RSTER, type = "l", lwd = 2.8, col = PAL[["good"]],
         xlab = "일", ylab = "RSTER", main = "(c) 멸균비 — 1을 넘느냐가 전부다")
    abline(h = 1, col = PAL[["bad"]], lty = 2, lwd = 1.8)
    vlines(b$tsurg, b$tstop)

    ts <- pmin(d$TSTER, 600)
    plot(d$time, ts, type = "l", lwd = 2.4, col = PAL[["bfp"]],
         xlab = "일", ylab = "T_ster (일)", main = "(d) 남은 멸균 소요시간")
    abline(h = c(42, 84), col = PAL[["grey"]], lty = 3)
    vlines(b$tsurg, b$tstop)
  })

  output$killexp <- renderUI(HTML(
    "<p style='color:#475569'>(a) 에서 파란 선 + 보라 선이 빨간 점선 위로 올라가야 바이오필름이 줄어든다.
     (b) 는 이 모델의 가장 임상적인 그림이다: 지속체를 죽이는 것이 <b>약</b>인지
     <b>면역</b>인지 <b>단지 기다림(KWAKE)</b>인지를 분리해서 보여준다. 리팜핀 없는 요법에서는
     파란 선이 회색 점선 아래로 깔리고, 그때 완치 여부는 사실상 소파술이 남긴 P0 하나로 정해진다.</p>"))

  ## ---- ⑥ bone & perfusion -------------------------------------------------
  output$bone <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(2, 2)); panelplot()

    plot(d$time, d$BM, type = "l", lwd = 2.4, col = PAL[["bone"]], ylim = c(0, 105),
         xlab = "일", ylab = "상대값", main = "(a) 생존 피질골과 격리골")
    par(new = TRUE)
    plot(d$time, d$SEQ, type = "l", lwd = 2.4, col = PAL[["seq"]], axes = FALSE,
         xlab = "", ylab = "")
    axis(4); mtext("SEQ (cm3)", side = 4, line = 2.2, cex = 0.85, las = 0)
    vlines(b$tsurg, b$tstop)
    legend("right", c("생존 골(좌축)", "격리골(우축)"),
           col = c(PAL[["bone"]], PAL[["seq"]]), lwd = 2.2, bty = "n", cex = 0.8)

    plot(d$time, d$PERF, type = "l", lwd = 2.6, col = PAL[["perf"]], ylim = c(0, 1),
         xlab = "일", ylab = "PERF", main = "(b) 관류 — 약과 면역의 공통 관문")
    abline(h = as.numeric(param(b$mod)$PERFCRIT), col = PAL[["bad"]], lty = 2)
    text(max(d$time) * 0.6, as.numeric(param(b$mod)$PERFCRIT) + 0.03,
         "이 아래에서 골이 죽는다", col = PAL[["bad"]], cex = 0.8)
    vlines(b$tsurg, b$tstop)

    plot(d$time, d$OC, type = "l", lwd = 2.2, col = PAL[["bad"]],
         xlab = "일", ylab = "상대값", main = "(c) 골파괴세포 · 골형성세포 · 농양",
         ylim = c(0, max(d$OC, d$OB, d$PUS, na.rm = TRUE)))
    lines(d$time, d$OB, lwd = 2.2, col = PAL[["good"]])
    lines(d$time, d$PUS, lwd = 2.2, col = PAL[["ic"]])
    vlines(b$tsurg, b$tstop)
    legend("topright", c("골파괴세포", "골형성세포", "농양(mL)"),
           col = c(PAL[["bad"]], PAL[["good"]], PAL[["ic"]]), lwd = 2.2, bty = "n", cex = 0.8)

    plot(d$SEQ, d$ETAA, type = "l", lwd = 2.6, col = PAL[["drugA"]],
         xlab = "격리골 SEQ (cm3)", ylab = "침투효율 eta",
         main = "(d) 악순환의 단면: SEQ 가 커지면 약이 못 들어간다")
    points(d$SEQ[1], d$ETAA[1], pch = 19, col = PAL[["bad"]], cex = 1.3)
    points(d$SEQ[nrow(d)], d$ETAA[nrow(d)], pch = 17, col = PAL[["good"]], cex = 1.3)
    legend("topright", c("시작", "종료"), pch = c(19, 17),
           col = c(PAL[["bad"]], PAL[["good"]]), bty = "n", cex = 0.85)
  })

  output$boneexp <- renderUI(HTML(
    "<p style='color:#475569'>(d) 가 이 모델의 핵심 고리를 한 장에 담은 그림이다.
     감염이 뼈를 죽이고(격리골 증가) 죽은 뼈는 확산거리를 늘려 약의 침투를 떨어뜨리며,
     떨어진 침투는 감염을 지속시킨다. 소파술은 이 곡선 위에서 점을 <b>왼쪽 위로</b> 되돌리는
     유일한 개입이고, 그 점의 위치가 같은 약의 효과를 몇 배로 바꾼다.</p>"))

  ## ---- ⑦ clinical --------------------------------------------------------
  output$clin <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(2, 2)); panelplot()

    plot(d$time, d$CRP, type = "l", lwd = 2.6, col = PAL[["crp"]],
         ylim = c(0, max(d$CRP, d$ESR, na.rm = TRUE) * 1.05),
         xlab = "일", ylab = "CRP (mg/L) / ESR (mm/h)",
         main = "(a) 같은 사건, 다른 시계: CRP 19h vs ESR 10d")
    lines(d$time, d$ESR, lwd = 2.6, col = PAL[["esr"]], lty = 2)
    vlines(b$tsurg, b$tstop)
    legend("topright", c("CRP", "ESR"), col = c(PAL[["crp"]], PAL[["esr"]]),
           lty = c(1, 2), lwd = 2.4, bty = "n", cex = 0.85)

    plot(d$time, d$CRP / max(d$CRP), type = "l", lwd = 2.4, col = PAL[["crp"]],
         ylim = c(0, 1.05), xlab = "일", ylab = "정규화 (0-1)",
         main = "(b) 바이오마커는 저수지를 보지 못한다")
    lines(d$time, d$RESV / max(d$RESV), lwd = 2.6, col = PAL[["bfp"]])
    vlines(b$tsurg, b$tstop)
    legend("topright", c("CRP (정규화)", "저수지 (정규화)"),
           col = c(PAL[["crp"]], PAL[["bfp"]]), lwd = 2.4, bty = "n", cex = 0.85)

    plot(d$time, d$PAIN, type = "l", lwd = 2.6, col = PAL[["pain"]], ylim = c(0, 10),
         xlab = "일", ylab = "NRS", main = "(c) 통증")
    vlines(b$tsurg, b$tstop)

    plot(d$time, d$PRLPS, type = "l", lwd = 2.8, col = PAL[["bad"]], ylim = c(0, 1),
         xlab = "일", ylab = "확률", main = "(d) 저수지로부터의 재발확률")
    vlines(b$tsurg, b$tstop)
  })

  output$clinexp <- renderUI(HTML(
    "<p style='color:#475569'>(b) 가 임상에서 가장 자주 속는 지점이다. CRP 는 부유 세균과
     사이토카인을 따라가므로 <b>바이오필름이 그대로 남아 있어도 정상화될 수 있다</b>.
     '염증 수치가 떨어졌으니 끊어도 된다'가 위험한 이유가 이 두 곡선의 어긋남이다.
     ESR 은 반대 방향으로 틀린다 — 실제로는 나았는데도 수 주 동안 계속 높게 남는다.</p>"))

  ## ---- ⑧ toxicity --------------------------------------------------------
  output$tox <- renderPlot({
    b <- built(); d <- b$sim
    par(mfrow = c(2, 2)); panelplot()
    plot(d$time, d$EGFRC, type = "l", lwd = 2.4, col = PAL[["bad"]],
         xlab = "일", ylab = "mL/min/1.73m2", main = "(a) eGFR (반코마이신 AKI)")
    vlines(b$tsurg, b$tstop)
    plot(d$time, d$AUC24A, type = "l", lwd = 2.4, col = PAL[["drugA"]],
         xlab = "일", ylab = "mg*h/L", main = "(b) 백본 AUC24 (목표대 400-600)")
    abline(h = c(400, 600), col = PAL[["grey"]], lty = c(3, 2))
    vlines(b$tsurg, b$tstop)
    plot(d$time, d$PLT, type = "l", lwd = 2.4, col = PAL[["ic"]],
         xlab = "일", ylab = "10^9/L", main = "(c) 혈소판 (리네졸리드, 기간 의존)")
    abline(h = 100, col = PAL[["bad"]], lty = 2); abline(v = 28, col = PAL[["grey"]], lty = 3)
    vlines(b$tsurg, b$tstop)
    plot(d$time, d$ALT, type = "l", lwd = 2.4, col = PAL[["esr"]],
         xlab = "일", ylab = "U/L", main = "(d) ALT (리팜핀)")
    abline(h = 3 * 25, col = PAL[["bad"]], lty = 2)
    vlines(b$tsurg, b$tstop)
  })

  output$toxexp <- renderUI(HTML(
    "<p style='color:#475569'>효능은 기간을 <b>요구</b>하고 독성은 기간이 <b>만든다</b>.
     그래서 최적 기간은 T* = argmax [P_cure(T) - sum w_k Tox_k(T)] 이고,
     소파술은 P_cure(T) 곡선을 왼쪽으로 밀어 T* 를 줄인다 — 수술이 곧 독성 감축이다.
     반코마이신은 신기능을 깎고 그 신기능 저하가 다시 노출을 올리는 되먹임을 갖는다((a)+(b)).</p>"))

  ## ---- ⑨ scenarios -------------------------------------------------------
  scn_res <- eventReactive(input$runscn, {
    scns <- list(scn01(), scn02(), scn03(), scn04(), scn04b(), scn05(),
                 scn06(FALSE), scn06(TRUE), scn07(), scn08(FALSE), scn08(TRUE),
                 scn10(), scn11(2), scn11(120))
    lapply(scns, summarise_scn, verbose = FALSE)
  })

  output$scntab <- renderTable({
    r <- scn_res()
    tb <- do.call(rbind, lapply(r, function(x) x$summary))
    tb[, c("scenario", "RSTER_d7", "KILLBF_d7", "KILLPS_d7", "ETAA_d7",
           "logB_tx", "logPersist", "logRes_tx", "reservoir", "P_relapse",
           "logB_end", "CRP_tx", "SEQ_end", "PERF_end")]
  }, digits = 3, striped = TRUE, hover = TRUE, width = "100%")

  output$scnplot <- renderPlot({
    r <- scn_res()
    panelplot(); par(mfrow = c(1, 2))
    tb <- do.call(rbind, lapply(r, function(x) x$summary))
    cols <- ifelse(tb$P_relapse < 0.5, PAL[["good"]], PAL[["bad"]])
    barplot(tb$RSTER_d7, horiz = TRUE, col = cols, border = NA,
            names.arg = sub(" .*", "", tb$scenario), las = 1, cex.names = 0.7,
            xlab = "RSTER (7일)", main = "요법별 멸균 가능성")
    abline(v = 1, col = "#111827", lty = 2, lwd = 1.8)
    plot(tb$RSTER_d7, pmax(tb$reservoir, 1e-8), log = "y", pch = 19, col = cols, cex = 1.5,
         xlab = "RSTER (7일)", ylab = "치료종료 저수지 (CFU, 로그축)",
         main = "RSTER 이 1을 넘는지가 결과를 가른다")
    abline(v = 1, col = "#111827", lty = 2); abline(h = 1, col = PAL[["bad"]], lty = 2)
    text(tb$RSTER_d7, pmax(tb$reservoir, 1e-8), sub(" .*", "", tb$scenario),
         pos = 4, cex = 0.65, col = "#475569")
  })

  ## ---- ⑩ duration sweep --------------------------------------------------
  sw_res <- eventReactive(input$runsweep, {
    out <- list()
    for (w in c(2, 4, 6, 12, 26)) for (rif in c(TRUE, FALSE)) {
      s <- scn09(w, with_rif = rif, deblog = 3.0, impkeep = TRUE)
      x <- summarise_scn(s, verbose = FALSE)$summary
      x$weeks <- w; x$rif <- rif
      out[[length(out) + 1]] <- x
    }
    do.call(rbind, out)
  })

  output$swtab <- renderTable({
    sw <- sw_res()
    sw$regimen <- ifelse(sw$rif, "\ub808\ubcf4\ud50c\ub85d\uc0ac\uc2e0 + \ub9ac\ud31c\ud54c", "\ub808\ubcf4\ud50c\ub85d\uc0ac\uc2e0 \ub2e8\ub3c5")
    sw[, c("regimen", "weeks", "RSTER_d7", "KILLPS_d7", "logB_tx", "logPersist",
           "reservoir", "P_relapse", "logB_end")]
  }, digits = 3, striped = TRUE, hover = TRUE, width = "100%")

  output$swplot <- renderPlot({
    sw <- sw_res(); panelplot(); par(mfrow = c(1, 2))
    a <- sw[sw$rif, ]; bq <- sw[!sw$rif, ]
    plot(a$weeks, pmax(a$reservoir, 1e-8), type = "b", log = "y", pch = 19, lwd = 2.4,
         col = PAL[["good"]], ylim = c(1e-8, 1e10),
         xlab = "투여 기간 (주)", ylab = "치료종료 저수지 (CFU, 로그축)",
         main = "기간은 언제 도움이 되는가")
    lines(bq$weeks, pmax(bq$reservoir, 1e-8), type = "b", pch = 17, lwd = 2.4, col = PAL[["bad"]])
    abline(h = 1, col = "#111827", lty = 2, lwd = 1.6)
    legend("topright", c("리팜핀 병용 (RSTER > 1)", "단독 (RSTER < 1)"),
           col = c(PAL[["good"]], PAL[["bad"]]), pch = c(19, 17), lwd = 2.2, bty = "n", cex = 0.85)
    plot(a$weeks, a$P_relapse, type = "b", pch = 19, lwd = 2.4, col = PAL[["good"]],
         ylim = c(0, 1), xlab = "투여 기간 (주)", ylab = "재발확률",
         main = "재발확률 대 기간")
    lines(bq$weeks, bq$P_relapse, type = "b", pch = 17, lwd = 2.4, col = PAL[["bad"]])
  })

  output$swexp <- renderUI(HTML(
    "<p style='color:#475569'>이 탭이 이 모델 전체의 결론이다. 리팜핀을 병용한 쪽(초록)은
     <b>문턱</b>을 갖는다 — 2주는 부족하고 4주부터 저수지가 1 CFU 아래로 떨어지며 그 뒤로는
     더 늘려도 달라지지 않는다. 리팜핀이 없는 쪽(빨강)은 기간에 <b>단조롭게 반응하지 않는다</b>:
     2주에서 26주까지 저수지가 오히려 커진다. RSTER &le; 1 일 때 기간은 부족한 것이 아니라
     <b>무관</b>하다는 뜻이고, 그때 필요한 것은 더 긴 항생제가 아니라 다시 수술방이다.</p>"))
}

shinyApp(ui, server)
