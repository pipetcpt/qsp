## =====================================================================
##  Multiple System Atrophy (MSA) — QSP Shiny dashboard
##  ---------------------------------------------------------------------
##  10 tabs, driven by the 61-ODE model in msa_mrgsolve_model.R.
##
##  DESIGN NOTE — why the compute functions are top-level.
##  Every panel is produced by a plain function (`tab_*`) that takes a
##  settings list and returns data / draws a plot.  Nothing that computes
##  lives inside renderPlot().  That makes the whole dashboard testable
##  head-less: msa_shiny_selftest() below walks all 10 tabs on a null PNG
##  device and asserts that each one produces finite numbers.  Run it with
##      Rscript -e 'source("msa_shiny_app.R"); msa_shiny_selftest()'
##  which needs mrgsolve but NOT a browser and NOT shiny.
##
##  THE ORGANISING IDEA THE DASHBOARD IS BUILT AROUND.
##  Tab 1 does not show "disease severity".  It shows the SEVEN MULTIPLICANDS
##  of upright blood pressure separately, because the whole clinical logic of
##  MSA is that different diseases break different factors of the same
##  product, and that is what decides which drug works.  Every other tab is
##  a consequence of that decomposition.
##
##  Launch:  Rscript -e 'shiny::runApp("msa_shiny_app.R", port = 8080)'
##  Educational / research use only.  NOT for clinical decision-making.
## =====================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
})

## The model file defines: mod, pheno, pp, msa_history, msa_state_at,
## msa_tilt, tilt_metrics, msa_day, day_metrics, rx, rx_chronic,
## msa_history_rx.  Source it without triggering its own __main__ block.
.msa_src <- local({
  f <- "msa_mrgsolve_model.R"
  if (!file.exists(f)) f <- file.path(dirname(sys.frame(1)$ofile %||% "."), f)
  f
})
`%||%` <- function(a, b) if (is.null(a)) b else a
local({
  src <- readLines("msa_mrgsolve_model.R", warn = FALSE)
  cut <- grep('^if \\(!interactive\\(\\)', src)
  if (length(cut)) src <- src[seq_len(cut[1] - 1L)]
  eval(parse(text = paste(src, collapse = "\n")), envir = globalenv())
})

PH_LABEL <- c(MSA_P = "MSA-P (선조체흑질 우세)",
              MSA_C = "MSA-C (올리보폰토소뇌 우세)",
              PAF   = "PAF (순수 자율신경실패, 절후 병소)",
              PD    = "PD (파킨슨병, 후시냅스 온전)")

PAL <- c("#2b6cb0", "#c53030", "#2f855a", "#b7791f", "#6b46c1",
         "#0987a0", "#97266d", "#4a5568")

## ---------------------------------------------------------------------
##  settings -> simulations (cached per settings signature)
## ---------------------------------------------------------------------
.cache <- new.env(parent = emptyenv())

sim_key <- function(s) paste(unlist(s), collapse = "|")

drug_events <- function(s, base_time = 1) {
  e <- NULL
  add <- function(x) if (is.null(e)) x else c(e, x)
  if (isTRUE(s$mido))  e <- add(ev(amt = s$mido_mg, cmt = "MIDG",
                                   time = base_time, ii = 5, addl = 2))
  if (isTRUE(s$mido_hs)) e <- add(ev(amt = s$mido_mg, cmt = "MIDG",
                                     time = 22, ii = 24, addl = 6))
  if (isTRUE(s$drox))  e <- add(ev(amt = s$drox_mg, cmt = "DRXG",
                                   time = base_time, ii = 5, addl = 2))
  if (isTRUE(s$carbi)) e <- add(ev(amt = 50, cmt = "CBDG",
                                   time = base_time, ii = 5, addl = 2))
  if (isTRUE(s$atx))   e <- add(ev(amt = s$atx_mg, cmt = "ATXG",
                                   time = base_time, ii = 12, addl = 3))
  if (isTRUE(s$fludro))e <- add(ev(amt = 0.1, cmt = "FLUG",
                                   time = base_time, ii = 24, addl = 6))
  if (isTRUE(s$pyrido))e <- add(ev(amt = 60, cmt = "PYRG",
                                   time = base_time, ii = 5, addl = 2))
  if (isTRUE(s$ddavp)) e <- add(ev(amt = 2, cmt = "DDA",
                                   time = 22, ii = 24, addl = 6))
  if (isTRUE(s$octreo))e <- add(ev(amt = 50, cmt = "OCT", time = base_time + 6.5))
  if (isTRUE(s$water)) e <- add(ev(amt = 500, cmt = "WB", time = base_time + 6))
  if (isTRUE(s$ldopa)) e <- add(ev(amt = s$ldopa_mg, cmt = "LDG",
                                   time = base_time, ii = 5, addl = 2))
  e
}

extra_params <- function(s) {
  p <- list(COQ2F = s$coq2, FTON = s$fton, MEALAMP = if (isTRUE(s$fasting)) 0 else 1)
  if (isTRUE(s$binder))  p$BINDER <- 1
  if (isTRUE(s$counter)) p$COUNTERM <- 1
  if (isTRUE(s$headup))  p$HUTF <- 0.22
  p
}

get_history <- function(s) {
  k <- paste0("hist|", sim_key(s))
  if (!is.null(.cache[[k]])) return(.cache[[k]])
  d <- msa_history(14, s$pheno, p = extra_params(s), delta = 72)
  .cache[[k]] <- d
  d
}

get_state <- function(s) {
  k <- paste0("state|", sim_key(s))
  if (!is.null(.cache[[k]])) return(.cache[[k]])
  d <- msa_state_at(s$year, s$pheno, p = extra_params(s))
  .cache[[k]] <- d
  d
}

get_tilt <- function(s, with_drug = TRUE) {
  k <- paste0("tilt|", with_drug, "|", sim_key(s))
  if (!is.null(.cache[[k]])) return(.cache[[k]])
  st <- get_state(s)
  ev0 <- if (with_drug) drug_events(s, base_time = 0) else NULL
  d <- msa_tilt(st, s$pheno, p = extra_params(s), dose = ev0)
  .cache[[k]] <- d
  d
}

get_day <- function(s, with_drug = TRUE) {
  k <- paste0("day|", with_drug, "|", sim_key(s))
  if (!is.null(.cache[[k]])) return(.cache[[k]])
  st <- get_state(s)
  ev0 <- if (with_drug) drug_events(s, base_time = 1) else NULL
  d <- msa_day(st, s$pheno, p = extra_params(s), dose = ev0, days = 2, delta = 0.1)
  .cache[[k]] <- d
  d
}

## ---------------------------------------------------------------------
##  small plotting helpers (base graphics only — no extra dependencies)
## ---------------------------------------------------------------------
.par <- function(...) par(mar = c(4.2, 4.4, 2.6, 1.2), mgp = c(2.5, 0.7, 0),
                          las = 1, cex.axis = 0.9, cex.lab = 0.95, ...)

lineplot <- function(x, ys, labels, xlab, ylab, main, cols = PAL,
                     lty = 1, ylim = NULL, hlines = NULL, vshade = NULL) {
  ys <- if (is.list(ys)) ys else list(ys)
  if (is.null(ylim)) {
    r <- range(unlist(ys)[is.finite(unlist(ys))])
    if (!all(is.finite(r)) || diff(r) == 0) r <- r[1] + c(-1, 1)
    ylim <- r + c(-0.06, 0.10)*diff(r)
  }
  .par()
  plot(NA, xlim = range(x), ylim = ylim, xlab = xlab, ylab = ylab, main = main)
  if (!is.null(vshade)) for (i in seq_len(nrow(vshade)))
    rect(vshade[i, 1], ylim[1], vshade[i, 2], ylim[2],
         col = "#f0f3f8", border = NA)
  grid(col = "#e2e8f0", lty = 1)
  if (!is.null(hlines)) abline(h = hlines, col = "#a0aec0", lty = 3)
  for (i in seq_along(ys))
    lines(x, ys[[i]], col = cols[(i - 1) %% length(cols) + 1], lwd = 2.1,
          lty = if (length(lty) > 1) lty[i] else lty)
  if (!is.null(labels))
    legend("topleft", legend = labels, col = cols[seq_along(ys)], lwd = 2.1,
           lty = if (length(lty) > 1) lty[seq_along(ys)] else lty,
           bty = "n", cex = 0.82)
  box()
  invisible(TRUE)
}

## Dosing records add extra output rows, so a drug arm and its drug-free
## comparator do NOT share a time grid.  Interpolate the comparator onto the
## drug arm's grid instead of assuming equal lengths.
onto <- function(xt, d, col) approx(d$time, d[[col]], xout = xt, rule = 2)$y

night_bands <- function(tmax) {
  b <- NULL
  d <- 0
  while (d*24 < tmax) {
    b <- rbind(b, c(d*24, d*24 + 7), c(d*24 + 22, min(d*24 + 24, tmax)))
    d <- d + 1
  }
  b[b[, 1] < tmax, , drop = FALSE]
}

## =====================================================================
##  TAB 1 — Patient profile: the SEVEN MULTIPLICANDS, shown separately
## =====================================================================
tab1_gains <- function(s) {
  st <- get_state(s)
  gc <- st$NMED^2.20 * st$NIML^1.60
  data.frame(
    term = c("압수용체 구심로\nafferent", "중추 통합\nNTS/RVLM",
             "절전 출력\nIML", "절후 NE 저장\nPOSTG",
             "α1 수용체 밀도\nA1R", "유효 순환용적\nECF",
             "정맥 용량 조절\ncapacitance"),
    value = c(1.0, st$NMED^2.20, st$NIML^1.60, st$POSTG,
              min(st$A1R/1.6, 1), min(st$ECF/14, 1.05), 1.0),
    stringsAsFactors = FALSE
  ) -> d
  d$product_gain <- c(1, gc, 1, st$POSTG, st$A1R/1.6, st$ECF/14, 1)
  attr(d, "G_CENT") <- gc
  attr(d, "state") <- st
  d
}

tab1_plot <- function(s) {
  d <- tab1_gains(s)
  .par(mar = c(7.5, 4.4, 3.0, 1.2))
  bp <- barplot(d$value, names.arg = d$term, las = 2, ylim = c(0, 1.15),
                col = ifelse(d$value < 0.5, "#c53030",
                      ifelse(d$value < 0.8, "#b7791f", "#2f855a")),
                border = NA, ylab = "잔존 이득 (1 = 정상)",
                main = sprintf("%s · 발병 후 %.1f년 — 기립 혈압의 곱셈 항 분해",
                               PH_LABEL[[s$pheno]], s$year), cex.names = 0.72)
  abline(h = c(0.5, 0.8), col = "#a0aec0", lty = 3)
  text(bp, pmin(d$value + 0.06, 1.12), sprintf("%.2f", d$value), cex = 0.8)
  mtext(sprintf("중추 이득 G_CENT = %.3f   |   절후 온전성 POSTG = %.3f   →   %s",
                attr(d, "G_CENT"), attr(d, "state")$POSTG,
                if (attr(d, "state")$POSTG > 2.5*attr(d, "G_CENT"))
                  "중추 병소 우세: NET 억제제(atomoxetine)가 작동할 조건"
                else "절후 병소 우세: 직접 α1 작용제(midodrine)만 기대 가능"),
        side = 1, line = 6.1, cex = 0.78)
  invisible(d)
}

## =====================================================================
##  TAB 2 — Orthostatic challenge (head-up tilt / active stand)
## =====================================================================
tab2_plot <- function(s) {
  tt <- get_tilt(s, with_drug = TRUE)
  t0 <- get_tilt(s, with_drug = FALSE)
  x <- (tt$time - 1.0)*60          # minutes relative to tilt
  b <- function(col) onto(tt$time, t0, col)
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  lineplot(x, list(b("SBP"), tt$SBP, b("DBP"), tt$DBP),
           c("SBP 약물 없음", "SBP 약물", "DBP 약물 없음", "DBP 약물"),
           "기립 후 시간 (분)", "혈압 (mmHg)", "기립경사 혈압 반응",
           cols = c("#2b6cb0", "#c53030", "#90a4c8", "#e08a8a"),
           lty = c(2, 1, 2, 1), hlines = 0)
  abline(v = 0, col = "#4a5568", lty = 2)
  lineplot(x, list(b("HR"), tt$HR), c("HR 약물 없음", "HR 약물"),
           "기립 후 시간 (분)", "심박수 (bpm)", "심박수 반응 (신경성 지표)",
           lty = c(2, 1))
  abline(v = 0, col = "#4a5568", lty = 2)
  lineplot(x, list(b("NEPL"), tt$NEPL), c("혈장 NE 약물 없음", "혈장 NE 약물"),
           "기립 후 시간 (분)", "혈장 NE (pg/mL)",
           "혈장 노르에피네프린 — 앙와위 정상 / 기립 증가 둔화", lty = c(2, 1))
  abline(v = 0, col = "#4a5568", lty = 2)
  lineplot(x, list(b("CBFMARG"), tt$CBFMARG), c("약물 없음", "약물"),
           "기립 후 시간 (분)", "MAP − 자동조절 하한 (mmHg)",
           "뇌혈류 자동조절 여유", lty = c(2, 1), hlines = 0)
  abline(v = 0, col = "#4a5568", lty = 2)
  layout(1)
  invisible(TRUE)
}

tab2_table <- function(s) {
  m0 <- tilt_metrics(get_tilt(s, FALSE))
  m1 <- tilt_metrics(get_tilt(s, TRUE))
  ## NOTE: column NAMES are assigned as strings, never used as identifiers.
  ## Non-ASCII identifiers (e.g. `지표 = ...`) fail to parse outside a UTF-8
  ## locale, and Rscript defaults to LC_CTYPE=C in many containers.
  out <- data.frame(
    metric = c("앙와위 SBP (mmHg)", "3분 기립 SBP (mmHg)", "ΔSBP (mmHg)",
               "ΔHR (bpm)", "ΔHR/ΔSBP (신경성: <0.5)",
               "nOH 기준 충족 (ΔSBP ≤ −20)", "앙와위 혈장 NE (pg/mL)",
               "기립 NE 증가율 (%)"),
    no_drug = c(round(m0$SBP_supine, 1), round(m0$SBP_3min, 1),
                round(m0$dSBP, 1), round(m0$dHR, 1),
                round(m0$HR_SBP_ratio, 2), ifelse(m0$nOH, "예", "아니오"),
                round(m0$NE_supine, 0), round(m0$NE_rise_pct, 1)),
    on_drug = c(round(m1$SBP_supine, 1), round(m1$SBP_3min, 1),
                round(m1$dSBP, 1), round(m1$dHR, 1),
                round(m1$HR_SBP_ratio, 2), ifelse(m1$nOH, "예", "아니오"),
                round(m1$NE_supine, 0), round(m1$NE_rise_pct, 1)),
    stringsAsFactors = FALSE)
  names(out) <- c("지표", "약물 없음", "약물")
  out
}

## =====================================================================
##  TAB 3 — 24-h BP and THE NOCTURNAL LOOP
## =====================================================================
tab3_plot <- function(s) {
  dd <- get_day(s, TRUE)
  d0 <- get_day(s, FALSE)
  nb <- night_bands(max(dd$time))
  b <- function(col) onto(dd$time, d0, col)
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  lineplot(dd$time, list(b("SBP"), dd$SBP), c("약물 없음", "약물"),
           "시간 (h, 0 = 자정)", "SBP (mmHg)",
           "48시간 SBP — 회색 = 야간 앙와위",
           lty = c(2, 1), hlines = c(140, 90), vshade = nb)
  lineplot(dd$time, list(b("UNAV"), dd$UNAV), c("약물 없음", "약물"),
           "시간 (h)", "요 나트륨 배설 (mmol/h)",
           "압력 나트륨배설 — 앙와위 고혈압이 밤에 용적을 소비한다",
           lty = c(2, 1), vshade = nb)
  lineplot(dd$time, list(b("UVOL"), dd$UVOL), c("약물 없음", "약물"),
           "시간 (h)", "요량 (L/h)", "야간 다뇨 (desmopressin 표적)",
           lty = c(2, 1), vshade = nb)
  lineplot(dd$time, list(b("ECF"), dd$ECF), c("약물 없음", "약물"),
           "시간 (h)", "세포외액량 (L)",
           "밤사이 용적 손실 → 다음 날 아침 기립성 저혈압",
           lty = c(2, 1), vshade = nb)
  layout(1)
  invisible(TRUE)
}

tab3_table <- function(s) {
  a <- day_metrics(get_day(s, FALSE))
  b <- day_metrics(get_day(s, TRUE))
  out <- data.frame(
    metric = c("야간 앙와위 SBP 최고 (mmHg)", "아침 07:30 기립 SBP (mmHg)",
               "오후 기립 SBP (mmHg)", "야간 나트륨 배설 (mmol)",
               "야간 요량 (L)", "새벽 세포외액량 (L)", "OHSA 증상점수"),
    no_drug = round(unlist(a[c(1, 2, 3, 4, 5, 6, 7)]), 2),
    on_drug = round(unlist(b[c(1, 2, 3, 4, 5, 6, 7)]), 2),
    stringsAsFactors = FALSE, row.names = NULL)
  names(out) <- c("지표", "약물 없음", "약물")
  out
}

## =====================================================================
##  TAB 4 — Drug pharmacokinetics
## =====================================================================
tab4_plot <- function(s) {
  dd <- get_day(s, TRUE)
  series <- list(); labs <- character(0)
  push <- function(v, l) { if (max(v, na.rm = TRUE) > 1e-9) {
    series[[length(series) + 1]] <<- v; labs <<- c(labs, l) } }
  push(dd$CDMM, "desglymidodrine (ng/mL)")
  push(dd$CDRX/10, "droxidopa (×10 ng/mL)")
  push(dd$CATX, "atomoxetine (ng/mL)")
  push(dd$CFLU, "fludrocortisone (ng/mL)")
  push(dd$CPYR, "pyridostigmine (ng/mL)")
  push(dd$CLDC/10, "levodopa (×10 ng/mL)")
  if (!length(series)) { series <- list(rep(0, nrow(dd))); labs <- "약물 없음" }
  layout(matrix(1:2, 2, 1))
  lineplot(dd$time, series, labs, "시간 (h)", "농도",
           "약물 혈장 농도 (48시간)", vshade = night_bands(max(dd$time)))
  lineplot(dd$time, list(dd$OCC, dd$NETI, dd$CARBI),
           c("α1 점유율 OCC", "NET 차단율 (atomoxetine)",
             "말초 AADC 차단율 (carbidopa)"),
           "시간 (h)", "분율 (0–1)",
           "표적 점유율 — carbidopa와 droxidopa는 같은 효소를 두고 경쟁한다",
           ylim = c(0, 1), vshade = night_bands(max(dd$time)))
  layout(1)
  invisible(TRUE)
}

## =====================================================================
##  TAB 5 — Autonomic axis (SNA, NE, renin, aldosterone, AVP, AQP2, A1R)
## =====================================================================
tab5_plot <- function(s) {
  dd <- get_day(s, TRUE)
  h  <- get_history(s)
  nb <- night_bands(max(dd$time))
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  lineplot(dd$time, list(dd$SNA, dd$NES, dd$NEVES),
           c("교감신경 활성 SNA", "시냅스 NE", "소포 NE 저장"),
           "시간 (h)", "정규화 값", "교감신경 구동과 NE 저장", vshade = nb)
  lineplot(dd$time, list(dd$PRA, dd$ALD, dd$MR),
           c("혈장 레닌 활성", "알도스테론", "MR 활성화"),
           "시간 (h)", "정규화 값",
           "RAAS — MSA에서 기립 시 레닌이 오르지 않는다", vshade = nb)
  lineplot(dd$time, list(dd$AVP, dd$AQP2*10),
           c("혈장 AVP (pg/mL)", "AQP2 traffic (×10)"),
           "시간 (h)", "값", "AVP–AQP2 축 (압수용체 팔이 소실)", vshade = nb)
  lineplot(h$YEARS, list(h$A1R, h$G_CENT + 1, h$POSTG + 1),
           c("α1 수용체 밀도 A1R", "중추 이득 + 1", "절후 온전성 + 1"),
           "발병 후 연수", "정규화 값",
           "탈신경 과민성: 구동이 죽는 만큼 수용체가 늘어난다")
  abline(v = s$year, col = "#c53030", lty = 2)
  layout(1)
  invisible(TRUE)
}

## =====================================================================
##  TAB 6 — Neuropathology progression (the glial cascade)
## =====================================================================
tab6_plot <- function(s) {
  h <- get_history(s)
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  lineplot(h$YEARS, list(h$GCI, h$ASYNO, h$SEED*5, h$P25A),
           c("GCI 봉입체 부담", "α-syn 올리고머", "세포외 seed ×5",
             "p25α 세포체 이동"),
           "발병 후 연수", "임의단위", "올리고덴드로글리아 α-시누클레인 캐스케이드")
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$MYE, h$CQ, h$OXS/2, h$MGL/2),
           c("수초 온전성 MYE", "CoQ10 풀", "산화 스트레스 ÷2",
             "활성 미세아교세포 ÷2"),
           "발병 후 연수", "정규화 값", "수초 · 미토콘드리아 · 신경염증")
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$NIML, h$NMED, h$POSTG),
           c("IML 절전 뉴런", "뇌간 심혈관 뉴런", "절후 종말 (보존)"),
           "발병 후 연수", "생존 분율",
           "자율신경 축: 중추는 무너지고 절후는 남는다", ylim = c(0, 1.05))
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$NMSN, h$NSN, h$NOPC, h$NONUF, h$NRESP),
           c("선조체 MSN (후시냅스)", "흑질 DA (전시냅스)",
             "올리보폰토소뇌", "Onuf 핵", "호흡 뉴런"),
           "발병 후 연수", "생존 분율", "지역별 취약성", ylim = c(0, 1.05))
  abline(v = s$year, col = "#c53030", lty = 2)
  layout(1)
  invisible(TRUE)
}

## =====================================================================
##  TAB 7 — Motor symptoms and the levodopa question
## =====================================================================
tab7_curve <- function(s, years = c(2, 4, 6, 8, 10)) {
  do.call(rbind, lapply(years, function(y) {
    s2 <- s; s2$year <- y
    st <- msa_state_at(y, s$pheno, p = extra_params(s))
    stp<- msa_state_at(y, "PD",    p = extra_params(s))
    g <- function(state, ph, mg) {
      e <- if (mg > 0) ev(amt = mg, cmt = "LDG", time = 1, ii = 5, addl = 2) else NULL
      d <- msa_day(state, ph, p = extra_params(s), dose = e, days = 1, delta = 0.5)
      mean(d$PARK)
    }
    off_m <- g(st, s$pheno, 0); on_m <- g(st, s$pheno, s$ldopa_mg)
    off_p <- g(stp, "PD", 0);   on_p <- g(stp, "PD", s$ldopa_mg)
    data.frame(year = y,
               G_post = st$NMSN^1.15,
               park_off = off_m, park_on = on_m,
               gain_pct = 100*(off_m - on_m)/max(off_m, 1e-6),
               pd_gain_pct = 100*(off_p - on_p)/max(off_p, 1e-6))
  }))
}

tab7_plot <- function(s) {
  h <- get_history(s)
  cv <- tab7_curve(s)
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  lineplot(h$YEARS, list(h$PARK, h$ATAX, h$STRID*20),
           c("파킨슨증 (UMSARS-II점)", "소뇌 실조 (UMSARS-II점)",
             "협착음 지수 ×20"),
           "발병 후 연수", "점수", "운동 표현형의 구성")
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(cv$year, list(cv$gain_pct, cv$pd_gain_pct),
           c(paste0(PH_LABEL[[s$pheno]], " 반응"), "PD 비교군 반응"),
           "발병 후 연수", "레보도파에 의한 파킨슨증 감소 (%)",
           sprintf("레보도파 %d mg t.i.d. — 같은 노출, 다른 표적", s$ldopa_mg))
  lineplot(cv$year, list(cv$G_post), "후시냅스 이득 G_post",
           "발병 후 연수", "G_post = NMSN^1.15",
           "왜 반응이 사라지는가: 표적이 사라진다", ylim = c(0, 1.05))
  lineplot(h$YEARS, list(h$STRIAT, h$DA),
           c("선조체 출력 STRIAT", "선조체 도파민 DA"),
           "발병 후 연수", "정규화 값",
           "도파민은 남아도 출력이 회복되지 않는다")
  abline(v = s$year, col = "#c53030", lty = 2)
  layout(1)
  invisible(cv)
}

## =====================================================================
##  TAB 8 — Clinical endpoints and survival
## =====================================================================
tab8_plot <- function(s) {
  h <- get_history(s)
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  lineplot(h$YEARS, list(h$U1, h$U2, h$U1 + h$U2),
           c("UMSARS-I (ADL)", "UMSARS-II (운동검사)", "UMSARS I+II"),
           "발병 후 연수", "점수", "UMSARS 진행")
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$SCOPA, h$OHSA*4, h$PVRES/10),
           c("SCOPA-AUT", "OHSA 증상 ×4", "잔뇨 ÷10 (mL)"),
           "발병 후 연수", "점수", "자율신경 증상 부담")
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$U4), "UMSARS-IV 전반 장애도",
           "발병 후 연수", "1–5", "전반적 장애도", ylim = c(1, 5.2))
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$SURV), "생존 확률",
           "발병 후 연수", "S(t)", "모델 생존 곡선", ylim = c(0, 1.02),
           hlines = 0.5)
  abline(v = s$year, col = "#c53030", lty = 2)
  layout(1)
  invisible(TRUE)
}

tab8_table <- function(s) {
  h <- get_history(s)
  at <- function(y, col) h[[col]][which.min(abs(h$YEARS - y))]
  yrs <- c(2, 4, 6, 8, 10, 12)
  out <- data.frame(year = yrs,
             UMSARS_I = round(sapply(yrs, at, "U1"), 1),
             UMSARS_II = round(sapply(yrs, at, "U2"), 1),
             SCOPA_AUT = round(sapply(yrs, at, "SCOPA"), 1),
             OHSA = round(sapply(yrs, at, "OHSA"), 2),
             PVR_mL = round(sapply(yrs, at, "PVRES"), 0),
             survival = round(sapply(yrs, at, "SURV"), 3),
             stringsAsFactors = FALSE)
  names(out) <- c("연수", "UMSARS-I", "UMSARS-II", "SCOPA-AUT", "OHSA",
                  "잔뇨 (mL)", "생존확률")
  out
}

## =====================================================================
##  TAB 9 — Scenario comparison (the therapeutic trade-off table)
## =====================================================================
SCENARIOS <- list(
  "무치료"                        = list(),
  "비약물 (복대+대항동작+수분)"   = list(binder = TRUE, counter = TRUE, water = TRUE),
  "머리높임 수면"                 = list(headup = TRUE),
  "취침 desmopressin"             = list(ddavp = TRUE),
  "머리높임 + desmopressin"       = list(headup = TRUE, ddavp = TRUE),
  "midodrine 10 mg t.i.d."        = list(mido = TRUE, mido_mg = 10),
  "midodrine 취침 투여 (오류)"    = list(mido_hs = TRUE, mido_mg = 10),
  "droxidopa 300 mg t.i.d."       = list(drox = TRUE, drox_mg = 300),
  "droxidopa + carbidopa"         = list(drox = TRUE, drox_mg = 300, carbi = TRUE),
  "atomoxetine 18 mg b.i.d."      = list(atx = TRUE, atx_mg = 18),
  "fludrocortisone 0.1 mg"        = list(fludro = TRUE),
  "pyridostigmine 60 mg t.i.d."   = list(pyrido = TRUE),
  "midodrine + fludro + 비약물"   = list(mido = TRUE, mido_mg = 10, fludro = TRUE,
                                         binder = TRUE, counter = TRUE)
)

tab9_table <- function(s, which = names(SCENARIOS)) {
  do.call(rbind, lapply(which, function(nm) {
    s2 <- modifyList(s, SCENARIOS[[nm]])
    tm <- tilt_metrics(get_tilt(s2, TRUE))
    dm <- day_metrics(get_day(s2, TRUE))
    data.frame(scenario = nm,
               stand_SBP_3min = round(tm$SBP_3min, 1),
               dSBP = round(tm$dSBP, 1),
               night_supine_peak = round(dm$SBP_supine_peak_night, 1),
               morning_stand_SBP = round(dm$SBP_morning, 1),
               night_Na_mmol = round(dm$UNa_night_mmol, 1),
               OHSA = round(dm$OHSA, 2),
               stringsAsFactors = FALSE)
  }))
}

tab9_plot <- function(s, which = names(SCENARIOS)) {
  tb <- tab9_table(s, which)
  .par(mar = c(4.2, 14, 3.0, 1.2))
  tb <- tb[order(tb$morning_stand_SBP), ]
  xr <- range(c(tb$stand_SBP_3min, tb$night_supine_peak, tb$morning_stand_SBP))
  plot(NA, xlim = xr + c(-4, 8), ylim = c(0.5, nrow(tb) + 0.5),
       yaxt = "n", xlab = "SBP (mmHg)", ylab = "",
       main = "치료 상충관계: 기립 혈압 상승 대 앙와위 고혈압")
  grid(col = "#e2e8f0", lty = 1)
  axis(2, at = seq_len(nrow(tb)), labels = tb$scenario, las = 1, cex.axis = 0.8)
  abline(v = 140, col = "#c53030", lty = 2)
  for (i in seq_len(nrow(tb))) {
    segments(tb$stand_SBP_3min[i], i, tb$night_supine_peak[i], i,
             col = "#cbd5e0", lwd = 5)
    points(tb$stand_SBP_3min[i], i, pch = 19, col = "#2b6cb0", cex = 1.2)
    points(tb$morning_stand_SBP[i], i, pch = 17, col = "#2f855a", cex = 1.1)
    points(tb$night_supine_peak[i], i, pch = 19, col = "#c53030", cex = 1.2)
  }
  legend("bottomright", c("3분 기립 SBP", "아침 07:30 기립 SBP",
                          "야간 앙와위 최고 SBP", "앙와위 고혈압 기준 140"),
         pch = c(19, 17, 19, NA), lty = c(NA, NA, NA, 2),
         col = c("#2b6cb0", "#2f855a", "#c53030", "#c53030"),
         bty = "n", cex = 0.78)
  box()
  invisible(tb)
}

## =====================================================================
##  TAB 10 — Biomarkers and the differential diagnosis panel
## =====================================================================
tab10_plot <- function(s) {
  h <- get_history(s)
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  lineplot(h$YEARS, list(h$NFL), "혈장 NfL (pg/mL)",
           "발병 후 연수", "pg/mL",
           "신경필라멘트 경쇄 — 소실 속도의 지표")
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$SEED, h$GCI/10),
           c("세포외 seed (RT-QuIC 상관)", "GCI 부담 ÷10"),
           "발병 후 연수", "임의단위", "시드 증폭 검사의 모델 상관물")
  abline(v = s$year, col = "#c53030", lty = 2)
  lineplot(h$YEARS, list(h$POSTG, h$G_CENT),
           c("절후 온전성 (¹²³I-MIBG 상관)", "중추 이득 G_CENT"),
           "발병 후 연수", "분율",
           "MIBG가 정상인데 기립성 저혈압이 심하다 = MSA", ylim = c(0, 1.05))
  abline(v = s$year, col = "#c53030", lty = 2)
  ## differential-diagnosis panel across all four phenotypes at this year
  ph <- names(pheno)
  vals <- t(sapply(ph, function(p) {
    st <- msa_state_at(s$year, p, p = extra_params(s))
    c(G_CENT = st$NMED^2.20*st$NIML^1.60, POSTG = st$POSTG, NMSN = st$NMSN)
  }))
  ## vals is phenotypes x metrics; barplot() groups by COLUMN, so pass it as
  ## is to get one group per metric with one bar per phenotype.
  .par(mar = c(5.2, 4.4, 3.0, 1.2))
  cols <- c("#2b6cb0", "#c53030", "#2f855a", "#b7791f")
  barplot(vals, beside = TRUE, ylim = c(0, 1.25), border = NA,
          col = cols[seq_len(nrow(vals))],
          names.arg = c("중추 이득 G_CENT", "절후 온전성 POSTG", "선조체 MSN"),
          ylab = "잔존 분율",
          main = sprintf("표현형 감별 (발병 후 %.0f년)", s$year))
  legend("topright", rownames(vals), fill = cols[seq_len(nrow(vals))],
         bty = "n", cex = 0.78)
  layout(1)
  invisible(vals)
}

## =====================================================================
##  HEAD-LESS SELF-TEST — walks all 10 tabs with no browser and no shiny
## =====================================================================
default_settings <- function() list(
  pheno = "MSA_P", year = 7, coq2 = 1.0, fton = 0.22, fasting = FALSE,
  mido = FALSE, mido_hs = FALSE, mido_mg = 10,
  drox = FALSE, drox_mg = 300, carbi = FALSE,
  atx = FALSE, atx_mg = 18, fludro = FALSE, pyrido = FALSE,
  ddavp = FALSE, octreo = FALSE, water = FALSE,
  binder = FALSE, counter = FALSE, headup = FALSE,
  ldopa = FALSE, ldopa_mg = 200
)

msa_shiny_selftest <- function(quick = FALSE) {
  ok <- TRUE
  png(tempfile(fileext = ".png"), width = 1100, height = 800)
  on.exit({ dev.off(); layout(1) }, add = TRUE)
  s <- default_settings()
  s$year <- 7
  s$mido <- TRUE; s$fludro <- TRUE; s$ldopa <- TRUE
  chk <- function(name, expr) {
    r <- tryCatch({ v <- expr; TRUE }, error = function(e) { message("   ", conditionMessage(e)); FALSE })
    cat(sprintf("  [%s] %s\n", if (r) "PASS" else "FAIL", name))
    if (!r) ok <<- FALSE
    invisible(r)
  }
  cat("\n--- Shiny panel self-test (head-less) ---\n")
  chk("tab 1  환자 프로파일 / 곱셈 항 분해", tab1_plot(s))
  chk("tab 2  기립경사 검사 (plot)",        tab2_plot(s))
  chk("tab 2  기립경사 검사 (table)",       stopifnot(nrow(tab2_table(s)) == 8))
  chk("tab 3  24시간 혈압 · 야간 루프",     tab3_plot(s))
  chk("tab 3  야간 루프 지표표",            stopifnot(nrow(tab3_table(s)) == 7))
  chk("tab 4  약물 PK",                     tab4_plot(s))
  chk("tab 5  자율신경 축",                 tab5_plot(s))
  chk("tab 6  신경병리 진행",               tab6_plot(s))
  chk("tab 7  운동증상과 레보도파",         tab7_plot(s))
  chk("tab 8  임상 엔드포인트 (plot)",      tab8_plot(s))
  chk("tab 8  임상 엔드포인트 (table)",     stopifnot(nrow(tab8_table(s)) == 6))
  sub <- if (quick) names(SCENARIOS)[c(1, 6, 7, 8, 9, 10)] else names(SCENARIOS)
  chk("tab 9  시나리오 비교",               tab9_plot(s, sub))
  chk("tab 10 바이오마커 · 감별",           tab10_plot(s))
  cat(sprintf("--- Shiny panel self-test %s ---\n\n", if (ok) "PASSED" else "FAILED"))
  invisible(ok)
}

## =====================================================================
##  UI
## =====================================================================
if (requireNamespace("shiny", quietly = TRUE)) {

library(shiny)

ui <- fluidPage(
  titlePanel("다계통 위축 (Multiple System Atrophy) — QSP 대시보드"),
  tags$p(style = "color:#4a5568;font-size:13px;margin-top:-8px",
    "기립 혈압을 '압수용체 구심로 × 중추 통합 × 절전 출력 × 절후 NE 저장 × α1 수용체 밀도 × 유효 순환용적 × 정맥 용량'의 곱으로 기술한 61-ODE 모델. 교육·연구용이며 임상 의사결정에 사용하지 마십시오."),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("pheno", "표현형 (vulnerability만 바꿉니다)",
                  choices = setNames(names(PH_LABEL), PH_LABEL),
                  selected = "MSA_P"),
      sliderInput("year", "발병 후 연수", min = 1, max = 12, value = 7, step = 1),
      sliderInput("coq2", "CoQ10 생합성능 (COQ2 변이 시 <1)",
                  min = 0.4, max = 1.0, value = 1.0, step = 0.05),
      sliderInput("fton", "잔여 긴장성 교감 구동 FTON (높으면 앙와위 고혈압 표현형)",
                  min = 0.10, max = 0.45, value = 0.22, step = 0.01),
      checkboxInput("fasting", "금식 프로토콜 (식후 저혈압 제거)", FALSE),
      tags$hr(),
      tags$b("승압 약물"),
      checkboxInput("mido", "midodrine t.i.d.", FALSE),
      sliderInput("mido_mg", "midodrine 용량 (mg)", 2.5, 10, 10, 2.5),
      checkboxInput("mido_hs", "midodrine 취침 투여 (잘못된 타이밍)", FALSE),
      checkboxInput("drox", "droxidopa t.i.d.", FALSE),
      sliderInput("drox_mg", "droxidopa 용량 (mg)", 100, 600, 300, 100),
      checkboxInput("carbi", "carbidopa 병용 (AADC 억제)", FALSE),
      checkboxInput("atx", "atomoxetine b.i.d.", FALSE),
      sliderInput("atx_mg", "atomoxetine 용량 (mg)", 10, 18, 18, 2),
      checkboxInput("fludro", "fludrocortisone 0.1 mg", FALSE),
      checkboxInput("pyrido", "pyridostigmine 60 mg t.i.d.", FALSE),
      checkboxInput("ddavp", "취침 desmopressin", FALSE),
      checkboxInput("octreo", "octreotide (식전)", FALSE),
      tags$hr(),
      tags$b("비약물 · 도파민성"),
      checkboxInput("binder", "복대 / 압박스타킹", FALSE),
      checkboxInput("counter", "물리적 대항동작", FALSE),
      checkboxInput("headup", "머리높임 수면 (10–30°)", FALSE),
      checkboxInput("water", "500 mL 수분 볼루스", FALSE),
      checkboxInput("ldopa", "levodopa t.i.d.", FALSE),
      sliderInput("ldopa_mg", "levodopa 용량 (mg)", 100, 400, 200, 50)
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("1 환자 프로파일", plotOutput("p1", height = "560px"),
                 tags$p(style = "font-size:12px;color:#4a5568",
                   "각 막대는 기립 혈압을 만드는 독립적인 곱셈 항입니다. 어느 항이 무너졌는지가 어떤 약이 듣는지를 결정합니다.")),
        tabPanel("2 기립경사 검사", plotOutput("p2", height = "620px"),
                 tableOutput("t2")),
        tabPanel("3 24시간 혈압·야간 루프", plotOutput("p3", height = "620px"),
                 tableOutput("t3"),
                 tags$p(style = "font-size:12px;color:#4a5568",
                   "앙와위 고혈압 → 야간 압력 나트륨배설 → 밤사이 용적 손실 → 아침 기립성 저혈압 악화. 이 순환이 승압제의 천장을 결정합니다.")),
        tabPanel("4 약물 PK", plotOutput("p4", height = "620px")),
        tabPanel("5 자율신경 축", plotOutput("p5", height = "620px")),
        tabPanel("6 신경병리 진행", plotOutput("p6", height = "620px")),
        tabPanel("7 운동증상·레보도파", plotOutput("p7", height = "620px"),
                 tableOutput("t7")),
        tabPanel("8 임상 엔드포인트", plotOutput("p8", height = "620px"),
                 tableOutput("t8")),
        tabPanel("9 시나리오 비교", plotOutput("p9", height = "560px"),
                 tableOutput("t9")),
        tabPanel("10 바이오마커·감별", plotOutput("p10", height = "620px"))
      )
    )
  )
)

server <- function(input, output, session) {
  S <- reactive({
    s <- default_settings()
    for (nm in names(s)) if (!is.null(input[[nm]])) s[[nm]] <- input[[nm]]
    s
  })
  output$p1  <- renderPlot(tab1_plot(S()))
  output$p2  <- renderPlot(tab2_plot(S()))
  output$t2  <- renderTable(tab2_table(S()))
  output$p3  <- renderPlot(tab3_plot(S()))
  output$t3  <- renderTable(tab3_table(S()))
  output$p4  <- renderPlot(tab4_plot(S()))
  output$p5  <- renderPlot(tab5_plot(S()))
  output$p6  <- renderPlot(tab6_plot(S()))
  output$p7  <- renderPlot(tab7_plot(S()))
  output$t7  <- renderTable(tab7_curve(S()))
  output$p8  <- renderPlot(tab8_plot(S()))
  output$t8  <- renderTable(tab8_table(S()))
  output$p9  <- renderPlot(tab9_plot(S()))
  output$t9  <- renderTable(tab9_table(S()))
  output$p10 <- renderPlot(tab10_plot(S()))
}

if (identical(Sys.getenv("MSA_SHINY_RUN"), "1")) shinyApp(ui, server)

} else {
  message("shiny is not installed — the compute functions and ",
          "msa_shiny_selftest() still work without it.")
}
