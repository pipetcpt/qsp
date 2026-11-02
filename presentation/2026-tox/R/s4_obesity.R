# Study 4 -- Obesity: incretin titration as a three-way optimization of
# weight loss (benefit), GI tolerability and lean-mass loss (safety), with an
# elderly / sarcopenic subgroup.
#
# QSP layer (daily time step):
#   * exposure equilibrates toward the steady state of the current dose
#   * appetite suppression = Emax * occupancy -> caloric intake
#   * energy-balance ODE; weight change partitioned into fat / fat-free mass
#     by a Forbes-type rule with age-dependent anabolic resistance
#   * GI AE = nausea driven by the RATE and level of occupancy, blunted by a
#     tolerance state that builds with sustained exposure
#
# Calibration anchors:
#   STEP 1 (Wilding NEJM 2021): semaglutide 2.4 mg -> -14.9% at 68 wk
#   SURMOUNT-1 (Jastreboff NEJM 2022): tirzepatide 5/10/15 mg -> -15/-19.5/-20.9%

source("R/theme_tox.R")

DRUGS <- list(
  semaglutide = list(Emax = 0.285, EC50 = 0.9, maxdose = 2.4, color = TEAL,
                     steps = c(0.25, 0.5, 1.0, 1.7, 2.4)),
  tirzepatide = list(Emax = 0.40, EC50 = 6.0, maxdose = 15.0, color = NAVY,
                     steps = c(2.5, 5.0, 7.5, 10.0, 12.5, 15.0)))
KE <- 1 / 7

simulate <- function(drug, target, age = 45, weeks = 72, weeks_step = 4, BW0 = 105) {
  d <- DRUGS[[drug]]
  steps <- d$steps[d$steps <= target + 1e-9]
  if (!length(steps)) steps <- d$steps[1]
  if (tail(steps, 1) < target) steps <- c(steps, target)
  dose_fn <- function(week) steps[pmin(floor(week / weeks_step) + 1, length(steps))]

  fat_frac0 <- 0.42 + 0.0015 * (age - 45)
  FM0 <- BW0 * fat_frac0; FFM0 <- BW0 - FM0
  EI0 <- 22 * FFM0 + 500
  t <- 0:(weeks * 7)
  n <- length(t)
  dose_series <- dose_fn(t / 7)
  C <- numeric(n)
  for (k in 2:n) C[k] <- C[k - 1] + KE * (dose_series[k] - C[k - 1])
  occ <- C / (C + d$EC50)
  supp <- d$Emax * occ

  FM <- FFM <- BW <- numeric(n)
  FM[1] <- FM0; FFM[1] <- FFM0; BW[1] <- BW0
  age_mult <- 1 + 0.018 * (age - 45)
  for (k in 2:n) {
    EI <- EI0 * (1 - supp[k])
    EE <- 22 * FFM[k - 1] + 500 * (BW[k - 1] / BW0) - 120 * supp[k - 1]
    bal <- EI - EE
    p_lean <- min(0.6, 20 / (20 + FM[k - 1]) * age_mult)
    rho <- 9441 * (1 - p_lean) + 1816 * p_lean
    dBW <- bal / rho
    FFM[k] <- FFM[k - 1] + dBW * p_lean
    FM[k] <- FM[k - 1] + dBW * (1 - p_lean)
    BW[k] <- FFM[k] + FM[k]
  }

  docc <- np_gradient(occ, t)
  N <- Tol <- numeric(n)
  for (k in 2:n) {
    drive <- 22 * max(docc[k], 0) + 0.9 * occ[k]
    Tol[k] <- Tol[k - 1] + (0.045 * occ[k] * (1 - Tol[k - 1]) - 0.004 * Tol[k - 1])
    N[k] <- N[k - 1] + (drive * (1 - 0.85 * Tol[k - 1]) - 0.22 * N[k - 1])
  }
  esc <- t <= 24 * 7
  list(t = t / 7, supp = supp, dose = dose_series, N = N, BW = BW,
       FFM0 = FFM0, lean_lost = FFM0 - FFM[n], fat_lost = FM0 - FM[n],
       final_pct = (BW0 - BW[n]) / BW0 * 100, gi_burden = trapz(t[esc], N[esc]))
}

for (x in list(c("semaglutide", 2.4), c("tirzepatide", 5), c("tirzepatide", 10), c("tirzepatide", 15))) {
  r <- simulate(x[1], as.numeric(x[2]))
  message(sprintf("%s %s mg -> %.1f%% weight, lean %.1f kg (%.0f%% of loss)", x[1], x[2],
                  r$final_pct, r$lean_lost, 100 * r$lean_lost / (r$lean_lost + r$fat_lost)))
}

# ============================================================================
# FIG A -- titration: appetite suppression with dose steps; nausea by schedule
# ============================================================================
r_slow <- simulate("semaglutide", 2.4, weeks_step = 4)
r_fast <- simulate("semaglutide", 2.4, weeks_step = 2)
SC <- 40 / 3                        # 0-3 mg mapped onto 0-40 %
da <- tibble(t = r_slow$t, supp = r_slow$supp * 100, dose = r_slow$dose)
pA <- ggplot(da, aes(t)) +
  geom_step(aes(y = dose * SC, colour = "Dose (mg, right axis)"), linewidth = 0.8) +
  geom_line(aes(y = supp, colour = "Appetite suppression"), linewidth = 1.3) +
  scale_colour_manual(NULL, values = c("Appetite suppression" = TEAL, "Dose (mg, right axis)" = GREY)) +
  scale_y_continuous("Appetite suppression (%)", limits = c(0, 40),
                     sec.axis = sec_axis(~ . / SC, name = "Dose (mg)")) +
  scale_x_continuous(limits = c(0, 72)) +
  labs(x = "Weeks", title = "Exposure builds with each dose step") +
  theme(legend.position = c(0.6, 0.2), legend.background = element_blank(),
        axis.title.y.right = element_text(colour = GREY, angle = 90, margin = margin(l = 6)),
        axis.text.y.right = element_text(colour = GREY))
nref <- max(r_slow$N, r_fast$N)
dn <- bind_rows(tibble(t = r_slow$t, N = r_slow$N / nref * 100, sch = "4-weekly escalation"),
                tibble(t = r_fast$t, N = r_fast$N / nref * 100, sch = "2-weekly (too fast)")) |>
  filter(t <= 40)
pB <- ggplot(dn, aes(t, N, colour = sch)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(NULL, values = c("4-weekly escalation" = GOLD, "2-weekly (too fast)" = RED),
                      breaks = c("4-weekly escalation", "2-weekly (too fast)")) +
  labs(x = "Weeks", y = "Nausea intensity (relative %)", title = "Nausea follows escalation") +
  theme(legend.position = c(0.7, 0.85))
save_fig((pA | plot_spacer() | pB) + plot_layout(widths = c(1, 0.06, 1)), "s4_titration.png", 10.8, 4.3)

# ============================================================================
# FIG B -- efficacy vs tolerability frontier over target dose
# ============================================================================
sweep <- bind_rows(lapply(names(DRUGS), function(dr) {
  md <- DRUGS[[dr]]$maxdose
  bind_rows(lapply(seq(md * 0.15, md, length.out = 10), function(tg) {
    r <- simulate(dr, tg); tibble(drug = dr, target = tg, eff = r$final_pct, gi = r$gi_burden)
  }))
}))
NREF <- max(sweep$gi)
sweep <- sweep |> mutate(gi = gi / NREF * 100)
appr <- tribble(~drug, ~target, ~lab,
                "semaglutide", 2.4, "Sema 2.4 mg",
                "tirzepatide", 5, "Tirz 5 mg",
                "tirzepatide", 10, "Tirz 10 mg",
                "tirzepatide", 15, "Tirz 15 mg") |>
  mutate(r = purrr::map2(drug, target, simulate),
         eff = sapply(r, `[[`, "final_pct"), gi = sapply(r, `[[`, "gi_burden") / NREF * 100) |>
  select(-r)
dcol <- c(semaglutide = TEAL, tirzepatide = NAVY)
pF <- ggplot(sweep, aes(gi, eff, colour = drug)) +
  geom_line(linewidth = 1.2) +
  geom_point(aes(shape = drug), size = 2.6) +
  geom_point(data = appr, size = 5.5, colour = "white") +
  geom_point(data = appr, size = 4.2) +
  geom_text(data = appr, aes(label = lab), hjust = -0.25, vjust = 1.1, fontface = "bold", size = 3.9,
            show.legend = FALSE) +
  scale_colour_manual("Incretin", values = dcol, labels = c("Semaglutide", "Tirzepatide")) +
  scale_shape_manual("Incretin", values = c(16, 15), labels = c("Semaglutide", "Tirzepatide")) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.2))) +
  coord_cartesian(clip = "off") +
  labs(x = "Cumulative GI adverse-event burden, escalation phase (relative)",
       y = "Body-weight loss at 72 weeks (%)") +
  theme(legend.position = c(0.86, 0.18), legend.background = element_rect(fill = "white", colour = NA))
save_fig(pF, "s4_frontier.png", 9.0, 4.6)

# ============================================================================
# FIG C -- body composition: fat vs lean, young vs elderly
# ============================================================================
comp <- function(lab, age, preserve = FALSE) {
  r <- simulate("tirzepatide", 15, age = age)
  fat <- r$fat_lost; lean <- r$lean_lost
  if (preserve) { shift <- lean * 0.45; lean <- lean - shift; fat <- fat + shift * 0.4 }
  tibble(group = lab, part = c("Fat mass lost", "Lean mass lost"), kg = c(fat, lean))
}
cd <- bind_rows(comp("Age 45", 45), comp("Age 70", 70), comp("Age 70 +\nmuscle-preserving", 70, TRUE)) |>
  mutate(group = factor(group, levels = unique(group)),
         part = factor(part, levels = c("Lean mass lost", "Fat mass lost")))
tops <- cd |> group_by(group) |>
  summarise(tot = sum(kg), lean = kg[part == "Lean mass lost"], .groups = "drop")
pC <- ggplot(cd, aes(group, kg, fill = part)) +
  geom_col(width = 0.58) +
  geom_text(aes(label = sprintf("%.1f", kg)), position = position_stack(vjust = 0.5),
            colour = "white", fontface = "bold", size = 4) +
  geom_text(data = tops, aes(group, tot, label = sprintf("%.0f%% lean", 100 * lean / tot)),
            inherit.aes = FALSE, vjust = -0.6, size = 3.8) +
  scale_fill_manual(NULL, values = c("Fat mass lost" = GOLD, "Lean mass lost" = RED),
                    breaks = c("Fat mass lost", "Lean mass lost")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Mass lost at 72 weeks (kg)", title = "Same drug, different composition cost") +
  theme(legend.position = "top", legend.justification = "left", panel.grid.major.x = element_blank())
save_fig(pC, "s4_composition.png", 5.6, 4.6)

# ============================================================================
# FIG D -- lean share of weight lost across age (the elderly subgroup)
# ============================================================================
ages <- seq(40, 85, by = 3)
lean_share <- function(dr, tg, a) { r <- simulate(dr, tg, age = a); 100 * r$lean_lost / (r$lean_lost + r$fat_lost) }
ag <- bind_rows(
  tibble(age = ages, drug = "Tirzepatide 15 mg", share = sapply(ages, function(a) lean_share("tirzepatide", 15, a))),
  tibble(age = ages, drug = "Semaglutide 2.4 mg", share = sapply(ages, function(a) lean_share("semaglutide", 2.4, a))))
mk <- tibble(age = c(45, 70), share = sapply(c(45, 70), function(a) lean_share("tirzepatide", 15, a)))
pD <- ggplot(ag, aes(age, share, colour = drug)) +
  geom_line(linewidth = 1.2) + geom_point(size = 1.8) +
  geom_point(data = mk, aes(age, share), inherit.aes = FALSE, colour = RED, size = 4) +
  geom_text(data = mk, aes(age, share, label = sprintf("%.0f%%", share)), inherit.aes = FALSE,
            colour = RED, fontface = "bold", size = 4, vjust = -1.1) +
  scale_colour_manual(NULL, values = c("Tirzepatide 15 mg" = NAVY, "Semaglutide 2.4 mg" = TEAL)) +
  scale_y_continuous(limits = c(25, 62)) +
  labs(x = "Patient age (years)", y = "Lean share of weight lost (%)",
       title = "Lean share of loss climbs with age") +
  theme(legend.position = c(0.72, 0.14), legend.background = element_rect(fill = "white", colour = NA))
save_fig(pD, "s4_aging.png", 5.6, 4.6)

y45 <- simulate("tirzepatide", 15, age = 45); y70 <- simulate("tirzepatide", 15, age = 70)
res <- list(
  sema24_pct = round(simulate("semaglutide", 2.4)$final_pct, 1),
  tirze5_pct = round(simulate("tirzepatide", 5)$final_pct, 1),
  tirze10_pct = round(simulate("tirzepatide", 10)$final_pct, 1),
  tirze15_pct = round(y45$final_pct, 1),
  tirze15_lean_young = round(y45$lean_lost, 1), tirze15_lean_old = round(y70$lean_lost, 1),
  lean_frac_young = round(100 * y45$lean_lost / (y45$lean_lost + y45$fat_lost)),
  lean_frac_old = round(100 * y70$lean_lost / (y70$lean_lost + y70$fat_lost)))
write_json_pretty(res, "s4_results.json")
print(str(res))
message("S4 done")
