# Study 2 -- Multiple sclerosis: natalizumab dosing interval vs PML risk.
#
# Natalizumab (anti-alpha4-integrin) is among the most effective relapse
# suppressors but lifts CNS immune surveillance -> PML from JC-virus
# reactivation. Extended-interval dosing (EID, q6w) trades a small efficacy
# margin for a large PML-risk reduction -- a controllable handle.
#
# QSP layer:
#   * exposure -> alpha4 saturation over a steady-state dosing interval
#     (monoexponential decline, t1/2 11 d, Emax occupancy)
#   * time-averaged saturation -> relapse-rate suppression (ARR)
#   * stratified PML base hazard (Ho 2017) x EID hazard ratio (Ryerson 2019)
#
# Calibration anchors:
#   AFFIRM (Polman NEJM 2006): ARR 0.23 (q4w) vs 0.73 placebo (~68% reduction)
#   NOVA (Foley Lancet Neurol 2022): q6w EID retains efficacy
#   Ho (Lancet Neurol 2017), Bloomgren (NEJM 2012): PML risk table
#   Ryerson (Neurology 2019): EID PML hazard ratio (0.24 used)
#   FREEDOMS (Kappos NEJM 2010), OPERA (Hauser NEJM 2017): comparator DMTs

source("R/theme_tox.R")

NAT_THALF <- 11.0
ke <- log(2) / NAT_THALF
KD <- 2.0

sat_profile <- function(interval_d, weight_kg = 75) {
  kV <- 75 / weight_kg
  Cmax <- 30 * kV
  t <- seq(0, interval_d, length.out = 400)
  Css0 <- Cmax / (1 - exp(-ke * interval_d))   # steady-state start of interval
  C <- Css0 * exp(-ke * t)
  tibble(t = t, sat = C / (KD + C))
}
mean_sat <- function(interval_d, weight_kg = 75) mean(sat_profile(interval_d, weight_kg)$sat)

ARR0 <- 0.73
EMAX_ARR <- (1 - 0.23 / ARR0) / mean_sat(28)          # calibrate q4w -> ARR 0.23
arr_for <- function(interval_d) ARR0 * (1 - EMAX_ARR * mean_sat(interval_d))
red_for <- function(interval_d) (1 - arr_for(interval_d) / ARR0) * 100
arr_q4 <- arr_for(28); arr_q6 <- arr_for(42)
message(sprintf("EMAX_ARR=%.3f | ARR q4w=%.3f q6w=%.3f", EMAX_ARR, arr_q4, arr_q6))

# ---- PML risk stratification (per 1000 patients, by epoch of months) -------
PML <- tribble(
  ~stratum,                            ~e1,  ~e2,  ~e3,
  "Anti-JCV negative",                 0.09, 0.09, 0.09,
  "Anti-JCV+, index ≤ 0.9",       0.1,  0.3,  0.4,
  "Anti-JCV+, index 0.9–1.5",     0.3,  0.9,  3.0,
  "Anti-JCV+, index > 1.5",            1.0,  3.0,  10.0,
  "Anti-JCV+, prior immunosuppressant",1.2,  4.0,  8.0)

cum_pml <- function(stratum, horizon_m = 72, hr = 1) {
  inc <- unlist(PML[PML$stratum == stratum, c("e1", "e2", "e3")])
  hz <- -log(1 - inc / 1000) / 24                     # monthly hazard per epoch
  m <- 0:(horizon_m - 1)
  h <- hz[pmin(m %/% 24 + 1, 3)]
  1 - exp(-sum(h * hr))
}
EID_HR <- 0.24

# ============================================================================
# FIG A -- receptor desaturation + efficacy retention
# ============================================================================
prof <- bind_rows(
  sat_profile(28) |> mutate(reg = "q4w (standard)"),
  sat_profile(42) |> mutate(reg = "q6w (extended)"),
  sat_profile(56) |> mutate(reg = "q8w (over-extended)")) |>
  mutate(reg = factor(reg, levels = c("q4w (standard)", "q6w (extended)", "q8w (over-extended)")))
regcol <- c("q4w (standard)" = NAVY, "q6w (extended)" = GOLD, "q8w (over-extended)" = RED)
pA <- ggplot(prof, aes(t, sat * 100, colour = reg)) +
  geom_hline(yintercept = 50, colour = RED, linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = 1, y = 47, label = "surveillance-deficit zone", hjust = 0,
           vjust = 1, colour = RED, size = 3.5) +
  geom_line(linewidth = 1.2) +
  geom_point(data = prof |> group_by(reg) |> slice_tail(n = 1), size = 2.6) +
  scale_colour_manual(NULL, values = regcol) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Days since dose", y = "α4-integrin saturation (%)",
       title = "Receptor desaturation between doses") +
  theme(legend.position = c(0.3, 0.2))

ivs <- tibble(interval = 21:60) |> mutate(red = sapply(interval, red_for))
pts <- tibble(interval = c(28, 42), lab = c("q4w", "q6w"),
              red = sapply(c(28, 42), red_for), col = c(NAVY, GOLD))
pB <- ggplot(ivs, aes(interval, red)) +
  geom_line(colour = TEAL, linewidth = 1.2) +
  geom_point(data = pts, aes(colour = I(col)), size = 3.6) +
  geom_text(data = pts, aes(label = sprintf("%s  %.0f%%", lab, red), colour = I(col)),
            hjust = -0.1, vjust = -0.9, fontface = "bold", size = 4) +
  scale_y_continuous(limits = c(40, 80)) +
  labs(x = "Dosing interval (days)", y = "Relapse reduction vs placebo (%)",
       title = "Efficacy is retained out to q6w")
save_fig((pA | plot_spacer() | pB) + plot_layout(widths = c(1, 0.06, 1)), "s2_saturation.png", 10.6, 4.3)

# ============================================================================
# FIG B -- PML stratification heatmap
# ============================================================================
heat <- PML |> pivot_longer(e1:e3, names_to = "epoch", values_to = "risk") |>
  mutate(epoch = factor(epoch, labels = c("Months\n1–24", "Months\n25–48", "Months\n49–72")),
         stratum = factor(stratum, levels = rev(PML$stratum)))
pH <- ggplot(heat, aes(epoch, stratum, fill = risk)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = format(risk, drop0trailing = TRUE),
                colour = I(ifelse(risk > 4, "white", INK))),
            fontface = "bold", size = 4.6) +
  scale_fill_gradientn("PML / 1000", colours = c("#EAF4F0", "#BFE0D2", GOLD, RED),
                       limits = c(0, 10)) +
  scale_x_discrete(position = "top") +
  labs(x = "Cumulative time on natalizumab", y = NULL) +
  theme(panel.grid = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12),
        legend.title = element_text(margin = margin(b = 8)), legend.key.height = unit(22, "pt"))
save_fig(pH, "s2_pml_table.png", 10.0, 4.6)

# ============================================================================
# FIG C -- benefit-risk over a 6-year horizon (dumbbell)
# ============================================================================
br_strata <- PML$stratum[2:5]
br_lab <- c("Index ≤ 0.9", "Index 0.9–1.5", "Index > 1.5", "Prior immunosuppressant")
br <- tibble(stratum = br_strata, lab = factor(br_lab, levels = rev(br_lab))) |>
  mutate(sid = sapply(stratum, cum_pml, hr = 1) * 1000,
         eid = sapply(stratum, cum_pml, hr = EID_HR) * 1000,
         drop = (1 - eid / sid) * 100)
pC <- ggplot(br, aes(y = lab)) +
  geom_segment(aes(x = eid, xend = sid, yend = lab), colour = LGREY, linewidth = 4, lineend = "round") +
  geom_point(aes(x = sid, colour = "Standard q4w"), size = 4.8) +
  geom_point(aes(x = eid, colour = "Extended q6w"), size = 4.8) +
  geom_text(aes(x = sid, label = sprintf("%.1f", sid)), hjust = -0.6, colour = RED, fontface = "bold", size = 4) +
  geom_text(aes(x = eid, label = sprintf("%.1f", eid)), hjust = 1.6, colour = TEAL, fontface = "bold", size = 4) +
  geom_text(aes(x = (sid + eid) / 2, label = sprintf("−%.0f%%", drop)), vjust = -1.3,
            colour = NAVY, fontface = "italic", size = 3.7) +
  scale_colour_manual("Natalizumab", values = c("Standard q4w" = RED, "Extended q6w" = TEAL),
                      breaks = c("Standard q4w", "Extended q6w")) +
  scale_x_continuous(limits = c(-1.5, 17)) +
  labs(x = "6-year cumulative PML risk (per 1000 patients)", y = "Anti-JCV+ stratum") +
  theme(panel.grid.major.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 12), legend.position = "top",
        legend.justification = "right")
save_fig(pC, "s2_benefit_risk.png", 9.4, 4.4)

# ============================================================================
# FIG D -- cross-DMT benefit-risk map (illustrative composite safety index)
# ============================================================================
DMT <- tribble(
  ~name,               ~eff, ~saf, ~col,
  "Interferon-β / GA", 32,  2,   GREY,
  "Fingolimod",          54,  9,   BLUE,
  "Ocrelizumab",         68,  11,  PURPLE,
  "Natalizumab q4w",     68,  24,  RED,
  "Natalizumab q6w",     63,  6,   TEAL) |>
  mutate(lx = c(2.8, 9.8, 10.2, 25.2, 6.8), ly = c(32, 54, 70.4, 68, 60.2), hj = c(0, 0, 1, 0, 0))
pD <- ggplot(DMT, aes(saf, eff)) +
  annotate("segment", x = 23.2, y = 67.6, xend = 6.9, yend = 63.3, colour = NAVY, linewidth = 1,
           arrow = arrow(length = unit(9, "pt"), type = "closed")) +
  annotate("text", x = 17.3, y = 71, label = "Extended-interval dosing", colour = NAVY,
           fontface = "bold", size = 3.9) +
  geom_point(aes(colour = I(col)), size = 6.5) +
  geom_text(aes(lx, ly, label = name, colour = I(col), hjust = hj), fontface = "bold", size = 4) +
  annotate("text", x = 0.6, y = 79, label = "Better: more efficacy, less risk", hjust = 0,
           colour = GREEN, fontface = "italic", size = 3.6) +
  scale_x_continuous(limits = c(0, 34), breaks = c(0, 10, 20, 30)) + scale_y_continuous(limits = c(25, 80)) +
  labs(x = "Serious-safety-event index (illustrative, per 1000 patient-years)",
       y = "Relapse-rate reduction vs placebo (%)")
save_fig(pD, "s2_dmt_map.png", 9.0, 4.6)

res <- list(
  arr_placebo = ARR0, arr_q4w = round(arr_q4, 3), arr_q6w = round(arr_q6, 3),
  red_q4w_pct = round(red_for(28), 1), red_q6w_pct = round(red_for(42), 1),
  eid_pml_hr = EID_HR,
  relapse_cost_6y = round((arr_q6 - arr_q4) * 6, 2),
  pml_6y_highrisk_sid = round(cum_pml("Anti-JCV+, index > 1.5", 72, 1) * 1000, 1),
  pml_6y_highrisk_eid = round(cum_pml("Anti-JCV+, index > 1.5", 72, EID_HR) * 1000, 1),
  pml_6y_priorIS_sid = round(cum_pml("Anti-JCV+, prior immunosuppressant", 72, 1) * 1000, 1),
  pml_6y_priorIS_eid = round(cum_pml("Anti-JCV+, prior immunosuppressant", 72, EID_HR) * 1000, 1),
  pml_drop_pct = round((1 - cum_pml("Anti-JCV+, index > 1.5", 72, EID_HR) /
                          cum_pml("Anti-JCV+, index > 1.5", 72, 1)) * 100))
write_json_pretty(res, "s2_results.json")
print(str(res))
message("S2 done")
