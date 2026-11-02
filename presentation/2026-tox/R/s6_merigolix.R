# Flagship -- Merigolix (oral GnRH antagonist, endometriosis): the
# estrogen-threshold window as a dose-optimization problem.
#
# Data: outputs of the merigolix QSP model (PK 2-compartment + HPG-axis
# indirect-response PD with feedback; 24 states; virtual population, week 12),
# built on the Phase 1 SAD/MAD PK/PD model (n = 51; Bae, Kang et al., Clin PK,
# in revision) -- github.com/pipetcpt/merigolix-qsp.
#   Week-12 E2 distribution by dose, % in the 20-40 pg/mL window,
#   % below 20 pg/mL (full suppression -> BMD / vasomotor risk),
#   and the composite therapeutic index
#   TI = 0.4 x efficacy + 0.3 x safety + 0.2 x pain + 0.1 x (100 - hot flashes).

source("R/theme_tox.R")

e2 <- tribble(
  ~dose,     ~mean, ~median, ~p10, ~p90, ~in_window, ~below20,
  "PBO",     58.2,  52.1,    28.3, 94.8, 18.4,       2.1,
  "80",      47.3,  42.4,    24.1, 76.3, 31.2,       8.3,
  "120",     35.6,  31.1,    17.8, 56.4, 45.6,       22.4,
  "160",     27.4,  23.4,    12.1, 44.2, 42.1,       41.8,
  "240",     17.8,  14.2,    6.8,  30.1, 26.3,       68.4) |>
  mutate(dose = factor(dose, levels = dose))
ti <- tibble(dose = factor(c("80", "120", "160", "240", "320"), levels = c("80", "120", "160", "240", "320")),
             ti = c(58.4, 72.1, 78.3, 74.2, 68.6))

pA <- ggplot(e2, aes(dose)) +
  scale_x_discrete() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 20, ymax = 40, fill = PALEGREEN) +
  geom_linerange(aes(ymin = p10, ymax = p90), colour = NAVY, linewidth = 1, alpha = 0.6) +
  geom_point(aes(y = median), colour = NAVY, size = 3.6) +
  labs(x = "Merigolix (mg/day)", y = "Estradiol at week 12 (pg/mL)",
       title = "E2 at week 12", subtitle = "Median, P10–P90 · band = window") +
  theme(plot.subtitle = element_text(colour = GREEN, face = "bold"))

wb <- e2 |> select(dose, `In 20–40 window` = in_window, `Below 20 (full suppression)` = below20) |>
  pivot_longer(-dose)
pB <- ggplot(wb, aes(dose, value, colour = name, group = name)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3) +
  scale_colour_manual(NULL, values = c("In 20–40 window" = TEAL, "Below 20 (full suppression)" = RED)) +
  scale_y_continuous(limits = c(0, 75)) +
  labs(x = "Merigolix (mg/day)", y = "Virtual patients (%)", title = "In window vs over-suppressed") +
  theme(legend.position = "bottom", legend.direction = "vertical")

pC <- ggplot(ti, aes(dose, ti, fill = I(ifelse(dose == "160", NAVY, LGREY)))) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.1f", ti), colour = I(ifelse(dose == "160", NAVY, INK))),
            vjust = -0.5, fontface = "bold", size = 3.8) +
  scale_y_continuous(limits = c(0, 90), expand = expansion(mult = c(0, 0))) +
  labs(x = "Merigolix (mg/day)", y = "Therapeutic index score", title = "Therapeutic index", subtitle = "Composite benefit–risk score") +
  theme(panel.grid.major.x = element_blank())

save_fig((pA | plot_spacer() | pB | plot_spacer() | pC) +
           plot_layout(widths = c(1, 0.08, 1, 0.08, 1)), "s6_merigolix.png", 12.4, 4.7)

res <- list(best_dose = 160, best_ti = 78.3, in_window_120 = 45.6, in_window_160 = 42.1,
            below20_160 = 41.8, below20_240 = 68.4)
write_json_pretty(res, "s6_results.json")
message("S6 done")
