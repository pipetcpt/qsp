# Study 3 -- Checkpoint inhibitors: efficacy vs immune-related toxicity
# (metastatic melanoma, nivolumab + ipilimumab).
#
# Key biology: ipilimumab's contribution to anti-tumour PRIMING saturates at a
# low dose, while its contribution to autoimmune TOXICITY keeps rising
# super-linearly -- efficacy and toxicity read CTLA-4 occupancy through two
# different curves, so the nivolumab:ipilimumab dose RATIO is a handle on the
# benefit-risk frontier (what CheckMate 511 showed).
#
# Calibration:
#   CheckMate 067 (Larkin NEJM 2015/2019): monotherapy arms
#     ORR nivo 45%, ipi 19%; grade 3-4 TRAE nivo 23%, ipi 28%
#   CheckMate 511 (Lebbe JCO 2019): NIVO1/IPI3 vs NIVO3/IPI1
#     grade 3-4 48% vs 34%; ORR 51% vs 46% (n.s.)
#   Weber (JCO 2017): irAE onset timing by organ

source("R/theme_tox.R")
suppressPackageStartupMessages(library(deSolve))

nivo_occ <- function(d) d / (d + 0.3)
ipi_eff  <- function(d) d / (d + 0.30)
ipi_tox  <- function(d) d^1.8 / (d^1.8 + 2.6)

REG <- tribble(
  ~name,                 ~nivo, ~ipi, ~orr_t, ~g34_t, ~label,
  "nivolumab 3 (mono)",  3,     0,    45,     23,     "Nivolumab mono",
  "ipilimumab 3 (mono)", 0,     3,    19,     28,     "Ipilimumab mono",
  "NIVO1 / IPI3",        1,     3,    51,     48,     "NIVO1 / IPI3",
  "NIVO3 / IPI1",        3,     1,    46,     34,     "NIVO3 / IPI1")
no <- nivo_occ(REG$nivo); ie <- ipi_eff(REG$ipi); it <- ipi_tox(REG$ipi)

orr_fun <- function(p, n, e) { Tt <- p[1]*n + p[2]*e + p[3]*n*e; 75 * Tt / (p[4] + Tt) }
G34_EMAX <- 95
g34_fun <- function(p, n, x) { R <- p[1]*n + p[2]*x + p[3]*n*x; G34_EMAX * R / (p[4] + R) }

fit_o <- optim(c(1.0, 0.4, 1.2, 0.8),
               function(p) sum((orr_fun(p, no, ie) - REG$orr_t)^2),
               method = "L-BFGS-B", lower = c(0, 0, 0, 0.05), upper = c(5, 5, 5, 5))$par
w_g <- c(1, 1, 1.6, 1.6)          # weight the two combination arms (headline result)
fit_g <- optim(c(0.3, 1.0, 0.8, 0.6),
               function(p) sum(((g34_fun(p, no, it) - REG$g34_t) * w_g)^2),
               method = "L-BFGS-B", lower = c(0, 0, 0, 0.05), upper = c(5, 5, 5, 5))$par
message("ORR pred: ", paste(round(orr_fun(fit_o, no, ie), 1), collapse = " "))
message("G34 pred: ", paste(round(g34_fun(fit_g, no, it), 1), collapse = " "))

eval_reg <- function(nd, id) {
  n <- nivo_occ(nd)
  c(orr = orr_fun(fit_o, n, ipi_eff(id)), g34 = g34_fun(fit_g, n, ipi_tox(id)))
}

# ============================================================================
# FIG A -- efficacy / toxicity frontier over the dose ratio
# ============================================================================
grid <- expand.grid(nd = seq(0.3, 3, length.out = 16), id = seq(0, 3, length.out = 16)) |>
  as_tibble() |>
  mutate(v = purrr::map2(nd, id, eval_reg),
         orr = sapply(v, `[`, "orr"), g34 = sapply(v, `[`, "g34")) |> select(-v)
front <- grid |> arrange(g34) |> filter(orr > cummax(lag(orr, default = -Inf)))
regs <- REG |> mutate(v = purrr::map2(nivo, ipi, eval_reg),
                      orr = sapply(v, `[`, "orr"), g34 = sapply(v, `[`, "g34")) |> select(-v)
regcol <- c("Nivolumab mono" = BLUE, "Ipilimumab mono" = GOLD,
            "NIVO1 / IPI3" = RED, "NIVO3 / IPI1" = TEAL)
r13 <- regs |> filter(label == "NIVO1 / IPI3"); r31 <- regs |> filter(label == "NIVO3 / IPI1")
drop_pts <- r13$g34 - r31$g34

pA <- ggplot() +
  geom_point(data = grid, aes(g34, orr, fill = id), shape = 21, colour = "white", stroke = 0.2,
             size = 2.6, alpha = 0.85) +
  geom_line(data = front, aes(g34, orr, linetype = "Efficacy–safety frontier"),
            colour = NAVY, linewidth = 1.2) +
  annotate("curve", x = r13$g34 - 0.4, y = r13$orr - 0.9, xend = r31$g34 + 0.3, yend = r31$orr - 1.1,
           curvature = -0.25, colour = NAVY, linewidth = 1,
           arrow = arrow(length = unit(9, "pt"), type = "closed")) +
  annotate("label", x = 35, y = 31, hjust = 0.5,
           label = sprintf("CheckMate 511: lower ipilimumab\n\u2212%.0f pts grade 3\u20134 irAE, ORR held", drop_pts),
           colour = NAVY, fontface = "bold", size = 3.8, border.colour = NA, fill = alpha("white", 0.9),
           lineheight = 1.05) +
  geom_point(data = regs, aes(g34, orr, colour = label), size = 5.5) +
  geom_label(data = regs, aes(g34, orr, label = label, colour = label), fontface = "bold", size = 4,
             nudge_y = c(-2.4, 2.4, -2.8, 2.6), nudge_x = c(0, 0, -2.6, -1.5), border.colour = NA,
             fill = alpha("white", 0.85), label.padding = unit(1.5, "pt"), family = BASE_FAMILY) +
  scale_fill_gradientn("Ipilimumab\n(mg/kg)", colours = c("#EAF4F0", "#BFE0D2", GOLD, RED)) +
  scale_colour_manual(values = regcol, guide = "none") +
  scale_linetype_manual(NULL, values = "solid") +
  labs(x = "Grade 3–4 immune-related AE rate (%)", y = "Objective response rate (%)") +
  theme(legend.position = "right", legend.title = element_text(margin = margin(b = 8)))
save_fig(pA, "s3_pareto.png", 9.6, 4.7)

# ============================================================================
# FIG B -- dynamic tumour + effector time course (NIVO3/IPI1)
# ============================================================================
rhs <- function(t, y, p) {
  n <- nivo_occ(p[1]); i <- ipi_eff(p[2])
  drive_t <- 0.9*n + 0.4*i + 1.3*n*i
  drive_s <- 0.3*n + 1.0*i + 0.8*n*i
  dEt <- 0.20*drive_t - 0.08*y[2]
  dEs <- 0.12*drive_s - 0.05*y[3]
  dT <- 0.045*y[1]*(1 - y[1]/100) - 0.9*y[2]*y[1]/(y[1] + 20)
  list(c(dT, dEt, dEs))
}
sol <- as.data.frame(ode(c(T = 60, Et = 0, Es = 0), seq(0, 52, length.out = 400), rhs, c(3, 1)))
tc1 <- ggplot(sol, aes(time, T / 60 * 100)) +
  geom_hline(yintercept = 70, colour = GREY, linetype = "dotted") +
  annotate("text", x = 1, y = 56, label = "−30% = partial response", hjust = 0, colour = GREY, size = 3.5) +
  geom_line(colour = NAVY, linewidth = 1.2) +
  scale_y_continuous(limits = c(0, 105), breaks = c(0, 50, 100)) +
  labs(x = NULL, y = "Tumor burden\n(% baseline)", title = "Response builds first") +
  theme(axis.text.x = element_blank())
eff <- sol |> select(time, `Anti-tumor effector` = Et, `Self-reactive effector (irAE driver)` = Es) |>
  pivot_longer(-time)
tc2 <- ggplot(eff, aes(time, value, colour = name, linetype = name)) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(NULL, values = c("Anti-tumor effector" = TEAL,
                                       "Self-reactive effector (irAE driver)" = RED)) +
  scale_linetype_manual(NULL, values = c("Anti-tumor effector" = "22",
                                         "Self-reactive effector (irAE driver)" = "solid")) +
  scale_y_continuous(breaks = c(0, 2.5, 5)) +
  labs(x = "Weeks on therapy", y = "Effector pool\n(relative)", title = "Toxicity driver follows") +
  theme(legend.position = "bottom", legend.direction = "vertical")
save_fig(tc1 / tc2, "s3_timecourse.png", 6.6, 5.6)

# ============================================================================
# FIG C -- irAE onset timing by organ (Weber 2017 pattern)
# ============================================================================
org <- tribble(
  ~organ,                              ~med, ~sd, ~col,
  "Skin (rash, pruritus)",              3.5, 2.0, GOLD,
  "GI (colitis, diarrhea)",             6.5, 2.5, RED,
  "Hepatic (transaminitis)",            7.5, 2.5, PURPLE,
  "Pulmonary (pneumonitis)",            8.5, 3.5, BLUE,
  "Endocrine (thyroid, pituitary)",    10.5, 3.5, TEAL,
  "Renal (nephritis)",                 12.0, 3.5, NAVY) |>
  mutate(organ = factor(organ, levels = rev(organ)))
pC <- ggplot(org, aes(y = organ)) +
  geom_segment(aes(x = med - sd, xend = med + sd, yend = organ, colour = I(col)),
               linewidth = 7, alpha = 0.55) +
  geom_point(aes(x = med, colour = I(col)), size = 3.5) +
  geom_text(aes(x = med + sd + 0.3, label = sprintf("~%g wk", med)), hjust = 0, size = 3.6) +
  scale_x_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15)) +
  labs(x = "Weeks to irAE onset", y = NULL, title = "irAEs arrive organ by organ") +
  theme(panel.grid.major.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank())
save_fig(pC, "s3_organ_timing.png", 6.6, 5.0)

res <- list(
  regimens = setNames(lapply(seq_len(nrow(regs)), function(i)
    list(orr = round(regs$orr[i], 1), g34 = round(regs$g34[i], 1))), regs$label),
  cm511_tox_drop_pts = round(drop_pts, 1),
  cm511_orr_change_pts = round(r13$orr - r31$orr, 1))
write_json_pretty(res, "s3_results.json")
print(str(res))
message("S3 done")
