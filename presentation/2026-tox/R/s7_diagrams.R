# Schematic + NAM-bridge figures, drawn in ggplot2.
#   figures/architecture.png   -- the daily autonomous build loop
#   figures/glk_validation.png -- gut-liver-kidney MPS + whole-body PBPK:
#                                 per-drug AUC GMFE (6 drugs, 24 clinical
#                                 observations, no clinical PK used to fit)

source("R/theme_tox.R")

# ============================================================================
# Architecture: literature -> LLM agent -> 4 artifacts -> gates -> git
# ============================================================================
box <- function(x, y, w, h, label, fill, col, txt, size = 4.2, face = "bold") {
  list(rect = tibble(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2,
                     fill = fill, col = col),
       text = tibble(x = x, y = y, label = label, txt = txt, size = size, face = face))
}
B <- list(
  box(1.1, 3, 1.75, 1.25, "Disease list\n+ PubMed", "#F4F6F9", GREY, INK),
  box(3.55, 3, 2.0, 1.6, "Claude Code\nRoutine\n(LLM agent)", NAVY, NAVY, "white", 4.6),
  box(6.5, 5.1, 2.45, 0.86, "Mechanistic map\n≥ 100 nodes", "white", TEAL, NAVY, 3.9),
  box(6.5, 3.9, 2.45, 0.86, "ODE model\n≥ 15 states", "white", TEAL, NAVY, 3.9),
  box(6.5, 2.7, 2.45, 0.86, "Shiny app\n≥ 6 tabs", "white", TEAL, NAVY, 3.9),
  box(6.5, 1.5, 2.45, 0.86, "References\n≥ 30 PubMed", "white", TEAL, NAVY, 3.9),
  box(9.25, 3, 1.75, 1.25, "Quality\ngates", "#FFF6E6", GOLD, INK, 4.1),
  box(11.55, 3, 1.75, 1.25, "Git commit\n& push", "#F4F6F9", GOLD, NAVY, 4.1))
rects <- bind_rows(lapply(B, `[[`, "rect")); texts <- bind_rows(lapply(B, `[[`, "text"))
arrows <- tribble(
  ~x,    ~y,  ~xend, ~yend, ~col,
  1.98,  3,   2.52,  3,     TEAL,
  4.58,  3.3, 5.25,  5.05,  TEAL,
  4.58,  3.1, 5.25,  3.9,   TEAL,
  4.58,  2.9, 5.25,  2.72,  TEAL,
  4.58,  2.7, 5.25,  1.55,  TEAL,
  7.75,  5.05, 8.35, 3.35,  GOLD,
  7.75,  3.9, 8.35,  3.12,  GOLD,
  7.75,  2.72, 8.35, 2.9,   GOLD,
  7.75,  1.55, 8.35, 2.65,  GOLD,
  10.15, 3,   10.65, 3,     GOLD)
pArch <- ggplot() +
  geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend, colour = I(col)), linewidth = 0.9,
               arrow = arrow(length = unit(8, "pt"), type = "closed")) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                              fill = I(fill), colour = I(col)), linewidth = 0.9) +
  geom_text(data = texts, aes(x, y, label = label, colour = I(txt), size = I(size), fontface = face),
            lineheight = 1.05, family = BASE_FAMILY) +
  geom_path(data = tibble(x = c(11.55, 11.55, 3.55, 3.55), y = c(2.37, 0.72, 0.72, 2.18)),
            aes(x, y), colour = GREY, linewidth = 0.8, linetype = "22",
            arrow = arrow(length = unit(7, "pt"), type = "closed")) +
  annotate("label", x = 9.25, y = 0.72, label = "Next day: next disease", colour = INK,
           fontface = "italic", size = 3.9, family = BASE_FAMILY, fill = "white", border.colour = NA) +
  coord_cartesian(xlim = c(0.1, 12.5), ylim = c(0.25, 5.6), expand = FALSE) +
  theme_void(base_family = BASE_FAMILY) +
  theme(plot.background = element_rect(fill = "white", colour = NA))
save_fig(pArch, "architecture.png", 11.0, 4.6, scale = 0.85)

# ============================================================================
# GLK MPS + PBPK validation (per-drug AUC GMFE)
# ============================================================================
glk <- tribble(
  ~drug,              ~mech,                        ~gmfe, ~in2, ~in4,
  "Antipyrine",       "CYP1A2, low extraction",     1.59,  4,    4,
  "Benzylpenicillin", "OAT1/3 renal secretion",     1.23,  4,    4,
  "Testosterone",     "CYP3A4 high-E + SHBG",       1.98,  2,    4,
  "Crizotinib",       "CYP3A4 high-E + biliary",    2.61,  1,    4,
  "Gefitinib",        "CYP3A4 + biliary-dominant",  2.93,  0,    4,
  "Simvastatin",      "OATP1B1 + CYP3A4",           1.37,  4,    4) |>
  arrange(gmfe) |>
  mutate(lab = factor(drug, levels = rev(drug)),
         col = ifelse(gmfe <= 2, TEAL, GOLD))
pG1 <- ggplot(glk, aes(y = lab)) +
  annotate("rect", xmin = 1, xmax = 2, ymin = -Inf, ymax = Inf, fill = GREEN, alpha = 0.08) +
  annotate("rect", xmin = 2, xmax = 3.15, ymin = -Inf, ymax = Inf, fill = GOLD, alpha = 0.07) +
  geom_vline(xintercept = 1.85, colour = NAVY, linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 1.88, y = 0.45, label = "Overall 1.85×", hjust = 0, colour = NAVY,
           fontface = "bold", size = 3.6) +
  geom_segment(aes(x = 1, xend = gmfe, yend = lab, colour = I(col)), linewidth = 1.6) +
  geom_point(aes(x = gmfe, colour = I(col)), size = 4.6) +
  annotate("text", x = 1.03, y = 6.75, hjust = 0, label = "≤ 2-fold", colour = GREEN, size = 3.3, fontface = "italic") +
  annotate("text", x = 2.05, y = 6.75, hjust = 0, label = "2–4-fold", colour = GOLD, size = 3.3, fontface = "italic") +
  scale_x_continuous(limits = c(1, 3.15), breaks = 1:3, labels = paste0(1:3, "×"),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(add = c(1.0, 1.0))) +
  labs(x = "AUC GMFE (4 observations per drug)", y = NULL,
       title = "Human AUC predicted from chip data alone") +
  theme(panel.grid.major.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold", colour = NAVY))
pT <- ggplot(glk, aes(y = lab)) +
  geom_text(aes(x = 0, label = sprintf("%.2f×", gmfe)), hjust = 0, size = 3.8, fontface = "bold") +
  geom_text(aes(x = 1.15, label = paste0(in2, "/4")), hjust = 0, size = 3.8) +
  geom_text(aes(x = 2.25, label = mech), hjust = 0, size = 3.6, colour = "#46525E") +
  annotate("text", x = c(0, 1.15, 2.25), y = 6.75, label = c("GMFE", "2-fold", "Mechanism"),
           hjust = 0, size = 3.4, fontface = "bold", colour = NAVY) +
  scale_x_continuous(limits = c(0, 7.6), expand = expansion(0)) +
  scale_y_discrete(expand = expansion(add = c(1.0, 1.0))) +
  labs(title = " ") +
  theme_void(base_family = BASE_FAMILY) +
  theme(plot.title = element_text(size = 14))
corr <- tibble(rule = factor(c("Uniform\ncorrection", "Selective\nrule"),
                             levels = c("Uniform\ncorrection", "Selective\nrule")),
               gmfe = c(9.7, 1.85), col = c(RED, TEAL))
pG2 <- ggplot(corr, aes(rule, gmfe, fill = I(col))) +
  geom_col(width = 0.56) +
  geom_text(aes(label = paste0(c("9.7", "1.85"), "×")), vjust = -0.5, fontface = "bold", size = 4.4) +
  scale_y_continuous(limits = c(0, 11.5), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Overall AUC GMFE", title = "Scaling-rule choice") +
  theme(panel.grid.major.x = element_blank())
save_fig((pG1 | pT | plot_spacer() | pG2) + plot_layout(widths = c(1.15, 1.6, 0.06, 0.8)),
         "glk_validation.png", 12.2, 4.6)

message("S7 done")
