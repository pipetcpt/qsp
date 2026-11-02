# Shared ggplot2 style for the 2026-tox talk figures.
#
# Run every study script from presentation/2026-tox, e.g.
#   Rscript R/s2_ms.R
# Figures are written to figures/ at 300 dpi on a white canvas, sized to the
# physical box they occupy on the slide so on-slide text stays >= ~11 pt.
# All figure text is English only (deck language rule).

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(jsonlite)
  library(scales)
  library(patchwork)
  library(ggrepel)
})

if (!file.exists("R/theme_tox.R")) {
  stop("Run from presentation/2026-tox (e.g. Rscript R/s2_ms.R)")
}
dir.create("figures", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

# ---- typeface: S-Core Dream (matches the slide body font) ------------------
.font_dir <- path.expand("~/Library/Fonts")
.reg <- file.path(.font_dir, "SCDream4.otf")
.bold <- file.path(.font_dir, "SCDream6.otf")
if (file.exists(.reg) && file.exists(.bold)) {
  systemfonts::register_font("Dream", plain = .reg, bold = .bold,
                             italic = .reg, bolditalic = .bold)
  BASE_FAMILY <- "Dream"
} else {
  BASE_FAMILY <- "sans"
}

# ---- palette (benefit vs risk) ----------------------------------------------
NAVY   <- "#163A5F"   # structure / titles
BLUE   <- "#2B6CB0"
TEAL   <- "#1F8A8C"   # benefit / efficacy
GOLD   <- "#C8922A"   # caution
RED    <- "#C0392B"   # harm / toxicity
GREEN  <- "#2E8B6F"   # safe / target window
PURPLE <- "#6B4E8C"
GREY   <- "#6B7785"
LGREY  <- "#D7DEE5"
INK    <- "#1A222C"
PALEGREEN <- "#E4F1EB"
CYCLE  <- c(NAVY, TEAL, GOLD, RED, BLUE, GREEN, PURPLE, GREY)

theme_tox <- function(base_size = 13) {
  theme_minimal(base_size = base_size, base_family = BASE_FAMILY) %+replace%
    theme(
      plot.title = element_text(face = "bold", colour = NAVY, size = rel(1.08),
                                hjust = 0, margin = margin(b = 8)),
      plot.subtitle = element_text(colour = GREY, size = rel(0.9), hjust = 0,
                                   margin = margin(b = 8)),
      axis.title = element_text(colour = INK, size = rel(0.95)),
      axis.title.y = element_text(colour = INK, size = rel(0.95), angle = 90,
                                  margin = margin(r = 6)),
      axis.title.x = element_text(colour = INK, size = rel(0.95),
                                  margin = margin(t = 6)),
      axis.text = element_text(colour = INK, size = rel(0.85)),
      axis.line = element_line(colour = "#AEB9C4", linewidth = 0.5),
      axis.ticks = element_line(colour = "#AEB9C4", linewidth = 0.4),
      axis.ticks.length = unit(3, "pt"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E6EBF0", linewidth = 0.45),
      legend.title = element_text(colour = INK, size = rel(0.85), face = "bold"),
      legend.text = element_text(colour = INK, size = rel(0.85)),
      legend.key.width = unit(18, "pt"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(8, 14, 8, 8),
      strip.text = element_text(face = "bold", colour = NAVY, size = rel(0.95),
                                margin = margin(b = 6))
    )
}
theme_set(theme_tox())
update_geom_defaults("text", list(family = BASE_FAMILY, colour = INK))
update_geom_defaults("label", list(family = BASE_FAMILY))

# Figures are drawn on a canvas 0.8x the nominal size, so that when placed in
# their slide box the text renders ~1.25x larger (>= ~12 pt on the slide).
FIG_SCALE <- 0.8
save_fig <- function(p, name, w, h, scale = FIG_SCALE) {
  out <- file.path("figures", name)
  ggsave(out, p, width = w * scale, height = h * scale, units = "in", dpi = 300,
         device = ragg::agg_png, bg = "white")
  message("wrote ", out)
  invisible(out)
}

write_json_pretty <- function(x, name) {
  out <- file.path("data", name)
  writeLines(toJSON(x, auto_unbox = TRUE, pretty = TRUE, digits = NA), out)
  message("wrote ", out)
  invisible(out)
}

# numpy.gradient equivalent (central differences, one-sided ends)
np_gradient <- function(y, x) {
  n <- length(y)
  g <- numeric(n)
  g[1] <- (y[2] - y[1]) / (x[2] - x[1])
  g[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])
  if (n > 2) {
    i <- 2:(n - 1)
    g[i] <- (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
  }
  g
}

trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
