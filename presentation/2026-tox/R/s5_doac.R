# Study 5 -- Aging and organ function: a PBPK/PD workflow for direct oral
# anticoagulants (DOACs) in the patients Phase 1 excludes (elderly, renal
# impairment).
#
# age -> creatinine clearance -> exposure (Rowland renal-fraction scaling)
#     -> bleeding / stroke (exposure-response) -> fixed vs individualized dose
#
# Renal-clearance fractions: dabigatran 0.80, edoxaban 0.50, rivaroxaban 0.33,
# apixaban 0.27.
# Validation anchors (AUC fold vs normal renal function; US labels / PK studies;
# mild CrCl 50-80, moderate 30-50, severe <30):
#   dabigatran  mild 1.5,  moderate 3.2   (severe not formally studied)
#   edoxaban    mild 1.32, moderate 1.74, severe 1.72
#   rivaroxaban mild 1.44, moderate 1.52, severe 1.64
#   apixaban    mild 1.16, moderate 1.29, severe 1.44
# Mechanistic result: renal scaling under-predicts dabigatran (P-gp /
# absorption); a mechanism-based correction -- not a free fit -- fixes it.

source("R/theme_tox.R")
set.seed(7)

DOAC <- tribble(
  ~drug,         ~fr,  ~col,
  "dabigatran",  0.80, RED,
  "edoxaban",    0.50, GOLD,
  "rivaroxaban", 0.33, BLUE,
  "apixaban",    0.27, TEAL)
OBS <- tribble(
  ~drug,         ~cat,       ~obs,
  "dabigatran",  "mild",     1.5,
  "dabigatran",  "moderate", 3.2,
  "edoxaban",    "mild",     1.32,
  "edoxaban",    "moderate", 1.74,
  "edoxaban",    "severe",   1.72,
  "rivaroxaban", "mild",     1.44,
  "rivaroxaban", "moderate", 1.52,
  "rivaroxaban", "severe",   1.64,
  "apixaban",    "mild",     1.16,
  "apixaban",    "moderate", 1.29,
  "apixaban",    "severe",   1.44)
CAT_CRCL <- c(mild = 65, moderate = 40, severe = 20)

auc_fold <- function(drug, crcl, corrected = FALSE) {
  fr <- DOAC$fr[match(drug, DOAC$drug)]
  RF <- crcl / 100
  fold <- 1 / (1 - fr * (1 - RF))
  corr <- corrected & (drug == "dabigatran")   # mechanism-based P-gp / absorption term
  fold * (1 + corr * 0.9 * (1 - RF))
}
crcl_of_age <- function(a) pmin(pmax(95 - 1.175 * (a - 40), 18), 125)

# ============================================================================
# FIG A -- age -> renal function and the dose-reduction threshold
# ============================================================================
ages <- 20:95
cloud <- tibble(age = rep(seq(20, 95, by = 2), each = 40)) |>
  mutate(m = crcl_of_age(age), crcl = rnorm(n(), m, 0.18 * m), age = age + rnorm(n(), 0, 0.5))
cross <- ages[which(crcl_of_age(ages) <= 50)[1]]
pA <- ggplot() +
  geom_point(data = cloud, aes(age, crcl), colour = "#C3CDD7", size = 0.9, alpha = 0.7) +
  geom_hline(yintercept = 50, colour = GOLD, linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 30, colour = RED, linetype = "dashed", linewidth = 0.8) +
  annotate("label", x = 20, y = 55, label = "Dose-reduction threshold (CrCl 50)", hjust = 0,
           vjust = 0, colour = GOLD, fontface = "bold", size = 3.7, fill = alpha("white", 0.9),
           border.colour = NA, label.padding = unit(1.5, "pt")) +
  annotate("label", x = 20, y = 21, label = "Severe impairment (CrCl 30)", hjust = 0,
           vjust = 0, colour = RED, fontface = "bold", size = 3.7, fill = alpha("white", 0.9),
           border.colour = NA, label.padding = unit(1.5, "pt")) +
  geom_line(data = tibble(age = ages, crcl = crcl_of_age(ages)), aes(age, crcl,
            linetype = "Population mean CrCl"), colour = NAVY, linewidth = 1.4) +
  annotate("curve", x = 88, y = 92, xend = cross + 0.6, yend = 52, curvature = 0.25,
           colour = GOLD, linewidth = 0.8, arrow = arrow(length = unit(7, "pt"), type = "closed")) +
  annotate("text", x = 89, y = 96, label = sprintf("Mean crosses 50\nat age %d", cross),
           hjust = 0.5, vjust = 0, colour = GOLD, fontface = "bold", size = 3.9, lineheight = 1) +
  scale_linetype_manual(NULL, values = "solid") +
  coord_cartesian(ylim = c(15, 160)) +
  labs(x = "Age (years)", y = "Creatinine clearance (mL/min)") +
  theme(legend.position = c(0.84, 0.92))
save_fig(pA, "s5_age_crcl.png", 9.0, 4.6)

# ============================================================================
# FIG B -- validation: predicted vs observed AUC fold + dabigatran correction
# ============================================================================
val <- OBS |> mutate(crcl = CAT_CRCL[cat],
                     pred = auc_fold(drug, crcl),
                     pred_c = auc_fold(drug, crcl, corrected = TRUE),
                     drug = factor(drug, levels = DOAC$drug))
gmfe <- exp(mean(abs(log(val$pred / val$obs))))
gmfe_c <- exp(mean(abs(log(val$pred_c / val$obs))))
band <- tibble(x = c(1, 7)) |> mutate(lo = x / 2, hi = x * 2)
dab <- val |> filter(drug == "dabigatran")
brks <- c(1, 1.5, 2, 3, 5, 7)
pB <- ggplot() +
  geom_ribbon(data = band, aes(x, ymin = lo, ymax = hi), fill = GREEN, alpha = 0.08) +
  geom_abline(slope = 1, intercept = 0, colour = GREY) +
  geom_abline(slope = 1, intercept = log10(2), colour = GREY, linetype = "dotted") +
  geom_abline(slope = 1, intercept = -log10(2), colour = GREY, linetype = "dotted") +
  geom_segment(data = dab, aes(x = obs, xend = obs, y = pred * 1.04, yend = pred_c / 1.05),
               colour = RED, linewidth = 0.8, arrow = arrow(length = unit(7, "pt"), type = "closed")) +
  geom_point(data = dab, aes(obs, pred_c), shape = 21, colour = RED, fill = "white", stroke = 1.5, size = 4.5) +
  geom_point(data = val, aes(obs, pred, fill = drug), shape = 21, colour = "white", stroke = 0.8, size = 4.2) +
  annotate("text", x = 1.04, y = 6.7, hjust = 0, vjust = 1, colour = RED, size = 3.5, lineheight = 1.05,
           label = "Dabigatran:\nrenal scaling\nunder-predicts") +
  scale_fill_manual("Renal-fraction PBPK", values = setNames(DOAC$col, DOAC$drug),
                    labels = tools::toTitleCase(DOAC$drug)) +
  scale_x_log10(breaks = brks, labels = brks) + scale_y_log10(breaks = brks, labels = brks) +
  coord_equal(xlim = c(1, 7), ylim = c(1, 7), expand = FALSE) +
  labs(x = "Observed AUC fold-increase (label)", y = "Predicted AUC fold-increase (PBPK)") +
  theme(legend.position = "right")
save_fig(pB, "s5_validation.png", 7.6, 4.9)

# ============================================================================
# FIG C -- exposure-response: bleeding up, stroke down, age shifts both
# ============================================================================
bleed <- function(e, age) (0.9 + 0.025 * (age - 70)) * e^1.5
stroke <- function(e, age) (3.0 + 0.04 * (age - 70)) / (1 + 2.2 * e)
expo <- seq(0.3, 3.5, length.out = 200)
er <- bind_rows(lapply(c(55, 80), function(a) bind_rows(
  tibble(e = expo, rate = bleed(expo, a), what = "Major bleeding", age = paste("Age", a)),
  tibble(e = expo, rate = stroke(expo, a), what = "Stroke / SE", age = paste("Age", a)))))
pC <- ggplot(er, aes(e, rate, colour = what, linetype = age)) +
  annotate("rect", xmin = 0.8, xmax = 1.6, ymin = -Inf, ymax = Inf, fill = GREEN, alpha = 0.09) +
  annotate("text", x = 1.2, y = 7.6, label = "Target\nexposure", colour = GREEN, fontface = "bold",
           size = 3.6, lineheight = 1, vjust = 1) +
  geom_line(linewidth = 1.2) +
  scale_colour_manual(NULL, values = c("Major bleeding" = RED, "Stroke / SE" = TEAL)) +
  scale_linetype_manual(NULL, values = c("Age 55" = "solid", "Age 80" = "22")) +
  scale_y_continuous(limits = c(0, 8)) +
  labs(x = "DOAC exposure (relative to young-normal target)", y = "Events per 100 patient-years",
       title = "Two risks move in opposite directions") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "pt"))
save_fig(pC, "s5_exposure_response.png", 5.8, 4.8)

# ============================================================================
# FIG D -- dose individualization in a virtual elderly dabigatran cohort
# ============================================================================
N <- 4000
age_c <- runif(N, 75, 92)
crcl_c <- pmin(pmax(rnorm(N, crcl_of_age(age_c), 11), 15), 95)
expo_fixed <- auc_fold("dabigatran", crcl_c, corrected = TRUE)
reduce <- (age_c >= 80) | (crcl_c < 50)            # 150 -> 110 mg bid
expo_adj <- ifelse(reduce, expo_fixed * 0.73, expo_fixed)
CEIL <- auc_fold("dabigatran", 50, corrected = TRUE)
br_fixed <- mean(bleed(expo_fixed, age_c)); br_adj <- mean(bleed(expo_adj, age_c))
over_fixed <- 100 * mean(expo_fixed > CEIL); over_adj <- 100 * mean(expo_adj > CEIL)

coh <- bind_rows(tibble(e = expo_fixed, rule = "Fixed full dose (150 mg bid)"),
                 tibble(e = expo_adj, rule = "Adjusted per label rule (110 mg bid)"))
d1 <- ggplot(coh, aes(e, fill = rule)) +
  annotate("rect", xmin = 1, xmax = CEIL, ymin = -Inf, ymax = Inf, fill = GREEN, alpha = 0.08) +
  geom_histogram(binwidth = 4 / 29, boundary = 1, position = "identity", alpha = 0.6, colour = NA) +
  geom_vline(xintercept = CEIL, colour = NAVY, linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = CEIL + 0.06, y = Inf, label = "Ceiling", hjust = 0, vjust = 1.5,
           colour = NAVY, fontface = "bold", size = 3.6) +
  scale_fill_manual(NULL, values = c("Fixed full dose (150 mg bid)" = RED, "Adjusted per label rule (110 mg bid)" = TEAL),
                    breaks = c("Fixed full dose (150 mg bid)", "Adjusted per label rule (110 mg bid)")) +
  coord_cartesian(xlim = c(1, 5)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Predicted exposure (× young-normal)", y = "Virtual patients",
       title = "Exposure distribution") +
  theme(legend.position = "bottom", legend.direction = "vertical")
bars <- tibble(rule = factor(c("Fixed\ndose", "Adjusted\ndose"), levels = c("Fixed\ndose", "Adjusted\ndose")),
               pct = c(over_fixed, over_adj), col = c(RED, TEAL))
d2 <- ggplot(bars, aes(rule, pct, fill = I(col))) +
  geom_col(width = 0.58) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), vjust = -0.5, fontface = "bold", size = 4.6) +
  scale_y_continuous(limits = c(0, max(bars$pct) * 1.25), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Above ceiling (%)", title = "Over-exposure") +
  theme(panel.grid.major.x = element_blank())
save_fig((d1 | plot_spacer() | d2) + plot_layout(widths = c(1.25, 0.1, 1)), "s5_individualization.png", 6.8, 4.8)

res <- list(
  gmfe_renal_only = round(gmfe, 2), gmfe_corrected = round(gmfe_c, 2),
  age_cross_crcl50 = cross,
  overexposed_fixed_pct = round(over_fixed), overexposed_adj_pct = round(over_adj),
  bleed_fixed = round(br_fixed, 1), bleed_adj = round(br_adj, 1),
  dabi_moderate_obs = 3.2,
  dabi_moderate_pred_renal = round(auc_fold("dabigatran", 40), 2),
  dabi_moderate_pred_corr = round(auc_fold("dabigatran", 40, TRUE), 2),
  n_cohort = N)
write_json_pretty(res, "s5_results.json")
print(str(res))
message("S5 done")
