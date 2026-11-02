# Study 1 -- the QSP library as a safety atlas.
#
# Walks every disease directory of the repository, counts ODE compartments
# (unique dxdt_ assignments) in each mrgsolve model, nodes/clusters in each
# Graphviz map and PubMed references in each reference file, classifies each
# model into a therapeutic area from the README gallery table, and flags a
# model as carrying a "safety axis" when its model code, mechanistic map,
# directory README and references together contain
# >= 3 DISTINCT terms of a toxicology lexicon (word-bounded; acronyms such as
# QT, AKI, PML, DILI, irAE matched case-sensitively so that e.g. "making" does
# not count as AKI).
#
# Output: data/library_models.csv, data/library_summary.json,
#         figures/s1_areas.png, s1_quality.png, s1_atlas.png
# (No time-ordered growth figure, by design.)

source("R/theme_tox.R")

ROOT <- normalizePath(file.path("..", ".."))
readme <- readLines(file.path(ROOT, "README.md"), encoding = "UTF-8", warn = FALSE)

# ---- README gallery table: number, category, directory ---------------------
rows <- grep("^\\|\\s*[0-9]+\\s*\\|", readme, value = TRUE, perl = TRUE)
mm <- regmatches(rows, regexec(
  "^\\|\\s*([0-9]+)\\s*\\|\\s*([^|]+?)\\s*\\|.*?\\]\\(([a-z0-9\\-]+)/\\)",
  rows, perl = TRUE))
mm <- mm[lengths(mm) == 4]
gallery <- tibble(num = as.integer(sapply(mm, `[`, 2)),
                  category_ko = trimws(sapply(mm, `[`, 3)),
                  dir = sapply(mm, `[`, 4))
N_GALLERY <- nrow(gallery)
message("README gallery rows: ", N_GALLERY)

KO2EN <- c(
  "종양" = "Oncology", "신경" = "Neurology", "내분비·대사" = "Endocrine/metabolic",
  "소화기·간담도" = "GI/hepatobiliary", "혈액" = "Hematology",
  "심혈관" = "Cardiovascular", "신장·비뇨" = "Renal/urologic",
  "희귀·유전" = "Rare/genetic", "감염" = "Infectious", "호흡기" = "Respiratory",
  "자가면역·류마티스" = "Autoimmune/rheumatic", "정신·신경" = "Psychiatric",
  "피부" = "Dermatologic", "안과" = "Ophthalmic", "혈관염" = "Vasculitis",
  "근골격·신경근육" = "Musculoskeletal/NM", "부인·생식" = "Gynecologic/repro",
  "신생아·소아" = "Neonatal/pediatric", "중독·환경" = "Toxicology/environ.",
  "기타" = "Other", "면역결핍" = "Immunodeficiency", "알레르기" = "Allergy")

# ---- safety / toxicology lexicon ------------------------------------------
SAFETY_WORDS <- c("toxic", "adverse", "safety", "hepatotox", "nephrotox",
  "cardiotox", "torsade", "myelosuppress", "neutropeni", "cytokine release",
  "immune-related", "bleeding", "hemorrhag", "haemorrhag", "overdose",
  "poisoning", "drug-induced", "ototox", "neurotox", "bone loss", "hypoglyc",
  "lactic acidosis", "serotonin syndrome", "infection risk", "anaphyla",
  "teratogen", "withdrawal", "black box", "liver injury", "kidney injury",
  "thrombocytopeni", "cytopeni", "agranulocyt", "side effect", "side-effect",
  "dose-limiting", "hypersensitiv", "mucositis", "osteonecrosis", "tolerability")
SAFETY_ACR <- c("QT", "QTc", "AKI", "DILI", "PML", "irAE", "CRS", "SJS", "TEN",
  "TEAE", "DLT", "MTD", "CIPN")

safety_terms <- function(txt) {
  w <- vapply(SAFETY_WORDS, function(k)
    grepl(paste0("\\b", k), txt, ignore.case = TRUE, perl = TRUE), logical(1))
  a <- vapply(SAFETY_ACR, function(k)
    grepl(paste0("\\b", k, "\\b"), txt, perl = TRUE), logical(1))
  c(names(w)[w], names(a)[a])
}

read_txt <- function(f) {
  if (!length(f) || !file.exists(f[1])) return("")
  paste(readLines(f[1], encoding = "UTF-8", warn = FALSE), collapse = "\n")
}

count_odes <- function(rtext) {
  hits <- regmatches(rtext, gregexpr("dxdt_([A-Za-z0-9_]+)\\s*=", rtext, perl = TRUE))[[1]]
  n <- length(unique(sub("\\s*=$", "", hits)))
  if (n > 0) return(n)
  m <- regmatches(rtext, regexec("\\$CMT\\b(.*?)(\\$[A-Z]|$)", rtext, perl = TRUE))[[1]]
  if (length(m) >= 2) {
    return(length(regmatches(m[2], gregexpr("[A-Za-z][A-Za-z0-9_]*", m[2]))[[1]]))
  }
  0L
}

map_stats <- function(dot) {
  clusters <- lengths(regmatches(dot, gregexpr("subgraph\\s+cluster", dot, perl = TRUE)))
  nodes <- 0L
  for (s in trimws(strsplit(dot, "\n", fixed = TRUE)[[1]])) {
    if (grepl("->", s, fixed = TRUE) || grepl("^(//|#)", s)) next
    if (grepl("[", s, fixed = TRUE) && grepl("]", s, fixed = TRUE)) {
      after <- sub("^[^\\[]*\\[", "", s)
      if (!grepl("=", after, fixed = TRUE)) next
      head <- trimws(sub("\\[.*$", "", s))
      if (nzchar(head) &&
          !grepl("^(graph|node|edge|subgraph|rankdir|label|\\{|\\})", head)) {
        nodes <- nodes + 1L
      }
    }
  }
  c(nodes = nodes, clusters = clusters)
}

skip <- c("scripts", "docs", "presentation", "talks")
dirs <- list.dirs(ROOT, recursive = FALSE, full.names = TRUE)
dirs <- dirs[!grepl("^\\.", basename(dirs)) & !(basename(dirs) %in% skip)]

models <- list()
for (d in sort(dirs)) {
  rf <- c(list.files(d, "_mrgsolve_model\\.R$", full.names = TRUE),
          list.files(d, "mrgsolve.*\\.R$", full.names = TRUE))
  if (!length(rf)) rf <- setdiff(list.files(d, "_model\\.R$", full.names = TRUE),
                               list.files(d, "shiny", full.names = TRUE))
  if (!length(rf)) next
  rtext <- read_txt(rf)
  dotf <- list.files(d, "\\.dot$", full.names = TRUE)
  ms <- if (length(dotf)) map_stats(read_txt(dotf)) else c(nodes = 0L, clusters = 0L)
  reff <- c(list.files(d, "_references\\.md$", full.names = TRUE),
            list.files(d, "reference.*\\.md$", full.names = TRUE))
  reftxt <- read_txt(reff)
  refs <- lengths(regmatches(reftxt, gregexpr("PMID[:\\s]*\\d+", reftxt,
                                              perl = TRUE, ignore.case = TRUE)))
  if (refs == 0) refs <- lengths(regmatches(reftxt, gregexpr("https?://pubmed", reftxt)))
  # the safety axis may live in the code, the map, the README or the references
  dottxt <- if (length(dotf)) read_txt(dotf) else ""
  readmetxt <- read_txt(file.path(d, "README.md"))
  terms <- safety_terms(paste(rtext, dottxt, readmetxt, reftxt, sep = "\n"))
  models[[length(models) + 1]] <- tibble(
    dir = basename(d), odes = count_odes(rtext),
    nodes = ms[["nodes"]], clusters = ms[["clusters"]], refs = refs,
    safety_terms = length(terms), safety = length(terms) >= 3)
}
models <- bind_rows(models) |>
  left_join(gallery |> select(dir, category_ko), by = "dir") |>
  mutate(category_ko = coalesce(category_ko, "기타"),
         area = unname(KO2EN[category_ko]),
         area = coalesce(area, "Other"))
message("models mined: ", nrow(models))
readr::write_csv(models |> select(-category_ko), "data/library_models.csv")

by_area <- models |> count(area, name = "total") |>
  left_join(models |> filter(safety) |> count(area, name = "safety"), by = "area") |>
  mutate(safety = coalesce(safety, 0L)) |> arrange(desc(total))

summ <- list(
  n_gallery = N_GALLERY,
  n_mined = nrow(models),
  n_areas = n_distinct(gallery$category_ko),
  n_safety = sum(models$safety),
  pct_safety = round(100 * mean(models$safety)),
  total_odes = sum(models$odes),
  total_refs = sum(models$refs),
  median_odes = median(models$odes),
  min_odes = min(models$odes[models$odes > 0]),
  max_odes = max(models$odes),
  median_nodes = median(models$nodes),
  median_refs = median(models$refs),
  n_odes_gate = sum(models$odes >= 15),
  n_refs_gate = sum(models$refs >= 30),
  pct_refs_gate = round(100 * mean(models$refs >= 30))
)
write_json_pretty(summ, "library_summary.json")
print(str(summ))

# ============================================================================
# FIG 1 -- therapeutic-area coverage with safety-axis fraction
# ============================================================================
# areas with < 14 models are pooled so the bar labels stay legible
small <- by_area |> filter(total < 14)
ba <- bind_rows(by_area |> filter(total >= 14),
                tibble(area = sprintf("Other %d areas", nrow(small)),
                       total = sum(small$total), safety = sum(small$safety))) |>
  mutate(area = factor(area, levels = rev(area)))
p <- ggplot(ba, aes(y = area)) +
  geom_col(aes(x = total, fill = "Models in area"), width = 0.72) +
  geom_col(aes(x = safety, fill = "Carry an explicit safety / toxicity axis"), width = 0.72) +
  geom_text(aes(x = total, label = total), hjust = -0.35, size = 3.6, colour = INK) +
  scale_fill_manual(NULL, values = c("Models in area" = LGREY,
                                     "Carry an explicit safety / toxicity axis" = TEAL),
                    breaks = c("Carry an explicit safety / toxicity axis", "Models in area")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Number of disease models", y = NULL) +
  theme(panel.grid.major.y = element_blank(), axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top", legend.justification = "left")
save_fig(p, "s1_areas.png", 8.4, 4.9)

# ============================================================================
# FIG 2 -- quality-gate distributions (ODE states, references)
# ============================================================================
gate_hist <- function(x, binw, gate, fill, xlab, ttl, xmax) {
  df <- tibble(x = pmin(x, xmax))
  ggplot(df, aes(x)) +
    geom_histogram(binwidth = binw, boundary = gate, fill = fill, colour = "white", linewidth = 0.3) +
    geom_vline(xintercept = gate, colour = RED, linetype = "dashed", linewidth = 0.9) +
    annotate("text", x = Inf, y = Inf, label = paste0("Gate ≥ ", gate, " (dashed)"),
             hjust = 1.05, vjust = 1.6, colour = RED, fontface = "bold", size = 3.9) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = xlab, y = "Models", title = ttl)
}
p1 <- gate_hist(models$odes, 4, 15, BLUE, "ODE compartments per model",
                sprintf("Model size (median %d)", as.integer(summ$median_odes)), 90)
p2 <- gate_hist(models$refs, 12, 30, TEAL, "PubMed references per model",
                sprintf("Literature grounding (median %d)", as.integer(summ$median_refs)), 240)
save_fig(p1 + p2, "s1_quality.png", 10.4, 4.2)

# ============================================================================
# FIG 3 -- the atlas: ODE size x map size, marker area ~ references
# ============================================================================
top_areas <- head(by_area$area, 7)
atl <- models |>
  mutate(group = ifelse(area %in% top_areas, area, "Other areas"),
         group = factor(group, levels = c(top_areas, "Other areas"))) |>
  arrange(group == "Other areas")
p <- ggplot(atl |> arrange(desc(group == "Other areas")),
            aes(odes, nodes, colour = group, size = refs)) +
  geom_point(alpha = 0.62, stroke = 0) +
  scale_colour_manual(NULL, values = setNames(c(CYCLE[1:7], LGREY), c(top_areas, "Other areas"))) +
  scale_size_area(max_size = 9, guide = "none") +
  scale_x_continuous(limits = c(10, NA)) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), nrow = 2)) +
  labs(x = "ODE compartments (PK + disease dynamics)", y = "Mechanistic-map nodes",
       caption = "Each point is one disease model; marker area ~ reference count") +
  theme(legend.position = "bottom",
        plot.caption = element_text(colour = GREY, size = 11, hjust = 1))
save_fig(p, "s1_atlas.png", 10.2, 5.2)

message("S1 done")
