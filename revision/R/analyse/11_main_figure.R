# revision/R/analyse/11_main_figure.R ---------------------------------------------------
# Step 11: the replacement for manuscript Figure 3.
#
# WHAT CHANGED AND WHY, so this is not undone by accident.
#
# 1. THE TILE GRID IS GONE. The 48-row colour grid plotted exactly the same numbers as
#    the points on the violins above it: the submitted `scatter_plot.R` passes one data
#    frame to both `geom_tile` and the violin layer, so 192 values were drawn twice.
#    Reviewer 1 asked whether the lower panels were "just another way to display the raw
#    data included in the violin plots above". They were. Removing them is the answer to
#    that comment, and it frees the space used below.
#
# 2. THE PRIMARY-STUDY VIOLIN NOW SHOWS ALL 5,740 ESTIMATES, not 48 per-model medians.
#    The primary-study level is the paper's main result, and its distribution did not
#    appear anywhere in the submitted figure - only a summary of a summary did.
#
# 3. THE SUMMARY LINE IS THE NUMBER THE PAPER REPORTS. The submitted figure drew
#    `stat_summary(fun = mean)`, an unweighted arithmetic mean of the plotted values.
#    At the meta-analysis level that is 0.710 for uncorrected power, while the text
#    reports 0.822. A reader could not find the paper's own figure on its own figure.
#    Here the line is read from the canonical summary tables.
#
# 4. NO CENSORING AND NO WHITE-ON-WHITE. The submitted fill scale ran to white, so the
#    19 models at power >= 0.99 were invisible, and `limits = c(0, 20)` with no `oob`
#    turned 11 Type M tiles grey, reading as missing data. Neither failure mode can
#    occur here: identity is carried by position, colour and fill, and the script
#    asserts that nothing falls outside the plotted range.
#
# 5. A THIRD SERIES. The bias-robust estimate (Yang et al. 2024) is now a reported
#    sensitivity analysis, so it appears beside the uncorrected and Yang-2023 values.
#
# Colours are Okabe-Ito, assigned in fixed order and validated for colour-vision
# deficiency and contrast (worst adjacent pair deutan dE 11.0; all three above 3:1
# against the surface). Fill duplicates position, so identity is never colour-alone.

source(here::here("revision", "R", "analyse", "00_revision_functions.R"))
suppressMessages(library(ggplot2))

message("== 11: main figure (violins, tile grid removed) ==")

S  <- readRDS(file.path(REV_TMP, "original_estimates.rds"))
o  <- S$original; L <- all_datasets(S$dat)
BR <- readRDS(file.path(REV_TMP, "bias_robust.rds"))
stopifnot(identical(o$MA_model, BR$MA_model))

FIG_DIR <- file.path(REV_OUT, "figures"); dir.create(FIG_DIR, showWarnings = FALSE)

EFFECTS <- tibble::tribble(
  ~key,          ~label,                ~mu,                ~ma_se,
  "uncorrected", "Uncorrected",         o$beta0,            o$se_beta0,
  "yang2023",    "Yang 2023 corrected", o$beta0_c3,         o$se_beta0,
  "yang2024",    "Yang 2024 bias-robust", BR$FE_VCV_estimate, BR$FE_VCV_CRVE_SE_CR2
)
PAL <- setNames(c("#0072B2", "#D55E00", "#009E73"), EFFECTS$label)
# Points use a darker shade of the same hue than the violin fill. With 5,740 values
# crushed into the low end of the power axis, points drawn in the fill colour simply
# dissolve into it; the contrast is what makes the density readable.
PAL_PT <- setNames(c("#00436B", "#8F3E00", "#00614A"), EFFECTS$label)

# --- the plotted values -------------------------------------------------------
# meta-analysis level: one value per model (48). primary-study level: one value per
# effect-size estimate (5,740), mu constant within a meta-analysis.
sei_all <- unlist(lapply(L, function(x) x$sei), use.names = FALSE)
idx     <- rep(seq_along(L), vapply(L, nrow, integer(1)))
stopifnot(length(sei_all) == 5740L)

mk <- function(mu, se, lab, level) {
  tibble::tibble(level = level, effect = lab,
                 power = power_two_tailed_cf(mu, se),
                 type_M = type_M_cf(mu, se),
                 type_S = type_S_cf(mu, se))
}
d <- purrr::list_rbind(lapply(seq_len(nrow(EFFECTS)), function(i) {
  e <- EFFECTS[i, ]
  dplyr::bind_rows(
    mk(e$mu[[1]][idx], sei_all,     e$label, "Primary-study level"),
    mk(e$mu[[1]],      e$ma_se[[1]], e$label, "Meta-analysis level"))
})) |>
  tidyr::pivot_longer(c(power, type_M, type_S), names_to = "metric", values_to = "value") |>
  dplyr::mutate(
    level  = factor(level, levels = c("Primary-study level", "Meta-analysis level")),
    effect = factor(effect, levels = EFFECTS$label),
    metric = factor(metric, levels = c("power", "type_M", "type_S")))

# --- the summary lines, read from the canonical tables ------------------------
sm <- dplyr::bind_rows(
  readr::read_csv(file.path(REV_OUT, "primary_level_sensitivity.csv"), show_col_types = FALSE),
  readr::read_csv(file.path(REV_OUT, "meta_analysis_level_sensitivity.csv"), show_col_types = FALSE)) |>
  dplyr::filter(role %in% c("reference_uncorrected", "primary", "reported_sensitivity"),
                # the reported estimand at each level: the study-cluster random-intercept
                # summary at the primary-study level, the effect-size-count-weighted
                # geometric mean at the meta-analysis level. The two model-level
                # estimands are sensitivity analyses and are not drawn as the bar.
                weighting %in% c("study_cluster_random_intercept", "k_effect_sizes")) |>
  dplyr::transmute(
    level = factor(ifelse(aggregation == "primary_study_level",
                          "Primary-study level", "Meta-analysis level"),
                   levels = levels(d$level)),
    effect = factor(dplyr::recode(role, reference_uncorrected = "Uncorrected",
                                  primary = "Yang 2023 corrected",
                                  reported_sensitivity = "Yang 2024 bias-robust"),
                    levels = EFFECTS$label),
    metric = factor(metric, levels = levels(d$metric)),
    value = geometric_mean, raw_median, raw_q1, raw_q3) |>
  # TYPE S IS SUMMARISED BY THE MEDIAN AND QUARTILES AT BOTH LEVELS. The log-scale model
  # requires an offset of 0.025 because Type S can be exactly zero, and the fitted
  # summary is sensitive to that offset wherever the observed values are small relative
  # to it - which is every meta-analysis-level cell and the uncorrected primary-study
  # cell. The manuscript therefore reads Type S off the raw median and interquartile
  # range at both levels, and the bar on the figure has to be the number in the sentence.
  #
  # Until 2026-08-17 this line read `& level == "Meta-analysis level"`, so the
  # primary-study bars were the back-transformed model estimates (2.76 / 10.21 / 4.84%)
  # while the text reported the medians (1.83 / 13.58 / 5.24%). Power and Type M are
  # unaffected: they need no offset and keep the model-based summary.
  #
  # Nothing is recomputed here. `raw_median`, `raw_q1` and `raw_q3` are read from the
  # canonical summary tables, the same columns Table S1 is built from.
  dplyr::mutate(
    use_raw = metric == "type_S",
    value = ifelse(use_raw, raw_median, value),
    lo    = ifelse(use_raw, raw_q1, NA_real_),
    hi    = ifelse(use_raw, raw_q3, NA_real_))
stopifnot(nrow(sm) == 18L)

# --- panels -------------------------------------------------------------------
panel <- function(mt, title, yscale, ylab) {
  dd <- dplyr::filter(d, metric == mt); ss <- dplyr::filter(sm, metric == mt)
  ggplot(dd, aes(effect, value, fill = effect, colour = effect)) +
    geom_violin(width = 0.85, linewidth = 0.25, alpha = 0.22, colour = NA, scale = "width") +
    # Points are shown at BOTH levels. The two layers differ only in size and opacity,
    # because one holds 5,740 values and the other 48. Showing the primary-study points
    # matters: the violin's kernel density smooths across the hard bounds, so the real
    # accumulation of estimates at power = 1 and at Type S = 0.5 is otherwise invisible.
    geom_jitter(data = dplyr::filter(dd, level == "Primary-study level"),
                width = 0.36, height = 0, size = 0.22, alpha = 0.16, stroke = 0) +
    geom_jitter(data = dplyr::filter(dd, level == "Meta-analysis level"),
                width = 0.13, height = 0, size = 0.85, alpha = 0.8, stroke = 0) +
    geom_linerange(data = dplyr::filter(ss, !is.na(lo)),
                   aes(x = effect, ymin = lo, ymax = hi),
                   linewidth = 0.35, colour = "grey15", inherit.aes = FALSE) +
    geom_crossbar(data = ss, aes(x = effect, y = value, ymin = value, ymax = value),
                  width = 0.7, linewidth = 0.45, colour = "grey15",
                  inherit.aes = FALSE) +
    facet_wrap(~ level, nrow = 1) +
    scale_fill_manual(values = PAL, name = NULL) +
    scale_x_discrete(labels = c("Uncorrected", "Yang 2023\ncorrected", "Yang 2024\nbias-robust")) +
    scale_colour_manual(values = PAL_PT, name = NULL) +
    yscale +
    labs(x = NULL, y = ylab, title = title) +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linewidth = 0.15, colour = "grey92"),
          panel.border = element_rect(linewidth = 0.3, colour = "grey70"),
          strip.background = element_rect(fill = "grey96", colour = NA),
          strip.text = element_text(size = 7.2),
          axis.text.x = element_text(size = 6.6, colour = "grey20"),
          axis.ticks.x = element_line(linewidth = 0.2),
          plot.title = element_text(size = 8.5, face = "bold"),
          legend.position = "none")
}

g1 <- panel("power", "a. Statistical power",
            scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)), "Probability")
g2 <- panel("type_M", "b. Type M error (exaggeration ratio)",
            scale_y_log10(breaks = c(1, 3, 10, 100, 1000, 10000),
                          labels = c("1", "3", "10", "100", "1000", "10000")), "Ratio (log scale)")
g3 <- panel("type_S", "c. Type S error",
            scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.125)), "Probability")

g <- patchwork::wrap_plots(g1, g2, g3, ncol = 1)
# NO CAPTION IS DRAWN INTO THE IMAGE. An earlier version baked six lines of caption
# text in with `plot_annotation(caption = ...)`. Journals typeset figure captions
# separately from the artwork, so that text would have appeared twice in the
# typeset article. The caption now lives only in the manuscript, and the wording
# there is the authority. Anything the reader needs in order to decode the panels
# must therefore be in the manuscript caption, not here: which summary each
# horizontal bar shows, that the violins are scaled to equal width rather than
# equal area, and that the meta-analysis-level Type S bar is a raw median with an
# interquartile range while every other bar is model-based.
#
# No legend anywhere either: the three groups are named directly on the x axis, so a
# colour key would be redundant and, applied per panel, would repeat itself three times.

# --- checks that would have caught the submitted figure's defects -------------
oob <- sum(dplyr::filter(d, metric == "power")$value > 1 | dplyr::filter(d, metric == "power")$value < 0) +
       sum(dplyr::filter(d, metric == "type_S")$value > 0.5) +
       sum(dplyr::filter(d, metric == "type_M")$value < 1)
if (oob > 0) stop(sprintf("%d values fall outside the plotted range and would be dropped", oob))
message(sprintf("scale check: 0 of %d plotted values fall outside their axis limits", nrow(d)))
message("summary bars (read from the canonical tables, not recomputed):")
for (i in seq_len(nrow(sm))) if (sm$metric[i] == "power")
  message(sprintf("   %-20s %-22s %.4f", sm$level[i], sm$effect[i], sm$value[i]))

ggsave(file.path(FIG_DIR, "main_metrics.pdf"), g, width = 6.8, height = 8.6, device = "pdf")
ggsave(file.path(FIG_DIR, "main_metrics.png"), g, width = 6.8, height = 8.6, dpi = 200)
message("wrote figures/main_metrics.{pdf,png}")
