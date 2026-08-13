# R/revision/08_model_level_figure.R --------------------------------------------
# Step 8: the per-model values behind the meta-analysis-level summaries, as a table
# and as a figure.
#
# PURPOSE. The meta-analysis-level summary is a k-weighted geometric mean over 48
# models, and it is not representative of them: one model can move it by 45%. Showing
# all 48 values lets a reader see the spread and the leverage directly instead of
# taking the aggregate on trust.
#
# WRITES TO results/revision/figures/. The submitted `figure/` directory is frozen and
# is not touched.
#
# DESIGN NOTES, so a later editor does not undo them by accident:
#   - Models are ordered by k, descending, because k IS the weight. Leverage therefore
#     reads top to bottom, and each label carries its own k.
#   - Three panels, one per metric, with free x scales. NOT a dual axis.
#   - Type M is on a log10 scale (values span 1 to 4,260). Power and Type S are linear
#     on their natural bounds. Type S is NOT logged: it reaches 1e-90, which would
#     give ~90 useless decades.
#   - Colours are Okabe-Ito #0072B2 / #D55E00 / #009E73, assigned in fixed order and
#     validated: all adjacent pairs clear the CVD separation threshold (worst 11.0
#     deutan) and all three clear 3:1 contrast against the surface. Shape duplicates
#     colour, so identity is never carried by colour alone.
#   - The dashed vertical rules are the k-weighted summaries. They are drawn so the
#     gap between the summary and the models it summarises is visible.
#   - UWLS is omitted from the figure only because it is the supplementary
#     specification; its per-model values are in the CSV written below.
#
# The predecessor figure in the submitted manuscript has two defects this one avoids:
# a `scale_fill_gradient(limits = c(0, 20))` with no `oob`, which censors 11 tiles to
# NA so they render grey and read as missing data, and a ramp running to white, which
# hides the 19 models with power >= 0.99. Do not reintroduce either.

source(here::here("R", "revision", "00_revision_functions.R"))
suppressMessages(library(ggplot2))

message("== 08: model-level values, table and figure ==")
p <- readr::read_csv(file.path(REV_OUT, "per_meta_analysis_estimates.csv"),
                     show_col_types = FALSE)
stopifnot(nrow(p) == 48L)

FIG_DIR <- file.path(REV_OUT, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# --- the 48 values, long ------------------------------------------------------
specs <- tibble::tribble(
  ~effect_estimator,          ~label,                    ~mu,                ~se,                          ~se_source,  ~se_method,          ~role,
  "uncorrected_beta0",        "Uncorrected",             p$beta0,            p$se_beta0,                   "se_beta0",  NA_character_,       "reference_uncorrected",
  "yang2023_gated_beta0_c3",  "Yang 2023 (primary)",     p$beta0_c3,         p$se_beta0,                   "se_beta0",  NA_character_,       "primary",
  "yang2024_FE_VCV",          "Yang 2024 FE + VCV",      p$FE_VCV_estimate,  p$FE_VCV_CRVE_SE_CR2,         "own_CRVE",  "CR2_Satterthwaite", "reported_sensitivity",
  "yang2024_UWLS",            "Yang 2024 UWLS",          p$UWLS_estimate,    p$UWLS_CRVE_SE_CR2_naive_t,   "own_CRVE",  "CR2_naive_t",       "supplementary"
)

long <- purrr::list_rbind(lapply(seq_len(nrow(specs)), function(i) {
  s <- specs[i, ]; mu <- s$mu[[1]]; se <- s$se[[1]]
  tibble::tibble(
    MA_model = p$MA_model, effect_size_type = p$effect_size_type, k = p$k,
    n_study_cluster = p$n_study_cluster,
    effect_estimator = s$effect_estimator, estimator_label = s$label,
    se_source = s$se_source, se_method = s$se_method, role = s$role,
    crit_value_method = "z_1.96",
    assumed_effect = mu, standard_error = se,
    power  = power_two_tailed_cf(mu, se),
    type_M = type_M_cf(mu, se),
    type_S = type_S_cf(mu, se)
  )
}))
stopifnot(nrow(long) == 48L * 4L)
write_revision(long, "model_level_metrics.csv")

# --- figure -------------------------------------------------------------------
PAL <- c("Uncorrected" = "#0072B2", "Yang 2023 (primary)" = "#D55E00",
         "Yang 2024 FE + VCV" = "#009E73")
SHP <- c("Uncorrected" = 16, "Yang 2023 (primary)" = 17, "Yang 2024 FE + VCV" = 15)

ord  <- p$MA_model[order(p$k, decreasing = TRUE)]
labs <- setNames(sprintf("%s  (k = %s)", sub("[.]csv", "", p$MA_model),
                         format(p$k, big.mark = ",", trim = TRUE)), p$MA_model)

d <- long |>
  dplyr::filter(estimator_label %in% names(PAL)) |>
  dplyr::mutate(MA_model = factor(MA_model, levels = rev(ord)),
                estimator_label = factor(estimator_label, levels = names(PAL)))

# k-weighted summaries, drawn as reference rules so the gap between the summary and
# the models it summarises is visible rather than asserted.
S <- readr::read_csv(file.path(REV_OUT, "meta_analysis_level_sensitivity.csv"),
                     show_col_types = FALSE)
rule <- S |>
  dplyr::filter(weighting == "k_effect_sizes",
                role %in% c("reference_uncorrected", "primary", "reported_sensitivity")) |>
  dplyr::transmute(
    estimator_label = factor(dplyr::recode(role,
        reference_uncorrected = "Uncorrected", primary = "Yang 2023 (primary)",
        reported_sensitivity  = "Yang 2024 FE + VCV"), levels = names(PAL)),
    metric = metric, value = geometric_mean)

hi <- data.frame(y = which(levels(d$MA_model) == "MA09.csv"))

# One panel per metric, built separately so each can carry its own scale. Three
# independent x scales on three panels - never two scales on one panel.
panel <- function(col, title, xscale, show_y) {
  ggplot(d, aes(.data[[col]], MA_model, colour = estimator_label, shape = estimator_label)) +
    geom_hline(data = hi, aes(yintercept = y), colour = "grey55",
               linewidth = 3.4, alpha = 0.20, inherit.aes = FALSE) +
    geom_vline(data = dplyr::filter(rule, metric == col),
               aes(xintercept = value, colour = estimator_label),
               linetype = "dashed", linewidth = 0.4, alpha = 0.9,
               inherit.aes = FALSE, show.legend = FALSE) +
    geom_point(size = 1.9, stroke = 0, alpha = 0.95) +
    scale_colour_manual(values = PAL, name = NULL, drop = FALSE) +
    scale_shape_manual(values = SHP, name = NULL, drop = FALSE) +
    scale_y_discrete(labels = if (show_y) labs else NULL) +
    xscale +
    labs(x = NULL, y = NULL, title = title) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.15, colour = "grey92"),
      panel.border = element_rect(linewidth = 0.3, colour = "grey70"),
      axis.text.y = if (show_y) element_text(size = 5.4, family = "mono") else element_blank(),
      axis.ticks.y = if (show_y) element_line(linewidth = 0.2) else element_blank(),
      axis.ticks.x = element_line(linewidth = 0.2),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      legend.position = "bottom", legend.key.height = grid::unit(9, "pt")
    )
}

g1 <- panel("power", "Statistical power",
            scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)), TRUE)
g2 <- panel("type_M", "Type M (exaggeration ratio, log scale)",
            scale_x_log10(breaks = c(1, 3, 10, 100, 1000, 4000),
                          labels = c("1", "3", "10", "100", "1000", "4000")), FALSE)
g3 <- panel("type_S", "Type S (probability of wrong sign)",
            scale_x_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.125)), FALSE)

g <- patchwork::wrap_plots(g1, g2, g3, nrow = 1, widths = c(1.35, 1, 1), guides = "collect") +
  patchwork::plot_annotation(
    title = "Design-analysis metrics for each of the 48 meta-analysis models",
    subtitle = paste0(
      "Models are ordered by k, the number of effect-size estimates, which is also the weight used in the meta-analysis-level summary; ",
      "leverage therefore reads top to bottom.\nDashed rules mark those k-weighted summaries. MA09 (shaded) alone holds 22.6% of the weight, ",
      "and dropping it moves the Yang-2024 summary from 0.479 to 0.693.\nAll metrics use a common normal-theory two-sided 5% threshold (z = 1.96). ",
      "Type S is not logged because it reaches 1e-90; Type M is, because it reaches 4,260."),
    theme = theme(
      plot.title = element_text(face = "bold", size = 10.5),
      plot.subtitle = element_text(size = 6.5, colour = "grey30", lineheight = 1.3))
  ) & theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "model_level_metrics.pdf"), g, width = 10, height = 8.4, device = "pdf")
ggsave(file.path(FIG_DIR, "model_level_metrics.png"), g, width = 10, height = 8.4, dpi = 200)
message("wrote figures/model_level_metrics.{pdf,png}")

# Sanity checks that would have caught the defects in the submitted Figure 2: nothing
# may be dropped by a scale limit, and nothing may be invisible.
oob <- sum(d$power < 0 | d$power > 1) + sum(d$type_S < 0 | d$type_S > 0.5) + sum(d$type_M < 1)
if (oob > 0) stop(sprintf("%d points fall outside the plotted scales and would be silently dropped", oob))
message(sprintf("scale check: 0 of %d plotted points fall outside their axis limits", nrow(d) * 3))
