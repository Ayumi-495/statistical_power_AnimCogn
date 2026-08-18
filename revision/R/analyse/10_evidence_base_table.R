# revision/R/analyse/10_evidence_base_table.R -------------------------------------------
# Step 10: the characteristics of the 28 included meta-analytical papers.
#
# This table does not currently exist, and four reviewer comments depend on it:
# R1C17 (describe the taxonomic and demographic composition), R2C11 (characteristics
# of the included meta-analyses and primary studies), R2C15 (study-level information;
# none of the contributing studies are described) and R2C17 (readers cannot tell what
# portion of the field the results represent). The revised Results already cite it,
# and the response to reviewers already promises it.
#
# Source: the systematic map that supplied the corpus,
#   .../systematic_mapping_AnimCogn/data/S1mapping_Nov2024.xlsx
# keyed by `meta_id` (MA01 ... MA49), joined to our 48 analysis models through the
# `source_paper` column of revision/results/per_meta_analysis_estimates.csv.
#
# The map is read-only and belongs to another project. Nothing here writes to it.

source(here::here("revision", "R", "analyse", "00_revision_functions.R"))
suppressMessages(library(readxl))

message("== 10: evidence-base characteristics table ==")

MAP <- file.path("/Users/ayumi/Library/CloudStorage/GoogleDrive-amizuno@ualberta.ca",
                 "My Drive/systematic_mapping_AnimCogn/data/S1mapping_Nov2024.xlsx")
if (!file.exists(MAP))
  stop("The systematic-map dataset is not available at:\n  ", MAP,
       "\nThis script needs it; it belongs to the sibling project systematic_mapping_AnimCogn.")

sheet <- function(s) suppressMessages(readxl::read_excel(MAP, sheet = s))

p <- readr::read_csv(file.path(REV_OUT, "per_meta_analysis_estimates.csv"),
                     show_col_types = FALSE)
stopifnot(nrow(p) == 48L)

# --- what our 48 models contribute, per source paper --------------------------
ours <- p |>
  dplyr::group_by(meta_id = source_paper) |>
  dplyr::summarise(
    n_models_included   = dplyr::n(),
    effect_size_metrics = paste(sort(unique(effect_size_type)), collapse = ", "),
    k_effect_sizes      = sum(k),
    n_study_clusters    = sum(n_study_cluster),
    n_sign_reversals_reported = sum(reversal_c3),
    .groups = "drop")
message(sprintf("our corpus: %d source papers, %d models, %d effect sizes",
        nrow(ours), sum(ours$n_models_included), sum(ours$k_effect_sizes)))

# --- characteristics from the systematic map ----------------------------------
info <- sheet("summary_note") |>
  dplyr::select(meta_id, authors_id, doi, year) |>
  dplyr::left_join(dplyr::select(sheet("taxon_focus"), meta_id, taxonomic_focus), by = "meta_id") |>
  dplyr::left_join(dplyr::select(sheet("new_categorisation"), meta_id, journal,
                                 contents1, contents2, contents3,
                                 factors_1, factors_2, factors_3, factors_4,
                                 life_stage_reported, life_stage_included,
                                 sex_reported, sex_included), by = "meta_id") |>
  dplyr::left_join(dplyr::select(sheet("task"), meta_id, paradigm, task), by = "meta_id") |>
  dplyr::left_join(dplyr::select(sheet("manipulation"), meta_id, manipulation,
                                 manipulation_category), by = "meta_id") |>
  dplyr::left_join(dplyr::select(sheet("effect size_numbers"), meta_id,
                                 effect_size_statistics,
                                 map_number_studies      = number_studies,
                                 map_number_effect_sizes = number_effect_sizes,
                                 map_number_species      = number_species), by = "meta_id") |>
  dplyr::left_join(dplyr::select(sheet("species_focus"), meta_id, research_aim, species_focus),
                   by = "meta_id") |>
  # `cognition_category` holds more than one row per paper (a paper can span domains),
  # so collapse to one row per paper rather than joining and duplicating.
  dplyr::left_join(
    sheet("cognition_category") |>
      dplyr::group_by(meta_id) |>
      dplyr::summarise(cognition_domain = paste(sort(unique(cognition)), collapse = "; "),
                       cognition_subtopic = paste(sort(unique(stats::na.omit(sub_topic))), collapse = "; "),
                       .groups = "drop"),
    by = "meta_id")

paste_nonmissing <- function(...) {
  v <- c(...); v <- v[!is.na(v) & nzchar(trimws(v)) & !tolower(trimws(v)) %in% c("na", "none")]
  if (!length(v)) NA_character_ else paste(unique(trimws(v)), collapse = "; ")
}

tab <- ours |>
  dplyr::left_join(info, by = "meta_id") |>
  dplyr::rowwise() |>
  dplyr::mutate(
    cognitive_domain = paste_nonmissing(contents1, contents2, contents3),
    moderators_examined = paste_nonmissing(factors_1, factors_2, factors_3, factors_4)
  ) |>
  dplyr::ungroup() |>
  dplyr::select(
    meta_id, authors_id, year, journal, doi,
    taxonomic_focus, species_focus, research_aim, cognition_domain, cognition_subtopic,
    cognitive_domain_study_type = cognitive_domain, paradigm, task,
    manipulation, manipulation_category, moderators_examined,
    life_stage_reported, life_stage_included, sex_reported, sex_included,
    effect_size_statistics,
    map_number_studies, map_number_effect_sizes, map_number_species,
    n_models_included, effect_size_metrics, k_effect_sizes, n_study_clusters,
    n_sign_reversals_reported
  ) |>
  dplyr::arrange(meta_id)

# --- verification gates -------------------------------------------------------
miss <- tab$meta_id[is.na(tab$authors_id)]
if (length(miss))
  stop("No systematic-map record for: ", paste(miss, collapse = ", "))
stopifnot(nrow(tab) == 28L, sum(tab$k_effect_sizes) == 5740L,
          sum(tab$n_models_included) == 48L)

write_revision(tab, "evidence_base_characteristics.csv")

# --- the composition claims the Results paragraph makes, recomputed ------------
# The revised Results assert specific counts. They are checked here rather than
# transcribed, because the paragraph was written before this table existed.
rodent <- grepl("Rattus norvegicus", tab$species_focus) & grepl("Mus musculus", tab$species_focus) &
          !grepl("Bos taurus|Gallus|Macaca|Canis", tab$species_focus)
dom <- function(w) sum(grepl(w, tab$cognition_domain, ignore.case = TRUE))
message(sprintf("\nRESULTS-PARAGRAPH CHECKS (manuscript claim -> computed):"))
message(sprintf("  'Sixteen papers were restricted to rats and/or mice'  -> %d papers match a rodent taxonomic focus", sum(rodent)))
message(sprintf("  'learning ... 16 papers'                             -> %d", dom("learning")))
message(sprintf("  'memory ... 13 papers'                               -> %d", dom("memory")))
message(sprintf("  'Twelve syntheses examined ... manipulations'         -> %d papers record a manipulation",
        sum(!is.na(tab$manipulation) & !tolower(trimws(tab$manipulation)) %in% c("na", "none", ""))))
message(sprintf("  life stage reported in %d of 28; included as a moderator in %d",
        sum(tolower(trimws(tab$life_stage_reported)) %in% c("yes", "y", "true")),
        sum(tolower(trimws(tab$life_stage_included)) %in% c("yes", "y", "true"))))
message(sprintf("  sex reported in %d of 28; included as a moderator in %d",
        sum(tolower(trimws(tab$sex_reported)) %in% c("yes", "y", "true")),
        sum(tolower(trimws(tab$sex_included)) %in% c("yes", "y", "true"))))

# --- primary-study sample size: what is actually recoverable ------------------
# Reviewer 2 (R2C11) asked for the average sample size of the primary studies. It is
# NOT recoverable for the whole corpus. For Fisher's Zr it follows exactly from the
# sampling variance, n = 3 + 1/v. For SMD it can be approximated only under an
# equal-group-size assumption. For lnRR the sampling variance depends on the means and
# standard deviations as well as the sample sizes, so it cannot be inverted.
S <- readRDS(file.path(REV_TMP, "original_estimates.rds"))
L <- all_datasets(S$dat); o <- S$original
zr <- which(o$effect_size_type == "Zr")
n_zr <- unlist(lapply(zr, function(i) 3 + 1 / L[[i]]$var))
message(sprintf("\nR2C11 sample size: Zr exactly recoverable for %d effect sizes; median n = %.0f (IQR %.0f-%.0f)",
        length(n_zr), stats::median(n_zr), stats::quantile(n_zr, .25), stats::quantile(n_zr, .75)))
message(sprintf("  SMD (%d effect sizes) approximable only under equal group sizes; lnRR (%d) not recoverable from the variance",
        sum(o$k[o$effect_size_type == "SMD"]), sum(o$k[o$effect_size_type == "lnRR"])))
message("  -> report the Zr subset exactly, state the SMD assumption if used, and say plainly that lnRR cannot be reconstructed.")
