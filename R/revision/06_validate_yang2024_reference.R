# R/revision/06_validate_yang2024_reference.R -----------------------------------
# EXTERNAL VALIDATION of our Yang et al. (2024) implementation against the authors'
# own published worked example.
#
# Why this script exists. The FE + VCV point estimate has two independent derivations
# in this repository (a closed-form `V^-1 1` weighted mean and a `metafor::rma.mv`
# fit), but the CR2 standard error has only one implementation - clubSandwich's. A
# hand-rolled second CR2 would not be an independent derivation; it would be a
# re-implementation tuned to agree with the first, which is precisely the class of
# self-invented diagnostic that every error in this project's audit came from.
# So the second check is an EXTERNAL ANCHOR instead: run our code path on the
# authors' own example data and require it to reproduce their published output.
#
# NOT part of run_all.R, because it needs a 160 kB data file from another project's
# GitHub repository and run_all.R must work offline. Run it on its own:
#   Rscript R/revision/06_validate_yang2024_reference.R
# Its verified output is committed as results/revision/yang2024_reference_validation.csv
# so the result is durable without re-downloading.

source(here::here("R", "revision", "00_revision_functions.R"))

message("== 06: external validation against Yang et al. (2024) worked example ==")

# --- the authors' example data ------------------------------------------------
# Bird et al. (2019) Ecology Letters, redistributed in the tutorial repository named
# by the paper's data-availability statement.
DATA_URL <- paste0("https://raw.githubusercontent.com/Yefeng0920/",
                   "BiasRobustMA_tutorial/main/data/bird.et.al.2019.ecoletts.csv")
DATA_MD5 <- "1bcad390a96bdc6a8e07a81a2a31347e"
local <- file.path(REV_TMP, "bird.et.al.2019.ecoletts.csv")

if (!file.exists(local)) {
  message("downloading the tutorial's example data ...")
  ok <- tryCatch({ utils::download.file(DATA_URL, local, quiet = TRUE); TRUE },
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok || !file.exists(local))
    stop("Could not download the tutorial data. This script needs network access.\n",
         "  Source: ", DATA_URL)
}
got_md5 <- unname(tools::md5sum(local))
if (!identical(got_md5, DATA_MD5))
  stop("The example data does not match the verified file.\n",
       "  expected md5 ", DATA_MD5, "\n  got      md5 ", got_md5)
message("example data md5 verified: ", DATA_MD5)

# --- the tutorial's own data preparation, in its order ------------------------
# hands_on_R.Rmd: read.csv -> select(study, eff.size, var.eff.size) -> filter
# outliers to (-20, 20) -> obs <- 1:nrow(dat). The order matters: `obs` is assigned
# after the outlier filter, so it is a contiguous 1..k index on the filtered data.
raw <- utils::read.csv(local)
dat <- raw[, c("study", "eff.size", "var.eff.size")]
dat <- dat[dat$eff.size < 20 & dat$eff.size > -20, ]
dat$obs <- seq_len(nrow(dat))
message(sprintf("k = %d effect sizes, J = %d studies", nrow(dat),
                dplyr::n_distinct(dat$study)))

# --- our implementation, unmodified ------------------------------------------
f <- fit_fe_vcv_cr2(dat$eff.size, dat$var.eff.size, dat$study, rho = 0.5)
stopifnot(!is.null(f))

# the working (step one) model's own CI, for the step-one comparison
crit_work <- stats::qt(0.975, nrow(dat) - 1L)

# --- the authors' published values --------------------------------------------
# Taken from the RENDERED tutorial (BiasRobustMA_tutorial/R/hands_on_R.html), which
# contains the values their own run produced, quoted to 3 significant decimals:
#   step one: beta = 0.074, SE = 0.018, 95% CI = [0.039, 0.108], df = 1632
#   step two: SE = 0.053, p = 0.168, 95% CI = [-0.032, 0.18], df 1632 -> 52
# The paper's prose quotes beta = 0.075, SE = 0.054, t = 1.375, p = 0.175,
# CI = [-0.034, 0.184] for the same example; those come from the 448-model corpus
# pipeline with its own exclusions, agree to two decimals, and are recorded below as
# a second, non-binding reference rather than as the target.
checks <- tibble::tribble(
  ~step,      ~quantity,        ~ours,                                  ~published, ~tol,
  "step_one", "estimate",       f$beta,                                 0.074,      5e-4,
  "step_one", "se_working",     f$se_working,                           0.018,      5e-4,
  "step_one", "ci_lower",       f$beta - crit_work * f$se_working,      0.039,      5e-4,
  "step_one", "ci_upper",       f$beta + crit_work * f$se_working,      0.108,      5e-4,
  "step_two", "se_CR2",         f$se_cr2,                               0.053,      5e-4,
  "step_two", "df_Satterthwaite", round(f$df_satt),                     52,         0.5,
  "step_two", "p_value",        f$pval,                                 0.168,      5e-4,
  "step_two", "ci_lower",       f$ci_lb,                               -0.032,      5e-4,
  "step_two", "ci_upper",       f$ci_ub,                                0.180,      5e-4
) |>
  dplyr::mutate(abs_diff = abs(ours - published), matches = abs_diff < tol)

print(as.data.frame(checks), digits = 6)

message(sprintf("\ninternal cross-checks: vcalc vs build_vcv max|diff| = %.3g | ",
                f$vcv_max_abs_diff),
        sprintf("metafor::robust vs clubSandwich direct max|diff| = %.3g",
                f$wrapper_max_abs_diff))

if (!all(checks$matches))
  stop("VALIDATION FAILED against the published worked example - do not use the ",
       "CR2 results until this is resolved.")
if (f$vcv_max_abs_diff > 1e-10 || f$wrapper_max_abs_diff > 1e-10)
  stop("Internal cross-checks exceeded tolerance.")

message("\nVALIDATED: all ", nrow(checks),
        " published values reproduced at the precision they are reported to.")

out <- checks |>
  dplyr::mutate(
    source = "Yefeng0920/BiasRobustMA_tutorial R/hands_on_R.html (rendered tutorial)",
    example_dataset = "bird.et.al.2019.ecoletts.csv",
    example_data_md5 = DATA_MD5,
    specification = "FE+VCV (rho=0.5) then CR2 with Satterthwaite df, cluster = study",
    k = nrow(dat), n_study = dplyr::n_distinct(dat$study),
    vcalc_vs_build_vcv_max_abs_diff = f$vcv_max_abs_diff,
    metafor_vs_clubsandwich_max_abs_diff = f$wrapper_max_abs_diff,
    R_version = PKG_VERSIONS[["R"]], metafor_version = PKG_VERSIONS[["metafor"]],
    clubSandwich_version = PKG_VERSIONS[["clubSandwich"]]
  )
write_revision(out, "yang2024_reference_validation.csv")
