# 00_setup.R -----------------------------------------------------------------
# Biology Open revision. Shared setup: packages, functions, data loading.
#
# S2_v2.R is the immutable baseline and is never sourced, edited, or run by
# these scripts. Everything here re-implements its documented behaviour so the
# revision analyses can be run with the audit corrections applied.
#
# Corrections applied here relative to S2_v2.R, each documented at the line it affects:
#   - the Zr directory is read as "zr", not "Zr". S2_v2.R:98,103 use "Zr",
#     which resolves only on a case-insensitive filesystem; on Linux it silently
#     returns character(0) and the pipeline yields 37 models instead of 48.
#   - hard assertions on the data hierarchy (48 models, 28 papers, 5,740 rows)
#     so a silent load failure cannot pass as a valid run.
#   - a fixed seed, because error_M() is a Monte Carlo estimator and neither
#     S2_v2.R nor Yang et al. set one.

pacman::p_load(tidyverse, metafor, lme4, here)

SEED <- 20260810L
set.seed(SEED)

# expected hierarchy, verified against the data during the audit
EXPECTED_N_MODELS <- 48L
EXPECTED_N_PAPERS <- 28L
EXPECTED_N_ROWS   <- 5740L

# metric-level metrics ------------------------------------------------------
# Verbatim from S2_v2.R:15-44 (itself verbatim from Yang et al. Rmd:46-84),
# except that error_M() now takes an explicit alpha/N and is called under a
# seeded RNG. Not modified: the estimators themselves must stay identical for
# the reproduction check in 01_estimates.R to be meaningful.

power_two_tailed <- function(mu, se, alpha = 0.05) {
  2 - pnorm(qnorm(1 - alpha / 2) - abs(mu) / se) -
      pnorm(qnorm(1 - alpha / 2) + abs(mu) / se)
}

error_S <- function(mu, se, alpha = 0.05) {
  p.u <- 1 - pnorm(qnorm(1 - alpha / 2) - abs(mu) / se)
  p.l <- pnorm(-qnorm(1 - alpha / 2) - abs(mu) / se)
  p.l / (p.u + p.l)
}

error_M <- function(mu, se, alpha = 0.05, N = 10000) {
  est.random <- rnorm(n = N, mean = mu, sd = se)
  sig.index <- abs(est.random) > se * qnorm(1 - alpha / 2)
  overestimate <- mean(abs(est.random)[sig.index]) / abs(mu)
  round(abs(overestimate), 3)
}

# vectorised over se with a single mu, which is how both levels use it
error_M_vec <- function(mu, se, alpha = 0.05, N = 10000) {
  vapply(se, function(s) error_M(mu, s, alpha = alpha, N = N), numeric(1))
}

# data loading ---------------------------------------------------------------
# Mirrors S2_v2.R:94-228. The five SMD/des_stat datasets carry raw descriptive
# statistics rather than es/var, so SMD is computed with escalc() for those and
# gains the effective-sample-size moderators. Reviewer 2 asked (§3.3) how
# "retained in original metrics" squares with "calculated from raw descriptive
# statistics": these five datasets are the answer -- 43 of 48 supply es/var
# directly, 5 are computed. No metric is ever converted into another.

load_datasets <- function(root = here::here()) {
  read_dir <- function(path) {
    files <- list.files(file.path(root, path), pattern = "[.]csv$", full.names = TRUE)
    out <- lapply(files, readr::read_csv, show_col_types = FALSE, progress = FALSE)
    names(out) <- basename(files)
    out
  }

  lnRR <- read_dir("lnRR")
  smd_direct <- read_dir("SMD")
  zr <- read_dir("zr")          # audit fix: lower case, see header
  smd_raw <- read_dir(file.path("SMD", "des_stat"))

  # effective sample size for the raw-descriptive SMD datasets (S2_v2.R:120-147)
  smd_raw <- lapply(smd_raw, function(d) {
    d$ess <- (4 * d$C_n * d$T_n) / (d$C_n + d$T_n)
    d$ess.var <- 1 / d$C_n + 1 / d$T_n
    d$ess.sei <- sqrt(d$ess.var)
    d
  })
  smd_computed <- lapply(smd_raw, function(d) {
    e <- metafor::escalc(measure = "SMD", m1i = T_mean, m2i = C_mean,
                         sd1i = T_sd, sd2i = C_sd, n1i = T_n, n2i = C_n, data = d)
    names(e)[names(e) == "yi"] <- "es"
    names(e)[names(e) == "vi"] <- "var"
    e
  })

  # S2_v2.R:170 -- append() order matters downstream: the ess-moderator models
  # are matched to the last five SMD datasets by position.
  smd <- append(smd_direct, smd_computed)

  prep <- function(L) {
    lapply(L, function(d) {
      d <- d[!is.na(d$es) & !is.na(d$var) & d$var != 0 & !is.na(d$year), ]
      d$sei <- sqrt(d$var)
      d$year_pub.l <- as.vector(d$year - max(d$year))
      d
    })
  }

  list(
    lnRR = prep(lnRR),
    SMD  = prep(smd),
    Zr   = prep(zr),
    smd_uses_ess = names(smd) %in% names(smd_computed)
  )
}

# One flat table of dataset-level metadata, in the canonical order used by every
# downstream script: lnRR, then SMD (direct then computed), then Zr.
dataset_index <- function(dat) {
  tibble::tibble(
    case = c(names(dat$lnRR), names(dat$SMD), names(dat$Zr)),
    es_type = c(rep("lnRR", length(dat$lnRR)),
                rep("SMD",  length(dat$SMD)),
                rep("Zr",   length(dat$Zr))),
    uses_ess = c(rep(FALSE, length(dat$lnRR)), dat$smd_uses_ess,
                 rep(FALSE, length(dat$Zr)))
  ) |>
    dplyr::mutate(
      # source meta-analytical paper: the MA__ filename prefix. Verified during
      # the audit to give exactly 28 papers across the 48 models.
      paper = stringr::str_extract(case, "^[A-Za-z]+[0-9]+")
    )
}

all_datasets <- function(dat) c(dat$lnRR, dat$SMD, dat$Zr)

# Assert the hierarchy rather than trusting the load. A silent path failure is
# the specific defect this guards against.
check_hierarchy <- function(dat, idx) {
  L <- all_datasets(dat)
  n_models <- length(L)
  n_rows <- sum(vapply(L, nrow, integer(1)))
  n_papers <- dplyr::n_distinct(idx$paper)

  stopifnot(
    "expected 48 meta-analysis models"  = n_models == EXPECTED_N_MODELS,
    "expected 28 meta-analytical papers" = n_papers == EXPECTED_N_PAPERS,
    "expected 5,740 effect-size rows"    = n_rows == EXPECTED_N_ROWS,
    "dataset index out of step with data" = nrow(idx) == n_models
  )
  message(sprintf("hierarchy OK: %d models, %d papers, %d effect-size rows",
                  n_models, n_papers, n_rows))
  invisible(TRUE)
}

# optimizer table -----------------------------------------------------------
# S2_v2.R uses "optim" everywhere except two places: the intercept-only Zr fits
# (S2_v2.R:311-319, with a comment that "optim" failed there) and the SMD
# bias-detection model fitted with the `sei` moderator (S2_v2.R:350-358). The
# five SMD datasets whose detection model uses `ess.sei` instead
# (S2_v2.R:360-368) use "optim", as do all scenario-specific correction models.
# Reproducing the baseline requires mirroring this exactly, so it is encoded
# rather than harmonised -- getting it wrong shifts scenario assignment for at
# least one dataset and moves the corrected results.
optimizer_for <- function(es_type, stage, uses_ess = FALSE) {
  if (stage == "null"   && es_type == "Zr")  return("nlminb")
  if (stage == "detect" && es_type == "SMD" && !uses_ess) return("nlminb")
  "optim"
}

fit_rma <- function(data, mods = NULL, optimizer = "optim") {
  args <- list(
    yi = quote(es), V = quote(var),
    random = list(~ 1 | study_ID / obs_ID),
    method = "REML", test = "t", data = data, sparse = TRUE,
    control = list(optimizer = optimizer)
  )
  if (!is.null(mods)) args$mods <- mods
  tryCatch(do.call(metafor::rma.mv, args), error = function(e) NULL)
}

out_dir <- function(...) {
  p <- here::here("revision", "results", "reproduce", ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}
