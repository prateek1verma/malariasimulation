#!/usr/bin/env Rscript

# Simple smoke test for scheduled mosquito genotype releases.
# Run from the repository root, e.g.:
#   Rscript test-scripts/check_release_simple.R

# devtools::load_all()

load_pkg <- function() {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
    return(invisible(TRUE))
  }
  if (requireNamespace("malariasimulationGD", quietly = TRUE)) {
    library(malariasimulationGD)
    return(invisible(TRUE))
  }
  stop("Need either pkgload (preferred) or an installed malariasimulationGD package")
}

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) {
    stop(msg, call. = FALSE)
  }
}

assert_equal <- function(x, y, msg) {
  if (!isTRUE(all.equal(x, y, check.attributes = FALSE))) {
    stop(
      paste0(
        msg, "\n  expected: ", paste(capture.output(str(y)), collapse = " "),
        "\n  actual: ", paste(capture.output(str(x)), collapse = " ")
      ),
      call. = FALSE
    )
  }
}

load_pkg()

if (!requireNamespace("MGDrivE", quietly = TRUE)) {
  stop("MGDrivE is required for this script")
}

release_day <- 20L
release_n <- 50L
timesteps <- 25L

set.seed(20260225)

params <- get_parameters(list(
  individual_mosquitoes = TRUE,
  human_population = 50,
  total_M = 200,
  init_foim = 0,
  progress_bar = FALSE
))
params <- parameterise_total_M(params, params$total_M)

cube_3 <- MGDrivE::cubeMendelian(gtype = c("AA", "Aa", "aa"))
cube_3$releaseType <- "aa"
params$cube <- cube_3

# Intentionally omit releaseGenotype: set_releases() should default to cube$releaseType.
params <- set_releases(params, list(
  releasesStart = release_day,
  releasesNumber = 1,
  releaseCount = release_n,
  releaseSex = "M"
))

out <- run_resumable_simulation(timesteps, parameters = params)
df <- out$data
geno <- out$mosquito_genotypes
schedule <- attr(df, "mosquito_release_schedule")

expected_schedule <- data.frame(
  timestep = release_day,
  species = "gamb",
  sex = "M",
  genotype = "aa",
  count = release_n,
  stringsAsFactors = FALSE
)

assert_true(!is.null(schedule), "Expected mosquito_release_schedule attribute")
assert_equal(schedule, expected_schedule, "Release schedule metadata does not match expected single event")

assert_true("n_released_gamb" %in% names(df), "Expected n_released_gamb column in output")
assert_equal(as.integer(df$n_released_gamb[release_day]), release_n, "Release-day n_released_gamb mismatch")
assert_equal(as.integer(sum(df$n_released_gamb[-release_day])), 0L, "Expected no releases on non-release days")

# Male aa count should jump exactly on the release day in this single-release setup.
assert_equal(as.integer(unname(geno$male[release_day, "aa"])), release_n, "Male aa count on release day mismatch")

# On the release day, female non-WT counts should still be zero (male-only release).
female_non_wt_release_day <- geno$female[release_day, c("Aa", "aa"), drop = FALSE]
assert_true(all(female_non_wt_release_day == 0), "Female non-WT counts should be 0 on the male release day")

cat("PASS: simple scheduled male genotype release behaves as expected.\n")
cat(sprintf("Release day %d, n_released_gamb=%d, male aa=%d\n",
  release_day,
  as.integer(df$n_released_gamb[release_day]),
  as.integer(unname(geno$male[release_day, "aa"]))
))

