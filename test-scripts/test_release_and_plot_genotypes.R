#!/usr/bin/env Rscript

# test_release_and_plot_genotypes.R
# Run from repo root:
#   Rscript test-scripts/test_release_and_plot_genotypes.R

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(malariasimulationGD)
  }
  library(MGDrivE)
})

# ----------------------------
# Helpers
# ----------------------------
stopf <- function(...) stop(sprintf(...), call. = FALSE)

assert_true <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
}

# Compare data.frame *column contents* bitwise (ignores attributes/rownames differences)
assert_df_columns_identical <- function(a, b, labelA = "A", labelB = "B") {
  if (!identical(names(a), names(b))) {
    stopf("FAIL: column names/order differ between %s and %s", labelA, labelB)
  }
  if (nrow(a) != nrow(b)) stopf("FAIL: nrow differs between %s and %s", labelA, labelB)
  
  diffs <- names(a)[!vapply(names(a), function(nm) identical(a[[nm]], b[[nm]]), logical(1))]
  if (length(diffs) > 0) {
    stopf("FAIL: dataframe column content differs (%s vs %s). First differing column: %s",
          labelA, labelB, diffs[1])
  }
  invisible(TRUE)
}

get_geno_from_out <- function(out_obj) {
  # Works for both run_simulation (attributes on df) and run_resumable_simulation (list output)
  if (is.list(out_obj) && !is.null(out_obj$mosquito_genotypes)) {
    return(out_obj$mosquito_genotypes)
  }
  if (is.data.frame(out_obj)) {
    F <- attr(out_obj, "mosquito_genotype_counts_female")
    M <- attr(out_obj, "mosquito_genotype_counts_male")
    V <- attr(out_obj, "mosquito_genotype_V")
    Tad <- attr(out_obj, "mosquito_genotype_total_adults")
    if (is.null(F) || is.null(M)) return(NULL)
    return(list(female = F, male = M, V = V, total_adults = Tad))
  }
  NULL
}

# ----------------------------
# Simulation setup
# ----------------------------
year <- 365
timesteps <- 200L

seed <- 123

params_base <- get_parameters(list(
  individual_mosquitoes = TRUE,
  model_seasonality = FALSE,
  progress_bar = FALSE,
  human_population = 200,
  total_M = 2000,
  init_foim = 0
))
params_base <- parameterise_total_M(params_base, params_base$total_M)

# ----------------------------
# 1) Baseline vs WT-only cube exact check
# ----------------------------
cat("== (1) Baseline vs WT-only cube equality check ==\n")

cube_WT <- MGDrivE::cubeMendelian(gtype = c("AA"))

set.seed(seed)
out0 <- run_simulation(timesteps, params_base)

params_wt <- params_base
params_wt$cube <- cube_WT

set.seed(seed)
out_wt <- run_simulation(timesteps, params_wt)

assert_df_columns_identical(out0, out_wt, "baseline(no cube)", "WT cube(AA only)")
cat("PASS: baseline(no cube) and WT cube(AA) have bitwise-identical dataframe columns under same seed.\n\n")

# ----------------------------
# 2) Mendelian 3-genotype cube + repeated aa male releases
# ----------------------------
cat("== (2) 3-genotype cube + repeated aa male releases ==\n")

cube_3 <- MGDrivE::cubeMendelian(gtype = c("AA", "Aa", "aa"))
cube_3$releaseType <- "aa"

params_3 <- params_base
params_3$cube <- cube_3

# Release schedule: start day 50, every 7 days, 10 releases of 200 males each
releases <- list(
  releasesStart = 50L,
  releasesNumber = 10L,
  releasesInterval = 7L,
  releaseCount = 200L,
  releaseSex = "M"
)

params_3 <- set_releases(params_3, releases)

set.seed(seed)
out_rel <- run_resumable_simulation(timesteps, parameters = params_3)
df <- out_rel$data
geno <- out_rel$mosquito_genotypes

assert_true(!is.null(attr(df, "mosquito_release_schedule")), "Expected mosquito_release_schedule attribute")
assert_true("n_released_gamb" %in% names(df), "Expected n_released_gamb column in df output")
cat("PASS: releases metadata + n_released_gamb present.\n")

# quick sanity: sum releases equals expected total
expected_total_release <- releases$releaseCount * releases$releasesNumber
actual_total_release <- sum(df$n_released_gamb, na.rm = TRUE)
assert_true(actual_total_release == expected_total_release,
            sprintf("Total released mismatch: expected %d, got %d",
                    expected_total_release, actual_total_release))
cat(sprintf("PASS: total released = %d.\n\n", actual_total_release))

# ----------------------------
# 3) Plot female adults by genotype
# ----------------------------
cat("== (3) Plotting female adult genotypes over time ==\n")

F <- geno$female
assert_true(!is.null(F), "Missing geno$female in output")
assert_true(all(colnames(F) %in% c("AA","Aa","aa")), "Unexpected genotype columns in female matrix")

t <- df$timestep
# If timestep column is missing in resumable output (unlikely), fallback
if (is.null(t)) t <- seq_len(nrow(F))

female_total <- rowSums(F)
cat(sprintf("Final day female totals: %s\n", paste(sprintf("%s=%d", colnames(F), F[nrow(F),]), collapse = ", ")))
cat(sprintf("Final day female_total=%d\n\n", female_total[nrow(F)]))

# Base R plot (no extra deps)
png("genotype_female_adults.png", width = 1200, height = 650)
matplot(
  x = t,
  y = F,
  type = "l",
  lty = 1,
  lwd = 2,
  xlab = "Timestep (day)",
  ylab = "Adult female count",
  main = "Adult female mosquitoes by genotype (Mendelian cube, repeated aa male releases)"
)
legend("topright", legend = colnames(F), lty = 1, lwd = 2, bty = "n")

# Overlay release days as vertical dashed lines
rel_days <- attr(df, "mosquito_release_schedule")$timestep
abline(v = rel_days, lty = 3)
dev.off()

cat("Wrote plot: genotype_female_adults.png\n\n")

# Also write a small CSV for convenience
geno_df <- data.frame(timestep = t, F, n_released_gamb = df$n_released_gamb)
write.csv(geno_df, "genotype_female_adults_timeseries.csv", row.names = FALSE)
cat("Wrote CSV: genotype_female_adults_timeseries.csv\n")

# --- Save + plot MALE genotype time series (like females) ---

M <- geno$male
stopifnot(!is.null(M))
t <- df$timestep
if (is.null(t)) t <- seq_len(nrow(M))

# Save CSV
male_df <- data.frame(timestep = t, M, n_released_gamb = df$n_released_gamb)
write.csv(male_df, "genotype_male_adults_timeseries.csv", row.names = FALSE)
cat("Wrote CSV: genotype_male_adults_timeseries.csv\n")

# Plot (base R)
png("genotype_male_adults.png", width = 1200, height = 650)
matplot(
  x = t,
  y = M,
  type = "l",
  lty = 1,
  lwd = 2,
  xlab = "Timestep (day)",
  ylab = "Adult male count (pool)",
  main = "Adult male mosquitoes by genotype (Mendelian cube, repeated aa male releases)"
)
legend("topright", legend = colnames(M), lty = 1, lwd = 2, bty = "n")

# Optional: mark release days
rel_days <- attr(df, "mosquito_release_schedule")$timestep
abline(v = rel_days, lty = 3)
dev.off()

cat("Wrote plot: genotype_male_adults.png\n")

cat("\nALL DONE.\n")

