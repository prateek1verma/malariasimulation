#!/usr/bin/env Rscript
############################################################
## Regression tests for new "queue male + queue female emergence,
## commit once per timestep before rendering" implementation.
##
## Tests:
##  (A) No-cube baseline vs WT-only cube: shared columns identical (exact).
##  (B) Emergence timing (WT-cube run): no large early male-only spike;
##      male_total and female_total should be close for early timesteps.
##  (C) Releases still work: n_released_gamb totals match schedule and
##      released genotype appears in male counts on release day.
##
## Outputs (CSV):
##  - out_baseline_no_cube.csv
##  - out_wt_only_cube.csv
##  - out_diff_wt_minus_baseline_shared_numeric.csv
##  - adult_totals_wt_cube.csv
##  - out_wt_cube_with_releases.csv
##  - geno_female_wt_release.csv
##  - geno_male_wt_release.csv
############################################################

suppressPackageStartupMessages({
  library(malariasimulationGD)
  library(MGDrivE)
  library(dplyr)
})

# -------------------------
# Helpers
# -------------------------
stopf <- function(...) stop(sprintf(...), call.=FALSE)
assert_true <- function(x, msg) if (!isTRUE(x)) stop(msg, call.=FALSE)

detect_time_col <- function(df) {
  candidates <- c("timestep", "t", "time", "day")
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) stopf("Could not detect time column. Columns are: %s", paste(names(df), collapse=", "))
  hit[1]
}

# Compare only shared columns, exactly (bitwise), ignoring attributes
assert_shared_columns_identical <- function(a, b, labelA="A", labelB="B") {
  shared <- intersect(names(a), names(b))
  if (length(shared) == 0) stopf("No shared columns to compare (%s vs %s).", labelA, labelB)
  if (nrow(a) != nrow(b)) stopf("FAIL: nrow differs (%s vs %s)", labelA, labelB)
  
  diffs <- shared[!vapply(shared, function(nm) identical(a[[nm]], b[[nm]]), logical(1))]
  if (length(diffs) > 0) {
    stopf("FAIL: first differing shared column (%s vs %s): %s", labelA, labelB, diffs[1])
  }
  invisible(TRUE)
}

mean_tail <- function(x, n = 60) mean(tail(x, n), na.rm = TRUE)

write_csv_safe <- function(df, path) {
  utils::write.csv(df, path, row.names = FALSE)
  cat("Wrote:", path, "\n")
}

# Extract totals with your actual column naming:
# - females (no-cube): Sm_gamb_count + Pm_gamb_count + Im_gamb_count
# - males (cube): from genotype matrices
extract_adult_totals <- function(df, mosquito_genotypes = NULL) {
  time_col <- detect_time_col(df)
  
  # FEMALES
  if (!is.null(mosquito_genotypes) && !is.null(mosquito_genotypes$female)) {
    female_total <- rowSums(mosquito_genotypes$female)
  } else if (all(c("Sm_gamb_count","Pm_gamb_count","Im_gamb_count") %in% names(df))) {
    female_total <- df$Sm_gamb_count + df$Pm_gamb_count + df$Im_gamb_count
  } else {
    stopf("Could not extract female totals. Available columns: %s", paste(names(df), collapse=", "))
  }
  
  # MALES
  if (!is.null(mosquito_genotypes) && !is.null(mosquito_genotypes$male)) {
    male_total <- rowSums(mosquito_genotypes$male)
  } else if ("male_gamb_count" %in% names(df)) {
    male_total <- df$male_gamb_count
  } else if ("n_male_gamb" %in% names(df)) {
    male_total <- df$n_male_gamb
  } else {
    male_total <- NA_real_  # no-cube path usually does not expose males
  }
  
  data.frame(
    timestep = df[[time_col]],
    female_total = female_total,
    male_total = male_total
  )
}

# -------------------------
# Config
# -------------------------
Tsim <- 365L
seed <- 42L

params_base <- get_parameters(list(
  individual_mosquitoes = TRUE,
  model_seasonality     = FALSE,
  progress_bar          = FALSE,
  human_population      = 400,
  total_M               = 2000,
  init_foim             = 0
))
params_base <- parameterise_total_M(params_base, params_base$total_M)

cube_WT <- MGDrivE::cubeMendelian(gtype = c("AA"))

releases <- list(
  releasesStart    = 60L,
  releasesNumber   = 3L,
  releasesInterval = 7L,
  releaseCount     = 100L,
  releaseSex       = "M",
  releaseGenotype  = "AA"
)

# -------------------------
# (A) Baseline vs WT-only cube: shared-column identity
# -------------------------
cat("=== (A) Baseline vs WT-only cube: shared-column identity check ===\n")

set.seed(seed)
df0 <- run_simulation(Tsim, params_base)

params_wt <- params_base
params_wt$cube <- cube_WT
params_wt$vector_infectivity_g <- c(AA = 1)  # intended no-op

set.seed(seed)
df_wt <- run_simulation(Tsim, params_wt)

assert_shared_columns_identical(df0, df_wt, "no-cube", "WT-only cube")
cat("PASS: shared columns identical (no-cube vs WT-only cube).\n\n")

write_csv_safe(df0,  "out_baseline_no_cube.csv")
write_csv_safe(df_wt, "out_wt_only_cube.csv")

# Write numeric diffs for shared numeric columns only
shared <- intersect(names(df0), names(df_wt))
shared_num <- shared[vapply(shared, function(nm) is.numeric(df0[[nm]]) && is.numeric(df_wt[[nm]]), logical(1))]
diff_shared <- df_wt[, shared_num, drop=FALSE]
for (nm in shared_num) diff_shared[[nm]] <- df_wt[[nm]] - df0[[nm]]
diff_shared[[detect_time_col(df0)]] <- df0[[detect_time_col(df0)]]
write_csv_safe(diff_shared, "out_diff_wt_minus_baseline_shared_numeric.csv")

# -------------------------
# (B) Emergence timing check (only meaningful when male totals exist)
# -------------------------
cat("=== (B) Emergence timing check (WT-only cube resumable run) ===\n")

set.seed(seed)
run_wt2 <- run_resumable_simulation(Tsim, parameters = params_wt)

df_wt2 <- run_wt2$data
geno_wt2 <- run_wt2$mosquito_genotypes

totWT <- extract_adult_totals(df_wt2, geno_wt2)
write_csv_safe(totWT, "adult_totals_wt_cube.csv")

# Early timesteps diagnostic window
early <- totWT %>% filter(timestep <= 3)

if (any(is.na(early$male_total))) {
  stopf("Male totals missing in WT-cube run; cannot test emergence timing.")
}

early <- early %>%
  mutate(rel_diff = (male_total - female_total) / pmax(1, female_total))

cat("Early timesteps (t<=3):\n")
print(early)

# Allow small stochastic deviation; disallow large male-only spikes
assert_true(max(abs(early$rel_diff)) < 0.10,
            "FAIL: Large early male/female discrepancy suggests males are still applied earlier than females.")

cat("PASS: No large early male-only emergence spike (|male-female|/female < 10% for t<=3).\n\n")

# -------------------------
# (C) Releases still work under new queue/commit logic
# -------------------------
cat("=== (C) Release regression check (WT-only cube) ===\n")

params_rel <- params_base
params_rel$cube <- cube_WT
params_rel$vector_infectivity_g <- c(AA = 1)
params_rel <- set_releases(params_rel, releases)

set.seed(seed)
run_rel <- run_resumable_simulation(Tsim, parameters = params_rel)

df_rel <- run_rel$data
geno_rel <- run_rel$mosquito_genotypes

# 1) n_released_gamb exists and totals match
assert_true("n_released_gamb" %in% names(df_rel), "Missing n_released_gamb in release run output.")
expected_total_release <- releases$releaseCount * releases$releasesNumber
observed_total_release <- sum(df_rel$n_released_gamb, na.rm = TRUE)
cat(sprintf("Expected total released: %d | Observed: %d\n", expected_total_release, observed_total_release))
assert_true(observed_total_release == expected_total_release, "Total released mismatch.")

# 2) On first release day, male AA should jump by a material amount
time_col <- detect_time_col(df_rel)
t0 <- releases$releasesStart
idx <- which(df_rel[[time_col]] == t0)
assert_true(length(idx) == 1, "Could not locate releasesStart timestep row in output.")

if (idx > 1) {
  assert_true("AA" %in% colnames(geno_rel$male), "Genotype 'AA' not found in male genotype matrix.")
  male_AA_t0 <- geno_rel$male[idx, "AA"]
  male_AA_prev <- geno_rel$male[idx - 1, "AA"]
  jump <- male_AA_t0 - male_AA_prev
  cat(sprintf("Male AA jump at releaseStart (t=%d): %.1f\n", t0, jump))
  
  # Not necessarily full releaseCount due to same-day ordering/mortality;
  # require at least half as a robust regression signal.
  assert_true(jump >= 0.5 * releases$releaseCount,
              "Release did not produce expected increase in male AA (too small).")
}

cat("PASS: Releases still functioning (schedule + genotype presence).\n\n")

write_csv_safe(df_rel, "out_wt_cube_with_releases.csv")

# Write genotype matrices to CSV for external checking
geno_f <- as.data.frame(geno_rel$female, check.names = FALSE)
geno_m <- as.data.frame(geno_rel$male,   check.names = FALSE)
geno_f[[time_col]] <- df_rel[[time_col]]
geno_m[[time_col]] <- df_rel[[time_col]]

write_csv_safe(geno_f, "geno_female_wt_release.csv")
write_csv_safe(geno_m, "geno_male_wt_release.csv")

cat("ALL TESTS PASSED.\n")

