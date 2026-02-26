# check_WT_cube_bitwise_contents.R

suppressPackageStartupMessages({
  library(malariasimulationGD)
  library(MGDrivE)
})

# ---- Setup ----
year <- 365
sim_length <- 120
human_population <- 1000
starting_EIR <- 5
age_min <- seq(0, 80, 5) * year
age_max <- seq(5, 85, 5) * year

simparams <- get_parameters(list(
  human_population = human_population,
  age_group_rendering_min_ages = age_min,
  age_group_rendering_max_ages = age_max,
  model_seasonality = FALSE,
  individual_mosquitoes = TRUE
))
simparams <- set_equilibrium(simparams, starting_EIR)

cube_WT <- cubeMendelian(gtype = c("AA"))
seed <- 123

# ---- Runs ----
set.seed(seed)
out0 <- run_simulation(sim_length, simparams)

simparams_wt <- simparams
simparams_wt$cube <- cube_WT

set.seed(seed)
out_wt <- run_simulation(sim_length, simparams_wt)

# ---- Check: same columns and same column values (ignoring attributes/rownames) ----
if (!identical(names(out0), names(out_wt))) stop("FAIL: column names/order differ.", call. = FALSE)
if (!identical(nrow(out0), nrow(out_wt))) stop("FAIL: number of rows differ.", call. = FALSE)

# Compare each column with identical()
diffs <- names(out0)[!vapply(names(out0), function(nm) identical(out0[[nm]], out_wt[[nm]]), logical(1))]

if (length(diffs) > 0) {
  cat("First differing column:", diffs[1], "\n")
  stop("FAIL: At least one column differs between baseline and WT cube (same seed).", call. = FALSE)
}

# Also do a full-object check ignoring attributes (includes row.names etc.)
ok_all_equal <- isTRUE(all.equal(out0, out_wt, tolerance = 0, check.attributes = FALSE))
if (!ok_all_equal) {
  stop("FAIL: all.equal(..., check.attributes=FALSE) indicates mismatch (unexpected if all columns identical).", call. = FALSE)
}

cat("PASS: Baseline vs WT-only cube have bitwise-identical column contents (same seed).\n")

# Optional: confirm WT cube genotype attributes exist (baseline may not have them)
if (is.null(attr(out_wt, "mosquito_genotype_counts_female"))) stop("FAIL: WT cube genotype attributes missing.", call. = FALSE)
cat("PASS: WT cube genotype attributes present.\n")

