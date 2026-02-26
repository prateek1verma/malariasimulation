# test_genotype_effect_blackbox.R
suppressPackageStartupMessages({
  library(malariasimulationGD)
  library(MGDrivE)
})

year <- 365
sim_length <- 200
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

# Cubes
cube_3 <- cubeMendelian(gtype = c("AA", "Aa", "aa"))
cube_sterile <- cube_3
cube_sterile$tau[] <- 0  # everything non-viable => V(t) should go to 0 if implemented

# Runs
seed <- 123
set.seed(seed)
out0 <- run_simulation(sim_length, simparams)

simparams3 <- simparams; simparams3$cube <- cube_3
set.seed(seed)
out3 <- run_simulation(sim_length, simparams3)

simparamsS <- simparams; simparamsS$cube <- cube_sterile
set.seed(seed)
outS <- run_simulation(sim_length, simparamsS)

# Quick numeric checks (choose robust mosquito columns)
mosq_cols <- c("E_gamb_count","L_gamb_count","P_gamb_count","Sm_gamb_count","Pm_gamb_count","Im_gamb_count")
missing <- mosq_cols[!mosq_cols %in% names(out0)]
if (length(missing) > 0) stop("Missing expected columns: ", paste(missing, collapse=", "))

end <- nrow(out0)

cat("\nFinal day mosquito counts:\n")
print(rbind(
  baseline = unlist(out0[end, mosq_cols]),
  mendelian_cube = unlist(out3[end, mosq_cols]),
  sterile_tau0 = unlist(outS[end, mosq_cols])
))

# Sterile cube should have MUCH smaller counts if V(t) is applied
ratio_E <- (outS[end, "E_gamb_count"] + 1) / (out0[end, "E_gamb_count"] + 1)
ratio_Sm <- (outS[end, "Sm_gamb_count"] + 1) / (out0[end, "Sm_gamb_count"] + 1)
cat(sprintf("\nSterile/baseline ratios: E=%.4g, Sm=%.4g\n", ratio_E, ratio_Sm))

# Plots: baseline vs sterile
par(mfrow=c(2,3), mar=c(4,4,2,1))
for (nm in mosq_cols) {
  plot(out0$timestep, out0[[nm]], type="l", xlab="timestep", ylab=nm, main=nm)
  lines(outS$timestep, outS[[nm]], lty=2)
  legend("topright", legend=c("baseline","sterile tau=0"), lty=c(1,2), bty="n", cex=0.8)
}

cat("\nInterpretation:\n")
cat("- If sterile tau=0 collapses E/L/P and adults vs baseline => genotype->V(t) scaling is working.\n")
cat("- If sterile tau=0 looks identical to baseline => cube is not being used in dynamics (or V not applied).\n")






# test_genotype_patch.R
# Run: Rscript test_genotype_patch.R

suppressPackageStartupMessages({
  library(malariasimulationGD)
  library(MGDrivE)
})

cat("malariasimulationGD version:", as.character(packageVersion("malariasimulationGD")), "\n")

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
cube_3  <- cubeMendelian(gtype = c("AA", "Aa", "aa"))

seed <- 123

get_geno_attrs <- function(out) {
  list(
    F = attr(out, "mosquito_genotype_counts_female"),
    M = attr(out, "mosquito_genotype_counts_male"),
    V = attr(out, "mosquito_genotype_V"),
    Tad = attr(out, "mosquito_genotype_total_adults")
  )
}

check_geno_consistency <- function(out, cube, label) {
  a <- get_geno_attrs(out)
  if (any(vapply(a, is.null, logical(1)))) stop("Missing genotype attributes for: ", label, call. = FALSE)
  
  F <- a$F; M <- a$M; V <- a$V; Tad <- a$Tad
  
  if (!identical(colnames(F), cube$genotypesID)) stop(label, ": female genotype colnames != cube$genotypesID", call. = FALSE)
  if (!identical(colnames(M), cube$genotypesID)) stop(label, ": male genotype colnames != cube$genotypesID", call. = FALSE)
  
  # Totals
  tot_from_geno <- rowSums(F) + rowSums(M)
  if (!isTRUE(all.equal(as.numeric(tot_from_geno), as.numeric(Tad), tolerance = 0))) {
    stop(label, ": total adults != sum(F)+sum(M)", call. = FALSE)
  }
  
  # Female state totals (if present)
  if (all(c("Sm_gamb_count","Pm_gamb_count","Im_gamb_count") %in% names(out))) {
    female_states <- out$Sm_gamb_count + out$Pm_gamb_count + out$Im_gamb_count
    if (!isTRUE(all.equal(as.numeric(rowSums(F)), as.numeric(female_states), tolerance = 0))) {
      stop(label, ": rowSums(F) != Sm+Pm+Im", call. = FALSE)
    }
  }
  
  # Mendelian cube tau=1 => V(t)=1 (no deterministic viability effect)
  Vv <- if (is.matrix(V) || is.data.frame(V)) as.numeric(V[, 1]) else as.numeric(V)
  if (any(abs(Vv - 1) > 1e-12)) stop(label, ": expected V(t)=1 for Mendelian cube (tau=1)", call. = FALSE)
  
  invisible(TRUE)
}

check_all_WT <- function(out, cube, label) {
  a <- get_geno_attrs(out)
  F <- a$F; M <- a$M
  wt <- cube$wildType
  non_wt <- setdiff(colnames(F), wt)
  if (length(non_wt) > 0) {
    if (any(F[, non_wt, drop=FALSE] != 0) || any(M[, non_wt, drop=FALSE] != 0)) {
      stop(label, ": expected non-WT genotype counts to be 0 (no release), but found nonzero.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

compare_exact_df <- function(outA, outB, labelA, labelB) {
  common <- intersect(names(outA), names(outB))
  ok <- isTRUE(all.equal(outA[common], outB[common], tolerance = 0))
  if (!ok) {
    diffs <- common[!vapply(common, function(nm) isTRUE(all.equal(outA[[nm]], outB[[nm]], tolerance = 0)), logical(1))]
    stop("Exact outputs differ between ", labelA, " and ", labelB, ". First differing column: ", diffs[1], call. = FALSE)
  }
  invisible(TRUE)
}

compare_summary <- function(outA, outB, cols, tail_n = 30, max_rel = 0.05, labelA="A", labelB="B") {
  # Compare mean over last tail_n days; allow e.g. 5% relative difference by default
  n <- min(nrow(outA), nrow(outB))
  idx <- (n - tail_n + 1):n
  for (nm in cols) {
    if (!(nm %in% names(outA)) || !(nm %in% names(outB))) next
    a <- mean(outA[[nm]][idx], na.rm = TRUE)
    b <- mean(outB[[nm]][idx], na.rm = TRUE)
    denom <- max(1e-12, abs(a), abs(b))
    rel <- abs(a - b) / denom
    if (rel > max_rel) {
      stop(sprintf("Summary check failed for %s: mean(last %d) differs too much (%s=%.4g, %s=%.4g, rel=%.3g > %.3g)",
                   nm, tail_n, labelA, a, labelB, b, rel, max_rel), call. = FALSE)
    }
  }
  invisible(TRUE)
}

# ---- 1) Baseline ----
set.seed(seed)
out0 <- run_simulation(sim_length, simparams)
cat("Baseline OK.\n")

# ---- 2) WT cube (AA only): should be exact no-op (no extra RNG) ----
simparams_wt <- simparams; simparams_wt$cube <- cube_WT
set.seed(seed)
out_wt <- run_simulation(sim_length, simparams_wt)
cat("WT cube OK.\n")

compare_exact_df(out0, out_wt, "baseline", "WT cube")
check_geno_consistency(out_wt, cube_WT, "WT cube")
check_all_WT(out_wt, cube_WT, "WT cube")
cat("PASS: WT cube exact no-op + genotype attributes consistent.\n")

# ---- 3) 3-genotype Mendelian cube (no release): RNG differs; test invariants, not path equality ----
simparams_3 <- simparams; simparams_3$cube <- cube_3
set.seed(seed)
out_3 <- run_simulation(sim_length, simparams_3)
cat("3-genotype cube OK.\n")

check_geno_consistency(out_3, cube_3, "3-genotype Mendelian cube")
check_all_WT(out_3, cube_3, "3-genotype Mendelian cube (no release)")
cat("PASS: 3-genotype Mendelian cube has consistent genotype outputs; non-WT counts are zero (no release); V(t)=1.\n")

# Optional: sanity summary check that mosquito metrics are in the same ballpark (not identical)
# Summary sanity check: use only large-count compartments
mosq_cols_robust <- c("E_gamb_count","L_gamb_count","P_gamb_count","Sm_gamb_count")
compare_summary(out0, out_3, mosq_cols_robust, tail_n = 30, max_rel = 0.10,
                labelA="baseline", labelB="3-genotype")
cat("PASS: Robust mosquito compartments within 10% over last 30 days.\n")

# Optional: low-count compartments use a looser rule (absolute diff on mean)
mean_abs_check <- function(outA, outB, cols, tail_n = 30, max_abs = 20, labelA="A", labelB="B") {
  n <- min(nrow(outA), nrow(outB))
  idx <- (n - tail_n + 1):n
  for (nm in cols) {
    if (!(nm %in% names(outA)) || !(nm %in% names(outB))) next
    a <- mean(outA[[nm]][idx], na.rm = TRUE)
    b <- mean(outB[[nm]][idx], na.rm = TRUE)
    if (abs(a - b) > max_abs) {
      stop(sprintf("Mean abs check failed for %s: %s=%.4g, %s=%.4g, abs diff=%.4g > %.4g",
                   nm, labelA, a, labelB, b, abs(a-b), max_abs), call. = FALSE)
    }
  }
  invisible(TRUE)
}

# For Im/Pm counts, allow e.g. ±20 mean difference over last 30 days
mean_abs_check(out0, out_3, c("Pm_gamb_count","Im_gamb_count"), tail_n = 30, max_abs = 20,
               labelA="baseline", labelB="3-genotype")
cat("PASS: Low-count Pm/Im means within absolute tolerance.\n")
cat("\nAll tests PASSED.\n")
