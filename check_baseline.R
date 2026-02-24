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

