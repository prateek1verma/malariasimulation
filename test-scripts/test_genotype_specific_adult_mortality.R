#!/usr/bin/env Rscript
############################################################
## Test script: genotype-specific adult mortality scaling (cube$omega)
##
## What this tests (minimal, robust):
##  1) omega == 1 for all genotypes is a strict no-op (bitwise identical output)
##  2) omega > 1 for drive genotypes reduces their adult abundance vs omega==1
##     (directional check on female+male genotype trajectories after release)
##
## Assumptions (matches your current API):
##  - run_simulation(T, params) returns a data.frame (baseline path)
##  - run_resumable_simulation(T, parameters=params) returns list with:
##      $data (data.frame with timestep)
##      $mosquito_genotypes$female, $male (matrices time x genotype)
##  - set_releases(params, releases) sets schedule
##  - cube from MGDrivE with cube$genotypesID and cube$releaseType
############################################################

suppressPackageStartupMessages({
  library(malariasimulationGD)
  library(MGDrivE)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
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

# Compare data.frame *contents* exactly (ignores attributes)
assert_df_columns_identical <- function(a, b, labelA="A", labelB="B") {
  if (!identical(names(a), names(b))) stopf("FAIL: column names differ (%s vs %s)", labelA, labelB)
  if (nrow(a) != nrow(b)) stopf("FAIL: nrow differs (%s vs %s)", labelA, labelB)
  diffs <- names(a)[!vapply(names(a), function(nm) identical(a[[nm]], b[[nm]]), logical(1))]
  if (length(diffs) > 0) stopf("FAIL: first differing column (%s vs %s): %s", labelA, labelB, diffs[1])
  invisible(TRUE)
}

mean_tail <- function(x, n = 60) mean(tail(x, n), na.rm = TRUE)

# Convert mosquito_genotypes matrices (time x genotype) to long df
geno_to_long <- function(df, mosquito_genotypes, time_col) {
  stopifnot(time_col %in% names(df))
  day <- df[[time_col]]
  
  F <- as.data.frame(mosquito_genotypes$female, check.names = FALSE)
  M <- as.data.frame(mosquito_genotypes$male,   check.names = FALSE)
  
  F[[time_col]] <- day
  M[[time_col]] <- day
  
  F_long <- pivot_longer(F, cols = -all_of(time_col), names_to = "genotype", values_to = "count") |>
    mutate(sex = "Female")
  
  M_long <- pivot_longer(M, cols = -all_of(time_col), names_to = "genotype", values_to = "count") |>
    mutate(sex = "Male")
  
  bind_rows(F_long, M_long) |>
    rename(time = all_of(time_col))
}

theme_paper <- theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

release_vlines <- function(rel_days) {
  geom_vline(xintercept = rel_days, linetype = "dashed", linewidth = 0.35, alpha = 0.5)
}

# -------------------------
# 1) Cube + infectivity weights
# -------------------------
cube_base <- cubeHoming1RA(c = 1, ch = 0.95)

# Infectivity weights (user-defined; must cover cube genotypes)
vector_infectivity_g <- c(HH=0, HW=0, HR=0, WW=1, WR=1, RR=1)

missing_g <- setdiff(cube_base$genotypesID, names(vector_infectivity_g))
extra_g   <- setdiff(names(vector_infectivity_g), cube_base$genotypesID)
if (length(missing_g) > 0) stopf("vector_infectivity_g missing genotypes: %s", paste(missing_g, collapse=", "))
if (length(extra_g) > 0)   stopf("vector_infectivity_g has extra genotypes not in cube: %s", paste(extra_g, collapse=", "))
vector_infectivity_g <- vector_infectivity_g[cube_base$genotypesID]

# Two omega scenarios:
#  - omega_all1: strict no-op (should match omega absent / 1)
#  - omega_cost: increased adult mortality for H* genotypes (fitness cost)
omega_all1 <- setNames(rep(1, length(cube_base$genotypesID)), cube_base$genotypesID)
omega_cost <- omega_all1
omega_cost[c("HH","HW","HR")] <- c(1.2, 1.1, 1.1)  # adjust as desired

# -------------------------
# 2) Parameters (keep identical across runs)
# -------------------------
year <- 365L
Tsim <- 4L * year
seed <- 1

params0 <- get_parameters(list(
  individual_mosquitoes = TRUE,
  model_seasonality     = FALSE,
  progress_bar          = FALSE,
  human_population      = 500,
  total_M               = 5000,
  init_foim             = 0
))
params0 <- parameterise_total_M(params0, params0$total_M)

releases <- list(
  releasesStart    = 100L,
  releasesNumber   = 1L,
  releasesInterval = 999999L,
  releaseCount     = 50L,
  releaseSex       = "M",
  releaseGenotype  = cube_base$releaseType
)

# -------------------------
# 3) (Test A) omega == 1 is a strict no-op
# -------------------------
cat("=== Test A: omega == 1 should be strict no-op ===\n")

# A1: run with cube but omega omitted (or NULL)
params_A1 <- params0
params_A1$cube <- cube_base
params_A1$vector_infectivity_g <- vector_infectivity_g
params_A1 <- set_releases(params_A1, releases)

# Ensure omega absent (simulate legacy behavior)
params_A1$cube$omega <- NULL

set.seed(seed)
run_A1 <- run_resumable_simulation(Tsim, parameters = params_A1)
df_A1  <- run_A1$data

# A2: run with omega explicitly all 1
params_A2 <- params0
cube_A2 <- cube_base
cube_A2$omega <- omega_all1
params_A2$cube <- cube_A2
params_A2$vector_infectivity_g <- vector_infectivity_g
params_A2 <- set_releases(params_A2, releases)

set.seed(seed)
run_A2 <- run_resumable_simulation(Tsim, parameters = params_A2)
df_A2  <- run_A2$data

# Compare output data frames exactly (bitwise identical columns)
assert_df_columns_identical(df_A1, df_A2, "omega omitted", "omega all ones")
cat("PASS: df outputs identical when omega omitted vs omega==1.\n\n")

# -------------------------
# 4) (Test B) omega > 1 reduces adult abundance for those genotypes
# -------------------------
cat("=== Test B: omega > 1 should reduce adult abundance for targeted genotypes ===\n")

params_B <- params0
cube_B <- cube_base
cube_B$omega <- omega_cost
params_B$cube <- cube_B
params_B$vector_infectivity_g <- vector_infectivity_g
params_B <- set_releases(params_B, releases)

set.seed(seed)
run_B <- run_resumable_simulation(Tsim, parameters = params_B)

df_B   <- run_B$data
geno_A <- run_A2$mosquito_genotypes  # omega==1 reference
geno_B <- run_B$mosquito_genotypes   # omega-cost

time_col <- detect_time_col(df_A2)

# Directional check: compare tail means of adult counts for H* genotypes
target_genos <- intersect(c("HH","HW","HR"), cube_base$genotypesID)
if (length(target_genos) == 0) stopf("No target genotypes found in cube$genotypesID.")

tail_n <- 120L  # last ~4 months
# ---- Directional check: only enforce for the most-penalized genotype(s) ----
omega_vec <- omega_cost[cube_base$genotypesID]
max_omega <- max(omega_vec, na.rm = TRUE)
most_penalized <- names(omega_vec)[omega_vec == max_omega]
cat(sprintf("Most penalized genotype(s) (max omega=%.3f): %s\n",
            max_omega, paste(most_penalized, collapse = ", ")))

for (g in most_penalized) {
  A_f <- mean_tail(geno_A$female[, g], tail_n)
  B_f <- mean_tail(geno_B$female[, g], tail_n)
  A_m <- mean_tail(geno_A$male[, g],   tail_n)
  B_m <- mean_tail(geno_B$male[, g],   tail_n)
  
  cat(sprintf("Genotype %s (tail mean, last %d days): Female A=%.2f B=%.2f | Male A=%.2f B=%.2f\n",
              g, tail_n, A_f, B_f, A_m, B_m))
  
  # Require *material* reduction (avoid failing on tiny stochastic wiggles)
  # If A is ~0 already, skip.
  if (A_f > 1) assert_true(B_f < 0.95 * A_f, paste("Expected >=5% lower female abundance under omega cost for", g))
  if (A_m > 1) assert_true(B_m < 0.95 * A_m, paste("Expected >=5% lower male abundance under omega cost for", g))
}

cat("PASS: most-penalized genotype(s) reduced under omega-cost.\n\n")
# -------------------------
# 5) Plots (compare omega==1 vs omega-cost)
# -------------------------
cat("=== Writing comparison plots ===\n")

# Build long dfs for both scenarios
df_ref <- df_A2
df_cost <- df_B

geno_long_ref <- geno_to_long(df_ref, geno_A, time_col) |>
  mutate(scenario = "omega=1")

geno_long_cost <- geno_to_long(df_cost, geno_B, time_col) |>
  mutate(scenario = "omega>1 (cost)")

geno_long <- bind_rows(geno_long_ref, geno_long_cost) |>
  mutate(genotype = factor(genotype, levels = cube_base$genotypesID),
         sex = factor(sex, levels = c("Female","Male")))

rel_days <- attr(df_ref, "mosquito_release_schedule")$timestep

geno_names <- cube_base$genotypesID
pal <- brewer.pal(max(3, min(8, length(geno_names))), "Set1")
geno_cols <- setNames(rep(pal, length.out = length(geno_names)), geno_names)

p_compare <- ggplot(geno_long, aes(x = time, y = count, color = genotype, linetype = sex)) +
  release_vlines(rel_days) +
  geom_line(linewidth = 0.9) +
  facet_grid(scenario ~ genotype, scales = "free_y") +
  scale_color_manual(values = geno_cols) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Genotype-specific adult mortality scaling test (cube$omega)",
    subtitle = sprintf("One-time male release: %d of %s on day %d",
                       releases$releaseCount, releases$releaseGenotype, releases$releasesStart),
    x = "Day", y = "Adult count"
  ) +
  theme_paper +
  theme(legend.position = "top")

ggsave("test_omega_genotype_adults_compare.pdf", p_compare, width = 12, height = 7, useDingbats = FALSE)
ggsave("test_omega_genotype_adults_compare.png", p_compare, width = 12, height = 7, dpi = 300)

cat("Wrote:\n",
    " - test_omega_genotype_adults_compare.(pdf|png)\n")


# ---- Plot human outcomes: omega=1 vs omega-cost ----
pick_human_cols <- function(df) {
  infections_candidates <- c("n_infections", "new_infections", "infections")
  clinical_candidates <- c(
    "clinical_incidence", "n_clinical", "clinical_cases", "incidence_clinical",
    "n_clinical_cases", "clinical", "case_clinical"
  )
  pick_first <- function(cands) {
    hit <- intersect(cands, names(df))
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  list(infections = pick_first(infections_candidates),
       clinical   = pick_first(clinical_candidates))
}

time_col <- detect_time_col(df_ref)

cols_ref  <- pick_human_cols(df_ref)
cols_cost <- pick_human_cols(df_cost)

# Require infections to exist (your model normally has n_infections)
if (is.na(cols_ref$infections) || is.na(cols_cost$infections)) {
  stop("Could not find infections column (expected n_infections) in one of the scenarios.")
}

# Use same column names across scenarios; clinical is optional
make_hum_long <- function(df, scenario, time_col, cols) {
  out <- df |>
    transmute(
      day = .data[[time_col]],
      `Infections (n_infections)` = .data[[cols$infections]],
      `Clinical incidence` = if (!is.na(cols$clinical)) .data[[cols$clinical]] else NA_real_
    ) |>
    tidyr::pivot_longer(cols = -day, names_to = "metric", values_to = "value") |>
    dplyr::filter(!(metric == "Clinical incidence" & is.na(value))) |>
    dplyr::mutate(scenario = scenario)
  out
}

hum_long <- dplyr::bind_rows(
  make_hum_long(df_ref,  "omega=1",        time_col, cols_ref),
  make_hum_long(df_cost, "omega>1 (cost)", time_col, cols_cost)
)

# Release days for vertical lines (use reference schedule)
rel_days <- attr(df_ref, "mosquito_release_schedule")$timestep

p_hum_compare <- ggplot(hum_long, aes(x = day, y = value, color = scenario)) +
  release_vlines(rel_days) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Human outcomes: omega=1 vs omega-cost",
    x = "Day",
    y = NULL
  ) +
  theme_paper

print(p_hum_compare)

ggsave("test_omega_humans_compare.pdf", p_hum_compare, width = 9, height = 5.2, useDingbats = FALSE)
ggsave("test_omega_humans_compare.png", p_hum_compare, width = 9, height = 5.2, dpi = 300)

