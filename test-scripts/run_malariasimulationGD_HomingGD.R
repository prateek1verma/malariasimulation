#!/usr/bin/env Rscript
############################################################
## malariasimulationGD: Arbitrary-G MGDrivE cube + 2y run
## Correct usage for YOUR current API:
##  - run_resumable_simulation(T, parameters=params) for genotype matrices
##  - df uses time column: "timestep"
##  - genotype time series live in run$mosquito_genotypes$female / $male
##  - releases set via set_releases(params, releases)
##  - vector_infectivity_g named by cube$genotypesID (auto-checked)
##
## Outputs:
##  (1) Paper-quality genotype-wise adult females + males vs time (single figure)
##  (2) Infected humans vs time (uses n_infections as a robust default)
##  (3) Clinical incidence vs time (auto-detects a likely column if present)
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

detect_time_col <- function(df) {
  candidates <- c("timestep", "t", "time", "day")
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) stopf("Could not detect time column. Columns are: %s", paste(names(df), collapse=", "))
  hit[1]
}

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

# Robust human metric picker:
# - Always plots n_infections if present (your code uses it)
# - Also tries to find a "clinical incidence" column if present
pick_human_cols <- function(df) {
  # Most reliable from your current outputs
  infections_candidates <- c("n_infections", "new_infections", "infections")
  
  # Clinical-ish candidates (depends on your model naming)
  clinical_candidates <- c(
    "clinical_incidence", "n_clinical", "clinical_cases", "incidence_clinical",
    "n_clinical_cases", "clinical", "case_clinical"
  )
  
  pick_first <- function(cands) {
    hit <- intersect(cands, names(df))
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  
  list(
    infections = pick_first(infections_candidates),
    clinical   = pick_first(clinical_candidates)
  )
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
# 1) Cube (arbitrary G via MGDrivE)
# -------------------------
# Example cube: Homing 1RA (arbitrary genotype set comes from MGDrivE)
cube <- cubeHoming1RA(c = 1, ch = 0.95)

# -------------------------
# 2) Genotype-specific vector infectivity (user input)
# -------------------------
# Must be a named vector covering exactly cube$genotypesID (order can differ)
vector_infectivity_g <- c(HH=0, HW=0, HR=0, WW=1, WR=1, RR=1)

missing_g <- setdiff(cube$genotypesID, names(vector_infectivity_g))
extra_g   <- setdiff(names(vector_infectivity_g), cube$genotypesID)
if (length(missing_g) > 0) stopf("vector_infectivity_g missing genotypes: %s", paste(missing_g, collapse=", "))
if (length(extra_g) > 0)   stopf("vector_infectivity_g has extra genotypes not in cube: %s", paste(extra_g, collapse=", "))

# Reorder to cube genotype order (nice hygiene)
vector_infectivity_g <- vector_infectivity_g[cube$genotypesID]

# -------------------------
# 3) Parameters (2 years)
# -------------------------
year <- 365L
Tsim <- 2L * year

params <- get_parameters(list(
  individual_mosquitoes = TRUE,
  model_seasonality     = FALSE,
  progress_bar          = FALSE,
  human_population      = 500,
  total_M               = 5000,
  init_foim             = 0
))
params <- parameterise_total_M(params, params$total_M)

params$cube <- cube
params$vector_infectivity_g <- vector_infectivity_g

# -------------------------
# 4) One-time male-only release (small)
# -------------------------
releases <- list(
  releasesStart    = 100L,
  releasesNumber   = 1L,
  releasesInterval = 999999L,  # irrelevant since releasesNumber=1
  releaseCount     = 50L,
  releaseSex       = "M",
  releaseGenotype  = cube$releaseType  # default; can override e.g. "RR"
)

params <- set_releases(params, releases)

# -------------------------
# 5) Run (use resumable API so we ALWAYS get genotype matrices)
# -------------------------
set.seed(1)
run <- run_resumable_simulation(Tsim, parameters = params)

df   <- run$data
geno <- run$mosquito_genotypes

time_col <- detect_time_col(df)
rel_days <- attr(df, "mosquito_release_schedule")$timestep

# -------------------------
# 6) Prepare plotting data
# -------------------------
geno_long <- geno_to_long(df, geno, time_col = time_col) |>
  mutate(genotype = factor(genotype, levels = cube$genotypesID),
         sex      = factor(sex, levels = c("Female", "Male")))

# Human outcomes
cols <- pick_human_cols(df)
if (is.na(cols$infections)) stopf("Could not find infections column (expected n_infections). Available: %s",
                                  paste(names(df), collapse=", "))

hum_df <- df |>
  transmute(
    day = .data[[time_col]],
    `Infections (n_infections)` = .data[[cols$infections]],
    `Clinical incidence` = if (!is.na(cols$clinical)) .data[[cols$clinical]] else NA_real_
  ) |>
  pivot_longer(cols = -day, names_to = "metric", values_to = "value") |>
  filter(!(metric == "Clinical incidence" & is.na(value)))

# -------------------------
# 7) Colors (paper-friendly)
# -------------------------
geno_names <- cube$genotypesID
pal <- brewer.pal(max(3, min(8, length(geno_names))), "Set1")
geno_cols <- setNames(rep(pal, length.out = length(geno_names)), geno_names)

# -------------------------
# 8) Plot A: Adult mosquitoes by genotype (Female + Male)
# -------------------------
p_adults <- ggplot(geno_long, aes(x = time, y = count, color = genotype, linetype = sex)) +
  release_vlines(rel_days) +
  geom_line(linewidth = 0.95) +
  facet_wrap(~ genotype, scales = "free_y", ncol = 2) +
  scale_color_manual(values = geno_cols) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Adult mosquitoes by genotype (female & male)",
    x = "Day",
    y = "Adult count",
    subtitle = sprintf(
      "One-time male release: %d of genotype %s on day %d | vector_infectivity_g supplied",
      releases$releaseCount, releases$releaseGenotype, releases$releasesStart
    )
  ) +
  theme_paper +
  theme(legend.position = "top")

# -------------------------
# 9) Plot B: Human outcomes (infections + clinical incidence if present)
# -------------------------
p_humans <- ggplot(hum_df, aes(x = day, y = value)) +
  release_vlines(rel_days) +
  geom_line(linewidth = 1.0) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Human outcomes over time",
    x = "Day",
    y = NULL
  ) +
  theme_paper

# Print
print(p_adults)
print(p_humans)

# -------------------------
# 10) Save (paper-ready)
# -------------------------
ggsave("fig_adults_by_genotype.pdf", p_adults, width = 9, height = 6.5, useDingbats = FALSE)
ggsave("fig_human_outcomes.pdf",    p_humans, width = 9, height = 5.2, useDingbats = FALSE)

# Optional PNGs
ggsave("fig_adults_by_genotype.png", p_adults, width = 9, height = 6.5, dpi = 300)
ggsave("fig_human_outcomes.png",     p_humans, width = 9, height = 5.2, dpi = 300)

cat("Wrote:\n",
    " - fig_adults_by_genotype.(pdf|png)\n",
    " - fig_human_outcomes.(pdf|png)\n")

