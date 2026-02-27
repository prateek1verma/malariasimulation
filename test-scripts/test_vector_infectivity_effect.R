#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly=TRUE)) {
    pkgload::load_all(".", quiet=TRUE)
  } else {
    library(malariasimulationGD)
  }
  library(MGDrivE)
})

stopf <- function(...) stop(sprintf(...), call.=FALSE)
assert_true <- function(x, msg) if (!isTRUE(x)) stop(msg, call.=FALSE)

# Compare data.frame *column contents* exactly (ignores attributes)
assert_df_columns_identical <- function(a, b, labelA="A", labelB="B") {
  if (!identical(names(a), names(b))) stopf("FAIL: column names differ (%s vs %s)", labelA, labelB)
  if (nrow(a) != nrow(b)) stopf("FAIL: nrow differs (%s vs %s)", labelA, labelB)
  diffs <- names(a)[!vapply(names(a), function(nm) identical(a[[nm]], b[[nm]]), logical(1))]
  if (length(diffs) > 0) stopf("FAIL: first differing column (%s vs %s): %s", labelA, labelB, diffs[1])
  invisible(TRUE)
}

mean_tail <- function(x, n=60) mean(tail(x, n), na.rm=TRUE)

# -------------------------
# Setup
# -------------------------
T <- 365L
seed <- 42

params_base <- get_parameters(list(
  individual_mosquitoes = TRUE,
  model_seasonality = FALSE,
  progress_bar = FALSE,
  human_population = 400,
  total_M = 3000,
  init_foim = 0
))
params_base <- parameterise_total_M(params_base, params_base$total_M)

cube_WT <- MGDrivE::cubeMendelian(gtype=c("AA"))
cube_3  <- MGDrivE::cubeMendelian(gtype=c("AA","Aa","aa"))
cube_3$releaseType <- "aa"

releases <- list(
  releasesStart = 60L,
  releasesNumber = 12L,
  releasesInterval = 7L,
  releaseCount = 400L,
  releaseSex = "M"
)

cat("=== (1) Baseline vs WT-only cube no-op check ===\n")
set.seed(seed)
out0 <- run_simulation(T, params_base)

params_wt <- params_base
params_wt$cube <- cube_WT
params_wt$vector_infectivity_g <- c(AA=1)  # should be strict no-op path per implementation note

set.seed(seed)
out_wt <- run_simulation(T, params_wt)

# WT-only cube with AA=1 should not add diagnostics columns either (per your note).
assert_df_columns_identical(out0, out_wt, "baseline", "WT cube + vector_infectivity_g(AA=1)")
cat("PASS: baseline vs WT cube(AA) + AA=1 are bitwise identical in dataframe columns.\n\n")

cat("=== (2) Mendelian cube + releases: control (all genotypes infective) ===\n")
params_ctl <- params_base
params_ctl$cube <- cube_3
params_ctl <- set_releases(params_ctl, releases)
params_ctl$vector_infectivity_g <- c(AA=1, Aa=1, aa=1)

set.seed(seed)
run_ctl <- run_resumable_simulation(T, parameters=params_ctl)
df_ctl <- run_ctl$data
geno_ctl <- run_ctl$mosquito_genotypes

# Ensure feature columns exist when active
req_cols <- c("infectivity_weighted_I_gamb", "vector_infectivity_mean_gamb", "EIR_gamb", "n_infections")
for (nm in req_cols) assert_true(nm %in% names(df_ctl), paste("Missing expected column:", nm))
cat("PASS: feature diagnostics columns present in control run.\n\n")

cat("=== (3) Mendelian cube + releases: blocking (AA=1, Aa=0, aa=0) ===\n")
params_blk <- params_base
params_blk$cube <- cube_3
params_blk <- set_releases(params_blk, releases)
params_blk$vector_infectivity_g <- c(AA=1, Aa=0, aa=0)

set.seed(seed)
run_blk <- run_resumable_simulation(T, parameters=params_blk)
df_blk <- run_blk$data
geno_blk <- run_blk$mosquito_genotypes

for (nm in req_cols) assert_true(nm %in% names(df_blk), paste("Missing expected column:", nm))
cat("PASS: feature diagnostics columns present in blocking run.\n\n")

# -------------------------
# Checks: genotype dynamics and directional effects
# -------------------------
cat("=== (4) Directional checks ===\n")

# Sanity: releases occurred as expected
assert_true("n_released_gamb" %in% names(df_blk), "Missing n_released_gamb")
expected_total_release <- releases$releaseCount * releases$releasesNumber
assert_true(sum(df_blk$n_released_gamb, na.rm=TRUE) == expected_total_release,
            "Total released mismatch in blocking run")

# Genotype response: after repeated aa male releases, Aa or aa should appear in adults (often females too)
# (May depend on your aquatic timings; use a weak check: at least males have aa > 0 after first release)
first_rel <- releases$releasesStart
assert_true(geno_blk$male[first_rel, "aa"] > 0, "Expected male aa > 0 on first release day")

# Core effect: vector_infectivity_mean should drop below 1 after releases in blocking run
m_ctl <- mean_tail(df_ctl$vector_infectivity_mean_gamb, 60)
m_blk <- mean_tail(df_blk$vector_infectivity_mean_gamb, 60)
cat(sprintf("Mean vector_infectivity_mean_gamb (last 60d): control=%.3f, blocking=%.3f\n", m_ctl, m_blk))
assert_true(m_blk < m_ctl - 0.05, "Expected blocking run to reduce vector_infectivity_mean_gamb vs control")

# Effective infectious count should drop relative to control
Ieff_ctl <- mean_tail(df_ctl$infectivity_weighted_I_gamb, 60)
Ieff_blk <- mean_tail(df_blk$infectivity_weighted_I_gamb, 60)
cat(sprintf("Mean infectivity_weighted_I_gamb (last 60d): control=%.3f, blocking=%.3f\n", Ieff_ctl, Ieff_blk))
assert_true(Ieff_blk < Ieff_ctl, "Expected infectivity_weighted_I_gamb to be lower in blocking run")

# EIR should drop (since EIR storage uses weighted I)
EIR_ctl <- mean_tail(df_ctl$EIR_gamb, 60)
EIR_blk <- mean_tail(df_blk$EIR_gamb, 60)
cat(sprintf("Mean EIR_gamb (last 60d): control=%.3f, blocking=%.3f\n", EIR_ctl, EIR_blk))
assert_true(EIR_blk < EIR_ctl, "Expected EIR_gamb to be lower in blocking run")

# Human infections should drop (directional; stochastic so allow modest margin)
inf_ctl <- mean_tail(df_ctl$n_infections, 60)
inf_blk <- mean_tail(df_blk$n_infections, 60)
cat(sprintf("Mean n_infections (last 60d): control=%.3f, blocking=%.3f\n", inf_ctl, inf_blk))
assert_true(inf_blk < inf_ctl, "Expected n_infections to be lower in blocking run")

cat("PASS: directional effects consistent with AA=1, Aa=0, aa=0 infectivity.\n\n")

# -------------------------
# Optional plots (ggplot2 + RColorBrewer)
# -------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
})

cat("=== (5) Writing paper-quality ggplot diagnostics ===\n")
rel_days <- attr(df_blk, "mosquito_release_schedule")$timestep

# Colors (ColorBrewer Reds for comparison; use 2 darkest)
cmp_cols <- brewer.pal(9, "Reds")[c(7, 9)]
names(cmp_cols) <- c("control", "blocking")

theme_paper <- theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# Helper to add release-day vertical lines
release_vlines <- function(rel_days) {
  geom_vline(xintercept = rel_days, linetype = "dashed", linewidth = 0.4, alpha = 0.5)
}

# ---- (a) vector_infectivity_mean_gamb ----
df_vi <- rbind(
  data.frame(day = df_ctl$timestep, value = df_ctl$vector_infectivity_mean_gamb, scenario = "control"),
  data.frame(day = df_blk$timestep, value = df_blk$vector_infectivity_mean_gamb, scenario = "blocking")
)

p1 <- ggplot(df_vi, aes(day, value, color = scenario)) +
  release_vlines(rel_days) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = cmp_cols) +
  labs(
    title = "Genotype-weighted mean infectivity among infectious mosquitoes",
    x = "Day",
    y = "Mean infectivity (vector_infectivity_mean_gamb)"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_paper

ggsave("vector_infectivity_mean_gamb.png", p1, width = 9, height = 4.8, dpi = 300)
ggsave("vector_infectivity_mean_gamb.pdf", p1, width = 9, height = 4.8)

# ---- (b) EIR_gamb ----
df_eir <- rbind(
  data.frame(day = df_ctl$timestep, value = df_ctl$EIR_gamb, scenario = "control"),
  data.frame(day = df_blk$timestep, value = df_blk$EIR_gamb, scenario = "blocking")
)

p2 <- ggplot(df_eir, aes(day, value, color = scenario)) +
  release_vlines(rel_days) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = cmp_cols) +
  labs(
    title = "EIR trajectory (gambiae): control vs blocking",
    x = "Day",
    y = "EIR_gamb"
  ) +
  theme_paper

ggsave("EIR_gamb_compare.png", p2, width = 9, height = 4.8, dpi = 300)
ggsave("EIR_gamb_compare.pdf", p2, width = 9, height = 4.8)

# ---- (c) n_infections ----
df_inf <- rbind(
  data.frame(day = df_ctl$timestep, value = df_ctl$n_infections, scenario = "control"),
  data.frame(day = df_blk$timestep, value = df_blk$n_infections, scenario = "blocking")
)

p3 <- ggplot(df_inf, aes(day, value, color = scenario)) +
  release_vlines(rel_days) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = cmp_cols) +
  labs(
    title = "Human infections per day: control vs blocking",
    x = "Day",
    y = "n_infections"
  ) +
  theme_paper

ggsave("n_infections_compare.png", p3, width = 9, height = 4.8, dpi = 300)
ggsave("n_infections_compare.pdf", p3, width = 9, height = 4.8)

# ---- (d) Female genotype trajectories (blocking run) ----
F <- geno_blk$female
df_F <- data.frame(day = df_blk$timestep, F, check.names = FALSE)

# reshape to long without extra deps
geno_names <- setdiff(colnames(df_F), "day")
df_F_long <- do.call(rbind, lapply(geno_names, function(g) {
  data.frame(day = df_F$day, genotype = g, value = df_F[[g]])
}))

# Brewer palette for genotypes (Set1 works well)
geno_cols <- brewer.pal(max(3, length(geno_names)), "Set1")[seq_along(geno_names)]
names(geno_cols) <- geno_names

p4 <- ggplot(df_F_long, aes(day, value, color = genotype)) +
  release_vlines(rel_days) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = geno_cols) +
  labs(
    title = "Adult female mosquitoes by genotype (blocking run)",
    x = "Day",
    y = "Adult female count"
  ) +
  theme_paper

ggsave("female_genotypes_blocking.png", p4, width = 9, height = 4.8, dpi = 300)
ggsave("female_genotypes_blocking.pdf", p4, width = 9, height = 4.8)

cat("Wrote plots (PNG + PDF):\n",
    " - vector_infectivity_mean_gamb.(png|pdf)\n",
    " - EIR_gamb_compare.(png|pdf)\n",
    " - n_infections_compare.(png|pdf)\n",
    " - female_genotypes_blocking.(png|pdf)\n")

