test_that("hybrid mosquito genotype tracking is backward compatible and consistent", {
  skip_if_not_installed("MGDrivE")

  base_params <- malariasimulationGD::get_parameters(list(
    individual_mosquitoes = TRUE,
    human_population = 50,
    total_M = 200,
    init_foim = 0,
    progress_bar = FALSE
  ))
  base_params <- malariasimulationGD::parameterise_total_M(base_params, base_params$total_M)

  timesteps <- 20

  set.seed(123)
  baseline <- malariasimulationGD::run_resumable_simulation(timesteps, parameters = base_params)
  baseline_df <- baseline$data
  baseline_adults <- baseline_df$Sm_gamb_count + baseline_df$Pm_gamb_count + baseline_df$Im_gamb_count

  cube1 <- MGDrivE::cubeMendelian(gtype = c("AA"))
  params_cube1 <- base_params
  params_cube1$cube <- cube1

  set.seed(123)
  run_cube1 <- malariasimulationGD::run_resumable_simulation(timesteps, parameters = params_cube1)
  cube1_df <- run_cube1$data

  expect_equal(cube1_df$E_gamb_count, baseline_df$E_gamb_count, tolerance = 0)
  expect_equal(cube1_df$L_gamb_count, baseline_df$L_gamb_count, tolerance = 0)
  expect_equal(cube1_df$P_gamb_count, baseline_df$P_gamb_count, tolerance = 0)
  expect_equal(
    cube1_df$Sm_gamb_count + cube1_df$Pm_gamb_count + cube1_df$Im_gamb_count,
    baseline_adults,
    tolerance = 0
  )

  cube3 <- MGDrivE::cubeMendelian(gtype = c("AA", "Aa", "aa"))
  params_cube3 <- base_params
  params_cube3$cube <- cube3

  set.seed(123)
  run_cube3 <- malariasimulationGD::run_resumable_simulation(timesteps, parameters = params_cube3)
  geno <- run_cube3$mosquito_genotypes
  cube3_df <- run_cube3$data

  expect_false(is.null(geno))
  expect_identical(colnames(geno$female), cube3$genotypesID)
  expect_identical(colnames(geno$male), cube3$genotypesID)
  expect_equal(geno$V[, 1], rep(1, nrow(geno$V)), tolerance = 0)

  female_total <- rowSums(geno$female)
  male_total <- rowSums(geno$male)
  expect_equal(female_total + male_total, geno$total_adults, tolerance = 0)
  expect_equal(
    female_total,
    cube3_df$Sm_gamb_count + cube3_df$Pm_gamb_count + cube3_df$Im_gamb_count,
    tolerance = 0
  )
})
