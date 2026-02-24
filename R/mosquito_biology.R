#' @title Calculate equilibrium solution for vector counts
#' @description taken from
#' "Modelling the impact of vector control interventions on Anopheles gambiae
#' population dynamics"
#' @param parameters model parameters
#' @param species the index of the species to find the equilibrium for
#' @param foim equilibrium foim
#' @param m the total number of female adult mosquitos
#' @noRd
initial_mosquito_counts <- function(parameters, species, foim, m) {
  omega <- calculate_omega(parameters, species)
  mum <- parameters$mum[[species]]
  n_E <- 2 * omega * mum * parameters$dl * (
    1. + parameters$dpl * parameters$mup
  ) * m

  n_L <- 2 * mum * parameters$dl * (
    1. + parameters$dpl * parameters$mup
  ) * m

  n_P <- 2 * parameters$dpl * mum * m

  n_Sm <- m * mum / (foim + mum)

  incubation_survival <- exp(-mum * parameters$dem)

  n_Pm <- m * foim / (foim + mum) * (
    1. - incubation_survival
  )

  n_Im <- m * foim / (foim + mum) * incubation_survival

  c(n_E, n_L, n_P, n_Sm, n_Pm, n_Im)
}

#' @title Calculate omega value
#' @description useful value for calculating equilibrium solutions for vectors
#' taken from
#' "Modelling the impact of vector control interventions on Anopheles gambiae
#' population dynamics"
#' @param parameters model parameters
#' @param species the index of the species to calculate for
#' @noRd
calculate_omega <- function(parameters, species) {
  sub_omega <- parameters$gamma * parameters$ml / parameters$me - (
    parameters$del / parameters$dl
  ) + (
    (parameters$gamma - 1) * parameters$ml * parameters$del
  )

  mum <- parameters$mum[[species]]

  beta <- eggs_laid(
    parameters$beta,
    mum,
    parameters$blood_meal_rates[[species]]
  )

  -.5 * sub_omega + sqrt(
    .25 * sub_omega**2 +
      .5 * parameters$gamma * beta * parameters$ml * parameters$del /
      (parameters$me * mum * parameters$dl * (
        1. + parameters$dpl * parameters$mup
      ))
  )
}

#' @title Calculate the vector carrying capacity
#' @description taken from
#' "Modelling the impact of vector control interventions on Anopheles gambiae
#' population dynamics"
#' @param parameters model parameters
#' @param m number of adult mosquitoes
#' @param species index of the species to calculate for
calculate_carrying_capacity <- function(parameters, m, species) {
  omega <- calculate_omega(parameters, species)

  m * 2 * parameters$dl * parameters$mum[[species]] * (
    1. + parameters$dpl * parameters$mup
  ) * parameters$gamma * (omega + 1) / (
    omega / (parameters$ml * parameters$del) - (
      1. / (parameters$ml * parameters$dl)
    ) - 1.
  )
}

#' @title Calculate the mean rainfall throughout the year
#' @param parameters model parameters
#' @noRd
calculate_R_bar <- function(parameters) {
  mean(vnapply(1:365, function(t) rainfall(
		t,
    parameters$g0,
    parameters$g,
    parameters$h,
    parameters$rainfall_floor
	)))
}

#' @title Calculate equilibrium total_M from parameters
#'
#' @param parameters to work from
#' @param EIR equilibrium to use, bites per person per year
#' @importFrom stats weighted.mean
#' @noRd
equilibrium_total_M <- function(parameters, EIR) {
  if (EIR == 0) {
    return(0)
  }
  if (parameters$init_foim == 0) {
    stop('init_foim must be > 0 to calculate a non-zero equilibrium total_M')
  }
  mum <- weighted.mean(parameters$mum, parameters$species_proportions)
  total_daily_eir <- EIR * parameters$human_population / 365
  lifetime <- parameters$init_foim * exp(-mum * parameters$dem) / (
    parameters$init_foim + mum
  )
  total_daily_eir / sum(
    parameters$species_proportions * parameters$blood_meal_rates * parameters$Q0 * lifetime
  )
}

#' @title Calculate the yearly offset (in timesteps) for the peak mosquito
#' season
#'
#' @param parameters to work from
#' @export
peak_season_offset <- function(parameters) {
  if (!parameters$model_seasonality) {
    return(0)
  }
  which.max(vnapply(seq(365), function(t) {
    rainfall(
      t,
      parameters$g0,
      parameters$g,
      parameters$h,
      parameters$rainfall_floor
    )
  }))[[1]]
}

#' @title Calculate the death rate of mosquitoes given interventions
#'
#' @param f the feeding rate for this species of mosquito
#' @param W the mean probability that a mosquito feeds and survives
#' @param Z the mean probability that a mosquito is repelled
#' @param Z the mean probability that a mosquito is repelled
#' @noRd
death_rate <- function(f, W, Z, species, parameters) {
  mum <- parameters$mum[[species]]
  p1_0 <- exp(-mum * parameters$foraging_time[[species]])
  gonotrophic_cycle <- get_gonotrophic_cycle(species, parameters)
  p2 <- exp(-mum * gonotrophic_cycle)
  p1 <- p1_0 * W / (1 - Z * p1_0)
  -f * log(p1 * p2)
}

get_gonotrophic_cycle <- function(v, parameters) {
  f <- parameters$blood_meal_rates[[v]]
  gonotrophic_cycle <- 1 / f - parameters$foraging_time[[v]]
}

#' @title Update the individual mosquito model after biting
#'
#' @param variables a list of variables in this simulation
#' @param foim force of infection towards mosquitoes
#' @param events events in the simulation
#' @param species the index of the species to calculate for
#' @param susceptible_species the indices of susceptible mosquitos of the
#' species
#' @param adult_species the indices of adult mosquitos of the species
#' @param mu the death rate of the current species
#' @param parameters the model parameters
#' @param timestep the current timestep
#' @noRd
biting_effects_individual <- function(
    variables,
    foim,
    events,
    species,
    susceptible_species,
    adult_species,
    mu,
    parameters,
    timestep
  ) {
  # deal with mosquito infections
  target <- sample_bitset(susceptible_species, foim)
  variables$mosquito_state$queue_update('Pm', target)
  events$mosquito_infection$schedule(
    target,
    log_uniform(target$size(), parameters$dem)
  )

  # deal with mosquito deaths
  died <- sample_bitset(adult_species, mu)

  events$mosquito_death$schedule(died, 0)
}

cube_genotype_info <- function(cube) {
  if (is.null(cube)) {
    return(list(
      G = 1L,
      genotypesID = "WT",
      wild_type_index = 1L
    ))
  }

  if (is.null(cube$ih)) {
    stop("cube$ih must be provided")
  }
  ih_dim <- dim(cube$ih)
  if (length(ih_dim) != 3) {
    stop("cube$ih must be a 3D array")
  }
  G <- ih_dim[[3]]
  genotypes_id <- cube$genotypesID
  if (is.null(genotypes_id)) {
    genotypes_id <- as.character(seq_len(G))
  }
  if (length(genotypes_id) != G) {
    stop("length(cube$genotypesID) must match dim(cube$ih)[3]")
  }

  list(
    G = as.integer(G),
    genotypesID = genotypes_id,
    wild_type_index = cube_wild_type_index(cube)
  )
}

cube_wild_type_index <- function(cube) {
  if (is.null(cube) || is.null(cube$wildType)) {
    return(1L)
  }
  wt <- cube$wildType
  if (is.numeric(wt) && length(wt) >= 1) {
    return(as.integer(wt[[1]]))
  }
  if (is.character(wt) && length(wt) >= 1) {
    if (is.null(cube$genotypesID)) {
      stop("cube$genotypesID is required when cube$wildType is a character")
    }
    idx <- match(wt[[1]], cube$genotypesID)
    if (is.na(idx)) {
      stop("cube$wildType was not found in cube$genotypesID")
    }
    return(as.integer(idx))
  }
  1L
}

sample_genotype_counts <- function(n, p) {
  if (n <= 0) {
    return(integer(length(p)))
  }
  if (length(p) == 1L) {
    # exact and RNG-free in the trivial one-genotype case
    out <- rep.int(0L, length(p))
    out[[which.max(p)]] <- as.integer(n)
    return(out)
  }
  as.integer(stats::rmultinom(1, size = n, prob = p)[, 1])
}

sample_genotype_ids <- function(n, p) {
  counts <- sample_genotype_counts(n, p)
  rep.int(seq_along(counts), counts)
}

adult_female_genotype_counts_by_species <- function(variables, species_name, G) {
  adult_index <- variables$mosquito_state$get_index_of('NonExistent')$not(TRUE)
  species_index <- variables$species$get_index_of(species_name)$and(adult_index)
  if (species_index$size() == 0) {
    return(rep.int(0, G))
  }
  tabulate(variables$geno_id$get_values(species_index), nbins = G)
}

#' @title Calculate offspring genotype proportions and viability from cube
#' @param cube MGDrivE-style inheritance cube
#' @param female_counts adult female counts by genotype
#' @param male_counts adult male counts by genotype
#' @noRd
calc_pg_V_from_cube <- function(cube, female_counts, male_counts) {
  G <- length(female_counts)
  if (length(male_counts) != G) {
    stop("female_counts and male_counts must have the same length")
  }
  wt <- cube_wild_type_index(cube)
  p_fallback <- rep.int(0, G)
  p_fallback[[wt]] <- 1
  B_zero <- rep.int(0, G)

  if (is.null(cube) || is.null(cube$ih)) {
    return(list(p = p_fallback, V = 1, B = B_zero))
  }

  ih <- cube$ih
  ih_dim <- dim(ih)
  if (length(ih_dim) != 3 || any(ih_dim != c(G, G, G))) {
    stop("cube$ih dimensions must be G x G x G")
  }
  tau <- cube$tau
  if (is.null(tau)) {
    tau <- array(1, dim = ih_dim)
  } else if (!all(dim(tau) == ih_dim)) {
    stop("cube$tau dimensions must match cube$ih")
  }
  eta <- cube$eta
  if (is.null(eta)) {
    eta <- matrix(1, nrow = G, ncol = G)
  }
  if (!all(dim(eta) == c(G, G))) {
    stop("cube$eta dimensions must be G x G")
  }

  total_males <- sum(male_counts)
  if (total_males <= 0) {
    return(list(p = p_fallback, V = 1, B = B_zero))
  }

  Q <- outer(female_counts, male_counts / total_males)
  B <- vnapply(seq_len(G), function(g) {
    sum(Q * ih[, , g] * tau[, , g] * eta)
  })
  total_B <- sum(B)
  if (total_B > 0) {
    p <- B / total_B
  } else {
    p <- p_fallback
  }

  den <- sum(Q * apply(ih, c(1, 2), sum) * eta)
  V <- if (den > 0) total_B / den else 1

  list(p = p, V = V, B = B)
}

#' @title Mosquito emergence process
#' @description Move mosquitos from NonExistent to Sm in line with the number of
#' pupals in the ODE models
#'
#' @param solvers a list of solver objects for each species of mosquito
#' @param models mosquito model objects (used for genotype state when cube is set)
#' @param variables simulation variables (mosquito_state, species, geno_id)
#' @param parameters model parameters
#' @noRd
create_mosquito_emergence_process <- function(
  solvers,
  models,
  variables,
  parameters
  ) {
  state <- variables$mosquito_state
  species <- variables$species
  geno_id <- variables$geno_id
  species_names <- parameters$species
  rate <- .5 * 1 / parameters$dpl
  function(timestep) {
    p_counts <- vnapply(
      solvers,
      function(solver) {
        solver$get_states()[[ODE_INDICES[['P']]]]
      }
    )
    n <- sum(p_counts) * rate
    available <- state$get_size_of('NonExistent')
    if (n > available) {
      stop(paste0(
        'Not enough mosquitoes (short by ',
        n - available,
        '). Please raise parameters$mosquito_limit. ',
        'If you have used parameterise_mosquito_equilibrium,',
        'your seasonality parameters lead to more mosquitoes than expected.'
      ))
    }
    non_existent <- state$get_index_of('NonExistent')
    latest <- 1
    for (i in seq_along(species_names)) {
      to_hatch <- p_counts[[i]] * rate
      hatched <- bitset_at(non_existent, seq(latest, latest + to_hatch))
      n_hatched <- hatched$size()

      if (!is.null(models[[i]]$cube) && n_hatched > 0) {
        cube_info <- cube_genotype_info(models[[i]]$cube)
        female_counts <- adult_female_genotype_counts_by_species(
          variables,
          species_names[[i]],
          cube_info$G
        )
        male_counts <- models[[i]]$genotype_state$male_counts
        pgv <- calc_pg_V_from_cube(models[[i]]$cube, female_counts, male_counts)
        geno_id$queue_update(sample_genotype_ids(n_hatched, pgv$p), hatched)
        models[[i]]$genotype_state$male_counts <- male_counts + sample_genotype_counts(n_hatched, pgv$p)
      } else if (n_hatched > 0) {
        geno_id$queue_update(1L, hatched)
      }

      state$queue_update('Sm', hatched)
      species$queue_update(species_names[[i]], hatched)
      latest <- latest + to_hatch + 1
    }
  }
}
