#' @title Biting process
#' @description
#' This is the biting process. It results in human and mosquito infection and
#' mosquito death.
#' @param renderer the model renderer object
#' @param solvers mosquito ode solvers
#' @param models mosquito ode models
#' @param variables a list of all of the model variables
#' @param events a list of all of the model events
#' @param parameters model pararmeters
#' @param lagged_infectivity a list of LaggedValue objects with historical sums
#' of infectivity, one for every metapopulation
#' @param lagged_eir a LaggedValue class with historical EIRs
#' @param mixing_fn a function to retrieve the mixed EIR and infectivity based
#' on the other populations
#' @param mixing_index an index for this population's position in the
#' lagged_infectivity list (default: 1)
#' @param infection_outcome competing hazards object for infection rates
#' @param timestep the current timestep
#' @noRd
create_biting_process <- function(
  renderer,
  solvers,
  models,
  variables,
  events,
  parameters,
  lagged_infectivity,
  lagged_eir,
  mixing_fn = NULL,
  mixing_index = 1,
  infection_outcome
  ) {
  function(timestep) {
    # Calculate combined EIR
    age <- get_age(variables$birth$get_values(), timestep)
    bitten <- simulate_bites(
      renderer,
      solvers,
      models,
      variables,
      events,
      age,
      parameters,
      timestep,
      lagged_infectivity,
      lagged_eir,
      mixing_fn,
      mixing_index
    )
    
    simulate_infection(
      variables,
      events,
      bitten$bitten_humans,
      bitten$n_bites_per_person,
      age,
      parameters,
      timestep,
      renderer,
      infection_outcome
    )
  }
}

#' @importFrom stats rpois
simulate_bites <- function(
  renderer,
  solvers,
  models,
  variables,
  events,
  age,
  parameters,
  timestep,
  lagged_infectivity,
  lagged_eir,
  mixing_fn = NULL,
  mixing_index = 1
  ) {
  bitten_humans <- individual::Bitset$new(parameters$human_population)
  n_bites_per_person <- numeric(0)
  vector_infectivity_active <- !is.null(parameters$vector_infectivity_g_by_species)
  
  human_infectivity <- variables$infectivity$get_values()
  if (parameters$tbv) {
    human_infectivity <- account_for_tbv(
      timestep,
      human_infectivity,
      variables,
      parameters
    )
  }
  renderer$render('infectivity', mean(human_infectivity), timestep)
  
  # Calculate pi (the relative biting rate for each human)
  psi <- unique_biting_rate(age, parameters)
  zeta <- variables$zeta$get_values()
  .pi <- human_pi(zeta, psi)
  
  # Get some indices for later
  if (parameters$individual_mosquitoes) {
    infectious_index <- variables$mosquito_state$get_index_of('Im')
    susceptible_index <- variables$mosquito_state$get_index_of('Sm')
    adult_index <- variables$mosquito_state$get_index_of('NonExistent')$not(TRUE)
    genotype_tracking <- !is.null(parameters$cube) && !is.null(variables$geno_id)
    if (vector_infectivity_active && !genotype_tracking) {
      stop("Internal error: genotype-specific vector_infectivity_g requires genotype tracking state")
    }
    if (genotype_tracking) {
      cube_info <- cube_genotype_info(parameters$cube)
      female_geno_totals <- rep.int(0, cube_info$G)
      male_geno_totals <- rep.int(0, cube_info$G)
      mu_by_species <- rep(0, length(parameters$species))
      female_totals_by_species <- rep(0L, length(parameters$species))
      V_by_species <- rep(1, length(parameters$species))
      infectious_geno_totals <- rep.int(0, cube_info$G)
    }
  }
  if (parameters$individual_mosquitoes && genotype_tracking && genotype_debug_enabled(parameters, timestep)) {
    for (s_i in seq_along(parameters$species)) {
      genotype_debug_log_counts(
        parameters,
        timestep,
        "BITE_START",
        parameters$species[[s_i]],
        genotype_debug_species_counts(variables, models, parameters, s_i),
        extra = "visible state at start of biting (after release process; emergence females may still be queued)"
      )
    }
  }
  
  EIR <- 0
  
  for (s_i in seq_along(parameters$species)) {
    species_name <- parameters$species[[s_i]]
    solver_states <- solvers[[s_i]]$get_states()
    p_bitten <- prob_bitten(timestep, variables, s_i, parameters)
    Q0 <- parameters$Q0[[s_i]]
    W <- average_p_successful(p_bitten$prob_bitten_survives, .pi, Q0)
    Z <- average_p_repelled(p_bitten$prob_repelled, .pi, Q0)
    f <- blood_meal_rate(s_i, Z, parameters)
    a <- .human_blood_meal_rate(f, s_i, W, parameters)
    lambda <- effective_biting_rates(a, .pi, p_bitten)

    if (parameters$individual_mosquitoes) {
      species_index <- variables$species$get_index_of(
        parameters$species[[s_i]]
      )$and(adult_index)
      vector_infectivity_weights <- NULL
      if (vector_infectivity_active) {
        vector_infectivity_weights <- vector_infectivity_g_weights_for_species(
          parameters,
          species_name
        )
      }
      if (is.null(vector_infectivity_weights)) {
        n_infectious <- calculate_infectious_individual(
          s_i,
          variables,
          infectious_index,
          adult_index,
          species_index,
          parameters
        )
      } else {
        infectious_counts_g <- calculate_infectious_individual_genotype_counts(
          variables,
          infectious_index,
          species_index,
          cube_info$G
        )
        n_infectious_total <- sum(infectious_counts_g)
        n_infectious <- sum(infectious_counts_g * vector_infectivity_weights)
        infectious_geno_totals <- infectious_geno_totals + infectious_counts_g
        renderer$render(
          paste0("infectivity_weighted_I_", species_name),
          n_infectious,
          timestep
        )
        renderer$render(
          paste0("vector_infectivity_mean_", species_name),
          if (n_infectious_total > 0) n_infectious / n_infectious_total else NA_real_,
          timestep
        )
      }
    } else {
      n_infectious <- calculate_infectious_compartmental(solver_states)
    }
    
    # store the current population's EIR for later
    lagged_eir[[s_i]]$save(
      n_infectious * a,
      timestep
    )

    # lagged EIR
    if (is.null(mixing_fn)) {
      species_eir <- lagged_eir[[s_i]]$get(timestep - parameters$de)
    } else {
      species_eir <- mixing_fn(timestep=timestep)$eir[[mixing_index, s_i]]
    }

    renderer$render(paste0('EIR_', species_name), species_eir, timestep)
    EIR <- EIR + species_eir
    if(parameters$parasite == "falciparum"){
      # p.f model factors eir by psi
      expected_bites <- species_eir * mean(psi)
    } else if (parameters$parasite == "vivax"){
      # p.v model standardises biting rate het to eir
      expected_bites <- species_eir
    }
    
    if (expected_bites > 0) {
      n_bites <- rpois(1, expected_bites)
      if (n_bites > 0) {
        bitten <- fast_weighted_sample(n_bites, lambda)
        bitten_humans$insert(bitten)
        renderer$render('n_bitten', bitten_humans$size(), timestep)
        if(parameters$parasite == "vivax"){
          # p.v must pass through the number of bites per person
          n_bites_per_person <- tabulate(bitten, nbins = length(lambda))
        }
      }
    }

    lagged_infectivity$save(sum(human_infectivity * .pi), timestep)

    if (is.null(mixing_fn)) {
      infectivity <- lagged_infectivity$get(timestep - parameters$delay_gam)
    } else {
      infectivity <- mixing_fn(timestep=timestep)$inf[[mixing_index]]
    }

    foim <- calculate_foim(a, infectivity)
    renderer$render(paste0('FOIM_', species_name), foim, timestep)
    mu <- death_rate(f, W, Z, s_i, parameters)
    renderer$render(paste0('mu_', species_name), mu, timestep)
    
    if (parameters$individual_mosquitoes) {
      # update the ODE with stats for ovoposition calculations
      effective_total_M <- species_index$size()
      if (genotype_tracking && !is.null(models[[s_i]]$genotype_state)) {
        female_counts <- tabulate(
          variables$geno_id$get_values(species_index),
          nbins = cube_info$G
        )
        male_counts <- models[[s_i]]$genotype_state$male_counts
        pgv <- calc_pg_V_from_cube(models[[s_i]]$cube, female_counts, male_counts)
        V_by_species[[s_i]] <- pgv$V
        effective_total_M <- effective_total_M * pgv$V
        aquatic_mosquito_model_set_egg_proportions(models[[s_i]]$.model, as.numeric(pgv$p))
        female_geno_totals <- female_geno_totals + female_counts
        male_geno_totals <- male_geno_totals + male_counts
        mu_by_species[[s_i]] <- mu
        female_totals_by_species[[s_i]] <- species_index$size()
        models[[s_i]]$genotype_state$last_V <- pgv$V
      }

      aquatic_mosquito_model_update(
        models[[s_i]]$.model,
        effective_total_M,
        f,
        mu
      )
      
      # update the individual mosquitoes
      susceptible_species_index <- susceptible_index$copy()$and(species_index)
      
      biting_effects_individual(
        variables,
        foim,
        events,
        s_i,
        susceptible_species_index,
        species_index,
        mu,
        parameters,
        timestep
      )
    } else {
      adult_mosquito_model_update(
        models[[s_i]]$.model,
        mu,
        foim,
        solver_states[[ADULT_ODE_INDICES['Sm']]],
        f
      )
    }
  }

  if (parameters$individual_mosquitoes && vector_infectivity_active &&
      !is.null(parameters$mosquito_infectious_genotype_history)) {
    parameters$mosquito_infectious_genotype_history[timestep, ] <- infectious_geno_totals
  }

  if (parameters$individual_mosquitoes && genotype_tracking) {
    history <- parameters$mosquito_genotype_history
    if (!is.null(history)) {
      history$female[timestep, ] <- female_geno_totals
      history$male[timestep, ] <- male_geno_totals
      history$V[timestep, ] <- V_by_species
      history$total_adults[timestep] <- sum(female_geno_totals + male_geno_totals)
      if (genotype_debug_enabled(parameters, timestep)) {
        for (s_i in seq_along(parameters$species)) {
          genotype_debug_log(
            parameters,
            timestep,
            "HWRITE",
            parameters$species[[s_i]],
            sprintf(
              "row=%d writes visible counts before male death update: F{%s} M{%s}",
              timestep,
              genotype_debug_fmt_counts(female_geno_totals, names(parameters$mosquito_genotype_history$female[timestep, ])),
              genotype_debug_fmt_counts(male_geno_totals, names(parameters$mosquito_genotype_history$male[timestep, ]))
            )
          )
        }
      }
    }
    for (s_i in seq_along(parameters$species)) {
      if (is.null(models[[s_i]]$genotype_state)) {
        next
      }
      male_counts <- models[[s_i]]$genotype_state$male_counts
      if (length(male_counts) == 1L) {
        # Keep the trivial cube case RNG-free and aligned with the implicit 1:1 sex ratio.
        male_counts[[1]] <- female_totals_by_species[[s_i]]
      } else {
        death_prob <- min(mu_by_species[[s_i]], 1)
        male_counts <- stats::rbinom(length(male_counts), size = male_counts, prob = 1 - death_prob)
      }
      models[[s_i]]$genotype_state$male_counts <- male_counts
      if (genotype_debug_enabled(parameters, timestep)) {
        genotype_debug_log_counts(
          parameters,
          timestep,
          "AFTER_MALE_DEATH",
          parameters$species[[s_i]],
          genotype_debug_species_counts(variables, models, parameters, s_i),
          extra = "female deaths are scheduled via mosquito_death event (delay=0)"
        )
      }
    }
  }

  list(bitten_humans = bitten_humans, n_bites_per_person = n_bites_per_person)
}


# =================
# Utility functions
# =================

calculate_eir <- function(species, solvers, variables, parameters, timestep) {
  a <- human_blood_meal_rate(species, variables, parameters, timestep)
  infectious <- calculate_infectious(species, solvers, variables, parameters)
  infectious * a
}

effective_biting_rates <- function(a, .pi, p_bitten) {
  a * .pi * p_bitten$prob_bitten / sum(.pi * p_bitten$prob_bitten_survives)
}

calculate_infectious <- function(species, solvers, variables, parameters) {
  if (parameters$individual_mosquitoes) {
    adult_index <- variables$mosquito_state$get_index_of('NonExistent')$not(TRUE)
    species_name <- parameters$species[[species]]
    species_index <- variables$species$get_index_of(
      species_name
    )$and(adult_index)
    vector_infectivity_weights <- vector_infectivity_g_weights_for_species(
      parameters,
      species_name
    )
    if (!is.null(vector_infectivity_weights)) {
      if (is.null(parameters$cube) || is.null(variables$geno_id)) {
        stop("Internal error: genotype-specific vector_infectivity_g requires genotype tracking state")
      }
      cube_info <- cube_genotype_info(parameters$cube)
      counts_g <- calculate_infectious_individual_genotype_counts(
        variables,
        variables$mosquito_state$get_index_of('Im'),
        species_index,
        cube_info$G
      )
      return(sum(counts_g * vector_infectivity_weights))
    }
    return(
      calculate_infectious_individual(
        species,
        variables,
        variables$mosquito_state$get_index_of('Im'),
        adult_index,
        species_index,
        parameters
      )
    )
  }
  calculate_infectious_compartmental(solvers[[species]]$get_states())
}

calculate_infectious_individual_genotype_counts <- function(
  variables,
  infectious_index,
  species_index,
  G
  ) {
  infectious_species_index <- infectious_index$copy()$and(species_index)
  if (infectious_species_index$size() == 0) {
    return(rep.int(0, G))
  }
  tabulate(variables$geno_id$get_values(infectious_species_index), nbins = G)
}

calculate_infectious_individual <- function(
  species,
  variables,
  infectious_index,
  adult_index,
  species_index,
  parameters
  ) {
  infectious_index$copy()$and(species_index)$size()
}

calculate_infectious_compartmental <- function(solver_states) {
  max(solver_states[[ADULT_ODE_INDICES['Im']]], 0)
}

intervention_coefficient <- function(p_bitten) {
  p_bitten$prob_bitten / sum(p_bitten$prob_bitten_survives)
}

human_pi <- function(zeta, psi) {
  (zeta * psi) / sum(zeta * psi)
}

blood_meal_rate <- function(v, z, parameters) {
  gonotrophic_cycle <- get_gonotrophic_cycle(v, parameters)
  interrupted_foraging_time <- parameters$foraging_time[[v]] / (1 - z)
  1 / (interrupted_foraging_time + gonotrophic_cycle)
}

human_blood_meal_rate <- function(species, variables, parameters, timestep) {
  age <- get_age(variables$birth$get_values(), timestep)
  psi <- unique_biting_rate(age, parameters)
  zeta <- variables$zeta$get_values()
  p_bitten <- prob_bitten(timestep, variables, species, parameters)
  .pi <- human_pi(zeta, psi)
  Q0 <- parameters$Q0[[species]]
  W <- average_p_successful(p_bitten$prob_bitten_survives, .pi, Q0)
  Z <- average_p_repelled(p_bitten$prob_repelled, .pi, Q0)
  f <- blood_meal_rate(species, Z, parameters)
  .human_blood_meal_rate(f, species, W, parameters)
}

.human_blood_meal_rate <- function(f, v, W, parameters) {
  Q <- 1 - (1 - parameters$Q0[[v]]) / W
  Q * f
}

average_p_repelled <- function(p_repelled, .pi, Q0) {
  Q0 * sum(.pi * p_repelled)
}

average_p_successful <- function(prob_bitten_survives, .pi, Q0) {
  (1 - Q0) + Q0 * sum(.pi *  prob_bitten_survives)
}

# Unique biting rate (psi) for a human of a given age
unique_biting_rate <- function(age, parameters) {
  1 - parameters$rho * exp(- age / parameters$a0)
}

#' @title Calculate the force of infection towards mosquitoes
#'
#' @param a human blood meal rate
#' @param infectivity_sum the sum of each individual's infectivity 
#' @noRd
calculate_foim <- function(a, infectivity_sum) {
  a * infectivity_sum
}
