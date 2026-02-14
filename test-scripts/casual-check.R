# Load the requisite packages:
# library(malariasimulation)
# Set colour palette:
cols <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7","#F0E442", "#0072B2", "#D55E00")
# Set the timespan over which to simulate
year <- 365; years <- 5; sim_length <- year * years

# Set an initial human population and initial entomological inoculation rate (EIR)
human_population <- 1000
starting_EIR <- 5

# Set the age ranges (in days)
age_min <- seq(0, 80, 5) * 365
age_max <- seq(5, 85, 5) * 365

# Use the get_parameters() function to establish the default simulation parameters, specifying
# age categories from 0-85 in 5 year intervals.
simparams <- get_parameters(
  list(
    human_population = human_population,
    age_group_rendering_min_ages = age_min,
    age_group_rendering_max_ages = age_max
  )
)

# Use set_equilibrium to tune the human and mosquito populations to those required for the
# defined EIR
simparams <- set_equilibrium(simparams, starting_EIR)
# Copy the simulation parameters as demography parameters:
dem_params <- simparams

# We can set our own custom demography:
ages <- round(c(0.083333, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45,
                50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 200) * year)

# Set deathrates for each age group (divide annual values by 365:
deathrates <- c(
  .4014834, .0583379, .0380348, .0395061, .0347255, .0240849, .0300902,
  .0357914, .0443123, .0604932, .0466799, .0426199, .0268332, .049361,
  .0234852, .0988317, .046755, .1638875, .1148753, .3409079, .2239224,
  .8338688) / 365

# Update the parameter list with the custom ages and death rates through the in-built
# set_demography() function, instructing the custom demography to be implemented at the
# beginning of the model run (timestep = 0):
dem_params <- set_demography(
  dem_params,
  agegroups = ages,
  timesteps = 0,
  deathrates = matrix(deathrates, nrow = 1)
)

# Confirm that the custom demography has been set:
dem_params$custom_demography
# Run the simulation with the default demographic set-up:
exp_output <- run_simulation(sim_length, simparams)
exp_output$run <- 'exponential'

# Run the simulation for the custom demographic set-up:
custom_output <- run_simulation(sim_length, dem_params)
custom_output$run <- 'custom'
# Combine the two dataframes:
dem_output <- rbind(exp_output, custom_output)

# Subset the final day of the simulation for each of the two demography runs:
dem_output <- dem_output[dem_output$timestep == 5 * 365,]

# Extract the age variables and convert the dataframe to long format:
convert_to_long <- function(age_min, age_max, output) {
  output <- lapply(
    seq_along(age_min),
    function(i) {
      data.frame(
        age_lower = age_min[[i]],
        age_upper = age_max[[i]],
        n = output[,paste0('n_age_', age_min[[i]], '_',age_max[[i]])],
        age_plot = age_min[[i]]/365,
        run = output$run,
        timestep = output$timestep)
    }
  )
  output <- do.call("rbind", output)
}

# Convert the output for plotting:
dem_output <- convert_to_long(age_min, age_max, dem_output)

# Define the plotting window
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))

# a) Default/Exponentially-distributed demography
plot.new(); grid(lty = 2, col = "grey80", lwd = 0.5, ny = 5, nx = 6); par(new = TRUE)
barplot(height = dem_output[dem_output$run == "exponential", c("n", "age_plot")]$n,
        names = c(paste0(seq(0,75, by = 5),"-",seq(0,75, by = 5)+4), "80+"),
        axes = TRUE, space = 0, ylim = c(0, 250), xaxs = "i", yaxs = "i",
        main = "Default", xlab = "Age Group", ylab = "Individuals",
        cex.axis = 0.8, cex.names = 0.8, cex.lab = 1, cex.main = 1, las = 2,
        col = cols[2]); box()

# b) Custom demography
plot.new()
grid(lty = 2, col = "grey80", lwd = 0.5, ny = 5, nx = 6)
par(new = TRUE)
barplot(height = dem_output[dem_output$run == "custom", c("n", "age_plot")]$n,
        names = c(paste0(seq(0,75, by = 5),"-",seq(0,75, by = 5)+4), "80+"),
        axes = TRUE, space = 0, ylim = c(0, 250), xaxs = "i", yaxs = "i",
        main = "Custom", xlab = "Age Group",
        cex.axis = 0.8, cex.names = 0.8, cex.lab = 1, cex.main = 1, las = 2,
        col = cols[1]); box()

