
####################################################################
#####  Catch Location Refinement using Species Compositions    #####
#####                  J.A.Smith NSW DPI, 2025                 #####
#####               james.a.smith@dpi.nsw.gov.au               #####
####################################################################

## This code provides the inverse prediction function
## This also requires a custom likelihood function; two examples are
## also specified below


####################################################
## INVERT THE MODEL TO PREDICT LOCATION FROM SPECIES
####################################################

## This relies on having a 'true' location specified, e.g. it will place a 
## prior around that location to influence the predictions. But if using this
## for location identification, and no inaccurate location exists, add a default
## inaccurate location, but place a wide uninfoirmative prior over that location

inverse_predict_with_prior_mv_custom <- function(
  obs,                           # The inaccurate locations to estimate
  fit_mv,                        # The fitted SDM to use as a likelihood function
  species = c("sp1", "sp2", "sp3", "sp4", "sp5"),  # The species IDs
  x_var = "x_obs",               # Name of x variable for inaccurate data
  y_var = "y_obs",               # Name of y variable for inaccurate data
  model_x = "x_true",            # Name of x variable for accurate data
  model_y = "y_true",            # Name of x variable for accurate data
  grid_radius_x = 4,             # The search radius around the inaccurate observation
  grid_radius_y = 4,
  grid_size_x = 21,              # The resolution (grid cells) of the search grid
  grid_size_y = 21,
  prior_sd_x = 2.0,              # The SD of the prior around inaccurate observation
  prior_sd_y = 2.0,
  use_normal_prior = TRUE,       # A normal prior around the inaccurate observation
  use_grid_constraint = FALSE,   # Use a grid cell to constrain predictions (insetad of or in addition to normal)
  lik_fun = NULL                 # User-defined likelihood function; same as fitted SDM
) {
  
  # Observed location
  x_obs <- obs[[x_var]]
  y_obs <- obs[[y_var]]
  
  # Build search grid
  x_seq <- seq(x_obs - grid_radius_x, x_obs + grid_radius_x, length.out = grid_size_x)
  y_seq <- seq(y_obs - grid_radius_y, y_obs + grid_radius_y, length.out = grid_size_y)
  grid <- expand.grid(x = x_seq, y = y_seq)
  grid$monthf <- obs$monthf
  grid$log_lik <- 0
  
  # Rename for prediction
  newdata <- grid %>%
    rename(!!model_x := x, !!model_y := y)
  
  # Predict species means
  preds <- posterior_epred(fit_mv, newdata = newdata, ndraws = 50)
  mu_hat <- apply(preds, c(2, 3), median)  # rows × species
  
  # Check that lik_fun is provided
  if (is.null(lik_fun)) {
    stop("Provide a species-level log-likelihood function via the `lik_fun` argument.")
  }
  
  # Calculate joint log-likelihood for each grid point
  log_lik_matrix <- matrix(NA, nrow = nrow(grid), ncol = length(species))
  for (j in seq_along(species)) {
    obs_count <- obs[[species[j]]]
    log_lik_matrix[, j] <- lik_fun(obs_count, mu_hat[, j], species[j], fit_mv)
  }
  grid$log_lik <- rowSums(log_lik_matrix)
  
  # Priors
  grid$log_prior <- 0
  if (use_normal_prior && !use_grid_constraint) {
    grid$log_prior <- dnorm(grid$x, mean = x_obs, sd = prior_sd_x, log = TRUE) +
      dnorm(grid$y, mean = y_obs, sd = prior_sd_y, log = TRUE)
  } else if (!use_normal_prior && use_grid_constraint) {
    if (all(c("grid_id_x", "grid_id_y") %in% names(obs))) {
      grid$in_bounds <- with(grid,
                             x >= obs$grid_min_x & x <= obs$grid_max_x &
                               y >= obs$grid_min_y & y <= obs$grid_max_y)
      grid$log_prior <- ifelse(grid$in_bounds, 0, -Inf)
    }
  } else if (use_normal_prior && use_grid_constraint) {
    if (all(c("grid_id_x", "grid_id_y") %in% names(obs))) {
      grid$in_bounds <- with(grid,
                             x >= obs$grid_min_x & x <= obs$grid_max_x &
                               y >= obs$grid_min_y & y <= obs$grid_max_y)
      grid$log_prior <- -Inf
      grid$log_prior[grid$in_bounds] <-
        dnorm(grid$x[grid$in_bounds], mean = x_obs, sd = prior_sd_x, log = TRUE) +
        dnorm(grid$y[grid$in_bounds], mean = y_obs, sd = prior_sd_y, log = TRUE)
    }
  }
  
  # Posterior and normalization
  grid$log_posterior <- grid$log_lik + grid$log_prior
  grid$prob <- exp(grid$log_posterior - max(grid$log_posterior, na.rm = TRUE))  # Log-Sum-Exp trick for stability
  grid$prob[is.na(grid$prob)] <- 0
  grid$prob <- grid$prob / sum(grid$prob)
  
  best_idx <- which.max(grid$log_posterior)
  
  list(map = grid[best_idx, c("x", "y")],
       posterior = grid)
}



###################################
## SOME CUSTOM LIKELIHOOD FUNCTIONS
###################################

## Negative binomial function
lik_fun_nb <- function(y_obs, mu, species_name, fit_mv) {
  shape <- median(as.matrix(fit_mv, variable = paste0("shape_", species_name)))
  mu_hat <- pmax(mu, 1e-3)  # to prevent log(0) and instability
  #mu_hat <- pmin(pmax(mu_hat, 1e-3), 1e5)  #also prevent massive numbers if a problem
  dnbinom(y_obs, mu = mu_hat, size = shape, log = TRUE)
}
# When mu_hat ≈ 0 and y_obs > 0, the log-likelihood is extremely negative;
# just one such species can dominate the joint log-likelihood and collapse the posterior
# into one or two grid cells. Bounding mu prevents this.


## Poisson function
lik_fun_pois <- function(y_obs, mu, species_name, fit_mv) {
  dpois(y_obs, lambda = mu, log = TRUE)
}

# ## Tweedie function (untested)
# lik_fun_tweedie <- function(y_obs, mu, species_name, fit_mv) {
#   require(tweedie)
#   # Extract power and phi from the fitted model
#   p <- median(as.matrix(fit_mv, variable = paste0("power_", species_name)))
#   phi <- median(as.matrix(fit_mv, variable = paste0("phi_", species_name)))
#   dtweedie(y_obs, mu = mu, phi = phi, power = p, log = TRUE)
# }
