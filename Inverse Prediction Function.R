
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

inverse_predict_with_prior_mv_custom <- function(
  obs,
  fit_mv,
  species = c("sp1", "sp2", "sp3", "sp4", "sp5"),
  x_var = "x_obs",
  y_var = "y_obs",
  model_x = "x_true",
  model_y = "y_true",
  grid_radius_x = 4,
  grid_radius_y = 4,
  grid_size_x = 21,
  grid_size_y = 21,
  prior_sd_x = 2.0,
  prior_sd_y = 2.0,
  use_normal_prior = TRUE,
  use_grid_constraint = FALSE,
  lik_fun = NULL  #user-defined likelihood function, same as fitted SDM
) {
  # Observed location
  x_obs <- obs[[x_var]]
  y_obs <- obs[[y_var]]
  
  # Build grid
  x_seq <- seq(x_obs - grid_radius_x, x_obs + grid_radius_x, length.out = grid_size_x)
  y_seq <- seq(y_obs - grid_radius_y, y_obs + grid_radius_y, length.out = grid_size_y)
  grid <- expand.grid(x = x_seq, y = y_seq)
  grid$monthf <- obs$monthf
  grid$log_area_km2 <- obs$log_area_km2  #*** for OPT example only
  grid$log_lik <- 0
  
  # Rename for prediction
  newdata <- grid %>%
    rename(!!model_x := x, !!model_y := y)
  
  # Predict species means
  #preds <- posterior_predict(fit_mv, newdata = newdata, ndraws = 50)
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
  grid$prob <- exp(grid$log_posterior - max(grid$log_posterior, na.rm = TRUE))  #the Log-Sum-Exp trick for stability
  grid$prob[is.na(grid$prob)] <- 0
  grid$prob <- grid$prob / sum(grid$prob)
  
  best_idx <- which.max(grid$log_posterior)
  
  list(map = grid[best_idx, c("x", "y")],
       posterior = grid)
}

