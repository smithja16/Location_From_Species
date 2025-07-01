
####################################################################
#####  Catch Location Refinement using Species Compositions    #####
#####                  J.A.Smith NSW DPI, 2025                 #####
#####               james.a.smith@dpi.nsw.gov.au               #####
####################################################################

## This script provides various functions for data processing and visualisation


###################################################################
## EXTRACT ESTIMATED LOCATIONS AND TEST ACCURACY - Hierarchical GLM
###################################################################

summarize_accuracy_lv <- function(fish_data, loc_unknown_samples) {
  
  loc_unknown_means <- apply(loc_unknown_samples, c(2, 3), mean)
  
  ## Add the refined (estimated) locations
  fish_data$x_refined <- fish_data$x_obs
  fish_data$y_refined <- fish_data$y_obs
  for (i in 1:length(unknown_ids)) {
    idx <- unknown_ids[i]
    fish_data$x_refined[idx] <- loc_unknown_means[i, 1]
    fish_data$y_refined[idx] <- loc_unknown_means[i, 2]
  }
  
  ## Evaluate model performance using Euclidean RMSE
  inaccurate <- which(!fish_data$accurate)  # Calculate errors for inaccurate locations only
  obs_dist <- sqrt(
    (fish_data$x_obs[inaccurate] - fish_data$x_true[inaccurate])^2 +
      (fish_data$y_obs[inaccurate] - fish_data$y_true[inaccurate])^2 )  #Error of original observations
  ref_dist <- sqrt(
    (fish_data$x_refined[inaccurate] - fish_data$x_true[inaccurate])^2 +
      (fish_data$y_refined[inaccurate] - fish_data$y_true[inaccurate])^2 )  # Error of updated observations
  
  ## Performance metrics
  mean_obs_error <- mean(obs_dist)  # Error of original inaccurate locations
  mean_ref_error <- mean(ref_dist)  # Error of updated estimated locations
  mean_improvement <- mean(obs_dist - ref_dist)  # Mean improvement in this error
  percent_improved <- 100 * mean(ref_dist < obs_dist)  # Perc. of observations improved
  
  list(fish_data = fish_data,
       summary = data.frame(mean_obs_error,
                            mean_ref_error,
                            mean_improvement,
                            percent_improved))
}



#####################################################################
## PLOT LOCATION AND ERROR FOR EXAMPLE OBSERVATION - Hierarchical GLM
#####################################################################

plot_prediction_diagnostics_lv <- function(samples, fish_data) {
  samples_df <- as.data.frame(samples)
  colnames(samples_df) <- c("x", "y")
  
  # Compute mean and covariance and marginal 90% credible intervals
  mean_loc <- colMeans(samples_df)
  cov_loc <- cov(samples_df)
  x_ci <- quantile(samples_df$x, probs = c(0.05, 0.95))
  y_ci <- quantile(samples_df$y, probs = c(0.05, 0.95))
  
  # Data frames for ggplot
  mean_df <- data.frame(x = mean_loc[1],
                        y = mean_loc[2])
  obs_df <- data.frame(x = fish_data$x_obs[unknown_ids[i]],
                       y = fish_data$y_obs[unknown_ids[i]])
  true_df <- data.frame(x = fish_data$x_true[unknown_ids[i]],
                        y = fish_data$y_true[unknown_ids[i]])
  
  # Optional ellipse calculations (difficult to compare across methods)
  eig <- eigen(cov_loc)
  angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])
  
  # Plot
  plot <- ggplot(data = mean_df, aes(x = x, y = y)) +
    #geom_point(data = samples_df, aes(x = x, y = y),
               #alpha = 0.05, color = "blue") +  # Posterior samples
    geom_point(data = mean_df, aes(x = x, y = y),
               color = "black", shape = 3, size = 3) +  # Posterior mean
    geom_ellipse(aes(x0 = mean_loc[1], y0 = mean_loc[2],
                     a = sqrt(qchisq(0.90, df = 2) * eig$values[1]),
                     b = sqrt(qchisq(0.90, df = 2) * eig$values[2]),
                     angle = angle),
                 fill = NA, color = "blue", linewidth = 1) +  # 95% ellipse
    annotate("rect", xmin = x_ci[1], xmax = x_ci[2],
             ymin = y_ci[1], ymax = y_ci[2],
             color = "darkgreen", fill = NA, linetype = "dashed") +  # Marginal CI box
    geom_point(data = obs_df, aes(x = x, y = y), color = "red", size = 4) +  # Observed
    geom_point(data = true_df, aes(x = x, y = y), color = "orange", size = 4) +  # True
    xlim(-2, 12) + ylim(-2, 12) +
    coord_fixed() +
    labs(title = paste0("Hierachical GLM, i=", i),
         subtitle = "Ornage = True; Red = Observed; Black = Estimated") +
    theme_bw()
  
  return(plot)
}



################################################################
## EXTRACT ESTIMATED LOCATIONS AND TEST ACCURACY - Inverse Model
################################################################

summarize_accuracy_inv <- function(test_obs, inverse_results) {
  results <- test_obs %>%
    mutate(
      x_pred = sapply(inverse_results, function(r) r$map$x),
      y_pred = sapply(inverse_results, function(r) r$map$y),
      error_obs = sqrt((x_obs - x_true)^2 + (y_obs - y_true)^2),
      error_pred = sqrt((x_pred - x_true)^2 + (y_pred - y_true)^2),
      improvement = error_obs - error_pred )
  
  summary_stats <- results %>%
    summarise(
      mean_error_obs = mean(error_obs),
      mean_error_pred = mean(error_pred),
      mean_improvement = mean(improvement),
      percent_improved = mean(improvement > 0) * 100 )
  
  list(results = results, summary = summary_stats)
}


########################################
## CALCULATE UNCERTAINTY - Inverse Model
########################################

estimate_uncertainty_inv <- function(posterior_df) {
  stopifnot(all(c("x", "y", "prob") %in% names(posterior_df)))
  
  # Marginal means
  x_mean <- weighted.mean(posterior_df$x, posterior_df$prob)
  y_mean <- weighted.mean(posterior_df$y, posterior_df$prob)
  
  # Marginal standard deviations
  x_sd <- sqrt(weighted.mean((posterior_df$x - x_mean)^2, posterior_df$prob))
  y_sd <- sqrt(weighted.mean((posterior_df$y - y_mean)^2, posterior_df$prob))
  
  # Marginal quantiles via grouping
  x_marg <- posterior_df %>%
    group_by(x) %>%
    summarise(prob = sum(prob), .groups = "drop") %>%
    arrange(x) %>%
    mutate(cdf = cumsum(prob))
  
  y_marg <- posterior_df %>%
    group_by(y) %>%
    summarise(prob = sum(prob), .groups = "drop") %>%
    arrange(y) %>%
    mutate(cdf = cumsum(prob))
  
  # Find 5% and 95% quantile bounds (90% credible interval)
  x_lower <- x_marg$x[which.min(abs(x_marg$cdf - 0.05))]
  x_upper <- x_marg$x[which.min(abs(x_marg$cdf - 0.95))]
  
  y_lower <- y_marg$y[which.min(abs(y_marg$cdf - 0.05))]
  y_upper <- y_marg$y[which.min(abs(y_marg$cdf - 0.95))]
  
  list(x_mean = x_mean, x_sd = x_sd, x_lower = x_lower, x_upper = x_upper,
       y_mean = y_mean, y_sd = y_sd, y_lower = y_lower, y_upper = y_upper)
}



##################################################################
## PLOT LOCATION AND ERROR FOR EXAMPLE OBSERVATION - Inverse model
##################################################################

sample_posterior_points <- function(posterior_df, n = 1000) {
  posterior_df[sample(1:nrow(posterior_df), size = n, replace = TRUE, prob = posterior_df$prob),
               c("x", "y")]
}


plot_prediction_diagnostics_inv <- function(
    obs,
    result,
    n_samples = 1000) {
  
  posterior <- result$posterior
  map_point <- result$map
  
  # Extract marginal bounds
  uncertainty <- estimate_uncertainty_inv(posterior_df=posterior)
  x_lower <- uncertainty$x_lower
  x_upper <- uncertainty$x_upper
  y_lower <- uncertainty$y_lower
  y_upper <- uncertainty$y_upper
  
  # Normalize prior and likelihood for plotting
  posterior$prior_prob <- exp(posterior$log_prior - max(posterior$log_prior, na.rm = TRUE))
  posterior$prior_prob <- posterior$prior_prob / sum(posterior$prior_prob, na.rm = TRUE)
  
  posterior$likelihood_prob <- exp(posterior$log_lik - max(posterior$log_lik, na.rm = TRUE))
  posterior$likelihood_prob <- posterior$likelihood_prob / sum(posterior$likelihood_prob, na.rm = TRUE)
  
  # Sample from posterior
  samples <- sample_posterior_points(posterior, n_samples)
  names(samples) <- c("x", "y")
  
  # Ellipse from covariance
  cov_mat <- cov(samples)
  center <- colMeans(samples)
  eig <- eigen(cov_mat)
  angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])
  a <- sqrt(qchisq(0.90, df = 2) * eig$values[1])
  b <- sqrt(qchisq(0.90, df = 2) * eig$values[2])
  ellipse_df <- data.frame(x0 = center[1], y0 = center[2], a = a, b = b, angle = angle)
  
  base <- ggplot(posterior, aes(x = x, y = y))# + coord_fixed()
  
  p1 <- base +
    geom_tile(aes(fill = prior_prob)) +
    annotate("point", x = obs$x_obs, y = obs$y_obs, color = "red", size = 3) +
    #annotate("point", x = map_point$x, y = map_point$y, color = "orange", shape = 17, size = 4) +
    #annotate("point", x = map_point$x, y = map_point$y, color = "white", shape = 3, size = 3) +
    xlim(-2, 12) + ylim(-2, 12) +
    coord_fixed() +
    labs(title = "Inv. model: Prior Surface", fill = "Prior")
  
  p2 <- base +
    geom_tile(aes(fill = likelihood_prob)) +
    annotate("point", x = obs$x_obs, y = obs$y_obs, color = "red", size = 3) +
    annotate("point", x = map_point$x, y = map_point$y, color = "black", shape = 3, size = 3) +
    annotate("point", x = obs$x_true, y = obs$y_true, color = "orange", size = 3) +
    annotate("point", x = map_point$x, y = map_point$y, color = "white", shape = 3, size = 3) +
    xlim(-2, 12) + ylim(-2, 12) +
    coord_fixed() +
    labs(title = "Inv. model: Joint Likelihood Surface", fill = "Likelihood")
  
  p3 <- base +
    geom_tile(aes(fill = prob)) +
    annotate("point", x = obs$x_obs, y = obs$y_obs, color = "red", size = 3) +
    annotate("point", x = obs$x_true, y = obs$y_true, color = "orange", size = 3) +
    annotate("point", x = map_point$x, y = map_point$y, color = "white", shape = 3, size = 3) +
    xlim(-2, 12) + ylim(-2, 12) +
    coord_fixed() +
    labs(title = "Inv. model: Posterior Surface", fill = "Posterior")
  
  p4 <- ggplot(samples, aes(x = x, y = y)) +
    #geom_point(color = "grey40", alpha = 0.2, size = 0.7) +
    geom_ellipse(data = ellipse_df, aes(x0 = x0, y0 = y0, a = a, b = b, angle = angle),
                 inherit.aes = FALSE, color = "blue", linewidth = 1, fill = NA) +
    annotate("rect", xmin = x_lower, xmax = x_upper, ymin = y_lower, ymax = y_upper,
             color = "darkgreen", fill = NA, linewidth = 1, linetype = "dashed") +
    annotate("point", x = obs$x_obs, y = obs$y_obs, color = "red", size = 3) +
    annotate("point", x = obs$x_true, y = obs$y_true, color = "orange", size = 3) +
    annotate("point", x = map_point$x, y = map_point$y, color = "black", shape = 3, size = 3) +
    xlim(-2, 12) + ylim(-2, 12) +
    coord_fixed() +
    labs(title = paste0("Inv. model: Posterior Prediction, i=",i),
         subtitle = "True=orange, Obs.=red, Pred.=black") +
    theme_bw()
  
  list(prior = p1, likelihood = p2, posterior = p3, uncertainty = p4)
}



#######################################################################
## EXTRACT ESTIMATED LOCATIONS AND TEST ACCURACY - Location as Response
#######################################################################

summarize_accuracy_gam_rf <- function(test_obs, predictions) {
  results <- test_obs %>%
    mutate(
      x_pred = predictions[,1],
      y_pred = predictions[,2],
      error_obs = sqrt((x_obs - x_true)^2 + (y_obs - y_true)^2),
      error_pred = sqrt((x_pred - x_true)^2 + (y_pred - y_true)^2),
      improvement = error_obs - error_pred
    )
  
  summary_stats <- results %>%
    summarise(
      mean_error_obs = mean(error_obs),
      mean_error_pred = mean(error_pred),
      mean_improvement = mean(improvement),
      percent_improved = mean(improvement > 0) * 100
    )
  
  list(results = results, summary = summary_stats)
}



#############################################################################
## PLOT LOCATION AND ERROR FOR EXAMPLE OBSERVATION - Location as Response GAM
#############################################################################

plot_prediction_diagnostics_gam <- function(test_obs, model, i) {
  
  # Get predictions with unconditional standard errors (includes smoothing parameter uncertainty)
  pred_unconditional <- predict(model, newdata = test_obs[i,], se.fit = TRUE, unconditional = TRUE)
  
  # Extract residual variance from the model
  sigma_residual <- sqrt(model$sig2)
  #sigma_residual <- 0
  
  # Combine unconditional SE with residual variance for full prediction uncertainty
  pred_se_full <- sqrt(pred_unconditional$se.fit^2 + sigma_residual^2)
  
  # Extract point predictions and full standard errors
  x0 <- pred_unconditional$fit[1]
  y0 <- pred_unconditional$fit[2]
  sd_x <- pred_se_full[1]
  sd_y <- pred_se_full[2]
  
  # Extract residual correlation structure
  Sigma <- solve(crossprod(model$family$data$R))  # Residual covariance matrix
  rho <- cov2cor(Sigma)[1, 2]  # Residual correlation between x and y
  
  # Calculate ellipse rotation angle using analytical formula
  angle <- 0.5 * atan2(2 * rho * sd_x * sd_y, sd_x^2 - sd_y^2)
  
  # Scale for 90% confidence ellipse
  radius_mult <- sqrt(qchisq(0.90, df = 2))
  
  # Create the plot
  plot <- ggplot() +
    # 90% confidence ellipse
    geom_ellipse(aes(x0 = x0, y0 = y0,
                     a = sd_x * radius_mult,
                     b = sd_y * radius_mult,
                     angle = angle),
                 color = "blue", linewidth = 1, fill = NA) +
    
    # Marginal confidence intervals (rectangular box)
    annotate("rect", 
             xmin = (x0 - 1.645*sd_x),  # Approximately 90% CI
             xmax = (x0 + 1.645*sd_x),
             ymin = (y0 - 1.645*sd_y),
             ymax = (y0 + 1.645*sd_y),
             color = "darkgreen", fill = NA, linewidth = 1, linetype = "dashed") +
    
    # Reference points
    geom_point(aes(x = test_obs$x_obs[i], y = test_obs$y_obs[i]), 
               color = "red", size = 3) +     # Observed location
    geom_point(aes(x = test_obs$x_true[i], y = test_obs$y_true[i]), 
               color = "orange", size = 3) +  # True location
    geom_point(aes(x = x0, y = y0), 
               color = "black", size = 3, shape = 3) +  # Predicted location
    
    xlim(-2, 12) + ylim(-2, 12) +
    coord_fixed() +
    labs(title = paste0("Location as response - mv GAM, point = ", i),
         subtitle = "Orange = True; Red = Observed; Black = Estimated",
         caption = "Blue ellipse: 90% joint confidence region; Green box: marginal 90% intervals") +
    theme_bw()
  
  # Print uncertainty components for diagnostic purposes
  cat("Uncertainty decomposition for observation", i, ":\n")
  cat("  Unconditional SE (x):", round(pred_unconditional$se.fit[1], 4), "\n")
  cat("  Unconditional SE (y):", round(pred_unconditional$se.fit[2], 4), "\n")
  cat("  Full prediction SE (x):", round(sd_x, 4), "\n")
  cat("  Full prediction SE (y):", round(sd_y, 4), "\n")
  
  return(plot)
}



############################################################################
## PLOT LOCATION AND ERROR FOR EXAMPLE OBSERVATION - Location as Response RF
############################################################################

plot_prediction_diagnostics_rf <- function(test_obs, model, i) {
  
  # Loop over trees to get predictions for each
  tree_preds <- lapply(1:model$ntree, function(b) {
    predict(model, newdata = test_obs, get.tree = b)$regrOutput
  })
  
  # Extract per-tree predictions for x_true and y_true
  tree_preds_x <- sapply(tree_preds, function(tree) tree$x_true$predicted)
  tree_preds_y <- sapply(tree_preds, function(tree) tree$y_true$predicted)
  
  # Compute mean and SD across trees
  pred_mean_x <- rowMeans(tree_preds_x)
  pred_mean_y <- rowMeans(tree_preds_y)
  pred_sd_x <- apply(tree_preds_x, 1, sd)
  pred_sd_y <- apply(tree_preds_y, 1, sd)
  
  # Estimate correlation across trees per point
  pred_corr <- sapply(1:nrow(tree_preds_x), function(z) {
    cor(tree_preds_x[z, ], tree_preds_y[z, ])
  })
  
  # Combine results into a plotting dataframe
  plot_data <- data.frame(
    x_obs = test_obs$x_obs,
    y_obs = test_obs$y_obs,
    x_true = test_obs$x_true,
    y_true = test_obs$y_true,
    x_pred = pred_mean_x,
    y_pred = pred_mean_y,
    sd_x = pred_sd_x,
    sd_y = pred_sd_y,
    corr = pred_corr )
  
  # Create ellipse parameters
  # For each observation, construct proper ellipse parameters
  plot_data$angle <- 0.5 * atan2(2 * plot_data$corr * plot_data$sd_x * plot_data$sd_y, 
                                 plot_data$sd_x^2 - plot_data$sd_y^2)
  
  # Scale for 90% confidence ellipse
  radius_mult <- sqrt(qchisq(0.90, df = 2))
  plot_data$ellipse_a <- plot_data$sd_x * radius_mult
  plot_data$ellipse_b <- plot_data$sd_y * radius_mult
  
  # # More robust method?
  # # For each observation, construct covariance matrix and decompose
  # for(ii in 1:nrow(plot_data)) {
  #   # Construct covariance matrix
  #   cov_mat <- matrix(c(plot_data$sd_x[ii]^2,
  #                       plot_data$corr[ii] * plot_data$sd_x[ii] * plot_data$sd_y[ii],
  #                       plot_data$corr[ii] * plot_data$sd_x[ii] * plot_data$sd_y[ii],
  #                       plot_data$sd_y[ii]^2),
  #                     nrow = 2)
  # 
  #   # Eigendecomposition
  #   eig <- eigen(cov_mat)
  #   plot_data$angle[ii] <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])
  #   plot_data$ellipse_a[ii] <- sqrt(qchisq(0.95, df = 2) * eig$values[1])
  #   plot_data$ellipse_b[ii] <- sqrt(qchisq(0.95, df = 2) * eig$values[2])
  # }
  
  # Plot
  plot <- ggplot(plot_data[i, ]) +
    geom_point(aes(x = x_true, y = y_true), color = "orange", size = 3) +  # true
    geom_point(aes(x = x_obs, y = y_obs), color = "red", size = 3) +      # observed
    geom_point(aes(x = x_pred, y = y_pred), color = "black", shape = 3, size = 3) +  # predicted
    annotate("rect", 
             xmin = (plot_data$x_pred[i] - 1.645*plot_data$sd_x[i]),  #*2 seems overly pessimistic
             xmax = (plot_data$x_pred[i] + 1.645*plot_data$sd_x[i]),
             ymin = (plot_data$y_pred[i] - 1.645*plot_data$sd_y[i]),
             ymax = (plot_data$y_pred[i] + 1.645*plot_data$sd_y[i]),
             color = "darkgreen", fill = NA, linewidth = 1, linetype = "dashed") +
    geom_ellipse(
      aes(x0 = x_pred, y0 = y_pred, a = ellipse_a, b = ellipse_b, angle = angle),
      color = "blue", linewidth = 1, fill = NA ) +
    xlim(-2, 12) + ylim(-2, 12) +
    coord_fixed() +
    theme_bw() +
    labs(title = paste0("Location as response - mv RF, point = ", i),
         subtitle = "Orange = True; Red = Observed; Black = Estimated")
  
  return(list(results = plot_data, plot = plot))
}

