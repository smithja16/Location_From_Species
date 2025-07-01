
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
    geom_point(data = samples_df, aes(x = x, y = y),
               alpha = 0.05, color = "blue") +  # Posterior samples
    geom_point(data = mean_df, aes(x = x, y = y),
               color = "black", shape = 3, size = 4, stroke=2) +  # Posterior mean
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
