
####################################################################
#####  Catch Location Refinement using Species Compositions    #####
#####                  J.A.Smith NSW DPI, 2025                 #####
#####               james.a.smith@dpi.nsw.gov.au               #####
####################################################################

## This script provides the data simulation function and a function
## to visualise the simulated species distributions.


################################################
## SIMULATE SPECIES CATCH AND LOCATIONS FUNCTION
################################################

## This generates a set of locations and counts at those locations. It is also
## hard-coded to add a 'month' covariate, but could be adapted for any
## additional covariates.

## There are many 'extras' in this function, some of which were not used in the
## published analysis.

simulate_fishing_data <- function(
    n_accurate = 300,                  # Number of accurate location records
    n_approx = 300,                    # Number of approximate location records  
    n_species = 5,                     # Number of species
    loc_error_sd = 1.0,                # SD of location measurement error
    species_correlation = 0.3,         # Base correlation between species
    dist_family = "negbin",            # Distribution family: "poisson" or "negbin" only
    negbin_size = 2,                   # Size parameter for negative binomial
    zero_inflation = 0,                # Probability of structural zeros
    month_effect_sd = 0.1,             # SD of monthly effects
    x_range = c(0, 10),                # Range of x coordinates
    y_range = c(0, 10),                # Range of y coordinates
    spatial_pattern = "random",        # Options: "random", "specified", "mixed"
    specified_patterns = NULL,         # Matrix of x,y coefficients for "specified" option
    use_autocorrelation = FALSE,       # Whether to use spatial autocorrelation for "random" option
    spatial_range = 3,                 # Range parameter for spatial autocorrelation
    spatial_sigma = 1,                 # Sigma parameter for spatial autocorrelation
    spatial_strength = 1.0,            # Overall strength of spatial effects
    mean_abundance = 5,                # Target mean abundance for each species
    max_abundance = 30,                # Target maximum abundance for each species
    seed = NULL ) {                    # Random seed for reproducibility
  
  # "specified" uses a provided matrix of species spatial preferences
  # "mixed" uses a random sample of pre-made spatial preferences
  # "random" uses a random set of spatial preferences
  # Each will be consistent if a "seed" is specified
  
  # Input validation
  if(!dist_family %in% c("poisson", "negbin")) {
    stop("dist_family must be 'poisson' or 'negbin'")
  }
  if(zero_inflation < 0 || zero_inflation > 1) {
    stop("zero_inflation must be between 0 and 1")
  }
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Total number of records
  n_total <- n_accurate + n_approx
  
  # True locations (uniformly across specified ranges)
  true_locs <- cbind(
    runif(n_total, min = x_range[1], max = x_range[2]), 
    runif(n_total, min = y_range[1], max = y_range[2]) )
  colnames(true_locs) <- c("x", "y")
  
  # Standardize locations for stability
  true_locs_std <- scale(true_locs)
  
  # Measurement error for approximate locations
  if(n_approx > 0) {
    direction_bias <- runif(n_approx, 0, 2*pi)  #Random direction for each observation
    error_magnitude <- rlnorm(n_approx, log(loc_error_sd * 0.8), 0.5)  # Log-normal error magnitudes
    approx_locs <- true_locs[(n_accurate + 1):n_total, ] + 
      cbind(sin(direction_bias) * error_magnitude, 
            cos(direction_bias) * error_magnitude)
  } else {
    approx_locs <- matrix(nrow = 0, ncol = 2)
  }
  
  # Observed locations: accurate ones are true, others are noisy
  obs_locs <- rbind(true_locs[1:n_accurate, ], approx_locs)
  loc_accuracy <- c(rep(TRUE, n_accurate), rep(FALSE, n_approx))
  
  # Month covariate (1-12)
  month <- sample(1:12, n_total, replace = TRUE)
  
  # Create a data structure to store the species associatoins
  true_associations <- list()
  
  # Standard patterns that can be used
  standard_patterns <- list(
    north = c(0, 1),       # Preference for northern areas (high y)
    south = c(0, -1),      # Preference for southern areas (low y)
    east = c(1, 0),        # Preference for eastern areas (high x)
    west = c(-1, 0),       # Preference for western areas (low x)
    northeast = c(1, 1),   # Preference for northeastern areas
    northwest = c(-1, 1),  # Preference for northwestern areas
    southeast = c(1, -1),  # Preference for southeastern areas
    southwest = c(-1, -1), # Preference for southwestern areas
    center = c(-0.5, -0.5, 0.5, 0.5), # Preference for central areas (quadratic)
    edge = c(0.5, 0.5, 0.5, 0.5) )     # Preference for edge areas (quadratic)
  
  # Initiate matrix for spatial effects
  spatial_effects <- matrix(0, n_total, n_species)
  
  # Create spatial effects based on the specified pattern
  if(spatial_pattern == "specified" && !is.null(specified_patterns)) {
    # Use user-specified patterns
    if(is.matrix(specified_patterns) && nrow(specified_patterns) == n_species) {
      species_patterns <- specified_patterns
    } else {
      stop("specified_patterns must be a matrix with n_species rows")
    }
    
    # Apply the patterns to calculate spatial effects
    for(s in 1:n_species) {
      if(ncol(specified_patterns) == 2) {
        # Linear effects
        spatial_effects[, s] <- spatial_strength * (
          specified_patterns[s, 1] * true_locs_std[, 1] + 
            specified_patterns[s, 2] * true_locs_std[, 2])
      } else if(ncol(specified_patterns) == 4) {
        # Quadratic effects
        spatial_effects[, s] <- spatial_strength * (
          specified_patterns[s, 1] * true_locs_std[, 1] + 
            specified_patterns[s, 2] * true_locs_std[, 2] +
            specified_patterns[s, 3] * true_locs_std[, 1]^2 + 
            specified_patterns[s, 4] * true_locs_std[, 2]^2)
      } else {
        stop("specified_patterns must have either 2 columns (linear) or 4 columns (quadratic)")
      }
    }
    
    true_associations$pattern_type <- "specified"
    true_associations$patterns <- species_patterns
    
  } else if(spatial_pattern == "mixed") {
    # Use a mix of standard patterns
    available_patterns <- names(standard_patterns)
    selected_patterns <- sample(available_patterns, n_species, replace = TRUE)
    
    species_patterns <- matrix(NA, n_species, 4)  # Use 4 columns for both linear and quadratic
    
    for(s in 1:n_species) {
      pattern_name <- selected_patterns[s]
      pattern <- standard_patterns[[pattern_name]]
      
      if(length(pattern) == 2) {
        # Linear pattern
        species_patterns[s, 1:2] <- pattern
        species_patterns[s, 3:4] <- 0  # No quadratic terms
        
        spatial_effects[, s] <- spatial_strength * (
          pattern[1] * true_locs_std[, 1] + 
            pattern[2] * true_locs_std[, 2])
      } else {
        # Quadratic pattern
        species_patterns[s, ] <- pattern
        
        spatial_effects[, s] <- spatial_strength * (
          pattern[1] * true_locs_std[, 1] + 
            pattern[2] * true_locs_std[, 2] +
            pattern[3] * true_locs_std[, 1]^2 + 
            pattern[4] * true_locs_std[, 2]^2)
      }
    }
    
    true_associations$pattern_type <- "mixed"
    true_associations$patterns <- species_patterns
    true_associations$pattern_names <- selected_patterns
    
  } else {
    # Use random patterns
    if(use_autocorrelation) {
      # Spatial autocorrelation
      for(s in 1:n_species) {
        # Create distance matrix
        D <- as.matrix(dist(true_locs))
        # Exponential covariance function
        Sigma <- spatial_sigma^2 * exp(-D/spatial_range)
        # Ensure matrix is positive definite
        diag(Sigma) <- diag(Sigma) + 1e-8
        # Generate multivariate normal with spatial correlation
        spatial_effects[, s] <- spatial_strength * MASS::mvrnorm(1, rep(0, n_total), Sigma)
      }
      
      true_associations$pattern_type <- "autocorrelation"
      
    } else {
      # Random linear effects
      species_patterns <- matrix(rnorm(n_species * 2), n_species, 2)
      
      # Normalize patterns to have unit length
      for(s in 1:n_species) {
        species_patterns[s, ] <- species_patterns[s, ] / sqrt(sum(species_patterns[s, ]^2))
      }
      
      for(s in 1:n_species) {
        spatial_effects[, s] <- spatial_strength * (
          species_patterns[s, 1] * true_locs_std[, 1] + 
            species_patterns[s, 2] * true_locs_std[, 2])
      }
      
      true_associations$pattern_type <- "random_linear"
      true_associations$patterns <- species_patterns
    }
  }
  
  # Species intercepts and month effects
  intercepts <- rep(log(mean_abundance), n_species)  # Intercepts set to target mean abundance
  month_effects <- matrix(rnorm(n_species * 12, 0, month_effect_sd), n_species, 12)
  
  # Linear predictors before species correlations
  log_lambda <- matrix(NA, n_total, n_species)
  for(s in 1:n_species) {
    log_lambda[, s] <- intercepts[s] + spatial_effects[, s] + month_effects[s, month]
  }
  
  # Add optional species correlation structure
  if(species_correlation > 0) {
    # Create correlation matrix
    species_corr_matrix <- matrix(species_correlation, n_species, n_species)
    diag(species_corr_matrix) <- 1  # Diagonal is 1
    
    # Make sure the correlation matrix is positive definite
    eigen_decomp <- eigen(species_corr_matrix)
    if(any(eigen_decomp$values <= 0)) {
      species_corr_matrix <- species_corr_matrix + diag(0.01, n_species)
    }
    
    # Generate correlated noise
    chol_corr <- chol(species_corr_matrix)
    correlated_effects <- matrix(rnorm(n_total * n_species), n_total, n_species) %*% chol_corr
    
    # Scale and add the correlated effects
    correlated_effects <- correlated_effects * 0.2  # Smaller value to avoid overwhelming spatial effects
    log_lambda <- log_lambda + correlated_effects
    
    # Store the correlation matrix
    true_associations$species_correlation_matrix <- species_corr_matrix
  }
  
  # Generate counts based on the specified distribution
  counts <- matrix(NA, n_total, n_species)
  
  # Apply a limit to extreme values in log_lambda to avoid extremely large counts
  log_lambda_capped <- pmin(pmax(log_lambda, log(0.1)), log(max_abundance))
  
  if(dist_family == "poisson") {
    # Poisson counts
    lambda <- exp(log_lambda_capped)
    for(s in 1:n_species) {
      counts[, s] <- rpois(n_total, lambda = lambda[, s])
    }
  } else if(dist_family == "negbin") {
    # Negative binomial counts
    mu <- exp(log_lambda_capped)
    for(s in 1:n_species) {
      counts[, s] <- MASS::rnegbin(n_total, mu = mu[, s], theta = negbin_size)
    }
  }
  
  # Add optional zero-inflation
  if(zero_inflation > 0) {
    zero_mask <- matrix(runif(n_total * n_species) < zero_inflation, n_total, n_species)
    counts[zero_mask] <- 0
  }
  
  # Assemble into a data frame
  fish_data <- data.frame(
    id = 1:n_total,
    x_obs = obs_locs[, 1],
    y_obs = obs_locs[, 2],
    x_true = true_locs[, 1],
    y_true = true_locs[, 2],
    accurate = loc_accuracy,
    month = month,
    monthf = as.factor(month) )
  
  # Add species counts
  species_counts <- as.data.frame(counts)
  colnames(species_counts) <- paste0("sp", 1:n_species)
  fish_data <- cbind(fish_data, species_counts)
  
  # Add grid cell information
  grid_size <- 2.5  # size of each reporting grid cell (relative to x_range and y_range)
  fish_data$grid_id_x <- floor(fish_data$x_true / grid_size)
  fish_data$grid_id_y <- floor(fish_data$y_true / grid_size)
  fish_data$grid_id <- paste0(fish_data$grid_id_x, "_", fish_data$grid_id_y)
  fish_data$grid_min_x <- grid_size*fish_data$grid_id_x
  fish_data$grid_min_y <- grid_size*fish_data$grid_id_y
  fish_data$grid_max_x <- grid_size*(fish_data$grid_id_x+1)
  fish_data$grid_max_y <- grid_size*(fish_data$grid_id_y+1)
  
  # Create description of spatial patterns
  pattern_descriptions <- data.frame(
    species = paste0("sp", 1:n_species),
    stringsAsFactors = FALSE )
  
  if(true_associations$pattern_type %in% c("specified", "mixed", "random_linear")) {
    pattern_descriptions$x_effect <- ifelse(true_associations$patterns[, 1] > 0, "increases with x", 
                                            ifelse(true_associations$patterns[, 1] < 0, "decreases with x", "neutral to x"))
    pattern_descriptions$y_effect <- ifelse(true_associations$patterns[, 2] > 0, "increases with y", 
                                            ifelse(true_associations$patterns[, 2] < 0, "decreases with y", "neutral to y"))
    
    # Add direction description
    pattern_descriptions$spatial_pattern <- apply(true_associations$patterns[, 1:2], 1, function(pattern) {
      x_coef <- pattern[1]
      y_coef <- pattern[2]
      
      if(abs(x_coef) < 0.1 && abs(y_coef) < 0.1) return("No clear pattern")
      
      directions <- c()
      if(y_coef > 0.3) directions <- c(directions, "north")
      if(y_coef < -0.3) directions <- c(directions, "south")
      if(x_coef > 0.3) directions <- c(directions, "east")
      if(x_coef < -0.3) directions <- c(directions, "west")
      
      if(length(directions) == 0) {
        if(abs(x_coef) > abs(y_coef)) {
          if(x_coef > 0) return("weak east")
          else return("weak west")
        } else {
          if(y_coef > 0) return("weak north")
          else return("weak south")
        }
      } else {
        return(paste(directions, collapse = "-"))
      }
    })
    
    # Add quadratic terms if present
    if(ncol(true_associations$patterns) >= 4) {
      has_quad <- apply(true_associations$patterns[, 3:4], 1, function(quad) any(abs(quad) > 0.1))
      quad_desc <- apply(true_associations$patterns[, 3:4], 1, function(quad) {
        x_quad <- quad[1]
        y_quad <- quad[2]
        
        if(abs(x_quad) < 0.1 && abs(y_quad) < 0.1) return("")
        
        x_part <- if(abs(x_quad) > 0.1) {
          ifelse(x_quad > 0, "peaks at edges in x", "peaks in center in x")
        } else ""
        
        y_part <- if(abs(y_quad) > 0.1) {
          ifelse(y_quad > 0, "peaks at edges in y", "peaks in center in y")
        } else ""
        
        if(x_part != "" && y_part != "") {
          return(paste(x_part, "and", y_part))
        } else if(x_part != "") {
          return(x_part)
        } else {
          return(y_part)
        }
      })
      
      pattern_descriptions$quadratic_effect <- ifelse(has_quad, quad_desc, "none")
    }
    
    if("pattern_names" %in% names(true_associations)) {
      pattern_descriptions$named_pattern <- true_associations$pattern_names
    }
    
  } else if(true_associations$pattern_type == "autocorrelation") {
    # For autocorrelation, we need to analyze the resulting patterns
    for(s in 1:n_species) {
      # Calculate correlations with coordinates
      x_cor <- cor(true_locs[, 1], spatial_effects[, s])
      y_cor <- cor(true_locs[, 2], spatial_effects[, s])
      
      pattern_descriptions$x_effect[s] <- ifelse(x_cor > 0.2, "tends to increase with x", 
                                                 ifelse(x_cor < -0.2, "tends to decrease with x", "complex pattern with x"))
      pattern_descriptions$y_effect[s] <- ifelse(y_cor > 0.2, "tends to increase with y", 
                                                 ifelse(y_cor < -0.2, "tends to decrease with y", "complex pattern with y"))
      
      pattern_descriptions$spatial_pattern[s] <- "complex spatial autocorrelation"
    }
  }
  
  true_associations$pattern_descriptions <- pattern_descriptions
  
  # Return both the data and the true associations
  return(list(
    data = fish_data,
    true_associations = true_associations ))
}



#####################################################
## VISUALISE SIMULATED SPECIES DISTRIBUTIONS FUNCTION
#####################################################

verify_species_associations <- function(sim_result) {
  data <- sim_result$data
  true_assoc <- sim_result$true_associations
  
  species_cols <- grep("^sp", names(data), value = TRUE)
  
  # Calculate observed correlations
  observed_correlations <- data.frame(
    species = species_cols,
    x_correlation = sapply(species_cols, function(sp) {
      if(sd(data[[sp]]) > 0) cor(data$x_true, data[[sp]]) else NA
    }),
    y_correlation = sapply(species_cols, function(sp) {
      if(sd(data[[sp]]) > 0) cor(data$y_true, data[[sp]]) else NA
    }),
    stringsAsFactors = FALSE
  )
  
  # Fit simple models to quantify relationships 
  fit_models <- lapply(species_cols, function(sp) {
    if(all(data[[sp]] == 0)) {
      return(NULL)  # Skip if all zeros
    }
    
    if(all(data[[sp]] %in% c(0, 1))) {
      # Logistic regression for binary data
      fit <- glm(paste0(sp, " ~ x_true + y_true"), data = data, family = "binomial")
    } else {
      # Poisson regression for count data
      fit <- glm(paste0(sp, " ~ x_true + y_true"), data = data, family = "poisson")
    }
    return(fit)
  })
  
  # Compare true vs. observed associations
  comparison <- list(
    true_associations = true_assoc$pattern_descriptions,
    observed_correlations = observed_correlations )
  
  # Add visualization function
  comparison$visualize_patterns <- function() {
    
    # Prepare data for plotting
    plot_data <- data %>%
      dplyr::select(x_true, y_true, starts_with("sp")) %>%
      tidyr::pivot_longer(cols = starts_with("sp"),
                          names_to = "species",
                          values_to = "count")
    
    # Log-transform counts (with offset to avoid log(0))
    plot_data$log_count <- log1p(plot_data$count)
    
    # Create plot for each species
    p <- ggplot(plot_data, aes(x = x_true, y = y_true, color = log_count)) +
      geom_point(size = 2, alpha = 0.8) +
      scale_color_viridis_c(option = "plasma") +
      facet_wrap(~ species) +
      theme_minimal() +
      labs(title = "Spatial distribution of species (log-transformed counts)",
           x = "X coordinate", y = "Y coordinate", color = "log(count + 1)")
    
    return(p)
  }
  
  return(comparison)
}
