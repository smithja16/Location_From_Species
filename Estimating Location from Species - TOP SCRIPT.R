
####################################################################
#####  Catch Location Refinement using Species Compositions    #####
#####                  J.A.Smith NSW DPI, 2025                 #####
#####               james.a.smith@dpi.nsw.gov.au               #####
####################################################################

## NOTE: some code was edited and improved by AI: Claude-Sonnet-4 & ChatGPT-4.5

## The overall goal is to use the spatial signal inherent in species compositions
## to estimate/improve inaccurate survey locations. This was created mainly for
## improving fishing catch records with inaccurate or missing location data.

## This code tests 3 main approaches:

## 1) A hierarchical joint Bayesian GLM with location as a latent variable estimated
## for the inaccurate records; this gains skill from the subset of accurate records

## 2) A joint Bayesian GLM trained on accurate locations, which is then
## inverted (used as a likelihood function) to find the most likely locations
## given new data (the inaccurate location data)

## 3) Multivariate GAM and ranbdom forest which model location as a function of
## species composition, which flips the usual modelling order.

## The code below simulates catch data containing accurate and inaccurate
## location values, then tests out each modelling approach.


##############################
##  Load libraries          ##
##############################

library(rstan)
library(brms)
library(MASS)
library(mgcv)
library(randomForestSRC)

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(patchwork)


##############################
##  Load functions          ##
##############################

source("Simulate Catch Data Function.R")
source("Inverse Prediction Function.R")
source("Processing Functions.R")


##############################
##  1) Generate the data    ##
##############################

## Define specific spatial patterns for each species (here, n=5)
## This creates linear spatial patterns, which allows for simpler 
## brms and Stan coding and better comparability of approaches

species_patterns <- matrix(c(
  -1.5,  0.0,    # sp1: strong west gradient
  0.0,  1.5,     # sp2: strong north gradient
  1.0, -1.0,     # sp3: southeast gradient
  -1.0, -1.0,    # sp4: southwest gradient
  1.5,  1.5 ),   # sp5: strong northeast gradient
  ncol = 2, byrow = TRUE)


## Generate data with specified patterns
sim_result <- simulate_fishing_data(
  n_accurate = 300,                 # Generate 300 accurate points (training)
  n_approx = 300,                   # Generate 300 innacurate points (testing)
  n_species = 5,                    # Use 5 species
  spatial_pattern = "specified",    # Uses the species patterns matrix
  specified_patterns = species_patterns,
  loc_error_sd = 3.0,               # Error of uncertain locations from true
  spatial_strength = 2.0,           # Spatial signal for species
  month_effect_sd = 0.05,           # Small monthly variation
  mean_abundance = 5,               # Target mean count of 5
  dist_family = "negbin",           # Negative binomial or Poisson
  negbin_size = 2,                  # Size (dispersion) parameter for neg bin
  seed = 117 )

fish_data <- sim_result$data        # The simulated data for later use
sim_result$true_associations        # The simulated spatial associations
comparison <- verify_species_associations(sim_result)
comparison$visualize_patterns()  # plot the spatial distributions



############################################
##  2) Fit and plot the hierarchical GLM  ##
############################################

## First write and load the Stan model code
## This is hard-coded as a negative binomial GLM with linear terms

source("Save Stan model.R")  #writes and saves the Stan model

## Prepare data for Stan
## Identify accurate and inaccurate observations
known_ids <- which(fish_data$accurate)
unknown_ids <- which(!fish_data$accurate)

## Create data object which matches Stan model syntax
stan_data <- list(
  N = nrow(fish_data),
  S = 5,  # 5 species
  M = 12, # 12 months
  N_known = length(known_ids),
  N_unknown = length(unknown_ids),
  known_ids = known_ids,
  unknown_ids = unknown_ids,
  loc_obs = as.matrix(fish_data[, c("x_obs", "y_obs")]),
  month = fish_data$month,
  counts = as.matrix(fish_data[, paste0("sp", 1:5)]),
  loc_sd = 1.5 )  # Standard deviation for location prior

## Fit the model
options(mc.cores = 4)  # Cores for parallel processing (same as chains)
rstan_options(auto_write = TRUE)

fit <- stan(
  file = "location_estimation_model.stan",
  data = stan_data,
  iter = 4000,
  chains = 4,
  # seed = 117,  # If you want same output each time
  control = list(adapt_delta = 0.9) )

## Extract results and calculate accuracy of estimated locations
posterior <- rstan::extract(fit)
loc_unknown_samples <- posterior$loc_true_unknown

accuracy_summary_lv <- summarize_accuracy_lv(
  fish_data = fish_data,
  loc_unknown_samples = loc_unknown_samples )

accuracy_summary_lv$fish_data  # All results
lv_summary <- accuracy_summary_lv$summary  # Summary

## Calculate marginal credible intervals
cred_intervals <- apply(loc_unknown_samples, c(2,3), function(x) {
  quantile(x, probs = c(0.05, 0.95)) })

n_unknown <- dim(loc_unknown_samples)[2]
ci_df <- data.frame(
  point = rep(1:n_unknown, each = 2),
  coord = rep(c("x", "y"), times = n_unknown),
  lower = as.numeric(cred_intervals[1,,]),
  upper = as.numeric(cred_intervals[2,,]) )

## Visualise an example observation
i <- 100  # index of the inaccurate observations
samples <- loc_unknown_samples[, i, ]
p_lv <- plot_prediction_diagnostics_lv(
  samples = samples,
  fish_data = fish_data ) 



##################################################
##  2) Fit and plot the Inverse Prediction GLM  ##
##################################################

## Specify each species' model structure
bf_sp1 <- bf(sp1 ~ x_true + y_true + monthf)  # Linear spatial terms
bf_sp2 <- bf(sp2 ~ x_true + y_true + monthf)  # Month as a factor
bf_sp3 <- bf(sp3 ~ x_true + y_true + monthf)
bf_sp4 <- bf(sp4 ~ x_true + y_true + monthf)
bf_sp5 <- bf(sp5 ~ x_true + y_true + monthf)

multi_formula <- bf_sp1 + bf_sp2 + bf_sp3 + bf_sp4 + bf_sp5

## Fit the model
fit_mv <- brm(
  formula = multi_formula,
  data = fish_data %>% filter(accurate),  # Fit using accurate data only
  family = negbinomial(),
  chains = 4, cores = 4, iter = 4000,
  control = list(adapt_delta = 0.9) )

#pp_check(fit_mv, type="scatter_avg",resp="sp1")

## Invert the model to predict inaccurate locations
test_obs <- fish_data %>% filter(!accurate)

test_with_priors <- pbapply::pblapply(1:nrow(test_obs), function(i) {
  inverse_predict_with_prior_mv_custom(
    obs = test_obs[i, ],
    fit_mv = fit_mv,
    x_var = "x_obs",        # x and y are split up so that they can be different units (e.g. latitude and depth)
    y_var = "y_obs",
    model_x = "x_true",
    model_y = "y_true",
    grid_radius_x = 5,      # In map units; make large enough to cover likely locations
    grid_radius_y = 5,
    grid_size_x = 31,       # Odd number ensure that the observed point is at the grid center
    grid_size_y = 31,
    prior_sd_x = 1.5,       # Increase this as grid radius and resolution increases (otherwise grid limits are un-priored); grid_radius/2 is a good default
    prior_sd_y = 1.5,
    use_normal_prior = TRUE,
    use_grid_constraint = FALSE,
    lik_fun = lik_fun_nb )  # Negative binomial to match the SDM
} )

## Extract results and calculate accuracy of estimated locations
accuracy_summary_inv <- summarize_accuracy_inv(
  test_obs = test_obs,
  inverse_results = test_with_priors )

accuracy_summary_inv$results  # All results
inv_summary <- accuracy_summary_inv$summary  #Summary

## Estimate uncertainty
uncertainty_results <- lapply(test_with_priors, function(res) {
  estimate_uncertainty_inv(res$posterior)
})

test_obs_uncertainty <- test_obs %>%
  mutate(
    x_pred = sapply(test_with_priors, function(r) r$map$x),
    y_pred = sapply(test_with_priors, function(r) r$map$y),
    x_lower = sapply(uncertainty_results, function(r) r$x_lower),
    x_upper = sapply(uncertainty_results, function(r) r$x_upper),
    y_lower = sapply(uncertainty_results, function(r) r$y_lower),
    y_upper = sapply(uncertainty_results, function(r) r$y_upper),
    x_se = sapply(uncertainty_results, function(r) r$x_sd),
    y_se = sapply(uncertainty_results, function(r) r$y_sd) )

## Visualise an example observation
i <- 100  # pick a point to explore (of the updated locations)
obs <- test_obs[i, ]
result <- test_with_priors[[i]]
plots <- plot_prediction_diagnostics_inv(obs, result)
(plots$prior + plots$likelihood) / (plots$posterior + plots$uncertainty)
p_inv <- plots$uncertainty



#######################################################
##  3) Fit and plot the Location As Repsonse models  ##
#######################################################

#### 3a) The multivariate GAM

## Fit the model
## MVN can model horizontal, vertical, and diagonal spatial patterns
fit_gam <- gam(list(
  x_true ~ monthf + s(sp1) + s(sp2) + s(sp3) + s(sp4) + s(sp5),
  y_true ~ monthf + s(sp1) + s(sp2) + s(sp3) + s(sp4) + s(sp5)),
  data = fish_data %>% filter(accurate),
  family = mvn(d = 2))

# summary(fit_gam)

## Examine covariance of x and y
Sigma <- solve(crossprod(fit_gam$family$data$R))
colnames(Sigma) <- rownames(Sigma) <- c("x_true", "y_true")
Sigma  #covariance matrix
cor_mat <- cov2cor(Sigma)
cor_mat  #correlation matrix

## Predict inaccurate locations
test_obs <- fish_data %>% filter(!accurate)
pred_gam <- predict(fit_gam, newdata=test_obs, se.fit=T)

## Calculate accuracy of estimated locations
accuracy_summary_gam <- summarize_accuracy_gam_rf(
  test_obs = test_obs,
  predictions = pred_gam$fit )

accuracy_summary_gam$results  # All results
gam_summary <- accuracy_summary_gam$summary  # Summary

## Visualise an example observation
i <- 100
p_gam <- plot_prediction_diagnostics_gam(
  test_obs = test_obs,
  model = fit_gam,
  i = i )


#### 3b) The multivariate random forest

## Fit the model
formula_mrf <- as.formula("Multivar(x_true, y_true) ~ monthf + sp1 + sp2 + sp3 + sp4 + sp5") 

fit_mrf <- rfsrc(formula_mrf,
                 data = fish_data %>% filter(accurate),
                 ntree = 1000,
                 mtry = 2,
                 nodesize = 3)

## Predict inaccurate locations
pred_rf <- predict(fit_mrf, newdata = test_obs)
x_pred <- pred_rf$regrOutput$x_true$predicted
y_pred <- pred_rf$regrOutput$y_true$predicted

## Calculate accuracy of estimated locations
accuracy_summary_rf <- summarize_accuracy_gam_rf(
  test_obs = test_obs,
  predictions = cbind(x_pred, y_pred) )

accuracy_summary_rf$results  # All results
rf_summary <- accuracy_summary_rf$summary  # SUmmary

## Visualise an example observation
i <- 100
results <- plot_prediction_diagnostics_rf(
  test_obs = test_obs,
  model = fit_mrf,
  i = i )
p_rf <- results$plot



####################################
##  3) Compare the Approaches     ##
####################################

(p_lv + p_inv) / (p_gam + p_rf)

lv_summary
inv_summary
gam_summary
rf_summary


