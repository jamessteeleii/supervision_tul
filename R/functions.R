# Read in and prepare data -----
read_prepare_data <- function(file) {
  data <- read_csv(file) |>
    janitor::clean_names()|>
    mutate(tul_centre = tul-120)
}

# Model fitting functions -----
rstan_setup <- function() {
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores() - 1)
}

hurdle_student_t_setup <- function() {
  hurdle_student_t <- 
    # Create a custom family that is logit if y = 0, student_t if not
    custom_family("hurdle_student_t", 
                  dpars = c("nu", "mu", "sigma", "hu"),
                  links = c("logm1", "identity", "log", "logit"),
                  lb = c(0, NA, 0, NA),
                  type = "real")
}

stan_funs_setup <- function() {
  stan_funs <- "
  real hurdle_student_t_lpdf(real y, real nu, real mu, real sigma, real hu) { 
    if (y == 0) { 
      return bernoulli_lpmf(1 | hu); 
    } else { 
      return bernoulli_lpmf(0 | hu) +  
             student_t_lpdf(y | nu, mu, sigma); 
    } 
  }
"
}

stan_vars_setup <- function(stan_funs) {
  # Prepare Stan code for use in brm()
  stanvars <- stanvar(scode = stan_funs, block = "functions")
}

# posterior predict functions
posterior_predict_hurdle_student_t <- function(i, prep, ...) {
  nu <- brms::get_dpar(prep, "nu", i = i)
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  theta <- brms::get_dpar(prep, "hu", i = i)
  
  hu <- runif(prep$ndraws, 0, 1)
  ifelse(hu < theta, 0, brms::rstudent_t(prep$ndraws, nu, mu,sigma))
}


posterior_epred_hurdle_student_t <- function(prep) {
  with(prep$dpars, mu * (1 - hu))
}

# Fit model to prior data in order to get prior distributions for experimental models -----

fit_model_prior_sample_tul <- function(data, stan_vars, hurdle_student_t) {
  brm(
    bf(
      tul_centre ~ 1 + core_assisted + (1 | location) + (1 | id) + (1 | machine),
      hu ~ 1 + core_assisted + (1 | location) + (1 | id) + (1 | machine)
    ),
    family = hurdle_student_t,
    stanvars = stan_vars,
    # prior = c(
    #   # set wide but constrained priors for just the fixed student t effects to realistic range for posterior predictions
    #   set_prior("student_t(3, 0, 15)", class = "b", lb = -60, ub = 60)
    # ), 
    data = data,
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 4000
  )
}
  

fit_model_prior_sample_rpe <- function(data) {
  ordbetareg(
    bf(
      session_rpe ~ 1 + core_assisted + (1|location),
      phi ~ 1 + core_assisted + (1|location)
    ),
    phi_reg = TRUE,
    true_bounds = c(6,20), 
    data = data,
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 4000
  )
}
  

# Fit model to experimental dataset -----

# Time under load outcome

# set priors from model parameters
set_priors_tul <- function(prior_model) {
  nu <- mean(get_dpar(prepare_predictions(prior_model), dpar = "nu"))
  params <- tidy(prior_model)
  
  priors_tul <- c(
    # Priors for Intercepts 
    set_prior(paste("student_t(",nu,",", params$estimate[1],",", params$std.error[1],")"),
              class = "Intercept"),
    set_prior(paste("logistic(", params$estimate[2],",", params$std.error[2],")"),
              class = "Intercept", dpar = "hu"),
    
    # Priors for core-assisted
    set_prior(paste("student_t(",nu,",", params$estimate[3],",", params$std.error[3],")"),
              class = "b", coef = "core_assistedcore"),
    set_prior(paste("logistic(", params$estimate[4],",", params$std.error[4],")"),
              class = "b", coef = "core_assistedcore", dpar = "hu"),
    
    # Priors for random effects for student_t component
    set_prior(paste("student_t(3,", params$estimate[6],",", params$std.error[6],")"),
              class = "sd", coef = "Intercept", group = "location"),
    
    set_prior(paste("student_t(3,", params$estimate[5],",", params$std.error[5],")"),
              class = "sd", coef = "Intercept", group = "location:id"),
    
    set_prior(paste("student_t(3,", params$estimate[7],",", params$std.error[7],")"),
              class = "sd", coef = "Intercept", group = "location:id:machine"),
    
    # Priors for random effects for hurdle component
    set_prior(paste("student_t(3,", params$estimate[9],",", params$std.error[9],")"),
              class = "sd", coef = "Intercept", group = "location", dpar = "hu"),
    
    set_prior(paste("student_t(3,", params$estimate[8],",", params$std.error[8],")"),
              class = "sd", coef = "Intercept", group = "location:id", dpar = "hu"),
    
    set_prior(paste("student_t(3,", params$estimate[10],",", params$std.error[10],")"),
              class = "sd", coef = "Intercept", group = "location:id:machine", dpar = "hu"),
    
    # Priors for sigma
    set_prior(paste("student_t(3,", params$estimate[11],",", params$std.error[11],")"),
              class = "sigma")
  )
}

# fit model
fit_model_tul <- function(data, stan_vars, hurdle_student_t, priors)
  brm(
    bf(
      tul_centre ~ 1 + core_assisted + (1 | location/id/machine),
      hu ~ 1 + core_assisted + (1 | location/id/machine)
    ),
    family = hurdle_student_t,
    stanvars = stan_vars,
    prior = priors, 
    data = data,
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 4000
  )

# Rating of perceive effort outcome
# set_priors_rpe <- function(prior_model) {
#   nu <- mean(get_dpar(prepare_predictions(prior_model), dpar = "nu"))
#   params <- tidy(prior_model)
#   
#   priors_rpe <- c(
#     # Priors for core training 
#     set_prior(paste("student_t(",nu,",", params$estimate[2],",", params$std.error[2],")"),
#               class = "b", coef = "core_assistedcore"),
#     set_prior(paste("logistic(", params$estimate[4],",", params$std.error[4],")"),
#               class = "b", coef = "core_assistedcore", dpar = "hu"),
#     
#     # Priors for assisted training i.e., supervision = 1
#     set_prior(paste("student_t(",nu,",", params$estimate[1],",", params$std.error[1],")"),
#               class = "b", coef = "core_assistedassisted"),
#     set_prior(paste("logistic(", params$estimate[3],",", params$std.error[3],")"),
#               class = "b", coef = "core_assistedassisted", dpar = "hu"),
#     
#     # Priors for random effects for student_t component
#     set_prior(paste("student_t(3,", params$estimate[6],",", params$std.error[6],")"),
#               class = "sd", coef = "Intercept", group = "location"),
#     
#     set_prior(paste("student_t(3,", params$estimate[5],",", params$std.error[5],")"),
#               class = "sd", coef = "Intercept", group = "location:id"),
#     
#     set_prior(paste("student_t(3,", params$estimate[7],",", params$std.error[7],")"),
#               class = "sd", coef = "Intercept", group = "location:id:machine"),
#     
#     # Priors for random effects for hurdle component
#     set_prior(paste("student_t(3,", params$estimate[9],",", params$std.error[9],")"),
#               class = "sd", coef = "Intercept", group = "location", dpar = "hu"),
#     
#     set_prior(paste("student_t(3,", params$estimate[8],",", params$std.error[8],")"),
#               class = "sd", coef = "Intercept", group = "location:id", dpar = "hu"),
#     
#     set_prior(paste("student_t(3,", params$estimate[10],",", params$std.error[10],")"),
#               class = "sd", coef = "Intercept", group = "location:id:machine", dpar = "hu")
#   )
# }
# 
# # fit model
# fit_model_rpe <- function(data, stan_vars, hurdle_student_t, priors)
#   brm(
#     bf(
#       rpe_centre ~ 0 + core_assisted + (1 | location/id/machine),
#       hu ~ 0 + core_assisted + (1 | location/id/machine)
#     ),
#     family = hurdle_student_t,
#     stanvars = stan_vars,
#     prior = priors, 
#     data = data,
#     chains = 4,
#     cores = 4,
#     seed = 1988,
#     warmup = 2000,
#     iter = 4000
#   )