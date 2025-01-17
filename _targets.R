# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c(
    "tidyverse",
    "brms",
    "ordbetareg",
    "bayesplot",
    "tidybayes",
    "rstan",
    "patchwork",
    "here",
    "broom.mixed"
  ),
  memory = "transient",
  format = "qs",
  garbage_collection = TRUE,
  storage = "worker",
  retrieval = "worker"
)
   
# Run the R scripts in the R/ folder with your custom functions:
tar_source("R/functions.R")

list(
  # Read and prepare sample of prior client data
  tar_target(
    prior_data_tul_file,
    here("data", "prior_sample_data.csv"),
    format = "file"
  ),
  tar_target(
    prior_data_tul, 
    read_prepare_data(prior_data_tul_file)
  ),
  
  tar_target(
    prior_data_rpe_file,
    here("data", "prior_sample_rpe_data.csv"),
    format = "file"
  ),
  tar_target(
    prior_data_rpe, 
    read_csv(prior_data_rpe_file)
  ),
  
  # Read and prepare study data
  tar_target(
    data_file,
    here("data", "data.csv"),
    format = "file"
  ),
  tar_target(
    data, 
    read_prepare_data(data_file)
  ),
  
  # Setup rstan to run chains in parallel
  tar_target(
    rstan, 
    rstan_setup()
  ),
  
  # Setup custom hurdle_student_t family
  tar_target(
    hurdle_student_t, 
    hurdle_student_t_setup()
  ),
  tar_target(
    stan_funs, 
    stan_funs_setup()
  ),
  tar_target(
    stan_vars, 
    stan_vars_setup(stan_funs)
  ),
  
  # Fit prior models
  tar_target(
    model_prior_sample_tul,
    fit_model_prior_sample_tul(prior_data_tul, stan_vars, hurdle_student_t)
  ),
  tar_target(
    model_prior_sample_rpe,
    fit_model_prior_sample_rpe(prior_data_rpe)
  ),
  
  # Fit experimental models
  tar_target(
    priors,
    set_priors_tul(model_prior_sample_tul)
  ),
  tar_target(
    model_tul,
    fit_model_tul(data, stan_vars, hurdle_student_t, priors)
  )
  
)
