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
    "ggdist",
    "tidybayes",
    "marginaleffects",
    "lme4",
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
  # Read and prepare sample of prior client data -----
  tar_target(
    prior_data_tul_file,
    here("data", "prior_sample_data.csv"),
    format = "file"
  ),
  tar_target(
    prior_data_tul, 
    read_prepare_data(prior_data_tul_file) |>
      mutate(across(2:5, as.factor))
  ),
  
  tar_target(
    prior_data_rpe_file,
    here("data", "prior_sample_rpe_data.csv"),
    format = "file"
  ),
  tar_target(
    prior_data_rpe, 
    read_csv(prior_data_rpe_file) |>
      mutate(across(2:4, as.factor))
  ),
  
  # Read and prepare study data ----- 
  tar_target(
    data_file,
    here("data", "data.csv"),
    format = "file"
  ),
  tar_target(
    data, 
    read_prepare_data(data_file) |>
      mutate(across(2:5, as.factor))
  ),
  
  # Setup rstan to run chains in parallel -----
  tar_target(
    rstan, 
    rstan_setup()
  ),
  
  # Setup custom hurdle_student_t family -----
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
  
  # Fit prior models -----
  tar_target(
    model_prior_sample_tul,
    fit_model_prior_sample_tul(prior_data_tul, stan_vars, hurdle_student_t)
  ),
  tar_target(
    model_prior_sample_rpe,
    fit_model_prior_sample_rpe(prior_data_rpe)
  ),
  
  # Fit experimental models -----
  tar_target(
    priors_tul,
    set_priors_tul(model_prior_sample_tul)
  ),
  tar_target(
    model_tul,
    fit_model_tul(data, stan_vars, hurdle_student_t, priors_tul)
  ),
  
  tar_target(
    priors_rpe,
    set_priors_rpe(model_prior_sample_rpe)
  ),
  tar_target(
    model_rpe,
    fit_model_rpe(data, priors_rpe)
  ),
  
  tar_target(
    model_discomfort,
    fit_model_discomfort(data)
  ),
  
  
  
  # Model checks -----
  
  # rhat plots
  tar_target(
    rhat_prior_sample_tul,
    make_rhat_plot(model_prior_sample_tul)
  ),
  
  tar_target(
    rhat_prior_sample_rpe,
    make_rhat_plot(model_prior_sample_rpe)
  ),
  
  tar_target(
    rhat_tul,
    make_rhat_plot(model_tul)
  ),
  
  tar_target(
    rhat_rpe,
    make_rhat_plot(model_rpe)
  ),
  
  tar_target(
    rhat_discomfort,
    make_rhat_plot(model_discomfort)
  ),
  
  # trace plots
  tar_target(
    trace_plot_prior_sample_tul,
    make_trace_plot(model_prior_sample_tul)
  ),
  
  tar_target(
    trace_plot_prior_sample_rpe,
    make_trace_plot(model_prior_sample_rpe)
  ),
  
  tar_target(
    trace_plot_tul,
    make_trace_plot(model_tul)
  ),
  
  tar_target(
    trace_plot_rpe,
    make_trace_plot(model_rpe)
  ),
  
  tar_target(
    trace_plot_discomfort,
    make_trace_plot(model_discomfort)
  ),
  
  # pp check
  tar_target(
    pp_check_prior_sample_tul,
    make_pp_check(model_prior_sample_tul)
  ),
  
  tar_target(
    pp_check_prior_sample_rpe,
    make_pp_check(model_prior_sample_rpe)
  ),
  
  tar_target(
    pp_check_tul,
    make_pp_check(model_tul)
  ),
  
  tar_target(
    pp_check_rpe,
    make_pp_check(model_rpe)
  ),
  
  tar_target(
    pp_check_discomfort,
    make_pp_check(model_discomfort)
  ),
  
  # Model post processing -----
  tar_target(
    pred_draws_tul,
    get_pred_draws_tul(model_tul, model_prior_sample_tul)
  ),
  
  tar_target(
    pred_draws_tul_hu,
    get_pred_draws_tul_hu(model_tul, model_prior_sample_tul)
  ),
  
  tar_target(
    contrast_draws_tul,
    get_contrast_draws_tul(model_tul, model_prior_sample_tul)
  ),
  
  tar_target(
    contrast_draws_tul_hu,
    get_contrast_draws_tul_hu(model_tul, model_prior_sample_tul)
  ),
  
  tar_target(
    pred_draws_rpe,
    get_pred_draws_rpe(model_rpe, model_prior_sample_rpe)
  ),
  
  tar_target(
    contrast_draws_rpe,
    get_contrast_draws_rpe(model_rpe, model_prior_sample_rpe)
  ),
  
  tar_target(
    pred_draws_rpe_phi,
    get_pred_draws_rpe_phi(model_rpe, model_prior_sample_rpe)
  ),
  
  tar_target(
    pred_draws_discomfort,
    get_pred_draws_discomfort(model_discomfort)
  ),
  
  tar_target(
    contrast_draws_discomfort,
    get_contrast_draws_discomfort(model_discomfort)
  ),
  
  tar_target(
    pred_draws_discomfort_phi,
    get_pred_draws_discomfort_phi(model_discomfort)
  ),
  
  # Plotting data and results -----
  
  # Data plots
  tar_target(
    plot_prior_sample_tul,
    make_plot_prior_sample_tul(prior_data_tul)
  ),
  tar_target(
    plot_data_tul,
    make_plot_data_tul(data)
  ),
  tar_target(
    plot_prior_sample_rpe,
    make_plot_prior_sample_rpe(prior_data_rpe)
  ),
  tar_target(
    plot_data_rpe,
    make_plot_data_rpe(data)
  ),
  tar_target(
    plot_data_discomfort,
    make_plot_data_discomfort(data)
  ),
  
  tar_target(
    plot_combined_data,
    make_plot_combined_data(plot_prior_sample_tul, plot_prior_sample_rpe,
                            plot_data_tul, plot_data_rpe, plot_data_discomfort)
  ),
  
  # Model results plots
  tar_target(
    plot_preds_tul,
    make_plot_preds_tul(pred_draws_tul)
  ),
  tar_target(
    plot_contrasts_tul,
    make_plot_contrasts_tul(contrast_draws_tul)
  ),
  
  tar_target(
    plot_preds_tul_hu,
    make_plot_preds_tul_hu(pred_draws_tul_hu)
  ),
  tar_target(
    plot_contrasts_tul_hu,
    make_plot_contrasts_tul_hu(contrast_draws_tul_hu)
  ),
  
  tar_target(
    plot_combined_tul,
    make_plot_combined_tul(plot_preds_tul, plot_contrasts_tul,
                           plot_preds_tul_hu, plot_contrasts_tul_hu)
  ),
  
  tar_target(
    plot_preds_rpe,
    make_plot_preds_rpe(pred_draws_rpe)
  ),
  tar_target(
    plot_contrasts_rpe,
    make_plot_contrasts_rpe(contrast_draws_rpe)
  ),
  
  tar_target(
    plot_preds_discomfort,
    make_plot_preds_discomfort(pred_draws_discomfort)
  ),
  tar_target(
    plot_contrasts_discomfort,
    make_plot_contrasts_discomfort(contrast_draws_discomfort)
  ),
  
  tar_target(
    plot_combined_rpe_discomfort,
    make_plot_combined_rpe_discomfort(plot_preds_rpe, plot_contrasts_rpe,
                           plot_preds_discomfort, plot_contrasts_discomfort)
  ),
  
  # Make plots into tiffs
  tar_target(
    plot_combined_data_tiff,
    ggsave(
      "plots/plot_combined_data.tiff", 
      plot_combined_data, 
      device = "tiff",
      dpi = 300,
      w = 7.5,
      h = 7.5
    )
  ),
  
  tar_target(
    plot_combined_tul_tiff,
    ggsave(
      "plots/plot_combined_tul.tiff", 
      plot_combined_tul, 
      device = "tiff",
      dpi = 300,
      w = 7.5,
      h = 7
    )
  ),
  
  tar_target(
    plot_combined_rpe_discomfort_tiff,
    ggsave(
      "plots/plot_combined_rpe_discomfort.tiff", 
      plot_combined_rpe_discomfort, 
      device = "tiff",
      dpi = 300,
      w = 7.5,
      h = 7
    )
  ),
  
  # Additional descriptive model of load progression -----
  
  tar_target(
    prior_loadprog_data_file,
    here("data", "prior_sample_loadprog_data.csv"),
    format = "file"
  ),
  tar_target(
    prior_loadprog_data, 
    read_csv(prior_loadprog_data_file) |>
      janitor::clean_names()
  ),
  
  tar_target(
    model_loadprog,
    fit_loadprog_model(prior_loadprog_data)
  ),
  
  tar_target(
    preds_loadprog,
    get_preds_loadprog(model_loadprog)
  ),
  
  tar_target(
    slopes_loadprog,
    get_slopes_loadprog(model_loadprog)
  ),
  
  tar_target(
    plot_loadprog,
    make_plot_loadprog(preds_loadprog, slopes_loadprog)
  ),
  
  tar_target(
    plot_loadprog_tiff,
    ggsave(
      "plots/plot_loadprog.tiff", 
      plot_loadprog, 
      device = "tiff",
      dpi = 300,
      w = 7.5,
      h = 5
    )
  ),
  
  # Add flow chart for modelling strategy -----
  tar_target(
    flow_chart,
    make_flow_chart()
  ),
  
  tar_target(
    flow_chart_tiff,
    ggsave(
      "plots/flow_chart.tiff", 
      flow_chart, 
      device = "tiff",
      dpi = 300,
      w = 10,
      h = 7.5
    )
  )
  
)
