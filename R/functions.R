# Read in and prepare data -----
read_prepare_data <- function(file) {
  data <- read_csv(file) |>
    janitor::clean_names() |>
    mutate(tul_centre = tul - 120)
}

# Model fitting functions -----
rstan_setup <- function() {
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores() - 1)
}

hurdle_student_t_setup <- function() {
  hurdle_student_t <-
    # Create a custom family that is logit if y = 0, student_t if not
    custom_family(
      "hurdle_student_t",
      dpars = c("nu", "mu", "sigma", "hu"),
      links = c("logm1", "identity", "log", "logit"),
      lb = c(0, NA, 0, NA),
      type = "real"
    )
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
  ifelse(hu < theta, 0, brms::rstudent_t(prep$ndraws, nu, mu, sigma))
}


posterior_epred_hurdle_student_t <- function(prep) {
  with(prep$dpars, mu * (1 - hu))
}

# Fit model to prior data in order to get prior distributions for experimental models -----

fit_model_prior_sample_tul <-
  function(data, stan_vars, hurdle_student_t) {
    brm(
      bf(
        tul_centre ~ 1 + core_assisted + (1 |
                                            location) + (1 | id) + (1 | machine),
        hu ~ 1 + core_assisted + (1 |
                                    location) + (1 | id) + (1 | machine)
      ),
      family = hurdle_student_t,
      stanvars = stan_vars,
      prior = c(
        # set wide but constrained priors for just the fixed effects
        set_prior("student_t(3, 0, 15)", class = "b", lb = -60, ub = 60)
      ),
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
      session_rpe ~ 1 + core_assisted + (1 | location),
      phi ~ 1 + core_assisted + (1 | location)
    ),
    phi_reg = TRUE,
    true_bounds = c(6, 20),
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
    set_prior(
      paste(
        "student_t(",
        nu,
        ",",
        params$estimate[1],
        ",",
        params$std.error[1],
        ")"
      ),
      class = "Intercept"
    ),
    set_prior(
      paste(
        "student_t(",
        nu,
        ",",
        params$estimate[2],
        ",",
        params$std.error[2],
        ")"
      ),
      class = "Intercept",
      dpar = "hu"
    ),
    
    # Priors for core-assisted
    set_prior(
      paste(
        "student_t(",
        nu,
        ",",
        params$estimate[3],
        ",",
        params$std.error[3],
        ")"
      ),
      class = "b",
      coef = "core_assistedcore"
    ),
    set_prior(
      paste(
        "student_t(",
        nu,
        ",",
        params$estimate[4],
        ",",
        params$std.error[4],
        ")"
      ),
      class = "b",
      coef = "core_assistedcore",
      dpar = "hu"
    ),
    
    # Priors for random effects for student_t component
    set_prior(
      paste(
        "student_t(3,",
        params$estimate[6],
        ",",
        params$std.error[6],
        ")"
      ),
      class = "sd",
      coef = "Intercept",
      group = "location"
    ),
    
    set_prior(
      paste(
        "student_t(3,",
        params$estimate[5],
        ",",
        params$std.error[5],
        ")"
      ),
      class = "sd",
      coef = "Intercept",
      group = "location:id"
    ),
    
    set_prior(
      paste(
        "student_t(3,",
        params$estimate[7],
        ",",
        params$std.error[7],
        ")"
      ),
      class = "sd",
      coef = "Intercept",
      group = "location:id:machine"
    ),
    
    # Priors for random effects for hurdle component
    set_prior(
      paste(
        "student_t(3,",
        params$estimate[9],
        ",",
        params$std.error[9],
        ")"
      ),
      class = "sd",
      coef = "Intercept",
      group = "location",
      dpar = "hu"
    ),
    
    set_prior(
      paste(
        "student_t(3,",
        params$estimate[8],
        ",",
        params$std.error[8],
        ")"
      ),
      class = "sd",
      coef = "Intercept",
      group = "location:id",
      dpar = "hu"
    ),
    
    set_prior(
      paste(
        "student_t(3,",
        params$estimate[10],
        ",",
        params$std.error[10],
        ")"
      ),
      class = "sd",
      coef = "Intercept",
      group = "location:id:machine",
      dpar = "hu"
    ),
    
    # Priors for sigma
    set_prior(
      paste(
        "student_t(3,",
        params$estimate[11],
        ",",
        params$std.error[11],
        ")"
      ),
      class = "sigma"
    )
  )
}

# fit model
fit_model_tul <-
  function(data, stan_vars, hurdle_student_t, priors) {
    brm(
      bf(
        tul_centre ~ 1 + core_assisted + (1 | location / id / machine),
        hu ~ 1 + core_assisted + (1 | location / id / machine)
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
  }


# Rating of perceived effort outcome

# set priors from model parameters
set_priors_rpe <- function(prior_model) {
  params <- tidy(prior_model)
  
  priors_rpe <- c(
    # Priors for Intercepts
    set_prior(
      paste(
        "student_t(3,",
        params$estimate[1],
        ",",
        params$std.error[1],
        ")"
      ),
      class = "Intercept"
    ),
    set_prior(
      paste(
        "student_t(3,",
        params$estimate[2],
        ",",
        params$std.error[2],
        ")"
      ),
      class = "Intercept",
      dpar = "phi"
    ),
    
    # Priors for core-assisted
    set_prior(
      paste(
        "student_t(3,",
        params$estimate[3],
        ",",
        params$std.error[3],
        ")"
      ),
      class = "b",
      coef = "core_assistedcore"
    ),
    set_prior(
      paste(
        "student_t(3,",
        params$estimate[4],
        ",",
        params$std.error[4],
        ")"
      ),
      class = "b",
      coef = "core_assistedcore",
      dpar = "phi"
    ),
    
    # Priors for random effects only for location (others left weakly regularising)
    set_prior(
      paste(
        "student_t(3,",
        params$estimate[5],
        ",",
        params$std.error[5],
        ")"
      ),
      class = "sd",
      coef = "Intercept",
      group = "location"
    ),
    
    set_prior(
      paste(
        "student_t(3,",
        params$estimate[6],
        ",",
        params$std.error[6],
        ")"
      ),
      class = "sd",
      coef = "Intercept",
      group = "location",
      dpar = "phi"
    )
    
  )
}

# fit model
fit_model_rpe <- function(data, priors) {
  ordbetareg(
    bf(
      rpe_e ~ 1 + core_assisted + (1 | location / id / machine),
      phi ~ 1 + core_assisted + (1 | location / id / machine)
    ),
    phi_reg = TRUE,
    true_bounds = c(0, 10),
    extra_prior = priors,
    data = data,
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 4000
  )
}

# Rating of perceived discomfort outcome

# fit model - note weakly regularising uninformed priors used for this outcome
fit_model_discomfort <- function(data) {
  ordbetareg(
    bf(
      rpe_d ~ 1 + core_assisted + (1 | location / id / machine),
      phi ~ 1 + core_assisted + (1 | location / id / machine)
    ),
    phi_reg = TRUE,
    true_bounds = c(0, 10),
    data = data,
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 4000
  )
}

# Model post-processing -----

# Time under load outcome

get_pred_draws_tul <- function(model, prior_model) {
  pred_draws <- bind_rows(
    avg_predictions(prior_model, 
                variables = "core_assisted",
                re_formula = NA) |>
      get_draws() |>
      mutate(distribution = "prior",
             draw = draw + 120),
    avg_predictions(model, 
                variables = "core_assisted",
                re_formula = NA) |>
      get_draws() |>
      mutate(distribution = "posterior",
             draw = draw + 120)
  ) |>
    mutate(core_assisted = factor(
      case_when(
        core_assisted == "core" ~ "Core",
        core_assisted == "assisted" ~ "Assisted"
      ),
      levels = c("Core", "Assisted")
    ))
}

get_pred_draws_tul_hu <- function(model, prior_model) {
  pred_draws <- bind_rows(
    avg_predictions(prior_model, 
                    variables = "core_assisted",
                    re_formula = NA,
                dpar = "hu") |>
      get_draws() |>
      mutate(distribution = "prior"),
    avg_predictions(model, 
                    variables = "core_assisted",
                    re_formula = NA,
                dpar = "hu") |>
      get_draws() |>
      mutate(distribution = "posterior")
  ) |>
    mutate(core_assisted = factor(
      case_when(
        core_assisted == "core" ~ "Core",
        core_assisted == "assisted" ~ "Assisted"
      ),
      levels = c("Core", "Assisted")
    ))
}

get_contrast_draws_tul <- function(model, prior_model) {
  contrast_draws <- bind_rows(
    avg_comparisons(prior_model, 
                    variables = "core_assisted",
                    re_formula = NA) |>
      get_draws() |>
      mutate(distribution = "prior",
             draw = draw),
    avg_comparisons(model, 
                    variables = "core_assisted",
                    re_formula = NA) |>
      get_draws() |>
      mutate(distribution = "posterior",
             draw = draw)
  ) 
}

get_contrast_draws_tul_hu <- function(model, prior_model) {
  contrast_draws <- bind_rows(
    avg_comparisons(prior_model, 
                    variables = "core_assisted",
                    re_formula = NA,
                    dpar = "hu") |>
      get_draws() |>
      mutate(distribution = "prior"),
    avg_comparisons(model, 
                    variables = "core_assisted",
                    re_formula = NA,
                    dpar = "hu") |>
      get_draws() |>
      mutate(distribution = "posterior")
  )
}

# Rating of perceived effort outcome

get_pred_draws_rpe <- function(model, prior_model) {
  pred_draws <- bind_rows(
    avg_predictions(prior_model, 
                    variables = "core_assisted",
                    re_formula = NA) |>
      get_draws() |>
      mutate(distribution = "prior"),
    avg_predictions(model, 
                    variables = "core_assisted",
                    re_formula = NA) |>
      get_draws() |>
      mutate(distribution = "posterior")
  ) |>
    mutate(core_assisted = factor(
      case_when(
        core_assisted == "core" ~ "Core",
        core_assisted == "assisted" ~ "Assisted"
      ),
      levels = c("Core", "Assisted")
    ))
}

get_contrast_draws_rpe <- function(model, prior_model) {
  contrast_draws <- bind_rows(
    avg_comparisons(prior_model, 
                    variables = "core_assisted",
                    re_formula = NA) |>
      get_draws() |>
      mutate(distribution = "prior"),
    avg_comparisons(model, 
                    variables = "core_assisted",
                    re_formula = NA) |>
      get_draws() |>
      mutate(distribution = "posterior")
  )
}

get_pred_draws_rpe_phi <- function(model, prior_model) {
  pred_draws <- bind_rows(
    predictions(prior_model) |>
      get_draws() |>
      mutate(distribution = "prior"),
    predictions(model) |>
      get_draws() |>
      mutate(distribution = "posterior")
  ) |>
    mutate(core_assisted = factor(
      case_when(
        core_assisted == "core" ~ "Core",
        core_assisted == "assisted" ~ "Assisted"
      ),
      levels = c("Core", "Assisted")
    ))
}

# Rating of perceived discomfort outcome

get_pred_draws_discomfort <- function(model) {
  pred_draws <- avg_predictions(model, 
                    variables = "core_assisted",
                    re_formula = NA) |>
      get_draws() |>
    mutate(distribution = "prior") |>
    mutate(core_assisted = factor(
      case_when(
        core_assisted == "core" ~ "Core",
        core_assisted == "assisted" ~ "Assisted"
      ),
      levels = c("Core", "Assisted")
    ))
}

get_contrast_draws_discomfort <- function(model) {
  contrast_draws <-avg_comparisons(model, 
                    variables = "core_assisted",
                    re_formula = NA) |>
      get_draws()
}

get_pred_draws_discomfort_phi <- function(model) {
  pred_draws <- predictions(model) |>
    mutate(core_assisted = factor(
      case_when(
        core_assisted == "core" ~ "Core",
        core_assisted == "assisted" ~ "Assisted"
      ),
      levels = c("Core", "Assisted")
    ))
}

# Plot results -----

# Raw plots for both prior sample and study data

# Time under load outcome
make_plot_prior_sample_tul <- function(data) {
  data |>
    mutate(core_assisted = factor(
      case_when(
        core_assisted == "core" ~ "Core",
        core_assisted == "assisted" ~ "Assisted"
      ),
      levels = c("Core", "Assisted")
    )) |>
    ggplot(aes(x = tul)) +
    geom_vline(xintercept = c(90, 120), linetype = "dashed",
               linewidth = 0.25
    ) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    scale_x_continuous(breaks = seq(60,180,10)) +
    facet_grid(core_assisted ~ .) +
    labs(x = "Time Under Load (seconds)") +
    theme_classic(base_size = 8) +
    theme(
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
}

make_plot_data_tul <- function(data) {
  data |>
    mutate(core_assisted = factor(
      case_when(
        core_assisted == "core" ~ "Core",
        core_assisted == "assisted" ~ "Assisted"
      ),
      levels = c("Core", "Assisted")
    )) |>
    ggplot(aes(
      x = core_assisted,
      y = tul,
      group = interaction(id, machine)
    )) +
    geom_hline(
      yintercept = c(90, 120),
      linetype = "dashed",
      alpha = 0.5,
      linewidth = 0.25
    ) +
    geom_line(
      position = position_jitter(width = 0.05, height = 0.05),
      alpha = 0.5,
      linewidth = 0.25
    ) + 
    scale_y_continuous(breaks = seq(60,600,30)) +
    labs(x = "Condition",
         y = "Time Under Load (Seconds)") +
    theme_classic(base_size = 8)
}

# Rating of perceived effort outcome
make_plot_prior_sample_rpe <- function(data) {
  data |>
    mutate(core_assisted = factor(
      case_when(
        core_assisted == "core" ~ "Core",
        core_assisted == "assisted" ~ "Assisted"
      ),
      levels = c("Core", "Assisted")
    )) |>
    ggplot(aes(x = session_rpe)) +
    geom_histogram(binwidth = 1, alpha = 0.5) +
    facet_grid(core_assisted ~ .) +
    scale_x_continuous(limits = c(5.5,20.5), breaks = seq(6,20)) +
    labs(x = "Session Rating of Perceived Effort (6-20 AU)") +
    theme_classic(base_size = 8) +
    theme(
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
}

make_plot_data_rpe <- function(data) {
  data <- data |>
    mutate(core_assisted = factor(
      case_when(
        core_assisted == "core" ~ "Core",
        core_assisted == "assisted" ~ "Assisted"
      ),
      levels = c("Core", "Assisted")
    )) 
  
  data |>
    ggplot(aes(
      x = core_assisted,
      y = rpe_e
    )) +
    geom_line(
      aes(      group = interaction(id, machine)),
      position = position_jitter(width = 0.05, height = 0.05),
      alpha = 0.5,
      linewidth = 0.25
    ) +
    stat_slab(
      data = filter(data, core_assisted == "Core"),
      density = "histogram",
      breaks = breaks_fixed(width = 0.5),
      align = align_center(),
      position = position_nudge(x = -0.1),
      side = "left"
    ) +
    stat_slab(
      data = filter(data, core_assisted == "Assisted"),
      density = "histogram",
      breaks = breaks_fixed(width = 0.5),
      align = align_center(),
      position = position_nudge(x = 0.1),
      side = "right"
    ) +
    scale_y_continuous(limits = c(-0.5,10.5), breaks = seq(0,10)) +
    labs(x = "Condition",
         y = "Rating of Perceived Effort (0-10 AU)") +
    theme_classic(base_size = 8)
}

# Rating of perceived discomfort outcome
make_plot_data_discomfort <- function(data) {
  data <- data |>
    mutate(core_assisted = factor(
      case_when(
        core_assisted == "core" ~ "Core",
        core_assisted == "assisted" ~ "Assisted"
      ),
      levels = c("Core", "Assisted")
    )) 
  
  data |>
    ggplot(aes(
      x = core_assisted,
      y = rpe_d
    )) +
    geom_line(
      aes(      group = interaction(id, machine)),
      position = position_jitter(width = 0.05, height = 0.05),
      alpha = 0.5,
      linewidth = 0.25
    ) +
    stat_slab(
      data = filter(data, core_assisted == "Core"),
      density = "histogram",
      breaks = breaks_fixed(width = 0.5),
      align = align_center(),
      position = position_nudge(x = -0.1),
      side = "left"
    ) +
    stat_slab(
      data = filter(data, core_assisted == "Assisted"),
      density = "histogram",
      breaks = breaks_fixed(width = 0.5),
      align = align_center(),
      position = position_nudge(x = 0.1),
      side = "right"
    ) +
    scale_y_continuous(limits = c(-0.5,10.5), breaks = seq(0,10)) +
    labs(x = "Condition",
         y = "Rating of Perceived Discomfort (0-10 AU)") +
    theme_classic(base_size = 8)
}

make_plot_combined_data <- function(plot1, plot2, plot3, plot4, plot5) {
  prior_plots <- (
    (plot1) + theme(strip.text = element_blank()) | (plot2)
  ) +
    plot_annotation(title = "Prior Sample Data Distributions",
                    subtitle = "Time Under Load and Session Rating of Perceived Effort") +
    plot_layout(axes = "collect")
  
  data_plots <- (
    (plot3) | (plot4) | (plot5)
  ) +
    plot_annotation(title = "Current Experimental Data Distributions",
                    subtitle = "Time Under Load, Rating of Perceived Effort, and Rating of Perceived Discomfort") +
    plot_layout(axes = "collect_x")
  
  
  wrap_elements(prior_plots) / wrap_elements(data_plots) 
}

# Model predictions and contrasts
make_plot_preds_tul <- function(pred_draws) {
  ggplot() +
    geom_vline(
      xintercept = c(90, 120),
      linetype = "dashed",
      alpha = 0.5,
      linewidth = 0.25
    ) +
    stat_slabinterval(
      aes(
        y = core_assisted,
        x = draw,
        color = distribution,
        fill = distribution
      ),
      data = filter(pred_draws, core_assisted == "Core"),
      point_interval = "median_qi",
      .width = .95,
      slab_alpha = 0.5,
      position = position_dodge(w = -0.1),
      size = 0.1
    ) +
    stat_slabinterval(
      aes(
        y = core_assisted,
        x = draw,
        color = distribution,
        fill = distribution
      ),
      data = filter(pred_draws, core_assisted == "Assisted"),
      point_interval = "median_qi",
      .width = .95,
      slab_alpha = 0.5,
      position = position_dodge(w = -0.1),
      size = 0.1
    ) +
    scale_color_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
    scale_x_continuous(limits = c(60, 180), breaks = seq(60,180,10)) +
    labs(
      y = "Condition",
      x = "Time Under Load (Seconds)",
      color = "Distribution",
      fill = "Distribution"
    ) +
    theme_classic(base_size = 8)
}

make_plot_contrasts_tul <- function(contrast_draws) {
  contrast_draws |>
    ggplot() +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      alpha = 0.5,
      linewidth = 0.25
    ) +
    stat_slabinterval(
      aes(
        x = draw,
        color = distribution,
        fill = distribution
      ),
      point_interval = "median_qi",
      .width = .95,
      slab_alpha = 0.5,
      position = position_dodge(w = -0.1),
      size = 0.1,
      scale = 0.5
    ) +
    scale_color_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
    # scale_x_continuous(limits = c(0, 200)) +
    labs(
      y = "Contrast (Core minus Assisted)",
      x = "Difference in Time Under Load (Seconds)",
      color = "Distribution",
      fill = "Distribution"
    ) +
    theme_classic(base_size = 8) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

make_plot_preds_tul_hu <- function(pred_draws) {
  ggplot() +
    geom_vline(
      xintercept = c(90, 120),
      linetype = "dashed",
      alpha = 0.5,
      linewidth = 0.25
    ) +
    stat_slabinterval(
      aes(
        y = core_assisted,
        x = draw,
        color = distribution,
        fill = distribution
      ),
      data = filter(pred_draws, core_assisted == "Core"),
      point_interval = "median_qi",
      .width = .95,
      slab_alpha = 0.5,
      position = position_dodge(w = -0.1),
      size = 0.1
    ) +
    stat_slabinterval(
      aes(
        y = core_assisted,
        x = draw,
        color = distribution,
        fill = distribution
      ),
      data = filter(pred_draws, core_assisted == "Assisted"),
      point_interval = "median_qi",
      .width = .95,
      slab_alpha = 0.5,
      position = position_dodge(w = -0.1),
      size = 0.1
    ) +
    scale_color_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
    scale_x_continuous(labels = scales::percent, limits = c(0,0.3)) +
    labs(
      y = "Condition",
      x = "Probability of Stopping at 120 Seconds (%)",
      color = "Distribution",
      fill = "Distribution"
    ) +
    theme_classic(base_size = 8)
}

make_plot_contrasts_tul_hu <- function(contrast_draws) {
  contrast_draws |>
    ggplot() +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      alpha = 0.5,
      linewidth = 0.25
    ) +
    stat_slabinterval(
      aes(
        x = draw,
        color = distribution,
        fill = distribution
      ),
      point_interval = "median_qi",
      .width = .95,
      slab_alpha = 0.5,
      position = position_dodge(w = -0.1),
      size = 0.1,
      scale = 0.5
    ) +
    scale_color_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
    # scale_x_continuous(limits = c(0, 200)) +
    scale_x_continuous(labels = scales::percent) +
    labs(
      y = "Contrast (Core minus Assisted)",
      x = "Difference in Probability of Stopping at 120 Seconds (%)",
      color = "Distribution",
      fill = "Distribution"
    ) +
    theme_classic(base_size = 8) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

make_plot_combined_tul <- function(plot1, plot2, plot3, plot4) {
  
  tul_plots <- (
    (plot1 + guides(fill = "none", color = "none")) | (plot2 + guides(fill = "none", color = "none"))
  ) +
    plot_annotation(title = "Time Under Load",
                    subtitle = "Global Grand Means for Predictions and Contrasts")
  
  tul_hu_plots <- (
    (plot3 + guides(fill = "none", color = "none")) | (plot4)
  ) +
    plot_annotation(title = "Stopping at 120 seconds?",
                    subtitle = "Global Grand Means for Predictions and Contrasts") + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom")
  
  wrap_elements(tul_plots) / wrap_elements(tul_hu_plots) 
  
}

make_plot_preds_rpe <- function(pred_draws) {
  ggplot() +
    geom_vline(
      xintercept = c(90, 120),
      linetype = "dashed",
      alpha = 0.5,
      linewidth = 0.25
    ) +
    stat_slabinterval(
      aes(
        y = core_assisted,
        x = draw,
        color = distribution,
        fill = distribution
      ),
      data = filter(pred_draws, core_assisted == "Core"),
      point_interval = "median_qi",
      .width = .95,
      slab_alpha = 0.5,
      position = position_dodge(w = -0.1),
      size = 0.1
    ) +
    stat_slabinterval(
      aes(
        y = core_assisted,
        x = draw,
        color = distribution,
        fill = distribution
      ),
      data = filter(pred_draws, core_assisted == "Assisted"),
      point_interval = "median_qi",
      .width = .95,
      slab_alpha = 0.5,
      position = position_dodge(w = -0.1),
      size = 0.1
    ) +
    scale_color_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
    scale_x_continuous(limits = c(0,1), labels = scales::percent, breaks = seq(0,1,0.1)) +
    labs(
      y = "Condition",
      x = "Rating of Perceived Effort (%)",
      color = "Distribution",
      fill = "Distribution"
    ) +
    theme_classic(base_size = 8)
}

make_plot_contrasts_rpe <- function(contrast_draws) {
  contrast_draws |>
    ggplot() +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      alpha = 0.5,
      linewidth = 0.25
    ) +
    stat_slabinterval(
      aes(
        x = draw,
        color = distribution,
        fill = distribution
      ),
      point_interval = "median_qi",
      .width = .95,
      slab_alpha = 0.5,
      position = position_dodge(w = -0.1),
      size = 0.1,
      scale = 0.5
    ) +
    scale_color_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
    scale_fill_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
    scale_x_continuous(labels = scales::percent, breaks = seq(-1,1,0.05)) +
    labs(
      y = "Contrast (Core minus Assisted)",
      x = "Difference in Rating of Perceived Effort (%)",
      color = "Distribution",
      fill = "Distribution"
    ) +
    theme_classic(base_size = 8) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

make_plot_preds_discomfort <- function(pred_draws) {
  ggplot() +
    geom_vline(
      xintercept = c(90, 120),
      linetype = "dashed",
      alpha = 0.5,
      linewidth = 0.25
    ) +
    stat_slabinterval(
      aes(
        y = core_assisted,
        x = draw * 10
        ),
      data = filter(pred_draws, core_assisted == "Core"),
      point_interval = "median_qi",
      .width = .95,
      slab_alpha = 0.5,
      position = position_dodge(w = -0.1),
      size = 0.1
    ) +
    stat_slabinterval(
      aes(
        y = core_assisted,
        x = draw * 10
        ),
      data = filter(pred_draws, core_assisted == "Assisted"),
      point_interval = "median_qi",
      .width = .95,
      slab_alpha = 0.5,
      position = position_dodge(w = -0.1),
      size = 0.1
    ) +
    scale_x_continuous(limits = c(0,10), breaks = seq(0,10)) +
    labs(
      y = "Condition",
      x = "Rating of Perceived Discomfort (0-10 AU)"
    ) +
    theme_classic(base_size = 8)
}

make_plot_contrasts_discomfort <- function(contrast_draws) {
  contrast_draws |>
    ggplot() +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      alpha = 0.5,
      linewidth = 0.25
    ) +
    stat_slabinterval(
      aes(
        x = draw * 10
      ),
      point_interval = "median_qi",
      .width = .95,
      slab_alpha = 0.5,
      position = position_dodge(w = -0.1),
      size = 0.1,
      scale = 0.5
    ) +
    labs(
      y = "Contrast (Core minus Assisted)",
      x = "Difference in Rating of Perceived Discomfort (0-10 AU)",
      color = "Distribution",
      fill = "Distribution"
    ) +
    theme_classic(base_size = 8) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

make_plot_combined_rpe_discomfort <- function(plot1, plot2, plot3, plot4) {
  rpe_plots <- (
    (plot1 + guides(fill = "none", color = "none")) | (plot2)
  ) +
    plot_annotation(title = "Rating of Perceived Effort",
                    subtitle = "Global Grand Means for Predictions and Contrasts") + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom")
  
  discomfort_plots <- (
    (plot3) | (plot4)
  ) +
    plot_annotation(title = "Rating of Perceived Discomfort",
                    subtitle = "Global Grand Means for Predictions and Contrasts")
  
  
  wrap_elements(rpe_plots) / wrap_elements(discomfort_plots) 
  
}
