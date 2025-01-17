

library(tidyverse)
library(brms)
library(bayesplot)
library(tidybayes)
library(marginaleffects)
library(rstan)
library(patchwork)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 1)

kieser_theme <- theme(
  plot.background = element_rect(fill = "#0084CE"),
  panel.grid = element_line(colour = "#FFFFFF", size = 0.1),
  panel.background = element_rect(fill = "#7FC1E6"),
  axis.title = element_text(colour = "white"),
  axis.text = element_text(colour = "white"),
  strip.background = element_rect(fill = "#3FA2DA", colour = NA),
  strip.text = element_text(colour = "white"),
  plot.title = element_text(colour = "white"),
  plot.subtitle = element_text(colour = "white"),
  plot.caption = element_text(colour = "white")
)

# Read csv as data frame into environment - Note: change FILENAME
prior_data <- read_csv("data/prior_sample_data.csv") |>
  mutate(across(2:5, as.factor)) |>
  mutate(tul_centre = tul-120)

data <- read_csv("data/data.csv") |>
  janitor::clean_names() |>
  mutate(across(2:5, as.factor))


# data |>
#   group_by(supervision) |>
#   summarise(mean = mean(rpe_d, na.rm=TRUE),
#             sd = sd(rpe_d, na.rm=TRUE)) |>
#   ggplot(aes(x=supervision)) +
#   # geom_hline(yintercept = 120, linetype = "dashed") +
#   geom_line(aes(x=supervision, y=rpe_d, group=interaction(id, machine)),
#             data=data,
#             position = position_jitter(width=0.05, height=0.05),
#             alpha = 0.5) +
#   geom_pointrange(aes(y=mean, ymin=mean-sd, ymax=mean+sd),
#                   position = position_nudge(x=c(-0.1,0.1))) +
#   scale_x_discrete(labels = c("Core", "Assisted")) +
#   # scale_y_continuous(breaks = seq(0, 600, by = 50)) +
#   labs(
#     x = "Session Type",
#     y = "Time Under Load"
#   ) +
#   kieser_theme

  

# custom hurdle student_t

hurdle_student_t <- 
  # Create a custom family that is logit if y = 0, normal/student_t if not
  custom_family("hurdle_student_t", 
                dpars = c("nu", "mu", "sigma", "hu"),
                links = c("logm1", "identity", "log", "logit"),
                lb = c(0, NA, 0, NA),
                type = "real")

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

# Prepare Stan code for use in brm()
stanvars <- stanvar(scode = stan_funs, block = "functions")

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
  

small_sample_core_ids <- prior_data |> 
  group_by(core_assisted) |>
  slice_sample(n=100) |> 
  filter(core_assisted == "core") 

small_sample_core <- prior_data |>
  filter(id %in% small_sample_core_ids$id)

model_prior_sample_core_tul <-
  brm(
    bf(
      tul_centre ~ 1 + (1 | location) + (1 | id) + (1 | machine),
      hu ~ 1 + (1 | location) + (1 | id) + (1 | machine)
    ),
    family = hurdle_student_t,
    stanvars = stanvars,
    data = small_sample_core,
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 4000,
  )


plot(model_prior_sample_core_tul)

pp_check(model_prior_sample_core_tul, type = "hist", binwidth = 1) +
  scale_x_continuous(limits = c(-60,60))

pp_check(model_prior_sample_core_tul)

ppc_hist(model_prior_sample_core_tul, binwidth = 1)


pred_draws <- model_prior_sample_core_tul |> 
  predicted_draws(newdata = expand_grid(location = NA,
                                        id = NA,
                                        machine = NA),
                  re_formula = NA) |>
  mutate(tul_pred = .prediction + 120)


pred_plot <- pred_draws |>
  mutate(
    color = if_else(round(tul_pred) == 120, "red", "gray") 
  )|>
  ggplot(aes(x=tul_pred, color=color, fill=color)) +
  geom_vline(xintercept = c(90, 120), linetype = "dashed") +
  geom_histogram(binwidth = 1, alpha = 0.5) +
  scale_color_identity() +
  scale_fill_identity() +
  labs(x = "Predicted Time Under Load") +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  )


prior_sample_plot / pred_plot

get_prior(model_prior_sample_core_tul)

# fit intercept only models to each group from prior client sample data
get_prior(tul ~ 1 + (1 | location / id / machine),
          data = prior_data |> filter(core_assisted == "core"),
          family = "hurdle_lognormal")

get_prior(tul ~ 1 + (1 | location / id / machine),
          data = prior_data |> filter(core_assisted == "assisted"),
          family = "lognormal")

# sample priors
model_sample_core_tul_prior <-
  brm(
    tul_centre ~ 1 + (1 | location / id / machine),
    family = hurdl,
    data = prior_data |> filter(core_assisted == "core"),
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 8000,
    sample_prior = "only"
  )

pp_check(model_sample_core_tul_prior)

model_sample_assisted_tul_prior <-
  brm(
    tul ~ 1 + (1 | location / id / machine),
    family = "lognormal",
    data = prior_data |> filter(core_assisted == "assisted"),
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 8000,
    sample_prior = "only"
  )

# models
model_sample_core_tul <-
  brm(
    tul ~ 1 + (1 | location / id / machine),
    family = "lognormal",
    data = prior_data |> filter(core_assisted == "core"),
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 8000
  )

model_sample_assisted_tul <-
  brm(
    tul ~ 1 + (1 | location / id / machine),
    family = "lognormal",
    data = prior_data |> filter(core_assisted == "assisted"),
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 8000
  )



# set priors from model parameters


prior_tul <- c(
  # We set this prior based on the assumption that most core members tend to train until they hit the upper 120s TUL range without continuing to failure
  set_prior("lognormal(125, 15, -5)", class = "b", coef = "supervision0"),
  
  # We set the prior on the supervised group such that the difference in modes is roughly 0.5 of an SD if we were to assume tul is typically normally distributed
  set_prior("normal(127.5, 10)", class = "b", coef = "supervision1")
)


tibble(
  x = seq(0, 200, length = (160 - 50) * 1000),
  y1 = dskew_normal(
    x,
    xi = 125,
    omega = 15,
    alpha = -5
  ),
  y2 = dnorm(x, 127.5, 10)
) |>
  ggplot() +
  geom_vline(xintercept = c(90, 120), linetype = "dashed") +
  geom_line(aes(x, y1)) +
  geom_line(aes(x, y2)) +
  labs(x = "Time Under Load",
       y = "Density")

model_tul_prior <-
  brm(
    tul ~ 0 + supervision + (1 | location / id / machine),
    data = data,
    chains = 4,
    cores = 4,
    seed = 1988,
    warmup = 2000,
    iter = 8000,
    prior = prior_tul,
    # control = list(adapt_delta = 0.99),
    sample_prior = "only"
  )

model_tul <- brm(
  tul ~ 0 + supervision + (1 | location / id / machine),
  data = data,
  chains = 4,
  cores = 4,
  seed = 1988,
  warmup = 2000,
  iter = 8000,
  prior = prior_tul,
  # control = list(adapt_delta = 0.99),
  # sample_prior = "only"
)

plot(model_tul)
bayesplot::pp_check(model_tul, type = "dens_overlay_grouped", group = "supervision")


# draws <- bind_rows(
#   data |>
#     modelr::data_grid(core_assisted = c("core", "assisted")) |>
#     add_epred_draws(model_tul,
#                     re_formula = NA) |>
#         mutate(distribution = "posterior",
#                .value = .value+120),
#   prior_data |>
#     modelr::data_grid(core_assisted = c("core", "assisted")) |>
#     add_epred_draws(model_prior_sample_tul,
#                     re_formula = NA) |>
#     mutate(distribution = "prior",
#            .value = .value+120),
#   ) |>
#     mutate(core_assisted = factor(
#       case_when(
#         core_assisted == "core" ~ "Core",
#         core_assisted == "assisted" ~ "Assisted"
#       ),
#       levels = c("Core", "Assisted")
#     ))

# mu plots

draws <- bind_rows(
  model_prior_sample_tul |>
    gather_draws(b_core_assistedcore, b_core_assistedassisted) |>
    mutate(distribution = "prior",
           .value = .value+120),
  model_tul |>
    gather_draws(b_core_assistedcore, b_core_assistedassisted) |>
    mutate(distribution = "posterior",
           .value = .value+120)
) |>
  mutate(core_assisted = factor(
    case_when(
      .variable == "b_core_assistedcore" ~ "Core",
      .variable == "b_core_assistedassisted" ~ "Assisted"
    ),
    levels = c("Core", "Assisted")
  ))

contrasts <- bind_rows(
  model_prior_sample_tul |>
    spread_draws(b_core_assistedcore, b_core_assistedassisted) |>
    mutate(distribution = "prior",
           b_core_assistedcore = b_core_assistedcore+120,
           b_core_assistedassisted = b_core_assistedassisted+120,
           .contrast = b_core_assistedassisted - b_core_assistedcore),
  model_tul |>
    spread_draws(b_core_assistedcore, b_core_assistedassisted) |>
    mutate(distribution = "posterior",
           b_core_assistedcore = b_core_assistedcore+120,
           b_core_assistedassisted = b_core_assistedassisted+120,
           .contrast = b_core_assistedassisted - b_core_assistedcore)
) 


draws_hu <- bind_rows(
  model_prior_sample_tul |>
    gather_draws(b_hu_core_assistedcore, b_hu_core_assistedassisted) |>
    mutate(distribution = "prior",
           .value = plogis(.value)),
  model_tul |>
    gather_draws(b_hu_core_assistedcore, b_hu_core_assistedassisted) |>
    mutate(distribution = "posterior",
           .value = plogis(.value))
) |>
  mutate(core_assisted = factor(
    case_when(
      .variable == "b_hu_core_assistedcore" ~ "Core",
      .variable == "b_hu_core_assistedassisted" ~ "Assisted"
    ),
    levels = c("Core", "Assisted")
  ))

contrasts_hu <- bind_rows(
  model_prior_sample_tul |>
    spread_draws(b_hu_core_assistedcore, b_hu_core_assistedassisted) |>
    mutate(distribution = "prior",
           b_hu_core_assistedcore = plogis(b_hu_core_assistedcore),
           b_hu_core_assistedassisted = plogis(b_hu_core_assistedassisted),
           .contrast = b_hu_core_assistedassisted - b_hu_core_assistedcore),
  model_tul |>
    spread_draws(b_hu_core_assistedcore, b_hu_core_assistedassisted) |>
    mutate(distribution = "posterior",
           b_hu_core_assistedcore = plogis(b_hu_core_assistedcore),
           b_hu_core_assistedassisted = plogis(b_hu_core_assistedassisted),
           .contrast = b_hu_core_assistedassisted - b_hu_core_assistedcore),
  
) 

raw_plot <- data |>
  mutate(core_assisted = factor(
    case_when(core_assisted == "core" ~ "Core",
              core_assisted == "assisted" ~ "Assisted"),
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
  scale_y_continuous(limits = c(0, 600)) +
  labs(x = "Condition",
       y = "Time Under Load (Seconds)",
       title = "Raw data") +
  theme_classic(base_size = 8)


model_plot <- ggplot() +
  geom_hline(
    yintercept = c(90, 120),
    linetype = "dashed",
    alpha = 0.5,
    linewidth = 0.25
  ) +
  stat_slabinterval(
    aes(
      x = core_assisted,
      y = .value,
      color = distribution,
      fill = distribution
    ),
    data = filter(draws, core_assisted == "Core"),
    point_interval = "mean_qi",
    .width = .95,
    side = c("left"),
    slab_alpha = 0.5,
    position = position_dodge(w = -0.1),
    size = 0.1
  ) +
  stat_slabinterval(
    aes(
      x = core_assisted,
      y = .value,
      color = distribution,
      fill = distribution
    ),
    data = filter(draws, core_assisted == "Assisted"),
    point_interval = "mean_qi",
    .width = .95,
    side = c("right"),
    slab_alpha = 0.5,
    position = position_dodge(w = -0.1),
    size = 0.1
  ) +
  scale_color_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
  # scale_y_continuous(limits = c(0, 200)) +
  labs(
    x = "Condition",
    y = "Time Under Load (Seconds)",
    color = "Distribution",
    fill = "Distribution"
  ) +
  theme_classic(base_size = 8)

contrast_plot <- contrasts |>
  ggplot() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    alpha = 0.5,
    linewidth = 0.25
  ) +
  stat_slabinterval(
    aes(
      y = .contrast,
      color = distribution,
      fill = distribution
    ),
    point_interval = "mean_qi",
    .width = .95,
    slab_alpha = 0.5,
    position = position_dodge(w = -0.1),
    size = 0.1,
    scale = 0.5
  ) +
  scale_color_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
  # scale_y_continuous(limits = c(0, 200)) +
  labs(
    x = "Contrast (Assisted minus Core)",
    y = "Difference in Time Under Load (Seconds)"
  ) +
  guides(color = "none",
         fill = "none") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_blank(),
         axis.ticks.x = element_blank())
  

model_hu_plot <- ggplot() +
  stat_slabinterval(
    aes(
      x = core_assisted,
      y = .value,
      color = distribution,
      fill = distribution
    ),
    data = filter(draws_hu, core_assisted == "Core"),
    point_interval = "mean_qi",
    .width = .95,
    side = c("left"),
    slab_alpha = 0.5,
    position = position_dodge(w = -0.1),
    size = 0.1
  ) +
  stat_slabinterval(
    aes(
      x = core_assisted,
      y = .value,
      color = distribution,
      fill = distribution
    ),
    data = filter(draws_hu, core_assisted == "Assisted"),
    point_interval = "mean_qi",
    .width = .95,
    side = c("right"),
    slab_alpha = 0.5,
    position = position_dodge(w = -0.1),
    size = 0.1
  ) +
  scale_color_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
  scale_y_continuous(labels = scales::percent, limits = c(0,0.3)) +
  labs(
    x = "Condition",
    y = "Probability of Stopping at 120 Seconds (%)",
    color = "Distribution",
    fill = "Distribution"
  ) +
  theme_classic(base_size = 8)

contrast_hu_plot <- contrasts_hu |>
  ggplot() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    alpha = 0.5,
    linewidth = 0.25
  ) +
  stat_slabinterval(
    aes(
      y = .contrast,
      color = distribution,
      fill = distribution
    ),
    point_interval = "mean_qi",
    .width = .95,
    slab_alpha = 0.5,
    position = position_dodge(w = -0.1),
    size = 0.1,
    scale = 0.5
  ) +
  scale_color_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"), limits = rev) +
  scale_y_continuous(labels = scales::percent, limits = c(-0.25,0.05)) +
  labs(
    x = "Contrast (Assisted minus Core)",
    y = "Difference in Probability of Stopping at 120 Seconds (%)"
  ) +
  guides(color = "none",
         fill = "none") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

((model_plot + contrast_plot ) / (model_hu_plot + contrast_hu_plot)) +
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(legend.position = "bottom")

(
  (raw_plot) / ((contrast_plot + model_plot) + 
    plot_layout(axes = "collect",
                guides = "collect",
                widths = c(1, 2)) &
    theme(legend.position = "bottom"))
) 




# rpe

prior_data_rpe <- read_csv("data/prior_sample_rpe_data.csv")

model_rpe_prior <- ordbetareg(bf(session_rpe ~ 0 + core_assisted + (1|location),
                                 phi ~ 0 + core_assisted + (1|location)),
                              data = prior_data_rpe,
                              phi_reg = TRUE,
                              true_bounds = c(6,20))

pp_check(model_rpe_prior, type = "dens_overlay_grouped", group = "core_assisted")


(plogis(0.39) * 14) + 6

(plogis(0.08) * 14) + 6

get_prior(model_rpe_prior)

rpe_preds <- predictions(model_rpe_prior, variables = "core_assisted") |>
  get_draws()

slopes(model_rpe_prior)

prior_data_rpe |>
  ggplot(aes(x=session_rpe)) +
  geom_histogram(binwidth = 1, alpha = 0.5) +
  facet_grid(core_assisted~., scales = "free_y") +
  scale_x_continuous(limits = c(5.5,20.5), breaks = seq(6,20)) +
  labs(x = "Session RPE (6-20)") +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  )
