# load packages
if (!require(pacman)) install.packages("pacman")
pacman::p_load(here, tidyverse, pipebind, cmdstanr, loo, bayesplot, tidybayes, posterior, MetBrewer, patchwork)
options(mc.cores = 10)
if (file.exists(here("figs/my-theme.R"))) {
  source(here("figs/my-theme.R"))
  theme_set(my_theme())
}

# choose paired or full acoustic data and compile model
model <- c("pair", "full")[1]
mod <- cmdstan_model(here("analysis/mod.stan"))

# prepare data for Stan
if (model == "pair") {
  dat <- read_rds(here("data/dat.rds")) |> drop_na()
  N <- nrow(dat)
  B_max <- max(dat$bellows)
  times <- read_rds(here("data/times.rds"))[1:N, 1:(B_max + 1)]
  region_i <- dat |> 
    distinct(site, region) |> 
    pull(region) |> 
    as.numeric()
  region_j <- dat |> 
    distinct(date, region) |> 
    pull(region) |> 
    as.numeric()
  stan_data <- list(N = c(N, N),
                    I = nlevels(dat$site),
                    J = length(unique(dat$date)),
                    R = nlevels(dat$region),
                    X = 4,
                    B_max = B_max,
                    site = as.numeric(dat$site),
                    region_i = region_i,
                    region_j = region_j,
                    survey = as.numeric(factor(dat$date_fct)),
                    C = dat$count[!is.na(dat$count)],
                    B = dat$bellows,
                    Delta = as.numeric(dat$Delta),
                    area = dat$area_sqm / mean(dat |> distinct(site, .keep_all = T) |> pull(area_sqm)),
                    y = times,
                    x = dat |>
                      select(region, 
                             contains(c("temp", "windspeed"))) |> 
                      mutate(region = as.numeric(region) - 1.5,
                             windspeed_start = log(windspeed_start),
                             across(contains("start"), ~scale(.)[, 1])) |> 
                      rename_with(~str_remove(., "_start")) |> 
                      mutate(temperature_sq = temperature^2, 
                             .after = temperature) |> 
                      as.matrix() |> 
                      t(),
                    pred = 1000) |> 
    glimpse()
} else {
  dat <- read_rds(here("data/dat.rds"))
  times <- read_rds(here("data/times.rds"))
  region_i <- dat |> 
    distinct(site, region) |> 
    pull(region) |> 
    as.numeric()
  region_j <- dat |> 
    distinct(date, region) |> 
    pull(region) |> 
    as.numeric()
  stan_data <- list(N = c(length(which(!is.na(dat$count))), nrow(dat)),
                    I = nlevels(dat$site),
                    J = nlevels(dat$date_fct),
                    R = nlevels(dat$region),
                    X = 0,
                    B_max = max(dat$bellows),
                    site = as.numeric(dat$site),
                    region_i = region_i,
                    region_j = region_j,
                    survey = as.numeric(dat$date_fct),
                    C = dat$count[!is.na(dat$count)],
                    B = dat$bellows,
                    Delta = as.numeric(dat$Delta),
                    area = (dat |>
                      drop_na(area_sqm) |> 
                      pull(area_sqm)) / mean(dat |> distinct(site, .keep_all = T) |> pull(area_sqm)),
                    y = times,
                    pred = 1000)
  stan_data$x <- matrix(0, stan_data$X, stan_data$N[2])
  glimpse(stan_data)
}
write_rds(stan_data, here(str_c("data/stan-data-", model, ".rds")))

# fit different configurations
configs <- expand_grid(NB1 = 0:1, NB2 = 0:1, HP = 0:1) |> 
  filter(!(NB2 == 1 & HP == 1)) |> 
  mutate(config = row_number())
file <- here(str_c("analysis/fits/fits-", model, ".rds"))
if (file.exists(file)) {
  fits <- read_rds(file)
} else {
  fits <- map(1:nrow(configs), ~{
    stan_data$NB <- configs |> select(contains("NB")) |> slice(.) |> unlist()
    stan_data$HP <- configs |> select(HP) |> slice(.) |> unlist()
    fit <- mod$sample(stan_data, chains = 10, refresh = 0, 
                      iter_warmup = 500, iter_sampling = 500,
                      adapt_delta = 0.95, show_exceptions = F)
    fit$save_object(here(str_c("analysis/fits/", model, "-", ., ".rds")))
  })
  write_rds(fits, here(str_c("analysis/fits/fits-", model, ".rds")))
}

# check diagnostics
diags <- map(fits, ~.$diagnostic_summary() |> as_tibble()) |> 
  list_rbind(names_to = "config") |> 
  left_join(configs, by = "config") |> 
  summarise(num_divergent = sum(num_divergent), 
            num_max_treedepth = sum(num_max_treedepth),
            ebfmi = mean(ebfmi), 
            .by = config)
diags

# loo and top model
loos <- map(fits, ~.$loo())
loo_compare(loos)
loo_model_weights(loos)
fit <- fits[[2]]

# prepare data for figure 1
out <- fit |> 
  spread_rvars(mu[d, n], hp[j]) |> 
  pivot_wider(names_from = j, values_from = hp) |> 
  mutate(mu = if_else(d == 2, mu / (1 - `1`/`2`), mu)) |> 
  left_join(dat |> 
              mutate(n = row_number()), 
            by = "n") |> 
  summarise(mu = exp(rvar_mean(log(mu), na.rm = T)), 
            .by = c(d, region, site)) |> 
  median_hdi(mu) |> 
  mutate(order = max(mu, na.rm = T), .by = site)

# figure 1
out |> 
  ggplot(aes(x = mu, y = fct_reorder(site, order))) + 
  facet_wrap(~ factor(d, labels = c("Drone Counts", "Bellow Rate (hourly)")),
             scales = "free_x",
             nrow = 1) + 
  geom_vline(xintercept = 0, 
             linetype = "dashed", 
             colour = "#333333", linewidth = 0.25) + 
  geom_point(aes(x = value),
             data = dat |>
               mutate(rate = bellows / as.numeric(Delta)) |>
               pivot_longer(c(count, rate), names_to = "d") |>
               left_join(out |> distinct(site, order),
                         relationship = "many-to-many"),
             colour = "#333333",
             shape = 16,
             size = 1,
             alpha = 0.6,
             position = position_jitter(width = 0, height = 0.2)) +
  geom_pointrange(aes(xmin = .lower, xmax = .upper, colour = region),
                  shape = 16, 
                  size = 0.4, 
                  linewidth = 0.25,
                  position = position_nudge(y = 0.3)) + 
  scale_y_discrete(labels = stan_data$I:1) + 
  scale_colour_met_d("Kandinsky") + 
  theme(legend.position = "inside",
        legend.position.inside = c(0.114, 0.925)) + 
  labs(x = "Posterior", 
       y = "Site", 
       colour = NULL)
ggsave(here(str_c("figs/fig-mu-", model, ".png")), width = 8, height = 5, dpi = 600)

# figure 2
wrap_plots(
  fit |>
    spread_rvars(mu[d, n]) |> 
    left_join(dat |> 
                mutate(n = row_number()), 
              by = "n") |> 
    mutate(d = factor(d, labels = c("C", "B"))) |> 
    summarise(mu = rvar_mean(log(mu), na.rm = T), .by = c(d, region, site)) |> 
    pivot_wider(names_from = d, values_from = mu) |> 
    median_hdci(C, B) |> 
    ggplot(aes(C, B)) + 
    geom_pointrange(aes(xmin = C.lower, xmax = C.upper, colour = region),
                    alpha = 0.8,
                    shape = 16, 
                    size = 0.4, 
                    linewidth = 0.25) + 
    geom_pointrange(aes(ymin = B.lower, ymax = B.upper, colour = region),
                    alpha = 0.8,
                    shape = 16, 
                    size = 0.4, 
                    linewidth = 0.25) + 
    scale_colour_met_d("Kandinsky") + 
    labs(x = "Expected Drone Counts (log)",
         y = "Baseline Bellow Rates (log)",
         colour = NULL) + 
    theme(legend.position = "inside",
          legend.position.inside = c(0.005, 0.995),
          legend.justification = c("left", "top")) + 
    coord_fixed(),
  map(fits[c(1, 3, 2)], ~spread_rvars(., rho[r, d])) |> 
    list_rbind(names_to = "model") |>
    filter(d == 1) |> 
    median_hdci(rho) |> 
    ggplot(aes(x = rho, 
               y = factor(model, labels = c("Poisson", "Negative Binomial", "Hawkes Process")))) + 
    geom_vline(xintercept = 0, 
               linetype = "dashed", 
               colour = "#333333", linewidth = 0.25) +
    geom_pointrange(aes(xmin = .lower, xmax = .upper, 
                        colour = factor(r, labels = levels(dat$region))),
                    shape = 16, 
                    size = 0.5, 
                    linewidth = 0.25, 
                    position = position_dodge(width = 0.2)) + 
    scale_colour_met_d("Kandinsky") + 
    scale_x_continuous(breaks = seq(-0.8, 0.8, 0.4), limits = c(-1, 1), expand = c(0, 0)) + 
    labs(x = "Correlation Coefficient", 
         y = NULL, 
         colour = NULL) + 
    theme(legend.position = "none"),
  fits[[5]] |> 
    spread_draws(theta_pred[i, d], theta_rep[i]) |> 
    pivot_wider(names_from = d, values_from = theta_pred) |> 
    rename(B = `1`, C = `2`) |> 
    filter(B > 0) |>
    ggplot(aes(exp(B), exp(C))) + 
    stat_lineribbon(point_interval = mean_qi, .width = 0.9,
                    show.legend = F) + 
    scale_fill_met_d("Hokusai2") + 
    scale_x_continuous(limits = c(NA, 17.5), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Predicted Drone Counts",
         y = "Predicted Baseline Bellow Rates")
) + 
  plot_annotation(tag_levels = "A")
ggsave(here(str_c("figs/fig-rho-", model, ".png")), width = 12, height = 5, dpi = 600)

# PPC
wrap_plots(
  pp_check(stan_data$C,
           fit$draws("C_rep", format = "draws_matrix")[1:50, ],
           fun = ppc_rootogram) + 
    scale_x_continuous(breaks = seq(10, 30, 10), limits = c(0, 35), expand = c(0, 0)) + 
    scale_y_continuous(breaks = seq(1, 3), limits = c(0, 3.6), expand = c(0, 0)) + 
    labs(x = "Drone Counts") +
    theme(legend.position = "none"),
  pp_check(stan_data$B,
           fit$draws("B_rep", format = "draws_matrix")[1:50, ],
           fun = ppc_rootogram) + 
    scale_x_continuous(breaks = seq(20, 60, 20), limits = c(0, 70), expand = c(0, 0)) + 
    scale_y_continuous(breaks = seq(1, 4), limits = c(0, 4.5), expand = c(0, 0)) + 
    labs(x = "Bellow Counts") +
    theme(legend.position = "inside",
          legend.position.inside = c(0.995, 0.995),
          legend.justification = c("right", "top")),
  nrow = 1
) + 
  plot_annotation(tag_levels = "A")
ggsave(here(str_c("figs/fig-ppc-", model, ".png")), width = 8, height = 5, dpi = 600)
