# load packages
if (!require(pacman)) install.packages("pacman")
pacman::p_load(here, tidyverse, cmdstanr, loo, bayesplot, tidybayes, posterior, MetBrewer, patchwork)
options(mc.cores = 10)
source(here("figs/my-theme.R"))
theme_set(my_theme())

# choose paired or full acoustic data and compile model
model <- "pair" # "full"
mod <- cmdstan_model(here("analysis/mod.stan"))

# prepare data for Stan
if (model == "pair") {
  dat <- read_rds(here("data/dat.rds")) |> drop_na()
  N <- nrow(dat)
  B_max <- max(dat$bellows)
  times <- read_rds(here("data/times.rds"))[1:N, 1:(B_max + 1)]
  stan_data <- list(N = c(N, N),
                    I = nlevels(dat$site),
                    J = length(unique(dat$date)),
                    X = 4,
                    B_max = B_max,
                    site = as.numeric(dat$site),
                    survey = as.numeric(factor(dat$date_fct)),
                    C = dat$count[!is.na(dat$count)],
                    B = dat$bellows,
                    Delta = as.numeric(dat$Delta),
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
                      t()) |> 
    glimpse()
} else {
  dat <- read_rds(here("data/dat.rds"))
  times <- read_rds(here("data/times.rds"))
  stan_data <- list(N = c(length(which(!is.na(dat$count))), nrow(dat)),
                    I = nlevels(dat$site),
                    J = nlevels(dat$date_fct),
                    X = 0,
                    B_max = max(dat$bellows),
                    site = as.numeric(dat$site),
                    survey = as.numeric(dat$date_fct),
                    C = dat$count[!is.na(dat$count)],
                    B = dat$bellows,
                    Delta = as.numeric(dat$Delta),
                    y = times)
  stan_data$x <- matrix(0, stan_data$X, stan_data$N[2])
  glimpse(stan_data)
}
write_rds(stan_data, here(str_c("data/stan-data-", model, ".rds")))

# fit different configurations
configs <- expand_grid(NB1 = 0:1, NB2 = 0:1, HP = 0:1) |> 
  mutate(config = row_number())
file <- here(str_c("analysis/fits/fits-", model, ".rds"))
if (file.exists(file)) {
  fits <- read_rds(file)
} else {
  fits <- map(1:nrow(configs), \(x) {
    stan_data$NB <- configs |> select(contains("NB")) |> slice(x) |> unlist()
    stan_data$HP <- configs |> select(HP) |> slice(x) |> unlist()
    fit <- mod$sample(stan_data, chains = 10, refresh = 0, 
                      iter_warmup = 500, iter_sampling = 500,
                      adapt_delta = 0.95, show_exceptions = F)
    fit$save_object(here(str_c("analysis/fits/", model, "-", x, ".rds")))
    fit
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

# loo and model stacking (minus models with HP and NB)
loos <- map(fits, ~.$loo())
loo_compare(loos[-c(4, 8)])
loo_model_weights(loos[-c(4, 8)])

# PPC
pp_check(stan_data$C,
         fits[[2]]$draws("C_rep", format = "draws_matrix")[1:50, ],
         fun = ppc_dens_overlay)
pp_check(stan_data$B,
         fits[[2]]$draws("B_rep", format = "draws_matrix")[1:50, ],
         fun = ppc_dens_overlay)

# prepare data for figure 1
out <- fits[[2]] |> 
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
               left_join(out |> select(site, order),
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
  theme(legend.position = "bottom") + 
  labs(x = "Posterior", 
       y = "Site", 
       colour = NULL)
ggsave(here(str_c("figs/fig-mu-", model, ".png")), width = 12, height = 7, dpi = 600)

# figure 2
wrap_plots(
  fits[[2]] |> 
    spread_rvars(mu[d, n]) |> 
    left_join(dat |> 
                mutate(n = row_number()), 
              by = "n") |> 
    mutate(d = factor(d, labels = c("C", "B"))) |> 
    summarise(mu = rvar_mean(log(mu), na.rm = T), .by = c(d, region, site)) |> 
    pivot_wider(names_from = d, values_from = mu) |> 
    median_hdi(C, B) |> 
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
    labs(x = "Expected Counts (log)",
         y = "Baseline Bellow Rate (log)",
         colour = NULL) + 
    theme(legend.position = "inside",
          legend.position.inside = c(0.005, 0.995),
          legend.justification = c("left", "top")) + 
    coord_fixed(),
  map(fits[c(1, 3, 2)], ~spread_rvars(., rho[d])) |> 
    list_rbind(names_to = "model") |>
    median_hdi(rho) |> 
    ggplot(aes(x = rho, 
               y = factor(d, labels = c("Site", "Survey")) |> 
                 fct_rev())) + 
    geom_vline(xintercept = 0, 
               linetype = "dashed", 
               colour = "#333333", linewidth = 0.25) +
    geom_pointrange(aes(xmin = .lower, xmax = .upper, 
                        colour = factor(model, labels = c("Poisson", "Negative Binomial", "Hawkes Process"))),
                    shape = 16, 
                    size = 0.5, 
                    linewidth = 0.25, 
                    position = position_dodge(width = 0.2)) + 
    scale_colour_met_d("Java") + 
    scale_x_continuous(breaks = seq(-0.8, 0.8, 0.4), limits = c(-1, 1), expand = c(0, 0)) + 
    labs(x = "Correlation Coefficient", 
         y = NULL, 
         colour = NULL) + 
    theme(legend.position = "inside",
          legend.position.inside = c(0.002, 0.995),
          legend.justification = c("left", "top"))
) + 
  plot_annotation(tag_levels = "A")
ggsave(here(str_c("figs/fig-rho-", model, ".png")), width = 8, height = 5, dpi = 600)
