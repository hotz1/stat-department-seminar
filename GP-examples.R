# Functions to create simple, visual examples of Gaussian processes

# Load necessary packages
library(tidyverse)
library(here)
library(latex2exp)
library(ggpubr)
library(glue)

# Setup
set.seed(5639)
img_folder <- here("img/GP-example-plots")
  
# colour scheme (muted palette from grafify)
muted_pal <- c("#CC6677", "#332288", "#DDCC77", "#117733",
               "#88CCEE", "#882255", "#44AA99", "#AA4499")

# Fix some set of hyperparameters 
ls = 0.10
sigma_eps = 0.05
sigma_f = 1

# function to calculate Gaussian Process posterior mean + variance
GP_posterior <- function(nu, prior_mean, X_new, X_obs, y_obs,
                         ls = 1, sigma_eps = 0.01, sigma_f = 1){
  # X_new, X_obs are matrices with same number of columns
  X_new <- as.matrix(X_new)
  X_obs <- as.matrix(X_obs)
  y_obs <- as.matrix(y_obs)
  new_to_new <- matern_kernel(X_new, X_new, nu, ls, sigma_f) + sigma_eps^2 * diag(nrow(X_new))
  new_to_obs <- matern_kernel(X_new, X_obs, nu, ls, sigma_f)
  obs_to_obs <- matern_kernel(X_obs, X_obs, nu, ls, sigma_f) + sigma_eps^2 * diag(nrow(X_obs))
  
  # posterior mean and variance
  post_mean <- prior_mean(X_new) + new_to_obs %*% solve(obs_to_obs, y_obs - prior_mean(X_obs))
  post_cov <- new_to_new - new_to_obs %*% solve(obs_to_obs, t(new_to_obs))
  post_sd <- sqrt(diag(post_cov))
  
  return(tibble(x = X_new[,1], post_mean = post_mean[,1], post_sd = post_sd, 
                nu = nu, ls = ls, sigma_eps = sigma_eps, sigma_f = sigma_f))
}

# prior mean function which is just zero everywhere
zero_mean <- function(X){
  matrix(0, nrow = nrow(X), ncol = 1)
}

# scaled euclidean distance
euc_distance <- function(x1, x2, ls = 1){
  sqrt(sum(((x1 - x2)/ls)^2))
}

# squared exponential function
rbf_kernel <- function(X1, X2, ls = 1, sigma_f = 1){
  # pairwise distances between rows
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  dists <- outer(1:nrow(X1), 1:nrow(X2), 
    FUN = Vectorize(function(a,b) euc_distance(X1[a,], X2[b,], ls)))
  
  return(sigma_f^2 * exp(-1/2 * dists^2))
}

# matern covariance function
matern_kernel <- function(X1, X2, nu, ls = 1, sigma_f = 1){
  # pairwise distances between rows
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  dists <- outer(1:nrow(X1), 1:nrow(X2), 
                 FUN = Vectorize(function(a,b) euc_distance(X1[a,], X2[b,], ls)))
  
  if(nu == 1/2){
    return(sigma_f^2 * exp(-dists))
  }
  else if(nu == 3/2){
    return(sigma_f^2 * (1 + sqrt(3) * dists) * exp(-sqrt(3) * dists))
  }
  else if(nu == 3/2){
    return(sigma_f^2 * (1 + sqrt(5) * dists + 5/3 * dists^2) * exp(-sqrt(5) * dists))
  }
  else{
    # Otherwise, just do RBF kernel
    return(sigma_f^2 * exp(-1/2 * dists^2))
  }
}

# true function, true maximum is at x = 0.18.
fstar <- function(t){-80 * t^5 + 102 * t^4 + 27.84 * t^3 - 73.44 * t^2 + 21.7728 * t + 0.1186}

# Plot fstar(x)
true_func_plot <- ggplot() +
  geom_function(fun = fstar, n = 501, colour = "#522a00") +
  labs(x = "x", y = TeX("$\\tilde{f}(x)$"), 
       title = TeX("Plot of the true function $y = \\tilde{f}(x)$ defined on $[0, 1]$"),
       subtitle = TeX("True maximum occurs at $x = 0.18$")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  geom_point(aes(x = 0.18, y = fstar(0.18)), fill = "#a01d2b", colour = "black", 
             shape = 24, size = 4, stroke = 1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5), 
        plot.background = element_rect(colour = NA, fill = NA))

ggsave(filename = "blackbox-function-plot.png", path = img_folder,
       width = 1800, height = 1200, units = "px")

# Priors
x_test <- matrix(seq(0, 1, by = 2e-3), ncol = 1)
GP_prior <- tibble(x = x_test[,1], post_mean = zero_mean(x_test)[,1], post_sd = 1)

# Plot prior mean
GP_prior %>% 
  mutate(lower_credible = post_mean + qnorm(0.025) * post_sd, 
         upper_credible = post_mean + qnorm(0.975) * post_sd) %>%
  ggplot() +
  geom_line(aes(x = x, y = post_mean), colour = "black", 
            linetype = "solid", linewidth = 0.75) +
  geom_line(aes(x = x, y = lower_credible), colour = "black", 
            linetype = "dashed", linewidth = 0.75) +
  geom_line(aes(x = x, y = upper_credible), colour = "black", 
            linetype = "dashed", linewidth = 0.75) +
  geom_ribbon(aes(x = x, ymin = lower_credible, ymax = upper_credible), 
              fill = muted_pal[5], alpha = 0.2) +
  labs(x = "x", y = TeX("$\\tilde{f}(x)$"), 
       title = TeX("95% credible intervals for $\\tilde{f}(x)$"),
       subtitle = TeX(paste(
         "RBF Kernel, ", "$\u2113=$", ls, ", ", "$\\sigma_{\\epsilon}=$", 
         sigma_eps, ", ", "$\\sigma_{f}=$", sigma_f, sep = ""))
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5), 
        plot.background = element_rect(colour = NA, fill = NA))

ggsave(filename = "GP-prior.png", path = img_folder,
       width = 1800, height = 1200, units = "px")

# simulate data sampling procedure
x_init <- matrix(runif(n = 10, min = 0, max = 1), ncol = 1)
y_init <- fstar(x_init)

# function to make posterior plots
plot_posterior <- function(post_df, X_obs, y_obs, ls, sigma_eps, sigma_f){
  obs_df <- as_tibble(cbind(X_obs, y_obs))
  colnames(obs_df) <- c("X_obs", "y_obs")
  
  post_df %>% 
    mutate(lower_credible = post_mean + qnorm(0.025) * post_sd, 
           upper_credible = post_mean + qnorm(0.975) * post_sd) %>%
    ggplot() +
    geom_line(aes(x = x, y = post_mean), colour = "black", 
              linetype = "solid", linewidth = 0.75) +
    geom_line(aes(x = x, y = lower_credible), colour = "black", 
              linetype = "dashed", linewidth = 0.75) +
    geom_line(aes(x = x, y = upper_credible), colour = "black", 
              linetype = "dashed", linewidth = 0.75) +
    geom_ribbon(aes(x = x, ymin = lower_credible, ymax = upper_credible), 
                fill = muted_pal[5], alpha = 0.2) +
    geom_point(aes(x = X_obs, y = y_obs), colour = muted_pal[1], 
               size = 3, stroke = 0.75, data = obs_df) +
    labs(x = "x", y = TeX("$\\tilde{f}(x)$"), 
         title = TeX("95% credible intervals for $\\tilde{f}(x)$"),
         subtitle = TeX(paste(
           "RBF Kernel, ", "$\u2113=$", ls, ", ", "$\\sigma_{\\epsilon}=$", 
           sigma_eps, ", ", "$\\sigma_{f}=$", sigma_f, sep = ""))
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5), 
          plot.background = element_rect(colour = NA, fill = NA))
}

for(i in c(1, 2, 3, 4, 5)){
  GP_post <- GP_posterior(
    nu = Inf, zero_mean, x_test, x_init[1:i,], y_init[1:i],
    ls = ls, sigma_eps = sigma_eps, sigma_f = sigma_f)
  
  plot_posterior(GP_post, x_init[1:i,], y_init[1:i,],
                 ls = ls, sigma_eps = sigma_eps, sigma_f = sigma_f)
  
  ggsave(filename = glue("GP-posterior-{i}.png"), path = img_folder,
         width = 1800, height = 1200, units = "px")
}

obs_df <- as_tibble(cbind(x_init[1:5,], y_init[1:5,]))
colnames(obs_df) <- c("X_obs", "y_obs")

# Compare posterior distributions with different kernels
kernels_post <- rbind(
  GP_posterior(nu = 1/2, zero_mean, x_test, x_init[4:5,], y_init[4:5,],
               ls = ls, sigma_eps = sigma_eps, sigma_f = sigma_f),
  GP_posterior(nu = 3/2, zero_mean, x_test, x_init[4:5,], y_init[4:5,],
               ls = ls, sigma_eps = sigma_eps, sigma_f = sigma_f),
  GP_posterior(nu = 5/2, zero_mean, x_test, x_init[4:5,], y_init[4:5,],
               ls = ls, sigma_eps = sigma_eps, sigma_f = sigma_f),
  GP_posterior(nu = Inf, zero_mean, x_test, x_init[4:5,], y_init[4:5,],
               ls = ls, sigma_eps = sigma_eps, sigma_f = sigma_f)
)
kernels_post <- kernels_post %>%
  mutate(plotname = case_when(
    nu %in% c(1/2, 3/2, 5/2) ~ glue("Mat\u00E9rn({2 * nu}/2) Kernel"),
    .default = "RBF Kernel"))
kernels_post %>% 
  mutate(lower_credible = post_mean + qnorm(0.025) * post_sd, 
         upper_credible = post_mean + qnorm(0.975) * post_sd) %>%
  ggplot() +
  geom_line(aes(x = x, y = post_mean), colour = "black", 
            linetype = "solid", linewidth = 0.75) +
  geom_line(aes(x = x, y = lower_credible), colour = "black", 
            linetype = "dashed", linewidth = 0.75) +
  geom_line(aes(x = x, y = upper_credible), colour = "black", 
            linetype = "dashed", linewidth = 0.75) +
  geom_ribbon(aes(x = x, ymin = lower_credible, ymax = upper_credible), 
              fill = muted_pal[5], alpha = 0.2) +
  geom_point(aes(x = X_obs, y = y_obs), colour = muted_pal[1], 
             size = 3, stroke = 0.75, data = obs_df[4:5,]) +
  labs(x = "x", y = TeX("$\\tilde{f}(x)$"), 
       title = TeX("95% credible intervals for $\\tilde{f}(x)$ under different distance-based kernels"),
       subtitle = TeX(paste(
         "$\u2113=$", ls, ", ", "$\\sigma_{\\epsilon}=$", 
         sigma_eps, ", ", "$\\sigma_{f}=$", sigma_f, sep = ""))
  ) + 
  facet_wrap(~plotname, nrow = 1, ncol = 4) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5), 
        plot.background = element_rect(colour = NA, fill = NA),
        strip.background = element_rect(fill = "white"))

ggsave(filename = glue("GP-diff-kernels-wide.png"), path = img_folder,
       width = 4800, height = 1200, units = "px")


# Compare posterior distributions with different lengthscales
lenscale_post <- rbind(
  GP_posterior(nu = Inf, zero_mean, x_test, x_init[1:5,], y_init[1:5,],
               ls = 0.01, sigma_eps = sigma_eps, sigma_f = sigma_f),
  GP_posterior(nu = Inf, zero_mean, x_test, x_init[1:5,], y_init[1:5,],
               ls = 0.05, sigma_eps = sigma_eps, sigma_f = sigma_f),
  GP_posterior(nu = Inf, zero_mean, x_test, x_init[1:5,], y_init[1:5,],
               ls = 0.1, sigma_eps = sigma_eps, sigma_f = sigma_f),
  GP_posterior(nu = Inf, zero_mean, x_test, x_init[1:5,], y_init[1:5,],
               ls = 0.25, sigma_eps = sigma_eps, sigma_f = sigma_f),
  GP_posterior(nu = Inf, zero_mean, x_test, x_init[1:5,], y_init[1:5,],
               ls = 0.5, sigma_eps = sigma_eps, sigma_f = sigma_f),
  GP_posterior(nu = Inf, zero_mean, x_test, x_init[1:5,], y_init[1:5,],
               ls = 1, sigma_eps = sigma_eps, sigma_f = sigma_f)
)
lenscale_post <- lenscale_post %>%
  mutate(plotname = glue("\u2113 = {ls}"))
lenscale_post %>% 
  mutate(lower_credible = post_mean + qnorm(0.025) * post_sd, 
         upper_credible = post_mean + qnorm(0.975) * post_sd) %>%
  ggplot() +
  geom_line(aes(x = x, y = post_mean), colour = "black", 
            linetype = "solid", linewidth = 0.75) +
  geom_line(aes(x = x, y = lower_credible), colour = "black", 
            linetype = "dashed", linewidth = 0.75) +
  geom_line(aes(x = x, y = upper_credible), colour = "black", 
            linetype = "dashed", linewidth = 0.75) +
  geom_ribbon(aes(x = x, ymin = lower_credible, ymax = upper_credible), 
              fill = muted_pal[5], alpha = 0.2) +
  geom_point(aes(x = X_obs, y = y_obs), colour = muted_pal[1], 
             size = 3, stroke = 0.75, data = obs_df) +
  labs(x = "x", y = TeX("$\\tilde{f}(x)$"), 
       title = TeX("95% credible intervals for $\\tilde{f}(x)$ under different lengthscales"),
       subtitle = TeX(paste(
         "RBF Kernel, ", "$\\sigma_{\\epsilon}=$", 
         sigma_eps, ", ", "$\\sigma_{f}=$", sigma_f, sep = ""))
  ) + 
  facet_wrap(~plotname, nrow = 2, ncol = 3) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5), 
        plot.background = element_rect(colour = NA, fill = NA),
        strip.background = element_rect(fill = "white"))

ggsave(filename = glue("GP-diff-lenscales.png"), path = img_folder,
       width = 2400, height = 1600, units = "px")


# Compare posterior distributions with different lengthscales
sigma_eps_post <- rbind(
  GP_posterior(nu = Inf, zero_mean, x_test, x_init[1:5,], y_init[1:5,],
               ls = ls, sigma_eps = 0.01, sigma_f = sigma_f),
  GP_posterior(nu = Inf, zero_mean, x_test, x_init[1:5,], y_init[1:5,],
               ls = ls, sigma_eps = 0.05, sigma_f = sigma_f),
  GP_posterior(nu = Inf, zero_mean, x_test, x_init[1:5,], y_init[1:5,],
               ls = ls, sigma_eps = 0.25, sigma_f = sigma_f),
  GP_posterior(nu = Inf, zero_mean, x_test, x_init[1:5,], y_init[1:5,],
               ls = ls, sigma_eps = 1, sigma_f = sigma_f)
)
sigma_eps_post %>% 
  mutate(lower_credible = post_mean + qnorm(0.025) * post_sd, 
         upper_credible = post_mean + qnorm(0.975) * post_sd) %>%
  ggplot() +
  geom_line(aes(x = x, y = post_mean), colour = "black", 
            linetype = "solid", linewidth = 0.75) +
  geom_line(aes(x = x, y = lower_credible), colour = "black", 
            linetype = "dashed", linewidth = 0.75) +
  geom_line(aes(x = x, y = upper_credible), colour = "black", 
            linetype = "dashed", linewidth = 0.75) +
  geom_ribbon(aes(x = x, ymin = lower_credible, ymax = upper_credible), 
              fill = muted_pal[5], alpha = 0.2) +
  geom_point(aes(x = X_obs, y = y_obs), colour = muted_pal[1], 
             size = 3, stroke = 0.75, data = obs_df) +
  labs(x = "x", y = TeX("$\\tilde{f}(x)$"), 
       title = TeX("95% credible intervals for $\\tilde{f}(x)$ under different observation noises"),
       subtitle = TeX(paste(
         "RBF Kernel, ", "$\u2113=$", ls, ", ", 
         "$\\sigma_{f}=$", sigma_f, sep = ""))
  ) + 
  facet_wrap(~sigma_eps, nrow = 2, ncol = 2, 
             labeller = label_bquote(sigma[epsilon]~"="~.(sigma_eps))) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  scale_y_continuous(limits = c(-3.5, 3.5), breaks = seq(-3, 3, by = 1)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5), 
        plot.background = element_rect(colour = NA, fill = NA),
        strip.background = element_rect(fill = "white"))

ggsave(filename = glue("GP-diff-noises.png"), path = img_folder,
       width = 2400, height = 1600, units = "px")


# Compare posterior distributions with different outputscales
sigma_f_post <- rbind(
  GP_posterior(nu = Inf, zero_mean, x_test, x_init[1:5,], y_init[1:5,],
               ls = ls, sigma_eps = sigma_eps, sigma_f = 0.05),
  GP_posterior(nu = Inf, zero_mean, x_test, x_init[1:5,], y_init[1:5,],
               ls = ls, sigma_eps = sigma_eps, sigma_f = 0.25),
  GP_posterior(nu = Inf, zero_mean, x_test, x_init[1:5,], y_init[1:5,],
               ls = ls, sigma_eps = sigma_eps, sigma_f = 1),
  GP_posterior(nu = Inf, zero_mean, x_test, x_init[1:5,], y_init[1:5,],
               ls = ls, sigma_eps = sigma_eps, sigma_f = 2)
)
sigma_f_post %>% 
  mutate(lower_credible = post_mean + qnorm(0.025) * post_sd, 
         upper_credible = post_mean + qnorm(0.975) * post_sd) %>%
  ggplot() +
  geom_line(aes(x = x, y = post_mean), colour = "black", 
            linetype = "solid", linewidth = 0.75) +
  geom_line(aes(x = x, y = lower_credible), colour = "black", 
            linetype = "dashed", linewidth = 0.75) +
  geom_line(aes(x = x, y = upper_credible), colour = "black", 
            linetype = "dashed", linewidth = 0.75) +
  geom_ribbon(aes(x = x, ymin = lower_credible, ymax = upper_credible), 
              fill = muted_pal[5], alpha = 0.2) +
  geom_point(aes(x = X_obs, y = y_obs), colour = muted_pal[1], 
             size = 3, stroke = 0.75, data = obs_df) +
  labs(x = "x", y = TeX("$\\tilde{f}(x)$"), 
       title = TeX("95% credible intervals for $\\tilde{f}(x)$ under different output scales"),
       subtitle = TeX(paste(
         "RBF Kernel, ", "$\u2113=$", ls, ", ",
         "$\\sigma_{\\epsilon}=$", sigma_eps, sep = ""))
  ) + 
  facet_wrap(~sigma_f, nrow = 2, ncol = 2, 
             labeller = label_bquote(sigma[f]~"="~.(sigma_f))) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  scale_y_continuous(limits = c(-5, 5), breaks = seq(-4, 4, by = 2)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5), 
        plot.background = element_rect(colour = NA, fill = NA),
        strip.background = element_rect(fill = "white"))

ggsave(filename = glue("GP-diff-outputscales.png"), path = img_folder,
       width = 2400, height = 1600, units = "px")