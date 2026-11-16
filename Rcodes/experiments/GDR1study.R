#################### GDR1 data analysis #############################
#####################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(scales)
library(ggdensity)
library(patchwork)

gdr1 <- read.csv("~/Downloads/gdr1set01.csv")

# Convert mas to arcsec so distances are in parsecs
gdr1 <- gdr1 %>%
  mutate(
    parallax_arcsec       = 1e-3 * parallax,
    parallax_error_arcsec = 1e-3 * parallax_error
  ) %>%
  filter(
    is.finite(parallax_arcsec),
    is.finite(parallax_error_arcsec),
    parallax_error_arcsec > 0
  )

dat <- gdr1[1:100, ]

L <- 1000  


log_prior_gamma <- function(r, L = 1000) {
  # Gamma(shape = 3, scale = L): ED prior
  dgamma(r, shape = 3, scale = L, log = TRUE)
}

log_prior_inv_gamma <- function(r, L = 1000) {
  a <- 4
  beta <- L
  a * log(beta) - lgamma(a) -
    (a + 1) * log(r) - beta / r
}

log_prior_half_cauchy <- function(r, L = 1000) {
  # Half-Cauchy(0,L): 2/[pi L {1 + (r/L)^2}]
  log(2 / pi) - log(L) - log1p((r / L)^2)
}

log_prior_rg <- function(r, L = 1000) {
  0.5 * log(2 / pi) + log(L) -
    2 * log(r) - L^2 / (2 * r^2)
}

log_prior_weibull <- function(r, L = 1000) {
  # Weibull(shape = 0.5, scale = L)
  dweibull(r, shape = 0.5, scale = L, log = TRUE)
}

log_prior_phc <- function(r, L = 1000) {
  # Product Half-Cauchy with scale L:
  # If Z = X*Y, X,Y iid Half-Cauchy(0,1), and R = LZ,
  # f_R(r) = (1/L) * 4 log(r/L) / [pi^2 {(r/L)^2 - 1}]
  z <- r / L
  out <- numeric(length(z))
  
  close <- abs(z - 1) < 1e-8
  out[close] <- log(2 / pi^2) - log(L)
  
  not_close <- !close
  dens_z <- 4 * log(z[not_close]) / (pi^2 * (z[not_close]^2 - 1))
  dens_z <- pmax(dens_z, 1e-300)
  out[not_close] <- log(dens_z) - log(L)
  
  out
}

log_priors <- list(
  "Gamma (3,1000)"                 = function(r) log_prior_gamma(r, L),
  "Reciprocal-Gaussian (0,1000)"   = function(r) log_prior_rg(r, L),
  "Weibull (0.5,1000)"             = function(r) log_prior_weibull(r, L),
  "Inverse-Gamma (4,1000)"       = function(r) log_prior_inv_gamma(r, L),
  "Half-Cauchy (0,1000)"           = function(r) log_prior_half_cauchy(r, L),
  "Product Half-Cauchy (0,1000)"   = function(r) log_prior_phc(r, L)
)

light_priors <- c(
  "Gamma (3,1000)",
  "Reciprocal-Gaussian (0,1000)",
  "Weibull (0.5,1000)"
)

heavy_priors <- c(
  "Inverse-Gamma (4,1000)",
  "Half-Cauchy (0,1000)",
  "Product Half-Cauchy (0,1000)"
)

prior_levels <- c(light_priors, heavy_priors)

my_cols <- c(
  "Gamma (3,1000)"                 = "#CB4154",
  "Reciprocal-Gaussian (0,1000)"   = "orange",
  "Weibull (0.5,1000)"             = "purple",
  "Inverse-Gamma (4,1000)"       = "#4DA3FF",
  "Half-Cauchy (0,1000)"           = "#D8B400",
  "Product Half-Cauchy (0,1000)"   = "#3A9D23"
)

my_shapes <- c(
  "Gamma (3,1000)"                 = 21,
  "Reciprocal-Gaussian (0,1000)"   = 22,
  "Weibull (0.5,1000)"             = 23,
  "Inverse-Gamma (4,1000)"       = 24,
  "Half-Cauchy (0,1000)"           = 25,
  "Product Half-Cauchy (0,1000)"   = 8
)


log_likelihood_r <- function(r, omega, sigma) {
  -0.5 * ((omega - 1 / r) / sigma)^2
}

log_posterior_r <- function(r, omega, sigma, log_prior_fun) {
  out <- rep(-Inf, length(r))
  ok <- r > 0 & is.finite(r)
  out[ok] <- log_prior_fun(r[ok]) + log_likelihood_r(r[ok], omega, sigma)
  out
}


posterior_median_r <- function(omega,
                               sigma,
                               log_prior_fun,
                               r_lower = 1e-6,
                               r_upper = 1e8,
                               n_grid = 12000) {
  
  if (!is.finite(omega) || !is.finite(sigma) || sigma <= 0) {
    return(NA_real_)
  }
  
  r_grid <- exp(seq(log(r_lower), log(r_upper), length.out = n_grid))
  
  lp <- log_posterior_r(
    r = r_grid,
    omega = omega,
    sigma = sigma,
    log_prior_fun = log_prior_fun
  )
  
  if (all(!is.finite(lp))) {
    return(NA_real_)
  }
  
  lp_max <- max(lp[is.finite(lp)])
  dens <- exp(lp - lp_max)
  dens[!is.finite(dens)] <- 0
  
  dr <- diff(r_grid)
  area_segments <- 0.5 * (dens[-1] + dens[-length(dens)]) * dr
  total_area <- sum(area_segments)
  
  if (!is.finite(total_area) || total_area <= 0) {
    return(NA_real_)
  }
  
  cdf <- c(0, cumsum(area_segments)) / total_area
  
  med <- approx(
    x = cdf,
    y = r_grid,
    xout = 0.5,
    ties = "ordered"
  )$y
  
  as.numeric(med)
}



median_list <- lapply(names(log_priors), function(pname) {
  message("Computing posterior medians for: ", pname)
  
  lp_fun <- log_priors[[pname]]
  
  r_med <- mapply(
    FUN = function(w, s) posterior_median_r(w, s, lp_fun),
    w = dat$parallax_arcsec,
    s = dat$parallax_error_arcsec
  )
  
  dat %>%
    mutate(
      Prior = pname,
      r_med = as.numeric(r_med)
    )
})

res_df <- bind_rows(median_list) %>%
  mutate(
    Prior = factor(Prior, levels = prior_levels),
    r_naive = ifelse(parallax_arcsec > 0, 1 / parallax_arcsec, NA_real_),
    f_obs = parallax_error_arcsec / abs(parallax_arcsec),
    shrinkage_fraction = 1 - r_med / r_naive,
    log_ratio = log(r_med / r_naive)
  )

plot_hdr_panel <- function(df, priors_to_show, panel_title) {
  
  df_plot <- df %>%
    filter(
      Prior %in% priors_to_show,
      parallax_arcsec > 0,
      is.finite(r_naive),
      is.finite(r_med),
      r_naive > 0,
      r_med > 0
    ) %>%
    mutate(
      Prior = factor(Prior, levels = priors_to_show)
    ) %>%
    droplevels()
  
  ggplot(
    df_plot,
    aes(
      x = r_naive,
      y = r_med,
      fill = Prior,
      shape = Prior
    )
  ) +
    geom_hdr(
      xlim = c(0, 5000),
      ylim = c(0, 5000)
    ) +
    geom_point(
      alpha = 0.55,
      size = 3
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "red",
      linetype = "dashed",
      linewidth = 0.7
    ) +
    scale_fill_manual(
      name = "Prior",
      values = my_cols[priors_to_show],
      breaks = priors_to_show,
      drop = TRUE
    ) +
    scale_shape_manual(
      name = "Prior",
      values = my_shapes[priors_to_show],
      breaks = priors_to_show,
      drop = TRUE
    ) +
    coord_cartesian(
      xlim = c(0, 5000),
      ylim = c(0, 5000),
      expand = FALSE
    ) +
    labs(
      title = panel_title,
      x = "Inverse Parallax [pc]",
      y = "Median Distance [pc]"
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        ncol = 1,
        override.aes = list(
          shape = unname(my_shapes[priors_to_show]),
          size = 3,
          alpha = 1
        )
      ),
      
      shape = "none",
      
      
      alpha = guide_legend(
        title = "HDR probability",
        order = 2,
        ncol = 1
      )
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      
      plot.title = element_text(
        size = 18,
        hjust = 0.5,
        margin = margin(b = 8)
      ),
      
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "top",
      
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      legend.key.height = grid::unit(0.55, "cm"),
      legend.key.width = grid::unit(0.75, "cm"),
      
      panel.grid.minor = element_blank(),
      
      plot.margin = margin(
        t = 10,
        r = 20,
        b = 10,
        l = 10
      )
    )
}


p_light <- plot_hdr_panel(
  res_df,
  light_priors,
  panel_title = ""
)

p_heavy <- plot_hdr_panel(
  res_df,
  heavy_priors,
  panel_title = ""
)

print(p_light)
print(p_heavy)

####################################################### Shrinkage Plots ###########################################################################

shrink_df <- res_df %>%
  filter(
    parallax_arcsec > 0,
    is.finite(r_naive),
    is.finite(r_med),
    r_naive > 0,
    r_med > 0,
    is.finite(shrinkage_fraction)
  )

shrink_summary <- shrink_df %>%
  group_by(Prior) %>%
  summarise(
    n = n(),
    median_shrinkage = median(shrinkage_fraction, na.rm = TRUE),
    mean_shrinkage = mean(shrinkage_fraction, na.rm = TRUE),
    q25_shrinkage = quantile(shrinkage_fraction, 0.25, na.rm = TRUE),
    q75_shrinkage = quantile(shrinkage_fraction, 0.75, na.rm = TRUE),
    prop_downward_shrinkage = mean(shrinkage_fraction > 0, na.rm = TRUE),
    median_log_ratio = median(log_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    median_shrinkage_pct = 100 * median_shrinkage,
    mean_shrinkage_pct = 100 * mean_shrinkage,
    prop_downward_shrinkage_pct = 100 * prop_downward_shrinkage
  ) %>%
  arrange(desc(median_shrinkage))

print(shrink_summary)

shrink_df2 <- shrink_df %>%
  filter(
    is.finite(f_obs),
    f_obs > 0,
    is.finite(shrinkage_fraction)
  ) %>%
  mutate(log_f = log10(f_obs))

breaks <- unique(
  quantile(
    shrink_df2$log_f,
    probs = seq(0, 1, length.out = 16),
    na.rm = TRUE
  )
)

shrink_curve <- shrink_df2 %>%
  mutate(
    f_bin = cut(
      log_f,
      breaks = breaks,
      include.lowest = TRUE
    )
  ) %>%
  group_by(Prior, f_bin) %>%
  summarise(
    f_mid = median(f_obs, na.rm = TRUE),
    shrink_med = median(shrinkage_fraction, na.rm = TRUE),
    shrink_q25 = quantile(shrinkage_fraction, 0.25, na.rm = TRUE),
    shrink_q75 = quantile(shrinkage_fraction, 0.75, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(
    is.finite(f_mid),
    is.finite(shrink_med)
  )
x_lower <- 0.03
x_upper <- max(shrink_curve$f_mid, na.rm = TRUE)

candidate_breaks <- c(
  0.03, 0.05, 0.10, 0.30, 0.50, 1.00, 2.00
)

x_breaks <- candidate_breaks[
  candidate_breaks >= x_lower &
    candidate_breaks <= x_upper
]

p_shrink_f <- ggplot(
  shrink_curve,
  aes(
    x = f_mid,
    y = shrink_med,
    color = Prior,
    fill = Prior
  )
) +
  geom_hline(
    yintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_ribbon(
    aes(
      ymin = shrink_q25,
      ymax = shrink_q75
    ),
    alpha = 0.12,
    color = NA
  ) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 1.7) +
  scale_x_log10(
    breaks = x_breaks,
    labels = format(x_breaks, nsmall = 2),
    expand = expansion(mult = c(0, 0.025))
  ) +
  coord_cartesian(
    xlim = c(x_lower, x_upper * 1.03)
  ) +
  scale_color_manual(
    values = my_cols,
    drop = FALSE
  ) +
  scale_fill_manual(
    values = my_cols,
    drop = FALSE
  ) +
  labs(
    x = expression(
      f == sigma[omega] / abs(omega)
    ),
    y = expression(
      "Median shrinkage fraction " *
        S == 1 - hat(r)[med] / (1 / omega)
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
print(p_shrink_f)
