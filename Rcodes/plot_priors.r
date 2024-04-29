# Required packages
library(ggplot2)
library(purrr)
library(dplyr)

# Function definition
f <- function(r, omega) {
  1/sqrt(2*pi)*exp(-1/2*(1/r^2 - 2*omega/r + omega^2))
}

# Vector of r values
r <- seq(1e-15, 5, length.out = 1000) 

# Vector of omega values
omegas <- c(0.01, 0.1, 0.5, 1, 5, 10, 1e3)

# Generate data frame with all combinations of r and omega
df <- expand.grid(r = r, omega = omegas) %>% 
  mutate(y = map2_dbl(r, omega, f)) # apply function f to each pair of r and omega

# Plot
(plt.likelihood <- ggplot(df, aes(x = r, y = y, colour = as.factor(omega))) + 
  geom_line() + 
  labs(x = "r", y = "f(r)", colour = "Omega") +
  theme_minimal() +
  ggtitle("Plot of likelihood L(r | omega)"))

ggsave("likelihood_plot.pdf",plt.likelihood, width = 7)
