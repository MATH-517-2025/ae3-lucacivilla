# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# Create directory for plots
dir.create("plots1")

# Set seed for reproducibility
set.seed(123)

# Generate simulation data
generate_data <- function(n = 1000, alpha = 2, beta = 2, error_sd = 1) {
  X <- rbeta(n, alpha, beta)
  m <- function(x) sin((x/3 + 0.1)^(-1))
  Y <- m(X) + rnorm(n, 0, error_sd)
  data.frame(X = X, Y = Y)
}

# Create quantile-based blocks with equal number of points
create_quantile_blocks <- function(data, N) {
  data <- data[order(data$X), ]
  quantiles <- quantile(data$X, probs = seq(0, 1, length.out = N + 1))
  data$block <- cut(data$X, breaks = quantiles, include.lowest = TRUE, labels = FALSE)
  return(data)
}

# Fit polynomial models in blocks and compute estimates
estimate_parameters <- function(data, N) {
  n <- nrow(data)
  data <- create_quantile_blocks(data, N)
  
  predictions <- numeric(n)
  second_derivs <- numeric(n)
  
  for (j in 1:N) {
    block_data <- data[data$block == j, ]
    
    block_data$X2 <- block_data$X^2
    block_data$X3 <- block_data$X^3
    block_data$X4 <- block_data$X^4
    
    model <- lm(Y ~ X + X2 + X3 + X4, data = block_data)
    
    block_indices <- which(data$block == j)
    pred_data <- data.frame(
      X = data$X[block_indices],
      X2 = data$X[block_indices]^2,
      X3 = data$X[block_indices]^3,
      X4 = data$X[block_indices]^4
    )
    
    predictions[block_indices] <- predict(model, newdata = pred_data)
    
    coefs <- coef(model)
    b2 <- coefs[3]; b3 <- coefs[4]; b4 <- coefs[5]
    second_derivs[block_indices] <- 2*b2 + 6*b3*data$X[block_indices] + 12*b4*data$X[block_indices]^2
  }
  
  theta_22 <- mean(second_derivs^2, na.rm = TRUE)
  sigma_sq <- sum((data$Y - predictions)^2, na.rm = TRUE) / (n - 5*N)
  
  list(theta_22 = theta_22, sigma_sq = sigma_sq)
}

# Compute optimal bandwidth
compute_optimal_bandwidth <- function(n, sigma_sq, theta_22, support_length = 1) {
  if (theta_22 <= 0 || is.na(theta_22)) return(NA)
  n^(-1/5) * (35 * sigma_sq * support_length / theta_22)^(1/5)
}

# Compute RSS for a given block size N
compute_rss <- function(data, N) {
  data <- create_quantile_blocks(data, N)
  rss <- 0
  
  for (j in 1:N) {
    block_data <- data[data$block == j, ]
    block_data$X2 <- block_data$X^2
    block_data$X3 <- block_data$X^3
    block_data$X4 <- block_data$X^4
    
    model <- lm(Y ~ X + X2 + X3 + X4, data = block_data)
    
    block_indices <- which(data$block == j)
    pred_data <- data.frame(
      X = data$X[block_indices],
      X2 = data$X[block_indices]^2,
      X3 = data$X[block_indices]^3,
      X4 = data$X[block_indices]^4
    )
    
    block_predictions <- predict(model, newdata = pred_data)
    rss <- rss + sum((data$Y[block_indices] - block_predictions)^2, na.rm = TRUE)
  }
  return(rss)
}

# Compute Mallow's Cp for block size selection
compute_mallows_cp <- function(data, N, N_max) {
  n <- nrow(data)
  rss_N <- compute_rss(data, N)
  rss_Nmax <- compute_rss(data, N_max)
  cp <- (rss_N / (rss_Nmax / (n - 5 * N_max))) - (n - 10 * N)
  return(cp)
}

# Find optimal block size using Mallow's Cp
find_optimal_N <- function(data) {
  n <- nrow(data)
  N_max <- max(min(floor(n / 20), 5), 1)
  N_values <- 1:N_max
  cp_values <- map_dbl(N_values, ~ compute_mallows_cp(data, ., N_max))
  optimal_N <- N_values[which.min(cp_values)]
  return(optimal_N)
}

# Define all alpha-beta combinations to test
beta_params <- list(
  c(0.5, 0.5),  # U-shaped
  c(1, 1),      # Uniform
  c(2, 2),      # Bell-shaped
  c(5, 1),      # Right-skewed
  c(1, 5)       # Left-skewed
)

# 1. Effect of sample size n for all alpha-beta combinations
n_values <- c(100, 200, 500, 1000, 2000, 5000)
sample_size_results <- data.frame()

for (n in n_values) {
  for (params in beta_params) {
    alpha <- params[1]
    beta <- params[2]
    
    data <- generate_data(n = n, alpha = alpha, beta = beta)
    
    N_optimal <- find_optimal_N(data)
    estimates <- estimate_parameters(data, N_optimal)
    h_opt <- compute_optimal_bandwidth(n, estimates$sigma_sq, estimates$theta_22)
    
    sample_size_results <- rbind(sample_size_results, data.frame(
      n = n,
      alpha = alpha,
      beta = beta,
      distribution = paste0("Beta(", alpha, ",", beta, ")"),
      N_optimal = N_optimal,
      h_AMISE = h_opt,
      theta_22 = estimates$theta_22,
      sigma_sq = estimates$sigma_sq
    ))
  }
}

# Plot 1: Effect of sample size on h_AMISE for all distributions
p1 <- ggplot(sample_size_results, aes(x = n, y = h_AMISE, color = distribution, group = distribution)) +
  geom_point(size = 2) +
  geom_line() +
  scale_x_log10() +
  labs(title = "Effect of Sample Size on Optimal Bandwidth",
       x = "Sample Size (n)",
       y = "Optimal Bandwidth (h_AMISE)",
       subtitle = "All distributions show n^(-1/5) relationship, but with different scales",
       color = "Distribution") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p1)
ggsave("plots1/effect_sample_size_all.png", p1, width = 9, height = 6, dpi = 300)

# 2. Effect of block size N for all alpha-beta combinations
N_values <- 1:10
block_size_results <- data.frame()

for (params in beta_params) {
  alpha <- params[1]
  beta <- params[2]
  
  data <- generate_data(n = 1000, alpha = alpha, beta = beta)
  
  # Find optimal N for this distribution
  N_optimal <- find_optimal_N(data)
  
  for (N in N_values) {
    estimates <- estimate_parameters(data, N)
    h_opt <- compute_optimal_bandwidth(1000, estimates$sigma_sq, estimates$theta_22)
    
    block_size_results <- rbind(block_size_results, data.frame(
      alpha = alpha,
      beta = beta,
      distribution = paste0("Beta(", alpha, ",", beta, ")"),
      N = N,
      h_AMISE = h_opt,
      theta_22 = estimates$theta_22,
      sigma_sq = estimates$sigma_sq,
      is_optimal = (N == N_optimal)
    ))
  }
}

# Plot 2: Effect of block size N on h_AMISE for all distributions
p2 <- ggplot(block_size_results, aes(x = N, y = h_AMISE, color = distribution, group = distribution)) +
  geom_point(aes(shape = is_optimal), size = 2) +
  geom_line(alpha = 0.7) +
  scale_shape_manual(values = c(16, 17), guide = "none") +
  labs(title = "Effect of Block Size on Optimal Bandwidth",
       x = "Number of Blocks (N)",
       y = "Optimal Bandwidth (h_AMISE)",
       subtitle = "Triangle points show optimal N for each distribution\nAll distributions show bias-variance tradeoff",
       color = "Distribution") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p2)
ggsave("plots1/effect_block_size_all.png", p2, width = 9, height = 6, dpi = 300)

# 3. Effect of Beta distribution parameters (already includes all)
distribution_results <- sample_size_results %>%
  filter(n == 1000) %>%
  group_by(distribution) %>%
  slice(1)  # Take one representative per distribution

# Plot 3: Effect of distribution shape on h_AMISE
p3 <- ggplot(distribution_results, aes(x = distribution, y = h_AMISE, fill = distribution)) +
  geom_col() +
  geom_text(aes(label = paste("N =", N_optimal)), vjust = -0.5) +
  labs(title = "Effect of Covariate Distribution on Optimal Bandwidth",
       x = "Beta Distribution Parameters",
       y = "Optimal Bandwidth (h_AMISE)",
       subtitle = "Numbers show optimal N from Mallow's Cp") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p3)
ggsave("plots1/effect_distribution_all.png", p3, width = 8, height = 6, dpi = 300)