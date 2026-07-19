library(forestBalance)
breakdown <- do.call(rbind, lapply(
  list(c(1000, 500), c(1000,1000), c(1000,2000), c(2000, 500), c(2000,1000), c(2000,2000),c(5000, 500), c(5000, 1000),c(5000,2000),
       
       c(10000, 500), c(10000, 1000),c(10000,2000), c(25000, 500), c(25000, 1000),c(25000,2000)),
  function(cfg) {
    nn <- cfg[1]; B_val <- cfg[2]
    set.seed(123)
    dat <- simulate_data(n = nn, p = 10)
    mns <- max(20L, min(floor(nn / 200) + 10, floor(nn / 50)))
    
    # 1. Forest fitting
    t_fit <- system.time({
      resp <- scale(cbind(dat$A, dat$Y))
      forest <- grf::multi_regression_forest(dat$X, Y = resp,
                                             num.trees = B_val, min.node.size = mns)
    })["elapsed"]
    
    # 2. Leaf extraction (Rcpp)
    t_leaf <- system.time(
      leaf_mat <- get_leaf_node_matrix(forest, dat$X)
    )["elapsed"]
    
    # 3. Z construction (Rcpp)
    t_z <- system.time(
      Z <- leaf_node_kernel_Z(leaf_mat)
    )["elapsed"]
    
    # 4. Weight computation (kernel_balance)
    t_bal <- system.time(
      w <- kernel_balance(dat$A, Z = Z, num.trees = B_val, solver = "cg")
    )["elapsed"]
    
    data.frame(n = nn, trees = B_val, mns = mns,
               fit = t_fit, leaf = t_leaf, Z = t_z, balance = t_bal,
               total = t_fit + t_leaf + t_z + t_bal)
  }
))


runtime_df <- tribble(
  ~n, ~trees, ~fit, ~leaf, ~kernel, ~balance, ~total,
  1000, 500, 0.049, 0.019, 0.007, 0.127, 0.202,
  1000, 1000, 0.087, 0.036, 0.010, 0.259, 0.392,
  1000, 2000, 0.175, 0.073, 0.019, 0.528, 0.795,
  2000, 500, 0.107, 0.043, 0.009, 0.268, 0.427,
  2000, 1000, 0.197, 0.087, 0.018, 0.540, 0.842,
  2000, 2000, 0.401, 0.175, 0.036, 1.154, 1.766,
  5000, 500, 0.270, 0.114, 0.023, 1.271, 1.678,
  5000, 1000, 0.538, 0.229, 0.045, 2.456, 3.268,
  5000, 2000, 1.064, 0.469, 0.092, 5.054, 6.679,
  10000, 500, 0.660, 0.236, 0.047, 5.770, 6.713,
  10000, 1000, 1.100, 0.475, 0.091, 10.249, 11.915,
  10000, 2000, 2.232, 0.974, 0.214, 19.824, 23.244,
  25000, 500, 1.552, 0.620, 0.111, 54.002, 56.285,
  25000, 1000, 3.212, 1.242, 0.255, 87.113, 91.822,
  25000, 2000, 6.642, 2.589, 0.529, 165.036, 174.796
)

runtime_long <- runtime_df %>%
  pivot_longer(
    cols = c(fit, leaf, kernel, balance),
    names_to = "stage",
    values_to = "time"
  )


ggplot(runtime_long,
       aes(x = factor(n),
           y = time,
           fill = stage)) +
  geom_col(width = 0.8) +
  facet_wrap(~ trees, nrow = 1) +
  scale_y_continuous("Runtime (seconds)") +
  scale_x_discrete("Sample size") +
  labs(fill = NULL) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )




library(ggplot2)

ggplot(runtime_long,
       aes(x = factor(n),
           y = time,
           fill = stage)) +
  geom_col(width = 0.8) +
  facet_wrap(
    ~ trees,
    nrow = 1,
    labeller = labeller(
      trees = function(x) paste("Trees =", x)
    )
  ) +
  labs(
    x = "Sample Size",
    y = "Runtime (seconds)",
    fill = "Step"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


set.seed(12345)

tree_vals <- c(200, 500, 1000,1500, 2000)
n_test <- c(1000, 2000,5000,10000,25000)

tree_bench <- do.call(rbind, lapply(n_test, function(nn) {
  do.call(rbind, lapply(tree_vals, function(B) {
    set.seed(123)
    dat <- simulate_data(n = nn, p = 10)
    t <- system.time(
      fit <- forest_balance(dat$X, dat$A, dat$Y, num.trees = B)
    )["elapsed"]
    data.frame(n = nn, trees = B, time = t)
  }))
}))

library(tidyverse)

library(tidyverse)

runtime_df <- tribble(
  ~n, ~trees, ~time,
  1000, 200, 0.060,
  1000, 500, 0.113,
  1000, 1000, 0.201,
  1000, 1500, 0.297,
  1000, 2000, 0.382,
  2000, 200, 0.157,
  2000, 500, 0.289,
  2000, 1000, 0.488,
  2000, 1500, 0.705,
  2000, 2000, 0.973,
  5000, 200, 0.930,
  5000, 500, 1.341,
  5000, 1000, 2.310,
  5000, 1500, 3.373,
  5000, 2000, 4.441,
  10000, 200, 4.957,
  10000, 500, 6.802,
  10000, 1000, 9.737,
  10000, 1500, 12.913,
  10000, 2000, 16.219,
  25000, 200, 23.952,
  25000, 500, 58.701,
  25000, 1000, 96.187,
  25000, 1500, 253.694,
  25000, 2000, 169.169
)

runtime_df$n <- factor(
  runtime_df$n,
  levels = c(1000, 2000, 5000, 10000, 25000),
  labels = c(
    "n = 1,000",
    "n = 2,000",
    "n = 5,000",
    "n = 10,000",
    "n = 25,000"
  )
)

ggplot(runtime_df,
       aes(x = trees,
           y = time,
           color = n,
           group = n)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = c(200, 500, 1000, 1500, 2000)
  ) +
  labs(
    x = "Number of Trees",
    y = "Runtime (seconds)",
    color = "Sample Size"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold")
  )




library(tidyverse)

runtime_df <- tribble(
  ~n, ~trees, ~time,
  1000, 200, 0.060,
  1000, 500, 0.113,
  1000, 1000, 0.201,
  1000, 1500, 0.297,
  1000, 2000, 0.382,
  2000, 200, 0.157,
  2000, 500, 0.289,
  2000, 1000, 0.488,
  2000, 1500, 0.705,
  2000, 2000, 0.973,
  5000, 200, 0.930,
  5000, 500, 1.341,
  5000, 1000, 2.310,
  5000, 1500, 3.373,
  5000, 2000, 4.441,
  10000, 200, 4.957,
  10000, 500, 6.802,
  10000, 1000, 9.737,
  10000, 1500, 12.913,
  10000, 2000, 16.219
)

runtime_df$n <- factor(
  runtime_df$n,
  levels = c(1000, 2000, 5000, 10000),
  labels = c(
    "n = 1,000",
    "n = 2,000",
    "n = 5,000",
    "n = 10,000"
  )
)

p <- ggplot(
  runtime_df,
  aes(
    x = trees,
    y = time,
    color = n,
    group = n
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = c(200, 500, 1000, 1500, 2000)
  ) +
  labs(
    x = "Number of Trees",
    y = "Runtime (seconds)",
    color = "Sample Size"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(linewidth = 0.8),
    axis.title = element_text(face = "bold")
  )

p

ggsave(
  "runtime_vs_trees.png",
  p,
  width = 7,
  height = 5,
  dpi = 300
)