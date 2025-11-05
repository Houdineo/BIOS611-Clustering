library(cluster)
library(ggplot2)
source("R/functions.R")

set.seed(12345)

n_shells <- 4
k_per_shell <- 100
noise_sd <- 0.1
d_threshold <- 1        
max_radii <- 10:0

spec_fun <- spectral_clustering_gap_fun(d_threshold = d_threshold)

res_list <- list()

for (Rmax in max_radii) {
  dat <- generate_shell_clusters(
    n_shells = n_shells,
    k_per_shell = k_per_shell,
    max_radius = Rmax,
    noise_sd = noise_sd,
    inner_radius = 1
  )
  
  gap <- clusGap(
    dat,
    FUNcluster = spec_fun,
    K.max = n_shells + 3, 
    B = 20
  )
  
  est_k <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method = "firstSEmax")
  
  res_list[[length(res_list) + 1]] <- data.frame(
    max_radius = Rmax,
    estimated_k = est_k,
    true_k = n_shells
  )
  cat("Finished max_radius", Rmax, "->", est_k, "\n")
}

res_df <- do.call(rbind, res_list)

p <- ggplot(res_df, aes(x = max_radius, y = estimated_k)) +
  geom_line() +
  geom_point() +
  geom_hline(aes(yintercept = true_k), linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = max_radii) +
  labs(
    title = "Spectral clustering + gap statistic on concentric shells",
    x = "max_radius (shells closer as max_radius ↓)",
    y = "Estimated number of clusters"
  ) +
  theme_minimal()

if (!dir.exists("figs")) dir.create("figs", recursive = TRUE)
ggsave("figs/task2_gap.png", p, width = 8, height = 5, dpi = 300)
