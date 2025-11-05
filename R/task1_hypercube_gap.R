
library(cluster)
library(ggplot2)
source("R/functions.R")

set.seed(12345)

dims <- c(6, 5, 4, 3, 2)
side_lengths <- 10:1
k_per_cluster <- 100
noise_sd <- 1.0

results <- list()

for (d in dims) {
  for (L in side_lengths) {
    dat <- generate_hypercube_clusters(
      n = d,
      k = k_per_cluster,
      side_length = L,
      noise_sd = noise_sd
    )
    gap <- clusGap(
      dat,
      FUNcluster = kmeans_gap_fun,
      K.max = d + 3,
      B = 20
    )
    
    bestK <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], method = "firstSEmax")
    
    results[[length(results) + 1]] <- data.frame(
      dimension = d,
      side_length = L,
      estimated_k = bestK,
      true_k = d
    )
    cat("Finished dim", d, "L", L, "->", bestK, "\n")
  }
}

res_df <- do.call(rbind, results)

p <- ggplot(res_df, aes(x = side_length, y = estimated_k)) +
  geom_line() +
  geom_point() +
  geom_hline(aes(yintercept = true_k), linetype = "dashed", color = "red") +
  scale_x_reverse(breaks = side_lengths) +
  facet_wrap(~ dimension, scales = "free_y") +
  labs(
    title = "Gap statistic estimated number of clusters vs side length",
    x = "Side length (clusters closer → smaller L)",
    y = "Estimated number of clusters"
  ) +
  theme_minimal()

if (!dir.exists("figs")) dir.create("figs", recursive = TRUE)
ggsave("figs/task1_gap.png", p, width = 8, height = 5, dpi = 300)
