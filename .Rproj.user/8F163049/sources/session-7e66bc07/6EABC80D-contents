library(cluster)
library(ggplot2)

### generate hypercube clusters
generate_hypercube_clusters <- function(n, k, side_length, noise_sd = 1.0) {
  centers <- diag(side_length, n)
  out_list <- vector("list", n)
  for (i in seq_len(n)) {
    center_i <- centers[i, ]
    pts <- matrix(rnorm(k * n, mean = 0, sd = noise_sd), nrow = k, ncol = n)
    pts <- sweep(pts, 2, center_i, "+")
    out_list[[i]] <- pts
  }
  data <- do.call(rbind, out_list)
  rownames(data) <- NULL
  data
}

### kmeans

kmeans_gap_fun <- function(x, k) {
  km <- kmeans(x, centers = k, nstart = 20, iter.max = 50)
  list(cluster = km$cluster)
}

## generate shelll clusters

generate_shell_clusters <- function(n_shells, k_per_shell, max_radius, noise_sd = 0.1, inner_radius = 1) {
  
  if (max_radius < inner_radius) {
    max_radius <- inner_radius
  }
  
  radii <- seq(inner_radius, max_radius, length.out = n_shells)
  pts_list <- vector("list", n_shells)
  for (s in seq_len(n_shells)) {
    r0 <- radii[s]
    theta <- runif(k_per_shell, 0, 2 * pi)
    phi   <- acos(runif(k_per_shell, -1, 1))
    r <- r0 + rnorm(k_per_shell, 0, noise_sd)
    x <- r * sin(phi) * cos(theta)
    y <- r * sin(phi) * sin(theta)
    z <- r * cos(phi)
    pts_list[[s]] <- cbind(x, y, z)
  }
  data <- do.call(rbind, pts_list)
  rownames(data) <- NULL
  data
}

### spectral clustering

spectral_clustering_gap_fun <- function(d_threshold = 1) {
  function(x, k) {
    n <- nrow(x)
    dmat <- as.matrix(dist(x))
    A <- (dmat < d_threshold) * 1 
    diag(A) <- 0
    
    degs <- rowSums(A)
    degs[degs == 0] <- 1
    D <- diag(degs, nrow = n)
    
    L <- D - A
    D.inv.sqrt <- diag(1 / sqrt(diag(D)), nrow = n)
    L.sym <- D.inv.sqrt %*% L %*% D.inv.sqrt
    
    eig <- eigen(L.sym, symmetric = TRUE)
    idx <- order(eig$values)[1:k]
    U <- eig$vectors[, idx, drop = FALSE]
    
    row_norms <- sqrt(rowSums(U^2))
    row_norms[row_norms == 0] <- 1
    U.norm <- U / row_norms
    
    km <- kmeans(U.norm, centers = k, nstart = 20, iter.max = 50)
    list(cluster = km$cluster)
  }
}