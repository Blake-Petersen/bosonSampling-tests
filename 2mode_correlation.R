library(BosonSampling)
library(ggplot2)

set.seed(2345)

# Huang et. al. Statistical Analysis for Collision-free Boson Sampling

correlator <- function(i, j, amat) {
  mean(amat[i, ] * amat[j, ], na.rm = TRUE)
  - (mean(amat[i, ], na.rm = TRUE) * mean(amat[j, ], na.rm = TRUE))
}

coefficient_variation <- function(c) {
  sqrt(mean(c ^ 2) - (mean(c) ^ 2)) / mean(c)
}

skewness <- function(c) {
  ( (mean(c ^ 3) - 3 * mean(c) * mean(c ^ 2) + 2 * mean(c) ^ 3)
    / sqrt( (mean(c ^ 2) - mean(c) ^ 2) ^ 3) )
}

# Change to collision-free samples in fock basis
focker_c_free <- function(bosons, modes, sample_num, amat) {
  fock <- matrix(nrow = modes, ncol = sample_num)
  for (j in 1:sample_num) {
    # Select only collision-free events
    if (!anyDuplicated(amat[, j])) {
      fock[, j] <- 0
      for (i in 1:bosons) {
        fock[amat[i, j], j] <- 1
      }
    }
  }
  fock
}

find_2mode_correlations <- function (modes, outcomes) {
  correlations <- vector(mode = "numeric", length = choose(modes, 2))
  counter <- 1
  for (j in 2:modes) {
    for (i in 1:(j - 1)) {
      correlations[counter] <- correlator(i, j, outcomes)
      counter <- counter + 1
    }
  }
  correlations
}

uniform_sampler <- function (modes, bosons, sample_num) {
  outcomes <- matrix(nrow = bosons, ncol = sample_num)
  for (i in 1:sample_num) {
    outcomes[, i] <- sample(modes, bosons)
  }
  outcomes
}

bosons <- 5L # N
modes <- 32L  # M
sample_num <- 2000L # How many samples for each network
network_num <- 200L # How many Haar random networks to sample

types <- 2L
t_boson <- 1L
t_uniform <- 2L

co_var <- matrix(nrow = network_num, ncol = types)
skew <- matrix(nrow = network_num, ncol = types)
for (nwrk in 1:network_num) {
  network_unitary <- randomUnitary(modes)
  sampling_matrix <- network_unitary[, 1:bosons]
  for (t in 1:types) {
    outcomes <- switch(
      t,
      t_boson = bosonSampler(sampling_matrix, sample_num)$values,
      t_uniform = uniform_sampler(modes, bosons, sample_num)
    )
    # Change to collision-free fock basis
    outcomes <- focker_c_free(bosons, modes, sample_num, outcomes)
    # Calculate 2-mode correlations (all of them (probably overkill?))
    correlations <- find_2mode_correlations(modes, outcomes)
    
    co_var[nwrk, t] <- coefficient_variation(correlations)
    skew[nwrk, t] <- skewness(correlations)
  }
}


stats_b <- data.frame(Type = "Boson", 
                      co_var = co_var[, t_boson],
                      skew = skew[, t_boson]
)
stats_u <- data.frame(Type = "Uniform", 
                      co_var = co_var[, t_uniform],
                      skew = skew[, t_uniform]
)
tot_stats <- rbind(stats_b, stats_u)

centroids <- aggregate(cbind(tot_stats$co_var, tot_stats$skew)~tot_stats$Type, tot_stats, mean)

colnames(centroids) <- c("Type", "co_var", "skew")

p <- ggplot(tot_stats, aes(x = co_var, y = skew, shape = Type, colour = Type))
p <- p + xlim(-1.2, 0) + ylim(-6, 4)
p <- p + stat_ellipse(alpha = 0.5)
p <- p + geom_point(size = 1, alpha = 0.75)
p <- p + geom_point(data = centroids,
                    color = "black", fill = "white",
                    shape = 21, size = 1, stroke = 0.5)
p <- p + labs(x = "Coefficient of Variation", y = "Skewness",
              title = "Statistical Analysis of Boson Sampling",
              subtitle = paste("N =", bosons, " M =", modes, sep = " "))


pdf("2mode_correlation.pdf")
print(p)
dev.off()

