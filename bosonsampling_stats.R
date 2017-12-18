library(BosonSampling)
library(ggplot2)

set.seed(2345)

# Huang et. al. Statistical Analysis for Collision-free Boson Sampling

correlator <- function(i, j, amat) {
  mean(amat[i, ] * amat[j, ], na.rm = TRUE) 
  - (mean(amat[i, ], na.rm = TRUE) * mean(amat[j, ], na.rm = TRUE))
}

coefficient_variation <- function(c) {
  sqrt(mean(c^2) - (mean(c)^2)) / mean(c)
}

skewness <- function(c) {
  ((mean(c^3) - 3 * mean(c) * mean(c^2) + 2 * mean(c)^3) 
    / sqrt((mean(c^2) - mean(c)^2)^3))
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
    for (i in 1:(j-1)) {
      correlations[counter] <- correlator(i, j, outcomes_fock)
      counter <- counter + 1
    }
  }
  correlations
}

bosons <- 5L # N
modes <- 32L  # M
sample_num <- 2000L # How many samples for each network
network_num <- 200L # How many Haar random networks to sample

co_var <- vector(mode = "numeric", length = network_num)
skew <- vector(mode = "numeric", length = network_num)
for (nwrk in 1:network_num) {
  network_unitary <- randomUnitary(modes)
  sampling_matrix <- network_unitary[, 1:bosons]
  
  outcomes <- bosonSampler(sampling_matrix, sample_num)$values
  
  # Change to collision-free fock basis
  outcomes_fock <- focker_c_free(bosons, modes, sample_num, outcomes)
  
  # Calculate 2-mode correlations (all of them (probably overkill?))
  correlations <- find_2mode_correlations(modes, outcomes_fock)
  
  co_var[nwrk] <- coefficient_variation(correlations)
  skew[nwrk] <- skewness(correlations)
}

pdf("correlation_stats.pdf")
res <- data.frame(co_var, skew)
p <- ggplot(res, aes(x = co_var, y = skew)) 
p <- p + xlim(-1.2, 0) + ylim(-6,4)
#p <- p + stat_density_2d()
p <- p + geom_point(aes(x = mean(co_var), y = mean(skew)), color = "magenta", size = 4, alpha = 1)
p <- p + geom_point(size = 1, alpha = 0.5, color = "darkslateblue") 
p <- p + labs(x = "Coefficient of Variation", y = "Skewness", 
              title = "Statistical Analysis of Boson Sampling")
p

dev.off()

