library(BosonSampling)
set.seed(2345)

# Implements Shchesnovich's generalised bunching test

bosons <- 4L # N in the paper
co_bin <- 2L # L
modes <- 8L  # M

bunch_bin <- modes - co_bin # K in the paper
sample_num <- 500L # How many samples for each network
network_num <- 200L # How many Haar random networks to sample

# Indistinguishable boson probabilites in [, 1], observed in [, 2]
expect <- 1
obs <- 2
results <- matrix(nrow = network_num, ncol = 2)
for (u_n in 1:network_num) {
  network_unitary <- randomUnitary(modes)

  H <- matrix(data = 0, nrow = bosons, ncol = bosons)
  for (i in 1:bosons) {
    for (j in 1:bosons) {
      for (l in 1:bunch_bin) {
        H[i, j] <- H[i, j] + (network_unitary[i,l] * Conj(network_unitary[j,l]))
      }
    }
  }
  results[u_n, expect] <- Mod(cxPerm(H))

  sampling_matrix <- network_unitary[, 1:bosons]
  outcomes <- bosonSampler(sampling_matrix, sample_num)$values > bunch_bin

  unbunched <- vector(mode = "logical", length = sample_num)
  for (i in 1:sample_num) {
    unbunched[i] <- TRUE %in% outcomes[, i]
  }
  results[u_n, obs] <- 1 - mean(unbunched)
}

p_expect <- mean(results[, expect])
p_obs <- mean(results[, obs])
p_obs_err <- sd(results[, obs])

pdf("generalized_bunching.pdf")
plot(results[, expect], results[, obs],
     main = "Observed bunching compared to indistinguishable bosons",
     sub = paste("<p_b> =", format(round(p_expect, 3), nsmall = 2),
                 " <p_obs> =", format(round(p_obs, 3), nsmall = 2),
                 sep = " "
                 ),
     xlab = "p_expect", ylab = "p_obs",
     xlim = c(0, 1),
     ylim = c(0, 1))
abline(a = 0, b = 1)
dev.off()
