library(BosonSampling)
set.seed(2345)

bosons <- 4L
test <- 2L
modes <- 8L
bunch_bin <- modes - test
sample_num <- 500L

network_unitary <- randomUnitary(modes)

H <- matrix(data = 0, nrow = bunch_bin, ncol = bunch_bin)
for (i in 1:bunch_bin) {
  for (j in 1:bunch_bin) {
    for (l in 1:bunch_bin) {
      H[i, j] <- H[i, j] + (network_unitary[i,l] * Conj(network_unitary[j,l]))
    }
  }
}
p_b <- Mod(cxPerm(H))

sampling_matrix <- network_unitary[, 1:bosons]
outcomes <- bosonSampler(sampling_matrix, sample_num)$values > bunch_bin
unbunched <- vector(mode = "logical", length = sample_num)
for (i in 1:sample_num) {
  unbunched[i] <- TRUE %in% outcomes[, i]
}
result <- 1 - mean(unbunched)

