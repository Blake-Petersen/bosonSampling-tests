library(BosonSampling)
set.seed(2345)

bosons <- 5L
modes <- 13L
sample_num <- 500L
bunch_bin <- 10L
#

network_unitary <- randomUnitary(modes)

H <- network_unitary[1:bunch_bin, 1:bunch_bin]
p_b <- Mod(cxPerm(H))

sampling_matrix <- network_unitary[, 1:bosons]
outcomes <- bosonSampler(sampling_matrix, sample_num)$values > bunch_bin

unbunched <- vector(mode = "logical", length = sample_num)
for (i in 1:sample_num) {
  unbunched[i] <- TRUE %in% outcomes[, i]
}

result <- 1 - mean(unbunched)
