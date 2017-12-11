library(BosonSampling)

bosons <- 10
modes <- 1000
sample_num <- 100
set.seed(2345)

unitary_matrix_sub <- randomUnitary(modes)[, 1:bosons]
outcomes <- matrix(bosonSampler(unitary_matrix_sub, sample_num)$values, nrow = bosons, ncol = sample_num)
#outcomes <- data.frame(run = 1:sample_num, t(apply(outcomes, 2, FUN = sort)))

# Convert to Fock basis.
#outcomes.fock <- matrix(0, nrow = sample_num, ncol = bosons)
#for (i in 1:sample_num) {
#  for (j in 1:bosons) {
#    outcomes.fock[i, j] <- length(which(outcomes[, i] == j))
#  }
#}

xvals <- matrix(1:modes)
#yvals <- matrix(apply(outcomes, 1, tabulate, nbins = modes), nrow = modes, ncol = bosons)

yvals <- apply(apply(outcomes, 1, tabulate, nbins = modes), 1, sum)
matplot(xvals, yvals, pch = 20, t = "p", 
        xlab = "output mode", ylab = "frequency")

# sample from uniform distribution
bosons <- 3
modes <- 10
x <- sample(modes, bosons, replace = TRUE)
dist_u <- vector(mode = "numeric", length = modes)
dist_u[x] <- dist_u[x] + 1
