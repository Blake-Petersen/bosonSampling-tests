library(BosonSampling)
set.seed(2345)

# Aaronson and Arkhipov's estimator for distinguishing bosonSampling
# distribution from uniform dist
estimator <- function(outcome_matrix) {
  tmp <- apply(outcome_matrix, 1, function(x) {norm(x, type = "2") ^ 2})
  n <- nrow(outcome_matrix)
  r <- prod(tmp)
  r_star <- r / (n ^ n)
  r_star
}

bosons <- 20
modes <- 500
sample_num <- 100

# Using BosonSampling's inbuilt unitary generator
sampling_matrix <- randomUnitary(modes)[, 1:bosons]

# Using unitary generator that specifically mentions Haar randomness
smatrix2 <- pracma::randortho(modes)[, 1:bosons]

rstar_d <- vector(mode = "numeric", length = sample_num)
rstar_e <- vector(mode = "numeric", length = sample_num)
rstar_u <- vector(mode = "numeric", length = sample_num)
for (i in 1:sample_num) {
  # Sample from BosonSampling distribution
  samp_d <- bosonSampler(sampling_matrix, 1)$values
  A_d <- sampling_matrix[samp_d, ]
  rstar_d[i] <- estimator(A_d)

  # Sample from BosonSampling distribution
  samp_e <- bosonSampler(smatrix2, 1)$values
  A_e <- smatrix2[samp_e, ]
  rstar_e[i] <- estimator(A_e)

  # Sample from uniform distribution
  samp_u <- sample(modes, bosons, replace = TRUE)
  A_u <- sampling_matrix[samp_u, ]
  rstar_u[i] <- estimator(A_u)
}

# Print density plot to file
pdf("row_estimator.pdf")
# Uniform distribution
plot(density(rstar_u),
     col = "red",
     main = "Comparison of uniform and BosonSampling distributions",
     sub = paste("n =",
         bosons,
         " m =",
         modes,
         " samples =",
         sample_num,
         sep = " "
         ),
     ylab = "Probability Density",
     xlab = "R*"
     )

# BosonSampling distribution, pracma::randOrtho function
lines(density(rstar_e),
      col = "blue",
      lt = 2
      )

# BosonSampling distribution, bosonSampling::randomUnitary function
lines(density(rstar_d),
      col = "green",
      lt = 4
      )

legend("topright",
  legend = c("Uniform",
    "BosonSampling,
    randOrtho",
    "BosonSampling,
    randomUnitary"
  ),
  col = c("red", "blue", "green"),
  lt = c(1, 2, 4),
  cex = 0.6
  )

dev.off()
