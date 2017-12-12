library(BosonSampling)

compose_element <- function(unitary, element, m1, m2) {
  tmp <- diag(nrow(unitary))
  tmp[c(m1, m2), c(m1, m2)] <- element
  unitary <- tmp %*% unitary
  unitary
}

beamsplitter <- c(1, 1, 1, -1)

bosons <- 4L
modes <- 4L
sample_num <- 100L
set.seed(2345)

# Parallel beamsplitters for dual Hong-Ou-Mandel effect
U <- diag(modes)
U <- compose_element(U, beamsplitter, 1, 2)
U <- compose_element(U, beamsplitter, 3, 4)

# Can of course be any unitary eg:
# U <- randomUnitary(modes)

A <- U[, 1:bosons]
outcomes <- bosonSampler(A, sample_num)$values

# Change to fock basis and group the outcomes
# There's probably a more "R-like" way to do this, but it isn't obvious.
fock_chars <- vector(mode = "character", length = sample_num)
for (j in 1:sample_num) {
  tmp <- vector(mode = "numeric", length = modes)
  for (i in 1:bosons) {
    tmp[outcomes[i, j]] <- tmp[outcomes[i, j]] + 1
  }
  fock_chars[j] <- paste(tmp, collapse = " ")
}
outcomes_factor <- factor(fock_chars)
pdf("hong_ou_mandel.pdf")
plot(outcomes_factor,
     main = "Hong-Ou-Mandel effect with parallel beamsplitters",
     sub = paste("samples =", sample_num),
     xlab = "Output state (shows only > 0)",
     ylab = "Observations"
     )
dev.off()
