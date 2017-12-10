compose_element <- function(unitary, element, m1, m2) {
  tmp <- diag(nrow(unitary))
  tmp[c(m1, m2), c(m1, m2)] <- element
  unitary <- tmp %*% unitary
  unitary
}

beamsplitter <- c(1, 1, 1, -1)

bosons <- 4L
modes <- 100L
sample_num <- 1000L
set.seed(7)

U <- diag(modes)
U <- compose_element(U, beamsplitter, 1, 2)
U <- compose_element(U, beamsplitter, 3, 4)


A <- U[, 1:bosons]

outcomes <- BosonSampling::bosonSampler(A, sample_num)$values
#outcomes <- data.frame(run = 1:sample_num, t(apply(outcomes, 2, FUN = sort)))

# Convert to Fock basis. (incorrect)
outcomes.fock <- matrix(0, ncol = sample_num, nrow = bosons)
for (j in 1:sample_num) {
  for (i in 1:bosons) {
    outcomes.fock[i, j] <- length(which(outcomes[, j] == i))
  }
}


#bins.fock <- as.matrix(partitions::compositions(bosons, modes), nrow = modes)


#outcomes.fock <- data.frame(run = 1:sample_num, outcomes.fock)
#agg <- aggregate(outcomes.fock$run, outcomes.fock[-1], FUN = length)



stars_n_stripes <- choose(bosons + modes - 1, bosons)
print(stars_n_stripes)

#bins.chars <- apply(bins.fock, 2, paste, collapse = " ")
outcomes.chars <- data.frame(run = 1:sample_num, state = apply(outcomes.fock, 2, paste, collapse = " "))
#totals <- matrix(1:stars_n_stripes)

#for (i in 1:stars_n_stripes) {
#  totals[i] <- length(which(outcomes.chars == bins.chars[i]))
#}


agg <- aggregate(outcomes.chars$run, outcomes.chars[-1], FUN = length)

barplot(agg$x / sample_num, names.arg = agg$state)

