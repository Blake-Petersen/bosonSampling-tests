library(BosonSampling)
library(ggplot2)
set.seed(2345)

###############################################################################
# Calculate random network statistics first
#==============================================================================
# FUNCTIONS

correlator <- function(i, j, amat) {
  mean(amat[i, ] * amat[j, ], na.rm = TRUE)
  - (mean(amat[i, ], na.rm = TRUE) * mean(amat[j, ], na.rm = TRUE))
}

coefficient_variation <- function(c) {
  sd(c) / mean(c)
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

# Calculate 2-mode correlations (all of them (probably overkill?))
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

#==============================================================================
# "MAIN"

bosons <- 4L # N
modes <- 8L # M
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
    
    outcomes <- focker_c_free(bosons, modes, sample_num, outcomes)
    correlations <- find_2mode_correlations(modes, outcomes)
    
    co_var[nwrk, t] <- coefficient_variation(correlations)
    skew[nwrk, t] <- skewness(correlations)
  }
}

#==============================================================================
# OUTPUT

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

###############################################################################
# Calculate statistics for constant-depth circuit
#==============================================================================
# FUNCTIONS

correlator <- function(i, j, amat) {
  mean(amat[i, ] * amat[j, ], na.rm = TRUE)
  - (mean(amat[i, ], na.rm = TRUE) * mean(amat[j, ], na.rm = TRUE))
}

coefficient_variation <- function(c) {
  sd(c) / mean(c)
}

skewness <- function(c) {
  ( (mean(c ^ 3) - 3 * mean(c) * mean(c ^ 2) + 2 * mean(c) ^ 3)
    / sqrt( (mean(c ^ 2) - mean(c) ^ 2) ^ 3) )
}


# Calculate 2-mode correlations (all of them)
find_2mode_correlations <- function (modes, outcomes) {
  correlations <- vector(mode = "numeric", length = choose(modes - 4, 2))
  counter <- 1
  for (j in 6:modes) {
    for (i in 5:(j - 1)) {
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

compose_element <- function(unitary, element, ...) {
  tmp <- diag(nrow(unitary))
  tmp[c(...), c(...)] <- element
  unitary <- tmp %*% unitary
  unitary
}

# Change to fock basis
focker_pselect <- function(bosons, modes, sample_num, amat) {
  fock <- matrix(nrow = modes, ncol = sample_num)
  for (j in 1:sample_num) {
    # Postselect on successful teleportation
    if ( (3 %in% amat[, j]) && (4 %in% amat[, j]) && !(1 %in% amat[, j]) && !(2 %in% amat[, j]) && !anyDuplicated(amat[, j])) {
      fock[, j] <- 0
      for (i in 1:bosons) {
        fock[amat[i, j], j] <- 1
      }
    }
  }
  fock
}

#==============================================================================
# ELEMENTS

# Matrices fill colums first, so these are transpose of normal ordering

gate_beamsplitter <- c(1, 1, 
                       1, -1)

gate_phase <- function(t) {
  tmp <- c(1, 0, 
    0, exp(complex(imaginary = t)))
  tmp
}

gate_swap <- c(0, 1,
               1, 0)

gate_z <- c(1, 0,
            0, -1)

U <- diag(nrow = 8)
# Input set
U <- compose_element(U, gate_swap, 1, 6)
U <- compose_element(U, gate_swap, 3, 8)

# Initial beamsplitters
U <- compose_element(U, gate_beamsplitter, 1, 2)
U <- compose_element(U, gate_beamsplitter, 3, 4)
U <- compose_element(U, gate_beamsplitter, 5, 6)
U <- compose_element(U, gate_beamsplitter, 7, 8)

# CZ?????
U <- compose_element(U, gate_z, 1, 7)
U <- compose_element(U, gate_z, 4, 5)
U <- compose_element(U, gate_z, 2, 8)
# ??????

# Phase ???
U <- compose_element(U, gate_phase(pi), 4, 1)
U <- compose_element(U, gate_phase(pi), 5, 8)

# Measure
U <- compose_element(U, gate_beamsplitter, 1, 4)
U <- compose_element(U, gate_beamsplitter, 2, 3)
U <- compose_element(U, gate_beamsplitter, 5, 8)
U <- compose_element(U, gate_beamsplitter, 6, 7)

# Output
U <- compose_element(U, gate_swap, 6, 1)
U <- compose_element(U, gate_swap, 7, 4)

#post-select on 1,2 = 0, 3,4 = 1

# Seems to be only 1 possible output. Is this correct?
#==============================================================================
# "MAIN"



bosons <- 4L # N
modes <- 8L # M
sample_num <- 2000L # How many samples for each network
network_num <- 1L # How many Haar random networks to sample

types <- 1L
t_boson <- 1L
t_uniform <- 2L

co_var <- matrix(nrow = network_num, ncol = types)
skew <- matrix(nrow = network_num, ncol = types)
for (nwrk in 1:network_num) {
  network_unitary <- U
  sampling_matrix <- network_unitary[, 1:bosons]
  for (t in 1:types) {
    outcomes <- bosonSampler(sampling_matrix, sample_num)$values
    outcomes <- focker_pselect(bosons, modes, sample_num, outcomes)
    correlations <- find_2mode_correlations(modes, outcomes)
    
    co_var[nwrk, t] <- coefficient_variation(correlations)
    skew[nwrk, t] <- skewness(correlations)
  }
}

c_depth <- data.frame(Type = "Const Depth", 
                      co_var = co_var[, t_boson],
                      skew = skew[, t_boson]
)



#==============================================================================
# OUTPUT



p <- ggplot(tot_stats, aes(x = co_var, y = skew, shape = Type, colour = Type))
p <- p + xlim(-2.5, 0) + ylim(-2, 2)
p <- p + stat_ellipse(alpha = 0.5)
p <- p + geom_point(size = 1, alpha = 1)
p <- p + geom_point(data = centroids,
                    color = "black", fill = "white",
                    shape = 21, size = 1, stroke = 0.5)
p <- p + geom_point(data = c_depth, aes(shape = Type, colour = Type), size = 2)
p <- p + labs(x = "Coefficient of Variation", y = "Skewness",
              title = "Statistics for Constant-Depth BosonSampling with Postselection",
              subtitle = paste("N =", bosons, " M =", modes, sep = " "))


pdf("constant_depth.pdf")
print(p)
dev.off()
