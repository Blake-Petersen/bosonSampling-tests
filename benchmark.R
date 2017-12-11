library(microbenchmark)
library(ggplot2)
library(BosonSampling)

set.seed(2345)

# n Bosons
# m Modes
sample <- function(n = 1, m = 1) {
  A <- U[, 1:n]
  bosonSampler(A, sampleSize = 10)
}

wrap <- function(b, m) {
  bquote(sample(n = .(b), m = .(m)))
}


# Characterise time depencence on boson number
modes <- 30L
boson_max <- 20L
boson_min <- 2L
U <- randomUnitary(modes) # This is used in the sampling function

job <- lapply(seq(boson_min, boson_max, 2), FUN = wrap, m = modes)
res <- microbenchmark(list = job)
pdf('boson_benchmark.pdf')
autoplot(res)
dev.off()

# Characterise time depencence on mode number
bosons <- 4L
modes_max <- 30L
A <- lapply(seq(bosons, modes_max, 2), function(x) randomUnitary(x)[, 1:bosons])
wrap2 <- function(A) {
  bquote(bosonSampler(.(A), sampleSize = 10))
}

job2 <- lapply(A, FUN = wrap2)
res2 <- microbenchmark(list = job2)
pdf('mode_benchmark.pdf')
autoplot(res2) + theme(axis.text.y = element_blank())
dev.off()
