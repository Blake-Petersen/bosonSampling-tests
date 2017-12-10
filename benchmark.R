library(microbenchmark)
library(ggplot2)
library(BosonSampling)

set.seed(2345)

# n Bosons
# m Modes
sample <- function(n = 1, m = 1) {
  A <- U[, 1:n]
  bosonSampler(A, sampleSize = 1)
}

wrap <- function(b, m) {
  bquote(sample(n = .(b), m = .(m)))
}

modes <- 30L
boson_max <- 20L
boson_min <- 2L
U <- randomUnitary(modes) # This is used in the sampling function

job <- lapply(seq(boson_min, boson_max, 2), FUN = wrap, m = modes)
res <- microbenchmark(list = job)
pdf('benchmark.pdf')
autoplot(res)
dev.off()
