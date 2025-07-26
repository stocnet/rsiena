library(Rcpp)
library(RcppArmadillo)
library(data.table)
library(microbenchmark)
sourceCpp("src/rcpp_utilities.cpp") 

# Re-define the relevant functions for the test
softmax_grouped_data_table_arma_df <- function(data){
  dt <- as.data.table(data)
  dt[, changeProb := softmax_arma(changeUtil), by = .(sim, period, ego)]
  as.data.frame(dt)
}
softmax_arma_by_group_df <- function(data){
  group <- interaction(data$sim, data$period, data$ego,drop=TRUE)
  group_int <- as.integer(group)
  data$changeProb <- softmax_arma_by_group(data$changeUtil, group_int)
  data
}

# Large data shape
n_sim <- 50000
n_period <- 3
n_ego <- 50
n_choice <- 50
N <- n_sim * n_period * n_ego * n_choice
cat(sprintf("Simulating %d rows (~%.1f GB for doubles only)\n", N, 8 * N / 1e9))

# Simulate SIENA-style changeUtil data
set.seed(2023)
period <- rep(rep(1:n_period, each = n_ego * n_choice), times = n_sim)
ego    <- rep(rep(1:n_ego,   each = n_choice), times = n_sim * n_period)
choice <- rep(1:n_choice, times = n_sim * n_period * n_ego)
sim <- rep(1:n_sim, each = N / n_sim)
changeUtil <- rnorm(N)
df <- data.frame(sim = sim, period = period, ego = ego, choice = choice, changeUtil = changeUtil)

# BENCHMARK
mb1 <- microbenchmark(
  softmax_grouped_data_table_arma_df = softmax_grouped_data_table_arma_df(df),
  softmax_arma_by_group_df           = softmax_arma_by_group_df(df),
  times = 1
)
print(mb1)

cat("DONE.\n")
