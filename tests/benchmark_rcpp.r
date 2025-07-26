library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(data.table)
sourceCpp("src/rcpp_utilities.cpp")

## Recursive softmax formula, see https://rpubs.com/FJRubio/softmax.
softmax <- function(par = NULL, recursive = TRUE) {
  if (recursive == TRUE) {
    n.par <- length(par)
    par1 <- sort(par, decreasing = TRUE)
    Lk <- par1[1]
    for (k in 1:(n.par - 1)) {
      Lk <- max(par1[k + 1], Lk) + log1p(exp(-abs(par1[k + 1] - Lk)))
    }
    val <- exp(par - Lk)
  } else {
    val <- exp(par) / sum(exp(par))
    # max approcach more numerically stable but faster than this?
  }
  val
}

softmax_grouped_R_df <- function(data) {
  group <- interaction(data$sim, data$period,data$ego, drop = TRUE)
  data$changeProb <- ave(data$changeUtil, group, FUN=softmax, recursive = FALSE)
  data
}
softmax_grouped_R_arma_df <- function(data) {
  group <- interaction(data$sim, data$period,data$ego, drop = TRUE)  # group for Rcpp
  data$changeProb <- ave(data$changeUtil, group, FUN=softmax_arma)
  data
}
softmax_grouped_data_table_df <- function(data){
  dt <- as.data.table(data)
  dt[, changeProb := {
      softmax(changeUtil, recursive = FALSE)
    }, by = .(sim, period, ego)]
  as.data.frame(dt)
}

softmax_grouped_data_table_arma_df <- function(data){
  dt <- as.data.table(data)
  dt[, changeProb := {
      softmax_arma(changeUtil)
    }, by = .(sim, period, ego)]
  as.data.frame(dt)
}

softmax_arma_by_group_df <- function(data){
  group <- interaction(data$sim, data$period,data$ego, drop = TRUE)  # group for Rcpp
  group_int <- as.integer(group)
  data$changeProb <- softmax_arma_by_group(data$changeUtil,group_int)
  data
}

n_sim <- 1
n_period <- 3
n_ego <- 100
n_choice <- 100
N <- n_sim * n_period * n_ego * n_choice
set.seed(100)
period <- rep(rep(1:n_period, each = n_ego * n_choice), times = n_sim)
ego    <- rep(rep(1:n_ego,   each = n_choice), times = n_sim * n_period)
choice <- rep(1:n_choice, times = n_sim * n_period * n_ego)
changeUtil <- rnorm(N)
sim <- rep(1:n_sim, each = N / n_sim)
df <- data.frame(sim = sim, period = period, ego = ego, choice = choice, changeUtil = changeUtil)

mb1 <- microbenchmark(
  softmax_grouped_R_df         = softmax_grouped_R_df(df),
  softmax_grouped_R_arma_df  = softmax_grouped_R_arma_df(df),
  softmax_arma_by_group_df     = softmax_arma_by_group_df(df),
  softmax_grouped_data_table_df = softmax_grouped_data_table_df(df),
  softmax_grouped_data_table_arma_df = softmax_grouped_data_table_arma_df(df),
  times = 5
)


## Check for identical results

dt_softmax_r      <- softmax_grouped_R_df(df)
dt_softmax_r_arma <- softmax_grouped_R_arma_df(df)
dt_softmax_dt     <- softmax_grouped_data_table_df(df)
dt_softmax_dt_arma<- softmax_grouped_data_table_arma_df(df)
dt_softmax_cpp    <- softmax_arma_by_group_df(df)

out_list <- list(
  softmax_grouped_R_df          = dt_softmax_r,
  softmax_grouped_R_arma_df     = dt_softmax_r_arma,
  softmax_grouped_data_table_df = dt_softmax_dt,
  softmax_grouped_data_table_arma_df = dt_softmax_dt_arma,
  softmax_arma_by_group_df      = dt_softmax_cpp
)

print(mb1, signif=4)
print(summary(mb1)[,c("expr","mean")], row.names=FALSE)

compare_changeProb <- function(df_list, tol = 1e-10) {
  combos <- combn(names(df_list), 2, simplify = FALSE)
  for (pair in combos) {
    a <- df_list[[pair[1]]][, "changeProb"]
    b <- df_list[[pair[2]]][, "changeProb"]
    max_diff <- max(abs(a - b), na.rm = TRUE)
    pass <- max_diff < tol
    cat(sprintf("%-35s vs %-35s : max abs diff = %.2e : %s\n",
                pair[1], pair[2], max_diff,
                if (pass) "PASS" else "FAIL"))
  }
}

compare_changeProb(out_list, tol = 1e-10)

## Larger data with nsim 200 for example

n_sim <- 200
n_period <- 3
n_ego <- 100
n_choice <- 100
N <- n_sim * n_period * n_ego * n_choice
set.seed(100)
period <- rep(rep(1:n_period, each = n_ego * n_choice), times = n_sim)
ego    <- rep(rep(1:n_ego,   each = n_choice), times = n_sim * n_period)
choice <- rep(1:n_choice, times = n_sim * n_period * n_ego)
changeUtil <- rnorm(N)
sim <- rep(1:n_sim, each = N / n_sim)
df <- data.frame(sim = sim, period = period, ego = ego, choice = choice, changeUtil = changeUtil)

mb2 <- microbenchmark(
  softmax_grouped_R_df         = softmax_grouped_R_df(df),
  softmax_grouped_R_arma_df  = softmax_grouped_R_arma_df(df),
  softmax_arma_by_group_df     = softmax_arma_by_group_df(df),
  softmax_grouped_data_table_df = softmax_grouped_data_table_df(df),
  softmax_grouped_data_table_arma_df = softmax_grouped_data_table_arma_df(df),
  times = 5
)

print(mb2, signif=4)
print(summary(mb1)[,c("expr","mean")], row.names=FALSE)

dt_softmax_r      <- softmax_grouped_R_df(df)
dt_softmax_r_arma <- softmax_grouped_R_arma_df(df)
dt_softmax_dt     <- softmax_grouped_data_table_df(df)
dt_softmax_dt_arma<- softmax_grouped_data_table_arma_df(df)
dt_softmax_cpp    <- softmax_arma_by_group_df(df)

out_list <- list(
  softmax_grouped_R_df          = dt_softmax_r,
  softmax_grouped_R_arma_df     = dt_softmax_r_arma,
  softmax_grouped_data_table_df = dt_softmax_dt,
  softmax_grouped_data_table_arma_df = dt_softmax_dt_arma,
  softmax_arma_by_group_df      = dt_softmax_cpp
)

compare_changeProb <- function(df_list, tol = 1e-10) {
  combos <- combn(names(df_list), 2, simplify = FALSE)
  for (pair in combos) {
    a <- df_list[[pair[1]]][, "changeProb"]
    b <- df_list[[pair[2]]][, "changeProb"]
    max_diff <- max(abs(a - b), na.rm = TRUE)
    pass <- max_diff < tol
    cat(sprintf("%-35s vs %-35s : max abs diff = %.2e : %s\n",
                pair[1], pair[2], max_diff,
                if (pass) "PASS" else "FAIL"))
  }
}

compare_changeProb(out_list, tol = 1e-10)

