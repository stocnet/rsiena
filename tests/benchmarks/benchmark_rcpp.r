# testthat::skip_on_cran()

library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(data.table)
sourceCpp("src/rcppUtilities.cpp")

## Recursive softmax formula, adjusted from https://rpubs.com/FJRubio/softmax.
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
    val <- par-max(par)
    val <- exp(val) / sum(exp(val))
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


log_benchmark <- function(logfile, ...) {
  cat(..., "\n", file=logfile, append=TRUE)
}
compare_changeProb_log <- function(df_list, tol = 1e-10, outfile=NULL) {
  combos <- combn(names(df_list), 2, simplify = FALSE)
    for (pair in combos) {
    a <- df_list[[pair[1]]][, "changeProb"]
    b <- df_list[[pair[2]]][, "changeProb"]
    max_diff <- max(abs(a - b), na.rm = TRUE)
    pass <- max_diff < tol
    msg <- sprintf("%-35s vs %-35s : max abs diff = %.2e : %s",
                pair[1], pair[2], max_diff,
                if (pass) "PASS" else "FAIL")
    cat(msg, "\n")
    if (!is.null(outfile)) log_benchmark(outfile, msg)
  }
}

run_softmax_benchmarks <- function(
  n_sim, n_period, n_ego, n_choice,
  logfile = "tests/benchmarks/softmax_benchmark.txt",
  compare = TRUE,
  fn_names = c(
    "softmax_grouped_R_df",
    "softmax_grouped_R_arma_df",
    "softmax_arma_by_group_df",
    "softmax_grouped_data_table_df",
    "softmax_grouped_data_table_arma_df"
  )
) {
  # Create directory if needed
  if (!dir.exists(dirname(logfile))) dir.create(dirname(logfile), recursive=TRUE)

  log_benchmark(logfile, sprintf("\n===== RUN: n_sim=%d, n_period=%d, n_ego=%d, n_choice=%d =====", 
                                 n_sim, n_period, n_ego, n_choice))

  N <- n_sim * n_period * n_ego * n_choice
  period <- rep(rep(1:n_period, each = n_ego * n_choice), times = n_sim)
  ego    <- rep(rep(1:n_ego,   each = n_choice), times = n_sim * n_period)
  choice <- rep(1:n_choice, times = n_sim * n_period * n_ego)
  changeUtil <- rnorm(N)
  sim <- rep(1:n_sim, each = N / n_sim)
  df <- data.frame(sim = sim, period = period, ego = ego, choice = choice, changeUtil = changeUtil)

  # Prepare functions to benchmark
  fns <- mget(fn_names, mode = "function", inherits = TRUE)
  names(fns) <- fn_names
  
  # Compose the expression list for microbenchmark dynamically
  expr_list <- setNames(
    lapply(fns, function(fn) bquote(.(fn)(df))),
    fn_names
  )
  # Now, evaluate these in microbenchmark (parse works for names; bquote/quote for actual calls)
  mb <- do.call(microbenchmark, c(expr_list, times=5))
  
  # Log benchmarks
  log_benchmark(logfile, "Microbenchmark results:")
  capture.output(print(mb, signif=4), file=logfile, append=TRUE)
  log_benchmark(logfile, "\nMean time per method (ms):")
  ms <- summary(mb)[,c("expr","mean")]
  capture.output(print(ms, row.names=FALSE), file=logfile, append=TRUE)
  
  # Run all selected implementations for correctness check
  out_list <- setNames(
    lapply(fns, function(fn) fn(df)),
    nm = names(fns)
  )
  
  if(compare){
      log_benchmark(logfile, "\nCorrectness check (max abs. diff < 1e-10):")
      compare_changeProb_log(out_list, tol = 1e-10, outfile = logfile)
  }
}

run_softmax_benchmarks(
  n_sim=1, n_period=3, n_ego=100, n_choice=100, 
  logfile='tests/benchmarks/softmax_benchmark_nsim1.txt'
)

run_softmax_benchmarks(
  n_sim=200, n_period=3, n_ego=100, n_choice=100, 
  logfile='tests/benchmarks/softmax_benchmark_nsim200.txt'
)

run_softmax_benchmarks(
  n_sim=10000, n_period=3, n_ego=50, n_choice=50,
  logfile='tests/benchmarks/softmax_benchmark_nsim10000.txt',
    fn_names=c(
    "softmax_arma_by_group_df",
    "softmax_grouped_data_table_df",
    "softmax_grouped_data_table_arma_df"
  )
)