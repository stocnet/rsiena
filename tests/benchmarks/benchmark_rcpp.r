# testthat::skip_on_cran()

library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(data.table)
has_peakRAM <- requireNamespace("peakRAM", quietly = TRUE)
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

softmax_R_grouped_R_df <- function(data) {
  grp <- cumsum(c(TRUE, diff(data$ego) != 0 | 
                        diff(data$period) != 0 | 
                        diff(data$sim) != 0))
  data$changeProb <- ave(data$changeUtil, grp, 
                         FUN = function(x) softmax(x, recursive = FALSE))
  data
}

softmax_R_grouped_R_dt <- function(data){
  dt <- as.data.table(data)
  dt[, changeProb := {
      softmax(changeUtil, recursive = FALSE)
    }, by = .(sim, period, ego)]
  dt
}

softmax_arma_grouped_R_df <- function(data) {
  grp <- cumsum(c(TRUE, diff(data$ego) != 0 | 
                        diff(data$period) != 0 | 
                        diff(data$sim) != 0))
  data$changeProb <- ave(data$changeUtil, grp, FUN = softmax_arma)
  data
}

softmax_arma_grouped_R_dt <- function(data){
  dt <- as.data.table(data)
  dt[, changeProb := {
      softmax_arma(changeUtil)
    }, by = .(sim, period, ego)]
  dt
}

softmax_arma_by_group_df <- function(data){
  grp <- cumsum(c(TRUE, diff(data$ego) != 0 | 
                        diff(data$period) != 0 | 
                        diff(data$sim) != 0))
  data$changeProb <- softmax_arma_by_group(data$changeUtil, grp)
  data
}

softmax_arma_by_group_dt <- function(data) {
  grp <- cumsum(c(1L, diff(data$sim) != 0L | 
                      diff(data$period) != 0L | 
                      diff(data$ego) != 0L))
  set(data, j = "changeProb", value = softmax_arma_by_group(data$changeUtil, grp))
  data
}

# best df-native implementation — no data.table required
softmax_rcpp_grouped_df <- function(data) {
  data$changeProb <- softmax_rcpp_grouped(
    data$changeUtil, data$sim, data$period, data$ego)
  data
}

# dt variant with set() for completenessƒ
softmax_rcpp_grouped_dt <- function(data) {
  set(data, j = "changeProb", value = softmax_rcpp_grouped(
    data$changeUtil, data$sim, data$period, data$ego))
  data
}

softmax_rcpp_grouped_mat_df <- function(data) {
  G <- cbind(data[["sim"]], data[["period"]], data[["ego"]])
  data$changeProb <- softmax_rcpp_grouped_mat(data$changeUtil, G)
  data
}

softmax_rcpp_grouped_mat_dt <- function(data) {
  G <- cbind(data[["sim"]], data[["period"]], data[["ego"]])
  set(data, j = "changeProb", value = softmax_rcpp_grouped_mat(data$changeUtil, G))
  data
}

# List-based: zero-copy on R side (list wrapper only, no cbind allocation)
softmax_rcpp_grouped_lst_df <- function(data) {
  G <- list(data[["sim"]], data[["period"]], data[["ego"]])
  data$changeProb <- softmax_rcpp_grouped_lst(data$changeUtil, G)
  data
}

softmax_rcpp_grouped_lst_dt <- function(data) {
  G <- list(data[["sim"]], data[["period"]], data[["ego"]])
  set(data, j = "changeProb", value = softmax_rcpp_grouped_lst(data$changeUtil, G))
  data
}

# DataFrame-direct: C++ extracts column pointers itself — cleanest production API
softmax_rcpp_grouped_direct_df <- function(data) {
  data$changeProb <- softmax_rcpp_grouped_cols(data, "changeUtil", c("sim", "period", "ego"))
  data
}

softmax_rcpp_grouped_direct_dt <- function(data) {
  set(data, j = "changeProb",
      value = softmax_rcpp_grouped_cols(data, "changeUtil", c("sim", "period", "ego")))
  data
}

run_memory_profile <- function(fn_names, fns, df, dt, logfile) {
  log_benchmark(logfile, "\nMemory profile (peak RAM per call):")
  results <- lapply(fn_names, function(nm) {
    fn <- fns[[nm]]
    gc(verbose = FALSE)
    input <- if (grepl("_dt$", nm)) copy(dt) else df
    if (has_peakRAM) {
      pr <- peakRAM::peakRAM(fn(input))
      data.frame(method = nm,
                 peak_ram_mb = round(pr$Peak_RAM_Used_MiB, 2),
                 stringsAsFactors = FALSE)
    } else {
      fn(input)
      data.frame(method = nm,
                 peak_ram_mb = NA_real_,
                 stringsAsFactors = FALSE)
    }
  })
  res <- do.call(rbind, results)
  capture.output(print(res, row.names = FALSE), file = logfile, append = TRUE)
  if (!has_peakRAM)
    log_benchmark(logfile, "  (install 'peakRAM' for peak RAM measurements)")
}

log_benchmark <- function(logfile, ...) {
  cat(..., "\n", file=logfile, append=TRUE)
}
compare_changeProb_log <- function(df_list, tol = 1e-10, outfile=NULL) {
  combos <- combn(names(df_list), 2, simplify = FALSE)
    for (pair in combos) {
    a <- as.numeric(df_list[[pair[1]]]$changeProb)
    b <- as.numeric(df_list[[pair[2]]]$changeProb)
    max_diff <- max(abs(a - b), na.rm = TRUE)
    pass <- max_diff < tol
    msg <- sprintf("%-35s vs %-35s : max abs diff = %.2e : %s",
                pair[1], pair[2], max_diff,
                if (pass) "PASS" else "FAIL")
    cat(msg, "\n")
    if (!is.null(outfile)) log_benchmark(outfile, msg)
  }
}

best_opponents <- c(
  "softmax_arma_grouped_R_dt",
  "softmax_rcpp_grouped_df",
  "softmax_rcpp_grouped_dt",
  "softmax_rcpp_grouped_lst_df",
  "softmax_rcpp_grouped_lst_dt",
  "softmax_rcpp_grouped_direct_df",
  "softmax_rcpp_grouped_direct_dt"
)

run_softmax_benchmarks <- function(
  n_sim, n_period, n_ego, n_choice,
  logfile = "tests/benchmarks/results/softmax_benchmark.txt",
  compare = TRUE,
  fn_names = c(
    "softmax_R_grouped_R_df",
    "softmax_R_grouped_R_dt",
    "softmax_arma_grouped_R_df",
    "softmax_arma_grouped_R_dt",
    "softmax_arma_by_group_df",
    "softmax_arma_by_group_dt",
    "softmax_rcpp_grouped_df",
    "softmax_rcpp_grouped_dt",
    "softmax_rcpp_grouped_mat_df",
    "softmax_rcpp_grouped_mat_dt",
    "softmax_rcpp_grouped_lst_df",
    "softmax_rcpp_grouped_lst_dt",
    "softmax_rcpp_grouped_direct_df",
    "softmax_rcpp_grouped_direct_dt"
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
  dt <- as.data.table(df)
  # Prepare functions to benchmark
  fns <- mget(fn_names, mode = "function", inherits = TRUE)
  names(fns) <- fn_names
  
  expr_list <- setNames(
    lapply(fn_names, function(nm) {
      fn <- fns[[nm]]
      if (grepl("_dt$", nm)) bquote(.(fn)(dt)) else bquote(.(fn)(df))
    }),
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
  
  if (compare) {
    out_list <- setNames(
      lapply(fn_names, function(nm) {
        fn <- fns[[nm]]
        input <- if (grepl("_dt$", nm)) copy(dt) else df
        fn(input)
      }),
      nm = fn_names
    )
    log_benchmark(logfile, "\nCorrectness check (max abs. diff < 1e-10):")
    compare_changeProb_log(out_list, tol = 1e-10, outfile = logfile)
  }

  run_memory_profile(fn_names, fns, df, dt, logfile)
}

run_softmax_benchmarks(
  n_sim=1, n_period=3, n_ego=100, n_choice=100, 
  logfile='tests/benchmarks/results/softmax_benchmark_nsim1.txt'
)

run_softmax_benchmarks(
  n_sim=200, n_period=3, n_ego=100, n_choice=100, 
  logfile='tests/benchmarks/results/softmax_benchmark_nsim200.txt',
  fn_names=best_opponents, compare=FALSE
)

run_softmax_benchmarks(
  n_sim=10000, n_period=3, n_ego=50, n_choice=50,
  logfile='tests/benchmarks/results/softmax_benchmark_nsim10000.txt',
  fn_names=best_opponents, compare=FALSE
)
