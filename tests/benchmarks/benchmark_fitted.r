library(microbenchmark)
library(Matrix)
library(profvis)


## Define two versions of fitted: dense and sparse

fitted_dense <- function(object, ...) {
    if (!"sims" %in% names(object)) {
        stop("No simulated values found: run siena07 with returnDeps = TRUE.")
    }
    if (object[["f"]][["nGroup"]] > 1) {
        stop("Fitted values for multiple groups not yet implemented.")
    }
    # Group implementation might require 4-dimensional array or long array
    n_deps <- object[["nDependentVariables"]]
    if (n_deps > 1) {
        stop("Fitted values for multiple dependent variables not yet implemented.")
    }
    # Multiple dependent variables might require 4-dimensional array or long array
    depVars <- object[["f"]][[1]][["depvars"]]
    periods <- object[["periodNos"]]
    n_periods <- length(periods)
    # check attr(object[["f"]][[1]][["depvars"]][[1]], "type")
    n_chains <- length(object[["sims"]])
    # better way to get info about dependent variables?
    dep1_dims <- attr(object[["f"]][[1]][["depvars"]][[1]], "netdims")
    n_egos <- dep1_dims[1]
    n_choices <- dep1_dims[2]
    # what about missing values?
    arr <- array(0, dim = c(n_egos, n_choices, n_periods))
    # Sum over chains first
    for(c in 1:n_chains){
        # dep2 <- object$sims[[c]][[2]]
        # for oneMode etc. depvar do:
        dep1 <- object$sims[[c]][[1]][[1]]
        for (p in periods) {
            # Convert edge list to adjacency matrix
            # use more efficient methods?
            mat_new <- matrix(0, n_choices, n_choices)
            edges <- dep1[[p]]
            mat_new[edges[,1:2]] <- edges[,3]
            arr[,,p] <- arr[,,p] + mat_new
        }
    }
    # Average over chains
    arr <- arr / n_chains
    arr
}

fitted_sparse <- function(object, ...) {
    if (!"sims" %in% names(object)) {
        stop("No simulated values found: run siena07 with returnDeps = TRUE.")
    }
    if (object[["f"]][["nGroup"]] > 1) {
        stop("Fitted values for multiple groups not yet implemented.")
    }
    # Group implementation might require 4-dimensional array or long array
    n_deps <- object[["nDependentVariables"]]
    if (n_deps > 1) {
        stop("Fitted values for multiple dependent variables not yet implemented.")
    }
    # Multiple dependent variables might require 4-dimensional array or long array
    depVars <- object[["f"]][[1]][["depvars"]]
    periods <- object[["periodNos"]]
    n_periods <- length(periods)
    # check attr(object[["f"]][[1]][["depvars"]][[1]], "type")
    n_chains <- length(object[["sims"]])
    # better way to get info about dependent variables?
    dep1_dims <- attr(object[["f"]][[1]][["depvars"]][[1]], "netdims")
    n_egos <- dep1_dims[1]
    n_choices <- dep1_dims[2]
    # what about missing values?
    # arr <- array(0, dim = c(n_egos, n_choices, n_periods))
    # Sum over chains first
    # Convert into sparse matrices to save memory/time
    ans <- lapply(seq_len(n_chains), function(c) {
        dep1 <- object$sims[[c]][[1]][[1]]  # list of edge lists per period
        lapply(periods, function(p) {
            edges <- dep1[[p]]
            sparseMatrix(i = edges[,1], j = edges[,2], x = edges[,3], 
                        dims = c(n_choices, n_choices))
        })
    })
    ans <- lapply(seq_along(periods), function(ip) {
        chain_mats <- lapply(seq_len(n_chains), function(c) ans[[c]][[ip]])
        Reduce("+", chain_mats) / n_chains  # sparse + sparse → sparse
    })
    # Average over chains
    simplify2array(lapply(ans, as.matrix))
}

fitted_hybrid <- function(object) {
    if (!"sims" %in% names(object)) {
        stop("No simulated values found: run siena07 with returnDeps = TRUE.")
    }
    if (object[["f"]][["nGroup"]] > 1) {
        stop("Fitted values for multiple groups not yet implemented.")
    }
    # Group implementation might require 4-dimensional array or long array
    n_deps <- object[["nDependentVariables"]]
    if (n_deps > 1) {
        stop("Fitted values for multiple dependent variables not yet implemented.")
    }
    # Multiple dependent variables might require 4-dimensional array or long array
    depVars <- object[["f"]][[1]][["depvars"]]
    periods <- object[["periodNos"]]
    n_periods <- length(periods)
    # check attr(object[["f"]][[1]][["depvars"]][[1]], "type")
    n_chains <- length(object[["sims"]])
    # better way to get info about dependent variables?
    dep1_dims <- attr(object[["f"]][[1]][["depvars"]][[1]], "netdims")
    n_egos <- dep1_dims[1]
    n_choices <- dep1_dims[2]
  
    arr <- array(0, dim = c(n_choices, n_choices, length(periods)))
    
    for(c in 1:n_chains) {
        dep1 <- object$sims[[c]][[1]][[1]]
        for(ip in seq_along(periods)) {
        p <- periods[ip]
        edges <- dep1[[p]]
        mat_new <- matrix(0, n_choices, n_choices)
        i1 <- edges[,1]; i2 <- edges[,2]; vals <- edges[,3]
        mat_new[cbind(i1, i2)] <- vals
        arr[, , ip] <- arr[, , ip] + mat_new
        }
    }
    arr / n_chains
}


mynet <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata <- sienaDataCreate(mynet)
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, transTrip, name = "mynet")
mycontrols <- sienaAlgorithmCreate(projname = NULL, n3 = 1000, cond = FALSE)
ans <- siena07(
  mycontrols,
  data = mydata,
  effects = mymodel,
  returnDeps = TRUE,
  silent = TRUE)

block_diag_k <- function(mat, k) {
  mats <- replicate(k, mat, simplify = FALSE)
  as.matrix(bdiag(mats))
}

k <- 4  # 4 * 50 = 200 actors
s501_big <- block_diag_k(s501, k)
s502_big <- block_diag_k(s502, k)
s503_big <- block_diag_k(s503, k)

mynet_big <- sienaDependent(
  array(c(s501_big, s502_big, s503_big),
        dim = c(50 * k, 50 * k, 3))
)

mydata_big <- sienaDataCreate(mynet_big)
mymodel_big <- getEffects(mydata_big)
mymodel_big <- includeEffects(mymodel_big, transTrip, name = "mynet_big")
ans_big <- siena07(
  mycontrols,
  data = mydata_big,
  effects = mymodel_big,
  returnDeps = TRUE,
  silent = TRUE)

sparse_matrix <- function(mat, keep_prob) {
  nz <- which(mat != 0, arr.ind = TRUE)
  keep <- runif(nrow(nz)) < keep_prob
  mat_sparse <- mat
  mat_sparse[nz[!keep, , drop = FALSE]] <- 0
  mat_sparse
}

keep_prob <- 0.3  # keep ~30% of ties
s501_big_sparse <- sparse_matrix(s501_big, keep_prob)
s502_big_sparse <- sparse_matrix(s502_big, keep_prob)
s503_big_sparse <- sparse_matrix(s503_big, keep_prob)

mynet_big_sparse <- sienaDependent(
  array(c(s501_big_sparse, s502_big_sparse, s503_big_sparse),
        dim = c(50 * k, 50 * k, 3))
)

mydata_big_sparse <- sienaDataCreate(mynet_big_sparse)
mymodel_big_sparse <- getEffects(mydata_big_sparse)
mymodel_big_sparse <- includeEffects(mymodel_big_sparse, transTrip, name = "mynet_big_sparse")
ans_big_sparse <- siena07(
  mycontrols,
  data = mydata_big_sparse,
  effects = mymodel_big_sparse,
  returnDeps = TRUE,
  silent = TRUE)

# Time all versions
bm <- microbenchmark(
  dense  = fitted_dense(ans),
  sparse = fitted_sparse(ans),
  hybrid = fitted_hybrid(ans),
  dense_big       = fitted_dense(ans_big),
  sparse_big      = fitted_sparse(ans_big),
  hybrid_big      = fitted_hybrid(ans_big),
  dense_big_sparse  = fitted_dense(ans_big_sparse),
  sparse_big_sparse = fitted_sparse(ans_big_sparse),
  hybrid_big_sparse = fitted_hybrid(ans_big_sparse),
  times = 5
)
print(bm)

# Profile sparse hotspots interactively:
profvis({
  fitted_sparse(ans_big_sparse)
})

# The overhead of using sparse matrices is just to high for small - medium sized networks

# Or Rprof for precise timings:
Rprof("sparse.prof")
invisible(fitted_sparse(ans_big_sparse))
Rprof(NULL)
summaryRprof("sparse.prof")$by.self

# Stopifnot has a large self.time 
# due to repeated evaluations of the same expressions.


