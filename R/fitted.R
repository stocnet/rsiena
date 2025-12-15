fitted.sienaFit <- function(object, ...) {
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