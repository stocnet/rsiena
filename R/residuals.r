##@residuals.sienaFit
residuals.sienaFit <- function(object, 
    type = c("response", "pearson", "quantile"), # add other types later
    fitted = NULL,
    useSimulatedVar = FALSE, # only relevant for pearson residuals
    randomized = FALSE, # only relevant for quantile residuals
    seed = NULL, # only relevant for randomized quantile residuals
    ...) # should also allow to give simulated values directly
{
    type <- match.arg(type)
    if (object[["f"]][["nGroup"]] > 1) {
        stop("Residuals for multiple groups not yet implemented.")
    }
    n_deps <- object[["nDependentVariables"]]
    if (n_deps > 1) {
        stop("Residuals for multiple dependent variables not yet implemented.")
    }
    
    if(is.null(fitted)){
            fitted_values <- fitted(object)
        } else {
            fitted_values <- fitted
    }
    # add more safety nets
    # the following does not ALWAYS work
    observed_values <- object[["f"]][[1]][["depvars"]][[1]][,,2:3]

    if(type == "response"){
        residuals <- observed_values - fitted_values
    } else if(type == "pearson"){
        if(useSimulatedVar){
            stop("useSimulatedVar = TRUE is the same for binary networks and
                residuals are not implemented yet for non binary-networks")
            # this should work for non-binary outcomes but otherwise not make a difference
            # compute empirical variance from simulations (robust to non-binary outcomes)
            # if (!"sims" %in% names(object)) {
            #     stop("No simulated values found: run siena07 with returnDeps = TRUE.")
            # }
            # n_chains <- length(object[["sims"]])
            # n_periods <- length(object[["periodNos"]])
            # dep1_dims <- attr(object[["f"]][[1]][["depvars"]][[1]], "netdims")
            # n_egos <- dep1_dims[1]
            # n_choices <- dep1_dims[2]
            # sumsq_arr <- array(0, dim = c(n_egos, n_choices, n_periods))
            # # Collect simulated values
            # for(c in 1:n_chains){
            #     dep1 <- object$sims[[c]][[1]][[1]]
            #     for (p in object[["periodNos"]]) {
            #         mat_new <- matrix(0, n_egos, n_choices)
            #         edges <- dep1[[p]]
            #         mat_new[edges[,1:2]] <- edges[,3]
            #         sumsq_arr[,,p] <- sumsq_arr[,,p] + mat_new^2
            #     }
            # }
            # mean_sq_arr <- sumsq_arr / n_chains
            # var_arr <- mean_sq_arr - fitted_values^2  # empirical variance formula
            # residuals <- (observed_values - fitted_values) / sqrt(var_arr)
        } else {
        # this assumes that we can approximate as a Bernoulli process
            residuals <- (observed_values - fitted_values) / 
                sqrt(fitted_values * (1 - fitted_values))
        }
    }  else if(type == "quantile"){
        # add check to only allow for binary data
        if (!"sims" %in% names(object)) {
            stop("No simulated values found but required for quantile residuals:
            run siena07 with returnDeps = TRUE.")
        }

        n_periods <- length(object[["periodNos"]])
        n_chains <- length(object[["sims"]])        
        dep1_dims <- attr(object[["f"]][[1]][["depvars"]][[1]], "netdims")
        n_egos <- dep1_dims[1]
        n_choices <- dep1_dims[2]

        # Flatten observed and fitted to vectors
        observed_vec <- as.vector(observed_values)
        fitted_vec <- as.vector(fitted_values)
        n_obs <- length(observed_vec)            
        chain_mat <- matrix(0, nrow = n_obs, ncol = n_chains)
        len <- n_egos * n_choices
        for (c in 1:n_chains) {
            dep1 <- object$sims[[c]][[1]][[1]]
            sim_vec <- numeric(length = n_obs)
            idx <- 1
            for (p in object[["periodNos"]]) {
                edges <- dep1[[p]]
                mat_new <- matrix(0, n_egos, n_choices)
                if (nrow(edges) > 0) {
                    mat_new[edges[,1:2]] <- edges[,3]
                }
                sim_vec[idx:(idx + len - 1)] <- as.vector(mat_new)
                idx <- idx + len
            }
            chain_mat[, c] <- sim_vec
        }
        residuals_vec <- getQuantileResiduals(chain_mat, observed_vec, 
            integerResponse = TRUE, method = "PIT", randomized = randomized)
        residuals <- array(residuals_vec, dim = c(n_egos, n_choices, n_periods))
        # if(normalize){
        #     residuals <- qnorm(u)
        # }
    }
    return(residuals)
}

# add resid as alias for residuals
# based on https://github.com/florianhartig/DHARMa/blob/master/DHARMa/R/helper.R
#' Calculate Quantile Residuals from Simulations
#'
#' @param simulations Matrix: rows = observations, columns = simulations.
#' @param observed Vector: observed values.
#' @param integerResponse Logical: is response integer-valued?
#' @param method "PIT" or "traditional"
#' @return Vector of quantile residuals
getQuantileResiduals <- function(simulations, observed, 
    integerResponse = TRUE, method = c("PIT", "traditional"),
    randomized = TRUE) {
    method <- match.arg(method)
    n <- length(observed)
    nSim <- ncol(simulations)
    scaledResiduals <- rep(NA, n)
    if (method == "traditional") {
        for (i in 1:n) {
            if (integerResponse) {
                if (randomized) {
                    scaledResiduals[i] <- DHARMa.ecdf(simulations[i,] + runif(nSim, -0.5, 0.5))(
                        observed[i] + runif(1, -0.5, 0.5))
                } else {
                    scaledResiduals[i] <- DHARMa.ecdf(simulations[i,])(
                        observed[i])
                }
            } else {
                scaledResiduals[i] <- DHARMa.ecdf(simulations[i,])(observed[i])
            }
        }
    } else { # PIT
        for (i in 1:n) {
            minSim <- mean(simulations[i,] < observed[i])
            maxSim <- mean(simulations[i,] <= observed[i])
            if (minSim == maxSim) {
                scaledResiduals[i] <- minSim
            } else if (randomized) {
                scaledResiduals[i] <- runif(1, minSim, maxSim)
            } else {
                scaledResiduals[i] <- 0.5 * (minSim + maxSim)
            }
        }
    }
    scaledResiduals
}

#' Modified ECDF function.
#'
#' @details Ensures symmetric ECDF (standard ECDF is <), and 
#' that 0 / 1 values are only produced if the data is strictly < > 
#' than the observed data.
#'
#' @keywords internal
DHARMa.ecdf <- function (x)
{
    x <- sort(x)
    n <- length(x)
    if (n < 1)
    stop(paste("DHARMa.ecdf - length vector < 1", x))
    vals <- unique(x)
    rval <- approxfun(vals, cumsum(tabulate(match(x, vals)))/ (n +1),
    method = "linear", yleft = 0, yright = 1, ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    assign("nobs", n, envir = environment(rval))
    attr(rval, "call") <- sys.call()
    rval
}