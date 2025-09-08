##@sienaPredict make into method?
sienaPredict <- function(
    ans,
    data,
    useTieProb = TRUE,
    depvar = NULL,
    level = "period",
    condition = NULL,
    sum.fun = mean,
    na.rm = TRUE,
    uncertainty = TRUE,
    useCluster = FALSE,
    nbrNodes = 1,
    nsim = 1000,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batch_dir = "temp",
    prefix = "simBatch_b",
    combine_batch = TRUE,
    batch_size = 50,
    keep_batch = FALSE,
    verbose = TRUE
    ){
        if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
        staticContributions <- calculateContribution(ans, data, depvar)
        probName <- ifelse(useTieProb, "tieProb", "changeProb")
        probFun <- predictProbability
        predictArgs <- list(
            ans = ans,
            staticContributions = staticContributions,
            useTieProb = useTieProb
        )
        sienaPostestimate(
            predictFun = probFun,
            predictArgs = predictArgs,
            outcome = probName,
            level = level,
            condition = condition,
            sum.fun = sum.fun,
            na.rm = na.rm,
            theta_hat = ans[["theta"]],
            cov_theta = ans[["covtheta"]],
            uncertainty = uncertainty,
            nsim = nsim,
            useCluster = useCluster,
            nbrNodes = nbrNodes,
            clusterType = clusterType,
            cluster = cluster,
            batch_dir = batch_dir,
            prefix = prefix,
            combine_batch = combine_batch,
            batch_size = batch_size,
            keep_batch = keep_batch,
            verbose = verbose
        )
}

predictProbability <- function(ans, staticContributions, theta, useTieProb = TRUE) {
    effects <- ans[["effects"]] # provide effects instead?
    include <- effects[["include"]]
    includedEffects <- effects[include, ]
    noRateIncluded <- includedEffects[["type"]] != "rate"
    effectNames  <- includedEffects[["shortName"]][noRateIncluded]

    # Align theta by name
    if (!is.null(names(theta))) {
        thetaNoRate <- theta[effectNames]
    } else {
        thetaNoRate <- theta[seq_along(effectNames)]
    }


    df <- staticContributions
    df <- addUtilityColumn(df, effectNames, thetaNoRate)
    df <- addProbabilityColumn(df, group_vars = c("period", "ego"), useTieProb = useTieProb)
    df <- subset(df, df[["density"]] != 0)
    df <- conditionalReplace(df, df[["density"]] == -1, setdiff(effectNames, "density"), function(x) x * -1)
    df
}

sienaPredictDynamic <- function(
    ans,
    data,
    effects,
    algorithm,
    useTieProb = TRUE,
    depvar = NULL,
    level = "period",
    condition = NULL,
    sum.fun = mean,
    na.rm = TRUE,
    n3 = 1000,
    useChangeContributions = TRUE,
    uncertainty = TRUE,
    useCluster = FALSE,
    nbrNodes = 1,
    nsim = 100,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batch_dir = "temp",
    prefix = "simBatch_b",
    combine_batch = TRUE,
    batch_size = 10,
    keep_batch = FALSE,
    verbose = TRUE
    ){
        if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
        probName <- ifelse(useTieProb, "tieProb", "changeProb")
        predictDynamicArgs <- list(
            ans = ans,
            data = data,
            effects = effects,
            algorithm = algorithm,
            useTieProb = useTieProb,
            depvar = depvar,
            n3 = n3,
            useChangeContributions = FALSE
        )
        sienaPostestimate(
            predictFun = predictProbabilityDynamic,
            predictArgs = predictDynamicArgs,
            outcome = probName,
            level = level,
            condition = condition,
            sum.fun = sum.fun,
            na.rm = na.rm,
            theta_hat = ans[["theta"]],
            cov_theta = ans[["covtheta"]],
            uncertainty = uncertainty,
            nsim = nsim,
            useCluster = useCluster,
            nbrNodes = nbrNodes,
            clusterType = clusterType,
            cluster = cluster,
            batch_dir = batch_dir,
            prefix = prefix,
            combine_batch = combine_batch,
            batch_size = batch_size,
            keep_batch = keep_batch,
            verbose = verbose,
            useChangeContributions = useChangeContributions
        )
}

predictProbabilityDynamic <- function(ans, data, theta, algorithm, effects,
    useTieProb = TRUE, depvar, n3 = NULL, useChangeContributions = FALSE) {
    include <- effects[["include"]]
    includedEffects <- effects[include, ]
    noRateIncluded <- includedEffects[["type"]] != "rate"
    effectNames  <- includedEffects[["shortName"]][noRateIncluded]

    # Align theta by name
    if (!is.null(names(theta))) {
        thetaNoRate <- theta[effectNames]
    } else {
        thetaNoRate <- theta[seq_along(effectNames)]
    }

    df <- getChangeContributionsDynamic(
        ans = ans, 
        theta = theta, 
        data = data, 
        algorithm = algorithm,
        effects = effects, 
        depvar = depvar, 
        n3 = n3, 
        useChangeContributions = useChangeContributions, 
        returnDataFrame = TRUE
    )
    df <- widenContribution(df)
    df <- addUtilityColumn(df, effectNames, thetaNoRate)
    df <- addProbabilityColumn(df, group_vars = c("chain","period", "ministep"), useTieProb = useTieProb)
    df <- subset(df, df[["density"]] != 0)
    df <- conditionalReplace(df, df[["density"]] == -1, setdiff(effectNames, "density"), function(x) x * -1)
    df
}

sienaPostestimate <- function(
    predictFun,
    predictArgs,
    outcome,
    level = "period",
    condition = NULL,
    sum.fun = mean,
    na.rm = TRUE,
    theta_hat,
    cov_theta,
    uncertainty = TRUE,
    nsim = 1000,
    useCluster = FALSE,
    nbrNodes = 1,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batch_dir = "temp",
    prefix = "simBatch_b",
    combine_batch = TRUE,
    batch_size = 50,
    keep_batch = FALSE,
    verbose = TRUE,
    useChangeContributions = NULL
) {
    estimator <- makeEstimator(
        predictFun = predictFun,
        predictArgs = predictArgs,
        outcome = outcome,
        level = level,
        condition = condition,
        sum.fun = sum.fun,
        na.rm = na.rm
    )

    if (!is.null(useChangeContributions)) {
        expect <- estimator(theta_hat, useChangeContributions = useChangeContributions)
    } else {
        expect <- estimator(theta_hat)
    }
    if (!uncertainty) {
        return(expect)
    }
    if (!useCluster) {
        clusterType <- "FORK"
        nbrNodes <- 1
    }
    uncert <- drawSim(
        estimator = estimator,
        theta_hat = theta_hat,
        cov_theta = cov_theta,
        nbrNodes = nbrNodes,
        nsim = nsim,
        clusterType = clusterType,
        cluster = cluster,
        batch_dir = batch_dir,
        prefix = prefix,
        combine_batch = combine_batch,
        batch_size = batch_size,
        keep_batch = keep_batch,
        verbose = verbose
    )
    # is this efficient for data.table?
    uncert <- agg(outcome, uncert, level = level, condition = condition, sum.fun = summarizeValue)
    if (is.null(condition)) {
        cbind(expect, uncert)
    } else {
        merge(expect, uncert)
    }
}

makeEstimator <- function(predictFun, predictArgs, outcome,
    level = "period", condition = NULL, sum.fun = mean, na.rm = TRUE) {
    function(theta, useChangeContributions = FALSE) {
        if(!is.null(predictArgs[["useChangeContributions"]])){
          predictArgs[["useChangeContributions"]] <- useChangeContributions
        }
        callArgs <- c(list(theta = theta), predictArgs)
        unit_pred <- do.call(predictFun, callArgs)
        agg(
            ME = outcome,
            data = unit_pred,
            level = level,
            condition = condition,
            sum.fun = sum.fun,
            na.rm = na.rm
        )
    }
}

drawSim <- function(
    estimator,
    theta_hat,
    cov_theta,
    nbrNodes = 1,
    nsim = 1000,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batch_dir = "temp",
    prefix = "simBatch_b",
    batch_size = NULL,
    combine_batch = TRUE,
    keep_batch = FALSE,
    verbose = TRUE){

  if(nbrNodes > 1){
    clusterType <- match.arg(clusterType) ## maybe add explicit single core
  }
  if (!dir.exists(batch_dir)) dir.create(batch_dir, recursive = TRUE)
  
  if(is.null(batch_size)) batch_size <- floor(nsim/nbrNodes)
  ## not really set up correctly yet - true number of simulations is generally not nsim this way!
  nbatches <- ceiling(nsim / batch_size)
  
  # Cluster logic
  cluster_created <- FALSE
  cl <- cluster
  if(clusterType == "PSOCK"){
    if (is.null(cl)) {
      cl <- parallel::makeCluster(nbrNodes, type = clusterType)
      cluster_created <- TRUE
      if (verbose) message("Created new parallel cluster with ", nbrNodes, " cores.")
    } else {
      if (verbose) message("Using existing cluster.")
    }
    # Export all arguments to the cluster
    parallel::clusterExport(cl, names(sim_args), envir = list2env(sim_args, parent = .GlobalEnv))
    # Make sure that all clusters have all dependencies available 
    parallel::clusterEvalQ(cl, {
      library(data.table)
    })
  } else {
    if (.Platform$OS.type == "windows") {
      if (nbrNodes > 1){
        warning("FORK/mclapply: parallel execution not available on Windows; running serially.")
        nbrNodes <- 1
      }
    }
  }
  
  if (verbose) message("Starting batch simulations (", nbatches, " batches of up to ", batch_size, ") ...")
  
  for (b in seq_len(nbatches)) {
    first_sim <- 1 + (b - 1) * batch_size
    last_sim  <- min(b * batch_size, nsim)
    batch     <- first_sim:last_sim
    
    if (verbose) message(sprintf("[Batch %d/%d] Simulations %d-%d", b, nbatches, first_sim, last_sim))
    # Use the argument list to call estimator in the cluster
    if(clusterType == "PSOCK"){
      batch_result <- parallel::parLapply(cl, batch, function(i) {
        if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
        theta_sim <- MASS::mvrnorm(1, mu = theta_hat, Sigma = cov_theta)
        sim <- do.call(estimator, list(theta = theta_sim))
        sim[,"sim"] <- i
        sim
      })
    } else {
      batch_result <- parallel::mclapply(batch, function(i){
        if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
        theta_sim <- MASS::mvrnorm(1, mu = theta_hat, Sigma = cov_theta)
        sim <- do.call(estimator, list(theta = theta_sim))
        sim[,"sim"] <- i
        sim
      },
      mc.cores = nbrNodes
      )
    }
    file_name <- file.path(batch_dir, sprintf("%s%03d.rds", prefix, b))
    saveRDS(batch_result, file = file_name)
    rm(batch_result)
    if (verbose) message(sprintf("Saved batch file '%s' (%d%% complete)", file_name, round(100 * last_sim/nsim)))
  }
  
  if (cluster_created) {
    parallel::stopCluster(cl)
    if (verbose) message("Stopped internal cluster.")
  }
  
  if (verbose) message("All batches complete. Now combining and cleaning up...")
  
  if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads()
  
  if(combine_batch){
    # Combine all batches
    batch_files <- list.files(batch_dir, pattern = sprintf("^%s\\d{3}\\.rds$", prefix), full.names = TRUE)
    batch_files <- sort(batch_files)
    combined_results <- vector("list", nsim)
    
    for (b in seq_along(batch_files)) {
      file_name <- batch_files[b]
      results <- readRDS(file_name)
      first <- 1 + (b - 1) * batch_size
      last  <- min(b * batch_size, nsim)
      combined_results[first:last] <- results
      if(!keep_batch){
        file.remove(file_name)
        if (verbose) message(sprintf("Deleted batch %d (%d-%d)", b, first, last))
      }
    }
    
    if (verbose) message("Returning combined results.")
    if (requireNamespace("data.table", quietly = TRUE)) {
      return(data.table::rbindlist(combined_results))
    } else {
      out <- do.call(rbind, combined_results)
      rownames(out) <- NULL
      return(out)
    }
  }
}

agg <- function(ME,
                data,
                level = "none",
                condition = NULL,
                sum.fun = mean,
                na.rm = TRUE) {
  levels <- list(
    none = character(0),
    period = "period",
    ego = c("period", "ego"),
    chain = c("period", "chain"),
    ministep = c("period", "chain", "ministep")
  )
  group_vars <- c(levels[[level]], condition)

  # Helper for complex output (vector/list output from sum.fun)
  expand_summary <- function(val) {
    if (length(val) == 1 && (is.null(names(val)) || names(val) == "")) {
      setNames(list(val), ME)
    } else {
      as.list(val)
    }
  }

  # Use data.table if available
  if (requireNamespace("data.table", quietly = TRUE)) {
    if (!data.table::is.data.table(data)) {
      data <- data.table::as.data.table(data)
    }
    result <- data[, expand_summary(sum.fun(get(ME), na.rm = na.rm)), by = group_vars]
    return(result)
  }

  # ---- Fallback base R: ----
  if (length(group_vars) == 0) {
    output <- expand_summary(sum.fun(data[[ME]], na.rm = na.rm))
    output <- as.data.frame(output)
    return(output)
  }
  grouping <- interaction(data[, group_vars, drop = FALSE], drop = TRUE)
  split_list <- split(data[[ME]], grouping)
  sum_list <- lapply(split_list, sum.fun, na.rm = na.rm)
  # Turn to data.frame with group codes and result; handle multi-valued outputs
  group_levels <- unique(grouping)
  group_df <- do.call(rbind, lapply(strsplit(as.character(group_levels), ".", fixed = TRUE), as.data.frame.list))
  names(group_df) <- group_vars
  res_df <- as.data.frame(do.call(rbind, lapply(sum_list, expand_summary)))
  out <- cbind(group_df, res_df)
  rownames(out) <- NULL
  return(out)
}

summarizeValue <- function(x, na.rm = TRUE){
  if(na.rm){
    x <- x[! is.na(x)]
  }
  qu <- unname(quantile(x, probs = c(.025,0.5,.975)))
  mn <- mean(x)
  se <- sd(x)
  n <- length(x)
  as.list(c("Mean" = mn, "SE" = se, "Median" = qu[2], "q_025" = qu[1], "q_975" = qu[3], "cases" = n))
}

## Simple helper to calculate the utility given contribution and theta values
calculateUtility <- function(mat, theta) {
  stopifnot(is.matrix(mat))
  as.numeric(mat %*% theta)
}

addUtilityColumn <- function(df, effectNames, theta) {
  if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(df)) {
    df[, ("changeUtil") := calculateUtility(as.matrix(.SD), theta), .SDcols = effectNames]
  } else {
    df[["changeUtil"]] <- calculateUtility(as.matrix(df[, effectNames, drop = FALSE]), theta)
    return(df)
  }
}

addProbabilityColumn <- function(
    df,
    group_vars,
    useTieProb = FALSE
) {
  stopifnot(is.data.frame(df))
  if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(df)) {
    df[, ("changeProb") := softmax_arma(changeUtil), by = group_vars]
    if (useTieProb) {
      df[, tieProb := get("changeProb")]
      if ("density" %in% names(df)) {
        df[density == -1, tieProb := 1 - tieProb]
      }
    }
    df
  } else {
    grouping <- interaction(df[, group_vars, drop = FALSE], drop = TRUE)
    df[["changeProb"]] <- ave(df[["changeUtil"]], grouping, FUN = softmax_arma)
    if (useTieProb) {
      df[,"tieProb"] <- df[, "changeProb"]
      if ("density" %in% names(df)) {
        idx <- which(df[["density"]] == -1)
        df[idx, "tieProb"] <- 1 - df[idx, "tieProb"]
      }
    }
    df
  }
}

conditionalReplace <- function(df, row_ids, cols, fun) {
  if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(df)) {
    df[eval(substitute(row_ids)), (cols) := lapply(.SD, fun), .SDcols = cols]
  } else {
    df[row_ids, cols] <- fun(df[row_ids, cols])
  }
  df
}
