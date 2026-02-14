##@predict.sienaFit
predict.sienaFit <- function(
    object,
    newdata,
    type = c("changeProb", "tieProb"), #behavior not implemented *here*
    newParams = NULL,
    effects = NULL,
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
    batch_size = NULL,
    keep_batch = FALSE,
    verbose = TRUE,
    memory_scale = NULL,
    ...
  ){
      type <- match.arg(type)
      if (is.null(depvar)) depvar <- names(newdata[["depvars"]])[1]
      if (is.null(memory_scale)) {
        memory_scale <- compute_memory_scale(newdata, depvar, dynamic = FALSE)
      }
      staticContributions <- getStaticChangeContributions(
          ans = object,
          data = newdata,
          effects = effects,
          depvar = depvar,
          returnDataFrame = TRUE
      )
      probFun <- predictProbability
      predictArgs <- list(
          ans = object,
          staticContributions = staticContributions,
          type = type
      )
      sienaPostestimate(
          predictFun = probFun,
          predictArgs = predictArgs,
          outcome = type,
          level = level,
          condition = condition,
          sum.fun = sum.fun,
          na.rm = na.rm,
          theta_hat = object[["theta"]], # change into coef
          cov_theta = object[["covtheta"]], # change into cov
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
          memory_scale = memory_scale
      )
}

predictProbability <- function(ans, staticContributions, theta, type = "changeProb") {
    effects <- ans[["requestedEffects"]] # provide effects instead?
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
    df <- widenStaticContribution(staticContributions)
    df <- addUtilityColumn(df, effectNames, thetaNoRate)
    df <- addProbabilityColumn(df, group_vars = c("period", "ego"), type = type)
    df <- conditionalReplace(df, df[["density"]] == -1, setdiff(effectNames, "density"), function(x) x * -1)
    df
}

predictDynamic <- function(
    ans,
    newdata,
    effects,
    algorithm,
    type = c("changeProb", "tieProb"),
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
    batch_size = NULL,
    keep_batch = FALSE,
    verbose = TRUE,
    silent = TRUE,
    memory_scale = NULL
){
    type <- match.arg(type)
    if (is.null(depvar)) depvar <- names(newdata[["depvars"]])[1]

    if (is.null(memory_scale)) {
        memory_scale <- compute_memory_scale(
            data = newdata,
            depvar = depvar,
            dynamic = TRUE,
            n3 = n3
        )
    }

    predictFun <- predictProbabilityDynamic
    predictArgs <- list(
        ans = ans,
        data = newdata,
        algorithm = algorithm,
        effects = effects,
        type = type,
        depvar = depvar,
        n3 = n3
    )

    sienaPostestimate(
        predictFun = predictFun,
        predictArgs = predictArgs,
        outcome = c("changeProb", "tieProb"),
        level = level,
        condition = condition,
        sum.fun = sum.fun,
        na.rm = na.rm,
        theta_hat = ans$theta,
        cov_theta = ans$covtheta,
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
        useChangeContributions = useChangeContributions,
        memory_scale = memory_scale
    )
}

predictProbabilityDynamic <- function(ans, data, theta, algorithm, effects,
    type = "changeProb", depvar, n3 = NULL, useChangeContributions = FALSE, silent = TRUE) {
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

    df <- getDynamicChangeContributions(
        ans = ans, 
        theta = theta, 
        data = data, 
        algorithm = algorithm,
        effects = effects, 
        depvar = depvar, 
        n3 = n3, 
        useChangeContributions = useChangeContributions, 
        returnDataFrame = TRUE,
        silent = silent
    )
    df <- widenDynamicContribution(df)
    df <- addUtilityColumn(df, effectNames, thetaNoRate)
    df <- addProbabilityColumn(df, group_vars = c("chain","period", "ministep"), type = type)
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
    batch_size = NULL,
    keep_batch = FALSE,
    verbose = TRUE,
    useChangeContributions = NULL,
    gc_each_batch = TRUE,
    gc_each_sim = FALSE,
    memory_scale = NULL
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

    if (!uncertainty) return(expect)

    if (is.null(memory_scale)) memory_scale <- 1L

    uncert <- drawSim(
      estimator = estimator,
        theta_hat = theta_hat,
        cov_theta = cov_theta,
        nbrNodes = if (useCluster) nbrNodes else 1L,
        nsim = nsim,
        clusterType = clusterType,
        cluster = cluster,
        batch_dir = batch_dir,
        prefix = prefix,
        combine_batch = combine_batch,
        batch_size = batch_size,
        keep_batch = keep_batch,
        verbose = verbose,
        gc_each_batch = gc_each_batch,
        gc_each_sim = gc_each_sim,
        memory_scale = memory_scale
    )

    uncert <- agg(outcome, uncert, level = level, condition = condition, sum.fun = summarizeValue)
    mergeEstimates(expect, uncert, level = level, condition = condition)
}

makeEstimator <- function(predictFun, predictArgs, outcome,
    level = "period", condition = NULL, sum.fun = mean, na.rm = TRUE) {
    function(theta, useChangeContributions = FALSE) {
        if(!is.null(predictArgs[["useChangeContributions"]])){
          predictArgs[["useChangeContributions"]] <- useChangeContributions
        }
        callArgs <- c(list(theta = theta), predictArgs)
        unit_pred <- do.call(predictFun, callArgs)
        main_result <- agg(
            ME = outcome,
            data = unit_pred,
            level = level,
            condition = condition,
            sum.fun = sum.fun,
            na.rm = na.rm
        )
        if (isTRUE(predictArgs$details)) {
          details <- unit_pred
          return(list(summary = main_result, details = details))
        } else {
          return(main_result)
        }
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
    verbose = TRUE,
    gc_each_batch = TRUE,
    gc_each_sim = FALSE,
    memory_scale = 1
) {
    if (nbrNodes > 1) {
        clusterType <- match.arg(clusterType)
    } else {
        clusterType <- "FORK"
    }

    if (!dir.exists(batch_dir)) dir.create(batch_dir, recursive = TRUE)

    batch_size <- resolve_batch_size(
        nsim = nsim,
        nbrNodes = nbrNodes,
        batch_size = batch_size,
        memory_scale = memory_scale
    )
    nbatches <- ceiling(nsim / batch_size)
  
    # Cluster logic
    cluster_created <- FALSE
    cl <- cluster
    use_psock <- (clusterType == "PSOCK" && nbrNodes > 1)
    use_fork  <- (clusterType == "FORK"  && nbrNodes > 1 && .Platform$OS.type != "windows")
    serial_mode <- !(use_psock || use_fork)

    if (use_psock) {
      if (is.null(cl)) {
        cl <- parallel::makeCluster(nbrNodes, type = clusterType)
        cluster_created <- TRUE
        if (verbose) {
          message("Created new parallel cluster with ", nbrNodes, " cores.")
        }
      }
      parallel::clusterExport(cl, 
        c("estimator", "theta_hat", "cov_theta"), 
        envir = environment()
      )
      parallel::clusterEvalQ(cl, { 
          if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1) 
      })
    }
    
    if (verbose) message("Starting batch simulations (", nbatches, 
      " batches of up to ", batch_size, ") ...")
    
    for (b in seq_len(nbatches)) {
      first_sim <- 1 + (b - 1) * batch_size
      last_sim  <- min(b * batch_size, nsim)
      batch     <- first_sim:last_sim

      if (use_psock) {
        batch_result <- parallel::parLapply(cl, batch, function(i) {
          theta_sim <- MASS::mvrnorm(1, mu = theta_hat, Sigma = cov_theta)
          sim <- do.call(estimator, list(theta = theta_sim))
          sim[, "sim"] <- i
          sim
        })
      } else if (use_fork) {
        batch_result <- parallel::mclapply(batch, function(i) {
          if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
          theta_sim <- MASS::mvrnorm(1, mu = theta_hat, Sigma = cov_theta)
          sim <- do.call(estimator, list(theta = theta_sim))
          sim[, "sim"] <- i
          sim
        }, mc.cores = nbrNodes)
      } else {
        batch_result <- vector("list", length(batch))
        for (k in seq_along(batch)) {
          i <- batch[k]
          theta_sim <- MASS::mvrnorm(1, mu = theta_hat, Sigma = cov_theta)
          sim <- do.call(estimator, list(theta = theta_sim))
          sim[, "sim"] <- i
          batch_result[[k]] <- sim
          if (gc_each_sim) gc(verbose = FALSE)
        }
      }

      file_name <- file.path(batch_dir, sprintf("%s%03d.rds", prefix, b))
      saveRDS(batch_result, file = file_name)
      rm(batch_result)
      if (gc_each_batch) gc(verbose = FALSE)
    }
    
    if (cluster_created) {
      parallel::stopCluster(cl)
      if (verbose) message("Stopped internal cluster.")
    }
    
    if (verbose) message("All batches complete. Now combining and cleaning up...")
    
    if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads()
    if(combine_batch){
      # Combine all batches
      batch_files <- list.files(batch_dir, 
        pattern = sprintf("^%s\\d{3}\\.rds$", prefix), full.names = TRUE)
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

compute_memory_scale <- function(data, depvar, dynamic = FALSE, n3 = NULL,
                                 cell_ref = 5000L, n3_ref = 500L) {
    dv <- data$depvars[[depvar]]

    n_actor <- if (!is.null(dim(dv))) {
        as.integer(nrow(dv))
    } else {
        as.integer(length(dv))
    }

    n_choice <- if (!is.null(dim(dv)) && length(dim(dv)) >= 2) {
        as.integer(ncol(dv))
    } else {
        3L
    }

    size_scale <- max(1L, ceiling((n_actor * n_choice) / as.integer(cell_ref)))

    if (isTRUE(dynamic)) {
        n3_val <- if (is.null(n3)) 1L else as.integer(n3)
        n3_scale <- max(1L, ceiling(n3_val / as.integer(n3_ref)))
    } else {
        n3_scale <- 1L
    }

    as.integer(max(1L, size_scale * n3_scale))
}

resolve_batch_size <- function(nsim, nbrNodes = 1L, batch_size = NULL, memory_scale = 1L) {
    if (!is.null(batch_size)) {
        return(min(max(1L, as.integer(batch_size)), as.integer(nsim)))
    }

    base_batch <- floor(as.integer(nsim) / max(1L, as.integer(nbrNodes)))
    auto_batch <- max(1L, floor(base_batch / max(1L, as.integer(memory_scale))))
    min(auto_batch, as.integer(nsim))
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
    egoChoice = c("period", "ego", "choice"),
    chain = c("period", "chain"),
    ministep = c("period", "chain", "ministep"),
    ministepChoice = c("period", "chain", "ministep", "choice")
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
    result <- data[, expand_summary(sum.fun(get(ME), na.rm = na.rm)), 
      by = group_vars]
    return(result)
  }
  # ---- Fallback base R: ----
  if (length(group_vars) == 0) {
    output <- expand_summary(sum.fun(data[[ME]], na.rm = na.rm))
    output <- as.data.frame(output)
    return(output)
  }

  grouping <- interaction(data[, group_vars, drop = FALSE], drop = TRUE)
  split_data <- split(data, grouping)
  agg_list <- lapply(split_data, function(subdf) {
    vals <- subdf[[ME]]
    res <- expand_summary(sum.fun(vals, na.rm = na.rm))
    # Add group columns
    group_vals <- subdf[1, group_vars, drop = FALSE]
    cbind(group_vals, as.data.frame(res, stringsAsFactors = FALSE))
  })
  out <- do.call(rbind, agg_list)
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
  available_effects <- effectNames[effectNames %in% names(df)]
  if (!length(available_effects)) {
    stop("No requested effect columns found in data for utility calculation")
  }

  if (!is.null(names(theta))) {
    theta_use <- theta[available_effects]
  } else {
    theta_use <- theta[seq_along(available_effects)]
  }

  if (requireNamespace("data.table", quietly = TRUE) && 
      data.table::is.data.table(df)) {
    df[, ("changeUtil") := calculateUtility(
      as.matrix(.SD), 
      theta_use), 
      .SDcols = available_effects]
  } else {
    df[["changeUtil"]] <- calculateUtility(
      as.matrix(df[, available_effects, drop = FALSE]), 
      theta_use)
    return(df)
  }
}

addProbabilityColumn <- function(
    df,
    group_vars,
    type = "changeProb"
) {
  stopifnot(is.data.frame(df))
  if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(df)) {
    # For non standard evaluation
    changeUtil <- tieProb <- NULL
    df[, ("changeProb") := softmax_arma(changeUtil), by = group_vars]
    if (type == "tieProb") {
      df[, tieProb := get("changeProb")]
      if ("density" %in% names(df)) {
        df[density == -1, tieProb := 1 - tieProb]
        df[density == 0, tieProb := NA]
      }
    }
    df
  } else {
    grouping <- interaction(df[, group_vars, drop = FALSE], drop = TRUE)
    df[["changeProb"]] <- ave(df[["changeUtil"]], 
      grouping, 
      FUN = softmax_arma)
    if (is.list(df[["changeProb"]])) df[["changeProb"]] <- unlist(df[["changeProb"]])
    if (type == "tieProb") {
      df[,"tieProb"] <- df[, "changeProb"]
      if ("density" %in% names(df)) {
        idx <- which(df[["density"]] == -1)
        df[idx, "tieProb"] <- 1 - df[idx, "tieProb"]
        df[df[["density"]] == 0, "tieProb"] <- NA
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

mergeEstimates <- function(df1, df2, level = "none", condition = NULL) {
  levels_list <- list(
    none = character(0),
    period = "period",
    ego = c("period", "ego"),
    egoChoice = c("period", "ego", "choice"),
    chain = c("period", "chain"),
    ministep = c("period", "chain", "ministep"),
    ministepChoice = c("period", "chain", "ministep", "choice")
  )
  group_vars <- c(levels_list[[level]], condition)
  if (length(group_vars) > 0) {
    merge(df1, df2, by = group_vars, all = TRUE)
  } else {
    cbind(df1, df2)
  }
}