##@sienaPredict make into method?
sienaPredict <- function(ans,
                         data,
                         effectNames,
                         tieProb = TRUE,
                         depvar = NULL,
                         level = "period",
                         condition = NULL,
                         sum.fun = mean,
                         na.rm = TRUE,
                         keep = NULL,
                         useCluster = FALSE,
                         nbrNodes = 1,
                         uncertainty = TRUE,
                         nsim = 1000,
                         clusterType = c("PSOCK", "FORK"),
                         cluster = NULL,
                         batch_dir = "temp",
                         prefix = "simBatch_b",
                         combine_batch = TRUE,
                         keep_batch = FALSE,
                         verbose = TRUE){
  sim_fun <- simProb
  sim_args <- list(
    ans = ans,
    data = data,
    effectNames = effectNames,
    tieProb = tieProb,
    depvar = depvar,
    sim_theta = TRUE,
    level = level,
    condition = condition,
    sum.fun = sum.fun,
    na.rm = na.rm
  )
  ME <- ifelse(tieProb, "tieProb", "changeProb")
  
  if(!useCluster){
    # not nice, but works since mclapply reduces to lapply for nbrNodes = 1
    clusterType <- "FORK" 
    nbrNodes <- 1
  }
  
  # Point estimate
  sim_args_theta <- c(sim_args, keep = keep)
  sim_args_theta$sim_theta <- FALSE
  
  point_AME <- do.call(sim_fun, sim_args_theta)
  # point_AME <- agg(ME, point_AME, level = level, condition = condition, sum.fun = sum.fun)
  # prob_df <- calculateChoiceProbability(ans, data, tieProb = tieProb)
  # prob_df <- subset(prob_df, density != 0) # set to NA?
  # prob_df[prob_df[, "density"] == -1, setdiff(effectNames, "density")] <- prob_df[prob_df[, "density"] == -1, setdiff(effectNames, "density")] * -1
  # 
  # if(tieProb){
  #   prob_table <- agg("tieProb", data = prob_df, condition = condition)
  # } else {
  #   prob_table <- agg("changeProb", data = prob_df, condition = condition)
  # }
  # prob_table
  uncert_AME <- drawSim(
    sim_fun = sim_fun,
    sim_args = sim_args,
    nbrNodes = nbrNodes,
    nsim = nsim,
    clusterType = clusterType,
    cluster = cluster,
    batch_dir = batch_dir,
    prefix = prefix,
    combine_batch = combine_batch,
    keep_batch = keep_batch,
    verbose = verbose
  )
  
  uncert_AME <- agg(ME, uncert_AME, level = level, condition = condition, sum.fun = summarizeValue)
  
  if(is.null(condition)){
    cbind(point_AME, uncert_AME) # merge would be safer, but does not handle results without id vars well
  } else {
    merge(point_AME, uncert_AME)
  }
}

simProb <- function(ans, data,
                    effectNames,
                    tieProb = TRUE, 
                    depvar = NULL, # curently unused
                    sim_theta = TRUE,
                    aggregateValues = TRUE,
                    level = "period",
                    condition = NULL,
                    sum.fun = mean,
                    na.rm = TRUE,
                    keep = NULL){
  ## effectNames can just be extracted
  if(sim_theta){
    ## might not work with estimated rate effects
    theta <- MASS::mvrnorm(n=1,
                           mu = ans$theta,
                           Sigma = ans$covtheta)
    ## add option to change algorithm -> number of chains -> 
  }else{
    theta <- ans$theta
  }
  prob_sim <- calculateChoiceProbability(ans, data, tieProb = tieProb, theta = theta)
  
  # density == 0 is only relevant for the normalizing constant, not needed for difference calculation
  prob_sim <- subset(prob_sim, density != 0)
  
  # Transform contributions to change statistics for aggregation and output
  prob_sim[prob_sim[, "density"] == -1, setdiff(effectNames, "density")] <- prob_sim[prob_sim[, "density"] == -1, setdiff(effectNames, "density")] * -1

    if(aggregateValues) {
      if(tieProb){
        prob_sim <- agg("tieProb",
                          prob_sim,
                          level = level,
                          condition = condition,
                          sum.fun = sum.fun,
                          na.rm = TRUE,
                          keep = keep)
      } else {
        prob_sim <- agg("changeProb",
                          prob_sim,
                          level = level,
                          condition = condition,
                          sum.fun = sum.fun,
                          na.rm = TRUE,
                          keep = keep)
      }

  }
  prob_sim
}

##@sienaPredictDynamic make into method? own dynamicsfile?
sienaPredictDynamic <- function(ans,
                         data,
                         effectNames,
                         effects,
                         tieProb = TRUE,
                         depvar = NULL,
                         level = "period",
                         condition = NULL,
                         sum.fun = mean,
                         na.rm = TRUE,
                         keep = NULL,
                         uncertainty = TRUE,
                         nsim = 1000,
                         useCluster = FALSE, 
                         nbrNodes = 1,
                         clusterType = c("PSOCK", "FORK"),
                         cluster = NULL,
                         batch_dir = "temp",
                         prefix = "simBatch_b",
                         batch_size = 1,
                         combine_batch = TRUE,
                         keep_batch = FALSE,
                         verbose = TRUE){
  sim_fun <- simProbDynamic
  sim_args <- list(
    ans = ans,
    data = data,
    effectNames = effectNames,
    effects = effects,
    tieProb = tieProb,
    depvar = depvar,
    sim_theta = TRUE,
    level = level,
    condition = condition,
    sum.fun = sum.fun,
    na.rm = na.rm
  )
  ME <- ifelse(tieProb, "tieProb", "changeProb")

  if(!useCluster){
    # not nice, but works since mclapply reduces to lapply for nbrNodes = 1
    clusterType <- "FORK" 
    nbrNodes <- 1
  }
  
  # Point estimate
  sim_args_theta <- c(sim_args, keep = keep)
  sim_args_theta$sim_theta <- FALSE
  
  point_AME <- do.call(sim_fun, sim_args_theta)

  if(!uncertainty){
    return(point_AME)
  }else{
    uncert_AME <- drawSim(
      sim_fun = sim_fun,
      sim_args = sim_args,
      nbrNodes = nbrNodes,
      nsim = nsim,
      clusterType = clusterType,
      cluster = cluster,
      batch_dir = batch_dir,
      prefix = prefix,
      batch_size = batch_size,
      combine_batch = combine_batch,
      keep_batch = keep_batch,
      verbose = verbose
    )
    
    uncert_AME <- agg(ME, uncert_AME, level = level, condition = condition, sum.fun = summarizeValue)
    
    if(is.null(condition)){
      return(cbind(point_AME, uncert_AME)) # merge would be safer, but does not handle results without id vars well
    } else {
      return(merge(point_AME, uncert_AME))
    }
  }
}

simProbDynamic <- function(ans, data,
                            effectNames, 
                            effects,
                            tieProb = TRUE, 
                            depvar = NULL,
                            sim_theta = TRUE,
                            aggregateValues = TRUE,
                            level = "none",
                            condition = NULL,
                            sum.fun = mean,
                            na.rm = TRUE,
                            keep = NULL){
  changeProb <- changeUtil <- chain <- period <- ministep <- NULL # To resolve R CMD checks not understanding data.table syntax
  if(is.null(depvar)){
    depvar <- names(data$depvars)[1]
  }
  n_choices <- dim(data[["depvars"]][[depvar]])[2]
  
  if(sim_theta){
    ## might not work with estimated rate effects
    theta <- mvrnorm(n=1,
                           mu = ans$theta,
                           Sigma = ans$covtheta)
    ## add option to change algorithm -> number of chains -> 
  }else{
    theta <- ans$theta
  }
  
  cont_df <- calculateContributionsDynamic(
    data = data,  
    theta = c(ans$rate, theta),
    algorithm = ans$x,
    effects = effects,
    depvar = depvar
  )
  
  prob_sim <- as.data.frame(cont_df)
  
  prob_sim[,"changeUtil"] <- calculateUtility(prob_sim[,effectNames], theta)
  
  prob_sim <- as.data.table(prob_sim)
  prob_sim[, changeProb := {
    ex <- exp(changeUtil - max(changeUtil))
    ex / sum(ex)
  }, by = .(chain, period, ministep)] # does work but not accepted by R CMD check! TODO
  prob_sim <- as.data.frame(prob_sim)
  
  # prob_sim[,"changeProb"] <- with(prob_sim, ave(changeUtil, chain, period, ministep, FUN = softmax))
  # grp <- with(prob_sim, interaction(chain, period, ministep, drop=TRUE))
  # 
  # # Compute grouped softmax using lapply
  # prob_sim[,"changeProb"] <- with(prob_sim, unsplit(
  #   lapply(split(changeUtil, grp), softmax),
  #   grp
  #   )
  # )
  
  # density == 0 is only relevant for the normalizing constant, not needed for difference calculation
  prob_sim <- subset(prob_sim, density != 0)
  
  if(tieProb == TRUE){
    prob_sim[,"tieProb"] <- prob_sim[,"changeProb"]
    prob_sim[prob_sim[,"density"] == -1,"tieProb"] <- 1 - prob_sim[prob_sim[,"density"] == -1,"changeProb"]
  }
  
  # Transform contributions to change statistics for aggregation and output
  prob_sim[prob_sim[, "density"] == -1, setdiff(effectNames, "density")] <- prob_sim[prob_sim[, "density"] == -1, setdiff(effectNames, "density")] * -1
  
  if(aggregateValues) {
    if(tieProb){
      prob_sim <- agg("tieProb",
                        prob_sim,
                        level = level,
                        condition = condition,
                        sum.fun = sum.fun,
                        na.rm = TRUE,
                        keep = keep)
    } else {
      prob_sim <- agg("changeProb",
                        prob_sim,
                        level = level,
                        condition = condition,
                        sum.fun = sum.fun,
                        na.rm = TRUE,
                        keep = keep)
    }
    
  }
  prob_sim
}

## Simple helper to calculate the utility given contribution and theta values
calculateUtility <- function(conts, theta){
  as.matrix(conts) %*% theta
}

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
  }
  val
}

calculateChoiceProbability <- function(ans, data,
                                       depvar = NULL,
                                       algorithm = NULL, 
                                       # diff1 = NULL, effectName1 = NULL, range1 = NULL,
                                       # diff2 = NULL, effectName2 = NULL, range2 = NULL,
                                       tieProb = FALSE,
                                       theta = NULL) {
  ## theta should be possible to change?
  ## insert option to draw from coefficients and/or use delta method here? -> probably better outside function
  

  effects <- ans$effects
  noRate <- effects$type != "rate"
  effects <- effects[noRate, ]
  effectNames <- effects[effects$include == TRUE,"shortName"]
  
  if(is.null(theta)){
    theta <- ans$theta
  }

  if(length(theta) > length(effectNames)){
    theta <- theta[noRate]
  }
  
  df <- calculateContribution(ans, data, depvar, algorithm)
  
  df[,"changeUtil"] <- calculateUtility(df[,effectNames], theta)
  df[,"changeProb"] <- with(df, ave(changeUtil, period, ego, FUN = softmax))
  
  # ## Currently all utilites and probabilities are recalculated
  # ## Could only be done for those that change & use simple + theta * diff_cont
  
  if(tieProb == TRUE){
    df[,"tieProb"] <- df[,"changeProb"]
    df[df[,"density"] == -1,"tieProb"] <- 1 - df[df[,"density"] == -1,"changeProb"]
  }

  return(df)
}

## Wrapper function to run nsim chains drawing from the multivariate normal
## of the theta coefficients in parallel depending on setup
## Base parallel setup has to be done by the user
## export?
drawSim <- function(
    sim_fun,           # simulation function (e.g. simFirstDiffDynamic or simSecondDiffDynamic)
    sim_args,          # named list of arguments for sim_fun
    nbrNodes = 1,
    nsim = 1000,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batch_dir = "temp",
    prefix = "simBatch_b",
    batch_size = NULL,
    combine_batch = TRUE,
    keep_batch = FALSE,
    verbose = TRUE
) {
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
    # Use the argument list to call sim_fun in the cluster
    if(clusterType == "PSOCK"){
      batch_result <- parallel::parLapply(cl, batch, function(i) {
        setDTthreads(1)
        sim <- do.call(sim_fun, sim_args)
        sim[,"sim"] <- i
        sim
      })
    } else {
      batch_result <- parallel::mclapply(batch, function(i){
        if (nbrNodes != 1) setDTthreads(1)
        sim <- do.call(sim_fun, sim_args)
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
  
  setDTthreads()
  
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
    data.table::rbindlist(combined_results)
  }
}

agg <- function(ME,
                  data,
                  level = "none",
                  condition = NULL,
                  sum.fun = mean,
                  na.rm = TRUE,
                  keep = NULL){
  levels <- list(
    none = character(0),
    period = "period",
    ego = c("period", "ego"),
    chain = c("period", "chain"),
    ministep = c("period", "chain", "ministep")
  )
  group_vars <- c(levels[[level]], condition, keep)
  data <- as.data.table(data)

  ## extract?
  expand_summary <- function(val){
    if (length(val) == 1 && (is.null(names(val)) || names(val) == "")){
      setNames(list(val), ME)
    }else{
      as.list(val)
    }
  }
  data[, expand_summary(sum.fun(get(ME), na.rm = na.rm)), by = group_vars]
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

