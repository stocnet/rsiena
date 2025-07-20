#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena/
# *
# * File: sienaMargins.r
# *
# * Description: Calculates predicted edge probabilities and
# * (average) marginal effects
# *****************************************************************************/

# Currently uses getChangeContribution from sienaRI.R
# We might try to use RSiena:::getTargets() with returnStaticCHangeContribution = TRUE?
# See also RSiena:::getTheActorStatistics()

# The following calculations are only correct if the effect statistics are independent of each other
# Extension for explicit interactions are relatively simple but are missing so far

##@sienaAME Modified from ergmMargins::ergmAME. Use as RSiena:sienaAME.
# Current implementation is only sensible for alternative-varying effects,
# i.e. no covariate related acivitiy effects.
# Effects statistics are also assumed to be independent of each other
sienaAME<-function(expectedChangeProbabilities,
                   modelFit,
                   shortNames,
                   simSE = FALSE,
                   simN = 100){
  probs <- expectedChangeProbabilities
  theta <- modelFit$theta
  effects <- modelFit$effects
  effectnumbers <- which(effects$shortName %in% shortNames)
  effectnames <- effects$effectName[effectnumbers]
  vc <- modelFit$covtheta

  ##direct marginal effects with no interaction; 
  AME.fun <- function(theta) { ## also depends on P?
    ME.ergm <- directMarginalEffect(theta[effectnumbers], probs)
    colMeans(ME.ergm, na.rm = TRUE)
    # A probability weighted average might be more sensible
    # colMeans(probs * ME.ergm, na.rm = TRUE)
  }
# if(simSE == FALSE){
  ## I am not sure if this captures the uncertainty properly
  AME <- AME.fun(theta)
  Jac <- numDeriv::jacobian(AME.fun, theta)
  variance.ame <- Jac %*% vc %*% t(Jac)

  AME.se <- sqrt(diag(variance.ame))
  AME.z <- AME / AME.se
  P.AME <- 2 * (stats::pnorm(-abs(AME.z))) # should probably be t_{1-alpha, n} instead of 2

  result <- cbind(AME, AME.se, AME.z, P.AME)
  colnames(result) <- c("AME", "Delta SE", "Z", "probs")
  rownames(result) <- effectnames
  result <- signif(result, digits = 5) #should be made more flexible
  # }else{
  # ## Draw coefficients from multivariate normal
  # simCoef <- MASS::mvrnorm(n = simN, theta, vc)
  # ## loop over each vector of coefficients ... 
  # sapply(1:lengthof(simCoef), 
  # }
  result
}

##@directMarginalEffect Effect of a change in alternative-specific statistic on its own probability. Use as RSiena:::directMarginalEffect
directMarginalEffect <- function(thetas, expectedChangeProbabilities, interactions = NULL) {
  probs <- expectedChangeProbabilities
  sapply(thetas, function(theta) {
    if (!is.null(interactions)) {
      d_util <- theta1 + theta2 * cont
    }else{
      d_util <- theta
    }
    probs * (1 - probs) * theta
  })
}

# crossMarginalEffect Effect of a change in alternative-specific statistic on all other alternative's probabilities is not implemented yet
# The formula is P_ij * -P_ih * beta_ij and leads to n-1 cross-marginal effects for each potential choice i -> j



## this should be divided in "getChangeStatistics" (maybe include in contributions?) and "calculateChangeProbabilities"
##@expectedChangeProbabilities. Use as RSiena:::expectedChangeProbabilities
expectedChangeProbabilities <- function(conts, effects, theta, thedata = NULL,
                                        getChangeStatistics = FALSE, effectNames = NULL) {
  waves <- length(conts[[1]])
  effects <- effects[effects$include == TRUE, ]
  noRate <- effects$type != "rate"
  effects <- effects[noRate, ]
  if (sum(noRate) != length(theta)) {
    theta <- theta[noRate]
  }
  effectNa <- attr(conts, "effectNames")
  effectTypes <- attr(conts, "effectTypes")
  networkNames <- attr(conts, "networkNames")
  networkTypes <- attr(conts, "networkTypes")
  networkInteraction <- effects$interaction1
  effectIds <- paste(effectNa, effectTypes, networkInteraction, sep = ".")
  currentDepName <- ""
  depNumber <- 0
  for(eff in 1:length(effectIds)) { # seq_along throws an error?
    if(networkNames[eff] != currentDepName) {
      currentDepName <- networkNames[eff]
      actors <- length(conts[[1]][[1]][[1]])
      depNumber <- depNumber + 1
      currentDepEffs <- effects$name == currentDepName
      depNetwork <- thedata$depvars[[depNumber]]
      if (networkTypes[eff] == "oneMode") {
        choices <- actors
      } else if (networkTypes[eff] == "behavior") {
        choices <- 3
      } else if (networkTypes[eff] == "bipartite") {
        if (dim(depNetwork)[2] >= actors) {
          stop("does not work for bipartite networks with second mode >= first mode")
        }
        choices <- dim(depNetwork)[2] + 1
      } else {
        stop("does not work for dependent variables of type 'continuous'")
      }
      # impute for wave 1
      if (networkTypes[eff] %in% c("oneMode", "bipartite")) {
        depNetwork[, , 1][is.na(depNetwork[, , 1])] <- 0
      } else {
        depNetwork[, , 1][is.na(depNetwork[, , 1])] <- attr(depNetwork, "modes")[1]
      }
      # impute for next waves;
      # this may be undesirable for structurals immediately followed by NA...
      for (m in 2:dim(depNetwork)[3]) {
        depNetwork[, , m][is.na(depNetwork[, , m])] <- depNetwork[, , m - 1][is.na(depNetwork[, , m])]
      }
      # Make sure the diagonals are not treated as structurals
      if (networkTypes[eff] == "oneMode") {
        for (m in 1:(dim(depNetwork)[3])) {
          diag(depNetwork[, , m]) <- 0
        }
      }
      structurals <- (depNetwork >= 10)
      if (networkTypes[eff] == "oneMode") {
        if (attr(depNetwork, "symmetric")) {
          message("\nNote that for symmetric networks, effect sizes are for modelType 2 (forcing).")
        }
      }
      #			currentDepObjEffsNames <- paste(effects$shortName[currentDepEffs],
      #				effects$type[currentDepEffs],effects$interaction1[currentDepEffs],sep=".")
      #			otherObjEffsNames <- paste(effects$shortName[!currentDepEffs],
      #				effects$type[!currentDepEffs],effects$interaction1[!currentDepEffs],sep=".")
      changeStats <- list()
      if (networkTypes[eff] == "behavior") {
        toggleProbabilities <- array(0, dim = c(actors, 3, waves))
      } else {
        toggleProbabilities <- array(0, dim = c(actors, choices, waves))
      }
      for(w in 1:waves) {
        currentDepEffectContributions <- conts[[1]][[w]][currentDepEffs]
        if (networkTypes[eff] == "bipartite") {
          currentDepEffectContributions <- lapply(
            currentDepEffectContributions,
            function(x) lapply(x, function(xx) xx[1:choices])
          )
        }
        # conts[[1]] is periods by effects by actors by actors
        currentDepEffectContributions <-
          sapply(lapply(currentDepEffectContributions, unlist),
                 matrix, nrow = actors, ncol = choices, byrow = TRUE,
                 simplify = "array")
        cdec <- apply(currentDepEffectContributions, c(2, 1), as.matrix)
        # cdec is effects by actors (alters) by actors (egos)
        if (dim(currentDepEffectContributions)[3] <= 1) { # only one effect
          cdec <- array(cdec, dim = c(1, dim(cdec)))
        }
        rownames(cdec) <- effectNa[currentDepEffs]
        if (getChangeStatistics) {
          changeStats[[w]] <- cdec
        }
        # replace structural 0s and 1s by NA,
        # so they are omitted from calculation
        if (networkTypes[eff] == "oneMode") {
          #	structuralsw <- structurals[,,w]
          for (ff in 1:(dim(cdec)[1])) {
            cdec[ff, , ][t(structurals[, , w])] <- NA
          }
        }
        distributions <- apply(cdec, 3,
                               calculateChoiceProbability, theta[which(currentDepEffs)])
        # matrix manipulation should be more efficient solution
        # distributions is an array of actor by choices
        if (networkTypes[eff] == "behavior") { #makes no difference here?
          toggleProbabilities[, , w] <- t(distributions)
        } else {
          toggleProbabilities[, , w] <- t(distributions)
        }
      }
    }
  }
  # toggleProbabilities is an array of choices by actors by wave
  attr(toggleProbabilities, "version") <- packageDescription(pkgname, fields = "Version")
  toggleProbabilities
}

# The following is inefficient but does extract probabilities for chains

## TODo: factor out changeContribution calculation over chains (might be more efficient in one, but less understandable)

##@expectedChangeDynamics. Use as RSiena:::expectedChangeDynamics.
# Simulates sequences of micro-steps and calculates the predicted probabilities each time
expectedChangeDynamics <- function(data = NULL, theta = NULL, algorithm = NULL, effects = NULL, depvar = NULL,
                                   returnActorStatistics = NULL, returnChangeContributions = FALSE) {
  x <- algorithm
  currentNetName <- depvar
  z  <-  NULL
  z$FRAN <- getFromNamespace(x$FRANname, pkgname)
  x$cconditional <-  FALSE
  z$print <- FALSE
  z$Phase <- 3
  z <- initializeFRAN(z, x, data, effects, prevAns=NULL, initC=FALSE, returnDeps=FALSE)
  z$returnChangeContributions <- TRUE
  z$theta <- theta
  if (!is.null(x$randomSeed))
  {
    set.seed(x$randomSeed, kind="default")
  } else {
    if (exists(".Random.seed")) {
      rm(.Random.seed, pos = 1)
      RNGkind(kind = "default")
    }
  }
  chains <- x$n3
  periods <- data$observation-1
  effects <- effects[effects$include == TRUE,]
  noRate <- effects$type != "rate"
  thetaNoRate <- theta[noRate]
  #	networkName <- effects$name[noRate]
  currentNetObjEffs <- effects$name[noRate] == currentNetName
  output <- list()
  for (chain in (1:chains))
  {
    # cat("The following line leads to an error\n")
    # browser()
    output[[chain]] <- list() #probably superfluous
    ans <- z$FRAN(z, x)
    for(period in 1:periods) {
      output[[chain]][[period]] <- list()
      microSteps <- length(ans$changeContributions[[1]][[period]])
      for(microStep in 1:microSteps) {
        if(attr(ans$changeContributions[[1]][[period]][[microStep]],
                "networkName")==currentNetName) {
          cdec <- ans$changeContributions[[1]][[period]][[microStep]]
          distributions <- calculateChoiceProbability(
            cdec,
            thetaNoRate[currentNetObjEffs]
          )
          if (returnChangeContributions) {
            output[[chain]][[period]][[microStep]] <- cbind(t(distributions), t(cdec))
            ## add colnames
          } else {
            output[[chain]][[period]][[microStep]] <- t(distributions)
            ## add colname
          }
        }
      }
    }
  }
  output
}

##@calculateChoiceProbability. Use as RSiena:::calculateChoiceProbability (just a simplified calculateDistribution)
# Calculate the probability of each potential choice for the focal actor
calculateChoiceProbability <- function(effectContributions = NULL, theta = NULL, diff = NULL) {
  nchoices <- dim(effectContributions)[2]
  distributions <- array(NA, dim = c(1, nchoices)) # probably could also work with a simple vector
  the.choices <- !is.na(colSums(effectContributions))
  if (sum(the.choices) >= 2) { # should produce an error instead of an empty array?
    utility <- colSums(theta * effectContributions[, the.choices, drop = FALSE], na.rm = TRUE)
    # why not theta %*% effectContributions[,the.choices,drop = FALSE] ?
    distributions[1, the.choices] <- softmax(utility)
  }
  distributions # returns the choice probabilitiy for the focal actor
}

##@softmax Recursive softmax formula, see https://rpubs.com/FJRubio/softmax. Use as RSiena:::softmax
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
