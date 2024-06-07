#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena/
# *
# * File: sienaMargins.r
# *
# * Description: Used to calculate predicted edge probabilities, to be extended to 
# * calculate (average) marginal effects
# *****************************************************************************/

## Currently uses getChangeContribution from RIsiena.R
## We might try to use RSiena:::getTargets() with returnStaticCHangeContribution = TRUE?



##@softmax Recursive softmax formula, see https://rpubs.com/FJRubio/softmax. Use as RSiena:::softmax
softmax <- function(par, recursive = FALSE){
  if (recursive == TRUE) {
    n.par <- length(par)
    par1 <- sort(par, decreasing = TRUE)
    Lk <- par1[1]
    for (k in 1:(n.par-1)) {
      Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk)))
    }
    val <- exp(par - Lk)
    val
  }
  else {
  val <- exp(par) /
    sum(exp(par))
  }
}

##@calculateChoiceProbability. Use as RSiena:::calculateChoiceProbability (just a simplified calculateDistribution)
calculateChoiceProbability <- function(effectContributions = NULL, theta = NULL)
{
  nchoices <- dim(effectContributions)[2]
  distributions <- array(NA, dim = c(1,nchoices)) #probably could also work with a simple vector
  the.choices <- !is.na(colSums(effectContributions))
  if (sum(the.choices) >= 2) { ## should produce an error instead of an empty array?
    utilitiy <- colSums(theta*effectContributions[,the.choices,drop = FALSE], na.rm = TRUE) # why not theta %*% effectContributions[,the.choices,drop = FALSE] ?
    distributions[1,the.choices] <- softmax(utilitiy)
  }
  distributions # returns the choice probabilitiy for the focal actor
}

##@expectedChangeProbabilities. Use as RSiena:::expectedChangeProbabilities
expectedChangeProbabilities <- function(conts, effects, theta, thedata = NULL, 
                                        getChangeStatistics = FALSE, effectNames = NULL) {
  waves <- length(conts[[1]])
  effects <- effects[effects$include == TRUE,]
  noRate <- effects$type != "rate"
  effects <- effects[noRate,]
  if(sum(noRate) != length(theta)) {
    theta <- theta[noRate]
  }
  effectNa <- attr(conts,"effectNames")
  effectTypes <- attr(conts,"effectTypes")
  networkNames <- attr(conts,"networkNames")
  networkTypes <- attr(conts,"networkTypes")
  networkInteraction <- effects$interaction1
  effectIds <- paste(effectNa,effectTypes,networkInteraction, sep = ".")
  currentDepName <- ""
  depNumber <- 0
  for(eff in 1:length(effectIds)) {
    if(networkNames[eff] != currentDepName) {
      currentDepName <- networkNames[eff]
      actors <- length(conts[[1]][[1]][[1]])
      depNumber <- depNumber + 1
      currentDepEffs <- effects$name == currentDepName
      depNetwork <- thedata$depvars[[depNumber]]
      if (networkTypes[eff] == "oneMode") {
        choices <- actors
      }
      else if (networkTypes[eff] == "behavior") {
        choices <- 3
      }
      else if (networkTypes[eff] == "bipartite") {
        if (dim(depNetwork)[2] >= actors) {
          stop("does not work for bipartite networks with second mode >= first mode")
        }
        choices <- dim(depNetwork)[2] + 1
      }
      else {
        stop("does not work for dependent variables of type 'continuous'")
      }
      
      # impute for wave 1
      if (networkTypes[eff] %in% c("oneMode", "bipartite")) {
        depNetwork[,,1][is.na(depNetwork[,,1])] <- 0
      }
      else {
        depNetwork[,,1][is.na(depNetwork[,,1])] <- attr(depNetwork, 'modes')[1]
      }
      # impute for next waves;
      # this may be undesirable for structurals immediately followed by NA...
      for (m in 2:dim(depNetwork)[3]) {
        depNetwork[,,m][is.na(depNetwork[,,m])] <- depNetwork[,,m - 1][is.na(depNetwork[,,m])]
        }
      # Make sure the diagonals are not treated as structurals
      if (networkTypes[eff] == "oneMode") {
        for (m in 1:(dim(depNetwork)[3])) {
          diag(depNetwork[,,m]) <- 0
        }
      }
      structurals <- (depNetwork >= 10)
      if (networkTypes[eff] == "oneMode") {
        if (attr(depNetwork, 'symmetric')) {
          message('\nNote that for symmetric networks, effect sizes are for modelType 2 (forcing).')
        }
      }
      
      #			currentDepObjEffsNames <- paste(effects$shortName[currentDepEffs],
      #				effects$type[currentDepEffs],effects$interaction1[currentDepEffs],sep=".")
      #			otherObjEffsNames <- paste(effects$shortName[!currentDepEffs],
      #				effects$type[!currentDepEffs],effects$interaction1[!currentDepEffs],sep=".")
      
      changeStats <- list()
      if (networkTypes[eff] == "behavior") {
        toggleProbabilities <- array(0, dim = c(actors, 3, waves))
      }
      else {
        toggleProbabilities <- array(0, dim = c(actors, choices, waves))
      }
      for(w in 1:waves) {
        currentDepEffectContributions <- conts[[1]][[w]][currentDepEffs]
        if (networkTypes[eff] == "bipartite") {
          currentDepEffectContributions <- lapply(
            currentDepEffectContributions,
            function(x){lapply(x,function(xx){xx[1:choices]})}
            )
        }
        # conts[[1]] is periods by effects by actors by actors
        currentDepEffectContributions <-
          sapply(lapply(currentDepEffectContributions, unlist),
                 matrix, nrow = actors, ncol = choices, byrow = TRUE,
                 simplify = "array")
        cdec <- apply(currentDepEffectContributions, c(2,1), as.matrix)
        # cdec is effects by actors (alters) by actors (egos)
        if (dim(currentDepEffectContributions)[3] <= 1) { # only one effect
          cdec <- array(cdec, dim = c(1,dim(cdec)))
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
            cdec[ff,,][t(structurals[,,w])] <- NA
          }
        }
        distributions <- apply(cdec, 3,
                               calculateChoiceProbability, theta[which(currentDepEffs)]) # matrix manipulation should be more efficient solution
        if (networkTypes[eff] == "behavior") {
          toggleProbabilities[,,w] <- t(distributions)
        }
        else {
          toggleProbabilities[,,w] <- t(distributions)
        }
      }
    }
  }
  toggleProbabilities                                             
}