## /*****************************************************************************
##	* SIENA: Simulation Investigation for Empirical Network Analysis
##	*
##	* Web: https://www.stats.ox.ac.uk/~snijders/siena
##	*
##	* File: sienaGOF.r
##	*
##	* Description: This file contains the code to assess goodness of fit:
##	* the main function sienaGOF, the plot method,
##	* and auxiliary statistics and extractor functions.
##	* Written by Josh Lospinoso, modifications by Tom Snijders.
##	*
##	****************************************************************************/

##@sienaGOF siena07 Does test for goodness of fit
sienaGOF <- function(
		sienaFitObject,	auxiliaryFunction,
		period=NULL, verbose=FALSE, join=TRUE, twoTailed=FALSE,
		cluster=NULL, robust=FALSE,
		groupName="Data1", varName, tested=NULL, iterations=NULL,
		giveNAWarning=TRUE, ...)
	{
	## require(MASS)
	## require(Matrix)
	##	Check input
	# fitList indicates whether the sienaFitObject is a list of
	# sienaFit objects, or a single such object.
	fitList <- FALSE
	if (inherits(sienaFitObject, "sienaFit"))
	{
		sFO <- sienaFitObject
	}
	else
	{
		if (all(vapply(sienaFitObject, function(x){inherits(x, "sienaFit")}, FUN.VALUE=TRUE)))		
		{
			sFO <- sienaFitObject[[1]]
# If fitList, it is assumed that all elements of this list
# have the same specification, and the first one is used
# to deduce this information.
			fitList <- TRUE
			if (!inherits(sFO, "sienaFit"))
			{
				stop("The first parameter of sienaGOF must be a sienaFit object or a list of such objects")
			}	
			if (length(period) > 1)
			{
				stop("For operating on a list of sienaFit objects, only a single period can be used.")
			}
		}
	}
	
	if (sFO$maxlike)
	{
		stop(
	"sienaGOF can only operate on results from Method of Moments estimation.")
	}
	if (! sFO$returnDeps)
	{
		stop("You must instruct siena07 to return the simulated networks")
	}
	if (!is.null(sFO$sf2.byIteration))
	{
		if (!sFO$sf2.byIteration)
    	{
        	stop("sienaGOF needs sf2 by iterations (use lessMem=FALSE)")
    	}
	}
	if (is.null(iterations))
	{
		iterations <- length(sFO$sims)
	}
	else
	{
		iterations <- min(iterations, length(sFO$sims))
	}
	if (iterations < 1)
	{
		stop("You need at least one iteration.")
	}
	if (missing(varName))
	{
		stop("You need to supply the parameter <<varName>>.")
	}
	if (missing(auxiliaryFunction))
	{
		stop("You need to supply the parameter <<auxiliaryFunction>>.")
	}
# There might be more than one varName:
	if (is.null(sFO$f[[groupName]]$depvars[[varName[1]]]))
	{
		stop("There is a mismatch between the sienaFitObject and the groupName or varName.")
	}
	
	groups <- length(sFO$f$groupNames)
	if (fitList)
	{
		tested <- FALSE
	}
	else if (is.null(tested))
	{
		tested <- sFO$test
	}
	else
	{
		if (!inherits(tested, "logical"))
		{
			stop('tested should be a logical vector')
		}
		if ((length(tested) != length(sFO$test)) | (all(tested == FALSE)))
		{
			tested <- rep(FALSE, length(sFO$test))
		}
		else
		{
			tested <- (tested & sFO$test)
		}
	}
	if (verbose)
	{
		if (groups <= 1)
		{
			message("Detected ", iterations, " iterations and ", groups, " group.")
		}
		else
		{
			message("Detected ", iterations, " iterations and ", groups, " groups.")
		}
		if (fitList)
		{
			message("The data for analysis is a list of ", length(sienaFitObject),
						" sienaFit objects.")
		}
	}

	if (fitList)
	{
		if (is.null(period) )
		{
			period <- 1
		}
		per <- period # can only be a single number
		auxFunction <- function(i, sienaFitObject, j, groupName, varName, ...){
			auxiliaryFunction(i, sienaFitObject[[j]]$f,
						sienaFitObject[[j]]$sims, per, groupName, varName, ...)}
		period <- seq_along(sienaFitObject) 
		# from now on, period will denote the rank number of the sienaFit object.
	}	
	else
	{
		if (is.null(period) )
		{
			period <- 1:(attr(sienaFitObject$f[[1]]$depvars[[1]], "netdims")[3] - 1)
		}
		auxFunction <- function(i, sienaFitObject, j, groupName, varName, ...){
			auxiliaryFunction(i, sienaFitObject$f,
						sienaFitObject$sims, j, groupName, varName, ...)}
	}

	obsStatsByPeriod <- lapply(period, function (j) {
						matrix(
						auxFunction(NULL,
								sienaFitObject, j, groupName, varName, ...)
								, nrow=1)
				})

	if (join)
	{
		obsStats <- Reduce("+", obsStatsByPeriod)
		obsStats <- list(Joint=obsStats)
	}
	else
	{
		obsStats <- obsStatsByPeriod
		names(obsStats) <- paste("Period", period)
	}
	plotKey <- names(auxiliaryFunction(NULL, sFO$f,
				sFO$sims, 1, groupName, varName, ...))
	class(obsStats) <- "observedAuxiliaryStatistics"
	attr(obsStats,"auxiliaryStatisticName") <-
			deparse(substitute(auxiliaryFunction))
	attr(obsStats,"joint") <- join
	

	##	Calculate the simulated auxiliary statistics
	if (verbose)
	{
		if (length(period) <= 1)
		{
			cat("Calculating auxiliary statistics for period ", period, ".\n")
		}
		else
		{
			cat("Calculating auxiliary statistics for periods ", period, ".\n")
		}
	}
	if (!is.null(cluster))
	{
		ttcSimulation <- system.time(simStatsByPeriod <-
			lapply(period, function (j) {
				simStatsByPeriod <- parSapply(cluster, 1:iterations,
					function (i){auxFunction(i, sienaFitObject,
										j, groupName, varName, ...)})
				simStatsByPeriod <- matrix(simStatsByPeriod, ncol=iterations)
				dimnames(simStatsByPeriod)[[1]] <-	plotKey
				dimnames(simStatsByPeriod)[[2]] <-	1:iterations
				t(simStatsByPeriod)
				}))
	}
	else
	{
		ttcSimulation <- system.time( simStatsByPeriod <- lapply(period,
					function (j) {
						if (verbose)
						{
							message("  Period ", j, "\n")
							flush.console()
						}
						simStatsByPeriod <- sapply(1:iterations, function (i)
						{
							if (verbose && (i %% 100 == 0) )
								{
								cat("  > Completed ", i,
										" calculations\r")
								flush.console()
								}
								auxFunction(i, sienaFitObject,
										j, groupName, varName, ...)
						})
					if (verbose)
					{
						cat("  > Completed ", iterations, " calculations\n\n")
					}
					flush.console()
					simStatsByPeriod <-
							matrix(simStatsByPeriod, ncol=iterations)
					dimnames(simStatsByPeriod)[[1]] <-	plotKey
					dimnames(simStatsByPeriod)[[2]] <-	1:iterations
					t(simStatsByPeriod)
					})
	  )
	}

	## Give a warning in case of missings.
	nmissings <- vapply(simStatsByPeriod,
		function(sp){apply(sp, 2, function(x){sum(is.na(x))})},
					FUN.VALUE=rep(0,dim(simStatsByPeriod[[1]])[2]))
	rownames(nmissings) <- plotKey
	if ((sum(nmissings) > 0) & giveNAWarning)
	{
		cat("Number of missing values in the simulated functions:\n")
		print(t(nmissings))
		warning("Some simulated values are missing.")
	}

	## Aggregate by period if necessary to produce simStats
	if (join)
	{
			simStats <- Reduce("+", simStatsByPeriod)
		simStats <- list(Joint=simStats)
	}
	else
	{
		simStats <- simStatsByPeriod
		names(simStats) <- paste("Period",period)
	}
	class(simStats) <- "simulatedAuxiliaryStatistics"
	attr(simStats,"auxiliaryStatisticName") <-
			deparse(substitute(auxiliaryFunction))
	attr(simStats,"joint") <- join
	attr(simStats,"time") <- ttcSimulation

	applyTest <-  function (observed, simulated)
# Test using Mahalanobis distances
	{
		if (!inherits(simulated,"matrix"))
		{
			stop("Invalid input.")
		}
		if (!inherits(simulated,"matrix"))
		{
			observed <- matrix(observed,nrow=1)
		}
		if (!inherits(simulated,"matrix"))
		{
			stop("Observation must be a matrix.")
		}
		if (ncol(observed) != ncol(simulated))
		{
			stop("Dimensionality of function parameters do not match.")
		}
		observations <- nrow(observed)
	#	simulations<-nrow(simulated)
		variates<-ncol(simulated)
		if (robust) {
			a <- cov.rob(simulated)$cov
		}
		else
		{
			a <- cov(simulated, use="pairwise.complete.obs")
			a[is.na(a)] <- 0
		}
		ainv <- ginv(a)
		arank <- rankMatrix(a)
		expectation <- colMeans(simulated)
		centeredSimulations <- scale(simulated, scale=FALSE)
		if (variates==1)
		{
			centeredSimulations <- t(centeredSimulations)
		}
		mhd <- function(x) # Mahalanobis distance
		{
			x %*% ainv %*% x
		}
		simTestStat <- apply(centeredSimulations, 1, mhd)
		centeredObservations <- observed - expectation
		obsTestStat <- apply(centeredObservations, 1, mhd)
		if (twoTailed)
		{
			p <- sapply(1:observations, function (i)
						1 - abs(1 - 2 * sum(obsTestStat[i] <=
						simTestStat)/length(simTestStat)) )
		}
		else
		{
			p <- sapply(1:observations, function (i)
				sum(obsTestStat[i] <= simTestStat) /length(simTestStat))
		}

		ret <- list( p = p,
				SimulatedTestStat=simTestStat,
				ObservedTestStat=obsTestStat,
				TwoTailed=twoTailed,
				Simulations=simulated,
				Observations=observed,
				InvCovSimStats=a,
				Rank=arank)
		class(ret) <- "sienaGofTest"
		attr(ret, "version") <- packageDescription(pkgname, fields = "Version")
		attr(ret,"sienaFitName") <- deparse(substitute(sienaFitObject))
		attr(ret,"auxiliaryStatisticName") <-
				attr(obsStats,"auxiliaryStatisticName")
		attr(ret, "key") <- plotKey
		ret
	}

	res <- lapply(1:length(simStats),
					function (i) {
				 applyTest(obsStats[[i]], simStats[[i]]) })
	mhdTemplate <- rep(0, sum(tested))
	names(mhdTemplate) <- rep(0, sum(tested))

	JoinedOneStepMHD_old <- mhdTemplate
	OneStepMHD_old <- lapply(period, function(i) (mhdTemplate))
	JoinedOneStepMHD <- mhdTemplate
	OneStepMHD <- lapply(period, function(i) (mhdTemplate))

	obsMhd <- NULL

	ExpStat <-
		lapply(period, function(i) {colMeans(simStatsByPeriod[[i]])})
	simStatsByPeriod_tilde <-
		lapply(period, function(i) {
			t(apply(simStatsByPeriod[[i]],1, function(x){x - ExpStat[[i]]}))})

	OneStepSpecs <- matrix(0, ncol=sum(tested),
			nrow=length(sienaFitObject$theta))
	if (robust) {
		covInvByPeriod <- lapply(period, function(i) ginv(
							cov.rob(simStatsByPeriod[[i]]) ))
	}
	else
	{
		covInvByPeriod <- lapply(period, function(i){
				b <- cov(simStatsByPeriod[[i]], use="pairwise.complete.obs")
				b[is.na(b)] <- 0
				ginv(b)})
	}

	obsMhd <- sapply(period, function (i) {
				 (obsStatsByPeriod[[i]] - ExpStat[[i]])	 %*%
						covInvByPeriod[[i]] %*%
						t(obsStatsByPeriod[[i]] - ExpStat[[i]] )
			})

	if (sum(tested) > 0) {
		effectsObject <- sienaFitObject$requestedEffects
		for (i in period) {
			names(OneStepMHD_old[[i]]) <-
					effectsObject$effectName[tested]
			names(OneStepMHD[[i]]) <-
					effectsObject$effectName[tested]
		}
		names(JoinedOneStepMHD_old) <-
			effectsObject$effectName[tested]
		names(JoinedOneStepMHD) <-
				effectsObject$effectName[tested]

		rownames(OneStepSpecs) <- effectsObject$effectName
		colnames(OneStepSpecs) <- effectsObject$effectName[tested]
		counterTestEffects <- 0
		for(index in which(tested)) {
			if (verbose) {
				message("Estimating test statistic for model including ",
						effectsObject$effectName[index], "\n")
			}
			counterTestEffects <- counterTestEffects + 1
			effectsToInclude <- !tested
			effectsToInclude[index] <- TRUE
			theta0 <- sienaFitObject$theta
			names(theta0) <- effectsObject$effectName
			theta0 <- theta0[effectsToInclude]
			obsSuffStats <-
					t(sienaFitObject$targets2[effectsToInclude, , drop=FALSE])
			G <- sienaFitObject$sf2[, , effectsToInclude, drop=FALSE] -
					rep(obsSuffStats, each=iterations)
			sigma <- cov(apply(G, c(1, 3), sum), use="pairwise.complete.obs")
			SF <- sienaFitObject$ssc[ , , effectsToInclude, drop=FALSE]
			dimnames(SF)[[3]] <- effectsObject$effectName[effectsToInclude]
			dimnames(G) <- dimnames(SF)
			if (!(sienaFitObject$maxlike || sienaFitObject$FinDiff.method))
			{
				D <- derivativeFromScoresAndDeviations(SF, G, , , , TRUE, )
			}
			else
			{
				DF <- sienaFitObject$
						sdf2[ , , effectsToInclude, effectsToInclude,
						drop=FALSE]
				D <- t(apply(DF, c(3, 4), mean))
			}
			fra <- apply(G, 3, sum) / iterations
			doTests <- rep(FALSE, sum(effectsToInclude))
			names(doTests) <- effectsObject$effectName[effectsToInclude]
			doTests[effectsObject$effectName[index]] <- TRUE
			redundant <- rep(FALSE, length(doTests))
			mmThetaDelta <- as.numeric(ScoreTest(length(doTests), D,
							sigma, fra, doTests, redundant,
							maxlike=sienaFitObject$maxlike)$oneStep )

      # \mu'_\theta(X)	  
			JacobianExpStat_old <- lapply(period, function (i) {
				t(SF[1:iterations,i,]) %*% simStatsByPeriod[[i]]/ iterations	 })
			JacobianExpStat <- lapply(period, function (i) {
				t(SF[1:iterations,i,]) %*% simStatsByPeriod_tilde[[i]]/ iterations })
      # List structure: Period, effect index
      thetaIndices <- 1:sum(effectsToInclude)
	  # \Gamma_i(\theta)  i=period, j=parameter, k=replication
			ExpStatCovar_old <- lapply(period, function (i) {
            lapply(thetaIndices, function(j){
              Reduce("+", lapply(1:iterations,function(k){
                simStatsByPeriod[[i]][k,] %*% t(simStatsByPeriod[[i]][k,]) * SF[k,i,j]
              })) / iterations
				- JacobianExpStat[[i]][j,] %*%
			t(ExpStat[[i]]) - ExpStat[[i]] %*% t(JacobianExpStat[[i]][j,])
            })
        })
			ExpStatCovar <- lapply(period, function (i) {
				lapply(thetaIndices, function(j){
				Reduce("+", lapply(1:iterations,function(k){
			simStatsByPeriod_tilde[[i]][k,] %*%
				t(simStatsByPeriod_tilde[[i]][k,]) * SF[k,i,j] })) / iterations})})
      # \Xi_i(\theta)
			JacobianCovar_old <- lapply(period, function (i) {
				lapply(thetaIndices, function(j){
					-1 * covInvByPeriod[[i]] %*% ExpStatCovar_old[[i]][[j]] %*%
						covInvByPeriod[[i]] })
			})
      JacobianCovar <- lapply(period, function (i) {
        lapply(thetaIndices, function(j){
					-1 * covInvByPeriod[[i]] %*% ExpStatCovar[[i]][[j]] %*%
						covInvByPeriod[[i]] })
        })

			Gradient_old <- lapply(period, function(i) {
				sapply(thetaIndices, function(j){
					( obsStatsByPeriod[[i]] - ExpStat[[i]] ) %*%
						JacobianCovar_old[[i]][[j]] %*%
					t( obsStatsByPeriod[[i]] - ExpStat[[i]] )
					})
				-2 * JacobianExpStat_old[[i]] %*% covInvByPeriod[[i]] %*%
					t( obsStatsByPeriod[[i]] - ExpStat[[i]] )
				})
			Gradient <- lapply(period, function(i) {
          sapply(thetaIndices, function(j){
          ( obsStatsByPeriod[[i]] - ExpStat[[i]] ) %*%
            JacobianCovar[[i]][[j]] %*%
          t( obsStatsByPeriod[[i]] - ExpStat[[i]] )
          })
				-2 * JacobianExpStat[[i]] %*% covInvByPeriod[[i]] %*%
            t( obsStatsByPeriod[[i]] - ExpStat[[i]] )
					})

			OneStepSpecs[effectsToInclude,counterTestEffects] <-
								theta0 + mmThetaDelta

			for (i in 1:length(obsMhd)) {
				OneStepMHD_old[[i]][counterTestEffects] <-
					as.numeric(obsMhd[i] + mmThetaDelta %*% Gradient_old[[i]] )
      }
			for (i in 1:length(obsMhd)) {
				OneStepMHD[[i]][counterTestEffects] <-
					as.numeric(obsMhd[i] + mmThetaDelta %*% Gradient[[i]] )
		}
			JoinedOneStepMHD_old[counterTestEffects] <-
						Reduce("+",OneStepMHD_old)[counterTestEffects]
			JoinedOneStepMHD[counterTestEffects] <-
						Reduce("+",OneStepMHD)[counterTestEffects]
		} # end 'for index'
	}

	names(res) <- names(obsStats)
	class(res) <- "sienaGOF"
	attr(res, "version") <- packageDescription(pkgname, fields = "Version")
	attr(res, "scoreTest") <- (sum(tested) > 0)
	attr(res, "originalMahalanobisDistances") <- obsMhd
	attr(res, "oneStepMahalanobisDistances") <- OneStepMHD
	attr(res, "joinedOneStepMahalanobisDistances") <-
			JoinedOneStepMHD
	attr(res, "oneStepMahalanobisDistances_old") <- OneStepMHD_old
	attr(res, "joinedOneStepMahalanobisDistances_old") <-
			JoinedOneStepMHD_old
	attr(res, "oneStepSpecs") <- OneStepSpecs
	attr(res,"auxiliaryStatisticName") <-
			attr(obsStats,"auxiliaryStatisticName")
	attr(res, "simTime") <- attr(simStats,"time")
	attr(res, "twoTailed") <- twoTailed
	attr(res, "joined") <- join
	attr(res, "nmissings") <- nmissings
	res
}

##@print.sienaGOF siena07 Print method for sienaGOF
print.sienaGOF <- function (x, ...) {
	## require(Matrix)
	levls <- 1:length(x)
	pVals <- sapply(levls, function(i) x[[i]]$p)
	titleStr <- "Monte Carlo Mahalanobis distance test p-value: "

	if (! attr(x,"joined"))
	{
		cat("Siena Goodness of Fit (",
			attr(x,"auxiliaryStatisticName"),"),", length(levls)," periods\n=====\n")
		cat(" >",titleStr, "\n")
		for (i in 1:length(pVals))
		{
			cat(names(x)[i], ": ", round(pVals[i],3), "\n")
		}
		for (i in 1:length(pVals))
		{
			if (x[[i]]$Rank < dim(x[[i]]$Observations)[2])
			{
				cat(" * Note for", names(x)[i],
					": Only", x[[i]]$Rank, "statistics are",
					"necessary in the auxiliary function.\n")
			}
		}
	}
	else
	{
		cat("Siena Goodness of Fit (",
			attr(x,"auxiliaryStatisticName"),"), all periods\n=====\n")
		cat(titleStr, round(pVals[1],3), "\n")
		if (x[[1]]$Rank < dim(x[[1]]$Observations)[2])
			{
				cat("**Note: Only", x[[1]]$Rank, "statistics are",
				"necessary in the auxiliary function.\n")
			}
	}

	if ( attr(x, "twoTailed") )
	{
		cat("-----\nTwo tailed test used.")
	}
	else
	{
		cat("-----\nOne tailed test used ",
		"(i.e. estimated probability of greater distance than observation).\n")
	}
	originalMhd <- attr(x, "originalMahalanobisDistances")
	if (attr(x, "joined")) {
		cat("-----\nCalculated joint MHD = (",
				round(sum(originalMhd),2),") for current model.\n")
	}
	else
	{
		for (j in 1:length(originalMhd)) {
			cat("-----\nCalculated period ", j, " MHD = (",
					round(originalMhd[j],2),") for current model.\n")
		}
	}
	invisible(x)
}

##@summary.sienaGOF siena07 Summary method for sienaGOF
summary.sienaGOF <- function(object, ...) {
	x <- object
	print(x)
	if (sum(attr(x, "nmissings"))> 0){
		cat("\nThere were missing values in the simulated statistics.\n")
		cat("Their number (by period):\n")
		print(t(attr(x, "nmissings")))
	}
	if (attr(x, "scoreTest")) {
		oneStepSpecs <- attr(x, "oneStepSpecs")
		oneStepMhd <- attr(x, "oneStepMahalanobisDistances")
		joinedOneStepMhd <- attr(x, "joinedOneStepMahalanobisDistances")
		cat("\nOne-step estimates and predicted Mahalanobis distances")
		cat(" for modified models.\n")
		if (attr(x, "joined")) {
			for (i in 1:ncol(oneStepSpecs)) {
				a <- cbind(oneStepSpecs[,i, drop=FALSE] )
				b <- matrix( c(joinedOneStepMhd[i] ), ncol=1)
				rownames(b) <- c("MHD")
				a <- rbind(a, b)
				a <- round(a, 3)
				cat("\n**Model including", colnames(a)[1], "\n")
				colnames(a) <- "one-step"
				print(a)
			}
		}
		else
		{
			for (j in 1:length(oneStepMhd)) {
				for (i in 1:ncol(oneStepSpecs)) {
					a <- cbind( oneStepSpecs[,i, drop=FALSE] )
					b <- matrix( oneStepMhd[[j]][i], ncol=1 )
					rownames(b) <- c("MHD")
					a <- rbind(a, b)
					a <- round(a, 3)
					cat("\n**Model including", colnames(a)[1], "\n")
					colnames(a) <- c("one-step")
					print(a)
				}
			}
		}
		cat("\n-----")
	}
	cat("\nComputation time for auxiliary statistic calculations on simulations: ",
			attr(x, "simTime")["elapsed"] , "seconds.\n")
	invisible(x)
}

##@plot.sienaGOF siena07 Plot method for sienaGOF
plot.sienaGOF <- function (x, center=FALSE, scale=FALSE, violin=TRUE,
		key=NULL, perc=.05, period=1, position=4, fontsize=12, ...)
{
	## require(lattice)
	args <- list(...)
	if (is.null(args$main))
	{
		main=paste("Goodness of Fit of",
				attr(x,"auxiliaryStatisticName"))
		if (!attr(x,"joined"))
		{
			main = paste(main, "Period", period)
		}
	}
	else
	{
		main=args$main
	}

	if (attr(x,"joined"))
	{
		x <- x[[1]]
	}
	else
	{
		x <- x[[period]]
	}
	sims <- x$Simulations
	obs <- x$Observations
	itns <- nrow(sims)
#	vars <- ncol(sims)
	## Need to check for useless statistics here:
	n.obs <- nrow(obs)

	screen <- sapply(1:ncol(obs),function(i){
						(sum(is.nan(rbind(sims,obs)[,i])) == 0) }) &
				(diag(var(rbind(sims,obs)))!=0)

	if (any((diag(var(rbind(sims,obs)))==0)))
	{	cat("Note: some statistics are not plotted because their variance is 0.\n")
		cat("This holds for the statistic")
		if (sum(diag(var(rbind(sims,obs)))==0) > 1){cat("s")}
		cat(": ")
		cat(paste(attr(x,"key")[which(diag(var(rbind(sims,obs)))==0)], sep=", "))
		cat(".\n")
	}

	sims <- sims[,screen, drop=FALSE]
	obs <- obs[,screen, drop=FALSE]
	obsLabels <- round(x$Observations[,screen, drop=FALSE],3)

	sims.min <- apply(sims, 2, min)
	sims.max <- apply(sims, 2, max)
	sims.min <- pmin(sims.min, obs)
	sims.max <- pmax(sims.max, obs)

	if (center)
	{
		sims.median <- apply(sims, 2, median)
		sims <- sapply(1:ncol(sims), function(i)
					(sims[,i] - sims.median[i]) )
		obs <- matrix(sapply(1:ncol(sims), function(i)
							(obs[,i] - sims.median[i])), nrow=n.obs )
		sims.min <- sims.min - sims.median
		sims.max <- sims.max - sims.median
	}
	if (scale)
	{
		sims.range <- sims.max - sims.min + 1e-6
		sims <- sapply(1:ncol(sims), function(i) sims[,i]/(sims.range[i]))
		obs <- matrix(sapply(1:ncol(sims), function(i) obs[,i]/(sims.range[i]))
				, nrow=n.obs )
		sims.min <- sims.min/sims.range
		sims.max <- sims.max/sims.range
	}

	ymin <- 1.05*min(sims.min) - 0.05*max(sims.max)
	ymax <- -0.05*min(sims.min) + 1.05*max(sims.max)

	if (is.null(args$ylab))
	{
		ylabel = "Statistic"
		if (center && scale) {
			ylabel = "Statistic (centered and scaled)"
		}
		else if (scale)
		{
			ylabel = "Statistic (scaled)"
		}
		else if (center)
		{
			ylabel = "Statistic (center)"
		}
		else
		{
			ylabel = "Statistic"
		}
	}
	else
	{
		ylabel = args$ylab
	}

	if (is.null(args$xlab))
	{
		xlabel = paste( paste("p:", round(x$p, 3),
						collapse = " "), collapse = "\n")
	}
	else
	{
		xlabel = args$xlab
	}

	if (is.null(args$cex))
	{
		cexpar <- par("cex")
	}
	else
	{
		cexpar <- args$cex
	}

	if (is.null(args$cex.axis))
	{
		cexaxispar <- par("cex.axis")
	}
	else
	{
		cexaxispar <- args$cex.axis
	}

	if (is.null(args$cex.main))
	{
		cexmainpar <- par("cex.main")
	}
	else
	{
		cexmainpar <- args$cex.main
	}

	if (is.null(args$cex.lab))
	{
		cexlabpar <- par("cex.lab")
	}
	else
	{
		cexlabpar <- args$cex.lab
	}

	if (is.null(args$cex.sub))
	{
		cexsubpar <- par("cex.sub")
	}
	else
	{
		cexsubpar <- args$cex.sub
	}

	xAxis <- (1:sum(screen))

	if (is.null(key))
	{
		if (is.null(attr(x, "key")))
		{
			key=xAxis
		}
		else
		{
			key <- attr(x,"key")[screen]
		}
	}
	else
	{
		key <- key[screen] ## added 1.1-244
		if (length(key) != ncol(obs))
		{
			stop("Key length does not match the number of variates.")
		}
	}

	br <- trellis.par.get("box.rectangle")
	br$col <- 1
	trellis.par.set("box.rectangle", br)
	bu <- trellis.par.get("box.umbrella")
	bu$col <- 1
	trellis.par.set("box.umbrella", bu)
	plot.symbol <- trellis.par.get("plot.symbol")
	plot.symbol$col <- "black"
	plot.symbol$pch <- 4
	plot.symbol$cex <- cexpar  # default 1
	trellis.par.set("plot.symbol", plot.symbol)

#	plot.axis <- trellis.par.get("axis.text")
#	plot.axis$cex <- cexaxispar # default 1
	trellis.par.set("axis.text", list(cex=cexaxispar))

#	plot.xlab <- trellis.par.get("par.xlab.text")
#	plot.xlab$cex <- cexlabpar # default 1
	trellis.par.set("par.xlab.text", list(cex=cexlabpar))

#	plot.ylab <- trellis.par.get("par.ylab.text")
#	plot.ylab$cex <- cexlabpar # default 1
	trellis.par.set("par.ylab.text", list(cex=cexlabpar))

#	plot.main <- trellis.par.get("par.main.text")
#	plot.main$cex <- cexmainpar # default 1.2
	trellis.par.set("par.main.text", list(cex=cexmainpar))

#	plot.font <- trellis.par.get("fontsize")
#	plot.font$text <- fontsize
	trellis.par.set("fontsize", list(text=fontsize))

	panelFunction <- function(..., x=x, y=y, box.ratio){
		ind.lower <- max( round(itns * perc/2), 1)
		ind.upper <- round(itns * (1-perc/2))
		yperc.lower <- sapply(1:ncol(sims), function(i)
					sort(sims[,i])[ind.lower]  )
		yperc.upper <- sapply(1:ncol(sims), function(i)
					sort(sims[,i])[ind.upper]  )
		if (violin) {
			panel.violin(x, y, box.ratio=box.ratio, col = "transparent", ...)
		}
		panel.bwplot(x, y, box.ratio=.1, fill = "gray", ...)
		panel.xyplot(xAxis, yperc.lower, lty=3, col = "gray", lwd=3, type="l",
				...)
		panel.xyplot(xAxis, yperc.upper, lty=3, col = "gray", lwd=3, type="l",
				...)
		for(i in 1:nrow(obs))
		{
			panel.xyplot(xAxis, obs[i,],  col="red", type="l", lwd=1, ...)
			panel.xyplot(xAxis, obs[i,],  col="red", type="p", lwd=3, pch=19,
					...)
			panel.text(xAxis, obs[i,], labels=obsLabels[i,], pos=position)
		}
	}
	bwplot(as.numeric(sims)~rep(xAxis, each=itns), horizontal=FALSE,
			panel = panelFunction, xlab=xlabel, ylab=ylabel, ylim=c(ymin,ymax),
			scales=list(x=list(labels=key), y=list(draw=FALSE)),
			main=main)
}

##@descriptives.sienaGOF siena07 Gives numerical values in the plot.
descriptives.sienaGOF <- function (x, center=FALSE, scale=FALSE,
			perc=.05, key=NULL, period=1, showAll=FALSE)
{
# adapted excerpt from plot.sienaGOF
	if (attr(x,"joined"))
	{
		x <- x[[1]]
	}
	else
	{
		x <- x[[period]]
	}

	sims <- x$Simulations
	obs <- x$Observations
	itns <- nrow(sims)
	screen <- sapply(1:ncol(obs),function(i){
						(sum(is.nan(rbind(sims,obs)[,i])) == 0) })
	if (!showAll)
	{
		screen <- screen & (diag(var(rbind(sims,obs)))!=0)
	}
	sims <- sims[,screen, drop=FALSE]
	obs <- obs[,screen, drop=FALSE]
	## Need to check for useless statistics here:
	n.obs <- nrow(obs)

	if (is.null(key))
	{
		if (is.null(attr(x, "key")))
		{
			key=(1:sum(screen))
		}
		else
		{
			key <- attr(x,"key")[screen]
		}
	}
	else
	{
		if (length(key) != ncol(obs))
		{
			stop("Key length does not match the number of variates.")
		}
		key <- key[screen]
	}

	sims.themin <- apply(sims, 2, min, na.rm=TRUE)
	sims.themax <- apply(sims, 2, max, na.rm=TRUE)
	sims.mean <- apply(sims, 2, mean, na.rm=TRUE)
	sims.sd <- apply(sims, 2, sd, na.rm=TRUE)
	sims.min <- pmin(sims.themin, obs)
	sims.max <- pmax(sims.themax, obs)

	if (center)
	{
		sims.median <- apply(sims, 2, median)
		sims <- sapply(1:ncol(sims), function(i)
					(sims[,i] - sims.median[i]) )
		obs <- matrix(sapply(1:ncol(sims), function(i)
							(obs[,i] - sims.median[i])), nrow=n.obs )
		sims.mean <- sims.mean - sims.median
		sims.min <- sims.min - sims.median
		sims.max <- sims.max - sims.median
	}

	if (scale)
	{
		sims.range <- sims.max - sims.min + 1e-6
		sims <- sapply(1:ncol(sims), function(i) sims[,i]/(sims.range[i]))
		obs <- matrix(sapply(1:ncol(sims), function(i) obs[,i]/(sims.range[i]))
				, nrow=n.obs )
		sims.mean <- sims.mean/sims.range
		sims.min <- sims.min/sims.range
		sims.max <- sims.max/sims.range
		sims.sd <- sims.sd/sims.range
	}

	screen <- sapply(1:ncol(obs),function(i){
						(sum(is.nan(rbind(sims,obs)[,i])) == 0) })
	if (!showAll)
	{
		screen <- screen & (diag(var(rbind(sims,obs)))!=0)
	}
	sims <- sims[,screen, drop=FALSE]
	obs <- obs[,screen, drop=FALSE]
	sims.themin <- sims.themin[screen, drop=FALSE]
	sims.themax <- sims.themax[screen, drop=FALSE]

	ind.lower = max( round(itns * perc/2), 1)
	ind.upper = round(itns * (1-perc/2))
	ind.median = round(itns * 0.5)
	yperc.mid = sapply(1:ncol(sims), function(i)
				sort(sims[,i])[ind.median])
	yperc.lower = sapply(1:ncol(sims), function(i)
				sort(sims[,i])[ind.lower]  )
	yperc.upper = sapply(1:ncol(sims), function(i)
				sort(sims[,i])[ind.upper]  )
	ypg <- sapply(1:ncol(sims), function(i)	mean(sims[,i] > obs[1,i]))
	ypp <- sapply(1:ncol(sims), function(i)	mean(sims[,i] >= obs[1,i]))
    violins <- matrix(NA, 10, ncol(sims))
	violins[1,] <- sims.themax
	violins[2,] <- yperc.upper
	violins[3,] <- sims.mean
	violins[4,] <- yperc.mid
	violins[5,] <- yperc.lower
	violins[6,] <- sims.themin
	violins[7,] <- sims.sd
	violins[8,] <- obs
    violins[9, ] <- ypg
    violins[10, ] <- ypp
    rownames(violins) <- c("max", "perc.upper", "mean", "median",
        "perc.lower", "min", "sd", "obs", "p>", "p>=")
	colnames(violins) <- key
	violins
}

##@changeToStructural sienaGOF Utility to change
# values in X to structural values in S
# X must have values 0, 1.
# NA values in X will be 0 in the result.
changeToStructural <- function(X, S) {
	if (any(S >= 10, na.rm=TRUE))
		{
			S[is.na(S)] <- 0
			S0 <- (S==10)
			S1 <- (S==11)
# the 1* turns the logical into numeric
			X <- 1*((X - S0 + S1)>=1)
		}
	X[is.na(X)] <- 0
	drop0(X)
}

##@changeToNewStructural sienaGOF Utility to change
# values in X to structural values in SAfter
# for tie variables that have no structural values in SBefore.
# X must have values 0, 1.
# NA values in X or SBefore or SAfter will be 0 in the result.
changeToNewStructural <- function(X, SBefore, SAfter) {
		SB <- (SBefore>=10)
		SA <- (SAfter>=10)
		difAB <- (SA > SB)
		if (any(difAB, na.rm=TRUE))
		{
			S0 <- (difAB)*(SAfter==10)
			S1 <- (difAB)*(SAfter==11)
# the 1* turns the logical into numeric
			X <- 1*((X - S0 + S1)>=1)
		}
	X[is.na(X)] <- 0
	drop0(X)
}

##@sparseMatrixExtraction sienaGOF Extracts simulated networks
# This function returns the simulated network as a dgCMatrix;
# this is the "standard" class for sparse numeric matrices
# in the Matrix package. See the help file for "dgCMatrix-class".
# Ties for ordered pairs with a missing value for wave=period or period+1
# are zeroed;
# note that this also is done in RSiena for calculation of target statistics.
# To obtain equality between observed and simulated tie values
# in the case of structurally determined values, the following is done.
# The difficulty lies in the possibility
# that there is change in structural values.
# The reasoning is as follows:
# structural values affect the following period.
# Therefore the simulated values at the end of the period
# should be compared with an observation containing the structural values
# present at the beginning of the period.
# This implies that observations (wave=period+1) should be modified to contain
# the structural values of the preceding observation (wave=period).
# But if there are any tie variables with
# structural values for wave=period+1 and free values for wave=period,
# then there is no valid reference value for the simulations in this period,
# and the simulated tie values should be set to
# the observed (structural) values for wave=period+1.
# Concluding:
# For ties that have a structurally determined value at wave=period,
# this value is used for the observation at the end of the period.
# For ties that have a structurally determined value at the end of the period
# and a free value at the start,
# the structurally determined value at wave=period+1 is used
# for the simulations at the end of the period.
# TODO: Calculate the matrix of structurals and of missings outside
# of this procedure, doing it only once. Perhaps in sienaGOF.
sparseMatrixExtraction <-
	function(i, obsData, sims, period, groupName, varName){
	# require(Matrix)
	isBipartite <- "bipartite" == attr(obsData[[groupName]]$depvars[[varName]], "type")
	dimsOfDepVar<- attr(obsData[[groupName]]$depvars[[varName]], "netdims")
	if (attr(obsData[[groupName]]$depvars[[varName]], "missing"))
	{
		if (attr(obsData[[groupName]]$depvars[[varName]], "sparse"))
		{
			missings <-
				(is.na(obsData[[groupName]]$depvars[[varName]][[period]]) |
				is.na(obsData[[groupName]]$depvars[[varName]][[period+1]]))*1
		}
		else
		{
			missings <- Matrix(
				(is.na(obsData[[groupName]]$depvars[[varName]][,,period]) |
				is.na(obsData[[groupName]]$depvars[[varName]][,,period+1]))*1)
		}
	}
	if (is.null(i))
	{
		# sienaGOF wants the observation;
		# transform structurally fixed values into regular values
		# by "modulo 10" (%%10) operation
		# If preceding observation contains structural values
		# use these to replace the observations at period+1.
		if (attr(obsData[[groupName]]$depvars[[varName]], "sparse"))
		{
			extractedMatrix <- drop0(Matrix(
				obsData[[groupName]]$depvars[[varName]][[period+1]] %% 10))
			extractedMatrix[is.na(extractedMatrix)] <- 0
			if (attr(obsData[[groupName]]$depvars[[varName]], "structural"))
			{
				extractedMatrix <- changeToStructural(extractedMatrix,
					Matrix(obsData[[groupName]]$depvars[[varName]][[period]]))
			}
		}
		else # not sparse
		{
			extractedMatrix <-
			 Matrix(obsData[[groupName]]$depvars[[varName]][,,period+1] %% 10)
			extractedMatrix[is.na(extractedMatrix)] <- 0
			if (attr(obsData[[groupName]]$depvars[[varName]], "structural"))
			{
				extractedMatrix <- changeToStructural(extractedMatrix,
				Matrix(obsData[[groupName]]$depvars[[varName]][,,period]))
			}
		}
		if(!isBipartite){ diag(extractedMatrix) <- 0} # not guaranteed by data input
	}
	else
	{
		# sienaGOF wants the i-th simulation:
		extractedMatrix <- sparseMatrix(
				sims[[i]][[groupName]][[varName]][[period]][,1],
				sims[[i]][[groupName]][[varName]][[period]][,2],
				x=sims[[i]][[groupName]][[varName]][[period]][,3],
				dims=dimsOfDepVar[1:2], repr="T")
		if (attr(obsData[[groupName]]$depvars[[varName]], "structural"))
		{
		# If observation at end of period contains structural values
		# use these to replace the simulations.
			if (attr(obsData[[groupName]]$depvars[[varName]], "sparse"))
			{
				extractedMatrix <- changeToNewStructural(extractedMatrix,
				Matrix(obsData[[groupName]]$depvars[[varName]][[period]]),
				Matrix(obsData[[groupName]]$depvars[[varName]][[period+1]]))
			}
			else # not sparse
			{
				extractedMatrix <- changeToNewStructural(extractedMatrix,
				Matrix(obsData[[groupName]]$depvars[[varName]][,,period]),
				Matrix(obsData[[groupName]]$depvars[[varName]][,,period+1]))
			}
		}
	}
	## Zero missings (the 1* turns the logical into numeric):
	if (attr(obsData[[groupName]]$depvars[[varName]], "missing"))
	{
		extractedMatrix <- 1*drop0((extractedMatrix - missings) > 0)
	}
	extractedMatrix
}

##@sparseMatrixExtraction0 sienaGOF Extracts simulated networks
## simplified version of sparseMatrixExtraction if there are no missings or structurals
sparseMatrixExtraction0 <-
	function(i, obsData, sims, period, groupName, varName){
	# require(Matrix)
	isBipartite <- "bipartite" == attr(obsData[[groupName]]$depvars[[varName]], "type")
	dimsOfDepVar<- attr(obsData[[groupName]]$depvars[[varName]], "netdims")
	if (is.null(i))
	{
		# sienaGOF wants the observation:
		if (attr(obsData[[groupName]]$depvars[[varName]], "sparse"))
		{
			extractedValue <- Matrix(
				obsData[[groupName]]$depvars[[varName]][[period+1]])
			extractedValue[is.na(extractedValue)] <- 0
		}
		else # not sparse
		{
			extractedValue <-
			 Matrix(obsData[[groupName]]$depvars[[varName]][,,period+1])
			extractedValue[is.na(extractedValue)] <- 0
		}
		diag(extractedValue) <- 0 # not guaranteed by data input
		extractedValue <- as(drop0(extractedValue), "sparseMatrix")
	}
	else
	{
		# sienaGOF wants the i-th simulation:
		extractedValue <- sparseMatrix(
				sims[[i]][[groupName]][[varName]][[period]][,1],
				sims[[i]][[groupName]][[varName]][[period]][,2],
				x=1,
				dims=dimsOfDepVar[1:2], repr="T")
	}
	extractedValue
}

##@networkExtraction sienaGOF Extracts simulated networks
# This function provides a standard way of extracting simulated and observed
# networks from the results of a siena07 run.
# It returns the network as an edge list of class "network"
# according to the <network> package (used for package sna).
# Ties for ordered pairs with a missing value for wave=period or period+1
# are zeroed;
# note that this also is done in RSiena for calculation of target statistics.
# Structural values are treated as in sparseMatrixExtraction.
networkExtraction <- function (i, obsData, sims, period, groupName, varName){
	## suppressPackageStartupMessages(require(network))
	dimsOfDepVar<- attr(obsData[[groupName]]$depvars[[varName]], "netdims")
	isbipartite <- (attr(obsData[[groupName]]$depvars[[varName]], "type")
						=="bipartite")
	# For bipartite networks in package <network>,
	# the number of nodes is equal to
	# the number of actors (rows) plus the number of events (columns)
	# with all actors preceding all events.
	# Therefore the bipartiteOffset will come in handy:
	bipartiteOffset <- ifelse(isbipartite, dimsOfDepVar[1], 0)

	# Initialize empty networks:
	if (isbipartite)
	{
		emptyNetwork <- network::network.initialize(dimsOfDepVar[1]+dimsOfDepVar[2],
											bipartite=dimsOfDepVar[1])
	}
	else
	{
		emptyNetwork <- network::network.initialize(dimsOfDepVar[1], bipartite=NULL)
	}
	# Use what was defined in the function above:
	matrixNetwork <- sparseMatrixExtraction(i, obsData, sims,
						period, groupName, varName)
	if (sum(matrixNetwork) <= 0) # else network.edgelist() below will not work
	{
		extractedValue <- emptyNetwork
	}
	else
	{
		tripEV <- mat2triplet(matrixNetwork)
		extractedValue <- network::network.edgelist(
					cbind(tripEV$i, tripEV$j + bipartiteOffset, 1),
					emptyNetwork)
	}
	extractedValue
}

##@behaviorExtraction sienaGOF Extracts simulated behavioral variables.
# This function provides a standard way of extracting simulated and observed
# dependent behavior variables from the results of a siena07 run.
# The result is an integer vector.
# Values for actors with a missing value for wave=period or period+1 are
# transformed to NA.
behaviorExtraction <- function (i, obsData, sims, period, groupName, varName) {
  missings <- is.na(obsData[[groupName]]$depvars[[varName]][,,period]) |
	is.na(obsData[[groupName]]$depvars[[varName]][,,period+1])
  if (is.null(i))
	{
		# sienaGOF wants the observation:
		original <- obsData[[groupName]]$depvars[[varName]][,,period+1]
		original[missings] <- NA
		extractedValue <- original
	}
	else
	{
		#sienaGOF wants the i-th simulation:
		extractedValue <- sims[[i]][[groupName]][[varName]][[period]]
		extractedValue[missings] <- NA
	}
	extractedValue
}

##@OutdegreeDistribution sienaGOF Calculates Outdegree distribution
OutdegreeDistribution <- function(i, obsData, sims, period, groupName, varName,
						levls=0:8, cumulative=TRUE) {
	if (!(attr(obsData[[groupName]]$depvars[[varName]],'missing') |
	      attr(obsData[[groupName]]$depvars[[varName]],'structural')))
	{
		x <- sparseMatrixExtraction0(i, obsData, sims, period, groupName, varName)
	}
	else
	{
	x <- sparseMatrixExtraction(i, obsData, sims, period, groupName, varName)
	}

	a <- apply(x, 1, sum)
	if (cumulative)
	{
		oddi <- sapply(levls, function(i){ sum(a<=i) })
	}
	else
	{
		oddi <- sapply(levls, function(i){ sum(a==i) })
	}
	names(oddi) <- as.character(levls)
	oddi
}

##@IndegreeDistribution sienaGOF Calculates Indegree distribution
IndegreeDistribution <- function (i, obsData, sims, period, groupName, varName,
						levls=0:8, cumulative=TRUE){
	if (!(attr(obsData[[groupName]]$depvars[[varName]],'missing') |
	      attr(obsData[[groupName]]$depvars[[varName]],'structural')))
	{
		x <- sparseMatrixExtraction0(i, obsData, sims, period, groupName, varName)
	}
	else
	{
  x <- sparseMatrixExtraction(i, obsData, sims, period, groupName, varName)
	}
  a <- apply(x, 2, sum)
  if (cumulative)
  {
	iddi <- sapply(levls, function(i){ sum(a<=i) })
  }
  else
  {
	iddi <- sapply(levls, function(i){ sum(a==i) })
  }
  names(iddi) <- as.character(levls)
  iddi
}

##@BehaviorDistribution sienaGOF Calculates behavior distribution
BehaviorDistribution <- function (i, obsData, sims, period, groupName, varName,
							levls=NULL, cumulative=TRUE){
	x <- behaviorExtraction(i, obsData, sims, period, groupName, varName)
	if (is.null(levls))
	{
		levls <- attr(obsData[[groupName]]$depvars[[varName]],"behRange")[1]:
					attr(obsData[[groupName]]$depvars[[varName]],"behRange")[2]
	}
	if (cumulative)
	{
		bdi <- sapply(levls, function(i){ sum(x<=i, na.rm=TRUE) })
	}
	else
	{
	bdi <- sapply(levls, function(i){ sum(x==i, na.rm=TRUE) })
	}
	names(bdi) <- as.character(levls)
	bdi
}

##@mixedTriadCensus sienaGOF Calculates mixed triad census
# Contributed by Christoph Stadtfeld.
# For more details see
# Hollway, J., Lomi, A., Pallotti, F., & Stadtfeld, C. (2017)
# Multilevel social spaces: The network dynamics of organizational fields
# Network Science, 5(2), 187-212. doi:10.1017/nws.2017.8
#
# https://www.cambridge.org/core/journals/network-science/article/multilevel-social-spaces-the-network-dynamics-of-organizational-fields/602BB810A44497EBDE2A111A6C2771A3
#
# The function is called with two varName parameters, e.g.
# gof <- sienaGOF(res, mixedTriadCensus, varName = c("oneModeNet", "twoModeNet"))
#
mixedTriadCensus <- function (i, obsData, sims, period, groupName, varName) {
  if (length(varName) != 2) stop("mixedTriadCensus expects two varName parameters")
  varName1 <- varName[1]
  varName2 <- varName[2]

  # get matrices
  m1 <- as.matrix(sparseMatrixExtraction(i, obsData, sims, period, groupName,
											varName = varName1))
  m2 <- as.matrix(sparseMatrixExtraction(i, obsData, sims, period, groupName,
											varName = varName2))

  # check if the first network is one-mode, the second (potentially) bipartite
  if (dim(m1)[1] != dim(m1)[2]) {
	stop("Error: The first element in varName must be one-mode")
  }
  if (dim(m1)[1] != dim(m2)[1]) {
	stop("Error: Both elements in varName must have the same node sets")
  }
  # complement of a binary matrix
  cp <- function(m) (-m + 1)

  # all ties of two-paths in the one mode network
  onemode.reciprocal <- m1 * t(m1)
  onemode.forward <- m1 * cp(t(m1))
  onemode.backward <- cp(m1) * t(m1)
  onemode.null <- cp(m1) * cp(t(m1))
  diag(onemode.forward) <- 0
  diag(onemode.backward) <- 0
  diag(onemode.null) <- 0

  # one mode projections to the first mode
  bipartite.twopath <- m2 %*% t(m2)
  bipartite.null <- cp(m2) %*% cp(t(m2))
  bipartite.onestep1 <- m2 %*% cp(t(m2))
  bipartite.onestep2 <- cp(m2) %*% t(m2)
  diag(bipartite.twopath) <- 0
  diag(bipartite.null) <- 0
  diag(bipartite.onestep1) <- 0
  diag(bipartite.onestep2) <- 0

  # The coding is explained in the above referenced paper, pages 191-192
  # The first digit refers to the number of two-mode ties,
  # the second to the number of one-mode ties in triadic configurations
  # with two first-mode nodes, and one second-mode node.
  res <- c("22" =  sum(onemode.reciprocal * bipartite.twopath) / 2,
           "21" = (sum(onemode.forward * bipartite.twopath) +
						sum(onemode.backward * bipartite.twopath)) / 2,
           "20" =  sum(onemode.null * bipartite.twopath) / 2,
           "12" = (sum(onemode.reciprocal * bipartite.onestep1) +
						sum(onemode.reciprocal * bipartite.onestep2)) / 2,
           "11D" = (sum(onemode.forward * bipartite.onestep1) +
						sum(onemode.backward * bipartite.onestep2)) / 2,
           "11U" = (sum(onemode.forward * bipartite.onestep2) +
						sum(onemode.backward * bipartite.onestep1)) / 2,
           "10" = (sum(onemode.null * bipartite.onestep2) +
						sum(onemode.null * bipartite.onestep1)) / 2,
           "02" =  sum(onemode.reciprocal * bipartite.null) / 2,
           "01" = (sum(onemode.forward * bipartite.null) +
						sum(onemode.backward * bipartite.null)) / 2,
           "00" = sum(onemode.null * bipartite.null) / 2
  )

  # An ad-hoc check. Error should never be thrown.
  dim1 <- dim(m2)[1]
  dim2 <- dim(m2)[2]
  nTriads <- dim2 * dim1 * (dim1 - 1) / 2
  if (sum(res) != nTriads){
  stop(paste("Error in calculation. More than", nTriads,
						"triads counted (", sum(res), ")"))
  }
  res
}

##@TriadCensus sienaGOF Calculates triad census
# Contributed by Christoph Stadtfeld.
#
# Implementation of the Batagelj-Mrvar (Social Networks, 2001) algorithm
# based on the summary in the thesis of Sindhuja
#
TriadCensus <- function (i, obsData, sims, period, groupName, varName, levls = 1:16) {
  # get matrix and prepare data
  mat <- as.matrix(sparseMatrixExtraction(i, obsData, sims, period, groupName, varName = varName))
  N <- nrow(mat)
  # matrix with reciprocal ties
  matReciprocal <- mat + t(mat)
  # matrix with direction information for triad
  matDirected <- mat - t(mat)
  matDirected[matDirected == -1] <- 2
  matDirected[matReciprocal == 2] <- 3
  matDirected <- matDirected + 1
  # reciproal matrix with ties from higher to lower IDs
  matHigher <- matReciprocal
  matHigher[lower.tri(matHigher)] <- 0
  # neighbors lookup
  neighbors<- apply(matReciprocal, 1, function(x) which(x > 0))
  # neighbors with lower ids
  neighborsHigher <- apply(matHigher, 1, function(x) which(x > 0))

  # lookup table for 64 triad types
  # i->j, j->k, i->k
  # 1: empty, 2: forward, 3: backward, 4: reciprocal
  # order as in vector tc
  lookup <- array(NA, dim = rep(4,3))
  lookup[1,1,1] <- 1
  lookup[2,1,1] <- lookup[1,2,1] <- lookup[1,1,2] <- lookup[3,1,1] <- lookup[1,3,1] <- lookup[1,1,3] <- 2
  lookup[4,1,1] <- lookup[1,4,1] <- lookup[1,1,4] <- 3
  lookup[2,1,2] <- lookup[3,2,1] <- lookup[1,3,3] <- 4
  lookup[2,3,1] <- lookup[3,1,3] <- lookup[1,2,2] <- 5
  lookup[2,2,1] <- lookup[3,3,1] <- lookup[2,1,3] <- lookup[3,1,2] <- lookup[1,2,3] <- lookup[1,3,2] <- 6
  lookup[4,3,1] <- lookup[4,1,3] <- lookup[2,4,1] <- lookup[1,4,2] <- lookup[3,1,4] <- lookup[1,2,4] <- 7
  lookup[4,2,1] <- lookup[4,1,2] <- lookup[3,4,1] <- lookup[1,4,3] <- lookup[2,1,4] <- lookup[1,3,4] <- 8
  lookup[2,2,2] <- lookup[2,3,3] <- lookup[2,3,2] <- lookup[3,3,3] <- lookup[3,2,2] <- lookup[3,2,3] <- 9
  lookup[2,2,3] <- lookup[3,3,2] <- 10  # 3-cycle
  lookup[4,4,1] <- lookup[4,1,4] <- lookup[1,4,4] <- 11
  lookup[2,4,2] <- lookup[3,2,4] <- lookup[4,3,3] <- 12
  lookup[2,3,4] <- lookup[3,4,3] <- lookup[4,2,2] <- 13
  lookup[2,2,4] <- lookup[3,3,4] <- lookup[2,4,3] <- lookup[3,4,2] <- lookup[4,2,3] <- lookup[4,3,2]<- 14
  lookup[2,4,4] <- lookup[4,2,4] <- lookup[4,4,2] <- lookup[3,4,4] <- lookup[4,3,4] <- lookup[4,4,3] <- 15
  lookup[4,4,4] <- 16

  # initialize triad census
  tc <- c("003"  = 0,
          "012"  = 0,
          "102"  = 0,
          "021D" = 0,
          "021U" = 0,
          "021C" = 0,
          "111D" = 0,
          "111U" = 0,
          "030T" = 0,
          "030C" = 0,
          "201"  = 0,
          "120D" = 0,
          "120U" = 0,
          "120C" = 0,
          "210"  = 0,
          "300"  = 0)

  # iterate through all non-empty dyads (from lower to higher ID)
  if (length(neighborsHigher) > 0){ # else mat is the zero matrix
	for(ii in 1:N){
		for(j in neighborsHigher[[ii]]){
      # set of nodes that are linked to ii and j
		third <- setdiff( union(neighbors[[ii]], neighbors[[j]]),
                    c(ii, j) )
		# store triads with just one tie
		triadType <- ifelse(matReciprocal[ii,j] == 2, 3, 2)
		tc[triadType] <- tc[triadType] + N - length(third) - 2
		for (k in third){
        # only store triads once
			if(j < k || ( ii < k && k < j && !(k %in% neighbors[[ii]]) ) ){
				t1 <- matDirected[ii,j]
				t2 <- matDirected[j,k]
				t3 <- matDirected[ii,k]
				triadType <- lookup[t1, t2, t3]
				tc[triadType] <- tc[triadType] + 1
				}
			}
		}
	}
  }
  # assign residual to empty triad count
  tc[1] <- 1/6 * N*(N-1)*(N-2) - sum(tc[2:16])
  tc[levls]
}


##@dyadicCov sienaGOF Auxiliary variable for dyadic covariate
#
# An auxiliary function calculating the proportion of ties
# for subsets of ordered pairs corresponding to
# certain values of the categorical dyadic covariate dc.
# dc should be a matrix of the same dimensions as
# the dependent variable,
# or an array where the third dimension is time.
# Frequencies of ties with dc == 0 are not counted.
dyadicCov <-  function (i, obsData, sims, period, groupName, varName, dc){
	m <- sparseMatrixExtraction(i, obsData, sims, period, groupName, varName)
    # note that m*dc is a sparse matrix, too:
	if (length(dim(dc))==3)
	{
		tmdyv <- table((m * dc[,,period])@x, useNA = "no")
	}
	else
	{
		tmdyv <- table((m * dc)@x, useNA = "no")
	}
	values <- unique(as.vector(dc))
	tdyv <- sort(values[!is.na(values)])
	tdyv <- tdyv[-which(tdyv==0)] # if 0 is included, take it out
	# Now we want to construct the table of numbers of m*dyv;
	# and categories in dc not represented in m*dyv should get a 0.
	# First make a named vector of the correct length with 0s in place.
	ttmdyv <- 0*tdyv
	names(ttmdyv) <- tdyv
	dims <- dimnames(tmdyv)[[1]]
	ttmdyv[dims] <- tmdyv # The other entries remain 0
	ttmdyv
}
