##*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: https://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: sienaTimeTest.r
## *
## * Description: This file contains the sienaTimeTest code for testing the
## * significance of additional time dummies interacted with effects, and
## * sienaTimeFix, which is called to set up time dummy interacted effects.
## ****************************************************************************/
##@sienaTimeTest siena07 Does test for time homogeneity of effects
sienaTimeTest <- function (sienaFit, effects=NULL, excludedEffects=NULL,
	condition=FALSE)
{
	if (!inherits(sienaFit, "sienaFit"))
	{
		stop("not a legitimate sienaFit object")
	}
	if (!sienaFit$OK)
	{
		stop("Error end of estimation algorithm")
	}
	else if (sienaFit$termination == "UserInterrupt")
	{
		stop("User interrupted run, object possibly incomplete")
	}
	if (!is.null(sienaFit$gmm))
	{
		if (sienaFit$gmm)
		{
			stop("sienaTimeTest cannot be applied to results of GMoM estimation")
		}
	}
	if (!is.null(sienaFit$sf2.byIteration))
	{
		if (!sienaFit$sf2.byIteration)
		{
			stop("sienaTimeTest needs sf2 by iterations (use lessMem=FALSE)")
		}
	}
	waveNumbers <- attr(sienaFit$f, "periodNos")
	nWaves <- length(waveNumbers)

	fitEffects <- sienaFit$requestedEffects

	## There must be more than 2 periods (more than 1 wave)
	## to do a time test!
	if (nWaves < 2)
	{
		stop("You must have at least two time periods (three waves) to test ",
			"for non-heterogeneity across time.")
	}
	## if this is a maximum likelihood fit or finite differences we need
	## the optional objects sdf2
	if (sienaFit$maxlike || sienaFit$FinDiff.method)
	{
		if (is.null(sienaFit$sdf2[[1]][[1]]))
		{
			stop("rerun siena07 with the byWave option TRUE")
		}
	}
	## get the desired effects
	if (!is.null(effects))
	{
		## a little validation helps ensure we have at least one effect we want
		if (any(!is.numeric(effects)))
		{
			stop("non numeric effect number requested")
		}
		if (any(!effects %in% 1:nrow(fitEffects)))
		{
			stop("Some effect number requested not found")
		}
		use <- 1:nrow(fitEffects) %in% effects
		if (any(fitEffects$basicRate[use]))
		{
			stop("siena time tests are inappropriate for basic rates")
		}
		#        if (any(grepl("Dummy", fitEffects$effectName[use])))
		#        {
		#            stop("siena time tests are inappropriate for time dummy effects")
		#        }
	}
	else
	{
		use <- !fitEffects$basicRate
	}

	## if (sienaFit$maxlike || sienaFit$FinDiff.method)
	## {
	##     stop("Not yet implemented for finite differences or maxlike")
	## }
	## Identify the effects which will potentially be tested
	baseInFit <- use & !grepl("Dummy", fitEffects$effectName)

	nBaseEffects <- sum(baseInFit)
	if (nBaseEffects == 0)
	{
		stop("No effects available to test")
	}

	fixedDummies <- fitEffects$shortName=='egoX' &
		fitEffects$fix & grepl("Dummy", fitEffects$effectName)

	## establish effects for top left of derivative matrix D
	estimatedInFit <- use &  !fixedDummies
	topleftEffectNumbers <- fitEffects$effectNumber[estimatedInFit]

	nRowsToTest <- nBaseEffects * nWaves

	## construct the data frame that will be the basis for all calculations
	toTest <- data.frame(baseEffect=rep(1:nBaseEffects, each=nWaves),
		effectNumber=rep(fitEffects$effectNumber[baseInFit],
			each=nWaves),
		type=rep(fitEffects$type[baseInFit], each=nWaves),
		period = rep(waveNumbers, nBaseEffects),
		period1 = rep(1:nWaves, nBaseEffects),
		toTest=rep(TRUE, nRowsToTest),
		baseRowInD=rep(0, nRowsToTest),
		effectName=rep(fitEffects$effectName[baseInFit],
			each=nWaves),
		rowInD=rep(0, nRowsToTest))
	toTest$toTest[toTest$period == 1] <- FALSE
	toTest$rowInD <- match(toTest$effectNumber,
		topleftEffectNumbers, nomatch=0)
	toTest$baseRowInD <- toTest$rowInD
	toTest$rowInD[toTest$period > 1]  <- 0
	## Go through each effect which is an estimated time dummy, and
	## incorporate this information into the toTest data frame. i.e.
	## if a time dummy was estimated, set
	## its element in toTest equal to FALSE so that we do not time test it
	for (i in which(grepl("isDummy", fitEffects$timeDummy) & estimatedInFit))
	{
		tmp <- toString(fitEffects$timeDummy[i])
		tmp <- strsplit(tmp, split=",", fixed=TRUE)[[1]]
		if (any(tmp != ""))
		{
			## Dont test the dummy for the corresponding effect
			thisRow <- toTest$effectNumber == tmp[3] & toTest$period == tmp[2]
			toTest[thisRow, "toTest"] <- FALSE
			## We want to be able to reference this effect given an
			## index for the base effect and a time period, so store
			## this information in rowInD -- this is used
			## extensively in plot.sienaTimeTest
			toTest[thisRow, "rowInD"] <-
				match(fitEffects$effectNumber[i],
					topleftEffectNumbers)
		}
	}	
	##  nEffects, nSims, nameslist, nDummies convert commonly used ingredients
	##  from sienaFit into an easily accessed form based on the screens
	##  set up above
	nEffects <- sum(!toTest$toTest)
	## add row indices in D to the rows to be tested
	toTest$rowInD[toTest$toTest] <- nEffects + 1:sum(toTest$toTest)
	toTest <- toTest[order(toTest$rowInD), ]

	nSims <- sienaFit$Phase3nits
	nameslist <- list(Iteration=paste("it", 1:nSims, sep=""),
		Wave=paste("Wave", waveNumbers, sep=""),
		Effect=fitEffects$effectName[estimatedInFit]
		)
	nDummies <- sum(toTest$toTest)
	ttt <- toTest$toTest
	type <- ifelse(toTest$type == "eval", "",
		paste(" (", toTest$type, ")", sep=""))
	toTest$dummyNames[!ttt] <- paste(fitEffects$effectName[estimatedInFit],
		type[!ttt], sep="")
	toTest$dummyNames[ttt] <- paste("(*)Dummy", toTest$period[ttt], ":",
		toTest$effectName[ttt], type[ttt], sep="")

	## obsStats, moment, scores are the crucial ingredients from sienaFit which
	## screen for the base effects and make the rest of the code clean
	## be careful as the original logical vectors cannot be used once rows
	## have been extracted.
	obsStats <- t(sienaFit$targets2[estimatedInFit, , drop=FALSE])
	moment <- sienaFit$sf2[, , estimatedInFit, drop=FALSE] -
		rep(obsStats, each=nSims)
	G <- array(0, dim=c(nSims, nWaves, nEffects + nDummies))
	## Set the base effects G equal to the moments from sienaFit
	G[, , 1:nEffects] <- moment
	## copy over the others
	subs1 <- cbind(rep(1:nSims, nDummies), rep(toTest$period1[ttt], each=nSims),
		rep(toTest$rowInD[ttt], each=nSims))
	subs2 <- cbind(rep(1:nSims, nDummies), rep(toTest$period1[ttt], each=nSims),
		rep(toTest$baseRowInD[ttt], each=nSims))
	G[subs1] <- moment[subs2]
	## Put names onto G for easy? reference
	## Use dimnames(G) <- list(NULL, NULL) to view the data as a matrix!
	dimnames(G) <- list(nameslist$Iteration, nameslist$Wave, toTest$dummyNames)
	## Make the covariance matrix for the new moments
	sigma <- cov(apply(G, c(1, 3), sum))
	if (length(sigma) == 1)
	{
		sigma <- matrix(sigma,1,1)
	}
	if (!(sienaFit$maxlike || sienaFit$FinDiff.method))
	{
		scores <- sienaFit$ssc[ , , estimatedInFit, drop=FALSE]
		SF <- array(0, dim=c(nSims, nWaves, nEffects + nDummies))
		SF[, , 1:nEffects] <- scores
		SF[subs1] <- scores[subs2]
		## Put names onto SF for easy reference
		dimnames(SF) <- dimnames(G)
		D <- derivativeFromScoresAndDeviations(SF, G, , , , TRUE, )
	}
	else
	{
		if (sienaFit$FinDiff.method)
		{
			##sienaFit$sdf2 is a list of arrays, one for each iteration
			derivs <- t(sapply(sienaFit$sdf2, function(x)
					x[, estimatedInFit, estimatedInFit, drop=FALSE]))
		}
		else
		{
			##sienaFit$sdf2 is a list of lists of Matrix objects
			derivs <-
				t(sapply(sienaFit$sdf2, function(xx)
						array(t(sapply(xx, function(x)
									as.matrix(x[estimatedInFit,
										estimatedInFit,
										drop=FALSE]))),
						dim=c(nWaves, nEffects, nEffects))))
		}
		derivs <- array(derivs, dim=c(nSims, nWaves, nEffects, nEffects))
		DF <- array(0, dim=c(nSims, nWaves, nEffects + nDummies,
				nEffects + nDummies))
		DF[, , 1:nEffects, 1:nEffects] <- derivs
		for (wave in 2:nWaves)
		{
			thisWave <- toTest$period == wave & toTest$toTest
			subs1 <- (1: (nEffects + nDummies))[thisWave]
			subs2 <- toTest$baseRowInD[thisWave]
			DF[, wave, subs1, subs1] <- derivs[, wave, subs2, subs2]

			DF[, wave, subs1, 1:nEffects] <- derivs[, wave, subs2, 1:nEffects]
			DF[, wave, 1:nEffects, subs1] <- derivs[, wave, 1:nEffects, subs2]
		}

		D <- t(apply(DF, c(3, 4), mean))
	}

	## The following parts about excluding effects were added later
	## (preceding code is mainly from JL; parts about excluding are from TS).
	## This leads to some duplication in functionality;
	## but since the preceding construction works well,
	## I (TS) wanted to keep it mostly as it was/is.

	## Check for fixed effects and exclude them from testing.
	if (is.null(excludedEffects))
	{
		excludedEffects <- which(fitEffects$fix)
	}
	else
	{
		excludedEffects <-
			unique(c(excludedEffects, which(fitEffects$fix)))
	}

	## In case some effects are excluded from testing, this is done now:
	if ((!is.null(excludedEffects)) & (length(excludedEffects) > 0))
	{
		useRows <- rep(TRUE, dim(toTest)[1])
		for (i in seq_along(excludedEffects))
		{
			useRows[toTest$baseRowInD==excludedEffects[i]] <- FALSE
		}
		useRows[1:nBaseEffects] <- TRUE
		toTest <- toTest[useRows,]
		D <- D[useRows, useRows]
		sigma <- sigma[useRows, useRows, drop=FALSE]
		nDummies <- sum(toTest$toTest)
		G <- G[,,useRows]
		nRowsToTest <- dim(toTest)[1]
		exclusions <- 1:nBaseEffects %in% excludedEffects
	}
	else
	{
		exclusions <- rep(FALSE, nBaseEffects)
	}

	## Check if there are any collinearities in the effects and dummies;
	## note that time dummies added by sienaTimeFix were already excluded.
	rankSigma <- qr(sigma)$rank
	extraExclusions <- rep(FALSE, nBaseEffects)
	if (rankSigma < dim(sigma)[1])
	{
		es <- ifelse((nRowsToTest - nDummies == 1), '', 's')
		message("TimeTest constructed a null hypothesis with ",
							nRowsToTest - nDummies, " estimated parameter", es)
		message("and ", nDummies," dummy variables to be tested.")
		message("However, there are ", dim(sigma)[1] - rankSigma,
							" linear dependencies between these.\n")
		message("This may be because some of the parameters are already")
		message("interactions with time dummies or other time variables.\n")
		message("Look in the help page for the use of the <effects> parameter.\n")

		subset0 <- rep(FALSE, dim(sigma)[1])
		subset0[1:nBaseEffects] <- TRUE
		for (i in 1:nBaseEffects)
		{
			if (!(i %in% excludedEffects))
				## Check if the dummies for this effect
				## are dependent on the base effects.
			{
				subset <- subset0 | (toTest$baseRowInD==i)
				if (qr(sigma[subset,subset])$rank
					< dim(sigma[subset,subset])[1])
				{
					extraExclusions[i] <- TRUE
				}
			}
		}
		if (any(extraExclusions))
		{
			useRows <- rep(TRUE, dim(toTest)[1])
			for (i in 1:nBaseEffects)
			{
				if (extraExclusions[i])
				{
					useRows[toTest$baseRowInD==i] <- FALSE
				}
				useRows[1:nBaseEffects] <- TRUE
			}
			toTest <- toTest[useRows,, drop=FALSE]
			D <- D[useRows, useRows, drop=FALSE]
			sigma <- sigma[useRows, useRows, drop=FALSE]
			nDummies <- sum(toTest$toTest)
			G <- G[,,useRows, drop=FALSE]
			nRowsToTest <- dim(toTest)[1]
			excludedNumber <- sum(extraExclusions)
			es <- ifelse((excludedNumber == 1), '', 's')
			efNames <- fitEffects[baseInFit,"effectName"]
	cat("Automatic discovery of dependencies yielded the exclusion of effect", 
							es, sep="")
			if (excludedNumber <= 1)
			{			
				cat(" ", efNames[extraExclusions],".\n ")
			}
			else
			{
				cat(" ", paste(efNames[extraExclusions], collapse=", "),".\n ")
			}
			rankSigma <- qr(sigma)$rank
			if (rankSigma < dim(sigma)[1])
			{
				cat("After these exclusions, there still are linear dependencies.\n")
				cat("Advice: use sienaTimeTest with a smaller set of effects.\n")
				cat("Look in the help page for the use of the <effects> parameter.\n\n")
				stop("Please rerun sienaTimeTest with appropriate excluded effects.")
			}
		}
		else
		{
			cat("Automatic discovery of dependencies had no effect.\n")
			cat("Advice: use sienaTimeTest with a smaller set of effects.\n")
			cat("Look in the help page for the use of the <effects> parameter.\n\n")
			stop("Please rerun sienaTimeTest with appropriate excluded effects.")
		}
	}
	## We have now set up all of the ingredients properly, so we may proceed
	## with the score type test of Schweinberger (2007-2012)
	fra <- apply(G, 3, sum) / nSims
	doTests <- toTest$toTest
	redundant <- rep(FALSE, length(doTests))
	jointTest <- ScoreTest(nrow(toTest), D, sigma, fra, doTests, redundant,
		maxlike=sienaFit$maxlike)
	jointTestP <- 1 - pchisq(jointTest$testresOverall, nDummies)
	if (! condition)
	{
		individualTest <- jointTest$testresulto[1:nDummies]
	}
	else
	{
		individualTest <- sapply(1:nDummies, function (i)
			{
				doTests <- rep(FALSE, nEffects + nDummies)
				doTests[nEffects + i] <- TRUE
				test <- ScoreTest(nrow(toTest), D, sigma,
					fra, doTests, redundant, maxlike=sienaFit$maxlike)
				test$testresulto[1]
			}
			)
	}
	individualTestP <- 2 * (1- pnorm(abs(individualTest)))
	if (!is.na(jointTestP))
	{
		rownames(jointTestP) <- "Joint Significant Test"
		colnames(jointTestP) <- "p-Val"
	}
	thetaOneStep <- c(sienaFit$theta[estimatedInFit, drop=FALSE], rep(0, nDummies)) +
		jointTest$oneStep
	## Define the null hypothesis for the later tests (effect-wise and group-wise)
	## as the null hypothesis tested in the overall test:
	nullHyp <- 	!toTest$toTest
	effectTest <- as.vector(by(toTest, toTest$baseEffect, function (x)
				{
					doTests <- rep(FALSE, nEffects + nDummies)
					if (any(x$toTest))
					{
						if (!condition)
						{
							## Control for all estimated effects.
							doTests[toTest$baseEffect == (x$baseEffect[1]) &
								toTest$toTest] <- TRUE
							testresult <- partial.scoreTest(D, sigma, fra,
								doTests, nullHyp, maxlike=sienaFit$maxlike)
							testresult
						}
						else
						{
							## Control for everything.
							doTests[toTest$baseEffect == (x$baseEffect[1]) &
								toTest$toTest] <- TRUE
							redundant <- rep(FALSE, length(doTests))
							test <- ScoreTest(nEffects + nDummies, D, sigma,
								fra, doTests, redundant, maxlike=sienaFit$maxlike)
							test$testresOverall
						}
					}
					else
					{
						NA
					}
				}
				))
	dim(effectTest) <- c(length(effectTest), 1)
	effectTestP <- matrix(NA, dim(effectTest)[1], 3)
	effectTestP[,1] <- round(effectTest[,1],2)
	effectTestP[,2] <- tapply(toTest$toTest, toTest$baseEffect, sum)
	effectTestP[,3] <- round(1 - pchisq(effectTestP[,1], effectTestP[,2]), 3)
	rownames(effectTestP) <- toTest$dummyNames[toTest$period == 1]
	colnames(effectTestP) <- c("chi-sq.", "df", "p-value")
	pvalues <-
		round(c(2 * (1 -
					pnorm(abs(sienaFit$theta[estimatedInFit] /
							sqrt(diag(sienaFit$covtheta)[estimatedInFit])))),
				individualTestP), 3)
	thetaStar <- cbind(c(sienaFit$theta[estimatedInFit], rep(0, nDummies)),
		thetaOneStep, pvalues)
	colnames(thetaStar) <- c("Initial Est.", "One Step Est.", "p-Value")
	rownames(thetaStar) <- dimnames(G)[[3]]
	## put things on toTest to make plot easier
	toTest[, c("InitialEst", "OneStepEst", "p.value")] <- thetaStar
	toTest$effectTest <- NA
	toTest$effectTest <- effectTestP[toTest$baseEffect, 1]
	type <- ifelse(toTest$type =="eval", "", paste(" (",
			as.character(toTest$type), ")", sep=""))
	toTest$effectName <-
		factor(paste(toTest$effectName, type,
				" \n(chisq=", toTest$effectTest, ")", sep=""))
	toTest$valsplus <- toTest$OneStepEst +
		ifelse(toTest$period == 1, 0, toTest$OneStepEst[toTest$baseEffect])
	toTest$dummysd <- abs(toTest$OneStepEst / qnorm(1 - toTest$p.value / 2))
	toTest$dummysd[toTest$period == 1] <-
		sqrt(diag(sienaFit$covtheta))[estimatedInFit][toTest$period == 1]

	# Added by Tom: groupwise tests.
	cef <- attr(sienaFit$f,"condEffects")
	if (inherits(cef,"sienaGroupEffects"))
	{
		# This is a multi-group result
		groupsKind <- "Group"
		groupSizes <- attr(sienaFit$f,"groupPeriods") - 1
		nGroups <- sienaFit$f$nGroup
		if (nGroups != length(groupSizes)) {stop("No good group sizes 1.")}
		groups <- rep(1:nGroups, groupSizes)
		if (nWaves != length(groups)) {stop("No good group sizes 2.")}
		toTest$whichGroup <- groups[toTest$period1]
	}
	else
	{
		groupsKind <- "Period"
		nGroups <- sienaFit$observations-1
		toTest$whichGroup <- toTest$period
	}
	groupTest <-
		as.vector(by(toTest, toTest$whichGroup, function (x)
				{
					doTests <- rep(FALSE, nEffects + nDummies)
					if (any(x$toTest))
					{
						doTests[toTest$whichGroup == (x$whichGroup[1]) &
							toTest$toTest] <- TRUE
						testresult <- partial.scoreTest(D, sigma, fra, doTests,
							nullHyp, maxlike=sienaFit$maxlike)
						testresult
					}
					else
					{
						NA
					}
				}
				))
	# The preceding does not yield the test for the first group/period,
	# because it was used as the reference category.
	# This is a test for the joint contribution 
	# of the dummies for the other groups.
	# Therefore matrix A is created as the design matrix for this test.	
	remainingEffects <- unique(toTest$baseEffect[toTest$toTest])
	nRemainingEffects <- length(remainingEffects)
	if (nRemainingEffects >= 1)
	{
		A <- matrix(0, nRemainingEffects, nRowsToTest)
		for (i in seq_along(remainingEffects))
		{
			A[i, which((toTest$baseEffect==remainingEffects[i]) & (!nullHyp))] <- 1
		}		
		groupTest[1] <- transformed.scoreTest(D, sigma, fra, A,
			nullHyp, maxlike=sienaFit$maxlike)
		dim(groupTest) <- c(nGroups, 1)
		groupTestP <- matrix(NA, nGroups, 3)
		groupTestP[,1] <- round(groupTest[,1], 2)
		groupTestP[,2] <- tapply(toTest$toTest, toTest$whichGroup, sum)
		groupTestP[1,2] <- nRemainingEffects
		groupTestP[,3] <- round(1 - pchisq(groupTest[,1], groupTestP[,2]), 3)
		rownames(groupTestP) <- paste(groupsKind, 1:nGroups)
		colnames(groupTestP) <- c("chi-sq.", "df", "p-value")
	}
	else
	{
		groupTest[1] <- NA
		groupTestP <- NA
		cat('No first dummy tested ...\n')
	}
	returnObj <- list(
		sienaFitName=deparse(substitute(sienaFit)),
		JointTest=jointTestP,
		EffectTest=effectTestP,
		IndividualTest=thetaStar,
		JointTestStatistics=jointTest,
		EffectTestStatistics=effectTest,
		IndividualTestStatistics=individualTest,
		GroupTest=groupTestP,
		GroupTestStatistics=groupTest,
		GroupsKind=groupsKind,
		CovDummyEst=jointTest$covMatrix,
		Moments=G,
		Deriv=D,
		BaseRowInD=match(which(baseInFit),
			which(estimatedInFit)),
		Waves=dim(G)[2],
		Sims=dim(G)[1],
		Effects=dim(G)[3],
		DummyStdErr=sqrt(diag(jointTest$covMatrix)),
		OriginalEffects=nEffects,
		OriginalThetaStderr=
			sqrt(diag(sienaFit$covtheta))[estimatedInFit],
		ToTest=toTest,
		ScreenedEffects=which(!use),
		RequestedExclusions=exclusions,
		ExtraExclusions=extraExclusions,
		WaveNumbers=waveNumbers,
		IndividualTestsOrthogonalized=condition
		)
	class(returnObj) <- "sienaTimeTest"
	attr(returnObj, "version") <- packageDescription(pkgname, fields = "Version")
	returnObj
}
##@transformed.scoreTest siena07 score test for part of the additional effects
transformed.scoreTest <- function(dfra, msf, fra, A, nullHyp, maxlike)
	# score test of A %*% fra;
	# Boolean vector nullHyp defines vector of estimated effects
	# non-zero columns of A should not overlap with nullHyp
{
	test <- colSums(abs(A)) > 0
	if (any(test & nullHyp))
	{
		stop("Null hyp and test in transformed score test should not overlap")
	}
	if (dim(A)[1] + sum(nullHyp) < 2)
	{
		stop("Meaningless transformed score test")
	}
	A1 <- matrix(0, sum(nullHyp), dim(dfra)[1])
	A1[cbind(1:sum(nullHyp), which(nullHyp))] <- 1
	AA <- rbind(A1,A)
	D <- AA %*% dfra %*% t(AA)
	sig <- AA %*% msf %*% t(AA)
	g <- AA %*% fra
	tested <- rep(FALSE, dim(AA)[1])
	tested[(sum(nullHyp)+1):dim(AA)[1]] <- TRUE
	redundant <- rep(FALSE, dim(AA)[1])
	EvaluateTestStatistic(maxlike, tested, redundant, D, sig, g)$cvalue
}
##@partial.scoreTest siena07 score test for part of the additional effects
partial.scoreTest <- function(dfra, msf, fra, test, nullHyp, maxlike)
	# test gives vector of tested effects,
	# nullHyp gives vector of estimated effects
{
	if (any(test & nullHyp))
	{
		stop("Null hypothesis and test in partial score test should not overlap")
	}
	if (!any(test)){stop("Nothing to be tested in partial score test")}
	A <- matrix(0, sum(test), dim(dfra)[1])
	A[cbind(1:sum(test),which(test))] <- 1
	transformed.scoreTest(dfra, msf, fra, A, nullHyp, maxlike)
}
##@summary.sienaTimeTest siena07 summary method for sienaTimeTest objects
summary.sienaTimeTest <- function(object, ...)
{
	if (!inherits(object, "sienaTimeTest"))
	{
		stop("not a legitimate Siena time test object")
	}
	class(object) <- c("summary.sienaTimeTest", class(object))
	object
}
##@print.summary.sienaTimeTest siena07 print method for summary.sienaTimeTest
print.summary.sienaTimeTest <- function(x, ...)
{
	if (!inherits(x, "summary.sienaTimeTest"))
	{
		stop("not a legitimate Siena time test summary object")
	}
	print.sienaTimeTest(x)
	## Additional output to the print will go in here:
	## The first condition is merely for compatibility with older versions.
	if ((!is.null(x$RequestedExclusions))||(!is.null(x$ExtraExclusions)))
	{
		if (any(x$RequestedExclusions) || any(x$ExtraExclusions))
		{
			cat("\nSome effects were excluded from being tested.\n")
		}

		if (any(x$RequestedExclusions))
		{
			if (sum(x$RequestedExclusions) <= 1)
			{
				cat("The effect with number", which(x$RequestedExclusions))
				cat(" in", x$sienaFitName, "was excluded.\n")
			}
			else
			{
				cat("\nRequested exclusions were the effects numbered\n")
				cat(which(x$RequestedExclusions)[1],
					paste(", ",
						which(x$RequestedExclusions)[2:sum(x$RequestedExclusions)],
						sep=""),
					" in ", x$sienaFitName, ".\n", sep="")
			}
		}

		if (any(x$ExtraExclusions))
		{
			if (sum(x$ExtraExclusions) <= 1)
			{
				cat("The effect with number", which(x$ExtraExclusions),
					" in", x$sienaFitName,
					" was excluded automatically because of collinearity.\n")
			}
			else
			{
				cat("\nThe effects ")
				cat(which(x$ExtraExclusions)[1], paste(", ",
						which(x$ExtraExclusions)[2:sum(x$ExtraExclusions)],
						sep=""),
					" in ", x$sienaFitName, ".\n", sep="")
				cat("were excluded automatically because of collinearity.\n")
			}
		}
	}
	cat("\nIndividual significance tests and one-step estimators:\n")
	print(format(round(x$IndividualTest, digits=4)), quote=F)
	cat("\nEffect-wise joint significance tests\n")
	cat("(i.e. each effect across all dummies):\n")
	print(x$EffectTest)
	if (x$Waves <=2 && ! x$IndividualTestsOrthogonalized)
	{
		#		cat("\n\nNote that these effect-wise tests have a different form")
		#		cat("\nthan the individual tests, thus testing with 3 observations")
		#		cat("\nmay yield different individual and effect-wise values.\n\n")
	}
	else
	{
		if (x$IndividualTestsOrthogonalized)
		{
			cat("\nNote that the 1-df and effect-wise test statistics were",
				"\northogonalized with respect to each other (condition=TRUE).")
			cat("\n")
		}
	}
	if (!is.null(x$GroupTest))
	{
		cat(paste("\n",x$GroupsKind,"-wise joint significance tests\n", sep=""))
		cat(paste("(i.e. each", tolower(x$GroupsKind),
				"across all parameters):\n"))
		print(x$GroupTest)
	}

	tmp <- paste(" (", 1:length(x$BaseRowInD), ") ",
		rownames(x$IndividualTest)[x$BaseRowInD], "\n", sep="")
	cat("\nUse the following indices for plotting:\n",
		tmp[!is.na(x$EffectTestStatistics)])
	cat("\nIf you would like to fit time dummies to your model,")
	cat("\nuse the includeTimeDummy function.")
	cat("\nType \"?sienaTimeTest\" for more information on this output.\n")
	invisible(x)
}
##@print.sienaTimeTest siena07 print method for sienaTimeTest objects
print.sienaTimeTest <- function(x, ...)
{
	if (!inherits(x, "sienaTimeTest"))
	{
		stop("not a legitimate Siena time test object")
	}
	effectNames <- rownames(x$IndividualTest)
	dummies <- x$ToTest$toTest
	dummyIndex <- paste(" (", format(1:sum(dummies)), ") ",
		effectNames[dummies], "\n", sep="")
	pvalue <- ifelse((x$JointTest < 0.0001)," < 0.0001",
		paste("=",sprintf("%6.4f", x$JointTest)))
	cat("Joint significance test of time heterogeneity:\n",
		"chi-squared = ",
		sprintf("%6.2f", x$JointTestStatistics$testresOverall),
		", d.f. = ", sprintf("%1d", length(x$IndividualTestStatistics)),
		", p", pvalue,
		", \nwhere H0: The following parameters are zero:\n",
		dummyIndex, sep="")
	invisible(x)
}
##@plot.sienaTimeTest siena07 plot method for sienaTimeTest objects
plot.sienaTimeTest <- function(x, pairwise=FALSE, effects,
	scale=0.2, plevels=c(0.1, 0.05, 0.025), ...)
{
	## require(lattice)
	tmp <- paste(" (", 1:length(x$BaseRowInD), ") ",
		rownames(x$IndividualTest)[x$BaseRowInD], "\n", sep="")
	if (missing(effects))
	{
		effects<- (1:length(tmp))[!is.na(x$EffectTestStatistics)]
	}
	if (any(!effects %in% 1:length(tmp)))
	{
		cat("Detected an error with the effects included. For a
			parameter-plot, use the following indices:")
			stop("\nUse the following indices for plotting the pairwise
				moment correlations:\n", tmp)
	}
	if (pairwise)
	{
		## On a "pairwise" call, print a pairwise plot of moments
		if (length(effects) == 0)
		{
			x <- x$Moments
		}
		else
		{
			x <- x$Moments[, , effects, drop=FALSE]
		}
		panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
		{
#			usr <- par("usr"); # on.exit(par("usr" = usr))
			par(usr = c(0, 1, 0, 1))
			r <- abs(cor(x, y))
			txt <- format(c(r, 0.123456789), digits=digits)[1]
			txt <- paste(prefix, txt, sep="")
			if(missing(cex.cor)) cex.cor <- 0.8 / strwidth(txt)
			text(0.5, 0.5, txt, cex = cex.cor * r)
		}
		pairs(apply(x, c(1, 3), sum),
			lower.panel=panel.smooth,
			upper.panel=panel.cor,
			pch=5, ...) ##pch=5 is a diamond
	}
	else
	{
		# require(grid)
		## Otherwise, make the parameter plots:
		## create a vector to split the IndividualTests by
		toTest <- x$ToTest
		waveNumbers <- x$WaveNumbers
		nWaves <- x$Waves
		nLevels <- length(plevels)
		cis <- lapply(1:nLevels, function(i)
			{
				cbind(toTest$valsplus - abs(qnorm(plevels[i]  /  2,
							sd=toTest$dummysd)),
					toTest$valsplus  + abs(qnorm(plevels[i] / 2,
							sd=toTest$dummysd)))
			}
			)
		toTest[, paste(rep(c("lower", "upper"), nLevels),
			rep(1:nLevels, each=2), sep="")] <-
				do.call(cbind, cis)
			## subset the data frame
			toTest <- toTest[toTest$baseEffect %in% effects, ]
			toTest$effectName <- factor(toTest$effectName,
				levels=unique(toTest$effectName))
			if (nlevels(toTest$effectName) != length(unique(toTest$baseEffect)))
			{
				stop("non unique effect names: contact RSiena help")
			}
			#devAskNewPage(TRUE)
			xyplot(valsplus ~ period1 | effectName, data=toTest,
				type = "p", as.table=TRUE, bty="n",
				xlab="Wave", ylab="Parameter Value", auto.key=TRUE,
				scales=list(relation='free', x=list(at=0:(nWaves+1),
						labels=c(" ", waveNumbers, " "))), subscripts=TRUE,
				par.strip.text=list(lines=2.2), #layout=layout,
				xlim=c(0, nWaves + 1), prepanel=function(x, y, ...)
				{
					ymin <- min(y - scale * abs(y))
					ymax <- max(y + scale * abs(y))
					list(ylim=c(ymin, ymax))
				},
				panel=function(x, y, subscripts, ...)
				{
					sapply(1:length(x), function(i, subscripts)
						{ ## pch 45 is a hyphen, pch 20 a filled circle
							tmp <- if(toTest$toTest[subscripts][i]) "red" else "gray"
							panel.xyplot(c(x[i], x[i]),
								c(toTest$lower1[subscripts][i],
									toTest$upper1[subscripts][i]),
								col=tmp, alpha=.50,
								type="l", lend=1, lwd=10)
							panel.xyplot(c(x[i], x[i]),
								c(toTest$lower1[subscripts][i],
									toTest$upper1[subscripts][i]),
								col=tmp, alpha=.75,
								type="p", pch=45, cex=3)
							panel.xyplot(c(x[i], x[i]),
								c(toTest$lower2[subscripts][i],
									toTest$upper2[subscripts][i]),
								col=tmp, alpha=.50,
								type="l", lend=1, lwd=10)
							panel.xyplot(c(x[i], x[i]),
								c(toTest$lower2[subscripts][i],
									toTest$upper2[subscripts][i]),
								col=tmp, alpha=.75,
								type="p", pch=45, cex=3)
							panel.xyplot(c(x[i], x[i]),
								c(toTest$lower3[subscripts][i],
									toTest$upper3[subscripts][i]),
								col=tmp, alpha=.25,
								type="l", lend=1, lwd=10)
						}, subscripts=subscripts)
					panel.xyplot(x, y, type="s",
						col="black", alpha=.75, pch=2)
					panel.xyplot(x, y, type="p",  pch=20,
						col=1)
					panel.abline(a=toTest[subscripts[1], "initial"],
						reference=TRUE,
						col="black", lwd=2, alpha=.75)
					# require grid	
					# grid.text(paste("p = ", toTest$effectTest[subscripts[1]]),
					#           0.02, 0.1, just="left")
				}, ...)
	}
}

##@sienaTimeFix siena07 Adds time dummy terms to the effects object
sienaTimeFix <- function(effects, data=NULL, getDocumentation=FALSE)
{
	##@addEffect internal sienaTimeFix add one or more effects
	addEffect <- function(newEffects, i, newname, effectGroup, shortName,
		timeDummy, fix=FALSE, include=TRUE, yName=NULL)
	{
		tmprows <- createEffects(effectGroup, xName=newname, yName=yName,
			name=effects$name[i][1],
			groupName=effects$groupName[i][1],
			group=effects$group[i][1],
			netType=effects$netType[i][1])
		tmprows <- tmprows[tmprows$shortName==shortName &
			tmprows$type %in% effects$type[i], ]
		tmprows$fix <- fix
		tmprows$include <- include
		tmprows$effectNumber <- max(newEffects$effectNumber) +
			(1:nrow(tmprows))
		tmprows$timeDummy <- timeDummy
		rownames(tmprows) <- paste(newname, effects$type[i], sep=".")
		newEffects <- rbind(newEffects, tmprows)
		newEffects
	}

	##@addVarCovar internal sienaTimeFix add a varying covariate
	addVarCovar <- function(data, base, nodeSet, newname, dataSub)
	{
		base <- varCovar(base, nodeSet=nodeSet, warn=FALSE)
		base <- addAttributes.varCovar(base, name=dname)
		data[[dataSub]]$vCovars <- c(data[[dataSub]]$vCovars, list(base))
		n <- length(data[[dataSub]]$vCovars)
		names(data[[dataSub]]$vCovars)[n] <- dname
		data
	}
	if (getDocumentation)
	{
		return(getInternals())
	}
	nGroups <- max(unique(effects$group))
	# note that period is a character and must be changed to numeric
	groupPeriods <- as.numeric(tapply(as.numeric(effects$period), effects$group,
			max, na.rm=TRUE)) + 1
	intervals <- rep(1:length(groupPeriods), groupPeriods)
	localPeriodNos <- unlist(sapply(groupPeriods, function(x)1:x))
	nGroups <- length(groupPeriods)
	v1 <- do.call(c, lapply(groupPeriods, function(x)1:x))
	v2 <- rep(cumsum(c(0, groupPeriods))[1:nGroups], groupPeriods)
	v3 <- v1 + v2
	periodNos <- v3[!v3 %in% cumsum(groupPeriods)]

	if (!is.null(data))
	{
		atts <- attributes(data)
	}
	if ((length(periodNos) < 2) && any(effects$timeDummy != ","))
	{
		warning("Time dummies not relevant with only 2 periods")
		effects$timeDummy <- ","
	}

	structuralRate <- effects$type == "rate" & effects$rateType %in%
		"structural"
	if (any(effects$timeDummy[structuralRate] != ","))
	{
		warning("Time dummy effects are not implemented",
			" for structural rate effects.")
		effects$timeDummy[structuralRate] <- ","
	}
	# JAL: Implementing these 20-FEB-2011 in RSienaTest
	# TODO: Behavioral interactions need to be implemented for these to work.
	# Once they are, we can comment lines 575-577 and 717-720 out.
	##	behaviorNonRateX <- effects$netType =="behavior" & effects$type != "rate"
	##   if (any(effects$timeDummy[behaviorNonRateX] != ","))
	##  {
	## 		warning("Time dummy effects are not implemented",
	##                " for behavior effects of type eval or endow.")
	##        effects$timeDummy[behaviorNonRateX] <- ","
	##    }

	if (all(effects$timeDummy == ",") )
	{
		##		No time dummy interactions to add, so kick the inputs back.
		return(list(effects=effects, data=data))
	}
	else
	{
		effects$timeDummy[effects$timeDummy=="all"]  <-
			paste(periodNos[-1], collapse = ",")

		alreadyDummied <- grepl("isDummy", effects$timeDummy)
		if (length(alreadyDummied)  >  0)
		{
			## Just remove those effects that have already been dummied so as to
			## not messy things up. The assumption is that the user will retain
			## all of the previous dummied effects within the column.
			effects <- effects[!alreadyDummied, ]
		}

		timesd <- strsplit(effects$timeDummy, ",| ")

		opar <- options(warn=-1)
		vals <- sapply(timesd, function(x) any(x=="" | !is.na(as.numeric(x))))
		if (any(!vals))
		{
			stop("Invalid input in timeDummy column", ": ",
				paste(effects$timeDummy[!vals], collapse=" and "))
		}
		options(opar)

		timesd <- lapply(timesd, function(x)as.numeric(x[x %in% periodNos]))
		dummiedEffects <- sapply(timesd, function(x)length(x) > 0)

		baseType <- list(number=which(dummiedEffects),
			type=effects$netType[dummiedEffects],
			name=effects$name[dummiedEffects])

		rateXDummies <- effects$shortName == "RateX" & dummiedEffects

		newEffects <- effects

		##@getCovar internal sienaTimeFix find the covariate
		getCovar <- function(cdvind, vdvind, bdvind, p)
		{
			if (!is.na(cdvind))
			{
				tmp <- data[[dataSub]]$cCovars[[cdvind]]
				val <- tmp
			}
			else if (!is.na(vdvind))
			{
				tmp <- data[[dataSub]]$vCovars[[vdvind]]
				val <- tmp[, p]
			}
			else
			{
				tmp <- data[[dataSub]]$depvars[[bdvind]]
				val <- tmp[, , p]
			}
			list(covar=tmp, val=val)
		}
		for (i in which(rateXDummies))
		{
			effect <- newEffects[i, ]
			## Figure out the base values:
			if (!is.null(data))
			{
				cdvind <- match(effect$interaction1, atts$cCovars)
				vdvind <- match(effect$interaction1, atts$vCovars)
				bdvind <- match(effect$interaction1, atts$netnames)
				if (any((is.na(cdvind) & is.na(vdvind) &
					(is.na(bdvind))))) #|| atts$types[bdvind] != "behavior"))
				{
					stop("Having trouble finding the covariate for your rate ",
						"effect. Please contact the developers.")
				}
			}
			for (p in timesd[[i]])
			{
				dname <-  paste(effect$interaction1, "Dummy", p, ":",
					effect$name, sep="")
				dataSub <- intervals[p]
				pp <- localPeriodNos[p]
				if (!is.null(data))
				{
					nPer <- attr(data[[dataSub]]$depvars[[effect$name]],
						"netdims")[3]
					covar <- getCovar(cdvind, vdvind, bdvind, pp)
					## add a new varCovar:
					nodeSet <- attr(covar$covar, "nodeSet")
					base <- matrix(0, nrow=length(covar$val), ncol=nPer - 1)
					base[, pp] <- covar$val
					data <- addVarCovar(data, base, nodeSet, dname, dataSub)

					## also add to other data objects as all must have same set

					for (ii in 1:length(data))
					{
						if (ii != dataSub)
						{
							nActors <-
								length(data[[ii]]$nodeSets[[nodeSet]])
							nPer <-
								attr(data[[ii]]$depvars[[effect$name]],
									"netdims")[3]
							base <- matrix(0, nrow=nActors, ncol=nPer - 1)
							data <- addVarCovar(data, base, nodeSet, dname, ii)
						}
					}
				}
				## add the rate effect row:
				newEffects <-
					addEffect(newEffects, i, dname, "covarNonSymmetricRate",
						'RateX', paste('isDummy', p, i, sep=','))
			}
		}
		## now the non rateX effects

		## first add dummies for all requested periods for each
		## dependent variable

		## then construct the interaction effects with the other effects
		## the dummies are fixed unless we need dummies for the density effect

		## find which periods we need
		use <- sapply(timesd, length) > 0 & !rateXDummies
		timeslist <- split(timesd[use], effects$name[use])
		timeslist <- lapply(timeslist, function(x)sort(unique(unlist(x))))
		timesTypelist <- split(timesd[use],
			list(effects$name[use], effects$type[use]))
		timesTypelist <- lapply(timesTypelist,
			function(x)unique(unlist(x)))

		types <- unique(effects$type[use ])
		for (depvar in names(timeslist))
		{
			##if (!is.null(data) && atts$types[[depvar]] == "behavior")
			##{
			##    stop ("Function is not specified for behavior effects")
			##}
			for (p in timeslist[[depvar]])
			{
				## create the dummy covariate
				pp <- localPeriodNos[p]
				dname <- paste("Dummy", p, ":", depvar, sep="")
				dataSub <- intervals[p]
				if (!is.null(data))
				{
					dims <- attr(data[[dataSub]]$depvars[[depvar]],
						"netdims")
					nodeSet <- attr(data[[dataSub]]$depvars[[depvar]],
						"nodeSet")[1]
					nActors <- dims[1]
					nPer <- dims[3]
					tmp <- matrix(0, nActors, nPer - 1)
					tmp[, pp] <- 1
					data <- addVarCovar(data, tmp, nodeSet, dname, dataSub)

					## also add to other data objects as all must be the same

					for (i in 1:length(data))
					{
						if ( i != dataSub)
						{
							dims <- attr(data[[i]]$depvars[[depvar]],
								"netdims")
							nodeSet <-
								attr(data[[i]]$depvars[[depvar]],
									"nodeSet")[1]
							nActors <- dims[1]
							nPer <- dims[3]
							base <- matrix(0, nActors, nPer - 1)
							data <- addVarCovar(data, base, nodeSet, dname, i)
						}
					}
				}
				## This is where the behaviors are treated differently
				if (baseType$type[baseType$name==depvar][1] == "behavior")
				{
					##find the linear effect rows for this depvar
					i <- which(effects$name == depvar &
						effects$shortName=="linear" &
						effects$type %in% types)

					if (length(i) == 0)
					{
						stop("Cannot find 'effect from' effect(s) for ", depvar,
							" ", types)
					}

					## establish whether we want the effFroms
					## included, fixed or not
					fix <- sapply(timesd[i], function(x)!p %in% x)
					typesNames <- paste(depvar, effects$type[i], sep=".")
					includeTypes <- !sapply(timesTypelist[typesNames], is.null)

					## add one or more effects (depending on types requested)
					newEffects <- addEffect(newEffects, i, dname,
						"covarBehaviorObjective", "effFrom",
						paste('isDummy', p,
							effects$effectNumber[i],
							sep=','),
						fix=fix, include=includeTypes,
						yName=dname)
				}
				else
				{
					##find the density rows for this depvar
					i <- which(effects$name == depvar &
						effects$shortName=="density" &
						effects$type %in% types)

					if (length(i) == 0)
					{
						stop("Cannot find density effect(s) for ", depvar,
							" ", types)
					}

					## establish whether we want the egoXs included,
					## fixed or not
					fix <- sapply(timesd[i], function(x)!p %in% x)
					typesNames <- paste(depvar, effects$type[i], sep=".")
					includeTypes <- !sapply(timesTypelist[typesNames], is.null)

					## add one or more effects (depending on types requested)
					newEffects <-
						addEffect(newEffects, i, dname,
							"covarNonSymmetricObjective", "egoX",
							paste('isDummy', p,
								effects$effectNumber[i], sep=','),
							fix=fix, include=includeTypes)
				}
			}

			## now the interactions for any dummied effects for this
			## dependent variable
			for (j in
				which(dummiedEffects & !rateXDummies &
					effects$name==depvar &
					effects$shortName != "density" &
					effects$shortName != "linear"))
			{
				for (p in timesd[[j]])
				{
					dname <- paste("Dummy", p, ":", depvar, sep="")
					effect <- effects[j, ]
					if (baseType$type[baseType$name==depvar][1] == "behavior") {
						newEffects <-
							includeInteraction(newEffects,
								effect$shortName, "effFrom",
								character=TRUE,
								type=effect$type,
								interaction1= c(effect$interaction1,
									dname),
								interaction2=effect$interaction2,
								name=depvar, verbose=FALSE)
					} else {
						newEffects <-
							includeInteraction(newEffects,
								effect$shortName, "egoX",
								character=TRUE,
								type=effect$type,
								interaction1= c(effect$interaction1,
									dname),
								interaction2=effect$interaction2,
								name=depvar, verbose=FALSE)
					}
					## find the row altered
					newrow <- newEffects$effect1 == j &
						newEffects$effect2 ==
						newEffects[paste(dname, effect$type, sep="."),
					"effectNumber"]
						newEffects$timeDummy[newrow] <-
							paste('isDummy', p, j, sep=',')
				}
			}
		}
		## reorder the effects so Dummies come before interactions and
		## rates where they are expected
		e0 <- split(newEffects, newEffects$name)

		e0 <- lapply(e0, function(x)
			{
				e1 <- !grepl("unspecified interaction effect",
					x$effectName)
				e2 <- grepl("unspecified interaction effect",
					x$effectName)
				e3 <- x$type == "rate"
				rbind(x[e3, ], x[e1 & ! e3, ], x[e2 & ! e3, ])
			})
		## sort into order of names in original effect
		effNames <- unique(effects$name)
		order1 <- match(effNames, names(e0))
		e0 <- e0[order1]
		newEffects <- do.call(rbind, e0)
		## now reconstruct the group object to make sure all attributes are
		## correct after adding covariates
		if (!is.null(data))
		{
			data <- sienaGroupCreate(data, singleOK=TRUE)
		}
		list(effects=newEffects, data=data)
	}
}
##@includeTimeDummy DataCreate
includeTimeDummy <- function(myeff, ..., timeDummy="all", name=myeff$name[1],
	type="eval", interaction1="", interaction2="", include=TRUE,
	character=FALSE)
{
	if (character)
	{
		dots <- sapply(list(...), function(x)x)
	}
	else
	{
		dots <- substitute(list(...))[-1] ##first entry is the word 'list'
	}
	if (length(dots) == 0)
	{
		stop("need some effect short names")
	}
	if (!character)
	{
		effectNames <- sapply(dots, function(x)deparse(x))
	}
	else
	{
		effectNames <- dots
	}
	use <- myeff$shortName %in% effectNames &
		myeff$type==type &
		myeff$name==name &
		myeff$interaction1 == interaction1 &
		myeff$interaction2 == interaction2
	if (any(use))
	{
		myeff[use, "timeDummy"] <- timeDummy
		myeff[use, "include"] <- include
		print.data.frame(myeff[use,])
	}
	else
	{
		warning("No effects could be found with this specification.")
	}
	myeff
}
