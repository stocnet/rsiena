#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: Sienatest.r
# *
# * Description: This module contains the function for instability analysis and
# * score tests.
# *
# *****************************************************************************/

test_parameter <- function(x, ...) UseMethod("test_parameter", x)

##@test_parameter.sienaFit test_parameter method for all tests for sienaFit 
test_parameter.sienaFit <- function(x, method=NULL, 
										tested=NULL, tested2=NULL, ...)
{
	if (is.null(method))
	{
		if (inherits(tested,"matrix")) 
		{
			Wald.RSiena(tested, x)
		}
		else if (is.vector(tested)) 		
		{
			Multipar.RSiena(x, tested)
		}
		else if (tested=="all")
		{
			Multipar.RSiena(x, which(!x$test))
		}
		else
		{
			stop("keyword 'tested' is not given properly")
		}
	}
	else if (method=="score")
	{
		if (is.null(tested))
		{
			tested <- x$test
		}
		score.Test(x, tested=tested)
	}
	else if (method=="same")
	{
		testSame.RSiena(x, tested, tested2)
	}
	else
	{
		stop(paste("There is no test_parameter.sienaFit with method = ", method))
	}
}

##@TestOutput siena07 Print test report
TestOutput <- function(z, x)
{
	testn <- sum(z$test)
	## browser()
	if (testn)
	{
		if (x$maxlike)
		{
			Heading(2, outf,'Score test <c>')
		}
		else
		{
			Heading(2, outf, 'Generalised score test <c>')
		}
		Report('Testing the goodness-of-fit of the model restricted by\n', outf)
		j <- 0
		for (k in 1:z$pp)
		{
			if (z$test[k])
			{
				j <- j + 1
				Report(c(" (", j, ")   ",
						format(paste(z$requestedEffects$type[k], ":  ",
								z$requestedEffects$effectName[k],
								sep=""),
							width=50), " = ",
						sprintf("%8.4f", z$theta[k]), "\n"),
					sep = "", outf)
			}
		}
		Report('_________________________________________________\n',outf)
		Report('                ',outf)
		Report('   \n',outf)
		if (testn > 1)
			Report('Joint test:\n-----------\n',outf)
		Report(c('   c = ',sprintf("%8.4f", z$testresOverall),
				'   d.f. = ',j,'   p-value '),sep='',outf)
		pvalue <- 1-pchisq(z$testresOverall,j)
		if (!is.na(pvalue))
		{
			if (pvalue < 0.0001)
			{
				Report('< 0.0001',outf)
			}
			else
			{
				Report(c('= ',sprintf("%8.4f",pvalue)), sep = '', outf)
			}
		}
		else
		{
			Report('  NA  ',outf)
		}
		if (testn==1)
			Report(c('\n   one-sided (normal variate): ',
					sprintf("%8.4f",z$testresulto[1])), sep = '', outf)
		if (testn> 1)
		{
			Report('\n\n',outf)
			for (k in 1:j)
			{
				Report(c('(',k,') tested separately:\n'),sep='',outf)
				Report('-----------------------\n',outf)
				Report(' - two-sided:\n',outf)
				Report(c('  c = ', sprintf("%8.4f", z$testresult[k]),
						'   d.f. = 1  p-value '), sep = '', outf)
				pvalue<- 1-pchisq(z$testresult[k],1)
				if (!is.na(pvalue))
				{
					if (pvalue < 0.0001)
					{
						Report('< 0.0001\n',outf)
					}
					else
					{
						Report(c('= ', sprintf("%8.4f", pvalue), '\n'), sep = '',
							outf)
					}
				}
				else
				{
					Report('  NA  ',outf)
				}
				Report(c(' - one-sided (normal variate): ',
						sprintf("%8.4f", z$testresulto[k])), sep = '', outf)
				if (k < j)
				{
					Report('\n\n',outf)
				}
			}
		}
		Report('    \n_________________________________________________\n\n',
			outf)
		Report('One-step estimates: \n\n', outf)
		for (i in 1 : z$pp)
		{
			onestepest <- z$oneStep[i] + z$theta[i]
			Report(c(format(paste(z$requestedEffects$type[i], ':  ',
							z$requestedEffects$effectName[i], sep = ''),
						width=50),
					sprintf("%8.4f", onestepest), '\n'), sep = '', outf)
		}
		Report('\n',outf)
	}
}

##@ScoreTest siena07 Do score tests
ScoreTest<- function(pp, dfra, msf, fra, test, redundant, maxlike)
{
	testresult <- rep(NA, pp) ##for chisq per parameter
	testresulto <- rep(NA, pp) ##for one-sided tests per parameter
	##first the general one
	ans <- EvaluateTestStatistic(maxlike, test, redundant, dfra, msf, fra)
	testresOverall <- ans$cvalue
	covMatrix <- ans$covMatrix
	if (sum(test) == 1)
	{
		testresulto[1] <- ans$oneSided
	}
	else
	{
		## single df tests
		use <- !test
		k <- 0
		for (i in 1:pp)
		{
			if (test[i])
			{
				k <- k + 1
				use[i] <- TRUE
				ans <- EvaluateTestStatistic(maxlike, test[use], redundant[use],
					dfra[use, use], msf[use, use], fra[use])
				testresult[k] <- ans$cvalue
				testresulto[k] <- ans$oneSided
				use[i] <- FALSE
			}
		}
	}
	##onestep estimator
	if (maxlike)
		dfra2 <- dfra + msf
	else
		dfra2 <- dfra
	if (inherits(try(dinv2 <- solve(dfra2), silent=TRUE), "try-error"))
	{
		Report("Error message for inversion to get onestep estimator: \n", cf)
		dinv2 <- dfra2
		dinv2[] <- NA
		oneStep <- rep(NA, nrow(dfra2))
	}
	else
	{
		oneStep<- -dinv2 %*% fra
	}
	list(testresult=testresult, testresulto=testresulto,
		testresOverall=testresOverall, covMatrix=covMatrix,
		oneStep=oneStep, dinv2= dinv2, dfra2=dfra2)
}
##@EvaluateTestStatistic siena07 Calculate score test statistics
EvaluateTestStatistic<- function(maxlike, test, redundant, dfra, msf, fra)
{
	##uses local arrays set up in the calling procedure
	d11 <- dfra[!(test|redundant), !(test|redundant), drop=FALSE]
	d22 <- dfra[test, test, drop=FALSE]
	d21 <- dfra[test, !(test|redundant), drop=FALSE]
	d12 <- t(d21)
	sigma11 <- msf[!(test|redundant), !(test|redundant), drop=FALSE]
	sigma22<- msf[test, test,drop=FALSE]
	sigma12 <- msf[!(test|redundant), test, drop=FALSE]
	sigma21<- t(sigma12)
	z1 <- fra[!(test|redundant)]
	z2 <- fra[test]
	if (inherits(try(id11 <- solve(d11), silent=TRUE), "try-error"))
	{
		warning('Score test: Error for inversion of d11 \n')
		oneSided <- NA
		v9 <- d22
		v9[] <- NA
		cvalue <- matrix(NA, 1, 1)
	}
	else
	{
		rg <- d21 %*% id11
		if (!maxlike)
		{
			##orthogonalise deviation vector
			ov <- z2 - rg %*% z1
			##compute var(ov) = sigma22 - (d21 %*% id11) %*% sigma12 -
			##      sigma21 %*% t(id11)%*% t(d21) +
			##      d21%*%id11 %*% sigma11 %*% t(id11) %*% t(d21)
			v2 <- sigma21 - rg %*% sigma11
			v6 <- v2 %*% t(id11) %*% t(d21)
			v9 <- sigma22 -  rg %*% sigma12 - v6
		}
		else
		{
			ov <- -z2
			v9 <- d22 - rg %*% d12
		}
		if (inherits(try(vav <- solve(v9), silent=TRUE), "try-error"))
			## vav is the inverse variance matrix of ov
		{
			warning('Score test: Error for inversion of v9\n')
			vav <- v9
			vav[] <- NA
			cvalue <- NA
			oneSided <- NA
		}
		else
		{
			cvalue <- t(ov) %*% vav %*% ov
			if (cvalue < 0) cvalue <- 0
			if (sum(test) == 1)
			{
				if (vav > 0)
					oneSided <- ov * sqrt(vav)
				else
					oneSided <- 0
				if (!maxlike) oneSided <- - oneSided
				## change the sign for intuition for users
			}
			else
			{
				oneSided <- 0
			}
		}
	}
	list(cvalue=cvalue, oneSided=oneSided, covMatrix=v9)
}


##@theEfNames internal RSiena function getting estimated effect names
theEfNames <- function(ans){
	if(length(unique(ans$requestedEffects$name)) > 1)
	{
		res <- paste(ans$requestedEffects$name, ans$requestedEffects$effectName, sep=': ')
	}
	else
	{
		res <- ans$requestedEffects$effectName
	}
	if (any(ans$requestedEffects$type != "eval"))
	{
		res <- paste(res, ans$requestedEffects$type, sep=' ')
	}
	res
}

##@scoreTest Calculate score tested
score.Test <- function(x, tested=x$test)
	# use: x must be a sienaFit object;
	# tested must be a boolean vector with length equal to the number of parameters of x,
	# or a vector of integer numbers between 1 and x$pp.
{
	if ((is.numeric(tested)) || (is.integer(tested)))
	{
		if (max(tested) > x$pp)
		{
			stop(paste('The maximum requested coordinate is too high:',
					max(tested)))
		}
		tested <- (1:x$pp) %in% tested
	}
	if (sum(tested) <= 0) stop(paste('Something should be tested, but the total requested is',
			sum(tested)))
	if (length(tested) != x$pp) stop('Dimensions of tested must agree')
	if (any(tested & (!x$fix))) warning('Warning: some tested parameters were not fixed; do you know what you are doing??? \n')
	fra <- colMeans(x$sf, na.rm=TRUE)
	redundant <- (x$fix & (!tested))
	tests <- EvaluateTestStatistic(x$maxlike, tested, redundant, x$dfra, x$msf, fra)
	teststat <- tests$cvalue
	df <- sum(tested)
	if (df == 1)
	{
		onesided <- tests$oneSided
	}
	else
	{
		onesided <- NULL
	}
	pval <- 1 - pchisq(teststat, df)
	efnames <- theEfNames(x)[tested]
	t_x <- list(chisquare=teststat, df=df, pvalue=pval, onesided=onesided, efnames=efnames)
	class(t_x) <- "sienaTest"
	attr(t_x, "version") <- packageDescription(pkgname, fields = "Version")
	t_x
}

##@Wald.RSiena  Calculate Wald test statistics
Wald.RSiena <- function(A, x)
{
	if (is.vector(A)){A <- matrix(A, nrow=1)}
	if (dim(A)[2] != x$pp){stop(paste('A must have', x$pp, 'columns.'))}
	sigma <- x$covtheta
	if (any(is.na(sigma))) {
		# happens when some coordinates were fixed
		# in the call of siena07 leading to x;
		# then the non-used part of sigma,
		# which partially consists of NA,
		# is replaced by the identity matrix.
		zero.cols <- apply(A, 2, function(colum){all(colum==0)})
		sigma[zero.cols, ] <- 0
		sigma[, zero.cols] <- 0
		diag(sigma)[zero.cols] <- 1
	}
	th <- A %*% x$theta
	covmat <- A %*% sigma %*% t(A)
	chisq <- drop(t(th) %*% solve(covmat) %*% th)
	df <- nrow(A)
	pval <- 1 - pchisq(chisq, df)
	if (df == 1)
	{
		onesided <- sign(th) * sqrt(chisq)
	}
	else
	{
		onesided <- NULL
	}
	t_x <- list(chisquare=chisq, df=df, pvalue=pval, onesided=onesided)
	class(t_x) <- "sienaTest"
	attr(t_x, "version") <- packageDescription(pkgname, fields = "Version")
	t_x
}

##@Multipar.RSiena  Calculate Wald test statistic for hypothesis that subvector = 0.
Multipar.RSiena <- function(x, tested)
{
	p <- length(x$theta)
	efnames <- theEfNames(x)[tested]
	k <- length(tested)
	A <- matrix(0, nrow=k, ncol=p)
	A[cbind(1:k,tested)] <- 1
	t_x <- Wald.RSiena(A, x)
	t_x$efnames <- efnames
	t_x
}

##@testSame.RSiena Test that two subvectors of parameter are identical
testSame.RSiena <- function(x, e1, e2){
	if (length(e1) != length(e2))
	{
		stop('Tested parameters should have the same length')
	}
	p <- x$pp
	if ((min(e1) < 1) | (min(e2) < 1) | (max(e1) > p) | (max(e2) > p))
	{
		stop(paste('Tested parameters should have numbers between 1 and', p))
	}
	if (any(x$fixed[c(e1,e2)]))
	{
		stop('Fixed parameters cannot be tested')
	}
	t <- length(e1)
	efnames1 <- theEfNames(x)[e1]
	efnames2 <- theEfNames(x)[e2]
	if (t > 1)
	{
		max.width <- max(c(nchar(efnames1)))
		efnames1 <- format(efnames1, width=max.width)
	}
	efnames <- paste(efnames1, ' == ', efnames2)
	mat <- matrix(0,t,p)
	for (i in seq_along(1:t)){mat[i,e1[i]] <- 1}
	for (i in seq_along(1:t)){mat[i,e2[i]] <- -1}
	t_x <- Wald.RSiena(mat, x)
	t_x$efnames <- efnames
	t_x
}

##@print.sienaTest Methods
print.sienaTest <- function(x, ...)
{
	if (!inherits(x, "sienaTest"))
	{
		stop("not a legitimate sienaTest object")
	}
	if (any(is.na(unlist(x))))
	{
		stop("some estimates or standard errors are NA")
	}
	if (!is.null(x$efnames))
	{
		cat('Tested effects:\n ')
		cat(paste(x$efnames,'\n'))
	}
	cat(paste('chi-squared = ',
		format(round(x$chisquare, digits=2), nsmall = 2),
		', d.f. = ', x$df, '; ', sep=''))
	if ((x$df == 1) & (!is.null(x$onesided)))
	{
		cat(paste('one-sided Z = ',
			format(round(x$onesided, digits=2), nsmall = 2), '; ', sep=''))
		cat('two-sided')
	}
	if (x$pvalue < 0.001)
	{
		cat(' p < 0.001. \n')
	}
	else
	{
		cat(paste(' p = ',
			format(round(x$pvalue, digits=3), nsmall = 2), '. \n', sep=''))
	}
	invisible(x)
}
