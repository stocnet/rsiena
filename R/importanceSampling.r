
stdError <- function(z, x, subphase,...)
{
	z$addChainToStore <- TRUE
	z$nbrNodes <- 1
	z$storedChains <- FALSE
	
	theta <- z$thav/z$thavn
	
	N <- length(z$myloglik2)

	z$iterSeq <- (N:1)
	z$muLik <- getLikelihoods(theta, z, getScores=TRUE,
		iterSequence= z$iterSeq)

	z
}

entropy <- function(x)
{
	y<- x/sum(x)
	l <- length(y)
	
	-sum(y*log(y))/log(l)
}

E <- function(t, z, x)
{
	theta <- t*x$thetaEnd + (1-t)*z$theta
	cat(paste("\n",theta,"\n"))
	z$muLik1 <- getLikelihoods(theta, z, getScores=TRUE,
	iterSequence= z$iterSeq)
	z$w1 <- exp(z$muLik1$lik - z$myloglik)
	z$bridgeEntropy <- entropy(z$w1)
	
	z
}

bridge1 <- function(z, x, subphase,...)
{
	z$addChainToStore <- TRUE
	z$nbrNodes <- 1
	z$storedChains <- FALSE
	
	z$iterSeq <- (z$nit:1)
	
	converge <- TRUE
	
	if (ifelse(is.null(x$checkTstat), FALSE, x$checkTstat))
	{
		tstat <- apply(z$sf,2,mean)/apply(z$sf,2,sd)
		if (max(abs(tstat)) > 0.2) ### don't want to apply this if an effect is fixed?
		{
			converge <- FALSE
		}
	}
	
if (converge)
{
	c <- ifelse(is.null(x$c), 0.8, x$c)
	z$min1 <- 0
	z$max1 <- 1

	if (ifelse(is.null(x$finalBridge), TRUE, !x$finalBridge))
	{

		z <- E(1,z,x)

		if (z$bridgeEntropy > c)
		{
			z$max1 <- 1
		} else {
			
			for (i in 1:6)
			{
				z$t1 <- (z$min1 + z$max1)/2
				z <- E(z$t1,z,x)
				cat("t1: ");cat(z$t1);cat("\n")
				cat("entropy: ");cat(z$bridgeEntropy);cat("\n")
				if (is.na(z$bridgeEntropy)) {break}
		
				if (z$bridgeEntropy < c)
				{
					z$max1 <- 0.5*(z$min1 + z$max1)
				} else {
					z$min1 <- 0.5*(z$max1 + z$min1)
				}
			}
	
			z$max1 <- z$t1
		}

		z$thetaNext <- z$max1*x$thetaEnd + (1-z$max1)*z$theta

		if (!is.null(z$thetaNext))
		{
			newLik <- getLikelihoods(z$thetaNext, z, getScores=TRUE,
													iterSequence = z$iterSeq)
			z$loglikNext <- newLik$lik
		}
	}

	if (!is.null(x$thetaPrevious))
	{
		newLik <- getLikelihoods(x$thetaPrevious, z, getScores=TRUE,
													iterSequence = z$iterSeq)

		z$loglikPrevious <- newLik$lik
	}
	
}

	z
}

eThermo <- function(h,H, z, x)
{
	if (!is.null(z$likelihoodsList[[h+1]]))
	{
		cat("\nAlready have likelihoods\n")
		z$muLik1 <- z$likelihoodsList[[h+1]]
	} else {
		cat("\nNeed to calculate likelihoods\n")
		theta <- (h/H)*x$thetaEnd + (1-(h/H))*z$theta
		z$muLik1 <- getLikelihoods(theta, z, getScores=TRUE,
		iterSequence= z$iterSeq)
		z$likelihoodsList[[h+1]] <- z$muLik1
	}

	z$w1 <- exp(z$muLik1$lik - z$myloglik)
	z$bridgeEntropy <- entropy(z$w1)
	
	z
}

zThermo <- function(z, x)
{
	z$r <- rep(0, length(z$theta))
	
	for (h in z$hmin:z$hmax)
	{
		cat("h: ");cat(h);cat("\n")
		z$r <- z$r + rThermo(z,x,h)
		cat("z$r: ");cat(z$r);cat("\n")
	}
	
	sum((x$thetaStart - x$thetaEnd)*z$r)/(z$H + 1)
}

rThermo <- function(z,x,h)
{
	r <- rep(0, z$pp)
	z <- E(h/z$H,z,x)
	scores <- t(z$muLik1$sc)

	for (i in 1:z$pp)
	{
		r[i] <- mean(z$w1*scores[,i])/mean(z$w1)
	}
	
	r
}

thermo <- function(z, x, subphase,...)
{
	z$addChainToStore <- TRUE
	z$nbrNodes <- 1
	z$storedChains <- FALSE

	z$iterSeq <- (z$nit:1)
	
	#### Interval bisection ####
	
	c <- ifelse(is.null(x$c), 0.8, x$c)

	z$H <- ifelse(is.null(x$H), 10, x$H)
	z$hmin <- ifelse(is.null(x$hmin), 0, x$hmin)
	z$hmaxMax <- z$H
	z$hmaxMin <- z$hmin
	
	z$likelihoodsList <- vector("list",z$H+1)

	z <- eThermo(z$hmaxMax,z$H,z,x)
	z$likelihoodsList[[z$H+1]] <- z$muLik1
	
	if (z$bridgeEntropy > c)
	{
		z$hmax <- z$H
	} else {
			
		for (i in 1:6) 
		{
			if (z$hmaxMin < z$hmaxMax)
			{
				z$t1 <- floor((z$hmaxMin + z$hmaxMax)/2)
				z <- eThermo(z$t1,z$H,z,x)
				cat("t1: ");cat(z$t1);cat("\n")
				cat("entropy: ");cat(z$bridgeEntropy);cat("\n")
				if (is.na(z$bridgeEntropy)) {break}
			
				if (z$bridgeEntropy < c)
				{
					z$hmaxMax <- ceiling(0.5*(z$hmaxMin + z$hmaxMax))
				} else {
					z$hmaxMin <- floor(0.5*(z$hmaxMax + z$hmaxMin))
				}
				cat("hmaxMin:");cat(z$hmaxMin);cat("\n")
				cat("hmaxMax:");cat(z$hmaxMax);cat("\n")
			}
		}
		
		z$hmax <- z$t1
	}
	
	cat("hmin:");cat(z$hmin);cat("\n")
	cat("hmax:");cat(z$hmax);cat("\n")
	
	z$zThermo <- zThermo(z,x)
	
	z
}


