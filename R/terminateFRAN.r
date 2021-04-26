terminateFRAN <- function(z, x)
{
    if (z$cconditional)
    {
        z$rate<- colMeans(z$ntim, na.rm=TRUE)
        z$vrate <- apply(z$ntim, 2, sd, na.rm=TRUE)
        z$theta[z$posj] <- z$theta[z$posj] * z$rate
		if (!x$simOnly)
		{
		  if (!gmm(x))
		  {
		    z$covtheta[z$posj, ] <- z$covtheta[z$posj, ] * z$rate
		    z$covtheta[, z$posj] <- z$covtheta[,z$posj ] * z$rate
		  }
			else 
			{ 
			  z$covtheta[z$posj[-which(z$requestedEffects$type=="gmm")], ] <- z$covtheta[z$posj[-which(z$requestedEffects$type=="gmm")], ] * z$rate
			  z$covtheta[, z$posj[-which(z$requestedEffects$type=="gmm")]] <- z$covtheta[,z$posj[-which(z$requestedEffects$type=="gmm")] ] * z$rate  
			}
		}	
    }
    f <- FRANstore()
    f$pModel <- NULL
    f$pData <- NULL
    z$f$myeffects <- lapply(z$f$myeffects, function(x){x$effectPtr <- NULL;x})
    FRANstore(NULL) ## clear the stored object
    if (is.null(z$print) || z$print)
    {
        PrintReport(z, x)
    }
    if (sum(z$test))
    {
        z$fra <- colMeans(z$sf, na.rm=TRUE)
		z$redundant <- (z$fix & (!z$test))
        ans <- ScoreTest(z$pp, z$dfra, z$msf, z$fra, z$test, z$redundant, x$maxlike)
        z <- c(z, ans)
        TestOutput(z, x)
    }
    if (!is.null(z$dfra))
    {
        dimnames(z$dfra)[[1]] <- as.list(z$requestedEffects$shortName)
    }
    return(z)
}
