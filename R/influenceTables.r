##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: https://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: influenceTables.r
## *
## ****************************************************************************/

#############################################################################
# R script InfluenceTables.r                                                #
# written by Tom A.B. Snijders                                              #
# October 12, 2023                                                          #
# (Gratitude to Steffen Triebel and Rene Veenstra for corrections)          #
#############################################################################

#############################################################################
#                                                                           #
# These are functions for constructing and presenting influence tables      #
# for the interpretation of results for network and behavior dynamics       #
# obtained with the RSiena or multiSiena packages.                          #
# Also consult the manual, s                         #
#                                                                           #
# The main function to use is                                               #
# influenceMatrix <- function(x,xd,netname,behname,levls, levls.alt=levls)  #
# which creates a matrix containing the influence table for                 #
# siena data set xd, sienaFit or sienaMeta or sienaBayesFit object x,       #
# network variable netname (should be a character string),                  #
# dependent actor variable behname (also a character string),               #
# levels for ego levls, levels for alter levls.alt.                         #
# If levls is NULL (the default), it is taken as the integer range of       #
# the dependent actor variable.                                             #
# Mostly levls.alt will be the same as levls,                               #
# in which case it does not have to be specified.                           #
#                                                                           #
# For models utilizing creation and/or endowment effects,                   #
# by specifying include.endow=TRUE                                          #
# the sum of evaluation and endowment effects is used;                      #
# by specifying include.creation=TRUE                                       #
# the sum of evaluation and creation effects is used;                       #
# the first is relevant for maintaining, i.e., not decreasing,              #
# the values of the dependent variable,                                     #
# the second for increasing them.                                           #
#                                                                           #
# Other functions to use (see below):                                       #
# influenceTable.se                                                         #
# Computes standard error of a linear combination of elements               #
# of the influence table.                                                   #
# Not for sienaMeta or sienaBayesFit objects.                               #
#                                                                           #
# influenceTable.plot                                                       #
# Constructs a plot of the influence table using ggplot2.                   #
#                                                                           #
# If these functions are used for a sienaMeta object x,                     #
# (produced by siena08, as opposed to a sienaFit object produced by siena07)#
# xd should be one of the individual data sets used for creating x,         #
# (a single-group instead of multi-group data set)                          #
# with an average ('representative') value                                  #
# for the mean of variable <behname>.                                       #
#                                                                           #
#############################################################################

##@influenceTable.basis  Basis for influenceTable
influenceTable.basis <- function(x, xd, netname, behname,
                            levls=NULL, levls.alt=levls,
                            out.ego=1, silent=FALSE, nfirst=x$nwarm+1,
                            include.endow=FALSE, include.creation=FALSE){
# Creates basic materials for a influence table for
# siena data set xd, sienaFit or sienaBayesFit object x,
# network variable netname (should be a character string),
# dependent actor variable behname (also a character string),
# levels for ego levls, levels for alter levls.alt.
# For effect avAlt, levls.alt refers to average values of alter's behavior;
# for effect totAlt, levls.alt refers to total values of alter's behavior;
# For effects avSim and totSim, levls.alt refers to
#                    constant values of alter's behavior.
# For effects totSim and totAlt, out.ego is the presumed outdegree of ego
# for which the contributions to the objective function are calculated.
#
# This function can be used also for a sienaMeta object x.
# Then xd must be one of the data sets used for creating x,
# with an average ('representative') value for the mean of variable <behname>.
#
	if(!silent)
	{
		cat("Network",netname, "; dependent behavior",behname,".\n")
	}
    if (!(inherits(x, "sienaFit") | inherits(x, "sienaBayesFit")))
    {
        stop("x should be a sienaFit or sienaBayesFit object")
    }
    if (!inherits(xd, "siena"))
    {
        stop("xd should be a siena data set")
    }
    if (include.endow && include.creation)
    {
        warning(paste('It is not meaningful to include creation',
                ' and maintenance effects simultaneously.'))
    }
# Obtain estimate
    if (inherits(x, "sienaBayesFit"))
    {
        theta <- sienaFitThetaTable(x, fromBayes=TRUE,
                    groupOnly=0, nfirst=nfirst)$mydf$value
        theEffects <- x$requestedEffects
        if (length(theta) != dim(theEffects)[1])
        {
            stop('mismatch between theta and effect names')
        }
# Is this correct for all model specifications?
    }
    else
    {
        theta <- x$theta
        theEffects <- x$requestedEffects
    }
    theparm <- theEffects$parm
# Note that theta and theparm now have the same length.
# Obtain means
    if (inherits(xd, "sienaGroup"))
    {
        cat("A sienaGroup data object was given.\n")
        # Does the variable exist, and is it a dependent variable?
        thebeh <- xd[[1]]$depvars[[behname]]
        thenet <- xd[[1]]$depvars[[netname]]
    }
    else
    {
        thebeh <- xd$depvars[[behname]]
        thenet <- xd$depvars[[netname]]
    }
    if (is.null(thenet)){stop(paste('There is no network <',netname,'>.'))}
    if (is.null(thebeh)){
            stop(paste('There is no dependent behaviour variable <',behname,'>.'))}

    if (inherits(xd, "sienaGroup"))
    {
#       zmean  <- mean(sapply(xd, function(z){mean(z$depvars[[behname]], na.rm=TRUE)}))
        ztot   <- sum(sapply(xd, function(z){sum(z$depvars[[behname]], na.rm=TRUE)}))
        ztotN  <- sum(sapply(xd, function(z){sum(!is.na(z$depvars[[behname]]))}))
        zmean  <- ztot/ztotN
        zsmean <- attr(xd, "bSim")[[behname]]
        Delta  <- attr(xd, "behRange")[,behname]
        if (is.null(levls))
        {
            levls <- Delta[1] : Delta[2]
        }
        Delta  <- Delta[2] - Delta[1]
    }
    else
    {
        zmean  <- mean(colMeans(thebeh, na.rm=TRUE))
        zsmean <- attr(thebeh, 'simMean')
        Delta  <- attr(thebeh, 'range')
        if (is.null(levls))
        {
            range2 <- attr(thebeh, 'range2')
            levls <- range2[1] : range2[2]
        }
    }
    if (is.null(levls.alt))
    {
        levls.alt <- levls
# I'm not sure this works if the effect used is totAlt and out.ego is not 1.
    }
# perhaps attr(thebeh, 'simMean') fails if there are more than one dependent network?
# is this attribute then a vector of length >= 2?
# Note that zsmean is used only if the model includes avSim or totSim.
    replace0 <- function(k){ifelse(length(k)==0,0,k)}
    efNames <- c('linear','quad','avAlt','avSim','totAlt','totSim',
                    'avAttHigher', 'avAttLower',
                    'threshold', 'threshold2', 'threshold3', 'threshold4')
    zeff.eval <- sapply(efNames, function(s)
        {replace0(which(theEffects$name == behname &
                (theEffects$type == 'eval') &
                (theEffects$interaction1 %in% c(netname,'')) &
                theEffects$include &
                theEffects$shortName==s))})
    zeff.endow <- sapply(efNames, function(s)
        {replace0(which(theEffects$name == behname &
                (theEffects$type == 'endow') &
                (theEffects$interaction1 %in% c(netname,'')) &
                theEffects$include &
                theEffects$shortName==s))})
    zeff.creation <- sapply(efNames, function(s)
        {replace0(which(theEffects$name == behname &
                (theEffects$type == 'creation') &
                (theEffects$interaction1 %in% c(netname,'')) &
                theEffects$include &
                theEffects$shortName==s))})
    zeff.parm <- sapply(efNames, function(s)
        {replace0(which(theEffects$name == behname &
                (theEffects$interaction1 %in% c(netname,'')) &
                theEffects$include &
                theEffects$shortName==s))})
# zeff gives the indicators of the effects in x
    taketheta <- function(k){ifelse(k==0,0,theta[k])}
    takeparm <- function(k){ifelse(k==0,0,theparm[k])}
    ztheta <- sapply(zeff.eval,taketheta)
    zparm <- sapply(zeff.parm,takeparm)
    if (include.endow)
    {
        ztheta <- ztheta + sapply(zeff.endow,taketheta)
    }
    if (include.creation)
    {
        ztheta <- ztheta + sapply(zeff.creation,taketheta)
    }
    if (all(ztheta == 0))
    {
        cat('All parameters found for the effect of', netname, 'on',
                                                behname, 'are 0.\n')
    stop('there seems to be something wrong with x, xd, netname, or behname.')
    }
    if (!silent){
        cat('Parameters found are\n')
        print(noquote(format(round(ztheta[ztheta != 0.0],4))))
        if (include.endow)
        {
            cat('This includes evaluation and endowment (maintenance) effects.\n')
        }
        if (include.creation)
        {
            cat('This includes evaluation and creation effects.\n')
        }
    }
    found <- (ztheta != 0.0)
    if (!silent){
		if (any(found[9:12]))
		{
			cat("Threshold internal effect parameters are ",
                        zparm[intersect(which(found), 9:12)], '\n')
		}
		if (sum(found[3:8]) == 0)
		{
			cat("Note: no influence effect was found.\n")
		}
		if (sum(found[3:8]) >= 2)
		{
			cat("Note: more that one influence effect was found.\n")
        }
		if (found[5] | found[6])
		{
			cat("Note that influence in this model depends on ego's outdegree. \n")
			cat("Outdegree of ego is given as ", out.ego, "; can be changed.\n")
        }
		if (found[4] | found[6] | found[7] | found[8])
		{
			cat("Levels of alter refer to constant values of alter's behavior.\n\n")
        } 
		else if (found[3] | found[5])
		{
			cat("Levels of alter refer to average values of alter's behavior.\n\n")
        }
		cat('\n')
		flush.console()
	}
    quad <- !any(found[c(4,6:12)])
    # These effects do not lead to quadratic functions
# ztheta contains the parameter values in x
    K <- length(levls)
    KA <- length(levls.alt)
    zalter <- rep(levls.alt,each=KA)
    zego <- rep(levls,K)
    fact <- 1:K
    alter <- factor(rep(fact,each=KA))
    coeffs <- matrix(NA, K*KA, length(efNames))
    coeffs[,1] <- (zego - zmean)
    coeffs[,2] <- (zego - zmean)*(zego - zmean)
    coeffs[,3] <- (zego - zmean)*(zalter - zmean)
    coeffs[,4] <- (1-(abs(zalter - zego)/Delta)-zsmean)
    coeffs[,5] <- out.ego*(zego - zmean)*(zalter - zmean)
    coeffs[,6] <- out.ego*(1-(abs(zalter - zego)/Delta)-zsmean)
    coeffs[,7] <- (1-(pmax(zalter - zego,0)/Delta)-zsmean)
    coeffs[,8] <- (1-(pmax(zego - zalter,0)/Delta)-zsmean)
    coeffs[,9] <- 1*(zego >= zparm[9])
    coeffs[,10] <- 1*(zego >= zparm[10])
    coeffs[,11] <- 1*(zego >= zparm[11])
    coeffs[,12] <- 1*(zego >= zparm[12])
    select <- coeffs %*% ztheta
    df <- data.frame(alter, zalter, zego, select=(coeffs%*%ztheta))
    list(df=df, zeff.eval=zeff.eval, zeff.endow=zeff.endow,
         zeff.creation=zeff.creation, ztheta=ztheta, zmean=zmean, zsmean=zsmean,
         Delta=Delta, coeffs=coeffs, levls=levls, levls.alt=levls.alt, quad=quad)
}

##@influenceTable  influenceTable
influenceTable <- function(x, xd, netname, behname,
                      as.matrix=FALSE,
                       levls=NULL, levls.alt=levls, out.ego=1, 
						silent=FALSE, nfirst=x$nwarm+1,
                        include.endow=FALSE, include.creation=FALSE){
# Creates a data frame containing the influence table for
# siena data set xd, sienaFit object x,
# network netname (should be a character string),
# dependent actor variable behname (also a character string),
# levels for ego levls, levels for alter levls.alt.
    infl.t <- influenceTable.basis(x, xd, netname, behname,
                   levls=levls, levls.alt=levls.alt, out.ego=out.ego, 
				   silent=silent, nfirst=nfirst,
                   include.endow=include.endow,
					include.creation=include.creation)
    df <- infl.t$df
    df$alter <- as.character(df$zalter)
    df$zego <- as.numeric(as.character(df$zego))
    df$select <- as.numeric(as.character(df$select))
    attr(df, "quad") <- infl.t$quad
	class(df) <- c("influenceTable", class(df))
	attr(df, "netname") <- netname
	attr(df, "behname") <- behname
	attr(df, "levls") <- infl.t$levls
	attr(df, "levls.alt") <- infl.t$levls.alt
	if (as.matrix)
	{
	    {
			levls <- unique(df$zego)
		}
		if (is.null(levls.alt))
		{
			levls.alt <- unique(df$alter)
		}
		mat <- matrix(as.numeric(as.character(df$select)),
            length(unique(df$zalter)),
            length(unique(df$zego)), byrow=TRUE)
		colnames(mat) <- infl.t$levls
		rownames(mat) <- infl.t$levls.alt
		class(mat) <- c("influenceTable", class(mat))
		return(mat)
	}
	else
	{
		return(df)
	}
}

influenceTable.se <- function(x, xd, netname, behname, levls, ww, levls.alt=levls){
# Calculates the standard error for a linear combination
# of elements of the influence table for siena data set xd, sienaFit object x,
# dependent actor variable behname (also a character string),
# network variable netname (should be a character string),
# levels for ego levls, levels for alter levls.alt.
# Coefficients of the linear combination are in matrix ww,
# which is assumed to have rows corresponding to levls and columns to levls.alt.
# The linear combination for which the standard error is computed is
# sum_{h,k} ww[h,k] * influenceTable[h,k].
    if (!all(dim(ww) == c(length(levls), length(levls.alt)))){
        stop('Dimension of ww should be the lengths of levls and levls.alt.\n')
    }
    if (inherits(x, "sienaBayesFit")){
        stop('This function does not work for sienaBayesFit objects.\n')
    }
    cat("Requested cell weights (row = alter's, col = ego's behavior value):\n")
    print(cbind(which(ww != 0, arr.ind=TRUE),
            value=ww[which(ww != 0, arr.ind=TRUE)]))
    if (!all(rowSums(ww) == 0)){
        cat("Warning: not all row sums of cell weights are 0.\n Row sums are:\n")
        print(rowSums(ww))
        cat("Only contrasts are meaningful. Change weight matrix.\n\n")
    }
    inftb <- influenceTable.basis(x=x, xd=xd, netname=netname, 
				behname=behname, levls=levls, levls.alt=levls.alt,
                silent=TRUE)
    zeff.r <- inftb$zeff.eval[inftb$zeff.eval != 0] # effects included in x
    cth <- x$covtheta[zeff.r, zeff.r] # their covariance matrix
    cth[is.na(cth)] <- 0 # if any effects were fixed, they have NA in covtheta
    lincomb <- sum(as.vector(t(ww)) * inftb$df["select"])
    # the desired linear combination of cell values of the influence matrix
    wt <- colSums(as.vector(t(ww)) * inftb$coeff)
    # Vector of weights for the linear combination of parameters
    names(wt) <- names(inftb$zeff.eval)
    wt.r <- wt[names(zeff.r)] # should be restricted to what is in the model
    cat("Parameter estimates and their resulting weights are \n")
    print(rbind('param'=inftb$ztheta[inftb$zeff.eval != 0], 'weight'=wt.r))
    cat("Linear combination of cells of influence matrix", round(lincomb,4),
        "\nStandard error\n")
    print(sqrt((wt.r %*% cth %*% wt.r)[1,1]))
}


##@print.influenceTable Methods
print.influenceTable <- function(x, ...)
{
	if (!inherits(x, "influenceTable"))
	{
		stop("not a legitimate influenceTable object")
	}
	if (inherits(x, "matrix"))
	{
		mat <- x
		class(mat) <- "matrix"
		print(mat)
	}
	else if (inherits(x, "data.frame"))
	{		
		mat <- matrix(as.numeric(as.character(x$select)),
            length(unique(x$zalter)),
            length(unique(x$zego)), byrow=TRUE)
		colnames(mat) <- attr(x, "levls") 
		rownames(mat) <- attr(x, "levls.alt")
		print(mat)
	}
	invisible(x)
}
	
	