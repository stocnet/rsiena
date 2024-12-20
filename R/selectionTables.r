##/*****************************************************************************
## * SIENA: Simulation Investigation for Empirical Network Analysis
## *
## * Web: https://www.stats.ox.ac.uk/~snijders/siena
## *
## * File: selectionTables.r
## *
## *
## * Description: This file contains functions for making selection tables
## *
## ****************************************************************************/

########################################################################
# R script SelectionTables.r                                           #
# written by Tom A.B. Snijders                                         #
# August 1, 2024                                                       #
########################################################################
########################################################################

########################################################################
#                                                                      #
# These are functions for constructing and presenting selection tables #
# for the interpretation of results for network dynamics               #
# obtained with the RSiena program.                                    #
# Also consult the manual, Sections 13.1 and 13.3!                     #
#                                                                      #
# It is recommended that you download this file and read it;           #
# in particular, for the functions you wish to use,                    #
# read the comments that are written at the start of the functions.    #
#                                                                      #
# The main functions to use are                                        #
# selectionTable.plot <- function(x, xd, name, vname, levls,           #
# and many other arguments, as shown below;                            #
# which makes a selection plot;                                        #
# and                                                                  #
# selectionMatrix <- function(x,xd,name,vname,levls, levls.alt=levls)  #
# which creates a matrix containing the selection table                #
# for siena or sienaGroup data set xd,                                 #
# sienaFit, sienaBayesFit, or sienaMeta object x,                      #
# actor covariate vname (should be a character string),                #
# dependent network variable name (also a character string),           #
# levels for ego levls, levels for alter levls.alt.                    #
# Mostly levls.alt will be the same as levls,                          #
# in which case it does not have to be specified.                      #
#                                                                      #
# These functions use the parameter estimates given by x$theta.        #
#                                                                      #
# You can use selectionTable.plot directly,                            #
# without first using some of the other functions.                     #
#                                                                      #
# For models utilizing creation and/or endowment effects,              #
# by specifying include.endow=TRUE                                     #
# the sum of evaluation and endowment effects is used;                 #
# by specifying include.creation=TRUE                                  #
# the sum of evaluation and creation effects is used;                  #
# the first is relevant for maintaining existing ties,                 #
# the second for creating new ties.                                    #
#                                                                      #
# The data set is only used to get the means and similarity means      #
# which are subtracted somewhere in the effects.                       #
#                                                                      #
# It is assumed variable <vname> is centered;                          #
# if it is not, some things below need to be changed.                  #
#                                                                      #
# Other functions to use (see below):                                  #
# selectionTable.se (only for sienaFit objects)                        #
# Computes standard error of a linear combination of elements          #
# of the selection table.                                              #
# This function uses the covariance matrix of parameter estimates      #
# given by x$covtheta.                                                 #
#                                                                      #
# selectionTableWithMax                                                #
# Creates a data frame of the selection table together with            #
# the maximum per ego of the selection function over alters.           #
#                                                                      #
# selectionTable.norm  (only for sienaFit objects)                     #
# Computes the location of the social norm and its standard error      #
# for the model discussed in Snijders & Lomi (2019).                   #
# This also uses the covariance matrix of parameter estimates          #
# given by x$covtheta.                                                 #
#                                                                      #
# selectionTable.plot uses ggplot2.                                    #
#                                                                      #
# If these functions are used                                          #
# for a sienaMeta object x constructed by siena08                      #
# (as opposed to a sienaFit object constructed by siena07),            #
# xd should be one of the individual data sets used for creating x,    #
# (a single-group instead of multi-group data set)                     #
# with an average ('representative') value                             #
# for the mean of variable <name>.                                     #
#                                                                      #
# If these functions are used                                          #
# for a sienaBayesFit object x  as created by sienaBayes               #
# (as opposed to a sienaFit object),                                   #
# information about the mean and similarity mean for variable vname    #
# is taken from the first data set in the list xd.                     #
########################################################################


##@selectionTable.basis  Basis for selectionTable
selectionTable.basis <- function(x, xd, name, vname,
                    levls=NULL, levls.alt=levls, nfirst=x$nwarm+1,
                    multiplier=1,
                    include.endow=FALSE, include.creation=FALSE,
                    silent=FALSE){
# Creates basic materials for a selection table for
# siena data set xd, sienaFit object x,
# actor variable vname (should be a character string),
# dependent variable name (also a character string),
# levels for ego levls, levels for alter levls.alt.
# The multiplier is used in case
# the variable vname has a different natural scale,
# and the values for ego and alter as reported should be multiplied.
# The values given for levls and levls.alt are before this multiplication,
# i.e., they are the values as in the data set.
#
# The function can be used also for a sienaMeta or sienaBayesFit object x.
# Then xd must be one of the data sets used for creating x,
# with an average ('representative') value for the mean of variable <vname>.
# nfirst is used only for sienaBayesFit objects x.
#
    dsign <- function(d){
        0.5*(ifelse(d > 0, 1, 0) + ifelse(d >= 0, 1, 0))
    }
	if (!silent)
    {
		cat("Dependent network",name, "; actor variable",vname,".\n")
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
# Obtain means
    depvar <- FALSE
    if (inherits(xd, "sienaGroup"))
    {
        cat("A sienaGroup data object was given.\n")
        # Does the variable exist, and is it a dependent variable?
        xd1 <- xd[[1]]
        thevar <- xd1$cCovars[[vname]]
        if (is.null(thevar)) {thevar <- xd1$vCovars[[vname]]}
        if (is.null(thevar))
            {
                depvar <- TRUE
                thevar <- xd1$depvars[[vname]]
            }
        if (is.null(thevar)){stop(paste('There is no actor variable <',vname,'>.'))}
        if (depvar)
        {
#           vmean <- mean(sapply(xd, function(z){mean(z$depvars[[vname]], na.rm=TRUE)}))
            vtot <- sum(sapply(xd, function(z){sum(z$depvars[[vname]], na.rm=TRUE)}))
            vtotN <- sum(sapply(xd, function(z){sum(!is.na(z$depvars[[vname]]))}))
            vmean <- vtot/vtotN
            vsmean <- attr(xd, "bSim")[[vname]]
        }
        else
        { # constant covariates are transformed by sienaGroupCreate to varying covariates
            vmean <- attr(xd, "vCovarMean")[vname]
            vsmean <- attr(xd, "vCovarSim")[vname]
        }
    }
    else
    {
        thevar <- xd$cCovars[[vname]]
        if (is.null(thevar)) {thevar <- xd$vCovars[[vname]]}
        if (is.null(thevar))
            {
                depvar <- TRUE
                thevar <- xd$depvars[[vname]]
            }
        if (is.null(thevar)){stop(paste('There is no actor variable <',vname,'>.'))}
        if (depvar)
        {# Then the mean is not stored as an attribute
            means <- colMeans(thevar, na.rm=TRUE)
            vmean <- mean(means)
        }
        else
        {
            vmean <- attr(thevar, 'mean')
        }
        vsmean <- attr(thevar, 'simMean')
    }
    if (!depvar)
    {
        if (!attr(thevar, 'centered'))
        {
            cat('Warning: this variable is not centered. The mean is subtracted anyway; this may be wrong.\n')
        }
    }
    Delta  <- attr(thevar, 'range')
    Delta12  <- attr(thevar, 'range2')
    if (is.null(levls))
    {
        levls <- floor(Delta12[1]) : ceiling(Delta12[2])
    }
    if (is.null(levls.alt))
    {
        levls.alt <- levls
    }
    replace0 <- function(k){ifelse(length(k)==0,0,k)}
# If additional efNames are going to be included,
# adapt also the coeffs matrix below and function selectionTableWithMax.
    efNames <- c('altX','altSqX','egoX','egoSqX','egoXaltX',
                    'simX','diffX','diffSqX','higher','sameX','egoDiffX','egoPlusAltX')
    veff.eval <- sapply(efNames, function(s){replace0(which(theEffects$name == name &
                    theEffects$interaction1 == vname &
                    (theEffects$type == 'eval') &
                    theEffects$include &
                    theEffects$shortName==s))})
    veff.endow <- sapply(efNames, function(s){replace0(which(theEffects$name == name &
                    theEffects$interaction1 == vname &
                    (theEffects$type == 'endow') &
                    theEffects$include &
                    theEffects$shortName==s))})
    veff.creation <- sapply(efNames, function(s){replace0(which(theEffects$name == name &
                    theEffects$interaction1 == vname &
                    (theEffects$type == 'creation') &
                    theEffects$include &
                    theEffects$shortName==s))})

# veff gives the indicators of the effects in x
    taketheta <- function(k){ifelse(k==0,0,theta[k])}
    vtheta <- sapply(veff.eval, taketheta)
    if (include.endow)
    {
        vtheta <- vtheta + sapply(veff.endow,taketheta)
    }
    if (include.creation)
    {
        vtheta <- vtheta + sapply(veff.creation,taketheta)
    }
    if (all(vtheta == 0))
    {
        cat('All parameters found for the effect of', vname, 'on', name, 'are 0.\n')
        stop('there seems to be something wrong with x, xd, name, or vname.')
    }

    if (!silent){
        cat('Parameters found are\n')
        print(vtheta[vtheta != 0.0])
        cat('Check that these parameter values are correct!\n')
        if (include.endow)
        {
            cat('This includes evaluation and endowment (maintenance) effects.\n')
        }
        if (include.creation)
        {
            cat('This includes evaluation and creation effects.\n')
        }
    }
    if (!silent)
    {
        cat('The subtracted mean value of', vname, 'is', vmean,'.\n')
        if (vtheta["simX"] != 0)
		{
			cat('The subtracted similarity mean value of', vname, 'is', vsmean,
															'.\n')
		}
    }
    flush.console()
# vtheta contains the parameter values in x
    K <- length(levls)
    KA <- length(levls.alt)
    valter <- rep(levls.alt,K)
    vego <- rep(levls,each=KA)
# coeffs duplicates the information in f6 in selectionTableWithMax.
# This duplication is undesirable, but I do this anyway.
# When maintaining this script, coeffs and f6 in selectionTable.WithMax
# must stay in line.
    coeffs <- matrix(NA, K*KA, length(efNames))
    coeffs[,1] <- (valter - vmean)
    coeffs[,2] <- (valter - vmean)*(valter - vmean)
    coeffs[,3] <- (vego - vmean)
    coeffs[,4] <- (vego - vmean)*(vego - vmean)
    coeffs[,5] <- (vego - vmean)*(valter - vmean)
    coeffs[,6] <- (1-(abs(valter - vego)/Delta)-vsmean)
    coeffs[,7] <- (valter - vego)
    coeffs[,8] <- (valter - vego)*(valter - vego)
    coeffs[,9] <- dsign(vego - valter)
    coeffs[,10] <- 1*(vego == valter)
    coeffs[,11] <- vego*(valter - vego)
    coeffs[,12] <- (vego - vmean) + (valter - vmean)
    select <- coeffs %*% vtheta
    fact <- multiplier*(1:K)
    ego <- factor(rep(fact,each=KA), ordered=TRUE)
    vego <- multiplier*vego
    valter <- multiplier*valter
    df <- data.frame(ego,vego,valter,select)
    list(df=df, veff.eval=veff.eval, veff.endow=veff.endow,
         veff.creation=veff.creation, vtheta=vtheta, vmean=vmean, vsmean=vsmean,
         Delta=Delta, coeffs=coeffs, levls=levls, levls.alt=levls.alt)
}


##@selectionTable  selectionTable
selectionTable <- function(x, xd, name, vname,
                    as.matrix=FALSE,
                    levls=NULL, levls.alt=levls, nfirst=x$nwarm+1,
                    multiplier=1,
					include.endow=FALSE, include.creation=FALSE,
					silent=FALSE){
# Creates a data frame containing the selection table for
# siena data set xd, sienaFit object x,
# actor covariate vname (should be a character string),
# dependent variable name (also a character string),
# levels for ego levls, levels for alter levls.alt.
    sel.t <- selectionTable.basis(x=x, xd=xd, name=name, vname=vname,
				levls=levls, levls.alt=levls.alt, nfirst=nfirst,
                multiplier=multiplier,
				include.endow=include.endow,
				include.creation=include.creation,
				silent=silent)
	df <- sel.t$df
    vals <- as.character(unique(df$vego))
    df$ego <- factor(as.character(df$vego), levels=vals, ordered=TRUE)
    df$valter <- as.numeric(as.character(df$valter))
    df$select <- as.numeric(as.character(df$select))
	class(df) <- c("selectionTable", class(df))
	attr(df, "name") <- name
	attr(df, "vname") <- vname
	attr(df, "multiplier") <- multiplier
	attr(df, "levls") <- sel.t$levls
	attr(df, "levls.alt") <- sel.t$levls.alt
	if (as.matrix)
	{
		mat <- matrix(as.numeric(as.character(df$select)),
            length(unique(df$vego)),
            length(unique(df$valter)), byrow=TRUE)
		colnames(mat) <- sel.t$levls.alt
		rownames(mat) <- sel.t$levls
		class(mat) <- c("selectionTable" ,class(mat))
		return(mat)
	}
	else
	{
		return(df)
	}
}

selectionTableWithMax <- function(x, xd, name, vname,
                levls=NULL, levls.alt=levls, nfirst=x$nwarm+1,
                multiplier=1, discrete=TRUE, 
                include.endow=FALSE, include.creation=FALSE,
                silent=FALSE){
# Creates a data frame including the selection table for
# siena data set xd, sienaFit object x,
# actor covariate vname (should be a character string),
# dependent variable name (also a character string),
# levels for ego levls, levels for alter levls.alt,
# to which is appended the values for the maximum across alter,
# where alter ranges over levls.alt
# (a grid with 200 points between minimum and maximum is used).
# The variable named 'kind' is 1 for the selection table and 2 for the maximum.
    dsign <- function(d){
        0.5*(ifelse(d > 0, 1, 0) + ifelse(d >= 0, 1, 0))
    }
    st <- selectionTable.basis(x, xd, name, vname, levls, levls.alt, nfirst=nfirst,
                                    multiplier,
				include.endow=include.endow,
				include.creation=include.creation,
				silent=silent)
    df1 <- st$df
    df1$kind <- rep(1,dim(df1)[1])
    vtheta <- st$vtheta
    vmean <- st$vmean
    Delta <- st$Delta
    vsmean <- st$vsmean
# f6 duplicates the information in coeffs in selectionTable.basis.
# This duplication is undesirable as programming style, but I do this anyway.
# When maintaining this script, f6 and coeffs in selectionTable.basis
# must stay in line.
    f6 <- function(ve,va){
        contr <- vtheta[1]*(va-vmean) + vtheta[2]*(va-vmean)*(va-vmean) +
            vtheta[3]*(ve-vmean) +
            vtheta[4]*(ve-vmean)*(ve-vmean) + vtheta[5]*(ve-vmean)*(va-vmean) +
            vtheta[6]*(1-(abs(va-ve)/Delta)-vsmean) +
            vtheta[7]*(va-ve) + vtheta[8]*(va-ve)*(va-ve) +
            vtheta[9]*dsign(ve-va) +
            vtheta[10]*(1*(ve == va)) +
            vtheta[11]*ve*(va-ve) +
            vtheta[12]*((va-ve) + (va-ve))
        names(contr) <- 'contribution'
        contr}
# vtheta contains the parameter values in x
    K <- length(levls)
    KA <- length(levls.alt)
    valter <- rep(levls.alt,K)
    vego <- rep(levls,each=KA)
    fact <- 1:K
    ego <- factor(rep(fact,each=KA))
# Now calculate the maximum.
    minv <- min(levls.alt)
    maxv <- max(levls.alt)
    gridv <- minv + (0:200)*((maxv-minv)/200)
# vsmean is not important because it is an additive constant
    if (discrete){
        altm <- as.integer(round(minv)):as.integer(round(maxv))
        egom <- rep(levls[1], length(altm)) # dummy
        maxm <- sapply(altm, function(x){max(f6(x, gridv))})
    } else {
        egom <- rep(NA, K*4)
        altm <- rep(NA, K*4)
        maxm <- rep(NA, K*4)
        for (i in 1:K){
            for (j in 1:4){
            egom[4*(i-1) + j] <- levls[i]
            x1 <- ifelse(i <= 1, levls[1], 0.5*(levls[i-1]+levls[i]))
            x2 <- ifelse(i >= K, levls[K], 0.5*(levls[i+1]+levls[i]))
            dd <- x2-x1
            alterij <- ((j-1)/3)*x2 + (1 - (j-1)/3)*x1
            altm[4*(i-1) + j] <- alterij
#           maxm[4*(i-1) + j]  <- max(f6(alterij, gridv))
            maxm[4*(i-1) + j]  <- max(f6(levls[i], gridv))
            }
        }
        }
    df1$ego <- as.character(df1$vego)
    df2 <- data.frame(ego=egom,vego=egom,valter=altm,select=maxm,kind=2)
    df <- rbind(df1,df2)
    df$ego <- as.character(df$vego)
#   df$valter <- as.numeric(as.character(df$valter)) # superfluous
#   df$select <- as.numeric(as.character(df$select)) # superfluous
    df
}

selectionTable.se <- function(x, xd, name, vname,
                    levls=NULL, ww, levls.alt=levls, nfirst=x$nwarm+1,
                    multiplier=1){
# Calculates standard errors for a linear combination
# of elements of the selection table for siena data set xd, sienaFit object x,
# actor covariate vname (should be a character string),
# dependent variable name (also a character string),
# levels for ego levls, levels for alter levls.alt,
# Coefficients of the linear combination are in matrix ww,
# which is assumed to have rows corresponding to levls and columns to levls.alt.
# The linear combination for which the standard error is computed is
# sum_{h,k} ww[h,k] * selectionTable[h,k].
    if (inherits(x, "sienaBayesFit")){
        stop('This function does not work for sienaBayesFit objects.\n')
    }
    if (!all(dim(ww) == c(length(levls), length(levls.alt)))){
        stop('Dimension of ww should be the lengths of levls and levls.alt.\n')
    }
    cat("Requested cell weights (row = ego; col = alter):\n")
    print(cbind(which(ww != 0, arr.ind=TRUE),
            value=ww[which(ww != 0, arr.ind=TRUE)]))
    sb <- selectionTable.basis(x, xd, name, vname, levls, levls.alt,
                nfirst=nfirst, multiplier, silent=TRUE)
    veff.r <- sb$veff.eval[sb$veff.eval != 0] # effects included in x
    cth <- x$covtheta[veff.r, veff.r] # their covariance matrix
    cth[is.na(cth)] <- 0 # if any effects were fixed, they have NA in covtheta
    lincomb <- sum(as.vector(t(ww)) * sb$df["select"])
    # the desired linear combination of cell values of the selection matrix
    wt <- colSums(as.vector(t(ww)) * sb$coeff)
    # Vector of weights for the linear combination of parameters
    names(wt) <- names(sb$veff.eval)
    wt.r <- wt[names(veff.r)] # should be restricted to what is in the model
    cat("Parameter estimates and their resulting weights are \n")
    print(rbind('param'=sb$vtheta[sb$veff.eval != 0], 'weight'=wt.r))
    cat("Linear combination of cells of selection matrix", round(lincomb,4),
        "\nStandard error\n")
    print(sqrt((wt.r %*% cth %*% wt.r)[1,1]))
}

selectionTable.norm <- function(x, xd, name=attr(xd$depvars,"name")[1], vname,
                                    nfirst=x$nwarm+1, multiplier=1){
# Calculates the location of the social norm and its standard error
# in the attraction function for siena data set xd, sienaFit object x,
# actor covariate vname (should be a character string),
# dependent variable name (also a character string),
# if the attraction function with terms
# 'altX','altSqX','egoX','egoSqX','diffSqX' is used.
# Note that the attraction function should not include the 'egoXaltX' effect!
    stab <- selectionTable.basis(x,xd,name,vname,1:2, nfirst=nfirst, multiplier)
# the 1:2 is arbitrary, not used
    vtheta <- stab$vtheta
    vmean <- stab$vmean
    veff <- stab$veff
	veff.eval <- stab$veff.eval
    if (vtheta["egoXaltX"] != 0){
        cat('Warning: effect of egoXaltX is ',vtheta["egoXaltX"],
             ', not equal to 0; \n')
        cat(' this function does not apply to such a sienaFit object.\n')
        }
# calculate gradient
    grad <- matrix(0, length(x$theta), 1)
    grad[veff.eval["altX"],1] <- -1/(2*vtheta["altSqX"])
    grad[veff.eval["altSqX"],1] <- vtheta["altX"]/(2*vtheta["altSqX"]*vtheta["altSqX"])
    covtheta <- x$covtheta
    covtheta[is.na(covtheta)] <- 0
    if (vtheta["altSqX"] > 0){
        cat('Warning: the coefficient of alter squared is positive.\n')
        cat('This means the extremum of the attraction function is a minimum,\n')
        cat('and cannot be interpreted as a social norm.\n')
        flush.console()
        }
    list(vnorm = vmean - (vtheta["altX"]/(2*vtheta["altSqX"])),
        se.vnorm = sqrt(t(grad) %*% covtheta %*% grad))
}


##@print.selectionTable Methods
print.selectionTable <- function(x, ...)
{
	if (!inherits(x, "selectionTable"))
	{
		stop("not a legitimate selectionTable object")
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
            length(unique(x$vego)),
            length(unique(x$valter)), byrow=TRUE)
		colnames(mat) <- attr(x, "levls.alt")
		rownames(mat) <- attr(x, "levls") 
		print(mat)
	}
	invisible(x)
}
	
	