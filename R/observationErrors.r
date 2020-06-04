#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: observationErrors.r
# *
# * Description: This module contains various functions for parameter estimation
# * in models with observation errors
# * by maximum likelihood under normality,
# * extending Snijders & Baerveldt (2003), for use in siena08.
# *****************************************************************************/
##args:x: vector of observations
##     se: vector of same length as x: standard errors of observations
## returns updated z

##@deviance.observationErrors siena08: deviance under normality assumptions
deviance.observationErrors <- function(mu, sig, x, se)
{
    ## deviance.observationErrors: given vector observations x, assuming model
    ## x ~ normal(expected = mu, variance = sig^2 + se^2)
    ## required is length(x) = length(se)
    sig2 <- sig^2
    sum((x - mu)^2/(sig2 + se^2)) + sum(log(sig2 + se^2))
}

##@profdev.mu siena08 (profile deviance =) minus twice profile loglik for mu
profdev.mu <- function(mu,x,se)
{
    ## deviance.observationErrors maximized over sigma
    ra <- max(x) - min(x)  # easy upper bound for sigma
    optimize(function(sig)deviance.observationErrors(mu, sig, x, se),
             c(0, ra))
}


##@profdev.sig siena08 (profile deviance =) minus twice prof. loglik for sigma
profdev.sig <- function(sig,x,se)
{
    ## deviance.observationErrors maximized over mu
    optimize(function(mu)deviance.observationErrors(mu, sig, x, se),
             range(x))
}

##@maxlik siena08 maximum likelihood estimator
maxlik <- function(x,se)
{
    if (length(x) > 1)
    {
        ## deviance.observationErrors minimized over mu and sigma
        ## Minimize profile deviance for mu:
        opmu  <- optimize(function(mu)profdev.mu(mu, x, se)$objective,
                          range(x))
        ## MLE for mu:
        mu    <- opmu$minimum
        ## minimized deviance:
        dev   <- opmu$objective
        ## Location for minimum of deviance for this value of mu:
        sig   <- profdev.mu(mu, x , se)$minimum
        ## Standard error of MLE(mu):
        se.mu <- sqrt(1 / sum(1 / (se^2 + sig^2)))
    }
    else
    {
        mu <- x
        se.mu <- NA
        sig <- NA
        dev <- NA
    }
    return(list(mu = mu, se.mu = se.mu, sigma = sig, deviance = dev))
}

##@unisroot siena08 root finder when interval may be inadequate
unisroot <- function(f, interval, ..., left=TRUE)
{
    ## tries to solve f(x) = 0,
    ## first within interval,
    ## if endpoints do not have opposite signs
    ## then in interval extended to left or to right
    ## Return value is a list which includes root = location of the root,
    ## and f.root = function value at this location.
    if (f(interval[1]) * f(interval[2]) < 0)
    {
        uniroot(f, interval, ...)
    }
    else
    {
        x1 <- interval[1]
        x2 <- interval[2]
        ra <- x2 - x1
        for (it in (1:1000))
        {
            if (left)
            {
                x1 <- x1 - ra
            }
            else
            {
                x2 <- x2 + ra
            }
            if (f(x1) * f(x2) < 0)
            {
                break
            }
        }
        if (f(x1) * f(x2) < 0)
        {
            uniroot(f, c(x1, x2), ...)
        }
        else
        {
            list(root=NA, f.root=NA)
        }
    }
}

##@confint.mu siena08 confidence interval for mu
confint.mu <- function(x, se, alpha=0.05)
{
    ## confidence interval for mu in model
    ## x ~ normal(expected = mu, variance = sigma^2 + se^2)
    ## with length(x) = length(se).
    ## returns list consisting of bounds confidence interval and
    ## confidence level (1-alpha).
    ma <- max(x)
    mi <- min(x)
    maxli <- maxlik(x, se)        # ML estimators
    mindev <- maxli$deviance       # minimized deviance
    mlemu <- maxli$mu             # MLE for mu
    chidev <- qchisq(1 - alpha, 1) # critical value chi-squared distribution
    if (length(x) > 1)
    {
        ## The end points of the confidence interval are the values
        ## where the profile deviance for mu is equal to mindev + chidev.
        ## Now first the solution for the left side of the confidence interval:
        tmp1 <- unisroot(function(mu)
                     {
                         tmp <- profdev.mu(mu, x, se)
                         tmp$objective - mindev - chidev
                     },
                         c(mi, mlemu))
        mu1 <- tmp1$root
        ## and then the solution for the right side of the confidence interval:
        tmp1 <- unisroot(function(mu)
                     {
                         tmp <- profdev.mu(mu, x, se)
                         tmp$objective - mindev - chidev
                     },
                         c(mlemu, ma), left=FALSE)
        mu2 <- tmp1$root
    }
    else
    {
        mu1 <- NA
        mu2 <- NA
    }
    c(mu1, mu2, 1 - alpha)
}

##@confint.sig siena08 confidence interval for sigma
confint.sig <- function(x, se, alpha=0.05)
{
    ## confidence interval for sigma in model
    ## x ~ normal(expected = mu, variance = sigma^2 + se^2)
    ## with length(x) = length(se)
    ## returns list consisting of bounds confidence interval and
    ## confidence level (1 - alpha).     ma <- max(x)
    ra <- max(x) - min(x)          # easy upper bound for sigma
    maxli  <- maxlik(x, se)        # ML estimators
    mindev <- maxli$deviance       # minimized deviance
    mlesig <- maxli$sigma          # MLE for sigma
    chidev <- qchisq(1 - alpha, 1) # critical value chi-squared distribution
    ## The end points of the confidence interval are the values
    ## where the profile deviance for sigma is equal to mindev + chidev,
    ## unless the profile deviance in 0 is smaller than this.
    ## Now first the solution for the left side of the confidence interval;
    ## this may be 0.
    if (length(x) > 1)
    {
        if (profdev.sig(0, x, se)$objective <= mindev + chidev)
        {
            sig1 <- 0
        }
        else
        {
            sig1 <- uniroot(function(sig)
                            profdev.sig(sig, x, se)$objective - mindev - chidev,
                            c(0, mlesig))$root
        }
        ## and then the solution for the right side of the confidence interval:
        sig2 <- unisroot(function(sig)
                         profdev.sig(sig, x, se)$objective - mindev - chidev,
                         c(mlesig, ra), left=FALSE)$root
    }
    else
    {
        sig1 <- NA
        sig2 <- NA
    }
    return(c(sig1, sig2, 1 - alpha))
}


