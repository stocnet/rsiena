\name{simstats0c}
\alias{simstats0c}
\alias{maxlikec}
\alias{initializeFRAN}
\alias{terminateFRAN}
\title{Versions of FRAN}
\description{
 The functions to be called as "FRAN" by \code{\link{siena07}}.   
 They call compiled C++. Not for general users' use.
}
\usage{
simstats0c(z, x, data=NULL, effects=NULL, fromFiniteDiff=FALSE,
           returnDeps=FALSE, returnChains=FALSE, returnChangeContributions=FALSE,
           byWave=FALSE, returnDataFrame=FALSE, returnLoglik=FALSE)
maxlikec(z, x, data=NULL, effects=NULL,
           returnChains=FALSE, byGroup = FALSE, 
           byWave=FALSE, returnDataFrame=FALSE, returnLoglik=FALSE, onlyLoglik=FALSE)
initializeFRAN(z, x, data, effects, prevAns=NULL, initC,
           profileData=FALSE, returnDeps=FALSE, returnChains=FALSE,
           returnChangeContributions=FALSE, byGroup=FALSE, returnDataFrame=FALSE,
           byWave=FALSE, returnLoglik=FALSE, onlyLoglik=FALSE)
terminateFRAN(z, x)
}
\arguments{
  \item{z}{Control object, passed in automatically in  \code{\link{siena07}}.}
  \item{x}{A sienaAlgorithm object, passed in automatically in  \code{\link{siena07}}.}
  \item{data}{A sienaData object as returned by \code{\link{sienaDataCreate}}.}
  \item{effects}{A sienaEffects object as
    returned by \code{\link{getEffects}}.}
  \item{fromFiniteDiff}{Boolean used during calculation of derivatives
    by finite differences. Not for user use.}
  \item{returnDeps}{Boolean. Whether to return the simulated networks in
    Phase 3.}
  \item{returnChains}{Boolean. Whether to return the chains.}
  \item{returnChangeContributions}{Boolean. Whether to return the change contributions.}
  \item{byWave}{Boolean. Whether to return the finite difference or
    maximum likelihood
    derivatives by wave (uses a great deal of memory). Only necessary
    for \code{\link{sienaTimeTest}}}
  \item{byGroup}{Boolean. For internal use: allows different thetas for each
    group to be used in \code{sienaBayes}.}
  \item{returnDataFrame}{Boolean. Whether to return the chains as
    lists or data frames.}
  \item{returnLoglik}{Boolean. Whether to return the log likelihood of
    the simulated chain.}
   \item{onlyLoglik}{Boolean: whether to return just the likelihood for
    the simulated chain, plus details of steps accepted and rejected. }
  \item{prevAns}{An object of class "sienaFit" as returned by
   \code{\link{siena07}}, from which scaling
    information (derivative matrix and standard deviation of the
    deviations) will be extracted along with the latest version of the
    parameters which will be used as the initial values, unless the
    model requests the use of standard initial values. If the previous
    model is exactly the same as the current one, Phase 1 will be omitted. If
    not, any parameter estimates for effects which are included in
    the new model will be used as initial
    values, but phase 1 will still be carried out.
    If the results used as \code{prevAns} are a reasonable starting point,
    this will increase the efficiency of the algorithm. }
  \item{initC}{If TRUE, call is to setup the data and model in C++. For
    use with multiple processes only.}
 \item{profileData}{Boolean to force dumping of the data for profiling
    with \command{sienaProfile.exe}.}
}
\details{
  Not for general users' use.\cr
  The name of \code{simstats0c} or \code{maxlikec}
  should be used for the element FRAN of the model object, the former when using
  estimation by forward simulation, the latter for maximum likelihood estimation.
  The arguments with no defaults must be passed in on the call to
  \code{\link{siena07}}.
  \code{initializeFRAN} and \code{terminateFRAN} are called in both cases.
}
\value{\code{simstats0c} returns a list containing:
  \item{fra}{Simulated statistics.}
  \item{sc}{Scores with which to calculate the derivative (not phase 2
    or if using finite differences or maximum likelihood).}
  \item{dff}{Contributions to the derivative if finite differences}
  \item{ntim}{For conditional processing, time taken.}
  \item{feasible}{Currently set to TRUE.}
  \item{OK}{Could be set to FALSE if serious error has occurred.}
  \item{sims}{A list of simulation results, one for each period. Each
    list consists of a list for each data object, each of which consists of
    a list for each network, each of which consists of a list for each period,
    each component of which is an edgelist in matrix form (the columns
    are from, to, value) (or vector for behavior variables).
    Only if \code{returnDeps} is TRUE.}
  \code{maxlikec} returns a list containing:
  \item{fra}{Simulated scores.}
   \item{dff}{Simulated Hessians: stored as lower triangular matrices}
  \item{ntim}{NULL, compatibility only}
  \item{feasible}{Currently set to TRUE.}
  \item{OK}{Could be set to FALSE if serious error has occurred.}
  \item{dff}{Simulated Hessian}
  \item{sims}{NULL, for compatibility only}
  \item{chain}{A list of sampled chains, one for each period. Each
    list consists of a list for each data object, each of which consists of
    a list for each network, each of which consists of a list for each period,
    each component of which is a list or a data frame depending on the
    value of \code{returnDataFrame}.
    Only if \code{returnChains} is TRUE.}
  \item{changeContributions}{A list of sampled chains, one for
    each period. Each list consists a list for each dependent variable, each 
    of which consists of a list for each period, each component of which 
    is a matrix of effects times choices containing the changeContributions 
    (+- change statistics) for each ministep.
    Only if \code{returnChangeContributions} is TRUE.}
  \item{accepts}{Number of accepted MH steps by dependent variable
    (permute steps are counted under first dependent variable)}
  \item{rejects}{Number of rejected MH steps by dependent variable
    (permute steps are counted under first dependent variable)}
  \item{aborts}{Number of aborted MH steps counted under first dependent
    variable.}
  \item{loglik}{Loglikelihood of the simulations. Only if
    \code{returnLoglik} is TRUE. If \code{onlyLoglik} is TRUE, only
    \code{loglik}, \code{accepts},\code{rejects} and \code{aborts} are
    returned.}
  \code{initializeFRAN} and \code{terminateFRAN} return the control
    object \code{z}.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }
\seealso{\code{\link{siena07}} }
\examples{

mynet1 <- sienaNet(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip)
myalgorithm <- sienaAlgorithmCreate(fn=simstats0c, nsub=2, n3=100, projname=NULL)
ans <- siena07(myalgorithm, data=mydata, effects=myeff,  batch=TRUE)
}
\keyword{models}
