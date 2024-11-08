\name{siena07}
\alias{siena07}
\concept{estimation}
\title{Function to estimate parameters in a Siena model}
\description{
 Estimates parameters in a Siena model using Method of Moments,
 based on direct simulation, conditional or otherwise;
 or using Generalized Method of Moments;
 or using Maximum Likelihood by MCMC simulation.
 Estimation is done using a Robbins-Monro algorithm.
 Note that the data and particular model to be used
 must be passed in using named arguments as the \code{...},
 and the specification for the algorithm must be passed on
 as \code{x},
 which is a \code{\link{sienaAlgorithm}} object as produced by
 \code{\link{sienaAlgorithmCreate}} (see examples).
}
\usage{
siena07(x, batch=FALSE, verbose=FALSE, silent=FALSE,
        useCluster=FALSE, nbrNodes=2,
        thetaValues = NULL,
        returnThetas = FALSE,
        thetaBound = 50,
        targets = NULL,
        initC=TRUE,
        clusterString=rep("localhost", nbrNodes), tt=NULL,
        parallelTesting=FALSE, clusterIter=!x$maxlike,
        clusterType=c("PSOCK", "FORK"), cl=NULL, ...)
      }
\arguments{
  \item{x}{A control object, of class \code{\link{sienaAlgorithm}. }}
  \item{batch}{ Desired interface: \code{FALSE} gives a gui
     (graphical user interface implemented as a tcl/tk screen),
     \code{TRUE} gives a small (if \code{verbose=FALSE}) amount of printout
     to the console.}
  \item{verbose}{Produces various output to the console if \code{TRUE}.}
  \item{silent}{Produces no output to the console if \code{TRUE},
     even if batch mode.}
  \item{useCluster}{Boolean: whether to use a cluster of processes (useful if
  multiple processors are available).}
  \item{nbrNodes}{Number of processes to use if useCluster is \code{TRUE}.}
  \item{thetaValues}{If not \code{NULL}, this should be a matrix with parameter
    values to be used in Phase 3. The number of columns must be equal to the
    number of estimated parameters in the effects object (if conditional
    estimation is used, without the rate parameters for the conditioning
    dependent variable). Can only be used if \code{x$simOnly=TRUE}.}
  \item{returnThetas}{Boolean: whether to return theta values
   and generated estimation statistics of Phase 2 runs.}
  \item{thetaBound}{Numeric: if at some moment during estimation the maximum
   absolute value of the parameters exceeds \code{thetaBound}, estimation
   is interrupted. In interactive use, the user is requested to give
   a higher number; if non-numeric input is given, estimation stops.
   In non-interactive use, estimation stops anyway.
   The check is turned off by using \code{thetaBound=Inf}.}
  \item{targets}{Numeric vector of length equal to the number of estimated
    parameters, meant to supersede the targets calculated from the data set;
    see "Details". Not for regular use.}
  \item{initC}{Boolean: set to \code{TRUE} if the simulation will use C
   routines (currently always needed). Only for use if using
   multiple processors, to ensure all copies are initialised
   correctly. Ignored otherwise, so is set to \code{TRUE} by default.}
  \item{clusterString}{Definitions of clusters. Default set up to use
   the local machine only.}
  \item{tt}{A \code{tcltk} toplevel window. Used if called from the
   model options screen, if \code{tcltk} is available.}
  \item{parallelTesting}{Boolean. If \code{TRUE}, sets up random numbers to
   parallel those in Siena 3.}
  \item{clusterIter}{Boolean. If \code{TRUE}, multiple processes execute
   complete iterations at each call.
   If \code{FALSE}, multiple processes execute a single wave at each call.}
  \item{clusterType}{Either "PSOCK" or "FORK". On Windows, must be
  "PSOCK". On a single non-Windows machine may be "FORK", and
  subprocesses will be formed by forking. If "PSOCK", subprocesses are
  formed using \R scripts.}
  \item{cl}{An object of class c("SOCKcluster", "cluster") (see Details).}
  \item{\dots}{Arguments for the simulation function, see
    \code{\link{simstats0c}}:
    in any case, \code{data} and \code{effects}, as in the examples below;\cr
    possibly also \code{prevAns} if a previous reasonable
    provisional estimate was obtained for a similar model;\cr
    possibly also \code{returnDeps} if the simulated dependent variables
    (networks, behaviour) should be returned;\cr
    possibly also \code{returnChains} if the simulated sequences (chains)
    of ministeps should be returned; this may produce a very big file.}
}
\details{
  This is the main function and workhorse of RSiena.

  For use of \code{siena07}, it is necessary to specify parameters \code{data}
  (RSiena data set) and \code{effects} (effects object), which are
  required parameters in function \code{\link{simstats0c}}.
  (These parameters are inserted through '\dots'.) See the examples.

 \code{siena07} runs a Robbins-Monro algorithm for parameter estimation
  using the three-phase implementation described in Snijders (2001, 2017),
 with (if \code{x$findiff=FALSE})
 derivative estimation as in Schweinberger and Snijders (2007).
 The default is estimation according to the Method of Moments
 as in Snijders, Steglich and Schweinberger (2007).\cr
 If  \code{x$gmm=TRUE} and \code{myeff} contains one or more
 \code{gmm} statistics as included by function
 \code{\link{includeGMoMStatistics}}, the algorithm employs the
 Generalized Method of Moments as defined in Amati, Schoenenberger, and
 Snijders (2015, 2019).\cr
 For continuous behavior variables defined with \code{type="continuous"} in
 \code{\link{sienaDependent}}, estimation is done as
 described in Niezink and Snijders (2017).\cr
 If  \code{x$maxlike=TRUE}, estimation is done by Maximum Likelihood
 implemented as in Snijders, Koskinen and Schweinberger (2010).\cr
 Phase 1 does a few iterations to
 estimate the derivative matrix of the targets with respect to the
 parameter vector. Phase 2 does the estimation. Phase 3 runs a
 simulation to estimate standard errors and check convergence of the model.
 The simulation function is called once for each iteration in these phases
 and also once to initialise the model fitting and once to complete it.
 Unless in batch mode, a tcl/tk screen is displayed to allow interruption
 and to show progress.

 If \code{targets} is specified (which should be done only in special cases),
 and provided that estimation is by the Method of Moments, the data
 is not a multi-group data set and has exactly 2 waves,
 and if the length of the vector \code{targets} is equal to the number
 of estimated parameters (not counting the rate parameters estimated
 by conditional estimation), then the vector \code{targets} supersedes the
 targets calculated from the data set.

 It is necessary to check that convergence has been achieved.
 The rule of thumb is that the all t-ratios for convergence should be
 in absolute value less than 0.1 and
 the overall maximum convergence ratio should be less than 0.25.
 If this was not achieved, the result can be used to start
 another estimation run from the estimate obtained, using
 the parameter \code{prevAns} as illustrated in the example below.
 (This parameter is inserted through '\dots' into the function
 \code{\link{initializeFRAN}}.)

 For good estimation of standard errors, it is necessary that \code{x$n3} is
 large enough. More about this is in the manual. The default value \code{x$n3}
 set in \code{\link{sienaAlgorithmCreate}} is adequate for most explorative use,
 but for presentation in publications larger values are necessary, depending
 on the data set and model; e.g., \code{x$n3=3000} or larger.

 Parameters can be tested against zero by dividing the estimate by its
 standard error and using an approximate standard normal null distribution.
 Further, functions \code{\link{Wald.RSiena}} and
 \code{\link{Multipar.RSiena}} are available for multi-parameter testing.\cr
 Parameters specified in \code{\link{includeEffects}} or
 \code{\link{setEffect}} with \code{fix=TRUE, test=TRUE} will not be estimated;
 score tests of their hypothesized values are reported in the output file
 specified in the control (algorithm) object.
 These tests can be obtained also using \code{\link{score.Test}}.

 If \code{x$simOnly} is \code{TRUE}, which is meant to go together with
 \code{x$nsub=0}, the calculation of the standard errors and covariance
 matrix at the end of Pase 3 is skipped. No estimation is performed.
 If \code{thetaValues} is not \code{NULL}, the parameter values in
 the rows of this matrix will be used in the consecutive runs of Phase 3.
 If \code{x$n3} is larger than the number of rows times \code{nbrNodes}
 (see below), the last row of \code{thetaValues} will continue to be used.
 The parameter values actually used will be stored in the
 output matrix \code{thetaUsed}.

 Multiple processors are used for estimation by MoM to distribute each iteration
 in each subphase over the cluster of nodes. The number of iterations
 accordingly will be divided (approximately) by the number of nodes;
 for phase 2, unless \code{n2start} is specified.
 This implies that if multiple processors are used, think of dividing
 \code{n2start} by \code{nbrNodes}.\cr
 For estimation by ML, multiple processing is done per period. Therefore,
 for one period (two waves) and one group, this will have no effect.

 In the case of using multiple processors, there are two options for telling
 \code{siena07} to use them. By specifying the options \code{useCluster},
 \code{nbrNodes}, \code{clusterString} and \code{initC}, \code{siena07}
 will create a \code{cluster object} that will be used by the
 \code{parallel} package. After finishing the estimation procedure,
 \code{siena07} will automatically stop the cluster. Alternatively, instead of
 having the function to create a cluster, the user may provide its own by
 specifying the option \code{cl}, similar to what the boot function does in
 the \CRANpkg{boot} package. By using the option \code{cl} the user may be
 able to create more complex clusters (see examples below).

 If \code{thetaValues} is not \code{NULL} and \code{nbrNodes >= 2}, parameters
 in Phase 3 will be constant for each set of \code{nbrNodes} consecutive
 simulations. This must be noted in the interpretation, and will be visible
 in \code{thetaUsed} (see below).

 The keyword \code{thetaBound} is used because in practice, parameter values
 of Stochastic Actor-oriented Models will be relatively small,
 and for usual models values larger than 30 never occur, which means
 that when they occur this is regarded as a signal of divergence
 of the algorithm.\cr
 Note that covariates should have large enough standard deviations;
 see the manual.
}
\value{
  Returns an object of class \code{\link{sienaFit}}, some parts of which are:
  \item{OK}{Boolean indicating successful termination}
  \item{termination }{Character string, values: "OK", "Error", or
    "UserInterrupt". "UserInterrupt" indicates that the user asked for early
  termination before phase 3. }
  \item{f}{Various characteristics of the data and model definition.}
  \item{requestedEffects}{The included effects in the effects object.}
  \item{effects}{The included effects in the effects object to which are added
    the main effects of the requested interaction effects, if any.}
  \item{theta}{Estimated value of theta, if \code{x$simOnly=FALSE}.}
  \item{thetas}{Matrix, returned if \code{returnThetas} and \code{x$nsub >= 1}.
    First column is subphase; further columns are values of theta as generated
    during this subphase of Phase 2.}
  \item{sfs}{Matrix, returned if \code{returnThetas} and \code{x$nsub >= 1}.
    First column is subphase; further columns are deviations from targets
    generated during this subphase of Phase 2.}
  \item{covtheta}{Estimated covariance matrix of theta; this is not available
    if \code{x$simOnly=TRUE}.}
  \item{se}{Vector of standard errors of estimated theta,
    if \code{x$simOnly=FALSE}.}
  \item{dfra}{Matrix of estimated derivatives.}
  \item{sf}{Matrix of simulated deviations from targets in phase 3.}
  \item{sf2}{Array of periodwise deviations from simulations in phase 3.
    Not included if \code{x$lessMem=TRUE}.}
  \item{W}{If \code{x$gmm=TRUE}: Estimated optimal matrix of weights for
    estimation by the Generalized Method of Moments.}
  \item{B}{If \code{x$gmm=TRUE}: Row-normalized matrix of weights for
    equating the linear combination of estimation statistics to 0,
    for estimation by the Generalized Method of Moments.}
  \item{tconv}{t-statistics for convergence.}
  \item{tmax}{maximum absolute t-statistic for convergence for
    non-fixed parameters.}
  \item{tconv.max}{overall maximum convergence ratio.}
  \item{ac3}{If \code{x$maxlike=TRUE}:
                    autocorrelations of statistics in Phase 3.}
  \item{targets}{Observed statistics; for ML, zero vector.}
  \item{targets2}{Observed statistics by wave,
    starting with second wave; for ML, zero matrix.}
  \item{ssc}{Score function contributions for each wave for each
    simulation in phase 3. Not included if finite difference method is used
    or if \code{x$lessMem=TRUE}.}
  \item{scores}{Score functions, added over waves, for each
    simulation in phase 3. Only included if \code{x$lessMem=FALSE}.}
  \item{regrCoef}{If \code{x$dolby} and not \code{x$maxlike}:
    regression coefficients of estimation statistics on score functions.}
  \item{regrCor}{If \code{x$dolby} and not \code{x$maxlike}:
    correlations between estimation statistics and score functions.}
  \item{estMeans}{Estimated means of estimation statistics.}
  \item{estMeans.sem}{If \code{x$simOnly}: Standard errors of the
    estimated means of estimation statistics.}
 \item{sims}{If \code{returnDeps=TRUE}:
    list of simulated dependent variables (networks, behaviour).
    Networks are given as a list of edgelists, one for each period.\cr
    The structure of sims is a nested list:
    \code{sims[[run]][[group]][[dependent variable]][[period]]}.
    If \code{x$maxlike=TRUE} and there is only one group and one period,
    the structure is \code{[[run]][[dependent variable]]}.}
 \item{chain}{If \code{returnChains = TRUE}:
    list, or data frame, of simulated chains of ministeps.
    The chain has the structure
    \code{chain[[run]][[depvar]][[period]][[ministep]]}.
    }
   \item{Phase3nits}{Number of iterations actually performed in phase 3.}
   \item{thetaUsed}{If \code{thetaValues} is not \code{NULL}, the matrix
   of parameter values actually used in the simulations of Phase 3.}
  Writes text output to the file named "projname.txt", where projname
  is defined in the \code{\link{sienaAlgorithm}} object \code{x}.
}
\references{
    Amati, V., Schoenenberger, F., and Snijders, T.A.B. (2015),
    Estimation of stochastic actor-oriented models for the evolution of
    networks by generalized method of moments.
    \emph{Journal de la Societe Francaise de Statistique} \bold{156}, 140--165.

    Amati, V., Schoenenberger, F., and Snijders, T.A.B. (2019),
    Contemporaneous statistics for estimation in stochastic actor-oriented
    co-evolution models. \emph{Psychometrika} \bold{84}, 1068--1096.

    Greenan, C. (2015),
    \emph{Evolving Social Network Analysis: developments in statistical
    methodology for dynamic stochastic actor-oriented models}.
    DPhil dissertation, University of Oxford.

    Niezink, N.M.D., and Snijders, T.A.B. (2017),
    Co-evolution of Social Networks and Continuous Actor Attributes.
    \emph{The Annals of Applied Statistics} \bold{11}, 1948--1973.

    Schweinberger, M., and Snijders, T.A.B. (2007),
    Markov models for digraph panel data:
    Monte Carlo based derivative estimation.
    \emph{Computational Statistics and Data Analysis} \bold{51}, 4465--4483.

     Snijders, T.A.B. (2001),
    The statistical evaluation of social network dynamics.
    \emph{Sociological Methodology} \bold{31}, 361--395.

    Snijders, T.A.B. (2017), Stochastic Actor-Oriented Models for Network
   Dynamics. \emph{Annual Review of Statistics and Its Application}
    \bold{4}, 343--363.

    Snijders, T.A.B., Koskinen, J., and Schweinberger, M. (2010).
    Maximum likelihood estimation for social network dynamics.
    \emph{Annals of Applied Statistics} \bold{4}, 567--588.

    Snijders, T.A.B., Steglich, C.E.G., and Schweinberger, Michael (2007),
    Modeling the co-evolution of networks and behavior.
    Pp. 41--71 in \emph{Longitudinal models in the behavioral
    and related sciences},
    edited by van Montfort, K., Oud, H., and Satorra, A.; Lawrence Erlbaum.

    Steglich, C.E.G., Snijders, T.A.B., and Pearson, M.A. (2010),
    Dynamic networks and behavior:     Separating selection from influence.
    \emph{Sociological Methodology} \bold{40}, 329--393.
    Information about the implementation of the algorithm is in
     \url{https://www.stats.ox.ac.uk/~snijders/siena/Siena_algorithms.pdf}.
    Further see \url{https://www.stats.ox.ac.uk/~snijders/siena/} .
}
\author{Ruth Ripley, Tom Snijders, Viviana Amati, Felix Schoenenberger,
        Nynke Niezink}

\seealso{\code{\link{siena}}, \code{\link{sienaAlgorithmCreate}},
  \code{\link{sienaEffects}},
  \code{\link{Wald.RSiena}}, \code{\link{Multipar.RSiena}},
  \code{\link{score.Test}}.

  There are print, summary and xtable methods for \code{\link{sienaFit}}
  objects: \code{\link{xtable}}, \code{\link{print.sienaFit}}.}
\examples{
myalgorithm <- sienaAlgorithmCreate(nsub=2, n3=100, seed=1293)
# nsub=2, n3=100 is used here for having a brief computation, not for practice.
mynet1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE)

# or for non-conditional estimation --------------------------------------------
\dontrun{
model <- sienaAlgorithmCreate(nsub=2, n3=100, cond=FALSE, seed=1283)
ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE)
        }

# or if a previous "on track" result ans was obtained --------------------------
\dontrun{
ans1 <- siena07(myalgorithm, data=mydata, effects=myeff, prevAns=ans)
         }

# Running in multiple processors -----------------------------------------------
\dontrun{
# Not tested because dependent on presence of processors
# Find out how many processors there are
library(parallel)
(n.clus <- detectCores() - 1)
n.clus <- min(n.clus, 4)  # keep time for other processes
ans2 <- siena07(myalgorithm, data=mydata, effects=myeff,
                useCluster=TRUE, nbrNodes=n.clus, initC=TRUE)

# Suppose 8 processors are going to be used.
# Loading the parallel package and creating a cluster
# with 8 processors (this should be equivalent)

library(parallel)
cl <- makeCluster(n.clus)

ans3 <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE, cl = cl)

# Notice that now -siena07- perhaps won't stop the cluster for you.
# stopCluster(cl)

# You can create even more complex clusters using several computers. In this
# example we are creating a cluster with 3*8 = 24 processors on three
# different machines.
#cl <- makePSOCKcluster(
#    rep(c('localhost', 'machine2.website.com' , 'machine3.website.com'), 8),
#    user='myusername', rshcmd='ssh -p PORTNUMBER')

#ans4 <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE, cl = cl)
#stopCluster(cl)
}

# for a continuous behavior variable -------------------------------------------
# simulate behavior data according to dZ(t) = [-0.1 Z + 1] dt + 1 dW(t)
set.seed(123)
y1 <- rnorm(50, 0,3)
y2 <- exp(-0.1) * y1 + (1-exp(-0.1)) * 1/ -0.1 + rnorm(50, 0, (exp(-0.2)- 1) / -0.2 * 1^2)
friend <- sienaDependent(array(c(s501, s502), dim = c(50,50,2)))
behavior <- sienaDependent(matrix(c(y1,y2), 50,2), type = "continuous")
(mydata <- sienaDataCreate(friend, behavior))
(myeff <- getEffects(mydata, onePeriodSde = TRUE))
algorithmMoM <- sienaAlgorithmCreate(nsub=1, n3=20, seed=321)
(ans <- siena07(myalgorithm, data = mydata, effects = myeff, batch=TRUE))

# Accessing simulated networks for ML ------------------------------------------
# The following is an example for accessing the simulated networks for ML,
# which makes sense only if there are some missing tie variables;
# observed tie variables are identically simulated
# at the moment of observation,
# missing tie variable are imputed in a model-based way.
mat1 <- matrix(c(0,0,1,1,
                 1,0,0,0,
                 0,0,0,1,
                 0,1,0,0),4,4, byrow=TRUE)
mat2 <- matrix(c(0,1,1,1,
                 1,0,0,0,
                 0,0,0,1,
                 0,0,1,0),4,4, byrow=TRUE)
mat3 <- matrix(c(0,1,0,1,
                 1,0,0,0,
                 0,0,0,0,
                 NA,1,1,0),4,4, byrow=TRUE)
mats <- array(c(mat1,mat2,mat3), dim=c(4,4,3))
net <- sienaDependent(mats, allowOnly=FALSE)
sdat <- sienaDataCreate(net)
alg <- sienaAlgorithmCreate(maxlike=TRUE, nsub=3, n3=100, seed=12534)
effs <- getEffects(sdat)
(ans <- siena07(alg, data=sdat, effects=effs, returnDeps=TRUE, batch=TRUE))
# See manual Section 9.1 for information about the following functions
edges.to.adj <- function(x,n){
# create empty adjacency matrix
    adj <- matrix(0, n, n)
# put edge values in desired places
    adj[x[, 1:2]] <- x[, 3]
    adj
}
the.edge <- function(x,n,h,k){
    edges.to.adj(x,n)[h,k]
}
# Now show the results
n <- 4
ego <- rep.int(1:n,n)
alter <- rep(1:n, each=n)
# Get the average simulated adjacency matrices for wave 3 (period 2):
ones <- sapply(1:n^2, function(i)
    {mean(sapply(ans$sims,
           function(x){the.edge(x[[1]][[2]][[1]],n,ego[i],alter[i])}))})
# Note that for maximum likelihood estimation,
# if there is one group and one period,
# the nesting levels for group and period are dropped from ans$sims.
cbind(ego,alter,ones)
matrix(ones,n,n)
}
\keyword{models}
