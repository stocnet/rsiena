\name{siena}
\alias{siena}
\alias{siena.sienadata}
\alias{siena07}
\concept{estimation}
\title{Function to estimate parameters in a Siena model}
\description{
 Estimates parameters in a Siena model using Method of Moments,
 based on direct simulation, conditional or otherwise;
 or using Generalized Method of Moments;
 or using Maximum Likelihood by MCMC simulation.
 Estimation is done using a Robbins-Monro algorithm.\cr
  \code{siena()} and \code{siena07()}
  perform the same estimation. The first of these is preferred,
  the other is kept for backward compatibility.
}
\usage{
\method{siena}{sienadata}(x, effects=NULL,
        control_model=NULL, control_algo=NULL, control_out=NULL,
        thetaBound = 50,
        batch = FALSE, verbose = FALSE, silent=TRUE,
        initC=TRUE,
        nbrNodes = 1, clusterString=rep("localhost", nbrNodes),
        clusterType=c("PSOCK", "FORK"), cl=NULL, ...)

siena07(alg, batch=FALSE, verbose=FALSE, silent=FALSE,
        useCluster=FALSE, nbrNodes=2,
        thetaValues = NULL,
        returnThetas = FALSE,
        thetaBound = 50,
        targets = NULL,
        initC=TRUE, clusterString=rep("localhost", nbrNodes), tt=NULL,
        parallelTesting=FALSE, clusterIter=!x$maxlike,
        clusterType=c("PSOCK", "FORK"), cl=NULL, ...)
      }
\arguments{
  \item{x}{Siena data set of class \code{\link{sienadata}}. }
  \item{effects}{\code{sienaEffects} object specifying which effects are
    included the model.}
  \item{control_model}{Control object for the model,
                       of class \code{\link{sienaModel}}. }
  \item{control_algo}{Control object for the algorithm,
                       of class \code{\link{sienaAlgorithmSettings}}. }
  \item{control_out}{Control object for the output,
                       of class \code{\link{sienaOutputOptions}}. }
  \item{thetaBound}{Numeric: if at some moment during estimation the maximum
   absolute value of the parameters exceeds \code{thetaBound}, estimation
   is interrupted. In interactive use, the user is requested to give
   a higher number; if non-numeric input is given, estimation stops.
   In non-interactive use, estimation stops anyway.
   The check is turned off by using \code{thetaBound=Inf}.}
  \item{batch}{ Desired interface: \code{FALSE} gives a gui
     (graphical user interface implemented as a tcl/tk screen),
     \code{TRUE} gives a small (if \code{verbose=FALSE}) amount of printout
     to the console.}
  \item{verbose}{Produces various output to the console if \code{TRUE}.}
  \item{silent}{Produces no output to the console if \code{TRUE},
     even if batch mode.}
  \item{initC}{Boolean: set to \code{TRUE} if the simulation will use C
   routines (currently always needed). Only for use if using
   multiple processors, to ensure all copies are initialised
   correctly. Ignored otherwise, so is set to \code{TRUE} by default.}
  \item{nbrNodes}{Number of processes to use; when using \code{siena07},
        this is ignored if useCluster is \code{FALSE}.}
  \item{clusterString}{Definitions of clusters. Default set up to use
   the local machine only.}
  \item{clusterType}{Either "PSOCK" or "FORK". On Windows, must be
  "PSOCK". On a single non-Windows machine may be "FORK", and
  subprocesses will be formed by forking. If "PSOCK", subprocesses are
  formed using \R scripts.}
  \item{cl}{An object of class c("SOCKcluster", "cluster") (see Details).}
  \item{alg}{A control object, of class \code{\link{sienaAlgorithm}. }}
  \item{useCluster}{Boolean: whether to use a cluster of processes (useful if
  multiple processors are available).}
  \item{thetaValues}{If not \code{NULL}, this should be a matrix with parameter
    values to be used in Phase 3. The number of columns must be equal to the
    number of estimated parameters in the effects object (if conditional
    estimation is used, without the rate parameters for the conditioning
    dependent variable). Can only be used if \code{x$simOnly=TRUE}.}
  \item{returnThetas}{Boolean: whether to return theta values
   and generated estimation statistics of Phase 2 runs.}
  \item{targets}{Numeric vector of length equal to the number of estimated
    parameters, meant to supersede the targets calculated from the data set;
    see "Details". Not for regular use.}
  \item{tt}{A \code{tcltk} toplevel window. Used if called from the
   model options screen, if \code{tcltk} is available.}
  \item{parallelTesting}{Boolean. If \code{TRUE}, sets up random numbers to
   parallel those in Siena 3.}
  \item{clusterIter}{Boolean. If \code{TRUE}, multiple processes execute
   complete iterations at each call.
   If \code{FALSE}, multiple processes execute a single wave at each call.}
  \item{\dots}{Arguments for the simulation function,
    see \code{\link{simstats0c}}:\cr
    \code{prevAns} if a previous reasonable
    provisional estimate was obtained for a similar model;\cr
    the following are ignored by \code{siena}
    (because included in \code{control_out}) but an be mentioned when using \code{siena07}:\cr
    in any case, \code{data} and \code{effects}, as in the examples below;\cr
    possibly also \code{returnDeps} if the simulated dependent variables
    (networks, behaviour) should be returned;\cr
    possibly also \code{returnChains} or \code{returnChangeContributions}
    if the simulated sequences (chains) of ministeps or the simulated
    change statistics (contributions) should be returned;
    this may produce a very big file.}
}
\details{
  This is the main function and workhorse of RSiena.

  For use of \code{siena07}, it is necessary to specify parameters \code{data}
  (RSiena data set) and \code{effects} (effects object), which are
  required parameters in function \code{\link{simstats0c}}.
  (These parameters are inserted through '\dots'.) See the examples.

 \code{siena} runs a Robbins-Monro algorithm for parameter estimation
  using the three-phase implementation described in Snijders (2001, 2017),
 with (if \code{findiff=FALSE} in the algorithm object)
 derivative estimation as in Schweinberger and Snijders (2007).
 The default is estimation according to the Method of Moments
 as in Snijders, Steglich and Schweinberger (2007).\cr
 If  \code{gmm=TRUE} in the algorithm object and \code{myeff} contains one or more
 \code{gmm} statistics as included by function
 \code{\link{includeGMoMStatistics}}, the algorithm employs the
 Generalized Method of Moments as defined in Amati, Schoenenberger, and
 Snijders (2015, 2019).\cr
 For continuous behavior variables defined with \code{type="continuous"} in \cr
 \code{\link{as_dependent_rsiena}}, estimation is done as
 described in Niezink and Snijders (2017).\cr
 If  \code{maxlike=TRUE} in the algorithm object, estimation is done by Maximum Likelihood
 implemented as in Snijders, Koskinen and Schweinberger (2010).\cr
 Phase 1 does a few iterations to
 estimate the derivative matrix of the targets with respect to the
 parameter vector. Phase 2 does the estimation. Phase 3 runs a
 simulation to estimate standard errors and check convergence of the model.
 The simulation function is called once for each iteration in these phases
 and also once to initialise the model fitting and once to complete it.
 Unless in batch mode, a tcl/tk screen is displayed to allow interruption
 and to show progress.

 If \code{targets} is specified in the algorithm object for \code{siena}
 or in the call of \code{siena07} (which should be done only in special cases),
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

 For good estimation of standard errors, it is necessary that \code{n3}
 in the algorithm object is large enough. More about this is in the manual.
 The default value \code{n3} is adequate for most explorative use,
 but for presentation in publications larger values are necessary, depending
 on the data set and model; e.g., \code{n3=3000} or larger.

 Parameters can be tested against zero by dividing the estimate by its
 standard error and using an approximate standard normal null distribution.
 Further, function \code{\link{test_parameter}} is available for
 multi-parameter testing.\cr
 Parameters specified in \code{\link{set_effect}} with
 \code{fix=TRUE, test=TRUE} will not be estimated;
 score tests of their hypothesized values can be obtained
 using \code{\link{test_parameter}} with \code{method="score"}.

 If \code{simOnly} is \code{TRUE} in the algorithm object,
 which is meant to go together with  \code{nsub=0},
 the calculation of the standard errors and covariance
 matrix at the end of Pase 3 is skipped. No estimation is performed.
 If \code{thetaValues} in the algorithm object is not \code{NULL},
 the parameter values in  the rows of this matrix will be used
 in the consecutive runs of Phase 3.
 If \code{n3} is larger than the number of rows times \code{nbrNodes}
 (see below), the last row of \code{thetaValues} will continue to be used.
 The parameter values actually used will be stored in the
 output matrix \code{thetaUsed}.

 Multiple processors are used for estimation by MoM to distribute each iteration
 in each subphase over the cluster of nodes. The number of iterations
 accordingly will be divided (approximately) by the number of nodes;
 for phase 2, unless \code{n2start} is specified in the algorithm object.
 This implies that if multiple processors are used, think of dividing
 \code{n2start} by \code{nbrNodes}.\cr
 For estimation by ML, multiple processing is done per period. Therefore,
 for one period (two waves) and one group, this will have no effect.

 In the case of using multiple processors, there are two options for telling
 \code{siena} or \code{siena07} to use them.
 By specifying the options \code{useCluster},
 \code{nbrNodes}, \code{clusterString} and \code{initC}, the estimation
 will create a \code{cluster object} that will be used by the
 \code{parallel} package. After finishing the estimation procedure,
 the estimation will automatically stop the cluster. Alternatively, instead of
 having the function to create a cluster, the user may provide its own by
 specifying the option \code{cl}, similar to what the boot function does in
 the \CRANpkg{boot} package. By using the option \code{cl} the user may be
 able to create more complex clusters (see examples below).

 If \code{thetaValues} is not \code{NULL} and \code{nbrNodes >= 2}, parameters
 in Phase 3 will be constant for each set of \code{nbrNodes} consecutive
 simulations. This must be noted in the interpretation, and will be visible
 in \code{thetaUsed} (see below).

 The keyword \code{thetaBound} is used because in practice, parameter values
 of Stochastic Actor-oriented Models will be limited in magnitude,
 and for usual models values larger than 50 never occur, which means
 that when they occur this may be regarded as a signal of divergence
 of the algorithm.\cr
 Note that covariates should have large enough standard deviations
 for parameters to have moderate values; see the manual.
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
    list, or data frame, of simulated chains of ministeps. \cr
    The chain has the structure
    \code{chain[[run]][[depvar]][[period]][[ministep]]}.
    }
  \item{changeContributions}{If \code{returnChangeContributions = TRUE}:
    list, or data frame, of simulated change contributions. \cr
    The contributions have the structure
    \code{returnChangeContributions[[run]][[depvar]][[period]][[ministep]]}.
    }
   \item{Phase3nits}{Number of iterations actually performed in phase 3.}
   \item{thetaUsed}{If \code{thetaValues} is not \code{NULL}, the matrix
   of parameter values actually used in the simulations of Phase 3.}
  Writes text output to the file indicated for \code{siena}
  in \code{control_output} and for \code{siena07} in \code{alg}.
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
     \url{https://www.stats.ox.ac.uk/~snijders/siena/Siena_algorithms.pdf}.\cr
    Further see \url{https://www.stats.ox.ac.uk/~snijders/siena/} .
}
\author{Ruth Ripley, Tom Snijders, Viviana Amati, Felix Schoenenberger,
        Nynke Niezink}

\seealso{\code{\link{sienadata}}, \code{\link{set_algorithm_saom}},
  \code{\link{sienaEffects}},
  \code{\link{test_parameter}}.
  There are print, summary and xtable methods for \code{\link{sienaFit}}
  objects: \code{\link{xtable}}, \code{\link{print.sienaFit}}.}
\examples{
myalgo <- set_algorithm_saom(nsub=2, n3=100, seed=1293)
# nsub=2, n3=100 is used here for having a brief computation, not for practice.
mynet1 <- as_dependent_rsiena(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- make_data_rsiena(mynet1)
myeff <- make_specification(mydata)
ans <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE)
# or, when using the old siena07,
myalgorithm <- sienaAlgorithmCreate(nsub=2, n3=100, seed=1293)
ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE)
ans
coef(ans)
# or for non-conditional estimation ---------------------------------
myalgo <- set_algorithm_saom(nsub=2, n3=100, cond=FALSE, seed=1273)
ans <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE)
coef(ans)
coef(ans, dropRates=TRUE)
# or if a previous "on track" result ans was obtained ---------------
\dontrun{
ans1 <- siena(mydata, effects=myeff, control_algo=myalgo, prevAns=ans)
         }

# Running in multiple processors ------------------------------------
\dontrun{
# Not tested because dependent on presence of processors
# Find out how many processors there are
library(parallel)
(n.clus <- detectCores() - 1)
n.clus <- min(n.clus, 4)  # keep time for other processes
ans2 <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE,
                nbrNodes=n.clus)

# Suppose 8 processors are going to be used.
# Loading the parallel package and creating a cluster
# with 8 processors (this should be equivalent)

library(parallel)
cl <- makeCluster(n.clus)

ans3 <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE, cl = cl)
# I did not succeed in constructing an alternative with siena for this.
# Notice that now -siena- perhaps won't stop the cluster for you.
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

# for a continuous behavior variable --------------------------------
# simulate behavior data according to dZ(t) = [-0.1 Z + 1] dt + 1 dW(t)
set.seed(123)
y1 <- rnorm(50, 0,3)
y2 <- exp(-0.1) * y1 + (1-exp(-0.1)) * 1/ -0.1 + 
                       rnorm(50, 0, (exp(-0.2)- 1) / -0.2 * 1^2)
friend <- as_dependent_rsiena(array(c(s501, s502), dim = c(50,50,2)))
behavior <- as_dependent_rsiena(matrix(c(y1,y2), 50,2), type = "continuous")
(mydata <- make_data_rsiena(friend, behavior))
(myeff <- make_specification(mydata, onePeriodSde = TRUE))
algorithmMoM <- set_algorithm_saom(nsub=1, n3=20, seed=321)
(ans <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE))

# Accessing simulated networks for ML -------------------------------
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
net <- as_dependent_rsiena(mats, allowOnly=FALSE)
sdat <- make_data_rsiena(net)
alg <- set_algorithm_saom(maxlike=TRUE, nsub=3, n3=100, seed=12534)
effs <- make_specification(sdat)
(ans <- siena(sdat, effects=effs, control_algo=alg,
                              batch=TRUE, returnDeps=TRUE))
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
