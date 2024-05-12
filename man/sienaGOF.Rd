\name{sienaGOF}
\alias{sienaGOF}
\alias{plot.sienaGOF}
\alias{descriptives.sienaGOF}
\title{Functions to assess goodness of fit for SAOMs}
\description{
 The function \code{sienaGOF} assesses goodness of fit
 for a model specification as represented by an estimated
 \code{sienaFit} object created by \code{\link{siena07}}.
 This is done by simulations of auxiliary
 statistics, that differ from the statistics used for estimating
 the parameters. The auxiliary statistics must be given explicitly.

 The fit is good if the average values
 of the auxiliary statistics over many simulation runs are close
 to the values observed in the data.
 A Monte Carlo test based on the Mahalanobis
 distance is used to calculate frequentist \eqn{p}-values.

 Plotting functions can be used to diagnose bad fit. There are basic
 functions for calculating auxiliary statistics available out of the
 box, and the user is invited to create additional ones.
 }

\usage{
sienaGOF(sienaFitObject, auxiliaryFunction,
        period=NULL, verbose=FALSE, join=TRUE, twoTailed=FALSE,
        cluster=NULL, robust=FALSE, groupName="Data1",
        varName, tested=NULL, iterations=NULL, giveNAWarning=TRUE, \dots)
\method{plot}{sienaGOF}(x, center=FALSE, scale=FALSE, violin=TRUE, key=NULL,
        perc=.05, period=1, position=4, fontsize=12, \dots)
descriptives.sienaGOF(x, center=FALSE, scale=FALSE, perc=.05, key=NULL,
        period=1, showAll=FALSE)
}
\arguments{
  \item{sienaFitObject}{An object of class \code{\link{sienaFit}}, produced by a
    call to \code{\link{siena07}} with \code{returnDeps = TRUE} and
    \code{maxlike=FALSE} (the latter is the default, the former is not);
    or a list of such objects; if a list, then the first period
    of each \code{\link{sienaFit}} object will be used.\cr
    If this is a list of \code{\link{sienaFit}} objects,
    where \code{sienaFitObject} is mentioned below, it refers to the first
    element of this list.}
  \item{auxiliaryFunction}{Function to be used to calculate
    the auxiliary statistics;
    this can be a user-defined function, e.g. depending on the \code{sna}
    or \code{igraph} packages.

    See \code{Examples} and \code{\link{sienaGOF-auxiliary}}
    for more information on the signature of this function.
    The basic signature is\cr
    \code{function(index, data, sims, period, groupName, varName, \dots)},
    where \code{index} is the index of the simulated network,
    or \code{NULL} if the observed variable is needed;
    \code{data} is the observed data object from which
    the relevant variables are extracted; \code{sims} is the list of simulations
    returned from \code{siena07}; \code{period} is the index of the period;
    and \dots are further arguments (like \code{levls} in the examples below
    and in \code{\link{sienaGOF-auxiliary}}).}
  \item{period}{Vector of period(s) to be used (may run from 1 to
    number of waves - 1). Has an effect only if \code{join=FALSE}.\cr
    May be only a single number if \code{sienaFitObject} is a list of
    \code{\link{sienaFit}} objects.}
  \item{verbose}{Whether to print intermediate results. This may give some
       peace of mind to the user because calculations can take some time.}
  \item{join}{Boolean: should sienaGOF do tests on all of the periods
       individually (\code{FALSE}), or sum across periods (\code{TRUE})?}
  \item{twoTailed}{Whether to use two tails for calculating \eqn{p}-values on
        the Monte Carlo test. Recommended for advanced users only,
        as it is probably only applicable in rare cases.}
  \item{cluster}{Optionally, a \code{parallel} or \code{snow} cluster to execute
        the auxiliary function calculations on.}
  \item{robust}{ Whether to use robust estimation of the covariance matrix.}
  \item{groupName}{ Name of group; relevant for multi-group data sets.}
  \item{varName}{ Name of dependent variable.}
  \item{tested}{ A logical vector of length \code{sienaFitObject$pp}
          (number of parameters), indicating a subset of tested parameters;
          or \code{NULL}, indicating all tested parameters (see below);
          or \code{FALSE}, indicating nothing is to be tested.}
  \item{iterations}{ Number of iterations for the goodness of fit calculations.
        If \code{NULL}, the number of simulated data sets in
        \code{sienaFitObject}.  }
  \item{giveNAWarning}{If \code{TRUE}, a warning is given if any simulated
        values are missing.}
  \item{x}{ Result from a call to sienaGOF. }
  \item{center}{ Whether to center the statistics by median during plotting.}
  \item{scale}{ Whether to scale the statistics by range during plotting.
                scale=\code{TRUE} makes little sense without also
                center=\code{TRUE}.}
  \item{violin}{ Use violin plots (vs. box plots only)? }
  \item{key}{ Keys in the plot for the levels of the auxiliary statistic
        (as given by parameter \code{levls} in the examples). }
  \item{perc}{ 1 minus confidence level for the confidence bands (two sided). }
  \item{position}{ Position where the observed value is plotted: 1=under,
        2=to the left, 3=above, 4=to the right of the red dot. Can be a single
        number from 1 to 4, or a vector with positions for each statistic
        (possibly recycled).}
  \item{fontsize}{ Font size for the observed values plotted.}
  \item{\dots}{Other arguments; for \code{sienaGOF()}, e.g., \code{levls} as a
       parameter for the auxiliary statistic in
       \code{\link{sienaGOF-auxiliary}};\cr
       for \code{plot.sienaGOF()}, e.g., the usual plotting parameters
       \code{main}, \code{xlab}, \code{ylab}, \code{cex}, \code{cex.main},
       \code{cex.lab}, and \code{cex.axis}.}
  \item{showAll}{If \code{FALSE}, drops statistics with variance 0,
       like in the plot.}
}
\details{
 This function is used to assess the goodness of fit of an estimated stochastic
 actor-oriented model for an arbitrarily defined multidimensional auxiliary
 statistic. It operates basically by comparing the observed values, at the ends
 of the periods, with the simulated values for the ends of the periods.
 The differences are assessed by combining the components of the auxiliary
 statistic using the Mahalanobis distance.

 For \code{sienaFitObject}s that were made for a multi-group data set,
 if you are not sure about the \code{groupName}s to use, these can be retrieved
 by the command \code{"names(dataObject)"} (where \code{dataObject} is
 the data used to produce the \code{sienaFitObject}).
 Mostly they are \code{"Data1", "Data2"}, etc.

 To save computation time, \code{iterations} can be set to a lower number
 than what is available in \code{sienaFitObject}; this will yield
 a less precise result.

 The function does not work properly for data sets that include a
 \code{\link{sienaCompositionChange}} object. If you wish to test the
 fit for such a data set, you need (for the purpose of fit assessment only)
 to replace the data set by a data set where absent actors are represented
 by structural zeros, and estimate the same model for this data set with the
 corresponding effects object, and use \code{sienaGOF} for this
 \code{\link{sienaFit}} object.

 To achieve comparability between simulated and observed dependent variables,
 variables that are missing in the data at the start or end of a period
 are replaced by 0 (for tie variables) or NA (for behavior variables).\cr
 If there are any differences between structural values at the beginning and
 at the end of a period, these are dealt with as follows.
 For tie variables that have a structural value at the start of the period,
 this value is used to replace the observed value at the end of the period
 (for the goodness of fit assessment only).
 For tie variables that have a structural value at the end of the period but a
 free value value at the start of the period, the reference value for the
 simulated values is lacking; therefore, the simulated values at the end of
 the period then are replaced by the structural value at the end of the period
 (again, for the goodness of fit assessment only).

 The auxiliary statistics documented in \code{\link{sienaGOF-auxiliary}}
 are calculated for the simulated dependent
 variables in Phase 3 of the estimation algorithm,
 returned in \code{sienaFitObject} because of having used
 \code{returnDeps = TRUE} in the call to \code{\link{siena07}}.
 These statistics should be chosen to represent features of the network
 that are not explicitly fit by the estimation procedure but can be
 considered important properties that the model at hand should represent well.
 Some examples are:
 \itemize{
 \item Outdegree distribution
 \item Indegree distribution
 \item Distribution of the dependent behavior variable (if any).
 \item Distribution of geodesic distances
 \item Triad census
 \item Edgewise homophily counts
 \item Edgewise shared partner counts
 \item Statistics depending on the combination of network
       and behavioral variables.
}
The function is written so that the user can easily define other functions to
capture some other relevant aspects of the network, behaviors, etc.
This is further illustrated in the help page
\code{\link{sienaGOF-auxiliary}}.

We recommend the following heuristic approach to model checking:
\enumerate{
 \item Check convergence of the estimation.
 \item Assess goodness of fit (primarily using \code{join=TRUE}) on
  auxiliary statistics, and if necessary refine the model.
 \item Assess time heterogeneity by \code{\link{sienaTimeTest}}
   and if there is evidence for time heterogeneity
   either modify the base effects or include time dummy terms.
}
 No general rules can be given about whether time heterogeneity
 (\code{\link{sienaTimeTest}}) or goodness of fit using \code{sienaGOF}
 have precedence. This is an explorative issue.

 The \code{summary} function will display some useful information
 to help with model selection if some effects are set in the effects object
 to be fixed and tested. In that case, for all parameters indicated
 in the vector \code{tested}, a rough estimator is computed for the Mahalanobis
 distance that would be obtained at each proposed specification.
 This is then given in the \code{summary}. This can help guide model selection.
 This estimator is called the modified Mahalanobis distance (MMD).
 See Lospinoso and Snijders (2019) or the manual for more information.

 The following functions are pre-fabricated for ease of use, and can
 be passed in as the \code{auxiliaryFunction} with no extra effort;
 see \code{\link{sienaGOF-auxiliary}} and the examples below.
\itemize{
 \item \code{\link{IndegreeDistribution}}
 \item \code{\link{OutdegreeDistribution}}
 \item \code{\link{BehaviorDistribution}}
 \item \code{\link{TriadCensus}}
 \item \code{\link{mixedTriadCensus}}
 \item \code{\link{dyadicCov}}
}

}
\value{
  \code{sienaGOF} returns a result of class \code{sienaGOF};
    this is a list of elements of class \code{sienaGofTest};
    if \code{join=TRUE}, the list has length 1; if \code{join=FALSE},
    each list element corresponds to a period analyzed;
    the list elements are themselves lists again, including the
    following elements:
    \item{- sienaFitName}{The name of sienaFitObject.}
    \item{- auxiliaryStatisticName}{The name of auxiliaryFunction.}
    \item{- Observations}{The observed values for the auxiliary statistics.}
    \item{- Simulations}{The simulated auxiliary statistics.}
    \item{- ObservedTestStat}{The observed Mahalanobis distance in the data.}
    \item{- SimulatedTestStat}{The Mahalanobis distance for the simulations.}
    \item{- TwoTailed}{Whether the \eqn{p}-value corresponds to a one- or
                     two-tailed Monte Carlo test.}
    \item{- p}{The \eqn{p}-value for the observed Mahalanobis distance in the
             permutation distribution of the
             simulated Mahalanobis distances.}
    \item{- Rank}{ Rank of the covariance matrix of the simulated auxiliary
              statistics.}
    In addition there are several attributes which give, for model
    specifications with fixed-and-tested effects, approximations
    to the expected Mahalanobis distance for model specifications where
    each of these effects would be added. This is reported
    in the \code{summary} method. \cr
    The \code{plot} method makes violin plots or box plots, with superimposed
    confidence bands, for the simulated distributions of all elements
    of the \code{auxiliaryFunction},
    with the observed values indicated by red dots;
    but statistics with variance 0 are dropped.

    \code{descriptives.sienaGOF} returns a matrix giving numerical
    information about what is plotted in the \code{plot} method:
    maximum, upper percentile, mean, median, lower percentile, minimum,
    and standard deviation
    of the simulated distributions of the auxiliary statistics, the
    observed values, and the proportions of simulated values greater
    and greater-or-equal than the observed values.
    If \code{center=TRUE} the median is subtracted
    from mean, median, and percentiles; if \code{scale=TRUE}
    these numbers and the standard deviation
    are divided by (maximum - minimum). \cr
    If \code{showAll=FALSE}, statistics with variance 0 will be dropped.
}

\references{
 Lospinoso, J.A. and Snijders, T.A.B. (2019,
  Goodness of fit for stochastic actor-oriented models.
  \emph{Methodological Innovations}, \bold{12}:2059799119884282.

  Also see \url{https://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Josh Lospinoso, modifications by Ruth Ripley and Tom Snijders}
\seealso{\code{\link{siena07}}, \code{\link{sienaGOF-auxiliary}},
         \code{\link{sienaTimeTest}} }
\examples{
mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mybeh <- sienaDependent(s50a[,1:2], type="behavior")
mydata <- sienaDataCreate(mynet, mybeh)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip)
myeff <- setEffect(myeff, cycle3, fix=TRUE, test=TRUE)
myeff <- setEffect(myeff, transTies, fix=TRUE, test=TRUE)
myalgorithm <- sienaAlgorithmCreate(nsub=1, n3=10, projname=NULL)
# Shorter phases 2 and 3, just for example.
ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE, returnDeps=TRUE)
gofi <- sienaGOF(ans, IndegreeDistribution, verbose=TRUE, join=TRUE,
  varName="mynet")
summary(gofi)
plot(gofi)

# Illustration just for showing a case with two dependent networks;
# running time backwards is not meaningful!
mynet1 <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaDependent(array(c(s503, s501), dim=c(50, 50, 2)))
mybeh <- sienaDependent(s50a[,1:2], type="behavior")
mydata <- sienaDataCreate(mynet1, mynet2, mybeh)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip)
myeff <- includeEffects(myeff, recip, name="mynet2")
# Shorter phases 2 and 3, just for example.
ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE, returnDeps=TRUE)
gofi <- sienaGOF(ans, IndegreeDistribution, verbose=TRUE, join=TRUE,
  varName="mynet1")
summary(gofi)
plot(gofi)

\dontrun{
(gofi.nc <- sienaGOF(ans, IndegreeDistribution, cumulative=FALSE,
    varName="mynet1"))
# cumulative is an example of "...".
plot(gofi.nc)
descriptives.sienaGOF(gofi.nc)

(gofi2 <- sienaGOF(ans, IndegreeDistribution, varName="mynet2"))
plot(gofi2)

(gofb <- sienaGOF(ans, BehaviorDistribution, varName = "mybeh"))
plot(gofb)

(gofo <- sienaGOF(ans, OutdegreeDistribution, varName="mynet1",
    levls=0:6, cumulative=FALSE))
# levls is another example of "...".
plot(gofo)
}

## A demonstration of using multiple processes
\dontrun{
library(parallel)
(n.clus <- detectCores() - 1)
n.clus <- min(n.clus, 4)  # keep time for other processes
myalgorithm.c <- sienaAlgorithmCreate(nsub=4, n3=1000, seed=1265)
(ans.c <- siena07(myalgorithm.c, data=mydata, effects=myeff, batch=TRUE,
    returnDeps=TRUE, useCluster=TRUE, nbrNodes=n.clus))
gofi.1 <- sienaGOF(ans.c, TriadCensus, verbose=TRUE, varName="mynet1")
cl <- makeCluster(n.clus)
gofi.cl <- sienaGOF(ans.c, TriadCensus, varName="mynet1", cluster=cl)
cl2 <- makeCluster(2)
gofi.cl2 <- sienaGOF(ans.c, TriadCensus, varName="mynet1", cluster=cl2)
# compare simulation times
attr(gofi.1,"simTime")
attr(gofi.cl,"simTime")
attr(gofi.cl2,"simTime")
}

}
\keyword{models}
