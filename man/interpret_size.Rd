\name{interpret_size}
\alias{interpret_size}
\alias{interpret_size.sienaFit}
\alias{interpret_size.sienaEffects}
\alias{plot.interpret_size}
\alias{print.interpret_size}
\alias{plot.sienaRI}
\alias{print.sienaRI}
\alias{interpret_size.methods}
\title{Functions to assess the relative importance of effects at observation
       moments}
\description{
 The function \code{interpret_size} computes various effect sizes for
 Stochastic Actor-oriented Models.
 As its value it returns the relative importance of effects
 of a SAOM according to the measure of relative importance described in
 Section 3.1 of Indlekofer and Brandes (2013).
 This is based on the influence of effects on potential tie change or
 behavior change decisions of individual actors at the given observation
 moments (waves).  \cr
 Other measures useful for effect sizes are entropy-based
 effect sizes as in Snijders (2004) and the within-ego
 standard deviations of change statistics.
 If \code{getChangeStats=TRUE}, the arrays of change statistics are
 stored in the \code{sienaRI} object that is created.

 The \code{sienaFit} method delegates the core RI computation to
 \code{\link{relativeImportance}}, which can also be called directly
 when only the RI values (plus optional uncertainty) are needed.

 It takes the data as well as the complete model specification into account.
 Therefore, necessary arguments are the analysed data given as a
 \code{\link{sienadata}} data object as well as the complete model
 specification represented either by an estimated \code{sienaFit}
 object or by the combination of a
 suitable parameter vector \code{theta} and the corresponding
 \code{\link{sienaEffects}} object.
}
\usage{
\method{interpret_size}{sienaFit}(object, data, getChangeStats=FALSE, depvar=NULL,
                                  effects=NULL, dynamic=FALSE,
                                  distFun=c("L1D","KLD"), algorithm=NULL,
                                  n3=500, useChangeContributions=FALSE,
                                  uncertainty=FALSE, level="ego",
                                  sum_fun=mean, na.rm=TRUE, nsim=1000,
                                  uncertaintySd=TRUE, uncertaintyCi=TRUE,
                                  uncertaintyMean=FALSE, uncertaintyMedian=FALSE,
                                  uncertaintyProbs=c(0.025,0.5,0.975),
                                  useCluster=FALSE, nbrNodes=1,
                                  clusterType=c("PSOCK","FORK"),
                                  batchDir="temp", batchSize=NULL,
                                  keepBatch=FALSE, verbose=TRUE,
                                  memoryScale=NULL, batchUnitBudget=2.5e8,
                                  dynamicMinistepFactor=10, ...)
\method{interpret_size}{sienaEffects}(object, data, theta, getChangeStats=FALSE, ...)
\method{print}{sienaRI}(x, printSigma=FALSE,\dots)
\method{plot}{sienaRI}(x, actors = NULL, col = NULL, addPieChart = FALSE,
              radius = 1, width = NULL, height = NULL, legend = TRUE,
              legendColumns = NULL, legendHeight = NULL,
              cex.legend = NULL, cex.names = NULL,
              intervalsPerPeriod = 20, \dots)
}
\arguments{
  \item{object}{for \code{interpret_size.sienaFit}: a \code{sienaFit} object,
           resulting from a call of \code{\link{siena}};\cr
           for \code{interpret_size.sienaEffects}: a \code{sienaEffects} object
           that was used for calling \code{\link{siena}}.}
  \item{x}{for \code{print.sienaRI} and \code{plot.sienaRI}:
           a \code{sienaRI} object resulting from a call to \code{interpret_size}.}
  \item{data}{\code{siena} data object representing the analyzed data and
              resulting from a call to \code{\link{make_data_rsiena}}.}
  \item{theta}{Vector of parameter values of evaluation effects
               included in the model.
               Length of \code{theta} has to be equal to the number of
               included evaluation effects in the effects object \code{object}.}
  \item{getChangeStats}{Boolean: If \code{TRUE}, an array of change statistics
               is added to the \code{interpret_size} object.}
  \item{depvar}{Character: name of dependent variable to calculate RI for.
               If \code{NULL}, all dependent variables are used.}
  \item{effects}{Optional \code{sienaEffects} object. If \code{NULL},
               uses \code{object$requestedEffects}.}
  \item{dynamic}{Boolean: If \code{TRUE}, compute relative importance
               along simulated dynamic chains rather than at observation
               moments. This replaces \code{\link{interpret_size_dynamics}}
               which is now deprecated.}
  \item{distFun}{Character: distance function used (\code{"L1D"} or \code{"KLD"}).}
  \item{algorithm}{\code{sienaAlgorithm} object used for the dynamic
               simulation. Required when \code{dynamic=TRUE}.}
  \item{n3}{Integer: number of chains to simulate for dynamic RI
               (default 500).}
  \item{useChangeContributions}{Boolean: If \code{TRUE}, reuse change
               contributions stored in \code{object}.}
  \item{uncertainty}{Boolean: If \code{TRUE}, uncertainty quantification
               via theta-draws from \eqn{N(\hat\theta, \Sigma_\theta)} is computed.}
  \item{level}{Character: aggregation level (\code{"ego"}).}
  \item{sum_fun}{Function used to aggregate RI values across actors (default \code{mean}).}
  \item{na.rm}{Boolean: whether to remove \code{NA}s in aggregation (default \code{TRUE}).}
  \item{nsim}{Integer: number of theta draws for uncertainty quantification (default 1000).}
  \item{uncertaintySd}{Boolean: include posterior SD in uncertainty output.}
  \item{uncertaintyCi}{Boolean: include credible interval in uncertainty output.}
  \item{uncertaintyMean}{Boolean: include mean of draws (default \code{FALSE}).}
  \item{uncertaintyMedian}{Boolean: include median of draws (default \code{FALSE}).}
  \item{uncertaintyProbs}{Numeric vector: quantile probabilities for credible intervals.}
  \item{useCluster}{Boolean: if \code{TRUE}, use a parallel cluster.}
  \item{nbrNodes}{Integer: number of cluster nodes (default 1).}
  \item{clusterType}{Character: cluster type (\code{"PSOCK"} or \code{"FORK"}).}
  \item{batchDir}{Character: directory for batch simulation files (default \code{"temp"}).}
  \item{batchSize}{Integer or \code{NULL}: number of simulations per batch.}
  \item{keepBatch}{Boolean: keep batch files after completion (default \code{FALSE}).}
  \item{verbose}{Boolean: print progress messages (default \code{TRUE}).}
  \item{memoryScale}{Numeric or \code{NULL}: memory scaling factor.}
  \item{batchUnitBudget}{Numeric: per-unit memory budget (default \code{2.5e8}).}
  \item{dynamicMinistepFactor}{Numeric: scaling factor for the number
               of ministeps to simulate in dynamic mode (default 10).}
  \item{printSigma}{Boolean: If \code{TRUE},
       average within-actor standard deviations of change statistics
       ('\code{sigma}') are included in the print.}
  \item{actors}{vector of integers: set of actors to be included in the plot;
       if \code{NULL}, all actors.}
  \item{col}{Colors used in the plot. If \code{col=NULL} a default color
       scheme is used.}
  \item{addPieChart}{Boolean: If \code{TRUE}, pie charts of aggregated relative
       importances for the complete set of actors will be added to the plot.}
  \item{radius}{Radius of pie charts. Only effective if \code{addPieCharts = TRUE}.}
  \item{width}{Width of the plot. If \code{width=NULL} a default value is used.}
  \item{height}{Height of the plot. If \code{height=NULL} a default value is used.}
  \item{legend}{Boolean: if \code{TRUE} a legend is added to the plot.
        \code{x$effectNames} are used as labels.}
  \item{legendColumns}{Number of columns in legend. If \code{legendColumns=NULL}
        a default value is used. Only effective if \code{legend=TRUE}.}
  \item{legendHeight}{Height of legend. If \code{legendHeight=NULL} a default
        value is used. Only effective if \code{legend=TRUE}.}
  \item{cex.legend}{Specifies the relative font size of legend labels.}
  \item{cex.names}{Specifies the relative font size of bar graph labels.}
  \item{intervalsPerPeriod}{Integer: number of bins for time-series plots
        of dynamic RI. Default 20. Ignored for static RI.}
  \item{\dots}{Additional arguments passed to \code{\link{relativeImportance}}.}
}
\details{\code{interpret_size} works only for estimation by the Method of Moments
      using \code{\link{modelType}} 1 (directed networks) or 2
      ('forcing' model for non-directed networks).
      For dependent behavior variables, \code{behModelType=1} ('standard') is assumed.
      It does not yet work for endowment or creation effects
      (i.e., included effects have to be of type \code{eval}) or \code{rate},
    and also not for models with interaction effects.
    For two-mode (bipartite) networks as dependent variables,
    it works only if the number of second-mode nodes is less than the
    number of actors.

    If there are any missing tie values in the network data set, they are imputed by
    initial zeros and Last Observation Carried Forward. Structural zeros and ones
    are replaced by \code{NA} and treated as impossible choices in the probability
    vectors and ignored in the standard deviations; but the change statistics
    for these dyads still are given in \code{changeStatistics} (if requested).

    The averages reported in the components \code{sigmas} (average across actors)
    and \code{meansigmas} (average across waves) are obtained by averaging at the
    variance level and then taking square roots.
}
\value{
 If the model contains only one dependent variable, \code{interpret_size} returns
 an object of class \code{sienaRI}. Otherwise, it returns a list of objects
 of class \code{sienaRI}, each corresponding to one dependent variable.

 The \code{sienaFit} method returns a new-style \code{sienaRI} object with:
 \describe{
    \item{- \code{dependentVariable}}{name of the dependent variable.}
    \item{- \code{data}}{a flat data.frame (or data.table) with one row per
          ego-period, containing:
          columns \code{period}, \code{ego};
          \code{RI_<effectName>} columns with normalised relative importance;
          \code{absSum} (total L1 distance before normalisation);
          \code{sigma_<effectName>} (within-ego SD of change statistics).}
    \item{- \code{effectNames}}{names of the included evaluation effects.}
    \item{- \code{shortEffectNames}}{short (internal) effect names,
          matching the column-name suffixes in \code{data}.}
    \item{- \code{uncertainty}}{when \code{uncertainty = TRUE}, a
          data.frame with uncertainty estimates (SE, quantile columns)
          from \code{sienaPostestimate}; otherwise absent.}
 }
 The \code{print} and \code{plot} methods work with both the new format
 (from the \code{sienaFit} method) and the old nested-list format
 (from calling \code{sienaRI()} directly or via the \code{sienaEffects}
 method).

 The \code{sienaEffects} method returns the old nested-list format.
 A returned old-format \code{sienaRI} object stores the expected relative
 importances of effects of one dependent variable at observation moments as
 defined in Section 3.1 of Indlekofer and Brandes (2013).

Old-format \code{sienaRI} components:
    \describe{
    \item{- \code{dependentVariable}}{the name of the corresponding
          dependent variable.}
    \item{- \code{effectNames}}{the names of considered effects.}
    \item{- \code{RIActors}}{a list that contains the expected
          relative importances of effects
          for each potential actor decision at observation moments.
          This is equation (3) in Indlekofer and Brandes (2013).}
    \item{- \code{expectedRI}}{a list that contains the expected
          relative importances of effects aggregated over all actors for each
          network observation. These are the averages of the actor related
          values in \code{RIActors}.
          This is equation (4) in Indlekofer and Brandes (2013).}
    \item{- \code{IActors}}{a list that contains the expected
          importances of effects for each potential actor decision at
          observation moments. This is the numerator of equation (3)
          in Indlekofer and Brandes (2013).}
    \item{- \code{expectedI}}{a list that contains the expected
          importances of effects
          aggregated over all actors in each observation.
          More precisely, it contains the averages of the actor related values
          in \code{IActors}.}
    \item{- \code{absoluteSumActors}}{a list that contains the sum of
          the (unstandardized) L1-differences
          calculated for each potential actor decision at observation moments.
          This is the denominator of equation (3)
          in Indlekofer and Brandes (2013).}
    \item{- \code{RHActors}}{a list that contains the degree of certainty,
          also called degree of determination, in
          the potential ministep taken by an actor at the observation moments;
          this is R_H(i,x) of formula (6) in Snijders (2004).
          The mean over actors of these degrees of certainty, given by
          formula (7) in Snijders (2004), is printed by the
          \code{print} method for \code{interpret_size} objects.}
    \item{- \code{sigma}}{a list of effects by ego matrices of the values of
          the within-ego standard deviations of the change statistics.}
    \item{- \code{sigmas}}{effects by wave matrices of averages
         (see above for how this is done; across actors) of \code{sigma}.
          These are printed if \code{printSigma=TRUE}.}
    \item{- \code{meansigmas}}{average (see above for how this is done)
          over waves of \code{sigmas}.}
    \item{- \code{changeStatistics}}{a list of arrays
          (effects by choices by egos)
          containing, for each observation wave,
          the values of the change statistics for making this choice,
          where for one-mode networks the choice is defined by the alters
          with alter=ego referring to 'no change',
          for two-mode networks the choice is defined by the second-mode nodes
          with the last choice referring to 'no change',
          and for behavior the choice is defined as
          going down, staying constant, or going up;
          this output is produced only if \code{getChangeStats=TRUE}.}
    \item{- \code{toggleProbabilities}}{an array
         (egos by choices by waves), where "choices" are as directly above,
         giving the choice probabilities of ego in a ministep,
         when the data are as in this wave.}
    }
}

\references{
 Indlekofer, N. and Brandes, U. (2013),
 Relative Importance of Effects in Stochastic Actor-oriented Models.
 \emph{Network Science}, \bold{1}, 278--304.

 Snijders, T.A.B. (2004),
  Explained Variation in Dynamic Network Models.
  \emph{Mathematics and Social Sciences}, \bold{168}, 31--41.
}
\author{Natalie Indlekofer, Tom Snijders, Daniel Gotthardt}
\seealso{\code{\link{relativeImportance}}, \code{\link{interpret_size_dynamics}},
  \code{\link{sienaEntropy}}}
\examples{
myalgo <- set_algorithm_saom(nsub=1, n3=50, seed=1293)
mynet1 <- as_dependent_rsiena(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- make_data_rsiena(mynet1)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, list(density, recip, outAct, inPop))
(myeff <- set_effect(myeff, reciAct, parameter=1))
ans <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE)

## New format (via sienaFit method): data.frame-based RI
RI <- interpret_size(ans, mydata)
RI
print(RI, printSigma=TRUE)
# Access the underlying data.frame:
head(RI$data)
# RI columns per effect:
grep("^RI_", names(RI$data), value=TRUE)
# sigma columns per effect:
grep("^sigma_", names(RI$data), value=TRUE)
\dontrun{
plot(RI, addPieChart=TRUE)
plot(RI, actors=1:20, addPieChart=TRUE, radius=1.08)
}

## Legacy format (via sienaEffects method): list-based RI
theta.eval <- coef(ans, dropRates=TRUE)
RI_old <- interpret_size(myeff, data=mydata, theta=theta.eval)
RI_old
# average within-ego standard deviations of change statistics by wave:
RI_old$sigmas
# sigma averaged over waves:
RI_old$meansigmas
# to get the change statistics:
RI_old.cs <- interpret_size(myeff, data=mydata, theta=theta.eval,
                            getChangeStats=TRUE)
# For the network at the first wave:
dim(RI_old.cs$changeStatistics[[1]])
# semi-standardized parameters by wave:
theta.eval * RI_old$meansigmas

# For two dependent variables:
myalgo <- set_algorithm_saom(nsub=1, n3=50, seed=1293)
mynet2 <- as_dependent_rsiena(array(c(s502, s503), dim=c(50, 50, 2)))
mybeh <- as_dependent_rsiena(s50a[,2:3], type="behavior")
mydata2 <- make_data_rsiena(mynet2, mybeh)
myeff2 <- make_specification(mydata2)
myeff2 <- set_effect(myeff2, list(density, recip, transTies))
myeff2 <- set_effect(myeff2, avAlt, depvar="mybeh", covar1="mynet2")
ans2 <- siena(mydata2, effects=myeff2, control_algo=myalgo, batch=TRUE)
theta.eval2 <- coef(ans2, dropRates=TRUE)
RI2 <- interpret_size(myeff2, data=mydata2, theta=theta.eval2)
RI2[[1]]
RI2[[2]]
\dontrun{
plot(RI2[[2]], col = c("red", "green"), legend=FALSE)
plot(RI2[[1]], addPieChart = TRUE, legendColumns=2)
}
}
\keyword{models}
