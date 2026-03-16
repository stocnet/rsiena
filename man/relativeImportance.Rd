\name{relativeImportance}
\alias{relativeImportance}
\alias{relativeImportance.sienaFit}
\alias{relativeImportance.sienaEffects}
\alias{sienaRI}
\title{Relative importance of effects in a Stochastic Actor-Oriented Model}
\description{
  Computes the relative importance (RI) of effects for a fitted SAOM,
  following Section 3.1 of Indlekofer and Brandes (2013).

  For each actor, RI is derived from the L1 distance between the
  full-\eqn{\theta} choice distribution and the leave-one-out
  distribution obtained by zeroing each effect in turn.  Normalising
  across effects gives the relative importance.

  \code{relativeImportance} is the core computation with methods for
  both \code{sienaFit} (full model with optional uncertainty and dynamic
  RI) and \code{sienaEffects} (point-estimate only, requiring a
  user-supplied \eqn{\theta}).
  \code{\link{interpret_size}} is the user-facing orchestrator that
  calls \code{relativeImportance} and may attach additional measures
  (entropy, etc.).

  Both static (at observation moments) and dynamic (along simulated
  ministep chains) RI are supported through the \code{dynamic} argument.
  The latter replaces the deprecated \code{\link{interpret_size_dynamics}}.
}
\usage{
relativeImportance(object, \dots)

\method{relativeImportance}{sienaFit}(object, data,
    depvar = NULL, effects = NULL,
    distFun = c("L1D", "KLD"),
    uncertainty = FALSE,
    dynamic = FALSE, algorithm = NULL, n3 = 500,
    useChangeContributions = FALSE,
    dynamicMinistepFactor = 10,
    level = "ego", sum_fun = mean, na.rm = TRUE,
    nsim = 1000,
    uncertaintyMode = c("batch", "stream"),
    uncertaintySd = TRUE, uncertaintyCi = TRUE,
    uncertaintyProbs = c(0.025, 0.5, 0.975),
    useCluster = FALSE, nbrNodes = 1,
    clusterType = c("PSOCK", "FORK"),
    batchDir = "temp", batchSize = NULL,
    keepBatch = FALSE, verbose = TRUE,
    memoryScale = NULL, batchUnitBudget = 2.5e8,
    getChangeStats = FALSE, \dots)

\method{relativeImportance}{sienaEffects}(object, data,
    theta, depvar = NULL, getChangeStats = FALSE,
    distFun = c("L1D", "KLD"), \dots)
}
\arguments{
  \item{object}{A \code{sienaFit} object resulting from a call to
        \code{\link{siena}}, or a \code{sienaEffects} object.}
  \item{data}{\code{siena} data object representing the analysed data,
        resulting from a call to \code{\link{make_data_rsiena}}.}
  \item{depvar}{Character: name of the dependent variable.
        If \code{NULL}, all dependent variables are processed.}
  \item{theta}{Numeric vector of parameter values for evaluation effects
        (required for the \code{sienaEffects} method).  Length must equal
        the number of included evaluation effects.}
  \item{effects}{A \code{sienaEffects} object.
        If \code{NULL}, uses \code{object$requestedEffects}.}
  \item{distFun}{Character: distance function used to measure the
        difference between full and leave-one-out distributions.
        \code{"L1D"} (default) is the L1 distance;
        \code{"KLD"} is the Kullback-Leibler divergence.}
  \item{uncertainty}{Logical.
        If \code{TRUE}, uncertainty quantification via draws from
        \eqn{N(\hat\theta, \Sigma_\theta)} is computed using
        \code{sienaPostestimate}.}
  \item{dynamic}{Logical. If \code{TRUE}, compute RI along simulated
        dynamic ministep chains rather than at observation moments.
        Requires \code{algorithm} to be specified.
        This mode replaces the deprecated \code{\link{sienaRIDynamics}} and 
        the intermediate \code{\link{interpret_size_dynamics}}.}
  \item{algorithm}{\code{sienaAlgorithm} object used for the dynamic
        simulation. Required when \code{dynamic = TRUE}.}
  \item{n3}{Integer: number of chains to simulate for dynamic RI
        (default 500).}
  \item{useChangeContributions}{Logical: if \code{TRUE}, reuse
        change contributions already stored in \code{object} instead
        of recomputing them.}
  \item{dynamicMinistepFactor}{Numeric: scaling factor controlling the
        number of ministeps per chain in dynamic mode (default 10).}
  \item{level}{Aggregation level: \code{"ego"} (one row per ego-period,
        default), \code{"period"} (aggregated over egos), or
        \code{"none"} (all rows).}
  \item{sum_fun}{Aggregation function applied within groups when
        \code{level} aggregates (default \code{mean}).}
  \item{na.rm}{Logical: whether to remove \code{NA}s in aggregation.}
  \item{nsim}{Number of theta draws for uncertainty quantification
        (default 1000).}
  \item{uncertaintyMode}{Character: \code{"batch"} (default) or
        \code{"stream"}.}
  \item{uncertaintySd}{Logical: include Monte Carlo standard errors
        in uncertainty output?}
  \item{uncertaintyCi}{Logical: include confidence intervals?}
  \item{uncertaintyProbs}{Numeric vector of quantile probabilities for
        confidence intervals (default \code{c(0.025, 0.5, 0.975)}).}
  \item{useCluster}{Logical: use a parallel cluster for uncertainty
        computation?}
  \item{nbrNodes}{Number of cluster nodes.}
  \item{clusterType}{Type of cluster: \code{"PSOCK"} or \code{"FORK"}.}
  \item{batchDir}{Directory for temporary batch files.}
  \item{batchSize}{Batch size; \code{NULL} (default) for automatic sizing.}
  \item{keepBatch}{Logical: retain batch files after computation?}
  \item{verbose}{Logical: print progress messages?}
  \item{memoryScale}{Optional memory scaling factor.}
  \item{batchUnitBudget}{Memory budget for automatic batch sizing.}
  \item{getChangeStats}{Logical: if \code{TRUE}, attach the change
        contribution struct (wide format) as a \code{changeStatistics}
        component of the returned object.}
  \item{\dots}{Additional arguments passed to \code{sienaPostestimate}.}
}
\details{
  \bold{Method and applicability.}
  Works only for estimation by the Method of Moments using
  \code{\link{modelType}} 1 (directed networks) or 2 (forcing model for
  non-directed networks).  For dependent behavior variables,
  \code{behModelType = 1} (standard) is assumed.
  Does not work for endowment or creation effects (included effects must
  be of type \code{eval} or \code{rate}), or for models with interaction
  effects.  For two-mode (bipartite) networks as dependent variables,
  it works only if the number of second-mode nodes is less than the
  number of actors.

  \bold{Static RI (\code{dynamic = FALSE}).}
  Change contributions are computed analytically at each observation
  moment.  For each ego, the multinomial choice probabilities are
  computed under the full parameter vector and under each
  leave-one-out alternative (\eqn{\theta_k = 0}).  The L1 distance
  between the full and leave-one-out distributions gives the
  (unnormalised) importance \eqn{I_k(i)}.  Normalising:
  \deqn{RI_k(i) = \frac{I_k(i)}{\sum_{k'} I_{k'}(i)}}
  (equation (3) of Indlekofer and Brandes, 2013).
  The averages \eqn{\overline{RI}_k} over actors (equation (4)) are
  also stored.

  \bold{Dynamic RI (\code{dynamic = TRUE}).}
  Uses \code{\link{getDynamicChangeContributions}} to simulate ministep
  chains and computes per-ministep RI values, capturing how effect
  importance evolves over the inter-observation period rather than only
  at the endpoints.

  \bold{Missing values.}
  If there are missing tie values in the network data, they are imputed by
  initial zeros and Last Observation Carried Forward.  Structural zeros
  and ones are replaced by \code{NA} and treated as impossible choices
  in the probability vectors (ignored in the standard deviations);
  the change statistics for these dyads are still included in
  \code{changeStatistics} if \code{getChangeStats = TRUE}.

  \bold{Standard deviations.}
  The components \code{sigma} (average across actors) and
  \code{meansigmas} (average across waves) are obtained by averaging at
  the variance level and then taking square roots.  Multiplying
  \code{meansigmas} by the corresponding parameter values gives
  semi-standardised parameters.
}
\value{
  If the model contains one dependent variable, returns a \code{sienaRI}
  object.  Otherwise, a named list of \code{sienaRI} objects.

  \bold{New format} (from \code{relativeImportance.sienaFit} and
  \code{interpret_size.sienaFit}):
  \describe{
    \item{\code{dependentVariable}}{Name of the dependent variable.}
    \item{\code{data}}{A flat data.frame (or data.table) with one row per
          ego-period, containing columns \code{period}, \code{ego};
          \code{RI_<effectName>} columns with normalised relative
          importance; \code{absSum} (total L1 distance before
          normalisation); and \code{sigma_<effectName>} columns
          (within-ego SD of change statistics).}
    \item{\code{effectNames}}{Names of the included evaluation effects.}
    \item{\code{shortEffectNames}}{Short (internal) effect names,
          matching the column-name suffixes in \code{data}.}
    \item{\code{uncertainty}}{When \code{uncertainty = TRUE}, a data.frame
          with uncertainty estimates (SE, quantile columns) from
          \code{sienaPostestimate}; otherwise absent.}
    \item{\code{changeStatistics}}{When \code{getChangeStats = TRUE}, the
          change contribution wide struct.}
  }

  \bold{Old format} (from \code{interpret_size.sienaEffects} / legacy
  \code{sienaRI()} calls):
  \describe{
    \item{\code{dependentVariable}}{Name of the dependent variable.}
    \item{\code{effectNames}}{Names of considered effects.}
    \item{\code{RIActors}}{List of per-actor expected relative importances
          at observation moments (equation (3) of Indlekofer and Brandes,
          2013).}
    \item{\code{expectedRI}}{List of expected relative importances
          aggregated over all actors per observation
          (equation (4) of Indlekofer and Brandes, 2013).}
    \item{\code{IActors}}{List of (unnormalised) expected importances per
          actor (numerator of equation (3)).}
    \item{\code{expectedI}}{List of importances averaged over actors.}
    \item{\code{absoluteSumActors}}{List of total L1 distances per actor
          (denominator of equation (3)).}
    \item{\code{RHActors}}{List of degrees of certainty per actor;
          see \code{\link{entropy}} for definition and documentation.}
    \item{\code{sigma}}{List of effects-by-ego matrices of within-ego
          standard deviations of change statistics.}
    \item{\code{sigmas}}{Effects-by-wave matrices of averages of
          \code{sigma} across actors (printed when
          \code{printSigma = TRUE}).}
    \item{\code{meansigmas}}{Average of \code{sigmas} over waves.}
    \item{\code{changeStatistics}}{List of arrays
          (effects by choices by egos) of change statistics per
          observation wave; produced only when
          \code{getChangeStats = TRUE}.  For one-mode networks, choices
          are alters with alter\eqn{=}ego referring to no change; for
          two-mode networks, choices are second-mode nodes with the last
          entry referring to no change; for behavior, choices are
          down, stay, or up.}
    \item{\code{toggleProbabilities}}{Array (egos by choices by waves)
          of multinomial choice probabilities under the full parameter
          vector, where choices are as defined for
          \code{changeStatistics} above.}
  }
}
\references{
  Indlekofer, N. and Brandes, U. (2013),
  Relative Importance of Effects in Stochastic Actor-oriented Models.
  \emph{Network Science}, \bold{1}, 278--304.
}
\author{Natalie Indlekofer, Tom A.B. Snijders, Daniel Gotthardt}
\seealso{\code{\link{interpret_size}}, \code{\link{interpret_size_dynamics}},
  \code{\link{sienaEntropy}}}
\examples{
myalgo <- set_algorithm_saom(nsub=1, n3=50, seed=1293)
mynet1 <- as_dependent_rsiena(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- make_data_rsiena(mynet1)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, list(density, recip, outAct, inPop))
ans <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE)

## Core computation (same as interpret_size for sienaFit):
RI <- relativeImportance(ans, mydata)
RI
head(RI$data)
grep("^RI_", names(RI$data), value=TRUE)

## With uncertainty quantification:
\dontrun{
RI_unc <- relativeImportance(ans, mydata, uncertainty=TRUE, nsim=500)
RI_unc$uncertainty
}

## Dynamic RI (along simulated chains):
\dontrun{
myalgo2 <- set_algorithm_saom(nsub=1, n3=100, seed=42, cond=FALSE)
RI_dyn <- relativeImportance(ans, mydata, dynamic=TRUE, algorithm=myalgo2)
plot(RI_dyn)
}
}
\keyword{models}
