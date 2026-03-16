\name{marginalEffects}
\alias{marginalEffects}
\alias{marginalEffects.sienaFit}
\alias{marginalEffects.sienaEffects}
\title{Marginal Effects for SAOMs}
\description{
  \code{marginalEffects} computes marginal effects (MEs) for
  Stochastic Actor-oriented Models.  An ME is the expected change in
  actor-level choice probabilities that results from a hypothetical
  unit (or user-specified) change in an effect's change statistic, averaged
  across all actors and choice alternatives.

  \strong{First differences} (\code{second = FALSE}, the default)
  measure the effect of a change in one predictor on the probability
  scale, analogous to average marginal effects in logistic regression.
  For a change \eqn{\Delta x} in effect \eqn{k}:
  \deqn{\text{AME}_k = E_i \left[
    P(y_i = j \mid x_k + \Delta x) - P(y_i = j \mid x_k)
  \right]}

  \strong{Second differences} (\code{second = TRUE}) measure how the
  first difference of one effect changes when another effect is also
  modified, capturing interaction or moderation on the probability scale.

  Two probability types are supported:
  \describe{
    \item{\code{changeProb}}{Change probability: the probability of
      toggling/changing a tie or behaviour choice in a ministep.}
    \item{\code{tieProb}}{Tie probability: the probability that the tie
      is present \emph{after} the ministep (accounts for whether the
      current state is tied or untied).}
  }

  Behaviour dependent variables are not currently supported and will
  trigger an error.  See the design document
  (\code{docs/design_count_behavior.md}) for planned extensions.

  \code{marginalEffects} is a generic function that dispatches on the class
  of its first argument:
  \itemize{
    \item \code{marginalEffects.sienaFit}: Full functionality with
          uncertainty quantification via theta draws from
          \eqn{N(\hat\theta, \hat\Sigma_\theta)}.
    \item \code{marginalEffects.sienaEffects}: Takes a bare effects object
          with a user-supplied \code{theta} vector and optional
          \code{covtheta} for uncertainty.
  }
}
\usage{
\method{marginalEffects}{sienaFit}(object, data, effects = NULL, depvar = NULL,
    effectName1, diff1 = NULL, contrast1 = NULL,
    interaction1 = FALSE, intEffectNames1 = NULL,
    modEffectNames1 = NULL,
    effectName2 = NULL, diff2 = NULL, contrast2 = NULL,
    interaction2 = FALSE, intEffectNames2 = NULL,
    modEffectNames2 = NULL,
    perturbType1 = NULL, perturbType2 = NULL,
    type = c("changeProb", "tieProb"), second = FALSE,
    level = "period", condition = NULL, sum_fun = mean, na.rm = TRUE,
    dynamic = FALSE, algorithm = NULL, n3 = 500,
    useChangeContributions = FALSE, uncertainty = TRUE,
    uncertaintyMode = c("batch", "stream"), nsim = 1000,
    uncertaintySd = TRUE, uncertaintyCi = TRUE,
    uncertaintyProbs = c(0.025, 0.5, 0.975),
    uncertaintyMcse = FALSE, uncertaintymcseBatches = NULL,
    useCluster = FALSE, nbrNodes = 1,
    clusterType = c("PSOCK", "FORK"), cluster = NULL,
    batchDir = "temp", prefix = "simBatch_b",
    combineBatch = TRUE, batchSize = NULL,
    keepBatch = FALSE, verbose = TRUE,
    mainEffect = "riskDifference", details = FALSE,
    memoryScale = NULL, batchUnitBudget = 2.5e8,
    dynamicMinistepFactor = 10)

\method{marginalEffects}{sienaEffects}(object, data, theta,
    covtheta = NULL, depvar = NULL, effectName1,
    uncertainty = !is.null(covtheta), ...)
}
\arguments{
  \item{object}{For \code{marginalEffects.sienaFit}: a \code{sienaFit} object
           resulting from a call to \code{\link{siena}}.
           For \code{marginalEffects.sienaEffects}: a \code{sienaEffects}
           object.}
  \item{data}{\code{siena} data object.}
  \item{effects}{Optional \code{sienaEffects} object. If \code{NULL},
      taken from \code{x$requestedEffects}.}
  \item{theta}{Parameter vector (required for the \code{sienaEffects} method).}
  \item{covtheta}{Covariance matrix of the parameter estimates (optional
      for the \code{sienaEffects} method). When provided, \code{uncertainty}
      defaults to \code{TRUE} and uncertainty quantification is performed.}
  \item{depvar}{Name of dependent variable. Default: first depvar in data.}
  \item{effectName1}{Name of the effect for which the first difference
      is calculated.  This should match a short effect name as used
      internally (e.g.\ \code{"transTrip"}, \code{"recip"}).
      Partial matching is supported: the function matches against the
      composite effect names in the contributions
      (e.g.\ \code{"mynet_transTrip_eval"}).  It is sufficient to
      supply the short name.}
  \item{diff1}{Numeric: size of the hypothetical change in the effect
      variable.  Default \code{NULL} uses unit change (\eqn{\pm 1}).
      For continuous effects this controls the step size of the
      discrete difference.}
  \item{contrast1}{Numeric vector of length 2: contrast values for the
      first effect.  Used when you want to evaluate the first difference
      at specific values rather than a unit change.  The marginal effect
      is computed as \eqn{P(y \mid c_2) - P(y \mid c_1)}.}
  \item{interaction1}{Logical: whether \code{effectName1} is part of an
      interaction.}
  \item{intEffectNames1}{Character: names of interaction partner effects.}
  \item{modEffectNames1}{Character: names of moderation effects.}
  \item{effectName2}{For second differences: name of the second effect.}
  \item{diff2, contrast2, interaction2, intEffectNames2, modEffectNames2}{
      As above, for the second effect.}
  \item{perturbType1}{Perturbation mode for the first effect.
      \code{NULL} (default) auto-detects from the effect's
      \code{interactionType}: ego-type effects use \code{"ego"}
      (ego-wide renormalization), all others use \code{"alter"}
      (one-alternative logit update).  Can be set to \code{"alter"}
      or \code{"ego"} explicitly to override auto-detection.}
  \item{perturbType2}{As \code{perturbType1}, for the second effect
      (used only when \code{second = TRUE}).}
  \item{type}{Type of probability: \code{"changeProb"} (change probability)
      or \code{"tieProb"} (tie probability).}
  \item{second}{Logical: if \code{TRUE}, compute second differences.}
  \item{level}{Aggregation level for the output: \code{"none"}
      (unaggregated, one row per ego-choice), \code{"period"}
      (averaged across egos within period), \code{"ego"}
      (averaged across choices within ego), \code{"egoChoice"},
      \code{"chain"}, \code{"ministep"}, or \code{"ministepChoice"}.
      Default \code{"period"}.}
  \item{condition}{Optional conditioning variable for aggregation.}
  \item{sum_fun}{Aggregation function. Default: \code{mean}.}
  \item{na.rm}{Logical: remove NAs before aggregation.}
  \item{dynamic}{Logical: if \code{TRUE}, use dynamic contributions from
      simulated chains.}
  \item{algorithm}{Required when \code{dynamic = TRUE}: a
      \code{sienaAlgorithm} object.}
  \item{n3}{Number of chains for dynamic simulation.}
  \item{useChangeContributions}{Use pre-stored change contributions from
      \code{ans}.}
  \item{uncertainty}{Logical: whether to compute uncertainty via
      theta-draws.  For the \code{sienaEffects} method, defaults to
      \code{TRUE} when \code{covtheta} is supplied.}
  \item{uncertaintyMode}{\code{"batch"} or \code{"stream"}.}
  \item{nsim}{Number of theta-draws for uncertainty.}
  \item{uncertaintySd, uncertaintyCi, uncertaintyProbs, uncertaintyMcse,
      uncertaintymcseBatches}{Control output of uncertainty summaries.}
  \item{useCluster, nbrNodes, clusterType, cluster}{Parallelisation settings.}
  \item{batchDir, prefix, combineBatch, batchSize, keepBatch}{Batch mode
      settings.}
  \item{verbose}{Logical: print progress messages.}
  \item{mainEffect}{\code{"riskDifference"} or \code{"riskRatio"}.}
  \item{details}{Logical: if \code{TRUE}, include detailed intermediate
      quantities.}
  \item{memoryScale, batchUnitBudget, dynamicMinistepFactor}{Memory management.}
  \item{\dots}{Additional arguments passed to internal functions.}
}
\value{
  A data.frame (or data.table if available) with the marginal effect
  estimates.  The exact columns depend on \code{level}, \code{type},
  and \code{uncertainty}:

  \describe{
    \item{Coordinate columns}{Depending on \code{level}:
      \code{period}, \code{ego}, \code{choice} (static) or
      \code{chain}, \code{ministep} (dynamic).}
    \item{\code{firstDiff} or \code{secondDiff}}{The marginal effect
      on the probability scale (when
      \code{mainEffect = "riskDifference"}).}
    \item{\code{firstRiskRatio} or \code{secondRiskRatio}}{The marginal
      effect as a risk ratio (when \code{mainEffect = "riskRatio"}).}
    \item{Effect contribution columns}{One column per effect, giving
      the change statistic value for each ego-choice pair.}
  }

  When \code{uncertainty = TRUE}, the output is aggregated at the
  requested \code{level} with additional columns:
  \describe{
    \item{\code{Mean}}{Point estimate (mean across theta draws).}
    \item{\code{SE}}{Monte Carlo standard error of the mean.}
    \item{Quantile columns}{Named by the requested
      \code{uncertaintyProbs} (e.g.\ \code{2.5\%}, \code{50\%},
      \code{97.5\%}).}
  }
}
\details{
  \strong{Effect name resolution}: The \code{effectName1} (and
  \code{effectName2} for second differences) argument is matched
  against the composite effect names used internally, which have the
  form \code{depvar_shortName_type} (e.g.\
  \code{"mynet_transTrip_eval"}).  It is sufficient to supply just the
  short name (e.g.\ \code{"transTrip"}) — the function resolves it
  unambiguously.

  \strong{Interaction effects}: When the effect of interest is an
  interaction (e.g.\ specified via \code{set_interaction}), set
  \code{interaction1 = TRUE} and supply the names of the interacting
  effects in \code{intEffectNames1}.

  \strong{Perturbation mode (ego vs.\ alter)}: For most effects the
  change statistic varies across alters and the standard
  one-alternative logit update is appropriate (\code{perturbType =
  "alter"}).  For ego-type effects (e.g.\ \code{egoX}) the change
  statistic is constant across alters, so instead the ego-wide
  multinomial logit is renormalised:
  \deqn{p_j' = \frac{p_j \exp(\Delta u)}{\sum_k p_k \exp(\Delta u_k)}}
  This is auto-detected from the \code{interactionType} column of
  the effects object.  If the auto-detection gives an undesirable
  result, the user can override it via \code{perturbType1} /
  \code{perturbType2}.

  When the ego-wide perturbation mode is used, the output includes
  two additional columns:
  \describe{
    \item{\code{massCreation}}{Creation probability-mass contrast
      \eqn{\Delta P_i^+ = \sum_{j \in C_i} (p'_{ij} - p_{ij})},
      summing change-probability first differences over the creation
      risk set for each ego and period.}
    \item{\code{massDissolution}}{Dissolution probability-mass contrast
      \eqn{\Delta P_i^- = \sum_{j \in D_i} (p'_{ij} - p_{ij})},
      summing change-probability first differences over the dissolution
      risk set for each ego and period.}
  }
  These mass contrasts are always on the change-probability scale
  (toggle probabilities), regardless of the \code{type} argument.

  \strong{Static vs.\ dynamic}: With \code{dynamic = FALSE} (the
  default), change contributions are extracted at the observation
  moments and predictions are fast.  With \code{dynamic = TRUE},
  Markov chains are simulated between observations using
  \code{siena07}, and each theta draw triggers a fresh simulation.
  Dynamic predictions are therefore substantially slower.

  \strong{Aggregation}: The \code{level} parameter controls
  granularity.  For reporting purposes, \code{level = "period"} is
  typical: it gives one estimate per observation period, averaged
  across all egos and alternatives.

  Multi-group data (\code{sienaGroup}) is not yet supported.
}
\references{
  Gotthardt, D. and Steglich, C.E.G.. (2026).
  Marginal effects for the stochastic actor-oriented model.
}
\author{Daniel Gotthardt}
\seealso{
  \code{\link{predict.sienaFit}}, \code{\link{interpret_size}},
  \code{\link{relativeImportance}}, \code{\link{sienaEntropy}},
  \code{\link{getChangeStatistics}}
}
\examples{
myalgo <- set_algorithm_saom(nsub=1, n3=50, seed=1293)
mynet1 <- as_dependent_rsiena(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- make_data_rsiena(mynet1)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, list(density, recip, transTrip))
ans <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE)

## First difference of transTrip on tie probability, averaged by period
me <- marginalEffects(ans, data=mydata, effectName1="transTrip",
        diff1=1, type="tieProb", depvar="mynet1",
        level="period", uncertainty=FALSE)
me

## First difference on change probability, unaggregated
me_raw <- marginalEffects(ans, data=mydata, effectName1="transTrip",
            diff1=1, type="changeProb", depvar="mynet1",
            level="none", uncertainty=FALSE)
head(me_raw)

\dontrun{
## With uncertainty quantification
me_unc <- marginalEffects(ans, data=mydata, effectName1="transTrip",
            diff1=1, type="tieProb", depvar="mynet1",
            level="period", nsim=100, verbose=FALSE)
me_unc

## Second difference: how does the marginal effect of transTrip
## change when recip also changes?
me2 <- marginalEffects(ans, data=mydata,
         effectName1="transTrip", diff1=1,
         effectName2="recip", diff2=1,
         type="tieProb", depvar="mynet1",
         second=TRUE, level="period",
         uncertainty=FALSE)
me2
}
}
\keyword{models}
