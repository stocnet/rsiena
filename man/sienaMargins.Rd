\name{marginalEffects}
\alias{marginalEffects}
\alias{marginalEffects.sienaFit}
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
  }
}
\usage{
\method{marginalEffects}{sienaFit}(object, data, targets,
    control_uncertainty = set_postest_uncertainty_saom(),
    control_algo        = set_postest_algo_saom(),
    control_out         = set_postest_output_saom())
}
\arguments{
  \item{object}{A fitted \code{sienaFit} object, as returned by
      \code{\link{siena07}}.}
  \item{data}{The \code{siena} data object the model was fitted to.
      Multi-group data (\code{sienaGroup}) is not supported.}
  \item{targets}{A \code{sienaPostestTargets} object saying \emph{what} to
      compute: which effects to perturb, how, and the estimand they are
      evaluated under. Built with \code{\link{make_marginal_targets}} and
      refined with \code{\link{set_target}},
      \code{\link{set_second_diff}} and \code{\link{set_dependency}}.
      Required.}
  \item{control_uncertainty}{A \code{sienaPostestUncertainty} object from
      \code{\link{set_postest_uncertainty_saom}}: whether standard errors
      are computed and by which method.}
  \item{control_algo}{A \code{sienaPostestControl} object from
      \code{\link{set_postest_algo_saom}}: \emph{how} the quantities are
      computed -- chain simulation, batching, parallelism, memory, logging.}
  \item{control_out}{A \code{sienaPostestOutput} object from
      \code{\link{set_postest_output_saom}}: what the returned object
      looks like.}
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
    \item{Quantile columns}{Named by the requested \code{ciInterval}
      bounds (e.g.\ \code{q_025}, \code{q_975}), plus \code{Median} when
      \code{uncertaintyMedian = TRUE}.}
  }
}
\details{
  \strong{Effect name resolution}: targets are named by effect short name,
  as in \code{\link{set_effect}}. Internally these resolve to the
  composite names the engine uses, of the form
  \code{depvar_shortName[_covariate]_type} (e.g.\
  \code{"mynet_transTrip_eval"}). A short name that matches more than one
  effect -- \code{egoX} on several covariates, or two \code{unspInt}
  terms -- is an error rather than a silent choice of the first; see
  \code{\link{set_target}} for how \code{covar1}/\code{covar2}
  disambiguate.

  \strong{Interaction effects}: When perturbing an effect also moves
  another effect's change statistic -- an interaction term, say -- declare
  the relationship with \code{\link{set_dependency}}, e.g.
  \code{set_dependency(tg, egoXaltX ~ egoX:altX)}. The relation is
  evaluated at each counterfactual state, so it applies to both steps of a
  second difference and is not restricted to products of two effects.

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

## What to compute: one target, a unit step in transTrip's change statistic,
## reported as a tie probability averaged by period.
mytg <- make_marginal_targets(ans, effects = myeff, depvar = "mynet1",
                             type = "tieProb", level = "period",
                             includeDefaults = FALSE)
mytg <- set_target(mytg, transTrip, diff = 1)

me <- marginalEffects(ans, data = mydata, targets = mytg,
        control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE))
me

## The same perturbation on the change probability, unaggregated.  Only the
## estimand and its reporting change, so only the targets object does.
mytg_raw <- make_marginal_targets(ans, effects = myeff, depvar = "mynet1",
                                 type = "changeProb", level = "none",
                                 includeDefaults = FALSE)
mytg_raw <- set_target(mytg_raw, transTrip, diff = 1)
me_raw <- marginalEffects(ans, data = mydata, targets = mytg_raw,
            control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE))
head(me_raw)

\dontrun{
## With uncertainty: how precisely is a control_uncertainty question, so the
## targets object is unchanged.
me_unc <- marginalEffects(ans, data = mydata, targets = mytg,
            control_uncertainty = set_postest_uncertainty_saom(nsim = 100),
            control_algo = set_postest_algo_saom(verbose = FALSE))
me_unc

## Second difference: how does the marginal effect of transTrip change when
## recip also changes?  Each effect's perturbation is stated next to it.
mytg2 <- make_marginal_targets(ans, effects = myeff, depvar = "mynet1",
                              type = "tieProb", level = "period",
                              includeDefaults = FALSE)
mytg2 <- set_second_diff(mytg2, list(transTrip = list(diff = 1),
                                     recip     = list(diff = 1)))
me2 <- marginalEffects(ans, data = mydata, targets = mytg2,
         control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE))
me2
}
}
\keyword{models}
