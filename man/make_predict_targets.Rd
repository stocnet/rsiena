\name{make_predict_targets}
\alias{make_predict_targets}
\alias{make_predict_targets.sienaFit}
\alias{set_condition}
\alias{set_condition.sienaPredictTargets}
\alias{print.sienaPredictTargets}
\title{Prediction Targets for a Stochastic Actor-Oriented Model}
\description{
  Builds the object that says WHAT \code{\link{predict.sienaFit}} should
  compute.  It is the sibling of \code{\link{make_marginal_targets}} and is
  much smaller: a marginal-effect target names an effect to perturb and how,
  whereas a prediction has nothing to perturb, so a prediction target names
  only what to condition on and at what level.

  One row per requested prediction.  \code{set_condition} appends a row;
  several stratifications of the same prediction therefore come back from a
  single \code{predict} call.
}
\usage{
\method{make_predict_targets}{sienaFit}(x, data = NULL, effects = NULL,
    depvar = NULL, type = c("changeProb", "tieProb"),
    dynamic = FALSE, level = "period", condition = NULL,
    egoNormalize = TRUE, accumulated = FALSE, rateWeight = FALSE,
    na.rm = TRUE, includeDefaults = TRUE, \dots)

\method{set_condition}{sienaPredictTargets}(x, condition, level,
    egoNormalize, accumulated, rateWeight,
    name = NULL, verbose = TRUE, \dots)
}
\arguments{
  \item{x}{For \code{make_predict_targets}, a \code{sienaFit} object; for
      \code{set_condition}, a prediction-targets object.}
  \item{data}{Optional \code{siena} data object.}
  \item{effects}{The \code{sienaEffects} object the model was fitted with.
      Required: the fit does carry one, but conditional estimation and
      stripped effect objects make the implicit route unreliable, so the
      caller states it.}
  \item{depvar}{Name of the dependent variable.  Defaults to the first in
      \code{effects}.}
  \item{type}{\code{"changeProb"}, the probability that ego changes the tie,
      or \code{"tieProb"}, the probability that the tie exists after the
      ministep.  The two differ in sign for maintenance opportunities.}
  \item{dynamic}{Evaluate along simulated chains rather than at the observed
      waves.  Requires an \code{algorithm} on the algorithm control.}
  \item{level}{How far to aggregate: \code{"none"}, \code{"period"},
      \code{"ego"}, \code{"egoChoice"}, \code{"chain"}, \code{"chainEgo"},
      \code{"ministep"}, \code{"ministepChoice"}.  \code{"egoChoice"} and its
      dynamic counterpart \code{"ministepChoice"} keep one row per
      alternative, which is the ego's tie profile rather than a summary of
      it.}
  \item{condition}{Change statistic(s) to stratify by, given as bare short
      names.  Conditioning on \code{density} separates tie creation from tie
      maintenance.  This is a separate axis from \code{level}: a profile is a
      grain, not a stratification.}
  \item{egoNormalize}{Average within each realised ego decision before
      aggregating further.}
  \item{accumulated}{Sum over the ministeps of a chain rather than averaging
      over them.  Requires \code{dynamic = TRUE}.}
  \item{rateWeight}{Weight by the rate parameter.  Absorbed by the simulation
      when \code{dynamic = TRUE}, and a message says so.}
  \item{na.rm}{Passed to the aggregation.}
  \item{includeDefaults}{If \code{TRUE} (default) the object starts with one
      row carrying the model-level request.  \code{FALSE} starts empty, for
      when every prediction wanted is stated explicitly -- otherwise the
      unconditional row rides along and has to be ignored downstream.}
  \item{name}{Optional label for the row, used to name the corresponding
      element of \code{predict}'s output.  Defaults to
      \code{"by_<condition>"}, or \code{"unconditional"} when no condition is
      given.}
  \item{verbose}{Whether to report what was added.}
  \item{\dots}{Additional arguments (currently unused).}
}
\details{
  Settings given to \code{make_predict_targets} are model-level defaults;
  \code{set_condition} may override \code{level}, \code{egoNormalize},
  \code{accumulated} and \code{rateWeight} for a single row.  An argument
  that is not supplied inherits the default, which is distinct from
  supplying \code{NULL} -- an explicit \code{condition = NULL} asks for the
  unconditional prediction even when the model sets one.

  The estimand rules are shared with \code{\link{make_marginal_targets}}, so
  the two entry points cannot drift apart: \code{accumulated} requires
  \code{dynamic}, and \code{rateWeight} is absorbed by a dynamic run.
}
\value{
  An object of class \code{c("sienaPredictTargets", "data.frame")}, with one
  row per requested prediction.  The model-level defaults, the dependent
  variable and the \code{effects} object are carried as attributes.
}
\author{Daniel Gotthardt}
\seealso{
  \code{\link{predict.sienaFit}}, \code{\link{make_marginal_targets}},
  \code{\link{set_postest_uncertainty_saom}}
}
\examples{
myalgo <- set_algorithm_saom(nsub=1, n3=50, seed=1293)
mynet1 <- as_dependent_rsiena(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- make_data_rsiena(mynet1)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, list(density, recip, transTrip))
ans <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE)

## The ego tie profile, and tie creation versus tie maintenance, together.
ptg <- make_predict_targets(ans, effects=myeff, type="tieProb",
                            level="period", includeDefaults=FALSE)
ptg <- set_condition(ptg, level="egoChoice", name="profile")
ptg <- set_condition(ptg, density, name="by_context")
ptg
}
\keyword{models}
