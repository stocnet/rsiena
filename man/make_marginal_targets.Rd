\name{make_marginal_targets}
\alias{make_marginal_targets}
\alias{make_marginal_targets.sienaFit}
\alias{print.sienaPostestTargets}
\alias{sienaPostestTargets}
\title{Build a postestimation targets specification}
\description{
  \code{make_marginal_targets} builds a \code{sienaPostestTargets} object: a
  table listing which marginal effects to compute, and how each one is
  perturbed. It is built and then modified the same way an effects
  specification is, via \code{\link{set_target}}, \code{\link{set_second_diff}},
  and \code{\link{set_dependency}}, and is then passed to
  \code{\link{marginalEffects}} as its \code{targets} argument.
}
\usage{
\method{make_marginal_targets}{sienaFit}(x, data = NULL, effects = NULL,
    depvar = NULL,
    dynamic = FALSE,
    type = c("changeProb", "tieProb"),
    mainEffect = c("riskDifference", "riskRatio"),
    level = "period",
    condition = NULL,
    egoNormalize = TRUE,
    accumulated = FALSE,
    rateWeight = FALSE,
    massContrasts = NULL,
    dependencies = NULL,
    includeDefaults = TRUE, \dots)

\method{print}{sienaPostestTargets}(x, \dots)
}
\arguments{
  \item{x}{For \code{make_marginal_targets}, a fitted \code{sienaFit} object
      (as returned by \code{\link{siena07}}) to build targets for. For
      \code{print.sienaPostestTargets}, the \code{sienaPostestTargets}
      object to print, as returned by \code{make_marginal_targets} or by a
      later \code{\link{set_target}}/\code{\link{set_second_diff}}/
      \code{\link{set_dependency}} call.}
  \item{data}{Currently unused (reserved).}
  \item{effects}{The \code{sienaEffects} object the model was fitted with,
      e.g.\ \code{effects = mymodel}. Required: this is the full effects
      object, not \code{x$requestedEffects}, because interaction
      components used by \code{\link{set_target}} and
      \code{\link{set_dependency}} live in it as non-included rows.}
  \item{depvar}{Name of the dependent variable the targets are enumerated
      for. Defaults to the first dependent variable named in
      \code{effects}.}
  \item{dynamic}{Logical: if \code{TRUE}, the marginal effect is averaged
      over Markov chains simulated between observations (via
      \code{\link{siena07}}); if \code{FALSE} (the default) it is
      evaluated on the change contributions at the observation moments.
      These are different estimands, not two routes to the same number,
      which is why the choice is made here rather than among the
      computation settings in \code{\link{set_postest_algo_saom}} --
      although \code{dynamic = TRUE} does require an \code{algorithm}
      to be supplied there. \code{accumulated = TRUE} requires it.}
  \item{type}{Character: \code{"changeProb"} (change probability, the
      default) or \code{"tieProb"} (tie probability). Model-level default
      for every target; see \code{\link{marginalEffects}} for what the two
      types mean.}
  \item{mainEffect}{Character: \code{"riskDifference"} (the default) or
      \code{"riskRatio"} — the scale the marginal effect is reported on.
      Model-level default for every target.}
  \item{level}{Aggregation level for the output, e.g.\ \code{"period"}
      (the default). Model-level default for every target; see
      \code{\link{marginalEffects}} for the available levels. May be
      overridden per target via \code{\link{set_target}} (see Details).}
  \item{condition}{Optional character vector of effect short names to
      condition the aggregation on. Model-level default for every target;
      \code{NULL} (the default) means no conditioning. May be overridden
      per target via \code{\link{set_target}} (see Details).}
  \item{egoNormalize}{Logical: normalise contributions per ego (default
      \code{TRUE}). Model-level default for every target.}
  \item{accumulated}{Logical: accumulate predictions across ministeps
      (default \code{FALSE}). Model-level default for every target; may be
      overridden per target via \code{\link{set_target}} (see Details).}
  \item{rateWeight}{Logical: weight by the rate function (default
      \code{FALSE}). Model-level default for every target; may be
      overridden per target via \code{\link{set_target}} (see Details).}
  \item{massContrasts}{Logical or \code{NULL} (the default): whether to
      compute mass contrasts (\code{massCreation}, \code{massDissolution}).
      \code{NULL} auto-detects from the perturbation mode, as in
      \code{\link{marginalEffects}}. Model-level default for every target;
      may be overridden per target via \code{\link{set_target}} (see
      Details).}
  \item{dependencies}{Optional list of one-sided dependency formulas, of
      the kind added by \code{\link{set_dependency}} (e.g.\
      \code{list(transRecTrip ~ transTrip:recip)}). Ordinarily left
      \code{NULL} here and built up afterwards with
      \code{\link{set_dependency}}; see that page for the declaration
      syntax.}
  \item{includeDefaults}{Logical: if \code{TRUE} (the default), every
      non-rate effect in \code{effects} is included as a selected target
      with its default perturbation (\code{diff = 1}, or
      \code{contrast = c(-1, 1)} for \code{density}). If \code{FALSE},
      every effect is listed but none is selected
      (\code{include = FALSE} for every row); use this to build the
      target set up from nothing with \code{\link{set_target}} rather
      than narrow a full default set down.}
  \item{\dots}{Further arguments passed to methods (currently unused for
      \code{print.sienaPostestTargets}, present for S3 consistency with
      \code{\link{print}}).}
}
\value{
  \code{make_marginal_targets} returns an object of class
  \code{c("sienaPostestTargets", "data.frame")} with one row per candidate
  effect (rate effects are excluded outright, see Details). Columns
  include \code{name}, \code{effectName1}, \code{effectType},
  \code{include}, \code{second}, \code{effectName2}, \code{diff1},
  \code{contrast1}, \code{diff2}, \code{contrast2}, \code{perturbType1},
  \code{perturbType2}, and \code{.overrides} (a list-column of per-target
  overrides set via \code{\link{set_target}}). The model-level defaults,
  the dependent variable name, and the \code{effects} object used to build
  the targets are carried as attributes, along with a \code{"version"}
  attribute giving the RSiena package version used to construct it.

  \code{\link{set_target}}, \code{\link{set_second_diff}}, and
  \code{\link{set_dependency}} each return a modified \emph{copy} of a
  \code{sienaPostestTargets} object; none of them mutate their input.

  \code{print.sienaPostestTargets} returns \code{x}, invisibly. It prints
  the estimand and reporting settings the targets share, any declared
  dependencies, and then one line per selected target -- name, effect(s)
  perturbed, what is done to each change statistic, and whether it is a
  first or a second difference -- in the order the results are returned
  in. Settings a target overrides are shown in brackets on its line.
}
\details{
  Only non-rate effects are candidate targets: rate effects have no
  meaningful first difference in this framework and are dropped rather
  than listed-but-unselected. Creation/endow/eval variants of the same
  effect are collapsed to a single target, since perturbing them is one
  perturbation, not three.

  \strong{Default perturbation}: when \code{includeDefaults = TRUE},
  every candidate effect gets a default first-difference perturbation.
  Effects whose change statistic is not on a scale a unit step means
  anything on default to a contrast between its two possible values:
  \code{density}'s only ever takes \eqn{-1} and \eqn{+1}, so it defaults
  to \code{contrast = c(-1, 1)}, and \code{recip}'s records whether the
  reverse tie is present, so it defaults to \code{contrast = c(0, 1)}.
  Every other effect defaults to \code{diff = 1}, a unit step.

  \strong{Per-target overrides}: \code{level}, \code{condition},
  \code{accumulated}, \code{rateWeight}, and \code{massContrasts} are set
  once here as defaults for the whole model, but \code{\link{set_target}}
  can override any of them for a single target. This matters for
  evaluating one fitted model under several conditioning sets or
  aggregation levels in a single \code{\link{marginalEffects}} call: for a
  dynamic model, the simulated chains are shared across all targets within
  one call, so covering two conditioning sets this way costs one
  simulation, versus two simulations if the same comparison is made with
  two separate targets objects and two separate calls.
}
\seealso{
  \code{\link{set_target}} and \code{\link{set_second_diff}} for
  selecting, tuning, and adding targets, \code{\link{set_dependency}} for
  declaring how effects relate, \code{\link{marginalEffects}} for the
  function \code{targets} is passed to, and
  \code{\link{set_postest_uncertainty_saom}} /
  \code{\link{set_postest_algo_saom}} for the accompanying uncertainty and
  simulation control objects.
}
\examples{
\dontrun{
## ans is a fitted sienaFit; mymodel is the sienaEffects object it was
## fitted with; mydata is the corresponding siena data object.
mytg <- make_marginal_targets(ans, effects = mymodel, depvar = "mynet",
                             type = "tieProb", level = "period")
mytg <- set_target(mytg, transTrip, diff = 2)
mytg <- set_second_diff(mytg, list(density = NULL, recip = list(diff = 1)))
print(mytg)

me <- marginalEffects(ans, data = mydata, targets = mytg,
                      control_uncertainty = set_postest_uncertainty_saom(),
                      control_algo = set_postest_algo_saom())

## Start from nothing and build the target set up explicitly
mytg2 <- make_marginal_targets(ans, effects = mymodel, depvar = "mynet",
                              includeDefaults = FALSE)
mytg2 <- set_target(mytg2, c(density, transTrip))
}
}
\keyword{models}
