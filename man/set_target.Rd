\name{set_target}
\alias{set_target}
\alias{set_target.sienaPostestTargets}
\alias{set_second_diff}
\alias{set_second_diff.sienaPostestTargets}
\alias{set_dependency}
\alias{set_dependency.sienaPostestTargets}
\title{Select, tune, and relate postestimation targets}
\description{
  Three functions modify a \code{sienaPostestTargets} object, as built by
  \code{\link{make_postest_targets}}. Each returns a modified \emph{copy}
  of its input; none of them mutate the object in place.

  \code{set_target} selects or deselects one or more existing targets and
  tunes how each one is perturbed, optionally overriding the model-level
  settings from \code{\link{make_postest_targets}} for just that target.

  \code{set_second_diff} adds a new target crossing two effects: the first
  difference of one effect after also perturbing a second, i.e.\ a second
  difference (see \code{\link{marginalEffects}}).

  \code{set_dependency} declares that one effect's change statistic moves
  together with others', so that perturbing a component effect also
  perturbs the effects that depend on it.
}
\usage{
set_target(x, \dots)

\method{set_target}{sienaPostestTargets}(x, shortNames,
    diff = NULL, contrast = NULL,
    covar1 = NULL, covar2 = NULL,
    perturbType = NULL,
    level, condition, accumulated,
    rateWeight, massContrasts,
    include = TRUE, verbose = TRUE, \dots)

set_second_diff(x, \dots)

\method{set_second_diff}{sienaPostestTargets}(x, perturbations,
    returnDecisionDetails = FALSE,
    name = NULL,
    include = TRUE, verbose = TRUE, \dots)

set_dependency(x, \dots)

\method{set_dependency}{sienaPostestTargets}(x, \dots, verbose = TRUE)
}
\arguments{
  \item{x}{A \code{sienaPostestTargets} object, as returned by
      \code{\link{make_postest_targets}} or by an earlier
      \code{set_target}, \code{set_second_diff}, or \code{set_dependency}
      call.}
  \item{shortNames}{One or more effect short names identifying which
      target(s) to change, given as bare symbols or characters, as in
      \code{\link{set_effect}} (e.g.\ \code{transTrip},
      \code{"transTrip"}, or \code{c(density, recip)}). Any number of
      names may be given, and the same change is applied to every one of
      them.}
  \item{perturbations}{A named list saying what is done to each of the two
      effects being crossed: names are effect short names, values are lists
      using the same settings \code{set_target} takes for a single
      effect -- \code{diff}, \code{contrast}, \code{covar1},
      \code{covar2}, \code{perturbType}. For example

      \preformatted{list(transTrip = list(diff = 1),
     recip = list(contrast = c(0, 1)))}

      A \code{NULL} value, or an empty list, gives that effect its default
      perturbation (see \code{\link{make_postest_targets}}). Exactly two
      entries are required; the same effect may appear twice, which crosses
      it with itself.

      This replaced a numbered form -- two bare effect names plus
      \code{diff1}/\code{contrast2} and friends -- which is no longer
      accepted; a call in that shape is an error naming the list form
      rather than being quietly ignored.}
  \item{diff}{Single numeric step size for the target's first-difference
      perturbation. Mutually exclusive with \code{contrast}.}
  \item{contrast}{Numeric vector of length 2 giving the two change-statistic
      values to contrast directly, instead of stepping by a fixed amount.
      This is the default for effects whose change statistic is not on a
      scale a step means anything on: \code{c(-1, 1)} for \code{density},
      which only ever takes those two values, and \code{c(0, 1)} for
      \code{recip}, which records whether the reverse tie is present.
      Mutually exclusive with \code{diff}.}
  \item{covar1, covar2}{Optional covariate name(s) identifying \emph{which}
      target a short name refers to, named and used exactly as in
      \code{\link{set_effect}}. A short name does not always identify one
      target: a model may carry \code{egoX} on several covariates, or two
      \code{unspInt} terms. A short name matching more than one target is
      an error rather than a silent choice of the first, and
      \code{covar1}/\code{covar2} (or the qualified short name, e.g.\
      \code{egoX_gender} or \code{unspInt2}, which is accepted in place of
      \code{shortNames}) resolve it. Only meaningful when a single short
      name is given.

      \code{set_effect}'s other identifying arguments do not apply here:
      \code{depvar} is fixed for a whole targets object, and
      \code{type} likewise, since \code{eval}/\code{endow}/\code{creation}
      variants of an effect are deliberately a single target (perturbing
      them is one perturbation, not three).}
  \item{perturbType}{Optional override of the perturbation mode:
      \code{"alter"} (one-alternative logit update) or \code{"ego"}
      (ego-wide renormalisation). \code{NULL} (the default) leaves the
      auto-detection used by \code{\link{marginalEffects}} in place; see
      that page for details.}
  \item{level}{Optional per-target override of the \code{level} set in
      \code{\link{make_postest_targets}}, for this target only. If not
      supplied, the target inherits the model-level value (see Details).}
  \item{condition}{Optional per-target override of the \code{condition}
      set in \code{\link{make_postest_targets}}, for this target only. An
      explicit \code{condition = NULL} is different from not supplying
      \code{condition} at all: it switches conditioning off for this
      target even when the model sets one (see Details).}
  \item{accumulated}{Optional per-target override of the \code{accumulated}
      set in \code{\link{make_postest_targets}}, for this target only (see
      Details for how it combines with the model-level value).}
  \item{rateWeight}{Optional per-target override of the \code{rateWeight}
      set in \code{\link{make_postest_targets}}, for this target only (see
      Details for how it combines with the model-level value).}
  \item{massContrasts}{Optional per-target override of the
      \code{massContrasts} set in \code{\link{make_postest_targets}}, for
      this target only (see Details).}
  \item{include}{Logical, default \code{TRUE}: for \code{set_target},
      whether the named target(s) are selected for computation; for
      \code{set_second_diff}, whether the newly added target is selected.
      \code{FALSE} records the target without selecting it.}
  \item{verbose}{Logical: print a message describing what changed.}
  \item{name}{Optional character name for the new target row. Defaults to
      the two effects' names joined with \code{"_"} and suffixed
      \code{"_sd"}. An error is raised if a target with this name already
      exists.}
  \item{returnDecisionDetails}{Logical, default \code{FALSE}: return the
      per-decision detail columns for this target.}
  \item{\dots}{For \code{set_target} and \code{set_second_diff}, further
      arguments (currently unused, present for S3 consistency). For
      \code{set_dependency}, one or more one-sided dependency formulas of
      the form \code{target ~ a:b} (e.g.\
      \code{transRecTrip ~ transTrip:recip}); at least one must be
      given. See Details.}
}
\value{
  All three functions return a modified copy of \code{x}, of class
  \code{c("sienaPostestTargets", "data.frame")}.

  \code{set_target} changes the \code{include}, \code{diff1}/
  \code{contrast1}, \code{perturbType1} and \code{.overrides} columns of
  the rows named in \code{shortNames}; it never adds or removes rows.

  \code{set_second_diff} returns a copy of \code{x} with one row appended
  for the new cross-partial target.

  \code{set_dependency} returns a copy of \code{x} whose
  \code{"dependencies"} attribute has the new formula(s) appended; no rows
  are added, changed, or removed.

  Unless \code{verbose = FALSE}, all three also print a message describing
  the change.
}
\details{
  \strong{Per-target overrides} (\code{level}, \code{condition},
  \code{accumulated}, \code{rateWeight}, \code{massContrasts} in
  \code{set_target}) let one fitted model be evaluated at several
  conditioning sets or aggregation levels in a single
  \code{\link{marginalEffects}} call, rather than requiring one call per
  setting. This matters in particular for a dynamic model: the simulated
  chains are shared across all targets within one call, so two
  conditioning sets evaluated this way cost one simulation, versus two
  simulations if done via two separate targets objects and two separate
  calls.

  Only a field actually named in the \code{set_target} call becomes an
  override; an argument left missing means the target inherits the
  model-level value, which is different from explicitly passing
  \code{NULL}. This distinction matters most for \code{condition}: an
  explicit \code{condition = NULL} on a target switches conditioning off
  for that target even when \code{\link{make_postest_targets}} set one,
  whereas simply not supplying \code{condition} inherits the model-level
  setting. When the targets object is used, \code{level} and
  \code{condition} overrides replace the model-level value outright
  (including an explicit \code{NULL}); \code{accumulated} and
  \code{rateWeight} overrides are combined with the model-level value by
  logical \code{OR} (so a target can turn one of these on when the model
  default is off, but cannot turn it off when the model default is on);
  and \code{massContrasts} replaces the model-level value whenever the
  override is non-\code{NULL}.

  \strong{Second differences} (\code{set_second_diff}): each of the two
  components defaults its own perturbation independently, the same way a
  single \code{set_target} call would (see
  \code{\link{make_postest_targets}}), unless that component's list entry
  says otherwise.

  \strong{Dependencies} (\code{set_dependency}): a dependency is a
  one-sided formula \code{target ~ a:b} meaning that \code{target}'s
  change statistic is the product of \code{a}'s and \code{b}'s change
  statistics. Perturbing \code{a} by \code{delta} therefore also moves
  \code{target}'s change statistic by \code{delta * s_b} (and
  symmetrically when \code{b} is the one perturbed). Only multiplicative
  relationships of exactly this form are supported; anything else (e.g.\
  \code{+}, \code{==}, or a function call on the right-hand side) is
  rejected with an error rather than silently approximated.

  A dependency only helps when the component effect being perturbed (
  \code{a} or \code{b}) is itself an included effect in the fitted
  model's \code{effects} object, since it is resolved by reusing that
  effect's own change statistic; it cannot supply a change statistic for
  a component that was never estimated, and does not, for instance,
  reconstruct the components of a \code{gwesp}-type effect or of
  \code{transRecTrip} out of thin air.

  A dependency is model-level: once declared, it applies regardless of
  which target ends up perturbing \code{a} or \code{b}, and it applies to
  both steps of a second difference. Dependencies accumulate across
  repeated \code{set_dependency} calls.
}
\seealso{
  \code{\link{make_postest_targets}} for building the targets object,
  \code{\link{marginalEffects}} for the function \code{targets} is passed
  to, and \code{\link{set_postest_uncertainty_saom}} /
  \code{\link{set_postest_algo_saom}} for the accompanying uncertainty and
  simulation control objects.
}
\examples{
\dontrun{
## ans is a fitted sienaFit; mymodel is the sienaEffects object it was
## fitted with; mydata is the corresponding siena data object.
mytg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                             includeDefaults = FALSE)
mytg <- set_target(mytg, transTrip, diff = 2)
mytg <- set_target(mytg, c(density, recip))
mytg <- set_second_diff(mytg, list(density   = list(contrast = c(-1, 1)),
                                   transTrip = list(diff = 1)))
mytg <- set_dependency(mytg, transRecTrip ~ transTrip:recip)
print(mytg)

me <- marginalEffects(ans, data = mydata, targets = mytg,
                      control_uncertainty = set_postest_uncertainty_saom(),
                      control_algo = set_postest_algo_saom())

## Per-target overrides: evaluate two different effects, each conditioned
## on a different variable, in a single call. For a dynamic model the
## simulated chains are shared across both targets, instead of being
## simulated once per call for two separate calls.
mytg2 <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                              includeDefaults = FALSE)
mytg2 <- set_target(mytg2, transTrip, diff = 1, condition = "recip")
mytg2 <- set_target(mytg2, density, contrast = c(-1, 1), condition = NULL)
me2 <- marginalEffects(ans, data = mydata, targets = mytg2,
                       control_uncertainty = set_postest_uncertainty_saom(),
                       control_algo = set_postest_algo_saom())
}
}
\keyword{models}
