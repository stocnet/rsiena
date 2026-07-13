\name{utilitySpec}
\alias{utilitySpec}
\alias{print.utilitySpec}
\title{Spec-based utility framework for structured marginal effects}
\description{
  \code{utilitySpec} declares interaction decompositions between SAOM
  change statistics using R formula syntax, for use with the marginal
  effects / counterfactual utility machinery. Density is implicit
  (never declared in the spec) and always interacts with every other
  effect.
}
\usage{
utilitySpec(..., contrasts = list(), diffs = list())
}
\arguments{
  \item{\dots}{One or more two-sided formulas of the form
    \code{intEffectName ~ factor1 : factor2}, or a single list of such
    formulas. The left-hand side names the SIENA interaction effect;
    the right-hand side gives the change-statistic factors that
    combine to form it, separated by \code{:}. Wrapping a factor name
    in \code{factor()} marks it as a binary factor (default contrast
    \code{c(0, 1)}); bare names are treated as continuous (default
    perturbation \code{diff = 1}).}
  \item{contrasts}{Named list overriding the default \code{c(0, 1)}
    contrast for specific factor terms.}
  \item{diffs}{Named list overriding the default perturbation size
    (\code{1}) for specific continuous terms.}
}
\details{
  Change statistics follow the \code{attachContributions(flip = TRUE)}
  convention: the density column is \eqn{\pm 1} (the sign of the
  proposed tie change), and all other columns are unsigned change
  statistics \eqn{\delta_k}. The utility for dyad \code{ij} is
  \deqn{U_{ij} = \delta_{density} \times (\theta_{density} +
    \sum_{k \neq density} \delta_k \theta_k),}
  equivalently \code{contribMat \%*\% theta}, since
  \eqn{contrib_k = \delta_{density} \times \delta_k}.

  Because density is implicit, perturbing it directly gives
  \code{utilDiff = -2 * U_ij}; perturbing another effect \code{k} by
  \code{Delta} gives
  \code{utilDiff = delta_density * Delta * (theta_k + interaction terms)}.
  State updates are then simply \code{delta[k] <- delta[k] + Delta},
  with no density multiplication needed.

  For example, \code{utilitySpec(transRecTrip ~ factor(recip) : transTrip)}
  maps \eqn{\delta_{recip} \times \delta_{transTrip}} to the SIENA
  effect \code{transRecTrip}.
}
\value{
  An object of class \code{"utilitySpec"} describing the declared
  interaction decompositions, for use internally by the marginal
  effects / counterfactual prediction machinery
  (see \code{\link{marginalEffects}}).
}
\seealso{
  \code{\link{marginalEffects}}
}
\keyword{models}
