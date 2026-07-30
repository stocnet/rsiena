\name{set_postest_uncertainty_saom}
\alias{set_postest_uncertainty_saom}
\alias{print.sienaPostestUncertainty}
\alias{sienaPostestUncertainty}
\title{Configure uncertainty quantification for SAOM postestimation}
\description{
  \code{set_postest_uncertainty_saom} constructs a
  \code{sienaPostestUncertainty} object bundling the settings that control
  how \code{\link{marginalEffects}} (and, prospectively,
  \code{\link{predict.sienaFit}}) quantify uncertainty of postestimation
  quantities: whether uncertainty is computed at all, the method used
  (bootstrap-style theta draws vs.\ analytic delta method), how many draws
  to take, and which summaries (standard error, confidence interval, mean,
  median) to report.
}
\usage{
set_postest_uncertainty_saom(
    enabled    = TRUE,
    mode       = c("delta", "bootstrap", "deltaFull"),
    nsim       = 1000,
    sd         = TRUE,
    ci         = TRUE,
    ciInterval = c(0.025, 0.975),
    simMean    = FALSE,
    simMedian  = FALSE
)

\method{print}{sienaPostestUncertainty}(x, \dots)
}
\arguments{
  \item{enabled}{Logical: whether uncertainty is computed at all. If
      \code{FALSE}, only point estimates are returned.}
  \item{mode}{How standard errors are derived. \code{"delta"} (the default)
      is analytic and deterministic; \code{"bootstrap"} draws \code{nsim}
      parameter vectors; \code{"deltaFull"} adds the path-distribution term
      to the delta SE and requires a dynamic targets object.

      The default was \code{"bootstrap"} before RSiena's postestimation
      redesign. It changed because the out-of-the-box call was two orders of
      magnitude more expensive than it needed to be -- around 114 seconds
      against 0.57 for \code{"delta"} on a 50-node example -- and stochastic
      as well. Standard errors stay on by default; only the way they are
      obtained changed.

      For a dynamic model \code{"delta"} gives the conditional
      (frozen-chain) standard error, which omits the path-distribution term.
      That term is a covariance and may have either sign, so the conditional
      SE is not a bound on the full one; use \code{"deltaFull"} when the
      chain distribution is part of the question.}
  \item{nsim}{Single positive numeric, coerced to integer: the number of
      theta draws used when \code{mode = "bootstrap"}. Ignored (but still
      validated) for the delta modes.}
  \item{sd}{Logical: include the standard error / Monte Carlo SE of the
      estimate in the output.}
  \item{ci}{Logical: include a confidence interval in the output.}
  \item{ciInterval}{Numeric vector of length 2 giving the lower and upper
      quantile bounds of the confidence interval, e.g.\
      \code{c(0.025, 0.975)} for a 95\% CI. Both values must be strictly
      between 0 and 1, and \code{ciInterval[1] < ciInterval[2]}.}
  \item{simMean}{Logical: include the mean across simulation draws in the
      output. Only meaningful when \code{mode = "bootstrap"}; setting this
      to \code{TRUE} together with \code{mode = "delta"} or
      \code{mode = "deltaFull"} is an error, since those modes produce no
      draws to average.}
  \item{simMedian}{As \code{simMean}, for the median across simulation
      draws.}
  \item{x}{A \code{sienaPostestUncertainty} object, as returned by
      \code{set_postest_uncertainty_saom}.}
  \item{\dots}{Further arguments, currently ignored (present for S3
      consistency with \code{\link{print}}).}
}
\value{
  \code{set_postest_uncertainty_saom} returns a list of class
  \code{"sienaPostestUncertainty"} with elements \code{enabled},
  \code{mode}, \code{nsim}, \code{sd}, \code{ci}, \code{ciInterval},
  \code{simMean}, and \code{simMedian}, and a \code{"version"} attribute
  giving the RSiena package version used to construct it.

  \code{print.sienaPostestUncertainty} returns \code{x}, invisibly; it is
  called for its side effect of printing a human-readable summary.
}
\details{
  A \code{sienaPostestUncertainty} object is passed to
  \code{\link{marginalEffects}} to select the uncertainty method used for
  postestimation quantities and to control which summaries of that
  uncertainty are reported.

  \strong{Modes}: \code{mode = "bootstrap"} draws \code{nsim} values of
  theta from \eqn{N(\hat\theta, \hat\Sigma_\theta)} and re-evaluates the
  postestimation quantity for each draw, producing an empirical
  distribution of simulation draws. \code{mode = "delta"} and
  \code{mode = "deltaFull"} instead propagate \eqn{\hat\Sigma_\theta}
  analytically via the delta method; they are typically much cheaper than
  bootstrapping but produce no draws, so \code{simMean} and
  \code{simMedian} — which summarise the draw distribution — are not
  available under these modes, and requesting them is an error.

  \code{ciInterval} gives the lower and upper quantile bounds of the
  confidence interval (e.g.\ \code{c(0.025, 0.975)} for a 95\% CI); both
  bounds must lie strictly between 0 and 1 and be given in increasing
  order.
}
\seealso{
  \code{\link{set_postest_algo_saom}} for the accompanying simulation /
  parallelisation control object, and \code{\link{marginalEffects}} for
  the function these settings configure.
}
\examples{
## Defaults: bootstrap uncertainty with SE and 95\% CI
u1 <- set_postest_uncertainty_saom()
print(u1)

## Analytic delta method, SE only, no CI
u2 <- set_postest_uncertainty_saom(mode = "delta", ci = FALSE)
print(u2)

## No uncertainty quantification at all (point estimates only)
u3 <- set_postest_uncertainty_saom(enabled = FALSE)
print(u3)

## Bootstrap with a 90\% CI and both mean and median of the draws
u4 <- set_postest_uncertainty_saom(
    mode = "bootstrap", nsim = 500,
    ciInterval = c(0.05, 0.95),
    simMean = TRUE, simMedian = TRUE
)
print(u4)
}
\keyword{models}
