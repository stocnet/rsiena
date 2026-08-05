\name{predict.sienaFit}
\alias{predict.sienaFit}
\alias{fitted.sienaFit}
\alias{residuals.sienaFit}
\title{Predicted Probabilities from a Stochastic Actor-Oriented Model}
\description{
  Compute predicted probabilities of tie changes (or tie states) from a
  fitted SAOM.  The method uses the estimated (or user-supplied)
  parameters together with the change contributions extracted from the
  data to produce actor-level multinomial choice probabilities.

  Two probability types are available:
  \describe{
    \item{\code{changeProb}}{The probability that an actor changes the
      tie to each alter (or, for behaviour, makes each possible choice)
      in a single ministep.  These are the raw multinomial logit
      (softmax) probabilities.}
    \item{\code{tieProb}}{The probability of a tie being present
      \emph{after} the ministep.  For alters currently tied,
      \code{tieProb = 1 - changeProb} (the tie is kept unless the actor
      chooses to dissolve it); for alters currently not tied,
      \code{tieProb = changeProb}.}
  }

  Static predictions use change contributions computed at the
  observation moments.  Dynamic predictions simulate Markov chains
  between observations using \code{siena07} and produce predictions at
  each simulated ministep.

  When \code{uncertainty = TRUE} (the default), parameters are drawn
  from the asymptotic distribution
  \eqn{N(\hat\theta, \hat\Sigma_\theta)} and predictions are
  aggregated across draws to produce point estimates with Monte Carlo
  standard errors and confidence intervals.
}
\usage{
\method{predict}{sienaFit}(object, data, targets,
    control_uncertainty = set_postest_uncertainty_saom(),
    control_algo        = set_postest_algo_saom(),
    control_out         = set_postest_output_saom(),
    newParams = NULL, \dots)
}
\arguments{
  \item{object}{A \code{sienaFit} object.}
  \item{data}{A \code{siena} data object (typically the same data used
      for estimation).}
  \item{targets}{A prediction-targets object from
      \code{\link{make_predict_targets}}, optionally extended with
      \code{\link{set_condition}}.  It states WHAT to predict: the
      quantity (\code{type}), whether to evaluate at the observed waves or
      along simulated chains (\code{dynamic}), how far to aggregate
      (\code{level}), and what to stratify by (\code{condition}).  One
      prediction is returned per row.}
  \item{control_uncertainty}{From
      \code{\link{set_postest_uncertainty_saom}}.

      For a STATIC prediction the decision contexts are fixed, so
      \code{mode = "delta"} is the whole of the uncertainty and agrees with
      \code{"bootstrap"} to within a few percent at roughly one hundredth of
      the cost.

      For a DYNAMIC prediction the contexts are simulated and therefore
      themselves depend on the parameters.  \code{"delta"} then reports the
      CONDITIONAL (frozen-chain) standard error, which omits that channel and
      is correspondingly smaller -- three to four times smaller on the s50
      example.  \code{"deltaFull"} adds the path-distribution term and is the
      one comparable with \code{"bootstrap"}: on the same example the two
      agree to within Monte-Carlo error (ratios 1.15 and 0.86 at 300 draws),
      at 7 seconds against 31 minutes.  A \code{deltaFull} run reports both,
      as \code{SE} and \code{SE_conditional}.

      Both analytic modes are estimated from the simulated chains, so
      \code{n3PointEst} governs their accuracy: at very small chain counts
      the frozen-chain Jacobian is unreliable.}
  \item{control_algo}{From \code{\link{set_postest_algo_saom}}.  Only
      consulted when the targets object sets \code{dynamic = TRUE}, which
      requires an \code{algorithm}.}
  \item{control_out}{From \code{\link{set_postest_output_saom}}.  Also
      carries \code{returnDecisionDetails} and \code{returnComponents}.}
  \item{newParams}{Optional named numeric vector of parameter values.  If
      \code{NULL}, the estimates from \code{object} are used.}
  \item{\dots}{Additional arguments (currently unused).}
}
\details{
  \strong{Static predictions} (\code{dynamic = FALSE}, the default)
  extract change contributions at each observation moment from the
  data and C++ internals.  These contributions are computed once and
  reused across all theta draws, making static predictions fast.

  \strong{Dynamic predictions} (\code{dynamic = TRUE}) simulate
  Markov chains between observation moments using \code{siena07}.
  Each theta draw triggers a fresh set of chain simulations, so
  dynamic predictions with uncertainty are substantially slower.

  \strong{Level and condition are separate axes.}  \code{level} sets the
  grain: \code{"period"} averages within a period, \code{"ego"} gives one
  row per actor, and \code{"egoChoice"} (static) or
  \code{"ministepChoice"} (dynamic) keep one row per alternative -- the
  ego's whole tie profile rather than a summary of it.  \code{condition}
  stratifies within whatever grain is chosen, by the value of a change
  statistic; conditioning on density separates tie creation from tie
  maintenance.  A profile is a grain, not a stratification, and the two can
  be combined or used alone.

  Multi-group data (\code{sienaGroup}) is not yet supported.
}
\value{
  With one row in \code{targets}, a data.frame; with several, a named list
  of them, named by the targets' \code{name} column.  Each carries class
  \code{"sienaPrediction"}.

  Columns depend on \code{level} and \code{condition}: the grouping
  variables for the requested level (\code{period}, \code{ego},
  \code{choice}, \code{chain}, \code{ministep} as applicable), one column
  per conditioning statistic, and the predicted quantity named after
  \code{type} (\code{changeProb} or \code{tieProb}).

  When uncertainty is enabled the output additionally carries \code{SE} and
  the requested interval columns.  The standard error is always called
  \code{SE} whatever produced it; \code{attr(x, "uncertaintyMethod")}
  records the method, since two standard errors that hold different things
  fixed should not share a name.

  \code{returnDecisionDetails = TRUE} attaches the per-decision detail as an
  attribute rather than adding columns.
}
\section{Choosing an uncertainty mode}{
  Measured on the \code{s50} example; treat the numbers as indicative, not as
  general constants.

  \strong{Standard errors.}  \code{"delta"} and \code{"bootstrap"} agree
  closely: within 2\% on predicted probabilities and first differences, and
  within 9\% on a second difference, against a 4000-draw bootstrap.  If a
  standard error is what you report, the analytic modes are enough.

  \strong{Intervals are a different matter.}  \code{"delta"} reports
  \code{est} \eqn{\pm} 1.96 SE, which is symmetric by construction, whereas
  the bootstrap reports empirical quantiles.  For a SECOND DIFFERENCE the
  sampling distribution is noticeably skewed -- upper and lower half-widths in
  a ratio of about 0.53 to 0.71 on this example -- and the symmetric interval
  placed the upper bound 25--39\% too high even though the standard errors
  agreed.  Use \code{"bootstrap"} when an interval on a second difference is
  the reported quantity.

  \strong{Dynamic quantities.}  \code{"delta"} gives the conditional
  (frozen-chain) standard error, which omits the path-distribution channel;
  \code{"deltaFull"} adds it and is the one comparable with a bootstrap.  Two
  chain settings do different jobs and are easy to confuse:
  \code{n3PointEst} sizes the POINT ESTIMATE simulation and has no effect on
  the standard error; \code{n3} sizes the UNCERTAINTY simulation and governs
  its stability -- about 7\% relative variation across chain realisations at
  \code{n3 = 200}, 5\% at 1000, 3\% at 2000, with cost linear in
  \code{n3}.  For publishable dynamic results use \code{"deltaFull"} with
  \code{n3} of at least 500, and \code{n3PointEst} of 1000 or more for a
  stable point estimate.
}
\author{Daniel Gotthardt}
\seealso{
  \code{\link{make_predict_targets}}, \code{\link{set_condition}},
  \code{\link{marginalEffects}}, \code{\link{interpret_size}},
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

## One prediction, aggregated by period, no uncertainty.
ptg <- make_predict_targets(ans, effects=myeff, type="tieProb",
                            level="period")
predict(ans, mydata, targets=ptg,
        control_uncertainty=set_postest_uncertainty_saom(enabled=FALSE))

## Several predictions in one call: the ego tie profile, and the same
## quantity split into tie creation versus tie maintenance.
ptg <- make_predict_targets(ans, effects=myeff, type="tieProb",
                            level="period", includeDefaults=FALSE)
ptg <- set_condition(ptg, level="egoChoice", name="profile")
ptg <- set_condition(ptg, density, name="by_context")
ptg

out <- predict(ans, mydata, targets=ptg,
               control_uncertainty=set_postest_uncertainty_saom(enabled=FALSE))
head(out$profile)
out$by_context
}
\keyword{models}
