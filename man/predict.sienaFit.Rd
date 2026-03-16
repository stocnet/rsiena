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
\method{predict}{sienaFit}(object, newdata, type = c("changeProb", "tieProb"),
    newParams = NULL, effects = NULL, depvar = NULL,
    dynamic = FALSE, algorithm = NULL, n3 = 1000,
    useChangeContributions = FALSE,
    level = "period", condition = NULL, sum_fun = mean, na.rm = TRUE,
    uncertainty = TRUE, uncertaintyMode = c("batch", "stream"),
    useCluster = FALSE, nbrNodes = 1,
    nsim = 1000, uncertaintySd = TRUE, uncertaintyCi = TRUE,
    uncertaintyProbs = c(0.025, 0.5, 0.975),
    uncertaintyMcse = FALSE, uncertaintymcseBatches = NULL,
    clusterType = c("PSOCK", "FORK"), cluster = NULL,
    batchDir = "temp", prefix = "simBatch_b",
    combineBatch = TRUE, batch = TRUE, silent = NULL,
    batchSize = NULL, keepBatch = FALSE, verbose = TRUE,
    memoryScale = NULL, batchUnitBudget = 5e8,
    dynamicMinistepFactor = 10, \dots)
}
\arguments{
  \item{object}{A \code{sienaFit} object.}
  \item{newdata}{A \code{siena} data object (typically the same data
      used for estimation).}
  \item{type}{\code{"changeProb"} (default) or \code{"tieProb"}.}
  \item{newParams}{Optional named numeric vector of parameter values.
      If \code{NULL}, the estimates from \code{object} are used.}
  \item{effects}{Optional \code{sienaEffects} object.  Defaults to
      \code{object$requestedEffects}.}
  \item{depvar}{Name of the dependent variable.  Defaults to the first
      dependent variable in \code{newdata}.}
  \item{dynamic}{Logical.  If \code{TRUE}, probabilities are computed
      along simulated Markov chains between observations.}
  \item{algorithm}{A \code{sienaAlgorithm} object.  Required when
      \code{dynamic = TRUE}.}
  \item{n3}{Number of simulated chains for dynamic predictions
      (default 1000).}
  \item{useChangeContributions}{Logical.  If \code{TRUE}, reuse
      change contributions stored in \code{object} (from the
      estimation run) instead of re-extracting them.}
  \item{level}{Aggregation level: \code{"none"}, \code{"period"},
      \code{"ego"}, \code{"egoChoice"}, \code{"chain"},
      \code{"ministep"}, or \code{"ministepChoice"}.
      Default \code{"period"}.}
  \item{condition}{Optional conditioning variable for grouped
      aggregation.}
  \item{sum_fun}{Aggregation function (default \code{mean}).}
  \item{na.rm}{Logical: remove NAs before aggregation.}
  \item{uncertainty}{Logical (default \code{TRUE}).  If \code{TRUE},
      compute uncertainty via Monte Carlo theta draws.}
  \item{uncertaintyMode}{\code{"batch"} (default) or \code{"stream"}.}
  \item{useCluster}{Logical: parallelize theta draws across a cluster.}
  \item{nbrNodes}{Number of cluster nodes.}
  \item{nsim}{Number of theta draws (default 1000).}
  \item{uncertaintySd, uncertaintyCi, uncertaintyProbs,
      uncertaintyMcse, uncertaintymcseBatches}{Control output columns
      for SE, confidence intervals, quantiles, and Monte Carlo SE.}
  \item{clusterType}{\code{"PSOCK"} or \code{"FORK"}.}
  \item{cluster}{An existing cluster object (optional).}
  \item{batchDir, prefix, combineBatch, batchSize, keepBatch}{Batch
      mode settings for large-scale computations.}
  \item{batch}{Logical: run \code{siena07} in batch mode (default
      \code{TRUE}).  Applies only when \code{dynamic = TRUE}.}
  \item{silent}{Logical: suppress \code{siena07} output.  Defaults to
      the value of \code{batch}.}
  \item{verbose}{Logical: print progress messages.}
  \item{memoryScale, batchUnitBudget, dynamicMinistepFactor}{Memory
      management parameters for batch sizing.}
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

  \strong{Aggregation}: The \code{level} parameter controls how the
  raw per-ego-per-choice predictions are aggregated.  For example,
  \code{level = "period"} averages across egos within each period.
  The \code{condition} parameter allows further grouping
  (e.g.\ by a density (creation vs. dissolution) or reciprocity effect).

  Multi-group data (\code{sienaGroup}) is not yet supported.
}
\value{
  A data.frame (or data.table if available) with columns depending on
  \code{level} and \code{type}.

  When \code{uncertainty = FALSE}, the \code{level = "none"} output
  contains one row per ego-choice-period (static) or
  ego-choice-chain-ministep (dynamic), with columns including:
  \describe{
    \item{\code{period}, \code{ego}, \code{choice}}{Coordinates.}
    \item{\code{changeUtil}}{The multinomial logit utility for this
          choice.}
    \item{\code{changeProb}}{The multinomial choice probability.}
    \item{\code{tieProb}}{The tie probability (if
          \code{type = "tieProb"}).}
    \item{Effect contribution columns}{One column per effect with
          the change statistic value.}
  }

  When \code{uncertainty = TRUE}, the output is at the requested
  \code{level} of aggregation with columns \code{Mean}, \code{SE},
  quantiles, etc.
}
\author{Daniel Gotthardt}
\seealso{
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

## Static changeProb, no uncertainty, unaggregated
p_raw <- predict(ans, newdata=mydata, type="changeProb",
                 uncertainty=FALSE, level="none")
head(p_raw)

## Static tieProb, aggregated by period (with uncertainty)
p_tie <- predict(ans, newdata=mydata, type="tieProb",
                 level="period", nsim=100, verbose=FALSE)
head(p_tie)
}
\keyword{models}
