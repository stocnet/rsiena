\name{entropy}
\alias{entropy}
\alias{entropy.sienaFit}
\alias{entropy.default}
\alias{sienaEntropy}
\title{Entropy-based degree of certainty for SAOMs}
\description{
  \code{entropy} is a generic function.  The \code{sienaFit} method
  computes the normalised entropy (degree of certainty) \eqn{R_H}
  with optional uncertainty quantification via \code{sienaPostestimate}.

  \code{sienaEntropy} is a standalone function that computes \eqn{R_H}
  from predicted choice probabilities (the output of
  \code{\link{predict.sienaFit}}) or from raw probability vectors.

  The degree of certainty is defined as
  \eqn{R_H = 1 - H / \log(J)} where \eqn{H = -\sum_j p_j \log p_j}
  is the Shannon entropy and \eqn{J} is the number of available choices.
  See Snijders (2004), equation (6).
}
\usage{
\method{entropy}{sienaFit}(object, data, effects = NULL, depvar = NULL,
    dynamic = FALSE, algorithm = NULL, n3 = 500,
    useChangeContributions = FALSE,
    level = "ego", sum_fun = mean, na.rm = TRUE,
    uncertainty = TRUE, uncertaintyMode = c("batch", "stream"),
    nsim = 1000, uncertaintySd = TRUE, uncertaintyCi = TRUE,
    uncertaintyProbs = c(0.025, 0.5, 0.975),
    useCluster = FALSE, nbrNodes = 1,
    clusterType = c("PSOCK", "FORK"),
    batchDir = "temp", batchSize = NULL,
    keepBatch = FALSE, verbose = TRUE,
    memoryScale = NULL, batchUnitBudget = 2.5e8,
    dynamicMinistepFactor = 10, ...)

\method{entropy}{default}(object, ...)

sienaEntropy(x, group_id = NULL, level = "ego",
             sum_fun = mean, na.rm = TRUE)
}
\arguments{
  \item{object}{For \code{entropy.sienaFit}: a \code{sienaFit} object.
      For \code{entropy.default}: a numeric distribution vector.}
  \item{x}{For \code{sienaEntropy}: a data.frame with a \code{changeProb}
      column (typically from \code{\link{predict.sienaFit}}), or a
      numeric vector of choice probabilities.}
  \item{data}{\code{siena} data object.}
  \item{effects}{Optional \code{sienaEffects} object.}
  \item{depvar}{Name of dependent variable. Default: first depvar in data.}
  \item{dynamic}{Logical: if \code{TRUE}, use dynamic contributions.}
  \item{algorithm}{Required when \code{dynamic = TRUE}.}
  \item{n3}{Number of chains for dynamic simulation.}
  \item{useChangeContributions}{Use pre-stored change contributions.}
  \item{level}{Aggregation level (default \code{"ego"}).}
  \item{sum_fun}{Aggregation function (default \code{mean}).}
  \item{na.rm}{Logical: remove NAs.}
  \item{uncertainty}{Logical: compute uncertainty via theta-draws.}
  \item{uncertaintyMode}{\code{"batch"} or \code{"stream"}.}
  \item{nsim}{Number of theta-draws for uncertainty.}
  \item{uncertaintySd, uncertaintyCi, uncertaintyProbs}{Uncertainty
      summary controls.}
  \item{useCluster, nbrNodes, clusterType}{Parallelisation settings.}
  \item{batchDir, batchSize, keepBatch}{Batch mode settings.}
  \item{verbose}{Logical: print progress messages.}
  \item{memoryScale, batchUnitBudget, dynamicMinistepFactor}{Memory
      management.}
  \item{group_id}{Integer vector of group IDs (required when \code{x}
      is a numeric vector for \code{sienaEntropy}).}
  \item{\dots}{Additional arguments.}
}
\value{
  \code{entropy.sienaFit}: a data.frame with \code{R_entropy} and, when
  \code{uncertainty = TRUE}, additional uncertainty columns (SE,
  quantiles).

  \code{entropy.default}: a single numeric value (degree of certainty).

  \code{sienaEntropy}: a data.frame with columns \code{group} and
  \code{R_entropy}.
}
\references{
  Snijders, T.A.B. (2004),
  Explained Variation in Dynamic Network Models.
  \emph{Mathematics and Social Sciences}, \bold{168}, 31--41.
}
\author{Daniel Gotthardt, Tom A.B. Snijders}
\seealso{\code{\link{interpret_size}}, \code{\link{predict.sienaFit}}}
\examples{
## entropy.default: degree of certainty from a probability vector
p_uniform <- rep(1/4, 4)
entropy(p_uniform)  # 0 (maximum uncertainty)

p_skewed <- c(0.9, 0.05, 0.05)
entropy(p_skewed)   # high certainty

\dontrun{
## entropy.sienaFit: degree of certainty from a fitted model
myalgo <- set_algorithm_saom(nsub=1, n3=50, seed=1293)
mynet1 <- as_dependent_rsiena(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- make_data_rsiena(mynet1)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, list(density, recip, transTrip))
ans <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE)

ent <- entropy(ans, data=mydata, depvar="mynet1",
               nsim=20, verbose=FALSE)
ent
}
}
\keyword{models}
