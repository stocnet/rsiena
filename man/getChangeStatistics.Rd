\name{getChangeStatistics}
\alias{getChangeStatistics}
\title{Extract change statistics from a fitted SAOM}
\description{
  Extracts change statistics (sign-adjusted contributions) from a
  fitted Stochastic Actor-oriented Model as a data.frame.

  In the SAOM, contributions represent the change in a network statistic
 when switching from \dQuote{no change} to a specific tie variable.
  For the density effect, the \dQuote{no change} option has a
  contribution of \eqn{-1}; \code{getChangeStatistics} (with
  \code{flip = TRUE}) adjusts the sign of all other effects on those
  rows so that the values represent genuine change statistics.

  Both static (at observation moments) and dynamic (over simulated
  chains of ministeps) change statistics are supported.
}
\usage{
getChangeStatistics(ans, data, effects = NULL, depvar = NULL,
    dynamic = FALSE, algorithm = NULL,
    n3 = 500, useChangeContributions = FALSE,
    flip = TRUE)
}
\arguments{
  \item{ans}{A \code{sienaFit} object.}
  \item{data}{A \code{siena} / \code{sienadata} object.}
  \item{effects}{Optional \code{sienaEffects}; defaults to
      \code{ans$requestedEffects}.}
  \item{depvar}{Character: dependent variable name. Default: first depvar.}
  \item{dynamic}{Logical: if \code{TRUE}, simulate chains to obtain
      dynamic change statistics.}
  \item{algorithm}{Required when \code{dynamic = TRUE}: a
      \code{sienaAlgorithm} object.}
  \item{n3}{Number of chains for dynamic simulation.}
  \item{useChangeContributions}{Logical: reuse stored contributions
      from \code{ans}.}
  \item{flip}{Logical: if \code{TRUE} (default), apply sign adjustment
      so values represent change statistics rather than raw
      contributions.}
}
\value{
  A data.frame (or data.table if available) with coordinate columns
  (\code{group}, \code{period}, \code{ego}, \code{choice} for static;
  adds \code{chain}, \code{ministep} for dynamic) and one column per
  effect containing the (possibly sign-adjusted) change statistics.
}
\author{Daniel Gotthardt}
\seealso{
  \code{\link{marginalEffects}}, \code{\link{interpret_size}},
  \code{\link{predict.sienaFit}}
}
\examples{
myalgo <- set_algorithm_saom(nsub=1, n3=50, seed=1293)
mynet1 <- as_dependent_rsiena(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- make_data_rsiena(mynet1)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, list(density, recip, transTrip))
ans <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE)

## Static change statistics
cs <- getChangeStatistics(ans, data=mydata, depvar="mynet1")
head(cs)

## Raw contributions (no sign adjustment)
cs_raw <- getChangeStatistics(ans, data=mydata, depvar="mynet1",
                               flip=FALSE)
head(cs_raw)

\dontrun{
## Dynamic change statistics
cs_dyn <- getChangeStatistics(ans, data=mydata, depvar="mynet1",
                               dynamic=TRUE, algorithm=myalgo, n3=50)
head(cs_dyn)
}
}
\keyword{models}
