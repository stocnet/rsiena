\name{print.sienaTest}
\alias{print.sienaTest}
\alias{sienaTest.methods}
\title{
  Print method for Wald and score tests for RSiena results
}
\description{
  This method prints Wald-type and score-type tests for results   
  estimated by \code{\link{siena07}}.
}
\usage{
\method{print}{sienaTest}(x, \dots)
}
\arguments{
  \item{x}{An object of type \code{sienaTest}, produced by
    \code{\link{Wald.RSiena}}, \code{\link{Multipar.RSiena}},
    or \code{\link{score.Test}}.}
  \item{\dots}{Extra arguments (not used at present).}
}

\details{The functions \code{\link{Wald.RSiena}}, \code{\link{Multipar.RSiena}},
    and \code{\link{score.Test}} produce an object of type \code{\link{sienaTest}}.
    These can be printed by this method.}
\value{An object of type  \code{\link{sienaTest}}.
}
\author{Tom Snijders}
\seealso{ \code{\link{siena07}}, \code{\link{Wald.RSiena}},
    \code{\link{Multipar.RSiena}}, \code{\link{score.Test}}}

\examples{
mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet)
myeff <- getEffects(mydata)
myalgorithm <- sienaAlgorithmCreate(nsub=1, n3=40, seed=123, projname=NULL)
# nsub=1 and n3=40 is used here for having a brief computation,
# not for practice.
myeff <- includeEffects(myeff, transTrip, transTies)
myeff <- includeEffects(myeff, outAct, outPop, fix=TRUE, test=TRUE)
(ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE))
mprs <- Multipar.RSiena(ans, 3, 4)
print(mprs)
}
\keyword{ print }
