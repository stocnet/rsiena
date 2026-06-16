\name{print.sienaTest}
\alias{print.sienaTest}
\alias{sienaTest.methods}
\title{
  Print method for Wald and score tests for RSiena results
}
\description{
  This method prints Wald-type and score-type tests for results
  estimated by \code{\link{siena}}.
}
\usage{
\method{print}{sienaTest}(x, \dots)
}
\arguments{
  \item{x}{An object of type \code{sienaTest}, produced by
    \code{\link{test_parameter}}.}
  \item{\dots}{Extra arguments (not used at present).}
}

\details{The function \code{\link{test_parameter}} produces
    an object of type \code{\link{sienaTest}}.
    These can be printed by this method.}
\value{An object of type  \code{\link{sienaTest}}.
}
\author{Tom Snijders}
\seealso{ \code{\link{siena}},  \code{\link{test_parameter}}.}

\examples{
mynet <- as_dependent_rsiena(array(c(s501, s502), dim=c(50, 50, 2)))
mydata <- make_data_rsiena(mynet)
myeff <- make_specification(mydata)
myalgo <- set_algorithm_saom(nsub=2, n3=100, seed=1291)
# nsub=1, n3=40 is used here for having a brief computation, not for practice.
myeff <- set_effect(myeff, list(transTrip, gwespFF))
myeff <- set_effect(myeff, list(outAct, inPop, inAct))
myeff <- set_effect(myeff, list(inInAss, outOutAss), fix=TRUE, test=TRUE)
(ans <- siena(mydata, effects=myeff, control_algo=myalgo,
                                            silent=TRUE, batch=TRUE))
test_parameter(ans, method="score", tested=ans$test)
for (i in which(ans$test)){
   sct <- test_parameter(ans, method="score", tested=i)
   print(sct)}
test_parameter(ans, method="score", tested=TRUE, onestep=TRUE)
myeff1 <- set_effect(myeff, list(inInAss, outOutAss))
(ans1 <- siena(mydata, effects=myeff1, control_algo=myalgo,
                                   prevAns=ans, silent=TRUE, batch=TRUE))
}
\keyword{print}

