\name{test_parameter}
\alias{test_parameter}
\alias{test_parameter.sienaFit}
\alias{Wald}
\alias{Wald.RSiena}
\alias{Multipar.RSiena}
\alias{score.Test}
\alias{scoreTest}
\alias{sienaTest}
\alias{siena.Test}
\alias{scoretest}
\alias{sienatest}
\alias{testSame.RSiena}
\title{  Wald and score tests for RSiena results }
\description{
  These functions test parameters in RSiena results
  estimated by \code{\link{siena}}.
  Tests can be Wald-type (if the parameters were estimated)
  or score-type tests (if the parameters were fixed and tested).\cr
  The function \code{test_parameter} chooses among the
  functions \code{Wald.RSiena},  \code{Multipar.RSiena},
  \code{testSame.RSiena},   and \code{score.Test},
  based on the input given. The export of the latter three functions is superfluous,
  but they are kept for backward compatibility.
}
\usage{

\method{test_parameter}{sienaFit}(x, method=NULL, tested=NULL, tested2=NULL, ...)

Wald.RSiena(A, x)

Multipar.RSiena(x, tested)

testSame.RSiena(x, e1, e2)

score.Test(x, tested=x$test)

}
\arguments{
  \item{x}{An object of class \code{\link{sienaFit}},
    resulting from a call to \code{\link{siena}.}}
  \item{method}{\code{NULL}, \code{"all"}, \code{"same"}, or \code{"score"}.}
  \item{tested}{One number or a vector of numbers between 1 and \code{p},
      where \code{p = x$pp}, the number of parameters in \code{x}
      excluding the basic rate parameters used for
      conditional estimation; these indicate the tested parameters.
      Can also be a logical vector of length \code{p}; \cr
      or, for \code{test_parameter}, a \code{k} * \code{p} matrix. }
  \item{tested2}{A single number or vector as mentioned for \code{tested};
     these two parameters are tested for equality for \code{method = "same"}.}
  \item{A}{A \code{k} * \code{p} matrix, where \code{p = ans$pp},
  the number of parameters in \code{ans} excluding the basic rate parameters
  used for conditional estimation.}
  \item{e1,e2}{Each an integer number between 1 and \code{p}, or a vector
    of such numbers; the hypothesis tested is that the parameters for effects
    with number/s \code{e1} are equal to those in \code{e2}.}
 \item{...}{Additional arguments (currently not used.)}
}

\details{ The tests for \code{test_parameter}, \code{Wald.RSiena},
    \code{Multipar.RSiena}, and \cr
    \code{testSame.RSiena} are Wald-type tests
    (but see below for \code{test_parameter} with \code{method=score}.\cr
    For \code{test_parameter} and \code{Wald.RSiena},
    if \code{tested} is a matrix \code{A}, the hypothesis tested
   is \eqn{A\theta = 0}{A theta = 0}, where \eqn{\theta}{theta} is
   the parameter estimated in the process leading to \code{x}.\cr
    For \code{test_parameter} and \code{Multipar.RSiena}
    if \code{tested} is a number or vector, a multivariate
    (i.e., simultaneous) test  is given of the hypothesis that
    the  parameters mentioned in \code{tested} are 0.   \cr
    For \code{test_parameter} if \code{"all"},
    a test is given of all parameters in the model.\cr
    For \code{test_parameter}
    if \code{"same"}, the items \code{"tested"} and  \code{"tested2"}
    should be single numbers or vectors of equal length, and the test
    of equality of the corresponding two parameters is given.\cr
    For \code{test_parameter} if \code{"score"}, and for \code{score.Test},
    the parameters with numbers in \code{"tested"}
    should all be used in \code{x} as fixed-and-tested parameters,
    and the multivariate score test of these parameters is given
    of the hypothesis that the tested parameters have
    the value indicated in the effects object used for obtaining \code{x}.  \cr
    The numbering of parameters is as in \code{print(x)};
    if conditional estimation was used, numbered
    as the 'Other parameters').

   These tests should be carried out only when convergence is adequate
   (overall maximum convergence ratio less than 0.25
   and all \emph{t}-ratios for convergence less than 0.1 in absolute value).

   These functions have their own print method, see \code{\link{print.sienaTest}}.
}
\value{An object of class \code{sienaTest}, which is a list with elements:
  \item{chisquare: }{ The test statistic, assumed to have a chi-squared null
    distribution.}
  \item{df: }{ The degrees of freedom.}
  \item{pvalue: }{ The associated \emph{p}-value.}
  \item{onesided: }{For \code{df}=1, the onesided test statistic.}
  \item{efnames: }{For \code{Multipar.RSiena} and \code{score.Test}, the names
              of the tested effects.}
}
\references{See the manual and \url{https://www.stats.ox.ac.uk/~snijders/siena/}

M. Schweinberger (2012). Statistical modeling of network panel data:
Goodness-of-fit. \emph{British Journal of Statistical
and Mathematical Psychology} \bold{65}, 263--281.}

\author{Tom Snijders}

\seealso{ \code{\link{siena}}, \code{\link{print.sienaTest}} }

\examples{
myalgo <- set_algorithm_saom(nsub=1, n3=40, seed=1777)
# nsub=1 and n3=40 is used here for having a brief computation,
# not for practice.
mynet <- as_dependent_rsiena(array(c(s501, s502), dim=c(50, 50, 2)))
mydata <- make_data_rsiena(mynet)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, list(transTrip, transTies))
myeff <- set_effect(myeff, list(outAct, outPop), fix=TRUE, test=TRUE)
(ans <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE))
A <- matrix(0, 2, 6)
A[1, 3] <- 1
A[2, 4] <- 1
(wa <- test_parameter(ans, tested=A))
wa
# An alternative specification of the above is:
test_parameter(ans, tested=c(3, 4))
# The following two are also equivalent:
sct <- test_parameter(ans, method="score",
                tested=c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE))
sct <- test_parameter(ans, method="score", tested=6)
print(sct)
# Getting all 1-df score tests separately:
for (i in which(ans$test)){
   sct <- test_parameter(ans, method="score", tested=i)
   print(sct)}
# Testing that endowment and creation effects are identical:
myeff1 <- make_specification(mydata)
myeff1 <- set_effect(myeff1, list(transTrip, transTies))

myeff1 <- make_specification(mydata)
myeff1 <- set_effect(myeff1, recip, include=FALSE)
myeff1 <- set_effect(myeff1, recip, type='creation')
(myeff1 <- set_effect(myeff1, recip, type='endow'))
ans1 <- siena(mydata, effects=myeff1, control_algo=myalgo, batch=TRUE)
test_parameter(ans1, method="same", tested=2, tested2=3)
}
\keyword{tests}
