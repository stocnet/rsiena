\name{Wald}
\alias{Wald.RSiena}
\alias{Multipar.RSiena}
\alias{score.Test}
\alias{scoreTest}
\alias{sienaTest}
\alias{testSame.RSiena}
\title{  Wald and score tests for RSiena results }
\description{
  These functions test parameters in RSiena results
  estimated by \code{\link{siena07}}.
  Tests can be Wald-type (if the parameters were estimated)
  or score-type tests (if the parameters were fixed and tested).
}
\usage{
Wald.RSiena(A, ans)

Multipar.RSiena(ans, \dots)

score.Test(ans, test=ans$test)

testSame.RSiena(ans, e1, e2)
}
\arguments{
  \item{A}{A \code{k} * \code{p} matrix, where \code{p = ans$pp},
  the number of parameters in \code{ans} excluding the basic rate parameters
  used for conditional estimation.}
  \item{ans}{An object of class \code{\link{sienaFit}},
    resulting from a call to \code{\link{siena07}.}}
  \item{\dots}{One or more integer numbers between 1 and \code{p}, specifying
    the tested effects
    (numbered as in \code{print(ans)}; if conditional estimation was used, numbered
    as the 'Other parameters').}
  \item{test}{One or more integer numbers between 1 and \code{p},
    or a logical vector of length \code{p}; these should specify the tested
    effects (numbered as described for the \dots).}
  \item{e1,e2}{Each an integer number between 1 and \code{p}, or a vector
    of such numbers; the hypothesis tested is that the parameters for effects
    with number/s \code{e1} are equal to those in \code{e2}.}
}

\details{\code{Wald.RSiena} produces a Wald-type test,
   applicable to estimated parameters. \code{Multipar.RSiena} and
   \code{testSame.RSiena} are special cases of \code{Wald.RSiena}.
   The hypothesis tested by \code{Wald.RSiena}
   is \eqn{A\theta = 0}{A theta = 0}, where \eqn{\theta}{theta} is
   the parameter estimated in the process leading to \code{ans}.

   The hypothesis tested by \code{Multipar.RSiena} is that all
   parameters given in \eqn{\ldots}{...} are 0.

   The hypothesis tested by \code{testSame.RSiena} is that all
   parameters given in \code{e1} are equal to those in \code{e2}.

   \code{score.Test} produces a score-type test.
   The tested effects for \code{score.Test} should have been specified
   in \code{\link{includeEffects}} or \code{\link{setEffect}} with
   \code{fix=TRUE, test=TRUE}, i.e., they should not have been estimated.
   The hypothesis tested by \code{score.Test} is that the tested parameters have
   the value indicated in the effects object used for obtaining \code{ans}.

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

\seealso{ \code{\link{siena07}}, \code{\link{print.sienaTest}} }

\examples{
mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mydata <- sienaDataCreate(mynet)
myeff <- getEffects(mydata)
myalgorithm <- sienaAlgorithmCreate(nsub=1, n3=40, seed=1777, projname=NULL)
# nsub=1 and n3=40 is used here for having a brief computation,
# not for practice.
myeff <- includeEffects(myeff, transTrip, transTies)
myeff <- includeEffects(myeff, outAct, outPop, fix=TRUE, test=TRUE)
(ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE))
A <- matrix(0, 2, 6)
A[1, 3] <- 1
A[2, 4] <- 1
wa <- Wald.RSiena(A, ans)
wa
# A shortcut for the above is:
Multipar.RSiena(ans, 3, 4)
# The following two are equivalent:
sct <- score.Test(ans, c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE))
sct <- score.Test(ans,6)
print(sct)
# Getting all 1-df score tests separately:
for (i in which(ans$test)){
   sct <- score.Test(ans,i)
   print(sct)}
# Testing that endowment and creation effects are identical:
myeff1 <- getEffects(mydata)
myeff1 <- includeEffects(myeff1, recip, include=FALSE)
myeff1 <- includeEffects(myeff1, recip, type='creation')
(myeff1 <- includeEffects(myeff1, recip, type='endow'))
(ans1 <- siena07(myalgorithm, data=mydata, effects=myeff1, batch=TRUE))
testSame.RSiena(ans1, 2, 3)
}
\keyword{ tests }
