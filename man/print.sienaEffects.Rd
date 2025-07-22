\name{print.sienaEffects}
\alias{print.sienaEffects}
\alias{summary.sienaEffects}
\alias{print.summary.sienaEffects}
\title{
  Print methods for Siena effects objects
}
\description{
  Prints the major columns of the effects object.   
  Or all, with any non-atomic columns listed separately.
}
\usage{
\method{print}{sienaEffects}(x, fileName = NULL, includeOnly=TRUE,
   expandDummies = FALSE, includeRandoms = FALSE, dropRates=FALSE, 
   includeShortNames=FALSE, ...)
\method{summary}{sienaEffects}(object, fileName = NULL,
   includeOnly=TRUE, expandDummies = FALSE, ...)
\method{print}{summary.sienaEffects}(x, fileName = NULL, ...)
}
\arguments{
  \item{object}{An object of class \code{sienaEffects}.}
  \item{x}{An object of class \code{sienaEffects} or
    \code{summary.sienaEffects} as appropriate.
  }
  \item{fileName}{Character string denoting file name if file output desired.}
  \item{includeOnly}{Boolean. If \code{TRUE}, only effects with the include
    flag \code{TRUE} will be printed.}
  \item{expandDummies}{Interpret the \code{timeDummy} column and show
    any effects which would be added by \code{sienaTimeFix}.}
  \item{includeRandoms}{Boolean. If \code{TRUE}, also the
    \code{randomEffects} column will be printed.}
  \item{includeShortNames}{Boolean. If \code{TRUE}, also the
    \code{shortName} column will be printed.}
  \item{dropRates}{Boolean. If \code{TRUE}, do not print the rows for basic rate
    effects.}
  \item{\dots}{For extra arguments (none used at present).}
}
\value{
  The function \code{print.sienaEffects} prints details of the main
  columns of the selected rows of the effects object.\cr
  If the effects object includes statistics for the
  Generalized Method of Moments (GMoM), as included by function
 \code{\link{includeGMoMStatistics}} and for which \code{type=gmm},
  the print consists of two parts: the first consists of the
  included effects for the probability model, the second
  of the statistics used for GMoM estimation.

 The function \code{summary.sienaEffects} checks the rows for valid
 printing via \code{print.data.frame} and excludes any that will
 fail. The OK columns are printed first, followed by any others.

 Output from either can be directed to a file by using the argument
 \code{filename}.
}
\references{
  See \url{https://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Ruth Ripley, modifications by Tom Snijders and Viviana Amati.
}
\seealso{\code{\link{sienaEffects}},
    \code{\link{getEffects}}, \code{\link{includeEffects}},
    \code{\link{includeGMoMStatistics}},
    \code{\link{sienaTimeTest}},
    \code{\link{effectsDocumentation}}
}
\examples{
mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- sienaDependent(s50a, type="behavior")
mycovar <- coCovar(rnorm(50))
mydyadcovar <- coDyadCovar(matrix(as.numeric(rnorm(2500) > 2), nrow=50))
mydata <- sienaDataCreate(mynet1, mybeh, mycovar, mydyadcovar)
myeff <- getEffects(mydata)
myeff
print(myeff, includeShortNames=TRUE)
summary(myeff)
}
\keyword{print}
