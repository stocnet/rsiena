% file RSiena/man/summary.iwlsm.Rd
% based on file MASS/man/summary.rlm.Rd
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{summary.iwlsm}
\alias{summary.iwlsm}
\alias{print.summary.iwlsm}
\title{
Summary method for Iterative Weighted Least Squares Models
}
\description{
\code{summary} method for objects of class \code{"iwlsm"}
}
\usage{
\method{summary}{iwlsm}(object, method = c("XtX", "XtWX"),
 correlation = FALSE, \dots)
}
\arguments{
\item{object}{
the fitted model.
This is assumed to be the result of some fit that produces
an object inheriting from the class \code{iwlsm}, in the sense that
the components returned by the \code{iwlsm} function will be available.
}
\item{method}{
Should the weighted (by the IWLS weights) or unweighted cross-products
matrix be used?
}
\item{correlation}{
logical. Should correlations be computed (and printed)?
}
\item{\dots}{
arguments passed to or from other methods.
}}
\value{
If printing takes place, only a null value is returned.
Otherwise, a list is returned with the following components.
Printing always takes place if this function is invoked automatically
as a method for the \code{summary} function.

\item{correlation}{
The computed correlation coefficient matrix for the coefficients in the model.
}
\item{cov.unscaled}{
The unscaled covariance matrix; i.e, a matrix such that multiplying it by
an estimate of the error variance produces an estimated covariance matrix
for the coefficients.
}
\item{sigma}{
The scale estimate.
}
\item{stddev}{
A scale estimate used for the standard errors.
}
\item{df}{
The number of degrees of freedom for the model and for residuals.
}
\item{coefficients}{
A matrix with three columns, containing the coefficients, their standard errors
and the corresponding t statistic.
}
\item{terms}{
The terms object used in fitting this model.
}}
\details{
This function is a method for the generic function
\code{summary()} for class \code{"iwlsm"}.
It can be invoked by calling \code{summary(x)} for an
object \code{x} of the appropriate class, or directly by
calling \code{summary.iwlsm(x)} regardless of the
class of the object.
}
\references{
  Venables, W. N. and Ripley, B. D. (2002),
  \emph{Modern Applied Statistics with S.} Fourth edition.  Springer.
  See also  \url{https://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Adapted by Ruth Ripley}
\seealso{
\code{\link{summary}}
}
\examples{
\dontrun{
##not enough data here for a sensible example, but shows the idea.
myalgorithm <- sienaAlgorithmCreate(nsub=2, n3=100)
mynet1 <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
mydata1 <- sienaDataCreate(mynet1)
mydata2 <- sienaDataCreate(mynet2)
myeff1 <- getEffects(mydata1)
myeff2 <- getEffects(mydata2)
myeff1 <- setEffect(myeff1, transTrip, fix=TRUE, test=TRUE)
myeff2 <- setEffect(myeff2, transTrip, fix=TRUE, test=TRUE)
myeff1 <- setEffect(myeff1, cycle3, fix=TRUE, test=TRUE)
myeff2 <- setEffect(myeff2, cycle3, fix=TRUE, test=TRUE)
ans1 <- siena07(myalgorithm, data=mydata1, effects=myeff1, batch=TRUE)
ans2 <- siena07(myalgorithm, data=mydata2, effects=myeff2, batch=TRUE)
meta <- siena08(ans1, ans2)
metadf <- split(meta$thetadf, meta$thetadf$effects)[[1]]
metalm <- iwlsm(theta ~ tconv, metadf, ses=se^2)
summary(metalm)
}
}
\keyword{methods}
