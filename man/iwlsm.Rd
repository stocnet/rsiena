%file RSiena/man/iwlsm.Rd
% based on file MASS/man/rlm.Rd
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{iwlsm}
\alias{iwlsm}
\alias{iwlsm.default}
\alias{iwlsm.formula}
\alias{print.iwlsm}
\alias{predict.iwlsm}
\alias{psi.iwlsm}
\title{
Function to fit an iterated weighted least squares model.
}
\description{
  Fits an iterated weighted least squares model.
}
\usage{
iwlsm(x, ...)

\method{iwlsm}{formula}(formula, data, weights, ses, \dots, subset, na.action,
    method = c("M", "MM", "model.frame"),
    wt.method = c("inv.var", "case"),
    model = TRUE, x.ret = TRUE, y.ret = FALSE, contrasts = NULL)

\method{iwlsm}{default}(x, y, weights, ses, \dots, w = rep(1/nrow(x), nrow(x)),
    init = "ls", psi = psi.iwlsm,
    scale.est = c("MAD", "Huber", "proposal 2"), k2 = 1.345,
    method = c("M", "MM"), wt.method = c("inv.var", "case"),
    maxit = 20, acc = 1e-4, test.vec = "resid", lqs.control = NULL)

psi.iwlsm(u, k, deriv = 0, w, sj2, hh)
}
\arguments{
  \item{formula}{
    a formula of the form \code{y ~ x1 + x2 + \dots}.
  }
  \item{data}{
    data frame from which variables specified in \code{formula} are
    preferentially to be taken.
  }
  \item{weights}{
    a vector of prior weights for each case.
  }
  \item{subset}{
    An index vector specifying the cases to be used in fitting.
  }
  \item{ses}{Estimated variance of the responses. Will be paseed to
	\code{psi} as sj2
  }
  \item{na.action}{
    A function to specify the action to be taken if \code{NA}s are found.
    The \sQuote{factory-fresh} default action in \R is
    \code{\link{na.omit}}, and can be changed by
    \code{\link{options}(na.action=)}.
  }
  \item{x}{
    a matrix or data frame containing the explanatory variables.
  }
  \item{y}{
    the response: a vector of length the number of rows of \code{x}.
  }
  \item{method}{Must be "M". (argument not used here).
  }
  \item{wt.method}{
    are the weights case weights (giving the relative importance of case,
    so a weight of 2 means there are two of these) or the inverse of the
    variances, so a weight of two means this error is half as variable?
	This will not work at present.
  }
  \item{model}{
    should the model frame be returned in the object?
  }
  \item{x.ret}{
    should the model matrix be returned in the object?
  }
  \item{y.ret}{
    should the response be returned in the object?
  }
  \item{contrasts}{
    optional contrast specifications: se \code{\link{lm}}.
  }
  \item{w}{
    (optional) initial down-weighting for each case. Will not work at present.
  }
  \item{init}{
    (optional) initial values for the coefficients OR a method to find
    initial values OR the result of a fit with a \code{coef} component.  Known
    methods are \code{"ls"} (the default) for an initial least-squares fit
    using weights \code{w*weights}, and \code{"lts"} for an unweighted
    least-trimmed squares fit with 200 samples. Probably not functioning.
  }
  \item{psi}{
    the psi function is specified by this argument.  It must give
    (possibly by name) a function \code{g(x, \dots, deriv, w)} that for
    \code{deriv=0} returns psi(x)/x and for \code{deriv=1} returns
    some value.  Extra arguments may be passed in via \code{\dots}.
  }
  \item{scale.est}{
    method of scale estimation: re-scaled MAD of the residuals (default)
    or Huber"s proposal 2 (which can be selected by either \code{"Huber"}
    or \code{"proposal 2"}).
  }
  \item{k2}{
    tuning constant used for Huber proposal 2 scale estimation.
  }
  \item{maxit}{
    the limit on the number of IWLS iterations.
  }
  \item{acc}{
    the accuracy for the stopping criterion.
  }
  \item{test.vec}{
    the stopping criterion is based on changes in this vector.
  }
  \item{\dots}{
    additional arguments to be passed to \code{iwlsm.default} or to the \code{psi}
    function.
  }
  \item{lqs.control}{
    An optional list of control values for \code{\link[MASS]{lqs}}.
  }
  \item{u}{
    numeric vector of evaluation points.
  }
  \item{k}{    tuning constant. Not used.
  }
  \item{deriv}{
    \code{0} or \code{1}: compute values of the psi function or of its
    first derivative. (Latter not used).
  }
  \item{sj2}{Estimated variance of the responses}
  \item{hh}{Diagonal values of the hat matrix}
}
\details{This function is very slightly adapted from \code{rlm} in
  packages \code{MASS}. It alternates between weighted
  least squares and estimation of variance on the basis of a common
  variance. The function \code{psi.iwlsm} calculates the weights
  for the next iteration. Used by \code{siena08} to combine estimates
  from different \code{sienaFits}.
}
\value{
  An object of class \code{"iwlsm"} inheriting from \code{"lm"}.
  Note that the \code{df.residual} component is deliberately set to
  \code{NA} to avoid inappropriate estimation of the residual scale from
  the residual mean square by \code{"lm"} methods.

  The additional components not in an \code{lm} object are

  \item{s}{
    the robust scale estimate used
  }
  \item{w}{
    the weights used in the IWLS process
  }
  \item{psi}{
    the psi function with parameters substituted
  }
  \item{conv}{
    the convergence criteria at each iteration
  }
  \item{converged}{
    did the IWLS converge?
  }
  \item{wresid}{
    a working residual, weighted for \code{"inv.var"} weights only.
  }
}
\references{
  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.}  Fourth edition.  Springer.
  See also \url{https://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Ruth Ripley}
\note{The function has been changed as little as possible, but
   has only been used with default arguments. The other options have
    been retained just in
    case they may prove useful.
}
\seealso{
  \code{\link{siena08}}, \code{\link{sienaMeta}}, \code{\link{sienaFit}}
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
}
}
\keyword{models }
