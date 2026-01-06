\name{update_theta}
\alias{update_theta}
\alias{updateTheta}
\alias{update_theta.sienaEffects}
\alias{update_specification}
\alias{updateSpecification}
\alias{update_specification.sienaEffects}
\title{
A function to update the initial values of theta, and a function
to update an effects object.
}
\description{
  \code{update_theta} copies the final values of any matching selected effects
  from a \code{\link{sienaFit}} object to a Siena effects object.\cr
  \code{update_specification} includes in a Siena effects object
  a set of effects that are included in another effects object.\cr
  The functions \code{update_theta} and
  \code{updateTheta} are identical;\cr
  the same holds for
  \code{update_specification} and \code{updateSpecification}.
  The first of these names are preferred,
  the others are kept for backward compatibility.
}
\usage{
\method{update_theta}{sienaEffects}(x, prevAns, varName=NULL, ...)
\method{update_specification}{sienaEffects}(x, effects.from,
         effects.extra=NULL, name.to=NULL, name.from=NULL, ...)

updateTheta(effects, prevAns, varName=NULL)
updateSpecification(effects.to, effects.from,
         effects.extra=NULL, name.to=NULL, name.from=NULL)
}
\arguments{
  \item{x}{Object of class \code{\link{sienaEffects}}.}
  \item{prevAns}{Object of class \code{\link{sienaFit}}
    as returned by \code{\link{siena07}}.}
  \item{varName}{Character string or vector of character strings;
       is this is not \code{NULL}, the update will only applied to this
       dependent variable / these dependent variables.}
  \item{effects.to}{Object of class \code{\link{sienaEffects}}.}
  \item{effects.from}{Object of class \code{\link{sienaEffects}}.}
  \item{effects.extra}{Object of class \code{\link{sienaEffects}}.}
  \item{name.to}{Character string, name of dependent variable in \code{object.to}.}
  \item{name.from}{Character string, name of dependent variable in
       \code{object.from}.}
  \item{effects}{Object of class \code{\link{sienaEffects}}.}
  \item{...}{Additional arguments, not used now.}
}
\details{
  The initial values of any selected effects
  in the input effects object which match an effect estimated in
  \code{prevAns} will be updated by \code{update_theta}.
  If the previous run was conditional, the estimated rate parameters for
  the dependent variable on which the run was conditioned are added to
  the final value of theta. If \code{varName} is not \code{NULL},
  this update is restricted to effects for the dependent variable/s
  specified by \code{varName}.
  \cr
  By \code{update_specification}, the effects included in \code{effects.from}
  are also included in \code{x} (for \code{update_specification}) \code{effects.to}
  (for \code{updateSpecification}); if \code{name.to} and/or
  \code{name.from} is specified, this is restricted to effects
  for those dependent variables.\cr
  If \code{effects.from} contains interaction effects, the corresponding
  main effects will be looked for in \code{effects.from}; if they are not
  found there, they will be looked for in \code{effects.extra}.
  It is not guaranteed that this will be successful.
  For \code{effects.extra}, it is best to use an effects object constructed for
  the same data set as \code{effects.from}, and by the same version of RSiena.
}
\value{
  Updated effects object.
}
\note{Using \code{update_theta} explicitly before calling \code{\link{siena}}
  rather than using it via the argument
  \code{prevAns} of \code{\link{siena}} will not permit the use of
  the previous derivative matrix for option \code{Dolby}.
  In most cases, using \code{\link{siena}} with \code{prevAns}
  will therefore be more efficient;
  but sometimes the use of  \code{update_theta} gives better results.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley, Tom A.B. Snijders}
\seealso{
    \code{\link{siena}}, \code{\link{make_specification}}
}
\examples{
## For updateTheta:
mynet1 <- as_dependent_rsiena(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- make_data_rsiena(mynet1)
myeff1 <- make_specification(mydata)
myeff1 <- set_effect(myeff1, transTrip)
myeff1 <- set_interaction(myeff1, list(recip, inPop))
myalgo <- set_algorithm_saom(nsub=1, n3=100, seed=1291)
(ans <- siena(mydata, effects=myeff1, control_algo=myalgo,
                                          batch=TRUE))
ans$theta
(myeff <- update_theta(myeff1, ans))
##
## For updateSpecification:
myeff2 <- make_specification(mydata)
myeff2 <- set_effect(myeff2, outAct)
update_specification(myeff2, myeff1)
}
\keyword{models}
