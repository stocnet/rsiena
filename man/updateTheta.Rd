\name{updateTheta}
\alias{updateTheta}
\alias{updateSpecification}
\title{
A function to update the initial values of theta, and a function
to update an effects object.
}
\description{
  \code{updateTheta} copies the final values of any matching selected effects   
  from a \code{\link{sienaFit}} object to a Siena effects object.\cr   
  \code{updateSpecification} includes in a Siena effects object   
  a set of effects that are included in another effects object.  
}
\usage{
updateTheta(effects, prevAns, varName=NULL)
updateSpecification(effects.to, effects.from, name.to=NULL, name.from=NULL)
}
\arguments{
  \item{effects}{Object of class \code{\link{sienaEffects}}.}
  \item{prevAns}{Object of class \code{\link{sienaFit}}
    as returned by \code{\link{siena07}}.}
  \item{varName}{Character string or vector of character strings;
       is this is not \code{NULL}, the update will only applied to this
       dependent variable / these dependent variables.}
  \item{effects.to}{Object of class \code{\link{sienaEffects}}.}
  \item{effects.from}{Object of class \code{\link{sienaEffects}}.}
  \item{name.to}{Character string, name of dependent variable in \code{object.to}.}
  \item{name.from}{Character string, name of dependent variable in
       \code{object.from}.}
}
\details{
  The initial values of any selected effects
  in the input effects object which match an effect estimated in
  \code{prevAns} will be updated by \code{updateTheta}.
  If the previous run was conditional, the estimated rate parameters for
  the dependent variable on which the run was conditioned are added to
  the final value of theta. If \code{varName} is not \code{NULL},
  this update is restricted to effects for the dependent variable/s
  specified by \code{varName}.
  \cr
  By \code{updateSpecification}, the effects included in \code{effects.from}
  are also included in \code{effects.to}; if \code{name.to} and/or
  \code{name.from} is specified, this is restricted to effects
  for those dependent variables.
  If \code{name.to = "all"} (should then not be used as variable name!),
  the effects for all dependent variables will be updated. \cr
  Correspondence between effects is defined by \code{"name", "shortName"
  "type", "groupName", "interaction1", "interaction2", "period", "effect1",
  "effect2",} and \code{"effect3"}. This means that inclusion of user-defined
  interactions will be updated only if they were available (i.e., defined)
  already in \code{effects.to}.
}
\value{
  Updated effects object.
}
\note{Using \code{updateTheta} explicitly before calling \code{\link{siena07}}
  rather than using it via the argument
  \code{prevAns} of \code{\link{siena07}} will not permit the use of
  the previous derivative matrix.
  In most cases, using \code{\link{siena07}} with \code{prevAns} will be more
  efficient.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley, Tom A.B. Snijders}

\seealso{\code{\link{siena07}}, \code{\link{getEffects}}
}
\examples{
## For updateTheta:
mynet1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff1 <- getEffects(mydata)
myeff1 <- includeEffects(myeff1, transTrip)
myalgorithm <- sienaAlgorithmCreate(nsub=1, n3=100, projname=NULL)
ans <- siena07(myalgorithm, data=mydata, effects=myeff1, batch=TRUE)
ans$theta
(myeff <- updateTheta(myeff1, ans))
##
## For updateSpecification:
myeff2 <- getEffects(mydata)
myeff2 <- includeEffects(myeff2, inPop)
updateSpecification(myeff2, myeff1)
# Create (meaningless) two-dimensional dependent network
mynet1 <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mynet2 <- sienaDependent(array(c(s503, s501), dim=c(50, 50, 2)))
mydata12 <- sienaDataCreate(mynet1, mynet2)
myeff12 <- getEffects(mydata12)
myeff.new <- getEffects(mydata12)
(myeff12 <- includeEffects(myeff12, inPop, outPop, outAct))
# update myeff.new only for mynet1:
updateSpecification(myeff.new, myeff12)
# update myeff.new for all dependent networks:
(myeff.updated <- updateSpecification(myeff.new, myeff12, "all"))
# use multivariate effects object to update univariate effects object:
myeff1 <- getEffects(sienaDataCreate(mynet1))
updateSpecification(myeff1, myeff.updated)
}
\keyword{ models }
