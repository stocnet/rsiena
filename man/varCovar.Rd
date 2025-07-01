\name{varCovar}
\alias{varCovar}
\title{Function to create a changing covariate object.}
\description{
  This function creates a changing covariate object from a matrix.
}
\usage{
varCovar(val, centered=TRUE, nodeSet="Actors",  warn=TRUE,
  imputationValues=NULL)
}
\arguments{
  \item{val}{ Matrix of covariate values, one row for each actor, one
    column for each period. }
  \item{centered}{Boolean: if TRUE, then the overall mean value is subtracted.}
  \item{nodeSet}{Character string containing the name of the associated
    node set. If the entire data set contains more than one node set,
    then the node sets must be specified in all data objects.}
  \item{warn}{Logical: is a warning given if all values are \code{NA},
  or all non-missing values are the same.}
  \item{imputationValues}{Matrix of covariate values of same
   dimensions as \code{val}, to be used for imputation of
   \code{NA} values (if any) in \code{val}. Must not contain
   any \code{NA}.}
}
\details{
  When part of a Siena data object, the covariate is assumed to be
  associated with node set \code{nodeSet} of the Siena data object.
  In practice, the node set needs to be specified only in the case of
  the use of the covariate with a two-mode network.\cr
  If there are any \code{NA} values in \code{val}, and
  \code{imputationValues} is given, then the corresponding elements
  of \code{imputationValues} are used for imputation.
  If \code{imputationValues} is \code{NULL}, imputation
  is by the overall mean value. In both cases, cases with imputed values
  are not used for calculating target statistics (see the manual).\cr
  The value of the changing covariate for wave \code{m} is supposed
  in the simulations to be valid in the whole period from wave \code{m}
  to wave \code{m+1}.
  If the data set has \code{M} waves, this means that the values,
  if any, for wave \code{M} will not be used. Therefore, the number of columns
  can be \code{M} or \code{M-1}; if the former, the values in the last
  column will not be used.
}
\value{
  Returns the covariate as an object of class "varCovar", in which
  form it can be used as an argument to \code{\link{sienaDataCreate}}.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{sienaDataCreate}}, \code{\link{coCovar}},
\code{\link{coDyadCovar}}, \code{\link{varDyadCovar}},
\code{\link{sienaNodeSet}}}
\examples{

myvarCovar <- varCovar(s50a)
senders <- sienaNodeSet(50, nodeSetName="senders")
receivers <- sienaNodeSet(30, nodeSetName="receivers")
senders.covariate <- varCovar(s50a, nodeSet="senders")
receivers.covariate <- varCovar(s50s[1:30,], nodeSet="receivers")
}
\keyword{classes}
