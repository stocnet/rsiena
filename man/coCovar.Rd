\name{coCovar}
\alias{coCovar}
\title{Function to create a constant covariate object}
\description{
  This function creates a constant covariate object from a vector.  
  }
\usage{
coCovar(val, centered=TRUE, nodeSet="Actors", warn=TRUE, imputationValues=NULL)
}
\arguments{
  \item{val}{Vector of covariate values}
  \item{centered}{Boolean: if TRUE, then the mean value is subtracted.}
  \item{nodeSet}{Name of node set: character string.
  If the entire data set contains more than one node set,
  then the node sets must be specified in all data objects.}
  \item{warn}{Logical: is a warning given if all values are \code{NA},
  or all non-missing values are the same.}
  \item{imputationValues}{Vector of covariate values of same
  length as \code{val}, to be used for imputation of
  \code{NA} values (if any) in \code{val}. Must not contain
  any \code{NA}.}
}
\details{
  When part of a Siena data object,
  the covariate is associated with the node set
  \code{nodeSet} of the Siena data object.
  In practice, the node set needs to be specified only in the case of
  the use of the covariate with a two-mode network.\cr
  If there are any \code{NA} values in \code{val}, and
  \code{imputationValues} is given, then the corresponding elements
  of \code{imputationValues} are used for imputation.
  If \code{imputationValues} is \code{NULL}, imputation
  is by the mean value. In both cases, cases with imputed values
  are not used for calculating target statistics (see the manual).
}
\value{
  Returns the covariate as an object of class "coCovar", in which
  form it can be used as an argument to \code{\link{sienaDataCreate}}.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{sienaDataCreate}}, \code{\link{varCovar}},
\code{\link{coDyadCovar}}, \code{\link{varDyadCovar}},
\code{\link{sienaNodeSet}}}
\examples{

myconstCovar <- coCovar(s50a[,1])
senders <- sienaNodeSet(50, nodeSetName="senders")
receivers <- sienaNodeSet(30, nodeSetName="receivers")
senders.attribute <- coCovar(rep(1:10, each=5), nodeSet="senders")
receivers.attribute <- coCovar(rep(1:5, each=6), nodeSet="receivers")
}
\keyword{classes}
