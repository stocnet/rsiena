\name{as_nodeset_rsiena}
\alias{as_nodeset_rsiena}
\alias{sienaNodeSet}
\title{Function to create a node set }
\description{Creates a Siena node set which can be used as the nodes in a  
  Siena network.\cr
  The functions \code{as_nodeset_rsiena} and \code{sienaNodeSet}
  are identical. The first name is preferred, 
  the other is kept for backward compatibility.
}
\usage{
as_nodeset_rsiena(n, nodeSetName="Actors", names=NULL)

sienaNodeSet(n, nodeSetName="Actors", names=NULL)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{n}{integer, size of set.}
  \item{nodeSetName}{character string naming the node set.}
  \item{names}{optional character string vector of length n of the names of the
 nodes.}
}
\details{
  This function is important for data sets having more than one node set,
  but not otherwise.
}
\value{
  Returns a Siena node set, an integer vector, possibly with names, plus the
  attributes, class equal to "as_nodeset_rsiena", and nodeSetName equal to
  the argument nodeSetName.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/ }}
\author{Ruth Ripley }
\seealso{\code{\link{as_dependent_rsiena}}, \code{\link{make_data_rsiena}}}
\examples{
senders <- as_nodeset_rsiena(50, nodeSetName="senders")
receivers <- as_nodeset_rsiena(30, nodeSetName="receivers")
senders.attribute <- as_covariate_rsiena(rep(1:10, each=5), 
                        type="monadic", nodeSet="senders")
receivers.attribute <- as_covariate_rsiena(rep(1:5, each=6), 
                        type="monadic", nodeSet="receivers")
mynet <- as_dependent_rsiena(
      array(c(s501[,1:30], s502[,1:30]), dim=c(50, 30, 2)),
      nodeSet=c("senders", "receivers"))
(mydata <- make_data_rsiena(mynet, senders.attribute, receivers.attribute,
      nodeSets=list(senders, receivers)))
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{classes }
