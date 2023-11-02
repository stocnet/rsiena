\name{sienaNodeSet}
\alias{sienaNodeSet}
\title{Function to create a node set }
\description{Creates a Siena node set which can be used as the nodes in a  
  Siena network.
}
\usage{
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
  attributes, class equal to "sienaNodeSet", and nodeSetName equal to
  the argument nodeSetName.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/ }}
\author{Ruth Ripley }
\seealso{\code{\link{sienaDependent}}, \code{\link{sienaDataCreate}}}
\examples{
senders <- sienaNodeSet(50, nodeSetName="senders")
receivers <- sienaNodeSet(30, nodeSetName="receivers")
senders.attribute <- coCovar(rep(1:10, each=5), nodeSet="senders")
receivers.attribute <- coCovar(rep(1:5, each=6), nodeSet="receivers")
mynet <- sienaDependent(array(c(s501[,1:30], s502[,1:30]), dim=c(50, 30, 2)),
      nodeSet=c("senders", "receivers"))
(mydata <- sienaDataCreate(mynet, senders.attribute, receivers.attribute,
      nodeSets=list(senders, receivers)))
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{classes }
