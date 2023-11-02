\name{sienaDataConstraint}
\alias{sienaDataConstraint}
\title{
Function to change the values of the constraints between networks.
}
\description{This function allows the user to change the constraints   
  of "higher", "disjoint" and "atLeastOne" for a specified pair of  
  networks in a Siena data object.
}
\usage{
sienaDataConstraint(x, net1, net2,
          type = c("higher", "disjoint", "atLeastOne"), value = FALSE)
}
\arguments{
  \item{x}{Siena data object; maybe a group object?}
  \item{net1}{name of first network}
  \item{net2}{name of second network}
  \item{type}{one of "higher", "disjoint", "atleastOne". Default is "higher".}
  \item{value}{Boolean giving the value.}
}
\details{
  The value of the appropriate attribute is set to the value requested.
  Note that, for \code{value=TRUE}, the correspondence of this value 
  to the data is not checked.
}
\value{
  Updated Siena data object.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{sienaDataCreate}}, \code{\link{sienaGroupCreate}}}

\examples{
nowFriends <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
ever <- array(c(s501, s502, s503), dim=c(50, 50, 3))
ever[,,2] <- pmax(ever[,,1], ever[,,2])
ever[,,3] <- pmax(ever[,,2], ever[,,3])
everFriends <- sienaDependent(ever)
# Note: this data set serves to illustrate this function,
# but it is not an appropriate data set for estimation by siena07,
# because everFriends (for the three waves together) depends deterministically
# on nowFriends (for the three waves together).
nowOrEver <- sienaDataCreate(nowFriends, everFriends)
attr(nowOrEver, "higher")
nowOrEver
nowOrEver.unconstrained <-
   sienaDataConstraint(nowOrEver, everFriends, nowFriends, "higher", FALSE)
nowOrEver.unconstrained
attr(nowOrEver.unconstrained, "higher")
}
\keyword{ classes }
