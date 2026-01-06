\name{make_constraint}
\alias{make_constraint}
\alias{make_constraint.sienadata}
\alias{make_constraint.sienaGroup}
\alias{sienaDataConstraint}
\title{
Function to change the values of the constraints between networks
}
\description{This function allows the user to change the constraints   
  of "higher", "disjoint" and "atLeastOne" for a specified pair of  
  networks in a Siena data object.\cr  
  The functions \code{make_constraint} and \code{sienaDataConstraint}
  are identical. The first name is preferred,
  the other is kept for backward compatibility.\cr
  The functions \code{make_constraint} and 
  \code{sienaDataConstraint} are identical. 
  The first name is preferred, 
  the other is kept for backward compatibility.
}
\usage{
\method{make_constraint}{sienadata}(x, net1, net2, type, value = FALSE, ...)
\method{make_constraint}{sienaGroup}(x, net1, net2, type, value = FALSE, ...)

sienaDataConstraint(x, net1, net2, type, value = FALSE)
}
\arguments{
  \item{x}{\code{sienadata} object or a \code{sienaGroup} object}
  \item{net1}{name of first network}
  \item{net2}{name of second network}
  \item{type}{one of "higher", "disjoint", "atleastOne".}
  \item{value}{Boolean giving the value.}
  \item{...}{Additional arguments, not used now.}
}
\details{
  The value of the appropriate attribute is set to the value requested.
  Note that, for \code{value=TRUE}, the correspondence of this value 
  to the data is not checked.
}
\value{
  Updated \code{sienadata} object.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{make_data_rsiena}}, \code{\link{make_group_rsiena}}}
\examples{
nowFriends <- as_dependent_rsiena(
                array(c(s501, s502, s503), dim=c(50, 50, 3)))
ever <- array(c(s501, s502, s503), dim=c(50, 50, 3))
ever[,,2] <- pmax(ever[,,1], ever[,,2])
ever[,,3] <- pmax(ever[,,2], ever[,,3])
everFriends <- as_dependent_rsiena(ever)
# Note: this data set serves to illustrate this function,
# but it is not an appropriate data set for estimation by siena07,
# because everFriends (for the three waves together) 
# depends deterministically on nowFriends (for the three waves together).
nowOrEver <- make_data_rsiena(nowFriends, everFriends)
attr(nowOrEver, "higher")
nowOrEver
nowOrEver.unconstrained <-
 make_constraint(nowOrEver, everFriends, nowFriends, "higher", FALSE)
nowOrEver.unconstrained
attr(nowOrEver.unconstrained, "higher")
}
\keyword{methods}
