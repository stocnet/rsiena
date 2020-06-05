\name{print01Report}
\Rdversion{1.1}
\alias{print01Report}
\title{
Function to produce the Siena01 report from R objects}
\description{
Prints a report of a Siena data object and its default effects.
}
\usage{
print01Report(data, modelname = "Siena", getDocumentation=FALSE)
}
\arguments{
  \item{data}{a Siena data object
}
\item{modelname}{
  Character string used to name the output file "modelname.txt"
}
\item{getDocumentation}{Flag to allow documentation of internal
  functions, not for use by users.}
}
\details{
  First deletes any file of the name "modelname.txt", then prints a new one.
  }
\value{
  No value returned.
}
\references{See \url{http://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley }

\examples{
mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mydata <- sienaDataCreate(mynet1)
\dontrun{
print01Report(mydata, modelname="mydescription")
}
}
\keyword{print}
