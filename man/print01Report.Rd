\name{write_report}
\alias{write_report}
\alias{print01Report}
\alias{write_report.sienadata}
\title{
Function to produce a basic descriptive report of SIENA data sets}
\description{
Prints a report of a Siena data object and its default effects.\cr
  The functions \code{write_report} and \code{print01Report}
  are identical. \cr
  The first name is preferred, 
  the other is kept for backward compatibility.
}
\usage{
\method{write_report}{sienadata}(x, outputName=NULL, ...)

print01Report(data, modelname="Siena", getDocumentation=FALSE)
}
\arguments{
  \item{x}{a Siena data object of class \code{\link{sienadata}}
}
  \item{data}{a Siena data object of class \code{\link{sienadata}}
}
\item{outputName}{If \code{outputName = NULL} and the Siena data set
   has name \code{dataName}, output will be written by
   \code{write_report} to a file with name "\code{dataName}_report.txt"
   in the working directory.\cr
    If \code{outputName} is not \code{NULL}, it should be a character string
    without embedded spaces, and output will be written
    to a file with name "\code{outputName}_report.txt".}
\item{modelname}{
  Character string used to name the output file "\code{modelname}.txt"
}
\item{getDocumentation}{Flag to allow documentation of internal
  functions, not for use by users.}
\item{...}{
   Additional arguments, not used now.
}
}
\details{
  First deletes any file of the given name, then writes a new one.
  }
\value{
  No value returned.
}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley, Tom Snijders }

\examples{
mynet1 <- as_dependent_rsiena(
             array(c(s501, s502, s503), dim=c(50, 50, 3)))
mydata <- make_data_rsiena(mynet1)
write_report(mydata, outputName=tempfile())
}
\keyword{print}
