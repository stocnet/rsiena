\name{edit.sienaEffects}
\alias{edit.sienaEffects}
\title{
  Allow editing of a sienaEffects object if a gui is available.
}
\description{
  Interactive editor for an effects object. A wrapper to edit.data.frame.
}
\usage{
\method{edit}{sienaEffects}(name, ...)
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{name}{An object of class \code{sienaEffects}}
  \item{\dots}{For extra arguments (none used at present)
}}
\details{Will be invoked by \code{fix(name)} for an object of class
  \code{sienaEffects}.
}
\value{
  The updated object. There is no backup copy, and the edits cannot be undone.
}
\seealso{\code{\link{make_specification}}}
\references{
  See \url{https://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Ruth Ripley
}
\examples{
mynet1 <- as_dependent_rsiena(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- as_dependent_rsiena(s50a, type="behavior")
mycovar <- as_covariate_rsiena(rnorm(50))
mydyadcovar <- as_covariate_rsiena(
                  matrix(as.numeric(rnorm(2500) > 2), nrow=50))
(mydata <- make_data_rsiena(mynet1, mybeh, mycovar, mydyadcovar))
myeff <- make_specification(mydata)
\dontrun{
fix(myeff)
}
}
\keyword{methods}
