\name{xtable}
\alias{xtable}
\title{
Access xtable in package xtable}
\description{
Dummy function to allow access to xtable in package xtable
}
\usage{
xtable(x, ...)
}
\arguments{
  \item{x}{\code{\link{sienaFit}} object

}
\item{\dots}{
  Other arguments for \code{\link{xtable.sienaFit}}
}
}
\value{
  Value returned from \code{\link{xtable.sienaFit}}
}
\references{
\url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley}
\seealso{\code{\link{xtable.sienaFit}}}

\examples{

## The function is currently defined as
function (x, ...)
{
    xtable::xtable(x, ...)
  }
\dontrun{
myalgo <- set_algorithm_saom(nsub=2, n3=100, seed=1293)
mynet <- as_dependent_rsiena(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- make_data_rsiena(mynet)
myeff <- make_specification(mydata)
ans <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE)
ans
summary(ans)
xtable(ans, type="html", file="ans.html")
}
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{print}
