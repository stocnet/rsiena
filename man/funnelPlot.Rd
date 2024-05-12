\name{funnelPlot}
\alias{funnelPlot}
\title{Plot function for a list of sienaFit objects}
\usage{
funnelPlot(anslist, k, threshold=NULL, origin=TRUE,
           plotAboveThreshold=TRUE, verbose=TRUE, ...)
}
\description{
  Draws a funnel plot for a list of sienaFit objects that all  
  have estimated the same parameter.
}
\arguments{
  \item{anslist}{A list of object of class \code{sienaFit}.}
  \item{k}{The number of the parameter to be plotted.}
  \item{threshold}{threshold for standard errors:
     all estimations where the standard error for parameter
     \code{k} is larger than this threshold will be disregarded.}
  \item{origin}{Boolean: whether to include the origin in the plot,
    if all estimates have the same sign.}
  \item{plotAboveThreshold}{Boolean: whether to include the estimates
    for which the standard error is larger than \code{threshold},
    and plot them with an asterisk at \code{se=threshold}.}
  \item{verbose}{Boolean: whether to report in the console all estimates
    omitted, because either their standard error is larger than
    \code{threshold}, or they were fixed.}
  \item{\dots}{For extra arguments (passed to \code{plot}).}
}
\details{
 The function \code{funnelPlot} plots estimates against standard
 errors for a given effect \code{k}, with red
 reference lines added at the two-sided
 significance threshold 0.05. Effects for which a score test
 was requested are not plotted (and reported to the console
 if \code{verbose}).\cr
 If not all effects with number \code{k} are the same
 in all \code{sienaFit} objects, a warning is given.
 The effect name for the first object is used as the plot title.\cr
 Another funnel plot is available as \code{\link{print.sienaMeta}}.
}
\value{
 The two-column matrix of
 values of the plotted points is invisibly returned.\cr
}
\author{Tom Snijders
}
\seealso{\code{\link{siena08}}, \code{\link{print.sienaMeta}}}

\examples{
# A meta-analysis for three groups does not make much sense.
# But using three groups shows the idea.

Group1 <- sienaDependent(array(c(N3401, HN3401), dim=c(45, 45, 2)))
Group3 <- sienaDependent(array(c(N3403, HN3403), dim=c(37, 37, 2)))
Group4 <- sienaDependent(array(c(N3404, HN3404), dim=c(33, 33, 2)))
dataset.1 <- sienaDataCreate(Friends = Group1)
dataset.3 <- sienaDataCreate(Friends = Group3)
dataset.4 <- sienaDataCreate(Friends = Group4)
OneAlgorithm <- sienaAlgorithmCreate(projname = NULL, nsub=1, n3=50, seed=123)
effects.1 <- getEffects(dataset.1)
effects.3 <- getEffects(dataset.3)
effects.4 <- getEffects(dataset.4)
ans.1 <- siena07(OneAlgorithm, data=dataset.1, effects=effects.1, batch=TRUE)
ans.3 <- siena07(OneAlgorithm, data=dataset.3, effects=effects.3, batch=TRUE)
ans.4 <- siena07(OneAlgorithm, data=dataset.4, effects=effects.4, batch=TRUE)
funnelPlot(list(ans.1, ans.3, ans.4), k=2)
funnelPlot(list(ans.1, ans.3, ans.4), k=2, origin=FALSE)
}
\keyword{ plot }
