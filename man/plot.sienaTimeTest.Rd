\name{plot.sienaTimeTest}
\alias{plot.sienaTimeTest}
\title{Functions to plot assessment of time heterogeneity of parameters}
\description{
 Plot method for \code{\link{sienaTimeTest}} objects.
 }

\usage{
\method{plot}{sienaTimeTest}(x, pairwise=FALSE, effects,
     scale=.2, plevels=c(.1, .05, .025), ...)
}
\arguments{
  \item{x}{ A sienaTimeTest object returned by \code{sienaTimeTest}. }
  \item{pairwise}{ A Boolean value corresponding to whether the user
  would like a pairwise plot of the simulated statistics to assess
  correlation among the effects (\code{pairwise=TRUE}), or a plot of the
  estimates across waves in order to assess graphically the results of
  the score type test. }
\item{effects}{ A vector of integers corresponding to the indices given
  in the \code{sienaTimeTest} output for effects which are to be
  plotted.}
\item{scale}{ A positive number corresponding to the number of standard
  deviations on one step estimates to use for computing the maximum and
  minimum of the plotting range. We recommend experimenting with this
  number when the y-axes of the plots are not satisfactory. Smaller
  numbers shrink the axes.}
\item{plevels}{ A list of three decimals indicating the gradients at
  which to draw the confidence interval bars.}
 \item{\dots}{For extra arguments. The Lattice parameter \code{layout}
   can be used to control the layout of the graphs.}
}
\details{
 The \code{pairwise=TRUE} plot may be used to assess whether effects are
 highly correlated. This information may be important when considering
 forward-model selection, since highly correlated effects may have
 highly correlated one-step estimates, particularly since the individual
 score type tests are not orthogonalized against the scores and
 deviations of yet-unestimated dummies.  For example, reciprocity and
 outdegree may have highly correlated statistics as indicated by a
 strong, positive correlation coefficient.  When considering whether to
 include dummy terms, it may be a good idea to include, e.g., outdegree,
 estimate the parameter, and see whether reciprocity dummies remain
 significant after method of moments estimation of the updated model--as
 opposed to including both outdegree and reciprocity.

 The \code{pairwise=FALSE} plot displays the most of the information
 garnered from \code{sienaTimeTest} in a graphical fashion.  For each
 effect, the method of moments parameter estimate for the base period
 (i.e. wave 1) is given as a blue, horizontal reference line. One step
 estimates are given for all of the parameters by dots at each wave. The
 dots are colored black if the parameter has been included in the model
 already (i.e. has been estimated via method of moments), or red if they
 have not been included. Confidence intervals are given based on pivots
 given at \code{pvalues}. Evidence of time heterogeneity is suggested by
 points with confidence intervals not overlapping with the base period.
}
\value{None
}

\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}
  for general information on RSiena.

}
\author{Josh Lospinoso}
\seealso{\code{\link{siena07}}, \code{\link{sienaTimeTest}},
  \code{\link[lattice]{xyplot}}}
\examples{
\dontrun{
myalgorithm <- sienaAlgorithmCreate(nsub=2, n3=500)
# It makes no sense to put together the following data set,
# but just for demonstration:
mynet1 <- sienaDependent(array(c(s501, s502, s503, s501, s503, s502), dim=c(50, 50, 6)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip)
myeff <- includeTimeDummy(myeff, density, timeDummy="all")
myeff <- includeTimeDummy(myeff, recip, timeDummy="2,3,5")
myeff <- includeTimeDummy(myeff, transTrip, timeDummy="2,3")
(ansp <- siena07(myalgorithm, data=mydata, effects=myeff))
ttp <- sienaTimeTest(ansp)
summary(ttp)

## Pairwise plots show
plot(ttp, pairwise=TRUE)

## Time test plots show
plot(ttp, effects=1:3) ## default layout
plot(ttp, effects=1:3, layout=c(3,1))
}
}

\keyword{methods}
