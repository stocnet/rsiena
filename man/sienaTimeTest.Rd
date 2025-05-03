\name{sienaTimeTest}
\alias{sienaTimeTest}
\title{Functions to assess and account for time heterogeneity of parameters}
\description{
 Takes a \code{sienaFit} object estimated by Method of Moments,  
 and tests for time heterogeneity by the addition of interactions with  
 time dummy variables at waves \code{m=2...(M-1)}.  
 The test used is the score-type test of Schweinberger (2012).

 Tests for joint significance, parameter-wise significance,  
 period-wise significance, individual significance, and one-step   
 estimates of the unrestricted model parameters are returned in a list.
 }
\usage{
sienaTimeTest(sienaFit, effects=NULL, excludedEffects=NULL, condition=FALSE)
}
\arguments{
  \item{sienaFit}{ A \code{sienaFit} object returned by \code{siena07}. }
  \item{effects}{Optional vector of effect numbers to test. Use the
  numbering on the print of the sienaFit object.}
  \item{excludedEffects}{Optional vector of effect numbers
  for which time heterogeneity is not to be tested. Use the
  numbering on the print of the sienaFit object.}
  \item{condition}{ Whether to orthogonalize effect-wise score-type tests
  and individual significance tests against estimated effects
  and un-estimated dummy terms, or just against estimated effects.}
}
\details{
  This test follows the score type test of Schweinberger (2012) as
 elaborated by Lospinoso et al. (2011) by using statistics already
 calculated at each wave to obtain vectors of partitioned moment
 functions corresponding to a restricted model (the model in the
 \code{sienaFit} object; used as null hypothesis) and an
 unrestricted model (which contains dummies for waves
 \code{m=2...(M-1)}; used as alternative hypothesis).

 \code{condition=TRUE} leads to a rough-and-easy approximation to
 controlling the mentioned tests also for the unestimated effects.

 After assessing time heterogeneity, effects objects can be modified by
 adding numbers of all or some periods to the \code{timeDummy} column.
 This is facilitated by the \code{\link{includeTimeDummy}} function.
 For an effects object in which the \code{timeDummy} column
 of some of the included effects includes some or all period numbers,
 interactions of those effects  with time dummies for the indicated
 periods will also be estimated.

 An alternative to the use of \code{\link{includeTimeDummy}}
 is to define time-dependent actor covariates (dummy variables or other
 functions of wave number that are the same for all actors),
 include these in the data set through \code{\link{sienaAlgorithmCreate}},
 and include interactions of other effects with ego effects of these
 time-dependent actor covariates by \code{\link{includeInteraction}}.
 This is illustrated in an example below.
 Using \code{\link{includeTimeDummy}} is easier;
 using self-defined interactions with time-dependent variables
 gives more control.

 If you wish to use this function with \code{sienaFit} objects that use the
 finite differences method of derivative estimation, or which use maximum
 likelihood estimation, you must request the derivatives to be returned
 by wave using the \code{byWave=TRUE} option for \code{\link{siena07}}.

 Effects leading to dummy interactions that are collinear with the
 model originally fitted, after excluding the effects mentioned,
 will be automatically excluded from the time heterogeneity testing.

 If \code{sienaTimeTest} gives errors that there are too
 many collinear effects, run it with a smaller set of effects as specified
 by the \code{effects} parameter.
 For example, if the model has 40 effects of which the
 first 8 are rate parameters and therefore uninteresting,
 and there is such an error message, try \code{effects=9:30};
 if that still does not work, decrease the upper limit of 30, if it does work
 increase it, to find the largest possible set of effects for which
 heterogeneity assessment still is possible; then as a next step
 try the remaining effects in a similar way.

 Also if the execution is time-consuming, e.g., for a multi-group
 \code{sienaFit} object with many groups and many effects,
 it can be helpful to carry out the function in smaller subsets of effects.
}
\value{
  \code{sienaTimeTest} returns a list containing many items,
    including the following:
  \item{JointTest }{A chi-squared test for joint significance of the dummies.}
  \item{EffectTest }{A chi-squared test for joint significance
      across dummies for each separate effect.}
  \item{GroupTest }{A chi-squared test for joint significance
      across dummies; if \code{sienaFit} is a fit for a multi-group object
      then these refer to each group; else they refer to each period.}
  \item{IndividualTest }{A matrix displaying initial estimates, one-step
    estimates, and \eqn{p}-values for the individual interactions.}
}

\references{
J.A. Lospinoso, M. Schweinberger, T.A.B. Snijders, and R.M. Ripley
(2011), Assessing and Accounting for Time Heterogeneity in Stochastic
Actor Oriented Models. \emph{Advances in Data Analysis and Computation},
\bold{5}, 147--176.

M. Schweinberger (2012), Statistical modeling of network panel data:
Goodness-of-fit. \emph{British Journal of Statistical
and Mathematical Psychology} \bold{65}, 263--281.
}

\author{Josh Lospinoso, Tom Snijders}
\seealso{\code{\link{siena07}}, \code{\link{plot.sienaTimeTest}},
  \code{\link{includeTimeDummy}}}
\examples{
## Estimate a restricted model
myalgorithm <- sienaAlgorithmCreate(nsub=1, n3=50, projname=NULL)
# Short estimation not for practice, just for having a quick demonstration
mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip)
ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE)

## Conduct the score-type test to assess whether heterogeneity is present.
tt <- sienaTimeTest(ans)
summary(tt)

## Suppose that we wish to include time dummies.
## Add them in the following way:
myeff <- includeTimeDummy(myeff, recip, transTrip, timeDummy="2")
ans2 <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE)

## Re-assess the time heterogeneity
(tt2 <- sienaTimeTest(ans2))

## And so on..

## A demonstration of the plotting facilities, on a larger dataset:
## (Of course pasting these identical sets of three waves after each other
##  in a sequence of six is not really meaningful. It's just a demonstration.)

myalgorithm <- sienaAlgorithmCreate(nsub=1, n3=50, seed=654, projname=NULL)
mynet1 <- sienaDependent(array(c(s501, s502, s503, s501, s503, s502),
                               dim=c(50, 50, 6)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip, balance)
myeff <- includeTimeDummy(myeff, density, timeDummy="all")
myeff <- includeTimeDummy(myeff, recip, timeDummy="2,3,5")
myeff <- includeTimeDummy(myeff, balance, timeDummy="4")
\dontrun{
(ansp <- siena07(myalgorithm, data=mydata, effects=myeff))
ttp <- sienaTimeTest(ansp, effects=1:4)

## Pairwise plots show
plot(ttp, pairwise=TRUE)

## Time test plots show
plot(ttp, effects=1:4, dims=c(2,2))
}

## Instead of working with includeTimeDummy,
## you can also define time dummies explicitly;
## this may give more control and more clarity:
dum2 <- matrix(c(0,1,0,0,0), nrow=50, ncol=5, byrow=TRUE)
dum3 <- matrix(c(0,0,1,0,0), nrow=50, ncol=5, byrow=TRUE)
dum4 <- matrix(c(0,0,0,1,0), nrow=50, ncol=5, byrow=TRUE)
dum5 <- matrix(c(0,0,0,0,1), nrow=50, ncol=5, byrow=TRUE)
time2 <- varCovar(dum2)
time3 <- varCovar(dum3)
time4 <- varCovar(dum4)
time5 <- varCovar(dum5)
mydata <- sienaDataCreate(mynet1, time2, time3, time4, time5)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip, balance)
## corresponding to includeTimeDummy(myeff, density, timeDummy="all"):
myeff <- includeEffects(myeff, egoX, interaction1='time2')
myeff <- includeEffects(myeff, egoX, interaction1='time3')
myeff <- includeEffects(myeff, egoX, interaction1='time4')
myeff <- includeEffects(myeff, egoX, interaction1='time5')
## corresponding to myeff <- includeTimeDummy(myeff, recip, timeDummy="2,3,5"):
myeff <- includeInteraction(myeff, egoX, recip, interaction1=c('time2', ''))
myeff <- includeInteraction(myeff, egoX, recip, interaction1=c('time3', ''))
myeff <- includeInteraction(myeff, egoX, recip, interaction1=c('time5', ''))
## corresponding to myeff <- includeTimeDummy(myeff, balance, timeDummy="4"):
myeff <- includeInteraction(myeff, egoX, balance, interaction1=c('time4', ''))
\dontrun{
(anspp <- siena07(myalgorithm, data=mydata, effects=myeff))
## anspp contains identical results as ansp above.
}

## A demonstration of RateX heterogeneity.
\dontrun{
mynet1 <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
myccov <- coCovar(s50a[,1])
mydata <- sienaDataCreate(mynet1, myccov)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip, balance)
myeff <- includeTimeDummy(myeff, RateX, type="rate", interaction1="myccov")
(ans <- siena07(myalgorithm, data=mydata, effects=myeff))
}
}
\keyword{tests}
