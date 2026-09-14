\name{interpret_size}
\alias{interpret_size}
\alias{sienaRI}
\alias{siena_sizes}
\alias{interpret_size.sienaFit}
\alias{interpret_size.sienaEffects}
\alias{entropy}
\alias{plot.siena_sizes}
\alias{print.siena_sizes}
\alias{summary.siena_sizes}
\alias{plot.sienaRI}
\alias{print.sienaRI}
\alias{siena_sizes.methods}
\title{Function to assess  effect sizes for \code{siena} results}
\description{
 The function \code{interpret_size} computes various effect sizes for
 Stochastic Actor-oriented Models.
 It computes estimates of the semi-standardized parameters
 treated in Section 6.5.1 of Snijders and Steglich (2026), 
 based on the within-ego standard deviations of change statistics;
 the entropy-based degrees of certainty of  Snijders (2004),
 treated in Section 6.5.3 of Snijders and Steglich (2026); 
 and the relative importance of effects according to the measure 
 proposed in Section 3.1 of Indlekofer and Brandes (2013).\cr
 The return object is of class \code{siena_sizes}.
}
\usage{
\method{interpret_size}{sienaFit}(x, data, getChangeStats=FALSE, ...)
\method{print}{siena_sizes}(x, \dots)
\method{summary}{siena_sizes}(object, \dots)
\method{plot}{siena_sizes}(x, depvar=NULL, actors=NULL, col=NULL, 
              addPieChart=FALSE, radius=1, width=NULL, height=NULL, 
              legend=TRUE, legendColumns=NULL, legendHeight=NULL,
              cex.legend = NULL, cex.names = NULL,\dots)
}
\arguments{
  \item{x}{for \code{interpret_size.sienaFit}: a \code{sienaFit} object,
           resulting from a call of \code{\link{siena}};\cr
           for \code{print.siena_sizes} and \code{plot.siena_sizes}:
           a \code{siena_sizes} object resulting
           from a call to \code{interpret_size}.}
  \item{object}{a \code{siena_sizes} object resulting
           from a call to \code{interpret_size}.}
  \item{data}{\code{siena} data object representing the analyzed data and
              resulting from a call to \code{\link{make_data_rsiena}}.}
  \item{getChangeStats}{Boolean: If \code{TRUE}, an arrays of change statistics
               and a list of toggle probabilities
               are added to the \code{interpret_size} object.}
  \item{actors}{vector of integers: set of actors to be included in the plot;
       if \code{NULL}, all actors.}
  \item{depvar}{string: name of the dependent variable;
       if \code{NULL}, the first dependent variabe in the data set.}
  \item{col}{Colors used in the plot. If \code{col=NULL} a default color
       scheme is used.}
  \item{addPieChart}{Boolean: If \code{TRUE}, pie charts of aggregated relative
       importances for the complete set of actors will be added to the plot.}
  \item{radius}{Radius of pie charts. Only effective if \code{addPieCharts = TRUE}.}
  \item{width}{Width of the plot. If \code{width=NULL} a default value is used.}
  \item{height}{Height of the plot. If \code{height=NULL} a default value is used.}
  \item{legend}{Boolean: if \code{TRUE} a legend is added to the plot.
        \code{x$effectNames} are used as labels.}
  \item{legendColumns}{Number of columns in legend. If \code{legendColumns=NULL}
        a default value is used. Only effective if \code{legend=TRUE}.}
  \item{legendHeight}{Height of legend. If \code{legendHeight=NULL} a default
        value is used. Only effective if \code{legend=TRUE}.}
  \item{cex.legend}{Specifies the relative font size of legend labels.}
  \item{cex.names}{Specifies the relative font size of bar graph labels.}
  \item{\dots}{Other arguments.}
}
\details{\code{interpret_size} works only for estimation by the Method of Moments
      using \code{\link{modelType}} 1 (directed networks) or 2
      ('forcing' model for non-directed networks).
      For dependent behavior variables, \code{behModelType=1} ('standard') is assumed.
      It does not yet work for endowment or creation effects
      (i.e., included effects have to be of type \code{eval} or \code{rate})
    and also not for models with interaction effects.

    The averages reported in the components \code{sigmas} (average across actors)
    and \code{meansigmas} (average across waves) are obtained by averaging at the
    variance level and then taking square roots.
    
    If \code{getChangeStats=TRUE}, the arrays of change statistics 
    and toggle probabilities are stored in the created object.
    
    The \code{print} method prints a table of estimates with their 
    \emph{p}-values, the semi-standardized parameter estimates,
    and the degrees of certainty.\cr
    The \code{summary} method prints also the standard deviations 
    of the change statistics and the expected relative importances.\cr
    The \code{plot} method plots the expected relative importances.
}
\value{
 \code{interpret_size} returns an object of class \code{siena_sizes}.
 This is a list of which the first elements are each for one dependent variable. 

  Each component that is a list for one dependent variable 
  is itself a list with the following components.
For the components referred to below as lists, these are lists
corresponding to the observation moments.
    \describe{
    \item{- \code{dependentVariable}}{the name of the corresponding
          dependent variable.}
    \item{- \code{effectNames}}{the names of considered effects.}
    \item{- \code{RIActors}}{a list that contains the expected
          relative importances of effects
          for each potential actor decision at observation moments.
          This is equation (3) in Indlekofer and Brandes (2013).}
    \item{- \code{expectedRI}}{a list that contains the expected
          relative importances of effects aggregated over all actors for each
          network observation. These are the averages of the actor related
          values in \code{RIActors}.
          This is equation (4) in Indlekofer and Brandes (2013).}
    \item{- \code{IActors}}{a list that contains the expected
          importances of effects for each potential actor decision at
          observation moments. This is the numerator of equation (3)
          in Indlekofer and Brandes (2013).}
    \item{- \code{expectedI}}{a list that contains the expected
          importances of effects
          aggregated over all actors in each observation.
          More precisely, it contains the averages of the actor related values
          in \code{IActors}.}
    \item{- \code{absoluteSumActors}}{a list that contains the sum of
          the (unstandardized) L1-differences
          calculated for each potential actor decision at observation moments.
          This is the denominator of equation (3)
          in Indlekofer and Brandes (2013).}
    \item{- \code{RHActors}}{a list that contains the degree of certainty,
          also called degree of determination, in
          the potential ministep taken by an actor at the observation moments;
          this is R_H(i,beta,x) of formula (6.14) in Snijders and Steglich (2026)
          and R_H(i,x) of formula (6) in Snijders (2004).
          The mean over actors of these degrees of certainty, given by
          formula (6.15) in Snijders and Steglich (2026), is printed by the
          \code{print} method for \code{interpret_size} objects.}
    \item{- \code{sigma}}{a list of effects by ego matrices of the values of
          the within-ego standard deviations of the change statistics.
          These are the square roots of formula (6.8) 
          in Snijders and Steglich (2026).  }     
    \item{- \code{sigmas}}{effects by wave matrices of averages
         (see above for how this is done; across actors) of \code{sigma}.}
    \item{- \code{meansigmas}}{average (see above for how this is done)
          over waves of \code{sigmas}.}
    \item{- \code{changeStatistics}}{a list of arrays
          (effects by choices by egos)
          containing, for each observation wave,
          the values of the change statistics for making this choice,
          where for one-mode networks the choice is defined by the alters
          with alter=ego referring to 'no change',
          for two-mode networks the choice is defined by the second-mode nodes
          with the last choice referring to 'no change',
          and for behavior the choice is defined as
          going down, staying constant, or going up;
          this output is produced only if \code{getChangeStats=TRUE}.}
    \item{- \code{toggleProbabilities}}{an array
         (egos by choices by waves), where "choices" are as directly above,
         giving the choice probabilities of ego in a ministep,
         when the data are as in this wave;
         this output is produced only if \code{getChangeStats=TRUE}.}
    }
}

\references{
 Indlekofer, N. and Brandes, U. (2013),
 Relative Importance of Effects in Stochastic Actor-oriented Models.
 \emph{Network Science}, \bold{1}, 278--304.

 Snijders, T.A.B. (2004),
  Explained Variation in Dynamic Network Models.
  \emph{Mathematics and Social Sciences}, \bold{168}, 31--41.
  
  Snijders, Tom A.B., and Steglich, Christian (2026),
  Stochastic Actor-Oriented Models for Longitudinal Networks,
  \emph{Cambridge University Press}.
}
\author{Tom Snijders, Daniel Gotthardt, Natalie Indlekofer}
\examples{
myalgo <- set_algorithm_saom(nsub=1, n3=50, seed=1293)
mynet1 <- as_dependent_rsiena(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- make_data_rsiena(mynet1)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, list(density, recip, outAct, inPop))
(myeff <- set_effect(myeff, reciAct, parameter=1))
ans <- siena(mydata, effects=myeff, control_algo=myalgo, batch=TRUE)
RI <- interpret_size(ans, mydata)
RI
# average within-ego standard deviations of change statistics by wave:
RI[["mynet1"]]$sigmas
# sigma averaged over waves:
RI$meansigmas
# to get the change statistics:
RI.cs <- interpret_size(ans, mydata, getChangeStats=TRUE)
# For the network at the first wave, the dimensions of the change statistics:
dim(RI.cs[["mynet1"]]$changeStatistics[[1]])
# E.g., the change statistics for the first actor (note the 0 first column):
RI.cs[["mynet1"]]$changeStatistics[[1]][,1,]
# semi-standardized parameters by wave:
ans$theta * RI[["mynet1"]]$meansigmas
\dontrun{
plot(RI[["mynet1"]], addPieChart=TRUE)
plot(RI[["mynet1"]], actors=1:20, addPieChart=TRUE, radius=1.08)
}
# For two dependent variables:
myalgo <- set_algorithm_saom(nsub=1, n3=50, seed=1293)
mynet2 <- as_dependent_rsiena(array(c(s502, s503), dim=c(50, 50, 2)))
mybeh <- as_dependent_rsiena(s50a[,2:3], type="behavior")
mydata2 <- make_data_rsiena(mynet2, mybeh)
myeff2 <- make_specification(mydata2)
myeff2 <- set_effect(myeff2, list(density, recip, transTies))
myeff2 <- set_effect(myeff2, avAlt, depvar="mybeh", covar1="mynet2")
ans2 <- siena(mydata2, effects=myeff2, control_algo=myalgo, batch=TRUE)
# Use only the parameters for the evaluation function:
(theta.eval <- coef(ans2, dropRates=TRUE))
RI <- interpret_size(ans2, data=mydata2)
RI
\dontrun{
plot(RI[["mybeh"]], col = c("red", "green"), legend=FALSE)
plot(RI[["mynet2"]], addPieChart = TRUE, legendColumns=2)
}
}
\keyword{models}
