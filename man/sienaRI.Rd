\name{sienaRI}
\alias{sienaRI}
\alias{entropy}
\alias{plot.sienaRI}
\alias{print.sienaRI}
\alias{sienaRI.methods}
\title{Functions to assess the relative importance of effects at observation
       moments}
\description{
 The function \code{sienaRI} computes various effect sizes for   
 Stochastic Actor-oriented Models.    
 As its value it returns the relative importance of effects   
 of a SAOM according to the measure of relative importance described in   
 Section 3.1 of Indlekofer and Brandes (2013).   
 This is based on the influence of effects on potential tie change or   
 behavior change decisions of individual actors at the given observation  
 moments (waves).  

 It takes the data as well as the complete model specification into account.   
 Therefore, necessary arguments are the analysed data given as a  
 \code{\link{siena}} data object as well as the complete model  
 specification represented either by an estimated \code{sienaFit}  
 object or by the triple consisting of a  
 suitable parameter vector \code{theta} and the corresponding   
 \code{\link{sienaAlgorithm}} and \code{\link{sienaEffects}} objects.\cr  
 Other measures useful for effect sizes are entropy-based  
 effect sizes as in Snijders (2004) and the within-ego  
 standard deviations of change statistics.  
 If \code{getChangeStats=TRUE}, the arrays of change statistics are   
 stored in the \code{sienaRI} object.
}
\usage{
sienaRI(data, ans=NULL, theta=NULL, algorithm=NULL, effects=NULL,
        getChangeStats=FALSE)
\method{print}{sienaRI}(x, printSigma=FALSE,\dots)
\method{plot}{sienaRI}(x, actors = NULL, col = NULL, addPieChart = FALSE,
              radius = 1, width = NULL, height = NULL, legend = TRUE,
              legendColumns = NULL, legendHeight = NULL,
              cex.legend = NULL, cex.names = NULL,\dots)
}
\arguments{
  \item{data}{\code{siena} data object representing the analyzed data and
              resulting from a call to \code{\link{sienaDataCreate}}.}
  \item{ans}{\code{\link{sienaFit}} object resulting from a call to
              \code{\link{siena07}}. The \code{sienaFit} object contains
              all necessary information on the model specification,
              in particular, the vector of parameter values \code{ans$theta},
              the used algorithm for estimation \code{ans$x},
              and information on included model effects \code{ans$effects}.
              If \code{ans} is a valid \code{sienaFit} object, the calculations
              of relative importances are based on \code{ans$theta},
              \code{ans$x}, and \code{ans$effects}.
              Alternatively, the necessary information can be given directly
              as a suitable parameter vector \code{theta}, a
              \code{sienaAlgorithm} object, and a \code{sienaEffects} object.
              In this case, \code{ans} has to be unspecified
              (i.e., \code{ans=NULL}).}
  \item{theta}{Vector of parameter values of effects included in the model.
               Length of \code{theta} has to be equal to the number of
               included effects.}
  \item{algorithm}{\code{sienaAlgorithm} object as resulting from a call to
               \code{\link{sienaAlgorithmCreate}}. Works only for estimation
               by Method of Moments (i.e., if  \code{maxlike = FALSE}).}
  \item{effects}{\code{sienaEffects} object specifying which effects are
               included the model. Note that \code{sienaRI} does not yet work
               for endowment or creation effect, i.e., included effects have to
               be of type \code{eval} (or \code{rate}).}
  \item{getChangeStats}{Boolean: If \code{TRUE}, an array of change statistics
               is added to the \code{sienaRI} object.}
  \item{x}{\code{sienaRI} object resulting from a call to \code{sienaRI}.}
  \item{printSigma}{Boolean: If \code{TRUE},
       average within-actor standard deviations of change statistics
       ('\code{sigma}') are included in the print.}
  \item{actors}{vector of integers: set of actors to be included in the plot;
       if \code{NULL}, all actors.}
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
\details{
\code{sienaRI} takes the data as well as the complete model specification
 into account. Therefore, necessary arguments are the analyzed data given as a
 \code{siena} data object as well as the complete model specification represented
 either by an estimated \code{sienaFit} object or by the triple consisting of a
 suitable parameter vector \code{theta} and the corresponding
 \code{sienaAlgorithm} and \code{sienaEffects} objects.

A \code{sienaFit} object contains all necessary information on the model
specification, in particular, the vector of parameter values \code{ans$theta},
the used algorithm for estimation \code{ans$x}, and information on included
model effects \code{ans$effects}.

If \code{ans} is a valid \code{sienaFit} object, the calculations of relative
importances are based on \code{ans$theta}, \code{ans$x}, and
\code{ans$effects}. Alternatively, the necessary information can be given
directly as a suitable parameter vector \code{theta}, a \code{sienaAlgorithm}
object, and a \code{sienaEffects} object. In this case, \code{ans} has to be
unspecified, i.e., \code{ans=NULL}.

Note that \code{sienaRI} works only with Method of Moments
(i.e., for \code{sienaAlgorithm} objects with \code{maxlike = FALSE})
and that it does not yet work for endowment or creation effects
(i.e., included effects have to be of type \code{eval}),
and also not for models with interaction effects.
For two-mode (bipartite) networks as dependent variables,
it works only if the number of second-mode nodes is less than the
number of actors.
If the network is non-directed, the relative importances and entropy-based
'degrees of certainty' are calculated for \code{modelType=2}
('forcing'; see  \code{\link{sienaAlgorithmCreate}}).
For dependent behavior variables, \code{behModelType=1} ('standard') is assumed.

If there are any missing tie values in the network data set, they are imputed by
initial zeros and Last Observation Carried Forward. Structural zeros and ones
are replaced by \code{NA} and treated as impossible choices in the probability
vectors and ignored in the standard deviations; but the change statistics
for these dyads still are given in \code{changeStatistics} (if requested).

The averages reported in the components \code{sigmas} (average across actors)
and \code{meansigmas} (average across waves) are obtained by averaging at the
variance level and then taking square roots.
}
\value{
 If the model contains only one dependent variable, \code{sienaRI} returns
 an object of class \code{sienaRI}. Otherwise, it returns a list of objects
 of class \code{sienaRI}, each corresponding to one dependent variable.

 A returned \code{sienaRI} object stores the expected relative importances of
 effects of one dependent variable at observation moments as defined in
 Section 3.1 of Indlekofer and Brandes (2013).

A \code{sienaRI} object is a list with the following components.
For the components referred to as lists themselves, these are lists
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
          this is R_H(i,x) of formula (6) in Snijders (2004).
          The mean over actors of these degrees of certainty, given by
          formula (7) in Snijders (2004), is printed by the
          \code{print} method for \code{sienaRI} objects.}
    \item{- \code{sigma}}{a list of effects by ego matrices of the values of
          the within-ego standard deviations of the change statistics.}
    \item{- \code{sigmas}}{effects by wave matrices of averages
         (see above for how this is done; across actors) of \code{sigma}.
          These are printed if \code{printSigma=TRUE}.}
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
         when the data are as in this wave.}
    }
}

\references{
 Indlekofer, N. and Brandes, U. (2013),
 Relative Importance of Effects in Stochastic Actor-oriented Models.
 \emph{Network Science}, \bold{1}, 278--304.
 
 Snijders, T.A.B. (2004),
  Explained Variation in Dynamic Network Models.
  \emph{Mathematics and Social Sciences}, \bold{168}, 31--41.
}
\author{Natalie Indlekofer, additions by Tom Snijders}
\examples{
myalgorithm <- sienaAlgorithmCreate(nsub=1, n3=50, projname=NULL)
mynet1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, density, recip, outAct, inPop)
(myeff <- setEffect(myeff, reciAct, parameter=1))
ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE)
RI <- sienaRI(mydata, ans)
RI
print(RI, printSigma=TRUE)
# average within-ego standard deviations of change statistics by wave:
RI$sigmas
# sigma averaged over waves:
RI$meansigmas
# semi-standardized parameters by wave:
ans$theta * RI$meansigmas
\dontrun{
plot(RI, addPieChart=TRUE)
plot(RI, actors=1:20, addPieChart=TRUE, radius=1.08)
}
myalgorithm <- sienaAlgorithmCreate(nsub=1, n3=50, projname=NULL)
mynet2 <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
mybeh <- sienaDependent(s50a[,2:3], type="behavior")
mydata2 <- sienaDataCreate(mynet2, mybeh)
myeff2 <- getEffects(mydata2)
myeff2 <- includeEffects(myeff2, density, recip, transTies)
ans2 <- siena07(myalgorithm, data=mydata2, effects=myeff2, batch=TRUE)
# Use only the parameters for the evaluation function:
theta.eval <- ans2$theta[myeff2$type[myeff2$include]=="eval"]
RI <- sienaRI(mydata2, theta=theta.eval, algorithm=myalgorithm,
              effects = myeff2)
RI
\dontrun{
plot(RI[[2]], col = c("red", "green"), legend=FALSE)
plot(RI[[1]], addPieChart = TRUE, legendColumns=2)
}
}
\keyword{models}
