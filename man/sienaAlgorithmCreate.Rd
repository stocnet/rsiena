\name{sienaAlgorithmCreate}
\alias{sienaModelCreate}
\alias{sienaAlgorithmCreate}
\alias{sienaModel}
\alias{sienaAlgorithm}
\alias{model.create}
\title{Function to create an object containing the algorithm
  specifications for parameter estimation in RSiena}
\description{
  Creates an object with specifications for the algorithm for parameter
  estimation in RSiena.\cr
  \code{sienaAlgorithmCreate()} and \code{sienaModelCreate()}
  are identical functions; the second name was
  used from the start of the \code{RSiena} package, but the first name
  indicates more precisely the purpose of this function.
}
\usage{
sienaAlgorithmCreate(fn, projname = "Siena", MaxDegree = NULL, Offset = NULL,
     useStdInits = FALSE, n3 = 1000, nsub = 4, n2start = NULL,
     dolby=TRUE, maxlike = FALSE, diagonalize=0.2*!maxlike,
     condvarno = 0, condname = "", firstg = 0.2, reduceg = 0.5,
     cond = NA, findiff = FALSE, seed = NULL,
     pridg=0.05, prcdg=0.05, prper=0.2, pripr=0.3, prdpr=0.3,
     prirms=0.05, prdrms=0.05, maximumPermutationLength=40,
     minimumPermutationLength=2, initialPermutationLength=20,
     modelType=NULL, behModelType=NULL, mult=5, simOnly=FALSE, localML=FALSE,
     truncation=5, doubleAveraging=0, standardizeVar=(diagonalize<1),
     lessMem=FALSE)

sienaModelCreate(fn, projname = "Siena", MaxDegree = NULL, Offset = NULL,
     useStdInits = FALSE, n3 = 1000, nsub = 4, n2start = NULL,
     dolby=TRUE, maxlike = FALSE, diagonalize=0.2*!maxlike,
     condvarno = 0, condname = "", firstg = 0.2, reduceg = 0.5,
     cond = NA, findiff = FALSE, seed = NULL,
     pridg=0.05, prcdg=0.05, prper=0.2, pripr=0.3, prdpr=0.3,
     prirms=0.05, prdrms=0.05, maximumPermutationLength=40,
     minimumPermutationLength=2, initialPermutationLength=20,
     modelType=NULL, behModelType=NULL, mult=5, simOnly=FALSE, localML=FALSE,
     truncation=5, doubleAveraging=0, standardizeVar=(diagonalize<1),
     lessMem=FALSE)
}
\arguments{
  \item{fn}{Function to do one simulation in the Robbins-Monro
    algorithm. Not to be touched.}
  \item{projname}{Character string name of project; the output file will be
    called projname.txt. No embedded spaces!!!\cr
	If \code{projname=NULL}, output will be written to a file in the temporary
	session directory, created as \code{\link{tempfile}(Siena)}.}
  \item{MaxDegree}{Named vector of maximum degree values for
    corresponding networks. Allows to restrict the model to networks
    with degrees not higher than this maximum.
    Names should be the names of all dependent network variables,
    in the same order as in the Siena data set.\cr
    Default as well as value 0 imply no restrictions.}
  \item{Offset}{Named vector of offset values for symmetric networks with
    \code{modelType = 3} (M.1), and for universal setting in Settings model.
    Names should be the names of all dependent network variables,
    in the same order as in the Siena data set.
    Default \code{NULL} implies values 0.}
  \item{useStdInits}{Boolean. If TRUE, the initial values in the effects
    object will be ignored and default values used instead.
    If FALSE, the initial values in the effects object will be used.}
  \item{n3}{Number of iterations in phase 3. For regular use with the
    Method of Moments, n3=1000 mostly suffices.
    For use in publications and for Maximum Likelihood,
    at least n3=3000 is advised. Sometimes much higher values are required
    for stable estimation of standard errors.}
  \item{nsub}{Number of subphases in phase 2.}
  \item{n2start}{Minimum number of interations in subphase 1 of phase 2;
    default is \code{2.52*(p+7)},
    where \code{p} = number of estimated parameters.}
  \item{dolby}{Boolean. Should there be noise reduction by regression
    on augmented data score. In most cases dolby=TRUE yields better
    convergence, but takes some extra computing time;
    if convergence is problematic, however, dolby=FALSE may be tried.
    Just use whatever works best.}
  \item{maxlike}{Whether to use maximum likelihood method or
    Method of Moments estimation.}
  \item{diagonalize}{Number between 0 and 1 (bounds included),
    values outside this interval will be truncated;
    for diagonalize=0 the complete estimated derivative matrix will be used
    for updates in the Robbins-Monro procedure;
    for diagonalize=1 only the diagonal entries will be used;
    for values between 0 and 1, the weighted average will be used
    with weight diagonalize for the diagonalized matrix.
    Has no effect for ML estimation.\cr
    Higher values are more stable, lower values potentially more efficient.
    Default: for ML estimation, diagonalize=0; for MoM estimation,
    diagonalize = 1.0.}
  \item{condvarno}{If \code{cond} (conditional simulation), the
    sequential number of the network
    or behavior variable on which to condition.}
  \item{condname}{If conditional, the name of the dependent variable on
    which to condition. Use one or other of \code{condname} or
    \code{condvarno} to specify the variable.}
  \item{firstg}{Initial value of scaling ("gain") parameter for updates
    in the Robbins-Monro procedure.}
  \item{reduceg}{Reduction factor for scaling ("gain") parameter for updates
    in the Robbins-Monro procedure (MoM only).}
  \item{cond}{Boolean. Only relevant for Method of Moments
    simulation/estimation.
    If TRUE, use conditional simulation; if FALSE, unconditional simulation.
    If missing, decision is deferred until \code{\link{siena07}},
    when it is set to TRUE if there is only one dependent variable,
    FALSE otherwise.}
  \item{findiff}{ Boolean: If TRUE, estimate derivatives using finite
    differences. If FALSE, use scores.}
  \item{seed}{Integer. Starting value of random seed. Not used if
    parallel testing.}
  \item{pridg}{Real number. Probability used in Metropolis-Hastings
    routine in ML estimation. See Siena_Algorithms.pdf.}
  \item{prcdg}{Real number. Probability used in Metropolis-Hastings
    routine in ML estimation. See Siena_Algorithms.pdf.}
  \item{prper}{Real number. Probability used in Metropolis-Hastings
    routine in ML estimation. See Siena_Algorithms.pdf.}
  \item{pripr}{Real number. Probability used in Metropolis-Hastings
    routine in ML estimation. See Siena_Algorithms.pdf.}
  \item{prdpr}{Real number. Probability used in Metropolis-Hastings
    routine in ML estimation. See Siena_Algorithms.pdf.}
  \item{prirms}{Real number. Probability used in Metropolis-Hastings
    routine in ML estimation. See Siena_Algorithms.pdf.}
  \item{prdrms}{Real number. Probability used in Metropolis-Hastings
    routine in ML estimation. See Siena_Algorithms.pdf.}
  \item{maximumPermutationLength}{Maximum length of permutation in steps in ML
    estimation.}
  \item{minimumPermutationLength}{Minimum length of permutation in steps in ML
    estimation.}
  \item{initialPermutationLength}{Initial length of permutation in steps in ML
    estimation.}
  \item{modelType}{Named vector indicating the type of model to be fitted for
    dependent network variables. Possible values are:\cr
    1=directed, 2:6 for symmetric networks only: 2=dictatorial forcing (D.1),
    3=Initiative model with reciprocal confirmation (M.1),
    4=Pairwise dictatorial forcing model (D.2),
    5=Pairwise mutual model (M.2), 6=Pairwise joint model (C.2).\cr
    Names should be the names of all dependent network variables,
    in the same order as in the Siena data set.\cr
    See Snijders and Pickup (2016) for the meanings of these models.\cr
    Default \code{NULL} implies 1 for directed or two-mode, 2 for symmetric.}
  \item{behModelType}{Named vector indicating the type of model to be fitted for
    behavioral dependent variables. Possible values are:\cr
    1=standard (restricted), 2=absorbing.\cr
    Names should be the names of all dependent behavioral variables,
    in the same order as in the Siena data set.\cr
    Default \code{NULL} implies values 1.}
  \item{mult}{Multiplication factor for maximum likelihood and Bayes. Number of
    steps per iteration is set to this multiple of the total distance
    between the observations at start and finish of the wave.
    Decreasing \code{mult} below a certain value has no further effect.\cr
    This can be either a number (which needs to be positive) or a vector
    of numbers, of length equal to the number of basic rate parameters in the
    model, i.e., the number of periods times the number of dependent
    variables.}
  \item{simOnly}{Logical: If TRUE, then the calculation of the covariance
    matrix and standard errors of the estimates at the end of
    Phase 3 of the estimation algorithm in function siena07 is skipped.
    This is suitable if nsub=0 and \code{\link{siena07}} is used only for the
    purpose of simulation.}
  \item{localML}{Logical: If TRUE, and \code{maxlike}, then calculations are
    sped up for models with all local effects.}
  \item{truncation}{Used for step truncation in the Robbins Monro algorithm
    (applied to deviate/(standard deviation)).}
  \item{doubleAveraging}{subphase after which double averaging is used
    in the Robbins Monro algorithm,
    which probably increases algorithm efficiency.}
  \item{standardizeVar}{Logical: whether to limit deviations used in
    Robbins-Monro updates to unit variances.}
  \item{lessMem}{Logical: whether to reduce storage during operation of
    \code{\link{siena07}}, and of the object produced, by leaving out arrays
    by iteration and by period of simulated statistics \code{sf2} and scores
    \code{ssc}.\cr
    if \code{lessMem=TRUE}, it will be impossible to run
    \code{\link{sienaTimeTest}} or \code{\link{sienaGOF}} on
    the object produced by \code{\link{siena07}}.}
}
\details{Model specification is done via this object for
  \code{\link{siena07}}.
  This function creates an object with the elements required to control the
  Robbins-Monro algorithm. Those not
  available as arguments can be changed manually where desired.\cr
  Further information about the implementation of the algorithm is in \cr
  http://www.stats.ox.ac.uk/~snijders/siena/Siena_algorithms.pdf.\cr
  Some of the examples use \code{projname=NULL}; this is just for the sake of
  checking the examples, not necessarily intended for normal use.
}
\value{Returns an object of class \code{sienaAlgorithm} containing
  values implied by the parameters.}
\references{ \itemize{
\item For the model types: \cr 
   Tom A. B. Snijders and Mark Pickup,
   Stochastic Actor-Oriented Models for Network Dynamics.
   In: Jennifer N. Victor, Mark Lubell and Alexander H. Montgomery,
   Oxford Handbook of Political Networks. Oxford University Press, 2016.
\item Also see \url{http://www.stats.ox.ac.uk/~snijders/siena/}
}}
\author{Ruth Ripley and Tom A.B. Snijders}
\seealso{\code{\link{siena07}}, \code{\link{simstats0c}}.}
\examples{
myAlgorithm <- sienaAlgorithmCreate(projname="NetworkDyn")
StdAlgorithm <- sienaAlgorithmCreate(projname="NetworkDyn", useStdInits=TRUE)
CondAlgorithm <- sienaAlgorithmCreate(projname="NetworkDyn", condvarno=1, cond=TRUE)
Max10Algorithm <- sienaAlgorithmCreate(projname="NetworkDyn", MaxDegree=c(mynet=10),
     modelType=c(mynet=1))
Beh2Algorithm <- sienaAlgorithmCreate(projname="NetBehDyn", behModelType=c(mybeh=2))
# where mynet is the name of the network object created by sienaDependent(),
# and mybeh the name of the behavior object created by the same function.
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{ classes }
