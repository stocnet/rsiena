\name{set_algorithm}
\alias{set_algorithm_saom}
\alias{set_model_saom}
\alias{set_output_saom}
\alias{set_algorithm}
\alias{set_model}
\alias{set_output}
\alias{sienaModelCreate}
\alias{sienaAlgorithmCreate}
\alias{sienaModel}
\alias{sienaAlgorithm}
\alias{sienaModel}
\alias{sienaAlgorithmSettings}
\alias{sienaOutputOptions}
\alias{model.create}
\alias{modelType}
\alias{behModelType}

\title{Functions to create objects containing various aspects of
  algorithm specifications for parameter estimation in RSiena}
\description{
  Create objects containing various aspects of
  algorithm specifications for parameter estimation in RSiena.  \cr
  \code{set_model_saom} sets the specifications of the Stochastic
  Actor-oriented Model. \cr
  \code{set_algorithm_saom} sets the specifications of the algorithm
  for parameter estimation of the model. \cr
  \code{set_output_saom} sets the specifications of the generated output.

  These three aspects are combined in the older functions
  \code{sienaAlgorithmCreate} and \cr
  \code{sienaModelCreate},
  which are identical functions although having different names.
  The older functions are not preferred any more,
  but are kept for backward compatibility.
}
\usage{
  set_model_saom(modelType = NULL, behModelType = NULL,
                        MaxDegree = NULL, Offset = NULL)

  set_algorithm_saom(maxlike = FALSE, gmm = FALSE,
    cond = NA, condvarno = 0, condname = "",
    simOnly = FALSE,
    targets = NULL, thetaValues = NULL, seed = NULL,
    n3 = 1000, nsub = 4, n2start = NULL, firstg = 0.2, reduceg = 0.5,
    truncation = 5, doubleAveraging = 0,
    diagonalize = 0.2*!maxlike, standardizeVar = (diagonalize<1),
    dolby = TRUE, useStdInits = FALSE, findiff = FALSE,
    mult = 5, prML = 1, maximumPermutationLength = 40,
    minimumPermutationLength = 2, initialPermutationLength = 20,
    localML = FALSE)

set_output_saom(outputName = NULL, lessMem = FALSE,
    returnThetas = FALSE, 
    returnChains = FALSE, returnDataFrame = FALSE,
    returnChangeContributions = FALSE)

sienaAlgorithmCreate(fn, projname = "Siena", MaxDegree = NULL,
     Offset = NULL, useStdInits = FALSE, 
     n3 = 1000, nsub = 4, n2start = NULL,
     dolby = TRUE, maxlike = FALSE, gmm = FALSE, 
     diagonalize = 0.2*!maxlike,
     condvarno = 0, condname = "", firstg = 0.2, reduceg = 0.5,
     cond = NA, findiff = FALSE, seed = NULL,
     prML = 1, maximumPermutationLength = 40,
     minimumPermutationLength = 2, initialPermutationLength = 20,
     modelType = NULL, behModelType = NULL, mult = 5, simOnly = FALSE,
     localML = FALSE, truncation = 5, doubleAveraging = 0,
     standardizeVar = (diagonalize<1), lessMem = FALSE, silent = FALSE)

sienaModelCreate(fn, projname = "Siena", MaxDegree = NULL,
     Offset = NULL, useStdInits = FALSE, 
     n3 = 1000, nsub = 4, n2start = NULL,
     dolby = TRUE, maxlike = FALSE, gmm = FALSE, 
     diagonalize = 0.2*!maxlike,
     condvarno = 0, condname = "", firstg = 0.2, reduceg = 0.5,
     cond = NA, findiff = FALSE, seed = NULL,
     prML = 1, maximumPermutationLength = 40,
     minimumPermutationLength = 2, initialPermutationLength = 20,
     modelType = NULL, behModelType = NULL, mult = 5, simOnly = FALSE,
     localML = FALSE, truncation = 5, doubleAveraging = 0,
     standardizeVar = (diagonalize<1), lessMem = FALSE, silent = FALSE)
}
\arguments{
  For \code{set_model_saom}:
  \item{modelType}{Named vector indicating the type of model to be fitted for
    dependent network variables.
    (See the examples below for how to specify a named vector.)\cr
    Possible values are:\cr
    1 = directed standard,
    \cr
    2:6 for symmetric networks only: 2 = dictatorial forcing (D.1),
    3 = Initiative model with reciprocal confirmation (M.1),
    4 = Pairwise dictatorial forcing model (D.2),
    5 = Pairwise mutual model (M.2), 6 = Pairwise joint model (C.2),\cr
    7:10 for directed one-mode only:
    7 = Double Step model with double step probability 0.25,
    8 = Double Step model with double step probability 0.50,
    9 = Double Step model with double step probability 0.75,
    10 = Double Step model with double step probability 1.00,\cr
    11 = Contemporaneous evaluation statistics model.
    \cr
    Names should be the names of all dependent network variables,
    in the same order as in the Siena data set.\cr
    See Snijders and Pickup (2016) for the meanings of the various models
    for symmetric networks.\cr
    Default \code{NULL} implies 1 for directed or two-mode, 2 for symmetric.}
  \item{behModelType}{Named vector indicating the type of model to be fitted for
    behavioral dependent variables.
    (See the examples below for how to specify a named vector.)\cr
    Possible values are:\cr
    1 = standard (restricted), 2 = absorbing.\cr
    Names should be the names of all dependent behavioral variables,
    in the same order as in the Siena data set.\cr
    Default \code{NULL} implies values 1.}
  \item{MaxDegree}{Named vector of maximum degree values for
    corresponding networks. Allows to restrict the model to networks
    with degrees not higher than this maximum.
    Names should be the names of all dependent network variables,
    in the same order as in the Siena data set.\cr
    Default as well as value 0 imply no restrictions.\cr
    This option is not available for maximum likelihood estimation.}
  \item{Offset}{Named vector of offset values for symmetric networks with
    \code{modelType = 3} (M.1), and for universal setting in Settings model.
    Names should be the names of all dependent network variables,
    in the same order as in the Siena data set.
    Default \code{NULL} implies values 0.}
  For \code{set_algorithm_saom}:
  \item{maxlike}{Whether to use maximum likelihood method or
    Method of Moments estimation.}
  \item{gmm}{Whether to use the Generalized Method of Moments or the regular
    Method of Moments estimation.}
  \item{cond}{Boolean. Only relevant for Method of Moments
    simulation/estimation.
    If TRUE, use conditional simulation; if FALSE, unconditional simulation.
    If missing, decision is deferred until \code{\link{siena}},
    when it is set to TRUE if there is only one dependent variable,
    FALSE otherwise.}
  \item{condvarno}{If \code{cond} (conditional simulation), the
    sequential number of the network
    or behavior variable on which to condition.}
  \item{condname}{If conditional, the name of the dependent variable on
    which to condition. Use one or other of \code{condname} or
    \code{condvarno} to specify the variable.}
  \item{simOnly}{Logical: If TRUE, then the calculation of the covariance
    matrix and standard errors of the estimates at the end of
    Phase 3 of the estimation algorithm is skipped.
    This is suitable if nsub = 0 and \code{\link{siena}}
    is used only for the purpose of simulation.}
  \item{targets}{Numeric vector of length equal to the number of estimated
    parameters, meant to supersede the targets calculated from the data set;
    see "Details". Not for regular use.}
  \item{thetaValues}{If not \code{NULL}, this should be a matrix with parameter
    values to be used in Phase 3. The number of columns must be equal to the
    number of estimated parameters in the effects object (if conditional
    estimation is used, without the rate parameters for the conditioning
    dependent variable). Can only be used if \code{simOnly=TRUE}.}
  \item{seed}{Integer. Starting value of random seed.}
  \item{n3}{Number of iterations in phase 3. For regular use with the
    Method of Moments, n3 = 1000 mostly suffices.
    For use in publications and for Maximum Likelihood,
    at least n3 = 3000 is advised. Sometimes much higher values are required
    for stable estimation of standard errors.}
  \item{nsub}{Number of subphases in phase 2.}
  \item{n2start}{Minimum number of iterations in subphase 1 of phase 2;
    default is \code{2.52*(p+7)}, where \code{p} = number of estimated
    parameters; if \code{useCluster = TRUE} in the call of
    \code{\link{siena07}}, this is divided by \code{nbrNodes}.  }
  \item{firstg}{Initial value of scaling ("gain") parameter for updates
    in the Robbins-Monro procedure.}
  \item{reduceg}{Reduction factor for scaling ("gain") parameter for updates
    in the Robbins-Monro procedure (MoM only).}
  \item{truncation}{Used for step truncation in the Robbins Monro algorithm
    (applied to deviate/(standard deviation)).}
  \item{doubleAveraging}{subphase after which double averaging is used
    in the Robbins Monro algorithm,
    which probably increases algorithm efficiency.}
  \item{diagonalize}{Number between 0 and 1 (bounds included),
    values outside this interval will be truncated;
    for diagonalize = 0 the complete estimated derivative matrix will be used
    for updates in the Robbins-Monro procedure;
    for diagonalize = 1 only the diagonal entries will be used;
    for values between 0 and 1, the weighted average will be used
    with weight diagonalize for the diagonalized matrix.
    Has no effect for ML estimation.\cr
    Higher values are more stable, lower values potentially more efficient.
    Default: for ML estimation, diagonalize = 0; for MoM estimation,
    diagonalize = 0.2.}
  \item{standardizeVar}{Logical: whether to limit deviations used in
    Robbins-Monro updates to unit variances.}
  \item{dolby}{Boolean. Should there be noise reduction by regression
    on augmented data score. In most cases dolby = TRUE yields better
    convergence, but takes some extra computing time;
    if convergence is problematic, however, dolby = FALSE may be tried.
    Just use whatever works best.}
  \item{useStdInits}{Boolean. If TRUE, the initial values in the effects
    object will be ignored and default values used instead.
    If FALSE, the initial values in the effects object will be used.}
  \item{findiff}{ Boolean: If TRUE, estimate derivatives using finite
    differences. If FALSE, use scores.}
  \item{mult}{Multiplication factor for maximum likelihood and Bayes. Number of
    steps per iteration is set to this multiple of the total distance
    between the observations at start and finish of the wave (and rounded).
    Decreasing \code{mult} below a certain value has no further effect.\cr
     \code{mult} can be either a number (which needs to be positive) or a vector
    of numbers, of length equal to the total number of periods. Note that for
    multi-group data, the total number of periods is equal to the number
    of groups times the number of periods per group
    (if the latter is constant).}
  Algorithm parameters only relevant for maximum likelihood
    and Bayesian estimation:
  \item{prML}{Either one real number, or a vector of 7 numbers.
    Determines update probabilities used in Metropolis-Hastings
    routine in ML estimation. Should be nonnegative; if a vector,
    the sum should be < =  1. See \code{Details}.}
  \item{maximumPermutationLength}{Maximum length of permutation
    in steps in ML estimation.}
  \item{minimumPermutationLength}{Minimum length of permutation
    in steps in MLestimation.}
  \item{initialPermutationLength}{Initial length of permutation
    in steps in MLestimation.}
  \item{localML}{Logical: If TRUE, and \code{maxlike}, then calculations
    are sped up for models with all local effects.}
  For \code{set_output_saom}:
  \item{outputName}{If \code{outputName = NULL} and the Siena data set
   has name \code{dataName}, output will be written  by
   \code{\link{siena}} to a file with name \code{dataName_out.txt}
   in the working directory, appending if this file already exists.\cr
    If \code{outputName} is not \code{NULL}, it should be a character string
    without embedded spaces, and output will be written
    to a file with name \code{outputName_out.txt}.}
  \item{lessMem}{Logical: whether to reduce storage during operation of
    \code{\link{siena}}, and of the object produced, by leaving out
    arrays by iteration and by period of simulated statistics \code{sf2}
    and scores \code{ssc}.\cr
    if \code{lessMem = TRUE}, it will be impossible to run
    \code{\link{test_time}} or \code{\link{test_gof}} on
    the object produced by \code{\link{siena}}.}
  \item{returnThetas}{Boolean: whether to return theta values
   and generated estimation statistics of Phase 2 runs.}
  \item{returnChains}{Boolean. Whether to return the chains generated in
    the Phase 3 runs.}
  \item{returnDataFrame}{Boolean. Whether to return the chains as
    lists or data frames.}
  \item{returnChangeContributions}{Boolean. Whether to return the change contributions.
     See \code{\link{interpret_size_dynamics}}.}

  And additional, for \code{sienaAlgorithmeCreate}:
  \item{fn}{Function to do one simulation in the Robbins-Monro
    algorithm. Not to be touched.}
  \item{projname}{Character string name of project; the output file will be
    called projname.txt. No embedded spaces!!!\cr
    If \code{projname = NULL}, output will be written to a file in the temporary
    session directory, created as \code{\link{tempfile}(Siena)}.}
  \item{silent}{Logical: whether to give a note about the output file.}
}
\details{Model specification is done via this object for
  \code{\link{siena}}, which calls function \code{\link{simstats0c}}.
  This function creates an object with the elements required to control the
  Robbins-Monro algorithm. Those not
  available as arguments can be changed manually when desired.\cr
  The value \code{prML = 1} defines the defaults valid in RSiena up to version
  1.3.16. For ML estimation with only one dependent variable,
  \code{prML = 2} may be more efficient.\cr
  If \code{prML} is given as a vector of 7 probabilities, these are,
  consecutively: the probabilities of inserting a diagonal step, deleting a
  diagonal step, permuting, inserting a CCP, deleting a CCP, inserting random
  missing, deleting random missing; the residual (1 minus the sum)
  is the probability of a move step.\cr
  Further information about the implementation of the algorithm is in \cr
  \url{https://www.stats.ox.ac.uk/~snijders/siena/Siena_algorithms.pdf}.\cr
  Some of the examples use \code{projname = NULL}; this is just for the sake of
  checking the examples, not necessarily intended for normal use.
}
\value{\code{set_model_saom} returns an object of class \code{sienaModel},
  \code{set_algorithm_saom} returns an object
  of class \code{sienaAlgorithmSettings},
  \code{set_output_saom} returns an object
  of class \code{sienaOutputOptions},
  and \code{sienaAlgorithmCreate} returns an object of class \cr
  \code{sienaAlgorithm}, all being lists containing values implied by the parameters.
}
\references{
    For \code{modelType}:\cr
   Snijders, T.A.B., and Pickup, M. (2016),
   Stochastic Actor-Oriented Models for Network Dynamics.
   In: Victor, J.N., Lubell, M., and Montgomery, A.H.,
   \emph{Oxford Handbook of Political Networks}. Oxford University Press.

    Also see \url{https://www.stats.ox.ac.uk/~snijders/siena/}
}
\author{Ruth Ripley and Tom A.B. Snijders}
\seealso{\code{\link{siena}}, \code{\link{simstats0c}}.}
\examples{
CondAlgorithm <- set_algorithm_saom(cond = TRUE, condvarno = 1)
Max10Model <- set_model_saom(MaxDegree = c(mynet = 10),
                                modelType = c(mynet = 1))
Beh2Model <- set_model_saom(behModelType = c(mybeh = 2))
# where mynet is the name of the network object
# created by as_dependent_rsiena,
# and mybeh the name of the behavior object created by the same function.
Algorithm2 <- set_algorithm_saom(nsub = 2, n2start = 2000)
OutAlgorithm <- set_output_saom(outputName = "LaterEstimation",
                        returnThetas = TRUE)
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{classes}
