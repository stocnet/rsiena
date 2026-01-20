\name{RSiena-package}
\alias{RSiena-package}
\alias{RSiena}
\docType{package}
\title{
Simulation Investigation for Empirical Network Analysis
}
\description{
Fits statistical models to longitudinal sets of networks, and to   
longitudinal sets of networks and behavioral variables.   
Not only one-mode networks but also two-mode networks   
and multivariate networks are allowed.   
The models are stochastic actor-oriented models,   
described in Snijders (2017).

Recent versions of the package are distributed through GitHub,   
see \url{https://github.com/stocnet/rsiena/}.

Bug reports can be submitted at 
\url{https://github.com/stocnet/rsiena/issues}.
}
\details{
  The main flow of operations of this package is as follows.

  Data objects can be created from matrices and
  vectors using \code{\link{as_dependent_rsiena}}, \cr
  \code{\link{as_covariate_rsiena}},
  possibly \code{\link{as_composition_rsiena}}, \cr
  and finally \code{\link{make_data_rsiena}}.

  Effects are selected using an \code{\link{sienaEffects}} object,
  which can be created using \code{\link{make_specification}}
  and may be further specified by \code{\link{set_effect}},
  and \code{\link{set_interaction}}.

  Control of the model family, the estimation algorithm,
  and the output produced  requires control objects that
  define the settings of the algorithm,
  and which can be created by \code{\link{set_algorithm_saom}}, 
  \code{\link{set_model_saom}}, and \code{\link{set_output_saom}}.

  Function \code{\link{siena}} is used to fit a model.
  Function \code{\link{test_gof}} can be used for studying goodness of fit.
  
  Version 1.6 has established an overhaul of the function names used,
  while retaining the old function names for backward compatibility.

  A general introduction to the method is available in the tutorial paper
  Snijders, van de Bunt, and Steglich (2010).
  Next to the help pages, more detailed help is available in the manual
  (see below) and a lot of information is at the website (also see below).
  \tabular{ll}{
    Package: \tab RSiena\cr
    Type: \tab Package\cr
    Version: \tab 1.6.1\cr
    Date: \tab 2026-01-20\cr
    Depends: \tab R (>= 3.5.0)\cr
    Imports: \tab Matrix, lattice, parallel, MASS, methods, xtable\cr
    Suggests: \tab network, tools, codetools, tcltk\cr
    SystemRequirements: \tab GNU make\cr
    License: \tab GPL-2 | GPL-3 \cr
    LazyData: \tab yes\cr
    NeedsCompilation: \tab yes\cr
    BuildResaveData: \tab no\cr
  }
}
\author{
  Ruth Ripley, Krists Boitmanis, Tom Snijders, Felix Schoenenberger,
  Nynke Niezink, Christian Steglich, Viviana Amati.
  Contributions by Josh Lospinoso, Charlotte Greenan,
  Johan Koskinen, Mark Ortmann, Natalie Indlekofer, Mark Huisman,
  Christoph Stadtfeld, Per Block, Marion Hoffman, Michael Schweinberger,
  Robert Hellpap, Alvaro Uzaheta, Robert Krause, James Hollway, 
  Daniel Gotthardt, and Steffen Triebel.

  Maintainer: Christian Steglich <c.e.g.steglich@rug.nl>
}
\references{
    Amati, V., Schoenenberger, F., and Snijders, T.A.B. (2015),
    Estimation of stochastic actor-oriented models for the evolution of 
    networks by generalized method of moments.
    \emph{Journal de la Societe Francaise de Statistique} \bold{156}, 140--165.
    
    Amati, V., Schoenenberger, F., and Snijders, T.A.B. (2019),
    Contemporaneous statistics for estimation in stochastic actor-oriented
    co-evolution models. \emph{Psychometrika} \bold{84}, 1068--1096.
    
    Greenan, C. (2015),
    \emph{Evolving Social Network Analysis: developments in statistical
    methodology for dynamic stochastic actor-oriented models}.
    DPhil dissertation, University of Oxford.
    
    Niezink, N.M.D., and Snijders, T.A.B. (2017),
    Co-evolution of Social Networks and Continuous Actor Attributes.
    \emph{The Annals of Applied Statistics} \bold{11}, 1948--1973.
    
    Schweinberger, M., and Snijders, T.A.B. (2007),
    Markov models for digraph panel data:
    Monte Carlo based derivative estimation.
    \emph{Computational Statistics and Data Analysis} \bold{51}, 4465--4483.
    
     Snijders, T.A.B. (2001),
    The statistical evaluation of social network dynamics.
    \emph{Sociological Methodology} \bold{31}, 361--395.
    
    Snijders, T.A.B. (2017), Stochastic Actor-Oriented Models for Network
    Dynamics. \emph{Annual Review of Statistics and Its Application}
    \bold{4}, 343--363.
    
    Snijders, T.A.B., Koskinen, J., and Schweinberger, M. (2010).
    Maximum likelihood estimation for social network dynamics.
    \emph{Annals of Applied Statistics} \bold{4}, 567--588.
    
    Snijders, T.A.B., Steglich, C.E.G., and Schweinberger, Michael (2007), 
    Modeling the co-evolution of networks and behavior.
    Pp. 41--71 in \emph{Longitudinal models in the behavioral
    and related sciences},
    edited by van Montfort, K., Oud, H., and Satorra, A.; Lawrence Erlbaum.
    
    Steglich, C.E.G., Snijders, T.A.B., and Pearson, M.A. (2010),
    Dynamic networks and behavior:     Separating selection from influence.
    \emph{Sociological Methodology} \bold{40}, 329--393.
    Information about the implementation of the algorithm is in
     \url{https://www.stats.ox.ac.uk/~snijders/siena/Siena_algorithms.pdf}.
     
    Further see \url{https://www.stats.ox.ac.uk/~snijders/siena/} and \cr
    \url{https://github.com/stocnet/rsiena/wiki}.
}
\seealso{\code{\link{siena}}}

\examples{
mynet <- as_dependent_rsiena(array(c(s501, s502), dim=c(50, 50, 2)))
mybeh  <- as_dependent_rsiena(s50a[,1:2], type="behavior")
mydata <- make_data_rsiena(mynet, mybeh)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, avAlt, depvar="mybeh", covar1="mynet")
myalgo <- set_algorithm_saom(nsub=2, n3=100, seed=1291)
# nsub=2, n3=100 is used here for having a brief computation, not for practice.
(ans <- siena(mydata, effects=myeff, control_algo=myalgo,
                                            silent=TRUE, batch=TRUE))
summary(ans)
}
\keyword{package}
