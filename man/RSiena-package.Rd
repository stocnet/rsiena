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
see \url{https://github.com/snlab-nl/rsiena/ }.
}
\details{
  The main flow of operations of this package is as follows.

  Data objects can be created from matrices and
  vectors using \code{\link{sienaDependent}}, \code{\link{coCovar}},
  \code{\link{varCovar}}, \code{\link{coDyadCovar}}, etc.,
  and finally \code{\link{sienaDataCreate}}.

  Effects are selected using an \code{\link{sienaEffects}} object,
  which can be created using \code{\link{getEffects}}
  and may be further specified by \code{\link{includeEffects}},
  \code{\link{setEffect}}, and \code{\link{includeInteraction}}.

  Control of the estimation algorithm requires a
  \code{\link{sienaAlgorithm}} object that
  defines the settings (parameters) of the algorithm,
  and which can be created by \code{\link{sienaAlgorithmCreate}}.

  Function \code{\link{siena07}} is used to fit a model.
  Function \code{\link{sienaGOF}} can be used for studying goodness of fit.

  A general introduction to the method is available in the tutorial paper
  Snijders, van de Bunt, and Steglich (2010).
  Next to the help pages, more detailed help is available in the manual
  (see below) and a lot of information is at the website (also see below).
  \tabular{ll}{
    Package: \tab RSiena\cr
    Type: \tab Package\cr
    Version: \tab 1.4.2\cr
    Date: \tab 2023-12-14\cr
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
  Robert Hellpap, Alvaro Uzaheta, Robert Krause, James Hollway, and Steffen Triebel.

  Maintainer: Tom A.B. Snijders <tom.snijders@nuffield.ox.ac.uk>
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
    Further see \url{https://www.stats.ox.ac.uk/~snijders/siena/} .
}
\seealso{\code{\link{siena07}}}

\examples{
mynet1 <- sienaDependent(array(c(tmp3, tmp4), dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff <- getEffects(mydata)
myeff <- includeEffects(myeff, transTrip)
myeff
myalgorithm <- sienaAlgorithmCreate(nsub=3, n3=200)
ans <- siena07(myalgorithm, data=mydata, effects=myeff, batch=TRUE)
summary(ans)
}
\keyword{ package }
