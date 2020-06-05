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
The models are stochastic actor-oriented models.

Package \code{"RSienaTest"} is the development version, and
is distributed through R-Forge, see
\url{http://r-forge.r-project.org/R/?group_id=461}.
Package \code{"RSiena"} is the official release.
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
    Version: \tab 1.2-25\cr
    Date: \tab 2020-04-06\cr
    Depends: \tab R (>= 3.0.0)\cr
    Imports: \tab Matrix\cr
    Suggests: \tab tcltk, network, codetools, lattice, MASS, parallel,
    xtable, tools, utils\cr
    SystemRequirements: \tab GNU make, tcl/tk 8.5, Tktable\cr
    License: \tab GPL-2 \cr
    LazyData: \tab yes\cr
    NeedsCompilation: \tab yes\cr
    BuildResaveData: \tab no\cr
  }
}
\author{
  Ruth Ripley, Krists Boitmanis, Tom Snijders, Felix Schoenenberger,
  Nynke Niezink.
  Contributions by Josh Lospinoso, Charlotte Greenan, Christian Steglich,
  Johan Koskinen, Mark Ortmann, Natalie Indlekofer,
  Christoph Stadtfeld, Per Block, Marion Hoffman, and Robert Hellpap.

  Maintainer: Tom A.B. Snijders <tom.snijders@nuffield.ox.ac.uk>
}
\references{
  \itemize{
    \item Schweinberger, Michael, and Snijders, Tom A.B. (2007).
      Markov models for digraph panel data:
      Monte Carlo-based derivative estimation.
      \emph{Computational Statistics and Data Analysis} \bold{51}, 4465--4483.
    \item Snijders, Tom A.B. (2001).
      The statistical evaluation of social network dynamics.
      \emph{Sociological Methodology} \bold{31}, 361-395.
    \item Snijders, Tom A.B. (2017). Stochastic Actor-Oriented Models for Network
      Dynamics. 
      \emph{Annual Review of Statistics and Its Application} \bold{4}, 343--363.
    \item  Snijders, Tom A.B., van de Bunt, Gerhard G., and
      Steglich, Christian E.G. (2010).
      Introduction to actor-based models for network dynamics.
      \emph{Social Networks} \bold{32}, 44--60.
    \item Snijders, Tom A.B., Steglich, Christian E.G., and Schweinberger,
      Michael (2007). Modeling the co-evolution of networks and behavior.
      Pp. 41--71 in \emph{Longitudinal models in the behavioral
        and related sciences},
      edited by Kees van Montfort, Han Oud and Albert Satorra; Lawrence Erlbaum.
    \item Steglich, Christian E.G., Snijders, Tom A.B., and
      Pearson, Michael A. (2010). Dynamic networks and behavior:
      Separating selection from influence. \emph{Sociological Methodology}
      \bold{40}, 329--393.
    \item The manual: \url{http://www.stats.ox.ac.uk/~snijders/siena/RSiena_Manual.pdf}
    \item The website: \url{http://www.stats.ox.ac.uk/~snijders/siena/}.
  }
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
