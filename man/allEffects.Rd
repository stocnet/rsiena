\name{allEffects}
\Rdversion{1.1}
\alias{allEffects}
\docType{data}
\title{
Internal data frame used to construct effect objects.
}
\description{
This data frame is used internally to construct effect objects.  
}
\usage{data(allEffects)}
\format{
  A data frame with values for the following 23 variables.
  \describe{
    \item{\code{effectGroup}}{a character vector}
    \item{\code{effectName}}{a character vector}
    \item{\code{functionName}}{a character vector}
    \item{\code{shortName}}{a character vector}
    \item{\code{endowment}}{a logical vector}
    \item{\code{interaction1}}{a character vector}
    \item{\code{interaction2}}{a character vector}
    \item{\code{type}}{a character vector}
    \item{\code{basicRate}}{a logical vector}
    \item{\code{include}}{a logical vector}
    \item{\code{randomEffects}}{a logical vector}
    \item{\code{fix}}{a logical vector}
    \item{\code{test}}{a logical vector}
    \item{\code{timeDummy}}{a character vector, default ","}
    \item{\code{initialValue}}{a numeric vector}
    \item{\code{parm}}{a numeric vector}
    \item{\code{functionType}}{a character vector}
    \item{\code{period}}{a character vector}
    \item{\code{rateType}}{a character vector}
    \item{\code{untrimmedValue}}{a numeric vector}
    \item{\code{effect1}}{a logical vector}
    \item{\code{effect2}}{a logical vector}
    \item{\code{effect3}}{a logical vector}
    \item{\code{interactionType}}{a character vector}
    \item{\code{local}}{a logical vector}
    \item{\code{setting}}{Settings name: '' (no settings), 'primary', 'universal' or the name of the defining covariate.}
  }
}
\details{
Used to define effects. Not for general user use.}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\keyword{datasets}
