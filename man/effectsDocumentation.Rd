\name{effectsDocumentation}
\alias{effectsDocumentation}
\title{
Function to create a table of documentation of effect names, short names etc.
}
\description{
  Produces a table of the shortnames and other information for effects,  
  either in html or latex.
}
\usage{
effectsDocumentation(effects = NULL, type = "html", display = (type=="html"),
     filename = ifelse(is.null(effects), "effects", deparse(substitute(effects))))
}
%- maybe also "usage" for other objects documented here.
\arguments{
  \item{effects}{
    A Siena effects object, or NULL.
}
  \item{type}{
    Type of output required. Valid options are "html" or "latex".
}
\item{display}{
  Boolean: should the output be displayed after creation. Only
  applicable to html output.
}
\item{filename}{
  Character string denoting file name.
}
}
\details{
  If \code{effects=NULL}, the \code{allEffects} object is written to
  a table, either latex or html. This table presents all the available
  effects present in this version of RSiena, not delimited
  by a particular data set. The default file name is "effects.tex" or
  "effects.html", respectively.

  The table lists all effects, with their name,
  shortName, whether an endowment (and creation) effect exists,
  the value of an effect parameter - if any -,
  and the interactionType (which can be empty or: "ego" or "dyadic" for
  dependent network variables; "OK" for dependent behavior variables).
  The latter is important for knowing how the effects can be used
  in interaction effects. (See \code{\link{includeInteraction}}).

  If an existing effects object is specified for \code{effects}, then all
  available effects in this effects object are listed.
  This table lists the name (i.e., dependent variable), effect name,
  shortName, type (rate/evaluation/endowment/creation),
  the variables defined as \code{interaction1} and \code{interaction2}
  (see \code{\link{includeEffects}}) that specify this effect,
  the value of an effect parameter - if any -, and the interactionType.\cr
  The GMoM effects, which are those with \code{type=gmm}, are listed at the end.
  For these, the distinction between the fields \code{name} and
  \code{interaction1}, referring to the dependent and the explanatory roles
  of the variables, has no meaning.\cr
  The default root file name is the name of the input effects object.
}
\value{
  Nothing returned. Output files are created in the current working directory.
}
\seealso{\code{\link{getEffects}}, \code{\link{includeEffects}},
    \code{\link{summary.sienaEffects}}, \code{\link{includeInteraction}}.}
\references{See \url{https://www.stats.ox.ac.uk/~snijders/siena/}}
\author{Ruth Ripley, Tom A.B. Snijders}
\examples{
\dontrun{effectsDocumentation()}
}
% Add one or more standard keywords, see file "KEYWORDS" in the
% R documentation directory.
\keyword{ print }
