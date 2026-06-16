\name{transformScript}
\alias{transformScript}
\title{Function to transform \code{RSiena}  and \code{multiSiena} scripts }
\description{
   Function to transform a script with  'old' function names
     of \code{RSiena}  and \code{multiSiena} to the new names
     that are used since version \code{RSiena.1.6.0} and
     \code{multiSiena.1.2.36}.
}
\usage{
transformScript(theoldscript, fileName="newScript.R",
                    linelength=80, initialblanks=6, logfile="log.txt",
                    keepOldAlgorithms=FALSE)
}
\arguments{
  \item{theoldscript}{Character vector: an R script which uses \code{RSiena}
   and/or \code{multiSiena}.}
  \item{fileName}{Character string naming a file for the new script.}
  \item{linelength}{The maximum line length for the new script.}
  \item{initialblanks}{The number of blanks at the start of each follow-up
     line for a given transformed R command.}
  \item{logfile}{Character string naming a file, for the log file.}
  \item{keepOldAlgorithms}{Boolean: should the old calls of
        \code{sienaAlgorithmCreate} be retained in the new script,
        even though they are redundant.}
}
\details{This function transforms \code{theoldscript}, supposedly a script with
     'old' function names of \code{RSiena} and \code{multiSiena},
     and transforms it to \code{fileName} in which the old
     names are replaced with the new names
     that are used since version \code{RSiena.1.6.0}; the use of the functions
     is adapted correspondingly.\cr
     The calls for constructing the algorithms required for
     the  estimation functions \code{siena} and \code{multi_siena}
     are positioned immediately before the calls
     of the respective estimation functions. This may involve duplication.\cr
     For creation of a time-constant dyadic covariate ("\code{coDyadCovar}")
     the argument \code{type=oneMode} is always used;
     the user should check whether, instead,
     it should be \code{type=bipartite}.\cr
     This function is not guaranteed to be correct,
     and the result should be checked!\cr
     Hash signs (\code{#}) in strings that are names of objects
     in \code{RSiena} functions will lead to errors, because it is assumed
     that everything following a has sign in a line 
     is not part of the command.\cr
     Functions utilizing other \code{RSiena} functions 
     may also lead to errors.\cr
     If an error is encountered, the end of the log file will show the line 
     number to which the transformation proceeded without errors.\cr
     The use of \code{keepOldAlgorithms} is just for the sake of
     possible information or checking. These calls of
     \code{sienaAlgorithmCreate} can be dropped from the new script
     without any consequences for its functioning.\cr\cr
     Note that using \code{""} for a filename will write to the console.\cr\cr
     This function transforms calls of \code{RSiena} and
     \code{multiSiena} functions to calls utilizing the new names.
     The only \code{RSiena} functions in \code{theoldscript}
     that are executed are calls of \code{\link{sienaAlgorithmCreate}};
     their results are used for transforming
     calls of \code{\link{siena07}} and \code{sienaBayes}.
}
\value{Invisibly, the new script will be returned as a character vector.
}
\references{See the Chapter on "New function names" in the manual
   available from \url{https://www.stats.ox.ac.uk/~snijders/siena/ }.}
\author{Tom Snijders}
\examples{
oldscript <- c("mynet <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))",
    "# construct actor covariate",
"drink <- coCovar(s50a[1,], centered=FALSE)",
"mydata <- sienaDataCreate(mynet, drink)",
"mymodel <- getEffects(mydata)",
"mymodel <- setEffect(mymodel, transTrip)",
"mymodel <- includeEffects(mymodel, simX ,egoX, altX, interaction1='drink')",
"mycontrols0 <- sienaAlgorithmCreate(projname='ans.txt')",
"(ans <- siena07(mycontrols0, data=mydata, effects=mymodel))",
"mycontrols <- sienaAlgorithmCreate(projname='ans.txt', nsub=1, ",
"                                   n2start=2000, n3=2000, seed=1234)",
"(ans1 <- siena07(mycontrols, data=mydata, effects=mymodel, prevAns=ans))",
"(mtest <- Multipar.RSiena(ans1, tested=3:4))")
transformScript(oldscript, fileName="", logfile=tempfile())
}
\keyword{classes}
