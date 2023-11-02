\name{s50s}
\alias{s50s}
\docType{data}
\title{Smoking data: excerpt from "Teenage Friends and Lifestyle Study" data}
\description{
Smoking data from    
an excerpt of the network, alcohol consumption, and smoking data     
for 50 randomly chosen girls from the
Teenage Friends and Lifestyle Study data set. Useful as a small    
example of network and behaviour, for which models can be    
fitted quickly, and for which there are no missing values.
}
\format{
 A matrix of variables relating to the smoking habits for the actors in
 the network. Three columns, one for each time point. Coding is 1--3:
 1 = no smoking, 2 = moderate smoking, 3 = serious smoking.
}
\source{
\url{https://www.stats.ox.ac.uk/~snijders/siena/s50_data.zip}
}
\references{
West, P. and Sweeting, H. (1995), \emph{Background Rationale and Design of the
West of Scotland 11-16 Study}.
Working Paper No. 52. MRC Medical Sociology Unit Glasgow.

See \url{https://www.stats.ox.ac.uk/~snijders/siena/s50_data.htm}
}
\seealso{\code{\link{s501}}, \code{\link{s502}}, \code{\link{s503}}, \code{\link{s50a}}
}
\examples{
mynet <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
myvar <- varCovar(s50s)
mydata <- sienaDataCreate(mynet, myvar)
mydata
}
\keyword{datasets}