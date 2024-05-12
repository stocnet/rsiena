\name{s50a}
\alias{s50a}
\docType{data}
\title{Alcohol use data: excerpt from "Teenage Friends and Lifestyle Study" data}
\description{
Alcohol use data from an excerpt of 50 girls from    
an excerpt of the network, alcohol consumption, and smoking data     
for 50 randomly chosen girls from the
Teenage Friends and Lifestyle Study data set. Useful as a small    
example of network and behaviour, for which models can be    
fitted quickly, and for which there are no missing values.
}
\format{
 A matrix of variables relating to the use of alcohol for the actors in
 the network. Three columns, one for each time point. Coding is 1--5, high
 values indicating higher consumption.
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
\seealso{\code{\link{s501}}, \code{\link{s502}}, \code{\link{s503}}, \code{\link{s50s}}
}
\examples{
mynet <- sienaDependent(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mybeh <- sienaDependent(s50a, type="behavior")
mydata <- sienaDataCreate(mynet, mybeh)
mydata
}
\keyword{datasets}