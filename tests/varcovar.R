## A changing covariate with a period that has only missing values
## (for example, a variable that was not measured at the first wave)
## must not give warnings, and its range for that period is NA.
library(RSiena)

alc <- s50a
alc[, 1] <- NA
mynet <- as_dependent_rsiena(array(c(s501, s502, s503), dim = c(50, 50, 3)))
alcohol <- varCovar(alc)

mydata <- withCallingHandlers(
	make_data_rsiena(mynet, alcohol),
	warning = function(w) stop("unexpected warning: ", conditionMessage(w)))

rangep <- attr(mydata$vCovars$alcohol, "rangep")
stopifnot(all(is.na(rangep[, 1])))
stopifnot(identical(as.numeric(rangep[, 2]), as.numeric(range(s50a[, 2]))))

## The descriptive report must still run.
invisible(capture.output(write_report(mydata, modelname = tempfile())))
