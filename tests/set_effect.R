## set_effect with a list of effects must respect include = FALSE.
library(RSiena)

mynet <- as_dependent_rsiena(array(c(s501, s502, s503), dim = c(50, 50, 3)))
alcohol <- as_covariate_rsiena(s50a)
mydata <- make_data_rsiena(mynet, alcohol)
myeff <- make_specification(mydata)
myeff <- set_effect(myeff, list(egoX, altX, egoXaltX), covar1 = "alcohol",
                    verbose = FALSE)
included <- function(eff, sn) {
	eff$include[eff$shortName %in% sn & eff$type == "eval" &
		eff$interaction1 == "alcohol"]
}
stopifnot(all(included(myeff, c("egoX", "altX", "egoXaltX"))))

myeff2 <- set_effect(myeff, list(egoXaltX, altX), covar1 = "alcohol",
                     include = FALSE, verbose = FALSE)
stopifnot(!any(included(myeff2, c("egoXaltX", "altX"))))
stopifnot(all(included(myeff2, "egoX")))
