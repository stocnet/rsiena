#' @importFrom generics tidy
#' @export
generics::tidy
# tidy <- function(x) UseMethod("tidy") # just for testing, don't use because overwrites use in other packages

#' @export
tidy.sienaFit <- function(x, conf.int = FALSE, conf.level = 0.95, fromBayes = FALSE, ...) {
  
  terms <- 
  
  result <- tibble::tibble(term = x$requestedEffects$effectName,
                           estimate = x$theta,
                           std.error = x$se,
                           statistic = x$theta/x$se,
                           p.value = 2 * (1 - pnorm(abs(x$theta/x$se))))
  
  # if (conf.int) {
  #   ci <- confint(x, level = conf.level)
  #   result <- dplyr::left_join(result, ci, by = "term")
  # }
  
  result
}

#' #' @importFrom generics glance
#' #' @export
#' generics::glance
#' # glance <- function(x) UseMethod("glance") # just for testing, don't use because overwrites use in other packages
#' 
#' glance.sienaFit <- function(x, ...) {
#'   with(
#'     summary(x),
#'     tibble::tibble(
#'       # r.squared = r.squared,
#'       # adj.r.squared = adj.r.squared,
#'       # sigma = sigma,
#'       # statistic = fstatistic["value"],
#'       # p.value = pf(
#'       #   fstatistic["value"],
#'       #   fstatistic["numdf"],
#'       #   fstatistic["dendf"],
#'       #   lower.tail = FALSE
#'       # ),
#'       # df = fstatistic["numdf"],
#'       logLik = as.numeric(stats::logLik(x)),
#'       AIC = stats::AIC(x),
#'       BIC = stats::BIC(x),
#'       # deviance = stats::deviance(x),
#'       # df.residual = df.residual(x),
#'       nobs = x$nEvents
#'     )
#'   )
#' }
#' 
#' # to be implemented...
#' # #' @importFrom generics augment
#' # #' @export
#' # generics::augment