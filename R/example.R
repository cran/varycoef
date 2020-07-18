#' Sample Function for SVCs
#'
#' @description Samples SVCs on a regular quadratic (Cartesian) grid. The SVCs
#' have all mean 0 and an exponential covariance function is used.
#'
#' @param m   (\code{numeric(1)}) \cr
#'    Number of observations in one dimension, i.i, the square root number of
#'    total number of observation locations \eqn{n = m^2}.
#' @param p   (\code{numeric(1)}) \cr
#'    Number of SVCs.
#' @param cov_pars (\code{data.frame(p, 2)}) \cr
#'    Contains the covariance parameters of SVCs. The two columns must have the
#'    names \code{"var"} and \code{"scale"}. These covariance parameters are
#'    then used for sampling the respective SVCs.
#' @param nugget  (\code{numeric(1)}) \cr
#'    Variance of the nugget / error term.
#' @param seed  (\code{numeric(1)}) \cr
#'    Seed set within the function for sampling.
#' @param given.locs (\code{NULL} or \code{data.frame(n, 2)}) \cr
#'    If \code{NULL}, the observations locations are sampled from a regular grid,
#'    Otherwise, the \code{data.frame} contains the observation locations.
#'    The data frame must have two columns of name \code{"x"} and \code{"y"}.
#'    The number of observations is then the number of rows \code{n}.
#'
#' @return \code{SpatialPointsDataFrame} \cr
#'    (see \code{\link[sp]{SpatialPointsDataFrame-class}}) of the sampled SVC
#'    including the nugget.
#'
#' @examples
#' # number of SVC
#' p <- 3
#' # sqrt of total number of observations
#' m <- 20
#' # covariance parameters
#' (pars <- data.frame(var = c(0.1, 0.2, 0.3),
#'                     scale = c(0.3, 0.1, 0.2)))
#' nugget.var <- 0.05
#'
#' # function to sample SVCs
#' sp.SVC <- fullSVC_reggrid(m = m, p = p,
#'                           cov_pars = pars,
#'                           nugget = nugget.var)
#'
#' library(sp)
#' # visualization of sampled SVC
#' spplot(sp.SVC, colorkey = TRUE)
#'
#' @importFrom RandomFields RMnugget RFsimulate RMexp
#' @importFrom sp SpatialPointsDataFrame
#' @export
fullSVC_reggrid <- function(m, p, cov_pars, nugget, seed = 123, given.locs = NULL) {

  if (is.null(given.locs)) {
    # number of observations
    n <- as.integer(m)^2

    # regular grid locations
    locs <- expand.grid(x = seq(0, 1, length.out = as.integer(m)),
                        y = seq(0, 1, length.out = as.integer(m)))
  } else {
    # take given locations
    locs <- given.locs
  }


  set.seed(seed)

  # SVC model
  model <- apply(cov_pars, 1, function(x) {
    RandomFields::RFsimulate(
      RandomFields::RMexp(x["var"], x["scale"]),
                          x = locs[, "x"], y = locs[, "y"])
  })

  model[[p+1]] <- RandomFields::RFsimulate(
    RandomFields::RMnugget(var = nugget),
                           x = locs[, "x"], y = locs[, "y"])
  sp.SVC <- Reduce(cbind, model)
  sp.SVC <- sp::SpatialPointsDataFrame(coords = sp.SVC@coords,
                                       data = sp.SVC@data)
  colnames(sp.SVC@data) <- c(paste0("SVC_", 1:p), "nugget")

  return(sp.SVC)
}
