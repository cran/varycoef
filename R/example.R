
#' Sample Function for GP-based SVC Model for Given Locations
#'
#' @description Samples SVC data at given locations. The SVCs parameters and the
#' covariance function have to be provided. The sampled model matrix can be 
#' provided or it is sampled. The SVCs are sampled according to their given parametrization and at
#' respective observation locations. The error vector is sampled from a nugget
#' effect. Finally, the response vector is computed. Please note that the
#' function is not optimized for sampling large data sets.
#'
#' @param df.pars (\code{data.frame(p, 3)}) \cr
#'    Contains the mean and covariance parameters of SVCs. The three columns
#'    must have the names \code{"mean"}, \code{"var"}, and \code{"scale"}.
#' @param nugget.sd  (\code{numeric(1)}) \cr
#'    Standard deviation of the nugget / error term.
#' @param cov.name (\code{character}(1)) \cr
#'    Character defining the covariance function, c.f. \code{\link{SVC_mle_control}}.
#' @param locs (\code{numeric(n)} or \code{matrix(n, d)}) \cr
#'    The numeric vector or matrix contains the observation locations and
#'    therefore defines the number of observations to be \code{n}. For a vector,
#'    we assume locations on the real line, i.e., \eqn{d=1}.
#' @param X (\code{NULL} or \code{matrix(n, p)}) \cr
#'    If \code{NULL}, the covariates are sampled, where the first column contains 
#'    only ones to model an intercept and further columns are sampled from a 
#'    standard normal. If it is provided as a \code{matrix}, then the dimensions
#'    must match the number of locations in \code{locs} (\code{n}) and the number of SVCs
#'    defined by the number of rows in \code{df.pars} (\code{p}).
#'
#' @return \code{list} \cr
#'    Returns a list with the response \code{y}, model matrix
#'    \code{X}, a matrix \code{beta} containing the sampled SVC at given
#'    locations, a vector \code{eps} containing the error, and a matrix
#'    \code{locs} containing the original locations. The \code{true_pars}
#'    contains the data frame of covariance parameters that were used to
#'    sample the GP-based SVCs. The nugget variance has been added to the 
#'    original argument of the function with its respective variance, but 
#'    \code{NA} for \code{"mean"} and \code{"scale"}.
#'
#' @details The parameters of the model can be chosen such that we obtain data
#'    from a not full model, i.e., not all covariates are associated with a 
#'    fixed and a random effect. Using \code{var = 0} for instance yields a 
#'    constant beta coefficient for respective covariate. Note that in that 
#'    case the \code{scale} value is neglected.
#'
#' @examples
#' set.seed(123)
#' # SVC parameters
#' (df.pars <- data.frame(
#'    var = c(2, 1),
#'    scale = c(3, 1),
#'    mean = c(1, 2)))
#' # nugget standard deviation
#' tau <- 0.5
#'
#' # sample locations
#' s <- sort(runif(500, min = 0, max = 10))
#' SVCdata <- sample_SVCdata(
#'   df.pars = df.pars, nugget.sd = tau, locs = s, cov.name = "mat32"
#' )
#' @importFrom spam rmvnorm cov.exp cov.mat cov.sph cov.wend1 cov.wend2
#' @importFrom stats rnorm
#' @export
sample_SVCdata <- function(
    df.pars, nugget.sd, locs,
    cov.name = c("exp", "sph", "mat32", "mat52", "wend1", "wend2"),
    X = NULL
) {
  # transform to matrix for further computations
  if (is.vector(locs)) {
    locs <- matrix(locs, ncol = 1)
  }
  # check covariance parameters and locations
  stopifnot(
    is.data.frame(df.pars),
    all(df.pars$var >= 0),
    all(df.pars$scale > 0),
    nugget.sd > 0,
    is.matrix(locs)
  )
  # dimensions
  d <- dim(locs)[2]
  n <- dim(locs)[1]
  p <- nrow(df.pars)
  
  ## build SVC models depending on covariance function, i.e., Sigma_y
  D <- as.matrix(dist(locs, diag = TRUE, upper = TRUE))
  
  ## covariance functions
  cov_fun <- function(theta) {
    do.call(
      what = MLE.cov.func(cov.name),
      args = list(h = D, theta = theta)
    )
  }
    
  
  ## sample SVCs (including mean effect)
  beta <- apply(df.pars, 1, function(x) {
    if (x["var"] == 0) {
      rep(x["mean"], n)
    } else {
      spam::rmvnorm(
        n = 1,
        mu = rep(x["mean"], n), 
        Sigma = cov_fun(theta = x[c("scale", "var")])
      )
    }
  })
  # nugget
  eps <- rnorm(n, sd = nugget.sd)
  
  # data
  if (is.null(X)) {
    X <- cbind(1, matrix(rnorm(n*(p-1)), ncol = p-1))
  } else {
    stopifnot(
      is.matrix(X),
      dim(X)[1] == n,
      dim(X)[2] == p
    )
  }
  y <- apply(beta*X, 1, sum) + eps
  
  list(
    y = y, X = X, beta = beta, eps = eps, locs = locs, 
    true_pars = rbind(
      df.pars, 
      data.frame(
        var = nugget.sd^2, 
        scale = NA, mean = NA
      )
    )
  )
}
