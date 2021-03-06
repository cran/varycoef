% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SVC_mle.R
\name{predict.SVC_mle}
\alias{predict.SVC_mle}
\title{Prediction of SVCs (and response variable)}
\usage{
\method{predict}{SVC_mle}(
  object,
  newlocs = NULL,
  newX = NULL,
  newW = NULL,
  compute.y.var = FALSE,
  ...
)
}
\arguments{
\item{object}{(\code{SVC_mle}) \cr
Model obtained from \code{\link{SVC_mle}} function call.}

\item{newlocs}{(\code{NULL} or \code{matrix(n.new, 2)}) \cr
If \code{NULL}, then function uses observed locations of model to estimate
SVCs. Otherwise, these are the new locations the SVCs are predicted for.}

\item{newX}{(\code{NULL} or \code{matrix(n.new, q)}) \cr
If provided (together with \code{newW}), the function also returns the
predicted response variable.}

\item{newW}{(\code{NULL} or \code{matrix(n.new, p)}) \cr
If provided (together with \code{newX}), the function also returns the
predicted response variable.}

\item{compute.y.var}{(\code{logical(1)}) \cr
If \code{TRUE} and the response is being estimated, the predictive
variance of each estimate will be computed.}

\item{...}{further arguments}
}
\value{
The function returns a data frame of \code{n.new} rows and with
columns
\itemize{
  \item \code{SVC_1, ..., SVC_p}: the predicted SVC at locations \code{newlocs}.
  \item \code{y.pred}, if \code{newX} and \code{newW} are provided
  \item \code{y.var}, if \code{newX} and \code{newW} are provided and
  \code{compute.y.var} is set to \code{TRUE}.
  \item \code{loc_x, loc_y}, the locations of the predictions
}
}
\description{
Prediction of SVCs (and response variable)
}
\examples{
## ---- toy example ----
## sample data
# setting seed for reproducibility
set.seed(123)
m <- 7
# number of observations
n <- m*m
# number of SVC
p <- 3
# sample data
y <- rnorm(n)
X <- matrix(rnorm(n*p), ncol = p)
# locations on a regular m-by-m-grid
locs <- expand.grid(seq(0, 1, length.out = m),
                    seq(0, 1, length.out = m))

## preparing for maximum likelihood estimation (MLE)
# controls specific to MLE
control <- SVC_mle_control(
  # initial values of optimization
  init = rep(0.1, 2*p+1),
  # lower bound
  lower = rep(1e-6, 2*p+1),
  # using profile likelihood
  profileLik = TRUE
)

# controls specific to optimization procedure, see help(optim)
opt.control <- list(
  # number of iterations (set to one for demonstration sake)
  maxit = 1,
  # tracing information
  trace = 6
)

## starting MLE
fit <- SVC_mle(y = y, X = X, locs = locs,
               control = control,
               optim.control = opt.control)

## output: convergence code equal to 1, since maxit was only 1
summary(fit)

## prediction
# new location
newlocs <- matrix(0.5, ncol = 2, nrow = 2)

# new data
X.new <- matrix(rnorm(2*p), ncol = p)

# predicting SVCs
predict(fit, newlocs = newlocs)

# predicting SVCs and calculating response
predict(fit, newlocs = newlocs,
        newX = X.new, newW = X.new)

# predicting SVCs, calculating response and predictive variance
predict(fit, newlocs = newlocs,
        newX = X.new, newW = X.new,
        compute.y.var = TRUE)

}
\references{
Dambon, J. A., Sigrist, F., Furrer, R. (2021)
   \emph{Maximum likelihood estimation of spatially varying coefficient
   models for large data with an application to real estate price prediction},
   Spatial Statistics \doi{10.1016/j.spasta.2020.100470}
}
\seealso{
\code{\link{SVC_mle}}
}
\author{
Jakob Dambon
}
