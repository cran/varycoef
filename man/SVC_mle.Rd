% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SVC_mle.R
\name{SVC_mle}
\alias{SVC_mle}
\alias{SVC_mle.default}
\alias{SVC_mle.formula}
\title{MLE of SVC model}
\usage{
SVC_mle(...)

\method{SVC_mle}{default}(y, X, locs, W = NULL, control = NULL, optim.control = list(), ...)

\method{SVC_mle}{formula}(
  formula,
  data,
  RE_formula = NULL,
  locs,
  control,
  optim.control = list(),
  ...
)
}
\arguments{
\item{...}{further arguments}

\item{y}{(\code{numeric(n)}) \cr
Response vector.}

\item{X}{(\code{matrix(n, q)}) \cr
Design matrix. Intercept has to be added manually.}

\item{locs}{(\code{matrix(n, d)}) \cr
Locations in a \eqn{d}-dimensional space. May contain multiple
observations at single location.}

\item{W}{(\code{NULL} or \code{matrix(n, p)}) \cr
If \code{NULL}, the same matrix as provided in \code{X} is used. This
fits a full SVC model, i.e., each covariate effect is modeled with a mean
and an SVC. In this case we have \eqn{p = q}. If optional matrix \code{W}
is provided, SVCs are only modeled for covariates within matrix \code{W}.}

\item{control}{(\code{list}) \cr
Control paramaters given by \code{\link{SVC_mle_control}}.}

\item{optim.control}{(\code{list}) \cr
Control arguments for optimization function, see Details in
\code{\link{optim}}.}

\item{formula}{Formula describing the fixed effects in SVC model. The response, i.e. LHS of the formula, is not allowed to have functions such as \code{sqrt()} or \code{log()}.}

\item{data}{data frame containing the observations}

\item{RE_formula}{Formula describing the random effects in SVC model. Only RHS is considered. If \code{NULL}, the same RHS of argument \code{formula} for fixed effects is used.}
}
\value{
Object of class \code{SVC_mle} if \code{control$extract_fun = FALSE},
meaning that a MLE has been conducted. Otherwise, if \code{control$extract_fun = TRUE},
the function returns a list with two entries:
\itemize{
   \item \code{obj_fun}: the objective function used in the optimization
   \item \code{args}: the arguments to evaluate the objective function.
}
For further detials, see description of \code{\link{SVC_mle_control}}.
}
\description{
Conducts a maximum likelihood estimation (MLE) for a Gaussian
  process-based SVC model as described in Dambon et al. (2021)
  \doi{10.1016/j.spasta.2020.100470}. More specifially, the model is
  defined as:

\deqn{y(s) = X \mu + W \eta (s) + \epsilon(s)}

where:
\itemize{
  \item \eqn{y} is the response (vector of length \eqn{n})
  \item \eqn{X} is the data matrix for the fixed effects covariates. The
  dimensions are \eqn{n} times \eqn{q}. This leads to \eqn{q} fixed effects.
  \item \eqn{\mu} is the vector containing the fixed effects
  \item W is the data matrix for the SVCs modeled by GPs. The dimensions are
  \eqn{n} times \eqn{p}. This lead to \eqn{p} SVCs in the model.
  \item \eqn{\eta} are the SVCs represented by a GP.
  \item \eqn{\epsilon} is the nugget effect
}

The MLE is an numeric optimization that runs \code{\link[stats]{optim}} or
(if parallelized) \code{\link[optimParallel]{optimParallel}}.
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
class(fit)

## output: convergence code equal to 1, since maxit was only 1
summary(fit)

## extract the optimization arguments, including objective function
control$extract_fun <- TRUE
opt <- SVC_mle(y = y, X = X, locs = locs,
               control = control)

# objective function and its arguments of optimization
class(opt$obj_fun)
class(opt$args)

# single evaluation with initial value
do.call(opt$obj_fun,
        c(list(x = control$init), opt$args))

\donttest{
## ---- real data example ----
require(sp)
## get data set
data("meuse", package = "sp")

# construct data matrix and response, scale locations
y <- log(meuse$cadmium)
X <- model.matrix(~1+dist+lime+elev, data = meuse)
locs <- as.matrix(meuse[, 1:2])/1000


## starting MLE
# the next call takes a couple of seconds
fit <- SVC_mle(y = y, X = X, locs = locs,
               # has 4 fixed effects, but only 3 random effects (SVC)
               # elev is missing in SVC
               W = X[, 1:3],
               control = SVC_mle_control(
                 # inital values for 3 SVC
                 # 7 = (3 * 2 covariance parameters + nugget)
                 init = c(rep(c(0.4, 0.2), 3), 0.2),
                 profileLik = TRUE
               ))

## summary and residual output
summary(fit)
plot(fit)

## predict
# new locations
newlocs <- expand.grid(
  x = seq(min(locs[, 1]), max(locs[, 1]), length.out = 30),
  y = seq(min(locs[, 2]), max(locs[, 2]), length.out = 30))
# predict SVC for new locations
SVC <- predict(fit, newlocs = as.matrix(newlocs))
# visualization
sp.SVC <- SVC
coordinates(sp.SVC) <- ~loc_x+loc_y
spplot(sp.SVC, colorkey = TRUE)
}
}
\references{
Dambon, J. A., Sigrist, F., Furrer, R. (2021)
   \emph{Maximum likelihood estimation of spatially varying coefficient
   models for large data with an application to real estate price prediction},
   Spatial Statistics \doi{10.1016/j.spasta.2020.100470}
}
\seealso{
\code{\link{predict.SVC_mle}}
}
\author{
Jakob Dambon
}
