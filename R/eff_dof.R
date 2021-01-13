tr <- function(A) {
  # computes the trace of a matrix
  sum(diag(A))
}


eff_dof <- function(cov.par, X, cov_func, outer.W, taper) {

  n <- nrow(X)
  p <- length(outer.W)
  nug.var <- cov.par[length(cov.par)]
  Sigma <- Sigma_y(cov.par, p, cov_func, outer.W, taper)

  iSigma <- solve(Sigma)

  XtiS <- crossprod(X, iSigma)

  # trace of hat matrix
  as.numeric(tr(X %*% solve(XtiS %*% X) %*% XtiS))

  # as.numeric(nug.var * tr(solve(XtiS %*% X) %*% XtiS%*% iSigma %*% X) +
  #   n + nug.var * tr(iSigma))

}

tr_Sigma <- function(cov.par, X, cov_func, outer.W, taper) {

  n <- nrow(X)
  nug.var <- tail(cov.par, 1)
  p <- length(outer.W)


  tr(Sigma_y(cov.par, p, cov_func, outer.W, taper)) - nug.var*n
}

