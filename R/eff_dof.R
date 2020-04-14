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

  as.numeric(nug.var * tr(solve(XtiS %*% X) %*% XtiS%*% iSigma %*% X) +
    n + nug.var * tr(iSigma))

}
