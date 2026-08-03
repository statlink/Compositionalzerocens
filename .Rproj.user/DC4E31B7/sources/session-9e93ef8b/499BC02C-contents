prob.cens <- function(mu, sigma) {

  D <- length(mu) + 1
  setup <- .precompute_faces(D)
  tH <- t(setup$H)
  z_mean <- as.numeric(tH %*% mu) + 1
  z_cov  <- tH %*% sigma %*% t(tH)     # D x D, singular, rank D-1

  1 - as.numeric(TruncatedNormal::pmvnorm(mu = z_mean, sigma = z_cov, lb = rep(0, D), ub = rep(Inf, D)))
}
