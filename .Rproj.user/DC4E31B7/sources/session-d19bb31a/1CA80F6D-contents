gof.zerocens <- function(x, mu, sigma, B = 999, theoretical = TRUE, nsim = 1e6) {

  n <- dim(x)[1]  ;  D <- dim(x)[2]
  neg <- x == 0
  n_neg <- Rfast::rowsums(neg)                      # 0 = interior, 1 = recorded zero
  obs_counts <- c( sum(n_neg == 0), Rfast::colsums(neg) )
  names(obs_counts) <- c( "interior", paste0("X", 1:D) )
  p_hat <- Compositionalzerocens::prob.zerocens(mu, sigma, theoretical = theoretical, n = nsim)
  p_full <- c(1 - sum(p_hat), p_hat)
  E <- n * p_full
  X2_obs <- sum( (obs_counts - E)^2 / E )
  X2_sim <- numeric(B)
  for ( b in 1:B ) {
    y_b <- Compositionalzerocens::rzerocens(n, mu, sigma)
    neg_b <- y_b == 0
    nneg_b <- Rfast::rowsums(neg_b)
    counts_b <- c( sum(nneg_b == 0), Rfast::colsums(neg_b) )
    X2_sim[b] <- sum( (counts_b - E)^2 / E )
  }

  ( sum(X2_sim >= X2_obs) + 1 ) / (B + 1)

}
