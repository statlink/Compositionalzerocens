zerocens.em <- function(x, tol = 1e-6, maxit = 1000) {

  z <- Compositional::alfa(x, 1)$aff
  ind <- which( Rfast::rowsums(x == 0) == 0 )
  y1 <- z[ind, ] 
  sy1 <- Rfast::colsums(y1) 
  sy12 <- crossprod(y1)

  w2 <- z[-ind, , drop = FALSE]
  n1 <- dim(y1)[1] ;  n2 <- dim(w2)[1]
  D <- dim(x)[2]  ;  d  <- D - 1
  n <- n1 + n2
  e1 <- c(1, rep(0, d - 1))

  con1 <- n2 * (d - 1) * log(2 * pi)
  con2 <- n1 * d * log(2 * pi)

  a_vec <- sqrt( Rfast::rowsums(w2^2) )
  ## --- precompute B_j ONCE, outside the EM loop ---
  B_list <- vector("list", n2)
  for ( j in 1:n2 )  B_list[[j]] <- far::orthonormalization(cbind(w2[j, ], e1))

  loglik_val <- function(mu, Sigma) {
    se_inv <- solve(Sigma)
    log_det_se <- as.numeric( determinant(Sigma, logarithm = TRUE)$modulus )
    m1 <-  - 0.5 * n1 * log_det_se - 0.5 * sum( Rfast::mahala(y1, mu, Sigma) )
    lam <- numeric(n2)
    for ( j in 1:n2 ) {
      B <- B_list[[ j ]]                 # <-- reused, not recomputed
      a <- a_vec[j]
      mt <- mu %*% B
      ts <- crossprod(B, Sigma) %*% B
      my <- mt[-1]  ;  s12 <- ts[1, -1]  ;  s22 <- ts[-1, -1]
      s22inv <- solve(s22)
      com <- s12 %*% s22inv
      mx <- as.numeric(mt[1] - com %*% my)
      sx <- as.numeric(ts[1, 1] - com %*% s12)
      lam[j] <-  -0.5 * as.numeric(determinant(as.matrix(s22), logarithm = TRUE)$modulus) -
                  0.5 * my %*% s22inv %*% my + pnorm(a, mx, sqrt(sx), lower.tail = FALSE, log.p = TRUE)
    }
    m1 + sum(lam) - 0.5 * con1 - 0.5 * con2 + (n * d + 0.5 * n) * log(D)
  }

  mu <- Rfast::colmeans(z)
  Sigma <- cov(z)
  ll_hist <- numeric(maxit)

  Ey <- matrix(0, n2, d)

  for ( it in 1:maxit ) {
    Eyy_sum <- 0
    for ( j in 1:n2 ) {
      v <- w2[j, ]
      a <- a_vec[j]
      B <- B_list[[ j ]]                  
      mt <- mu %*% B
      ts <- crossprod(B, Sigma) %*% B
      my  <- mt[-1]  ;  s12 <- ts[1, -1]  ;  s22 <- ts[-1, -1]
      s22inv <- solve(s22)
      com <- s12 %*% s22inv
      mx <- as.numeric( mt[1] - com %*% my )
      sx <- as.numeric( ts[1, 1] - com %*% s12 )
      sdx <- sqrt(sx)

      alpha_i <- (a - mx) / sdx
      log_phi <- dnorm(alpha_i, log = TRUE)
      log_Sbar <- pnorm(alpha_i, lower.tail = FALSE, log.p = TRUE)
      lambda_i <- exp(log_phi - log_Sbar)
      zhat1 <- mx + sdx * lambda_i
      var_trunc <- sx * ( 1 - lambda_i * (lambda_i - alpha_i) )
      Ez1sq <- var_trunc + zhat1^2

      r1 <- zhat1 / a
      r2 <- Ez1sq / a^2
      Ey[j, ] <- r1 * v
      Eyy_sum <- Eyy_sum + r2 * tcrossprod(v) ## ( v %*% t(v) )
    }

    mu <- ( sy1 + Rfast::colsums(Ey) ) / n
    S_total <- sy12 + Eyy_sum
    Sigma <- S_total / n - mu %*% t(mu)
    ll_hist[it] <- loglik_val(mu, Sigma)

    if ( it > 1 && abs(ll_hist[it] - ll_hist[it - 1]) < tol ) {
      ll_hist <- ll_hist[1:it]
      break
    }
  }

  list(loglik = ll_hist[it], iters = it, mu = mu, sigma = Sigma)
}
