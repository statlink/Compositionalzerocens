zerocens.mle <- function(x) {

  z <- Compositional::alfa(x, 1)$aff
  ind <- which( Rfast::rowsums(x == 0) == 0 )
  y1 <- z[ind, , drop = FALSE]
  n1 <- dim(y1)[1]
  w2 <- z[-ind, , drop = FALSE]
  n2 <- dim(w2)[1]
  D <- dim(x)[2]  ;  d <- D - 1
  n <- n1 + n2
  con1 <- n2 * (d - 1) * log(2 * pi)
  con2 <- n1 * d * log(2 * pi)

  ## pack/unpack mu and Sigma via Cholesky factor
  unpack <- function(par) {
    mu <- par[1:d]
    lvec <- par[(d + 1):length(par)]
    L <- matrix(0, d, d)
    L[lower.tri(L, diag = TRUE)] <- lvec
    diag(L) <- exp(diag(L))
    list(mu = mu, Sigma = crossprod(L))
  }

  mu0 <- Rfast::colmeans(z)
  S0 <- cov(z)
  L0 <- t( chol(S0) )
  diag(L0) <- log( diag(L0) )
  par0 <- c( mu0, L0[lower.tri(L0, diag = TRUE)] )

  loglik <- function(par, z, y1, d, D, w2, n2, n1, n) {
    PENALTY <- 1e10

    up <- unpack(par)
    me <- up$mu
    se <- up$Sigma

    se_inv <- tryCatch( solve(se), error = function(e) NULL )
    if ( is.null(se_inv) ) return(PENALTY)
    log_det_se <- tryCatch( as.numeric(determinant(se, logarithm = TRUE)$modulus),
                            error = function(e) NA )
    if ( !is.finite(log_det_se) ) return(PENALTY)

    ## --- Interior (uncensored) points
    m1 <- -0.5 * n1 * log_det_se - 0.5 * sum( Rfast::mahala(y1, me, se) )

    ## --- Censored (boundary) points
    lam <- numeric(n2)
    for ( j in 1:n2 ) {
      v <- w2[j, ]

      M <- cbind(v, diag(d))
      B <- qr.Q(qr(M))[, 1:d, drop = FALSE]
      if ( sum(B[, 1] * v) < 0 ) B[, 1] <- -B[, 1]
      a <- sum(v * B[, 1])
      mt <- as.numeric(crossprod(B, me))
      ts <- crossprod(B, se) %*% B

      my <- mt[-1]
      s12 <- ts[1, -1, drop = FALSE]
      s22 <- ts[-1, -1, drop = FALSE]

      s22inv <- tryCatch( solve(s22), error = function(e) NULL )
      if ( is.null(s22inv) ) return(PENALTY)
      log_det_s22 <- tryCatch( as.numeric(determinant(s22, logarithm = TRUE)$modulus),
                               error = function(e) NA )
      if ( !is.finite(log_det_s22) ) return(PENALTY)

      com <- s12 %*% s22inv
      mx  <- mt[1] - com %*% my
      sx  <- ts[1, 1] - com %*% t(s12)
      if ( !is.finite(sx) || sx <= 0 ) return(PENALTY)

      lam[j] <-  -0.5 * log_det_s22 - 0.5 * as.numeric(my %*% s22inv %*% my) +
                 pnorm( a, mx, sqrt(sx), lower.tail = FALSE, log.p = TRUE )
    }

    val <-  -m1 - sum(lam)
    if ( !is.finite(val) ) return(PENALTY)
    val
  }

  qa <- nlm( loglik, par0, z = z, y1 = y1, d = d, D = D, w2 = w2,
             n2 = n2, n1 = n1, n = n, iterlim = 10000 )
  qa <- optim( qa$estimate, loglik, z = z, y1 = y1, d = d, D = D, w2 = w2,
               n2 = n2, n1 = n1, n = n, control = list( maxit = 10000 ) )

  up <- unpack(qa$par)
  lik <-  -qa$value - 0.5 * con1 - 0.5 * con2 + (n * d + 0.5 * n) * log(D)

  list( loglik = lik, mu = up$mu, sigma = up$Sigma )
}

