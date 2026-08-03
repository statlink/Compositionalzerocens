prob.zerocens <- function(mu, sigma, n = 1e+6) {

   x <- Compositionalzerocens::rzerocens(n, mu, sigma)
   sum(x == 0) / n

}
