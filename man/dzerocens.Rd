\name{dzerocens}
\alias{dzerocens}

\title{
Density of the zero-censored model
}

\description{
Density of the zero-censored model.
}

\usage{
dzerocens(x, mu, sigma, logdens = FALSE)
}

\arguments{
\item{x}{
A vector or a matrix with compositional data.
}
\item{mu}{
The mean vector in \eqn{R^{D-1}}.
}
\item{sigma}{
The covariance matrix in \eqn{R^{D-1}}.
}
\item{logdens}{
If you want the log of the density set this TRUE, otherwise leave it FALSE.
}
}

\details{
The function computes the density values of the zero-censored model.
}

\value{
The density values at the given compositional data x.
}

\references{
Tsagris M. (2026). Modelling structural zeros in compositional data via
a zero-censored multivariate normal model.

https://arxiv.org/pdf/2208.13073
}

\author{
R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\examples{
mu <- c(0.325, 0.121)
sigma <- matrix(c(0.149, -0.200,
                  -0.200,  0.323), 2, 2)
x <- rzerocens(100, mu, sigma)
mod <- zerocens.em(x)
mu <- mod$mu
sigma <- mod$sigma
f <- dzerocens(x, mu, sigma)
}


