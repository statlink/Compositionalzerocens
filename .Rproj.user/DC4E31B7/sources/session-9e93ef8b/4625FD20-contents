\name{rzerocens}
\alias{rzerocens}

\title{
Maximum likelihood estimation of the zero-censored model
}

\description{
Maximum likelihood estimation of the zero-censored model.
}

\usage{
rzerocens(n, mu, sigma)
}

\arguments{
\item{n}{
The sample size.
}
\item{mu}{
The mean vector in \eqn{R^{D-1}}.
}
\item{sigma}{
The covariance matrix in \eqn{R^{D-1}}.
}
}

\details{
The function generates compositional data from the cero-censored model (Tsagris, 2026).
The drawback is that only 1 zero, at most, is allowed in each compositional vector.
}

\value{
A numerical matrix with compositional data.
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
x <- rzerocens(1000, mu, sigma)
}


