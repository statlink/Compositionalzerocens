\name{gof.zerocens}
\alias{gof.zerocens}

\title{
Goodness--of--fit test for the zero-censored model
}

\description{
Goodness--of--fit test for the zero-censored model.
}

\usage{
gof.zerocens(mu, sigma, B = 999, theoretical = TRUE, nsim = 1e+6)
}

\arguments{
\item{mu}{
The mean vector in \eqn{R^{D-1}}.
}
\item{sigma}{
The covariance matrix in \eqn{R^{D-1}}.
}
\item{theoretical}{
Do you want the theoretical probability?
}
\item{n}{
The number of draws to use in the Monte Carlo estimation of the zero-censoring probability.
}
}

\details{
The function computes the probability of a zero-censoring, in either face of the simplex,
using Monte-Carlo (Tsagris and Alharbi, 2026).
}

\value{
The (theoretical) probability of zero-censoring at each component.
}

\references{
Tsagris M. and Alharbi N. (2026). Modelling structural zeros in compositional data via
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
prob.zerocens(mu, sigma, theoretical = TRUE)
}


