\name{zerocens.em}
\alias{zerocens.em}
\alias{zerocens.mle}

\title{
Maximum likelihood estimation of the zero-censored model
}

\description{
Maximum likelihood estimation of the zero-censored model.
}

\usage{
zerocens.em(x, tol = 1e-6, maxit = 1000)
zerocens.mle(x)
}

\arguments{
\item{x}{
A numerical matrix with compositional data. Only one zero value is allowed in each row.
}
\item{tol}{
The tolerance value to terminate the EM algorithm.
}
\item{maxit}{
The maximum number of iterations allowed for the EM algorithm.
}
}

\details{
The function fits the cero-censored model (Tsagris, 2026) to compositional data with zero values.
The drawback is that only 1 zero, at most, is allowed in each compositional vector.
The zerocens.em() function fits the model using the EM algorithm and it is quite efficient, whereas the
second function, zerocens.mle() uses optim() and is quite slow.
}

\value{
A list including:
\item{loglik}{
The log-likelihood value.
}
\item{iters}{
The number of iterations required by the EM algorithm.
}
\item{mu}{
The estimated mean vector.
}
\item{sigma}{
The estimated covariance matrix.
}
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
zerocens.em(x)
zerocens.mle(x)
}


