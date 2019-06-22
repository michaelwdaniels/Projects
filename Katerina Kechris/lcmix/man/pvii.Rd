\name{pvii}

\alias{PVII}
\alias{pvii}
\alias{dpvii}
\alias{rpvii}

\title{
The Pearson Type VII (PVII) Distribution
}

\description{
Density and random generation functions for the Pearson Type VII (PVII) distribution.
}

\usage{
dpvii(x, mean, scale, shape, log=FALSE)
rpvii(n, mean, scale, shape)
}

\arguments{
\item{x}{a numeric vector of observations.}
\item{mean}{mean parameter (real number).}
\item{scale}{scale parameter (positive real number).}
\item{shape}{shape parameter (positive real number).}
\item{log}{logical; if \code{TRUE}, density is given as the log-density.}
\item{n}{number of observations to simulate.}
}

\details{
The PVII distribution, a generalization of the t-distribution, arises when the inverse of the variance of a normal random variable is itself a random variable with the gamma distribution.  See Sun et al. (2010) for details.  See \code{\link{mvpvii}} for the multivariate density which easily reduces to the univariate version.
}

\value{
For \code{dpvii}, a vector of densities.  For \code{rpvii}, a vector of length \code{n} representing a sample from the PVII distribution with the specified parameters.
}

\references{
Sun, J., Kaban, A., and Garibaldi, J.M.  (2010)  Robust mixture clustering using Pearson Type VII distribution.  \emph{Pattern Recognition Letters} \bold{31}(16), 2447--2454.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{mvpvii}} for the multivariate version; \code{\link{thetahat}} for parameter estimation.
}

\examples{
set.seed(123)
mu <- 1
Sigma <- 3
alpha <- 2
x <- rpvii(5, mu, Sigma, alpha)
x
# [1]  1.8452349 -0.2040292 -2.1309664  0.4013371  2.0181986
dpvii(x, mu, Sigma, alpha)
# [1] 0.23112248 0.17824559 0.02719683 0.26484673 0.20555914
}

\keyword{distribution}
\keyword{univar}
