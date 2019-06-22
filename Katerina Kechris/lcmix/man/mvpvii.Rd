\name{mvpvii}

\alias{MVPVII}
\alias{mvpvii}
\alias{dmvpvii}
\alias{rmvpvii}

\title{
The Multivariate Pearson Type VII (PVII) Distribution
}

\description{
Density and random generation functions for the multivariate Pearson Type VII (PVII) distribution.
}

\usage{
dmvpvii(x, mean, scale, shape, log=FALSE)
rmvpvii(n, mean, scale, shape)
}

\arguments{
\item{x}{a numeric matrix of which each row represents an observation.}
\item{mean}{a vector of mean parameters for the columns of \code{x}.  Let \eqn{D}{D} = \code{ncol(x)}, and \code{length(mean)} should be equal to \eqn{D}{D}.}
\item{scale}{a positive definite \eqn{D \times D}{D-by-D} matrix.}
\item{shape}{a positive scalar.}
\item{log}{logical; if \code{TRUE}, density is given as the log-density.}
\item{n}{number of vectors to simulate.}
}

\details{
The multivariate PVII distribution, a generalization of the multivariate t-distribution, arises when the inverse of the covariance of a multivariate normal random variable is itself a random variable with the Wishart distribution.  See Sun et al. (2010) for details.  As parameterized here, the density of a multivariate PVII random variable with mean \eqn{\mu}{mu}, scale \eqn{\Sigma}{Sigma}, and shape \eqn{\alpha}{alpha} is

\deqn{f(x) = (2\pi)^{-D/2} \ |\Sigma|^{-1/2} \ \Gamma(\alpha)^{-1} \ \Gamma\left(\alpha + \frac{D}{2}\right) \ \left\{1 + \frac{1}{2}(x-\mu)^{\mbox{\tiny T}} \Sigma^{-1} (x-\mu)\right\}^{-\alpha + \frac{D}{2}}}{f(x) = (2 pi)^(-D/2) |Sigma|^(-1/2) gamma(alpha)^(-1) gamma(alpha + D/2) {1 + (1/2) t(x-mu) Sigma^(-1) (x-mu)}^(alpha + D/2) }

where \eqn{\Gamma}{\dQuote(gamma)} denotes the gamma function.
}

\value{
For \code{dmvpvii}, a vector of densities.  For \code{rmvpvii}, a vector with \code{n} rows and \eqn{D}{D} columns representing a sample from the multivariate PVII distribution with the specified parameters.
}

\references{
Sun, J., Kaban, A., and Garibaldi, J.M.  (2010)  Robust mixture clustering using Pearson Type VII distribution.  \emph{Pattern Recognition Letters} \bold{31}(16), 2447--2454.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{pvii}} for the univariate version; \code{\link{mvnorm}} for a related distribution; \code{\link{thetahat}} for parameter estimation.
}

\examples{
set.seed(123)

mu <- -1:1

Sigma <- mleCov(matrix(rnorm(30), ncol=3))
Sigma
#            [,1]       [,2]       [,3]
# [1,]  0.8187336  0.5147059 -0.3243663
# [2,]  0.5147059  0.9698367 -0.4933797
# [3,] -0.3243663 -0.4933797  0.7797652

alpha <- 2

x <- rmvpvii(5, mu, Sigma, alpha)
x
#            [,1]       [,2]     [,3]
# [1,] -0.7774939 -0.5824543 2.139000
# [2,] -0.3941455  0.5651861 1.157972
# [3,] -0.3595201  0.2209538 0.588348
# [4,] -1.4053874 -0.5759132 1.055372
# [5,] -1.6673451  1.3083343 1.625087

dmvpvii(x, mu, Sigma, alpha)
# [1] 0.024770962 0.103843268 0.147672666 0.189416293 0.001638714
}

\keyword{distribution}
\keyword{multivariate}
