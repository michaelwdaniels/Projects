\name{mvexp}

\alias{mvexp}
\alias{dmvexp}
\alias{rmvexp}

\title{
The Multivariate Exponential Distribution
}

\description{
Density and random generation functions for the multivariate exponential distribution constructed using a normal (Gaussian) copula.
}

\usage{
dmvexp(x, rate, corr=diag(ncol(x)), log=FALSE)
rmvexp(n, rate=1, corr=diag(length(rate)))
}

\arguments{
\item{x}{a numeric matrix of which each row represents an observation.}
\item{rate}{a vector of rate parameters for the marginal distributions of the columns of \code{x}.  If \code{length(rate) < ncol(x)}, the elements of \code{rate} will be recycled.  If \code{length(rate) > ncol(x)}, the \code{rate} vector will be truncated and a warning given.}
\item{corr}{the correlation matrix.  See Details.}
\item{log}{logical; if \code{TRUE}, density is given as the log-density.}
\item{n}{number of vectors to simulate.}
}

\details{
The construction of multivariate distributions from univariate marginal distributions using normal copulas is discussed in Song (2000).
Briefly, given univariate marginal densities and the corresponding distribution functions (here, the exponential distribution), the standard normal quantiles of the values of the distribution
functions follow a multivariate standard normal distribution, that is, a multivariate normal distribution with marginal means of 0 and marginal variances of 1.
Thus the covariance matrix is referred to as the correlation matrix in this context.
} % line breaks inside paragraphs prevent Rd2txt from truncating lines

\value{
For \code{dmvexp}, a vector of densities.  For \code{rmvexp}, a vector with \code{n} rows and \code{ncol(corr)} columns representing a sample from the multivariate exp distribution with the specified parameters.
}

\references{
Song, P.  (2000)  Multivariate dispersion models generated from Gaussian copula.  \emph{Scandinavian Journal of Statistics} \bold{27}, 305--320.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{mvnorm}}, \code{\link{mvweisd}} for related distributions; \code{\link{thetahat}} for parameter estimation.
}

\examples{
set.seed(123)

r <- 2:3
rho <- matrix(c(1, 0.5, 0.5, 1), ncol=2)

x <- rmvexp(5, r, rho)
x
#           [,1]        [,2]
#[1,] 0.05569265 0.025983192
#[2,] 0.02328782 0.009049892
#[3,] 0.10527587 0.003011938
#[4,] 0.11576801 0.227582536
#[5,] 0.10297173 0.330607553


dmvexp(x, r, rho)
# [1] 10.389920 18.910258  7.690541  2.501572  1.478549
}

\keyword{distribution}
\keyword{multivariate}
