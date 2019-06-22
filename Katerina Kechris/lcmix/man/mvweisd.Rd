\name{mvweisd}

\alias{mvweisd}
\alias{dmvweisd}
\alias{rmvweisd}

\title{
The Multivariate Weibull (Shape-Decay) Distribution
}

\description{
Density and random generation functions for the multivariate Weibull (Shape-Decay) distribution constructed using a normal (Gaussian) copula.
}

\usage{
dmvweisd(x, shape, decay, corr=diag(ncol(x)), log=FALSE)
rmvweisd(n, shape=1, decay=1, corr=diag(length(shape)))
}

\arguments{
\item{x}{a numeric matrix of which each row represents an observation.}
\item{shape}{a vector of shape parameters for the marginal distributions of the columns of \code{x}.  If \code{length(shape) < ncol(x)}, the elements of \code{shape} will be recycled.  If \code{length(shape) > ncol(x)}, the \code{shape} vector will be truncated and a warning given.}
\item{decay}{a vector of decay parameters for the marginal distributions of the columns of \code{x}.  If \code{length(decay) < ncol(x)}, the elements of \code{decay} will be recycled.  If \code{length(decay) > ncol(x)}, the \code{decay} vector will be truncated and a warning given.}
\item{corr}{the correlation matrix.  See Details.}
\item{log}{logical; if \code{TRUE}, density is given as the log-density.}
\item{n}{number of vectors to simulate.}
}

\details{
The construction of multivariate distributions from univariate marginal
distributions using normal copulas is discussed in Song (2000).  Briefly, given
univariate marginal densities and the corresponding distribution functions
(here, the shape-decay parameterization of the Weibull distribution), the
standard normal quantiles of the values of the distribution functions follow a
multivariate standard normal distribution, that is, a multivariate normal
distribution with marginal means of 0 and marginal variances of 1.  Thus the
covariance matrix is referred to as the correlation matrix in this context.
} % line breaks inside paragraphs prevent Rd2txt from truncating lines

\value{
For \code{dmvweisd}, a vector of densities.  For \code{rmvweisd}, a vector with \code{n} rows and \code{ncol(corr)} columns representing a sample from the multivariate Weibull (shape-decay) distribution with the specified parameters.
}

\references{
Song, P.  (2000)  Multivariate dispersion models generated from Gaussian copula.  \emph{Scandinavian Journal of Statistics} \bold{27}, 305--320.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{weisd}} for the underlying univariate distribution; \code{\link{mvnorm}}, \code{\link{mvgamma}} for related distributions; \code{\link{thetahat}} for parameter estimation.
}

\examples{
set.seed(123)

s <- 1:2
d <- 2:3
rho <- matrix(c(1, 0.5, 0.5, 1), ncol=2)

x <- rmvweisd(5, s, d, rho)
print(x)
#           [,1]      [,2]
# [1,] 0.1600585 0.3834426
# [2,] 1.3762076 0.6174464
# [3,] 0.6280634 1.0148760
# [4,] 0.3958020 0.2199443
# [5,] 0.1229824 0.3249533


dmvweisd(x, s, d, rho)
# [1] 2.6471540 0.1836727 0.2035865 0.8923407 2.9891143
}

\keyword{distribution}
\keyword{multivariate}
