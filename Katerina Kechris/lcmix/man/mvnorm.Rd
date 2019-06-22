\name{mvnorm}

\alias{mvnorm}
\alias{dmvnorm}
\alias{rmvnorm}

\title{
The Multivariate Normal Distribution
}

\description{
Density and random generation functions for the multivariate normal distribution.
}

\usage{
dmvnorm(x, mean=rep(0, ncol(x)), cov=diag(ncol(x)), log=FALSE)
rmvnorm(n, mean=NULL, cov=NULL) 
}

\arguments{
\item{x}{a numeric matrix of which each row represents an observation.}
\item{mean}{a vector of column means.}
\item{cov}{a covariance matrix.}
\item{log}{logical; if \code{TRUE}, density is given as the log-density.}
\item{n}{number of vectors to simulate.}
}

\details{
See Evans et al. (2000) for one of many references on the multivariate normal density.

For \code{rmvnorm}, if the \code{mean} argument is \code{NULL}, then the scalar 0 will be used by default, unless \code{cov} is not \code{NULL}, in which case \code{mean = rep(0, nrow(cov))} will be used.  If the \code{cov} argument is \code{NULL}, then \code{diag(length(mean))} will be used by default.  Thus \code{rmvnorm(n)} is equivalent to \code{rnorm(n)}.
}

\value{
For \code{dmvnorm}, a vector of densities.  For \code{rmvnorm}, a vector with \code{n} rows and \code{length(mean)} columns representing a sample from the multivariate normal distribution with the specified parameters.
}

\references{
Evans, M., Hastings, N., and Peacock, B.  (2000)  \emph{Statistical Distributions, third edition}, John Wiley & Sons.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{mvgamma}}, \code{\link{mvweisd}} for related distributions; \code{Normal} in package \code{stats}; \code{dmvnorm} and \code{rmvnorm} in package \code{mixtools}.
}

\examples{
set.seed(123)

m = c(-3, 0, 3)
v = cov(matrix(rnorm(90), ncol=3))

x = rmvnorm(5, m, v)
print(x)
#            [,1]       [,2]     [,3]
# [1,] -2.0498637  0.3599609 3.218045
# [2,] -3.7479240  1.2302183 2.290916
# [3,] -0.9852752  1.1433559 2.790147
# [4,] -3.9353966 -0.5451158 3.223321
# [5,] -3.2864769 -0.1672326 2.201353

dmvnorm(x, m, v)
# 0.048079901 0.025591976 0.002587824 0.041810685 0.054688032
}

\keyword{distribution}
\keyword{multivariate}
