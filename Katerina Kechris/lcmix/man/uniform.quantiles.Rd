\name{uniformQuantiles}

\alias{uniformQuantiles}

\title{
Evenly Spaced Sample Quantiles
}

\description{
Quantiles corresponding to the uniform distribution of the distribution function of a vector, or of the columns of a matrix or data frame.
}

\usage{
uniformQuantiles(X, K, decreasing=TRUE)
}

\arguments{
\item{X}{a numeric vector, matrix, or data frame.}
\item{K}{the number of quantiles.}
\item{decreasing}{logical; if \code{TRUE}, the quantiles are given in descending order.}
}

\details{
For a vector \eqn{x}, the \eqn{k}th uniform quantile, for \eqn{k = 1,\dots,K}, is the sample quantile at \eqn{q = \left(k-\frac{1}{2}\right)/K}{q = (k-0.5)/K}.  In other words, \code{uniformQuantiles(X, K, FALSE)[k] = quantile(X, (k-0.5)/K)} if \code{X} is a vector.  If \code{X} is a matrix or data frame, the function is applied to the columns of \code{X}.
}

\value{
A vector of quantiles if \code{X} is a vector, a matrix of quantiles of the columns of \code{X} if \code{X} is a matrix or data frame.
}

\author{
Daniel Dvorkin
}

\note{
\code{uniformQuantiles} is used to initialize component means in \code{\link{mixmod}} and \code{\link{mdmixmod}}.
}

\seealso{
\code{\link{mixmod}}, \code{\link{mdmixmod}}; \code{quantile} in package \code{stats}.
}

\examples{
set.seed(123)

X = matrix(rnorm(300), ncol=3)
uniformQuantiles(X, 5)
#             [,1]       [,2]        [,3]
# [1,]  1.26449867  1.0584323  1.55267383
# [2,]  0.51301442  0.3039516  0.56563159
# [3,]  0.06175631 -0.2258300  0.03591471
# [4,] -0.38719515 -0.6200279 -0.45546412
# [5,] -1.06822046 -1.2894076 -1.02776368
}

\keyword{univar}
\keyword{multivariate}
