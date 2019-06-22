\name{covariance etc.}

\alias{mleCov}
\alias{mleSd}
\alias{mleVar}

\title{
Covariance and Standard Deviation MLE
}

\description{
Calcuate the maximum likelihood estimators of covariance and standard deviation.
}

\usage{
mleCov(x, mean=colMeans(as.matrix(x)), simplify=TRUE)
mleSd(x, mean=colMeans(as.matrix(x)))
}

\arguments{
\item{x}{a numeric vector, matrix, or data frame.}
\item{mean}{the mean, or vector of means, with respect to which covariance should be calculated.}
\item{simplify}{logical; if \code{TRUE}, the return value will be reduced to the simplest possible form, for example, a single-element matrix will be reduced to a scalar.}
}

\details{
The standard \R functions \code{cov}, \code{var}, and \code{sd} return unbiased estimators of (co)variance and standard deviation, which are not equal to the maximum likelihood estimator (MLE).  The functions given here return MLEs instead.  For example, if \code{x} is a vector consisting of elements \eqn{x_1, \dots, x_N} with mean \eqn{\mu}{mu}, then the MLE of the variance is

\deqn{\hat{\sigma^2} = \frac{1}{N} \sum_{n=1}^N (x_n - \mu)^2}{sigmasq.hat = (1/N) sum_n (x_n - mu)^2}

while the unbiased estimator is

\deqn{\hat{\sigma^2} = \frac{1}{N-1} \sum_{n=1}^N (x_n - \mu)^2.}{sigmasq.hat = {1/(N-1)} sum_n (x_n - mu)^2.}

\code{mleVar} is an alias for \code{mleCov}.
}

\value{
For \code{mleCov}, the estimated variance if \code{x} is a vector, or the esitmated covariance matrix of the columns of \code{x} if \code{x} is a matrix or data frame.  If \code{simplify} is \code{FALSE}, then a matrix will always be returned even if \code{x} is a vector.

For \code{mleSd}, the estimated standard deviation if \code{x} is a vector, or the vector of estimated column standard deviations if \code{x} is a matrix or a data frame.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{thetahat.norm}}, \code{\link{thetahat.mvnorm}} for weighted parameter estimation; \code{cov}, \code{var}, \code{sd} in package \code{stats}.
}

\examples{
set.seed(123)

x <- rnorm(10)
mleCov(x)
# [1] 0.8187336
var(x)
# [1] 0.909704
mleSd(x)
# [1] 0.904839
sd(x)
# [1] 0.953784

x <- matrix(rnorm(30), ncol=3)
mleCov(x)
#            [,1]       [,2]       [,3]
# [1,]  0.9698367 -0.4933797  0.1336627
# [2,] -0.4933797  0.7797652 -0.1931557
# [3,]  0.1336627 -0.1931557  0.2502430
cov(x)
#            [,1]       [,2]       [,3]
# [1,]  1.0775964 -0.5481997  0.1485141
# [2,] -0.5481997  0.8664058 -0.2146175
# [3,]  0.1485141 -0.2146175  0.2780478
mleSd(x)
# [1] 0.9848029 0.8830432 0.5002430
apply(x, 2, sd) # sd(x) is deprecated
# [1] 1.0380734 0.9308092 0.5273024
}

\keyword{matrix}