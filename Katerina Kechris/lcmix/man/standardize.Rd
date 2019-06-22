\name{standardize}

\alias{standardize}
\alias{standardise}

\title{
Standardize Data
}

\description{
Standardize data by mean or standard deviation, with optional centering.
}

\usage{
standardize(X, method=c("sd","mean"), mle=TRUE, center=FALSE)
standardise(X, method=c("sd","mean"), mle=TRUE, center=FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{X}{A numeric vector, matrix, or data frame.}
\item{method}{a string giving the standardization method to use, either \code{"sd"} (the default) or \code{"mean"}.  Partial matches are allowed.}
\item{mle}{logical; if \code{TRUE}, then the maximum likelihood estimator of standard deviation rather than the unbiased estimator will be used.  Has no effect if \code{method} is \code{"mean"}.}
\item{center}{logical; if \code{TRUE}, then after standardization, data will be centered to have mean (or column means) 0.}
}

\details{
\code{standardize} is an alternative to the base \R function \code{scale}, offering more options.  Data on a \eqn{(-\infty, \infty)}{(-Inf, Inf)} scale, such as normally distributed data, are generally standardized by dividing by the standard deviation.  However, it is often more appropriate to standardize data on a \eqn{[0, \infty)}{[0, Inf)} scale, such as Weibull or gamma distributed data, by dividing by the mean.

So far as the numerical values of the data are concerned, \code{standardize(X, "sd", FALSE, TRUE)} is equivalent to \code{scale(X)}.  See the note on attributes in Value, below.

\code{standardise} is an alias for \code{standardize}.
}

\value{
A vector or matrix of the same dimension as \code{X}, with the attribute \code{standardized}, which is a list with elements \code{scale} and \code{center}.  These have the same meanings as \code{scaled:scale} and \code{scaled:center} in the return value attributes of \code{scale}.  If \code{center} is \code{FALSE}, then \code{standardized@center} will be \code{NULL}.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{scale} in package \code{base}.
}

\examples{
set.seed(123)

X = rmvweisd(5, 1:3, 3:1)
standardize(X, "mean")
#           [,1]      [,2]      [,3]
# [1,] 0.4600131 0.8392135 1.3481948
# [2,] 1.0188317 1.0361609 1.3975700
# [3,] 1.5355348 0.3813876 0.6260400
# [4,] 0.5390977 1.7176610 0.9614062
# [5,] 1.4465227 1.0255769 0.6667891
# attr(,"standardized")
# attr(,"standardized")$scale
# [1] 0.2457062 0.6110324 1.0480777
# 
# attr(,"standardized")$center
# NULL

colMeans(standardize(X, "mean"))
# [1] 1 1 1
}

\keyword{univar}
\keyword{multivariate}
