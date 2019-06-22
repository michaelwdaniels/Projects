\name{distance functions}

\alias{euclidean}
\alias{manhattan}

\title{
Euclidean and Manhattan Distance
}

\description{
Returns the squared Euclidean distance or the Manhattan distance between all rows in \code{x} and the vector \code{center}.
}

\usage{
euclidean(x, center, ...)
manhattan(x, center, ...)
}

\arguments{
\item{x}{a matrix or data frame, or a vector to be interpreted as a one-row matrix.}
\item{center}{a vector such that \code{length(center) = ncol(x)}.}
\item{\dots}{no effect; included for signature compatibility with \code{mahalanobis}.}
}

\value{
A vector of which the \code{i}th element is the \emph{squared} Euclidean distance or the Manhattan distance between \code{x[i,]} and \code{center}.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{mahalanobis} in package \code{stats}.
}

\examples{
set.seed(123)
x <- matrix(rnorm(12), ncol=3)
euclidean(x, rep(0,3))
# [1] 0.8026151 3.1930442 4.1403916 1.7348174
manhattan(x, rep(0,3))
# [1] 1.376616 2.390904 3.243706 1.695383
}

\keyword{array}
\keyword{arith}
