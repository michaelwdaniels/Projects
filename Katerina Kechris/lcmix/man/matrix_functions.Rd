\name{matrix functions}

\alias{allDiv}
\alias{colDiv}
\alias{rowDiv}

\alias{asSymmetric}

\alias{colSample}
\alias{rowSample}

\alias{colSort}
\alias{rowSort}

\alias{isRowMax}
\alias{isRowMin}
\alias{matrix2design}

\alias{labsDet}

\alias{matrixMaxIdx}
\alias{matrixMinIdx}

\alias{matrixTrace}

\title{
Convenience Functions For Matrices
}

\description{
Convenience functions for manipulating and reporting on matrices.
}

\usage{
allDiv(x, log=FALSE)
colDiv(x, log=FALSE)
rowDiv(x, log=FALSE)

asSymmetric(x, keep.lower=TRUE)

colSample(x, size, replace=FALSE, prob=NULL)
rowSample(x, size, replace=FALSE, prob=NULL)

colSort(x, decreasing=FALSE, na.last=NA)
rowSort(x, decreasing=FALSE, na.last=NA)

isRowMax(x)
isRowMin(x)
matrix2design(x)

labsDet(x)

matrixMaxIdx(x)
matrixMinIdx(x)

matrixTrace(x)
}

\arguments{
\item{x}{a numeric matrix or data frame (square for \code{asSymmetric}, \code{labsDet}, and \code{matrixTrace}; of any dimension for the other functions).}
\item{log}{logical; if \code{TRUE}, the natural logarithm of the result will be returned.}
\item{keep.lower}{logical; if \code{TRUE}, the lower triangular portion of \code{x} will be used, otherwise the upper triangular portion.}
\item{size, replace, prob}{as the equivalent arguments to \code{sample}.}
\item{decreasing, na.last}{as the equivalent arguments to \code{sort}.}
}

\details{
\code{allDiv} returns \code{x} divided by the sum of all its elements, while \code{colDiv} returns the columns of \code{x} divided by their sums, and \code{rowDiv} returns the rows of \code{x} divided by their sums.  Because calculations are performed by subtraction of logarithms when \code{log} is \code{TRUE}, this option may be useful for precision, e.g., \code{exp(rowDiv(x, TRUE))} may be more precise than \code{rowDiv(x)}.

If \code{x} is square, \code{asSymmetric} mirrors the lower triangular portion (by default) of \code{x} to the upper portion and returns a symmetric matrix.  \code{matrixTrace} returns the trace of the matrix, that is, the sum of the elements along the diagonal.  No checking is done, so it is the user's responsibility to ensure that \code{x} is square; the function will \dQuote{work} with non-square matrices as well!

\code{colSample} and \code{rowSample} return a matrix consisting of samples from the columns or rows, respectively, of \code{x}.

\code{colSort} and \code{rowSort} return \code{x} with its columns sorted by the sums of its columns or its rows sorted by the sums of its rows, respectively.

\code{isRowMax} returns a matrix of the same dimension as \code{x} of which the \eqn{(n,k)}{(n,k)}th element is \code{TRUE} if the \eqn{(n,k)}{(n,k)}th element of \code{x} is the maximum value in the \eqn{k}{k}th row, otherwise \code{FALSE}.  \code{isRowMin} is the same for row minima instead of maxima.  For \code{matrix2design}, the result is the same as \code{isRowMax} except that numerical values (1 or 0) are used instead of logical values.

\code{labsDet} returns the logarithm of the absolute value of the determinant of \code{x}.  Like the built-in function \code{det}, this is a convenience wrapper around \code{determinant}.

\code{matrixMaxIdx} and \code{matrixMinIdx} return the row and column indices of the maximum and minimum elements, respectively, of \code{x}.
}

\value{
For \code{allDiv}, \code{colDiv}, \code{rowDiv}, \code{asSymmetric}, \code{isRowMax}, \code{isRowMin}, \code{matrix2design}, \code{colSort}, and \code{rowSort}, a matrix of the same dimension as \code{x}.  For \code{colSample} and \code{rowSample}, a matrix of the specified dimension.  For \code{matrixMaxIdx} and \code{matrixMinIdx}, a two-element vector.  For \code{matrixTrace}, a scalar.
}

\author{
Daniel Dvorkin
}

\examples{
x <- matrix(1:12, ncol=3)
x
#     [,1] [,2] [,3]
#[1,]    1    5    9
#[2,]    2    6   10
#[3,]    3    7   11
#[4,]    4    8   12
rowDiv(x)
#           [,1]      [,2]      [,3]
#[1,] 0.06666667 0.3333333 0.6000000
#[2,] 0.11111111 0.3333333 0.5555556
#[3,] 0.14285714 0.3333333 0.5238095
#[4,] 0.16666667 0.3333333 0.5000000
rowDiv(x, log=TRUE)
#          [,1]      [,2]       [,3]
#[1,] -2.708050 -1.098612 -0.5108256
#[2,] -2.197225 -1.098612 -0.5877867
#[3,] -1.945910 -1.098612 -0.6466272
#[4,] -1.791759 -1.098612 -0.6931472
matrixMaxIdx(x)
# [1] 4 3
}

\keyword{matrix}