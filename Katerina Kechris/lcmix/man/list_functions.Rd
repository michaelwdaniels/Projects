\name{list functions}

\alias{listGrid}

\alias{listMean}
\alias{listSum}

\alias{namedDataFrame}
\alias{namedList}
\alias{namedVector}

\title{
Convenience Functions For Lists
}

\description{
Convenience structures for creating, manipulating, and reporting on lists.
}

\usage{
listGrid(x, KEEP.OUT.ATTRS=TRUE, stringsAsFactors=FALSE)

listMean(x, w=NULL)
listSum(x)

namedDataFrame(..., stringsAsFactors=FALSE) 
namedList(...)
namedVector(...)
}

\arguments{
\item{x}{for \code{listGrid}, a list of which each element is a vector; for \code{listMean} and \code{listSum}, a list of numeric scalars, vectors, matrices, or data frames all having the same dimension.}
\item{KEEP.OUT.ATTRS, stringsAsFactors}{as the equivalent arguments to \code{expand.grid}.}
\item{w}{optional vector of weights.}
\item{...}{any number of objects.}
}

\details{
\code{listGrid} is a replacement for \code{expand.grid} in which the last factor, rather than the first as in \code{expand.grid}, varies fastest.

\code{listMean} and \code{listSum} return the (optionally weighted) mean and sum, respectively, of a list of numeric elements.  All elements of \code{x} must have the same dimension.

The \code{named*} functions return named objects of the appropriate type containing the objects in the \code{...} argument.  For example, instead of writing \code{x <- list(y=y, z=z)}, one may write \code{x <- namedList(x, y)}.
}

\value{
For \code{listGrid}, a data frame.  For \code{listMean} and \code{listSum}, an object of the same dimension as the elements of \code{x}.  For \code{named*}, an object of the appropriate class.
}

\author{
Daniel Dvorkin
}

\examples{
species <- c("cat", "dog")
status <- c("spayed", "neutered", "intact")
age <- c(2, 5, 7)
listGrid(namedList(species, status, age))
#    species   status age
# 1      cat   spayed   2
# 2      cat   spayed   5
# 3      cat   spayed   7
# 4      cat neutered   2
# 5      cat neutered   5
# 6      cat neutered   7
# 7      cat   intact   2
# 8      cat   intact   5
# 9      cat   intact   7
# 10     dog   spayed   2
# 11     dog   spayed   5
# 12     dog   spayed   7
# 13     dog neutered   2
# 14     dog neutered   5
# 15     dog neutered   7
# 16     dog   intact   2
# 17     dog   intact   5
# 18     dog   intact   7

set.seed(123)
x <- replicate(3, matrix(rnorm(6), ncol=3), simplify=FALSE)
x
# [[1]]
#            [,1]       [,2]      [,3]
# [1,] -0.5604756 1.55870831 0.1292877
# [2,] -0.2301775 0.07050839 1.7150650
# 
# [[2]]
#            [,1]       [,2]      [,3]
# [1,]  0.4609162 -0.6868529 1.2240818
# [2,] -1.2650612 -0.4456620 0.3598138
# 
# [[3]]
#           [,1]       [,2]       [,3]
# [1,] 0.4007715 -0.5558411  0.4978505
# [2,] 0.1106827  1.7869131 -1.9666172
listSum(x)
#           [,1]      [,2]      [,3]
# [1,]  0.301212 0.3160143 1.8512200
# [2,] -1.384556 1.4117596 0.1082617
listMean(x, 1:3)
#            [,1]       [,2]       [,3]
# [1,]  0.2606119 -0.2470868  0.6785005
# [2,] -0.4047086  0.7566540 -0.5775265

namedList(species, status, age)
# $species
# [1] "cat" "dog"
# 
# $status
# [1] "spayed"   "neutered" "intact"  
# 
# $age
# [1] 2 5 7
}

\keyword{list}