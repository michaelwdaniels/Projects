\name{rlapply}

\alias{rlapply}

\title{
Recursive List Apply
}

\description{
Recursively apply function \code{f} to a list \code{object} containing, at some level, objects of a class appearing in \code{classes}.  This is a simpler alternative to the base \R function \code{rapply}.
}

\usage{
rlapply(object, f, classes="ANY", ...)
}

\arguments{
\item{object, f, classes}{as equivalent arguments to \code{rapply}.}
\item{...}{further arguments to \code{f}.}
}

\value{A list with the same structure as \code{object}.}

\author{
Daniel Dvorkin
}

\examples{
set.seed(123)
x <- lapply(2:3, function(D)
        lapply(c(10, 20), function(N)
            matrix(rnorm(N*D), ncol=D)))
rlapply(x, dim, "matrix")
# [[1]]
# [[1]][[1]]
# [1] 10  2
# 
# [[1]][[2]]
# [1] 20  2
# 
# 
# [[2]]
# [[2]][[1]]
# [1] 10  3
# 
# [[2]][[2]]
# [1] 20  3
}

