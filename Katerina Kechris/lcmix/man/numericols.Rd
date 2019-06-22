\name{numericols}

\alias{numericols}

\title{
Numeric Columns Of Data Frame
}

\description{
Extract numerical columns from a data frame.
}

\usage{
numericols(x, include.logical=TRUE)
}

\arguments{
\item{x}{a data frame with at least one numerical column.}
\item{include.logical}{logical; if \code{TRUE}, logical columns will be considered numerical.}
}

\value{A data frame containing the numerical columns of \code{x}.}

\author{
Daniel Dvorkin
}

\examples{
x = data.frame(number=1:3, species=c("cat", "dog", "rat"),
        forsale=c(TRUE, FALSE, TRUE))
numericols(x)
#   number forsale
# 1      1    TRUE
# 2      2   FALSE
# 3      3    TRUE
numericols(x, FALSE)
#   number
# 1      1
# 2      2
# 3      3
}
