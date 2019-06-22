\name{forceRange}

\alias{forceRange}

\title{
Force Within Range
}

\description{
	Force all elements of the input value to be within the specified range.
}

\usage{
forceRange(x, min=-Inf, max=Inf)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{A numeric vector, matrix, or data frame, possibly containing values less that \code{min} or greater than \code{max}.}
\item{min}{The minimum value.}
\item{max}{The maximum value.}
}

\details{
\code{forceRange} will replace all elements of \code{x} less than \code{min} with \code{min}, and all elements of \code{x} greater than \code{max} with \code{max}.  If the default parameters are used, \code{x} will be returned unchanged.
}

\value{
A matrix, vector, or data frame.
}

\examples{
set.seed(123)

x = rnorm(12)
y = matrix(x, ncol=4)
y
#            [,1]       [,2]       [,3]       [,4]
# [1,] -0.5604756 0.07050839  0.4609162 -0.4456620
# [2,] -0.2301775 0.12928774 -1.2650612  1.2240818
# [3,]  1.5587083 1.71506499 -0.6868529  0.3598138

forceRange(x, 0)
# [1] 0.00000000 0.00000000 1.55870831 0.07050839 0.12928774 1.71506499
# [7] 0.46091621 0.00000000 0.00000000 0.00000000 1.22408180 0.35981383

forceRange(y, -1, 1)
#            [,1]       [,2]       [,3]       [,4]
# [1,] -0.5604756 0.07050839  0.4609162 -0.4456620
# [2,] -0.2301775 0.12928774 -1.0000000  1.0000000
# [3,]  1.0000000 1.00000000 -0.6868529  0.3598138
}