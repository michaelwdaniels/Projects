\name{isMDF}

\alias{isMDF}

\title{
Test For Matrix Or Data Frame
}

\description{
See Value.
}

\usage{
isMDF(x)
}

\arguments{
\item{x}{an object.}
}

\value{\code{TRUE} if \code{x} is a matrix or a data frame, otherwise \code{FALSE}.
}

\author{
Daniel Dvorkin
}

\examples{
x = 1:3
y = cbind(1:3, 2:4)
z = as.data.frame(y)
isMDF(x) # FALSE
isMDF(y) # TRUE
isMDF(z) # TRUE
isMDF(z$V1) # FALSE
}
