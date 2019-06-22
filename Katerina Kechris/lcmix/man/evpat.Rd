\name{evpat}

\alias{evpat}

\title{
Evaluate Parsable Text
}

\description{
See Details.
}

\usage{
evpat(x)
}

\arguments{
\item{x}{a string containing a parsable expression.}
}

\details{
This function attempts to evaluate \code{x} as a command.
}

\value{The return value of the expression in \code{x}.}

\author{
Daniel Dvorkin
}

\examples{
evpat("c(1, 2, 3)") + 1 # [1] 2 3 4
}
