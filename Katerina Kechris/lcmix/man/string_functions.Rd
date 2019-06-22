\name{string functions}

\alias{collapseVec}
\alias{notification}
\alias{upperFirst}

\title{
Functions for Creating or Manipulating Strings
}

\description{
See Details for a description of the individual functions.
}

\usage{
collapseVec(x, space=TRUE, quote=FALSE, as.expr=FALSE)
notification(x, depth=0)
upperFirst(x)
}

\arguments{
\item{x}{a vector for \code{collapseVec}, otherwise a string.}
\item{space}{logical; if \code{TRUE}, items will be separated by spaces after commas.}
\item{quote}{logical; if \code{TRUE}, items will be quoted.}
\item{as.expr}{logical; if \code{TRUE}, return a string which can be parsed as an expression.}
\item{depth}{a non-negative integer.}
}

\details{
\code{collapseVec} returns a string representing a vector, intended for printing or for entry as a value.  The elements of the vector (\dQuote{items}) are separated by commas and enclosed by parentheses.

\code{notification} returns \code{x} with the time and date appended, and \code{depth} tabs prepended.

\code{upperFirst} returns a version of \code{x} with the first letter (if any) in the string in upppercase.
}

\value{A string.}
\author{
Daniel Dvorkin
}

\examples{
collapseVec(1:3) # "(1, 2, 3)"
collapseVec(1:3, as.expr=TRUE) # "c(1, 2, 3)"
evpat(collapseVec(1:3, as.expr=TRUE))^2 # [1] 1 4 9

notification("hello, world")
# "hello, world @ 2012-10-08 19:10:21 MDT"

message(notification("goodbye, cruel world", 2))
#         goodbye, cruel world @ 2012-10-08 19:11:05 MDT

upperFirst(collapseVec(c("spam", "eggs"))) # "(Spam, eggs)"
}

\keyword{manip}
\keyword{character}
