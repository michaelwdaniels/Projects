\name{model reporting functions}

\alias{reporting}

\alias{assignment}
\alias{assignment.mixmod}
\alias{design}
\alias{design.mixmod}
\alias{iter}
\alias{iter.mixmod}
\alias{posterior}
\alias{posterior.mixmod}
\alias{stats}
\alias{stats.mixmod}

\title{
Mixture Model Reporting
}

\description{
See Details for a description of the individual functions.
}

\usage{
\method{assignment}{mixmod}(x, ...)
\method{design}{mixmod}(x, ...)
\method{iter}{mixmod}(x, ...)
\method{posterior}{mixmod}(x, ...)
\method{stats}{mixmod}(x, ...)
}

\arguments{
\item{x}{an object of class \code{\link{mixmod}} or \code{\link{mdmixmod}}.}
\item{...}{currently unused.}
}

\details{
These functions provide convenience reporting for objects of class \code{\link{mixmod}} or \code{\link{mdmixmod}}:  \code{assignment} returns the vector of most probable component assignments; \code{design} returns the design matrix, that is, a binary matrix of which the \eqn{(n,k)}{(n,k)}th element is 1 if the \eqn{n}{n}th datum is assigned to component \eqn{k}{k}, otherwise 0; \code{iter} returns the number of iterations used to fit the model; \code{posterior} returns the posterior probability matrix; \code{stats} returns a vector of model statistics.  Note that in the case of \code{\link{mdmixmod}} objects, \dQuote{components} refers here to the top-level classification, that is, \eqn{\mathcal{Y}_0}{Y_0}.
}

\value{
A numeric matrix, vector, or scalar.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{mixmod}} and \code{\link{mdmixmod}} for model details; \code{\link{likelihood}} for likelihood-specific model reporting; \code{\link{rocinfo}} for ROC curves on models.
}

\examples{\dontrun{
data(CiData)
fit <- mixmod(CiData$expression, 3)
design(fit)[1:10,]
#      [,1] [,2] [,3]
# [1,]    0    1    0
# [2,]    0    0    1
# [3,]    0    1    0
# [4,]    1    0    0
# [5,]    0    1    0
# [6,]    0    1    0
# [7,]    1    0    0
# [8,]    1    0    0
# [9,]    0    0    1
#[10,]    0    0    1
}}

\keyword{models}
\keyword{classif}