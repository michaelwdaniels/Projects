\name{convergencePlot}

\alias{convergencePlot}

\concept{convergence}

\title{
Show EM Algorithm Convergence
}

\description{
Plot convergence of the log-likelihood to a maximum during execution of the EM algorithm for fitting mixture models.  Optionally also show the Q-value, that is, the expected complete data log-likelihood.
}
\usage{
convergencePlot(x, show.qval=FALSE, main="")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{an object of class \code{mixmod} or \code{mdmixmod}.}
\item{show.qval}{logical; if \code{TRUE}, the Q-value will be plotted as well as the log-likelihood.}
\item{main}{plot title.}
}

\details{
Log-likelihood is plotted vs. iteration in black.  If \code{show.qval} is \code{TRUE}, Q-value is plotted in red.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{mixmod}}, \code{\link{mdmixmod}}, \code{\link{rocinfo}}; \code{plot} in package \code{graphics}.
}

\examples{\dontrun{
data(CiData)
fit <- mixmod(CiData$expression, 3)
convergencePlot(fit, TRUE, main="log-likelihood and Q-value")
}}

\keyword{hplot}
