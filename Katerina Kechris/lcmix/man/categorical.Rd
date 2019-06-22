\name{categorical}

\alias{categorical}
\alias{dcategorical}
\alias{rcategorical}

\title{
The Categorical Distribution
}

\description{
Density and random generation functions for the categorical distribution.
}

\usage{
dcategorical(x, prob, log=FALSE)
rcategorical(n, prob)
}

\arguments{
\item{x}{vector of integers between 1 and some integer \eqn{K}, inclusive.}
\item{prob}{vector of probabilities such that \code{sum(prob) = 1} and \code{length(prob)} = \eqn{K}.}
\item{log}{logical; if \code{TRUE}, density is given as the log-density.}
\item{n}{number of observations.}
}

\details{
The categorical distribution is a special case of the multinomial distribution with size 1.  Values may be between 1 and some integer \eqn{K}, inclusive.  If the random variable \eqn{\mathcal X}{X} is categorical with probability vector \eqn{p} of length \eqn{K}, then \eqn{P(\mathcal X = x) = p_x}{P(X=x) = p_x}.  The density is

\deqn{f(x) \ = \ \prod_k p_k^{I(x=k)} \ = \ p_x}{
f(x) = prod_k p_k^I(x=k) = p_x}

where \eqn{I(.)} is the indicator function and \eqn{k = 1,\dots,K}.
}

\value{
\code{dcategorical} gives the density and \code{rcategorical} generates random deviates.
}

\author{
Daniel Dvorkin
}

\seealso{
\code{Multinom} in package \code{stats}.
}

\keyword{distribution}
\keyword{univar}
