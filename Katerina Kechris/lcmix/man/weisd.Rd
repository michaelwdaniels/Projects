\name{weisd}

\alias{weisd}
\alias{dweisd}
\alias{pweisd}
\alias{qweisd}
\alias{rweisd}
\alias{eweisd}
\alias{vweisd}
\alias{medweisd}
\alias{rmomweisd}

\title{
The Weibull (Shape-Decay) Distribution
}

\description{
Density, distribution function, quantile function, random generation, expected value, variance, median, and raw moments for the shape-decay parameterization of the Weibull distribution.
}

\usage{
dweisd(x, shape, decay, log=FALSE)
pweisd(q, shape, decay, lower.tail=TRUE, log.p=FALSE)
qweisd(p, shape, decay, lower.tail=TRUE, log.p=FALSE)
rweisd(n, shape, decay)
eweisd(shape, decay)
vweisd(shape, decay)
medweisd(shape, decay)
rmomweisd(r, shape, decay)
}

\arguments{
\item{x,q}{vector of quantiles.}
\item{p}{vector of probabilities.}
\item{n}{number of observations. If \code{length(n) > 1}, the length is taken to be the number required.}
\item{shape, decay}{shape and decay parameters.  See Details.}
\item{log, log.p}{logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log p}{log(p)}.}
\item{lower.tail}{logical; if \code{TRUE} (default), probabilities are \eqn{P(X \leq x)}{P[X <= x]}, otherwise, \eqn{P(X > x)}{P[X > x]}.}
\item{r}{the \code{r}th raw moment is returned.  See Value.}
}

\details{
The \code{weisd} distribution with \code{shape} parameter \eqn{\sigma}{s} and \code{decay} parameter \eqn{\delta}{d} has density

\deqn{f(x) = \sigma \delta x^{\sigma-1} \exp(-\delta x^\sigma)}{f(x) = s*d*x^(s-1) * exp(-d*x^s)}

for \eqn{x \geq 0; \  \sigma, \delta > 0}{x >= 0; s, d > 0}.  Compared to the usual shape-scale parameterization used in \code{dweibull} etc., with \code{shape} parameter \eqn{\alpha}{a} and \code{scale} parameter \eqn{\beta}{b}, the relationship between the parameters is given by \eqn{\alpha = \sigma}{a = s} and \eqn{\beta = \delta^{-1/\sigma}}{b = d^(-1/s)}, or equivalently, \eqn{\sigma = \alpha}{s = a} and \eqn{\delta = \beta^{-\alpha}}{d = b^(-a)}.
}

\value{
\code{dweisd} gives the density, \code{pweisd} gives the distribution function, \code{qweisd} gives the quantile function, and \code{rweisd} generates random deviates.  \code{eweisd} gives the expected value, \code{vweisd} gives the variance, and \code{medweisd} gives the median for  the specified parameters.  \code{rmomweisd} gives the \code{r}th raw moment for the specified parameters; for example, \code{rmomweisd(1, 2, 3)} is equivalent to \code{eweisd(2, 3)}.
}

\author{
Daniel Dvorkin
}

\note{
The distribution function, using the parameters given in Details, has the very simple form

\deqn{F(x) = 1 - \exp(-\delta x^\sigma)}{F(x) = 1 - exp(-d*x^s)}

from which other functions such as the hazard function are easily derived.
}

\seealso{
\code{\link{mvweisd}} for the multivariate version; \code{\link{thetahat}} for parameter estimation; \code{Weibull} in package \code{stats}.
}

\keyword{distribution}
\keyword{univar}
