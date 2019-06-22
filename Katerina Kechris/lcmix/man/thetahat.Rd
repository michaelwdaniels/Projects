\name{thetahat}

\alias{thetahat}
\alias{thetahat.norm}
\alias{thetahat.mvnorm}
\alias{thetahat.pvii}
\alias{thetahat.mvpvii}
\alias{thetahat.weisd}
\alias{thetahat.mvweisd}
\alias{thetahat.gamma}
\alias{thetahat.mvgamma}
\alias{thetahat.exp}
\alias{thetahat.mvexp}

\concept{parameter estimation}
\concept{maximum likelihood estimation}
\concept{maximum likelihood estimator}

\title{
Parameter Estimation for Common Distributions
}

\description{
Functions to estimate the parameters of distributions used in the \code{lcmix} package by maximum likelihood estimation.
}

\usage{
thetahat.norm(x, w=1, aslist=TRUE)
thetahat.mvnorm(x, w=1, aslist=TRUE)
thetahat.pvii(x, w=1, aslist=TRUE, iter.max=LC_ITER_MAX, theta=NULL)
thetahat.mvpvii(x, w=1, aslist=TRUE, iter.max=LC_ITER_MAX, theta=NULL)
thetahat.weisd(x, w=1, aslist=TRUE, shape.guess=NULL)
thetahat.mvweisd(x, w=1, aslist=TRUE, shape.guess=NULL)
thetahat.gamma(x, w=1, aslist=TRUE, shape.guess=NULL)
thetahat.mvgamma(x, w=1, aslist=TRUE, shape.guess=NULL)
thetahat.exp(x, w=1, aslist=TRUE)
thetahat.mvexp(x, w=1, aslist=TRUE)
}

\arguments{
\item{x}{for univariate, a numeric vector, with each element representing an observation; for multivariate, a numeric matrix or data frame, with each row representing an observation.}
\item{w}{a numeric vector of non-negative values; weights to be applied to the elements (univariate) or rows (multivariate) of \code{x} in parameter estimation.}
\item{aslist}{logical; if \code{TRUE}, estimated parameters will be returned as a named list.}
\item{iter.max}{maximum iterations for the \dQuote{internal} EM algorithm used in estimation of parameters for these distributions.  If \code{iter.max=0}, then only a 1-step estimation will be carried out, which is sufficient inside an \dQuote{external} EM algorithm such as used by \code{\link{mixmod}} and \code{\link{mdmixmod}}.}
\item{theta}{optional list containing initial guesses for parameters.  If not \code{NULL}, then for \code{thetahat.pvii} and \code{thetahat.mvpvii} the elements of \code{theta} should be named \code{mean}, \code{scale}, and \code{shape}.}
\item{shape.guess}{optional initial guess at the shape parameter:  a scalar for univariate, a vector of the appropriate length for multivariate.}
}

\details{
By default, estimation is unweighted.  If weights are given, the length of the weights should generally be the same as the length of the data, that is, \code{length(w) = length(x)} (univariate) or \code{length(w) = nrow(x)} (multivariate).  If weights are shorter than the data, they will be recycled.  If longer than the data, weights will be truncated and a warning will be given.

If \code{aslist} is \code{FALSE}, parameters will be returned as a vector with named elements if possible.  However, multivariate distributions generally have matrix parameters and using \code{aslist=FALSE} for these distributions may produce results which are difficult to interpret.
}

\value{
A named list, a vector, or a matrix, depending on the distribution and options used, of estimated parameters.
}

\references{
Lawless, J.F.  (1982)  \emph{Statistical Models and Methods for Lifetime Data}, John Wiley & Sons.

Joe, H. and Xu, J.  (1996)  The estimation method of inference functions for margins for multivariate models.  Technical Report 166, University of British Columbia, Department of Statistics.
}

\author{
Daniel Dvorkin
}

\note{
For distributions without closed-form MLEs, a method-of-moments approach is used to obtain initial parameter estimates and then numerical optimization of the likelihood function is performed using \code{nlminb} in package \code{stats}.  See Lawless (1982) for issues in parameter estimation for the Weibull and gamma distributions.

For multivariate distributions other than the normal and PVII, the correlation matrix for the copula is estimated using the inference function for margins procedure of Joe and Xu (1996).
}

\seealso{
\code{\link{mvnorm}}, \code{\link{pvii}}, \code{\link{mvpvii}}, \code{\link{weisd}}, \code{\link{mvweisd}}, \code{\link{mvgamma}}, \code{\link{mvexp}}; \code{nlminb}, \code{Normal}, \code{Weibull}, \code{Gamma}, \code{Exponential} in package \code{stats}; \code{fitdistr} in package \code{MASS} which performs parameter estimation using different numerical optimization functions.
}

\examples{
set.seed(123)
x <- rweisd(1000, 2, 3)
thetahat.weisd(x, aslist=FALSE)
#    shape    decay
# 2.015511 2.996586
}

\keyword{distribution}
\keyword{univar}
\keyword{multivariate}