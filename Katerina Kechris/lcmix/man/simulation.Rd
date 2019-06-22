\name{simulation}
\alias{simulation}
\alias{simulateMixdata}
\alias{simulateMdmixdata}
\alias{simulateFromFit}
\alias{simulateFromFit.mixmod}
\alias{simulateFromFit.mdmixmod}
\alias{resampleFromFit}
\alias{resampleFromFit.mixmod}
\alias{resampleFromFit.mdmixmod}
\alias{simulateCiTypeData}

\title{
Generate Simulated Data
}

\description{
Simulation functions allow the generation of simulated data suitable for model fitting with \code{\link{mixmod}} or \code{\link{mdmixmod}} from either a \emph{de novo} set of parameters or the parameters of a fitted model.
}

\usage{
simulateMixdata(n, distn, params)
simulateMdmixdata(n, distn, params, topology=LC_TOPOLOGY)
\method{simulateFromFit}{mixmod}(x, n=x$N, ...)
\method{simulateFromFit}{mdmixmod}(x, n=x$N, ...)
\method{resampleFromFit}{mixmod}(x, n=x$N,
    replace=TRUE, hidpar=x$params$hidden, ...)
\method{resampleFromFit}{mdmixmod}(x, n=x$N,
    replace=TRUE, hidpar=x$params$hidden, topology=x$topology, ...)
simulateCiTypeData(n=10000, par=LC_SIMPAR, topology=LC_TOPOLOGY)
}

\arguments{
\item{n}{the length of the data to be simulated, that is, the number of independent simulated observations.}
\item{distn}{the name of the distribution for \code{simulateMixdata}, or the vector of names of distributions for \code{simulateMdmixdata}, to be used in generating the observed portion of the data.  Must be one of the distribution names (\emph{not} family names) in \code{\link{LC_FAMILY}}.}
\item{params}{the parameters of the complete data (hidden and observed) distributions; a list formatted as the \code{params} element of an object of class \code{mixmod} or \code{mdmixmod}.}
\item{topology}{one of the model topologies in \code{\link{LC_TOPOLOGY}},  either \code{"layered"} or \code{"chained"}.}
\item{x}{an object of class \code{mixmod} or \code{mdmixmod}.}
\item{replace}{logical; if \code{TRUE}, sample with replacement.}
\item{hidpar}{hidden data simulation parameters.}
\item{par}{marginal simulation parameters with the same structure as \code{\link{LC_SIMPAR}}.}
\item{...}{currently unused.}
}

\details{
\code{simulateMixdata} and \code{simulateMdmixdata} are for generating data from \emph{de novo} parameters, while \code{simulateFromFit} generates data from the parameters of a fitted model, and \code{resampleFromFit} generates hidden data from the parameters of a fitted model and uses the data to which that model was fitted to generate the observed data.  The \code{params} argument must be a list formatted in the same way as the \code{params} element of a \code{mixmod} or \code{mdmixmod} object, such as is returned by \code{\link{mixmod}} or \code{\link{mdmixmod}}, respectively.  However, in the case of simulating a multiple data model from a layered topology, that is, \code{simulateMdmixdata(..., topology="layered")}, the \code{probz} and \code{rprob} elements of \code{params$hidden} are not required.

\code{simulateCiTypeData} is a convenience function to generate simulated data with performance characteristics similar to those of \code{\link{CiData}}.  The distribution is considerably simplified, but fitted models should show similar results, with either the layered or chained topology being favored depending on the topology chosen for the simulation.
}

\value{
For \code{simulateMixdata} and \code{simulateFromFit.mixmod}, a list with elements:
\item{X}{observed data.}
\item{Y}{hidden data.}

For \code{simulateMdmixdata}, \code{simulateFromFit.mdmixmod}, and \code{simulateCiTypeData}, a list with elements:
\item{X}{observed data.}
\item{Y}{intermediate-level hidden data, that is, observations on for \eqn{\mathcal Y_1, \ldots, Y_Z}{Y_1, ..., Y_Z}.}
\item{Y0}{top-level hidden data, that is, observations on \eqn{Y_0}{Y_0}.}
}

\author{
Daniel Dvorkin
}

\seealso{
\code{\link{mixmod}}, \code{\link{mdmixmod}} for specifications of the parameters; \code{\link{constants}} for distribution names and topologies; \code{\link{mvnorm}}, \code{\link{mvweisd}}, \code{\link{weisd}}, \code{\link{mvgamma}} for relevant distributions; \code{\link{rocinfo}} for performance evaluations using complete data.
}

\examples{\dontrun{

set.seed(123)

SimMarginal <- simulateMixdata(1000, distn="norm",
    params=list(hidden=list(prob=c(0.1, 0.9)),
                observed=list(mean=c(1, -1), sd=c(1, 2))))
SimMarginalFit <- mixmod(SimMarginal$X, 2)
SimMarginalFit
# Normal mixture model ('norm')
# Data 'SimMarginal$X' of length 1000 fitted to 2 components
# Model statistics:
#      iter      llik      qval       bic    iclbic 
#   406.000 -2115.747 -2302.603 -4266.032 -4639.745 
plot(rocinfo(SimMarginalFit, SimMarginal$Y==1))
                 

SimJoint <- simulateCiTypeData(10000) # layered topology
SimJointFits <- lapply(namedList("layered", "chained"), function(top)
    mdmixmod(SimJoint$X, c(2,3,2), topology=top))
SimJointFits
# $layered
# Layered (normal, normal, normal) mixture model ('norm', 'mvnorm', 'norm')
# Data 'SimJoint$X' of size 10000-by-(1,3,1) fitted to 2 (2,3,2) components
# Model statistics:
#       iter       llik       qval        bic     iclbic 
#      55.00  -72885.59  -85368.43 -146176.44 -171142.12 
# 
# $chained
# Chained (normal, normal, normal) mixture model ('norm', 'mvnorm', 'norm')
# Data 'SimJoint$X' of size 10000-by-(1,3,1) fitted to 2 (2,3,2) components
# Model statistics:
#       iter       llik       qval        bic     iclbic 
#      35.00  -72887.69  -82999.68 -146189.84 -166413.83 
SimJointMarginalFits <- marginals(SimJointFits$layered)
SimJointAllFits <- c(SimJointFits, SimJointMarginalFits)
plot(multiroc(SimJointAllFits, SimJoint$Y0==1))
}}

\keyword{datasets}
\keyword{datagen}
